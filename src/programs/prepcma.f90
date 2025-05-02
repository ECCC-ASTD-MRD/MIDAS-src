program midas_prepcma
  !
  !:Purpose: Read the observation files (usually after output by the background check) 
  !          and apply further quality control and thinning for use by the ``LETKF`` program 
  !
  !          ---
  !
  !:Algorithm: After reading the observation files that have been processed by the background check,
  !            the ``prepcma`` program rejects more observations based on the
  !            various flags and conditions. It also performs further thinnings
  !            on the data types of aircraft (AI), scatterometer (SC),
  !            satellite winds (SW) and some radiance (TO). The rejection and thinning
  !            are controlled by the options in the namelist of ``NAMPREPCMA``.
  !            
  !            ---
  !
  !:File I/O: The required input files and produced output files are listed as follows. 
  !
  !          --
  !
  !============================================== ==============================================================
  ! Input and Output Files                         Description of file
  !============================================== ==============================================================
  ! ``flnml``                                      In - Main namelist file with parameters user may modify
  ! ``flnml_static``                               In - The "static" namelist that should not be modified
  ! ``obserr``                                     In - Observation error statistics
  ! ``obsfiles_$FAM/obs$FAM_$NNNN_$NNNN``          In - Observation file for each "family"
  ! ``stats_tovs``                                 In - Satellite radiance observation errors
  ! ``rtcoef_$PLATFORM_$SENSOR.dat``               In - RTTOV coefficient files 
  ! ``obsfiles_$FAM.updated/obs$FAM_$NNNN_$NNNN``  Out - final observation file for each family
  !============================================== ==============================================================
  !
  !          ---
  !
  !:Synopsis: Below is a summary of the ``prepcma`` program calling sequence:
  !
  !           - **Initial setups:**
  !            
  !             - Read the NAMPREPCMA namelist and check/modify some values.            
  !
  !             - ``filt_setup``: set up list of elements to be assimilated and flags for rejection
  !                              
  !             - ``obsf_setup``: get observation file names and datestamp
  !
  !           - **Computation:**
  !
  !             - ``obsf_readFiles``: get the observations 
  !
  !             - ``filt_suprep``: select the elements to assimilate and apply rejection flags 
  !
  !             - ``oer_setObsErrors``: initialize obs error covariances and set flag 
  !
  !             - ``oti_setup``: reject any observations outside the data assimilation window
  !
  !             - ``enkf_rejectHighLatIR``: reject all IR radiance observation in arctic and antarctic
  !
  !             - ``enkf_modifyAmsubObsError``: modify the obs error stddev for AMSUB in the tropics 
  !
  !             - ``thinning_fam``: perform thinning for aircraft (AI), scatterometer (SC),
  !                                 satellite winds (SW) and some radiance (TO)
  !
  !           - **Final steps:**
  !
  !             - ``obsf_writeFiles``: write to burp/sqlite files 
  !
  !             - ``obsf_printFiles``: print to ascci file and to unformatted files
  !
  !          ---
  !
  !:Options: `List of namelist blocks <../namelists_in_each_program.html#prepcma>`_
  !          that can affect the ``prepcma`` program.
  !
  !          * Some of the relevant namelist blocks used to configure the
  !            prepcma are listed in the following table:
  !
  !          --
  !   
  !=================== ====================== ===========================================
  ! Program/Module     Namelist               Description of what is controlled
  !=================== ====================== ===========================================
  ! ``midas_prepcma``  ``NAMPREPCMA``         parameters for CMA format and others 
  !                                           to modify, reject and thinning some  
  !                                           observation data 
  ! ``timeCoord_mod``  ``NAMTIME``            assimilation time window length, temporal
  !                                           resolution of the background state and the 
  !                                           analysis
  ! ``tovs_nl_mod``    ``NAMTOV``             The list of satellite and instrument
  !=================== ====================== ===========================================
  !
  !
  use version_mod
  use obsSpaceData_mod
  use obsFiles_mod
  use obsFilter_mod
  use obsTimeInterp_mod
  use obsErrors_mod
  use tovs_mod
  use timeCoord_mod
  use enkf_mod
  use utilities_mod
  use midasMpi_mod
  use ramDisk_mod
  use regions_mod
  use burpRead_mod
  use gridStateVectorFileIO_mod
  use codtyp_mod
  use rttov_const, only : inst_name, platform_name

  implicit none

  integer :: fnom, ierr, dateStampFromObs
  type(struct_obs), target  :: obsSpaceData
  type(struct_oti), pointer :: oti => null()
  real(kind=8) :: hx_dummy(1,1)
  integer :: ncmahdr, ncmahx, ncmabdy, ncmadim, nobsout, nbrpform
  integer :: numHeader, numBody
  logical :: qcvar
  character(len=7) :: resumeType

  ! number of pressure ranges used for the thinning of aircraft (and other) data:
  integer, parameter :: npres_ai = 5
  integer, parameter :: npres_sw = 2
  integer, parameter :: nai_target = 10
  integer, parameter :: nsc_target = 10
  integer, parameter :: nsw_target = 6
  integer, parameter :: nto_target = 6 
  integer :: numTovsInstNameList, sensorIndex
  integer :: maxNumHeaderPerInst
  real(8) :: nai_pmax(npres_ai) = (/ 25000.0, 40000.0, 60000.0, 80000.0, 110000.0/)
  real(8) :: nsw_pmax(npres_sw) = (/ 60000.0, 110000.0/)
  ! For a scalar array, no layer selection will be done
  real(8) :: nsc_pmax(1) = (/ 0.0 /)
  real(8) :: nto_pmax(1) = (/ 0.0 /)
  character(len=codtyp_name_length) :: tovsInstName
  character(len=codtyp_name_length), allocatable :: tovsInstNameList(:)

  ! Namelist variables:
  integer :: maxNumHeadersForTovsInst(tvs_maxNumberOfSensors) ! max number of headers for each TOVS inst
  character(len=codtyp_name_length) :: tovsInstNamesWithMaxNumHeaders(tvs_maxNumberOfSensors) ! List of TOVS inst names 
                                      ! to prescribe max number of headers
  character(len=256) :: cmahdr        ! should not be used anymore
  character(len=256) :: cmabdy        ! should not be used anymore
  character(len=256) :: cmadim        ! should not be used anymore
  character(len=256) :: obsout        ! file name for ascii output
  character(len=256) :: brpform       ! should not be used anymore
  logical :: suprep                   ! choose to execute 'suprep' obs filtering
  logical :: rejectOutsideTimeWindow  ! choose to reject obs outside time window
  logical :: thinning                 ! choose to apply 'extra' thinning of some obs types
  logical :: thinningConv             ! choose to apply 'extra' thinning of conventional (AI+SW+SC) obs types
  logical :: thinningRadiance         ! choose to apply 'extra' thinning of radiance (TO) obs types
  logical :: applySatUtil             ! choose to reject satellite obs based on 'util' column of stats_tovs
  logical :: modifyAmsubObsError      ! choose to modify the obs error stddev for AMSUB/MHS in the tropics
  logical :: rejectHighLatIR          ! choose to reject IR data in high latitudes
  logical :: obsClean                 ! choose to remove rejected observations from files
  logical :: writeObsFiles            ! choose to update the (burp or sqlite) observation files
  logical :: writeAsciiCmaFiles       ! choose to write ascii output
  logical :: thinTovsPerInst          ! choose to thin tovs per instrument

  NAMELIST /NAMPREPCMA/ cmahdr, cmabdy, cmadim, obsout, brpform,  &
                        suprep, rejectOutsideTimeWindow, &
                        thinning, thinningConv, thinningRadiance, &
                        applySatUtil, modifyAmsubObsError, rejectHighLatIR, &
                        obsClean, writeObsFiles, writeAsciiCmaFiles, &
                        thinTovsPerInst, tovsInstNamesWithMaxNumHeaders, &
                        maxNumHeadersForTovsInst

  call ver_printNameAndVersion('prepcma','Prepare observations for LETKF')

  !- 1.0 mpi
  call mmpi_initialize

  !- 1.1 timings
  call tmg_init(mmpi_myid, 'TMG_INFO')
  call utl_tmg_start(0,'Main')
  call utl_printTime()

  if ( mmpi_myid == 0 ) call utl_writeStatus('PREPCMA_BEG')

  call utl_readNml()

  !- Specify default values for namelist variables
  cmahdr        = 'NOT_DEFINED'
  cmabdy        = 'NOT_DEFINED'
  cmadim        = 'NOT_DEFINED'
  obsout        = 'NOT_DEFINED'
  brpform       = 'brpform'
  suprep                  = .true.
  rejectOutsideTimeWindow = .true.
  thinning                = .true.
  thinningConv            = .true.
  thinningRadiance        = .true.
  applySatUtil            = .true.
  modifyAmsubObsError     = .true.
  rejectHighLatIR         = .true.
  obsClean                = .true.
  writeObsFiles           = .false.
  writeAsciiCmaFiles      = .false.
  thinTovsPerInst         = .false.
  tovsInstNamesWithMaxNumHeaders(:) = 'NOT_DEFINED'
  maxNumHeadersForTovsInst(:) = nto_target

  call utl_tmg_start(181,'low-level--readNML')
  read(utl_flnml, nml=namprepcma, iostat=ierr)
  if (ierr /= 0) call utl_abort('midas-prepcma: Error reading namelist')
  if (mmpi_myid == 0) write(*,nml=namprepcma)
  call utl_tmg_stop(181)

  !- RAM disk usage
  call ram_setup

  ! Setup the format of the output RPN standard files to 'XDF' or 'RSF'
  call gio_setup

  call utl_tmg_start(10,'--Observations')

  !- Set up list of elements to be assimilated and flags for rejection (from namelist)
  call filt_setup('prepcma')

  !- Observation file names and get datestamp
  call obsf_setup(dateStampFromObs, 'prepcma')

  !- Allocate obsSpaceData
  call obs_class_initialize('ENKFMIDAS')
  call obs_initialize( obsSpaceData, mpi_local_opt=obsf_filesSplit() )

  !- Read observations
  call utl_tmg_start(11,'----ReadObsFiles')
  call obsf_readFiles( obsSpaceData )
  call utl_tmg_stop(11)

  numHeader = obs_numheader(obsSpaceData)
  numBody   = obs_numbody(obsSpaceData)
  write(*,*) 'midas-prepcma: obs_numheader =', numheader
  write(*,*) 'midas-prepcma: obs_numbody   =', numbody

  !- Determine if qcvar flag is expected to be present
  resumeType = brpr_getTypeResume() 
  write(*,*) 'midas_prepcma: RESUME type =', resumeType
  qcvar = (resumeType == 'POSTALT')
  if (qcvar) then 
    write(*,*) 'midas_prepcma: The input file is a postalt file'
  else
    write(*,*) 'midas_prepcma: The input file is NOT a postalt file'
  end if

  !- Initialize TOVS processing
  if (obs_famExist(obsSpaceData,'TO')) then
    call tvs_setup
    call checkTovsInstNamesInNml()
  end if

  !- Select the elements to assimilate and apply rejection flags
  if (suprep) call filt_suprep(obsSpaceData)

  !- Allocation for TOVS
  if (obs_famExist(obsSpaceData,'TO')) call tvs_setupAlloc(obsSpaceData)

  !- Initialize obs error covariances and set flag using 'util' column of stats_tovs
  call oer_setObsErrors(obsSpaceData, 'analysis', useTovsUtil_opt=applySatUtil) ! IN

  !- Call suprep again to 'black list' channels according to 'util' column of stats_tovs
  if (applySatUtil) call filt_suprep(obsSpaceData)

  call utl_tmg_stop(10)

  !- Setup timeCoord module and set dateStamp from env variable
  call tim_setup()
  if (tim_getDateStamp() == 0) then
    if (dateStampFromObs > 0) then
      ! use dateStamp from obs if not already set by env variable
      call tim_setDateStamp(dateStampFromObs)
    else
      call utl_abort('midas-prepcma: DateStamp was not set')
    end if
  end if

  !- Reject any observation outside the data assimilation window
  if (rejectOutsideTimeWindow) then
    call oti_setup( oti, obsSpaceData, numStep=1, &
                    headerIndexBeg=1, headerIndexEnd=obs_numheader(obsSpaceData), &
                    flagObsOutside_opt=.true. )
  end if
  
  !- Reject all IR radiance observation in arctic and antarctic (.i.e |lat|>60. )
  if (rejectHighLatIR) call enkf_rejectHighLatIR(obsSpaceData)

  !- Modify the obs error stddev for AMSUB in the tropics
  if (modifyAmsubObsError) call enkf_modifyAmsubObsError(obsSpaceData)

  !- Perform thinning for several observation types
  if (thinning) then
    ! perform thinning for aircraft observations
    if (thinningConv)     call thinning_fam(obsSpaceData, nai_pmax, nai_target, 'AI')
    ! perform thinning for satwind observations
    if (thinningConv)     call thinning_fam(obsSpaceData, nsw_pmax, nsw_target, 'SW')
    ! perform thinning for scatterometer observations
    if (thinningConv)     call thinning_fam(obsSpaceData, nsc_pmax, nsc_target, 'SC')
    ! perform thinning for radiance observations
    if (thinningRadiance) then
      ! thinning per instrument
      if (thinTovsPerInst) then
        call getTovsInstNameList(numTovsInstNameList,tovsInstNameList)
        write(*,*) 'midas-prepcma: numTovsInstNameList=', numTovsInstNameList, &
                    'tovsInstNameList(1:numTovsInstNameList)=', tovsInstNameList(1:numTovsInstNameList)

        do sensorIndex = 1, numTovsInstNameList
          tovsInstName = trim(tovsInstNameList(sensorIndex))

          maxNumHeaderPerInst = getMaxNumHeadersForTovsInst(tovsInstName)

          write(*,*) 'midas-prepcma: thinning for TO inst=',  trim(tovsInstName)

          call thinning_fam(obsSpaceData, nto_pmax, maxNumHeaderPerInst, 'TO', &
                            codtyp_opt=codtyp_get_codtyp(tovsInstName))
        end do

      else ! if (thinTovsPerInst)
        call thinning_fam(obsSpaceData, nto_pmax, nto_target, 'TO')
      end if ! thinTovsPerInst

    end if  ! thinningRadiance
  end if

  !- Write the results
  write(*,*)
  write(*,*) '> midas-prepcma: writing to files'

  !- Write to burp/sqlite files if requested
  if (writeObsFiles) then
    call obsf_writeFiles(obsSpaceData)
    if ( obsClean ) call obsf_cleanObsFiles()
  end if

  if (writeAsciiCmaFiles) then

    !- Remove all observations from obsSpaceData that will not be assimilated
    !- But, unlike the EnKF program, do not check value of OBS_ZHA
    if (obsClean) then
      call obs_clean(obsSpaceData, hx_dummy, 0, -1, qcvar, checkZha_opt=.false.)
    end if

    if (mmpi_nprocs > 1) then
      call obs_expandToMpiGlobal(obsSpaceData)
    end if

    if (mmpi_myid == 0) then
      !- Open file for ascii output
      nobsout = 0
      ierr = fnom(nobsout, obsout, 'FMT+SEQ+R/W', 0)
      call obs_print(obsSpaceData,nobsout)
      close(nobsout)

      !- Write the results in CMA format
      ncmahdr = 0
      ierr = fnom(ncmahdr, cmahdr, 'FTN+SEQ+UNF+R/W', 0)
      ncmabdy = 0
      ierr = fnom(ncmabdy, cmabdy, 'FTN+SEQ+UNF+R/W', 0)
      ncmadim = 0
      ierr = fnom(ncmadim, cmadim, 'FTN+SEQ+R/W', 0)
      ncmahx  = -1
      call obs_write(obsSpaceData, hx_dummy, 0, ncmahdr, ncmabdy, ncmahx, ncmadim)
      close(ncmahdr)
      close(ncmabdy)
      close(ncmadim)

      !! This used to contain a .true. or .false. value indicating if observations passed the QCVar
      !! Since, this is not the case, we can write .false.
      nbrpform = 0
      ierr = fnom(nbrpform, brpform, 'FTN+SEQ+R/W', 0)
      write(nbrpform,*) .false.
      close(nbrpform)

    end if

  end if

  !
  !- 3.  Ending
  !
  write(*,*)
  write(*,*) '> midas-prepcma: Ending'
  call obs_finalize(obsSpaceData) ! deallocate obsSpaceData

  call utl_printTime()
  call utl_tmg_stop(0)
  call tmg_terminate(mmpi_myid, 'TMG_INFO')

  call mmpi_finalize

  if ( mmpi_myid == 0 ) then
    call utl_writeStatus('PREPCMA_END')
  end if

contains

  subroutine thinning_fam(obsSpaceData, n_pmax, n_target, cfam, codtyp_opt)
    !
    ! :Purpose: thin the observations of the selected family
    !
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsSpaceData  ! the data in observation space
    real(8),           intent(in)    :: n_pmax(:)   ! pressure levels that separate vertical layers for the thinning
    integer,           intent(in)    :: n_target    ! maximum desired amount of data per 3-D box
    character(len=2),  intent(in)    :: cfam        ! family type
    integer, optional, intent(in)    :: codtyp_opt  ! optional supplied codtyp

    ! Locals:
    type(struct_reg) :: lsc
    integer :: idist, n_count_thin, iai, iseed(4)
    integer :: nobs_count, nobs_count_mpiGlobal, nobs_count_thin, nobs_count_thin_mpiGlobal
    integer :: nrep_count, nrep_count_mpiGlobal, nrep_count_thin, nrep_count_thin_mpiGlobal
    integer :: headerIndex, bodyIndex, bodyIndexBeg, bodyIndexEnd
    integer :: iblock, codtyp, ilat, incr, ipres, nblocksum, npres, nsize, numHeader
    logical :: count_obs, allRejected, beSilent
    real(4) :: lat_r4, lon_r4
    real(8) :: pressure, rannum
    real(8), allocatable :: latcenter(:), latmin(:), latmax(:), ranvals(:)
    integer, allocatable :: nblockoffset(:), nlonblock(:)
    integer, allocatable :: ai_indices(:,:), nstation(:,:), nstationMpiGlobal(:,:)
    real(8), allocatable :: keep_ai(:,:)
    integer :: sensorIndex
    integer, allocatable :: numHeaderPerTovsInstBeforeThin(:,:), numHeaderPerTovsInstAfterThin(:,:)
    integer, allocatable :: numHeaderPerTovsInstBeforeThin_mpiGlobal(:,:), numHeaderPerTovsInstAfterThin_mpiGlobal(:,:)
    logical :: thinTovsPerInst, headerFoundForInst
    character(len=codtyp_name_length) :: instName
    logical, save :: firstCall = .true.

    ! box size that is used for observation thinning 
    ! (the numerator is an approximate distance in km)
    real(8), parameter :: r0_count_km = 200.0/(2**0.5)
    ! next two parameters are not used in this program
    real(8), parameter :: r1_dum = 1.0
    real(8), parameter :: rz_dum = 1.0

    beSilent = (.not. firstCall) 
    if (firstCall) firstCall = .false.

    numHeader = obs_numheader(obsSpaceData)
    npres = size(n_pmax,1) 
    write(*,*) 'Start thinning for ', cfam, ' data'
    ! at this stage we still have many radiance channels 
    ! that will be rejected at a later stage.
    if (cfam == 'XX') then
      ! never getting here
      write(*,*) 'count individual observations'   
      count_obs = .true.
    else
      write(*,*) 'count the number of reports'
      count_obs = .false.
    end if

    call reg_init_struct(lsc, r0_count_km, r1_dum, rz_dum)
    if (mmpi_myid == 0) write(*,*) 'number of latitude bands: ', lsc%nlatband
    nsize = lsc%nlatband
    allocate(latmin(nsize))
    allocate(latmax(nsize))
    allocate(latcenter(nsize))
    allocate(nlonblock(nsize))
    allocate(nblockoffset(nsize))

    call reg_getlatitude(lsc%r0_rad, lsc%nlatband, latmin, latcenter, latmax, beSilent_opt=beSilent)
    if (mmpi_myid == 0) write(*,*) 'number of latitude bands: ',lsc%nlatband
    do ilat = 1, lsc%nlatband
      if (.not. beSilent .and. mmpi_myid == 0) write(*,*) ' band: ', ilat, ' latitude between ', latmin(ilat), latmax(ilat)
    end do
    call reg_getblock(lsc%nlatband, lsc%r0_rad, latmin, latmax, nlonblock)
    nblocksum = 0
    do ilat = 1, lsc%nlatband
      nblockoffset(ilat) = nblocksum
      nblocksum = nblocksum + nlonblock(ilat)
      if (.not. beSilent .and. mmpi_myid == 0) write(*,*) 'latband: ', ilat, ' no of blocks: ', nlonblock(ilat)
    end do 
    if (mmpi_myid == 0) write(*,*) 'total number of blocks: ', nblocksum

    nrep_count = 0
    nobs_count = 0

    if (cfam == 'TO') then
      allocate(numHeaderPerTovsInstBeforeThin(nblocksum,tvs_nsensors))
      allocate(numHeaderPerTovsInstBeforeThin_mpiGlobal(nblocksum,tvs_nsensors))
      allocate(numHeaderPerTovsInstAfterThin(nblocksum,tvs_nsensors))
      allocate(numHeaderPerTovsInstAfterThin_mpiGlobal(nblocksum,tvs_nsensors))
      numHeaderPerTovsInstBeforeThin(:,:) = 0
      numHeaderPerTovsInstBeforeThin_mpiGlobal(:,:) = 0
      numHeaderPerTovsInstAfterThin(:,:) = 0
      numHeaderPerTovsInstAfterThin_mpiGlobal(:,:) = 0
    end if

    allocate(nstation(nblocksum, npres))
    allocate(keep_ai(nblocksum, npres))

    nstation=0
    ! keep_ai = 1 corresponds to no thinning
    keep_ai(:,:) = 1.0

    allocate(ai_indices(numHeader,3))

    if (cfam /= 'TO') then
      if (present(codtyp_opt)) then
        call utl_abort('thining_fam: codtyp_opt argument only allowed for TO family')
      end if
    end if

    thinTovsPerInst = .false.
    instName = ''
    if (present(codtyp_opt)) then
      thinTovsPerInst = .true.
      instName = codtyp_get_name(codtyp_opt)
    end if

    header_loop: do headerIndex = 1, numHeader
      codtyp = obs_headElem_i(obsSpaceData, obs_ity, headerIndex)
      if ((cfam == 'AI' .and. (codtyp == 42  .or. codtyp == 128 .or. &
                               codtyp == 157 .or. codtyp == 177)) .or. &
          (cfam == 'SC' .and. codtyp == 254) .or. &
          (cfam == 'SW' .and. (codtyp == 88  .or. codtyp == 188)) .or. &
          (cfam == 'TO' .and. tvs_isIdBurpTovs(codtyp))) then

        ! skip this header if does not match the supplied codtyp(s)
        if (present(codtyp_opt)) then
          if (codtyp /= codtyp_opt) cycle header_loop
        end if

        ! skip this header if all observations already rejected
        bodyIndexBeg = obs_headElem_i(obsSpaceData, obs_rln, headerIndex)
        bodyIndexEnd = obs_headElem_i(obsSpaceData, obs_nlv, headerIndex) + bodyIndexBeg - 1
        allRejected = .true.
        body_loop: do bodyIndex = bodyIndexBeg, bodyIndexEnd
          if (obs_bodyElem_i(obsSpaceData, obs_ass, bodyIndex) == obs_assimilated) then
            allRejected = .false.
            exit body_loop
          end if
        end do body_loop
        if (allRejected) cycle header_loop

        lat_r4 = obs_headElem_r(obsSpaceData, obs_lat, headerIndex)
        lon_r4 = obs_headElem_r(obsSpaceData, obs_lon, headerIndex)
        call reg_locatestn(lsc%r0_rad, lat_r4, lon_r4, &
                           lsc%nlatband, nlonblock, &
                           nblockoffset, iblock)
        ! note that all data from the aircraft are at the same pressure
        if (npres == 1) then
          ipres = 1
        else
          pressure= obs_bodyElem_r(obsSpaceData, obs_ppp, bodyIndexBeg)
          ipres = 1
          do
            if ((n_pmax(ipres) > pressure) .or. (ipres > npres)) exit
            ipres = ipres + 1 
          end do
        end if  
        if (ipres <= npres) then
          if (count_obs) then
            incr = obs_headElem_i(obsSpaceData, obs_nlv, headerIndex)
          else
            incr = 1
          end if
          nstation(iblock,ipres) = nstation(iblock,ipres) + incr
          nobs_count = nobs_count + obs_headElem_i(obsSpaceData, obs_nlv, headerIndex)
          nrep_count = nrep_count + 1
          ai_indices(nrep_count, 1) = headerIndex
          ai_indices(nrep_count, 2) = iblock
          ai_indices(nrep_count, 3) = ipres 
        end if

        if (cfam=='TO' .and. tvs_isIdBurpTovs(obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex))) then
          sensorIndex = tvs_lsensor(headerIndex)
          numHeaderPerTovsInstBeforeThin(iblock,sensorIndex) = numHeaderPerTovsInstBeforeThin(iblock,sensorIndex) + 1
        end if

      end if
    end do header_loop

    ! do mpi communication of the accumulators
    allocate(nstationMpiGlobal(nblocksum, npres))
    call mmpi_allReduce(nstation, nstationMpiGlobal, 'mpi_sum')
    call mmpi_allReduce(nrep_count, nrep_count_mpiGlobal, 'mpi_sum')
    call mmpi_allReduce(nobs_count, nobs_count_mpiGlobal, 'mpi_sum')

    write(*,*) 'total number of ', cfam, ' reports (local and mpiglobal): ',  &
                nrep_count, nrep_count_mpiGlobal
    allocate(ranvals(nrep_count))
    write(*,*) 'total number of ', cfam, ' observations (local and mpiglobal): ',  &
                nobs_count, nobs_count_mpiGlobal

    if (cfam == 'TO') then
      call mmpi_allReduce(numHeaderPerTovsInstBeforeThin, numHeaderPerTovsInstBeforeThin_mpiGlobal, "mpi_sum")

      do sensorIndex = 1, tvs_nsensors
        write(*,*) 'total number of ', cfam, ' headers (local and mpiglobal) for ', &
                    inst_name(tvs_instruments(sensorIndex)), platform_name(tvs_platforms(sensorIndex)), tvs_satellites(sensorIndex), ':', &
                    sum(numHeaderPerTovsInstBeforeThin(:,sensorIndex)), sum(numHeaderPerTovsInstBeforeThin_mpiGlobal(:,sensorIndex))
      end do

      call printToListingsForTovs(instName, nblocksum, numHeaderPerTovsInstBeforeThin, &
                                  numHeaderPerTovsInstBeforeThin_mpiGlobal, thinTovsPerInst, &
                                  headerFoundForInst)

      if (.not. headerFoundForInst) return

    end if ! cfam=='TO'

    n_count_thin = 0
    do iblock = 1, nblocksum
      do ipres = 1, npres
        if (nstationMpiGlobal(iblock,ipres) >= 1) then
          if (.not. beSilent .and. mmpi_myid == 0) then
            write(*,*) 'block ipres and count: ',iblock, ipres, nstationMpiGlobal(iblock,ipres)
          end if

          if (nstationMpiGlobal(iblock,ipres) > n_target) then
            keep_ai(iblock,ipres) = dble(n_target) / dble(nstationMpiGlobal(iblock,ipres))
            n_count_thin = n_count_thin + n_target
          else
            n_count_thin = n_count_thin + nstationMpiGlobal(iblock,ipres)
          end if
        end if
      end do
    end do

    if (count_obs) then
      write(*,*) 'Estimated remaining number of ', cfam, ' observations (mpiGlobal): ', n_count_thin
    else 
      write(*,*) 'Estimated remaining number of ', cfam, ' reports (mpiGlobal): ', n_count_thin
    end if

    nrep_count_thin = 0
    nobs_count_thin = 0

    idist = 1
    iseed(1) = 1
    iseed(2) = 5
    iseed(3) = 9
    iseed(4) = 11
    call dlarnv(idist,iseed,nrep_count,ranvals)
    do iai = 1, nrep_count
      headerIndex = ai_indices(iai,1)
      iblock = ai_indices(iai,2)
      ipres  = ai_indices(iai,3)
      rannum = ranvals(iai)
      if (rannum <= keep_ai(iblock,ipres)) then
        nrep_count_thin = nrep_count_thin + 1
        nobs_count_thin = nobs_count_thin + obs_headElem_i(obsSpaceData, obs_nlv, headerIndex)

        if (cfam=='TO' .and. tvs_isIdBurpTovs(obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex))) then
          sensorIndex = tvs_lsensor(headerIndex)
          numHeaderPerTovsInstAfterThin(iblock,sensorIndex) = numHeaderPerTovsInstAfterThin(iblock,sensorIndex) + 1
        end if
      else
        ! reject the profile
        bodyIndexBeg = obs_headElem_i(obsSpaceData, obs_rln, headerIndex)
        bodyIndexEnd = obs_headElem_i(obsSpaceData, obs_nlv, headerIndex) + bodyIndexBeg - 1
        do bodyIndex = bodyIndexBeg, bodyIndexEnd
          call obs_bodySet_i(obsSpaceData, obs_ass, bodyIndex, obs_notAssimilated)
          ! also set the 'rejected by selection process' flag (bit 11)
          call obs_bodySet_i( obsSpaceData, obs_flg, bodyIndex,  &
                              ibset( obs_bodyElem_i( obsSpaceData, obs_flg, bodyIndex ), 11) )
        end do       
      end if
    end do

    ! mpi communication of accumulators
    call mmpi_allReduce(nrep_count_thin, nrep_count_thin_mpiGlobal, 'mpi_sum')
    call mmpi_allReduce(nobs_count_thin, nobs_count_thin_mpiGlobal, 'mpi_sum')
    
    write(*,*) 'True remaining number of ', cfam, ' reports (local, mpiGlobal): ',  &
          nrep_count_thin, nrep_count_thin_mpiGlobal
    write(*,*) 'True remaining number of ', cfam, ' observations (local, mpiGlobal): ',  &
          nobs_count_thin, nobs_count_thin_mpiGlobal

    if (cfam == 'TO') then
      call mmpi_allReduce(numHeaderPerTovsInstAfterThin, numHeaderPerTovsInstAfterThin_mpiGlobal, "mpi_sum")

      do sensorIndex = 1, tvs_nsensors
        write(*,*) 'True remaining number of ', cfam, ' headers (local and mpiglobal) for ', &
                    inst_name(tvs_instruments(sensorIndex)), platform_name(tvs_platforms(sensorIndex)), tvs_satellites(sensorIndex), ':', &
                    sum(numHeaderPerTovsInstAfterThin(:,sensorIndex)), sum(numHeaderPerTovsInstAfterThin_mpiGlobal(:,sensorIndex))
      end do

      call printToListingsForTovs(instName, nblocksum, numHeaderPerTovsInstAfterThin, &
                                  numHeaderPerTovsInstAfterThin_mpiGlobal, thinTovsPerInst, &
                                  headerFoundForInst, ifAfterThinning_opt=.true.)

      if (.not. headerFoundForInst) then
        call utl_abort('thinning_fam: should have returned before thinning was performed')
      end if

    end if ! cfam=='TO'

    deallocate(ranvals)
    deallocate(latmin)
    deallocate(latmax)
    deallocate(latcenter)
    deallocate(nlonblock)
    deallocate(nblockoffset)
    deallocate(nstation)
    deallocate(nstationMpiGlobal)
    deallocate(keep_ai)
    deallocate(ai_indices)

    if (cfam == 'TO') then
      deallocate(numHeaderPerTovsInstBeforeThin)
      deallocate(numHeaderPerTovsInstBeforeThin_mpiGlobal)
      deallocate(numHeaderPerTovsInstAfterThin)
      deallocate(numHeaderPerTovsInstAfterThin_mpiGlobal)
    end if

  end subroutine thinning_fam

  !----------------------------------------------------------------------
  ! getTovsInstNameList
  !----------------------------------------------------------------------
  subroutine getTovsInstNameList(numInstNameUniqueList, instNameUniqueList)
    !
    ! :Purpose: Create a unique list of tovs instrument names.
    !
    implicit none

    ! Arguments:
    integer,                       intent(out)   :: numInstNameUniqueList ! Number of tovs instrument names in the list
    character(len=*), allocatable, intent(inout) :: instNameUniqueList(:) ! Unique list of tovs instrument names

    ! Locals:
    integer :: sensorIndex, sensorIndex2
    character(len=codtyp_name_length) :: instNameList(tvs_nsensors)

    ! replace all amsub with mhs in the original list
    do sensorIndex = 1, tvs_nsensors
      if (trim(inst_name(tvs_instruments(sensorIndex))) == 'amsub') then
        instNameList(sensorIndex) = 'mhs'
      else
        instNameList(sensorIndex) = trim(inst_name(tvs_instruments(sensorIndex)))
      end if
    end do

    if (mmpi_myid == 0) then
      write(*,*) 'getTovsInstNameList: all occurances of amsub are replaced by mhs in instNameUniqueList'
    end if

    allocate(instNameUniqueList(tvs_nsensors))
    instNameUniqueList(:) = ''
    numInstNameUniqueList = 1
    instNameUniqueList(numInstNameUniqueList) = trim(instNameList(1))

    loopSensor3: do sensorIndex = 1, tvs_nsensors
      do sensorIndex2 = 1, numInstNameUniqueList
        if (trim(instNameUniqueList(sensorIndex2)) == trim(instNameList(sensorIndex))) then
          cycle loopSensor3
        end if
      end do

      numInstNameUniqueList = numInstNameUniqueList + 1
      instNameUniqueList(numInstNameUniqueList) = trim(instNameList(sensorIndex))
    end do loopSensor3

  end subroutine getTovsInstNameList

  !----------------------------------------------------------------------
  ! checkTovsInstNamesInNml
  !----------------------------------------------------------------------
  subroutine checkTovsInstNamesInNml()
    !
    ! :Purpose: Perform the following checks on namelist variable tovsInstNamesWithMaxNumHeaders: 
    !           1- no duplicate in the list; 2- amsub is not present in the list. 
    !
    implicit none

    ! Locals:
    integer :: sensorIndex, sensorIndex2

    ! check no duplicate in instrument names in nml
    loopSensor4: do sensorIndex = 1, tvs_nsensors
      do sensorIndex2 = sensorIndex + 1, tvs_nsensors
        if (trim(tovsInstNamesWithMaxNumHeaders(sensorIndex)) /= 'NOT_DEFINED') then
          if (trim(tovsInstNamesWithMaxNumHeaders(sensorIndex)) == trim(tovsInstNamesWithMaxNumHeaders(sensorIndex2))) then
            write(*,*) trim(tovsInstNamesWithMaxNumHeaders(sensorIndex))
            call utl_abort('checkTovsInstNamesInNml: duplicate instrument names exist in tovsInstNamesWithMaxNumHeaders')
          end if
        else
          exit loopSensor4
        end if
      end do
    end do loopSensor4

    do sensorIndex = 1, tvs_nsensors
      if (trim(tovsInstNamesWithMaxNumHeaders(sensorIndex)) == 'amsub') then
        call utl_abort('checkTovsInstNamesInNml: amsub exist in tovsInstNamesWithMaxNumHeaders, replace with mhs.')        
      end if
    end do

  end subroutine checkTovsInstNamesInNml

  !----------------------------------------------------------------------
  ! getMaxNumHeadersForTovsInst
  !----------------------------------------------------------------------
  function getMaxNumHeadersForTovsInst(tovsInstName) result(maxValue)
    !
    ! :Purpose: Return max number of headers for tovs instrument.
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: tovsInstName ! tovs instrument name
    ! Result:
    integer                      :: maxValue ! max value of the headers in each region for tovs instrument name

    ! Locals:
    integer :: instrumentIndex
    logical :: instFound

    instFound = .false.
    inst_loop: do instrumentIndex = 1, tvs_nsensors
      if (trim(tovsInstNamesWithMaxNumHeaders(instrumentIndex)) == trim(tovsInstName)) then
        maxValue = maxNumHeadersForTovsInst(instrumentIndex)
        instFound = .true.

        exit inst_loop
      end if
    end do inst_loop

    if (instFound) then
      if (mmpi_myid == 0) write(*,*) 'for instrument, max number of headers:', trim(tovsInstName), maxValue
    else
      maxValue = nto_target
      if (mmpi_myid == 0) write(*,*) 'for instrument, using default max number of headers:', trim(tovsInstName)
    end if

  end function getMaxNumHeadersForTovsInst

  !----------------------------------------------------------------------
  ! printToListingsForTovs
  !----------------------------------------------------------------------
  subroutine printToListingsForTovs(instName, nblocksum, numHeaderPerTovsInst, &
                                    numHeaderPerTovsInst_mpiGlobal, thinTovsPerInst, &
                                    headerFoundForInst, ifAfterThinning_opt)
    !
    ! :Purpose: Print to the listings for TOVS instrument.
    !
    implicit none

    ! Arguments:
    character(len=codtyp_name_length), intent(in) :: instName                  ! instrument name
    integer,                           intent(in) :: nblocksum                 ! sum of number of regions over all TOVS
    integer,                           intent(in) :: numHeaderPerTovsInst(:,:) ! local number of headers per instrument
    integer,                           intent(in) :: numHeaderPerTovsInst_mpiGlobal(:,:) ! global number of headers per instrument
    logical,                           intent(in) :: thinTovsPerInst           ! do thinning for tovs per instrument independently
    logical,                           intent(out) :: headerFoundForInst       ! if headers found for instrument
    logical, optional,                 intent(in) :: ifAfterThinning_opt       ! if subroutine called after thinning is performed

    ! Locals:
    integer :: sensorIndex, sensorIndex2, numInstNameUniqueListWithHeader, numMatchFound, matchFoundIndex
    integer :: iblock, numHeadersFound, numHeadersFound_mpiGlobal
    integer, allocatable :: matchIndexList(:), numHeadersFoundInBlock(:,:), numHeadersFoundInBlock_mpiGlobal(:,:)
    logical :: ifAfterThinning
    character(len=codtyp_name_length) :: instNameUniqueListWithHeader(tvs_nsensors)

    if (present(ifAfterThinning_opt)) then
      ifAfterThinning = ifAfterThinning_opt
    else
      ifAfterThinning = .false.
    end if

    headerFoundForInst = .true.
 
    ! find the first inst with non-zero number of headers
    instNameUniqueListWithHeader(:) = ''
    numInstNameUniqueListWithHeader = 1
    loopSensor1: do sensorIndex = 1, tvs_nsensors
      if (sum(numHeaderPerTovsInst_mpiGlobal(:,sensorIndex)) > 0) then
        instNameUniqueListWithHeader(numInstNameUniqueListWithHeader) = trim(inst_name(tvs_instruments(sensorIndex)))
        exit loopSensor1
      end if
    end do loopSensor1

    ! find rest of the inst with non-zero number of headers
    loopSensor2: do sensorIndex = 1, tvs_nsensors
      do sensorIndex2 = 1, numInstNameUniqueListWithHeader
        if (trim(instNameUniqueListWithHeader(sensorIndex2)) == trim(inst_name(tvs_instruments(sensorIndex))) .or. &
            sum(numHeaderPerTovsInst_mpiGlobal(:,sensorIndex)) == 0) then
          cycle loopSensor2
        end if
      end do

      numInstNameUniqueListWithHeader = numInstNameUniqueListWithHeader + 1
      instNameUniqueListWithHeader(numInstNameUniqueListWithHeader) = trim(inst_name(tvs_instruments(sensorIndex)))
    end do loopSensor2

    if (mmpi_myid == 0) then
      write(*,*) 'numInstNameUniqueListWithHeader=', numInstNameUniqueListWithHeader, &
                 ', instNameUniqueListWithHeader(:)=', instNameUniqueListWithHeader(:)
    end if

    do sensorIndex2 = 1, numInstNameUniqueListWithHeader
      matchIndexList = utl_findlocs(inst_name(tvs_instruments(:)),instNameUniqueListWithHeader(sensorIndex2))
      if (matchIndexList(1) > 0) then
        numMatchFound = size(matchIndexList)
        numHeadersFound = 0
        numHeadersFound_mpiGlobal = 0
        do matchFoundIndex = 1, numMatchFound
          numHeadersFound = numHeadersFound + sum(numHeaderPerTovsInst(:,matchIndexList(matchFoundIndex)))
          numHeadersFound_mpiGlobal = numHeadersFound_mpiGlobal + sum(numHeaderPerTovsInst_mpiGlobal(:,matchIndexList(matchFoundIndex)))
        end do

        if (mmpi_myid == 0 .and. .not. ifAfterThinning) then
          write(*,*) 'total number of TO headers (local and mpiglobal) for all ', &
                      trim(instNameUniqueListWithHeader(sensorIndex2)), ':', &
                      numHeadersFound, numHeadersFound_mpiGlobal
        end if

        if (thinTovsPerInst) then
          if (numInstNameUniqueListWithHeader == 1) then
            if (trim(instNameUniqueListWithHeader(1)) /= trim(instName)) then
              write(*,*) 'instName=', instName
              call utl_abort('printToListingsForTovs: instNameUniqueListWithHeader(1)) /= instName')
            end if
          else
            call utl_abort('printToListingsForTovs: numInstNameUniqueListWithHeader > 1')
          end if
        end if ! thinTovsPerInst

      else
        ! return because there is no header for this instrument
        if (thinTovsPerInst) then
          write(*,*) 'no headers found for ', trim(instNameUniqueListWithHeader(sensorIndex2)), ', returning ...'
          headerFoundForInst = .false.
          return
        end if

      end if ! matchIndexList(1) > 0
    end do ! sensorIndex2

    ! compute counts per block per sensor
    allocate(numHeadersFoundInBlock(nblocksum,numInstNameUniqueListWithHeader))
    allocate(numHeadersFoundInBlock_mpiGlobal(nblocksum,numInstNameUniqueListWithHeader))
    numHeadersFoundInBlock(:,:) = 0
    numHeadersFoundInBlock_mpiGlobal(:,:) = 0
    do sensorIndex2 = 1, numInstNameUniqueListWithHeader
      matchIndexList = utl_findlocs(inst_name(tvs_instruments(:)),instNameUniqueListWithHeader(sensorIndex2))
      if (matchIndexList(1) > 0) then
        numMatchFound = size(matchIndexList)

        do matchFoundIndex = 1, numMatchFound
          numHeadersFoundInBlock(:,sensorIndex2) = numHeadersFoundInBlock(:,sensorIndex2) + &
                                                   numHeaderPerTovsInst(:,matchIndexList(matchFoundIndex))
          numHeadersFoundInBlock_mpiGlobal(:,sensorIndex2) = numHeadersFoundInBlock_mpiGlobal(:,sensorIndex2) + &
                                                             numHeaderPerTovsInst_mpiGlobal(:,matchIndexList(matchFoundIndex))
        end do
      end if
    end do ! sensorIndex2

    ! check the sum over all blocks match the counts per sensor 
    do sensorIndex2 = 1, numInstNameUniqueListWithHeader
      matchIndexList = utl_findlocs(inst_name(tvs_instruments(:)),instNameUniqueListWithHeader(sensorIndex2))
      if (matchIndexList(1) > 0) then
        numMatchFound = size(matchIndexList)

        numHeadersFound_mpiGlobal = 0
        do matchFoundIndex = 1, numMatchFound
          numHeadersFound_mpiGlobal = numHeadersFound_mpiGlobal + &
                                      sum(numHeaderPerTovsInst_mpiGlobal(:,matchIndexList(matchFoundIndex)))        
        end do

        if (sum(numHeadersFoundInBlock_mpiGlobal(:,sensorIndex2)) /= numHeadersFound_mpiGlobal) then
          write(*,*) 'sensorIndex2=', sensorIndex2, ', instName=', trim(instNameUniqueListWithHeader(sensorIndex2)), &
                      ', numHeadersFound_mpiGlobal=', numHeadersFound_mpiGlobal, &
                      ', sum=', sum(numHeadersFoundInBlock_mpiGlobal(:,sensorIndex2))

          if (.not. ifAfterThinning) then
            call utl_abort('printToListingsForTovs: the sums do not match before thinning')
          else
            call utl_abort('printToListingsForTovs: the sums do not match after thinning')
          end if

        end if
      end if
    end do

    ! print counts per block per sensor
    if (mmpi_myid == 0) then
      if (thinTovsPerInst) then
        if (.not. ifAfterThinning) then
          write(*,'(a,1x,a)') 'before thinning global one instr:', trim(instNameUniqueListWithHeader(1))
        else
          write(*,'(a,1x,a)') 'after thinning global one instr:', trim(instNameUniqueListWithHeader(1))
        end if

        do iblock = 1, nblocksum
          if (.not. ifAfterThinning) then
            write(*,'(a,1x,a,i6,a,i8)') 'before thinning global', trim(instNameUniqueListWithHeader(1)), &
                                      iblock, ':', numHeadersFoundInBlock_mpiGlobal(iblock,1)
          else
            write(*,'(a,1x,a,i6,a,i8)') 'after thinning global', trim(instNameUniqueListWithHeader(1)), &
                                      iblock, ':', numHeadersFoundInBlock_mpiGlobal(iblock,1)
          end if
        end do ! iblock
      else
        if (.not. ifAfterThinning) then
          write(*,'(a39)',advance='no') 'before thinning global all instr       '
        else
          write(*,'(a39)',advance='no') 'after thinning global all instr        '
        end if
        do sensorIndex2 = 1, numInstNameUniqueListWithHeader
          write(*,'(a8,1x)',advance='no') trim(instNameUniqueListWithHeader(sensorIndex2))
        end do
        write(*,*)

        do iblock = 1, nblocksum
          if (.not. ifAfterThinning) then
            write(*,'(a32,i6,a)',advance='no') 'before thinning global all instr', iblock, ':'
          else
            write(*,'(a32,i6,a)',advance='no') 'after thinning global all instr ', iblock, ':'
          end if

          do sensorIndex2 = 1, numInstNameUniqueListWithHeader
            write(*,'(i8,1x)',advance='no') numHeadersFoundInBlock_mpiGlobal(iblock,sensorIndex2)
          end do
          write(*,*)
          
        end do ! iblock

        ! print max of all block for sensor
        if (.not. ifAfterThinning) then
          write(*,'(a39)',advance='no') 'before thinning global all instr MAX   '
        else
          write(*,'(a39)',advance='no') 'after thinning global all instr MAX    '
        end if

        do sensorIndex2 = 1, numInstNameUniqueListWithHeader
          write(*,'(i8,1x)',advance='no') maxval(numHeadersFoundInBlock_mpiGlobal(:,sensorIndex2))
        end do
        write(*,*)
      end if ! thinTovsPerInst
    end if ! mmpi_myid

    deallocate(numHeadersFoundInBlock)
    deallocate(numHeadersFoundInBlock_mpiGlobal)

  end subroutine printToListingsForTovs

end program
