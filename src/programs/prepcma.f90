!-------------------------------------- LICENCE BEGIN ------------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------

program midas_prepcma
  !
  ! :Purpose: Read observation files that are in bgckalt format (i.e. output from
  !           the background check) and transform them to the ObsSpaceData
  !           format. All observations will collected in a single ObsSpaceData
  !           structure. The structure is output in binary format.
  !
  use version_mod
  use obsSpaceData_mod
  use obsFiles_mod
  use obsFilter_mod
  use obsTimeInterp_mod
  use obsErrors_mod
  use tovs_nl_mod
  use timeCoord_mod
  use enkf_mod
  use utilities_mod
  use midasMpi_mod
  use ramDisk_mod
  use regions_mod
  use burpRead_mod

  implicit none

  integer :: fnom, fclos, nulnam, ierr, dateStampFromObs
  type(struct_obs), target  :: obsSpaceData
  type(struct_oti), pointer :: oti => null()
  real(kind=8) :: hx_dummy(1,1)
  integer :: ncmahdr, ncmahx, ncmabdy, ncmadim, nobsout, nbrpform
  logical :: qcvar, numHeader, numBody
  character(len=7) :: resumeType

  ! number of pressure ranges used for the thinning of aircraft (and other) data:
  integer, parameter :: npres_ai = 5
  integer, parameter :: npres_sw = 2
  integer, parameter :: nai_target = 10
  integer, parameter :: nsc_target = 10
  integer, parameter :: nsw_target = 6
  integer, parameter :: nto_target = 6 
  real(8) :: nai_pmax(npres_ai) = (/ 25000.0, 40000.0, 60000.0, 80000.0, 110000.0/)
  real(8) :: nsw_pmax(npres_sw) = (/ 60000.0, 110000.0/)
  ! For a scalar array, no layer selection will be done
  real(8) :: nsc_pmax(1) = (/ 0.0 /)
  real(8) :: nto_pmax(1) = (/ 0.0 /)

  ! Namelist variables:
  character(len=256) :: cmahdr
  character(len=256) :: cmabdy
  character(len=256) :: cmadim
  character(len=256) :: obsout
  character(len=256) :: brpform
  logical :: suprep
  logical :: rejectOutsideTimeWindow
  logical :: thinning
  logical :: applySatUtil
  logical :: modifyAmsubObsError
  logical :: rejectHighLatIR
  logical :: obsClean
  logical :: writeObsFiles
  logical :: writeAsciiCmaFiles

  NAMELIST /NAMPREPCMA/ cmahdr, cmabdy, cmadim, obsout, brpform,  &
                        suprep, rejectOutsideTimeWindow, thinning, &
                        applySatUtil, modifyAmsubObsError, rejectHighLatIR, &
                        obsClean, writeObsFiles, writeAsciiCmaFiles

  call ver_printNameAndVersion('prepcma','Prepare observations for LETKF')

  !- 1.0 mpi
  call mmpi_initialize

  !- 1.1 timings
  call tmg_init(mmpi_myid, 'TMG_INFO')
  call utl_tmg_start(0,'Main')

  if ( mmpi_myid == 0 ) call utl_writeStatus('PREPCMA_BEG')

  !- Specify default values for namelist variables
  cmahdr        = 'NOT_DEFINED'
  cmabdy        = 'NOT_DEFINED'
  cmadim        = 'NOT_DEFINED'
  obsout        = 'NOT_DEFINED'
  suprep                  = .true.
  rejectOutsideTimeWindow = .true.
  thinning                = .true.
  applySatUtil            = .true.
  modifyAmsubObsError     = .true.
  rejectHighLatIR         = .true.
  obsClean                = .true.
  writeObsFiles           = .false.
  writeAsciiCmaFiles       = .false.

  nulnam = 0
  ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read(nulnam,nml=namprepcma,iostat=ierr)
  if (ierr /= 0) call utl_abort('midas-prepcma: Error reading namelist')
  if (mmpi_myid == 0) write(*,nml=namprepcma)
  ierr = fclos(nulnam)

  !- RAM disk usage
  call ram_setup

  call utl_tmg_start(10,'--Observations')

  !- Set up list of elements to be assimilated and flags for rejection (from namelist)
  call filt_setup('prepcma')

  !- Observation file names and get datestamp
  call obsf_setup(dateStampFromObs, 'prepcma' )

  !- Allocate obsSpaceData
  call obs_class_initialize('ENKFMIDAS')
  call obs_initialize( obsSpaceData, mpi_local=obsf_filesSplit() )

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
  if (obs_famExist(obsSpaceData,'TO')) call tvs_setup

  !- Select the elements to assimilate and apply rejection flags
  if (suprep) call filt_suprep(obsSpaceData)

  !- Allocation for TOVS
  if (obs_famExist(obsSpaceData,'TO')) call tvs_setupAlloc(obsSpaceData)

  !- Initialize obs error covariances and set flag using 'util' column of stats_tovs
  call oer_setObsErrors(obsSpaceData, 'analysis', useTovsUtil_opt=applySatUtil) ! IN

  !- Call suprep again to 'black list' channels according to 'util' column of stats_tovs
  if (applySatUtil) call filt_suprep(obsSpaceData)

  call utl_tmg_stop(10)

  !- Setup timeCoord module
  call tim_setup()
  call tim_setDateStamp(dateStampFromObs)

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
    call thinning_fam(obsSpaceData, nai_pmax, nai_target, 'AI')
    ! perform thinning for scatterometer observations
    call thinning_fam(obsSpaceData, nsc_pmax, nsc_target, 'SC')
    ! perform thinning for radiance observations
    call thinning_fam(obsSpaceData, nto_pmax, nto_target, 'TO')
    ! perform thinning for satwind observations
    call thinning_fam(obsSpaceData, nsw_pmax, nsw_target, 'SW')
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

  call utl_tmg_stop(0)
  call tmg_terminate(mmpi_myid, 'TMG_INFO')

  call rpn_comm_finalize(ierr)

  if ( mmpi_myid == 0 ) then
    call utl_writeStatus('PREPCMA_END')
  end if

contains

  subroutine thinning_fam(obsSpaceData, n_pmax, n_target, cfam)
    !
    ! :Purpose: thin the observations of the selected family
    !
    implicit none

    ! arguments:
    type (struct_obs) :: obsSpaceData  ! the data in observation space
    real(8), intent(in) :: n_pmax(:)   ! pressure levels that separate vertical layers for the thinning
    integer, intent(in) :: n_target    ! maximum desired amount of data per 3-D box
    character(len=2) :: cfam           ! family type

    ! locals:
    type(struct_reg) :: lsc
    integer :: idist, n_count_thin, iai, iseed(4)
    integer :: nobs_count, nobs_count_mpiGlobal, nobs_count_thin, nobs_count_thin_mpiGlobal
    integer :: nrep_count, nrep_count_mpiGlobal, nrep_count_thin, nrep_count_thin_mpiGlobal
    integer :: headerIndex, bodyIndex, bodyIndexBeg, bodyIndexEnd
    integer :: iblock, codeType, ilat, incr, ipres, nblocksum, npres, nsize, num_stn
    logical :: count_obs, allRejected
    real(4) :: lat_r4, lon_r4
    real(8) :: pressure, rannum
    real(8), allocatable :: latcenter(:), latmin(:), latmax(:), ranvals(:)
    integer, allocatable :: nblockoffset(:), nlonblock(:)
    integer, allocatable :: ai_indices(:,:), nstation(:,:), nstationMpiGlobal(:,:)
    real(8), allocatable :: keep_ai(:,:)

    ! box size that is used for observation thinning 
    ! (the numerator is an approximate distance in km)
    real(8), parameter :: r0_count_km = 200.0/(2**0.5)
    ! next two parameters are not used in this program
    real(8), parameter :: r1_dum = 1.0
    real(8), parameter :: rz_dum = 1.0
   
    num_stn = obs_numheader(obsSpaceData)
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

    call reg_getlatitude(lsc%r0_rad, lsc%nlatband, latmin, latcenter, latmax)
    if (mmpi_myid == 0) write(*,*) 'number of latitude bands: ',lsc%nlatband
    do ilat = 1, lsc%nlatband
      if (mmpi_myid == 0) write(*,*) ' band: ', ilat, ' latitude between ', latmin(ilat), latmax(ilat)
    end do
    call reg_getblock(lsc%nlatband, lsc%r0_rad, latmin, latmax, nlonblock)
    nblocksum = 0
    do ilat = 1, lsc%nlatband
      nblockoffset(ilat) = nblocksum
      nblocksum = nblocksum + nlonblock(ilat)
      if (mmpi_myid == 0) write(*,*) 'latband: ', ilat, ' no of blocks: ', nlonblock(ilat)
    end do 
    if (mmpi_myid == 0) write(*,*) 'total number of blocks: ', nblocksum

    nrep_count = 0
    nobs_count = 0

    allocate(nstation(nblocksum, npres))
    allocate(keep_ai(nblocksum, npres))
 
    nstation=0
    ! keep_ai = 1 corresponds to no thinning
    keep_ai(:,:) = 1.0

    allocate(ai_indices(num_stn,3))

    header_loop: do headerIndex = 1, num_stn
      codeType= obs_headElem_i(obsSpaceData, obs_ity, headerIndex)
      if ( (cfam=='AI' .and. (codeType==42  .or. codeType==128 .or. &
                              codeType==157 .or. codeType==177)) .or. &
           (cfam=='SC' .and. codeType==254) .or. &
           (cfam=='SW' .and. (codeType==88  .or. codeType==188)) .or. &
           (cfam=='TO' .and. tvs_isIdBurpTovs(codeType)) ) then

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
      end if
    end do header_loop

    ! do mpi communication of the accumulators
    nsize = nblocksum * npres
    allocate(nstationMpiGlobal(nblocksum, npres))
    call rpn_comm_allreduce(nstation, nstationMpiGlobal, nsize,  &
                            'mpi_integer','mpi_sum', 'GRID', ierr)
    call rpn_comm_allreduce(nrep_count, nrep_count_mpiGlobal, 1,  &
                            'mpi_integer','mpi_sum', 'GRID', ierr)
    call rpn_comm_allreduce(nobs_count, nobs_count_mpiGlobal, 1,  &
                            'mpi_integer','mpi_sum', 'GRID', ierr)

    write(*,*) 'total number of ', cfam, ' reports (local and mpiglobal): ',  &
         nrep_count, nrep_count_mpiGlobal
    allocate(ranvals(nrep_count))
    write(*,*) 'total number of ', cfam, ' observations (local and mpiglobal): ',  &
         nobs_count, nobs_count_mpiGlobal

    n_count_thin = 0
    do iblock = 1, nblocksum
      do ipres = 1, npres
        if (nstationMpiGlobal(iblock,ipres) .ge. 1) then
          if (mmpi_myid == 0) write(*,*) 'block ipres and count: ',iblock,ipres, &
                     nstationMpiGlobal(iblock,ipres)
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
    call rpn_comm_allreduce(nrep_count_thin, nrep_count_thin_mpiGlobal, 1,  &
                            'mpi_integer','mpi_sum', 'GRID', ierr)
    call rpn_comm_allreduce(nobs_count_thin, nobs_count_thin_mpiGlobal, 1,  &
                            'mpi_integer','mpi_sum', 'GRID', ierr)
    
    write(*,*) 'True remaining number of ', cfam, ' reports (local, mpiGlobal): ',  &
         nrep_count_thin, nrep_count_thin_mpiGlobal
    write(*,*) 'True remaining number of ', cfam, ' observations (local, mpiGlobal): ',  &
         nobs_count_thin, nobs_count_thin_mpiGlobal

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

  end subroutine thinning_fam

end program
