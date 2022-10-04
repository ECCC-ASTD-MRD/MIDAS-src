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

module obsErrors_mod
  ! MODULE obsErrors_mod (prefix='oer' category='2. B and R matrices')
  !
  ! :Purpose: Subroutines to set up the observation-error standard deviations.
  !
  use midasMpi_mod
  use mathPhysConstants_mod
  use obsSpaceData_mod
  use obsSubSpaceData_mod
  use tovs_nl_mod
  use codtyp_mod
  use bufr_mod
  use utilities_mod
  use earthConstants_mod
  use gps_mod
  use columnData_mod
  use rmatrix_mod
  use varNameList_mod
  use obsfiles_mod
  use burp_module
  use rttov_const, only: surftype_sea
  implicit none
  save
  private

  ! public procedures
  public :: oer_setObsErrors, oer_SETERRGPSGB, oer_SETERRGPSRO, oer_setErrBackScatAnisIce, oer_sw
  public :: oer_setInterchanCorr, oer_computeInflatedStateDepSigmaObs
  
  ! public functions
  public :: oer_getSSTdataParam_char, oer_getSSTdataParam_int, oer_getSSTdataParam_R8

  ! public variables (parameters)
  public :: oer_ascatAnisOpenWater, oer_ascatAnisIce
  
 ! Temporary arrays for QC purpose
  public :: oer_toverrst, oer_clwThreshArr, oer_tovutil
  public :: oer_sigmaObsErr, oer_useStateDepSigmaObs 
  ! TOVS OBS ERRORS
  real(8) :: toverrst(tvs_maxChannelNumber,tvs_maxNumberOfSensors)
  real(8) :: clwThreshArr(tvs_maxChannelNumber,tvs_maxNumberOfSensors,2)
  real(8) :: sigmaObsErr(tvs_maxChannelNumber,tvs_maxNumberOfSensors,2)
  integer :: tovutil(tvs_maxChannelNumber,tvs_maxNumberOfSensors)
  logical :: useStateDepSigmaObs(tvs_maxChannelNumber,tvs_maxNumberOfSensors)
 ! Temporary arrays for QC purpose
  real(8) :: oer_toverrst(tvs_maxChannelNumber,tvs_maxNumberOfSensors)
  real(8) :: oer_clwThreshArr(tvs_maxChannelNumber,tvs_maxNumberOfSensors,2)
  real(8) :: oer_sigmaObsErr(tvs_maxChannelNumber,tvs_maxNumberOfSensors,2)
  real(8) :: clearClwThresholdSigmaObsInflation(tvs_maxChannelNumber,tvs_maxNumberOfSensors)
  real(8) :: stateDepSigmaObsInflationCoeff(tvs_maxNumberOfSensors)
  integer :: oer_tovutil(tvs_maxChannelNumber,tvs_maxNumberOfSensors)
  logical :: oer_useStateDepSigmaObs(tvs_maxChannelNumber,tvs_maxNumberOfSensors)

  ! SST data 
  integer, parameter :: maxNumberSSTDatasets = 15
  integer :: numberSSTDatasets = 0 ! number of SST datasets in namelist
  type SSTdataParamsType
    character(len=20) :: dataType   = '' ! type of data: insitu, satellite, pseudo
    character(len=20) :: instrument = '' ! instrument: drifts, bouys, ships, AVHRR, VIIRS, AMSR2
    character(len=20) :: sensor     = '' ! sensor of satellite data: NOAA19, NOAA20,...
    character(len=20) :: sensorType = '' ! type of satellite data sensors: infrared, microwave,..
    integer           :: codeType   = MPC_missingValue_INT ! data codtype
    real(8)           :: dayError   = MPC_missingValue_R8  ! data error for daytime 
    real(8)           :: nightError = MPC_missingValue_R8  ! data error for nighttime
  end type SSTdataParamsType
  type(SSTdataParamsType) :: SSTdataParams(maxNumberSSTDatasets)

  ! CONVENTIONAL OBS ERRORS
  real(8) :: xstd_ua_ai_sw(20,11)
  real(8) :: xstd_sf(9,6)
  real(8) :: xstd_pr(2)
  real(8) :: xstd_sc(1)
  real(8) :: LVLVALUE(9), HGT_ERR(200,9)

  ! CONSTITUENT/CHEMISTRY OBS ERROR STD DEV
  ! Associated below to routines beginning with the prefix 'chm'
  type :: struct_chm_std
     !
     ! Structure containing information retrieved from auxiliary obs data file 
     ! holding observation std dev information for constituent obs
     !
     !  Variable               Description
     !  --------               -----------
     !  n_stnid                Number of sub-families (identified via STNIDs)
     !  stnids                 Sub-families (STNIDs; * are wild cards)
     !  element                BUFR element in data block 
     !  source                 0: Set entirely from the auxiliary file being read. No 
     !                            initial values read from observation files
     !                         1: Initial values in observation files 
     !                            (may be adjusted after input)
     !                         2: Initial values in observation files for variable number
     !                            of vertical levels (for error std deviations only)
     !  std_type               Index of setup approach (used in combination with source)
     !                         For source value 0 or 1, 
     !                         0: std1 or observation file values (sigma)
     !                         1: max(std3,std2*ZVAL)  if source=0
     !                            max(std3,std2*sigma) otherwise
     !                         2: sqrt(std3**2+(std2*ZVAL)**2))  if source=0
     !                            sqrt(std3**2+(std2*sigma)**2)) otherwise
     !                         3: min(std3,max(std2,std1_chm*ZVAL)) if source=0
     !                            min(std3,max(std2,std1_chm*sigma))  otherwise
     !                         4: sqrt(std2**2+(std1*ZVAL)**2))  if source=0 
     !                            sqrt(std2**2+(std1*sigma)**2)) otherwise
     !  ibegin                 Position index of start of data for given
     !                         sub-family in the arrays std1,levels,lat
     !  n_lvl                  Number of vertical levels (max number when source=2)
     !  levels                 Vertical levels (in coordinate of sub-family data)
     !  n_lat                  Number of latitudes
     !  lat                    Latitudes (degrees; ordered in increasing size)
     !  std1                   See std_type for usage (dependent on vertical level)
     !  std2                   See std_type for usage (independent of vertical level)
     !  std3                   See std_type for usage (independent of vertical level)

     integer ::  n_stnid
     character(len=12), allocatable :: stnids(:)
     integer, allocatable :: element(:),std_type(:),n_lat(:)
     integer, allocatable :: source(:),ibegin(:),n_lvl(:)
     real(8), allocatable :: std1(:),std2(:),std3(:)
     real(8), allocatable :: levels(:),lat(:)

     ! Array to hold std dev's read from auxiliary obs data/info file
     type(struct_oss_obsdata), allocatable :: obsStdDev(:)
     
  end type struct_chm_std

  type(struct_chm_std)  :: chm_std

  ! Sea Ice Concentration obs-error standard deviation
  real(8) :: xstd_sic(9)
  ! tiepoint standard deviation for ASCAT backscatter anisotropy
  integer, parameter :: ncells = 21
  real(8) :: ascatAnisSigmaOpenWater(ncells,12), ascatAnisSigmaIce(ncells,12)
  real    :: oer_ascatAnisOpenWater(ncells,12), oer_ascatAnisIce(ncells,12)

  ! Hydrology
  real(8) :: xstd_hydro(1)

  integer :: n_sat_type, n_categorie
  integer :: tbl_m(200), tbl_h(200), tbl_t(200), tbl_g(200)
  integer :: surfaceObsTypeNumber

  character(len=9) :: SAT_AMV(200,10), SAT_LIST(200), MET_LIST(200)
  character(len=9) :: HTM_LIST(200), TMG_LIST(200), NSW_LIST(200)

  logical :: new_oer_sw, obsfile_oer_sw, visAndGustAdded, useTovsUtil
  logical :: mwAllskyInflateByOmp, mwAllskyInflateByClwDiff
  real(8) :: amsuaClearClwThresholdSigmaObsInflation(5)
  real(8) :: amsuaStateDepSigmaObsInflationCoeff

  logical :: readOldSymmetricObsErrFile

  character(len=48) :: obserrorMode

contains

  !--------------------------------------------------------------------------
  ! oer_setInterchanCorr
  !--------------------------------------------------------------------------
  subroutine oer_setInterchanCorr()
    !
    !  :Purpose: Setup of interchannel observation errors correlations
    !    
    use rmatrix_mod
    IMPLICIT NONE

    INTEGER ::  ISENS

    if (tvs_nsensors == 0) then
      write(*,*) 'oer_setInterchanCorr: tvs_nsensors is zero, skipping setup'
      return
    end if

! Initialization of the correlation matrices
    call rmat_init(tvs_nsensors,tvs_nobtov)
    if (rmat_lnondiagr) then
      do isens = 1, tvs_nsensors
        if (tvs_isReallyPresent(isens)) call rmat_readCMatrix(tvs_listSensors(:,isens), isens, tvs_ichan(1:tvs_nchan(isens),isens))
      end do
    end if

  END SUBROUTINE oer_setInterchanCorr

  !--------------------------------------------------------------------------
  ! oer_setObsErrors
  !--------------------------------------------------------------------------
  subroutine oer_setObsErrors(obsSpaceData, obserrorMode_in, useTovsUtil_opt)
    !
    ! :Purpose: read and set observation errors (from former sucovo subroutine).
    !
    type(struct_obs)             :: obsSpaceData
    character(len=*), intent(in) :: obserrorMode_in
    logical, optional            :: useTovsUtil_opt

    namelist /namoer/ new_oer_sw, obsfile_oer_sw, visAndGustAdded
    namelist /namoer/ mwAllskyInflateByOmp, mwAllskyInflateByClwDiff
    namelist /namoer/ amsuaClearClwThresholdSigmaObsInflation
    namelist /namoer/ amsuaStateDepSigmaObsInflationCoeff
    namelist /namoer/ readOldSymmetricObsErrFile
    integer :: fnom, fclos, nulnam, ierr

    !
    !- 1.  Setup Mode
    !
    obserrorMode = obserrorMode_in

    ! Additional key to allow the use of 'util' column in stats_tovs file
    if (present(useTovsUtil_opt)) then
      useTovsUtil = useTovsUtil_opt
    else
      useTovsUtil = .false.
    end if

    ! read namelist namoer
    new_oer_sw = .false.
    obsfile_oer_sw  = .false.
    visAndGustAdded = .false.
    mwAllskyInflateByOmp = .false.
    mwAllskyInflateByClwDiff = .false.
    amsuaClearClwThresholdSigmaObsInflation(:) = 0.03D0
    amsuaClearClwThresholdSigmaObsInflation(1) = 0.05D0
    amsuaClearClwThresholdSigmaObsInflation(4) = 0.02D0
    amsuaStateDepSigmaObsInflationCoeff = 13.0D0
    readOldSymmetricObsErrFile = .true.

    if (utl_isNamelistPresent('namoer','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read (nulnam, nml = NAMOER, iostat = ierr)
      if (ierr /= 0) call utl_abort('oer_setObsErrors: Error reading namelist')
      if (mmpi_myid == 0) write(*,nml=namoer)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'oer_setObsErrors: namoer is missing in the namelist. The default value will be taken.'
    end if

    !
    !- 2.  Read in the observation std dev errors
    !

    !- 2.1 Radiance data
    if (obs_famexist(obsSpaceData,'TO')) then
      call oer_readObsErrorsTOVS
    else 
      write(*,*) "oer_setObsErrors: No brightness temperature observations found."
    end if
    
    !- 2.2 Conventional data
    if (obs_famExist(obsSpaceData, 'UA') .or. obs_famExist(obsSpaceData, 'AI') .or. obs_famExist(obsSpaceData, 'SW') .or. &
        obs_famExist(obsSpaceData, 'SF') .or. obs_famExist(obsSpaceData, 'GP') .or. obs_famExist(obsSpaceData, 'SC') .or. &
        obs_famExist(obsSpaceData, 'PR')) then

      call oer_readObsErrorsCONV()

    else

      write(*,*) "oer_setObsErrors: No conventional weather observations found."

    end if

    !- 2.3 Constituent data
    if (obs_famexist(obsSpaceData,'CH')) then
      call chm_read_obs_err_stddev
    else
      write(*,*) "oer_setObsErrors: No CH observations found."
    end if

    !- 2.4 Sea ice concentration
    if (obs_famexist(obsSpaceData,'GL')) then
      call oer_readObsErrorsIce
    else
      write(*,*) "oer_setObsErrors: No GL observations found."
    end if

    !- 2.5 SST
    if (obs_famexist(obsSpaceData,'TM')) then    
      call oer_readObsErrorsSST
    else
      write(*,*) "oer_setObsErrors: No TM observations found."
    end if

    !- 2.6 Hydrology
    if (obs_famexist(obsSpaceData,'HY')) then
      call oer_readObsErrorsHydro
    else
      write(*,*) "oer_setObsErrors: No HY observations found."
    end if

    !
    !- 3.  Set obs error information in obsSpaceData object
    !
    call oer_fillObsErrors(obsSpaceData)

    !
    !- 4.  Deallocate temporary storage
    !
    if (obs_famExist(obsSpaceData,'CH')) call chm_dealloc_obs_err_stddev

  end subroutine oer_setObsErrors

  !--------------------------------------------------------------------------
  ! oer_readObsErrorsTOVS
  !--------------------------------------------------------------------------
  subroutine oer_readObsErrorsTOVS
    !
    ! :Purpose: Read the observation error statistics and
    !           utilization flag for TOVS processing. This information
    !           resides on an ASCII file and is read using a free format.
    !
    implicit none

    integer, parameter :: bgckColumnIndex = 1
    integer, parameter :: analysisColumnIndex = 2

    integer,external  :: FNOM, FCLOS
    integer :: IER, ILUTOV, ILUTOV2, JI, obsErrorColumnIndex, JL, JM 
    integer :: INUMSAT, INUMSAT2, ISAT, IPLF
    integer :: IPLATFORM(tvs_maxNumberOfSensors), ISATID(tvs_maxNumberOfSensors)
    integer :: IINSTRUMENT(tvs_maxNumberOfSensors), NUMCHN(tvs_maxNumberOfSensors)
    integer :: NUMCHNIN(tvs_maxNumberOfSensors)
    integer :: IPLATFORM2, ISATID2, IINSTRUMENT2, NUMCHNIN2
    integer :: IUTILST(tvs_maxChannelNumber,tvs_maxNumberOfSensors)
    integer :: ICHN(tvs_maxChannelNumber,tvs_maxNumberOfSensors)
    integer :: ICHNIN(tvs_maxChannelNumber,tvs_maxNumberOfSensors)
    integer :: ICHNIN2(tvs_maxChannelNumber)

    real(8) :: TOVERRIN(tvs_maxChannelNumber,2,tvs_maxNumberOfSensors)
    real(8) :: clwThreshArrInput(tvs_maxChannelNumber,tvs_maxNumberOfSensors,2)
    real(8) :: sigmaObsErrInput(tvs_maxChannelNumber,tvs_maxNumberOfSensors,2)
    real(8) :: tovsObsInflation(tvs_maxChannelNumber,tvs_maxNumberOfSensors)
    real(8) :: clearClwThresholdSigmaObsInflationInput(tvs_maxChannelNumber,tvs_maxNumberOfSensors)
    real(8) :: stateDepSigmaObsInflationCoeffInput(tvs_maxNumberOfSensors)
    integer :: useStateDepSigmaObsInput(tvs_maxChannelNumber,tvs_maxNumberOfSensors)
    integer :: amsuaChannelOffset, amsuaChannelNum  

    character (len=132) :: CLDUM,CPLATF,CINSTR

    write(*,*) 'oer_readObsErrorsTOVS: reading observation error statistics required for TOVS processing'

    !
    !    1. Initialize
    !       ----------
    !
    TOVERRST(:,:) = 0.0D0
    TOVERRIN(:,:,:) = 0.0D0
    clwThreshArr(:,:,:) = 0.0d0
    clwThreshArrInput(:,:,:) = 0.0d0
    sigmaObsErr(:,:,:) = 0.0d0
    sigmaObsErrInput(:,:,:) = 0.0d0
    tovsObsInflation(:,:) = 0.0d0
    IUTILST(:,:) = 0
    useStateDepSigmaObs(:,:) = .false.
    useStateDepSigmaObsInput(:,:) = 0
    clearClwThresholdSigmaObsInflation(:,:) = 0.0d0
    clearClwThresholdSigmaObsInflationInput(:,:) = 0.0d0
    stateDepSigmaObsInflationCoeff(:) = 0.0d0
    stateDepSigmaObsInflationCoeffInput(:) = 0.0d0

    IPLATFORM(:) = 0
    NUMCHN(:) = 0
    NUMCHNIN(:) = 0
    ICHN(:,:) = 0
    ICHNIN(:,:) = 0
    ICHNIN2(:) = 0

    if (tvs_nobtov == 0) return

    !
    !     2. Open the file
    !        -------------
    !
    ilutov = 0
    IER =  FNOM(ILUTOV,'stats_tovs','SEQ+FMT',0)
    if (IER < 0) call utl_abort ('oer_readObsErrorsTOVS: Problem opening file stats_tovs')

    !
    !     3. Read number of satellites
    !        -------------------------
    !
    read(ILUTOV,*)
    read(ILUTOV,*) INUMSAT
    read(ILUTOV,*)

    if (inumsat > tvs_maxNumberOfSensors) then
      write(*,'(A)') ' Number of sensors in stats_tovs file is greater than maximum allowed (tvs_maxNumberOfSensors)'
      call utl_abort('oer_readObsErrorsTOVS')
    end if

    !
    !     4. Read the satellite identification, the number of channels,
    !        the observation errors and the utilization flags
    !        ----------------------------------------------------------
    !
    write(*,*) 'oer_readObsErrorsTOVS: Reading stats_tovs file.'

    DO JL = 1, INUMSAT

      read(ILUTOV,*)
      read(ILUTOV,'(A)') CLDUM
      write(*,'(A)') CLDUM
      CINSTR=CLDUM
      call split(CINSTR," ",CPLATF)
      write(*,*) "CINSTR: ",CINSTR
      write(*,*) "CPLATF: ",CPLATF
      read(ILUTOV,*)
      read(ILUTOV,*) ISATID(JL), NUMCHNIN(JL)

      do JI = 1, 3
        read(ILUTOV,*)
      end do

      if (CPLATF == 'FY-3C') THEN 
         CPLATF = 'FY3-3'
         CINSTR = 'MWHS2'
      end if

      IPLATFORM(JL) =  tvs_getPlatformId(CPLATF)

      if (IPLATFORM(JL) == -1) call utl_abort ('oer_readObsErrorsTOVS: Unknown platform!')

      IINSTRUMENT(JL) = tvs_getInstrumentId(CINSTR)
      
      if (IINSTRUMENT(JL) == -1) call utl_abort ('oer_readObsErrorsTOVS: Unknown instrument!')

      do JI = 1, NUMCHNIN(JL)
        read(ILUTOV,*) ICHNIN(JI,JL), TOVERRIN(ICHNIN(JI,JL),1,JL), TOVERRIN(ICHNIN(JI,JL),2,JL), IUTILST(ICHNIN(JI,JL),JL), tovsObsInflation(ICHNIN(JI,JL),JL)
      end do
      read(ILUTOV,*)

    end do

    ! read in the parameters to define the user-defined symmetric TOVS errors
    if (tvs_mwAllskyAssim) then
      ilutov2 = 10
      IER =  FNOM(ILUTOV2,'stats_tovs_symmetricObsErr','SEQ+FMT',0)
      if (IER < 0) call utl_abort ('oer_readObsErrorsTOVS: Problem opening symmetricObsErr file.')

      read(ILUTOV2,*)
      read(ILUTOV2,*) INUMSAT2
      read(ILUTOV2,*)

      write(*,*) 'oer_readObsErrorsTOVS: Reading symmetricObsErr file.'

      DO JL = 1, INUMSAT2

        read(ILUTOV2,*)
        read(ILUTOV2,'(A)') CLDUM
        write(*,'(A)') CLDUM
        CINSTR = CLDUM
        call split(CINSTR," ",CPLATF)
        write(*,*) "CINSTR: ",CINSTR
        write(*,*) "CPLATF: ",CPLATF
        read(ILUTOV2,*)

        ! If reading the old style stats_tovs_symmetricObsErr, the the all-sky parameters are available only for AMSUA.
        if (readOldSymmetricObsErrFile) then
          read(ILUTOV2,*) ISATID2, NUMCHNIN2
          if (CINSTR == "AMSUA") stateDepSigmaObsInflationCoeffInput(JL) = amsuaStateDepSigmaObsInflationCoeff
        else
          read(ILUTOV2,*) ISATID2, NUMCHNIN2, stateDepSigmaObsInflationCoeffInput(JL)
        end if

        if (ISATID2 /= ISATID(JL) .or. NUMCHNIN2 /= NUMCHNIN(JL)) &
          call utl_abort ('oer_readObsErrorsTOVS: problem with ISATID2, NUMCHNIN2 in symmetricObsErr')

        do JI = 1, 3
          read(ILUTOV2,*)
        end do

        IPLATFORM2 = tvs_getPlatformId(CPLATF)
        IINSTRUMENT2 = tvs_getInstrumentId(CINSTR)
        if (IPLATFORM2 /= IPLATFORM(JL) .or. IINSTRUMENT2 /= IINSTRUMENT(JL)) & 
          call utl_abort ('oer_readObsErrorsTOVS: problem with IPLATFORM2, IINSTRUMENT2 in symmetricObsErr')

        do JI = 1, NUMCHNIN2
          ! If reading the old style stats_tovs_symmetricObsErr, the the all-sky parameters are available only for AMSUA.
          if (readOldSymmetricObsErrFile) then
            read(ILUTOV2,*) ICHNIN2(JI), &
                  clwThreshArrInput(ICHNIN2(JI),JL,1), clwThreshArrInput(ICHNIN2(JI),JL,2), &
                  sigmaObsErrInput(ICHNIN2(JI),JL,1), sigmaObsErrInput(ICHNIN2(JI),JL,2), &
                  useStateDepSigmaObsInput(ICHNIN2(JI),JL)
            if (CINSTR == "AMSUA") then
              amsuaChannelOffset = 27
              amsuaChannelNum = ICHNIN2(JI) - amsuaChannelOffset
              if (amsuaChannelNum >= 1 .and.  amsuaChannelNum <= 5) then
                clearClwThresholdSigmaObsInflationInput(ICHNIN2(JI),JL) = &
                        amsuaClearClwThresholdSigmaObsInflation(amsuaChannelNum)
              end if
            end if
          else
            read(ILUTOV2,*) ICHNIN2(JI), &
                  clwThreshArrInput(ICHNIN2(JI),JL,1), clwThreshArrInput(ICHNIN2(JI),JL,2), &
                  sigmaObsErrInput(ICHNIN2(JI),JL,1), sigmaObsErrInput(ICHNIN2(JI),JL,2), &
                  clearClwThresholdSigmaObsInflationInput(ICHNIN2(JI),JL), &
                  useStateDepSigmaObsInput(ICHNIN2(JI),JL)
          end if

          if (ICHNIN2(JI) /= ICHNIN(JI,JL)) & 
            call utl_abort ('oer_readObsErrorsTOVS: problem with ICHNIN2 in symmetricObsErr')

        end do
        read(ILUTOV2,*)

      end do

      IER = FCLOS(ILUTOV2)
      if (IER /= 0) call utl_abort ('oer_readObsErrorsTOVS')

    end if

    !
    !   Select input error to use: if ANAL mode, use ERRANAL (obsErrorColumnIndex=2);
    !   otherwise use ERRBGCK (obsErrorColumnIndex=1)
    !
    if (trim(obserrorMode) == 'analysis' .or. trim(obserrorMode) == 'FSO') THEN
      obsErrorColumnIndex = analysisColumnIndex
    ELSE
      obsErrorColumnIndex = bgckColumnIndex
    end if

    !
    !   Fill the observation error array TOVERRST
    !
    write(*,*) 'oer_readObsErrorsTOVS: Fill error array TOVERRST.'
    do JM = 1, INUMSAT
      do JL = 1, tvs_nsensors
        if (tvs_platforms (JL) == IPLATFORM(JM) .AND. tvs_satellites(JL) == ISATID(JM)) THEN
          if (tvs_instruments (JL) == IINSTRUMENT(JM)) THEN
            NUMCHN(JL)=NUMCHNIN(JM)
            do JI = 1, tvs_maxChannelNumber
              TOVERRST(JI,JL) = TOVERRIN(JI,obsErrorColumnIndex,JM)
              ICHN(JI,JL) = ICHNIN(JI,JM)

              if (tvs_mwAllskyAssim) then
                clwThreshArr(JI,JL,:) = clwThreshArrInput(JI,JM,:)
                sigmaObsErr(JI,JL,:) = sigmaObsErrInput(JI,JM,:)
                clearClwThresholdSigmaObsInflation(JI,JL) = &
                        clearClwThresholdSigmaObsInflationInput(JI,JM)
                useStateDepSigmaObs(JI,JL) = &
                        (useStateDepSigmaObsInput(JI,JM) == 1)

                if (JI == 1) stateDepSigmaObsInflationCoeff(JL) = &
                                stateDepSigmaObsInflationCoeffInput(JM)

                ! inflate the sigmaObsErr in analysis mode
                if (obsErrorColumnIndex == analysisColumnIndex) then
                  sigmaObsErr(JI,JL,1) = sigmaObsErr(JI,JL,1) * tovsObsInflation(JI,JM)
                  sigmaObsErr(JI,JL,2) = sigmaObsErr(JI,JL,2) * tovsObsInflation(JI,JM)
                end if
              end if

            end do
            if (trim(obserrorMode) == 'bgck' .or. useTovsUtil) THEN
              do JI = 1, tvs_maxChannelNumber
                tovutil(JI,JL) =  IUTILST(JI,JM)
              end do
            end if
          end if
        end if
      end do
    end do

    !
    !  Check that oberservation error statistics have been defined for
    !  all the satellites specified in the namelist.
    !
    do JL = 1, tvs_nsensors
      IPLF = utl_findArrayIndex(IPLATFORM  , INUMSAT, tvs_platforms  (JL))
      ISAT =  utl_findArrayIndex(ISATID     , INUMSAT, tvs_satellites (JL))
      if (IPLF == 0 .OR. ISAT == 0) THEN
        write(*,'(A,I3)') 'oer_readObsErrorsTOVS: Observation errors not defined for sensor #', JL
        call utl_abort ('oer_readObsErrorsTOVS')
      end if
      if (NUMCHN(JL) == 0) THEN 
        write(*,'(A,I3)') 'oer_readObsErrorsTOVS: Problem setting errors for sensor #', JL
        call utl_abort ('oer_readObsErrorsTOVS')
      end if
    end do

    !
    !    5. Print out observation errors for each sensor
    !       --------------------------------------------
    !
    if (mmpi_myid == 0) THEN
      write(*,*) 'Radiance observation errors read from file'
      write(*,*) '------------------------------------------'
      do JL = 1, tvs_nsensors
        write(*,'(A,I2,4(A))') 'SENSOR #', JL, ', Platform: ', tvs_satelliteName(JL), &
                                ', Instrument: ',tvs_instrumentName(JL)
        if (tvs_mwAllskyAssim .and. any(useStateDepSigmaObs(ICHN(1:NUMCHN(JL),JL),JL))) then
          write(*,'(A,4(2X,A8),(1X,A9),(2X,A3))') 'Channel','clw1','clw2','sigmaO1','sigmaO2','anlErrInf','use'
          do JI = 1, NUMCHN(JL)
            write(*,'(I7,5(2X,F8.4),(2X,L3))') ICHN(JI,JL), &
              clwThreshArr(ICHN(JI,JL),JL,1), clwThreshArr(ICHN(JI,JL),JL,2), &
              sigmaObsErr(ICHN(JI,JL),JL,1), sigmaObsErr(ICHN(JI,JL),JL,2), &
              clearClwThresholdSigmaObsInflation(ICHN(JI,JL),JL), &
              useStateDepSigmaObs(ICHN(JI,JL),JL)
          end do
          write(*,'(A,(2X,F8.4))') 'stateDepSigmaObsInflationCoeff=', &
                stateDepSigmaObsInflationCoeff(JL) 
        else
          write(*,'(A,2X,A8)') 'Channel','sigmaO'
          do JI = 1, NUMCHN(JL)
            write(*,'(I7,1(2X,F8.4))') ICHN(JI,JL),TOVERRST(ICHN(JI,JL),JL)
          end do
        end if
      end do
    end if

    !
    !    6. Close the file
    !       --------------
    !
    IER = FCLOS(ILUTOV)
    if (IER /= 0) call utl_abort ('oer_readObsErrorsTOVS')

    !
    !    7. Temporary: Filll the publics variables for QC purpose
    !       --------------
    oer_toverrst(:,:) = toverrst(:,:)
    oer_tovutil (:,:) = tovutil(:,:)
    oer_clwThreshArr(:,:,:) = clwThreshArr(:,:,:)
    oer_sigmaObsErr(:,:,:) = sigmaObsErr(:,:,:)
    oer_useStateDepSigmaObs(:,:) = useStateDepSigmaObs(:,:)

  contains

    subroutine compact(str)
      ! Code from Benthien's module: http://www.gbenthien.net/strings/index.html
      ! Converts multiple spaces and tabs to single spaces; deletes control characters;
      ! removes initial spaces.

      character(len=*):: str
      character(len=1):: ch
      character(len=len_trim(str)):: outstr
      integer isp,k,lenstr,i,ich

      str=adjustl(str)
      lenstr=len_trim(str)
      outstr=' '
      isp=0
      k=0

      do i=1,lenstr
        ch=str(i:i)
        ich=iachar(ch)

        select case(ich)
        case(9,32)    ! space or tab character         
          if(isp==0) then
            k=k+1
            outstr(k:k)=' '
          end if
          isp=1
        case(33:)              ! not a space, quote, or control character
          k=k+1
          outstr(k:k)=ch
          isp=0
        end select

      end do

      str=adjustl(outstr)

    end subroutine compact

    subroutine split(str, delims, before)
      !
      ! :Comment:
      ! Code extracted from Benthien's module: http://www.gbenthien.net/strings/index.html
      ! Routine finds the first instance of a character from 'delims' in the
      ! the string 'str'. The characters before the found delimiter are
      ! output in 'before'. The characters after the found delimiter are
      ! output in 'str'. 
      !

      character(len=*) :: str
      character(len=*) :: delims
      character(len=*) :: before

      character :: ch,cha
      integer lenstr,i,k,ipos,iposa
      str=adjustl(str)
      call compact(str)
      lenstr=len_trim(str)

      if(lenstr == 0) return ! string str is empty
      k=0
      before=' '
      do i=1,lenstr
        ch=str(i:i)
        
        ipos=index(delims,ch)

        if(ipos == 0) then ! character is not a delimiter
          k=k+1
          before(k:k)=ch
          cycle
        end if
        if(ch /= ' ') then ! character is a delimiter that is not a space
          str=str(i+1:)
          exit
        end if

        cha=str(i+1:i+1)  ! character is a space delimiter
        iposa=index(delims,cha)
        if(iposa > 0) then   ! next character is a delimiter 
          str=str(i+2:)
          exit
        else
          str=str(i+1:)
          exit
        end if
      end do
      if(i >= lenstr) str=''

      str=adjustl(str) ! remove initial spaces

    end subroutine split

  end subroutine oer_readObsErrorsTOVS

  !--------------------------------------------------------------------------
  ! oer_readObsErrorsCONV
  !--------------------------------------------------------------------------
  subroutine oer_readObsErrorsCONV()
    ! 
    ! :Purpose: read observation errors (modification of former readcovo subroutine) of conventional data 
    !
    implicit none

    integer :: fnom, fclos, ierr, jlev, jelm, jcat, icodtyp, nulstat
    logical             :: LnewExists
    character (len=128) :: ligne

    if (visAndGustAdded) then
      surfaceObsTypeNumber = 6
    else
      surfaceObsTypeNumber = 4
    end if

    ! CHECK THE EXISTENCE OF THE NEW FILE WITH STATISTICS
    inquire(file = 'obserr', exist = LnewExists)
    if (LnewExists) then
      write(*,*) '--------------------------------------------------------'
      write(*,*) 'oer_readObsErrorsCONV: reads observation errors in obserr'
      write(*,*) '--------------------------------------------------------'
    else
      call utl_abort('oer_readObsErrorsCONV: NO OBSERVATION STAT FILE FOUND!!')     
    end if

    ! Read observation errors from file obserr for conventional data
    nulstat=0
    ierr=fnom(nulstat, 'obserr', 'SEQ', 0)
    if (ierr == 0) then
      write(*,*) 'oer_readObsErrorsCONV: File =  ./obserr'
      write(*,*) ' opened as unit file ',nulstat
      open(unit=nulstat, file='obserr', status='OLD')
    else
      call utl_abort('oer_readObsErrorsCONV: COULD NOT OPEN FILE obserr!!!')
    end if

    write(*, '(A)') ' '

    do jlev = 1,3
      read(nulstat, '(A)') ligne
      write(*, '(A)') ligne
    end do

    do jlev = 1, 19
      read(nulstat, *) (xstd_ua_ai_sw(jlev,jelm), jelm=1,11)
      write(*, '(f6.0,10f6.1)')  (xstd_ua_ai_sw(jlev,jelm), jelm=1,11)
    end do

    do jlev = 1,5
      read(nulstat, '(A)') ligne
      write(*, '(A)') ligne
    end do

    read(nulstat, *) xstd_pr(1),xstd_pr(2)
    write(*, '(2f6.1)')  xstd_pr(1),xstd_pr(2)

    do jlev = 1,5
      read(nulstat, '(A)') ligne
      write(*, '(A)') ligne
    end do

    read(nulstat, *) xstd_sc(1)
    write(*, '(f8.3)')  xstd_sc(1)

    read(nulstat, '(A)') ligne
    write(*, '(A)') ligne

    do icodtyp = 1,9
      do jlev = 1,4
        read(nulstat, '(A)') ligne
        write(*, '(A)') ligne
      end do
      read(nulstat, *) (xstd_sf(icodtyp,jelm), jelm=1,surfaceObsTypeNumber)
      write(*, '(f6.2,2f6.1,f8.3)')  (xstd_sf(icodtyp,jelm), jelm=1,surfaceObsTypeNumber)
    end do

    !
    !     Read height assignment errors
    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !

    if(new_oer_sw) then

      do jlev = 1,5
        read(nulstat, '(A)') ligne
        write(*, '(A)') ligne
      end do

      read(nulstat, *) (LVLVALUE(jelm), jelm=1,9)
      write(*, '(9f6.1)')  (LVLVALUE(jelm), jelm=1,9)

      do jlev = 1,2
        read(nulstat, '(A)') ligne
        write(*, '(A)') ligne
      end do

      read(nulstat, '(i10)') n_sat_type
      write(*,'(i10)') n_sat_type
      do jlev = 1,n_sat_type
        read(nulstat, '(10A10)') (SAT_AMV(jlev,jelm), jelm=1,10)
        write(*, '(10A10)') (SAT_AMV(jlev,jelm), jelm=1,10)
      end do

      read(nulstat, '(A)') ligne
      write(*, '(A)') ligne

      tbl_m(:) = 0
      tbl_h(:) = 0
      tbl_t(:) = 0
      tbl_g(:) = 0
      read(nulstat, '(i10)') n_categorie
      write(*,'(i10)') n_categorie
      do jcat = 1,n_categorie
        read(nulstat, *) SAT_LIST(jcat),MET_LIST(jcat),HTM_LIST(jcat),TMG_LIST(jcat),NSW_LIST(jcat),(HGT_ERR(jcat,jelm), jelm=1,9)
        write(*, '(5A10,9f6.1)') SAT_LIST(jcat),MET_LIST(jcat),HTM_LIST(jcat),TMG_LIST(jcat),NSW_LIST(jcat),(HGT_ERR(jcat,jelm), jelm=1,9)
        if(trim(MET_LIST(jcat)) == 'ir')   tbl_m(jcat) = 1
        if(trim(MET_LIST(jcat)) == 'vis')  tbl_m(jcat) = 2
        if(trim(MET_LIST(jcat)) == 'wv')   tbl_m(jcat) = 3
        if(trim(HTM_LIST(jcat)) == 'ebbt') tbl_h(jcat) = 1
        if(trim(HTM_LIST(jcat)) == 'co2')  tbl_h(jcat) = 4
        if(trim(TMG_LIST(jcat)) == 'sea')  tbl_t(jcat) = 1
        if(trim(TMG_LIST(jcat)) == 'land') tbl_t(jcat) = 2 
        if(trim(TMG_LIST(jcat)) == 'ice')  tbl_t(jcat) = 3
        if(trim(NSW_LIST(jcat)) == 'NH')   tbl_g(jcat) = 1
        if(trim(NSW_LIST(jcat)) == 'SH')   tbl_g(jcat) = 2
      end do

    end if

    write(*, '(A)') ' '

    close(unit = nulstat)
    ierr = fclos(nulstat)

  end subroutine oer_readObsErrorsCONV

  !--------------------------------------------------------------------------
  ! oer_readObsErrorsIce
  !--------------------------------------------------------------------------
  subroutine oer_readObsErrorsIce
    !
    ! :Purpose: read observation errors for sea ice concentration analysis
    !
    implicit none

    external fnom, fclos
    integer :: fnom, fclos, ierr, jlev, jelm, nulstat
    logical            :: fileExists
    character(len=128) :: ligne
    character(len=15), parameter :: fileName = 'sea_ice_obs-err'

    ! Variables for the ASCAT backscatter anisotropy
    integer :: jcell_no, icell_no, imonth
    real :: tiePoint12, tiePoint13, tiePoint23

    ! CHECK THE EXISTENCE OF THE FILE WITH STATISTICS
    inquire(file = fileName, exist = fileExists)
    if (fileExists) then
      write(*,*) '--------------------------------------------------------'
      write(*,*) 'oer_readObsErrorsICE: reads observation errors in ',fileName
      write(*,*) '--------------------------------------------------------'
    else
      call utl_abort('oer_readObsErrorsICE: NO OBSERVATION STAT FILE FOUND!!')     
    end if

    ! Read observation errors from file
    nulstat=0
    ierr=fnom(nulstat, fileName, 'SEQ', 0)
    if (ierr == 0) then
      write(*,*) 'oer_readObsErrorsICE: File = ./',fileName
      write(*,*) ' opened as unit file ',nulstat
      open(unit=nulstat, file=fileName, status='OLD')
    else
      call utl_abort('oer_readObsErrorsICE: COULD NOT OPEN FILE '//fileName//'!!!')
    end if

    if (mmpi_myid == 0) write(*, '(A)') ' '

    do jelm = 1, 9

      do jlev = 1, 3
        read(nulstat, '(A)') ligne
        if (mmpi_myid == 0) write(*, '(A)') ligne
      end do

      read(nulstat, *) xstd_sic(jelm)
      if (mmpi_myid == 0) write(*,*) xstd_sic(jelm)

      do jlev = 1, 2
        read(nulstat, '(A)') ligne
        if (mmpi_myid == 0) write(*, '(A)') ligne
      end do

    end do

    ! Read coefficients for ASCAT backscatter anisotropy
    do jlev = 1, 5
      read(nulstat, '(A)') ligne
      if (mmpi_myid == 0) write(*, '(A)') ligne
    end do

    do imonth = 1, 12
      do jlev = 1, 3
        read(nulstat, '(A)') ligne
        if (mmpi_myid == 0) write(*, '(A)') ligne
      end do
      do jcell_no = 1, ncells
         read(nulstat, *) icell_no, tiePoint12, tiePoint13, tiePoint23, &
                           oer_ascatAnisIce(jcell_no,imonth), oer_ascatAnisOpenWater(jcell_no,imonth), &
                           ascatAnisSigmaIce(jcell_no,imonth), ascatAnisSigmaOpenWater(jcell_no,imonth)
         if (mmpi_myid == 0) write(*,*) icell_no, tiePoint12, tiePoint13, tiePoint23, &
                           oer_ascatAnisIce(jcell_no,imonth), oer_ascatAnisOpenWater(jcell_no,imonth), &
                           ascatAnisSigmaIce(jcell_no,imonth), ascatAnisSigmaOpenWater(jcell_no,imonth)
      end do
    end do

    if (mmpi_myid == 0) write(*, '(A)') ' '

    close(unit = nulstat)
    ierr = fclos(nulstat)

  end subroutine oer_readObsErrorsIce

  !--------------------------------------------------------------------------
  ! oer_readObsErrorsSST
  !--------------------------------------------------------------------------
  subroutine oer_readObsErrorsSST
    !
    ! :Purpose: read observation errors for SST data
    !
    implicit none

    integer :: fnom, fclos, nulnam, ierr, indexData
    namelist /namSSTObsErrors/ numberSSTDatasets, SSTdataParams
    
    if (utl_isNamelistPresent('namSSTObsErrors','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
      read (nulnam, nml = namSSTObsErrors, iostat = ierr)
      if (ierr /= 0) call utl_abort('oer_readObsErrorsSST: Error reading namelist')
      if (mmpi_myid == 0) write(*,nml=namSSTObsErrors)
      ierr = fclos(nulnam)
    else
       call utl_abort('oer_readObsErrorsSST: namSSTObsErrors is missing in the namelist.')
    end if
    
  end subroutine oer_readObsErrorsSST

  !--------------------------------------------------------------------------
  ! oer_readObsErrorsHydro
  !--------------------------------------------------------------------------
  subroutine oer_readObsErrorsHydro
    implicit none

    external fnom, fclos

    integer                      :: fnom, fclos, ierr, nulstat
    logical                      :: fileExists
    character(len=15), parameter :: fileName = 'obserr_hydro'    
    character(len=*) , parameter :: myName   = 'oer_readObsErrorsHydro'

    inquire(file = fileName, exist = fileExists)
    if (fileExists) then
      write(*,*) '--------------------------------------------------------'
      write(*,*) myName//': reads observation errors in ', fileName
      write(*,*) '--------------------------------------------------------'
    else
      call utl_abort(myName//': NO OBSERVATION STAT FILE FOUND!!')     
    end if

    nulstat=0
    ierr=fnom(nulstat, fileName, 'SEQ', 0)
    if (ierr == 0) then
      write(*,*) myName//': File = ./',fileName
      write(*,*) myName//' opened as unit file ',nulstat
      open(unit=nulstat, file=fileName, status='OLD')
    else
      call utl_abort(myName//': COULD NOT OPEN FILE '//fileName//'!!!')
    end if

    read(nulstat, *)    xstd_hydro(1)
    write(*, '(f8.3)')  xstd_hydro(1)

    close(unit = nulstat)
    ierr = fclos(nulstat)

  end subroutine oer_readObsErrorsHydro

  !--------------------------------------------------------------------------
  ! oer_fillObsErrors
  !--------------------------------------------------------------------------
  subroutine oer_fillObsErrors(obsSpaceData)
    !
    ! :Purpose: fill observation errors in obsSpaceData
    !
    implicit none

    ! arguments
    type(struct_obs) :: obsSpaceData

    !  locals
    integer :: jn, JI, bodyIndex, headerIndex, ityp, iass, idata, idatend, codeType
    integer :: sensorIndex 
    integer :: isat, channelNumber, iplatf, instr, iplatform, instrum
    integer :: ilev, nlev, idate, itime
    integer :: ielem, icodtyp, header_prev, indexDataset, indexSensor

    real(8) :: zlat, zlon, zlev, zval, zwb, zwt, obs_err_stddev, solarZenith
    real(8) :: obsValue, obsStdDevError, obsPPP, obsOER
    real(8) :: clwThresh1, clwThresh2, clw_avg, clwObs, clwFG
    real(8) :: sigmaThresh1, sigmaThresh2, sigmaObsErrUsed
    real(8), parameter :: minRetrievableClwValue = 0.0D0
    real(8), parameter :: maxRetrievableClwValue = 3.0D0

    logical :: ifirst, surfTypeIsWater, unsupportedCodeType, unsupportedSensor

    character(len=2)  :: cfam
    character(len=12) :: cstnid

    write(*,'(10X, "***********************************")')
    write(*,'(10X, "oer_fillObsErrors: starting", /)')
    write(*,'(10X, "***********************************")')

    header_prev=-1

    ! SET STANDARD DEVIATION ERRORS FOR EACH DATA FAMILY
    HEADER: do headerIndex = 1, obs_numheader(obsSpaceData)

      idata    = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
      idatend  = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) + idata - 1
      cfam     = obs_getFamily (obsSpaceData, headerIndex)
      zlat     = obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex)
      zlon     = obs_headElem_r(obsSpaceData, OBS_LON, headerIndex)
      codeType = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      iplatf   = obs_headElem_i(obsSpaceData, OBS_SAT, headerIndex)
      instr    = obs_headElem_i(obsSpaceData, OBS_INS, headerIndex)
      cstnid   = obs_elem_c    (obsSpaceData, 'STID' , headerIndex)
      idate    = obs_headElem_i(obsSpaceData, OBS_DAT, headerIndex) 
      itime    = obs_headElem_i(obsSpaceData, OBS_ETM, headerIndex)

      surfTypeIsWater = (tvs_ChangedStypValue(obsSpaceData,headerIndex) == surftype_sea)

      nlev = idatend - idata + 1

      BODY: do bodyIndex  = idata, idatend

        ityp  = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex)
        iass  = obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex)
        zval  = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex)

        if (iass == obs_assimilated) then

             !***********************************************************************
             !                           TOVS DATA
             !***********************************************************************

          if (cfam == 'TO') then

            if (ityp == BUFR_NBT1 .or. &
                ityp == BUFR_NBT2 .or. &
                ityp == BUFR_NBT3) then

              channelNumber = NINT(obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex))
              call tvs_mapSat(iplatf, iplatform, isat)
              call tvs_mapInstrum(instr, instrum)

              do sensorIndex = 1, tvs_nsensors
                if (iplatform ==  tvs_platforms(sensorIndex)  .and. &
                    isat      ==  tvs_satellites(sensorIndex) .and. &
                    instrum   == tvs_instruments(sensorIndex)) then

                  ! decide whether or not use the state dependent sigmaObsErrUsed for OBS_OER
                  if (tvs_mwAllskyAssim .and. &
                      useStateDepSigmaObs(channelNumber,sensorIndex) .and. &
                      surfTypeIsWater) then

                    ! set dummy value for OBS_OER in bgck mode
                    if (trim(obserrorMode) == 'bgck') then
                      sigmaObsErrUsed = 1.0d0
                    else
                      clwThresh1 = clwThreshArr(channelNumber,sensorIndex,1)
                      clwThresh2 = clwThreshArr(channelNumber,sensorIndex,2)
                      sigmaThresh1 = sigmaObsErr(channelNumber,sensorIndex,1)
                      sigmaThresh2 = sigmaObsErr(channelNumber,sensorIndex,2)
                      clwObs = obs_headElem_r(obsSpaceData, OBS_CLWO, headerIndex)
                      clwFG  = obs_headElem_r(obsSpaceData, OBS_CLWB, headerIndex)
                      clw_avg = 0.5D0 * (clwObs + clwFG)

                      ! check to ensure CLW is retrieved and properly set
                      if (clw_avg < minRetrievableClwValue .or. &
                          clw_avg > maxRetrievableClwValue) then
                        write(*,*) 'This observation should have been rejected ', &
                                  'in all-sky mode at background check!' 
                        write(*,*) 'oer_fillObsErrors: clwObs=', clwObs, &
                                  ', clwFG=', clwFG
                        call utl_abort('oer_fillObsErrors: CLW is not usable to define obs error')
                      end if

                      sigmaObsErrUsed = calcStateDepObsErr(clwThresh1, clwThresh2, &
                                sigmaThresh1,sigmaThresh2,clw_avg)
                    end if
                  else
                    sigmaObsErrUsed = TOVERRST(channelNumber, sensorIndex)
                  end if
                  call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, sigmaObsErrUsed)

                  !   Utilization flag for AIRS,IASI and CrIS channels (bgck mode only)
                  if (trim(obserrorMode) == 'bgck' .or. useTovsUtil) then
                    if  (tovutil(channelNumber, sensorIndex) == 0) &
                      call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex, ibset(obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex), 8))
                  end if
                end if
              end do

            end if

                !***********************************************************************
                !                      RADIOSONDE DATA
                !***********************************************************************

          else if (cfam == 'UA') then

            zlev = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex)

            if ((ityp == BUFR_NEUS) .or. (ityp == BUFR_NEVS)) then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sf(1,4))
            else if (ityp == BUFR_NETS) then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sf(1,2))
            else if (ityp == BUFR_NESS) then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sf(1,3))
            else if (ityp == BUFR_NEPS) then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sf(1,1))
            else if (ityp == BUFR_NEPN) then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sf(1,1))
            else

              if ((ityp == BUFR_NEUU) .or. (ityp == BUFR_NEVV)) then
                ielem = 4
              else if (ityp == BUFR_NETT) then
                ielem = 2
              else if (ityp == BUFR_NEES) then
                ielem = 3
              else if (ityp == BUFR_NEGZ) then
                ielem = 5
              end if

              if ((zlev * MPC_MBAR_PER_PA_R8) >= xstd_ua_ai_sw(1,1)) then

                call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_ua_ai_sw(1, ielem))

              else if ((zlev * MPC_MBAR_PER_PA_R8) <= xstd_ua_ai_sw(19, 1)) then

                call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_ua_ai_sw(19, ielem))

              else

                do jn = 1, 18
                  if ((zlev * MPC_MBAR_PER_PA_R8) >= xstd_ua_ai_sw(jn + 1, 1)) exit
                end do

                zwb = log((zlev * MPC_MBAR_PER_PA_R8) / xstd_ua_ai_sw(jn, 1)) / log(xstd_ua_ai_sw(jn + 1, 1) / xstd_ua_ai_sw(jn, 1))
                zwt = 1.0D0 - zwb
                
                call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex,  &
                                    zwt * xstd_ua_ai_sw(jn, ielem) +  &
                                    zwb * xstd_ua_ai_sw(jn + 1, ielem))

              end if
                   
            end if

                !***********************************************************************
                !                          AMV, AIREP, AMDAR DATA
                !***********************************************************************

          else if (cfam == 'AI' .or. cfam == 'SW') then

            zlev=obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex)

            if (codeType == codtyp_get_codtyp('windsbufr')) then ! AMV
          
              if ((ityp == BUFR_NEUU) .or. (ityp == BUFR_NEVV)) ielem = 11
       
            else if (codeType == codtyp_get_codtyp('airep')) then ! AIREP
              if ((ityp == BUFR_NEUU) .or. (ityp == BUFR_NEVV)) then
                ielem = 7
              else if (ityp == BUFR_NETT) then
                ielem = 6
              end if
              
            else if (codeType == codtyp_get_codtyp('amdar') .or. &
                     codeType == codtyp_get_codtyp('acars') .or. &
                     codeType == codtyp_get_codtyp('ads')) then
              if (ityp == BUFR_NEUU .or. ityp == BUFR_NEVV) then
                ielem = 10
              else if (ityp == BUFR_NETT) then
                ielem = 8
              else if (ityp == BUFR_NEES) then
                ielem = 9
              end if
            end if

            if ((zlev * MPC_MBAR_PER_PA_R8) >= xstd_ua_ai_sw(1, 1)) then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_ua_ai_sw(1, ielem))
            else if ((zlev * MPC_MBAR_PER_PA_R8) <= xstd_ua_ai_sw(19, 1)) then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_ua_ai_sw(19, ielem))
            else
              do jn = 1, 18
                if ((zlev * MPC_MBAR_PER_PA_R8) >= xstd_ua_ai_sw(jn + 1, 1)) exit
              end do

              zwb = log((zlev * MPC_MBAR_PER_PA_R8) / xstd_ua_ai_sw(jn, 1)) / log(xstd_ua_ai_sw(jn + 1, 1) / xstd_ua_ai_sw(jn, 1))
              zwt = 1.0D0 - zwb 

              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex,  &
                                  zwt * xstd_ua_ai_sw(jn, ielem) +  &
                                  zwb * xstd_ua_ai_sw(jn + 1, ielem))

            end if

                !***********************************************************************
                !                         SURFACE DATA
                !***********************************************************************

          else if (cfam == 'SF') then
            icodtyp = 1   ! Default values
            if (codeType == codtyp_get_codtyp('synopnonauto')) icodtyp = 2 ! SYNOP
            if (codeType == codtyp_get_codtyp('shipnonauto'))  icodtyp = 3 ! SHIP NON-AUTOMATIQUE
            if (codeType == codtyp_get_codtyp('synopmobil'))   icodtyp = 4 ! DRIBU
            if (codeType == codtyp_get_codtyp('drifter'))      icodtyp = 5 ! DRifTER
            if (codeType == codtyp_get_codtyp('synoppatrol'))  icodtyp = 6 ! STATION AUTOMATIQUE
            if (codeType == codtyp_get_codtyp('asynopauto'))   icodtyp = 7 ! ASYNOP
            if (codeType == codtyp_get_codtyp('ashipauto'))    icodtyp = 8 ! ASHIP
            if (ityp == BUFR_NEUU .or. ityp == BUFR_NEVV .or. &
                ityp == BUFR_NEGZ .or. ityp == BUFR_NETT .or. ityp == BUFR_NEES) icodtyp = 9  ! Others
            if (ityp == BUFR_NEUS .or. ityp == BUFR_NEVS) then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sf(icodtyp, 4))
            else if (ityp == BUFR_NETS) then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sf(icodtyp, 2))
            else if (ityp == BUFR_NESS) then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sf(icodtyp, 3))
            else if (ityp == BUFR_NEPS) then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sf(icodtyp, 1))
            else if (ityp == BUFR_NEPN) then
              if (icodtyp == 2  .or. icodtyp == 7) then
                call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sf(1, 1))
              else
                call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sf(icodtyp, 1))
              end if
            else if ((ityp == BUFR_NEUU) .or. (ityp == BUFR_NEVV))then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sf(icodtyp, 4))
            else if (ityp == BUFR_NEGZ) then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sf(icodtyp, 1))
            else if (ityp == BUFR_NETT) then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sf(icodtyp, 2))
            else if (ityp == BUFR_NEES) then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sf(icodtyp, 3))
            else if (ityp == bufr_vis .or. ityp == bufr_logVis) then
              if (surfaceObsTypeNumber >= 5) then
                call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sf(icodtyp, 5))
              else
                call utl_abort('oer_fillObsErrors: observation error missing for visibility')
              end if
            else if (ityp == bufr_gust) then
              if (surfaceObsTypeNumber >= 6) then
                call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sf(icodtyp, 6))
              else
                call utl_abort('oer_fillObsErrors: observation error missing for wind gust')
              end if
            else if (ityp == BUFR_NEFS) then ! SAR wind speed
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, 2.0d0)
            end if

                !***********************************************************************
                !                             GPS RO DATA
                !***********************************************************************

          else if (cfam == 'RO') then
            ! Process only refractivity data (codtyp 169)
            if (codeType == codtyp_get_codtyp('ro')) then
              if (ityp == BUFR_NEPS) then
                call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, 50.D0)
              end if
              if (ityp == BUFR_NETT) then
                call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, 10.D0)
              end if
              if (ityp == BUFR_NERF) then
                call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, 500.D0)
              end if
              if (ityp == BUFR_NEBD) then
                call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, 0.08D0)
              end if
            end if

          !***********************************************************************
          !                   GB-GPS SFC MET DATA
          !***********************************************************************

          !              ERRORS ARE SET TO SYNO SFC OBS ERRORS FROM S/R SUCOVO
          !              AND WEIGHTED BY FACTOR YSFERRWGT FOR 3D-VAR FGAT OR 4D-VAR ASSIM.
          !              OF TIME-SERIES (YSFERRWGT = 1.0 FOR 3D THINNING) 
          !
          else if (cfam == 'GP') then

            if (ityp == BUFR_NEPS) then ! Psfc Error (Pa)
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sf(2, 1))
            end if
            if (ityp == BUFR_NETS) then ! Tsfc Error (K)
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sf(2, 2))
            end if
            if (ityp == BUFR_NESS) then ! T-Td Error (K)
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sf(2, 3))
            end if
                ! ZTD Error (m) (value is formal error, real error set later in s/r seterrgpsgb)
                ! If error is missing, set to dummy value (1 m).
            if (ityp == BUFR_NEZD) then
              if (obs_bodyElem_r(obsSpaceData, OBS_OER, bodyIndex) <=  0.0D0) then
                call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, 1.0D0)
              end if
            end if

                !***********************************************************************
                !        SCATTEROMETER, WIND PROFILER DATA
                !***********************************************************************

          else if (cfam == 'SC') then

            call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sc(1))

          else if (cfam == 'PR') then

            zlev= obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex) - obs_headElem_r(obsSpaceData, OBS_ALT, headerIndex)

            if (zlev  >=  6000.) then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_pr(2))
            else
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_pr(1))
            end if

                !***********************************************************************
                !        ALADIN HORIZONTAL LINE-OF-SIGHT WIND DATA
                !***********************************************************************

          else if (cfam == 'AL') then

            ! TEMPORARILY, hard-code the observation error of AL winds to 2 m/s
            call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, 2.0d0)
              
                !***********************************************************************
                !               CONSTITUENT DATA (OZONE AND OTHER CHEMICALS)
                !***********************************************************************

          else if (cfam == 'CH') then

                !        Process only retrieved constituent data
                !        Also, exclude BUFR_SCALE_EXPONENT element as a data value!

            if (codeType == codtyp_get_codtyp('CHEMREMOTE') .or. codeType == codtyp_get_codtyp('CHEMINSITU')) then

              ifirst = headerIndex /= header_prev
              if (ifirst) then

                header_prev = headerIndex

                ! Subtract the number of occurences of code BUFR_SCALE_EXPONENT from number of levels
                do ji = 0, nlev - 1
                  if (obs_bodyElem_i(obsSpaceData, OBS_VNM, idata + ji) == BUFR_SCALE_EXPONENT) nlev = nlev-1
                end do
              end if

              zlev   = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)

              ilev = 0
              do ji = idata, bodyIndex
                if (obs_bodyElem_i(obsSpaceData, OBS_VNM, ji) /= BUFR_SCALE_EXPONENT) ilev = ilev+1
              end do

              obs_err_stddev = chm_get_obs_err_stddev(cstnid, nlev, ityp, zlat, zlon, idate, itime, zval, zlev, ilev, ifirst)
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, obs_err_stddev)

            end if     

                !***********************************************************************
                !               Sea Surface Temperature
                !***********************************************************************

          else if (cfam == 'TM') then

            if (obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex) /= bufr_sst) cycle BODY

            if (codeType == codtyp_get_codtyp('drifter')   .or. codeType == codtyp_get_codtyp('shipnonauto') .or. &
                codeType == codtyp_get_codtyp('ashipauto') .or. codeType == codtyp_get_codtyp('pseudosfc')        ) then

              unsupportedCodeType = .true.
              dataset_loop: do indexDataset = 1, numberSSTDatasets
                if(codeType == SSTdataParams(indexDataset)%codeType) then
                  unsupportedCodeType = .false.
                  call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, SSTdataParams(indexDataset)%dayError)
                  exit dataset_loop
                end if
              end do dataset_loop

              if(unsupportedCodeType) then
                write(*,'(a,i5,2a)') 'oer_fillObsErrors: unsupported SST data, codtype ', codeType,', cstnid ', cstnid,' found in dataset_loop!' 
                call utl_abort('oer_fillObsErrors: unsupported codeType!')
              end if

            else if (codeType == codtyp_get_codtyp('satob')) then

              solarZenith = obs_headElem_r(obsSpaceData, obs_sun, headerIndex)
              if (solarZenith == MPC_missingValue_R8) then
                write(*,*) 'oer_fillObsErrors: Solar zenith value is missing for satellite SST data ', &
                           cstnid, ', headerIndex: ', headerIndex,', bodyIndex: ', bodyIndex
                call utl_abort('oer_fillObsErrors: Solar zenith value is missing!')
              end if
              
              unsupportedSensor = .true.
              sensor_loop: do indexSensor = 1, numberSSTDatasets
                if (cstnid == trim(SSTdataParams(indexSensor)%sensor)) then
                  unsupportedSensor = .false.
                  if (solarZenith <= 90.) then ! day
                    call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, SSTdataParams(indexSensor)%dayError)
                  else ! night
                    call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, SSTdataParams(indexSensor)%nightError)
                  end if
                  exit sensor_loop
                end if
              end do sensor_loop

              if(unsupportedSensor) then
                write(*,'(3a)') 'oer_fillObsErrors: unsupported satellite SST data sensor ', cstnid,' found in sensor_loop!' 
                call utl_abort('oer_fillObsErrors: unsupported satellite SST data sensor!')
              end if

            else

              write(*,'(a,i5,2a)') 'oer_fillObsErrors: unsupported SST data, codtype ', codeType,', cstnid ', cstnid,' found!' 
              call utl_abort('oer_fillObsErrors: unsupported codeType!')

            end if

                !***********************************************************************
                !               Sea Ice Concentration
                !***********************************************************************

          else if (cfam == 'GL') then

            if (index(cstnid,'DMSP') == 1) then
              select case(cstnid)
              case('DMSP15')
                call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sic(1))
              case('DMSP16','DMSP17','DMSP18')
                call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sic(2))
              case DEFAULT
                call utl_abort('oer_fillObsErrors: UNKNOWN station id: '//cstnid)
              end select
            else if (cstnid == 'GCOM-W1') then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sic(3))
            else if (cstnid(1:6) == 'METOP-') then
              if (ityp == BUFR_ICEP) then
                 call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sic(4))
              else if (ityp == BUFR_ICES) then
                 ! This is backscatter anisotropy, obs-error will be set in oer_setErrBackScatAnisIce
                 ! because the need of background ice concentration.
                 ! For now just put a large positive value
                 call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, 1000.0)
              else
                 call utl_abort('oer_fillObsErrors: UNKNOWN varno: '// trim(utl_str(ityp)) // 'station id: '//cstnid)
              end if
            else if (cstnid == 'noaa-19') then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sic(5))
            else if (cstnid == 'CIS_DAILY') then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sic(6))
            else if (cstnid == 'RS1_IMG') then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sic(7))
            else if (codeType ==  178) then ! lake ice
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sic(8))
            else if (cstnid == 'CIS_REGIONAL') then
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, xstd_sic(9))
            else
              call utl_abort('oer_fillObsErrors: UNKNOWN station id: '//cstnid)
            end if
            
            !***********************************************************************
            !               Hydrology
            !***********************************************************************
          else if (cfam == 'HY') then

            if (obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex) /= bufr_riverFlow) cycle BODY

            if (codeType == 12) then ! pseudo-SYNOP
              obsValue = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex)
              obsStdDevError = obsValue * xstd_hydro(1)
              write(*,*) 'Hydro observation std dev error: ', bodyIndex, obsValue, xstd_hydro(1), obsStdDevError
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, obsStdDevError)
            else
              write(*,*) 'oer_fillObsErrors: unsupported codeType for hydro data found in the observations: ', codeType 
              call utl_abort('oer_fillObsErrors: unsupported codeType')
            end if  

          else if (cfam == 'RA') then

            if (ityp == BUFR_radarPrecip .or. ityp == BUFR_logRadarPrecip) then
              ! Temporary hardcoded value for log-transformed radar precipitation
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, 1.0D0)
            else if (ityp == bufr_radvel) then
              ! Temporary hardcoded value for radar Doppler velocity
              call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, 1.0D0)
            else
              write(*,*) 'varnum = ', ityp
              call utl_abort('oer_fillObsErrors: unknown varnum for RA family')
            end if

          else

            write(*,*) 'oer_fillObsErrors: UNKNOWN DATA FAMILY:',cfam

          end if

          !***********************************************************************
          !              Check for case where error should have been set but was
          !              not. 3dvar will abort in this case.
          !***********************************************************************

          if (obs_bodyElem_r(obsSpaceData, OBS_OER, bodyIndex)  <=  0.0D0) then

            write(*,*)'  PROBLEM OBSERR VARIANCE FAM= ',cfam

            write(*,'(1X,"STNID= ",A10,"codeType= ",I5," LAT= ",F10.2," LON = ",F10.2)') &
                 cstnid,    &
                 codeType,                                           &
                 zlat*MPC_DEGREES_PER_RADIAN_R8,                   &
                 zlon*MPC_DEGREES_PER_RADIAN_R8

            obsPPP = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
            obsOER = obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex)
            write(*,'(1X,"ELEMENT= ",I6," LEVEL= ",F10.2," OBSERR = ",E10.2)')         &
                 ityp, obsPPP, obsOER

            call utl_abort('oer_fillObsErrors: PROBLEM OBSERR VARIANCE.')

          end if

        end if ! end of iass == obs_assimilated

      end do BODY

    end do HEADER

    write(*,'(10X,"Fill_obs_errors: Done")')
    write(*,'(10X,"---------------",/)')

  contains

    function calcStateDepObsErr(clwThresh1,clwThresh2,sigmaThresh1,sigmaThresh2,clw_avg) result(sigmaObsErrUsed)
      implicit none
      real(8) :: clwThresh1
      real(8) :: clwThresh2
      real(8) :: sigmaThresh1
      real(8) :: sigmaThresh2
      real(8) :: clw_avg
      real(8) :: sigmaObsErrUsed

      if (clw_avg <= clwThresh1) then
        sigmaObsErrUsed = sigmaThresh1
      else if (clw_avg >  clwThresh1 .and. & 
                    clw_avg <= clwThresh2) then
        sigmaObsErrUsed = sigmaThresh1 + &
                        (sigmaThresh2 - sigmaThresh1) / &
                        (clwThresh2 - clwThresh1) * &
                        (clw_avg - clwThresh1) 
      else
        sigmaObsErrUsed = sigmaThresh2
      end if

    end function calcStateDepObsErr

  end subroutine oer_fillObsErrors

  subroutine oer_computeInflatedStateDepSigmaObs(obsSpaceData, headerIndex, bodyIndex, &
                                                 ompOmaObsColumn, beSilent_opt)
    !
    ! :Purpose: Update OBS_OER with inflated state dependant observation error
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData
    integer,          intent(in)    :: headerIndex
    integer,          intent(in)    :: bodyIndex
    integer,          intent(in)    :: ompOmaObsColumn  ! obsSpaceData OBS_OMP or OBS_OMA column
    logical, intent(in), optional   :: beSilent_opt

    ! Locals:
    integer :: channelNumber_withOffset
    integer :: channelNumber, channelIndex
    integer :: tovsIndex, sensorIndex
    logical :: surfTypeIsWater 
    real(8) :: clwObs
    real(8) :: clwFG
    real(8) :: deltaE1
    real(8) :: deltaE2
    real(8) :: ompValue
    real(8) :: sigmaObsBeforeInflation
    real(8) :: sigmaObsAfterInflation
    logical :: beSilent

    if (present(beSilent_opt)) then
      beSilent = beSilent_opt
    else
      beSilent = .true.
    end if

    tovsIndex = tvs_tovsIndex(headerIndex)
    sensorIndex = tvs_lsensor(tovsIndex)

    call tvs_getChannelNumIndexFromPPP(obsSpaceData, headerIndex, bodyIndex, &
                                        channelNumber, channelIndex)
    channelNumber_withOffset = channelNumber + tvs_channelOffset(sensorIndex)

    surfTypeIsWater = (obs_headElem_i(obsSpaceData, OBS_STYP, headerIndex) == surftype_sea)

    if (.not. tvs_mwAllskyAssim .or. &
         .not. useStateDepSigmaObs(channelNumber_withOffset,sensorIndex) .or. &
         .not. surfTypeIsWater .or. &
         (.not. mwAllskyInflateByOmp .and. .not. mwAllskyInflateByClwDiff)) return

    if (.not. beSilent) then
      write(*,*) 'oer_computeInflatedStateDepSigmaObs: headerIndex=', headerIndex, &
                        ', sensorIndex=', sensorIndex, &
                        ', chan_noOff=', channelNumber_withOffset, &
                        ', chan_no=', channelNumber
    end if

    clwObs  = obs_headElem_r(obsSpaceData, OBS_CLWO, headerIndex)
    clwFG  = obs_headElem_r(obsSpaceData, OBS_CLWB, headerIndex)

    sigmaObsBeforeInflation = obs_bodyElem_r(obsSpaceData, OBS_OER, bodyIndex)
    ompValue                = obs_bodyElem_r(obsSpaceData, ompOmaObsColumn, bodyIndex)

    ! error inflation for cloud placement 
    deltaE1 = 0.0D0
    if (mwAllskyInflateByOmp .and. &
         ((clwObs - clearClwThresholdSigmaObsInflation(channelNumber_withOffset,sensorIndex)) * &
          (clwFG  - clearClwThresholdSigmaObsInflation(channelNumber_withOffset,sensorIndex)) < 0) .and. &
         abs(clwObs - clwFG) >= 0.005) then
      deltaE1 = abs(ompValue)
    end if

    ! error inflation due to cloud liquid water difference
    deltaE2 = 0.0D0
    if (mwAllskyInflateByClwDiff) then
      deltaE2 = stateDepSigmaObsInflationCoeff(sensorIndex) * abs(clwObs - clwFG) * &
                      sigmaObsBeforeInflation
    end if
    deltaE2 = min(deltaE2,3.5D0 * sigmaObsBeforeInflation)

    sigmaObsAfterInflation = sqrt(sigmaObsBeforeInflation ** 2 + &
                              (deltaE1 + deltaE2) ** 2)

    if (.not. beSilent) then
      write(*,*) 'oer_computeInflatedStateDepSigmaObs: clwObs=', clwObs, &
                        ', clwFG=', clwFG, ', OMP=', ompValue
      write(*,*) 'oer_computeInflatedStateDepSigmaObs: deltaE1=', deltaE1, &
                        ', deltaE2=', deltaE2
      write(*,*) 'oer_computeInflatedStateDepSigmaObs: sigmaObs=', &
          sigmaObsBeforeInflation, ', sigmaObsInflated=', sigmaObsAfterInflation
    end if

    call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, sigmaObsAfterInflation)

  end subroutine oer_computeInflatedStateDepSigmaObs

  !--------------------------------------------------------------------------
  ! readOerFromObsFileForSW
  !--------------------------------------------------------------------------
  subroutine readOerFromObsFileForSW(obsSpaceData)
    !
    ! :Purpose: Read observation errors for AMVs from the obsFiles so as to use
    !           the values computed by a previous MIDAS execution (e.g. GDPS).
    !
    implicit none
    ! arguments
    type(struct_obs) :: obsSpaceData

    ! locals
    character(len=1060) :: filename
    type(BURP_FILE)  :: fileIn
    type(BURP_BLOCK) :: blkoer
    type(BURP_RPT)   :: report
    character(len=9)      :: stnid
    integer(kind=int_def) :: error, ref_rpt
    integer  :: numLevels, numValues, numReports, obsCount
    logical  :: fileFound
    integer  :: levelIndex, reportIndex, obsIndex
    integer  :: uuIndex, vvIndex, headerIndex, bodyIndex, blockIndex, g_btyp_oer
    integer  :: bodyIndexBeg, bodyIndexEnd
    real(8), allocatable :: uu_oer(:), vv_oer(:)

    filename = obsf_getFileName('SW',fileFound)
    if (fileFound) then
      write(*,*) 'oer_readOerFromObsFileForSW: reading OER from the file: ', trim(filename)
    else
      write(*,*) 'oer_readOerFromObsFileForSW: no obsfile with SW family, returning'
      return
    end if

    g_btyp_oer  = 1134

    call burp_init(report)
    call burp_init(blkoer)
    call burp_init(fileIn)
    call burp_new(fileIn, FILENAME=filename, MODE=FILE_ACC_READ, IOSTAT=error)

    !
    ! Scans reports
    !
    call burp_get_property(fileIn, NRPTS = numReports, FILENAME = filename, IOSTAT = error)
    ref_rpt  = 0
    obsCount = 0

    do reportIndex = 1, numReports
      ref_rpt = burp_find_report(fileIn, REPORT=report, SEARCH_FROM=ref_rpt, IOSTAT=error)
      call burp_get_property(report, ELEV=numLevels, IOSTAT=error)
      obsCount = obsCount + numLevels
    end do

    write(*,*) 'readOerFromObsFileForSW: numReports, obsCount  = ',numReports, obsCount

    allocate (uu_oer(obsCount))
    allocate (vv_oer(obsCount))

    !
    ! Read observation errors from file
    !
    ref_rpt  = 0
    obsIndex  = 0

    records_in: do

      ref_rpt = burp_find_report(fileIn, REPORT=report, SEARCH_FROM=ref_rpt, IOSTAT=error)
      if (ref_rpt <0) exit
      
      call burp_get_property(report, STNID=stnid,ELEV=numLevels, IOSTAT=error)

      if (stnid(1:2) == ">>") cycle records_in

      blockIndex = 0
      blockIndex = burp_find_block(report, block=blkoer, search_from=blockIndex, btyp=g_btyp_oer, bfam=10, convert=.false., iostat=error)    

      call burp_get_property(blkoer, NVAL=numValues, IOSTAT=error) 

      uuIndex  = burp_find_element(blkoer, ELEMENT=BUFR_NEUU, IOSTAT=error) 
      vvIndex  = burp_find_element(blkoer, ELEMENT=BUFR_NEVV, IOSTAT=error) 

      if (uuIndex == -1) then
        write(*,*) 'readOerFromObsFileForSW: WARNING: wind element not found, skipping report'
        cycle records_in
      end if
      
      do levelIndex = 1, numLevels

        obsIndex = obsIndex + 1
        if (obsIndex > obsCount) call abort('readOerFromObsFileForSW: Something went wrong')
        uu_oer(obsIndex) = (burp_get_tblval(blkoer, NELE_IND=uuIndex, NVAL_IND=numValues, NT_IND=levelIndex, IOSTAT=error) - 4096)/10.0
        vv_oer(obsIndex) = (burp_get_tblval(blkoer, NELE_IND=vvIndex, NVAL_IND=numValues, NT_IND=levelIndex, IOSTAT=error) - 4096)/10.0

      end do

    end do records_in

    !
    ! Clean-up
    !
    call burp_free(report)
    call burp_free(blkoer)
    call burp_free(fileIn)

    !
    ! Fill obsSpaceData with observation errors from SW file
    !
    obsIndex = 0
    HEADER_UU: do headerIndex = 1, obs_numheader(obsSpaceData)
      if (obs_getFamily(obsSpaceData,headerIndex) /= 'SW') cycle HEADER_UU
      bodyIndexBeg   = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
      bodyIndexEnd = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + bodyIndexBeg - 1
      do bodyIndex = bodyIndexBeg, bodyIndexEnd
        if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) == BUFR_NEUU) then
          obsIndex = obsIndex + 1
          if (obsIndex > obsCount) call abort('readOerFromObsFileForSW: Something went wrong')
          call obs_bodySet_r(obsSpaceData,OBS_OER,bodyIndex,uu_oer(obsIndex))
        end if
      end do
    end do HEADER_UU

    obsIndex = 0
    HEADER_VV: do headerIndex = 1, obs_numheader(obsSpaceData)
      if (obs_getFamily(obsSpaceData,headerIndex) /= 'SW') cycle HEADER_VV
      bodyIndexBeg   = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
      bodyIndexEnd = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + bodyIndexBeg - 1
      do bodyIndex = bodyIndexBeg, bodyIndexEnd
        if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) == BUFR_NEVV) then
          obsIndex = obsIndex + 1
          if (obsIndex > obsCount) call abort('readOerFromObsFileForSW: Something went wrong')
          call obs_bodySet_r(obsSpaceData,OBS_OER,bodyIndex,vv_oer(obsIndex))
        end if
      end do
    end do HEADER_VV

    deallocate (uu_oer)
    deallocate (vv_oer)

  end subroutine readOerFromObsFileForSW

  !--------------------------------------------------------------------------
  ! oer_sw
  !--------------------------------------------------------------------------
  subroutine oer_sw(columnTrlOnTrlLev,obsSpaceData)
    !
    ! :Purpose: Calculate observation errors for AMVs according to the Met-Office 
    !           situation dependant approach.
    !
    implicit none

    type(struct_columnData) :: columnTrlOnTrlLev
    type(struct_obs) :: obsSpaceData

    integer :: headerIndex,bodyIndex,ilyr,jlev
    integer :: iass,ixtr,ivco,ivnm,iqiv,iqiv1,iqiv2,imet,ilsv,igav,ihav,itrn,J_SAT
    real(8) :: zvar,zoer
    real(8) :: zwb,zwt,ZOTR,ZMOD
    real(8) :: zlat,zlon,zlev,zpt,zpb,zpc
    real(8) :: SP_WGH,TO_WGH,TO_DSP,E_VHGT,E_DRIFT,E_HEIGHT
    character(len=4) :: varName
    character(len=4) :: varLevel
    character(len=9) :: cstnid
    real(8), pointer :: col_ptr_uv(:)
    logical :: passe_once, valeurs_defaut, print_debug
    logical, save :: firstCall=.true.

    ! If requested, just read oer from the burp file (only 1st time)
    if(obsfile_oer_sw) then
      if (firstCall) then
        call readOerFromObsFileForSW(obsSpaceData)
        firstCall = .false.
      end if
      return
    end if
    
    if(.not. new_oer_sw) return

    valeurs_defaut = .false.
    passe_once     = .true.
    print_debug    = .false.

    if (firstCall) write(*,*) "Entering subroutine oer_sw"
    firstCall = .false.

    call obs_set_current_body_list(obsSpaceData, 'SW')
    BODY: do
      bodyIndex = obs_getBodyIndex(obsSpaceData)
      if (bodyIndex < 0) exit BODY

      ! Only process pressure level observations flagged to be assimilated
      iass=obs_bodyElem_i (obsSpaceData,OBS_ASS,bodyIndex)
      ivco=obs_bodyElem_i (obsSpaceData,OBS_VCO,bodyIndex)
      if(iass /= obs_assimilated .or. ivco /= 2) cycle BODY

      ixtr = obs_bodyElem_i (obsSpaceData,OBS_XTR,bodyIndex)
      ivnm = obs_bodyElem_i (obsSpaceData,OBS_VNM,bodyIndex)
      zvar = obs_bodyElem_r (obsSpaceData,OBS_VAR,bodyIndex)
      zlev = obs_bodyElem_r (obsSpaceData,OBS_PPP,bodyIndex)
      zoer = obs_bodyElem_r (obsSpaceData,OBS_OER,bodyIndex)
      headerIndex = obs_bodyElem_i (obsSpaceData,OBS_HIND,bodyIndex)

      iqiv1 = obs_headElem_i (obsSpaceData,OBS_SWQ1,headerIndex)
      iqiv2 = obs_headElem_i (obsSpaceData,OBS_SWQ2,headerIndex)
      imet =  obs_headElem_i (obsSpaceData,OBS_SWMT,headerIndex)
      ilsv =  obs_headElem_i (obsSpaceData,OBS_SWLS,headerIndex)
      igav =  obs_headElem_i (obsSpaceData,OBS_SWGA,headerIndex)
      ihav =  obs_headElem_i (obsSpaceData,OBS_SWHA,headerIndex)
      zlat =  obs_headElem_r (obsSpaceData,OBS_LAT,headerIndex)*MPC_DEGREES_PER_RADIAN_R8
      zlon =  obs_headElem_r (obsSpaceData,OBS_LON,headerIndex)*MPC_DEGREES_PER_RADIAN_R8
      cstnid = obs_elem_c (obsSpaceData,'STID' ,headerIndex)

      itrn = ilsv
      if(igav == 1) itrn = 2
      if(imet >  3) imet = 3
      if(ihav /= 4) ihav = 1
      ! Use the quality score qiv2, but if it is missing then use qiv1
      iqiv = iqiv2
      if(iqiv < 0) iqiv = iqiv1

      if(valeurs_defaut) then
        E_DRIFT  = 2.5
        E_HEIGHT = 8000.0
      else
        E_DRIFT = 7.5 - 0.05*iqiv
        call get_height_error(cstnid,imet,itrn,ihav,zlat,zlev,E_HEIGHT,J_SAT)
      end if

      varLevel = vnl_varLevelFromVarnum(ivnm)
      varName  = vnl_varnameFromVarnum(ivnm)

      col_ptr_uv=>col_getColumn(columnTrlOnTrlLev,headerIndex,varName)
      ilyr = obs_bodyElem_i (obsSpaceData,OBS_LYR,bodyIndex)
      ZPT  = col_getPressure(columnTrlOnTrlLev,ilyr  ,headerIndex,varLevel)
      ZPB  = col_getPressure(columnTrlOnTrlLev,ilyr+1,headerIndex,varLevel)
      ZWB  = LOG(zlev/ZPT)/LOG(ZPB/ZPT)
      ZWT  = 1. - ZWB

      ZMOD = ZWB*col_ptr_uv(ilyr+1) + ZWT*col_ptr_uv(ilyr)

      TO_WGH = 0.0D0
      TO_DSP = 0.0D0
      if(zlev > 20000. .and. zlev < 30000. .and. passe_once .and. print_debug) then
        write(*,'(3a10,2i10,4f12.3)') 'stlchka',cstnid,varName,iqiv,imet,zlev,ZMOD,ZVAR,E_DRIFT
      end if
      do jlev = 2, col_getNumLev(columnTrlOnTrlLev,'MM') - 1
        ZPT= col_getPressure(columnTrlOnTrlLev,jlev-1,headerIndex,varLevel)
        ZPC= col_getPressure(columnTrlOnTrlLev,jlev  ,headerIndex,varLevel)
        ZPB= col_getPressure(columnTrlOnTrlLev,jlev+1,headerIndex,varLevel)
        ZOTR = col_ptr_uv(jlev)
        SP_WGH = exp(-0.5*((ZPC - zlev)**2)/(E_HEIGHT**2))*((ZPB - ZPT)/2)
        TO_DSP = TO_DSP + SP_WGH*((ZOTR - ZMOD)**2)
        TO_WGH = TO_WGH + SP_WGH 
        if(zlev > 20000. .and. zlev < 30000. .and. passe_once .and. print_debug) then
          write(*,'(a10,i10,4f12.3)') 'stlchk',jlev,ZPT,ZPC,ZPB,ZOTR
        end if

      end do

      E_VHGT = sqrt(TO_DSP/TO_WGH)
      zoer = sqrt(E_VHGT**2 + E_DRIFT**2)
       
      if(zlev > 20000. .and. zlev < 30000. .and. passe_once .and. print_debug) then
        write(*,'(a10,4f10.2)') 'stlchkb',zoer,E_VHGT,E_DRIFT,E_HEIGHT
        passe_once = .false.
      end if

      call obs_bodySet_r(obsSpaceData,OBS_OER,bodyIndex,zoer)

      if(print_debug) write(*,'(2a10,6f12.3,4i10)') 'hgterr',cstnid,zlat,zlon,zlev/100.,E_HEIGHT/100.0,E_DRIFT,zoer,imet,itrn,ihav,J_SAT
      
    end do BODY

  end subroutine oer_sw

  !--------------------------------------------------------------------------
  ! get_height_error
  !--------------------------------------------------------------------------
  subroutine get_height_error(stnid, methode, terrain, htasmet, zlat, zlev, E_HEIGHT, J_SAT)
    character(len=9), intent(in)     :: stnid
    integer,          intent(in)     :: methode
    integer,          intent(in)     :: terrain
    integer,          intent(in)     :: htasmet
    integer,          intent(out)    :: J_SAT 
    real(8),          intent(in)     :: zlat
    real(8),          intent(in)     :: zlev
    real(8),          intent(out)    :: E_HEIGHT

    integer, parameter :: NLVLVALUE=9

    integer :: I_HGT, jelm, jlev, jcat, i_typ, hemisphere
    real(8) :: zlev_hpa, ZWB, ZWT

    logical :: interpole

    if(zlat >= 0) hemisphere = 1
    if(zlat <  0) hemisphere = 2


    J_SAT = 1 ! default value

    i_typ = 0
    list1: do jlev = 1,n_sat_type
      do jelm = 2, 10
        if(trim(stnid(2:9)) == trim(SAT_AMV(jlev,jelm))) then
          i_typ = jlev
          exit list1
        end if
      end do
    end do list1

    if (i_typ /= 0) then
      list2:   do jcat = 1,n_categorie
        if(trim(SAT_AMV(i_typ,1)) == trim(SAT_LIST(jcat))) then
          if(tbl_m(jcat) == 0 .or. tbl_m(jcat) == methode) then
            if(tbl_h(jcat) == 0 .or. tbl_h(jcat) == htasmet) then
              if(tbl_t(jcat) == 0 .or. tbl_t(jcat) == terrain+1) then
                if(tbl_g(jcat) == 0 .or. tbl_g(jcat) == hemisphere) then
                  J_SAT = jcat
                  exit list2
                end if
              end if
            end if
          end if
        end if
      end do list2
    end if

    zlev_hpa  = zlev/100.
    interpole = .true.
    do I_HGT=1,NLVLVALUE
      if(zlev_hpa >= LVLVALUE(I_HGT)) exit
    end do
    if(I_HGT == 1)            interpole = .false.
    if(I_HGT == NLVLVALUE + 1) interpole = .false.
    if(I_HGT == NLVLVALUE + 1) I_HGT     = NLVLVALUE

    if(interpole) then

      ZWB      = LOG(zlev_hpa/LVLVALUE(I_HGT-1))/LOG(LVLVALUE(I_HGT)/LVLVALUE(I_HGT-1))
      ZWT      = 1. - ZWB
      E_HEIGHT = ZWT*HGT_ERR(J_SAT,I_HGT-1) + ZWB*HGT_ERR(J_SAT,I_HGT)

    else

      E_HEIGHT = HGT_ERR(J_SAT,I_HGT)
      
    end if

!
!   Retourne l'erreur en (Pa)
!
    E_HEIGHT = E_HEIGHT*100.0

  end subroutine get_height_error

  !--------------------------------------------------------------------------
  ! oer_SETERRGPSRO
  !--------------------------------------------------------------------------
  SUBROUTINE oer_SETERRGPSRO(columnTrlOnTrlLev, obsSpaceData, beSilent)
    !
    ! :Purpose: Compute estimated errors for GPSRO observations
    !
    IMPLICIT NONE
    !
    type(struct_columnData) :: columnTrlOnTrlLev
    type(struct_obs)        :: obsSpaceData
    logical                 :: beSilent
    !
    integer headerIndex, IDATYP, bodyIndex, iProfile, varNum
    REAL*8 zLat, Lat, sLat
    REAL*8 zLon, Lon
    REAL*8 zAzm
    REAL*8, allocatable :: ZPP(:)
    REAL*8, allocatable :: ZTT(:)
    REAL*8, allocatable :: ZHU(:)
    REAL*8, allocatable :: zHeight(:)
    REAL*8, allocatable :: ZUU(:)
    REAL*8, allocatable :: ZVV(:)
    !
    REAL*8 DH,DDH
    integer JL, isat, JH, NGPSLEV, NWNDLEV
    REAL*8 zMT, Rad, Geo, zP0
    REAL*8 HNH1, HJH, SUM0, SUM1, ZMIN, WFGPS, H1, F2, F3, F4
    !
    LOGICAL  ASSIM, L1, L2, L3

    integer NH, NH1
    TYPE(GPS_PROFILE)           :: PRF
    REAL*8       , allocatable :: H   (:),AZMV(:)
    REAL*8       , allocatable :: ZOBS(:),ZREF(:),ZOFF(:),ZERR(:), ZMHX(:)
    TYPE(GPS_DIFF), allocatable :: RSTV(:)

    if (.not. beSilent) write(*,*)'ENTER SETERRGPSRO'
    !
    !     * 1.  Initializations
    !     *     ---------------
    !
    NGPSLEV=col_getNumLev(columnTrlOnTrlLev,'TH')
    NWNDLEV=col_getNumLev(columnTrlOnTrlLev,'MM')
    allocate(ZPP (NGPSLEV))
    allocate(ZTT (NGPSLEV))
    allocate(ZHU (NGPSLEV))
    allocate(zHeight (NGPSLEV))
    allocate(ZUU (NGPSLEV))
    allocate(ZVV (NGPSLEV))
    !
    allocate(H    (gps_RO_MAXPRFSIZE))
    allocate(AZMV (gps_RO_MAXPRFSIZE))
    allocate(ZOBS (gps_RO_MAXPRFSIZE))
    allocate(ZREF (gps_RO_MAXPRFSIZE))
    allocate(ZOFF (gps_RO_MAXPRFSIZE))
    allocate(ZERR (gps_RO_MAXPRFSIZE))
    allocate(RSTV (gps_RO_MAXPRFSIZE))
    allocate(ZMHX (gps_RO_MAXPRFSIZE))
    !
    !     Loop over all header indices of the 'RO' family:
    !
    call obs_set_current_header_list(obsSpaceData,'RO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER
       !
       !     *  Process only refractivity data (codtyp 169)
       !
      IDATYP = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
      IF (IDATYP == 169) then
          !
          !     *     Scan for requested data values of the profile, and count them
          !
        ASSIM = .FALSE.
        NH = 0
        call obs_set_current_body_list(obsSpaceData, headerIndex)
        BODY: do 
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY
          IF (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) == obs_assimilated) then
            ASSIM = .TRUE.
            NH = NH + 1
          end if
        end do BODY
          !
          !     *     If assimilations are requested, prepare and apply the observation operator
          !
        IF (ASSIM) then
          iProfile=gps_iprofile_from_index(headerIndex)
          varNum = gps_vRO_IndexPrf(iProfile, 2)
             !
             !     *        Basic geometric variables of the profile:
             !
          isat = obs_headElem_i(obsSpaceData,OBS_SAT,headerIndex)
          Rad  = obs_headElem_r(obsSpaceData,OBS_TRAD,headerIndex)
          Geo  = obs_headElem_r(obsSpaceData,OBS_GEOI,headerIndex)
          zAzm = obs_headElem_r(obsSpaceData,OBS_AZA,headerIndex) / MPC_DEGREES_PER_RADIAN_R8
          zMT  = col_getHeight(columnTrlOnTrlLev,0,headerIndex,'SF')
             !     
             !     *        Profile at the observation location:
             !
          zLat = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)
          zLon = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex)
          Lat  = zLat * MPC_DEGREES_PER_RADIAN_R8
          Lon  = zLon * MPC_DEGREES_PER_RADIAN_R8
          sLat = sin(zLat)
          zMT  = zMT * ec_rg / gps_gravitysrf(sLat)
          zP0  = col_getElem(columnTrlOnTrlLev,1,headerIndex,'P0')
          DO JL = 1, NGPSLEV
                !
                !     *           Profile x
                !
            ZPP(JL) = col_getPressure(columnTrlOnTrlLev,JL,headerIndex,'TH')
            ZTT(JL) = col_getElem(columnTrlOnTrlLev,JL,headerIndex,'TT') - MPC_K_C_DEGREE_OFFSET_R8
            ZHU(JL) = col_getElem(columnTrlOnTrlLev,JL,headerIndex,'HU')
            ZUU(JL) = 0.d0
            ZVV(JL) = 0.d0
            zHeight(jl) = col_getHeight(columnTrlOnTrlLev,jl,headerIndex,'TH')
          end do

          if((col_getPressure(columnTrlOnTrlLev,1,headerIndex,'TH') + 1.0d-4)  <  &
               col_getPressure(columnTrlOnTrlLev,1,headerIndex,'MM')) then
                ! case with top thermo level above top momentum level (Vcode=5002)
            do jl = 1, nwndlev
              zuu(jl) = col_getElem(columnTrlOnTrlLev,jl, headerIndex, 'UU')
              zvv(jl) = col_getElem(columnTrlOnTrlLev,jl, headerIndex, 'VV')
            end do
          else
                ! case without top thermo above top momentum level or unstaggered (Vcode=5001/4/5)
            do jl = 1, nwndlev - 1
              zuu(jl) = col_getElem(columnTrlOnTrlLev, jl + 1, headerIndex, 'UU')
              zvv(jl) = col_getElem(columnTrlOnTrlLev, jl + 1, headerIndex, 'VV')
            end do
            zuu(nwndlev) = zuu(nwndlev - 1)
            zvv(nwndlev) = zuu(nwndlev - 1)
          end if
          zuu(ngpslev) = zuu(nwndlev)
          zvv(ngpslev) = zuu(nwndlev)
             !     
             !     *        GPS profile structure:
             !
          call gps_struct1sw_v2(ngpslev,zLat,zLon,zAzm,zMT,Rad,geo,zP0,zPP,zTT,zHU,zUU,zVV,zHeight,prf)
             !
             !     *        Prepare the vector of all the observations:
             !
          NH1 = 0
             !
             !     *        Loop over all body indices for this index_header:
             !     *        (start at the beginning of the list)
             !
          call obs_set_current_body_list(obsSpaceData, headerIndex)
          BODY_2: do 
            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if (bodyIndex < 0) exit BODY_2
            IF (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) == obs_assimilated) then
              NH1      = NH1 + 1
              H(NH1)   = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex)
              AZMV(NH1)= zAzm
              ZOBS(NH1)= obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex)
                   !     *              Reference value:
              IF (varNum == bufr_nebd) then
                ZREF(NH1) = 0.025d0*exp(-(H(NH1)-Rad)/6500.d0)
              ELSE
                ZREF(NH1) = 300.d0*exp(-H(NH1)/6500.d0)
              end if
            end if
          end do BODY_2
             !
             !     *        Apply the observation operator:
             !
          IF (varNum == bufr_nebd) then
            call GPS_BNDOPV1(H, AZMV, NH, PRF, RSTV)
          ELSE
            call GPS_REFOPV (H,       NH, PRF, RSTV)
          end if
             !
             !     *        Perform the (H(x)-Y)/R operation:
             !
          DO NH1 = 1, NH
            ZMHX(NH1) = RSTV(NH1)%VAR
                !
                !     *           Normalized offset:
                !
            IF (.NOT.gps_roBNorm) then
              ZOFF(NH1) = (ZOBS(NH1) - ZMHX(NH1)) / ZREF(NH1)
            ELSE
              ZOFF(NH1) = (ZOBS(NH1) - ZMHX(NH1)) / ZMHX(NH1)
            end if
          end do
             !
             !     *        The procedure below is well tested to collectively
             !     *        create error profiles from the data profile, and
             !     *        intended to be used for these data.
             !
          DH = 5000.d0
          if (varNum == bufr_nebd) then
            ZMIN=0.01D0
          else
            ZMIN=0.002D0
          end if

          if (varNum == bufr_nerf) then
            if (trim(gps_roError) == 'DYNAMIC') then
              do NH1 = 1, NH
                SUM0=1.d-30
                SUM1=0.d0
                do JH = 1, NH
                  if (H(JH) <= gps_HtpMaxEr) then
                    DDH=H(JH)-H(NH1)
                    SUM0=SUM0+EXP(-(DDH/DH)**2)
                    SUM1=SUM1+EXP(-(DDH/DH)**2)*ZOFF(JH)**2
                  end if
                end do
                ZERR(NH1)=(SUM1/SUM0)**0.5D0
                if (ZERR(NH1) < ZMIN) ZERR(NH1) = ZMIN
              end do
            else if (trim(gps_roError) == 'STATIC_2018') then
              ! this was introduced by Maziar in late 2018 on advice by Josep
              do NH1 = 1, NH
                HNH1 = H(NH1)
                ZERR(NH1) = 0.05d0
                L1=(HNH1 <= 10000.d0)
                L2=(HNH1 > 10000.d0 .and. HNH1 < 30000.d0)
                L3=(HNH1 > 30000.d0)
                IF (L1) ZERR(NH1)=0.005d0+0.020d0*(10000.d0-HNH1)/10000.d0
                IF (L2) ZERR(NH1)=0.005d0
                IF (L3) ZERR(NH1)=0.005d0+0.030d0*(HNH1-30000.d0)/30000.d0
                if (ZERR(NH1) < ZMIN) ZERR(NH1) = ZMIN
              end do
            else if (trim(gps_roError) == 'STATIC_2014') then
              ! recipe used in EnKF from Josep by email on February 25 2014 
              do NH1 = 1, NH
                HNH1 = H(NH1)
                select case (nint(hnh1))
                case(:10000)
                  ZERR(NH1) = 0.005 + 0.015 * ((10000.0-hnh1)/10000.0)
                case (10001:30000)
                  ZERR(NH1) = 0.005
                case (30001:)
                  ZERR(NH1) = 0.005 + 0.010 * ((hnh1-30000.0)/10000.0)
                end select
              end do
            else
              call utl_abort('oer_setErrGPSro: Invalid value for gps_roError')
            endif
          else
            if (trim(gps_roError) == 'DYNAMIC') then
              ZMIN = 0.01D0
              do NH1 = 1, NH
                HNH1=H(NH1)-Rad
                SUM0=1.d-30
                SUM1=0.d0
                DH = 1000.d0 + 0.1d0 * HNH1
                do JH = 1, NH
                  HJH=H(JH)-Rad
                  if (HJH <= gps_HtpMaxEr) then
                    DDH=HJH-HNH1
                    SUM0=SUM0+EXP(-(DDH/DH)**2+(DDH/DH))
                    SUM1=SUM1+EXP(-(DDH/DH)**2+(DDH/DH))*ZOFF(JH)**2
                  end if
                end do
                ZERR(NH1)=(SUM1/SUM0)**0.5D0
                if (ZERR(NH1) < ZMIN) ZERR(NH1) = ZMIN
                if (H(NH1) < PRF%ast(ngpslev)%Var) ZERR(NH1)=0.08
              end do
            else if (trim(gps_roError) == 'STATIC_2018') then
              do NH1 = 1, NH
                ZERR(NH1)=0.05d0
                HNH1=H(NH1)-Rad
                L1=(HNH1 <= 10000.d0)
                L2=(HNH1 > 10000.d0 .and. HNH1 < 30000.d0)
                L3=(HNH1 > 30000.d0)
                IF (L1) ZERR(NH1)=0.02d0+0.08d0*(10000.d0-HNH1)/10000.d0
                IF (L2) ZERR(NH1)=0.02d0
                IF (L3) ZERR(NH1)=0.02d0+0.13d0*(HNH1-30000.d0)/30000.d0
                IF (isat==3 .or. isat==4 .or. isat==5) ZERR(NH1) = 2*ZERR(NH1)
                IF (ZERR(NH1) < ZMIN) ZERR(NH1) = ZMIN
              end do
            end if
          end if

          NH1 = 0
             !
             !     *        Loop over all body indices for this index_header:
             !     *        (start at the beginning of the list)
             !
          call obs_set_current_body_list(obsSpaceData, headerIndex)
          BODY_4: do 
            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if (bodyIndex < 0) exit BODY_4
            IF (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) == obs_assimilated) then
              NH1 = NH1 + 1
              H1 = H(NH1)
              if (varNum == bufr_nebd) then
                 H1 = H1 - Rad
                 F2 = 0.5d0*(erf((H1-22000.d0)/5000.d0)+1.d0)
              else
                 F2 = exp(0.5d0*H1/6500.d0)
              endif
              F3 = exp(-((H1- 6500.d0)/6500.d0)**2)
              F4 = exp(-((H1-14500.d0)/6500.d0)**2)
              WFGPS = gps_WGPS(ISAT,1) + gps_WGPS(ISAT,2) * F2 + gps_WGPS(ISAT,3) * F3 + gps_WGPS(ISAT,4) * F4
              IF (WFGPS < 1.D0) WFGPS = 1.D0
              WFGPS = sqrt(WFGPS)
              !
                   !     *              Observation error    S
                   !
              if (.NOT.gps_roBNorm) then
                call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, ZERR(NH1)*ZREF(NH1)*WFGPS)
              else
                call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, ZERR(NH1)*ZMHX(NH1)*WFGPS)
              end if
            end if
          end do BODY_4
        end if
      end if
    end do HEADER

    deallocate(RSTV)
    deallocate(ZERR)
    deallocate(ZOFF)
    deallocate(ZREF)
    deallocate(ZOBS)
    deallocate(AZMV)
    deallocate(H)
    deallocate(ZMHX)

    deallocate(zVV)
    deallocate(zUU)
    deallocate(zHU)
    deallocate(zHeight)
    deallocate(zTT)
    deallocate(zPP)

    if (.not.beSilent) write(*,*)'oer_setErrGPSRO: done'

  END SUBROUTINE OER_SETERRGPSRO

  !--------------------------------------------------------------------------
  ! oer_SETERRGPSGB
  !--------------------------------------------------------------------------
  SUBROUTINE oer_SETERRGPSGB(columnTrlOnTrlLev, obsSpaceData, beSilent, ldata, analysisMode)
    !
    ! :Purpose:
    !
    !             - Set the observation errors [OBS_OER] and Std(O-P) [OBS_HPHT] for GB-GPS ZTD data.
    !               (GPS surface met obs errors are set before in observation_erreurs_mod.ftn90.
    !               The ZTD error is also initialized to the "formal error" or to 1.0 if missing.)
    !             - Returns ldata=.false. if there are no GPS ZTD data to assimilate
    !               and also sets the modgpsztd_mod variable gps_gb_numZTD = 0.
    !
    IMPLICIT NONE
    !!
    !! NOTE: gps_gb_YZDERRWGT IN modgpsztd_mod (FROM NML FILE) IS USED FOR ERROR WEIGHTING
    !!       OF TIME SERIES (FGAT) GPS ZTD OBSERVATIONS TO ACCOUNT FOR TEMPORAL ERROR
    !!       CORRELATIONS.
    !!
    type(struct_columnData) :: columnTrlOnTrlLev
    type(struct_obs)        :: obsSpaceData
    logical                 :: beSilent
    logical                 :: ldata
    logical                 :: analysisMode

    integer bodyIndex, headerIndex, ityp, iass, IZTDJ, NBRPDATE, ICOUNT, ICOUNT2
    integer nlev_T

    LOGICAL LLCZTDE, LLFER, LLFZTDE, LLZTD, LLRZTDE, ASSIM, ERRSET, DEBUG, LESTP
    LOGICAL LLZWD
    character(len=12) :: cstnid

    !
    REAL*8  ZTDERR, ZZTD, ZMINZDE, ZPSFC, ZHD, ZWD, ZTDOER, zlev, zval, ZZWD
    REAL*8  ZBTSFC, ZBPSFC, ZBZSFC, ZDZ, ZSTDOMP

    !
    !     ZZDERMIN = MIN ZTD OER VALUE (M)
    !     ZZDERMAX = MAX ZTD OER VALUE (M)
    !     ZTDERFAC = MULTIPLICATION FACTOR FOR FORMAL ZTD MEASUREMENT ERROR
    !     ZOPEFAC  = FRACTION OF REGRESSION EQUATION SD(O-P) TO USE AS ZTD OBSERVATION ERROR
    !     ----------------------------------------------------------------------------------
    !
    REAL*8 ZZDERMIN, ZZDERMAX, ZTDERFAC, ZOPEFAC
    DATA ZZDERMIN /0.004D0/
    DATA ZZDERMAX /0.030D0/
    DATA ZTDERFAC /3.0D0/
    DATA ZOPEFAC  /1.0D0/
    !
    !     FOR ESTIMATION OF PSFC (IF MISSING)
    !       ZGAMMA = (NEG. OF) TEMPERATURE LAPSE RATE (K/M)
    !
    REAL*8 ZGAMMA
    DATA ZGAMMA /0.0065D0/

    !     ----------------------------------------------------------------------------------
    !     LINEAR REGRESSION EQUATION CONSTANTS AND COEFFS FOR ZTD ERROR AND STD(O-P):

    !     ZRCONST, ZRCOEFF:
    !       - From linear regression of Desroziers error estimates binned by observed ZWD.
    !       - Gives ZTDerror (mm) as function of ZWD (m):
    !            ZTDerror(mm) = ZRCONST + ZRCOEFF*ZWD(m)
    !     ZRCONST2, ZRCOEFF2:
    !       - From linear regression of Std(O-P) binned by observed ZWD.
    !       - Gives Std(O-P) (mm) as function of ZWD (m):
    !            Std(O-P)(mm) = ZRCONST2 + ZRCOEFF2*ZWD(m)
    !     ----------------------------------------------------------------------------------
    REAL*8 ZRCONST, ZRCOEFF, ZRCONST2, ZRCOEFF2
    DATA  ZRCONST  /5.12D0/
    DATA  ZRCOEFF  /26.4D0/
    DATA  ZRCONST2 /6.67D0/
    DATA  ZRCOEFF2 /42.6D0/
    !
    !
    if (.not.beSilent) write(*,*) 'ENTER SETERRGPSGB'
    !
    DEBUG = .FALSE.

    LLCZTDE = .FALSE.
    LLRZTDE = .FALSE.
    LLFZTDE = .FALSE.
    IF (gps_gb_YZTDERR  <  0.0D0) then
      LLFZTDE = .TRUE.
    else if (gps_gb_YZTDERR  >  0.0D0) then
      LLCZTDE = .TRUE.
    ELSE
      LLRZTDE = .TRUE.
    end if

    nlev_T = col_getNumLev(columnTrlOnTrlLev,'TH')

    ldata = .false.
    ICOUNT  = 0
    ICOUNT2 = 0
    !
    !     Loop over all header indices of the 'GP' family:
    !
    call obs_set_current_header_list(obsSpaceData,'GP')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER
      NBRPDATE  = obs_headElem_i(obsSpaceData, OBS_DAT, headerIndex)
      LLZTD     = .FALSE.
      LLFER     = .FALSE.
      LLZWD     = .FALSE.
      ASSIM     = .FALSE.
      ERRSET    = .FALSE.
      ZZTD      = -1.0D0
      ZZWD      = -1.0D0
      ZPSFC     = -1.0D0
      LESTP     = .FALSE.
      ZSTDOMP   = 15.0D0*0.001D0

       !   Get Psfc (Pa), Tsfc (K) and model surface height (m) from background profile

      ZBPSFC = col_getElem(columnTrlOnTrlLev, 1, headerIndex, 'P0')
      ZBTSFC = col_getElem(columnTrlOnTrlLev, nlev_T, headerIndex, 'TT')
      ZBZSFC = col_getHeight(columnTrlOnTrlLev, nlev_T, headerIndex, 'TH')
       !
       !    Loop over all body indices of current report; Set the ZTD error if
       !    constant value specified (LLCZTDE=true). Get GPS height and Psfc obs (if any).
       !    Get ZTD obs, ZTD formal error and ZWD observation.
       !
      call obs_set_current_body_list(obsSpaceData, headerIndex)
      BODY: do
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if (bodyIndex < 0) exit BODY
        ityp   = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex)
        iass   = obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex)
        zval   = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex)
          !         Store Psfc
        IF (ityp == BUFR_NEPS) then
          IF (zval  >  0.0D0) ZPSFC = zval
        end if
          !         Set ZTDOER to constant value (if LLCZTDE); get value of ZTD, 
          !         ZTD formal error (OBS_OER) and antenna height (OBS_PPP).
        IF (ityp == BUFR_NEZD) then
          IF (LLCZTDE) then
            ZTDOER = gps_gb_YZTDERR
            ERRSET = .TRUE.
          end if
          zlev   = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex)
          ZTDERR = obs_bodyElem_r(obsSpaceData, OBS_OER, bodyIndex)
          IF (ZTDERR  /=  1.0D0) LLFER = .TRUE.
          IZTDJ = bodyIndex
          IF (zval  >  0.0D0) then
            ZZTD = zval
            LLZTD = .TRUE.
          end if
          IF (iass == obs_assimilated) ASSIM = .TRUE.
        end if
        IF (ityp == BUFR_NEZW) then
          IF (zval  >  0.0D0) then
            ZZWD = zval
            LLZWD = .TRUE.
          end if
        end if
      end do BODY

       !      Initialize Std(O-P) to 15 mm  for ZTD observation (for bgck mode)
      IF (LLZTD .and. .NOT.analysisMode) &
           call obs_bodySet_r(obsSpaceData,OBS_HPHT,IZTDJ,ZSTDOMP)

       !      Replace formal ZTD error with real error for all ZTD to be assimilated.
       !      Set Std(O-P) as function of ZWD for ZTD observation and store in OBS_HPHT. 

      if (ASSIM) then
        if (LLZTD) then
          ldata = .true.
          ICOUNT = ICOUNT + 1
          if (LLZWD) then
            ZWD = ZZWD
          else
                !               If Psfc obs is missing, estimate the pressure from background state
            if (ZPSFC  <  0.0D0) then
              LESTP = .TRUE.
              ZDZ = zlev - ZBZSFC
              ZPSFC  = ZBPSFC * &
                   (1.0D0-(ZGAMMA/ZBTSFC)*ZDZ)**(ec_rg/(MPC_RGAS_DRY_AIR_R8*ZGAMMA))
              ICOUNT2 = ICOUNT2 + 1
            end if
                !                Compute the hydrostatic delay ZHD (m) from Psfc (Pa)
            ZHD = 2.2766D-05 * ZPSFC
                !               Compute the wet delay (m) from ZTD and ZHD. Avoid negative ZWD.
            if (ZHD  >  ZZTD) then
              ZWD = 0.0D0
            else
              ZWD = ZZTD - ZHD
            end if
          end if
             !              Std(O-P) for background check. Limit to 30 mm in case ZTD obs is bad (too high).             
          ZSTDOMP = (ZRCONST2 + ZRCOEFF2*ZWD)*0.001D0
          ZSTDOMP = MIN(ZZDERMAX, ZSTDOMP)
             !             Compute ZTD error as a function of ZWD using regression coeff (SD(O-P) vs ZWD).
             !             Take fraction ZOPEFAC of computed error and convert from mm to m.
             !             Ensure error is > ZZDERMIN and < ZZDERMAX
          IF (.NOT. ERRSET) then 
            ZMINZDE = ZRCONST + ZRCOEFF*ZWD
            ZMINZDE = ZMINZDE * ZOPEFAC * 0.001D0
            IF (LLRZTDE) then
              ZTDOER = MAX(ZZDERMIN, ZMINZDE)
              ZTDOER = MIN(ZZDERMAX, ZTDOER)
            ELSE
              IF (LLFER) then
                ZTDOER = MAX(ZZDERMIN, ZTDERR*ZTDERFAC)
              ELSE
                ZTDOER = MAX(ZZDERMIN, ZMINZDE)
                ZTDOER = MIN(ZZDERMAX, ZTDOER)
              end if
            end if
                !  Ensure that error is not less than formal error ZTDERR
            IF (LLFER) then
              IF (DEBUG .and. ICOUNT  <=  50) then
                write(*,*) cstnid, &
                     ' FORMAL ERR, OBS ERROR (mm) = ', &
                     ZTDERR*1000.D0, ZTDOER*1000.D0
              end if
              ZTDOER = MAX(ZTDOER, ZTDERR)
            end if
          end if
             !  *** APPLY TIME-SERIES WEIGHTING FACTOR TO OBSERVATION ERROR (gps_gb_YZDERRWGT=1 FOR 3D THINNING)
          call obs_bodySet_r(obsSpaceData,OBS_OER,IZTDJ, ZTDOER*gps_gb_YZDERRWGT)
          IF (.NOT.analysisMode) call obs_bodySet_r(obsSpaceData,OBS_HPHT,IZTDJ, ZSTDOMP)
          IF (DEBUG .and. (ICOUNT2  <=  50) .and. LESTP) then
            write(*,*) 'TAG    SITE    ZTD    ERROR    ELEV    PSFC    ZWD     STDOMP'
            write(*,*) 'ERRDEBUG ', cstnid, &
                 ZZTD*1000.D0, ZTDOER*1000.D0, zlev, ZPSFC/100.D0, ZWD*1000.D0, ZSTDOMP*1000.D0
          end if
        ELSE
          call utl_abort('SETERRGPSGB: ERROR:NEGATIVE ZTD VALUE!')
        end if
      end if

       !
    end do  HEADER

    !      IF (DEBUG) call utl_abort('******DEBUG STOP*******')

    IF (.not.ldata) gps_gb_numZTD = 0

    IF (ldata .and. .not.beSilent) write(*,*) ' gps_gb_numZTD = ', ICOUNT

    if (.not.beSilent) write(*,*) 'EXIT SETERRGPSGB'

  END SUBROUTINE OER_SETERRGPSGB

  !--------------------------------------------------------------------------
  ! oer_setErrBackScatAnisIce
  !--------------------------------------------------------------------------
  subroutine oer_setErrBackScatAnisIce(columnTrlOnTrlLev, obsSpaceData, beSilent)
    !
    ! :Purpose: Compute estimated errors for ASCAT backscatter anisotropy observations
    !
    implicit none

    type(struct_columnData), intent(in) :: columnTrlOnTrlLev
    type(struct_obs)                    :: obsSpaceData
    logical,                 intent(in) :: beSilent

    ! locals
    integer :: headerIndex, bodyIndex
    integer :: idate, imonth, varno
    integer :: trackCellNum

    real(8) :: conc, obsErrStdDev

    character(len=*), parameter :: myName = 'oer_setErrBackScatAnisIce'
    character(len=8)            :: ccyymmdd

    if (.not. beSilent) write(*,*) 'ENTER '//myName

    !
    !     Loop over all header indices of the 'GL' family:
    !
    call obs_set_current_header_list(obsSpaceData,'GL')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER
      idate = obs_headElem_i(obsSpaceData, OBS_DAT, headerIndex)
      trackCellNum = obs_headElem_i(obsSpaceData, OBS_FOV, headerIndex)
      call obs_set_current_body_list(obsSpaceData, headerIndex)
      BODY: do 
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if (bodyIndex < 0) exit BODY
        !
        !     *  Process only ASCAT backscatter anisotropy observations
        !
        varno = obs_bodyElem_i(obsSpaceData, OBS_VNM , bodyIndex)
        if (varno == BUFR_ICES) then
           write(ccyymmdd, FMT='(i8.8)') idate
           read(ccyymmdd(5:6), FMT='(i2)') imonth
           conc = col_getElem(columnTrlOnTrlLev,1,headerIndex,'GL')
           obsErrStdDev = SQRT(((1.0-conc)*ascatAnisSigmaOpenWater(trackCellNum,imonth))**2 + &
                                       (conc*ascatAnisSigmaIce(trackCellNum,imonth))**2)

           call obs_bodySet_r(obsSpaceData, OBS_OER, bodyIndex, obsErrStdDev)

        end if
      end do BODY
    end do HEADER

    if (.not. beSilent) write(*,*) myName//': done'

  end subroutine oer_setErrBackScatAnisIce

  !--------------------------------------------------------------------------
  ! chm_read_obs_err_stddev
  !--------------------------------------------------------------------------
  subroutine chm_read_obs_err_stddev
    !
    !:Purpose: To read observation errors from auxiliary file or observation
    !          file.
    !
    implicit none

    integer, parameter :: ndim=1
    integer :: istnid

    ! read the error std. dev. information from the auxiliary file
    call chm_read_obs_err_stddev_file

    ! set size of observation sub-space std array
    allocate(chm_std%obsStdDev(chm_std%n_stnid))

    ! read from the observation file if requested
    do istnid=1,chm_std%n_stnid
       if (chm_std%source(istnid).ge.1) then
          
          ! retrieve data from stats blocks (with bkstp=14 and block_type='DATA')
          chm_std%obsStdDev(istnid) = obsf_obsSub_read('CH',chm_std%stnids(istnid),chm_std%element(istnid), &
                                 chm_std%n_lvl(istnid),ndim,bkstp_opt=14,block_opt='DATA', &
                                 match_nlev_opt=chm_std%source(istnid).eq.1)

       end if
    end do

  end subroutine chm_read_obs_err_stddev

  !--------------------------------------------------------------------------
  ! chm_read_obs_err_stddev_file
  !--------------------------------------------------------------------------
  subroutine chm_read_obs_err_stddev_file
    !
    !:Purpose:  To read and store observation error std. dev. as needed for CH
    !           family obs.
    !
  implicit none

  integer :: FNOM, FCLOS
  integer :: IERR, JLEV, JELM, nulstat, ios, isize, icount
  logical :: LnewExists
  character(len=11) :: chemAuxObsDataFile = 'obsinfo_chm'
  
  character (len=128) :: ligne

  EXTERNAL FNOM,FCLOS

  ! Initialization

  chm_std%n_stnid=0
  
  ! CHECK THE EXISTENCE OF THE NEW FILE WITH STATISTICS

  INQUIRE(FILE=trim(chemAuxObsDataFile),EXIST=LnewExists)
  IF (.not.LnewExists) then
    WRITE(*,*) '---------------------------------------------------------------'
    WRITE(*,*) 'WARNING! chm_read_obs_err_stddev: auxiliary file ' // trim(chemAuxObsDataFile) 
    WRITE(*,*) 'WARNING! not available. Default CH family stddev to be applied if needed.'
    WRITE(*,*) '---------------------------------------------------------------'
    return
  ENDIF

  ! Read observation error std dev. from auxiliary file for constituent data

  NULSTAT=0
  IERR=FNOM(NULSTAT,trim(chemAuxObsDataFile),'SEQ',0)
  IF (IERR .EQ. 0) THEN
    open(unit=nulstat, file=trim(chemAuxObsDataFile), status='OLD')
  ELSE
    CALL utl_abort('chm_read_obs_err_stddev_file: COULD NOT OPEN AUXILIARY FILE ' // trim(chemAuxObsDataFile))
  ENDIF

  ! Read error standard deviations for constituents if available.
  ! (CH family; ozone and others)
  
  ios=0
  read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne
  do while (trim(adjustl(ligne(1:12))).ne.'SECTION I:') 
      read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne
  end do    
  
  ! Read number of observation set sub-families (STNIDs and ...) and allocate space

  read(nulstat,*,iostat=ios,err=10,end=10) chm_std%n_stnid
  read(nulstat,*,iostat=ios,err=10,end=10) isize
  
  allocate(chm_std%stnids(chm_std%n_stnid))
  allocate(chm_std%std_type(chm_std%n_stnid),chm_std%n_lat(chm_std%n_stnid))
  allocate(chm_std%source(chm_std%n_stnid),chm_std%ibegin(chm_std%n_stnid))
  allocate(chm_std%element(chm_std%n_stnid),chm_std%n_lvl(chm_std%n_stnid))
  allocate(chm_std%std1(isize),chm_std%std2(chm_std%n_stnid),chm_std%std3(chm_std%n_stnid))
  allocate(chm_std%levels(isize),chm_std%lat(isize))
 
  chm_std%element(:)=0
  chm_std%source(:)=0
  chm_std%std_type(:)=0
  chm_std%n_lvl(:)=1
  chm_std%n_lat(:)=1

  ! Begin reading for each sub-family
  ! Important: Combination of STNID, BUFR element and number of vertical levels
  !            to determine association to the observations.

  icount=0
  STNIDLOOP: do jelm=1,chm_std%n_stnid
    chm_std%ibegin(jelm)=icount+1

    ! disregard line of dashes
    read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne

    ! Read STNID (* as wildcard)    
    read(nulstat,'(2X,A9)',iostat=ios,err=10,end=10) chm_std%stnids(jelm) 
    
    !   Read (1) BUFR element,
    !        (2) Flag indication if EOR provided from this auxiliary file or
    !            to be read from the observation file,
    !        (3) Index specifying OER setup method,
    !        (4) Number of vertical levels
    !        (5) Number of latitudes
    !
    !   Important: Combination of STNID, BUFR element and number of vertical levels
    !              to determine association to the observations.
    
    read(nulstat,*,iostat=ios,err=10,end=10) chm_std%element(jelm),chm_std%source(jelm),  &
       chm_std%std_type(jelm), chm_std%n_lvl(jelm), chm_std%n_lat(jelm),  &
       chm_std%std2(jelm), chm_std%std3(jelm)

    if (chm_std%n_lvl(jelm).lt.1) chm_std%n_lvl(jelm)=1
    if (chm_std%n_lat(jelm).lt.1) chm_std%n_lat(jelm)=1
    
    if (icount+chm_std%n_lvl(jelm)*chm_std%n_lat(jelm).gt.isize) then
       write(*,'(10X,"Max array size exceeded: ",I6)') isize
       CALL utl_abort('chm_read_obs_err_stddev_file: PROBLEM READING OBSERR STD DEV.')    
    end if

    ! disregard line of dashes
    read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne
    
    ! disregard data section if not needed
    if (chm_std%std_type(jelm).eq.1.or.chm_std%std_type(jelm).eq.2.or.(chm_std%source(jelm).ge.1.and.chm_std%std_type(jelm).eq.0)) &
         cycle STNIDLOOP 

    if (chm_std%n_lvl(jelm).eq.1.and.chm_std%n_lat(jelm).eq.1) then
    
       ! Read one value only (independent of level and latitude)
       
       icount=icount+1
       read(nulstat,*,iostat=ios,err=10,end=10) chm_std%std1(icount)

    else if (chm_std%n_lvl(jelm).eq.1.and.chm_std%n_lat(jelm).gt.1) then
    
!      Value dependent on latitude only
       
       ! Read reference latitudes (must be in order of increasing size)
       
       read(nulstat,*,iostat=ios,err=10,end=10)                      &
              chm_std%lat(icount+1:icount+chm_std%n_lat(jelm))
      
       ! Read OER-related values
  
       read(nulstat,*,iostat=ios,err=10,end=10)                 &
                 chm_std%std1(icount+1:icount+chm_std%n_lat(jelm))

       icount=icount+chm_std%n_lat(jelm)

    else if (chm_std%n_lvl(jelm).gt.1.and.chm_std%n_lat(jelm).eq.1) then
    
       ! Value dependent on vertical level only
      
       do jlev=1,chm_std%n_lvl(jelm)
          icount=icount+1
          
          ! Read vertical level and OER-related value.
          
          read(nulstat,*,iostat=ios,err=10,end=10)                 &
                 chm_std%levels(icount),chm_std%std1(icount)

       end do
   
    else if (chm_std%n_lvl(jelm).gt.1.and.chm_std%n_lat(jelm).gt.1) then
    
       ! Value dependent on vertical level and latitude 
       
       ! Read reference latitudes (must be in order of increasing size)
       read(nulstat,*,iostat=ios,err=10,end=10)                      &
              chm_std%lat(icount+1:icount+chm_std%n_lat(jelm))
       ! write(*, '(10X,20F9.3)') chm_std%lat(icount+1:icount+chm_std%n_lat(jelm))
      
       do jlev=1,chm_std%n_lvl(jelm)
          
          ! Read vertical level and OER-related lat-dependent values.
          
          read(nulstat,*,iostat=ios,err=10,end=10)                   &
                 chm_std%levels(icount+jlev),                           &
                 chm_std%std1(icount+(jlev-1)*chm_std%n_lat(jelm)+1:icount+jlev*chm_std%n_lat(jelm))

       end do
       icount=icount+chm_std%n_lat(jelm)*chm_std%n_lvl(jelm)
    end if
 end do STNIDLOOP
   
 10 if (ios.gt.0) then
       WRITE(*,*) 'File read error message number: ',ios
       CALL utl_abort('chm_read_obs_err_stddev_file: PROBLEM READING OBSERR STD DEV.')    
    end if
 
 11 CLOSE(UNIT=NULSTAT)
    IERR=FCLOS(NULSTAT)    

  end subroutine chm_read_obs_err_stddev_file

  !--------------------------------------------------------------------------
  ! chm_obs_err_stddev_index
  !--------------------------------------------------------------------------
  subroutine chm_obs_err_stddev_index(CSTNID,NLEV,VARNO,ZLAT,ISTNID,JINT)
    ! 
    !:Purpose: To return the station ID and latitude indices corresponding to a
    !          measurement.
    !
    implicit none

    character(len=*), intent(in)  :: CSTNID
    integer, intent(in)           :: NLEV,VARNO
    real(8), intent(in)           :: ZLAT
    integer, intent(out)          :: ISTNID,JINT
    integer                       :: JN,ibegin
    real(8)                       :: lat

    !  Important: Combination of STNID, BUFR element and number of vertical levels
    !             to determine association to the observations.

    !             Find stnid with same number of vertical levels and same BUFR element.
    !             Note: * in chm_std%stnids stands for a wildcard
     
    ISTNID=0
    DO JN=1,chm_std%n_stnid

       ! First compare STNID values allowing for * and blanks in 
       ! chm_std%stnids(JN) as wildcards

       IF (utl_stnid_equal(chm_std%stnids(JN),CSTNID)) THEN
          IF ((NLEV.EQ.chm_std%n_lvl(JN) .OR. chm_std%source(JN).eq.2) .AND. VARNO.EQ.chm_std%element(JN)) THEN
             ISTNID=JN
             exit
          END IF
       END IF

    END DO

    IF (ISTNID.EQ.0) THEN
       write(*,*) 'chm_obs_err_stddev_index: Error std. dev. is unavailable for STNID ' // trim(CSTNID) // & 
                  ' and NLEV = ' // trim(utl_str(NLEV)) // ' and VARNO = ' // trim(utl_str(VARNO))
       write(*,*)
       write(*,*) ' Contents of chm_std (n_stnid = ' // trim(utl_str(chm_std%n_stnid)) // '):'
       if (chm_std%n_stnid.gt.0) then
          write(*,'(A)') '  STNID                 BUFR     NLEVELS'
          do jn=1,chm_std%n_stnid
             write(*,'(2X,A,2X,I12,2X,I10)') chm_std%stnids(jn),chm_std%element(jn),chm_std%n_lvl(jn)
          end do
       end if
       write(*,*)
       call utl_abort('chm_obs_err_stddev_index: Check section I of the auxiliary file.')
    ELSE

       IF (chm_std%n_lat(ISTNID) .GT. 1) THEN

          ! Find latitude index for interpolation.
          ! Assuming increasing latitudes in chm_std%lat

          lat = zlat / MPC_RADIANS_PER_DEGREE_R8  ! radians to degrees

          ibegin=chm_std%ibegin(ISTNID)-1
          IF (lat .GE. chm_std%lat(ibegin+chm_std%n_lat(ISTNID))) THEN
             JINT=chm_std%n_lat(ISTNID)+1
          ELSE
             DO JINT=1,chm_std%n_lat(ISTNID)
                IF (lat .LE. chm_std%lat(ibegin+JINT)) exit
             END DO
          END IF
                                           
       END IF       
    END IF         

  end subroutine chm_obs_err_stddev_index

  !--------------------------------------------------------------------------
  ! chm_get_obs_err_stddev
  !--------------------------------------------------------------------------
  function chm_get_obs_err_stddev(cstnid,nlev,varno,zlat,zlon,idate,itime,zval,&
                                  zlev,ilev,ifirst) result(obs_err_stddev) 
    ! 
    !:Purpose: To return the observational error std dev for a CH family
    !          measurement
    implicit none
    real(8)  :: obs_err_stddev 
   
    ! Arguments:
    character(len=*), intent(in) :: CSTNID ! station ID
    integer, intent(in) :: NLEV ! number of levels
    integer, intent(in) :: VARNO ! BUFR number
    real(8), intent(in) :: ZLAT ! latitude (radians)
    real(8), intent(in) :: ZLON ! longitude (radians)
    integer, intent(in) :: IDATE ! date in YYYYMMDD format
    integer, intent(in) :: ITIME ! time in HHMM format
    real(8), intent(in) :: ZVAL  ! observation values
    real(8), intent(in) :: ZLEV  ! vertical coordinate value
    integer, intent(in) :: ILEV  ! observation number in the profile
    logical, intent(in) :: IFIRST! true:  first call for a profile

    ! Locals:
    real(8) :: wgt,zwb,sigma
    integer :: ibegin,JLEV,JN,stat

    integer, save :: ISTNID,JINT
    
    ! If this call is for the first level for this measurement, get
    ! the station ID and latitude indices corresponding to this measurement
    if (ifirst) call chm_obs_err_stddev_index(CSTNID,NLEV,VARNO,ZLAT,ISTNID,JINT)                  
            
    ! Get weighting of error std. dev. if required

    if (chm_std%std_type(ISTNID).gt.2 .or. &
       (chm_std%source(ISTNID).eq.0 .and. chm_std%std_type(ISTNID).eq.0)) then

       IF (chm_std%n_lvl(ISTNID) .GT. 1) THEN
                 
          ! Find nearest vertical level (no interpolation)
                 
          zwb=1.E10
          ibegin=chm_std%ibegin(ISTNID)-1
          DO JN=1,chm_std%n_lvl(ISTNID)
             IF (zwb .GT. abs(ZLEV-chm_std%levels(ibegin+JN))) THEN
                JLEV=JN
                zwb=abs(ZLEV-chm_std%levels(ibegin+JN))
             END IF
          END DO
          JLEV=ibegin+(JLEV-1)*chm_std%n_lat(ISTNID)+1
       ELSE
          JLEV=chm_std%ibegin(ISTNID)
       END IF

       IF (chm_std%n_lat(ISTNID) .GT. 1) THEN
                
          ! Apply interpolation

          JLEV=JLEV+JINT-1
          ibegin=chm_std%ibegin(ISTNID)-1
          IF (JINT.EQ.1.OR.JINT.GT.chm_std%n_lat(ISTNID)) THEN
             wgt=chm_std%std1(JLEV)
          ELSE
             wgt=(chm_std%std1(JLEV-1)*(chm_std%lat(ibegin+JINT)-ZLAT)+ &
                  chm_std%std1(JLEV)*(ZLAT-chm_std%lat(ibegin+JINT-1)))/ &
                  (chm_std%lat(ibegin+JINT)-chm_std%lat(ibegin+JINT-1))
          END IF
       ELSE
          wgt=chm_std%std1(JLEV)             
       END IF
         
    end if
             
    ! Set the error std. dev.
                   
    IF (chm_std%source(ISTNID).EQ.0) THEN
               
       ! Set error standard deviations from scratch using content of
       ! previously read content of the auxiliary file.
                
       select case(chm_std%std_type(ISTNID))
       case(0)
          obs_err_stddev = wgt
       case(1)
          obs_err_stddev = max(chm_std%std3(ISTNID),chm_std%std2(ISTNID)*ZVAL)
       case(2)
          obs_err_stddev = sqrt(chm_std%std3(ISTNID)**2+(chm_std%std2(ISTNID)*ZVAL)**2)
       case(3)
          obs_err_stddev = min(chm_std%std3(ISTNID),max(chm_std%std2(ISTNID),wgt*ZVAL))
       case(4)
          obs_err_stddev = sqrt(chm_std%std2(ISTNID)**2+(wgt*ZVAL)**2)
       case default
          call utl_abort('chm_get_obs_err_stddev: std_type = ' // trim(utl_str(chm_std%std_type(ISTNID))) // &
               ' for STNID = ' // trim(CSTNID) // ' is not recognized.')
       end select

    ELSE

       ! Adjust error standard deviations read from observation file if requested.

       sigma = oss_obsdata_get_element(chm_std%obsStdDev(istnid), oss_obsdata_get_header_code(zlon,zlat,idate,itime,cstnid), ilev, stat_opt=stat)

       select case(stat)
       case(1)
          call utl_abort("chm_get_obs_err_stddev: No reports available for STNID = " // trim(cstnid) // &
                       ", nlev = " // trim(utl_str(nlev)) // ", varno = " // trim(utl_str(varno)))
       case(2)
          call utl_abort("chm_get_obs_err_stddev: Report not found for STNID = " // trim(cstnid) // &
                       ", nlev = " // trim(utl_str(nlev)) // ", varno = " // trim(utl_str(varno)))
       end select

       select case(chm_std%std_type(ISTNID))
       case(0)
          obs_err_stddev = sigma
       case(1)
          obs_err_stddev = max(chm_std%std3(ISTNID),chm_std%std2(ISTNID)*sigma)
       case(2)
          obs_err_stddev = sqrt(chm_std%std3(ISTNID)**2+(chm_std%std2(ISTNID)*sigma)**2)
       case(3)
          obs_err_stddev = min(chm_std%std3(ISTNID),max(chm_std%std2(ISTNID),wgt*sigma))
       case(4)
          obs_err_stddev = sqrt(chm_std%std2(ISTNID)**2+(wgt*sigma)**2)
       case default
          call utl_abort('chm_get_obs_err_stddev: std_type = ' // trim(utl_str(chm_std%std_type(ISTNID))) // &
               ' for STNID = ' // trim(CSTNID) // ' is not recognized.')
       end select
       
    END IF
    
  end function chm_get_obs_err_stddev

  !--------------------------------------------------------------------------
  ! chm_dealloc_obs_err_stddev
  !--------------------------------------------------------------------------
  subroutine chm_dealloc_obs_err_stddev
    ! 
    !:Purpose: To deallocate temporary storage space used for observation errors
    !          for the CH family.
    !
    implicit none

    integer :: istnid

    if (chm_std%n_stnid.eq.0) return
    
    if (allocated(chm_std%obsStdDev)) then
       do istnid=1,chm_std%n_stnid
          if (chm_std%source(istnid).ge.1) call oss_obsdata_dealloc(chm_std%obsStdDev(istnid))
       end do
       deallocate(chm_std%obsStdDev)       
    end if

    if (allocated(chm_std%stnids))   deallocate(chm_std%stnids)
    if (allocated(chm_std%n_lvl))    deallocate(chm_std%n_lvl)
    if (allocated(chm_std%std_type)) deallocate(chm_std%std_type)
    if (allocated(chm_std%ibegin))   deallocate(chm_std%ibegin)
    if (allocated(chm_std%element))  deallocate(chm_std%element)
    if (allocated(chm_std%source))   deallocate(chm_std%source)
    if (allocated(chm_std%n_lat))    deallocate(chm_std%n_lat)
    if (allocated(chm_std%std1))     deallocate(chm_std%std1)
    if (allocated(chm_std%std2))     deallocate(chm_std%std2)
    if (allocated(chm_std%std3))     deallocate(chm_std%std3)
    if (allocated(chm_std%levels))   deallocate(chm_std%levels)
    if (allocated(chm_std%lat))      deallocate(chm_std%lat)

  end subroutine chm_dealloc_obs_err_stddev

  !--------------------------------------------------------------------------
  ! oer_getSSTdataParam_char
  !--------------------------------------------------------------------------
  function oer_getSSTdataParam_char(item, itemIndex) result(value)
    !
    !:Purpose: get character item value from SSTdataParams derived type
    !
    implicit none
    
    character(len=20) :: value 

    ! Arguments:
    character(len=*), intent(in) :: item
    integer         , intent(in) :: itemIndex

    select case(trim(item))      
      case('dataType')
        value = SSTdataParams(itemIndex)%dataType
      case('instrument')
        value = SSTdataParams(itemIndex)%instrument
      case('sensor')
        value = SSTdataParams(itemIndex)%sensor
      case('sensorType')
        value = SSTdataParams(itemIndex)%sensorType
      case default
        call utl_abort('oer_getSSTdataParam_char: invalid item '//(trim(item)))
    end select

  end function oer_getSSTdataParam_char

  !--------------------------------------------------------------------------
  ! oer_getSSTdataParam_int
  !--------------------------------------------------------------------------
  function oer_getSSTdataParam_int(item, itemIndex_opt) result(value)
    !
    !:Purpose: get integer item value from SSTdataParams derived type
    !
    implicit none
    
    integer :: value 

    ! Arguments:
    character(len=*), intent(in)           :: item
    integer         , intent(in), optional :: itemIndex_opt

    if (present(itemIndex_opt)) then
      select case(trim(item))      
        case('codeType')
          value = SSTdataParams(itemIndex_opt)%codeType
        case default
          call utl_abort('oer_getSSTdataParam_int: invalid item '//(trim(item)))
      end select
    else
      select case(trim(item))
        case('maxNumberSSTDatasets')
          value = maxNumberSSTDatasets
        case('numberSSTDatasets')
          value = numberSSTDatasets
        case default
          call utl_abort('oer_getSSTdataParam_int: invalid item '//(trim(item)))
      end select

    end if

  end function oer_getSSTdataParam_int

  !--------------------------------------------------------------------------
  ! oer_getSSTdataParam_R8
  !--------------------------------------------------------------------------
  function oer_getSSTdataParam_R8(item, itemIndex) result(value)
    !
    !:Purpose: get real(8) item value from SSTdataParams derived type
    !
    implicit none
    
    real(8) :: value 

    ! Arguments:
    character(len=*), intent(in) :: item
    integer         , intent(in) :: itemIndex

    select case(trim(item))      
      case('dayError')
        value = SSTdataParams(itemIndex)%dayError
      case('nightError')
        value = SSTdataParams(itemIndex)%nightError
      case default
        call utl_abort('oer_getSSTdataParam_R8: invalid item '//(trim(item)))
    end select

  end function oer_getSSTdataParam_R8

end module obsErrors_mod
