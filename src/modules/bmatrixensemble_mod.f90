!--------------------------------------- LICENCE BEGIN -----------------------------------
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

module BmatrixEnsemble_mod
  ! MODULE BmatrixEnsemble_mod (prefix='ben' category='2. B and R matrices')
  !
  ! :Purpose: Performs transformation from control vector to analysis increment 
  !           using the spatially localized ensemble covariance matrix. This 
  !           module works for both global and limited-area applications.
  !
  use ramDisk_mod
  use midasMpi_mod
  use fileNames_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use ensembleStateVector_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use localization_mod
  use mathPhysConstants_mod
  use earthConstants_mod
  use gridVariableTransforms_mod
  use utilities_mod
  use globalSpectralTransform_mod
  use lamSpectralTransform_mod
  use spectralFilter_mod
  use varNameList_mod
  use advection_mod
  use gridBinning_mod
  use humidityLimits_mod
  use analysisGrid_mod
  implicit none
  save
  private

  ! public procedures
  public :: ben_Setup, ben_BSqrt, ben_BSqrtAd, ben_writeAmplitude
  public :: ben_reduceToMPILocal, ben_reduceToMPILocal_r4, ben_expandToMPIGlobal, ben_expandToMPIGlobal_r4
  public :: ben_getScaleFactor, ben_getnEns, ben_getPerturbation, ben_getEnsMean, ben_Finalize
  public :: ben_setFsoLeadTime, ben_getNumStepAmplitudeAssimWindow, ben_getAmplitudeAssimWindow
  public :: ben_getAmp3dStepIndexAssimWindow, ben_getNumLoc, ben_getLoc, ben_getNumInstance

  integer, parameter   :: maxNumLevels      = 200
  integer, parameter   :: maxNumLocalLength =  20
  integer, parameter   :: nInstanceMax      =  10

  type :: struct_bEns
    logical             :: initialized = .false.

    real(8),allocatable :: scaleFactor_M(:), scaleFactor_T(:), scaleFactor_DP(:)
    real(8)             :: scaleFactor_SF

    integer             :: ni, nj, myLonBeg, myLonEnd, myLatBeg, myLatEnd
    integer             :: nLevInc_M, nLevInc_T, nLevInc_DP, nLevEns_M, nLevEns_T, nLevEns_DP
    integer             :: topLevIndex_M, topLevIndex_T, topLevIndex_DP
    integer             :: nEnsOverDimension
    integer             :: cvDim_mpilocal, cvDim_mpiglobal
    integer             :: numStep, numStepAssimWindow
    integer             :: numStepAmplitudeFSOFcst, numStepAmplitudeAssimWindow, numStepAdvectAssimWindow
    integer             :: numSubEns
    integer,allocatable :: dateStampList(:)
    integer,allocatable :: dateStampListAdvectedFields(:)
    
    integer             :: numIncludeAnlVar
    
    ! FSO
    real(8)             :: fsoLeadTime = -1.0D0
    integer             :: numStepAdvectFSOFcst
    
    ! Localizations
    integer             :: nWaveBand
    integer             :: nWaveBandForFiltering = 0
    
    ! Ensemble perturbations
    type(struct_ens), allocatable :: ensPerts(:)
    
    ! Ensemble amplitude (only used in diagnostic mode)
    type(struct_ens)    :: ensAmplitudeStorage
    character(len=4)    :: varNameALFA(1)
    
    ! Localization
    type(struct_loc), pointer :: locStorage(:)
    
    ! The HU LQ mess
    logical :: gsvHUcontainsLQ
    logical :: ensShouldNotContainLQvarName 
    
    ! Vertical grid
    type(struct_vco), pointer :: vco_anl, vco_ens, vco_file => null()

    ! Horizontal grid
    type(struct_hco), pointer :: hco_ens  ! Ensemble   horizontal grid parameters
    type(struct_hco), pointer :: hco_core ! Core grid for limited area EnVar
    type(struct_hco), pointer :: hco_file ! Input file horizontal grid parameters
    
    ! Amplitude parameters
    type(struct_adv)          :: adv_amplitudeFSOFcst
    type(struct_adv), pointer :: adv_amplitudeAssimWindow
    type(struct_adv)          :: adv_ensPerts
    type(struct_adv)          :: adv_analInc

    integer           :: amp3dStepIndexAssimWindow
    integer           :: amp3dStepIndexFSOFcst
    
    ! Variance smoothing
    logical           :: ensPertsNormalized
    
    type(struct_gsv)  :: statevector_ensStdDev
    
    ! Optimization
    logical             :: useSaveAmp

    ! Namelist variables
    integer             :: nEns ! number of ensemble members
    real(8)             :: scaleFactor(maxNumLevels)
    real(8)             :: scaleFactorHumidity(maxNumLevels)
    real(8)             :: advectFactorFSOFcst(maxNumLevels)
    real(8)             :: advectFactorAssimWindow(maxNumLevels)
    integer             :: ntrunc
    character(len=256)  :: enspathname
    real(8)             :: hLocalize(maxNumLocalLength)
    real(8)             :: vLocalize(maxNumLocalLength)
    character(len=256)  :: LocalizationType
    integer             :: waveBandPeaks(maxNumLocalLength)
    integer             :: waveBandIndexSelected
    logical             :: ensDiagnostic
    logical             :: advDiagnostic
    character(len=2)    :: ctrlVarHumidity
    character(len=32)   :: advectTypeAssimWindow
    character(len=32)   :: advectStartTimeIndexAssimWindow
    logical             :: advectAmplitudeFSOFcst
    logical             :: advectAmplitudeAssimWindow = .false.
    logical             :: advectEnsPertAnlInc        = .false.
    logical             :: removeSubEnsMeans
    logical             :: keepAmplitude
    character(len=4)    :: IncludeAnlVar(vnl_numvarmax)
    logical             :: ensContainsFullField
    character(len=24)   :: varianceSmoothing
    real(8)             :: footprintRadius
    real(8)             :: footprintTopoThreshold
    logical             :: useCmatrixOnly
    integer             :: ensDateOfValidity
    character(len=20)   :: transformVarKindCH
    real(8)             :: huMinValue
    character(len=12)   :: hInterpolationDegree
  end type struct_bEns

  integer :: nInstance = 0 ! The number of Bens instances

  type(struct_bEns) :: bEns(nInstanceMax)

  type(struct_ens), target :: ensAmplitudeSave(nInstanceMax) ! Save this to allow early allocation 
                                                               ! for better efficiency

  character(len=15) :: ben_mode

  character(len=4), parameter  :: varNameALFAatmMM(1) = (/ 'ALFA' /)
  character(len=4), parameter  :: varNameALFAatmTH(1) = (/ 'ALFT' /)
  character(len=4), parameter  :: varNameALFAsfc(1)   = (/ 'ALFS' /)
  character(len=4), parameter  :: varNameALFAocean(1) = (/ 'ALFO' /)

  logical, parameter :: verbose = .false. ! Control parameter for the level of listing output

  integer, external    :: get_max_rss, omp_get_thread_num

  character(len=20)   :: transformVarKindCH                                   

CONTAINS

  !--------------------------------------------------------------------------
  ! ben_setup
  !--------------------------------------------------------------------------
  subroutine ben_setup(hco_anl_in, hco_core_in, vco_anl_in, cvDimPerInstance, &
                       mode_opt)
    !
    !:Purpose: To configure the ensemble B matrix
    !
    implicit none

    ! Arguments:
    type(struct_hco), pointer, intent(in) :: hco_anl_in
    type(struct_hco), pointer, intent(in) :: hco_core_in
    type(struct_vco), pointer, intent(in) :: vco_anl_in

    character(len=*), intent(in), optional :: mode_opt

    integer, allocatable, intent(out)      :: cvDimPerInstance(:)

    ! Locals:
    integer        :: fnom, fclos, ierr
    integer        :: cvDimStorage(nInstanceMax)
    integer        :: nulnam = 0

    ! Namelist variables
    integer             :: nEns
    real(8)             :: scaleFactor(maxNumLevels)
    real(8)             :: scaleFactorHumidity(maxNumLevels)
    real(8)             :: advectFactorFSOFcst(maxNumLevels)
    real(8)             :: advectFactorAssimWindow(maxNumLevels)
    integer             :: ntrunc
    character(len=256)  :: enspathname
    real(8)             :: hLocalize(maxNumLocalLength)
    real(8)             :: vLocalize(maxNumLocalLength)
    character(len=256)  :: localizationType
    integer             :: waveBandPeaks(maxNumLocalLength)
    integer             :: waveBandIndexSelected
    logical             :: ensDiagnostic
    logical             :: advDiagnostic
    character(len=2)    :: ctrlVarHumidity
    character(len=32)   :: advectTypeAssimWindow
    character(len=32)   :: advectStartTimeIndexAssimWindow
    logical             :: advectAmplitudeFSOFcst     = .false.
    logical             :: advectAmplitudeAssimWindow = .false.
    logical             :: advectEnsPertAnlInc        = .false.
    logical             :: removeSubEnsMeans
    logical             :: keepAmplitude
    character(len=4)    :: includeAnlVar(vnl_numvarmax)
    logical             :: ensContainsFullField
    character(len=24)   :: varianceSmoothing
    real(8)             :: footprintRadius
    real(8)             :: footprintTopoThreshold
    logical             :: useCmatrixOnly
    integer             :: ensDateOfValidity
    real(8)             :: huMinValue
    character(len=12)   :: hInterpolationDegree ! select degree of horizontal interpolation (if needed)

    ! Namelist
    NAMELIST /NAMBEN/nEns, scaleFactor, scaleFactorHumidity, ntrunc, enspathname,             &
         hLocalize, vLocalize, LocalizationType, waveBandPeaks, ensDiagnostic, advDiagnostic, &
         ctrlVarHumidity, advectFactorFSOFcst, advectFactorAssimWindow, removeSubEnsMeans,    &
         keepAmplitude, advectTypeAssimWindow, advectStartTimeIndexAssimWindow, IncludeAnlVar,&
         ensContainsFullField, varianceSmoothing, footprintRadius, footprintTopoThreshold,    &
         useCmatrixOnly, waveBandIndexSelected, ensDateOfValidity, transformVarKindCH,        &
         huMinValue, hInterpolationDegree

    if (verbose) write(*,*) 'Entering ben_Setup'

    call utl_tmg_start(54,'----B_ENS_Setup')

    !- Set the module mode
    if ( present(mode_opt) ) then
      if ( trim(mode_opt) == 'Analysis' .or. trim(mode_opt) == 'BackgroundCheck') then
        ben_mode = trim(mode_opt)
        if (mmpi_myid == 0) write(*,*)
        if (mmpi_myid == 0) write(*,*) 'ben_setup: Mode activated = ', trim(ben_mode)
      else
        write(*,*)
        write(*,*) 'mode = ', trim(mode_opt)
        call utl_abort('ben_setup: unknown mode')
      end if
    else
      ben_mode = 'Analysis'
      if (mmpi_myid == 0) write(*,*)
      if (mmpi_myid == 0) write(*,*) 'ben_setup: Analysis mode activated (by default)'
    end if

    !- Open the namelist and loop through it
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
 
    instanceLoop: do
 
      !- Set the default values for the namalist parameters
      scaleFactor(:)        =    0.0d0
      scaleFactorHumidity(:)=    1.0d0
      nEns                  =   10
      ntrunc                =   30
      enspathname           = 'ensemble'
      localizationType      = 'LevelDependent'
      waveBandPeaks(:)      =   -1.0d0
      waveBandIndexSelected =   -1
      ensDiagnostic         = .false.
      advDiagnostic         = .false.
      hLocalize(:)          =   -1.0d0
      hLocalize(1)          = 2800.0d0
      vLocalize(:)          =   -1.0d0
      vLocalize(1)          =    2.0d0
      ctrlVarHumidity       = 'LQ'
      advectFactorFSOFcst(:)=   0.0D0
      advectTypeAssimWindow = 'amplitude'
      advectStartTimeIndexAssimWindow = 'first'
      advectFactorAssimWindow(:) = 0.0D0
      removeSubEnsMeans     = .false.
      keepAmplitude         = .false.
      ensContainsFullField  = .true.
      includeAnlVar(:)      = ''
      varianceSmoothing     = 'none'
      footprintRadius        =  250.0d3 ! 250km
      footprintTopoThreshold =  200.0d0 ! 200 m
      useCmatrixOnly        = .false.
      ensDateOfValidity     = MPC_missingValue_INT ! i.e. undefined
      transformVarKindCH    = ''
      huMinValue            = MPC_missingValue_R8  ! i.e. undefined
      hInterpolationDegree  = 'LINEAR' ! or 'CUBIC' or 'NEAREST'
      
      !- Read the namelist
      read(nulnam,nml=namben,iostat=ierr)
      if (ierr /= 0) then
        if (nInstance >= 1) then
          exit instanceLoop
        else
          call utl_abort('ben_setup: Error reading the first instance namelist')
        end if
      end if

      if (mmpi_myid == 0) write(*,nml=namben)
 
      ! We have found a valid instance
      nInstance = nInstance + 1

      !- Adjust some namelist-dependent variables
      
      ! If zero weight, skip rest of setup
      if ( sum(scaleFactor(:)) == 0.0d0 ) then
        if (mmpi_myid == 0) write(*,*) 'ben_setup: scaleFactor=0, skipping rest of setup and exit instance loop'
        cvDimStorage(nInstance) = 0
        bEns(nInstance)%initialized = .false.
        exit instanceLoop
      end if

      if (nInstance > nInstanceMax) then
        call utl_abort('ben_setup: the number of instance exceed the maximum currently allowed')
      end if

      if (trim(ben_mode) == 'BackgroundCheck' .and. nInstance > 1) then
        call utl_abort('ben_setup: the background check mode is not compatible with multiple instance')
      end if

      if ( (huMinValue == MPC_missingValue_R8) .and. &
           gsv_varExist(varName='HU') .and. &
           (ctrlVarHumidity == 'LQ') ) then
        call utl_abort('ben_setup: the value of huMinValue must be specified in namelist NAMBEN')
      end if

      !- Transfer the info to the structure
      bEns(nInstance)%nEns                       = nEns
      bEns(nInstance)%scaleFactor(:)             = scaleFactor(:)
      bEns(nInstance)%scaleFactorHumidity(:)     = scaleFactorHumidity(:)
      bEns(nInstance)%nTrunc                     = nTrunc
      bEns(nInstance)%ensPathName                = ensPathName
      bEns(nInstance)%hLocalize(:)               = hLocalize(:)
      bEns(nInstance)%vLocalize(:)               = vLocalize(:)
      bEns(nInstance)%localizationType           = localizationType
      bEns(nInstance)%waveBandPeaks(:)           = waveBandPeaks(:)
      bEns(nInstance)%waveBandIndexSelected      = waveBandIndexSelected
      bEns(nInstance)%ensDiagnostic              = ensDiagnostic
      bEns(nInstance)%advDiagnostic              = advDiagnostic
      bEns(nInstance)%ctrlVarHumidity            = ctrlVarHumidity
      bEns(nInstance)%advectTypeAssimWindow      = advectTypeAssimWindow
      bEns(nInstance)%advectStartTimeIndexAssimWindow = advectStartTimeIndexAssimWindow
      bEns(nInstance)%advectFactorFSOFcst(:)     = advectFactorFSOFcst(:)
      bEns(nInstance)%advectFactorAssimWindow(:) = advectFactorAssimWindow(:)
      bEns(nInstance)%advectAmplitudeFSOFcst     = advectAmplitudeFSOFcst
      bEns(nInstance)%advectAmplitudeAssimWindow = advectAmplitudeAssimWindow
      bEns(nInstance)%advectEnsPertAnlInc        = advectEnsPertAnlInc
      bEns(nInstance)%removeSubEnsMeans          = removeSubEnsMeans
      bEns(nInstance)%keepAmplitude              = keepAmplitude
      bEns(nInstance)%includeAnlVar(:)           = includeAnlVar(:)
      bEns(nInstance)%ensContainsFullField       = ensContainsFullField
      bEns(nInstance)%varianceSmoothing          = varianceSmoothing
      bEns(nInstance)%footprintRadius            = footprintRadius
      bEns(nInstance)%footprintTopoThreshold     = footprintTopoThreshold
      bEns(nInstance)%useCmatrixOnly             = useCmatrixOnly
      bEns(nInstance)%ensDateOfValidity          = ensDateOfValidity
      bEns(nInstance)%transformVarKindCH         = transformVarKindCH
      bEns(nInstance)%huMinValue                 = huMinValue
      bEns(nInstance)%hInterpolationDegree       = hInterpolationDegree
      
      bEns(nInstance)%hco_ens  => hco_anl_in
      bEns(nInstance)%hco_core => hco_core_in
      bEns(nInstance)%vco_anl  => vco_anl_in

      !- Setup the LAM analysis grid metrics
      call agd_SetupFromHCO( hco_anl_in, hco_core_in ) ! IN
      
      !- Set the instance
      call ben_setupOneInstance(nInstance,               & ! IN
                                cvDimStorage(nInstance))   ! OUT

      bEns(nInstance)%initialized = .true.
      
    end do instanceLoop

    !- Close the namelist
    ierr = fclos(nulnam)

    !- Set the output control variable dimensions array
    allocate(cvDimPerInstance(nInstance))
    cvDimPerInstance(:) = cvDimStorage(1:nInstance)

    call utl_tmg_stop(54)

  end subroutine ben_setup

  !--------------------------------------------------------------------------
  ! ben_setupOneInstance
  !--------------------------------------------------------------------------
  subroutine ben_setupOneInstance(instanceIndex, cvDim)
    !
    !:Purpose: To configure a single instance of the ensemble B matrix
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: instanceIndex
    integer, intent(out) :: cvDim

    ! Locals:
    type(struct_gsv) :: statevector_ensMean4D, statevector_oneEnsPert4D

    type(struct_gbi) :: gbi_horizontalMean, gbi_landSeaTopo

    real(8) :: pSurfRef, delT_hour

    real(8), allocatable :: advectFactorFSOFcst_M(:),advectFactorAssimWindow_M(:)

    real(8),pointer :: vertLocationEns(:), vertLocationFile(:), vertLocationInc(:)

    real(4), pointer :: bin2d(:,:,:)
    real(8), pointer :: HeightSfc(:,:)

    integer        :: lonPerPE, latPerPE, lonPerPEmax, latPerPEmax
    integer        :: myMemberBeg, myMemberEnd, myMemberCount, maxMyMemberCount
    integer        :: levIndex, jvar, status, ierr
    integer        :: waveBandIndex, stepIndex
    character(len=256) :: ensFileName
    integer        :: dateStampFSO, ensDateStampOfValidity, idate, itime, newdate

    logical        :: EnsTopMatchesAnlTop, useAnlLevelsOnly

    character(len=32)   :: direction, directionEnsPerts, directionAnlInc

    if (verbose) write(*,*) 'Entering ben_SetupOneInstance'
    
    write(*,*) 'ben_setupOneInstance: enspathname = ', trim(bEns(instanceIndex)%ensPathName)

    !
    !- 1.  B matrix configuration
    !

    !- 1.1 Number of time step bins
    bEns(instanceIndex)%numStep = tim_nstepobsinc
    if (bEns(instanceIndex)%numStep /= 1.and.bEns(instanceIndex)%numStep /= 3.and.bEns(instanceIndex)%numStep /= 5.and.bEns(instanceIndex)%numStep /= 7) then
      call utl_abort('ben_setupOneInstance: Invalid value for numStep (choose 1 or 3 or 5 or 7)!')
    end if

    !- 1.2 FSO-related options
    bEns(instanceIndex)%numStepAssimWindow = bEns(instanceIndex)%numStep
    if (bEns(instanceIndex)%fsoLeadTime > 0.0D0) then
      bEns(instanceIndex)%numStep = bEns(instanceIndex)%numStep + 1
      call incdatr(dateStampFSO, tim_getDatestamp(), bEns(instanceIndex)%fsoLeadTime)
    end if

    allocate(bEns(instanceIndex)%dateStampList(bEns(instanceIndex)%numStep))
    if (bEns(instanceIndex)%fsoLeadTime > 0.0D0) then
      call tim_getstamplist(bEns(instanceIndex)%dateStampList,bEns(instanceIndex)%numStep-1,tim_getDatestamp())
      bEns(instanceIndex)%dateStampList(bEns(instanceIndex)%numStep) = dateStampFSO
    else
      if (bEns(instanceIndex)%ensDateOfValidity == MPC_missingValue_INT) then
        call tim_getstamplist(bEns(instanceIndex)%dateStampList,bEns(instanceIndex)%numStep,tim_getDatestamp())
      else
        if (bEns(instanceIndex)%numStep == 1) then
          if (bEns(instanceIndex)%ensDateOfValidity == -1) then
            ensDateStampOfValidity = bEns(instanceIndex)%ensDateOfValidity
          else
            idate = bEns(instanceIndex)%ensDateOfValidity/100
            itime = (bEns(instanceIndex)%ensDateOfValidity-idate*100)*1000000
            ierr = newdate(ensDateStampOfValidity, idate, itime, 3)
          end if
          bEns(instanceIndex)%dateStampList(:) = ensDateStampOfValidity
        else
          call utl_abort('ben_setupOneInstance: A single date of validity cannot be specified for numStep > 1')
        end if
      end if
    end if

    !- 1.3 Horizontal grid
    bEns(instanceIndex)%ni = bEns(instanceIndex)%hco_ens%ni
    bEns(instanceIndex)%nj = bEns(instanceIndex)%hco_ens%nj
    if (bEns(instanceIndex)%hco_ens%global) then
      if (mmpi_myid == 0) write(*,*)
      if (mmpi_myid == 0) write(*,*) 'ben_setupOneInstance: GLOBAL mode activated'
    else
      if (mmpi_myid == 0) write(*,*)
      if (mmpi_myid == 0) write(*,*) 'ben_setupOneInstance: LAM mode activated'
    end if

    !- 1.4 Vertical levels
    if ( mmpi_myid == 0 ) then
      call fln_ensfileName(ensFileName, bEns(instanceIndex)%ensPathName, memberIndex_opt=1)
      call vco_SetupFromFile(bEns(instanceIndex)%vco_file, ensFileName)
    end if
    call vco_mpiBcast(bEns(instanceIndex)%vco_file)

    !- Do we need to read all the vertical levels from the ensemble?
    useAnlLevelsOnly = vco_subsetOrNot(bEns(instanceIndex)%vco_anl, bEns(instanceIndex)%vco_file)
    if ( useAnlLevelsOnly ) then
      write(*,*)
      write(*,*) 'ben_setupOneInstance: only the analysis levels will be read in the ensemble '
      bEns(instanceIndex)%vco_ens  => bEns(instanceIndex)%vco_anl ! the ensemble target grid is the analysis grid
      call vco_deallocate(bEns(instanceIndex)%vco_file)
      bEns(instanceIndex)%vco_file => bEns(instanceIndex)%vco_anl ! only the analysis levels will be read in the ensemble
      EnsTopMatchesAnlTop = .true.
    else
      write(*,*)
      write(*,*) 'ben_setupOneInstance: all the vertical levels will be read in the ensemble '
      if ( bEns(instanceIndex)%vco_anl%nLev_M > 0 .and. bEns(instanceIndex)%vco_anl%vgridPresent ) then
        pSurfRef = 101000.D0
        nullify(vertLocationInc)
        status = vgd_levels( bEns(instanceIndex)%vco_anl%vgrid, &
                             ip1_list=bEns(instanceIndex)%vco_anl%ip1_M, &
                             levels=vertLocationInc, &
                             sfc_field=pSurfRef, in_log=.false. )
        if (status /= VGD_OK) call utl_abort('ben_setupOneInstance: ERROR from vgd_levels')
        nullify(vertLocationFile)
        status = vgd_levels( bEns(instanceIndex)%vco_file%vgrid, &
                             ip1_list=bEns(instanceIndex)%vco_file%ip1_M, &
                             levels=vertLocationFile, &
                             sfc_field=pSurfRef, in_log=.false.)
        if (status /= VGD_OK) call utl_abort('ben_setupOneInstance: ERROR from vgd_levels')
      
        do levIndex = 1, bEns(instanceIndex)%vco_anl%nLev_M
          vertLocationInc(levIndex) = log(vertLocationInc(levIndex))
        end do
        do levIndex = 1, bEns(instanceIndex)%vco_file%nLev_M
          vertLocationFile(levIndex) = log(vertLocationFile(levIndex))
        end do

        EnsTopMatchesAnlTop = abs( vertLocationFile(1) - vertLocationInc(1) ) < 0.1d0
        write(*,*) 'ben_setupOneInstance: EnsTopMatchesAnlTop, presEns, presInc = ', &
             EnsTopMatchesAnlTop, vertLocationFile(1), vertLocationInc(1)
        deallocate(vertLocationFile)
        deallocate(vertLocationInc)
      else
        ! not sure what this mean when no MM levels
        write(*,*) 'ben_setupOneInstance: nLev_M       = ', bEns(instanceIndex)%vco_anl%nLev_M
        write(*,*) 'ben_setupOneInstance: vgridPresent = ', bEns(instanceIndex)%vco_anl%vgridPresent
        EnsTopMatchesAnlTop = .true.
      end if

      if ( EnsTopMatchesAnlTop ) then
        if ( mmpi_myid == 0 ) write(*,*) 'ben_setupOneInstance: top level of ensemble member and analysis grid match'
        bEns(instanceIndex)%vco_ens => bEns(instanceIndex)%vco_anl  ! IMPORTANT: top levels DO match, therefore safe
        ! to force members to be on analysis vertical levels
      else
        if ( mmpi_myid == 0 ) write(*,*) 'ben_setupOneInstance: top level of ensemble member and analysis grid are different, therefore'
        if ( mmpi_myid == 0 ) write(*,*) '                      assume member is already be on correct levels - NO CHECKING IS DONE'
        bEns(instanceIndex)%vco_ens => bEns(instanceIndex)%vco_file ! IMPORTANT: top levels do not match, therefore must
        ! assume file is already on correct vertical levels
      end if
    end if
    
    if (bEns(instanceIndex)%vco_anl%Vcode /= bEns(instanceIndex)%vco_ens%Vcode) then
      write(*,*) 'ben_setupOneInstance: vco_anl%Vcode = ', bEns(instanceIndex)%vco_anl%Vcode, ', vco_ens%Vcode = ', bEns(instanceIndex)%vco_ens%Vcode
      call utl_abort('ben_setupOneInstance: vertical levels of ensemble not compatible with analysis grid')
    end if
    bEns(instanceIndex)%nLevEns_M  = bEns(instanceIndex)%vco_ens%nLev_M
    bEns(instanceIndex)%nLevEns_T  = bEns(instanceIndex)%vco_ens%nLev_T
    bEns(instanceIndex)%nLevEns_DP = bEns(instanceIndex)%vco_ens%nLev_Depth
    bEns(instanceIndex)%nLevInc_M  = bEns(instanceIndex)%vco_anl%nLev_M
    bEns(instanceIndex)%nLevInc_T  = bEns(instanceIndex)%vco_anl%nLev_T
    bEns(instanceIndex)%nLevInc_DP = bEns(instanceIndex)%vco_anl%nLev_Depth
    bEns(instanceIndex)%topLevIndex_M  = bEns(instanceIndex)%nLevInc_M  - bEns(instanceIndex)%nLevEns_M+1
    bEns(instanceIndex)%topLevIndex_T  = bEns(instanceIndex)%nLevInc_T  - bEns(instanceIndex)%nLevEns_T+1
    bEns(instanceIndex)%topLevIndex_DP = bEns(instanceIndex)%nLevInc_DP - bEns(instanceIndex)%nLevEns_DP+1

    if (bEns(instanceIndex)%vco_anl%Vcode == 5002) then
      if ( (bEns(instanceIndex)%nLevEns_T /= (bEns(instanceIndex)%nLevEns_M+1)) .and. (bEns(instanceIndex)%nLevEns_T /= 1 .or. bEns(instanceIndex)%nLevEns_M /= 1) ) then
        write(*,*) 'ben_setupOneInstance: nLevEns_T, nLevEns_M = ',bEns(instanceIndex)%nLevEns_T,bEns(instanceIndex)%nLevEns_M
        call utl_abort('ben_setupOneInstance: Vcode=5002, nLevEns_T must equal nLevEns_M+1!')
      end if
    else if (bEns(instanceIndex)%vco_anl%Vcode == 5005) then
      if ( bEns(instanceIndex)%nLevEns_T /= bEns(instanceIndex)%nLevEns_M .and. &
           bEns(instanceIndex)%nLevEns_T /= 0 .and. &
           bEns(instanceIndex)%nLevEns_M /= 0 ) then
        write(*,*) 'ben_setup: nLevEns_T, nLevEns_M = ',bEns(instanceIndex)%nLevEns_T,bEns(instanceIndex)%nLevEns_M
        call utl_abort('ben_setupOneInstance: Vcode=5005, nLevEns_T must equal nLevEns_M!')
      end if
    else if (bEns(instanceIndex)%vco_anl%Vcode == 0) then
      if ( bEns(instanceIndex)%nLevEns_T /= 0 .and. bEns(instanceIndex)%nLevEns_M /= 0 ) then
        write(*,*) 'ben_setup: nLevEns_T, nLevEns_M = ',bEns(instanceIndex)%nLevEns_T, bEns(instanceIndex)%nLevEns_M
        call utl_abort('ben_setupOneInstance: surface-only case (Vcode=0), bEns(instanceIndex)%nLevEns_T and nLevEns_M must equal 0!')
      end if
    else
      write(*,*) 'vco_anl%Vcode = ',bEns(instanceIndex)%vco_anl%Vcode
      call utl_abort('ben_setupOneInstance: unknown vertical coordinate type!')
    end if

    if (bEns(instanceIndex)%nLevEns_M.gt.bEns(instanceIndex)%nLevInc_M) then
      call utl_abort('ben_setupOneInstance: ensemble has more levels than increment - not allowed!')
    end if

    if (bEns(instanceIndex)%nLevEns_M.lt.bEns(instanceIndex)%nLevInc_M) then
      if (mmpi_myid == 0) write(*,*) 'ben_setupOneInstance: ensemble has less levels than increment'
      if (mmpi_myid == 0) write(*,*) '                      some levels near top will have zero increment'
    end if

    !- 1.5 Bmatrix Weight
    if (bEns(instanceIndex)%vco_anl%Vcode == 5002 .or. bEns(instanceIndex)%vco_anl%Vcode == 5005) then
      if (bEns(instanceIndex)%nLevEns_M > 0) then
        ! Multi-level or momentum-level-only analysis
        bEns(instanceIndex)%varNameALFA(:) = varNameALFAatmMM(:)
      else
        ! Thermo-level-only analysis
        bEns(instanceIndex)%varNameALFA(:) = varNameALFAatmTH(:)
      end if
      allocate(bEns(instanceIndex)%scaleFactor_M(bEns(instanceIndex)%nLevEns_M))
      allocate(bEns(instanceIndex)%scaleFactor_T(bEns(instanceIndex)%nLevEns_T))
      do levIndex = 1, bEns(instanceIndex)%nLevEns_T
        if (bEns(instanceIndex)%scaleFactor(levIndex) > 0.0d0) then 
          bEns(instanceIndex)%scaleFactor(levIndex) = sqrt(bEns(instanceIndex)%scaleFactor(levIndex))
        else
          bEns(instanceIndex)%scaleFactor(levIndex) = 0.0d0
        end if
      end do
      bEns(instanceIndex)%scaleFactor_T(1:bEns(instanceIndex)%nLevEns_T) = bEns(instanceIndex)%scaleFactor(1:bEns(instanceIndex)%nLevEns_T)
      if (bEns(instanceIndex)%vco_anl%Vcode == 5002) then
        bEns(instanceIndex)%scaleFactor_M(1:bEns(instanceIndex)%nLevEns_M) = bEns(instanceIndex)%scaleFactor(2:(bEns(instanceIndex)%nLevEns_M+1))
      else
        bEns(instanceIndex)%scaleFactor_M(1:bEns(instanceIndex)%nLevEns_M) = bEns(instanceIndex)%scaleFactor(1:bEns(instanceIndex)%nLevEns_M)
      end if

      do levIndex = 1, bEns(instanceIndex)%nLevEns_T
        if (bEns(instanceIndex)%scaleFactorHumidity(levIndex) > 0.0d0) then 
          bEns(instanceIndex)%scaleFactorHumidity(levIndex) = sqrt(bEns(instanceIndex)%scaleFactorHumidity(levIndex))
        else
          bEns(instanceIndex)%scaleFactorHumidity(levIndex) = 0.0d0
        end if
      end do
      
      bEns(instanceIndex)%scaleFactor_SF = bEns(instanceIndex)%scaleFactor_T(bEns(instanceIndex)%nLevEns_T)

    else if (bEns(instanceIndex)%nLevEns_DP > 0) then
      ! Ocean variables on depth levels
      write(*,*) 'nlev_Depth=', bEns(instanceIndex)%nLevEns_DP
      bEns(instanceIndex)%varNameALFA(:) = varNameALFAocean(:)
      allocate(bEns(instanceIndex)%scaleFactor_DP(bEns(instanceIndex)%nLevEns_DP))
      do levIndex = 1, bEns(instanceIndex)%nLevEns_DP
        if (bEns(instanceIndex)%scaleFactor(levIndex) > 0.0d0) then 
          bEns(instanceIndex)%scaleFactor(levIndex) = sqrt(bEns(instanceIndex)%scaleFactor(levIndex))
        else
          bEns(instanceIndex)%scaleFactor(levIndex) = 0.0d0
        end if
      end do
      bEns(instanceIndex)%scaleFactor_DP(1:bEns(instanceIndex)%nLevEns_DP) = bEns(instanceIndex)%scaleFactor(1:bEns(instanceIndex)%nLevEns_DP)
    else
      ! 2D surface variables
      bEns(instanceIndex)%varNameALFA(:) = varNameALFAsfc(:)
      if (bEns(instanceIndex)%scaleFactor(1) > 0.0d0) then 
        bEns(instanceIndex)%scaleFactor_SF = sqrt(bEns(instanceIndex)%scaleFactor(1))
      else
        call utl_abort('ben_setupOneInstance: with vCode == 0, the scale factor should never be equal to 0')
      end if
    end if

    !- 1.5 Domain Partionning
    call mmpi_setup_latbands(bEns(instanceIndex)%nj, latPerPE, latPerPEmax, bEns(instanceIndex)%myLatBeg, bEns(instanceIndex)%myLatEnd)
    call mmpi_setup_lonbands(bEns(instanceIndex)%ni, lonPerPE, lonPerPEmax, bEns(instanceIndex)%myLonBeg, bEns(instanceIndex)%myLonEnd)

    !- 1.6 Localization
    if ( trim(ben_mode) == 'Analysis' ) then

      call mmpi_setup_levels(bEns(instanceIndex)%nEns,myMemberBeg,myMemberEnd,myMemberCount)
      call rpn_comm_allreduce(myMemberCount, maxMyMemberCount, &
           1,"MPI_INTEGER","MPI_MAX","GRID",ierr)
      bEns(instanceIndex)%nEnsOverDimension = mmpi_npex * maxMyMemberCount

      select case(trim(bEns(instanceIndex)%localizationType))
      case('LevelDependent')
        if (mmpi_myid == 0) write(*,*)
        if (mmpi_myid == 0) write(*,*) 'ben_setupOneInstance: Level-Dependent (Standard) localization will be used'
        bEns(instanceIndex)%nWaveBand = 1

      case('ScaleDependent','ScaleDependentWithSpectralLoc')
        if (mmpi_myid == 0) write(*,*)
        if (trim(bEns(instanceIndex)%localizationType) == 'ScaleDependent') then
          if (mmpi_myid == 0) write(*,*) 'ben_setupOneInstance: Scale-Dependent localization (SDL) will be used'
          bEns(instanceIndex)%nWaveBand             = count(bEns(instanceIndex)%waveBandPeaks >= 0)
          bEns(instanceIndex)%nWaveBandForFiltering = bEns(instanceIndex)%nWaveBand
        else
          if (mmpi_myid == 0) write(*,*) 'ben_setupOneInstance: Scale-Dependent localization with Spectral localization (SDLwSL) will be used'
          bEns(instanceIndex)%nWaveBand             = 1
          bEns(instanceIndex)%nWaveBandForFiltering = count(bEns(instanceIndex)%waveBandPeaks >= 0)
        end if

        if ( bEns(instanceIndex)%nWaveBandForFiltering <= 1 ) then
          call utl_abort('ben_setupOneInstance: nWaveBandForFiltering <= 1')
        end if
        ! You must provide nWaveBand wavenumbers in decreasing order
        ! e.g. For a 3 wave bands decomposition...
        !      wavenumber #1 = where the response function for wave band 1 (hgh res) reaches 1 
        !                      and stays at 1 for higher wavenumbers
        !      wavenumber #2 = where the response function for wave band 2 reaches 1
        !      wavenumber #3 = where the response function for wave band 3 (low res) reaches 1 
        !                      and stays at 1 for lower wavenumbers
        ! See FilterResponseFunction for further info...

        ! Make sure that the wavenumbers are in the correct (decreasing) order
        do waveBandIndex = 1, bEns(instanceIndex)%nWaveBandForFiltering-1
          if ( bEns(instanceIndex)%waveBandPeaks(waveBandIndex)-bEns(instanceIndex)%waveBandPeaks(waveBandIndex+1) <= 0 ) then
            call utl_abort('ben_setupOneInstance: waveBandPeaks are not in decreasing wavenumber order') 
          end if
        end do

        if (trim(bEns(instanceIndex)%localizationType) == 'ScaleDependent') then
          ! Make sure that we have valid localization length scales for each wave bands
          do  waveBandIndex = 1, bEns(instanceIndex)%nWaveBand
            if ( bEns(instanceIndex)%hLocalize(waveBandIndex) <= 0.0d0 ) then
              call utl_abort('ben_setupOneInstance: Invalid HORIZONTAL localization length scale')
            end if
            if ( bEns(instanceIndex)%vLocalize(waveBandIndex) <= 0.0d0 .and. (bEns(instanceIndex)%nLevInc_M > 1 .or. bEns(instanceIndex)%nLevInc_T > 1) ) then
              call utl_abort('ben_setupOneInstance: Invalid VERTICAL localization length scale')
            end if
          end do

          ! Make sure the truncation is compatible with the waveBandPeaks
          if ( bEns(instanceIndex)%nTrunc < bEns(instanceIndex)%waveBandPeaks(1) ) then
            call utl_abort('ben_setupOneInstance: The truncation is not compatible with the your scale-dependent localization')
          end if

        else
          ! Do we have a valid selected waveBand index?
          if (bEns(instanceIndex)%waveBandIndexSelected < 1                              .and. &
              bEns(instanceIndex)%waveBandIndexSelected > bEns(instanceIndex)%nWaveBand ) then
            write(*,*) 'ben_setupOneInstance: waveBandIndexSelected = ', bEns(instanceIndex)%waveBandIndexSelected
            write(*,*) 'ben_setupOneInstance: nWaveBand             = ', bEns(instanceIndex)%nWaveBand
            call utl_abort('ben_setupOneInstance: The selected waveBand index is not valid')
          end if

          ! Make sure we have only ONE localization length scales for the selected waveBand index
          if ( bEns(instanceIndex)%hLocalize(2) > 0.0d0 .or. bEns(instanceIndex)%vLocalize(2) > 0.0d0 ) then
            call utl_abort('ben_setupOneInstance: only a single localization length scale must be provided with SDLwSL')
          end if
        end if

      case default
        call utl_abort('ben_setupOneInstance: Invalid mode for LocalizationType')
      end select

      ! Setup the localization
      if ( bEns(instanceIndex)%vco_anl%Vcode == 5002 .or. bEns(instanceIndex)%vco_anl%Vcode == 5005 ) then
        pSurfRef = 101000.D0
        nullify(vertLocationInc)
        status = vgd_levels( bEns(instanceIndex)%vco_anl%vgrid, &
                             ip1_list=bEns(instanceIndex)%vco_anl%ip1_M, &
                             levels=vertLocationInc, &
                             sfc_field=pSurfRef, in_log=.false.)
        if (status /= VGD_OK)then
          call utl_abort('ben_setupOneInstance: ERROR from vgd_levels')
        end if

        allocate(vertLocationEns(bEns(instanceIndex)%nLevEns_M))
        do levIndex = 1, bEns(instanceIndex)%nLevEns_M
          vertLocationEns(levIndex) = log(vertLocationInc(levIndex+bEns(instanceIndex)%topLevIndex_M-1))
        end do
        deallocate(vertLocationInc)
      else if ( bEns(instanceIndex)%vco_anl%nLev_depth > 0 ) then
        allocate(vertLocationEns(bEns(instanceIndex)%vco_anl%nLev_depth))
        vertLocationEns(:) = bEns(instanceIndex)%vco_anl%depths(:)
      else
        pSurfRef = 101000.D0
        allocate(vertLocationEns(1))
        vertLocationEns(:) = pSurfRef
      end if

      allocate(bEns(instanceIndex)%locStorage(bEns(instanceIndex)%nWaveBand))
      do waveBandIndex = 1, bEns(instanceIndex)%nWaveBand
        call loc_setup(bEns(instanceIndex)%locStorage(waveBandIndex), bEns(instanceIndex)%cvDim_mpilocal,          & ! OUT
                       bEns(instanceIndex)%hco_ens, bEns(instanceIndex)%vco_ens, bEns(instanceIndex)%nEns,         & ! IN
                       vertLocationEns, bEns(instanceIndex)%nTrunc, 'spectral',                               & ! IN
                       bEns(instanceIndex)%localizationType, bEns(instanceIndex)%hLocalize(waveBandIndex),         & ! IN
                       bEns(instanceIndex)%hLocalize(waveBandIndex+1), bEns(instanceIndex)%vLocalize(waveBandIndex)) ! IN
      end do

      cvDim = bEns(instanceIndex)%cvDim_mpilocal
      deallocate(vertLocationEns)

    end if

    !- 1.7 Control variables
    if      ( bEns(instanceIndex)%ctrlVarHumidity == 'LQ' ) then
      write(*,*)
      write(*,*) 'ben_setupOneInstance: Humidity control variable = ', bEns(instanceIndex)%ctrlVarHumidity
      bEns(instanceIndex)%gsvHUcontainsLQ = .true.
    else if ( bEns(instanceIndex)%ctrlVarHumidity == 'HU' ) then
      write(*,*)
      write(*,*) 'ben_setupOneInstance: Humidity control variable = ', bEns(instanceIndex)%ctrlVarHumidity
      bEns(instanceIndex)%gsvHUcontainsLQ = .false.
    else
      write(*,*)
      write(*,*) 'Unknown humidity control variable'
      write(*,*) 'Should be LQ or LU, found = ', bEns(instanceIndex)%ctrlVarHumidity
      call utl_abort('ben_setupOneInstance')
    end if

    !
    !- 2.  Read/Process the Ensemble
    !
    
    !- 2.1 Identify set of variables for which ensembles are required    
    do jvar = 1, vnl_numvarmax
      if (trim(bEns(instanceIndex)%includeAnlVar(jvar)) == '') exit
      if (.not.gsv_varExist(varName = trim(bEns(instanceIndex)%includeAnlVar(jvar)))) then
        write(*,*) 'ben_setupOneInstance: This variable is not a member of ANLVAR: ', trim(bEns(instanceIndex)%includeAnlVar(jvar))
        call utl_abort('ben_setupOneInstance: Invalid variable in includeAnlVar')
      else
        bEns(instanceIndex)%numIncludeAnlVar = bEns(instanceIndex)%numIncludeAnlVar+1
      end if
    end do
    if (bEns(instanceIndex)%numIncludeAnlVar == 0) then
      do jvar = 1, vnl_numvarmax
        if (.not.gsv_varExist(varName = trim(vnl_varNamelist(jvar)))) cycle
        bEns(instanceIndex)%numIncludeAnlVar = bEns(instanceIndex)%numIncludeAnlVar+1
        bEns(instanceIndex)%includeAnlVar(bEns(instanceIndex)%numIncludeAnlVar) = vnl_varNamelist(jvar)
      end do
    end if

    if (bEns(instanceIndex)%numIncludeAnlVar == 0) call utl_abort('ben_setupOneInstance: Ensembles not being requested for any variable')

    bEns(instanceIndex)%ensShouldNotContainLQvarName=.false.
    if (bEns(instanceIndex)%ctrlVarHumidity == 'LQ' .and. .not. bEns(instanceIndex)%ensContainsFullField) then
      ! In this particular case, we must force readEnsemble to contains the LQ varName 
      ! to be able to read LQ perturbations
      do jvar = 1, bEns(instanceIndex)%numIncludeAnlVar
        if (bEns(instanceIndex)%includeAnlVar(jvar) == 'LQ')  then
          call utl_abort('ben_setup: LQ must not be present in ANLVAR in this case')
        end if
        if (bEns(instanceIndex)%includeAnlVar(jvar) == 'HU')  then 
          bEns(instanceIndex)%includeAnlVar(jvar) = 'LQ'
          bEns(instanceIndex)%ensShouldNotContainLQvarName=.true.
        end if
      end do
    end if

    !- 2.2 Read the ensemble data
    call setupEnsemble(instanceIndex)

    if ( trim(ben_mode) /= 'Analysis' ) then
      cvDim = 9999 ! Dummy value > 0 to indicate to the background check (s/r ose_compute_HBHT_ensemble)
      return
    end if

    !- 2.3 Convert into a C matrix
    if (bEns(instanceIndex)%useCmatrixOnly .or. trim(bEns(instanceIndex)%varianceSmoothing) /= 'none') then
      call ens_computeStdDev(bEns(instanceIndex)%ensPerts(1), containsScaledPerts_opt=.true.)
      call ens_normalize(bEns(instanceIndex)%ensPerts(1))
      bEns(instanceIndex)%ensPertsNormalized = .true.
    else
      bEns(instanceIndex)%ensPertsNormalized = .false.
    end if

    !- 2.4 Variance smoothing
    if (trim(bEns(instanceIndex)%varianceSmoothing) /= 'none' .and. .not. bEns(instanceIndex)%useCmatrixOnly) then
      if (mmpi_myid == 0) write(*,*) 'ben_setupOneInstance: variance smoothing will be performed'

      call ens_copyEnsStdDev(bEns(instanceIndex)%ensPerts(1), bEns(instanceIndex)%statevector_ensStdDev)
      if (bEns(instanceIndex)%ensDiagnostic) then
        call gio_writeToFile(bEns(instanceIndex)%statevector_ensStdDev,'./ens_stddev.fst',       & ! IN
                             'STDDEV_RAW', HUcontainsLQ_opt=bEns(instanceIndex)%gsvHUcontainsLQ)  ! IN
      end if

      call gsv_power(bEns(instanceIndex)%statevector_ensStdDev, 2.d0 ) ! StdDev -> Variance

      if (trim(bEns(instanceIndex)%varianceSmoothing) == 'horizMean') then
        if (mmpi_myid == 0) write(*,*) 'ben_setup: variance smoothing type = horizMean'
        call gbi_setup(gbi_horizontalMean, 'HorizontalMean', &
                       bEns(instanceIndex)%statevector_ensStdDev, bEns(instanceIndex)%hco_core)
        call gbi_mean(gbi_horizontalMean, bEns(instanceIndex)%statevector_ensStdDev, & ! IN
                      bEns(instanceIndex)%statevector_ensStdDev)                       ! OUT
        call gbi_deallocate(gbi_horizontalMean)
      else if (trim(bEns(instanceIndex)%varianceSmoothing) == 'footprint') then
        call gsv_smoothHorizontal(bEns(instanceIndex)%statevector_ensStdDev, & ! INOUT
                                  bEns(instanceIndex)%footprintRadius)         ! IN
      else if (trim(bEns(instanceIndex)%varianceSmoothing) == 'footprintLandSeaTopo') then
        call gbi_setup(gbi_landSeaTopo, 'landSeaTopo', bEns(instanceIndex)%statevector_ensStdDev, &
                       bEns(instanceIndex)%hco_core, &
                       mpi_distribution_opt='None', writeBinsToFile_opt=bEns(instanceIndex)%ensDiagnostic)
        call gsv_getField(gbi_landSeaTopo%statevector_bin2d,bin2d)
        HeightSfc => gsv_getHeightSfc(gbi_landSeaTopo%statevector_bin2d)
        call gsv_smoothHorizontal(bEns(instanceIndex)%statevector_ensStdDev,                                       & ! INOUT
                                  bEns(instanceIndex)%footprintRadius, binInteger_opt=bin2d, binReal_opt=HeightSfc,& ! IN
                                  binRealThreshold_opt=bEns(instanceIndex)%footprintTopoThreshold)                   ! IN
        call gbi_deallocate(gbi_landSeaTopo)
      else
        call utl_abort('ben_setupOneInstance: Invalid variance smoothing type = '//trim(bEns(instanceIndex)%varianceSmoothing))
      end if

      call gsv_power(bEns(instanceIndex)%statevector_ensStdDev, 0.5d0) ! Variance -> StdDev

      if (bEns(instanceIndex)%ensDiagnostic) then
        call gio_writeToFile(bEns(instanceIndex)%statevector_ensStdDev,'./ens_stddev.fst',& ! IN
                             'STDDEV_SMOOT', HUcontainsLQ_opt=bEns(instanceIndex)%gsvHUcontainsLQ)   ! IN
      end if
    end if

    !- 2.5 Pre-compute everything for advection in FSO mode
    if (bEns(instanceIndex)%fsoLeadTime > 0.0D0) then
      bEns(instanceIndex)%amp3dStepIndexFSOFcst = 1
      if ( sum(bEns(instanceIndex)%advectFactorFSOFcst(:)) == 0.0D0 .or. bEns(instanceIndex)%numStep == 1) then
        if (mmpi_myid == 0) write(*,*) 'ben_setupOneInstance: advection not activated for FSO'
        bEns(instanceIndex)%advectAmplitudeFSOFcst = .false.
        bEns(instanceIndex)%numStepAmplitudeFSOFcst = 1
      else
        if (mmpi_myid == 0) write(*,*) 'ben_setupOneInstance: advection activated in FSO mode'
        bEns(instanceIndex)%advectAmplitudeFSOFcst = .true.
        bEns(instanceIndex)%numStepAmplitudeFSOFcst = 2
        bEns(instanceIndex)%numStepAdvectFSOFcst = nint(bEns(instanceIndex)%fsoLeadTime/6.0D0) + 1
        allocate(advectFactorFSOFcst_M(bEns(instanceIndex)%vco_ens%nLev_M))
        advectFactorFSOFcst_M(:) = bEns(instanceIndex)%advectFactorFSOFcst(1:bEns(instanceIndex)%vco_ens%nLev_M)
        allocate(bEns(instanceIndex)%dateStampListAdvectedFields(bEns(instanceIndex)%numStepAmplitudeFSOFcst))
        bEns(instanceIndex)%dateStampListAdvectedFields(1) = tim_getDatestamp()
        bEns(instanceIndex)%dateStampListAdvectedFields(2) = bEns(instanceIndex)%dateStampList(bEns(instanceIndex)%numStep)
        delT_hour = bEns(instanceIndex)%fsoLeadTime/real(bEns(instanceIndex)%numStepAdvectFSOFcst-1,8) ! time between winds
        call utl_tmg_start(55,'------B_ENS_SetupAdvecFSO')
        call adv_setup( bEns(instanceIndex)%adv_amplitudeFSOFcst,                                   & ! OUT
                        'fromFirstTimeIndex', bEns(instanceIndex)%hco_ens, bEns(instanceIndex)%vco_ens,                 & ! IN
                        bEns(instanceIndex)%numStepAmplitudeFSOFcst, bEns(instanceIndex)%dateStampListAdvectedFields,   & ! IN
                        bEns(instanceIndex)%numStepAdvectFSOFcst, delT_hour, advectFactorFSOFcst_M, & ! IN
                        'MMLevsOnly',                                           & ! IN
                        steeringFlowFilename_opt=trim(bEns(instanceIndex)%ensPathName)//'/forecast_for_advection' ) ! IN
        call utl_tmg_stop(55)
        deallocate(advectFactorFSOFcst_M)
      end if
    else
      bEns(instanceIndex)%numStepAmplitudeFSOFcst = 0
    end if

    !- 2.6 Pre-compute everything for advection in ANALYSIS mode
    if ( sum(bEns(instanceIndex)%advectFactorAssimWindow(:)) == 0.0D0 .or. bEns(instanceIndex)%numStep == 1) then
      if (mmpi_myid == 0) write(*,*) 'ben_setupOneInstance: advection not activated in ANALYSIS mode'

      bEns(instanceIndex)%advectAmplitudeAssimWindow = .false.
      bEns(instanceIndex)%numStepAmplitudeAssimWindow = 1
      bEns(instanceIndex)%amp3dStepIndexAssimWindow   = 1

    else
      if (mmpi_myid == 0) write(*,*) 'ben_setupOneInstance: advection activated in ANALYSIS mode'

      delT_hour                 = tim_dstepobsinc
      allocate(advectFactorAssimWindow_M(bEns(instanceIndex)%vco_ens%nLev_M))
      advectFactorAssimWindow_M(:) = bEns(instanceIndex)%advectFactorAssimWindow(1:bEns(instanceIndex)%vco_ens%nLev_M)
      allocate(bEns(instanceIndex)%dateStampListAdvectedFields(bEns(instanceIndex)%numStep))
      bEns(instanceIndex)%dateStampListAdvectedFields(:) = bEns(instanceIndex)%dateStampList(:)
      call gsv_allocate(statevector_ensMean4D, bEns(instanceIndex)%numStep, bEns(instanceIndex)%hco_ens, bEns(instanceIndex)%vco_ens, &
                        datestampList_opt=bEns(instanceIndex)%dateStampListAdvectedFields,     &
                        mpi_local_opt=.true., &
                        allocHeight_opt=.false., allocPressure_opt=.false.)
      call ens_copyEnsMean(bEns(instanceIndex)%ensPerts(1), & ! IN
                           statevector_ensMean4D  )   ! OUT

      call utl_tmg_start(56,'------B_ENS_SetupAdvecAnl')

      select case(trim(bEns(instanceIndex)%advectTypeAssimWindow))
      case ('amplitude')
        if (mmpi_myid == 0) write(*,*) '         amplitude fields will be advected'
        bEns(instanceIndex)%advectAmplitudeAssimWindow  = .true.
        bEns(instanceIndex)%numStepAmplitudeAssimWindow = bEns(instanceIndex)%numStep
        bEns(instanceIndex)%numStepAdvectAssimWindow    = bEns(instanceIndex)%numStep

        select case(trim(bEns(instanceIndex)%advectStartTimeIndexAssimWindow))
        case ('first')
          direction='fromFirstTimeIndex'
          bEns(instanceIndex)%amp3dStepIndexAssimWindow = 1
        case ('middle')
          direction='fromMiddleTimeIndex'
          bEns(instanceIndex)%amp3dStepIndexAssimWindow = (bEns(instanceIndex)%numStepAmplitudeAssimWindow+1)/2
        case default
          write(*,*)
          write(*,*) 'Unsupported starting timeIndex : ', trim(bEns(instanceIndex)%advectStartTimeIndexAssimWindow)
          call utl_abort('ben_setupOneInstance')
        end select
        call adv_setup( bEns(instanceIndex)%adv_amplitudeAssimWindow,                                          & ! OUT
                        direction, bEns(instanceIndex)%hco_ens, bEns(instanceIndex)%vco_ens,                                       & ! IN
                        bEns(instanceIndex)%numStepAmplitudeAssimWindow, bEns(instanceIndex)%dateStampListAdvectedFields,          & ! IN
                        bEns(instanceIndex)%numStepAdvectAssimWindow, delT_hour, advectFactorAssimWindow_M,    & ! IN
                        'MMLevsOnly', statevector_steeringFlow_opt = statevector_ensMean4D ) ! IN

      case('ensPertAnlInc')
        if (mmpi_myid == 0) write(*,*) '         ensPerts and AnalInc will be advected'

        if (.not. EnsTopMatchesAnlTop) then
          call utl_abort('ben_setupOneInstance: for advectTypeAssimWindow=ensPertAnlInc, ensTop and anlTop must match!')
        end if

        bEns(instanceIndex)%advectEnsPertAnlInc         = .true.
        bEns(instanceIndex)%amp3dStepIndexAssimWindow   = 1
        bEns(instanceIndex)%numStepAmplitudeAssimWindow = 1
        bEns(instanceIndex)%numStepAdvectAssimWindow    = bEns(instanceIndex)%numStep
        
        select case(trim(bEns(instanceIndex)%advectStartTimeIndexAssimWindow))
        case ('first')
          directionEnsPerts='towardFirstTimeIndex'
          directionAnlInc  ='towardFirstTimeIndexInverse'
        case ('middle')
          directionEnsPerts='towardMiddleTimeIndex'
          directionAnlInc  ='towardMiddleTimeIndexInverse'
        case default
          write(*,*)
          write(*,*) 'Unsupported starting timeIndex : ', trim(bEns(instanceIndex)%advectStartTimeIndexAssimWindow)
          call utl_abort('ben_setupOneInstance')
        end select

        call adv_setup( bEns(instanceIndex)%adv_ensPerts,                                                              & ! OUT
                        directionEnsPerts, bEns(instanceIndex)%hco_ens, bEns(instanceIndex)%vco_ens,                   & ! IN
                        bEns(instanceIndex)%numStepAdvectAssimWindow, bEns(instanceIndex)%dateStampListAdvectedFields, & ! IN
                        bEns(instanceIndex)%numStepAdvectAssimWindow, delT_hour, advectFactorAssimWindow_M,            & ! IN
                        'allLevs', statevector_steeringFlow_opt=statevector_ensMean4D )   ! IN

        call adv_setup( bEns(instanceIndex)%adv_analInc,                                                               & ! OUT
                        directionAnlInc, bEns(instanceIndex)%hco_ens, bEns(instanceIndex)%vco_ens,                     & ! IN
                        bEns(instanceIndex)%numStepAdvectAssimWindow, bEns(instanceIndex)%dateStampListAdvectedFields, & ! IN
                        bEns(instanceIndex)%numStepAdvectAssimWindow, delT_hour, advectFactorAssimWindow_M,            & ! IN
                        'allLevs', statevector_steeringFlow_opt=statevector_ensMean4D )   ! IN

      case default
        write(*,*)
        write(*,*) 'Unsupported advectTypeAssimWindow : ', trim(bEns(instanceIndex)%advectTypeAssimWindow)
        call utl_abort('ben_setupOneInstance')
      end select

      call utl_tmg_stop(56)

      deallocate(advectFactorAssimWindow_M)

      !- If wanted, write the ensemble mean
      if (bEns(instanceIndex)%advDiagnostic) then
        do stepIndex = 1, tim_nstepobsinc
          call gio_writeToFile(statevector_ensMean4D,'./ens_mean.fst','ENSMEAN4D',        & ! IN
                               stepIndex_opt=stepIndex, HUcontainsLQ_opt=bEns(instanceIndex)%gsvHUcontainsLQ ) ! IN
        end do
      end if

      call gsv_deallocate(statevector_ensMean4D)

    end if

    !- 2.7 Ensemble perturbations advection
    if ( bEns(instanceIndex)%advectEnsPertAnlInc ) then

      !- If wanted, write the original ensemble perturbations for member #1
      if (bEns(instanceIndex)%advDiagnostic) then
        call ens_copyMember(bEns(instanceIndex)%ensPerts(1), statevector_oneEnsPert4D, 1)
        do stepIndex = 1, tim_nstepobsinc
          call gio_writeToFile(statevector_oneEnsPert4D,'./ens_pert1.fst','ORIGINAL', & ! IN
               stepIndex_opt=stepIndex, HUcontainsLQ_opt=bEns(instanceIndex)%gsvHUcontainsLQ )             ! IN
        end do
      end if

      !- Do the advection of all the members
      call adv_ensemble_tl( bEns(instanceIndex)%ensPerts(1), &       ! INOUT
                            bEns(instanceIndex)%adv_ensPerts, bEns(instanceIndex)%nEns ) ! IN

      !- If wanted, write the advected ensemble perturbations for member #1
      if (bEns(instanceIndex)%advDiagnostic) then
        call ens_copyMember(bEns(instanceIndex)%ensPerts(1), statevector_oneEnsPert4D, 1)
        do stepIndex = 1, tim_nstepobsinc
          call gio_writeToFile(statevector_oneEnsPert4D,'./ens_pert1_advected.fst','ADVECTED', & ! IN
               stepIndex_opt=stepIndex,HUcontainsLQ_opt=bEns(instanceIndex)%gsvHUcontainsLQ )                       ! IN
        end do
      end if

    end if

    !- 2.8 Compute and write Std. Dev.
    if (bEns(instanceIndex)%ensDiagnostic) call ensembleDiagnostic(instanceIndex,'FullPerturbations')

    !- 2.9 Partitioned the ensemble perturbations into wave bands
    if (trim(bEns(instanceIndex)%localizationType) == 'ScaleDependent'                .or. &
        trim(bEns(instanceIndex)%localizationType) == 'ScaleDependentWithSpectralLoc') then
      call ensembleScaleDecomposition(instanceIndex)
      if (bEns(instanceIndex)%ensDiagnostic) call ensembleDiagnostic(instanceIndex,'WaveBandPerturbations')
    end if

    ! Allocate main ensemble Amplitude arrays for improved efficiency - do not know why!!!
    if (bEns(instanceIndex)%numStepAmplitudeAssimWindow > 1 .or. bEns(instanceIndex)%numStepAmplitudeFSOFcst > 1 ) then
      write(*,*) 'ben_setupOneInstance: WARNING: due to advection being activated, cannot currently used saved ensemble amplitudes!'
      bEns(instanceIndex)%useSaveAmp = .false.
    else
      write(*,*) 'ben_setupOneInstance: using saved memory for ensemble amplitudes for improved efficiency'
      bEns(instanceIndex)%useSaveAmp = .true.
      call ens_allocate(ensAmplitudeSave(instanceIndex),                 &
                        bEns(instanceIndex)%nEnsOverDimension,           &
                        bEns(instanceIndex)%numStepAmplitudeAssimWindow, &
                        bEns(instanceIndex)%hco_ens,                     &
                        bEns(instanceIndex)%vco_ens,                     &
                        bEns(instanceIndex)%dateStampList,               &
                        hco_core_opt=bEns(instanceIndex)%hco_core,       &
                        varNames_opt=bEns(instanceIndex)%varNameALFA, dataKind_opt=8)
      call ens_zero(ensAmplitudeSave(instanceIndex))
    end if

    !- 2.10 Setup en ensGridStateVector to store the amplitude fields (for writing)
    if (bEns(instanceIndex)%keepAmplitude) then
      write(*,*)
      write(*,*) 'ben_setupOneInstance: ensAmplitude fields will be store for potential write to file'
      call ens_allocate(bEns(instanceIndex)%ensAmplitudeStorage, bEns(instanceIndex)%nEns, &
                        bEns(instanceIndex)%numStepAmplitudeAssimWindow, &
                        bEns(instanceIndex)%hco_ens, bEns(instanceIndex)%vco_ens, &
                        bEns(instanceIndex)%dateStampList, &
                        hco_core_opt=bEns(instanceIndex)%hco_core, &
                        varNames_opt=bEns(instanceIndex)%varNameALFA, dataKind_opt=8)
    end if

  end subroutine ben_setupOneInstance

  !--------------------------------------------------------------------------
  ! ben_finalize
  !--------------------------------------------------------------------------
  subroutine ben_finalize()
    implicit none
    integer :: waveBandIndex
    integer :: instanceIndex

    if (verbose) write(*,*) 'Entering ben_Finalize'

    do instanceIndex = 1, nInstance
      if (bEns(instanceIndex)%initialized) then
        write(*,*) 'ben_finalize: deallocating B_ensemble arrays for instance #', instanceIndex
        do waveBandIndex = 1, bEns(instanceIndex)%nWaveBand
          call ens_deallocate(bEns(instanceIndex)%ensPerts(waveBandIndex))
          call loc_finalize(bEns(instanceIndex)%locStorage(waveBandIndex))
        end do
        deallocate(bEns(instanceIndex)%ensPerts)
        if (bEns(instanceIndex)%keepAmplitude) call ens_deallocate(bEns(instanceIndex)%ensAmplitudeStorage)
      end if
    end do
  
  end subroutine ben_finalize

  !--------------------------------------------------------------------------
  ! ben_getScaleFactor
  !--------------------------------------------------------------------------
  subroutine ben_getScaleFactor(scaleFactor_out, instanceIndex_opt)
    implicit none
    real(8), intent(out)          :: scaleFactor_out(:)
    integer, optional, intent(in) :: instanceIndex_opt

    integer :: levIndex, instanceIndex

    if (verbose) write(*,*) 'Entering ben_getScaleFactor'

    instanceIndex = ben_setInstanceIndex(instanceIndex_opt)

    ! return value of 0 above highest level of ensemble
    do levIndex = 1, (bEns(instanceIndex)%topLevIndex_T - 1)
      scaleFactor_out(levIndex) = 0.0d0
    end do
    ! return scale factor for thermo levels
    do levIndex = bEns(instanceIndex)%topLevIndex_T, bEns(instanceIndex)%nLevInc_T
      scaleFactor_out(levIndex) = bEns(instanceIndex)%scaleFactor_T(levIndex-bEns(instanceIndex)%topLevIndex_T+1)
    end do

  end subroutine ben_getScaleFactor

  !--------------------------------------------------------------------------
  ! ben_setInstanceIndex
  !--------------------------------------------------------------------------
  function ben_setInstanceIndex(instanceIndex_opt) result(instanceIndex)
    !
    !:Purpose: To return the appropriate instance index
    !
    implicit none
    integer, optional :: instanceIndex_opt
    integer :: instanceIndex

    if (present(instanceIndex_opt)) then
      instanceIndex = instanceIndex_opt
    else
      instanceIndex = 1
    end if

  end function ben_setInstanceIndex

  !--------------------------------------------------------------------------
  ! ben_getnEns
  !--------------------------------------------------------------------------
  integer function ben_getnEns(instanceIndex_opt)
    !
    !:Purpose: To return the number of ensemble members
    !
    implicit none
    integer, optional, intent(in) :: instanceIndex_opt
    integer :: instanceIndex

    instanceIndex = ben_setInstanceIndex(instanceIndex_opt)

    ben_getnEns = bEns(instanceIndex)%nEns

  end function ben_getnEns

  !--------------------------------------------------------------------------
  ! setupEnsemble
  !--------------------------------------------------------------------------
  subroutine setupEnsemble(instanceIndex)
    implicit none

    integer, intent(in) :: instanceIndex

    real(4), pointer     :: ptr4d_r4(:,:,:,:)
    real(8) :: multFactor
    integer :: stepIndex,levIndex,lev,waveBandIndex,memberIndex,varIndex
    logical :: makeBiPeriodic
    character(len=4)  :: varName
    character(len=30) :: transform

    write(*,*) 'setupEnsemble: Start'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !- 1. Memory allocation
    allocate(bEns(instanceIndex)%ensPerts(bEns(instanceIndex)%nWaveBand))    
    do waveBandIndex = 1, bEns(instanceIndex)%nWaveBand
      call ens_allocate(bEns(instanceIndex)%ensPerts(waveBandIndex),  &
                        bEns(instanceIndex)%nEns, bEns(instanceIndex)%numStep, &
                        bEns(instanceIndex)%hco_ens,  &
                        bEns(instanceIndex)%vco_ens, bEns(instanceIndex)%dateStampList, &
                        hco_core_opt = bEns(instanceIndex)%hco_core, &
                        varNames_opt = bEns(instanceIndex)%includeAnlVar(1:bEns(instanceIndex)%numIncludeAnlVar), &
                        hInterpolateDegree_opt = bEns(instanceIndex)%hInterpolationDegree)
    end do

    !- 2. Read ensemble
    makeBiPeriodic = (trim(bEns(instanceIndex)%localizationType) == 'ScaleDependent' .or. &
                      trim(bEns(instanceIndex)%localizationType) == 'ScaleDependentWithSpectralLoc' .or. &
                      trim(bEns(instanceIndex)%varianceSmoothing) /= 'none')
    call ens_readEnsemble(bEns(instanceIndex)%ensPerts(1), bEns(instanceIndex)%ensPathName, makeBiPeriodic,         &
                          vco_file_opt = bEns(instanceIndex)%vco_file,                          &
                          varNames_opt = bEns(instanceIndex)%includeAnlVar(1:bEns(instanceIndex)%numIncludeAnlVar), & 
                          containsFullField_opt=bEns(instanceIndex)%ensContainsFullField)

    if ( bEns(instanceIndex)%ctrlVarHumidity == 'LQ' .and. ens_varExist(bEns(instanceIndex)%ensPerts(1),'HU') .and. &
         bEns(instanceIndex)%ensContainsFullField ) then
      call gvt_transform(bEns(instanceIndex)%ensPerts(1),'HUtoLQ',huMinValue_opt=bEns(instanceIndex)%huMinValue)
    else if ( bEns(instanceIndex)%ctrlVarHumidity == 'HU' .and. ens_varExist(bEns(instanceIndex)%ensPerts(1),'HU') .and. &
         bEns(instanceIndex)%ensContainsFullField .and. bEns(instanceIndex)%huMinValue /= MPC_missingValue_R8 ) then
      call qlim_setMin(bEns(instanceIndex)%ensPerts(1), bEns(instanceIndex)%huMinValue)
    else if ( trim(bEns(instanceIndex)%transformVarKindCH) /= '' ) then
      do varIndex = 1, bEns(instanceIndex)%numIncludeAnlVar
        if ( vnl_varKindFromVarname(bEns(instanceIndex)%includeAnlVar(varIndex)) /= 'CH' ) cycle            

        transform = trim(bens(instanceIndex)%transformVarKindCH)//'CH'
        call gvt_transform( bEns(instanceIndex)%ensPerts(1), trim(transform), &          
                            varName_opt=bEns(instanceIndex)%includeAnlVar(varIndex) ) 
      end do
    end if

    !- 3. From ensemble FORECASTS to ensemble PERTURBATIONS

    !- 3.1 remove mean
    call ens_computeMean( bEns(instanceIndex)%ensPerts(1), bEns(instanceIndex)%removeSubEnsMeans, numSubEns_opt=bEns(instanceIndex)%numSubEns )
    call ens_removeMean( bEns(instanceIndex)%ensPerts(1) )

    !- 3.2 normalize and apply scale factors
    !$OMP PARALLEL DO PRIVATE (levIndex,varName,lev,ptr4d_r4,stepIndex,memberIndex,multFactor)
    do levIndex = 1, ens_getNumK(bEns(instanceIndex)%ensPerts(1))
      varName = ens_getVarNameFromK(bEns(instanceIndex)%ensPerts(1),levIndex)
      lev = ens_getLevFromK(bEns(instanceIndex)%ensPerts(1),levIndex)

      if ( .not. ens_varExist(bEns(instanceIndex)%ensPerts(1), varName) ) cycle 

      ptr4d_r4 => ens_getOneLev_r4(bEns(instanceIndex)%ensPerts(1),levIndex)

      do stepIndex = 1, bEns(instanceIndex)%numStep
        do memberIndex = 1, bEns(instanceIndex)%nEns

          if ( vnl_varLevelFromVarname(varName) == 'MM' ) then
            multFactor = bEns(instanceIndex)%scaleFactor_M(lev)
          else if ( vnl_varLevelFromVarname(varName) == 'TH' ) then
            multFactor = bEns(instanceIndex)%scaleFactor_T(lev)
          else  if ( vnl_varLevelFromVarname(varName) == 'SF' ) then
            multFactor = bEns(instanceIndex)%scaleFactor_SF
          else if ( vnl_varLevelFromVarname(varName) == 'DP' ) then
            multFactor = bEns(instanceIndex)%scaleFactor_DP(lev)
          else if ( vnl_varLevelFromVarname(varName) == 'SFDP' ) then
            multFactor = bEns(instanceIndex)%scaleFactor_DP(1)
          else
            write(*,*) 'varName = ', varName, ', varLevel = ', vnl_varLevelFromVarname(varName)
            call utl_abort('setupEnsemble: unknown varLevel')
          end if

          multFactor = multFactor/sqrt(1.0d0*dble(bEns(instanceIndex)%nEns-bEns(instanceIndex)%numSubEns))

          if (trim(varName) == 'HU') then
            multFactor = multFactor*bEns(instanceIndex)%scaleFactorHumidity(lev)
          end if

          ptr4d_r4(memberIndex,stepIndex,:,:) = real( real(ptr4d_r4(memberIndex,stepIndex,:,:),8)*multFactor, 4 )

        end do ! memberIndex
      end do ! stepIndex

    end do ! levIndex
    !$OMP END PARALLEL DO

    write(*,*) 'ben_setupEnsemble: finished adjusting ensemble members...'

  end subroutine setupEnsemble

  !--------------------------------------------------------------------------
  ! ben_getPerturbation
  !--------------------------------------------------------------------------
  subroutine ben_getPerturbation(statevector, memberIndexWanted,  &
                                 upwardExtrapolationMethod, waveBandIndexWanted_opt, &
                                 undoNormalization_opt, instanceIndex_opt)
    implicit none

    type(struct_gsv) :: statevector
    integer,          intent(in) :: memberIndexWanted
    character(len=*), intent(in) :: upwardExtrapolationMethod
    integer, optional, intent(in):: waveBandIndexWanted_opt
    logical, optional :: undoNormalization_opt
    integer, optional, intent(in) :: instanceIndex_opt

    integer :: instanceIndex

    real(8), pointer :: ptr4d_r8(:,:,:,:)
    real(4), pointer :: ensOneLev_r4(:,:,:,:)
    real(8) :: dnens2, scaleFactor_MT
    logical :: undoNormalization
    integer :: waveBandIndex
    integer :: lonIndex,latIndex,stepIndex,levIndex,lev,levInc,topLevOffset
    character(len=4) :: varName

    if (verbose) write(*,*) 'Entering ben_getPerturbation'

    instanceIndex = ben_setInstanceIndex(instanceIndex_opt)

    if ( trim(upwardExtrapolationMethod) /= "ConstantValue" ) then
      call utl_abort('ben_getPerturbation : Invalid value for upwardExtrapolationMethod')
    end if

    if ( present(waveBandIndexWanted_opt) ) then
      waveBandIndex = waveBandIndexWanted_opt
    else
      waveBandIndex = 1
    end if

    ! set default value for optional argument undoNormalization
    if ( present(undoNormalization_opt) ) then
      undoNormalization = undoNormalization_opt
    else
      undoNormalization = .false.
    end if

    do levIndex = 1, ens_getNumK(bEns(instanceIndex)%ensPerts(1))
      varName = ens_getVarNameFromK(bEns(instanceIndex)%ensPerts(1),levIndex)
      lev = ens_getLevFromK(bEns(instanceIndex)%ensPerts(1),levIndex)

      if (varName == 'LQ' .and. bEns(instanceIndex)%ensShouldNotContainLQvarName) then
        call gsv_getField(statevector, ptr4d_r8, 'HU')
      else
        call gsv_getField(statevector, ptr4d_r8, varName)
      end if
      ensOneLev_r4 => ens_getOneLev_r4(bEns(instanceIndex)%ensPerts(waveBandIndex),levIndex)

      !$OMP PARALLEL DO PRIVATE(stepIndex,topLevOffset,scaleFactor_MT,levInc,dnens2,latIndex,lonIndex)
      do stepIndex = 1, bEns(instanceIndex)%numStep

        if ( vnl_varLevelFromVarname(varName) == 'MM' ) then
          topLevOffset = bEns(instanceIndex)%topLevIndex_M - 1
          scaleFactor_MT = bEns(instanceIndex)%scaleFactor_M(lev)
        else if ( vnl_varLevelFromVarname(varName) == 'TH' ) then
          topLevOffset = bEns(instanceIndex)%topLevIndex_T - 1
          scaleFactor_MT = bEns(instanceIndex)%scaleFactor_T(lev)
        else ! SF
          topLevOffset = 0
          scaleFactor_MT = bEns(instanceIndex)%scaleFactor_SF
        end if

        levInc = lev + topLevOffset

        ! undo the normalization (optional)
        if (undoNormalization) then
          if (scaleFactor_MT > 0.0d0) then
            dnens2 = sqrt(1.0d0*dble(bEns(instanceIndex)%nEns-1))/scaleFactor_MT
          else
            if (stepIndex == 1) then 
              write(*,*) 'scalefactor not positive, cannot undo normalization!'
              write(*,*) varName,scaleFactor_MT,lev
            end if
            dnens2 = 0.0d0
          end if
          if (varName == 'HU  ') then
            if (bEns(instanceIndex)%scaleFactorHumidity(lev).gt.0.0d0) then
              dnens2 = dnens2/bEns(instanceIndex)%scaleFactorHumidity(lev)
            else
              if (stepIndex == 1) then
                write(*,*) 'Humidity scalefactor not positive, cannot undo normalization!'
                write(*,*) varName,bEns(instanceIndex)%scaleFactorHumidity(lev),lev
              end if
            end if
          end if
        else
          dnens2 = 1.0d0
        end if

        do latIndex = bEns(instanceIndex)%myLatBeg, bEns(instanceIndex)%myLatEnd
          do lonIndex = bEns(instanceIndex)%myLonBeg, bEns(instanceIndex)%myLonEnd
            ptr4d_r8(lonIndex,latIndex,levInc,stepIndex) =   &
                 dnens2*dble(ensOneLev_r4(memberIndexWanted,stepIndex,lonIndex,latIndex))
          end do
        end do

        if ( topLevOffset > 0 .and. lev == 1) then
          ! Fill the gap between the ensemble lid and the analysis lid

          ! undo the normalization (optional)
          if (undoNormalization) then
            if (bEns(instanceIndex)%scaleFactor(1) > 0.0d0) then
              dnens2 = sqrt(1.0d0*dble(bEns(instanceIndex)%nEns-1))/bEns(instanceIndex)%scaleFactor(1)
            else
              if (stepIndex == 1) then
                write(*,*) 'scalefactor(top) not positive, cannot undo normalization!'
                write(*,*) varName,bEns(instanceIndex)%scaleFactor(1)
              end if
              dnens2 = 0.0d0
            end if
            if (varName == 'HU  ') then
              if (bEns(instanceIndex)%scaleFactorHumidity(1) > 0.0d0) then
                dnens2 = dnens2/bEns(instanceIndex)%scaleFactorHumidity(1)
              else
                if (stepIndex == 1) then
                  write(*,*) 'Humidity scalefactor(top) not positive, cannot undo normalization!'
                  write(*,*) varName,bEns(instanceIndex)%scaleFactorHumidity(1)
                end if
              end if
            end if
          else
            dnens2 = 1.0d0
          end if

          do levInc = 1, topLevOffset
            ! using a constant value
            do latIndex = bEns(instanceIndex)%myLatBeg, bEns(instanceIndex)%myLatEnd
              do lonIndex = bEns(instanceIndex)%myLonBeg, bEns(instanceIndex)%myLonEnd
                ptr4d_r8(lonIndex,latIndex,levInc,stepIndex) = dnens2 *  &
                     dble(ensOneLev_r4(memberIndexWanted,stepIndex,lonIndex,latIndex))
              end do
            end do
          end do

        end if ! topLevOffset > 0

      end do ! stepIndex
      !$OMP END PARALLEL DO

    end do ! levIndex

  end subroutine ben_getPerturbation

  !--------------------------------------------------------------------------
  ! ben_getEnsMean
  !--------------------------------------------------------------------------
  subroutine ben_getEnsMean(statevector, upwardExtrapolationMethod, &
                            instanceIndex_opt)
    implicit none

    type(struct_gsv) :: statevector
    character(len=*), intent(in) :: upwardExtrapolationMethod
    integer, optional, intent(in) :: instanceIndex_opt

    real(8), pointer :: ptr4d_out(:,:,:,:)
    real(8), pointer :: ensOneLev_mean(:,:,:)
    integer :: instanceIndex
    integer :: lonIndex,latIndex,stepIndex,levIndex,lev,levInc,topLevOffset
    character(len=4) :: varName

    instanceIndex = ben_setInstanceIndex(instanceIndex_opt)

    if (.not. bEns(instanceIndex)%initialized) then
      if (mmpi_myid == 0) write(*,*) 'ben_getEnsMean: bMatrixEnsemble not initialized, returning zero vector'
      call gsv_zero(statevector)
      return
    end if

    if ( trim(upwardExtrapolationMethod) /= "ConstantValue" ) then
      call utl_abort('ben_getEnsMean : Invalid value for upwardExtrapolationMethod')
    end if

    do levIndex = 1, ens_getNumK(bEns(instanceIndex)%ensPerts(1))
      varName = ens_getVarNameFromK(bEns(instanceIndex)%ensPerts(1),levIndex)
      lev = ens_getLevFromK(bEns(instanceIndex)%ensPerts(1),levIndex)

      if (varName == 'LQ' .and. bEns(instanceIndex)%ensShouldNotContainLQvarName) then
        call gsv_getField(statevector, ptr4d_out, 'HU')
      else
        call gsv_getField(statevector, ptr4d_out, varName)
      end if
      ensOneLev_mean => ens_getOneLevMean_r8(bEns(instanceIndex)%ensPerts(1), 1, levIndex)

      !$OMP PARALLEL DO PRIVATE(stepIndex,topLevOffset,levInc,latIndex,lonIndex)
      do stepIndex = 1, bEns(instanceIndex)%numStep

        if ( vnl_varLevelFromVarname(varName) == 'MM' ) then
          topLevOffset = bEns(instanceIndex)%topLevIndex_M - 1
        else if ( vnl_varLevelFromVarname(varName) == 'TH' ) then
          topLevOffset = bEns(instanceIndex)%topLevIndex_T - 1
        else ! SF
          topLevOffset = 0
        end if

        levInc = lev + topLevOffset

        do latIndex = bEns(instanceIndex)%myLatBeg, bEns(instanceIndex)%myLatEnd
          do lonIndex = bEns(instanceIndex)%myLonBeg, bEns(instanceIndex)%myLonEnd
            ptr4d_out(lonIndex,latIndex,levInc,stepIndex) = ensOneLev_mean(stepIndex,lonIndex,latIndex)
          end do
        end do

        if ( topLevOffset > 0 .and. lev == 1 ) then
          ! Fill the gap between the ensemble lid and the analysis lid

          do levInc = 1, topLevOffset
            ! using a constant value
            do latIndex = bEns(instanceIndex)%myLatBeg, bEns(instanceIndex)%myLatEnd
              do lonIndex = bEns(instanceIndex)%myLonBeg, bEns(instanceIndex)%myLonEnd
                ptr4d_out(lonIndex,latIndex,levInc,stepIndex) = ensOneLev_mean(stepIndex,lonIndex,latIndex)
              end do
            end do
          end do

        end if ! topLevOffset > 0

      end do ! stepIndex
      !$OMP END PARALLEL DO

    end do ! levIndex

  end subroutine ben_getEnsMean

  !--------------------------------------------------------------------------
  ! ensembleScaleDecomposition
  !--------------------------------------------------------------------------
  subroutine ensembleScaleDecomposition(instanceIndex)
    implicit none

    integer, intent(in) :: instanceIndex

    integer :: waveBandIndex, memberindex, stepIndex, levIndex, latIndex, lonIndex
    integer :: ila_filter, p, nla_filter, nphase_filter

    real(8), allocatable :: ResponseFunction(:,:)

    real(8), allocatable :: bandSum(:,:)
    real(8) :: totwvnb_r8

    real(8), allocatable :: ensPertSP(:,:,:)
    real(8), allocatable :: ensPertSPfiltered(:,:,:)
    real(8), allocatable :: ensPertGD(:,:,:)
    real(4), pointer     :: ptr4d_r4(:,:,:,:)

    integer, allocatable :: nIndex_vec(:)

    integer :: nTrunc
    integer :: gstFilterID, mIndex, nIndex, mymBeg, mymEnd, mynBeg, mynEnd, mymSkip, mynSkip
    integer :: mymCount, mynCount
    integer :: waveBandIndexStart, waveBandIndexEnd
    integer :: waveBandIndexLoopStart, waveBandIndexLoopEnd, waveBandIndexLoopDirection

    type(struct_lst)    :: lst_ben_filter ! Spectral transform Parameters for filtering

    character(len=19)   :: kind

    !
    ! --> For SDL
    !
    ! --- Ensemble Perturbation Data at the Start  ---
    ! ensPerts(1          ,:) contains the full perturbations
    ! ensPerts(2:nWaveBand,:) already allocated but empty
    !
    ! --- Ensemble Perturbation Data at the End    ---
    ! ensPerts(nWaveBand,:) contains the largest scales
    ! ...
    ! ensPerts(1        ,:) contains the smallest scales
    !

    !
    ! --> For SDLwSL
    !
    ! --- Ensemble Perturbation Data at the Start  ---
    ! ensPerts(1,:) contains the full perturbations
    !
    ! --- Ensemble Perturbation Data at the End    ---
    ! ensPerts(1,:) contains the selected scales
    !

    if ( mmpi_myid == 0 ) then
      write(*,*)
      write(*,*) 'Scale decomposition of the ensemble perturbations'
      write(*,*) '   number of WaveBands for filtering = ', bEns(instanceIndex)%nWaveBandForFiltering
      write(*,*) '   WaveBand Peaks (total wavenumber)...'
      do waveBandIndex = 1, bEns(instanceIndex)%nWaveBandForFiltering
        write(*,*) waveBandIndex, bEns(instanceIndex)%waveBandPeaks(waveBandIndex)
      end do
    end if

    !
    !- Setup a spectral transform for filtering (nk = nEnsOverDimension)
    !
    if (trim(bEns(instanceIndex)%localizationType) == 'ScaleDependent') then
      nTrunc = bEns(instanceIndex)%waveBandPeaks(1)
    else ! ScaleDependentWithSpectralLoc
      if (bEns(instanceIndex)%waveBandIndexSelected == 1) then
        nTrunc = -1 ! no truncation needed to extract the smallest scales
      else
        nTrunc = bEns(instanceIndex)%waveBandPeaks(bEns(instanceIndex)%waveBandIndexSelected-1)
      end if
    end if

    if (bEns(instanceIndex)%hco_ens%global) then
      ! Global mode
      gstFilterID = gst_setup(bEns(instanceIndex)%ni,bEns(instanceIndex)%nj,nTrunc,bEns(instanceIndex)%nEnsOverDimension)
      if (mmpi_myid == 0) write(*,*) 'ben : returned value of gstFilterID = ',gstFilterID

      nla_filter = gst_getNla(gstFilterID)
      nphase_filter = 2

      allocate(nIndex_vec(nla_filter))
      call mmpi_setup_m(nTrunc,mymBeg,mymEnd,mymSkip,mymCount)
      call mmpi_setup_n(nTrunc,mynBeg,mynEnd,mynSkip,mynCount)
      ila_filter = 0
      do mIndex = mymBeg, mymEnd, mymSkip
        do nIndex = mynBeg, mynEnd, mynSkip
          if (mIndex.le.nIndex) then
            ila_filter = ila_filter + 1
            nIndex_vec(ila_filter) = nIndex
          end if
        end do
      end do

    else
      ! LAM mode
      call lst_Setup(lst_ben_filter,                                                                           & ! OUT
                     bEns(instanceIndex)%ni, bEns(instanceIndex)%nj, bEns(instanceIndex)%hco_ens%dlon, nTrunc, & ! IN
                     'LatLonMN', maxlevels_opt=bEns(instanceIndex)%nEnsOverDimension, gridDataOrder_opt='kij' )  ! IN

      nla_filter    = lst_ben_filter%nla
      nphase_filter = lst_ben_filter%nphase
    end if

    !
    !- 1.  Scale decomposition for every wave band except for wave band #1
    !
    if (trim(bEns(instanceIndex)%localizationType) == 'ScaleDependent') then
      waveBandIndexStart     = 2                ! Skip the smallest scales
      waveBandIndexEnd       = bEns(instanceIndex)%nWaveBand
      waveBandIndexLoopStart = waveBandIndexEnd ! Start with the largest scales
      waveBandIndexLoopEnd   = waveBandIndexStart
      waveBandIndexLoopDirection = -1 
    else ! ScaleDependentWithSpectralLoc
      waveBandIndexStart     = bEns(instanceIndex)%waveBandIndexSelected
      waveBandIndexEnd       = bEns(instanceIndex)%waveBandIndexSelected
      waveBandIndexLoopStart = waveBandIndexStart
      waveBandIndexLoopEnd   = waveBandIndexEnd 
      waveBandIndexLoopDirection = 1 
    end if

    allocate(ResponseFunction(nla_filter,waveBandIndexStart:waveBandIndexEnd))
    allocate(ensPertSP(nla_filter,nphase_filter,bEns(instanceIndex)%nEnsOverDimension))
    allocate(ensPertSPfiltered(nla_filter,nphase_filter,bEns(instanceIndex)%nEnsOverDimension))
    allocate(ensPertGD(bEns(instanceIndex)%nEnsOverDimension,bEns(instanceIndex)%myLonBeg:bEns(instanceIndex)%myLonEnd,bEns(instanceIndex)%myLatBeg:bEns(instanceIndex)%myLatEnd))

    ensPertSP        (:,:,:) = 0.0d0
    ensPertSPfiltered(:,:,:) = 0.0d0

    !- 1.1 Pre-compute the response function
    do waveBandIndex = waveBandIndexLoopStart, waveBandIndexLoopEnd, waveBandIndexLoopDirection
      do ila_filter = 1, nla_filter
        if (bEns(instanceIndex)%hco_ens%global) then
          totwvnb_r8 = real(nIndex_vec(ila_filter),8)
        else
          totwvnb_r8 = lst_ben_filter%k_r8(ila_filter)
        end if
        ResponseFunction(ila_filter,waveBandIndex) = spf_FilterResponseFunction(totwvnb_r8,waveBandIndex, bEns(instanceIndex)%waveBandPeaks, &
                                                                                bEns(instanceIndex)%nWaveBandForFiltering)
        if (trim(bEns(instanceIndex)%localizationType) == 'ScaleDependentWithSpectralLoc') then
          ResponseFunction(ila_filter,waveBandIndex) = sqrt(ResponseFunction(ila_filter,waveBandIndex)) 
        end if
        write(*,*) totwvnb_r8, ResponseFunction(ila_filter,waveBandIndex)
      end do
    end do
    if (bEns(instanceIndex)%hco_ens%global) deallocate(nIndex_vec)

    do stepIndex = 1, bEns(instanceIndex)%numStep ! Loop on ensemble time bin
      do levIndex = 1, ens_getNumK(bEns(instanceIndex)%ensPerts(1)) ! Loop on variables and vertical levels

        ptr4d_r4 => ens_getOneLev_r4(bEns(instanceIndex)%ensPerts(1),levIndex)
          
        !- 1.2 GridPoint space -> Spectral Space
        !$OMP PARALLEL DO PRIVATE (latIndex)
        do latIndex = bEns(instanceIndex)%myLatBeg, bEns(instanceIndex)%myLatEnd
          ensPertGD(:,:,latIndex) = 0.0d0
        end do
        !$OMP END PARALLEL DO
        !$OMP PARALLEL DO PRIVATE (memberIndex,latIndex,lonIndex)
        do latIndex = bEns(instanceIndex)%myLatBeg, bEns(instanceIndex)%myLatEnd
          do lonIndex = bEns(instanceIndex)%myLonBeg, bEns(instanceIndex)%myLonEnd
            do memberIndex = 1, bEns(instanceIndex)%nEns
              ensPertGD(memberIndex,lonIndex,latIndex) = dble(ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex))
            end do
          end do
        end do
        !$OMP END PARALLEL DO
        if (bEns(instanceIndex)%hco_ens%global) then
          ! Global Mode
          call gst_setID(gstFilterID) ! IN
          call gst_reespe_kij(ensPertSP, & ! OUT
                              ensPertGD)   ! IN
        else
          ! LAM mode
          kind = 'GridPointToSpectral'
          call lst_VarTransform(lst_ben_filter,         & ! IN
                                ensPertSP,              & ! OUT
                                ensPertGD,              & ! IN 
                                kind, bEns(instanceIndex)%nEnsOverDimension ) ! IN
        end if

        !- 1.3 Filtering and transformation back to grid point space 
        do waveBandIndex = waveBandIndexLoopStart, waveBandIndexLoopEnd, waveBandIndexLoopDirection

          ! Filtering
          !$OMP PARALLEL DO PRIVATE (memberIndex,p,ila_filter)
          do memberIndex = 1, bEns(instanceIndex)%nEns
            do p = 1, nphase_filter
              do ila_filter = 1, nla_filter
                ensPertSPfiltered(ila_filter,p,memberIndex) = &
                     ensPertSP(ila_filter,p,memberIndex) * ResponseFunction(ila_filter,waveBandIndex)
              end do
            end do
          end do
          !$OMP END PARALLEL DO

          ! Spectral Space -> GridPoint space
          if (bEns(instanceIndex)%hco_ens%global) then
            ! Global Mode
            call gst_setID(gstFilterID) ! IN
            call gst_speree_kij(ensPertSPfiltered, & ! IN
                                ensPertGD)           ! OUT
          else
            ! LAM mode
            kind = 'SpectralToGridPoint'
            call lst_VarTransform(lst_ben_filter,         & ! IN
                                  ensPertSPfiltered,      & ! IN
                                  ensPertGD,              & ! OUT
                                  kind, bEns(instanceIndex)%nEnsOverDimension ) ! IN
          end if

          if (trim(bEns(instanceIndex)%localizationType) == 'ScaleDependent') then
            ptr4d_r4 => ens_getOneLev_r4(bEns(instanceIndex)%ensPerts(waveBandIndex),levIndex)
          else
            ptr4d_r4 => ens_getOneLev_r4(bEns(instanceIndex)%ensPerts(1),levIndex)
          end if
          !$OMP PARALLEL DO PRIVATE (memberIndex,latIndex,lonIndex)
          do latIndex = bEns(instanceIndex)%myLatBeg, bEns(instanceIndex)%myLatEnd
            do lonIndex = bEns(instanceIndex)%myLonBeg, bEns(instanceIndex)%myLonEnd
              do memberIndex = 1, bEns(instanceIndex)%nEns
                ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex) = sngl(ensPertGD(memberIndex,lonIndex,latIndex))
              end do
            end do
          end do
          !$OMP END PARALLEL DO

        end do ! waveBandIndex
      end do ! time bins
    end do ! variables&levels

    deallocate(ensPertGD)
    deallocate(ResponseFunction)
    deallocate(ensPertSP)
    deallocate(ensPertSPfiltered)

    !
    !- 2.  For SDL only, Isolate the smallest scales in waveBandIndex = 1 by difference in grid point space
    !
    if (trim(bEns(instanceIndex)%localizationType) == 'ScaleDependent') then

      allocate(bandSum(bEns(instanceIndex)%myLonBeg:bEns(instanceIndex)%myLonEnd,bEns(instanceIndex)%myLatBeg:bEns(instanceIndex)%myLatEnd))
      do stepIndex = 1, bEns(instanceIndex)%numStep
        !$OMP PARALLEL DO PRIVATE (memberIndex,levIndex,latIndex,lonIndex,waveBandIndex,bandsum,ptr4d_r4)
        do levIndex = 1, ens_getNumK(bEns(instanceIndex)%ensPerts(1))
          do memberIndex = 1, bEns(instanceIndex)%nEns
            bandSum(:,:) = 0.d0
            do waveBandIndex = 2, bEns(instanceIndex)%nWaveBand
              ptr4d_r4 => ens_getOneLev_r4(bEns(instanceIndex)%ensPerts(waveBandIndex),levIndex)
              do latIndex = bEns(instanceIndex)%myLatBeg, bEns(instanceIndex)%myLatEnd
                do lonIndex = bEns(instanceIndex)%myLonBeg, bEns(instanceIndex)%myLonEnd
                  bandSum(lonIndex,latIndex) = bandSum(lonIndex,latIndex) + dble(ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex))
                end do
              end do
            end do
            ptr4d_r4 => ens_getOneLev_r4(bEns(instanceIndex)%ensPerts(1),levIndex)
            do latIndex = bEns(instanceIndex)%myLatBeg, bEns(instanceIndex)%myLatEnd
              do lonIndex = bEns(instanceIndex)%myLonBeg, bEns(instanceIndex)%myLonEnd
                ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex) = sngl(dble(ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex)) - bandSum(lonIndex,latIndex))
              end do
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do
      deallocate(bandSum)

    end if

  end subroutine ensembleScaleDecomposition

  !--------------------------------------------------------------------------
  ! ben_reduceToMPILocal
  !--------------------------------------------------------------------------
  subroutine ben_reduceToMPILocal(cv_mpilocal,cv_mpiglobal,instanceIndex)
    implicit none
    integer, intent(in)  :: instanceIndex
    real(8), intent(out) :: cv_mpilocal(bEns(instanceIndex)%cvDim_mpilocal)
    real(8), intent(in)  :: cv_mpiglobal(:)

    if (verbose) write(*,*) 'Entering ben_reduceToMPILocal'

    call loc_reduceToMPILocal(bEns(instanceIndex)%locStorage(1),cv_mpilocal,cv_mpiglobal)
    
  end subroutine ben_reduceToMPILocal

  !--------------------------------------------------------------------------
  ! ben_reduceToMPILocal_r4
  !--------------------------------------------------------------------------
  subroutine ben_reduceToMPILocal_r4(cv_mpilocal,cv_mpiglobal,instanceIndex)
    implicit none
    integer, intent(in)  :: instanceIndex
    real(4), intent(out) :: cv_mpilocal(bEns(instanceIndex)%cvDim_mpilocal)
    real(4), intent(in)  :: cv_mpiglobal(:)

    if (verbose) write(*,*) 'Entering reduceToMPILocal_r4'

    call loc_reduceToMPILocal_r4(bEns(instanceIndex)%locStorage(1),cv_mpilocal,cv_mpiglobal) ! IN

  end subroutine ben_reduceToMPILocal_r4
  
  !--------------------------------------------------------------------------
  ! ben_expandToMPIGlobal
  !--------------------------------------------------------------------------
  subroutine ben_expandToMPIGlobal(cv_mpilocal,cv_mpiglobal,instanceIndex)
    implicit none
    integer, intent(in)  :: instanceIndex
    real(8), intent(in)  :: cv_mpilocal(bEns(instanceIndex)%cvDim_mpilocal)
    real(8), intent(out) :: cv_mpiglobal(:)

    if (verbose) write(*,*) 'Entering ben_expandToMPIGlobal'

    call loc_expandToMPIGlobal(bEns(instanceIndex)%locStorage(1), cv_mpilocal,  & ! IN
                               cv_mpiglobal)                                      ! OUT  

  end subroutine ben_expandToMPIGlobal

  !--------------------------------------------------------------------------
  ! ben_expandToMPIGlobal_r4
  !--------------------------------------------------------------------------
  subroutine ben_expandToMPIGlobal_r4(cv_mpilocal,cv_mpiglobal,instanceIndex)
    implicit none
    integer, intent(in)  :: instanceIndex
    real(4), intent(in)  :: cv_mpilocal(bEns(instanceIndex)%cvDim_mpilocal)
    real(4), intent(out) :: cv_mpiglobal(:)

    if (verbose) write(*,*) 'Entering ben_expandToMPIGlobal_r4'

    call loc_expandToMPIGlobal_r4(bEns(instanceIndex)%locStorage(1), cv_mpilocal,  & ! IN
                                  cv_mpiglobal)                                      ! OUT

  end subroutine ben_expandToMPIGlobal_r4

  !--------------------------------------------------------------------------
  ! ben_BSqrt
  !--------------------------------------------------------------------------
  subroutine ben_BSqrt(instanceIndex, controlVector_in, statevector,  &
                       useFSOFcst_opt, stateVectorRef_opt)
    implicit none

    ! Arguments
    integer, intent(in) :: instanceIndex
    real(8)          :: controlVector_in(bEns(instanceIndex)%cvDim_mpilocal) 
    type(struct_gsv) :: statevector
    logical,optional :: useFSOFcst_opt
    type(struct_gsv), optional :: statevectorRef_opt

    ! Locals
    type(struct_ens), target  :: ensAmplitude
    type(struct_ens), pointer :: ensAmplitude_ptr

    integer   :: ierr, waveBandIndex
    integer   :: numStepAmplitude, amp3dStepIndex
    logical   :: immediateReturn
    logical   :: useFSOFcst

    if (mmpi_doBarrier) call rpn_comm_barrier('GRID',ierr)

    !
    !- 1.  Tests
    !
    if (.not. bEns(instanceIndex)%initialized) then
      if (mmpi_myid == 0) write(*,*) 'ben_bsqrt: bMatrixEnsemble not initialized'
      return
    end if

    if (sum(bEns(instanceIndex)%scaleFactor) == 0.0d0) then
      if (mmpi_myid == 0) write(*,*) 'ben_bsqrt: scaleFactor=0, skipping bSqrt'
      return
    end if

    ! only check controlVector on proc 0, since may be zero-length on some procs
    if (mmpi_myid == 0) then
      immediateReturn = .false.
      if (maxval(controlVector_in) == 0.0d0 .and. minval(controlVector_in) == 0.0d0) then
        write(*,*) 'ben_bsqrt: controlVector=0, skipping bSqrt'
        immediateReturn = .true.
      end if
    end if
    call rpn_comm_bcast(immediateReturn, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)
    if (immediateReturn) return

    if (mmpi_myid == 0) write(*,*) 'ben_bsqrt: starting'
    if (mmpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if (present(useFSOFcst_opt)) then
      useFSOFcst = useFSOFcst_opt
    else
      useFSOFcst = .false.
    end if

    !
    !- 2.  Compute the analysis increment from Bens
    !
    if (verbose) write(*,*) 'ben_bsqrt: allocating ensAmplitude'
    if (useFSOFcst) then
      numStepAmplitude = bEns(instanceIndex)%numStepAmplitudeFSOFcst
      amp3dStepIndex   = bEns(instanceIndex)%amp3dStepIndexFSOFcst
    else
      numStepAmplitude = bEns(instanceIndex)%numStepAmplitudeAssimWindow
      amp3dStepIndex   = bEns(instanceIndex)%amp3dStepIndexAssimWindow
    end if
    if (bEns(instanceIndex)%useSaveAmp) then
      ensAmplitude_ptr => ensAmplitudeSave(instanceIndex)
    else
      call ens_allocate(ensAmplitude, bEns(instanceIndex)%nEnsOverDimension, numStepAmplitude,                     &
                        bEns(instanceIndex)%hco_ens,   &
                        bEns(instanceIndex)%vco_ens, bEns(instanceIndex)%dateStampList, &
                        hco_core_opt=bEns(instanceIndex)%hco_core, &
                        varNames_opt=bEns(instanceIndex)%varNameALFA, dataKind_opt=8)
      ensAmplitude_ptr => ensAmplitude
    end if
    call gsv_zero(statevector)

    do waveBandIndex = 1, bEns(instanceIndex)%nWaveBand !  Loop on WaveBand (for ScaleDependent Localization)

      ! 2.1 Compute the ensemble amplitudes
      call utl_tmg_start(60,'------LocSpectral_TL')
      call loc_Lsqrt(bEns(instanceIndex)%locStorage(waveBandIndex),controlVector_in, & ! IN
                     ensAmplitude_ptr,                                               & ! OUT
                     amp3dStepIndex)                                                   ! IN
      call utl_tmg_stop(60)

      ! 2.2 Advect the amplitudes
      if      (bEns(instanceIndex)%advectAmplitudeFSOFcst   .and. useFSOFcst) then
        call adv_ensemble_tl( ensAmplitude_ptr,                                                  & ! INOUT
                              bEns(instanceIndex)%adv_amplitudeFSOFcst, bEns(instanceIndex)%nEns )   ! IN
      else if (bEns(instanceIndex)%advectAmplitudeAssimWindow .and. .not. useFSOFcst) then
        call adv_ensemble_tl( ensAmplitude_ptr,                                                    & ! INOUT
                              bEns(instanceIndex)%adv_amplitudeAssimWindow, bEns(instanceIndex)%nEns ) ! IN
      end if

      if ( bEns(instanceIndex)%keepAmplitude .and. waveBandIndex == 1 ) then
        call ens_copy(ensAmplitude_ptr, bEns(instanceIndex)%ensAmplitudeStorage)
      end if

      ! 2.3 Compute increment by multiplying amplitudes by member perturbations
      call addEnsMember( ensAmplitude_ptr, statevector,          & ! INOUT 
                         instanceIndex, waveBandIndex, useFSOFcst )  ! IN

    end do ! Loop on WaveBand

    if (.not. bEns(instanceIndex)%useSaveAmp) call ens_deallocate(ensAmplitude)

    ! 2.4 Apply the Std. Dev (if needed)
    if (bEns(instanceIndex)%ensPertsNormalized .and. .not. bEns(instanceIndex)%useCmatrixOnly) then
      call gsv_schurProduct(bEns(instanceIndex)%statevector_ensStdDev, & ! IN
                            statevector)                                 ! INOUT 
    end if

    ! 2.5 Advect Increments
    if ( bEns(instanceIndex)%advectEnsPertAnlInc ) then
      call adv_statevector_tl( statevector,                      & ! INOUT
                               bEns(instanceIndex)%adv_analInc )   ! IN
    end if

    !
    !- 3.  Variable transforms
    !
    if ( bEns(instanceIndex)%ctrlVarHumidity == 'LQ' .and. gsv_varExist(varName='HU') ) then
       call gvt_transform( statevector,   &                        ! INOUT
                           'LQtoHU_tlm',  &                        ! IN
                           stateVectorRef_opt=stateVectorRef_opt ) ! IN
    end if

    if (mmpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if (mmpi_myid == 0) write(*,*) 'ben_bsqrt: done'

  end subroutine ben_BSqrt

  !--------------------------------------------------------------------------
  ! ben_BSqrtAd
  !--------------------------------------------------------------------------
  subroutine ben_BSqrtAd(instanceIndex, statevector, controlVector_out,  &
                         useFSOFcst_opt, stateVectorRef_opt)
    implicit none

    ! Arguments
    integer, intent(in)        :: instanceIndex
    real(8)                    :: controlVector_out(bEns(instanceIndex)%cvDim_mpilocal) 
    type(struct_gsv)           :: statevector
    logical, optional          :: useFSOFcst_opt
    type(struct_gsv), optional :: statevectorRef_opt

    ! Locals
    type(struct_ens), target  :: ensAmplitude
    type(struct_ens), pointer :: ensAmplitude_ptr

    integer           :: ierr, waveBandIndex
    integer           :: numStepAmplitude,amp3dStepIndex
    logical           :: useFSOFcst

    !
    !- 1.  Tests
    !
    if (mmpi_doBarrier) call rpn_comm_barrier('GRID',ierr)

    if (.not. bEns(instanceIndex)%initialized) then
      if (mmpi_myid == 0) write(*,*) 'ben_bsqrtad: bMatrixEnsemble not initialized'
      return
    end if

    if (sum(bEns(instanceIndex)%scaleFactor) == 0.0d0) then
      if (mmpi_myid == 0) write(*,*) 'ben_bsqrtad: scaleFactor=0, skipping bSqrtAd'
      return
    end if

    if (mmpi_myid == 0) write(*,*) 'ben_bsqrtad: starting'
    if (mmpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if (present(useFSOFcst_opt)) then
      useFSOFcst = useFSOFcst_opt
    else
      useFSOFcst = .false.
    end if

    !
    !- 3.  Variable transforms
    !
    if ( bEns(instanceIndex)%ctrlVarHumidity == 'LQ' .and. gsv_varExist(varName='HU') ) then
      call gvt_transform( statevector,  &                         ! INOUT
                          'LQtoHU_ad',  &                         ! IN
                          stateVectorRef_opt=stateVectorRef_opt ) ! IN
    end if

    !
    !- 2.  Compute the analysis increment from Bens
    !

    ! 2.5 Advect Increments
    if ( bEns(instanceIndex)%advectEnsPertAnlInc ) then
      call adv_statevector_ad( statevector,  & ! INOUT
                               bEns(instanceIndex)%adv_analInc )   ! IN
    end if

    ! 2.4 Apply the Std. Dev (if needed)
    if (bEns(instanceIndex)%ensPertsNormalized .and. .not. bEns(instanceIndex)%useCmatrixOnly) then
      call gsv_schurProduct(bEns(instanceIndex)%statevector_ensStdDev, & ! IN
                            statevector)                                 ! INOUT 
    end if

    if (verbose) write(*,*) 'ben_bsqrtAd: allocating ensAmplitude'
    if (useFSOFcst) then
      numStepAmplitude = bEns(instanceIndex)%numStepAmplitudeFSOFcst
      amp3dStepIndex   = bEns(instanceIndex)%amp3dStepIndexFSOFcst
    else
      numStepAmplitude = bEns(instanceIndex)%numStepAmplitudeAssimWindow
      amp3dStepIndex   = bEns(instanceIndex)%amp3dStepIndexAssimWindow
    end if
    if (bEns(instanceIndex)%useSaveAmp) then
      ensAmplitude_ptr => ensAmplitudeSave(instanceIndex)
    else
      call ens_allocate(ensAmplitude, bEns(instanceIndex)%nEnsOverDimension, numStepAmplitude,                     &
                        bEns(instanceIndex)%hco_ens,  &
                        bEns(instanceIndex)%vco_ens, bEns(instanceIndex)%dateStampList, &
                        hco_core_opt=bEns(instanceIndex)%hco_core,  &
                        varNames_opt=bEns(instanceIndex)%varNameALFA, dataKind_opt=8)
      ensAmplitude_ptr => ensAmplitude
    end if

    do waveBandIndex = 1, bEns(instanceIndex)%nWaveBand !  Loop on WaveBand (for ScaleDependent Localization)

      ! 2.3 Compute increment by multiplying amplitudes by member perturbations
      call addEnsMemberAd( statevector, ensAmplitude_ptr,         & ! INOUT
                           instanceIndex, waveBandIndex, useFSOFcst)  ! IN

      ! 2.2 Advect the  amplitudes
      if      (bEns(instanceIndex)%advectAmplitudeFSOFcst .and. useFSOFcst) then
        call adv_ensemble_ad( ensAmplitude_ptr,          & ! INOUT
                              bEns(instanceIndex)%adv_amplitudeFSOFcst, bEns(instanceIndex)%nEns )   ! IN
      else if (bEns(instanceIndex)%advectAmplitudeAssimWindow .and. .not. useFSOFcst) then
        call adv_ensemble_ad( ensAmplitude_ptr,                                                    & ! INOUT
                              bEns(instanceIndex)%adv_amplitudeAssimWindow, bEns(instanceIndex)%nEns ) ! IN
      end if

      ! 2.1 Compute the ensemble amplitudes
      call utl_tmg_start(64,'------LocSpectral_AD')
      call loc_LsqrtAd(bEns(instanceIndex)%locStorage(waveBandIndex),ensAmplitude_ptr, & ! IN
                       controlVector_out,                                              & ! OUT
                       amp3dStepIndex)                                                   ! IN
      call utl_tmg_stop(64)

    end do ! Loop on WaveBand

    if (.not. bEns(instanceIndex)%useSaveAmp) call ens_deallocate(ensAmplitude)

    if (mmpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if (mmpi_myid == 0) write(*,*) 'ben_bsqrtAd: done'

  end subroutine ben_BSqrtAd

  !--------------------------------------------------------------------------
  ! addEnsMember
  !--------------------------------------------------------------------------
  subroutine addEnsMember(ensAmplitude, statevector_out, &
                          instanceIndex, waveBandIndex, useFSOFcst_opt)
    implicit none

    type(struct_ens)    :: ensAmplitude
    type(struct_gsv)    :: statevector_out
    integer, intent(in) :: instanceIndex
    integer, intent(in) :: waveBandIndex
    logical, optional   :: useFSOFcst_opt

    real(8), pointer    :: ensAmplitude_oneLev(:,:,:,:)
    real(8), pointer    :: ensAmplitude_oneLevM1(:,:,:,:), ensAmplitude_oneLevP1(:,:,:,:)
    real(8), allocatable, target :: ensAmplitude_MT(:,:,:,:)
    real(8), pointer     :: ensAmplitude_MT_ptr(:,:,:,:)
    real(4), pointer     :: increment_out_r4(:,:,:,:)
    real(8), pointer     :: increment_out_r8(:,:,:,:)
    real(8), allocatable :: increment_out2(:,:,:)
    real(4), pointer     :: ensMemberAll_r4(:,:,:,:)
    integer     :: lev, lev2, levIndex, stepIndex, stepIndex_amp, latIndex, lonIndex, topLevOffset, memberIndex
    character(len=4)     :: varName

    logical             :: useFSOFcst
    integer             :: stepIndex2, stepBeg, stepEnd, numStepAmplitude

    if (verbose) write(*,*) 'Entering ben_addEnsMember'

    call utl_tmg_start(58,'------AddMem_TL')

    if (present(useFSOFcst_opt)) then
      useFSOFcst = useFSOFcst_opt
    else
      useFSOFcst = .false.
    end if
    if (useFSOFcst .and. bEns(instanceIndex)%fsoLeadTime > 0.0d0) then
      stepBeg = bEns(instanceIndex)%numStep
      stepEnd = stepBeg
      if (mmpi_myid == 0) write(*,*) 'ben_bsqrtad: using forecast ensemble stored at timestep ',stepEnd
    else
      stepBeg = 1
      stepEnd = bEns(instanceIndex)%numStepAssimWindow
    end if
    if (useFSOFcst) then
      numStepAmplitude =  bEns(instanceIndex)%numStepAmplitudeFSOFcst
    else
      numStepAmplitude =  bEns(instanceIndex)%numStepAmplitudeAssimWindow
    end if
    
    allocate(ensAmplitude_MT(bEns(instanceIndex)%nEns,numStepAmplitude,bEns(instanceIndex)%myLonBeg:bEns(instanceIndex)%myLonEnd,bEns(instanceIndex)%myLatBeg:bEns(instanceIndex)%myLatEnd))
    allocate(increment_out2(bEns(instanceIndex)%numStep,bEns(instanceIndex)%myLonBeg:bEns(instanceIndex)%myLonEnd,bEns(instanceIndex)%myLatBeg:bEns(instanceIndex)%myLatEnd))

    do levIndex = 1, ens_getNumK(bEns(instanceIndex)%ensPerts(waveBandIndex))

      lev = ens_getLevFromK(bEns(instanceIndex)%ensPerts(1),levIndex)
      varName = ens_getVarNameFromK(bEns(instanceIndex)%ensPerts(1),levIndex)

      !$OMP PARALLEL DO PRIVATE (latIndex)
      do latIndex = bEns(instanceIndex)%myLatBeg, bEns(instanceIndex)%myLatEnd
        increment_out2(:,:,latIndex) = 0.0d0
      end do
      !$OMP END PARALLEL DO

      if ( vnl_varLevelFromVarname(varName) /= 'SF' .and. &
           vnl_varLevelFromVarname(varName) == vnl_varLevelFromVarname(bEns(instanceIndex)%varNameALFA(1)) ) then
        ! The non-surface variable varName is on the same levels than the amplitude field
        ensAmplitude_oneLev => ens_getOneLev_r8(ensAmplitude,lev)
        ensAmplitude_MT_ptr(1:,1:,bEns(instanceIndex)%myLonBeg:,bEns(instanceIndex)%myLatBeg:) => ensAmplitude_oneLev(1:bEns(instanceIndex)%nEns,:,:,:)

      else if (vnl_varLevelFromVarname(varName) == 'TH') then
        ! The non-surface variable varName is on TH levels whereas the amplitude field is on MM levels

        if (bEns(instanceIndex)%vco_anl%Vcode == 5002) then

          if (lev == 1) then
            ! use top momentum level amplitudes for top thermo level
            ensAmplitude_oneLev => ens_getOneLev_r8(ensAmplitude,lev)
            ensAmplitude_MT_ptr(1:,1:,bEns(instanceIndex)%myLonBeg:,bEns(instanceIndex)%myLatBeg:) => ensAmplitude_oneLev(1:bEns(instanceIndex)%nEns,:,:,:)
          else if (lev == bEns(instanceIndex)%nLevEns_T) then
            ! use surface momentum level amplitudes for surface thermo level
            ensAmplitude_oneLev => ens_getOneLev_r8(ensAmplitude,bEns(instanceIndex)%nLevEns_M)
            ensAmplitude_MT_ptr(1:,1:,bEns(instanceIndex)%myLonBeg:,bEns(instanceIndex)%myLatBeg:) => ensAmplitude_oneLev(1:bEns(instanceIndex)%nEns,:,:,:)
          else
            ! for other levels, interpolate momentum weights to get thermo amplitudes
            ensAmplitude_oneLev   => ens_getOneLev_r8(ensAmplitude,lev)
            ensAmplitude_oneLevM1 => ens_getOneLev_r8(ensAmplitude,lev-1)
            !$OMP PARALLEL DO PRIVATE (latIndex)
            do latIndex = bEns(instanceIndex)%myLatBeg, bEns(instanceIndex)%myLatEnd
              ensAmplitude_MT(:,:,:,latIndex) = 0.5d0*( ensAmplitude_oneLevM1(1:bEns(instanceIndex)%nEns,:,:,latIndex) +   &
                                                        ensAmplitude_oneLev(1:bEns(instanceIndex)%nEns,:,:,latIndex) )
            end do
            !$OMP END PARALLEL DO
            ensAmplitude_MT_ptr(1:,1:,bEns(instanceIndex)%myLonBeg:,bEns(instanceIndex)%myLatBeg:) => ensAmplitude_MT(:,:,:,:)
          end if

        else if (bEns(instanceIndex)%vco_anl%Vcode == 5005) then

          if (lev == bEns(instanceIndex)%nLevEns_T) then
            ! use surface momentum level amplitudes for surface thermo level
            ensAmplitude_oneLev => ens_getOneLev_r8(ensAmplitude,bEns(instanceIndex)%nLevEns_M)
            ensAmplitude_MT_ptr(1:,1:,bEns(instanceIndex)%myLonBeg:,bEns(instanceIndex)%myLatBeg:) => ensAmplitude_oneLev(1:bEns(instanceIndex)%nEns,:,:,:)
          else
            ! for other levels, interpolate momentum weights to get thermo amplitudes
            ensAmplitude_oneLev   => ens_getOneLev_r8(ensAmplitude,lev)
            ensAmplitude_oneLevP1 => ens_getOneLev_r8(ensAmplitude,lev+1)
            !$OMP PARALLEL DO PRIVATE (latIndex)
            do latIndex = bEns(instanceIndex)%myLatBeg, bEns(instanceIndex)%myLatEnd
              ensAmplitude_MT(:,:,:,latIndex) = 0.5d0*( ensAmplitude_oneLevP1(1:bEns(instanceIndex)%nEns,:,:,latIndex) +   &
                                                        ensAmplitude_oneLev(1:bEns(instanceIndex)%nEns,:,:,latIndex) )
            end do
            !$OMP END PARALLEL DO
            ensAmplitude_MT_ptr(1:,1:,bEns(instanceIndex)%myLonBeg:,bEns(instanceIndex)%myLatBeg:) => ensAmplitude_MT(:,:,:,:)
          end if

        else
          call utl_abort('ben_addEnsMember: incompatible vcode')
        end if

      else if (vnl_varLevelFromVarname(varName) == 'SF') then

        ! Surface variable cases (atmosphere or surface only)
        if (bEns(instanceIndex)%vco_anl%Vcode == 5002 .or. bEns(instanceIndex)%vco_anl%Vcode == 5005) then
          ensAmplitude_oneLev   => ens_getOneLev_r8(ensAmplitude,bEns(instanceIndex)%nLevEns_M)
        else ! vco_anl%Vcode == 0
          ensAmplitude_oneLev   => ens_getOneLev_r8(ensAmplitude,1)
        end if
        ensAmplitude_MT_ptr(1:,1:,bEns(instanceIndex)%myLonBeg:,bEns(instanceIndex)%myLatBeg:) => ensAmplitude_oneLev(1:bEns(instanceIndex)%nEns,:,:,:)

      else if (vnl_varLevelFromVarname(varName) == 'SFDP') then

        ! Surface variable cases (ocean)
        ensAmplitude_oneLev   => ens_getOneLev_r8(ensAmplitude,1)

      else

        write(*,*) 'variable name = ', varName
        write(*,*) 'varLevel      = ', vnl_varLevelFromVarname(varName)
        call utl_abort('ben_addEnsMember: unknown value of varLevel')

      end if

      call utl_tmg_start(59,'--------AddMemInner_TL')

      ensMemberAll_r4 => ens_getOneLev_r4(bEns(instanceIndex)%ensPerts(waveBandIndex),levIndex)
      !$OMP PARALLEL DO PRIVATE (latIndex,lonIndex,stepIndex,stepIndex2,stepIndex_amp,memberIndex)
      do latIndex = bEns(instanceIndex)%myLatBeg, bEns(instanceIndex)%myLatEnd
        do lonIndex = bEns(instanceIndex)%myLonBeg, bEns(instanceIndex)%myLonEnd
          do stepIndex = StepBeg, StepEnd
            stepIndex2 = stepIndex - stepBeg + 1
            if      (bEns(instanceIndex)%advectAmplitudeFSOFcst   .and. useFSOFcst) then
              stepIndex_amp = 2
            else if (bEns(instanceIndex)%advectAmplitudeAssimWindow .and. .not. useFSOFcst) then
              stepIndex_amp = stepIndex
            else
              stepIndex_amp = 1
            end if
            do memberIndex = 1, bEns(instanceIndex)%nEns
              increment_out2(stepIndex2,lonIndex,latIndex) = increment_out2(stepIndex2,lonIndex,latIndex) +   &
                ensAmplitude_MT_ptr(memberIndex,stepIndex_amp,lonIndex,latIndex) *  &
                dble(ensMemberAll_r4(memberIndex,stepIndex,lonIndex,latIndex))
            end do ! memberIndex
          end do ! stepIndex
        end do ! lonIndex
      end do ! latIndex
      !$OMP END PARALLEL DO

      call utl_tmg_stop(59)

      ! compute increment level from amplitude/member level
      if (vnl_varLevelFromVarname(varName) == 'SF') then
        topLevOffset = 1
      else if (vnl_varLevelFromVarname(varName) == 'MM') then
        topLevOffset = bEns(instanceIndex)%topLevIndex_M
      else if (vnl_varLevelFromVarname(varName) == 'TH') then
        topLevOffset = bEns(instanceIndex)%topLevIndex_T
      else if (vnl_varLevelFromVarname(varName) == 'SFDP') then
        topLevOffset = 1
      else if (vnl_varLevelFromVarname(varName) == 'DP') then
        topLevOffset = 1
      else
        call utl_abort('ben_addEnsMember: unknown value of varLevel')
      end if
      lev2 = lev - 1 + topLevOffset

      if (varName == 'LQ' .and. bEns(instanceIndex)%ensShouldNotContainLQvarName) then
        if (gsv_getDataKind(statevector_out) == 4) then
          call gsv_getField(statevector_out, increment_out_r4, 'HU')
        else
          call gsv_getField(statevector_out, increment_out_r8, 'HU')
        end if
      else
        if (gsv_getDataKind(statevector_out) == 4) then
          call gsv_getField(statevector_out, increment_out_r4, varName)
        else
          call gsv_getField(statevector_out, increment_out_r8, varName)
        end if
      end if
      !$OMP PARALLEL DO PRIVATE (stepIndex, stepIndex2)
      do stepIndex = StepBeg, StepEnd
        stepIndex2 = stepIndex - StepBeg + 1
        if (gsv_getDataKind(statevector_out) == 4) then
          increment_out_r4(:,:,lev2,stepIndex2) = increment_out_r4(:,:,lev2,stepIndex2) + increment_out2(stepIndex2,:,:)
        else
          increment_out_r8(:,:,lev2,stepIndex2) = increment_out_r8(:,:,lev2,stepIndex2) + increment_out2(stepIndex2,:,:)
        end if
      end do
      !$OMP END PARALLEL DO

    end do ! levIndex

    deallocate(ensAmplitude_MT)
    deallocate(increment_out2)

    call utl_tmg_stop(58)

  end subroutine addEnsMember

  !--------------------------------------------------------------------------
  ! addEnsMemberAd
  !--------------------------------------------------------------------------
  subroutine addEnsMemberAd(statevector_in, ensAmplitude, &
                            instanceIndex, waveBandIndex, useFSOFcst_opt)
    implicit none

    type(struct_ens)   :: ensAmplitude
    type(struct_gsv)   :: statevector_in
    integer,intent(in) :: instanceIndex
    integer,intent(in) :: waveBandIndex
    logical,optional   :: useFSOFcst_opt

    real(8), pointer    :: ensAmplitude_oneLev(:,:,:,:)
    real(8), pointer    :: ensAmplitude_oneLevM1(:,:,:,:), ensAmplitude_oneLevP1(:,:,:,:)
    real(8), allocatable :: ensAmplitude_MT(:,:)
    real(4), pointer     :: increment_in_r4(:,:,:,:)
    real(8), pointer     :: increment_in_r8(:,:,:,:)
    real(8), allocatable :: increment_in2(:,:,:)
    real(4), pointer :: ensMemberAll_r4(:,:,:,:)
    integer          :: levIndex, lev, lev2, stepIndex, stepIndex_amp, latIndex, lonIndex, topLevOffset, memberIndex
    character(len=4) :: varName
    character(len=4) :: varLevel, varLevelAlfa
    integer          :: stepBeg, stepEnd, stepIndex2, numStepAmplitude
    logical          :: useFSOFcst

    if (verbose) write(*,*) 'Entering ben_addEnsMemberAd'

    call utl_tmg_start(62,'------AddMem_AD')

    if (present(useFSOFcst_opt)) then
      useFSOFcst = useFSOFcst_opt
    else
      useFSOFcst = .false.
    end if

    if (useFSOFcst .and. bEns(instanceIndex)%fsoLeadTime > 0.0d0) then
      stepBeg = bEns(instanceIndex)%numStep
      stepEnd = stepBeg
      if (mmpi_myid == 0) write(*,*) 'ben_addEnsMemberAd: using forecast ensemble stored at timestep ',stepEnd
    else
      stepBeg = 1
      stepEnd = bEns(instanceIndex)%numStepAssimWindow
    end if
    if (useFSOFcst) then
      numStepAmplitude =  bEns(instanceIndex)%numStepAmplitudeFSOFcst
    else
      numStepAmplitude =  bEns(instanceIndex)%numStepAmplitudeAssimWindow
    end if

    allocate(ensAmplitude_MT(bEns(instanceIndex)%nEns,numStepAmplitude))
    allocate(increment_in2(bEns(instanceIndex)%numStep,bEns(instanceIndex)%myLonBeg:bEns(instanceIndex)%myLonEnd,bEns(instanceIndex)%myLatBeg:bEns(instanceIndex)%myLatEnd))

    ! set output ensemble Amplitude to zero
    !$OMP PARALLEL DO PRIVATE (levIndex, ensAmplitude_oneLev)
    do levIndex = 1, ens_getNumLev(ensAmplitude,vnl_varLevelFromVarname(bEns(instanceIndex)%varNameALFA(1)))
      ensAmplitude_oneLev => ens_getOneLev_r8(ensAmplitude,levIndex)
      ensAmplitude_oneLev(:,:,:,:) = 0.0d0
    end do
    !$OMP END PARALLEL DO

    do levIndex = 1, ens_getNumK(bEns(instanceIndex)%ensPerts(waveBandIndex))

      lev = ens_getLevFromK(bEns(instanceIndex)%ensPerts(1),levIndex)
      varName = ens_getVarNameFromK(bEns(instanceIndex)%ensPerts(1),levIndex)
      varLevel     = vnl_varLevelFromVarname(varName)
      varLevelAlfa = vnl_varLevelFromVarname(bEns(instanceIndex)%varNameALFA(1))

      ! compute increment level from amplitude/member level
      if (varLevel == 'SF') then
        topLevOffset = 1
      else if (varLevel == 'MM') then
        topLevOffset = bEns(instanceIndex)%topLevIndex_M
      else if (varLevel == 'TH') then
        topLevOffset = bEns(instanceIndex)%topLevIndex_T
      else
        call utl_abort('ben_addEnsMemberAd: unknown value of varLevel')
      end if
      lev2 = lev - 1 + topLevOffset

      if (varName == 'LQ' .and. bEns(instanceIndex)%ensShouldNotContainLQvarName) then
        if (gsv_getDataKind(statevector_in) == 4) then
          call gsv_getField(statevector_in, increment_in_r4, 'HU')
        else
          call gsv_getField(statevector_in, increment_in_r8, 'HU')
        end if
      else
        if (gsv_getDataKind(statevector_in) == 4) then
          call gsv_getField(statevector_in, increment_in_r4, varName)
        else
          call gsv_getField(statevector_in, increment_in_r8, varName)
        end if
      end if
      !$OMP PARALLEL DO PRIVATE (stepIndex, stepIndex2)
      do stepIndex = stepBeg, stepEnd
        stepIndex2 = stepIndex - stepBeg + 1
        if (gsv_getDataKind(statevector_in) == 4) then
          increment_in2(stepIndex2,:,:) = increment_in_r4(:,:,lev2,stepIndex2)
        else
          increment_in2(stepIndex2,:,:) = increment_in_r8(:,:,lev2,stepIndex2)
        end if
      end do
      !$OMP END PARALLEL DO

      ensMemberAll_r4 => ens_getOneLev_r4(bEns(instanceIndex)%ensPerts(waveBandIndex),levIndex)
      !$OMP PARALLEL DO PRIVATE (latIndex,lonIndex,stepIndex, stepIndex2, stepIndex_amp, &
           memberIndex,ensAmplitude_oneLev, ensAmplitude_oneLevM1, &
           ensAmplitude_oneLevP1, ensAmplitude_MT)
      do latIndex = bEns(instanceIndex)%myLatBeg, bEns(instanceIndex)%myLatEnd
        do lonIndex = bEns(instanceIndex)%myLonBeg, bEns(instanceIndex)%myLonEnd

          if (omp_get_thread_num() == 0) call utl_tmg_start(63,'--------AddMemInner_AD')
          ensAmplitude_MT(:,:) = 0.0d0
          do stepIndex = stepBeg, stepEnd
            stepIndex2 = stepIndex - stepBeg + 1
            if      (bEns(instanceIndex)%advectAmplitudeFSOFcst   .and. useFSOFcst) then
              stepIndex_amp = 2
            else if (bEns(instanceIndex)%advectAmplitudeAssimWindow .and. .not. useFSOFcst) then
              stepIndex_amp = stepIndex
            else
              stepIndex_amp = 1
            end if
            do memberIndex = 1, bEns(instanceIndex)%nEns
              ensAmplitude_MT(memberIndex,stepIndex_amp) = ensAmplitude_MT(memberIndex,stepIndex_amp) +  &
                increment_in2(stepIndex2,lonIndex,latIndex) * dble(ensMemberAll_r4(memberIndex,stepIndex,lonIndex,latIndex))
            end do ! memberIndex
          end do ! stepIndex
          if (omp_get_thread_num() == 0) call utl_tmg_stop(63)

          ! transform thermo/momentum level amplitude sensitivites appropriately

          if ( varLevel /= 'SF' .and. &
               varLevel == varLevelAlfa ) then
            ! The non-surface variable varName is on the same levels than the amplitude field

            ensAmplitude_oneLev => ens_getOneLev_r8(ensAmplitude,lev)
            ensAmplitude_oneLev(1:bEns(instanceIndex)%nEns,:,lonIndex,latIndex) = &
                 ensAmplitude_oneLev(1:bEns(instanceIndex)%nEns,:,lonIndex,latIndex) + ensAmplitude_MT(:,:)

          else if (varLevel == 'TH') then
            ! The non-surface variable varName is on TH levels whereas the amplitude field is on MM levels

            if (bEns(instanceIndex)%vco_anl%Vcode == 5002) then

              if (lev == 1) then
                ! use top momentum level amplitudes for top thermo level
                ensAmplitude_oneLev => ens_getOneLev_r8(ensAmplitude,lev)
                ensAmplitude_oneLev(1:bEns(instanceIndex)%nEns,:,lonIndex,latIndex) = &
                   ensAmplitude_oneLev(1:bEns(instanceIndex)%nEns,:,lonIndex,latIndex) + ensAmplitude_MT(:,:)
              else if (lev == bEns(instanceIndex)%nLevEns_T) then
                ! use surface momentum level amplitudes for surface thermo level
                ensAmplitude_oneLev => ens_getOneLev_r8(ensAmplitude,bEns(instanceIndex)%nLevEns_M)
                ensAmplitude_oneLev(1:bEns(instanceIndex)%nEns,:,lonIndex,latIndex) = &
                     ensAmplitude_oneLev(1:bEns(instanceIndex)%nEns,:,lonIndex,latIndex) + ensAmplitude_MT(:,:)
              else
                ! for other levels, interpolate momentum weights to get thermo amplitudes
                ensAmplitude_oneLev   => ens_getOneLev_r8(ensAmplitude,lev)
                ensAmplitude_oneLevM1 => ens_getOneLev_r8(ensAmplitude,lev-1)
                ensAmplitude_oneLev(1:bEns(instanceIndex)%nEns,:,lonIndex,latIndex)   = &
                     ensAmplitude_oneLev(1:bEns(instanceIndex)%nEns,:,lonIndex,latIndex)   + 0.5d0*ensAmplitude_MT(:,:)
                ensAmplitude_oneLevM1(1:bEns(instanceIndex)%nEns,:,lonIndex,latIndex) = &
                     ensAmplitude_oneLevM1(1:bEns(instanceIndex)%nEns,:,lonIndex,latIndex) + 0.5d0*ensAmplitude_MT(:,:)
              end if

            else if (bEns(instanceIndex)%vco_anl%Vcode == 5005) then

              if (lev == bEns(instanceIndex)%nLevEns_T) then
                ! use surface momentum level amplitudes for surface thermo level
                ensAmplitude_oneLev => ens_getOneLev_r8(ensAmplitude,bEns(instanceIndex)%nLevEns_M)
                ensAmplitude_oneLev(1:bEns(instanceIndex)%nEns,:,lonIndex,latIndex) = &
                     ensAmplitude_oneLev(1:bEns(instanceIndex)%nEns,:,lonIndex,latIndex) + ensAmplitude_MT(:,:)
              else
                ! for other levels, interpolate momentum weights to get thermo amplitudes
                ensAmplitude_oneLev   => ens_getOneLev_r8(ensAmplitude,lev)
                ensAmplitude_oneLevP1 => ens_getOneLev_r8(ensAmplitude,lev+1)
                ensAmplitude_oneLev(1:bEns(instanceIndex)%nEns,:,lonIndex,latIndex)   = &
                     ensAmplitude_oneLev(1:bEns(instanceIndex)%nEns,:,lonIndex,latIndex)   + 0.5d0*ensAmplitude_MT(:,:)
                ensAmplitude_oneLevP1(1:bEns(instanceIndex)%nEns,:,lonIndex,latIndex) = &
                     ensAmplitude_oneLevP1(1:bEns(instanceIndex)%nEns,:,lonIndex,latIndex) + 0.5d0*ensAmplitude_MT(:,:)
              end if

            else
              call utl_abort('ben_addEnsMemberAd: incompatible vcode')
            end if

          else if (varLevel == 'SF') then
            ! Surface variable cases

            if (bEns(instanceIndex)%vco_anl%Vcode == 5002 .or. bEns(instanceIndex)%vco_anl%Vcode == 5005) then
              ensAmplitude_oneLev   => ens_getOneLev_r8(ensAmplitude,bEns(instanceIndex)%nLevEns_M)
            else ! vco_anl%Vcode == 0
              ensAmplitude_oneLev   => ens_getOneLev_r8(ensAmplitude,1)
            end if
            ensAmplitude_oneLev (1:bEns(instanceIndex)%nEns,:,lonIndex,latIndex) = &
                 ensAmplitude_oneLev(1:bEns(instanceIndex)%nEns,:,lonIndex,latIndex) + ensAmplitude_MT(:,:)

          else
            call utl_abort('ben_addEnsMemberAd: unknown value of varLevel')
          end if

        end do ! lonIndex
      end do ! latIndex
      !$OMP END PARALLEL DO

    end do ! levIndex

    deallocate(ensAmplitude_MT)
    deallocate(increment_in2)

    call utl_tmg_stop(62)

  end subroutine addEnsMemberAd

  !--------------------------------------------------------------------------
  ! ensembleDiagnostic
  !--------------------------------------------------------------------------
  subroutine ensembleDiagnostic(instanceIndex, mode)
    implicit none

    integer, intent(in) :: instanceIndex
    character(len=*), intent(in) :: mode

    type(struct_gsv) :: statevector, statevector_temp

    integer :: nWaveBandToDiagnose, waveBandIndex, memberIndex

    real(8) :: dnens2

    character(len=48):: fileName
    character(len=12):: etiket
    character(len=2) :: waveBandNumber
    character(len=2) :: instanceNumber

    if ( trim(mode) == 'FullPerturbations') then
      nWaveBandToDiagnose = 1
    else if ( trim(mode) == 'WaveBandPerturbations' ) then
      if (trim(bEns(instanceIndex)%localizationType) == 'ScaleDependent') then
        nWaveBandToDiagnose = bEns(instanceIndex)%nWaveBand
      else
        nWaveBandToDiagnose = 1
      end if
    else
      write(*,*)
      write(*,*) 'mode = ', trim(mode)
      call utl_abort('EnsembleDiagnostic: unknown mode')
    end if

    if ( mmpi_myid == 0 ) write(*,*)
    if ( mmpi_myid == 0 ) write(*,*) 'EnsembleDiagnostic in mode: ', mode

    !
    !- Write each wave band for a selected member
    !
    if ( mmpi_myid == 0 ) write(*,*) '   writing perturbations for member 001'
    memberIndex = 1
    dnens2 = sqrt(1.0d0*dble(bEns(instanceIndex)%nEns-1))
    do waveBandIndex = 1, nWaveBandToDiagnose
      if ( mmpi_myid == 0 ) write(*,*) '     waveBandIndex = ', waveBandIndex
      call gsv_allocate(statevector, tim_nstepobsinc, bEns(instanceIndex)%hco_ens, bEns(instanceIndex)%vco_anl, &
                        datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true., &
                        allocHeight_opt=.false., allocPressure_opt=.false.)
      call ben_getPerturbation( statevector,    & ! OUT
                                memberIndex,    & ! IN
                                'ConstantValue', waveBandIndex, instanceIndex_opt = instanceIndex) ! IN
      if ( trim(mode) == 'FullPerturbations') then
        etiket = 'PERT001_FULL'
      else
        write(waveBandNumber,'(I2.2)') waveBandIndex
        etiket = 'PERT001_WB' // trim(waveBandNumber)
      end if
      write(instanceNumber,'(I2.2)') instanceIndex
      fileName = './ens_pert001_i' // trim(instanceNumber) // '.fst'

      call gio_writeToFile(statevector,fileName,etiket, &                               ! IN
                           scaleFactor_opt=dnens2, &                                    ! IN
                           HUcontainsLQ_opt=bEns(instanceIndex)%gsvHUcontainsLQ )       ! IN
      call gsv_deallocate(statevector)
    end do

    !
    !- Compute the standard deviations for each wave band
    !
    if ( mmpi_myid == 0 ) write(*,*) '   computing Std.Dev.'
    call gsv_allocate(statevector_temp, tim_nstepobsinc, bEns(instanceIndex)%hco_ens, bEns(instanceIndex)%vco_anl, &
                      mpi_local_opt=.true., dataKind_opt=ens_getDataKind(bEns(instanceIndex)%ensPerts(1)), &
                      allocHeight_opt=.false., allocPressure_opt=.false.)

    do waveBandIndex = 1, nWaveBandToDiagnose
       if ( mmpi_myid == 0 ) write(*,*) '     waveBandIndex = ', waveBandIndex
       call gsv_allocate(statevector, tim_nstepobsinc, bEns(instanceIndex)%hco_ens, bEns(instanceIndex)%vco_anl, &
                         datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true., &
                         dataKind_opt=ens_getDataKind(bEns(instanceIndex)%ensPerts(1)), &
                         allocHeight_opt=.false., allocPressure_opt=.false.)
       call gsv_zero(statevector)
       do memberIndex = 1, bEns(instanceIndex)%nEns
          !- Get normalized perturbations
          call ens_copyMember(bEns(instanceIndex)%ensPerts(waveBandIndex), & ! IN
                              statevector_temp,        & ! OUT
                              memberIndex)               ! IN

          !- Square
          call gsv_power(statevector_temp, & ! INOUT
                         2.d0)               ! IN
          !- Sum square values, result in statevector
          call gsv_add(statevector_temp, & ! IN
                       statevector)        ! INOUT
       end do

       !- Convert to StdDev
       call gsv_power(statevector, & ! INOUT
                      0.5d0)         ! IN

       !- Write to file
       if ( trim(mode) == 'FullPerturbations') then
          etiket = 'STDDEV_FULL'
       else
          write(waveBandNumber,'(I2.2)') waveBandIndex
          etiket = 'STDDEV_WB' // trim(waveBandNumber)
       end if
       write(instanceNumber,'(I2.2)') instanceIndex
       fileName = './ens_stddev_i' // trim(instanceNumber) // '.fst'

       call gio_writeToFile(statevector,fileName,etiket,                         & ! IN
                            HUcontainsLQ_opt=bEns(instanceIndex)%gsvHUcontainsLQ)  ! IN
       call gsv_deallocate(statevector)
    end do

    call gsv_deallocate(statevector_temp)

  end subroutine ensembleDiagnostic

  !--------------------------------------------------------------------------
  ! ben_writeAmplitude
  !--------------------------------------------------------------------------
  subroutine ben_writeAmplitude(ensPathName, ensFileNamePrefix, ip3, instanceIndex_opt)
    implicit none

    character(len=*),  intent(in) :: ensPathName
    character(len=*),  intent(in) :: ensFileNamePrefix
    integer,           intent(in) :: ip3
    integer, optional, intent(in) :: instanceIndex_opt

    integer :: instanceIndex

    instanceIndex = ben_setInstanceIndex(instanceIndex_opt)

    if (bEns(instanceIndex)%initialized .and. bEns(instanceIndex)%keepAmplitude) then
      if ( mmpi_myid == 0 ) write(*,*)
      if ( mmpi_myid == 0 ) write(*,*) 'bmatrixEnsemble_mod: Writing the amplitude field'
      call ens_writeEnsemble(bEns(instanceIndex)%ensAmplitudeStorage, ensPathName, ensFileNamePrefix, &
                             'FROM_BENS', 'R',varNames_opt=bEns(instanceIndex)%varNameALFA, ip3_opt=ip3)
    end if

  end subroutine ben_writeAmplitude

  !--------------------------------------------------------------------------
  ! ben_setFsoLeadTime
  !--------------------------------------------------------------------------
  subroutine ben_setFsoLeadTime(fsoLeadTime_in, instanceIndex_opt)
    implicit none
    real(8)  :: fsoLeadTime_in
    integer, optional, intent(in) :: instanceIndex_opt
    integer :: instanceIndex

    instanceIndex = ben_setInstanceIndex(instanceIndex_opt)

    bEns(instanceIndex)%fsoLeadTime = fsoLeadTime_in

  end subroutine ben_setFsoLeadTime

  !--------------------------------------------------------------------------
  ! ben_getNumStepAmplitudeAssimWindow
  !--------------------------------------------------------------------------
  function ben_getNumStepAmplitudeAssimWindow(instanceIndex_opt) result(numStepAmplitude)
    implicit none
    integer, optional, intent(in) :: instanceIndex_opt
    integer :: instanceIndex
    integer numStepAmplitude

    instanceIndex = ben_setInstanceIndex(instanceIndex_opt)

    numStepAmplitude = bEns(instanceIndex)%numStepAmplitudeAssimWindow

  end function ben_getNumStepAmplitudeAssimWindow

  !--------------------------------------------------------------------------
  ! ben_getAmplitudeAssimWindow
  !--------------------------------------------------------------------------
  function ben_getAmplitudeAssimWindow(instanceIndex_opt) result(adv_amplitude)
    implicit none
    integer, optional, intent(in) :: instanceIndex_opt
    integer :: instanceIndex
    type(struct_adv), pointer  :: adv_amplitude

    instanceIndex = ben_setInstanceIndex(instanceIndex_opt)

    adv_amplitude => bEns(instanceIndex)%adv_amplitudeAssimWindow

  end function ben_getAmplitudeAssimWindow

  !--------------------------------------------------------------------------
  ! ben_getAmp3dStepIndexAssimWindow
  !--------------------------------------------------------------------------
  function ben_getAmp3dStepIndexAssimWindow(instanceIndex_opt) result(stepIndex)
    implicit none
    integer, optional, intent(in) :: instanceIndex_opt
    integer :: instanceIndex
    integer  :: stepIndex

    instanceIndex = ben_setInstanceIndex(instanceIndex_opt)

    stepIndex = bEns(instanceIndex)%amp3dStepIndexAssimWindow

  end function ben_getAmp3dStepIndexAssimWindow

  !--------------------------------------------------------------------------
  ! ben_getNumInstance
  !--------------------------------------------------------------------------
  function ben_getNumInstance() result(numInstance)
    implicit none
    integer  :: numInstance

    numInstance = nInstance

  end function ben_getNumInstance

  !--------------------------------------------------------------------------
  ! ben_getNumLoc
  !--------------------------------------------------------------------------
  function ben_getNumLoc(instanceIndex_opt) result(numLoc)
    implicit none
    integer, optional, intent(in) :: instanceIndex_opt
    integer :: instanceIndex
    integer  :: numLoc

    instanceIndex = ben_setInstanceIndex(instanceIndex_opt)

    numLoc = bEns(instanceIndex)%nWaveBand

  end function ben_getNumLoc

  !--------------------------------------------------------------------------
  ! ben_getLoc
  !--------------------------------------------------------------------------
  function ben_getLoc(locIndex,instanceIndex_opt) result(loc)
    implicit none
    integer, intent(in)  :: locIndex
    integer, optional, intent(in) :: instanceIndex_opt
    integer :: instanceIndex
    type(struct_loc), pointer :: loc

    instanceIndex = ben_setInstanceIndex(instanceIndex_opt)

    loc => bEns(instanceIndex)%locStorage(locIndex)

  end function ben_getLoc

end module BMatrixEnsemble_mod
