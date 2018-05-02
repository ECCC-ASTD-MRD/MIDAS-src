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

!--------------------------------------------------------------------------
!! MODULE BmatrixEnsemble (prefix="ben")
!!
!! *Purpose*: Performs transformation from control vector to analysis increment 
!!            using the spatially localized ensemble covariance matrix. This 
!!            module works for both global and limited-area applications.
!!
!--------------------------------------------------------------------------
MODULE BmatrixEnsemble_mod
  use ramDisk_mod
  use mpivar_mod
  use fileNames_mod
  use gridStateVector_mod
  use ensembleStateVector_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use localization_mod
  use mathPhysConstants_mod
  use earthConstants_mod
  use variableTransforms_mod
  use utilities_mod
  use globalSpectralTransform_mod
  use lamSpectralTransform_mod
  use spectralFilter_mod
  use varNameList_mod
  use advection_mod
  implicit none
  save
  private

  ! public procedures
  public :: ben_Setup, ben_BSqrt, ben_BSqrtAd, ben_writeAmplitude
  public :: ben_reduceToMPILocal, ben_reduceToMPILocal_r4, ben_expandToMPIGlobal, ben_expandToMPIGlobal_r4
  public :: ben_getScaleFactor, ben_getnEns, ben_getPerturbation, ben_getEnsMean, ben_Finalize
  public :: ben_setFsoLeadTime, ben_getNumStepAmplitudeAssimWindow, ben_getAmplitudeAssimWindow
  public :: ben_getAmp3dStepIndexAssimWindow

  logical             :: initialized = .false.

  integer,parameter   :: maxNumLevels=200

  real(8),allocatable :: scaleFactor_M(:), scaleFactor_T(:)

  integer             :: nj,ni,lonPerPE,lonPerPEmax,myLonBeg,myLonEnd,latPerPE,latPerPEmax,myLatBeg,myLatEnd
  integer,allocatable :: allLonBeg(:), allLatBeg(:)
  integer             :: nLevInc_M,nLevInc_T,nLevEns_M,nLevEns_T
  integer             :: topLevIndex_M,topLevIndex_T
  integer             :: nEnsOverDimension
  integer             :: cvDim_mpilocal,cvDim_mpiglobal
  integer             :: numStep, numStepAssimWindow
  integer             :: numStepAmplitudeFSOFcst, numStepAmplitudeAssimWindow, numStepAdvectAssimWindow
  integer             :: numSubEns
  integer,allocatable :: dateStampList(:)
  integer,allocatable :: dateStampListAdvectedFields(:)

  character(len=32)   :: direction, directionEnsPerts, directionAnlInc

  integer,external    :: get_max_rss, omp_get_thread_num
  integer             :: numIncludeAnlVar

  ! FSO
  real(8)             :: fsoLeadTime = -1.0D0
  integer             :: numStepAdvectFSOFcst

  ! Localizations
  integer, parameter  :: maxNumLocalLength = 20
  integer             :: nWaveBand

  ! Ensemble perturbations
  type(struct_ens), allocatable :: ensPerts(:)

  ! Ensemble amplitude (only used in diagnostic mode)
  type(struct_ens)    :: ensAmplitudeStorage
  character(len=4), parameter  :: varNameALFA(1) = (/ 'ALFA' /)

  ! Localization
  integer, allocatable :: locIDs(:)

  logical            :: HUcontainsLQ_gsv

  ! Vertical grid
  type(struct_vco),pointer :: vco_anl, vco_ens, vco_file => null()

  ! Horizontal grid
  type(struct_hco), pointer :: hco_anl  ! Analysis   horizontal grid parameters
  type(struct_hco), pointer :: hco_ens  ! Ensemble   horizontal grid parameters
  type(struct_hco), pointer :: hco_file ! Input file horizontal grid parameters

  ! Amplitude parameters
  character(len=32) :: advectTypeAssimWindow
  character(len=32) :: advectStartTimeIndexAssimWindow

  type(struct_adv)          :: adv_amplitudeFSOFcst
  type(struct_adv), target  :: adv_amplitudeAssimWindow
  type(struct_adv)          :: adv_ensPerts
  type(struct_adv)          :: adv_analInc

  integer           :: amp3dStepIndexAssimWindow
  integer           :: amp3dStepIndexFSOFcst

  ! Namelist variables
  integer             :: nEns ! number of ensemble members
  real(8)             :: scaleFactor(maxNumLevels)
  real(8)             :: scaleFactorHumidity(maxNumLevels)
  integer             :: ntrunc
  character(len=256)  :: enspathname
  real(8)             :: hLocalize(maxNumLocalLength)
  real(8)             :: vLocalize(maxNumLocalLength)
  character(len=256)  :: LocalizationType
  integer             :: waveBandPeaks(maxNumLocalLength)
  logical             :: ensDiagnostic
  logical             :: advDiagnostic
  character(len=2)    :: ctrlVarHumidity
  logical             :: advectAmplitudeFSOFcst
  logical             :: advectAmplitudeAssimWindow = .false.
  logical             :: advectEnsPertAnlInc        = .false.
  logical             :: removeSubEnsMeans
  logical             :: keepAmplitude
  character(len=4)    :: IncludeAnlVar(vnl_numvarmax)
  
  ! Control parameter for the level of listing output
  logical, parameter :: verbose = .false.

CONTAINS

  !--------------------------------------------------------------------------
  ! ben_setup
  !--------------------------------------------------------------------------
  SUBROUTINE ben_setup(hco_anl_in, vco_anl_in, cvDim_out, &
                       mode_opt)
    implicit none

    type(struct_hco), pointer, intent(in) :: hco_anl_in
    type(struct_vco), pointer, intent(in) :: vco_anl_in

    character(len=*), intent(in), optional :: mode_opt

    character(len=15) :: ben_mode

    type(struct_gsv) :: statevector_ensMean4D, statevector_oneEnsPert4D

    real(8) :: pSurfRef, delT_hour
    real(8) :: advectFactorFSOFcst
    real(8) :: advectFactorAssimWindow

    real(8),pointer :: pressureProfileEns_M(:), pressureProfileFile_M(:), pressureProfileInc_M(:)

    integer        :: cvDim_out, myMemberBeg,myMemberEnd,myMemberCount,maxMyMemberCount
    integer        :: levIndex,nIndex,mIndex,jvar,ila,return_code,status
    integer        :: fnom,fclos,ierr,nulnam
    integer        :: waveBandIndex,locID, stepIndex
    integer        :: stamp_last,newdate,ndate,ntime
    character(len=256) :: ensFileName
    integer        :: dateStampFSO

    logical        :: EnsTopMatchesAnlTop, useAnlLevelsOnly
    logical        :: lExists

    !namelist
    NAMELIST /NAMBEN/nEns,scaleFactor,scaleFactorHumidity,ntrunc,enspathname, &
         hLocalize,vLocalize,LocalizationType,waveBandPeaks, &
         ensDiagnostic,advDiagnostic,ctrlVarHumidity,advectFactorFSOFcst,advectFactorAssimWindow,&
         removeSubEnsMeans, keepAmplitude, advectTypeAssimWindow, advectStartTimeIndexAssimWindow, &
         IncludeAnlVar

    if (verbose) write(*,*) 'Entering ben_Setup'

    call tmg_start(12,'BEN_SETUP')

    !
    !- 1.1  Read namelist-dependent options
    !

    ! parameters from namelist
    scaleFactor(:)        =    0.0d0
    scaleFactorHumidity(:)=    1.0d0
    nEns                  =   10
    ntrunc                =   30
    enspathname           = 'ensemble'
    LocalizationType      = 'LevelDependent'
    waveBandPeaks(:)      =   -1.0d0
    ensDiagnostic         = .false.
    advDiagnostic         = .false.
    hLocalize(:)          =   -1.0d0
    hLocalize(1)          = 2800.0d0
    vLocalize(:)          =   -1.0d0
    vLocalize(1)          =    2.0d0
    ctrlVarHumidity       = 'LQ'
    advectFactorFSOFcst   =   0.0D0
    advectTypeAssimWindow   = 'amplitude'
    advectStartTimeIndexAssimWindow = 'first'
    advectFactorAssimWindow =   0.0D0
    removeSubEnsMeans     = .false.
    keepAmplitude         = .false. 

    includeAnlVar(:) = ''
    numIncludeAnlVar = 0

    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=namben,iostat=ierr)
    if (ierr.ne.0) call utl_abort('ben_setup: Error reading namelist')
    if (mpi_myid == 0) write(*,nml=namben)
    ierr = fclos(nulnam)
 
    ! If zero weight, skip rest of setup
    if ( sum(scaleFactor(:)) == 0.0d0 ) then
      if (mpi_myid == 0) write(*,*) 'ben_setup: scaleFactor=0, skipping rest of setup'
      cvDim_out = 0
      return
    end if

    write(*,*) 'ben_setup: enspathname = ', trim(enspathname)

    !
    !- 2.  Settings
    !

    !- 2.1 Mode
    if ( present(mode_opt) ) then
      if ( trim(mode_opt) == 'Analysis' .or. trim(mode_opt) == 'BackgroundCheck') then
        ben_mode = trim(mode_opt)
        if (mpi_myid == 0) write(*,*)
        if (mpi_myid == 0) write(*,*) 'ben_setup: Mode activated = ', trim(ben_mode)
      else
        write(*,*)
        write(*,*) 'mode = ', trim(mode_opt)
        call utl_abort('ben_setup: unknown mode')
      end if
    else
      ben_mode = 'Analysis'
      if (mpi_myid == 0) write(*,*)
      if (mpi_myid == 0) write(*,*) 'ben_setup: Analysis mode activated (by default)'
    end if

    !- 2.2 Number of time step bins
    numStep = tim_nstepobsinc
    if (numStep /= 1.and.numStep /= 3.and.numStep /= 5.and.numStep /= 7) then
      call utl_abort('ben_setup: Invalid value for numStep (choose 1 or 3 or 5 or 7)!')
    end if

    !- for FSO
    numStepAssimWindow = numStep
    if (fsoLeadTime > 0.0D0) then
      numStep = numStep + 1
      call incdatr(dateStampFSO, tim_getDatestamp(), fsoLeadTime)
    end if

    allocate(dateStampList(numStep))
    if (fsoLeadTime > 0.0D0) then
      call tim_getstamplist(dateStampList,numStep-1,tim_getDatestamp())
      dateStampList(numStep) = dateStampFSO
    else
      call tim_getstamplist(dateStampList,numStep,tim_getDatestamp())
    end if

    !- 2.3 Horizontal grid
    hco_anl => hco_anl_in
    hco_ens => hco_anl ! ensemble members must be on analysis grid
    ni = hco_ens%ni
    nj = hco_ens%nj
    if (hco_ens%global) then
      if (mpi_myid == 0) write(*,*)
      if (mpi_myid == 0) write(*,*) 'ben_setup: GLOBAL mode activated'
    else
      if (mpi_myid == 0) write(*,*)
      if (mpi_myid == 0) write(*,*) 'ben_setup: LAM mode activated'
    end if

    !- 2.4 Vertical levels
    vco_anl => vco_anl_in

    if ( mpi_myid == 0 ) then
      call fln_ensfileName(ensFileName, ensPathName, 1)
      call vco_SetupFromFile(vco_file, ensFileName)
    end if
    call vco_mpiBcast(vco_file)

    !- Do we need to read all the vertical levels from the ensemble?
    useAnlLevelsOnly = vco_subsetOrNot(vco_anl, vco_file)
    if ( useAnlLevelsOnly ) then
      write(*,*)
      write(*,*) 'ben_setup: only the analysis levels will be read in the ensemble '
      vco_ens  => vco_anl ! the ensemble target grid is the analysis grid
      call vco_deallocate(vco_file)
      vco_file => vco_anl ! only the analysis levels will be read in the ensemble
      EnsTopMatchesAnlTop = .true.
    else
      write(*,*)
      write(*,*) 'ben_setup: all the vertical levels will be read in the ensemble '
      pSurfRef = 101000.D0
      nullify(pressureProfileInc_M)
      status = vgd_levels( vco_anl%vgrid, ip1_list=vco_anl%ip1_M, levels=pressureProfileInc_M, &
           sfc_field=pSurfRef, in_log=.false.)
      if (status /= VGD_OK) call utl_abort('ben_setup: ERROR from vgd_levels')
      nullify(pressureProfileFile_M)
      status = vgd_levels( vco_file%vgrid, ip1_list=vco_file%ip1_M, levels=pressureProfileFile_M, &
           sfc_field=pSurfRef, in_log=.false.)
      if (status /= VGD_OK) call utl_abort('ben_setup: ERROR from vgd_levels')
      
      EnsTopMatchesAnlTop = abs( log(pressureProfileFile_M(1)) - log(pressureProfileInc_M(1)) ) < 0.1d0
      write(*,*) 'EnsTopMatchesAnlTop: EnsTopMatchesAnlTop, presEns, presInc = ', &
           EnsTopMatchesAnlTop, pressureProfileFile_M(1), pressureProfileInc_M(1)
      deallocate(pressureProfileFile_M)
      deallocate(pressureProfileInc_M)

      if ( EnsTopMatchesAnlTop ) then
        if ( mpi_myid == 0 ) write(*,*) 'ben_setup: top level of ensemble member and analysis grid match'
        vco_ens => vco_anl  ! IMPORTANT: top levels DO match, therefore safe
        ! to force members to be on analysis vertical levels
      else
        if ( mpi_myid == 0 ) write(*,*) 'ben_setup: top level of ensemble member and analysis grid are different, therefore'
        if ( mpi_myid == 0 ) write(*,*) '           assume member is already be on correct levels - NO CHECKING IS DONE'
        vco_ens => vco_file ! IMPORTANT: top levels do not match, therefore must
        ! assume file is already on correct vertical levels
      end if
    end if
    
    if (vco_anl%Vcode.ne.vco_ens%Vcode) then
      write(*,*) 'ben_setup: vco_anl%Vcode = ', vco_anl%Vcode, ', vco_ens%Vcode = ', vco_ens%Vcode
      call utl_abort('ben_setup: vertical levels of ensemble not compatible with analysis grid')
    end if
    nLevEns_M = vco_ens%nLev_M
    nLevEns_T = vco_ens%nLev_T
    nLevInc_M = vco_anl%nLev_M
    nLevInc_T = vco_anl%nLev_T
    topLevIndex_M = nLevInc_M-nLevEns_M+1
    topLevIndex_T = nLevInc_T-nLevEns_T+1

    if (vco_anl%Vcode == 5002) then
      if ( (nLevEns_T /= (nLevEns_M+1)) .and. (nLevEns_T /= 1 .or. nLevEns_M /= 1) ) then
        write(*,*) 'ben_setup: nLevEns_T, nLevEns_M = ',nLevEns_T,nLevEns_M
        call utl_abort('ben_setup: Vcode=5002, nLevEns_T must equal nLevEns_M+1!')
      end if
    else if (vco_anl%Vcode == 5005) then
      if (nLevEns_T.ne.nLevEns_M) then
        write(*,*) 'ben_setup: nLevEns_T, nLevEns_M = ',nLevEns_T,nLevEns_M
        call utl_abort('ben_setup: Vcode=5005, nLevEns_T must equal nLevEns_M!')
      end if
    else
      write(*,*) 'vco_anl%Vcode = ',vco_anl%Vcode
      call utl_abort('ben_setup: unknown vertical coordinate type!')
    end if

    if (nLevEns_M.gt.nLevInc_M) then
      call utl_abort('ben_setup: ensemble has more levels than increment - not allowed!')
    end if

    if (nLevEns_M.lt.nLevInc_M) then
      if (mpi_myid == 0) write(*,*) 'ben_setup: ensemble has less levels than increment'
      if (mpi_myid == 0) write(*,*) '           some levels near top will have zero increment'
    end if

    !- 2.5 Bmatrix Weight
    allocate(scaleFactor_M(nLevEns_M))
    allocate(scaleFactor_T(nLevEns_T))
    do levIndex = 1, nLevEns_T
      if (scaleFactor(levIndex).gt.0.0d0) then 
        scaleFactor(levIndex) = sqrt(scaleFactor(levIndex))
      else
        scaleFactor(levIndex) = 0.0d0
      end if
    end do
    scaleFactor_T(1:nLevEns_T) = scaleFactor(1:nLevEns_T)
    if (vco_anl%Vcode == 5002) then
      scaleFactor_M(1:nLevEns_M) = scaleFactor(2:(nLevEns_M+1))
    else
      scaleFactor_M(1:nLevEns_M) = scaleFactor(1:nLevEns_M)
    end if

    do levIndex = 1, nLevEns_T
      if (scaleFactorHumidity(levIndex).gt.0.0d0) then 
        scaleFactorHumidity(levIndex) = sqrt(scaleFactorHumidity(levIndex))
      else
        scaleFactorHumidity(levIndex) = 0.0d0
      end if
    end do

    !- 2.5 Domain Partionning
    call mpivar_setup_latbands(nj, latPerPE, latPerPEmax, myLatBeg, myLatEnd)
    call mpivar_setup_lonbands(ni, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd)
    allocate(allLonBeg(mpi_npex))
    CALL rpn_comm_allgather(myLonBeg,1,"mpi_integer",       &
         allLonBeg,1,"mpi_integer","EW",ierr)
    allocate(allLatBeg(mpi_npey))
    CALL rpn_comm_allgather(myLatBeg,1,"mpi_integer",       &
         allLatBeg,1,"mpi_integer","NS",ierr)

    !- 2.6 Localization
    if ( trim(ben_mode) == 'Analysis' ) then

      call mpivar_setup_levels_npex(nEns,myMemberBeg,myMemberEnd,myMemberCount)
      call rpn_comm_allreduce(myMemberCount, maxMyMemberCount, &
           1,"MPI_INTEGER","MPI_MAX","GRID",ierr)
      nEnsOverDimension = mpi_npex * maxMyMemberCount

      if (trim(LocalizationType) == 'LevelDependent') then
        if (mpi_myid == 0) write(*,*)
        if (mpi_myid == 0) write(*,*) 'ben_setup: Level-Dependent (Standard) localization will be used'
        nWaveBand = 1
      else if (trim(LocalizationType) == 'ScaleDependent') then
        if (mpi_myid == 0) write(*,*)
        if (mpi_myid == 0) write(*,*) 'ben_setup: Scale-Dependent localization will be used'
        nWaveBand = count(waveBandPeaks .ge. 0)
        if ( nWaveBand <= 1 ) then
          call utl_abort('ben_setup: nWaveBand <= 1')
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
        do waveBandIndex = 1, nWaveBand-1
          if ( waveBandPeaks(waveBandIndex)-waveBandPeaks(waveBandIndex+1) <= 0 ) then
            call utl_abort('ben_setup: waveBandPeaks are not in decreasing wavenumber order') 
          end if
        end do

        ! Make sure that we have valid localization length scales for each wave bands
        do  waveBandIndex = 1, nWaveBand
          if ( hLocalize(waveBandIndex) <= 0.0d0 ) then
            call utl_abort('ben_setup: Invalid HORIZONTAL localization length scale')
          end if
          if ( vLocalize(waveBandIndex) <= 0.0d0 ) then
            call utl_abort('ben_setup: Invalid VERTICAL localization length scale')
          end if
        end do

        ! Make sure the truncation is compatible with the waveBandPeaks
        if ( ntrunc < waveBandPeaks(1) ) then
          call utl_abort('ben_setup: The truncation is not compatible with the your scale-dependent localization')
        end if

      else
        call utl_abort('ben_setup: Invalid mode for LocalizationType')
      end if

      pSurfRef = 101000.D0
      nullify(pressureProfileInc_M)
      status = vgd_levels( vco_anl%vgrid, ip1_list=vco_anl%ip1_M, levels=pressureProfileInc_M, &
           sfc_field=pSurfRef, in_log=.false.)
      if (status /= VGD_OK)then
        call utl_abort('ben_setup: ERROR from vgd_levels')
      end if

      allocate(pressureProfileEns_M(nLevEns_M))
      pressureProfileEns_M(1:nLevEns_M) = pressureProfileInc_M(topLevIndex_M:nLevInc_M)

      allocate(locIDs(nWaveBand))
      do waveBandIndex = 1, nWaveBand
        call loc_setup(hco_ens, vco_ens, nEns, pressureProfileEns_M, nTrunc, 'spectral', & ! IN
             LocalizationType, hLocalize(waveBandIndex), hLocalize(waveBandIndex+1),     & ! IN
             vLocalize(waveBandIndex),                                                   & ! IN
             cvDim_mpilocal, locID)                                                        ! OUT
        locIDs(waveBandIndex) = locID
      end do

      cvDim_out = cvDim_mpilocal

      deallocate(pressureProfileEns_M)
      deallocate(pressureProfileInc_M)
    end if

    !- 2.7 Control variables
    if      ( ctrlVarHumidity == 'LQ' ) then
      write(*,*)
      write(*,*) 'ben_setup: Humidity control variable = ', ctrlVarHumidity
      HUcontainsLQ_gsv = .true.
    else if ( ctrlVarHumidity == 'HU' ) then
      write(*,*)
      write(*,*) 'ben_setup: Humidity control variable = ', ctrlVarHumidity
      HUcontainsLQ_gsv = .false.
    else
      write(*,*)
      write(*,*) 'Unknown humidity control variable'
      write(*,*) 'Should be LQ or LU, found = ', ctrlVarHumidity
      call utl_abort('ben_setup')
    end if

    !
    !- 3.  Read/Process the Ensemble
    !
    
    !- 3.1 Identify set of variables for which ensembles are required    
    do jvar = 1, vnl_numvarmax
      if (trim(includeAnlVar(jvar)) == '') exit
      if (.not.gsv_varExist(varName = trim(includeAnlVar(jvar)))) then
        write(*,*) 'ben_setup: This variable is not a member of ANLVAR: ', trim(includeAnlVar(jvar))
        call utl_abort('ben_setup: Invalid variable in includeAnlVar')
      else
        numIncludeAnlVar = numIncludeAnlVar+1
      end if
    end do
    if (numIncludeAnlVar == 0) then
      do jvar = 1, vnl_numvarmax
        if (.not.gsv_varExist(varName = trim(vnl_varNamelist(jvar)))) cycle
        numIncludeAnlVar = numIncludeAnlVar+1
        includeAnlVar(numIncludeAnlVar) = vnl_varNamelist(jvar)
      end do
    end if
    if (numIncludeAnlVar == 0) call utl_abort('ben_setup: Ensembles not being requested for any variable')
     
    !- 3.2 Read the ensemble data
    call setupEnsemble()

    if ( trim(ben_mode) /= 'Analysis' ) then
      cvDim_out = 9999 ! Dummy value > 0 to indicate to the background check (s/r compute_HBHT_ensemble)
      initialized = .true.
      call tmg_stop(12)
      return
    end if

    !- 3.3 Pre-compute everything for advection in FSO mode
    if (fsoLeadTime > 0.0D0) then
      amp3dStepIndexFSOFcst = 1
      if ( advectFactorFSOFcst == 0.0D0 .or. numStep == 1) then
        if (mpi_myid == 0) write(*,*) 'ben_setup: advection not activated for FSO'
        advectAmplitudeFSOFcst = .false.
        numStepAmplitudeFSOFcst = 1
      else
        if (mpi_myid == 0) write(*,*) 'ben_setup: advection activated in FSO mode'
        advectAmplitudeFSOFcst = .true.
        numStepAmplitudeFSOFcst = 2
        numStepAdvectFSOFcst = nint(fsoLeadTime/6.0D0) + 1
        allocate(dateStampListAdvectedFields(numStepAmplitudeFSOFcst))
        dateStampListAdvectedFields(1) = tim_getDatestamp()
        dateStampListAdvectedFields(2) = dateStampList(numStep)
        delT_hour = fsoLeadTime/real(numStepAdvectFSOFcst-1,8) ! time between winds
        call tmg_start(135,'BEN_SETUP_ADVEC_FSO')
        call adv_setup( adv_amplitudeFSOFcst,                                   & ! OUT
                        'fromFirstTimeIndex', hco_ens, vco_ens,                 & ! IN
                        numStepAmplitudeFSOFcst, dateStampListAdvectedFields,   & ! IN
                        numStepAdvectFSOFcst, delT_hour, advectFactorFSOFcst,   & ! IN
                        'MMLevsOnly',                                           & ! IN
                        steeringFlowFilename_opt=trim(enspathname)//'/forecast_for_advection' ) ! IN
        call tmg_stop(135)
      end if
    end if

    !- 3.4 Pre-compute everything for advection in ANALYSIS mode
    if ( advectFactorAssimWindow == 0.0D0 .or. numStep == 1) then
      if (mpi_myid == 0) write(*,*) 'ben_setup: advection not activated in ANALYSIS mode'

      advectAmplitudeAssimWindow = .false.
      numStepAmplitudeAssimWindow = 1
      amp3dStepIndexAssimWindow   = 1

    else
      if (mpi_myid == 0) write(*,*) 'ben_setup: advection activated in ANALYSIS mode'

      delT_hour                 = tim_dstepobsinc
      allocate(dateStampListAdvectedFields(numStep))
      dateStampListAdvectedFields(:) = dateStampList(:)
      call gsv_allocate(statevector_ensMean4D, numStep, hco_ens, vco_ens, &
                        datestampList_opt=dateStampListAdvectedFields,     &
                        mpi_local_opt=.true.)
      call ens_copyEnsMean(ensPerts(1), & ! IN
                           statevector_ensMean4D  )   ! OUT

      call tmg_start(136,'BEN_SETUP_ADVEC_ANL')

      select case(trim(advectTypeAssimWindow))
      case ('amplitude')
        if (mpi_myid == 0) write(*,*) '         amplitude fields will be advected'
        advectAmplitudeAssimWindow  = .true.
        numStepAmplitudeAssimWindow = numStep
        numStepAdvectAssimWindow    = numStep

        select case(trim(advectStartTimeIndexAssimWindow))
        case ('first')
          direction='fromFirstTimeIndex'
          amp3dStepIndexAssimWindow = 1
        case ('middle')
          direction='fromMiddleTimeIndex'
          amp3dStepIndexAssimWindow = (numStepAmplitudeAssimWindow+1)/2
        case default
          write(*,*)
          write(*,*) 'Unsupported starting timeIndex : ', trim(advectStartTimeIndexAssimWindow)
          call utl_abort('ben_setup')
        end select
        call adv_setup( adv_amplitudeAssimWindow,                                     & ! OUT
                        direction, hco_ens, vco_ens,                                  & ! IN
                        numStepAmplitudeAssimWindow, dateStampListAdvectedFields,     & ! IN
                        numStepAdvectAssimWindow, delT_hour, advectFactorAssimWindow, & ! IN
                        'MMLevsOnly', statevector_steeringFlow_opt = statevector_ensMean4D )       ! IN

      case('ensPertAnlInc')
        if (mpi_myid == 0) write(*,*) '         ensPerts and AnalInc will be advected'

        if (.not. EnsTopMatchesAnlTop) then
          call utl_abort('ben_setup: for advectTypeAssimWindow=ensPertAnlInc, ensTop and anlTop must match!')
        end if

        advectEnsPertAnlInc         = .true.
        amp3dStepIndexAssimWindow   = 1
        numStepAmplitudeAssimWindow = 1
        numStepAdvectAssimWindow    = numStep
        
        select case(trim(advectStartTimeIndexAssimWindow))
        case ('first')
          directionEnsPerts='towardFirstTimeIndex'
          directionAnlInc  ='towardFirstTimeIndexInverse'
        case ('middle')
          directionEnsPerts='towardMiddleTimeIndex'
          directionAnlInc  ='towardMiddleTimeIndexInverse'
        case default
          write(*,*)
          write(*,*) 'Unsupported starting timeIndex : ', trim(advectStartTimeIndexAssimWindow)
          call utl_abort('ben_setup')
        end select

        call adv_setup( adv_ensPerts,                                                 & ! OUT
                        directionEnsPerts, hco_ens, vco_ens,                          & ! IN
                        numStepAdvectAssimWindow, dateStampListAdvectedFields,        & ! IN
                        numStepAdvectAssimWindow, delT_hour, advectFactorAssimWindow, & ! IN
                        'allLevs', statevector_steeringFlow_opt=statevector_ensMean4D ) ! IN

        call adv_setup( adv_analInc,                                                  & ! OUT
                        directionAnlInc, hco_ens, vco_ens,                            & ! IN
                        numStepAdvectAssimWindow, dateStampListAdvectedFields,        & ! IN
                        numStepAdvectAssimWindow, delT_hour, advectFactorAssimWindow, & ! IN
                        'allLevs', statevector_steeringFlow_opt=statevector_ensMean4D ) ! IN

      case default
        write(*,*)
        write(*,*) 'Unsupported advectTypeAssimWindow : ', trim(advectTypeAssimWindow)
        call utl_abort('ben_setup')
      end select

      call tmg_stop(136)

      !- If wanted, write the ensemble mean
      if (advDiagnostic) then
        do stepIndex = 1, tim_nstepobsinc
          call gsv_writeToFile(statevector_ensMean4D,'./ens_mean.fst','ENSMEAN4D',        & ! IN
                               stepIndex_opt=stepIndex, HUcontainsLQ_opt=HUcontainsLQ_gsv ) ! IN
        end do
      end if

      call gsv_deallocate(statevector_ensMean4D)

    end if

    !- 3.5 Compute and write Std. Dev.
    if (ensDiagnostic) call EnsembleDiagnostic('FullPerturbations')

    !- 3.6 Ensemble perturbations advection
    if ( advectEnsPertAnlInc ) then

      !- If wanted, write the original ensemble perturbations for member #1
      if (advDiagnostic) then
        call ens_copyMember(ensPerts(1), statevector_oneEnsPert4D, 1)
        do stepIndex = 1, tim_nstepobsinc
          call gsv_writeToFile(statevector_oneEnsPert4D,'./ens_pert1.fst','ORIGINAL', & ! IN
               stepIndex_opt=stepIndex, HUcontainsLQ_opt=HUcontainsLQ_gsv )             ! IN
        end do
      end if

      !- Do the advection of all the members
      call tmg_start(137,'BEN_ADVEC_ENSPERT_TL')
      call adv_ensemble_tl( ensPerts(1), &       ! INOUT
                            adv_ensPerts, nEns ) ! IN
      call tmg_stop(137)

      !- If wanted, write the advected ensemble perturbations for member #1
      if (advDiagnostic) then
        call ens_copyMember(ensPerts(1), statevector_oneEnsPert4D, 1)
        do stepIndex = 1, tim_nstepobsinc
          call gsv_writeToFile(statevector_oneEnsPert4D,'./ens_pert1_advected.fst','ADVECTED', & ! IN
               stepIndex_opt=stepIndex,HUcontainsLQ_opt=HUcontainsLQ_gsv )                       ! IN
        end do
      end if

    end if

    !- 3.7 Compute and write Std. Dev.
    if (ensDiagnostic) call EnsembleDiagnostic('FullPerturbations')

    !- 3.8 Partitioned the ensemble perturbations into wave bands
    if (trim(LocalizationType) == 'ScaleDependent') then
      call EnsembleScaleDecomposition()
      if (ensDiagnostic) call EnsembleDiagnostic('WaveBandPerturbations')
    end if

    !- 3.9 Setup en ensGridStateVector to store the amplitude fields (for writing)
    if (keepAmplitude) then
      write(*,*)
      write(*,*) 'ben_setup: ensAmplitude fields will be store for potential write to file'
      call ens_allocate(ensAmplitudeStorage, nEns, numStepAmplitudeAssimWindow, hco_ens, vco_ens, &
                        dateStampList, varNames_opt=varNameALFA, dataKind_opt=8)
    end if

    !
    !- 4.  Ending
    !
    initialized = .true.

    call tmg_stop(12)

  END SUBROUTINE ben_setup

  !--------------------------------------------------------------------------
  ! ben_finalize
  !--------------------------------------------------------------------------
  SUBROUTINE ben_finalize()
    implicit none
    integer :: memberIndex, waveBandIndex, subEnsIndex

    if (verbose) write(*,*) 'Entering ben_Finalize'

    if (initialized) then
      write(*,*) 'ben_finalize: deallocating B_ensemble arrays'
      do waveBandIndex = 1, nWaveBand
        call ens_deallocate(ensPerts(waveBandIndex))
        call loc_finalize(locIDs(waveBandIndex))
      end do
      deallocate(ensPerts)
      if (keepAmplitude) call ens_deallocate(ensAmplitudeStorage)
    end if

  END SUBROUTINE ben_finalize

  !--------------------------------------------------------------------------
  ! ben_getScaleFactor
  !--------------------------------------------------------------------------
  subroutine ben_getScaleFactor(scaleFactor_out)
    implicit none
    real(8) :: scaleFactor_out(:)
    integer :: levIndex

    if (verbose) write(*,*) 'Entering ben_getScaleFactor'

    ! return value of 0 above highest level of ensemble
    do levIndex = 1, (topLevIndex_T - 1)
      scaleFactor_out(levIndex) = 0.0d0
    end do
    ! return scale factor for thermo levels
    do levIndex = topLevIndex_T, nLevInc_T
      scaleFactor_out(levIndex) = scaleFactor_T(levIndex-topLevIndex_T+1)
    end do

  end subroutine ben_getScaleFactor

  !--------------------------------------------------------------------------
  ! ben_getnEns
  !--------------------------------------------------------------------------
  integer function ben_getnEns()
    !func getnEns - returns the number ensemble member
    implicit none
    ben_getnEns = nEns
  end function ben_getnEns

  !--------------------------------------------------------------------------
  ! setupEnsemble
  !--------------------------------------------------------------------------
  SUBROUTINE setupEnsemble()
    implicit none

    real(4), pointer     :: ptr4d_r4(:,:,:,:)
    real(8) :: multFactor
    integer :: stepIndex,levIndex,lev,waveBandIndex,memberIndex
    integer :: ierr,varIndex
    logical :: makeBiPeriodic
    character(len=4) :: varName
    character(len=256) :: ensFileName

    write(*,*) 'setupEnsemble: Start'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !- 1. Memory allocation
    allocate(ensPerts(nWaveBand))    
    do waveBandIndex = 1, nWaveBand
      call ens_allocate(ensPerts(waveBandIndex), nEns, numStep, hco_ens, vco_ens, dateStampList, &
                         varNames_opt = includeAnlVar(1:numIncludeAnlVar))
    end do
    
    !- 2. Read ensemble
    makeBiPeriodic = (trim(LocalizationType) == 'ScaleDependent')
    call ens_readEnsemble(ensPerts(1), ensPathName, makeBiPeriodic, &
                          ctrlVarHumidity, vco_file_opt = vco_file, &
                          varNames_opt = includeAnlVar(1:numIncludeAnlVar))

    !- 3. From ensemble FORECASTS to ensemble PERTURBATIONS

    !- 3.1 remove mean
    call ens_computeMean( ensPerts(1), removeSubEnsMeans, numSubEns_opt=numSubEns )
    call ens_removeMean( ensPerts(1) )

    !- 3.2 normalize and apply scale factors
    !$OMP PARALLEL DO PRIVATE (levIndex,varName,lev,ptr4d_r4,stepIndex,memberIndex,multFactor)
    do levIndex = 1, ens_getNumK(ensPerts(1))
      varName = ens_getVarNameFromK(ensPerts(1),levIndex)
      lev = ens_getLevFromK(ensPerts(1),levIndex)

      if ( .not. ens_varExist(ensPerts(1), varName) ) cycle 

      ptr4d_r4 => ens_getOneLev_r4(ensPerts(1),levIndex)

      do stepIndex = 1, numStep
        do memberIndex = 1, nEns

          if ( vnl_varLevelFromVarname(varName) == 'MM' ) then
            multFactor = scaleFactor_M(lev)
          else if ( vnl_varLevelFromVarname(varName) == 'TH' ) then
            multFactor = scaleFactor_T(lev)
          else ! SF
            multFactor = scaleFactor_T(nLevEns_T)
          end if

          multFactor = multFactor/sqrt(1.0d0*dble(nEns-numSubEns))

          if (trim(varName) == 'HU') then
            multFactor = multFactor*scaleFactorHumidity(lev)
          end if

          ptr4d_r4(memberIndex,stepIndex,:,:) = real( real(ptr4d_r4(memberIndex,stepIndex,:,:),8)*multFactor, 4 )

        end do ! memberIndex
      end do ! stepIndex

    end do ! levIndex
    !$OMP END PARALLEL DO

    write(*,*) 'ben_setupEnsemble: finished adjusting ensemble members...'

  END SUBROUTINE setupEnsemble

  !--------------------------------------------------------------------------
  ! ben_getPerturbation
  !--------------------------------------------------------------------------
  SUBROUTINE ben_getPerturbation(statevector, memberIndexWanted,  &
       upwardExtrapolationMethod, waveBandIndexWanted_opt, &
       undoNormalization_opt)
    implicit none

    type(struct_gsv) :: statevector
    integer,          intent(in) :: memberIndexWanted
    character(len=*), intent(in) :: upwardExtrapolationMethod
    integer, optional, intent(in):: waveBandIndexWanted_opt
    logical, optional :: undoNormalization_opt

    real(8), pointer :: ptr4d_r8(:,:,:,:)
    real(4), pointer :: ensOneLev_r4(:,:,:,:)
    real(8) :: dnens2, scaleFactor_MT
    logical :: undoNormalization
    integer :: waveBandIndex
    integer :: lonIndex,latIndex,stepIndex,levIndex,lev,levInc,topLevOffset
    character(len=4) :: varName

    if (verbose) write(*,*) 'Entering ben_getPerturbation'

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

    do levIndex = 1, ens_getNumK(ensPerts(1))
      varName = ens_getVarNameFromK(ensPerts(1),levIndex)
      lev = ens_getLevFromK(ensPerts(1),levIndex)

      ptr4d_r8 => gsv_getField_r8(statevector, varName)
      ensOneLev_r4 => ens_getOneLev_r4(ensPerts(waveBandIndex),levIndex)

      !$OMP PARALLEL DO PRIVATE(stepIndex,topLevOffset,scaleFactor_MT,levInc,dnens2,latIndex,lonIndex)
      do stepIndex = 1, numStep

        if ( vnl_varLevelFromVarname(varName) == 'MM' ) then
          topLevOffset = topLevIndex_M - 1
          scaleFactor_MT = scaleFactor_M(lev)
        else if ( vnl_varLevelFromVarname(varName) == 'TH' ) then
          topLevOffset = topLevIndex_T - 1
          scaleFactor_MT = scaleFactor_T(lev)
        else ! SF
          topLevOffset = 0
          scaleFactor_MT = scaleFactor_T(nLevEns_T)
        end if

        levInc = lev + topLevOffset

        ! undo the normalization (optional)
        if (undoNormalization) then
          if (scaleFactor_MT > 0.0d0) then
            dnens2 = sqrt(1.0d0*dble(nEns-1))/scaleFactor_MT
          else
            if (stepIndex == 1) then 
              write(*,*) 'scalefactor not positive, cannot undo normalization!'
              write(*,*) varName,scaleFactor_MT,lev
            end if
            dnens2 = 0.0d0
          end if
          if (varName == 'HU  ') then
            if (scaleFactorHumidity(lev).gt.0.0d0) then
              dnens2 = dnens2/scaleFactorHumidity(lev)
            else
              if (stepIndex == 1) then
                write(*,*) 'Humidity scalefactor not positive, cannot undo normalization!'
                write(*,*) varName,scaleFactorHumidity(lev),lev
              end if
            end if
          end if
        else
          dnens2 = 1.0d0
        end if

        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            ptr4d_r8(lonIndex,latIndex,levInc,stepIndex) =   &
                 dnens2*dble(ensOneLev_r4(memberIndexWanted,stepIndex,lonIndex,latIndex))
          end do
        end do

        if ( topLevOffset > 0 .and. lev == 1) then
          ! Fill the gap between the ensemble lid and the analysis lid

          ! undo the normalization (optional)
          if (undoNormalization) then
            if (scaleFactor(1) > 0.0d0) then
              dnens2 = sqrt(1.0d0*dble(nEns-1))/scaleFactor(1)
            else
              if (stepIndex == 1) then
                write(*,*) 'scalefactor(top) not positive, cannot undo normalization!'
                write(*,*) varName,scaleFactor(1)
              end if
              dnens2 = 0.0d0
            end if
            if (varName == 'HU  ') then
              if (scaleFactorHumidity(1) > 0.0d0) then
                dnens2 = dnens2/scaleFactorHumidity(1)
              else
                if (stepIndex == 1) then
                  write(*,*) 'Humidity scalefactor(top) not positive, cannot undo normalization!'
                  write(*,*) varName,scaleFactorHumidity(1)
                end if
              end if
            end if
          else
            dnens2 = 1.0d0
          end if

          do levInc = 1, topLevOffset
            ! using a constant value
            do latIndex = myLatBeg, myLatEnd
              do lonIndex = myLonBeg, myLonEnd
                ptr4d_r8(lonIndex,latIndex,levInc,stepIndex) = dnens2 *  &
                     dble(ensOneLev_r4(memberIndexWanted,stepIndex,lonIndex,latIndex))
              end do
            end do
          end do

        end if ! topLevOffset > 0

      end do ! stepIndex
      !$OMP END PARALLEL DO

    end do ! levIndex

  END SUBROUTINE ben_getPerturbation

  !--------------------------------------------------------------------------
  ! ben_getEnsMean
  !--------------------------------------------------------------------------
  SUBROUTINE ben_getEnsMean(statevector, upwardExtrapolationMethod)
    implicit none

    type(struct_gsv) :: statevector
    character(len=*), intent(in) :: upwardExtrapolationMethod

    real(8), pointer :: ptr4d_out(:,:,:,:)
    real(8), pointer :: ensOneLev_mean(:,:,:)
    integer :: lonIndex,latIndex,stepIndex,levIndex,lev,levInc,topLevOffset
    character(len=4) :: varName

    if (.not. initialized) then
      if (mpi_myid == 0) write(*,*) 'ben_getEnsMean: bMatrixEnsemble not initialized, returning zero vector'
      call gsv_zero(statevector)
      return
    end if

    if ( trim(upwardExtrapolationMethod) /= "ConstantValue" ) then
      call utl_abort('ben_getEnsMean : Invalid value for upwardExtrapolationMethod')
    end if

    do levIndex = 1, ens_getNumK(ensPerts(1))
      varName = ens_getVarNameFromK(ensPerts(1),levIndex)
      lev = ens_getLevFromK(ensPerts(1),levIndex)

      ptr4d_out => gsv_getField_r8(statevector, varName)
      ensOneLev_mean => ens_getOneLevMean_r8(ensPerts(1), 1, levIndex)

      !$OMP PARALLEL DO PRIVATE(stepIndex,topLevOffset,levInc,latIndex,lonIndex)
      do stepIndex = 1, numStep

        if ( vnl_varLevelFromVarname(varName) == 'MM' ) then
          topLevOffset = topLevIndex_M - 1
        else if ( vnl_varLevelFromVarname(varName) == 'TH' ) then
          topLevOffset = topLevIndex_T - 1
        else ! SF
          topLevOffset = 0
        end if

        levInc = lev + topLevOffset

        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            ptr4d_out(lonIndex,latIndex,levInc,stepIndex) = ensOneLev_mean(stepIndex,lonIndex,latIndex)
          end do
        end do

        if ( topLevOffset > 0 .and. lev == 1 ) then
          ! Fill the gap between the ensemble lid and the analysis lid

          do levInc = 1, topLevOffset
            ! using a constant value
            do latIndex = myLatBeg, myLatEnd
              do lonIndex = myLonBeg, myLonEnd
                ptr4d_out(lonIndex,latIndex,levInc,stepIndex) = ensOneLev_mean(stepIndex,lonIndex,latIndex)
              end do
            end do
          end do

        end if ! topLevOffset > 0

      end do ! stepIndex
      !$OMP END PARALLEL DO

    end do ! levIndex

  END SUBROUTINE ben_getEnsMean

  !--------------------------------------------------------------------------
  ! EnsembleScaleDecomposition
  !--------------------------------------------------------------------------
  SUBROUTINE EnsembleScaleDecomposition()
    implicit none

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

    integer :: gstFilterID, mIndex, nIndex, mymBeg, mymEnd, mynBeg, mynEnd, mymSkip, mynSkip
    integer :: mymCount, mynCount, ierr
    integer :: myMemberBeg, myMemberEnd, myMemberCount, maxMyMemberCount, nEnsOverDimension

    type(struct_lst)    :: lst_ben_filter ! Spectral transform Parameters for filtering

    character(len=19)   :: kind

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
    if ( mpi_myid == 0 ) then
      write(*,*)
      write(*,*) 'Scale decomposition of the ensemble perturbations'
      write(*,*) '   number of WaveBands = ', nWaveBand
      write(*,*) '   WaveBand Peaks (total wavenumber)...'
      do waveBandIndex = 1, nWaveBand
        write(*,*) waveBandIndex, waveBandPeaks(waveBandIndex)
      end do
    end if

    !
    !- Setup a spectral transform for filtering (nk = nEnsOverDimension)
    !

    call mpivar_setup_levels_npex(nEns,                                  & ! IN
                                  myMemberBeg,myMemberEnd,myMemberCount)   ! OUT
    call rpn_comm_allreduce(myMemberCount, maxMyMemberCount, &
                            1,"MPI_INTEGER","mpi_max","GRID",ierr)
    nEnsOverDimension  = mpi_npex * maxMyMemberCount

    if (hco_ens%global) then
      ! Global mode
      gstFilterID = gst_setup(ni,nj,ntrunc,nEnsOverDimension)
      if (mpi_myid == 0) write(*,*) 'ben : returned value of gstFilterID = ',gstFilterID

      nla_filter = gst_getNla(gstFilterID)
      nphase_filter = 2

      allocate(nIndex_vec(nla_filter))
      call mpivar_setup_m(ntrunc,mymBeg,mymEnd,mymSkip,mymCount)
      call mpivar_setup_n(ntrunc,mynBeg,mynEnd,mynSkip,mynCount)
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
      call lst_Setup( lst_ben_filter,                   & ! OUT
           ni, nj, hco_ens%dlon, ntrunc,                & ! IN
           'LatLonMN', maxlevels_opt=nEnsOverDimension, & ! IN
           gridDataOrder_opt='kij' )                      ! IN

      nla_filter = lst_ben_filter%nla
      nphase_filter = lst_ben_filter%nphase
    end if

    !
    !- 1.  Scale decomposition for every wave band except for wave band #1
    !
    allocate(ResponseFunction(nla_filter,2:nWaveBand))
    allocate(ensPertSP(nla_filter,nphase_filter,nEnsOverDimension))
    allocate(ensPertSPfiltered(nla_filter,nphase_filter,nEnsOverDimension))
    allocate(ensPertGD(nEnsOverDimension,myLonBeg:myLonEnd,myLatBeg:myLatEnd))

    ensPertSP        (:,:,:) = 0.0d0
    ensPertSPfiltered(:,:,:) = 0.0d0

    !- 1.1 Pre-compute the response function for each wave band except for wave band #1
    do waveBandIndex = nWaveBand, 2, -1 ! Start with the largest scales
      do ila_filter = 1, nla_filter
        if (hco_ens%global) then
          totwvnb_r8 = real(nIndex_vec(ila_filter),8)
        else
          totwvnb_r8 = lst_ben_filter%k_r8(ila_filter)
        end if
        ResponseFunction(ila_filter,waveBandIndex) = spf_FilterResponseFunction(totwvnb_r8,waveBandIndex, waveBandPeaks, nWaveBand)
        write(*,*) totwvnb_r8, ResponseFunction(ila_filter,waveBandIndex)
      end do
    end do
    if (hco_ens%global) deallocate(nIndex_vec)

    do stepIndex = 1, numStep ! Loop on ensemble time bin
      do levIndex = 1, ens_getNumK(ensPerts(waveBandIndex)) ! Loop on variables and vertical levels
        ptr4d_r4 => ens_getOneLev_r4(ensPerts(1),levIndex)

        !- 1.2 GridPoint space -> Spectral Space
 !$OMP PARALLEL DO PRIVATE (latIndex)
        do latIndex = myLatBeg, myLatEnd
          ensPertGD(:,:,latIndex) = 0.0d0
        end do
 !$OMP END PARALLEL DO
 !$OMP PARALLEL DO PRIVATE (memberIndex,latIndex,lonIndex)
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            do memberIndex = 1, nEns
              ensPertGD(memberIndex,lonIndex,latIndex) = dble(ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex))
            end do
          end do
        end do
 !$OMP END PARALLEL DO
        if (hco_ens%global) then
          ! Global Mode
          call gst_setID(gstFilterID) ! IN
          call gst_reespe_kij(ensPertSP, & ! OUT
                              ensPertGD)   ! IN
        else
          ! LAM mode
          kind = 'GridPointToSpectral'
          call lst_VarTransform( lst_ben_filter%id,      & ! IN
               ensPertSP,              & ! OUT
               ensPertGD,              & ! IN 
               kind, nEnsOverDimension ) ! IN
        end if

        !- 1.3 Filtering and transformation back to grid point space 
        do waveBandIndex = nWaveBand, 2, -1 ! Start with the largest scales
          ! Filtering
!$OMP PARALLEL DO PRIVATE (memberIndex,p,ila_filter)
          do memberIndex = 1, nEns
            do p = 1, nphase_filter
              do ila_filter = 1, nla_filter
                ensPertSPfiltered(ila_filter,p,memberIndex) = &
                     ensPertSP(ila_filter,p,memberIndex) * ResponseFunction(ila_filter,waveBandIndex)
              end do
            end do
          end do
 !$OMP END PARALLEL DO

          ! Spectral Space -> GridPoint space
          if (hco_ens%global) then
            ! Global Mode
            call gst_setID(gstFilterID) ! IN
            call gst_speree_kij(ensPertSPfiltered, & ! IN
                                ensPertGD)           ! OUT
          else
            ! LAM mode
            kind = 'SpectralToGridPoint'
            call lst_VarTransform( lst_ben_filter%id,      & ! IN
                                   ensPertSPfiltered,      & ! IN
                                   ensPertGD,              & ! OUT
                                   kind, nEnsOverDimension ) ! IN
          end if
          ptr4d_r4 => ens_getOneLev_r4(ensPerts(waveBandIndex),levIndex)
!$OMP PARALLEL DO PRIVATE (memberIndex,latIndex,lonIndex)
          do latIndex = myLatBeg, myLatEnd
            do lonIndex = myLonBeg, myLonEnd
              do memberIndex = 1, nEns
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
    !- 2.  Isolate the smallest scales in waveBandIndex = 1 by difference in grid point space
    !
    allocate(bandSum(myLonBeg:myLonEnd,myLatBeg:myLatEnd))
    do stepIndex = 1, numStep
!$OMP PARALLEL DO PRIVATE (memberIndex,levIndex,latIndex,lonIndex,waveBandIndex,bandsum,ptr4d_r4)
      do levIndex = 1, ens_getNumK(ensPerts(1))
        do memberIndex = 1, nEns
          bandSum(:,:) = 0.d0
          do waveBandIndex = 2, nWaveBand
            ptr4d_r4 => ens_getOneLev_r4(ensPerts(waveBandIndex),levIndex)
            do latIndex = myLatBeg, myLatEnd
              do lonIndex = myLonBeg, myLonEnd
                bandSum(lonIndex,latIndex) = bandSum(lonIndex,latIndex) + dble(ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex))
              end do
            end do
          end do
          ptr4d_r4 => ens_getOneLev_r4(ensPerts(1),levIndex)
          do latIndex = myLatBeg, myLatEnd
            do lonIndex = myLonBeg, myLonEnd
              ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex) = sngl(dble(ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex)) - bandSum(lonIndex,latIndex))
            end do
          end do
        end do
      end do
!$OMP END PARALLEL DO
    end do
    deallocate(bandSum)

  END SUBROUTINE EnsembleScaleDecomposition

!--------------------------------------------------------------------------
! ben_reduceToMPILocal
!--------------------------------------------------------------------------
  SUBROUTINE ben_reduceToMPILocal(cv_mpilocal,cv_mpiglobal)
    implicit none
    real(8), intent(out) :: cv_mpilocal(cvDim_mpilocal)
    real(8), intent(in)  :: cv_mpiglobal(:)

    if (verbose) write(*,*) 'Entering ben_reduceToMPILocal'

    call loc_reduceToMPILocal(locIDs(1),cv_mpilocal,cv_mpiglobal)

 END SUBROUTINE ben_reduceToMPILocal

!--------------------------------------------------------------------------
! ben_reduceToMPILocal_r4
!--------------------------------------------------------------------------
  SUBROUTINE ben_reduceToMPILocal_r4(cv_mpilocal,cv_mpiglobal)
    implicit none
    real(4), intent(out) :: cv_mpilocal(cvDim_mpilocal)
    real(4), intent(in)  :: cv_mpiglobal(:)

    if (verbose) write(*,*) 'Entering reduceToMPILocal_r4'

    call loc_reduceToMPILocal_r4(locIDs(1),cv_mpilocal,cv_mpiglobal) ! IN

 END SUBROUTINE ben_reduceToMPILocal_r4

!--------------------------------------------------------------------------
! ben_expandToMPIGlobal
!--------------------------------------------------------------------------
  SUBROUTINE ben_expandToMPIGlobal(cv_mpilocal,cv_mpiglobal)
    implicit none

    real(8), intent(in)  :: cv_mpilocal(cvDim_mpilocal)
    real(8), intent(out) :: cv_mpiglobal(:)

    if (verbose) write(*,*) 'Entering ben_expandToMPIGlobal'

    call loc_expandToMPIGlobal(locIDs(1), cv_mpilocal,  & ! IN
                               cv_mpiglobal)              ! OUT  

  end SUBROUTINE ben_expandToMPIGlobal

!--------------------------------------------------------------------------
! ben_expandToMPIGlobal_r4
!--------------------------------------------------------------------------
  SUBROUTINE ben_expandToMPIGlobal_r4(cv_mpilocal,cv_mpiglobal)
    implicit none

    real(4), intent(in)  :: cv_mpilocal(cvDim_mpilocal)
    real(4), intent(out) :: cv_mpiglobal(:)

    if (verbose) write(*,*) 'Entering ben_expandToMPIGlobal_r4'

    call loc_expandToMPIGlobal_r4(locIDs(1), cv_mpilocal,  & ! IN
                                  cv_mpiglobal)              ! OUT

  end SUBROUTINE ben_expandToMPIGlobal_r4

!--------------------------------------------------------------------------
! ben_BSqrt
!--------------------------------------------------------------------------
  SUBROUTINE ben_BSqrt(controlVector_in,statevector,useFSOFcst_opt)
    implicit none

    real(8)          :: controlVector_in(cvDim_mpilocal) 
    type(struct_gsv) :: statevector
    logical,optional :: useFSOFcst_opt

    type(struct_ens)    :: ensAmplitude_M
    
    integer   :: ierr, levIndex, latIndex, memberIndex, waveBandIndex
    integer   :: numStepAmplitude, amp3dStepIndex
    logical   :: immediateReturn
    logical   :: useFSOFcst

    call tmg_start(67,'BEN_BARR')
    if (mpi_doBarrier) call rpn_comm_barrier('GRID',ierr)
    call tmg_stop(67)

    !
    !- 1.  Tests
    !
    if (.not. initialized) then
      if (mpi_myid == 0) write(*,*) 'ben_bsqrt: bMatrixEnsemble not initialized'
      return
    end if

    if (sum(scaleFactor) == 0.0d0) then
      if (mpi_myid == 0) write(*,*) 'ben_bsqrt: scaleFactor=0, skipping bSqrt'
      return
    end if

    ! only check controlVector on proc 0, since may be zero-length on some procs
    if (mpi_myid == 0) then
      immediateReturn = .false.
      if (maxval(controlVector_in) == 0.0d0 .and. minval(controlVector_in) == 0.0d0) then
        write(*,*) 'ben_bsqrt: controlVector=0, skipping bSqrt'
        immediateReturn = .true.
      end if
    end if
    call rpn_comm_bcast(immediateReturn, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)
    if (immediateReturn) return

    if (mpi_myid == 0) write(*,*) 'ben_bsqrt: starting'
    if (mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if (present(useFSOFcst_opt)) then
      useFSOFcst = useFSOFcst_opt
    else
      useFSOFcst = .false.
    end if

    !
    !- 2.  Compute the analysis increment from Bens
    !
    if (verbose) write(*,*) 'ben_bsqrt: allocating ensAmplitude_M'
    if (useFSOFcst) then
      numStepAmplitude = numStepAmplitudeFSOFcst
      amp3dStepIndex   = amp3dStepIndexFSOFcst
    else
      numStepAmplitude =  numStepAmplitudeAssimWindow
      amp3dStepIndex   = amp3dStepIndexAssimWindow
    end if
    call ens_allocate(ensAmplitude_M, nEnsOverDimension, numStepAmplitude, hco_ens, vco_ens, dateStampList, &
                      varNames_opt=varNameALFA, dataKind_opt=8)

    do waveBandIndex = 1, nWaveBand !  Loop on WaveBand (for ScaleDependent Localization)

      ! 2.1 Compute the ensemble amplitudes
      call loc_Lsqrt( locIDs(waveBandIndex),controlVector_in, & ! IN
                      ensAmplitude_M,                         & ! OUT
                      amp3dStepIndex)                           ! IN

      ! 2.2 Advect the amplitudes
      if      (advectAmplitudeFSOFcst   .and. useFSOFcst) then
        call tmg_start(131,'BEN_ADVEC_AMP_FSO_TL')
        call adv_ensemble_tl( ensAmplitude_M,              & ! INOUT
                              adv_amplitudeFSOFcst, nEns )   ! IN
        call tmg_stop(131)
      else if (advectAmplitudeAssimWindow .and. .not. useFSOFcst) then
        call tmg_start(133,'BEN_ADVEC_AMP_ANL_TL')
        call adv_ensemble_tl( ensAmplitude_M,              & ! INOUT
                              adv_amplitudeAssimWindow, nEns ) ! IN
        call tmg_stop(133)
      end if

      if ( keepAmplitude .and. waveBandIndex == 1 ) call ens_copy(ensAmplitude_M, ensAmplitudeStorage)

      ! 2.3 Compute increment by multiplying amplitudes by member perturbations
      call addEnsMember( ensAmplitude_M, statevector,  & ! INOUT 
                         waveBandIndex, useFSOFcst )     ! IN

    end do ! Loop on WaveBand

    call ens_deallocate(ensAmplitude_M)

    ! 2.4 Advect Increments
    if ( advectEnsPertAnlInc ) then
      call tmg_start(138,'BEN_ADVEC_ANLINC_TL')
      call adv_statevector_tl( statevector,  & ! INOUT
                               adv_analInc )   ! IN
      call tmg_stop(138)
    end if

    !
    !- 3.  Variable transforms
    !
    if ( ctrlVarHumidity == 'LQ') then
       call vtr_transform( statevector, & ! INOUT
                           'LQtoHU_tlm' ) ! IN
    end if

    if (mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if (mpi_myid == 0) write(*,*) 'ben_bsqrt: done'

  END SUBROUTINE ben_BSqrt

!--------------------------------------------------------------------------
! ben_BSqrtAd
!--------------------------------------------------------------------------
  SUBROUTINE ben_BSqrtAd(statevector,controlVector_out,useFSOFcst_opt)
    implicit none

    real(8)           :: controlVector_out(cvDim_mpilocal) 
    type(struct_gsv)  :: statevector
    logical, optional :: useFSOFcst_opt

    type(struct_ens)   :: ensAmplitude_M

    integer           :: ierr, levIndex, latIndex, memberIndex, waveBandIndex
    integer           :: numStepAmplitude,amp3dStepIndex
    logical           :: useFSOFcst

    !
    !- 1.  Tests
    !
    call tmg_start(67,'BEN_BARR')
    if (mpi_doBarrier) call rpn_comm_barrier('GRID',ierr)
    call tmg_stop(67)

    if (.not. initialized) then
      if (mpi_myid == 0) write(*,*) 'ben_bsqrtad: bMatrixEnsemble not initialized'
      return
    end if

    if (sum(scaleFactor) == 0.0d0) then
      if (mpi_myid == 0) write(*,*) 'ben_bsqrtad: scaleFactor=0, skipping bSqrtAd'
      return
    end if

    if (mpi_myid == 0) write(*,*) 'ben_bsqrtad: starting'
    if (mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if (present(useFSOFcst_opt)) then
      useFSOFcst = useFSOFcst_opt
    else
      useFSOFcst = .false.
    end if

    !
    !- 3.  Variable transforms
    !
    if ( ctrlVarHumidity == 'LQ') then
      call vtr_transform( statevector, & ! INOUT
                          'LQtoHU_tlm' ) ! IN
    end if

    !
    !- 2.  Compute the analysis increment from Bens
    !

    ! 2.4 Advect Increments
    if ( advectEnsPertAnlInc ) then
      call tmg_start(139,'BEN_ADVEC_ANLINC_AD')
      call adv_statevector_ad( statevector,  & ! INOUT
                               adv_analInc )   ! IN
      call tmg_stop(139)
    end if

    if (verbose) write(*,*) 'ben_bsqrtAd: allocating ensAmplitude_M'
    if (useFSOFcst) then
      numStepAmplitude = numStepAmplitudeFSOFcst
      amp3dStepIndex   = amp3dStepIndexFSOFcst
    else
      numStepAmplitude = numStepAmplitudeAssimWindow
      amp3dStepIndex   = amp3dStepIndexAssimWindow
    end if
    call ens_allocate(ensAmplitude_M, nEnsOverDimension, numStepAmplitude, hco_ens, vco_ens, dateStampList, &
                      varNames_opt=varNameALFA, dataKind_opt=8)

    do waveBandIndex = 1, nWaveBand !  Loop on WaveBand (for ScaleDependent Localization)

      ! 2.3 Compute increment by multiplying amplitudes by member perturbations
      call addEnsMemberAd( statevector, ensAmplitude_M,  & ! INOUT
                           waveBandIndex, useFSOFcst)      ! IN
      ! 2.2 Advect the  amplitudes
      if      (advectAmplitudeFSOFcst   .and. useFSOFcst) then
        call tmg_start(132,'BEN_ADVEC_AMP_FSO_AD')
        call adv_ensemble_ad( ensAmplitude_M,              & ! INOUT
                              adv_amplitudeFSOFcst, nEns )   ! IN
        call tmg_stop(132)
      else if (advectAmplitudeAssimWindow .and. .not. useFSOFcst) then
        call tmg_start(134,'BEN_ADVEC_AMP_ANL_AD')
        call adv_ensemble_ad( ensAmplitude_M,              & ! INOUT
                              adv_amplitudeAssimWindow, nEns ) ! IN
        call tmg_stop(134)
      end if

      ! 2.1 Compute the ensemble amplitudes
      call loc_LsqrtAd( locIDs(waveBandIndex),ensAmplitude_M, & ! IN
                        controlVector_out,                    & ! OUT
                        amp3dStepIndex)                         ! IN

    end do ! Loop on WaveBand

    call ens_deallocate(ensAmplitude_M)

    if (mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if (mpi_myid == 0) write(*,*) 'ben_bsqrtAd: done'

  END SUBROUTINE ben_BSqrtAd

!--------------------------------------------------------------------------
! addEnsMember
!--------------------------------------------------------------------------
  SUBROUTINE addEnsMember(ensAmplitude_M, statevector_out, &
                          waveBandIndex, useFSOFcst_opt)
    implicit none

    type(struct_ens)    :: ensAmplitude_M
    type(struct_gsv)    :: statevector_out
    integer, intent(in) :: waveBandIndex
    logical, optional   :: useFSOFcst_opt

    real(8), pointer    :: ensAmplitude_M_oneLev(:,:,:,:), ensAmplitude_M_oneLevM1(:,:,:,:)
    real(8), allocatable, target :: ensAmplitude_MT(:,:,:,:)
    real(8), pointer     :: ensAmplitude_MT_ptr(:,:,:,:)
    real(8), pointer     :: increment_out(:,:,:,:)
    real(8), allocatable :: increment_out2(:,:,:)
    real(4), pointer     :: ensMemberAll_r4(:,:,:,:)
    integer     :: lev, lev2, levIndex, stepIndex, stepIndex_amp, latIndex, lonIndex, topLevOffset, numLev, memberIndex
    character(len=4)     :: varName

    logical             :: useFSOFcst
    integer             :: stepIndex2, stepBeg, stepEnd, numStepAmplitude

    
    if (verbose) write(*,*) 'Entering ben_addEnsMember'

    if (vco_anl%Vcode /= 5002 .and. (vco_anl%nlev_T > 1 .or. vco_anl%nlev_M > 1) ) then
      call utl_abort('addEnsMemberAd: Only 5002 supported in 3D mode for now!')
    end if

    call tmg_start(62,'ADDMEM')

    if (present(useFSOFcst_opt)) then
      useFSOFcst = useFSOFcst_opt
    else
      useFSOFcst = .false.
    end if
    if (useFSOFcst .and. fsoLeadTime > 0.0d0) then
      stepBeg = numStep
      stepEnd = stepBeg
      if (mpi_myid == 0) write(*,*) 'ben_bsqrtad: using forecast ensemble stored at timestep ',stepEnd
    else
      stepBeg = 1
      stepEnd = numStepAssimWindow
    end if
    if (useFSOFcst) then
      numStepAmplitude =  numStepAmplitudeFSOFcst
    else
      numStepAmplitude =  numStepAmplitudeAssimWindow
    end if
    
    allocate(ensAmplitude_MT(nEns,numStepAmplitude,myLonBeg:myLonEnd,myLatBeg:myLatEnd))
    allocate(increment_out2(numStep,myLonBeg:myLonEnd,myLatBeg:myLatEnd))

    do levIndex = 1, ens_getNumK(ensPerts(waveBandIndex))

      lev = ens_getLevFromK(ensPerts(1),levIndex)
      varName = ens_getVarNameFromK(ensPerts(1),levIndex)

!$OMP PARALLEL DO PRIVATE (latIndex)
      do latIndex = myLatBeg, myLatEnd
        increment_out2(:,:,latIndex) = 0.0d0
      end do
!$OMP END PARALLEL DO

      call tmg_start(66,'ADDMEM_PREPAMP')

      if (vnl_varLevelFromVarname(varName) == 'MM') then

        ensAmplitude_M_oneLev => ens_getOneLev_r8(ensAmplitude_M,lev)
        ensAmplitude_MT_ptr(1:,1:,myLonBeg:,myLatBeg:) => ensAmplitude_M_oneLev(1:nEns,:,:,:)

      else if (vnl_varLevelFromVarname(varName) == 'TH') then

        if (lev == 1) then
          ! use top momentum level amplitudes for top thermo level
          ensAmplitude_M_oneLev => ens_getOneLev_r8(ensAmplitude_M,lev)
          ensAmplitude_MT_ptr(1:,1:,myLonBeg:,myLatBeg:) => ensAmplitude_M_oneLev(1:nEns,:,:,:)
        else if (lev == nLevEns_T) then
          ! use surface momentum level amplitudes for surface thermo level
          ensAmplitude_M_oneLev => ens_getOneLev_r8(ensAmplitude_M,nLevEns_M)
          ensAmplitude_MT_ptr(1:,1:,myLonBeg:,myLatBeg:) => ensAmplitude_M_oneLev(1:nEns,:,:,:)
        else
          ! for other levels, interpolate momentum weights to get thermo amplitudes
!$OMP PARALLEL DO PRIVATE (latIndex, ensAmplitude_M_oneLev, ensAmplitude_M_oneLevM1)
          do latIndex = myLatBeg, myLatEnd
            ensAmplitude_M_oneLev   => ens_getOneLev_r8(ensAmplitude_M,lev)
            ensAmplitude_M_oneLevM1 => ens_getOneLev_r8(ensAmplitude_M,lev-1)
            ensAmplitude_MT(:,:,:,latIndex) = 0.5d0*( ensAmplitude_M_oneLevM1(1:nEns,:,:,latIndex) +   &
                                                   ensAmplitude_M_oneLev(1:nEns,:,:,latIndex) )
          end do
!$OMP END PARALLEL DO
          ensAmplitude_MT_ptr(1:,1:,myLonBeg:,myLatBeg:) => ensAmplitude_MT(:,:,:,:)
        end if

      else if (vnl_varLevelFromVarname(varName) == 'SF') then

        ! surface variable
        ensAmplitude_M_oneLev   => ens_getOneLev_r8(ensAmplitude_M,nLevEns_M)
        ensAmplitude_MT_ptr(1:,1:,myLonBeg:,myLatBeg:) => ensAmplitude_M_oneLev(1:nEns,:,:,:)
      end if
      call tmg_stop(66)

      call tmg_start(77,'ADDMEM_INNER')

      ensMemberAll_r4 => ens_getOneLev_r4(ensPerts(waveBandIndex),levIndex)
!$OMP PARALLEL DO PRIVATE (latIndex,lonIndex,stepIndex,stepIndex2,stepIndex_amp,memberIndex)
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          do stepIndex = StepBeg, StepEnd
            stepIndex2 = stepIndex - stepBeg + 1
            if      (advectAmplitudeFSOFcst   .and. useFSOFcst) then
              stepIndex_amp = 2
            else if (advectAmplitudeAssimWindow .and. .not. useFSOFcst) then
              stepIndex_amp = stepIndex
            else
              stepIndex_amp = 1
            end if
            do memberIndex = 1, nEns
              increment_out2(stepIndex2,lonIndex,latIndex) = increment_out2(stepIndex2,lonIndex,latIndex) +   &
                ensAmplitude_MT_ptr(memberIndex,stepIndex_amp,lonIndex,latIndex) *  &
                dble(ensMemberAll_r4(memberIndex,stepIndex,lonIndex,latIndex))
            end do ! memberIndex
          end do ! stepIndex
        end do ! lonIndex
      end do ! latIndex
!$OMP END PARALLEL DO

      call tmg_stop(77)

      ! compute increment level from amplitude/member level
      if (vnl_varLevelFromVarname(varName) == 'SF') then
        topLevOffset = 1
      else if (vnl_varLevelFromVarname(varName) == 'MM') then
        topLevOffset = topLevIndex_M
      else
        topLevOffset = topLevIndex_T
      end if
      lev2 = lev - 1 + topLevOffset

      increment_out => gsv_getField_r8(statevector_out, varName)
!$OMP PARALLEL DO PRIVATE (stepIndex, stepIndex2)
      do stepIndex = StepBeg, StepEnd
        stepIndex2 = stepIndex - StepBeg + 1
        increment_out(:,:,lev2,stepIndex2) = increment_out(:,:,lev2,stepIndex2) + increment_out2(stepIndex2,:,:)
      end do
!$OMP END PARALLEL DO

    end do ! levIndex

    deallocate(ensAmplitude_MT)
    deallocate(increment_out2)

    call tmg_stop(62)

  END SUBROUTINE addEnsMember

!--------------------------------------------------------------------------
! addEnsMemberAd
!--------------------------------------------------------------------------
  SUBROUTINE addEnsMemberAd(statevector_in, ensAmplitude_M, &
                                   waveBandIndex, useFSOFcst_opt)
    implicit none

    type(struct_ens)    :: ensAmplitude_M
    type(struct_gsv)   :: statevector_in
    integer,intent(in) :: waveBandIndex
    logical,optional   :: useFSOFcst_opt

    real(8), pointer    :: ensAmplitude_M_oneLev(:,:,:,:), ensAmplitude_M_oneLevM1(:,:,:,:)
    real(8), allocatable :: ensAmplitude_MT(:,:)
    real(8), pointer :: increment_in(:,:,:,:)
    real(8), allocatable :: increment_in2(:,:,:)
    real(4), pointer :: ensMemberAll_r4(:,:,:,:)
    integer          :: levIndex, lev, lev2, stepIndex, stepIndex_amp, latIndex, lonIndex, topLevOffset, numLev, memberIndex
    character(len=4)     :: varName
    integer     ::  stepBeg, stepEnd, stepIndex2, numStepAmplitude
    logical          :: useFSOFcst

    if (verbose) write(*,*) 'Entering ben_addEnsMemberAd'

    if (vco_anl%Vcode /= 5002 .and. (vco_anl%nlev_T > 1 .or. vco_anl%nlev_M > 1) ) then
      call utl_abort('addEnsMemberAd: Only 5002 supported in 3D mode for now!')
    end if

    call tmg_start(63,'ADDMEMAD')

    if (present(useFSOFcst_opt)) then
      useFSOFcst = useFSOFcst_opt
    else
      useFSOFcst = .false.
    end if

    if (useFSOFcst .and. fsoLeadTime > 0.0d0) then
      stepBeg = numStep
      stepEnd = stepBeg
      if (mpi_myid == 0) write(*,*) 'ben_bsqrtad: using forecast ensemble stored at timestep ',stepEnd
    else
      stepBeg = 1
      stepEnd = numStepAssimWindow
    end if
    if (useFSOFcst) then
      numStepAmplitude =  numStepAmplitudeFSOFcst
    else
      numStepAmplitude =  numStepAmplitudeAssimWindow
    end if

    allocate(ensAmplitude_MT(nEns,numStepAmplitude))
    allocate(increment_in2(numStep,myLonBeg:myLonEnd,myLatBeg:myLatEnd))

    ! set output ensemble Amplitude to zero
    call tmg_start(69,'ADDMEMAD_ZERO')
!$OMP PARALLEL DO PRIVATE (levIndex, ensAmplitude_M_oneLev)
    do levIndex = 1, nLevEns_M
      ensAmplitude_M_oneLev => ens_getOneLev_r8(ensAmplitude_M,levIndex)
      ensAmplitude_M_oneLev(:,:,:,:) = 0.0d0
    end do
!$OMP END PARALLEL DO
    call tmg_stop(69)

    do levIndex = 1, ens_getNumK(ensPerts(waveBandIndex))

      lev = ens_getLevFromK(ensPerts(1),levIndex)
      varName = ens_getVarNameFromK(ensPerts(1),levIndex)

      ! compute increment level from amplitude/member level
      if (vnl_varLevelFromVarname(varName) == 'SF') then
        topLevOffset = 1
      else if (vnl_varLevelFromVarname(varName) == 'MM') then
        topLevOffset = topLevIndex_M
      else
        topLevOffset = topLevIndex_T
      end if
      lev2 = lev - 1 + topLevOffset

      call tmg_start(65,'ADDMEMAD_SHUFFLE')
      increment_in => gsv_getField_r8(statevector_in, varName)
!$OMP PARALLEL DO PRIVATE (stepIndex, stepIndex2)
      do stepIndex = stepBeg, stepEnd
        stepIndex2 = stepIndex - stepBeg + 1
        increment_in2(stepIndex2,:,:) = increment_in(:,:,lev2,stepIndex2)
      end do
!$OMP END PARALLEL DO
      call tmg_stop(65)

      !ensAmpZeroed(:,:,:) = .false.
      ensMemberAll_r4 => ens_getOneLev_r4(ensPerts(waveBandIndex),levIndex)
!$OMP PARALLEL DO PRIVATE (latIndex,lonIndex,stepIndex, stepIndex2, stepIndex_amp,memberIndex,ensAmplitude_M_oneLev, ensAmplitude_M_oneLevM1, ensAmplitude_MT)
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd

          if (omp_get_thread_num() == 0) call tmg_start(78,'ADDMEMAD_INNER')
          ensAmplitude_MT(:,:) = 0.0d0
          do stepIndex = StepBeg, StepEnd
            stepIndex2 = stepIndex-StepBeg+1
            if      (advectAmplitudeFSOFcst   .and. useFSOFcst) then
              stepIndex_amp = 2
            else if (advectAmplitudeAssimWindow .and. .not. useFSOFcst) then
              stepIndex_amp = stepIndex
            else
              stepIndex_amp = 1
            end if
            do memberIndex = 1, nEns
              ensAmplitude_MT(memberIndex,stepIndex_amp) = ensAmplitude_MT(memberIndex,stepIndex_amp) +  &
                increment_in2(stepIndex2,lonIndex,latIndex) * dble(ensMemberAll_r4(memberIndex,stepIndex,lonIndex,latIndex))
            end do ! memberIndex
          end do ! stepIndex
          if (omp_get_thread_num() == 0) call tmg_stop(78)

          ! transform thermo/momentum level amplitude sensitivites appropriately

          if (omp_get_thread_num() == 0) call tmg_start(68,'ADDMEMAD_PREPAMP')
          if (vnl_varLevelFromVarname(varName) == 'MM') then

            ensAmplitude_M_oneLev => ens_getOneLev_r8(ensAmplitude_M,lev)
            ensAmplitude_M_oneLev(1:nEns,:,lonIndex,latIndex) = &
                 ensAmplitude_M_oneLev(1:nEns,:,lonIndex,latIndex) + ensAmplitude_MT(:,:)

          else if (vnl_varLevelFromVarname(varName) == 'TH') then

            if (lev == 1) then
              ! use top momentum level amplitudes for top thermo level
              ensAmplitude_M_oneLev => ens_getOneLev_r8(ensAmplitude_M,lev)
              ensAmplitude_M_oneLev(1:nEns,:,lonIndex,latIndex) = &
                 ensAmplitude_M_oneLev(1:nEns,:,lonIndex,latIndex) + ensAmplitude_MT(:,:)
            else if (lev == nLevEns_T) then
              ! use surface momentum level amplitudes for surface thermo level
              ensAmplitude_M_oneLev => ens_getOneLev_r8(ensAmplitude_M,nLevEns_M)
              ensAmplitude_M_oneLev(1:nEns,:,lonIndex,latIndex) = &
                   ensAmplitude_M_oneLev(1:nEns,:,lonIndex,latIndex) + ensAmplitude_MT(:,:)
            else
              ! for other levels, interpolate momentum weights to get thermo amplitudes
              ensAmplitude_M_oneLev   => ens_getOneLev_r8(ensAmplitude_M,lev)
              ensAmplitude_M_oneLevM1 => ens_getOneLev_r8(ensAmplitude_M,lev-1)
              ensAmplitude_M_oneLev(1:nEns,:,lonIndex,latIndex)   = &
                   ensAmplitude_M_oneLev(1:nEns,:,lonIndex,latIndex)   + 0.5d0*ensAmplitude_MT(:,:)
              ensAmplitude_M_oneLevM1(1:nEns,:,lonIndex,latIndex) = &
                   ensAmplitude_M_oneLevM1(1:nEns,:,lonIndex,latIndex) + 0.5d0*ensAmplitude_MT(:,:)
            end if

          else if (vnl_varLevelFromVarname(varName) == 'SF') then

            ! surface variable
            ensAmplitude_M_oneLev   => ens_getOneLev_r8(ensAmplitude_M,nLevEns_M)
            ensAmplitude_M_oneLev (1:nEns,:,lonIndex,latIndex) = &
                 ensAmplitude_M_oneLev(1:nEns,:,lonIndex,latIndex) + ensAmplitude_MT(:,:)

          end if
          if (omp_get_thread_num() == 0) call tmg_stop(68)

        end do ! lonIndex
      end do ! latIndex
!$OMP END PARALLEL DO

    end do ! levIndex

    deallocate(ensAmplitude_MT)
    deallocate(increment_in2)

    call tmg_stop(63)

  END SUBROUTINE addEnsMemberAd

!--------------------------------------------------------------------------
! EnsembleDiagnostic
!--------------------------------------------------------------------------
  SUBROUTINE EnsembleDiagnostic(mode)
    implicit none

    character(len=*), intent(in) :: mode

    type(struct_gsv) :: statevector, statevector_temp

    integer :: nWaveBandToDiagnose, waveBandIndex, memberIndex

    real(8) :: dnens2

    character(len=12):: etiket, modeEtiket
    character(len=2) :: wbnum

    if ( trim(mode) == 'FullPerturbations') then
       nWaveBandToDiagnose = 1
    else if ( trim(mode) == 'WaveBandPerturbations' ) then
       nWaveBandToDiagnose = nWaveBand
    else
       write(*,*)
       write(*,*) 'mode = ', trim(mode)
       call utl_abort('EnsembleDiagnostic: unknown mode')
    end if

    if ( mpi_myid == 0 ) write(*,*)
    if ( mpi_myid == 0 ) write(*,*) 'EnsembleDiagnostic in mode: ', mode

    !
    !- Write each wave band for a selected member
    !
    if (trim(LocalizationType) == 'ScaleDependent') then
       if ( mpi_myid == 0 ) write(*,*) '   writing perturbations for member 001'
       memberIndex = 1
       dnens2 = sqrt(1.0d0*dble(nEns-1))
       do waveBandIndex = 1, nWaveBandToDiagnose
          if ( mpi_myid == 0 ) write(*,*) '     waveBandIndex = ', waveBandIndex
          call gsv_allocate(statevector, tim_nstepobsinc, hco_ens, vco_anl, &
                            datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true.)
          call ben_getPerturbation( statevector,    & ! OUT
                                    memberIndex,    & ! IN
                                    'ConstantValue', waveBandIndex ) ! IN
          if ( trim(mode) == 'FullPerturbations') then
             etiket = 'PERT001_FULL'
          else
             write(wbnum,'(I2.2)') waveBandIndex
             etiket = 'PERT001_WB' // trim(wbnum)
          end if
          call gsv_writeToFile(statevector,'./ens_pert001.fst',etiket, & ! IN
                               dnens2,HUcontainsLQ_opt=HUcontainsLQ_gsv )    ! IN
          call gsv_deallocate(statevector)
       end do
    end if

    !
    !- Compute the standard deviations for each wave band
    !
    if ( mpi_myid == 0 ) write(*,*) '   computing Std.Dev.'
    call gsv_allocate(statevector_temp, tim_nstepobsinc, hco_ens, vco_anl, &
         mpi_local_opt=.true.)

    do waveBandIndex = 1, nWaveBandToDiagnose
       if ( mpi_myid == 0 ) write(*,*) '     waveBandIndex = ', waveBandIndex
       call gsv_allocate(statevector, tim_nstepobsinc, hco_ens, vco_anl, &
                         datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true.)
       call gsv_zero(statevector)
       do memberIndex = 1, nEns
          !- Get normalized perturbations
          call ens_copyMember(ensPerts(waveBandIndex), & ! IN
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
          write(wbnum,'(I2.2)') waveBandIndex
          etiket = 'STDDEV_WB' // trim(wbnum)
       end if
       call gsv_writeToFile(statevector,'./ens_stddev.fst',etiket, & ! IN
                            HUcontainsLQ_opt=HUcontainsLQ_gsv)       ! IN
       call gsv_deallocate(statevector)
    end do

    call gsv_deallocate(statevector_temp)

  END SUBROUTINE EnsembleDiagnostic

!--------------------------------------------------------------------------
! ben_writeAmplitude
!--------------------------------------------------------------------------
  SUBROUTINE ben_writeAmplitude(ensPathName, ensFileNamePrefix, ip3)
    implicit none

    character(len=*), intent(in) :: ensPathName
    character(len=*), intent(in) :: ensFileNamePrefix
    integer,          intent(in) :: ip3

    if (initialized .and. keepAmplitude) then
      if ( mpi_myid == 0 ) write(*,*)
      if ( mpi_myid == 0 ) write(*,*) 'bmatrixEnsemble_mod: Writing the amplitude field'
      call ens_writeEnsemble(ensAmplitudeStorage, ensPathName, ensFileNamePrefix, &
                             'LQ', 'FROM_BENS', 'R',varNames_opt=varNameALFA, ip3_opt=ip3)
    end if

  END SUBROUTINE ben_writeAmplitude

!--------------------------------------------------------------------------
! ben_setFsoLeadTime
!--------------------------------------------------------------------------
  SUBROUTINE ben_setFsoLeadTime(fsoLeadTime_in)
    implicit none
    real(8)  :: fsoLeadTime_in

    fsoLeadTime = fsoLeadTime_in

  END SUBROUTINE ben_setFsoLeadTime

  !--------------------------------------------------------------------------
  ! ben_getNumStepAmplitudeAssimWindow
  !--------------------------------------------------------------------------
  function ben_getNumStepAmplitudeAssimWindow() result(numStepAmplitude)
    implicit none
    integer numStepAmplitude

    numStepAmplitude = numStepAmplitudeAssimWindow

  end function ben_getNumStepAmplitudeAssimWindow

  !--------------------------------------------------------------------------
  ! ben_getAmplitudeAssimWindow
  !--------------------------------------------------------------------------
  function ben_getAmplitudeAssimWindow() result(adv_amplitude)
    implicit none
    type(struct_adv), pointer  :: adv_amplitude

    adv_amplitude => adv_amplitudeAssimWindow

  end function ben_getAmplitudeAssimWindow

  !--------------------------------------------------------------------------
  ! ben_getAmp3dStepIndexAssimWindow
  !--------------------------------------------------------------------------
  function ben_getAmp3dStepIndexAssimWindow() result(stepIndex)
    implicit none
    integer  :: stepIndex

    stepIndex = amp3dStepIndexAssimWindow

  end function ben_getAmp3dStepIndexAssimWindow

END MODULE BMatrixEnsemble_mod
