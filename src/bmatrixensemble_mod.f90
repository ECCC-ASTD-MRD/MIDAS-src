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
  use mpivar_mod
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
  implicit none
  save
  private

  ! public procedures
  public :: ben_Setup, ben_BSqrt, ben_BSqrtAd, ben_writeAmplitude
  public :: ben_reduceToMPILocal, ben_reduceToMPILocal_r4, ben_expandToMPIGlobal, ben_expandToMPIGlobal_r4
  public :: ben_getScaleFactor, ben_getnEns, ben_getPerturbation, ben_getEnsMean, ben_Finalize

  logical             :: initialized = .false.

  integer,parameter   :: maxNumLevels=200

  real(8),allocatable :: scaleFactor_M(:), scaleFactor_T(:)

  integer             :: nj,ni,lonPerPE,myLonBeg,myLonEnd,latPerPE,myLatBeg,myLatEnd
  integer,allocatable :: allLonBeg(:), allLatBeg(:)
  integer             :: nLevInc_M,nLevInc_T,nkgdimInc,nLevEns_M,nLevEns_T,nkgdimEns
  integer             :: topLevIndex_M,topLevIndex_T
  integer             :: nEnsOverDimension
  integer             :: cvDim_mpilocal,cvDim_mpiglobal
  integer             :: numStep, numStepAmplitude
  integer             :: numSubEns

  integer,external    :: get_max_rss, omp_get_thread_num

  ! Localizations
  integer, parameter  :: maxNumLocalLength = 20
  integer             :: nWaveBand

  ! Ensemble perturbations
  type(struct_ens), allocatable :: ensPerts(:)

  ! Ensemble amplitude (only used in diagnostic mode)
  type(struct_ens)    :: ensAmplitudeStorage

  ! Localization
  integer, allocatable :: locIDs(:)

  logical            :: HUcontainsLQ_gsv

  ! Vertical grid
  type(struct_vco),pointer :: vco_anl, vco_ens, vco_file => null()
  integer                  :: Vcode_anl, Vcode_ens

  ! Horizontal grid
  type(struct_hco), pointer :: hco_anl  ! Analysis   horizontal grid parameters
  type(struct_hco), pointer :: hco_ens  ! Ensemble   horizontal grid parameters
  type(struct_hco), pointer :: hco_file ! Input file horizontal grid parameters

  ! Amplitude advection
  real(8)              :: advectAmplitudeFactor
  integer, allocatable :: lonIndexAdvect(:,:,:,:)
  integer, allocatable :: latIndexAdvect(:,:,:,:)
  real(8), allocatable :: interpWeightAdvect_BL(:,:,:,:)
  real(8), allocatable :: interpWeightAdvect_BR(:,:,:,:)
  real(8), allocatable :: interpWeightAdvect_TL(:,:,:,:)
  real(8), allocatable :: interpWeightAdvect_TR(:,:,:,:)

  ! Namelist variables
  integer             :: nEns ! number of ensemble members
  real(8)             :: scaleFactor(maxNumLevels)
  real(8)             :: scaleFactorHumidity(maxNumLevels)
  integer             :: ntrunc
  character(len=256)  :: enspathname,ensfilebasename
  real(8)             :: hLocalize(maxNumLocalLength)
  real(8)             :: vLocalize(maxNumLocalLength)
  logical             :: TweakTG
  character(len=256)  :: LocalizationType
  integer             :: waveBandPeaks(maxNumLocalLength)
  logical             :: diagnostic
  character(len=2)    :: ctrlVarHumidity
  logical             :: advectAmplitude
  logical             :: removeSubEnsMeans
  logical             :: keepAmplitude

CONTAINS

!--------------------------------------------------------------------------
! ben_setup
!--------------------------------------------------------------------------
  SUBROUTINE ben_setup(hco_anl_in,vco_anl_in,cvDim_out,&
                       mode)
    implicit none
  
    type(struct_hco), pointer, intent(in) :: hco_anl_in
    type(struct_vco), pointer, intent(in) :: vco_anl_in

    character(len=*), intent(in), optional :: mode

    character(len=15) :: ben_mode

    real(8) :: zps

    real(8),pointer :: pressureProfileEns_M(:), pressureProfileFile_M(:), pressureProfileInc_M(:)

    integer        :: cvDim_out, myMemberBeg,myMemberEnd,myMemberCount,maxMyMemberCount
    integer        :: levIndex,nIndex,mIndex,jvar,ila,return_code,status
    integer        :: fnom,fclos,ierr,nulnam
    integer        :: waveBandIndex,locID
    integer        :: stamp_last,newdate,ndate,ntime
    character(len=256) :: ensFileName

    logical        :: EnsTopMatchesAnlTop
    logical        :: lExists

    !namelist
    NAMELIST /NAMBEN/nEns,scaleFactor,scaleFactorHumidity,ntrunc,enspathname,ensfilebasename, &
                     hLocalize,vLocalize,tweakTG,LocalizationType,waveBandPeaks, &
                     diagnostic,ctrlVarHumidity,advectAmplitudeFactor,removeSubEnsMeans, &
                     keepAmplitude

    call tmg_start(12,'BEN_SETUP')

    !
    !- 1.1  Read namelist-dependent options
    !

    ! parameters from namelist
    scaleFactor(:)   = 0.0d0
    scaleFactorHumidity(:) = 1.0d0
    nEns             = 10
    ntrunc           = 30
    enspathname      = 'ensemble'
    ensfilebasename  = ''
    LocalizationType = 'LevelDependent'
    waveBandPeaks(:) =   -1.0d0
    diagnostic       = .false.
    hLocalize(:)     =   -1.0d0
    hLocalize(1)     = 2800.0d0
    vLocalize(:)     =   -1.0d0
    vLocalize(1)     =    2.0d0
    tweakTG          = .false.
    ctrlVarHumidity  = 'LQ'
    advectAmplitudeFactor = 0.0D0
    removeSubEnsMeans= .false.
    keepAmplitude    = .false. 

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
    if ( present(mode) ) then
       if ( trim(mode) == 'Analysis' .or. trim(mode) == 'BackgroundCheck') then
         ben_mode = trim(mode)
         if (mpi_myid == 0) write(*,*)
         if (mpi_myid == 0) write(*,*) 'ben_setup: Mode activated = ', trim(ben_mode)
       else
          write(*,*)
          write(*,*) 'mode = ', trim(mode)
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
      call ens_fileName(ensFileName, ensPathName, ensFileBaseName, 1)
      call vco_SetupFromFile(vco_file, ensFileName)
    end if
    call vco_mpiBcast(vco_file)

    zps = 101000.D0
    nullify(pressureProfileInc_M)
    status = vgd_levels( vco_anl%vgrid, ip1_list=vco_anl%ip1_M, levels=pressureProfileInc_M, &
                         sfc_field=zps, in_log=.false.)
    if (status /= VGD_OK) call utl_abort('ben_setup: ERROR from vgd_levels')
    nullify(pressureProfileFile_M)
    status = vgd_levels( vco_file%vgrid, ip1_list=vco_file%ip1_M, levels=pressureProfileFile_M, &
                         sfc_field=zps, in_log=.false.)
    if (status /= VGD_OK) call utl_abort('ben_setup: ERROR from vgd_levels')

    EnsTopMatchesAnlTop = abs( log(pressureProfileFile_M(1)) - log(pressureProfileInc_M(1)) ) < 0.1d0
    write(*,*) 'EnsTopMatchesAnlTop: EnsTopMatchesAnlTop, presEns, presInc = ', &
               EnsTopMatchesAnlTop, pressureProfileFile_M(1), pressureProfileInc_M(1)
    deallocate(pressureProfileFile_M)
    deallocate(pressureProfileInc_M)

    if( EnsTopMatchesAnlTop ) then
      if( mpi_myid == 0 ) write(*,*) 'ben_setup: top level of ensemble member and analysis grid match'
      vco_ens => vco_anl  ! IMPORTANT: top levels DO match, therefore safe
                          ! to force members to be on analysis vertical levels
    else
      if( mpi_myid == 0 ) write(*,*) 'ben_setup: top level of ensemble member and analysis grid are different, therefore'
      if( mpi_myid == 0 ) write(*,*) '           assume member is already be on correct levels - NO CHECKING IS DONE'
      vco_ens => vco_file ! IMPORTANT: top levels do not match, therefore must
                          ! assume file is already on correct vertical levels
    end if
    status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)
    status = vgd_get(vco_ens%vgrid,key='ig_1 - vertical coord code',value=Vcode_ens)
    if (Vcode_anl.ne.Vcode_ens) then
      write(*,*) 'ben_setup: vcode_anl = ', vcode_anl, ', vcode_ens = ', vcode_ens
      call utl_abort('ben_setup: vertical levels of ensemble not compatible with analysis grid')
    end if
    nLevEns_M = vco_ens%nLev_M
    nLevEns_T = vco_ens%nLev_T
    nLevInc_M = vco_anl%nLev_M
    nLevInc_T = vco_anl%nLev_T
    nkgdimInc = 2*nLevInc_M+2*nLevInc_T+2  ! assume 4 3d and 2 2d variables
    nkgdimEns = 2*nLevEns_M+2*nLevEns_T+2  ! assume 4 3d and 2 2d variables
    topLevIndex_M = nLevInc_M-nLevEns_M+1
    topLevIndex_T = nLevInc_T-nLevEns_T+1

    if (Vcode_anl == 5002) then
      if (nLevEns_T.ne.(nLevEns_M+1)) then
        write(*,*) 'ben_setup: nLevEns_T, nLevEns_M = ',nLevEns_T,nLevEns_M
        call utl_abort('ben_setup: Vcode=5002, nLevEns_T must equal nLevEns_M+1!')
      end if
    else if (Vcode_anl == 5005) then
      if (nLevEns_T.ne.nLevEns_M) then
        write(*,*) 'ben_setup: nLevEns_T, nLevEns_M = ',nLevEns_T,nLevEns_M
        call utl_abort('ben_setup: Vcode=5005, nLevEns_T must equal nLevEns_M!')
      end if
    else
      write(*,*) 'Vcode_anl = ',Vcode_anl
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
    if (Vcode_anl == 5002) then
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
    call mpivar_setup_latbands(nj,latPerPE,myLatBeg,myLatEnd)
    call mpivar_setup_lonbands(ni,lonPerPE,myLonBeg,myLonEnd)
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

      zps = 101000.D0
      nullify(pressureProfileInc_M)
      status = vgd_levels( vco_anl%vgrid, ip1_list=vco_anl%ip1_M, levels=pressureProfileInc_M, &
                           sfc_field=zps, in_log=.false.)
      if (status /= VGD_OK)then
        call utl_abort('ben_setup: ERROR from vgd_levels')
      end if

      allocate(pressureProfileEns_M(nLevEns_M))
      pressureProfileEns_M(1:nLevEns_M) = pressureProfileInc_M(topLevIndex_M:nLevInc_M)

      allocate(locIDs(nWaveBand))
      do waveBandIndex = 1, nWaveBand
         call loc_setup(hco_ens, vco_ens, nEns, pressureProfileEns_M, nTrunc, 'spectral',       & ! IN
                        LocalizationType, hLocalize(waveBandIndex), hLocalize(waveBandIndex+1), & ! IN
                        vLocalize(waveBandIndex),                                               & ! IN
                        cvDim_mpilocal, locID)                                                    ! OUT
         locIDs(waveBandIndex) = locID
      end do

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
    
    ! Read the ensemble data
    call setupEnsemble()

    ! Pre-compute everything for advectAmplitude (only if numStep > 1)
    if(advectAmplitudeFactor == 0.0D0 .or. numStep == 1) then
      if(mpi_myid == 0) write(*,*) 'ben_setup: advection not activated'
      advectAmplitude = .false.
      numStepAmplitude = 1
    else
      if(mpi_myid == 0) write(*,*) 'ben_setup: advection activated'
      if(advectAmplitudeFactor < 0.0D0) advectAmplitudeFactor = 0.0D0 ! FOR TESTING
      advectAmplitude = .true.
      numStepAmplitude = numStep
      call setupAdvectAmplitude
    end if

    ! Special treatment for TG (skin temperature) ensemble perturbations
    if (tweakTG) call AdjustTGOverOpenWater() 

    ! Compute and write Std. Dev.
    if (diagnostic) call EnsembleDiagnostic('FullPerturbations')

    if ( trim(ben_mode) == 'Analysis' ) then

      ! Partitioned the ensemble perturbations into wave bands
      if (trim(LocalizationType) == 'ScaleDependent') then
        call EnsembleScaleDecomposition()
        if (diagnostic) call EnsembleDiagnostic('WaveBandPerturbations')
      end if

      cvDim_out = cvDim_mpilocal
    else
      cvDim_out = 9999 ! Dummy value > 0 to indicate to the background check (s/r compute_HBHT_ensemble) 
                       ! that Bens is used
    end if

    !- Setup en ensGridStateVector to store the amplitude fields (for writing)
    if (keepAmplitude) then
      write(*,*)
      write(*,*) 'ben_setup: ensAmplitude fields will be store for potential write to file'
      call ens_allocate(ensAmplitudeStorage, nEns, numStepAmplitude, hco_ens, vco_ens, &
                        varName='ALFA')
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
    logical :: makeBiPeriodic
    character(len=4) :: varName

    write(*,*) 'setupEnsemble: Start'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !- 1. Memory allocation
    allocate(ensPerts(nWaveBand))
    do waveBandIndex = 1, nWaveBand
      call ens_allocate(ensPerts(waveBandIndex), nEns, numStep, hco_ens, vco_ens)
    end do

    !- 2. Read ensemble
    makeBiPeriodic = (trim(LocalizationType) == 'ScaleDependent')
    call ens_readEnsemble(ensPerts(1), ensPathName, ensFileBaseName, makeBiPeriodic, ctrlVarHumidity)

    !- 3. From ensemble FORECASTS to ensemble PERTURBATIONS

    !- 3.1 remove mean
    call ens_computeMean( ensPerts(1), removeSubEnsMeans, numSubEns=numSubEns )
    call ens_removeMean( ensPerts(1) )

    !- 3.2 normalize and apply scale factors
!$OMP PARALLEL DO PRIVATE (levIndex,varName,lev,ptr4d_r4,stepIndex,memberIndex,multFactor)
    do levIndex = 1, ens_getNumK(ensPerts(1))
      varName = ens_getVarNameFromK(ensPerts(1),levIndex)
      lev = ens_getLevFromK(ensPerts(1),levIndex)

      if ( .not. ens_varExist(ensPerts(1), varName) ) cycle 

      ptr4d_r4 => ens_getRepack_r4(ensPerts(1),levIndex)

      do stepIndex = 1, numStep
        do memberIndex = 1, nEns

          if ( vnl_varLevelFromVarname(varName) == 'MM' ) then
            multFactor = scaleFactor_M(lev)
          elseif ( vnl_varLevelFromVarname(varName) == 'TH' ) then
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
                                 upwardExtrapolationMethod, waveBandIndexWanted, &
                                 undoNormalization)
    implicit none

    type(struct_gsv) :: statevector
    integer,          intent(in) :: memberIndexWanted
    character(len=*), intent(in) :: upwardExtrapolationMethod
    integer, optional, intent(in):: waveBandIndexWanted
    logical, optional :: undoNormalization

    real(8), pointer :: ptr4d_r8(:,:,:,:)
    real(4), pointer :: repack_r4(:,:,:,:)
    real(8) :: dnens2, scaleFactor_MT
    logical :: undoNormalization2
    integer :: waveBandIndex
    integer :: lonIndex,latIndex,stepIndex,levIndex,lev,levInc,topLevOffset
    character(len=4) :: varName

    if ( trim(upwardExtrapolationMethod) /= "ConstantValue" ) then
       call utl_abort('ben_getPerturbation : Invalid value for upwardExtrapolationMethod')
    end if

    if ( present(waveBandIndexWanted) ) then
       waveBandIndex = waveBandIndexWanted
    else
       waveBandIndex = 1
    end if

    ! set default value for optional argument undoNormalization
    if ( present(undoNormalization) ) then
      undoNormalization2 = undoNormalization
    else
      undoNormalization2 = .false.
    end if

    do levIndex = 1, ens_getNumK(ensPerts(1))
      varName = ens_getVarNameFromK(ensPerts(1),levIndex)
      lev = ens_getLevFromK(ensPerts(1),levIndex)

      ptr4d_r8 => gsv_getField_r8(statevector, varName)
      repack_r4 => ens_getRepack_r4(ensPerts(waveBandIndex),levIndex)

!$OMP PARALLEL DO PRIVATE(stepIndex,topLevOffset,scaleFactor_MT,levInc,dnens2,latIndex,lonIndex)
      do stepIndex = 1, numStep

        if ( vnl_varLevelFromVarname(varName) == 'MM' ) then
          topLevOffset = topLevIndex_M - 1
          scaleFactor_MT = scaleFactor_M(lev)
        elseif ( vnl_varLevelFromVarname(varName) == 'TH' ) then
          topLevOffset = topLevIndex_T - 1
          scaleFactor_MT = scaleFactor_T(lev)
        else ! SF
          topLevOffset = 0
          scaleFactor_MT = scaleFactor_T(nLevEns_T)
        end if

        levInc = lev + topLevOffset

        ! undo the normalization (optional)
        if (undoNormalization2) then
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
              dnens2*dble(repack_r4(memberIndexWanted,stepIndex,lonIndex,latIndex))
          end do
        end do

        if ( topLevOffset > 0 .and. lev == 1) then
          ! Fill the gap between the ensemble lid and the analysis lid

          ! undo the normalization (optional)
          if (undoNormalization2) then
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
                  dble(repack_r4(memberIndexWanted,stepIndex,lonIndex,latIndex))
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
    real(8), pointer :: repack_mean(:,:,:)
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
      repack_mean => ens_getRepackMean_r8(ensPerts(1), 1, levIndex)

!$OMP PARALLEL DO PRIVATE(stepIndex,topLevOffset,levInc,latIndex,lonIndex)
      do stepIndex = 1, numStep

        if ( vnl_varLevelFromVarname(varName) == 'MM' ) then
          topLevOffset = topLevIndex_M - 1
        elseif ( vnl_varLevelFromVarname(varName) == 'TH' ) then
          topLevOffset = topLevIndex_T - 1
        else ! SF
          topLevOffset = 0
        end if

        levInc = lev + topLevOffset

        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            ptr4d_out(lonIndex,latIndex,levInc,stepIndex) = repack_mean(stepIndex,lonIndex,latIndex)
          end do
        end do

        if ( topLevOffset > 0 .and. lev == 1 ) then
          ! Fill the gap between the ensemble lid and the analysis lid

          do levInc = 1, topLevOffset
            ! using a constant value
            do latIndex = myLatBeg, myLatEnd
              do lonIndex = myLonBeg, myLonEnd
                ptr4d_out(lonIndex,latIndex,levInc,stepIndex) = repack_mean(stepIndex,lonIndex,latIndex)
              end do
            end do
          end do

        end if ! topLevOffset > 0

      end do ! stepIndex
!$OMP END PARALLEL DO

    end do ! levIndex

  END SUBROUTINE ben_getEnsMean

!--------------------------------------------------------------------------
! AdjustTGOverOpenWater
!--------------------------------------------------------------------------
  SUBROUTINE AdjustTGOverOpenWater()
    implicit none

    real(4), allocatable :: TrialLandSeaMask(:,:), TrialSeaIceMask(:,:)
    real(4), allocatable :: AnalLandSeaMask(:,:), AnalSeaIceMask(:,:)

    real(4), allocatable :: tgstdbg_tmp(:,:)
    real(4), pointer     :: repack_TG(:,:,:,:), repack_TT2m(:,:,:,:)

    integer :: latIndex, lonIndex, stepIndex, memberIndex

    integer :: ini,inj,ink, inpas, inbits, idatyp, ideet
    integer :: ip1,ip2,ip3,ig1,ig2,ig3,ig4
    integer :: ierr,ipas,ntrials
    integer :: idateo
    integer :: fstprm,fstinf,iultg,fnom,fclos,fstouv,fstfrm
    integer :: TrlmNumberWanted
    integer :: fstlir, fstecr, key, fstinl, nultrl, ni_trial, nj_trial
    integer :: deet, npas, nbits, datyp
    integer :: swa, lng, dltf, ubc
    integer :: extra1, extra2, extra3
    
    integer :: ezdefset, ezqkdef, ezsint, ezsetopt, TrialGridID

    character(len=2)  :: typvar
    character(len=1)  :: grtyp
    character(len=4)  :: nomvar
    character(len=12) :: etiket

    character(len=2) :: flnum
    character(len=128) :: trialfile

    logical :: trialExists

    integer :: levTG, levTT2m

    if ( mpi_myid == 0 ) write(*,*)
    if ( mpi_myid == 0 ) write(*,*) 'Adjusting TG ensemble perturbations over Open Water'

    !- Read MG and GL in the middle of the assimilation time window
    if ( tim_nStepObs == 1 ) then
       TrlmNumberWanted = 1
    else
       TrlmNumberWanted = nint( (tim_nStepObs + 1.d0) / 2.d0)
    end if
    
    write(flnum,'(I2.2)') TrlmNumberWanted
    trialfile='./trlm_'//trim(flnum)
    inquire(file=trim(trialfile),exist=trialExists)

    if ( .not. trialExists ) then
       if ( mpi_myid == 0 ) write(*,*)
       if ( mpi_myid == 0 ) write(*,*) 'Trial file not found = ', trialfile
       if ( mpi_myid == 0 ) write(*,*) 'Look for an ensemble of trial files '
       
       trialfile='./trlm_'//trim(flnum)//'_0001'
       inquire(file=trim(trialfile),exist=trialExists)
       if ( .not. trialExists ) then
          if ( mpi_myid == 0 ) write(*,*) 'Ensemble trial file not found = ', trialfile
          call utl_abort('BMatrixEnsemble : DID NOT FIND A TRIAL FIELD FILE')
       end if
    end if

    nultrl = 0
    ierr = fnom(nultrl,trim(trialfile),'RND+OLD+R/O',0)
    ierr = fstouv(nultrl,'RND+OLD')

    !- Determine grid size and EZSCINT ID
    idateo    = -1
    etiket = ' '
    ip1      = -1
    ip2      = -1
    ip3      = -1
    typvar = ' '
    nomvar = 'MG'

    key = fstinf( nultrl,                                             & ! IN
                  ni_trial, nj_trial, ink,                            & ! OUT
                  idateo, etiket, ip1, ip2, ip3, typvar, nomvar ) ! IN

    if (key < 0) then
       write(*,*)
       write(*,*) 'bMatrixEnsemble: Unable to find trial field = ',nomvar
       call utl_abort('BMatrixEnsemble')
    end if

    ierr = fstprm( key,                                                 & ! IN
                   idateo, deet, npas, ni_trial, nj_trial, ink, nbits,  & ! OUT
                   datyp, ip1, ip2, ip3, typvar, nomvar, etiket,  & ! OUT
                   grtyp, ig1, ig2, ig3,                              & ! OUT
                   ig4, swa, lng, dltf, ubc, extra1, extra2, extra3 )     ! OUT

    allocate(TrialLandSeaMask(ni_trial, nj_trial))
    allocate(TrialSeaIceMask(ni_trial, nj_trial))

    idateo   = -1
    etiket = ' '
    ip1      = -1
    ip2      = -1
    ip3      = -1
    typvar = ' '
    nomvar = 'MG'
    ierr = fstlir(TrialLandSeaMask, nultrl, ini, inj, ink,  &
                  idateo ,etiket, ip1, ip2, ip3, typvar, nomvar)
    if ( ierr < 0 ) then
       write(*,*)
       write(*,*) 'bMatrixEnsemble: Unable to read trial field = ',nomvar
       call utl_abort('BMatrixEnsemble : fstlir failed')
    end if

    if (ini /= ni_trial .or. inj /= nj_trial) then
       write(*,*)
       write(*,*) 'bMatrixEnsemble: Invalid dimensions for ...'
       write(*,*) 'nomvar      =', trim(nomvar)
       write(*,*) 'etiket      =', trim(etiket)
       write(*,*) 'ip1         =', ip1
       write(*,*) 'Found ni,nj =', ini, inj 
       write(*,*) 'Should be   =', ni_trial, nj_trial
       call utl_abort('bMatrixEnsemble')
    end if

    nomvar = 'GL'
    ierr = fstlir(TrialSeaIceMask, nultrl, ini, inj, ink,  &
                    idateo ,etiket, ip1, ip2, ip3, typvar, nomvar)
    if ( ierr < 0 ) then
       write(*,*)
       write(*,*) 'bMatrixEnsemble: Unable to read trial field = ',nomvar
       call utl_abort('BMatrixEnsemble : fstlir failed')
    end if

    if (ini /= ni_trial .or. inj /= nj_trial) then
       write(*,*)
       write(*,*) 'bMatrixEnsemble: Invalid dimensions for ...'
       write(*,*) 'nomvar      =', trim(nomvar)
       write(*,*) 'etiket      =', trim(etiket)
       write(*,*) 'ip1         =', ip1
       write(*,*) 'Found ni,nj =', ini, inj 
       write(*,*) 'Should be   =', ni_trial, nj_trial
       call utl_abort('bMatrixEnsemble')
    end if

    TrialGridID  = ezqkdef( ni_trial, nj_trial, grtyp, ig1, ig2, ig3, ig4, nultrl )   ! IN
 
    ierr = fstfrm(nultrl)  
    ierr = fclos(nultrl)

    !- Interpolate to the Analysis Grid
    allocate(AnalLandSeaMask(hco_ens%ni, hco_ens%nj))
    allocate(AnalSeaIceMask(hco_ens%ni, hco_ens%nj))

    ierr = ezdefset( hco_ens%EZscintID, TrialGridID ) ! IN,  IN
    ierr = ezsetopt('INTERP_DEGREE', 'NEAREST') ! Nearest-neighbor interpolation

    ierr = ezsint  (AnalLandSeaMask, TrialLandSeaMask) ! OUT, IN
    ierr = ezsint  (AnalSeaIceMask , TrialSeaIceMask ) ! OUT, IN

    deallocate(TrialLandSeaMask)
    deallocate(TrialSeaIceMask)

    !- Modify the TG ensemble perturbations
    levTG   = 2 + 2*nLevEns_M + 2*nLevEns_T
    levTT2m = 2*nLevEns_M + nLevEns_T

    repack_TG => ens_getRepack_r4(ensPerts(1), levTG)
    repack_TT2m => ens_getRepack_r4(ensPerts(1), levTT2m)
    do latIndex = myLatBeg, myLatEnd
       do lonIndex = myLonBeg, myLonEnd

          if ( AnalLandSeaMask(lonIndex,latIndex) <= 0.1 .and. AnalSeaIceMask(lonIndex,latIndex) <= 0.2 ) then
             ! We have an open water point. Replace TG perturbations by 2m T perturbations.
             ! This is done because the EnKF does not perturbed the SST => TG ens pert = 0 over open water.
             do stepIndex = 1, numStep
                do memberIndex = 1, nens
                   repack_TG(memberIndex,stepIndex,lonIndex,latIndex) = repack_TT2m(memberIndex,stepIndex,lonIndex,latIndex)
                end do
             end do
          end if
       end do
    end do

    !- Write the modified TG

    !allocate(tgstdbg_tmp(hco_ens%ni,hco_ens%nj))
    !tgstdbg_tmp(:,:) = 0.0
    !ptr4d_r4 => ens_getField_r4(ensPerts(1),1)
    !do latIndex = myLatBeg, myLatEnd
    !   do lonIndex = myLonBeg, myLonEnd
    !      tgstdbg_tmp(lonIndex,latIndex) = ptr4d_r4(lonIndex,latIndex,levTG,1)
    !   end do
    !end do
    !call rpn_comm_allreduce(tgstdbg_tmp,tgstdbg_tmp,hco_ens%ni*hco_ens%nj,"mpi_real","mpi_sum","GRID",ierr)
    !
    !if ( mpi_myid == 0 ) then
    !
    !    iultg = 0
    !    ierr = fnom(iultg,'tg_enspert_tweaked.fst','RND',0)
    !    ierr = fstouv(iultg,'RND')
    !
    !    ini = hco_ens%ni
    !    inj = hco_ens%nj
    !    ink = 1
    !    ip1 = 0
    !    ip2 = 0
    !    ip3 = 0
    !    idateo = 0
    !    typvar = 'E'
    !    nomvar = 'TG'
    !    etiket = 'MOD_ENS_PERT'
    !    grtyp = 'G'
    !    ig1 = 0
    !    ig2 = 0
    !    ig3 = 0
    !    ig4 = 0
    !    idatyp = 1
    !
    !    ierr = fstecr(tgstdbg_tmp, tgstdbg_tmp, -32, iultg, idateo,         &
    !                   0, 0, ini, inj, ink, ip1, ip2, ip3, typvar,         &
    !                   nomvar,etiket,grtyp,ig1, ig2, ig3, ig4, idatyp, &
    !                   .true.)
    !
    ! end if
    !
    ! tgstdbg_tmp(:,:) = 0.0
    ! ptr4d_r4 => ens_getField_r4(ensPerts(1),1)
    ! do latIndex = myLatBeg, myLatEnd
    !    do lonIndex = myLonBeg, myLonEnd
    !       tgstdbg_tmp(lonIndex,latIndex) = ptr4d_r4(lonIndex,latIndex,levTT2m,1)
    !    end do
    ! end do
    ! call rpn_comm_allreduce(tgstdbg_tmp,tgstdbg_tmp,hco_ens%ni*hco_ens%nj,"mpi_real","mpi_sum","GRID",ierr)
    !
    ! if ( mpi_myid == 0 ) then
    !    typvar = 'E'
    !    nomvar = 'TT'
    !    etiket = 'RAW_ENS_PERT'
    !    ierr = fstecr(tgstdbg_tmp, tgstdbg_tmp, -32, iultg, idateo,         &
    !                   0, 0, ini, inj, ink, 12000, ip2, ip3, typvar,         &
    !                   nomvar,etiket,grtyp,ig1, ig2, ig3, ig4, idatyp, &
    !                   .true.)
    !
    !    do latIndex = 1, hco_ens%nj
    !      do lonIndex = 1,hco_ens%ni
    !         tgstdbg_tmp(lonIndex,latIndex) = AnalLandSeaMask(lonIndex,latIndex)
    !      end do
    !    end do
    !    typvar = 'P'
    !    nomvar = 'MG'
    !    etiket = 'TRIAL2ANAL'
    !    ierr = fstecr(tgstdbg_tmp, tgstdbg_tmp, -32, iultg, idateo,         &
    !                   0, 0, ini, inj, ink, ip1, ip2, ip3, typvar,         &
    !                   nomvar,etiket,grtyp,ig1, ig2, ig3, ig4, idatyp, &
    !                   .true.)
    !
    !    do latIndex = 1, hco_ens%nj
    !      do lonIndex = 1,hco_ens%ni
    !         tgstdbg_tmp(lonIndex,latIndex) = AnalSeaIceMask(lonIndex,latIndex)
    !      end do
    !    end do       
    !    typvar = 'P'
    !    nomvar = 'GL'
    !    etiket = 'TRIAL2ANAL'
    !    ierr = fstecr(tgstdbg_tmp, tgstdbg_tmp, -32, iultg, idateo,         &
    !                   0, 0, ini, inj, ink, ip1, ip2, ip3, typvar,         &
    !                   nomvar,etiket,grtyp,ig1, ig2, ig3, ig4, idatyp, &
    !                   .true.)
    !
    !    ierr = fstfrm(iultg)
    !    ierr = fclos(iultg)
    ! end if
    !
    ! deallocate(tgstdbg_tmp)

    deallocate(AnalLandSeaMask)
    deallocate(AnalSeaIceMask)

  END SUBROUTINE AdjustTGOverOpenWater

!--------------------------------------------------------------------------
! CheckEnsDim
!--------------------------------------------------------------------------
  SUBROUTINE CheckEnsDim(niEns,njEns,nkEns,nomvar)
    implicit none

    integer,      intent(in) :: niEns,njEns,nkEns
    character(*), intent(in) :: nomvar

    if ( niEns /= lonPerPE  .or. &
         njEns /= latPerPE  .or. &
         nkEns /= 1 ) then

       write(*,*) 'Variable :', trim(nomvar)
       write(*,*) 'i-dim = ', niEns, lonPerPE
       write(*,*) 'j-dim = ', njEns, latPerPE
       write(*,*) 'k-dim = ', nkEns, 1
       call utl_abort('Ensemble dimensions are incompatible with the topology and/or the analysis grid')
    end if

  END SUBROUTINE CheckEnsDim

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

    integer :: gstID2, mIndex, nIndex, mymBeg, mymEnd, mynBeg, mynEnd, mymSkip, mynSkip

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

    call utl_abort('EnsembleScaleDecomposition: this subroutine is not implemented yet!')

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
    !- Setup a spectral transform for filtering (nk = nkgdimEns)
    !
    if (hco_ens%global) then
       ! Global mode
       gstID2 = gst_setup(ni,nj,ntrunc,nkgdimEns)
       if (mpi_myid == 0) write(*,*) 'ben : returned value of gstID2 = ',gstID2
       nla_filter = (ntrunc+1)*(ntrunc+2)/2 !nla_mpilocal ! should be recomputed here based on gstID2 info
       nphase_filter = 2
       allocate(nIndex_vec(nla_filter))
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
        call lst_Setup( lst_ben_filter,                    & ! OUT
                        ni, nj, hco_ens%dlon, ntrunc,      & ! IN
                        'LatLonMN', maxlevels_in=nkgdimEns ) ! IN
        nla_filter = lst_ben_filter%nla
        nphase_filter = lst_ben_filter%nphase
    end if

    !
    !- 1.  Scale decomposition for every wave band except for wave band #1
    !
    allocate(ResponseFunction(nla_filter,2:nWaveBand))
    allocate(ensPertSP(nla_filter,nphase_filter,nkgdimEns))
    allocate(ensPertSPfiltered(nla_filter,nphase_filter,nkgdimEns))
    allocate(ensPertGD(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdimEns))

    ensPertSP(:,:,:) = 0.0d0
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
       end do
    end do
    if (hco_ens%global) deallocate(nIndex_vec)

    do stepIndex = 1, numStep ! Loop on ensemble time bin
       do memberIndex = 1, nEns ! Loop on ensemble member
!NOT IMPLEMENTED YET          ptr4d_r4 => ens_getField_r4(ensPerts(1),memberIndex)

          !- 1.2 GridPoint space -> Spectral Space
!$OMP PARALLEL DO PRIVATE (levIndex,latIndex,lonIndex)
          do levIndex = 1, nkgdimEns
             do latIndex = myLatBeg, myLatEnd
                do lonIndex = myLonBeg, myLonEnd
                   ensPertGD(lonIndex,latIndex,levIndex) = dble(ptr4d_r4(lonIndex,latIndex,levIndex,stepIndex))
                end do
             end do
          end do
!$OMP END PARALLEL DO
          if (hco_ens%global) then
             ! Global Mode
             call gst_setID(gstID2) ! IN
             call gst_reespe(ensPertSP, & ! OUT
                             ensPertGD)   ! IN
          else
             ! LAM mode
             kind = 'GridPointToSpectral'
             call lst_VarTransform( lst_ben_filter%id, & ! IN
                                    ensPertSP,           & ! OUT
                                    ensPertGD,           & ! IN 
                                    kind, nkgdimEns )      ! IN
          end if

          !- 1.3 Filtering and transformation back to grid point space 
          do waveBandIndex = nWaveBand, 2, -1 ! Start with the largest scales
             ! Filtering
!$OMP PARALLEL DO PRIVATE (levIndex,p,ila_filter)
             do levIndex = 1, nkgdimEns
                do p = 1, nphase_filter
                   do ila_filter = 1, nla_filter
                      ensPertSPfiltered(ila_filter,p,levIndex) = &
                           ensPertSP(ila_filter,p,levIndex) * ResponseFunction(ila_filter,waveBandIndex)
                   end do
                end do
             end do
!$OMP END PARALLEL DO

             ! Spectral Space -> GridPoint space
             if (hco_ens%global) then
                ! Global Mode
                call gst_setID(gstID2) ! IN
                call gst_speree(ensPertSPfiltered, & ! IN
                                ensPertGD)           ! OUT
             else
                ! LAM mode
                kind = 'SpectralToGridPoint'
                call lst_VarTransform( lst_ben_filter%id, & ! IN
                                       ensPertSPfiltered, & ! IN
                                       ensPertGD,         & ! OUT
                                       kind, nkgdimEns )    ! IN
             end if
!NOT IMPLEMENTED YET             ptr4d_r4 => ens_getField_r4(ensPerts(waveBandIndex),memberIndex)
!$OMP PARALLEL DO PRIVATE (levIndex,latIndex,lonIndex)
             do levIndex = 1, nkgdimEns
                do latIndex = myLatBeg, myLatEnd
                   do lonIndex = myLonBeg, myLonEnd
                      ptr4d_r4(lonIndex,latIndex,levIndex,stepIndex) = sngl(ensPertGD(lonIndex,latIndex,levIndex))
                   end do
                end do
             end do
!$OMP END PARALLEL DO

          end do ! waveBandIndex
       end do ! ensemble member
    end do ! time bins

    deallocate(ensPertGD)
    deallocate(ResponseFunction)
    deallocate(ensPertSP)
    deallocate(ensPertSPfiltered)

    !
    !- 2.  Isolate the smallest scales in waveBandIndex = 1 by difference in grid point space
    !
    allocate(bandSum(myLonBeg:myLonEnd,myLatBeg:myLatEnd))
    do stepIndex = 1, numStep
!$OMP PARALLEL DO PRIVATE (memberindex,levIndex,latIndex,lonIndex,waveBandIndex,bandsum,ptr4d_r4)
       do memberIndex = 1, nEns
          do levIndex = 1, nkgdimEns
             bandSum(:,:) = 0.d0
             do waveBandIndex = 2, nWaveBand
!NOT IMPLEMENTED YET                ptr4d_r4 => ens_getField_r4(ensPerts(waveBandIndex),memberIndex)
                do latIndex = myLatBeg, myLatEnd
                   do lonIndex = myLonBeg, myLonEnd
                     bandSum(lonIndex,latIndex) = bandSum(lonIndex,latIndex) + dble(ptr4d_r4(lonIndex,latIndex,levIndex,stepIndex))
                   end do
                end do
             end do
!NOT IMPLEMENTED YET             ptr4d_r4 => ens_getField_r4(ensPerts(1),memberIndex)
             do latIndex = myLatBeg, myLatEnd
                do lonIndex = myLonBeg, myLonEnd
                   ptr4d_r4(lonIndex,latIndex,levIndex,stepIndex) = sngl(dble(ptr4d_r4(lonIndex,latIndex,levIndex,stepIndex)) - bandSum(lonIndex,latIndex))
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
  SUBROUTINE ben_reduceToMPILocal(cv_mpilocal,cv_mpiglobal,cvDim_mpilocal_out)
    implicit none
    real(8), intent(out) :: cv_mpilocal(cvDim_mpilocal)
    real(8), intent(in)  :: cv_mpiglobal(:)
    integer, intent(out) :: cvDim_mpilocal_out

    call loc_reduceToMPILocal(locIDs(1),cv_mpilocal,cv_mpiglobal, & ! IN
                              cvDim_mpilocal_out)         ! OUT

   
 END SUBROUTINE ben_reduceToMPILocal

!--------------------------------------------------------------------------
! ben_reduceToMPILocal_r4
!--------------------------------------------------------------------------
  SUBROUTINE ben_reduceToMPILocal_r4(cv_mpilocal,cv_mpiglobal,cvDim_mpilocal_out)
    implicit none
    real(4), intent(out) :: cv_mpilocal(cvDim_mpilocal)
    real(4), intent(in)  :: cv_mpiglobal(:)
    integer, intent(out) :: cvDim_mpilocal_out

    call loc_reduceToMPILocal_r4(locIDs(1),cv_mpilocal,cv_mpiglobal, & ! IN
                                 cvDim_mpilocal_out)         ! OUT

 END SUBROUTINE ben_reduceToMPILocal_r4

!--------------------------------------------------------------------------
! ben_expandToMPIGlobal
!--------------------------------------------------------------------------
  SUBROUTINE ben_expandToMPIGlobal(cv_mpilocal,cv_mpiglobal,cvDim_mpiglobal_out)
    implicit none

    real(8), intent(in)  :: cv_mpilocal(cvDim_mpilocal)
    real(8), intent(out) :: cv_mpiglobal(:)
    integer, intent(out) :: cvDim_mpiglobal_out

    call loc_expandToMPIGlobal(locIDs(1), cv_mpilocal,          & ! IN
                               cv_mpiglobal,cvDim_mpiglobal_out)  ! OUT  

  end SUBROUTINE ben_expandToMPIGlobal

!--------------------------------------------------------------------------
! ben_expandToMPIGlobal_r4
!--------------------------------------------------------------------------
  SUBROUTINE ben_expandToMPIGlobal_r4(cv_mpilocal,cv_mpiglobal,cvDim_mpiglobal_out)
    implicit none

    real(4), intent(in)  :: cv_mpilocal(cvDim_mpilocal)
    real(4), intent(out) :: cv_mpiglobal(:)
    integer, intent(out) :: cvDim_mpiglobal_out

    call loc_expandToMPIGlobal_r4(locIDs(1), cv_mpilocal ,         & ! IN
                                  cv_mpiglobal,cvDim_mpiglobal_out)  ! OUT

  end SUBROUTINE ben_expandToMPIGlobal_r4

!--------------------------------------------------------------------------
! ben_BSqrt
!--------------------------------------------------------------------------
  SUBROUTINE ben_BSqrt(controlVector_in,statevector)
    implicit none

    real(8)          :: controlVector_in(cvDim_mpilocal) 
    type(struct_gsv) :: statevector

    real(8), pointer :: ensAmplitudeAll_M(:,:,:,:,:)
    integer   :: ierr, levIndex, latIndex, memberIndex, waveBandIndex
    logical   :: immediateReturn

    call tmg_start(67,'BEN_BARR')
    if(mpi_doBarrier) call rpn_comm_barrier('GRID',ierr)
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

    !
    !- 2.  Compute the analysis increment from Bens
    !
    allocate(ensAmplitudeAll_M(nEnsOverDimension,numStepAmplitude,myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLevEns_M))
    call gsv_zero(statevector)

    do waveBandIndex = 1, nWaveBand !  Loop on WaveBand (for ScaleDependent Localization)

      ! 2.1 Compute the ensemble amplitudes
      call loc_Lsqrt( locIDs(waveBandIndex),controlVector_in, & ! IN
                      ensAmplitudeAll_M(:,1,:,:,:) )            ! OUT

      ! 2.2 Advect the initial time amplitudes
      if ( advectAmplitude ) call advectAmplitude_tl( ensAmplitudeAll_M ) ! INOUT

      if ( keepAmplitude .and. waveBandIndex == 1 ) call copyAmplitude(ensAmplitudeAll_M) ! IN

      ! 2.3 Compute increment by multiplying amplitudes by member perturbations
      call addEnsMember_repack( ensAmplitudeAll_M, statevector,  & ! INOUT 
                                waveBandIndex )        ! IN

    end do ! Loop on WaveBand

    deallocate(ensAmplitudeAll_M)

    !
    !- 3.  Variable transforms
    !
    if ( ctrlVarHumidity == 'HU') then
       ! convert HU to LQ
       call vtr_transform( statevector, & ! INOUT
                           'HUtoLQ_tlm' ) ! IN
    end if

    if (mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if (mpi_myid == 0) write(*,*) 'ben_bsqrt: done'

  END SUBROUTINE ben_BSqrt

!--------------------------------------------------------------------------
! ben_BSqrtAd
!--------------------------------------------------------------------------
  SUBROUTINE ben_BSqrtAd(statevector,controlVector_out)
    implicit none

    real(8)           :: controlVector_out(cvDim_mpilocal) 
    type(struct_gsv)  :: statevector

    real(8), pointer  :: ensAmplitudeAll_M(:,:,:,:,:)
    integer           :: ierr, levIndex, latIndex, memberIndex, waveBandIndex

    !
    !- 1.  Tests
    !
    call tmg_start(67,'BEN_BARR')
    if(mpi_doBarrier) call rpn_comm_barrier('GRID',ierr)
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

    allocate(ensAmplitudeAll_M(nEnsOverDimension,numStepAmplitude,myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLevEns_M))

    !
    !- 3.  Variable transforms
    !
    if ( ctrlVarHumidity == 'HU') then
       ! convert HU to LQ
       call vtr_transform( statevector, & ! INOUT
                           'HUtoLQ_tlm' ) ! IN
    end if

    !
    !- 2.  Compute the analysis increment from Bens
    !
    allocate(ensAmplitudeAll_M(nEnsOverDimension,numStepAmplitude,myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLevEns_M))

    do waveBandIndex = 1, nWaveBand !  Loop on WaveBand (for ScaleDependent Localization)

      ! 2.3 Compute increment by multiplying amplitudes by member perturbations
      call addEnsMemberAd_repack( statevector, ensAmplitudeAll_M,  & ! INOUT
                                  waveBandIndex )                    ! IN
      ! 2.2 Advect the initial time amplitudes
      if ( advectAmplitude ) call advectAmplitude_ad( ensAmplitudeAll_M ) ! INOUT

      ! 2.1 Compute the ensemble amplitudes
      call loc_LsqrtAd( locIDs(waveBandIndex),ensAmplitudeAll_M(:,1,:,:,:), & ! IN
                        controlVector_out )                        ! OUT

    end do ! Loop on WaveBand

    deallocate(ensAmplitudeAll_M)

    if (mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if (mpi_myid == 0) write(*,*) 'ben_bsqrtAd: done'

  END SUBROUTINE ben_BSqrtAd

!--------------------------------------------------------------------------
! setupAdvectAmplitude
!--------------------------------------------------------------------------
  SUBROUTINE setupAdvectAmplitude
    implicit none
    integer :: latIndex0, lonIndex0, latIndex, lonIndex, levIndex, jsubStep, stepIndex, stepIndex0, ierr, gdxyfll
    integer :: nsize, latIndex_mpiglobal, lonIndex_mpiglobal
    integer :: procID, procIDx, procIDy
    integer, allocatable :: numSubSteps(:)
    real(8) :: uu, vv, delT, subDelT, lonAdvect, latAdvect, delx, dely, sumWeight
    real(8) :: uu_p, vv_p, lonAdvect_p, latAdvect_p, Gcoef, Scoef
    real(4) :: lonAdvect_deg_r4, latAdvect_deg_r4, xpos_r4, ypos_r4
    real(8), allocatable :: uu_mpiglobal_tiles(:,:,:,:), uu_mpiglobal(:,:,:)
    real(8), allocatable :: vv_mpiglobal_tiles(:,:,:,:), vv_mpiglobal(:,:,:)
    logical :: verbose
    real(8) :: numGridPts
    real(8) :: latitudePatch
    real(8), pointer  :: ptr3d_r8(:,:,:)
    character(len=64) :: filename
    character(len=3)  :: filenumber
    type(struct_gsv) :: statevector_ensMean1

    !- Set some important values
    verbose = .false.
    numGridPts = 1.0d0 ! used to compute numSubStep
    latitudePatch = 80.0d0 ! this defines latitude where rotated grid used
    delT = 6.0D0*3600.0D0/real(numStep-1,8) ! time between winds (assume 6h window)

    !- Write out the ensemble mean to a standard file
    call ens_copyEnsMean(ensPerts(1),statevector_ensMean1,1)
    do stepIndex = 1, numStep
      write(filenumber,'(i3.3)') stepIndex
      filename = './ensmean_' // filenumber // '.fst'
      call gsv_writeToFileMPI(statevector_ensMean1,filename,'ENSMEAN',indexStep_in=stepIndex)
    end do
    call gsv_deallocate(statevector_ensMean1)

    if(mpi_myid == 0) write(*,*) 'setupAdvectAmplitude: starting'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    allocate(numSubSteps(nj))

    allocate(lonIndexAdvect(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLevEns_M,numStep))
    allocate(latIndexAdvect(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLevEns_M,numStep))
    allocate(interpWeightAdvect_BL(numStep,myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLevEns_M))
    allocate(interpWeightAdvect_BR(numStep,myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLevEns_M))
    allocate(interpWeightAdvect_TL(numStep,myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLevEns_M))
    allocate(interpWeightAdvect_TR(numStep,myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLevEns_M))
    interpWeightAdvect_BL(:,:,:,:) = 1.0d0
    interpWeightAdvect_BR(:,:,:,:) = 0.0d0
    interpWeightAdvect_TL(:,:,:,:) = 0.0d0
    interpWeightAdvect_TR(:,:,:,:) = 0.0d0

    allocate(uu_mpiglobal_tiles(numStep, lonPerPE, latPerPE, mpi_nprocs))
    allocate(uu_mpiglobal(numStep, ni, nj))
    allocate(vv_mpiglobal_tiles(numStep, lonPerPE, latPerPE, mpi_nprocs))
    allocate(vv_mpiglobal(numStep, ni, nj))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    do levIndex = 1, nLevEns_M ! loop over levels in amplitude field 

      if(mpi_myid == 0) write(*,*) 'setupAdvectAmplitude: levIndex = ', levIndex
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      ! gather the global winds for this level
      nsize = lonPerPE*latPerPE*numStep
 
      ptr3d_r8 => ens_getRepackMean_r8(ensPerts(1),1,levIndex + ens_getOffsetFromVarName(ensPerts(1),'UU'))
      call rpn_comm_allgather(ptr3d_r8(:,:,:), nsize, "mpi_double_precision",  &
                              uu_mpiglobal_tiles(:,:,:,:),     nsize, "mpi_double_precision",  &
                              "GRID", ierr )
      ptr3d_r8 => ens_getRepackMean_r8(ensPerts(1),1,levIndex + ens_getOffsetFromVarName(ensPerts(1),'VV'))
      call rpn_comm_allgather(ptr3d_r8(:,:,:), nsize, "mpi_double_precision",  &
                              vv_mpiglobal_tiles(:,:,:,:),     nsize, "mpi_double_precision",  &
                              "GRID", ierr )

      ! rearrange gathered winds for convenience
      do procIDy = 0, (mpi_npey-1)
        do procIDx = 0, (mpi_npex-1)
          procID = procIDx + procIDy*mpi_npex

          do latIndex = 1, latPerPE
            latIndex_mpiglobal = latIndex + allLatBeg(procIDy+1) - 1
            do lonIndex = 1, lonPerPE
              lonIndex_mpiglobal = lonIndex + allLonBeg(procIDx+1) - 1
              uu_mpiglobal(:, lonIndex_mpiglobal, latIndex_mpiglobal) = uu_mpiglobal_tiles(:, lonIndex, latIndex, procID+1)
              vv_mpiglobal(:, lonIndex_mpiglobal, latIndex_mpiglobal) = vv_mpiglobal_tiles(:, lonIndex, latIndex, procID+1)
            end do
          end do

        end do
      end do

      ! determine the number of time steps required as a function of latitude
      do latIndex0 = 1, nj
        latAdvect = hco_ens%lat(latIndex0)
        if(abs(latAdvect) < latitudePatch*MPC_RADIANS_PER_DEGREE_R8) then
          uu = maxval(abs(uu_mpiglobal(:,:,latIndex0) /(RA*cos(latAdvect)))) ! in rad/s
          vv = maxval(abs(vv_mpiglobal(:,:,latIndex0) / RA)) ! in rad/s
        else
          uu = maxval(abs(uu_mpiglobal(:,:,latIndex0) / RA)) ! in rad/s
          vv = maxval(abs(vv_mpiglobal(:,:,latIndex0) / RA)) ! in rad/s
        end if
        numSubSteps(latIndex0) = max( 1,  &
                               nint( (delT * advectAmplitudeFactor * uu) / (numGridPts*(hco_ens%lon(2)-hco_ens%lon(1))) ),  &
                               nint( (delT * advectAmplitudeFactor * vv) / (numGridPts*(hco_ens%lat(2)-hco_ens%lat(1))) ) )
      end do
      if(mpi_myid == 0) write(*,*) 'min and max of numSubSteps',minval(numSubSteps(:)),maxval(numSubSteps(:))

      ! loop over all initial grid points within tile for determining back trajectories
      do latIndex0 = myLatBeg, myLatEnd
        do lonIndex0 = myLonBeg, myLonEnd

          subDelT = delT/real(numSubSteps(latIndex0),8)  ! in seconds

          ! compute backtrajectories starting at time stepIndex0
          do stepIndex0  = 1, numStep

            ! position at the initial time of back trajectory
            lonAdvect = hco_ens%lon(lonIndex0)  ! in radians
            latAdvect = hco_ens%lat(latIndex0)
            lonIndex = lonIndex0  ! index
            latIndex = latIndex0
            xpos_r4 = real(lonIndex,4)
            ypos_r4 = real(latIndex,4)

            ! initial positions in rotated coordinate system
            lonAdvect_p = 0.0d0
            latAdvect_p = 0.0d0
          
          if(mpi_myid == 0 .and. verbose) &
          write(*,*) 'final lonAdvect,latAdvect=', &
                       lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                       latAdvect*MPC_DEGREES_PER_RADIAN_R8
          if(mpi_myid == 0 .and. verbose) &
          write(*,*) 'numSubSteps=', numSubSteps(latIndex0)

            ! time step of back trajectory, stepping backwards
            do stepIndex = (stepIndex0-1), 1, -1

              if(mpi_myid == 0 .and. verbose) &
              write(*,*) 'stepIndex0,stepIndex,lonIndex,latIndex=',stepIndex0,stepIndex,lonIndex,latIndex

              do jsubStep = 1, numSubSteps(latIndex0)

                ! perform one timestep of back trajectory
                if(abs(hco_ens%lat(latIndex0)) < latitudePatch*MPC_RADIANS_PER_DEGREE_R8) then
                  ! points away from pole, handled normally

                  ! determine wind at current location (now at BL point)
                  uu = uu_mpiglobal(stepIndex,lonIndex,latIndex)/(RA*cos(hco_ens%lat(latIndex))) ! in rad/s
                  vv = vv_mpiglobal(stepIndex,lonIndex,latIndex)/RA
                  ! apply user-specified scale factor to advecting winds
                  uu = advectAmplitudeFactor * uu
                  vv = advectAmplitudeFactor * vv

                  ! compute next position
                  lonAdvect = lonAdvect - subDelT*uu  ! in radians
                  latAdvect = latAdvect - subDelT*vv

                  if(mpi_myid==0 .and. verbose) &
                  write(*,*) 'not near pole, lonAdvect,latAdvect,uu,vv=', &
                               lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                               latAdvect*MPC_DEGREES_PER_RADIAN_R8,uu,vv
                else
                  ! points near pole, handled in a special way

                  ! determine wind at current location (now at BL point)
                  uu = uu_mpiglobal(stepIndex,lonIndex,latIndex) ! in m/s
                  vv = vv_mpiglobal(stepIndex,lonIndex,latIndex)

                  ! transform wind vector into rotated coordinate system
                  Gcoef = ( cos(latAdvect)*cos(hco_ens%lat(latIndex0)) + &
                            sin(latAdvect)*sin(hco_ens%lat(latIndex0))*cos(lonAdvect-hco_ens%lon(lonIndex0)) ) / &
                          cos(latAdvect_p)
                  Scoef = ( sin(hco_ens%lat(latIndex0))*sin(lonAdvect-hco_ens%lon(lonIndex0)) ) / &
                          cos(latAdvect_p)
                  uu_p = Gcoef * uu - Scoef * vv ! in m/s
                  vv_p = Scoef * uu + Gcoef * vv 

                  ! apply user-specified scale factor to advecting winds
                  uu_p = advectAmplitudeFactor * uu_p ! in m/s
                  vv_p = advectAmplitudeFactor * vv_p

                  ! compute next position (in rotated coord system)
                  lonAdvect_p = lonAdvect_p - subDelT*uu_p/(RA*cos(latAdvect_p))  ! in radians
                  latAdvect_p = latAdvect_p - subDelT*vv_p/RA

                  if(mpi_myid == 0 .and. verbose) &
                  write(*,*) '    near pole, uu_p,vv_p,Gcoef,Scoef=', &
                               uu_p, vv_p, Gcoef, Scoef

                  if(mpi_myid == 0 .and. verbose) &
                  write(*,*) '    near pole, lonAdvect_p,latAdvect_p=', &
                               lonAdvect_p*MPC_DEGREES_PER_RADIAN_R8, &
                               latAdvect_p*MPC_DEGREES_PER_RADIAN_R8

                  ! compute lon/lat in original coordinate system
                  lonAdvect = hco_ens%lon(lonIndex0) +                                                  &
                              atan2( cos(latAdvect_p)*sin(lonAdvect_p) ,                            &
                                     ( cos(latAdvect_p)*cos(lonAdvect_p)*cos(hco_ens%lat(latIndex0)) -  &
                                       sin(latAdvect_p)*sin(hco_ens%lat(latIndex0)) ) )
                  latAdvect = asin( cos(latAdvect_p)*cos(lonAdvect_p)*sin(hco_ens%lat(latIndex0)) + &
                                    sin(latAdvect_p)*cos(hco_ens%lat(latIndex0)) )

                  if(mpi_myid == 0 .and. verbose) &
                  write(*,*) '    near pole, lonAdvect,latAdvect=', &
                               lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                               latAdvect*MPC_DEGREES_PER_RADIAN_R8
                end if

                ! convert lon/lat position into index
                lonAdvect_deg_r4 = real(lonAdvect,4)* MPC_DEGREES_PER_RADIAN_R4
                latAdvect_deg_r4 = real(latAdvect,4)* MPC_DEGREES_PER_RADIAN_R4
                ierr = gdxyfll(hco_ens%EZscintID, xpos_r4, ypos_r4, &
                               latAdvect_deg_r4, lonAdvect_deg_r4, 1)

                ! determine the bottom-left grid point
                lonIndex = floor(xpos_r4)
                latIndex = floor(ypos_r4)

                ! check if position is east of the grid
                if(floor(xpos_r4) > ni) then
                  if(mpi_myid == 0 .and. verbose) &
                  write(*,*) 'stepIndex0,lonIndex0,latIndex0,lon,lat,stepIndex,x,y xpos_r4 > ni :', &
                              stepIndex0,lonIndex0,latIndex0,lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                              latAdvect*MPC_DEGREES_PER_RADIAN_R8,stepIndex,xpos_r4,ypos_r4
                  ! add 10*epsilon(real*4) to ensure do not go too far due to limited precision
                  lonAdvect = lonAdvect - 2.0D0*MPC_PI_R8 + 10.0D0*real(epsilon(1.0),8)
                  lonAdvect_deg_r4 = real(lonAdvect,4)* MPC_DEGREES_PER_RADIAN_R4
                  latAdvect_deg_r4 = real(latAdvect,4)* MPC_DEGREES_PER_RADIAN_R4
                  ierr = gdxyfll(hco_ens%EZscintID, xpos_r4, ypos_r4, &
                                 latAdvect_deg_r4, lonAdvect_deg_r4, 1)
                  if(mpi_myid == 0 .and. verbose) &
                  write(*,*) 'new                            xpos_r4 > ni :', &
                              stepIndex0,lonIndex0,latIndex0,lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                              latAdvect*MPC_DEGREES_PER_RADIAN_R8,stepIndex,xpos_r4,ypos_r4
                end if

                ! check if position is west of the grid
                if(floor(xpos_r4) < 1) then
                  if(mpi_myid == 0 .and. verbose) &
                  write(*,*) 'stepIndex0,lonIndex0,latIndex0,lon,lat,stepIndex,x,y xpos_r4 <  1 :', &
                              stepIndex0,lonIndex0,latIndex0,lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                              latAdvect*MPC_DEGREES_PER_RADIAN_R8,stepIndex,xpos_r4,ypos_r4
                  ! subtract 10*epsilon(real*4) to ensure do not go too far due to limited precision
                  lonAdvect = lonAdvect + 2.0D0*MPC_PI_R8 - 10.0D0*real(epsilon(1.0),8)
                  lonAdvect_deg_r4 = real(lonAdvect,4)* MPC_DEGREES_PER_RADIAN_R4
                  latAdvect_deg_r4 = real(latAdvect,4)* MPC_DEGREES_PER_RADIAN_R4
                  ierr = gdxyfll(hco_ens%EZscintID, xpos_r4, ypos_r4, &
                                 latAdvect_deg_r4, lonAdvect_deg_r4, 1)
                  if(mpi_myid == 0 .and. verbose) &
                  write(*,*) 'new                            xpos_r4 <  1 :', &
                              stepIndex0,lonIndex0,latIndex0,lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                              latAdvect*MPC_DEGREES_PER_RADIAN_R8,stepIndex,xpos_r4,ypos_r4
                end if

                ! longitude is still outside grid - should not happen!
                if(floor(xpos_r4) > ni) then 
                  write(*,*) '***still outside lonIndex > ni: stepIndex0,stepIndex,jsubStep,lonIndex0,latIndex0,levIndex,x,y,uu=', &
                                                          stepIndex0,stepIndex,jsubStep,lonIndex0,latIndex0,levIndex,xpos_r4,ypos_r4,uu
                  xpos_r4 = real(ni)
                  lonAdvect = hco_ens%lon(ni)
                end if
                if(floor(xpos_r4) <  1) then 
                  write(*,*) '***still outside lonIndex < 1 : stepIndex0,stepIndex,jsubStep,lonIndex0,latIndex0,levIndex,x,y,uu=', &
                                                      stepIndex0,stepIndex,jsubStep,lonIndex0,latIndex0,levIndex,xpos_r4,ypos_r4,uu
                  xpos_r4 = 1.0
                  lonAdvect = hco_ens%lon(1)
                end if

                ! if position is poleward of last lat circle, ensure valid lat index
                if(latIndex > nj) then
                  if(verbose) &
                  write(*,*) 'stepIndex0,lonIndex0,latIndex0,lon,lat,stepIndex,x,y ypos_r4 > nj :', &
                              stepIndex0,lonIndex0,latIndex0,lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                              latAdvect*MPC_DEGREES_PER_RADIAN_R8,stepIndex,xpos_r4,ypos_r4
                  ypos_r4 = real(nj)
                  latAdvect = hco_ens%lat(nj)
                end if

                ! if position is poleward of first lat circle, ensure valid lat index
                if(latIndex < 1) then
                  if(verbose) &
                  write(*,*) 'stepIndex0,lonIndex0,latIndex0,lon,lat,stepIndex,x,y ypos_r4 <  1 :', &
                              stepIndex0,lonIndex0,latIndex0,lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                              latAdvect*MPC_DEGREES_PER_RADIAN_R8,stepIndex,xpos_r4,ypos_r4
                  ypos_r4 = 1.0
                  latAdvect = hco_ens%lat(1)
                end if

                ! determine bottom left grid point again after possible adjustments
                lonIndex = floor(xpos_r4)
                latIndex = floor(ypos_r4)

              end do ! jsubStep

            end do ! stepIndex

            ! store the final position of the back trajectory and interp weights
            lonIndexAdvect(lonIndex0,latIndex0,levIndex,stepIndex0) = lonIndex
            latIndexAdvect(lonIndex0,latIndex0,levIndex,stepIndex0) = latIndex

            if(mpi_myid == 0 .and. verbose) &
            write(*,*) 'final, initial lonIndex,latIndex', lonIndex0,latIndex0,lonIndex,latIndex

            delx = real(xpos_r4,8) - real(lonIndex,8)
            dely = real(ypos_r4,8) - real(latIndex,8)

            interpWeightAdvect_BL(stepIndex0,lonIndex0,latIndex0,levIndex) = min(max( (1.d0-delx) * (1.d0-dely), 0.0d0), 1.0d0)
            interpWeightAdvect_BR(stepIndex0,lonIndex0,latIndex0,levIndex) = min(max(       delx  * (1.d0-dely), 0.0d0), 1.0d0)
            interpWeightAdvect_TL(stepIndex0,lonIndex0,latIndex0,levIndex) = min(max( (1.d0-delx) *       dely , 0.0d0), 1.0d0)
            interpWeightAdvect_TR(stepIndex0,lonIndex0,latIndex0,levIndex) = min(max(       delx  *       dely , 0.0d0), 1.0d0)

            sumWeight = interpWeightAdvect_BL(stepIndex0,lonIndex0,latIndex0,levIndex) + &
                        interpWeightAdvect_BR(stepIndex0,lonIndex0,latIndex0,levIndex) + &
                        interpWeightAdvect_TL(stepIndex0,lonIndex0,latIndex0,levIndex) + &
                        interpWeightAdvect_TR(stepIndex0,lonIndex0,latIndex0,levIndex)
            if(sumWeight > 1.1d0) then
              write(*,*) 'sumWeight > 1.1 : ', sumWeight
              write(*,*) '          BL, BR, TL, TR=',interpWeightAdvect_BL(stepIndex0,lonIndex0,latIndex0,levIndex), &
                                                     interpWeightAdvect_BR(stepIndex0,lonIndex0,latIndex0,levIndex), &
                                                     interpWeightAdvect_TL(stepIndex0,lonIndex0,latIndex0,levIndex), &
                                                     interpWeightAdvect_TR(stepIndex0,lonIndex0,latIndex0,levIndex)
              write(*,*) '          stepIndex0, levIndex, lonIndex0, latIndex0, lonIndex, latIndex, delx, dely  =', &
                                    stepIndex0, levIndex, lonIndex0, latIndex0, lonIndex, latIndex, delx, dely

              interpWeightAdvect_BL(stepIndex0,lonIndex0,latIndex0,levIndex) = 0.25d0
              interpWeightAdvect_BR(stepIndex0,lonIndex0,latIndex0,levIndex) = 0.25d0
              interpWeightAdvect_TL(stepIndex0,lonIndex0,latIndex0,levIndex) = 0.25d0
              interpWeightAdvect_TR(stepIndex0,lonIndex0,latIndex0,levIndex) = 0.25d0

            end if

          end do ! stepIndex0

        end do ! lonIndex0
      end do ! latIndex0
    end do ! levIndex

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if(mpi_myid == 0) write(*,*) 'setupAdvectAmplitude: done'

  end SUBROUTINE setupAdvectAmplitude

!--------------------------------------------------------------------------
! advectAmplitude_tl
!--------------------------------------------------------------------------
  SUBROUTINE advectAmplitude_tl( ensAmplitudeAll_M )
    implicit none
    real(8)              :: ensAmplitudeAll_M(nEnsOverDimension,numStepAmplitude,myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLevEns_M)
    real(8), allocatable :: ensAmplitude1_mpiglobal_tiles(:,:,:,:)
    real(8), allocatable :: ensAmplitude1_mpiglobal(:,:,:)
    integer :: memberIndex, stepIndex, levIndex, lonIndex, latIndex, lonIndex2, latIndex2, lonIndex2_p1, latIndex2_p1, nsize, ierr
    integer :: procID, procIDx, procIDy, lonIndex_mpiglobal, latIndex_mpiglobal

    allocate(ensAmplitude1_mpiglobal_tiles(nEns,lonPerPE,latPerPE,mpi_nprocs))
    allocate(ensAmplitude1_mpiglobal(nEns,ni,nj))

    do levIndex = 1, nLevEns_M

      ! gather the global initial amplitude field on all tasks
      nsize = nEns*lonPerPE*latPerPE
      call rpn_comm_allgather(ensAmplitudeAll_M(1:nEns,1,:,:,levIndex),        nsize, "mpi_double_precision",  &
                              ensAmplitude1_mpiglobal_tiles(:,:,:,:), nsize, "mpi_double_precision",  &
                              "GRID", ierr )

      ! rearrange gathered amplitudes for convenience
!$OMP PARALLEL DO PRIVATE (procIDy,procIDx,procID,latIndex,lonIndex,latIndex_mpiglobal,lonIndex_mpiglobal,memberIndex)
      do procIDy = 0, (mpi_npey-1)
        do procIDx = 0, (mpi_npex-1)
          procID = procIDx + procIDy*mpi_npex
          do latIndex = 1, latPerPE
            latIndex_mpiglobal = latIndex + allLatBeg(procIDy+1) - 1
            do lonIndex = 1, lonPerPE
              lonIndex_mpiglobal = lonIndex + allLonBeg(procIDx+1) - 1
              do memberIndex = 1, nEns
                ensAmplitude1_mpiglobal(memberIndex,lonIndex_mpiglobal, latIndex_mpiglobal) = &
                  ensAmplitude1_mpiglobal_tiles(memberIndex, lonIndex, latIndex, procID+1)
              end do ! memberIndex
            end do ! lonIndex
          end do ! latIndex
        end do ! procIDx
      end do ! procIDy
!$OMP END PARALLEL DO
      
!$OMP PARALLEL DO PRIVATE (latIndex,lonIndex,lonIndex2,latIndex2,lonIndex2_p1,latIndex2_p1,stepIndex,memberIndex)
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          ! this is the bottom-left grid point
          lonIndex2 = lonIndexAdvect(lonIndex,latIndex,levIndex,stepIndex)
          latIndex2 = latIndexAdvect(lonIndex,latIndex,levIndex,stepIndex)
          lonIndex2_p1 = mod(lonIndex2,ni)+1 ! assume periodic
          latIndex2_p1 = min(latIndex2+1,nj)
          do stepIndex = 2, numStep
            do memberIndex = 1, nEns
              ensAmplitudeAll_M(memberIndex,stepIndex,lonIndex,latIndex,levIndex) =   &
                interpWeightAdvect_BL(stepIndex,lonIndex,latIndex,levIndex)*ensAmplitude1_mpiglobal(memberIndex, lonIndex2   ,latIndex2) +  &
                interpWeightAdvect_BR(stepIndex,lonIndex,latIndex,levIndex)*ensAmplitude1_mpiglobal(memberIndex, lonIndex2_p1,latIndex2) +  &
                interpWeightAdvect_TL(stepIndex,lonIndex,latIndex,levIndex)*ensAmplitude1_mpiglobal(memberIndex, lonIndex2   ,latIndex2_p1) +  &
                interpWeightAdvect_TR(stepIndex,lonIndex,latIndex,levIndex)*ensAmplitude1_mpiglobal(memberIndex, lonIndex2_p1,latIndex2_p1)
            end do ! memberIndex
          end do ! stepIndex
        end do ! lonIndex
      end do ! latIndex
!$OMP END PARALLEL DO

    end do ! levIndex

    deallocate(ensAmplitude1_mpiglobal_tiles)
    deallocate(ensAmplitude1_mpiglobal)

  END SUBROUTINE advectAmplitude_tl

!--------------------------------------------------------------------------
! advectAmplitude_ad
!--------------------------------------------------------------------------
  SUBROUTINE advectAmplitude_ad( ensAmplitudeAll_M )
    implicit none
    real(8)              :: ensAmplitudeAll_M(nEnsOverDimension,numStepAmplitude,myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLevEns_M)
    real(8), allocatable :: ensAmplitude1_mpiglobal(:,:,:)
    real(8), allocatable :: ensAmplitude1_mpiglobal_tiles(:,:,:,:)
    real(8), allocatable :: ensAmplitude1_mpiglobal_tiles2(:,:,:,:)
    integer :: memberIndex, stepIndex, levIndex, lonIndex, latIndex, lonIndex2, latIndex2, lonIndex2_p1, latIndex2_p1
    integer :: lonIndex_mpiglobal, latIndex_mpiglobal, procID, procIDx, procIDy, ierr, nsize

    allocate(ensAmplitude1_mpiglobal(nEns,ni,nj))
    allocate(ensAmplitude1_mpiglobal_tiles(nEns,lonPerPE,latPerPE,mpi_nprocs))
    allocate(ensAmplitude1_mpiglobal_tiles2(nEns,lonPerPE,latPerPE,mpi_nprocs))

    do levIndex = 1, nLevEns_M
      ensAmplitude1_mpiglobal(:,:,:) = 0.0d0

!$OMP PARALLEL DO PRIVATE (latIndex,lonIndex,latIndex2,lonIndex2,latIndex2_p1,lonIndex2_p1)
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          ! this is the bottom-left grid point
          lonIndex2 = lonIndexAdvect(lonIndex,latIndex,levIndex,stepIndex)
          latIndex2 = latIndexAdvect(lonIndex,latIndex,levIndex,stepIndex)
          lonIndex2_p1 = mod(lonIndex2,ni)+1 ! assume periodic
          latIndex2_p1 = min(latIndex2+1,nj)
          do stepIndex = 2, numStep
            do memberIndex = 1, nEns
              ensAmplitude1_mpiglobal(memberIndex, lonIndex2   ,latIndex2) = ensAmplitude1_mpiglobal(memberIndex, lonIndex2   ,latIndex2) +  &
                interpWeightAdvect_BL(stepIndex,lonIndex,latIndex,levIndex)*ensAmplitudeAll_M(memberIndex,stepIndex,lonIndex,latIndex,levIndex)
              ensAmplitude1_mpiglobal(memberIndex, lonIndex2_p1,latIndex2) = ensAmplitude1_mpiglobal(memberIndex, lonIndex2_p1,latIndex2) +  &
                interpWeightAdvect_BR(stepIndex,lonIndex,latIndex,levIndex)*ensAmplitudeAll_M(memberIndex,stepIndex,lonIndex,latIndex,levIndex)
              ensAmplitude1_mpiglobal(memberIndex, lonIndex2   ,latIndex2_p1) = ensAmplitude1_mpiglobal(memberIndex, lonIndex2   ,latIndex2_p1) +  &
                interpWeightAdvect_TL(stepIndex,lonIndex,latIndex,levIndex)*ensAmplitudeAll_M(memberIndex,stepIndex,lonIndex,latIndex,levIndex)
              ensAmplitude1_mpiglobal(memberIndex, lonIndex2_p1,latIndex2_p1) = ensAmplitude1_mpiglobal(memberIndex, lonIndex2_p1,latIndex2_p1) +  &
                interpWeightAdvect_TR(stepIndex,lonIndex,latIndex,levIndex)*ensAmplitudeAll_M(memberIndex,stepIndex,lonIndex,latIndex,levIndex)
            end do ! memberIndex
          end do ! stepIndex
        end do ! lonIndex
      end do ! latIndex
!$OMP END PARALLEL DO

      ! redistribute the global initial time field across mpi tasks by tiles
!$OMP PARALLEL DO PRIVATE(procIDy,procIDx,procID,levIndex,latIndex,latIndex_mpiglobal,lonIndex,lonIndex_mpiglobal)
      do procIDy = 0, (mpi_npey-1)
        do procIDx = 0, (mpi_npex-1)
          procID = procIDx + procIDy*mpi_npex
            do latIndex = 1, latPerPE
              latIndex_mpiglobal = latIndex + allLatBeg(procIDy+1) - 1
              do lonIndex = 1, lonPerPE
                lonIndex_mpiglobal = lonIndex + allLonBeg(procIDx+1) - 1
                do memberIndex = 1, nEns
                  ensAmplitude1_mpiglobal_tiles(memberIndex, lonIndex, latIndex, procID+1) =  &
                    ensAmplitude1_mpiglobal(memberIndex, lonIndex_mpiglobal, latIndex_mpiglobal)
                end do ! memberIndex
              end do ! lonIndex
            end do ! latIndex
          end do ! procIDx
        end do ! procIDy
!$OMP END PARALLEL DO

      nsize = nEns*lonPerPE*latPerPE
      if(mpi_nprocs > 1) then
        call rpn_comm_alltoall(ensAmplitude1_mpiglobal_tiles, nsize,"mpi_double_precision",  &
                               ensAmplitude1_mpiglobal_tiles2,nsize,"mpi_double_precision","GRID",ierr)
      else
        ensAmplitude1_mpiglobal_tiles2(:,:,:,1) = ensAmplitude1_mpiglobal_tiles(:,:,:,1)
      end if

      do procID = 0, (mpi_nprocs-1)
!$OMP PARALLEL DO PRIVATE(levIndex,latIndex,latIndex2,lonIndex,lonIndex2)
        do latIndex = 1, latPerPE
          latIndex2= latIndex + myLatBeg - 1
          do lonIndex = 1, lonPerPE
            lonIndex2 = lonIndex + myLonBeg - 1
            do memberIndex = 1, nEns
              ensAmplitudeAll_M(memberIndex, 1, lonIndex2, latIndex2, levIndex) = ensAmplitudeAll_M(memberIndex, 1, lonIndex2, latIndex2, levIndex) +  &
                 ensAmplitude1_mpiglobal_tiles2(memberIndex, lonIndex, latIndex, procID+1)
            end do ! memberIndex
          end do ! lonIndex
        end do ! latIndex
!$OMP END PARALLEL DO
      end do ! procID

    end do ! levIndex

    deallocate(ensAmplitude1_mpiglobal)
    deallocate(ensAmplitude1_mpiglobal_tiles)
    deallocate(ensAmplitude1_mpiglobal_tiles2)

  END SUBROUTINE advectAmplitude_ad

!--------------------------------------------------------------------------
! addEnsMember
!--------------------------------------------------------------------------
  SUBROUTINE addEnsMember(statevector_ensAmp, statevector_out, &
                          waveBandIndex)
    implicit none

    type(struct_gsv)    :: statevector_ensAmp(:)
    type(struct_gsv)    :: statevector_out
    integer, intent(in) :: waveBandIndex

    real(8), pointer    :: ensAmplitude_M(:,:,:,:)
    real(8), pointer    :: ensAmplitude_MT(:,:,:,:)
    real(8), target     :: ensAmplitude_T(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLevEns_T,numStepAmplitude)
    real(8), pointer    :: increment_out(:,:,:,:)
    real(4), pointer    :: ensMember_r4(:,:,:,:)

    integer     :: jvar, levIndex, levIndex2, stepIndex, stepIndex_amp, latIndex, lonIndex, numVar, topLevOffset, numLev, memberIndex

    call utl_abort('ERROR: this subroutine is not used anymore')

    call tmg_start(62,'ADDMEM')

    do memberIndex = 1, nEns

    ensAmplitude_M => gsv_getField_r8(statevector_ensAmp(memberIndex))

    ! compute thermo level amplitudes from momentum level amplitudes
    if (Vcode_anl == 5002) then

!$OMP PARALLEL DO PRIVATE (levIndex)
      do levIndex = 1, nLevEns_T
        if (levIndex == 1) then
          ! use top momentum level amplitudes for top thermo level
          ensAmplitude_T(:,:,levIndex,:) = ensAmplitude_M(:,:,levIndex,:)
        else if (levIndex == nLevEns_T) then
          ! use surface momentum level amplitudes for surface thermo level
          ensAmplitude_T(:,:,levIndex,:) = ensAmplitude_M(:,:,nLevEns_M,:)
        else
          ! for other levels, interpolate momentum weights to get thermo amplitudes
          ensAmplitude_T(:,:,levIndex,:) = 0.5d0*( ensAmplitude_M(:,:,levIndex-1,:) +   &
                                               ensAmplitude_M(:,:,levIndex,:) )
        end if
      end do
!$OMP END PARALLEL DO

    else if (Vcode_anl == 5005) then

!$OMP PARALLEL DO PRIVATE (levIndex)
      do levIndex = 1, nLevEns_T
        if (levIndex == nLevEns_T) then
          ! use surface momentum level amplitudes for surface thermo level (this is not correct, 
          ! but consistent with how it has been done until now)
          ensAmplitude_T(:,:,levIndex,:) = ensAmplitude_M(:,:,levIndex,:)
        else
          ! for other levels, interpolate momentum weights to get thermo amplitudes
          ensAmplitude_T(:,:,levIndex,:) = 0.5d0*( ensAmplitude_M(:,:,levIndex,:) +   &
                                               ensAmplitude_M(:,:,levIndex+1,:) )
        end if
      end do
!$OMP END PARALLEL DO

    end if

    if (Vcode_anl == 5002 .or. Vcode_anl == 5005) then

      do jvar = 1, vnl_numvarmax 
        if (.not. gsv_varExist(statevector_out, vnl_varNameList(jvar))) cycle
        increment_out => gsv_getField_r8(statevector_out, vnl_varNameList(jvar))
!NOT IMPLEMENTED YET        ensMember_r4 => ens_getField_r4(ensPerts(waveBandIndex),memberIndex, &
!                                        vnl_varNameList(jvar))
        if (vnl_varLevelFromVarname(vnl_varNameList(jvar)) == 'SF') then
          numLev = 1
          topLevOffset = 1
          ensAmplitude_MT(myLonBeg:, myLatBeg:, 1:        , 1:) =>   &
           ensAmplitude_M(myLonBeg:, myLatBeg:, nLevEns_M:, 1:)
        else if (vnl_varLevelFromVarname(vnl_varNameList(jvar)) == 'MM') then
          numLev = nlevEns_M
          topLevOffset = topLevIndex_M
          ensAmplitude_MT => ensAmplitude_M
        else
          numLev = nlevEns_T
          topLevOffset = topLevIndex_T
          ensAmplitude_MT => ensAmplitude_T
        end if

        call tmg_start(77,'ADDMEM_INNER')

        do stepIndex = 1, numStep
          if(advectAmplitude) then
            stepIndex_amp = stepIndex
          else
            stepIndex_amp = 1
          end if
!$OMP PARALLEL DO PRIVATE (levIndex,levIndex2,latIndex,lonIndex)
          do levIndex = 1, numLev                 ! levels for the amplitude and ensemble member
            levIndex2 = levIndex - 1 + topLevOffset   ! levels for the 1 variable increment field
            do latIndex = myLatBeg, myLatEnd
              do lonIndex = myLonBeg, myLonEnd
                increment_out(lonIndex,latIndex,levIndex2,stepIndex) = increment_out(lonIndex,latIndex,levIndex2,stepIndex) +   &
                  ensAmplitude_MT(lonIndex,latIndex,levIndex,stepIndex_amp) * dble(ensMember_r4(lonIndex,latIndex,levIndex,stepIndex))
              end do
            end do
          end do
!$OMP END PARALLEL DO
        end do ! stepIndex

        call tmg_stop(77)

      end do ! jvar

    else ! other Vcode

      write(*,*) 'addEnsMember: this Vcode not supported =',Vcode_anl
      call utl_abort('addEnsMember')

    end if

    end do ! memberIndex

    call tmg_stop(62)

  END SUBROUTINE addEnsMember

!--------------------------------------------------------------------------
! addEnsMember_repack
!--------------------------------------------------------------------------
  SUBROUTINE addEnsMember_repack(ensAmplitudeAll_M, statevector_out, &
                                 waveBandIndex)
    implicit none

    real(8), target    :: ensAmplitudeAll_M(nEnsOverDimension,numStepAmplitude,myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLevEns_M)
    type(struct_gsv)    :: statevector_out
    integer, intent(in) :: waveBandIndex

    real(8), allocatable, target :: ensAmplitudeAll_MT(:,:,:,:)
    real(8), pointer     :: ensAmplitudeAll_MT_ptr(:,:,:,:)
    real(8), pointer     :: increment_out(:,:,:,:)
    real(8), allocatable :: increment_out2(:,:,:)
    real(4), pointer     :: ensMemberAll_r4(:,:,:,:)
    integer     :: lev, lev2, levIndex, stepIndex, stepIndex_amp, latIndex, lonIndex, topLevOffset, numLev, memberIndex
    character(len=4)     :: varName

    if (Vcode_anl /= 5002) call utl_abort('addEnsMember_repack: Only 5002 supported for now!')

    call tmg_start(62,'ADDMEM')

    allocate(ensAmplitudeAll_MT(nEns,numStepAmplitude,myLonBeg:myLonEnd,myLatBeg:myLatEnd))
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

        ensAmplitudeAll_MT_ptr(1:,1:,myLonBeg:,myLatBeg:) => ensAmplitudeAll_M(1:nEns,:,:,:,lev)

      else if (vnl_varLevelFromVarname(varName) == 'TH') then

        if (lev == 1) then
          ! use top momentum level amplitudes for top thermo level
          ensAmplitudeAll_MT_ptr(1:,1:,myLonBeg:,myLatBeg:) => ensAmplitudeAll_M(1:nEns,:,:,:,lev)
        else if (lev == nLevEns_T) then
          ! use surface momentum level amplitudes for surface thermo level
          ensAmplitudeAll_MT_ptr(1:,1:,myLonBeg:,myLatBeg:) => ensAmplitudeAll_M(1:nEns,:,:,:,nLevEns_M)
        else
          ! for other levels, interpolate momentum weights to get thermo amplitudes
!$OMP PARALLEL DO PRIVATE (latIndex)
          do latIndex = myLatBeg, myLatEnd
            ensAmplitudeAll_MT(:,:,:,latIndex) = 0.5d0*( ensAmplitudeAll_M(1:nEns,:,:,latIndex,lev-1) +   &
                                                   ensAmplitudeAll_M(1:nEns,:,:,latIndex,lev) )
          end do
!$OMP END PARALLEL DO
          ensAmplitudeAll_MT_ptr(1:,1:,myLonBeg:,myLatBeg:) => ensAmplitudeAll_MT(:,:,:,:)
        end if

      else if (vnl_varLevelFromVarname(varName) == 'SF') then

        ! surface variable
        ensAmplitudeAll_MT_ptr(1:,1:,myLonBeg:,myLatBeg:) => ensAmplitudeAll_M(1:nEns,:,:,:,nLevEns_M)

      end if
      call tmg_stop(66)

      call tmg_start(77,'ADDMEM_INNER')

      ensMemberAll_r4 => ens_getRepack_r4(ensPerts(waveBandIndex),levIndex)
!$OMP PARALLEL DO PRIVATE (latIndex,lonIndex,stepIndex,stepIndex_amp,memberIndex)
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          do stepIndex = 1, numStep
            if(advectAmplitude) then
              stepIndex_amp = stepIndex
            else
              stepIndex_amp = 1
            end if
            do memberIndex = 1, nEns
              increment_out2(stepIndex,lonIndex,latIndex) = increment_out2(stepIndex,lonIndex,latIndex) +   &
                ensAmplitudeAll_MT_ptr(memberIndex,stepIndex_amp,lonIndex,latIndex) *  &
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
!$OMP PARALLEL DO PRIVATE (stepIndex)
      do stepIndex = 1, numStep
        increment_out(:,:,lev2,stepIndex) = increment_out2(stepIndex,:,:)
      end do
!$OMP END PARALLEL DO

    end do ! levIndex

    deallocate(ensAmplitudeAll_MT)
    deallocate(increment_out2)

    call tmg_stop(62)

  END SUBROUTINE addEnsMember_repack

!--------------------------------------------------------------------------
! addEnsMemberAd
!--------------------------------------------------------------------------
  SUBROUTINE addEnsMemberAd(statevector_in, statevector_ensamp, &
                            memberIndex, waveBandIndex)
    implicit none

    type(struct_gsv)   :: statevector_ensamp
    type(struct_gsv)   :: statevector_in
    integer,intent(in) :: memberIndex, waveBandIndex

    real(8), pointer :: ensAmplitude_M(:,:,:,:)
    real(8), pointer :: ensAmplitude_MT(:,:,:,:)
    real(8), target  :: ensAmplitude_T(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLevEns_T,numStepAmplitude)
    real(8), pointer :: increment_in(:,:,:,:)
    real(4), pointer :: ensMember_r4(:,:,:,:)

    integer     :: jvar, levIndex, levIndex2, stepIndex, stepIndex_amp, latIndex, lonIndex, numVar, topLevOffset, numLev

    call utl_abort('addEnsMemberAd: this subroutine is not used anymore!')

    call tmg_start(63,'ADDMEMAD')

    ensAmplitude_M => gsv_getField_r8(statevector_ensAmp)

    ensAmplitude_M(:,:,:,:) = 0.0d0
    ensAmplitude_T(:,:,:,:) = 0.0d0

    if (Vcode_anl == 5002 .or. Vcode_anl == 5005) then


      do jvar = 1, vnl_numvarmax
        if (.not. gsv_varExist(statevector_in, vnl_varNameList(jvar))) cycle
        increment_in => gsv_getField_r8(statevector_in, vnl_varNameList(jvar))
!NOT IMPLEMENTED YET        ensMember_r4 => ens_getField_r4(ensPerts(waveBandIndex),memberIndex, &
!                                        vnl_varNameList(jvar))

        if (vnl_varLevelFromVarname(vnl_varNameList(jvar)) == 'SF') then
          numLev = 1
          topLevOffset = 1
          ensAmplitude_MT(myLonBeg:, myLatBeg:, 1:        , 1:) =>   &
           ensAmplitude_M(myLonBeg:, myLatBeg:, nLevEns_M:, 1:)
        else if (vnl_varLevelFromVarname(vnl_varNameList(jvar)) == 'MM') then
          numLev = nlevEns_M
          topLevOffset = topLevIndex_M
          ensAmplitude_MT => ensAmplitude_M
        else
          numLev = nlevEns_T
          topLevOffset = topLevIndex_T
          ensAmplitude_MT => ensAmplitude_T
        end if

        call tmg_start(78,'ADDMEMAD_INNER')
        do stepIndex = 1, numStep
          if(advectAmplitude) then
            stepIndex_amp = stepIndex
          else
            stepIndex_amp = 1
          end if
!$OMP PARALLEL DO PRIVATE (levIndex,levIndex2,latIndex,lonIndex)
          do levIndex = 1, numLev                 ! levels for the amplitude and ensemble member
            levIndex2 = levIndex - 1 + topLevOffset   ! levels for the 1 variable increment field
            do latIndex = myLatBeg, myLatEnd
              do lonIndex = myLonBeg, myLonEnd
                ensAmplitude_MT(lonIndex,latIndex,levIndex,stepIndex_amp) = ensAmplitude_MT(lonIndex,latIndex,levIndex,stepIndex_amp) +   &
                  increment_in(lonIndex,latIndex,levIndex2,stepIndex) * dble(ensMember_r4(lonIndex,latIndex,levIndex,stepIndex))
              end do
            end do
          end do
!$OMP END PARALLEL DO
        end do ! stepIndex
        call tmg_stop(78)

      end do ! jvar

    else ! other Vcode

      write(*,*) 'addEnsMember: this Vcode not supported =',Vcode_anl
      call utl_abort('addEnsMember')

    end if

    ! combine thermo and momentum level amplitude sensitivites
    if (Vcode_anl == 5002) then

!$OMP PARALLEL DO PRIVATE (levIndex)
      do levIndex = 1, nLevEns_T
        if (levIndex == 1) then
          ! use top momentum level amplitudes for top thermo level
          ensAmplitude_M(:,:,levIndex,:) = ensAmplitude_M(:,:,levIndex,:) + ensAmplitude_T(:,:,levIndex,:)
        else if (levIndex == nLevEns_T) then
          ! use surface momentum level amplitudes for surface thermo level
          ensAmplitude_M(:,:,nLevEns_M,:) = ensAmplitude_M(:,:,nLevEns_M,:) + ensAmplitude_T(:,:,levIndex,:)
        else
          ! for other levels, interpolate momentum weights to get thermo amplitudes
          ensAmplitude_M(:,:,levIndex,:)   = ensAmplitude_M(:,:,levIndex,:)   + 0.5d0*ensAmplitude_T(:,:,levIndex,:)
        end if
      end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE (levIndex)
      do levIndex = 2, (nLevEns_T-1)
        ensAmplitude_M(:,:,levIndex-1,:) = ensAmplitude_M(:,:,levIndex-1,:) + 0.5d0*ensAmplitude_T(:,:,levIndex,:)
      end do
!$OMP END PARALLEL DO

    else if (Vcode_anl == 5005) then

!$OMP PARALLEL DO PRIVATE (levIndex)
      do levIndex = 1, nLevEns_T
        if (levIndex == nLevEns_T) then
          ! use surface momentum level amplitudes for surface thermo level (not correct, but
          ! consistent with how it has been done until now)
          ensAmplitude_M(:,:,levIndex,:) = ensAmplitude_M(:,:,levIndex,:) + ensAmplitude_T(:,:,levIndex,:)
        else
          ! for other levels, interpolate momentum weights to get thermo amplitudes
          ensAmplitude_M(:,:,levIndex,:)   = ensAmplitude_M(:,:,levIndex,:)   + 0.5d0*ensAmplitude_T(:,:,levIndex,:)
        end if
      end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE (levIndex)
      do levIndex = 1, (nLevEns_T-1)
        ensAmplitude_M(:,:,levIndex+1,:) = ensAmplitude_M(:,:,levIndex+1,:) + 0.5d0*ensAmplitude_T(:,:,levIndex,:)
      end do
!$OMP END PARALLEL DO

    end if

    call tmg_stop(63)

  END SUBROUTINE addEnsMemberAd

!--------------------------------------------------------------------------
! addEnsMemberAd_repack
!--------------------------------------------------------------------------
  SUBROUTINE addEnsMemberAd_repack(statevector_in, ensAmplitudeAll_M, &
                                   waveBandIndex)
    implicit none

    real(8)            :: ensAmplitudeAll_M(nEnsOverDimension,numStepAmplitude,myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLevEns_M)
    type(struct_gsv)   :: statevector_in
    integer,intent(in) :: waveBandIndex

    real(8), allocatable :: ensAmplitudeAll_MT(:,:)
    real(8), pointer :: increment_in(:,:,:,:)
    real(8), allocatable :: increment_in2(:,:,:)
    real(4), pointer :: ensMemberAll_r4(:,:,:,:)
    integer          :: levIndex, lev, lev2, stepIndex, stepIndex_amp, latIndex, lonIndex, topLevOffset, numLev, memberIndex
    character(len=4)     :: varName
    !logical :: ensAmpZeroed(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLevEns_M)

    if (Vcode_anl /= 5002) call utl_abort('addEnsMemberAd_repack: Only 5002 supported for now!')

    call tmg_start(63,'ADDMEMAD')

    allocate(ensAmplitudeAll_MT(nEns,numStepAmplitude))
    allocate(increment_in2(numStep,myLonBeg:myLonEnd,myLatBeg:myLatEnd))

    ! set output ensemble Amplitude to zero
    call tmg_start(69,'ADDMEMAD_ZERO')
!$OMP PARALLEL DO PRIVATE (levIndex)
    do levIndex = 1, nLevEns_M
      ensAmplitudeAll_M(:,:,:,:,levIndex) = 0.0d0
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
!$OMP PARALLEL DO PRIVATE (stepIndex)
      do stepIndex = 1, numStep
        increment_in2(stepIndex,:,:) = increment_in(:,:,lev2,stepIndex)
      end do
!$OMP END PARALLEL DO
      call tmg_stop(65)

      !ensAmpZeroed(:,:,:) = .false.
      ensMemberAll_r4 => ens_getRepack_r4(ensPerts(waveBandIndex),levIndex)
!$OMP PARALLEL DO PRIVATE (latIndex,lonIndex,stepIndex,stepIndex_amp,memberIndex,ensAmplitudeAll_MT)
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd

          if(omp_get_thread_num() == 0) call tmg_start(78,'ADDMEMAD_INNER')
          ensAmplitudeAll_MT(:,:) = 0.0d0
          do stepIndex = 1, numStep
            if(advectAmplitude) then
              stepIndex_amp = stepIndex
            else
              stepIndex_amp = 1
            end if
            do memberIndex = 1, nEns
              ensAmplitudeAll_MT(memberIndex,stepIndex_amp) = ensAmplitudeAll_MT(memberIndex,stepIndex_amp) +  &
                increment_in2(stepIndex,lonIndex,latIndex) * dble(ensMemberAll_r4(memberIndex,stepIndex,lonIndex,latIndex))
            end do ! memberIndex
          end do ! stepIndex
          if(omp_get_thread_num() == 0) call tmg_stop(78)

          ! transform thermo/momentum level amplitude sensitivites appropriately

          if(omp_get_thread_num() == 0) call tmg_start(68,'ADDMEMAD_PREPAMP')
          if (vnl_varLevelFromVarname(varName) == 'MM') then

            ensAmplitudeAll_M(1:nEns,:,lonIndex,latIndex,lev) = ensAmplitudeAll_M(1:nEns,:,lonIndex,latIndex,lev) + ensAmplitudeAll_MT(:,:)

          else if (vnl_varLevelFromVarname(varName) == 'TH') then

            if (lev == 1) then
              ! use top momentum level amplitudes for top thermo level
              ensAmplitudeAll_M(1:nEns,:,lonIndex,latIndex,lev) = ensAmplitudeAll_M(1:nEns,:,lonIndex,latIndex,lev) + ensAmplitudeAll_MT(:,:)
            else if (lev == nLevEns_T) then
              ! use surface momentum level amplitudes for surface thermo level
              ensAmplitudeAll_M(1:nEns,:,lonIndex,latIndex,nLevEns_M) = ensAmplitudeAll_M(1:nEns,:,lonIndex,latIndex,nLevEns_M) + ensAmplitudeAll_MT(:,:)
            else
              ! for other levels, interpolate momentum weights to get thermo amplitudes
              ensAmplitudeAll_M(1:nEns,:,lonIndex,latIndex,lev)   = ensAmplitudeAll_M(1:nEns,:,lonIndex,latIndex,lev)   + 0.5d0*ensAmplitudeAll_MT(:,:)
              ensAmplitudeAll_M(1:nEns,:,lonIndex,latIndex,lev-1) = ensAmplitudeAll_M(1:nEns,:,lonIndex,latIndex,lev-1) + 0.5d0*ensAmplitudeAll_MT(:,:)
            end if

          else if (vnl_varLevelFromVarname(varName) == 'SF') then

            ! surface variable
            ensAmplitudeAll_M(1:nEns,:,lonIndex,latIndex,nLevEns_M) = ensAmplitudeAll_M(1:nEns,:,lonIndex,latIndex,nLevEns_M) + ensAmplitudeAll_MT(:,:)

          end if
          if(omp_get_thread_num() == 0) call tmg_stop(68)

        end do ! lonIndex
      end do ! latIndex
!$OMP END PARALLEL DO


    end do ! levIndex

    deallocate(ensAmplitudeAll_MT)
    deallocate(increment_in2)

    call tmg_stop(63)

  END SUBROUTINE addEnsMemberAd_repack

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
                            datestamp=tim_getDatestamp(), mpi_local=.true.)
          call ben_getPerturbation( statevector,    & ! OUT
                                    memberIndex,    & ! IN
                                    'ConstantValue', waveBandIndex ) ! IN
          if ( trim(mode) == 'FullPerturbations') then
             etiket = 'PERT001_FULL'
          else
             write(wbnum,'(I2.2)') waveBandIndex
             etiket = 'PERT001_WB' // trim(wbnum)
          end if
          call gsv_writeToFileMPI(statevector,'./ens_pert001.fst',etiket, & ! IN
                                  dnens2,HUcontainsLQ=HUcontainsLQ_gsv )    ! IN
          call gsv_deallocate(statevector)
       end do
    end if

    !
    !- Compute the standard deviations for each wave band
    !
    if ( mpi_myid == 0 ) write(*,*) '   computing Std.Dev.'
    call gsv_allocate(statevector_temp, tim_nstepobsinc, hco_ens, vco_anl, mpi_local=.true.)

    do waveBandIndex = 1, nWaveBandToDiagnose
       if ( mpi_myid == 0 ) write(*,*) '     waveBandIndex = ', waveBandIndex
       call gsv_allocate(statevector, tim_nstepobsinc, hco_ens, vco_anl, &
                         datestamp=tim_getDatestamp(), mpi_local=.true.)
       call gsv_zero(statevector)
       do memberIndex = 1, nEns
          !- Get normalized perturbations
          call ben_getPerturbation( statevector_temp, & ! OUT
                                    memberIndex,      & ! IN
                                    'ConstantValue', waveBandIndex ) ! IN
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
       call gsv_writeToFileMPI(statevector,'./ens_stddev.fst',etiket, & ! IN
                               HUcontainsLQ=HUcontainsLQ_gsv)           ! IN
       call gsv_deallocate(statevector)
    end do

    call gsv_deallocate(statevector_temp)

  END SUBROUTINE EnsembleDiagnostic

!--------------------------------------------------------------------------
! copyAmplitude
!--------------------------------------------------------------------------
  SUBROUTINE copyAmplitude(ensAmplitude)
    implicit none

    real(8), intent(in) :: ensAmplitude(nEnsOverDimension,numStepAmplitude,myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLevEns_M)

    integer :: ensIndex, stepIndex, lonIndex, latIndex, levIndex

    real(4), pointer :: ensAmpAll_r4(:,:,:,:)

    do levIndex = 1, nLevEns_M
      ensAmpAll_r4 => ens_getRepack_r4(ensAmplitudeStorage,levIndex)
!$OMP PARALLEL DO PRIVATE(ensIndex, latIndex,lonIndex,stepIndex)
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          do stepIndex = 1, numStepAmplitude
            do ensIndex = 1, nEns
              ensAmpAll_r4(ensIndex,stepIndex,lonIndex,latIndex) = &
                   real(ensAmplitude(ensIndex,stepIndex,lonIndex,latIndex,levIndex),4)  
            end do
          end do
        end do
      end do
!$OMP END PARALLEL DO
    end do

  END SUBROUTINE copyAmplitude

!--------------------------------------------------------------------------
! ben_writeAmplitude
!--------------------------------------------------------------------------
  SUBROUTINE ben_writeAmplitude(ensPathName, ensFileBaseName, ip3)
    implicit none

    character(len=*), intent(in) :: ensPathName
    character(len=*), intent(in) :: ensFileBaseName
    integer,          intent(in) :: ip3

    if (initialized .and. keepAmplitude) then
      if ( mpi_myid == 0 ) write(*,*)
      if ( mpi_myid == 0 ) write(*,*) 'bmatrixEnsemble_mod: Writing the amplitude field'
      call ens_writeEnsemble(ensAmplitudeStorage, ensPathName, ensFileBaseName, &
                             'LQ', 'FROM_BENS', 'R',varName='ALFA', ip3_in=ip3)
    end if

  END SUBROUTINE ben_writeAmplitude

END MODULE BMatrixEnsemble_mod
