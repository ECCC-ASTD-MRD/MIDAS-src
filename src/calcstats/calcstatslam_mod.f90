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
!! MODULE CalcStatsLam (prefix="csl")
!!
!! *Purpose*: Compute homogeneous and isotropic background error covariances 
!!            from forecast error estimate in model variable space (limited-area version).
!!
!--------------------------------------------------------------------------
module calcstatslam_mod
  use gridStateVector_mod
  use LamSpectralTransform_mod
  use analysisGrid_mod
  use HorizontalCoord_mod
  use localizationFunction_mod
  use utilities_mod
  use menetrierDiag_mod
  implicit none
  save
  private

  ! Public Subroutines
  public :: csl_setup, csl_computeStats
  public :: csl_toolbox

  type(struct_hco), pointer :: hco_ens => null() ! Ensemble horizontal grid parameters
  type(struct_hco), pointer :: hco_bhi => null() ! B matrix horizontal grid parameters
  type(struct_vco), pointer :: vco_bhi => null() ! B matrix vertical grid parameters
  type(struct_lst)          :: lst_bhi ! Spectral transform Parameters

  character(len=256), allocatable :: cflensin(:)

  integer,external    :: get_max_rss

  integer, parameter :: cv_model = 1
  integer, parameter :: cv_bhi   = 2
  type  :: calcb_cv
    character(len=4)      :: NomVar(2)
    character(len=2)      :: GridType ! TH=Thermo, MM=Momentum, NS=Non-staggered
    integer               :: nlev
    integer               :: kDimStart
    integer               :: kDimEnd
    integer, allocatable  :: ip1(:)
  end type calcb_cv

  character(len=4)      :: momentumControlVar(2)

  integer,parameter      :: nMaxControlVar = 10
  type(calcb_cv)         :: ControlVariable(nMaxControlVar)

  integer                :: nControlVariable
  integer :: nens
  integer :: ntrunc
  integer :: nkgdim
  integer :: ip2_ens

  logical :: initialized = .false.
  logical :: NormByStdDev, SetTGtoZero, writeEnsPert
  logical :: vertLoc, horizLoc, correlationLQ

  character(len=12) :: WindTransform
  character(len=12) :: SpectralWeights
  character(len=2)  :: spatialDimensions

  real(8), pointer     :: pressureProfile_M(:), pressureProfile_T(:)

  real(8) :: vLocalize_wind, vlocalize_mass, vlocalize_humidity ! vertical length scale (in units of ln(Pressure))
  real(8) :: hlocalize_wind, hlocalize_mass, hlocalize_humidity ! vertical length scale (in km)

  contains

!--------------------------------------------------------------------------
! CSL_SETUP
!--------------------------------------------------------------------------
    subroutine csl_setup( nens_in, cflens_in, hco_ens_in, vco_ens_in, ip2_in)
    use vGrid_Descriptors , only: vgrid_descriptor, vgd_levels, VGD_OK  
    implicit none

    integer,                   intent(in)   :: nens_in
    character(len=*),          intent(in)   :: cflens_in(nens_in)
    type(struct_vco), pointer, intent(in)   :: vco_ens_in
    type(struct_hco), pointer, intent(in)   :: hco_ens_in
    integer,                   intent(in)   :: ip2_in

    integer :: nulnam, ier, status
    integer :: fclos, fnom, fstouv, fstfrm
    integer :: grd_ext_x, grd_ext_y
    integer :: var, k

    real(8) :: SurfacePressure

    character(len=4) :: cgneed(nMaxControlVar)

    logical :: mask(nMaxControlVar)

    NAMELIST /NAMCALCSTATS_LAM/ntrunc,grd_ext_x,grd_ext_y,WindTransform,NormByStdDev,&
         SetTGtoZero,writeEnsPert,SpectralWeights,&
         vLocalize_wind,vlocalize_mass,vlocalize_humidity,&
         hLocalize_wind,hlocalize_mass,hlocalize_humidity,correlationLQ
    NAMELIST /NAMSTATE/CGNEED

    write(*,*)
    write(*,*) 'csl_setup: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- 1. Initialized the info on the ensemble
    !
    nens=nens_in

    allocate(cflensin(nens))
    cflensin(:)=cflens_in(:)

    hco_ens => hco_ens_in
    vco_bhi => vco_ens_in

    if ( vco_bhi%nlev_T == 1 .and. vco_bhi%nlev_T == 1 ) then
       write(*,*)
       write(*,*) 'Spatial dimensions = 2D'
       spatialDimensions = '2D'
    else
       write(*,*)
       write(*,*) 'Spatial dimensions = 3D'
       spatialDimensions = '3D'
    end if

    ip2_ens = ip2_in

    !
    !- 2. Read namelist NAMCALCSTATS_LAM
    !
    ntrunc        = 75       ! Default value
    grd_ext_x     = 10       ! Default value
    grd_ext_y     = 10       ! Default value
    WindTransform = 'PsiChi' ! Default value
    NormByStdDev  = .true.   ! Default value
    SetTGtoZero   = .false.  ! Default value
    writeEnsPert  = .false.  ! Default value
    correlationLQ = .false.  ! Default value
    SpectralWeights = 'lst'  ! Default value
    vLocalize_wind      = -1.0d0 ! Default value (no vloc)
    vLocalize_mass      = -1.0d0 ! Default value (no vloc)
    vLocalize_humidity  = -1.0d0 ! Default value (no vloc)
    hLocalize_wind      = -1.0d0 ! Default value (no hloc)
    hLocalize_mass      = -1.0d0 ! Default value (no hloc)
    hLocalize_humidity  = -1.0d0 ! Default value (no hloc)

    nulnam = 0

    ier=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read (nulnam,nml=namcalcstats_lam)
    write(*     ,nml=namcalcstats_lam)
    ier=fclos(nulnam)

    write(*,*)
    write(*,*) 'Truncation = ', ntrunc

    select case(trim(SpectralWeights))
    case ('lst')
       write(*,*)
       write(*,*) 'Using spectral weights from the lamSpectralTransform module'
    case ('legacy')
       write(*,*)
       write(*,*) 'Using spectral weights from the OLD VAR code in LAM mode'
    case default
       write(*,*)
       write(*,*) 'Unknown spectral weights TYPE : ', trim(SpectralWeights)
       call utl_abort('csl_setup')
    end select

    !
    !- 3. Initialized the extended (bi-periodic) grid
    !

    !- 3.1 Create the extended and non-extended grid prototype file
    call CreateLamTemplateGrids('./analysisgrid',grd_ext_x,grd_ext_y) ! IN

    !- 3.2 Setup the Extended B_HI grid
    call hco_SetupFromFile( hco_bhi,'./analysisgrid', 'ANALYSIS', 'BHI' ) ! IN

    !- 3.3 Setup the LAM analysis grid metrics
    call agd_SetupFromHCO( hco_bhi, hco_ens ) ! IN
    !- 3.4 Setup the LAM spectral transform
    call lst_Setup( lst_bhi,                & ! OUT
                    hco_bhi%ni, hco_bhi%nj, & ! IN
                    hco_bhi%dlon, ntrunc,   & ! IN
                    'NoMpi' )                 ! IN

    !
    !- 4.  Setup the control variables (model space and B_hi space)
    !

    !- 4.1 Read NAMELIST NAMSTATE to find which fields are needed
    cgneed(:) = 'NONE'
    ier=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read (nulnam,nml=namstate)
    write(*     ,nml=namstate)
    ier=fclos(nulnam)

    mask = cgneed .ne. 'NONE'
    nControlVariable = count(mask) 
    write(*,*)
    write(*,*) 'Number of Control Variables = ', nControlVariable

    !- 4.2 Set ControlVariable structure
    nkgdim = 0
    momentumControlVar(:) = 'NULL' 
    do var = 1, nControlVariable

       !- Set variable name
       ControlVariable(var)%nomvar(cv_model)  = trim(cgneed(var))
       if ( ControlVariable(var)%nomvar(cv_model) == 'UU' ) then
          if ( trim(WindTransform) == 'PsiChi') then
             write(*,*)
             write(*,*) '--- Momentum Control Variables = Psi-Chi ---'
             momentumControlVar(1) = 'PP'
          else if ( trim(WindTransform) == 'VortDiv') then
             write(*,*)
             write(*,*) '--- Momentum Control Variables = Vort-Div ---'
             momentumControlVar(1) = 'QR'
          else if ( trim(WindTransform) == 'UV') then
             write(*,*)
             write(*,*) '--- Momentum Control Variables = U-V ---'
             momentumControlVar(1) = 'UU'
          else
             write(*,*)
             write(*,*) 'Wind Transform not available ', trim((WindTransform))
             call utl_abort('csl_setup')
          end if
          ControlVariable(var)%nomvar(cv_bhi) = trim(momentumControlVar(1))
       else if ( ControlVariable(var)%nomvar(cv_model) == 'VV' ) then
          if ( trim(WindTransform) == 'PsiChi') then
             momentumControlVar(2) = 'CC'
          else if ( trim(WindTransform) == 'VortDiv') then
             momentumControlVar(2) = 'DD'
          else if ( trim(WindTransform) == 'UV') then
             momentumControlVar(2) = 'VV'
          else
             write(*,*)
             write(*,*) 'Wind Transform not available ', trim((WindTransform))
             call utl_abort('csl_setup')
          end if
          ControlVariable(var)%nomvar(cv_bhi) = momentumControlVar(2)
       else
          ControlVariable(var)%nomvar(cv_bhi)    = ControlVariable(var)%nomvar(cv_model)
       end if

       !- Set Level info
       if ( ControlVariable(var)%nomvar(cv_model) /= 'TG' .and. spatialDimensions /= '2D' ) then
          call VarLevInfo(ControlVariable(var)%nlev,             & ! OUT
                          ControlVariable(var)%GridType,         & ! OUT
                          ControlVariable(var)%nomvar(cv_model), & ! IN 
                          cflensin(1) )                            ! IN
       else
          ControlVariable(var)%nlev     = 1
          ControlVariable(var)%GridType = 'NS'
       end if

       write(*,*)
       write(*,*) 'Control Variable Name ', ControlVariable(var)%nomvar(cv_model)
       write(*,*) '   Number of Levels = ', ControlVariable(var)%nlev
       write(*,*) '   Type   of Levels = ', ControlVariable(var)%GridType
       
       allocate( ControlVariable(var)%ip1 (ControlVariable(var)%nlev) )

       if ( spatialDimensions /= '2D' ) then
          if (ControlVariable(var)%nlev /= 1) then
             if (ControlVariable(var)%GridType == 'TH') then
                ControlVariable(var)%ip1(:) = vco_bhi%ip1_T(:)
             else
                ControlVariable(var)%ip1(:) = vco_bhi%ip1_M(:)
             end if
          else
             ControlVariable(var)%ip1(:) = 0
          end if
       else
          ControlVariable(var)%ip1(:) = vco_bhi%ip1_T(:) ! which is identical to ip1_M(:)
       end if

       ControlVariable(var)%kDimStart = nkgdim + 1
       nkgdim = nkgdim + ControlVariable(var)%nlev
       ControlVariable(var)%kDimEnd   = nkgdim

    end do

    !
    !- 5.  Setup localization 
    !
    
    !- 5.1 Setup horizontal localization
    if ( hLocalize_wind < 0.d0 .and. hLocalize_mass < 0.d0 .and. hLocalize_humidity < 0.d0 ) then
       write(*,*) 
       write(*,*) 'csl_setup: NO horizontal correlation localization will be performed'
       horizLoc=.false.
    else
       write(*,*) 
       write(*,*) 'csl_setup: horizontal correlation localization WILL BE performed'
       horizLoc=.true.
    end if

    !- 5.2 Setup vertical localization
    if ( vLocalize_wind < 0.d0 .and. vLocalize_mass < 0.d0 .and. vLocalize_humidity < 0.d0 ) then
       write(*,*) 
       write(*,*) 'csl_setup: NO vertical correlation localization will be performed'
       vertLoc=.false.
    else
       write(*,*) 
       write(*,*) 'csl_setup: vertical correlation localization WILL BE performed'
       vertLoc=.true.
    end if

    !
    !- 6.  Setup pressure profile for vertical localization
    !
    SurfacePressure = 101000.D0

    status = vgd_levels( vco_bhi%vgrid, ip1_list=vco_bhi%ip1_M,     & ! IN
                         levels=pressureProfile_M,                  & ! OUT
                         sfc_field=SurfacePressure, in_log=.false.)   ! IN

    if ( status /= VGD_OK ) then
      write(*,*)
      write(*,*) 'csl_setup: ERROR with vgd_levels for MOMENTUM levels '
      call utl_abort('csl_setup')
    else
      write(*,*)
      write(*,*) 'Pressure profile...'
      do k = 1, vco_bhi%nlev_M
        write(*,*) k, pressureProfile_M(k) / 100.d0, ' hPa'
      end do
    endif
    
    status = vgd_levels( vco_bhi%vgrid, ip1_list=vco_bhi%ip1_T,     & ! IN
                         levels=pressureProfile_T,                  & ! OUT
                         sfc_field=SurfacePressure, in_log=.false.)   ! IN
    
    if ( status /= VGD_OK ) then
      write(*,*)
      write(*,*) 'csl_setup: ERROR with vgd_levels for THERMO levels '
      call utl_abort('csl_setup')
    else
      write(*,*)
      write(*,*) 'Pressure profile...'
      do k = 1, vco_bhi%nlev_T
        write(*,*) k, pressureProfile_T(k) / 100.d0, ' hPa'
      end do
    endif

    !
    !- 7.  Ending
    !
    initialized = .true.

    write(*,*)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'csl_setup: Done!'

  end subroutine csl_setup

!--------------------------------------------------------------------------
! CSL_ComputeStats
!--------------------------------------------------------------------------
  subroutine csl_computeStats
    implicit none

    integer :: ier

    integer :: NumBins2d, NumBins2dGridPoint

    real(4),allocatable :: ensPerturbations(:,:,:,:)

    real(8),allocatable :: ensMean3d(:,:,:)
    real(8),allocatable :: StdDev3d(:,:,:)
    real(8),allocatable :: StdDev3dGridPoint(:,:,:)
    real(8),allocatable :: SpVertCorrel(:,:,:)
    real(8),allocatable :: TotVertCorrel(:,:)
    real(8),allocatable :: NormB(:,:,:)
    real(8),allocatable :: NormBsqrt(:,:,:)
    real(8),allocatable :: PowerSpectrum(:,:)
    real(8),allocatable :: HorizScale(:)
    
    integer,allocatable :: Bin2d(:,:)
    integer,allocatable :: Bin2dGridPoint(:,:)

    write(*,*)
    write(*,*) 'csl_computeStats: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    allocate(ensPerturbations(hco_bhi%ni,hco_bhi%nj,nkgdim,nens),stat=ier)
    if(ier /= 0) then
      write(*,*) 'Problem allocating ensPerturbations memory!',ier
      write(*,*)  hco_bhi%ni,hco_bhi%nj,nkgdim,nens
      call utl_abort('csl_computeStats') 
    endif
    allocate(ensMean3d(hco_bhi%ni,hco_bhi%nj,nkgdim),stat=ier)
    if (ier /= 0) then
      write(*,*) 'Problem allocating ensMean3d memory!',ier
      call utl_abort('csl_computeStats')
    end if
    allocate(Bin2d(hco_bhi%ni,hco_bhi%nj),stat=ier)
    if (ier /= 0) then
      write(*,*) 'Problem allocating Bin2d memory!',ier
      call utl_abort('csl_computeStats')
    end if
    allocate(Bin2dGridPoint(hco_bhi%ni,hco_bhi%nj),stat=ier)
    if (ier /= 0) then
      write(*,*) 'Problem allocating Bin2d memory!',ier
      call utl_abort('csl_computeStats')
    end if
    allocate(StdDev3d(hco_bhi%ni,hco_bhi%nj,nkgdim),stat=ier)
    if (ier /= 0) then
      write(*,*) 'Problem allocating StdDev3d memory!',ier
      call utl_abort('csl_computeStats')
    end if
    allocate(StdDev3dGridPoint(hco_bhi%ni,hco_bhi%nj,nkgdim),stat=ier)
    if (ier /= 0) then
      write(*,*) 'Problem allocating StdDev3dGridPoint memory!',ier
      call utl_abort('csl_computeStats')
    end if
    allocate(SpVertCorrel(nkgdim,nkgdim,0:ntrunc),stat=ier)
    if(ier.ne.0) then
      write(*,*) 'Problem allocating SpVertCorrel memory!',ier
      call utl_abort('csl_computeStats')
    endif
    allocate(TotVertCorrel(nkgdim,nkgdim),stat=ier)
    if(ier.ne.0) then
      write(*,*) 'Problem allocating TotVertCorrel memory!',ier
      call utl_abort('csl_computeStats')
    endif
    allocate(PowerSpectrum(nkgdim,0:ntrunc),stat=ier)
    if(ier.ne.0) then
      write(*,*) 'Problem allocating PowerSpectrum memory!',ier
      call utl_abort('csl_computeStats')
    endif
    allocate(NormB(nkgdim,nkgdim,0:ntrunc),stat=ier)
    if(ier.ne.0) then
      write(*,*) 'Problem allocating Bsqrt memory!',ier
      call utl_abort('csl_computeStats')
    endif
    allocate(NormBsqrt(nkgdim,nkgdim,0:ntrunc),stat=ier)
    if(ier.ne.0) then
      write(*,*) 'Problem allocating Bsqrt memory!',ier
      call utl_abort('csl_computeStats')
    endif
    allocate(HorizScale(nkgdim),stat=ier)
    if(ier.ne.0) then
      write(*,*) 'Problem allocating HorizScale memory!',ier
      call utl_abort('csl_computeStats')
    endif

    ensPerturbations(:,:,:,:) = 0.d0

    !
    !- 1.  Read forecast error estimates (all variables, levels, ens. members)
    !
    call ReadEnsemble(ensPerturbations) ! INOUT

    !
    !- 2.  Create periodic fields in X and Y directions
    !
    call BiPeriodization(ensPerturbations) ! INOUT

    if (writeEnsPert) call WriteEnsemble(ensPerturbations, './EnsPert_cv_model', cv_model ) ! IN

    !
    !- 3.  Transform u-wind and v-wind to control variables 
    !
    if ( momentumControlVar(1) /= 'NULL' .and. momentumControlVar(2) /= 'NULL' ) then
       call UVToCtrlVar(ensPerturbations) ! INOUT
    end if

    if (writeEnsPert) call WriteEnsemble(ensPerturbations, './EnsPert_cv_bhi', cv_bhi ) ! IN

    !
    !- 4.  Remove the ensemble mean
    !
    call RemoveMean(ensPerturbations, & ! INOUT
                    ensMean3d)          ! OUT

    !
    !- 5.  Calculate the Standard Deviations in grid point space
    !

    !- 5.1 Binned Std. Dev. to be used in the analysis step
    call CreateBins(Bin2d, NumBins2d, & ! OUT
                    'HorizontalMean')   ! IN

    call CalcStdDev3d(ensPerturbations, & ! INOUT
                      StdDev3d,         & ! OUT
                      Bin2d, NumBins2d)   ! IN

    !- 5.2 Non-binned Std. Dev. to be used for normalization and/or
    !      diagnostics
    call CreateBins(Bin2dGridPoint, NumBins2dGridPoint, & ! OUT
                    'GridPoint')                      ! IN

    call CalcStdDev3d(ensPerturbations,             & ! INOUT
                      StdDev3dGridPoint,              & ! OUT
                      Bin2dGridPoint, NumBins2dGridPoint) ! IN

    !- 5.3 Normalization
    if ( NormByStdDev ) then
      call Normalize3d(ensPerturbations, & ! INOUT
                        StdDev3dGridPoint)    ! IN
    end if

    !
    !- 6.  Covariance statistics in Spectral Space
    !

    !- 6.1 Vertical correlations and Power Spectra
    call CalcSpectralStats(ensPerturbations,              & ! IN
                           SpVertCorrel, PowerSpectrum,   & ! OUT
                           NormB)                           ! OUT

    !- 6.2 Calculate the horiontal correlation lenght scales
    call CalcHorizScale(HorizScale, & ! OUT
                        NormB)        ! IN

    !- 6.3 Calculate the total vertical correlation matrix
    call CalcTotVertCorrel(TotVertCorrel,             & ! OUT
                           SpVertCorrel, PowerSpectrum) ! IN

    !- 6.4 Set cross-correlations
    call SetSpVertCorrel(NormB) ! INOUT

    !- 6.5 Calculate the square-root of the correlation-based B matrix
    call CalcBsqrt(NormBsqrt, & ! OUT
                   NormB   )    ! IN

    !
    !- 7.  Writing statistics to files
    !

    !- 7.1 Statistics needed by VAR in analysis mode
    call WriteVarStats(NormBsqrt, StdDev3d) ! IN

    !- 7.2 Diagnostics fields
    call WriteDiagStats(NormB, SpVertCorrel, TotVertCorrel, EnsMean3d, & ! 
                        StdDev3dGridPoint, PowerSpectrum, HorizScale)    ! IN

    deallocate(Bin2d)
    deallocate(Bin2dGridPoint)
    deallocate(ensPerturbations)
    deallocate(StdDev3d)
    deallocate(StdDev3dGridPoint)
    deallocate(ensMean3d)
    deallocate(SpVertCorrel)
    deallocate(TotVertCorrel)
    deallocate(PowerSpectrum)
    deallocate(NormB)
    deallocate(NormBsqrt)
    deallocate(HorizScale)

    write(*,*)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'csl_computeStats: Done!'

  end subroutine csl_computeStats

!--------------------------------------------------------------------------
! CSL_TOOLBOOX
!--------------------------------------------------------------------------
  subroutine csl_toolbox
    implicit none

    integer :: nulnam, ier, fstouv, fnom, fstfrm, fclos
    integer :: iunstats

    integer :: NumBins2d

    real(4),allocatable :: ensPerturbations(:,:,:,:)

    real(8),allocatable :: ensMean3d(:,:,:)
    real(8),allocatable :: StdDev3dGridPoint(:,:,:)
    real(8),allocatable :: VertCorrel(:,:,:)

    real(8),allocatable :: SpVertCorrel(:,:,:)
    real(8),allocatable :: NormB(:,:,:)
    real(8),allocatable :: PowerSpectrum(:,:)

    integer,allocatable :: Bin2d(:,:)

    character(len=4), allocatable :: nomvar3d(:,:),nomvar2d(:,:)

    integer :: varLevOffset(6), nvar3d, nvar2d, var3d, var2d, var

    character(len=60) :: tool

    NAMELIST /NAMTOOLBOX/tool

    write(*,*)
    write(*,*) 'csl_toolbox: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Tool selection
    !
    nulnam = 0
    ier = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=NAMTOOLBOX)
    write(*,nml=NAMTOOLBOX)
    ier = fclos(nulnam)

    select case(trim(tool))
    case ('VERTCORREL_GRIDPOINT')
       write(*,*)
       write(*,*) 'Computing vertical correlation in grid point space'
    case ('HORIZCORREL_FUNCTION')
       write(*,*)
       write(*,*) 'Computing horizontal correlation functions'
    case ('STDDEV')
       write(*,*)
       write(*,*) 'Computing Standard-Deviations'
    case ('LOCALIZATIONRADII')
       write(6,*)
       write(6,*) 'Estimating the optimal covariance localization radii'

       nVar3d = 0
       nVar2d = 0
       do var = 1, nControlVariable
         if ( ControlVariable(var)%nlev == 1 ) then
           nVar2d = nVar2d + 1
         else
           nVar3d = nVar3d + 1
         end if
       end do
       allocate(nomvar3d(nVar3d,3))
       allocate(nomvar2d(nVar2d,3))

       var3d = 0
       var2d = 0
       do var = 1, nControlVariable
         if ( ControlVariable(var)%nlev == 1 ) then
           var2d = var2d + 1
           nomvar2d(var2d,1) = ControlVariable(var)%nomvar(cv_model)
           nomvar2d(var2d,2) = ControlVariable(var)%nomvar(cv_bhi)
         else
           var3d = var3d + 1
           nomvar3d(var3d,1) = ControlVariable(var)%nomvar(cv_model)
           nomvar3d(var3d,2) = ControlVariable(var)%nomvar(cv_bhi)
         end if
         varLevOffset(var) = ControlVariable(var)%kDimStart - 1
       end do
       call bmd_setup( hco_ens, nens, vco_bhi%nlev_M, vco_bhi%nlev_T,    & ! IN
                       nkgdim, pressureProfile_M, pressureProfile_T,     & ! IN
                       nvar3d, nvar2d, varLevOffset, nomvar3d, nomvar2d, & ! IN
                       1)
    case default
       write(*,*)
       write(*,*) 'Unknown TOOL in csl_toolbox : ', trim(tool)
       call utl_abort('calbmatrix_lam')
    end select

    !
    !- Horizontal and vertical correlation diagnostics
    !
    allocate(ensPerturbations(hco_bhi%ni,hco_bhi%nj,nkgdim,nens),stat=ier)
    if(ier /= 0) then
      write(*,*) 'Problem allocating ensPerturbations memory!',ier
      write(*,*)  hco_bhi%ni,hco_bhi%nj,nkgdim,nens
      call utl_abort('csl_computeStats') 
    endif
    allocate(ensMean3d(hco_bhi%ni,hco_bhi%nj,nkgdim),stat=ier)
    if (ier /= 0) then
      write(*,*) 'Problem allocating ensMean3d memory!',ier
      call utl_abort('csl_computeStats')
    end if
    allocate(StdDev3dGridPoint(hco_bhi%ni,hco_bhi%nj,nkgdim),stat=ier)
    if (ier /= 0) then
      write(*,*) 'Problem allocating StdDev3dGridPoint memory!',ier
      call utl_abort('csl_computeStats')
    end if
    allocate(Bin2d(hco_bhi%ni,hco_bhi%nj),stat=ier)
    if (ier /= 0) then
      write(*,*) 'Problem allocating Bin2d memory!',ier
      call utl_abort('csl_computeStats')
    end if

    ensPerturbations(:,:,:,:) = 0.d0
    
    !
    !- 1.  Read forecast error estimates (all variables, levels, ens. members)
    !
    call ReadEnsemble(ensPerturbations) ! INOUT

    !
    !- 2.  Create periodic fields in X and Y directions
    !
    call BiPeriodization(ensPerturbations) ! INOUT

    !
    !- 3.  Transform u-wind and v-wind to control variables 
    !
    if ( momentumControlVar(1) /= 'NULL' .and. momentumControlVar(2) /= 'NULL' ) then
       call UVToCtrlVar(ensPerturbations) ! INOUT
    end if

    !
    !- 4.  Remove the ensemble mean
    !
    call RemoveMean(ensPerturbations, & ! INOUT
                    ensMean3d)          ! OUT

    !
    !- 5.  Calculate the Standard Deviations in grid point space
    !
    
    !- 5.1 Non-binned Std. Dev. to be used for normalization
    call CreateBins(Bin2d, NumBins2d, & ! OUT
                    'GridPoint')        ! IN

    call CalcStdDev3d(ensPerturbations,   & ! INOUT
                      StdDev3dGridPoint,  & ! OUT
                      Bin2d, NumBins2d)     ! IN

    if ( trim(tool) == 'STDDEV' ) then
       iunstats = 0
       ier    = fnom(iunstats,'./stddev.fst','RND',0)
       ier    = fstouv(iunstats,'RND')
       call Write3d(StdDev3dGridPoint,iunstats,'STDDEV_GRIDP',cv_bhi) ! IN
       call WriteTicTacToc(iunstats) ! IN
       ier =  fstfrm(iunstats)
       ier =  fclos (iunstats)
       return
    end if

    if ( trim(tool) == 'VERTCORREL_GRIDPOINT' ) then
       !
       !- 6.  VERTCORREL_GRIDPOINT
       !
       
       !- 6.0 Normalization
       call Normalize3d(ensPerturbations, & ! INOUT
                        StdDev3dGridPoint)  ! IN

       !- 6.1  Remove the (2D) domain mean
       call removeDomainMean(ensPerturbations) ! INOUT

       !- 6.2  Calculate the vertical correlations
       iunstats = 0
       ier    = fnom(iunstats,'./vertCorrel.fst','RND',0)
       ier    = fstouv(iunstats,'RND')
       
       !- 6.2.1 Stats from the whole horizontal domain
       call CreateBins(Bin2d, NumBins2d, & ! OUT
                      'FullDomain')        ! IN

       allocate(VertCorrel(nkgdim,nkgdim,NumBins2d))

       call calcVertCorrel(VertCorrel,       &  ! OUT
                           ensPerturbations, &  ! IN
                           Bin2d, NumBins2d)    ! IN

       call WriteTotVertCorrel(VertCorrel(:,:,1),iunstats,'ZT','GPVCOR_FULL') ! IN

       deallocate(VertCorrel)

       !- 6.2.2 Stats from the core and the extension zone
       call CreateBins(Bin2d, NumBins2d, & ! OUT
                       'CoreExt')          ! IN

       allocate(VertCorrel(nkgdim,nkgdim,NumBins2d))

       call calcVertCorrel(VertCorrel,       &  ! OUT
                           ensPerturbations, &  ! IN
                           Bin2d, NumBins2d)    ! IN

       call WriteTotVertCorrel(VertCorrel(:,:,1),iunstats,'ZT','GPVCOR_CORE') ! IN
       call WriteTotVertCorrel(VertCorrel(:,:,2),iunstats,'ZT','GPVCOR_EXTN') ! IN

       deallocate(VertCorrel)

       ier =  fstfrm(iunstats)
       ier =  fclos (iunstats)

     else if (trim(tool) == 'HORIZCORREL_FUNCTION') then
       
       !
       !- 7.  HORIZCORREL_FUNCTION
       !
       allocate(SpVertCorrel(nkgdim,nkgdim,0:ntrunc))
       allocate(PowerSpectrum(nkgdim,0:ntrunc))
       allocate(NormB(nkgdim,nkgdim,0:ntrunc))

       !- 7.0 Normalization
       call Normalize3d(ensPerturbations, & ! INOUT
                        StdDev3dGridPoint)  ! IN

       !- 7.1 Vertical correlations and Power Spectra
       call CalcSpectralStats(ensPerturbations,              & ! IN
                              SpVertCorrel, PowerSpectrum,   & ! OUT
                              NormB)                           ! OUT

       !- 7.2 Compute the horizontal correlation functions
       call horizCorrelFunction(NormB) ! IN

       deallocate(SpVertCorrel)
       deallocate(PowerSpectrum)
       deallocate(NormB)

    else
      call bmd_localizationRadii(ensPerturbations, StdDev3dGridPoint, cv_bhi, waveBandIndex=1) ! IN
    endif

    !
    !- Write the estimated pressure profiles
    !
    call writePressureProfiles

    deallocate(Bin2d)
    deallocate(ensPerturbations)
    deallocate(StdDev3dGridPoint)
    deallocate(ensMean3d)

    write(*,*)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'csl_toolbox: Done!'

  end subroutine csl_toolbox

!--------------------------------------------------------------------------
! CreateLamTemplateGrids
!--------------------------------------------------------------------------
  subroutine CreateLamTemplateGrids(TemplateFileName,grd_ext_x,grd_ext_y)
    use MathPhysConstants_mod, only : MPC_DEGREES_PER_RADIAN_R8
    implicit none

    character(len=*), intent(in) :: TemplateFileName
    integer         , intent(in) :: grd_ext_x
    integer         , intent(in) :: grd_ext_y

    integer :: ni_ext, nj_ext, i, j
    integer :: iun = 0
    integer :: ier, fnom, fstouv, fstfrm, fclos

    real(8), allocatable :: Field2d(:,:)

    real(8), allocatable :: lat_ext(:)
    real(8), allocatable :: lon_ext(:)

    real(8) :: dlat, dlon

    integer :: dateo,npak
    integer :: ip1,ip2,ip3,deet,npas,datyp,ig1,ig2,ig3,ig4
    integer :: ig1_tictac,ig2_tictac,ig3_tictac,ig4_tictac

    character(len=1)  :: grtyp
    character(len=2)  :: typvar
    character(len=12) :: etiket

    !
    !- 1.  Opening the output template file
    !
    ier = fnom(iun, trim(TemplateFileName), 'RND', 0)
    ier = fstouv(iun, 'RND')

    npak     = -32

    !
    !- 2.  Writing the core grid (Ensemble) template
    !

    !- 2.1 Tic-Tac
    deet     =  0
    ip1      =  hco_ens%ig1
    ip2      =  hco_ens%ig2
    ip3      =  0
    npas     =  0
    datyp    =  1
    grtyp    = 'E'
    typvar   = 'X'
    etiket   = 'COREGRID'
    dateo =  0

    call cxgaig ( grtyp,                                          & ! IN
                  ig1_tictac, ig2_tictac, ig3_tictac, ig4_tictac, & ! OUT
                  real(hco_ens%xlat1), real(hco_ens%xlon1),   & ! IN
                  real(hco_ens%xlat2), real(hco_ens%xlon2)  )   ! IN

    ig1      =  ig1_tictac
    ig2      =  ig2_tictac
    ig3      =  ig3_tictac
    ig4      =  ig4_tictac

    ier = utl_fstecr(hco_ens%lon*MPC_DEGREES_PER_RADIAN_R8, npak, &
                     iun, dateo, deet, npas, hco_ens%ni, 1, 1, ip1,    &
                     ip2, ip3, typvar, '>>', etiket, grtyp, ig1,          &
                     ig2, ig3, ig4, datyp, .true.)

    ier = utl_fstecr(hco_ens%lat*MPC_DEGREES_PER_RADIAN_R8, npak, &
                     iun, dateo, deet, npas, 1, hco_ens%nj, 1, ip1,    &
                     ip2, ip3, typvar, '^^', etiket, grtyp, ig1,          &
                     ig2, ig3, ig4, datyp, .true.)

    !- 2.2 2D Field
    allocate(Field2d(hco_ens%ni,hco_ens%nj))
    Field2d(:,:) = 10.d0

    deet      =  0
    ip1       =  0
    ip2       =  0
    ip3       =  0
    npas      =  0
    datyp     =  1
    grtyp     =  hco_ens%grtyp
    typvar    = 'A'
    etiket    = 'COREGRID'
    dateo  =  0
    ig1       =  hco_ens%ig1
    ig2       =  hco_ens%ig2
    ig3       =  hco_ens%ig3
    ig4       =  hco_ens%ig4

    ier = utl_fstecr(Field2d, npak,                                    &
                     iun, dateo, deet, npas, hco_ens%ni, hco_ens%nj, 1, ip1, &
                     ip2, ip3, typvar, 'P0', etiket, grtyp, ig1,                &
                     ig2, ig3, ig4, datyp, .true.)

    deallocate(Field2d)

    !
    !- 3.  Create and Write the extended grid (Analysis) template
    !
    ni_ext = hco_ens%ni + grd_ext_x
    nj_ext = hco_ens%nj + grd_ext_y

    !- 3.1 Tic-Tac
    allocate(lon_ext(ni_ext))
    allocate(lat_ext(nj_ext))

    !- Copy core grid info
    lon_ext(1:hco_ens%ni) = hco_ens%lon(:) 
    lat_ext(1:hco_ens%nj) = hco_ens%lat(:)

    !- Extend the lat lon
    dlon = hco_ens%lon(2) - hco_ens%lon(1) 
    do i = hco_ens%ni + 1, ni_ext
       lon_ext(i) = lon_ext(hco_ens%ni) + (i - hco_ens%ni) * dlon
    end do

    dlat = hco_ens%lat(2) - hco_ens%lat(1) 
    do j = hco_ens%nj + 1, nj_ext
       lat_ext(j) = lat_ext(hco_ens%nj) + (j - hco_ens%nj) * dlat
    end do

    !- Write
    deet     =  0
    ip1      =  hco_ens%ig1 + 100 ! Must be different from the core grid
    ip2      =  hco_ens%ig2 + 100 ! Must be different from the core grid
    ip3      =  0
    npas     =  0
    datyp    =  1
    grtyp    = 'E'
    typvar   = 'X'
    etiket   = 'ANALYSIS'
    dateo =  0
    ig1      =  ig1_tictac
    ig2      =  ig2_tictac
    ig3      =  ig3_tictac
    ig4      =  ig4_tictac

    ier = utl_fstecr(lon_ext*MPC_DEGREES_PER_RADIAN_R8, npak, &
                     iun, dateo, deet, npas, ni_ext, 1, 1, ip1,  &
                     ip2, ip3, typvar, '>>', etiket, grtyp, ig1,    &
                     ig2, ig3, ig4, datyp, .true.)

    ier = utl_fstecr(lat_ext*MPC_DEGREES_PER_RADIAN_R8, npak, &
                     iun, dateo, deet, npas, 1, nj_ext, 1, ip1,  &
                     ip2, ip3, typvar, '^^', etiket, grtyp, ig1,    &
                     ig2, ig3, ig4, datyp, .true.)

    deallocate(lon_ext)
    deallocate(lat_ext)

    !- 3.2 2D Field
    allocate(Field2d(ni_ext,nj_ext))
    Field2d(:,:) = 10.d0

    deet      =  0
    ip1       =  0
    ip2       =  0
    ip3       =  0
    npas      =  0
    datyp     =  1
    grtyp     =  hco_ens%grtyp
    typvar    = 'A'
    etiket    = 'ANALYSIS'
    dateo  =  0
    ig1       =  hco_ens%ig1 + 100 ! Must be different from the core grid
    ig2       =  hco_ens%ig2 + 100 ! Must be different from the core grid
    ig3       =  0
    ig4       =  0

    ier = utl_fstecr(Field2d, npak,                                  &
                     iun, dateo, deet, npas, ni_ext, nj_ext, 1, ip1, &
                     ip2, ip3, typvar, 'P0', etiket, grtyp, ig1,     &
                     ig2, ig3, ig4, datyp, .true.)

    deallocate(Field2d)

    !
    !- 4.  Closing the output template file
    !
    ier = fstfrm(iun)
    ier = fclos (iun)

  end subroutine CreateLamTemplateGrids

!--------------------------------------------------------------------------
! VarLevInfo
!--------------------------------------------------------------------------
  subroutine VarLevInfo( nlev, GridType, VarName, infile )
    implicit none

    character(len=4), intent(in)  :: VarName
    character(len=*), intent(in)  :: infile

    integer         , intent(out) :: nlev
    character(len=2), intent(out) :: GridType

    integer :: fstouv, fnom, fstfrm, fclos, iunens
    integer :: key, fstinf, fstinl, fstprm, infon
    integer, parameter :: nmax=2000
    integer :: liste(nmax)

    integer           :: ip1, ip2, ip3, ier
    integer           :: ni_t, nj_t, nlev_t, dateo, deet, npas, nbits, datyp
    integer           :: ig1, ig2, ig3, ig4, extra1, extra2, extra3, swa, lng, dltf, ubc
    character(len=4 ) :: nomvar
    character(len=2 ) :: typvar
    character(len=1 ) :: grtyp
    character(len=12) :: etiket

    !
    !- 1.  Open the input file
    !
    iunens = 0
    ier    = fnom(iunens,infile,'RND+OLD+R/O',0)
    ier    = fstouv(iunens,'RND+OLD')

    !
    !- 2.  Find the number of records
    !
    dateo  = -1
    etiket = ' '
    ip1    = -1
    ip2    = ip2_ens
    ip3    = -1
    typvar = ' '

    key = fstinl(iunens,                                       & ! IN
                ni_t, nj_t, nlev_t,                            & ! OUT
                dateo, etiket, ip1, ip2, ip3, typvar, VarName, & ! IN
                liste, infon,                                  & ! OUT
                nmax )                                           ! IN

    if (key < 0) then
      write(*,*)
      write(*,*) 'VarLevInfo: Cannot find variable ', VarName 
      call utl_abort('VarLevInfo') 
    end if

    !
    !- 3.  Find de the Number of verical levels (ip1)
    !
    
    ! We assume that the number of levels = number of records (this could be improve...)
    nlev = infon 
    if ( nlev /= 1 .and. nlev /= vco_bhi%nlev_T .and. nlev /= vco_bhi%nlev_M ) then
       write(*,*)
       write(*,*) 'The number of levels found does not match with the vco structure !!!'
       write(*,*) 'nlev found = ', nlev
       write(*,*) 'vco nlevs  = ', vco_bhi%nlev_T, vco_bhi%nlev_M
       call utl_abort('VarLevInfo')
    end if

    !
    !- 4.  Find the type of vertical grid
    !
    if ( nlev /=1 ) then
      ier = fstprm( liste(1),                                            & ! IN
                    dateo, deet, npas, ni_t, nj_t, nlev_t,               & ! OUT
                    nbits, datyp, ip1, ip2, ip3, typvar, nomvar, etiket, & ! OUT
                    grtyp, ig1, ig2, ig3, ig4, swa, lng, dltf, ubc,      & ! OUT
                    extra1, extra2, extra3 )                               ! OUT

      if      ( ANY(vco_bhi%ip1_M == ip1) ) then
         GridType = 'MM'
      else if ( ANY(vco_bhi%ip1_T == ip1) ) then
         GridType = 'TH'
      else
         write(*,*)
         write(*,*) 'The ip1 found does not match with the vco structure !!!', ip1
         call utl_abort('VarLevInfo')
      end if

    else
       GridType = 'NS'
    end if

    !
    !- 5.  Close the input file
    !
    ier =  fstfrm(iunens)
    ier =  fclos (iunens)

  end subroutine VarLevInfo

!--------------------------------------------------------------------------
! ReadEnsemble
!--------------------------------------------------------------------------
  subroutine ReadEnsemble(ensPerturbations)
    use MathPhysConstants_mod, only: MPC_M_PER_S_PER_KNOT_R4
    implicit none

    real(4), intent(inout) :: ensPerturbations(hco_bhi%ni,hco_bhi%nj,nkgdim,nens)

    real(4), allocatable :: work2d(:,:)

    real(4)   :: factor

    integer   :: ier, fstouv, fnom, fstfrm, fclos, fstecr, fstlir
    integer   :: iunens, ens, var, k, kgdim
    integer   :: ip1, ip2, ip3
    integer   :: ni_t, nj_t, nlev_t, dateo  

    character(len=4 )           :: nomvar
    character(len=1 )           :: typvar
    character(len=12)           :: etiket

    write(*,*)
    write(*,*) 'ReadEnsemble: Starting...'

    allocate(work2d(hco_ens%ni, hco_ens%nj))

    !- Loop over the Ensemble members
    do ens = 1, nens

      write(*,*)
      write(*,*) 'Reading ensemble member: ', ens, trim(cflensin(ens))
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      iunens = 0
      ier    = fnom(iunens,cflensin(ens),'RND+OLD+R/O',0)
      ier    = fstouv(iunens,'RND+OLD')

      !- Loop over Control Variables
      do var = 1, nControlVariable

        !- Loop over vertical Levels
        do k = 1, ControlVariable(var)%nlev

          dateo  = -1
          ip1    = ControlVariable(var)%ip1(k)
          ip2    = ip2_ens
          ip3    = -1
          typvar = ' '
          nomvar = trim(ControlVariable(var)%nomvar(cv_model))
          etiket = ' '
          if ( trim(nomvar) == 'UU' .or. trim(nomvar) == 'VV') then
            factor = MPC_M_PER_S_PER_KNOT_R4 ! knots -> m/s
          else if ( trim(nomvar) == 'P0' ) then
            factor = 100.0 ! hPa -> Pa
          else if ( trim(nomvar) == 'TG' .and. SetTGtoZero ) then
            nomvar = 'TT'
            ip1    = 93423264
            factor = 0.0
          else
            factor = 1.0
          end if

          !- Reading 
          ier = fstlir( work2d,                                      & ! OUT 
                        iunens,                                      & ! IN
                        ni_t, nj_t, nlev_t,                          & ! OUT
                        dateo, etiket, ip1, ip2, ip3, typvar, nomvar)  ! IN 

          if (ier < 0) then
            write(*,*)
            write(*,*) 'ReadEnsemble: Cannot find field ...'
            write(*,*) 'nomvar =', trim(nomvar)
            write(*,*) 'etiket =', trim(etiket)
            write(*,*) 'ip1    =', ip1
            write(*,*) 'ip2    =', ip2
            call utl_abort('ReadEnsemble')
          end if

          if (ni_t /= hco_ens%ni .or. nj_t /= hco_ens%nj) then
            write(*,*)
            write(*,*) 'ReadEnsemble: Invalid dimensions for ...'
            write(*,*) 'nomvar      =', trim(nomvar)
            write(*,*) 'etiket      =', trim(etiket)
            write(*,*) 'ip1         =', ip1
            write(*,*) 'Found ni,nj =', ni_t, nj_t 
            write(*,*) 'Should be   =', hco_ens%ni, hco_ens%nj
            call utl_abort('ReadEnsemble')
          end if

          !- Insert in EnsPerturbations
          kgdim = ControlVariable(var)%kDimStart + k - 1
          EnsPerturbations(1:hco_ens%ni, 1:hco_ens%nj,kgdim,ens) = factor * work2d(:,:) 

        end do ! Vertical Levels

      end do ! Variables
    
      ier =  fstfrm(iunens)
      ier =  fclos (iunens)

    end do ! Ensemble members

    deallocate(work2d)

    write(*,*)
    write(*,*) 'ReadEnsemble: Done!'

  end subroutine ReadEnsemble

!--------------------------------------------------------------------------
! UVToCtrlVar
!--------------------------------------------------------------------------
  subroutine UVToCtrlVar(ensPerturbations)
    implicit none

    real(4), intent(inout) :: ensPerturbations(hco_bhi%ni,hco_bhi%nj,nkgdim,nens)

    real(8), allocatable :: uwind_3d(:,:,:)
    real(8), allocatable :: vwind_3d(:,:,:)
    real(8), allocatable :: CtrlVar1_3d(:,:,:)
    real(8), allocatable :: CtrlVar2_3d(:,:,:)

    integer :: nlev, ens, kStart_u, kStart_v, kEnd_u, kEnd_v

    write(*,*)
    write(*,*) 'UVToCtrlVar: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !- Error traps
    if ( ControlVariable(1)%nlev /= ControlVariable(2)%nlev ) then
      write(*,*)
      write(*,*) 'UVToCtrlVar: Error in Wind Field vertical grid !!!'
      call utl_abort('UVToCtrlVar')
    end if
    if ( ControlVariable(1)%nomvar(cv_model) /= 'UU' .or. &
         ControlVariable(2)%nomvar(cv_model) /= 'VV' ) then
      write(*,*)
      write(*,*) 'UVToCtrlVar: Wind field(s) missing on input !!!'
      write(*,*) 'nomvar(1) = ', ControlVariable(1)%nomvar(cv_model)
      write(*,*) 'nomvar(2) = ', ControlVariable(2)%nomvar(cv_model)
      call utl_abort('UVToCtrlVar')
    end if

    if (trim(WindTransform) == 'UV') then
       write(*,*)
       write(*,*) 'WindTransform = UV. Nothing to do!'
       return
    end if

    !- Allocation
    nlev = ControlVariable(1)%nlev
    allocate(uwind_3d   (hco_bhi%ni, hco_bhi%nj, nlev))
    allocate(vwind_3d   (hco_bhi%ni, hco_bhi%nj, nlev))
    allocate(CtrlVar1_3d(hco_bhi%ni, hco_bhi%nj, nlev))
    allocate(CtrlVar2_3d(hco_bhi%ni, hco_bhi%nj, nlev))

    !- Loop over all ensemble members
    do ens   = 1, nens

      !- Wind position in EnsPerturbations
      kStart_u = ControlVariable(1)%kDimStart
      kStart_v = ControlVariable(2)%kDimStart
      kEnd_u   = ControlVariable(1)%kDimEnd
      kEnd_v   = ControlVariable(2)%kDimEnd

      !- Extract from EnsPerturbations
      uwind_3d(:,:,:) = real(ensPerturbations(:,:,kStart_u:kEnd_u,ens),8)
      vwind_3d(:,:,:) = real(ensPerturbations(:,:,kStart_v:kEnd_v,ens),8)

      !- U-wind,V-wind -> Vorticity,Divergence
      call agd_UVToVortDiv(CtrlVar1_3d, CtrlVar2_3d, & ! OUT
                           uwind_3d   , vwind_3d   , & ! IN
                           nlev                      ) ! IN

      if (trim(WindTransform) == 'PsiChi') then
        !- Vorticity,Divergence -> Stream Function,Velocity Potential
        call VortDivToPsiChi(CtrlVar1_3d, CtrlVar2_3d, & ! INOUT
                             nlev)                       ! IN
      end if

      !- Insert results in EnsPerturbations
      ensPerturbations(:,:,kStart_u:kEnd_u,ens) = real(CtrlVar1_3d(:,:,:),4)
      ensPerturbations(:,:,kStart_v:kEnd_v,ens) = real(CtrlVar2_3d(:,:,:),4)

    end do

    deallocate(uwind_3d   )
    deallocate(vwind_3d   )
    deallocate(CtrlVar1_3d)
    deallocate(CtrlVar2_3d)

    write(*,*)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'UVToCtrlVar: Done!'

  end subroutine UVToCtrlVar

!--------------------------------------------------------------------------
! VortDivToPsiChi
!--------------------------------------------------------------------------
  subroutine VortDivToPsiChi(VortPsi,DivChi,nk)
    implicit none

    integer, intent(in)    :: nk
    real(8), intent(inout) :: VortPsi(hco_bhi%ni, hco_bhi%nj,nk)
    real(8), intent(inout) :: DivChi (hco_bhi%ni, hco_bhi%nj,nk)

    ! Vort -> Psi
    call lst_Laplacian( lst_bhi%id,     & ! IN
                        VortPsi,        & ! INOUT
                        'Inverse', nk )   ! IN    

    ! Div -> Chi
    call lst_Laplacian( lst_bhi%id,     & ! IN
                        DivChi,         & ! INOUT
                        'Inverse', nk )   ! IN

  end subroutine VortDivToPsiChi

!--------------------------------------------------------------------------
! RemoveMean
!--------------------------------------------------------------------------
  subroutine RemoveMean(ensPerturbations,ensMean3d)
    implicit none

    real(4), intent(inout)  :: ensPerturbations(hco_bhi%ni,hco_bhi%nj,nkgdim,nens)
    real(8), intent(out)    :: ensMean3d(hco_bhi%ni,hco_bhi%nj,nkgdim)

    real(8) :: inens
    integer :: i,j,kgdim,ens

    write(*,*)
    write(*,*) 'RemoveMean: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    inens=1.0d0/real(nens,8)

    ensMean3d(:,:,:) = 0.0d0

!$OMP PARALLEL
!$OMP DO PRIVATE (kgdim,ens,j,i)
    do kgdim = 1, nkgdim

      !- Sum
      do ens = 1, nens
        do j = 1 , hco_bhi%nj
          do i = 1, hco_bhi%ni
            ensMean3d(i,j,kgdim) = ensMean3d(i,j,kgdim) + real(ensPerturbations(i,j,kgdim,ens),8)
          end do
        end do
      end do

      !- Mean
      do j = 1, hco_bhi%nj
        do i = 1, hco_bhi%ni
          ensMean3d(i,j,kgdim) = ensMean3d(i,j,kgdim) * inens
        end do
      end do

      !- Remove Mean
      do ens = 1, nens
        do j = 1, hco_bhi%nj
          do i = 1, hco_bhi%ni
            ensPerturbations(i,j,kgdim,ens) = ensPerturbations(i,j,kgdim,ens) - &
                                                real(ensMean3d(i,j,kgdim),4)
          end do
        end do
      end do

    end do
!$OMP END DO
!$OMP END PARALLEL

    write(*,*)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'RemoveMean: Done!'

  end subroutine removeMean

!--------------------------------------------------------------------------
! REMOVEDOMAINMEAN
!--------------------------------------------------------------------------
  subroutine removeDomainMean(ensPerturbations)
    implicit none

    real(4), intent(inout)  :: ensPerturbations(hco_bhi%ni,hco_bhi%nj,nkgdim,nens)

    real(8) :: domainMean, iSize

    integer :: i,j,kgdim,ens

    write(*,*) 'RemoveDomainMean: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    iSize= 1.d0/(dble(hco_bhi%ni)*dble(hco_bhi%nj))

!$OMP PARALLEL
!$OMP DO PRIVATE (kgdim,ens,j,i,domainMean)
    do ens = 1,nens
      do kgdim = 1,nkgdim
        domainMean=0.0d0
        do j = 1,hco_bhi%nj
          do i = 1,hco_bhi%ni
            domainMean = domainMean + real(ensPerturbations(i,j,kgdim,ens),8)
          end do
        end do
        domainMean = domainMean * iSize
        do j = 1,hco_bhi%nj
          do i = 1,hco_bhi%ni
            ensPerturbations(i,j,kgdim,ens) = ensPerturbations(i,j,kgdim,ens) - real(domainMean,4)
          end do
        end do
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    write(*,*)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'RemoveDomainMean: Done!'

  end subroutine removeDomainMean

!--------------------------------------------------------------------------
! CalcStdDev3d
!--------------------------------------------------------------------------
  subroutine CalcStdDev3d(ensPerturbations,StdDev3d,Bin2d,NumBins2d)
    implicit none

    real(4), intent(inout)  :: ensPerturbations(hco_bhi%ni,hco_bhi%nj,nkgdim,nens)
    real(8), intent(out)    :: StdDev3d(hco_bhi%ni,hco_bhi%nj,nkgdim)
    integer, intent(in)     :: Bin2d(hco_bhi%ni,hco_bhi%nj)
    integer, intent(in)     :: NumBins2d

    integer :: BinCount (NumBins2d)
    real(8) :: StdDevBin(NumBins2d)

    real(8) :: inens
    integer :: i,j,kgdim,ens,b

    write(*,*)
    write(*,*) 'CalcStdDev3d: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

!$OMP PARALLEL
!$OMP DO PRIVATE (kgdim,BinCount,StdDevBin,b,ens,j,i)
    do kgdim = 1, nkgdim

      !- Sum of Squares per bin
      BinCount (:) = 0
      StdDevBin(:) = 0.d0
      do ens = 1, nens
        do j = 1, hco_bhi%nj
          do i = 1, hco_bhi%ni
            b = Bin2d(i,j)
            BinCount(b)  = BinCount(b) + 1
            StdDevBin(b) = StdDevBin(b) + real(ensPerturbations(i,j,kgdim,ens),8)**2
          end do
        end do
      end do

      !- Convert to Standard Deviation
      StdDevBin(:) = sqrt( StdDevBin(:) / real((BinCount(:)-1),8) )

      !- Distribute bin values at each grid point
      do j = 1, hco_bhi%nj
        do i = 1, hco_bhi%ni
          b = Bin2d(i,j)
          StdDev3d(i,j,kgdim) = StdDevBin(b)
        end do
      end do

    end do
!$OMP END DO
!$OMP END PARALLEL

    write(*,*)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'CalcStdDev3d: Done!'

  end subroutine CalcStdDev3d

!--------------------------------------------------------------------------
! Normalize3d
!--------------------------------------------------------------------------
  subroutine normalize3d(ensPerturbations,StdDev3d)
    implicit none

    real(4), intent(inout) :: ensPerturbations(hco_bhi%ni,hco_bhi%nj,nkgdim,nens)
    real(8), intent(in)    :: StdDev3d(hco_bhi%ni,hco_bhi%nj,nkgdim)

    integer :: ens, kgdim, j, i

    real(4) :: fact

    write(*,*)
    write(*,*) 'normalize3d: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

!$OMP PARALLEL
!$OMP DO PRIVATE (ens, kgdim, j, i, fact)
    do kgdim = 1, nkgdim
      do j = 1, hco_bhi%nj
         do i = 1, hco_bhi%ni
           if (StdDev3d(i,j,kgdim) > 0.0d0 ) then
             fact = real(1.0d0/StdDev3d(i,j,kgdim),4)
           else
             fact = 0.0
           endif
           do ens = 1, nens
              ensPerturbations(i,j,kgdim,ens) = ensPerturbations(i,j,kgdim,ens) &
                                                * fact
           end do
         end do
       end do
     end do
!$OMP END DO
!$OMP END PARALLEL

    write(*,*)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'normalize3d: Done!'
  
  end subroutine normalize3d

!--------------------------------------------------------------------------
! CalcSpectralStats
!--------------------------------------------------------------------------
  subroutine CalcSpectralStats(ensPerturbations,SpVertCorrel,PowerSpectrum, &
                               NormB)
    implicit none

    real(4), intent(in)     :: ensPerturbations(hco_bhi%ni,hco_bhi%nj,nkgdim,nens)
    real(8), intent(out)    :: SpVertCorrel(nkgdim,nkgdim,0:ntrunc)
    real(8), intent(out)    :: PowerSpectrum(nkgdim,0:ntrunc)
    real(8), intent(out)    :: NormB(nkgdim,nkgdim,0:ntrunc)

    real(8), allocatable    :: NormPowerSpectrum(:,:)
    real(8), allocatable    :: SpectralStateVar(:,:,:)
    real(8), allocatable    :: GridState(:,:,:)
    real(8), allocatable    :: SumWeight(:)

    real(8)           :: weight, sum

    integer           :: i, j, k1, k2, ens, b, e, ila, p, k, totwvnb

    character(len=24) :: kind

    write(*,*)
    write(*,*) 'CalcSpectralStats: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    SpVertCorrel(:,:,:) = 0.d0

    !
    !- 1.  Calculate the Vertical Covariances in Spectral Space
    !
    allocate( SpectralStateVar(lst_bhi%nla,lst_bhi%nphase,nkgdim) )
    allocate( GridState(hco_bhi%ni, hco_bhi%nj, nkgdim) )

    allocate(SumWeight(0:ntrunc))
    SumWeight(:) = 0.d0

    do ens = 1, nens

      !- 1.1 Extract fields from ensPerturbations
      GridState(:,:,:) = real(ensPerturbations(:,:,:,ens),8)

      !- 1.2 Grid Point Space -> Spectral Space
      kind = 'GridPointToSpectral'
      call lst_VarTransform( lst_bhi%id,          & ! IN
                             SpectralStateVar,      & ! OUT
                             GridState,             & ! IN
                             kind, nkgdim     )       ! IN

      !- 1.3 Compute the covariances
!$OMP PARALLEL
!$OMP DO PRIVATE (totwvnb,weight,e,ila,p,k2,k1)
      do totwvnb = 0, ntrunc
        do e = 1, lst_bhi%nePerK(totwvnb)
          ila = lst_bhi%ilaFromEK(e,totwvnb)
          do p = 1, lst_bhi%nphase
            if (trim(SpectralWeights) == 'lst') then
               weight = lst_bhi%Weight(ila)
               SumWeight(totwvnb) = SumWeight(totwvnb) + weight
            else
               weight = 2.0d0
            end if
            do k2 = 1, nkgdim
              do k1 = 1, nkgdim
                SpVertCorrel(k1,k2,totwvnb) = SpVertCorrel(k1,k2,totwvnb) &
                    + weight * (SpectralStateVar(ila,p,k1) * SpectralStateVar(ila,p,k2))
              end do
            end do
          end do
        end do
      end do
!$OMP END DO
!$OMP END PARALLEL

    end do ! Loop in Ensemble

    deallocate(SpectralStateVar)
    deallocate(GridState)

    !- 1.4 Compute the weighted COVARIANCES for each total wavenumber
    do totwvnb = 0, ntrunc
      if (trim(SpectralWeights) == 'legacy') then
         if (totwvnb /= 0 ) then
            SumWeight(totwvnb) = 2.d0 * real(lst_bhi % nphase * lst_bhi % nePerK(totwvnb),8) - 2.d0
         else
            SumWeight(totwvnb) = 1.d0
         end if
         SumWeight(totwvnb) = SumWeight(totwvnb)*nens
      end if

      if ( SumWeight(totwvnb) /= 0.d0 ) then 
         SpVertCorrel(:,:,totwvnb) = SpVertCorrel(:,:,totwvnb) / SumWeight(totwvnb)
      else
         SpVertCorrel(:,:,totwvnb) = 0.d0
      end if

    end do

    deallocate(SumWeight)

    !- 1.5 Extract the power spectrum (the variances on the diagonal elements)
    do k = 1, nkgdim
      PowerSpectrum(k,:) = SpVertCorrel(k,k,:)
    end do

    !
    !- 2.  Calculate the Vertical Correlations in Spectral Space
    !
!$OMP PARALLEL
!$OMP DO PRIVATE (totwvnb,k2,k1)
    do totwvnb = 0, ntrunc
      do k2 = 1, nkgdim
        do k1 = 1, nkgdim 
          if ( PowerSpectrum(k1,totwvnb) /= 0.d0 .and. &
               PowerSpectrum(k2,totwvnb) /= 0.d0 ) then
            SpVertCorrel(k1,k2,totwvnb) = SpVertCorrel(k1,k2,totwvnb) / &
               sqrt( PowerSpectrum(k1,totwvnb) * PowerSpectrum(k2,totwvnb) )
          else
            SpVertCorrel(k1,k2,totwvnb) = 0.d0
          end if
        end do
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    ! Apply vertical localization (if wanted)
    if (vertLoc) then                                                                 
       call applyVertLoc(SpVertCorrel) ! INOUT                                        
    end if

    !
    !- 3.  Normalize the power spectrum (i.e. build normalised spectral densities of the variance)
    !
    allocate(NormPowerSpectrum(nkgdim,0:ntrunc))

    call NormalizePowerSpectrum(PowerSpectrum,     & ! IN
                                NormPowerSpectrum)   ! OUT

    ! Apply horizontal localization (if wanted)
    if (horizLoc) then                                                                 
       call applyHorizLoc(NormPowerSpectrum) ! INOUT                                        
    end if

    !
    !- 4.  Normalize the spectral vertical correlation matrix to ensure correlations in horizontal
    !

!$OMP PARALLEL
!$OMP DO PRIVATE (totwvnb,k2,k1)
    do totwvnb = 0, ntrunc
      do k2 = 1, nkgdim
        do k1 = 1, nkgdim
          NormB(k1,k2,totwvnb) = SpVertCorrel(k1,k2,totwvnb) * &
                sqrt( NormPowerSpectrum(k1,totwvnb) * NormPowerSpectrum(k2,totwvnb) )
        end do
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    deallocate(NormPowerSpectrum)

    write(*,*)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'CalcSpectralStats: Done!'

  end subroutine CalcSpectralStats

!--------------------------------------------------------------------------
! NormalizePowerSpectrum
!--------------------------------------------------------------------------
  subroutine NormalizePowerSpectrum(PowerSpectrum, NormPowerSpectrum)
    implicit none

    real(8), intent(in)    :: PowerSpectrum(nkgdim,0:ntrunc)
    real(8), intent(out)   :: NormPowerSpectrum(nkgdim,0:ntrunc)

    real(8), allocatable   :: SpectralStateVar(:,:,:)
    real(8), allocatable   :: GridState(:,:,:)

    real(8)           :: sum

    integer           :: i, j, e, ila, p, k, totwvnb

    character(len=24) :: kind

    write(*,*)
    write(*,*) 'NormalizePowerSpectrum: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- 1.  Normalize the power spectrum (i.e. build normalised spectral densities of the variance)
    !

    !- 1.1 Part 1

!$OMP PARALLEL
!$OMP DO PRIVATE (totwvnb,k,sum)
    do k = 1, nkgdim
      sum = 0.0d0
      do totwvnb = 0, ntrunc
        sum = sum + real(totwvnb,8) * PowerSpectrum(k,totwvnb)
      end do
      do totwvnb = 0, ntrunc
        if ( sum /= 0.0d0 ) then
          NormPowerSpectrum(k,totwvnb) = PowerSpectrum(k,totwvnb) / sum
        else
          NormPowerSpectrum(k,totwvnb) = 0.d0
        end if
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    !- 1.2 Part 2
    allocate( SpectralStateVar(lst_bhi%nla,lst_bhi%nphase,nkgdim) )
    allocate( GridState(hco_bhi%ni, hco_bhi%nj, nkgdim) )

    !- 1.2.1 Spectral transform of a delta function (at the center of the domain)
    GridState(:,:,:) = 0.d0
    GridState(hco_bhi%ni/2,hco_bhi%nj/2,:) = 1.d0

    kind = 'GridPointToSpectral'
    call lst_VarTransform( lst_bhi%id,          & ! IN
                           SpectralStateVar,      & ! OUT
                           GridState,             & ! IN
                           kind, nkgdim     )       ! IN

    !- 1.2.2 Apply the horizontal correlation function
!$OMP PARALLEL
!$OMP DO PRIVATE (totwvnb,e,ila,p,k)
    do totwvnb = 0, ntrunc
       do e = 1, lst_bhi%nePerK(totwvnb)
          ila = lst_bhi%ilaFromEK(e,totwvnb)
          do p = 1, lst_bhi%nphase
             do k = 1, nkgdim
                SpectralStateVar(ila,p,k) = SpectralStateVar(ila,p,k) * NormPowerSpectrum(k,totwvnb) * &
                              lst_bhi%NormFactor(ila,p) * lst_bhi%NormFactorAd(ila,p)
             end do
          end do
       end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    !- 1.2.3 Move back to physical space
    kind = 'SpectralToGridPoint'
    call lst_VarTransform( lst_bhi%id,      & ! IN
                           SpectralStateVar,  & ! IN
                           GridState,         & ! OUT
                           kind, nkgdim )       ! IN

    !- 1.2.4 Normalize to 1
    do k = 1, nkgdim
       if ( GridState(hco_bhi%ni/2,hco_bhi%nj/2,k) < 0.d0 ) then
          write(*,*) 'NormalizePowerSpectrum: Problem in normalization ', k, GridState(hco_bhi%ni/2,hco_bhi%nj/2,k)
          call utl_abort('aborting in NormalizePowerSpectrum')
       end if

       if ( GridState(hco_bhi%ni/2,hco_bhi%nj/2,k) /= 0.d0 ) then
          write(*,*) 'Normalization factor = ', k, GridState(hco_bhi%ni/2,hco_bhi%nj/2,k), 1.d0 / GridState(hco_bhi%ni/2,hco_bhi%nj/2,k)
          NormPowerSpectrum(k,:) = NormPowerSpectrum(k,:) / GridState(hco_bhi%ni/2,hco_bhi%nj/2,k)
       else
          write(*,*) 'Setting NormPowerSpectrum to zero = ', k
          NormPowerSpectrum(k,:) = 0.d0
       end if
    end do

    deallocate(SpectralStateVar)
    deallocate(GridState)

    write(*,*)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'NormalizePowerSpectrum: Done!'

  end subroutine NormalizePowerSpectrum

!--------------------------------------------------------------------------
! CalcHorizScale
!--------------------------------------------------------------------------
  subroutine CalcHorizScale(HorizScale,SpCovariance)
    use MathPhysConstants_mod, only: MPC_PI_R8
    use EarthConstants_mod, only: RA
    implicit none

    real(8), intent(out) :: HorizScale(nkgdim)
    real(8), intent(in)  :: SpCovariance(nkgdim,nkgdim,0:ntrunc)

    real(8) :: circ_eq, cur_circ_eq, un_deg_lon, dx, dist
    real(8) :: a, b, beta

    integer :: totwvnb, k, var

    write(*,*)
    write(*,*) 'CalcHorizScale: Starting...'

    !
    !- Computing distance-related variables
    !

    ! Grid spacing in meters
    dx = hco_bhi%dlon * RA
    write(*,*)
    write(*,*) 'grid spacing (m) =', dx

    dist = max(hco_bhi%ni, hco_bhi %nj) * dx
    beta = (dist/(2.d0*MPC_PI_R8))**2

    !
    !- Estimate horizontal correlation scales based on the power spectra
    !
    do k = 1, nkgdim
      a = 0.d0
      b = 0.d0
      do totwvnb = 0, ntrunc
        a = a + SpCovariance(k,k,totwvnb) * totwvnb
        b = b + SpCovariance(k,k,totwvnb) * totwvnb**3
      end do
      if (b <= 0.d0) then 
        HorizScale(k) = 0.d0
      else
        HorizScale(k) = sqrt(2.d0*a*beta/b)
      end if
    end do

    do var = 1, nControlVariable
      write(*,*)
      write(*,*) ControlVariable(var)%nomvar(cv_bhi)
      do k = ControlVariable(var)%kDimStart, ControlVariable(var)%kDimEnd
        write(*,'(i3,2X,f9.2,2X,a2)') k, HorizScale(k)/1000.d0, 'km'
      end do
    end do

    write(*,*)
    write(*,*) 'CalcHorizScale: Done!'

  end subroutine CalcHorizScale

!--------------------------------------------------------------------------
! CalcTotVertCorrel
!--------------------------------------------------------------------------
  subroutine CalcTotVertCorrel(TotVertCorrel, SpVertCorrel, PowerSpectrum)
    implicit none

    real(8), intent(out)    :: TotVertCorrel(nkgdim,nkgdim)
    real(8), intent(in)     :: SpVertCorrel(nkgdim,nkgdim,0:ntrunc)
    real(8), intent(in)     :: PowerSpectrum(nkgdim,0:ntrunc)

    real(8), allocatable    :: TotVertCov(:,:)

    integer           :: k1, k2, totwvnb

    write(*,*)
    write(*,*) 'CalcTotVertCorrel: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- 1.  Calculate the total Normalized Covariance Matrix
    !
    allocate(TotVertCov(nkgdim,nkgdim))
    TotVertCov(:,:) = 0.d0

    do k2 = 1, nkgdim
      do k1 = 1, nkgdim
        do totwvnb = 0, ntrunc
          TotVertCov(k1,k2) = TotVertCov(k1,k2) + &
                              real(totwvnb,8) * SpVertCorrel(k1,k2,totwvnb) *    &
                              sqrt(PowerSpectrum(k1,totwvnb) * PowerSpectrum(k2,totwvnb))
        end do
      end do
    end do

    !
    !- 2.  Transform into correlations
    !
    do k2 = 1, nkgdim
      do k1 = 1, nkgdim
        if ( TotVertCov(k1,k1) /= 0.d0 .and. &
             TotVertCov(k2,k2) /= 0.d0 ) then
          TotVertCorrel(k1,k2) = TotVertCov(k1,k2) / &
               ( sqrt(TotVertCov(k1,k1)) * sqrt(TotVertCov(k2,k2)) )
        else
          TotVertCorrel(k1,k2) = 0.d0
        end if
      end do
    end do

    deallocate(TotVertCov)

    write(*,*)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'CalcTotVertCorrel: Done!'

  end subroutine CalcTotVertCorrel

!--------------------------------------------------------------------------
! CalcBsqrt
!--------------------------------------------------------------------------
  subroutine CalcBsqrt(Bsqrt,B)
    implicit none
    !
    !  - produce matrix B^0.5 = V D^0.5 D^t where V and D are the
    !    eigenvectors and eigenvalues of B
    !
    real(8), intent(out)   :: Bsqrt(nkgdim,nkgdim,0:ntrunc)
    real(8), intent(in)    :: B    (nkgdim,nkgdim,0:ntrunc)

    real(8), allocatable :: EigenValues(:)
    real(8), allocatable :: Work(:)

    real(8), allocatable :: EigenVectors(:,:)

    integer :: sizework, info, totwvnb, k, k1, k2 

    write(*,*)
    write(*,*) 'CalcBsqrt: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    sizework = 64 * nkgdim
    allocate(work(sizework))

    allocate(EigenValues(nkgdim))
    allocate(EigenVectors(nkgdim,nkgdim))

    !
    !-  Calculate B^0.5 for each total wave number
    !
    do totwvnb = 0, ntrunc

       EigenVectors(:,:) = B(:,:,totwvnb)

       !- Calculate EigenVectors (V) and EigenValues (D) of B matrix
       call dsyev('V','U',nkgdim,  & ! IN
                   EigenVectors,   & ! INOUT
                   nkgdim,         & ! IN
                   EigenValues,    & ! OUT
                   work, sizework, & ! IN
                   info )            ! OUT
       
       if ( info /= 0 ) then
         write(*,*)
         write(*,*) 'CalcBsqrt: DSYEV failed !!! ', totwvnb, info
         call utl_abort('CalcBsqrt')
       end if

       !- Calculate B^0.5 = V D^0.5 V^t
       where(EigenValues < 0.d0)
         EigenValues = 0.d0
       end where

       do k1 = 1, nkgdim
         do k2 = 1, nkgdim
           Bsqrt(k1,k2,totwvnb) = sum ( EigenVectors (k1,1:nkgdim)   &
                                   *    EigenVectors (k2,1:nkgdim)   &
                                   *    sqrt(EigenValues(1:nkgdim)) )
         end do
       end do

    end do  ! total wave number

    deallocate(EigenVectors)
    deallocate(EigenValues)
    deallocate(work)

    write(*,*)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'CalcBsqrt: Done!'

  end subroutine CalcBsqrt

!--------------------------------------------------------------------------
! SetSpVertCorrel
!--------------------------------------------------------------------------
  subroutine SetSpVertCorrel(SpVertCorrel)
    implicit none

    real(8), intent(inout) :: SpVertCorrel(nkgdim,nkgdim,0:ntrunc)

    real(8), allocatable :: KeepOrDiscard(:,:)

    integer :: totwvnb, var1, var2, k1, k2 

    write(*,*)
    write(*,*) 'SetSpVertCorrel: Starting...'

    !
    !-  Determine which bloc of the correlation matrix to keep/discard
    !
    allocate(KeepOrDiscard(nControlVariable,nControlVariable))

    !- Calculate upper half
    do var2 = 1, nControlVariable
      do var1 = 1, var2
         if (var1 == var2) then
           KeepOrDiscard(var1,var2) = 1.d0 ! Keep Auto-Correlations
         elseif( (ControlVariable(var2)%nomvar(cv_bhi) == momentumControlVar(2) .and. & ! e.g. PP-CC
                  ControlVariable(var1)%nomvar(cv_bhi) == momentumControlVar(1)) .or. &
                 (ControlVariable(var2)%nomvar(cv_bhi) == 'TT'                  .and. & ! TT-(PP, QR or UU)
                  ControlVariable(var1)%nomvar(cv_bhi) == momentumControlVar(1)) .or. &
                 (ControlVariable(var2)%nomvar(cv_bhi) == 'LQ'                  .and. & ! LQ-(PP, QR or UU)
                  ControlVariable(var1)%nomvar(cv_bhi) == momentumControlVar(1) .and. &
                  correlationLQ ) .or. &
                 (ControlVariable(var2)%nomvar(cv_bhi) == 'P0'                  .and. & ! P0-(PP, QR or UU)
                  ControlVariable(var1)%nomvar(cv_bhi) == momentumControlVar(1)) .or. &
                 (ControlVariable(var2)%nomvar(cv_bhi) == 'P0'                  .and. & ! P0-TT
                  ControlVariable(var1)%nomvar(cv_bhi) == 'TT')                  .or. &
                 (ControlVariable(var2)%nomvar(cv_bhi) == 'P0'                  .and. & ! P0-LQ
                  ControlVariable(var1)%nomvar(cv_bhi) == 'LQ'                  .and. &
                  correlationLQ ) .or. &
                 (ControlVariable(var2)%nomvar(cv_bhi) == 'LQ'                  .and. & ! LQ-TT
                  ControlVariable(var1)%nomvar(cv_bhi) == 'TT'                  .and. &
                  correlationLQ ) .or. &
                 (ControlVariable(var2)%nomvar(cv_bhi) == 'TT'                  .and. & ! TT-(CC, DD or VV)
                  ControlVariable(var1)%nomvar(cv_bhi) == momentumControlVar(2)) .or. &
                 (ControlVariable(var2)%nomvar(cv_bhi) == 'LQ'                  .and. & ! LQ-VV
                  ControlVariable(var1)%nomvar(cv_bhi) == momentumControlVar(2) .and. &
                  correlationLQ ) .or. &
                 (ControlVariable(var2)%nomvar(cv_bhi) == 'P0'                  .and. & ! P0-(CC, DD or VV)
                  ControlVariable(var1)%nomvar(cv_bhi) == momentumControlVar(2)) ) then
           KeepOrDiscard(var1,var2) = 1.d0 ! Keep these Cross-Correlations
         else
           KeepOrDiscard(var1,var2) = 0.d0 ! Discard these Cross-Correlation
         end if
         write(*,*) var1, var2, ControlVariable(var1)%nomvar(cv_bhi), ControlVariable(var2)%nomvar(cv_bhi), KeepOrDiscard(var1,var2) 
      end do
    end do

    ! Symmetrize
    do var1 = 2, nControlVariable
      do var2 = 1, var1-1
        KeepOrDiscard(var1,var2) = KeepOrDiscard(var2,var1)
      end do
    end do

    !
    !- Modify the Vertical Correlation Matrix
    !
    do totwvnb = 0, ntrunc

      do var2 = 1, nControlVariable
        do var1 = 1, nControlVariable

          do k2 = ControlVariable(var2)%kDimStart, ControlVariable(var2)%kDimEnd
            do k1 = ControlVariable(var1)%kDimStart, ControlVariable(var1)%kDimEnd
              SpVertCorrel(k1,k2,totwvnb) = KeepOrDiscard(var1,var2) * &
                                            SpVertCorrel(k1,k2,totwvnb)
            end do
          end do

        end do
      end do

    end do  ! total wave number

    deallocate(KeepOrDiscard)

    write(*,*)
    write(*,*) 'SetSpVertCorrel: Done!'

  end subroutine SetSpVertCorrel

!--------------------------------------------------------------------------
! CreateBins
!--------------------------------------------------------------------------
  subroutine CreateBins(Bin2d,NumBins2d,BinningStrategy)
    implicit none

    integer, intent(out)     :: Bin2d(hco_bhi%ni,hco_bhi%nj)
    integer, intent(out)     :: NumBins2d
    character(*), intent(in) :: BinningStrategy

    real(8) :: BinCat
    integer :: i,j,kgdim,ens,b

    select case(trim(BinningStrategy))
    case ('GridPoint')
      write(*,*) ' BIN_TYPE : No horizontal averaging '
      write(*,*) '          = One bin per horizontal grid point'

      NumBins2d = hco_bhi%ni * hco_bhi%nj
      BinCat    = 0
      do j = 1, hco_bhi%nj
        do i = 1, hco_bhi%ni
          BinCat = BinCat + 1
          bin2d(i,j) = BinCat
        end do
      end do

    case ('YrowBand')
      write(*,*)
      write(*,*) ' BIN_TYPE : One bin per Y row'

      NumBins2d = hco_bhi%nj
      BinCat    = 0
      do j = 1, hco_bhi%nj
        BinCat = BinCat + 1
        bin2d(:,j) = BinCat
      end do

    case ('HorizontalMean','FullDomain')
      write(*,*)
      write(*,*) ' BIN_TYPE : Average over all horizontal points'

      NumBins2d   = 1
      bin2d(:, :) = 1

    case ('CoreExt')
      write(*,*)
      write(*,*) ' BIN_TYPE : Partition between CORE and EXTENSION zone'

      NumBins2d = 2
      do j = 1, hco_bhi%nj
        do i = 1, hco_bhi%ni
           if ( i <= hco_ens%ni .and. j <= hco_ens%nj ) then
              bin2d(i,j) = 1 ! Core
           else
              bin2d(i,j) = 2 ! Extension
           end if
        end do
      end do

    case default
      write(*,*)
      write(*,*) 'Invalid Binning Strategy : ',trim(BinningStrategy)
      call utl_abort('CreateBins')
    end select

  end subroutine CreateBins

!--------------------------------------------------------------------------
! BiPeriodization
!--------------------------------------------------------------------------
  subroutine BiPeriodization(ensPerturbations)
    implicit none

    real(4), intent(inout) :: ensPerturbations(hco_bhi%ni,hco_bhi%nj,nkgdim,nens)

    real(8), allocatable :: work3d(:,:,:)

    integer :: kgdim, ens

    write(*,*)
    write(*,*) 'BiPeriodization: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    allocate(work3d(hco_bhi%ni, hco_bhi%nj, nkgdim))

    !- Loop over all variables / levels / ensemble members
    do ens   = 1, nens

        work3d(:,:,:) = real(ensPerturbations(:,:,:,ens),8)

        call agd_mach(work3d,                             & ! INOUT
                      hco_bhi%ni, hco_bhi%nj, nkgdim)   ! IN

        ensPerturbations(:,:,:,ens) = real(work3d(:,:,:),4)

    end do

    deallocate(work3d)

    write(*,*)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'BiPeriodization: Done!'

  end subroutine BiPeriodization

!--------------------------------------------------------------------------
! CalcVertCorrel
!--------------------------------------------------------------------------
  subroutine CalcVertCorrel(VertCorrel,ensPerturbations,Bin2d,NumBins2d)
    implicit none

    real(4), intent(in)     :: ensPerturbations(hco_bhi%ni,hco_bhi%nj,nkgdim,nens)
    real(8), intent(out)    :: VertCorrel(nkgdim,nkgdim,NumBins2d)
    integer, intent(in)     :: Bin2d(hco_bhi%ni,hco_bhi%nj)
    integer, intent(in)     :: NumBins2d

    integer :: BinCount (NumBins2d)
    integer :: i, j, k1, k2, ens, b

    write(*,*)
    write(*,*) 'CalcVertCorrel: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    VertCorrel(:,:,:) = 0.d0
    BinCount (:) = 0

    !
    !- Calculate the Vertical Correlation in GridPoint Space
    !
    
    ! ... we assume that the ensemble grid point mean was removed and that
    !     the ensemble values were divided by the grid point std dev.

    do ens = 1, nens

       do j = 1, hco_bhi%nj
          do i = 1, hco_bhi%ni

             b = Bin2d(i,j)
             if (ens == 1) BinCount(b) = BinCount(b) + 1

!$OMP PARALLEL
!$OMP DO PRIVATE (k1,k2)
             do k2 = 1, nkgdim
                do k1 = 1, nkgdim
                   VertCorrel(k1,k2,b) = VertCorrel(k1,k2,b) &
                        + real(ensPerturbations(i,j,k1,ens),8)*real(ensPerturbations(i,j,k2,ens),8)
                end do
             end do
!$OMP END DO
!$OMP END PARALLEL

          end do
       end do

    end do ! Loop in Ensemble

    do b = 1, NumBins2d
       VertCorrel(:,:,b) = VertCorrel(:,:,b) / real((nens-1)*BinCount(b),8)
    end do

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'CalcVertCorrel: Done!'

  end subroutine CalcVertCorrel

!--------------------------------------------------------------------------
! HORIZCORRELFUNCTION
!--------------------------------------------------------------------------
  subroutine horizCorrelFunction(NormB)
    implicit none

    real(8), intent(in)    :: NormB(nkgdim,nkgdim,0:ntrunc)

    real(8), allocatable    :: SpectralStateVar(:,:,:)
    real(8), allocatable    :: GridState(:,:,:)

    integer   :: i, j, e, ila, p, k, totwvnb

    integer   :: ier, fstouv, fnom, fstfrm, fclos
    integer   :: iunstats

    character(len=24) :: kind

    allocate( SpectralStateVar(lst_bhi%nla,lst_bhi%nphase,nkgdim) )
    allocate( GridState(hco_bhi%ni, hco_bhi%nj, nkgdim) )

    !- 3.2.1 Spectral transform of a delta function (at the center of the domain)
    GridState(:,:,:) = 0.d0
    GridState(hco_bhi%ni/2,hco_bhi%nj/2,:) = 1.d0

    kind = 'GridPointToSpectral'
    call lst_VarTransform( lst_bhi%id,            & ! IN
                           SpectralStateVar,      & ! OUT
                           GridState,             & ! IN
                           kind, nkgdim     )       ! IN

    !- 3.2.2 Apply the horizontal correlation function
!$OMP PARALLEL
!$OMP DO PRIVATE (totwvnb,e,ila,p,k)
    do totwvnb = 0, ntrunc
       do e = 1, lst_bhi%nePerK(totwvnb)
          ila = lst_bhi%ilaFromEK(e,totwvnb)
          do p = 1, lst_bhi%nphase
             do k = 1, nkgdim
                SpectralStateVar(ila,p,k) = SpectralStateVar(ila,p,k) * NormB(k,k,totwvnb) * &
                              lst_bhi%NormFactor(ila,p) * lst_bhi%NormFactorAd(ila,p)
             end do
          end do
       end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    !- 3.2.3 Move back to physical space
    kind = 'SpectralToGridPoint'
    call lst_VarTransform( lst_bhi%id,        & ! IN
                           SpectralStateVar,  & ! IN
                           GridState,         & ! OUT
                           kind, nkgdim )       ! IN

    !
    !- 4.  Write to file
    !
    iunstats = 0
    ier    = fnom(iunstats,'./horizCorrel.fst','RND',0)
    ier    = fstouv(iunstats,'RND')
    call write3d(GridState,iunstats,'HORIZCORFUNC',cv_bhi) ! IN
    call WriteTicTacToc(iunstats) ! IN
    ier =  fstfrm(iunstats)
    ier =  fclos (iunstats)

    deallocate(SpectralStateVar)
    deallocate(GridState)

  end subroutine horizCorrelFunction

!--------------------------------------------------------------------------
! ApplyHorizLoc
!--------------------------------------------------------------------------
  subroutine applyHorizLoc(NormPowerSpectrum)
    implicit none

    real(8), intent(inout) :: NormPowerSpectrum(nkgdim,0:ntrunc)

    real(8), allocatable   :: SpectralStateVar(:,:,:)
    real(8), allocatable   :: GridState(:,:,:)
    real(8), allocatable   :: GridStateLoc(:,:,:)
    real(8), allocatable   :: PowerSpectrum(:,:)
    real(8), allocatable   :: SumWeight(:)
    real(8), allocatable   :: local_length(:)

    integer :: totwvnb, var, k, i, j, e, ila, p

    real(8)  :: dist, fact, hlocalize

    character(len=24) :: kind

    write(*,*)
    write(*,*) 'applyHorizLoc: Starting...'

    allocate( SpectralStateVar(lst_bhi%nla,lst_bhi%nphase,nkgdim) )
    allocate( GridState(hco_bhi%ni, hco_bhi%nj, nkgdim) )

    !
    !- 1. Get the original gridpoint horizontal correlations
    !

    !- 1.1 Spectral transform of a delta function (at the lower-left of the domain)
    GridState(:,:,:) = 0.d0
    GridState(1,1,:) = 1.d0

    kind = 'GridPointToSpectral'
    call lst_VarTransform( lst_bhi%id,            & ! IN
                           SpectralStateVar,      & ! OUT
                           GridState,             & ! IN
                           kind, nkgdim     )       ! IN

    !- 1.2 Apply the horizontal correlation function
!$OMP PARALLEL
!$OMP DO PRIVATE (totwvnb,e,ila,p,k)
    do totwvnb = 0, ntrunc
       do e = 1, lst_bhi%nePerK(totwvnb)
          ila = lst_bhi%ilaFromEK(e,totwvnb)
          do p = 1, lst_bhi%nphase
             do k = 1, nkgdim
                SpectralStateVar(ila,p,k) = SpectralStateVar(ila,p,k) * NormPowerSpectrum(k,totwvnb) * &
                              lst_bhi%NormFactor(ila,p) * lst_bhi%NormFactorAd(ila,p)
             end do
          end do
       end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    !- 1.3 Move back to physical space
    kind = 'SpectralToGridPoint'
    call lst_VarTransform( lst_bhi%id,        & ! IN
                           SpectralStateVar,  & ! IN
                           GridState,         & ! OUT
                           kind, nkgdim )       ! IN

    !
    !- 2.  Create and apply the localization function in gridpoint space
    !
    allocate(local_length(nkgdim))
    allocate(GridStateLoc(hco_bhi%ni, hco_bhi%nj, nkgdim) )

    do var = 1, nControlVariable
       if      ( ControlVariable(var)%nomvar(cv_bhi) == 'CC' .or. &
                 ControlVariable(var)%nomvar(cv_bhi) == 'PP' .or. &
                 ControlVariable(var)%nomvar(cv_bhi) == 'QR' .or. &
                 ControlVariable(var)%nomvar(cv_bhi) == 'DD' ) then
          hLocalize = hLocalize_wind
       else if ( ControlVariable(var)%nomvar(cv_bhi) == 'TT' .or. &
                 ControlVariable(var)%nomvar(cv_bhi) == 'TG' .or. &
                 ControlVariable(var)%nomvar(cv_bhi) == 'P0' ) then
          hLocalize = hLocalize_mass
       else if ( ControlVariable(var)%nomvar(cv_bhi) == 'LQ' ) then
          hLocalize = hLocalize_humidity
       else
          call utl_abort('applyHorizLoc')
       end if

       do k = ControlVariable(var)%kDimStart, ControlVariable(var)%kDimEnd
          local_length(k) = hLocalize
       end do

    end do

    call lfn_setup('FifthOrder') ! IN
    call lfn_CreateBiPerFunction( GridStateLoc,                  & ! OUT
                                  local_length, hco_bhi%dlon,    & ! IN
                                  hco_bhi%ni, hco_bhi%nj, nkgdim)  ! IN

    GridState(:,:,:) = GridStateLoc(:,:,:) * GridState(:,:,:)

    deallocate(GridStateLoc)
    deallocate(local_length)

    !
    !- 3. Create the localized normalized power spectrum
    !    
    allocate(PowerSpectrum(nkgdim,0:ntrunc))

    !- 3.1 Transform to spectral space
    kind = 'GridPointToSpectral'
    call lst_VarTransform( lst_bhi%id   ,     & ! IN
                           SpectralStateVar,  & ! OUT
                           GridState,         & ! IN
                           kind, nkgdim )       ! IN
 
    !- 3.2 Compute band mean
    allocate(SumWeight(0:ntrunc))
    SumWeight(:) = 0.d0

    PowerSpectrum(:,:) = 0.d0
    do totwvnb = 0, ntrunc
       do e = 1, lst_bhi%nePerK(totwvnb)
          ila = lst_bhi%ilaFromEK(e,totwvnb)
          do p = 1, lst_bhi%nphase
             SumWeight(totwvnb) = SumWeight(totwvnb) + lst_bhi%Weight(ila)
             do k = 1, nkgdim
                PowerSpectrum(k,totwvnb) = PowerSpectrum(k,totwvnb) + &
                     lst_bhi%Weight(ila) * abs(SpectralStateVar(ila,p,k))
             end do
          end do
       end do
    end do
    
    do totwvnb = 0, ntrunc
       if (SumWeight(totwvnb) /= 0.d0) then
          PowerSpectrum(:,totwvnb) = PowerSpectrum(:,totwvnb) / SumWeight(totwvnb)
       else
          PowerSpectrum(:,totwvnb) = 0.d0
       endif
    end do

    deallocate(SumWeight)

    !- 3.3 Normalize
    call NormalizePowerSpectrum(PowerSpectrum,     & ! IN
                                NormPowerSpectrum)   ! OUT


    deallocate(SpectralStateVar)
    deallocate(GridState)
    deallocate(PowerSpectrum)

    write(*,*)
    write(*,*) 'applyHorizLoc: Done!'

  end subroutine applyHorizLoc

!--------------------------------------------------------------------------
! ApplyVertLoc
!--------------------------------------------------------------------------
  subroutine applyVertLoc(SpVertCorrel)
    implicit none

    real(8), intent(inout) :: SpVertCorrel(nkgdim,nkgdim,0:ntrunc)

    integer :: totwvnb, var1, var2, k1, k2, lev1, lev2

    real(8)  :: dist, fact, vLocalize, pres1, pres2, vLocalize1, vLocalize2

    write(*,*)
    write(*,*) 'applyVertLoc: Starting...'

    !
    !- 1.  Select the localization function
    !
    call lfn_setup('FifthOrder') ! IN 

    !
    !- 2.  Apply localization to the spectral vertical correlations
    !

    !- 2.1 Loop on control variables
    do var2 = 1, nControlVariable
       do var1 = 1, nControlVariable

          !-  2.2 Set the vertical length scale

          !-  2.2.1 Select the length scale for control variable 1
          if      ( ControlVariable(var1)%nomvar(cv_bhi) == 'CC' .or. &
                    ControlVariable(var1)%nomvar(cv_bhi) == 'PP' .or. &
                    ControlVariable(var1)%nomvar(cv_bhi) == 'QR' .or. &
                    ControlVariable(var1)%nomvar(cv_bhi) == 'DD') then
             vLocalize1 = vLocalize_wind
          else if ( ControlVariable(var1)%nomvar(cv_bhi) == 'TT' .or. &
                    ControlVariable(var1)%nomvar(cv_bhi) == 'TG' .or. &
                    ControlVariable(var1)%nomvar(cv_bhi) == 'P0' ) then
             vLocalize1 = vLocalize_mass
          else if ( ControlVariable(var1)%nomvar(cv_bhi) == 'LQ' ) then
             vLocalize1 = vLocalize_humidity
          else
             call utl_abort('applyVertLoc')
          end if

          !-  2.2.2 Select the length scale for control variable 2
          if      ( ControlVariable(var2)%nomvar(cv_bhi) == 'CC' .or. &
                    ControlVariable(var2)%nomvar(cv_bhi) == 'PP' .or. &
                    ControlVariable(var2)%nomvar(cv_bhi) == 'QR' .or. &
                    ControlVariable(var2)%nomvar(cv_bhi) == 'DD') then
             vLocalize2 = vLocalize_wind
          else if ( ControlVariable(var2)%nomvar(cv_bhi) == 'TT' .or. &
                    ControlVariable(var2)%nomvar(cv_bhi) == 'TG' .or. &
                    ControlVariable(var2)%nomvar(cv_bhi) == 'P0' ) then
             vLocalize2 = vLocalize_mass
          else if ( ControlVariable(var2)%nomvar(cv_bhi) == 'LQ' ) then
             vLocalize2 = vLocalize_humidity
          else
             call utl_abort('applyVertLoc')
          end if
          
          !- 2.2.3 Length scale to be use for var1-var2 correlation
          vLocalize = (vLocalize1+vLocalize2)/2.d0

          !- 2.3 Loop on vertical levels
          do k2 = ControlVariable(var2)%kDimStart, ControlVariable(var2)%kDimEnd
             do k1 = ControlVariable(var1)%kDimStart, ControlVariable(var1)%kDimEnd

                !- 2.4 Set the pressure values

                !- 2.4.1 Pressure for control-variable-1 level
                if (ControlVariable(var1)%nlev /= 1) then ! variable 3D
                   lev1 = k1 - ControlVariable(var1)%kDimStart + 1
                   if (ControlVariable(var1)%GridType == 'TH') then
                      pres1 = pressureProfile_T(lev1)
                   else
                      pres1 = pressureProfile_M(lev1)
                   end if
                else
                   pres1 = pressureProfile_M(vco_bhi%nlev_M) ! variable 2D
                end if

                !- 2.4.2 Pressure for control-variable-2 level
                if (ControlVariable(var2)%nlev /= 1) then ! variable 3D
                   lev2 = k2 - ControlVariable(var2)%kDimStart + 1
                   if (ControlVariable(var2)%GridType == 'TH') then
                      pres2 = pressureProfile_T(lev2)
                   else
                      pres2 = pressureProfile_M(lev2)
                   end if
                else
                   pres2 = pressureProfile_M(vco_bhi%nlev_M) ! variable 2D
                end if

                !- 2.5 Compute the localization factor
                dist = abs(log(pres2) - log(pres1))
                fact = lfn_response(dist,vLocalize)

                !- 2.6 Localize each total wavenumber (not scale-dependent!)
                do totwvnb = 0, ntrunc
                   SpVertCorrel(k1,k2,totwvnb) = fact * SpVertCorrel(k1,k2,totwvnb)
                end do

             end do
          end do

       end do
    end do

    deallocate(pressureProfile_M)
    deallocate(pressureProfile_T)

    write(*,*)
    write(*,*) 'applyVertLoc: Done!'

  end subroutine applyVertLoc

!--------------------------------------------------------------------------
! WriteVarStats
!--------------------------------------------------------------------------
  subroutine WriteVarStats(Bsqrt,StdDev3d)
    implicit none

    real(8), intent(in) :: Bsqrt(nkgdim,nkgdim,0:ntrunc)
    real(8), intent(in) :: StdDev3d(hco_bhi%ni,hco_bhi%nj,nkgdim)

    integer   :: ier, fstouv, fnom, fstfrm, fclos
    integer   :: iunstats

    write(*,*)
    write(*,*) 'Writing covariance statistics for VAR'

    !
    !- Opening Output file
    !
    iunstats = 0
    ier    = fnom(iunstats,'./bgcov.fst','RND',0)
    ier    = fstouv(iunstats,'RND')

    !
    !- Add Tic-Tac and Toc-Toc
    !
    call WriteTicTacToc(iunstats)

    !
    !- Add Control Variable Info
    !
    call WriteControlVarInfo(iunstats)

    !
    !- Bsqrt
    !
    call WriteSpVertCorrel(Bsqrt,iunstats,'ZN','B_SQUAREROOT') ! IN
 
    !
    !- 3D Standard Deviations
    !
    call Write3d(StdDev3d,iunstats,'STDDEV',cv_bhi) ! IN

    !
    !- Closing output file
    !
    ier =  fstfrm(iunstats)
    ier =  fclos (iunstats)

  end subroutine WriteVarStats

!--------------------------------------------------------------------------
! WriteDiagStats
!--------------------------------------------------------------------------
  subroutine WriteDiagStats(NormB,SpVertCorrel,TotVertCorrel,EnsMean3d, &
                            StdDev3dGridPoint,PowerSpectrum,HorizScale)
    implicit none

    real(8), intent(in) :: NormB(nkgdim,nkgdim,0:ntrunc)
    real(8), intent(in) :: SpVertCorrel(nkgdim,nkgdim,0:ntrunc)
    real(8), intent(in) :: TotVertCorrel(nkgdim,nkgdim)
    real(8), intent(in) :: EnsMean3d(hco_bhi%ni,hco_bhi%nj,nkgdim)
    real(8), intent(in) :: StdDev3dGridPoint(hco_bhi%ni,hco_bhi%nj,nkgdim)
    real(8), intent(in) :: PowerSpectrum(nkgdim,0:ntrunc)
    real(8), intent(in) :: HorizScale(nkgdim)

    integer   :: ier, fstouv, fnom, fstfrm, fclos
    integer   :: iunstats

    write(*,*)
    write(*,*) 'Writing Diagnostics'

    !
    !- Opening Output file
    !
    iunstats = 0
    ier    = fnom(iunstats,'./bgcov_diag.fst','RND',0)
    ier    = fstouv(iunstats,'RND')

    !
    !- Add Tic-Tac and Toc-Toc
    !
    call WriteTicTacToc(iunstats)

    !
    !- Spectral Vertical correlations
    !

    !- Spectral Vertical Correlations
    call WriteSpVertCorrel(SpVertCorrel,iunstats,'ZZ','SPVERTCORREL') ! IN
 
    !- Total Vertical Correlations
    call WriteTotVertCorrel(TotVertCorrel,iunstats,'ZT','TTVERTCORREL') ! IN

    !- Normalized Vertical Correlations
    call WriteSpVertCorrel(NormB,iunstats,'ZN','NRVERTCORREL') ! IN

    !
    !- 3D Grid Point Ensemble Mean
    !
    call Write3d(EnsMean3d,iunstats,'ENSMEAN',cv_bhi) ! IN

    !
    !- 3D Grid Point Standard Deviation
    !
    call Write3d(StdDev3dGridPoint,iunstats,'STDDEV_GRIDP',cv_bhi) ! IN

    !
    !- Power Spectrum
    !
    call WritePowerSpectrum(PowerSpectrum,iunstats,'POWERSPECT',cv_bhi) ! IN

    !
    !- Horizontal Correlation Length scale
    !
    call WriteHorizScale(HorizScale,iunstats,'HORIZSCALE',cv_bhi) ! IN

    !
    !- Closing output file
    !
    ier =  fstfrm(iunstats)
    ier =  fclos (iunstats)

  end subroutine WriteDiagStats

!--------------------------------------------------------------------------
! WriteSpVertCorrel
!--------------------------------------------------------------------------
  subroutine WriteSpVertCorrel(SpVertCorrel,iun,nomvar_in,etiket_in)
    implicit none

    real(8), intent(in) :: SpVertCorrel(nkgdim,nkgdim,0:ntrunc)
    integer, intent(in) :: iun
    character(len=*), intent(in) :: nomvar_in
    character(len=*), intent(in) :: etiket_in

    real(4), allocatable :: work2d(:,:)

    real(4) :: work

    integer   :: ier, fstecr
    integer   :: k, kgdim, totwvnb

    integer :: dateo, npak, ni, nj, nk
    integer :: ip1, ip2, ip3, deet, npas, datyp
    integer :: ig1 ,ig2 ,ig3 ,ig4

    character(len=1)  :: grtyp
    character(len=4)  :: nomvar
    character(len=2)  :: typvar
    character(len=12) :: etiket

    allocate(work2d(nkgdim, nkgdim))

    !- Loop over Total Wavenumbers
    do totwvnb = 0, ntrunc

      npak   = -32
      dateo  = 0
      deet   = 0
      npas   = 0
      ni     = nkgdim
      nj     = nkgdim
      nk     = 1
      ip1    = 0
      ip2    = totwvnb
      ip3    = nens
      typvar = 'XX'
      nomvar = nomvar_in
      etiket = etiket_in
      grtyp  = 'X'
      ig1    = 0
      ig2    = 0
      ig3    = 0
      ig4    = 0
      datyp  = 5

      !- Extract from full Matrix
      work2d(:,:) = real(SpVertCorrel(:,:,totwvnb),4)

      !- Writing 
      ier = fstecr(work2d, work, npak, iun, dateo, deet, npas, ni, nj, &
                   nk, ip1, ip2, ip3, typvar, nomvar, etiket, grtyp,        &
                   ig1, ig2, ig3, ig4, datyp, .true.)

    end do ! Total Wavenumbers

    deallocate(work2d)

  end subroutine WriteSpVertCorrel

!--------------------------------------------------------------------------
! WriteTotVertCorrel
!--------------------------------------------------------------------------
  subroutine WriteTotVertCorrel(TotVertCorrel,iun,nomvar_in,etiket_in)
    implicit none

    real(8), intent(in) :: TotVertCorrel(nkgdim,nkgdim)
    integer, intent(in) :: iun
    character(len=*), intent(in) :: nomvar_in
    character(len=*), intent(in) :: etiket_in

    real(4), allocatable :: workecr(:,:)

    real(4)   :: work

    integer   :: ier, fstecr

    integer :: dateo, npak, ni, nj, nk
    integer :: ip1, ip2, ip3, deet, npas, datyp
    integer :: ig1 ,ig2 ,ig3 ,ig4

    character(len=1)  :: grtyp
    character(len=4)  :: nomvar
    character(len=2)  :: typvar
    character(len=12) :: etiket

    allocate(workecr(nkgdim, nkgdim))

    npak   = -32
    dateo  = 0
    deet   = 0
    npas   = 0
    ni     = nkgdim
    nj     = nkgdim
    nk     = 1
    ip1    = 0
    ip2    = 0
    ip3    = nens
    typvar = 'XX'
    nomvar = nomvar_in
    etiket = etiket_in
    grtyp  = 'X'
    ig1    = 0
    ig2    = 0
    ig3    = 0
    ig4    = 0
    datyp  = 5

    !- Covert to real 4
    workecr(:,:) = real(TotVertCorrel(:,:),4)

    !- Writing 
    ier = fstecr(workecr, work, npak, iun, dateo, deet, npas, ni, nj, &
                 nk, ip1, ip2, ip3, typvar, nomvar, etiket, grtyp,        &
                 ig1, ig2, ig3, ig4, datyp, .true.)

    deallocate(workecr)

  end subroutine WriteTotVertCorrel

!--------------------------------------------------------------------------
! WriteEnsemble
!--------------------------------------------------------------------------
  subroutine WriteEnsemble(ensPerturbations, OutputBaseName, cv_type)
    use MathPhysConstants_mod, only: MPC_KNOTS_PER_M_PER_S_R4
    implicit none

    real(4), intent(in) :: ensPerturbations(hco_bhi%ni,hco_bhi%nj,nkgdim,nens)
    integer, intent(in) :: cv_type
    character(len=*), intent(in) :: OutputBaseName

    real(4), allocatable :: work2d(:,:)

    real(4)   :: factor, work

    integer   :: ier, fstouv, fnom, fstfrm, fclos, fstecr
    integer   :: iunens, ens, var, k, kgdim

    integer :: dateo, npak, ni, nj, nk
    integer :: ip1, ip2, ip3, deet, npas, datyp
    integer :: ig1 ,ig2 ,ig3 ,ig4

    character(len=1)  :: grtyp
    character(len=4)  :: nomvar
    character(len=2)  :: typvar
    character(len=12) :: etiket
    character(len=128):: OutputFileName
    character(len=3)  :: cens

    allocate(work2d(hco_bhi%ni, hco_bhi%nj))

    !- Loop over the Ensemble members
    do ens = 1, nens

      if (ens < 10) then
        write( cens, '(i1)' ) ens
        cens='00'//cens
      else if (ens < 100) then
        write( cens, '(i2)' ) ens
        cens='0'//cens
      else
        write( cens, '(i3)' ) ens
      end if
      OutputFileName= trim(OutputBaseName) // '_' // cens // '.fst'

      write(*,*)
      write(*,*) 'Writing ensemble member: ', ens, trim(OutputFileName)

      iunens = 0
      ier    = fnom(iunens,trim(OutputFileName),'RND',0)
      ier    = fstouv(iunens,'RND')

      !- Add Tic-Tac and Toc-Toc
      call WriteTicTacToc(iunens)

      !- Loop over Control Variables
      do var = 1, nControlVariable

        !- Loop over vertical Levels
        do k = 1, ControlVariable(var)%nlev

          npak   = -32
          dateo  = 0
          deet   = 0
          npas   = 0
          ni     = hco_bhi%ni
          nj     = hco_bhi%nj
          nk     = 1
          ip1    = ControlVariable(var)%ip1(k)
          ip2    = 0
          ip3    = 0
          typvar = 'E'
          nomvar = trim(ControlVariable(var)%nomvar(cv_type))
          etiket = 'DEBUG'
          grtyp  = hco_bhi%grtyp
          ig1    = hco_bhi%ig1
          ig2    = hco_bhi%ig2
          ig3    = hco_bhi%ig3
          ig4    = hco_bhi%ig4
          datyp  = 1

          if ( trim(nomvar) == 'UU' .or. trim(nomvar) == 'VV') then
            factor = MPC_KNOTS_PER_M_PER_S_R4 ! m/s -> knots
          else if ( trim(nomvar) == 'P0' ) then
            factor = 0.01 ! Pa -> hPa
          else
            factor = 1.0
          end if

          !- Extract from EnsPerturbations
          kgdim = ControlVariable(var)%kDimStart + k - 1
          work2d(:,:) = factor * EnsPerturbations(:,:,kgdim,ens)

          !- Writing 
          ier = fstecr(work2d, work, npak, iunens, dateo, deet, npas, ni, nj, &
                       nk, ip1, ip2, ip3, typvar, nomvar, etiket, grtyp,      &
                       ig1, ig2, ig3, ig4, datyp, .true.)

        end do ! Vertical Levels

      end do ! Variables
    
      ier =  fstfrm(iunens)
      ier =  fclos (iunens)

    end do ! Ensemble members

    deallocate(work2d)

  end subroutine WriteEnsemble

!--------------------------------------------------------------------------
! Write3D
!--------------------------------------------------------------------------
  subroutine Write3D(Field3d, iun, Etiket_in, cv_type)
    use MathPhysConstants_mod, only: MPC_KNOTS_PER_M_PER_S_R4
    implicit none

    real(8), intent(in) :: Field3d(hco_bhi%ni,hco_bhi%nj,nkgdim)
    integer, intent(in) :: iun
    integer, intent(in) :: cv_type
    character(len=*), intent(in) :: Etiket_in

    real(4), allocatable :: work2d(:,:)

    real(4)   :: factor, work

    integer   :: ier, fstecr
    integer   :: var, k, kgdim

    integer :: dateo, npak, ni, nj, nk
    integer :: ip1, ip2, ip3, deet, npas, datyp
    integer :: ig1 ,ig2 ,ig3 ,ig4

    character(len=1)  :: grtyp
    character(len=4)  :: nomvar
    character(len=2)  :: typvar
    character(len=12) :: etiket
    character(len=3)  :: cens

    allocate(work2d(hco_bhi%ni, hco_bhi%nj))
    
    !- Loop over Control Variables
    do var = 1, nControlVariable

      !- Loop over vertical Levels
      do k = 1, ControlVariable(var)%nlev
         
        npak   = -32
        dateo  = 0
        deet   = 0
        npas   = 0
        ni     = hco_bhi%ni
        nj     = hco_bhi%nj
        nk     = 1
        ip1    = ControlVariable(var)%ip1(k)
        ip2    = 0
        ip3    = 0
        typvar = 'E'
        nomvar = trim(ControlVariable(var)%nomvar(cv_type))
        etiket = trim(Etiket_in)
        grtyp  = hco_bhi%grtyp
        ig1    = hco_bhi%ig1
        ig2    = hco_bhi%ig2
        ig3    = hco_bhi%ig3
        ig4    = hco_bhi%ig4
        datyp  = 1

        if ( trim(nomvar) == 'UU' .or. trim(nomvar) == 'VV') then
          factor = MPC_KNOTS_PER_M_PER_S_R4 ! m/s -> knots
        else if ( trim(nomvar) == 'P0' ) then
          factor = 0.01 ! Pa -> hPa
        else
          factor = 1.0
        end if

        !- Extract from EnsPerturbations
        kgdim = ControlVariable(var)%kDimStart + k - 1
        work2d(:,:) = factor * real(Field3d(:,:,kgdim),4)

        !- Writing 
        ier = fstecr(work2d, work, npak, iun, dateo, deet, npas, ni, nj, &
                     nk, ip1, ip2, ip3, typvar, nomvar, etiket, grtyp,      &
                     ig1, ig2, ig3, ig4, datyp, .true.)

      end do ! Vertical Levels

    end do ! Variables

    deallocate(work2d)

  end subroutine Write3D

!--------------------------------------------------------------------------
! WritePowerSpectrum
!--------------------------------------------------------------------------
  subroutine WritePowerSpectrum(PowerSpectrum,iun,etiket_in,cv_type)
    implicit none

    real(8), intent(in) :: PowerSpectrum(nkgdim,0:ntrunc)
    integer, intent(in) :: iun
    integer, intent(in) :: cv_type
    character(len=*), intent(in) :: Etiket_in

    real(4), allocatable :: workecr(:,:)

    real(4)   :: factor, work

    integer   :: ier, fstecr
    integer   :: var, k, kgdim

    integer :: dateo, npak, ni, nj, nk
    integer :: ip1, ip2, ip3, deet, npas, datyp
    integer :: ig1 ,ig2 ,ig3 ,ig4

    character(len=1)  :: grtyp
    character(len=4)  :: nomvar
    character(len=2)  :: typvar
    character(len=12) :: etiket
    character(len=3)  :: cens

    allocate(workecr(ntrunc+1, 1))
    
    !- Loop over Control Variables
    do var = 1, nControlVariable

      !- Loop over vertical Levels
      do k = 1, ControlVariable(var)%nlev

        npak   = -32
        dateo  = 0
        deet   = 0
        npas   = 0
        ni     = ntrunc + 1
        nj     = 1
        nk     = 1
        ip1    = ControlVariable(var)%ip1(k)
        ip2    = 0
        ip3    = 0
        typvar = 'E'
        nomvar = trim(ControlVariable(var)%nomvar(cv_type))
        etiket = trim(Etiket_in)
        grtyp  = 'X'
        ig1    = 0
        ig2    = 0
        ig3    = 0 
        ig4    = 0
        datyp  = 1

        !- Extract
        kgdim = ControlVariable(var)%kDimStart + k - 1
        workecr(:,1) = real(PowerSpectrum(kgdim,:),4)

        !- Writing 
        ier = fstecr(workecr, work, npak, iun, dateo, deet, npas, ni, nj, &
                     nk, ip1, ip2, ip3, typvar, nomvar, etiket, grtyp,   &
                     ig1, ig2, ig3, ig4, datyp, .true.)

      end do ! Vertical Levels

    end do ! Variables

    deallocate(workecr)

  end subroutine WritePowerSpectrum

!--------------------------------------------------------------------------
! WriteHorizScale
!--------------------------------------------------------------------------
  subroutine WriteHorizScale(HorizScale,iun,etiket_in,cv_type)
    implicit none

    real(8), intent(in) :: HorizScale(nkgdim)
    integer, intent(in) :: iun
    integer, intent(in) :: cv_type
    character(len=*), intent(in) :: Etiket_in

    real(4), allocatable :: workecr(:,:,:)

    real(4)   :: factor, work

    integer   :: ier, fstecr
    integer   :: var, k, kgdim

    integer :: dateo, npak, ni, nj, nk
    integer :: ip1, ip2, ip3, deet, npas, datyp
    integer :: ig1 ,ig2 ,ig3 ,ig4

    character(len=1)  :: grtyp
    character(len=4)  :: nomvar
    character(len=2)  :: typvar
    character(len=12) :: etiket
    character(len=3)  :: cens
    
    !- Loop over Control Variables
    do var = 1, nControlVariable

      allocate(workecr(1,1,ControlVariable(var)%nlev))

      npak   = -32
      dateo  = 0
      deet   = 0
      npas   = 0
      ni     = 1
      nj     = 1
      nk     = ControlVariable(var)%nlev
      ip1    = 0
      ip2    = 0
      ip3    = 0
      typvar = 'E'
      nomvar = trim(ControlVariable(var)%nomvar(cv_type))
      etiket = trim(Etiket_in)
      grtyp  = 'X'
      ig1    = 0
      ig2    = 0
      ig3    = 0
      ig4    = 0
      datyp  = 1

      !- Extract
      workecr(1,1,:) = real(HorizScale(ControlVariable(var)%kDimStart:ControlVariable(var)%kDimEnd),4)

      !- Writing 
      ier = fstecr(workecr, work, npak, iun, dateo, deet, npas, ni, nj, &
                   nk, ip1, ip2, ip3, typvar, nomvar, etiket, grtyp,   &
                   ig1, ig2, ig3, ig4, datyp, .true.)

      deallocate(workecr)
      
    end do ! Variables

  end subroutine WriteHorizScale

!--------------------------------------------------------------------------
! WriteTicTacToc
!--------------------------------------------------------------------------
  subroutine WriteTicTacToc(iun)
    use vGrid_Descriptors , only: vgrid_descriptor, vgd_write, VGD_OK
    use MathPhysConstants_mod, only : MPC_DEGREES_PER_RADIAN_R8
    implicit none

    integer, intent(in) :: iun

    integer :: ier
    integer :: dateo,npak, status
    integer :: ip1,ip2,ip3,deet,npas,datyp,ig1,ig2,ig3,ig4
    integer :: ig1_tictac,ig2_tictac,ig3_tictac,ig4_tictac

    character(len=1)  :: grtyp
    character(len=2)  :: typvar
    character(len=12) :: etiket

    !
    !- 1.  Writing Tic-Tac
    !
    npak     = -32
    deet     =  0
    ip1      =  hco_bhi%ig1
    ip2      =  hco_bhi%ig2
    ip3      =  0
    npas     =  0
    datyp    =  1
    grtyp    = 'E'
    typvar   = 'X'
    etiket   = 'EXTENDEDGRID'
    dateo =  0

    call cxgaig ( grtyp,                                          & ! IN
                  ig1_tictac, ig2_tictac, ig3_tictac, ig4_tictac, & ! OUT
                  real(hco_bhi%xlat1), real(hco_bhi%xlon1),   & ! IN
                  real(hco_bhi%xlat2), real(hco_bhi%xlon2)  )   ! IN

    ig1      =  ig1_tictac
    ig2      =  ig2_tictac
    ig3      =  ig3_tictac
    ig4      =  ig4_tictac

    ier = utl_fstecr(hco_bhi%lon*MPC_DEGREES_PER_RADIAN_R8, npak, &
                     iun, dateo, deet, npas, hco_bhi%ni, 1, 1, ip1,    &
                     ip2, ip3, typvar, '>>', etiket, grtyp, ig1,          &
                     ig2, ig3, ig4, datyp, .true.)

    ier = utl_fstecr(hco_bhi%lat*MPC_DEGREES_PER_RADIAN_R8, npak, &
                     iun, dateo, deet, npas, 1, hco_bhi%nj, 1, ip1,    &
                     ip2, ip3, typvar, '^^', etiket, grtyp, ig1,          &
                     ig2, ig3, ig4, datyp, .true.)

    !
    !- Writing Toc-Toc
    !
    if ( spatialDimensions /= '2D' ) then
       status = vgd_write(vco_bhi%vgrid,iun,'fst')

       if ( status /= VGD_OK ) then
          write(*,*)
          write(*,*) 'WriteTicTacToc: ERROR with vgd_write '
          call utl_abort('WriteTicTacToc')
       end if
    end if

  end subroutine WriteTicTacToc

!--------------------------------------------------------------------------
! WriteControlVarInfo
!--------------------------------------------------------------------------
  subroutine WriteControlVarInfo(iun)
    implicit none

    integer, intent(in) :: iun

    integer :: ier, fstecr, fstecr_s

    real(8) :: work

    integer :: npak, var, dateo, ni, nj
    integer :: ip1,ip2,ip3,deet,npas,datyp,ig1,ig2,ig3,ig4
    integer :: ig1_tictac,ig2_tictac,ig3_tictac,ig4_tictac

    character(len=1)  :: grtyp
    character(len=2)  :: typvar
    character(len=4)  :: nomvar
    character(len=12) :: etiket

    character(len=4)  :: ControlModelVarnameList(nControlVariable)
    character(len=4)  :: ControlBhiVarnameList  (nControlVariable)
    character(len=2)  :: ControlVarGridTypeList (nControlVariable)
    integer           :: ControlVarNlevList     (nControlVariable)

    !
    !- 1. Gathering the info
    !

    !- Loop over Control Variables
    do var = 1, nControlVariable
       ControlModelVarnameList(var) = trim(ControlVariable(var)%nomvar(cv_model))
       ControlBhiVarnameList(var)   = trim(ControlVariable(var)%nomvar(cv_bhi))
       ControlVarNlevList(var)      = ControlVariable(var)%nlev
       ControlVarGridTypeList(var)  = ControlVariable(var)%GridType
    end do

    write(*,*)
    write(*,*) 'ControlModelVarnameList = ',ControlModelVarnameList(:)
    write(*,*) 'ControlBhiVarnameList   = ',ControlBhiVarnameList(:)
    write(*,*) 'ControlVarNlevList      = ',ControlVarNlevList(:)
    write(*,*) 'ControlVarGridTypeList  = ',ControlVarGridTypeList(:)

    !
    !- 2.  Writing the list of control variables and number of vertical levels
    !
    npak     = -32
    dateo    =  0
    deet     =  0
    ip1      =  0
    ip2      =  0
    ip3      =  0
    npas     =  0
    grtyp    = 'X'
    typvar   = 'X'
    ig1      =  0
    ig2      =  0
    ig3      =  0
    ig4      =  0

    nomvar   = 'CVN'
    ni       =  4  ! 4 Characters
    nj       =  nControlVariable
    datyp    =  7 ! Character

    ier = fstecr_s(ControlModelVarnameList, work, npak, &
                 iun, dateo, deet, npas, ni, nj, 1, ip1,    &
                 ip2, ip3, typvar, nomvar, 'MODEL', grtyp, ig1, &
                 ig2, ig3, ig4, datyp, .true.)

    ier = fstecr_s(ControlBhiVarnameList, work, npak, &
                 iun, dateo, deet, npas, ni, nj, 1, ip1,    &
                 ip2, ip3, typvar, nomvar, 'B_HI', grtyp, ig1, &
                 ig2, ig3, ig4, datyp, .true.)

    nomvar   = 'CVL'
    ni       =  2  ! 2 Characters
    ier = fstecr_s(ControlVarGridTypeList, work, npak, &
                 iun, dateo, deet, npas, ni, nj, 1, ip1,    &
                 ip2, ip3, typvar, nomvar, 'LEVTYPE', grtyp, ig1, &
                 ig2, ig3, ig4, datyp, .true.)
 
    datyp    =  2 ! Integer
    ni       =  nControlVariable
    nj       =  1 
    ier = fstecr(ControlVarNlevList, work, npak, &
                 iun, dateo, deet, npas, ni, nj, 1, ip1,    &
                 ip2, ip3, typvar, nomvar, 'NLEV', grtyp, ig1, &
                 ig2, ig3, ig4, datyp, .true.)

  end subroutine WriteControlVarInfo

  !--------------------------------------------------------------------------
  ! WRITEPRESSUREPROFILES
  !--------------------------------------------------------------------------
  subroutine writePressureProfiles
    implicit none

    character(len=128) :: outfilename

    integer :: jk

    outfilename = "./pressureProfile_M.txt"
    open (unit=99,file=outfilename,action="write",status="new")
    do jk = 1, vco_bhi%nlev_M
       write(99,'(I3,2X,F6.1)') jk, pressureProfile_M(jk)/100.d0
    end do
    close(unit=99)
       
    outfilename = "./pressureProfile_T.txt"
    open (unit=99,file=outfilename,action="write",status="new")
    do jk = 1, vco_bhi%nlev_T
       write(99,'(I3,2X,F6.1)') jk, pressureProfile_T(jk)/100.d0
    end do
    close(unit=99)

    write(6,*) 'finished writing pressure profiles...'
    call flush(6)

  end subroutine writePressureProfiles

end module calcstatslam_mod
