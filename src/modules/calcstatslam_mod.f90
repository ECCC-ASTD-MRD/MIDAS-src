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
!! MODULE calcStatsLam (prefix="csl" category='1. High-level functionality')
!!
!! *Purpose*: Compute homogeneous and isotropic background error covariances 
!!            from forecast error estimate in model variable space (limited-area version).
!!
!--------------------------------------------------------------------------
module calcstatslam_mod
  use gridStateVector_mod
  use lamSpectralTransform_mod
  use analysisGrid_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use localizationFunction_mod
  use utilities_mod
  use menetrierDiag_mod
  use ensemblestatevector_mod
  use variableTransforms_mod
  use varNameList_mod
  use gridBinning_mod
  use timeCoord_mod
  implicit none
  save
  private

  ! Public Subroutines
  public :: csl_setup, csl_computeBhi
  public :: csl_toolbox

  type(struct_hco), pointer :: hco_ens => null() ! Ensemble horizontal grid parameters
  type(struct_hco), pointer :: hco_bhi => null() ! B matrix horizontal grid parameters
  type(struct_vco), pointer :: vco_bhi => null() ! B matrix vertical grid parameters
  type(struct_lst)          :: lst_bhi ! Spectral transform Parameters

  integer,external   :: get_max_rss

  integer, parameter :: cv_model = 1
  integer, parameter :: cv_bhi   = 2
  integer, parameter :: nMaxControlVar = 10
  
  type  :: struct_cv
    character(len=4)     :: NomVar(2)
    integer              :: varLevIndexStart
    integer              :: varLevIndexEnd
    integer              :: nlev
    character(len=2)     :: GridType
    integer, allocatable :: ip1(:)
  end type struct_cv

  type  :: struct_bhi
    type(struct_cv)  :: controlVariable(nMaxControlVar)
    character(len=4) :: momentumControlVar(2)
    integer          :: nControlVariable
    integer          :: nVarLev
  end type struct_bhi

  type(struct_bhi)  :: bhi

  integer :: nens
  integer :: ntrunc
  integer :: ip2_ens

  logical :: initialized = .false.
  logical :: NormByStdDev, SetTGtoZero, writeEnsPert
  logical :: vertLoc, horizLoc, correlationLQ
  logical :: ensContainsFullField

  character(len=12) :: WindTransform
  character(len=12) :: SpectralWeights
  character(len=2)  :: ctrlVarHumidity

  type(struct_ens)  :: ensPerts

  real(8), pointer  :: pressureProfile_M(:), pressureProfile_T(:)

  real(8) :: vLocalize_wind, vlocalize_mass, vlocalize_humidity ! vertical length scale (in units of ln(Pressure))
  real(8) :: hlocalize_wind, hlocalize_mass, hlocalize_humidity ! vertical length scale (in km)

contains
  
  !--------------------------------------------------------------------------
  ! csl_setup
  !--------------------------------------------------------------------------
  subroutine csl_setup(nens_in, ensFileName, hco_ens_in, vco_ens_in, ip2_in)
    use vGrid_Descriptors , only: vgrid_descriptor, vgd_levels, VGD_OK  
    implicit none

    integer,                   intent(in)   :: nens_in
    type(struct_vco), pointer, intent(in)   :: vco_ens_in
    type(struct_hco), pointer, intent(in)   :: hco_ens_in
    integer,                   intent(in)   :: ip2_in
    character(len=*), intent(in)   :: ensFileName

    integer :: nulnam, ier, status
    integer :: fclos, fnom, fstouv, fstfrm
    integer :: grd_ext_x, grd_ext_y
    integer :: varIndex, k

    integer :: numStep
    integer, allocatable :: dateStampList(:)

    real(8) :: SurfacePressure

    character(len=256)  :: enspathname
    logical :: makeBiPeriodic

    character(len=4), pointer :: controlVarNames(:)

    NAMELIST /NAMCALCSTATS_LAM/ntrunc,grd_ext_x,grd_ext_y,WindTransform,NormByStdDev, &
                               SetTGtoZero,writeEnsPert,SpectralWeights,              &
                               vLocalize_wind,vlocalize_mass,vlocalize_humidity,      &
                               hLocalize_wind,hlocalize_mass,hlocalize_humidity,      &
                               correlationLQ

    write(*,*)
    write(*,*) 'csl_setup: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- 1. Initialized the info on the ensemble
    !
    nens=nens_in

    hco_ens => hco_ens_in
    vco_bhi => vco_ens_in

    if ( vco_bhi%nlev_T == 1 .and. vco_bhi%nlev_T == 1 ) then
      write(*,*)
      write(*,*) 'Spatial dimensions = 2D'
    else
      write(*,*)
      write(*,*) 'Spatial dimensions = 3D'
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
    call createLamTemplateGrids('./analysisgrid',grd_ext_x,grd_ext_y) ! IN

    !- 3.2 Setup the Extended B_HI grid
    call hco_setupFromFile( hco_bhi,'./analysisgrid', 'ANALYSIS', 'BHI' ) ! IN

    !- 3.3 Setup the LAM analysis grid metrics
    call agd_setupFromHCO( hco_bhi, hco_ens ) ! IN

    !- 3.4 Setup the LAM spectral transform
    call lst_setup( lst_bhi,     & ! OUT
         hco_bhi%ni, hco_bhi%nj, & ! IN
         hco_bhi%dlon, ntrunc,   & ! IN
         'NoMpi' )                 ! IN

    !
    !- 4. Read the ensemble
    !
    numStep = 1
    allocate(dateStampList(numStep))
    dateStampList(:)  = -1
    call ens_allocate(ensPerts, nEns, numStep, hco_bhi, vco_bhi, dateStampList)

    ensContainsFullField = .false.
    ctrlVarHumidity = 'LQ'
    ensPathName    = './ensemble'
    makeBiPeriodic = .true.
    call ens_readEnsemble(ensPerts, ensPathName, makeBiPeriodic, &
                          containsFullField_opt=ensContainsFullField)

    if ( ctrlVarHumidity == 'LQ' .and. ens_varExist(ensPerts,'HU') .and. &
         ensContainsFullField ) then
      call vtr_transform(ensPerts,'HUtoLQ')
    end if

    !
    !- 5.  Setup the control variables (model space and B_hi space)
    !
    nullify(controlVarNames)
    call ens_varNamesList(ensPerts,controlVarNames)

    bhi%nControlVariable = size(controlVarNames) !count(mask) 
    write(*,*)
    write(*,*) 'Number of Control Variables = ', bhi%nControlVariable

    !- 5.1 Set ControlVariable structure
    bhi%momentumControlVar(:) = 'NULL' 
    do varIndex = 1, bhi%nControlVariable

      !- Set variable name
      bhi%controlVariable(varIndex)%nomvar(cv_model)  = trim(controlVarNames(varIndex))
      if ( bhi%controlVariable(varIndex)%nomvar(cv_model) == 'UU' ) then
        if ( trim(WindTransform) == 'PsiChi') then
          write(*,*)
          write(*,*) '--- Momentum Control Variables = Psi-Chi ---'
          bhi%momentumControlVar(1) = 'PP'
        else if ( trim(WindTransform) == 'VortDiv') then
          write(*,*)
          write(*,*) '--- Momentum Control Variables = Vort-Div ---'
          bhi%momentumControlVar(1) = 'QR'
          call utl_abort('Momentum Control Variables = Vort-Div not yest availble')
        else if ( trim(WindTransform) == 'UV') then
          write(*,*)
          write(*,*) '--- Momentum Control Variables = U-V ---'
          bhi%momentumControlVar(1) = 'UU'
        else
          write(*,*)
          write(*,*) 'Wind Transform not available ', trim((WindTransform))
          call utl_abort('csl_setup')
        end if
        bhi%controlVariable(varIndex)%nomvar(cv_bhi) = trim(bhi%momentumControlVar(1))
      else if ( bhi%controlVariable(varIndex)%nomvar(cv_model) == 'VV' ) then
        if ( trim(WindTransform) == 'PsiChi') then
          bhi%momentumControlVar(2) = 'CC'
        else if ( trim(WindTransform) == 'VortDiv') then
          bhi%momentumControlVar(2) = 'DD'
        else if ( trim(WindTransform) == 'UV') then
          bhi%momentumControlVar(2) = 'VV'
        else
          write(*,*)
          write(*,*) 'Wind Transform not available ', trim((WindTransform))
          call utl_abort('csl_setup')
        end if
        bhi%controlVariable(varIndex)%nomvar(cv_bhi) = bhi%momentumControlVar(2)
      else
        bhi%controlVariable(varIndex)%nomvar(cv_bhi) = bhi%controlVariable(varIndex)%nomvar(cv_model)
      end if

      !- 5.2 Set Level info
      bhi%controlVariable(varIndex)%GridType =  vnl_varLevelFromVarname(bhi%controlVariable(varIndex)%nomvar(cv_model))
      bhi%controlVariable(varIndex)%nlev     =  ens_getNumLev(ensPerts,bhi%controlVariable(varIndex)%GridType)

      write(*,*)
      write(*,*) 'Control Variable Name ', bhi%controlVariable(varIndex)%nomvar(cv_model)
      write(*,*) '   Number of Levels = ', bhi%controlVariable(varIndex)%nlev
      write(*,*) '   Type   of Levels = ', bhi%controlVariable(varIndex)%GridType

      allocate( bhi%controlVariable(varIndex)%ip1 (bhi%controlVariable(varIndex)%nlev) )

      if (bhi%controlVariable(varIndex)%GridType == 'TH') then
        bhi%controlVariable(varIndex)%ip1(:) = vco_bhi%ip1_T(:)
      else if (bhi%controlVariable(varIndex)%GridType == 'MM') then
        bhi%controlVariable(varIndex)%ip1(:) = vco_bhi%ip1_M(:)
      else if (bhi%controlVariable(varIndex)%GridType == 'SF') then
        bhi%controlVariable(varIndex)%ip1(:) = 0
      end if

      bhi%controlVariable(varIndex)%varLevIndexStart =  &
           ens_getOffsetFromVarName(ensPerts,bhi%controlVariable(varIndex)%nomvar(cv_model)) + 1 
      bhi%controlVariable(varIndex)%varLevIndexEnd   =  &
           ens_getOffsetFromVarName(ensPerts,bhi%controlVariable(varIndex)%nomvar(cv_model)) + &
           bhi%controlVariable(varIndex)%nlev

    end do

    bhi%nVarLev = ens_getNumK(ensPerts)

    !
    !- 6.  Transform u-wind and v-wind to control variables 
    !
    if (writeEnsPert) then
      call ens_writeEnsemble(ensPerts, './', 'MODELVAR_', ctrlVarHumidity, 'MODELVAR', 'E', &
                             containsFullField_opt = ensContainsFullField)
    end if

    if      ( bhi%momentumControlVar(1) == 'PP' .and. bhi%momentumControlVar(2) == 'CC' ) then
      call vtr_transform(ensPerts,'UVtoPsiChi')
      call ens_modifyVarName(ensPerts, 'UU', 'PP')
      call ens_modifyVarName(ensPerts, 'VV', 'CC')
    else if ( bhi%momentumControlVar(1) == 'QR' .and. bhi%momentumControlVar(2) == 'DD' ) then
      call vtr_transform(ensPerts,'UVtoVortDiv')
      call ens_modifyVarName(ensPerts, 'UU', 'QR')
      call ens_modifyVarName(ensPerts, 'VV', 'DD')
    end if

    if (writeEnsPert) then
      call ens_writeEnsemble(ensPerts, './', 'CTRLVAR', ctrlVarHumidity, 'CTRLVAR', 'E', &
           containsFullField_opt = ensContainsFullField)
    end if

    !
    !- 7.  Compute and remove the ensemble mean; compute the stdDev
    !
    call ens_computeMean(ensPerts)
    call ens_removeMean (ensPerts)
    call ens_computeStdDev(ensPerts)

    !
    !- 8.  Setup localization 
    !

    !- 8.1 Setup horizontal localization
    if ( hLocalize_wind < 0.d0 .and. hLocalize_mass < 0.d0 .and. hLocalize_humidity < 0.d0 ) then
      write(*,*) 
      write(*,*) 'csl_setup: NO horizontal correlation localization will be performed'
      horizLoc=.false.
    else
      write(*,*) 
      write(*,*) 'csl_setup: horizontal correlation localization WILL BE performed'
      horizLoc=.true.
    end if

    !- 8.2 Setup vertical localization
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
    !- 9.  Setup pressure profile for vertical localization
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
    !- 10.  Ending
    !
    initialized = .true.

    write(*,*)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'csl_setup: Done!'

  end subroutine csl_setup

  !--------------------------------------------------------------------------
  ! csl_computeBhi
  !--------------------------------------------------------------------------
  subroutine csl_computeBhi
    implicit none

    integer :: ier

    real(8),allocatable :: SpVertCorrel(:,:,:)
    real(8),allocatable :: TotVertCorrel(:,:)
    real(8),allocatable :: NormB(:,:,:)
    real(8),allocatable :: NormBsqrt(:,:,:)
    real(8),allocatable :: PowerSpectrum(:,:)
    real(8),allocatable :: HorizScale(:)

    character(len=4), pointer :: varNamesList(:)

    type(struct_gbi) :: gbi_horizontalMean
    type(struct_gsv) :: statevector_stdDev
    type(struct_gsv) :: statevector_mean, statevector_stdDevGridPoint

    write(*,*)
    write(*,*) 'csl_computeBhi: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    allocate(SpVertCorrel(bhi%nVarLev,bhi%nVarLev,0:ntrunc))
    allocate(TotVertCorrel(bhi%nVarLev,bhi%nVarLev))
    allocate(PowerSpectrum(bhi%nVarLev,0:ntrunc))
    allocate(NormB(bhi%nVarLev,bhi%nVarLev,0:ntrunc))
    allocate(NormBsqrt(bhi%nVarLev,bhi%nVarLev,0:ntrunc))
    allocate(HorizScale(bhi%nVarLev))

    !
    !- 1.  Calculate the gridded binned Std. Dev. to be used in the analysis step
    !
    nullify(varNamesList)
    call ens_varNamesList(ensPerts,varNamesList) 
    call gsv_allocate(statevector_stdDev, ens_getNumStep(ensPerts),                            &
                      ens_getHco(ensPerts), ens_getVco(ensPerts), varNames_opt=varNamesList, &
                      datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true.,      &
                      dataKind_opt=8 )

    call gbi_setup(gbi_horizontalMean, 'HorizontalMean', statevector_stdDev)

    call gbi_stdDev(gbi_horizontalMean, ensPerts, & ! IN
                    statevector_stdDev)             ! OUT

    !
    !- 2.  Normalization of the ensemble perturbations
    !
    if ( NormByStdDev ) then
      call ens_normalize(ensPerts)
    end if

    !
    !- 3.  Covariance statistics in Spectral Space
    !

    !- 3.1 Vertical correlations and Power Spectra
    call calcSpectralStats(ensPerts,                      & ! IN
                           SpVertCorrel, PowerSpectrum,   & ! OUT
                           NormB)                           ! OUT

    !- 3.2 Calculate the horiontal correlation lenght scales
    call calcHorizScale(HorizScale, & ! OUT
                        NormB)        ! IN

    !- 3.3 Calculate the total vertical correlation matrix
    call calcTotVertCorrel(TotVertCorrel,             & ! OUT
                           SpVertCorrel, PowerSpectrum) ! IN

    !- 3.4 Set cross-correlations
    call setSpVertCorrel(NormB) ! INOUT

    !- 3.5 Calculate the square-root of the correlation-based B matrix
    call calcBsqrt(NormBsqrt, & ! OUT
                   NormB   )    ! IN

    !
    !- 4.  Writing statistics to files
    !

    !- 4.1 Statistics needed by VAR in analysis mode
    call writeVarStats(NormBsqrt, statevector_stdDev) ! IN

    !- 4.2 Diagnostics fields
    call ens_copyEnsMean(ensPerts, statevector_mean)
    call ens_copyEnsStdDev(ensPerts, statevector_stdDevGridPoint)

    call writeDiagStats(NormB, SpVertCorrel, TotVertCorrel, statevector_mean, & ! IN 
                        statevector_stdDevGridPoint, PowerSpectrum, HorizScale) ! IN

    !
    !- 5.  Ending
    !
    call ens_deallocate(ensPerts)
    call gsv_deallocate(statevector_mean)
    call gsv_deallocate(statevector_stdDev)
    call gsv_deallocate(statevector_stdDevGridPoint)

    deallocate(SpVertCorrel)
    deallocate(TotVertCorrel)
    deallocate(PowerSpectrum)
    deallocate(NormB)
    deallocate(NormBsqrt)
    deallocate(HorizScale)

    write(*,*)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'csl_computeBhi: Done!'

  end subroutine csl_computeBhi

  !--------------------------------------------------------------------------
  ! csl_toolbox
  !--------------------------------------------------------------------------
  subroutine csl_toolbox
    implicit none

!!$    integer :: nulnam, ier, fstouv, fnom, fstfrm, fclos
!!$    integer :: iunstats
!!$
!!$    integer :: NumBins2d
!!$
!!$    real(4),allocatable :: ensPerturbations(:,:,:,:)
!!$
!!$    real(8),allocatable :: ensMean3d(:,:,:)
!!$    real(8),allocatable :: StdDev3dGridPoint(:,:,:)
!!$    real(8),allocatable :: VertCorrel(:,:,:)
!!$
!!$    real(8),allocatable :: SpVertCorrel(:,:,:)
!!$    real(8),allocatable :: NormB(:,:,:)
!!$    real(8),allocatable :: PowerSpectrum(:,:)
!!$
!!$    integer,allocatable :: Bin2d(:,:)
!!$
!!$    character(len=4), allocatable :: nomvar3d(:,:),nomvar2d(:,:)
!!$
!!$    integer :: varLevOffset(6), nvar3d, nvar2d, var3d, var2d, var
!!$
!!$    character(len=60) :: tool
!!$
!!$    NAMELIST /NAMTOOLBOX/tool
!!$
!!$    write(*,*)
!!$    write(*,*) 'csl_toolbox: Starting...'
!!$    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
!!$
!!$    !
!!$    !- Tool selection
!!$    !
!!$    nulnam = 0
!!$    ier = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
!!$    read(nulnam,nml=NAMTOOLBOX)
!!$    write(*,nml=NAMTOOLBOX)
!!$    ier = fclos(nulnam)
!!$
!!$    select case(trim(tool))
!!$    case ('VERTCORREL_GRIDPOINT')
!!$       write(*,*)
!!$       write(*,*) 'Computing vertical correlation in grid point space'
!!$    case ('HORIZCORREL_FUNCTION')
!!$       write(*,*)
!!$       write(*,*) 'Computing horizontal correlation functions'
!!$    case ('STDDEV')
!!$       write(*,*)
!!$       write(*,*) 'Computing Standard-Deviations'
!!$    case ('HVCORREL_LOCAL')
!!$       write(6,*)
!!$       write(6,*) 'Computing Local Correlation'
!!$    case ('LOCALIZATIONRADII')
!!$       write(6,*)
!!$       write(6,*) 'Estimating the optimal covariance localization radii'
!!$
!!$       nVar3d = 0
!!$       nVar2d = 0
!!$       do var = 1, bhi%nControlVariable
!!$         if ( bhi%controlVariable(var)%nlev == 1 ) then
!!$           nVar2d = nVar2d + 1
!!$         else
!!$           nVar3d = nVar3d + 1
!!$         end if
!!$       end do
!!$       allocate(nomvar3d(nVar3d,3))
!!$       allocate(nomvar2d(nVar2d,3))
!!$
!!$       var3d = 0
!!$       var2d = 0
!!$       do var = 1, bhi%nControlVariable
!!$         if ( bhi%controlVariable(var)%nlev == 1 ) then
!!$           var2d = var2d + 1
!!$           nomvar2d(var2d,1) = bhi%controlVariable(var)%nomvar(cv_model)
!!$           nomvar2d(var2d,2) = bhi%controlVariable(var)%nomvar(cv_bhi)
!!$         else
!!$           var3d = var3d + 1
!!$           nomvar3d(var3d,1) = bhi%controlVariable(var)%nomvar(cv_model)
!!$           nomvar3d(var3d,2) = bhi%controlVariable(var)%nomvar(cv_bhi)
!!$         end if
!!$         varLevOffset(var) = bhi%controlVariable(var)%varLevIndexStart - 1
!!$       end do
!!$       call bmd_setup( hco_ens, nens, vco_bhi%nlev_M, vco_bhi%nlev_T,    & ! IN
!!$                       bhi%nVarLev, pressureProfile_M, pressureProfile_T,     & ! IN
!!$                       nvar3d, nvar2d, varLevOffset, nomvar3d, nomvar2d, & ! IN
!!$                       1)
!!$    case default
!!$       write(*,*)
!!$       write(*,*) 'Unknown TOOL in csl_toolbox : ', trim(tool)
!!$       call utl_abort('calbmatrix_lam')
!!$    end select
!!$
!!$    !
!!$    !- Horizontal and vertical correlation diagnostics
!!$    !
!!$    allocate(ensPerturbations(hco_bhi%ni,hco_bhi%nj,bhi%nVarLev,nens),stat=ier)
!!$    if(ier /= 0) then
!!$      write(*,*) 'Problem allocating ensPerturbations memory!',ier
!!$      write(*,*)  hco_bhi%ni,hco_bhi%nj,bhi%nVarLev,nens
!!$      call utl_abort('csl_computeStats') 
!!$    endif
!!$    allocate(ensMean3d(hco_bhi%ni,hco_bhi%nj,bhi%nVarLev),stat=ier)
!!$    if (ier /= 0) then
!!$      write(*,*) 'Problem allocating ensMean3d memory!',ier
!!$      call utl_abort('csl_computeStats')
!!$    end if
!!$    allocate(StdDev3dGridPoint(hco_bhi%ni,hco_bhi%nj,bhi%nVarLev),stat=ier)
!!$    if (ier /= 0) then
!!$      write(*,*) 'Problem allocating StdDev3dGridPoint memory!',ier
!!$      call utl_abort('csl_computeStats')
!!$    end if
!!$    allocate(Bin2d(hco_bhi%ni,hco_bhi%nj),stat=ier)
!!$    if (ier /= 0) then
!!$      write(*,*) 'Problem allocating Bin2d memory!',ier
!!$      call utl_abort('csl_computeStats')
!!$    end if
!!$
!!$
!!$    if ( trim(tool) == 'STDDEV' ) then
!!$       iunstats = 0
!!$       ier    = fnom(iunstats,'./stddev.fst','RND',0)
!!$       ier    = fstouv(iunstats,'RND')
!!$       call Write3d(StdDev3dGridPoint,iunstats,'STDDEV_GRIDP',cv_bhi) ! IN
!!$       call WriteTicTacToc(iunstats) ! IN
!!$       ier =  fstfrm(iunstats)
!!$       ier =  fclos (iunstats)
!!$       return
!!$    else if ( trim(tool) == 'VERTCORREL_GRIDPOINT' ) then
!!$       !
!!$       !- 6.  VERTCORREL_GRIDPOINT
!!$       !
!!$       
!!$       !- 6.0 Normalization
!!$       call Normalize3d(ensPerturbations, & ! INOUT
!!$                        StdDev3dGridPoint)  ! IN
!!$
!!$       !- 6.1  Remove the (2D) domain mean
!!$       call removeDomainMean(ensPerturbations) ! INOUT
!!$
!!$       !- 6.2  Calculate the vertical correlations
!!$       iunstats = 0
!!$       ier    = fnom(iunstats,'./vertCorrel.fst','RND',0)
!!$       ier    = fstouv(iunstats,'RND')
!!$       
!!$       !- 6.2.1 Stats from the whole horizontal domain
!!$       call CreateBins(Bin2d, NumBins2d, & ! OUT
!!$                      'FullDomain')        ! IN
!!$
!!$       allocate(VertCorrel(bhi%nVarLev,bhi%nVarLev,NumBins2d))
!!$
!!$       call calcVertCorrel(VertCorrel,       &  ! OUT
!!$                           ensPerturbations, &  ! IN
!!$                           Bin2d, NumBins2d)    ! IN
!!$
!!$       call WriteTotVertCorrel(VertCorrel(:,:,1),iunstats,'ZT','GPVCOR_FULL') ! IN
!!$
!!$       deallocate(VertCorrel)
!!$
!!$       !- 6.2.2 Stats from the core and the extension zone
!!$       call CreateBins(Bin2d, NumBins2d, & ! OUT
!!$                       'CoreExt')          ! IN
!!$
!!$       allocate(VertCorrel(bhi%nVarLev,bhi%nVarLev,NumBins2d))
!!$
!!$       call calcVertCorrel(VertCorrel,       &  ! OUT
!!$                           ensPerturbations, &  ! IN
!!$                           Bin2d, NumBins2d)    ! IN
!!$
!!$       call WriteTotVertCorrel(VertCorrel(:,:,1),iunstats,'ZT','GPVCOR_CORE') ! IN
!!$       call WriteTotVertCorrel(VertCorrel(:,:,2),iunstats,'ZT','GPVCOR_EXTN') ! IN
!!$
!!$       deallocate(VertCorrel)
!!$
!!$       ier =  fstfrm(iunstats)
!!$       ier =  fclos (iunstats)
!!$
!!$     else if (trim(tool) == 'HORIZCORREL_FUNCTION') then
!!$       
!!$       !
!!$       !- 7.  HORIZCORREL_FUNCTION
!!$       !
!!$       allocate(SpVertCorrel(bhi%nVarLev,bhi%nVarLev,0:ntrunc))
!!$       allocate(PowerSpectrum(bhi%nVarLev,0:ntrunc))
!!$       allocate(NormB(bhi%nVarLev,bhi%nVarLev,0:ntrunc))
!!$
!!$       !- 7.0 Normalization
!!$       call Normalize3d(ensPerturbations, & ! INOUT
!!$                        StdDev3dGridPoint)  ! IN
!!$
!!$       !- 7.1 Vertical correlations and Power Spectra
!!$       !call CalcSpectralStats(ensPerturbations,              & ! IN
!!$       !                       SpVertCorrel, PowerSpectrum,   & ! OUT
!!$       !                       NormB)                           ! OUT
!!$
!!$       !- 7.2 Compute the horizontal correlation functions
!!$       call horizCorrelFunction(NormB) ! IN
!!$
!!$       deallocate(SpVertCorrel)
!!$       deallocate(PowerSpectrum)
!!$       deallocate(NormB)
!!$
!!$     else if (trim(tool) == 'HVCORREL_LOCAL') then
!!$       !
!!$       !- 8.  
!!$       !
!!$       call Normalize3d( ensPerturbations, & ! INOUT
!!$                         StdDev3dGridPoint)  ! IN
!!$       call calcLocalCorrelations(ensPerturbations) ! IN
!!$    else
!!$      
!!$      call bmd_localizationRadii(ensPerturbations, StdDev3dGridPoint, cv_bhi, waveBandIndex_opt=1) ! IN
!!$    endif
!!$
!!$    !
!!$    !- Write the estimated pressure profiles
!!$    !
!!$    call writePressureProfiles
!!$
!!$    deallocate(Bin2d)
!!$    deallocate(ensPerturbations)
!!$    deallocate(StdDev3dGridPoint)
!!$    deallocate(ensMean3d)

    write(*,*)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'csl_toolbox: Done!'

  end subroutine csl_toolbox

  !--------------------------------------------------------------------------
  ! createLamTemplateGrids
  !--------------------------------------------------------------------------
  subroutine createLamTemplateGrids(TemplateFileName,grd_ext_x,grd_ext_y)
    use MathPhysConstants_mod, only : MPC_DEGREES_PER_RADIAN_R8
    implicit none

    character(len=*), intent(in) :: TemplateFileName
    integer         , intent(in) :: grd_ext_x
    integer         , intent(in) :: grd_ext_y

    integer :: ni_ext, nj_ext, i, j, lev, ni, nj, nk
    integer :: iun = 0
    integer :: ier, fnom, fstouv, fstfrm, fclos, fstecr

    real(8), allocatable :: Field2d(:,:)
    real(8), allocatable :: lat_ext(:)
    real(8), allocatable :: lon_ext(:)

    real(4), allocatable :: dummy2D(:,:)

    real(8) :: dlat, dlon
    real(4) :: work

    integer :: dateo,npak,status
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
    ip3      =  hco_ens%ig3
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
    !- 4. Write the vertical grid description
    !

    !- 4.1 Write the toc-toc
    status = vgd_write(vco_bhi%vgrid,iun,'fst')

    if ( status /= VGD_OK ) then
      call utl_abort('createLamTemplateGrids: ERROR with vgd_write')
    end if

    !- 4.2 Write a dummy 2D field for each MM and TH levels
    npak   = -12
    dateo  = 0
    deet   = 0
    npas   = 0
    ni     = 4
    nj     = 2
    nk     = 1
    ip2    = 0
    ip3    = 0
    typvar = 'A'
    etiket = 'VERTICALGRID'
    grtyp  = 'G'
    ig1    = 0
    ig2    = 0
    ig3    = 0
    ig4    = 0
    datyp  = 1

    allocate(dummy2D(ni,nj))
    dummy2D(:,:) = 0.0

    do lev = 1, vco_bhi%nlev_M
      ip1 = vco_bhi%ip1_M(lev)
      ier = fstecr(dummy2D, work, npak, iun, dateo, deet, npas, ni, nj, &
           nk, ip1, ip2, ip3, typvar, 'MM', etiket, grtyp,              &
           ig1, ig2, ig3, ig4, datyp, .true.)
    end do
    do lev = 1, vco_bhi%nlev_T
      ip1 = vco_bhi%ip1_T(lev)
      ier = fstecr(dummy2D, work, npak, iun, dateo, deet, npas, ni, nj, &
           nk, ip1, ip2, ip3, typvar, 'TH', etiket, grtyp,              &
           ig1, ig2, ig3, ig4, datyp, .true.)
    end do

    deallocate(dummy2D)

    !
    !- 5.  Closing the output template file
    !
    ier = fstfrm(iun)
    ier = fclos (iun)

  end subroutine createLamTemplateGrids

  !--------------------------------------------------------------------------
  ! REMOVEDOMAINMEAN
  !--------------------------------------------------------------------------
  subroutine removeDomainMean(ensPerturbations)
    implicit none

    real(4), intent(inout)  :: ensPerturbations(hco_bhi%ni,hco_bhi%nj,bhi%nVarLev,nens)

    real(8) :: domainMean, iSize

    integer :: i,j,kgdim,ens

    write(*,*) 'RemoveDomainMean: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    iSize= 1.d0/(dble(hco_bhi%ni)*dble(hco_bhi%nj))

    !$OMP PARALLEL
    !$OMP DO PRIVATE (kgdim,ens,j,i,domainMean)
    do ens = 1,nens
      do kgdim = 1,bhi%nVarLev
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
  ! calcSpectralStats
  !--------------------------------------------------------------------------
  subroutine calcSpectralStats(ensPerts,SpVertCorrel,PowerSpectrum, &
                               NormB)
    implicit none

    type(struct_ens) :: ensPerts
    real(8), intent(out)    :: SpVertCorrel(bhi%nVarLev,bhi%nVarLev,0:ntrunc)
    real(8), intent(out)    :: PowerSpectrum(bhi%nVarLev,0:ntrunc)
    real(8), intent(out)    :: NormB(bhi%nVarLev,bhi%nVarLev,0:ntrunc)

    real(8), allocatable    :: NormPowerSpectrum(:,:)
    real(8), allocatable    :: SpectralStateVar(:,:,:)
    real(8), allocatable    :: GridState(:,:,:)
    real(8), allocatable    :: SumWeight(:)

    real(4), pointer     :: ptr4d_r4(:,:,:,:)

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
    allocate( SpectralStateVar(lst_bhi%nla,lst_bhi%nphase,bhi%nVarLev) )
    allocate( GridState(hco_bhi%ni, hco_bhi%nj, bhi%nVarLev) )

    allocate(SumWeight(0:ntrunc))
    SumWeight(:) = 0.d0

    do ens = 1, nens

      !- 1.1 Extract fields from ensPerturbations
      do k = 1, bhi%nVarLev
        ptr4d_r4 => ens_getOneLev_r4(ensPerts,k)
        GridState(:,:,k) = real(ptr4d_r4(ens,1,:,:),8)
      end do

      !- 1.2 Grid Point Space -> Spectral Space
      kind = 'GridPointToSpectral'
      call lst_VarTransform( lst_bhi%id,            & ! IN
                             SpectralStateVar,      & ! OUT
                             GridState,             & ! IN
                             kind, bhi%nVarLev     )       ! IN

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
            do k2 = 1, bhi%nVarLev
              do k1 = 1, bhi%nVarLev
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
    do k = 1, bhi%nVarLev
      PowerSpectrum(k,:) = SpVertCorrel(k,k,:)
    end do

    !
    !- 2.  Calculate the Vertical Correlations in Spectral Space
    !
    !$OMP PARALLEL
    !$OMP DO PRIVATE (totwvnb,k2,k1)
    do totwvnb = 0, ntrunc
      do k2 = 1, bhi%nVarLev
        do k1 = 1, bhi%nVarLev 
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
    allocate(NormPowerSpectrum(bhi%nVarLev,0:ntrunc))

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
      do k2 = 1, bhi%nVarLev
        do k1 = 1, bhi%nVarLev
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

  end subroutine calcSpectralStats

  !--------------------------------------------------------------------------
  ! normalizePowerSpectrum
  !--------------------------------------------------------------------------
  subroutine normalizePowerSpectrum(PowerSpectrum, NormPowerSpectrum)
    implicit none

    real(8), intent(in)    :: PowerSpectrum(bhi%nVarLev,0:ntrunc)
    real(8), intent(out)   :: NormPowerSpectrum(bhi%nVarLev,0:ntrunc)

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
    do k = 1, bhi%nVarLev
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
    allocate( SpectralStateVar(lst_bhi%nla,lst_bhi%nphase,bhi%nVarLev) )
    allocate( GridState(hco_bhi%ni, hco_bhi%nj, bhi%nVarLev) )

    !- 1.2.1 Spectral transform of a delta function (at the center of the domain)
    GridState(:,:,:) = 0.d0
    GridState(hco_bhi%ni/2,hco_bhi%nj/2,:) = 1.d0

    kind = 'GridPointToSpectral'
    call lst_VarTransform( lst_bhi%id,          & ! IN
         SpectralStateVar,      & ! OUT
         GridState,             & ! IN
         kind, bhi%nVarLev     )       ! IN

    !- 1.2.2 Apply the horizontal correlation function
    !$OMP PARALLEL
    !$OMP DO PRIVATE (totwvnb,e,ila,p,k)
    do totwvnb = 0, ntrunc
      do e = 1, lst_bhi%nePerK(totwvnb)
        ila = lst_bhi%ilaFromEK(e,totwvnb)
        do p = 1, lst_bhi%nphase
          do k = 1, bhi%nVarLev
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
         kind, bhi%nVarLev )       ! IN

    !- 1.2.4 Normalize to 1
    do k = 1, bhi%nVarLev
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
  ! calcHorizScale
  !--------------------------------------------------------------------------
  subroutine calcHorizScale(HorizScale,SpCovariance)
    use MathPhysConstants_mod, only: MPC_PI_R8
    use EarthConstants_mod, only: RA
    implicit none

    real(8), intent(out) :: HorizScale(bhi%nVarLev)
    real(8), intent(in)  :: SpCovariance(bhi%nVarLev,bhi%nVarLev,0:ntrunc)

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
    do k = 1, bhi%nVarLev
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

    do var = 1, bhi%nControlVariable
      write(*,*)
      write(*,*) bhi%controlVariable(var)%nomvar(cv_bhi)
      do k = bhi%controlVariable(var)%varLevIndexStart, bhi%controlVariable(var)%varLevIndexEnd
        write(*,'(i3,2X,f9.2,2X,a2)') k, HorizScale(k)/1000.d0, 'km'
      end do
    end do

    write(*,*)
    write(*,*) 'calcHorizScale: Done!'

  end subroutine calcHorizScale

  !--------------------------------------------------------------------------
  ! calcTotVertCorrel
  !--------------------------------------------------------------------------
  subroutine calcTotVertCorrel(TotVertCorrel, SpVertCorrel, PowerSpectrum)
    implicit none

    real(8), intent(out)    :: TotVertCorrel(bhi%nVarLev,bhi%nVarLev)
    real(8), intent(in)     :: SpVertCorrel(bhi%nVarLev,bhi%nVarLev,0:ntrunc)
    real(8), intent(in)     :: PowerSpectrum(bhi%nVarLev,0:ntrunc)

    real(8), allocatable    :: TotVertCov(:,:)

    integer           :: k1, k2, totwvnb

    write(*,*)
    write(*,*) 'CalcTotVertCorrel: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- 1.  Calculate the total Normalized Covariance Matrix
    !
    allocate(TotVertCov(bhi%nVarLev,bhi%nVarLev))
    TotVertCov(:,:) = 0.d0

    do k2 = 1, bhi%nVarLev
      do k1 = 1, bhi%nVarLev
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
    do k2 = 1, bhi%nVarLev
      do k1 = 1, bhi%nVarLev
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

  end subroutine calcTotVertCorrel

  !--------------------------------------------------------------------------
  ! calcBsqrt
  !--------------------------------------------------------------------------
  subroutine calcBsqrt(Bsqrt,B)
    implicit none
    !
    !  - produce matrix B^0.5 = V D^0.5 D^t where V and D are the
    !    eigenvectors and eigenvalues of B
    !
    real(8), intent(out)   :: Bsqrt(bhi%nVarLev,bhi%nVarLev,0:ntrunc)
    real(8), intent(in)    :: B    (bhi%nVarLev,bhi%nVarLev,0:ntrunc)

    real(8), allocatable :: EigenValues(:)
    real(8), allocatable :: Work(:)

    real(8), allocatable :: EigenVectors(:,:)

    integer :: sizework, info, totwvnb, k, k1, k2 

    write(*,*)
    write(*,*) 'CalcBsqrt: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    sizework = 64 * bhi%nVarLev
    allocate(work(sizework))

    allocate(EigenValues(bhi%nVarLev))
    allocate(EigenVectors(bhi%nVarLev,bhi%nVarLev))

    !
    !-  Calculate B^0.5 for each total wave number
    !
    do totwvnb = 0, ntrunc

      EigenVectors(:,:) = B(:,:,totwvnb)

      !- Calculate EigenVectors (V) and EigenValues (D) of B matrix
      call dsyev('V','U',bhi%nVarLev,  & ! IN
           EigenVectors,   & ! INOUT
           bhi%nVarLev,         & ! IN
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

      do k1 = 1, bhi%nVarLev
        do k2 = 1, bhi%nVarLev
          Bsqrt(k1,k2,totwvnb) = sum ( EigenVectors (k1,1:bhi%nVarLev)   &
               *    EigenVectors (k2,1:bhi%nVarLev)   &
               *    sqrt(EigenValues(1:bhi%nVarLev)) )
        end do
      end do

    end do  ! total wave number

    deallocate(EigenVectors)
    deallocate(EigenValues)
    deallocate(work)

    write(*,*)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'CalcBsqrt: Done!'

  end subroutine calcBsqrt

  !--------------------------------------------------------------------------
  ! setSpVertCorrel
  !--------------------------------------------------------------------------
  subroutine setSpVertCorrel(SpVertCorrel)
    implicit none

    real(8), intent(inout) :: SpVertCorrel(bhi%nVarLev,bhi%nVarLev,0:ntrunc)

    real(8), allocatable :: KeepOrDiscard(:,:)

    integer :: totwvnb, var1, var2, k1, k2 

    write(*,*)
    write(*,*) 'SetSpVertCorrel: Starting...'

    !
    !-  Determine which bloc of the correlation matrix to keep/discard
    !
    allocate(KeepOrDiscard(bhi%nControlVariable,bhi%nControlVariable))

    !- Calculate upper half
    do var2 = 1, bhi%nControlVariable
      do var1 = 1, var2
        if (var1 == var2) then
          KeepOrDiscard(var1,var2) = 1.d0 ! Keep Auto-Correlations
        elseif( (bhi%controlVariable(var2)%nomvar(cv_bhi) == bhi%momentumControlVar(2) .and. & ! e.g. PP-CC
             bhi%controlVariable(var1)%nomvar(cv_bhi) == bhi%momentumControlVar(1)) .or. &
             (bhi%controlVariable(var2)%nomvar(cv_bhi) == 'TT'                  .and. & ! TT-(PP, QR or UU)
             bhi%controlVariable(var1)%nomvar(cv_bhi) == bhi%momentumControlVar(1)) .or. &
             (bhi%controlVariable(var2)%nomvar(cv_bhi) == 'LQ'                  .and. & ! LQ-(PP, QR or UU)
             bhi%controlVariable(var1)%nomvar(cv_bhi) == bhi%momentumControlVar(1) .and. &
             correlationLQ ) .or. &
             (bhi%controlVariable(var2)%nomvar(cv_bhi) == 'P0'                  .and. & ! P0-(PP, QR or UU)
             bhi%controlVariable(var1)%nomvar(cv_bhi) == bhi%momentumControlVar(1)) .or. &
             (bhi%controlVariable(var2)%nomvar(cv_bhi) == 'P0'                  .and. & ! P0-TT
             bhi%controlVariable(var1)%nomvar(cv_bhi) == 'TT')                  .or. &
             (bhi%controlVariable(var2)%nomvar(cv_bhi) == 'P0'                  .and. & ! P0-LQ
             bhi%controlVariable(var1)%nomvar(cv_bhi) == 'LQ'                  .and. &
             correlationLQ ) .or. &
             (bhi%controlVariable(var2)%nomvar(cv_bhi) == 'LQ'                  .and. & ! LQ-TT
             bhi%controlVariable(var1)%nomvar(cv_bhi) == 'TT'                  .and. &
             correlationLQ ) .or. &
             (bhi%controlVariable(var2)%nomvar(cv_bhi) == 'TT'                  .and. & ! TT-(CC, DD or VV)
             bhi%controlVariable(var1)%nomvar(cv_bhi) == bhi%momentumControlVar(2)) .or. &
             (bhi%controlVariable(var2)%nomvar(cv_bhi) == 'LQ'                  .and. & ! LQ-VV
             bhi%controlVariable(var1)%nomvar(cv_bhi) == bhi%momentumControlVar(2) .and. &
             correlationLQ ) .or. &
             (bhi%controlVariable(var2)%nomvar(cv_bhi) == 'P0'                  .and. & ! P0-(CC, DD or VV)
             bhi%controlVariable(var1)%nomvar(cv_bhi) == bhi%momentumControlVar(2)) ) then
          KeepOrDiscard(var1,var2) = 1.d0 ! Keep these Cross-Correlations
        else
          KeepOrDiscard(var1,var2) = 0.d0 ! Discard these Cross-Correlation
        end if
        write(*,*) var1, var2, bhi%controlVariable(var1)%nomvar(cv_bhi), bhi%controlVariable(var2)%nomvar(cv_bhi), KeepOrDiscard(var1,var2) 
      end do
    end do

    ! Symmetrize
    do var1 = 2, bhi%nControlVariable
      do var2 = 1, var1-1
        KeepOrDiscard(var1,var2) = KeepOrDiscard(var2,var1)
      end do
    end do

    !
    !- Modify the Vertical Correlation Matrix
    !
    do totwvnb = 0, ntrunc

      do var2 = 1, bhi%nControlVariable
        do var1 = 1, bhi%nControlVariable

          do k2 = bhi%controlVariable(var2)%varLevIndexStart, bhi%controlVariable(var2)%varLevIndexEnd
            do k1 = bhi%controlVariable(var1)%varLevIndexStart, bhi%controlVariable(var1)%varLevIndexEnd
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

  end subroutine setSpVertCorrel

  !--------------------------------------------------------------------------
  ! calcVertCorrel
  !--------------------------------------------------------------------------
  subroutine calcVertCorrel(VertCorrel,ensPerturbations,Bin2d,NumBins2d)
    implicit none

    real(4), intent(in)     :: ensPerturbations(hco_bhi%ni,hco_bhi%nj,bhi%nVarLev,nens)
    real(8), intent(out)    :: VertCorrel(bhi%nVarLev,bhi%nVarLev,NumBins2d)
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
          do k2 = 1, bhi%nVarLev
            do k1 = 1, bhi%nVarLev
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

  end subroutine calcVertCorrel

  !--------------------------------------------------------------------------
  ! horizCorrelFunction
  !--------------------------------------------------------------------------
  subroutine horizCorrelFunction(NormB)
    implicit none

    real(8), intent(in)    :: NormB(bhi%nVarLev,bhi%nVarLev,0:ntrunc)

    real(8), allocatable    :: SpectralStateVar(:,:,:)
    real(8), allocatable    :: GridState(:,:,:)

    integer   :: i, j, e, ila, p, k, totwvnb

    integer   :: ier, fstouv, fnom, fstfrm, fclos
    integer   :: iunstats

    character(len=24) :: kind

    allocate( SpectralStateVar(lst_bhi%nla,lst_bhi%nphase,bhi%nVarLev) )
    allocate( GridState(hco_bhi%ni, hco_bhi%nj, bhi%nVarLev) )

    !- 3.2.1 Spectral transform of a delta function (at the center of the domain)
    GridState(:,:,:) = 0.d0
    GridState(hco_bhi%ni/2,hco_bhi%nj/2,:) = 1.d0

    kind = 'GridPointToSpectral'
    call lst_VarTransform( lst_bhi%id,            & ! IN
         SpectralStateVar,      & ! OUT
         GridState,             & ! IN
         kind, bhi%nVarLev     )       ! IN

    !- 3.2.2 Apply the horizontal correlation function
    !$OMP PARALLEL
    !$OMP DO PRIVATE (totwvnb,e,ila,p,k)
    do totwvnb = 0, ntrunc
      do e = 1, lst_bhi%nePerK(totwvnb)
        ila = lst_bhi%ilaFromEK(e,totwvnb)
        do p = 1, lst_bhi%nphase
          do k = 1, bhi%nVarLev
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
         kind, bhi%nVarLev )       ! IN

    !
    !- 4.  Write to file
    !
    !iunstats = 0
    !ier    = fnom(iunstats,'./horizCorrel.fst','RND',0)
    !ier    = fstouv(iunstats,'RND')
    !call write3d(GridState,iunstats,'HORIZCORFUNC',cv_bhi) ! IN
    !call WriteTicTacToc(iunstats) ! IN
    !ier =  fstfrm(iunstats)
    !ier =  fclos (iunstats)

    deallocate(SpectralStateVar)
    deallocate(GridState)

  end subroutine horizCorrelFunction

  !--------------------------------------------------------------------------
  ! applyHorizLoc
  !--------------------------------------------------------------------------
  subroutine applyHorizLoc(NormPowerSpectrum)
    implicit none

    real(8), intent(inout) :: NormPowerSpectrum(bhi%nVarLev,0:ntrunc)

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

    allocate( SpectralStateVar(lst_bhi%nla,lst_bhi%nphase,bhi%nVarLev) )
    allocate( GridState(hco_bhi%ni, hco_bhi%nj, bhi%nVarLev) )

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
         kind, bhi%nVarLev     )       ! IN

    !- 1.2 Apply the horizontal correlation function
    !$OMP PARALLEL
    !$OMP DO PRIVATE (totwvnb,e,ila,p,k)
    do totwvnb = 0, ntrunc
      do e = 1, lst_bhi%nePerK(totwvnb)
        ila = lst_bhi%ilaFromEK(e,totwvnb)
        do p = 1, lst_bhi%nphase
          do k = 1, bhi%nVarLev
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
         kind, bhi%nVarLev )       ! IN

    !
    !- 2.  Create and apply the localization function in gridpoint space
    !
    allocate(local_length(bhi%nVarLev))
    allocate(GridStateLoc(hco_bhi%ni, hco_bhi%nj, bhi%nVarLev) )

    do var = 1, bhi%nControlVariable
      if      ( bhi%controlVariable(var)%nomvar(cv_bhi) == 'CC' .or. &
           bhi%controlVariable(var)%nomvar(cv_bhi) == 'PP' .or. &
           bhi%controlVariable(var)%nomvar(cv_bhi) == 'QR' .or. &
           bhi%controlVariable(var)%nomvar(cv_bhi) == 'DD' ) then
        hLocalize = hLocalize_wind
      else if ( bhi%controlVariable(var)%nomvar(cv_bhi) == 'TT' .or. &
           bhi%controlVariable(var)%nomvar(cv_bhi) == 'TG' .or. &
           bhi%controlVariable(var)%nomvar(cv_bhi) == 'P0' ) then
        hLocalize = hLocalize_mass
      else if ( bhi%controlVariable(var)%nomvar(cv_bhi) == 'LQ' ) then
        hLocalize = hLocalize_humidity
      else
        call utl_abort('applyHorizLoc')
      end if

      do k = bhi%controlVariable(var)%varLevIndexStart, bhi%controlVariable(var)%varLevIndexEnd
        local_length(k) = hLocalize
      end do

    end do

    call lfn_setup('FifthOrder') ! IN
    call lfn_CreateBiPerFunction( GridStateLoc,                  & ! OUT
         local_length, hco_bhi%dlon,    & ! IN
         hco_bhi%ni, hco_bhi%nj, bhi%nVarLev)  ! IN

    GridState(:,:,:) = GridStateLoc(:,:,:) * GridState(:,:,:)

    deallocate(GridStateLoc)
    deallocate(local_length)

    !
    !- 3. Create the localized normalized power spectrum
    !    
    allocate(PowerSpectrum(bhi%nVarLev,0:ntrunc))

    !- 3.1 Transform to spectral space
    kind = 'GridPointToSpectral'
    call lst_VarTransform( lst_bhi%id   ,     & ! IN
         SpectralStateVar,  & ! OUT
         GridState,         & ! IN
         kind, bhi%nVarLev )       ! IN

    !- 3.2 Compute band mean
    allocate(SumWeight(0:ntrunc))
    SumWeight(:) = 0.d0

    PowerSpectrum(:,:) = 0.d0
    do totwvnb = 0, ntrunc
      do e = 1, lst_bhi%nePerK(totwvnb)
        ila = lst_bhi%ilaFromEK(e,totwvnb)
        do p = 1, lst_bhi%nphase
          SumWeight(totwvnb) = SumWeight(totwvnb) + lst_bhi%Weight(ila)
          do k = 1, bhi%nVarLev
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
  ! applyVertLoc
  !--------------------------------------------------------------------------
  subroutine applyVertLoc(SpVertCorrel)
    implicit none

    real(8), intent(inout) :: SpVertCorrel(bhi%nVarLev,bhi%nVarLev,0:ntrunc)

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
    do var2 = 1, bhi%nControlVariable
      do var1 = 1, bhi%nControlVariable

        !-  2.2 Set the vertical length scale

        !-  2.2.1 Select the length scale for control variable 1
        if      ( bhi%controlVariable(var1)%nomvar(cv_bhi) == 'CC' .or. &
             bhi%controlVariable(var1)%nomvar(cv_bhi) == 'PP' .or. &
             bhi%controlVariable(var1)%nomvar(cv_bhi) == 'QR' .or. &
             bhi%controlVariable(var1)%nomvar(cv_bhi) == 'DD') then
          vLocalize1 = vLocalize_wind
        else if ( bhi%controlVariable(var1)%nomvar(cv_bhi) == 'TT' .or. &
             bhi%controlVariable(var1)%nomvar(cv_bhi) == 'TG' .or. &
             bhi%controlVariable(var1)%nomvar(cv_bhi) == 'P0' ) then
          vLocalize1 = vLocalize_mass
        else if ( bhi%controlVariable(var1)%nomvar(cv_bhi) == 'LQ' ) then
          vLocalize1 = vLocalize_humidity
        else
          call utl_abort('applyVertLoc')
        end if

        !-  2.2.2 Select the length scale for control variable 2
        if      ( bhi%controlVariable(var2)%nomvar(cv_bhi) == 'CC' .or. &
             bhi%controlVariable(var2)%nomvar(cv_bhi) == 'PP' .or. &
             bhi%controlVariable(var2)%nomvar(cv_bhi) == 'QR' .or. &
             bhi%controlVariable(var2)%nomvar(cv_bhi) == 'DD') then
          vLocalize2 = vLocalize_wind
        else if ( bhi%controlVariable(var2)%nomvar(cv_bhi) == 'TT' .or. &
             bhi%controlVariable(var2)%nomvar(cv_bhi) == 'TG' .or. &
             bhi%controlVariable(var2)%nomvar(cv_bhi) == 'P0' ) then
          vLocalize2 = vLocalize_mass
        else if ( bhi%controlVariable(var2)%nomvar(cv_bhi) == 'LQ' ) then
          vLocalize2 = vLocalize_humidity
        else
          call utl_abort('applyVertLoc')
        end if

        !- 2.2.3 Length scale to be use for var1-var2 correlation
        vLocalize = (vLocalize1+vLocalize2)/2.d0

        !- 2.3 Loop on vertical levels
        do k2 = bhi%controlVariable(var2)%varLevIndexStart, bhi%controlVariable(var2)%varLevIndexEnd
          do k1 = bhi%controlVariable(var1)%varLevIndexStart, bhi%controlVariable(var1)%varLevIndexEnd

            !- 2.4 Set the pressure values

            !- 2.4.1 Pressure for control-variable-1 level
            if (bhi%controlVariable(var1)%nlev /= 1) then ! variable 3D
              lev1 = k1 - bhi%controlVariable(var1)%varLevIndexStart + 1
              if (bhi%controlVariable(var1)%GridType == 'TH') then
                pres1 = pressureProfile_T(lev1)
              else
                pres1 = pressureProfile_M(lev1)
              end if
            else
              pres1 = pressureProfile_M(vco_bhi%nlev_M) ! variable 2D
            end if

            !- 2.4.2 Pressure for control-variable-2 level
            if (bhi%controlVariable(var2)%nlev /= 1) then ! variable 3D
              lev2 = k2 - bhi%controlVariable(var2)%varLevIndexStart + 1
              if (bhi%controlVariable(var2)%GridType == 'TH') then
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
  ! writeVarStats
  !--------------------------------------------------------------------------
  subroutine writeVarStats(Bsqrt,statevector_stdDev)
    implicit none

    real(8), intent(in) :: Bsqrt(bhi%nVarLev,bhi%nVarLev,0:ntrunc)
    type(struct_gsv)    :: statevector_stdDev

    integer   :: ier, fstouv, fnom, fstfrm, fclos
    integer   :: iunstats

    character(len=24) :: fileName = './bgcov.fst'

    write(*,*)
    write(*,*) 'Writing covariance statistics for VAR'

    !
    !- 1. Add the gridded std dev (and the Tic-Tac-Toc)
    !
    call gsv_writeToFile(statevector_stdDev, trim(fileName), 'STDDEV', &
                         typvar_opt = 'E', numBits_opt = 32)

    !
    !- 2. Add C^1/2
    !

    !- Opening Output file
    iunstats = 0
    ier    = fnom(iunstats,trim(fileName),'RND',0)
    ier    = fstouv(iunstats,'RND')

    !- Add Control Variable Info
    call WriteControlVarInfo(iunstats)

    !- Bsqrt
    call WriteSpVertCorrel(Bsqrt,iunstats,'ZN','B_SQUAREROOT') ! IN

    !- Closing output file
    ier =  fstfrm(iunstats)
    ier =  fclos (iunstats)

  end subroutine writeVarStats

  !--------------------------------------------------------------------------
  ! writeDiagStats
  !--------------------------------------------------------------------------
  subroutine writeDiagStats(NormB,SpVertCorrel,TotVertCorrel,statevector_mean, &
                            statevector_stdDevGridPoint,PowerSpectrum,HorizScale)
    implicit none

    real(8), intent(in) :: NormB(bhi%nVarLev,bhi%nVarLev,0:ntrunc)
    real(8), intent(in) :: SpVertCorrel(bhi%nVarLev,bhi%nVarLev,0:ntrunc)
    real(8), intent(in) :: TotVertCorrel(bhi%nVarLev,bhi%nVarLev)
    type(struct_gsv)    :: statevector_mean
    type(struct_gsv)    :: statevector_stdDevGridPoint
    real(8), intent(in) :: PowerSpectrum(bhi%nVarLev,0:ntrunc)
    real(8), intent(in) :: HorizScale(bhi%nVarLev)

    integer   :: ier, fstouv, fnom, fstfrm, fclos
    integer   :: iunstats

    character(len=24) :: fileName = './bgcov_diag.fst'

    write(*,*)
    write(*,*) 'Writing Diagnostics'

    !
    !- 1. Add the gridded mean and std dev (and the Tic-Tac-Toc)
    !
    call gsv_writeToFile(statevector_mean, trim(fileName), 'ENSMEAN', &
                         typvar_opt = 'E', numBits_opt = 32)
    call gsv_writeToFile(statevector_stdDevGridPoint, trim(fileName), 'STDDEV_GRIDP', &
                         typvar_opt = 'E', numBits_opt = 32)

    !
    !- 2. Add stats in spectral space
    !

    !- Opening Output file
    iunstats = 0
    ier    = fnom(iunstats,trim(fileName),'RND',0)
    ier    = fstouv(iunstats,'RND')

    !- Spectral Vertical Correlations
    call writeSpVertCorrel(SpVertCorrel,iunstats,'ZZ','SPVERTCORREL') ! IN

    !- Total Vertical Correlations
    call writeTotVertCorrel(TotVertCorrel,iunstats,'ZT','TTVERTCORREL') ! IN

    !- Normalized Vertical Correlations
    call writeSpVertCorrel(NormB,iunstats,'ZN','NRVERTCORREL') ! IN

    !- Power Spectrum
    call writePowerSpectrum(PowerSpectrum,iunstats,'POWERSPECT',cv_bhi) ! IN

    !- Horizontal Correlation Length scale
    call writeHorizScale(HorizScale,iunstats,'HORIZSCALE',cv_bhi) ! IN

    !- Closing output file
    ier =  fstfrm(iunstats)
    ier =  fclos (iunstats)

  end subroutine writeDiagStats

  !--------------------------------------------------------------------------
  ! writeSpVertCorrel
  !--------------------------------------------------------------------------
  subroutine writeSpVertCorrel(SpVertCorrel,iun,nomvar_in,etiket_in)
    implicit none

    real(8), intent(in) :: SpVertCorrel(bhi%nVarLev,bhi%nVarLev,0:ntrunc)
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

    allocate(work2d(bhi%nVarLev, bhi%nVarLev))

    !- Loop over Total Wavenumbers
    do totwvnb = 0, ntrunc

      npak   = -32
      dateo  = 0
      deet   = 0
      npas   = 0
      ni     = bhi%nVarLev
      nj     = bhi%nVarLev
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

  end subroutine writeSpVertCorrel

  !--------------------------------------------------------------------------
  ! writeTotVertCorrel
  !--------------------------------------------------------------------------
  subroutine writeTotVertCorrel(TotVertCorrel,iun,nomvar_in,etiket_in)
    implicit none

    real(8), intent(in) :: TotVertCorrel(bhi%nVarLev,bhi%nVarLev)
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

    allocate(workecr(bhi%nVarLev, bhi%nVarLev))

    npak   = -32
    dateo  = 0
    deet   = 0
    npas   = 0
    ni     = bhi%nVarLev
    nj     = bhi%nVarLev
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

  end subroutine writeTotVertCorrel

  !--------------------------------------------------------------------------
  ! writePowerSpectrum
  !--------------------------------------------------------------------------
  subroutine writePowerSpectrum(PowerSpectrum,iun,etiket_in,cv_type)
    implicit none

    real(8), intent(in) :: PowerSpectrum(bhi%nVarLev,0:ntrunc)
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
    do var = 1, bhi%nControlVariable

      !- Loop over vertical Levels
      do k = 1, bhi%controlVariable(var)%nlev

        npak   = -32
        dateo  = 0
        deet   = 0
        npas   = 0
        ni     = ntrunc + 1
        nj     = 1
        nk     = 1
        ip1    = bhi%controlVariable(var)%ip1(k)
        ip2    = 0
        ip3    = 0
        typvar = 'E'
        nomvar = trim(bhi%controlVariable(var)%nomvar(cv_type))
        etiket = trim(Etiket_in)
        grtyp  = 'X'
        ig1    = 0
        ig2    = 0
        ig3    = 0 
        ig4    = 0
        datyp  = 1

        !- Extract
        kgdim = bhi%controlVariable(var)%varLevIndexStart + k - 1
        workecr(:,1) = real(PowerSpectrum(kgdim,:),4)

        !- Writing 
        ier = fstecr(workecr, work, npak, iun, dateo, deet, npas, ni, nj, &
             nk, ip1, ip2, ip3, typvar, nomvar, etiket, grtyp,   &
             ig1, ig2, ig3, ig4, datyp, .true.)

      end do ! Vertical Levels

    end do ! Variables

    deallocate(workecr)

  end subroutine writePowerSpectrum

  !--------------------------------------------------------------------------
  ! writeHorizScale
  !--------------------------------------------------------------------------
  subroutine writeHorizScale(HorizScale,iun,etiket_in,cv_type)
    implicit none

    real(8), intent(in) :: HorizScale(bhi%nVarLev)
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
    do var = 1, bhi%nControlVariable

      allocate(workecr(1,1,bhi%controlVariable(var)%nlev))

      npak   = -32
      dateo  = 0
      deet   = 0
      npas   = 0
      ni     = 1
      nj     = 1
      nk     = bhi%controlVariable(var)%nlev
      ip1    = 0
      ip2    = 0
      ip3    = 0
      typvar = 'E'
      nomvar = trim(bhi%controlVariable(var)%nomvar(cv_type))
      etiket = trim(Etiket_in)
      grtyp  = 'X'
      ig1    = 0
      ig2    = 0
      ig3    = 0
      ig4    = 0
      datyp  = 1

      !- Extract
      workecr(1,1,:) = real(HorizScale(bhi%controlVariable(var)%varLevIndexStart:bhi%controlVariable(var)%varLevIndexEnd),4)

      !- Writing 
      ier = fstecr(workecr, work, npak, iun, dateo, deet, npas, ni, nj, &
           nk, ip1, ip2, ip3, typvar, nomvar, etiket, grtyp,   &
           ig1, ig2, ig3, ig4, datyp, .true.)

      deallocate(workecr)

    end do ! Variables

  end subroutine writeHorizScale

  !--------------------------------------------------------------------------
  ! writeControlVarInfo
  !--------------------------------------------------------------------------
  subroutine writeControlVarInfo(iun)
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

    character(len=4)  :: ControlModelVarnameList(bhi%nControlVariable)
    character(len=4)  :: ControlBhiVarnameList  (bhi%nControlVariable)
    character(len=2)  :: ControlVarGridTypeList (bhi%nControlVariable)
    integer           :: ControlVarNlevList     (bhi%nControlVariable)

    !
    !- 1. Gathering the info
    !

    !- Loop over Control Variables
    do var = 1, bhi%nControlVariable
      ControlModelVarnameList(var) = trim(bhi%controlVariable(var)%nomvar(cv_model))
      ControlBhiVarnameList(var)   = trim(bhi%controlVariable(var)%nomvar(cv_bhi))
      ControlVarNlevList(var)      = bhi%controlVariable(var)%nlev
      ControlVarGridTypeList(var)  = bhi%controlVariable(var)%GridType
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
    nj       =  bhi%nControlVariable
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
    ni       =  bhi%nControlVariable
    nj       =  1 
    ier = fstecr(ControlVarNlevList, work, npak, &
         iun, dateo, deet, npas, ni, nj, 1, ip1,    &
         ip2, ip3, typvar, nomvar, 'NLEV', grtyp, ig1, &
         ig2, ig3, ig4, datyp, .true.)

  end subroutine writeControlVarInfo

  !--------------------------------------------------------------------------
  ! writePressureProfiles
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

  !--------------------------------------------------------------------------
  ! calcLocalCorrelations
  !--------------------------------------------------------------------------
  subroutine calcLocalCorrelations(ensPerturbations)
    implicit none

    real(4), intent(in) :: ensPerturbations(:,:,:,:)

    real(8), allocatable :: localHorizCorrel(:,:,:)

    real(8) :: dnens

    integer :: ier
    integer :: i, j, k, ens
    integer :: blocklength_x, blocklength_y, blockpadding, nirefpoint, njrefpoint
    integer :: iref_id, jref_id, iref, jref
    integer :: imin, imax, jmin, jmax

    integer :: nulstats, ierr, fclos, fnom, nulnam, iunstats, fstouv, fstfrm

    NAMELIST /NAMHVCORREL_LOCAL/nirefpoint, njrefpoint, blockpadding

    !
    ! To compute the local horizontal correlation for some 'reference' grid point
    ! ... we assume that the ensemble grid point mean was removed and that
    !     the ensemble values were divided by the grid point std dev.
    !

    nirefpoint = 4 ! Number of reference grid point in x
    njrefpoint = 2 ! Number of reference grid point in y
    blockpadding = 4  ! Number of grid point padding between blocks (to set correlation to 0 between each block)

    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=NAMHVCORREL_LOCAL)
    write(*,nml=NAMHVCORREL_LOCAL)
    ierr = fclos(nulnam)

    blocklength_x = hco_ens%ni / nirefpoint ! Horizontal correlation will be compute blocklength x blocklength gridpoint
    ! around each reference point
    blocklength_y = hco_ens%nj / njrefpoint ! Horizontal correlation will be compute blocklength x blocklength gridpoint
    ! around each reference point

    allocate(localHorizCorrel(hco_bhi%ni,hco_bhi%nj,bhi%nVarLev))

    localHorizCorrel(:,:,:)=0.0d0

    dnens = 1.0d0/dble(nens-1)

    !$OMP PARALLEL DO PRIVATE (k,jref_id,iref_id,iref,jref,jmin,jmax,imin,imax,j,i,ens)
    do k = 1, bhi%nVarLev

      do ens = 1, nens
        do jref_id = 1, njrefpoint
          do iref_id = 1, nirefpoint
            iref = (2*iref_id-1)*blocklength_x/2
            jref = (2*jref_id-1)*blocklength_y/2
            jmin = max(jref-(blocklength_y-blockpadding)/2,1)
            jmax = min(jref+(blocklength_y-blockpadding)/2,hco_ens%nj)
            imin = max(iref-(blocklength_x-blockpadding)/2,1)
            imax = min(iref+(blocklength_x-blockpadding)/2,hco_ens%ni)
            do j = jmin, jmax
              do i = imin, imax
                localHorizCorrel(i,j,k)=localHorizCorrel(i,j,k) + &
                     ensPerturbations(i,j,k,ens) * ensPerturbations(iref,jref,k,ens)
              end do
            end do
          end do
        end do
      end do

      do j = 1, hco_ens%nj
        do i = 1, hco_ens%ni
          localHorizCorrel(i,j,k) = localHorizCorrel(i,j,k)*dnens
        end do
      end do

    end do
    !$OMP END PARALLEL DO

    write(6,*) 'finished computing the local horizontal correlations...'
    call flush(6)

    !
    !- 4.  Write to file
    !
    !iunstats = 0
    !ier    = fnom(iunstats,'./horizCorrelLocal.fst','RND',0)
    !ier    = fstouv(iunstats,'RND')
    !call write3d(localHorizCorrel,iunstats,'HCORREL_LOC',cv_bhi)
    !call WriteTicTacToc(iunstats) ! IN
    !ier =  fstfrm(iunstats)
    !ier =  fclos (iunstats)

    deallocate(localHorizCorrel)

  end subroutine calcLocalCorrelations

end module calcstatslam_mod
