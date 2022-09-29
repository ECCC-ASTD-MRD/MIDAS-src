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

module calcStatsLam_mod
  ! MODULE calcStatsLam_mod (prefix='csl' category='1. High-level functionality')
  !
  ! :Purpose: To compute homogeneous and isotropic background error covariances 
  !           from forecast error estimate in model variable space (limited-area
  !           version).
  !
  use MathPhysConstants_mod
  use earthConstants_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use lamSpectralTransform_mod
  use analysisGrid_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use localizationFunction_mod
  use utilities_mod
  use menetrierDiag_mod
  use ensemblestatevector_mod
  use gridVariableTransforms_mod
  use varNameList_mod
  use gridBinning_mod
  use timeCoord_mod
  use midasMpi_mod
  implicit none
  save
  private

  ! Public Subroutines
  public :: csl_setup, csl_computeBhi
  public :: csl_toolbox

  type(struct_hco), pointer :: hco_ens => null() ! Ensemble horizontal grid parameters
  type(struct_hco), pointer :: hco_bhi => null() ! B matrix horizontal grid parameters
  type(struct_vco), pointer :: vco_bhi => null() ! B matrix vertical grid parameters

  integer,external   :: get_max_rss

  integer, parameter :: cv_model = 1
  integer, parameter :: cv_bhi   = 2
  integer, parameter :: nMaxControlVar = 10
  integer, parameter :: maxNumLevels   = 200
  
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

  integer :: nEns
  integer :: ip2_ens

  logical :: initialized = .false.
  logical :: vertLoc, horizLoc, stdDevScaling
  logical :: ensContainsFullField

  character(len=2)  :: ctrlVarHumidity

  type(struct_ens)  :: ensPerts

  real(8), pointer  :: pressureProfile_M(:), pressureProfile_T(:)


  real(8),allocatable :: scaleFactor_M(:), scaleFactor_T(:)
  real(8)             :: scaleFactor_SF

  ! Namelist variables
  integer :: nTrunc
  character(len=12) :: WindTransform
  logical :: NormByStdDev
  logical :: SetTGtoZero
  logical :: writeEnsPert
  character(len=12) :: SpectralWeights
  real(8) :: vLocalize_wind     ! vertical length scale (in units of ln(Pressure))
  real(8) :: vlocalize_mass     ! vertical length scale (in units of ln(Pressure))
  real(8) :: vlocalize_humidity ! vertical length scale (in units of ln(Pressure))
  real(8) :: vlocalize_other    ! vertical length scale (in units of ln(Pressure))
  real(8) :: hlocalize_wind     ! horizontal length scale (in km)
  real(8) :: hlocalize_mass     ! horizontal length scale (in km)
  real(8) :: hlocalize_humidity ! horizontal length scale (in km)
  real(8) :: hlocalize_other    ! horizontal length scale (in km)
  character(len=4)  :: correlatedVariables(vnl_numvarmax)

contains
  
  !--------------------------------------------------------------------------
  ! csl_setup
  !--------------------------------------------------------------------------
  subroutine csl_setup(nEns_in, hco_ens_in, vco_ens_in, ip2_in)
    !
    ! :Purpose: To initialize this module
    !
    use vGrid_Descriptors , only: vgrid_descriptor, vgd_levels, VGD_OK  
    implicit none

    ! arguments
    integer,                   intent(in)   :: nEns_in
    type(struct_vco), pointer, intent(in)   :: vco_ens_in
    type(struct_hco), pointer, intent(in)   :: hco_ens_in
    integer,                   intent(in)   :: ip2_in

    ! locals
    integer :: nulnam, ier, status
    integer :: fclos, fnom
    integer :: varIndex, levIndex, k
    integer :: numStep
    integer, allocatable :: dateStampList(:)
    real(8) :: SurfacePressure
    character(len=256)  :: enspathname
    logical :: makeBiPeriodic
    character(len=4), pointer :: controlVarNames(:)

    ! Namelist variables (local)
    real(8) :: scaleFactor(maxNumLevels)
    integer :: grd_ext_x
    integer :: grd_ext_y

    NAMELIST /NAMCALCSTATS_LAM/nTrunc,grd_ext_x,grd_ext_y,WindTransform,NormByStdDev, &
                               SetTGtoZero,writeEnsPert,SpectralWeights,              &
                               vLocalize_wind,vlocalize_mass,vlocalize_humidity,      &
                               hLocalize_wind,hlocalize_mass,hlocalize_humidity,      &
                               hLocalize_other,vlocalize_other,                       &
                               correlatedVariables,scaleFactor

    write(*,*)
    write(*,*) 'csl_setup: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- 1. Initialized the info on the ensemble
    !
    nEns=nEns_in

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
    nTrunc        = 75       ! Default value
    grd_ext_x     = 10       ! Default value
    grd_ext_y     = 10       ! Default value
    WindTransform = 'PsiChi' ! Default value
    NormByStdDev  = .true.   ! Default value
    SetTGtoZero   = .false.  ! Default value
    writeEnsPert  = .false.  ! Default value
    correlatedVariables = '    '
    SpectralWeights = 'lst'  ! Default value
    vLocalize_wind      = -1.0d0 ! Default value (no vloc)
    vLocalize_mass      = -1.0d0 ! Default value (no vloc)
    vLocalize_humidity  = -1.0d0 ! Default value (no vloc)
    vLocalize_other     = -1.0d0 ! Default value (no vloc)
    hLocalize_wind      = -1.0d0 ! Default value (no hloc)
    hLocalize_mass      = -1.0d0 ! Default value (no hloc)
    hLocalize_humidity  = -1.0d0 ! Default value (no hloc)
    hLocalize_other     = -1.0d0 ! Default value (no hloc)
    scaleFactor(:)      =  1.0d0 ! Default value (no scaling)
    
    nulnam = 0

    ier=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read (nulnam,nml=namcalcstats_lam)
    write(*     ,nml=namcalcstats_lam)
    ier=fclos(nulnam)

    write(*,*)
    write(*,*) 'Truncation = ', nTrunc

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
    if (mmpi_myid == 0) then
      call agd_createLamTemplateGrids('./analysisgrid', hco_ens, vco_bhi, & ! IN
                                      grd_ext_x, grd_ext_y)                 ! IN
    end if
    call rpn_comm_barrier("GRID",ier)

    !- 3.2 Setup the Extended B_HI grid
    call hco_setupFromFile( hco_bhi,'./analysisgrid', 'ANALYSIS', 'BHI' ) ! IN

    !- 3.3 Setup the LAM analysis grid metrics
    call agd_setupFromHCO( hco_bhi, hco_ens ) ! IN

    !
    !- 4. Read the ensemble
    !
    numStep = 1
    allocate(dateStampList(numStep))
    dateStampList(:)  = -1
    call ens_allocate(ensPerts, nEns, numStep, hco_bhi, vco_bhi, dateStampList, &
                      hco_core_opt=hco_ens)

    ensContainsFullField = .false.
    ctrlVarHumidity = 'LQ'
    ensPathName    = './ensemble'
    makeBiPeriodic = .true.
    call ens_readEnsemble(ensPerts, ensPathName, makeBiPeriodic, &
                          containsFullField_opt=ensContainsFullField)

    if ( ctrlVarHumidity == 'LQ' .and. ens_varExist(ensPerts,'HU') .and. &
         ensContainsFullField ) then
      call gvt_transform(ensPerts,'HUtoLQ')
    end if

    !
    !- 5.  Setup the control variables (model space and B_hi space)
    !
    nullify(controlVarNames)
    call ens_varNamesList(controlVarNames,ensPerts)

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
      call ens_writeEnsemble(ensPerts, './', 'MODELVAR_', 'MODELVAR', 'E', &
                             containsFullField_opt = ensContainsFullField)
    end if

    if      ( bhi%momentumControlVar(1) == 'PP' .and. bhi%momentumControlVar(2) == 'CC' ) then
      call gvt_transform(ensPerts,'UVtoPsiChi')
      call ens_modifyVarName(ensPerts, 'UU', 'PP')
      call ens_modifyVarName(ensPerts, 'VV', 'CC')
    else if ( bhi%momentumControlVar(1) == 'QR' .and. bhi%momentumControlVar(2) == 'DD' ) then
      call gvt_transform(ensPerts,'UVtoVortDiv')
      call ens_modifyVarName(ensPerts, 'UU', 'QR')
      call ens_modifyVarName(ensPerts, 'VV', 'DD')
    end if

    if (writeEnsPert) then
      call ens_writeEnsemble(ensPerts, './', 'CTRLVAR', 'CTRLVAR', 'E', &
                             containsFullField_opt = ensContainsFullField)
    end if

    !
    !- 7.  Compute and remove the ensemble mean; compute the stdDev
    !
    call ens_computeMean(ensPerts)
    call ens_removeMean (ensPerts)
    call ens_computeStdDev(ensPerts)

    !
    !- 8.  Setup the localization 
    !

    !- 8.1 Setup horizontal localization
    if ( hLocalize_wind     < 0.d0 .and. hLocalize_mass  < 0.d0 .and. &
         hLocalize_humidity < 0.d0 .and. hLocalize_other < 0.d0) then
      write(*,*) 
      write(*,*) 'csl_setup: NO horizontal correlation localization will be performed'
      horizLoc=.false.
    else
      write(*,*) 
      write(*,*) 'csl_setup: horizontal correlation localization WILL BE performed'
      horizLoc=.true.
    end if

    !- 8.2 Setup vertical localization
    if ( vLocalize_wind     < 0.d0 .and. vLocalize_mass  < 0.d0 .and. &
         vLocalize_humidity < 0.d0 .and. vLocalize_other < 0.d0 ) then
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
    if (vco_bhi%vgridPresent) then
      SurfacePressure = 101000.D0

      status = vgd_levels(vco_bhi%vgrid, ip1_list=vco_bhi%ip1_M,     & ! IN
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
      end if
      
      status = vgd_levels(vco_bhi%vgrid, ip1_list=vco_bhi%ip1_T,     & ! IN
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
      end if

    end if

    !
    !- 10.  Setup the scaling
    !
    if ( all(scaleFactor(:) == 1.d0) ) then
      write(*,*) 
      write(*,*) 'csl_setup: NO scaling of the StdDev will be performed'
      stdDevScaling=.false.

    else
      write(*,*) 
      write(*,*) 'csl_setup: scaling of the StdDev WILL BE performed'
      stdDevScaling=.true.

      allocate(scaleFactor_M(vco_bhi%nlev_M))
      allocate(scaleFactor_T(vco_bhi%nlev_T))
      do levIndex = 1, vco_bhi%nlev_T
        if (scaleFactor(levIndex) > 0.0d0) then 
          scaleFactor_T(levIndex) = sqrt(scaleFactor(levIndex))
        else
          scaleFactor_T(levIndex) = 0.0d0
        end if
      end do
      scaleFactor_M(1:vco_bhi%nlev_M) = scaleFactor_T(1:vco_bhi%nlev_M)
      scaleFactor_SF = scaleFactor_T(vco_bhi%nlev_T)

    end if

    !
    !- 11.  Ending
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
    !
    ! :Purpose: To compute an homogeneous and isotopic B matrix
    !
    implicit none

    ! locals
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

    integer :: ier

    write(*,*)
    write(*,*) 'csl_computeBhi: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    allocate(SpVertCorrel(bhi%nVarLev,bhi%nVarLev,0:nTrunc))
    allocate(TotVertCorrel(bhi%nVarLev,bhi%nVarLev))
    allocate(PowerSpectrum(bhi%nVarLev,0:nTrunc))
    allocate(NormB(bhi%nVarLev,bhi%nVarLev,0:nTrunc))
    allocate(NormBsqrt(bhi%nVarLev,bhi%nVarLev,0:nTrunc))
    allocate(HorizScale(bhi%nVarLev))

    !
    !- 1.  Calculate the gridded binned Std. Dev. to be used in the analysis step
    !
    nullify(varNamesList)
    call ens_varNamesList(varNamesList,ensPerts) 
    call gsv_allocate(statevector_stdDev, ens_getNumStep(ensPerts),                            &
                      ens_getHco(ensPerts), ens_getVco(ensPerts), varNames_opt=varNamesList, &
                      datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true.,      &
                      dataKind_opt=8 )

    call gbi_setup(gbi_horizontalMean, 'HorizontalMean', statevector_stdDev, hco_ens)

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

     if (mmpi_myid == 0) then
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
     end if

     !- 3.6 Apply scaling
     if (stdDevScaling) then
       call scaleStdDev(statevector_stdDev) ! INOUT
     end if
     
     call rpn_comm_barrier("GRID",ier)

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
    !
    ! :Purpose: High-level control of various diagnostic tools
    !
    implicit none

    ! locals
    real(8),allocatable :: SpVertCorrel(:,:,:)
    real(8),allocatable :: NormB(:,:,:)
    real(8),allocatable :: PowerSpectrum(:,:)

    integer :: nulnam, ier, fnom, fclos

    type(struct_gsv) :: statevector_stdDev
    type(struct_gsv) :: statevector_template
    character(len=4), pointer :: varNamesList(:)

    character(len=60) :: tool

    NAMELIST /NAMTOOLBOX/tool

    write(*,*)
    write(*,*) 'csl_toolbox: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- 1.  Tool selection
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

       call ens_normalize(ensPerts)

       call calcVertCorrel(ensPerts)

    case ('HORIZCORREL_FUNCTION')
       write(*,*)
       write(*,*) 'Computing horizontal correlation functions'

       allocate(SpVertCorrel(bhi%nVarLev,bhi%nVarLev,0:nTrunc))
       allocate(PowerSpectrum(bhi%nVarLev,0:nTrunc))
       allocate(NormB(bhi%nVarLev,bhi%nVarLev,0:nTrunc))

       call ens_normalize(ensPerts)

       call calcSpectralStats(ensPerts,                      & ! IN
                              SpVertCorrel, PowerSpectrum,   & ! OUT
                              NormB)                           ! OUT

       if (mmpi_myid == 0) then
         call horizCorrelFunction(NormB) ! IN
       end if

       deallocate(NormB)
       deallocate(PowerSpectrum)
       deallocate(SpVertCorrel)

    case ('STDDEV')
       write(*,*)
       write(*,*) 'Computing Standard-Deviations'
       
       nullify(varNamesList)
       call ens_varNamesList(varNamesList,ensPerts) 
       call gsv_allocate(statevector_stdDev, ens_getNumStep(ensPerts),                            &
                         ens_getHco(ensPerts), ens_getVco(ensPerts), varNames_opt=varNamesList, &
                         datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true.,      &
                         dataKind_opt=8 )

       call ens_copyEnsStdDev(ensPerts, statevector_stdDev)

       call gio_writeToFile(statevector_stdDev, './stddev.fst', 'STDDEV_GRIDP', &
                            typvar_opt = 'E', numBits_opt = 32)

       call gsv_deallocate(statevector_stdDev)

    case ('HVCORREL_LOCAL')
       write(*,*)
       write(*,*) 'Computing Local Correlation'
       call ens_normalize(ensPerts)
       call calcLocalCorrelations(ensPerts)

     case ('LOCALIZATIONRADII')
       write(6,*)
       write(6,*) 'Estimating the optimal covariance localization radii'

       call ens_copyEnsStdDev(ensPerts, statevector_template)

       call bmd_setup(statevector_template, hco_ens, nEns, pressureProfile_M, pressureProfile_T, 1)

       call bmd_localizationRadii(ensPerts, waveBandIndex_opt=1) ! IN

    case default
       call utl_abort('csl_toolbox: Unknown TOOL '// trim(tool))
    end select

    !
    !- 2.  Write the estimated pressure profiles
    !
    if (mmpi_myid == 0 .and. vco_bhi%vgridPresent) then
      call writePressureProfiles
    end if

    !
    !- 3.  Ending
    !
    call ens_deallocate(ensPerts)

    write(*,*)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'csl_toolbox: Done!'

  end subroutine csl_toolbox

  !--------------------------------------------------------------------------
  ! calcSpectralStats
  !--------------------------------------------------------------------------
  subroutine calcSpectralStats(ensPerts,SpVertCorrel,PowerSpectrum, &
                               NormB)
    !
    ! :Purpose: To compute background-error covariances in spectral space
    !           from an ensemble of gridded data
    !
    implicit none

    ! arguments
    type(struct_ens)        :: ensPerts
    real(8), intent(out)    :: SpVertCorrel(bhi%nVarLev,bhi%nVarLev,0:nTrunc)
    real(8), intent(out)    :: PowerSpectrum(bhi%nVarLev,0:nTrunc)
    real(8), intent(out)    :: NormB(bhi%nVarLev,bhi%nVarLev,0:nTrunc)

    ! locals
    real(8), allocatable    :: NormPowerSpectrum(:,:)
    real(8), allocatable    :: SpectralStateVar(:,:,:)
    real(8), allocatable    :: SpVertCorrel_local(:,:,:)
    real(8), allocatable    :: GridState(:,:,:)
    real(8), allocatable    :: SumWeight(:)
    real(8), allocatable    :: SumWeight_local(:)

    type(struct_lst)  :: lst_bhi ! Spectral transform Parameters

    real(4), pointer  :: ptr4d_r4(:,:,:,:)

    real(8)           :: weight

    integer           :: k1, k2, ens, e, ila, p, k, totwvnb
    integer           :: myLonBeg, myLonEnd, myLatBeg, myLatEnd
    integer           :: nSize, ier

    character(len=24) :: kind

    write(*,*)
    write(*,*) 'CalcSpectralStats: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- 1. Setup the (LAM) spectral transform
    !

    call lst_setup(lst_bhi,                             & ! OUT
                   hco_bhi%ni, hco_bhi%nj,              & ! IN
                   hco_bhi%dlon, nTrunc,                & ! IN
                   'LatLonMN', maxlevels_opt=bhi%nVarLev) ! IN

    !
    !- 2.  Calculate the Vertical Covariances in Spectral Space
    !
    call ens_getLatLonBounds(ensPerts, myLonBeg, myLonEnd, myLatBeg, myLatEnd)

    allocate( SpVertCorrel_local(bhi%nVarLev, bhi%nVarLev, 0:nTrunc) )
    SpVertCorrel_local(:,:,:) = 0.d0

    allocate( SpectralStateVar(lst_bhi%nla,lst_bhi%nphase,bhi%nVarLev) )
    allocate( GridState(myLonBeg:myLonEnd, myLatBeg:myLatEnd, bhi%nVarLev) )

    allocate(SumWeight_local(0:nTrunc))
    SumWeight_local(:) = 0.d0

    do ens = 1, nEns

      write(6,*) 'ens = ', ens
      flush(6)

      !- 2.1 Extract fields from ensPerturbations
      do k = 1, bhi%nVarLev
        ptr4d_r4 => ens_getOneLev_r4(ensPerts,k)
        GridState(:,:,k) = real(ptr4d_r4(ens,1,:,:),8)
      end do

      !- 2.2 Grid Point Space -> Spectral Space
      kind = 'GridPointToSpectral'
      call lst_VarTransform(lst_bhi,               & ! IN
                            SpectralStateVar,      & ! OUT
                            GridState,             & ! IN
                            kind, bhi%nVarLev)       ! IN

      !- 2.3 Compute the covariances
      !$OMP PARALLEL DO PRIVATE(totwvnb,weight,e,ila,p,k2,k1)
      do totwvnb = 0, nTrunc
        do e = 1, lst_bhi%nePerK(totwvnb)
          ila = lst_bhi%ilaFromEK(e,totwvnb)
          do p = 1, lst_bhi%nphase
            if (trim(SpectralWeights) == 'lst') then
              weight = lst_bhi%Weight(ila)
              SumWeight_local(totwvnb) = SumWeight_local(totwvnb) + weight
            else
              weight = 2.0d0
            end if
            do k2 = 1, bhi%nVarLev
              do k1 = 1, bhi%nVarLev
                SpVertCorrel_local(k1,k2,totwvnb) = SpVertCorrel_local(k1,k2,totwvnb) &
                     + weight * (SpectralStateVar(ila,p,k1) * SpectralStateVar(ila,p,k2))
              end do
            end do
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    end do ! Loop in Ensemble

    deallocate(SpectralStateVar)
    deallocate(GridState)

    ! Gather the all the info in processor 0
    SpVertCorrel(:,:,:) = 0.d0
    nSize = bhi%nVarLev * bhi%nVarLev * (nTrunc + 1)
    call rpn_comm_reduce(SpVertCorrel_local,SpVertCorrel,nsize,"mpi_double_precision","mpi_sum",0,"GRID",ier)

    allocate(SumWeight(0:nTrunc))
    SumWeight(:) = 0.d0
    nSize = nTrunc + 1
    call rpn_comm_reduce(SumWeight_local,   SumWeight,   nsize,"mpi_double_precision","mpi_sum",0,"GRID",ier)

    deallocate(SumWeight_local)
    deallocate(SpVertCorrel_local)

    if (mmpi_myid == 0) then 

      !- 2.4 Compute the weighted COVARIANCES for each total wavenumber
      do totwvnb = 0, nTrunc
        if (trim(SpectralWeights) == 'legacy') then
          if (totwvnb /= 0 ) then
            SumWeight(totwvnb) = 2.d0 * real(lst_bhi%nphase * lst_bhi%nePerKglobal(totwvnb),8) - 2.d0
          else
            SumWeight(totwvnb) = 1.d0
          end if
          SumWeight(totwvnb) = SumWeight(totwvnb)*nEns
        end if
        
        if ( SumWeight(totwvnb) /= 0.d0 ) then 
          SpVertCorrel(:,:,totwvnb) = SpVertCorrel(:,:,totwvnb) / SumWeight(totwvnb)
        else
          SpVertCorrel(:,:,totwvnb) = 0.d0
        end if
        
      end do

      deallocate(SumWeight)
      
      !- 2.5 Extract the power spectrum (the variances on the diagonal elements)
      do k = 1, bhi%nVarLev
        PowerSpectrum(k,:) = SpVertCorrel(k,k,:)
      end do

      !
      !- 3.  Calculate the Vertical Correlations in Spectral Space
      !
      !$OMP PARALLEL DO PRIVATE (totwvnb,k2,k1)
      do totwvnb = 0, nTrunc
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
      !$OMP END PARALLEL DO
      
      ! Apply vertical localization (if wanted)
      if (vertLoc) then                                                                 
        call applyVertLoc(SpVertCorrel) ! INOUT                                        
      end if
      
      !
      !- 4.  Normalize the power spectrum (i.e. build normalised spectral densities of the variance)
      !
      allocate(NormPowerSpectrum(bhi%nVarLev,0:nTrunc))
      
      call NormalizePowerSpectrum(PowerSpectrum,     & ! IN
           NormPowerSpectrum)   ! OUT
      
      ! Apply horizontal localization (if wanted)
      if (horizLoc) then                                                                 
        call applyHorizLoc(NormPowerSpectrum) ! INOUT                                        
      end if
      
      !
      !- 5.  Normalize the spectral vertical correlation matrix to ensure correlations in horizontal
      !
      
      !$OMP PARALLEL DO PRIVATE (totwvnb,k2,k1)
      do totwvnb = 0, nTrunc
        do k2 = 1, bhi%nVarLev
          do k1 = 1, bhi%nVarLev
            NormB(k1,k2,totwvnb) = SpVertCorrel(k1,k2,totwvnb) * &
                 sqrt( NormPowerSpectrum(k1,totwvnb) * NormPowerSpectrum(k2,totwvnb) )
          end do
        end do
      end do
      !$OMP END PARALLEL DO
      
      deallocate(NormPowerSpectrum)

    end if ! mmpi_myid == 0

    write(*,*)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'CalcSpectralStats: Done!'

  end subroutine calcSpectralStats

  !--------------------------------------------------------------------------
  ! normalizePowerSpectrum
  !--------------------------------------------------------------------------
  subroutine normalizePowerSpectrum(PowerSpectrum, NormPowerSpectrum)
    !
    ! :Purpose: To convert spectral variances into spectral correlations 
    !
    implicit none

    ! arguments
    real(8), intent(in)    :: PowerSpectrum(bhi%nVarLev,0:nTrunc)
    real(8), intent(out)   :: NormPowerSpectrum(bhi%nVarLev,0:nTrunc)

    ! locals
    real(8), allocatable   :: SpectralStateVar(:,:,:)
    real(8), allocatable   :: GridState(:,:,:)

    real(8)           :: sum

    integer           :: e, ila, p, k, totwvnb

    character(len=24) :: kind

    type(struct_lst)  :: lst_norm ! Spectral transform Parameters

    write(*,*)
    write(*,*) 'NormalizePowerSpectrum: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- 1. Setup the (LAM) spectral transform
    !
    call lst_setup(lst_norm,                            & ! OUT
                   hco_bhi%ni, hco_bhi%nj,              & ! IN
                   hco_bhi%dlon, nTrunc,                & ! IN
                   'NoMpi')

    !
    !- 1.  Normalize the power spectrum (i.e. build normalised spectral densities of the variance)
    !

    !- 1.1 Part 1

    !$OMP PARALLEL DO PRIVATE (totwvnb,k,sum)
    do k = 1, bhi%nVarLev
      sum = 0.0d0
      do totwvnb = 0, nTrunc
        sum = sum + real(totwvnb,8) * PowerSpectrum(k,totwvnb)
      end do
      do totwvnb = 0, nTrunc
        if ( sum /= 0.0d0 ) then
          NormPowerSpectrum(k,totwvnb) = PowerSpectrum(k,totwvnb) / sum
        else
          NormPowerSpectrum(k,totwvnb) = 0.d0
        end if
      end do
    end do
    !$OMP END PARALLEL DO

    !- 1.2 Part 2
    allocate( SpectralStateVar(lst_norm%nla,lst_norm%nphase,bhi%nVarLev) )
    allocate( GridState(hco_bhi%ni, hco_bhi%nj, bhi%nVarLev) )

    !- 1.2.1 Spectral transform of a delta function (at the center of the domain)
    GridState(:,:,:) = 0.d0
    GridState(hco_bhi%ni/2,hco_bhi%nj/2,:) = 1.d0

    kind = 'GridPointToSpectral'
    call lst_VarTransform(lst_norm,           & ! IN
                          SpectralStateVar,   & ! OUT
                          GridState,          & ! IN
                          kind, bhi%nVarLev)    ! IN

    !- 1.2.2 Apply the horizontal correlation function
    !$OMP PARALLEL DO PRIVATE (totwvnb,e,ila,p,k)
    do totwvnb = 0, nTrunc
      do e = 1, lst_norm%nePerK(totwvnb)
        ila = lst_norm%ilaFromEK(e,totwvnb)
        do p = 1, lst_norm%nphase
          do k = 1, bhi%nVarLev
            SpectralStateVar(ila,p,k) = SpectralStateVar(ila,p,k) * NormPowerSpectrum(k,totwvnb) * &
                 lst_norm%NormFactor(ila,p) * lst_norm%NormFactorAd(ila,p)
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    !- 1.2.3 Move back to physical space
    kind = 'SpectralToGridPoint'
    call lst_VarTransform(lst_norm,          & ! IN
                          SpectralStateVar,  & ! IN
                          GridState,         & ! OUT
                          kind, bhi%nVarLev)   ! IN

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

  end subroutine normalizePowerSpectrum

  !--------------------------------------------------------------------------
  ! calcHorizScale
  !--------------------------------------------------------------------------
  subroutine calcHorizScale(HorizScale,SpCovariance)
    !
    ! :Purpose: To compute horizontal lenght scales based on the power spectra
    !
    implicit none

    ! arguments
    real(8), intent(out) :: HorizScale(bhi%nVarLev)
    real(8), intent(in)  :: SpCovariance(bhi%nVarLev,bhi%nVarLev,0:nTrunc)

    ! locals
    real(8) :: a, b, beta, dx, dist

    integer :: totwvnb, k, var

    write(*,*)
    write(*,*) 'CalcHorizScale: Starting...'

    !
    !- Computing distance-related variables
    !

    ! Grid spacing in meters
    dx = hco_bhi%dlon * ec_ra
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
      do totwvnb = 0, nTrunc
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
    !
    ! :Purpose: To compute the total vertical correlations (i.e. gridpoint equivalent)
    !
    implicit none

    ! arguments
    real(8), intent(out)    :: TotVertCorrel(bhi%nVarLev,bhi%nVarLev)
    real(8), intent(in)     :: SpVertCorrel(bhi%nVarLev,bhi%nVarLev,0:nTrunc)
    real(8), intent(in)     :: PowerSpectrum(bhi%nVarLev,0:nTrunc)

    ! locals
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
        do totwvnb = 0, nTrunc
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
    !
    ! :Purpose: To compute the sqare-root of B
    !
    implicit none

    ! arguments
    real(8), intent(out)   :: Bsqrt(bhi%nVarLev,bhi%nVarLev,0:nTrunc)
    real(8), intent(in)    :: B    (bhi%nVarLev,bhi%nVarLev,0:nTrunc)

    ! locals
    integer :: totwvnb

    !
    !-  Calculate B^0.5 for each total wave number
    !
    Bsqrt(:,:,:) = B(:,:,:)

    do totwvnb = 0, nTrunc
      call utl_matSqrt(Bsqrt(:,:,totwvnb),bhi%nVarLev,1.0d0,.true.)
    end do

  end subroutine calcBsqrt

  !--------------------------------------------------------------------------
  ! setSpVertCorrel
  !--------------------------------------------------------------------------
  subroutine setSpVertCorrel(SpVertCorrel)
    !
    ! :Purpose: To discard some user-defined cross-correlations
    !
    implicit none

    ! arguments
    real(8), intent(inout) :: SpVertCorrel(bhi%nVarLev,bhi%nVarLev,0:nTrunc)

    ! locals
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
          KeepOrDiscard(var1,var2) = 1.d0 ! Keep the Auto-Correlations
        elseif( any(correlatedVariables == bhi%controlVariable(var1)%nomvar(cv_bhi)) .and. &
                any(correlatedVariables == bhi%controlVariable(var2)%nomvar(cv_bhi)) ) then
          KeepOrDiscard(var1,var2) = 1.d0 ! Keep these Cross-Correlations
        else
          KeepOrDiscard(var1,var2) = 0.d0 ! Discard these Cross-Correlations
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
    do totwvnb = 0, nTrunc

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
  subroutine calcVertCorrel(ensPerts)
    !
    ! :Purpose: To compute vertical correlations from an ensemble of gridded data
    !
    implicit none

    ! arguments
    type(struct_ens)     :: ensPerts

    ! locals
    real(8), allocatable :: vertCorrel(:,:)

    real(4), pointer  :: ptr4d_k1_r4(:,:,:,:)
    real(4), pointer  :: ptr4d_k2_r4(:,:,:,:)

    real(8), allocatable :: vertCorrel_local(:,:)

    integer :: lonIndex, latIndex, k1, k2, memberIndex
    integer :: myLonBeg, myLonEnd, myLatBeg, myLatEnd, nSize, ier
    integer :: fstouv, fnom, fstfrm, fclos, iunstats

    write(*,*)
    write(*,*) 'CalcVertCorrel: Starting...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    allocate(vertCorrel(bhi%nVarLev,bhi%nVarLev))
    vertCorrel(:,:) = 0.d0

    allocate(vertCorrel_local(bhi%nVarLev,bhi%nVarLev))
    vertCorrel_local(:,:) = 0.d0

    !
    !- Calculate the Vertical Correlation in GridPoint Space
    !  ... we assume that the ensemble grid point mean was removed and that
    !      the ensemble values were divided by the grid point std dev.

    call ens_getLatLonBounds(ensPerts, myLonBeg, myLonEnd, myLatBeg, myLatEnd)

    do memberIndex = 1, nEns

      !$OMP PARALLEL DO PRIVATE (k1,k2,ptr4d_k1_r4,ptr4d_k2_r4,latIndex,lonIndex)
      do k2 = 1, bhi%nVarLev
        do k1 = 1, bhi%nVarLev
          
          ptr4d_k1_r4 => ens_getOneLev_r4(ensPerts,k1)
          ptr4d_k2_r4 => ens_getOneLev_r4(ensPerts,k2)

          do latIndex = myLatBeg, myLatEnd
            do lonIndex = myLonBeg, myLonEnd

              if (lonIndex > hco_ens%ni .or. latIndex > hco_ens%nj) cycle ! do not use data in the extension zone

              VertCorrel_local(k1,k2) = VertCorrel_local(k1,k2)            &
                   + real(ptr4d_k1_r4(memberIndex,1,lonIndex,latIndex),8)  &
                   * real(ptr4d_k2_r4(memberIndex,1,lonIndex,latIndex),8)
            end do
          end do

        end do
      end do
      !$OMP END PARALLEL DO

    end do ! Loop in Ensemble

    !- Communication
    nSize = bhi%nVarLev * bhi%nVarLev
    call rpn_comm_reduce(vertCorrel_local,vertCorrel,nsize,"mpi_double_precision","mpi_sum",0,"GRID",ier)

    deallocate(vertCorrel_local)

    !- Conversion to correlation
    if (mmpi_myid == 0) then
      vertCorrel(:,:) = vertCorrel(:,:) / real((nEns-1)*hco_ens%nj*hco_ens%ni,8)
    end if

    !- Output
    if (mmpi_myid == 0) then
      iunstats = 0
      ier    = fnom(iunstats,'./vertCorrel.fst','RND',0)
      ier    = fstouv(iunstats,'RND')
      call WriteTotVertCorrel(VertCorrel,iunstats,'ZT','GPVCOR_CORE') ! IN
      ier =  fstfrm(iunstats)
      ier =  fclos (iunstats)
    end if

    deallocate(VertCorrel)

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'CalcVertCorrel: Done!'

  end subroutine calcVertCorrel

  !--------------------------------------------------------------------------
  ! horizCorrelFunction
  !--------------------------------------------------------------------------
  subroutine horizCorrelFunction(NormB)
    !
    ! :Purpose: To compute and write the horizontal correlation function of
    !           every variables and levels in the correlation formulation of
    !           the B matrix (i.e. C matrix)
    !
    implicit none

    ! arguments
    real(8), intent(in)    :: NormB(bhi%nVarLev,bhi%nVarLev,0:nTrunc)

    ! locals
    real(8), allocatable   :: SpectralStateVar(:,:,:)
    real(8), allocatable   :: GridState(:,:,:)

    type(struct_gsv) :: statevector
    real(8), pointer :: ptr3d_r8(:,:,:)

    character(len=4), pointer :: varNamesList(:)

    integer   :: e, ila, p, k, totwvnb

    type(struct_lst)  :: lst_cor ! Spectral transform Parameters

    character(len=24) :: kind

    nullify(varNamesList)
    call ens_varNamesList(varNamesList,ensPerts) 
    call gsv_allocate(statevector, ens_getNumStep(ensPerts),                    &
                      hco_bhi, ens_getVco(ensPerts), varNames_opt=varNamesList, &
                      datestamp_opt=tim_getDatestamp(), mpi_local_opt=.false.,  &
                      dataKind_opt=8 )

    !
    !- 1. Setup the (LAM) spectral transform
    !
    call lst_setup(lst_cor,                             & ! OUT
                   hco_bhi%ni, hco_bhi%nj,              & ! IN
                   hco_bhi%dlon, nTrunc,                & ! IN
                   'NoMpi')

    allocate( SpectralStateVar(lst_cor%nla,lst_cor%nphase,bhi%nVarLev) )
    allocate( GridState(hco_bhi%ni, hco_bhi%nj, bhi%nVarLev) )

    !- 3.2.1 Spectral transform of a delta function (at the center of the domain)
    GridState(:,:,:) = 0.d0
    GridState(hco_bhi%ni/2,hco_bhi%nj/2,:) = 1.d0

    kind = 'GridPointToSpectral'
    call lst_VarTransform(lst_cor,            & ! IN
                          SpectralStateVar,   & ! OUT
                          GridState,          & ! IN
                          kind, bhi%nVarLev)    ! IN

    !- 3.2.2 Apply the horizontal correlation function
    !$OMP PARALLEL DO PRIVATE (totwvnb,e,ila,p,k)
    do totwvnb = 0, nTrunc
      do e = 1, lst_cor%nePerK(totwvnb)
        ila = lst_cor%ilaFromEK(e,totwvnb)
        do p = 1, lst_cor%nphase
          do k = 1, bhi%nVarLev
            SpectralStateVar(ila,p,k) = SpectralStateVar(ila,p,k) * NormB(k,k,totwvnb) * &
                 lst_cor%NormFactor(ila,p) * lst_cor%NormFactorAd(ila,p)
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    !- 3.2.3 Move back to physical space
    kind = 'SpectralToGridPoint'
    call lst_VarTransform(lst_cor,          & ! IN
                          SpectralStateVar, & ! IN
                          GridState,        & ! OUT
                          kind, bhi%nVarLev)  ! IN

    !
    !- 4.  Write to file
    !
    call gsv_getField(statevector,ptr3d_r8)
    ptr3d_r8(:,:,:) = GridState(:,:,:)
    call gio_writeToFile(statevector, './horizCorrel.fst', 'CORRELFUNC')

    deallocate(SpectralStateVar)
    deallocate(GridState)
    call gsv_deallocate(statevector)

  end subroutine horizCorrelFunction

  !--------------------------------------------------------------------------
  ! applyHorizLoc
  !--------------------------------------------------------------------------
  subroutine applyHorizLoc(NormPowerSpectrum)
    !
    ! :Purpose: To apply horizontal localization to the  spectral correlations
    !
    implicit none

    ! arguments
    real(8), intent(inout) :: NormPowerSpectrum(bhi%nVarLev,0:nTrunc)

    ! locals
    real(8), allocatable   :: SpectralStateVar(:,:,:)
    real(8), allocatable   :: GridState(:,:,:)
    real(8), allocatable   :: GridStateLoc(:,:,:)
    real(8), allocatable   :: PowerSpectrum(:,:)
    real(8), allocatable   :: SumWeight(:)
    real(8), allocatable   :: local_length(:)

    type(struct_lst)  :: lst_hloc ! Spectral transform Parameters

    integer :: totwvnb, var, k, e, ila, p

    real(8)  :: hlocalize

    character(len=24) :: kind

    write(*,*)
    write(*,*) 'applyHorizLoc: Starting...'

    !
    !- 1. Setup the (LAM) spectral transform
    !
    call lst_setup(lst_hloc,                            & ! OUT
                   hco_bhi%ni, hco_bhi%nj,              & ! IN
                   hco_bhi%dlon, nTrunc,                & ! IN
                   'NoMpi')

    !
    !- 1. Get the original gridpoint horizontal correlations
    !

    allocate( SpectralStateVar(lst_hloc%nla,lst_hloc%nphase,bhi%nVarLev) )
    allocate( GridState(hco_bhi%ni, hco_bhi%nj, bhi%nVarLev) )

    !- 1.1 Spectral transform of a delta function (at the lower-left of the domain)
    GridState(:,:,:) = 0.d0
    GridState(1,1,:) = 1.d0

    kind = 'GridPointToSpectral'
    call lst_VarTransform(lst_hloc,            & ! IN
                          SpectralStateVar,    & ! OUT
                          GridState,           & ! IN
                          kind, bhi%nVarLev)     ! IN

    !- 1.2 Apply the horizontal correlation function
    !$OMP PARALLEL DO PRIVATE (totwvnb,e,ila,p,k)
    do totwvnb = 0, nTrunc
      do e = 1, lst_hloc%nePerK(totwvnb)
        ila = lst_hloc%ilaFromEK(e,totwvnb)
        do p = 1, lst_hloc%nphase
          do k = 1, bhi%nVarLev
            SpectralStateVar(ila,p,k) = SpectralStateVar(ila,p,k) * NormPowerSpectrum(k,totwvnb) * &
                 lst_hloc%NormFactor(ila,p) * lst_hloc%NormFactorAd(ila,p)
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    !- 1.3 Move back to physical space
    kind = 'SpectralToGridPoint'
    call lst_VarTransform(lst_hloc,          & ! IN
                          SpectralStateVar,  & ! IN
                          GridState,         & ! OUT
                          kind, bhi%nVarLev)   ! IN

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
        hLocalize = hLocalize_other
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
    allocate(PowerSpectrum(bhi%nVarLev,0:nTrunc))

    !- 3.1 Transform to spectral space
    kind = 'GridPointToSpectral'
    call lst_VarTransform(lst_hloc,          & ! IN
                          SpectralStateVar,  & ! OUT
                          GridState,         & ! IN
                          kind, bhi%nVarLev)   ! IN

    !- 3.2 Compute band mean
    allocate(SumWeight(0:nTrunc))
    SumWeight(:) = 0.d0

    PowerSpectrum(:,:) = 0.d0
    do totwvnb = 0, nTrunc
      do e = 1, lst_hloc%nePerK(totwvnb)
        ila = lst_hloc%ilaFromEK(e,totwvnb)
        do p = 1, lst_hloc%nphase
          SumWeight(totwvnb) = SumWeight(totwvnb) + lst_hloc%Weight(ila)
          do k = 1, bhi%nVarLev
            PowerSpectrum(k,totwvnb) = PowerSpectrum(k,totwvnb) + &
                 lst_hloc%Weight(ila) * abs(SpectralStateVar(ila,p,k))
          end do
        end do
      end do
    end do

    do totwvnb = 0, nTrunc
      if (SumWeight(totwvnb) /= 0.d0) then
        PowerSpectrum(:,totwvnb) = PowerSpectrum(:,totwvnb) / SumWeight(totwvnb)
      else
        PowerSpectrum(:,totwvnb) = 0.d0
      endif
    end do

    deallocate(SumWeight)

    !- 3.3 Normalize
    call normalizePowerSpectrum(PowerSpectrum,     & ! IN
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
    !
    ! :Purpose: To apply vertical localization to the spectral correlations
    !
    implicit none

    ! arguments
    real(8), intent(inout) :: SpVertCorrel(bhi%nVarLev,bhi%nVarLev,0:nTrunc)

    ! locals
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
          vLocalize1 = vLocalize_other
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
          vLocalize2 = vLocalize_other
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
            do totwvnb = 0, nTrunc
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
  ! scaleStdDev
  !--------------------------------------------------------------------------
  subroutine scaleStdDev(statevector_stdDev)
    !
    ! :Purpose: To scale the gridpoint background-error standard deviations
    !
    implicit none
    
    ! arguments
    type(struct_gsv), intent(inout) :: statevector_stdDev

    ! locals
    real(8), pointer :: ptr3d_r8(:,:,:)
    real(8) :: multFactor
    integer :: nVarLev, varLevIndex, levIndex
    character(len=4) :: varName

    write(*,*)
    write(*,*) 'scaleStdDev: Starting...'
    
    nVarLev = gsv_getNumK(statevector_stdDev)
 
    call gsv_getField(statevector_stdDev,ptr3d_r8)

    do varLevIndex = 1, nVarLev
      varName = gsv_getVarNameFromK(statevector_stdDev,varLevIndex)
      levIndex = gsv_getLevFromK(statevector_stdDev,varLevIndex)

      if ( vnl_varLevelFromVarname(varName) == 'MM' ) then
        multFactor = scaleFactor_M(levIndex)
      else if ( vnl_varLevelFromVarname(varName) == 'TH' ) then
        multFactor = scaleFactor_T(levIndex)
      else ! SF
        multFactor = scaleFactor_SF
      end if

      ptr3d_r8(:,:,varLevIndex) = multFactor * ptr3d_r8(:,:,varLevIndex)
    end do

    write(*,*)
    write(*,*) 'scaleStdDev: Done...'

  end subroutine scaleStdDev

  !--------------------------------------------------------------------------
  ! writeVarStats
  !--------------------------------------------------------------------------
  subroutine writeVarStats(Bsqrt,statevector_stdDev)
    !
    ! :Purpose: To write data needed for VAR applications, i.e C^1/2 and stdDev
    !
    implicit none

    ! arguments
    real(8), intent(in) :: Bsqrt(bhi%nVarLev,bhi%nVarLev,0:nTrunc)
    type(struct_gsv)    :: statevector_stdDev

    ! locals
    integer   :: ier, fstouv, fnom, fstfrm, fclos
    integer   :: iunstats

    character(len=24) :: fileName = './bgcov.fst'

    write(*,*)
    write(*,*) 'Writing covariance statistics for VAR'

    !
    !- 1. Add the gridded std dev (and the Tic-Tac-Toc)
    !
    call gio_writeToFile(statevector_stdDev, trim(fileName), 'STDDEV', &
                         typvar_opt = 'E', numBits_opt = 32)

    !
    !- 2. Add C^1/2
    !
    if (mmpi_myid == 0) then
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
    end if

  end subroutine writeVarStats

  !--------------------------------------------------------------------------
  ! writeDiagStats
  !--------------------------------------------------------------------------
  subroutine writeDiagStats(NormB,SpVertCorrel,TotVertCorrel,statevector_mean, &
                            statevector_stdDevGridPoint,PowerSpectrum,HorizScale)
    !
    ! :Purpose: To write other relevant data computed during the
    !           calculation of Bhi that are not needed for VAR applications
    !
    implicit none

    ! arguments
    real(8), intent(in) :: NormB(bhi%nVarLev,bhi%nVarLev,0:nTrunc)
    real(8), intent(in) :: SpVertCorrel(bhi%nVarLev,bhi%nVarLev,0:nTrunc)
    real(8), intent(in) :: TotVertCorrel(bhi%nVarLev,bhi%nVarLev)
    type(struct_gsv)    :: statevector_mean
    type(struct_gsv)    :: statevector_stdDevGridPoint
    real(8), intent(in) :: PowerSpectrum(bhi%nVarLev,0:nTrunc)
    real(8), intent(in) :: HorizScale(bhi%nVarLev)

    ! locals
    integer   :: ier, fstouv, fnom, fstfrm, fclos
    integer   :: iunstats

    character(len=24) :: fileName = './bgcov_diag.fst'

    write(*,*)
    write(*,*) 'Writing Diagnostics'

    !
    !- 1. Add the gridded mean and std dev (and the Tic-Tac-Toc)
    !
    call gio_writeToFile(statevector_mean, trim(fileName), 'ENSMEAN', &
                         typvar_opt = 'E', numBits_opt = 32)
    call gio_writeToFile(statevector_stdDevGridPoint, trim(fileName), 'STDDEV_GRIDP', &
                         typvar_opt = 'E', numBits_opt = 32)

    !
    !- 2. Add stats in spectral space
    !
    if (mmpi_myid == 0) then
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
    end if

  end subroutine writeDiagStats

  !--------------------------------------------------------------------------
  ! writeSpVertCorrel
  !--------------------------------------------------------------------------
  subroutine writeSpVertCorrel(SpVertCorrel,iun,nomvar_in,etiket_in)
    !
    ! :Purpose: To write vertical correlations in spectral space
    !
    implicit none

    ! arguments
    real(8), intent(in) :: SpVertCorrel(bhi%nVarLev,bhi%nVarLev,0:nTrunc)
    integer, intent(in) :: iun
    character(len=*), intent(in) :: nomvar_in
    character(len=*), intent(in) :: etiket_in

    ! locals
    real(4), allocatable :: work2d(:,:)

    real(4) :: work

    integer   :: ier, fstecr, totwvnb

    integer :: dateo, npak, ni, nj, nk
    integer :: ip1, ip2, ip3, deet, npas, datyp
    integer :: ig1 ,ig2 ,ig3 ,ig4

    character(len=1)  :: grtyp
    character(len=4)  :: nomvar
    character(len=2)  :: typvar
    character(len=12) :: etiket

    allocate(work2d(bhi%nVarLev, bhi%nVarLev))

    !- Loop over Total Wavenumbers
    do totwvnb = 0, nTrunc

      npak   = -32
      dateo  = 0
      deet   = 0
      npas   = 0
      ni     = bhi%nVarLev
      nj     = bhi%nVarLev
      nk     = 1
      ip1    = 0
      ip2    = totwvnb
      ip3    = nEns
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
    !
    ! :Purpose: To write the total vertical correlations
    !           (i.e. gridpoint equivalent)
    !
    implicit none

    ! arguments
    real(8), intent(in) :: TotVertCorrel(bhi%nVarLev,bhi%nVarLev)
    integer, intent(in) :: iun
    character(len=*), intent(in) :: nomvar_in
    character(len=*), intent(in) :: etiket_in

    ! locals
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
    ip3    = nEns
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
    !
    ! :Purpose: To write the power spectrum 
    !
    implicit none

    ! arguments
    real(8), intent(in) :: PowerSpectrum(bhi%nVarLev,0:nTrunc)
    integer, intent(in) :: iun
    integer, intent(in) :: cv_type
    character(len=*), intent(in) :: Etiket_in

    ! locals
    real(4), allocatable :: workecr(:,:)

    real(4)   :: work

    integer   :: ier, fstecr
    integer   :: var, k, kgdim

    integer :: dateo, npak, ni, nj, nk
    integer :: ip1, ip2, ip3, deet, npas, datyp
    integer :: ig1 ,ig2 ,ig3 ,ig4

    character(len=1)  :: grtyp
    character(len=4)  :: nomvar
    character(len=2)  :: typvar
    character(len=12) :: etiket

    allocate(workecr(nTrunc+1, 1))

    !- Loop over Control Variables
    do var = 1, bhi%nControlVariable

      !- Loop over vertical Levels
      do k = 1, bhi%controlVariable(var)%nlev

        npak   = -32
        dateo  = 0
        deet   = 0
        npas   = 0
        ni     = nTrunc + 1
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
    !
    ! :Purpose: To write the horizontal lenght scales
    !
    implicit none

    real(8), intent(in) :: HorizScale(bhi%nVarLev)
    integer, intent(in) :: iun
    integer, intent(in) :: cv_type
    character(len=*), intent(in) :: Etiket_in

    real(4), allocatable :: workecr(:,:,:)

    real(4)   :: work

    integer   :: ier, fstecr, var

    integer :: dateo, npak, ni, nj, nk
    integer :: ip1, ip2, ip3, deet, npas, datyp
    integer :: ig1 ,ig2 ,ig3 ,ig4

    character(len=1)  :: grtyp
    character(len=4)  :: nomvar
    character(len=2)  :: typvar
    character(len=12) :: etiket

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
    !
    ! :Purpose: To write the control variable related info
    !
    implicit none

    ! arguments
    integer, intent(in) :: iun

    ! locals
    integer :: ier, fstecr, fstecr_s

    real(8) :: work

    integer :: npak, var, dateo, ni, nj
    integer :: ip1,ip2,ip3,deet,npas,datyp,ig1,ig2,ig3,ig4

    character(len=1)  :: grtyp
    character(len=2)  :: typvar
    character(len=4)  :: nomvar

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
    !
    ! :Purpose: To write the MM and TH pressure profiles used for vertical localization
    !
    implicit none

    ! locals
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
  subroutine calcLocalCorrelations(ensPerts)
    !
    ! :Purpose: To compute the local horizontal correlation for some 'reference' grid points
    !
    implicit none

    ! arguments
    type(struct_ens) :: ensPerts

    ! locals
    type(struct_gsv) :: statevector_locHorizCor
    type(struct_gsv) :: statevector_oneMember
    type(struct_gsv) :: statevector_oneMemberTiles

    real(8), pointer :: ptr3d_r8(:,:,:)
    real(8), pointer :: ptr3d_r8_oneMember(:,:,:)

    real(8) :: dnEns

    integer :: i, j, k, ens
    integer :: blocklength_x, blocklength_y, blockpadding, nirefpoint, njrefpoint
    integer :: iref_id, jref_id, iref, jref
    integer :: imin, imax, jmin, jmax

    character(len=4), pointer :: varNamesList(:)

    integer :: ier, fclos, fnom, nulnam

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
    ier = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=NAMHVCORREL_LOCAL)
    write(*,nml=NAMHVCORREL_LOCAL)
    ier = fclos(nulnam)

    blocklength_x = hco_ens%ni / nirefpoint ! Horizontal correlation will be compute blocklength x blocklength gridpoint
    ! around each reference point
    blocklength_y = hco_ens%nj / njrefpoint ! Horizontal correlation will be compute blocklength x blocklength gridpoint
    ! around each reference point

    nullify(varNamesList)
    call ens_varNamesList(varNamesList,ensPerts) 

    call gsv_allocate(statevector_locHorizCor, ens_getNumStep(ensPerts),                     &
                      ens_getHco(ensPerts), ens_getVco(ensPerts), varNames_opt=varNamesList, &
                      datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true.,                &
                      mpi_distribution_opt='VarsLevs', dataKind_opt=8 )

    call gsv_allocate(statevector_oneMemberTiles, ens_getNumStep(ensPerts),                  &
                      ens_getHco(ensPerts), ens_getVco(ensPerts), varNames_opt=varNamesList, &
                      datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true.,                &
                      mpi_distribution_opt='Tiles', dataKind_opt=8 )

    call gsv_allocate(statevector_oneMember, ens_getNumStep(ensPerts),                       &
                      ens_getHco(ensPerts), ens_getVco(ensPerts), varNames_opt=varNamesList, &
                      datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true.,                &
                      mpi_distribution_opt='VarsLevs', dataKind_opt=8 )

    call gsv_zero(statevector_locHorizCor)

    dnEns = 1.0d0/dble(nEns-1)

    call gsv_getField(statevector_locHorizCor,ptr3d_r8)

    do ens = 1, nEns
      call ens_copyMember(ensPerts, statevector_oneMemberTiles, ens)
      call gsv_transposeTilesToVarsLevs(statevector_oneMemberTiles, statevector_oneMember)
      call gsv_getField(statevector_oneMember,ptr3d_r8_oneMember)

      do k = statevector_locHorizCor%mykBeg, statevector_locHorizCor%mykEnd
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
                ptr3d_r8(i,j,k) = ptr3d_r8(i,j,k) + &
                     ptr3d_r8_oneMember(i,j,k)*ptr3d_r8_oneMember(iref,jref,k)
              end do
            end do
          end do
        end do
      end do

    end do

    call gsv_scale(statevector_locHorizCor,dnEns)

    write(*,*) 'finished computing the local horizontal correlations...'

    !
    !- 4.  Write to file
    !
    call gio_writeToFile(statevector_locHorizCor, './horizCorrelLocal.fst', 'HCORREL_LOC', &
                         typvar_opt = 'E', numBits_opt = 32)

    call gsv_deallocate(statevector_locHorizCor)
    call gsv_deallocate(statevector_oneMember)
    call gsv_deallocate(statevector_oneMemberTiles)

  end subroutine calcLocalCorrelations

end module calcStatsLam_mod
