
module bMatrix1DVar_mod
  ! MODULE bMatrix1DVar_mod (prefix='bmat1D' category='2. B and R matrices')
  !
  !:Purpose: contains all 1Dvar B matrices.
  !
  use mathPhysConstants_mod
  use columnData_mod
  use columnVariableTransforms_mod
  use controlVector_mod
  use gridStatevector_mod
  use gridstatevectorFileIO_mod
  use horizontalCoord_mod
  use midasMpi_mod 
  use obsSpaceData_mod
  use timeCoord_mod
  use utilities_mod
  use verticalCoord_mod
  use tovsNL_mod
  use var1D_mod
  use filenames_mod
  use localizationFunction_mod
  use varNameList_mod
  use ensembleStateVector_mod
  use stateToColumn_mod
  use calcHeightAndPressure_mod
  implicit none
  save
  private

  ! public procedures
  public :: bmat1D_bsetup
  public :: bmat1D_sqrtB, bmat1D_sqrtBT
  public :: bmat1D_finalize, bmat1D_get1DVarIncrement

  ! public variables
  public :: bmat1D_includeAnlVar, bmat1D_numIncludeAnlVar

  integer                       :: bmat1D_numIncludeAnlVar
  character(len=4), allocatable :: bmat1D_includeAnlVar(:)

  type(struct_hco), pointer :: hco_yGrid
  logical             :: initialized = .false.
  integer             :: nkgdim
  integer             :: cvDim_mpilocal
 
  real(8), allocatable :: bSqrtLand(:,:,:), bSqrtSea(:,:,:)
  real(8), allocatable :: bSqrtEns(:,:,:)
  real(4), allocatable :: latLand(:), lonLand(:), latSea(:), lonSea(:)
  integer              :: nLonLatPosLand, nLonLatPosSea
  integer, external    :: get_max_rss
  integer,          parameter :: numMasterBmat = 2
  character(len=4), parameter :: masterBmatTypeList (numMasterBmat) = (/ 'HI', 'ENS' /)
  character(len=8), parameter :: masterBmatLabelList(numMasterBmat) = (/'B_HI', 'B_ENS' /)
  logical,          parameter :: masterbmatIs3dList (numMasterBmat) = (/.true., .true. /) 
  integer            :: numBmat
  integer, parameter :: numBmatMax = 10
  character(len=4) :: bmatTypeList  (numBmatMax)
  character(len=9) :: bmatLabelList (numBmatMax)
  logical          :: bmatIs3dList  (numBmatMax)
  logical          :: bmatActive    (numBmatMax)
  type(struct_columnData), allocatable :: ensColumns(:)
  type(struct_columnData) :: meanColumn
  type(struct_ens) :: ensembles

  ! Namelist variables
  integer          :: nEns                             ! ensemble size
  real(8)          :: vlocalize                        ! vertical localization length scale
  character(len=4) :: includeAnlVar(vnl_numvarmax)     ! list of variable names to include in B matrix
  integer :: numIncludeAnlVar                          ! MUST NOT BE INCLUDED IN NAMELIST!
  real(8) :: scaleFactorHI(vco_maxNumLevels)           ! scaling factors for HI variances
  real(8) :: scaleFactorHIHumidity(vco_maxNumLevels)   ! scaling factors for HI humidity variances
  real(8) :: scaleFactorHITG                     ! scaling factors for HI skin temperature variances
  real(8) :: scaleFactorEns(vco_maxNumLevels)          ! scaling factors for Ens variances
  real(8) :: scaleFactorEnsHumidity(vco_maxNumLevels)  ! scaling factors for Ens humidity variances
  real(8) :: scaleFactorEnsTG                    ! scaling factors for Ens skin temperature variances
  real(8) :: scaleFactorEnsTGCorrelation                  ! scaling factors for corrleation between Ens skin temperature error and other variable error 
  logical :: dumpBmatrixTofile                         ! flag to control output of B matrices to Bmatrix.bin binary file
  logical :: doAveraging                               ! flag to control output the average instead of the invidual B matrices
  real(8) :: latMin                                    ! minimum latitude of the Bmatrix latitude-longitude output box
  real(8) :: latMax                                    ! maximum latitude of the Bmatrix latitude-longitude output box
  real(8) :: lonMin                                    ! minimum longitude of the Bmatrix latitude-longitude output box
  real(8) :: lonMax                                    ! maximum longitude of the Bmatrix latitude-longitude output box
  NAMELIST /NAMBMAT1D/ scaleFactorHI, scaleFactorHIHumidity, scaleFactorHITG, &
      scaleFactorENs, scaleFactorEnsHumidity, scaleFactorEnsTG, scaleFactorEnsTGCorrelation, &
      nEns, vLocalize, includeAnlVar, numIncludeAnlVar, &
      dumpBmatrixTofile, latMin, latMax, lonMin, lonMax, doAveraging

contains

  !--------------------------------------------------------------------------
  ! bmat1D_bsetup
  !--------------------------------------------------------------------------
  subroutine bmat1D_bsetup(vco_in, hco_in, obsSpaceData)
    !
    !:Purpose: To initialize the 1Dvar analysis Background term.
    !
    implicit none

    ! Arguments:
    type(struct_vco), pointer, intent(in)    :: vco_in
    type(struct_hco), pointer, intent(in)    :: hco_in
    type (struct_obs),         intent(inout) :: obsSpaceData

    ! Locals:
    integer :: cvdim
    integer :: masterBmatIndex, bmatIndex
    logical :: active
    integer :: nulnam, ierr
    integer, external ::  fnom, fclos
    integer :: varIndex

    call utl_tmg_start(50, '--Bmatrix')

    ! default values for namelist variables
    scaleFactorHI(:) = 0.d0
    scaleFactorHIHumidity(:) = 1.d0
    scaleFactorHITG = 1.d0
    scaleFactorEns(:) = 0.d0
    scaleFactorEnsHumidity(:) = 1.d0
    scaleFactorEnsTG = 1.d0
    scaleFactorEnsTGCorrelation = 1.d0
    
    nEns = -1
    vLocalize = -1.d0
    includeAnlVar(:)= ''
    numIncludeAnlVar = MPC_missingValue_INT
    dumpBmatrixTofile = .false.
    doAveraging = .false.
    latMin = -1000.d0
    latMax = 1000.d0
    lonMin = -1000.d0
    lonMax = 1000.d0
    
    nulnam = 0
    ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    read(nulnam, nml=nambmat1D, iostat=ierr)
    if ( ierr /= 0 ) call utl_abort( 'bmat1D_bsetup: Error reading namelist' )
    if ( mmpi_myid == 0 ) write( *, nml = nambmat1D )
    ierr = fclos( nulnam )
    if (numIncludeAnlVar /= MPC_missingValue_INT) then
      call utl_abort('bmat1D_bsetup: check NAMBMAT1D namelist section: numIncludeAnlVar should be removed')
    end if
    numIncludeAnlVar = 0
    do varIndex = 1, vnl_numvarmax
      if (trim(includeAnlVar(varIndex)) == '') exit
      numIncludeAnlVar = numIncludeAnlVar + 1
    end do
    bmat1D_numIncludeAnlVar = numIncludeAnlVar
    allocate( bmat1D_includeAnlVar(bmat1D_numIncludeAnlVar) )
    bmat1D_includeAnlVar(1:bmat1D_numIncludeAnlVar) = includeAnlVar(1:numIncludeAnlVar)

    !
    !- 1.  Setup the B matrices
    !
    do masterBmatIndex = 1, numMasterBmat

      select case( trim(masterBmatTypeList(masterBmatIndex)) )
      case ('HI')
        !- 1.1 Time-Mean Homogeneous and Isotropic...
        write(*,*) 'bmat1D_bsetup: Setting up the modular GLOBAL HI 1D covariances...'
        call utl_tmg_start(51, '----B_HI_Setup')
        call bmat1D_SetupBHi(vco_in, obsSpaceData, cvdim)
        call utl_tmg_stop(51)
        write(*,*) ' bmat1D_bsetup: cvdim= ', cvdim
      case ('ENS')
        !- 1.2 ensemble based
        write(*,*) 'bmat1D_bsetup: Setting up the ensemble based 1D matrix.'
        call utl_tmg_start(54, '----B_ENS_Setup')
        call bmat1D_SetupBEns(vco_in, hco_in, obsSpaceData, cvdim)
        call utl_tmg_stop(54)
        write(*,*) ' bmat1D_bsetup: cvdim= ', cvdim
      case default
        call utl_abort('bmat1D_bSetup: requested bmatrix type does not exist ' // trim(masterBmatTypeList(masterBmatIndex)))
      end select

      !- 1.2 Append the info to the B matrix info arrays and setup the proper control sub-vectors
      numBmat = numBmat + 1
      bmatLabelList (numBmat) = trim(masterbmatLabelList(masterBmatIndex))
      bmatTypeList  (numBmat) = masterBmatTypeList(masterBmatIndex)
      bmatIs3dList  (numBmat) = masterbmatIs3dList(masterBmatIndex)
      call cvm_setupSubVector(bmatLabelList(numBmat), bmatTypeList(numBmat), cvdim)
    end do

    !
    !- 2. Print a summary and set the active B matrices array
    !
    write(*,*)
    write(*,*) 'bmat1D_bsetup SUMMARY, number of B matrices found = ', numBmat
    do bmatIndex = 1, numBmat
      write(*,*) '  B matrix #', bmatIndex
      active = cvm_subVectorExists(bmatLabelList(bmatIndex))
      if (active) then
        write(*,*) '   ACTIVE'
      else
        write(*,*) '   NOT USED'
      end if
      write(*,*) '     -> label       = ', bmatLabelList (bmatIndex)
      write(*,*) '     -> type        = ', bmatTypeList  (bmatIndex)
      if (active) then
        write(*,*) '     -> is 3D       = ', bmatIs3dList  (bmatIndex)
      end if
      bmatActive(bmatIndex) = active
    end do

    call utl_tmg_stop(50)

  end subroutine bmat1D_bsetup

  !--------------------------------------------------------------------------
  !  bmat1D_setupBHi
  !--------------------------------------------------------------------------
  subroutine bmat1D_setupBHi(vco_in, obsSpaceData, cvDim_out)
    !
    ! :Purpose: to setup bmat1D module
    !
    implicit none

    ! Arguments:
    type(struct_vco), pointer, intent(in)  :: vco_in
    type (struct_obs)        , intent(in)  :: obsSpaceData
    integer                  , intent(out) :: cvDim_out

    ! Locals:
    integer :: levelIndex, ierr
    integer, external ::  fnom, fclos
    integer :: Vcode_anl
    logical :: fileExists
    integer :: nulbgst=0
    type(struct_vco), pointer :: vco_file => null()
    type(struct_vco), target  :: vco_1Dvar
    type(struct_vco), pointer :: vco_anl
    character(len=18) :: oneDBmatLand = './Bmatrix_land.bin'
    character(len=17) :: oneDBmatSea  = './Bmatrix_sea.bin'
    character(len=4), allocatable :: includeAnlVarHi(:)
    integer :: extractDate, locationIndex, varIndex, numIncludeAnlVarHi
    integer :: shiftLevel, varLevIndex1, varLevIndex2
    integer :: varLevIndexBmat
    logical, save :: firstCall=.true.
    real(8), allocatable :: bMatrix(:,:), multFactor(:)

    if (.not. (gsv_varExist(varName='TT') .and.  &
               gsv_varExist(varName='UU') .and.  &
               gsv_varExist(varName='VV') .and.  &
               (gsv_varExist(varName='HU').or.gsv_varExist(varName='LQ')) .and.  &
               gsv_varExist(varName='P0')) ) then
      call utl_abort('bmat1D_setupBHi: Some or all weather fields are missing. If it is desired to deactivate&
           & the weather assimilation, then all entries of the array SCALEFACTORHI in the namelist NAMVAR1D&
           & should be set to zero.')
    end if

    if (.not. gsv_varExist(varName='TG')) then
      write(*,*) 'bmat1D_setupBHi: WARNING: The TG field is missing. This must be present when assimilating'
      write(*,*) 'radiance observations.'
    end if

    if (firstCall) then
      call var1D_setup(obsSpaceData)
      firstCall = .false.
    end if

    if (mmpi_myid == 0) write(*,*) 'bmat1D_setupBHi: Starting'
    if (mmpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

    do levelIndex = 1, vco_maxNumLevels
      if( scaleFactorHI(levelIndex) > 0.0d0 ) then 
        scaleFactorHI(levelIndex) = sqrt( scaleFactorHI(levelIndex))
      else
        scaleFactorHI(levelIndex) = 0.0d0
      end if
    end do
   
    do levelIndex = 1, vco_maxNumLevels
      if(scaleFactorHIHumidity(levelIndex) > 0.0d0) then 
        scaleFactorHIHumidity(levelIndex) = sqrt(scaleFactorHIHumidity(levelIndex))
      else
        scaleFactorHIHumidity(levelIndex) = 0.0d0
      end if
    end do

    if(scaleFactorHITG > 0.0d0) then 
      scaleFactorHITG = sqrt(scaleFactorHITG)
    else
      scaleFactorHITG = 0.0d0
    end if

    if ( sum(scaleFactorHI(1:vco_maxNumLevels)) == 0.0d0 ) then
      if ( mmpi_myid == 0 ) write(*,*) 'bmat1D_setupBHi: scaleFactorHI=0, skipping rest of setup'
      cvDim_out = 0
      return
    end if

    if (mmpi_myid == 0) write(*,*) 'bmat1D_setupBHi: Read 1DVar background statistics'
    inquire(file=trim(oneDBmatLand), exist=fileExists)
    if ( fileExists ) then
      ierr = fnom(nulbgst, trim(oneDBmatLand), 'FTN+SEQ+UNF+OLD+R/O', 0)
    else
      call utl_abort('bmat1D_setupBHi: No 1DVar BACKGROUND STAT FILE ' // trim(oneDBmatLand))
    end if

    read(nulbgst) extractDate, vco_1Dvar%nLev_T, vco_1Dvar%nLev_M, vco_1Dvar%Vcode, &
         vco_1Dvar%ip1_sfc, vco_1Dvar%ip1_T_2m, vco_1Dvar%ip1_M_10m, numIncludeAnlVarHi, nkgdim, nLonLatPosLand
    allocate( vco_1Dvar%ip1_T(vco_1Dvar%nLev_T), vco_1Dvar%ip1_M(vco_1Dvar%nLev_M) )
    if (numIncludeAnlVarHi /= bmat1D_numIncludeAnlVar) then
      write(*,*) 'numIncludeAnlVarHi, bmat1D_numIncludeAnlVar= ', numIncludeAnlVarHi, bmat1D_numIncludeAnlVar
      call utl_abort('bmat1D_setupBHi: incompatible number of 1DVar analyzed variables in ' // trim(oneDBmatLand))
    end if

    allocate (multFactor(nkgdim))
    if(vco_1Dvar%Vcode == 5002) then
      shiftLevel = 1
    else
      shiftLevel = 0
    end if
    varLevIndexBmat = 0
    do varIndex = 1, bmat1D_numIncludeAnlVar
      select case(bmat1D_includeAnlVar(varIndex))
        case('TT')
          do levelIndex = 1, vco_1Dvar%nLev_T
            varLevIndexBmat = varLevIndexBmat + 1
            multFactor(varLevIndexBmat) = scaleFactorHI(levelIndex)
          end do
        case('HU')
          do levelIndex = 1, vco_1Dvar%nLev_T
            varLevIndexBmat = varLevIndexBmat + 1
            multFactor(varLevIndexBmat)= scaleFactorHIHumidity(levelIndex) * scaleFactorHI(levelIndex)
          end do
        case('UU','VV')
          do levelIndex = 1, vco_1Dvar%nLev_M
            varLevIndexBmat = varLevIndexBmat + 1
            multFactor(varLevIndexBmat) = scaleFactorHI(levelIndex+shiftLevel)
          end do
        case('TG')
          varLevIndexBmat = varLevIndexBmat + 1
          multFactor(varLevIndexBmat) = scaleFactorHI(max(vco_1Dvar%nLev_T,vco_1Dvar%nLev_M))* scaleFactorHITG
        case('P0')
          varLevIndexBmat = varLevIndexBmat + 1
          multFactor(varLevIndexBmat) = scaleFactorHI(max(vco_1Dvar%nLev_T,vco_1Dvar%nLev_M))
        case default
          call utl_abort('bmat1D_setupBHi: unsupported variable ' // bmat1D_includeAnlVar(varIndex))
        end select
    end do

    allocate( includeAnlVarHi(bmat1D_numIncludeAnlVar) )
    allocate( bMatrix(nkgdim,nkgdim) )
    allocate( latLand(nLonLatPosLand), lonLand(nLonLatPosLand))       
    allocate( bSqrtLand(nLonLatPosLand, nkgdim, nkgdim) )
    read(nulbgst) vco_1Dvar%ip1_T(:), vco_1Dvar%ip1_M(:), includeAnlVarHi(:)
    if (any(includeAnlVarHi /= bmat1D_includeAnlVar)) then
      do varIndex = 1, bmat1D_numIncludeAnlVar
        write(*,*) varIndex, includeAnlVarHi(varIndex), bmat1D_includeAnlVar(varIndex)
      end do
      call utl_abort('bmat1D_setupBHi: incompatible 1DVar analyzed variable list in ' // trim(oneDBmatLand))
    end if
    deallocate(includeAnlVarHi)

    do locationIndex = 1, nLonLatPosLand
      read(nulbgst) latLand(locationIndex), lonLand(locationIndex), bMatrix(:,:)
      !application of the scaling factor
      do varLevIndex2 = 1, nkgdim
        do varLevIndex1 = 1, nkgdim
          bSqrtLand(locationIndex, varLevIndex1, varLevIndex2) = &
               bMatrix(varLevIndex1, varLevIndex2) *  multFactor(varLevIndex1) * multFactor(varLevIndex2)
        end do
      end do
      call utl_matsqrt(bSqrtLand(locationIndex, :, :), nkgdim, 1.d0, printInformation_opt=.false. )
    end do
    ierr = fclos(nulbgst)

    inquire(file=trim(oneDBmatSea), exist=fileExists)
    if ( fileExists ) then
      ierr = fnom(nulbgst, trim(oneDBmatSea), 'FTN+SEQ+UNF+OLD+R/O', 0)
    else
      call utl_abort('bmat1D_setupBHi: No 1DVar BACKGROUND STAT FILE ' // trim(oneDBmatSea))
    end if
    read(nulbgst) extractDate, vco_1Dvar%nLev_T, vco_1Dvar%nLev_M, vco_1Dvar%Vcode, &
         vco_1Dvar%ip1_sfc, vco_1Dvar%ip1_T_2m, vco_1Dvar%ip1_M_10m, numIncludeAnlVarHi, nkgdim, nLonLatPosSea
    if (numIncludeAnlVarHi /= bmat1D_numIncludeAnlVar) then
      write(*,*) 'numIncludeAnlVarHi, bmat1D_numIncludeAnlVar= ', numIncludeAnlVarHi, bmat1D_numIncludeAnlVar
      call utl_abort('bmat1D_setupBHi: incompatible number of 1DVar analyzed variables in ' // trim(oneDBmatSea))
    end if
    allocate( bSqrtSea(nLonLatPosSea, nkgdim, nkgdim) )
    allocate( latSea(nLonLatPosSea), lonSea(nLonLatPosSea))
    allocate( includeAnlVarHi(bmat1D_numIncludeAnlVar) )
    read(nulbgst) vco_1Dvar%ip1_T(:), vco_1Dvar%ip1_M(:), includeAnlVarHi(:)
    if (any(includeAnlVarHi /= bmat1D_includeAnlVar)) then
      do varIndex = 1, bmat1D_numIncludeAnlVar
        write(*,*) varIndex, includeAnlVarHi(varIndex), bmat1D_includeAnlVar(varIndex)
      end do
      call utl_abort('bmat1D_setupBHi: incompatible 1DVar analyzed variable list in ' // trim(oneDBmatSea))
    end if
    deallocate(includeAnlVarHi)
    do locationIndex = 1, nLonLatPosSea
      read(nulbgst) latSea(locationIndex), lonSea(locationIndex), bMatrix(:,:)
      !application of the scaling factor
      do varLevIndex2 = 1, nkgdim
        do varLevIndex1 = 1, nkgdim
          bSqrtSea(locationIndex, varLevIndex1, varLevIndex2) = &
               bMatrix(varLevIndex1, varLevIndex2) *  multFactor(varLevIndex1) * multFactor(varLevIndex2)
        end do
      end do
      call utl_matsqrt(bSqrtSea(locationIndex,:,:), nkgdim, 1.d0, printInformation_opt=.false. )
    end do
    ierr = fclos(nulbgst)
  
    deallocate(bMatrix, multFactor)

    vco_1Dvar%initialized = .true.
    vco_1Dvar%vGridPresent = .false.
    vco_file => vco_1Dvar
    vco_anl => vco_in
    if (.not. vco_equal(vco_anl,vco_file)) then
      call utl_abort('bmat1D_setupBHi: vco from analysisgrid and cov file do not match')
    end if
    if (mmpi_myid == 0) write(*,*) 'bmat1D_setupBHi: nLev_M, nLev_T=', vco_1Dvar%nLev_M, vco_1Dvar%nLev_T
    Vcode_anl = vco_anl%vCode
    if(Vcode_anl /= 5002 .and. Vcode_anl /= 5005) then
      write(*,*) 'Vcode_anl = ',Vcode_anl
      call utl_abort('bmat1D_setupBHi: unknown vertical coordinate type!')
    end if
    
    cvDim_out = nkgdim * var1D_validHeaderCount
    cvDim_mpilocal = cvDim_out
    initialized = .true.

    if (mmpi_myid == 0) write(*,*) 'bmat1D_setupBHi: Exiting'
    if (mmpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  end subroutine bmat1D_setupBHi

  !--------------------------------------------------------------------------
  !  bmat1D_setupBEns
  !--------------------------------------------------------------------------
  subroutine bmat1D_setupBEns(vco_in, hco_in, obsSpaceData, cvDim_out)
    !
    ! :Purpose: to setup bmat1D module
    !
    implicit none

    ! Arguments
    type(struct_vco), pointer, intent(in)    :: vco_in
    type(struct_hco), pointer, intent(in)    :: hco_in
    type (struct_obs)        , intent(inout) :: obsSpaceData
    integer                  , intent(out)   :: cvDim_out

    ! Locals:
    character(len=256) :: ensPathName = 'ensemble'
    character(len=256) :: ensFileName
    type(struct_vco), pointer :: vco_file => null()
    type(struct_vco), pointer :: vco_ens => null()
    type(struct_hco), pointer :: hco_ens => null()
    type(struct_gsv)          :: stateVector, stateVectorMean
    integer, allocatable :: dateStampList(:)
    character(len=12) :: hInterpolationDegree='LINEAR' ! select degree of horizontal interpolation (if needed)
    integer :: memberIndex, columnIndex, headerIndex, varIndex, levIndex
    integer :: varLevIndex1, varLevIndex2
    integer :: varLevIndexBmat, varLevIndexCol
    integer :: numStep, levIndexColumn
    real(8), allocatable :: scaleFactor_M(:), scaleFactor_T(:)
    real(8) :: scaleFactor_SF, zr
    logical :: useAnlLevelsOnly, EnsTopMatchesAnlTop
    real(8), pointer :: pressureProfileFile_M(:), pressureProfileInc_M(:)
    real(8) :: pSurfRef
    integer :: nj, latPerPE, latPerPEmax, myLatBeg, myLatEnd
    integer :: ni, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd
    integer :: nLevEns_M, nLevEns_T
    integer :: nLevInc_M, nLevInc_T
    real(8) :: logP1, logP2
    real(8), pointer :: currentProfile(:), meanProfile(:)
    real(8), allocatable :: lineVector(:,:), meanPressureProfile(:), multFactor(:)
    integer, allocatable :: varLevColFromVarLevBmat(:)
    character(len=4), allocatable :: varNameFromVarLevIndexBmat(:)
    character(len=2) :: varLevel

    if (mmpi_myid == 0) write(*,*) 'bmat1D_setupBEns: Starting'
    if (mmpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

    if (nEns <= 0) then
      if (mmpi_myid == 0) write(*,*) 'bmat1D_setupBEns: no Ensemble members, skipping rest of setup'
      cvdim_out = 0
      return
    end if
    
    !- 1.1 Number of time step bins
    numStep = tim_nstepobsinc
    if (numStep /= 1 .and. numStep /= 3.and. numStep /= 5 .and. numStep /= 7) then
      call utl_abort('bmat1D_setupBEns: Invalid value for numStep (choose 1 or 3 or 5 or 7)!')
    end if
    allocate(dateStampList(numStep))
    call tim_getstamplist(dateStampList,numStep,tim_getDatestamp())
    
    hco_ens => hco_in

    !- 1.2 Horizontal grid
    ni = hco_ens%ni
    nj = hco_ens%nj
    if (hco_ens%global) then
      if (mmpi_myid == 0) write(*,*)
      if (mmpi_myid == 0) write(*,*) 'bmat1D_setupBEns: GLOBAL mode activated'
    else
      if (mmpi_myid == 0) write(*,*)
      if (mmpi_myid == 0) write(*,*) 'bmat1D_setupBEns: LAM mode activated'
    end if

    !- 1.3 Vertical levels
    call fln_ensfileName(ensFileName, ensPathName, memberIndex_opt=1)
    call vco_SetupFromFile(vco_file, ensFileName)

    !- Do we need to read all the vertical levels from the ensemble?
    useAnlLevelsOnly = vco_subsetOrNot(vco_in, vco_file)
    if (useAnlLevelsOnly) then
      write(*,*)
      write(*,*) 'bmat1D_setupBEns: only the analysis levels will be read in the ensemble '
      vco_ens  => vco_in ! the ensemble target grid is the analysis grid
      call vco_deallocate(vco_file)
      vco_file => vco_in ! only the analysis levels will be read in the ensemble
      EnsTopMatchesAnlTop = .true.
    else
      write(*,*)
      write(*,*) 'bmat1D_setupBEns: all the vertical levels will be read in the ensemble '
      if ( vco_in%nLev_M > 0 .and. vco_in%vgridPresent ) then
        pSurfRef = 101000.D0
        call czp_fetch1DLevels(vco_in, pSurfRef, profM_opt=pressureProfileInc_M)
        call czp_fetch1DLevels(vco_in, pSurfRef, profM_opt=pressureProfileFile_M)
      
        EnsTopMatchesAnlTop = abs( log(pressureProfileFile_M(1)) - log(pressureProfileInc_M(1)) ) < 0.1d0
        write(*,*) 'bmat1D_setupBEns: EnsTopMatchesAnlTop, presEns, presInc = ', &
             EnsTopMatchesAnlTop, pressureProfileFile_M(1), pressureProfileInc_M(1)
        deallocate(pressureProfileFile_M)
        deallocate(pressureProfileInc_M)
      else
        ! not sure what this mean when no MM levels
        write(*,*) 'bmat1D_setupBEns: nLev_M       = ', vco_in%nLev_M
        write(*,*) 'bmat1D_setupBEns: vgridPresent = ', vco_in%vgridPresent
        EnsTopMatchesAnlTop = .true.
      end if

      if ( EnsTopMatchesAnlTop ) then
        if ( mmpi_myid == 0 ) write(*,*) 'bmat1D_setupBEns: top level of ensemble member and analysis grid match'
        vco_ens => vco_in  ! IMPORTANT: top levels DO match, therefore safe
                           ! to force members to be on analysis vertical levels
      else
        if ( mmpi_myid == 0 ) write(*,*) 'bmat1D_setupBEns: top level of ensemble member and analysis grid are different, therefore'
        if ( mmpi_myid == 0 ) write(*,*) '                      assume member is already be on correct levels - NO CHECKING IS DONE'
        vco_ens => vco_file ! IMPORTANT: top levels do not match, therefore must
                            ! assume file is already on correct vertical levels
      end if
    end if

    if (vco_in%Vcode /= vco_ens%Vcode) then
      write(*,*) 'bmat1D_setupBEns: vco_in%Vcode = ', vco_in%Vcode, ', vco_ens%Vcode = ', vco_ens%Vcode
      call utl_abort('bmat1D_setupBEns: vertical levels of ensemble not compatible with analysis grid')
    end if
    nLevEns_M = vco_ens%nLev_M
    nLevEns_T = vco_ens%nLev_T
    nLevInc_M = vco_in%nLev_M
    nLevInc_T = vco_in%nLev_T

    if (vco_in%Vcode == 5002) then
      if ( (nLevEns_T /= (nLevEns_M+1)) .and. (nLevEns_T /= 1 .or. nLevEns_M /= 1) ) then
        write(*,*) 'bmat1D_setubBEns: nLevEns_T, nLevEns_M = ',nLevEns_T,nLevEns_M
        call utl_abort('bmat1D_setubBEns: Vcode=5002, nLevEns_T must equal nLevEns_M+1!')
      end if
    else if (vco_in%Vcode == 5005) then
      if ( nLevEns_T /= nLevEns_M .and. &
           nLevEns_T /= 0 .and. &
           nLevEns_M /= 0 ) then
        write(*,*) 'bmat1D_setubBEns: nLevEns_T, nLevEns_M = ',nLevEns_T,nLevEns_M
        call utl_abort('bmat1D_setubBEns: Vcode=5005, nLevEns_T must equal nLevEns_M!')
      end if
    else if (vco_in%Vcode == 0) then
      if ( nLevEns_T /= 0 .and. nLevEns_M /= 0 ) then
        write(*,*) 'bmat1D_setubBEns: nLevEns_T, nLevEns_M = ',nLevEns_T, nLevEns_M
        call utl_abort('bmat1D_setubBEns: surface-only case (Vcode=0), nLevEns_T and nLevEns_M must equal 0!')
      end if
    else
      write(*,*) 'bmat1D_setubBEns: vco_in%Vcode = ',vco_in%Vcode
      call utl_abort('bmat1D_setubBEns: unknown vertical coordinate type!')
    end if

    if (nLevEns_M > nLevInc_M) then
      call utl_abort('bmat1D_setubBEns: ensemble has more levels than increment - not allowed!')
    end if

    if (nLevEns_M < nLevInc_M) then
      if (mmpi_myid == 0) write(*,*) 'bmat1D_setubBEns: ensemble has less levels than increment'
      if (mmpi_myid == 0) write(*,*) '                  some levels near top will have zero increment'
    end if

    !- 1.4 Bmatrix Weight
    if (vco_in%Vcode == 5002 .or. vco_in%Vcode == 5005) then
      allocate(scaleFactor_M(nLevEns_M))
      allocate(scaleFactor_T(nLevEns_T))
      do levIndex = 1, nLevEns_T
        if (scaleFactorEns(levIndex) > 0.0d0) then 
          scaleFactorEns(levIndex) = sqrt(scaleFactorEns(levIndex))
        else
          scaleFactorEns(levIndex) = 0.0d0
        end if
      end do
      scaleFactor_T(1:nLevEns_T) = scaleFactorEns(1:nLevEns_T)
      if (vco_in%Vcode == 5002) then
        scaleFactor_M(1:nLevEns_M) = scaleFactorEns(2:(nLevEns_M+1))
      else
        scaleFactor_M(1:nLevEns_M) = scaleFactorEns(1:nLevEns_M)
      end if

      do levIndex = 1, nLevEns_T
        if (scaleFactorEnsHumidity(levIndex) > 0.0d0) then 
          scaleFactorEnsHumidity(levIndex) = sqrt(scaleFactorEnsHumidity(levIndex))
        else
          scaleFactorEnsHumidity(levIndex) = 0.0d0
        end if
      end do
      
      scaleFactor_SF = scaleFactor_T(nLevEns_T)

    else ! vco_in%Vcode == 0
      if (scaleFactorEns(1) > 0.0d0) then 
        scaleFactor_SF = sqrt(scaleFactorEns(1))
      else
        call utl_abort('bmat1D_setubBEns: with vCode == 0, the scale factor should never be equal to 0')
      end if
    end if

    if (scaleFactorEnsTG > 0.0d0) then 
      scaleFactorEnsTG = sqrt(scaleFactorEnsTG)
    else 
      scaleFactorEnsTG = 0.0d0
    end if 

    !- 1.5 Domain Partionning
    call mmpi_setup_latbands(nj, latPerPE, latPerPEmax, myLatBeg, myLatEnd)
    call mmpi_setup_lonbands(ni, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd)

    !- 1.6 Localization
    if ( vLocalize <= 0.0d0 .and. (nLevInc_M > 1 .or. nLevInc_T > 1) ) then
      call utl_abort('bmat1D_setubBEns: Invalid VERTICAL localization length scale')
    end if     
  
    call ens_allocate(ensembles,  &
         nEns, numStep, &
         hco_ens,  &
         vco_ens, dateStampList, &
         hco_core_opt = hco_in, &
         varNames_opt = bmat1D_includeAnlVar(1:bmat1D_numIncludeAnlVar), &
         hInterpolateDegree_opt = hInterpolationDegree)
    write(*,*) 'Read ensemble members'
    call ens_readEnsemble(ensembles, ensPathName, biPeriodic=.false., &
                          varNames_opt = bmat1D_includeAnlVar(1:bmat1D_numIncludeAnlVar))

    allocate(ensColumns(nEns))
    call gsv_allocate(stateVector, numstep, hco_ens, vco_ens, &
                     dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeightSfc_opt=.true.)
    
    call gsv_allocate(stateVectorMean, numstep, hco_ens, vco_ens, &
                     dateStamp_opt=tim_getDateStamp(),  &
                     mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                     dataKind_opt=4, allocHeightSfc_opt=.true.)
    call gsv_zero(stateVectorMean)
    do memberIndex = 1, nEns
      write(*,*) 'Copy member ', memberIndex
      call gsv_zero(stateVector)
      call ens_copyMember(ensembles, stateVector, memberIndex)
      write(*,*) 'interpolate member ', memberIndex
      call col_setVco(ensColumns(memberIndex), vco_ens)
      call col_allocate(ensColumns(memberIndex), obs_numheader(obsSpaceData), &
                        setToZero_opt=.true.)
      call s2c_nl(stateVector, obsSpaceData, ensColumns(memberIndex), hco_in, &
                  timeInterpType='NEAREST', dealloc_opt=.false. )
      call gsv_add(statevector, statevectorMean, scaleFactor_opt=(1.d0/nEns))
    end do

    call col_setVco(meanColumn, vco_ens)
    call col_allocate(meanColumn, obs_numheader(obsSpaceData), &
                      setToZero_opt=.true.)
    call s2c_nl(stateVectorMean, obsSpaceData, meanColumn, hco_in, &
                timeInterpType='NEAREST' )

    call gsv_deallocate(stateVector)
    call gsv_deallocate(stateVectorMean)

    nkgdim = 0
    do varIndex = 1, bmat1D_numIncludeAnlVar
      currentProfile => col_getColumn(meanColumn, var1D_validHeaderIndex(1), varName_opt=bmat1D_includeAnlVar(varIndex))
      nkgdim = nkgdim + size(currentProfile)
    end do
    write(*,*) 'bmat1D_setupBEns: nkgdim', nkgdim
    cvDim_out = nkgdim * var1D_validHeaderCount

    currentProfile => col_getColumn(meanColumn, var1D_validHeaderIndex(1))
    allocate (varLevColFromVarLevBmat(nkgdim))
    allocate (varNameFromVarLevIndexBmat(nkgdim))
    allocate (multFactor(nkgdim))
    varLevIndexBmat = 0
    do varIndex = 1, bmat1D_numIncludeAnlVar
      levIndex = 0
      do varLevIndexCol = 1, size(currentProfile)
        if ( trim( col_getVarNameFromK(meanColumn,varLevIndexCol) ) == trim( bmat1D_includeAnlVar(varIndex) ) ) then
          varLevIndexBmat = varLevIndexBmat + 1
          levIndex = levIndex + 1
          varLevColFromVarLevBmat(varLevIndexBmat) = varLevIndexCol
          varNameFromVarLevIndexBmat(varLevIndexBmat) = trim( bmat1D_includeAnlVar(varIndex) )
          varLevel = vnl_varLevelFromVarname(varNameFromVarLevIndexBmat(varLevIndexBmat))
          if ( varLevel == 'MM' ) then      ! Momentum
            multFactor(varLevIndexBmat) = scaleFactor_M(levIndex)
          else if ( varLevel == 'TH' ) then ! Thermo
            multFactor(varLevIndexBmat) = scaleFactor_T(levIndex)
          else                              ! SF
            multFactor(varLevIndexBmat) = scaleFactor_SF
          end if

          if (varNameFromVarLevIndexBmat(varLevIndexBmat) == 'HU') then
            multFactor(varLevIndexBmat) = multFactor(varLevIndexBmat) * scaleFactorEnsHumidity(levIndex)
          end if

          if (varNameFromVarLevIndexBmat(varLevIndexBmat) == 'TG') then
            multFactor(varLevIndexBmat) = multFactor(varLevIndexBmat) * scaleFactorEnsTG
          end if

          if (mmpi_myid == 0) write(*,*) 'bmat1D_setupBEns:  bmat1D_includeAnlVar ', bmat1D_includeAnlVar(varIndex), varLevIndexBmat, levIndex
        end if
      end do
    end do
    
    call ens_deallocate(ensembles)
    allocate(bSqrtEns(var1D_validHeaderCount,nkgdim,nkgdim))
    bSqrtEns(:,:,:) = 0.d0
    allocate(lineVector(1,nkgdim))

    !$OMP PARALLEL DO PRIVATE (columnIndex,headerIndex,memberIndex,meanProfile,currentProfile,lineVector)
    do columnIndex = 1, var1D_validHeaderCount
      headerIndex = var1D_validHeaderIndex(columnIndex)
      meanProfile => col_getColumn(meanColumn, headerIndex)
      do memberIndex = 1, nEns
        currentProfile => col_getColumn(ensColumns(memberIndex), headerIndex)
        lineVector(1,:) = currentProfile(varLevColFromVarLevBmat(:)) - meanProfile(varLevColFromVarLevBmat(:))
        lineVector(1,:) = lineVector(1,:) * multFactor(:)
        bSqrtEns(columnIndex,:,:) = bSqrtEns(columnIndex,:,:) + &
            matmul(transpose(lineVector),lineVector)
      end do
    end do
    !$OMP END PARALLEL DO

    do varLevIndexBmat =1, nkgdim
      write(*,*) 'bmat1D_setupBEns: variance = ', varLevIndexBmat, &
                 sum(bSqrtEns(:,varLevIndexBmat,varLevIndexBmat))/real((nEns-1)*var1D_validHeaderCount)
    end do

    deallocate(lineVector)
    deallocate(multFactor)
    allocate(meanPressureProfile(nkgdim))
    call lfn_Setup('FifthOrder')

    !$OMP PARALLEL DO PRIVATE (columnIndex,headerIndex,meanPressureProfile,varLevIndexBmat,varLevIndexCol,varLevIndex1,varLevIndex2,varLevel,levIndexColumn,logP1,logP2,zr)
    do columnIndex = 1, var1D_validHeaderCount
      headerIndex = var1D_validHeaderIndex(columnIndex)
      bSqrtEns(columnIndex,:,:) = bSqrtEns(columnIndex,:,:) / (nEns - 1)
      
      do varLevIndexBmat = 1, nkgdim
        varLevIndexCol = varLevColFromVarLevBmat(varLevIndexBmat)
        varLevel = vnl_varLevelFromVarname(varNameFromVarLevIndexBmat(varLevIndexBmat))
        if (varLevel=='SF') then
          meanPressureProfile(varLevIndexBmat) = col_getElem(meanColumn, 1, headerIndex, 'P0')
        else
          levIndexColumn = col_getLevIndexFromVarLevIndex(meanColumn, varLevIndexCol)
          meanPressureProfile(varLevIndexBmat) = col_getPressure(meanColumn, levIndexColumn, headerIndex, varLevel)
        end if
      end do

      ! Apply vertical localization
      if (vLocalize > 0.0d0) then
        do varLevIndex1 = 1, nkgdim 
          logP1 = log(meanPressureProfile(varLevIndex1))
          do varLevIndex2 = 1, nkgdim
            !-  do Schurr product with 5'th order function
            logP2 = log(meanPressureProfile(varLevIndex2))
            zr = abs(logP2 - logP1)
            bSqrtEns(columnIndex, varLevIndex2, varLevIndex1) = &
                bSqrtEns(columnIndex, varLevIndex2, varLevIndex1) * lfn_response(zr, vLocalize)
          end do
        end do
      end if

      ! Apply background error correlation scale factor between TG and other variables
      if (scaleFactorEnsTGCorrelation /= 1.0d0) then
        do varLevIndex1 = 1, nkgdim 
          do varLevIndex2 = 1, nkgdim
            if ((varNameFromVarLevIndexBmat(varLevIndex1) == 'TG' .and. varNameFromVarLevIndexBmat(varLevIndex2) /= 'TG') .or. &
                (varNameFromVarLevIndexBmat(varLevIndex1) /= 'TG' .and. varNameFromVarLevIndexBmat(varLevIndex2) == 'TG')) then 
                bSqrtEns(columnIndex, varLevIndex2, varLevIndex1) = &
                    bSqrtEns(columnIndex, varLevIndex2, varLevIndex1) * scaleFactorEnsTGCorrelation
            end if
          end do
        end do
      end if
    end do
    !$OMP END PARALLEL DO

    do varLevIndexBmat =1, nkgdim
      write(*,*) 'bmat1D_setupBEns: variance (after localization) = ', varLevIndexBmat, &
          sum(bSqrtEns(:,varLevIndexBmat,varLevIndexBmat))/real(var1D_validHeaderCount)
    end do
    
    if (dumpBmatrixTofile) then
      call dumpBmatrices(vco_in, obsSpaceData, dateStampList, bSqrtEns)
    end if
    
    write(*,*) 'bmat1D_setupBEns: computing matrix sqrt'
    do columnIndex = 1, var1D_validHeaderCount
      call utl_matsqrt(bSqrtEns(columnIndex, :, :), nkgdim, 1.d0, printInformation_opt=.false. )
    end do
    
    deallocate(varLevColFromVarLevBmat) 
    deallocate(varNameFromVarLevIndexBmat)
    deallocate(meanPressureProfile)
    call col_deallocate(meanColumn)
    do memberIndex = 1, nEns
      call col_deallocate(ensColumns(memberIndex))
    end do
    cvDim_mpilocal = cvDim_out
    initialized = .true.
    
    if (mmpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
    if (mmpi_myid == 0) write(*,*) 'bmat1D_setupBEns: Exiting'
    
  end subroutine bmat1D_setupBEns
  
  !--------------------------------------------------------------------------
  ! dumpBmatrices
  !--------------------------------------------------------------------------
  subroutine dumpBmatrices(vco_in, obsSpaceData, dateStampList, bEns)
    !
    ! :Purpose: output B matrices or their average to binary file Bmatrix.bin
    !
    implicit none
    
    ! Arguments:
    type(struct_vco), pointer, intent(in)    :: vco_in           ! Analysis vertical grid descriptor
    type (struct_obs)        , intent(inout) :: obsSpaceData     ! ObsSpaceData structure
    integer                  , intent(in)    :: dateStamplist(:) ! List of datestamps
    real(8)                  , intent(in)    :: bEns(:,:,:)      ! B matrices (first dimension is location index)

    ! Locals:
    integer :: nulmat, ierr
    integer :: yyyymmdd, hhmm, countDumped, countDumpedMax, countDumpedMpiGlobal
    integer :: globalDumpedIndex, countDumpedOut
    integer :: columnIndex, headerIndex
    integer :: taskIndex, dumpedIndex
    integer :: tag, status
    integer :: numstep, nkgdim
    integer, external ::  fnom, fclos, newdate
    integer              :: obsOffset(0:mmpi_nprocs-1)
    integer              :: countDumpedAllTasks(mmpi_nprocs)
    integer, allocatable :: listColumnDumped(:)
    real(8) :: latitude, longitude
    real(8), allocatable :: tempoBmatrix(:,:)
    real(8), allocatable :: outLats(:),outLons(:), outBmatrix(:,:,:)

    if (mmpi_myid == 0) write(*,*) 'dumpBmatrices: Starting'
    if (mmpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
    
    countDumped = 0
    do columnIndex = 1, var1D_validHeaderCount 
      headerIndex = var1D_validHeaderIndex(columnIndex)
      latitude = obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex) * MPC_DEGREES_PER_RADIAN_R8
      longitude = obs_headElem_r(obsSpaceData, OBS_LON, headerIndex) * MPC_DEGREES_PER_RADIAN_R8
      if (latitude  >= latMin .and. &
          latitude  <= latMax .and. &
          longitude >= lonMin .and. &
          longitude <= lonMax) then
        countDumped = countDumped + 1
      end if
    end do
    call rpn_comm_barrier('GRID', ierr)
    
    call rpn_comm_allgather(countDumped, 1, 'MPI_INTEGER', countDumpedAllTasks, 1,'MPI_INTEGER', 'GRID', ierr)
    countDumpedMpiGlobal = sum( countDumpedAllTasks(:) )
    countDumpedMax = maxval( countDumpedAllTasks(:) )
    obsOffset(0) = 0
    do taskIndex = 1, mmpi_nprocs - 1
      obsOffset(taskIndex) = obsOffset(taskIndex - 1) + countDumpedAllTasks(taskIndex)
    end do
    write(*,*) 'dumpBmatrices: obsOffset: ', obsOffset(:)
    
    if (countDumpedMax > 0) then
      allocate(listColumnDumped(countDumpedMax))
      listColumnDumped(:) = -1
      countDumped = 0
      do columnIndex = 1, var1D_validHeaderCount 
        headerIndex = var1D_validHeaderIndex(columnIndex)
        latitude = obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex) * MPC_DEGREES_PER_RADIAN_R8
        longitude = obs_headElem_r(obsSpaceData, OBS_LON, headerIndex) * MPC_DEGREES_PER_RADIAN_R8
        if (latitude  >= latMin .and. &
            latitude  <= latMax .and. &
            longitude >= lonMin .and. &
            longitude <= lonMax) then
          countDumped = countDumped + 1
          listColumnDumped(countDumped) = columnIndex
        end if
      end do
      nkgdim = size(bEns, dim=2)
      if (mmpi_myId == 0) then
        nulmat = 0
        ierr = fnom(nulmat, './Bmatrix.bin', 'FTN+SEQ+UNF', 0)
        numstep = size(dateStampList)
        ierr = newdate(dateStampList(1 + numstep / 2), yyyymmdd, hhmm, -3)
        countDumpedOut = countDumped
        if (doAveraging) countDumpedOut = 1
        write(nulmat) yyyymmdd * 100 + nint(hhmm/100.), vco_in%nlev_T, vco_in%nlev_M, vco_in%Vcode, &
            vco_in%ip1_sfc, vco_in%ip1_T_2m, vco_in%ip1_M_10m, bmat1D_numIncludeAnlVar, nkgdim, countDumpedOut
        write(nulmat) vco_in%ip1_T(:), vco_in%ip1_M(:), bmat1D_includeAnlVar(1:bmat1D_numIncludeAnlVar)
        allocate(outLats(countDumpedMpiGlobal), outLons(countDumpedMpiGlobal), OutBmatrix(countDumpedMpiGlobal, nkgdim, nkgdim))
        allocate(tempoBmatrix(nkgdim, nkgdim))
      end if
    end if
    do dumpedIndex = 1, countDumpedMax
      columnIndex = listColumnDumped(dumpedIndex)
      if (columnIndex > 0) then
        headerIndex = var1D_validHeaderIndex(columnIndex)
        latitude = obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex)
        longitude = obs_headElem_r(obsSpaceData, OBS_LON, headerIndex)
      else
        latitude = MPC_missingValue_R8
        longitude = MPC_missingValue_R8
      end if
      if (mmpi_myId == 0) then
        if (latitude /= MPC_missingValue_R8 .and. longitude /= MPC_missingValue_R8) then 
          outBmatrix(dumpedIndex, :, :) = bEns(columnIndex, :, :)
          outLats(dumpedIndex) = latitude
          outLons(dumpedIndex) = longitude
        end if
        do taskIndex = 1,  mmpi_nprocs - 1
          tag = 3 * taskIndex
          call rpn_comm_recv(latitude, 1, 'mpi_real8', taskIndex, tag, 'GRID', status, ierr)
          call rpn_comm_recv(longitude, 1, 'mpi_real8', taskIndex, tag+1, 'GRID', status, ierr)
          if (latitude /= MPC_missingValue_R8 .and. longitude /= MPC_missingValue_R8) then
            call rpn_comm_recv(tempoBmatrix(:,:), nkgdim*nkgdim, 'mpi_real8', taskIndex, tag + 2, 'GRID', status, ierr)
            globalDumpedIndex = dumpedIndex + obsOffset(taskIndex)
            outBmatrix(globalDumpedIndex, :, :) = tempoBmatrix(:, :)
            outLats(globalDumpedIndex) = latitude
            outLons(globalDumpedIndex) = longitude
          end if
        end do
      else
        tag = 3 * mmpi_myID
        call rpn_comm_send(latitude, 1, 'mpi_real8', 0, tag, 'GRID', ierr)
        call rpn_comm_send(longitude, 1, 'mpi_real8', 0, tag + 1, 'GRID', ierr)
        if (columnIndex > 0) then
          call rpn_comm_send(bEns(columnIndex, :, :), nkgdim*nkgdim, 'mpi_real8', 0, tag + 2, 'GRID', ierr)
        end if
      end if
    end do

    if (mmpi_myId == 0) then
      if (doAveraging) then
        write(nulmat) real(sum(outLats(:))/countDumpedMpiGlobal, kind=4), real(sum(outLons(:))/countDumpedMpiGlobal, kind=4), &
            sum(outBmatrix(:, :, :),dim=1)/countDumpedMpiGlobal
      else
        do dumpedIndex = 1, countDumpedMpiGlobal          
          write(nulmat) real(outLats(dumpedIndex), kind=4), real(outLons(dumpedIndex), kind=4), outBmatrix(dumpedIndex, :, :)
        end do
      end if
      ierr = fclos(nulmat)
      deallocate(outLats)
      deallocate(outLons)
      deallocate(outBmatrix)
      deallocate(tempoBmatrix)
    end if
    if (allocated(listColumnDumped)) deallocate(listColumnDumped)
  end subroutine dumpBmatrices

  !--------------------------------------------------------------------------
  ! bmat1D_bSqrtHi
  !--------------------------------------------------------------------------
  subroutine bmat1D_bSqrtHi(controlVector_in, column, obsSpaceData)
    !
    ! :Purpose: HI component of B square root in 1DVar mode
    !
    implicit none

    ! Arguments:
    real(8),                 intent(in)    :: controlVector_in(cvDim_mpilocal)
    type(struct_columnData), intent(inout) :: column
    type(struct_obs),        intent(in)    :: obsSpaceData

    ! Locals:
    integer :: headerIndex, latitudeBandIndex(1), varIndex, columnIndex
    real(8), pointer :: currentColumn(:)
    real(8), allocatable ::  oneDProfile(:)
    real(8) :: latitude
    integer :: surfaceType, offset

    if (mmpi_myid == 0) write(*,*) 'bmat1D_bsqrtHi: starting'
    if (mmpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

    if (.not. initialized) then
      if (mmpi_myid == 0) write(*,*) 'bmat1D_bsqrtHi: 1Dvar B matrix not initialized'
      return
    end if
    allocate(oneDProfile(nkgdim))

    !$OMP PARALLEL DO PRIVATE (columnIndex,headerIndex,oneDProfile,offset,varIndex,currentColumn,latitude,surfaceType,latitudeBandIndex)
    do columnIndex = 1, var1D_validHeaderCount 
      headerIndex = var1D_validHeaderIndex(columnIndex)
      latitude = obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex) !radian 
      surfaceType = tvs_ChangedStypValue(obsSpaceData, headerIndex)
      if (surfaceType == 1) then !Sea
        latitudeBandIndex = minloc( abs( latitude - latSea(:)) )
        oneDProfile(:) = matmul(bSqrtSea(latitudeBandIndex(1), :, :), controlVector_in(1+(columnIndex-1)*nkgdim:columnIndex*nkgdim))
      else ! Land or Sea Ice
        latitudeBandIndex = minloc( abs( latitude - latLand(:)) )
        oneDProfile(:) = matmul(bSqrtLand(latitudeBandIndex(1), :, :), controlVector_in(1+(columnIndex-1)*nkgdim:columnIndex*nkgdim))
      end if
      offset = 0
      do varIndex = 1, bmat1D_numIncludeAnlVar
        currentColumn => col_getColumn(column, headerIndex, varName_opt=bmat1D_includeAnlVar(varIndex))
        currentColumn(:) = currentColumn(:) + oneDProfile(offset+1:offset+size(currentColumn))
        offset = offset + size(currentColumn)
      end do
      if (offset /= nkgdim) then
        write(*,*) 'bmat1D_bsqrtHi: offset, nkgdim', offset, nkgdim
        call utl_abort('bmat1D_bSqrtHi: inconsistency between Bmatrix and statevector size')
      end if
    end do
    !$OMP END PARALLEL DO

    deallocate(oneDProfile)
    if (mmpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
    if (mmpi_myid == 0) write(*,*) 'bmat1D_bSqrtHi: done'

  end subroutine bmat1D_bSqrtHi

  !--------------------------------------------------------------------------
  ! bmat1D_bSqrtHiAd
  !--------------------------------------------------------------------------
  subroutine bmat1D_bSqrtHiAd(controlVector_in, column, obsSpaceData)
    !
    ! :Purpose: HI component of B square root adjoint in 1DVar mode
    !
    implicit none

    ! Arguments:
    real(8),                 intent(inout) :: controlVector_in(cvDim_mpilocal)
    type(struct_columnData), intent(inout) :: column
    type (struct_obs),       intent(in)    :: obsSpaceData

    ! Locals:
    integer :: headerIndex, latitudeBandIndex(1), varIndex, columnIndex
    real(8), pointer :: currentColumn(:)
    real(8), allocatable ::  oneDProfile(:)
    real(8) :: latitude
    integer :: surfaceType, offset

    if (mmpi_myid == 0) write(*,*) 'bmat1D_bSqrtHiAd: starting'
    if (mmpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
    if (.not. initialized) then
      if (mmpi_myid == 0) write(*,*) 'bmat1D_bSqrtHiAd: 1dvar Bmatrix not initialized'
      return
    end if
    allocate(oneDProfile(nkgdim))

    controlVector_in(:) = 0.d0

    !$OMP PARALLEL DO PRIVATE (columnIndex,headerIndex,oneDProfile,offset,varIndex,currentColumn,latitude,surfaceType,latitudeBandIndex)
    do columnIndex = 1, var1D_validHeaderCount
      headerIndex = var1D_validHeaderIndex(columnIndex)
      offset = 0
      do varIndex = 1, bmat1D_numIncludeAnlVar
        currentColumn => col_getColumn(column, headerIndex, varName_opt=bmat1D_includeAnlVar(varIndex))
        oneDProfile(offset+1:offset+size(currentColumn)) = currentColumn(:)
        offset = offset + size(currentColumn)
      end do
      if (offset /= nkgdim) then
        write(*,*) 'bmat1D_bSqrtHiAd: offset, nkgdim', offset, nkgdim
        call utl_abort('bmat1D_bSqrtHiAd: inconsistency between Bmatrix and statevector size')
      end if
      latitude = obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex) !radians
      surfaceType =  tvs_ChangedStypValue(obsSpaceData, headerIndex)
      if (surfaceType == 1) then !Sea
        latitudeBandIndex = minloc( abs( latitude - latSea(:)) )
        controlVector_in(1+(columnIndex-1)*nkgdim:columnIndex*nkgdim) =  &
             controlVector_in(1+(columnIndex-1)*nkgdim:columnIndex*nkgdim) + &
             matmul(bSqrtSea(latitudeBandIndex(1),:,:), oneDProfile)
      else ! Land or Sea Ice
        latitudeBandIndex = minloc( abs( latitude - latLand(:)) )
        controlVector_in(1+(columnIndex-1)*nkgdim:columnIndex*nkgdim) =  &
             controlVector_in(1+(columnIndex-1)*nkgdim:columnIndex*nkgdim) + &
             matmul(bSqrtLand(latitudeBandIndex(1),:,:), oneDProfile)
      end if
    end do
    !$OMP END PARALLEL DO

    deallocate(oneDProfile)
    if (mmpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
    if (mmpi_myid == 0) write(*,*) 'bmat1D_bSqrtHiAd: done'

  end subroutine bmat1D_bSqrtHiAd

  !--------------------------------------------------------------------------
  ! bmat1D_bSqrtEns
  !--------------------------------------------------------------------------
  subroutine bmat1D_bSqrtEns(controlVector_in, column)
    !
    ! :Purpose: Ensemble component of B square root in 1DVar mode
    !
    implicit none

    ! Arguments:
    real(8),                 intent(in)    :: controlVector_in(cvDim_mpilocal)
    type(struct_columnData), intent(inout) :: column

    ! Locals:
    integer :: headerIndex, varIndex, columnIndex
    real(8), pointer :: currentColumn(:)
    real(8), allocatable ::  oneDProfile(:)
    integer :: offset

    allocate(oneDProfile(nkgdim))

    !$OMP PARALLEL DO PRIVATE (columnIndex,headerIndex,oneDProfile,offset,varIndex,currentColumn)
    do columnIndex = 1, var1D_validHeaderCount 
      headerIndex = var1D_validHeaderIndex(columnIndex)
      oneDProfile(:) = matmul(bSqrtEns(columnIndex, :, :), controlVector_in(1+(columnIndex-1)*nkgdim:columnIndex*nkgdim))
      offset = 0
      do varIndex = 1, bmat1D_numIncludeAnlVar
        currentColumn => col_getColumn(column, headerIndex, varName_opt=bmat1D_includeAnlVar(varIndex))
        currentColumn(:) = currentColumn(:) + oneDProfile(offset+1:offset+size(currentColumn))
        offset = offset + size(currentColumn)
      end do
      if (offset /= nkgdim) then
        write(*,*) 'bmat1D_bsqrtEns: offset, nkgdim', offset, nkgdim
        call utl_abort('bmat1D_bSqrtEns: inconsistency between Bmatrix and statevector size')
      end if
    end do
    !$OMP END PARALLEL DO

    deallocate(oneDProfile)

  end subroutine bmat1D_bSqrtEns

  !--------------------------------------------------------------------------
  ! bmat1D_bSqrtEnsAd
  !--------------------------------------------------------------------------
  subroutine bmat1D_bSqrtEnsAd(controlVector_in, column)
    !
    ! :Purpose: Ensemble component of B square root in 1DVar mode
    !
    implicit none

    ! Arguments:
    real(8),                 intent(inout) :: controlVector_in(cvDim_mpilocal)
    type(struct_columnData), intent(inout) :: column

    ! Locals:
    integer :: headerIndex, varIndex, columnIndex
    real(8), pointer :: currentColumn(:)
    real(8), allocatable ::  oneDProfile(:)
    integer :: offset

    if (mmpi_myid == 0) write(*,*) 'bmat1D_bSqrtEnsAd: starting'
    if (mmpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
    if (.not. initialized) then
      if (mmpi_myid == 0) write(*,*) 'bmat1D_bSqrtEnsAd: 1dvar Bmatrix not initialized'
      return
    end if
    allocate(oneDProfile(nkgdim))

    controlVector_in(:) = 0.d0

    !$OMP PARALLEL DO PRIVATE (columnIndex,headerIndex,oneDProfile,offset,varIndex,currentColumn)
    do columnIndex = 1, var1D_validHeaderCount
      headerIndex = var1D_validHeaderIndex(columnIndex)
      offset = 0
      do varIndex = 1, bmat1D_numIncludeAnlVar
        currentColumn => col_getColumn(column, headerIndex, varName_opt=bmat1D_includeAnlVar(varIndex))
        oneDProfile(offset+1:offset+size(currentColumn)) = currentColumn(:)
        offset = offset + size(currentColumn)
      end do
      if (offset /= nkgdim) then
        write(*,*) 'bmat1D_bSqrtHiAd: offset, nkgdim', offset, nkgdim
        call utl_abort('bmat1D_bSqrtEnsAd: inconsistency between Bmatrix and statevector size')
      end if
      controlVector_in(1+(columnIndex-1)*nkgdim:columnIndex*nkgdim) =  &
             controlVector_in(1+(columnIndex-1)*nkgdim:columnIndex*nkgdim) + &
             matmul(bSqrtEns(columnIndex,:,:), oneDProfile)
    end do
    !$OMP END PARALLEL DO

    deallocate(oneDProfile)
    if (mmpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
    if (mmpi_myid == 0) write(*,*) 'bmat1D_bSqrtEnsAd: done'

  end subroutine bmat1D_bSqrtEnsAd

  !--------------------------------------------------------------------------
  ! bmat1D_sqrtB
  !-------------------------------------------------------------------------- 
  subroutine bmat1D_sqrtB(controlVector, cvdim, column, obsSpaceData)
    !
    !:Purpose: To transform model state from control-vector space to grid-point
    !          space.    
    !
    implicit none

    ! Arguments:
    integer,                 intent(in)    :: cvdim
    real(8),                 intent(in)    :: controlVector(cvdim)
    type(struct_columnData), intent(inout) :: column
    type(struct_obs),        intent(in)    :: obsSpaceData

    ! Locals:
    integer :: bmatIndex
    real(8), pointer :: subVector(:)

    !
    !- 1.  Compute the analysis increment
    !
    call col_zero(column)

    bmat_loop: do bmatIndex = 1, numBmat
      if ( .not. bmatActive(bmatIndex) ) cycle bmat_loop
      subVector => cvm_getSubVector( controlVector, bmatLabelList(bmatIndex) )

      call utl_tmg_start(50, '--Bmatrix')
      select case( trim(bmatTypeList(bmatIndex)) )
      case ('HI')
        !- 1.1 Time-Mean Homogeneous and Isotropic...
        call utl_tmg_start(52, '----B_HI_TL')
        call bmat1D_bsqrtHi(subVector,   & ! IN
                            column,      & ! OUT
                            obsspacedata ) ! IN
        call utl_tmg_stop(52)
      case ('ENS')
        !- 1.2 Ensemble based
        call utl_tmg_start(57, '----B_ENS_TL')
        call bmat1D_bsqrtEns(subVector, &  ! IN
                              column)      ! OUT
        call utl_tmg_stop(57)
      case default
        call utl_abort( 'bmat1D_sqrtB: requested bmatrix type does not exist ' // trim(bmatTypeList(bmatIndex)) )
      end select
      call utl_tmg_stop(50)

    end do bmat_loop

  end subroutine bmat1D_sqrtB

  !--------------------------------------------------------------------------
  ! bmat1D_sqrtBT
  !--------------------------------------------------------------------------
  subroutine bmat1D_sqrtBT(controlVector, cvdim, column, obsSpaceData)
    !
    !:Purpose: To transform model state from grid-point space to
    !          error-covariance space.
    !
    implicit none

    ! Arguments:
    integer,                 intent(in)    :: cvdim
    real(8),                 intent(in)    :: controlVector(cvdim)
    type(struct_columnData), intent(inout) :: column
    type(struct_obs),        intent(in)    :: obsSpaceData

    ! Locals:
    integer :: bmatIndex
    real(8), pointer :: subVector(:)

    ! Process components in opposite order as forward calculation
    bmat_loop: do bmatIndex = numBmat, 1, -1
      if ( .not. bmatActive(bmatIndex) ) cycle bmat_loop
      subVector => cvm_getSubVector( controlVector, bmatLabelList(bmatIndex) )

      call utl_tmg_start(50, '--Bmatrix')
      select case( trim(bmatTypeList(bmatIndex)) )
      case ('HI')
        !- Time-Mean Homogeneous and Isotropic...
        call utl_tmg_start(53, '----B_HI_AD')
        call bmat1D_bsqrtHiAd(subvector,  &  ! IN
                              column,     &  ! OUT
                              obSSpaceData ) ! IN
        call utl_tmg_stop(53)
      case ('ENS')
        !- Ensemble based
        call utl_tmg_start(61, '----B_ENS_AD')
        call bmat1D_bsqrtEnsAd(subvector, &  ! IN
                                column )     ! OUT
        call utl_tmg_stop(61)
      case default
        call utl_abort( 'bmat1D_sqrtBT: requested bmatrix type does not exist ' // trim(bmatTypeList(bmatIndex)) )
      end select
      call utl_tmg_stop(50)

    end do bmat_loop

  end subroutine bmat1D_sqrtBT

  !--------------------------------------------------------------------------
  ! bmat1D_get1DVarIncrement
  !--------------------------------------------------------------------------
  subroutine bmat1D_get1DVarIncrement(incr_cv, column, columnTrlOnAnlIncLev, &
                                     obsSpaceData, nvadim_mpilocal)
    !
    ! :Purpose: to compute 1Dvar increment from control vector
    !
    implicit none

    ! Arguments:
    real(8),                 intent(in)    :: incr_cv(:)
    type(struct_columnData), intent(inout) :: column
    type(struct_columnData), intent(in)    :: columnTrlOnAnlIncLev
    type(struct_obs),        intent(in)    :: obsSpaceData
    integer,                 intent(in)    :: nvadim_mpilocal

    ! compute increment from control vector (multiply by B^1/2)
    call bmat1D_sqrtB(incr_cv, nvadim_mpilocal, column, obsSpaceData)
    call cvt_transform(column, 'ZandP_tl', columnTrlOnAnlIncLev)

  end subroutine bmat1D_get1DVarIncrement

  !--------------------------------------------------------------------------
  ! bmat1D_Finalize
  !--------------------------------------------------------------------------
  subroutine bmat1D_Finalize()
    !
    ! :Purpose: to deallocate memory used by internal module structures
    !
    implicit none

    if (initialized) then
       if (allocated(bSqrtLand)) then
         deallocate( bSqrtLand )
         deallocate( bSqrtSea )
         deallocate( latLand, lonLand, latSea, lonSea )
       end if
       if (allocated(bSqrtEns)) deallocate( bSqrtEns )
       call var1D_finalize()
    end if

  end subroutine bmat1D_Finalize

end module bMatrix1DVar_mod
