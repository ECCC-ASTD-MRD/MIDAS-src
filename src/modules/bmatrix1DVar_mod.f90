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

module bmatrix1DVar_mod
  ! MODULE bmatrix1DVar_mod (prefix='bmat1D' category='3. High-level transformations')
  !
  ! :Purpose: contains all 1Dvar B matrices.
  !
  use columnData_mod
  use columnVariableTransforms_mod
  use controlVector_mod
  use gridStatevector_mod
  use horizontalCoord_mod
  use mpi_mod 
  use obsSpaceData_mod
  use timeCoord_mod
  use utilities_mod
  use verticalCoord_mod
  use codeprecision_mod
  use tovs_nl_mod
  use mathphysconstants_mod
  use var1D_mod
!  use localizationFunction_mod
!  use mpivar_mod

  implicit none
  save
  private

  ! public procedures
  public :: bmat1D_bsetup
  public :: bmat1D_sqrtB, bmat1D_sqrtBT
  public :: bmat1D_finalize, bmat1D_get1DVarIncrement

  ! public variables
  public :: bmat1D_varList

  type(struct_hco), pointer :: hco_yGrid
  logical             :: initialized = .false.
  integer             :: nkgdim
  integer             :: cvDim_mpilocal
  integer, parameter   :: maxNumLevels=200
 
  real(8), allocatable :: bMatrix(:,:)
  real(8), allocatable :: bSqrtLand(:,:,:), bSqrtSea(:,:,:)
  real(4), allocatable :: latLand(:), lonLand(:), latSea(:), lonSea(:)
  integer              :: nLonLatPosLand, nLonLatPosSea, bmat1D_varCount
  character(len=4), allocatable :: bmat1D_varList(:)
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
  integer          :: nEns
  real(8)          :: vlocalize
  real(8),allocatable :: LvertSqrt(:,:)
  !logical          :: useBHi, useBEns
  !Namelist variables
  real(8)             :: scaleFactorHI(maxNumLevels)    ! scaling factors for HI variances
  real(8)             :: scaleFactorHILQ(maxNumLevels)  ! scaling factors for HI LQ variances
  real(8)             :: scaleFactorEns(maxNumLevels)   ! scaling factors for Ens variances
  real(8)             :: scaleFactorEnsHumidity(maxNumLevels) ! scaling factors for Ens LQ variances
  NAMELIST /NAMBMAT1D/ scaleFactorHI, scaleFactorHILQ, scaleFactorENs, scaleFactorEnsHumidity, nEns, &
       vLocalize

contains

  !--------------------------------------------------------------------------
  ! bmat1D_bsetup
  !--------------------------------------------------------------------------
  subroutine bmat1D_bsetup(vco_in, obsSpaceData)
    !
    !:Purpose: To initialize the 1Dvar analysis Background term.
    !
    implicit none
    ! arguments:
    type(struct_vco), pointer, intent(in) :: vco_in
    type (struct_obs), intent(in)         :: obsSpaceData
    ! locals:
    integer :: cvdim
    integer :: masterBmatIndex, bmatIndex
    logical :: active
    integer :: nulnam, ierr
    integer, external ::  fnom, fclos
    ! default values for namelist variables
    scaleFactorHI(:) = 1.d0
    scaleFactorHILQ(:) = 1.d0
    scaleFactorEns(:) = 1.d0
    scaleFactorEnsHumidity(:) = 1.d0
    nEns = -1
    vLocalize = -1.d0
    nulnam = 0
    ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    read(nulnam, nml=nambmat1D, iostat=ierr)
    if ( ierr /= 0 ) call utl_abort( 'bmat1D_bsetup: Error reading namelist' )
    if ( mpi_myid == 0 ) write( *, nml = nambmat1D )
    ierr = fclos( nulnam )
    !
    !- 1.  Setup the B matrices
    !
    do masterBmatIndex = 1, numMasterBmat
      select case( trim(masterBmatTypeList(masterBmatIndex)) )
      case ('HI')
        !- 1.1 Time-Mean Homogeneous and Isotropic...
        write(*,*) 'bmat1D_bsetup: Setting up the modular GLOBAL HI 1D covariances...'
        call bmat1D_SetupBHi(vco_in, obsSpaceData, cvdim)
        write(*,*) " bmat1D_bsetup: cvdim= ", cvdim
      case ('ENS')
        !- 1.2 ensemble based
         write(*,*) 'bmat1D_bsetup: Setting up the ensemble based 1D matrix.'
        call bmat1D_SetupBEns(vco_in, obsSpaceData, cvdim)
        write(*,*) " bmat1D_bsetup: cvdim= ", cvdim
      case default
        call utl_abort( 'bmat1D_bSetup: requested bmatrix type does not exist ' // trim(masterBmatTypeList(masterBmatIndex)) )
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
    write(*,*) "bmat1D_bsetup SUMMARY, number of B matrices found = ", numBmat
    do bmatIndex = 1, numBmat
      write(*,*) "  B matrix #", bmatIndex
      active = cvm_subVectorExists(bmatLabelList(bmatIndex))
      if (active) then
        write(*,*) "   ACTIVE"
      else
        write(*,*) "   NOT USED"
      end if
      write(*,*) "     -> label       = ", bmatLabelList (bmatIndex)
      write(*,*) "     -> type        = ", bmatTypeList  (bmatIndex)
      if (active) then
        write(*,*) "     -> is 3D       = ", bmatIs3dList  (bmatIndex)
      end if
      bmatActive(bmatIndex) = active
    end do

  end subroutine bmat1D_bsetup

  !--------------------------------------------------------------------------
  !  bmat1D_setupBHi
  !--------------------------------------------------------------------------
  subroutine bmat1D_setupBHi(vco_in, obsSpaceData, cvDim_out)
    !
    ! :Purpose: to setup bmat1D module
    !
    implicit none
    ! arguments:
    type(struct_vco), pointer, intent(in):: vco_in
    type (struct_obs), intent(in)        :: obsSpaceData
    integer, intent(out)                 :: cvDim_out
    ! locals:
    integer :: levelIndex, nulnam, ierr
    integer, external ::  fnom, fclos
    integer :: status, Vcode_anl
    logical :: fileExists
    integer :: nulbgst=0
    type(struct_vco), pointer :: vco_file => null()
    type(struct_vco), target  :: vco_1Dvar
    type(struct_vco), pointer :: vco_anl
    character(len=18) :: oneDBmatLand = './Bmatrix_land.bin'
    character(len=17) :: oneDBmatSea  = './Bmatrix_sea.bin'
    integer :: extractDate, locationIndex, countGood, headerIndex
    integer :: bodyStart, bodyEnd, bodyIndex
    logical,save :: firstCall=.true.

    if (firstCall) then
      call var1D_setup(vco_in, obsSpaceData)
      firstCall = .false.
    end if

    call tmg_start(15,'BHI1D_SETUP')
    if(mpi_myid == 0) write(*,*) 'bmat1D_setupBHi: Starting'
    if(mpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

    do levelIndex = 1, maxNumLevels
      if( scaleFactorHI( levelIndex ) > 0.0d0 ) then 
        scaleFactorHI( levelIndex ) = sqrt( scaleFactorHI( levelIndex ))
      else
        scaleFactorHI( levelIndex ) = 0.0d0
      end if
    end do
    
    do levelIndex = 1, maxNumLevels
      if(scaleFactorHILQ(levelIndex) > 0.0d0) then 
        scaleFactorHILQ(levelIndex) = sqrt(scaleFactorHILQ(levelIndex))
      else
        scaleFactorHILQ(levelIndex) = 0.0d0
      end if
    end do
    if ( sum( scaleFactorHI( 1 : maxNumLevels ) ) == 0.0d0 ) then
      if ( mpi_myid == 0 ) write(*,*) 'bmat1D_setupBHi: scaleFactorHI=0, skipping rest of setup'
      cvDim_out = 0
      call tmg_stop(15)
      return
    end if
    if (mpi_myid == 0) write(*,*) 'bmat1D_setupBHi: Read 1DVar background statistics'
    inquire(file=trim(oneDBmatLand), exist=fileExists)
    if ( fileExists ) then
      ierr = fnom(nulbgst, trim(oneDBmatLand), 'FTN+SEQ+UNF+OLD+R/O', 0)
    else
      call utl_abort('bmat1D_setupBHi: No 1DVar BACKGROUND STAT FILE ' // trim(oneDBmatLand))
    end if
    read(nulbgst) extractDate, vco_1Dvar%nLev_T, vco_1Dvar%nLev_M, vco_1Dvar%Vcode, &
         vco_1Dvar%ip1_sfc, vco_1Dvar%ip1_T_2m, vco_1Dvar%ip1_M_10m, bmat1D_varCount, nkgdim, nLonLatPosLand
    allocate( vco_1Dvar%ip1_T(vco_1Dvar%nLev_T), vco_1Dvar%ip1_M(vco_1Dvar%nLev_M) )
    allocate( bmat1D_varList(bmat1D_varCount) )
    allocate( bMatrix(nkgdim,nkgdim) )
    allocate( latLand(nLonLatPosLand), lonLand(nLonLatPosLand))       
    allocate( bSqrtLand(nLonLatPosLand, nkgdim, nkgdim) )
    read(nulbgst) vco_1Dvar%ip1_T(:), vco_1Dvar%ip1_M(:), bmat1D_varList(:)
    do locationIndex = 1, nLonLatPosLand
      read(nulbgst) latLand(locationIndex), lonLand(locationIndex), bMatrix(:,:)      
      bSqrtLand(locationIndex, :, :) = bMatrix(:, :)
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
         vco_1Dvar%ip1_sfc, vco_1Dvar%ip1_T_2m, vco_1Dvar%ip1_M_10m, bmat1D_varCount, nkgdim, nLonLatPosSea
    allocate( bSqrtSea(nLonLatPosSea, nkgdim, nkgdim) )
    allocate( latSea(nLonLatPosSea), lonSea(nLonLatPosSea))
    read(nulbgst) vco_1Dvar%ip1_T(:), vco_1Dvar%ip1_M(:), bmat1D_varList(:)
    do locationIndex = 1, nLonLatPosSea
      read(nulbgst) latSea(locationIndex), lonSea(locationIndex), bMatrix(:,:)
      bSqrtSea(locationIndex, :, :) = bMatrix(:, :)
      call utl_matsqrt(bSqrtSea(locationIndex,:,:), nkgdim, 1.d0, printInformation_opt=.false. )
    end do
    ierr = fclos(nulbgst)

    vco_1Dvar%initialized = .true.
    vco_1Dvar%vGridPresent = .false.
    vco_file => vco_1Dvar
    vco_anl => vco_in
    if (.not. vco_equal(vco_anl,vco_file)) then
      call utl_abort('bmat1D_setupBHi: vco from analysisgrid and cov file do not match')
    end if
    if (mpi_myid == 0) write(*,*) 'bmat1D_setupBHi: nLev_M, nLev_T=', vco_1Dvar%nLev_M, vco_1Dvar%nLev_T
    status = vgd_get(vco_anl%vgrid, key='ig_1 - vertical coord code', value=Vcode_anl)
    if(Vcode_anl /= 5002 .and. Vcode_anl /= 5005) then
      write(*,*) 'Vcode_anl = ',Vcode_anl
      call utl_abort('bmat1D_setupBHi: unknown vertical coordinate type!')
    end if
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

    cvDim_out = nkgdim * var1D_validHeaderCount
    cvDim_mpilocal = cvDim_out
    initialized = .true.

    if(mpi_myid == 0) write(*,*) 'bmat1D_setupBHi: Exiting'
    if(mpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
    call tmg_stop(15)

  end subroutine bmat1D_setupBHi

  !--------------------------------------------------------------------------
  !  bmat1D_setupBEns
  !--------------------------------------------------------------------------
  subroutine bmat1D_setupBEns(vco_in, obsSpaceData, cvDim_out)
    !
    ! :Purpose: to setup bmat1D module
    !
    implicit none
    ! arguments:
    type(struct_vco), pointer, intent(in):: vco_in
    type (struct_obs), intent(in)        :: obsSpaceData
    integer, intent(out)                 :: cvDim_out
    ! locals:
    

    call tmg_start(15,'BENS1D_SETUP')
    if(mpi_myid == 0) write(*,*) 'bmat1D_setupBEns: Starting'
    if(mpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

    cvDim_out = 0
    cvDim_mpilocal = cvDim_out
    initialized = .true.

  end subroutine bmat1D_setupBEns
  
  !--------------------------------------------------------------------------
  ! bmat1D_bSqrtHi
  !--------------------------------------------------------------------------
  subroutine bmat1D_bSqrtHi(controlVector_in, column, obsSpaceData)
    !
    ! :Purpose: HI component of B square root in 1DVar mode
    !
    implicit none
    ! arguments:
    real(8), intent(in)                    :: controlVector_in(cvDim_mpilocal)
    type(struct_columnData), intent(inout) :: column
    type(struct_obs), intent(in)           :: obsSpaceData
    ! locals:
    integer :: headerIndex, latitudeBandIndex(1), varIndex, columnIndex
    real(8), pointer :: currentColumn(:)
    real(8), allocatable ::  oneDProfile(:)
    real(8) :: latitude
    integer :: surfaceType, offset

    if (mpi_myid == 0) write(*,*) 'bmat1D_bsqrtHi: starting'
    if (mpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

    if (.not. initialized) then
      if (mpi_myid == 0) write(*,*) 'bmat1D_bsqrtHi: 1Dvar B matrix not initialized'
      return
    end if
    allocate(oneDProfile(nkgdim))
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
      do varIndex = 1, bmat1D_varCount
        currentColumn => col_getColumn(column, headerIndex, varName_opt=bmat1D_varList(varIndex))
        currentColumn(:) = 0.d0
        currentColumn(:) = oneDProfile(offset+1:offset+size(currentColumn))
        offset = offset + size(currentColumn)
      end do
      if (offset /= nkgdim) then
        write(*,*) 'bmat1D_bsqrtHi: offset, nkgdim', offset, nkgdim
        call utl_abort('bmat1D_bSqrtHi: inconsistency between Bmatrix and statevector size')
      end if
    end do
    deallocate(oneDProfile)
    if (mpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
    if (mpi_myid == 0) write(*,*) 'bmat1D_bSqrtHi: done'

  end subroutine bmat1D_bSqrtHi

  !--------------------------------------------------------------------------
  ! bmat1D_bSqrtHiAd
  !--------------------------------------------------------------------------
  subroutine bmat1D_bSqrtHiAd(controlVector_in, column, obsSpaceData)
    !
    ! :Purpose: HI component of B square root adjoint in 1DVar mode
    !
    implicit none
    ! arguments:
    real(8), intent(inout)                 :: controlVector_in(cvDim_mpilocal)
    type(struct_columnData), intent(inout) :: column
    type (struct_obs), intent(in)          :: obsSpaceData
    ! locals:
    integer :: headerIndex, latitudeBandIndex(1), varIndex, columnIndex
    real(8), pointer :: currentColumn(:)
    real(8), allocatable ::  oneDProfile(:)
    real(8) :: latitude
    integer :: surfaceType, offset

    if (mpi_myid == 0) write(*,*) 'bmat1D_bSqrtHiAd: starting'
    if (mpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
    if (.not. initialized) then
      if (mpi_myid == 0) write(*,*) 'bmat1D_bSqrtHiAd: 1dvar Bmatrix not initialized'
      return
    end if
    allocate(oneDProfile(nkgdim))
    do columnIndex = 1, var1D_validHeaderCount
      headerIndex = var1D_validHeaderIndex(columnIndex)
      offset = 0
      do varIndex = 1, bmat1D_varCount
        currentColumn => col_getColumn(column, headerIndex, varName_opt=bmat1D_varList(varIndex))
        oneDProfile(offset+1:offset+size(currentColumn)) = currentColumn(:)
        offset = offset + size(currentColumn)
      end do
      if (offset /= nkgdim) then
        write(*,*) 'bmat1D_bSqrtHiAd: offset, nkgdim', offset, nkgdim
        call utl_abort('bmat1D_bSqrtHiAd: inconsistency between Bmatrix and statevector size')
      end if
      latitude = obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex) !radian
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
    deallocate(oneDProfile)
    if (mpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
    if (mpi_myid == 0) write(*,*) 'bmat1D_bSqrtHiAd: done'

  end subroutine bmat1D_bSqrtHiAd


  !--------------------------------------------------------------------------
  ! bmat1D_bSqrtEns
  !--------------------------------------------------------------------------
  subroutine bmat1D_bSqrtEns(controlVector_in, column, obsSpaceData)
    !
    ! :Purpose: Ensemble component of B square root in 1DVar mode
    !
    implicit none
    ! arguments:
    real(8), intent(in)                    :: controlVector_in(cvDim_mpilocal)
    type(struct_columnData), intent(inout) :: column
    type(struct_obs), intent(in)           :: obsSpaceData
    ! locals:
    integer :: headerIndex, varIndex, columnIndex
    real(8), pointer :: currentColumn(:)
    real(8), allocatable ::  oneDProfile(:)
    integer :: offset

  end subroutine bmat1D_bSqrtEns

  !--------------------------------------------------------------------------
  ! bmat1D_bSqrtEnsAd
  !--------------------------------------------------------------------------
  subroutine bmat1D_bSqrtEnsAd(controlVector_in, column, obsSpaceData)
    !
    ! :Purpose: Ensemble component of B square root in 1DVar mode
    !
    implicit none
    ! arguments:
    real(8), intent(in)                    :: controlVector_in(cvDim_mpilocal)
    type(struct_columnData), intent(inout) :: column
    type(struct_obs), intent(in)           :: obsSpaceData
    ! locals:
    integer :: headerIndex, varIndex, columnIndex
    real(8), pointer :: currentColumn(:)
    real(8), allocatable ::  oneDProfile(:)
    integer :: offset

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
    ! arguments:
    integer, intent(in)                    :: cvdim
    real(8), intent(in)                    :: controlVector(cvdim)
    type(struct_columnData), intent(inout) :: column
    type(struct_obs), intent(in)           :: obsSpaceData
    ! locals:
    integer :: bmatIndex
    real(8), pointer :: subVector(:)
    !
    !- 1.  Compute the analysis increment
    !
    bmat_loop: do bmatIndex = 1, numBmat
      if ( .not. bmatActive(bmatIndex) ) cycle bmat_loop
      subVector => cvm_getSubVector( controlVector, bmatLabelList(bmatIndex) )
      select case( trim(bmatTypeList(bmatIndex)) )
      case ('HI')
        !- 1.1 Time-Mean Homogeneous and Isotropic...
        call tmg_start(50,'B_HI')
        call bmat1D_bsqrtHi( subVector,   &  ! IN
                            column,      &  ! OUT
                            obsspacedata )  ! IN
        call tmg_stop(50)
      case ('ENS')
        !- 1.2 Ensemble based
        call tmg_start(50,'B_ENS')
        call bmat1D_bsqrtEns( subVector,   &  ! IN
                             column,      &  ! OUT
                             obsspacedata )  ! IN
        call tmg_stop(50)
      case default
        call utl_abort( 'bmat1D_sqrtB: requested bmatrix type does not exist ' // trim(bmatTypeList(bmatIndex)) )
      end select
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
    ! arguments:
    integer, intent(in)                    :: cvdim
    real(8), intent(in)                    :: controlVector(cvdim)
    type(struct_columnData), intent(inout) :: column
    type(struct_obs), intent(in)           :: obsSpaceData
    ! locals:
    integer :: bmatIndex
    real(8), pointer :: subVector(:)
    ! Process components in opposite order as forward calculation
    bmat_loop: do bmatIndex = numBmat, 1, -1
      if ( .not. bmatActive(bmatIndex) ) cycle bmat_loop
      subVector => cvm_getSubVector( controlVector, bmatLabelList(bmatIndex) )
      select case( trim(bmatTypeList(bmatIndex)) )
      case ('HI')
        !- Time-Mean Homogeneous and Isotropic...
        call tmg_start(51,'B_HI_T')
        call bmat1D_bsqrtHiAd( subvector, &  ! IN
                              column,    &  ! OUT
                              obSSpaceData )     ! IN
        call tmg_stop(51)
      case ('ENS')
        !- Ensemble based
        call tmg_start(51,'B_ENS_T')
        call bmat1D_bsqrtEnsAd( subvector, &  ! IN
                               column,    &  ! OUT
                               obSSpaceData )     ! IN
        call tmg_stop(51)
      case default
        call utl_abort( 'bmat1D_sqrtBT: requested bmatrix type does not exist ' // trim(bmatTypeList(bmatIndex)) )
      end select
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
    ! arguments:
    real(8), intent(in)                    :: incr_cv(:)
    type(struct_columnData), intent(inout) :: column
    type(struct_columnData), intent(in)    :: columnTrlOnAnlIncLev
    type(struct_obs),        intent(in)    :: obsSpaceData
    integer,                 intent(in)    :: nvadim_mpilocal
    ! compute increment from control vector (multiply by B^1/2)
    call bmat1D_sqrtB(incr_cv, nvadim_mpilocal, column, obsSpaceData)
    call cvt_transform(column, columnTrlOnAnlIncLev, 'PsfcToP_tl')

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
       deallocate( bMatrix)
       deallocate( bSqrtLand )
       deallocate( bSqrtSea )
       deallocate( latLand, lonLand, latSea, lonSea )
       !if (allocated(bSqrtEns)) deallocate( bSqrtEns )
       call var1D_finalize()
    end if

  end subroutine bmat1D_Finalize

end module bmatrix1DVar_mod
