
module bMatrix_mod
  ! MODULE bMatrix_mod (prefix='bmat' category='2. B and R matrices')
  !
  ! :Purpose: A higher-level module that takes care of calling subroutines 
  !           in the lower-level modules bmatrixHI/lambmatrixHI and
  !           bmatrixEnsemble
  !
  ! :Comments:
  !           - Considerations for ensemble-based and regional static
  !             covariances for constituents are not yet included.
  !
  use midasMpi_mod
  use bMatrixHI_mod
  use bMatrixEnsemble_mod
  use bMatrixChem_mod
  use bMatrixDiff_mod
  use controlVector_mod
  use verticalCoord_mod
  use gridStateVector_mod
  use LAMbMatrixHI_mod
  use horizontalCoord_mod
  use utilities_mod
  implicit none
  save
  private

  ! public procedures
  public :: bmat_setup, bmat_finalize, bmat_sqrtB, bmat_sqrtBT
  public :: bmat_reduceToMPILocal, bmat_reduceToMPILocal_r4, bmat_expandToMPIGlobal, bmat_expandToMPIGlobal_r4

  logical :: globalGrid = .true.

  integer,          parameter :: numMasterBmat = 4
  character(len=4), parameter :: masterBmatTypeList (numMasterBmat) = (/'HI'  , 'ENS'  , 'CHM'  , 'DIFF'  /)
  character(len=8), parameter :: masterBmatLabelList(numMasterBmat) = (/'B_HI', 'B_ENS', 'B_CHM', 'B_DIFF'/)
  logical,          parameter :: masterbmatIs3dList (numMasterBmat) = (/.true., .false., .true. , .true.  /)

  integer            :: numBmat
  integer, parameter :: numBmatMax = 50

  character(len=4) :: bmatTypeList  (numBmatMax)
  character(len=9) :: bmatLabelList (numBmatMax)
  integer          :: bmatInstanceID(numBmatMax)
  logical          :: bmatIs3dList  (numBmatMax)
  logical          :: bmatActive    (numBmatMax)

contains

  !--------------------------------------------------------------------------
  ! bmat_setup
  !--------------------------------------------------------------------------
  subroutine bmat_setup(hco_anl, hco_core, vco_anl)
    !
    !:Purpose: To initialize the analysis Background term for the specific
    !          analysis configuration used.
    !
    implicit none

    ! Arguments:
    type(struct_vco), pointer, intent(in) :: vco_anl
    type(struct_hco), pointer, intent(in) :: hco_anl
    type(struct_hco), pointer, intent(in) :: hco_core

    ! Locals:
    integer, allocatable :: cvDimPerInstance(:)
    integer :: cvdim
    integer :: masterBmatIndex, bMatInstanceIndex, nBmatInstance, bmatIndex
    character(len=2) :: bMatInstanceIndexString
    character(len=3) :: bMatExtraLabel
    logical :: active

    call utl_tmg_start(50,'--Bmatrix')

    !
    !- 1.  Get/Check the analysis grid info
    !

    !- 1.1 Horizontal Grid info
    globalGrid = hco_anl%global

    !
    !- 2.  Setup the B matrices
    !
    numBmat = 0

    do masterBmatIndex = 1, numMasterBmat

      select case( trim(masterBmatTypeList(masterBmatIndex)) )
      case ('HI')

        !- 2.1 Time-Mean Homogeneous and Isotropic...
        nBmatInstance = 1 ! hardwired
        allocate(cvdimPerInstance(nBmatInstance))

        if ( globalGrid ) then
          write(*,*)
          write(*,*) 'Setting up the modular GLOBAL HI covariances...'
          call bhi_Setup( hco_anl, vco_anl, & ! IN
                          cvdim )             ! OUT
        else
          write(*,*)
          write(*,*) 'Setting up the modular LAM HI covariances...'
          call lbhi_Setup( hco_anl, hco_core, vco_anl, & ! IN
                           cvdim )                       ! OUT
        end if

        cvdimPerInstance(1) = cvdim

      case ('ENS')

        !- 2.2 Flow-dependent Ensemble-Based
        write(*,*)
        write(*,*) 'Setting up the modular ENSEMBLE covariances...'
        call ben_Setup( hco_anl, hco_core, vco_anl, & ! IN
                        cvdimPerInstance )  ! OUT

        nBmatInstance = size(cvdimPerInstance)

      case ('CHM')

        !- 2.3  Static (Time-Mean Homogeneous and Isotropic) covariances for constituents
        nBmatInstance = 1 ! hardwired
        allocate(cvdimPerInstance(nBmatInstance))

        if ( globalGrid ) then
          write(*,*)
          write(*,*) 'Setting up the modular GLOBAL HI-chm covariances...'
          call bchm_SetupCH( hco_anl, vco_anl, & ! IN
                             cvdim )             ! OUT
        else
          cvdim=0
        end if

        cvdimPerInstance(1) = cvdim

      case ('DIFF')

        !- 2.4 Covariances modelled using a diffusion operator.
        nBmatInstance = 1 ! hardwired
        allocate(cvdimPerInstance(nBmatInstance))

        write(*,*)
        write(*,*) 'Setting up the modular DIFFUSION covariances...'
        call bdiff_Setup( hco_anl, vco_anl, & ! IN
                          cvdim )             ! OUT

        cvdimPerInstance(1) = cvdim

      case default

        call utl_abort( 'bmat_setup: requested bmatrix type does not exist ' // trim(masterBmatTypeList(masterBmatIndex)) )

      end select

      !- 2.5 Append the info to the B matrix info arrays and setup the proper control sub-vectors
      do bMatInstanceIndex = 1, nBmatInstance

        numBmat = numBmat + 1

        if (nBmatInstance == 1) then
          bMatExtraLabel= ""
        else
          write(bMatInstanceIndexString,'(I2.2)') bMatInstanceIndex 
          bMatExtraLabel="_"//trim(bMatInstanceIndexString)
        end if
        bmatLabelList (numBmat) = trim(masterbmatLabelList(masterBmatIndex))//trim(bMatExtraLabel)

        bmatTypeList  (numBmat) = masterBmatTypeList(masterBmatIndex)
        bmatIs3dList  (numBmat) = masterbmatIs3dList(masterBmatIndex)
        bmatInstanceID(numBmat) = bMatInstanceIndex

        call cvm_setupSubVector(bmatLabelList(numBmat), bmatTypeList(numBmat), cvdimPerInstance(bMatInstanceIndex))

      end do

      deallocate(cvdimPerInstance)

    end do

    !
    !- 3. Print a summary and set the active B matrices array
    !
    write(*,*)
    write(*,*) " bmat_setup SUMMARY, number of B matrices found = ", numBmat
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
        write(*,*) "     -> instance ID = ", bmatInstanceID(bmatIndex)
      end if
      bmatActive(bmatIndex) = active
    end do

    !
    !- 4. Check if the control vector has zero length (i.e. no B matrix is active)
    !
    if (cvm_nvadim_mpiglobal == 0) then
      call utl_abort('bmat_setup: The control vector has zero length. No B matrix has been set up.')
    end if

    call utl_tmg_stop(50)

  end subroutine bmat_setup

  !--------------------------------------------------------------------------
  ! bmat_sqrtB
  !-------------------------------------------------------------------------- 
  subroutine bmat_sqrtB(controlVector, cvdim, statevector,  &
                        useFSOFcst_opt, stateVectorRef_opt)
    !
    !:Purpose: To transform model state from control-vector space to grid-point
    !          space.
    !          
    !
    implicit none

    ! Arguments:
    integer                   , intent(in)    :: cvdim
    real(8)                   , intent(in)    :: controlVector(cvdim)
    type(struct_gsv)          , intent(inout) :: statevector
    logical, optional         , intent(in)    :: useFSOFcst_opt
    type(struct_gsv), optional, intent(in)    :: stateVectorRef_opt

    ! Locals:
    integer :: bmatIndex
    real(8),pointer :: subVector(:)
    type(struct_gsv) :: statevector_temp
    character(len=4), pointer :: varNames(:)

    call utl_tmg_start(50,'--Bmatrix')

    !
    !- 1.  Set analysis increment to zero and allocate a temporary statevector
    !
    nullify(varNames)
    call gsv_zero( statevector )
    call gsv_varNamesList(varNames, statevector)
    call gsv_allocate( statevector_temp, statevector%numStep,            &
                       gsv_getHco(statevector), gsv_getVco(statevector), &
                       dataKind_opt=gsv_getDataKind(statevector),        &
                       mpi_local_opt=.true., varNames_opt=varNames )
    deallocate(varNames)

    !
    !- 2.  Compute the analysis increment
    !
    bmat_loop: do bmatIndex = 1, numBmat

      if ( .not. bmatActive(bmatIndex) ) cycle bmat_loop

      subVector => cvm_getSubVector( controlVector, bmatLabelList(bmatIndex) )
      call gsv_zero( statevector_temp )

      select case( trim(bmatTypeList(bmatIndex)) )
      case ('HI')

        !- 2.1 Time-Mean Homogeneous and Isotropic...
        call utl_tmg_start(52,'----B_HI_TL')
        if ( globalGrid ) then
          call bhi_bsqrt( subVector,        &  ! IN
                          statevector_temp, &  ! OUT
                          stateVectorRef_opt ) ! IN
        else
          call lbhi_bSqrt( subVector,        &  ! IN
                           statevector_temp, &  ! OUT
                           stateVectorRef_opt ) ! IN
        end if
        call utl_tmg_stop(52)

      case ('CHM')

        !- 2.2  Static (Time-Mean Homogeneous and Isotropic) covariances for constituents
        call utl_tmg_start(68,'----B_CHM_TL')
        if ( globalGrid ) then
          call bchm_bsqrt( subVector,        &  ! IN
                           statevector_temp, &  ! OUT
                           stateVectorRef_opt ) ! IN
        end if
        call utl_tmg_stop(68)

      case ('DIFF')

        !- 2.3 Covariances modelled using a diffusion operator.
        call utl_tmg_start(66,'----B_DIFF_TL')
        call bdiff_bsqrt( subVector,       & ! IN
                          statevector_temp ) ! OUT
        call utl_tmg_stop(66)

      case ('ENS')

        !- 2.4 Flow-dependent Ensemble-Based
        call utl_tmg_start(57,'----B_ENS_TL')
        call ben_bsqrt( bmatInstanceID(bmatIndex), subVector, & ! IN
                        statevector_temp,                     & ! OUT
                        useFSOFcst_opt, stateVectorRef_opt )    ! IN
        call utl_tmg_stop(57)

      end select

      ! Make latest increment contribution 4D, if necessary
      if ( bmatIs3dList(bmatIndex) ) call gsv_3dto4d( statevector_temp )

      ! Add latest contribution to total increment in statevector
      call gsv_add( statevector_temp, statevector )

    end do bmat_loop

    call gsv_deallocate( statevector_temp )

    call utl_tmg_stop(50)

  end subroutine bmat_sqrtB

  !--------------------------------------------------------------------------
  ! bmat_sqrtBT
  !--------------------------------------------------------------------------
  subroutine bmat_sqrtBT(controlVector, cvdim, statevector,  &
                         useFSOFcst_opt, stateVectorRef_opt)
    !
    !:Purpose: To transform model state from grid-point space to
    !          error-covariance space.
    !
    implicit none

    ! Arguments:
    integer ,                   intent(in)    :: cvdim
    real(8) ,                   intent(inout) :: controlVector(cvdim)
    type(struct_gsv),           intent(in)    :: statevector
    logical,optional,           intent(in)    :: useFSOFcst_opt
    type(struct_gsv), optional, intent(in)    :: stateVectorRef_opt

    ! Locals:
    integer :: bmatIndex
    real(8),pointer :: subVector(:)
    type(struct_gsv) :: statevector_temp
    character(len=4), pointer :: varNames(:)

    call utl_tmg_start(50,'--Bmatrix')

    nullify(varNames)
    call gsv_varNamesList(varNames, statevector)
    call gsv_allocate( statevector_temp, statevector%numStep,            &
                       gsv_getHco(statevector), gsv_getVco(statevector), &
                       dataKind_opt=gsv_getDataKind(statevector),        &
                       mpi_local_opt=.true., varNames_opt=varNames )
    deallocate(varNames)

    ! Process components in opposite order as forward calculation
    bmat_loop: do bmatIndex = numBmat, 1, -1

      if ( .not. bmatActive(bmatIndex) ) cycle bmat_loop

      subVector => cvm_getSubVector( controlVector, bmatLabelList(bmatIndex) )
      subVector(:) = 0.0d0

      ! Adjoint of converting statevector from 3D to 4D
      call gsv_copy( statevector, statevector_temp )
      if ( bmatIs3dList(bmatIndex) ) call gsv_3dto4dAdj( statevector_temp )

      select case( trim(bmatTypeList(bmatIndex)) )
      case ('ENS')

        !- 2.1 Flow-dependent Ensemble-Based
        call utl_tmg_start(61,'----B_ENS_AD')

        call ben_bsqrtad( bmatInstanceID(bmatIndex), statevector_temp, &  ! IN
                          subVector,                                   &  ! OUT
                          useFSOFcst_opt, stateVectorRef_opt )            ! IN
        call utl_tmg_stop(61)

      case ('DIFF')

        !- 2.2 Covariances modelled using a diffusion operator.
        call utl_tmg_start(67,'----B_DIFF_AD')
        call bdiff_bsqrtad( statevector_temp, & ! IN
                            subVector )         ! OUT
        call utl_tmg_stop(67)
 
      case ('CHM')

        !- 2.3  Static (Time-Mean Homogeneous and Isotropic) covariances for constituents
        call utl_tmg_start(69,'----B_CHM_AD')
        if ( globalGrid ) then
          call bchm_bsqrtad( statevector_temp, &  ! IN
                             subVector,        &  ! OUT
                             stateVectorRef_opt ) ! IN
        end if
        call utl_tmg_stop(69)

      case ('HI')

        !- 2.4 Time-Mean Homogeneous and Isotropic...
        call utl_tmg_start(53,'----B_HI_AD')
        if ( globalGrid ) then
          call bhi_bsqrtad( statevector_temp, &  ! IN
                            subVector,        &  ! OUT
                            stateVectorRef_opt ) ! IN
        else
          call lbhi_bSqrtAdj( statevector_temp, &  ! IN
                              subVector,        &  ! OUT
                              stateVectorRef_opt ) ! IN
        end if
        call utl_tmg_stop(53)

      end select

    end do bmat_loop

    call gsv_deallocate( statevector_temp )

    call utl_tmg_stop(50)

  end subroutine bmat_sqrtBT

  !--------------------------------------------------------------------------
  ! bmat_finalize
  !--------------------------------------------------------------------------
  subroutine bmat_finalize()   
    !
    !:Purpose: To release memory used by B matrices.
    !
    implicit none 

    call bhi_finalize()
    call ben_finalize()
    call bchm_finalize()
    call lbhi_finalize()
    call bdiff_finalize()

  end subroutine bmat_finalize

  !--------------------------------------------------------------------------
  ! bmat_reduceToMPILocal
  !--------------------------------------------------------------------------
  subroutine bmat_reduceToMPILocal(cv_mpilocal,cv_mpiglobal)    
    !
    !:Purpose: To distribute MPI_global control vector from task 0 to all tasks
    !          where the arguments are real(8)'s
    !
    implicit none

    ! Arguments:
    real(8), intent(inout) :: cv_mpilocal(:)
    real(8), intent(in)    :: cv_mpiglobal(:)

    ! Locals:
    integer :: bmatIndex
    real(8), pointer :: subVector_mpilocal(:), subVector_mpiglobal(:)
    real(8), target  :: dummyVector(1)

    bmat_loop: do bmatIndex = 1, numBmat

      if ( .not. bmatActive(bmatIndex) ) cycle bmat_loop

      subVector_mpilocal => cvm_getSubVector( cv_mpilocal, bmatLabelList(bmatIndex) )
      if ( mmpi_myid == 0 ) then
         subVector_mpiglobal => cvm_getSubVector_mpiglobal( cv_mpiglobal, bmatLabelList(bmatIndex) )
      else
         subVector_mpiglobal => dummyVector
      end if

      select case( trim(bmatTypeList(bmatIndex)) )
      case ('HI')

        !- 2.1 Time-Mean Homogeneous and Isotropic...
        if ( globalGrid ) then
          call bhi_reduceToMPILocal( subVector_mpilocal,subVector_mpiglobal )
        else
          call lbhi_reduceToMPILocal( subVector_mpilocal,subVector_mpiglobal )
        end if

      case ('CHM')

        !- 2.2  Static (Time-Mean Homogeneous and Isotropic) covariances for constituents
        if ( globalGrid ) then
          call bchm_reduceToMPILocal( subVector_mpilocal, subVector_mpiglobal )
        end if

      case ('DIFF')

        !- 2.3 Covariances modelled using a diffusion operator.
        call bdiff_reduceToMPILocal( subVector_mpilocal, subVector_mpiglobal )

      case ('ENS')

        !- 2.4 Flow-dependent Ensemble-Based
        call ben_reduceToMPILocal(subVector_mpilocal, subVector_mpiglobal, bmatInstanceID(bmatIndex))

      end select

    end do bmat_loop

  end subroutine bmat_reduceToMPILocal

  !--------------------------------------------------------------------------
  ! bmat_reduceToMPILocal_r4
  !--------------------------------------------------------------------------
  subroutine bmat_reduceToMPILocal_r4(cv_mpilocal,cv_mpiglobal)
    !
    !:Purpose: To distribute MPI_global control vector from task 0 to all tasks
    !          where the arguments are real(4)'s.
    !
    implicit none

    ! Arguments:
    real(4), intent(inout) :: cv_mpilocal(:)
    real(4), intent(in)    :: cv_mpiglobal(:)

    ! Locals:
    integer :: bmatIndex
    real(4), pointer :: subVector_mpilocal(:), subVector_mpiglobal(:)
    real(4), target  :: dummyVector_r4(1)

    bmat_loop: do bmatIndex = 1, numBmat

      if ( .not. bmatActive(bmatIndex) ) cycle bmat_loop

      subVector_mpilocal => cvm_getSubVector_r4( cv_mpilocal, bmatLabelList(bmatIndex) )
      if ( mmpi_myid == 0 ) then
         subVector_mpiglobal => cvm_getSubVector_mpiglobal_r4( cv_mpiglobal, bmatLabelList(bmatIndex) )
      else
         subVector_mpiglobal => dummyVector_r4
      end if

      select case( trim(bmatTypeList(bmatIndex)) )
      case ('HI')

        !- 2.1 Time-Mean Homogeneous and Isotropic...
        if ( globalGrid ) then
          call bhi_reduceToMPILocal_r4( subVector_mpilocal,subVector_mpiglobal )
        else
          call lbhi_reduceToMPILocal_r4( subVector_mpilocal,subVector_mpiglobal )
        end if

      case ('CHM')

        !- 2.2  Static (Time-Mean Homogeneous and Isotropic) covariances for constituents
        if ( globalGrid ) then
          call bchm_reduceToMPILocal_r4( subVector_mpilocal, subVector_mpiglobal )
        end if

      case ('DIFF')

        !- 2.3 Covariances modelled using a diffusion operator.
        call utl_abort('bmat_reduceToMPILocal_r4: not yet implemented for bMatrixDiff')
        !call bdiff_reduceToMPILocal_r4( subVector_mpilocal, subVector_mpiglobal )

      case ('ENS')

        !- 2.4 Flow-dependent Ensemble-Based
        call ben_reduceToMPILocal_r4(subVector_mpilocal, subVector_mpiglobal, bmatInstanceID(bmatIndex))

      end select

    end do bmat_loop

  end subroutine bmat_reduceToMPILocal_r4

  !--------------------------------------------------------------------------
  ! bmat_expandToMPIGlobal
  !--------------------------------------------------------------------------
  subroutine bmat_expandToMPIGlobal(cv_mpilocal,cv_mpiglobal)
    !
    !:Purpose: To gather control vector from all tasks to task 0 where the
    !          arguments are real(8)'s.
    !
    implicit none

    ! Arguments:
    real(8), intent(in)    :: cv_mpilocal(:)
    real(8), intent(inout) :: cv_mpiglobal(:)

    ! Locals:
    integer :: bmatIndex
    real(8), pointer :: subVector_mpilocal(:), subVector_mpiglobal(:)
    real(8), target  :: dummyVector(1)

    bmat_loop: do bmatIndex = 1, numBmat

      if ( .not. bmatActive(bmatIndex) ) cycle bmat_loop

      subVector_mpilocal => cvm_getSubVector( cv_mpilocal, bmatLabelList(bmatIndex) )
      if ( mmpi_myid == 0 ) then
         subVector_mpiglobal => cvm_getSubVector_mpiglobal( cv_mpiglobal, bmatLabelList(bmatIndex) )
      else
         subVector_mpiglobal => dummyVector
      end if

      select case( trim(bmatTypeList(bmatIndex)) )
      case ('HI')

        !- 2.1 Time-Mean Homogeneous and Isotropic...
        if ( globalGrid ) then
          call bhi_expandToMPIGlobal( subVector_mpilocal,subVector_mpiglobal )
        else
          call lbhi_expandToMPIGlobal( subVector_mpilocal,subVector_mpiglobal )
        end if

      case ('CHM')

        !- 2.2  Static (Time-Mean Homogeneous and Isotropic) covariances for constituents
        if ( globalGrid ) then
          call bchm_expandToMPIGlobal( subVector_mpilocal, subVector_mpiglobal )
        end if

      case ('DIFF')

        !- 2.3 Covariances modelled using a diffusion operator.
        !call bdiff_expandToMPIGlobal( subVector_mpilocal, subVector_mpiglobal )

      case ('ENS')

        !- 2.4 Flow-dependent Ensemble-Based
        call ben_expandToMPIGlobal(subVector_mpilocal, subVector_mpiglobal, bmatInstanceID(bmatIndex) )

      end select

    end do bmat_loop

  end subroutine bmat_expandToMPIGlobal

  !--------------------------------------------------------------------------
  ! bmat_expandToMPIGlobal_r4
  !--------------------------------------------------------------------------
  subroutine bmat_expandToMPIGlobal_r4(cv_mpilocal,cv_mpiglobal)
    !
    !:Purpose: To gather control vector from all tasks to task 0 where the
    !          arguments are real(4)'s.
    !
    implicit none

    ! Arguments:
    real(4), intent(in)    :: cv_mpilocal(:)
    real(4), intent(inout) :: cv_mpiglobal(:)

    ! Locals:
    integer :: bmatIndex
    real(4), pointer :: subVector_mpilocal(:), subVector_mpiglobal(:)
    real(4), target  :: dummyVector_r4(1)

    bmat_loop: do bmatIndex = 1, numBmat

      if ( .not. bmatActive(bmatIndex) ) cycle bmat_loop

      subVector_mpilocal => cvm_getSubVector_r4( cv_mpilocal, bmatLabelList(bmatIndex) )
      if ( mmpi_myid == 0 ) then
         subVector_mpiglobal => cvm_getSubVector_mpiglobal_r4( cv_mpiglobal, bmatLabelList(bmatIndex) )
      else
         subVector_mpiglobal => dummyVector_r4
      end if

      select case( trim(bmatTypeList(bmatIndex)) )
      case ('HI')

        !- 2.1 Time-Mean Homogeneous and Isotropic...
        if ( globalGrid ) then
          call bhi_expandToMPIGlobal_r4( subVector_mpilocal,subVector_mpiglobal )
        else
          call lbhi_expandToMPIGlobal_r4( subVector_mpilocal,subVector_mpiglobal )
        end if

      case ('CHM')

        !- 2.2  Static (Time-Mean Homogeneous and Isotropic) covariances for constituents
        if ( globalGrid ) then
          call bchm_expandToMPIGlobal_r4( subVector_mpilocal, subVector_mpiglobal )
        end if

      case ('DIFF')

        !- 2.3 Covariances modelled using a diffusion operator.
        !call bdiff_expandToMPIGlobal_r4( subVector_mpilocal, subVector_mpiglobal )

      case ('ENS')

        !- 2.4 Flow-dependent Ensemble-Based
        call ben_expandToMPIGlobal_r4(subVector_mpilocal, subVector_mpiglobal, bmatInstanceID(bmatIndex) )

      end select

    end do bmat_loop

  end subroutine bmat_expandToMPIGlobal_r4

end module bMatrix_mod
