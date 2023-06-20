program midas_adjointTest
  !
  !:Purpose: Main program for adjoint test applications, i.e., testing if the identity
  !          :math:`<x,L(y)> = <L^T(x),y>` is respected.
  !
  !:Algorithm: Initialize a random grid-point space statevector :math:`x` and a random
  !            control vector :math:`y`.  Compute :math:`L(y)` by the action of a linear
  !            operator and compute the inner product :math:`<x,L(y)>`.  Then compute
  !            :math:`L^T(x)` through the action of the corresponding adjoint operator
  !            and the inner product :math:`<L^T(x),y>`.  Verify if both inner products
  !            are equal.
  !
  !            --
  !
  !:File I/O: The required input files vary according to the application (see options below).
  !
  !================================================= ==============================================================
  ! Input and Output Files (square root covariances)  Description of file
  !================================================= ==============================================================
  ! ``flnml``                                         In - Main namelist file with parameters user may modify
  ! ``trlm_01``                                       In - Background state (a.k.a. trial) file
  ! ``analysisgrid``                                  In - File defining grid for computing the analysis increment
  ! ``bgcov``                                         In - Static (i.e. NMC) B matrix file for NWP fields
  ! ``ensemble/$YYYYMMDDHH_006_$NNNN``                In - Ensemble member files defining ensemble B matrix
  ! ``innerProd.txt``                                 Out - Results of inner products and their difference
  !================================================= ==============================================================
  !
  !:Synopsis: There are several flavors of tests that can be run. The choice is configured
  !           in the namelist block ``NAMADT`` with the variable ``test`` (see below).
  !           All tests follow this general synopsis:
  !
  !           - initialize with gaussian noise (``rng_gaussian()``) a statevector
  !             :math:`x` and a control vector :math:`y`.
  !
  !           - apply the corresponding square root operator on :math:`y` to obtain :math:`Ly`
  !
  !           - compute the inner product :math:`<x ,Ly>`
  !
  !           - apply the corresponding square root adjoint operator on :math:`x` to
  !             obtain :math:`L^T(x)`
  !
  !           - compute the inner product :math:`<L^Tx,y>`
  !
  !           - print :math:`<x,L(y)>`, :math:`<L^T(x),y>` and their relative difference
  !
  !:Options: `List of namelist blocks <../namelists_in_each_program.html#adjointtest>`_
  !          that can affect the ``adjointTest`` program.
  !
  !          * As described in the algorithm section, four different tests are implemented.
  !            The test that is conducted is chosen through the namelist variable
  !            `&NAMADT Test` that can be either
  !
  !             * **'Bhi'** : test the adjoint of the square root of the homogeneous and
  !               isotropic covariance matrix.  The namelist block `&NAMBHI` will
  !               configure the covariance matrix properties.
  !
  !             * **'Bens'** : test the adjoint of the square root of the ensemble covariance
  !               matrix.  The namelist block `&NAMBEN` will configure the covariance
  !               matrix properties.
  !
  !             * **'loc'** : test the adjoint of the square root of localized covariance
  !               matrix.
  !
  !             * **'advEns'** : test the adjoint of the advection operator on ensemble
  !
  !             * **'advGSV'** : test the adjoint of the advection operator on a single
  !               statevector
  !
  use version_mod
  use codePrecision_mod
  use ramDisk_mod
  use utilities_mod
  use midasMpi_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use gridStateVector_mod
  use gridVariableTransforms_mod
  use bMatrixHI_mod
  use bMatrixEnsemble_mod
  use randomNumber_mod
  use advection_mod
  use ensembleStateVector_mod
  use localization_mod
  use lamBmatrixHI_mod
  implicit none

  type(struct_vco),       pointer :: vco_anl  => null()
  type(struct_hco),       pointer :: hco_anl  => null()
  type(struct_hco),       pointer :: hco_core => null()

  type(struct_gsv) :: statevector_x
  type(struct_gsv) :: statevector_y
  type(struct_gsv) :: statevector_Ly
  type(struct_gsv) :: statevector_LTx

  real(8) :: innerProduct1_local, innerProduct2_local, innerProduct1_global, innerProduct2_global

  real(8), allocatable ::  controlVector1(:)
  real(8), allocatable ::  controlVector2(:)

  integer :: get_max_rss, ierr, cvDim, nulnam, fnom, fclos
  character(len=20) ::  test ! adjoint test type ('Bhi','Bens','advEns','advGSV','loc')
  NAMELIST /NAMADT/test

  call ver_printNameAndVersion('adjointTest','Various tests of adjoint codes')

  !
  !- 1.  Settings and module initializations
  !
  write(*,*)
  write(*,*) '> midas-adjointTest: setup - START'

  !- 1.1 mpi
  call mmpi_initialize

  !- 1.2 timings
  call tmg_init(mmpi_myid, 'TMG_INFO')
  call utl_tmg_start(0,'Main')

  !- 1.3 RAM disk usage
  call ram_setup

  !- 1.4 Temporal grid and set dateStamp from env variable
  call tim_setup()
  if (tim_getDateStamp()==0) then
    call utl_abort('midas-adjointTest: dateStamp was not set')
  end if

  !- 1.6 Constants
  if ( mmpi_myid == 0 ) then
    call mpc_printConstants(6)
    call pre_printPrecisions
  end if

  !- 1.8 Variables of the model states
  call gsv_setup

  !- 1.9 Set the horizontal domain
  call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS') ! IN
  if ( hco_anl % global ) then
    hco_core => hco_anl
  else
    !- Iniatilized the core (Non-Exteded) analysis grid
    call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID') ! IN
  end if

  !- 1.10 Initialize the vertical coordinate
  call vco_SetupFromFile( vco_anl,        & ! OUT
                          './analysisgrid') ! IN

  write(*,*)
  write(*,*) '> midas-adjointTest: setup - END'
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !- 1.11 Variable transforms
  call gvt_Setup(hco_anl, hco_core, vco_anl)

  !- 1.12 Test selection
  test = 'Bhi' ! default test 

  if ( .not. utl_isNamelistPresent('NAMADT','./flnml') ) then
    if (mmpi_myid == 0) then
      write(*,*) 'midas-adjointTest: namadt is missing in the namelist. '&
                 //'The default values will be taken.'
    end if
  else
    nulnam=0
    ierr=fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    read(nulnam, nml=namadt, iostat=ierr)
    if(ierr.ne.0) call utl_abort('midas-adjointTest: Error reading namelist')
    if( mmpi_myid == 0 ) write(*,nml=namadt)
    ierr=fclos(nulnam)
  end if

  !
  !- 2.  The tests
  !
 
  write(*,*)
  write(*,*) '> midas-adjointTest: '//test
  if ( test == 'Bhi') then
    !- 2.1 Bhi
    call check_bhi
  else if ( test == 'Bens') then  
    !- 2.2 Bens
    call check_bens
  else if ( test == 'advEns') then
    !- 2.3 AdvectionENS
    call check_advectionENS
  else if ( test == 'advGSV') then 
    !- 2.4 AdvectionGSV
    call check_advectionGSV
  else if ( test == 'loc') then
    !- 2.5 Localization
    call check_loc 
  !else if ( test == 'addMem') then
  !  !- 2.6 Localization
  !  call check_addmem
  else
    call utl_abort('midas-adjointTest: inexistant test label ('//test//')')
  end if

  !
  !- 3.  Ending
  !
  write(*,*)
  write(*,*) '> midas-adjointTest: Ending'
  call utl_tmg_stop(0)
  call tmg_terminate(mmpi_myid, 'TMG_INFO')

  call rpn_comm_finalize(ierr) 

contains

  !--------------------------------------------------------------------------
  !- check Bhi
  !--------------------------------------------------------------------------
  subroutine check_bhi()
    implicit none

    integer :: seed, kIndex, stepIndex, latIndex, lonIndex, cvIndex
    real(8), pointer :: field4d_r8(:,:,:,:), field3d_Ly_r8(:,:,:), field3d_x_r8(:,:,:)

    call gvt_setupRefFromTrialFiles('HU')
    if (hco_anl%global) then
      call bhi_Setup( hco_anl, vco_anl, & ! IN
                      cvdim )             ! OUT
    else
      call lbhi_Setup( hco_anl, hco_core, vco_anl, & ! IN
                      cvdim )                        ! OUT
    end if

    call gsv_allocate(statevector_x  , tim_nstepobsinc, hco_anl, vco_anl, &
                      mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)

    call gsv_allocate(statevector_Ly , tim_nstepobsinc, hco_anl, vco_anl, &
                      mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)

    allocate ( controlVector1(cvDim) )
    
    ! x
    !statevector_x%gd3d_r8(:,:,:) = 13.3d0
    call gsv_getField(statevector_x,field4d_r8)
    seed=1
    call rng_setup(abs(seed+mmpi_myid))
    do kIndex = statevector_x%mykBeg, statevector_x%mykEnd
      do stepIndex = 1, statevector_x%numStep
        do latIndex = statevector_x%myLatBeg, statevector_x%myLatEnd
          do lonIndex = statevector_x%myLonBeg, statevector_x%myLonEnd
            field4d_r8(lonIndex,latIndex,kIndex,stepIndex) = rng_gaussian()
          end do
        end do
      end do
    end do

    ! y
    !controlVector1(:) = 2.4d0
    do cvIndex = 1, cvDim
      controlVector1(cvIndex) = rng_gaussian()
    end do

    ! y
    !controlVector1 = 0.d0
    !if (hco_anl%global) then
    !  call bhi_bSqrtAd( statevector_x, &  ! IN
    !                    controlVector1)  ! OUT
    !else
    !  call lbhi_bSqrtAdj( statevector_x, &  ! IN
    !                      controlVector1)  ! OUT
    !end if

    ! Ly
    if (hco_anl%global) then
      call bhi_bSqrt( controlVector1, & ! IN
                      statevector_Ly )  ! OUT
    else
      call lbhi_bSqrt( controlVector1, & ! IN
                       statevector_Ly )  ! OUT
    end if


    ! <x ,L(y)>
    innerProduct1_local = 0.d0
    call gsv_getField(statevector_Ly,field3d_Ly_r8)
    call gsv_getField(statevector_x, field3d_x_r8 )
    call euclid(innerProduct1_local, &
         field3d_x_r8(statevector_Ly%myLonBeg:statevector_Ly%myLonEnd,statevector_Ly%myLatBeg:statevector_Ly%myLatEnd,:), &
         field3d_Ly_r8(statevector_Ly%myLonBeg:statevector_Ly%myLonEnd,statevector_Ly%myLatBeg:statevector_Ly%myLatEnd,:), &
         statevector_Ly%myLonBeg, statevector_Ly%myLonEnd, statevector_Ly%myLatBeg, statevector_Ly%myLatEnd, statevector_Ly%nk, 1)
    write(*,*) "<x     ,L(y)> local = ",innerProduct1_local
    call rpn_comm_allreduce(innerProduct1_local,innerProduct1_global,1,"mpi_double_precision","mpi_sum","GRID",ierr)
    write(*,*) "<x     ,L(y)> global= ",innerProduct1_global
    
    ! L_T(x)
    allocate ( controlVector2(cvDim) )
    controlVector2 = 0.d0
    if (hco_anl%global) then
      call bhi_bSqrtAd( statevector_x, &  ! IN
                        controlVector2 ) ! OUT
    else
      call lbhi_bSqrtAdj( statevector_x, &  ! IN
                         controlVector2 ) ! OUT
    end if

    ! <L_T(x),y>
    innerProduct2_local = 0.d0
    call euclid(innerProduct2_local, controlVector2, controlVector1, 1, cvDim, 1, 1, 1, 1)
    print*,"<Lt(x) ,y   > local = ",innerProduct2_local
    call rpn_comm_allreduce(innerProduct2_local,innerProduct2_global,1,"mpi_double_precision","mpi_sum","GRID",ierr)
    write(*,*) "<Lt(x) ,y   > global= ",innerProduct2_global
    
    ! Results
    call checkAndOutputInnerProd
    
    deallocate(controlVector2)
    deallocate(controlVector1)

    call gsv_deallocate(statevector_x)
    call gsv_deallocate(statevector_Ly)

  end subroutine check_bhi

  !--------------------------------------------------------------------------
  !- check Bens
  !--------------------------------------------------------------------------
  subroutine check_bens()
    implicit none

    integer :: seed, kIndex, stepIndex, latIndex, lonIndex

    integer, allocatable :: cvDimPerInstance(:)
    real(8), pointer :: field4d_Ly_r8(:,:,:,:), field4d_x_r8(:,:,:,:)

    call gvt_setupRefFromTrialFiles('HU')
    call ben_Setup( hco_anl, hco_core, vco_anl, & ! IN
                    cvdimPerInstance )  ! OUT

    cvDim = cvdimPerInstance(1)
    call gsv_allocate(statevector_x  , tim_nstepobsinc, hco_anl, vco_anl, &
                      mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)
    call gsv_allocate(statevector_Ly , tim_nstepobsinc, hco_anl, vco_anl, &
                      mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)

    allocate ( controlVector1(cvdim) )

    ! x
    seed=1
    call rng_setup(abs(seed+mmpi_myid))
    call gsv_getField(statevector_x,field4d_x_r8)
    do kIndex = statevector_x%mykBeg, statevector_x%mykEnd
      do stepIndex = 1, statevector_x%numStep
        do latIndex = statevector_x%myLatBeg, statevector_x%myLatEnd
          do lonIndex = statevector_x%myLonBeg, statevector_x%myLonEnd
            field4d_x_r8(lonIndex,latIndex,kIndex,stepIndex) = rng_gaussian()
          end do
        end do
      end do
    end do

    ! y
    controlVector1 = 0.d0
    call ben_bSqrtAd( 1, statevector_x, &  ! IN
                      controlVector1)  ! OUT

    ! Ly
    call ben_bSqrt( 1, controlVector1, & ! IN
                    statevector_Ly )  ! OUT

    ! <x ,L(y)>
    innerProduct1_local = 0.d0
    call gsv_getField(statevector_Ly,field4d_Ly_r8)
    call gsv_getField(statevector_x, field4d_x_r8 )
    call euclid(innerProduct1_local, &
         field4d_x_r8(statevector_Ly%myLonBeg:statevector_Ly%myLonEnd,statevector_Ly%myLatBeg:statevector_Ly%myLatEnd,:,:), &
         field4d_Ly_r8(statevector_Ly%myLonBeg:statevector_Ly%myLonEnd,statevector_Ly%myLatBeg:statevector_Ly%myLatEnd,:,:), &
         statevector_Ly%myLonBeg, statevector_Ly%myLonEnd, statevector_Ly%myLatBeg, statevector_Ly%myLatEnd, statevector_Ly%nk, statevector_Ly%numStep)
    write(*,*) "<x     ,L(y)> local = ",innerProduct1_local
    call rpn_comm_allreduce(innerProduct1_local,innerProduct1_global,1,"mpi_double_precision","mpi_sum","GRID",ierr)
    write(*,*) "<x     ,L(y)> global= ",innerProduct1_global
    
    ! L_T(x)
    allocate ( controlVector2(cvDim) )
    controlVector2 = 0.d0
    call ben_bSqrtAd( 1, statevector_x, &  ! IN
                      controlVector2 )  ! OUT
    
    ! <L_T(x),y>
    innerProduct2_local = 0.d0
    call euclid(innerProduct2_local, controlVector2, controlVector1, 1, cvDim, 1, 1, 1, 1)
    print*,"<Lt(x) ,y   > local = ",innerProduct2_local
    call rpn_comm_allreduce(innerProduct2_local,innerProduct2_global,1,"mpi_double_precision","mpi_sum","GRID",ierr)
    write(*,*) "<Lt(x) ,y   > global= ",innerProduct2_global

    ! Results
    call checkAndOutputInnerProd
    
    deallocate(controlVector2)
    deallocate(controlVector1)

    call gsv_deallocate(statevector_x)
    call gsv_deallocate(statevector_Ly)

  end subroutine check_bens

  !--------------------------------------------------------------------------
  !- check Loc
  !--------------------------------------------------------------------------
  subroutine check_loc()
    implicit none

    integer :: seed, kIndex, stepIndex, latIndex, lonIndex, dateStamp
    integer :: numStepAmplitude, amp3dStepIndex, memberIndex

    type(struct_ens) :: ensAmplitude_x
    type(struct_ens) :: ensAmplitude_Ly

    type(struct_loc), pointer :: loc => null()

    real(8), pointer     :: ens_oneLev(:,:,:,:)
    real(8), pointer :: field4d_Ly_r8(:,:,:,:), field4d_x_r8(:,:,:,:)

    character(len=4), parameter  :: varNameALFAatm(1) = (/ 'ALFA' /)
    character(len=4), parameter  :: varNameALFAsfc(1) = (/ 'ALFS' /)
    character(len=4)             :: varNameALFA(1)

    integer,allocatable :: dateStampList(:)

    integer, allocatable :: cvDimPerInstance(:)

    allocate(dateStampList(tim_nstepobsinc))
    call tim_getstamplist(dateStampList,tim_nstepobsinc,tim_getDatestamp())

    dateStamp = tim_getDatestamp()
    write(*,*) 'check_loc: tim_getDatestamp = ', dateStamp
    write(*,*) 'check_loc: dateStampList = ', dateStampList(:)

    call ben_Setup( hco_anl, hco_core, vco_anl, & ! IN
                    cvDimPerInstance )             ! OUT

    cvDim = cvDimPerInstance(1)

    loc => ben_getLoc(1)
    if ( cvDim /= loc%cvDim ) then
      call utl_abort('check_loc: cvDim /= loc%cvDim')
    end if
    numStepAmplitude = ben_getNumStepAmplitudeAssimWindow()
    if ( numStepAmplitude /= 1 ) then
      call utl_abort('check_loc: not yet adapted for localization advection')
    end if
    amp3dStepIndex   = ben_getAmp3dStepIndexAssimWindow()

    if (loc%vco%Vcode == 5002 .or. loc%vco%Vcode == 5005) then
      varNameALFA(:) = varNameALFAatm(:)
    else ! vco_anl%Vcode == -1
      varNameALFA(:) = varNameALFAsfc(:)
    end if

    call gsv_allocate(statevector_x  , numStepAmplitude, loc%hco, loc%vco, &
                      mpi_local_opt=.true., varNames_opt=varNameALFA, dataKind_opt=8)
    call gsv_allocate(statevector_Ly , numStepAmplitude, loc%hco, loc%vco, &
                      mpi_local_opt=.true., varNames_opt=varNameALFA, dataKind_opt=8)

    call ens_allocate(ensAmplitude_x, loc%nEnsOverDimension, numStepAmplitude, loc%hco, loc%vco, &
                      datestampList=dateStampList, varNames_opt=varNameALFA, dataKind_opt=8)    
    call ens_allocate(ensAmplitude_Ly, loc%nEnsOverDimension, numStepAmplitude, loc%hco, loc%vco, &
                      datestampList=dateStampList, varNames_opt=varNameALFA, dataKind_opt=8)

    allocate ( controlVector1(cvDim) )

    ! x
    seed=1
    call rng_setup(abs(seed+mmpi_myid))
    do kIndex = 1, ens_getNumK(ensAmplitude_x)
      ens_oneLev => ens_getOneLev_r8(ensAmplitude_x,kIndex)
      do memberIndex = 1, loc%nEnsOverDimension
        do stepIndex = 1,numStepAmplitude
          do latIndex = statevector_x%myLatBeg, statevector_x%myLatEnd
            do lonIndex = statevector_x%myLonBeg, statevector_x%myLonEnd
              ens_oneLev(memberIndex,stepIndex,lonIndex,latIndex) = rng_gaussian()
            end do
          end do
        end do
      end do
    end do

    ! y
    controlVector1 = 0.d0
    call loc_LsqrtAd(loc,             & ! IN
                     ensAmplitude_x,  & ! IN
                     controlVector1,  & ! OUT
                     amp3dStepIndex)    ! IN

    ! Ly
    call loc_Lsqrt  (loc,           & ! IN
                     controlVector1, & ! IN
                     ensAmplitude_Ly,  & ! OUT
                     amp3dStepIndex)  ! IN

    ! <x ,L(y)>
    innerProduct1_local = 0.d0
    call gsv_getField(statevector_Ly,field4d_Ly_r8)
    call gsv_getField(statevector_x, field4d_x_r8 )
    do memberIndex = 1, loc%nEnsOverDimension
      call ens_copyMember(ensAmplitude_x , statevector_x , memberIndex)
      call ens_copyMember(ensAmplitude_Ly, statevector_Ly, memberIndex)
      call euclid(innerProduct1_local, &
           field4d_x_r8(statevector_Ly%myLonBeg:statevector_Ly%myLonEnd,statevector_Ly%myLatBeg:statevector_Ly%myLatEnd,:,:), &
           field4d_Ly_r8(statevector_Ly%myLonBeg:statevector_Ly%myLonEnd,statevector_Ly%myLatBeg:statevector_Ly%myLatEnd,:,:), &
           statevector_Ly%myLonBeg, statevector_Ly%myLonEnd, statevector_Ly%myLatBeg, statevector_Ly%myLatEnd, statevector_Ly%nk, statevector_Ly%numStep)
    end do
    write(*,*) "<x     ,L(y)> local = ",innerProduct1_local
    call rpn_comm_allreduce(innerProduct1_local,innerProduct1_global,1,"mpi_double_precision","mpi_sum","GRID",ierr)
    write(*,*) "<x     ,L(y)> global= ",innerProduct1_global
    
    ! L_T(x)
    allocate ( controlVector2(cvDim) )
    controlVector2 = 0.d0
    call loc_LsqrtAd(loc,             & ! IN
                     ensAmplitude_x,  & ! IN
                     controlVector2,  & ! OUT
                     amp3dStepIndex)    ! IN

    ! <L_T(x),y>
    innerProduct2_local = 0.d0
    call euclid(innerProduct2_local, controlVector2, controlVector1, 1, cvDim, 1, 1, 1, 1)
    print*,"<Lt(x) ,y   > local = ",innerProduct2_local
    call rpn_comm_allreduce(innerProduct2_local,innerProduct2_global,1,"mpi_double_precision","mpi_sum","GRID",ierr)
    write(*,*) "<Lt(x) ,y   > global= ",innerProduct2_global
    
    ! Results
    call checkAndOutputInnerProd
    
    deallocate(controlVector2)
    deallocate(controlVector1)

    call gsv_deallocate(statevector_x)
    call gsv_deallocate(statevector_Ly)

    call ens_deallocate(ensAmplitude_x)
    call ens_deallocate(ensAmplitude_Ly)

  end subroutine check_loc

!!$  !--------------------------------------------------------------------------
!!$  !- check addMem
!!$  !--------------------------------------------------------------------------
!!$  subroutine check_addmem()
!!$    implicit none
!!$
!!$    integer :: seed, kIndex, stepIndex, latIndex, lonIndex
!!$    integer :: numStepAmplitude, amp3dStepIndex, memberIndex
!!$
!!$    type(struct_ens) :: ensAmplitude_LTx
!!$    type(struct_ens) :: ensAmplitude_y
!!$
!!$    type(struct_loc), pointer :: loc => null()
!!$
!!$    real(8), pointer     :: ens_oneLev(:,:,:,:)
!!$
!!$    character(len=4), parameter  :: varNameALFAatm(1) = (/ 'ALFA' /)
!!$    character(len=4), parameter  :: varNameALFAsfc(1) = (/ 'ALFS' /)
!!$    character(len=4)             :: varNameALFA(1)
!!$
!!$    integer,allocatable :: dateStampList(:)
!!$
!!$    allocate(dateStampList(tim_nstepobsinc))
!!$    call tim_getstamplist(dateStampList,tim_nstepobsinc,tim_getDatestamp())
!!$
!!$    write(*,*) 'JFC tim_getDatestamp = ', tim_getDatestamp()
!!$    write(*,*) 'JFC dateStampList = ', dateStampList(:)
!!$
!!$    call ben_Setup( hco_anl, hco_core, vco_anl, & ! IN
!!$                    cvdim )             ! OUT
!!$
!!$    write(*,*) 'JFC ben_Setup done '
!!$    
!!$    
!!$    loc => ben_getLoc(1)
!!$    if ( cvDim /= loc%cvDim ) then
!!$      call utl_abort('check_loc: cvDim /= loc%cvDim')
!!$    end if
!!$    numStepAmplitude = ben_getNumStepAmplitudeAssimWindow()
!!$    if ( numStepAmplitude /= 1 ) then
!!$      call utl_abort('check_loc: not yet adapted for localization advection')
!!$    end if
!!$    amp3dStepIndex   = ben_getAmp3dStepIndexAssimWindow()
!!$
!!$    if (loc%vco%Vcode == 5002 .or. loc%vco%Vcode == 5005) then
!!$      varNameALFA(:) = varNameALFAatm(:)
!!$    else ! vco_anl%Vcode == -1
!!$      varNameALFA(:) = varNameALFAsfc(:)
!!$    end if
!!$
!!$    call gsv_allocate(statevector_x  , tim_nstepobsinc, hco_anl, vco_anl, &
!!$                      mpi_local_opt=.true., &
!!$                      allocHeight_opt=.false., allocPressure_opt=.false.)
!!$    call gsv_allocate(statevector_Ly , tim_nstepobsinc, hco_anl, vco_anl, &
!!$                      mpi_local_opt=.true., &
!!$                      allocHeight_opt=.false., allocPressure_opt=.false.)
!!$
!!$    call gsv_allocate(statevector_LTx  , numStepAmplitude, loc%hco, loc%vco, &
!!$                      mpi_local_opt=.true., varNames_opt=varNameALFA, dataKind_opt=8)
!!$    call gsv_allocate(statevector_y , numStepAmplitude, loc%hco, loc%vco, &
!!$                      mpi_local_opt=.true., varNames_opt=varNameALFA, dataKind_opt=8)
!!$
!!$    call ens_allocate(ensAmplitude_LTx, loc%nEnsOverDimension, numStepAmplitude, loc%hco, loc%vco, &
!!$                        datestampList=dateStampList, varNames_opt=varNameALFA, dataKind_opt=8)    
!!$    call ens_allocate(ensAmplitude_y, loc%nEnsOverDimension, numStepAmplitude, loc%hco, loc%vco, &
!!$                        datestampList=dateStampList, varNames_opt=varNameALFA, dataKind_opt=8)
!!$
!!$    ! x
!!$    seed=1
!!$    call rng_setup(abs(seed+mmpi_myid))
!!$    do kIndex = statevector_x%mykBeg, statevector_x%mykEnd
!!$      do stepIndex = 1, statevector_x%numStep
!!$        do latIndex = statevector_x%myLatBeg, statevector_x%myLatEnd
!!$          do lonIndex = statevector_x%myLonBeg, statevector_x%myLonEnd
!!$            statevector_x%gd_r8(lonIndex,latIndex,kIndex,stepIndex) = rng_gaussian()
!!$          end do
!!$        end do
!!$      end do
!!$    end do
!!$
!!$    ! y
!!$    do kIndex = 1, ens_getNumK(ensAmplitude_y)
!!$      ens_oneLev => ens_getOneLev_r8(ensAmplitude_y,kIndex)
!!$      do memberIndex = 1, loc%nEnsOverDimension
!!$        do stepIndex = 1,tim_nstepobsinc
!!$          do latIndex = statevector_x%myLatBeg, statevector_x%myLatEnd
!!$            do lonIndex = statevector_x%myLonBeg, statevector_x%myLonEnd
!!$              ens_oneLev(memberIndex,stepIndex,lonIndex,latIndex) = rng_gaussian()
!!$            end do
!!$          end do
!!$        end do
!!$      end do
!!$    end do
!!$
!!$    ! Ly
!!$    call addEnsMember( ensAmplitude_y,    & ! IN
!!$                       statevector_Ly,    & ! OUT 
!!$                       1, .false. )         ! IN
!!$
!!$    ! <x ,L(y)>
!!$    innerProduct1_local = 0.d0
!!$    call euclid(innerProduct1_local, &
!!$         statevector_x %gd_r8(statevector_Ly%myLonBeg:statevector_Ly%myLonEnd,statevector_Ly%myLatBeg:statevector_Ly%myLatEnd,:,:), &
!!$         statevector_Ly%gd_r8(statevector_Ly%myLonBeg:statevector_Ly%myLonEnd,statevector_Ly%myLatBeg:statevector_Ly%myLatEnd,:,:), &
!!$         statevector_Ly%myLonBeg, statevector_Ly%myLonEnd, statevector_Ly%myLatBeg, statevector_Ly%myLatEnd, statevector_Ly%nk, statevector_Ly%numStep)
!!$    write(*,*) "<x     ,L(y)> local = ",innerProduct1_local
!!$    call rpn_comm_allreduce(innerProduct1_local,innerProduct1_global,1,"mpi_double_precision","mpi_sum","GRID",ierr)
!!$    write(*,*) "<x     ,L(y)> global= ",innerProduct1_global
!!$    
!!$    ! L_T(x)
!!$    call addEnsMemberad( statevector_x,   & ! IN 
!!$                         ensAmplitude_LTx,  & ! OUT
!!$                         1, .false. )       ! IN
!!$
!!$    ! <L_T(x),y>
!!$    innerProduct2_local = 0.d0
!!$    do memberIndex = 1, loc%nEnsOverDimension
!!$      call ens_copyMember(ensAmplitude_LTx, statevector_LTx, memberIndex)
!!$      call ens_copyMember(ensAmplitude_y  , statevector_y  , memberIndex)
!!$      call euclid(innerProduct2_local, &
!!$           statevector_LTx %gd_r8(statevector_y%myLonBeg:statevector_y%myLonEnd,statevector_y%myLatBeg:statevector_y%myLatEnd,:,:), &
!!$           statevector_y%gd_r8(statevector_y%myLonBeg:statevector_y%myLonEnd,statevector_y%myLatBeg:statevector_y%myLatEnd,:,:), &
!!$           statevector_y%myLonBeg, statevector_y%myLonEnd, statevector_y%myLatBeg, statevector_y%myLatEnd, statevector_y%nk, statevector_y%numStep)
!!$    end do
!!$    print*,"<Lt(x) ,y   > local = ",innerProduct2_local
!!$    call rpn_comm_allreduce(innerProduct2_local,innerProduct2_global,1,"mpi_double_precision","mpi_sum","GRID",ierr)
!!$    write(*,*) "<Lt(x) ,y   > global= ",innerProduct2_global
!!$    
!!$    ! Results
!!$    call checkAndOutputInnerProd
!!$    
!!$    call gsv_deallocate(statevector_LTx)
!!$    call gsv_deallocate(statevector_y)
!!$    call gsv_deallocate(statevector_x)
!!$    call gsv_deallocate(statevector_Ly)
!!$
!!$    call ens_deallocate(ensAmplitude_LTx)
!!$    call ens_deallocate(ensAmplitude_y)
!!$
!!$  end subroutine check_addmem

  !--------------------------------------------------------------------------
  !- check AdvectionENS
  !--------------------------------------------------------------------------
  subroutine check_advectionENS()
    implicit none

    integer :: seed, kIndex, stepIndex, latIndex, lonIndex, dateStamp

    type(struct_adv)  :: adv_analInc

    type(struct_ens) :: ens_x, ens_Ly, ens_y, ens_LTx

    character(len=32)   :: directionAnlInc
    character(len=4), parameter  :: varNameALFA(1) = (/ 'ALFA' /)

    real(8), pointer     :: ens_oneLev(:,:,:,:)
    real(8), pointer     :: field4d_Ly_r8(:,:,:,:), field4d_x_r8(:,:,:,:)
    real(8), pointer     :: field4d_LTx_r8(:,:,:,:), field4d_y_r8(:,:,:,:)

    real(8) :: delT_hour

    real(8), allocatable :: advectFactor(:)

    integer :: numStepAdvect, numStepReferenceFlow, nEns, memberIndex

    integer,allocatable :: dateStampList(:)

    allocate(dateStampList(tim_nstepobsinc))
    call tim_getstamplist(dateStampList,tim_nstepobsinc,tim_getDatestamp())

    dateStamp = tim_getDatestamp()
    write(*,*) 'check_advectionENS: tim_getDatestamp = ', dateStamp
    write(*,*) 'check_advectionENS: dateStampList = ', dateStampList(:)

    directionAnlInc = 'fromFirstTimeIndex' !'towardFirstTimeIndexInverse'
    delT_hour = 1.0d0 !tim_dstepobsinc
    numStepAdvect             = tim_nstepobsinc
    numStepReferenceFlow      = 7
    allocate(advectFactor(vco_anl%nLev_M))
    advectFactor(:) = 0.75D0
    nEns = 1

    call adv_setup( adv_analInc,                                             & ! OUT
                    directionAnlInc, hco_anl, vco_anl,                       & ! IN
                    numStepAdvect, dateStampList,              & ! IN
                    numStepReferenceFlow, delT_hour, advectFactor,           & ! IN
                    'MMLevsOnly', steeringFlowFilename_opt='ensemble/forecast_for_advection' ) ! IN

    deallocate(advectFactor)
    
    call ens_allocate(ens_x, nEns, numStepAdvect, hco_anl, vco_anl, dateStampList, &
                      hco_core_opt=hco_core, varNames_opt=varNameALFA, dataKind_opt=8)
    call ens_allocate(ens_Ly, nEns, numStepAdvect, hco_anl, vco_anl, dateStampList, &
                      hco_core_opt=hco_core, varNames_opt=varNameALFA, dataKind_opt=8)
    call ens_allocate(ens_y, nEns, numStepAdvect, hco_anl, vco_anl, dateStampList, &
                      hco_core_opt=hco_core, varNames_opt=varNameALFA, dataKind_opt=8)
    call ens_allocate(ens_LTx, nEns, numStepAdvect, hco_anl, vco_anl, dateStampList, &
                      hco_core_opt=hco_core, varNames_opt=varNameALFA, dataKind_opt=8)

    call gsv_allocate(statevector_x  , tim_nstepobsinc, hco_anl, vco_anl, &
                      mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)
    call gsv_allocate(statevector_Ly , tim_nstepobsinc, hco_anl, vco_anl, &
                      mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)
    call gsv_allocate(statevector_y  , tim_nstepobsinc, hco_anl, vco_anl, &
                      mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)
    call gsv_allocate(statevector_LTx , tim_nstepobsinc, hco_anl, vco_anl, &
                      mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)

    ! x
    seed=1
    call rng_setup(abs(seed+mmpi_myid))
    do kIndex = 1, ens_getNumK(ens_x)
      ens_oneLev => ens_getOneLev_r8(ens_x,kIndex)
      do memberIndex = 1, nEns
        do stepIndex = 1,tim_nstepobsinc
          do latIndex = statevector_x%myLatBeg, statevector_x%myLatEnd
            do lonIndex = statevector_x%myLonBeg, statevector_x%myLonEnd
              ens_oneLev(memberIndex,stepIndex,lonIndex,latIndex) = rng_gaussian()
            end do
          end do
        end do
      end do
    end do

    ! y
    call ens_copy(ens_x, & ! IN
                  ens_y)  ! OUT
    call adv_ensemble_ad( ens_y, & ! INOUT
                          adv_analInc, nEns )      ! IN
    ! Ly
    call ens_copy(ens_y, ens_Ly)
    call adv_ensemble_tl( ens_Ly, & ! INOUT
                          adv_analInc, nEns )     ! IN

    ! <x ,L(y)>
    innerProduct1_local = 0.d0
    call gsv_getField(statevector_Ly,field4d_Ly_r8)
    call gsv_getField(statevector_x, field4d_x_r8 )
    do memberIndex = 1, nEns
      call ens_copyMember(ens_x , statevector_x , memberIndex)
      call ens_copyMember(ens_Ly, statevector_Ly, memberIndex)
      call euclid(innerProduct1_local, &
           field4d_x_r8(statevector_Ly%myLonBeg:statevector_Ly%myLonEnd,statevector_Ly%myLatBeg:statevector_Ly%myLatEnd,:,:), &
           field4d_Ly_r8(statevector_Ly%myLonBeg:statevector_Ly%myLonEnd,statevector_Ly%myLatBeg:statevector_Ly%myLatEnd,:,:), &
           statevector_Ly%myLonBeg, statevector_Ly%myLonEnd, statevector_Ly%myLatBeg, statevector_Ly%myLatEnd, statevector_Ly%nk, statevector_Ly%numStep)
    end do
    write(*,*) "<x     ,L(y)> local = ",innerProduct1_local
    call rpn_comm_allreduce(innerProduct1_local,innerProduct1_global,1,"mpi_double_precision","mpi_sum","GRID",ierr)
    write(*,*) "<x     ,L(y)> global= ",innerProduct1_global
    
    ! L_T(x)
    call ens_copy(ens_x, & ! IN
                  ens_LTx)  ! OUT
    call adv_ensemble_ad( ens_LTx, & ! INOUT
                          adv_analInc, nEns )      ! IN
    
    ! <L_T(x),y>
    innerProduct2_local = 0.d0
    call gsv_getField(statevector_LTx,field4d_LTx_r8)
    call gsv_getField(statevector_y,  field4d_y_r8 )
    do memberIndex = 1, nEns
      call ens_copyMember(ens_LTx, statevector_LTx, memberIndex)
      call ens_copyMember(ens_y  , statevector_y  , memberIndex)
      call euclid(innerProduct2_local, &
           field4d_LTx_r8(statevector_y%myLonBeg:statevector_y%myLonEnd,statevector_y%myLatBeg:statevector_y%myLatEnd,:,:), &
           field4d_y_r8(statevector_y%myLonBeg:statevector_y%myLonEnd,statevector_y%myLatBeg:statevector_y%myLatEnd,:,:), &
           statevector_y%myLonBeg, statevector_y%myLonEnd, statevector_y%myLatBeg, statevector_y%myLatEnd, statevector_y%nk, statevector_y%numStep)
    end do
    print*,"<Lt(x) ,y   > local = ",innerProduct2_local
    call rpn_comm_allreduce(innerProduct2_local,innerProduct2_global,1,"mpi_double_precision","mpi_sum","GRID",ierr)
    write(*,*) "<Lt(x) ,y   > global= ",innerProduct2_global
    
    ! Results
    call checkAndOutputInnerProd

    call gsv_deallocate(statevector_x)
    call gsv_deallocate(statevector_Ly)
    call gsv_deallocate(statevector_LTx)
    call gsv_deallocate(statevector_y)

    call ens_deallocate(ens_x)
    call ens_deallocate(ens_Ly)
    call ens_deallocate(ens_LTx)
    call ens_deallocate(ens_y)

  end subroutine check_advectionENS

  !--------------------------------------------------------------------------
  !- check AdvectionGSV
  !--------------------------------------------------------------------------
  subroutine check_advectionGSV()
    implicit none

    integer :: seed, kIndex, stepIndex, latIndex, lonIndex, dateStamp

    type(struct_adv)  :: adv_analInc

    character(len=32)   :: directionAnlInc

    real(8) :: delT_hour

    real(8), allocatable :: advectFactor(:)
    real(8), pointer     :: field4d_x_r8(:,:,:,:), field4d_y_r8(:,:,:,:)
    real(8), pointer     :: field4d_LTx_r8(:,:,:,:), field4d_Ly_r8(:,:,:,:)

    integer :: numStepAdvect, numStepReferenceFlow

    integer,allocatable :: dateStampList(:)

    allocate(dateStampList(tim_nstepobsinc))
    call tim_getstamplist(dateStampList,tim_nstepobsinc,tim_getDatestamp())

    dateStamp = tim_getDatestamp()
    write(*,*) 'check_advectionGSV: tim_getDatestamp = ', dateStamp
    write(*,*) 'check_advectionGSV: dateStampList = ', dateStampList(:)

    directionAnlInc = 'fromFirstTimeIndex' !'towardFirstTimeIndexInverse'
    delT_hour = 1.0d0 !tim_dstepobsinc
    numStepAdvect             = tim_nstepobsinc
    numStepReferenceFlow      = 7
    advectFactor = 0.75D0

    call adv_setup( adv_analInc,                                             & ! OUT
                    directionAnlInc, hco_anl, vco_anl,                       & ! IN
                    numStepAdvect, dateStampList,              & ! IN
                    numStepReferenceFlow, delT_hour, advectFactor,           & ! IN
                    'allLevs', steeringFlowFilename_opt='ensemble/forecast_for_advection' ) ! IN

    call gsv_allocate(statevector_x  , tim_nstepobsinc, hco_anl, vco_anl, &
                      mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)
    call gsv_allocate(statevector_Ly , tim_nstepobsinc, hco_anl, vco_anl, &
                      mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)
    call gsv_allocate(statevector_y  , tim_nstepobsinc, hco_anl, vco_anl, &
                      mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)
    call gsv_allocate(statevector_LTx , tim_nstepobsinc, hco_anl, vco_anl, &
                      mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)

    ! x
    seed=1
    call rng_setup(abs(seed+mmpi_myid))
    call gsv_getField(statevector_x,  field4d_x_r8 )
    do kIndex = statevector_x%mykBeg, statevector_x%mykEnd
      do stepIndex = 1, statevector_x%numStep
        do latIndex = statevector_x%myLatBeg, statevector_x%myLatEnd
          do lonIndex = statevector_x%myLonBeg, statevector_x%myLonEnd
            field4d_x_r8(lonIndex,latIndex,kIndex,stepIndex) = rng_gaussian()
          end do
        end do
      end do
    end do

    ! y
    call gsv_getField(statevector_y,  field4d_y_r8 )
    do kIndex = statevector_y%mykBeg, statevector_y%mykEnd
      do stepIndex = 1, statevector_y%numStep
        do latIndex = statevector_y%myLatBeg, statevector_y%myLatEnd
          do lonIndex = statevector_y%myLonBeg, statevector_y%myLonEnd
            field4d_y_r8(lonIndex,latIndex,kIndex,stepIndex) = rng_gaussian()
          end do
        end do
      end do
    end do

    ! Ly
    call gsv_copy(statevector_y, & ! IN
                  statevector_Ly)  ! OUT
    call adv_statevector_tl( statevector_Ly, & ! INOUT
                             adv_analInc )     ! IN

    ! <x ,L(y)>
    innerProduct1_local = 0.d0
    call gsv_getField(statevector_x,  field4d_x_r8 )
    call gsv_getField(statevector_Ly, field4d_Ly_r8 )
    call euclid(innerProduct1_local, &
         field4d_x_r8(statevector_Ly%myLonBeg:statevector_Ly%myLonEnd,statevector_Ly%myLatBeg:statevector_Ly%myLatEnd,:,:), &
         field4d_Ly_r8(statevector_Ly%myLonBeg:statevector_Ly%myLonEnd,statevector_Ly%myLatBeg:statevector_Ly%myLatEnd,:,:), &
         statevector_Ly%myLonBeg, statevector_Ly%myLonEnd, statevector_Ly%myLatBeg, statevector_Ly%myLatEnd, statevector_Ly%nk, statevector_Ly%numStep)
    write(*,*) "<x     ,L(y)> local = ",innerProduct1_local
    call rpn_comm_allreduce(innerProduct1_local,innerProduct1_global,1,"mpi_double_precision","mpi_sum","GRID",ierr)
    write(*,*) "<x     ,L(y)> global= ",innerProduct1_global
    
    ! L_T(x)
    call gsv_copy(statevector_x, & ! IN
                  statevector_LTx)  ! OUT
    call adv_statevector_ad( statevector_LTx, & ! INOUT
                             adv_analInc )      ! IN
    
    ! <L_T(x),y>
    innerProduct2_local = 0.d0
    call gsv_getField(statevector_LTx, field4d_LTx_r8 )
    call gsv_getField(statevector_y,   field4d_y_r8 )
    call euclid(innerProduct2_local, &
         field4d_LTx_r8(statevector_y%myLonBeg:statevector_y%myLonEnd,statevector_y%myLatBeg:statevector_y%myLatEnd,:,:), &
         field4d_y_r8(statevector_y%myLonBeg:statevector_y%myLonEnd,statevector_y%myLatBeg:statevector_y%myLatEnd,:,:), &
         statevector_y%myLonBeg, statevector_y%myLonEnd, statevector_y%myLatBeg, statevector_y%myLatEnd, statevector_y%nk, statevector_y%numStep)
!    call euclid(innerProduct2_local, controlVector2, controlVector1, 1, cvDim, 1, 1, 1, 1)
    print*,"<Lt(x) ,y   > local = ",innerProduct2_local
    call rpn_comm_allreduce(innerProduct2_local,innerProduct2_global,1,"mpi_double_precision","mpi_sum","GRID",ierr)
    write(*,*) "<Lt(x) ,y   > global= ",innerProduct2_global
    
    ! Results
    call checkAndOutputInnerProd

    call gsv_deallocate(statevector_x)
    call gsv_deallocate(statevector_Ly)
    call gsv_deallocate(statevector_LTx)
    call gsv_deallocate(statevector_y)

  end subroutine check_advectionGSV
  
  !--------------------------------------------------------------------------
  !- Inner product computation
  !--------------------------------------------------------------------------
  subroutine euclid (pvalue,px,py,nis,nie,njs,nje,nk,nStep)
    implicit none
    integer ::  nis,nie,njs,nje,nk,nStep
    real(8) ::  px(nis:nie,njs:nje,nk,nStep), py(nis:nie,njs:nje,nk,nStep)
    real(8) ::  pvalue
    
    integer i,j,k,s

    do s = 1, nStep
      do k = 1,nk
        do j = njs, nje
          do i = nis, nie
            pvalue = pvalue + px(i,j,k,s) * py(i,j,k,s)
          end do
        end do
      end do
    end do

  end subroutine euclid

  !--------------------------------------------------------------------------
  !- Inner product comparison
  !--------------------------------------------------------------------------
  subroutine checkAndOutputInnerProd()
    implicit none

    integer :: fun

    if ( mmpi_myid == 0 ) then
      open(newunit=fun, file="innerProd.txt", status="new", action="write")
      write(fun,'(A20)') test
      write(fun, '(G23.16)') innerProduct1_global
      write(fun, '(G23.16)') innerProduct2_global
      write(fun, '(G23.16)') innerProduct2_global-innerProduct1_global
      write(*,*)
      if ( innerProduct2_global + innerProduct1_global /= 0.d0 ) then
        write(fun, '(G23.16,A1)') 100.d0 * (abs(innerProduct2_global-innerProduct1_global) &
                                          /(0.5d0*(innerProduct2_global+innerProduct1_global)) ), &
                                  '%'
        write(*,'(A20,1X,A,1X,G23.16)') test, ": InnerProd Difference(%) = ", &
             100.d0 * (abs(innerProduct2_global-innerProduct1_global) / &
             (0.5d0*(innerProduct2_global+innerProduct1_global)) )
      else
        write(fun, '(A23)') 'InnerProduct = 0 ERROR'
        write(*,*) 'InnerProduct = 0 !!! Obviously, something went wrong...'
      end if
      close(fun)
    end if

  end subroutine checkAndOutputInnerProd

end program midas_adjointTest
