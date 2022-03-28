!-------------------------------------- LICENCE BEGIN ------------------------------------
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

program midas_adjointTest
  !
  ! :Purpose: Main program for adjoint test applications: 
  !           <x,L(y)> = <L^T(x),y>
  !
  use version_mod
  use codePrecision_mod
  use ramDisk_mod
  use utilities_mod
  use mpi_mod
  use mpiVar_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use gridStateVector_mod
  use gridVariableTransforms_mod
  use bmatrixhi_mod
  use bmatrixensemble_mod
  use randomNumber_mod
  use advection_mod
  use ensembleStateVector_mod
  use localization_mod
  use lamBMatrixHI_mod
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

  integer :: get_max_rss, ierr, cvDim

  call ver_printNameAndVersion('adjointTest','Various tests of adjoint codes')

  !
  !- 1.  Settings and module initializations
  !
  write(*,*)
  write(*,*) '> midas-adjointTest: setup - START'

  !- 1.1 mpi
  call mpi_initialize

  !- 1.2 timings
  call tmg_init(mpi_myid, 'TMG_INFO')
  call utl_tmg_start(0,'MAIN')

  !- 1.3 RAM disk usage
  call ram_setup

  !- 1.4 Temporal grid
  call tim_setup()

  !- 1.6 Constants
  if ( mpi_myid == 0 ) then
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

  !
  !- 2.  The tests
  !
 
  !- 2.1 Bhi
  !write(*,*)
  !write(*,*) '> midas-adjointTest: Bhi'
  !call check_bhi

  !- 2.2 Bens
  write(*,*)
  write(*,*) '> midas-adjointTest: Bens'
  call check_bens

  !- 2.3 AdvectionENS
  !write(*,*)
  !write(*,*) '> midas-adjointTest: AdvectionENS'
  !call check_advectionENS

  !- 2.4 AdvectionGSV
  !write(*,*)
  !write(*,*) '> midas-adjointTest: AdvectionGSV'
  !call check_advectionGSV

  !- 2.5 Localization
  !write(*,*)
  !write(*,*) '> midas-adjointTest: Localization'
  !call check_loc

  !- 2.6 Localization
  !write(*,*)
  !write(*,*) '> midas-adjointTest: Addmem'
  !call check_addmem

  !
  !- 3.  Ending
  !
  write(*,*)
  write(*,*) '> midas-adjointTest: Ending'
  call tmg_stop(0)
  call tmg_terminate(mpi_myid, 'TMG_INFO')

  call rpn_comm_finalize(ierr) 

contains

  !--------------------------------------------------------------------------
  !- check Bhi
  !--------------------------------------------------------------------------
  subroutine check_bhi()
    implicit none

    integer :: seed, kIndex, stepIndex, latIndex, lonIndex, cvIndex
    real(8), pointer :: field4d_r8(:,:,:,:), field3d_Ly_r8(:,:,:), field3d_x_r8(:,:,:)

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
    call rng_setup(abs(seed+mpi_myid))
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
    call gsv_getField3d(statevector_Ly,field3d_Ly_r8)
    call gsv_getField3d(statevector_x, field3d_x_r8 )
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
    call checkInnerProd ('Bhi')
    
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
    call rng_setup(abs(seed+mpi_myid))
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
    call checkInnerProd ('Bens')
    
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
    call rng_setup(abs(seed+mpi_myid))
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
    call checkInnerProd ('Loc')
    
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
!!$    call rng_setup(abs(seed+mpi_myid))
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
!!$    call checkInnerProd ('addMem')
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
    call rng_setup(abs(seed+mpi_myid))
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
    call checkInnerProd ('advection')

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
    call rng_setup(abs(seed+mpi_myid))
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
    call checkInnerProd ('advection')

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
  subroutine checkInnerProd (testName)
    implicit none
    character(len=*) :: testName

    if ( mpi_myid == 0 ) then
      write(*,*)
      if ( innerProduct2_global + innerProduct1_global /= 0.d0 ) then
        write(*,'(A20,1X,A,1X,G23.16)') testName, ": InnerProd Difference(%) = ", &
             100.d0 * (abs(innerProduct2_global-innerProduct1_global) / &
             (0.5d0*(innerProduct2_global+innerProduct1_global)) )
      else
        write(*,*) 'InnerProduct = 0 !!! Obviously, something went wrong...'
      end if
    end if

  end subroutine checkInnerProd

end program midas_adjointTest
