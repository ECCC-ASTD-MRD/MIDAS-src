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

!--------------------------------------------------------------------------
!!
!! *Purpose*: Main program for adjoint test applications: 
!!            <x,L(y)> = <L^T(x),y>
!!
!--------------------------------------------------------------------------
program midas_adjointTest
  use ramDisk_mod
  use utilities_mod
  use mpiVar_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use gridStateVector_mod
  use analysisGrid_mod
  use variabletransforms_mod
  use bmatrixhi_mod
  use bmatrixensemble_mod
  implicit none

  type(struct_vco),       pointer :: vco_anl  => null()
  type(struct_vco),       pointer :: vco_trl  => null()
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

  write(*,*) " --------------------------------------------"
  write(*,*) " ---  START OF MAIN PROGRAM ajointTest    ---"
  write(*,*) " ---  Adjoint Test                        ---"
  write(*,*) " --------------------------------------------"

  !
  !- 1.  Settings and module initializations
  !
  write(*,*)
  write(*,*) '> midas-adjointTest: setup - START'

  !- 1.1 mpi
  call mpi_initialize  

  !- 1.2 timings
  call tmg_init(mpi_myid, 'TMG_ADJOINTTEST' )
  call tmg_start(1,'MAIN')

  !- 1.3 RAM disk usage
  call ram_setup

  !- 1.4 Temporal grid
  call tim_setup(fileNameForDate_opt='./trlm_01')

  !- 1.6 Constants
  if (mpi_myid == 0) call mpc_printConstants(6)

  !- 1.8 Variables of the model states
  call gsv_setup

  !- 1.9 Set the horizontal domain
  call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS') ! IN
  if ( hco_anl % global ) then
    call agd_SetupFromHCO( hco_anl ) ! IN
  else
    !- Iniatilized the core (Non-Exteded) analysis grid
    call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID') ! IN
    !- Setup the LAM analysis grid metrics
    call agd_SetupFromHCO( hco_anl, hco_core ) ! IN
  end if

  !- 1.10 Initialize the vertical coordinate
  call vco_SetupFromFile( vco_anl,        & ! OUT
                          './analysisgrid') ! IN

  write(*,*)
  write(*,*) '> midas-adjointTest: setup - END'
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !- 1.11 Variable transforms
  call vtr_Setup(hco_anl, vco_anl)

  !
  !- 2.  The tests
  !
 
  !- 2.1 LQ to HU conversion
  write(*,*)
  write(*,*) '> midas-adjointTest: LQ to HU'
  call check_LQtoHU

  !- 2.2 Bhi
  write(*,*)
  write(*,*) '> midas-adjointTest: Bhi'
  if (hco_anl%global) then
    call check_bhi_global
  !!else
  !  ! call check_bhi_lam
  end if

  !- 2.3 Bens
  write(*,*)
  write(*,*) '> midas-adjointTest: Bens'
  call check_bens

  !
  !- 3.  Ending
  !
  write(*,*)
  write(*,*) '> midas-adjointTest: Ending'
  call tmg_stop(1)
  call tmg_terminate(mpi_myid, 'TMG_ADJOINTTEST' )

  call rpn_comm_finalize(ierr) 

contains
  !
  !- check LQtoHU
  !
  subroutine check_LQtoHU ()
    implicit none

    ! x = HU
    ! y = LQ
    ! L = LQtoHU_tlm
    call gsv_allocate(statevector_x  , tim_nstepobsinc, hco_anl, vco_anl, &
                      varNames_opt=(/'HU'/), mpi_local_opt=.true.)
    call gsv_allocate(statevector_y  , tim_nstepobsinc, hco_anl, vco_anl, &
                      varNames_opt=(/'HU'/), mpi_local_opt=.true.)
    call gsv_allocate(statevector_Ly , tim_nstepobsinc, hco_anl, vco_anl, &
                      varNames_opt=(/'HU'/), mpi_local_opt=.true.)
    call gsv_allocate(statevector_LTx, tim_nstepobsinc, hco_anl, vco_anl, &
                      varNames_opt=(/'HU'/), mpi_local_opt=.true.)
    
    statevector_x%gd3d_r8(:,:,:) = 0.0133d0
    statevector_y%gd3d_r8(:,:,:) = log(0.0133d0)
    
    ! <x,L(y)>
    call gsv_copy(statevector_y, & ! IN
                  statevector_Ly)  ! OUT
    call vtr_transform( statevector_Ly, & ! INOUT
         'LQtoHU_tlm' ) ! IN
    
    innerProduct1_local = 0.d0
    call euclid(innerProduct1_local, &
       statevector_x %gd3d_r8(statevector_Ly%myLonBeg:statevector_Ly%myLonEnd,statevector_Ly%myLatBeg:statevector_Ly%myLatEnd,:), &
       statevector_Ly%gd3d_r8(statevector_Ly%myLonBeg:statevector_Ly%myLonEnd,statevector_Ly%myLatBeg:statevector_Ly%myLatEnd,:), &
       statevector_Ly%myLonBeg, statevector_Ly%myLonEnd, statevector_Ly%myLatBeg, statevector_Ly%myLatEnd, statevector_Ly%nk)
    write(*,*) "<x     ,L(y)> local = ",innerProduct1_local
    call rpn_comm_allreduce(innerProduct1_local,innerProduct1_global,1,"mpi_double_precision","mpi_sum","GRID",ierr)
    write(*,*) "<x     ,L(y)> global= ",innerProduct1_global
    
    ! <LT(x),y>
    call gsv_copy(statevector_x, & ! IN
         statevector_LTx)  ! OUT
    call vtr_transform( statevector_LTx, & ! INOUT
         'LQtoHU_tlm' ) ! IN
    
    innerProduct2_local = 0.d0
    call euclid(innerProduct2_local, &
         statevector_LTx%gd3d_r8(statevector_Ly%myLonBeg:statevector_Ly%myLonEnd,statevector_Ly%myLatBeg:statevector_Ly%myLatEnd,:), &
         statevector_y  %gd3d_r8(statevector_Ly%myLonBeg:statevector_Ly%myLonEnd,statevector_Ly%myLatBeg:statevector_Ly%myLatEnd,:), &
         statevector_Ly%myLonBeg, statevector_Ly%myLonEnd, statevector_Ly%myLatBeg, statevector_Ly%myLatEnd, statevector_Ly%nk)
    print*,"<Lt(x) ,y   > local = ",innerProduct2_local
    call rpn_comm_allreduce(innerProduct2_local,innerProduct2_global,1,"mpi_double_precision","mpi_sum","GRID",ierr)
    write(*,*) "<Lt(x) ,y   > global= ",innerProduct2_global
    
    ! Results
    call checkInnerProd (innerProduct1_global, innerProduct2_global, 'LQtoHU')
  
    call gsv_deallocate(statevector_x  )
    call gsv_deallocate(statevector_y  )
    call gsv_deallocate(statevector_Ly )
    call gsv_deallocate(statevector_LTx)

  end subroutine check_LQtoHU

  !
  !- check Bhi global
  !
  subroutine check_bhi_global ()
    implicit none

    call bhi_Setup( hco_anl, vco_anl, & ! IN
                  cvdim )             ! OUT

    call gsv_allocate(statevector_x  , tim_nstepobsinc, hco_anl, vco_anl, &
                      mpi_local_opt=.true.)
    call gsv_allocate(statevector_Ly , tim_nstepobsinc, hco_anl, vco_anl, &
                      mpi_local_opt=.true.)

    allocate ( controlVector1(cvDim) )

    ! y
    controlVector1(:) = 2.4d0
    
    ! x
    statevector_x%gd3d_r8(:,:,:) = 13.3d0
    
    ! Ly
    call bhi_bSqrt( controlVector1, & ! IN
         statevector_Ly )  ! OUT
    
    ! <x ,L(y)>
    innerProduct1_local = 0.d0
    call euclid(innerProduct1_local, &
         statevector_x %gd3d_r8(statevector_Ly%myLonBeg:statevector_Ly%myLonEnd,statevector_Ly%myLatBeg:statevector_Ly%myLatEnd,:), &
         statevector_Ly%gd3d_r8(statevector_Ly%myLonBeg:statevector_Ly%myLonEnd,statevector_Ly%myLatBeg:statevector_Ly%myLatEnd,:), &
         statevector_Ly%myLonBeg, statevector_Ly%myLonEnd, statevector_Ly%myLatBeg, statevector_Ly%myLatEnd, statevector_Ly%nk)
    write(*,*) "<x     ,L(y)> local = ",innerProduct1_local
    call rpn_comm_allreduce(innerProduct1_local,innerProduct1_global,1,"mpi_double_precision","mpi_sum","GRID",ierr)
    write(*,*) "<x     ,L(y)> global= ",innerProduct1_global
    
    ! L_T(x)
    allocate ( controlVector2(cvDim) )
    controlVector2 = 0.d0
    call bhi_bSqrtAd( statevector_x, &  ! IN
         controlVector2 ) ! OUT
    
    ! <L_T(x),y>
    innerProduct2_local = 0.d0
    call euclid(innerProduct2_local, controlVector2, controlVector1, 1, cvDim, 1, 1, 1)
    print*,"<Lt(x) ,y   > local = ",innerProduct2_local
    call rpn_comm_allreduce(innerProduct2_local,innerProduct2_global,1,"mpi_double_precision","mpi_sum","GRID",ierr)
    write(*,*) "<Lt(x) ,y   > global= ",innerProduct2_global
    
    ! Results
    call checkInnerProd (innerProduct1_global, innerProduct2_global, 'Bhi-global')
    
    deallocate(controlVector2)
    deallocate(controlVector1)

    call gsv_deallocate(statevector_x)
    call gsv_deallocate(statevector_Ly)

  end subroutine check_bhi_global

  !
  !- check Bens
  !
  subroutine check_bens()
    implicit none

    call ben_Setup( hco_anl, vco_anl, & ! IN
                    cvdim )             ! OUT

    call gsv_allocate(statevector_x  , tim_nstepobsinc, hco_anl, vco_anl, &
                      mpi_local_opt=.true.)
    call gsv_allocate(statevector_Ly , tim_nstepobsinc, hco_anl, vco_anl, &
                      mpi_local_opt=.true.)

    allocate ( controlVector1(cvDim) )

    ! x
    statevector_x%gd3d_r8(:,:,:) = 13.3d0
    
    ! y
    !controlVector1(:) = 2.4d0
    controlVector1 = 0.d0
    call ben_bSqrtAd( statevector_x, &  ! IN
                      controlVector1)  ! OUT

    ! Ly
    call ben_bSqrt( controlVector1, & ! IN
                    statevector_Ly )  ! OUT
    
    ! <x ,L(y)>
    innerProduct1_local = 0.d0
    call euclid(innerProduct1_local, &
         statevector_x %gd3d_r8(statevector_Ly%myLonBeg:statevector_Ly%myLonEnd,statevector_Ly%myLatBeg:statevector_Ly%myLatEnd,:), &
         statevector_Ly%gd3d_r8(statevector_Ly%myLonBeg:statevector_Ly%myLonEnd,statevector_Ly%myLatBeg:statevector_Ly%myLatEnd,:), &
         statevector_Ly%myLonBeg, statevector_Ly%myLonEnd, statevector_Ly%myLatBeg, statevector_Ly%myLatEnd, statevector_Ly%nk)
    write(*,*) "<x     ,L(y)> local = ",innerProduct1_local
    call rpn_comm_allreduce(innerProduct1_local,innerProduct1_global,1,"mpi_double_precision","mpi_sum","GRID",ierr)
    write(*,*) "<x     ,L(y)> global= ",innerProduct1_global
    
    ! L_T(x)
    allocate ( controlVector2(cvDim) )
    controlVector2 = 0.d0
    call ben_bSqrtAd( statevector_x, &  ! IN
                      controlVector2 )  ! OUT
    
    ! <L_T(x),y>
    innerProduct2_local = 0.d0
    call euclid(innerProduct2_local, controlVector2, controlVector1, 1, cvDim, 1, 1, 1)
    print*,"<Lt(x) ,y   > local = ",innerProduct2_local
    call rpn_comm_allreduce(innerProduct2_local,innerProduct2_global,1,"mpi_double_precision","mpi_sum","GRID",ierr)
    write(*,*) "<Lt(x) ,y   > global= ",innerProduct2_global
    
    ! Results
    call checkInnerProd (innerProduct1_global, innerProduct2_global, 'Bens')
    
    deallocate(controlVector2)
    deallocate(controlVector1)

    call gsv_deallocate(statevector_x)
    call gsv_deallocate(statevector_Ly)

  end subroutine check_bens

  !
  !- Inner product computation
  !
  subroutine euclid (pvalue,px,py,nis,nie,njs,nje,nk)
    implicit none
    integer ::  nis,nie,njs,nje,nk
    real(8) ::  px(nis:nie,njs:nje,nk), py(nis:nie,njs:nje,nk)
    real(8) ::  pvalue
    
    integer i,j,k
    
    do k = 1,nk
      do j = njs, nje
        do i = nis, nie
          pvalue = pvalue + px(i,j,k) * py(i,j,k)
        end do
      end do
    end do
    
  end subroutine euclid

  !
  !- Inner product comparison
  !
  subroutine checkInnerProd (innerProduct1, innerProduct2, testName)
    implicit none
    real(8) ::  innerProduct1, innerProduct2
    character(len=*) :: testName

    if ( mpi_myid == 0 ) then
      write(*,*)
      write(*,'(A20,1X,A,1X,G23.16)') testName, ": InnerProd Difference(%) = ", &
           100.d0 * (abs(innerProduct2_global-innerProduct1_global) / &
           (0.5d0*(innerProduct2_global+innerProduct1_global)) )
    end if

  end subroutine checkInnerProd

end program midas_adjointTest
