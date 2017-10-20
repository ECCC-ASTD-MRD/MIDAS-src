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
!! MODULE variableTransforms (prefix="vtr")
!!
!! *Purpose*: To store various functions for variable transforms using inputs from
!!            gridStateVector(s). Outputs are also placed in a GridStateVector.
!!
!--------------------------------------------------------------------------
MODULE variableTransforms_mod
  use mpivar_mod
  use mathPhysConstants_mod
  use timeCoord_mod
  use gridStateVector_mod
  use lamSpectralTransform_mod
  use analysisGrid_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use utilities_mod
  implicit none
  save
  private

  ! public procedures
  public :: vtr_Setup, vtr_Transform

  logical             :: initialized = .false.

  type(struct_gsv)    :: statevector_trial
  integer             :: myLonBeg, myLonEnd, myLatBeg, myLatEnd

CONTAINS

  !--------------------------------------------------------------------------
  ! vtr_Setup
  !--------------------------------------------------------------------------
  subroutine vtr_setup(hco_anl,vco_anl)
    implicit none

    type(struct_hco), pointer :: hco_anl
    type(struct_vco), pointer :: vco_anl
    
    if (initialized) return

    write(*,*) 'vtr_setup: starting'

    initialized  = .true.

    call vtr_setupTrials(hco_anl,vco_anl) ! in

    write(*,*) 'vtr_setup: done'

  end subroutine vtr_setup

  !--------------------------------------------------------------------------
  ! vtr_setupTrials
  !--------------------------------------------------------------------------
  subroutine vtr_setupTrials(hco_anl,vco_anl)
    implicit none

    type(struct_hco), pointer :: hco_anl
    type(struct_vco), pointer :: vco_anl    
 
    call tmg_start(92,'VTR_READTRIALS')
    call gsv_readTrials(hco_anl,vco_anl,statevector_trial)
    call tmg_stop(92)
    
  end subroutine vtr_setupTrials

  !--------------------------------------------------------------------------
  ! vtr_transform
  !--------------------------------------------------------------------------
  subroutine vtr_transform(statevector,transform)
    implicit none
   
    type(struct_gsv) :: statevector
 
    character(len=*), intent(in) :: transform

    integer :: indexStep

    if ( .not. initialized ) then
      write(*,*)
      write(*,*) 'The variable transforms module was NOT initialized'
      call utl_abort('vtr_transform')
    endif

    myLonBeg = statevector%myLonBeg
    myLonEnd = statevector%myLonEnd
    myLatBeg = statevector%myLatBeg
    myLatEnd = statevector%myLatEnd

    select case(trim(transform))

    case ('UVtoVortDiv')
       call UVtoVortDiv(statevector)
    case ('VortDivToPsiChi')
       if ( .not. gsv_varExist(statevector,'QR') .or. .not. gsv_varExist(statevector,'DD') ) then
         write(*,*)
         write(*,*) 'for VortDivToPsiChi, variables QR and DD must be allocated in gridstatevector'
         call utl_abort('vtr_transform')
       end if
       call VortDivToPsiChi(statevector)
    case ('LQtoHU')
       call LQtoHU(statevector)
    case ('HUtoLQ')
       call HUtoLQ(statevector)
    case ('LQtoHU_tlm')
       call LQtoHU_tlm(statevector)
    case ('HUtoLQ_tlm')
       call HUtoLQ_tlm(statevector)
    case default
       write(*,*)
       write(*,*) 'Unsupported function : ', trim(transform)
       call utl_abort('vtr_transform')
    end select

  end subroutine vtr_transform

!--------------------------------------------------------------------------
! LQtoHU
!--------------------------------------------------------------------------
  subroutine LQtoHU(statevector)
    implicit none

    type(struct_gsv) :: statevector
    integer :: i,j,k,indexStep

    real(8), pointer :: hu_ptr(:,:,:,:), lq_ptr(:,:,:,:)

    hu_ptr => gsv_getField_r8(statevector,'HU')
    lq_ptr => gsv_getField_r8(statevector,'HU')

    do indexStep = 1, statevector%numStep
!$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname('HU'))
          do j = statevector%myLatBeg, statevector%myLatEnd
             do i = statevector%myLonBeg, statevector%myLonEnd
                hu_ptr(i,j,k,indexStep) = exp(lq_ptr(i,j,k,indexStep))
             end do
          end do
       end do
!$OMP END PARALLEL DO
   end do

  end subroutine LQtoHU

!--------------------------------------------------------------------------
! HUtoLQ
!--------------------------------------------------------------------------
  subroutine HUtoLQ(statevector)
    implicit none

    type(struct_gsv)    :: statevector
    integer :: i,j,k,indexStep

    real(8), pointer :: hu_ptr(:,:,:,:), lq_ptr(:,:,:,:)

    hu_ptr   => gsv_getField_r8(statevector,'HU')
    lq_ptr   => gsv_getField_r8(statevector,'HU')

    do indexStep = 1, statevector%numStep
!$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname('HU'))
          do j = statevector%myLatBeg, statevector%myLatEnd
             do i = statevector%myLonBeg, statevector%myLonEnd
                lq_ptr(i,j,k,indexStep) = log(max(hu_ptr(i,j,k,indexStep),MPC_MINIMUM_HU_R8))
             end do
          end do
       end do
!$OMP END PARALLEL DO
    end do

  end subroutine HUtoLQ

!--------------------------------------------------------------------------
! LQtoHU_tlm
!--------------------------------------------------------------------------
  subroutine LQtoHU_tlm(statevector)
    implicit none

    type(struct_gsv)    :: statevector
    integer :: i,j,k,indexStep

    real(8), pointer :: hu_ptr(:,:,:,:), lq_ptr(:,:,:,:), lq_trial(:,:,:,:)

    lq_trial => gsv_getField_r8(statevector_trial,'HU')
    hu_ptr   => gsv_getField_r8(statevector      ,'HU')
    lq_ptr   => gsv_getField_r8(statevector      ,'HU')

    do indexStep = 1, statevector%numStep
!$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname('HU'))
          do j = statevector%myLatBeg, statevector%myLatEnd
             do i = statevector%myLonBeg, statevector%myLonEnd       
                hu_ptr(i,j,k,indexStep) = lq_ptr(i,j,k,indexStep)*exp(lq_trial(i,j,k,indexStep))
             end do
          end do
       end do
!$OMP END PARALLEL DO
    end do

  end subroutine LQtoHU_tlm

!--------------------------------------------------------------------------
! HUtoLQ_tlm
!--------------------------------------------------------------------------
  subroutine HUtoLQ_tlm(statevector)
    implicit none

    type(struct_gsv)    :: statevector
    integer :: i,j,k,indexStep

    real(8), pointer :: hu_ptr(:,:,:,:), lq_ptr(:,:,:,:), lq_trial(:,:,:,:)

    lq_trial => gsv_getField_r8(statevector_trial,'HU')
    hu_ptr   => gsv_getField_r8(statevector      ,'HU')
    lq_ptr   => gsv_getField_r8(statevector      ,'HU')

    do indexStep = 1, statevector%numStep
!$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname('HU'))
          do j = statevector%myLatBeg, statevector%myLatEnd
             do i = statevector%myLonBeg, statevector%myLonEnd
                lq_ptr(i,j,k,indexStep) = hu_ptr(i,j,k,indexStep)/exp(lq_trial(i,j,k,indexStep))
             end do
          end do
       end do
!$OMP END PARALLEL DO
    end do

  end subroutine HUtoLQ_tlm

!--------------------------------------------------------------------------
! UVtoVortDiv
!--------------------------------------------------------------------------
  subroutine UVtoVortDiv(statevector)
    implicit none
   
    type(struct_gsv) :: statevector
    integer :: indexStep
 
    real(8), pointer :: uu_ptr(:,:,:,:), vv_ptr(:,:,:,:)
    real(8), pointer :: qr_ptr(:,:,:,:), dd_ptr(:,:,:,:)

    integer :: nlev_M

    uu_ptr => gsv_getField_r8(statevector,'UU')
    vv_ptr => gsv_getField_r8(statevector,'VV')
    qr_ptr => gsv_getField_r8(statevector,'QR')
    dd_ptr => gsv_getField_r8(statevector,'DD')
    nlev_M =  gsv_getNumLev  (statevector,'MM')
    
    if (  statevector%hco%global ) then
       call utl_abort('UVtoVortDiv:global mode not available')
    else

       do indexStep = 1, statevector%numStep
          
          call agd_UVToVortDiv(qr_ptr(:,:,:,indexStep), dd_ptr(:,:,:,indexStep), & ! OUT
                               uu_ptr(:,:,:,indexStep), vv_ptr(:,:,:,indexStep), & ! IN
                               nlev_M )      ! IN
          
       end do
          
    end if

  end subroutine UVtoVortDiv

!--------------------------------------------------------------------------
! VortDivToPsiChi
!--------------------------------------------------------------------------
  subroutine VortDivToPsiChi(statevector)
    implicit none
   
    type(struct_gsv) :: statevector
    integer :: indexStep

    real(8), pointer :: qr_ptr(:,:,:,:), dd_ptr(:,:,:,:)
    real(8), pointer :: pp_ptr(:,:,:,:), cc_ptr(:,:,:,:)

    type(struct_lst)     :: lst_lapl   ! Spectral transform Parameters for Vort/Div -> Psi/Chi
    integer :: nlev_M
    integer :: trunc=125 !!!!!********** UGLY *********!!!!!!!

    qr_ptr => gsv_getField_r8(statevector,'QR')
    dd_ptr => gsv_getField_r8(statevector,'DD')
    pp_ptr => gsv_getField_r8(statevector,'PP')
    cc_ptr => gsv_getField_r8(statevector,'CC')
    nlev_M =  gsv_getNumLev  (statevector,'MM')
    
    if (  statevector%hco%global ) then
       call utl_abort('VortDivToPsiChi:global mode not available')
    else
       call lst_Setup( lst_lapl,                                & ! OUT
                       statevector%hco%ni, statevector%hco%nj,  & ! IN
                       statevector%hco%dlon, trunc,             & ! IN
                      'LatLonMN', nlev_M )                        ! IN

       do indexStep = 1, statevector%numStep
          
          pp_ptr(:,:,:,indexStep) = qr_ptr(:,:,:,indexStep)
          cc_ptr(:,:,:,indexStep) = dd_ptr(:,:,:,indexStep)
         
          ! Vort -> Psi
          call lst_Laplacian( lst_lapl%id,             & ! IN
                              pp_ptr(:,:,:,indexStep), & ! INOUT
                              'Inverse', nlev_M )        ! IN

          ! Div -> Chi
          call lst_Laplacian( lst_lapl%id,             & ! IN
                              cc_ptr(:,:,:,indexStep), & ! INOUT
                              'Inverse', nlev_M )        ! IN

       end do
       
    end if

  end subroutine VortDivToPsiChi

END MODULE variableTransforms_mod