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
  use EarthConstants_mod
  use timeCoord_mod
  use gridStateVector_mod
  use lamSpectralTransform_mod
  use globalSpectralTransform_mod
  use analysisGrid_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use utilities_mod
  implicit none
  save
  private

  ! public procedures
  public :: vtr_Setup, vtr_Transform

  logical             :: trialsInitialized = .false.
  type(struct_hco), pointer :: hco_anl => null()
  type(struct_vco), pointer :: vco_anl => null()

  type(struct_gsv)    :: statevector_trial

CONTAINS

  !--------------------------------------------------------------------------
  ! vtr_Setup
  !--------------------------------------------------------------------------
  subroutine vtr_setup(hco_in,vco_in)
    implicit none

    type(struct_hco), pointer :: hco_in
    type(struct_vco), pointer :: vco_in
    
    if (trialsInitialized) return

    write(*,*) 'vtr_setup: starting'

    hco_anl => hco_in
    vco_anl => vco_in

    write(*,*) 'vtr_setup: done'

  end subroutine vtr_setup

  !--------------------------------------------------------------------------
  ! vtr_setupTrials
  !--------------------------------------------------------------------------
  subroutine vtr_setupTrials()
    implicit none

    call tmg_start(92,'VTR_READTRIALS')
    call gsv_readTrials( hco_anl, vco_anl,        & ! IN
                         statevector_trial,       & ! OUT
                         HUcontainsLQ_opt=.false. )  !IN
    call tmg_stop(92)

    trialsInitialized = .true.
    
  end subroutine vtr_setupTrials

  !--------------------------------------------------------------------------
  ! vtr_transform
  !--------------------------------------------------------------------------
  subroutine vtr_transform(statevector,transform)
    implicit none
   
    type(struct_gsv) :: statevector
 
    character(len=*), intent(in) :: transform

    integer :: stepIndex

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
    case ('UVtoPsiChi')
       if ( .not. gsv_varExist(statevector,'PP') .or. .not. gsv_varExist(statevector,'CC') ) then
         write(*,*)
         write(*,*) 'for UVToPsiChi, variables PP and CC must be allocated in gridstatevector'
         call utl_abort('vtr_transform')
       end if
       call UVtoPsiChi(statevector)
    case ('LQtoHU')
      if ( .not. gsv_varExist(statevector,'LQ') .and. gsv_varExist(statevector,'HU')  ) then
        call LQtoHU(statevector)
      else
        write(*,*) 'vtr_transform: skipping LQtoHU'
      end if
    case ('HUtoLQ')
       call HUtoLQ(statevector)
    case ('LQtoHU_tlm')
      if ( .not. gsv_varExist(statevector,'LQ') .and. gsv_varExist(statevector,'HU')  ) then
        call LQtoHU_tlm(statevector)
      else
        write(*,*) 'vtr_transform: skipping LQtoHU_tlm'
      end if
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
    integer :: i,j,k,stepIndex

    real(8), pointer :: hu_ptr(:,:,:,:), lq_ptr(:,:,:,:)

    hu_ptr => gsv_getField_r8(statevector,'HU')
    lq_ptr => gsv_getField_r8(statevector,'HU')

    do stepIndex = 1, statevector%numStep
!$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname('HU'))
          do j = statevector%myLatBeg, statevector%myLatEnd
             do i = statevector%myLonBeg, statevector%myLonEnd
                hu_ptr(i,j,k,stepIndex) = exp(lq_ptr(i,j,k,stepIndex))
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
    integer :: i,j,k,stepIndex

    real(8), pointer :: hu_ptr(:,:,:,:), lq_ptr(:,:,:,:)

    hu_ptr   => gsv_getField_r8(statevector,'HU')
    lq_ptr   => gsv_getField_r8(statevector,'HU')

    do stepIndex = 1, statevector%numStep
!$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname('HU'))
          do j = statevector%myLatBeg, statevector%myLatEnd
             do i = statevector%myLonBeg, statevector%myLonEnd
                lq_ptr(i,j,k,stepIndex) = log(max(hu_ptr(i,j,k,stepIndex),MPC_MINIMUM_HU_R8))
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
    integer :: i,j,k,stepIndex

    real(8), pointer :: hu_ptr(:,:,:,:), lq_ptr(:,:,:,:), hu_trial(:,:,:,:)

    if ( .not. trialsInitialized ) call vtr_setupTrials()

    hu_trial => gsv_getField_r8(statevector_trial,'HU')
    hu_ptr   => gsv_getField_r8(statevector      ,'HU')
    lq_ptr   => gsv_getField_r8(statevector      ,'HU')

    do stepIndex = 1, statevector%numStep
!$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname('HU'))
          do j = statevector%myLatBeg, statevector%myLatEnd
             do i = statevector%myLonBeg, statevector%myLonEnd       
               hu_ptr(i,j,k,stepIndex) =  lq_ptr(i,j,k,stepIndex)*hu_trial(i,j,k,stepIndex)
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
    integer :: i,j,k,stepIndex

    real(8), pointer :: hu_ptr(:,:,:,:), lq_ptr(:,:,:,:), hu_trial(:,:,:,:)

    if ( .not. trialsInitialized ) call vtr_setupTrials()

    hu_trial => gsv_getField_r8(statevector_trial,'HU')
    hu_ptr   => gsv_getField_r8(statevector      ,'HU')
    lq_ptr   => gsv_getField_r8(statevector      ,'HU')

    do stepIndex = 1, statevector%numStep
!$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = 1, gsv_getNumLev(statevector,vnl_varLevelFromVarname('HU'))
          do j = statevector%myLatBeg, statevector%myLatEnd
             do i = statevector%myLonBeg, statevector%myLonEnd
               lq_ptr(i,j,k,stepIndex) = hu_ptr(i,j,k,stepIndex)/hu_trial(i,j,k,stepIndex)
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
    integer :: stepIndex
 
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

       do stepIndex = 1, statevector%numStep
          
          call agd_UVToVortDiv(qr_ptr(:,:,:,stepIndex), dd_ptr(:,:,:,stepIndex), & ! OUT
                               uu_ptr(:,:,:,stepIndex), vv_ptr(:,:,:,stepIndex), & ! IN
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
    integer :: stepIndex

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

       do stepIndex = 1, statevector%numStep
          
          pp_ptr(:,:,:,stepIndex) = qr_ptr(:,:,:,stepIndex)
          cc_ptr(:,:,:,stepIndex) = dd_ptr(:,:,:,stepIndex)
         
          ! Vort -> Psi
          call lst_Laplacian( lst_lapl%id,             & ! IN
                              pp_ptr(:,:,:,stepIndex), & ! INOUT
                              'Inverse', nlev_M )        ! IN

          ! Div -> Chi
          call lst_Laplacian( lst_lapl%id,             & ! IN
                              cc_ptr(:,:,:,stepIndex), & ! INOUT
                              'Inverse', nlev_M )        ! IN

       end do
       
    end if

  end subroutine VortDivToPsiChi

!--------------------------------------------------------------------------
! UVtoPsiChi
!--------------------------------------------------------------------------
  subroutine UVtoPsiChi(statevector)
    implicit none
   
    type(struct_gsv) :: statevector

    integer :: stepIndex 
    real(8), pointer :: uu_ptr(:,:,:,:), vv_ptr(:,:,:,:)
    real(8), pointer :: psi_ptr(:,:,:,:), chi_ptr(:,:,:,:)
    real(8), allocatable :: gridState(:,:,:), spectralState(:,:,:)
    real(8) :: dla2
    integer :: nlev_M, levIndex
    integer :: ila_mpiglobal, ila_mpilocal

    ! spectral transform configuration (saved)
    integer, save :: gstID = -1
    integer, save :: nla_mpilocal, maxMyNla, ntrunc
    integer, save :: mymBeg,mymEnd,mymSkip,mymCount
    integer, save :: mynBeg,mynEnd,mynSkip,mynCount
    integer, save, pointer :: ilaList_mpiglobal(:), ilaList_mpilocal(:)

    write(*,*) 'UVtoPsiChi: starting'
    call flush(6)

    uu_ptr  => gsv_getField_r8(statevector,'UU')
    vv_ptr  => gsv_getField_r8(statevector,'VV')
    psi_ptr => gsv_getField_r8(statevector,'PP')
    chi_ptr => gsv_getField_r8(statevector,'CC')
    nlev_M = gsv_getNumLev(statevector,'MM')
    
    if ( .not. statevector%hco%global ) then
      call utl_abort('UVtoPsiChi:only global is available')
    else

      if ( gstID < 0 ) then
        !ntrunc = statevector%nj
        ntrunc = 180
        gstID = gst_setup(statevector%ni,statevector%nj,ntrunc,2*nlev_M)
        call mpivar_setup_m(ntrunc,mymBeg,mymEnd,mymSkip,mymCount)
        call mpivar_setup_n(ntrunc,mynBeg,mynEnd,mynSkip,mynCount)
        call gst_ilaList_mpiglobal(ilaList_mpiglobal,nla_mpilocal,maxMyNla,gstID,mymBeg,mymEnd,mymSkip,mynBeg,mynEnd,mynSkip)
        call gst_ilaList_mpilocal(ilaList_mpilocal,gstID,mymBeg,mymEnd,mymSkip,mynBeg,mynEnd,mynSkip)
      end if

      dla2   = dble(ra)*dble(ra)
      allocate(gridState(statevector%lonPerPE,statevector%latPerPE,2*nlev_M))
      allocate(spectralState(nla_mpilocal,2,2*nlev_M))

      do stepIndex = 1, statevector%numStep

        write(*,*) 'copying into gridstate'
        call flush(6)

        gridState(:,:,1:nlev_M)            = uu_ptr(:,:,:,stepIndex)
        gridState(:,:,(nlev_M+1):2*nlev_M) = vv_ptr(:,:,:,stepIndex)

        write(*,*) 'first tranform'; call flush(6)

        call gst_setID(gstID)
        call gst_gdsp(spectralState,gridState,nlev_M)

        write(*,*) 'scale spectral state'; call flush(6)
        do levIndex = 1, 2*nlev_M
          do ila_mpilocal = 1, nla_mpilocal
            ila_mpiglobal = ilaList_mpiglobal(ila_mpilocal)
            spectralState(ila_mpilocal,:,levIndex) =   &
                             spectralState(ila_mpilocal,:,levIndex) *   &
                             dla2*gst_getR1snp1(ila_mpiglobal)
          enddo
        enddo

        write(*,*) 'second transform'; call flush(6)
        call gst_speree(spectralState,gridState)

        write(*,*) 'copying back to psi/chi'
        call flush(6)

        psi_ptr(:,:,:,stepIndex) = gridState(:,:,1:nlev_M)
        chi_ptr(:,:,:,stepIndex) = gridState(:,:,(nlev_M+1):2*nlev_M)

      end do

      write(*,*) 'deallocate'; call flush(6)
      deallocate(gridState)
      deallocate(spectralState)

    end if

    write(*,*) 'UVtoPsiChi: finished'
    call flush(6)

  end subroutine UVtoPsiChi

END MODULE variableTransforms_mod
