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

module calcPressure_mod
  !
  ! MODULE calcPressure (prefix='clp' category='2. High-level data objects')
  !
  ! **Purpose**: Calculate pressure on the gridstatevector (nl/tl/ad).
  !
  ! @Author M. Bani Shahabadi, Jan 2019
  !
  use gridStateVector_mod
  use verticalCoord_mod
  implicit none
  save
  private

  ! public structure definition

  ! public subroutines and functions
  public :: clp_calcPressure_nl, clp_calcPressure_tl, clp_calcPressure_ad

  contains

  subroutine clp_calcPressure_nl(statevector, beSilent_opt)

    implicit none
    type(struct_gsv), intent(inout) :: statevector
    logical, optional :: beSilent_opt

    real(kind=8), allocatable   :: Psfc(:,:)
    real(kind=8), pointer       :: Pressure_out(:,:,:) 
    real(kind=8), pointer       :: dP_dPsfc_out(:,:,:)
    real(kind=8), pointer       :: field_Psfc(:,:,:,:)
    integer                     :: jobs, status, stepIndex
    logical                     :: beSilent
    real(8), pointer            :: P_T(:,:,:,:) => null()
    real(8), pointer            :: P_M(:,:,:,:) => null()

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .false.
    end if

    if ( .not. beSilent ) write(*,*) 'clp_calcPressure_nl: computing pressure on staggered or UNstaggered levels'

    if ( .not. gsv_varExist(statevector,'P_T') .or. .not. gsv_varExist(statevector,'P_M') .or. .not. gsv_varExist(statevector,'P0')) then
      call utl_abort('clp_calcPressure_nl: P_T/P_M/P0 do not exist in statevector!')
    end if

    P_T => gsv_getField_r8(statevector,'P_T')
    P_M => gsv_getField_r8(statevector,'P_M')

    allocate(Psfc(statevector%myLonBeg:statevector%myLonEnd, &
                  statevector%myLatBeg:statevector%myLatEnd))

    do stepIndex = 1, statevector%numStep

      field_Psfc => gsv_getField_r8(statevector,'P0')
      Psfc(:,:) = field_Psfc(:,:,1,stepIndex)

      ! P_T
      nullify(Pressure_out)
      status = vgd_levels(statevector%vco%vgrid, &
                        ip1_list=statevector%vco%ip1_M, &
                        levels=Pressure_out, &
                        sfc_field=Psfc, &
                        in_log=.false.)
      if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_levels')
      P_M(:,:,:,stepIndex) = Pressure_out(:,:,:)
      deallocate(Pressure_out)

      ! P_M
      nullify(Pressure_out)
      status = vgd_levels(statevector%vco%vgrid, &
                        ip1_list=statevector%vco%ip1_T, &
                        levels=Pressure_out, &
                        sfc_field=Psfc, &
                        in_log=.false.)
      if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_levels')
      P_T(:,:,:,stepIndex) = Pressure_out(:,:,:)
      deallocate(Pressure_out)

      deallocate(Psfc)

      if ( .not. beSilent .and. stepIndex == 1 ) then
        write(*,*) 'stepIndex=',stepIndex, ',P_M='
        write(*,*) P_M(statevector%myLonBeg,statevector%myLatBeg,:,stepIndex)
        write(*,*) 'stepIndex=',stepIndex, ',P_T='
        write(*,*) P_T(statevector%myLonBeg,statevector%myLatBeg,:,stepIndex)
      end if

    end do

  end subroutine clp_calcPressure_nl


  subroutine clp_calcPressure_tl(statevector_trial, statevector, beSilent_opt)

    implicit none
    type(struct_gsv), intent(inout) :: statevector, statevector_trial
    logical, optional :: beSilent_opt

    real(kind=8), allocatable   :: Psfc(:,:)
    real(kind=8), pointer       :: delPsfc(:,:,:,:) => null()
    real(kind=8), pointer       :: field_Psfc(:,:,:,:) => null()
    real(8), pointer            :: delP_T(:,:,:,:) => null()
    real(8), pointer            :: delP_M(:,:,:,:) => null()
    real(8), pointer            :: dP_dPsfc_T(:,:,:) => null()
    real(8), pointer            :: dP_dPsfc_M(:,:,:) => null()
    integer                     :: jobs, status, stepIndex,lonIndex,latIndex
    integer                     :: lev_M, lev_T, nlev_T, nlev_M, numStep
    logical                     :: beSilent

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .false.
    end if

    if ( .not. beSilent ) write(*,*) 'clp_calcPressure_tl: computing delP_T/delP_M on the gridstatevector'

    if ( .not. gsv_varExist(statevector,'P_T') .or. .not. gsv_varExist(statevector,'P_M') .or. .not. gsv_varExist(statevector,'P0') ) then
      call utl_abort('clp_calcPressure_tl: P_T/P_M/P0 do not exist in statevector!')
    end if

    delP_T => gsv_getField_r8(statevector,'P_T')
    delP_M => gsv_getField_r8(statevector,'P_M')
    delPsfc => gsv_getField_r8(statevector,'P0')
    field_Psfc => gsv_getField_r8(statevector_trial,'P0')

    nlev_T = gsv_getNumLev(statevector_trial,'TH')
    nlev_M = gsv_getNumLev(statevector_trial,'MM')
    numStep = statevector_trial%numstep

    allocate(Psfc(statevector%myLonBeg:statevector%myLonEnd, &
                  statevector%myLatBeg:statevector%myLatEnd))

    do stepIndex = 1, numStep

      Psfc(:,:) = field_Psfc(:,:,1,stepIndex)

      ! dP_dPsfc_M
      nullify(dP_dPsfc_M)
      status = vgd_dpidpis(statevector%vco%vgrid, &
                           statevector%vco%ip1_M, &
                           dP_dPsfc_M, &
                           Psfc)
      if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_dpidpis')
      ! calculate delP_M
      do lev_M = 1, nlev_M
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do lonIndex = statevector%myLonBeg, statevector%myLonEnd
            delP_M(lonIndex,latIndex,lev_M,stepIndex) = dP_dPsfc_M(lonIndex,latIndex,lev_M) * delPsfc(lonIndex,latIndex,1,stepIndex)
          end do
        end do
      end do
      deallocate(dP_dPsfc_M)

      ! dP_dPsfc_T
      nullify(dP_dPsfc_T)
      status = vgd_dpidpis(statevector%vco%vgrid, &
                           statevector%vco%ip1_T, &
                           dP_dPsfc_T, &
                           Psfc)
      if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_dpidpis')
      ! calculate delP_T
      do lev_T = 1, nlev_T
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do lonIndex = statevector%myLonBeg, statevector%myLonEnd
            delP_T(lonIndex,latIndex,lev_T,stepIndex) = dP_dPsfc_T(lonIndex,latIndex,lev_T) * delPsfc(lonIndex,latIndex,1,stepIndex)
          end do
        end do
      end do
      deallocate(dP_dPsfc_T)

      if ( .not. beSilent .and. stepIndex == 1 ) then
        write(*,*) 'stepIndex=',stepIndex, ',delP_M='
        write(*,*) delP_M(statevector%myLonBeg,statevector%myLatBeg,:,stepIndex)
        write(*,*) 'stepIndex=',stepIndex, ',delP_T='
        write(*,*) delP_T(statevector%myLonBeg,statevector%myLatBeg,:,stepIndex)
      end if

    end do

    deallocate(Psfc)

  end subroutine clp_calcPressure_tl


  subroutine clp_calcPressure_ad(statevector_trial, statevector, beSilent_opt)

    implicit none
    type(struct_gsv), intent(inout) :: statevector, statevector_trial
    logical, optional :: beSilent_opt

    real(kind=8), allocatable   :: Psfc(:,:)
    real(kind=8), pointer       :: delPsfc(:,:,:,:) => null()
    real(kind=8), pointer       :: field_Psfc(:,:,:,:) => null()
    real(8), pointer            :: delP_T(:,:,:,:) => null()
    real(8), pointer            :: delP_M(:,:,:,:) => null()
    real(8), pointer            :: dP_dPsfc_T(:,:,:) => null()
    real(8), pointer            :: dP_dPsfc_M(:,:,:) => null()
    integer                     :: jobs, status, stepIndex,lonIndex,latIndex
    integer                     :: lev_M, lev_T, nlev_T, nlev_M, numStep
    logical                     :: beSilent

    ! Add calculation of the GZ/P  

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .false.
    end if

    if ( .not. beSilent ) write(*,*) 'clp_calcPressure_ad: computing delP_T/delP_M on the gridstatevector'

    if ( .not. gsv_varExist(statevector,'P_T') .or. .not. gsv_varExist(statevector,'P_M') .or. .not. gsv_varExist(statevector,'P0') ) then
      call utl_abort('clp_calcPressure_ad: P_T/P_M/P0 do not exist in statevector!')
    end if

    delP_T => gsv_getField_r8(statevector,'P_T')
    delP_M => gsv_getField_r8(statevector,'P_M')
    delPsfc => gsv_getField_r8(statevector,'P0')
    field_Psfc => gsv_getField_r8(statevector_trial,'P0')

    nlev_T = gsv_getNumLev(statevector_trial,'TH')
    nlev_M = gsv_getNumLev(statevector_trial,'MM')
    numStep = statevector_trial%numstep

    allocate(Psfc(statevector%myLonBeg:statevector%myLonEnd, &
                  statevector%myLatBeg:statevector%myLatEnd))

    do stepIndex = 1, numStep

      Psfc(:,:) = field_Psfc(:,:,1,stepIndex)

      ! dP_dPsfc_M
      nullify(dP_dPsfc_M)
      status = vgd_dpidpis(statevector%vco%vgrid, &
                           statevector%vco%ip1_M, &
                           dP_dPsfc_M, &
                           Psfc)
      if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_dpidpis')
      ! calculate delP_M
      do lev_M = 1, nlev_M
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do lonIndex = statevector%myLonBeg, statevector%myLonEnd
            delPsfc(lonIndex,latIndex,1,stepIndex) = delPsfc(lonIndex,latIndex,1,stepIndex) + dP_dPsfc_M(lonIndex,latIndex,lev_M) * delP_M(lonIndex,latIndex,lev_M,stepIndex)
          end do
        end do
      end do
      deallocate(dP_dPsfc_M)

      ! dP_dPsfc_T
      nullify(dP_dPsfc_T)
      status = vgd_dpidpis(statevector%vco%vgrid, &
                           statevector%vco%ip1_T, &
                           dP_dPsfc_T, &
                           Psfc)
      if( status .ne. VGD_OK ) call utl_abort('ERROR with vgd_dpidpis')
      ! calculate delP_T
      do lev_T = 1, nlev_T
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do lonIndex = statevector%myLonBeg, statevector%myLonEnd
            delPsfc(lonIndex,latIndex,1,stepIndex) = delPsfc(lonIndex,latIndex,1,stepIndex) + dP_dPsfc_T(lonIndex,latIndex,lev_T) * delP_T(lonIndex,latIndex,lev_T,stepIndex)
          end do
        end do
      end do
      deallocate(dP_dPsfc_T)

      if ( .not. beSilent .and. stepIndex == 1 ) then
        write(*,*) 'stepIndex=',stepIndex, ',delPsfc='
        write(*,*) delPsfc(statevector%myLonBeg,statevector%myLatBeg,1,stepIndex)
      end if

    end do

    deallocate(Psfc)

  end subroutine clp_calcPressure_ad


end module calcPressure_mod

