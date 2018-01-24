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
!! MODULE tt2phi (prefix="tt2phi")
!!
!! *Purpose*: Subroutines for computing GZ from TT, HU and P0. Nonlinear, tangent-
!!            linear and adjoint versions of this transformation are included in separate
!!            subroutines.
!!
!--------------------------------------------------------------------------
module tt2phi_mod
  use mpi_mod
  use mathPhysConstants_mod
  use physicsFunctions_mod
  use earthConstants_mod
  use columnData_mod
  use utilities_mod
  implicit none
  save
  private

  ! public procedures
  public :: tt2phi, tt2phi_tl, tt2phi_ad

contains

subroutine tt2phi(columnghr,beSilent_opt)
  !
  !**s/r tt2phi - Temperature to geopotential transformation on GEM4 staggered levels
  !               NOTE: we assume 
  !                     1) nlev_T = nlev_M+1 
  !                     2) GZ_T(nlev_T) = GZ_M(nlev_M), both at the surface
  !                     3) a thermo level exists at the top, higher than the highest momentum level
  !                     4) the placement of the thermo levels means that GZ_T is the average of 2 nearest GZ_M
  !                        (according to Ron and Claude)
  !
  !Author  : M. Buehner, February 2014
  !
  implicit none

  type(struct_columnData) :: columnghr
  logical, optional       :: beSilent_opt

  integer :: columnIndex,lev_M,lev_T,nlev_M,nlev_T,status,Vcode
  real(8) :: hu,tt,ratioP
  real(8), allocatable :: tv(:)
  real(8), pointer     :: gz_T(:),gz_M(:)
  real                 :: gz_sfcOffset_T_r4, gz_sfcOffset_M_r4
  type(struct_vco), pointer :: vco_ghr
  logical                   :: beSilent

  if ( present(beSilent_opt) ) then
    beSilent = beSilent_opt
  else
    beSilent = .false.
  end if

  vco_ghr => col_getVco(columnghr)
  status = vgd_get(vco_ghr%vgrid,key='ig_1 - vertical coord code',value=Vcode)

  nlev_T = col_getNumLev(columnghr,'TH')
  nlev_M = col_getNumLev(columnghr,'MM')
  if(Vcode .eq. 5002 .and. nlev_T .ne. nlev_M+1) call utl_abort('tt2phi: nlev_T is not equal to nlev_M+1!')
  if(Vcode .eq. 5005 .and. nlev_T .ne. nlev_M)   call utl_abort('tt2phi: nlev_T is not equal to nlev_M!')

  allocate(tv(nlev_T))

  if(Vcode .eq. 5002) then

    ! loop over all columns
    do columnIndex = 1, col_getNumCol(columnghr)

      gz_M => col_getColumn(columnghr,columnIndex,'GZ','MM')
      gz_T => col_getColumn(columnghr,columnIndex,'GZ','TH')

      ! set the surface height
      gz_M(nlev_M) = col_getGZsfc(columnghr,columnIndex)
      gz_T(nlev_T) = col_getGZsfc(columnghr,columnIndex)

      ! initialize the rest to zero
      gz_M(1:(nlev_M-1)) = 0.0d0
      gz_T(1:(nlev_T-1)) = 0.0d0

      ! compute virtual temperature on thermo levels
      do lev_T = 1, nlev_T
        hu = exp(col_getElem(columnghr,lev_T,columnIndex,'HU'))
        tt = col_getElem(columnghr,lev_T,columnIndex,'TT')
        tv(lev_T) = fotvt8(tt,hu)
      enddo
    
      ! compute GZ on momentum levels
      do lev_M = (nlev_M-1), 1, -1
        lev_T = lev_M+1 ! thermo level just below momentum level being computed
        if(col_getPressure(columnghr,lev_M,columnIndex,'MM').eq.0.0d0) then
          write(*,*) 'tt2phi: pressure is zero, lev_m, columnIndex=',lev_m, columnIndex
          call utl_abort('tt2phi')
        endif
        ratioP = col_getPressure(columnghr,lev_M+1,columnIndex,'MM') / col_getPressure(columnghr,lev_M,columnIndex,'MM')
        gz_M(lev_M) = gz_M(lev_M+1) + MPC_RGAS_DRY_AIR_R8*tv(lev_T)*log(ratioP)
      enddo

      ! compute GZ on top thermo level (from top momentum level)
      ratioP = col_getPressure(columnghr,1,columnIndex,'MM') / col_getPressure(columnghr,1,columnIndex,'TH')
      gz_T(1) = gz_M(1) + MPC_RGAS_DRY_AIR_R8*tv(1)*log(ratioP)

      ! compute GZ on remaining thermo levels by simple averaging
      do lev_T = 2, (nlev_T-1)
        lev_M = lev_T ! momentum level just below thermo level being computed
        gz_T(lev_T) = 0.5d0*( gz_M(lev_M-1) + gz_M(lev_M) )
      enddo

      !if(columnIndex.eq.1) then
      !  do lev_M = 1, nlev_M
      !    write(*,*) 'tt2phi: pres_M,gz_M=',lev_M,col_getPressure(columnghr,lev_M,columnIndex,'MM'),gz_M(lev_M)
      !  enddo
      !  do lev_T = 1, nlev_T
      !    write(*,*) 'tt2phi: pres_T,gz_T=',lev_T,col_getPressure(columnghr,lev_T,columnIndex,'TH'),gz_T(lev_T)
      !  enddo
      !endif

    enddo

  elseif(Vcode .eq. 5005) then

    status = vgd_get(columnghr%vco%vgrid,key='DHM - height of the diagnostic level (m)',value=gz_sfcOffset_M_r4)
    status = vgd_get(columnghr%vco%vgrid,key='DHT - height of the diagnostic level (t)',value=gz_sfcOffset_T_r4)
    if(mpi_myid == 0 .and. .not.beSilent ) then
      write(*,*) 'col_fillmvo: height offset for near-sfc momentum level is: ', gz_sfcOffset_M_r4, ' metres'
      write(*,*) 'col_fillmvo: height offset for near-sfc thermo level is:   ', gz_sfcOffset_T_r4, ' metres'
    end if

    ! loop over all columns
    do columnIndex = 1, col_getNumCol(columnghr)

      gz_M => col_getColumn(columnghr,columnIndex,'GZ','MM')
      gz_T => col_getColumn(columnghr,columnIndex,'GZ','TH')

      ! set the surface height (this is the true surface, not the lowest UU/TT level)
      gz_M(nlev_M) = col_getGZsfc(columnghr,columnIndex) + real(gz_sfcOffset_M_r4,8) * RG
      gz_T(nlev_T) = col_getGZsfc(columnghr,columnIndex) + real(gz_sfcOffset_T_r4,8) * RG

      ! initialize the rest to zero
      gz_M(1:(nlev_M-1)) = 0.0d0
      gz_T(1:(nlev_T-1)) = 0.0d0

      ! compute virtual temperature on thermo levels
      do lev_T = 1, nlev_T
        hu = exp(col_getElem(columnghr,lev_T,columnIndex,'HU'))
        tt = col_getElem(columnghr,lev_T,columnIndex,'TT')
        tv(lev_T) = fotvt8(tt,hu)
      enddo
    
      ! compute GZ on momentum levels
      do lev_M = (nlev_M-1), 1, -1
        lev_T = lev_M ! thermo level just below momentum level being computed
        if(col_getPressure(columnghr,lev_M,columnIndex,'MM').eq.0.0d0) then
          write(*,*) 'tt2phi: pressure is zero, lev_m, columnIndex=',lev_m, columnIndex
          call utl_abort('tt2phi')
        endif
        ratioP = col_getPressure(columnghr,lev_M+1,columnIndex,'MM') / col_getPressure(columnghr,lev_M,columnIndex,'MM')
        gz_M(lev_M) = gz_M(lev_M+1) + MPC_RGAS_DRY_AIR_R8*tv(lev_T)*log(ratioP)
      enddo

      ! compute GZ on thermo levels by simple averaging
      do lev_T = 1, (nlev_T-1)
        lev_M = lev_T+1 ! momentum level just below thermo level being computed
        gz_T(lev_T) = 0.5d0*( gz_M(lev_M-1) + gz_M(lev_M) )
      enddo

      !if(columnIndex.eq.1) then
      !  do lev_M = 1, nlev_M
      !    write(*,*) 'tt2phi: pres_M,gz_M=',lev_M,col_getPressure(columnghr,lev_M,columnIndex,'MM'),gz_M(lev_M)
      !  enddo
      !  do lev_T = 1, nlev_T
      !    write(*,*) 'tt2phi: pres_T,gz_T=',lev_T,col_getPressure(columnghr,lev_T,columnIndex,'TH'),gz_T(lev_T)
      !  enddo
      !endif

    enddo

  endif

  deallocate(tv)

end subroutine tt2phi


subroutine tt2phi_tl(column,columng)
  !
  !**s/r tt2phi_tl - Temperature to geopotential transformation on GEM4 staggered levels
  !               NOTE: we assume 
  !                     1) nlev_T = nlev_M+1 
  !                     2) GZ_T(nlev_T) = GZ_M(nlev_M), both at the surface
  !                     3) a thermo level exists at the top, higher than the highest momentum level
  !                     4) the placement of the thermo levels means that GZ_T is the average of 2 nearest GZ_M
  !                        (according to Ron and Claude)
  !
  ! NOTE: after revision 680 removed code for vcode=5005 when rewriting for increased efficiency (M. Buehner)
  !
  !Author  : M. Buehner, February 2014
  !
  implicit none

  type(struct_columnData) :: column,columng

  integer :: columnIndex,lev_M,lev_T,nlev_M,nlev_T,status,Vcode_anl
  real(8) :: hu,tt,ratioP1,delLnP_M1,delLnP_T1
  real(8), allocatable :: ratioP(:), delThick(:)
  real(8), allocatable :: delLnP_M(:),delLnP_T(:)
  real(8), pointer     :: delGz_M(:),delGz_T(:),delTT(:),delLQ(:),delP0(:)
  type(struct_vco), pointer :: vco_anl

  real(8), allocatable, save :: coeff_M_TT(:,:), coeff_M_HU(:,:), coeff_M_P0(:,:), &
                                coeff_T_TT(:), coeff_T_HU(:), coeff_T_P0(:)
  logical, save :: firstTime = .true.

  vco_anl => col_getVco(columng)
  status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)

  nlev_T = col_getNumLev(columng,'TH')
  nlev_M = col_getNumLev(columng,'MM')

  allocate(ratioP(nlev_T))
  allocate(delThick(nlev_T))
  allocate(delLnP_M(nlev_M))
  allocate(delLnP_T(nlev_T))

  if(Vcode_anl .eq. 5002) then

    if(firstTime) then
      ! initialize and save coefficients for increased efficiency (assumes no relinearization)
      firstTime = .false.

      allocate(coeff_M_TT(nlev_T,col_getNumCol(columng)))
      allocate(coeff_M_HU(nlev_T,col_getNumCol(columng)))
      allocate(coeff_M_P0(nlev_T,col_getNumCol(columng)))
      allocate(coeff_T_TT(col_getNumCol(columng)))
      allocate(coeff_T_HU(col_getNumCol(columng)))
      allocate(coeff_T_P0(col_getNumCol(columng)))

      do columnIndex = 1, col_getNumCol(columng)

        ! compute coefficients for GZ on momentum levels
        do lev_M = 1, nlev_M
          delLnP_M(lev_M) = col_getPressureDeriv(columng,lev_M,columnIndex,'MM')/  &
                            col_getPressure(columng,lev_M,columnIndex,'MM')
        enddo
        do lev_T = 2, (nlev_T-1)
          ratioP1 = log( col_getPressure(columng,lev_T  ,columnIndex,'MM') /  &
                         col_getPressure(columng,lev_T-1,columnIndex,'MM') )
          hu = exp(col_getElem(columng,lev_T,columnIndex,'HU'))
          tt = col_getElem(columng,lev_T,columnIndex,'TT')
          coeff_M_TT(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * fottva(hu,1.0d0)
          coeff_M_HU(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * folnqva(hu,tt,1.0d0)
          coeff_M_P0(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * (delLnP_M(lev_T)-delLnP_M(lev_T-1)) * fotvt8(tt,hu)
        enddo

        ! compute coefficients for GZ on only the top thermo level
        ratioP1 = log( col_getPressure(columng,1,columnIndex,'MM') /  &
                       col_getPressure(columng,1,columnIndex,'TH') )
        delLnP_M1 = col_getPressureDeriv(columng,1,columnIndex,'MM')/  &
                    col_getPressure(columng,1,columnIndex,'MM')
        delLnP_T1 = col_getPressureDeriv(columng,1,columnIndex,'TH')/  &
                    col_getPressure(columng,1,columnIndex,'TH')
        hu = exp(col_getElem(columng,1,columnIndex,'HU'))
        tt = col_getElem(columng,1,columnIndex,'TT')
        coeff_T_TT(columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * fottva(hu,1.0d0)
        coeff_T_HU(columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * folnqva(hu,tt,1.0d0)
        coeff_T_P0(columnIndex) = MPC_RGAS_DRY_AIR_R8 * (delLnP_M1-delLnP_T1) * fotvt8(tt,hu)

      enddo

    endif

    ! loop over all columns

!$OMP PARALLEL DO PRIVATE(columnIndex,delGz_M,delGz_T,delThick,delTT,delLQ,delP0,lev_M,lev_T)
    do columnIndex = 1, col_getNumCol(columng)

      delGz_M => col_getColumn(column,columnIndex,'GZ','MM')
      delGz_T => col_getColumn(column,columnIndex,'GZ','TH')
      delTT   => col_getColumn(column,columnIndex,'TT')
      delLQ   => col_getColumn(column,columnIndex,'HU')
      delP0   => col_getColumn(column,columnIndex,'P0')

      ! ensure increment at sfc is zero (fixed height level)
      delGz_M(nlev_M) = 0.0d0
      delGz_T(nlev_T) = 0.0d0

      ! compute increment to thickness for each layer
      do lev_T = 2, (nlev_T-1)
        delThick(lev_T) = coeff_M_TT(lev_T,columnIndex) * delTT(lev_T) + &
                          coeff_M_HU(lev_T,columnIndex) * delLQ(lev_T) + &
                          coeff_M_P0(lev_T,columnIndex) * delP0(1)
      enddo

      ! compute GZ increment on momentum levels above the surface
      do lev_M = (nlev_M-1), 1, -1
        lev_T = lev_M+1 ! thermo level just below momentum level being computed
        delGz_M(lev_M) = delGz_M(lev_M+1) + delThick(lev_T) 
      enddo

      ! compute GZ increment on thermo levels by simple averaging, except for top level
      do lev_T = 2, (nlev_T-1)
        lev_M = lev_T ! momentum level just below thermo level being computed
        delGz_T(lev_T) = 0.5d0*( delGz_M(lev_M-1) + delGz_M(lev_M) )
      enddo

      ! compute GZ increment for top thermo level (from top momentum level)
      delGz_T(1) = delGz_M(1) +  &
                   coeff_T_TT(columnIndex) * delTT(1) + &
                   coeff_T_HU(columnIndex) * delLQ(1) + &
                   coeff_T_P0(columnIndex) * delP0(1)

      !if(columnIndex.eq.1) then
      !  do lev_M = 1, nlev_M
      !    write(*,*) 'tt2phi_tl: delGz_M=',lev_M,delGz_M(lev_M)
      !  enddo
      !  do lev_T = 1, nlev_T
      !    write(*,*) 'tt2phi_tl: delGz_T=',lev_T,delGz_T(lev_T)
      !  enddo
      !endif

    enddo
!$OMP END PARALLEL DO

  elseif(Vcode_anl .eq. 5005) then

    if(firstTime) then
      ! initialize and save coefficients for increased efficiency (assumes no relinearization)
      firstTime = .false.

      allocate(coeff_M_TT(nlev_T,col_getNumCol(columng)))
      allocate(coeff_M_HU(nlev_T,col_getNumCol(columng)))
      allocate(coeff_M_P0(nlev_T,col_getNumCol(columng)))

      do columnIndex = 1, col_getNumCol(columng)

        ! compute coefficients for GZ on momentum levels
        do lev_M = 1, nlev_M
          delLnP_M(lev_M) = col_getPressureDeriv(columng,lev_M,columnIndex,'MM')/  &
                            col_getPressure(columng,lev_M,columnIndex,'MM')
        enddo
        do lev_T = 1, (nlev_T-1)
          ratioP1 = log( col_getPressure(columng,lev_T+1,columnIndex,'MM') /  &
                         col_getPressure(columng,lev_T  ,columnIndex,'MM') )
          hu = exp(col_getElem(columng,lev_T,columnIndex,'HU'))
          tt = col_getElem(columng,lev_T,columnIndex,'TT')
          coeff_M_TT(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * fottva(hu,1.0d0)
          coeff_M_HU(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * folnqva(hu,tt,1.0d0)
          coeff_M_P0(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * (delLnP_M(lev_T+1)-delLnP_M(lev_T)) * fotvt8(tt,hu)
        enddo

      enddo

    endif

    ! loop over all columns

!$OMP PARALLEL DO PRIVATE(columnIndex,delGz_M,delGz_T,delThick,delTT,delLQ,delP0,lev_M,lev_T)
    do columnIndex = 1, col_getNumCol(columng)

      delGz_M => col_getColumn(column,columnIndex,'GZ','MM')
      delGz_T => col_getColumn(column,columnIndex,'GZ','TH')
      delTT   => col_getColumn(column,columnIndex,'TT')
      delLQ   => col_getColumn(column,columnIndex,'HU')
      delP0   => col_getColumn(column,columnIndex,'P0')

      ! ensure increment at sfc is zero (fixed height level)
      delGz_M(nlev_M) = 0.0d0
      delGz_T(nlev_T) = 0.0d0

      ! compute increment to thickness for each layer
      do lev_T = 1, (nlev_T-1)
        delThick(lev_T) = coeff_M_TT(lev_T,columnIndex) * delTT(lev_T) + &
                          coeff_M_HU(lev_T,columnIndex) * delLQ(lev_T) + &
                          coeff_M_P0(lev_T,columnIndex) * delP0(1)
      enddo

      ! compute GZ increment on momentum levels above the surface
      do lev_M = (nlev_M-1), 1, -1
        lev_T = lev_M ! thermo level just below momentum level being computed
        delGz_M(lev_M) = delGz_M(lev_M+1) + delThick(lev_T) 
      enddo

      ! compute GZ increment on thermo levels by simple averaging, except for top level
      do lev_T = 1, (nlev_T-1)
        lev_M = lev_T+1 ! momentum level just below thermo level being computed
        delGz_T(lev_T) = 0.5d0*( delGz_M(lev_M-1) + delGz_M(lev_M) )
      enddo

      !if(columnIndex.eq.1) then
      !  do lev_M = 1, nlev_M
      !    write(*,*) 'tt2phi_tl: delGz_M=',lev_M,delGz_M(lev_M)
      !  enddo
      !  do lev_T = 1, nlev_T
      !    write(*,*) 'tt2phi_tl: delGz_T=',lev_T,delGz_T(lev_T)
      !  enddo
      !endif

    enddo
!$OMP END PARALLEL DO

  endif

  deallocate(ratioP)
  deallocate(delThick)
  deallocate(delLnP_M)
  deallocate(delLnP_T)

end subroutine tt2phi_tl


subroutine tt2phi_ad(column,columng)
  !
  !**s/r tt2phi_ad- Adjoint of temperature to geopotential transformation on GEM4 staggered levels
  !               NOTE: we assume 
  !                     1) nlev_T = nlev_M+1 
  !                     2) GZ_T(nlev_T) = GZ_M(nlev_M), both at the surface
  !                     3) a thermo level exists at the top, higher than the highest momentum level
  !                     4) the placement of the thermo levels means that GZ_T is the average of 2 nearest GZ_M
  !                        (according to Ron and Claude)
  !
  ! NOTE: after revision 680 removed code for vcode=5005 when rewriting for increased efficiency (M. Buehner)
  !
  !Author  : M. Buehner, February 2014
  !
  implicit none

  type(struct_columnData) :: column,columng

  integer :: columnIndex,lev_M,lev_T,nlev_M,nlev_T,status,Vcode_anl
  real(8) :: hu,tt,ratioP1,delLnP_M1,delLnP_T1
  real(8), allocatable :: ratioP(:),sumGz_T(:)
  real(8), allocatable :: delLnP_M(:),delLnP_T(:),delGz_M(:),delGz_T(:)
  real(8), pointer     :: delGz_M_in(:),delGz_T_in(:),delTT(:),delLQ(:),delP0(:)
  type(struct_vco), pointer :: vco_anl

  real(8), allocatable, save :: coeff_M_TT(:,:), coeff_M_HU(:,:), coeff_M_P0(:,:), &
                                coeff_T_TT(:), coeff_T_HU(:), coeff_T_P0(:)
  logical, save :: firstTime = .true.

  vco_anl => col_getVco(columng)
  status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)

  nlev_T = col_getNumLev(columng,'TH')
  nlev_M = col_getNumLev(columng,'MM')

  allocate(ratioP(nlev_M))
  allocate(delLnP_M(nlev_M))
  allocate(delLnP_T(nlev_T))
  allocate(delGz_M(nlev_M))
  allocate(delGz_T(nlev_T))
  allocate(sumGz_T(0:nlev_T))

  if(Vcode_anl .eq. 5002) then

    if(firstTime) then
      ! initialize and save coefficients for increased efficiency (assumes no relinearization)
      firstTime = .false.

      allocate(coeff_M_TT(nlev_T,col_getNumCol(columng)))
      allocate(coeff_M_HU(nlev_T,col_getNumCol(columng)))
      allocate(coeff_M_P0(nlev_T,col_getNumCol(columng)))
      allocate(coeff_T_TT(col_getNumCol(columng)))
      allocate(coeff_T_HU(col_getNumCol(columng)))
      allocate(coeff_T_P0(col_getNumCol(columng)))

      do columnIndex = 1, col_getNumCol(columng)

        ! compute coefficients for GZ on momentum levels
        do lev_M = 1, nlev_M
          delLnP_M(lev_M) = col_getPressureDeriv(columng,lev_M,columnIndex,'MM')/  &
                            col_getPressure(columng,lev_M,columnIndex,'MM')
        enddo
        do lev_T = 2, (nlev_T-1)
          ratioP1 = log( col_getPressure(columng,lev_T  ,columnIndex,'MM') /  &
                         col_getPressure(columng,lev_T-1,columnIndex,'MM') )
          hu = exp(col_getElem(columng,lev_T,columnIndex,'HU'))
          tt = col_getElem(columng,lev_T,columnIndex,'TT')
          coeff_M_TT(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * fottva(hu,1.0d0)
          coeff_M_HU(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * folnqva(hu,tt,1.0d0)
          coeff_M_P0(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * (delLnP_M(lev_T)-delLnP_M(lev_T-1)) * fotvt8(tt,hu)
        enddo

        ! compute coefficients for GZ on only the top thermo level
        ratioP1 = log( col_getPressure(columng,1,columnIndex,'MM') /  &
                       col_getPressure(columng,1,columnIndex,'TH') )
        delLnP_M1 = col_getPressureDeriv(columng,1,columnIndex,'MM')/  &
                    col_getPressure(columng,1,columnIndex,'MM')
        delLnP_T1 = col_getPressureDeriv(columng,1,columnIndex,'TH')/  &
                    col_getPressure(columng,1,columnIndex,'TH')
        hu = exp(col_getElem(columng,1,columnIndex,'HU'))
        tt = col_getElem(columng,1,columnIndex,'TT')
        coeff_T_TT(columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * fottva(hu,1.0d0)
        coeff_T_HU(columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * folnqva(hu,tt,1.0d0)
        coeff_T_P0(columnIndex) = MPC_RGAS_DRY_AIR_R8 * (delLnP_M1-delLnP_T1) * fotvt8(tt,hu)

      enddo

    endif

    ! loop over all columns

!$OMP PARALLEL DO PRIVATE(columnIndex,delGz_M_in,delGz_T_in,delTT,  &
!$OMP delLQ,delP0,lev_M,lev_T,sumGz_T,delGz_M,delGz_T)
    do columnIndex = 1, col_getNumCol(columng)

      delGz_M_in => col_getColumn(column,columnIndex,'GZ','MM')
      delGz_T_in => col_getColumn(column,columnIndex,'GZ','TH')
      delTT      => col_getColumn(column,columnIndex,'TT')
      delLQ      => col_getColumn(column,columnIndex,'HU')
      delP0      => col_getColumn(column,columnIndex,'P0')

      !if(columnIndex.eq.1) then
      !  do lev_M = 1, nlev_M
      !    write(*,*) 'tt2phi_ad: gradient wrt GZ_M=',lev_M,delGZ_M_in(lev_M)
      !  enddo
      !  do lev_T = 1, nlev_T
      !    write(*,*) 'tt2phi_ad: gradient wrt GZ_T=',lev_T,delGZ_T_in(lev_T)
      !  enddo
      !endif

      delGz_M(:) = delGz_M_in(:)
      delGz_T(:) = delGz_T_in(:)

      ! adjoint of compute GZ increment on remaining thermo levels by simple averaging
      do lev_T = 2, (nlev_T-1)
        lev_M = lev_T ! momentum level just below thermo level being computed
        !delGz_T(lev_T) = 0.5d0*( delGz_M(lev_M-1) + delGz_M(lev_M) )
        delGz_M(lev_M-1) = delGz_M(lev_M-1) + 0.5d0*delGz_T(lev_T)
        delGz_M(lev_M)   = delGz_M(lev_M)   + 0.5d0*delGz_T(lev_T)
      enddo

      ! adjoint of compute GZ increment on top thermo level (from top momentum level)
      delGz_M(1)  = delGz_M(1)  + delGz_T(1)
      delTT(1) = delTT(1) + coeff_T_TT(columnIndex)*delGz_T(1)
      delLQ(1) = delLQ(1) + coeff_T_HU(columnIndex)*delGz_T(1)
      delP0(1) = delP0(1) + coeff_T_P0(columnIndex)*delGz_T(1)

      ! adjoint of compute GZ increment on momentum levels
      sumGz_T(0:1) = 0.0d0
      do lev_M = 1, (nlev_M-1)
        lev_T = lev_M+1 ! thermo level just below momentum level being computed
        sumGz_T(lev_T) = sumGz_T(lev_T-1) + delGz_M(lev_M)
      enddo
      do lev_T = 2, nlev_T-1
        delTT(lev_T) = delTT(lev_T) + coeff_M_TT(lev_T,columnIndex)*sumGz_T(lev_T)
        delLQ(lev_T) = delLQ(lev_T) + coeff_M_HU(lev_T,columnIndex)*sumGz_T(lev_T)
        delP0(1)     = delP0(1)     + coeff_M_P0(lev_T,columnIndex)*sumGz_T(lev_T)
      enddo

      !if(columnIndex.eq.1) then
      !  do lev_T = 1, nlev_T
      !    write(*,*) 'tt2phi_ad: gradient wrt TT, LQ=',lev_T,delTT(lev_T),delLQ(lev_T)
      !  enddo
      !  write(*,*) 'tt2phi_ad: gradient wrt P0=',delP0(1)
      !endif

    enddo
!$OMP END PARALLEL DO

  elseif(Vcode_anl .eq. 5005) then

    if(firstTime) then
      ! initialize and save coefficients for increased efficiency (assumes no relinearization)
      firstTime = .false.

      allocate(coeff_M_TT(nlev_T,col_getNumCol(columng)))
      allocate(coeff_M_HU(nlev_T,col_getNumCol(columng)))
      allocate(coeff_M_P0(nlev_T,col_getNumCol(columng)))

      do columnIndex = 1, col_getNumCol(columng)

        ! compute coefficients for GZ on momentum levels
        do lev_M = 1, nlev_M
          delLnP_M(lev_M) = col_getPressureDeriv(columng,lev_M,columnIndex,'MM')/  &
                            col_getPressure(columng,lev_M,columnIndex,'MM')
        enddo
        do lev_T = 1, (nlev_T-1)
          ratioP1 = log( col_getPressure(columng,lev_T+1,columnIndex,'MM') /  &
                         col_getPressure(columng,lev_T  ,columnIndex,'MM') )
          hu = exp(col_getElem(columng,lev_T,columnIndex,'HU'))
          tt = col_getElem(columng,lev_T,columnIndex,'TT')
          coeff_M_TT(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * fottva(hu,1.0d0)
          coeff_M_HU(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * folnqva(hu,tt,1.0d0)
          coeff_M_P0(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * (delLnP_M(lev_T+1)-delLnP_M(lev_T)) * fotvt8(tt,hu)
        enddo

      enddo

    endif

    ! loop over all columns

!$OMP PARALLEL DO PRIVATE(columnIndex,delGz_M_in,delGz_T_in,delTT,  &
!$OMP delLQ,delP0,lev_M,lev_T,sumGz_T,delGz_M,delGz_T)
    do columnIndex = 1, col_getNumCol(columng)

      delGz_M_in => col_getColumn(column,columnIndex,'GZ','MM')
      delGz_T_in => col_getColumn(column,columnIndex,'GZ','TH')
      delTT      => col_getColumn(column,columnIndex,'TT')
      delLQ      => col_getColumn(column,columnIndex,'HU')
      delP0      => col_getColumn(column,columnIndex,'P0')

      !if(columnIndex.eq.1) then
      !  do lev_M = 1, nlev_M
      !    write(*,*) 'tt2phi_ad: gradient wrt GZ_M=',lev_M,delGZ_M_in(lev_M)
      !  enddo
      !  do lev_T = 1, nlev_T
      !    write(*,*) 'tt2phi_ad: gradient wrt GZ_T=',lev_T,delGZ_T_in(lev_T)
      !  enddo
      !endif

      delGz_M(:) = delGz_M_in(:)
      delGz_T(:) = delGz_T_in(:)

      ! adjoint of compute GZ increment on remaining thermo levels by simple averaging
      do lev_T = 1, (nlev_T-1)
        lev_M = lev_T+1 ! momentum level just below thermo level being computed
        !delGz_T(lev_T) = 0.5d0*( delGz_M(lev_M-1) + delGz_M(lev_M) )
        delGz_M(lev_M-1) = delGz_M(lev_M-1) + 0.5d0*delGz_T(lev_T)
        delGz_M(lev_M)   = delGz_M(lev_M)   + 0.5d0*delGz_T(lev_T)
      enddo

      ! adjoint of compute GZ increment on momentum levels
      sumGz_T(0) = 0.0d0
      do lev_M = 1, (nlev_M-1)
        lev_T = lev_M ! thermo level just below momentum level being computed
        sumGz_T(lev_T) = sumGz_T(lev_T-1) + delGz_M(lev_M)
      enddo
      do lev_T = 1, nlev_T-1
        delTT(lev_T) = delTT(lev_T) + coeff_M_TT(lev_T,columnIndex)*sumGz_T(lev_T)
        delLQ(lev_T) = delLQ(lev_T) + coeff_M_HU(lev_T,columnIndex)*sumGz_T(lev_T)
        delP0(1)     = delP0(1)     + coeff_M_P0(lev_T,columnIndex)*sumGz_T(lev_T)
      enddo

      !if(columnIndex.eq.1) then
      !  do lev_T = 1, nlev_T
      !    write(*,*) 'tt2phi_ad: gradient wrt TT, LQ=',lev_T,delTT(lev_T),delLQ(lev_T)
      !  enddo
      !  write(*,*) 'tt2phi_ad: gradient wrt P0=',delP0(1)
      !endif

    enddo
!$OMP END PARALLEL DO

  endif

  deallocate(ratioP)
  deallocate(delLnP_M)
  deallocate(delLnP_T)
  deallocate(delGz_M)
  deallocate(delGz_T)

end subroutine tt2phi_ad


end module tt2phi_mod