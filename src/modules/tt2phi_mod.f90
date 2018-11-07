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
!! MODULE tt2phi (prefix="tt2phi" category='3. High-level transformations')
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
  use verticalCoord_mod
  use utilities_mod
  implicit none
  save
  private

  ! public procedures
  public :: tt2phi           , tt2phi_tl           , tt2phi_ad
  public :: tt2phi_gpsro_tl     
  public :: tt2phi_gpsro_1Col, tt2phi_gpsro_tl_1Col
  public :: tt2phi_gpsro_ad

  ! constants from gps_mod {{{
  real(8), parameter :: delta = 0.6077686814144D0

  ! Avogadro constant:
  real(8), parameter :: p_Avog = 6.02214129D23        ! From CODATA

  ! Boltzmann constant:
  real(8), parameter :: p_Boltz = 1.3806488D-23        ! From CODATA

  ! Air properties:
  real(8), parameter :: p_md = 28.965516D0            ! From Aparicio(2011)
  real(8), parameter :: p_mw = 18.015254D0            ! From Aparicio(2011)
  real(8), parameter :: p_wa = p_md/p_mw
  real(8), parameter :: p_wb = (p_md-p_mw)/p_mw

  real(8), parameter :: p_Rd = p_Avog*p_Boltz/(1.D-3*p_md)   ! per air mass

  ! Angular velocity of the Earth (omegaPrime) (radians/s).
  ! Standard Earth, rotating with a constant angular velocity (IAU, GRS67).
  real(8), parameter :: WGS_OmegaPrime = 7292115.1467D-11

  ! Units and scales:
  real(8), parameter :: p_TC    = 273.15D0
  real(8), parameter :: p_knot  = 0.514444D0

  ! Semimajor axis (a) (m)                             [*Defining constant*]
  real(8), parameter :: WGS_a = 6378137.0D0

  ! Theoretical (Normal) Gravity Formula Constant:
  real(8), parameter :: WGS_TNGk = 0.00193185265241D0

  ! First eccentricity squared:
  real(8), parameter :: WGS_e2 = 6.69437999014D-3

  ! Theoretical (Normal) Gravity at the equator (m/s2):
  real(8), parameter :: WGS_GammaE = 9.7803253359D0

  ! Flattening (f)                                     [*Defining constant*]
  real(8), parameter :: WGS_f = 1.D0 / 298.257223563D0

  ! m = omega^2 a^2 b / GM
  real(8), parameter :: WGS_m = 0.00344978650684D0
  ! }}}

contains
!
!subroutine tt2phi(columnghr,beSilent_opt) !{{{
!  !
!  !**s/r tt2phi - Temperature to geopotential transformation on GEM4 staggered levels
!  !               NOTE: we assume 
!  !                     1) nlev_T = nlev_M+1 
!  !                     2) GZ_T(nlev_T) = GZ_M(nlev_M), both at the surface
!  !                     3) a thermo level exists at the top, higher than the highest momentum level
!  !                     4) the placement of the thermo levels means that GZ_T is the average of 2 nearest GZ_M
!  !                        (according to Ron and Claude)
!  !
!  !Author  : M. Buehner, February 2014
!  !
!  implicit none
!
!  type(struct_columnData) :: columnghr
!  logical, optional       :: beSilent_opt
!
!  integer :: columnIndex,lev_M,lev_T,nlev_M,nlev_T,status,Vcode
!  real(8) :: hu,tt,ratioP
!  real(8), allocatable :: tv(:)
!  real(8), pointer     :: gz_T(:),gz_M(:)
!  real                 :: gz_sfcOffset_T_r4, gz_sfcOffset_M_r4
!  type(struct_vco), pointer :: vco_ghr
!  logical                   :: beSilent
!
!  if ( present(beSilent_opt) ) then
!    beSilent = beSilent_opt
!  else
!    beSilent = .false.
!  end if
!
!  vco_ghr => col_getVco(columnghr)
!  status = vgd_get(vco_ghr%vgrid,key='ig_1 - vertical coord code',value=Vcode)
!
!  nlev_T = col_getNumLev(columnghr,'TH')
!  nlev_M = col_getNumLev(columnghr,'MM')
!  if(Vcode .eq. 5002 .and. nlev_T .ne. nlev_M+1) call utl_abort('tt2phi: nlev_T is not equal to nlev_M+1!')
!  if(Vcode .eq. 5005 .and. nlev_T .ne. nlev_M)   call utl_abort('tt2phi: nlev_T is not equal to nlev_M!')
!
!  allocate(tv(nlev_T))
!
!  write(*,*) 'MAZIAR: tt2phi, Vcode=', Vcode
!
!  if(Vcode .eq. 5002) then
!
!    ! loop over all columns
!    do columnIndex = 1, col_getNumCol(columnghr)
!
!      gz_M => col_getColumn(columnghr,columnIndex,'GZ','MM')
!      gz_T => col_getColumn(columnghr,columnIndex,'GZ','TH')
!
!      ! set the surface height
!      gz_M(nlev_M) = col_getGZsfc(columnghr,columnIndex)
!      gz_T(nlev_T) = col_getGZsfc(columnghr,columnIndex)
!
!      ! initialize the rest to zero
!      gz_M(1:(nlev_M-1)) = 0.0d0
!      gz_T(1:(nlev_T-1)) = 0.0d0
!
!      ! compute virtual temperature on thermo levels
!      do lev_T = 1, nlev_T
!        hu = col_getElem(columnghr,lev_T,columnIndex,'HU')
!        tt = col_getElem(columnghr,lev_T,columnIndex,'TT')
!        tv(lev_T) = fotvt8(tt,hu)
!      enddo
!    
!      ! compute GZ on momentum levels
!      do lev_M = (nlev_M-1), 1, -1
!        lev_T = lev_M+1 ! thermo level just below momentum level being computed
!        if(col_getPressure(columnghr,lev_M,columnIndex,'MM').eq.0.0d0) then
!          write(*,*) 'tt2phi: pressure is zero, lev_m, columnIndex=',lev_m, columnIndex
!          call utl_abort('tt2phi')
!        endif
!        ratioP = col_getPressure(columnghr,lev_M+1,columnIndex,'MM') / col_getPressure(columnghr,lev_M,columnIndex,'MM')
!        gz_M(lev_M) = gz_M(lev_M+1) + MPC_RGAS_DRY_AIR_R8*tv(lev_T)*log(ratioP)
!      enddo
!
!      ! compute GZ on top thermo level (from top momentum level)
!      ratioP = col_getPressure(columnghr,1,columnIndex,'MM') / col_getPressure(columnghr,1,columnIndex,'TH')
!      gz_T(1) = gz_M(1) + MPC_RGAS_DRY_AIR_R8*tv(1)*log(ratioP)
!
!      ! compute GZ on remaining thermo levels by simple averaging
!      do lev_T = 2, (nlev_T-1)
!        lev_M = lev_T ! momentum level just below thermo level being computed
!        gz_T(lev_T) = 0.5d0*( gz_M(lev_M-1) + gz_M(lev_M) )
!      enddo
!
!      if ( columnIndex == 1 ) write(*,*) 'MAZIAR: tt2phi, gz_T=', gz_T
!      if ( columnIndex == 1 ) write(*,*) 'MAZIAR: tt2phi, gz_M=', gz_M
!
!    enddo
!
!  elseif(Vcode .eq. 5005) then
!
!    status = vgd_get(columnghr%vco%vgrid,key='DHM - height of the diagnostic level (m)',value=gz_sfcOffset_M_r4)
!    status = vgd_get(columnghr%vco%vgrid,key='DHT - height of the diagnostic level (t)',value=gz_sfcOffset_T_r4)
!    if(mpi_myid == 0 .and. .not.beSilent ) then
!      write(*,*) 'col_fillmvo: height offset for near-sfc momentum level is: ', gz_sfcOffset_M_r4, ' metres'
!      write(*,*) 'col_fillmvo: height offset for near-sfc thermo level is:   ', gz_sfcOffset_T_r4, ' metres'
!    end if
!
!    ! loop over all columns
!    do columnIndex = 1, col_getNumCol(columnghr)
!
!      gz_M => col_getColumn(columnghr,columnIndex,'GZ','MM')
!      gz_T => col_getColumn(columnghr,columnIndex,'GZ','TH')
!
!      ! set the surface height (this is the true surface, not the lowest UU/TT level)
!      gz_M(nlev_M) = col_getGZsfc(columnghr,columnIndex) + real(gz_sfcOffset_M_r4,8) * RG
!      gz_T(nlev_T) = col_getGZsfc(columnghr,columnIndex) + real(gz_sfcOffset_T_r4,8) * RG
!
!      ! initialize the rest to zero
!      gz_M(1:(nlev_M-1)) = 0.0d0
!      gz_T(1:(nlev_T-1)) = 0.0d0
!
!      ! compute virtual temperature on thermo levels
!      do lev_T = 1, nlev_T
!        hu = col_getElem(columnghr,lev_T,columnIndex,'HU')
!        tt = col_getElem(columnghr,lev_T,columnIndex,'TT')
!        tv(lev_T) = fotvt8(tt,hu)
!      enddo
!    
!      ! compute GZ on momentum levels
!      do lev_M = (nlev_M-1), 1, -1
!        lev_T = lev_M ! thermo level just below momentum level being computed
!        if(col_getPressure(columnghr,lev_M,columnIndex,'MM').eq.0.0d0) then
!          write(*,*) 'tt2phi: pressure is zero, lev_m, columnIndex=',lev_m, columnIndex
!          call utl_abort('tt2phi')
!        endif
!        ratioP = col_getPressure(columnghr,lev_M+1,columnIndex,'MM') / col_getPressure(columnghr,lev_M,columnIndex,'MM')
!        gz_M(lev_M) = gz_M(lev_M+1) + MPC_RGAS_DRY_AIR_R8*tv(lev_T)*log(ratioP)
!      enddo
!
!      ! compute GZ on thermo levels by simple averaging
!      do lev_T = 1, (nlev_T-1)
!        lev_M = lev_T+1 ! momentum level just below thermo level being computed
!        gz_T(lev_T) = 0.5d0*( gz_M(lev_M-1) + gz_M(lev_M) )
!      enddo
!
!      if ( columnIndex == 1 ) write(*,*) 'MAZIAR: tt2phi, gz_T=', gz_T
!      if ( columnIndex == 1 ) write(*,*) 'MAZIAR: tt2phi, gz_M=', gz_M
!
!    enddo
!
!  endif
!
!  deallocate(tv)
!
!end subroutine tt2phi !}}}
!

subroutine tt2phi_tl(column,columng) !{{{
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
  real(8), pointer     :: delGz_M(:),delGz_T(:),delTT(:),delHU(:),delP0(:)
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
          hu = col_getElem(columng,lev_T,columnIndex,'HU')
          tt = col_getElem(columng,lev_T,columnIndex,'TT')
          coeff_M_TT(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * fottva(hu,1.0d0)
          coeff_M_HU(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * folnqva(hu,tt,1.0d0) / hu
          coeff_M_P0(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * (delLnP_M(lev_T)-delLnP_M(lev_T-1)) * fotvt8(tt,hu)
        enddo

        ! compute coefficients for GZ on only the top thermo level
        ratioP1 = log( col_getPressure(columng,1,columnIndex,'MM') /  &
                       col_getPressure(columng,1,columnIndex,'TH') )
        delLnP_M1 = col_getPressureDeriv(columng,1,columnIndex,'MM')/  &
                    col_getPressure(columng,1,columnIndex,'MM')
        delLnP_T1 = col_getPressureDeriv(columng,1,columnIndex,'TH')/  &
                    col_getPressure(columng,1,columnIndex,'TH')
        hu = col_getElem(columng,1,columnIndex,'HU')
        tt = col_getElem(columng,1,columnIndex,'TT')
        coeff_T_TT(columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * fottva(hu,1.0d0)
        coeff_T_HU(columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * folnqva(hu,tt,1.0d0) / hu
        coeff_T_P0(columnIndex) = MPC_RGAS_DRY_AIR_R8 * (delLnP_M1-delLnP_T1) * fotvt8(tt,hu)

      enddo

    endif

    ! loop over all columns

!$OMP PARALLEL DO PRIVATE(columnIndex,delGz_M,delGz_T,delThick,delTT,delHU,delP0,lev_M,lev_T)
    do columnIndex = 1, col_getNumCol(columng)

      delGz_M => col_getColumn(column,columnIndex,'GZ','MM')
      delGz_T => col_getColumn(column,columnIndex,'GZ','TH')
      delTT   => col_getColumn(column,columnIndex,'TT')
      delHU   => col_getColumn(column,columnIndex,'HU')
      delP0   => col_getColumn(column,columnIndex,'P0')

      ! ensure increment at sfc is zero (fixed height level)
      delGz_M(nlev_M) = 0.0d0
      delGz_T(nlev_T) = 0.0d0

      ! compute increment to thickness for each layer
      do lev_T = 2, (nlev_T-1)
        delThick(lev_T) = coeff_M_TT(lev_T,columnIndex) * delTT(lev_T) + &
                          coeff_M_HU(lev_T,columnIndex) * delHU(lev_T) + &
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
                   coeff_T_HU(columnIndex) * delHU(1) + &
                   coeff_T_P0(columnIndex) * delP0(1)

      if ( columnIndex == 1 ) then
        write(*,*) 'MAZIAR: tt2phi_tl, delGZ_M/delGZ_T (columnIndex=1):'
        write(*,*) delGz_M(1:nlev_M)
        write(*,*) delGz_T(1:nlev_T)
      endif

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
          hu = col_getElem(columng,lev_T,columnIndex,'HU')
          tt = col_getElem(columng,lev_T,columnIndex,'TT')
          coeff_M_TT(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * fottva(hu,1.0d0)
          coeff_M_HU(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * folnqva(hu,tt,1.0d0) / hu
          coeff_M_P0(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * (delLnP_M(lev_T+1)-delLnP_M(lev_T)) * fotvt8(tt,hu)
        enddo

      enddo

    endif

    ! loop over all columns

!$OMP PARALLEL DO PRIVATE(columnIndex,delGz_M,delGz_T,delThick,delTT,delHU,delP0,lev_M,lev_T)
    do columnIndex = 1, col_getNumCol(columng)

      delGz_M => col_getColumn(column,columnIndex,'GZ','MM')
      delGz_T => col_getColumn(column,columnIndex,'GZ','TH')
      delTT   => col_getColumn(column,columnIndex,'TT')
      delHU   => col_getColumn(column,columnIndex,'HU')
      delP0   => col_getColumn(column,columnIndex,'P0')

      ! ensure increment at sfc is zero (fixed height level)
      delGz_M(nlev_M) = 0.0d0
      delGz_T(nlev_T) = 0.0d0

      ! compute increment to thickness for each layer
      do lev_T = 1, (nlev_T-1)
        delThick(lev_T) = coeff_M_TT(lev_T,columnIndex) * delTT(lev_T) + &
                          coeff_M_HU(lev_T,columnIndex) * delHU(lev_T) + &
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

      if ( columnIndex == 1 ) then
        write(*,*) 'MAZIAR: tt2phi_tl, delGZ_M/delGZ_T (columnIndex=1):'
        write(*,*) delGz_M(1:nlev_M)
        write(*,*) delGz_T(1:nlev_T)
      endif

    enddo
!$OMP END PARALLEL DO

  endif

  deallocate(ratioP)
  deallocate(delThick)
  deallocate(delLnP_M)
  deallocate(delLnP_T)

end subroutine tt2phi_tl !}}}


subroutine tt2phi_ad(column,columng) !{{{
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
  real(8), pointer     :: delGz_M_in(:),delGz_T_in(:),delTT(:),delHU(:),delP0(:)
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
          hu = col_getElem(columng,lev_T,columnIndex,'HU')
          tt = col_getElem(columng,lev_T,columnIndex,'TT')
          coeff_M_TT(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * fottva(hu,1.0d0)
          coeff_M_HU(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * folnqva(hu,tt,1.0d0) / hu
          coeff_M_P0(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * (delLnP_M(lev_T)-delLnP_M(lev_T-1)) * fotvt8(tt,hu)
        enddo

        ! compute coefficients for GZ on only the top thermo level
        ratioP1 = log( col_getPressure(columng,1,columnIndex,'MM') /  &
                       col_getPressure(columng,1,columnIndex,'TH') )
        delLnP_M1 = col_getPressureDeriv(columng,1,columnIndex,'MM')/  &
                    col_getPressure(columng,1,columnIndex,'MM')
        delLnP_T1 = col_getPressureDeriv(columng,1,columnIndex,'TH')/  &
                    col_getPressure(columng,1,columnIndex,'TH')
        hu = col_getElem(columng,1,columnIndex,'HU')
        tt = col_getElem(columng,1,columnIndex,'TT')
        coeff_T_TT(columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * fottva(hu,1.0d0)
        coeff_T_HU(columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * folnqva(hu,tt,1.0d0) / hu
        coeff_T_P0(columnIndex) = MPC_RGAS_DRY_AIR_R8 * (delLnP_M1-delLnP_T1) * fotvt8(tt,hu)

      enddo

    endif

    ! loop over all columns

!$OMP PARALLEL DO PRIVATE(columnIndex,delGz_M_in,delGz_T_in,delTT,  &
!$OMP delHU,delP0,lev_M,lev_T,sumGz_T,delGz_M,delGz_T)
    do columnIndex = 1, col_getNumCol(columng)

      delGz_M_in => col_getColumn(column,columnIndex,'GZ','MM')
      delGz_T_in => col_getColumn(column,columnIndex,'GZ','TH')
      delTT      => col_getColumn(column,columnIndex,'TT')
      delHU      => col_getColumn(column,columnIndex,'HU')
      delP0      => col_getColumn(column,columnIndex,'P0')

      delGz_M(:) = delGz_M_in(:)
      delGz_T(:) = delGz_T_in(:)

      ! adjoint of compute GZ increment on remaining thermo levels by simple averaging
      do lev_T = 2, (nlev_T-1)
        lev_M = lev_T ! momentum level just below thermo level being computed
        delGz_M(lev_M-1) = delGz_M(lev_M-1) + 0.5d0*delGz_T(lev_T)
        delGz_M(lev_M)   = delGz_M(lev_M)   + 0.5d0*delGz_T(lev_T)
      enddo

      ! adjoint of compute GZ increment on top thermo level (from top momentum level)
      delGz_M(1)  = delGz_M(1)  + delGz_T(1)
      delTT(1) = delTT(1) + coeff_T_TT(columnIndex)*delGz_T(1)
      delHU(1) = delHU(1) + coeff_T_HU(columnIndex)*delGz_T(1)
      delP0(1) = delP0(1) + coeff_T_P0(columnIndex)*delGz_T(1)

      ! adjoint of compute GZ increment on momentum levels
      sumGz_T(0:1) = 0.0d0
      do lev_M = 1, (nlev_M-1)
        lev_T = lev_M+1 ! thermo level just below momentum level being computed
        sumGz_T(lev_T) = sumGz_T(lev_T-1) + delGz_M(lev_M)
      enddo
      do lev_T = 2, nlev_T-1
        delTT(lev_T) = delTT(lev_T) + coeff_M_TT(lev_T,columnIndex)*sumGz_T(lev_T)
        delHU(lev_T) = delHU(lev_T) + coeff_M_HU(lev_T,columnIndex)*sumGz_T(lev_T)
        delP0(1)     = delP0(1)     + coeff_M_P0(lev_T,columnIndex)*sumGz_T(lev_T)
      enddo

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
          hu = col_getElem(columng,lev_T,columnIndex,'HU')
          tt = col_getElem(columng,lev_T,columnIndex,'TT')
          coeff_M_TT(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * fottva(hu,1.0d0)
          coeff_M_HU(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * ratioP1 * folnqva(hu,tt,1.0d0) / hu
          coeff_M_P0(lev_T,columnIndex) = MPC_RGAS_DRY_AIR_R8 * (delLnP_M(lev_T+1)-delLnP_M(lev_T)) * fotvt8(tt,hu)
        enddo

      enddo

    endif

    ! loop over all columns

!$OMP PARALLEL DO PRIVATE(columnIndex,delGz_M_in,delGz_T_in,delTT,  &
!$OMP delHU,delP0,lev_M,lev_T,sumGz_T,delGz_M,delGz_T)
    do columnIndex = 1, col_getNumCol(columng)

      delGz_M_in => col_getColumn(column,columnIndex,'GZ','MM')
      delGz_T_in => col_getColumn(column,columnIndex,'GZ','TH')
      delTT      => col_getColumn(column,columnIndex,'TT')
      delHU   => col_getColumn(column,columnIndex,'HU')
      delP0      => col_getColumn(column,columnIndex,'P0')

      delGz_M(:) = delGz_M_in(:)
      delGz_T(:) = delGz_T_in(:)

      ! adjoint of compute GZ increment on remaining thermo levels by simple averaging
      do lev_T = 1, (nlev_T-1)
        lev_M = lev_T+1 ! momentum level just below thermo level being computed
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
        delHU(lev_T) = delHU(lev_T) + coeff_M_HU(lev_T,columnIndex)*sumGz_T(lev_T)
        delP0(1)     = delP0(1)     + coeff_M_P0(lev_T,columnIndex)*sumGz_T(lev_T)
      enddo

    enddo
!$OMP END PARALLEL DO

  endif

  deallocate(ratioP)
  deallocate(delLnP_M)
  deallocate(delLnP_T)
  deallocate(delGz_M)
  deallocate(delGz_T)

end subroutine tt2phi_ad !}}}


! tt2phi used (tt2phi_gpsro)
subroutine tt2phi(columnghr,beSilent_opt) !{{{
  !
  !**s/r tt2phi_gpsro - Temperature to geopotential transformation on GEM4 staggered levels
  !               NOTE: we assume 
  !                     1) nlev_T = nlev_M+1 
  !                     2) GZ_T(nlev_T) = GZ_M(nlev_M), both at the surface
  !                     3) a thermo level exists at the top, higher than the highest momentum level
  !                     4) the placement of the thermo levels means that GZ_T is the average of 2 nearest GZ_M
  !                        (according to Ron and Claude)
  !
  !Author  : M. Buehner, February 2014
  !
  !Revision 001 : M. Bani Shahabadi, October 2018
  !          - adaptation of GPSRO calculation of GZ
  !          - Gravity acceleration include 2nd-order Eotvos effect for non-linear operator.
  !
  implicit none

  type(struct_columnData) :: columnghr
  logical, optional       :: beSilent_opt

  integer :: columnIndex,lev_M,lev_T,nlev_M,nlev_T,status,Vcode
  real(8) :: hu, tt, Pr, cmp, delThick, tvm, ratioP, ScaleFactorBottom, ScaleFactorTop
  real(8), allocatable :: tv(:), AL_T(:), AL_M(:) 
  real(8), pointer     :: gz_T(:),gz_M(:)
  real                 :: gz_sfcOffset_T_r4, gz_sfcOffset_M_r4
  real(8) :: rLat, rLon, latrot, lonrot, xpos, ypos, rMT, rUU, rVV
  real(8) :: h0, dh, Rgh, Eot, Eot2, sLat, cLat
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
  if (Vcode == 5002 .and. nlev_T /= nlev_M+1) call utl_abort('tt2phi_gpsro: nlev_T is not equal to nlev_M+1!')
  if (Vcode == 5005 .and. nlev_T /= nlev_M)   call utl_abort('tt2phi_gpsro: nlev_T is not equal to nlev_M!')

  allocate(tv(nlev_T))
  allocate(AL_T(nlev_T))
  allocate(AL_M(nlev_M))

  write(*,*) 'MAZIAR: tt2phi_gpsro, Vcode=', Vcode

  ! loop over all columns
  do columnIndex = 1, col_getNumCol(columnghr)

    gz_M => col_getColumn(columnghr,columnIndex,'GZ','MM')
    gz_T => col_getColumn(columnghr,columnIndex,'GZ','TH')

    ! initialize the GZ/AL to zero
    gz_M(:) = 0.0d0
    gz_T(:) = 0.0d0
    AL_T(1:nlev_T) = 0.0D0
    AL_M(1:nlev_M) = 0.0D0

    ! latitude/longitude
    call col_getLatLon(columnghr, columnIndex,                  & ! IN
                        rLat, rLon, ypos, xpos, LatRot, LonRot )  ! OUT
    sLat = sin(rLat)
    cLat = cos(rLat)

    ! compute virtual temperature on thermo levels (corrected of compressibility)
    do lev_T = 1, nlev_T
      hu = col_getElem(columnghr,lev_T,columnIndex,'HU')
      tt = col_getElem(columnghr,lev_T,columnIndex,'TT')
      Pr = col_getPressure(columnghr,lev_T,columnIndex,'TH')

      cmp = gpscompressibility(Pr,tt,hu)

      tv(lev_T) = fotvt8(tt,hu) * cmp 
    enddo

    ! compute altitude on bottom thermo level (Vcode=5005)
    if (Vcode == 5002) then
      Rgh = gpsgravitysrf(sLat)
      rMT = col_getGZsfc(columnghr,columnIndex) / Rgh
      AL_T(nlev_T) = rMT
    elseif (Vcode == 5005) then
      ratioP  = log(col_getPressure(columnghr,nlev_T,columnIndex,'TH') / &
                col_getElem(columnghr,1,columnIndex,'P0') )
      Rgh = gpsgravitysrf(sLat)
      delThick = (-p_Rd / Rgh) * tv(nlev_T) * ratioP
      rMT = col_getGZsfc(columnghr,columnIndex) / Rgh
      AL_T(nlev_T) = rMT + delThick
    endif

    ! compute altitude on rest of thermo levels
    do lev_T = nlev_T-1, 1, -1
      ratioP = log(col_getPressure(columnghr,lev_T  ,columnIndex,'TH') / &
                   col_getPressure(columnghr,lev_T+1,columnIndex,'TH'))
      tvm = 0.5D0 * (tv(lev_T) + tv(lev_T+1))
      
      ! UU/VV
      if (Vcode == 5002) then
        lev_M = lev_T
      elseif (Vcode == 5005) then
        lev_M = lev_T + 1
      endif
      rUU = col_getElem(columnghr,lev_M,columnIndex,'UU') * p_knot
      rVV = col_getElem(columnghr,lev_M,columnIndex,'VV') * p_knot

      ! Gravity acceleration 
      h0  = AL_T(lev_T+1)
      Eot = 2 * WGS_OmegaPrime * cLat * p_knot * rUU
      Eot2= ((p_knot*rUU) ** 2 + (p_knot*rVV) ** 2) / WGS_a
      Rgh = gpsgravityalt(sLat, h0) - Eot - Eot2
      dh  = (-p_Rd / Rgh) * tvm * ratioP
      Rgh = gpsgravityalt(sLat, h0+0.5D0*dh) - Eot - Eot2

      delThick   = (-p_Rd / Rgh) * tvm * ratioP
      AL_T(lev_T) = AL_T(lev_T+1) + delThick
    enddo

    ! compute Altitude on momentum levels
    if (Vcode == 5002) then
      AL_M(nlev_M) = AL_T(nlev_T)

      do lev_M = 1, nlev_M-1
        lev_T = lev_M + 1 ! thermo level just below momentum level being computed

        ScaleFactorBottom = log( col_getPressure(columnghr,lev_M  ,columnIndex,'MM')   / &
                                 col_getPressure(columnghr,lev_T-1,columnIndex,'TH') ) / &
                            log( col_getPressure(columnghr,lev_T  ,columnIndex,'TH')   / &
                                 col_getPressure(columnghr,lev_T-1,columnIndex,'TH') )
        ScaleFactorTop    = 1 - ScaleFactorBottom

        AL_M(lev_M) = ScaleFactorBottom * AL_T(lev_T) + ScaleFactorTop * AL_T(lev_T-1)
      end do

    elseif (Vcode == 5005) then

      do lev_M = 2, nlev_M
        lev_T = lev_M ! thermo level just below momentum level being computed

        ScaleFactorBottom = log( col_getPressure(columnghr,lev_M  ,columnIndex,'MM')   / &
                                 col_getPressure(columnghr,lev_T-1,columnIndex,'TH') ) / &
                            log( col_getPressure(columnghr,lev_T  ,columnIndex,'TH')   / &
                                 col_getPressure(columnghr,lev_T-1,columnIndex,'TH') )
        ScaleFactorTop    = 1 - ScaleFactorBottom

        AL_M(lev_M) = ScaleFactorBottom * AL_T(lev_T) + ScaleFactorTop * AL_T(lev_T-1)
      end do

      ! compute altitude on top momentum level
      ratioP = log(col_getPressure(columnghr,1,columnIndex,'MM') / &
                   col_getPressure(columnghr,1,columnIndex,'TH'))
      tvm = tv(1)
      rUU = col_getElem(columnghr,1,columnIndex,'UU') * p_knot
      rVV = col_getElem(columnghr,1,columnIndex,'VV') * p_knot

      ! Gravity acceleration 
      h0  = AL_T(1)
      Eot = 2 * WGS_OmegaPrime * cLat * p_knot * rUU
      Eot2= ((p_knot*rUU) ** 2 + (p_knot*rVV) ** 2) / WGS_a
      Rgh = gpsgravityalt(sLat, h0) - Eot - Eot2
      dh  = (-p_Rd / Rgh) * tvm * ratioP
      Rgh = gpsgravityalt(sLat, h0+0.5D0*dh) - Eot - Eot2

      delThick   = (-p_Rd / Rgh) * tvm * ratioP
      AL_M(1) = AL_T(1) + delThick

    endif

    ! compute GZ 
    gz_T(1:nlev_T) = AL_T(1:nlev_T) * RG
    gz_M(1:nlev_M) = AL_M(1:nlev_M) * RG

    if ( columnIndex == 1 ) then

      write(*,*) 'MAZIAR: tt2phi_gpsro, gz_T=', gz_T
      write(*,*) 'MAZIAR: tt2phi_gpsro, gz_M=', gz_M

      write(*,*) 'MAZIAR: tt2phi_gpsro, AL_M/AL_T='
      write(*,*) AL_M
      write(*,*) AL_T
    endif

  enddo

  deallocate(tv)
  deallocate(AL_T)
  deallocate(AL_M)

end subroutine tt2phi !}}}


subroutine tt2phi_gpsro_1Col(columnghr,columnIndex,beSilent_opt) !{{{
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
  !Revision 001 : M. Bani Shahabadi, October 2018
  !          - adaptation of GPSRO calculation of GZ
  !          - Gravity acceleration include 2nd-order Eotvos effect for non-linear operator.
  !
  implicit none

  type(struct_columnData) :: columnghr
  logical, optional       :: beSilent_opt

  integer :: columnIndex,lev_M,lev_T,nlev_M,nlev_T,status,Vcode
  real(8) :: hu, tt, Pr, cmp, delThick, tvm, ratioP, ScaleFactorBottom, ScaleFactorTop
  real(8), allocatable :: tv(:), AL_T(:), AL_M(:) 
  real(8), pointer     :: gz_T(:),gz_M(:)
  real                 :: gz_sfcOffset_T_r4, gz_sfcOffset_M_r4
  real(8) :: rLat, rLon, latrot, lonrot, xpos, ypos, rMT, rUU, rVV
  real(8) :: h0, dh, Rgh, Eot, Eot2, sLat, cLat
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
  if (Vcode == 5002 .and. nlev_T /= nlev_M+1) call utl_abort('tt2phi_gpsro_1Col: nlev_T is not equal to nlev_M+1!')
  if (Vcode == 5005 .and. nlev_T /= nlev_M)   call utl_abort('tt2phi_gpsro_1Col: nlev_T is not equal to nlev_M!')

  allocate(tv(nlev_T))
  allocate(AL_T(nlev_T))
  allocate(AL_M(nlev_M))
  AL_T(1:nlev_T) = 0.0D0
  AL_M(1:nlev_M) = 0.0D0

  write(*,*) 'MAZIAR: tt2phi_gpsro_1Col, Vcode=', Vcode

  ! latitude/longitude
  call col_getLatLon(columnghr, columnIndex,                  & ! IN
                      rLat, rLon, ypos, xpos, LatRot, LonRot )  ! OUT
  sLat = sin(rLat)
  cLat = cos(rLat)

  ! compute virtual temperature on thermo levels (corrected of compressibility)
  do lev_T = 1, nlev_T
    hu = col_getElem(columnghr,lev_T,columnIndex,'HU')
    tt = col_getElem(columnghr,lev_T,columnIndex,'TT')
    Pr = col_getPressure(columnghr,lev_T,columnIndex,'TH')

    cmp = gpscompressibility(Pr,tt,hu)

    tv(lev_T) = fotvt8(tt,hu) * cmp 
  enddo

  ! compute altitude on bottom thermo level (Vcode=5005)
  if (Vcode == 5005) then
    ! first thermo level for Vcode=5005
    ratioP  = log(col_getPressure(columnghr,nlev_T,columnIndex,'TH') / &
              col_getElem(columnghr,1,columnIndex,'P0') )
    Rgh = gpsgravitysrf(sLat)
    delThick = (-p_Rd / Rgh) * tv(nlev_T) * ratioP
    rMT = col_getGZsfc(columnghr,columnIndex) / Rgh
    AL_T(nlev_T) = rMT + delThick
  endif

  ! compute altitude on rest of thermo levels
  do lev_T = nlev_T-1, 1, -1
    ratioP = log(col_getPressure(columnghr,lev_T  ,columnIndex,'TH') / &
                 col_getPressure(columnghr,lev_T+1,columnIndex,'TH'))
    tvm = 0.5D0 * (tv(lev_T) + tv(lev_T+1))
    
    ! UU/VV
    if (Vcode == 5002) then
      lev_M = lev_T
    elseif (Vcode == 5005) then
      lev_M = lev_T + 1
    endif
    rUU = col_getElem(columnghr,lev_M,columnIndex,'UU') * p_knot
    rVV = col_getElem(columnghr,lev_M,columnIndex,'VV') * p_knot

    ! Gravity acceleration 
    h0  = AL_T(lev_T+1)
    Eot = 2 * WGS_OmegaPrime * cLat * rUU
    Eot2= (rUU ** 2 + rVV ** 2) / WGS_a
    Rgh = gpsgravityalt(sLat, h0) - Eot - Eot2
    dh  = (-p_Rd / Rgh) * tvm * ratioP
    Rgh = gpsgravityalt(sLat, h0+0.5D0*dh) - Eot - Eot2

    delThick   = (-p_Rd / Rgh) * tvm * ratioP
    AL_T(lev_T) = AL_T(lev_T+1) + delThick
  enddo

  ! compute Altitude on momentum levels
  if (Vcode == 5002) then

    do lev_M = 1, nlev_M
      lev_T = lev_M + 1 ! thermo level just below momentum level being computed

      ScaleFactorBottom = log( col_getPressure(columnghr,lev_M  ,columnIndex,'MM')   / &
                               col_getPressure(columnghr,lev_T-1,columnIndex,'TH') ) / &
                          log( col_getPressure(columnghr,lev_T  ,columnIndex,'TH')   / &
                               col_getPressure(columnghr,lev_T-1,columnIndex,'TH') )
      ScaleFactorTop    = 1 - ScaleFactorBottom

      AL_M(lev_M) = ScaleFactorBottom * AL_T(lev_T) + ScaleFactorTop * AL_T(lev_T-1)
    end do

  elseif (Vcode == 5005) then

    do lev_M = 2, nlev_M
      lev_T = lev_M ! thermo level just below momentum level being computed

      ScaleFactorBottom = log( col_getPressure(columnghr,lev_M  ,columnIndex,'MM')   / &
                               col_getPressure(columnghr,lev_T-1,columnIndex,'TH') ) / &
                          log( col_getPressure(columnghr,lev_T  ,columnIndex,'TH')   / &
                               col_getPressure(columnghr,lev_T-1,columnIndex,'TH') )
      ScaleFactorTop    = 1 - ScaleFactorBottom

      AL_M(lev_M) = ScaleFactorBottom * AL_T(lev_T) + ScaleFactorTop * AL_T(lev_T-1)
    end do

    ! compute altitude on top momentum level
    ratioP = log(col_getPressure(columnghr,1,columnIndex,'MM') / &
                 col_getPressure(columnghr,1,columnIndex,'TH'))
    tvm = tv(1)
    rUU = col_getElem(columnghr,1,columnIndex,'UU') * p_knot
    rVV = col_getElem(columnghr,1,columnIndex,'VV') * p_knot

    ! Gravity acceleration 
    h0  = AL_T(1)
    Eot = 2 * WGS_OmegaPrime * cLat * rUU
    Eot2= (rUU ** 2 + rVV ** 2) / WGS_a
    Rgh = gpsgravityalt(sLat, h0) - Eot - Eot2
    dh  = (-p_Rd / Rgh) * tvm * ratioP
    Rgh = gpsgravityalt(sLat, h0+0.5D0*dh) - Eot - Eot2

    delThick   = (-p_Rd / Rgh) * tvm * ratioP
    AL_M(1) = AL_T(1) + delThick

  endif

  write(*,*) 'MAZIAR: tt2phi_gpsro_1Col, AL_M/AL_T='
  write(*,*) AL_M
  write(*,*) AL_T

  deallocate(tv)
  deallocate(AL_T)
  deallocate(AL_M)

end subroutine tt2phi_gpsro_1Col !}}}


subroutine tt2phi_gpsro_tl(column,columng) !{{{
  !
  !**s/r tt2phi_gpsro_tl - Temperature to geopotential transformation on GEM4 staggered levels
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
  !Revision 001 : M. Bani Shahabadi, October 2018
  !          - adaptation of GPSRO calculation of delGZ
  !
  implicit none

  type(struct_columnData) :: column, columng

  integer :: columnIndex,lev_M,lev_T,nlev_M,nlev_T,status,Vcode_anl
  integer :: lev_gps, levCompare
  real(8) :: hu,tt,Pr,cmp,cmp_TT,cmp_HU,cmp_P0,delLnP_M1,delLnP_T1, ScaleFactorBottom, ScaleFactorTop
  real(8), allocatable :: delThick(:)
  real(8), allocatable :: delLnP_T(:)
  real(8)              :: delLnPsfc 
  real(8), pointer     :: gz_T(:),gz_M(:)
  real(8), pointer     :: delGz_M(:),delGz_T(:),delTT(:),delHU(:),delP0(:)
  real(8) :: rLat, rLon, latrot, lonrot, xpos, ypos
  real(8) :: h0, Rgh, sLat, cLat
  type(struct_vco), pointer :: vco_anl

  real(8), allocatable, save :: coeff_T_TT(:,:), coeff_T_HU(:,:), coeff_T_P0(:,:), &
                                ratioP(:,:), ratioP1(:), &
                                diff_delLnP_T(:,:), diff_delLnP_M1(:), & 
                                coeff_T_P0_dP(:,:)

  logical, save :: firstTime = .true.

  vco_anl => col_getVco(columng)
  status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)

  nlev_T = col_getNumLev(columng,'TH')
  nlev_M = col_getNumLev(columng,'MM')

  allocate(delThick(nlev_T-1))
  allocate(delLnP_T(nlev_T))

  if(firstTime) then
    ! initialize and save coefficients for increased efficiency (assumes no relinearization)
    firstTime = .false.

    allocate(coeff_T_TT   (nlev_T,col_getNumCol(columng)))
    allocate(coeff_T_HU   (nlev_T,col_getNumCol(columng)))
    allocate(coeff_T_P0   (nlev_T,col_getNumCol(columng)))
    allocate(coeff_T_P0_dP(nlev_T,col_getNumCol(columng)))

    allocate(ratioP       (nlev_T-1,col_getNumCol(columng)))
    allocate(diff_delLnP_T(nlev_T-1,col_getNumCol(columng)))

    allocate(ratioP1(col_getNumCol(columng)))
    allocate(diff_delLnP_M1(col_getNumCol(columng)))

    coeff_T_TT(:,:) = 0.0D0
    coeff_T_HU(:,:) = 0.0D0
    coeff_T_P0(:,:) = 0.0D0
    coeff_T_P0_dP(:,:) = 0.0D0
    ratioP(:,:) = 0.0D0
    ratioP1(:) = 0.0D0
    diff_delLnP_T(:,:) = 0.0D0
    diff_delLnP_M1(:) = 0.0D0

    write(*,*) 'MAZIAR: tt2phi_gpsro_tl, Vcode_anl=',Vcode_anl,',coef done.'

    do columnIndex = 1, col_getNumCol(columng)

      gz_T => col_getColumn(columng,columnIndex,'GZ','TH')

      ! latitude/longitude
      call col_getLatLon(columng, columnIndex,                  & ! IN
                          rLat, rLon, ypos, xpos, LatRot, LonRot )! OUT
      sLat = sin(rLat)
      cLat = cos(rLat)

      ! compute thermo level properties
      do lev_T = 1,nlev_T
        delLnP_T(lev_T) = col_getPressureDeriv(columng,lev_T,columnIndex,'TH')/  &
                         col_getPressure(columng,lev_T,columnIndex,'TH')

        hu = col_getElem(columng,lev_T,columnIndex,'HU')
        tt = col_getElem(columng,lev_T,columnIndex,'TT')
        Pr = col_getPressure(columng,lev_T,columnIndex,'TH')

        cmp = gpscompressibility(Pr,tt,hu)
        cmp_TT = gpscompressibility_TT(Pr,tt,hu)
        cmp_HU = gpscompressibility_HU(Pr,tt,hu)
        cmp_P0 = gpscompressibility_P0(Pr,tt,hu,col_getPressureDeriv(columng,lev_T,columnIndex,'TH'))

        ! Gravity acceleration 
        h0  = gz_T(lev_T) / RG
        Rgh = gpsgravityalt(sLat, h0)

        coeff_T_TT(lev_T,columnIndex) = (p_Rd / Rgh) * (fottva(hu,1.0D0) * cmp + fotvt8(tt,hu) * cmp_TT)
        coeff_T_HU(lev_T,columnIndex) = (p_Rd / Rgh) * (folnqva(hu,tt,1.0d0) / hu * cmp + fotvt8(tt,hu) * cmp_HU)
        coeff_T_P0(lev_T,columnIndex) = (p_Rd / Rgh) * fotvt8(tt,hu) * cmp
        coeff_T_P0_dP(lev_T,columnIndex) = (p_Rd / Rgh) * fotvt8(tt,hu) * cmp_P0
      enddo

      ! compute layer properties (between two adjacent thermo levels)
      do lev_T = 1,nlev_T-1
        ratioP(lev_T,columnIndex) = log( col_getPressure(columng,lev_T+1,columnIndex,'TH') /  &
                                         col_getPressure(columng,lev_T  ,columnIndex,'TH') )
        diff_delLnP_T(lev_T,columnIndex) = delLnP_T(lev_T+1) - delLnP_T(lev_T)
      enddo

      ! compute the property of the top layer (between first momentum and first thermo level for Vcode=5005)
      if (Vcode_anl == 5005) then

        ratioP1(columnIndex)      = log( col_getPressure(columng,1,columnIndex,'TH') /  &
                                         col_getPressure(columng,1,columnIndex,'MM') )

        delLnP_M1 = col_getPressureDeriv(columng,1,columnIndex,'MM') / &
                         col_getPressure(columng,1,columnIndex,'MM')
        delLnP_T1 = col_getPressureDeriv(columng,1,columnIndex,'TH') / &
                         col_getPressure(columng,1,columnIndex,'TH')
        diff_delLnP_M1(columnIndex) = delLnP_T1 - delLnP_M1

      endif

    enddo

  endif

!$OMP PARALLEL DO PRIVATE(columnIndex,gz_M,gz_T,delGz_M,delGz_T,delThick,delTT,delHU,delP0,lev_M,lev_T,ScaleFactorBottom,ScaleFactorTop)
  do columnIndex = 1, col_getNumCol(columng)

    gz_M => col_getColumn(columng,columnIndex,'GZ','MM')
    gz_T => col_getColumn(columng,columnIndex,'GZ','TH')

    delGz_M => col_getColumn(column,columnIndex,'GZ','MM')
    delGz_T => col_getColumn(column,columnIndex,'GZ','TH')
    delTT   => col_getColumn(column,columnIndex,'TT')
    delHU   => col_getColumn(column,columnIndex,'HU')
    delP0   => col_getColumn(column,columnIndex,'P0')

    ! ensure increment at sfc/fixed height level is zero
    delGz_M(nlev_M) = 0.0d0
    delGz_T(nlev_T) = 0.0d0

    ! compute increment to thickness for each layer between the two thermo levels
    do lev_T = nlev_T-1, 1, -1
      delThick(lev_T) = 0.5D0 * coeff_T_TT(lev_T  ,columnIndex) * delTT(lev_T  ) * ratioP(lev_T,columnIndex) + &
                        0.5D0 * coeff_T_TT(lev_T+1,columnIndex) * delTT(lev_T+1) * ratioP(lev_T,columnIndex) + &
                        0.5D0 * coeff_T_HU(lev_T  ,columnIndex) * delHU(lev_T  ) * ratioP(lev_T,columnIndex) + &
                        0.5D0 * coeff_T_HU(lev_T+1,columnIndex) * delHU(lev_T+1) * ratioP(lev_T,columnIndex) + &
                        0.5D0 * (coeff_T_P0   (lev_T,columnIndex) + coeff_T_P0   (lev_T+1,columnIndex)) * &
                                delP0(1) * diff_delLnP_T(lev_T,columnIndex) + & 
                        0.5D0 * (coeff_T_P0_dP(lev_T,columnIndex) + coeff_T_P0_dP(lev_T+1,columnIndex)) * &
                                delP0(1) * ratioP(lev_T,columnIndex)

    enddo

    ! compute GZ increment on thermo levels above the surface
    do lev_T = nlev_T-1, 1, -1
      delGz_T(lev_T) = delGz_T(lev_T+1) + delThick(lev_T) * RG
    enddo

    ! compute GZ increment on momentum levels using weighted average of GZ increment of thermo levels
    if (Vcode_anl == 5002) then

      do lev_M = 1, nlev_M-1

        lev_T = lev_M + 1 ! thermo level just below momentum level being computed

        ScaleFactorBottom = (gz_M(lev_M) - gz_T(lev_T-1)) / &
                            (gz_T(lev_T) - gz_T(lev_T-1))
        ScaleFactorTop    = 1 - ScaleFactorBottom

        delGz_M(lev_M) = ScaleFactorBottom * delGz_T(lev_T) + ScaleFactorTop * delGz_T(lev_T-1)
      enddo

    elseif (Vcode_anl == 5005) then

      do lev_M = 2, nlev_M-1

        lev_T = lev_M     ! thermo level just below momentum level being computed

        ScaleFactorBottom = (gz_M(lev_M) - gz_T(lev_T-1)) / &
                            (gz_T(lev_T) - gz_T(lev_T-1))
        ScaleFactorTop    = 1 - ScaleFactorBottom

        delGz_M(lev_M) = ScaleFactorBottom * delGz_T(lev_T) + ScaleFactorTop * delGz_T(lev_T-1)
      enddo

      ! compute GZ increment for top momentum level (from top thermo level)
      delGz_M(1) = delGz_T(1) +  RG * ( &
                    0.5D0 * coeff_T_TT(1,columnIndex) * delTT(1) * ratioP1(columnIndex) + &
                    0.5D0 * coeff_T_HU(1,columnIndex) * delHU(1) * ratioP1(columnIndex) + &
                    0.5D0 * coeff_T_P0(1,columnIndex) * delP0(1) * diff_delLnP_M1(columnIndex) + & 
                    0.5D0 * coeff_T_P0_dP(1,columnIndex) * delP0(1) * ratioP1(columnIndex) )

    endif

    ! print del altitude
    if ( columnIndex == 1 ) then
      write(*,*) 'MAZIAR: tt2phi_gpsro_tl, delGZ_M/delGZ_T (columnIndex=1):'
      write(*,*) delGz_M(1:nlev_M)
      write(*,*) delGz_T(1:nlev_T)
    endif

  enddo
!$OMP END PARALLEL DO

  deallocate(delThick)
  deallocate(delLnP_T)

end subroutine tt2phi_gpsro_tl !}}}


subroutine tt2phi_gpsro_tl_1Col(column,columng,columnIndex) !{{{
  !
  !**s/r tt2phi_gpsro_tl_1Col - Temperature to geopotential transformation on GEM4 staggered levels
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
  !Revision 001 : M. Bani Shahabadi, October 2018
  !          - adaptation of GPSRO calculation of delGZ
  !
  implicit none

  type(struct_columnData) :: column, columng

  integer :: columnIndex,lev_M,lev_T,nlev_M,nlev_T,status,Vcode_anl
  integer :: lev_gps, levCompare
  real(8) :: hu,tt,Pr,cmp,cmp_TT,cmp_HU,cmp_P0,delLnP_M1,delLnP_T1, ScaleFactorBottom, ScaleFactorTop
  real(8), allocatable :: delThick(:)
  real(8), allocatable :: delLnP_T(:)
  real(8)              :: delLnPsfc 
  real(8), pointer     :: gz_T(:),gz_M(:)
  real(8), allocatable :: delGz_M(:),delGz_T(:)
  real(8), pointer     :: delTT(:),delHU(:),delP0(:)
  real(8) :: rLat, rLon, latrot, lonrot, xpos, ypos
  real(8) :: h0, Rgh, sLat, cLat
  type(struct_vco), pointer :: vco_anl

  real(8), allocatable, save :: coeff_T_TT(:,:), coeff_T_HU(:,:), coeff_T_P0(:,:), &
                                ratioP(:,:), ratioP1(:), &
                                diff_delLnP_T(:,:), diff_delLnP_M1(:), & 
                                coeff_T_P0_dP(:,:)

  logical, save :: firstTime = .true.

  vco_anl => col_getVco(columng)
  status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)

  nlev_T = col_getNumLev(columng,'TH')
  nlev_M = col_getNumLev(columng,'MM')

  allocate(delThick(nlev_T-1))
  allocate(delLnP_T(nlev_T))

  allocate(delGz_M(nlev_M))
  allocate(delGz_T(nlev_T))
  delGz_M(1:nlev_M) = 0.0D0
  delGz_T(1:nlev_T) = 0.0D0 

  if(firstTime) then
    ! initialize and save coefficients for increased efficiency (assumes no relinearization)
    firstTime = .false.

    allocate(coeff_T_TT   (nlev_T,col_getNumCol(columng)))
    allocate(coeff_T_HU   (nlev_T,col_getNumCol(columng)))
    allocate(coeff_T_P0   (nlev_T,col_getNumCol(columng)))
    allocate(coeff_T_P0_dP(nlev_T,col_getNumCol(columng)))

    allocate(ratioP       (nlev_T-1,col_getNumCol(columng)))
    allocate(diff_delLnP_T(nlev_T-1,col_getNumCol(columng)))

    allocate(ratioP1(col_getNumCol(columng)))
    allocate(diff_delLnP_M1(col_getNumCol(columng)))

    coeff_T_TT(:,:) = 0.0D0
    coeff_T_HU(:,:) = 0.0D0
    coeff_T_P0(:,:) = 0.0D0
    coeff_T_P0_dP(:,:) = 0.0D0
    ratioP(:,:) = 0.0D0
    ratioP1(:) = 0.0D0
    diff_delLnP_T(:,:) = 0.0D0
    diff_delLnP_M1(:) = 0.0D0

!    do columnIndex = 1, col_getNumCol(columng)

      gz_T => col_getColumn(columng,columnIndex,'GZ','TH')

      ! latitude/longitude
      call col_getLatLon(columng, columnIndex,                  & ! IN
                          rLat, rLon, ypos, xpos, LatRot, LonRot )! OUT
      sLat = sin(rLat)
      cLat = cos(rLat)

      ! compute thermo level properties
      do lev_T = 1,nlev_T
        delLnP_T(lev_T) = col_getPressureDeriv(columng,lev_T,columnIndex,'TH')/  &
                         col_getPressure(columng,lev_T,columnIndex,'TH')

        hu = col_getElem(columng,lev_T,columnIndex,'HU')
        tt = col_getElem(columng,lev_T,columnIndex,'TT')
        Pr = col_getPressure(columng,lev_T,columnIndex,'TH')

        cmp = gpscompressibility(Pr,tt,hu)
        cmp_TT = gpscompressibility_TT(Pr,tt,hu)
        cmp_HU = gpscompressibility_HU(Pr,tt,hu)
        cmp_P0 = gpscompressibility_P0(Pr,tt,hu,col_getPressureDeriv(columng,lev_T,columnIndex,'TH'))

        ! Gravity acceleration 
        h0  = gz_T(lev_T) / RG
        Rgh = gpsgravityalt(sLat, h0)

        coeff_T_TT(lev_T,columnIndex) = (p_Rd / Rgh) * (fottva(hu,1.0D0) * cmp + fotvt8(tt,hu) * cmp_TT)
        coeff_T_HU(lev_T,columnIndex) = (p_Rd / Rgh) * (folnqva(hu,tt,1.0d0) / hu * cmp + fotvt8(tt,hu) * cmp_HU)
        coeff_T_P0(lev_T,columnIndex) = (p_Rd / Rgh) * fotvt8(tt,hu) * cmp
        coeff_T_P0_dP(lev_T,columnIndex) = (p_Rd / Rgh) * fotvt8(tt,hu) * cmp_P0
      enddo

      ! compute layer properties (between two adjacent thermo levels)
      do lev_T = 1,nlev_T-1
        ratioP(lev_T,columnIndex) = log( col_getPressure(columng,lev_T+1,columnIndex,'TH') /  &
                                         col_getPressure(columng,lev_T  ,columnIndex,'TH') )
        diff_delLnP_T(lev_T,columnIndex) = delLnP_T(lev_T+1) - delLnP_T(lev_T)
      enddo

      ! compute the property of the top layer (between first momentum and first thermo level for Vcode=5005)
      if (Vcode_anl == 5005) then

        ratioP1(columnIndex)      = log( col_getPressure(columng,1,columnIndex,'TH') /  &
                                         col_getPressure(columng,1,columnIndex,'MM') )

        delLnP_M1 = col_getPressureDeriv(columng,1,columnIndex,'MM') / &
                         col_getPressure(columng,1,columnIndex,'MM')
        delLnP_T1 = col_getPressureDeriv(columng,1,columnIndex,'TH') / &
                         col_getPressure(columng,1,columnIndex,'TH')
        diff_delLnP_M1(columnIndex) = delLnP_T1 - delLnP_M1

      endif

!    enddo

  endif

  write(*,*) 'MAZIAR: tt2phi_gpsro_tl_1Col, Vcode_anl=',Vcode_anl,',coef done.'

!!$OMP PARALLEL DO PRIVATE(columnIndex,delGz_M,delGz_T,delThick,delTT,delHU,delP0,lev_M,lev_T)
!  do columnIndex = 1, col_getNumCol(columng)

    gz_M => col_getColumn(columng,columnIndex,'GZ','MM')
    gz_T => col_getColumn(columng,columnIndex,'GZ','TH')

    !delGz_M => col_getColumn(column,columnIndex,'GZ','MM')
    !delGz_T => col_getColumn(column,columnIndex,'GZ','TH')
    delTT   => col_getColumn(column,columnIndex,'TT')
    delHU   => col_getColumn(column,columnIndex,'HU')
    delP0   => col_getColumn(column,columnIndex,'P0')

    ! ensure increment at sfc is zero (fixed height level)
    delGz_M(nlev_M) = 0.0d0
    delGz_T(nlev_T) = 0.0d0

    ! compute increment to thickness for each layer between the two thermo levels
    do lev_T = nlev_T-1, 1, -1
      delThick(lev_T) = 0.5D0 * coeff_T_TT(lev_T  ,columnIndex) * delTT(lev_T  ) * ratioP(lev_T,columnIndex) + &
                        0.5D0 * coeff_T_TT(lev_T+1,columnIndex) * delTT(lev_T+1) * ratioP(lev_T,columnIndex) + &
                        0.5D0 * coeff_T_HU(lev_T  ,columnIndex) * delHU(lev_T  ) * ratioP(lev_T,columnIndex) + &
                        0.5D0 * coeff_T_HU(lev_T+1,columnIndex) * delHU(lev_T+1) * ratioP(lev_T,columnIndex) + &
                        0.5D0 * (coeff_T_P0   (lev_T,columnIndex) + coeff_T_P0   (lev_T+1,columnIndex)) * &
                                delP0(1) * diff_delLnP_T(lev_T,columnIndex) + & 
                        0.5D0 * (coeff_T_P0_dP(lev_T,columnIndex) + coeff_T_P0_dP(lev_T+1,columnIndex)) * &
                                delP0(1) * ratioP(lev_T,columnIndex)

    enddo

    ! compute GZ increment on thermo levels above the surface
    do lev_T = nlev_T-1, 1, -1
      delGz_T(lev_T) = delGz_T(lev_T+1) + delThick(lev_T) * RG
    enddo

    ! compute GZ increment on momentum levels using weighted average of GZ increment of thermo levels
    if (Vcode_anl == 5002) then

      do lev_M = 1, nlev_M-1

        lev_T = lev_M + 1 ! thermo level just below momentum level being computed

        ScaleFactorBottom = (gz_M(lev_M) - gz_T(lev_T-1)) / &
                            (gz_T(lev_T) - gz_T(lev_T-1))
        ScaleFactorTop    = 1 - ScaleFactorBottom

        delGz_M(lev_M) = ScaleFactorBottom * delGz_T(lev_T) + ScaleFactorTop * delGz_T(lev_T-1)
      enddo

    elseif (Vcode_anl == 5005) then

      do lev_M = 2, nlev_M-1

        lev_T = lev_M     ! thermo level just below momentum level being computed

        ScaleFactorBottom = (gz_M(lev_M) - gz_T(lev_T-1)) / &
                            (gz_T(lev_T) - gz_T(lev_T-1))
        ScaleFactorTop    = 1 - ScaleFactorBottom

        delGz_M(lev_M) = ScaleFactorBottom * delGz_T(lev_T) + ScaleFactorTop * delGz_T(lev_T-1)
      enddo

      ! compute GZ increment for top momentum level (from top thermo level)
      delGz_M(1) = delGz_T(1) +  &
                    0.5D0 * coeff_T_TT(1,columnIndex) * delTT(1) * ratioP1(columnIndex) + &
                    0.5D0 * coeff_T_HU(1,columnIndex) * delHU(1) * ratioP1(columnIndex) + &
                    0.5D0 * coeff_T_P0(1,columnIndex) * delP0(1) * diff_delLnP_M1(columnIndex) + & 
                    0.5D0 * coeff_T_P0_dP(1,columnIndex) * delP0(1) * ratioP1(columnIndex)

    endif

!  enddo
!!$OMP END PARALLEL DO

  ! print del altitude
  write(*,*) 'MAZIAR: tt2phi_gpsro_tl_1Col, delAL_M/delAL_T:'
  write(*,*) delGz_M(1:nlev_M) / RG
  write(*,*) delGz_T(1:nlev_T) / RG

  deallocate(delThick)
  deallocate(delLnP_T)

end subroutine tt2phi_gpsro_tl_1Col !}}}


subroutine tt2phi_gpsro_ad(column,columng) !{{{
  !
  !**s/r tt2phi_gpsro_ad - Adjoint of temperature to geopotential transformation on GEM4 staggered levels
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
  !Revision 001 : M. Bani Shahabadi, October 2018
  !          - adaptation of GPSRO calculation of delGZ
  !
  implicit none

  type(struct_columnData) :: column,columng

  integer :: columnIndex,lev_M,lev_T,nlev_M,nlev_T,status,Vcode_anl
  integer :: lev_gps, levCompare
  real(8) :: hu,tt,Pr,cmp,cmp_TT,cmp_HU,cmp_P0,delLnP_M1,delLnP_T1, ScaleFactorBottom, ScaleFactorTop
  real(8), allocatable :: delThick(:)
  real(8), allocatable :: delLnP_T(:)
  real(8)              :: delLnPsfc 
  real(8), pointer     :: gz_T(:),gz_M(:)
  real(8), allocatable :: delGz_M(:),delGz_T(:)
  real(8), pointer     :: delGz_M_in(:),delGz_T_in(:),delTT(:),delHU(:),delP0(:)
  real(8) :: rLat, rLon, latrot, lonrot, xpos, ypos
  real(8) :: h0, Rgh, sLat, cLat
  type(struct_vco), pointer :: vco_anl

  real(8), allocatable, save :: coeff_T_TT(:,:), coeff_T_HU(:,:), coeff_T_P0(:,:), &
                                ratioP(:,:), ratioP1(:), &
                                diff_delLnP_T(:,:), diff_delLnP_M1(:), & 
                                coeff_T_P0_dP(:,:)

  logical, save :: firstTime = .true.

  vco_anl => col_getVco(columng)
  status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)

  nlev_T = col_getNumLev(columng,'TH')
  nlev_M = col_getNumLev(columng,'MM')

  allocate(delThick(0:nlev_T-1))
  allocate(delLnP_T(nlev_T))

  allocate(delGz_M(nlev_M))
  allocate(delGz_T(nlev_T))

  if(firstTime) then
    ! initialize and save coefficients for increased efficiency (assumes no relinearization)
    firstTime = .false.

    allocate(coeff_T_TT   (nlev_T,col_getNumCol(columng)))
    allocate(coeff_T_HU   (nlev_T,col_getNumCol(columng)))
    allocate(coeff_T_P0   (nlev_T,col_getNumCol(columng)))
    allocate(coeff_T_P0_dP(nlev_T,col_getNumCol(columng)))

    allocate(ratioP       (nlev_T-1,col_getNumCol(columng)))
    allocate(diff_delLnP_T(nlev_T-1,col_getNumCol(columng)))

    allocate(ratioP1(col_getNumCol(columng)))
    allocate(diff_delLnP_M1(col_getNumCol(columng)))

    coeff_T_TT(:,:) = 0.0D0
    coeff_T_HU(:,:) = 0.0D0
    coeff_T_P0(:,:) = 0.0D0
    coeff_T_P0_dP(:,:) = 0.0D0
    ratioP(:,:) = 0.0D0
    ratioP1(:) = 0.0D0
    diff_delLnP_T(:,:) = 0.0D0
    diff_delLnP_M1(:) = 0.0D0

    do columnIndex = 1, col_getNumCol(columng)

      gz_T => col_getColumn(columng,columnIndex,'GZ','TH')

      ! latitude/longitude
      call col_getLatLon(columng, columnIndex,                  & ! IN
                          rLat, rLon, ypos, xpos, LatRot, LonRot )! OUT
      sLat = sin(rLat)
      cLat = cos(rLat)

      ! compute thermo level properties
      do lev_T = 1,nlev_T
        delLnP_T(lev_T) = col_getPressureDeriv(columng,lev_T,columnIndex,'TH')/  &
                         col_getPressure(columng,lev_T,columnIndex,'TH')

        hu = col_getElem(columng,lev_T,columnIndex,'HU')
        tt = col_getElem(columng,lev_T,columnIndex,'TT')
        Pr = col_getPressure(columng,lev_T,columnIndex,'TH')

        cmp = gpscompressibility(Pr,tt,hu)
        cmp_TT = gpscompressibility_TT(Pr,tt,hu)
        cmp_HU = gpscompressibility_HU(Pr,tt,hu)
        cmp_P0 = gpscompressibility_P0(Pr,tt,hu,col_getPressureDeriv(columng,lev_T,columnIndex,'TH'))

        ! Gravity acceleration 
        h0  = gz_T(lev_T) / RG
        Rgh = gpsgravityalt(sLat, h0)

        coeff_T_TT(lev_T,columnIndex) = (p_Rd / Rgh) * (fottva(hu,1.0D0) * cmp + fotvt8(tt,hu) * cmp_TT)
        coeff_T_HU(lev_T,columnIndex) = (p_Rd / Rgh) * (folnqva(hu,tt,1.0d0) / hu * cmp + fotvt8(tt,hu) * cmp_HU)
        coeff_T_P0(lev_T,columnIndex) = (p_Rd / Rgh) * fotvt8(tt,hu) * cmp
        coeff_T_P0_dP(lev_T,columnIndex) = (p_Rd / Rgh) * fotvt8(tt,hu) * cmp_P0
      enddo

      ! compute layer properties (between two adjacent thermo levels)
      do lev_T = 1,nlev_T-1
        ratioP(lev_T,columnIndex) = log( col_getPressure(columng,lev_T+1,columnIndex,'TH') /  &
                                         col_getPressure(columng,lev_T  ,columnIndex,'TH') )
        diff_delLnP_T(lev_T,columnIndex) = delLnP_T(lev_T+1) - delLnP_T(lev_T)
      enddo

      ! compute the property of the top layer (between first momentum and first thermo level for Vcode=5005)
      if (Vcode_anl == 5005) then

        ratioP1(columnIndex)      = log( col_getPressure(columng,1,columnIndex,'TH') /  &
                                         col_getPressure(columng,1,columnIndex,'MM') )

        delLnP_M1 = col_getPressureDeriv(columng,1,columnIndex,'MM') / &
                         col_getPressure(columng,1,columnIndex,'MM')
        delLnP_T1 = col_getPressureDeriv(columng,1,columnIndex,'TH') / &
                         col_getPressure(columng,1,columnIndex,'TH')
        diff_delLnP_M1(columnIndex) = delLnP_T1 - delLnP_M1

      endif

    enddo

  endif

!$OMP PARALLEL DO PRIVATE(columnIndex,delGz_M,delGz_T,delThick,delTT,delHU,delP0,lev_M,lev_T,delGz_M_in,delGz_T_in,gz_M,gz_T,ScaleFactorBottom,ScaleFactorTop)
  do columnIndex = 1, col_getNumCol(columng)

    gz_M => col_getColumn(columng,columnIndex,'GZ','MM')
    gz_T => col_getColumn(columng,columnIndex,'GZ','TH')

    delGz_M_in => col_getColumn(column,columnIndex,'GZ','MM')
    delGz_T_in => col_getColumn(column,columnIndex,'GZ','TH')
    delTT      => col_getColumn(column,columnIndex,'TT')
    delHU      => col_getColumn(column,columnIndex,'HU')
    delP0      => col_getColumn(column,columnIndex,'P0')

    delGz_M(:) = delGz_M_in(:)
    delGz_T(:) = delGz_T_in(:)

    ! adjoint of compute GZ increment on momentum levels using weighted average of GZ increment of thermo levels
    if (Vcode_anl == 5002) then

      do lev_M = 1, nlev_M-1

        lev_T = lev_M + 1 ! thermo level just below momentum level being computed

        ScaleFactorBottom = (gz_M(lev_M) - gz_T(lev_T-1)) / &
                            (gz_T(lev_T) - gz_T(lev_T-1))
        ScaleFactorTop    = 1 - ScaleFactorBottom

        delGz_T(lev_T-1) = delGz_T(lev_T-1) + ScaleFactorTop     * delGz_M(lev_M)
        delGz_T(lev_T  ) = delGz_T(lev_T  ) + ScaleFactorBottom  * delGz_M(lev_M)

      enddo

    elseif (Vcode_anl == 5005) then

      do lev_M = 2, nlev_M-1

        lev_T = lev_M     ! thermo level just below momentum level being computed

        ScaleFactorBottom = (gz_M(lev_M) - gz_T(lev_T-1)) / &
                            (gz_T(lev_T) - gz_T(lev_T-1))
        ScaleFactorTop    = 1 - ScaleFactorBottom

        delGz_T(lev_T-1) = delGz_T(lev_T-1) + ScaleFactorTop     * delGz_M(lev_M)
        delGz_T(lev_T  ) = delGz_T(lev_T  ) + ScaleFactorBottom  * delGz_M(lev_M)

      enddo

      ! adjoint of compute GZ increment for top momentum level (from top thermo level)
      delGz_T(1) = delGz_T(1) + delGz_M(1)
      delTT(1) = delTT(1) + 0.5D0 * coeff_T_TT(1,columnIndex) * delGz_M(1) * ratioP1(columnIndex) * RG
      delHU(1) = delHU(1) + 0.5D0 * coeff_T_HU(1,columnIndex) * delGz_M(1) * ratioP1(columnIndex) * RG
      delP0(1) = delP0(1) + 0.5D0 * coeff_T_P0(1,columnIndex) * delGz_M(1) * diff_delLnP_M1(columnIndex) * RG + & 
                            0.5D0 * coeff_T_P0_dP(1,columnIndex) * delGz_M(1) * ratioP1(columnIndex) * RG

    endif

    ! adjoint of compute GZ increment on thermo levels above the surface
    delThick(0) = 0.0D0
    do lev_T = 1, nlev_T-1
      delThick(lev_T) = delThick(lev_T-1) + delGz_T(lev_T) * RG
    enddo

    ! adjoint of compute increment to thickness for each layer between the two thermo levels
    do lev_T = nlev_T-1, 1, -1

      delTT(lev_T  ) = delTT(lev_T  ) + 0.5D0 * coeff_T_TT(lev_T  ,columnIndex) * delThick(lev_T) * ratioP(lev_T,columnIndex)
      delTT(lev_T+1) = delTT(lev_T+1) + 0.5D0 * coeff_T_TT(lev_T+1,columnIndex) * delThick(lev_T) * ratioP(lev_T,columnIndex)
      delHU(lev_T  ) = delHU(lev_T  ) + 0.5D0 * coeff_T_HU(lev_T  ,columnIndex) * delThick(lev_T) * ratioP(lev_T,columnIndex)
      delHU(lev_T+1) = delHU(lev_T+1) + 0.5D0 * coeff_T_HU(lev_T+1,columnIndex) * delThick(lev_T) * ratioP(lev_T,columnIndex)
      delP0(1) = delP0(1) +  0.5D0 * (coeff_T_P0   (lev_T,columnIndex) + coeff_T_P0   (lev_T+1,columnIndex)) * delThick(lev_T) * diff_delLnP_T(lev_T,columnIndex) + &
                             0.5D0 * (coeff_T_P0_dP(lev_T,columnIndex) + coeff_T_P0_dP(lev_T+1,columnIndex)) * delThick(lev_T) * ratioP(lev_T,columnIndex)

    enddo

    !! print delTT/delHU/delP0
    !if ( columnIndex == 1 ) then
    !  write(*,*) 'MAZIAR: tt2phi_gpsro_ad, delTT/delHU/delP0 (columnIndex=1):'
    !  write(*,*) delTT(1:nlev_T)
    !  write(*,*) delHU(1:nlev_T)
    !  write(*,*) delP0(1:nlev_T)
    !endif

  enddo
!$OMP END PARALLEL DO

  deallocate(delThick)
  deallocate(delLnP_T)
  deallocate(delGz_M)
  deallocate(delGz_T)

end subroutine tt2phi_gpsro_ad !}}}


function gpscompressibility(p,t,q)
  real(8), intent(in)  :: p,t,q
  real(8)              :: gpscompressibility

  real(8), parameter   :: a0= 1.58123D-6
  real(8), parameter   :: a1=-2.9331D-8
  real(8), parameter   :: a2= 1.1043D-10
  real(8), parameter   :: b0= 5.707D-6
  real(8), parameter   :: b1=-2.051D-8
  real(8), parameter   :: c0= 1.9898D-4
  real(8), parameter   :: c1=-2.376D-6
  real(8), parameter   :: d = 1.83D-11
  real(8), parameter   :: e =-0.765D-8

  real(8)         :: x,tc,pt,tc2,x2

  x  = p_wa * q / (1.D0 + p_wb * q)
  ! Estimate, from CIPM, Picard (2008)
  tc = t - p_TC
  pt = p / t
  tc2= tc * tc
  x2 = x * x
  gpscompressibility = 1.D0 - pt * (a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2) + pt*pt*(d+e*x2)
end function gpscompressibility

function gpscompressibility_TT(p,t,q)
  real(8), intent(in)  :: p,t,q
  real(8)              :: gpscompressibility_TT

  real(8), parameter   :: a0= 1.58123D-6
  real(8), parameter   :: a1=-2.9331D-8
  real(8), parameter   :: a2= 1.1043D-10
  real(8), parameter   :: b0= 5.707D-6
  real(8), parameter   :: b1=-2.051D-8
  real(8), parameter   :: c0= 1.9898D-4
  real(8), parameter   :: c1=-2.376D-6
  real(8), parameter   :: d = 1.83D-11
  real(8), parameter   :: e =-0.765D-8

  real(8)         :: x,tc,pt,tc2,x2
  real(8)         :: d_x,d_tc,d_pt,d_tc2,d_x2

  x  = p_wa * q / (1.D0 + p_wb * q)
  d_x  = 0.0D0
  ! Estimate, from CIPM, Picard (2008)
  tc = t - p_TC
  d_tc = 1.D0
  pt = p / t
  d_pt = - p / t**2
  tc2= tc * tc
  d_tc2= 2 * tc * d_tc
  x2 = x * x
  d_x2 = 2 * x * d_x
  gpscompressibility_TT = -d_pt * (a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2) + &
          -pt * (a1*d_tc + a2*d_tc2 + b1*d_tc*x + (b0+b1*tc)*d_x + c1*d_tc*x2 + (c0+c1*tc)*d_x2) + &
          2*pt*d_pt*(d+e*x2) + pt*pt*e*d_x2
end function gpscompressibility_TT

function gpscompressibility_HU(p,t,q)
  real(8), intent(in)  :: p,t,q
  real(8)              :: gpscompressibility_HU

  real(8), parameter   :: a0= 1.58123D-6
  real(8), parameter   :: a1=-2.9331D-8
  real(8), parameter   :: a2= 1.1043D-10
  real(8), parameter   :: b0= 5.707D-6
  real(8), parameter   :: b1=-2.051D-8
  real(8), parameter   :: c0= 1.9898D-4
  real(8), parameter   :: c1=-2.376D-6
  real(8), parameter   :: d = 1.83D-11
  real(8), parameter   :: e =-0.765D-8

  real(8)         :: x,tc,pt,tc2,x2
  real(8)         :: d_x,d_tc,d_pt,d_tc2,d_x2

  x  = p_wa * q / (1.D0+p_wb*q)
  d_x  = p_wa * (1.0D0 / (1.D0+p_wb*q) - q / (1.D0+p_wb*q)**2 * p_wb * 1.0D0)
  ! Estimate, from CIPM, Picard (2008)
  tc = t - p_TC
  d_tc = 0.0D0
  pt = p / t
  d_pt = 0.0D0
  tc2= tc * tc
  d_tc2= 2 * tc * d_tc
  x2 = x * x
  d_x2 = 2 * x * d_x
  gpscompressibility_HU = -d_pt * (a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2) + &
          -pt * (a1*d_tc + a2*d_tc2 + b1*d_tc*x + (b0+b1*tc)*d_x + c1*d_tc*x2 + (c0+c1*tc)*d_x2) + &
          2*pt*d_pt*(d+e*x2) + pt*pt*e*d_x2
end function gpscompressibility_HU

function gpscompressibility_P0(p,t,q,dpdp0)
  real(8), intent(in)  :: p,t,q,dpdp0
  real(8)              :: gpscompressibility_P0

  real(8), parameter   :: a0= 1.58123D-6
  real(8), parameter   :: a1=-2.9331D-8
  real(8), parameter   :: a2= 1.1043D-10
  real(8), parameter   :: b0= 5.707D-6
  real(8), parameter   :: b1=-2.051D-8
  real(8), parameter   :: c0= 1.9898D-4
  real(8), parameter   :: c1=-2.376D-6
  real(8), parameter   :: d = 1.83D-11
  real(8), parameter   :: e =-0.765D-8

  real(8)         :: x,tc,pt,tc2,x2
  real(8)         :: d_x,d_tc,d_pt,d_tc2,d_x2

  x  = p_wa * q / (1.D0+p_wb*q)
  d_x  = 0.0D0
  ! Estimate, from CIPM, Picard (2008)
  tc = t - p_TC
  d_tc = 0.0D0
  pt = p / t
  d_pt = dpdp0 / t
  tc2= tc * tc
  d_tc2= 2 * tc * d_tc
  x2 = x * x
  d_x2 = 2 * x * d_x
  gpscompressibility_P0 = -d_pt * (a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2) + &
          -pt * (a1*d_tc + a2*d_tc2 + b1*d_tc*x + (b0+b1*tc)*d_x + c1*d_tc*x2 + (c0+c1*tc)*d_x2) + &
          2*pt*d_pt*(d+e*x2) + pt*pt*e*d_x2
end function gpscompressibility_P0

function gpsgravitysrf(sLat)
  !  Normal gravity on ellipsoidal surface:
  !  Input:  Latitude
  !          sin(Latitude)
  !
  !  Output: Normal gravity
  !          gpsgravitysrf         : m/s2
  !
  real(8), intent(in)  :: sLat
  real(8)              :: gpsgravitysrf
  
  real(8)              :: ks2
  real(8)              :: e2s

  ks2 = WGS_TNGk * sLat*sLat
  e2s = 1.D0 - WGS_e2 * sLat*sLat
  gpsgravitysrf = WGS_GammaE * (1.D0 + ks2) / sqrt(e2s)
end function gpsgravitysrf

function gpsgravityalt(sLat, Altitude)
  ! Normal gravity above the ellipsoidal surface:
  ! Input:  Latitude, altitude
  !         sin(Latitude)
  !         Altitude               : m
  !
  ! Output: Normal gravity
  !         gpsgravityalt          : m/s2
  !
  real(8), intent(in)  :: sLat
  real(8), intent(in)  :: Altitude
  real(8)              :: gpsgravityalt

  real(8)              :: C1
  real(8)              :: C2

  C1 =-2.D0/WGS_a*(1.D0+WGS_f+WGS_m-2*WGS_f*sLat*sLat)
  C2 = 3.D0/WGS_a**2
  gpsgravityalt = gpsgravitysrf(sLat)*                                   &
       (1.D0 + C1 * Altitude + C2 * Altitude**2)
end function gpsgravityalt

function gpsgeopotential(Latitude, Altitude)
  ! Geopotential energy at a given point.
  ! Result is based on the WGS84 approximate expression for the
  ! gravity acceleration as a function of latitude and altitude,
  ! integrated with the trapezoidal rule.
  ! Input:  Latitude, altitude
  !         Latitude               : rad
  !         Altitude               : m
  !
  ! Output: Geopotential
  !         gpsgeopotential                              : m2/s2
  !
  real(8), intent(in)  :: Latitude
  real(8), intent(in)  :: Altitude
  real(8)              :: gpsgeopotential

  real(8)              :: dh, sLat
  integer               :: n, i
  real(8), allocatable :: hi(:)
  real(8), allocatable :: gi(:)
  
  dh = 500.D0
  n = 1 + int(Altitude/dh)

  allocate(hi(0:n))
  allocate(gi(0:n))

  sLat=sin(Latitude)

  do i = 0, n-1
     hi(i) = i * dh
     gi(i) = gpsgravityalt(sLat, hi(i))
  enddo
  hi(n) = Altitude
  gi(n) = gpsgravityalt(sLat, hi(n))

  gpsgeopotential = 0.D0
  do i = 1, n
     gpsgeopotential = gpsgeopotential + 0.5D0 * (gi(i)+gi(i-1)) * (hi(i)-hi(i-1))
  enddo

  deallocate(hi)
  deallocate(gi)
end function gpsgeopotential

end module tt2phi_mod
