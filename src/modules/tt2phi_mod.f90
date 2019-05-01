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
!! *Purpose*: Subroutines for computing height from TT, HU and P0. Nonlinear, tangent-
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
  use obsSpaceData_mod
  implicit none
  save
  private

  ! public procedures
  public :: tt2phi           , tt2phi_tl           , tt2phi_ad

  ! constants from gps_mod
  ! Air properties:
  real(8), parameter :: p_md = 28.965516D0            ! From Aparicio(2011)
  real(8), parameter :: p_mw = 18.015254D0            ! From Aparicio(2011)
  real(8), parameter :: p_wa = p_md/p_mw
  real(8), parameter :: p_wb = (p_md-p_mw)/p_mw
  ! Angular velocity of the Earth (omegaPrime) (radians/s).
  ! Standard Earth, rotating with a constant angular velocity (IAU, GRS67).
  real(8), parameter :: WGS_OmegaPrime = 7292115.1467D-11

  ! private module variables
  real(8), allocatable, save :: coeff_M_TT(:,:), coeff_M_HU(:,:)
  real(8), allocatable, save :: coeff_T_TT(:), coeff_T_HU(:)
  real(8), allocatable, save :: coeff_M_P0(:,:), coeff_M_P0_dP(:,:)
  real(8), allocatable, save :: coeff_T_P0(:), coeff_T_P0_dP(:)

contains

subroutine tt2phi(columnghr,obsSpaceData,beSilent_opt)
  !
  !**s/r tt2phi - Temperature to geopotential transformation on GEM4 staggered levels
  !               NOTE: we assume 
  !                     1) nlev_T = nlev_M+1 
  !                     2) alt_T(nlev_T) = alt_M(nlev_M), both at the surface
  !                     3) a thermo level exists at the top, higher than the highest momentum level
  !                     4) the placement of the thermo levels means that alt_T is the average of 2 nearest alt_M
  !                        (according to Ron and Claude)
  !
  !Author  : M. Buehner, February 2014
  !
  !Revision 001 : M. Bani Shahabadi, October 2018
  !          - adaptation of GPSRO calculation of height
  !
  implicit none

  type(struct_columnData) :: columnghr
  type(struct_obs)        :: obsSpaceData
  logical, optional       :: beSilent_opt

  integer :: columnIndex,lev_M,lev_T,nlev_M,nlev_T,status,Vcode
  real(8) :: hu, tt, Pr, cmp, delThick, tvm, ratioP, ScaleFactorBottom, ScaleFactorTop
  real(8), allocatable :: tv(:), alt_T(:), alt_M(:) 
  real(8), pointer     :: alt_T_ptr(:),alt_M_ptr(:)
  real                 :: alt_sfcOffset_T_r4, alt_sfcOffset_M_r4
  real(8) :: rLat, rMT
  real(8) :: h0, dh, Rgh, sLat, cLat
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
  if (Vcode == 5002 .and. nlev_T /= nlev_M+1) call utl_abort('tt2phi: nlev_T is not equal to nlev_M+1!')
  if (Vcode == 5005 .and. nlev_T /= nlev_M)   call utl_abort('tt2phi: nlev_T is not equal to nlev_M!')

  if (Vcode == 5005) then
    status = vgd_get(columnghr%vco%vgrid,key='DHM - height of the diagnostic level (m)',value=alt_sfcOffset_M_r4)
    status = vgd_get(columnghr%vco%vgrid,key='DHT - height of the diagnostic level (t)',value=alt_sfcOffset_T_r4)
    if ( mpi_myid == 0 .and. .not.beSilent ) then
      write(*,*) 'tt2phi: height offset for near-sfc momentum level is: ', alt_sfcOffset_M_r4, ' metres'
      write(*,*) 'tt2phi: height offset for near-sfc thermo level is:   ', alt_sfcOffset_T_r4, ' metres'
      if ( .not.columnghr%addGZsfcOffset ) then
        write(*,*) '----------------------------------------------------------------------------------'
        write(*,*) 'tt2phi: BUT HEIGHT OFFSET REMOVED FOR DIAGNOSTIC LEVELS FOR BACKWARD COMPATIBILITY'
        write(*,*) '----------------------------------------------------------------------------------'
      end if
    end if
  end if

  allocate(tv(nlev_T))
  allocate(alt_T(nlev_T))
  allocate(alt_M(nlev_M))

  ! loop over all columns
  do columnIndex = 1, col_getNumCol(columnghr)

    alt_M_ptr => col_getColumn(columnghr,columnIndex,'GZ','MM')
    alt_T_ptr => col_getColumn(columnghr,columnIndex,'GZ','TH')

    ! initialize the ALT_ptr/ALT to zero
    alt_M_ptr(:) = 0.0d0
    alt_T_ptr(:) = 0.0d0
    alt_T(1:nlev_T) = 0.0D0
    alt_M(1:nlev_M) = 0.0D0

    ! latitude
    rLat = obs_headElem_r(obsSpaceData,OBS_LAT,columnIndex)
    sLat = sin(rLat)
    cLat = cos(rLat)

    ! compute virtual temperature on thermo levels (corrected of compressibility)
    do lev_T = 1, nlev_T
      hu = col_getElem(columnghr,lev_T,columnIndex,'HU')
      tt = col_getElem(columnghr,lev_T,columnIndex,'TT')
      Pr = col_getPressure(columnghr,lev_T,columnIndex,'TH')

      cmp = gpscompressibility(Pr,tt,hu)

      tv(lev_T) = fotvt8(tt,hu) * cmp 
    end do

    rMT = col_getHeight(columnghr,0,columnIndex,'SF')

    ! compute altitude on bottom momentum level
    if (Vcode == 5002) then
      alt_M(nlev_M) = rMT
    else if (Vcode == 5005) then
      alt_M(nlev_M) = rMT + alt_sfcOffset_M_r4
    end if

    ratioP  = log(col_getPressure(columnghr,nlev_M-1,columnIndex,'MM') / &
              col_getElem(columnghr,1,columnIndex,'P0') )

    ! Gravity acceleration 
    h0  = rMT
    Rgh = phf_gravityalt(sLat,h0)
    dh  = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(nlev_T-1) * ratioP
    Rgh = phf_gravityalt(sLat, h0+0.5D0*dh)

    delThick = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(nlev_T-1) * ratioP
    alt_M(nlev_M-1) = rMT + delThick

    ! compute altitude on rest of momentum levels
    do lev_M = nlev_M-2, 1, -1
      ratioP = log(col_getPressure(columnghr,lev_M  ,columnIndex,'MM') / &
                   col_getPressure(columnghr,lev_M+1,columnIndex,'MM'))
      
      if (Vcode == 5002) then
        lev_T = lev_M + 1
      elseif (Vcode == 5005) then
        lev_T = lev_M
      end if

      ! Gravity acceleration 
      h0  = alt_M(lev_M+1)
      Rgh = phf_gravityalt(sLat,h0)
      dh  = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(lev_T) * ratioP
      Rgh = phf_gravityalt(sLat, h0+0.5D0*dh)

      delThick   = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(lev_T) * ratioP
      alt_M(lev_M) = alt_M(lev_M+1) + delThick
    end do

    ! compute Altitude on thermo levels
    if (Vcode == 5002) then
      alt_T(nlev_T) = alt_M(nlev_M)

      do lev_T = 2, nlev_T-1
        lev_M = lev_T ! momentum level just below thermo level being computed
        ScaleFactorBottom = log( col_getPressure(columnghr,lev_T  ,columnIndex,'TH')   / &
                                 col_getPressure(columnghr,lev_M-1,columnIndex,'MM') ) / &
                            log( col_getPressure(columnghr,lev_M  ,columnIndex,'MM')   / &
                                 col_getPressure(columnghr,lev_M-1,columnIndex,'MM') )
        ScaleFactorTop    = 1 - ScaleFactorBottom
        alt_T(lev_T) = ScaleFactorBottom * alt_M(lev_M) + ScaleFactorTop * alt_M(lev_M-1)
      end do

      ! compute altitude on top thermo level
      ratioP = log(col_getPressure(columnghr,1,columnIndex,'TH') / &
                   col_getPressure(columnghr,1,columnIndex,'MM'))

      ! Gravity acceleration 
      h0  = alt_M(1)
      Rgh = phf_gravityalt(sLat, h0)
      dh  = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(1) * ratioP
      Rgh = phf_gravityalt(sLat, h0+0.5D0*dh)

      delThick   = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(1) * ratioP
      alt_T(1) = alt_M(1) + delThick

    elseif (Vcode == 5005) then
      alt_T(nlev_T) = rMT + alt_sfcOffset_T_r4
      do lev_T = 1, nlev_T-2
        lev_M = lev_T + 1  ! momentum level just below thermo level being computed
        ScaleFactorBottom = log( col_getPressure(columnghr,lev_T  ,columnIndex,'TH')   / &
                                 col_getPressure(columnghr,lev_M-1,columnIndex,'MM') ) / &
                            log( col_getPressure(columnghr,lev_M  ,columnIndex,'MM')   / &
                                 col_getPressure(columnghr,lev_M-1,columnIndex,'MM') )
        ScaleFactorTop    = 1 - ScaleFactorBottom
        alt_T(lev_T) = ScaleFactorBottom * alt_M(lev_M) + ScaleFactorTop * alt_M(lev_M-1)
      end do

      ! compute altitude on next to bottom thermo level
      ratioP  = log(col_getPressure(columnghr,nlev_T-1,columnIndex,'TH') / &
                col_getElem(columnghr,1,columnIndex,'P0') )

      h0  = rMT
      Rgh = phf_gravityalt(sLat,h0)
      dh  = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(nlev_T-1) * ratioP
      Rgh = phf_gravityalt(sLat, h0+0.5D0*dh)

      delThick = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(nlev_T-1) * ratioP
      alt_T(nlev_T-1) = rMT + delThick
    end if

    ! fill the height array
    alt_T_ptr(1:nlev_T) = alt_T(1:nlev_T)
    alt_M_ptr(1:nlev_M) = alt_M(1:nlev_M)

    ! remove the height offset for the diagnostic levels for backward compatibility only
    if ( .not.columnghr%addGZsfcOffset ) then
      alt_T_ptr(nlev_T) = rMT
      alt_M_ptr(nlev_M) = rMT
    end if

  end do

  deallocate(tv)
  deallocate(alt_T)
  deallocate(alt_M)

end subroutine tt2phi


subroutine tt2phi_tl(column,columng,obsSpaceData)
  !
  !**s/r tt2phi_tl - Temperature to geopotential transformation on GEM4 staggered levels
  !               NOTE: we assume 
  !                     1) nlev_T = nlev_M+1 
  !                     2) alt_T(nlev_T) = alt_M(nlev_M), both at the surface
  !                     3) a thermo level exists at the top, higher than the highest momentum level
  !                     4) the placement of the thermo levels means that alt_T is the average of 2 nearest alt_M
  !                        (according to Ron and Claude)
  !
  ! NOTE: after revision 680 removed code for vcode=5005 when rewriting for increased efficiency (M. Buehner)
  !
  !Author  : M. Buehner, February 2014
  !
  !Revision 001 : M. Bani Shahabadi, October 2018
  !          - adaptation of GPSRO calculation of delAL
  !
  implicit none

  type(struct_columnData) :: column, columng
  type(struct_obs)        :: obsSpaceData

  integer :: columnIndex,lev_M,lev_T,nlev_M,nlev_T,status,Vcode_anl
  real(8) :: ScaleFactorBottom, ScaleFactorTop
  real(8), allocatable :: delThick(:)
  real(8), pointer     :: alt_T_ptr(:),alt_M_ptr(:)
  real(8), pointer     :: delAlt_M_ptr(:),delAlt_T_ptr(:),delTT(:),delHU(:),delP0(:)
  type(struct_vco), pointer :: vco_anl

  vco_anl => col_getVco(columng)
  status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)

  nlev_T = col_getNumLev(columng,'TH')
  nlev_M = col_getNumLev(columng,'MM')

  allocate(delThick(nlev_T))

  ! generate the height coefficients
  call calcAltitudeCoeff(columng,obsSpaceData)

!$OMP PARALLEL DO PRIVATE(columnIndex,alt_M_ptr,alt_T_ptr,delAlt_M_ptr,delAlt_T_ptr,delThick,delTT,delHU,delP0,lev_M,lev_T,ScaleFactorBottom,ScaleFactorTop)
  do columnIndex = 1, col_getNumCol(columng)

    alt_M_ptr => col_getColumn(columng,columnIndex,'GZ','MM')
    alt_T_ptr => col_getColumn(columng,columnIndex,'GZ','TH')

    delAlt_M_ptr => col_getColumn(column,columnIndex,'GZ','MM')
    delAlt_T_ptr => col_getColumn(column,columnIndex,'GZ','TH')
    delTT   => col_getColumn(column,columnIndex,'TT')
    delHU   => col_getColumn(column,columnIndex,'HU')
    delP0   => col_getColumn(column,columnIndex,'P0')

    ! ensure increment at sfc/fixed height level is zero
    delAlt_M_ptr(nlev_M) = 0.0d0
    delAlt_T_ptr(nlev_T) = 0.0d0

    if (Vcode_anl == 5002) then

      ! compute increment to thickness for each layer between the two momentum levels
      do lev_T = 2, (nlev_T-1)
        delThick(lev_T) = coeff_M_TT   (lev_T,columnIndex) * delTT(lev_T) + &
                          coeff_M_HU   (lev_T,columnIndex) * delHU(lev_T) + &
                          coeff_M_P0   (lev_T,columnIndex) * delP0(1)     + & 
                          coeff_M_P0_dP(lev_T,columnIndex) * delP0(1) 

      end do

      ! compute height increment on momentum levels above the surface
      do lev_M = nlev_M-1, 1, -1
        lev_T = lev_M + 1 ! thermo level just below momentum level being computed
        delAlt_M_ptr(lev_M) = delAlt_M_ptr(lev_M+1) + delThick(lev_T)
      end do

      ! compute height increment on thermo levels using weighted average of height increment of momentum levels
      do lev_T = 2, (nlev_T-1)
        lev_M = lev_T ! momentum level just below thermo level being computed
        ScaleFactorBottom = (alt_T_ptr(lev_T) - alt_M_ptr(lev_M-1)) / &
                            (alt_M_ptr(lev_M) - alt_M_ptr(lev_M-1))
        ScaleFactorTop    = 1 - ScaleFactorBottom
        delAlt_T_ptr(lev_T) = ScaleFactorBottom * delAlt_M_ptr(lev_M) + ScaleFactorTop * delAlt_M_ptr(lev_M-1)
      end do

      ! compute height increment for top thermo level (from top momentum level)
      delAlt_T_ptr(1) = delAlt_M_ptr(1) + &
                    coeff_T_TT(columnIndex) * delTT(1) + &
                    coeff_T_HU(columnIndex) * delHU(1) + &
                    coeff_T_P0(columnIndex) * delP0(1) + & 
                    coeff_T_P0_dP(columnIndex) * delP0(1)

    elseif (Vcode_anl == 5005) then

      ! compute increment to thickness for each layer between the two momentum levels
      do lev_T = 1, (nlev_T-1)
        delThick(lev_T) = coeff_M_TT   (lev_T,columnIndex) * delTT(lev_T) + &
                          coeff_M_HU   (lev_T,columnIndex) * delHU(lev_T) + &
                          coeff_M_P0   (lev_T,columnIndex) * delP0(1)     + & 
                          coeff_M_P0_dP(lev_T,columnIndex) * delP0(1) 

      end do

      ! compute height increment on momentum levels above the surface
      do lev_M = nlev_M-1, 1, -1
        lev_T = lev_M ! thermo level just below momentum level being computed
        delAlt_M_ptr(lev_M) = delAlt_M_ptr(lev_M+1) + delThick(lev_T)
      end do

      ! compute height increment on thermo levels using weighted average of height increment of momentum levels
      do lev_T = 1, (nlev_T-1)
        lev_M = lev_T + 1 ! momentum level just below thermo level being computed
        ScaleFactorBottom = (alt_T_ptr(lev_T) - alt_M_ptr(lev_M-1)) / &
                            (alt_M_ptr(lev_M) - alt_M_ptr(lev_M-1))
        ScaleFactorTop    = 1 - ScaleFactorBottom
        delAlt_T_ptr(lev_T) = ScaleFactorBottom * delAlt_M_ptr(lev_M) + ScaleFactorTop * delAlt_M_ptr(lev_M-1)
      end do

    end if

  end do
!$OMP END PARALLEL DO

  deallocate(delThick)

end subroutine tt2phi_tl


subroutine tt2phi_ad(column,columng,obsSpaceData)
  !
  !**s/r tt2phi_ad - Adjoint of temperature to geopotential transformation on GEM4 staggered levels
  !               NOTE: we assume 
  !                     1) nlev_T = nlev_M+1 
  !                     2) alt_T(nlev_T) = alt_M(nlev_M), both at the surface
  !                     3) a thermo level exists at the top, higher than the highest momentum level
  !                     4) the placement of the thermo levels means that alt_T is the average of 2 nearest alt_M
  !                        (according to Ron and Claude)
  !
  ! NOTE: after revision 680 removed code for vcode=5005 when rewriting for increased efficiency (M. Buehner)
  !
  !Author  : M. Buehner, February 2014
  !
  !Revision 001 : M. Bani Shahabadi, October 2018
  !          - adaptation of GPSRO calculation of delAL
  !
  implicit none

  type(struct_columnData) :: column,columng
  type(struct_obs)        :: obsSpaceData

  integer :: columnIndex,lev_M,lev_T,nlev_M,nlev_T,status,Vcode_anl
  real(8) :: ScaleFactorBottom, ScaleFactorTop
  real(8), allocatable :: delThick(:)
  real(8), pointer     :: alt_T_ptr(:),alt_M_ptr(:)
  real(8), allocatable :: delAlt_M(:),delAlt_T(:)
  real(8), pointer     :: delAlt_M_ptr(:),delAlt_T_ptr(:),delTT(:),delHU(:),delP0(:)
  type(struct_vco), pointer :: vco_anl

  vco_anl => col_getVco(columng)
  status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)

  nlev_T = col_getNumLev(columng,'TH')
  nlev_M = col_getNumLev(columng,'MM')

  allocate(delThick(0:nlev_T))
  allocate(delAlt_M(nlev_M))
  allocate(delAlt_T(nlev_T))

  ! generate the height coefficients
  call calcAltitudeCoeff(columng,obsSpaceData)

!$OMP PARALLEL DO PRIVATE(columnIndex,delAlt_M,delAlt_T,delThick,delTT,delHU,delP0,lev_M,lev_T,delAlt_M_ptr,delAlt_T_ptr,alt_M_ptr,alt_T_ptr,ScaleFactorBottom,ScaleFactorTop)
  do columnIndex = 1, col_getNumCol(columng)

    alt_M_ptr => col_getColumn(columng,columnIndex,'GZ','MM')
    alt_T_ptr => col_getColumn(columng,columnIndex,'GZ','TH')

    delAlt_M_ptr => col_getColumn(column,columnIndex,'GZ','MM')
    delAlt_T_ptr => col_getColumn(column,columnIndex,'GZ','TH')
    delTT      => col_getColumn(column,columnIndex,'TT')
    delHU      => col_getColumn(column,columnIndex,'HU')
    delP0      => col_getColumn(column,columnIndex,'P0')

    delAlt_M(:) = delAlt_M_ptr(:)
    delAlt_T(:) = delAlt_T_ptr(:)

    ! adjoint of compute height increment on momentum levels using weighted average of height increment of thermo levels
    if (Vcode_anl == 5002) then

      do lev_T = 2, (nlev_T-1)
        lev_M = lev_T ! momentum level just below thermo level being computed
        ScaleFactorBottom = (alt_T_ptr(lev_T) - alt_M_ptr(lev_M-1)) / &
                            (alt_M_ptr(lev_M) - alt_M_ptr(lev_M-1))
        ScaleFactorTop    = 1 - ScaleFactorBottom
        delAlt_M(lev_M-1) = delAlt_M(lev_M-1) + ScaleFactorTop     * delAlt_T(lev_T)
        delAlt_M(lev_M  ) = delAlt_M(lev_M  ) + ScaleFactorBottom  * delAlt_T(lev_T)
      end do

      ! adjoint of compute height increment for top thermo level (from top momentum level)
      delAlt_M(1) = delAlt_M(1) + delAlt_T(1)
      delTT(1) = delTT(1) + coeff_T_TT   (columnIndex) * delAlt_T(1)
      delHU(1) = delHU(1) + coeff_T_HU   (columnIndex) * delAlt_T(1)
      delP0(1) = delP0(1) + coeff_T_P0   (columnIndex) * delAlt_T(1) + & 
                            coeff_T_P0_dP(columnIndex) * delAlt_T(1)

      ! adjoint of compute height increment on momentum levels above the surface
      delThick(0:1) = 0.0D0
      do lev_M = 1, (nlev_M-1)
        lev_T = lev_M + 1 ! thermo level just below momentum level being computed
        delThick(lev_T) = delThick(lev_T-1) + delAlt_M(lev_M)
      end do

      ! adjoint of compute increment to thickness for each layer between the two momentum levels
      do lev_T = 2, nlev_T-1
        delTT(lev_T) = delTT(lev_T) + coeff_M_TT   (lev_T,columnIndex) * delThick(lev_T)
        delHU(lev_T) = delHU(lev_T) + coeff_M_HU   (lev_T,columnIndex) * delThick(lev_T)
        delP0(1)     = delP0(1)     + coeff_M_P0   (lev_T,columnIndex) * delThick(lev_T) + &
                                      coeff_M_P0_dP(lev_T,columnIndex) * delThick(lev_T)
      end do

    elseif (Vcode_anl == 5005) then

      do lev_T = 1, nlev_T-1
        lev_M = lev_T + 1 ! momentum level just below thermo level being computed
        ScaleFactorBottom = (alt_T_ptr(lev_T) - alt_M_ptr(lev_M-1)) / &
                            (alt_M_ptr(lev_M) - alt_M_ptr(lev_M-1))
        ScaleFactorTop    = 1 - ScaleFactorBottom
        delAlt_M(lev_M-1) = delAlt_M(lev_M-1) + ScaleFactorTop     * delAlt_T(lev_T)
        delAlt_M(lev_M  ) = delAlt_M(lev_M  ) + ScaleFactorBottom  * delAlt_T(lev_T)
      end do

      ! adjoint of compute height increment on momentum levels above the surface
      delThick(0) = 0.0D0
      do lev_M = 1, nlev_M-1
        lev_T = lev_M ! thermo level just below momentum level being computed
        delThick(lev_T) = delThick(lev_T-1) + delAlt_M(lev_M)
      end do

      ! adjoint of compute increment to thickness for each layer between the two thermo levels
      do lev_T = 1, nlev_T-1
        delTT(lev_T) = delTT(lev_T) + coeff_M_TT   (lev_T,columnIndex) * delThick(lev_T)
        delHU(lev_T) = delHU(lev_T) + coeff_M_HU   (lev_T,columnIndex) * delThick(lev_T)
        delP0(1)     = delP0(1)     + coeff_M_P0   (lev_T,columnIndex) * delThick(lev_T) + &
                                      coeff_M_P0_dP(lev_T,columnIndex) * delThick(lev_T)
      end do

    end if

  end do
!$OMP END PARALLEL DO

  deallocate(delThick)
  deallocate(delAlt_M)
  deallocate(delAlt_T)

end subroutine tt2phi_ad


subroutine calcAltitudeCoeff(columng,obsSpaceData)
  !
  !**s/r calcAltitudeCoeff - Calculating the coefficients of height for tt2phi_tl/tt2phi_ad
  !
  !Author  : M. Bani Shahabadi, Oct 2018
  !          - based on the original tt2phi_tl/tt2phi_ad by Mark Buehner 
  !
  implicit none

  type(struct_columnData) :: columng
  type(struct_obs)        :: obsSpaceData

  integer :: columnIndex,lev_M,lev_T,nlev_M,nlev_T,status,Vcode_anl
  real(8) :: hu,tt,Pr,cmp,cmp_TT,cmp_HU,cmp_P0,delLnP_M1,delLnP_T1, ScaleFactorBottom, ScaleFactorTop, ratioP1
  real(8), allocatable :: delLnP_M(:)
  real(8), pointer     :: alt_T_ptr(:),alt_M_ptr(:)
  real(8), pointer     :: delAlt_M_ptr(:),delAlt_T_ptr(:),delTT(:),delHU(:),delP0(:)
  real(8) :: rLat
  real(8) :: h0, Rgh, sLat, cLat
  type(struct_vco), pointer :: vco_anl

  logical, save :: firstTimeAltCoeff = .true.

  if ( .not. firstTimeAltCoeff ) return

  Write(*,*) "Entering subroutine calcAltitudeCoeff"

  ! initialize and save coefficients for increased efficiency (assumes no relinearization)
  firstTimeAltCoeff = .false.

  vco_anl => col_getVco(columng)
  status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)

  nlev_T = col_getNumLev(columng,'TH')
  nlev_M = col_getNumLev(columng,'MM')

  allocate(delLnP_M(nlev_M))

  ! saved arrays
  allocate(coeff_M_TT   (nlev_T,col_getNumCol(columng)))
  allocate(coeff_M_HU   (nlev_T,col_getNumCol(columng)))
  allocate(coeff_M_P0   (nlev_T,col_getNumCol(columng)))
  allocate(coeff_M_P0_dP(nlev_T,col_getNumCol(columng)))
  allocate(coeff_T_TT   (col_getNumCol(columng)))
  allocate(coeff_T_HU   (col_getNumCol(columng)))
  allocate(coeff_T_P0   (col_getNumCol(columng)))
  allocate(coeff_T_P0_dP(col_getNumCol(columng)))

  coeff_M_TT(:,:) = 0.0D0
  coeff_M_HU(:,:) = 0.0D0
  coeff_M_P0(:,:) = 0.0D0
  coeff_M_P0_dP(:,:) = 0.0D0
  coeff_T_TT(:) = 0.0D0  
  coeff_T_HU(:) = 0.0D0  
  coeff_T_P0(:) = 0.0D0
  coeff_T_P0_dP(:) = 0.0D0  

  do columnIndex = 1, col_getNumCol(columng)

    alt_T_ptr => col_getColumn(columng,columnIndex,'GZ','TH')

    ! latitude
    rLat = obs_headElem_r(obsSpaceData,OBS_LAT,columnIndex)
    sLat = sin(rLat)
    cLat = cos(rLat)

    ! compute momentum level properties
    do lev_M = 1,nlev_M
      delLnP_M(lev_M) = col_getPressureDeriv(columng,lev_M,columnIndex,'MM')/  &
                        col_getPressure(columng,lev_M,columnIndex,'MM')
    end do

    if (Vcode_anl == 5002) then

      do lev_T = 2, (nlev_T-1)
        ratioP1 = log( col_getPressure(columng,lev_T  ,columnIndex,'MM') /  &
                       col_getPressure(columng,lev_T-1,columnIndex,'MM') )

        hu = col_getElem(columng,lev_T,columnIndex,'HU')
        tt = col_getElem(columng,lev_T,columnIndex,'TT')
        Pr = col_getPressure(columng,lev_T,columnIndex,'TH')

        cmp = gpscompressibility(Pr,tt,hu)
        cmp_TT = gpscompressibility_TT(Pr,tt,hu)
        cmp_HU = gpscompressibility_HU(Pr,tt,hu)
        cmp_P0 = gpscompressibility_P0(Pr,tt,hu,col_getPressureDeriv(columng,lev_T,columnIndex,'TH'))

        ! Gravity acceleration 
        h0  = alt_T_ptr(lev_T)
        Rgh = phf_gravityalt(sLat, h0)

        coeff_M_TT(lev_T,columnIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * (fottva(hu,1.0D0) * cmp + fotvt8(tt,hu) * cmp_TT) * ratioP1 
        coeff_M_HU(lev_T,columnIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * (folnqva(hu,tt,1.0d0) / hu * cmp + fotvt8(tt,hu) * cmp_HU) * ratioP1
        coeff_M_P0(lev_T,columnIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp * (delLnP_M(lev_T) - delLnP_M(lev_T-1))
        coeff_M_P0_dP(lev_T,columnIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0 * ratioP1
      end do

      ! compute the property of the top layer (between first momentum and first thermo level for Vcode=5002)
      ratioP1 = log( col_getPressure(columng,1,columnIndex,'MM') /  &
                     col_getPressure(columng,1,columnIndex,'TH') )

      hu = col_getElem(columng,1,columnIndex,'HU')
      tt = col_getElem(columng,1,columnIndex,'TT')
      Pr = col_getPressure(columng,1,columnIndex,'TH')

      cmp = gpscompressibility(Pr,tt,hu)
      cmp_TT = gpscompressibility_TT(Pr,tt,hu)
      cmp_HU = gpscompressibility_HU(Pr,tt,hu)
      cmp_P0 = gpscompressibility_P0(Pr,tt,hu,col_getPressureDeriv(columng,1,columnIndex,'TH'))

      delLnP_M1 = col_getPressureDeriv(columng,1,columnIndex,'MM') / &
                       col_getPressure(columng,1,columnIndex,'MM')
      delLnP_T1 = col_getPressureDeriv(columng,1,columnIndex,'TH') / &
                       col_getPressure(columng,1,columnIndex,'TH')

      coeff_T_TT   (columnIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * (fottva(hu,1.0D0) * cmp + fotvt8(tt,hu) * cmp_TT) * ratioP1
      coeff_T_HU   (columnIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * (folnqva(hu,tt,1.0d0) / hu * cmp + fotvt8(tt,hu) * cmp_HU) * ratioP1
      coeff_T_P0   (columnIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp * (delLnP_M1 - delLnP_T1)
      coeff_T_P0_dP(columnIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0 * ratioP1

    elseif (Vcode_anl == 5005) then

      do lev_T = 1, (nlev_T-1)
        ratioP1 = log( col_getPressure(columng,lev_T+1,columnIndex,'MM') /  &
                       col_getPressure(columng,lev_T  ,columnIndex,'MM') )

        hu = col_getElem(columng,lev_T,columnIndex,'HU')
        tt = col_getElem(columng,lev_T,columnIndex,'TT')
        Pr = col_getPressure(columng,lev_T,columnIndex,'TH')

        cmp = gpscompressibility(Pr,tt,hu)
        cmp_TT = gpscompressibility_TT(Pr,tt,hu)
        cmp_HU = gpscompressibility_HU(Pr,tt,hu)
        cmp_P0 = gpscompressibility_P0(Pr,tt,hu,col_getPressureDeriv(columng,lev_T,columnIndex,'TH'))

        ! Gravity acceleration 
        h0  = alt_T_ptr(lev_T)
        Rgh = phf_gravityalt(sLat, h0)

        coeff_M_TT(lev_T,columnIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * (fottva(hu,1.0D0) * cmp + fotvt8(tt,hu) * cmp_TT) * ratioP1 
        coeff_M_HU(lev_T,columnIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * (folnqva(hu,tt,1.0d0) / hu * cmp + fotvt8(tt,hu) * cmp_HU) * ratioP1
        coeff_M_P0(lev_T,columnIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp * (delLnP_M(lev_T+1) - delLnP_M(lev_T))
        coeff_M_P0_dP(lev_T,columnIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0 * ratioP1
      end do

    end if

  end do

  deallocate(delLnP_M)

  Write(*,*) "Exit subroutine calcAltitudeCoeff"

end subroutine calcAltitudeCoeff


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
  tc = t - MPC_K_C_DEGREE_OFFSET_R8
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
  tc = t - MPC_K_C_DEGREE_OFFSET_R8
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
  tc = t - MPC_K_C_DEGREE_OFFSET_R8
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
  tc = t - MPC_K_C_DEGREE_OFFSET_R8
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


end module tt2phi_mod
