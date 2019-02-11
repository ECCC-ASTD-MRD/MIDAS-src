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
  use gridstatevector_mod
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

  real(8), allocatable, save :: coeff_M_TT_gsv(:,:,:,:), coeff_M_HU_gsv(:,:,:,:)
  real(8), allocatable, save :: coeff_T_TT_gsv(:,:,:),   coeff_T_HU_gsv(:,:,:)
  real(8), allocatable, save :: coeff_M_P0_delPM(:,:,:,:), coeff_M_P0_dP_delPT(:,:,:,:), coeff_M_P0_dP_delP0(:,:,:,:)
  real(8), allocatable, save :: coeff_T_P0_delP1(:,:,:),   coeff_T_P0_dP_delPT(:,:,:),   coeff_T_P0_dP_delP0(:,:,:)

  ! interface for computing GZ in column/gridstatevector_trial
  interface tt2phi
    module procedure tt2phi_col
    module procedure tt2phi_gsv
  end interface tt2phi

  ! interface for computing GZ increment in column/gridstatevector
  interface tt2phi_tl
    module procedure tt2phi_tl_col
    module procedure tt2phi_tl_gsv
  end interface tt2phi_tl

  ! interface for adjoint of computing GZ increment in column/gridstatevector
  interface tt2phi_ad
    module procedure tt2phi_ad_col
    module procedure tt2phi_ad_gsv
  end interface tt2phi_ad

contains

subroutine tt2phi_col(columnghr,obsSpaceData,beSilent_opt)
  !
  !**s/r tt2phi_col - Temperature to geopotential transformation on GEM4 staggered levels
  !                   NOTE: we assume 
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
  if (Vcode == 5002 .and. nlev_T /= nlev_M+1) call utl_abort('tt2phi_col: nlev_T is not equal to nlev_M+1!')
  if (Vcode == 5005 .and. nlev_T /= nlev_M)   call utl_abort('tt2phi_col: nlev_T is not equal to nlev_M!')

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

    alt_M_ptr => col_getColumn(columnghr,columnIndex,'GZ_M')
    alt_T_ptr => col_getColumn(columnghr,columnIndex,'GZ_T')

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

    if (nlev_M > 1) then
      ratioP  = log(col_getPressure(columnghr,nlev_M-1,columnIndex,'MM') / &
                col_getElem(columnghr,1,columnIndex,'P0') )

      uu = col_getElem(columnghr,nlev_M,columnIndex,'UU')
      vv = col_getElem(columnghr,nlev_M,columnIndex,'VV')
      ! averaged wind speed in the layer
      aveUU = 0.5D0 * uu
      avevv = 0.5D0 * vv

      ! Gravity acceleration 
      h0  = rMT
      Rgh = phf_gravityalt(sLat,h0)
      dh  = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(nlev_T-1) * ratioP
      Rgh = phf_gravityalt(sLat, h0+0.5D0*dh)

      delThick = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(nlev_T-1) * ratioP
      alt_M(nlev_M-1) = rMT + delThick
    end if

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
      Rgh = phf_gravityalt(sLat, h0+0.5D0*dh) - Eot - Eot2

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
      aveUU = col_getElem(columnghr,1,columnIndex,'UU')
      aveVV = col_getElem(columnghr,1,columnIndex,'VV')

      ! Gravity acceleration 
      h0  = alt_M(1)
      Rgh = phf_gravityalt(sLat, h0)
      dh  = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(1) * ratioP
      Rgh = phf_gravityalt(sLat, h0+0.5D0*dh) - Eot - Eot2

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
      if (nlev_T > 1) then
        ratioP  = log(col_getPressure(columnghr,nlev_T-1,columnIndex,'TH') / &
                  col_getElem(columnghr,1,columnIndex,'P0') )

        h0  = rMT
        Rgh = phf_gravityalt(sLat,h0)
        dh  = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(nlev_T-1) * ratioP
        Rgh = phf_gravityalt(sLat, h0+0.5D0*dh)

        delThick = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(nlev_T-1) * ratioP
        alt_T(nlev_T-1) = rMT + delThick
      end if
    end if

    ! fill the height array
    alt_T_ptr(1:nlev_T) = alt_T(1:nlev_T)
    alt_M_ptr(1:nlev_M) = alt_M(1:nlev_M)

    ! remove the height offset for the diagnostic levels for backward compatibility only
    if ( .not.columnghr%addGZsfcOffset ) then
      alt_T_ptr(nlev_T) = rMT
      alt_M_ptr(nlev_M) = rMT
    end if

    !if ( columnIndex == 1 ) then
    !  write(*,*) 'MAZIAR: tt2phi_col, GZ_T='
    !  write(*,*) gz_T(:)
    !  write(*,*) 'MAZIAR: tt2phi_col, GZ_M='
    !  write(*,*) gz_M(:)
    !end if

  enddo

  deallocate(tv)
  deallocate(AL_T)
  deallocate(AL_M)

end subroutine tt2phi_col


subroutine tt2phi_gsv(statevector_trial)
  !
  !**s/r tt2phi_gsv - Temperature to geopotential transformation on GEM4 staggered levels
  !                   NOTE: we assume 
  !                     1) nlev_T = nlev_M+1 
  !                     2) GZ_T(nlev_T) = GZ_M(nlev_M), both at the surface
  !                     3) a thermo level exists at the top, higher than the highest momentum level
  !                     4) the placement of the thermo levels means that GZ_T is the average of 2 nearest GZ_M
  !                        (according to Ron and Claude)
  !
  !Author  : M. Bani Shahabadi, January 2019
  !
  implicit none

  type(struct_gsv) :: statevector_trial

  integer :: lev_M,lev_T,nlev_M,nlev_T,status,Vcode,numStep,stepIndex,latIndex,lonIndex
  real(8) :: hu, tt, Pr, cmp, delThick, tvm, ratioP, ScaleFactorBottom, ScaleFactorTop
  real(8), allocatable :: tv(:), AL_T(:), AL_M(:) 
  real(8), pointer     :: GZ_T_ptr_r8(:,:,:,:),GZ_M_ptr_r8(:,:,:,:)
  real(8), pointer     :: hu_ptr_r8(:,:,:,:),tt_ptr_r8(:,:,:,:)
  real(8), pointer     :: P_T_ptr_r8(:,:,:,:),P_M_ptr_r8(:,:,:,:)
  real(8), pointer     :: P0_ptr_r8(:,:,:,:)
  real(8), pointer     :: GZsfc_ptr_r8(:,:)
  real(8), pointer     :: UU_ptr_r8(:,:,:,:),VV_ptr_r8(:,:,:,:)
  real(4), pointer     :: GZ_T_ptr_r4(:,:,:,:),GZ_M_ptr_r4(:,:,:,:)
  real(4), pointer     :: hu_ptr_r4(:,:,:,:),tt_ptr_r4(:,:,:,:)
  real(4), pointer     :: P_T_ptr_r4(:,:,:,:),P_M_ptr_r4(:,:,:,:)
  real(4), pointer     :: P0_ptr_r4(:,:,:,:)
  real(4), pointer     :: UU_ptr_r4(:,:,:,:),VV_ptr_r4(:,:,:,:)
  real                 :: gz_sfcOffset_T_r4, gz_sfcOffset_M_r4
  real(4) :: lat_4
  real(8) :: lat_8, rMT, uu, vv, uuM1, vvM1, aveUU, aveVV
  real(8) :: h0, dh, Rgh, Eot, Eot2, sLat, cLat
  type(struct_vco), pointer :: vco_ghr

  real(8) :: P_M, P_M1, P_Mm1, P_T, P0

  call tmg_start(204,'tt2phi_gsv')

  write(*,*) 'entering tt2phi_gsv'

  vco_ghr => gsv_getVco(statevector_trial)
  status = vgd_get(vco_ghr%vgrid,key='ig_1 - vertical coord code',value=Vcode)

  nlev_T = gsv_getNumLev(statevector_trial,'TH')
  nlev_M = gsv_getNumLev(statevector_trial,'MM')
  numStep = statevector_trial%numstep

  if (Vcode == 5002 .and. nlev_T /= nlev_M+1) call utl_abort('tt2phi_gsv: nlev_T is not equal to nlev_M+1!')
  if (Vcode == 5005 .and. nlev_T /= nlev_M)   call utl_abort('tt2phi_gsv: nlev_T is not equal to nlev_M!')

  allocate(tv(nlev_T))
  allocate(AL_T(nlev_T))
  allocate(AL_M(nlev_M))

  if ( statevector_trial%dataKind == 4 ) then
    GZ_M_ptr_r4 => gsv_getField_r4(statevector_trial,'GZ_M')
    GZ_T_ptr_r4 => gsv_getField_r4(statevector_trial,'GZ_T')

    ! initialize the GZ/AL to zero
    GZ_M_ptr_r4(:,:,:,:) = 0.0
    GZ_T_ptr_r4(:,:,:,:) = 0.0
  else
    GZ_M_ptr_r8 => gsv_getField_r8(statevector_trial,'GZ_M')
    GZ_T_ptr_r8 => gsv_getField_r8(statevector_trial,'GZ_T')

    ! initialize the GZ/AL to zero
    GZ_M_ptr_r8(:,:,:,:) = 0.0d0
    GZ_T_ptr_r8(:,:,:,:) = 0.0d0
  end if


  if ( statevector_trial%dataKind == 4 ) then
    hu_ptr_r4 => gsv_getField_r4(statevector_trial,'HU')
    tt_ptr_r4 => gsv_getField_r4(statevector_trial,'TT')
    P_T_ptr_r4 => gsv_getField_r4(statevector_trial,'P_T')
    P_M_ptr_r4 => gsv_getField_r4(statevector_trial,'P_M')
    P0_ptr_r4 => gsv_getField_r4(statevector_trial,'P0')
    UU_ptr_r4 => gsv_getField_r4(statevector_trial,'UU')
    VV_ptr_r4 => gsv_getField_r4(statevector_trial,'VV')
  else
    hu_ptr_r8 => gsv_getField_r8(statevector_trial,'HU')
    tt_ptr_r8 => gsv_getField_r8(statevector_trial,'TT')
    P_T_ptr_r8 => gsv_getField_r8(statevector_trial,'P_T')
    P_M_ptr_r8 => gsv_getField_r8(statevector_trial,'P_M')
    P0_ptr_r8 => gsv_getField_r8(statevector_trial,'P0')
    UU_ptr_r8 => gsv_getField_r8(statevector_trial,'UU')
    VV_ptr_r8 => gsv_getField_r8(statevector_trial,'VV')
  end if
  GZsfc_ptr_r8 => gsv_getGZsfc(statevector_trial)


  ! compute virtual temperature on thermo levels (corrected of compressibility)
  do stepIndex = 1, numStep
    do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
      do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd

        AL_T(:) = 0.0D0
        AL_M(:) = 0.0D0

        ! latitude
        lat_4 = statevector_trial%hco%lat2d_4(lonIndex,latIndex)
        lat_8 = real(lat_4,8)
        sLat = sin(lat_8)
        cLat = cos(lat_8)

        do lev_T = 1, nlev_T
          if ( statevector_trial%dataKind == 4 ) then
            hu = real(hu_ptr_r4(lonIndex,latIndex,lev_T,stepIndex),8)
            tt = real(tt_ptr_r4(lonIndex,latIndex,lev_T,stepIndex),8)
            Pr = real(P_T_ptr_r4(lonIndex,latIndex,lev_T,stepIndex),8)
          else
            hu = hu_ptr_r8(lonIndex,latIndex,lev_T,stepIndex)
            tt = tt_ptr_r8(lonIndex,latIndex,lev_T,stepIndex)
            Pr = P_T_ptr_r8(lonIndex,latIndex,lev_T,stepIndex)
          end if
          cmp = gpscompressibility(Pr,tt,hu)

          tv(lev_T) = fotvt8(tt,hu) * cmp 
        enddo

        rMT = GZsfc_ptr_r8(lonIndex,latIndex)

        ! compute altitude on bottom momentum level
        if (Vcode == 5002) then
          AL_M(nlev_M) = rMT

          uuM1 = 0.0D0
          vvM1 = 0.0D0
        elseif (Vcode == 5005) then
          if ( statevector_trial%dataKind == 4 ) then
            P_M = real(P_M_ptr_r4(lonIndex,latIndex,nlev_M,stepIndex),8)
            P0 = real(P0_ptr_r4(lonIndex,latIndex,1,stepIndex),8)

            uu = real(UU_ptr_r4(lonIndex,latIndex,nlev_M,stepIndex),8)
            vv = real(VV_ptr_r4(lonIndex,latIndex,nlev_M,stepIndex),8)
          else
            P_M = P_M_ptr_r8(lonIndex,latIndex,nlev_M,stepIndex)
            P0 = P0_ptr_r8(lonIndex,latIndex,1,stepIndex)

            uu = UU_ptr_r8(lonIndex,latIndex,nlev_M,stepIndex)
            vv = VV_ptr_r8(lonIndex,latIndex,nlev_M,stepIndex)
          end if

          ratioP  = log( P_M / P0 )

          ! averaged wind speed in the layer
          aveUU = 0.5D0 * uu
          avevv = 0.5D0 * vv

          ! Gravity acceleration 
          h0  = rMT
          Eot = 2 * WGS_OmegaPrime * cLat * aveUU
          Eot2= (aveUU ** 2 + aveVV ** 2) / WGS_a
          Rgh = phf_gravityalt(sLat,h0) - Eot - Eot2
          dh  = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(nlev_T) * ratioP
          Rgh = phf_gravityalt(sLat, h0+0.5D0*dh) - Eot - Eot2

          delThick = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(nlev_T) * ratioP
          AL_M(nlev_M) = rMT + delThick

          uuM1 = uu 
          vvM1 = vv
        endif ! VCode 

        ! compute altitude on rest of momentum levels
        do lev_M = nlev_M-1, 1, -1
          if ( statevector_trial%dataKind == 4 ) then
            P_M = real(P_M_ptr_r4(lonIndex,latIndex,lev_M,stepIndex),8)
            P_M1 = real(P_M_ptr_r4(lonIndex,latIndex,lev_M+1,stepIndex),8)

            uu = real(UU_ptr_r4(lonIndex,latIndex,lev_M,stepIndex),8)
            vv = real(VV_ptr_r4(lonIndex,latIndex,lev_M,stepIndex),8)
          else
            P_M = P_M_ptr_r8(lonIndex,latIndex,lev_M,stepIndex)
            P_M1 = P_M_ptr_r8(lonIndex,latIndex,lev_M+1,stepIndex)

            uu = UU_ptr_r8(lonIndex,latIndex,lev_M,stepIndex)
            vv = VV_ptr_r8(lonIndex,latIndex,lev_M,stepIndex)
          end if

          ratioP  = log( P_M / P_M1 )

          ! averaged wind speed in the layer
          aveUU = 0.5D0 * (uuM1 + uu)
          avevv = 0.5D0 * (vvM1 + vv)
          
          if (Vcode == 5002) then
            lev_T = lev_M + 1
          elseif (Vcode == 5005) then
            lev_T = lev_M
          endif

          ! Gravity acceleration 
          h0  = AL_M(lev_M+1)
          Eot = 2 * WGS_OmegaPrime * cLat * aveUU
          Eot2= (aveUU ** 2 + aveVV ** 2) / WGS_a
          Rgh = phf_gravityalt(sLat,h0) - Eot - Eot2
          dh  = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(lev_T) * ratioP
          Rgh = phf_gravityalt(sLat, h0+0.5D0*dh) - Eot - Eot2

          delThick   = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(lev_T) * ratioP
          AL_M(lev_M) = AL_M(lev_M+1) + delThick

          uuM1 = uu
          vvM1 = vv
        enddo ! lev_M

        ! compute Altitude on thermo levels
        if (Vcode == 5002) then
          AL_T(nlev_T) = AL_M(nlev_M)

          do lev_T = 2, nlev_T-1
            lev_M = lev_T ! momentum level just below thermo level being computed

            if ( statevector_trial%dataKind == 4 ) then
              P_T   = real(P_T_ptr_r4(lonIndex,latIndex,lev_T  ,stepIndex),8)
              P_M   = real(P_M_ptr_r4(lonIndex,latIndex,lev_M  ,stepIndex),8)
              P_Mm1 = real(P_M_ptr_r4(lonIndex,latIndex,lev_M-1,stepIndex),8)
            else
              P_T   = P_T_ptr_r8(lonIndex,latIndex,lev_T  ,stepIndex)
              P_M   = P_M_ptr_r8(lonIndex,latIndex,lev_M  ,stepIndex)
              P_Mm1 = P_M_ptr_r8(lonIndex,latIndex,lev_M-1,stepIndex)
            end if

            ScaleFactorBottom = log( P_T / P_Mm1 ) / log( P_M / P_Mm1 )
            ScaleFactorTop    = 1 - ScaleFactorBottom
            AL_T(lev_T) = ScaleFactorBottom * AL_M(lev_M) + ScaleFactorTop * AL_M(lev_M-1)
          end do

          ! compute altitude on top thermo level
          if ( statevector_trial%dataKind == 4 ) then
            P_T   = real(P_T_ptr_r4(lonIndex,latIndex,1,stepIndex),8)
            P_M   = real(P_M_ptr_r4(lonIndex,latIndex,1,stepIndex),8)

            aveUU = real(UU_ptr_r4(lonIndex,latIndex,1,stepIndex),8)
            aveVV = real(VV_ptr_r4(lonIndex,latIndex,1,stepIndex),8)
          else
            P_T   = P_T_ptr_r8(lonIndex,latIndex,1,stepIndex)
            P_M   = P_M_ptr_r8(lonIndex,latIndex,1,stepIndex)

            aveUU = UU_ptr_r8(lonIndex,latIndex,1,stepIndex)
            aveVV = VV_ptr_r8(lonIndex,latIndex,1,stepIndex)
          end if

          ratioP = log( P_T / P_M )

          ! Gravity acceleration 
          h0  = AL_M(1)
          Eot = 2 * WGS_OmegaPrime * cLat * aveUU
          Eot2= (aveUU ** 2 + aveVV ** 2) / WGS_a
          Rgh = phf_gravityalt(sLat, h0) - Eot - Eot2
          dh  = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(1) * ratioP
          Rgh = phf_gravityalt(sLat, h0+0.5D0*dh) - Eot - Eot2

          delThick   = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(1) * ratioP
          AL_T(1) = AL_M(1) + delThick

        elseif (Vcode == 5005) then

          do lev_T = 1, nlev_T-1
            lev_M = lev_T + 1  ! momentum level just below thermo level being computed

            if ( statevector_trial%dataKind == 4 ) then
              P_T   = real(P_T_ptr_r4(lonIndex,latIndex,lev_T  ,stepIndex),8)
              P_M   = real(P_M_ptr_r4(lonIndex,latIndex,lev_M  ,stepIndex),8)
              P_Mm1 = real(P_M_ptr_r4(lonIndex,latIndex,lev_M-1,stepIndex),8)
            else
              P_T   = P_T_ptr_r8(lonIndex,latIndex,lev_T  ,stepIndex)
              P_M   = P_M_ptr_r8(lonIndex,latIndex,lev_M  ,stepIndex)
              P_Mm1 = P_M_ptr_r8(lonIndex,latIndex,lev_M-1,stepIndex)
            end if

            ScaleFactorBottom = log( P_T / P_Mm1 ) / log( P_M / P_Mm1 )
            ScaleFactorTop    = 1 - ScaleFactorBottom
            AL_T(lev_T) = ScaleFactorBottom * AL_M(lev_M) + ScaleFactorTop * AL_M(lev_M-1)
          end do


          ! compute altitude on bottom thermo level
          if ( statevector_trial%dataKind == 4 ) then
            P_T = real(P_T_ptr_r4(lonIndex,latIndex,nlev_T,stepIndex),8)
            P0  = real(P0_ptr_r4 (lonIndex,latIndex,1     ,stepIndex),8)
          else
            P_T = P_T_ptr_r8(lonIndex,latIndex,nlev_T,stepIndex)
            P0  = P0_ptr_r8 (lonIndex,latIndex,1     ,stepIndex)
          end if

          ratioP = log( P_T / P0 )

          h0  = rMT
          Rgh = phf_gravityalt(sLat,h0)
          dh  = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(nlev_T) * ratioP
          Rgh = phf_gravityalt(sLat, h0+0.5D0*dh)

          delThick = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(nlev_T) * ratioP
          AL_T(nlev_T) = rMT + delThick
        endif

        ! compute GZ 
        if ( statevector_trial%dataKind == 4 ) then
          do lev_T = 1, nlev_T
            GZ_T_ptr_r4(lonIndex,latIndex,lev_T,stepIndex) = real(AL_T(lev_T),4)
          end do
          do lev_M = 1, nlev_M
            GZ_M_ptr_r4(lonIndex,latIndex,lev_M,stepIndex) = real(AL_M(lev_M),4)
          end do
        else
          GZ_T_ptr_r8(lonIndex,latIndex,1:nlev_T,stepIndex) = AL_T(1:nlev_T)
          GZ_M_ptr_r8(lonIndex,latIndex,1:nlev_M,stepIndex) = AL_M(1:nlev_M)
        end if

        ! remove the height offset for the diagnostic levels for backward compatibility only
        if ( .not.statevector_trial%addGZsfcOffset ) then
          if ( statevector_trial%dataKind == 4 ) then
            GZ_T_ptr_r4(lonIndex,latIndex,nlev_T,stepIndex) = real(rMT,4)
            GZ_M_ptr_r4(lonIndex,latIndex,nlev_M,stepIndex) = real(rMT,4)
          else
            GZ_T_ptr_r8(lonIndex,latIndex,nlev_T,stepIndex) = rMT
            GZ_M_ptr_r8(lonIndex,latIndex,nlev_M,stepIndex) = rMT
          end if
        end if

      enddo
    enddo
  enddo

  deallocate(tv)
  deallocate(alt_T)
  deallocate(alt_M)

  if ( statevector_trial%dataKind == 4 ) then
    write(*,*) 'tt2phi_gsv, GZ_T='
    write(*,*) GZ_T_ptr_r4(statevector_trial%myLonBeg,statevector_trial%myLatBeg,:,1)
    write(*,*) 'tt2phi_gsv, GZ_M='
    write(*,*) GZ_M_ptr_r4(statevector_trial%myLonBeg,statevector_trial%myLatBeg,:,1)
  else
    write(*,*) 'tt2phi_gsv, GZ_T='
    write(*,*) GZ_T_ptr_r8(statevector_trial%myLonBeg,statevector_trial%myLatBeg,:,1)
    write(*,*) 'tt2phi_gsv, GZ_M='
    write(*,*) GZ_M_ptr_r8(statevector_trial%myLonBeg,statevector_trial%myLatBeg,:,1)
  end if

  write(*,*) 'tt2phi_gsv: statevector_trial%addGZsfcOffset=', statevector_trial%addGZsfcOffset 
  write(*,*) 'exiting tt2phi_gsv'

  call tmg_stop(204)

end subroutine tt2phi_gsv


subroutine tt2phi_tl_col(column,columng,obsSpaceData)
  !
  !**s/r tt2phi_tl_col - Temperature to geopotential transformation on GEM4 staggered levels
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

    alt_M_ptr => col_getColumn(columng,columnIndex,'GZ_M')
    alt_T_ptr => col_getColumn(columng,columnIndex,'GZ_T')

    delAlt_M_ptr => col_getColumn(column,columnIndex,'GZ_M')
    delAlt_T_ptr => col_getColumn(column,columnIndex,'GZ_T')
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

end subroutine tt2phi_tl_col


subroutine tt2phi_tl_gsv(statevector,statevector_trial)
  !
  !**s/r tt2phi_tl_gsv- temperature to geopotential transformation on gridstatevector
  !
  !Author  : M. Bani Shahabadi, September 2018
  !
  implicit none

  type(struct_gsv) :: statevector,statevector_trial

  integer :: lev_M,lev_T,nlev_M,nlev_T,status,Vcode_anl,numStep,stepIndex,latIndex,lonIndex
  real(8) :: ScaleFactorBottom, ScaleFactorTop
  real(8), allocatable :: delThick(:,:,:,:)
  real(8), pointer     :: delGz_M(:,:,:,:),delGz_T(:,:,:,:),delTT(:,:,:,:),delHU(:,:,:,:),delP0(:,:,:,:)
  real(8), pointer     :: gz_T(:,:,:,:),gz_M(:,:,:,:)
  real(8), pointer     :: P_T(:,:,:,:), P_M(:,:,:,:)
  real(8), pointer     :: delP_T(:,:,:,:), delP_M(:,:,:,:)
  type(struct_vco), pointer :: vco_anl

  call tmg_start(201,'tt2phi_tl_gsv')

  write(*,*) 'entering tt2phi_tl_gsv'

  vco_anl => gsv_getVco(statevector_trial)
  status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)

  nlev_T = gsv_getNumLev(statevector_trial,'TH')
  nlev_M = gsv_getNumLev(statevector_trial,'MM')
  numStep = statevector_trial%numstep

  allocate(delThick(statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
                    statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
                    nlev_T,numStep))

  ! generate the height coefficients on the grid
  call calcAltitudeCoeff_gsv(statevector_trial)

  ! loop over all lat/lon/step
!!$OMP PARALLEL DO PRIVATE(columnIndex,delGz_M,delGz_T,delThick,delTT,delHU,delP0,lev_M,lev_T)

  gz_M => gsv_getField_r8(statevector_trial,'GZ_M')
  gz_T => gsv_getField_r8(statevector_trial,'GZ_T')
  P_T => gsv_getField_r8(statevector_trial,'P_T')
  P_M => gsv_getField_r8(statevector_trial,'P_M')

  delGz_M => gsv_getField_r8(statevector,'GZ_M')
  delGz_T => gsv_getField_r8(statevector,'GZ_T')
  delTT => gsv_getField_r8(statevector,'TT')
  delHU => gsv_getField_r8(statevector,'HU')
  delP0 => gsv_getField_r8(statevector,'P0')
  delP_T => gsv_getField_r8(statevector,'P_T')
  delP_M => gsv_getField_r8(statevector,'P_M')

  ! ensure increment at sfc is zero (fixed height level)
  delGz_M(:,:,nlev_M,:) = 0.0d0
  delGz_T(:,:,nlev_T,:) = 0.0d0

  if (Vcode_anl .eq. 5002) then

    ! compute increment to thickness for each layer between the two momentum levels
    do stepIndex = 1, numStep
      do lev_T = 2, (nlev_T-1)
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd

            delThick(lonIndex,latIndex,lev_T,stepIndex) = coeff_M_TT_gsv(lonIndex,latIndex,lev_T,stepIndex) * delTT(lonIndex,latIndex,lev_T,stepIndex) + &
                                                          coeff_M_HU_gsv(lonIndex,latIndex,lev_T,stepIndex) * delHU(lonIndex,latIndex,lev_T,stepIndex) + &

                                                          coeff_M_P0_delPM (lonIndex,latIndex,lev_T,stepIndex) * &
                                                          ( delP_M(lonIndex,latIndex,lev_T  ,stepIndex) / P_M(lonIndex,latIndex,lev_T  ,stepIndex) - &
                                                            delP_M(lonIndex,latIndex,lev_T-1,stepIndex) / P_M(lonIndex,latIndex,lev_T-1,stepIndex) ) + &

                                                          coeff_M_P0_dP_delPT(lonIndex,latIndex,lev_T,stepIndex) * delP_T(lonIndex,latIndex,lev_T,stepIndex) + &
                                                          coeff_M_P0_dP_delP0(lonIndex,latIndex,lev_T,stepIndex) * delP0(lonIndex,latIndex,1,stepIndex)
          enddo
        enddo
      enddo
    enddo

    ! compute GZ increment on momentum levels above the surface
    do stepIndex = 1, numStep
      do lev_M = (nlev_M-1), 1, -1
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            lev_T = lev_M + 1 ! thermo level just below momentum level being computed
            delGz_M(lonIndex,latIndex,lev_M,stepIndex) = delGz_M(lonIndex,latIndex,lev_M+1,stepIndex) + delThick(lonIndex,latIndex,lev_T,stepIndex)
          enddo
        enddo
      enddo
    enddo

    ! compute GZ increment on thermo levels using weighted average of GZ increment of momentum levels
    do stepIndex = 1, numStep
      do lev_T = 1, (nlev_T-1)
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            if ( lev_T == 1) then
              ! compute GZ increment for top thermo level (from top momentum level)
              delGz_T(lonIndex,latIndex,1,stepIndex) = delGz_M(lonIndex,latIndex,1,stepIndex) +  &
                                 coeff_T_TT_gsv(lonIndex,latIndex,stepIndex) * delTT(lonIndex,latIndex,1,stepIndex) + &
                                 coeff_T_HU_gsv(lonIndex,latIndex,stepIndex) * delHU(lonIndex,latIndex,1,stepIndex) + &

                                 coeff_T_P0_delP1(lonIndex,latIndex,stepIndex) * &
                                 ( delP_M(lonIndex,latIndex,1,stepIndex) / P_M(lonIndex,latIndex,1,stepIndex) - &
                                   delP_T(lonIndex,latIndex,1,stepIndex) / P_T(lonIndex,latIndex,1,stepIndex) ) + &

                                 coeff_T_P0_dP_delPT(lonIndex,latIndex,stepIndex) * delP_T(lonIndex,latIndex,1,stepIndex) + &
                                 coeff_T_P0_dP_delP0(lonIndex,latIndex,stepIndex) * delP0(lonIndex,latIndex,1,stepIndex)

            else
              lev_M = lev_T ! momentum level just below thermo level being computed
              ScaleFactorBottom = (gz_T(lonIndex,latIndex,lev_T,stepIndex) - gz_M(lonIndex,latIndex,lev_M-1,stepIndex)) / &
                                  (gz_M(lonIndex,latIndex,lev_M,stepIndex) - gz_M(lonIndex,latIndex,lev_M-1,stepIndex))
              ScaleFactorTop    = 1 - ScaleFactorBottom
              delGz_T(lonIndex,latIndex,lev_T,stepIndex) = ScaleFactorBottom * delGz_M(lonIndex,latIndex,lev_M  ,stepIndex) + &
                                                              ScaleFactorTop * delGz_M(lonIndex,latIndex,lev_M-1,stepIndex)
            endif

          enddo
        enddo
      enddo
    enddo

  elseif(Vcode_anl .eq. 5005) then

    ! compute increment to thickness for each layer between the two momentum levels
    do stepIndex = 1, numStep
      do lev_T = 1, (nlev_T-1)
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd

            delThick(lonIndex,latIndex,lev_T,stepIndex) = coeff_M_TT_gsv(lonIndex,latIndex,lev_T,stepIndex) * delTT(lonIndex,latIndex,lev_T,stepIndex) + &
                                                          coeff_M_HU_gsv(lonIndex,latIndex,lev_T,stepIndex) * delHU(lonIndex,latIndex,lev_T,stepIndex) + &

                                                          coeff_M_P0_delPM (lonIndex,latIndex,lev_T,stepIndex) * &
                                                          ( delP_M(lonIndex,latIndex,lev_T+1,stepIndex) / P_M(lonIndex,latIndex,lev_T+1,stepIndex) - &
                                                            delP_M(lonIndex,latIndex,lev_T  ,stepIndex) / P_M(lonIndex,latIndex,lev_T  ,stepIndex) ) + &

                                                          coeff_M_P0_dP_delPT(lonIndex,latIndex,lev_T,stepIndex) * delP_T(lonIndex,latIndex,lev_T,stepIndex) + &
                                                          coeff_M_P0_dP_delP0(lonIndex,latIndex,lev_T,stepIndex) * delP0(lonIndex,latIndex,1,stepIndex)
          enddo
        enddo
      enddo
    enddo

    ! compute GZ increment on momentum levels above the surface
    do stepIndex = 1, numStep
      do lev_M = (nlev_M-1), 1, -1
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            lev_T = lev_M ! thermo level just below momentum level being computed
            delGz_M(lonIndex,latIndex,lev_M,stepIndex) = delGz_M(lonIndex,latIndex,lev_M+1,stepIndex) + delThick(lonIndex,latIndex,lev_T,stepIndex)
          enddo
        enddo
      enddo
    enddo

    ! compute GZ increment on thermo levels using weighted average of GZ increment of momentum levels
    do stepIndex = 1, numStep
      do lev_T = 1, (nlev_T-1)
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            lev_M = lev_T + 1 ! momentum level just below thermo level being computed
            ScaleFactorBottom = (gz_T(lonIndex,latIndex,lev_T,stepIndex) - gz_M(lonIndex,latIndex,lev_M-1,stepIndex)) / &
                                (gz_M(lonIndex,latIndex,lev_M,stepIndex) - gz_M(lonIndex,latIndex,lev_M-1,stepIndex))
            ScaleFactorTop    = 1 - ScaleFactorBottom
            delGz_T(lonIndex,latIndex,lev_T,stepIndex) = ScaleFactorBottom * delGz_M(lonIndex,latIndex,lev_M  ,stepIndex) + &
                                                            ScaleFactorTop * delGz_M(lonIndex,latIndex,lev_M-1,stepIndex)
          enddo
        enddo
      enddo
    enddo

  endif

!!$OMP END PARALLEL DO

  deallocate(delThick)

  write(*,*) 'tt2phi_tl_gsv, delGz_T='
  write(*,*) delGz_T(statevector_trial%myLonBeg,statevector_trial%myLatBeg,:,1)
  write(*,*) 'tt2phi_tl_gsv, delGz_M='
  write(*,*) delGz_M(statevector_trial%myLonBeg,statevector_trial%myLatBeg,:,1)

  write(*,*) 'exiting tt2phi_tl_gsv'

  call tmg_stop(201)

end subroutine tt2phi_tl_gsv


subroutine tt2phi_ad_col(column,columng,obsSpaceData)
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

    alt_M_ptr => col_getColumn(columng,columnIndex,'GZ_M')
    alt_T_ptr => col_getColumn(columng,columnIndex,'GZ_T')

    delAlt_M_ptr => col_getColumn(column,columnIndex,'GZ_M')
    delAlt_T_ptr => col_getColumn(column,columnIndex,'GZ_T')
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

end subroutine tt2phi_ad_col


subroutine tt2phi_ad_gsv(statevector,statevector_trial)
  !
  !**s/r tt2phi_ad_gsv- Adjoint of temperature to geopotential transformation on gridstatevector
  !
  !Author  : M. Bani Shahabadi, September 2018
  !
  implicit none

  type(struct_gsv) :: statevector,statevector_trial

  integer :: lev_M,lev_T,nlev_M,nlev_T,status,Vcode_anl,numStep,stepIndex,latIndex,lonIndex
  real(8) :: ScaleFactorBottom, ScaleFactorTop
  real(8), allocatable :: delThick(:,:,:,:)
  real(8), pointer     :: delGz_M_in(:,:,:,:),delGz_T_in(:,:,:,:),delTT(:,:,:,:),delHU(:,:,:,:),delP0(:,:,:,:)
  real(8), pointer     :: gz_M(:,:,:,:),gz_T(:,:,:,:)
  real(8), allocatable :: delGz_M(:,:,:,:),delGz_T(:,:,:,:)
  real(8), pointer     :: P_M(:,:,:,:),P_T(:,:,:,:)
  real(8), pointer     :: delP_M(:,:,:,:),delP_T(:,:,:,:)
  type(struct_vco), pointer :: vco_anl


  call tmg_start(202,'tt2phi_ad_gsv MAZIAR')

  write(*,*) 'MAZIAR: entering tt2phi_ad_gsv'

  vco_anl => gsv_getVco(statevector_trial)
  status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)

  nlev_T = gsv_getNumLev(statevector_trial,'TH')
  nlev_M = gsv_getNumLev(statevector_trial,'MM')
  numStep = statevector_trial%numstep

  allocate(delGz_M(statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
                   statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
                   nlev_M,numStep))
  allocate(delGz_T(statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
                   statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
                   nlev_T,numStep))
  allocate(delThick(statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
                   statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
                   0:nlev_T,numStep))

  ! generate the height coefficients on the grid
  call calcAltitudeCoeff_gsv(statevector_trial)

  ! loop over all lat/lon/step
!!$OMP PARALLEL DO PRIVATE(columnIndex,delGz_M,delGz_T,delThick,delTT,delHU,delP0,lev_M,lev_T)

  gz_M => gsv_getField_r8(statevector_trial,'GZ_M')
  gz_T => gsv_getField_r8(statevector_trial,'GZ_T')
  P_T => gsv_getField_r8(statevector_trial,'P_T')
  P_M => gsv_getField_r8(statevector_trial,'P_M')

  delGz_M_in => gsv_getField_r8(statevector,'GZ_M')
  delGz_T_in => gsv_getField_r8(statevector,'GZ_T')

  delTT => gsv_getField_r8(statevector,'TT')
  delHU => gsv_getField_r8(statevector,'HU')
  delP0 => gsv_getField_r8(statevector,'P0')
  delP_T => gsv_getField_r8(statevector,'P_T')
  delP_M => gsv_getField_r8(statevector,'P_M')

  delGz_M(:,:,:,:) = delGz_M_in(:,:,:,:)
  delGz_T(:,:,:,:) = delGz_T_in(:,:,:,:)

  if(Vcode_anl .eq. 5002) then

    ! adjoint of compute GZ increment on thermo levels by simple averaging
    do stepIndex = 1, numStep
      do lev_T = 1, (nlev_T-1)
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            lev_M = lev_T ! momentum level just below thermo level being computed

            ! adjoint of compute GZ increment on top thermo level (from top momentum level)
            if (lev_T == 1) then
              delGz_M(lonIndex,latIndex,1,stepIndex)  = delGz_M(lonIndex,latIndex,1,stepIndex) + &
                                                        delGz_T(lonIndex,latIndex,1,stepIndex)

              delTT(lonIndex,latIndex,1,stepIndex) = delTT(lonIndex,latIndex,1,stepIndex) + &
                                                     coeff_T_TT_gsv   (lonIndex,latIndex,stepIndex) * delGz_T(lonIndex,latIndex,1,stepIndex)
              delHU(lonIndex,latIndex,1,stepIndex) = delHU(lonIndex,latIndex,1,stepIndex) + &
                                                     coeff_T_HU_gsv   (lonIndex,latIndex,stepIndex) * delGz_T(lonIndex,latIndex,1,stepIndex)

              delP_M(lonIndex,latIndex,1,stepIndex) = delP_M(lonIndex,latIndex,1,stepIndex) + &
                                                      coeff_T_P0_delP1(lonIndex,latIndex,stepIndex) / P_M(lonIndex,latIndex,1,stepIndex) * &
                                                      delGz_T(lonIndex,latIndex,1,stepIndex)

              delP_T(lonIndex,latIndex,1,stepIndex) = delP_T(lonIndex,latIndex,1,stepIndex) - &
                                                      coeff_T_P0_delP1(lonIndex,latIndex,stepIndex) / P_T(lonIndex,latIndex,1,stepIndex) * &
                                                      delGz_T(lonIndex,latIndex,1,stepIndex)

              delP_T(lonIndex,latIndex,1,stepIndex) = delP_T(lonIndex,latIndex,1,stepIndex) + &
                                                     coeff_T_P0_dP_delPT(lonIndex,latIndex,stepIndex) * delGz_T(lonIndex,latIndex,1,stepIndex)

              delP0(lonIndex,latIndex,1,stepIndex) = delP0(lonIndex,latIndex,1,stepIndex) + &
                                                     coeff_T_P0_dp_delP0(lonIndex,latIndex,stepIndex) * delGz_T(lonIndex,latIndex,1,stepIndex)
            else
              ScaleFactorBottom = (gz_T(lonIndex,latIndex,lev_T,stepIndex) - gz_M(lonIndex,latIndex,lev_M-1,stepIndex)) / &
                                  (gz_M(lonIndex,latIndex,lev_M,stepIndex) - gz_M(lonIndex,latIndex,lev_M-1,stepIndex))
              ScaleFactorTop    = 1 - ScaleFactorBottom
              delGz_M(lonIndex,latIndex,lev_M-1,stepIndex) = delGz_M(lonIndex,latIndex,lev_M-1,stepIndex) + &
                                            ScaleFactorTop * delGz_T(lonIndex,latIndex,lev_T  ,stepIndex)
              delGz_M(lonIndex,latIndex,lev_M  ,stepIndex) = delGz_M(lonIndex,latIndex,lev_M  ,stepIndex) + &
                                         ScaleFactorBottom * delGz_T(lonIndex,latIndex,lev_T  ,stepIndex)
            end if
          enddo
        enddo
      enddo
    enddo

    ! adjoint of compute GZ increment on momentum levels above the surface
    delThick(:,:,0:1,:) = 0.0d0
    do stepIndex = 1, numStep
      do lev_M = 1, (nlev_M-1)
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            lev_T = lev_M + 1 ! thermo level just below momentum level being computed
            delThick(lonIndex,latIndex,lev_T,stepIndex) = delThick(lonIndex,latIndex,lev_T-1,stepIndex) + &
                                                           delGz_M(lonIndex,latIndex,lev_M  ,stepIndex)
          end do
        end do
      end do
    end do

    ! adjoint of compute increment to thickness for each layer between the two momentum levels
    do stepIndex = 1, numStep
      do lev_T = 2, nlev_T-1
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            delTT(lonIndex,latIndex,lev_T,stepIndex) = delTT            (lonIndex,latIndex,lev_T,stepIndex) + &
                                                       coeff_M_TT_gsv   (lonIndex,latIndex,lev_T,stepIndex) * delThick(lonIndex,latIndex,lev_T,stepIndex)
            delHU(lonIndex,latIndex,lev_T,stepIndex) = delHU            (lonIndex,latIndex,lev_T,stepIndex) + &
                                                       coeff_M_HU_gsv   (lonIndex,latIndex,lev_T,stepIndex) * delThick(lonIndex,latIndex,lev_T,stepIndex)

            delP_M(lonIndex,latIndex,lev_T,stepIndex)= delP_M(lonIndex,latIndex,lev_T,stepIndex) + &
                                                       coeff_M_P0_delPM (lonIndex,latIndex,lev_T,stepIndex) / P_M(lonIndex,latIndex,lev_T,stepIndex) * &
                                                       delThick(lonIndex,latIndex,lev_T,stepIndex)

            delP_M(lonIndex,latIndex,lev_T-1,stepIndex) = delP_M(lonIndex,latIndex,lev_T-1,stepIndex) - &
                                                          coeff_M_P0_delPM (lonIndex,latIndex,lev_T,stepIndex) / P_M(lonIndex,latIndex,lev_T-1,stepIndex) * &
                                                          delThick(lonIndex,latIndex,lev_T,stepIndex)

            delP_T(lonIndex,latIndex,lev_T,stepIndex) = delP_T(lonIndex,latIndex,lev_T,stepIndex) + &
                                                        coeff_M_P0_dP_delPT(lonIndex,latIndex,lev_T,stepIndex) * delThick(lonIndex,latIndex,lev_T,stepIndex)

            delP0(lonIndex,latIndex,1,stepIndex) = delP0(lonIndex,latIndex,1,stepIndex) + &
                                                   coeff_M_P0_dP_delP0(lonIndex,latIndex,lev_T,stepIndex) * delThick(lonIndex,latIndex,lev_T,stepIndex)
          enddo
        enddo
      enddo
    enddo

  elseif(Vcode_anl .eq. 5005) then

    ! adjoint of compute GZ increment on thermo levels by simple averaging
    do stepIndex = 1, numStep
      do lev_T = 1, (nlev_T-1)
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            lev_M = lev_T+1 ! momentum level just below thermo level being computed
            ScaleFactorBottom = (gz_T(lonIndex,latIndex,lev_T,stepIndex) - gz_M(lonIndex,latIndex,lev_M-1,stepIndex)) / &
                                (gz_M(lonIndex,latIndex,lev_M,stepIndex) - gz_M(lonIndex,latIndex,lev_M-1,stepIndex))
            ScaleFactorTop    = 1 - ScaleFactorBottom
            delGz_M(lonIndex,latIndex,lev_M-1,stepIndex) = delGz_M(lonIndex,latIndex,lev_M-1,stepIndex) + &
                                          ScaleFactorTop * delGz_T(lonIndex,latIndex,lev_T  ,stepIndex)
            delGz_M(lonIndex,latIndex,lev_M,stepIndex)   = delGz_M(lonIndex,latIndex,lev_M  ,stepIndex) + &
                                       ScaleFactorBottom * delGz_T(lonIndex,latIndex,lev_T  ,stepIndex)
          enddo
        enddo
      enddo
    enddo

    ! adjoint of compute GZ increment on momentum levels
    delThick(:,:,0,:) = 0.0d0
    do stepIndex = 1, numStep
      do lev_M = 1, (nlev_M-1)
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            lev_T = lev_M ! thermo level just below momentum level being computed
            delThick(lonIndex,latIndex,lev_T,stepIndex) = delThick(lonIndex,latIndex,lev_T-1,stepIndex) + &
                                                          delGz_M (lonIndex,latIndex,lev_M  ,stepIndex)
          end do
        end do
      end do
    end do

    do stepIndex = 1, numStep
      do lev_T = 1, nlev_T-1
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            delTT(lonIndex,latIndex,lev_T,stepIndex) = delTT(lonIndex,latIndex,lev_T,stepIndex) + &
                                                       coeff_M_TT_gsv(lonIndex,latIndex,lev_T,stepIndex) * delThick(lonIndex,latIndex,lev_T,stepIndex)
            delHU(lonIndex,latIndex,lev_T,stepIndex) = delHU(lonIndex,latIndex,lev_T,stepIndex) + &
                                                       coeff_M_HU_gsv(lonIndex,latIndex,lev_T,stepIndex) * delThick(lonIndex,latIndex,lev_T,stepIndex)

            delP_M(lonIndex,latIndex,lev_T+1,stepIndex) = delP_M(lonIndex,latIndex,lev_T+1,stepIndex) + &
                                                          coeff_M_P0_delPM(lonIndex,latIndex,lev_T,stepIndex) / P_M(lonIndex,latIndex,lev_T+1,stepIndex) * &
                                                          delThick(lonIndex,latIndex,lev_T,stepIndex)

            delP_M(lonIndex,latIndex,lev_T  ,stepIndex) = delP_M(lonIndex,latIndex,lev_T  ,stepIndex) - &
                                                          coeff_M_P0_delPM(lonIndex,latIndex,lev_T,stepIndex) / P_M(lonIndex,latIndex,lev_T  ,stepIndex) * &
                                                          delThick(lonIndex,latIndex,lev_T,stepIndex)

            delP_T(lonIndex,latIndex,lev_T  ,stepIndex) = delP_T(lonIndex,latIndex,lev_T  ,stepIndex) + &
                                                          coeff_M_P0_dP_delPT(lonIndex,latIndex,lev_T,stepIndex) * delThick(lonIndex,latIndex,lev_T,stepIndex)

            delP0(lonIndex,latIndex,1,stepIndex)     = delP0(lonIndex,latIndex,1,stepIndex) + &
                                                       coeff_M_P0_dP_delP0(lonIndex,latIndex,lev_T,stepIndex) * delThick(lonIndex,latIndex,lev_T,stepIndex)
          enddo
        enddo
      enddo
    enddo

  endif

!!$OMP END PARALLEL DO

  !write(*,*) 'MAZIAR: tt2phi_ad_gsv, delGz_T='
  !write(*,*) delGz_T(statevector_trial%myLonBeg,statevector_trial%myLatBeg,:,1)
  !write(*,*) 'MAZIAR: tt2phi_ad_gsv, delGz_M='
  !write(*,*) delGz_M(statevector_trial%myLonBeg,statevector_trial%myLatBeg,:,1)

  deallocate(delThick)
  deallocate(delGz_M)
  deallocate(delGz_T)

  write(*,*) 'MAZIAR: exiting tt2phi_ad_gsv'

  call tmg_stop(202)

end subroutine tt2phi_ad_gsv


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

    alt_T_ptr => col_getColumn(columng,columnIndex,'GZ_T')

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


subroutine calcAltitudeCoeff_gsv(statevector_trial)
  !
  !**s/r calcAltitudeCoeff_gsv - Calculating the coefficients of height for tt2phi_tl_gsv/tt2phi_ad_gsv
  !
  !Author  : M. Bani Shahabadi, Jan 2019
  !          - based on the original calcAltitudeCoef by M. Bani Shahabadi
  !
  implicit none

  type(struct_gsv) :: statevector_trial

  integer :: lev_M,lev_T,nlev_M,nlev_T,status,Vcode_an,numStep,stepIndex,latIndex,lonIndex,Vcode_anl
  real(8) :: hu,tt,Pr,gz,cmp,cmp_TT,cmp_HU,cmp_P0_1,cmp_P0_2,ratioP1
  real(4) :: lat_4
  real(8) :: Rgh, sLat, lat_8
  real(8), pointer :: hu_ptr(:,:,:,:),tt_ptr(:,:,:,:)
  real(8), pointer :: P_T_ptr(:,:,:,:),P_M_ptr(:,:,:,:)
  real(8), pointer :: gz_T_ptr(:,:,:,:)
  type(struct_vco), pointer :: vco_anl

  logical, save :: firstTimeAltCoeff_gsv = .true.

  if ( .not. firstTimeAltCoeff_gsv ) return

  Write(*,*) "Entering subroutine calcAltitudeCoeff_gsv"

  ! initialize and save coefficients for increased efficiency (assumes no relinearization)
  firstTimeAltCoeff_gsv = .false.

  vco_anl => gsv_getVco(statevector_trial)
  status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)

  nlev_T = gsv_getNumLev(statevector_trial,'TH')
  nlev_M = gsv_getNumLev(statevector_trial,'MM')
  numStep = statevector_trial%numstep

  ! saved arrays
  allocate(coeff_M_TT_gsv     (statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
                               statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
                               nlev_T,numStep))
  allocate(coeff_M_HU_gsv     (statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
                               statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
                               nlev_T,numStep))
  allocate(coeff_M_P0_delPM   (statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
                               statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
                               nlev_T,numStep))
  allocate(coeff_M_P0_dP_delPT(statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
                               statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
                               nlev_T,numStep))
  allocate(coeff_M_P0_dP_delP0(statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
                               statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
                               nlev_T,numStep))

  allocate(coeff_T_TT_gsv     (statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
                               statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
                               numStep))
  allocate(coeff_T_HU_gsv     (statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
                               statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
                               numStep))
  allocate(coeff_T_P0_delP1   (statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
                               statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
                               numStep))
  allocate(coeff_T_P0_dP_delPT(statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
                               statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
                               numStep))
  allocate(coeff_T_P0_dP_delP0(statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
                               statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
                               numStep))

  coeff_M_TT_gsv(:,:,:,:) = 0.0D0
  coeff_M_HU_gsv(:,:,:,:) = 0.0D0

  coeff_M_P0_delPM(:,:,:,:) = 0.0D0

  coeff_M_P0_dP_delPT(:,:,:,:) = 0.0D0
  coeff_M_P0_dP_delP0(:,:,:,:) = 0.0D0

  coeff_T_TT_gsv(:,:,:) = 0.0D0  
  coeff_T_HU_gsv(:,:,:) = 0.0D0  

  coeff_T_P0_delP1(:,:,:) = 0.0D0

  coeff_T_P0_dP_delPT(:,:,:) = 0.0D0  
  coeff_T_P0_dP_delP0(:,:,:) = 0.0D0  

!!$OMP PARALLEL DO PRIVATE(columnIndex,delGz_M,delGz_T,delThick,delTT,delHU,delP0,lev_M,lev_T)
  hu_ptr => gsv_getField_r8(statevector_trial,'HU')
  tt_ptr => gsv_getField_r8(statevector_trial,'TT')
  P_T_ptr => gsv_getField_r8(statevector_trial,'P_T')
  P_M_ptr => gsv_getField_r8(statevector_trial,'P_M')
  gz_T_ptr => gsv_getField_r8(statevector_trial,'GZ_T')

  if (Vcode_anl == 5002) then

    do stepIndex = 1, numStep
      do lev_T = 1, (nlev_T-1)
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            if ( lev_T == 1 ) then 
              ! compute coefficients for GZ on only the top thermo level
              ratioP1 = log( P_M_ptr(lonIndex,latIndex,1,stepIndex) / &
                             P_T_ptr(lonIndex,latIndex,1,stepIndex) ) 
              hu = max(hu_ptr(lonIndex,latIndex,1,stepIndex),MPC_MINIMUM_HU_R8)
              tt = tt_ptr(lonIndex,latIndex,1,stepIndex)
              Pr = P_T_ptr(lonIndex,latIndex,1,stepIndex)
              gz = gz_T_ptr(lonIndex,latIndex,1,stepIndex)

              cmp = gpscompressibility(Pr,tt,hu)
              cmp_TT = gpscompressibility_TT(Pr,tt,hu)
              cmp_HU = gpscompressibility_HU(Pr,tt,hu)
              cmp_P0_1 = gpscompressibility_P0_1(Pr,tt,hu,1.0d0)
              cmp_P0_2 = gpscompressibility_P0_2(Pr,tt,hu)

              ! Gravity acceleration 
              lat_4 = statevector_trial%hco%lat2d_4(lonIndex,latIndex)
              lat_8 = real(lat_4,8)
              sLat = sin(lat_8)
              Rgh = phf_gravityalt(sLat, gz)

              coeff_T_TT_gsv     (lonIndex,latIndex,stepIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * (fottva(hu,1.0D0) * cmp + fotvt8(tt,hu) * cmp_TT) * ratioP1
              coeff_T_HU_gsv     (lonIndex,latIndex,stepIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * (folnqva(hu,tt,1.0d0) / hu * cmp + fotvt8(tt,hu) * cmp_HU) * ratioP1 

              coeff_T_P0_delP1  (lonIndex,latIndex,stepIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp

              coeff_T_P0_dP_delPT(lonIndex,latIndex,stepIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0_1 * ratioP1
              coeff_T_P0_dP_delP0(lonIndex,latIndex,stepIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0_2 * ratioP1
            else
              ! compute coefficients for GZ on momentum levels
              ratioP1 = log( P_M_ptr(lonIndex,latIndex,lev_T  ,stepIndex) / &
                             P_M_ptr(lonIndex,latIndex,lev_T-1,stepIndex) )
              hu = max(hu_ptr(lonIndex,latIndex,lev_T,stepIndex),MPC_MINIMUM_HU_R8)
              tt = tt_ptr(lonIndex,latIndex,lev_T,stepIndex)
              Pr = P_T_ptr(lonIndex,latIndex,lev_T,stepIndex)
              gz = gz_T_ptr(lonIndex,latIndex,lev_T,stepIndex)

              cmp = gpscompressibility(Pr,tt,hu)
              cmp_TT = gpscompressibility_TT(Pr,tt,hu)
              cmp_HU = gpscompressibility_HU(Pr,tt,hu)
              cmp_P0_1 = gpscompressibility_P0_1(Pr,tt,hu,1.0d0)
              cmp_P0_2 = gpscompressibility_P0_2(Pr,tt,hu)

              ! Gravity acceleration 
              lat_4 = statevector_trial%hco%lat2d_4(lonIndex,latIndex)
              lat_8 = real(lat_4,8)
              sLat = sin(lat_8)
              Rgh = phf_gravityalt(sLat, gz)

              coeff_M_TT_gsv     (lonIndex,latIndex,lev_T,stepIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * (fottva(hu,1.0D0) * cmp + fotvt8(tt,hu) * cmp_TT) * ratioP1
              coeff_M_HU_gsv     (lonIndex,latIndex,lev_T,stepIndex) = (MPC_RGAS_DRY_AIR_R8 / RgH) * (folnqva(hu,tt,1.0d0) / hu * cmp + fotvt8(tt,hu) * cmp_HU) * ratioP1 

              coeff_M_P0_delPM   (lonIndex,latIndex,lev_T,stepIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp 

              coeff_M_P0_dP_delPT(lonIndex,latIndex,lev_T,stepIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0_1 * ratioP1
              coeff_M_P0_dP_delP0(lonIndex,latIndex,lev_T,stepIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0_2 * ratioP1
            endif
          enddo
        enddo
      enddo
    enddo

  elseif (Vcode_anl == 5005) then

    do stepIndex = 1, numStep
      do lev_T = 1, (nlev_T-1)
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            ! compute coefficients for GZ on momentum levels
            ratioP1 = log( P_M_ptr(lonIndex,latIndex,lev_T+1,stepIndex) / &
                           P_M_ptr(lonIndex,latIndex,lev_T  ,stepIndex) )
            hu = max(hu_ptr(lonIndex,latIndex,lev_T,stepIndex),MPC_MINIMUM_HU_R8)
            tt = tt_ptr(lonIndex,latIndex,lev_T,stepIndex)
            Pr = P_T_ptr(lonIndex,latIndex,lev_T,stepIndex)
            gz = gz_T_ptr(lonIndex,latIndex,lev_T,stepIndex)

            cmp = gpscompressibility(Pr,tt,hu)
            cmp_TT = gpscompressibility_TT(Pr,tt,hu)
            cmp_HU = gpscompressibility_HU(Pr,tt,hu)
            cmp_P0_1 = gpscompressibility_P0_1(Pr,tt,hu,1.0d0)
            cmp_P0_2 = gpscompressibility_P0_2(Pr,tt,hu)

            ! Gravity acceleration 
            lat_4 = statevector_trial%hco%lat2d_4(lonIndex,latIndex)
            lat_8 = real(lat_4,8)
            sLat = sin(lat_8)
            Rgh = phf_gravityalt(sLat, gz)

            coeff_M_TT_gsv     (lonIndex,latIndex,lev_T,stepIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * (fottva(hu,1.0D0) * cmp + fotvt8(tt,hu) * cmp_TT) * ratioP1
            coeff_M_HU_gsv     (lonIndex,latIndex,lev_T,stepIndex) = (MPC_RGAS_DRY_AIR_R8 / RgH) * (folnqva(hu,tt,1.0d0) / hu * cmp + fotvt8(tt,hu) * cmp_HU) * ratioP1 

            coeff_M_P0_delPM   (lonIndex,latIndex,lev_T,stepIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp

            coeff_M_P0_dP_delPT(lonIndex,latIndex,lev_T,stepIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0_1 * ratioP1
            coeff_M_P0_dP_delP0(lonIndex,latIndex,lev_T,stepIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0_2 * ratioP1
          enddo
        enddo
      enddo
    enddo

  end if

!!$OMP END PARALLEL DO

  Write(*,*) "Exit subroutine calcAltitudeCoeff_gsv"

end subroutine calcAltitudeCoeff_gsv


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

! gpscompressibility_P0_1 has dpdp0 dependency
function gpscompressibility_P0_1(p,t,q,dpdp0)
  real(8), intent(in)  :: p,t,q,dpdp0
  real(8)              :: gpscompressibility_P0_1

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
  gpscompressibility_P0_1 = -d_pt * (a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2) + &
          2*pt*d_pt*(d+e*x2)
end function gpscompressibility_P0_1


! gpscompressibility_P0_2 has NO dpdp0 dependency
function gpscompressibility_P0_2(p,t,q)
  real(8), intent(in)  :: p,t,q
  real(8)              :: gpscompressibility_P0_2

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
  real(8)         :: d_x,d_tc,d_tc2,d_x2

  x  = p_wa * q / (1.D0+p_wb*q)
  d_x  = 0.0D0
  ! Estimate, from CIPM, Picard (2008)
  tc = t - MPC_K_C_DEGREE_OFFSET_R8
  d_tc = 0.0D0
  pt = p / t
  tc2= tc * tc
  d_tc2= 2 * tc * d_tc
  x2 = x * x
  d_x2 = 2 * x * d_x
  gpscompressibility_P0_2 = -pt * (a1*d_tc + a2*d_tc2 + b1*d_tc*x + (b0+b1*tc)*d_x + c1*d_tc*x2 + (c0+c1*tc)*d_x2) + &
          + pt*pt*e*d_x2
end function gpscompressibility_P0_2


end module tt2phi_mod
