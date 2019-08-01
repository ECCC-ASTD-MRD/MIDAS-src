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

module tt2phi_mod
  ! MODULE tt2phi (prefix='tt2phi' category='3. High-level transformations')
  !
  ! :Purpose: Subroutines for computing height from TT, HU and P0. Nonlinear,
  !           tangent-linear and adjoint versions of this transformation are
  !           included in separate subroutines.
  !
  use mpi_mod
  use mathPhysConstants_mod
  use physicsFunctions_mod
  use earthConstants_mod
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
  real(8), allocatable, save :: coeff_M_TT_gsv(:,:,:,:), coeff_M_HU_gsv(:,:,:,:)
  real(8), allocatable, save :: coeff_T_TT_gsv(:,:,:),   coeff_T_HU_gsv(:,:,:)
  real(8), allocatable, save :: coeff_M_P0_delPM(:,:,:,:), coeff_M_P0_dP_delPT(:,:,:,:), coeff_M_P0_dP_delP0(:,:,:,:)
  real(8), allocatable, save :: coeff_T_P0_delP1(:,:,:),   coeff_T_P0_dP_delPT(:,:,:),   coeff_T_P0_dP_delP0(:,:,:)

contains


subroutine tt2phi(statevector_trial,beSilent_opt)
  !
  ! :Purpose: Temperature to geopotential transformation on GEM4 staggered levels
  !           NOTE: we assume 
  !           1) nlev_T = nlev_M+1 
  !           2) alt_T(nlev_T) = alt_M(nlev_M), both at the surface
  !           3) a thermo level exists at the top, higher than the highest momentum level
  !           4) the placement of the thermo levels means that alt_T is the average of 2 nearest alt_M
  !           (according to Ron and Claude)
  !
  implicit none

  type(struct_gsv) :: statevector_trial
  logical, optional :: beSilent_opt

  integer :: lev_M,lev_T,nlev_M,nlev_T,status,Vcode,numStep,stepIndex,latIndex,lonIndex
  real(8) :: hu, tt, Pr, cmp, delThick, ratioP, ScaleFactorBottom, ScaleFactorTop
  real(8), allocatable :: tv(:), height_T(:), height_M(:) 
  real(8), pointer     :: height_T_ptr_r8(:,:,:,:),height_M_ptr_r8(:,:,:,:)
  real(8), pointer     :: hu_ptr_r8(:,:,:,:),tt_ptr_r8(:,:,:,:)
  real(8), pointer     :: P_T_ptr_r8(:,:,:,:),P_M_ptr_r8(:,:,:,:)
  real(8), pointer     :: P0_ptr_r8(:,:,:,:)
  real(8), pointer     :: HeightSfc_ptr_r8(:,:)
  real(4), pointer     :: height_T_ptr_r4(:,:,:,:),height_M_ptr_r4(:,:,:,:)
  real(4), pointer     :: hu_ptr_r4(:,:,:,:),tt_ptr_r4(:,:,:,:)
  real(4), pointer     :: P_T_ptr_r4(:,:,:,:),P_M_ptr_r4(:,:,:,:)
  real(4), pointer     :: P0_ptr_r4(:,:,:,:)
  real                 :: heightSfcOffset_T_r4, heightSfcOffset_M_r4
  real(4) :: lat_4
  real(8) :: lat_8, rMT
  real(8) :: h0, dh, Rgh, sLat, cLat
  type(struct_vco), pointer :: vco_ghr
  logical                   :: beSilent

  real(8) :: P_M, P_M1, P_Mm1, P_T, P0

  if ( present(beSilent_opt) ) then
    beSilent = beSilent_opt
  else
    beSilent = .true.
  end if

  call tmg_start(192,'tt2phi')

  if (.not.beSilent) write(*,*) 'tt2phi: START'

  vco_ghr => gsv_getVco(statevector_trial)
  Vcode = vco_ghr%vcode

  nlev_T = gsv_getNumLev(statevector_trial,'TH')
  nlev_M = gsv_getNumLev(statevector_trial,'MM')
  numStep = statevector_trial%numstep

  if (Vcode == 5002 .and. nlev_T /= nlev_M+1) call utl_abort('tt2phi: nlev_T is not equal to nlev_M+1!')
  if (Vcode == 5005 .and. nlev_T /= nlev_M)   call utl_abort('tt2phi: nlev_T is not equal to nlev_M!')

  if (Vcode == 5005) then
    status = vgd_get(statevector_trial%vco%vgrid,key='DHM - height of the diagnostic level (m)',value=heightSfcOffset_M_r4)
    status = vgd_get(statevector_trial%vco%vgrid,key='DHT - height of the diagnostic level (t)',value=heightSfcOffset_T_r4)
    if ( mpi_myid == 0 .and. .not.beSilent ) then
      write(*,*) 'tt2phi: height offset for near-sfc momentum level is: ', heightSfcOffset_M_r4, ' metres'
      write(*,*) 'tt2phi: height offset for near-sfc thermo level is:   ', heightSfcOffset_T_r4, ' metres'
      if ( .not.statevector_trial%addHeightSfcOffset ) then
        write(*,*) '----------------------------------------------------------------------------------'
        write(*,*) 'tt2phi: BUT HEIGHT OFFSET REMOVED FOR DIAGNOSTIC LEVELS FOR BACKWARD COMPATIBILITY'
        write(*,*) '----------------------------------------------------------------------------------'
      end if
    end if
  end if

  allocate(tv(nlev_T))
  allocate(height_T(nlev_T))
  allocate(height_M(nlev_M))

  if ( statevector_trial%dataKind == 4 ) then
    height_M_ptr_r4 => gsv_getField_r4(statevector_trial,'Z_M')
    height_T_ptr_r4 => gsv_getField_r4(statevector_trial,'Z_T')

    ! initialize the height pointer to zero
    height_M_ptr_r4(:,:,:,:) = 0.0
    height_T_ptr_r4(:,:,:,:) = 0.0
  else
    height_M_ptr_r8 => gsv_getField_r8(statevector_trial,'Z_M')
    height_T_ptr_r8 => gsv_getField_r8(statevector_trial,'Z_T')

    ! initialize the height pointer to zero
    height_M_ptr_r8(:,:,:,:) = 0.0d0
    height_T_ptr_r8(:,:,:,:) = 0.0d0
  end if

  if ( statevector_trial%dataKind == 4 ) then
    hu_ptr_r4 => gsv_getField_r4(statevector_trial,'HU')
    tt_ptr_r4 => gsv_getField_r4(statevector_trial,'TT')
    P_T_ptr_r4 => gsv_getField_r4(statevector_trial,'P_T')
    P_M_ptr_r4 => gsv_getField_r4(statevector_trial,'P_M')
    P0_ptr_r4 => gsv_getField_r4(statevector_trial,'P0')
  else
    hu_ptr_r8 => gsv_getField_r8(statevector_trial,'HU')
    tt_ptr_r8 => gsv_getField_r8(statevector_trial,'TT')
    P_T_ptr_r8 => gsv_getField_r8(statevector_trial,'P_T')
    P_M_ptr_r8 => gsv_getField_r8(statevector_trial,'P_M')
    P0_ptr_r8 => gsv_getField_r8(statevector_trial,'P0')
  end if
  HeightSfc_ptr_r8 => gsv_getHeightSfc(statevector_trial)

  ! compute virtual temperature on thermo levels (corrected of compressibility)
  do stepIndex = 1, numStep
    do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
      do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd

        height_T(:) = 0.0D0
        height_M(:) = 0.0D0

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
        end do

        rMT = HeightSfc_ptr_r8(lonIndex,latIndex)

        ! compute altitude on bottom momentum level
        if (Vcode == 5002) then
          height_M(nlev_M) = rMT
        elseif (Vcode == 5005) then
          height_M(nlev_M) = rMT + heightSfcOffset_M_r4
        end if 

        ! compute altitude on 2nd momentum level
        if (nlev_M > 1) then
          if ( statevector_trial%dataKind == 4 ) then
            P_M = real(P_M_ptr_r4(lonIndex,latIndex,nlev_M-1,stepIndex),8)
            P0  = real(P0_ptr_r4(lonIndex,latIndex,1,stepIndex),8)
          else
            P_M = P_M_ptr_r8(lonIndex,latIndex,nlev_M-1,stepIndex)
            P0  = P0_ptr_r8(lonIndex,latIndex,1,stepIndex)
          end if

          ratioP  = log( P_M / P0 )

          ! Gravity acceleration 
          h0  = rMT
          Rgh = phf_gravityalt(sLat,h0)
          dh  = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(nlev_T-1) * ratioP
          Rgh = phf_gravityalt(sLat, h0+0.5D0*dh)

          delThick   = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(nlev_T-1) * ratioP
          height_M(nlev_M-1) = rMT + delThick
        end if

        ! compute altitude on rest of momentum levels
        do lev_M = nlev_M-2, 1, -1
          if ( statevector_trial%dataKind == 4 ) then
            P_M = real(P_M_ptr_r4(lonIndex,latIndex,lev_M,stepIndex),8)
            P_M1 = real(P_M_ptr_r4(lonIndex,latIndex,lev_M+1,stepIndex),8)
          else
            P_M = P_M_ptr_r8(lonIndex,latIndex,lev_M,stepIndex)
            P_M1 = P_M_ptr_r8(lonIndex,latIndex,lev_M+1,stepIndex)
          end if

          ratioP  = log( P_M / P_M1 )

          if (Vcode == 5002) then
            lev_T = lev_M + 1
          else if (Vcode == 5005) then
            lev_T = lev_M
          end if

          ! Gravity acceleration 
          h0  = height_M(lev_M+1)
          Rgh = phf_gravityalt(sLat,h0)
          dh  = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(lev_T) * ratioP
          Rgh = phf_gravityalt(sLat, h0+0.5D0*dh)

          delThick   = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(lev_T) * ratioP
          height_M(lev_M) = height_M(lev_M+1) + delThick
        end do

        ! compute Altitude on thermo levels
        if (Vcode == 5002) then
          height_T(nlev_T) = height_M(nlev_M)

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
            height_T(lev_T) = ScaleFactorBottom * height_M(lev_M) + ScaleFactorTop * height_M(lev_M-1)
          end do

          ! compute altitude on top thermo level
          if ( statevector_trial%dataKind == 4 ) then
            P_T = real(P_T_ptr_r4(lonIndex,latIndex,1,stepIndex),8)
            P_M = real(P_M_ptr_r4(lonIndex,latIndex,1,stepIndex),8)
          else
            P_T = P_T_ptr_r8(lonIndex,latIndex,1,stepIndex)
            P_M = P_M_ptr_r8(lonIndex,latIndex,1,stepIndex)
          end if

          ratioP = log( P_T / P_M )

          ! Gravity acceleration 
          h0  = height_M(1)
          Rgh = phf_gravityalt(sLat, h0)
          dh  = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(1) * ratioP
          Rgh = phf_gravityalt(sLat, h0+0.5D0*dh)

          delThick   = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(1) * ratioP
          height_T(1) = height_M(1) + delThick

        else if (Vcode == 5005) then
          height_T(nlev_T) = rMT + heightSfcOffset_T_r4

          do lev_T = 1, nlev_T-2
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
            height_T(lev_T) = ScaleFactorBottom * height_M(lev_M) + ScaleFactorTop * height_M(lev_M-1)
          end do

          ! compute altitude on next to bottom thermo level
          if (nlev_T > 1) then
            if ( statevector_trial%dataKind == 4 ) then
              P_T = real(P_T_ptr_r4(lonIndex,latIndex,nlev_T-1,stepIndex),8)
              P0  = real(P0_ptr_r4(lonIndex,latIndex,1,stepIndex),8)
            else
              P_T = P_T_ptr_r8(lonIndex,latIndex,nlev_T-1,stepIndex)
              P0  = P0_ptr_r8(lonIndex,latIndex,1,stepIndex)
            end if

            ratioP = log( P_T / P0 )

            h0  = rMT
            Rgh = phf_gravityalt(sLat,h0)
            dh  = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(nlev_T-1) * ratioP
            Rgh = phf_gravityalt(sLat, h0+0.5D0*dh)

            delThick = (-MPC_RGAS_DRY_AIR_R8 / Rgh) * tv(nlev_T-1) * ratioP
            height_T(nlev_T-1) = rMT + delThick
          end if
        end if

        ! fill the height array
        if ( statevector_trial%dataKind == 4 ) then
          do lev_T = 1, nlev_T
            height_T_ptr_r4(lonIndex,latIndex,lev_T,stepIndex) = real(height_T(lev_T),4)
          end do
          do lev_M = 1, nlev_M
            height_M_ptr_r4(lonIndex,latIndex,lev_M,stepIndex) = real(height_M(lev_M),4)
          end do
        else
          height_T_ptr_r8(lonIndex,latIndex,1:nlev_T,stepIndex) = height_T(1:nlev_T)
          height_M_ptr_r8(lonIndex,latIndex,1:nlev_M,stepIndex) = height_M(1:nlev_M)
        end if

        ! remove the height offset for the diagnostic levels for backward compatibility only
        if ( .not. statevector_trial%addHeightSfcOffset ) then
          if ( statevector_trial%dataKind == 4 ) then
            height_T_ptr_r4(lonIndex,latIndex,nlev_T,stepIndex) = real(rMT,4)
            height_M_ptr_r4(lonIndex,latIndex,nlev_M,stepIndex) = real(rMT,4)
          else
            height_T_ptr_r8(lonIndex,latIndex,nlev_T,stepIndex) = rMT
            height_M_ptr_r8(lonIndex,latIndex,nlev_M,stepIndex) = rMT
          end if
        end if

      enddo
    enddo
  enddo

  deallocate(height_M)
  deallocate(height_T)
  deallocate(tv)

  if ( .not.beSilent ) then
    if ( statevector_trial%dataKind == 4 ) then
      write(*,*) 'tt2phi, Z_T='
      write(*,*) height_T_ptr_r4(statevector_trial%myLonBeg,statevector_trial%myLatBeg,:,1)
      write(*,*) 'tt2phi, Z_M='
      write(*,*) height_M_ptr_r4(statevector_trial%myLonBeg,statevector_trial%myLatBeg,:,1)
    else
      write(*,*) 'tt2phi, Z_T='
      write(*,*) height_T_ptr_r8(statevector_trial%myLonBeg,statevector_trial%myLatBeg,:,1)
      write(*,*) 'tt2phi, Z_M='
      write(*,*) height_M_ptr_r8(statevector_trial%myLonBeg,statevector_trial%myLatBeg,:,1)
    end if
    write(*,*) 'tt2phi: statevector_trial%addHeightSfcOffset=', statevector_trial%addHeightSfcOffset 
  end if

  if (.not.beSilent) write(*,*) 'tt2phi: END'

  call tmg_stop(192)

end subroutine tt2phi


subroutine tt2phi_tl(statevector,statevector_trial)
  !
  ! :Purpose: Temperature to geopotential transformation on gridstatevector
  !
  !
  implicit none

  type(struct_gsv) :: statevector,statevector_trial

  integer :: lev_M,lev_T,nlev_M,nlev_T,Vcode_anl,numStep,stepIndex,latIndex,lonIndex
  real(8) :: ScaleFactorBottom, ScaleFactorTop
  real(8), allocatable :: delThick(:,:,:,:)
  real(8), pointer     :: delHeight_M_ptr(:,:,:,:),delHeight_T_ptr(:,:,:,:),delTT(:,:,:,:),delHU(:,:,:,:),delP0(:,:,:,:)
  real(8), pointer     :: height_T_ptr(:,:,:,:),height_M_ptr(:,:,:,:)
  real(8), pointer     :: P_T(:,:,:,:), P_M(:,:,:,:)
  real(8), pointer     :: delP_T(:,:,:,:), delP_M(:,:,:,:)
  type(struct_vco), pointer :: vco_anl

  call tmg_start(193,'tt2phi_tl')

  write(*,*) 'tt2phi_tl: START'

  vco_anl => gsv_getVco(statevector_trial)
  Vcode_anl = vco_anl%vcode

  nlev_T = gsv_getNumLev(statevector_trial,'TH')
  nlev_M = gsv_getNumLev(statevector_trial,'MM')
  numStep = statevector_trial%numstep

  allocate(delThick(statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
                    statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
                    nlev_T,numStep))

  ! generate the height coefficients on the grid
  call calcHeightCoeff_gsv(statevector_trial)

  ! loop over all lat/lon/step
!!$OMP PARALLEL DO PRIVATE(columnIndex,delHeight_M_ptr,delHeight_T_ptr,delThick,delTT,delHU,delP0,lev_M,lev_T)

  height_M_ptr => gsv_getField_r8(statevector_trial,'Z_M')
  height_T_ptr => gsv_getField_r8(statevector_trial,'Z_T')
  P_T => gsv_getField_r8(statevector_trial,'P_T')
  P_M => gsv_getField_r8(statevector_trial,'P_M')

  delHeight_M_ptr => gsv_getField_r8(statevector,'Z_M')
  delHeight_T_ptr => gsv_getField_r8(statevector,'Z_T')
  delTT => gsv_getField_r8(statevector,'TT')
  delHU => gsv_getField_r8(statevector,'HU')
  delP0 => gsv_getField_r8(statevector,'P0')
  delP_T => gsv_getField_r8(statevector,'P_T')
  delP_M => gsv_getField_r8(statevector,'P_M')

  ! ensure increment at sfc is zero (fixed height level)
  delHeight_M_ptr(:,:,nlev_M,:) = 0.0d0
  delHeight_T_ptr(:,:,nlev_T,:) = 0.0d0

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

    ! compute height increment on momentum levels above the surface
    do stepIndex = 1, numStep
      do lev_M = (nlev_M-1), 1, -1
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            lev_T = lev_M + 1 ! thermo level just below momentum level being computed
            delHeight_M_ptr(lonIndex,latIndex,lev_M,stepIndex) = delHeight_M_ptr(lonIndex,latIndex,lev_M+1,stepIndex) + delThick(lonIndex,latIndex,lev_T,stepIndex)
          enddo
        enddo
      enddo
    enddo

    ! compute height increment on thermo levels using weighted average of height increment of momentum levels
    do stepIndex = 1, numStep
      do lev_T = 1, (nlev_T-1)
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            if ( lev_T == 1) then
              ! compute height increment for top thermo level (from top momentum level)
              delHeight_T_ptr(lonIndex,latIndex,1,stepIndex) = delHeight_M_ptr(lonIndex,latIndex,1,stepIndex) +  &
                                 coeff_T_TT_gsv(lonIndex,latIndex,stepIndex) * delTT(lonIndex,latIndex,1,stepIndex) + &
                                 coeff_T_HU_gsv(lonIndex,latIndex,stepIndex) * delHU(lonIndex,latIndex,1,stepIndex) + &

                                 coeff_T_P0_delP1(lonIndex,latIndex,stepIndex) * &
                                 ( delP_M(lonIndex,latIndex,1,stepIndex) / P_M(lonIndex,latIndex,1,stepIndex) - &
                                   delP_T(lonIndex,latIndex,1,stepIndex) / P_T(lonIndex,latIndex,1,stepIndex) ) + &

                                 coeff_T_P0_dP_delPT(lonIndex,latIndex,stepIndex) * delP_T(lonIndex,latIndex,1,stepIndex) + &
                                 coeff_T_P0_dP_delP0(lonIndex,latIndex,stepIndex) * delP0(lonIndex,latIndex,1,stepIndex)

            else
              lev_M = lev_T ! momentum level just below thermo level being computed
              ScaleFactorBottom = (height_T_ptr(lonIndex,latIndex,lev_T,stepIndex) - height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex)) / &
                                  (height_M_ptr(lonIndex,latIndex,lev_M,stepIndex) - height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex))
              ScaleFactorTop    = 1 - ScaleFactorBottom
              delHeight_T_ptr(lonIndex,latIndex,lev_T,stepIndex) = ScaleFactorBottom * delHeight_M_ptr(lonIndex,latIndex,lev_M  ,stepIndex) + &
                                                              ScaleFactorTop * delHeight_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex)
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

    ! compute height increment on momentum levels above the surface
    do stepIndex = 1, numStep
      do lev_M = (nlev_M-1), 1, -1
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            lev_T = lev_M ! thermo level just below momentum level being computed
            delHeight_M_ptr(lonIndex,latIndex,lev_M,stepIndex) = delHeight_M_ptr(lonIndex,latIndex,lev_M+1,stepIndex) + delThick(lonIndex,latIndex,lev_T,stepIndex)
          enddo
        enddo
      enddo
    enddo

    ! compute height increment on thermo levels using weighted average of height increment of momentum levels
    do stepIndex = 1, numStep
      do lev_T = 1, (nlev_T-1)
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            lev_M = lev_T + 1 ! momentum level just below thermo level being computed
            ScaleFactorBottom = (height_T_ptr(lonIndex,latIndex,lev_T,stepIndex) - height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex)) / &
                                (height_M_ptr(lonIndex,latIndex,lev_M,stepIndex) - height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex))
            ScaleFactorTop    = 1 - ScaleFactorBottom
            delHeight_T_ptr(lonIndex,latIndex,lev_T,stepIndex) = ScaleFactorBottom * delHeight_M_ptr(lonIndex,latIndex,lev_M  ,stepIndex) + &
                                                            ScaleFactorTop * delHeight_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex)
          enddo
        enddo
      enddo
    enddo

  endif

!!$OMP END PARALLEL DO

  deallocate(delThick)

  write(*,*) 'tt2phi_tl: END'

  call tmg_stop(193)

end subroutine tt2phi_tl


subroutine tt2phi_ad(statevector,statevector_trial)
  !
  !:Purpose: Adjoint of temperature to geopotential transformation on
  !          gridstatevector
  !
  !
  implicit none

  type(struct_gsv) :: statevector,statevector_trial

  integer :: lev_M,lev_T,nlev_M,nlev_T,Vcode_anl,numStep,stepIndex,latIndex,lonIndex
  real(8) :: ScaleFactorBottom, ScaleFactorTop
  real(8), allocatable :: delThick(:,:,:,:)
  real(8), pointer     :: delHeight_M_ptr(:,:,:,:),delHeight_T_ptr(:,:,:,:),delTT(:,:,:,:),delHU(:,:,:,:),delP0(:,:,:,:)
  real(8), pointer     :: height_M_ptr(:,:,:,:),height_T_ptr(:,:,:,:)
  real(8), allocatable :: delHeight_M(:,:,:,:),delHeight_T(:,:,:,:)
  real(8), pointer     :: P_M(:,:,:,:),P_T(:,:,:,:)
  real(8), pointer     :: delP_M(:,:,:,:),delP_T(:,:,:,:)
  type(struct_vco), pointer :: vco_anl


  call tmg_start(194,'tt2phi_ad')

  write(*,*) 'tt2phi_ad: START'

  vco_anl => gsv_getVco(statevector_trial)
  Vcode_anl = vco_anl%vcode

  nlev_T = gsv_getNumLev(statevector_trial,'TH')
  nlev_M = gsv_getNumLev(statevector_trial,'MM')
  numStep = statevector_trial%numstep

  allocate(delHeight_M(statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
                   statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
                   nlev_M,numStep))
  allocate(delHeight_T(statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
                   statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
                   nlev_T,numStep))
  allocate(delThick(statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
                   statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
                   0:nlev_T,numStep))

  ! generate the height coefficients on the grid
  call calcHeightCoeff_gsv(statevector_trial)

  ! loop over all lat/lon/step
!!$OMP PARALLEL DO PRIVATE(columnIndex,delHeight_M,delHeight_T,delThick,delTT,delHU,delP0,lev_M,lev_T)

  height_M_ptr => gsv_getField_r8(statevector_trial,'Z_M')
  height_T_ptr => gsv_getField_r8(statevector_trial,'Z_T')
  P_T => gsv_getField_r8(statevector_trial,'P_T')
  P_M => gsv_getField_r8(statevector_trial,'P_M')

  delHeight_M_ptr => gsv_getField_r8(statevector,'Z_M')
  delHeight_T_ptr => gsv_getField_r8(statevector,'Z_T')

  delTT => gsv_getField_r8(statevector,'TT')
  delHU => gsv_getField_r8(statevector,'HU')
  delP0 => gsv_getField_r8(statevector,'P0')
  delP_T => gsv_getField_r8(statevector,'P_T')
  delP_M => gsv_getField_r8(statevector,'P_M')

  delHeight_M(:,:,:,:) = delHeight_M_ptr(:,:,:,:)
  delHeight_T(:,:,:,:) = delHeight_T_ptr(:,:,:,:)

  if(Vcode_anl .eq. 5002) then

    ! adjoint of compute height increment on thermo levels by simple averaging
    do stepIndex = 1, numStep
      do lev_T = 1, (nlev_T-1)
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            lev_M = lev_T ! momentum level just below thermo level being computed

            ! adjoint of compute height increment on top thermo level (from top momentum level)
            if (lev_T == 1) then
              delHeight_M(lonIndex,latIndex,1,stepIndex)  = delHeight_M(lonIndex,latIndex,1,stepIndex) + &
                                                        delHeight_T(lonIndex,latIndex,1,stepIndex)

              delTT(lonIndex,latIndex,1,stepIndex) = delTT(lonIndex,latIndex,1,stepIndex) + &
                                                     coeff_T_TT_gsv   (lonIndex,latIndex,stepIndex) * delHeight_T(lonIndex,latIndex,1,stepIndex)
              delHU(lonIndex,latIndex,1,stepIndex) = delHU(lonIndex,latIndex,1,stepIndex) + &
                                                     coeff_T_HU_gsv   (lonIndex,latIndex,stepIndex) * delHeight_T(lonIndex,latIndex,1,stepIndex)

              delP_M(lonIndex,latIndex,1,stepIndex) = delP_M(lonIndex,latIndex,1,stepIndex) + &
                                                      coeff_T_P0_delP1(lonIndex,latIndex,stepIndex) / P_M(lonIndex,latIndex,1,stepIndex) * &
                                                      delHeight_T(lonIndex,latIndex,1,stepIndex)

              delP_T(lonIndex,latIndex,1,stepIndex) = delP_T(lonIndex,latIndex,1,stepIndex) - &
                                                      coeff_T_P0_delP1(lonIndex,latIndex,stepIndex) / P_T(lonIndex,latIndex,1,stepIndex) * &
                                                      delHeight_T(lonIndex,latIndex,1,stepIndex)

              delP_T(lonIndex,latIndex,1,stepIndex) = delP_T(lonIndex,latIndex,1,stepIndex) + &
                                                     coeff_T_P0_dP_delPT(lonIndex,latIndex,stepIndex) * delHeight_T(lonIndex,latIndex,1,stepIndex)

              delP0(lonIndex,latIndex,1,stepIndex) = delP0(lonIndex,latIndex,1,stepIndex) + &
                                                     coeff_T_P0_dp_delP0(lonIndex,latIndex,stepIndex) * delHeight_T(lonIndex,latIndex,1,stepIndex)
            else
              ScaleFactorBottom = (height_T_ptr(lonIndex,latIndex,lev_T,stepIndex) - height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex)) / &
                                  (height_M_ptr(lonIndex,latIndex,lev_M,stepIndex) - height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex))
              ScaleFactorTop    = 1 - ScaleFactorBottom
              delHeight_M(lonIndex,latIndex,lev_M-1,stepIndex) = delHeight_M(lonIndex,latIndex,lev_M-1,stepIndex) + &
                                            ScaleFactorTop * delHeight_T(lonIndex,latIndex,lev_T  ,stepIndex)
              delHeight_M(lonIndex,latIndex,lev_M  ,stepIndex) = delHeight_M(lonIndex,latIndex,lev_M  ,stepIndex) + &
                                         ScaleFactorBottom * delHeight_T(lonIndex,latIndex,lev_T  ,stepIndex)
            end if
          enddo
        enddo
      enddo
    enddo

    ! adjoint of compute height increment on momentum levels above the surface
    delThick(:,:,0:1,:) = 0.0d0
    do stepIndex = 1, numStep
      do lev_M = 1, (nlev_M-1)
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            lev_T = lev_M + 1 ! thermo level just below momentum level being computed
            delThick(lonIndex,latIndex,lev_T,stepIndex) = delThick(lonIndex,latIndex,lev_T-1,stepIndex) + &
                                                           delHeight_M(lonIndex,latIndex,lev_M  ,stepIndex)
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

    ! adjoint of compute height increment on thermo levels by simple averaging
    do stepIndex = 1, numStep
      do lev_T = 1, (nlev_T-1)
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            lev_M = lev_T+1 ! momentum level just below thermo level being computed
            ScaleFactorBottom = (height_T_ptr(lonIndex,latIndex,lev_T,stepIndex) - height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex)) / &
                                (height_M_ptr(lonIndex,latIndex,lev_M,stepIndex) - height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex))
            ScaleFactorTop    = 1 - ScaleFactorBottom
            delHeight_M(lonIndex,latIndex,lev_M-1,stepIndex) = delHeight_M(lonIndex,latIndex,lev_M-1,stepIndex) + &
                                          ScaleFactorTop * delHeight_T(lonIndex,latIndex,lev_T  ,stepIndex)
            delHeight_M(lonIndex,latIndex,lev_M,stepIndex)   = delHeight_M(lonIndex,latIndex,lev_M  ,stepIndex) + &
                                       ScaleFactorBottom * delHeight_T(lonIndex,latIndex,lev_T  ,stepIndex)
          enddo
        enddo
      enddo
    enddo

    ! adjoint of compute height increment on momentum levels
    delThick(:,:,0,:) = 0.0d0
    do stepIndex = 1, numStep
      do lev_M = 1, (nlev_M-1)
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            lev_T = lev_M ! thermo level just below momentum level being computed
            delThick(lonIndex,latIndex,lev_T,stepIndex) = delThick(lonIndex,latIndex,lev_T-1,stepIndex) + &
                                                          delHeight_M (lonIndex,latIndex,lev_M  ,stepIndex)
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

  deallocate(delThick)
  deallocate(delHeight_M)
  deallocate(delHeight_T)

  write(*,*) 'tt2phi_ad: END'

  call tmg_stop(194)

end subroutine tt2phi_ad


subroutine calcHeightCoeff_gsv(statevector_trial)
  !
  ! :Purpose: Calculating the coefficients of height for tt2phi_tl/tt2phi_ad
  !
  implicit none

  type(struct_gsv) :: statevector_trial

  integer :: lev_T,nlev_M,nlev_T,numStep,stepIndex,latIndex,lonIndex,Vcode_anl
  real(8) :: hu,tt,Pr,height_T,cmp,cmp_TT,cmp_HU,cmp_P0_1,cmp_P0_2,ratioP1
  real(4) :: lat_4
  real(8) :: Rgh, sLat, lat_8
  real(8), pointer :: hu_ptr(:,:,:,:),tt_ptr(:,:,:,:)
  real(8), pointer :: P_T_ptr(:,:,:,:),P_M_ptr(:,:,:,:)
  real(8), pointer :: height_T_ptr(:,:,:,:)
  type(struct_vco), pointer :: vco_anl

  logical, save :: firstTimeHeightCoeff_gsv = .true.

  if ( .not. firstTimeHeightCoeff_gsv ) return

  Write(*,*) "calcHeightCoeff_gsv: START"

  ! initialize and save coefficients for increased efficiency (assumes no relinearization)
  firstTimeHeightCoeff_gsv = .false.

  vco_anl => gsv_getVco(statevector_trial)
  Vcode_anl = vco_anl%vcode 

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

!!$OMP PARALLEL DO PRIVATE(columnIndex,delThick,delTT,delHU,delP0,lev_M,lev_T)
  hu_ptr => gsv_getField_r8(statevector_trial,'HU')
  tt_ptr => gsv_getField_r8(statevector_trial,'TT')
  P_T_ptr => gsv_getField_r8(statevector_trial,'P_T')
  P_M_ptr => gsv_getField_r8(statevector_trial,'P_M')
  height_T_ptr => gsv_getField_r8(statevector_trial,'Z_T')

  if (Vcode_anl == 5002) then

    do stepIndex = 1, numStep
      do lev_T = 1, (nlev_T-1)
        do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
          do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
            if ( lev_T == 1 ) then 
              ! compute height coefficients on only the top thermo level
              ratioP1 = log( P_M_ptr(lonIndex,latIndex,1,stepIndex) / &
                             P_T_ptr(lonIndex,latIndex,1,stepIndex) ) 
              hu = max(hu_ptr(lonIndex,latIndex,1,stepIndex),MPC_MINIMUM_HU_R8)
              tt = tt_ptr(lonIndex,latIndex,1,stepIndex)
              Pr = P_T_ptr(lonIndex,latIndex,1,stepIndex)
              height_T = height_T_ptr(lonIndex,latIndex,1,stepIndex)

              cmp = gpscompressibility(Pr,tt,hu)
              cmp_TT = gpscompressibility_TT(Pr,tt,hu)
              cmp_HU = gpscompressibility_HU(Pr,tt,hu)
              cmp_P0_1 = gpscompressibility_P0_1(Pr,tt,hu,1.0d0)
              cmp_P0_2 = gpscompressibility_P0_2(Pr,tt,hu)

              ! Gravity acceleration 
              lat_4 = statevector_trial%hco%lat2d_4(lonIndex,latIndex)
              lat_8 = real(lat_4,8)
              sLat = sin(lat_8)
              Rgh = phf_gravityalt(sLat, height_T)

              coeff_T_TT_gsv     (lonIndex,latIndex,stepIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * (fottva(hu,1.0D0) * cmp + fotvt8(tt,hu) * cmp_TT) * ratioP1
              coeff_T_HU_gsv     (lonIndex,latIndex,stepIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * (folnqva(hu,tt,1.0d0) / hu * cmp + fotvt8(tt,hu) * cmp_HU) * ratioP1 

              coeff_T_P0_delP1  (lonIndex,latIndex,stepIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp

              coeff_T_P0_dP_delPT(lonIndex,latIndex,stepIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0_1 * ratioP1
              coeff_T_P0_dP_delP0(lonIndex,latIndex,stepIndex) = (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0_2 * ratioP1
            else
              ! compute height coefficients on momentum levels
              ratioP1 = log( P_M_ptr(lonIndex,latIndex,lev_T  ,stepIndex) / &
                             P_M_ptr(lonIndex,latIndex,lev_T-1,stepIndex) )
              hu = max(hu_ptr(lonIndex,latIndex,lev_T,stepIndex),MPC_MINIMUM_HU_R8)
              tt = tt_ptr(lonIndex,latIndex,lev_T,stepIndex)
              Pr = P_T_ptr(lonIndex,latIndex,lev_T,stepIndex)
              height_T = height_T_ptr(lonIndex,latIndex,lev_T,stepIndex)

              cmp = gpscompressibility(Pr,tt,hu)
              cmp_TT = gpscompressibility_TT(Pr,tt,hu)
              cmp_HU = gpscompressibility_HU(Pr,tt,hu)
              cmp_P0_1 = gpscompressibility_P0_1(Pr,tt,hu,1.0d0)
              cmp_P0_2 = gpscompressibility_P0_2(Pr,tt,hu)

              ! Gravity acceleration 
              lat_4 = statevector_trial%hco%lat2d_4(lonIndex,latIndex)
              lat_8 = real(lat_4,8)
              sLat = sin(lat_8)
              Rgh = phf_gravityalt(sLat, height_T)

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
            ! compute height coefficients on momentum levels
            ratioP1 = log( P_M_ptr(lonIndex,latIndex,lev_T+1,stepIndex) / &
                           P_M_ptr(lonIndex,latIndex,lev_T  ,stepIndex) )
            hu = max(hu_ptr(lonIndex,latIndex,lev_T,stepIndex),MPC_MINIMUM_HU_R8)
            tt = tt_ptr(lonIndex,latIndex,lev_T,stepIndex)
            Pr = P_T_ptr(lonIndex,latIndex,lev_T,stepIndex)
            height_T = height_T_ptr(lonIndex,latIndex,lev_T,stepIndex)

            cmp = gpscompressibility(Pr,tt,hu)
            cmp_TT = gpscompressibility_TT(Pr,tt,hu)
            cmp_HU = gpscompressibility_HU(Pr,tt,hu)
            cmp_P0_1 = gpscompressibility_P0_1(Pr,tt,hu,1.0d0)
            cmp_P0_2 = gpscompressibility_P0_2(Pr,tt,hu)

            ! Gravity acceleration 
            lat_4 = statevector_trial%hco%lat2d_4(lonIndex,latIndex)
            lat_8 = real(lat_4,8)
            sLat = sin(lat_8)
            Rgh = phf_gravityalt(sLat, height_T)

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

  Write(*,*) "calcHeightCoeff_gsv: END"

end subroutine calcHeightCoeff_gsv


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
