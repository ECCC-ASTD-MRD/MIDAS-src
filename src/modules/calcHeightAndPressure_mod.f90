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

module calcHeightAndPressure_mod
  ! MODULE czp_calcHeightAndPressure (prefix='czp' category='3. High-level transformations')
  !
  ! :Purpose: Subroutines for computing height or pressure from TT, HU, P0 or 
  !           orography. 
  !           Nonlinear, tangent-linear and adjoint versions of this 
  !           transformation are included in separate subroutines.
  !
  use codePrecision_mod
  use mpi_mod
  use mathPhysConstants_mod
  use physicsFunctions_mod
  use earthConstants_mod
  use verticalCoord_mod
  use gridstatevector_mod
  use columnData_mod
  use utilities_mod
  use varnamelist_mod
  implicit none
  save
  private

  ! public procedures
  public :: czp_tt2phi, czp_tt2phi_tl, czp_tt2phi_ad, &
            czp_calcGridPressure_nl, czp_calcGridPressure_tl, &
            czp_calcGridPressure_ad, & 
            czp_calcColumnZandP_tl, czp_calcColumnZandP_ad, &
            czp_calcColumnPressure_nl, czp_calcColumnPressure_tl, &
            czp_calcColumnPressure_ad 

  ! (mad001) @TODO : interface for calc(Grid|Pressure)_*

  interface czp_tt2phi_tl
    module procedure czp_tt2phi_gsv_tl
    module procedure czp_tt2phi_col_tl
  end interface czp_tt2phi_tl
  interface czp_tt2phi_ad
    module procedure czp_tt2phi_gsv_ad
    module procedure czp_tt2phi_col_ad
  end interface czp_tt2phi_ad
  
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
  real(8), allocatable :: coeff_M_TT_gsv(:,:,:,:), coeff_M_HU_gsv(:,:,:,:)
  real(8), allocatable :: coeff_T_TT_gsv(:,:,:),   coeff_T_HU_gsv(:,:,:)
  real(8), allocatable :: coeff_M_P0_delPM_gsv(:,:,:,:), coeff_M_P0_dP_delPT_gsv(:,:,:,:)
  real(8), allocatable :: coeff_M_P0_dP_delP0_gsv(:,:,:,:)
  real(8), allocatable :: coeff_T_P0_delP1_gsv(:,:,:),   coeff_T_P0_dP_delPT_gsv(:,:,:)
  real(8), allocatable :: coeff_T_P0_dP_delP0_gsv(:,:,:)

  real(8), allocatable :: coeff_M_TT_col(:,:), coeff_M_HU_col(:,:)
  real(8), allocatable :: coeff_T_TT_col(:),   coeff_T_HU_col(:)
  real(8), allocatable :: coeff_M_P0_delPM_col(:,:), coeff_M_P0_dP_delPT_col(:,:)
  real(8), allocatable :: coeff_M_P0_dP_delP0_col(:,:)
  real(8), allocatable :: coeff_T_P0_delP1_col(:),   coeff_T_P0_dP_delPT_col(:)
  real(8), allocatable :: coeff_T_P0_dP_delP0_col(:)

contains

  !---------------------------------------------------------
  ! czp_tt2phi
  !---------------------------------------------------------
  subroutine czp_tt2phi(statevector_trial,beSilent_opt)
    !
    ! :Purpose: Temperature to geopotential transformation on GEM4 staggered 
    !           levels
    !           NOTE: we assume 
    !           1) nlev_T = nlev_M+1 
    !           2) alt_T(nlev_T) = alt_M(nlev_M), both at the surface
    !           3) a thermo level exists at the top, higher than the highest 
    !              momentum level
    !           4) the placement of the thermo levels means that alt_T is the 
    !              average of 2 nearest alt_M
    !              (according to Ron and Claude)
    !
    implicit none
  
    type(struct_gsv) :: statevector_trial
    logical, optional :: beSilent_opt
  
    integer ::  lev_M,lev_T,nlev_M,nlev_T,status,Vcode,numStep,stepIndex, &
                latIndex,lonIndex
    real(8) ::  hu, tt, Pr, cmp, delThick, ratioP, ScaleFactorBottom, &
                ScaleFactorTop
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
    real(4)              :: heightSfcOffset_T_r4, heightSfcOffset_M_r4
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
  
    call tmg_start(192,'czp_tt2phi')
  
    if (.not.beSilent) write(*,*) 'czp_tt2phi: START'
  
    vco_ghr => gsv_getVco(statevector_trial)
    Vcode = vco_ghr%vcode
  
    nlev_T = gsv_getNumLev(statevector_trial,'TH')
    nlev_M = gsv_getNumLev(statevector_trial,'MM')
    numStep = statevector_trial%numstep
  
    if (Vcode == 5002 .and. nlev_T /= nlev_M+1) then
      call utl_abort('czp_tt2phi: nlev_T is not equal to nlev_M+1!')
    end if
    if (Vcode == 5005 .and. nlev_T /= nlev_M) then
      call utl_abort('czp_tt2phi: nlev_T is not equal to nlev_M!')
    end if

    if (Vcode == 5005) then
      status = vgd_get(&
          statevector_trial%vco%vgrid,key='DHM - height of the diagnostic level (m)', &
          value=heightSfcOffset_M_r4)
      status = vgd_get(&
            statevector_trial%vco%vgrid,key='DHT - height of the diagnostic level (t)', &
            value=heightSfcOffset_T_r4)
      if ( mpi_myid == 0 .and. .not.beSilent ) then
        write(*,*) 'czp_tt2phi: height offset for near-sfc momentum level is: ', &
              heightSfcOffset_M_r4, ' metres'
        write(*,*) 'czp_tt2phi: height offset for near-sfc thermo level is:   ', &
              heightSfcOffset_T_r4, ' metres'
        if ( .not.statevector_trial%addHeightSfcOffset ) then
          write(*,*) '----------------------------------------------------------------------------------'
          write(*,*) 'czp_tt2phi: BUT HEIGHT OFFSET REMOVED FOR DIAGNOSTIC LEVELS FOR BACKWARD COMPATIBILITY'
          write(*,*) '----------------------------------------------------------------------------------'
        end if
      end if
    end if
  
    allocate(tv(nlev_T))
    allocate(height_T(nlev_T))
    allocate(height_M(nlev_M))
  
    if ( statevector_trial%dataKind == 4 ) then
      call gsv_getField(statevector_trial,height_M_ptr_r4,'Z_M')
      call gsv_getField(statevector_trial,height_T_ptr_r4,'Z_T')
  
      ! initialize the height pointer to zero
      height_M_ptr_r4(:,:,:,:) = 0.0
      height_T_ptr_r4(:,:,:,:) = 0.0
    else
      call gsv_getField(statevector_trial,height_M_ptr_r8,'Z_M')
      call gsv_getField(statevector_trial,height_T_ptr_r8,'Z_T')
  
      ! initialize the height pointer to zero
      height_M_ptr_r8(:,:,:,:) = 0.0d0
      height_T_ptr_r8(:,:,:,:) = 0.0d0
    end if
  
    if ( statevector_trial%dataKind == 4 ) then
      call gsv_getField(statevector_trial,hu_ptr_r4,'HU')
      call gsv_getField(statevector_trial,tt_ptr_r4,'TT')
      call gsv_getField(statevector_trial,P_T_ptr_r4,'P_T')
      call gsv_getField(statevector_trial,P_M_ptr_r4,'P_M')
      call gsv_getField(statevector_trial,P0_ptr_r4,'P0')
    else
      call gsv_getField(statevector_trial,hu_ptr_r8,'HU')
      call gsv_getField(statevector_trial,tt_ptr_r8,'TT')
      call gsv_getField(statevector_trial,P_T_ptr_r8,'P_T')
      call gsv_getField(statevector_trial,P_M_ptr_r8,'P_M')
      call gsv_getField(statevector_trial,P0_ptr_r8,'P0')
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
              height_T(lev_T) = ScaleFactorBottom * height_M(lev_M) &
                + ScaleFactorTop * height_M(lev_M-1)
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
              height_T(lev_T) = ScaleFactorBottom * height_M(lev_M) &
                + ScaleFactorTop * height_M(lev_M-1)
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
              height_T_ptr_r4(lonIndex,latIndex,lev_T,stepIndex) &
                  = real(height_T(lev_T),4)
            end do
            do lev_M = 1, nlev_M
              height_M_ptr_r4(lonIndex,latIndex,lev_M,stepIndex) &
                  = real(height_M(lev_M),4)
            end do
          else
            height_T_ptr_r8(lonIndex,latIndex,1:nlev_T,stepIndex) &
                = height_T(1:nlev_T)
            height_M_ptr_r8(lonIndex,latIndex,1:nlev_M,stepIndex) &
                = height_M(1:nlev_M)
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
        write(*,*) 'czp_tt2phi, Z_T='
        write(*,*) height_T_ptr_r4(statevector_trial%myLonBeg, &
            statevector_trial%myLatBeg,:,1)
        write(*,*) 'czp_tt2phi, Z_M='
        write(*,*) height_M_ptr_r4(statevector_trial%myLonBeg, &
            statevector_trial%myLatBeg,:,1)
      else
        write(*,*) 'czp_tt2phi, Z_T='
        write(*,*) height_T_ptr_r8(statevector_trial%myLonBeg, &
            statevector_trial%myLatBeg,:,1)
        write(*,*) 'czp_tt2phi, Z_M='
        write(*,*) height_M_ptr_r8(statevector_trial%myLonBeg, &
            statevector_trial%myLatBeg,:,1)
      end if
      write(*,*) 'czp_tt2phi: statevector_trial%addHeightSfcOffset=', &
          statevector_trial%addHeightSfcOffset 
    end if
  
    if (.not.beSilent) write(*,*) 'czp_tt2phi: END'
  
    call tmg_stop(192)
  
  end subroutine czp_tt2phi
  
  !---------------------------------------------------------
  ! czp_tt2phi_gsv_tl
  !---------------------------------------------------------
  subroutine czp_tt2phi_gsv_tl(statevector,statevector_trial)
    !
    ! :Purpose: Temperature to geopotential transformation on gridstatevector
    !
    !
    implicit none
  
    type(struct_gsv) :: statevector,statevector_trial
  
    integer ::  lev_M,lev_T,nlev_M,nlev_T,Vcode_anl,numStep,stepIndex, &
                latIndex,lonIndex
    real(8) :: ScaleFactorBottom, ScaleFactorTop
    real(8), allocatable :: delThick(:,:,:,:)
    real(8), pointer     :: height_T_ptr(:,:,:,:),height_M_ptr(:,:,:,:)
    real(8), pointer     :: P_T(:,:,:,:), P_M(:,:,:,:)
    real(pre_incrReal), pointer ::  delHeight_M_ptr_r48(:,:,:,:), &
                                    delHeight_T_ptr_r48(:,:,:,:)
    real(pre_incrReal), pointer ::  delTT_r48(:,:,:,:), delHU_r48(:,:,:,:), &
                                    delP0_r48(:,:,:,:)
    real(pre_incrReal), pointer ::  delP_T_r48(:,:,:,:), delP_M_r48(:,:,:,:)
    type(struct_vco),   pointer ::  vco_anl
  
    call tmg_start(193,'czp_tt2phi_tl')
  
    write(*,*) 'czp_tt2phi_gsv_tl: START'
  
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
  
    call gsv_getField(statevector_trial,height_M_ptr,'Z_M')
    call gsv_getField(statevector_trial,height_T_ptr,'Z_T')
    call gsv_getField(statevector_trial,P_T,'P_T')
    call gsv_getField(statevector_trial,P_M,'P_M')
  
    call gsv_getField(statevector,delHeight_M_ptr_r48,'Z_M')
    call gsv_getField(statevector,delHeight_T_ptr_r48,'Z_T')
    call gsv_getField(statevector,delTT_r48,'TT')
    call gsv_getField(statevector,delHU_r48,'HU')
    call gsv_getField(statevector,delP0_r48,'P0')
    call gsv_getField(statevector,delP_T_r48,'P_T')
    call gsv_getField(statevector,delP_M_r48,'P_M')
    ! ensure increment at sfc is zero (fixed height level)
    delHeight_M_ptr_r48(:,:,nlev_M,:) = 0.0d0
    delHeight_T_ptr_r48(:,:,nlev_T,:) = 0.0d0
  
    if (Vcode_anl == 5002) then
  
      ! compute increment to thickness for each layer between the two momentum levels
      do stepIndex = 1, numStep
        do lev_T = 2, (nlev_T-1)
          do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
            do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
              delThick(lonIndex,latIndex,lev_T,stepIndex) = &
                  coeff_M_TT_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                  delTT_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                  coeff_M_HU_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                  delHU_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                  coeff_M_P0_delPM_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                  ( delP_M_r48(lonIndex,latIndex,lev_T  ,stepIndex) / &
                    P_M(lonIndex,latIndex,lev_T  ,stepIndex) - &
                    delP_M_r48(lonIndex,latIndex,lev_T-1,stepIndex) / &
                    P_M(lonIndex,latIndex,lev_T-1,stepIndex) ) + &
                  coeff_M_P0_dP_delPT_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                  delP_T_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                  coeff_M_P0_dP_delP0_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                  delP0_r48(lonIndex,latIndex,1,stepIndex)
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
              delHeight_M_ptr_r48(lonIndex,latIndex,lev_M,stepIndex) =  &
                   delHeight_M_ptr_r48(lonIndex,latIndex,lev_M+1,stepIndex) &
                   + delThick(lonIndex,latIndex,lev_T,stepIndex)
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
                delHeight_T_ptr_r48(lonIndex,latIndex,1,stepIndex) = &
                delHeight_M_ptr_r48(lonIndex,latIndex,1,stepIndex) + &
                coeff_T_TT_gsv(lonIndex,latIndex,stepIndex) * &
                delTT_r48(lonIndex,latIndex,1,stepIndex) + &
                coeff_T_HU_gsv(lonIndex,latIndex,stepIndex) * &
                delHU_r48(lonIndex,latIndex,1,stepIndex) + &
                coeff_T_P0_delP1_gsv(lonIndex,latIndex,stepIndex) * &
                ( delP_M_r48(lonIndex,latIndex,1,stepIndex) / &
                  P_M(lonIndex,latIndex,1,stepIndex) - &
                  delP_T_r48(lonIndex,latIndex,1,stepIndex) / &
                  P_T(lonIndex,latIndex,1,stepIndex) ) + &
                coeff_T_P0_dP_delPT_gsv(lonIndex,latIndex,stepIndex) * &
                delP_T_r48(lonIndex,latIndex,1,stepIndex) + &
                coeff_T_P0_dP_delP0_gsv(lonIndex,latIndex,stepIndex) * &
                delP0_r48(lonIndex,latIndex,1,stepIndex)
              else
                lev_M = lev_T ! momentum level just below thermo level being computed
                ScaleFactorBottom = &
                    (height_T_ptr(lonIndex,latIndex,lev_T,stepIndex) - &
                      height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex)) / &
                    (height_M_ptr(lonIndex,latIndex,lev_M,stepIndex) - &
                      height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex))
                ScaleFactorTop    = 1 - ScaleFactorBottom
                delHeight_T_ptr_r48(lonIndex,latIndex,lev_T,stepIndex) = &
                    ScaleFactorBottom * &
                    delHeight_M_ptr_r48(lonIndex,latIndex,lev_M  ,stepIndex) + &
                    ScaleFactorTop * &
                    delHeight_M_ptr_r48(lonIndex,latIndex,lev_M-1,stepIndex)
              endif
            enddo
          enddo
        enddo
      enddo
  
    elseif(Vcode_anl == 5005) then
  
      ! compute increment to thickness for each layer between the two momentum levels
      do stepIndex = 1, numStep
        do lev_T = 1, (nlev_T-1)
          do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
            do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
              delThick(lonIndex,latIndex,lev_T,stepIndex) = &
                  coeff_M_TT_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                  delTT_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                  coeff_M_HU_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                  delHU_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                  coeff_M_P0_delPM_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                  ( delP_M_r48(lonIndex,latIndex,lev_T+1,stepIndex) / &
                    P_M(lonIndex,latIndex,lev_T+1,stepIndex) - &
                    delP_M_r48(lonIndex,latIndex,lev_T  ,stepIndex) / &
                    P_M(lonIndex,latIndex,lev_T  ,stepIndex) ) + &
                  coeff_M_P0_dP_delPT_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                  delP_T_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                  coeff_M_P0_dP_delP0_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                  delP0_r48(lonIndex,latIndex,1,stepIndex)
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
              delHeight_M_ptr_r48(lonIndex,latIndex,lev_M,stepIndex) = &
              delHeight_M_ptr_r48(lonIndex,latIndex,lev_M+1,stepIndex) + &
              delThick(lonIndex,latIndex,lev_T,stepIndex)
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
              ScaleFactorBottom = &
                  (height_T_ptr(lonIndex,latIndex,lev_T,stepIndex) - &
                    height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex)) / &
                  (height_M_ptr(lonIndex,latIndex,lev_M,stepIndex) - &
                    height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex))
              ScaleFactorTop    = 1 - ScaleFactorBottom
              delHeight_T_ptr_r48(lonIndex,latIndex,lev_T,stepIndex) = &
                  ScaleFactorBottom * &
                  delHeight_M_ptr_r48(lonIndex,latIndex,lev_M  ,stepIndex) + &
                  ScaleFactorTop * &
                  delHeight_M_ptr_r48(lonIndex,latIndex,lev_M-1,stepIndex)
            enddo
          enddo
        enddo
      enddo
  
    endif
  
    deallocate(delThick)
  
    write(*,*) 'czp_tt2phi_gsv_tl: END'
  
    call tmg_stop(193)
  
  end subroutine czp_tt2phi_gsv_tl
  
  !---------------------------------------------------------
  ! czp_tt2phi_col_tl
  !---------------------------------------------------------
  subroutine czp_tt2phi_col_tl(columnInc,columnTrlOnAnlInc)
    !
    ! :Purpose: Temperature to geopotential transformation on gridstatevector
    !
    !
    implicit none
  
    type(struct_columnData) :: columnInc,columnTrlOnAnlInc
  
    integer :: lev_M,lev_T,nlev_M,nlev_T,Vcode_anl,colIndex,numColumns
    real(8) :: ScaleFactorBottom, ScaleFactorTop
    real(8), allocatable :: delThick(:,:)
    real(8), pointer     :: height_T_ptr(:,:),height_M_ptr(:,:)
    real(8), pointer     :: P_T(:,:), P_M(:,:)
    real(8), pointer  :: delHeight_M_ptr(:,:),delHeight_T_ptr(:,:)
    real(8), pointer  :: delTT(:,:),delHU(:,:),delP0(:,:)
    real(8), pointer  :: delP_T(:,:), delP_M(:,:)
    type(struct_vco), pointer :: vco_anl
  
    call tmg_start(193,'czp_tt2phi_tl')
  
    write(*,*) 'czp_tt2phi_col_tl: START'
  
    vco_anl => col_getVco(columnTrlOnAnlInc)
    Vcode_anl = vco_anl%vcode
  
    nlev_T = col_getNumLev(columnTrlOnAnlInc,'TH')
    nlev_M = col_getNumLev(columnTrlOnAnlInc,'MM')
  
    numColumns = col_getNumCol(columnInc)
    
    allocate(delThick(nlev_T,numColumns))
    delThick(:,:) = 0.0d0
  
    ! generate the height coefficients on the grid
    call calcHeightCoeff_col(columnTrlOnAnlInc)
  
    ! loop over all lat/lon/step
  
    height_M_ptr => col_getAllColumns(columnTrlOnAnlInc,'Z_M')
    height_T_ptr => col_getAllColumns(columnTrlOnAnlInc,'Z_T')
    P_M          => col_getAllColumns(columnTrlOnAnlInc,'P_M')
    P_T          => col_getAllColumns(columnTrlOnAnlInc,'P_T')
  
    delHeight_M_ptr => col_getAllColumns(columnInc,'Z_M')
    delHeight_T_ptr => col_getAllColumns(columnInc,'Z_T')
    delTT           => col_getAllColumns(columnInc,'TT')
    delHU           => col_getAllColumns(columnInc,'HU')
    delP0           => col_getAllColumns(columnInc,'P0')
    delP_M          => col_getAllColumns(columnInc,'P_M')
    delP_T          => col_getAllColumns(columnInc,'P_T')
  
    ! ensure increment at sfc is zero (fixed height level)
    delHeight_M_ptr(nlev_M,:) = 0.0d0
    delHeight_T_ptr(nlev_T,:) = 0.0d0
  
    if (Vcode_anl == 5002) then
  
      ! compute increment to thickness for each layer between the two momentum levels
      do colIndex = 1, numColumns
        do lev_T = 2, (nlev_T-1)
          delThick(lev_T,colIndex) =  &
                coeff_M_TT_col(lev_T,colIndex) * delTT(lev_T,colIndex) + &
                coeff_M_HU_col(lev_T,colIndex) * delHU(lev_T,colIndex) + &
                coeff_M_P0_delPM_col(lev_T,colIndex) * &
                ( delP_M(lev_T  ,colIndex) / P_M(lev_T  ,colIndex) - &
                  delP_M(lev_T-1,colIndex) / P_M(lev_T-1,colIndex) ) + &
                coeff_M_P0_dP_delPT_col(lev_T,colIndex) * &
                delP_T(lev_T,colIndex) + &
                coeff_M_P0_dP_delP0_col(lev_T,colIndex) * delP0(1,colIndex)
        enddo
      enddo
  
      ! compute height increment on momentum levels above the surface
      do colIndex = 1, numColumns
        do lev_M = (nlev_M-1), 1, -1
          lev_T = lev_M + 1 ! thermo level just below momentum level being computed
          delHeight_M_ptr(lev_M,colIndex) =  &
               delHeight_M_ptr(lev_M+1,colIndex) + delThick(lev_T,colIndex)
        enddo
      enddo
  
      ! compute height increment on thermo levels using weighted average of height increment of momentum levels
      do colIndex = 1, numColumns
        do lev_T = 1, (nlev_T-1)
          if ( lev_T == 1) then
            ! compute height increment for top thermo level (from top momentum level)
            delHeight_T_ptr(1,colIndex) = delHeight_M_ptr(1,colIndex) +  &
                 coeff_T_TT_col(colIndex) * delTT(1,colIndex) + &
                 coeff_T_HU_col(colIndex) * delHU(1,colIndex) + &
                 coeff_T_P0_delP1_col(colIndex) * &
                 ( delP_M(1,colIndex) / P_M(1,colIndex) - &
                   delP_T(1,colIndex) / P_T(1,colIndex) ) + &
                 coeff_T_P0_dP_delPT_col(colIndex) * delP_T(1,colIndex) + &
                 coeff_T_P0_dP_delP0_col(colIndex) * delP0(1,colIndex)
          else
            lev_M = lev_T ! momentum level just below thermo level being computed
            ScaleFactorBottom = (height_T_ptr(lev_T,colIndex) - &
                height_M_ptr(lev_M-1,colIndex)) / &
                (height_M_ptr(lev_M,colIndex) - height_M_ptr(lev_M-1,colIndex))
            ScaleFactorTop    = 1 - ScaleFactorBottom
            delHeight_T_ptr(lev_T,colIndex) =  &
                 ScaleFactorBottom * delHeight_M_ptr(lev_M  ,colIndex) + &
                 ScaleFactorTop * delHeight_M_ptr(lev_M-1,colIndex)
          endif
        enddo
      enddo
  
    elseif(Vcode_anl == 5005) then
  
      ! compute increment to thickness for each layer between the two momentum levels
      do colIndex = 1, numColumns
        do lev_T = 1, (nlev_T-1)
          delThick(lev_T,colIndex) =  &
               coeff_M_TT_col(lev_T,colIndex) * delTT(lev_T,colIndex) + &
               coeff_M_HU_col(lev_T,colIndex) * delHU(lev_T,colIndex) + &
               coeff_M_P0_delPM_col(lev_T,colIndex) * &
               ( delP_M(lev_T+1,colIndex) / P_M(lev_T+1,colIndex) - &
                 delP_M(lev_T  ,colIndex) / P_M(lev_T  ,colIndex) ) + &
               coeff_M_P0_dP_delPT_col(lev_T,colIndex) * &
               delP_T(lev_T,colIndex) + &
               coeff_M_P0_dP_delP0_col(lev_T,colIndex) * delP0(1,colIndex)
        enddo
      enddo
  
      ! compute height increment on momentum levels above the surface
      do colIndex = 1, numColumns
        do lev_M = (nlev_M-1), 1, -1
          lev_T = lev_M ! thermo level just below momentum level being computed
          delHeight_M_ptr(lev_M,colIndex) =  &
               delHeight_M_ptr(lev_M+1,colIndex) + delThick(lev_T,colIndex)
        enddo
      enddo
  
      ! compute height increment on thermo levels using weighted average of height increment of momentum levels
      do colIndex = 1, numColumns
        do lev_T = 1, (nlev_T-1)
          lev_M = lev_T + 1 ! momentum level just below thermo level being computed
          ScaleFactorBottom =  &
               (height_T_ptr(lev_T,colIndex) - height_M_ptr(lev_M-1,colIndex)) / &
               (height_M_ptr(lev_M,colIndex) - height_M_ptr(lev_M-1,colIndex))
          ScaleFactorTop    = 1 - ScaleFactorBottom
          delHeight_T_ptr(lev_T,colIndex) =  &
               ScaleFactorBottom * delHeight_M_ptr(lev_M  ,colIndex) + &
               ScaleFactorTop * delHeight_M_ptr(lev_M-1,colIndex)
        enddo
      enddo
  
    endif
  
    deallocate(delThick)
  
    write(*,*) 'czp_tt2phi_col_tl: END'
  
    call tmg_stop(193)
  
  end subroutine czp_tt2phi_col_tl
  
  !---------------------------------------------------------
  ! czp_tt2phi_gsv_ad
  !---------------------------------------------------------
  subroutine czp_tt2phi_gsv_ad(statevector,statevector_trial)
    !
    !:Purpose: Adjoint of temperature to geopotential transformation on
    !          gridstatevector
    !
    !
    implicit none
  
    type(struct_gsv) :: statevector,statevector_trial
  
    integer ::  lev_M,lev_T,nlev_M,nlev_T,Vcode_anl,numStep,stepIndex, &
                latIndex,lonIndex
    real(8) :: ScaleFactorBottom, ScaleFactorTop
    real(8), allocatable :: delThick(:,:,:,:)
    real(8), pointer     :: height_M_ptr(:,:,:,:),height_T_ptr(:,:,:,:)
    real(8), allocatable :: delHeight_M(:,:,:,:),delHeight_T(:,:,:,:)
    real(8), pointer     :: P_M(:,:,:,:),P_T(:,:,:,:)
    real(pre_incrReal), pointer :: delHeight_M_ptr_r48(:,:,:,:),delHeight_T_ptr_r48(:,:,:,:)
    real(pre_incrReal), pointer :: delTT_r48(:,:,:,:),delHU_r48(:,:,:,:),delP0_r48(:,:,:,:)
    real(pre_incrReal), pointer :: delP_M_r48(:,:,:,:),delP_T_r48(:,:,:,:)
    type(struct_vco),   pointer :: vco_anl
  
    call tmg_start(194,'czp_tt2phi_ad')
  
    write(*,*) 'czp_tt2phi_gsv_ad: START'
  
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
  
    call gsv_getField(statevector_trial,height_M_ptr,'Z_M')
    call gsv_getField(statevector_trial,height_T_ptr,'Z_T')
    call gsv_getField(statevector_trial,P_T,'P_T')
    call gsv_getField(statevector_trial,P_M,'P_M')
  
    call gsv_getField(statevector,delHeight_M_ptr_r48,'Z_M')
    call gsv_getField(statevector,delHeight_T_ptr_r48,'Z_T')
    call gsv_getField(statevector,delTT_r48,'TT')
    call gsv_getField(statevector,delHU_r48,'HU')
    call gsv_getField(statevector,delP0_r48,'P0')
    call gsv_getField(statevector,delP_T_r48,'P_T')
    call gsv_getField(statevector,delP_M_r48,'P_M')
    delHeight_M(:,:,:,:) = delHeight_M_ptr_r48(:,:,:,:)
    delHeight_T(:,:,:,:) = delHeight_T_ptr_r48(:,:,:,:)
  
    if(Vcode_anl == 5002) then
  
      ! adjoint of compute height increment on thermo levels by simple averaging
      do stepIndex = 1, numStep
        do lev_T = 1, (nlev_T-1)
          do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
            do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
              lev_M = lev_T ! momentum level just below thermo level being computed
  
              ! adjoint of compute height increment on top thermo level 
              ! (from top momentum level)
              if (lev_T == 1) then
                delHeight_M(lonIndex,latIndex,1,stepIndex)  =  &
                    delHeight_M(lonIndex,latIndex,1,stepIndex) + &
                    delHeight_T(lonIndex,latIndex,1,stepIndex)
                  
                delTT_r48(lonIndex,latIndex,1,stepIndex) =  &
                    delTT_r48(lonIndex,latIndex,1,stepIndex) + &
                    coeff_T_TT_gsv(lonIndex,latIndex,stepIndex) * &
                    delHeight_T(lonIndex,latIndex,1,stepIndex)
  
                delHU_r48(lonIndex,latIndex,1,stepIndex) =  &
                    delHU_r48(lonIndex,latIndex,1,stepIndex) + &
                    coeff_T_HU_gsv   (lonIndex,latIndex,stepIndex) * &
                    delHeight_T(lonIndex,latIndex,1,stepIndex)
  
                delP_M_r48(lonIndex,latIndex,1,stepIndex) =  &
                    delP_M_r48(lonIndex,latIndex,1,stepIndex) + &
                    coeff_T_P0_delP1_gsv(lonIndex,latIndex,stepIndex) / &
                    P_M(lonIndex,latIndex,1,stepIndex) * &
                    delHeight_T(lonIndex,latIndex,1,stepIndex)
  
                delP_T_r48(lonIndex,latIndex,1,stepIndex) =  &
                    delP_T_r48(lonIndex,latIndex,1,stepIndex) - &
                    coeff_T_P0_delP1_gsv(lonIndex,latIndex,stepIndex) / &
                    P_T(lonIndex,latIndex,1,stepIndex) * &
                    delHeight_T(lonIndex,latIndex,1,stepIndex)
  
                delP_T_r48(lonIndex,latIndex,1,stepIndex) =  &
                    delP_T_r48(lonIndex,latIndex,1,stepIndex) + &
                    coeff_T_P0_dP_delPT_gsv(lonIndex,latIndex,stepIndex) * &
                    delHeight_T(lonIndex,latIndex,1,stepIndex)
  
                delP0_r48(lonIndex,latIndex,1,stepIndex) =  &
                    delP0_r48(lonIndex,latIndex,1,stepIndex) + &
                    coeff_T_P0_dp_delP0_gsv(lonIndex,latIndex,stepIndex) * &
                    delHeight_T(lonIndex,latIndex,1,stepIndex)
              else
                ScaleFactorBottom =  &
                    (height_T_ptr(lonIndex,latIndex,lev_T,stepIndex) - &
                      height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex)) / &
                    (height_M_ptr(lonIndex,latIndex,lev_M,stepIndex) - &
                      height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex))
                ScaleFactorTop    = 1 - ScaleFactorBottom
  
                delHeight_M(lonIndex,latIndex,lev_M-1,stepIndex) =  &
                    delHeight_M(lonIndex,latIndex,lev_M-1,stepIndex) + &
                    ScaleFactorTop * &
                    delHeight_T(lonIndex,latIndex,lev_T,stepIndex)
  
                delHeight_M(lonIndex,latIndex,lev_M,stepIndex) =  &
                    delHeight_M(lonIndex,latIndex,lev_M  ,stepIndex) + &
                    ScaleFactorBottom * &
                    delHeight_T(lonIndex,latIndex,lev_T,stepIndex)
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
              delThick(lonIndex,latIndex,lev_T,stepIndex) =  &
                   delThick(lonIndex,latIndex,lev_T-1,stepIndex) + &
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
              delTT_r48(lonIndex,latIndex,lev_T,stepIndex) =  &
                  delTT_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                  coeff_M_TT_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                  delThick(lonIndex,latIndex,lev_T,stepIndex)
  
              delHU_r48(lonIndex,latIndex,lev_T,stepIndex) =  &
                  delHU_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                  coeff_M_HU_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                  delThick(lonIndex,latIndex,lev_T,stepIndex)
  
              delP_M_r48(lonIndex,latIndex,lev_T,stepIndex)=  &
                  delP_M_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                  coeff_M_P0_delPM_gsv(lonIndex,latIndex,lev_T,stepIndex) / &
                  P_M(lonIndex,latIndex,lev_T,stepIndex) * &
                  delThick(lonIndex,latIndex,lev_T,stepIndex)
  
              delP_M_r48(lonIndex,latIndex,lev_T-1,stepIndex) =  &
                  delP_M_r48(lonIndex,latIndex,lev_T-1,stepIndex) - &
                  coeff_M_P0_delPM_gsv(lonIndex,latIndex,lev_T,stepIndex) / &
                  P_M(lonIndex,latIndex,lev_T-1,stepIndex) * &
                  delThick(lonIndex,latIndex,lev_T,stepIndex)
  
              delP_T_r48(lonIndex,latIndex,lev_T,stepIndex) =  &
                  delP_T_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                  coeff_M_P0_dP_delPT_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                  delThick(lonIndex,latIndex,lev_T,stepIndex)
  
              delP0_r48(lonIndex,latIndex,1,stepIndex) =  &
                  delP0_r48(lonIndex,latIndex,1,stepIndex) + &
                  coeff_M_P0_dP_delP0_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                  delThick(lonIndex,latIndex,lev_T,stepIndex)
            enddo
          enddo
        enddo
      enddo
  
    elseif(Vcode_anl == 5005) then
  
      ! adjoint of compute height increment on thermo levels by simple averaging
      do stepIndex = 1, numStep
        do lev_T = 1, (nlev_T-1)
          do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
            do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
              lev_M = lev_T+1 ! momentum level just below thermo level being computed
              ScaleFactorBottom = &
                  (height_T_ptr(lonIndex,latIndex,lev_T,stepIndex) - &
                    height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex)) / &
                  (height_M_ptr(lonIndex,latIndex,lev_M,stepIndex) -&
                    height_M_ptr(lonIndex,latIndex,lev_M-1,stepIndex))
              ScaleFactorTop = 1 - ScaleFactorBottom
              delHeight_M(lonIndex,latIndex,lev_M-1,stepIndex) = &
                  delHeight_M(lonIndex,latIndex,lev_M-1,stepIndex) + &
                  ScaleFactorTop * &
                  delHeight_T(lonIndex,latIndex,lev_T,stepIndex)
              delHeight_M(lonIndex,latIndex,lev_M,stepIndex) = &
                  delHeight_M(lonIndex,latIndex,lev_M  ,stepIndex) + &
                  ScaleFactorBottom * &
                  delHeight_T(lonIndex,latIndex,lev_T,stepIndex)
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
              delThick(lonIndex,latIndex,lev_T,stepIndex) = &
                  delThick(lonIndex,latIndex,lev_T-1,stepIndex) + &
                  delHeight_M (lonIndex,latIndex,lev_M,stepIndex)
            end do
          end do
        end do
      end do
  
      do stepIndex = 1, numStep
        do lev_T = 1, nlev_T-1
          do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
            do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
              delTT_r48(lonIndex,latIndex,lev_T,stepIndex) =  &
                  delTT_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                  coeff_M_TT_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                  delThick(lonIndex,latIndex,lev_T,stepIndex)
  
              delHU_r48(lonIndex,latIndex,lev_T,stepIndex) =  &
                  delHU_r48(lonIndex,latIndex,lev_T,stepIndex) + &
                  coeff_M_HU_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                  delThick(lonIndex,latIndex,lev_T,stepIndex)
                
              delP_M_r48(lonIndex,latIndex,lev_T+1,stepIndex) =  &
                  delP_M_r48(lonIndex,latIndex,lev_T+1,stepIndex) + &
                  coeff_M_P0_delPM_gsv(lonIndex,latIndex,lev_T,stepIndex) / &
                  P_M(lonIndex,latIndex,lev_T+1,stepIndex) * &
                  delThick(lonIndex,latIndex,lev_T,stepIndex)
  
              delP_M_r48(lonIndex,latIndex,lev_T  ,stepIndex) =  &
                  delP_M_r48(lonIndex,latIndex,lev_T  ,stepIndex) - &
                  coeff_M_P0_delPM_gsv(lonIndex,latIndex,lev_T,stepIndex) / &
                  P_M(lonIndex,latIndex,lev_T  ,stepIndex) * &
                  delThick(lonIndex,latIndex,lev_T,stepIndex)
  
              delP_T_r48(lonIndex,latIndex,lev_T  ,stepIndex) =  &
                  delP_T_r48(lonIndex,latIndex,lev_T  ,stepIndex) + &
                  coeff_M_P0_dP_delPT_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                  delThick(lonIndex,latIndex,lev_T,stepIndex)
  
              delP0_r48(lonIndex,latIndex,1,stepIndex)     =  &
                  delP0_r48(lonIndex,latIndex,1,stepIndex) + &
                  coeff_M_P0_dP_delP0_gsv(lonIndex,latIndex,lev_T,stepIndex) * &
                  delThick(lonIndex,latIndex,lev_T,stepIndex)
            enddo
          enddo
        enddo
      enddo
  
    endif
  
    deallocate(delThick)
    deallocate(delHeight_M)
    deallocate(delHeight_T)
  
    write(*,*) 'czp_tt2phi_gsv_ad: END'
  
    call tmg_stop(194)
  
  end subroutine czp_tt2phi_gsv_ad
  
  !---------------------------------------------------------
  ! czp_tt2phi_col_ad
  !---------------------------------------------------------
  subroutine czp_tt2phi_col_ad(columnInc,columnTrlOnAnlInc)
    !
    !:Purpose: Adjoint of temperature to geopotential transformation on
    !          columnData
    !
    !
    implicit none
  
    type(struct_columnData) :: columnInc,columnTrlOnAnlInc
  
    integer :: lev_M,lev_T,nlev_M,nlev_T,Vcode_anl,numColumns,colIndex
    real(8) :: ScaleFactorBottom, ScaleFactorTop
    real(8), allocatable :: delThick(:,:)
    real(8), pointer     :: height_M_ptr(:,:),height_T_ptr(:,:)
    real(8), allocatable :: delHeight_M(:,:),delHeight_T(:,:)
    real(8), pointer     :: P_M(:,:),P_T(:,:)
    real(8), pointer     :: delHeight_M_ptr(:,:),delHeight_T_ptr(:,:)
    real(8), pointer     :: delTT(:,:),delHU(:,:),delP0(:,:)
    real(8), pointer     :: delP_M(:,:),delP_T(:,:)
    type(struct_vco), pointer :: vco_anl
  
    call tmg_start(194,'czp_tt2phi_ad')
  
    write(*,*) 'czp_tt2phi_col_ad: START'
  
    vco_anl => col_getVco(columnTrlOnAnlInc)
    Vcode_anl = vco_anl%vcode
  
    nlev_T = col_getNumLev(columnTrlOnAnlInc,'TH')
    nlev_M = col_getNumLev(columnTrlOnAnlInc,'MM')
    numColumns = col_getNumCol(columnTrlOnAnlInc)
  
    allocate(delHeight_M(nlev_M,numColumns))
    allocate(delHeight_T(nlev_T,numColumns))
    allocate(delThick(0:nlev_T,numColumns))
  
    ! generate the height coefficients on the grid
    call calcHeightCoeff_col(columnTrlOnAnlInc)
  
    height_M_ptr => col_getAllColumns(columnTrlOnAnlInc,'Z_M')
    height_T_ptr => col_getAllColumns(columnTrlOnAnlInc,'Z_T')
    P_M          => col_getAllColumns(columnTrlOnAnlInc,'P_M')
    P_T          => col_getAllColumns(columnTrlOnAnlInc,'P_T')
  
    delHeight_M_ptr => col_getAllColumns(columnInc,'Z_M')
    delHeight_T_ptr => col_getAllColumns(columnInc,'Z_T')
    delTT           => col_getAllColumns(columnInc,'TT')
    delHU           => col_getAllColumns(columnInc,'HU')
    delP0           => col_getAllColumns(columnInc,'P0')
    delP_M          => col_getAllColumns(columnInc,'P_M')
    delP_T          => col_getAllColumns(columnInc,'P_T')
  
    delHeight_M(:,:) = delHeight_M_ptr(:,:)
    delHeight_T(:,:) = delHeight_T_ptr(:,:)
  
    if(Vcode_anl == 5002) then
  
      ! adjoint of compute height increment on thermo levels by simple averaging
      do colIndex = 1, numColumns
        do lev_T = 1, (nlev_T-1)
          lev_M = lev_T ! momentum level just below thermo level being computed
  
          ! adjoint of compute height increment on top thermo level (from top momentum level)
          if (lev_T == 1) then
            delHeight_M(1,colIndex)  =  &
                 delHeight_M(1,colIndex) + &
                 delHeight_T(1,colIndex)
                  
            delTT(1,colIndex) =  &
                 delTT(1,colIndex) + &
                 coeff_T_TT_col(colIndex) * delHeight_T(1,colIndex)
  
            delHU(1,colIndex) =  &
                 delHU(1,colIndex) + &
                 coeff_T_HU_col   (colIndex) * delHeight_T(1,colIndex)
  
            delP_M(1,colIndex) =  &
                 delP_M(1,colIndex) + &
                 coeff_T_P0_delP1_col(colIndex) / P_M(1,colIndex) * &
                 delHeight_T(1,colIndex)
  
            delP_T(1,colIndex) =  &
                 delP_T(1,colIndex) - &
                 coeff_T_P0_delP1_col(colIndex) / P_T(1,colIndex) * &
                 delHeight_T(1,colIndex)
  
            delP_T(1,colIndex) =  &
                 delP_T(1,colIndex) + &
                 coeff_T_P0_dP_delPT_col(colIndex) * delHeight_T(1,colIndex)
  
            delP0(1,colIndex) =  &
                 delP0(1,colIndex) + &
                 coeff_T_P0_dp_delP0_col(colIndex) * delHeight_T(1,colIndex)
          else
            ScaleFactorBottom =  &
                (height_T_ptr(lev_T,colIndex) - & 
                  height_M_ptr(lev_M-1,colIndex)) / &
                (height_M_ptr(lev_M,colIndex) - height_M_ptr(lev_M-1,colIndex))
            ScaleFactorTop    = 1 - ScaleFactorBottom
  
            delHeight_M(lev_M-1,colIndex) =  &
                 delHeight_M(lev_M-1,colIndex) + &
                 ScaleFactorTop * delHeight_T(lev_T,colIndex)
  
            delHeight_M(lev_M,colIndex) =  &
                 delHeight_M(lev_M  ,colIndex) + &
                 ScaleFactorBottom * delHeight_T(lev_T,colIndex)
          end if
        enddo
      enddo
  
      ! adjoint of compute height increment on momentum levels above the surface
      delThick(0:1,:) = 0.0d0
      do colIndex = 1, numColumns
        do lev_M = 1, (nlev_M-1)
          lev_T = lev_M + 1 ! thermo level just below momentum level being computed
          delThick(lev_T,colIndex) =  &
               delThick(lev_T-1,colIndex) + &
               delHeight_M(lev_M  ,colIndex)
        end do
      end do
  
      ! adjoint of compute increment to thickness for each layer between the two momentum levels
      do colIndex = 1, numColumns
        do lev_T = 2, nlev_T-1
          delTT(lev_T,colIndex) =  &
              delTT            (lev_T,colIndex) + &
              coeff_M_TT_col   (lev_T,colIndex) * delThick(lev_T,colIndex)
  
          delHU(lev_T,colIndex) =  &
              delHU            (lev_T,colIndex) + &
              coeff_M_HU_col   (lev_T,colIndex) * delThick(lev_T,colIndex)
  
          delP_M(lev_T,colIndex)=  &
              delP_M(lev_T,colIndex) + &
              coeff_M_P0_delPM_col(lev_T,colIndex) / P_M(lev_T,colIndex) * &
              delThick(lev_T,colIndex)
  
          delP_M(lev_T-1,colIndex) =  &
              delP_M(lev_T-1,colIndex) - &
              coeff_M_P0_delPM_col(lev_T,colIndex) / P_M(lev_T-1,colIndex) * &
              delThick(lev_T,colIndex)
  
          delP_T(lev_T,colIndex) =  &
              delP_T(lev_T,colIndex) + &
              coeff_M_P0_dP_delPT_col(lev_T,colIndex) * &
              delThick(lev_T,colIndex)
  
          delP0(1,colIndex) =  &
              delP0(1,colIndex) + &
              coeff_M_P0_dP_delP0_col(lev_T,colIndex) * &
              delThick(lev_T,colIndex)
        enddo
      enddo
  
    elseif(Vcode_anl == 5005) then
  
      ! adjoint of compute height increment on thermo levels by simple averaging
      do colIndex = 1, numColumns
        do lev_T = 1, (nlev_T-1)
          lev_M = lev_T+1 ! momentum level just below thermo level being computed
          ScaleFactorBottom = (height_T_ptr(lev_T,colIndex) - &
                height_M_ptr(lev_M-1,colIndex)) / &
              (height_M_ptr(lev_M,colIndex) - height_M_ptr(lev_M-1,colIndex))
          ScaleFactorTop    = 1 - ScaleFactorBottom
          delHeight_M(lev_M-1,colIndex) = delHeight_M(lev_M-1,colIndex) + &
              ScaleFactorTop * delHeight_T(lev_T  ,colIndex)
          delHeight_M(lev_M,colIndex)   = delHeight_M(lev_M  ,colIndex) + &
              ScaleFactorBottom * delHeight_T(lev_T  ,colIndex)
        enddo
      enddo
  
      ! adjoint of compute height increment on momentum levels
      delThick(0,:) = 0.0d0
      do colIndex = 1, numColumns
        do lev_M = 1, (nlev_M-1)
          lev_T = lev_M ! thermo level just below momentum level being computed
          delThick(lev_T,colIndex) = delThick(lev_T-1,colIndex) + &
                                     delHeight_M (lev_M  ,colIndex)
        end do
      end do
  
      do colIndex = 1, numColumns
        do lev_T = 1, nlev_T-1
          delTT(lev_T,colIndex) =  &
              delTT(lev_T,colIndex) + &
              coeff_M_TT_col(lev_T,colIndex) * delThick(lev_T,colIndex)
  
          delHU(lev_T,colIndex) =  &
              delHU(lev_T,colIndex) + &
              coeff_M_HU_col(lev_T,colIndex) * delThick(lev_T,colIndex)
                
          delP_M(lev_T+1,colIndex) =  &
              delP_M(lev_T+1,colIndex) + &
              coeff_M_P0_delPM_col(lev_T,colIndex) / P_M(lev_T+1,colIndex) * &
              delThick(lev_T,colIndex)
  
          delP_M(lev_T  ,colIndex) =  &
              delP_M(lev_T  ,colIndex) - &
              coeff_M_P0_delPM_col(lev_T,colIndex) / P_M(lev_T  ,colIndex) * &
              delThick(lev_T,colIndex)
  
          delP_T(lev_T  ,colIndex) =  &
              delP_T(lev_T  ,colIndex) + &
              coeff_M_P0_dP_delPT_col(lev_T,colIndex) * delThick(lev_T,colIndex)
  
          delP0(1,colIndex)     =  &
              delP0(1,colIndex) + &
              coeff_M_P0_dP_delP0_col(lev_T,colIndex) * delThick(lev_T,colIndex)
        enddo
      enddo
  
    endif
  
    deallocate(delThick)
    deallocate(delHeight_M)
    deallocate(delHeight_T)
  
    write(*,*) 'czp_tt2phi_col_ad: END'
  
    call tmg_stop(194)
  
  end subroutine czp_tt2phi_col_ad
  
  !---------------------------------------------------------
  ! calcHeightCoeff_gsv
  !---------------------------------------------------------
  subroutine calcHeightCoeff_gsv(statevector_trial)
    !
    ! :Purpose: Calculating the coefficients of height for czp_tt2phi_tl/czp_tt2phi_ad
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
  
    ! initialize and save coefficients for increased efficiency 
    !(assumes no relinearization)
    firstTimeHeightCoeff_gsv = .false.
  
    vco_anl => gsv_getVco(statevector_trial)
    Vcode_anl = vco_anl%vcode 
  
    nlev_T = gsv_getNumLev(statevector_trial,'TH')
    nlev_M = gsv_getNumLev(statevector_trial,'MM')
    numStep = statevector_trial%numstep
  
    ! saved arrays
    allocate(coeff_M_TT_gsv(&
        statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
        statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
        nlev_T,numStep))
    allocate(coeff_M_HU_gsv(&
        statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
        statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
        nlev_T,numStep))
    allocate(coeff_M_P0_delPM_gsv(&
        statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
        statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
        nlev_T,numStep))
    allocate(coeff_M_P0_dP_delPT_gsv(&
        statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
        statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
        nlev_T,numStep))
    allocate(coeff_M_P0_dP_delP0_gsv(&
        statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
        statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
        nlev_T,numStep))
  
    allocate(coeff_T_TT_gsv(&
        statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
        statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
        numStep))
    allocate(coeff_T_HU_gsv(&
        statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
        statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
        numStep))
    allocate(coeff_T_P0_delP1_gsv(&
        statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
        statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
        numStep))
    allocate(coeff_T_P0_dP_delPT_gsv(&
        statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
        statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
        numStep))
    allocate(coeff_T_P0_dP_delP0_gsv(&
        statevector_trial%myLonBeg:statevector_trial%myLonEnd, &
        statevector_trial%myLatBeg:statevector_trial%myLatEnd, &
        numStep))
  
    coeff_M_TT_gsv(:,:,:,:) = 0.0D0
    coeff_M_HU_gsv(:,:,:,:) = 0.0D0
  
    coeff_M_P0_delPM_gsv(:,:,:,:) = 0.0D0
  
    coeff_M_P0_dP_delPT_gsv(:,:,:,:) = 0.0D0
    coeff_M_P0_dP_delP0_gsv(:,:,:,:) = 0.0D0
  
    coeff_T_TT_gsv(:,:,:) = 0.0D0  
    coeff_T_HU_gsv(:,:,:) = 0.0D0  
  
    coeff_T_P0_delP1_gsv(:,:,:) = 0.0D0
  
    coeff_T_P0_dP_delPT_gsv(:,:,:) = 0.0D0  
    coeff_T_P0_dP_delP0_gsv(:,:,:) = 0.0D0  
  
    call gsv_getField(statevector_trial,hu_ptr,'HU')
    call gsv_getField(statevector_trial,tt_ptr,'TT')
    call gsv_getField(statevector_trial,P_T_ptr,'P_T')
    call gsv_getField(statevector_trial,P_M_ptr,'P_M')
    call gsv_getField(statevector_trial,height_T_ptr,'Z_T')
  
    if (Vcode_anl == 5002) then
  
      do stepIndex = 1, numStep
        do lev_T = 1, (nlev_T-1)
          do latIndex = statevector_trial%myLatBeg, statevector_trial%myLatEnd
            do lonIndex = statevector_trial%myLonBeg, statevector_trial%myLonEnd
              if ( lev_T == 1 ) then 
                ! compute height coefficients on only the top thermo level
                ratioP1 = log( P_M_ptr(lonIndex,latIndex,1,stepIndex) / &
                               P_T_ptr(lonIndex,latIndex,1,stepIndex) ) 
                hu = max(hu_ptr(lonIndex,latIndex,1,stepIndex),&
                                MPC_MINIMUM_HU_R8)
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
  
                coeff_T_TT_gsv(lonIndex,latIndex,stepIndex) = &
                    (MPC_RGAS_DRY_AIR_R8 / Rgh) * (fottva(hu,1.0D0) * cmp + &
                    fotvt8(tt,hu) * cmp_TT) * ratioP1
                
                coeff_T_HU_gsv(lonIndex,latIndex,stepIndex) = &
                    (MPC_RGAS_DRY_AIR_R8 / Rgh) * (folnqva(hu,tt,1.0d0) / & 
                    hu * cmp + fotvt8(tt,hu) * cmp_HU) * ratioP1 
  
                coeff_T_P0_delP1_gsv(lonIndex,latIndex,stepIndex) = &
                    (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp
  
                coeff_T_P0_dP_delPT_gsv(lonIndex,latIndex,stepIndex) = &
                    (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0_1 * &
                    ratioP1
  
                coeff_T_P0_dP_delP0_gsv(lonIndex,latIndex,stepIndex) = &
                    (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0_2 * &
                    ratioP1
              else
                ! compute height coefficients on momentum levels
                ratioP1 = log( P_M_ptr(lonIndex,latIndex,lev_T  ,stepIndex) / &
                               P_M_ptr(lonIndex,latIndex,lev_T-1,stepIndex) )
                hu = max( hu_ptr(lonIndex,latIndex,lev_T,stepIndex),&
                          MPC_MINIMUM_HU_R8)
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
  
                coeff_M_TT_gsv(lonIndex,latIndex,lev_T,stepIndex) = &
                    (MPC_RGAS_DRY_AIR_R8 / Rgh) * (fottva(hu,1.0D0) * cmp + &
                    fotvt8(tt,hu) * cmp_TT) * ratioP1
  
                coeff_M_HU_gsv(lonIndex,latIndex,lev_T,stepIndex) = &
                    (MPC_RGAS_DRY_AIR_R8 / RgH) * (folnqva(hu,tt,1.0d0) / hu * &
                    cmp + fotvt8(tt,hu) * cmp_HU) * ratioP1 
  
                coeff_M_P0_delPM_gsv   (lonIndex,latIndex,lev_T,stepIndex) = &
                    (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp 
  
                coeff_M_P0_dP_delPT_gsv(lonIndex,latIndex,lev_T,stepIndex) = &
                    (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0_1 * &
                    ratioP1
                
                coeff_M_P0_dP_delP0_gsv(lonIndex,latIndex,lev_T,stepIndex) = &
                    (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0_2 * &
                    ratioP1
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
              hu = max( hu_ptr(lonIndex,latIndex,lev_T,stepIndex),&
                        MPC_MINIMUM_HU_R8)
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
  
              coeff_M_TT_gsv(lonIndex,latIndex,lev_T,stepIndex) = &
                  (MPC_RGAS_DRY_AIR_R8 / Rgh) * (fottva(hu,1.0D0) * cmp + &
                  fotvt8(tt,hu) * cmp_TT) * ratioP1
  
              coeff_M_HU_gsv(lonIndex,latIndex,lev_T,stepIndex) = &
                  (MPC_RGAS_DRY_AIR_R8 / RgH) * (folnqva(hu,tt,1.0d0) / hu * &
                  cmp + fotvt8(tt,hu) * cmp_HU) * ratioP1 
  
              coeff_M_P0_delPM_gsv(lonIndex,latIndex,lev_T,stepIndex) = &
                  (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp
  
              coeff_M_P0_dP_delPT_gsv(lonIndex,latIndex,lev_T,stepIndex) = &
                  (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0_1 * &
                  ratioP1
  
              coeff_M_P0_dP_delP0_gsv(lonIndex,latIndex,lev_T,stepIndex) = &
                  (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0_2 * &
                  ratioP1
            enddo
          enddo
        enddo
      enddo
  
    end if
  
    write(*,*) "calcHeightCoeff_gsv: END"
  
  end subroutine calcHeightCoeff_gsv
  
  !---------------------------------------------------------
  ! calcHeightCoeff_col
  !---------------------------------------------------------
  subroutine calcHeightCoeff_col(columnTrlOnAnlInc)
    !
    ! :Purpose: Calculating the coefficients of height for 
    !           czp_tt2phi_tl/czp_tt2phi_ad
    !
    implicit none
  
    type(struct_columnData) :: columnTrlOnAnlInc
  
    integer :: lev_T,nlev_M,nlev_T,numColumns,colIndex,Vcode_anl
    real(8) :: hu,tt,Pr,height_T,cmp,cmp_TT,cmp_HU,cmp_P0_1,cmp_P0_2,ratioP1
    real(8) :: Rgh, sLat, lat_8
    real(8), pointer :: hu_ptr(:,:),tt_ptr(:,:)
    real(8), pointer :: P_T_ptr(:,:),P_M_ptr(:,:)
    real(8), pointer :: height_T_ptr(:,:)
    type(struct_vco), pointer :: vco_anl
  
    logical, save :: firstTimeHeightCoeff_col = .true.
  
    if ( .not. firstTimeHeightCoeff_col ) return
  
    write(*,*) "calcHeightCoeff_col: START"
  
    ! initialize and save coefficients for increased efficiency 
    ! (assumes no relinearization)
    firstTimeHeightCoeff_col = .false.
  
    vco_anl => col_getVco(columnTrlOnAnlInc)
    Vcode_anl = vco_anl%vcode 
  
    nlev_T = col_getNumLev(columnTrlOnAnlInc,'TH')
    nlev_M = col_getNumLev(columnTrlOnAnlInc,'MM')
    numColumns = col_getNumCol(columnTrlOnAnlInc)
  
    ! saved arrays
    allocate(coeff_M_TT_col         (nlev_T,numColumns))
    allocate(coeff_M_HU_col         (nlev_T,numColumns))
    allocate(coeff_M_P0_delPM_col   (nlev_T,numColumns))
    allocate(coeff_M_P0_dP_delPT_col(nlev_T,numColumns))
    allocate(coeff_M_P0_dP_delP0_col(nlev_T,numColumns))
  
    allocate(coeff_T_TT_col         (numColumns))
    allocate(coeff_T_HU_col         (numColumns))
    allocate(coeff_T_P0_delP1_col   (numColumns))
    allocate(coeff_T_P0_dP_delPT_col(numColumns))
    allocate(coeff_T_P0_dP_delP0_col(numColumns))
  
    coeff_M_TT_col(:,:) = 0.0D0
    coeff_M_HU_col(:,:) = 0.0D0
  
    coeff_M_P0_delPM_col(:,:) = 0.0D0
  
    coeff_M_P0_dP_delPT_col(:,:) = 0.0D0
    coeff_M_P0_dP_delP0_col(:,:) = 0.0D0
  
    coeff_T_TT_col(:) = 0.0D0  
    coeff_T_HU_col(:) = 0.0D0  
  
    coeff_T_P0_delP1_col(:) = 0.0D0
  
    coeff_T_P0_dP_delPT_col(:) = 0.0D0  
    coeff_T_P0_dP_delP0_col(:) = 0.0D0  
  
    hu_ptr       => col_getAllColumns(columnTrlOnAnlInc,'HU')
    tt_ptr       => col_getAllColumns(columnTrlOnAnlInc,'TT')
    P_T_ptr      => col_getAllColumns(columnTrlOnAnlInc,'P_T')
    P_M_ptr      => col_getAllColumns(columnTrlOnAnlInc,'P_M')
    height_T_ptr => col_getAllColumns(columnTrlOnAnlInc,'Z_T')
  
    if (Vcode_anl == 5002) then
  
      do colIndex = 1, numColumns
        do lev_T = 1, (nlev_T-1)
          if ( lev_T == 1 ) then 
            ! compute height coefficients on only the top thermo level
            ratioP1 = log( P_M_ptr(1,colIndex) / &
                           P_T_ptr(1,colIndex) ) 
            hu = max(hu_ptr(1,colIndex),MPC_MINIMUM_HU_R8)
            tt = tt_ptr(1,colIndex)
            Pr = P_T_ptr(1,colIndex)
            height_T = height_T_ptr(1,colIndex)
  
            cmp = gpscompressibility(Pr,tt,hu)
            cmp_TT = gpscompressibility_TT(Pr,tt,hu)
            cmp_HU = gpscompressibility_HU(Pr,tt,hu)
            cmp_P0_1 = gpscompressibility_P0_1(Pr,tt,hu,1.0d0)
            cmp_P0_2 = gpscompressibility_P0_2(Pr,tt,hu)
  
            ! Gravity acceleration 
            lat_8 = columnTrlOnAnlInc%lat(colIndex)
            sLat = sin(lat_8)
            Rgh = phf_gravityalt(sLat, height_T)
  
            coeff_T_TT_col(colIndex) = &
                (MPC_RGAS_DRY_AIR_R8 / Rgh) * (fottva(hu,1.0D0) * cmp + &
                fotvt8(tt,hu) * cmp_TT) * ratioP1
                
            coeff_T_HU_col(colIndex) = &
                (MPC_RGAS_DRY_AIR_R8 / Rgh) * (folnqva(hu,tt,1.0d0) / hu * &
                cmp + fotvt8(tt,hu) * cmp_HU) * ratioP1 
  
            coeff_T_P0_delP1_col   (colIndex) = &
                (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp
  
            coeff_T_P0_dP_delPT_col(colIndex) = &
                (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0_1 * ratioP1
  
            coeff_T_P0_dP_delP0_col(colIndex) = &
                (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0_2 * ratioP1
          else
            ! compute height coefficients on momentum levels
            ratioP1 = log( P_M_ptr(lev_T  ,colIndex) / &
                           P_M_ptr(lev_T-1,colIndex) )
            hu = max(hu_ptr(lev_T,colIndex),MPC_MINIMUM_HU_R8)
            tt = tt_ptr(lev_T,colIndex)
            Pr = P_T_ptr(lev_T,colIndex)
            height_T = height_T_ptr(lev_T,colIndex)
  
            cmp = gpscompressibility(Pr,tt,hu)
            cmp_TT = gpscompressibility_TT(Pr,tt,hu)
            cmp_HU = gpscompressibility_HU(Pr,tt,hu)
            cmp_P0_1 = gpscompressibility_P0_1(Pr,tt,hu,1.0d0)
            cmp_P0_2 = gpscompressibility_P0_2(Pr,tt,hu)
  
            ! Gravity acceleration 
            lat_8 = columnTrlOnAnlInc%lat(colIndex)
            sLat = sin(lat_8)
            Rgh = phf_gravityalt(sLat, height_T)
  
            coeff_M_TT_col(lev_T,colIndex) = &
                (MPC_RGAS_DRY_AIR_R8 / Rgh) * (fottva(hu,1.0D0) * cmp + &
                fotvt8(tt,hu) * cmp_TT) * ratioP1
  
            coeff_M_HU_col(lev_T,colIndex) = &
                (MPC_RGAS_DRY_AIR_R8 / RgH) * (folnqva(hu,tt,1.0d0) / hu * &
                cmp + fotvt8(tt,hu) * cmp_HU) * ratioP1 
  
            coeff_M_P0_delPM_col(lev_T,colIndex) = &
                (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp 
  
            coeff_M_P0_dP_delPT_col(lev_T,colIndex) = &
                (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0_1 * ratioP1
                
            coeff_M_P0_dP_delP0_col(lev_T,colIndex) = &
                (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0_2 * ratioP1
          endif
        enddo
      enddo
  
    elseif (Vcode_anl == 5005) then
  
      do colIndex = 1, numColumns
        do lev_T = 1, (nlev_T-1)
          ! compute height coefficients on momentum levels
          ratioP1 = log( P_M_ptr(lev_T+1,colIndex) / &
                         P_M_ptr(lev_T  ,colIndex) )
          hu = max(hu_ptr(lev_T,colIndex),MPC_MINIMUM_HU_R8)
          tt = tt_ptr(lev_T,colIndex)
          Pr = P_T_ptr(lev_T,colIndex)
          height_T = height_T_ptr(lev_T,colIndex)
  
          cmp = gpscompressibility(Pr,tt,hu)
          cmp_TT = gpscompressibility_TT(Pr,tt,hu)
          cmp_HU = gpscompressibility_HU(Pr,tt,hu)
          cmp_P0_1 = gpscompressibility_P0_1(Pr,tt,hu,1.0d0)
          cmp_P0_2 = gpscompressibility_P0_2(Pr,tt,hu)
  
          ! Gravity acceleration 
          lat_8 = columnTrlOnAnlInc%lat(colIndex)
          sLat = sin(lat_8)
          Rgh = phf_gravityalt(sLat, height_T)
  
          coeff_M_TT_col(lev_T,colIndex) = &
              (MPC_RGAS_DRY_AIR_R8 / Rgh) * (fottva(hu,1.0D0) * cmp + &
              fotvt8(tt,hu) * cmp_TT) * ratioP1
  
          coeff_M_HU_col(lev_T,colIndex) = &
              (MPC_RGAS_DRY_AIR_R8 / RgH) * (folnqva(hu,tt,1.0d0) / hu * cmp + &
              fotvt8(tt,hu) * cmp_HU) * ratioP1 
  
          coeff_M_P0_delPM_col   (lev_T,colIndex) = &
              (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp
  
          coeff_M_P0_dP_delPT_col(lev_T,colIndex) = &
              (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0_1 * ratioP1
  
          coeff_M_P0_dP_delP0_col(lev_T,colIndex) = &
              (MPC_RGAS_DRY_AIR_R8 / Rgh) * fotvt8(tt,hu) * cmp_P0_2 * ratioP1
        enddo
      enddo
  
    end if
  
    write(*,*) "calcHeightCoeff_col: END"
  
  end subroutine calcHeightCoeff_col
  
  !---------------------------------------------------------
  ! gpscompressibility
  !---------------------------------------------------------
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
    gpscompressibility = 1.D0 - pt * &
        (a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2) + pt*pt*(d+e*x2)
  end function gpscompressibility
  
  !---------------------------------------------------------
  ! gpscompressibility_TT
  !---------------------------------------------------------
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
    gpscompressibility_TT = &
        -d_pt * (a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2) - pt * &
        (a1*d_tc + a2*d_tc2 + b1*d_tc*x + (b0+b1*tc)*d_x + c1*d_tc*x2 + &
        (c0+c1*tc)*d_x2) + 2*pt*d_pt*(d+e*x2) + pt*pt*e*d_x2
  end function gpscompressibility_TT
  
  !---------------------------------------------------------
  ! gpscompressibility_HU
  !---------------------------------------------------------
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
    gpscompressibility_HU = &
        -d_pt * (a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2) - &
        pt * (a1*d_tc + a2*d_tc2 + b1*d_tc*x + (b0+b1*tc)*d_x + c1*d_tc*x2 + &
        (c0+c1*tc)*d_x2) + 2*pt*d_pt*(d+e*x2) + pt*pt*e*d_x2
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
    gpscompressibility_P0 = &
        -d_pt * (a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2) - &
        pt * (a1*d_tc + a2*d_tc2 + b1*d_tc*x + (b0+b1*tc)*d_x + c1*d_tc*x2 + &
        (c0+c1*tc)*d_x2) + 2*pt*d_pt*(d+e*x2) + pt*pt*e*d_x2
  end function gpscompressibility_P0
  
  !---------------------------------------------------------
  ! gpscompressibility_P0_1
  !---------------------------------------------------------
  function gpscompressibility_P0_1(p,t,q,dpdp0)
    ! gpscompressibility_P0_1 has dpdp0 dependency
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
    gpscompressibility_P0_1 = &
        -d_pt * (a0+a1*tc+a2*tc2+(b0+b1*tc)*x+(c0+c1*tc)*x2) + &
        2*pt*d_pt*(d+e*x2)
  end function gpscompressibility_P0_1
  
  !---------------------------------------------------------
  ! gpscompressibility_P0_2
  !---------------------------------------------------------
  function gpscompressibility_P0_2(p,t,q)
    ! gpscompressibility_P0_2 has NO dpdp0 dependency
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
    gpscompressibility_P0_2 = &
        -pt * (a1*d_tc + a2*d_tc2 + b1*d_tc*x + (b0+b1*tc)*d_x + c1*d_tc*x2 &
              + (c0+c1*tc)*d_x2) + pt*pt*e*d_x2
  end function gpscompressibility_P0_2

  !---------------------------------------------------------
  ! czp_calcGridPressure_nl
  !---------------------------------------------------------
  subroutine czp_calcGridPressure_nl(statevector, beSilent_opt)
    !
    ! :Purpose: calculation of the pressure on the grid subroutine
    !
    implicit none

    type(struct_gsv)    :: statevector
    logical, optional   :: beSilent_opt

    if ( gsv_getDataKind(statevector) == 4 ) then
      call calcGridPressure_nl_r4(statevector, beSilent_opt=beSilent_opt)
    else
      call calcGridPressure_nl_r8(statevector, beSilent_opt=beSilent_opt)
    end if

  end subroutine czp_calcGridPressure_nl


  !---------------------------------------------------------
  ! calcGridPressure_nl_r8
  !---------------------------------------------------------
  subroutine calcGridPressure_nl_r8(statevector, beSilent_opt)
    !
    !:Purpose: double-precision calculation of the pressure on the grid.
    !
    implicit none
  
    ! Arguments:
    type(struct_gsv), intent(inout) :: statevector ! statevector that will contain P_T/P_M
    logical, optional               :: beSilent_opt
  
    ! Locals:
    real(kind=8), allocatable   :: Psfc(:,:)
    real(kind=8), pointer       :: Pressure_out(:,:,:) 
    real(kind=8), pointer       :: field_Psfc(:,:,:,:)
    integer                     :: status, stepIndex, numStep
    logical                     :: beSilent
    real(8), pointer            :: P_T(:,:,:,:)
    real(8), pointer            :: P_M(:,:,:,:)
  
    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .true.
    end if
  
    if ( .not. beSilent ) write(*,*) 'calcGridPressure_nl_r8: computing Pressure on staggered or UNstaggered levels'
  
    if ( .not. gsv_varExist(statevector,'P_T') .or. &
        .not. gsv_varExist(statevector,'P_M') .or. &
        .not. gsv_varExist(statevector,'P0')) then
      call utl_abort('calcGridPressure_nl_r8: P_T/P_M/P0 do not exist in statevector!')
    end if
  
    nullify(P_T)
    nullify(P_M)
  
    call gsv_getField(statevector,P_T,'P_T')
    call gsv_getField(statevector,P_M,'P_M')
  
    allocate(Psfc(statevector%myLonBeg:statevector%myLonEnd, &
                  statevector%myLatBeg:statevector%myLatEnd))
    call gsv_getField(statevector,field_Psfc,'P0')
    numStep = statevector%numStep
  
    do stepIndex = 1, numStep
      Psfc(:,:) = field_Psfc(:,:,1,stepIndex)
  
      ! P_M
      nullify(Pressure_out)
      status = vgd_levels(statevector%vco%vgrid, &
                        ip1_list=statevector%vco%ip1_M, &
                        levels=Pressure_out, &
                        sfc_field=Psfc, &
                        in_log=.false.)
      if( status .ne. VGD_OK ) then
          call utl_abort('calcGridPressure_nl_r8: ERROR with vgd_levels')
      end if
      P_M(:,:,:,stepIndex) = Pressure_out(:,:,:)
      deallocate(Pressure_out)
  
      ! P_T
      nullify(Pressure_out)
      status = vgd_levels(statevector%vco%vgrid, &
                        ip1_list=statevector%vco%ip1_T, &
                        levels=Pressure_out, &
                        sfc_field=Psfc, &
                        in_log=.false.)
      if( status .ne. VGD_OK ) then
          call utl_abort('calcGridPressure_nl_r8: ERROR with vgd_levels')
      end if
      P_T(:,:,:,stepIndex) = Pressure_out(:,:,:)
      deallocate(Pressure_out)
  
      if ( .not. beSilent .and. stepIndex == 1 ) then
        write(*,*) 'stepIndex=',stepIndex, ',P_M='
        write(*,*) P_M(statevector%myLonBeg,statevector%myLatBeg,:,stepIndex)
        write(*,*) 'stepIndex=',stepIndex, ',P_T='
        write(*,*) P_T(statevector%myLonBeg,statevector%myLatBeg,:,stepIndex)
      end if
  
    end do

    deallocate(Psfc)

    if ( .not. beSilent ) write(*,*) 'calcGridPressure_nl_r8: END'
  
  end subroutine calcGridPressure_nl_r8

  !---------------------------------------------------------
  ! calcGridPressure_nl_r4
  !---------------------------------------------------------
  subroutine calcGridPressure_nl_r4(statevector, beSilent_opt)
    !
    !:Purpose: single-precision calculation of the pressure on the grid.
    !
    implicit none
  
    ! Arguments:
    type(struct_gsv), intent(inout) :: statevector ! statevector that will contain P_T/P_M
    logical, optional               :: beSilent_opt
  
    ! Locals:
    real(kind=4), allocatable   :: Psfc(:,:)
    real(kind=4), pointer       :: Pressure_out(:,:,:) 
    real(kind=4), pointer       :: field_Psfc(:,:,:,:)
    integer                     :: status, stepIndex, numStep
    logical                     :: beSilent
    real(4), pointer            :: P_T(:,:,:,:)
    real(4), pointer            :: P_M(:,:,:,:)
  
    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .true.
    end if
  
    if ( .not. beSilent ) write(*,*) 'calcGridPressure_nl_r4: computing Pressure on staggered or UNstaggered levels'
  
    if ( .not. gsv_varExist(statevector,'P_T') .or. &
        .not. gsv_varExist(statevector,'P_M') .or. &
        .not. gsv_varExist(statevector,'P0')) then
      call utl_abort('calcGridPressure_nl_r4: P_T/P_M/P0 do not exist in statevector!')
    end if
  
    nullify(P_T)
    nullify(P_M)
  
    call gsv_getField(statevector,P_T,'P_T')
    call gsv_getField(statevector,P_M,'P_M')
  
    allocate(Psfc(statevector%myLonBeg:statevector%myLonEnd, &
                  statevector%myLatBeg:statevector%myLatEnd))
    call gsv_getField(statevector,field_Psfc,'P0')
    numStep = statevector%numStep
  
    do stepIndex = 1, numStep
      Psfc(:,:) = field_Psfc(:,:,1,stepIndex)
  
      ! P_T
      nullify(Pressure_out)
      status = vgd_levels(statevector%vco%vgrid, &
                          ip1_list=statevector%vco%ip1_M, &
                          levels=Pressure_out, &
                          sfc_field=Psfc, &
                          in_log=.false.)
      if( status .ne. VGD_OK ) then
          call utl_abort('calcGridPressure_nl_r4: ERROR with vgd_levels')
      end if
      P_M(:,:,:,stepIndex) = Pressure_out(:,:,:)
      deallocate(Pressure_out)
  
      ! P_M
      nullify(Pressure_out)
      status = vgd_levels(statevector%vco%vgrid, &
                          ip1_list=statevector%vco%ip1_T, &
                          levels=Pressure_out, &
                          sfc_field=Psfc, &
                          in_log=.false.)
      if( status .ne. VGD_OK ) then
          call utl_abort('calcGridPressure_nl_r4: ERROR with vgd_levels')
      end if
      P_T(:,:,:,stepIndex) = Pressure_out(:,:,:)
      deallocate(Pressure_out)
  
      if ( .not. beSilent .and. stepIndex == 1 ) then
        write(*,*) 'stepIndex=',stepIndex, ',P_M='
        write(*,*) P_M(statevector%myLonBeg,statevector%myLatBeg,:,stepIndex)
        write(*,*) 'stepIndex=',stepIndex, ',P_T='
        write(*,*) P_T(statevector%myLonBeg,statevector%myLatBeg,:,stepIndex)
      end if
  
    end do
  
    deallocate(Psfc)

    if ( .not. beSilent ) write(*,*) 'calcGridPressure_nl_r4: END'
  
  end subroutine calcGridPressure_nl_r4

  !---------------------------------------------------------
  ! czp_calcGridPressure_tl
  !---------------------------------------------------------
  subroutine czp_calcGridPressure_tl( statevector, statevector_trial, &
                                      beSilent_opt)
    !
    !:Purpose: calculation of the Pressure increment on the grid.
    !
    implicit none
  
    ! Arguments:
    type(struct_gsv), intent(inout) :: statevector       ! statevector that will contain the P_T/P_M increments
    type(struct_gsv), intent(in)    :: statevector_trial ! statevector that has the Psfc
    logical, optional               :: beSilent_opt
  
    ! Locals:
    real(8), allocatable  :: Psfc(:,:)
    real(4), pointer      :: delPsfc_r4(:,:,:,:)
    real(8), pointer      :: delPsfc_r8(:,:,:,:)
    real(8), pointer      :: field_Psfc(:,:,:,:)
    real(4), pointer      :: delP_T_r4(:,:,:,:)
    real(8), pointer      :: delP_T_r8(:,:,:,:)
    real(4), pointer      :: delP_M_r4(:,:,:,:)
    real(8), pointer      :: delP_M_r8(:,:,:,:)
    real(8), pointer      :: dP_dPsfc_T(:,:,:)
    real(8), pointer      :: dP_dPsfc_M(:,:,:)
    integer               :: status, stepIndex,lonIndex,latIndex
    integer               :: lev_M, lev_T, nlev_T, nlev_M, numStep
    logical               :: beSilent
  
    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .true.
    end if
  
    if ( .not. beSilent ) write(*,*) 'czp_calcGridPressure_tl: computing delP_T/delP_M on the gridstatevector'
  
    nullify(dP_dPsfc_T)
    nullify(dP_dPsfc_M)
    nullify(delPsfc_r4,delPsfc_r8)
    nullify(delP_T_r4,delP_T_r8)
    nullify(delP_M_r4,delP_M_r8)
  
    if (gsv_getDataKind(statevector) == 4) then
      call gsv_getField(statevector,delP_T_r4,'P_T')
      call gsv_getField(statevector,delP_M_r4,'P_M')
      call gsv_getField(statevector,delPsfc_r4,'P0')
    else
      call gsv_getField(statevector,delP_T_r8,'P_T')
      call gsv_getField(statevector,delP_M_r8,'P_M')
      call gsv_getField(statevector,delPsfc_r8,'P0')
    end if
  
    nullify(field_Psfc)
    call gsv_getField(statevector_trial,field_Psfc,'P0')
  
    nlev_T = gsv_getNumLev(statevector,'TH')
    nlev_M = gsv_getNumLev(statevector,'MM')
    numStep = statevector%numstep
  
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
      if( status .ne. VGD_OK ) then
          call utl_abort('czp_calcGridPressure_tl: ERROR with vgd_dpidpis')
      end if
      ! calculate delP_M
      if (gsv_getDataKind(statevector) == 4) then
        do lev_M = 1, nlev_M
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd
              delP_M_r4(lonIndex,latIndex,lev_M,stepIndex) =  &
                   dP_dPsfc_M(lonIndex-statevector%myLonBeg+1,&
                              latIndex-statevector%myLatBeg+1,lev_M) * &
                   delPsfc_r4(lonIndex,latIndex,1,stepIndex)
            end do
          end do
        end do
      else
        do lev_M = 1, nlev_M
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd
              delP_M_r8(lonIndex,latIndex,lev_M,stepIndex) =  &
                   dP_dPsfc_M(lonIndex-statevector%myLonBeg+1,&
                              latIndex-statevector%myLatBeg+1,lev_M) * &
                   delPsfc_r8(lonIndex,latIndex,1,stepIndex)
            end do
          end do
        end do
      end if
      deallocate(dP_dPsfc_M)
  
      ! dP_dPsfc_T
      nullify(dP_dPsfc_T)
      status = vgd_dpidpis(statevector%vco%vgrid, &
                           statevector%vco%ip1_T, &
                           dP_dPsfc_T, &
                           Psfc)
      if( status .ne. VGD_OK ) then
          call utl_abort('czp_calcGridPressure_tl: ERROR with vgd_dpidpis')
      end if
      ! calculate delP_T
      if (gsv_getDataKind(statevector) == 4) then
        do lev_T = 1, nlev_T
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd
              delP_T_r4(lonIndex,latIndex,lev_T,stepIndex) =  &
                   dP_dPsfc_T(lonIndex-statevector%myLonBeg+1,&
                              latIndex-statevector%myLatBeg+1,lev_T) * &
                   delPsfc_r4(lonIndex,latIndex,1,stepIndex)
            end do
          end do
        end do
      else
        do lev_T = 1, nlev_T
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd
              delP_T_r8(lonIndex,latIndex,lev_T,stepIndex) =  &
                   dP_dPsfc_T(lonIndex-statevector%myLonBeg+1,&
                              latIndex-statevector%myLatBeg+1,lev_T) * &
                   delPsfc_r8(lonIndex,latIndex,1,stepIndex)
            end do
          end do
        end do
      end if
      deallocate(dP_dPsfc_T)
  
    end do
  
    deallocate(Psfc)
  
  end subroutine czp_calcGridPressure_tl

  !---------------------------------------------------------
  ! czp_calcGridPressure_ad
  !---------------------------------------------------------
  subroutine czp_calcGridPressure_ad( statevector, statevector_trial, &
                                      beSilent_opt)
    !
    !:Purpose: adjoint of calculation of the Pressure on the grid.
    !
    implicit none
  
    ! Arguments:
    type(struct_gsv), intent(inout) :: statevector       ! statevector that will contain increment of Psfc.
    type(struct_gsv), intent(in)    :: statevector_trial ! statevector that has the Psfc
    logical, optional               :: beSilent_opt
  
    ! Locals:
    real(8), allocatable     :: Psfc(:,:)
    real(4), pointer         :: delPsfc_r4(:,:,:,:)
    real(8), pointer         :: delPsfc_r8(:,:,:,:)
    real(8), pointer         :: field_Psfc(:,:,:,:)
    real(4), pointer         :: delP_T_r4(:,:,:,:)
    real(8), pointer         :: delP_T_r8(:,:,:,:)
    real(4), pointer         :: delP_M_r4(:,:,:,:)
    real(8), pointer         :: delP_M_r8(:,:,:,:)
    real(8), pointer         :: dP_dPsfc_T(:,:,:)
    real(8), pointer         :: dP_dPsfc_M(:,:,:)
    integer                  :: status, stepIndex,lonIndex,latIndex
    integer                  :: lev_M, lev_T, nlev_T, nlev_M, numStep
    logical                  :: beSilent
  
    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .true.
    end if
  
    if ( .not. beSilent ) write(*,*) 'czp_calcGridPressure_ad: computing delP_T/delP_M on the gridstatevector'
  
    nullify(delPsfc_r4, delPsfc_r8)
    nullify(field_Psfc)
    nullify(delP_T_r4, delP_T_r8)
    nullify(delP_M_r4, delP_M_r8)
    nullify(dP_dPsfc_T)
    nullify(dP_dPsfc_M)
  
    if (gsv_getDataKind(statevector) == 4) then
      call gsv_getField(statevector,delP_T_r4,'P_T')
      call gsv_getField(statevector,delP_M_r4,'P_M')
      call gsv_getField(statevector,delPsfc_r4,'P0')
    else
      call gsv_getField(statevector,delP_T_r8,'P_T')
      call gsv_getField(statevector,delP_M_r8,'P_M')
      call gsv_getField(statevector,delPsfc_r8,'P0')
    end if
    call gsv_getField(statevector_trial,field_Psfc,'P0')
  
    nlev_T = gsv_getNumLev(statevector,'TH')
    nlev_M = gsv_getNumLev(statevector,'MM')
    numStep = statevector%numstep
  
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
      if( status .ne. VGD_OK ) call utl_abort('czp_calcGridPressure_ad: ERROR with vgd_dpidpis')
      ! calculate delP_M
      if (gsv_getDataKind(statevector) == 4) then
        do lev_M = 1, nlev_M
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd
              delPsfc_r4(lonIndex,latIndex,1,stepIndex) =  &
                   delPsfc_r4(lonIndex,latIndex,1,stepIndex) + &
                   dP_dPsfc_M(lonIndex-statevector%myLonBeg+1,&
                              latIndex-statevector%myLatBeg+1,lev_M) * &
                   delP_M_r4(lonIndex,latIndex,lev_M,stepIndex)
            end do
          end do
        end do
      else
        do lev_M = 1, nlev_M
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd
              delPsfc_r8(lonIndex,latIndex,1,stepIndex) =  &
                   delPsfc_r8(lonIndex,latIndex,1,stepIndex) + &
                   dP_dPsfc_M(lonIndex-statevector%myLonBeg+1,&
                              latIndex-statevector%myLatBeg+1,lev_M) * &
                   delP_M_r8(lonIndex,latIndex,lev_M,stepIndex)
            end do
          end do
        end do
      end if
      deallocate(dP_dPsfc_M)
  
      ! dP_dPsfc_T
      nullify(dP_dPsfc_T)
      status = vgd_dpidpis(statevector%vco%vgrid, &
                           statevector%vco%ip1_T, &
                           dP_dPsfc_T, &
                           Psfc)
      if( status .ne. VGD_OK ) then
          call utl_abort('czp_calcGridPressure_ad: ERROR with vgd_dpidpis')
      end if
      ! calculate delP_T
      if (gsv_getDataKind(statevector) == 4) then
        do lev_T = 1, nlev_T
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd
              delPsfc_r4(lonIndex,latIndex,1,stepIndex) =  &
                   delPsfc_r4(lonIndex,latIndex,1,stepIndex) + &
                   dP_dPsfc_T(lonIndex-statevector%myLonBeg+1,&
                              latIndex-statevector%myLatBeg+1,lev_T) * &
                   delP_T_r4(lonIndex,latIndex,lev_T,stepIndex)
            end do
          end do
        end do
      else
        do lev_T = 1, nlev_T
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd
              delPsfc_r8(lonIndex,latIndex,1,stepIndex) =  &
                   delPsfc_r8(lonIndex,latIndex,1,stepIndex) + &
                   dP_dPsfc_T(lonIndex-statevector%myLonBeg+1,&
                              latIndex-statevector%myLatBeg+1,lev_T) * &
                   delP_T_r8(lonIndex,latIndex,lev_T,stepIndex)
            end do
          end do
        end do
      end if
      deallocate(dP_dPsfc_T)
  
    end do
  
    deallocate(Psfc)
  
  end subroutine czp_calcGridPressure_ad

  !---------------------------------------------------------
  ! czp_calcColumnZandP_tl
  !---------------------------------------------------------
  subroutine czp_calcColumnZandP_tl(columnInc, columnRefOnIncLev, beSilent_opt)
    !
    ! :Purpose: compute pressure and height increment in the column in proper 
    !           order depending on the vgrid kind
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(inout) :: columnInc    ! column that will contain the P_T/P_M increments
    type(struct_columnData), intent(in)    :: columnRefOnIncLev ! column that has the Psfc
    logical, intent(in), optional          :: beSilent_opt

    ! Locals
    integer   :: Vcode

    Vcode = columnInc%vco%vcode
    if (Vcode == 5002 .or. Vcode == 5005) then
      if (col_varExist(columnInc,'P_T') .and. &
          col_varExist(columnInc,'P_M') ) then
        call czp_calcColumnPressure_tl( columnInc, columnRefOnIncLev, &
                                        beSilent_opt=beSilent_opt)
      end if
      if (col_varExist(columnInc,'Z_T') .and. &
          col_varExist(columnInc,'Z_M') ) then
        call czp_tt2phi_col_tl(columnInc, columnRefOnIncLev)
      end if
    else if (Vcode == 21001) then
      if (col_varExist(columnInc,'Z_T') .and. &
          col_varExist(columnInc,'Z_M') ) then
        call czp_tt2phi_col_tl(columnInc, columnRefOnIncLev)
      end if
      if (col_varExist(columnInc,'P_T') .and. &
          col_varExist(columnInc,'P_M') ) then
        call czp_calcColumnPressure_tl( columnInc, columnRefOnIncLev, &
                                        beSilent_opt=beSilent_opt)
      end if
    else
      call utl_abort('czp_calcColumnZandP_tl: invalid Vcode')
    end if

  end subroutine czp_calcColumnZandP_tl

  !---------------------------------------------------------
  ! czp_calcColumnZandP_ad
  !---------------------------------------------------------
  subroutine czp_calcColumnZandP_ad(columnInc, columnRefOnIncLev, beSilent_opt)
    !
    ! :Purpose: adjoint of pressure and height increment computation in the 
    !           column in proper order depending on the vgrid kind
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(inout) :: columnInc    ! column that will contain the P_T/P_M increments
    type(struct_columnData), intent(in)    :: columnRefOnIncLev ! column that has the Psfc
    logical, intent(in), optional          :: beSilent_opt

    ! Locals
    integer   :: Vcode

    Vcode = columnInc%vco%vcode
    if (Vcode == 5002 .or. Vcode == 5005) then
      if (col_varExist(columnInc,'Z_T') .and. &
          col_varExist(columnInc,'Z_M') ) then
        call czp_tt2phi_col_ad(columnInc, columnRefOnIncLev)
      end if
      if (col_varExist(columnInc,'P_T') .and. &
          col_varExist(columnInc,'P_M') ) then
        call czp_calcColumnPressure_ad( columnInc, columnRefOnIncLev, &
                                      beSilent_opt=beSilent_opt)
      end if
    else if (Vcode == 21001) then
      if (col_varExist(columnInc,'P_T') .and. &
          col_varExist(columnInc,'P_M') ) then
        call czp_calcColumnPressure_ad( columnInc, columnRefOnIncLev, &
                                        beSilent_opt=beSilent_opt)
      end if
      if (col_varExist(columnInc,'Z_T') .and. &
          col_varExist(columnInc,'Z_M') ) then
        call czp_tt2phi_col_ad(columnInc, columnRefOnIncLev)
      end if
    else
      call utl_abort('czp_calcColumnZandP_ad: invalid Vcode')
    end if

  end subroutine czp_calcColumnZandP_ad

  !---------------------------------------------------------
  ! czp_calcColumnPressure_nl
  !---------------------------------------------------------
  subroutine czp_calcColumnPressure_nl(column, beSilent_opt)
    !
    !:Purpose: calculation of the Pressure in the column.
    !
    implicit none
    type(struct_columnData), intent(inout) :: column
    logical, optional :: beSilent_opt
  
    real(kind=8), allocatable :: Psfc(:,:),zppobs2(:,:)
    real(kind=8), pointer     :: zppobs1(:,:,:) => null()
    integer :: headerIndex, status, ilev1, ilev2
    logical                   :: beSilent
  
    if ( col_getNumCol(column) <= 0 ) return
  
    if (.not.col_varExist(column,'P0')) then
      call utl_abort('czp_calcColumnPressure_nl: P0 must be present as an analysis variable!')
    end if
  
    allocate(Psfc(1,col_getNumCol(column)))
    do headerIndex = 1,col_getNumCol(column)
      Psfc(1,headerIndex) = col_getElem(column,1,headerIndex,'P0')
    end do
  
    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .false.
    end if
  
    if ( .not.beSilent ) write(*,*) 'czp_calcColumnPressure_nl: computing Pressure on staggered or UNstaggered levels'
  
    status=vgd_levels(column%vco%vgrid,ip1_list=column%vco%ip1_M,  &
                      levels=zppobs1,sfc_field=Psfc,in_log=.false.)
    if(status.ne.VGD_OK) call utl_abort('ERROR with vgd_levels')
    allocate(zppobs2(col_getNumLev(column,'MM'),col_getNumCol(column)))
    zppobs2 = transpose(zppobs1(1,:,:))
    ilev1 = 1 + column%varOffset(vnl_varListIndex('P_M'))
    ilev2 = ilev1 - 1 + column%varNumLev(vnl_varListIndex('P_M'))
    column%all(ilev1:ilev2,:) = zppobs2(:,:)
    if (associated(zppobs1))  deallocate(zppobs1)
    deallocate(zppobs2)
  
    status=vgd_levels(column%vco%vgrid,ip1_list=column%vco%ip1_T,  &
                      levels=zppobs1,sfc_field=Psfc,in_log=.false.)
    if(status.ne.VGD_OK) call utl_abort('ERROR with vgd_levels')
    allocate(zppobs2(col_getNumLev(column,'TH'),col_getNumCol(column)))
    zppobs2 = transpose(zppobs1(1,:,:))
    ilev1 = 1 + column%varOffset(vnl_varListIndex('P_T'))
    ilev2 = ilev1 - 1 + column%varNumLev(vnl_varListIndex('P_T'))
    column%all(ilev1:ilev2,:) = zppobs2(:,:)
    if (associated(zppobs1)) deallocate(zppobs1)
    deallocate(zppobs2)
  
    deallocate(Psfc)
  
  end subroutine czp_calcColumnPressure_nl

  !---------------------------------------------------------
  ! czp_calcColumnPressure_tl
  !---------------------------------------------------------
  subroutine czp_calcColumnPressure_tl( columnInc, columnRefOnIncLev, &
                                        beSilent_opt)
    !
    !:Purpose: calculation of the Pressure increment in the column.
    !
    implicit none
  
    ! Arguments:
    type(struct_columnData), intent(inout) :: columnInc    ! column that will contain the P_T/P_M increments
    type(struct_columnData), intent(in)    :: columnRefOnIncLev ! column that has the Psfc
    logical, intent(in), optional          :: beSilent_opt
  
    ! Locals:
    real(8)          :: Psfc
    real(8), pointer :: delPsfc(:,:), PsfcRef(:,:)
    real(8), pointer :: delP_T(:,:), delP_M(:,:)
    real(8), pointer :: dP_dPsfc_T(:), dP_dPsfc_M(:)
    integer          :: status, colIndex
    integer          :: lev_M, lev_T, nlev_T, nlev_M, numColumns
    logical          :: beSilent
  
    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .true.
    end if
  
    if ( .not. beSilent ) write(*,*) 'czp_calcColumnPressure_tl: computing delP_T/delP_M on the column'
  
    nullify(dP_dPsfc_T)
    nullify(dP_dPsfc_M)
    nullify(delPsfc)
    nullify(delP_T)
    nullify(delP_M)
  
    delP_M  => col_getAllColumns(columnInc,'P_M')
    delP_T  => col_getAllColumns(columnInc,'P_T')
    delPsfc => col_getAllColumns(columnInc,'P0')
    PsfcRef => col_getAllColumns(columnRefOnIncLev,'P0')
  
    nlev_T = col_getNumLev(columnInc,'TH')
    nlev_M = col_getNumLev(columnInc,'MM')
    numColumns = col_getNumCol(columnInc)
  
    do colIndex = 1, numColumns
  
      Psfc = PsfcRef(1,colIndex)
  
      ! dP_dPsfc_M
      nullify(dP_dPsfc_M)
      status = vgd_dpidpis(columnInc%vco%vgrid, &
                           columnInc%vco%ip1_M, &
                           dP_dPsfc_M, &
                           Psfc)
      if( status .ne. VGD_OK ) then
          call utl_abort('czp_calcColumnPressure_tl: ERROR with vgd_dpidpis')
      end if
      ! calculate delP_M
      do lev_M = 1, nlev_M
        delP_M(lev_M,colIndex) = dP_dPsfc_M(lev_M) * delPsfc(1,colIndex)
      end do
      deallocate(dP_dPsfc_M)
  
      ! dP_dPsfc_T
      nullify(dP_dPsfc_T)
      status = vgd_dpidpis(columnInc%vco%vgrid, &
                           columnInc%vco%ip1_T, &
                           dP_dPsfc_T, &
                           Psfc)
      if( status .ne. VGD_OK ) then
          call utl_abort('czp_calcColumnPressure_tl: ERROR with vgd_dpidpis')
      end if
      ! calculate delP_T
      do lev_T = 1, nlev_T
        delP_T(lev_T,colIndex) = dP_dPsfc_T(lev_T) * delPsfc(1,colIndex)
      end do
      deallocate(dP_dPsfc_T)
  
    end do
  
  end subroutine czp_calcColumnPressure_tl

  !---------------------------------------------------------
  ! czp_calcColumnPressure_ad
  !---------------------------------------------------------
  subroutine czp_calcColumnPressure_ad( columnInc, columnRefOnIncLev, &
                                        beSilent_opt)
    !
    !:Purpose: adjoint of calculation of the Pressure in the column.
    !
    implicit none
  
    ! Arguments:
    type(struct_columnData), intent(inout) :: columnInc    ! column that will contain increments of Psfc.
    type(struct_columnData), intent(in)    :: columnRefOnIncLev ! column that has the Psfc
    logical, optional                      :: beSilent_opt
  
    ! Locals:
    real(8)          :: Psfc
    real(8), pointer :: delPsfc(:,:), PsfcRef(:,:)
    real(8), pointer :: delP_T(:,:), delP_M(:,:)
    real(8), pointer :: dP_dPsfc_T(:), dP_dPsfc_M(:)
    integer          :: status, colIndex
    integer          :: lev_M, lev_T, nlev_T, nlev_M, numColumns
    logical          :: beSilent
  
    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .true.
    end if
  
    if ( .not. beSilent ) write(*,*) 'czp_calcColumnPressure_ad: computing delP_T/delP_M on the column'
  
    nullify(delPsfc)
    nullify(PsfcRef)
    nullify(delP_T)
    nullify(delP_M)
    nullify(dP_dPsfc_T)
    nullify(dP_dPsfc_M)
  
    delP_M  => col_getAllColumns(columnInc,'P_M')
    delP_T  => col_getAllColumns(columnInc,'P_T')
    delPsfc => col_getAllColumns(columnInc,'P0')
    PsfcRef => col_getAllColumns(columnRefOnIncLev,'P0')
  
    nlev_T = col_getNumLev(columnInc,'TH')
    nlev_M = col_getNumLev(columnInc,'MM')
    numColumns = col_getNumCol(columnInc)
  
    do colIndex = 1, numColumns
  
      Psfc = PsfcRef(1,colIndex)
  
      ! dP_dPsfc_M
      nullify(dP_dPsfc_M)
      status = vgd_dpidpis(columnInc%vco%vgrid, &
                           columnInc%vco%ip1_M, &
                           dP_dPsfc_M, &
                           Psfc)
      if( status .ne. VGD_OK ) then
          call utl_abort('czp_calcColumnPressure_ad: ERROR with vgd_dpidpis')
      end if
      ! calculate delP_M
      do lev_M = 1, nlev_M
        delPsfc(1,colIndex) = delPsfc(1,colIndex) + &
             dP_dPsfc_M(lev_M) * delP_M(lev_M,colIndex)
      end do
      deallocate(dP_dPsfc_M)
  
      ! dP_dPsfc_T
      nullify(dP_dPsfc_T)
      status = vgd_dpidpis(columnInc%vco%vgrid, &
                           columnInc%vco%ip1_T, &
                           dP_dPsfc_T, &
                           Psfc)
      if( status .ne. VGD_OK ) call utl_abort('czp_calcColumnPressure_ad: ERROR with vgd_dpidpis')
      ! calculate delP_T
      do lev_T = 1, nlev_T
        delPsfc(1,colIndex) = &
            delPsfc(1,colIndex) + dP_dPsfc_T(lev_T) * delP_T(lev_T,colIndex)
      end do
      deallocate(dP_dPsfc_T)
  
    end do
  
  end subroutine czp_calcColumnPressure_ad

end module calcHeightAndPressure_mod
