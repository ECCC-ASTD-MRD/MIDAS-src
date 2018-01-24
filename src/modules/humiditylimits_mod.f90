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
!! MODULE humidityLimits_mod (prefix="qlim")
!!
!! *Purpose*: Various manipulations of humidity-related quantities.
!!
!--------------------------------------------------------------------------
module humidityLimits_mod
  use mpi_mod
  use utilities_mod
  use mathPhysConstants_mod
  use physicsFunctions_mod
  use verticalCoord_mod
  use gridStateVector_mod
  implicit none
  save
  private
  
  ! public procedures
  public :: qlim_gsvSaturationLimit, qlim_gsvRttovLimit

  real(8), parameter :: mixratio_to_ppmv = 1.60771704d+6

  contains

  subroutine qlim_gsvSaturationLimit(statevector, HUcontainsLQ_opt)
    !
    ! Purpose:
    !          impose saturation limit on humidity variable of a statevector
    !
    implicit none
    type(struct_gsv) :: statevector
    logical, optional :: HUcontainsLQ_opt

    type(struct_vco), pointer :: vco_ptr
    logical :: HUcontainsLQ
    real(8), pointer :: lq_ptr(:,:,:,:), hu_ptr(:,:,:,:), tt_ptr(:,:,:,:), psfc_ptr(:,:,:,:)
    real(8), pointer :: pressure(:,:,:)
    real(8)          :: hu, husat, hu_modified, tt
    real(8), allocatable :: psfc(:,:)
    integer          :: lon1, lon2, lat1, lat2, lev1, lev2, ierr
    integer          :: lonIndex, latIndex, levIndex, stepIndex

    write(*,*) 'qlim_gsvSaturationLimit: STARTING'

    if( .not. gsv_varExist(statevector,'HU') ) then
      if( mpi_myid == 0 ) write(*,*) 'qlim_gsvSaturationLimit: statevector does not ' // &
                                     'contain humidity ... doing nothing'
      return
    end if

    if( statevector%dataKind /= 8 ) then
      call utl_abort('qlim_gsvSaturationLimit: only compatible with double precision ' // &
                     'data.')
    end if

    if( present(HUcontainsLQ_opt) ) then
      HUcontainsLQ = HUcontainsLQ_opt
    else
      HUcontainsLQ = .true.
    end if

    vco_ptr => gsv_getVco(statevector)
    lq_ptr => gsv_getField_r8(statevector,'HU')
    hu_ptr => gsv_getField_r8(statevector,'HU')
    tt_ptr => gsv_getField_r8(statevector,'TT')

    lon1 = statevector%myLonBeg
    lon2 = statevector%myLonEnd
    lat1 = statevector%myLatBeg
    lat2 = statevector%myLatEnd
    lev1 = 1
    lev2 = vco_getNumLev(vco_ptr,'TH')

    allocate(psfc(lon2-lon1+1,lat2-lat1+1))
    do stepIndex = 1, statevector%numStep
      psfc_ptr => gsv_getField_r8(statevector,'P0')
      psfc(:,:) = psfc_ptr(:,:,1,stepIndex)
      nullify(pressure)
      ierr = vgd_levels(vco_ptr%vgrid,           &
                        ip1_list=vco_ptr%ip1_T,  &
                        levels=pressure,         &
                        sfc_field=psfc, &
                        in_log=.false.)
!$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex, hu, tt, husat, hu_modified)
      do levIndex = lev1, lev2
        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2

            ! Obtain HU from LQ, if necessary
            if( HUcontainsLQ ) then
              hu = exp(lq_ptr(lonIndex,latIndex,levIndex,stepIndex))
            else
              hu = hu_ptr(lonIndex,latIndex,levIndex,stepIndex)
            end if
            tt = tt_ptr(lonIndex,latIndex,levIndex,stepIndex)

            ! get the saturated vapor pressure from HU
            husat = foqst8(tt + MPC_K_C_DEGREE_OFFSET_R8, pressure(lonIndex-lon1+1,latIndex-lat1+1,levIndex) )

            ! limit the humidity to the saturated humidity
            hu_modified = min(husat, hu)

            ! Recompute LQ from HU, if necessary
            if( HUcontainsLQ ) then
              if( hu /= hu_modified ) then
                lq_ptr(lonIndex,latIndex,levIndex,stepIndex) = log(max(hu_modified,MPC_MINIMUM_HU_R8))
              end if
            else
              hu_ptr(lonIndex,latIndex,levIndex,stepIndex) = hu_modified
            end if

          end do ! lonIndex
        end do ! latIndex
      end do ! levIndex
!$OMP END PARALLEL DO

      deallocate(pressure)

    end do ! stepIndex

  end subroutine qlim_gsvSaturationLimit


  subroutine qlim_gsvRttovLimit(statevector, HUcontainsLQ_opt)
    !
    ! Purpose:
    !          impose RTTOV limits on humidity
    !
    implicit none
    type(struct_gsv) :: statevector
    logical, optional :: HUcontainsLQ_opt
    type(struct_vco), pointer :: vco_ptr
    logical :: HUcontainsLQ
    real(8), pointer :: lq_ptr(:,:,:,:), hu_ptr(:,:,:,:), psfc_ptr(:,:,:,:)
    real(8), pointer :: pressure(:,:,:)
    real(8)          :: hu, hu_modified
    real(8), allocatable :: psfc(:,:)
    integer          :: lon1, lon2, lat1, lat2, lev1, lev2, ierr, nulfile
    integer          :: lonIndex, latIndex, levIndex, stepIndex, numLev_rttov
    integer          :: ni, nj, numLev
    real(8), allocatable :: press_rttov(:), qmin_rttov(:), qmax_rttov(:)
    real(8), allocatable :: qmin3D_rttov(:,:,:), qmax3D_rttov(:,:,:)
    character(len=256) :: fileName
    integer :: fnom, fclos

    write(*,*) 'qlim_gsvRttovLimit: STARTING'

    if( .not. gsv_varExist(statevector,'HU') ) then
      if( mpi_myid == 0 ) write(*,*) 'qlim_gsvRttovLimit: statevector does not ' // &
                                     'contain humidity ... doing nothing'
      return
    end if

    if( statevector%dataKind /= 8 ) then
      call utl_abort('qlim_gsvRttovLimit: only compatible with double precision ' // &
                     'data.')
    end if

    if( present(HUcontainsLQ_opt) ) then
      HUcontainsLQ = HUcontainsLQ_opt
    else
      HUcontainsLQ = .true.
    end if

    ! Read in RTTOV humidity limits
    fileName = "rttov_h2o_limits.dat"
    nulfile = 0
    ierr = fnom(nulfile, fileName, "FMT+OLD+R/O", 0)
    if( ierr /= 0 ) then
      write(*,*) 'fileName = ', fileName
      call utl_abort('qlim_gsvRttovLimit: error opening the humidity limits file')
    end if

    read(nulfile,*) numLev_rttov
    write(*,*) 'qlim_gsvRttovLimit: rttov number of levels = ', numLev_rttov
    allocate(press_rttov(numLev_rttov))
    allocate(qmin_rttov(numLev_rttov))
    allocate(qmax_rttov(numLev_rttov))
    do levIndex = 1, numLev_rttov
      read(nulfile,*) press_rttov(levIndex), qmax_rttov(levIndex), qmin_rttov(levIndex)
    end do
    ierr = fclos(nulfile)
    press_rttov(:) = press_rttov(:) * mpc_pa_per_mbar_r8
    qmin_rttov(:) = qmin_rttov(:) / mixratio_to_ppmv
    qmax_rttov(:) = qmax_rttov(:) / mixratio_to_ppmv

    write(*,*) ' '
    do levIndex = 1, numLev_rttov
      write(*,fmt='(" qlim_gsvRttovLimit:   LEVEL = ",I4,", PRES = ",F9.0,", HUMIN = ",E10.2,", HUMAX = ",E10.2)') &
        levIndex, press_rttov(levIndex), qmin_rttov(levIndex), qmax_rttov(levIndex)
    end do

    vco_ptr => gsv_getVco(statevector)
    lq_ptr => gsv_getField_r8(statevector,'HU')
    hu_ptr => gsv_getField_r8(statevector,'HU')

    lon1 = statevector%myLonBeg
    lon2 = statevector%myLonEnd
    lat1 = statevector%myLatBeg
    lat2 = statevector%myLatEnd
    lev1 = 1
    lev2 = vco_getNumLev(vco_ptr,'TH')

    ni = lon2 - lon1 + 1
    nj = lat2 - lat1 + 1
    numLev = lev2 - lev1 + 1
    allocate( qmin3D_rttov(ni,nj,numLev) )
    allocate( qmax3D_rttov(ni,nj,numLev) )
    allocate( psfc(lon2-lon1+1,lat2-lat1+1) )

    do stepIndex = 1, statevector%numStep
      psfc_ptr => gsv_getField_r8(statevector,'P0')
      psfc(:,:) = psfc_ptr(:,:,1,stepIndex)
      nullify(pressure)
      ierr = vgd_levels(vco_ptr%vgrid,           &
                        ip1_list=vco_ptr%ip1_T,  &
                        levels=pressure,         &
                        sfc_field=psfc, &
                        in_log=.false.)

      ! Interpolate RTTOV limits onto model levels
      call lintv_minmax(press_rttov, qmin_rttov, qmax_rttov, numLev_rttov, &
                        ni, nj, numLev, pressure, qmin3D_rttov, qmax3D_rttov)

!$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex, hu, hu_modified)
      do levIndex = lev1, lev2
        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2

            ! Obtain HU from LQ, if necessary
            if( HUcontainsLQ ) then
              hu = exp(lq_ptr(lonIndex,latIndex,levIndex,stepIndex))
            else
              hu = hu_ptr(lonIndex,latIndex,levIndex,stepIndex)
            end if

            ! limit the humidity according to the rttov limits
            hu_modified = max(hu, qmin3D_rttov(lonIndex - lon1 + 1, latIndex - lat1 + 1, levIndex) )
            hu_modified = min(hu_modified, qmax3D_rttov(lonIndex - lon1 + 1, latIndex - lat1 + 1, levIndex) )

            ! Recompute LQ from HU, if necessary
            if( HUcontainsLQ ) then
              if( hu /= hu_modified ) then
                lq_ptr(lonIndex,latIndex,levIndex,stepIndex) = log(max(hu_modified,MPC_MINIMUM_HU_R8))
              end if
            else
              hu_ptr(lonIndex,latIndex,levIndex,stepIndex) = hu_modified
            end if

          end do ! lonIndex
        end do ! latIndex
      end do ! levIndex
!$OMP END PARALLEL DO

      deallocate(pressure)

    end do ! stepIndex

    deallocate( psfc )
    deallocate( qmin3D_rttov )
    deallocate( qmax3D_rttov )

  end subroutine qlim_gsvRttovLimit


  subroutine lintv_minmax(press_src, qmin_src, qmax_src, numLev_src, &
                          ni_dest, nj_dest, numLev_dest, press_dest, qmin_dest, qmax_dest)
    !Arguments
    !     i   press_src(numLev_src)    : Vertical levels, pressure (source)
    !     i   qmin/max_src(numLev_src)   : Vectors to be interpolated (source)
    !     i   numLev_src               : Number of input levels (source)
    !     i   ni_dest,nj_dest          : Number of profiles
    !     i   numLev_dest              : Number of output levels (destination)
    !     i   press_dest(ni_dest,nj_dest,numLev_dest) : Vertical levels, pressure (destination)
    !     o   qmin/max_dest(ni_dest,nj_dest,numLev_dest) : Interpolated profiles (destination)
    !
    !    -------------------
    !*    Purpose: Performs the vertical interpolation in log of pressure
    !              and constant value extrapolation of one-dimensional vectors.
    implicit none

    integer, intent(in) :: numLev_src, ni_dest, nj_dest, numLev_dest
    real(8), intent(in) :: press_src(numLev_src)
    real(8), intent(in) :: qmin_src(numLev_src) ,qmax_src(numLev_src)
    real(8)             :: press_dest(:,:,:), qmin_dest(:,:,:), qmax_dest(:,:,:)

    integer :: ji, jk, jo, ii, jj, ik, iorder
    integer :: ilen, ierr
    real(8) :: zpo(numLev_dest)
    integer :: il(numLev_dest)
    real(8) :: zpi(0:numLev_src+1)
    real(8) :: zqmin_src(0:numLev_src+1), zqmax_src(0:numLev_src+1)
    real(8) :: zw1, zw2, zp, xi, zrt, zp1, zp2

    zpi(0)=200000.d0
    zpi(numLev_src+1)=200000.d0

    ! Determine if input pressure levels are in ascending or
    ! descending order.
    if ( press_src(1) < press_src(numLev_src) ) then
      iorder = 1
    else
      iorder = -1
    end if

    ! Source levels
    do jk = 1, numLev_src
      zpi(jk) = press_src(jk)
      zqmin_src(jk) = qmin_src(jk)
      zqmax_src(jk) = qmax_src(jk)
    enddo
    zqmin_src(0) = qmin_src(1)
    zqmin_src(numLev_src+1) = qmin_src(numLev_src)
    zqmax_src(0) = qmax_src(1)
    zqmax_src(numLev_src+1) = qmax_src(numLev_src)

    do jj = 1, nj_dest
      do ii = 1, ni_dest
        zpo(:) = press_dest(ii,jj,:)

        ! Interpolate in log of pressure or extrapolate with constant value
        ! for each destination pressure level

        ! Find the adjacent level below
        il(:) = 0
        do ji = 1, numLev_src
          do jo = 1, numLev_dest
            zrt = zpo(jo)
            zp = zpi(ji)
            xi = sign(1.0d0,iorder*(zrt-zp))
            il(jo) = il(jo) + max(0.0d0,xi)
          end do
        end do

        ! Interpolation/extrapolation
        do jo = 1, numLev_dest
          ik = il(jo)
          zp = zpo(jo)
          zp1 = zpi(ik)
          zp2 = zpi(ik+1)
          zw1 = log(zp/zp2)/log(zp1/zp2)
          zw2 = 1.d0 - zw1
          qmin_dest(ii,jj,jo) = zw1*zqmin_src(ik) +  zw2*zqmin_src(ik+1)
          qmax_dest(ii,jj,jo) = zw1*zqmax_src(ik) +  zw2*zqmax_src(ik+1)
        end do

      end do ! ii
    end do ! jj

  end subroutine lintv_minmax

end module humidityLimits_mod