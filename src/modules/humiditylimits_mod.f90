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

module humidityLimits_mod
  ! MODULE humidityLimits_mod (prefix='qlim' category='1. High-level functionality')
  !
  ! :Purpose: Various manipulations of humidity-related quantities.
  !
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

  !--------------------------------------------------------------------------
  ! qlim_gsvSaturationLimit
  !--------------------------------------------------------------------------
  subroutine qlim_gsvSaturationLimit(statevector)
    !
    !:Purpose: To impose saturation limit on humidity variable of a statevector
    !
    implicit none

    ! Arguments:
    type(struct_gsv) :: statevector

    write(*,*) 'qlim_gsvSaturationLimit: STARTING'

    if( .not. gsv_varExist(statevector,'HU') ) then
      if( mpi_myid == 0 ) write(*,*) 'qlim_gsvSaturationLimit: statevector does not ' // &
           'contain humidity ... doing nothing'
      return
    end if

    if ( statevector%dataKind == 8 ) then
      call qlim_gsvSaturationLimit_r8(statevector)
    else if ( statevector%dataKind == 4 ) then
      call qlim_gsvSaturationLimit_r4(statevector)
    else
      call utl_abort('qlim_gsvSaturationLimit: only compatible with single or double precision ' // &
           'data.')
    end if

  end subroutine qlim_gsvSaturationLimit

  !--------------------------------------------------------------------------
  ! qlim_gsvSaturationLimit_r8
  !--------------------------------------------------------------------------
  subroutine qlim_gsvSaturationLimit_r8(statevector)
    !
    !:Purpose: To impose saturation limit on humidity variable of a r8
    !          statevector
    !
    implicit none

    ! Arguments:
    type(struct_gsv) :: statevector

    ! Locals:
    type(struct_vco), pointer :: vco_ptr
    real(8), pointer :: lq_ptr(:,:,:,:), hu_ptr(:,:,:,:), tt_ptr(:,:,:,:), psfc_ptr(:,:,:,:)
    real(8), pointer :: pressure(:,:,:)
    real(8)          :: hu, husat, hu_modified, tt
    real(8), allocatable :: psfc(:,:)
    integer          :: lon1, lon2, lat1, lat2, lev1, lev2, ierr
    integer          :: lonIndex, latIndex, levIndex, stepIndex

    if ( statevector%dataKind /= 8 ) then
      call utl_abort('qlim_gsvSaturationLimit_r8: only compatible with double precision ' // &
           'data.')
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
            hu = hu_ptr(lonIndex,latIndex,levIndex,stepIndex)
            tt = tt_ptr(lonIndex,latIndex,levIndex,stepIndex)

            ! get the saturated vapor pressure from HU
            husat = foqst8(tt, pressure(lonIndex-lon1+1,latIndex-lat1+1,levIndex) )

            ! limit the humidity to the saturated humidity
            hu_modified = min(husat, hu)
            hu_ptr(lonIndex,latIndex,levIndex,stepIndex) = hu_modified

          end do ! lonIndex
        end do ! latIndex
      end do ! levIndex
      !$OMP END PARALLEL DO

      deallocate(pressure)

    end do ! stepIndex

  end subroutine qlim_gsvSaturationLimit_r8

  !--------------------------------------------------------------------------
  ! qlim_gsvSaturationLimit_r4
  !--------------------------------------------------------------------------
  subroutine qlim_gsvSaturationLimit_r4(statevector)
    !
    !:Purpose: To impose saturation limit on humidity variable of a r4
    !          statevector
    !
    implicit none

    ! Arguments:
    type(struct_gsv) :: statevector

    ! Locals:
    type(struct_vco), pointer :: vco_ptr
    real(4), pointer :: lq_ptr(:,:,:,:), hu_ptr(:,:,:,:), tt_ptr(:,:,:,:), psfc_ptr(:,:,:,:)
    real(8), pointer :: pressure(:,:,:)
    real(8)          :: hu, husat, hu_modified, tt
    real(8), allocatable :: psfc(:,:)
    integer          :: lon1, lon2, lat1, lat2, lev1, lev2, ierr
    integer          :: lonIndex, latIndex, levIndex, stepIndex

    if ( statevector%dataKind /= 4 ) then
      call utl_abort('qlim_gsvSaturationLimit_r4: only compatible with single precision ' // &
           'data.')
    end if

    vco_ptr => gsv_getVco(statevector)
    lq_ptr => gsv_getField_r4(statevector,'HU')
    hu_ptr => gsv_getField_r4(statevector,'HU')
    tt_ptr => gsv_getField_r4(statevector,'TT')

    lon1 = statevector%myLonBeg
    lon2 = statevector%myLonEnd
    lat1 = statevector%myLatBeg
    lat2 = statevector%myLatEnd
    lev1 = 1
    lev2 = vco_getNumLev(vco_ptr,'TH')

    allocate(psfc(lon2-lon1+1,lat2-lat1+1))
    do stepIndex = 1, statevector%numStep
      psfc_ptr => gsv_getField_r4(statevector,'P0')
      psfc(:,:) = real(psfc_ptr(:,:,1,stepIndex),8)
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
            hu = real(hu_ptr(lonIndex,latIndex,levIndex,stepIndex),8)
            tt = real(tt_ptr(lonIndex,latIndex,levIndex,stepIndex),8)

            ! get the saturated vapor pressure from HU
            husat = foqst8(tt, pressure(lonIndex-lon1+1,latIndex-lat1+1,levIndex) )

            ! limit the humidity to the saturated humidity
            hu_modified = min(husat, hu)
            hu_ptr(lonIndex,latIndex,levIndex,stepIndex) = real(hu_modified,4)

          end do ! lonIndex
        end do ! latIndex
      end do ! levIndex
      !$OMP END PARALLEL DO

      deallocate(pressure)

    end do ! stepIndex

  end subroutine qlim_gsvSaturationLimit_r4

  !--------------------------------------------------------------------------
  ! qlim_gsvRttovLimit
  !--------------------------------------------------------------------------
  subroutine qlim_gsvRttovLimit(statevector)
    !
    !:Purpose: To impose RTTOV limits on humidity
    !
    implicit none

    ! Arguments:
    type(struct_gsv) :: statevector

    ! Locals:
    integer          :: levIndex, numLev_rttov
    real(8), allocatable :: press_rttov(:), qmin_rttov(:), qmax_rttov(:)
    character(len=256) :: fileName
    integer :: fnom, fclos, ierr, nulfile

    write(*,*) 'qlim_gsvRttovLimit: STARTING'

    if ( .not. gsv_varExist(statevector,'HU') ) then
      if ( mpi_myid == 0 ) write(*,*) 'qlim_gsvRttovLimit: statevector does not ' // &
           'contain humidity ... doing nothing'
      return
    end if

    ! Read in RTTOV humidity limits
    fileName = "rttov_h2o_limits.dat"
    nulfile = 0
    ierr = fnom(nulfile, fileName, "FMT+OLD+R/O", 0)
    if( ierr /= 0 ) then
      if ( mpi_myid == 0 ) write(*,*) 'fileName = ', fileName
      call utl_abort('qlim_gsvRttovLimit: error opening the humidity limits file')
    end if

    read(nulfile,*) numLev_rttov
    if ( mpi_myid == 0 ) write(*,*) 'qlim_gsvRttovLimit: rttov number of levels = ', numLev_rttov
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
      if ( mpi_myid == 0 ) write(*,fmt='(" qlim_gsvRttovLimit:   LEVEL = ",I4,", PRES = ",F9.0,", HUMIN = ",E10.2,", HUMAX = ",E10.2)') &
           levIndex, press_rttov(levIndex), qmin_rttov(levIndex), qmax_rttov(levIndex)
    end do

    if ( statevector%dataKind == 8 ) then
      call qlim_gsvRttovLimit_r8(statevector, press_rttov, qmin_rttov, qmax_rttov, numLev_rttov)
    else if ( statevector%dataKind == 4 ) then
      call qlim_gsvRttovLimit_r4(statevector, press_rttov, qmin_rttov, qmax_rttov, numLev_rttov)
    else
      call utl_abort('qlim_gsvSaturationLimit: only compatible with single or double precision ' // &
           'data.')
    end if

    deallocate(qmax_rttov)
    deallocate(qmin_rttov)
    deallocate(press_rttov)

  end subroutine qlim_gsvRttovLimit

  !--------------------------------------------------------------------------
  ! qlim_gsvRttovLimit_r8
  !--------------------------------------------------------------------------
  subroutine qlim_gsvRttovLimit_r8(statevector, press_rttov, qmin_rttov, &
                                   qmax_rttov, numLev_rttov)
    !
    !:Purpose: To impose RTTOV limits on humidity
    !
    implicit none

    ! Arguments:
    type(struct_gsv) :: statevector
    integer, intent(in) :: numLev_rttov
    real(8), intent(in) :: press_rttov(numLev_rttov), qmin_rttov(numLev_rttov), qmax_rttov(numLev_rttov)

    ! Locals:
    type(struct_vco), pointer :: vco_ptr
    real(8), pointer :: lq_ptr(:,:,:,:), hu_ptr(:,:,:,:), psfc_ptr(:,:,:,:)
    real(8), pointer :: pressure(:,:,:)
    real(8)          :: hu, hu_modified
    real(8), allocatable :: psfc(:,:)
    integer          :: lon1, lon2, lat1, lat2, lev1, lev2
    integer          :: lonIndex, latIndex, levIndex, stepIndex
    integer          :: ni, nj, numLev, ierr
    real(8), allocatable :: qmin3D_rttov(:,:,:), qmax3D_rttov(:,:,:)

    if( statevector%dataKind /= 8 ) then
      call utl_abort('qlim_gsvRttovLimit_r8: only compatible with double precision ' // &
           'data.')
    end if

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
      call qlim_lintv_minmax(press_rttov, qmin_rttov, qmax_rttov, numLev_rttov, &
           ni, nj, numLev, pressure, qmin3D_rttov, qmax3D_rttov)

      !$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex, hu, hu_modified)
      do levIndex = lev1, lev2
        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2
            hu = hu_ptr(lonIndex,latIndex,levIndex,stepIndex)

            ! limit the humidity according to the rttov limits
            hu_modified = max(hu, qmin3D_rttov(lonIndex - lon1 + 1, latIndex - lat1 + 1, levIndex) )
            hu_modified = min(hu_modified, qmax3D_rttov(lonIndex - lon1 + 1, latIndex - lat1 + 1, levIndex) )
            hu_ptr(lonIndex,latIndex,levIndex,stepIndex) = hu_modified

          end do ! lonIndex
        end do ! latIndex
      end do ! levIndex
      !$OMP END PARALLEL DO

      deallocate(pressure)

    end do ! stepIndex

    deallocate( psfc )
    deallocate( qmin3D_rttov )
    deallocate( qmax3D_rttov )

  end subroutine qlim_gsvRttovLimit_r8

  !--------------------------------------------------------------------------
  ! qlim_gsvRttovLimit_r4
  !--------------------------------------------------------------------------
  subroutine qlim_gsvRttovLimit_r4(statevector, press_rttov, qmin_rttov, &
                                   qmax_rttov, numLev_rttov)
    !
    !:Purpose: To impose RTTOV limits on humidity
    !
    implicit none

    ! Arguments:
    type(struct_gsv) :: statevector
    integer, intent(in) :: numLev_rttov
    real(8), intent(in) :: press_rttov(numLev_rttov), qmin_rttov(numLev_rttov), qmax_rttov(numLev_rttov)

    ! Locals:
    type(struct_vco), pointer :: vco_ptr
    real(4), pointer :: lq_ptr(:,:,:,:), hu_ptr(:,:,:,:), psfc_ptr(:,:,:,:)
    real(8), pointer :: pressure(:,:,:)
    real(8)          :: hu, hu_modified
    real(8), allocatable :: psfc(:,:)
    integer          :: lon1, lon2, lat1, lat2, lev1, lev2
    integer          :: lonIndex, latIndex, levIndex, stepIndex
    integer          :: ni, nj, numLev, ierr
    real(8), allocatable :: qmin3D_rttov(:,:,:), qmax3D_rttov(:,:,:)

    if( statevector%dataKind /= 4 ) then
      call utl_abort('qlim_gsvRttovLimit_r4: only compatible with single precision ' // &
           'data.')
    end if

    vco_ptr => gsv_getVco(statevector)
    lq_ptr => gsv_getField_r4(statevector,'HU')
    hu_ptr => gsv_getField_r4(statevector,'HU')

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
      psfc_ptr => gsv_getField_r4(statevector,'P0')
      psfc(:,:) = real(psfc_ptr(:,:,1,stepIndex),8)
      nullify(pressure)
      ierr = vgd_levels(vco_ptr%vgrid,           &
           ip1_list=vco_ptr%ip1_T,  &
           levels=pressure,         &
           sfc_field=psfc, &
           in_log=.false.)

      ! Interpolate RTTOV limits onto model levels
      call qlim_lintv_minmax(press_rttov, qmin_rttov, qmax_rttov, numLev_rttov, &
           ni, nj, numLev, pressure, qmin3D_rttov, qmax3D_rttov)

      !$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex, hu, hu_modified)
      do levIndex = lev1, lev2
        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2
            hu = real(hu_ptr(lonIndex,latIndex,levIndex,stepIndex),8)

            ! limit the humidity according to the rttov limits
            hu_modified = max(hu, qmin3D_rttov(lonIndex - lon1 + 1, latIndex - lat1 + 1, levIndex) )
            hu_modified = min(hu_modified, qmax3D_rttov(lonIndex - lon1 + 1, latIndex - lat1 + 1, levIndex) )
            hu_ptr(lonIndex,latIndex,levIndex,stepIndex) = real(hu_modified,4)

          end do ! lonIndex
        end do ! latIndex
      end do ! levIndex
      !$OMP END PARALLEL DO

      deallocate(pressure)

    end do ! stepIndex

    deallocate( psfc )
    deallocate( qmin3D_rttov )
    deallocate( qmax3D_rttov )

  end subroutine qlim_gsvRttovLimit_r4

  !--------------------------------------------------------------------------
  ! qlim_lintv_minmax
  !--------------------------------------------------------------------------
  subroutine qlim_lintv_minmax(press_src, qmin_src, qmax_src, numLev_src, &
                               ni_dest, nj_dest, numLev_dest, press_dest, &
                               qmin_dest, qmax_dest)
    !
    !:Purpose: To perform the vertical interpolation in log of pressure and
    !          and constant value extrapolation of one-dimensional vectors.
    implicit none

    ! Arguments:
    real(8), intent(in) :: press_src(numLev_src) ! Vertical levels, pressure (source)
    real(8), intent(in) :: qmin_src(numLev_src)  ! Vectors to be interpolated (source)
    real(8), intent(in) :: qmax_src(numLev_src)  ! Vectors to be interpolated (source)
    integer, intent(in) :: numLev_src ! Number of input levels (source)
    integer, intent(in) :: ni_dest ! Number of profiles
    integer, intent(in) :: nj_dest ! Number of profiles
    integer, intent(in) :: numLev_dest ! Number of output levels (destination)
    real(8)             :: press_dest(:,:,:) ! Vertical levels, pressure (destination)
    real(8)             :: qmin_dest(:,:,:)  ! Interpolated profiles (destination)
    real(8)             :: qmax_dest(:,:,:)  ! Interpolated profiles (destination)

    ! Locals:
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

  end subroutine qlim_lintv_minmax

end module humidityLimits_mod
