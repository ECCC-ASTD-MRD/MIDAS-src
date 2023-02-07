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
  ! MODULE humidityLimits_mod (prefix='qlim' category='4. Data Object transformations')
  !
  ! :Purpose: Various manipulations of humidity-related quantities.
  !
  use midasMpi_mod
  use utilities_mod
  use mathPhysConstants_mod
  use varNameList_mod
  use physicsFunctions_mod
  use verticalCoord_mod
  use gridStateVector_mod
  use ensembleStateVector_mod
  use calcHeightAndPressure_mod
  implicit none
  save
  private

  ! public procedures
  public :: qlim_saturationLimit, qlim_rttovLimit, qlim_setMin
  public :: qlim_getMinValueCloud, qlim_getMaxValueCloud

  real(8), parameter :: mixratio_to_ppmv = 1.60771704d+6
  real(8) :: qlim_minValueLWCR, qlim_minValueIWCR, qlim_minValueRF, qlim_minValueSF, qlim_minValueCLDR
  real(8) :: qlim_maxValueLWCR, qlim_maxValueIWCR, qlim_maxValueRF, qlim_maxValueSF, qlim_maxValueCLDR

  ! interface for qlim_saturationLimit
  interface qlim_saturationLimit
    module procedure qlim_saturationLimit_gsv
    module procedure qlim_saturationLimit_ens
  end interface qlim_saturationLimit
  
  ! interface for qlim_rttovLimit
  interface qlim_rttovLimit
    module procedure qlim_rttovLimit_gsv
    module procedure qlim_rttovLimit_ens
  end interface qlim_rttovLimit

  ! interface for qlim_setMin
  interface qlim_setMin
    module procedure qlim_setMin_ens
  end interface qlim_setMin
  
contains

  !--------------------------------------------------------------------------
  ! readNameList
  !--------------------------------------------------------------------------
  subroutine readNameList
    !
    ! :Purpose: Reading NAMQLIM namelist by any subroutines in humidityLimits_mod module.
    !
    implicit none

    ! Locals:
    integer :: nulnam, ierr
    integer, external :: fnom, fclos
    logical, save :: nmlAlreadyRead = .false.

    ! Namelist variables:
    real(8) :: minValueLWCR ! minimum LWCR value
    real(8) :: maxValueLWCR ! maximum LWCR value
    real(8) :: minValueIWCR ! minimum IWCR value
    real(8) :: maxValueIWCR ! maximum IWCR value
    real(8) :: minValueRF   ! minimum   RF value
    real(8) :: maxValueRF   ! maximum   RF value
    real(8) :: minValueSF   ! minimum   SF value
    real(8) :: maxValueSF   ! maximum   SF value
    real(8) :: minValueCLDR ! minimum CLDR value
    real(8) :: maxValueCLDR ! maximum CLDR value

    NAMELIST /NAMQLIM/ minValueLWCR, maxValueLWCR, minValueIWCR, maxValueIWCR, &
                       minValueRF, maxValueRF, minValueSF, maxValueSF, minValueCLDR, maxValueCLDR

    if ( nmlAlreadyRead ) return

    nmlAlreadyRead = .true.

    !- Setting default values
    minValueLWCR = 1.0d-9
    maxValueLWCR = 1.0d0

    minValueIWCR = 1.0d-9
    maxValueIWCR = 1.0d0
    
    minValueRF = 0.0d0
    maxValueRF = 1.0d0
    
    minValueSF = 0.0d0
    maxValueSF = 1.0d0

    minValueCLDR = 0.0d0
    maxValueCLDR = 1.0d0

    if ( .not. utl_isNamelistPresent('NAMQLIM','./flnml') ) then
      if ( mmpi_myid == 0 ) then
        write(*,*) 'NAMQLIM is missing in the namelist. The default values will be taken.'
      end if

    else
      ! Reading the namelist
      nulnam = 0
      ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
      read(nulnam, nml=namqlim, iostat=ierr)
      if ( ierr /= 0) call utl_abort('humidityLimits_mod: Error reading namelist')
      ierr = fclos(nulnam)

    end if
    if ( mmpi_myid == 0 ) write(*,nml=namqlim)

    ! Transfer namelist variables to module variables.
    qlim_minValueLWCR = minValueLWCR
    qlim_maxValueLWCR = maxValueLWCR

    qlim_minValueIWCR = minValueIWCR
    qlim_maxValueIWCR = maxValueIWCR
    
    qlim_minValueRF   = minValueRF
    qlim_maxValueRF   = maxValueRF
    
    qlim_minValueSF   = minValueSF
    qlim_maxValueSF   = maxValueSF

    qlim_minValueCLDR = minValueCLDR
    qlim_maxValueCLDR = maxValueCLDR

  end subroutine readNameList

  !--------------------------------------------------------------------------
  ! qlim_saturationLimit_gsv
  !--------------------------------------------------------------------------
  subroutine qlim_saturationLimit_gsv(statevector)
    !
    !:Purpose: To impose saturation limit on humidity variable of a statevector
    !
    implicit none

    ! Arguments:
    type(struct_gsv) :: statevector

    ! Locals:
    type(struct_vco), pointer :: vco_ptr
    real(4), pointer :: hu_ptr_r4(:,:,:,:), tt_ptr_r4(:,:,:,:), psfc_ptr_r4(:,:,:,:)
    real(8), pointer :: hu_ptr_r8(:,:,:,:), tt_ptr_r8(:,:,:,:), psfc_ptr_r8(:,:,:,:)
    real(8), pointer :: pressure(:,:,:)
    real(8)          :: hu, husat, hu_modified, tt
    real(8), allocatable :: psfc(:,:)
    integer          :: lon1, lon2, lat1, lat2, lev1, lev2, ierr
    integer          :: lonIndex, latIndex, levIndex, stepIndex

    if (mmpi_myid == 0) write(*,*) 'qlim_saturationLimit_gsv: STARTING'

    if( .not. gsv_varExist(statevector,'HU') ) then
      if( mmpi_myid == 0 ) write(*,*) 'qlim_saturationLimit_gsv: statevector does not ' // &
           'contain humidity ... doing nothing'
      return
    end if

    vco_ptr => gsv_getVco(statevector)
    if (stateVector%dataKind == 8) then
      call gsv_getField(statevector,hu_ptr_r8,'HU')
      call gsv_getField(statevector,tt_ptr_r8,'TT')
    else
      call gsv_getField(statevector,hu_ptr_r4,'HU')
      call gsv_getField(statevector,tt_ptr_r4,'TT')
    end if

    lon1 = statevector%myLonBeg
    lon2 = statevector%myLonEnd
    lat1 = statevector%myLatBeg
    lat2 = statevector%myLatEnd
    lev1 = 1
    lev2 = gsv_getNumLev(statevector,'TH')

    allocate(psfc(lon2-lon1+1,lat2-lat1+1))
    do stepIndex = 1, statevector%numStep
      if (stateVector%dataKind == 8) then
        call gsv_getField(statevector,psfc_ptr_r8,'P0')
        psfc(:,:) = psfc_ptr_r8(:,:,1,stepIndex)
      else
        call gsv_getField(statevector,psfc_ptr_r4,'P0')
        psfc(:,:) = psfc_ptr_r4(:,:,1,stepIndex)
      end if
      call czp_fetch3DField(vco_ptr, psfc, fldT_opt=pressure)
      if (stateVector%dataKind == 8) then
        !$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex, hu, tt, husat, hu_modified)
        do levIndex = lev1, lev2
          do latIndex = lat1, lat2
            do lonIndex = lon1, lon2
              hu = hu_ptr_r8(lonIndex,latIndex,levIndex,stepIndex)
              tt = tt_ptr_r8(lonIndex,latIndex,levIndex,stepIndex)

              ! get the saturated vapor pressure from HU
              husat = phf_foqst8(tt, pressure(lonIndex-lon1+1,latIndex-lat1+1,levIndex) )

              ! limit the humidity to the saturated humidity
              hu_modified = min(husat, hu)
              hu_ptr_r8(lonIndex,latIndex,levIndex,stepIndex) = hu_modified

            end do ! lonIndex
          end do ! latIndex
        end do ! levIndex
        !$OMP END PARALLEL DO
      else
        !$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex, hu, tt, husat, hu_modified)
        do levIndex = lev1, lev2
          do latIndex = lat1, lat2
            do lonIndex = lon1, lon2
              hu = hu_ptr_r4(lonIndex,latIndex,levIndex,stepIndex)
              tt = tt_ptr_r4(lonIndex,latIndex,levIndex,stepIndex)

              ! get the saturated vapor pressure from HU
              husat = phf_foqst8(tt, pressure(lonIndex-lon1+1,latIndex-lat1+1,levIndex) )

              ! limit the humidity to the saturated humidity
              hu_modified = min(husat, hu)
              hu_ptr_r4(lonIndex,latIndex,levIndex,stepIndex) = hu_modified

            end do ! lonIndex
          end do ! latIndex
        end do ! levIndex
        !$OMP END PARALLEL DO
      end if

      deallocate(pressure)

    end do ! stepIndex

    deallocate(psfc)

  end subroutine qlim_saturationLimit_gsv

  !--------------------------------------------------------------------------
  ! qlim_saturationLimit_ens
  !--------------------------------------------------------------------------
  subroutine qlim_saturationLimit_ens(ensemble)
    !
    !:Purpose: To impose saturation limit on humidity variable of an ensemble
    !
    implicit none

    ! Arguments:
    type(struct_ens) :: ensemble

    ! Locals:
    type(struct_vco), pointer :: vco_ptr
    real(4), pointer :: hu_ptr_r4(:,:,:,:), tt_ptr_r4(:,:,:,:), psfc_ptr_r4(:,:,:,:)
    real(8), pointer :: pressure(:,:,:)
    real(8)          :: hu, husat, hu_modified, tt
    real(8), allocatable :: psfc(:,:)
    integer          :: lon1, lon2, lat1, lat2, numLev, ierr
    integer          :: lonIndex, latIndex, levIndex, stepIndex, memberIndex, varLevIndex

    if (mmpi_myid == 0) write(*,*) 'qlim_saturationLimit_ens: STARTING'

    if (ens_getDataKind(ensemble) == 8) then
      call utl_abort('qlim_saturationLimit_ens: Not compatible with dataKind = 8')
    end if

    if( .not. ens_varExist(ensemble,'HU') ) then
      if( mmpi_myid == 0 ) write(*,*) 'qlim_saturationLimit_ens: ensemble does not ' // &
           'contain humidity ... doing nothing'
      return
    end if

    vco_ptr => ens_getVco(ensemble)
    numLev = ens_getNumLev(ensemble,'TH')
    call ens_getLatLonBounds(ensemble, lon1, lon2, lat1, lat2)
    allocate(psfc(ens_getNumMembers(ensemble),ens_getNumStep(ensemble)))

    do latIndex = lat1, lat2
      do lonIndex = lon1, lon2

        ! compute pressure for all members and steps
        varLevIndex = ens_getKFromLevVarName(ensemble, 1, 'P0')
        psfc_ptr_r4 => ens_getOneLev_r4(ensemble,varLevIndex)
        psfc(:,:) = psfc_ptr_r4(:,:,lonIndex,latIndex)
        call czp_fetch3DField(vco_ptr, psfc, fldT_opt=pressure)

        do levIndex = 1, numLev
          varLevIndex = ens_getKFromLevVarName(ensemble, levIndex, 'HU')
          hu_ptr_r4 => ens_getOneLev_r4(ensemble,varLevIndex)
          varLevIndex = ens_getKFromLevVarName(ensemble, levIndex, 'TT')
          tt_ptr_r4 => ens_getOneLev_r4(ensemble,varLevIndex)

          !$OMP PARALLEL DO PRIVATE (stepIndex, memberIndex, hu, tt, husat, hu_modified)
          do stepIndex = 1, ens_getNumStep(ensemble)
            do memberIndex = 1, ens_getNumMembers(ensemble)
              hu = hu_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex)
              tt = tt_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex)

              ! get the saturated vapor pressure from HU
              husat = phf_foqst8(tt, pressure(memberIndex,stepIndex,levIndex) )

              ! limit the humidity to the saturated humidity
              hu_modified = min(husat, hu)
              hu_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) = hu_modified

            end do ! memberIndex
          end do ! stepIndex
          !$OMP END PARALLEL DO

        end do ! levIndex
        deallocate(pressure)

      end do ! lonIndex
    end do ! latIndex

    deallocate(psfc)

  end subroutine qlim_saturationLimit_ens

  !--------------------------------------------------------------------------
  ! qlim_rttovLimit_gsv
  !--------------------------------------------------------------------------
  subroutine qlim_rttovLimit_gsv(statevector, varName_opt, applyLimitToCloud_opt)
    !
    !:Purpose: To impose RTTOV limits on humidity/cloud
    !
    implicit none

    ! Arguments:
    type(struct_gsv),        intent(inout) :: statevector
    character(len=*), optional, intent(in) :: varName_opt
    logical,          optional, intent(in) :: applyLimitToCloud_opt

    ! Locals:
    type(struct_vco), pointer :: vco_ptr
    real(8), allocatable :: press_rttov(:), qmin_rttov(:), qmax_rttov(:)
    real(8), allocatable :: psfc(:,:)
    real(8), allocatable :: qmin3D_rttov(:,:,:), qmax3D_rttov(:,:,:)
    real(8), pointer     :: hu_ptr_r8(:,:,:,:), psfc_ptr_r8(:,:,:,:)
    real(8), pointer     :: cld_ptr_r8(:,:,:,:)
    real(4), pointer     :: hu_ptr_r4(:,:,:,:), psfc_ptr_r4(:,:,:,:)
    real(4), pointer     :: cld_ptr_r4(:,:,:,:)
    real(8), pointer     :: pressure(:,:,:)
    real(8)              :: hu, hu_modified
    real(8)              :: cld, cld_modified
    real(8)              :: minValueCld, maxValueCld
    integer              :: lon1, lon2, lat1, lat2, lev1, lev2, varNameIndex
    integer              :: lonIndex, latIndex, levIndex, stepIndex
    integer              :: ni, nj, numLev, numLev_rttov
    integer              :: fnom, fclos, ierr, nulfile
    character(len=256)   :: fileName
    character(len=4)     :: varName
    logical, save        :: firstTime=.true.
    logical              :: applyLimitToAllVarname
    logical              :: applyLimitToCloud

    if (mmpi_myid == 0) write(*,*) 'qlim_rttovLimit_gsv: STARTING'

    if (present(varName_opt)) then
      applyLimitToAllVarname = .false.
      varName = varName_opt
    else
      applyLimitToAllVarname = .true.
      varName = 'XXXX'  
    end if

    ! for cloud, limits are applied to ALL cloud variables.
    if (present(applyLimitToCloud_opt)) then
      applyLimitToCloud = applyLimitToCloud_opt
      if (.not. applyLimitToCloud) call utl_abort('qlim_rttovLimit_gsv: remove applyLimitToCloud_opt argument')
    else
      if (vnl_isCloudVar(varName)) then
        if (mmpi_myid == 0) write(*,*) 'qlim_rttovLimit_gsv: limits are applied to ALL cloud variables'
        applyLimitToCloud = .true.
      else
        applyLimitToCloud = .false.
      end if
    end if

    if (applyLimitToCloud) applyLimitToAllVarname = .false.

    if ((applyLimitToAllVarname .or. trim(varName) == 'HU') .and. &
         gsv_varExist(statevector,'HU')) then

      if (mmpi_myid == 0) write(*,*) 'qlim_rttovLimit_gsv: applying limits to HU.'

      ! Read in RTTOV humidity limits
      fileName = "rttov_h2o_limits.dat"
      nulfile = 0
      ierr = fnom(nulfile, fileName, "FMT+OLD+R/O", 0)
      if( ierr /= 0 ) then
        if ( mmpi_myid == 0 ) write(*,*) 'fileName = ', fileName
        call utl_abort('qlim_rttovLimit_gsv: error opening the humidity limits file')
      end if

      read(nulfile,*) numLev_rttov
      if ( mmpi_myid == 0 .and. firstTime ) write(*,*) 'qlim_rttovLimit_gsv: rttov number of levels = ', numLev_rttov
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

      if (firstTime) then
        write(*,*) ' '
        do levIndex = 1, numLev_rttov
          if ( mmpi_myid == 0 ) write(*,fmt='(" qlim_rttovLimit_gsv:   LEVEL = ",I4,", PRES = ",F9.0,", HUMIN = ",E10.2,", HUMAX = ",E10.2)') &
              levIndex, press_rttov(levIndex), qmin_rttov(levIndex), qmax_rttov(levIndex)
        end do
        firstTime = .false.
      end if

      vco_ptr => gsv_getVco(statevector)
      if (statevector%dataKind == 8) then
        call gsv_getField(statevector,hu_ptr_r8,'HU')
      else
        call gsv_getField(statevector,hu_ptr_r4,'HU')
      end if

      lon1 = statevector%myLonBeg
      lon2 = statevector%myLonEnd
      lat1 = statevector%myLatBeg
      lat2 = statevector%myLatEnd
      lev1 = 1
      lev2 = gsv_getNumLev(statevector,'TH')

      ni = lon2 - lon1 + 1
      nj = lat2 - lat1 + 1
      numLev = lev2 - lev1 + 1
      allocate( qmin3D_rttov(ni,nj,numLev) )
      allocate( qmax3D_rttov(ni,nj,numLev) )
      allocate( psfc(lon2-lon1+1,lat2-lat1+1) )

      do stepIndex = 1, statevector%numStep
        if (statevector%dataKind == 8) then
          call gsv_getField(statevector,psfc_ptr_r8,'P0')
          psfc(:,:) = psfc_ptr_r8(:,:,1,stepIndex)
        else
          call gsv_getField(statevector,psfc_ptr_r4,'P0')
          psfc(:,:) = real(psfc_ptr_r4(:,:,1,stepIndex),8)
        end if
        call czp_fetch3DField(vco_ptr, psfc, fldT_opt=pressure)

        ! Interpolate RTTOV limits onto model levels
        call qlim_lintv_minmax(press_rttov, qmin_rttov, qmax_rttov, numLev_rttov, &
            ni, nj, numLev, pressure, qmin3D_rttov, qmax3D_rttov)

        !$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex, hu, hu_modified)
        do levIndex = lev1, lev2
          do latIndex = lat1, lat2
            do lonIndex = lon1, lon2
              if (statevector%dataKind == 8) then
                hu = hu_ptr_r8(lonIndex,latIndex,levIndex,stepIndex)
              else
                hu = real(hu_ptr_r4(lonIndex,latIndex,levIndex,stepIndex),8)
              end if

              ! limit the humidity according to the rttov limits
              hu_modified = max(hu, qmin3D_rttov(lonIndex - lon1 + 1, latIndex - lat1 + 1, levIndex) )
              hu_modified = min(hu_modified, qmax3D_rttov(lonIndex - lon1 + 1, latIndex - lat1 + 1, levIndex) )
              if (statevector%dataKind == 8) then
                hu_ptr_r8(lonIndex,latIndex,levIndex,stepIndex) = hu_modified
              else
                hu_ptr_r4(lonIndex,latIndex,levIndex,stepIndex) = real(hu_modified,4)
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

      deallocate(qmax_rttov)
      deallocate(qmin_rttov)
      deallocate(press_rttov)

    end if

    ! apply limits to ALL available cloud variables
    if ((applyLimitToAllVarname .or. applyLimitToCloud) .and. cloudExistInStateVector(statevector)) then

      do varNameIndex = 1, vnl_numvarmaxCloud
        if (.not. gsv_varExist(statevector, vnl_varNameListCloud(varNameIndex))) cycle

        if (mmpi_myid == 0) write(*,*) 'qlim_rttovLimit_gsv: applying limits to ', &
                                        vnl_varNameListCloud(varNameIndex)

        if (statevector%dataKind == 8) then
          call gsv_getField(statevector, cld_ptr_r8, vnl_varNameListCloud(varNameIndex))
        else
          call gsv_getField(statevector, cld_ptr_r4, vnl_varNameListCloud(varNameIndex))
        end if

        minValueCld = qlim_getMinValueCloud(vnl_varNameListCloud(varNameIndex))
        maxValueCld = qlim_getMaxValueCloud(vnl_varNameListCloud(varNameIndex))

        lon1 = statevector%myLonBeg
        lon2 = statevector%myLonEnd
        lat1 = statevector%myLatBeg
        lat2 = statevector%myLatEnd
        lev1 = 1
        lev2 = gsv_getNumLev(statevector,'TH')

        do stepIndex = 1, statevector%numStep

          !$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex, cld, cld_modified)
          do levIndex = lev1, lev2
            do latIndex = lat1, lat2
              do lonIndex = lon1, lon2
                if (statevector%dataKind == 8) then
                  cld = cld_ptr_r8(lonIndex,latIndex,levIndex,stepIndex)
                else
                  cld = real(cld_ptr_r4(lonIndex,latIndex,levIndex,stepIndex),8)
                end if

                cld_modified = max(cld,minValueCld)
                cld_modified = min(cld_modified,maxValueCld)
                if (statevector%dataKind == 8) then
                  cld_ptr_r8(lonIndex,latIndex,levIndex,stepIndex) = cld_modified
                else
                  cld_ptr_r4(lonIndex,latIndex,levIndex,stepIndex) = real(cld_modified,4)
                end if

              end do ! lonIndex
            end do ! latIndex
          end do ! levIndex
          !$OMP END PARALLEL DO

        end do ! stepIndex
      end do ! varNameIndex

    end if

  end subroutine qlim_rttovLimit_gsv

  !--------------------------------------------------------------------------
  ! qlim_rttovLimit_ens
  !--------------------------------------------------------------------------
  subroutine qlim_rttovLimit_ens(ensemble, varName_opt, applyLimitToCloud_opt)
    !
    !:Purpose: To impose RTTOV limits on humidity/cloud
    !
    implicit none

    ! Arguments:
    type(struct_ens),        intent(inout) :: ensemble
    character(len=*), optional, intent(in) :: varName_opt
    logical,          optional, intent(in) :: applyLimitToCloud_opt

    ! Locals:
    type(struct_vco), pointer :: vco_ptr
    real(8), allocatable :: press_rttov(:), qmin_rttov(:), qmax_rttov(:)
    real(8), allocatable :: psfc(:,:)
    real(8), allocatable :: qmin3D_rttov(:,:,:), qmax3D_rttov(:,:,:)
    real(4), pointer     :: hu_ptr_r4(:,:,:,:), psfc_ptr_r4(:,:,:,:)
    real(4), pointer     :: cld_ptr_r4(:,:,:,:)
    real(8), pointer     :: pressure(:,:,:)
    real(8)              :: hu, hu_modified
    real(8)              :: cld, cld_modified
    real(8)              :: minValueCld, maxValueCld
    integer              :: lon1, lon2, lat1, lat2, varNameIndex
    integer              :: lonIndex, latIndex, levIndex, stepIndex, varLevIndex, memberIndex
    integer              :: numMember, numStep, numLev, numLev_rttov
    integer              :: fnom, fclos, ierr, nulfile
    character(len=256)   :: fileName
    character(len=4)     :: varName
    logical, save        :: firstTime=.true.
    logical              :: applyLimitToAllVarname
    logical              :: applyLimitToCloud

    if (mmpi_myid == 0) write(*,*) 'qlim_rttovLimit_ens: STARTING'

    if (present(varName_opt)) then
      applyLimitToAllVarname = .false.
      varName = varName_opt
    else
      applyLimitToAllVarname = .true.
      varName = 'XXXX'    
    end if

    ! for cloud, limits are applied to ALL cloud variables.
    if (present(applyLimitToCloud_opt)) then
      applyLimitToCloud = applyLimitToCloud_opt
      if (.not. applyLimitToCloud) call utl_abort('qlim_rttovLimit_ens: remove applyLimitToCloud_opt argument')
    else
      if (vnl_isCloudVar(varName)) then
        if (mmpi_myid == 0) write(*,*) 'qlim_rttovLimit_ens: limits are applied to ALL cloud variables'
        applyLimitToCloud = .true.
      else
        applyLimitToCloud = .false.
      end if
    end if

    if (applyLimitToCloud) applyLimitToAllVarname = .false.

    if ((applyLimitToAllVarname .or. trim(varName) == 'HU') .and. &
         ens_varExist(ensemble,'HU')) then

      if ( mmpi_myid == 0 ) write(*,*) 'qlim_rttovLimit_ens:  applying limits to HU.'

      if (ens_getDataKind(ensemble) == 8) then
        call utl_abort('qlim_rttovLimit_ens: Not compatible with dataKind = 8')
      end if

      ! Read in RTTOV humidity limits
      fileName = "rttov_h2o_limits.dat"
      nulfile = 0
      ierr = fnom(nulfile, fileName, "FMT+OLD+R/O", 0)
      if( ierr /= 0 ) then
        if ( mmpi_myid == 0 ) write(*,*) 'fileName = ', fileName
        call utl_abort('qlim_rttovLimit_ens: error opening the humidity limits file')
      end if

      read(nulfile,*) numLev_rttov
      if ( mmpi_myid == 0 .and. firstTime ) write(*,*) 'qlim_rttovLimit_ens: rttov number of levels = ', numLev_rttov
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

      if (firstTime) then
        write(*,*) ' '
        do levIndex = 1, numLev_rttov
          if ( mmpi_myid == 0 ) write(*,fmt='(" qlim_rttovLimit_ens:   LEVEL = ",I4,", PRES = ",F9.0,", HUMIN = ",E10.2,", HUMAX = ",E10.2)') &
              levIndex, press_rttov(levIndex), qmin_rttov(levIndex), qmax_rttov(levIndex)
        end do
        firstTime = .false.
      end if

      vco_ptr => ens_getVco(ensemble)
      numLev = ens_getNumLev(ensemble,'TH')
      numMember = ens_getNumMembers(ensemble)
      numStep = ens_getNumStep(ensemble)
      call ens_getLatLonBounds(ensemble, lon1, lon2, lat1, lat2)

      allocate( psfc(numMember,numStep) )
      allocate( qmin3D_rttov(numMember,numStep,numLev) )
      allocate( qmax3D_rttov(numMember,numStep,numLev) )

      do latIndex = lat1, lat2
        do lonIndex = lon1, lon2

          varLevIndex = ens_getKFromLevVarName(ensemble, 1, 'P0')
          psfc_ptr_r4 => ens_getOneLev_r4(ensemble,varLevIndex)
          psfc(:,:) = real(psfc_ptr_r4(:,:,lonIndex,latIndex),8)
          call czp_fetch3DField(vco_ptr, psfc, fldT_opt=pressure)

          ! Interpolate RTTOV limits onto model levels
          call qlim_lintv_minmax(press_rttov, qmin_rttov, qmax_rttov, numLev_rttov, &
                                numMember, numStep, numLev, pressure,  &
                                qmin3D_rttov, qmax3D_rttov)

          do levIndex = 1, numLev

            varLevIndex = ens_getKFromLevVarName(ensemble, levIndex, 'HU')
            hu_ptr_r4 => ens_getOneLev_r4(ensemble,varLevIndex)

            !$OMP PARALLEL DO PRIVATE (stepIndex, memberIndex, hu, hu_modified)
            do stepIndex = 1, numStep
              do memberIndex = 1, numMember

                hu = real(hu_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex),8)

                ! limit the humidity according to the rttov limits
                hu_modified = max(hu, qmin3D_rttov(memberIndex, stepIndex, levIndex) )
                hu_modified = min(hu_modified, qmax3D_rttov(memberIndex, stepIndex, levIndex) )
                hu_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) = real(hu_modified,4)

              end do ! memberIndex
            end do ! stepIndex
            !$OMP END PARALLEL DO
          end do ! levIndex

          deallocate(pressure)

        end do ! lonIndex
      end do ! latIndex

      deallocate( psfc )
      deallocate( qmin3D_rttov )
      deallocate( qmax3D_rttov )

      deallocate(qmax_rttov)
      deallocate(qmin_rttov)
      deallocate(press_rttov)
    end if

    ! apply limits to ALL available cloud variables
    if ((applyLimitToAllVarname .or. applyLimitToCloud) .and. cloudExistInEnsemble(ensemble)) then

      vco_ptr => ens_getVco(ensemble)
      numLev = ens_getNumLev(ensemble,'TH')
      numMember = ens_getNumMembers(ensemble)
      numStep = ens_getNumStep(ensemble)
      call ens_getLatLonBounds(ensemble, lon1, lon2, lat1, lat2)

      do varNameIndex = 1, vnl_numvarmaxCloud
        if (.not. ens_varExist(ensemble, vnl_varNameListCloud(varNameIndex))) cycle

        if (mmpi_myid == 0) write(*,*) 'qlim_rttovLimit_ens:  applying limits to ', &
                                        vnl_varNameListCloud(varNameIndex)

        minValueCld = qlim_getMinValueCloud(vnl_varNameListCloud(varNameIndex))
        maxValueCld = qlim_getMaxValueCloud(vnl_varNameListCloud(varNameIndex))

        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2

            do levIndex = 1, numLev
                varLevIndex = ens_getKFromLevVarName(ensemble, levIndex, vnl_varNameListCloud(varNameIndex))
                cld_ptr_r4 => ens_getOneLev_r4(ensemble,varLevIndex)

                !$OMP PARALLEL DO PRIVATE (stepIndex, memberIndex, cld, cld_modified)
                do stepIndex = 1, numStep
                  do memberIndex = 1, numMember

                      cld = real(cld_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex),8)

                      cld_modified = max(cld,minValueCld)
                      cld_modified = min(cld_modified,maxValueCld)
                      cld_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) = real(cld_modified,4)

                  end do ! memberIndex
                end do ! stepIndex
                !$OMP END PARALLEL DO

            end do ! levIndex

          end do ! lonIndex
        end do ! latIndex
      end do ! varNameIndex

    end if

  end subroutine qlim_rttovLimit_ens

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

  !--------------------------------------------------------------------------
  ! qlim_setMin_ens
  !--------------------------------------------------------------------------
  subroutine qlim_setMin_ens(ensemble,huMinValue)
    !
    !:Purpose: To impose lower limit on humidity variable of an ensemble
    !
    implicit none

    ! Arguments:
    type(struct_ens), intent(inout) :: ensemble
    real(8),          intent(in)    :: huMinValue

    ! Locals:
    real(4), pointer :: hu_ptr_r4(:,:,:,:)
    real(4)          :: hu, hu_modified
    integer          :: lon1, lon2, lat1, lat2, numLev
    integer          :: lonIndex, latIndex, levIndex, stepIndex, memberIndex, varLevIndex

    if (mmpi_myid == 0) write(*,*) 'qlim_setMin_ens: STARTING'

    if (ens_getDataKind(ensemble) == 8) then
      call utl_abort('qlim_setMin_ens: Not compatible with dataKind = 8')
    end if

    if( .not. ens_varExist(ensemble,'HU') ) then
      if( mmpi_myid == 0 ) write(*,*) 'qlim_setMin_ens: ensemble does not ' // &
           'contain humidity ... doing nothing'
      return
    end if

    numLev = ens_getNumLev(ensemble,'TH')
    call ens_getLatLonBounds(ensemble, lon1, lon2, lat1, lat2)

    do latIndex = lat1, lat2
      do lonIndex = lon1, lon2
        do levIndex = 1, numLev
          varLevIndex = ens_getKFromLevVarName(ensemble, levIndex, 'HU')
          hu_ptr_r4 => ens_getOneLev_r4(ensemble,varLevIndex)

          !$OMP PARALLEL DO PRIVATE (stepIndex, memberIndex, hu, hu_modified)
          do stepIndex = 1, ens_getNumStep(ensemble)
            do memberIndex = 1, ens_getNumMembers(ensemble)
              hu = hu_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex)
              hu_modified = max(hu, real(huMinValue,4))
              hu_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) = hu_modified
            end do ! memberIndex
          end do ! stepIndex
          !$OMP END PARALLEL DO

        end do ! levIndex
      end do ! lonIndex
    end do ! latIndex

  end subroutine qlim_setMin_ens
  
  !-----------------------------------------------------------------------
  ! qlim_getMinValueCloud
  !----------------------------------------------------------------------
  function qlim_getMinValueCloud(varName) result(minValue)
    !
    ! :Purpose: Return the minValue for the hydrometeor.
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: varName
    real(8)                      :: minValue

    ! readNameList runs one time during program execution
    call readNameList

    select case (trim(varName))
    case ('LWCR')
      minValue = qlim_minValueLWCR
    case ('IWCR')
      minValue = qlim_minValueIWCR
    case ('RF')
      minValue = qlim_minValueRF
    case ('SF')
      minValue = qlim_minValueSF
    case ('CLDR')
      minValue = qlim_minValueCLDR      
    case default
      write(*,*)
      write(*,*) 'ERROR unknown varName: ', trim(varName)
      call utl_abort('qlim_getMinValueCloud')
   end select

  end function qlim_getMinValueCloud

  !-----------------------------------------------------------------------
  ! qlim_getMaxValueCloud
  !----------------------------------------------------------------------
  function qlim_getMaxValueCloud(varName) result(maxValue)
    !
    ! :Purpose: Return the maxValue for the hydrometeor.
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: varName
    real(8)                      :: maxValue

    ! readNameList runs one time during program execution
    call readNameList

    select case (trim(varName))
    case ('LWCR')
      maxValue = qlim_maxValueLWCR
    case ('IWCR')
      maxValue = qlim_maxValueIWCR
    case ('RF')
      maxValue = qlim_maxValueRF
    case ('SF')
      maxValue = qlim_maxValueSF
    case ('CLDR')
      maxValue = qlim_maxValueCLDR      
    case default
      write(*,*)
      write(*,*) 'ERROR unknown varName: ', trim(varName)
      call utl_abort('qlim_getMaxValueCloud')
   end select

  end function qlim_getMaxValueCloud

  !-----------------------------------------------------------------------
  ! cloudExistInEnsemble
  !----------------------------------------------------------------------
  function cloudExistInEnsemble(ensemble) result(cloudExist)
    !
    ! :Purpose: determine if any cloud variable exists in the ensemble.
    !
    implicit none

    ! Arguments:
    type(struct_ens), intent(in) :: ensemble
    logical                      :: cloudExist
    
    ! Locals:
    integer :: varNameIndex

    cloudExist = .false.
    do varNameIndex = 1, vnl_numvarmaxCloud
      if (ens_varExist(ensemble, vnl_varNameListCloud(varNameIndex))) then
        cloudExist = .true.
        return
      end if
    end do

  end function cloudExistInEnsemble

  !-----------------------------------------------------------------------
  ! cloudExistInStateVector
  !----------------------------------------------------------------------
  function cloudExistInStateVector(stateVector) result(cloudExist)
    !
    ! :Purpose: determine if any cloud variable exists in the stateVector.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in) :: stateVector
    logical                      :: cloudExist
    
    ! Locals:
    integer :: varNameIndex

    cloudExist = .false.
    do varNameIndex = 1, vnl_numvarmaxCloud
      if (gsv_varExist(stateVector, vnl_varNameListCloud(varNameIndex))) then
        cloudExist = .true.
        return
      end if
    end do

  end function cloudExistInStateVector

end module humidityLimits_mod
