
module humidityLimits_mod
  ! MODULE humidityLimits_mod (prefix='qlim' category='4. Data Object transformations')
  !
  !:Purpose: Various manipulations of humidity-related quantities.
  !
  use rmn_fnom
  use midasMpi_mod
  use utilities_mod
  use runtimeInfo_mod
  use mathPhysConstants_mod
  use varNameList_mod
  use physicsFunctions_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use gridStateVector_mod
  use ensembleStateVector_mod
  use calcHeightAndPressure_mod
  use columnData_mod
  use message_mod

  implicit none
  save
  private

  ! Public procedures
  public :: qlim_saturationLimit, qlim_rttovLimit, qlim_setMin, qlim_applyQcLimit
  public :: qlim_getMinValueCloud, qlim_getMaxValueCloud

  real(8), parameter :: mixratio_to_ppmv = 1.60771704d+6
  real(8) :: qlim_minValueLWCR, qlim_minValueIWCR, qlim_minValueRF, qlim_minValueSF, qlim_minValueCLDR
  real(8) :: qlim_maxValueLWCR, qlim_maxValueIWCR, qlim_maxValueRF, qlim_maxValueSF, qlim_maxValueCLDR
  real(8) :: qlim_qcThresh

  ! interface for qlim_saturationLimit
  interface qlim_saturationLimit
    module procedure qlim_saturationLimit_gsv
    module procedure qlim_saturationLimit_ens
  end interface qlim_saturationLimit

  ! interface for qlim_rttovLimit
  interface qlim_rttovLimit
    module procedure qlim_rttovLimit_gsv
    module procedure qlim_rttovLimit_ens
    module procedure qlim_rttovLimit_col
  end interface qlim_rttovLimit

  ! interface for qlim_applyQcLimit
  interface qlim_applyQcLimit
    module procedure qlim_applyQcLimit_gsv
    module procedure qlim_applyQcLimit_ens
  end interface qlim_applyQcLimit

  ! interface for qlim_setMin
  interface qlim_setMin
    module procedure qlim_setMin_ens
  end interface qlim_setMin

contains

  !--------------------------------------------------------------------------
  ! readNameList
  !--------------------------------------------------------------------------
  subroutine readNameList()
    !
    ! :Purpose: Reading NAMQLIM namelist by any subroutines in humidityLimits_mod module.
    !
    implicit none

    ! Locals:
    integer :: ierr
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
    real(8) :: qcThresh     ! threshold for QC


    NAMELIST /NAMQLIM/ minValueLWCR, maxValueLWCR, minValueIWCR, maxValueIWCR, &
                       minValueRF, maxValueRF, minValueSF, maxValueSF, minValueCLDR, maxValueCLDR, &
                       qcThresh

    if (nmlAlreadyRead) return

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

    qcThresh = 3.0d-8

    if (.not. utl_isNamelistPresent('NAMQLIM','./flnml')) then
      if (mmpi_myid == 0) then
        write(*,*) 'NAMQLIM is missing in the namelist. The default values will be taken.'
      end if

    else
      ! Reading the namelist
      call rti_tmg_start(181,'low-level--readNML')
      read(utl_flnml, nml=namqlim, iostat=ierr)
      if (ierr /= 0) call rti_abort('humidityLimits_mod: Error reading namelist')
      call rti_tmg_stop(181)

    end if
    if (mmpi_myid == 0) write(*,nml=namqlim)

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

    qlim_qcThresh = qcThresh

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
    type(struct_gsv), intent(inout) :: statevector ! state vector to modify

    ! Locals:
    type(struct_vco), pointer :: vco_ptr
    real(4), pointer :: hu_ptr_r4(:,:,:,:), tt_ptr_r4(:,:,:,:)
    real(8), pointer :: hu_ptr_r8(:,:,:,:), tt_ptr_r8(:,:,:,:)
    real(8), pointer :: pressure4D_T_r8(:,:,:,:), pressure4D_M_r8(:,:,:,:)
    real(4), pointer :: pressure4D_T_r4(:,:,:,:), pressure4D_M_r4(:,:,:,:)
    real(8), pointer :: height4D_T_r8(:,:,:,:), height4D_M_r8(:,:,:,:)
    real(4), pointer :: height4D_T_r4(:,:,:,:), height4D_M_r4(:,:,:,:)
    real(8)          :: hu, husat, hu_modified, tt
    integer          :: lonIndex, latIndex, levIndex, stepIndex

    if (mmpi_myid == 0) write(*,*) 'qlim_saturationLimit_gsv: STARTING'

    if (.not. gsv_varExist(statevector,'HU')) then
      if (mmpi_myid == 0) write(*,*) 'qlim_saturationLimit_gsv: statevector does not ' // &
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

    !
    !- Compute pressure (4D)
    !
    if (stateVector%dataKind == 8) then

      allocate(pressure4D_T_r8(statevector%myLonBeg:statevector%myLonEnd, &
                               statevector%myLatBeg:statevector%myLatEnd, &
                               gsv_getNumLev(statevector,'TH'), statevector%numStep))
      allocate(pressure4D_M_r8(statevector%myLonBeg:statevector%myLonEnd, &
                               statevector%myLatBeg:statevector%myLatEnd, &
                               gsv_getNumLev(statevector,'MM'), statevector%numStep))
      if (vco_ptr%vcode == 5002 .or. vco_ptr%vcode == 5005) then
        call czp_calcReturnPressure_gsv_nl(statevector,                  &
                                           PTout_r8_opt=pressure4D_T_r8, &
                                           PMout_r8_opt=pressure4D_M_r8)
      else if (vco_ptr%vcode == 21001) then
        allocate(height4D_T_r8(statevector%myLonBeg:statevector%myLonEnd, &
                               statevector%myLatBeg:statevector%myLatEnd, &
                               gsv_getNumLev(statevector,'TH'), statevector%numStep))
        allocate(height4D_M_r8(statevector%myLonBeg:statevector%myLonEnd, &
                               statevector%myLatBeg:statevector%myLatEnd, &
                               gsv_getNumLev(statevector,'MM'), statevector%numStep))
        call czp_calcReturnHeight_gsv_nl(statevector,                &
                                         ZTout_r8_opt=height4D_T_r8, &
                                         ZMout_r8_opt=height4D_M_r8)
        call czp_calcReturnPressure_gsv_nl(statevector,                  &
                                           ZTin_r8_opt=height4D_T_r8,    &
                                           ZMin_r8_opt=height4D_M_r8,    &
                                           PTout_r8_opt=pressure4D_T_r8, &
                                           PMout_r8_opt=pressure4D_M_r8)
        deallocate(height4D_M_r8)
        deallocate(height4D_T_r8)
      else
        call rti_abort('qlim_saturationLimit_gsv: Not compatible with this vCode: '//str(vco_ptr%vcode))
      end if
      deallocate(pressure4D_M_r8)

    else ! real 4

      allocate(pressure4D_T_r4(statevector%myLonBeg:statevector%myLonEnd, &
                               statevector%myLatBeg:statevector%myLatEnd, &
                               gsv_getNumLev(statevector,'TH'), statevector%numStep))
      allocate(pressure4D_M_r4(statevector%myLonBeg:statevector%myLonEnd, &
                               statevector%myLatBeg:statevector%myLatEnd, &
                               gsv_getNumLev(statevector,'MM'), statevector%numStep))
      allocate(pressure4D_T_r8(statevector%myLonBeg:statevector%myLonEnd, &
                               statevector%myLatBeg:statevector%myLatEnd, &
                               gsv_getNumLev(statevector,'TH'), statevector%numStep))
      if (vco_ptr%vcode == 5002 .or. vco_ptr%vcode == 5005) then
        call czp_calcReturnPressure_gsv_nl(statevector,                  &
                                           PTout_r4_opt=pressure4D_T_r4, &
                                           PMout_r4_opt=pressure4D_M_r4)
      else if (vco_ptr%vcode == 21001) then
        allocate(height4D_T_r4(statevector%myLonBeg:statevector%myLonEnd, &
                                 statevector%myLatBeg:statevector%myLatEnd, &
                                 gsv_getNumLev(statevector,'TH'), statevector%numStep))
        allocate(height4D_M_r4(statevector%myLonBeg:statevector%myLonEnd, &
                                 statevector%myLatBeg:statevector%myLatEnd, &
                                 gsv_getNumLev(statevector,'MM'), statevector%numStep))
        call czp_calcReturnHeight_gsv_nl(statevector,                &
                                         ZTout_r4_opt=height4D_T_r4, &
                                         ZMout_r4_opt=height4D_M_r4)
        call czp_calcReturnPressure_gsv_nl(statevector,                  &
                                           ZTin_r4_opt=height4D_T_r4,    &
                                           ZMin_r4_opt=height4D_M_r4,    &
                                           PTout_r4_opt=pressure4D_T_r4, &
                                           PMout_r4_opt=pressure4D_M_r4)
        deallocate(height4D_M_r4)
        deallocate(height4D_T_r4)
      else
        call rti_abort('qlim_saturationLimit_gsv: Not compatible with this vCode: '//str(vco_ptr%vcode))
      end if
      pressure4D_T_r8(:,:,:,:) = real(pressure4D_T_r4(:,:,:,:),8)
      deallocate(pressure4D_M_r4)
      deallocate(pressure4D_T_r4)

    end if ! stateVector%dataKind

    !
    !- Cap specific humidity (HU) one time step at a time
    !
    do stepIndex = 1, statevector%numStep

      if (stateVector%dataKind == 8) then
        !$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex, hu, tt, husat, hu_modified)
        do levIndex = 1, gsv_getNumLev(statevector,'TH')
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd
              hu = hu_ptr_r8(lonIndex,latIndex,levIndex,stepIndex)
              tt = tt_ptr_r8(lonIndex,latIndex,levIndex,stepIndex)
              ! get the saturated vapor pressure from HU
              husat = phf_foqst8(tt, pressure4D_T_r8(lonIndex,latIndex,levIndex,stepIndex))
              ! limit the humidity to the saturated humidity
              hu_modified = min(husat, hu)
              hu_ptr_r8(lonIndex,latIndex,levIndex,stepIndex) = hu_modified
            end do ! lonIndex
          end do ! latIndex
        end do ! levIndex
        !$OMP END PARALLEL DO
      else
        !$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex, hu, tt, husat, hu_modified)
        do levIndex = 1, gsv_getNumLev(statevector,'TH')
          do latIndex = statevector%myLatBeg, statevector%myLatEnd
            do lonIndex = statevector%myLonBeg, statevector%myLonEnd
              hu = hu_ptr_r4(lonIndex,latIndex,levIndex,stepIndex)
              tt = tt_ptr_r4(lonIndex,latIndex,levIndex,stepIndex)
              ! get the saturated vapor pressure from HU
              husat = phf_foqst8(tt, pressure4D_T_r8(lonIndex,latIndex,levIndex,stepIndex))
              ! limit the humidity to the saturated humidity
              hu_modified = min(husat, hu)
              hu_ptr_r4(lonIndex,latIndex,levIndex,stepIndex) =real( hu_modified,4)
            end do ! lonIndex
          end do ! latIndex
        end do ! levIndex
        !$OMP END PARALLEL DO
      end if

    end do ! stepIndex

    deallocate(pressure4D_T_r8)

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
    type(struct_ens), intent(inout) :: ensemble ! ensemble to modify

    ! Locals:
    type(struct_gsv)          :: stateVector
    type(struct_vco), pointer :: vco_ptr
    type(struct_hco), pointer :: hco_ptr
    real(4), pointer :: hu_ptr_r4(:,:,:,:), tt_ptr_r4(:,:,:,:), psfc_ptr_r4(:,:,:,:)
    real(4), pointer :: pressure_ptr_r4(:,:,:,:)
    real(8), pointer :: pressureEns(:,:,:)
    real(8)          :: hu, husat, hu_modified, tt
    integer          :: lon1, lon2, lat1, lat2, numLev, numStep, numMember
    real(8), allocatable :: psfc(:,:), pressure(:,:,:,:)
    integer          :: lonIndex, latIndex, levIndex, stepIndex, memberIndex, varLevIndex

    if (mmpi_myid == 0) write(*,*) 'qlim_saturationLimit_ens: STARTING'

    if (ens_getDataKind(ensemble) == 8) then
      call rti_abort('qlim_saturationLimit_ens: Not compatible with dataKind = 8')
    end if

    if (.not. ens_varExist(ensemble,'HU')) then
      if (mmpi_myid == 0) write(*,*) 'qlim_saturationLimit_ens: ensemble does not ' // &
           'contain humidity ... doing nothing'
      return
    end if

    vco_ptr => ens_getVco(ensemble)
    hco_ptr => ens_getHco(ensemble)

    numLev = ens_getNumLev(ensemble,'TH')
    numMember = ens_getNumMembers(ensemble)
    numStep = ens_getNumStep(ensemble)
    call ens_getLatLonBounds(ensemble, lon1, lon2, lat1, lat2)
    allocate(psfc(lon1:lon2,lat1:lat2))
    allocate(pressure(lon1:lon2,lat1:lat2,numLev,numStep))

    call gsv_allocate(stateVector, numStep,  &
                      hco_ptr, vco_ptr,  &
                      mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                      dataKind_opt=4, allocHeightSfc_opt=.true.,  &
                      varNames_opt=(/'P0 ','P_M','P_T','Z_M','Z_T','TT ','HU '/))
    call gsv_zero(stateVector)

    do memberIndex = 1, numMember

      if (vco_ptr%vcode == 21001) then

        call ens_copyMember(ensemble, stateVector, memberIndex)
        call czp_calcZandP_nl(stateVector)
        call gsv_getField(stateVector,pressure_ptr_r4,'P_T')
        pressure(:,:,:,:) = real(pressure_ptr_r4(:,:,:,:), 8)

      else

        varLevIndex = ens_getKFromLevVarName(ensemble, 1, 'P0')
        psfc_ptr_r4 => ens_getOneLev_r4(ensemble,varLevIndex)

        do stepIndex = 1, numStep
          psfc(:,:) = psfc_ptr_r4(memberIndex,stepIndex,:,:)
          nullify(pressureEns)
          if (vco_ptr%vcode == 5002 .or. vco_ptr%vcode == 5005) then
            call czp_fetch3DLevels(vco_ptr, psfc, fldT_opt=pressureEns)
          else
            write(*,*) 'vcode = ', vco_ptr%vcode
            call rti_abort('qlim_saturationLimit_ens: Unknown vcode value')
          end if
          pressure(:,:,:,stepIndex) = pressureEns(:,:,:)
          deallocate(pressureEns)
        end do ! stepIndex

      end if

      !$OMP PARALLEL DO PRIVATE (latIndex, lonIndex, levIndex, varLevIndex, hu_ptr_r4, tt_ptr_r4, stepIndex, hu, tt, husat, hu_modified)
      do latIndex = lat1, lat2
        do lonIndex = lon1, lon2

          do levIndex = 1, numLev
            varLevIndex = ens_getKFromLevVarName(ensemble, levIndex, 'HU')
            hu_ptr_r4 => ens_getOneLev_r4(ensemble,varLevIndex)
            varLevIndex = ens_getKFromLevVarName(ensemble, levIndex, 'TT')
            tt_ptr_r4 => ens_getOneLev_r4(ensemble,varLevIndex)

            do stepIndex = 1, numStep
              hu = hu_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex)
              tt = tt_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex)

              ! get the saturated vapor pressure from HU
              husat = phf_foqst8(tt, pressure(lonIndex,latIndex,levIndex,stepIndex))

              ! limit the humidity to the saturated humidity
              hu_modified = min(husat, hu)
              hu_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) = real(hu_modified,4)

            end do ! stepIndex

          end do ! levIndex

        end do ! lonIndex
      end do ! latIndex
      !$OMP END PARALLEL DO

    end do ! memberIndex

    deallocate(psfc)

  end subroutine qlim_saturationLimit_ens

  !--------------------------------------------------------------------------
  ! qlim_rttovLimit_gsv
  !--------------------------------------------------------------------------
  subroutine qlim_rttovLimit_gsv(statevector, applyLimitToHumidity_opt)
    !
    !:Purpose: To impose RTTOV limits on humidity/cloud
    !
    implicit none

    ! Arguments:
    type(struct_gsv),  intent(inout) :: statevector              ! state vector to be modified
    logical, optional, intent(in)    :: applyLimitToHumidity_opt ! apply limits to humidity variable

    ! Locals:
    type(struct_vco), pointer :: vco_ptr
    real(8), allocatable :: press_rttov(:), qmin_rttov(:), qmax_rttov(:)
    real(8), allocatable :: psfc(:,:)
    real(8), allocatable :: qmin3D_rttov(:,:,:), qmax3D_rttov(:,:,:)
    real(8), pointer     :: hu_ptr_r8(:,:,:,:), psfc_ptr_r8(:,:,:,:)
    real(8), pointer     :: cld_ptr_r8(:,:,:,:)
    real(4), pointer     :: hu_ptr_r4(:,:,:,:), psfc_ptr_r4(:,:,:,:)
    real(4), pointer     :: cld_ptr_r4(:,:,:,:)
    real(8), pointer     :: height4D_T_r8(:,:,:,:), height4D_M_r8(:,:,:,:)
    real(4), pointer     :: height4D_T_r4(:,:,:,:), height4D_M_r4(:,:,:,:)
    real(8), pointer     :: pressure3D(:,:,:)
    real(8), pointer     :: pressure4D_T_r8(:,:,:,:), pressure4D_M_r8(:,:,:,:)
    real(4), pointer     :: pressure4D_T_r4(:,:,:,:), pressure4D_M_r4(:,:,:,:)
    real(8)              :: hu, hu_modified, cld, cld_modified
    real(8)              :: minValueCld, maxValueCld
    integer              :: lon1, lon2, lat1, lat2, numStep, varNameIndex
    integer              :: lonIndex, latIndex, levIndex, stepIndex
    integer              :: ni, nj, numLev_T, numLev_M, numLev_rttov
    logical              :: applyLimitToHumidity, applyLimitToCloud

    if (mmpi_myid == 0) write(*,*) 'qlim_rttovLimit_gsv: STARTING'

    ! User can choose whether or not to apply limits to humidity variable
    if (present(applyLimitToHumidity_opt)) then
      applyLimitToHumidity = applyLimitToHumidity_opt
    else
      applyLimitToHumidity = .true.
    end if

    ! Always apply limits to cloud variables (if they are present)
    applyLimitToCloud = .true.

    ! Initialize some convenient variables
    vco_ptr => gsv_getVco(stateVector)
    numLev_T = gsv_getNumLev(stateVector,'TH')
    numLev_M = gsv_getNumLev(stateVector,'MM')
    numStep  = stateVector%numStep
    lon1 = statevector%myLonBeg
    lon2 = statevector%myLonEnd
    lat1 = statevector%myLatBeg
    lat2 = statevector%myLatEnd
    ni = lon2 - lon1 + 1
    nj = lat2 - lat1 + 1

    ! Apply limits to humidity
    if (applyLimitToHumidity .and. gsv_varExist(statevector,'HU')) then

      if (mmpi_myid == 0) write(*,*) 'qlim_rttovLimit_gsv: applying limits to HU.'

      ! Read in RTTOV humidity limits
      call readRttovLimitsFile(press_rttov, qmin_rttov, qmax_rttov, numLev_rttov)

      if (statevector%dataKind == 8) then
        call gsv_getField(statevector,hu_ptr_r8,'HU')
      else
        call gsv_getField(statevector,hu_ptr_r4,'HU')
      end if

      allocate(qmin3D_rttov(ni,nj,numLev_T))
      allocate(qmax3D_rttov(ni,nj,numLev_T))

      !
      !- Compute pressure (4D)
      !
      allocate(pressure4D_T_r8(lon1:lon2,lat1:lat2,numLev_T,numStep))

      if (vco_ptr%vcode == 5002 .or. vco_ptr%vcode == 5005) then
        !
        !- For GEM-P coordinate
        !
        allocate(psfc(ni,nj))

        do stepIndex = 1, numStep

          if (statevector%dataKind == 8) then
            call gsv_getField(statevector,psfc_ptr_r8,'P0')
            psfc(:,:) = psfc_ptr_r8(:,:,1,stepIndex)
          else
            call gsv_getField(statevector,psfc_ptr_r4,'P0')
            psfc(:,:) = real(psfc_ptr_r4(:,:,1,stepIndex),8)
          end if
          call czp_fetch3DLevels(vco_ptr, psfc, fldT_opt=pressure3d)

          pressure4D_T_r8(:,:,:,stepIndex) = pressure3d(:,:,:)
          deallocate(pressure3d)

        end do ! stepIndex

        deallocate(psfc)

      else if (vco_ptr%vcode == 21001) then
        !
        !- For GEM-H coordinate
        !
        if (statevector%dataKind == 8) then
          allocate(height4D_T_r8(lon1:lon2, lat1:lat2, numLev_T, numStep))
          allocate(height4D_M_r8(lon1:lon2, lat1:lat2, numLev_M, numStep))
          allocate(pressure4D_M_r8(lon1:lon2, lat1:lat2, numLev_M, numStep))
          call czp_calcReturnHeight_gsv_nl(statevector,                &
                                           ZTout_r8_opt=height4D_T_r8, &
                                           ZMout_r8_opt=height4D_M_r8)
          call czp_calcReturnPressure_gsv_nl(statevector,                  &
                                             ZTin_r8_opt=height4D_T_r8,    &
                                             ZMin_r8_opt=height4D_M_r8,    &
                                             PTout_r8_opt=pressure4D_T_r8, &
                                             PMout_r8_opt=pressure4D_M_r8)
          deallocate(height4D_M_r8)
          deallocate(height4D_T_r8)
          deallocate(pressure4D_M_r8)
        else
          allocate(height4D_T_r4(lon1:lon2, lat1:lat2, numLev_T, numStep))
          allocate(height4D_M_r4(lon1:lon2, lat1:lat2, numLev_M, numStep))
          allocate(pressure4D_T_r4(lon1:lon2, lat1:lat2, numLev_M, numStep))
          allocate(pressure4D_M_r4(lon1:lon2, lat1:lat2, numLev_M, numStep))
          call czp_calcReturnHeight_gsv_nl(statevector,                &
                                           ZTout_r4_opt=height4D_T_r4, &
                                           ZMout_r4_opt=height4D_M_r4)
          call czp_calcReturnPressure_gsv_nl(statevector,                  &
                                             ZTin_r4_opt=height4D_T_r4,    &
                                             ZMin_r4_opt=height4D_M_r4,    &
                                             PTout_r4_opt=pressure4D_T_r4, &
                                             PMout_r4_opt=pressure4D_M_r4)
          pressure4D_T_r8(:,:,:,:) = real(pressure4D_T_r4(:,:,:,:),8)
          deallocate(height4D_M_r4)
          deallocate(height4D_T_r4)
          deallocate(pressure4D_M_r4)
          deallocate(pressure4D_T_r4)
        end if

      else

        write(*,*) 'vcode = ', vco_ptr%vcode
        call rti_abort('qlim_rttovLimit_gsv: unknown vcode')

      end if

      ! Apply limits
      do stepIndex = 1, numStep

        ! Interpolate RTTOV limits onto model levels
        call qlim_lintv_minmax(press_rttov, qmin_rttov, qmax_rttov, numLev_rttov, &
                               ni, nj, numLev_T, pressure4D_T_r8(:,:,:,stepIndex), &
                               qmin3D_rttov, qmax3D_rttov)

        !$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex, hu, hu_modified)
        do levIndex = 1, numLev_T
          do latIndex = lat1, lat2
            do lonIndex = lon1, lon2
              if (statevector%dataKind == 8) then
                hu = hu_ptr_r8(lonIndex,latIndex,levIndex,stepIndex)
              else
                hu = real(hu_ptr_r4(lonIndex,latIndex,levIndex,stepIndex),8)
              end if

              ! limit the humidity according to the rttov limits
              hu_modified = max(hu, qmin3D_rttov(lonIndex - lon1 + 1, latIndex - lat1 + 1, levIndex))
              hu_modified = min(hu_modified, qmax3D_rttov(lonIndex - lon1 + 1, latIndex - lat1 + 1, levIndex))
              if (statevector%dataKind == 8) then
                hu_ptr_r8(lonIndex,latIndex,levIndex,stepIndex) = hu_modified
              else
                hu_ptr_r4(lonIndex,latIndex,levIndex,stepIndex) = real(hu_modified,4)
              end if

            end do ! lonIndex
          end do ! latIndex
        end do ! levIndex
        !$OMP END PARALLEL DO

      end do ! stepIndex

      deallocate(qmin3D_rttov)
      deallocate(qmax3D_rttov)

      deallocate(qmax_rttov)
      deallocate(qmin_rttov)
      deallocate(press_rttov)

      deallocate(pressure4D_T_r8)

    end if

    ! Apply limits to ALL available cloud variables
    if (applyLimitToCloud .and. cloudExistInStateVector(statevector)) then

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

        do stepIndex = 1, numStep

          !$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex, cld, cld_modified)
          do levIndex = 1, numLev_T
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

    if (mmpi_myid == 0) write(*,*) 'qlim_rttovLimit_gsv: FINISHED'

  end subroutine qlim_rttovLimit_gsv

  !--------------------------------------------------------------------------
  ! qlim_rttovLimit_ens
  !--------------------------------------------------------------------------
  subroutine qlim_rttovLimit_ens(ensemble, applyLimitToHumidity_opt)
    !
    !:Purpose: To impose RTTOV limits on humidity/cloud
    !
    implicit none

    ! Arguments:
    type(struct_ens),  intent(inout) :: ensemble                 ! ensemble that will be modified
    logical, optional, intent(in)    :: applyLimitToHumidity_opt ! apply limits to humidity variable

    ! Locals:
    type(struct_gsv)          :: stateVector
    type(struct_vco), pointer :: vco_ptr
    type(struct_hco), pointer :: hco_ptr
    real(8), allocatable :: press_rttov(:), qmin_rttov(:), qmax_rttov(:)
    real(8), allocatable :: psfc(:,:)
    real(8), allocatable :: qmin3D_rttov(:,:,:), qmax3D_rttov(:,:,:)
    real(4), pointer     :: pressure4D_ptr_r4(:,:,:,:)
    real(8), pointer     :: pressure3D_ptr_r8(:,:,:)
    real(4), pointer     :: hu_ptr_r4(:,:,:,:), psfc_ptr_r4(:,:,:,:)
    real(4), pointer     :: cld_ptr_r4(:,:,:,:)
    real(8), allocatable :: pressure4D_r8(:,:,:,:)
    real(8)              :: hu, hu_modified, cld, cld_modified
    real(8)              :: minValueCld, maxValueCld
    integer              :: lon1, lon2, lat1, lat2, ni, nj, varNameIndex
    integer              :: lonIndex, latIndex, levIndex, stepIndex, varLevIndex, memberIndex
    integer              :: numMember, numStep, numLev_M, numLev_T, numLev_rttov
    logical              :: applyLimitToHumidity, applyLimitToCloud

    if (mmpi_myid == 0) write(*,*) 'qlim_rttovLimit_ens: STARTING'

    if (ens_getDataKind(ensemble) == 8) then
      call rti_abort('qlim_rttovLimit_ens: Not compatible with dataKind = 8')
    end if

    ! User can choose whether or not to apply limits to humidity variable
    if (present(applyLimitToHumidity_opt)) then
      applyLimitToHumidity = applyLimitToHumidity_opt
    else
      applyLimitToHumidity = .true.
    end if

    ! Always apply limits to cloud variables (if they are present)
    applyLimitToCloud = .true.

    ! Initialize some convenient variables
    vco_ptr => ens_getVco(ensemble)
    hco_ptr => ens_getHco(ensemble)

    numLev_T = ens_getNumLev(ensemble,'TH')
    numLev_M = ens_getNumLev(ensemble,'MM')
    numMember = ens_getNumMembers(ensemble)
    numStep = ens_getNumStep(ensemble)
    call ens_getLatLonBounds(ensemble, lon1, lon2, lat1, lat2)
    ni = lon2 - lon1 + 1
    nj = lat2 - lat1 + 1

    ! Apply limits to humidity
    if (applyLimitToHumidity .and. ens_varExist(ensemble,'HU')) then

      if (mmpi_myid == 0) write(*,*) 'qlim_rttovLimit_ens:  applying limits to HU.'

      ! Read in RTTOV humidity limits
      call readRttovLimitsFile(press_rttov, qmin_rttov, qmax_rttov, numLev_rttov)

      call gsv_allocate(stateVector, numStep,  &
                        hco_ptr, vco_ptr,  &
                        mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                        dataKind_opt=4, allocHeightSfc_opt=.true.,  &
                        varNames_opt=(/'P0 ','P_M','P_T','Z_M','Z_T','TT ','HU '/))
      call gsv_zero(stateVector)

      allocate(qmin3D_rttov(ni,nj,numLev_T))
      allocate(qmax3D_rttov(ni,nj,numLev_T))
      allocate(pressure4D_r8(lon1:lon2,lat1:lat2,numLev_T,numStep))

      do memberIndex = 1, numMember

        if (vco_ptr%vcode == 21001) then
          !
          !- GEM-H
          !
          call ens_copyMember(ensemble, stateVector, memberIndex)
          call czp_calcZandP_nl(stateVector)
          call gsv_getField(stateVector,pressure4D_ptr_r4,'P_T')
          pressure4D_r8(:,:,:,:) = real(pressure4D_ptr_r4(:,:,:,:), 8)

        else
          !
          !- GEM-P
          !
          varLevIndex = ens_getKFromLevVarName(ensemble, 1, 'P0')
          psfc_ptr_r4 => ens_getOneLev_r4(ensemble,varLevIndex)

          do stepIndex = 1, numStep
            nullify(pressure3D_ptr_r8)
            if (.not. allocated(psfc)) allocate(psfc(ni,nj))
            psfc(:,:) = psfc_ptr_r4(memberIndex,stepIndex,:,:)
            write(*,*) 'psfc min/max = ', minval(psfc), maxval(psfc)
            if (vco_ptr%vcode == 5002 .or. vco_ptr%vcode == 5005) then
              call czp_fetch3DLevels(vco_ptr, psfc, fldT_opt=pressure3D_ptr_r8)
            else
              write(*,*) 'vcode = ', vco_ptr%vcode
              call rti_abort('qlim_rttovLimit_ens: Unknown vcode value')
            end if
            pressure4D_r8(:,:,:,stepIndex) = pressure3D_ptr_r8(:,:,:)
            deallocate(pressure3D_ptr_r8)
          end do ! stepIndex

        end if

        do stepIndex = 1, numStep

          ! Interpolate RTTOV limits onto model levels
          call qlim_lintv_minmax(press_rttov, qmin_rttov, qmax_rttov, numLev_rttov, &
                                 ni, nj, numLev_T, pressure4D_r8(:,:,:,stepIndex),  &
                                 qmin3D_rttov, qmax3D_rttov)

          do levIndex = 1, numLev_T

            varLevIndex = ens_getKFromLevVarName(ensemble, levIndex, 'HU')
            hu_ptr_r4 => ens_getOneLev_r4(ensemble,varLevIndex)

            !$OMP PARALLEL DO PRIVATE (latIndex, lonIndex, hu, hu_modified)
            do latIndex = lat1, lat2
              do lonIndex = lon1, lon2

                hu = real(hu_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex),8)

                ! limit the humidity according to the rttov limits
                hu_modified = max(hu, qmin3D_rttov(lonIndex - lon1 + 1, latIndex - lat1 + 1, levIndex))
                hu_modified = min(hu_modified, qmax3D_rttov(lonIndex - lon1 + 1, latIndex - lat1 + 1, levIndex))
                hu_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) = real(hu_modified,4)

              end do ! lonIndex
            end do ! latIndex
            !$OMP END PARALLEL DO

          end do ! levIndex

        end do ! stepIndex

      end do ! memberIndex

      if (allocated(psfc)) deallocate(psfc)
      deallocate(qmin3D_rttov)
      deallocate(qmax3D_rttov)
      deallocate(pressure4D_r8)

      deallocate(qmax_rttov)
      deallocate(qmin_rttov)
      deallocate(press_rttov)

    end if

    ! Apply limits to ALL available cloud variables
    if (applyLimitToCloud .and. cloudExistInEnsemble(ensemble)) then

      do varNameIndex = 1, vnl_numvarmaxCloud
        if (.not. ens_varExist(ensemble, vnl_varNameListCloud(varNameIndex))) cycle

        if (mmpi_myid == 0) write(*,*) 'qlim_rttovLimit_ens:  applying limits to ', &
                                        vnl_varNameListCloud(varNameIndex)

        minValueCld = qlim_getMinValueCloud(vnl_varNameListCloud(varNameIndex))
        maxValueCld = qlim_getMaxValueCloud(vnl_varNameListCloud(varNameIndex))

        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2

            do levIndex = 1, numLev_T
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
  ! qlim_applyQcLimit_gsv
  !--------------------------------------------------------------------------
  subroutine qlim_applyQcLimit_gsv(statevector)
    !
    !:Purpose: To impose limit on QC
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: statevector ! statevector object to be modified

    ! Locals:
    real(8), pointer :: qc_ptr_r8(:,:,:,:)
    real(4), pointer :: qc_ptr_r4(:,:,:,:)
    real(8)          :: qc, qc_modified
    integer          :: lon1, lon2, lat1, lat2, numLev_T
    integer          :: lonIndex, latIndex, levIndex, stepIndex
    character(len=4) :: varName

    varName = 'QC'
    if (.not. gsv_varExist(stateVector,varName)) return

    if (mmpi_myid == 0) write(*,*) 'qlim_applyQcLimit_gsv: applying limits to QC'

    if (statevector%dataKind == 8) then
      call gsv_getField(statevector, qc_ptr_r8, varName)
    else
      call gsv_getField(statevector, qc_ptr_r4, varName)
    end if

    lon1 = statevector%myLonBeg
    lon2 = statevector%myLonEnd
    lat1 = statevector%myLatBeg
    lat2 = statevector%myLatEnd
    numLev_T = gsv_getNumLev(statevector,'TH')

    do stepIndex = 1, statevector%numStep

      !$OMP PARALLEL DO PRIVATE (levIndex, latIndex, lonIndex, qc, qc_modified)
      do levIndex = 1, numLev_T
        do latIndex = lat1, lat2
          do lonIndex = lon1, lon2
            if (statevector%dataKind == 8) then
              qc = qc_ptr_r8(lonIndex,latIndex,levIndex,stepIndex)
            else
              qc = real(qc_ptr_r4(lonIndex,latIndex,levIndex,stepIndex),8)
            end if

            qc_modified = qc
            if (qc <  qlim_qcThresh) qc_modified = 0.0d0

            if (statevector%dataKind == 8) then
              qc_ptr_r8(lonIndex,latIndex,levIndex,stepIndex) = qc_modified
            else
              qc_ptr_r4(lonIndex,latIndex,levIndex,stepIndex) = real(qc_modified,4)
            end if

          end do ! lonIndex
        end do ! latIndex
      end do ! levIndex
      !$OMP END PARALLEL DO

    end do ! stepIndex

  end subroutine qlim_applyQcLimit_gsv

  !--------------------------------------------------------------------------
  ! qlim_applyQcLimit_ens
  !--------------------------------------------------------------------------
  subroutine qlim_applyQcLimit_ens(ensemble)
    !
    !:Purpose: To impose limit on QC
    !
    implicit none

    ! Arguments:
    type(struct_ens), intent(inout) :: ensemble ! ensemble object to modify

    ! Locals:
    real(4), pointer :: qc_ptr_r4(:,:,:,:)
    real(8)          :: qc, qc_modified
    integer          :: lon1, lon2, lat1, lat2
    integer          :: numLev, numMember, numStep, varLevIndex
    integer          :: lonIndex, latIndex, levIndex, stepIndex, memberIndex
    character(len=4) :: varName

    varName = 'QC'
    if (.not. ens_varExist(ensemble,varName)) return

    if (mmpi_myid == 0) write(*,*) 'qlim_applyQcLimit_ens: applying limits to QC'

    numLev = ens_getNumLev(ensemble,'TH')
    numMember = ens_getNumMembers(ensemble)
    numStep = ens_getNumStep(ensemble)
    call ens_getLatLonBounds(ensemble, lon1, lon2, lat1, lat2)

    do latIndex = lat1, lat2
      do lonIndex = lon1, lon2

        do levIndex = 1, numLev
            varLevIndex = ens_getKFromLevVarName(ensemble, levIndex, varName)
            qc_ptr_r4 => ens_getOneLev_r4(ensemble,varLevIndex)

            !$OMP PARALLEL DO PRIVATE (stepIndex, memberIndex, qc, qc_modified)
            do stepIndex = 1, numStep
              do memberIndex = 1, numMember

                  qc = real(qc_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex),8)

                  qc_modified = qc
                  if (qc <  qlim_qcThresh) qc_modified = 0.0d0

                  qc_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) = real(qc_modified,4)

              end do ! memberIndex
            end do ! stepIndex
            !$OMP END PARALLEL DO

        end do ! levIndex

      end do ! lonIndex
    end do ! latIndex

  end subroutine qlim_applyQcLimit_ens

  !--------------------------------------------------------------------------
  ! qlim_lintv_minmax
  !--------------------------------------------------------------------------
  subroutine qlim_lintv_minmax(press_src, qmin_src, qmax_src, numLev_src, &
                               ni_dest, nj_dest, numLev_dest, press_dest, &
                               qmin_dest, qmax_dest)
    !
    !:Purpose: To perform the vertical interpolation in log of pressure and
    !          and constant value extrapolation of one-dimensional vectors.
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: numLev_src            ! Number of input levels (source)
    real(8), intent(in)  :: press_src(numLev_src) ! Vertical levels, pressure (source)
    real(8), intent(in)  :: qmin_src(numLev_src)  ! Vectors to be interpolated (source)
    real(8), intent(in)  :: qmax_src(numLev_src)  ! Vectors to be interpolated (source)
    integer, intent(in)  :: ni_dest               ! Number of profiles in i dimension
    integer, intent(in)  :: nj_dest               ! Number of profiles in j dimension
    integer, intent(in)  :: numLev_dest           ! Number of output levels (destination)
    real(8), intent(in)  :: press_dest(:,:,:)     ! Vertical levels, pressure (destination)
    real(8), intent(out) :: qmin_dest(:,:,:)      ! Interpolated profiles (destination)
    real(8), intent(out) :: qmax_dest(:,:,:)      ! Interpolated profiles (destination)

    ! Locals:
    integer :: ji, jk, jo, ii, jj, ik, iorder, xi
    real(8) :: zpo(numLev_dest)
    integer :: il(numLev_dest)
    real(8) :: zpi(0:numLev_src+1)
    real(8) :: zqmin_src(0:numLev_src+1), zqmax_src(0:numLev_src+1)
    real(8) :: zw1, zw2, zp, zrt, zp1, zp2

    zpi(0)=200000.d0
    zpi(numLev_src+1)=200000.d0

    ! Determine if input pressure levels are in ascending or
    ! descending order.
    if (press_src(1) < press_src(numLev_src)) then
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
            xi = int(sign(1.0d0,iorder*(zrt-zp)))
            il(jo) = il(jo) + max(0, xi)
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
    type(struct_ens), intent(inout) :: ensemble   ! ensemble to be modified
    real(8),          intent(in)    :: huMinValue ! lower limit on HU value to impose

    ! Locals:
    real(4), pointer :: hu_ptr_r4(:,:,:,:)
    real(4)          :: hu, hu_modified
    integer          :: lon1, lon2, lat1, lat2, numLev
    integer          :: lonIndex, latIndex, levIndex, stepIndex, memberIndex, varLevIndex

    if (mmpi_myid == 0) write(*,*) 'qlim_setMin_ens: STARTING'

    if (ens_getDataKind(ensemble) == 8) then
      call rti_abort('qlim_setMin_ens: Not compatible with dataKind = 8')
    end if

    if (.not. ens_varExist(ensemble,'HU')) then
      if (mmpi_myid == 0) write(*,*) 'qlim_setMin_ens: ensemble does not ' // &
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
    character(len=*), intent(in) :: varName  ! variable name
    ! Result:
    real(8)                      :: minValue ! minimum value for this variable

    ! readNameList runs one time during program execution
    call readNameList()

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
      call rti_abort('qlim_getMinValueCloud')
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
    character(len=*), intent(in) :: varName  ! variable name
    ! Result:
    real(8)                      :: maxValue ! maximum value for this variable

    ! readNameList runs one time during program execution
    call readNameList()

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
      call rti_abort('qlim_getMaxValueCloud')
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
    type(struct_ens), intent(in) :: ensemble   ! ensemble to be examined
    ! Result:
    logical                      :: cloudExist ! indicate if any cloud variable exists in ensemble

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
    type(struct_gsv), intent(in) :: stateVector ! state vector to be examined
    ! Result:
    logical                      :: cloudExist  ! indicate if any clound variable exists in state

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

  !--------------------------------------------------------------------------
  ! qlim_rttovLimit_col
  !--------------------------------------------------------------------------
  subroutine qlim_rttovLimit_col(column)
    !
    !:Purpose: To impose RTTOV limits on humidity/LWCR for column data objects
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(inout) :: column  ! column object that will be modified

    ! Locals:
    integer              :: numLev_rttov, numCol, numLev_T
    real(8), allocatable :: press_rttov(:), qmin_rttov(:), qmax_rttov(:)
    real(8), allocatable :: qmin2D_rttov(:,:), qmax2D_rttov(:,:)
    integer              :: levIndex, columnIndex
    real(8), pointer     :: huColPtr(:,:), pressTColPtr(:,:)
    real(8)              :: hu, hu_modified

    if (mmpi_myid == 0) write(*,*) 'qlim_rttovLimit_col: STARTING'

    if (.not. col_varExist(column,'HU')) then
      if (mmpi_myid == 0) write(*,*) 'qlim_rttovLimit_col: column does not ' // &
           'contain humidity ... doing nothing'
      return
    end if

    ! Read in RTTOV humidity limits
    call readRttovLimitsFile(press_rttov, qmin_rttov, qmax_rttov, numLev_rttov)

    huColPtr => col_getAllColumns(column,'HU')
    pressTColPtr => col_getAllColumns(column,'P_T')

    numCol = col_getNumCol(column)
    numLev_T = col_getNumLev(column,'TH')

    write(*,*) 'numLev', numLev_T
    allocate(qmin2D_rttov(numCol,numLev_T))
    allocate(qmax2D_rttov(numCol,numLev_T))

    call qlim_lintv_minmax_col(press_rttov, qmin_rttov, qmax_rttov, numLev_rttov, &
                               numCol, numLev_T, pressTColPtr, qmin2D_rttov, qmax2D_rttov)

    !$OMP PARALLEL DO PRIVATE (levIndex, columnIndex, hu, hu_modified)
    do columnIndex = 1, numCol
      do levIndex = 1, numLev_T
        hu = huColPtr(levIndex, columnIndex)

        ! limit the humidity according to the rttov limits
        hu_modified = max(hu, qmin2D_rttov(columnIndex, levIndex))
        hu_modified = min(hu_modified, qmax2D_rttov(columnIndex, levIndex))

        huColPtr(levIndex, columnIndex) = hu_modified
      end do ! levIndex
    end do ! columnIndex
    !$OMP END PARALLEL DO

    deallocate(qmin2D_rttov)
    deallocate(qmax2D_rttov)
    deallocate(qmax_rttov)
    deallocate(qmin_rttov)
    deallocate(press_rttov)

  end subroutine qlim_rttovLimit_col

  !--------------------------------------------------------------------------
  ! readRttovLimitsFile
  !--------------------------------------------------------------------------
  subroutine readRttovLimitsFile(press_rttov, qmin_rttov, qmax_rttov, numLev_rttov)
    !
    !:Purpose: Read the contents of the file containing the RTTOV
    !          limits on pressure levels
    !
    implicit none

    ! Arguments:
    real(8), allocatable, intent(out) :: press_rttov(:) ! pressure levels for rttov limits
    real(8), allocatable, intent(out) :: qmin_rttov(:)  ! minimum values for each level
    real(8), allocatable, intent(out) :: qmax_rttov(:)  ! maximum values for each level
    integer,              intent(out) :: numLev_rttov   ! number of pressure levels in the file

    ! Locals:
    character(len=*), parameter :: fileName = 'rttov_h2o_limits.dat'
    integer :: nulfile, ierr, levIndex
    logical, save :: firstTime = .true.

    ! Open the file
    nulfile = 0
    ierr = fnom(nulfile, fileName, "FMT+OLD+R/O", 0)
    if (ierr /= 0) then
      if (mmpi_myid == 0) write(*,*) 'fileName = ', fileName
      call rti_abort('qlim_rttovLimit_col: error opening the humidity limits file')
    end if

    ! Read the contents
    read(nulfile,*) numLev_rttov
    if (mmpi_myid == 0 .and. firstTime) write(*,*) 'qlim_rttovLimit_col: rttov number of levels = ', numLev_rttov
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

    ! Print the file contents to the listing
    if (firstTime) then
      write(*,*) ' '
      do levIndex = 1, numLev_rttov
        if (mmpi_myid == 0) write(*,fmt='(" qlim_rttovLimit_col:   LEVEL = ",I4,", PRES = ",F9.0,", HUMIN = ",E10.2,", HUMAX = ",E10.2)') &
             levIndex, press_rttov(levIndex), qmin_rttov(levIndex), qmax_rttov(levIndex)
      end do
      firstTime = .false.
    end if

  end subroutine readRttovLimitsFile

  !--------------------------------------------------------------------------
  ! qlim_lintv_minmax_col
  !--------------------------------------------------------------------------
  subroutine qlim_lintv_minmax_col(press_src, qmin_src, qmax_src, numLev_src, &
                                   numColumn_dest, numLev_dest, press_dest, &
                                   qmin_dest, qmax_dest)

    !
    !:Purpose: To perform the vertical interpolation in log of pressure and
    !          and constant value extrapolation of one-dimensional column.
    implicit none

    ! Arguments:
    integer, intent(in)  :: numLev_src            ! Number of input levels (source)
    real(8), intent(in)  :: press_src(numLev_src) ! Vertical levels, pressure (source)
    real(8), intent(in)  :: qmin_src(numLev_src)  ! Vectors to be interpolated (source)
    real(8), intent(in)  :: qmax_src(numLev_src)  ! Vectors to be interpolated (source)
    integer, intent(in)  :: numColumn_dest        ! Number of profiles
    integer, intent(in)  :: numLev_dest           ! Number of output levels (destination)
    real(8), intent(in)  :: press_dest(:,:)       ! Vertical levels, pressure (destination)
    real(8), intent(out) :: qmin_dest(:,:)        ! Interpolated profiles (destination)
    real(8), intent(out) :: qmax_dest(:,:)        ! Interpolated profiles (destination)

    ! Locals:
    integer :: ji, jk, jo, ii, ik, iorder, xi
    real(8) :: zpo(numLev_dest)
    integer :: il(numLev_dest)
    real(8) :: zpi(0:numLev_src+1)
    real(8) :: zqmin_src(0:numLev_src+1), zqmax_src(0:numLev_src+1)
    real(8) :: zw1, zw2, zp, zrt, zp1, zp2

    zpi(0)=200000.d0
    zpi(numLev_src+1)=200000.d0

    ! Determine if input pressure levels are in ascending or
    ! descending order.
    if (press_src(1) < press_src(numLev_src)) then
      iorder = 1
    else
      iorder = -1
    end if

    ! Source levels
    do jk = 1, numLev_src
      zpi(jk) = press_src(jk)
      zqmin_src(jk) = qmin_src(jk)
      zqmax_src(jk) = qmax_src(jk)
    end do

    zqmin_src(0) = qmin_src(1)
    zqmin_src(numLev_src+1) = qmin_src(numLev_src)
    zqmax_src(0) = qmax_src(1)
    zqmax_src(numLev_src+1) = qmax_src(numLev_src)

    do ii = 1, numColumn_dest
      zpo(:) = press_dest(:,ii)

      ! Find the adjacent level below
      il(:) = 0
      do ji = 1, numLev_src
        do jo = 1, numLev_dest
          zrt = zpo(jo)
          zp = zpi(ji)
          xi = int(sign(1.0d0,iorder*(zrt-zp)))
          il(jo) = il(jo) + max(0, xi)
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
        qmin_dest(ii,jo) = zw1*zqmin_src(ik) +  zw2*zqmin_src(ik+1)
        qmax_dest(ii,jo) = zw1*zqmax_src(ik) +  zw2*zqmax_src(ik+1)
      end do

    end do
  end subroutine qlim_lintv_minmax_col

end module humidityLimits_mod
