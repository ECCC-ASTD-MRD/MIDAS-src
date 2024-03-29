
module costFunction_mod
  ! MODULE costFunction_mod, (prefix='cfn' category='5. Observation operators')
  !
  !:Purpose: To compute Jo term.
  !
  use midasMpi_mod
  use obsSpaceData_mod
  use rttov_const, only : inst_name, platform_name
  use tovsNL_mod
  use utilities_mod
  use obserrors_mod
  use codtyp_mod

  implicit none

  save
  private

  ! public procedures

  public :: cfn_calcJo, cfn_sumJo

  integer,           allocatable :: channelNumberList(:,:)
  character(len=15), allocatable :: sensorNameList(:)

contains

  !--------------------------------------------------------------------------
  ! cfn_calcJo
  !--------------------------------------------------------------------------
  subroutine cfn_calcJo(lobsSpaceData)
    !
    !:Purpose: To compute JO contribution of each assimilated and diagnosed
    !          datum, and to store the result in OBS_JOBS
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: lobsSpaceData

    ! Locals:
    integer :: bodyIndex

    !$OMP PARALLEL DO PRIVATE(bodyIndex)
    do bodyIndex=1,obs_numbody(lobsSpaceData)

      if ( obs_bodyElem_i( lobsSpaceData, OBS_ASS, bodyIndex) == obs_assimilated) then
        call obs_bodySet_r( lobsSpaceData, OBS_JOBS, bodyIndex, &
             ( obs_bodyElem_r( lobsSpaceData, OBS_WORK, bodyIndex ) &
             * obs_bodyElem_r( lobsSpaceData, OBS_WORK, bodyIndex ) &
             ) / 2.d0 )
      else
        call obs_bodySet_r(lobsSpaceData,OBS_JOBS,bodyIndex, 0.0d0)
      end if

    end do
    !$OMP END PARALLEL DO

  end subroutine cfn_calcJo

  !--------------------------------------------------------------------------
  ! cfn_sumJo
  !--------------------------------------------------------------------------
  subroutine cfn_sumJo( lobsSpaceData, pjo, beSilent_opt )
    !
    !:Purpose: To compute the sum of Jo contributions saved in OBS_JOBS. Also,
    !          to compute contribution of each family of observation (for
    !          diagnostic purposes)
    implicit none

    ! Arguments:
    type(struct_obs),    intent(in)  :: lobsSpaceData
    real(8),             intent(out) :: pjo ! Total observation cost function
    logical, optional,   intent(in) :: beSilent_opt

    ! Locals:
    integer :: bodyIndex, tovsIndex, sensorIndex, headerIndex, bodyIndexBeg, bodyIndexEnd
    integer :: channelNumber, channelNumberIndexInListFound, channelIndex
    integer :: sensorIndexInList, sensorIndexInListFound
    logical :: beSilent

    real(8) :: dljoraob, dljoairep, dljosatwind, dljoscat, dljosurfc, dljotov, dljosst, dljoice
    real(8) :: dljoprof, dljogpsro, dljogpsztd, dljochm, pjo_1, dljoaladin, dljohydro, dljoradar
    real(8) :: dljotov_sensors( tvs_nsensors )
    real(8) :: joTovsPerChannelSensor(tvs_maxNumberOfChannels,tvs_nsensors)
    character(len=15) :: lowerCaseName

    logical :: printJoTovsPerChannelSensor
  
    real(8), allocatable :: joSSTInstrument(:)
    integer, allocatable :: nobsInstrument(:), nobsInstrumentGlob(:)
    integer :: SSTdatasetIndex, codeType, ierr

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .false.
    end if

    if (.not. allocated(channelNumberList)) then
      allocate(channelNumberList(tvs_maxNumberOfChannels,tvs_nsensors))
    end if
    if (.not. allocated(sensorNameList)) then
      allocate(sensorNameList(tvs_nsensors))
    end if
    
    if(oer_getSSTdataParam_int('numberSSTDatasets') > 0) then
      allocate(joSSTInstrument(oer_getSSTdataParam_int('numberSSTDatasets')))
      allocate(nobsInstrument(oer_getSSTdataParam_int('numberSSTDatasets')))
      allocate(nobsInstrumentGlob(oer_getSSTdataParam_int('numberSSTDatasets')))
    end if

    call readNameList
    printJoTovsPerChannelSensor = .false.
    if (any(sensorNameList(:) /= '') .and. any(channelNumberList(:,:) /= 0)) then
      printJoTovsPerChannelSensor = .true.
    end if

    dljogpsztd = 0.d0
    dljoraob = 0.d0
    dljoairep = 0.d0
    dljosatwind = 0.d0
    dljosurfc = 0.d0
    dljoscat = 0.d0
    dljotov = 0.d0
    dljogpsro = 0.d0
    dljoprof = 0.d0
    dljochm = 0.d0
    dljosst = 0.0d0
    dljoaladin = 0.d0
    dljoice = 0.0d0
    dljotov_sensors(:) = 0.d0
    joTovsPerChannelSensor(:,:) = 0.0d0
    dljohydro = 0.0d0
    dljoradar = 0.0d0
    joSSTInstrument(:) = 0.0d0
    nobsInstrumentGlob(:) = 0
    nobsInstrument(:) = 0

    pjo = 0.0d0

    do bodyIndex = 1, obs_numbody(lobsSpaceData)

      pjo_1 = obs_bodyElem_r(lobsSpaceData, OBS_JOBS, bodyIndex)

      ! total observation cost function
      pjo   = pjo + pjo_1

      ! subcomponents of observation cost function (diagnostic only)
      select case(obs_getFamily(lobsSpaceData, bodyIndex_opt = bodyIndex))
      case('UA')
        dljoraob    = dljoraob    + pjo_1
      case('AI')
        dljoairep   = dljoairep   + pjo_1
      case('SW')
        dljosatwind = dljosatwind + pjo_1
      case('SF')
        dljosurfc   = dljosurfc   + pjo_1
      case('SC')
        dljoscat    = dljoscat    + pjo_1
      case('TO')
        dljotov     = dljotov     + pjo_1
      case('RO')
        dljogpsro   = dljogpsro   + pjo_1
      case('PR')
        dljoprof    = dljoprof    + pjo_1
      case('GP')
        dljogpsztd  = dljogpsztd  + pjo_1
      case('CH')
        dljochm     = dljochm     + pjo_1
      case('TM')
        dljosst     = dljosst     + pjo_1
      case('AL')
        dljoaladin  = dljoaladin  + pjo_1
      case('GL')
        dljoice     = dljoice     + pjo_1
      case('HY')
        dljohydro   = dljohydro   + pjo_1
      case('RA')
        dljoradar   = dljoradar   + pjo_1
      end select
    enddo

    do tovsIndex = 1, tvs_nobtov
      headerIndex = tvs_headerIndex( tovsIndex )
      if ( headerIndex > 0 ) then
        bodyIndexBeg = obs_headElem_i(lobsSpaceData, OBS_RLN, headerIndex)
        bodyIndexEnd = obs_headElem_i(lobsSpaceData, OBS_NLV, headerIndex) + bodyIndexBeg - 1
        sensorIndex = tvs_lsensor (tovsIndex)

        sensorIndexInListFound = 0
        if ( printJoTovsPerChannelSensor ) then
          loopSensor1: do sensorIndexInList = 1, tvs_nsensors
            call up2low(sensorNameList(sensorIndexInList),lowerCaseName)

            if ( trim(lowerCaseName) == trim(inst_name(tvs_instruments(sensorIndex))) ) then
              sensorIndexInListFound = sensorIndexInList
              exit loopSensor1
            end if

          end do loopSensor1
        end if

        do bodyIndex = bodyIndexBeg, bodyIndexEnd
          pjo_1 = obs_bodyElem_r(lobsSpaceData, OBS_JOBS, bodyIndex)
          dljotov_sensors(sensorIndex) =  dljotov_sensors(sensorIndex) + pjo_1

          if ( printJoTovsPerChannelSensor .and. &
               sensorIndexInListFound > 0 ) then
            call tvs_getChannelNumIndexFromPPP( lobsSpaceData, headerIndex, bodyIndex, &
                                                channelNumber, channelIndex )
            channelNumberIndexInListFound = utl_findloc(channelNumberList(:,sensorIndexInListFound), &
                                                        channelNumber)

            if ( channelNumberIndexInListFound > 0 ) then
              joTovsPerChannelSensor(channelNumberIndexInListFound,sensorIndexInListFound) = &
                        joTovsPerChannelSensor(channelNumberIndexInListFound,sensorIndexInListFound) + &
                        pjo_1
            end if

          end if

        end do
      end if
    end do

    if(oer_getSSTdataParam_int('numberSSTDatasets') > 0) then
      do headerIndex = 1, obs_numheader(lobsSpaceData)
        codeType     = obs_headElem_i(lobsSpaceData, OBS_ITY, headerIndex)
        bodyIndexBeg = obs_headElem_i(lobsSpaceData, OBS_RLN, headerIndex)
        bodyIndexEnd = obs_headElem_i(lobsSpaceData, OBS_NLV, headerIndex) + bodyIndexBeg - 1
        do bodyIndex = bodyIndexBeg, bodyIndexEnd
          pjo_1 = obs_bodyElem_r(lobsSpaceData, OBS_JOBS, bodyIndex)
          dataset_loop: do SSTdatasetIndex = 1, oer_getSSTdataParam_int('numberSSTDatasets')
            if (codeType == oer_getSSTdataParam_int('codeType', SSTdatasetIndex) .and. codeType /= codtyp_get_codtyp('satob')) then
              joSSTInstrument(SSTdatasetIndex) = joSSTInstrument(SSTdatasetIndex) + pjo_1
              nobsInstrument(SSTdatasetIndex) = nobsInstrument(SSTdatasetIndex) + 1
              exit dataset_loop
            else
              if (obs_elem_c(lobsSpaceData, 'STID', headerIndex) == oer_getSSTdataParam_char('sensor', SSTdatasetIndex)) then
                joSSTInstrument(SSTdatasetIndex) = joSSTInstrument(SSTdatasetIndex) + pjo_1
                nobsInstrument(SSTdatasetIndex) = nobsInstrument(SSTdatasetIndex) + 1
                exit dataset_loop
              end if
            end if
          end do dataset_loop
        end do
      end do
    end if

    call mmpi_allreduce_sumreal8scalar( pjo, "GRID" )
    call mmpi_allreduce_sumreal8scalar( dljoraob, "GRID" )
    call mmpi_allreduce_sumreal8scalar( dljoairep, "GRID" )
    call mmpi_allreduce_sumreal8scalar( dljosatwind, "GRID" )
    call mmpi_allreduce_sumreal8scalar( dljosurfc, "GRID" )
    call mmpi_allreduce_sumreal8scalar( dljoscat, "GRID" )
    call mmpi_allreduce_sumreal8scalar( dljotov, "GRID" )
    call mmpi_allreduce_sumreal8scalar( dljogpsro, "GRID" )
    call mmpi_allreduce_sumreal8scalar( dljoprof, "GRID" )
    call mmpi_allreduce_sumreal8scalar( dljogpsztd, "GRID" )
    call mmpi_allreduce_sumreal8scalar( dljochm, "GRID" )
    call mmpi_allreduce_sumreal8scalar( dljosst, "GRID" )
    call mmpi_allreduce_sumreal8scalar( dljoaladin,"GRID")
    call mmpi_allreduce_sumreal8scalar( dljoice, "GRID" )
    call mmpi_allreduce_sumreal8scalar( dljohydro, "GRID" )
    call mmpi_allreduce_sumreal8scalar( dljoradar, "GRID" )
    do sensorIndex = 1, tvs_nsensors
      call mmpi_allreduce_sumreal8scalar(dljotov_sensors(sensorIndex), "GRID")
    end do
    if (printJoTovsPerChannelSensor) then
      loopSensor2: do sensorIndex = 1, tvs_nsensors
        if (trim(sensorNameList(sensorIndex)) == '') cycle loopSensor2

        call mmpi_allreduce_sumR8_1d(joTovsPerChannelSensor(:,sensorIndex), "GRID")
      end do loopSensor2
    end if

    ! SST data per instrument
    do SSTdatasetIndex = 1, oer_getSSTdataParam_int('numberSSTDatasets')
      call mmpi_allreduce_sumreal8scalar(joSSTInstrument(SSTdatasetIndex), "grid")
      call rpn_comm_allreduce(nobsInstrument(SSTdatasetIndex), nobsInstrumentGlob(SSTdatasetIndex), &
                              1, "mpi_integer", "mpi_sum", "grid", ierr)
    end do

    if ( mmpi_myid == 0 .and. .not. beSilent ) then
      write(*,'(a15,f30.17)') 'Jo(UA)   = ', dljoraob
      write(*,'(a15,f30.17)') 'Jo(AI)   = ', dljoairep
      write(*,'(a15,f30.17)') 'Jo(SF)   = ', dljosurfc
      write(*,'(a15,f30.17)') 'Jo(SC)   = ', dljoscat
      write(*,'(a15,f30.17)') 'Jo(TO)   = ', dljotov
      write(*,'(a15,f30.17)') 'Jo(SW)   = ', dljosatwind
      write(*,'(a15,f30.17)') 'Jo(PR)   = ', dljoprof
      write(*,'(a15,f30.17)') 'Jo(RO)   = ', dljogpsro
      write(*,'(a15,f30.17)') 'Jo(GP)   = ', dljogpsztd
      write(*,'(a15,f30.17)') 'Jo(CH)   = ', dljochm
      write(*,'(a15,f30.17)') 'Jo(TM)   = ', dljosst
      write(*,'(a15,f30.17)') 'Jo(AL)   = ', dljoaladin
      write(*,'(a15,f30.17)') 'Jo(GL)   = ', dljoice
      write(*,'(a15,f30.17)') 'Jo(HY)   = ', dljohydro
      write(*,'(a15,f30.17)') 'Jo(RA)   = ', dljoradar
      write(*,*) ' '
      if ( tvs_nsensors > 0 ) then
        write(*,'(1x,a)') 'For TOVS decomposition by sensor:'
        write(*,'(1x,a)') '#  plt sat ins    Jo'
        do sensorIndex = 1, tvs_nsensors
          write(*,'(i2,1x,a,1x,a,1x,i2,1x,f30.17)') sensorIndex, inst_name(tvs_instruments(sensorIndex)), &
                                                    platform_name(tvs_platforms(sensorIndex)), &
                                                    tvs_satellites(sensorIndex), &
                                                    dljotov_sensors(sensorIndex)
        end do
        write(*,*) ' '
      end if

      ! print per channel information
      if ( tvs_nsensors > 0 .and. printJoTovsPerChannelSensor ) then
        write(*,'(1x,a)') 'For TOVS decomposition by sensor/channel:'
        write(*,'(1x,a)') 'index  sensorName  channel  Jo'
        loopSensor3: do sensorIndex = 1, tvs_nsensors
        if ( trim(sensorNameList(sensorIndex)) == '' ) cycle loopSensor3

          loopChannel: do channelIndex = 1, tvs_maxNumberOfChannels
            if ( channelNumberList(channelIndex,sensorIndex) == 0 ) cycle loopChannel

            write(*,'(i2,1x,a,1x,i4,1x,f30.17)') sensorIndex, &
                                                 sensorNameList(sensorIndex), &
                                                 channelNumberList(channelIndex,sensorIndex), &
                                                 joTovsPerChannelSensor(channelIndex,sensorIndex)
          end do loopChannel

        end do loopSensor3
        write(*,*) ' '
      end if

      ! print SST data per instrument
      if(oer_getSSTdataParam_int('numberSSTDatasets') > 0) then
        write(*,*) 'cfn_sumJo: SST data by data type:'
        write(*,'(a5,a15,a10,a30,a10, a15)') 'index', ' instrument', ' sensor', ' Jo', ' nobs', ' Jo/nobs'
        do SSTdatasetIndex = 1, oer_getSSTdataParam_int('numberSSTDatasets')
          if (nobsInstrumentGlob(SSTdatasetIndex) > 0) then
            write(*,'(i5,a15,a10,f30.17,i15,f10.5)') SSTdatasetIndex, &
                                                     oer_getSSTdataParam_char('instrument', SSTdatasetIndex),&
                                                     oer_getSSTdataParam_char('sensor', SSTdatasetIndex), &
                                                     joSSTInstrument(SSTdatasetIndex), &
                                                     nobsInstrumentGlob(SSTdatasetIndex),&
                                                     joSSTInstrument(SSTdatasetIndex) / &
                                                     real(nobsInstrumentGlob(SSTdatasetIndex))
          end if    
        end do
      end if

    end if
    
    if(oer_getSSTdataParam_int('numberSSTDatasets') > 0) then
      deallocate(joSSTInstrument)
      deallocate(nobsInstrument)
      deallocate(nobsInstrumentGlob)
    end if

  end subroutine cfn_sumJo

  !--------------------------------------------------------------------------
  ! readNameList
  !--------------------------------------------------------------------------
  subroutine readNameList
    !
    ! :Purpose: Reading NAMCFN namelist by any subroutines in costfunction_mod module.
    !
    implicit none

    ! Locals:
    integer :: nulnam, ierr
    integer, external :: fnom, fclos
    logical, save :: nmlAlreadyRead = .false.
    NAMELIST /NAMCFN/ sensorNameList, channelNumberList

    if ( .not. nmlAlreadyRead ) then
      nmlAlreadyRead = .true.

      !- Setting default values
      sensorNameList(:) = ''
      channelNumberList(:,:) = 0

      if ( .not. utl_isNamelistPresent('NAMCFN','./flnml') ) then
        if ( mmpi_myid == 0 ) then
          write(*,*) 'NAMCFN is missing in the namelist. The default values will be taken.'
        end if

      else
        ! Reading the namelist
        nulnam = 0
        ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
        read(nulnam, nml=namcfn, iostat=ierr)
        if ( ierr /= 0) call utl_abort('costfunction_mod: Error reading namelist')
        ierr = fclos(nulnam)

        call sortChannelNumbersInNml
      end if
      if ( mmpi_myid == 0 ) write(*,nml=namcfn)
    end if

  end subroutine readNameList

  !--------------------------------------------------------------------------
  ! sortChannelNumbersInNml
  !--------------------------------------------------------------------------
  subroutine sortChannelNumbersInNml
    !
    !:Purpose: Sort channelNumbers in NAMCFN namelist. This involves removing
    !          the duplicates and combine channelNumbers of same sensor
    !          prescribed on different lines.
    !
    implicit none

    ! Locals:
    integer :: channelIndex
    integer :: channelIndex1, channelIndex2
    integer :: sensorIndexInList, sensorIndexInList1, sensorIndexInList2

    character(len=15) :: sensorName1LowerCase, sensorName2LowerCase

    if ( mmpi_myid == 0 ) then
      write(*,*) 'costfunction_mod: sortChannelNumbersInNml START'
    end if

    ! set duplicate channelNumber for each sensor to zero
    loopSensor3: do sensorIndexInList = 1, tvs_nsensors
      if ( trim(sensorNameList(sensorIndexInList)) == '' ) cycle loopSensor3

      loopChannel2: do channelIndex = tvs_maxNumberOfChannels, 2, -1
        if ( channelNumberList(channelIndex,sensorIndexInList) == 0 ) cycle loopChannel2

        if ( any(channelNumberList(1:channelIndex-1,sensorIndexInList) == &
                 channelNumberList(channelIndex    ,sensorIndexInList)) ) then
          channelNumberList(channelIndex,sensorIndexInList) = 0
        end if
      end do loopChannel2
    end do loopSensor3

    ! remove later duplicates of sensor/channel pairs
    loopSensor4: do sensorIndexInList1 = 1, tvs_nsensors-1
      if ( trim(sensorNameList(sensorIndexInList1)) == '' ) cycle loopSensor4

      call up2low(sensorNameList(sensorIndexInList1),sensorName1LowerCase)

      loopSensor5: do sensorIndexInList2 = sensorIndexInList1+1, tvs_nsensors
        if ( trim(sensorNameList(sensorIndexInList2)) == '' ) cycle loopSensor5

        call up2low(sensorNameList(sensorIndexInList2),sensorName2LowerCase)

        if ( trim(sensorName2LowerCase) == trim(sensorName1LowerCase) ) then
          loopChannel3: do channelIndex2 = 1, tvs_maxNumberOfChannels
            if ( channelNumberList(channelIndex2,sensorIndexInList2) == 0 ) cycle loopChannel3

            if ( any(channelNumberList(:           ,sensorIndexInList1) == &
                     channelNumberList(channelIndex2,sensorIndexInList2)) ) then
              ! set second occurance of channelNumber to zero
              channelNumberList(channelIndex2,sensorIndexInList2) = 0
            else

              ! replace the first zero channelNumber of first sensor with new non-zero channelNumber
              do channelIndex1 = 1, tvs_maxNumberOfChannels
                if ( channelNumberList(channelIndex1,sensorIndexInList1) == 0 ) then
                  channelNumberList(channelIndex1,sensorIndexInList1) = &
                        channelNumberList(channelIndex2,sensorIndexInList2)
                  channelNumberList(channelIndex2,sensorIndexInList2) = 0
                  exit
                end if
              end do

            end if

          end do loopChannel3
        end if
      end do loopSensor5
    end do loopSensor4

    ! if all entries for a sensor are zero, remove the sensor from namelist.
    ! otherwise sort in ascending order
    loopSensor6: do sensorIndexInList = 1, tvs_nsensors
      if ( trim(sensorNameList(sensorIndexInList)) == '' ) cycle loopSensor6

      if ( all(channelNumberList(:,sensorIndexInList) == 0) ) then
        sensorNameList(sensorIndexInList) = ''
      else
        call isort(channelNumberList(:,sensorIndexInList),tvs_maxNumberOfChannels)
      end if
    end do loopSensor6

    if ( mmpi_myid == 0 ) then
      write(*,*) 'costfunction_mod: sortChannelNumbersInNml END'
    end if

  end subroutine sortChannelNumbersInNml

end module costFunction_mod
