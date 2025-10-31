
module costFunction_mod
  ! MODULE costFunction_mod, (prefix='cfn' category='5. Observation operators')
  !
  !:Purpose: To compute Jo term.
  !
  use midasMpi_mod
  use obsSpaceData_mod
  use utilities_mod
  use obserrors_mod
  use codtyp_mod

  implicit none
  save
  private

  ! Public procedures
  public :: cfn_calcJo, cfn_sumJo

  integer,           allocatable :: channelNumberList(:,:)
  character(len=15), allocatable :: sensorNameList(:)
  logical                        :: allReduceForward

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
    integer :: bodyIndex, sensorIndex, headerIndex, bodyIndexBeg, bodyIndexEnd
    integer :: channelNumber, channelNumberIndexInListFound, channelIndex
    integer :: sensorIndexInList, sensorIndexInListFound
    logical :: beSilent

    real(8) :: dljoraob, dljoairep, dljosatwind, dljoscat, dljosurfc, dljotov, dljosst, dljoice
    real(8) :: dljoprof, dljogpsro, dljogpsztd, dljochm, pjo_1, dljoaladin, dljohydro, dljoradar
    character(len=15) :: lowerCaseName

    logical :: printJoTovsPerChannelSensor

    real(8), allocatable :: joSSTInstrument(:)
    integer, allocatable :: nobsInstrument(:), nobsInstrumentGlob(:)
    integer :: SSTdatasetIndex, codeType

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .false.
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
    if(oer_getSSTdataParam_int('numberSSTDatasets') > 0) then
      joSSTInstrument(:) = 0.0d0
      nobsInstrumentGlob(:) = 0
      nobsInstrument(:) = 0
    end if

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

    call mmpi_allreduce_sumreal8scalar(pjo,         allReduceForward_opt = allReduceForward)
    call mmpi_allreduce_sumreal8scalar(dljoraob,    allReduceForward_opt = allReduceForward)
    call mmpi_allreduce_sumreal8scalar(dljoairep,   allReduceForward_opt = allReduceForward)
    call mmpi_allreduce_sumreal8scalar(dljosatwind, allReduceForward_opt = allReduceForward)
    call mmpi_allreduce_sumreal8scalar(dljosurfc,   allReduceForward_opt = allReduceForward)
    call mmpi_allreduce_sumreal8scalar(dljoscat,    allReduceForward_opt = allReduceForward)
    call mmpi_allreduce_sumreal8scalar(dljotov,     allReduceForward_opt = allReduceForward)
    call mmpi_allreduce_sumreal8scalar(dljogpsro,   allReduceForward_opt = allReduceForward)
    call mmpi_allreduce_sumreal8scalar(dljoprof,    allReduceForward_opt = allReduceForward)
    call mmpi_allreduce_sumreal8scalar(dljogpsztd,  allReduceForward_opt = allReduceForward)
    call mmpi_allreduce_sumreal8scalar(dljochm,     allReduceForward_opt = allReduceForward)
    call mmpi_allreduce_sumreal8scalar(dljosst,     allReduceForward_opt = allReduceForward)
    call mmpi_allreduce_sumreal8scalar(dljoaladin,  allReduceForward_opt = allReduceForward)
    call mmpi_allreduce_sumreal8scalar(dljoice,     allReduceForward_opt = allReduceForward)
    call mmpi_allreduce_sumreal8scalar(dljohydro,   allReduceForward_opt = allReduceForward)
    call mmpi_allreduce_sumreal8scalar(dljoradar,   allReduceForward_opt = allReduceForward)

    ! SST data per instrument
    do SSTdatasetIndex = 1, oer_getSSTdataParam_int('numberSSTDatasets')
      call mmpi_allreduce_sumreal8scalar(joSSTInstrument(SSTdatasetIndex), allReduceForward_opt = allReduceForward)
      call mmpi_allReduce(nobsInstrument(SSTdatasetIndex), nobsInstrumentGlob(SSTdatasetIndex), mmpi_sum)
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

      ! print SST data per instrument
      if(oer_getSSTdataParam_int('numberSSTDatasets') > 0) then
        write(*,*) 'cfn_sumJo: SST data by data type:'
        write(*,'(a10, a15, a10, a30, a20, a20)') 'index', ' instrument', ' sensor', ' Jo', ' nobs', ' Jo/nobs'
        do SSTdatasetIndex = 1, oer_getSSTdataParam_int('numberSSTDatasets')
          if (nobsInstrumentGlob(SSTdatasetIndex) > 0) then
            write(*,'(i10, a15, a10, f30.17,i20, f20.5)') SSTdatasetIndex, &
                                                          '      '//oer_getSSTdataParam_char('instrument', SSTdatasetIndex),&
                                                          '    '//oer_getSSTdataParam_char('sensor', SSTdatasetIndex), &
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
    integer :: ierr
    logical, save :: nmlAlreadyRead = .false.
    NAMELIST /NAMCFN/ sensorNameList, channelNumberList, allReduceForward

    if ( .not. nmlAlreadyRead ) then
      nmlAlreadyRead = .true.

      !- Setting default values
      sensorNameList(:) = ''
      channelNumberList(:,:) = 0
      allReduceForward = .true.

      if ( .not. utl_isNamelistPresent('NAMCFN','./flnml') ) then
        if ( mmpi_myid == 0 ) then
          write(*,*) 'NAMCFN is missing in the namelist. The default values will be taken.'
        end if

      else
        ! Reading the namelist
        call utl_tmg_start(181,'low-level--readNML')
        read(utl_flnml, nml=namcfn, iostat=ierr)
        if ( ierr /= 0) call utl_abort('costfunction_mod: Error reading namelist')
        call utl_tmg_stop(181)

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

    if ( mmpi_myid == 0 ) then
      write(*,*) 'costfunction_mod: sortChannelNumbersInNml END'
    end if

  end subroutine sortChannelNumbersInNml

end module costFunction_mod
