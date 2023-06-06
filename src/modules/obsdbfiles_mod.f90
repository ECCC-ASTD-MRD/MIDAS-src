
module obsdbFiles_mod
  ! MODULE obsdbFiles (prefix='odbf' category='3. Observation input/output')
  !
  ! :Purpose: To read and update sqlite files that are in the new 'obsDB' format.
  !
  use midasMpi_mod
  use codePrecision_mod
  use mathPhysConstants_mod
  use bufr_mod
  use codtyp_mod
  use obsSpaceData_mod
  use fSQLite
  use utilities_mod
  use clib_interfaces_mod
  use obsUtil_mod
  use obsVariableTransforms_mod
  use sqliteUtilities_mod
  use timeCoord_mod

  implicit none
  save
  private

  ! Public subroutines and functions:
  public :: odbf_getDateStamp, odbf_readFile, odbf_updateFile, obdf_clean

  ! Arrays used to match obsDB column names with obsSpaceData column names

  integer, parameter :: lenSqlName    = 60
  integer, parameter :: sqlColIndex   = 1
  integer, parameter :: obsColIndex   = 2
  integer, parameter :: varNoColIndex = 2
  integer, parameter :: maxElementnumber = 100
  
  character(len=lenSqlName) :: headTableName = ' '
  character(len=lenSqlName) :: bodyTableName = ' '
  character(len=lenSqlName) :: midasBodyTableName = ' '
  character(len=lenSqlName) :: midasHeadTableName = ' '

  ! ...for the header table

  integer :: numHeadMatch 
  character(len=lenSqlName) :: obsHeadKeySqlName = ' '
  character(len=lenSqlName) :: headDateSqlName = ' ' 
  character(len=lenSqlName), allocatable :: headMatchList(:,:) 
  character(len=6), allocatable :: headBufrList(:)

  ! ...for the body table
  
  integer :: numBodyMatch
  character(len=lenSqlName) :: obsBodyKeySqlName = ' '
  character(len=lenSqlName), allocatable :: bodyMatchList(:,:)
  character(len=6), allocatable :: bodyBufrList(:)

  ! Dictionary of 'varno' value for each obsDB observation value column
  integer :: numVarNo
  character(len=lenSqlName), allocatable :: varNoList(:,:)

  ! Column names for the MIDAS Header table and corresponding obsSpace names
  character(len=lenSqlName) :: midasHeadKeySqlName
  integer :: numMidasHeadMatch
  character(len=lenSqlName), allocatable :: midasHeadNamesList(:,:)

  ! Column names for the MIDAS Body table and corresponding obsSpace names
  character(len=lenSqlName) :: midasBodyKeySqlName
  integer :: numMidasBodyMatch
  character(len=lenSqlName), allocatable :: midasBodyNamesList(:,:)  
  integer :: numBodyMidasTableRequired
   
  ! Other constants
  logical, parameter :: setObsFlagZero = .true.
  character(len=20)  :: combinedTableName = 'combinedTable'

  ! NAMELIST variables
  integer :: numElemIdList                ! MUST NOT BE INCLUDED IN NAMELIST!
  integer :: elemIdList(maxElementNumber) ! list of bufr element IDs to read from file

contains

  !--------------------------------------------------------------------------
  ! readNml
  !--------------------------------------------------------------------------
  subroutine readNml()
    !
    ! :Purpose: Read the namelist for obsDB files
    !
    implicit none

    ! locals
    integer            :: nulfile, ierr
    integer, external  :: fnom, fclos
    logical, save      :: alreadyRead = .false.
    integer            :: elementIndex

    namelist /namobsdb/ numElemIdList, elemIdList

    if ( alreadyRead ) return

    alreadyRead = .true.

    ! default values
    numElemIdList = MPC_missingValue_INT
    elemIdList(:) = 0

    if ( .not. utl_isNamelistPresent('NAMOBSDB','./flnml') ) then
      if ( mmpi_myid == 0 ) then
        write(*,*) 'readNml (odbf): namObsDB is missing in the namelist.'
        write(*,*) '                The default values will be taken.'
      end if
    else
      ! reading namelist variables
      nulfile = 0
      ierr = fnom(nulfile,'./flnml','FTN+SEQ+R/O',0)
      read(nulfile, nml=namobsdb, iostat=ierr)
      if ( ierr /= 0 ) call utl_abort('readNml (odbf): Error reading namelist')
      ierr = fclos(nulfile)
      if (numElemIdList /= MPC_missingValue_INT) then
        call utl_abort('readNml (odbf): check namobsdb namelist section: numElemIdList should be removed')
      end if
      numElemIdList = 0
      do elementIndex = 1, maxElementNumber
        if (elemIdList(elementIndex) == 0) exit
        numElemIdList = numElemIdList + 1
      end do
    end if
    if ( mmpi_myid == 0 ) write(*, nml=namObsDb)

    if (numElemIdList==0) then
      call utl_abort('odbf_setup: element list is empty')
    end if

  end subroutine readNml

  !--------------------------------------------------------------------------
  ! odbf_setup
  !--------------------------------------------------------------------------
  subroutine odbf_setup()
    !
    ! :Purpose: Read the namelist for obsDB files and read the obsDB column table
    !
    implicit none

    ! locals
    
    integer            :: nulfile, ierr
    integer, external  :: fnom, fclos
    logical, save      :: alreadyRead = .false.
    
    character(len=512), parameter :: obsDbColumnFile = 'ObsDBColumnTable.dat' 

    character(len=512) :: readLine
    character(len=lenSqlName) :: readDBColumn, readObsSpaceColumn, readBufrColumn  ! Strings to temporary assign 
                                                                                ! read header and body column names
    
    integer            :: headerTableRow, bodyTableRow, midasBodyTableRow, midasHeadTableRow ! Counters to determine the size of 
                                                                         ! header and body table
    
    integer            :: countRow, countMatchRow, countVarRow                        
    integer            :: stringFoundIndex
    logical            :: obsColumnIsValid
                                                                       
    if ( alreadyRead ) return

    alreadyRead = .true.

    call readNml()

    ! initialize obsDb columns names to be consistent with MIDAS obsSpaceData Report column names

    nulfile = 0 
    ierr = fnom(nulfile,trim(obsDbColumnFile),'FTN+SEQ+R/O',0)
    if ( ierr /= 0 ) call utl_abort('odbf_setup: Error reading ObsDBColumnTable file') 
  
    ! Get number of rows in ObsDBColumnTable file
    headerTableRow = 0
    bodyTableRow = 0
    midasBodyTableRow = 0
    midasHeadTableRow = 0
    numHeadMatch = 0
    numBodyMatch = 0 
    numMidasBodyMatch = 0
    numMidasHeadMatch = 0
    numVarNo = 0
    do while(ierr == 0)
      read(nulfile, '(A)' ,iostat = ierr) readline

      if (index(trim(readline),'OBS_HEADER TABLE INFO BEGINS') > 0) then
        ! Found the start of header table in the file
        do while (ierr == 0)
          read(nulfile,'(A)' ,iostat = ierr) readline
          ! Stop search at the end of header table
          if (index(trim(readline),'OBS_HEADER TABLE INFO ENDS') > 0) exit
          headerTableRow = headerTableRow +1
          if (index(trim(readline),'[') == 0 .and. index(trim(readline),']') == 0) then 
            ! Count obsSpaceData header columns  
            numHeadMatch = numHeadMatch + 1
          end if
        end do

      else if (index(trim(readline),'OBS_BODY TABLE INFO BEGINS') > 0) then  
        ! Found the start of the body table
        do while (ierr == 0)
          read(nulfile,'(A)' ,iostat = ierr) readline
          ! Stop search at the end of body table
          if (index(trim(readline),'OBS_BODY TABLE INFO ENDS') > 0) exit
          bodyTableRow = bodyTableRow + 1
          if (index(trim(readline),'[') == 0 .and. index(trim(readline),']') == 0) then 
            ! Count obsSpaceData body columns
            numBodyMatch = numBodyMatch + 1 
            if (index(trim(readline),'VAR') /= 0) then
              ! Count VarNo columns
              numVarNo = numVarNo + 1
            end if
          end if
        end do

      else if (index(trim(readline),'MIDAS_HEADER TABLE INFO BEGINS') > 0) then
        ! Found the start of the MIDAS_HEADER table
        do while (ierr == 0)
          read(nulfile,'(A)' ,iostat = ierr) readline
          ! Stop search at the end of MIDAS_HEADER table
          if (index(trim(readline),'MIDAS_HEADER TABLE INFO ENDS') > 0) exit
          midasHeadTableRow = midasHeadTableRow + 1
          ! Count obsSpaceData  MIDAS_HEADER columns
          if (index(trim(readline),'[') == 0 .and. index(trim(readline),']') == 0) then 
            numMidasHeadMatch = numMidasHeadMatch + 1
          end if
        end do

      else if (index(trim(readline),'MIDAS_BODY TABLE INFO BEGINS') > 0) then
        ! Found the start of the MIDAS_BODY table
        do while (ierr == 0)
          read(nulfile,'(A)' ,iostat = ierr) readline
          ! Stop search at the end of MIDAS_BODY table
          if (index(trim(readline),'MIDAS_BODY TABLE INFO ENDS') > 0) exit
          midasBodyTableRow = midasBodyTableRow + 1
          ! Count obsSpaceData MIDAS_BODY columns
          if (index(trim(readline),'[') == 0 .and. index(trim(readline),']') == 0) then 
            numMidasBodyMatch = numMidasBodyMatch + 1
          end if
        end do
      end if ! if (index(trim(readline)
    end do ! do while(ierr == 0)

    write(*,*) 'odbf_setup: Number of obs header columns in ObsDBColumnTable', headerTableRow
    write(*,*) 'odbf_setup: Number of obs body columns in ObsDBColumnTable', bodyTableRow
    write(*,*) 'odbf_setup: Number of midas header columns in ObsDBColumnTable', midasHeadTableRow
    write(*,*) 'odbf_setup: Number of midas body columns in ObsDBColumnTable', midasBodyTableRow
    
    write(*,*) 'odbf_setup: numHeadMatch', numHeadMatch
    write(*,*) 'odbf_setup: numBodyMatch', numBodyMatch
    write(*,*) 'odbf_setup: numMidasHeadMatch',numMidasHeadMatch
    write(*,*) 'odbf_setup: numMidasBodyMatch', numMidasBodyMatch

    ! Read Report obsdb and MIDAS obsSpaceData columns name
    allocate(headMatchList(2,numHeadMatch))
    allocate(headBufrList(numHeadMatch))
    allocate(bodyMatchList(2,numBodyMatch))
    allocate(bodyBufrList(numBodyMatch))
    allocate(midasBodyNamesList(2,numMidasBodyMatch))
    allocate(midasHeadNamesList(2,numMidasHeadMatch))
    allocate(varNoList(2,numVarNo))

    rewind(nulfile)

    ierr = 0
    do while (ierr == 0)
      read(nulfile, '(A)' ,iostat = ierr) readline

      if (index(trim(readline),'OBS_HEADER TABLE INFO BEGINS') > 0) then
        ! Found the start of Obs Header table in the file  
        countMatchRow = 0
        do countRow = 1, headerTableRow
          ! Read all header table rows
          read(nulfile,*,iostat = ierr) readDBColumn, readObsSpaceColumn, readBufrColumn      
          select case (trim(readObsSpaceColumn))
            ! Assign read values to appropriate variable
            case ('[headTableName]')
              headTableName = trim(readDBColumn)
              ! Remove the brackets from column name string 
              headTableName = headTableName(1:len(headTableName)-1)
            case ('[headPrimaryKey]')
              obsHeadKeySqlName = trim(readDBColumn)
              ! Remove the brackets from column name string 
              obsHeadKeySqlName = obsHeadKeySqlName(1:len(obsHeadKeySqlName)-1)
            case ('[headDateSqlName]')
              headDateSqlName = trim(readDBColumn)
              ! Remove the brackets from column name string  
              headDateSqlName = headDateSqlName(1:len(headDateSqlName)-1)
            case default
              countMatchRow = countMatchRow + 1
              headMatchList(1,countMatchRow) = trim(readDBColumn)
              headMatchList(2,countMatchRow) = trim(readObsSpaceColumn)
              headBufrList(countMatchRow) = trim(readBufrColumn)
          end select
        end do ! do countRow

      else if (index(trim(readline),'OBS_BODY TABLE INFO BEGINS') > 0) then
        ! Found the start of Obs BODY table in the file  
        countMatchRow = 0
        countVarRow = 0 
        do countRow = 1, bodyTableRow
          ! Read all body table rows 
          read(nulfile,*,iostat = ierr) readDBColumn, readObsSpaceColumn, readBufrColumn      
          select case (trim(readObsSpaceColumn))
            ! Assign read values to appropriate variable
            case ('[bodyTableName]')
              bodyTableName = trim(readDBColumn)
              ! Remove the brackets from column name string
              bodyTableName = bodyTableName(1:len(bodyTableName)-1)
            case ('[bodyPrimaryKey]')
              obsBodyKeySqlName = trim(readDBColumn)
              ! Remove the brackets from column name string
              obsBodyKeySqlName = obsBodyKeySqlName(1:len(obsBodyKeySqlName)-1)
            case default
              countMatchRow = countMatchRow + 1
              bodyMatchList(1,countMatchRow) = trim(readDBColumn)
              bodyMatchList(2,countMatchRow) = trim(readObsSpaceColumn)
              bodyBufrList(countMatchRow) = trim(readBufrColumn)

              ! Assign varNoList if appropriate
              if (bodyMatchList(2,countMatchRow) == 'VAR') then
                countVarRow = countVarRow + 1 
                varNoList(1, countVarRow) = bodyMatchList(1,countMatchRow) 
                varNoList(2, countVarRow) = bodyBufrList(countMatchRow)
              end if
          end select
        end do ! do countRow

      else if (index(trim(readline),'MIDAS_HEADER TABLE INFO BEGINS') > 0) then
        ! Found the start of MIDAS Header table in the file
        countMatchRow = 0
        do countRow = 1, midasHeadTableRow
          !Read all header table rows
          read(nulfile,*,iostat = ierr) readDBColumn, readObsSpaceColumn, readBufrColumn
          select case (trim(readObsSpaceColumn))
            case ('[midasHeadTableName]')
              midasHeadTableName = trim(readDBColumn)
              midasHeadTableName = midasHeadTableName(1:len(midasHeadTableName)-1)
            case ('[midasHeadPrimaryKey]')
              midasHeadKeySqlName = trim(readDBColumn)
              ! Remove the brackets from column name string
              midasHeadKeySqlName = midasHeadKeySqlName(1:len(midasHeadKeySqlName)-1)
            case default
              countMatchRow = countMatchRow + 1
              midasHeadNamesList(1,countMatchRow) = trim(readDBColumn)              
              midasHeadNamesList(2,countMatchRow) = trim(readObsSpaceColumn)
          end select
        end do ! do countRow

      else if (index(trim(readline),'MIDAS_BODY TABLE INFO BEGINS') > 0) then
        ! Found the start of MIDAS Body table in the file
        numBodyMidasTableRequired = 0
        countMatchRow = 0
        do countRow = 1, midasBodyTableRow
          ! Read all body table rows
          read(nulfile,*,iostat = ierr) readDBColumn, readObsSpaceColumn, readBufrColumn      
          select case (trim(readObsSpaceColumn))
            ! Assign read values to appropriate variable
            case ('[midasBodyTableName]')
              midasBodyTableName = trim(readDBColumn)
              ! Remove the brackets from column name string
              midasBodyTableName = midasBodyTableName(1:len(midasBodyTableName)-1)
            case ('[midasBodyPrimaryKey]')
              midasBodyKeySqlName = trim(readDBColumn)
              ! Remove the brackets from column name string
              midasBodyKeySqlName = midasBodyKeySqlName(1:len(midasBodyKeySqlName)-1)
            case default
              countMatchRow = countMatchRow + 1
              midasBodyNamesList(1,countMatchRow) = trim(readDBColumn)              
              midasBodyNamesList(2,countMatchRow) = trim(readObsSpaceColumn)
             
              ! Check if the ObsSpaceData column name is mandatory 
              stringFoundIndex = index(midasBodyNamesList(2,countMatchRow),'*')
              if (stringFoundIndex > 0) then
                midasBodyNamesList(2,countMatchRow) = trim(midasBodyNamesList(2,countMatchRow)(1:stringFoundIndex-1))
                numBodyMidasTableRequired = numBodyMidasTableRequired + 1
              end if
          end select
        end do ! do countRow
      end if ! if (index(trim(readline)
    end do ! do while (ierr == 0)

    ierr = fclos(nulfile)
 
    ! Check if the header, body and MIDAS table column names are read correctly 
    if (len(trim(headTableName)) == 0) call utl_abort('odbf_setup: headTableName is incorrectly defined or missing')
    if (len(trim(obsHeadKeySqlName)) == 0) call utl_abort('odbf_setup: obsHeadKeySqlName is incorrectly defined or missing')
    if (len(trim(headDateSqlName)) == 0) call utl_abort('odbf_setup: headDateSqlName is incorrectly defined or missing')
    if (len(trim(bodyTableName)) == 0) call utl_abort('odbf_setup: bodyTableName is incorrectly defined or missing')
    if (len(trim(obsBodyKeySqlName)) == 0) call utl_abort('odbf_setup: obsBodyKeySqlName is incorrectly defined or missing')
    if (len(trim(midasBodyTableName)) == 0) call utl_abort('odbf_setup: midasBodyTableName is incorrectly defined or missing')
    if (len(trim(midasBodyKeySqlName)) == 0) call utl_abort('odbf_setup: midasBodyKeySqlName is incorrectly defined or missing')
    if (len(trim(midasHeadKeySqlName)) == 0) call utl_abort('odbf_setup: midasHeadKeySqlName is incorrectly defined or missing')

    do countRow = 1, numHeadMatch
      obsColumnIsValid = obs_isColumnNameValid(trim(headMatchList(2,countRow)))
      if (.not. obsColumnIsValid) then
        call utl_abort('odbf_setup: Column in Header Table does not exist in ObsSpaceData: ' // &
                       trim(headMatchList(2,countRow)))
      end if
    end do
      
    do countRow = 1, numBodyMatch 
      obsColumnIsValid = obs_isColumnNameValid(trim(bodyMatchList(2,countRow)))
      if (.not. obsColumnIsValid) then
        call utl_abort('odbf_setup: Column in Body Table does not exist in ObsSpaceData: ' // &
                       trim(bodyMatchList(2,countRow)))
      end if
    end do

    do countRow = 1, numMidasBodyMatch 
      obsColumnIsValid = obs_isColumnNameValid(trim(midasBodyNamesList(2,countRow)))
      if (.not. obsColumnIsValid) then
        call utl_abort('odbf_setup: Column in MIDAS Body Table does not exist in ObsSpaceData: ' // &
                       trim(midasBodyNamesList(2,countRow)))
      end if
    end do

    do countRow = 1, numMidasHeadMatch 
      obsColumnIsValid = obs_isColumnNameValid(trim(midasHeadNamesList(2,countRow)))
      if (.not. obsColumnIsValid) then
        call utl_abort('odbf_setup: Column in MIDAS Header Table does not exist in ObsSpaceData: ' // &
                       trim(midasHeadNamesList(2,countRow)))
      end if
    end do

  end subroutine odbf_setup

  !--------------------------------------------------------------------------
  ! odbf_getDateStamp
  !--------------------------------------------------------------------------
  subroutine odbf_getDateStamp(dateStamp, fileName)
    !
    ! Purpose: get dateStamp from an obsDB file
    !
    implicit none
    
    ! arguments
    integer         , intent(out) :: dateStamp
    character(len=*), intent(in)  :: fileName
    
    ! locals
    integer,    allocatable :: headDateValues(:), headTimeValues(:)
    integer                 :: ier, imode, validTime, validDate, validDateRecv, validTimeRecv
    integer                 :: newdate

    call odbf_setup()

    call sqlu_getColumnValuesDateStr(headDateValues, headTimeValues, fileName=trim(fileName), &
                                     tableName=headTableName, sqlColumnName=headDateSqlName)

    validDate = MPC_missingValue_INT 
    validTime = MPC_missingValue_INT 

    call tim_getValidDateTimeFromList(headDateValues, headTimeValues, validDate, validTime)

    ! Make sure all mpi tasks have a valid date (important for split sqlite files)
    call rpn_comm_allreduce(validDate, validDateRecv, 1, "MPI_INTEGER", "MPI_MAX", "GRID", ier)
    call rpn_comm_allreduce(validTime, validTimeRecv, 1, "MPI_INTEGER", "MPI_MAX", "GRID", ier)
    
    if (validDateRecv == MPC_missingValue_INT .or. validTimeRecv == MPC_missingValue_INT) then
      call utl_abort('odbf_getDateStamp: Error in getting valid date and time!')
    end if
    
    ! printable to stamp, validTime must be multiplied with 1e6 to make newdate work
    imode = 3
    ier = newdate(dateStamp, validDateRecv, validTimeRecv * 1000000, imode)
    write(*,*)'odbf_getDateStamp: obsDB files valid date (YYYYMMDD): ', validDateRecv
    write(*,*)'odbf_getDateStamp: obsDB files valid time       (HH): ', validTimeRecv
    write(*,*)'odbf_getDateStamp: obsDB files dateStamp            : ', datestamp

    deallocate(headDateValues)
    deallocate(headTimeValues)

  end subroutine odbf_getDateStamp

  !--------------------------------------------------------------------------
  ! odbf_readFile
  !--------------------------------------------------------------------------
  subroutine odbf_readFile(obsdat, fileName, familyType, fileIndex)
    !
    ! :Purpose: Read the contents of an obsDB file and put in obsSpaceData
    !
    implicit none

    ! arguments:
    type (struct_obs), intent(inout) :: obsdat
    character(len=*)                 :: fileName
    character(len=*)                 :: familyType
    integer                          :: fileIndex

    ! locals:
    integer :: bodyIndex, bodyIndexBegin, bodyIndexEnd, headIndexBegin
    integer :: headIndexEnd, headIndex, obsRln
    integer :: numBody, numHead, columnIndex
    integer :: headTableIndex, numRowsHeadTable, bodyTableIndex, numRowsBodyTable
    character(len=lenSqlName), allocatable :: headCharSqlNames(:)
    character(len=lenSqlName), allocatable :: headSqlNames(:),    bodySqlNames(:)
    character(len=50),         allocatable :: headCharValues(:,:)
    integer,                   allocatable :: headDateValues(:), headTimeValues(:)
    real(8),                   allocatable :: headValues(:,:), bodyValues(:,:)
    integer(8),                allocatable :: headPrimaryKey(:)
    integer(8),                allocatable :: bodyPrimaryKey(:), bodyHeadKey(:)

    write(*,*)
    write(*,*) 'odbf_readFile: Starting'
    write(*,*)
    write(*,*) 'odbf_readFile: FileName   : ', trim(FileName)
    write(*,*) 'odbf_readFile: FamilyType : ', FamilyType

    !- 0.0 Some initialization
    call odbf_setup()
    call ovt_setup(elemIdList(1:numElemIdList))

    !- 1.0 Determine names of columns present in obsDB file

    call sqlu_getSqlColumnNames(headCharSqlNames, fileName=trim(fileName), &
                                tableName=headTableName, dataType='varchar')
    call sqlu_getSqlColumnNames(headSqlNames, fileName=trim(fileName), &
                                tableName=headTableName, dataType='numeric' )
    call sqlu_getSqlColumnNames(bodySqlNames, fileName=trim(fileName), &
                                tableName=bodyTableName, dataType='numeric')

    ! Print all of the column names to the listing
    do columnIndex = 1, size(headCharSqlNames)
      write(*,*) 'odbf_readFile: headCharSqlNames =', columnIndex, &
                 trim(headCharSqlNames(columnIndex))
    end do
    do columnIndex = 1, size(headSqlNames)
      write(*,*) 'odbf_readFile: headSqlNames =', columnIndex, &
                 trim(headSqlNames(columnIndex))
    end do
    do columnIndex = 1, size(bodySqlNames)
      write(*,*) 'odbf_readFile: bodySqlNames   =', columnIndex, &
                 trim(bodySqlNames(columnIndex))
    end do

    !- 1.1 Read most of the contents of the file into local tables

    call odbf_getPrimaryKeys(headPrimaryKey, bodyPrimaryKey, bodyHeadKey, &
                             fileName=trim(fileName))

    call sqlu_getColumnValuesDateStr(headDateValues, headTimeValues, fileName=trim(fileName), &
                                     tableName=headTableName, sqlColumnName=headDateSqlName)
    call sqlu_getColumnValuesChar(headCharValues, fileName=trim(fileName), &
                                  tableName=headTableName, sqlColumnNames=headCharSqlNames)
    call sqlu_getColumnValuesNum (headValues, fileName=trim(fileName), &
                                  tableName=headTableName, sqlColumnNames=headSqlNames)
    call sqlu_getColumnValuesNum (bodyValues, fileName=trim(fileName), &
                                  tableName=bodyTableName, sqlColumnNames=bodySqlNames)
    numRowsBodyTable = size(bodyValues,1)
    numRowsHeadTable = size(headValues,1)
    
    ! For debugging, print first 10 rows of each local table to the listing
    do headTableIndex = 1, min(10, numRowsHeadTable)
      write(*,*) 'odbf_readFile: headKeyValues  = ', headPrimaryKey(headTableIndex)
    end do
    do headTableIndex = 1, min(10, numRowsHeadTable)
      write(*,*) 'odbf_readFile: headCharValues = ', headCharValues(headTableIndex,:)
    end do
    do headTableIndex = 1, min(10, numRowsHeadTable)
      write(*,*) 'odbf_readFile: headValues = ', headValues(headTableIndex,:)
    end do
    do bodyTableIndex = 1, min(10, numRowsBodyTable)
      write(*,*) 'odbf_readFile: bodyKeyValues  = ', bodyPrimaryKey(bodyTableIndex)
    end do
    do bodyTableIndex = 1, min(10, numRowsBodyTable)
      write(*,*) 'odbf_readFile: bodyHeadKeyValues  = ', bodyHeadKey(bodyTableIndex)
    end do
    do bodyTableIndex = 1, min(10, numRowsBodyTable)
      write(*,*) 'odbf_readFile: bodyValues = ', bodyValues(bodyTableIndex,:)
    end do

    

    ! Starting point for adding rows to obsSpaceData
    bodyIndexBegin = obs_numBody(obsdat) + 1
    headIndexBegin = obs_numHeader(obsdat) + 1

    !- 1.2 Initialize header and body columns as 'missing'

    ! Set a bunch of obsSpaceData header columns to 'missing'
    do headTableIndex = 1, numRowsHeadTable
      headIndex = headTableIndex + headIndexBegin - 1
      call obs_headSet_i(obsdat, OBS_SEN, headIndex, nint(MPC_missingValue_R8))
      call obs_headSet_i(obsdat, OBS_INS, headIndex, nint(MPC_missingValue_R8))
    end do

    ! Set a bunch of obsSpaceData body columns to 'missing'
    do bodyTableIndex = 1, numRowsBodyTable
      bodyIndex = bodyTableIndex + bodyIndexBegin -1
      if (obs_columnActive_RB(obsdat, OBS_OMA)) call obs_bodySet_r(obsdat, OBS_OMA , bodyIndex, obs_missingValue_R)
      if (obs_columnActive_RB(obsdat, OBS_OMA0)) call obs_bodySet_r(obsdat, OBS_OMA0, bodyIndex, obs_missingValue_R)
      if (obs_columnActive_RB(obsdat, OBS_OMP)) call obs_bodySet_r(obsdat, OBS_OMP , bodyIndex, obs_missingValue_R)
      if (obs_columnActive_RB(obsdat, OBS_OMP6)) call obs_bodySet_r(obsdat, OBS_OMP6, bodyIndex, obs_missingValue_R)
      if (obs_columnActive_RB(obsdat, OBS_OER)) call obs_bodySet_r(obsdat, OBS_OER , bodyIndex, obs_missingValue_R)
      if (obs_columnActive_RB(obsdat, OBS_HPHT)) call obs_bodySet_r(obsdat, OBS_HPHT, bodyIndex, obs_missingValue_R)
      if (obs_columnActive_RB(obsdat, OBS_HAHT)) call obs_bodySet_r(obsdat, OBS_HAHT, bodyIndex, obs_missingValue_R)
      if (obs_columnActive_RB(obsdat, OBS_WORK)) call obs_bodySet_r(obsdat, OBS_WORK, bodyIndex, obs_missingValue_R)
      if (obs_columnActive_RB(obsdat, OBS_SIGI)) call obs_bodySet_r(obsdat, OBS_SIGI, bodyIndex, obs_missingValue_R)
      if (obs_columnActive_RB(obsdat, OBS_SIGO)) call obs_bodySet_r(obsdat, OBS_SIGO, bodyIndex, obs_missingValue_R)
      if (obs_columnActive_RB(obsdat, OBS_ZHA )) call obs_bodySet_r(obsdat, OBS_ZHA , bodyIndex, obs_missingValue_R)
      if (obs_columnActive_RB(obsdat, OBS_SEM )) call obs_bodySet_r(obsdat, OBS_SEM , bodyIndex, obs_missingValue_R)
      if (obs_columnActive_RB(obsdat, OBS_BCOR)) call obs_bodySet_r(obsdat, OBS_BCOR, bodyIndex, obs_missingValue_R)
    end do
    
    !- 1.3 Copy values from local tables into obsSpaceData

    ! Set the columns related to surface type
    call odbf_setSurfaceType(obsdat, headIndexBegin, numRowsHeadTable, &
                             fileName=trim(fileName), tableName=headTableName)

    ! Header date/time values
    call odbf_copyToObsSpaceHeadDate(obsdat, headDateValues, headTimeValues, &
                                     headIndexBegin)

    ! Header-Character values
    call odbf_copyToObsSpaceHeadChar(obsdat, headCharSqlNames, headCharValues, &
                                     headIndexBegin)

    ! Header-Numeric values
    call odbf_copyToObsSpaceHead(obsdat, headSqlNames, headPrimaryKey, &
                                 headValues, headIndexBegin)

    ! Body-Numeric values
    call odbf_copyToObsSpaceBody(obsdat, bodySqlNames, bodyPrimaryKey, bodyHeadKey, &
                                 bodyValues, bodyIndexBegin, headIndexBegin)

    ! Get indexes of last rows added to obsSpaceData
    bodyIndexEnd = obs_numBody(obsdat)
    headIndexEnd = obs_numHeader(obsdat)

    !- 1.4 Set some other quantities in obsSpaceData Header table

    do headIndex = headIndexBegin, headIndexEnd
      call obs_headSet_i(obsdat, OBS_ONM, headIndex, headIndex)
      call obs_headSet_i(obsdat, OBS_IDF, headIndex, fileIndex)
      call obs_setFamily(obsdat, trim(familyType), headIndex)
      if ( headIndex == 1 ) then
        call obs_headSet_i(obsdat, OBS_RLN, headIndex, 1 )
      else
        obsRln = obs_headElem_i(obsdat, OBS_RLN, headIndex-1) + &
                 obs_headElem_i(obsdat, OBS_NLV, headIndex-1)
        call obs_headSet_i(obsdat, OBS_RLN, headIndex, obsRln)
      end if
    end do

    !- 1.5 Read values written during previous MIDAS program executions

    call odbf_readMidasBodyTable(obsdat, trim(fileName), familyType, fileIndex)

    !- 2.0 Additional changes to values after they are in obsSpaceData

    call odbf_adjustValues(obsdat, headIndexBegin, headIndexEnd)

    if ( trim(familyType) /= 'TO' ) then
      call ovt_transformObsValues      (obsdat, headIndexBegin, headIndexEnd )
      call ovt_adjustHumGZ             (obsdat, headIndexBegin, headIndexEnd )
      call obsu_computeVertCoordSurfObs(obsdat, headIndexBegin, headIndexEnd )
    end if

    ! For GP family, initialize OBS_OER to element 15032 (ZTD formal error) 
    ! for all ZTD observations (element 15031)
    if ( trim(familyType) == 'GP' ) then
      write(*,*) 'odbf_readFile: Initializing OBS_OER for GB-GPS ZTD to formal error (ele 15032)'
      call obsu_setGbgpsError( obsdat, headIndexBegin, headIndexEnd )
    end if

    numHead = obs_numHeader(obsdat)
    numBody = obs_numBody(obsdat)
    write(*,*) 'odbf_readFile: after reading file, obs_numHeader = ', numHead
    write(*,*) 'odbf_readFile: after reading file, obs_numBody   = ', numBody
    write(*,*)
    write(*,*) 'odbf_readFile: finished'
    write(*,*)

  end subroutine odbf_readFile

  !--------------------------------------------------------------------------
  ! odbf_readMidasBodyTable
  !--------------------------------------------------------------------------
  subroutine odbf_readMidasBodyTable(obsdat, fileName, familyType, fileIndex)
    !
    ! :Purpose: Read values from any column found the MIDAS table, if it
    !           already exists in the file. This will replace any existing
    !           values read from the original obs-DB tables (e.g. the obs
    !           value).
    !
    implicit none

    ! arguments:
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in)    :: fileName   
    character(len=*), intent(in)    :: familyType
    integer,          intent(in)    :: fileIndex

    ! locals:
    type(fSQL_STATUS)    :: stat ! sqlite error status
    type(fSQL_DATABASE)  :: db   ! sqlite file handle
    type(fSQL_STATEMENT) :: stmt ! precompiled sqlite statements
    integer(8)           :: obsIdd
    integer              :: obsIdf, obsVarNo, sqlNameIndex, matdata_int(1,1)
    integer              :: headIndex, bodyIndex, bodyIndexBegin, bodyIndexEnd
    integer              :: obsSpaceColIndex, ierr, numRows, numColumns
    real(8)              :: obsPPP, matdata_r8(1,1)
    character(len=4)     :: obsSpaceColumnName
    character(len=lenSqlName) :: sqlColumnName, vnmSqlName, pppSqlName, varSqlName
    character(len=3000)  :: query
    logical              :: midasBodyTableExists
    logical, allocatable :: midasColumnExists(:)

    write(*,*)
    write(*,*) 'odbf_readMidasBodyTable: Starting'
    write(*,*)
    write(*,*) 'odbf_readMidasBodyTable: FileName   : ', trim(FileName)
    write(*,*) 'odbf_readMidasBodyTable: FamilyType : ', FamilyType

    ! check if midasTable already exists in the file
    midasBodyTableExists = sqlu_sqlTableExists(fileName, midasBodyTableName)

    if (.not. midasBodyTableExists) then
      write(*,*) 'odbf_readMidasBodyTable: MIDAS table not present in file'
      return
    else
      write(*,*) 'odbf_readMidasBodyTable: MIDAS table present in file, will read contents'
    end if

    ! some sql column names
    vnmSqlName = odbf_midasTabColFromObsSpaceName('VNM', midasBodyNamesList)
    pppSqlName = odbf_midasTabColFromObsSpaceName('PPP', midasBodyNamesList)
    varSqlName = odbf_midasTabColFromObsSpaceName('VAR', midasBodyNamesList)

    ! check which columns exist in the MIDAS output table
    allocate(midasColumnExists(numMidasBodyMatch))
    do sqlNameIndex = 1, numMidasBodyMatch
      sqlColumnName = midasBodyNamesList(1,sqlNameIndex)
      midasColumnExists(sqlNameIndex) = sqlu_sqlColumnExists(fileName, midasBodyTableName, sqlColumnName)
    end do

    ! open the obsDB file
    call fSQL_open(db, trim(fileName), stat)
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_readMidasBodyTable: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort('odbf_readMidasBodyTable: fSQL_open')
    end if
    
    ! read the contents of the MIDAS table, one column at a time
    SQLNAME: do sqlNameIndex = 1, numMidasBodyMatch

      ! skip this sql column name if it is not present in the file
      if (.not. midasColumnExists(sqlNameIndex)) cycle SQLNAME

      ! get obsSpaceData column name and index corresponding to sql column
      sqlColumnName      = midasBodyNamesList(1,sqlNameIndex)
      obsSpaceColumnName = midasBodyNamesList(2,sqlNameIndex)
      ierr = clib_toUpper(obsSpaceColumnName)
      obsSpaceColIndex = obs_columnIndexFromName(trim(obsSpaceColumnName))

      ! skip this sql column name if the related obsSpaceData column is not active
      if (obs_columnDataType(obsSpaceColIndex) == 'real') then
        if (.not. obs_columnActive_RB(obsdat, obsSpaceColIndex)) cycle SQLNAME
      else
        if (.not. obs_columnActive_IB(obsdat, obsSpaceColIndex)) cycle SQLNAME
      end if

      write(*,*) 'odbf_readMidasBodyTable: reading midasTable column: ', &
                 trim(sqlColumnName)
      write(*,*) 'odbf_readMidasBodyTable: to update obsSpaceData column: ', &
                 trim(obsSpaceColumnName)

      ! prepare sql update query
      query = 'select ' // trim(sqlColumnName) // &
              ' from ' // trim(midasBodyTableName) // &
              ' where ' // &
              trim(obsBodyKeySqlName) // ' = ? and '   // &
              trim(vnmSqlName)     // ' = ? and '   // &
              trim(pppSqlName)     // ' = ? ;'
      write(*,*) 'odbf_readMidasBodyTable: query ---> ', trim(query)

      call fSQL_prepare(db, query , stmt, stat)
      if ( fSQL_error(stat) /= FSQL_OK ) then
        write(*,*) 'odbf_readMidasBodyTable: fSQL_prepare: ', fSQL_errmsg(stat)
        call utl_abort('odbf_readMidasBodyTable: fSQL_prepare')
      end if

      call fSQL_begin(db)
      HEADER2: do headIndex = 1, obs_numHeader(obsdat)

        obsIdf = obs_headElem_i( obsdat, OBS_IDF, headIndex )
        if ( obsIdf /= fileIndex ) cycle HEADER2

        bodyIndexBegin = obs_headElem_i( obsdat, OBS_RLN, headIndex )
        bodyIndexEnd = bodyIndexBegin + &
                       obs_headElem_i( obsdat, OBS_NLV, headIndex ) - 1

        BODY2: do bodyIndex = bodyIndexBegin, bodyIndexEnd

          ! execute the sql query to select the desired value
          obsIdd  = obs_bodyPrimaryKey(obsdat, bodyIndex)
          call fSQL_bind_param(stmt, PARAM_INDEX=1, INT8_VAR=obsIdd)
          obsVarNo  = obs_bodyElem_i(obsdat, obs_vnm, bodyIndex)
          call fSQL_bind_param(stmt, PARAM_INDEX=2, INT_VAR=obsVarNo)
          obsPPP  = obs_bodyElem_r(obsdat, obs_ppp, bodyIndex)
          call fSQL_bind_param(stmt, PARAM_INDEX=3, REAL8_VAR=obsPPP)
          call fSQL_exec_stmt(stmt)

          ! read the real or integer value
          if (obs_columnDataType(obsSpaceColIndex) == 'real') then
            call fSQL_get_many(stmt, nrows=numRows, ncols=numColumns, &
                               mode=FSQL_REAL8, status=stat)
            if ( fSQL_error(stat) /= FSQL_OK ) then
              write(*,*) 'odbf_readMidasBodyTable: fSQL_get_many: ', fSQL_errmsg(stat)
              call utl_abort('odbf_readMidasBodyTable: problem with fSQL_get_many')
            end if
            if (numRows /= 1 .or. numColumns /= 1) then
              write(*,*) 'odbf_readMidasBodyTable: numRows, numColumns =', numRows, numColumns
              call utl_abort('odbf_readMidasBodyTable: sql query did not return 1 value')
            end if
            call fSQL_fill_matrix(stmt, matdata_r8)
            call obs_bodySet_r(obsdat,obsSpaceColIndex,bodyIndex,matdata_r8(1,1))
          else
            call fSQL_get_many(stmt, nrows=numRows, ncols=numColumns, &
                               mode=FSQL_INT, status=stat)
            if ( fSQL_error(stat) /= FSQL_OK ) then
              write(*,*) 'odbf_readMidasBodyTable: fSQL_get_many: ', fSQL_errmsg(stat)
              call utl_abort('odbf_readMidasBodyTable: problem with fSQL_get_many')
            end if
            if (numRows /= 1 .or. numColumns /= 1) then
              write(*,*) 'odbf_readMidasBodyTable: numRows, numColumns =', numRows, numColumns
              call utl_abort('odbf_readMidasBodyTable: sql query did not return 1 value')
            end if
            call fSQL_fill_matrix(stmt, matdata_int)
            call obs_bodySet_i(obsdat,obsSpaceColIndex,bodyIndex,matdata_int(1,1))
          end if ! if (obs_columnDataType(obsSpaceColIndex) == 'real')

          call fSQL_free_mem(stmt)

        end do BODY2

      end do HEADER2

      call fSQL_finalize(stmt)
      call fSQL_commit(db)

    end do SQLNAME

    ! close the obsDB file
    call fSQL_close(db, stat) 

    deallocate(midasColumnExists)

    write(*,*)
    write(*,*) 'odbf_readMidasBodyTable: finished'
    write(*,*)

  end subroutine odbf_readMidasBodyTable

  !--------------------------------------------------------------------------
  ! odbf_getPrimaryKeys
  !--------------------------------------------------------------------------
  subroutine odbf_getPrimaryKeys(headPrimaryKey, bodyPrimaryKey, bodyHeadKey, &
                                 fileName)
    !
    ! :Purpose: Read the values from obsDB file for the head and body table
    !           primary keys.
    !
    implicit none

    ! arguments:
    integer(8), allocatable, intent(out) :: headPrimaryKey(:)
    integer(8), allocatable, intent(out) :: bodyPrimaryKey(:)
    integer(8), allocatable, intent(out) :: bodyHeadKey(:)
    character(len=*),        intent(in)  :: fileName

    ! locals:
    integer :: numRows, numColumns
    integer(8), allocatable   :: tempHeadKey(:,:), tempBodyKey(:,:)
    character(len=3000)       :: query
    type(fSQL_STATUS)         :: stat ! sqlite error status
    type(fSQL_DATABASE)       :: db   ! sqlite file handle
    type(fSQL_STATEMENT)      :: stmt ! precompiled sqlite statements

    ! open the obsDB file
    call fSQL_open(db, trim(fileName), status=stat)
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getPrimaryKeys: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort('odbf_getPrimaryKeys: fSQL_open')
    end if

    ! build the sqlite query for the HEADER primary key
    query = 'select ' // trim(obsHeadKeySqlName) // ' from ' // &
            trim(headTableName) // ';'
    write(*,*) 'odbf_getPrimaryKeys: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare(db, trim(query), stmt, status=stat)
    ! note: "status" not set when getting integers
    call fSQL_get_many(stmt, nrows=numRows, ncols=numColumns, &
                       mode=FSQL_INT8)
    write(*,*) 'odbf_getPrimaryKeys: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( headPrimaryKey(numRows) )
    allocate( tempHeadKey(numRows,1) )
    call fSQL_fill_matrix(stmt, tempHeadKey)
    headPrimaryKey(:) = tempHeadKey(:,1)
    deallocate(tempHeadKey)

    ! build the sqlite query for the BODY primary key
    query = 'select ' // trim(obsBodyKeySqlName) // ' from ' // &
            trim(bodyTableName) // ';'
    write(*,*) 'odbf_getPrimaryKeys: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare(db, trim(query) , stmt, status=stat)
    ! note: "status" not set when getting integers
    call fSQL_get_many(stmt, nrows=numRows, ncols=numColumns, &
                       mode=FSQL_INT8)
    write(*,*) 'odbf_getPrimaryKeys: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( bodyPrimaryKey(numRows) )
    allocate( tempBodyKey(numRows,1) )
    call fSQL_fill_matrix(stmt, tempBodyKey)
    bodyPrimaryKey(:) = tempBodyKey(:,1)
    deallocate(tempBodyKey)

    ! build the sqlite query for the BODY-HEAD key
    query = 'select ' // trim(obsHeadKeySqlName) // ' from ' // &
            trim(bodyTableName) // ';'
    write(*,*) 'odbf_getPrimaryKeys: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare(db, trim(query) , stmt, status=stat)
    ! note: "status" not set when getting integers
    call fSQL_get_many(stmt, nrows=numRows, ncols=numColumns, &
                       mode=FSQL_INT8)
    write(*,*) 'odbf_getPrimaryKeys: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( bodyHeadKey(numRows) )
    allocate( tempBodyKey(numRows,1) )
    call fSQL_fill_matrix(stmt, tempBodyKey)
    bodyHeadKey(:) = tempBodyKey(:,1)
    deallocate(tempBodyKey)

    ! close the obsDB file
    call fSQL_free_mem(stmt)
    call fSQL_finalize(stmt)
    call fSQL_close(db, stat) 

  end subroutine odbf_getPrimaryKeys

  !--------------------------------------------------------------------------
  ! odbf_setSurfaceType
  !--------------------------------------------------------------------------
  subroutine odbf_setSurfaceType(obsdat, headIndexBegin, numRowsHeadTable, &
                                 fileName, tableName)
    !
    ! :Purpose: Set the surface type based on lat-lon and some external mask files.
    !
    implicit none

    ! arguments:
    type(struct_obs), intent(inout) :: obsdat
    integer,          intent(in)    :: headIndexBegin
    integer,          intent(in)    :: numRowsHeadTable
    character(len=*), intent(in)    :: fileName
    character(len=*), intent(in)    :: tableName


    ! locals:
    integer              :: numRows, numColumns, headTableIndex, headIndex
    integer, allocatable :: columnValues(:,:)
    character(len=3000)  :: query
    type(fSQL_STATUS)    :: stat ! sqlite error status
    type(fSQL_DATABASE)  :: db   ! sqlite file handle
    type(fSQL_STATEMENT) :: stmt ! precompiled sqlite statements

    ! open the obsDB file
    call fSQL_open(db, trim(fileName), status=stat)
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_setSurfaceType: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort('odbf_setSurfaceType: fSQL_open')
    end if

    ! build the sqlite query
    query = 'select mask_mer(lat,lon), mask_glace_clim(lat,lon) from '&
             // trim(tableName) // ';'
    write(*,*) 'odbf_setSurfaceType: query ---> ', trim(query)

    ! read the values from the query result
    call fSQL_prepare(db, trim(query), stmt, status=stat)
    call fSQL_get_many(stmt, nrows=numRows, ncols=numColumns, &
                       mode=FSQL_INT, status=stat)
    write(*,*) 'odbf_setSurfaceType: numRows = ', numRows, ', numColumns = ', numColumns
    if (numRows /= numRowsHeadTable) then
      write(*,*) 'odbf_setSurfaceType: numRows = ', numRows, &
                 ', numRowsHeadTable = ', numRowsHeadTable      
      call utl_abort('odbf_setSurfaceType: Number of rows found in mask query is &
                       not equal to total number of rows in head table')
    end if      
    allocate( columnValues(numRows, numColumns) )
    call fSQL_fill_matrix(stmt, columnValues)

    ! set the values of STYP and TTYP
    do headTableIndex = 1, numRows
      headIndex = headTableIndex + headIndexBegin - 1
      call obs_headSet_i(obsdat, OBS_STYP, headIndex, columnValues(headTableIndex,1))
      
      if (columnValues(headTableIndex,1) == 1 .and. columnValues(headTableIndex,2) == 1) then
        ! set terrain type to 0 over water with sea ice
        call obs_headSet_i(obsdat, OBS_TTYP, headIndex, 0)
      else
        ! otherwise set terrain type to -1
        call obs_headSet_i(obsdat, OBS_TTYP, headIndex, -1)
      end if

    end do

    ! close the obsDB file
    call fSQL_free_mem(stmt)
    call fSQL_finalize(stmt)
    call fSQL_close(db, stat) 

  end subroutine odbf_setSurfaceType

  !--------------------------------------------------------------------------
  ! odbf_copyToObsSpaceHeadChar
  !--------------------------------------------------------------------------
  subroutine odbf_copyToObsSpaceHeadChar(obsdat, headCharSqlNames, &
                                         headCharValues, headIndexBegin)
    !
    ! :Purpose: Copy character string values from a local table into
    !           obsSpaceData header rows. Currently, only the STATION ID
    !           and OBS_ITY (i.e. codeType from character sql column
    !           containing obs type name).
    !
    implicit none

    ! arguments:
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in)    :: headCharSqlNames(:)
    character(len=*), intent(in)    :: headCharValues(:,:)
    integer,          intent(in)    :: headIndexBegin

    ! locals:
    character(len=lenSqlName), allocatable :: stIdSqlName(:), codeTypeSqlName(:)
    integer :: columnIndex, headTableIndex, headIndex
    integer :: numRowsHeadTable, codeType

    numRowsHeadTable = size(headCharValues,1)

    ! Set the STATION ID
    stIdSqlName = odbf_sqlNameFromObsSpaceName('STID')
    columnIndex = utl_findloc(headCharSqlNames(:), trim(stIdSqlName(1)))
    if (columnIndex == 0) then
      call utl_abort('odbf_copyToObsSpaceHeadChar: Station ID column not found in sql table')
    end if
    do headTableIndex = 1, numRowsHeadTable
      headIndex = headTableIndex + headIndexBegin - 1
      if (headTableIndex == 1) then
        write(*,*) 'odbf_copyToObsSpaceHeadChar: set header char column   : ', &
                   trim(headCharSqlNames(columnIndex))
      end if
      call obs_set_c(obsdat, 'STID', &
                     headIndex, headCharValues(headTableIndex,columnIndex))
    end do

    ! Set the codeType (obs_ity) from a character column containing the obs type
    codeTypeSqlName = odbf_sqlNameFromObsSpaceName('ITY')
    columnIndex = utl_findloc(headCharSqlNames(:), trim(codeTypeSqlName(1)))
    if (columnIndex == 0) then
      call utl_abort('odbf_copyToObsSpaceHeadChar: Obs type column not found in sql table')
    end if
    do headTableIndex = 1, numRowsHeadTable
      headIndex = headTableIndex + headIndexBegin - 1
      if (headTableIndex == 1) then
        write(*,*) 'odbf_copyToObsSpaceHeadChar: set header char column   : ', &
                   trim(headCharSqlNames(columnIndex))
      end if
      codeType = codtyp_get_codtyp(trim(headCharValues(headTableIndex,columnIndex)))
      if (codeType == -1) then
        write(*,*) 'odbf_copyToObsSpaceHeadChar: obs type =', &
                   trim(headCharValues(headTableIndex,columnIndex))
        call utl_abort('odbf_copyToObsSpaceHeadChar: codtyp for this obs type not found') 
      end if
      call obs_headSet_i(obsdat, OBS_ITY, headIndex, codeType)
    end do ! do headTableIndex

  end subroutine odbf_copyToObsSpaceHeadChar

  !--------------------------------------------------------------------------
  ! odbf_copyToObsSpaceHeadDate
  !--------------------------------------------------------------------------
  subroutine odbf_copyToObsSpaceHeadDate(obsdat, headDateValues, headTimeValues, headIndexBegin)
    !
    ! :Purpose: Copy date values from a local table into
    !           obsSpaceData header rows.
    !
    implicit none

    ! arguments:
    type(struct_obs), intent(inout) :: obsdat
    integer,          intent(in)    :: headDateValues(:)
    integer,          intent(in)    :: headTimeValues(:)
    integer,          intent(in)    :: headIndexBegin

    ! locals:
    integer :: headTableIndex, headIndex
    integer :: numRowsHeadTable

    numRowsHeadTable = size(headDateValues,1)

    do headTableIndex = 1, numRowsHeadTable
      headIndex = headTableIndex + headIndexBegin - 1
      if (headTableIndex == 1) then
        write(*,*) 'odbf_copyToObsSpaceHeadDate: set header date/time column'
      end if
      call obs_headSet_i(obsdat, obs_dat, headIndex, headDateValues(headTableIndex))
      call obs_headSet_i(obsdat, obs_etm, headIndex, headTimeValues(headTableIndex))
    end do

  end subroutine odbf_copyToObsSpaceHeadDate

  !--------------------------------------------------------------------------
  ! odbf_copyToObsSpaceHead
  !--------------------------------------------------------------------------
  subroutine odbf_copyToObsSpaceHead(obsdat, headSqlNames, headPrimaryKey, &
                                     headValues, headIndexBegin)
    !
    ! :Purpose: Copy real and integer values from a local table into
    !           obsSpaceData header rows.
    !
    implicit none

    ! arguments:
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in)    :: headSqlNames(:)
    integer(8),       intent(in)    :: headPrimaryKey(:)
    real(8),          intent(in)    :: headValues(:,:)
    integer,          intent(in)    :: headIndexBegin

    ! locals:
    integer :: columnIndex, matchIndex, headTableIndex, headIndex
    integer :: numRowsHeadTable, obsColumnIndex

    numRowsHeadTable = size(headValues,1)

    write(*,*) 'odbf_copyToObsSpaceHead: set header primary key'
    do headTableIndex = 1, numRowsHeadTable
      headIndex = headTableIndex + headIndexBegin - 1
      call obs_setHeadPrimaryKey(obsdat, headIndex, headPrimaryKey(headTableIndex))
    end do

    column_loop: do columnIndex = 1, size(headSqlNames)
      matchIndex = utl_findloc(headMatchList(sqlColIndex,:), headSqlNames(columnIndex))

      if (matchIndex == 0) then

        write(*,*) 'odbf_copyToObsSpaceHead: unknown column name      : ', &
                   trim(headSqlNames(columnIndex))

      else

        obsColumnIndex = obs_columnIndexFromName(trim(headMatchList(obsColIndex,matchIndex)))

        headTable_loop: do headTableIndex = 1, numRowsHeadTable
          headIndex = headTableIndex + headIndexBegin - 1

          if (obs_columnDataType(obsColumnIndex) == 'real') then
            if ( obs_columnActive_RH(obsdat, obsColumnIndex) ) then
              if (headTableIndex == 1) then
                write(*,*) 'odbf_copyToObsSpaceHead: set header real column   : ', &
                           trim(headSqlNames(columnIndex))
              end if
              call obs_headSet_r(obsdat, obsColumnIndex, headIndex, &
                                 real(headValues(headTableIndex,columnIndex),&
                                 pre_obsReal))
            end if
          else if (obs_columnDataType(obsColumnIndex) == 'integer') then
            if ( obs_columnActive_IH(obsdat, obsColumnIndex) ) then
              if (headTableIndex == 1) then
                write(*,*) 'odbf_copyToObsSpaceHead: set header integer column: ', &
                           trim(headSqlNames(columnIndex))
              end if
              call obs_headSet_i(obsdat, obsColumnIndex, &
                                 headIndex, nint(headValues(headTableIndex,columnIndex)))
            end if
          else
            call utl_abort('odbf_copyToObsSpaceHead: unknown data type for obs header column')
          end if ! if (obs_columnDataType(obsColumnIndex) == 'real')

        end do headTable_loop

      end if ! if (matchIndex == 0)

    end do column_loop

  end subroutine odbf_copyToObsSpaceHead

  !--------------------------------------------------------------------------
  ! odbf_copyToObsSpaceBody
  !--------------------------------------------------------------------------
  subroutine odbf_copyToObsSpaceBody(obsdat, bodySqlNames, &
                                     bodyPrimaryKey, bodyHeadKey, bodyValues, &
                                     bodyIndexBegin, headIndexBegin)
    !
    ! :Purpose: Copy real and integer values from a local table into
    !           obsSpaceData body rows. Note: this version currently
    !           assumes that only 1 observed quantity is present for
    !           each row of the sqlite table. This is likely only valid
    !           for radiance observation types and therefore modifications
    !           will be required for other observation types.
    !
    implicit none

    ! arguments:
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in)    :: bodySqlNames(:)
    integer(8),       intent(in)    :: bodyPrimaryKey(:)
    integer(8),       intent(in)    :: bodyHeadKey(:)
    real(8),          intent(in)    :: bodyValues(:,:)
    integer,          intent(in)    :: bodyIndexBegin
    integer,          intent(in)    :: headIndexBegin

    ! locals:
    character(len=lenSqlName), allocatable :: obsValueSqlNames(:)
    integer :: columnIndex, matchIndex, bodyTableIndex, bodyIndex, headIndex
    integer :: numRowsBodyTable, obsNlv, obsValueIndex, numObsValues
    integer(8) :: lastHeadKey
    integer, allocatable :: bodyColumnIndexObsValueList(:)
    integer, allocatable :: obsVarNoList(:)
    integer, allocatable :: matchIndexVec(:)
    integer, allocatable :: obsColumnIndex(:)
    logical              :: firstHead

    numRowsBodyTable = size(bodyValues,1)

    ! initialize some arrays to save lots of time in the main loops
    allocate(matchIndexVec(size(bodySqlNames)))
    do columnIndex = 1, size(bodySqlNames)
      matchIndexVec(columnIndex) = utl_findloc(bodyMatchList(sqlColIndex,:), &
                                               bodySqlNames(columnIndex))
    end do
    allocate(obsColumnIndex(numBodyMatch))
    do matchIndex = 1, numBodyMatch
      obsColumnIndex(matchIndex) = obs_columnIndexFromName(bodyMatchList(obsColIndex,matchIndex))
    end do

    ! figure out column indexes for observation values (OBS_VAR)
    obsValueSqlNames = odbf_sqlNameFromObsSpaceName('VAR')
    numObsValues = size(obsValueSqlNames)
    allocate(obsVarNoList(numObsValues))
    allocate(bodyColumnIndexObsValueList(numObsValues))
    do obsValueIndex = 1, numObsValues
      bodyColumnIndexObsValueList(obsValueIndex) = &
           utl_findloc(bodySqlNames(:), obsValueSqlNames(obsValueIndex))
      if (bodyColumnIndexObsValueList(obsValueIndex) == 0) then
        write(*,*) 'odbf_copyToObsSpaceBody: obsValueSqlName = ', &
                   trim(obsValueSqlNames(obsValueIndex))
        call utl_abort('odbf_copyToObsSpaceBody: column with obs value not present')
      end if
      ! determine varNo for the observation value
      obsVarNoList(obsValueIndex) = odbf_varNoFromSqlName(obsValueSqlNames(obsValueIndex))
      write(*,*) 'odbf_copyToObsSpaceBody: obsVarNo = ', obsVarNoList(obsValueIndex)
    end do

    lastHeadKey = 0
    firstHead = .true.
    headIndex = headIndexBegin - 1
    bodyIndex = bodyIndexBegin - 1

    bodyIndex_loop: do bodyTableIndex = 1, numRowsBodyTable

      obsValueIndex_loop: do obsValueIndex = 1, numObsValues

        ! initialize count of number of body rows for each header row (OBS_NLV)
        if ( firstHead .or. (bodyHeadKey(bodyTableIndex) /= lastHeadKey) ) then
          headIndex = headIndex + 1
          call obs_headSet_i(obsdat, OBS_NLV, headIndex, 0)
          lastHeadKey = bodyHeadKey(bodyTableIndex)
          firstHead = .false.
        end if

        ! check that the primary key for header table matches the value in the body table
        if ( obs_headPrimaryKey( obsdat, headIndex ) /= &
             bodyHeadKey(bodyTableIndex) ) then
          write(*,*) 'odbf_copyToObsSpaceBody: primary key in HEADER table = ', &
                     obs_headPrimaryKey( obsdat, headIndex )
          write(*,*) 'odbf_copyToObsSpaceBody: same key in BODY table      = ', &
                     bodyHeadKey(bodyTableIndex)
          call utl_abort('odbf_copyToObsSpaceBody: Primary key of HEADER table ' // &
                         'not equal to value in BODY table')
        end if

        ! check if obs value is null/missing
        if ( bodyValues(bodyTableIndex,bodyColumnIndexObsValueList(obsValueIndex)) == &
             MPC_missingValue_R8 ) then
          cycle obsValueIndex_loop
        end if

        ! check if element id is in list
        if ( utl_findloc(elemIdList(1:numElemIdList),obsVarNoList(obsValueIndex)) == 0 ) then
          cycle obsValueIndex_loop
        end if

        ! add to count of number of body rows for each header row (OBS_NLV)
        obsNLV = obs_headElem_i(obsdat, OBS_NLV, headIndex)
        call obs_headSet_i(obsdat, OBS_NLV, headIndex, obsNLV + 1)

        bodyIndex = bodyIndex + 1

        ! copy body primary key to obsSpaceData
        if (bodyTableIndex == 1) then
          write(*,*) 'odbf_copyToObsSpaceBody: set body primary key'
        end if
        call obs_setBodyPrimaryKey(obsdat, bodyIndex, bodyPrimaryKey(bodyTableIndex))

        ! set the varNo for this obsValue
        call obs_bodySet_i(obsdat, OBS_VNM, bodyIndex, obsVarNoList(obsValueIndex))

        ! copy real and integer values into obsSpaceData
        columnIndex_loop: do columnIndex = 1, size(bodySqlNames)
          matchIndex = matchIndexVec(columnIndex)

          if (matchIndex == 0) then

            if (bodyTableIndex == 1) then
              write(*,*) 'odbf_copyToObsSpaceBody: unknown column name    : ', &
                         trim(bodySqlNames(columnIndex))
            end if

          else

            ! if this column corresponds to the obs value, then check if it is the one we want
            if ( obsColumnIndex(matchIndex) == OBS_VAR ) then
              if ( columnIndex /= bodyColumnIndexObsValueList(obsValueIndex) ) then
                ! skip this column
                cycle columnIndex_loop
                if (bodyTableIndex == 1) then
                  write(*,*) 'odbf_copyToObsSpaceBody: skip obs body column   : ', &
                             trim(bodySqlNames(columnIndex)), &
                             ' for obsValueIndex = ', obsValueIndex
                end if
              end if
            end if

            if (obs_columnDataType(obsColumnIndex(matchIndex)) == 'real') then
              ! real values
              if ( obs_columnActive_RB(obsdat, obsColumnIndex(matchIndex)) ) then
                if (bodyTableIndex == 1) then
                  write(*,*) 'odbf_copyToObsSpaceBody: set body real column   : ', &
                             trim(bodySqlNames(columnIndex))
                end if
                call obs_bodySet_r(obsdat, obsColumnIndex(matchIndex), bodyIndex, &
                                   bodyValues(bodyTableIndex,columnIndex))
              end if
            else if (obs_columnDataType(obsColumnIndex(matchIndex)) == 'integer') then
              ! integer values
              if ( obs_columnActive_IB(obsdat, obsColumnIndex(matchIndex)) ) then
                if (bodyTableIndex == 1) then
                  write(*,*) 'odbf_copyToObsSpaceBody: set body integer column: ', &
                             trim(bodySqlNames(columnIndex))
                end if
                call obs_bodySet_i(obsdat, obsColumnIndex(matchIndex), &
                                   bodyIndex, nint(bodyValues(bodyTableIndex,columnIndex)))
              end if
            else
              call utl_abort('odbf_copyToObsSpaceBody: unknown data type for obs body column')
            end if ! if (obs_columnDataType(obsColumnIndex(matchIndex))

          end if ! if (matchIndex == 0)

        end do columnIndex_loop

      end do obsValueIndex_loop

    end do bodyIndex_loop

    deallocate(matchIndexVec)
    deallocate(obsColumnIndex)
    deallocate(obsVarNoList)
    deallocate(bodyColumnIndexObsValueList)

  end subroutine odbf_copyToObsSpaceBody

  !--------------------------------------------------------------------------
  ! odbf_adjustValues
  !--------------------------------------------------------------------------
  subroutine odbf_adjustValues(obsdat, headIndexBegin, headIndexEnd)
    !
    ! :Purpose: Adjust units and other minor modifications of some
    !           obsSpaceData columns after transfer from sqlite files
    !
    implicit none

    ! arguments:
    type(struct_obs), intent(inout) :: obsdat
    integer,          intent(in)    :: headIndexBegin
    integer,          intent(in)    :: headIndexEnd

    ! locals:
    integer :: headIndex, bodyIndexStart, bodyIndexEnd, bodyIndex
    integer :: instrument, obsSat, codeType, sensor
    real(8) :: obsLon, obsLat, surfEmiss
    character(len=2) :: obsFamily

    do headIndex = headIndexBegin, headIndexEnd

      obsFamily = obs_getfamily( obsdat, headIndex )
      bodyIndexStart = obs_headElem_i(obsdat, OBS_RLN, headIndex)
      bodyIndexEnd   = obs_headElem_i(obsdat, OBS_NLV, headIndex) + &
                       bodyIndexStart - 1

      ! Convert lon-lat from degrees to radians

      obsLon = obs_headElem_r( obsdat, OBS_LON, headIndex )
      obsLat = obs_headElem_r( obsdat, OBS_LAT, headIndex )

      if ( obsLon < 0.0D0 ) obsLon = obsLon + 360.0D0
      obsLon = obsLon * MPC_RADIANS_PER_DEGREE_R8
      obsLat = obsLat * MPC_RADIANS_PER_DEGREE_R8

      call obs_headSet_r(obsdat, OBS_LON, headIndex, real(obsLon,pre_obsReal))
      call obs_headSet_r(obsdat, OBS_LAT, headIndex, real(obsLat,pre_obsReal))

      ! Set global and observation flags to zero if specified

      if (setObsFlagZero) then
        call obs_headSet_i(obsdat, OBS_ST1, headIndex, 0)
        do bodyIndex = bodyIndexStart, bodyIndexEnd
          call obs_bodySet_i(obsdat, OBS_FLG, bodyIndex, 0)
        end do
      end if

      ! Various adjustment for radiance observations

      if ( obsFamily == 'TO' ) then

        instrument = obs_headElem_i( obsdat, OBS_INS, headIndex )
        obsSat     = obs_headElem_i( obsdat, OBS_SAT, headIndex )
        codeType   = obs_headElem_i( obsdat, OBS_ITY, headIndex )
        sensor     = obs_headElem_i( obsdat, OBS_SEN, headIndex )

        ! set sensor to missing if not amsua/b, mhs or atms

        if ( codeType /= codtyp_get_codtyp('amsua') .and. &
             codeType /= codtyp_get_codtyp('amsub') .and. &
             codeType /= codtyp_get_codtyp('mhs') .and. &
             codeType /= codtyp_get_codtyp('atms') ) then
          sensor = nint(MPC_missingValue_R8)
        end if

        ! modify OBS_SAT, OBS_INS and OBS_SEN

        if ( instrument == 420 ) obsSat  = 784
        if ( codeType == 202 .and. instrument == 620 ) instrument = 2046
        if ( sensor == nint(MPC_missingValue_R8) ) then
          sensor = 0
          if (instrument == nint(MPC_missingValue_R8) ) instrument = 0
        else
          instrument = obsu_cvt_obs_instrum(sensor)
        end if
        call obs_headSet_i(obsdat, OBS_INS, headIndex, instrument)
        call obs_headSet_i(obsdat, OBS_SAT, headIndex, obsSat)
        call obs_headSet_i(obsdat, OBS_SEN, headIndex, sensor)

        ! change units for surface emissivity

        if (obs_columnActive_RB(obsdat, OBS_SEM)) then
          do bodyIndex = bodyIndexStart, bodyIndexEnd
            surfEmiss = obs_bodyElem_r(obsdat, OBS_SEM, bodyIndex)
            surfEmiss = surfEmiss * 0.01D0
            call obs_bodySet_r(obsdat, OBS_SEM, bodyIndex, surfEmiss)
          end do
        end if

      end if ! obsFamily = 'TO'

    end do ! do headIndex

  end subroutine odbf_adjustValues

  !--------------------------------------------------------------------------
  ! odbf_sqlNameFromObsSpaceName
  !--------------------------------------------------------------------------
  function odbf_sqlNameFromObsSpaceName(obsSpaceName) result(sqlName)
    !
    ! :Purpose: Return the corresponding sql file column name for a
    !           given obsSpaceData column name from the matching
    !           tables.
    !
    implicit none

    ! arguments:
    character(len=*), intent(in)           :: obsSpaceName
    character(len=lenSqlName), allocatable :: sqlName(:)

    ! locals:
    integer                   :: numMatchFound, matchFoundIndex
    integer, allocatable      :: matchIndexList(:)

    if (allocated(sqlName)) deallocate(sqlName)

    ! first try the body matching list
    matchIndexList = utl_findlocs(bodyMatchList(obsColIndex,:), trim(obsSpaceName))
    if (matchIndexList(1) > 0) then
      numMatchFound = size(matchIndexList)
      allocate(sqlName(numMatchFound))
      do matchFoundIndex = 1, numMatchFound
        sqlName(matchFoundIndex) = bodyMatchList(sqlColIndex,matchIndexList(matchFoundIndex))
      end do
      return
    end if

    ! now try the header matching list
    matchIndexList = utl_findlocs(headMatchList(obsColIndex,:), trim(obsSpaceName))
    if (matchIndexList(1) > 0) then
      numMatchFound = size(matchIndexList)
      allocate(sqlName(numMatchFound))
      do matchFoundIndex = 1, numMatchFound
        sqlName(matchFoundIndex) = headMatchList(sqlColIndex,matchIndexList(matchFoundIndex))
      end do
      return
    end if

    ! not found in either list, abort
    write(*,*) 'odbf_sqlNameFromObsSpaceName: requested obsSpace name = ', trim(obsSpaceName)
    call utl_abort('odbf_sqlNameFromObsSpaceName: obsSpace name not found in matching list')
    
  end function odbf_sqlNameFromObsSpaceName

  !--------------------------------------------------------------------------
  ! odbf_midasTabColFromObsSpaceName
  !--------------------------------------------------------------------------
  function odbf_midasTabColFromObsSpaceName(obsSpaceName, midasSQLColumnList) result(sqlColName)
    !
    ! :Purpose: Return the corresponding sql file column name for a
    !           given obsSpaceData column name from the midas table
    !           matching list.
    !
    implicit none

    ! arguments:
    character(len=*), intent(in) :: obsSpaceName
    character(len=lenSqlName), intent(in)    :: midasSQLColumnList(:,:)
    character(len=lenSqlName)    :: sqlColName
    ! locals:
    integer :: matchIndex

    ! look in the midas table matching list
    matchIndex = utl_findloc(midasSQLColumnList(obsColIndex,:), trim(obsSpaceName))
    if (matchIndex > 0) then
      sqlColName = midasSQLColumnList(sqlColIndex,matchIndex)
      return
    end if

    ! not found, abort
    write(*,*) 'odbf_midasTabColFromObsSpaceName: requested obsSpace name = ', trim(obsSpaceName)
    call utl_abort('odbf_midasTabColFromObsSpaceName: obsSpace name not found in midasSQLNamesList')
    
  end function odbf_midasTabColFromObsSpaceName

  !--------------------------------------------------------------------------
  ! odbf_varNoFromSqlName
  !--------------------------------------------------------------------------
  function odbf_varNoFromSqlName(sqlName) result(varNo)
    !
    ! :Purpose: Return the bufr element id number from the corresponding
    !           sql file column name of an observed value.
    !
    implicit none

    ! arguments:
    character(len=*), intent(in) :: sqlName
    integer                      :: varNo

    ! locals:
    integer           :: matchIndex
    character(len=10) :: varNoStr
  
    matchIndex = utl_findloc(varNoList(sqlColIndex,:), trim(sqlName))
    if (matchIndex > 0) then
      varNoStr = varNoList(varNoColIndex,matchIndex)
      read(varNoStr,*) varNo
      return
    end if

    write(*,*) 'odbf_varNoFromSqlName: requested sqlName = ', trim(sqlName)
    call utl_abort('odbf_varNoFromSqlName: not found in varNo list')
    
  end function odbf_varNoFromSqlName

  !--------------------------------------------------------------------------
  ! odbf_updateFile
  !--------------------------------------------------------------------------
  subroutine odbf_updateFile(obsdat, fileName, familyType, fileIndex)
    !
    ! :Purpose: Call subroutines to update MIDAS Header and Body tables
    !
    implicit none
    
    ! arguments
    type (struct_obs), intent(inout) :: obsdat
    character(len=*),  intent(in)    :: fileName, familyType
    integer,           intent(in)    :: fileIndex

    ! Locals:
    logical                          :: fileExists

    call utl_tmg_start(14,'----UpdateObsDBfile')
    
    ! prepare to create obsDB from scratch if the file does not exist
    inquire(file=trim(fileName), exist=fileExists)
    if (.not. fileExists) call odbf_setup()

    ! Check if the Midas Header Table needs to be updated, specified from namelist
    if ( .not. utl_isNamelistPresent('namObsDbMIDASHeaderUpdate','./flnml') ) then
      if ( mmpi_myid == 0 ) then
        write(*,*) 'odbf_updateFile: namObsDbMIDASHeaderUpdate is missing in the namelist.'
        write(*,*) '                 The MIDAS Header Output Table will not be updated.'
      end if
    else
      ! Update the MIDAS Header Output Table
      call odbf_insertInMidasHeaderTable(obsdat, fileIndex, fileName, familyType)
    end if

    ! Check if the Midas Body Table needs to be updated, specified from namelist
    if ( .not. utl_isNamelistPresent('namObsDbMIDASBodyUpdate','./flnml') ) then
      if ( mmpi_myid == 0 ) then
        write(*,*) 'odbf_updateFile: namObsDbMIDASBodyUpdate is missing in the namelist.'
        write(*,*) '                 The MIDAS Body Output Table will not be updated.'
      end if
    else
      ! Update the MIDAS Body Output Table
      call odbf_insertInMidasBodyTable(obsdat, fileIndex, fileName, familyType)
    end if

    call utl_tmg_stop(14)

  end subroutine odbf_updateFile

  !--------------------------------------------------------------------------
  ! odbf_insertInMidasHeaderTable
  !--------------------------------------------------------------------------
  subroutine odbf_insertInMidasHeaderTable(obsdat, fileIndex, fileName, familyType)
    !
    ! :Purpose: Insert selected columns in the MIDAS Header Output table using
    !           values from obsSpaceData. If the MIDAS Header table does not already
    !           exist, it is created by copying the observation table.
    !           A single table is created that contains all quantities being
    !           updated. Unlike the observation table, each observed variable
    !           is stored in a separate row and all quantities are in columns
    !           (e.g. ETOP, VTOP, ECF,...).
    !
    implicit none
    
    ! arguments:
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in)    :: fileName
    character(len=*), intent(in)    :: familyType 
    integer,          intent(in)    :: fileIndex
    
    ! locals:
    type(fSQL_STATUS)    :: stat ! sqlite error status
    type(fSQL_DATABASE)  :: db   ! sqlite file handle
    type(fSQL_STATEMENT) :: stmt ! precompiled sqlite statements
    integer(8)           :: obsIdo
    integer              :: columnParamIndex
    integer              :: midasKey, obsIdf, updateItemIndex, updateValue_i
    integer              :: headIndex, maxNumHeader
    integer              :: obsSpaceColIndexSource, fnom, fclos, nulnam, ierr
    integer              :: maxNumHeaderAllMpi(mmpi_nprocs), obsSpaceColIndexSourceArr(15)
    real(8)              :: updateValue_r
    character(len=3000)  :: query, queryCreateTable, queryInsertInTable, queryForValues
    character(len=5000)  :: tableInsertColumnList
    character(len=20)    :: sqlDataType
    character(len=lenSqlName) :: sqlColumnName
    character(len=4)     :: obsSpaceColumnName
    logical              :: midasTableExists
    logical, save        :: nmlAlreadyRead = .false.
    character(len=6), parameter     :: midasTableType = 'header' !Define the type of MIDAS table: header/body
    integer, parameter   :: maxItemNumber = 15
    ! namelist variables
    integer,          save :: numberUpdateItems             ! MUST NOT BE INCLUDED IN NAMELIST!
    character(len=4), save :: updateItemList(maxItemNumber) ! obsSpace column names used to update the file

    namelist/namObsDbMIDASHeaderUpdate/ numberUpdateItems, updateItemList

    write(*,*)
    write(*,*) 'odbf_insertInMidasHeaderTable: Starting'
    write(*,*)
    write(*,*) 'odbf_insertInMidasHeaderTable: FileName   : ', trim(FileName)
    write(*,*) 'odbf_insertInMidasHeaderTable: FamilyType : ', FamilyType
    write(*,*) 'odbf_insertInMidasHeaderTable: fileIndex   : ', fileIndex

    if (.not. nmlAlreadyRead) then
      nmlAlreadyRead = .true.

      ! set default values of namelist variables
      updateItemList(:) = ''
      numberUpdateItems = MPC_missingValue_INT

      ! Read the namelist for directives
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam, nml=namObsDbMIDASHeaderUpdate, iostat=ierr)
      if ( ierr /= 0 ) call utl_abort('odbf_insertInMidasHeaderTable: Error reading namelist')
      ierr = fclos(nulnam)
      if ( numberUpdateItems /= MPC_missingValue_INT) then
        call utl_abort('odbf_insertInMidasHeaderTable: check namObsDbMIDASHeaderUpdate namelist section: numberUpdateItems should be removed')
      end if
      numberUpdateItems = 0
      do updateItemIndex = 1, maxItemNumber
        if (trim(updateItemList(updateItemIndex)) == '') exit
        numberUpdateItems = numberUpdateItems + 1
      end do
      if ( mmpi_myid == 0 ) then
        write(*, nml=namObsDbMIDASHeaderUpdate)
      end if
    end if ! not nmlAlreadyRead

    if (numberUpdateItems == 0) then
      write(*,*) 'odbf_insertInMidasHeaderTable: numberUpdateItems=0. ' // &
                 'MIDAS Header Output Table will not be updated.'
      return
    end if

    ! check if midasTable already exists in the file
    midasTableExists = sqlu_sqlTableExists(fileName, midasHeadTableName)

    if (.not. midasTableExists) then
      ! create midasTable by copying rearranging contents of observation table
      call odbf_createMidasHeaderTable(fileName)
    
      ! open the obsDB file
      call fSQL_open(db, trim(fileName), stat)
      if ( fSQL_error(stat) /= FSQL_OK ) then
        write(*,*) 'odbf_insertInMidasHeaderTable: fSQL_open: ', fSQL_errmsg(stat)
        call utl_abort('odbf_insertInMidasHeaderTable: fSQL_open')
      end if

      ! Obtain the max number of header rows per mpi task
      maxNumHeader = 0
      HEADER1: do headIndex = 1, obs_numHeader(obsdat)
        obsIdf = obs_headElem_i( obsdat, OBS_IDF, headIndex )
        if ( obsIdf /= fileIndex ) cycle HEADER1
        maxNumHeader = maxNumHeader + 1
      end do HEADER1

      call rpn_comm_allGather(maxNumHeader,       1, 'mpi_integer',  &
                              maxNumHeaderAllMpi, 1, 'mpi_integer', &
                              'GRID', ierr)

      ! Set the midasKey to start counting based on the latest value from the previous 
      ! mpi task
      if( mmpi_myid == 0 ) then
        midasKey = 0
      else
        midasKey = sum(maxNumHeaderAllMpi(1:mmpi_myid))
      end if

      ! set the primary key, keys to main obsDB tables and other basic info
      query = 'insert into ' // trim(midasHeadTableName) // '(' // &
              trim(midasHeadKeySqlName) // ',' //trim(obsHeadKeySqlName) // &
              ') values(?, ?);'

      write(*,*) 'odbf_insertInMidasHeaderTable: query = ', trim(query)
      call fSQL_prepare(db, query, stmt, stat)
      call fSQL_begin(db)

      HEADER: do headIndex = 1, obs_numHeader(obsdat)

        obsIdf = obs_headElem_i( obsdat, OBS_IDF, headIndex )
        if ( obsIdf /= fileIndex ) cycle HEADER

        midasKey = midasKey +1
        call fSQL_bind_param(stmt, PARAM_INDEX=1, INT_VAR=midasKey)
        obsIdo = obs_headPrimaryKey(obsdat, headIndex)
        call fSQL_bind_param(stmt, PARAM_INDEX = 2, INT8_VAR  = obsIdo)

        call fSQL_exec_stmt (stmt)
      
      end do HEADER

      call fSQL_finalize(stmt)
      call fSQL_commit(db)

      ! close the obsDB file
      call fSQL_close(db, stat) 
    else
      write(*,*) 'odbf_insertInMidasHeaderTable: the midas header output table already exists, ' // &
                 'proceed to join with temporary table that has new columns.'
    end if ! .not.midasTableExists

    ! open the obsDB file
    call fSQL_open(db, trim(fileName), stat)
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_insertInMidasHeaderTable: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort('odbf_insertInMidasHeaderTable: fSQL_open')
    end if
    
    tableInsertColumnList = ''
    obsSpaceColIndexSourceArr(:) = mpc_missingValue_int
    call getCreateTableInsertQueries(numberUpdateItems, updateItemList, midasTableType, &
                                     queryCreateTable, queryInsertInTable, &
                                     tableInsertColumnList, obsSpaceColIndexSourceArr(:))

    ! Create a temporary table with new columns and values from obsSpaceData
    write(*,*) 'odbf_insertInMidasHeaderTable: queryCreateTable   -->', trim(queryCreateTable)
    call fSQL_do_many(db, queryCreateTable, stat)
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'fSQL_do_many: ', fSQL_errmsg(stat)
      call utl_abort('odbf_insertInMidasHeaderTable: Problem with fSQL_do_many')
    end if 

    ! Prepare to insert into the table
    write(*,*) 'odbf_insertInMidasHeaderTable: queryInsertInTable -->', trim(queryInsertInTable)
    call fSQL_prepare(db, queryInsertInTable, stmt, stat)
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_insertInMidasHeaderTable: fSQL_prepare: ', fSQL_errmsg(stat)
      call utl_abort('odbf_insertInMidasHeaderTable: fSQL_prepare')
    end if

    call fSQL_begin(db)
    HEADER2: do headIndex = 1, obs_numHeader(obsdat)
      obsIdf = obs_headElem_i( obsdat, OBS_IDF, headIndex )
      if ( obsIdf /= fileIndex ) cycle HEADER2

      obsIdo  = obs_headPrimaryKey( obsdat, headIndex )
      call fSQL_bind_param(stmt, PARAM_INDEX=1, INT8_VAR=obsIdo)

      do updateItemIndex = 1, numberUpdateItems

        obsSpaceColIndexSource = obsSpaceColIndexSourceArr(updateItemIndex)
        columnParamIndex = updateItemIndex + 1

        ! update the value, but set to null if it is missing
        if (obs_columnDataType(obsSpaceColIndexSource) == 'real') then
          updateValue_r = obs_headElem_r(obsdat, obsSpaceColIndexSource, headIndex)

          if ( updateValue_r == obs_missingValue_R ) then
            call fSQL_bind_param(stmt, PARAM_INDEX=columnParamIndex)  ! sql null values
          else
            call fSQL_bind_param(stmt, PARAM_INDEX=columnParamIndex, REAL8_VAR=updateValue_r)
          end if
        else
          updateValue_i = obs_headElem_i(obsdat, obsSpaceColIndexSource, headIndex)
          if ( updateValue_i == mpc_missingValue_int ) then
            call fSQL_bind_param(stmt, PARAM_INDEX=columnParamIndex)  ! sql null values
          else
            call fSQL_bind_param(stmt, PARAM_INDEX=columnParamIndex, INT_VAR=updateValue_i)
          end if
        end if ! if (obs_columnDataType(obsSpaceColIndexSource) == 'real')
      end do ! do updateItemIndex

      call fSQL_exec_stmt(stmt)

    end do HEADER2

    call fSQL_finalize(stmt)
    call fSQL_commit(db)

    ! close the obsDB file
    call fSQL_close(db, stat)

    call mergeTableInMidasTables(fileName, midasHeadTableName, obsHeadKeySqlName, &
                                 tableInsertColumnList)

    write(*,*)
    write(*,*) 'odbf_insertInMidasHeaderTable: finished'
    write(*,*)

  end subroutine odbf_insertInMidasHeaderTable

  !--------------------------------------------------------------------------
  ! odbf_insertInMidasBodyTable
  !--------------------------------------------------------------------------
  subroutine odbf_insertInMidasBodyTable(obsdat, fileIndex, fileName, familyType)
    !
    ! :Purpose: Insert selected columns in the MIDAS body table using
    !           values from obsSpaceData. If the MIDAS Body table does not already 
    !           exist, it is created by copying the observation table.
    !           A single table is created that contains all quantities being
    !           updated. Unlike the observation table, each observed variable
    !           is stored in a separate row and all quantities are in columns
    !           (e.g. obsValue, OMP, OMA,...).
    !
    implicit none

    ! arguments:
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in)    :: fileName   
    character(len=*), intent(in)    :: familyType
    integer,          intent(in)    :: fileIndex

    ! locals:
    type(fSQL_STATUS)    :: stat ! sqlite error status
    type(fSQL_DATABASE)  :: db   ! sqlite file handle
    type(fSQL_STATEMENT) :: stmt ! precompiled sqlite statements
    integer(8)           :: obsIdo, obsIdd
    integer              :: columnParamIndex, columnIndex
    integer              :: obsIdf, obsVarNo, midasKey, updateItemIndex, updateValue_i
    integer              :: headIndex, bodyIndex, bodyIndexBegin, bodyIndexEnd, maxNumBody
    integer              :: obsSpaceColIndexSource, fnom, fclos, nulnam, ierr
    integer              :: maxNumBodyAllMpi(mmpi_nprocs), obsSpaceColIndexSourceArr(15)
    real(8)              :: updateValue_r, obsValue, obsPPP, obsVAR
    character(len=4)     :: obsSpaceColumnName
    character(len=3000)  :: query, queryCreateTable, queryInsertInTable, queryForValues
    character(len=5000)  :: tableInsertColumnList
    character(len=20)        :: sqlDataType
    logical              :: midasTableExists
    logical, save        :: nmlAlreadyRead = .false.
    character(len=6), parameter  :: midasTableType='body' ! Define the type of MIDAS table: header/body
    integer, parameter   :: maxItemNumber = 15
    ! namelist variables
    integer,          save :: numberUpdateItems             ! MUST NOT BE INCLUDED IN NAMELIST!
    character(len=4), save :: updateItemList(maxItemNumber) ! obsSpace column names used to update the file

    namelist/namObsDbMIDASBodyUpdate/ numberUpdateItems, updateItemList

    write(*,*)
    write(*,*) 'odbf_insertInMidasBodyTable: Starting'
    write(*,*)
    write(*,*) 'odbf_insertInMidasBodyTable: FileName   : ', trim(FileName)
    write(*,*) 'odbf_insertInMidasBodyTable: FamilyType : ', FamilyType
    write(*,*) 'odbf_insertInMidasBodyTable: fileIndex : ', fileIndex

    if (.not. nmlAlreadyRead) then
      nmlAlreadyRead = .true.

      ! set default values of namelist variables
      updateItemList(:) = ''
      numberUpdateItems = MPC_missingValue_INT

      ! Read the namelist for directives
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam, nml=namObsDbMIDASBodyUpdate, iostat=ierr)
      if ( ierr /= 0 ) call utl_abort('odbf_insertInMidasBodyTable: Error reading namelist')
      ierr = fclos(nulnam)
      if (numberUpdateItems /=  MPC_missingValue_INT) then
        call utl_abort('odbf_insertInMidasBodyTable: check namObsDbMIDASBodyUpdate namelist section: numberUpdateItems should be removed')
      end if
      numberUpdateItems = 0
      do updateItemIndex = 1, maxItemNumber
        if (trim(updateItemList(updateItemIndex)) == '') exit
        numberUpdateItems = numberUpdateItems + 1
      end do
      ! Add VNM/PPP/VAR/FLG to the updateItemList to ensure they are included in the sqlite table with values from obsSpaceData
      do columnIndex = 1, numBodyMidasTableRequired
        numberUpdateItems = numberUpdateItems + 1
        updateItemList(numberUpdateitems) = trim(midasBodyNamesList(2,columnIndex))
      end do
         
      numberUpdateItems = numberUpdateItems + 1
      updateItemList(numberUpdateItems) = 'FLG'      

      if ( mmpi_myid == 0 ) then
        write(*,*) 'odbf_insertInMidasBodyTable: NOTE: the FLG/VNM/PPP/VAR columns are always added to update list'
        write(*, nml=namObsDbMIDASBodyUpdate)
      end if
    end if ! not nmlAlreadyRead

    if (numberUpdateItems == 0) then
      write(*,*) 'odbf_insertInMidasBodyTable: numberUpdateItems=0. ' // &
                 'MIDAS Body Output Table will not be updated.'
      return
    end if

    ! check if midasTable already exists in the file
    midasTableExists = sqlu_sqlTableExists(fileName, midasBodyTableName)

    if (.not. midasTableExists) then
      ! create midasTable by copying rearranging contents of observation table
      call odbf_createMidasBodyTable(fileName)

      ! open the obsDB file
      call fSQL_open(db, trim(fileName), stat)
      if ( fSQL_error(stat) /= FSQL_OK ) then
        write(*,*) 'odbf_insertInMidasBodyTable: fSQL_open: ', fSQL_errmsg(stat)
        call utl_abort('odbf_insertInMidasBodyTable: fSQL_open')
      end if

      ! Obtain the max number of body rows per mpi task
      maxNumBody = 0
      HEADER1: do headIndex = 1, obs_numHeader(obsdat)
        obsIdf = obs_headElem_i( obsdat, OBS_IDF, headIndex )
        if ( obsIdf /= fileIndex ) cycle HEADER1

        bodyIndexBegin = obs_headElem_i( obsdat, OBS_RLN, headIndex )
        bodyIndexEnd = bodyIndexBegin + &
                       obs_headElem_i( obsdat, OBS_NLV, headIndex ) - 1
        
        BODY1: do bodyIndex = bodyIndexBegin, bodyIndexEnd
          obsValue = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex)
          if ( obsValue == obs_missingValue_R ) cycle BODY1
          maxNumBody = maxNumBody + 1
        end do BODY1
      end do HEADER1
      
      call rpn_comm_allGather(maxNumBody,       1, 'mpi_integer',  &
                              maxNumBodyAllMpi, 1, 'mpi_integer', &
                            'GRID', ierr )

      ! Set the midasKey to start counting based on the latest value from the previous 
      ! mpi task
      if( mmpi_myid == 0 ) then
        midasKey = 0
      else
        midasKey = sum(maxNumBodyAllMpi(1:mmpi_myid))
      end if
                      
      ! set the primary key, keys to main obsDB tables and other basic info
      query = 'insert into ' // trim(midasBodyTableName) // '(' // &
              trim(midasBodyKeySqlName) // ',' // trim(obsHeadKeySqlName) // ',' // &
              trim(obsBodyKeySqlName)  // ') values(?,?,?);'
      write(*,*) 'odbf_insertInMidasBodyTable: query = ', trim(query)
      call fSQL_prepare(db, query, stmt, stat)
      call fSQL_begin(db)
      HEADER: do headIndex = 1, obs_numHeader(obsdat)

        obsIdf = obs_headElem_i( obsdat, OBS_IDF, headIndex )
        if ( obsIdf /= fileIndex ) cycle HEADER

        bodyIndexBegin = obs_headElem_i( obsdat, OBS_RLN, headIndex )
        bodyIndexEnd = bodyIndexBegin + &
                       obs_headElem_i( obsdat, OBS_NLV, headIndex ) - 1

        BODY: do bodyIndex = bodyIndexBegin, bodyIndexEnd

          ! do not try to update if the observed value is missing
          obsValue = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex)
          if ( obsValue == obs_missingValue_R ) cycle BODY

          midasKey = midasKey + 1
          call fSQL_bind_param(stmt, PARAM_INDEX=1, INT_VAR=midasKey)
          obsIdo  = obs_headPrimaryKey( obsdat, headIndex )
          call fSQL_bind_param(stmt, PARAM_INDEX=2, INT8_VAR=obsIdo)
          obsIdd  = obs_bodyPrimaryKey(obsdat, bodyIndex)
          call fSQL_bind_param(stmt, PARAM_INDEX=3, INT8_VAR=obsIdd)

          call fSQL_exec_stmt(stmt)

        end do BODY

      end do HEADER

      call fSQL_finalize(stmt)
      call fSQL_commit(db)

      ! close the obsDB file
      call fSQL_close(db, stat) 

    else
      write(*,*) 'odbf_insertInMidasBodyTable: the midas body output table already exists, ' // &
                 'proceed to join with temporary table that has new columns.'
    end if ! .not.midasTableExists

    ! open the obsDB file
    call fSQL_open(db, trim(fileName), status=stat)
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_insertInMidasBodyTable: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort('odbf_insertInMidasBodyTable: fSQL_open '//fSQL_errmsg(stat) )
    end if

    tableInsertColumnList = ''
    obsSpaceColIndexSourceArr(:) = mpc_missingValue_int
    call getCreateTableInsertQueries(numberUpdateItems, updateItemList, midasTableType, &
                                     queryCreateTable, queryInsertInTable, &
                                     tableInsertColumnList, obsSpaceColIndexSourceArr(:))

    ! Create a temporary table with new columns and values from obsSpaceData
    write(*,*) 'odbf_insertInMidasBodyTable: queryCreateTable   -->', trim(queryCreateTable)
    call fSQL_do_many(db, queryCreateTable, stat)
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'fSQL_do_many: ', fSQL_errmsg(stat)
      call utl_abort('odbf_insertInMidasBodyTable: Problem with fSQL_do_many')
    end if 

    ! Prepare to insert into the table
    write(*,*) 'odbf_insertInMidasBodyTable: queryInsertInTable -->', trim(queryInsertInTable)
    call fSQL_prepare(db, queryInsertInTable, stmt, stat)
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_insertInMidasBodyTable: fSQL_prepare: ', fSQL_errmsg(stat)
      call utl_abort('odbf_insertInMidasBodyTable: fSQL_prepare')
    end if

    ! updateItemList(:) columns are updated within the BODY2 loop
    call fSQL_begin(db)
    HEADER2: do headIndex = 1, obs_numHeader(obsdat)

      obsIdf = obs_headElem_i( obsdat, OBS_IDF, headIndex )
      if ( obsIdf /= fileIndex ) cycle HEADER2

      bodyIndexBegin = obs_headElem_i( obsdat, OBS_RLN, headIndex )
      bodyIndexEnd = bodyIndexBegin + &
                     obs_headElem_i( obsdat, OBS_NLV, headIndex ) - 1

      BODY2: do bodyIndex = bodyIndexBegin, bodyIndexEnd
      
        ! do not try to update if the observed value is missing
        obsValue = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex)
        if ( obsValue == obs_missingValue_R ) cycle BODY2

        obsIdd = obs_bodyPrimaryKey(obsdat, bodyIndex)
        call fSQL_bind_param(stmt, PARAM_INDEX=1, INT8_VAR=obsIdd)

        do updateItemIndex = 1, numberUpdateItems

          obsSpaceColIndexSource = obsSpaceColIndexSourceArr(updateItemIndex)
          columnParamIndex = updateItemIndex + 1

          ! update the value, but set to null if it is missing
          if (obs_columnDataType(obsSpaceColIndexSource) == 'real') then
            updateValue_r = obs_bodyElem_r(obsdat, obsSpaceColIndexSource, bodyIndex)

            ! change units for surface emissivity
            if (obsSpaceColIndexSource == OBS_SEM) then
              updateValue_r =updateValue_r * 100.0D0
            end if

            if ( updateValue_r == obs_missingValue_R ) then
              call fSQL_bind_param(stmt, PARAM_INDEX=columnParamIndex)  ! sql null values
            else
              call fSQL_bind_param(stmt, PARAM_INDEX=columnParamIndex, REAL8_VAR=updateValue_r)
            end if
          else
            updateValue_i = obs_bodyElem_i(obsdat, obsSpaceColIndexSource, bodyIndex)
            if ( updateValue_i == mpc_missingValue_int ) then
              call fSQL_bind_param(stmt, PARAM_INDEX=columnParamIndex)  ! sql null values
            else
              call fSQL_bind_param(stmt, PARAM_INDEX=columnParamIndex, INT_VAR=updateValue_i)
            end if
          end if ! if (obs_columnDataType(obsSpaceColIndexSource) == 'real')
        end do ! do updateItemIndex

        call fSQL_exec_stmt(stmt)

      end do BODY2

    end do HEADER2

    call fSQL_finalize(stmt)
    call fSQL_commit(db)

    ! close the obsDB file
    call fSQL_close(db, stat)

    call mergeTableInMidasTables(fileName, midasBodyTableName, obsBodyKeySqlName, &
                                 tableInsertColumnList)

    write(*,*)
    write(*,*) 'odbf_insertInMidasBodyTable: finished'
    write(*,*)

  end subroutine odbf_insertInMidasBodyTable

  !--------------------------------------------------------------------------
  ! odbf_createMidasHeaderTable
  !--------------------------------------------------------------------------
  subroutine odbf_createMidasHeaderTable(fileName)
    !
    ! :Purpose: Create the midasOutput Header table that stores all quantities computed
    !           in MIDAS at the level of the obsSpaceData Header table (e.g. ETOP, VTOP, ECF).
    !
    implicit none

    ! arguments:
    character(len=*),              intent(in)  :: fileName

    ! locals:
    character(len=3000)         :: query
    type(fSQL_STATUS)           :: stat ! sqlite error status
    type(fSQL_DATABASE)         :: db   ! sqlite file handle
    character(len=*), parameter :: myName = 'odbf_createMidasHeaderTable'

    ! open the obsDB file
    call fSQL_open(db, trim(fileName), status=stat)
    if ( fSQL_error(stat) /= FSQL_OK ) then
      call utl_abort(myName//': fSQL_open '//fSQL_errmsg(stat))
    end if

    ! create the new MIDAS table
    query = 'create table ' // trim(midasHeadTableName) // ' (' // new_line('A') // &
            '  ' // trim(midasHeadKeySqlName) // ' integer primary key,' // new_line('A') // &
            '  ' // trim(obsHeadKeySqlName) // ' integer' // new_line('A')

    query = trim(query) // ');'
    write(*,*) myName//': query = ', trim(query)
    call fSQL_do_many(db, query, stat)
    if ( fSQL_error(stat) /= FSQL_OK ) then
      call utl_abort(myName//': Problem with fSQL_do_many '//fSQL_errmsg(stat))
    end if

    ! close the obsDB file
    call fSQL_close(db, stat) 

  end subroutine odbf_createMidasHeaderTable

  !--------------------------------------------------------------------------
  ! odbf_createMidasBodyTable
  !--------------------------------------------------------------------------
  subroutine odbf_createMidasBodyTable(fileName)
    !
    ! :Purpose: Create the midasOutput table that stores all quantities computed
    !           in MIDAS at the level of the obsSpaceData Body table (e.g. OMP, OMA, FLG).
    !
    implicit none

    ! arguments:
    character(len=*),              intent(in)  :: fileName

    ! locals:
    integer :: columnIndex, obsColumnIndex
    character(len=3000)      :: query
    character(len=20)        :: sqlDataType
    type(fSQL_STATUS)        :: stat ! sqlite error status
    type(fSQL_DATABASE)      :: db   ! sqlite file handle
    character(len=*), parameter :: myName = 'odbf_createMidasBodyTable'

    ! open the obsDB file
    call fSQL_open(db, trim(fileName), status=stat)
    if ( fSQL_error(stat) /= FSQL_OK ) then
      call utl_abort(myName//': fSQL_open '//fSQL_errmsg(stat))
    end if

    ! create the new MIDAS table
    query = 'create table ' // trim(midasBodyTableName) // ' (' // new_line('A') // &
            '  ' // trim(midasBodyKeySqlName) // ' integer primary key,' // new_line('A') // &
            '  ' // trim(obsHeadKeySqlName) // ' integer,' // new_line('A') // &
            '  ' // trim(obsBodyKeySqlName) // ' integer );'

    write(*,*) myName//': query = ', trim(query)
    call fSQL_do_many(db, query, stat)
    if ( fSQL_error(stat) /= FSQL_OK ) then
      call utl_abort(myName//': Problem with fSQL_do_many '//fSQL_errmsg(stat))
    end if

    ! close the obsDB file
    call fSQL_close(db, stat) 

  end subroutine odbf_createMidasBodyTable

  !--------------------------------------------------------------------------
  ! obdf_clean
  !--------------------------------------------------------------------------
  subroutine obdf_clean(fileName, familyType)

    ! :Purpose: After the observational thinning procedure, this subroutine removes
    !           rows that are flagged as thinned in MIDAS_BODY_OUTPUT Table 
    !           the rows in the Report, Observation and MIDAS_HEADER_OUTPUT with corresponding 
    !           ID_Report and ID_Observation are also removed. 

    implicit none

    ! arguments
    character(len=*),  intent(in) :: fileName
    character(len=*),  intent(in) :: familyType

    ! locals:
    character(len = 512)        :: query
    type(fSQL_STATUS)           :: stat ! sqlite error status
    type(fSQL_DATABASE)         :: db   ! sqlite file handle
    type(fSQL_STATEMENT)        :: stmt ! precompiled sqlite statements
    integer                     :: nulnam , ierr, fnom, fclos
    character(len = lenSqlName) :: flgSqlName
    
    ! namelist variables
    logical, save               :: useVacuum ! choose to 'vacuum' the file after cleaning to reduce file size

    namelist/namObsDbClean/ useVacuum

    ! default value
    useVacuum = .false.

    if ( .not. utl_isNamelistPresent('namObsDbClean','./flnml') ) then
      if ( mmpi_myid == 0 ) then
        write(*,*) 'odbf_clean: namObsDbClean is missing in the namelist.'
        write(*,*) '            The default values will be taken.'
      end if
    else
      ! reading namelist variables
      nulnam  = 0
      ierr = fnom(nulnam ,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam , nml=namObsDbClean, iostat=ierr)
      if ( ierr /= 0 ) call utl_abort('obdf_clean: Error reading namelist')
      ierr = fclos(nulnam )
    end if
    if ( mmpi_myid == 0 ) write(*, nml=namObsDbClean)

    write(*,*)
    write(*,*) 'obdf_clean: Starting'
    write(*,*)
    write(*,*) 'obdf_clean: FileName   : ', trim(FileName)
    write(*,*) 'obdf_clean: FamilyType : ', FamilyType

    ! open the obsDB file
    call fSQL_open(db, trim(fileName), stat)
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'obdf_clean: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort('obdf_clean: fSQL_open')
    end if
    
    flgSqlName = odbf_midasTabColFromObsSpaceName('FLG', midasBodyNamesList)

    ! Mark for deletion all records with bit 11 (2048) set
    query = ' delete from '// trim(midasBodyTableName) //' where '// trim(flgSqlName) //' & 2048 =2048;'

    call fSQL_prepare(db, query, stmt, stat)
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'obdf_clean: fSQL_prepare: ', fSQL_errmsg(stat)
      call utl_abort('obdf_clean: fSQL_prepare')
    end if

    call fSQL_exec_stmt(stmt)
  
    query = 'create temporary table good_headers as select distinct '// trim(obsHeadKeySqlName) //' from '// trim(midasBodyTableName) //';'
    write(*,*) 'obdf_clean: query = ', trim(query)
    call fSQL_do_many(db, query, stat)

    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'fSQL_do_many: ', fSQL_errmsg(stat)
      call utl_abort('obdf_clean: Problem with fSQL_do_many')
    end if

    query = 'delete from '// trim(midasHeadTableName) //' where '// trim(obsHeadKeySqlName) // &
            ' not in ( select '// trim(obsHeadKeySqlName) //' from good_headers );'
    write(*,*) 'obdf_clean: query = ', trim(query)
    call fSQL_prepare(db, query, stmt, stat)
    
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'obdf_clean: fSQL_prepare: ', fSQL_errmsg(stat)
      call utl_abort('obdf_clean: fSQL_prepare')
    end if

    call fSQL_exec_stmt(stmt)

    query = 'delete from '// trim(bodyTableName) //' where '// trim(obsHeadKeySqlName) // &
            ' not in ( select '// trim(obsHeadKeySqlName) //'  from good_headers );'
    write(*,*) 'obdf_clean: query = ', trim(query)
    call fSQL_prepare(db, query, stmt, stat)

    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'obdf_clean: fSQL_prepare: ', fSQL_errmsg(stat)
      call utl_abort('obdf_clean: fSQL_prepare')
    end if

    call fSQL_exec_stmt(stmt)

    query = 'delete from '// trim(headTableName) //' where '// trim(obsHeadKeySqlName) // &
            ' not in ( select '// trim(obsHeadKeySqlName) //'  from good_headers );'
    write(*,*) 'obdf_clean: query = ', trim(query)
    call fSQL_prepare(db, query, stmt, stat)

    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'obdf_clean: fSQL_prepare: ', fSQL_errmsg(stat)
      call utl_abort('obdf_clean: fSQL_prepare')
    end if

    call fSQL_exec_stmt(stmt)
    call fSQL_finalize(stmt)

    ! Reduces the size of SQL file 
    if ( useVacuum ) then
      query = 'vacuum;'
      write(*,*) 'obdf_clean: query = ', trim(query)
      call fSQL_do_many(db, query, stat)

      if ( fSQL_error(stat) /= FSQL_OK ) then
        write(*,*) 'fSQL_do_many: ', fSQL_errmsg(stat)
        call utl_abort('obdf_clean: Problem with fSQL_do_many')
      end if
    end if

    call fSQL_close(db, stat)

  end subroutine obdf_clean

  !--------------------------------------------------------------------------
  ! getCreateTableInsertQueries
  !--------------------------------------------------------------------------
  subroutine getCreateTableInsertQueries(numberUpdateItems, updateItemList, midasTableType, &
                                         queryCreateTable, queryInsertInTable, &
                                         tableInsertColumnList, obsSpaceColIndexSourceArr)
    !
    ! :Purpose: Generate the queries for creating the table and insert columns into it.
    !
    implicit none

    ! arguments:
    integer, intent(in)          :: numberUpdateItems               ! number of items in update list
    character(len=4), intent(in) :: updateItemList(:)               ! update list
    character(len=*), intent(in) :: midasTableType                  ! table type: 'header' or 'body'
    character(len=*), intent(inout) :: queryCreateTable             ! query to create table
    character(len=*), intent(inout) :: queryInsertInTable           ! query to insert new columns in the table
    character(len=*), intent(inout) :: tableInsertColumnList        ! char of "combinedTableName.column1, combinedTableName.column2, .."
    integer         , intent(inout) :: obsSpaceColIndexSourceArr(:) ! list of obsSpaceData columnIndex for items in update list

    ! locals:
    integer :: updateItemIndex, ierr
    integer :: obsSpaceColIndexSource
    character(len=4)          :: obsSpaceColumnName
    character(len=20)         :: sqlDataType
    character(len=lenSqlName) :: sqlColumnName
    character(len=3000)       :: queryForValues


    do updateItemIndex = 1, numberUpdateItems

      ! get obsSpaceData column index for source of updated sql column
      obsSpaceColumnName = updateItemList(updateItemIndex)
      ierr = clib_toUpper(obsSpaceColumnName)
      obsSpaceColIndexSource = obs_columnIndexFromName(trim(obsSpaceColumnName))      

      if (midasTableType == 'header') then
        sqlColumnName = odbf_midasTabColFromObsSpaceName(updateItemList(updateItemIndex), midasHeadNamesList)
      else if (midasTableType == 'body') then
        sqlColumnName = odbf_midasTabColFromObsSpaceName(updateItemList(updateItemIndex), midasBodyNamesList)
      end if
      write(*,*) 'getCreateTableInsertQueries: updating midasTable column: ', trim(sqlColumnName)
      write(*,*) 'getCreateTableInsertQueries: with contents of obsSpaceData column: ', &
                 trim(obsSpaceColumnName)

      if (updateItemIndex == 1) then
        if (midasTableType == 'header') then
          queryCreateTable = 'create table '// trim(combinedTableName) // '(' // new_line('A') // &
                             '  ' // trim(obsHeadKeySqlName) // ' integer ' // new_line('A')

          queryInsertInTable = 'insert into '// trim(combinedTableName) // '(' // new_line('A') // &
                               '  ' // trim(obsHeadKeySqlName) // new_line('A')

        else if (midasTableType == 'body') then
          queryCreateTable = 'create table '// trim(combinedTableName) // '(' // new_line('A') // &
                             '  ' // trim(obsBodyKeySqlName) // ' integer ' // new_line('A')

          queryInsertInTable = 'insert into '// trim(combinedTableName) // '(' // new_line('A') // &
                               '  ' // trim(obsBodyKeySqlName) // new_line('A')                              
        end if ! if (midasTableType == 'header')

        queryForValues = 'values(?'
      end if ! if (updateItemIndex == 1)

      if (obs_columnDataType(obsSpaceColIndexSource) == 'real') then
        sqlDataType = 'double'
      else
        sqlDataType = 'integer'
      end if
      queryCreateTable = trim(queryCreateTable) // ', ' // trim(sqlColumnName) // ' ' // trim(sqlDataType) // new_line('A')
      queryInsertInTable = trim(queryInsertInTable) // ', ' // trim(sqlColumnName) // new_line('A')
      queryForValues = trim(queryForValues) // ',?'

      if (updateItemIndex == numberUpdateItems) then
        queryCreateTable = trim(queryCreateTable) // ');'
        queryForValues = trim(queryForValues) // ')'        
        queryInsertInTable = trim(queryInsertInTable) // ') ' // trim(queryForValues) // ';'
      end if

      tableInsertColumnList = trim(tableInsertColumnList) // ', '// &
                              trim(combinedTableName) // '.' // trim(sqlColumnName) // new_line('A')
      obsSpaceColIndexSourceArr(updateItemIndex) = obsSpaceColIndexSource
    end do ! do updateItemIndex

  end subroutine getCreateTableInsertQueries

  !--------------------------------------------------------------------------
  ! mergeTableInMidasTables
  !--------------------------------------------------------------------------
  subroutine mergeTableInMidasTables(fileName, midasTableName, jointColumnName, &
                                     tableInsertColumnList)
    !
    ! :Purpose: In a series of join/drop/alter merge input table and midasTable
    !           to create a new midasTable which contains the original columns plus 
    !           columns from the input table.   
    !
    implicit none

    ! arguments:
    character(len=*), intent(in) :: fileName              ! obsDB filename
    character(len=*), intent(in) :: midasTableName        ! name of original midas table to add column to
    character(len=*), intent(in) :: jointColumnName       ! name of column used to match original midas table and temporary table
    character(len=*), intent(in) :: tableInsertColumnList ! char of "combinedTableName.column1, combinedTableName.column2, .." to add to original midas table

    ! locals:
    type(fSQL_STATUS)   :: stat ! sqlite error status
    type(fSQL_DATABASE) :: db   ! sqlite file handle
    integer             :: columnIndex
    character(len=3000) :: query

    write(*,*)
    write(*,*) 'mergeTableInMidasTables: START'
    write(*,*)

    ! open the obsDB file
    call fSQL_open(db, trim(fileName), stat)
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'mergeTableInMidasTables: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort('mergeTableInMidasTables: fSQL_open')
    end if

    ! Combine midasTable + inputTable -> table_tmp
    query = 'create table table_tmp as select ' // trim(midasTableName) // '.*' // new_line('A')

    query = trim(query) // trim(tableInsertColumnList) // '  from ' // new_line('A') // &
            '  ' // trim(combinedTableName) // ' inner join ' // trim(midasTableName) // ' on ' // new_line('A') // &
            '  ' // trim(combinedTableName) // '.' // trim(jointColumnName) // '=' // &
            trim(midasTableName) // '.' // trim(jointColumnName) // ';'
   
    write(*,*) 'mergeTableInMidasTables: query ---> ', trim(query)
    call fSQL_do_many(db, query, stat)
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'fSQL_do_many: ', fSQL_errmsg(stat)
      call utl_abort('mergeTableInMidasTables: Problem with fSQL_do_many')
    end if

    ! Drop the midasTable
    query = 'drop table ' // trim(midasTableName) // ';'
      
    write(*,*) 'mergeTableInMidasTables: query ---> ', trim(query)
    call fSQL_do_many(db, query, stat)
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'fSQL_do_many: ', fSQL_errmsg(stat)
      call utl_abort('mergeTableInMidasTables: Problem with fSQL_do_many')
    end if

    ! Drop the inputTable
    query = 'drop table ' // trim(combinedTableName) // ';'
    write(*,*) 'mergeTableInMidasTables: query ---> ', trim(query)
    call fSQL_do_many(db, query, stat)
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'fSQL_do_many: ', fSQL_errmsg(stat)
      call utl_abort('mergeTableInMidasTables: Problem with fSQL_do_many')
    end if

    ! Rename table_tmp -> midasTable
    query = 'alter table table_tmp rename to ' // trim(midasTableName) // ';'

    write(*,*) 'mergeTableInMidasTables: query ---> ', trim(query)
    call fSQL_do_many(db, query, stat)
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'fSQL_do_many: ', fSQL_errmsg(stat)
      call utl_abort('mergeTableInMidasTables: Problem with fSQL_do_many')
    end if

    ! close the obsDB file
    call fSQL_close(db, stat) 

    write(*,*)
    write(*,*) 'mergeTableInMidasTables: END'
    write(*,*)

  end subroutine mergeTableInMidasTables

end module obsdbFiles_mod
