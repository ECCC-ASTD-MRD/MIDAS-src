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

module obsdbFiles_mod
  ! MODULE obsdbFiles (prefix='odbf' category='6. Observation input/output')
  !
  ! :Purpose: To read and update sqlite files that are in the new 'obsDB' format.
  !
  use mpi_mod
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

  implicit none
  save
  private

  ! Public subroutines and functions:
  public :: odbf_isActive, odbf_readFile, odbf_updateFile

 
  ! Arrays used to match obsDB column names with obsSpaceData column names

  integer, parameter :: lenSqlName    = 60
  integer, parameter :: sqlColIndex   = 1
  integer, parameter :: obsColIndex   = 2
  integer, parameter :: varNoColIndex = 2
  
  character(len=lenSqlName) :: headTableName
  character(len=lenSqlName) :: bodyTableName

  character(len=lenSqlName) :: midasTableName = 'midasOutput'

  ! ...for the header table

  integer :: numHeadMatch 
  character(len=lenSqlName) :: headKeySqlName
  character(len=lenSqlName) :: headDateSqlName 
  character(len=lenSqlName), allocatable :: headMatchList(:,:) 
  character(len=6), allocatable :: headBufrList(:)

  ! ...for the body table
  
  integer :: numBodyMatch
  character(len=lenSqlName) :: bodyKeySqlName
  character(len=lenSqlName), allocatable :: bodyMatchList(:,:)
  character(len=6), allocatable :: bodyBufrList(:)

  ! Dictionary of 'varno' value for each obsDB observation value column
  integer, parameter :: numVarNo = 1
  character(len=lenSqlName) :: varNoList(2,numVarNo)
  
  ! Column names for the MIDAS output table and corresponding obsSpace names
  character(len=lenSqlName) :: midasKeySqlName = 'ID_MIDAS'
  integer, parameter :: numColMidasTable = 11
  integer, parameter :: numColMidasTableRequired = 3
  character(len=lenSqlName) :: midasOutputNamesList(2,numColMidasTable) = (/ &
       'vcoord',             'PPP',  &
       'varNo',              'VNM',  &
       'obsValue',           'VAR',  &
       'flag',               'FLG',  &
       'obsMinusBackground', 'OMP',  &
       'obsMinusAnalysis',   'OMA',  &
       'obsError',           'OER',  &
       'backgroundError',    'HPHT', &
       'fso',                'FSO',  &
       'biasCorrection',     'BCOR', &
       'sfcEmissivity',      'SEM'/)

  ! Other constants
  logical, parameter :: setObsFlagZero = .true.

  ! NAMELIST variables
  logical :: obsDbActive
  integer :: numElemIdList
  integer :: elemIdList(100)

contains

  !--------------------------------------------------------------------------
  ! odbf_setup
  !--------------------------------------------------------------------------
  subroutine odbf_setup()
    !
    ! :Purpose: Read the namelist for obsDB files
    !
    implicit none

    ! locals
    integer            :: nulnam, ierr
    integer, external  :: fnom, fclos
    logical, save      :: nmlAlreadyRead = .false.
    
    !  to be changed - used for testing only
    character(len=lenSqlName) :: obsdb_column_file                              ! Temporary Directory to 
                                                                         ! ObsDBColumnTable.dat' 

    character(len=lenSqlName) :: readLine
    character(len=lenSqlName) :: readDBColumn, readMidasColumn, readBufrColumn  ! Strings to temporary assign 
                                                                         ! read header and body column names
    
    integer            :: headerTableRow, bodyTableRow                   ! Counters to determine the size of 
                                                                         ! header and body table
    
    integer            :: countRow, countHeadMatchRow, countBodyMatchRow ! Counters to write to 
                                                                         ! headMatchList and bodyMatchList

    namelist /namobsdb/ obsDbActive, numElemIdList, elemIdList

    if ( nmlAlreadyRead ) return

    nmlAlreadyRead = .true.

    ! default values
    obsDbActive = .false.
    numElemIdList = 0
    elemIdList(:) = 0

    if ( .not. utl_isNamelistPresent('NAMOBSDB','./flnml') ) then
      if ( mpi_myid == 0 ) then
        write(*,*) 'odbf_setup: namObsDB is missing in the namelist.'
        write(*,*) '            The default values will be taken.'
      end if
    else
      ! reading namelist variables
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam, nml=namobsdb, iostat=ierr)
      if ( ierr /= 0 ) call utl_abort('odbf_setup: Error reading namelist')
      ierr = fclos(nulnam)
    end if
    if ( mpi_myid == 0 ) write(*, nml=namObsDb)

    if (obsDbActive .and. numElemIdList==0) then
      call utl_abort('odbf_setup: element list is empty')
    end if

  !initialize obsDb columns names to be consistent with MIDAS obsSpaceData Report column names
  nulnam = 0
    
  obsdb_column_file='/home/zqw002/midas-work/src/modules/ObsDBColumnTable.dat'
   
  ierr = fnom(nulnam,trim(obsdb_column_file),'FTN+SEQ+R/O',0)
  if ( ierr /= 0 ) call utl_abort('odbf_setup: Error reading ObsDBColumnTable file') 
  
  !get number of rows in obsdb coloumn summary file(obsdb_column_file)
  headerTableRow = 0
  bodyTableRow = 0
  do while(ierr==0)
    read(nulnam, '(A)' ,iostat=ierr) readline

    !Search for the start of header table in the file
    if(index(trim(readline),'HEADER TABLE INFO BEGINS')>0)then

      !read the header table  
      do while(ierr==0)
        read(nulnam,'(A)' ,iostat=ierr) readline
        !Stop search at the end of header table
        if(index(trim(readline),'HEADER TABLE INFO ENDS')>0) exit
        headerTableRow = headerTableRow +1
      end do

    !Search for start of the body table
    else if(index(trim(readline),'BODY TABLE INFO BEGINS')>0)then

      !read the body table 
      do while(ierr==0)
        read(nulnam,'(A)' ,iostat=ierr) readline
        !Stop search at the end of body table
        if(index(trim(readline),'BODY TABLE INFO ENDS')>0) exit
        bodyTableRow = bodyTableRow + 1
      end do
    end if
  end do
  

  numHeadMatch = headerTableRow -3 ! 3 rows in header table are not included in the numHeadMatch
  numBodyMatch = bodyTableRow -2 ! 2 rows in body table are not included in the numHeadMatch

  ! Read Report obsdb and MIDAS obsSpaceData columns name
  allocate(headMatchList(2,numHeadMatch))
  allocate(headBufrList(numHeadMatch))
  allocate(bodyMatchList(2,numBodyMatch))
  allocate(bodyBufrList(numBodyMatch))

  rewind(nulnam)

  ierr = 0
  do while(ierr==0)
    read(nulnam, '(A)' ,iostat=ierr) readline

    !Search for the start of HEADER table in the file
    if(index(trim(readline),'HEADER TABLE INFO BEGINS')>0)then
          
      countHeadMatchRow = 0
      !Read all header table rows
      do countRow = 1, headerTableRow
        read(nulnam,*,iostat=ierr) readDBColumn, readMidasColumn, readBufrColumn      
        ! Assign read values to appropriate variable
        select case(trim(readMidasColumn))
          case('headTableName')
            headTableName = trim(readDBColumn)
          case('headPrimaryKey')
            headKeySqlName = trim(readDBColumn)
          case('DAT;ETM')
            headDateSqlName = trim(readDBColumn)
          case default
            countHeadMatchRow = countHeadMatchRow +1
            headMatchList(1,countHeadMatchRow) = trim(readDBColumn)
            headMatchList(2,countHeadMatchRow) = trim(readMidasColumn)
            headBufrList(countHeadMatchRow) = trim(readBufrColumn)
        end select
      end do

    !Search for the start of BODY table in the file
    else if(index(trim(readline),'BODY TABLE INFO BEGINS')>0)then

      countBodyMatchRow = 0
      !Read all body table rows
      do countRow = 1, bodyTableRow
        read(nulnam,*,iostat=ierr) readDBColumn, readMidasColumn, readBufrColumn      
        ! Assign read values to appropriate variable
        select case(trim(readMidasColumn))
          case('bodyTableName')
            bodyTableName = trim(readDBColumn)
          case('bodyPrimaryKey')
            bodyKeySqlName = trim(readDBColumn)
          case default
            countBodyMatchRow = countBodyMatchRow +1
            bodyMatchList(1,countBodyMatchRow) = trim(readDBColumn)
            bodyMatchList(2,countBodyMatchRow) = trim(readMidasColumn)
            bodyBufrList(countBodyMatchRow) = trim(readBufrColumn)

            !Assign varNoList if appropriate
            if(bodyMatchList(2,countBodyMatchRow)=='VAR') then
              varNoList(1, 1) = bodyMatchList(1,countBodyMatchRow) 
              varNoList(2, 1) = bodyBufrList(countBodyMatchRow)
            end if
        end select
      end do 
    end if
  end do

  ierr = fclos(nulnam)

  end subroutine odbf_setup

  !--------------------------------------------------------------------------
  ! odbf_isActive
  !--------------------------------------------------------------------------
  function odbf_isActive() result(isActive)
    !
    ! :Purpose: Tell the caller if the namelist indicates obsDB is active
    !
    implicit none

    ! return value
    logical :: isActive

    call odbf_setup()

    isActive = obsDbActive

  end function odbf_isActive

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

    call odbf_getColumnValuesDate(headDateValues, headTimeValues, fileName=trim(fileName), &
                                  tableName=headTableName, sqlColumnName=headDateSqlName)
    call odbf_getColumnValuesChar(headCharValues, fileName=trim(fileName), &
                                  tableName=headTableName, sqlColumnNames=headCharSqlNames)
    call odbf_getColumnValuesNum (headValues, fileName=trim(fileName), &
                                  tableName=headTableName, sqlColumnNames=headSqlNames)
    call odbf_getColumnValuesNum (bodyValues, fileName=trim(fileName), &
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

    !- 1.2 Copy values from local tables into obsSpaceData

    ! Starting point for adding rows to obsSpaceData
    bodyIndexBegin = obs_numBody(obsdat) + 1
    headIndexBegin = obs_numHeader(obsdat) + 1

    ! Set the columns related to surface type
    call odbf_setSurfaceType(obsdat, headIndexBegin, fileName=trim(fileName), &
                             tableName=headTableName)

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

    !- 1.3 Set some other quantities in obsSpaceData Header table

    do headIndex = headIndexBegin, headIndexEnd
      call obs_headSet_i(obsdat, OBS_SEN, headIndex, nint(MPC_missingValue_R8))
      call obs_headSet_i(obsdat, OBS_ONM, headIndex, headIndex)
      call obs_headSet_i(obsdat, OBS_OTP, headIndex, fileIndex)
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

    !- 1.4 Read values written during previous MIDAS program executions

    call odbf_readMidasTable(obsdat, trim(fileName), familyType, fileIndex)

    !- 2.0 Additional changes to values after they are in obsSpaceData

    call odbf_adjustValues(obsdat, headIndexBegin, headIndexEnd)

    if ( trim(familyType) /= 'TO' ) then
      call ovt_transformObsValues      (obsdat, headIndexBegin, headIndexEnd )
      call ovt_adjustHumGZ             (obsdat, headIndexBegin, headIndexEnd )
      call obsu_computeVertCoordSurfObs(obsdat, headIndexBegin, headIndexEnd )
    end if

    ! Set a bunch of obsSpaceData body columns to 'missing'
    do bodyIndex = bodyIndexBegin, bodyIndexEnd
      if ( obs_columnActive_RB(obsdat, OBS_OMA) )  call obs_bodySet_r(obsdat, OBS_OMA , bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_OMA0))  call obs_bodySet_r(obsdat, OBS_OMA0, bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_OMP) )  call obs_bodySet_r(obsdat, OBS_OMP , bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_OMP6))  call obs_bodySet_r(obsdat, OBS_OMP6, bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_OER) )  call obs_bodySet_r(obsdat, OBS_OER , bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_HPHT) ) call obs_bodySet_r(obsdat, OBS_HPHT, bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_HAHT))  call obs_bodySet_r(obsdat, OBS_HAHT, bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_WORK) ) call obs_bodySet_r(obsdat, OBS_WORK, bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_SIGI))  call obs_bodySet_r(obsdat, OBS_SIGI, bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_SIGO))  call obs_bodySet_r(obsdat, OBS_SIGO, bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_ZHA ))  call obs_bodySet_r(obsdat, OBS_ZHA , bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_SEM ))  call obs_bodySet_r(obsdat, OBS_SEM , bodyIndex, obs_missingValue_R)
      if ( obs_columnActive_RB(obsdat, OBS_BCOR))  call obs_bodySet_r(obsdat, OBS_BCOR, bodyIndex, obs_missingValue_R)
    end do

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
  ! odbf_readMidasTable
  !--------------------------------------------------------------------------
  subroutine odbf_readMidasTable(obsdat, fileName, familyType, fileIndex)
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
    logical              :: midasTableExists
    logical, allocatable :: midasColumnExists(:)

    write(*,*)
    write(*,*) 'odbf_readMidasTable: Starting'
    write(*,*)
    write(*,*) 'odbf_readMidasTable: FileName   : ', trim(FileName)
    write(*,*) 'odbf_readMidasTable: FamilyType : ', FamilyType

    ! check if midasTable already exists in the file
    midasTableExists = sqlu_sqlTableExists(fileName, midasTableName)

    if (.not. midasTableExists) then
      write(*,*) 'odbf_readMidasTable: MIDAS table not present in file'
      return
    else
      write(*,*) 'odbf_readMidasTable: MIDAS table present in file, will read contents'
    end if

    ! some sql column names
    vnmSqlName = odbf_midasTabColFromObsSpaceName('VNM')
    pppSqlName = odbf_midasTabColFromObsSpaceName('PPP')
    varSqlName = odbf_midasTabColFromObsSpaceName('VAR')

    ! check which columns exist in the MIDAS output table
    allocate(midasColumnExists(numColMidasTable))
    do sqlNameIndex = 1, numColMidasTable
      sqlColumnName = midasOutputNamesList(1,sqlNameIndex)
      midasColumnExists(sqlNameIndex) = &
           sqlu_sqlColumnExists(fileName, midasTableName, sqlColumnName)
    end do

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_readMidasTable: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_readMidasTable: fSQL_open' )
    end if
    
    ! read the contents of the MIDAS table, one column at a time
    SQLNAME: do sqlNameIndex = 1, numColMidasTable

      ! skip this sql column name if it is not present in the file
      if (.not. midasColumnExists(sqlNameIndex)) cycle SQLNAME

      ! get obsSpaceData column name and index corresponding to sql column
      sqlColumnName      = midasOutputNamesList(1,sqlNameIndex)
      obsSpaceColumnName = midasOutputNamesList(2,sqlNameIndex)
      ierr = clib_toUpper(obsSpaceColumnName)
      obsSpaceColIndex = obs_columnIndexFromName(trim(obsSpaceColumnName))

      ! skip this sql column name if the related obsSpaceData column is not active
      if (obs_columnDataType(obsSpaceColIndex) == 'real') then
        if (.not. obs_columnActive_RB(obsdat, obsSpaceColIndex)) cycle SQLNAME
      else
        if (.not. obs_columnActive_IB(obsdat, obsSpaceColIndex)) cycle SQLNAME
      end if

      write(*,*) 'odbf_readMidasTable: reading midasTable column: ', &
                 trim(sqlColumnName)
      write(*,*) 'odbf_readMidasTable: to update obsSpaceData column: ', &
                 trim(obsSpaceColumnName)

      ! prepare sql update query
      query = 'select ' // trim(sqlColumnName) // &
              ' from ' // trim(midasTableName) // &
              ' where ' // &
              trim(bodyKeySqlName) // ' = ? and '   // &
              trim(vnmSqlName)     // ' = ? and '   // &
              trim(pppSqlName)     // ' = ? ;'
      write(*,*) 'odbf_readMidasTable: query ---> ', trim(query)

      call fSQL_prepare( db, query , stmt, stat )
      if ( fSQL_error(stat) /= FSQL_OK ) then
        write(*,*) 'odbf_readMidasTable: fSQL_prepare: ', fSQL_errmsg(stat)
        call utl_abort( 'odbf_readMidasTable: fSQL_prepare' )
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
          obsIdd  = obs_bodyPrimaryKey( obsdat, bodyIndex )
          call fSQL_bind_param(stmt, PARAM_INDEX=1, INT8_VAR=obsIdd)
          obsVarNo  = obs_bodyElem_i( obsdat, obs_vnm, bodyIndex )
          call fSQL_bind_param(stmt, PARAM_INDEX=2, INT_VAR=obsVarNo)
          obsPPP  = obs_bodyElem_r( obsdat, obs_ppp, bodyIndex )
          call fSQL_bind_param(stmt, PARAM_INDEX=3, REAL8_VAR=obsPPP)
          call fSQL_exec_stmt(stmt)

          ! read the real or integer value
          if (obs_columnDataType(obsSpaceColIndex) == 'real') then
            call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                 mode=FSQL_REAL8, status=stat )
            if ( fSQL_error(stat) /= FSQL_OK ) then
              write(*,*) 'odbf_readMidasTable: fSQL_get_many: ', fSQL_errmsg(stat)
              call utl_abort('odbf_readMidasTable: problem with fSQL_get_many')
            end if
            if (numRows /= 1 .or. numColumns /= 1) then
              write(*,*) 'odbf_readMidasTable: numRows, numColumns =', numRows, numColumns
              call utl_abort('odbf_readMidasTable: sql query did not return 1 value')
            end if
            call fSQL_fill_matrix ( stmt, matdata_r8 )
            call obs_bodySet_r(obsdat,obsSpaceColIndex,bodyIndex,matdata_r8(1,1))
          else
            call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                                mode=FSQL_INT, status=stat )
            if ( fSQL_error(stat) /= FSQL_OK ) then
              write(*,*) 'odbf_readMidasTable: fSQL_get_many: ', fSQL_errmsg(stat)
              call utl_abort('odbf_readMidasTable: problem with fSQL_get_many')
            end if
            if (numRows /= 1 .or. numColumns /= 1) then
              write(*,*) 'odbf_readMidasTable: numRows, numColumns =', numRows, numColumns
              call utl_abort('odbf_readMidasTable: sql query did not return 1 value')
            end if
            call fSQL_fill_matrix ( stmt, matdata_int )
            call obs_bodySet_i(obsdat,obsSpaceColIndex,bodyIndex,matdata_int(1,1))
          end if

          call fSQL_free_mem( stmt )

        end do BODY2

      end do HEADER2

      call fSQL_finalize( stmt )
      call fSQL_commit( db )

    end do SQLNAME

    ! close the obsDB file
    call fSQL_close( db, stat ) 

    deallocate(midasColumnExists)

    write(*,*)
    write(*,*) 'odbf_readMidasTable: finished'
    write(*,*)

  end subroutine odbf_readMidasTable

  !--------------------------------------------------------------------------
  ! odbf_updateFile
  !--------------------------------------------------------------------------
  subroutine odbf_updateFile(obsdat, fileName, familyType, fileIndex)
    !
    ! :Purpose: Update the selected quantities in an obsDB file using
    !           values from obsSpaceData. If the MIDAS table does not already
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
    integer              :: obsIdf, obsVarNo, midasKey, updateItemIndex, updateValue_i
    integer              :: headIndex, bodyIndex, bodyIndexBegin, bodyIndexEnd
    integer              :: obsSpaceColIndexSource, fnom, fclos, nulnam, ierr
    real(8)              :: updateValue_r, obsValue, obsPPP, obsVAR
    character(len=4)     :: obsSpaceColumnName
    character(len=lenSqlName) :: sqlColumnName, vnmSqlName, pppSqlName, varSqlName
    character(len=3000)  :: query
    character(len=20)    :: sqlDataType
    logical              :: midasTableExists
    logical, allocatable :: midasColumnExists(:)
    logical, save        :: nmlAlreadyRead = .false.

    ! namelist variables
    integer,          save :: numberUpdateItems  ! number of items to use from the list
    character(len=4), save :: updateItemList(15) ! obsSpace column names used to update the file

    namelist/namObsDbUpdate/ numberUpdateItems, updateItemList

    call tmg_start(97,'obdf_updateFile')

    write(*,*)
    write(*,*) 'odbf_updateFile: Starting'
    write(*,*)
    write(*,*) 'odbf_updateFile: FileName   : ', trim(FileName)
    write(*,*) 'odbf_updateFile: FamilyType : ', FamilyType

    if (.not. nmlAlreadyRead) then
      nmlAlreadyRead = .true.

      ! set default values of namelist variables
      updateItemList(:) = ''
      numberUpdateItems = 0

      ! Read the namelist for directives
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam, nml=namObsDbUpdate, iostat=ierr)
      if ( ierr /= 0 ) call utl_abort('odbf_updateFile: Error reading namelist')
      ierr = fclos(nulnam)

      ! Add "FLG" to the updateItemList to ensure it is always updated
      numberUpdateItems = numberUpdateItems + 1
      updateItemList(numberUpdateItems) = 'FLG'

      if ( mpi_myid == 0 ) then
        write(*,*) 'odbf_updateFile: NOTE: the FLG column is always added to update list'
        write(*, nml=namObsDbUpdate)
      end if
    end if ! not nmlAlreadyRead

    ! some sql column names
    vnmSqlName = odbf_midasTabColFromObsSpaceName('VNM')
    pppSqlName = odbf_midasTabColFromObsSpaceName('PPP')
    varSqlName = odbf_midasTabColFromObsSpaceName('VAR')

    ! check if midasTable already exists in the file
    midasTableExists = sqlu_sqlTableExists(fileName, midasTableName)

    if (.not. midasTableExists) then
      ! create midasTable by copying rearranging contents of observation table
      call odbf_createMidasTable(fileName)

      ! open the obsDB file
      call fSQL_open( db, trim(fileName), stat )
      if ( fSQL_error(stat) /= FSQL_OK ) then
        write(*,*) 'odbf_updateFile: fSQL_open: ', fSQL_errmsg(stat)
        call utl_abort( 'odbf_updateFile: fSQL_open' )
      end if

      ! set the primary key, keys to main obsDB tables and other basic info
      midasKey = 0
      query = 'insert into ' // trim(midasTableName) // '(' // &
              trim(midasKeySqlName) // ',' // trim(headKeySqlName) // ',' // &
              trim(bodyKeySqlName)  // ',' // trim(vnmSqlname)     // ',' // &
              trim(pppSqlName)      // ',' // trim(varSqlname)     // &
              ') values(?,?,?,?,?,?);'
      write(*,*) 'odbf_updateFile: query = ', trim(query)
      call fSQL_prepare( db, query, stmt, stat )
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
          obsIdd  = obs_bodyPrimaryKey( obsdat, bodyIndex )
          call fSQL_bind_param(stmt, PARAM_INDEX=3, INT8_VAR=obsIdd)
          obsVarNo  = obs_bodyElem_i( obsdat, obs_vnm, bodyIndex )
          call fSQL_bind_param(stmt, PARAM_INDEX=4, INT_VAR=obsVarNo)
          obsPPP  = obs_bodyElem_r( obsdat, obs_ppp, bodyIndex )
          call fSQL_bind_param(stmt, PARAM_INDEX=5, REAL8_VAR=obsPPP)
          obsVAR  = obs_bodyElem_r( obsdat, obs_var, bodyIndex )
          call fSQL_bind_param(stmt, PARAM_INDEX=6, REAL8_VAR=obsVAR)

          call fSQL_exec_stmt(stmt)

        end do BODY

      end do HEADER

      call fSQL_finalize( stmt )
      call fSQL_commit( db )

      ! create an index for the new table - necessary to speed up the update
      query = 'create index idx_midasTable on ' // &
              trim(midasTableName) // &
              '(' // trim(bodyKeySqlName) // ',' // trim(vnmSqlName) // ');'
      write(*,*) 'odbf_updateFile: query = ', trim(query)
      call fSQL_do_many( db, query, stat )
      if ( fSQL_error(stat) /= FSQL_OK ) then
        write(*,*) 'fSQL_do_many: ', fSQL_errmsg(stat)
        call utl_abort('odbf_updateFile: Problem with fSQL_do_many')
      end if

      ! close the obsDB file
      call fSQL_close( db, stat ) 

    else

      write(*,*) 'odbf_updateFile: the midas output table already exists, ' // &
                 'will just update its values'

    end if ! .not.midasTableExists

    ! check which columns to be updated already exist in the table
    allocate(midasColumnExists(numberUpdateItems))
    do updateItemIndex = 1, numberUpdateItems
      sqlColumnName = odbf_midasTabColFromObsSpaceName(updateItemList(updateItemIndex))
      midasColumnExists(updateItemIndex) = &
           sqlu_sqlColumnExists(fileName, midasTableName, sqlColumnName)
    end do

    ! now that the table exists, we can update the selected columns

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_updateFile: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_updateFile: fSQL_open' )
    end if
    
    ! updating the contents of the MIDAS table, one column at a time
    do updateItemIndex = 1, numberUpdateItems

      ! get obsSpaceData column index for source of updated sql column
      obsSpaceColumnName = updateItemList(updateItemIndex)
      ierr = clib_toUpper(obsSpaceColumnName)
      obsSpaceColIndexSource = obs_columnIndexFromName(trim(obsSpaceColumnName))

      sqlColumnName = odbf_midasTabColFromObsSpaceName(updateItemList(updateItemIndex))
      write(*,*) 'odbf_updateFile: updating midasTable column: ', trim(sqlColumnName)
      write(*,*) 'odbf_updateFile: with contents of obsSpaceData column: ', &
                 trim(obsSpaceColumnName)

      ! add column to sqlite table
      if (.not. midasColumnExists(updateItemIndex)) then
        if (obs_columnDataType(obsSpaceColIndexSource) == 'real') then
          sqlDataType = 'double'
        else
          sqlDataType = 'integer'
        end if
        query = 'alter table ' // trim(midasTableName) // ' add column ' // &
                trim(sqlColumnName) // ' ' // trim(sqlDataType) // ';'
        write(*,*) 'odbf_updateFile: query ---> ', trim(query)
        call fSQL_do_many( db, query, stat )
        if ( fSQL_error(stat) /= FSQL_OK ) then
          write(*,*) 'fSQL_do_many: ', fSQL_errmsg(stat)
          call utl_abort('odbf_updateFile: Problem with fSQL_do_many')
        end if
      end if

      ! prepare sql update query
      query = 'update ' // trim(midasTableName) // ' set ' // &
              trim(sqlColumnName)  // ' = ? where ' // &
              trim(bodyKeySqlName) // ' = ? and '   // &
              trim(vnmSqlName)     // ' = ? and '   // &
              trim(pppSqlName)     // ' = ? ;'
      write(*,*) 'odbf_updateFile: query ---> ', trim(query)

      call fSQL_prepare( db, query , stmt, stat )
      if ( fSQL_error(stat) /= FSQL_OK ) then
        write(*,*) 'odbf_updateFile: fSQL_prepare: ', fSQL_errmsg(stat)
        call utl_abort( 'odbf_updateFile: fSQL_prepare' )
      end if

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

          ! update the value, but set to null if it is missing
          if (obs_columnDataType(obsSpaceColIndexSource) == 'real') then
            updateValue_r = obs_bodyElem_r(obsdat, obsSpaceColIndexSource, bodyIndex)

            ! change units for surface emissivity
            if (obsSpaceColIndexSource == OBS_SEM) then
              updateValue_r =updateValue_r * 100.0D0
            end if

            if ( updateValue_r == obs_missingValue_R ) then
              call fSQL_bind_param(stmt, PARAM_INDEX=1)  ! sql null values
            else
              call fSQL_bind_param(stmt, PARAM_INDEX=1, REAL8_VAR=updateValue_r)
            end if
          else
            updateValue_i = obs_bodyElem_i(obsdat, obsSpaceColIndexSource, bodyIndex)
            if ( updateValue_i == mpc_missingValue_int ) then
              call fSQL_bind_param(stmt, PARAM_INDEX=1)  ! sql null values
            else
              call fSQL_bind_param(stmt, PARAM_INDEX=1, INT_VAR=updateValue_i)
            end if
          end if

          obsIdd  = obs_bodyPrimaryKey( obsdat, bodyIndex )
          call fSQL_bind_param(stmt, PARAM_INDEX=2, INT8_VAR=obsIdd)
          obsVarNo  = obs_bodyElem_i( obsdat, obs_vnm, bodyIndex )
          call fSQL_bind_param(stmt, PARAM_INDEX=3, INT_VAR=obsVarNo)
          obsPPP  = obs_bodyElem_r( obsdat, obs_ppp, bodyIndex )
          call fSQL_bind_param(stmt, PARAM_INDEX=4, REAL8_VAR=obsPPP)

          call fSQL_exec_stmt(stmt)

        end do BODY2

      end do HEADER2

      call fSQL_finalize( stmt )
      call fSQL_commit( db )
    end do

    ! close the obsDB file
    call fSQL_close( db, stat ) 

    deallocate(midasColumnExists)

    write(*,*)
    write(*,*) 'odbf_updateFile: finished'
    write(*,*)

    call tmg_stop(97)

  end subroutine odbf_updateFile

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
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getPrimaryKeys: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_getPrimaryKeys: fSQL_open' )
    end if

    ! build the sqlite query for the HEADER primary key
    query = 'select ' // trim(headKeySqlName) // ' from ' // &
            trim(headTableName) // ';'
    write(*,*) 'odbf_getPrimaryKeys: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query), stmt, status=stat )
    ! note: "status" not set when getting integers
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_INT8 )
    write(*,*) 'odbf_getPrimaryKeys: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( headPrimaryKey(numRows) )
    allocate( tempHeadKey(numRows,1) )
    call fSQL_fill_matrix( stmt, tempHeadKey )
    headPrimaryKey(:) = tempHeadKey(:,1)
    deallocate(tempHeadKey)

    ! build the sqlite query for the BODY primary key
    query = 'select ' // trim(bodyKeySqlName) // ' from ' // &
            trim(bodyTableName) // ';'
    write(*,*) 'odbf_getPrimaryKeys: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query) , stmt, status=stat )
    ! note: "status" not set when getting integers
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_INT8 )
    write(*,*) 'odbf_getPrimaryKeys: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( bodyPrimaryKey(numRows) )
    allocate( tempBodyKey(numRows,1) )
    call fSQL_fill_matrix( stmt, tempBodyKey )
    bodyPrimaryKey(:) = tempBodyKey(:,1)
    deallocate(tempBodyKey)

    ! build the sqlite query for the BODY-HEAD key
    query = 'select ' // trim(headKeySqlName) // ' from ' // &
            trim(bodyTableName) // ';'
    write(*,*) 'odbf_getPrimaryKeys: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query) , stmt, status=stat )
    ! note: "status" not set when getting integers
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_INT8 )
    write(*,*) 'odbf_getPrimaryKeys: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( bodyHeadKey(numRows) )
    allocate( tempBodyKey(numRows,1) )
    call fSQL_fill_matrix( stmt, tempBodyKey )
    bodyHeadKey(:) = tempBodyKey(:,1)
    deallocate(tempBodyKey)

    ! close the obsDB file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine odbf_getPrimaryKeys

  !--------------------------------------------------------------------------
  ! odbf_getColumnValuesDate
  !--------------------------------------------------------------------------
  subroutine odbf_getColumnValuesDate(columnDateValues, columnTimeValues, fileName, &
                                      tableName, sqlColumnName)
    !
    ! :Purpose: Read the column values from obsDB file for the specified table
    !           and column names.
    !
    implicit none

    ! arguments:
    integer, allocatable, intent(out) :: columnDateValues(:)
    integer, allocatable, intent(out) :: columnTimeValues(:)
    character(len=*),     intent(in)  :: sqlColumnName
    character(len=*),     intent(in)  :: fileName
    character(len=*),     intent(in)  :: tableName

    ! locals:
    integer              :: numRows, numColumns, rowIndex
    character(len=20), allocatable :: columnValuesStr(:,:)
    character(len=3000)  :: query
    type(fSQL_STATUS)    :: stat ! sqlite error status
    type(fSQL_DATABASE)  :: db   ! sqlite file handle
    type(fSQL_STATEMENT) :: stmt ! precompiled sqlite statements

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getColumnValuesDate: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_getColumnValuesDate: fSQL_open' )
    end if

    ! Get the date and time

    ! build the sqlite query
    query = "select strftime('%Y%m%d'," // trim(sqlColumnName) // &
            "), strftime('%H%M'," // trim(sqlColumnName) // ") " // &
            "from " // trim(tableName) // ";"
    write(*,*) 'odbf_getColumnValuesDate: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query), stmt, status=stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_CHAR, status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getColumnValuesDate: fSQL_get_many: ', fSQL_errmsg(stat)
      call utl_abort('odbf_getColumnValuesDate: problem with fSQL_get_many')
    end if
    write(*,*) 'odbf_getColumnValuesDate: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( columnValuesStr(numRows,2) )
    call fSQL_fill_matrix( stmt, columnValuesStr )
    allocate( columnDateValues(numRows) )
    allocate( columnTimeValues(numRows) )
    do rowIndex = 1, numRows
      read(columnValuesStr(rowIndex,1),*) columnDateValues(rowIndex)
      read(columnValuesStr(rowIndex,2),*) columnTimeValues(rowIndex)
    end do

    deallocate(columnValuesStr)

    ! close the obsDB file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine odbf_getColumnValuesDate

  !--------------------------------------------------------------------------
  ! odbf_setSurfaceType
  !--------------------------------------------------------------------------
  subroutine odbf_setSurfaceType(obsdat, headIndexBegin, fileName, tableName)
    !
    ! :Purpose: Set the surface type based on lat-lon and some external mask files.
    !
    implicit none

    ! arguments:
    type(struct_obs), intent(inout) :: obsdat
    integer,          intent(in)    :: headIndexBegin
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
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_setSurfaceType: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_setSurfaceType: fSQL_open' )
    end if

    ! build the sqlite query
    query = 'select mask_mer(lat,lon) from ' // trim(tableName) // ';'
    write(*,*) 'odbf_setSurfaceType: query ---> ', trim(query)

    ! read the values from the query result
    call fSQL_prepare( db, trim(query), stmt, status=stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_INT, status=stat )
    write(*,*) 'odbf_setSurfaceType: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( columnValues(numRows, numColumns) )
    call fSQL_fill_matrix( stmt, columnValues )

    ! set the values of STYP and TTYP
    do headTableIndex = 1, numRows
      headIndex = headTableIndex + headIndexBegin - 1
      call obs_headSet_i(obsdat, OBS_STYP, headIndex, columnValues(headTableIndex,1))
      call obs_headSet_i(obsdat, OBS_TTYP, headIndex, -1) ! Not sea ice
    end do

    ! close the obsDB file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine odbf_setSurfaceType

  !--------------------------------------------------------------------------
  ! odbf_getColumnValuesChar
  !--------------------------------------------------------------------------
  subroutine odbf_getColumnValuesChar(columnValues, fileName, tableName, &
                                      sqlColumnNames)
    !
    ! :Purpose: Read the column values from obsDB file for the specified table
    !           and column names.
    !
    implicit none

    ! arguments:
    character(len=50), allocatable, intent(out) :: columnValues(:,:)
    character(len=*),               intent(in)  :: sqlColumnNames(:)
    character(len=*),               intent(in)  :: fileName
    character(len=*),               intent(in)  :: tableName

    ! locals:
    integer :: numRows, numColumns, columnIndex
    character(len=3000)      :: query
    type(fSQL_STATUS)        :: stat ! sqlite error status
    type(fSQL_DATABASE)      :: db   ! sqlite file handle
    type(fSQL_STATEMENT)     :: stmt ! precompiled sqlite statements

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getColumnValuesChar: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_getColumnValuesChar: fSQL_open' )
    end if

    ! build the sqlite query
    query = 'select'
    numColumns = size(sqlColumnNames)
    do columnIndex = 1, numColumns
      query = trim(query) // ' ' // trim(sqlColumnNames(columnIndex))
      if (columnIndex < numColumns) query = trim(query) // ','
    end do
    query = trim(query) // ' from ' // trim(tableName) // ';'
    write(*,*) 'odbf_getColumnValuesChar: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query), stmt, status=stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_CHAR, status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getColumnValuesChar: fSQL_get_many: ', fSQL_errmsg(stat)
      call utl_abort('odbf_getColumnValuesChar: problem with fSQL_get_many')
    end if
    write(*,*) 'odbf_getColumnValuesChar: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( columnValues(numRows, numColumns) )
    call fSQL_fill_matrix( stmt, columnValues )

    ! close the obsDB file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine odbf_getColumnValuesChar

  !--------------------------------------------------------------------------
  ! odbf_getColumnValuesNum
  !--------------------------------------------------------------------------
  subroutine odbf_getColumnValuesNum(columnValues, fileName, tableName, &
                                     sqlColumnNames)
    !
    ! :Purpose: Read the column values from obsDB file for the specified table
    !           and column names.
    !
    implicit none

    ! arguments:
    real(8), allocatable, intent(out) :: columnValues(:,:)
    character(len=*),     intent(in)  :: sqlColumnNames(:)
    character(len=*),     intent(in)  :: fileName
    character(len=*),     intent(in)  :: tableName

    ! locals:
    integer :: numRows, numColumns, columnIndex
    character(len=3000)       :: query
    type(fSQL_STATUS)         :: stat ! sqlite error status
    type(fSQL_DATABASE)       :: db   ! sqlite file handle
    type(fSQL_STATEMENT)      :: stmt ! precompiled sqlite statements

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getColumnValuesNum: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_getColumnValuesNum: fSQL_open' )
    end if

    ! build the sqlite query
    query = 'select'
    numColumns = size(sqlColumnNames)
    do columnIndex = 1, numColumns
      query = trim(query) // ' ' // trim(sqlColumnNames(columnIndex))
      if (columnIndex < numColumns) query = trim(query) // ','
    end do
    query = trim(query) // ' from ' // trim(tableName) // ';'
    write(*,*) 'odbf_getColumnValuesNum: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query) , stmt, status=stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, mode=FSQL_REAL8, &
                        real8_missing=MPC_missingValue_R8, status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getColumnValuesNum: fSQL_get_many: ', fSQL_errmsg(stat)
      call utl_abort('odbf_getColumnValuesNum: problem with fSQL_get_many')
    end if
    write(*,*) 'odbf_getColumnValuesNum: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( columnValues(numRows, numColumns) )
    call fSQL_fill_matrix( stmt, columnValues )

    ! close the obsDB file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine odbf_getColumnValuesNum

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
    end do

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
      call obs_setHeadPrimaryKey( obsdat,  headIndex, headPrimaryKey(headTableIndex))
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
          end if

        end do headTable_loop

      end if

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
            end if

          end if

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
            surfEmiss = obs_bodyElem_r( obsdat, OBS_SEM, bodyIndex )
            surfEmiss = surfEmiss * 0.01D0
            call obs_bodySet_r(obsdat, OBS_SEM, bodyIndex, surfEmiss)
          end do
        end if

      end if ! obsFamily = 'TO'

    end do

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
  function odbf_midasTabColFromObsSpaceName(obsSpaceName) result(tableName)
    !
    ! :Purpose: Return the corresponding sql file column name for a
    !           given obsSpaceData column name from the midas table
    !           matching list.
    !
    implicit none

    ! arguments:
    character(len=*), intent(in) :: obsSpaceName
    character(len=lenSqlName)    :: tableName

    ! locals:
    integer :: matchIndex

    ! look in the midas table matching list
    matchIndex = utl_findloc(midasOutputNamesList(obsColIndex,:), trim(obsSpaceName))
    if (matchIndex > 0) then
      tableName = midasOutputNamesList(sqlColIndex,matchIndex)
      return
    end if

    ! not found, abort
    write(*,*) 'odbf_midasTabColFromObsSpaceName: requested obsSpace name = ', trim(obsSpaceName)
    call utl_abort('odbf_midasTabColFromObsSpaceName: obsSpace name not found in midasOutputNamesList')
    
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
  ! odbf_createMidasTable
  !--------------------------------------------------------------------------
  subroutine odbf_createMidasTable(fileName)
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
    character(len=*), parameter :: myName = 'odbf_createMidasTable'

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      call utl_abort( myName//': fSQL_open '//fSQL_errmsg(stat) )
    end if

    ! create the new MIDAS table
    query = 'create table ' // trim(midasTableName) // ' (' // new_line('A') // &
            '  ' // trim(midasKeySqlName) // ' integer primary key,' // new_line('A') // &
            '  ' // trim(headKeySqlName) // ' integer,' // new_line('A') // &
            '  ' // trim(bodyKeySqlName) // ' integer,' // new_line('A')
    do columnIndex = 1, numColMidasTableRequired
      obsColumnIndex = obs_columnIndexFromName(trim(midasOutputNamesList(2,columnIndex)))
      if (obs_columnDataType(obsColumnIndex) == 'real') then
        sqlDataType = 'double'
     else
        sqlDataType = 'integer'
      end if
      query = trim(query) // '  ' // trim(midasOutputNamesList(1,columnIndex)) // &
              ' ' // trim(sqlDataType)
      if (columnIndex < numColMidasTableRequired) query = trim(query) // ', '
      query = trim(query) // new_line('A')
    end do
    query = trim(query) // ');'
    write(*,*) myName//': query = ', trim(query)
    call fSQL_do_many( db, query, stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      call utl_abort( myName//': Problem with fSQL_do_many '//fSQL_errmsg(stat) )
    end if

    ! close the obsDB file
    call fSQL_close( db, stat ) 

  end subroutine odbf_createMidasTable
  
end module obsdbFiles_mod
