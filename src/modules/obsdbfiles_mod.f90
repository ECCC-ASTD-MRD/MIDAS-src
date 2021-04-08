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

  implicit none
  save
  private

  ! Public subroutines and functions:
  public :: odbf_setup, odbf_isActive, odbf_readFile, odbf_updateFile

 
  ! Arrays used to match obsDB column names with obsSpaceData column names

  integer, parameter :: lenSqlName = 60
  integer, parameter :: sqlColIndex = 1
  integer, parameter :: obsColIndex = 2

  character(len=lenSqlName) :: headTableName = 'header'
  character(len=lenSqlName) :: bodyTableName = 'data'

  ! ...for the header table
  integer, parameter :: numHeadMatch = 18
  character(len=lenSqlName) :: headKeySqlName = 'ID_OBS'
  character(len=lenSqlName) :: headMatchList(2,numHeadMatch) = (/ &
       'ID_STN',        'STNID', &
       'LAT',           'LAT ', &
       'LON',           'LON ', &
       'CODTYP',        'ITY ', &
       'DATE',          'DAT ', &
       'TIME',          'ETM ', &
       'STATUS',        'ST1 ', &
       'ELEV',          'ALT ', &
       'ID_SAT',        'SAT ', &
       'INSTRUMENT',    'INS ', &
       'LAND_SEA',      'STYP', &
       'ZENITH',        'SZA ', &
       'SOLAR_ZENITH',  'SUN ', &
       'AZIMUTH',       'AZA ', &
       'TERRAIN_TYPE',  'TTYP', &
       'SENSOR',        'SEN ', &
       'SOLAR_AZIMUTH', 'SAZ ', &
       'CLOUD_COVER',   'CLF ' /)

  ! ...for the body table
  integer, parameter :: numBodyMatch = 10
  character(len=lenSqlName) :: bodyKeySqlName = 'ID_DATA'
  character(len=lenSqlName) :: bodyMatchList(2,numBodyMatch) = (/ &
       'VCOORD',        'PPP ', &
       'VARNO',         'VNM ', &
       'OBSVALUE',      'VAR ', &
       'FLAG',          'FLG ', &
       'OMP',           'OMP ', &
       'OMA',           'OMA ', &
       'FSO',           'FSO ', &
       'FG_ERROR',      'HPHT', &
       'OBS_ERROR',     'OER ', &
       'SURF_EMISS',    'SEM ' /)

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
    integer           :: nulnam, ierr
    integer, external :: fnom, fclos
    logical, save     :: nmlAlreadyRead = .false.

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
    real(8),                   allocatable :: headValues(:,:), bodyValues(:,:)
    integer(8),                allocatable :: headPrimaryKey(:)
    integer(8),                allocatable :: bodyPrimaryKey(:), bodyHeadKey(:)

    write(*,*)
    write(*,*) 'odbf_readFile: Starting'
    write(*,*)
    write(*,*) 'odbf_readFile: FileName   : ', trim(FileName)
    write(*,*) 'odbf_readFile: FamilyType : ', FamilyType

    !- 0.0 Some initialization
    call ovt_setup(elemIdList(1:numElemIdList))

    !- 1.0 Determine names of columns present in obsDB file

    call odbf_getSqlColumnNames(headCharSqlNames, fileName=trim(fileName), &
                                tableName=headTableName, dataType='varchar')
    call odbf_getSqlColumnNames(headSqlNames, fileName=trim(fileName), &
                                tableName=headTableName, dataType='numeric' )
    call odbf_getSqlColumnNames(bodySqlNames, fileName=trim(fileName), &
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

    !- 1.1 Read the contents of the file into local tables

    call odbf_getPrimaryKeys(headPrimaryKey, bodyPrimaryKey, bodyHeadKey, fileName=trim(fileName))

    call odbf_getColumnValues_char(headCharValues, fileName=trim(fileName), &
                                   tableName=headTableName, sqlColumnNames=headCharSqlNames)
    call odbf_getColumnValues_num(headValues, fileName=trim(fileName), &
                                  tableName=headTableName, sqlColumnNames=headSqlNames)
    call odbf_getColumnValues_num(bodyValues, fileName=trim(fileName), &
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

    !- 1.3 Set some other quantities in obsSpaceData

    do headIndex = headIndexBegin, headIndexEnd
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
    end do

    ! For GP family, initialize OBS_OER to element 15032 (ZTD formal error) 
    ! for all ZTD observations (element 15031)
    if ( trim(familyType) == 'GP') then
      write(*,*)'odbf_readFile: Initializing OBS_OER for GB-GPS ZTD to formal error (ele 15032)'
      call obsu_setGbgpsError(obsdat, headIndexBegin, headIndexEnd )
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
  ! odbf_updateFile
  !--------------------------------------------------------------------------
  subroutine odbf_updateFile(obsdat, fileName, familyType, fileIndex)
    !
    ! :Purpose: Update the selected columns in an obsDB file using
    !           values from obsSpaceData
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
    integer(8)           :: obsIdd, obsIdo
    integer              :: obsIdf, obsStatus
    integer              :: columnIndex, updateItemIndex, matchIndex
    integer              :: headIndex, bodyIndex, bodyIndexBegin, bodyIndexEnd
    integer              :: updateValue_i, updateList(20), fnom, fclos, nulnam, ierr
    real(8)              :: updateValue_r, obsValue
    logical              :: headFlagPresent
    character(len=4)     :: obsSpaceColumnName
    character(len=lenSqlName) :: sqlColumnName, headFlagSqlName
    character(len=3000)  :: query
    character(len=lenSqlName), allocatable :: headSqlNames(:)
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
      if ( mpi_myid == 0 ) write(*, nml=namObsDbUpdate)
      ierr = fclos(nulnam)

      ! add the observation flag (OBS_FLG) to the list of items being updated
      numberUpdateItems = numberUpdateItems + 1
      updateItemList(numberUpdateItems) = 'FLG'

    end if

    ! later, we need to know if the header-level flag (OBS_ST1) is in sql file
    call odbf_getSqlColumnNames(headSqlNames, fileName=trim(fileName), &
                                tableName=headTableName, dataType='numeric' )
    headFlagPresent = .false.
    do columnIndex = 1, size(headSqlNames)
      matchIndex = utl_findloc(headMatchList(sqlColIndex,:), headSqlNames(columnIndex))

      ! if this column name is not in the match list, skip it
      if (matchIndex == 0) cycle

      ! if this column corresponds with the header-level flag (OBS_ST1), save it
      if (trim(headMatchList(obsColIndex,matchIndex)) == 'ST1') then
        headFlagPresent = .true.
        headFlagSqlName = trim(headSqlNames(columnIndex))
      end if
    end do
    
    ! open the obsDB file
    call fSQL_open( db, trim(fileName), stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_updateFile: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_updateFile: fSQL_open' )
    end if

    ! Create the update sql query
    query = 'update ' // trim(bodyTableName) // ' set'
    do updateItemIndex = 1, numberUpdateItems
      obsSpaceColumnName = updateItemList(updateItemIndex)
      ierr = clib_toUpper(obsSpaceColumnName)

      ! get obsSpaceData column index from the name
      updateList(updateItemIndex) = obs_columnIndexFromName(trim(obsSpaceColumnName))

      ! get the sql column name from the obsSpaceData column index
      matchIndex = utl_findloc(bodyMatchList(obsColIndex,:), obsSpaceColumnName)
      if (matchIndex == 0) then
        call utl_abort( 'odbf_updateFile: Invalid obsSpaceData column name ' // &
                        'for file update: ' // trim(obsSpaceColumnName) )
      end if
      sqlColumnName = bodyMatchList(sqlColIndex,matchIndex)

      write(*,*) 'odbf_updateFile: Updating ', updateItemIndex, &
                 ', obsSpace name = ', trim(obsSpaceColumnName), &
                 ', sql column name = ', trim(sqlColumnName)
      query = trim(query) // ' ' // trim(sqlColumnName) // ' = ?'
      if (updateItemIndex < numberUpdateItems) then
        query = trim(query) // ', '
      end if
    end do

    query = trim(query)//' where ' // trim(bodyKeySqlName) // ' = ?  ;'
    write(*,*) 'odbf_updateFile: query ---> ', trim(query)
    call fSQL_prepare( db, query , stmt, stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_updateFile: fSQL_prepare: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_updateFile: fSQL_prepare' )
    end if

    call fSQL_begin(db)
    HEADER: do headIndex = 1, obs_numHeader(obsdat)
 
      obsIdf = obs_headElem_i( obsdat,OBS_IDF, headIndex )
 
      if ( obsIdf /= fileIndex ) cycle HEADER
      bodyIndexBegin = obs_headElem_i( obsdat, OBS_RLN, headIndex )
      bodyIndexEnd = bodyIndexBegin + &
                     obs_headElem_i( obsdat, OBS_NLV, headIndex ) - 1

      BODY: do bodyIndex = bodyIndexBegin, bodyIndexEnd

        obsIdd  = obs_bodyPrimaryKey( obsdat, bodyIndex )
        call fSQL_bind_param(stmt, PARAM_INDEX=numberUpdateItems+1, INT8_VAR=obsIdd)

        ITEMS: do updateItemIndex = 1, numberUpdateItems

          obsValue = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex)
          if ( obsValue == obs_missingValue_R ) cycle BODY

          if (obs_columnDataType(updateList(updateItemIndex)) == 'real') then

            updateValue_r = obs_bodyElem_r(obsdat, updateList(updateItemIndex), bodyIndex)
            if ( updateValue_r == obs_missingValue_R ) then
              call fSQL_bind_param(stmt, PARAM_INDEX=updateItemIndex)  ! sql null values
            else
              call fSQL_bind_param(stmt, PARAM_INDEX=updateItemIndex, REAL8_VAR=updateValue_r)
            end if

          else if (obs_columnDataType(updateList(updateItemIndex)) == 'integer') then

            updateValue_i = obs_bodyElem_i(obsdat, updateList(updateItemIndex), bodyIndex)
            if ( updateValue_i == nint(MPC_missingValue_R8) ) then
              call fSQL_bind_param(stmt, PARAM_INDEX=updateItemIndex)  ! sql null values
            else
              call fSQL_bind_param(stmt, PARAM_INDEX=updateItemIndex, INT_VAR=updateValue_i)
            end if

          else
            call utl_abort('odbf_updateFile: unknown data type for obs header column')
          end if

        end do ITEMS

        call fSQL_exec_stmt(stmt)

      end do BODY

    end do HEADER

    call fSQL_finalize( stmt )

    ! Update the header-level flag, if it is present in file

    if ( headFlagPresent ) then

      query = 'update ' // trim(headTableName) // ' set ' // &
              trim(headFlagSqlName) // ' = ? where ' // &
              trim(headKeySqlName) // ' = ? ;'
      write(*,*) 'odbf_updateFile: query ---> ', trim(query)
      call fSQL_prepare( db, query , stmt, status=stat )
      if ( fSQL_error(stat) /= FSQL_OK ) then
        write(*,*) 'odbf_updateFile: fSQL_prepare: ', fSQL_errmsg(stat)
        call utl_abort( 'odbf_updateFile: fSQL_prepare' )
      end if

      HEADER2: do headIndex = 1,obs_numHeader(obsdat)

        obsIdf = obs_headElem_i(obsdat, OBS_IDF, headIndex)
        if ( obsIdf /= fileIndex ) cycle HEADER2

        obsIdo    = obs_headPrimaryKey( obsdat, headIndex )
        obsStatus = obs_headElem_i(obsdat, OBS_ST1, headIndex)
        call fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR  = obsStatus )
        call fSQL_bind_param( stmt, PARAM_INDEX = 2, INT8_VAR = obsIdo )
        call fSQL_exec_stmt ( stmt )

      end do HEADER2
    
      call fSQL_finalize( stmt )

    end if

    call fSQL_commit( db )

    ! close the obsDB file
    call fSQL_close( db, stat ) 

    write(*,*)
    write(*,*) 'odbf_updateFile: finished'
    write(*,*)

    call tmg_stop(97)

  end subroutine odbf_updateFile

  !--------------------------------------------------------------------------
  ! odbf_getSqlColumnNames
  !--------------------------------------------------------------------------
  subroutine odbf_getSqlColumnNames(sqlColumnNames, fileName, tableName, dataType)
    !
    ! :Purpose: Read the column names in the obsDB file for the specified table.
    !
    implicit none

    ! arguments:
    character(len=*), allocatable, intent(out) :: sqlColumnNames(:)
    character(len=*),              intent(in)  :: fileName
    character(len=*),              intent(in)  :: tableName
    character(len=*),              intent(in)  :: dataType

    ! locals:
    integer :: numRows, numColumns, rowIndex, ierr
    character(len=100), allocatable :: charValues(:,:)
    character(len=100)       :: dataTypeCriteria
    character(len=3000)      :: query
    type(fSQL_STATUS)        :: stat ! sqlite error status
    type(fSQL_DATABASE)      :: db   ! sqlite file handle
    type(fSQL_STATEMENT)     :: stmt ! precompiled sqlite statements

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getSqlColumnNames: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_getSqlColumnNames: fSQL_open' )
    end if

    ! read the column names
    if (trim(dataType) == 'varchar') then
      dataTypeCriteria = 'substr(type,1,7)="varchar"'
    else if (trim(dataType) == 'numeric') then
      dataTypeCriteria = 'type="real" or type="integer"'
    else
      call utl_abort('odbf_getSqlColumnNames: invalid dataType = ' // trim(dataType))
    end if
    query = 'select name from pragma_table_info("' // trim(tableName) // &
            '") where ' // trim(dataTypeCriteria) // ' ;'
    call fSQL_prepare( db, trim(query) , stmt, stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_CHAR, status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getSqlColumnNames: fSQL_get_many: ', fSQL_errmsg(stat)
      call utl_abort('odbf_getSqlColumnNames: problem with fSQL_get_many')
    end if
    allocate( charValues(numRows, numColumns) )
    call fSQL_fill_matrix( stmt, charValues )

    ! copy to output array and ensure they are upper case
    allocate( sqlColumnNames(numRows) )
    do rowIndex = 1, numRows
      sqlColumnNames(rowIndex) = charValues(rowIndex,1)
      ierr = clib_toUpper(sqlColumnNames(rowIndex))
    end do
    deallocate( charValues )

    ! clean up and close the obsDB file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine odbf_getSqlColumnNames

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
  ! odbf_getColumnValues_char
  !--------------------------------------------------------------------------
  subroutine odbf_getColumnValues_char(columnValues, fileName, tableName, &
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
      write(*,*) 'odbf_getColumnValues_char: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_getColumnValues_char: fSQL_open' )
    end if

    ! build the sqlite query
    query = 'select'
    numColumns = size(sqlColumnNames)
    do columnIndex = 1, numColumns
      query = trim(query) // ' ' // trim(sqlColumnNames(columnIndex))
      if (columnIndex < numColumns) query = trim(query) // ','
    end do
    query = trim(query) // ' from ' // trim(tableName) // ';'
    write(*,*) 'odbf_getColumnValues_char: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query), stmt, status=stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_CHAR, status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getColumnValues_char: fSQL_get_many: ', fSQL_errmsg(stat)
      call utl_abort('odbf_getColumnValues_char: problem with fSQL_get_many')
    end if
    write(*,*) 'odbf_getColumnValues_char: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( columnValues(numRows, numColumns) )
    call fSQL_fill_matrix( stmt, columnValues )

    ! close the obsDB file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine odbf_getColumnValues_char

  !--------------------------------------------------------------------------
  ! odbf_getColumnValues_num
  !--------------------------------------------------------------------------
  subroutine odbf_getColumnValues_num(columnValues, fileName, tableName, &
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
    integer :: numRows, numColumns, columnIndex, elemIdIndex
    character(len=lenSqlName) :: elemIdSqlName
    character(len=3000)       :: query
    character(len=10)         :: elemIdStr
    logical                   :: elemIdPresent
    type(fSQL_STATUS)         :: stat ! sqlite error status
    type(fSQL_DATABASE)       :: db   ! sqlite file handle
    type(fSQL_STATEMENT)      :: stmt ! precompiled sqlite statements

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getColumnValues_num: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_getColumnValues_num: fSQL_open' )
    end if

    ! figure out sql column name for element ID (OBS_VNM)
    elemIdSqlName = odbf_sqlNameFromObsSpaceName('VNM')

    ! build the sqlite query
    elemIdPresent = .false.
    query = 'select'
    numColumns = size(sqlColumnNames)
    do columnIndex = 1, numColumns
      if (trim(sqlColumnNames(columnIndex)) == trim(elemIdSqlName)) elemIdPresent = .true.
      query = trim(query) // ' ' // trim(sqlColumnNames(columnIndex))
      if (columnIndex < numColumns) query = trim(query) // ','
    end do
    query = trim(query) // ' from ' // trim(tableName)
    if (elemIdPresent) then
      write(*,*) 'odbf_getColumnValues_num: selection only these element ids:', &
                 elemIdList(1:numElemIdList)
      query = trim(query) // ' where ' // trim(elemIdSqlName) // ' in ('
      do elemIdIndex = 1, numElemIdList
        write(elemIdStr,'(i6)') elemIdList(elemIdIndex)
        query = trim(query) // elemIdStr
        if (elemIdIndex < numElemIdList) then
          query = trim(query) // ','
        end if
      end do
      query = trim(query) // ')'
    end if
    query = trim(query) // ';'
    write(*,*) 'odbf_getColumnValues_num: query ---> ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query) , stmt, status=stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, mode=FSQL_REAL8, &
                        real8_missing=MPC_missingValue_R8, status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getColumnValues_num: fSQL_get_many: ', fSQL_errmsg(stat)
      call utl_abort('odbf_getColumnValues_num: problem with fSQL_get_many')
    end if
    write(*,*) 'odbf_getColumnValues_num: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( columnValues(numRows, numColumns) )
    call fSQL_fill_matrix( stmt, columnValues )

    ! close the obsDB file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine odbf_getColumnValues_num

  !--------------------------------------------------------------------------
  ! odbf_copyToObsSpaceHeadChar
  !--------------------------------------------------------------------------
  subroutine odbf_copyToObsSpaceHeadChar(obsdat, headCharSqlNames, &
                                         headCharValues, headIndexBegin)
    !
    ! :Purpose: Copy character string values from a local table into
    !           obsSpaceData header rows.
    !
    implicit none

    ! arguments:
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in)    :: headCharSqlNames(:)
    character(len=*), intent(in)    :: headCharValues(:,:)
    integer,          intent(in)    :: headIndexBegin

    ! locals:
    character(len=lenSqlName) :: stnIdSqlName
    integer :: columnIndex, matchIndex, headTableIndex, headIndex
    integer :: numRowsHeadTable

    numRowsHeadTable = size(headCharValues,1)

    do columnIndex = 1, size(headCharSqlNames)
      matchIndex = utl_findloc(headMatchList(sqlColIndex,:), headCharSqlNames(columnIndex))

      if (matchIndex == 0) then

        write(*,*) 'odbf_copyToObsSpaceHeadChar: name not in list of known ' // &
                   'header column names = ', trim(headCharSqlNames(columnIndex))

      else
        stnIdSqlName = odbf_sqlNameFromObsSpaceName('STNID')
        if (headCharSqlNames(columnIndex) /= trim(stnIdSqlName)) then
          call utl_abort('odbf_copyToObsSpaceHeadChar: only valid char column is STNID')
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

      end if

    end do

  end subroutine odbf_copyToObsSpaceHeadChar

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
                write(*,*) 'odbf_copyToObsSpaceHead: set header real column   : ', trim(headSqlNames(columnIndex))
              end if
              call obs_headSet_r(obsdat, obsColumnIndex, &
                                 headIndex, real(headValues(headTableIndex,columnIndex),pre_obsReal))
            end if
          else if (obs_columnDataType(obsColumnIndex) == 'integer') then
            if ( obs_columnActive_IH(obsdat, obsColumnIndex) ) then
              if (headTableIndex == 1) then
                write(*,*) 'odbf_copyToObsSpaceHead: set header integer column: ', trim(headSqlNames(columnIndex))
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
    !           obsSpaceData body rows.
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
    character(len=lenSqlName) :: obsValueSqlName
    integer :: columnIndex, matchIndex, bodyTableIndex, bodyIndex, headIndex
    integer :: numRowsBodyTable
    integer :: bodyColumnIndexObsValue
    integer :: lastHeadKey, obsNlv
    integer, allocatable :: matchIndexVec(:)
    integer, allocatable :: obsColumnIndex(:)

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

    ! figure out column index for observation value (OBS_VAR)
    obsValueSqlName = odbf_sqlNameFromObsSpaceName('VAR')
    bodyColumnIndexObsValue = utl_findloc(bodySqlNames(:), obsValueSqlName)

    lastHeadKey = 0
    headIndex = headIndexBegin - 1
    bodyIndex = bodyIndexBegin - 1

    bodyIndex_loop: do bodyTableIndex = 1, numRowsBodyTable
      bodyIndex = bodyIndex + 1

      ! check if obs value is null/missing
      if ( bodyValues(bodyTableIndex,bodyColumnIndexObsValue) == MPC_missingValue_R8 ) then
        bodyIndex = bodyIndex - 1
        cycle bodyIndex_loop
      end if

      ! count number of body rows for each header row (OBS_NLV)
      if ( bodyHeadKey(bodyTableIndex) /= lastHeadKey ) then
        headIndex = headIndex + 1
        call obs_headSet_i(obsdat, OBS_NLV, headIndex, 0)
        lastHeadKey = bodyHeadKey(bodyTableIndex)
      end if
      obsNLV = obs_headElem_i(obsdat, OBS_NLV, headIndex)
      call obs_headSet_i(obsdat, OBS_NLV, headIndex, obsNLV + 1)

      ! check that the primary key for header table matches the value in the body table
      if ( obs_headPrimaryKey( obsdat, headIndex ) /= &
           bodyHeadKey(bodyTableIndex) ) then
        write(*,*) 'odbf_copyToObsSpaceBody: primary key in HEADER table = ', &
                   obs_headPrimaryKey( obsdat, headIndex )
        write(*,*) 'odbf_copyToObsSpaceBody: same key in BODY table      = ', &
                   bodyHeadKey(bodyTableIndex)
        call utl_abort('odbf_copyToObsSpaceBody: Primary key of HEADER table not equal ' // &
                       'to value in BODY table')
      end if

      ! copy body primary key to obsSpaceData
      if (bodyTableIndex == 1) then
        write(*,*) 'odbf_copyToObsSpaceBody: set body primary key'
      end if
      call obs_setBodyPrimaryKey(obsdat, bodyIndex, bodyPrimaryKey(bodyTableIndex))

      ! copy real and integer values into obsSpaceData
      columnIndex_loop: do columnIndex = 1, size(bodySqlNames)
        matchIndex = matchIndexVec(columnIndex)

        if (matchIndex == 0) then
          if (bodyTableIndex == 1) then
            write(*,*) 'odbf_copyToObsSpaceBody: unknown column name    : ', &
                       trim(bodySqlNames(columnIndex))
          end if
        else
          if (obs_columnDataType(obsColumnIndex(matchIndex)) == 'real') then
            ! real values
            if ( obs_columnActive_RB(obsdat, obsColumnIndex(matchIndex)) ) then
              if (bodyTableIndex == 1) then
                write(*,*) 'odbf_copyToObsSpaceBody: set body real column   : ', trim(bodySqlNames(columnIndex))
              end if
              call obs_bodySet_r(obsdat, obsColumnIndex(matchIndex), &
                                 bodyIndex, real(bodyValues(bodyTableIndex,columnIndex),pre_obsReal))
            end if
          else if (obs_columnDataType(obsColumnIndex(matchIndex)) == 'integer') then
            ! integer values
            if ( obs_columnActive_IB(obsdat, obsColumnIndex(matchIndex)) ) then
              if (bodyTableIndex == 1) then
                write(*,*) 'odbf_copyToObsSpaceBody: set body integer column: ', trim(bodySqlNames(columnIndex))
              end if
              call obs_bodySet_i(obsdat, obsColumnIndex(matchIndex), &
                                 bodyIndex, nint(bodyValues(bodyTableIndex,columnIndex)))
            end if
          else
            call utl_abort('odbf_copyToObsSpaceBody: unknown data type for obs body column')
          end if
        end if

      end do columnIndex_loop

    end do bodyIndex_loop

    deallocate(matchIndexVec)
    deallocate(obsColumnIndex)

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
    integer :: obsTime, instrument, obsSat, codeType, sensor
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

      ! Time

      obsTime = obs_headElem_i( obsdat, OBS_ETM, headIndex )
      obsTime = obsTime/100
      call obs_headSet_i(obsdat, OBS_ETM, headIndex, obsTime)

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

        do bodyIndex = bodyIndexStart, bodyIndexEnd
          surfEmiss = obs_bodyElem_r( obsdat, OBS_SEM, bodyIndex )
          surfEmiss = surfEmiss * 0.01D0
          call obs_bodySet_r(obsdat, OBS_SEM, bodyIndex, real(surfEmiss, pre_obsReal))
        end do
        
      end if
      
    end do

  end subroutine odbf_adjustValues

  !--------------------------------------------------------------------------
  ! odbf_adjustValues
  !--------------------------------------------------------------------------
  function odbf_sqlNameFromObsSpaceName(obsSpaceName) result(sqlName)
    !
    ! :Purpose: Return the corresponding sql file column name for a
    !           given obsSpaceData column name from the matching
    !           tables.
    !
    implicit none

    ! arguments:
    character(len=*), intent(in) :: obsSpaceName
    character(len=lenSqlName)    :: sqlName

    ! locals:
    integer                   :: matchIndex

    ! first try the body matching list
    matchIndex = utl_findloc(bodyMatchList(obsColIndex,:), trim(obsSpaceName))
    if (matchIndex > 0) then
      sqlName = bodyMatchList(sqlColIndex,matchIndex)
      return
    end if

    ! now try the header matching list
    matchIndex = utl_findloc(headMatchList(obsColIndex,:), trim(obsSpaceName))
    if (matchIndex > 0) then
      sqlName = headMatchList(sqlColIndex,matchIndex)
      return
    end if

    ! not found in either list, abort
    write(*,*) 'odbf_sqlNameFromObsSpaceName: requested obsSpace name = ', trim(obsSpaceName)
    call utl_abort('odbf_sqlNameFromObsSpaceName: obsSpace name not found in matching list')
    
  end function odbf_sqlNameFromObsSpaceName

end module obsdbFiles_mod
