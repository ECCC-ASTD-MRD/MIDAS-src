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

 
  ! Arrays used to match obsDB column names with obsSpaceData column indexes

  ! ...for the header table
  integer, parameter :: numHeadMatch = 19
  character(len=100) :: headSqlNameMatch(numHeadMatch) = &
       (/ 'ID_STN', 'ID_OBS', 'LAT',   'LON',   'CODTYP', 'DATE',  'TIME', 'STATUS', 'ELEV', 'ID_SAT', 'INSTRUMENT', 'LAND_SEA', 'ZENITH', 'SOLAR_ZENITH', 'AZIMUTH', 'TERRAIN_TYPE', 'SENSOR', 'SOLAR_AZIMUTH', 'CLOUD_COVER' /)
  integer            :: headObsColMatch(numHeadMatch)  = &
       (/ 0,        OBS_IDO,  OBS_LAT, OBS_LON, OBS_ITY,  OBS_DAT, OBS_ETM, OBS_ST1, OBS_ALT, OBS_SAT, OBS_INS,      OBS_STYP,   OBS_SZA,  OBS_SUN,        OBS_AZA,   OBS_TTYP,       OBS_SEN,  OBS_SAZ,         OBS_CLF       /)

  ! ...for the body table
  integer, parameter :: numBodyMatch = 11
  character(len=100) :: bodySqlNameMatch(numBodyMatch) = &
       (/ 'ID_DATA', 'VCOORD', 'VARNO','OBSVALUE','FLAG',  'OMP',  'OMA',    'FSO',   'FG_ERROR', 'OBS_ERROR', 'SURF_EMISS'  /)
  integer            :: bodyObsColMatch(numBodyMatch)  = &
       (/ OBS_IDD,   OBS_PPP,  OBS_VNM, OBS_VAR,  OBS_FLG, OBS_OMP, OBS_OMA, OBS_FSO, OBS_HPHT,   OBS_OER,     OBS_SEM       /)

  ! NAMELIST variables
  logical :: obsDbActive
  integer :: numElemIdList
  integer :: elemIdList(100)

  ! Some important column names in obsDB files
  character(len=100) :: keyHeadSqlName = 'ID_OBS'
  character(len=100) :: keyDataSqlName = 'ID_DATA'
  character(len=100) :: elemIdSqlName  = 'VARNO'

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
    numElemIdList = 2
    elemIdList(:) = 0
    elemIdList(1) = BUFR_SST
    elemIdList(2) = BUFR_SOZ

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

  end subroutine odbf_setup

  !--------------------------------------------------------------------------
  ! odbf_isActive
  !--------------------------------------------------------------------------
  function odbf_isActive() result(isActive)
    !
    ! :Purpose: Tell the caller if the namelist indicates obsDB are active
    !
    implicit none

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

    ! arguments
    type (struct_obs), intent(inout) :: obsdat
    character(len=*)                 :: fileName
    character(len=*)                 :: familyType
    integer                          :: fileIndex

    ! locals
    integer :: bodyIndex, bodyIndexBegin, bodyIndexEnd, headIndexBegin
    integer :: headIndexEnd, headIndex, obsRln
    integer :: numBody, numHead, columnIndex
    integer :: headTableIndex, numRowsHeadTable, dataTableIndex, numRowsDataTable
    character(len=100), allocatable :: headCharSqlNames(:)
    character(len=100), allocatable :: headSqlNames(:),    dataSqlNames(:)
    character(len=50),  allocatable :: headCharValues(:,:)
    real(8),            allocatable :: headValues(:,:), dataValues(:,:)

    write(*,*)
    write(*,*) 'odbf_readFile: Starting'
    write(*,*)
    write(*,*) 'odbf_readFile: FileName   : ', trim(FileName)
    write(*,*) 'odbf_readFile: FamilyType : ', FamilyType

    !- 0.0 Some initialization
    call ovt_setup(elemIdList(1:numElemIdList))

    !- 1.0 Determine names of columns present in obsDB file

    call odbf_getSqlColumnNames(headCharSqlNames, fileName=trim(fileName), &
                                tableName='header', dataType='varchar')
    call odbf_getSqlColumnNames(headSqlNames, fileName=trim(fileName), &
                                tableName='header', dataType='numeric' )
    call odbf_getSqlColumnNames(dataSqlNames, fileName=trim(fileName), &
                                tableName='data', dataType='numeric')

    ! Print all of the column names to the listing
    do columnIndex = 1, size(headCharSqlNames)
      write(*,*) 'odbf_readFile: headCharSqlNames =', columnIndex, &
                 trim(headCharSqlNames(columnIndex))
    end do
    do columnIndex = 1, size(headSqlNames)
      write(*,*) 'odbf_readFile: headSqlNames =', columnIndex, &
                 trim(headSqlNames(columnIndex))
    end do
    do columnIndex = 1, size(dataSqlNames)
      write(*,*) 'odbf_readFile: dataSqlNames   =', columnIndex, &
                 trim(dataSqlNames(columnIndex))
    end do

    !- 1.1 Read the contents of the file into local tables

    call odbf_getColumnValues_char(headCharValues, fileName=trim(fileName), &
                                   tableName='header', sqlColumnNames=headCharSqlNames)
    call odbf_getColumnValues_num(headValues, fileName=trim(fileName), &
                                  tableName='header', sqlColumnNames=headSqlNames)
    call odbf_getColumnValues_num(dataValues, fileName=trim(fileName), &
                                  tableName='data', sqlColumnNames=dataSqlNames)
    numRowsDataTable = size(dataValues,1)
    numRowsHeadTable = size(headValues,1)

    ! For debugging, print first 10 rows of each local table to the listing
    do headTableIndex = 1, min(10, numRowsHeadTable)
      write(*,*) 'odbf_readFile: headCharValues = ', headCharValues(headTableIndex,:)
    end do
    do headTableIndex = 1, min(10, numRowsHeadTable)
      write(*,*) 'odbf_readFile: headValues = ', headValues(headTableIndex,:)
    end do
    do dataTableIndex = 1, min(10, numRowsDataTable)
      write(*,*) 'odbf_readFile: dataValues = ', dataValues(dataTableIndex,:)
    end do

    !- 1.2 Copy values from local tables into obsSpaceData

    ! Starting point for adding data to obsSpaceData
    bodyIndexBegin = obs_numBody(obsdat) + 1
    headIndexBegin = obs_numHeader(obsdat) + 1

    ! Header-Character values
    call odbf_copyToObsSpaceHeadChar(obsdat, headCharSqlNames, headCharValues, &
                                     headIndexBegin)

    ! Header-Numeric values
    call odbf_copyToObsSpaceHead(obsdat, headSqlNames, headValues, &
                                 headIndexBegin)

    ! Body-Numeric values
    call odbf_copyToObsSpaceBody(obsdat, dataSqlNames, dataValues, &
                                 bodyIndexBegin, headIndexBegin)

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

    !- 2.0 Additional changes to data after they are in obsSpaceData

    call odbf_adjustValues(obsdat, headIndexBegin, headIndexEnd)

    if ( trim(familyType) /= 'TO' ) then
      call ovt_transformObsValues      (obsdat, headIndexBegin, headIndexEnd )
      call ovt_adjustHumGZ             (obsdat, headIndexBegin, headIndexEnd )
      call obsu_computeVertCoordSurfObs(obsdat, headIndexBegin, headIndexEnd )
    end if

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
    ! for all ZTD data (element 15031)
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

    ! arguments
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in)    :: fileName   
    character(len=*), intent(in)    :: familyType
    integer,          intent(in)    :: fileIndex

    ! locals
    type(fSQL_STATUS)    :: stat ! sqlite error status
    type(fSQL_DATABASE)  :: db   ! sqlite file handle
    type(fSQL_STATEMENT) :: stmt ! precompiled sqlite statements
    integer              :: obsIdf, obsIdd, obsIdo, obsStatus
    integer              :: columnIndex, updateItemIndex, matchIndex
    integer              :: headIndex, bodyIndex, bodyIndexBegin, bodyIndexEnd
    integer              :: updateValue_i, updateList(20), fnom, fclos, nulnam, ierr
    real(8)              :: updateValue_r, obsValue
    logical              :: headFlagPresent
    character(len=4)     :: obsSpaceColumnName
    character(len=50)    :: sqlColumnName, headFlagSqlName
    character(len=3000)  :: query
    character(len=100), allocatable :: headSqlNames(:)
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
                                tableName='header', dataType='numeric' )
    headFlagPresent = .false.
    do columnIndex = 1, size(headSqlNames)
      matchIndex = findloc(headSqlNameMatch(:), trim(headSqlNames(columnIndex)), 1)

      ! if this column name is not in the match list, skip it
      if (matchIndex == 0) cycle

      ! if this column corresponds with the header-level flag (OBS_ST1), save it
      if (headObsColMatch(matchIndex) == OBS_ST1) then
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
    query = ' update data set'
    do updateItemIndex = 1, numberUpdateItems
      obsSpaceColumnName = updateItemList(updateItemIndex)
      ierr = clib_toUpper(obsSpaceColumnName)

      ! get obsSpaceData column index from the name
      updateList(updateItemIndex) = obs_columnIndexFromName(trim(obsSpaceColumnName))

      ! get the sql column name from the obsSpaceData column index
      matchIndex = findloc(bodyObsColMatch(:), updateList(updateItemIndex), 1)
      if (matchIndex == 0) then
        call utl_abort( 'odbf_updateFile: Invalid obsSpaceData column name ' // &
                        'for file update: ' // trim(obsSpaceColumnName) )
      end if
      sqlColumnName = bodySqlNameMatch(matchIndex)

      write(*,*) 'odbf_updateFile: Updating ', updateItemIndex, &
                 ', obsSpace name = ', trim(obsSpaceColumnName), &
                 ', sql column name = ', trim(sqlColumnName)
      query = trim(query) // ' ' // trim(sqlColumnName) // ' = ?'
      if (updateItemIndex < numberUpdateItems) then
        query = trim(query) // ', '
      end if
    end do

    query = trim(query)//' where id_data = ?  ;'
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
      obsIdo = obs_headElem_i( obsdat, OBS_IDO, headIndex )
      bodyIndexBegin = obs_headElem_i( obsdat, OBS_RLN, headIndex )
      bodyIndexEnd = bodyIndexBegin + &
                     obs_headElem_i( obsdat, OBS_NLV, headIndex ) - 1

      BODY: do bodyIndex = bodyIndexBegin, bodyIndexEnd

        obsIdd  = obs_bodyElem_i( obsdat, OBS_IDD, bodyIndex )
        call fSQL_bind_param(stmt, PARAM_INDEX=numberUpdateItems+1, INT_VAR=obsIdd)

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

      query = ' update header set ' // trim(headFlagSqlName) // ' = ? where id_obs = ? ;'
      write(*,*) 'odbf_updateFile: query ---> ', trim(query)
      call fSQL_prepare( db, query , stmt, stat)
      if ( fSQL_error(stat) /= FSQL_OK ) then
        write(*,*) 'odbf_updateFile: fSQL_prepare: ', fSQL_errmsg(stat)
        call utl_abort( 'odbf_updateFile: fSQL_prepare' )
      end if

      HEADER2: do headIndex = 1,obs_numHeader(obsdat)

        obsIdf = obs_headElem_i(obsdat, OBS_IDF, headIndex)
        if ( obsIdf /= fileIndex ) cycle HEADER2

        obsIdo    = obs_headElem_i(obsdat, OBS_IDO, headIndex)
        obsStatus = obs_headElem_i(obsdat, OBS_ST1, headIndex)
        call fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR  = obsStatus )
        call fSQL_bind_param( stmt, PARAM_INDEX = 2, INT_VAR  = obsIdo )
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

    ! arguments
    character(len=*), allocatable, intent(out) :: sqlColumnNames(:)
    character(len=*),              intent(in)  :: fileName
    character(len=*),              intent(in)  :: tableName
    character(len=*),              intent(in)  :: dataType

    ! locals
    integer :: numRows, numColumns, rowIndex, ierr
    character(len=100), allocatable :: charData(:,:)
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
    allocate( charData(numRows, numColumns) )
    call fSQL_fill_matrix( stmt, charData )

    ! copy to output array and ensure they are upper case
    allocate( sqlColumnNames(numRows) )
    do rowIndex = 1, numRows
      sqlColumnNames(rowIndex) = charData(rowIndex,1)
      ierr = clib_toUpper(sqlColumnNames(rowIndex))
    end do
    deallocate( charData )

    ! clean up and close the obsDB file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine odbf_getSqlColumnNames

  !--------------------------------------------------------------------------
  ! odbf_getColumnValues_char
  !--------------------------------------------------------------------------
  subroutine odbf_getColumnValues_char(columnData, fileName, tableName, &
                                       sqlColumnNames)
    !
    ! :Purpose: Read the column data from obsDB file for the specified table
    !           and column names.
    !
    implicit none

    ! arguments
    character(len=50), allocatable, intent(out) :: columnData(:,:)
    character(len=*),               intent(in)  :: sqlColumnNames(:)
    character(len=*),               intent(in)  :: fileName
    character(len=*),               intent(in)  :: tableName

    ! locals
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
    write(*,*) 'odbf_getColumnValues_char: query = ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query) , stmt, status=stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, &
                        mode=FSQL_CHAR, status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getColumnValues_char: fSQL_get_many: ', fSQL_errmsg(stat)
      call utl_abort('odbf_getColumnValues_char: problem with fSQL_get_many')
    end if
    write(*,*) 'odbf_getColumnValues_char: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( columnData(numRows, numColumns) )
    call fSQL_fill_matrix( stmt, columnData )

    ! close the obsDB file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine odbf_getColumnValues_char

  !--------------------------------------------------------------------------
  ! odbf_getColumnValues_num
  !--------------------------------------------------------------------------
  subroutine odbf_getColumnValues_num(columnData, fileName, tableName, &
                                      sqlColumnNames)
    !
    ! :Purpose: Read the column data from obsDB file for the specified table
    !           and column names.
    !
    implicit none

    ! arguments
    real(8), allocatable, intent(out) :: columnData(:,:)
    character(len=*),     intent(in)  :: sqlColumnNames(:)
    character(len=*),     intent(in)  :: fileName
    character(len=*),     intent(in)  :: tableName

    ! locals
    integer :: numRows, numColumns, columnIndex, elemIdIndex
    character(len=3000)      :: query
    character(len=10)        :: elemIdStr
    logical                  :: elemIdPresent
    type(fSQL_STATUS)        :: stat ! sqlite error status
    type(fSQL_DATABASE)      :: db   ! sqlite file handle
    type(fSQL_STATEMENT)     :: stmt ! precompiled sqlite statements

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getColumnValues_num: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_getColumnValues_num: fSQL_open' )
    end if

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
    write(*,*) 'odbf_getColumnValues_num: query = ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query) , stmt, status=stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, mode=FSQL_REAL8, &
                        real8_missing=MPC_missingValue_R8, status=stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getColumnValues_num: fSQL_get_many: ', fSQL_errmsg(stat)
      call utl_abort('odbf_getColumnValues_num: problem with fSQL_get_many')
    end if
    write(*,*) 'odbf_getColumnValues_num: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( columnData(numRows, numColumns) )
    call fSQL_fill_matrix( stmt, columnData )

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
    integer :: columnIndex, matchIndex, headTableIndex, headIndex
    integer :: numRowsHeadTable

    numRowsHeadTable = size(headCharValues,1)
    
    do columnIndex = 1, size(headCharSqlNames)
      matchIndex = findloc(headSqlNameMatch(:), headCharSqlNames(columnIndex), 1)

      if (matchIndex == 0) then

        write(*,*) 'odbf_copyToObsSpaceHeadChar: name not in list of known header column names = ', &
                   trim(headCharSqlNames(columnIndex))

      else

        if (headCharSqlNames(columnIndex) /= 'ID_STN') then
          call utl_abort('odbf_copyToObsSpaceHeadChar: only valid char column is ID_STN')
        end if
        do headTableIndex = 1, numRowsHeadTable
          headIndex = headTableIndex + headIndexBegin - 1
          call obs_set_c(obsdat, 'STID', &
                         headIndex, headCharValues(headTableIndex,columnIndex))
        end do

      end if

    end do

  end subroutine odbf_copyToObsSpaceHeadChar

  !--------------------------------------------------------------------------
  ! odbf_copyToObsSpaceHead
  !--------------------------------------------------------------------------
  subroutine odbf_copyToObsSpaceHead(obsdat, headSqlNames, &
                                     headValues, headIndexBegin)
    !
    ! :Purpose: Copy real and integer values from a local table into
    !           obsSpaceData header rows.
    !
    implicit none

    ! arguments:
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in)    :: headSqlNames(:)
    real(8),          intent(in)    :: headValues(:,:)
    integer,          intent(in)    :: headIndexBegin

    ! locals:
    integer :: columnIndex, matchIndex, headTableIndex, headIndex
    integer :: numRowsHeadTable

    numRowsHeadTable = size(headValues,1)

    column_loop: do columnIndex = 1, size(headSqlNames)
      matchIndex = findloc(headSqlNameMatch(:), headSqlNames(columnIndex), 1)

      if (matchIndex == 0) then

        write(*,*) 'odbf_copyToObsSpaceHead: unknown column name      : ', &
                   trim(headSqlNames(columnIndex))

      else

        headTable_loop: do headTableIndex = 1, numRowsHeadTable
          headIndex = headTableIndex + headIndexBegin - 1

          if (obs_columnDataType(headObsColMatch(matchIndex)) == 'real') then
            if ( obs_columnActive_RH(obsdat, headObsColMatch(matchIndex)) ) then
              if (headTableIndex == 1) then
                write(*,*) 'odbf_copyToObsSpaceHead: set header real column   : ', trim(headSqlNames(columnIndex))
              end if
              call obs_headSet_r(obsdat, headObsColMatch(matchIndex), &
                                 headIndex, real(headValues(headTableIndex,columnIndex),pre_obsReal))
            end if
          else if (obs_columnDataType(headObsColMatch(matchIndex)) == 'integer') then
            if ( obs_columnActive_IH(obsdat, headObsColMatch(matchIndex)) ) then
              if (headTableIndex == 1) then
                write(*,*) 'odbf_copyToObsSpaceHead: set header integer column: ', trim(headSqlNames(columnIndex))
              end if
              call obs_headSet_i(obsdat, headObsColMatch(matchIndex), &
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
  subroutine odbf_copyToObsSpaceBody(obsdat, dataSqlNames, &
                                     dataValues, bodyIndexBegin, headIndexBegin)
    !
    ! :Purpose: Copy real and integer values from a local table into
    !           obsSpaceData body rows.
    !
    implicit none

    ! arguments:
    type(struct_obs), intent(inout) :: obsdat
    character(len=*), intent(in)    :: dataSqlNames(:)
    real(8),          intent(in)    :: dataValues(:,:)
    integer,          intent(in)    :: bodyIndexBegin
    integer,          intent(in)    :: headIndexBegin

    ! locals:
    integer :: columnIndex, matchIndex, dataTableIndex, bodyIndex, headIndex
    integer :: numRowsDataTable
    integer :: dataColumnIndexObsValue, dataColumnIndexHeadKey
    integer :: lastHeadKey, obsNlv

    numRowsDataTable = size(dataValues,1)

    dataColumnIndexObsValue = findloc(dataSqlNames(:), 'OBSVALUE', 1)
    dataColumnIndexHeadKey  = findloc(dataSqlNames(:), keyHeadSqlName, 1)
    lastHeadKey = 0
    headIndex = headIndexBegin - 1
    bodyIndex = bodyIndexBegin - 1

    dataIndex_loop: do dataTableIndex = 1, numRowsDataTable
      bodyIndex = bodyIndex + 1

      ! check if obs value is null/missing
      if ( dataValues(dataTableIndex,dataColumnIndexObsValue) == MPC_missingValue_R8 ) then
        write(*,*) 'odbf_copyToObsSpaceBody: obs value missing, skip row - dataTableIndex = ', dataTableIndex
        bodyIndex = bodyIndex - 1
        cycle dataIndex_loop
      end if

      ! count number of body rows for each header row (OBS_NLV)
      if ( dataValues(dataTableIndex,dataColumnIndexHeadKey) /= lastHeadKey ) then
        headIndex = headIndex + 1
        call obs_headSet_i(obsdat, OBS_NLV, headIndex, 0)
        lastHeadKey = dataValues(dataTableIndex,dataColumnIndexHeadKey)
      end if
      obsNLV = obs_headElem_i(obsdat, OBS_NLV, headIndex)
      call obs_headSet_i(obsdat, OBS_NLV, headIndex, obsNLV + 1)

      ! check that the primary key for header table matches the value in the data table
      if (obs_headElem_i(obsdat, OBS_IDO, headIndex) /= &
          dataValues(dataTableIndex,dataColumnIndexHeadKey)) then
        write(*,*) 'odbf_copyToObsSpaceBody: primary key in HEADER table = ', &
                   obs_headElem_i(obsdat, OBS_IDO, headIndex)
        write(*,*) 'odbf_copyToObsSpaceBody: same key in DATA table      = ', &
                   dataValues(dataTableIndex,dataColumnIndexHeadKey)
        call utl_abort('odbf_copyToObsSpaceBody: Primary key of HEADER table not equal ' // &
                       'to value in DATA table')
      end if

      ! copy real and integer values into obsSpaceData
      columnIndex_loop: do columnIndex = 1, size(dataSqlNames)
        matchIndex = findloc(bodySqlNameMatch(:), dataSqlNames(columnIndex), 1)

        if (matchIndex == 0) then
          if (dataTableIndex == 1) then
            write(*,*) 'odbf_copyToObsSpaceBody: unknown column name    : ', &
                       trim(dataSqlNames(columnIndex))
          end if
        else
          if (obs_columnDataType(bodyObsColMatch(matchIndex)) == 'real') then
            ! real values
            if ( obs_columnActive_RB(obsdat, bodyObsColMatch(matchIndex)) ) then
              if (dataTableIndex == 1) then
                write(*,*) 'odbf_copyToObsSpaceBody: set body real column   : ', trim(dataSqlNames(columnIndex))
              end if
              call obs_bodySet_r(obsdat, bodyObsColMatch(matchIndex), &
                                 bodyIndex, real(dataValues(dataTableIndex,columnIndex),pre_obsReal))
            end if
          else if (obs_columnDataType(bodyObsColMatch(matchIndex)) == 'integer') then
            ! integer values
            if ( obs_columnActive_IB(obsdat, bodyObsColMatch(matchIndex)) ) then
              if (dataTableIndex == 1) then
                write(*,*) 'odbf_copyToObsSpaceBody: set body integer column: ', trim(dataSqlNames(columnIndex))
              end if
              call obs_bodySet_i(obsdat, bodyObsColMatch(matchIndex), &
                                 bodyIndex, nint(dataValues(dataTableIndex,columnIndex)))
            end if
          else
            call utl_abort('odbf_copyToObsSpaceBody: unknown data type for obs body column')
          end if
        end if

      end do columnIndex_loop

    end do dataIndex_loop

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

end module obsdbFiles_mod
