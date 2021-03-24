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
  ! :Purpose: To store the filenames of the sqlite observation files and call
  !           subroutines in readObsdb to read and update sqlite files that are in
  !           the new "obsDB" format.
  !
  use codePrecision_mod
  use mathPhysConstants_mod
  use obsSpaceData_mod
  use fSQLite
  use utilities_mod
  use clib_interfaces_mod
  use obsUtil_mod
  use obsVariableTransforms_mod

  implicit none
  save
  private
  public :: odbf_readFile

  ! Arrays used to match obsDB column names with obsSpaceData column indexes
  ! for the header table
  integer, parameter :: numHeadMatch = 7
  character(len=100) :: headSqlNameMatch(numHeadMatch) = &
       (/ 'ID_STN', 'ID_OBS', 'LAT',   'LON',   'CODTYP', 'DATE',  'TIME'  /)
  integer            :: headObsColMatch(numHeadMatch)  = &
       (/ 0,        OBS_IDO,  OBS_LAT, OBS_LON, OBS_ITY,  OBS_DAT, OBS_ETM /)

  ! for the body table
  integer, parameter :: numBodyMatch = 5
  character(len=100) :: bodySqlNameMatch(numBodyMatch) = &
       (/ 'ID_DATA', 'VCOORD', 'VARNO','OBSVALUE','FLAG'  /)
  integer            :: bodyObsColMatch(numBodyMatch)  = &
       (/ OBS_IDD,   OBS_PPP,  OBS_VNM, OBS_VAR,  OBS_FLG /)
contains

  !--------------------------------------------------------------------------
  ! odbf_readFile
  !--------------------------------------------------------------------------
  subroutine odbf_readFile(obsdat, fileName, familyType, fileIndex )
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
    integer :: headIndexEnd, headIndex
    integer :: numBody, numHead, columnIndex, matchIndex
    integer :: localHeadIndex, numLocalHead, localDataIndex, numLocalData
    real(pre_obsReal)               :: missingValue
    character(len=100), allocatable :: headCharSqlNames(:)
    character(len=100), allocatable :: headSqlNames(:),    dataSqlNames(:)
    character(len=50),  allocatable :: headCharValues(:,:)
    real(8),            allocatable :: headValues(:,:), dataValues(:,:)

    write(*,*)
    write(*,*) 'odbf_readFile: Starting'
    write(*,*)
    write(*,*) 'odbf_readFile: FileName   : ', trim(FileName)
    write(*,*) 'odbf_readFile: FamilyType : ', FamilyType

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

    !- 1.1 Read the contents of the file into load tables

    call odbf_getColumnData_char(headCharValues, fileName=trim(fileName), &
                                 tableName='header', sqlColumnNames=headCharSqlNames)
    call odbf_getColumnData_num(headValues, fileName=trim(fileName), &
                               tableName='header', sqlColumnNames=headSqlNames)
    call odbf_getColumnData_num(dataValues, fileName=trim(fileName), &
                                tableName='data', sqlColumnNames=dataSqlNames)
    numLocalData = size(dataValues,1)
    numLocalHead = size(headValues,1)

    ! Print first 10 rows of each local table to the listing
    do localHeadIndex = 1, min(10, numLocalHead)
      write(*,*) 'odbf_readFile: headCharValues = ', headCharValues(localHeadIndex,:)
    end do
    do localHeadIndex = 1, min(10, numLocalHead)
      write(*,*) 'odbf_readFile: headValues = ', headValues(localHeadIndex,:)
    end do
    do localDataIndex = 1, min(10, numLocalData)
      write(*,*) 'odbf_readFile: dataValues = ', dataValues(localDataIndex,:)
    end do

    !- 1.2 Copy values from local tables into obsSpaceData

    ! Starting point for adding data to obsSpaceData
    bodyIndexBegin = obs_numBody(obsdat) + 1
    headIndexBegin = obs_numHeader(obsdat) + 1

    ! Header-Character values
    do columnIndex = 1, size(headCharSqlNames)
      matchIndex = findloc(headSqlNameMatch(:), headCharSqlNames(columnIndex), 1)
      if (matchIndex == 0) then
        write(*,*) 'odbf_readFile: name not in list of known column names = ', &
                   trim(headCharSqlNames(columnIndex))
      else
        if (headCharSqlNames(columnIndex) /= 'ID_STN') then
          call utl_abort('odbf_readFile: only valid char column is ID_STN')
        end if
        do localHeadIndex = 1, numLocalHead
          headIndex = localHeadIndex + headIndexBegin - 1
          call obs_set_c(obsdat, 'STID', &
                         headIndex, headCharValues(localHeadIndex,columnIndex))
        end do
      end if
    end do

    ! Header-Numeric values
    do columnIndex = 1, size(headSqlNames)
      matchIndex = findloc(headSqlNameMatch(:), headSqlNames(columnIndex), 1)
      if (matchIndex == 0) then
        write(*,*) 'odbf_readFile: name not in list of known column names = ', &
                   trim(headSqlNames(columnIndex))
      else
        do localHeadIndex = 1, numLocalHead
          headIndex = localHeadIndex + headIndexBegin - 1
          if (obs_columnDataType(headObsColMatch(matchIndex)) == 'real') then
            if ( obs_columnActive_RH(obsdat, headObsColMatch(matchIndex)) ) then
              call obs_headSet_r(obsdat, headObsColMatch(matchIndex), &
                                 headIndex, headValues(localHeadIndex,columnIndex))
            end if
          else if (obs_columnDataType(headObsColMatch(matchIndex)) == 'integer') then
            if ( obs_columnActive_IH(obsdat, headObsColMatch(matchIndex)) ) then
              call obs_headSet_i(obsdat, headObsColMatch(matchIndex), &
                                 headIndex, nint(headValues(localHeadIndex,columnIndex)))
            end if
          else
            call utl_abort('odbf_readFile: unknown data type for obs header column')
          end if
        end do
      end if
    end do

    ! Body-Numeric values
    do columnIndex = 1, size(dataSqlNames)
      matchIndex = findloc(bodySqlNameMatch(:), dataSqlNames(columnIndex), 1)
      if (matchIndex == 0) then
        write(*,*) 'odbf_readFile: name not in list of known column names = ', &
                   trim(dataSqlNames(columnIndex))
      else
        do localDataIndex = 1, numLocalData
          bodyIndex = localDataIndex + bodyIndexBegin - 1
          if (obs_columnDataType(bodyObsColMatch(matchIndex)) == 'real') then
            if ( obs_columnActive_RB(obsdat, headObsColMatch(matchIndex)) ) then
              call obs_bodySet_r(obsdat, bodyObsColMatch(matchIndex), &
                                 bodyIndex, dataValues(localDataIndex,columnIndex))
            end if
          else if (obs_columnDataType(bodyObsColMatch(matchIndex)) == 'integer') then
            if ( obs_columnActive_IB(obsdat, headObsColMatch(matchIndex)) ) then
              call obs_bodySet_i(obsdat, bodyObsColMatch(matchIndex), &
                                 bodyIndex, nint(dataValues(localDataIndex,columnIndex)))
            end if
          else
            call utl_abort('odbf_readFile: unknown data type for obs body column')
          end if
        end do
      end if
    end do

    ! last rows added to obsSpaceData
    bodyIndexEnd = obs_numBody(obsdat)
    headIndexEnd = obs_numHeader(obsdat)

    !- 2.0 Additional changes to data after they are in obsSpaceData

    if ( trim(familyType) /= 'TO' ) then
      call ovt_transformObsValues      (obsdat, headIndexBegin, headIndexEnd )
      call ovt_adjustHumGZ             (obsdat, headIndexBegin, headIndexEnd )
      call obsu_computeVertCoordSurfObs(obsdat, headIndexBegin, headIndexEnd )
    end if

    do headIndex = headIndexBegin, headIndexEnd
      call obs_headSet_i(obsdat, OBS_OTP, headIndex, fileIndex)
      call obs_headSet_i(obsdat, OBS_IDF, headIndex, fileIndex)
      call obs_setFamily(obsdat, trim(familyType), headIndex)
    end do

    missingValue = real(MPC_missingValue_R8,pre_obsReal)
    do bodyIndex = bodyIndexBegin, bodyIndexEnd
      if ( obs_columnActive_RB(obsdat, OBS_OMA) )  call obs_bodySet_r(obsdat, OBS_OMA , bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_OMA0))  call obs_bodySet_r(obsdat, OBS_OMA0, bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_OMP) )  call obs_bodySet_r(obsdat, OBS_OMP , bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_OMP6))  call obs_bodySet_r(obsdat, OBS_OMP6, bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_OER) )  call obs_bodySet_r(obsdat, OBS_OER , bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_HPHT) ) call obs_bodySet_r(obsdat, OBS_HPHT, bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_HAHT))  call obs_bodySet_r(obsdat, OBS_HAHT, bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_WORK) ) call obs_bodySet_r(obsdat, OBS_WORK, bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_SIGI))  call obs_bodySet_r(obsdat, OBS_SIGI, bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_SIGO))  call obs_bodySet_r(obsdat, OBS_SIGO, bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_ZHA ))  call obs_bodySet_r(obsdat, OBS_ZHA , bodyIndex, missingValue )
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
    character(len=512)       :: query
    type(fSQL_STATUS)        :: stat ! sqlite error status
    type(fSQL_DATABASE)      :: db   ! sqlite file handle
    type(fSQL_STATEMENT)     :: stmt ! precompiled sqlite statements

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_readFile: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_readFile: fSQL_open' )
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
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, mode=FSQL_CHAR )
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
  ! odbf_getColumnData_char
  !--------------------------------------------------------------------------
  subroutine odbf_getColumnData_char(columnData, fileName, tableName, sqlColumnNames)
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
    character(len=512)       :: query
    type(fSQL_STATUS)        :: stat ! sqlite error status
    type(fSQL_DATABASE)      :: db   ! sqlite file handle
    type(fSQL_STATEMENT)     :: stmt ! precompiled sqlite statements

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getColumnData_char: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_getColumnData_char: fSQL_open' )
    end if

    ! build the sqlite query
    query = 'select'
    numColumns = size(sqlColumnNames)
    do columnIndex = 1, numColumns
      query = trim(query) // ' ' // trim(sqlColumnNames(columnIndex))
      if (columnIndex < numColumns) query = trim(query) // ','
    end do
    query = trim(query) // ' from ' // trim(tableName) // ';'
    write(*,*) 'odbf_getColumnData_char: query = ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query) , stmt, stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, mode=FSQL_CHAR )
    write(*,*) 'odbf_getColumnData_char: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( columnData(numRows, numColumns) )
    call fSQL_fill_matrix( stmt, columnData )

    ! close the obsDB file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine odbf_getColumnData_char

  !--------------------------------------------------------------------------
  ! odbf_getColumnData_num
  !--------------------------------------------------------------------------
  subroutine odbf_getColumnData_num(columnData, fileName, tableName, sqlColumnNames)
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
    integer :: numRows, numColumns, columnIndex
    character(len=512)       :: query
    type(fSQL_STATUS)        :: stat ! sqlite error status
    type(fSQL_DATABASE)      :: db   ! sqlite file handle
    type(fSQL_STATEMENT)     :: stmt ! precompiled sqlite statements

    ! open the obsDB file
    call fSQL_open( db, trim(fileName), stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) 'odbf_getColumnData_num: fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( 'odbf_getColumnData_num: fSQL_open' )
    end if

    ! build the sqlite query
    query = 'select'
    numColumns = size(sqlColumnNames)
    do columnIndex = 1, numColumns
      query = trim(query) // ' ' // trim(sqlColumnNames(columnIndex))
      if (columnIndex < numColumns) query = trim(query) // ','
    end do
    query = trim(query) // ' from ' // trim(tableName) // ';'
    write(*,*) 'odbf_getColumnData_num: query = ', trim(query)

    ! read the values from the file
    call fSQL_prepare( db, trim(query) , stmt, stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, mode=FSQL_REAL8 )
    write(*,*) 'odbf_getColumnData_num: numRows = ', numRows, ', numColumns = ', numColumns
    allocate( columnData(numRows, numColumns) )
    call fSQL_fill_matrix( stmt, columnData )

    ! close the obsDB file
    call fSQL_free_mem( stmt )
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 

  end subroutine odbf_getColumnData_num

end module obsdbFiles_mod
