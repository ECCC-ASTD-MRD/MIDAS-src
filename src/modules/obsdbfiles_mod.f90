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
  
  contains

  !--------------------------------------------------------------------------
  ! odbf_getColumnNames
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
    integer :: bodyIndex, bodyIndexBegin, bodyIndexEnd, headerIndexBegin, headerIndexEnd, headerIndex
    integer :: numBody, numHeader, columnIndex
    real(pre_obsReal)               :: missingValue
    character(len=100), allocatable :: headerColumnNames(:), dataColumnNames(:)

    write(*,*)
    write(*,*) 'odbf_readFile: Starting'
    write(*,*)
    write(*,*) 'odbf_readFile: FileName   : ', trim(FileName)
    write(*,*) 'odbf_readFile: FamilyType : ', FamilyType

    !- 1.0 Determine columns present in obsDB file
    call odbf_getColumnNames(headerColumnNames, fileName=trim(fileName), tableName='header')
    do columnIndex = 1, size(headerColumnNames)
      write(*,*) 'odbf_readFile: headerColumnNames =', columnIndex, trim(headerColumnNames(columnIndex))
    end do
    call odbf_getColumnNames(dataColumnNames, fileName=trim(fileName), tableName='data')
    do columnIndex = 1, size(headerColumnNames)
      write(*,*) 'odbf_readFile: dataColumnNames   =', columnIndex, trim(dataColumnNames(columnIndex))
    end do

    !- 1.1 Determine active columns in obsSpaceData

    !- 1.2 Figure out what columns are read from file and their destination

    ! read the contents of the file into obsSpaceData
    bodyIndexBegin   = obs_numbody(obsdat) + 1
    headerIndexBegin = obs_numheader(obsdat) + 1
    !call odbr_readSqlite(obsdat, trim(familyType), trim(fileName) )
    bodyIndexEnd   = obs_numbody(obsdat)
    headerIndexEnd = obs_numheader(obsdat)

    if ( trim(familyType) /= 'TO' ) then
      call ovt_transformObsValues      (obsdat, headerIndexBegin, headerIndexEnd )
      call ovt_adjustHumGZ             (obsdat, headerIndexBegin, headerIndexEnd )
      call obsu_computeVertCoordSurfObs(obsdat, headerIndexBegin, headerIndexEnd )
    end if

    do headerIndex = headerIndexBegin, headerIndexEnd
      call obs_headSet_i(obsdat, OBS_OTP, headerIndex, fileIndex)
      call obs_headSet_i(obsdat, OBS_IDF, headerIndex, fileIndex)
      call obs_setFamily(obsdat, trim(familyType), headerIndex)
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
      call obsu_setGbgpsError(obsdat, headerIndexBegin, headerIndexEnd )
    end if

    numHeader = obs_numHeader(obsdat)
    numBody   = obs_numBody(obsdat)
    write(*,*) 'odbf_readFile: after reading file, obs_numheader = ', numHeader
    write(*,*) 'odbf_readFile: after reading file, obs_numbody   = ', numBody
    write(*,*)
    write(*,*) 'odbf_readFile: finished'
    write(*,*)

  end subroutine odbf_readFile

  !--------------------------------------------------------------------------
  ! odbf_getColumnNames
  !--------------------------------------------------------------------------
  subroutine odbf_getColumnNames(columnNames, fileName, tableName)
    !
    ! :Purpose: Read the column names in the obsDB file for the specified table.
    !
    implicit none

    ! arguments
    character(len=*), allocatable, intent(out) :: columnNames(:)
    character(len=*),              intent(in)  :: fileName
    character(len=*),              intent(in)  :: tableName

    ! locals
    integer :: numRows, numColumns, rowIndex, ierr
    character(len=100), allocatable :: charData(:,:)
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
    query = 'select name from pragma_table_info("' // trim(tableName) // '");'
    call fSQL_prepare( db, trim(query) , stmt, stat )
    call fSQL_get_many( stmt, nrows=numRows, ncols=numColumns, mode=FSQL_CHAR )
    allocate( charData(numRows, numColumns) )
    call fSQL_fill_matrix_char( stmt, charData )

    ! copy to output array and ensure they are upper case
    allocate( columnNames(numRows) )
    do rowIndex = 1, numRows
      columnNames(rowIndex) = charData(rowIndex,1)
      ierr = clib_toUpper(columnNames(rowIndex))
    end do
    deallocate( charData )

    ! close the obsDB file
    call fSQL_close( db, stat ) 

  end subroutine odbf_getColumnNames

end module obsdbFiles_mod
