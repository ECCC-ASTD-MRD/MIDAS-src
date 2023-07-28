
module sqliteFiles_mod
  ! MODULE sqliteFiles_mod (prefix='sqlf' category='3. Observation input/output')
  !
  !:Purpose:  To store the filenames of the sqlite observation files and call
  !           subroutines in readSqlite to read and update sqlite files.
  !
  use mathPhysConstants_mod
  use sqliteRead_mod
  use obsSpaceData_mod
  use fSQLite
  use utilities_mod
  use codePrecision_mod
  use obsUtil_mod
  use obsVariableTransforms_mod
  use timeCoord_mod

  implicit none
  save
  private
  public :: sqlf_getDateStamp, sqlf_updateFile, sqlf_readFile, sqlf_cleanFile
  public :: sqlf_addCloudParametersandEmissivity
  
  type(fSQL_DATABASE) :: db         ! type for SQLIte  file handle
  type(FSQL_STATUS)   :: statusSqlite

  contains

  !--------------------------------------------------------------------------
  ! sqlf_getDateStamp
  !--------------------------------------------------------------------------
  subroutine sqlf_getDateStamp(dateStamp, sqliteFileName)
    !
    ! Purpose: get dateStamp from an SQLite file
    !
  
    implicit none
    
    ! Arguments:
    integer         , intent(out) :: dateStamp
    character(len=*), intent(in)  :: sqliteFileName
    
    ! Locals:
    logical              :: fileExists 
    integer              :: ier, imode, validTime, validDate, validDateRecv, validTimeRecv
    integer              :: newdate
    integer, allocatable :: headDateValues(:), headTimeValues(:)

    validDate = MPC_missingValue_INT 
    validTime = MPC_missingValue_INT 
    
    inquire(file = trim(sqliteFileName), exist = fileExists)

    if (fileExists) then
      call sqlr_getColumnValuesDate(headDateValues, headTimeValues, fileName=trim(sqliteFileName))
      call tim_getValidDateTimeFromList(headDateValues, headTimeValues, validDate, validTime)
    end if

    ! Make sure all mpi tasks have a valid date (important for split sqlite files)
    call rpn_comm_allreduce(validDate, validDateRecv, 1, "MPI_INTEGER", "MPI_MAX", "GRID", ier)
    call rpn_comm_allreduce(validTime, validTimeRecv, 1, "MPI_INTEGER", "MPI_MAX", "GRID", ier)
    
    if (validDateRecv == MPC_missingValue_INT .or. validTimeRecv == MPC_missingValue_INT) then
      write(*,*) 'sqlf_getDateStamp: WARNING: Error in getting valid date and time!'
      dateStamp = 0
    else    
      ! printable to stamp, validTime must be multiplied with 1e6 to make newdate work
      imode = 3
      ier = newdate(dateStamp, validDateRecv, validTimeRecv * 1000000, imode)
      write(*,*)'sqlf_getDateStamp: SQLite files valid date (YYYYMMDD): ', validDateRecv
      write(*,*)'sqlf_getDateStamp: SQLite files valid time       (HH): ', validTimeRecv
      write(*,*)'sqlf_getDateStamp: SQLite files dateStamp            : ', datestamp
    end if

  end subroutine sqlf_getDateStamp

  !--------------------------------------------------------------------------
  ! sqlf_readFile
  !--------------------------------------------------------------------------
  subroutine sqlf_readFile(obsdat, fileName, familyType, fileIndex)
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsdat
    character(len=*),  intent(in)    :: fileName
    character(len=*),  intent(in)    :: familyType
    integer         ,  intent(in)    :: fileIndex

    ! Locals:
    integer :: bodyIndex, bodyIndexBegin, bodyIndexEnd, headerIndexBegin, headerIndexEnd, headerIndex
    integer :: numBody, numHeader
    character(len=*), parameter :: my_name = 'sqlf_readFile'
    character(len=*), parameter :: my_warning = '****** '// my_name //' WARNING: '
    character(len=*), parameter :: my_error   = '******** '// my_name //' ERROR: '
    real(pre_obsReal)           :: missingValue

    write(*,*)' '
    write(*,*)'                '//trim(my_name)//': Starting          '
    write(*,*)' '
    missingValue = real(MPC_missingValue_R8,pre_obsReal)
    write(*,*) my_name//': FileName   : ', trim(FileName)
    write(*,*) my_name//': FamilyType : ', FamilyType

    bodyIndexBegin   = obs_numbody(obsdat) + 1
    headerIndexBegin = obs_numheader(obsdat) + 1
    call sqlr_readSqlite(obsdat, trim(familyType), trim(fileName))
    bodyIndexEnd   = obs_numbody(obsdat)
    headerIndexEnd = obs_numheader(obsdat)
    if (trim(familyType) == 'TO') then
      call sqlr_readSqlite_avhrr(obsdat, trim(fileName), headerIndexBegin, headerIndexEnd)
    end if

    if (trim(familyType) /= 'TO') then
      call ovt_transformObsValues      (obsdat, headerIndexBegin, headerIndexEnd)
      call ovt_adjustHumGZ             (obsdat, headerIndexBegin, headerIndexEnd)
      call obsu_computeVertCoordSurfObs(obsdat, headerIndexBegin, headerIndexEnd)
    end if

    do headerIndex = headerIndexBegin, headerIndexEnd
      call obs_headSet_i(obsdat, OBS_IDF, headerIndex, fileIndex)
      call obs_setFamily(obsdat, trim(familyType), headerIndex)
    end do

    ! For GP family, initialize OBS_OER to element 15032 (ZTD formal error) 
    ! for all ZTD data (element 15031)
    if (trim(familyType) == 'GP') then
      write(*,*)' Initializing OBS_OER for GB-GPS ZTD to formal error (ele 15032)'
      call obsu_setGbgpsError(obsdat, headerIndexBegin, headerIndexEnd)
    end if

    numHeader = obs_numHeader(obsdat)
    numBody   = obs_numBody(obsdat)
    write(*,*) my_name//': obs_numheader', trim(familyType), numHeader
    write(*,*) my_name//': obs_numbody  ', trim(familyType), numBody
    write(*,*)' '
    write(*,*)'      '//trim(my_name)//'     END                   '
    write(*,*)' '

  end subroutine  sqlf_readFile

  !--------------------------------------------------------------------------
  ! sqlf_updateFile
  !--------------------------------------------------------------------------
  subroutine sqlf_updateFile(obsSpaceData, fileName, familyType, fileIndex)
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsSpaceData
    character(len=*),  intent(in)    :: fileName
    character(len=*),  intent(in)    :: familyType
    integer         ,  intent(in)    :: fileIndex

    ! Locals:
    character(len=*), parameter :: myName = 'sqlf_updateFile'
    character(len=*), parameter :: myWarning = '****** '// myName //' WARNING: '
    character(len=*), parameter :: myError   = '******** '// myName //' ERROR: '

    call utl_tmg_start(13,'----UpdateSqliteFile')
    write(*,*) myName//' Starting'
    write(*,*) myName//': FileName   : ',trim(fileName)
    write(*,*) myName//': FamilyType : ',FamilyType

    call fSQL_open(db, fileName, statusSqlite)
    if (fSQL_error(statusSqlite) /= FSQL_OK) then
      write(*,*) 'fSQL_open: ', fSQL_errmsg(statusSqlite)
      write(*,*) myError, fSQL_errmsg(statusSqlite)
    end if

    call sqlr_updateSqlite(db, obsSpaceData, familyType, fileName, fileIndex)
    call sqlr_insertSqlite(db, obsSpaceData, familyType, fileName, fileIndex)

    write(*,*)'  closed database -->', trim(FileName)
    call fSQL_close(db, statusSqlite)
    write(*,*)' '
    write(*,*)'================================================='
    write(*,*)'                '//trim(myName)//'    END               '
    write(*,*)'================================================='
    write(*,*)' '
    call utl_tmg_stop(13)

  end subroutine sqlf_updateFile

  !--------------------------------------------------------------------------
  ! sqlf_cleanFile
  !--------------------------------------------------------------------------
  subroutine sqlf_cleanFile(fileName, familyType)
    !
    ! :Purpose: to reduce the number of observation data in an SQLite file
    !
    implicit none

    ! Arguments:
    character(len=*),  intent(in) :: fileName
    character(len=*),  intent(in) :: familyType

    ! Locals:
    character(len=*), parameter :: myName = 'sqlf_cleanFile'

    write(*,*) myName//': Starting'
    write(*,*) myName//': FileName   : ',trim(fileName)
    write(*,*) myName//': FamilyType : ',FamilyType

    call sqlr_cleanSqlite(db, fileName)

    write(*,*)myName//': Finished'

  end subroutine sqlf_cleanFile

  !--------------------------------------------------------------------------
  ! sqlf_addCloudParametersandEmissivity
  !--------------------------------------------------------------------------
  subroutine sqlf_addCloudParametersandEmissivity(obsSpaceData, fileIndex, fileName)
    !
    ! :Purpose: To insert cloud parameters in obsspace data into sqlite file
    !
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsSpaceData
    character(len=*),  intent(in)    :: fileName
    integer,           intent(in)    :: fileIndex

    ! Locals:
    character(len=*), parameter :: myName  = 'sqlf_addCloudParametersandEmissivity'
    character(len=*), parameter :: myError = '******** '// myName //' ERROR: '

    call fSQL_open(db, fileName, statusSqlite)
    if (fSQL_error(statusSqlite) /= FSQL_OK) then
      write(*,*) 'fSQL_open: ', fSQL_errmsg(statusSqlite)
      write(*,*) myError, fSQL_errmsg(statusSqlite)
    end if
    call sqlr_addCloudParametersandEmissivity(db, obsSpaceData,fileIndex)
    call fSQL_close(db, statusSqlite)
  end subroutine sqlf_addCloudParametersandEmissivity

end module sqliteFiles_mod
