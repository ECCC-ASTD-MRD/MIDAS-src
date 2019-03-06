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

!--------------------------------------------------------------------------
!! MODULE sqliteFiles (prefix="sqlite" category='6. Observation input/output')
!!
!! *Purpose*: To store the filenames of the sqlite observation files and call
!!            subroutines in readSqlite to read and update sqlite files.
!!
!--------------------------------------------------------------------------
module sqliteFiles_mod
  
  use mathPhysConstants_mod
  use mpiVar_mod
  use sqliteRead_mod
  use obsSpaceData_mod
  use fSQLite
  use utilities_mod
  use codePrecision_mod
  use obsUtil_mod

  implicit none
  save
  private
  public :: sqlf_getDateStamp, sqlf_updateFile, sqlf_readFile, sqlf_thinFile, sqlf_writeSqlDiagFiles
  
  type(fSQL_DATABASE) :: db         ! type for SQLIte  file handle
  type(FSQL_STATUS)   :: statusSqlite

  contains

  subroutine sqlf_getDateStamp(datestamp, sqliteFileName)
    implicit none
    ! arguments
    integer                      :: dateStamp
    character(len=*), intent(in) :: sqliteFileName
    ! locals
    logical             :: isExistLogical 
    integer             :: ier, ktime, kdate, ivals, kdate_recv, ktime_recv
    integer             :: inrecs, mrfopc, dateSqlite , timeSqlite
    real(8)             :: delhh
    integer             :: nbrpdate, nbrphh, istampobs, inewhh, newdate
    character(len=128)  :: querySqlite
    character(len=9)    :: datetimeSqliteCharacter

    ier = mrfopc('MSGLVL','FATAL')
    ivals = 8
    kdate = MPC_missingValue_INT 
    ktime = MPC_missingValue_INT 
    inquire(file = trim(sqliteFileName), exist = isExistLogical )

    if ( isExistLogical ) then
      write(*,*)' Open File : ', trim(sqliteFileName)
      call fSQL_open( db, trim(sqliteFileName), statusSqlite )
      querySqlite = "select date from resume;"
      datetimeSqliteCharacter = sqlr_query(db, trim(querySqlite) )
      read(datetimeSqliteCharacter,*)  dateSqlite
      kdate = dateSqlite
      querySqlite = "select time from resume;"
      datetimeSqliteCharacter = sqlr_query( db, trim( querySqlite ) )
      read( datetimeSqliteCharacter, * ) timeSqlite 
      ktime = timeSqlite
      write(*,*) ' DATE and TIME in Sqlite file = ', dateSqlite, timeSqlite, ' (kdate, ktime) (',kdate,ktime,')'
      call fSQL_close( db, statusSqlite )
    end if

    ! Make sure all mpi tasks have a valid date (important for split sqlite files)
    call rpn_comm_allreduce(kdate, kdate_recv, 1, "MPI_INTEGER", "MPI_MAX", "GRID", ier)
    call rpn_comm_allreduce(ktime, ktime_recv, 1, "MPI_INTEGER", "MPI_MAX", "GRID", ier)

    kdate = kdate_recv
    ktime = ktime_recv
    ier = newdate(istampobs, kdate, ktime, 3)
    delhh = 3.0d0
    call INCDATR (datestamp, istampobs, delhh)
    ier = newdate(datestamp, nbrpdate, inewhh, -3)
    nbrphh = ktime
    if (nbrphh >= 21 .or. nbrphh < 3 ) then
      nbrphh = 0
    else if(nbrphh >= 3 .and. nbrphh < 9 ) then
      nbrphh = 6
    else if(nbrphh >= 9 .and. nbrphh < 15 ) then
      nbrphh = 12
    else
      nbrphh = 18
    end if
    ier = newdate(datestamp, nbrpdate, nbrphh * 1000000, 3)
    write(*,*)' SQLITE FILES VALID DATE (YYYYMMDD) : ', nbrpdate
    write(*,*)' SQLITE FILES VALID TIME       (HH) : ', nbrphh
    write(*,*)' SQLITE FILES DATESTAMP             : ', datestamp

  end subroutine sqlf_getDateStamp


  subroutine sqlf_readFile(obsdat, fileName, familyType, fileIndex )
    implicit none
    ! arguments
    type (struct_obs), intent(inout) :: obsdat
    character(len=*)                 :: fileName
    character(len=*)                 :: familyType
    integer                          :: fileINdex
    ! locals
    integer :: bodyIndex, bodyIndexBegin, bodyIndexEnd, headerIndexBegin, headerIndexEnd, headerIndex
    logical :: obsFull
    character(len=*), parameter :: my_name = 'sqlf_readFile'
    character(len=*), parameter :: my_warning = '****** '// my_name //' WARNING: '
    character(len=*), parameter :: my_error   = '******** '// my_name //' ERROR: '
    real(obs_real)  :: missingValue

    write(*,*)' '
    write(*,*)'                '//trim(my_name)//': Starting          '
    write(*,*)' '
    missingValue = real(MPC_missingValue_R8,OBS_REAL)
    write(*,*) my_name//': FileName   : ', trim(FileName)
    write(*,*) my_name//': FamilyType : ', FamilyType

    bodyIndexBegin   = obs_numbody(obsdat) + 1
    headerIndexBegin = obs_numheader(obsdat) + 1
    call sqlr_readSqlite(obsdat, trim(familyType), trim(fileName), fileIndex )
    bodyIndexEnd   = obs_numbody(obsdat)
    headerIndexEnd = obs_numheader(obsdat)

    if ( trim(familyType) /= 'TO' ) then
      call obsu_windDirectionToUV     (obsdat, headerIndexBegin, headerIndexEnd, MPC_missingValue_R4 )
      call obsu_adjustHumGZ             (obsdat, headerIndexBegin, headerIndexEnd )
      call obsu_computeVertCoordSurfObs (obsdat, headerIndexBegin, headerIndexEnd )
    end if

    do headerIndex = headerIndexBegin, headerIndexEnd
      call obs_headSet_i(obsdat, OBS_OTP, headerIndex, fileIndex)
      call obs_headSet_i(obsdat, OBS_IDF, headerIndex, fileIndex)
      call obs_setFamily(obsdat, trim(familyType), headerIndex)
    end do

    do bodyIndex = bodyIndexBegin, bodyIndexEnd
      if ( obs_columnActive_RB(obsdat, OBS_OMA) )  call obs_bodySet_r(obsdat, OBS_OMA , bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_OMP) )  call obs_bodySet_r(obsdat, OBS_OMP , bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_OER) )  call obs_bodySet_r(obsdat, OBS_OER , bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_HPHT) ) call obs_bodySet_r(obsdat, OBS_HPHT, bodyIndex, missingValue )
      if ( obs_columnActive_RB(obsdat, OBS_WORK) ) call obs_bodySet_r(obsdat, OBS_WORK, bodyIndex, missingValue )
    end do

    ! For GP family, initialize OBS_OER to element 15032 (ZTD formal error) 
    ! for all ZTD data (element 15031)
    if ( trim(familyType) == 'GP') then
      write(*,*)' Initializing OBS_OER for GB-GPS ZTD to formal error (ele 15032)'
      call obsu_setGbgpsError(obsdat, headerIndexBegin, headerIndexEnd )
    end if

    write(*,*) my_name//': obs_numheader(obsdat)', trim(familyType), obs_numheader(obsdat)
    write(*,*) my_name//': obs_numbody(obsdat)  ', trim(familyType), obs_numbody  (obsdat)
    write(*,*)' '
    write(*,*)'      '//trim(my_name)//'     END                   '
    write(*,*)' '

  end subroutine  sqlf_readFile

  
  subroutine sqlf_updateFile(obsSpaceData, fileName, familyType, fileIndex )
    implicit none
    ! arguments
    type (struct_obs), intent(inout) :: obsSpaceData
    character(len=*)                 :: fileName
    character(len=*)                 :: familyType
    integer                          :: fileIndex
    ! locals
    integer :: headerIndex
    character(len=*), parameter :: myName = 'sqlf_updateFile'
    character(len=*), parameter :: myWarning = '****** '// myName //' WARNING: '
    character(len=*), parameter :: myError   = '******** '// myName //' ERROR: '

    call tmg_start(97,'POST_UPDATESQL')
    write(*,*) myName//' Starting'
    write(*,*) myName//': FileName   : ',trim(fileName)
    write(*,*) myName//': FamilyType : ',FamilyType

    call fSQL_open( db, fileName, statusSqlite )
    if ( fSQL_error(statusSqlite) /= FSQL_OK ) then
      write(*,*) 'fSQL_open: ', fSQL_errmsg(statusSqlite )
      write(*,*) myError, fSQL_errmsg(statusSqlite )
    end if

    call sqlr_updateSqlite(db, obsSpaceData, familyType, fileName, fileIndex )
    call sqlr_insertSqlite(db, obsSpaceData, familyType, fileName, fileIndex )

    write(*,*)'  closed database -->', trim(FileName)
    call fSQL_close( db, statusSqlite )
    write(*,*)' '
    write(*,*)'================================================='
    write(*,*)'                '//trim(myName)//'    END               '
    write(*,*)'================================================='
    write(*,*)' '
    call tmg_stop(97)

  end subroutine sqlf_updateFile


  !--------------------------------------------------------------------------
  !!
  !! *Purpose*: to reduce the number of observation data in an SQL file
  !!
  !--------------------------------------------------------------------------
  subroutine sqlf_thinFile(obsSpaceData, fileName, familyType, fileIndex)
    implicit none
    ! arguments
    type (struct_obs), intent(inout) :: obsSpaceData
    character(len=*),  intent(in) :: fileName
    character(len=*),  intent(in) :: familyType
    integer,           intent(in) :: fileIndex

    ! locals
    character(len=*), parameter :: myName = 'sqlf_thinFile'

    call tmg_start(96,'POST_THINSQL')
    write(*,*) myName//' Starting'
    write(*,*) myName//': FileName   : ',trim(fileName)
    write(*,*) myName//': FamilyType : ',FamilyType

    call sqlr_thinSqlite(db, obsSpaceData, familyType, fileName, fileIndex )

    write(*,*)' '
    write(*,*)'================================================='
    write(*,*)'                '//trim(myName)//'    END               '
    write(*,*)'================================================='
    write(*,*)' '
    call tmg_stop(96)
  end subroutine sqlf_thinFile


  subroutine sqlf_writeSqlDiagFiles( obsSpaceData )
    implicit none
    ! arguments
    type (struct_obs), intent(inout) :: obsSpaceData

    call tmg_start(99,'sqlf_writeSqlDiagFiles: TOTAL: ')
    
    call sqlr_writeAllSqlDiagFiles( obsSpaceData )
    
    call tmg_stop(99)

  end subroutine sqlf_writeSqlDiagFiles


end module sqliteFiles_mod
