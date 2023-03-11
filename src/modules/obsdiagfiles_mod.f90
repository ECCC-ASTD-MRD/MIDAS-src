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

module obsDiagFiles_mod
  ! MODULE sqliteRead (prefix='diaf' category='3. Observation input/output')
  !
  ! :Purpose: To write the "diag" format SQLITE observation files. Data is stored in 
  !           obsSpaceData object.
  !
  use obsSpaceData_mod
  use midasMpi_mod
  use fSQLite
  use mathPhysConstants_mod
  use utilities_mod
  use ramDisk_mod
  use tovs_nl_mod
  use rttov_const, only : ninst
  use codtyp_mod
  use sqliteUtilities_mod
  use ensembleObservations_mod

  implicit none

  save

  private

  public :: diaf_writeAllSqlDiagFiles

  contains

  !--------------------------------------------------------------------------
  ! diaf_writeAllSqlDiagFiles
  !--------------------------------------------------------------------------
  subroutine diaf_writeAllSqlDiagFiles(obsdat, sfFileName, onlyAssimObs, addFSOdiag, ensObs_opt)
    !
    ! :Purpose: To prepare the writing of obsSpaceData content into SQLite format files
    !  
    implicit none

    ! Arguments:
    type(struct_obs)                :: obsdat         ! obsSpaceData object
    character(len=*)                :: sfFileName     ! fileName acronym used for surface obs file
    logical                         :: onlyAssimObs   ! only write assimilated obs
    logical                         :: addFSOdiag     ! include FSO column in body table
    type(struct_eob), optional      :: ensObs_opt     ! ensObs object
    
    ! Locals:
    integer                :: familyIndex, codeTypeIndex, fileIndex
    character(len=2)       :: obsFamilyList(50)
    integer                :: obsFamilyListSize
    integer                :: tovsAllCodeTypeListSize, tovsAllCodeTypeList(ninst)
    integer                :: tovsCodeTypeListSize, tovsCodeTypeList(10)
    integer                :: tovsFileNameListSize
    character(len=20)      :: tovsFileNameList(30)
    character(len=20)      :: fileName

    ! ensure all mpi tasks have same list of common obs family names
    call diaf_getObsFamilyListMpiGlobal(obsdat, obsFamilyListSize, obsFamilyList)

    ! get list of all possible tovs codetype values and unique list of corresponding filenames
    call tvs_getAllIdBurpTovs(tovsAllCodeTypeListSize, tovsAllCodeTypeList)
    write(*,*) 'tovsAllCodeTypeListSize = ', tovsAllCodeTypeListSize
    write(*,*) 'tovsAllCodeTypeList = ', tovsAllCodeTypeList(1:tovsAllCodeTypeListSize)
    
    tovsFileNameListSize = 0
    tovsFileNameList(:) = 'XXXXX'
    do codeTypeIndex = 1, tovsAllCodeTypeListSize
      fileName = diaf_getObsFileName('TO', codeType_opt=tovsAllCodeTypeList(codeTypeIndex))
      if (all(tovsFileNameList(:) /= fileName)) then
        tovsFileNameListSize = tovsFileNameListSize + 1
        tovsFileNameList(tovsFileNameListSize) = fileName
      end if
    end do
    write(*,*) 'tovsFileNameListSize = ', tovsFileNameListSize
    write(*,*) 'tovsFileNameList = ', tovsFileNameList(1:tovsFileNameListSize)
    
    do familyIndex = 1, obsFamilyListSize

      write(*,*) 'diaf_writeAllSqlDiagFiles: Family = ', familyIndex, obsFamilyList(familyIndex)

      if (obsFamilyList(familyIndex) == 'TO') then

        do fileIndex = 1, tovsFileNameListSize
          fileName = tovsFileNameList(fileIndex)
          write(*,*) 'tovs filename = ', fileName

          ! get list of codetypes associated with this filename
          tovsCodeTypeListSize = 0
          tovsCodeTypeList(:) = MPC_missingValue_INT
          do codeTypeIndex = 1, tovsAllCodeTypeListSize
            if (fileName == diaf_getObsFileName('TO', codeType_opt=tovsAllCodeTypeList(codeTypeIndex))) then
              tovsCodeTypeListSize = tovsCodeTypeListSize + 1
              tovsCodeTypeList(tovsCodeTypeListSize) = tovsAllCodeTypeList(codeTypeIndex)
            end if
          end do

          write(*,*) 'tovsCodeTypeListSize = ', tovsCodeTypeListSize
          write(*,*) 'tovsCodeTypeList = ', tovsCodeTypeList(1:tovsCodeTypeListSize) 
          call diaf_writeSqlDiagFile(obsdat, 'TO', onlyAssimObs, addFSOdiag, &
                                     tovsFileNameList(fileIndex), &
                                     tovsCodeTypeList(1:tovsCodeTypeListSize), & 
                                     ensObs_opt=ensObs_opt ) 
        end do

      else

        fileName = diaf_getObsFileName(obsFamilyList(familyIndex), sfFileName_opt=sfFileName)
        call diaf_writeSqlDiagFile(obsdat, obsFamilyList(familyIndex), &
                                   onlyAssimObs, addFSOdiag, fileName, & 
                                   ensObs_opt=ensObs_opt ) 

      end if   
      
    end do

  end subroutine diaf_writeAllSqlDiagFiles
 
  !--------------------------------------------------------------------------
  ! diaf_writeSqlDiagFile
  !--------------------------------------------------------------------------
  subroutine diaf_writeSqlDiagFile(obsdat, obsFamily, onlyAssimObs, addFSOdiag, instrumentFileName, codeTypeList_opt, ensObs_opt)
    !
    ! :Purpose: To write the obsSpaceData content into SQLite format files
    !
    implicit none

    ! Arguments:
    type(struct_obs)           , intent(inout) :: obsdat
    character(len=*)           , intent(in)    :: obsFamily
    logical                    , intent(in)    :: onlyAssimObs
    logical                    , intent(in)    :: addFSOdiag
    character(len=*)           , intent(in)    :: instrumentFileName
    integer          , optional, intent(in)    :: codeTypeList_opt(:)
    type(struct_eob) , optional                :: ensObs_opt 

    ! Locals:
    type(fSQL_DATABASE)    :: db                                   ! type for SQLIte  file handle
    type(fSQL_STATEMENT)   :: stmtData, stmtHeader, stmtEnsObs     ! type for precompiled SQLite statements
    type(fSQL_STATUS)      :: stat                                 ! type for error status
    integer                :: obsVarno, obsFlag, ASS, vertCoordType, codeType, date, time, idObs, idData, memberIndex
    real                   :: obsValue, OMA, OMP, OER, FGE, PPP, lon, lat, altitude, ENSOBSTRL, ENSOBSANL
    real                   :: latData, lonData
    real                   :: ensInnovStdDev, ensObsErrStdDev, zhad, fso
    integer                :: numberInsertions, numHeaders, headerIndex, bodyIndex, obsNlv, obsRln
    character(len = 512)   :: queryData, queryHeader, queryCreate, queryCreateEnsObs 
    character(len = 12)    :: idStation
    character(len=30)      :: fileNameExtention
    character(len=256)     :: fileName, fileNameDir
    character(len=4)       :: cmyidx, cmyidy
    logical                :: writeHeader
        
    ! determine initial idData,idObs to ensure unique values across mpi tasks
    call sqlu_getInitialIdObsData(obsDat, obsFamily, idObs, idData, codeTypeList_opt)

    ! return if this mpi task does not have any observations of this family
    if (trim(instrumentFileName) == 'XXXXX') return

    ! check if any obs exist for this file, if not return
    numHeaders = 0
    call obs_set_current_header_list(obsdat, obsFamily)
    HEADERCOUNT: do
      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADERCOUNT
      if (present(codeTypeList_opt)) then
        codeType  = obs_headElem_i(obsdat, OBS_ITY, headerIndex)
        if (all(codeTypeList_opt(:) /= codeType)) cycle HEADERCOUNT
      end if
      numHeaders = numHeaders + 1
    end do HEADERCOUNT
    if (numHeaders == 0) return

    fileNameDir = trim(ram_getRamDiskDir())
    if (fileNameDir == ' ') &
    write(*,*) 'diaf_writeSqlDiagFile: WARNING! The program may be slow creating many sqlite files in the same directory.'
    write(*,*) 'diaf_writeSqlDiagFile: WARNING! Please, use the ram disk option prior to MIDAS run!'

    if (obs_mpiLocal(obsdat)) then
      write(cmyidy,'(I4.4)') (mmpi_myidy + 1)
      write(cmyidx,'(I4.4)') (mmpi_myidx + 1)
      fileNameExtention  = trim(cmyidx) // '_' // trim(cmyidy)
    else
      if (mmpi_myid > 0) return
      fileNameExtention = ' '
    end if
    
    fileName = trim(fileNameDir) // 'obs/dia' // trim(instrumentFileName) // '_' // trim(fileNameExtention)

    write(*,*) 'diaf_writeSqlDiagFile: Creating file: ', trim(fileName)
    call fSQL_open(db, fileName, stat)
    if (fSQL_error(stat) /= FSQL_OK) write(*,*) 'diaf_writeSqlDiagFile: fSQL_open: ', fSQL_errmsg(stat),' filename: '//trim(fileName)

    ! Create the tables HEADER and DATA
    if (addFSOdiag) then
      queryCreate = 'create table header (id_obs integer primary key, id_stn varchar(50), lat real, lon real, &
                    &codtyp integer, date integer, time integer, elev real); &
                    &create table data (id_data integer primary key, id_obs integer, varno integer, vcoord real, &
                    &vcoord_type integer, obsvalue real, flag integer, oma real, ompt real, oma0 real, omp real, &
                    &an_error real, fg_error real, obs_error real, sigi real, sigo real, zhad real, lat real, lon real, &
                    &fso real);'
    else
      queryCreate = 'create table header (id_obs integer primary key, id_stn varchar(50), lat real, lon real, &
                    &codtyp integer, date integer, time integer, elev real); &
                    &create table data (id_data integer primary key, id_obs integer, varno integer, vcoord real, &
                    &vcoord_type integer, obsvalue real, flag integer, oma real, ompt real, oma0 real, omp real, &
                    &an_error real, fg_error real, obs_error real, sigi real, sigo real, zhad real, lat real, lon real);'
    end if

    call fSQL_do_many(db, queryCreate, stat)
    if (fSQL_error(stat) /= FSQL_OK) call sqlu_handleError(stat, 'diaf_writeSqlDiagFile: fSQL_do_many with query: '//trim(queryCreate))
    
    ! If the analysis members in obs space are allocated, make table queries for trial members and analysis members
    if ( present( ensObs_opt ) ) then
      ! Create
      queryCreateEnsObs = 'create table ensobs (id_data integer, id_obs integer, id_member integer, obstrl real, obsanl real);'
      call fSQL_do_many(db, queryCreateEnsObs, stat)
      if (fSQL_error(stat) /= FSQL_OK) call sqlu_handleError(stat, 'diaf_writeSqlDiagFile: fSQL_do_many with query: '//trim(queryCreateEnsObs))
    end if
    
    if (addFSOdiag) then
      queryData = 'insert into data (id_data, id_obs, varno, vcoord, vcoord_type, obsvalue, flag, oma, oma0, ompt, fg_error, &
                   &obs_error, sigi, sigo, zhad, lat, lon, fso) values(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);'
    else
      queryData = 'insert into data (id_data, id_obs, varno, vcoord, vcoord_type, obsvalue, flag, oma, oma0, ompt, fg_error, &
                  &obs_error, sigi, sigo, zhad, lat, lon) values(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);'
    end if
    
    queryHeader = 'insert into header (id_obs, id_stn, lat, lon, date, time, codtyp, elev) values(?,?,?,?,?,?,?,?); '

    write(*,*) 'diaf_writeSqlDiagFile: Insert query Data   = ', trim(queryData)
    write(*,*) 'diaf_writeSqlDiagFile: Insert query Header = ', trim(queryHeader)

    call fSQL_begin(db)
    call fSQL_prepare(db, queryData, stmtData, stat)
    if (fSQL_error(stat) /= FSQL_OK) call sqlu_handleError(stat, 'diaf_writeSqlDiagFile: fSQL_prepare: ')
    call fSQL_prepare(db, queryHeader, stmtHeader, stat)
    if (fSQL_error(stat) /= FSQL_OK) call sqlu_handleError(stat, 'diaf_writeSqlDiagFile: fSQL_prepare: ')
    
    if ( present( ensObs_opt ) ) then
      ! Insert
      queryCreateEnsObs = 'insert into ensobs (id_data, id_obs, id_member, obstrl, obsanl) values(?,?,?,?,?);'
      write(*,*) 'diaf_writeSqlDiagFile: Insert query EnsObs   = ', trim(queryCreateEnsObs)

      call fSQL_prepare(db, queryCreateEnsObs, stmtEnsObs, stat)
      
      if (fSQL_error(stat) /= FSQL_OK) call sqlu_handleError(stat, 'diaf_writeSqlDiagFile: fSQL_prepare: ')
    end if
    
    numberInsertions = 0
    call obs_set_current_header_list(obsdat, obsFamily)
    HEADER: do

      headerIndex = obs_getHeaderIndex(obsdat)
      if (headerIndex < 0) exit HEADER
        
      codeType  = obs_headElem_i(obsdat, OBS_ITY, headerIndex)
      if (present(codeTypeList_opt)) then
        if (all(codeTypeList_opt(:) /= codeType)) cycle HEADER
      end if

      obsRln    = obs_headElem_i(obsdat, OBS_RLN, headerIndex)
      obsNlv    = obs_headElem_i(obsdat, OBS_NLV, headerIndex)
      idStation = obs_elem_c    (obsdat, 'STID' , headerIndex) 
      altitude  = obs_headElem_r(obsdat, OBS_ALT, headerIndex)      
      lon       = obs_headElem_r(obsdat, OBS_LON, headerIndex) * MPC_DEGREES_PER_RADIAN_R8
      lat       = obs_headElem_r(obsdat, OBS_LAT, headerIndex) * MPC_DEGREES_PER_RADIAN_R8
      if (lon > 180.) lon = lon - 360.
      date      = obs_headElem_i(obsdat, OBS_DAT, headerIndex)
      time      = obs_headElem_i(obsdat, OBS_ETM, headerIndex) * 100.

      ! check if at least one observation will be written, if not skip this header
      if (onlyAssimObs) then
        writeHeader = .false.
        BODYCHECK: do bodyIndex = obsRln, obsNlv + obsRln -1
          OMP = obs_bodyElem_r(obsdat, OBS_OMP , bodyIndex)
          FGE = obs_bodyElem_r(obsdat, OBS_HPHT, bodyIndex)
          ASS = obs_bodyElem_i(obsdat, OBS_ASS , bodyIndex)
          if ((ASS == obs_assimilated) .and.     &
              (OMP /= obs_missingValue_R) .and.  &
              (FGE /= obs_missingValue_R)) writeHeader = .true.
        end do BODYCHECK
        if (.not. writeHeader) cycle HEADER
      end if

      idObs = idObs + 1
      call fSQL_bind_param(stmtHeader, param_index = 1, int_var  = idObs)
      call fSQL_bind_param(stmtHeader, param_index = 2, char_var = idStation)
      call fSQL_bind_param(stmtHeader, param_index = 3, real_var = lat) 
      call fSQL_bind_param(stmtHeader, param_index = 4, real_var = lon) 
      call fSQL_bind_param(stmtHeader, param_index = 5, int_var  = date) 
      call fSQL_bind_param(stmtHeader, param_index = 6, int_var  = time) 
      call fSQL_bind_param(stmtHeader, param_index = 7, int_var  = codeType) 
      call fSQL_bind_param(stmtHeader, param_index = 8, real_var = altitude)
      call fSQL_exec_stmt (stmtHeader)

      BODY: do bodyIndex = obsRln, obsNlv + obsRln -1
         
        obsVarno      = obs_bodyElem_i(obsdat, OBS_VNM , bodyIndex)
        obsFlag       = obs_bodyElem_i(obsdat, OBS_FLG , bodyIndex)
        vertCoordType = obs_bodyElem_i(obsdat, OBS_VCO , bodyIndex)
        obsValue      = obs_bodyElem_r(obsdat, OBS_VAR , bodyIndex)
        OMA           = obs_bodyElem_r(obsdat, OBS_OMA , bodyIndex)
        OMP           = obs_bodyElem_r(obsdat, OBS_OMP , bodyIndex)
        OER           = obs_bodyElem_r(obsdat, OBS_OER , bodyIndex)
        FGE           = obs_bodyElem_r(obsdat, OBS_HPHT, bodyIndex)
        PPP           = obs_bodyElem_r(obsdat, OBS_PPP , bodyIndex)
        ASS           = obs_bodyElem_i(obsdat, OBS_ASS , bodyIndex)
        latData       = obs_bodyElem_r(obsdat, OBS_LATD, bodyIndex)
        lonData       = obs_bodyElem_r(obsdat, OBS_LOND, bodyIndex)

        ! skip obs if it was not assimilated
        if (onlyAssimObs) then
          if ((ASS /= obs_assimilated) .or.     &
              (OMP == obs_missingValue_R) .or.  &
              (FGE == obs_missingValue_R)) cycle BODY
        end if

        if (obs_columnActive_RB(obsdat, OBS_SIGI)) then
          ensInnovStdDev = obs_bodyElem_r(obsdat, OBS_SIGI, bodyIndex)
        else
          ensInnovStdDev = obs_missingValue_R
        end if
        if (obs_columnActive_RB(obsdat, OBS_SIGO)) then
          ensObsErrStdDev = obs_bodyElem_r(obsdat, OBS_SIGO, bodyIndex)
        else
          ensObsErrStdDev = obs_missingValue_R
        end if
        if (obs_columnActive_RB(obsdat, OBS_ZHA)) then
          zhad = obs_bodyElem_r(obsdat, OBS_ZHA, bodyIndex)
        else
          zhad = obs_missingValue_R
        end if
        if (addFSOdiag) then
          if (obs_columnActive_RB(obsdat, OBS_FSO)) then
            fso = obs_bodyElem_r(obsdat, OBS_FSO, bodyIndex)
          else
            fso = obs_missingValue_R
          end if
        end if

        select case(obsFamily)
          case ('UA', 'AI', 'SW')
            if (vertCoordType == 2) vertCoordType = 7004
          case ('RO')
            vertCoordType = 7007
          case ('RA')
            vertCoordType = 7007
          case ('PR')
            vertCoordType = 7006
            PPP = PPP - altitude
          case ('TO')
            vertCoordType = 5042
            if(codeType == codtyp_get_codtyp('amsua') .or. &
               codeType == codtyp_get_codtyp('amsub') .or. &
               codeType == codtyp_get_codtyp('mhs')) vertCoordType = 2150
          case ('SF', 'SC', 'GP')
            vertCoordType = MPC_missingValue_INT 
        end select

        ! insert order: id_data, id_obs, varno, vcoord, vcoord_type, obsvalue, flag, oma, oma0, ompt, fg_error, obs_error, sigi, sigo, zhad, fso
        idData = idData + 1
        call fSQL_bind_param(stmtData, param_index = 1, int_var  = idData)
        call fSQL_bind_param(stmtData, param_index = 2, int_var  = idObs)
        call fSQL_bind_param(stmtData, param_index = 3, int_var  = obsVarno)
        call fSQL_bind_param(stmtData, param_index = 4, real_var = PPP)
        if (vertCoordType == MPC_missingValue_INT) then
          call fSQL_bind_param(stmtData, param_index = 5) 
        else
          call fSQL_bind_param(stmtData, param_index = 5, int_var  = vertCoordType) 
        end if
        call fSQL_bind_param(stmtData, param_index = 6, real_var = obsValue) 
        call fSQL_bind_param(stmtData, param_index = 7, int_var  = obsFlag)
        if (OMA == obs_missingValue_R) then
          call fSQL_bind_param(stmtData, param_index = 8) 
          call fSQL_bind_param(stmtData, param_index = 9) 
        else
          call fSQL_bind_param(stmtData, param_index = 8, real_var = OMA)
          call fSQL_bind_param(stmtData, param_index = 9, real_var = OMA)
        end if
        if (OMP == obs_missingValue_R) then
          call fSQL_bind_param(stmtData, param_index = 10) 
        else
          call fSQL_bind_param(stmtData, param_index = 10, real_var = OMP)
        end if
        if (FGE == obs_missingValue_R) then
          call fSQL_bind_param(stmtData, param_index = 11) 
        else
          call fSQL_bind_param(stmtData, param_index = 11, real_var = FGE)
        end if
        if (OER == obs_missingValue_R) then
          call fSQL_bind_param(stmtData, param_index = 12) 
        else
          call fSQL_bind_param(stmtData, param_index = 12, real_var = OER)
        end if 
        if (ensInnovStdDev == obs_missingValue_R) then
          call fSQL_bind_param(stmtData, param_index = 13) 
        else
          call fSQL_bind_param(stmtData, param_index = 13, real_var = ensInnovStdDev)
        end if 
        if (ensObsErrStdDev == obs_missingValue_R) then
          call fSQL_bind_param(stmtData, param_index = 14) 
        else
          call fSQL_bind_param(stmtData, param_index = 14, real_var = ensObsErrStdDev)
        end if 
        if (zhad == obs_missingValue_R) then
          call fSQL_bind_param(stmtData, param_index = 15) 
        else
          call fSQL_bind_param(stmtData, param_index = 15, real_var = zhad)
        end if 
        if (latData == obs_missingValue_R) then
          call fSQL_bind_param(stmtData, param_index = 16) 
        else
          latData = latData * MPC_DEGREES_PER_RADIAN_R8
          call fSQL_bind_param(stmtData, param_index = 16, real_var = latData)
        end if 
        if (lonData == obs_missingValue_R) then
          call fSQL_bind_param(stmtData, param_index = 17) 
        else
          lonData = lonData * MPC_DEGREES_PER_RADIAN_R8
          call fSQL_bind_param(stmtData, param_index = 17, real_var = lonData)
        end if 
        if (addFSOdiag) then
          if (fso == obs_missingValue_R) then
            call fSQL_bind_param(stmtData, param_index = 18)
          else
            call fSQL_bind_param(stmtData, param_index = 18, real_var = fso)
          end if
        end if

        call fSQL_exec_stmt (stmtData)

        numberInsertions = numberInsertions + 1

        if ( present( ensObs_opt ) ) then
          if ( .not. allocated(ensObs_opt%Yb_r4) ) then
            call utl_abort('diaf_writeSqlDiagFile: ensObs%Yb_r4 must be allocated and it is not')
          end if
          ! Loop over members. insert order: id_data, id_obs, id_member, obstrl, obsanl
          do memberIndex = 1, ensObs_opt%numMembers
            ENSOBSTRL = ensObs_opt%Yb_r4(memberIndex,bodyIndex)
            if (ensObs_opt%meanRemoved) then
              ENSOBSTRL = ENSOBSTRL + ensObs_opt%meanYb(bodyIndex) ! Yb_r4 has mean removed, so add back
            end if
            call fSQL_bind_param(stmtEnsObs, param_index = 1, int_var  = idData)
            call fSQL_bind_param(stmtEnsObs, param_index = 2, int_var  = idObs)
            call fSQL_bind_param(stmtEnsObs, param_index = 3, int_var  = memberIndex)
            call fSQL_bind_param(stmtEnsObs, param_index = 4, real_var = ENSOBSTRL)
            if ( allocated(ensObs_opt%Ya_r4) ) then
              ENSOBSANL = ensObs_opt%Ya_r4(memberIndex,bodyIndex)
              call fSQL_bind_param(stmtEnsObs, param_index = 5, real_var = ENSOBSANL)
            else
              call fSQL_bind_param(stmtEnsObs, param_index = 5)
            end if
            call fSQL_exec_stmt (stmtEnsObs)
          end do
        end if  
          
      end do BODY
     
    end do HEADER
    
    call fSQL_finalize (stmtData)

    write(*,*) 'diaf_writeSqlDiagFile: Observation Family: ', obsFamily,', number of insertions: ', numberInsertions

    call fSQL_commit(db)
    call fSQL_close(db, stat)

  end subroutine diaf_writeSqlDiagFile

  !--------------------------------------------------------------------------
  ! diaf_getObsFileName
  !--------------------------------------------------------------------------
  function diaf_getObsFileName(obsFamily, sfFileName_opt, codetype_opt) result(fileName)
    !
    ! :Purpose: Return the part of the observation file name associated
    !           with the type of observation it contains.
    !
    implicit none

    ! Arguments:
    character(len=*)           :: obsFamily
    character(len=*), optional :: sfFileName_opt ! fileName acronym used for surface obs file
    integer, optional          :: codetype_opt
    character(len=20)          :: fileName

    if (obsFamily == 'TO') then
      if (.not. present(codetype_opt)) then
        call utl_abort('diaf_getObsFileName: codetype_opt must be specified for TO family')
      end if

      if (codtyp_get_name(codeType_opt) == 'radianceclear') then
        fileName  = 'csr'
      else if (codtyp_get_name(codeType_opt) == 'mhs' .or. codtyp_get_name(codeType_opt) == 'amsub') then
        if (tvs_isInstrumAllskyHuAssim(tvs_getInstrumentId(codtyp_get_name(codeType_opt)))) then
          fileName = 'to_amsub_allsky'
        else
          fileName = 'to_amsub'
        end if
      else if (codtyp_get_name(codeType_opt) == 'amsua') then
        if (tvs_isInstrumAllskyTtAssim(tvs_getInstrumentId(codtyp_get_name(codeType_opt))) .or. &
            tvs_isInstrumAllskyTtHuAssim(tvs_getInstrumentId(codtyp_get_name(codeType_opt)))) then
          fileName = 'to_amsua_allsky'
        else
          fileName = 'to_amsua'
        end if
      else if (codtyp_get_name(codeType_opt) == 'ssmi') then
        fileName = 'ssmis'
      else if (codtyp_get_name(codeType_opt) == 'crisfsr') then
        fileName = 'cris'
      else if (codtyp_get_name(codeType_opt) == 'atms') then
        if (tvs_isInstrumAllskyTtAssim(tvs_getInstrumentId(codtyp_get_name(codeType_opt))) .or. &
            tvs_isInstrumAllskyHuAssim(tvs_getInstrumentId(codtyp_get_name(codeType_opt)))) then
          fileName = 'atms_allsky'
        else
          fileName = 'atms'        
        end if
      else
        fileName = codtyp_get_name(codeType_opt)
      end if
    else
      if (.not. present(sfFileName_opt)) then
        call utl_abort('diaf_getObsFileName: sfFileName_opt must be specified')
      end if
      call up2low(obsFamily, fileName)
      if (fileName == 'ra') fileName = 'radar'
      if (fileName == 'sf') then
        ! use either 'sf' or 'sfc' for filename with surface obs
        fileName = sfFileName_opt
      end if
    end if

  end function diaf_getObsFileName

  !--------------------------------------------------------------------------
  ! diaf_getObsFamilyListMpiGlobal
  !--------------------------------------------------------------------------
  subroutine diaf_getObsFamilyListMpiGlobal(obsdat, obsFamilyListSizeCommon,  &
                                            obsFamilyListCommon)
    !
    ! :Purpose: Obtain a common set of obs family names over all mpi tasks
    !
    implicit none
      
    ! Arguments:
    type(struct_obs) :: obsdat
    integer          :: obsFamilyListSizeCommon
    character(len=*) :: obsFamilyListCommon(:)

    ! Locals:
    integer                       :: headerIndex, familyIndex, charIndex, procIndex, nsize, ierr
    integer                       :: obsFamilyListSizeMpiLocal, obsFamilyListSizeMaxMpiLocal, obsFamilyListSizeMax
    character(len=2), allocatable :: obsFamilyListMpiLocal(:)
    character(len=2), allocatable :: obsFamilyListMpiGlobal(:,:)
    character(len=2)              :: currentObsFamily
    integer, allocatable          :: intObsFamilyListMpiLocal(:,:)
    integer, allocatable          :: intObsFamilyListMpiGlobal(:,:,:)
    integer, allocatable          :: allObsFamilyListSizeMpiLocal(:)

    obsFamilyListSizeMax = size(obsFamilyListCommon)
    write(*,*) 'obsFamilyListSizeMax =', obsFamilyListSizeMax

    ! get family list for this mpi task
    obsFamilyListSizeMpiLocal = 0
    allocate(obsFamilyListMpiLocal(obsFamilyListSizeMax))
    obsFamilyListMpiLocal(:) = 'XX'
    HEADER: do headerIndex = 1, obs_numHeader(obsdat)
      currentObsFamily = obs_getFamily(obsdat, headerIndex) 
      if (any(obsFamilyListMpiLocal(:) == currentObsFamily)) cycle HEADER
      obsFamilyListSizeMpiLocal = obsFamilyListSizeMpiLocal + 1
      obsFamilyListMpiLocal(obsFamilyListSizeMpiLocal) = currentObsFamily
      write(*,*) 'add the family: ', currentObsFamily
    end do HEADER
    write(*,*) 'obsFamilyListSizeMpiLocal =', obsFamilyListSizeMpiLocal
    write(*,*) 'obsFamilyListMpiLocal = ', obsFamilyListMpiLocal(1:obsFamilyListSizeMpiLocal)

    allocate(allObsFamilyListSizeMpiLocal(mmpi_nprocs))
    call rpn_comm_allgather(obsFamilyListSizeMpiLocal,    1, 'mpi_integer',  &
                            allObsFamilyListSizeMpiLocal, 1, 'mpi_integer', 'GRID', ierr)
    call rpn_comm_allreduce(obsFamilyListSizeMpiLocal, obsFamilyListSizeMaxMpiLocal,1,'mpi_integer','mpi_max','GRID',ierr)

    ! convert local family list from characters to integers
    allocate(intObsFamilyListMpiLocal(len(currentObsFamily),obsFamilyListSizeMaxMpiLocal))
    intObsFamilyListMpiLocal(:,:)=0
    do familyIndex = 1, obsFamilyListSizeMpiLocal
      do charIndex = 1, len(currentObsFamily)
        intObsFamilyListMpiLocal(charIndex,familyIndex) =  &
             iachar(obsFamilyListMpiLocal(familyIndex)(charIndex:charIndex))
      end do
    end do

    ! communicate obs family list to all mpi tasks as integers
    allocate(intObsFamilyListMpiGlobal(len(currentObsFamily),obsFamilyListSizeMaxMpiLocal,mmpi_nprocs))
    nsize = size(intObsFamilyListMpiLocal)
    call rpn_comm_allgather(intObsFamilyListMpiLocal,  nsize, 'mpi_integer',  &
                            intObsFamilyListMpiGlobal, nsize, 'mpi_integer', 'GRID', ierr)

    ! convert global family lists from integers to characters
    allocate(obsFamilyListMpiGlobal(obsFamilyListSizeMaxMpiLocal,mmpi_nprocs))
    obsFamilyListMpiGlobal(:,:) = 'XX'
    do procIndex = 1, mmpi_nprocs
      do familyIndex = 1, allObsFamilyListSizeMpiLocal(procIndex)
        do charIndex=1,len(currentObsFamily)
          obsFamilyListMpiGlobal(familyIndex,procIndex)(charIndex:charIndex) =  &
               achar(intObsFamilyListMpiGlobal(charIndex,familyIndex,procIndex))
        end do
      end do
      write(*,*) 'obsFamilyListMpiGlobal = ', procIndex,  &
           obsFamilyListMpiGlobal(1:allObsFamilyListSizeMpiLocal(procIndex),procIndex)
    end do

    ! construct single common list of families to be used for all mpi tasks
    obsFamilyListCommon(:) = 'YY'
    obsFamilyListSizeCommon = 0
    do procIndex = 1, mmpi_nprocs
      FAMILY: do familyIndex = 1, obsFamilyListSizeMaxMpiLocal
        if (obsFamilyListMpiGlobal(familyIndex,procIndex) == 'XX') cycle FAMILY
        if (any(obsFamilyListCommon(:) == obsFamilyListMpiGlobal(familyIndex,procIndex))) cycle FAMILY
        obsFamilyListSizeCommon = obsFamilyListSizeCommon + 1
        obsFamilyListCommon(obsFamilyListSizeCommon) = obsFamilyListMpiGlobal(familyIndex,procIndex)
      end do FAMILY
    end do
    write(*,*) 'obsFamilyListSizeCommon = ', obsFamilyListSizeCommon
    write(*,*) 'obsFamilyListCommon = ', obsFamilyListCommon(1:obsFamilyListSizeCommon)

    deallocate(allObsFamilyListSizeMpiLocal)
    deallocate(obsFamilyListMpiGlobal)
    deallocate(intObsFamilyListMpiGlobal)
    deallocate(intObsFamilyListMpiLocal)
    deallocate(obsFamilyListMpiLocal)

  end subroutine diaf_getObsFamilyListMpiGlobal

end module obsDiagFiles_mod
