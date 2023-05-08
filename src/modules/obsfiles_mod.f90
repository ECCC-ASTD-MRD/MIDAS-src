
module obsFiles_mod
  ! MODULE obsFiles_mod (prefix='obsf' category='3. Observation input/output')
  !
  ! :Purpose: High-level module to handle reading/writing of observations that
  !           can be stored in one of several different formats. Currently, the
  !           only supported formats are:
  !
  !              1. BURP
  !              2. SQLITE (burp2rdb format)
  !              3. SQLITE (obsDB format)
  !
  use midasMpi_mod
  use ramdisk_mod
  use utilities_mod
  use obsSpaceData_mod
  use burpFiles_mod
  use sqliteFiles_mod
  use sqliteUtilities_mod
  use obsdbFiles_mod
  use obsDiagFiles_mod
  use bufr_mod
  use obsSubSpaceData_mod
  use obsUtil_mod
  use obsVariableTransforms_mod
  use burpread_mod
  use biasCorrectionConv_mod
  use clib_interfaces_mod
  use tovs_nl_mod 
  use ensembleobservations_mod

  implicit none
  save
  private

  ! public variables
  public :: obsf_nfiles, obsf_fileName

  ! public procedures
  public :: obsf_setup, obsf_filesSplit, obsf_determineFileType, obsf_determineSplitFileType
  public :: obsf_readFiles, obsf_writeFiles, obsf_obsSub_read, obsf_obsSub_update
  public :: obsf_addCloudParametersAndEmissivity, obsf_getFileName, obsf_copyObsDirectory
  public :: obsf_updateMissingObsFlags, obsf_cleanObsFiles
  logical           :: obsFilesSplit
  logical           :: initialized = .false.

  integer, parameter :: maxNumObsfiles = 150
  integer, parameter :: maxLengthFilename = 1060
  integer, parameter :: familyTypeLen = 2
  integer :: obsf_nfiles, obsf_numMpiUniqueList
  character(len=maxLengthFilename) :: obsf_fileName(maxNumObsfiles)
  character(len=familyTypeLen)     :: obsf_familyType(maxNumObsfiles)
  character(len=48)  :: obsFileMode
  character(len=maxLengthFilename) :: obsf_baseFileNameMpiUniqueList(maxNumObsfiles)
  character(len=familyTypeLen)     :: obsf_familyTypeMpiUniqueList(maxNumObsfiles)
  character(len=9) :: obsf_myIdExt

contains

  !--------------------------------------------------------------------------
  ! obsf_setup
  !--------------------------------------------------------------------------
  subroutine obsf_setup(dateStamp_out, obsFileMode_in)

    implicit none

    ! arguments
    integer         , intent(out)           :: dateStamp_out
    character(len=*), intent(in)            :: obsFileMode_in

    ! locals
    character(len=10)                       :: obsFileType

    obsFileMode = trim(obsFileMode_in)

    !
    ! Initialize file names and the file type
    !
    call obsf_setupFileNames()
    call obsf_determineFileType(obsFileType)

    !
    ! Determine if obsFiles are split
    !
    if ( obsFileType == 'OBSDB' .or. obsFileType == 'BURP' .or.  &
         obsFileType == 'SQLITE') then
      obsFilesSplit = .true.
    else
      call utl_abort('obsf_setup: invalid observation file type: ' // trim(obsFileType))
    end if

    !
    ! Do some setup of observation files
    !
    if ( obsFileType == 'BURP' ) then
      call brpf_getDateStamp( dateStamp_out, obsf_fileName(1) )
    else if ( obsFileType == 'OBSDB' ) then
      call odbf_getDateStamp( dateStamp_out, obsf_fileName(1) )
    else if ( obsFileType == 'SQLITE' ) then
      call sqlf_getDateStamp( dateStamp_out, obsf_fileName(1) )
    else
      dateStamp_out = -1
    end if

    initialized = .true.

  end subroutine obsf_setup

  !--------------------------------------------------------------------------
  ! obsf_filesSplit
  !--------------------------------------------------------------------------
  function obsf_filesSplit() result(obsFilesSplit_out)
    implicit none

    ! arguments
    logical :: obsFilesSplit_out

    if ( .not.initialized ) call utl_abort('obsf_filesSplit: obsFiles_mod not initialized!')

    obsFilesSplit_out = obsFilesSplit
        
  end function obsf_filesSplit

  !--------------------------------------------------------------------------
  ! obsf_readFiles
  !--------------------------------------------------------------------------
  subroutine obsf_readFiles( obsSpaceData )

    implicit none

    ! arguments
    type(struct_obs)  :: obsSpaceData

    ! locals
    integer :: fileIndex, fileIndexMpiLocal, numHeaders, numBodies
    integer :: numHeaderBefore, numBodyBefore, numHeaderRead, numBodyRead
    character(len=10) :: obsFileType
    character(len=maxLengthFilename) :: fileName
    character(len=256) :: fileNamefull
    character(len=familyTypeLen) :: obsFamilyType
    logical :: fileExists

    if ( .not.initialized ) call utl_abort('obsf_readFiles: obsFiles_mod not initialized!')

    call obsf_determineFileType(obsFileType)

    ! for every splitted file, the file type is defined separately 
    do fileIndex = 1, obsf_numMpiUniqueList

      fileName = trim(obsf_baseFileNameMpiUniqueList(fileIndex)) // '_' // trim(obsf_myIdExt)
      obsFamilyType = obsf_familyTypeMpiUniqueList(fileIndex)
      fileNameFull = ram_fullWorkingPath(fileName,noAbort_opt=.true.)
      inquire(file=trim(fileNameFull),exist=fileExists)

      if (fileExists) then
        ! get fileIndex on mpi local
        do fileIndexMpiLocal = 1, obsf_nfiles
          if (trim(obsf_fileName(fileIndexMpiLocal)) == fileNamefull .and. &
              trim(obsf_familyType(fileIndexMpiLocal)) == obsFamilyType) then
            exit
          end if
        end do

        call obsf_determineSplitFileType( obsFileType, fileNameFull )
        if ( obsFileType == 'BURP' )   then
          ! Add extra bias correction elements to conventional and TO files
          ! Bias correction elements for AI are added at the derivate file stage
          if ( obsFamilyType == 'TO' ) then
            call brpr_addElementsToBurp(fileNameFull, obsFamilyType, beSilent_opt=.false.)
          else 
            if ( bcc_biasActive(obsFamilyType) ) then
              call brpr_addElementsToBurp(fileNameFull, obsFamilyType, beSilent_opt=.false.)
            end if
          end if

          numHeaderBefore = obs_numHeader(obsSpaceData)
          numBodyBefore = obs_numBody(obsSpaceData)
          call brpf_readFile( obsSpaceData, fileNameFull, obsFamilyType, fileIndexMpiLocal )
          numHeaders = obs_numHeader(obsSpaceData)
          numBodies = obs_numBody(obsSpaceData)
          numHeaderRead = numHeaders - numHeaderBefore
          numBodyRead = numBodies - numBodyBefore
        end if
        if ( obsFileType == 'SQLITE' ) then
          call sqlf_readFile( obsSpaceData, fileNameFull, obsFamilyType, fileIndexMpiLocal )
        end if
        if ( obsFileType == 'OBSDB' ) then
          call odbf_readFile( obsSpaceData, fileNameFull, obsFamilyType, fileIndexMpiLocal )
        end if
      else
        numHeaderRead = 0
        numBodyRead = 0
      end if ! if (fileExists)

      if (obsFileType == 'BURP') then
        call setHeadBodyPrimaryKeyColumns(obsSpaceData, numHeaderRead, numBodyRead)
      end if
      
    end do

    ! abort if NAMTOV does not exist but there are radiance observation files
    if ( .not. utl_isNamelistPresent('NAMTOV','./flnml') .and. &
        any(obsf_familyType(:) == 'TO') ) then
      call utl_abort('obsf_readFiles: Namelist block NAMTOV is missing but there are radiance observation files')
    end if

    ! initialize OBS_HIND for each observation
    call obs_sethind(obsSpaceData)
  
  end subroutine obsf_readFiles

  !--------------------------------------------------------------------------
  ! obsf_writeFiles
  !--------------------------------------------------------------------------
  subroutine obsf_writeFiles( obsSpaceData, HXens_mpiglobal_opt, asciDumpObs_opt, writeDiagFiles_opt, ensObs_opt)
    implicit none

    ! arguments
    type(struct_obs)           :: obsSpaceData
    real(8), optional          :: HXens_mpiglobal_opt(:,:)
    logical, optional          :: asciDumpObs_opt
    logical, optional          :: writeDiagFiles_opt
    type(struct_eob), optional :: ensObs_opt          
  
    ! locals
    integer           :: fileIndex, fnom, fclos, nulnam, ierr
    integer           :: status, baseNameIndexBeg
    character(len=200) :: baseName, fullName, fileNameDir
    character(len=10) :: obsFileType, sfFileName
    character(len=*), parameter :: myName = 'obsf_writeFiles'
    character(len=*), parameter :: myWarning = myName //' WARNING: '

    ! namelist variables
    logical :: lwritediagsql ! choose to write 'diag' sqlite observation files
    logical :: onlyAssimObs  ! choose to not include unassimilated obs in 'diag' sqlite files
    logical :: addFSOdiag    ! choose to include FSO column in body table
    logical :: writeObsDb    ! write obDB file from scratch

    namelist /namwritediag/ lwritediagsql, onlyAssimObs, addFSOdiag, writeObsDb

    call utl_tmg_start(10,'--Observations')

    if ( .not.initialized ) call utl_abort('obsf_writeFiles: obsFiles_mod not initialized!')
 
    call obsf_determineFileType(obsFileType)

    nulnam=0
    lwritediagsql = .false.
    onlyAssimObs = .false.
    addFSOdiag = .false.
    writeObsDb = .false.
    ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=namwritediag,iostat=ierr)
    if (ierr /= 0) write(*,*) myWarning//' namwritediag is missing in the namelist. The default value will be taken.'
    if (mmpi_myid == 0) write(*,nml = namwritediag)
    if (present(writeDiagFiles_opt)) then
      lwritediagsql = lwritediagsql .and. writeDiagFiles_opt
    end if
    ierr=fclos(nulnam)

    if ( obsFileType == 'BURP' .or. obsFileType == 'OBSDB' .or. &
         obsFileType == 'SQLITE' ) then

      if (trim(obsFileMode) == 'analysis') call ovt_transformResiduals(obsSpaceData, obs_oma)
      if (trim(obsFileMode) /= 'prepcma' .and. trim(obsFileMode) /= 'thinning') then
        call ovt_transformResiduals(obsSpaceData, obs_omp)
      end if
      if (trim(obsFileMode) == 'analysis' .or. trim(obsFileMode) == 'FSO') call obsu_setassflg(obsSpaceData)
      if (trim(obsFileMode) /= 'prepcma' .and. trim(obsFileMode) /= 'thinning') then
        call obsu_updateSourceVariablesFlag(obsSpaceData)
      end if
      if (trim(obsFileMode) /= 'prepcma') call ovt_transformResiduals(obsSpaceData, obs_omp)
      if (trim(obsFileMode) /= 'prepcma') call obsu_updateSourceVariablesFlag(obsSpaceData)

      do fileIndex = 1, obsf_nfiles
        call obsf_determineSplitFileType( obsFileType, obsf_fileName(fileIndex) )

        if ( obsFileType == 'BURP'   ) then
          call brpf_updateFile( obsSpaceData, obsf_fileName(fileIndex), obsf_familyType(fileIndex), &
                                fileIndex )
        else if ( obsFileType == 'OBSDB' ) then
          call odbf_updateFile( obsSpaceData, obsf_fileName(fileIndex), &
                                obsf_familyType(fileIndex), fileIndex )
        else if ( obsFileType == 'SQLITE' ) then
          call sqlf_updateFile( obsSpaceData, obsf_fileName(fileIndex), &
                                obsf_familyType(fileIndex), fileIndex )
        end if
      end do

      if ( present(HXens_mpiglobal_opt) .and. mmpi_myid == 0 ) then
        call obsf_writeHX(obsSpaceData, HXens_mpiglobal_opt)
      end if

      if (obsFileType /= 'OBSDB' .and. writeObsDb) then
        ! Create destination directory
        fileNameDir = trim(ram_getRamDiskDir())
        if (fileNameDir == ' ') then
          write(*,*) 'obsf_writeFiles: WARNING! Writing obsDB to current working path ' // &
                     'instead of ramdisk'
        end if

        if (mmpi_myid == 0) status = clib_mkdir_r(trim(fileNameDir)//'obsDB')
        if (obsf_filesSplit()) call rpn_comm_barrier('GRID',status)

        ! update obsDB files
        do fileIndex = 1, obsf_nfiles
          fullName = trim(obsf_fileName(fileIndex))
          baseNameIndexBeg = index(fullName,'/',back=.true.)
          baseName = fullName(baseNameIndexBeg+1:)

          call odbf_updateFile(obsSpaceData, trim(fileNameDir)//'obsDB/'//trim(baseName), &
                               obsf_familyType(fileIndex), fileIndex)
        end do        
      end if

    end if

    if ( index(obsf_getFileName('SF'), 'sfc') > 0 ) then
      sfFileName = 'sfc'
    else
      sfFileName = 'sf'
    end if

    if (lwritediagsql) then 
      call diaf_writeAllSqlDiagFiles(obsSpaceData, sfFileName, onlyAssimObs, &
                                     addFSOdiag, ensObs_opt=ensObs_opt)
    end if

    if ( present(asciDumpObs_opt) ) then
      if ( asciDumpObs_opt ) then
        if ( obsFileType == 'BURP' .or. obsFileType == 'OBSDB' .or.  &
             obsFileType == 'SQLITE' .or. mmpi_myid == 0   ) then
          ! all processors write to files only for BURP, OBSDB and SQLITE    
          call obsf_writeAsciDump(obsSpaceData)
        end if
      end if
    end if

    call utl_tmg_stop(10)

  end subroutine obsf_writeFiles

  !--------------------------------------------------------------------------
  ! obsf_cleanObsFiles
  !--------------------------------------------------------------------------
  subroutine obsf_cleanObsFiles()
    implicit none

    ! locals
    integer           :: fileIndex
    character(len=10) :: obsFileType

    if ( .not.initialized ) call utl_abort('obsf_cleanObsFiles: obsFiles_mod not initialized!')
   
    call obsf_determineFileType(obsFileType)

    if ( obsFileType /= 'BURP' .and. obsFileType /= 'OBSDB' .and.  &
         obsFileType /= 'SQLITE' ) then
      write(*,*) 'obsf_cleanObsFiles: obsFileType=', obsFileType, &
                ' is not BURP, OBSDB nor SQLITE. Return.' 
      return
    end if
    
    call utl_tmg_start(23, '----ObsFileClean')

    do fileIndex = 1, obsf_nfiles
      call obsf_determineSplitFileType( obsFileType, obsf_fileName(fileIndex) )

      if ( obsFileType == 'BURP' ) then
        call brpr_burpClean( obsf_fileName(fileIndex), obsf_familyType(fileIndex) )
      else if ( obsFileType == 'OBSDB' ) then
        call obdf_clean( obsf_fileName(fileIndex), obsf_familyType(fileIndex) )
      else if ( obsFileType == 'SQLITE' ) then
        call sqlf_cleanFile( obsf_fileName(fileIndex), obsf_familyType(fileIndex) )
      end if
    end do

    call utl_tmg_stop(23)

  end subroutine obsf_cleanObsFiles 

  !--------------------------------------------------------------------------
  ! obsf_writeHX
  !--------------------------------------------------------------------------
  subroutine obsf_writeHX(obsSpaceData, HXens_mpiglobal)
    implicit none

    ! arguments
    type(struct_obs) :: obsSpaceData
    real(8)          :: HXens_mpiglobal(:,:)

    ! locals
    integer :: unitHX, ierr, headerIndex, fnom, fclos
    character(len=10) :: fileNameHX

    write(*,*) 'obsf_writeHX: Starting'

    fileNameHX     = 'cmahxout'
    unitHX = 0
    ierr = fnom(unitHX, fileNameHX, 'FTN+SEQ+UNF+R/W', 0)

    do headerIndex = 1, obs_numHeader(obsSpaceData)

      call obs_write_hx(obsSpaceData, HXens_mpiglobal, headerIndex, unitHX)

    enddo
 
    ierr = fclos(unitHX)

  end subroutine obsf_writeHX

  !--------------------------------------------------------------------------
  ! obsf_writeAsciDump
  !--------------------------------------------------------------------------
  subroutine obsf_writeAsciDump(obsSpaceData)

    implicit none

    ! arguments
    type(struct_obs) :: obsSpaceData

    ! locals
    character(len=25) :: fileNameAsciDump
    integer :: unitAsciDump, ierr, fnom, fclos
    character(len=4)    :: cmyidx, cmyidy
    character(len=9)    :: cmyid

    write(*,*) 'obsf_writeAsciDump: Starting'

    ! determine the file name depending on if obs data is mpi local or global
    if ( obs_mpiLocal(obsSpaceData) ) then
      ! separate file per mpi task
      write(cmyidy,'(I4.4)') (mmpi_myidy + 1)
      write(cmyidx,'(I4.4)') (mmpi_myidx + 1)
      cmyid  = trim(cmyidx) // '_' // trim(cmyidy)
      fileNameAsciDump = 'obsout_asci_' // trim(cmyid)
    else
      ! only task 0 write the global data
      if ( mmpi_myid > 0 ) return
      fileNameAsciDump = 'obsout_asci'
    end if

    write(*,*) 'obsf_writeAsciDump: writing to file : ', fileNameAsciDump
    unitAsciDump = 0
    ierr = fnom(unitAsciDump, fileNameAsciDump, 'FTN+SEQ+FMT+R/W', 0)
    call obs_print(obsSpaceData,unitAsciDump)
    ierr = fclos(unitAsciDump)

  end subroutine obsf_writeAsciDump

  !--------------------------------------------------------------------------
  ! obsf_setupFileNames
  !--------------------------------------------------------------------------
  subroutine obsf_setupFileNames()

    implicit none

    ! locals
    character(len=20) :: namePrefix(maxNumObsfiles)
    character(len=2) :: familyName(maxNumObsfiles)
    character(len=4) :: cmyidx, cmyidy
    character(len=256):: obsDirectory
    character(len=maxLengthFilename) :: fileName, baseFileName  ! the length should be more than 
                                                                      ! len(obsDirectory)+1+len(namePrefix)+1+len(obsf_myIdExt)
    character(len=maxLengthFilename) :: baseFileNameList(maxNumObsfiles)

    character(len=256)               :: fileNamefull
    logical :: fileExists
    integer :: fileIndex

    write(cmyidy,'(I4.4)') (mmpi_myidy + 1)
    write(cmyidx,'(I4.4)') (mmpi_myidx + 1)
    obsf_myIdExt  = trim(cmyidx) // '_' // trim(cmyidy)

    namePrefix(:)=''
    ! file names only for burp
    namePrefix( 1) = 'brpuan'
    namePrefix( 2) = 'brpuas'
    namePrefix( 3) = 'brpai'
    namePrefix( 4) = 'brpain'
    namePrefix( 5) = 'brpais'
    namePrefix( 6) = 'brpaie'
    namePrefix( 7) = 'brpaiw'
    namePrefix( 8) = 'brpsfc'
    namePrefix( 9) = 'brpsf'
    namePrefix(10) = 'brptov'
    namePrefix(11) = 'brpssmis'
    namePrefix(12) = 'brpairs'
    namePrefix(13) = 'brpto_amsua'
    namePrefix(14) = 'brpto_amsua_allsky'
    namePrefix(15) = 'brpto_amsub'
    namePrefix(16) = 'brpto_amsub_allsky'
    namePrefix(17) = 'brpcsr'
    namePrefix(18) = 'brpiasi'
    namePrefix(19) = 'brpatms'
    namePrefix(20) = 'brpatms_allsky'
    namePrefix(21) = 'brpcris'
    namePrefix(22) = 'brpcrisfsr'
    namePrefix(23) = 'brpcrisfsr1'
    namePrefix(24) = 'brpcrisfsr2'
    namePrefix(25) = 'brpcrisfsr3'
    namePrefix(26) = 'brpcrisfsr4'
    namePrefix(27) = 'brpsw'
    namePrefix(28) = 'brpswgoes9'
    namePrefix(29) = 'brpswgoese'
    namePrefix(30) = 'brpswgoesw'
    namePrefix(31) = 'brpswmodis'
    namePrefix(32) = 'brpswmtsate'
    namePrefix(33) = 'brpswmtsatw'
    namePrefix(34) = 'brpgo'
    namePrefix(35) = 'brpsc'
    namePrefix(36) = 'brppr'
    namePrefix(37) = 'brpro'
    namePrefix(38) = 'brphum'
    namePrefix(39) = 'brpsat'
    namePrefix(40) = 'brpssm'
    namePrefix(41) = 'brpgp'
    namePrefix(42) = 'brpch'
    namePrefix(43) = 'brpua'
    ! new general file names for burp and sqlite
    namePrefix(44) = 'obsua'
    namePrefix(45) = 'obsuan'
    namePrefix(46) = 'obsuas'
    namePrefix(47) = 'obsai'
    namePrefix(48) = 'obsain'
    namePrefix(49) = 'obsais'
    namePrefix(50) = 'obsaie'
    namePrefix(51) = 'obsaiw'
    namePrefix(52) = 'obssfc'
    namePrefix(53) = 'obssf'
    namePrefix(54) = 'obstov'
    namePrefix(55) = 'obsssmis'
    namePrefix(56) = 'obsairs'
    namePrefix(57) = 'obsto_amsua'
    namePrefix(58) = 'obsto_amsua_allsky'
    namePrefix(59) = 'obsto_amsub'
    namePrefix(60) = 'obsto_amsub_allsky'
    namePrefix(61) = 'obscsr'
    namePrefix(62) = 'obsiasi'
    namePrefix(63) = 'obsatms'
    namePrefix(64) = 'obsatms_allsky'
    namePrefix(65) = 'obscris'
    namePrefix(66) = 'obscrisfsr'
    namePrefix(67) = 'obscrisfsr1'
    namePrefix(68) = 'obscrisfsr2'
    namePrefix(69) = 'obscrisfsr3'
    namePrefix(70) = 'obscrisfsr4'
    namePrefix(71) = 'obssw'
    namePrefix(72) = 'obsswgoes9'
    namePrefix(73) = 'obsswgoese'
    namePrefix(74) = 'obsswgoesw'
    namePrefix(75) = 'obsswmodis'
    namePrefix(76) = 'obsswmtsate'
    namePrefix(77) = 'obsswmtsatw'
    namePrefix(78) = 'obsgo'
    namePrefix(79) = 'obssc'
    namePrefix(80) = 'obspr'
    namePrefix(81) = 'obsro'
    namePrefix(82) = 'obshum'
    namePrefix(83) = 'obssat'
    namePrefix(84) = 'obsssm'
    namePrefix(85) = 'obsgp'
    namePrefix(86) = 'obsch'
    namePrefix(87) = 'obsgl_ssmi'
    namePrefix(88) = 'obsgl_ssmis'
    namePrefix(89) = 'obsgl_amsr2'
    namePrefix(90) = 'obsgl_ascat'
    namePrefix(91) = 'obsgl_avhrr'
    namePrefix(92) = 'obsgl_cisA'
    namePrefix(93) = 'obsgl_cisI'
    namePrefix(94) = 'obsgl_cisL'
    namePrefix(95) = 'obsgl_cisR'
    namePrefix(96) = 'cmaheader' ! file name for CMA format used by EnKF
    namePrefix(97) = 'brpsst'
    namePrefix(98) = 'obsal'
    namePrefix(99) = 'obsradar'
    namePrefix(100) = 'obssst_insitu'
    namePrefix(101) = 'obshydro'
    namePrefix(102)= 'obsmwhs2'
    namePrefix(103)= 'brpmwhs2'
    namePrefix(104)= 'obssarwinds'
    namePrefix(105)= 'obssst_avhrr'
    namePrefix(106)= 'obssst_amsr2'
    namePrefix(107)= 'obssst_viirs'
    namePrefix(108)= 'obssst_pseudo'
    namePrefix(109)= 'obssst'

    familyName(:)   = ''
    familyName( 1)  = 'UA'
    familyName( 2)  = 'UA'
    familyName( 3)  = 'AI'
    familyName( 4)  = 'AI'
    familyName( 5)  = 'AI'
    familyName( 6)  = 'AI'
    familyName( 7)  = 'AI'
    familyName( 8)  = 'SF'
    familyName( 9)  = 'SF'
    familyName(10)  = 'TO'
    familyName(11)  = 'TO'
    familyName(12)  = 'TO'
    familyName(13)  = 'TO'
    familyName(14)  = 'TO'
    familyName(15)  = 'TO'
    familyName(16)  = 'TO'
    familyName(17)  = 'TO'
    familyName(18)  = 'TO'
    familyName(19)  = 'TO'
    familyName(20)  = 'TO'
    familyName(21)  = 'TO'
    familyName(22)  = 'TO'
    familyName(23)  = 'TO'
    familyName(24)  = 'TO'
    familyName(25)  = 'TO'
    familyName(26)  = 'TO'
    familyName(27)  = 'SW'
    familyName(28)  = 'SW'
    familyName(29)  = 'SW'
    familyName(30)  = 'SW'
    familyName(31)  = 'SW'
    familyName(32)  = 'SW'
    familyName(33)  = 'SW'
    familyName(34)  = 'GO'
    familyName(35)  = 'SC'
    familyName(36)  = 'PR'
    familyName(37)  = 'RO'
    familyName(38)  = 'HU'
    familyName(39)  = 'ST'
    familyName(40)  = 'MI'
    familyName(41)  = 'GP'
    familyName(42)  = 'CH'
    familyName(43)  = 'UA'
    familyName(44)  = 'UA'
    familyName(45)  = 'UA'
    familyName(46)  = 'UA'
    familyName(47)  = 'AI'
    familyName(48)  = 'AI'
    familyName(49)  = 'AI'
    familyName(50)  = 'AI'
    familyName(51)  = 'AI'
    familyName(52)  = 'SF'
    familyName(53)  = 'SF'
    familyName(54)  = 'TO'
    familyName(55)  = 'TO'
    familyName(56)  = 'TO'
    familyName(57)  = 'TO'
    familyName(58)  = 'TO'
    familyName(59)  = 'TO'
    familyName(60)  = 'TO'
    familyName(61)  = 'TO'
    familyName(62)  = 'TO'
    familyName(63)  = 'TO'
    familyName(64)  = 'TO'
    familyName(65)  = 'TO'
    familyName(66)  = 'TO'
    familyName(67)  = 'TO'
    familyName(68)  = 'TO'
    familyName(69)  = 'TO'
    familyName(70)  = 'TO'
    familyName(71)  = 'SW'
    familyName(72)  = 'SW'
    familyName(73)  = 'SW'
    familyName(74)  = 'SW'
    familyName(75)  = 'SW'
    familyName(76)  = 'SW'
    familyName(77)  = 'SW'
    familyName(78)  = 'GO'
    familyName(79)  = 'SC'
    familyName(80)  = 'PR'
    familyName(81)  = 'RO'
    familyName(82)  = 'HU'
    familyName(83)  = 'ST'
    familyName(84)  = 'MI'
    familyName(85)  = 'GP'
    familyName(86)  = 'CH'
    familyName(87)  = 'GL'
    familyName(88)  = 'GL'
    familyName(89)  = 'GL'
    familyName(90)  = 'GL'
    familyName(91)  = 'GL'
    familyName(92)  = 'GL'
    familyName(93)  = 'GL'
    familyName(94)  = 'GL'
    familyName(95)  = 'GL'
    familyName(96)  = 'XX' ! dummy family type for CMA, since it contains all families
    familyName(97)  = 'TM'
    familyName(98)  = 'AL'
    familyName(99)  = 'RA'
    familyName(100)  = 'TM'
    familyName(101)  = 'HY'
    familyName(102) = 'TO'
    familyName(103) = 'TO'
    familyName(104) = 'SF'
    familyName(105) = 'TM'
    familyName(106) = 'TM'
    familyName(107) = 'TM'
    familyName(108) = 'TM'
    familyName(109) = 'TM'

    obsDirectory = 'obs'

    obsf_nfiles = 0
    obsf_fileName(1) = 'DUMMY_FILE_NAME'
    baseFileNameList(:) = ''

    do fileIndex = 1, maxNumObsfiles 

      if(namePrefix(fileIndex) == '') exit

      baseFileName = trim(obsDirectory) // '/' // trim(namePrefix(fileIndex))
      fileName = trim(baseFileName) // '_' // trim(obsf_myIdExt)
      fileNameFull = ram_fullWorkingPath(fileName,noAbort_opt=.true.)
      inquire(file=trim(fileNameFull),exist=fileExists)

      if (.not. fileExists ) then
        fileName=trim(baseFileName)
        fileNameFull = ram_fullWorkingPath(fileName, noAbort_opt=.true.)
        inquire(file=trim(fileNameFull), exist=fileExists)
      end if

      if ( fileExists ) then
        obsf_nfiles=obsf_nfiles + 1
        baseFileNameList(obsf_nfiles) = trim(baseFileName)
        obsf_fileName(obsf_nfiles) = fileNameFull
        obsf_familyType(obsf_nfiles) = familyName(fileIndex)
      end if

    end do

    call setObsFilesMpiUniqueList(baseFileNameList)
    
    write(*,*) ' '
    write(*,*)'obsf_setupFileNames: Number of observation files is :', obsf_nfiles
    write(*,*)'Type  Name '
    write(*,*)'----  ---- '
    do fileIndex = 1, obsf_nfiles
      write(*,'(1X,A2,1X,A60)' ) obsf_familyType(fileIndex), trim(obsf_fileName(fileIndex))
    end do

  end subroutine obsf_setupFileNames

  !--------------------------------------------------------------------------
  ! setObsFilesMpiUniqueList
  !--------------------------------------------------------------------------
  subroutine setObsFilesMpiUniqueList(baseFileNameList)
    !
    ! :Purpose: Create a unique list of obs filenames/familyTypes across all mpi tasks.
    !
    implicit none

    ! Arguments
    character(len=*), intent(in) :: baseFileNameList(:)

    ! Locals:
    integer :: fileIndex, fileIndex2, procIndex, ierr
    character(len=maxLengthFilename), allocatable :: baseFileNameListAllMpi(:,:)
    character(len=familyTypeLen), allocatable :: familyTypeListAllMpi(:,:)

    ! Communicate filenames across all mpi tasks
    allocate(baseFileNameListAllMpi(maxNumObsfiles,mmpi_nprocs))
    baseFileNameListAllMpi(:,:) = ''
    call mmpi_allgather_string(baseFileNameList, baseFileNameListAllMpi, &
                               maxNumObsfiles, maxLengthFilename, mmpi_nprocs, &
                               "GRID", ierr)

    ! Communicate familyTypes across all mpi tasks
    allocate(familyTypeListAllMpi(maxNumObsfiles,mmpi_nprocs))
    familyTypeListAllMpi(:,:) = ''
    call mmpi_allgather_string(obsf_familyType, familyTypeListAllMpi, &
                               maxNumObsfiles, familyTypeLen, mmpi_nprocs, &
                               "GRID", ierr)

    ! Create a unique list of obs filenames/familytype across all mpi tasks without duplicates
    obsf_baseFileNameMpiUniqueList(:) = ''
    obsf_familyTypeMpiUniqueList(:) = ''
    obsf_numMpiUniqueList = 1
    obsf_baseFileNameMpiUniqueList(obsf_numMpiUniqueList) = baseFileNameListAllMpi(1,1)
    obsf_familyTypeMpiUniqueList(obsf_numMpiUniqueList) = familyTypeListAllMpi(1,1)
    do procIndex = 1, mmpi_nprocs
      loopFilename: do fileIndex = 1, maxNumObsfiles 
        if (trim((baseFileNameListAllMpi(fileIndex,procIndex))) == '') cycle loopFilename

        ! cycle if filename already exists in the unique list
        do fileIndex2 = 1, obsf_numMpiUniqueList
          if (trim(obsf_baseFileNameMpiUniqueList(fileIndex2)) == &
              trim(baseFileNameListAllMpi(fileIndex,procIndex))) then
            cycle loopFilename
          end if
        end do
        
        ! add the filename to the unique list
        obsf_numMpiUniqueList = obsf_numMpiUniqueList + 1
        obsf_baseFileNameMpiUniqueList(obsf_numMpiUniqueList) = baseFileNameListAllMpi(fileIndex,procIndex)
        obsf_familyTypeMpiUniqueList(obsf_numMpiUniqueList) = familyTypeListAllMpi(fileIndex,procIndex)

      end do loopFilename
    end do

    write(*,*) 'setObsFilesMpiUniqueList: obsf_numMpiUniqueList=', obsf_numMpiUniqueList
    write(*,*) 'setObsFilesMpiUniqueList: familyType/filename in unique list:'
    write(*,*) 'Type  Name '
    write(*,*) '----  ---- '
    do fileIndex = 1, obsf_numMpiUniqueList
      write(*,'(1X,A2,1X,A60)' ) trim(obsf_familyTypeMpiUniqueList(fileIndex)), &
                                 trim(obsf_baseFileNameMpiUniqueList(fileIndex))
    end do

  end subroutine setObsFilesMpiUniqueList
    
  !--------------------------------------------------------------------------
  ! obsf_determineFileType
  !--------------------------------------------------------------------------
  subroutine obsf_determineFileType( obsFileType )

    implicit none

    ! arguments
    character(len=*)                       :: obsFileType

    ! locals
    integer :: ierr, procID, all_nfiles(0:(mmpi_nprocs-1))
    logical :: fileExists

    call rpn_comm_allgather( obsf_nfiles, 1, 'MPI_INTEGER', &
                             all_nfiles,  1, 'MPI_INTEGER', 'GRID', ierr )
    fileExists = .false.
    procid_loop: do procID = 0, (mmpi_nprocs-1)
      if ( all_nfiles(procID) > 0 ) then
        fileExists = .true.
        exit procid_loop
      end if
    end do procid_loop

    if ( .not.fileExists ) then
      call utl_abort('obsf_determineFileType: No observation files found')
    end if

    write(*,*) 'obsf_determineFileType: read obs file that exists on mpi task id: ', procID

    if ( mmpi_myid == procID ) call obsf_determineSplitFileType( obsFileType, obsf_fileName(1) )

    call rpn_comm_bcastc(obsFileType , len(obsFileType), 'MPI_CHARACTER', procID, 'GRID', ierr)
    write(*,*) 'obsf_determineFileType: obsFileType = ', obsFileType

  end subroutine obsf_determineFileType

  !--------------------------------------------------------------------------
  ! obsf_determineSplitFileType
  !--------------------------------------------------------------------------
  subroutine obsf_determineSplitFileType( obsFileType, fileName )

    implicit none

    ! arguments
    character(len=*),intent(out) :: obsFileType
    character(len=*),intent(in)  :: fileName

    ! locals
    integer           :: ierr, unitFile
    integer           :: fnom, fclos, wkoffit
    character(len=20) :: fileStart
    character(len=*), parameter :: obsDbTableName = 'Report'

    write(*,*) 'obsf_determineSplitFileType: read obs file: ', trim(fileName)
   
    ierr = wkoffit(trim(fileName))

    if (ierr.eq.6) then
      obsFiletype = 'BURP'
    else

      unitFile = 0
      ierr = fnom(unitFile, fileName, 'FTN+SEQ+FMT+R/O', 0)
      read(unitFile,'(A)') fileStart
      ierr = fclos(unitFile)

      if ( index( fileStart, 'SQLite format 3' ) > 0 ) then
        if (sqlu_sqlTableExists(fileName, obsDbTableName)) then
          obsFileType = 'OBSDB'
        else
          obsFileType = 'SQLITE'
        end if
      else
        call utl_abort('obsf_determineSplitFileType: unknown obs file type')
      end if

    end if

    write(*,*) 'obsf_determineSplitFileType: obsFileType = ', obsFileType

  end subroutine obsf_determineSplitFileType

  !--------------------------------------------------------------------------
  ! obsf_getFileName
  !--------------------------------------------------------------------------
  function obsf_getFileName(obsfam,fileFound_opt) result(filename)
    !
    ! :Purpose: Returns the observations file name assigned to the calling processor.
    !           If the input family has more than one file, the first file found will
    !           be returned.
    !
    ! :Arguments:
    !           :obsfam: observation family name
    !           :found_opt: logical indicating if a file could be found for the family (optional)
    !
    implicit none
    ! arguments:
    character(len=2), intent(in) :: obsfam
    logical, intent(out), optional :: fileFound_opt
    character(len=maxLengthFilename) :: filename ! file name of associated observations file
    ! locals:
    integer :: numFound, ifile

    filename = ""
    numFound = 0
   
    do ifile=1,obsf_nfiles
       if (obsfam == obsf_familyType(ifile)) then
          filename = obsf_fileName(ifile)
          numFound = numFound + 1
          exit
       end if
    end do

    if (numFound == 0) then
      write(*,*) "obsf_getFileName: WARNING: File not found for obs family '" // trim(obsfam) // "'"
    end if
    if (numFound > 1) then
      write(*,*) "obsf_getFileName: WARNING: Multiple files found for obs family '" // trim(obsfam) // "'"
    end if

    if (present(fileFound_opt)) fileFound_opt = (numFound > 0)

  end function obsf_getFileName


  function obsf_obsSub_read( obsfam, stnid, varno, nlev, ndim, numColumns_opt, &
                             bkstp_opt, block_opt, match_nlev_opt, codtyp_opt ) result(obsdata)
    !
    ! :Purpose: Retrieves information for observations from observation files and returns the data
    !           in a struct_oss_obsdata object. Data will be retrieved for all nodes that have valid
    !           filenames for the specied observational family and combined into one struct_oss_obsdata
    !           if the observational files are split.
    !
    ! :Arguments:
    !           :obsfam:          observation family name
    !           :stnid:           station ID of observation
    !           :varno:           BUFR code (if <=0, to search through all codes to obtain first
    !                             between 10000 and 16000)
    !           :nlev:            number of levels in the observation
    !           :ndim:            number of dimensions for the retrieved data in
    !                             each report (e.g. ndim=1 for std, ndim=2 for
    !                             averaging kernel matrices) 
    !           :numColumns_opt:  Number of columns (if different from nlev and for ndim=2)
    !           :bkstp_opt:       bkstp number of requested block if BURP file type (optional)
    !           :block_opt:       block type of requested block if BURP file type (optional)
    !                             Valid values are 'DATA', 'INFO', '3-D', and 'MRQR', indicated
    !                             by the two rightmost bits of bknat.
    !           :match_nlev_opt:  determines if the report matching criteria includes checking
    !                             if the report number of levels is the same as the input
    !                             argument nlev (optional)
    !           :codtyp_opt:      optional CODTYP list for search (optional)
    !
    implicit none

    character(len=*), intent(in)           :: obsfam
    character(len=9), intent(in)           :: stnid
    integer         , intent(in)           :: varno
    integer         , intent(in)           :: nlev
    integer         , intent(in)           :: ndim
    integer         , intent(in), optional :: numColumns_opt ! Number of columns (if different from nlev and for ndim=2)
    integer         , intent(in), optional :: bkstp_opt
    integer         , intent(in), optional :: codtyp_opt(:)
    logical         , intent(in), optional :: match_nlev_opt
    character(len=4), intent(in), optional :: block_opt
    type(struct_oss_obsdata)               :: obsdata ! struct_oss_obsdata object

    character(len=maxLengthFilename) :: filename
    logical :: fileFound
    character(len=10) :: obsFileType

    filename = obsf_getFileName(obsfam,fileFound)

    if (fileFound) then
       call obsf_determineSplitFileType( obsFileType, filename )
       if (obsFileType=='BURP') then
          if (.not.present(block_opt)) &
               call utl_abort("obsf_obsSub_read: optional variable 'block_opt' is required for BURP observational files.")
          if (.not.present(numColumns_opt)) then
            obsdata = brpf_obsSub_read(filename,stnid,varno,nlev,ndim,block_opt,bkstp_opt=bkstp_opt, &
                                   match_nlev_opt=match_nlev_opt,codtyp_opt=codtyp_opt)
          else
            obsdata = brpf_obsSub_read(filename,stnid,varno,nlev,ndim,block_opt,bkstp_opt=bkstp_opt, &
                                       match_nlev_opt=match_nlev_opt,codtyp_opt=codtyp_opt,numColumns_opt=numColumns_opt)
          end if 
       else
          call utl_abort("obsf_obsSub_read: Only BURP observational files currently supported.")
       end if

    else

       write(*,*) "obsf_obsSub_read: No observational files found for family '" // trim(obsfam) // "' for this node."

       if (obsf_filesSplit()) then
          ! Must allocate obsdata so that it is available from ALL processors when
          ! requiring of rpn_comm_allgather via oss_obsdata_MPIallgather.
          if (ndim == 1) then
             call oss_obsdata_alloc(obsdata,1,dim1=nlev)
          else
             call oss_obsdata_alloc(obsdata,1,dim1=nlev,dim2_opt=nlev)
          end if
          obsdata%nrep=0
          write(*,*) "obsf_obsSub_read: Setting empty struct_oss_obsdata object for this node."
       else
          call utl_abort("obsf_obsSub_read: Abort since files are not split.")
       end if
    end if

    if (obsf_filesSplit()) call oss_obsdata_MPIallgather(obsdata)

  end function obsf_obsSub_read


  function obsf_obsSub_update( obsdata, obsfam, varno, bkstp_opt, block_opt, multi_opt ) &
                           result(nrep_modified)
    ! :Purpose: Add or modify data in observational files from data stored
    !           in a struct_oss_obsdata object.
    !
    ! :Arguments:
    !           :obsdata: Input struct_oss_obsdata object for varno.
    !           :obsfam:  observation family name
    !           :varno:   BUFR descriptors. Number of elements must be 
    !                     max(1,obsdata%dim2)
    !           :bkstp_opt: bkstp number of requested block if BURP file type (optional)
    !           :block_opt: block type of requested block if BURP file type (optional)
    !                       Valid values are 'DATA', 'INFO', '3-D', and 'MRQR', indicated
    !                       by the two rightmost bits of bknat.
    !           :multi_opt: Indicates if intended report are for 'UNI' or 'MULTI' level data (optional)
    !
    implicit none

    type(struct_oss_obsdata), intent(inout) :: obsdata
    character(len=*), intent(in) :: obsfam
    integer, intent(in) :: varno(:)
    integer, intent(in), optional :: bkstp_opt
    character(len=4), intent(in), optional :: block_opt
    character(len=*), intent(in), optional :: multi_opt
    integer :: nrep_modified    ! Number of modified reports

    integer :: ierr,nrep_modified_global
    character(len=maxLengthFilename) :: filename
    logical :: fileFound
    character(len=10) :: obsFileType

    filename = obsf_getFileName(obsfam,fileFound)

    if (fileFound) then
       if (obsf_filesSplit() .or. mmpi_myid == 0) then
          call obsf_determineSplitFileType( obsFileType, filename )
          if (obsFileType=='BURP') then
             if (.not.present(block_opt)) &
                  call utl_abort("obsf_obsSub_update: optional varaible 'block_opt' is required for BURP observational files.")
             nrep_modified = brpf_obsSub_update(obsdata,filename,varno,block_opt,bkstp_opt=bkstp_opt,multi_opt=multi_opt)
          else
             call utl_abort("obsf_obsSub_update: Only BURP observational files currently supported.")
          end if
       end if
    else
       write(*,*) "obsf_obsSub_update: No observational files found for family '" // trim(obsfam) // "' for this node."
    end if

    if (obsf_filesSplit()) then
       call rpn_comm_allreduce(nrep_modified,nrep_modified_global,1,"MPI_INTEGER","MPI_SUM","GRID",ierr)
       nrep_modified = nrep_modified_global
    end if

  end function obsf_obsSub_update

  !--------------------------------------------------------------------------
  ! obsf_addCloudParametersAndEmissivity
  !--------------------------------------------------------------------------
  subroutine obsf_addCloudParametersAndEmissivity(obsSpaceData)
    !
    ! :Purpose: Loop on observation files to add cloud parameters and emissivity
    !
    implicit none
    ! Arguments:
    type(struct_obs),  intent(inout)  :: obsSpaceData

    ! Locals:
    integer           :: fileIndex
    character(len=10) :: obsFileType

    ! If obs files not split and I am not task 0, then return
    if ( .not.obsf_filesSplit() .and. mmpi_myid /= 0 ) return

    do fileIndex = 1, obsf_nfiles

      call obsf_determineSplitFileType( obsFileType, obsf_fileName(fileIndex) )
      if ( trim(obsFileType) == 'SQLITE' )  then
        call sqlf_addCloudParametersandEmissivity(obsSpaceData, fileIndex, obsf_fileName(fileIndex))
      else if ( trim(obsFileType) == 'BURP' ) then 
        call brpr_addCloudParametersandEmissivity(obsSpaceData, fileIndex, trim( obsf_fileName(fileIndex) ) )
      else if ( trim(obsFileType) == 'OBSDB' ) then
        ! The variables updated here for other file types are added/updated in ObsDb files using a
        ! seperate subroutine called odbf_updateMidasHeaderTable which is called when writing everything else.
        write(*,*) 'obsf_addCloudParametersAndEmissivity: obsDB files handled separately'
      else  
        write(*,*) ' UNKNOWN FileType=',obsFileType
        call utl_abort("obsf_addCloudParametersAndEmissivity: Only BURP, OBSDB and SQLITE observational files supported.")
      end if
    end do

  end subroutine obsf_addCloudParametersAndEmissivity

  !--------------------------------------------------------------------------
  ! obsf_updateMissingObsFlags
  !--------------------------------------------------------------------------
  subroutine obsf_updateMissingObsFlags(obsSpaceData)
    !
    ! :Purpose: Loop on observation files to set missing observation flags to 2048
    !           For now, this is done for only ATMS and AMSUA
    !
    implicit none
    ! Arguments:
    type(struct_obs)  :: obsSpaceData

    ! Locals:
    integer           :: fileIndex
    character(len=10) :: obsFileType
    logical           :: toDataPresent
    integer           :: headerIndex
    integer           :: codtyp 


    toDataPresent = .false.
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
    if (headerIndex < 0) exit HEADER
      codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if ( (tvs_isIdBurpInst(codtyp,'atms' )) .or. &
           (tvs_isIdBurpInst(codtyp,'amsua')) .or. &
           (tvs_isIdBurpInst(codtyp,'amsub')) .or. &
           (tvs_isIdBurpInst(codtyp,'mhs'  )) .or. & 
           (tvs_isIdBurpInst(codtyp,'radianceclear'  )) ) then
        toDataPresent = .true.
      end if
    end do HEADER

    if ( .not. toDataPresent ) then
      write(*,*) 'WARNING: WILL NOT RUN obsf_updateMissingObsFlags since no ATMS or AMSUA, CSR, AMSUB, MHS'
      return
    end if

    ! If obs files not split and I am not task 0, then return
    if ( .not.obsf_filesSplit() .and. mmpi_myid /= 0 ) return

    FILELOOP: do fileIndex = 1, obsf_nfiles
      if ( obsf_familyType(fileIndex) /= 'TO' ) cycle FILELOOP
      write(*,*) 'INPUT FILE TO  obsf_updateMissingObsFlags = ', trim( obsf_fileName(fileIndex) )
      call obsf_determineSplitFileType( obsFileType, obsf_fileName(fileIndex) )
      if ( trim(obsFileType) /= 'BURP' ) then
        write(*,*) 'obsFileType = ',obsFileType
        write(*,*) 'obsf_updateMissingObsFlags: WARNING this s/r is currently only compatible with BURP files'
      else
        call brpr_updateMissingObsFlags( trim( obsf_fileName(fileIndex) ) )
      end if
    end do FILELOOP

  end subroutine obsf_updateMissingObsFlags

  !--------------------------------------------------------------------------
  ! obsf_copyObsDirectory
  !--------------------------------------------------------------------------
  subroutine obsf_copyObsDirectory(directoryInOut, direction)
    !
    ! :Purpose: Loop on observation files and copy each to and from
    !           the specified directory
    !
    implicit none
    ! Arguments:
    character(len=*) :: directoryInOut
    character(len=*) :: direction

    ! Locals:
    integer            :: status, fileIndex, baseNameIndexBeg
    character(len=200) :: baseName, fullName

    if (trim(direction) == 'TO') then
      ! Copy files TO the specified directory

      ! Create destination directory
      if (mmpi_myid == 0) status = clib_mkdir_r(trim(directoryInOut))
      if (obsf_filesSplit()) call rpn_comm_barrier('GRID',status)

      ! If obs files not split and I am not task 0, then return
      if ( .not.obsf_filesSplit() .and. mmpi_myid /= 0 ) return

      do fileIndex = 1, obsf_nfiles
        fullName = trim( obsf_fileName(fileIndex) )
        baseNameIndexBeg = index(fullName,'/',back=.true.)
        baseName = fullName(baseNameIndexBeg+1:)
        write(*,*) 'obsf_copyObsDirectory: Copying file ', trim(baseName)
        status = utl_copyFile(fullName,trim(directoryInOut)//'/'//trim(baseName))
      end do

    else if (trim(direction) == 'FROM') then
      ! Copy files FROM the specified directory

      ! If obs files not split and I am not task 0, then return
      if ( .not.obsf_filesSplit() .and. mmpi_myid /= 0 ) return

      do fileIndex = 1, obsf_nfiles
        fullName = trim( obsf_fileName(fileIndex) )
        baseNameIndexBeg = index(fullName,'/',back=.true.)
        baseName = fullName(baseNameIndexBeg+1:)
        write(*,*) 'obsf_copyObsDirectory: Copying file ', trim(baseName)
        status = utl_copyFile(trim(directoryInOut)//'/'//trim(baseName),fullName)
        status = clib_remove(trim(directoryInOut)//'/'//trim(baseName))
      end do

      ! Remove the directory
      if (obsf_filesSplit()) call rpn_comm_barrier('GRID',status)
      if (mmpi_myid == 0) status = clib_remove(trim(directoryInOut))
      
    else

      call utl_abort('obsf_copyObsDirectory: invalid value for direction')

    end if

  end subroutine obsf_copyObsDirectory

  !--------------------------------------------------------------------------
  ! setHeadBodyPrimaryKeyColumns
  !--------------------------------------------------------------------------
  subroutine setHeadBodyPrimaryKeyColumns(obsDat, numHeaderRead, numBodyRead)
    !
    ! :Purpose: Set header/body primary keys in obsSpaceData that 
    !           will ensure unique values over all mpi tasks.
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsDat
    integer, intent(in) :: numHeaderRead
    integer, intent(in) :: numBodyRead

    ! locals:
    integer    :: initialHeaderindex, initialBodyindex
    integer    :: numHeaders, numBodies
    integer    :: headerIndex, bodyIndex
    integer    :: headerIndexBegin, headerIndexEnd, bodyIndexBegin, bodyIndexEnd
    integer    :: ierr
    integer(8) :: headerPrimaryKey, bodyPrimaryKey
    integer, allocatable :: allNumHeaderRead(:), allNumBodyRead(:)

    write(*,*) 'setHeadBodyPrimaryKeyColumns: start'
    numHeaders = obs_numHeader(obsSpaceData)
    numBodies = obs_numBody(obsSpaceData)
     
    write(*,*) 'setHeadBodyPrimaryKeyColumns: numHeaders=', numHeaders, ', numBodies=', numBodies
    write(*,*) 'setHeadBodyPrimaryKeyColumns: numHeaderRead=', numHeaderRead, &
                ', numBodyRead=', numBodyRead

    allocate(allNumHeaderRead(mmpi_nprocs))
    allocate(allNumBodyRead(mmpi_nprocs))
    call rpn_comm_allgather(numHeaderRead,1,'mpi_integer',       &
                            allNumHeaderRead,1,'mpi_integer','GRID',ierr)
    call rpn_comm_allgather(numBodyRead,1,'mpi_integer',       &
                            allNumBodyRead,1,'mpi_integer','GRID',ierr)
    if (mmpi_myid > 0) then
      initialHeaderindex = sum(allNumHeaderRead(1:mmpi_myid))
      initialBodyindex = sum(allNumBodyRead(1:mmpi_myid))
    else
      initialHeaderindex = 0
      initialBodyindex = 0
    end if
    deallocate(allNumHeaderRead)
    deallocate(allNumBodyRead)

    write(*,*) 'setHeadBodyPrimaryKeyColumns: initialHeaderIndex=', initialHeaderindex , &
                ', initialBodyindex=', initialBodyindex

    headerIndexBegin = numHeaders - numHeaderRead + 1
    headerIndexEnd = numHeaders
    headerPrimaryKey = initialHeaderindex
    do headerIndex = headerIndexBegin, headerIndexEnd
      headerPrimaryKey = headerPrimaryKey + 1
      call obs_setHeadPrimaryKey(obsdat, headerIndex, headerPrimaryKey)
    end do
        
    bodyIndexBegin = numBodies - numBodyRead + 1
    bodyIndexEnd = numBodies
    bodyPrimaryKey = initialBodyindex
    do bodyIndex = bodyIndexBegin, bodyIndexEnd
      bodyPrimaryKey = bodyPrimaryKey + 1
      call obs_setBodyPrimaryKey(obsdat, bodyIndex, bodyPrimaryKey)
    end do

    write(*,*) 'setHeadBodyPrimaryKeyColumns: end'

  end subroutine setHeadBodyPrimaryKeyColumns

end module obsFiles_mod
