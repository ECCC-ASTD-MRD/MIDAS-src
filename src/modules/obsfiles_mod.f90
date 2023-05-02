
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
  public :: obsf_nfiles, obsf_cfilnam

  ! public procedures
  public :: obsf_setup, obsf_filesSplit, obsf_determineFileType, obsf_determineSplitFileType
  public :: obsf_readFiles, obsf_writeFiles, obsf_obsSub_read, obsf_obsSub_update
  public :: obsf_addCloudParametersAndEmissivity, obsf_getFileName, obsf_copyObsDirectory
  public :: obsf_updateMissingObsFlags, obsf_cleanObsFiles
  logical           :: obsFilesSplit
  logical           :: initialized = .false.

  integer, parameter :: jpfiles = 150
  integer, parameter :: maxLengthFilename = 1060
  integer, parameter :: familyTypeLen = 2
  integer :: obsf_nfiles, obsf_numMpiUniqueList
  character(len=maxLengthFilename) :: obsf_cfilnam(jpfiles)
  character(len=familyTypeLen)     :: obsf_cfamtyp(jpfiles)
  character(len=48)  :: obsFileMode
  character(len=maxLengthFilename) :: obsf_baseFileNameMpiUniqueList(jpfiles)
  character(len=familyTypeLen)     :: obsf_familyTypeMpiUniqueList(jpfiles)
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
      call brpf_getDateStamp( dateStamp_out, obsf_cfilnam(1) )
    else if ( obsFileType == 'OBSDB' ) then
      call odbf_getDateStamp( dateStamp_out, obsf_cfilnam(1) )
    else if ( obsFileType == 'SQLITE' ) then
      call sqlf_getDateStamp( dateStamp_out, obsf_cfilnam(1) )
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
          if (trim(obsf_cfilnam(fileIndexMpiLocal)) == fileNamefull) exit
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

          call getNumHeadersBodies(obsSpaceData, numHeaderBefore, numBodyBefore)
          call brpf_readFile( obsSpaceData, fileNameFull, obsFamilyType, fileIndexMpiLocal )
          call getNumHeadersBodies(obsSpaceData, numHeaders, numBodies)
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

      if (obsFileType == 'BURP') call brpr_setHeadBodyPrimaryKeyColumns(obsSpaceData, obsFamilyType)
      
    end do

    ! abort if NAMTOV does not exist but there are radiance observation files
    if ( .not. utl_isNamelistPresent('NAMTOV','./flnml') .and. &
        any(obsf_cfamtyp(:) == 'TO') ) then
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
        call obsf_determineSplitFileType( obsFileType, obsf_cfilnam(fileIndex) )

        if ( obsFileType == 'BURP'   ) then
          call brpf_updateFile( obsSpaceData, obsf_cfilnam(fileIndex), obsf_cfamtyp(fileIndex), &
                                fileIndex )
        else if ( obsFileType == 'OBSDB' ) then
          call odbf_updateFile( obsSpaceData, obsf_cfilnam(fileIndex), &
                                obsf_cfamtyp(fileIndex), fileIndex )
        else if ( obsFileType == 'SQLITE' ) then
          call sqlf_updateFile( obsSpaceData, obsf_cfilnam(fileIndex), &
                                obsf_cfamtyp(fileIndex), fileIndex )
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
          fullName = trim(obsf_cfilnam(fileIndex))
          baseNameIndexBeg = index(fullName,'/',back=.true.)
          baseName = fullName(baseNameIndexBeg+1:)

          call odbf_updateFile(obsSpaceData, trim(fileNameDir)//'obsDB/'//trim(baseName), &
                               obsf_cfamtyp(fileIndex), fileIndex)
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
      call obsf_determineSplitFileType( obsFileType, obsf_cfilnam(fileIndex) )

      if ( obsFileType == 'BURP' ) then
        call brpr_burpClean( obsf_cfilnam(fileIndex), obsf_cfamtyp(fileIndex) )
      else if ( obsFileType == 'OBSDB' ) then
        call obdf_clean( obsf_cfilnam(fileIndex), obsf_cfamtyp(fileIndex) )
      else if ( obsFileType == 'SQLITE' ) then
        call sqlf_cleanFile( obsf_cfilnam(fileIndex), obsf_cfamtyp(fileIndex) )
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
    character(len=20) :: clvalu(jpfiles)
    character(len=2)    :: cfami(jpfiles)
    character(len=4)    :: cmyidx, cmyidy
    character(len=256)  :: obsDirectory
    character(len=maxLengthFilename) :: fileName, baseFileNameNoMyId  ! the length should be more than 
                                                                      ! len(obsDirectory)+1+len(clvalu)+1+len(obsf_myIdExt)
    character(len=maxLengthFilename) :: baseFileName(jpfiles)

    character(len=256)               :: fileNamefull
    logical :: fileExists
    integer :: fileIndex

    write(cmyidy,'(I4.4)') (mmpi_myidy + 1)
    write(cmyidx,'(I4.4)') (mmpi_myidx + 1)
    obsf_myIdExt  = trim(cmyidx) // '_' // trim(cmyidy)

    clvalu(:)=''
    ! file names only for burp
    clvalu( 1) = 'brpuan'
    clvalu( 2) = 'brpuas'
    clvalu( 3) = 'brpai'
    clvalu( 4) = 'brpain'
    clvalu( 5) = 'brpais'
    clvalu( 6) = 'brpaie'
    clvalu( 7) = 'brpaiw'
    clvalu( 8) = 'brpsfc'
    clvalu( 9) = 'brpsf'
    clvalu(10) = 'brptov'
    clvalu(11) = 'brpssmis'
    clvalu(12) = 'brpairs'
    clvalu(13) = 'brpto_amsua'
    clvalu(14) = 'brpto_amsua_allsky'
    clvalu(15) = 'brpto_amsub'
    clvalu(16) = 'brpto_amsub_allsky'
    clvalu(17) = 'brpcsr'
    clvalu(18) = 'brpiasi'
    clvalu(19) = 'brpatms'
    clvalu(20) = 'brpatms_allsky'
    clvalu(21) = 'brpcris'
    clvalu(22) = 'brpcrisfsr'
    clvalu(23) = 'brpcrisfsr1'
    clvalu(24) = 'brpcrisfsr2'
    clvalu(25) = 'brpcrisfsr3'
    clvalu(26) = 'brpcrisfsr4'
    clvalu(27) = 'brpsw'
    clvalu(28) = 'brpswgoes9'
    clvalu(29) = 'brpswgoese'
    clvalu(30) = 'brpswgoesw'
    clvalu(31) = 'brpswmodis'
    clvalu(32) = 'brpswmtsate'
    clvalu(33) = 'brpswmtsatw'
    clvalu(34) = 'brpgo'
    clvalu(35) = 'brpsc'
    clvalu(36) = 'brppr'
    clvalu(37) = 'brpro'
    clvalu(38) = 'brphum'
    clvalu(39) = 'brpsat'
    clvalu(40) = 'brpssm'
    clvalu(41) = 'brpgp'
    clvalu(42) = 'brpch'
    clvalu(43) = 'brpua'
    ! new general file names for burp and sqlite
    clvalu(44) = 'obsua'
    clvalu(45) = 'obsuan'
    clvalu(46) = 'obsuas'
    clvalu(47) = 'obsai'
    clvalu(48) = 'obsain'
    clvalu(49) = 'obsais'
    clvalu(50) = 'obsaie'
    clvalu(51) = 'obsaiw'
    clvalu(52) = 'obssfc'
    clvalu(53) = 'obssf'
    clvalu(54) = 'obstov'
    clvalu(55) = 'obsssmis'
    clvalu(56) = 'obsairs'
    clvalu(57) = 'obsto_amsua'
    clvalu(58) = 'obsto_amsua_allsky'
    clvalu(59) = 'obsto_amsub'
    clvalu(60) = 'obsto_amsub_allsky'
    clvalu(61) = 'obscsr'
    clvalu(62) = 'obsiasi'
    clvalu(63) = 'obsatms'
    clvalu(64) = 'obsatms_allsky'
    clvalu(65) = 'obscris'
    clvalu(66) = 'obscrisfsr'
    clvalu(67) = 'obscrisfsr1'
    clvalu(68) = 'obscrisfsr2'
    clvalu(69) = 'obscrisfsr3'
    clvalu(70) = 'obscrisfsr4'
    clvalu(71) = 'obssw'
    clvalu(72) = 'obsswgoes9'
    clvalu(73) = 'obsswgoese'
    clvalu(74) = 'obsswgoesw'
    clvalu(75) = 'obsswmodis'
    clvalu(76) = 'obsswmtsate'
    clvalu(77) = 'obsswmtsatw'
    clvalu(78) = 'obsgo'
    clvalu(79) = 'obssc'
    clvalu(80) = 'obspr'
    clvalu(81) = 'obsro'
    clvalu(82) = 'obshum'
    clvalu(83) = 'obssat'
    clvalu(84) = 'obsssm'
    clvalu(85) = 'obsgp'
    clvalu(86) = 'obsch'
    clvalu(87) = 'obsgl_ssmi'
    clvalu(88) = 'obsgl_ssmis'
    clvalu(89) = 'obsgl_amsr2'
    clvalu(90) = 'obsgl_ascat'
    clvalu(91) = 'obsgl_avhrr'
    clvalu(92) = 'obsgl_cisA'
    clvalu(93) = 'obsgl_cisI'
    clvalu(94) = 'obsgl_cisL'
    clvalu(95) = 'obsgl_cisR'
    clvalu(96) = 'cmaheader' ! file name for CMA format used by EnKF
    clvalu(97) = 'brpsst'
    clvalu(98) = 'obsal'
    clvalu(99) = 'obsradar'
    clvalu(100)= 'obssst_insitu'
    clvalu(101)= 'obshydro'
    clvalu(102)= 'obsmwhs2'
    clvalu(103)= 'brpmwhs2'
    clvalu(104)= 'obssarwinds'
    clvalu(105)= 'obssst_avhrr'
    clvalu(106)= 'obssst_amsr2'
    clvalu(107)= 'obssst_viirs'
    clvalu(108)= 'obssst_pseudo'
    clvalu(109)= 'obssst'

    cfami(:)   = ''
    cfami( 1)  = 'UA'
    cfami( 2)  = 'UA'
    cfami( 3)  = 'AI'
    cfami( 4)  = 'AI'
    cfami( 5)  = 'AI'
    cfami( 6)  = 'AI'
    cfami( 7)  = 'AI'
    cfami( 8)  = 'SF'
    cfami( 9)  = 'SF'
    cfami(10)  = 'TO'
    cfami(11)  = 'TO'
    cfami(12)  = 'TO'
    cfami(13)  = 'TO'
    cfami(14)  = 'TO'
    cfami(15)  = 'TO'
    cfami(16)  = 'TO'
    cfami(17)  = 'TO'
    cfami(18)  = 'TO'
    cfami(19)  = 'TO'
    cfami(20)  = 'TO'
    cfami(21)  = 'TO'
    cfami(22)  = 'TO'
    cfami(23)  = 'TO'
    cfami(24)  = 'TO'
    cfami(25)  = 'TO'
    cfami(26)  = 'TO'
    cfami(27)  = 'SW'
    cfami(28)  = 'SW'
    cfami(29)  = 'SW'
    cfami(30)  = 'SW'
    cfami(31)  = 'SW'
    cfami(32)  = 'SW'
    cfami(33)  = 'SW'
    cfami(34)  = 'GO'
    cfami(35)  = 'SC'
    cfami(36)  = 'PR'
    cfami(37)  = 'RO'
    cfami(38)  = 'HU'
    cfami(39)  = 'ST'
    cfami(40)  = 'MI'
    cfami(41)  = 'GP'
    cfami(42)  = 'CH'
    cfami(43)  = 'UA'
    cfami(44)  = 'UA'
    cfami(45)  = 'UA'
    cfami(46)  = 'UA'
    cfami(47)  = 'AI'
    cfami(48)  = 'AI'
    cfami(49)  = 'AI'
    cfami(50)  = 'AI'
    cfami(51)  = 'AI'
    cfami(52)  = 'SF'
    cfami(53)  = 'SF'
    cfami(54)  = 'TO'
    cfami(55)  = 'TO'
    cfami(56)  = 'TO'
    cfami(57)  = 'TO'
    cfami(58)  = 'TO'
    cfami(59)  = 'TO'
    cfami(60)  = 'TO'
    cfami(61)  = 'TO'
    cfami(62)  = 'TO'
    cfami(63)  = 'TO'
    cfami(64)  = 'TO'
    cfami(65)  = 'TO'
    cfami(66)  = 'TO'
    cfami(67)  = 'TO'
    cfami(68)  = 'TO'
    cfami(69)  = 'TO'
    cfami(70)  = 'TO'
    cfami(71)  = 'SW'
    cfami(72)  = 'SW'
    cfami(73)  = 'SW'
    cfami(74)  = 'SW'
    cfami(75)  = 'SW'
    cfami(76)  = 'SW'
    cfami(77)  = 'SW'
    cfami(78)  = 'GO'
    cfami(79)  = 'SC'
    cfami(80)  = 'PR'
    cfami(81)  = 'RO'
    cfami(82)  = 'HU'
    cfami(83)  = 'ST'
    cfami(84)  = 'MI'
    cfami(85)  = 'GP'
    cfami(86)  = 'CH'
    cfami(87)  = 'GL'
    cfami(88)  = 'GL'
    cfami(89)  = 'GL'
    cfami(90)  = 'GL'
    cfami(91)  = 'GL'
    cfami(92)  = 'GL'
    cfami(93)  = 'GL'
    cfami(94)  = 'GL'
    cfami(95)  = 'GL'
    cfami(96)  = 'XX' ! dummy family type for CMA, since it contains all families
    cfami(97)  = 'TM'
    cfami(98)  = 'AL'
    cfami(99)  = 'RA'
    cfami(100) = 'TM'
    cfami(101) = 'HY'
    cfami(102) = 'TO'
    cfami(103) = 'TO'
    cfami(104) = 'SF'
    cfami(105) = 'TM'
    cfami(106) = 'TM'
    cfami(107) = 'TM'
    cfami(108) = 'TM'
    cfami(109) = 'TM'

    obsDirectory = 'obs'

    obsf_nfiles = 0
    obsf_cfilnam(1) = 'DUMMY_FILE_NAME'
    baseFileName(:) = ''

    do fileIndex = 1, jpfiles 

      if(clvalu(fileIndex) == '') exit

      baseFileNameNoMyId = trim(obsDirectory) // '/' // trim(clvalu(fileIndex))
      fileName = trim(baseFileNameNoMyId) // '_' // trim(obsf_myIdExt)
      fileNameFull = ram_fullWorkingPath(fileName,noAbort_opt=.true.)
      inquire(file=trim(fileNameFull),exist=fileExists)

      if (.not. fileExists ) then
        fileName=trim(baseFileNameNoMyId)
        fileNameFull = ram_fullWorkingPath(fileName, noAbort_opt=.true.)
        inquire(file=trim(fileNameFull), exist=fileExists)
      end if

      if ( fileExists ) then
        obsf_nfiles=obsf_nfiles + 1
        baseFileName(obsf_nfiles) = trim(baseFileNameNoMyId)
        obsf_cfilnam(obsf_nfiles) = fileNameFull
        obsf_cfamtyp(obsf_nfiles) = cfami(fileIndex)
      end if

    end do

    call setObsFilesMpiUniqueList(baseFileName)
    
    write(*,*) ' '
    write(*,*)'obsf_setupFileNames: Number of observation files is :', obsf_nfiles
    write(*,*)'Type  Name '
    write(*,*)'----  ---- '
    do fileIndex = 1, obsf_nfiles
      write(*,'(1X,A2,1X,A60)' ) obsf_cfamtyp(fileIndex), trim(obsf_cfilnam(fileIndex))
    end do

  end subroutine obsf_setupFileNames

  !--------------------------------------------------------------------------
  ! setObsFilesMpiUniqueList
  !--------------------------------------------------------------------------
  subroutine setObsFilesMpiUniqueList(baseFileName)
    !
    ! :Purpose: Create a unique list of obs filenames/familyTypes across all mpi tasks.
    !
    implicit none

    ! Arguments
    character(len=*), intent(in) :: baseFileName(:)

    ! Locals:
    integer :: fileIndex, fileIndex2, procIndex, ierr
    character(len=maxLengthFilename), allocatable :: baseFileNameAllMpi(:,:)
    character(len=familyTypeLen), allocatable :: familyTypeAllMpi(:,:)

    ! Communicate filenames across all mpi tasks
    allocate(baseFileNameAllMpi(jpfiles,mmpi_nprocs))
    baseFileNameAllMpi(:,:) = ''
    call mmpi_allgather_string(baseFileName, baseFileNameAllMpi, &
                                jpfiles, maxLengthFilename, mmpi_nprocs, &
                                "GRID", ierr)

    ! Communicate familyTypes across all mpi tasks
    allocate(familyTypeAllMpi(jpfiles,mmpi_nprocs))
    familyTypeAllMpi(:,:) = ''
    call mmpi_allgather_string(obsf_cfamtyp, familyTypeAllMpi, &
                                jpfiles, familyTypeLen, mmpi_nprocs, &
                                "GRID", ierr)

    ! Create a unique list of obs filenames/familytype across all mpi tasks without duplicates
    obsf_baseFileNameMpiUniqueList(:) = ''
    obsf_familyTypeMpiUniqueList(:) = ''
    obsf_numMpiUniqueList = 1
    obsf_baseFileNameMpiUniqueList(obsf_numMpiUniqueList) = baseFileNameAllMpi(1,1)
    obsf_familyTypeMpiUniqueList(obsf_numMpiUniqueList) = familyTypeAllMpi(1,1)
    do procIndex = 1, mmpi_nprocs
      loopFilename: do fileIndex = 1, jpfiles 
        if (trim((baseFileNameAllMpi(fileIndex,procIndex))) == '') cycle loopFilename

        ! cycle if filename already exists in the unique list
        do fileIndex2 = 1, obsf_numMpiUniqueList
          if (trim(obsf_baseFileNameMpiUniqueList(fileIndex2)) == &
              trim(baseFileNameAllMpi(fileIndex,procIndex))) then
            cycle loopFilename
          end if
        end do
        
        ! add the filename to the unique list
        obsf_numMpiUniqueList = obsf_numMpiUniqueList + 1
        obsf_baseFileNameMpiUniqueList(obsf_numMpiUniqueList) = baseFileNameAllMpi(fileIndex,procIndex)
        obsf_familyTypeMpiUniqueList(obsf_numMpiUniqueList) = familyTypeAllMpi(fileIndex,procIndex)

      end do loopFilename
    end do

    if (mmpi_myid == 0) then
      write(*,*) 'setObsFilesMpiUniqueList: all mpi familType/filename before creating unique list:' 
      write(*,*)'Type  Name '
      write(*,*)'----  ---- '
      do fileIndex = 1, jpfiles
        if (trim((baseFileNameAllMpi(fileIndex,1))) == '') cycle

        write(*,'(1X,A2,1X,A60)' ) trim(familyTypeAllMpi(fileIndex,1)), &
                                    trim(baseFileNameAllMpi(fileIndex,1))
      end do
    end if

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

    if ( mmpi_myid == procID ) call obsf_determineSplitFileType( obsFileType, obsf_cfilnam(1) )

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
       if (obsfam == obsf_cfamtyp(ifile)) then
          filename = obsf_cfilnam(ifile)
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

      call obsf_determineSplitFileType( obsFileType, obsf_cfilnam(fileIndex) )
      if ( trim(obsFileType) == 'SQLITE' )  then
        call sqlf_addCloudParametersandEmissivity(obsSpaceData, fileIndex, obsf_cfilnam(fileIndex))
      else if ( trim(obsFileType) == 'BURP' ) then 
        call brpr_addCloudParametersandEmissivity(obsSpaceData, fileIndex, trim( obsf_cfilnam(fileIndex) ) )
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
      if ( obsf_cfamtyp(fileIndex) /= 'TO' ) cycle FILELOOP
      write(*,*) 'INPUT FILE TO  obsf_updateMissingObsFlags = ', trim( obsf_cfilnam(fileIndex) )
      call obsf_determineSplitFileType( obsFileType, obsf_cfilnam(fileIndex) )
      if ( trim(obsFileType) /= 'BURP' ) then
        write(*,*) 'obsFileType = ',obsFileType
        write(*,*) 'obsf_updateMissingObsFlags: WARNING this s/r is currently only compatible with BURP files'
      else
        call brpr_updateMissingObsFlags( trim( obsf_cfilnam(fileIndex) ) )
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
        fullName = trim( obsf_cfilnam(fileIndex) )
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
        fullName = trim( obsf_cfilnam(fileIndex) )
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
  ! getNumHeadersBodies
  !--------------------------------------------------------------------------
  subroutine getNumHeadersBodies(obsSpaceData, numHeaders, numBodies)
    !
    ! :Purpose: Get number of local headers/bodies from obsSpaceData.
    !
    implicit none

    ! arguments
    type(struct_obs), intent(in) :: obsSpaceData
    integer, intent(out) :: numHeaders, numBodies

    numHeaders = obs_numHeader(obsSpaceData)
    numBodies = obs_numBody(obsSpaceData)

  end subroutine getNumHeadersBodies

end module obsFiles_mod
