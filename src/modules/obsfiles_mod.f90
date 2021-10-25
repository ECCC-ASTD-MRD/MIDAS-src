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

module obsFiles_mod
  ! MODULE obsFiles_mod (prefix='obsf' category='6. Observation input/output')
  !
  ! :Purpose: High-level module to handle reading/writing of observations that
  !           can be stored in one of several different formats. Currently, the
  !           only supported formats are:
  !
  !              1. BURP
  !              2. CMA (binary format of obsSpaceData contents)
  !              3. SQLITE (burp2rdb format)
  !              4. SQLITE (obsDB format)
  !
  use mpi_mod
  use ramdisk_mod
  use mathPhysConstants_mod
  use earthConstants_mod
  use utilities_mod
  use obsSpaceData_mod
  use burpFiles_mod
  use sqliteFiles_mod
  use obsdbFiles_mod
  use cmaFiles_mod
  use bufr_mod
  use obsSubSpaceData_mod
  use obsUtil_mod
  use obsVariableTransforms_mod
  use burpread_mod
  use biasCorrectionConv_mod
  use clib_interfaces_mod
  use tovs_nl_mod 
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

  integer, parameter :: jpfiles=150
  integer, parameter :: maxLengthFilename=1060
  integer :: obsf_nfiles
  character(len=maxLengthFilename) :: obsf_cfilnam(jpfiles)
  character(len=2)   :: obsf_cfamtyp(jpfiles)

  character(len=48)  :: obsFileMode

contains

  subroutine obsf_setup( dateStamp_out, obsFileMode_in, obsFileType_opt )

    implicit none

    ! arguments
    integer         , intent(out)           :: dateStamp_out
    character(len=*), intent(in)            :: obsFileMode_in
    character(len=*), intent(out), optional :: obsFileType_opt

    ! locals
    character(len=10)                       :: obsFileType

    obsFileMode = trim(obsFileMode_in)

    !
    ! Initialize file names and the file type
    !
    call obsf_setupFileNames()
    call obsf_determineFileType(obsFileType)

    if ( present(obsFileType_opt) ) obsFileType_opt = obsFileType

    !
    ! Determine if obsFiles are split
    !
    if ( obsFileType == 'BURP' .or. obsFileType == 'SQLITE') then
      obsFilesSplit = .true.
    else if ( obsFileType == 'CMA' ) then
      obsFilesSplit = .false.
    else
      call utl_abort('obsf_setup: invalid observation file type: ' // trim(obsFileType))
    end if

    !
    ! Do some setup of observation files
    !
    if ( obsFileType == 'BURP' ) then
      call brpf_getDateStamp( dateStamp_out, obsf_cfilnam(1) )
    else if ( obsFileType == 'SQLITE' ) then
      call sqlf_getDateStamp( dateStamp_out, obsf_cfilnam(1) )
    else
      dateStamp_out = -1
    end if

    initialized = .true.

  end subroutine obsf_setup


  function obsf_filesSplit() result(obsFilesSplit_out)
    implicit none

    ! arguments
    logical :: obsFilesSplit_out

    if ( .not.initialized ) call utl_abort('obsf_filesSplit: obsFiles_mod not initialized!')

    obsFilesSplit_out = obsFilesSplit
        
  end function obsf_filesSplit


  subroutine obsf_readFiles( obsSpaceData )

    implicit none

    ! arguments
    type(struct_obs)  :: obsSpaceData

    ! locals
    integer           :: fileIndex
    character(len=10) :: obsFileType

    if ( .not.initialized ) call utl_abort('obsf_readFiles: obsFiles_mod not initialized!')

    call obsf_determineFileType(obsFileType)

    if ( obsFileType == 'CMA' ) then

      ! read same global CMA file on all mpi tasks
      call cma_readFiles(obsSpaceData)

    else 

      ! for every splitted file, the file type is defined separately 
      do fileIndex = 1, obsf_nfiles

        call obsf_determineSplitFileType( obsFileType, obsf_cfilnam(fileIndex) )
        if ( obsFileType == 'BURP' )   then
          ! Add extra bias correction element to GP and TO files
          if ( bcc_biasActive( obsf_cfamtyp(fileIndex) ) .or. ( obsf_cfamtyp(fileIndex) == 'TO' ) ) &
               call brpr_addElementsToBurp(obsf_cfilnam(fileIndex),  obsf_cfamtyp(fileIndex), beSilent_opt=.false.)
          call brpf_readFile( obsSpaceData, obsf_cfilnam(fileIndex), obsf_cfamtyp(fileIndex), fileIndex )
        end if
        if ( obsFileType == 'SQLITE' ) then
          if (odbf_isActive()) then
            call odbf_readFile( obsSpaceData, obsf_cfilnam(fileIndex), obsf_cfamtyp(fileIndex), fileIndex )
          else
            call sqlf_readFile( obsSpaceData, obsf_cfilnam(fileIndex), obsf_cfamtyp(fileIndex), fileIndex )
          end if
        end if

      end do

    end if

    ! abort if NAMTOV does not exist but there are radiance observation files
    if ( .not. utl_isNamelistPresent('NAMTOV','./flnml') .and. &
        any(obsf_cfamtyp(:) == 'TO') ) then
      call utl_abort('obsf_readFiles: Namelist block NAMTOV is missing but there are radiance observation files')
    end if

    ! initialize OBS_HIND for each observation
    call obs_sethind(obsSpaceData)
  end subroutine obsf_readFiles


  subroutine obsf_writeFiles( obsSpaceData, HXens_mpiglobal_opt, asciDumpObs_opt, writeDiagFiles_opt)
  implicit none

  ! arguments
  type(struct_obs)  :: obsSpaceData
  real(8), optional :: HXens_mpiglobal_opt(:,:)
  logical, optional :: asciDumpObs_opt
  logical, optional :: writeDiagFiles_opt

  ! locals
  integer           :: fileIndex, fnom, fclos, nulnam, ierr
  character(len=10) :: obsFileType, sfFileName
  character(len=*), parameter :: myName = 'obsf_writeFiles'
  character(len=*), parameter :: myWarning = myName //' WARNING: '

  ! namelist variables
  logical :: lwritediagsql
  logical :: onlyAssimObs
  logical :: addFSOdiag

  namelist /namwritediag/lwritediagsql,onlyAssimObs,addFSOdiag

  if ( .not.initialized ) call utl_abort('obsf_writeFiles: obsFiles_mod not initialized!')
 
  call obsf_determineFileType(obsFileType)

  nulnam=0
  lwritediagsql = .false.
  onlyAssimObs = .false.
  addFSOdiag = .false.
  ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read(nulnam,nml=namwritediag,iostat=ierr)
  if (ierr /= 0) write(*,*) myWarning//' namwritediag is missing in the namelist. The default value will be taken.'
  if (mpi_myid == 0) write(*,nml = namwritediag)
  if (present(writeDiagFiles_opt)) then
    lwritediagsql = lwritediagsql .and. writeDiagFiles_opt
  end if
  ierr=fclos(nulnam)

  if ( obsFileType == 'BURP' .or. obsFileType == 'SQLITE' ) then

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
      else if ( obsFileType == 'SQLITE' ) then
        if (odbf_isActive()) then
          call odbf_updateFile( obsSpaceData, obsf_cfilnam(fileIndex), &
                                obsf_cfamtyp(fileIndex), fileIndex )
        else
          call sqlf_updateFile( obsSpaceData, obsf_cfilnam(fileIndex), &
                                obsf_cfamtyp(fileIndex), fileIndex )
        end if
      end if
    end do

    if ( present(HXens_mpiglobal_opt) .and. mpi_myid == 0 ) then
      call obsf_writeHX(obsSpaceData, HXens_mpiglobal_opt)
    end if

  else if ( obsFileType == 'CMA' ) then

    ! only 1 mpi task should do the writing
    if( mpi_myid == 0 ) call cma_writeFiles( obsSpaceData, HXens_mpiglobal_opt )

  end if

  if ( index(obsf_getFileName('SF'), 'sfc') > 0 ) then
    sfFileName = 'sfc'
  else
    sfFileName = 'sf'
  end if

  if (lwritediagsql) call sqlf_writeSqlDiagFiles( obsSpaceData, sfFileName, onlyAssimObs, addFSOdiag )

  if ( present(asciDumpObs_opt) ) then
    if ( asciDumpObs_opt ) then
      if ( obsFileType == 'BURP' .or. obsFileType == 'SQLITE' .or. mpi_myid == 0   ) then
        ! all processors write to files only for BURP and SQLITE    
        call obsf_writeAsciDump(obsSpaceData)
      end if
    end if

  end if
  end subroutine obsf_writeFiles


  subroutine obsf_cleanObsFiles()
    implicit none

    ! locals
    integer           :: fileIndex
    character(len=10) :: obsFileType

    if ( .not.initialized ) call utl_abort('obsf_cleanObsFiles: obsFiles_mod not initialized!')
   
    call obsf_determineFileType(obsFileType)

    if ( obsFileType /= 'BURP' .and. obsFileType /= 'SQLITE' ) then
      write(*,*) 'obsf_cleanObsFiles: obsFileType=', obsFileType, &
                ' is not BURP nor SQLITE. Return.' 
      return
    end if

    do fileIndex = 1, obsf_nfiles
      call obsf_determineSplitFileType( obsFileType, obsf_cfilnam(fileIndex) )

      if ( obsFileType == 'BURP' ) then
        call brpr_burpClean( obsf_cfilnam(fileIndex), obsf_cfamtyp(fileIndex) )
      else if ( obsFileType == 'SQLITE' ) then
        call sqlf_cleanFile( obsf_cfilnam(fileIndex), obsf_cfamtyp(fileIndex) )
      end if
    end do

  end subroutine obsf_cleanObsFiles 


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
      write(cmyidy,'(I4.4)') (mpi_myidy + 1)
      write(cmyidx,'(I4.4)') (mpi_myidx + 1)
      cmyid  = trim(cmyidx) // '_' // trim(cmyidy)
      fileNameAsciDump = 'obsout_asci_' // trim(cmyid)
    else
      ! only task 0 write the global data
      if ( mpi_myid > 0 ) return
      fileNameAsciDump = 'obsout_asci'
    end if

    write(*,*) 'obsf_writeAsciDump: writing to file : ', fileNameAsciDump
    unitAsciDump = 0
    ierr = fnom(unitAsciDump, fileNameAsciDump, 'FTN+SEQ+FMT+R/W', 0)
    call obs_print(obsSpaceData,unitAsciDump)
    ierr = fclos(unitAsciDump)

  end subroutine obsf_writeAsciDump


  subroutine obsf_setupFileNames()

    implicit none

    ! locals
    character(len=20)   :: clvalu(jpfiles)
    character(len=2)    :: cfami(jpfiles)
    character(len=4)    :: cmyidx, cmyidy
    character(len=9)    :: cmyid
    character(len=256)  :: obsDirectory
    character(len=maxLengthFilename) :: fileName   !! the length should be more than len(obsDirectory)+1+len(clvalu)+1+len(cmyid)
    character(len=256)               :: fileNamefull
    logical :: fileExists
    integer :: fileIndex

    write(cmyidy,'(I4.4)') (mpi_myidy + 1)
    write(cmyidx,'(I4.4)') (mpi_myidx + 1)
    cmyid  = trim(cmyidx) // '_' // trim(cmyidy)

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
    clvalu(16) = 'brpcsr'
    clvalu(17) = 'brpiasi'
    clvalu(18) = 'brpatms'
    clvalu(19) = 'brpatms_allsky'
    clvalu(20) = 'brpcris'
    clvalu(21) = 'brpcrisfsr'
    clvalu(22) = 'brpcrisfsr1'
    clvalu(23) = 'brpcrisfsr2'
    clvalu(24) = 'brpcrisfsr3'
    clvalu(25) = 'brpcrisfsr4'
    clvalu(26) = 'brpsw'
    clvalu(27) = 'brpswgoes9'
    clvalu(28) = 'brpswgoese'
    clvalu(29) = 'brpswgoesw'
    clvalu(30) = 'brpswmodis'
    clvalu(31) = 'brpswmtsate'
    clvalu(32) = 'brpswmtsatw'
    clvalu(33) = 'brpgo'
    clvalu(34) = 'brpsc'
    clvalu(35) = 'brppr'
    clvalu(36) = 'brpro'
    clvalu(37) = 'brphum'
    clvalu(38) = 'brpsat'
    clvalu(39) = 'brpssm'
    clvalu(40) = 'brpgp'
    clvalu(41) = 'brpch'
    clvalu(42) = 'brpua'
    ! new general file names for burp and sqlite
    clvalu(43) = 'obsua'
    clvalu(44) = 'obsuan'
    clvalu(45) = 'obsuas'
    clvalu(46) = 'obsai'
    clvalu(47) = 'obsain'
    clvalu(48) = 'obsais'
    clvalu(49) = 'obsaie'
    clvalu(50) = 'obsaiw'
    clvalu(51) = 'obssfc'
    clvalu(52) = 'obssf'
    clvalu(53) = 'obstov'
    clvalu(54) = 'obsssmis'
    clvalu(55) = 'obsairs'
    clvalu(56) = 'obsto_amsua'
    clvalu(57) = 'obsto_amsua_allsky'
    clvalu(58) = 'obsto_amsub'
    clvalu(59) = 'obscsr'
    clvalu(60) = 'obsiasi'
    clvalu(61) = 'obsatms'
    clvalu(62) = 'obsatms_allsky'
    clvalu(63) = 'obscris'
    clvalu(64) = 'obscrisfsr'
    clvalu(65) = 'obscrisfsr1'
    clvalu(66) = 'obscrisfsr2'
    clvalu(67) = 'obscrisfsr3'
    clvalu(68) = 'obscrisfsr4'
    clvalu(69) = 'obssw'
    clvalu(70) = 'obsswgoes9'
    clvalu(71) = 'obsswgoese'
    clvalu(72) = 'obsswgoesw'
    clvalu(73) = 'obsswmodis'
    clvalu(74) = 'obsswmtsate'
    clvalu(75) = 'obsswmtsatw'
    clvalu(76) = 'obsgo'
    clvalu(77) = 'obssc'
    clvalu(78) = 'obspr'
    clvalu(79) = 'obsro'
    clvalu(80) = 'obshum'
    clvalu(81) = 'obssat'
    clvalu(82) = 'obsssm'
    clvalu(83) = 'obsgp'
    clvalu(84) = 'obsch'
    clvalu(85) = 'obsgl_ssmi'
    clvalu(86) = 'obsgl_ssmis'
    clvalu(87) = 'obsgl_amsr2'
    clvalu(88) = 'obsgl_ascat'
    clvalu(89) = 'obsgl_avhrr'
    clvalu(90) = 'obsgl_cisA'
    clvalu(91) = 'obsgl_cisI'
    clvalu(92) = 'obsgl_cisL'
    clvalu(93) = 'obsgl_cisR'
    clvalu(94) = 'cmaheader' ! file name for CMA format used by EnKF
    clvalu(95) = 'brpsst'
    clvalu(96) = 'obsal'
    clvalu(97) = 'obsradar'
    clvalu(98) = 'obssst_insitu'
    clvalu(99) = 'obshydro'
    clvalu(100)= 'obsmwhs2'
    clvalu(101)= 'brpmwhs2'
    clvalu(102)= 'obssarwinds'
    clvalu(103)= 'obssst_avhrr'
    clvalu(104)= 'obssst_amsr2'
    clvalu(105)= 'obssst_viirs'
    clvalu(106)= 'obssst_pseudo'
    clvalu(107)= 'obssst'

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
    cfami(26)  = 'SW'
    cfami(27)  = 'SW'
    cfami(28)  = 'SW'
    cfami(29)  = 'SW'
    cfami(30)  = 'SW'
    cfami(31)  = 'SW'
    cfami(32)  = 'SW'
    cfami(33)  = 'GO'
    cfami(34)  = 'SC'
    cfami(35)  = 'PR'
    cfami(36)  = 'RO'
    cfami(37)  = 'HU'
    cfami(38)  = 'ST'
    cfami(39)  = 'MI'
    cfami(40)  = 'GP'
    cfami(41)  = 'CH'
    cfami(42)  = 'UA'
    cfami(43)  = 'UA'
    cfami(44)  = 'UA'
    cfami(45)  = 'UA'
    cfami(46)  = 'AI'
    cfami(47)  = 'AI'
    cfami(48)  = 'AI'
    cfami(49)  = 'AI'
    cfami(50)  = 'AI'
    cfami(51)  = 'SF'
    cfami(52)  = 'SF'
    cfami(53)  = 'TO'
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
    cfami(69)  = 'SW'
    cfami(70)  = 'SW'
    cfami(71)  = 'SW'
    cfami(72)  = 'SW'
    cfami(73)  = 'SW'
    cfami(74)  = 'SW'
    cfami(75)  = 'SW'
    cfami(76)  = 'GO'
    cfami(77)  = 'SC'
    cfami(78)  = 'PR'
    cfami(79)  = 'RO'
    cfami(80)  = 'HU'
    cfami(81)  = 'ST'
    cfami(82)  = 'MI'
    cfami(83)  = 'GP'
    cfami(84)  = 'CH'
    cfami(85)  = 'GL'
    cfami(86)  = 'GL'
    cfami(87)  = 'GL'
    cfami(88)  = 'GL'
    cfami(89)  = 'GL'
    cfami(90)  = 'GL'
    cfami(91)  = 'GL'
    cfami(92)  = 'GL'
    cfami(93)  = 'GL'
    cfami(94)  = 'XX' ! dummy family type for CMA, since it contains all families
    cfami(95)  = 'TM'
    cfami(96)  = 'AL'
    cfami(97)  = 'RA'
    cfami(98)  = 'TM'
    cfami(99)  = 'HY'
    cfami(100) = 'TO'
    cfami(101) = 'TO'
    cfami(102) = 'SF'
    cfami(103) = 'TM'
    cfami(104) = 'TM'
    cfami(105) = 'TM'
    cfami(106) = 'TM'
    cfami(107) = 'TM'

    obsDirectory = 'obs'

    obsf_nfiles = 0
    obsf_cfilnam(1) = 'DUMMY_FILE_NAME'

    do fileIndex = 1, jpfiles 

      if(clvalu(fileIndex) == '') exit
      fileName = trim(obsDirectory) // '/' // trim(clvalu(fileIndex)) // '_' // trim(cmyid)
      fileNameFull = ram_fullWorkingPath(fileName,noAbort_opt=.true.)

      inquire(file=trim(fileNameFull),exist=fileExists)
      if (.not. fileExists ) then
        fileName=trim(obsDirectory)//'/'//trim(clvalu(fileIndex))
        fileNameFull = ram_fullWorkingPath(fileName, noAbort_opt=.true.)
        inquire(file=trim(fileNameFull), exist=fileExists)
      end if

      if ( fileExists ) then
        obsf_nfiles=obsf_nfiles + 1
        obsf_cfilnam(obsf_nfiles) = fileNameFull
        obsf_cfamtyp(obsf_nfiles) = cfami(fileIndex)
      end if

    end do

    write(*,*) ' '
    write(*,*)'obsf_setupFileNames: Number of observation files is :', obsf_nfiles
    write(*,*)'Type  Name '
    write(*,*)'----  ---- '
    do fileIndex = 1, obsf_nfiles
      write(*,'(1X,A2,1X,A60)' ) obsf_cfamtyp(fileIndex), trim(obsf_cfilnam(fileIndex))
    end do

  end subroutine obsf_setupFileNames


  subroutine obsf_determineFileType( obsFileType )

    implicit none

    ! arguments
    character(len=*)                       :: obsFileType

    ! locals
    integer :: ierr, procID, all_nfiles(0:(mpi_nprocs-1))
    logical :: fileExists

    call rpn_comm_allgather( obsf_nfiles, 1, 'MPI_INTEGER', &
                             all_nfiles,  1, 'MPI_INTEGER', 'GRID', ierr )
    fileExists = .false.
    procid_loop: do procID = 0, (mpi_nprocs-1)
      if ( all_nfiles(procID) > 0 ) then
        fileExists = .true.
        exit procid_loop
      end if
    end do procid_loop

    if ( .not.fileExists ) then
      call utl_abort('obsf_determineFileType: No observation files found')
    end if

    write(*,*) 'obsf_determineFileType: read obs file that exists on mpi task id: ', procID

    if ( mpi_myid == procID ) call obsf_determineSplitFileType( obsFileType, obsf_cfilnam(1) )

    call rpn_comm_bcastc(obsFileType , len(obsFileType), 'MPI_CHARACTER', procID, 'GRID', ierr)
    write(*,*) 'obsf_determineFileType: obsFileType = ', obsFileType

  end subroutine obsf_determineFileType

  subroutine obsf_determineSplitFileType( obsFileType, fileName )

    implicit none

    ! arguments
    character(len=*),intent(out) :: obsFileType
    character(len=*),intent(in)  :: fileName

    ! locals
    integer           :: ierr, unitFile
    integer           :: fnom, fclos, wkoffit
    character(len=20) :: fileStart

    write(*,*) 'obsf_determineSplitFileType: read obs file: ', trim(fileName)
   
    if ( index(trim(fileName), 'cma' ) > 0 ) then
      obsFileType = 'CMA'
    else

      ierr = wkoffit(trim(fileName))

      if (ierr.eq.6) then
         obsFiletype = 'BURP'
      else

         unitFile = 0
         ierr = fnom(unitFile, fileName, 'FTN+SEQ+FMT+R/O', 0)
         read(unitFile,'(A)') fileStart
         ierr = fclos(unitFile)

         if ( index( fileStart, 'SQLite format 3' ) > 0 ) then
            obsFileType = 'SQLITE'
         else
            call utl_abort('obsf_determineSplitFileType: unknown obs file type')
         end if

      end if
    end if

    write(*,*) 'obsf_determineSplitFileType: obsFileType = ', obsFileType

  end subroutine obsf_determineSplitFileType


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
      write(*,*) "obsf_getFileName: File not found for obs family '" // trim(obsfam) // "'"
    end if
    if (numFound > 1) then
      write(*,*) "obsf_getFileName: WARNING: Multiple files found for obs family '" // trim(obsfam) // "'"
    end if

    if (present(fileFound_opt)) fileFound_opt = (numFound > 0)

  end function obsf_getFileName


  function obsf_obsSub_read( obsfam, stnid, varno, nlev, ndim, bkstp_opt, block_opt, match_nlev_opt, &
                              codtyp_opt ) result(obsdata)
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
    !                             averagine kernels) 
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
               call utl_abort("obsf_obsSub_read: optional varaible 'block_opt' is required for BURP observational files.")
          obsdata = brpf_obsSub_read(filename,stnid,varno,nlev,ndim,block_opt,bkstp_opt=bkstp_opt, &
                                   match_nlev_opt=match_nlev_opt,codtyp_opt=codtyp_opt)
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
       if (obsf_filesSplit() .or. mpi_myid == 0) then
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
    if ( .not.obsf_filesSplit() .and. mpi_myid /= 0 ) return

    do fileIndex = 1, obsf_nfiles

      call obsf_determineSplitFileType( obsFileType, obsf_cfilnam(fileIndex) )
      if ( trim(obsFileType) == 'SQLITE' ) then
        call sqlf_addCloudParametersandEmissivity(obsSpaceData, fileIndex, obsf_cfilnam(fileIndex))
      else if ( trim(obsFileType) == 'BURP' ) then 
        call brpr_addCloudParametersandEmissivity(obsSpaceData, fileIndex, trim( obsf_cfilnam(fileIndex) ) )
      else  
        write(*,*) ' UNKNOWN FileType=',obsFileType
        call utl_abort("obsf_addCloudParametersAndEmissivity: Only BURP or SQLITE observational files supported.")
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
    if ( .not.obsf_filesSplit() .and. mpi_myid /= 0 ) return

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
      if (mpi_myid == 0) status = clib_mkdir_r(trim(directoryInOut))
      if (obsf_filesSplit()) call rpn_comm_barrier('GRID',status)

      ! If obs files not split and I am not task 0, then return
      if ( .not.obsf_filesSplit() .and. mpi_myid /= 0 ) return

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
      if ( .not.obsf_filesSplit() .and. mpi_myid /= 0 ) return

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
      if (mpi_myid == 0) status = clib_remove(trim(directoryInOut))
      
    else

      call utl_abort('obsf_copyObsDirectory: invalid value for direction')

    end if

  end subroutine obsf_copyObsDirectory

end module obsFiles_mod
