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
!! MODULE obsFiles_mod (prefix="obsf" category='6. Observation input/output')
!!
!! *Purpose*: High-level module to handle reading/writing of observations that
!!            can be stored in one of several different formats. Currently, the
!!            only supported formats are:
!!
!!            1. BURP
!!            2. CMA (binary format of obsSpaceData contents)
!!            3. SQLITE
!!
!--------------------------------------------------------------------------
module obsFiles_mod
  use mpi_mod
  use ramdisk_mod
  use mathPhysConstants_mod
  use earthConstants_mod
  use utilities_mod
  use obsSpaceData_mod
  use burpFiles_mod
  use sqliteFiles_mod
  use cmaFiles_mod
  use bufr_mod
  use obsSubSpaceData_mod
  use obsUtil_mod

  implicit none
  save
  private

  ! public variables
  public :: obsf_nfiles, obsf_cfilnam

  ! public procedures
  public :: obsf_setup, obsf_filesSplit, obsf_determineFileType, obsf_determineSplitFileType
  public :: obsf_readFiles, obsf_writeFiles, obsf_obsSub_read, obsf_obsSub_update
  public :: obsf_thinFiles

  logical           :: obsFilesSplit
  logical           :: initialized = .false.

  integer, parameter :: jpfiles=100
  integer, parameter :: maxLengthFilename=1060
  integer :: obsf_nfiles
  character(len=maxLengthFilename) :: obsf_cfilnam(jpfiles)
  character(len=2)   :: obsf_cfamtyp(jpfiles)

  character(len=48)  :: obsFileMode

contains

  subroutine obsf_setup( dateStamp_out, obsFileMode_in, obsFileType_opt )
    implicit none
    ! arguments
    integer                                 :: dateStamp_out
    character(len=*)                        :: obsFileMode_in
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
        if ( obsFileType == 'BURP' )   call brpf_readFile( obsSpaceData, obsf_cfilnam(fileIndex), obsf_cfamtyp(fileIndex), fileIndex )
        if ( obsFileType == 'SQLITE' ) call sqlf_readFile( obsSpaceData, obsf_cfilnam(fileIndex), obsf_cfamtyp(fileIndex), fileIndex )

      end do

    end if

  end subroutine obsf_readFiles


  subroutine obsf_writeFiles( obsSpaceData, HXensT_mpiglobal_opt, asciDumpObs_opt )
  implicit none
  ! arguments
  type(struct_obs)                       :: obsSpaceData
  real(8),             pointer, optional :: HXensT_mpiglobal_opt(:,:)
  logical,                      optional :: asciDumpObs_opt
  ! locals
  integer           :: fileIndex
  character(len=10) :: obsFileType

  if ( .not.initialized ) call utl_abort('obsf_writeFiles: obsFiles_mod not initialized!')
 
  call obsf_determineFileType(obsFileType)

  if ( obsFileType == 'BURP' .or. obsFileType == 'SQLITE' ) then

    if (trim(obsFileMode) == 'analysis') call obsu_computeDirectionSpeedResiduals(obs_oma,obsSpaceData)
    call obsu_computeDirectionSpeedResiduals(obs_omp,obsSpaceData)
    if (trim(obsFileMode) == 'analysis' .or. trim(obsFileMode) == 'FSO') call obsu_setassflg(obsSpaceData)
    call obsu_updateFlagWindDirectionSpeed(obsSpaceData)
    ! Put the scale factor for FSO
    if (trim(obsFileMode) == 'FSO') call obsu_scaleFSO(obsSpaceData)

    do fileIndex = 1, obsf_nfiles

      call obsf_determineSplitFileType( obsFileType, obsf_cfilnam(fileIndex) )
      if ( obsFileType == 'BURP'   ) call brpf_updateFile( obsSpaceData, obsf_cfilnam(fileIndex), obsf_cfamtyp(fileIndex), fileIndex )
      if ( obsFileType == 'SQLITE' ) call sqlf_updateFile( obsSpaceData, obsf_cfilnam(fileIndex), obsf_cfamtyp(fileIndex), fileIndex )

    end do

    if ( present(HXensT_mpiglobal_opt) .and. mpi_myid == 0 ) then
      call obsf_writeHX(obsSpaceData, HXensT_mpiglobal_opt)
    end if

  else if ( obsFileType == 'CMA' ) then

    ! only 1 mpi task should do the writing
    if( mpi_myid == 0 ) call cma_writeFiles( obsSpaceData, HXensT_mpiglobal_opt )

  end if

  if ( present(asciDumpObs_opt) ) then
    if ( asciDumpObs_opt ) then
      if ( obsFileType == 'BURP' .or. obsFileType == 'SQLITE' .or. mpi_myid == 0   ) then
        ! all processors write to files only for BURP and SQLITE    
        call obsf_writeAsciDump(obsSpaceData)
      end if
    end if

  end if

  end subroutine obsf_writeFiles


  !--------------------------------------------------------------------------
  !!
  !! *Purpose*: to reduce the number of observation data
  !!
  !! *Note*:    operates only on SQL files. Issues a warning for other file types
  !!
  !--------------------------------------------------------------------------
  subroutine obsf_thinFiles(obsSpaceData)
    implicit none
    ! arguments
    type(struct_obs), intent(inout) :: obsSpaceData

    ! locals
    integer :: fileIndex
    character(len=10) :: obsFileType

    if(.not.initialized) call utl_abort( &
                                'obsf_thinFiles: obsFiles_mod not initialized!')

    call obsf_determineFileType(obsFileType)
    if ( obsFileType /= 'SQLITE' ) then
      write(*,*)"WARNING:  observation thinning cannot be done with a BURP file."
      write(*,*)"          No observation has been removed from the data base."
      return
    end if

    do fileIndex = 1, obsf_nfiles
      call sqlf_thinFile(obsSpaceData, obsf_cfilnam(fileIndex), &
                         obsf_cfamtyp(fileIndex), fileIndex)

    end do
  end subroutine obsf_thinFiles


  subroutine obsf_writeHX(obsSpaceData, HXensT_mpiglobal)
    implicit none
    ! arguments
    type(struct_obs) :: obsSpaceData
    real(8), pointer :: HXensT_mpiglobal(:,:)
    ! locals
    integer :: unitHX, ierr, headerIndex, fnom, fclos
    character(len=10) :: fileNameHX

    write(*,*) 'obsf_writeHX: Starting'

    fileNameHX     = 'cmahxout'
    unitHX = 0
    ierr = fnom(unitHX, fileNameHX, 'FTN+SEQ+UNF+R/W', 0)

    do headerIndex = 1, obs_numHeader(obsSpaceData)

      call obs_write_hx(obsSpaceData, HXensT_mpiglobal, headerIndex, unitHX)

    enddo
 
    ierr = fclos(unitHX)

  endsubroutine obsf_writeHX


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
    clvalu(14) = 'brpto_amsub'
    clvalu(15) = 'brpcsr'
    clvalu(16) = 'brpiasi'
    clvalu(17) = 'brpatms'
    clvalu(18) = 'brpcris'
    clvalu(19) = 'brpsw'
    clvalu(20) = 'brpswgoes9'
    clvalu(21) = 'brpswgoese'
    clvalu(22) = 'brpswgoesw'
    clvalu(23) = 'brpswmodis'
    clvalu(24) = 'brpswmtsate'
    clvalu(25) = 'brpswmtsatw'
    clvalu(26) = 'brpgo'
    clvalu(27) = 'brpsc'
    clvalu(28) = 'brppr'
    clvalu(29) = 'brpro'
    clvalu(30) = 'brphum'
    clvalu(31) = 'brpsat'
    clvalu(32) = 'brpssm'
    clvalu(33) = 'brpgp'
    clvalu(34) = 'brpch'
    clvalu(35) = 'brpua'
    ! new general file names for burp and sqlite
    clvalu(36) = 'obsua'
    clvalu(37) = 'obsuan'
    clvalu(38) = 'obsuas'
    clvalu(39) = 'obsai'
    clvalu(40) = 'obsain'
    clvalu(41) = 'obsais'
    clvalu(42) = 'obsaie'
    clvalu(43) = 'obsaiw'
    clvalu(44) = 'obssfc'
    clvalu(45) = 'obssf'
    clvalu(46) = 'obstov'
    clvalu(47) = 'obsssmis'
    clvalu(48) = 'obsairs'
    clvalu(49) = 'obsto_amsua'
    clvalu(50) = 'obsto_amsub'
    clvalu(51) = 'obscsr'
    clvalu(52) = 'obsiasi'
    clvalu(53) = 'obsatms'
    clvalu(54) = 'obscris'
    clvalu(55) = 'obssw'
    clvalu(56) = 'obsswgoes9'
    clvalu(57) = 'obsswgoese'
    clvalu(58) = 'obsswgoesw'
    clvalu(59) = 'obsswmodis'
    clvalu(60) = 'obsswmtsate'
    clvalu(61) = 'obsswmtsatw'
    clvalu(62) = 'obsgo'
    clvalu(63) = 'obssc'
    clvalu(64) = 'obspr'
    clvalu(65) = 'obsro'
    clvalu(66) = 'obshum'
    clvalu(67) = 'obssat'
    clvalu(68) = 'obsssm'
    clvalu(69) = 'obsgp'
    clvalu(70) = 'obsch'
    clvalu(71) = 'obsgl_ssmi'
    clvalu(72) = 'obsgl_ssmis'
    clvalu(73) = 'obsgl_amsr2'
    clvalu(74) = 'obsgl_ascat'
    clvalu(75) = 'obsgl_avhrr'
    clvalu(76) = 'obsgl_cisA'
    clvalu(77) = 'obsgl_cisI'
    clvalu(78) = 'obsgl_cisL'
    clvalu(79) = 'obsgl_cisR'
    ! file name for CMA format used by EnKF
    clvalu(80) = 'cmaheader'
    ! Sea Surface Temperature data file name
    clvalu(81) = 'brpsst'
    clvalu(82) = 'obsal'
    clvalu(83) = 'obsradar'

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
    cfami(19)  = 'SW'
    cfami(20)  = 'SW'
    cfami(21)  = 'SW'
    cfami(22)  = 'SW'
    cfami(23)  = 'SW'
    cfami(24)  = 'SW'
    cfami(25)  = 'SW'
    cfami(26)  = 'GO'
    cfami(27)  = 'SC'
    cfami(28)  = 'PR'
    cfami(29)  = 'RO'
    cfami(30)  = 'HU'
    cfami(31)  = 'ST'
    cfami(32)  = 'MI'
    cfami(33)  = 'GP'
    cfami(34)  = 'CH'
    cfami(35)  = 'UA'
    cfami(36)  = 'UA'
    cfami(37)  = 'UA'
    cfami(38)  = 'UA'
    cfami(39)  = 'AI'
    cfami(40)  = 'AI'
    cfami(41)  = 'AI'
    cfami(42)  = 'AI'
    cfami(43)  = 'AI'
    cfami(44)  = 'SF'
    cfami(45)  = 'SF'
    cfami(46)  = 'TO'
    cfami(47)  = 'TO'
    cfami(48)  = 'TO'
    cfami(49)  = 'TO'
    cfami(50)  = 'TO'
    cfami(51)  = 'TO'
    cfami(52)  = 'TO'
    cfami(53)  = 'TO'
    cfami(54)  = 'TO'
    cfami(55)  = 'SW'
    cfami(56)  = 'SW'
    cfami(57)  = 'SW'
    cfami(58)  = 'SW'
    cfami(59)  = 'SW'
    cfami(60)  = 'SW'
    cfami(61)  = 'SW'
    cfami(62)  = 'GO'
    cfami(63)  = 'SC'
    cfami(64)  = 'PR'
    cfami(65)  = 'RO'
    cfami(66)  = 'HU'
    cfami(67)  = 'ST'
    cfami(68)  = 'MI'
    cfami(69)  = 'GP'
    cfami(70)  = 'CH'
    cfami(71)  = 'GL'
    cfami(72)  = 'GL'
    cfami(73)  = 'GL'
    cfami(74)  = 'GL'
    cfami(75)  = 'GL'
    cfami(76)  = 'GL'
    cfami(77)  = 'GL'
    cfami(78)  = 'GL'
    cfami(79)  = 'GL'
    ! dummy family type for CMA, since it contains all families
    cfami(80)  = 'XX'
    cfami(81)  = 'TM'
    cfami(82)  = 'AL'
    cfami(83)  = 'RA'

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
    integer :: ierr, procID, unitFile, all_nfiles(0:(mpi_nprocs-1))
    integer :: fnom, fclos
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
    character(len=*) ,intent(in) :: fileName
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

!--------------------------------------------------------------------------
!! *Purpose*: Returns the observations file name assigned to the calling processor.
!!            If the input family has more than one file, the first file found will
!!            be returned.
!!
!! @author M. Sitwell  Sept 2016
!!
!! Input:
!!v    obsfam            observation family name
!!
!! Output:
!!v    filename       file name of associated observations file
!!v    found_opt      logical indicating if a file could be found for the family (optional)
!--------------------------------------------------------------------------
  function obsf_get_filename(obsfam,found_opt) result(filename)

    implicit none

    character(len=2), intent(in) :: obsfam
    logical, intent(out), optional :: found_opt
    character(len=maxLengthFilename) :: filename
 
    logical :: found
    integer :: ifile

    filename = ""
    found = .false.
       
    do ifile=1,obsf_nfiles
       if (obsfam == obsf_cfamtyp(ifile)) then
          filename = obsf_cfilnam(ifile)
          found = .true.
          exit
       end if
    end do

    if (.not.found) write(*,*) "obsf_get_filename: File not found for observation family '" // trim(obsfam) // "'"

    if (present(found_opt)) found_opt = found

  end function obsf_get_filename

!--------------------------------------------------------------------------
!! *Purpose*: Retrieves information for observations from observation files and returns the data
!!            in a struct_oss_obsdata object. Data will be retrieved for all nodes that have valid
!!            filenames for the specied observational family and combined into one struct_oss_obsdata
!!            if the observational files are split.
!!
!! @author M. Sitwell, ARQI/AQRD, Sept 2016
!!
!! Revisions:
!!v       Y. Rochon, ARQI/AQRD, Nov 2016
!!v         - Added optional input argument codtyp_opt and option of varno <=0 (in brpf_chem_read)
!!
!! Input:
!!v           obsfam          observation family name
!!v           stnid           station ID of observation
!!v           varno           BUFR code (if <=0, to search through all codes to obtain first
!!v                           between 10000 and 16000)
!!v           nlev            number of levels in the observation
!!v           ndim            number of dimensions for the retrieved data in
!!v                           each report (e.g. ndim=1 for std, ndim=2 for
!!v                           averagine kernels) 
!!v           bkstp_opt       bkstp number of requested block if BURP file type (optional)
!!v           block_opt       block type of requested block if BURP file type (optional)
!!v                           Valid values are 'DATA', 'INFO', '3-D', and 'MRQR', indicated
!!v                           by the two rightmost bits of bknat.
!!v           match_nlev_opt  determines if the report matching criteria includes checking
!!v                           if the report number of levels is the same as the input
!!v                           argument nlev (optional)
!!v           codtyp_opt      optional CODTYP list for search (optional)
!!
!! Output:
!!v           obsdata       struct_oss_obsdata object
!!v
!--------------------------------------------------------------------------
  function obsf_obsSub_read(obsfam,stnid,varno,nlev,ndim,bkstp_opt,block_opt,match_nlev_opt, &
                              codtyp_opt) result(obsdata)

    implicit none

    character(len=*), intent(in)  :: obsfam
    character(len=9), intent(in)  :: stnid
    integer, intent(in)           :: varno,nlev,ndim
    integer, intent(in), optional :: bkstp_opt,codtyp_opt(:)
    logical, intent(in), optional :: match_nlev_opt
    character(len=4), intent(in), optional :: block_opt
    type(struct_oss_obsdata) :: obsdata

    character(len=maxLengthFilename) :: filename
    logical :: found
    character(len=10) :: obsFileType

    filename = obsf_get_filename(obsfam,found)

    if (found) then
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


!--------------------------------------------------------------------------
!! *Purpose*: Add or modify data in observational files from data stored
!!            in a struct_oss_obsdata object.
!!
!! @author M. Sitwell, ARQI/AQRD, Sept 2016
!!
!! Input:
!!v           obsdata       Input struct_oss_obsdata object for varno.
!!v           obsfam        observation family name
!!v           varno         BUFR descriptors. Number of elements must be 
!!v                         max(1,obsdata%dim2)
!!v           bkstp_opt     bkstp number of requested block if BURP file type (optional)
!!v           block_opt     block type of requested block if BURP file type (optional)
!!v                         Valid values are 'DATA', 'INFO', '3-D', and 'MRQR', indicated
!!v                         by the two rightmost bits of bknat.
!!v           multi_opt     Indicates if intended report are for 'UNI' or 'MULTI' level data (optional)
!!v
!!
!! Output: 
!!v           nrep_modified Number of modified reports
!!
!--------------------------------------------------------------------------
  function obsf_obsSub_update(obsdata,obsfam,varno,bkstp_opt,block_opt,multi_opt) &
                           result(nrep_modified)

    implicit none

    type(struct_oss_obsdata), intent(inout) :: obsdata
    character(len=*), intent(in) :: obsfam
    integer, intent(in) :: varno(:)
    integer, intent(in), optional :: bkstp_opt
    character(len=4), intent(in), optional :: block_opt
    character(len=*), intent(in), optional :: multi_opt
    integer :: nrep_modified

    integer :: ierr,nrep_modified_global
    character(len=maxLengthFilename) :: filename
    logical :: found
    character(len=10) :: obsFileType

    filename = obsf_get_filename(obsfam,found)

    if (found) then
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

end module obsFiles_mod
