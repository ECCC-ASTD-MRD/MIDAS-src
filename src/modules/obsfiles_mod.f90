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
!! MODULE obsFiles_mod (prefix="obsf")
!!
!! *Purpose*: High-level module to handle reading/writing of observations that
!!            can be stored in one of several different formats. Currently, the
!!            only supported formats are:
!!
!!            1. BURP
!!            2. CMA (binary format of obsSpaceData contents)
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
  use cmaFiles_mod
  use bufr_mod

  implicit none
  save
  private

  ! public variables
  public :: obsf_nfiles, obsf_cfilnam

  ! public procedures
  public :: obsf_setup, obsf_filesSplit, obsf_getFileType, obsf_fileTypeIsBurp
  public :: obsf_readFiles, obsf_writeFiles

  character(len=10) :: obsFileType
  logical           :: obsFilesSplit
  logical           :: initialized = .false.

  integer, parameter :: jpfiles=100
  integer, parameter :: maxLengthFilename=1060
  integer :: obsf_nfiles
  character(len=maxLengthFilename) :: obsf_cfilnam(jpfiles)
  character(len=2)   :: obsf_cfamtyp(jpfiles)

  character(len=48)  :: obsFileMode

contains

  subroutine obsf_setup(dateStamp_out,obsFileMode_in)
    implicit none

    ! arguments
    integer :: dateStamp_out
    character(len=*) :: obsFileMode_in

    obsFileMode = trim(obsFileMode_in)

    !
    ! Initialize file names and the file type
    !
    call obsf_setupFileNames()
    call obsf_determineFileType(obsFileType)

    !
    ! Determine if obsFiles are split
    !
    if ( obsFileType == 'BURP' ) then
      obsFilesSplit = .true.
    else if ( obsFileType == 'CMA' ) then
      obsFilesSplit = .false.
    else if ( obsFileType == 'SQLITE' ) then
      call utl_abort('obsf_setup: SQLITE observation file type not yet implemented')
    else
      call utl_abort('obsf_setup: invalid observation file type: ' // trim(obsFileType))
    end if

    !
    ! Do some setup of observation files
    !
    if ( obsFileType == 'BURP' ) then
      call brpf_getDateStamp( dateStamp_out, obsf_cfilnam(1) )
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


  subroutine obsf_getFileType(obsFileType_out)
    implicit none

    ! arguments
    character(len=*) :: obsFileType_out

    if ( .not.initialized ) call utl_abort('obsf_getFileType: obsFiles_mod not initialized!')

    obsFileType_out = obsFileType
        
  end subroutine obsf_getFileType


  function obsf_fileTypeIsBurp() result(fileTypeIsBurp)
    implicit none
 
    ! arguments
    logical :: fileTypeIsBurp

    if ( .not.initialized ) call utl_abort('obsf_readFiles: obsFiles_mod not initialized!')

    fileTypeIsBurp = ( trim(obsFileType) == 'BURP' )

  end function obsf_fileTypeIsBurp


  subroutine obsf_readFiles(obsSpaceData)
    implicit none

    ! arguments
    type(struct_obs) :: obsSpaceData

    ! locals
    integer :: fileIndex

    if ( .not.initialized ) call utl_abort('obsf_readFiles: obsFiles_mod not initialized!')

    if ( obsFileType == 'BURP' ) then

      do fileIndex = 1, obsf_nfiles
        call brpf_readFile(obsSpaceData, obsf_cfilnam(fileIndex), obsf_cfamtyp(fileIndex), fileIndex)
      end do

    else if ( obsFileType == 'SQLITE' ) then

      do fileIndex = 1, obsf_nfiles
        !call sql_readFile(obsSpaceData, obsf_cfilnam(fileIndex), obsf_cfamtyp(fileIndex), fileIndex)
      end do

    else if ( obsFileType == 'CMA' ) then

      ! read same global CMA file on all mpi tasks
      call cma_readFiles(obsSpaceData)

    end if

  end subroutine obsf_readFiles


  subroutine obsf_writeFiles(obsSpaceData,HXensT_mpiglobal_opt,asciDumpObs_opt)
    implicit none

    ! arguments
    type(struct_obs)           :: obsSpaceData
    real(8), pointer, optional :: HXensT_mpiglobal_opt(:,:)
    logical, optional          :: asciDumpObs_opt

    ! locals
    integer :: fileIndex

    if ( .not.initialized ) call utl_abort('obsf_writeFiles: obsFiles_mod not initialized!')

    if ( obsFileType == 'BURP' ) then

      if (trim(obsFileMode) == 'analysis') call vint3dfd(obs_oma,obsSpaceData)
      call vint3dfd(obs_omp,obsSpaceData)
      if (trim(obsFileMode) == 'analysis' .or. trim(obsFileMode) == 'FSO') call setassflg(obsSpaceData)
      call flaguvtofd_obsdat(obsSpaceData)

      do fileIndex = 1, obsf_nfiles
        call brpf_updateFile(obsSpaceData, obsf_cfilnam(fileIndex), obsf_cfamtyp(fileIndex), fileIndex)
      end do

      if ( present(HXensT_mpiglobal_opt) .and. mpi_myid == 0 ) then
        call obsf_writeHX(obsSpaceData, HXensT_mpiglobal_opt)
      end if

    else if ( obsFileType == 'SQLITE' ) then

      if (trim(obsFileMode) == 'analysis') call vint3dfd(obs_oma,obsSpaceData)
      call vint3dfd(obs_omp,obsSpaceData)
      if (trim(obsFileMode) == 'analysis' .or. trim(obsFileMode) == 'FSO') call setassflg(obsSpaceData)
      call flaguvtofd_obsdat(obsSpaceData)

      do fileIndex = 1, obsf_nfiles
        !call sql_updateFile(obsSpaceData, obsf_cfilnam(fileIndex), obsf_cfamtyp(fileIndex), fileIndex)
      end do

      if ( present(HXensT_mpiglobal_opt) .and. mpi_myid == 0 ) then
        call obsf_writeHX(obsSpaceData, HXensT_mpiglobal_opt)
      end if

    else if ( obsFileType == 'CMA' ) then

      ! only 1 mpi task should do the writing
      call cma_writeFiles(obsSpaceData,HXensT_mpiglobal_opt)

    end if

    if ( present(asciDumpObs_opt) ) then
      if ( asciDumpObs_opt ) call obsf_writeAsciDump(obsSpaceData)
    end if

  end subroutine obsf_writeFiles


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
      write(cmyidy,'(I4.4)') (mpi_npey - mpi_myidy)
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
    !
    ! obsf_setupFileNames - initialze obs file names
    !
    ! locals
    character(len=20)   :: clvalu(jpfiles)
    character(len=2)    :: cfami(jpfiles)
    character(len=4)    :: cmyidx, cmyidy
    character(len=9)    :: cmyid
    character(len=256)  :: obsDirectory
    character(len=maxLengthFilename) :: fileName   !! the length should be more than len(obsDirectory)+1+len(clvalu)+1+len(cmyid)
    character(len=256)               :: fileNamefull
    logical :: isExist_L 
    integer :: fileIndex

    write(cmyidy,'(I4.4)') (mpi_npey - mpi_myidy)
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
    ! file name for CMA format used by EnKF
    clvalu(71) = 'cmaheader'

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
    ! dummy family type for CMA, since it contains all families
    cfami(71)  = 'XX'

    obsDirectory = 'obs'

    obsf_nfiles = 0
    obsf_cfilnam(1) = 'DUMMY_FILE_NAME'
    do fileIndex = 1, jpfiles 
      if(clvalu(fileIndex) == '') exit
      fileName = trim(obsDirectory) // '/' // trim(clvalu(fileIndex)) // '_' // trim(cmyid)
      fileNameFull = ram_fullWorkingPath(fileName,noAbort_opt=.true.)

      inquire(file=trim(fileNameFull),exist=isExist_L)
      if (.not. isExist_L ) then
        fileName=trim(obsDirectory)//'/'//trim(clvalu(fileIndex))
        fileNameFull = ram_fullWorkingPath(fileName, noAbort_opt=.true.)
        inquire(file=trim(fileNameFull), exist=isExist_L)
      end if

      if ( isExist_L ) then
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


  subroutine obsf_determineFileType(obsFileType_out)
    implicit none

    ! arguements
    character(len=10) :: obsFileType_out

    ! locals
    integer :: ierr, procID, unitFile, all_nfiles(0:(mpi_nprocs-1))
    integer :: fnom, fclos
    logical :: fileExists
    character(len=20) :: fileStart

    call rpn_comm_allgather(obsf_nfiles, 1, 'MPI_INTEGER', &
                            all_nfiles,  1, 'MPI_INTEGER', &
                            'GRID',ierr)
    fileExists = .false.
    PROCID_LOOP: do procID = 0, (mpi_nprocs-1)
      if ( all_nfiles(procID) > 0 ) then
        fileExists = .true.
        exit PROCID_LOOP
      end if
    end do PROCID_LOOP

    if ( .not.fileExists ) then
      call utl_abort('obsf_determineFileType: No observation files found')
    end if

    write(*,*) 'obsf_setupFileNames: read obs file that exists on mpi task id: ', procID
    if ( mpi_myid == procID ) then
      if ( index(obsf_cfilnam(1), 'cma' ) > 0 ) then
        obsFileType = 'CMA'
      else
        ! check first characters of file
        unitFile = 0
        ierr = fnom(unitFile, obsf_cfilnam(1), 'FTN+SEQ+FMT+R/O', 0)
        read(unitFile,'(A)') fileStart
        ierr = fclos(unitFile)
        !write(*,*) 'obsf_determineFileType: fileStart = |||', trim(fileStart), '|||'
        if ( index( fileStart, 'SQLite format 3' ) > 0 ) then
          obsFileType = 'SQLITE'
        else if ( index( fileStart, 'XDF0BRP0' ) > 0 ) then
          obsFiletype = 'BURP'
        else
          call utl_abort('obsf_determineFileType: unknown obs file type')
        end if
      end if
    end if

    call rpn_comm_bcastc(obsFileType , 10, 'MPI_CHARACTER', procID, 'GRID', ierr)
    write(*,*) 'obsf_setupFileNames: obsFileType = ', obsFileType

  end subroutine obsf_determineFileType


  SUBROUTINE FLAGUVTOFD_OBSDAT(obsSpaceData)
!
!**s/r FLAGUVTOFD_OBSDAT  - Update WIND DIRECTION AND SPEED FLAGS
!
!
!Author  : P. Koclas *CMC/CMDA  April 2013
!
!
!Arguments
!
      IMPLICIT NONE
!
      type(struct_obs) :: obsSpaceData
      INTEGER :: IUU,IVV,IFF,IDD
      INTEGER :: FLAGU,FLAGV,NEWFLAG
      INTEGER :: INDEX_HEADER,ISTART,IEND,jwintyp
      INTEGER :: INDEX_BODY,INDEX_BODY2
      REAL*8  :: ZLEVU
      LOGICAL ::  LLOK
      CHARACTER*9 :: STID
!-----------------------------------------------------------------------
!
      WIND_TYPE: do jwintyp=1,2

         if (jwintyp == 1) then
            IUU=BUFR_NEUU
            IVV=BUFR_NEVV
            IDD=BUFR_NEDD
            IFF=BUFR_NEFF
         else
            IUU=BUFR_NEUS
            IVV=BUFR_NEVS
            IDD=BUFR_NEDS
            IFF=BUFR_NEFS
         end if
!
!
!
         BODY: DO INDEX_BODY=1,obs_numBody(obsSpaceData)

            LLOK= ( obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY) == IUU)

             FLAGU=-1
    !----------------
            IF ( LLOK ) THEN
    !----------------
               INDEX_HEADER = obs_bodyElem_i(obsSpaceData,OBS_HIND,INDEX_BODY)
               ISTART       = obs_headElem_i(obsSpaceData,OBS_RLN,INDEX_HEADER)
               IEND=obs_headElem_i(obsSpaceData,OBS_NLV,INDEX_HEADER) +ISTART-1
       STID=obs_elem_c(obsSpaceData,'STID',INDEX_HEADER)


               ZLEVU = obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY)
!
!****************************************************************************
!  GET FLAG OF U COMPONENT
!***********************************************************************
!
       FLAGU=obs_bodyElem_i(obsSpaceData,OBS_FLG,INDEX_BODY)

               BODY_2: DO INDEX_BODY2=ISTART,IEND
                  IF ( ( obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY2) == IVV) &
                 .AND. ( obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY2) == ZLEVU) ) THEN
!
!****************************************************************************
!  GET FLAG OF V COMPONENT
!***********************************************************************
!
                     FLAGV= obs_bodyElem_i(obsSpaceData,OBS_FLG,INDEX_BODY2)
                     NEWFLAG =IOR(FLAGU,FLAGV)
!   
                  END IF
               END DO BODY_2
!
!***********************************************************************
!                UPDATE FLAGS OF DIRECTION AN SPEED
!***********************************************************************
!
               BODY_2_2: DO INDEX_BODY2=ISTART,IEND
       !===============================================
                  IF ((obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY2) == IDD) &
                 .AND. obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY2) == ZLEVU ) THEN

                     NEWFLAG =IOR(FLAGU,FLAGV)
                     call obs_bodySet_i(obsSpaceData, OBS_FLG, INDEX_BODY2, NEWFLAG) 

                  END IF
	       !===============================================

       !===============================================
                  IF ((obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY2) == IFF) &
                 .AND. obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY2) == ZLEVU ) THEN

                     NEWFLAG =IOR(FLAGU,FLAGV)
                     call obs_bodySet_i(obsSpaceData,OBS_FLG,INDEX_BODY2, NEWFLAG)
                  END IF
	       !===============================================
               END DO BODY_2_2

	    !----------------
            END IF
	    !----------------

         END DO BODY

      END DO WIND_TYPE

  END SUBROUTINE FLAGUVTOFD_OBSDAT


  SUBROUTINE VINT3DFD(elem_i,obsSpaceData)
      !
      ! s/r VINT3DFD  - Computation of DIRECTION AND SPEED RESIDUALS
      !
      ! Author  : P. Koclas *CMC/AES  September 1999
      ! Revision:
      !     1.0  P. Koclas CMC :  September 2000
      !                 -remove quality control flag and (ff dd) component initializtions
      !          JM Belanger CMDA/SMC  Jan 2001
      !                   . 32 bits conversion
      !
      !     Purpose:  -Compute direction and speed residuals from u and
      !                v residuals.
      !
      implicit none

      type(struct_obs) :: obsSpaceData
      integer, intent(in) :: elem_i
      INTEGER IUU,IVV,IFF,IDD
      INTEGER INDEX_HEADER,ISTART,IEND,jwintyp
      INTEGER INDEX_BODY,INDEX_BODY2
      REAL*8 ZLEVU
      REAL*8 MODUL,ANG,UU,VV
      LOGICAL LLOK

      WIND_TYPE: do jwintyp=1,2

         if (jwintyp == 1) then
            IUU=BUFR_NEUU
            IVV=BUFR_NEVV
            IDD=BUFR_NEDD
            IFF=BUFR_NEFF
         else
            IUU=BUFR_NEUS
            IVV=BUFR_NEVS
            IDD=BUFR_NEDS
            IFF=BUFR_NEFS
         end if

         ! Process all data within the domain of the model

         BODY: DO INDEX_BODY=1,obs_numBody(obsSpaceData)
            LLOK= (obs_bodyElem_i(obsSpaceData,OBS_ASS,INDEX_BODY) == 1)  &
            .AND. (obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY) == IUU)
            IF ( LLOK ) THEN
               INDEX_HEADER = obs_bodyElem_i(obsSpaceData,OBS_HIND,INDEX_BODY)
               ISTART=obs_headElem_i(obsSpaceData,OBS_RLN,INDEX_HEADER)
               IEND=obs_headElem_i(obsSpaceData,OBS_NLV,INDEX_HEADER) +ISTART-1
               ZLEVU = obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY)
               UU=-obs_bodyElem_r(obsSpaceData,elem_i,INDEX_BODY) +  &
                   obs_bodyElem_r(obsSpaceData,OBS_VAR,INDEX_BODY)
               BODY_2: DO INDEX_BODY2=ISTART,IEND
                  IF ((obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY2) == IVV)  &
                 .AND.(obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY2) == ZLEVU)) THEN
                   VV=-obs_bodyElem_r(obsSpaceData,elem_i,INDEX_BODY2) +  &
                       obs_bodyElem_r(obsSpaceData,OBS_VAR,INDEX_BODY2)

                     ! 1-calculate angle

                     MODUL=SQRT((UU**2)+(VV**2))
                     IF (MODUL == 0.) THEN
                        ANG=0.0D0
                     ELSE
                        ANG=ATAN2(VV,UU)
                        ANG= (270.0D0 - ANG  * MPC_DEGREES_PER_RADIAN_R8 )

                        ! 2-Change to meteorological definition of wind direction.

                        IF (ANG > 360.0D0) ANG=ANG-360.0D0
                        IF (ANG <= 0.0D0)   ANG=ANG+360.0D0
                     END IF
   
                  END IF
               END DO BODY_2

               ! insert resduals into obsSpaceData

               BODY_2_2: DO INDEX_BODY2=ISTART,IEND
                  IF ((obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY2) == IDD)  &
                 .AND. obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY2) == ZLEVU ) THEN

                     call obs_bodySet_r(obsSpaceData, elem_i, INDEX_BODY2,    &
                          obs_bodyElem_r(obsSpaceData,OBS_VAR,INDEX_BODY2) - ANG )

                     IF ( obs_bodyElem_r(obsSpaceData,elem_i,INDEX_BODY2) >  180.0d0)  &
                        call obs_bodySet_r(obsSpaceData, elem_i, INDEX_BODY2,   &
                                       obs_bodyElem_r(obsSpaceData,elem_i,INDEX_BODY2)-360.0d0)
                     IF ( obs_bodyElem_r(obsSpaceData,elem_i,INDEX_BODY2) <= -180.0d0)  &
                        call obs_bodySet_r(obsSpaceData, elem_i, INDEX_BODY2,  &
                                       obs_bodyElem_r(obsSpaceData,elem_i,INDEX_BODY2)+360.0d0)

                      call obs_bodySet_r(obsSpaceData, elem_i, INDEX_BODY2, -1.0d0*  &
                                         obs_bodyElem_r(obsSpaceData,elem_i,INDEX_BODY2))

                      call obs_bodySet_r(obsSpaceData,OBS_OER,INDEX_BODY2,1.0d0)
                      call obs_bodySet_i(obsSpaceData,OBS_ASS,INDEX_BODY2, 1)
                      call obs_bodySet_i(obsSpaceData,OBS_FLG,INDEX_BODY2, 0)
                  END IF
                  IF ((obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY2) == IFF)  &
                 .AND. obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY2) == ZLEVU ) THEN
                     call obs_bodySet_r(obsSpaceData,elem_i, INDEX_BODY2,   &
                          obs_bodyElem_r(obsSpaceData,OBS_VAR,INDEX_BODY2) - MODUL)
                     call obs_bodySet_r(obsSpaceData,OBS_OER,INDEX_BODY2,1.0d0)
                     call obs_bodySet_i(obsSpaceData,OBS_ASS,INDEX_BODY2, 1)
                     call obs_bodySet_i(obsSpaceData,OBS_FLG,INDEX_BODY2, 0)
                  END IF
               END DO BODY_2_2
            END IF

         END DO BODY

      END DO WIND_TYPE

  END SUBROUTINE VINT3DFD

  subroutine setassflg(obsSpaceData)
    ! Purpose:  Set banco quality control bit #12 for all data assimilated
    !           by current analysis.
    implicit none

    type(struct_obs) :: obsSpaceData
    integer :: index_body

    ! Process all data
    do index_body=1,obs_numBody(obsSpaceData)
      if (obs_bodyElem_i(obsSpaceData,OBS_ASS,index_body) == 1)  then
        call obs_bodySet_i(obsSpaceData,OBS_FLG,index_body,ibset( obs_bodyElem_i(obsSpaceData,OBS_FLG,index_body), 12 ))
      end if
    end do

  end subroutine setassflg

end module obsFiles_mod
