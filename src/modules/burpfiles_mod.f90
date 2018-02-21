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
!! MODULE burpFiles (prefix="burp")
!!
!! *Purpose*: To store the filenames of the burp observation files and call
!!            subroutines in readBurp to read and update burp files.
!!
!--------------------------------------------------------------------------
module burpFiles_mod
  use codePrecision_mod
  use MathPhysConstants_mod
  use EarthConstants_mod
  use utilities_mod
  use mpivar_mod
  use obsSpaceData_mod
  use burpread_mod
  use ramDisk_mod
  use bufr_mod
  use utilities_mod
  use obsSubSpaceData_mod
  use burp_module
  use mpi_mod, only: mpi_myid
  
  implicit none
  save
  private

  ! public variables
  public :: burp_nfiles, burp_cfilnam

  ! public procedures
  public :: burp_setupfiles, burp_readfiles, burp_updatefiles
  public :: burp_chem_read_all, burp_chem_update_all, burp_split, burp_splitMode

  integer, parameter :: jpfiles=64
  integer :: burp_nfiles
  integer, parameter :: jmaxburpfilename=1060

  character(len=jmaxburpfilename) :: burp_cfilnam(jpfiles)
  character(len=2)   :: burp_cfamtyp(jpfiles)
  character(len=256) :: burp_split_mode
  character(len=48)  :: burpFileMode

  logical :: burp_split_L

contains

  SUBROUTINE burp_setupfiles(datestamp, burpFileMode_in)
    implicit none
    !s/r  burp_setupfiles -INITIALZE BURP FILE NAMES and return datestamp
    !
    !Revisions:
    !          Ping Du (Fall 2014)
    !          -- Added CVALUE for brpch and CFAMI for CH.
    !          Yves Rochon (Jan 2015)
    !          -- Removed CVALUEs for brpo3 and brpoz
    !          -- Removed CFAMIs for OZ
    !          Mike Sitwell (Jan 2015)
    !          -- Initialized CVALUEs and CFAMI to ''. This prevents 
    !             the possibility of related issues with statements such as
    !                  IF (CVALUE(JJ) == '') EXIT
    !          Y. Rochon ARQI/AQRD, July 2016
    !          -- Allows nearest reference time to differ from synoptic hour when the
    !             resume record is present.
    !             The resume record first encountered must contain the desired reference time.
    !             

    character(len=*), intent(in) :: burpFileMode_in

    INTEGER IER,INBLKS,nulburp,JJ
    INTEGER FNOM,FCLOS,NUMBLKS
    CHARACTER(len=20) :: CLVALU(JPFILES)
    CHARACTER(len=2) :: CFAMI(JPFILES)
    CHARACTER(len=4) :: cmyidx, cmyidy
    CHARACTER(len=9) :: cmyid
    CHARACTER(len=256) :: burp_directory
    CHARACTER(len=jmaxburpfilename) :: burpin   !! the length should be more than len(burp_directory)+1+len(clvalu)+1+len(cmyid)
    character(len=256)              :: burpinfull
    LOGICAL   isExist_L 

    INTEGER KTIME,KDATE,KDATE_RECV,KTIME_RECV
    INTEGER IHANDL,ILONG,DATESTAMP
    INTEGER ITIME,IFLGS,IDBURP,ILAT,ILON,IDX,IDY
    INTEGER IALT,IDELAY,IDATE,IRS,IRUNN,INBLK,ISUP,IXAUX
    INTEGER INSUP,INXAUX

    INTEGER, ALLOCATABLE :: IBUF(:)
    INTEGER INRECS
    INTEGER MRFCLS,MRFOPN,MRFOPC,MRBHDR,MRFLOC,MRFGET
    INTEGER MRFMXL, status, length_burp_split, length_burpDatestamp

    INTEGER ISTAMPOBS,INEWHH,NEWDATE,nresume
    REAL*8 DELHH
    INTEGER       IVALS
    CHARACTER*9   CLSTNID

    character(len=12) :: obs_yyyymmddhhmm_char
    integer :: length_obs_yyyymmddhhmm, obsCentralDate_yyyyddmm, obsCentralDate_hhmm

    logical :: datestampFromBurpFile

    EXTERNAL FCLOS,FNOM,MRFCLS,MRFOPN,MRFOPC,MRBHDR,MRFLOC,MRFGET,MRFMXL,NUMBLKS

    !
    !- Setup the mode
    !
    burpFileMode = burpFileMode_in

    !
    !- Determine if the observation files are already split by subdomain
    !
    status = 0
    call get_environment_variable('OAVAR_BURP_SPLIT',burp_split_mode,length_burp_split,status,.true.)

    if ( status == 1 ) then
      write(*,*) 'burp_setupfiles: The environment variable OAVAR_BURP_SPLIT has not been detected!'
      write(*,*) '                 The observation files are NOT split'
      burp_split_L = .false.  
      burp_split_mode = 'DEFAULT'
      ! At this point the code only supports split files
      write(*,*) 'burp_setupfiles: Expecting split burp files!'
      call utl_abort('burp_setupfiles')
    else
      write(*,*)
      write(*,*) 'burp_setupfiles: The environment variable OAVAR_BURP_SPLIT has correctly been detected'
      select case ( trim(burp_split_mode) )
      case ('yes')
        write(*,*) 'burp_setupfiles: The observation files ARE split. Assuming ROUNDROBIN strategy'
        burp_split_L = .true.
        burp_split_mode = 'ROUNDROBIN'
      case ('no')
        write(*,*) 'burp_setupfiles: The observation files are NOT split'
        burp_split_L = .false.
        burp_split_mode = 'DEFAULT'
        ! At this point the code only supports split files
        write(*,*) 'burp_setupfiles: Expecting split burp files!'
        call utl_abort('burp_setupfiles')
      case ('ROUNDROBIN')
        write(*,*) 'burp_setupfiles: The observation files ARE split using the ROUNDROBIN strategy'
        burp_split_L = .true.
      case ('LATLONTILES')
        write(*,*) 'burp_setupfiles: The observation files ARE split by LAT-LON tiles'
        write(*,*) 'burp_setupfiles: LAT-LON TILES NO LONGER PERMITTED, USE ROUND ROBIN!'
        call utl_abort('burp_setupfiles')
      case default
        write(*,*) 'burp_setupfiles: Unknown burp_split_mode ', trim(burp_split_mode)
        call utl_abort('burp_setupfiles')
      end select
    end if

    !
    !- Determine datestamp
    !
    status = 0
    call get_environment_variable('OAVAR_OBS_CENTRAL_DATE',obs_yyyymmddhhmm_char,length_obs_yyyymmddhhmm,status,.true.)

    if ( status == 1 ) then
      write(*,*) 'burp_setupfiles: The environment variable OAVAR_OBS_CENTRAL_DATESTAMP has not been detected!'
      write(*,*) '                 The datestamp will be obtained from the burp files'
      datestampFromBurpFile = .true.
    else
      write(*,*)
      write(*,*) 'burp_setupfiles: The environment variable OAVAR_OBS_CENTRAL_DATE has correctly been detected'
      write(*,*) '                 The date found is: ', trim(obs_yyyymmddhhmm_char)
      datestampFromBurpFile = .false.
      read(obs_yyyymmddhhmm_char(1: 8),'(i8)') obsCentralDate_yyyyddmm
      read(obs_yyyymmddhhmm_char(9:12),'(i4)') obsCentralDate_hhmm
    end if

    !
    !- Process the input burp files
    !
    write(cmyidy,'(I4.4)') (mpi_npey-mpi_myidy)
    write(cmyidx,'(I4.4)') (mpi_myidx+1)
    cmyid  = trim(cmyidx)//'_'//trim(cmyidy)

    CLVALU(:)=''
    CLVALU( 1) = 'brpuan'  
    CLVALU( 2) = 'brpuas'  
    CLVALU( 3) = 'brpai'  
    CLVALU( 4) = 'brpain'  
    CLVALU( 5) = 'brpais'  
    CLVALU( 6) = 'brpaie'  
    CLVALU( 7) = 'brpaiw'  
    CLVALU( 8) = 'brpsfc'  
    CLVALU( 9) = 'brpsf'  
    CLVALU(10) = 'brptov'  
    CLVALU(11) = 'brpssmis'
    CLVALU(12) = 'brpairs'
    CLVALU(13) = 'brpto_amsua'
    CLVALU(14) = 'brpto_amsub'
    CLVALU(15) = 'brpcsr'
    CLVALU(16) = 'brpiasi'
    CLVALU(17) = 'brpatms'
    CLVALU(18) = 'brpcris'
    CLVALU(19) = 'brpsw'  
    CLVALU(20) = 'brpswgoes9'  
    CLVALU(21) = 'brpswgoese'  
    CLVALU(22) = 'brpswgoesw'  
    CLVALU(23) = 'brpswmodis'  
    CLVALU(24) = 'brpswmtsate'  
    CLVALU(25) = 'brpswmtsatw'  
    CLVALU(26) = 'brpgo'  
    CLVALU(27) = 'brpsc'  
    CLVALU(28) = 'brppr'  
    CLVALU(29) = 'brpro'  
    CLVALU(30) = 'brphum'  
    CLVALU(31) = 'brpsat'  
    CLVALU(32) = 'brpssm'  
    CLVALU(33) = 'brpgp'  
    CLVALU(34) = 'brpch' 
    CLVALU(35) = 'brpua'  


    CFAMI(:)   = ''
    CFAMI( 1)  = 'UA' 
    CFAMI( 2)  = 'UA' 
    CFAMI( 3)  = 'AI' 
    CFAMI( 4)  = 'AI' 
    CFAMI( 5)  = 'AI'
    CFAMI( 6)  = 'AI'
    CFAMI( 7)  = 'AI'
    CFAMI( 8)  = 'SF' 
    CFAMI( 9)  = 'SF' 
    CFAMI(10)  = 'TO' 
    CFAMI(11)  = 'TO' 
    CFAMI(12)  = 'TO' 
    CFAMI(13)  = 'TO' 
    CFAMI(14)  = 'TO' 
    CFAMI(15)  = 'TO' 
    CFAMI(16)  = 'TO' 
    CFAMI(17)  = 'TO' 
    CFAMI(18)  = 'TO' 
    CFAMI(19)  = 'SW' 
    CFAMI(20)  = 'SW' 
    CFAMI(21)  = 'SW' 
    CFAMI(22)  = 'SW' 
    CFAMI(23)  = 'SW' 
    CFAMI(24)  = 'SW' 
    CFAMI(25)  = 'SW' 
    CFAMI(26)  = 'GO' 
    CFAMI(27)  = 'SC' 
    CFAMI(28)  = 'PR' 
    CFAMI(29)  = 'RO' 
    CFAMI(30)  = 'HU' 
    CFAMI(31)  = 'ST' 
    CFAMI(32)  = 'MI' 
    CFAMI(33)  = 'GP' 
    CFAMI(34)  = 'CH' 
    CFAMI(35)  = 'UA' 

    IER =MRFOPC('MSGLVL','FATAL')

    burp_directory = 'obs'

    IVALS=8
    KDATE=-9999
    KTIME=-9999
    nresume=0
    burp_nfiles=0
    DO JJ=1,JPFILES 
      IF(CLVALU(JJ) == '') EXIT
      nulburp=0
      burpin=trim(burp_directory)//'/'//trim(CLVALU(JJ))//'_'//trim(cmyid)
      burpinFull = ram_fullWorkingPath(burpin,noAbort_opt=.true.)

      INQUIRE(FILE=trim(burpinFull),EXIST=isExist_L)
      IF (.NOT. isExist_L )THEN
        burpin=trim(burp_directory)//'/'//trim(CLVALU(JJ))
        burpinFull = ram_fullWorkingPath(burpin,noAbort_opt=.true.)
        INQUIRE(FILE=trim(burpinFull),EXIST=isExist_L)
      END IF
      IF ( isExist_L )THEN
        IER=FNOM(nulburp,burpinFull,'RND+OLD',0)
        WRITE(*,*)' Open File : ',trim(burpinFull)
        IF ( IER == 0 ) THEN
          INBLKS= -1
          INBLKS=NUMBLKS(nulburp)
          IF ( INBLKS > 0 ) THEN
            INRECS=MRFOPN(NULBURP,'READ')
            ILONG =MRFMXL(NULBURP)
            ALLOCATE(IBUF(ILONG + 20))
            IBUF(1)=ILONG + 20
            IHANDL  =MRFLOC(NULBURP,0,'>>*******',-1,-1,-1,-1,-1,-1,0)
            IF ( IHANDL < 0 ) THEN
              IHANDL=MRFLOC(NULBURP,0,'*********',-1,-1,-1,-1,-1,-1,0)
            ELSE
              nresume=nresume+1
            END IF
            IF ( IHANDL < 0 ) THEN
              WRITE(*,*) 'AUCUN ENREGISTREMENT VALIDE DANS LE FICHIER BURP'
            ELSE
              burp_nfiles=burp_nfiles + 1
              burp_cfilnam(burp_nfiles)=burpinFull
              burp_cfamtyp(burp_nfiles)=CFAMI(JJ)
              if ((kdate < 0.and.ktime < 0).or.nresume == 1) then 
                INSUP=0
                INXAUX=0
                IER=MRFGET(IHANDL,IBUF)
                IER=MRBHDR(IBUF,ITIME,IFLGS,CLSTNID,IDBURP,ILAT,   &
                     ILON,IDX,IDY, IALT,IDELAY,IDATE,IRS,IRUNN,INBLK, &
                     ISUP,INSUP,IXAUX,INXAUX)
                KTIME=ITIME
                KDATE=IDATE
                if (nresume == 1) nresume=2
              end if
            END IF
            DEALLOCATE(IBUF)
            IER=MRFCLS(NULBURP)
          END IF
        END IF
        IER= FCLOS(nulburp)
      END IF
    END DO

    WRITE(*,*) ' '
    WRITE(*,*)' NUMBER OF BURP FILES IS :',burp_nfiles
    WRITE(*,*)'TYPE  NAME '
    WRITE(*,*)'----  ---- '
    DO JJ=1,burp_nfiles
      WRITE(*,'(1X,A2,1X,A128)' ) burp_cfamtyp(JJ),trim(burp_cfilnam(JJ))
    END DO

    !
    !- Set reference datestamp
    !
    if ( datestampFromBurpFile ) then
      ! Make sure all mpi tasks have a valid date (important for split burp files)
      call rpn_comm_allreduce(kdate,kdate_recv,1,"MPI_INTEGER","MPI_MAX","GRID",ier)
      call rpn_comm_allreduce(ktime,ktime_recv,1,"MPI_INTEGER","MPI_MAX","GRID",ier)
      kdate = kdate_recv
      ktime = ktime_recv
       if (nresume >= 1 ) then  
        ier = newdate(datestamp,kdate,ktime*10000,3)
      else
        ! Assumes 6-hour windows with reference times being synoptic times.
        ! Does not require kdate and ktime to be from a resume record.
        ier = newdate(istampobs,kdate,ktime*10000,3)
        delhh = 3.0d0
        call INCDATR (datestamp, istampobs, delhh)
        ier = newdate(datestamp,kdate,inewhh,-3)
        ktime=KTIME/100
        if (ktime >= 21 .or. ktime < 3) then
          ktime = 0
        else if(ktime >= 3 .and. ktime < 9) then
          ktime = 6
        else if(ktime >= 9 .and. ktime < 15) then
          ktime = 12
        else
          ktime = 18
        end if
        ier = newdate(datestamp,kdate,ktime*1000000,3)
        ktime=ktime*100
      end if

    else
      kdate = obsCentralDate_yyyyddmm
      ktime = obsCentralDate_hhmm
      ier = newdate(datestamp,kdate,ktime*10000,3)
    end if

    WRITE(*, *)' BURP FILES VALID DATE (YYYYMMDD) : ',kdate
    WRITE(*, *)' BURP FILES VALID TIME     (HHMM) : ',ktime
    write(*, *)' BURP FILES DATESTAMP             : ',datestamp

  END SUBROUTINE burp_setupfiles

!--------------------------------------------------------------------------
! burp_splitMode
!--------------------------------------------------------------------------
    function burp_SplitMode() result(BurpSplitMode)
      implicit none
      character(len=48) :: BurpSplitMode

      BurpSplitMode = trim(burp_split_mode)
        
    end function burp_SplitMode

!--------------------------------------------------------------------------
! burp_split
!--------------------------------------------------------------------------
    function burp_split() result(BurpSplit)
      implicit none
      logical :: BurpSplit

      BurpSplit = burp_split_L
        
    end function burp_split

    SUBROUTINE burp_readFiles(obsdat)
!
!---------------------------------------------------------------------------------
!      PURPOSE: READ CMC BURP FILES FILL UP OBSSPACE DATA FILE
!
!    ARGUMENTS:
!                        obsdat   - obsdat-file object
!
!       AUTHOR: P. KOCLAS(CMC CMDA)
!
!       Revisions:
!         1 Oct 2013:  S. Macpherson ARMA
!                      Bug fix: Put CALL SET_ERR_GBGPS *after* OBS_OER initialized to 0
!         Dec 2014:    Ping Du, CMDA
!                      - Addition of CH family consideration for scaling
!                        via element BUFR_SCALE_EXPONENT
!
!     NOTE:
!     BURP FILES ARE ASSUMED TO BE PRESENT IN CURRENT WORKING DIRECTORY
!---------------------------------------------------------------------------------
!
      IMPLICIT NONE
      type (struct_obs), intent(inout) :: obsdat

      INTEGER :: IBEG, IEND, NSTN1, NSTN2
      logical :: obs_full,burp_chem

      INTEGER :: J,JO
      REAL(OBS_REAL)  :: MISG

      WRITE(*,*)' '
      WRITE(*,*)'================================================='
      WRITE(*,*)'                burp_readFiles BEGIN             '
      WRITE(*,*)'================================================='
      WRITE(*,*)' '
      MISG = real(MPC_missingValue_R8,OBS_REAL)

      IBEG=obs_numbody(obsdat)
      DO J =1,burp_nfiles

         IBEG=obs_numbody(obsdat) +1
         Nstn1=obs_numheader(obsdat)

         call READBURP(obsdat,burp_cfamtyp(J),burp_cfilnam(J),J)
         Nstn2=obs_numheader(obsdat)
         IEND=obs_numbody(obsdat)

         burp_chem = trim(burp_cfamtyp(J)) == 'CH'

         IF ( trim(burp_cfamtyp(J)) /= 'TO' .and. .not.burp_chem) THEN
            call FDTOUV_OBSDAT(  obsdat,Nstn1+1,Nstn2,MPC_missingValue_R4)
            call ADJUST_HUM_GZ(  obsdat,Nstn1+1,Nstn2)
            call ADJUST_SFVCOORD(obsdat,Nstn1+1,Nstn2)
         END IF
         DO JO=nstn1+1,nstn2
           call obs_headSet_i(obsdat,OBS_OTP,JO,J)
           call obs_setFamily(obsdat,trim(burp_cfamtyp(J)),JO)
            ! For CH family, apply scaling from the element BUFR_SCALE_EXPONENT when present.
            if (burp_chem) call set_scale_chm(obsdat,JO,forward=.true.)
         END DO
         !    initializations
         DO JO=IBEG,IEND
            if ( obs_columnActive_RB(obsdat,OBS_OMA) )  call obs_bodySet_r(obsdat,OBS_OMA ,JO,MISG)
            if ( obs_columnActive_RB(obsdat,OBS_OMP) )  call obs_bodySet_r(obsdat,OBS_OMP ,JO,MISG)
            if ( obs_columnActive_RB(obsdat,OBS_OER) )  call obs_bodySet_r(obsdat,OBS_OER ,JO,MISG)
            if ( obs_columnActive_RB(obsdat,OBS_HPHT) ) call obs_bodySet_r(obsdat,OBS_HPHT,JO,MISG)
            if ( obs_columnActive_RB(obsdat,OBS_WORK) ) call obs_bodySet_r(obsdat,OBS_WORK,JO,MISG)
         END DO

         ! For GP family, initialize OBS_OER to element 15032 (ZTD formal error) 
         ! for all ZTD data (element 15031)
         IF ( trim(burp_cfamtyp(J)) == 'GP') THEN
           print * ,' Initializing OBS_OER for GB-GPS ZTD to formal error (ele 15032)'
           CALL SET_ERR_GBGPS(obsdat,Nstn1+1,Nstn2)
         END IF

      END DO

      WRITE(*,*) '  readburp obs_numheader(obsdat)', obs_numheader(obsdat)
      WRITE(*,*) '  readburp obs_numbody(obsdat)  ', obs_numbody  (obsdat)

      WRITE(*,*)' '
      WRITE(*,*)'================================================='
      WRITE(*,*)'                burp_readFiles END               '
      WRITE(*,*)'================================================='
      WRITE(*,*)' '
      RETURN

end subroutine burp_readFiles


subroutine burp_updateFiles(obsSpaceData)
!
!     PURPOSE: READ OBSDAT AND UPDATE CMC BURP FILES
!
!     ARGUMENTS:
!                   obsSpaceData   - obsdat-file object
!
!       AUTHOR: P. KOCLAS(CMC CMDA)
!
!       Revision:
!                Ping Du, CMDA, Feb-Mar 2015
!                - Added scaling for CH family data
!                  from use of BUFR_SCALE_EXPONENT element when required.
!
!     NOTE:
!     BURP FILES ARE ASSUMED TO BE PRESENT IN CURRENT WORKING DIRECTORY
!
      IMPLICIT NONE
      type (struct_obs), intent(inout) :: obsSpaceData

      INTEGER J
      INTEGER FILENUMB,IBRP1,IER,INRECS,ISTAT,LNMX
      INTEGER  FNOM,FCLOS,MRFCLS,MRFOPN,MRFMXL
      EXTERNAL FNOM,FCLOS,MRFCLS,MRFOPN,MRFMXL
      integer nstn1,nstn2,headerIndex

      call tmg_start(93,'POST_UPDATEBRP')

      if (trim(burpFileMode) == 'analysis') call vint3dfd(obs_oma,obsSpaceData)
      call vint3dfd(obs_omp,obsSpaceData)
      if (trim(burpFileMode) == 'analysis' .or. trim(burpFileMode) == 'FSO') call setassflg(obsSpaceData)
      call flaguvtofd_obsdat(obsSpaceData)
!
!  ------NOTE----------
! currently supported families of data 'UA' 'AI' 'SC' 'SF' 'SW' 'TO' 'CH'
!
!     READ DATA FROM FILES CONTAINED IN ARRAY CLVAL.
!
      WRITE(*,*)' '
      WRITE(*,*)'================================================='
      WRITE(*,*)'                burp_updateFiles BEGIN           '
      WRITE(*,*)'================================================='
      WRITE(*,*)' '

      if ( .not. burp_split_L ) then 
         write(*,*) 'burp_updateFiles: We read/write global observation files'
         CALL obs_expandToMpiGlobal(obsSpaceData)
         IF(mpi_myid /= 0) then
           call tmg_stop(93)
           return
         end if
      else
         ! redistribute obs data to how it was just after reading the files
         call obs_MpiRedistribute(obsSpaceData,OBS_IPF)
      end if

      
      ! CH family: Scaling of the obs related values to be stored in the BURP files

      call obs_set_current_header_list(obsSpaceData,'CH')
      HEADER: do
         headerIndex = obs_getHeaderIndex(obsSpaceData)
         if (headerIndex < 0) exit HEADER
         call set_scale_chm(obsSpaceData,headerIndex,forward=.false.)
      end do HEADER

      do J =1,burp_nfiles
         call update_burp(obsSpaceData,burp_cfamtyp(J),burp_cfilnam(J),J)
      end do

      WRITE(*,*)' '
      WRITE(*,*)'================================================='
      WRITE(*,*)'                burp_updatefiles    END          '
      WRITE(*,*)'================================================='
      WRITE(*,*)' '

      call tmg_stop(93)

END SUBROUTINE burp_updateFiles

  SUBROUTINE  ADJUST_HUM_GZ(obsdat,START,END)
!**s/r ADJUST_HUM_GZ  - Adjust  t-td and GZ in obsdat
!
!
!Author  : P. Koclas *CMC/CMDA  April 2013
!Revision:
!
!*    Purpose:  - Adjust  t-td values to zesmax=30. in obsdat
!                 set Z to GZ                       in obsdat
!
!
!Arguments
! 
!               INPUT:
!                  -OBSDAT    : instance of obsspace_data module object
!                  -START     : FIRST OBERVATION
!                  -END       : LAST  OBERVATION
!
!
      IMPLICIT NONE
      INTEGER  :: START,END

      INTEGER  :: J,JO,RLN,NLV
      INTEGER  :: VARNO
      type (struct_obs), intent(inout):: obsdat

      REAL(OBS_REAL)    :: ZESMAX,GZ,OBSV
      REAL              :: RMIN

      !-----------------------
      ZESMAX=30.0
      !-----------------------
!-----------------------------------------------------------------------
      WRITE(*,*)'   ADJUST_HUM_GZ '
!
!-----------------------------------
!      STN LOOP
!-----------------------------------
!     DO JO=1,obs_numheader(obsdat)
      DO JO=START,END
        RLN=obs_headElem_i(obsdat,OBS_RLN,JO)
        NLV=obs_headElem_i(obsdat,OBS_NLV,JO)
        !=================================
        ! DATA LOOP
        !=================================
        DO J = RLN, NLV + RLN -1

          VARNO=obs_bodyElem_i(obsdat,OBS_VNM,j)
          SELECT CASE(VARNO)
            CASE(BUFR_NEES,BUFR_NESS)
             OBSV=obs_bodyElem_r(obsdat,OBS_VAR,j)
             IF ( OBSV > ZESMAX) THEN
                OBSV=ZESMAX
             END IF
             call obs_bodySet_r(obsdat,OBS_VAR,j, OBSV )
            CASE(BUFR_NEGZ)
             OBSV=obs_bodyElem_r(obsdat,OBS_VAR,j)
             GZ=OBSV*GRAV
             call obs_bodySet_r(obsdat,OBS_VAR,j,GZ )
          END SELECT
!
        END DO
        !=================================
!
      END DO
!-----------------------------------
!
      WRITE(*,*)' DONE   ADJUST_HUM_GZ '
      RETURN
  END SUBROUTINE ADJUST_HUM_GZ


  SUBROUTINE  SET_ERR_GBGPS(obsdat,START,END)
!**s/r SET_ERR_GBGPS  - SET INITIAL ERROR FRO GROUND BASED GPS
!
!
!Author  : P. Koclas *CMC/CMDA  July 2013
!Revision:
!
!*    Purpose:  - PUT 15032 observation element as error of 15031 element  in obsdat
!
!
!Arguments
! 
!               INPUT:
!                  -OBSDAT    : instance of obsspace_data module object
!                  -START     : FIRST OBERVATION
!                  -END       : LAST  OBERVATION
!
      IMPLICIT NONE
      INTEGER  :: START,END

      INTEGER  :: J,JO,RLN,NLV
      INTEGER  :: VARNO
      type (struct_obs), intent(inout):: obsdat

      REAL(OBS_REAL)    :: OBSV

!-----------------------------------------------------------------------
      WRITE(*,*)'   SET_ERR_GBGPS '
!
!-----------------------------------
!      STN LOOP
!-----------------------------------
      DO JO=START,END
        RLN=obs_headElem_i(obsdat,OBS_RLN,JO)
        NLV=obs_headElem_i(obsdat,OBS_NLV,JO)
        !=================================
        ! DATA LOOP
        !=================================
        OBSV=real(MPC_missingValue_R8,OBS_REAL)
        DO J = RLN, NLV + RLN -1

          VARNO=obs_bodyElem_i(obsdat,OBS_VNM,j)
          IF ( VARNO == 15032 ) THEN
             OBSV=obs_bodyElem_r(obsdat,OBS_VAR,j)
             call obs_bodySet_i(obsdat,OBS_VNM,j,999 )
             EXIT
          END IF
!
        END DO
        DO J = RLN, NLV + RLN -1

          VARNO=obs_bodyElem_i(obsdat,OBS_VNM,j)
          IF ( VARNO == 15031 .and. OBSV /= real(MPC_missingValue_R8,OBS_REAL)) THEN
             call obs_bodySet_r(obsdat,OBS_OER,j,OBSV )
             EXIT
          END IF
!
        END DO
        !=================================
!
      END DO
!-----------------------------------
!
      WRITE(*,*)' DONE   SET_ERR_GBGPS '
      RETURN
  END SUBROUTINE SET_ERR_GBGPS

!--------------------------------------------------------------------------
!! *Purpose*:  Apply or unapply scaling to CH observations  by multiplying
!!             (or dividing) with 10^{exponent} where the exponent is from
!!             element BUFR_SCALE_EXPONENT if provided.           
!!
!! @author P. Du *CMC/CMDA Nov. 2014
!!
!! Revisions:
!!v     Y. Rochon, ARQI/AQRD, Feb 2015
!!v       - Generalized except for the assumption that if
!!v         exponents are present, then they must be present
!!v         for all obs and in the same sequential order.
!!v     Y. Rochon, ARQI/AQRD, March/April 2016
!!v       - Added consideration of HBHT (HPHT)
!!v       - Avoid abort when BUFR_SCALE_EXPONENT present but not the mantissa as it could 
!!v         have been filtered in READBURP.
!!v     M. Sitwell, ARQI/AQRD, March 2017
!!v       - Modified to do scaling for a single profile (i.e. one headerIndex)
!!v         instead of for all observations (loop is done outside of set_scale_chm) 
!!
!! Input:
!!v      obsdat         struct_obs instance
!!v      headerIndex    header index in obsdat to apply/unapply scaling to
!!v      forward        applies scaling if .true., unapplies scaling if .false.
!--------------------------------------------------------------------------
  SUBROUTINE  SET_SCALE_CHM(obsdat,headerIndex,forward)

      IMPLICIT NONE
      
      type (struct_obs), intent(inout):: obsdat
      integer, intent(in) :: headerIndex
      logical, intent(in) :: forward

      INTEGER  :: bodyIndex,RLN,NLV

      REAL(OBS_REAL) :: OBSV
      real(OBS_REAL) :: vomp, voma, voer, vhpht, scale
      integer        :: nexp,iobs,iexp

      real(OBS_REAL), allocatable :: expnt(:)

      RLN=obs_headElem_i(obsdat,OBS_RLN,headerIndex)
      NLV=obs_headElem_i(obsdat,OBS_NLV,headerIndex)

      allocate(expnt(nlv))

      !=======================================
      ! Count number of power of 10 exponents
      !=======================================

      nexp=0
      DO bodyIndex = RLN, NLV + RLN -1
         if (obs_bodyElem_i(obsdat,OBS_VNM,bodyIndex) == BUFR_SCALE_EXPONENT) then
            nexp=nexp+1
            expnt(nexp) = obs_bodyElem_r(obsdat,OBS_VAR,bodyIndex)
         end if
      END DO

      if (nexp == 0) then
         deallocate(expnt)
         return
      end if

      if (nexp /= nlv/2) then
         ! Skip over obs assuming mantissa was filtered out in READBURP 
         ! (not inserted in obsSpaceData) due to quality flags.
         ! Set exponent quality flag to that of a 'Suspicious element' 
         
         do bodyIndex = RLN, NLV + RLN -1
            call obs_bodySet_i(obsdat,OBS_FLG,bodyIndex, ibset(obs_bodyElem_i(obsdat,OBS_FLG,bodyIndex),04) )
         end do
              
         ! write(*,*) 'NLV =',nlv,' Nexp=',nexp    
         ! call utl_abort('set_scale_chm: Inconsistent number of exponents')
         deallocate(expnt)
         return
      end if

      if (forward) then
         
         !========================================
         ! Apply power of 10 exponents if present
         !========================================
         
         iobs=0
         DO bodyIndex = RLN, NLV + RLN -1
            IF (obs_bodyElem_i(obsdat,OBS_VNM,bodyIndex) /= BUFR_SCALE_EXPONENT) THEN
               iobs=iobs+1
               obsv=obs_bodyElem_r(obsdat,OBS_VAR,bodyIndex)
               call obs_bodySet_r(obsdat,OBS_VAR,bodyIndex,obsv*10**(expnt(iobs)) )
            end if
         END DO
      
      else
             
         !==========================================
         ! Unapply power of 10 exponents if present
         !==========================================

         iobs=0
         iexp=0
         DO bodyIndex = RLN, NLV + RLN -1
            if (obs_bodyElem_i(obsdat,OBS_VNM,bodyIndex) == BUFR_SCALE_EXPONENT) then
               ! Store scaling exponents
               iexp=iexp+1
               call obs_bodySet_r(obsdat,OBS_OMP,bodyIndex,expnt(iexp))
               call obs_bodySet_r(obsdat,OBS_OMA,bodyIndex,expnt(iexp))
               call obs_bodySet_r(obsdat,OBS_OER,bodyIndex,expnt(iexp))
               call obs_bodySet_r(obsdat,OBS_HPHT,bodyIndex,expnt(iexp))
            else
               iobs=iobs+1
               obsv=obs_bodyElem_r(obsdat,OBS_VAR,bodyIndex)
               vomp=obs_bodyElem_r(obsdat,OBS_OMP,bodyIndex)
               voma=obs_bodyElem_r(obsdat,OBS_OMA,bodyIndex)
               voer=obs_bodyElem_r(obsdat,OBS_OER,bodyIndex)
               vhpht=obs_bodyElem_r(obsdat,OBS_HPHT,bodyIndex)
               scale=10**(-expnt(iobs))
               call obs_bodySet_r(obsdat,OBS_VAR,bodyIndex,obsv*scale )
               call obs_bodySet_r(obsdat,OBS_OMP,bodyIndex,vomp*scale )
               call obs_bodySet_r(obsdat,OBS_OMA,bodyIndex,voma*scale )
               call obs_bodySet_r(obsdat,OBS_OER,bodyIndex,voer*scale )
               call obs_bodySet_r(obsdat,OBS_HPHT,bodyIndex,vhpht*scale )
            end if
         END DO
                 
      end if

      RETURN
  END SUBROUTINE SET_SCALE_CHM


  SUBROUTINE  ADJUST_SFVCOORD(obsdat,START,END)
!
!**s/r ADJUST_SFVCOORD  - Computation of HEIGHT ASSIGNED TO SURFACE OBSERVATIONS
!
!
!Author  : P. Koclas *CMC/CMDA  April 2013
!Revision:
!          S. Macpherson *ARMA  Oct 2013
!              -- add GB-GPS (GP family) element BUFR_NEZD (ele 15031)
!              -- NOTE that for GP data, ELEV = GPS Antenna Height so
!                 no adjustment is needed (SFC_VCO=0).
!
!*    Purpose:  -Compute  HEIGHT ASSIGNED TO SURFACE OBSERVATIONS
!                and INSERT INTO CMA.
!
!
!Arguments
!               INPUT:
!                  -OBSDAT    : instance of obsspace_data module object
!                  -START     : FIRST OBERVATION
!                  -END       : LAST  OBERVATION
!
      IMPLICIT NONE
      INTEGER  :: START,END
      INTEGER  :: J,JO,RLN,NLV
      INTEGER  :: VARNO,CODTYP,ITY
      REAL     :: SFC_VCO,ELEV
      REAL(OBS_REAL) :: PPP
      type (struct_obs), intent(inout):: obsdat
!-----------------------------------------------------------------------
      WRITE(*,*)'   ADJUST_SFVCOORD '
!
!-----------------------------------
!      STN LOOP
!-----------------------------------
!     DO JO=1,obs_numheader(obsdat)
      DO JO=START,END
        RLN=obs_headElem_i(obsdat,OBS_RLN,JO)
        NLV=obs_headElem_i(obsdat,OBS_NLV,JO)
        ITY=obs_headElem_i(obsdat,OBS_ITY,JO)
        CODTYP = ITY
        !=================================
        ! DATA LOOP
        !=================================
        ELEV=obs_headElem_r(obsdat,OBS_ALT,JO)
        DO J = RLN, NLV + RLN -1

          VARNO=obs_bodyElem_i(obsdat,OBS_VNM,j)
          SELECT CASE(VARNO)
            CASE(BUFR_NEDS,BUFR_NEFS,BUFR_NEUS,BUFR_NEVS,BUFR_NETS,BUFR_NESS,BUFR_NEPN,BUFR_NEPS,BUFR_NEHS,BUFR_NEZD)
!           CASE(11011,11012,11215,11216,12004,12203,10051,10004,13220,15031)
!
             SFC_VCO= SURFVCORD(VARNO,CODTYP)
             IF ( VARNO /= BUFR_NEPN) THEN
                PPP=ELEV  + SFC_VCO
                call obs_bodySet_r(obsdat,OBS_PPP,j,PPP)
                call obs_bodySet_i(obsdat,OBS_VCO,j,1)
             ELSE
                 PPP=0.
                call obs_bodySet_r(obsdat,OBS_PPP,j,PPP)
                call obs_bodySet_i(obsdat,OBS_VCO,j,1)
             END IF
          END SELECT
        END DO
        !=================================
!
      END DO
!-----------------------------------
!
      WRITE(*,*)' DONE   ADJUST_SFVCOORD '
      RETURN

  END SUBROUTINE ADJUST_SFVCOORD


  REAL FUNCTION SURFVCORD(ILEM,IDTYP)
!
      implicit none
      INTEGER ILEM,IDTYP,TYPE
!     REAL SURFVCORD
      REAL VCORDSF2
!
!***********************************************************************
!
!      PURPOSE: SEt vertical coordinate for surface data.
!
!       AUTHOR:   P. KOCLAS (CMC/CMDA) December 2011
!
!       Revision : 
!
!    ARGUMENTS:
!               INPUT:
!                      -ILEMP   : BURP ELEMENT NUMBER
!                      -IDTYP   : BURP CODETYPE
!
!               OUTPUT:
!                      -SURFVCORD
!
!
!***********************************************************************
!
!
!     GENERATE TABLES TO ADJUST VERTICAL COORDINATE OF SURFACE DATA
!
!     DEFAULT VALUE 
!    =====================
        vcordsf2=0.
!    =====================

       select case(IDTYP)
        case(135,136,137,138,32,34,35,37,38,159,160,161,162)
!      -----------------
!       UPPER AIR LAND
!      -----------------
        TYPE=3

        case(139,140,141,142,33,36)
!      -----------------
!       UPPER AIR SHIP
!      -----------------
        TYPE=4

        case(12,14,146)
!      -----------------
!       SYNOPS
!      -----------------
        TYPE=1

        case(13,18,145,147)
!      -----------------
!       SHIPS
!      -----------------
        TYPE=2

        case(254)
!      --------------------
!       SCATTEROMETER WINDS
!      --------------------
        TYPE=5

!      -----------------
        case default
!      -----------------
        TYPE=-99
      end select

!

       select case(TYPE)
!===================================================================
         case (1)
          select case(ilem)
!           case (11011,11012,11215,11216)
            case (BUFR_NEDS,BUFR_NEFS,BUFR_NEUS,BUFR_NEVS)

!            us,vs,ffs,dds
!           ==============
            vcordsf2=10.0
!           ==============

!           case (11051)
            case (BUFR_NEPN)
!           pnm
!           ==============
            vcordsf2=0.0
!           ==============

!           ps
!           case (10004)
            case (BUFR_NEPS)
!           ==============
            vcordsf2=0.0
!           ==============

!           ts
!           case (12004)
            case (BUFR_NETS)
!           ==============
            vcordsf2=1.5
!           ==============

!           t-td
!           case (12192,12203)
            case (BUFR_NEES,BUFR_NESS)
!           ==============
            vcordsf2=1.5
!           ==============

        end select
!===================================================================

!===================================================================
         case (2)
          select case(ilem)
!            us,vs,ffs,dds
!           case (11011,11012,11215,11216)
            case (BUFR_NEDS,BUFR_NEFS,BUFR_NEUS,BUFR_NEVS)
!           ==============
            vcordsf2=20.0
!           ==============

!           case (11051)
            case (BUFR_NEPN)
!           pnm
!           ==============
            vcordsf2=0.0
!           ==============

!           case (10004)
            case (BUFR_NEPS)
!           ps
!           ==============
            vcordsf2=0.0
!           ==============

!           case (12004)
            case (BUFR_NETS)
!           ts
!           ==============
            vcordsf2=11.5
!           ==============

!           case (12192,12203)
            case (BUFR_NEES,BUFR_NESS)
!           t-td
!           ==============
            vcordsf2=11.5
!           ==============
        end select
!===================================================================

!===================================================================
         case (3)
            select case(ilem)
!           case (11011,11012,11215,11216)
            case (BUFR_NEDS,BUFR_NEFS,BUFR_NEUS,BUFR_NEVS)
               vcordsf2=10.0

!              case (11051)
               case (BUFR_NEPN)
!              pnm
!           ===============
               vcordsf2=0.0
!           ===============

!              case (10004)
               case (BUFR_NEPS)
!              ps
!           ===============
               vcordsf2=0.0
!           ===============

!              case (12004)
               case (BUFR_NETS)
!              ts
!           ===============
               vcordsf2=1.5
!           ===============

!              case (12192)
               case (BUFR_NEES)
!              t-td
!           ===============
               vcordsf2=0.0
!           ===============

!              t-td(surf)
!              case (12203)
               case (BUFR_NESS)
!           ===============
               vcordsf2=1.5
!           ===============
        end select
!===================================================================

!===================================================================
         case (4)
            select case(ilem)
!           case (11011,11012,11215,11216)
            case (BUFR_NEDS,BUFR_NEFS,BUFR_NEUS,BUFR_NEVS)
!           ===============
            vcordsf2=20.0
!           ===============

!              case (11051)
               case (BUFR_NEPN)
!           pnm
!           ===============
            vcordsf2=0.0
!           ===============

!              case (10004)
               case (BUFR_NEPS)
!           ps
!           ===============
            vcordsf2=0.0
!           ===============

!              case (12004)
               case (BUFR_NETS)
!           ts
!           ===============
            vcordsf2=1.5
!           ===============

!           t-td
!              case (12192)
               case (BUFR_NEES)
!           ===============
            vcordsf2=0.0
!           ===============
!              case (12203)
               case (BUFR_NESS)
!           ===============
            vcordsf2=1.5
!           ===============
        end select
!===================================================================

!===================================================================
         case (5)
            select case(ilem)
!           case (11011,11012,11215,11216)
            case (BUFR_NEDS,BUFR_NEFS,BUFR_NEUS,BUFR_NEVS)
!           ===============
            vcordsf2=10.0
!           ===============
            end select
!===================================================================

         end select

!
!        *******************
         SURFVCORD=VCORDSF2
!        *******************
!
      RETURN
  END FUNCTION  SURFVCORD


  SUBROUTINE FDTOUV_OBSDAT(obsdat,START,END,PPMIS)
!
!---------------------------------------------------------------
!
! Author  : P. Koclas, CMC/CMDA December  2012
!           CONVERT DD , FF  WINDS TO
!            UU (est-west),  VV (north-south) COMPONENTS
!
!    ARGUMENTS:
!                 INPUT:
!
!                       -obsdat     : CMA_table INSTANCE 
!                       -START     : FIRST OBERVATION
!                       -END       : LAST  OBERVATION
!                       -PPMIS     : MISSING VALUE  
!
!        **************************************************
!         IT IS ASSUMED THAT CMA CONTAINS ENTRIES   FOR 
!          UU AND VV  with observed values = missing value
!        **************************************************
!
!---------------------------------------------------------------
!
    implicit none
    type (struct_obs), intent(inout) :: obsdat
    
    REAL*4          :: PPMIS
    INTEGER*4       :: START,END
    INTEGER*4       :: VARNO,VARNO2,VARNO4

    REAL*4          :: OBSUV
    INTEGER*4       :: JO,RLN,NLV,j,j2,j4,Jpos,ilem
    INTEGER*4       :: DDFLAG,FFFLAG,NEWFLAG,UUFLAG,VVFLAG
    INTEGER*4       :: ILEMF,ILEMU,ILEMV,INDU_MISG,INDV_MISG,INDUM,INDVM
    LOGICAL         :: LLMISDD,LLMISFF,LLMIS,LLUV_misg,LLU_misg,LLV_misg
    LOGICAL         :: LLUV_PRESENT,LLU_PRESENT,LLV_PRESENT

    INTEGER         :: NOBSOUT

    REAL(OBS_REAL)  :: UU,VV,DD,FF
    REAL(OBS_REAL)  :: LEVEL_DD,LEVEL4,LEVEL,LEVEL_UU

    NOBSOUT = 6
    FFFLAG = 0  ! bhe 

    !--------------------------------
    !   HEADER LOOP
    !--------------------------------
    HEADER1: do JO=START,END
        
      RLN=obs_headElem_i(obsdat,OBS_RLN,JO)
      NLV=obs_headElem_i(obsdat,OBS_NLV,JO)
      !--------------------------------
      ! TOP DATA LOOP
      !--------------------------------
      DO J = RLN, NLV + RLN -1
        DD = PPMIS
        FF = PPMIS

        VARNO = obs_bodyElem_i(obsdat,OBS_VNM,j)
        LLMISDD =.true.

 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SELECT case (VARNO)
 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        case (BUFR_NEDD,BUFR_NEDS)
          IF( VARNO == BUFR_NEDS) then
            ILEMF=BUFR_NEFS
            ILEMU=BUFR_NEUS
            ILEMV=BUFR_NEVS
          ELSE
            ILEMF=BUFR_NEFF
            ILEMU=BUFR_NEUU
            ILEMV=BUFR_NEVV
          END IF

          DD      =obs_bodyElem_r(obsdat,OBS_VAR,j)
          DDFLAG  =obs_bodyElem_i(obsdat,OBS_FLG,j)
          LEVEL_dd = obs_bodyElem_r(obsdat,OBS_PPP,j)
          
          LLU_misg = .false.
          LLV_misg = .false.
          LLU_PRESENT = .false.
          LLV_PRESENT = .false.
          INDUM = -1
          INDVM = -1
          ! FIND IF  U AND V ARE ALREADY IN CMA
          !-------------------------------------
          uvinobsdat: do J4 =J, NLV + RLN -1
            !-------------------------------------
            LEVEL4 = obs_bodyElem_r(obsdat, OBS_PPP,j4)
            IF (LEVEL4 == LEVEL_dd) then
              VARNO4=obs_bodyElem_i(obsdat, OBS_VNM,j4)
              SELECT case (VARNO4)
              case (11003,11004,11002,11001,11215,11216,11011,11012)

                OBSUV =obs_bodyElem_r(obsdat, OBS_VAR,j4)
                IF (  (VARNO4 == ILEMU)     .and.  (obsuv /= PPMIS) ) THEN
                  LLU_PRESENT=.true.
                  INDUM=J4
                ELSE IF ( (VARNO4 == ILEMV) .and. (obsuv /= PPMIS) ) THEN
                  LLV_PRESENT=.true.
                  INDVM=J4
                END IF
                
                IF (  (VARNO4 == ILEMU)     .and. (obsuv == PPMIS) ) THEN
                  LLU_misg=.true.
                  INDU_MISG=J4
                ELSE IF ( (VARNO4 == ILEMV)  .and. (obsuv == PPMIS) ) THEN
                  LLV_misg=.true.
                  INDV_MISG=J4
                END IF
                
              END SELECT
            END IF

            !-------------------------------------
          end do uvinobsdat
          !-------------------------------------

          LLUV_misg = (LLU_misg .and. LLV_misg)
          LLUV_PRESENT= (LLU_PRESENT .and. LLV_PRESENT)

          !      *******************************
          IF (   LLUV_misg) THEN
            !      *******************************

            !---------------------------------
            calcuv: do J2 =J, NLV + RLN -1
              !---------------------------------

              LLMISFF =.true.
              LLMISDD =.true.
              LLMIS   =.true.
              LEVEL=obs_bodyElem_r(obsdat,OBS_PPP,j2)
              if ( LEVEL /= LEVEL_dd) cycle
              VARNO2=obs_bodyElem_i(obsdat,OBS_VNM,j2)
              !==VARNO2=============================================
              IF (  (VARNO2) == ILEMF ) THEN

                FF   =obs_bodyElem_r(obsdat,OBS_VAR,j2)
                FFFLAG=obs_bodyElem_i(obsdat,OBS_FLG,j2)
                IF (  (DD == 0.  .AND. FF > 0.) .or. ( DD > 360. .OR. DD  < 0.) ) THEN
                  LLMISDD =.true.
                  LLMISFF =.true.
                ELSE IF ( DD == PPMIS .OR. FF == PPMIS)  THEN
                  LLMISDD =.true.
                  LLMISFF =.true.
                ELSE
                  LLMISDD=.false.
                  LLMISFF=.false.
                ENDIF
                !
                !             IF SPEED = 0 CALM WIND IS ASSUMED.
                !             ==================================
                IF (FF == 0.0) THEN
                  DD = 0.
                ENDIF
                   
                DD=DD + 180.
                IF ( DD > 360.) DD=DD-360.
                DD=DD*MPC_RADIANS_PER_DEGREE_R8
                
                !                U,V COMPONENTS ARE
                !==============================================
                UU =FF*SIN(DD)
                VV =FF*COS(DD)
                if  ( ( llmisdd .eqv. .true.) .or. ( llmisff .eqv. .true. ) ) then
                  llmis=.true.
                  if ( INDU_MISG > 0 .or. INDV_MISG > 0 ) then
                    call obs_bodySet_i(obsdat,OBS_VNM,INDU_MISG,-1)
                    call obs_bodySet_i(obsdat,OBS_VNM,INDV_MISG,-1)
                  end if
                else
                  llmis=.false.
                end if

              END IF
              NEWFLAG = IOR(DDFLAG,FFFLAG)

              if ( INDUM > 0 .or. INDVM > 0 ) then
                call obs_bodySet_i(obsdat,OBS_VNM,INDU_MISG,-1)
                call obs_bodySet_i(obsdat,OBS_VNM,INDV_MISG,-1)
              end if
              IF (llmis .eqv. .true.) THEN
                if ( INDUM > 0 .or. INDVM > 0 ) then
                  call obs_bodySet_i(obsdat,OBS_FLG,induM,NEWFLAG)
                  call obs_bodySet_i(obsdat,OBS_FLG,indvM,NEWFLAG)
                end if
              ELSE IF (llmis .eqv. .false.) THEN
                call obs_bodySet_r(obsdat,OBS_VAR,INDU_MISG,UU)
                call obs_bodySet_i(obsdat,OBS_FLG,INDU_MISG,NEWFLAG)

                call obs_bodySet_r(obsdat,OBS_VAR,INDV_MISG,VV)
                call obs_bodySet_i(obsdat,OBS_FLG,INDV_MISG,NEWFLAG)
              END IF
!
              !---------------------
            END DO calcuv
            !---------------------
            !      *******************************
          ELSE                       
            !      *******************************
            IF ( LLUV_PRESENT .eqv. .true. )  THEN
              call obs_bodySet_i(obsdat,OBS_VNM,INDU_MISG,-1)
              call obs_bodySet_i(obsdat,OBS_VNM,INDV_MISG,-1)
            ELSE
              if (indum > 0) then
                call obs_bodySet_i(obsdat,OBS_VNM,indum,-1)
              end if
              if (indvm > 0) then
                call obs_bodySet_i(obsdat,OBS_VNM,indvm,-1)
              end if
            END IF
            !      *******************************
          END IF
          !      *******************************

          !---------------------

          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        END SELECT
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        !--------------------------------
        ! END TOP DATA LOOP
        !--------------------------------
      end do

    end do HEADER1

!==========================================================================================
!  do JO=1,obs_numHeader(obsdat)
   do JO=START,END

    RLN=obs_headElem_i(obsdat,OBS_RLN,JO)
    NLV=obs_headElem_i(obsdat,OBS_NLV,JO)
    !--------------------------------
    !  DATA LOOP
    !--------------------------------
    DO J = RLN, NLV + RLN -1

      LLMISDD =.true.
      VARNO=obs_bodyElem_i(obsdat,OBS_VNM,J)
      LEVEL=obs_bodyElem_r(obsdat,OBS_PPP,J)

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SELECT case (VARNO)
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      case (BUFR_NEUU)
        ILEM=BUFR_NEVV
      case (BUFR_NEUS)
        ILEM=BUFR_NEVS
        
        !         case (BUFR_NEVV)
        !   ILEM=BUFR_NEUU
        !         case (BUFR_NEVS)
        !   ILEM=BUFR_NEUS
      case default
        cycle
        Jpos=0
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      END SELECT
 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      Jpos=-1
      !---------------------------------------------------------------------
      !TRANSFER THE FLAG BITS  FROM ONE WIND COMPONENT TO THE OTHER
      !---------------------------------------------------------------------
      DO J4 = RLN, NLV + RLN -1
        UU      =obs_bodyElem_r(obsdat,OBS_VAR,J4)
        LEVEL_UU=obs_bodyElem_r(obsdat,OBS_PPP,J4)
        Jpos=-1
        
!
        if ( LEVEL_UU == LEVEL .and. UU == PPMIS ) then
          call obs_bodySet_i(obsdat,OBS_VNM,J4,-1)
        end if
!

        if ( LEVEL_UU == LEVEL .and. UU /= PPMIS ) then
          UUFLAG  =obs_bodyElem_i(obsdat,OBS_FLG,J4)
          VARNO2  =obs_bodyElem_i(obsdat,OBS_VNM,J4)
          !            SELECT case (VARNO2)
          !              case (BUFR_NEUU,BUFR_NEUS,BUFR_NEVV,BUFR_NEVS)
          !============================================================
          IF ( (ILEM == VARNO2)  ) THEN
            VVFLAG  =obs_bodyElem_i(obsdat,OBS_FLG,J)
            NEWFLAG =IOR(UUFLAG,VVFLAG)
            call obs_bodySet_i(obsdat,OBS_FLG,J, NEWFLAG)
            call obs_bodySet_i(obsdat,OBS_FLG,J4,NEWFLAG)
            Jpos=J4
            exit
          END IF
          !============================================================
          !            END SELECT

        end if
        !----------------------------------------------------------------
      END DO !J4
      !----------------------------------------------------------------

      !---------------------------------------------------------------------
      !ELIMINATE ENTRIES WHERE ONE COMPONENT OF WIND (UU OR VV) IS MISSING
      !---------------------------------------------------------------------
      if (Jpos < 0) then
        WRITE(*,*) ' eliminate winds for station : ',obs_elem_c(obsdat,'STID',JO),obs_bodyElem_i (obsdat,OBS_VNM,J),obs_bodyElem_r(obsdat,OBS_PPP,J)
        call obs_bodySet_i(obsdat,OBS_VNM,J,-1)
      end if

	!--------------------------------
        END DO !J
	!--------------------------------

        END DO !JO
!==========================================================================================

  END SUBROUTINE FDTOUV_OBSDAT


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

      RETURN
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

  SUBROUTINE SETASSFLG(lobsSpaceData)
!*    Purpose:  -Set BANCO QUALITY CONTROL BIT No 12 FOR ALL DATA ASSIMILATED
!                BY CURRENT ANALYSIS.
!
    IMPLICIT NONE
!
      type(struct_obs) :: lobsSpaceData
      INTEGER INDEX_BODY
!
!     Process all data
!
      DO INDEX_BODY=1,obs_numBody(lobsSpaceData)
         IF (obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY) == 1)  THEN
            call obs_bodySet_i(lobsSpaceData,OBS_FLG,INDEX_BODY,ibset( obs_bodyElem_i(lobsSpaceData,OBS_FLG,INDEX_BODY), 12 ))
         END IF
      END DO
!--------------------------------------------------------------------
    END SUBROUTINE SETASSFLG

!--------------------------------------------------------------------------
!! *Purpose*: Retrieves information for observations from all BURP files. Currently only used for
!!            chemical constituent but can potentially be used for any family. If the BURP files
!!            are not split, then the single BURP file is read. If the BURP files are split, then
!!            each split BURP file is read and the data combined into one struct_oss_obsdata
!!            object. See the function burp_chem_read for more details.
!!
!! @author M. Sitwell, ARQI/AQRD, Sept 2016
!!
!! Revisions:
!!v       Y. Rochon, ARQI/AQRD, Nov 2016
!!v         - Added optional input argument codtyplist and option of varno <=0 (in burp_chem_read)
!!
!! Input:
!!v           stnid         station ID of observation
!!v           varno         BUFR code 
!!v                         If <=0, to search through all codes to obtain first
!!v                         between 10000 and 16000.
!!v           nlev          number of levels in the observation
!!v           ndim          number of dimensions for the retrieved data in
!!v                         each report (e.g. ndim=1 for std, ndim=2 for
!!v                         averagine kernels) 
!!v           bkstp         bkstp number of requested block
!!v           block_type    block type indicated by the two rightmost bits
!!v                         of bknat. Valid values are 'DATA', 'INFO', '3-D',
!!v                         and 'MRQR'.
!!v           match_nlev    determines if the report matching criteria includes checking
!!v                         if the report number of levels is the same as the input
!!v                         argument nlev
!!v           codtyplist    optional CODTYP list for search (optional)
!!v           obsfam        observation family name
!!
!! Output:
!!v           burp_out      struct_oss_obsdata object
!!
!! Comment:
!!
!!v   This routine is general enough to be used by observation families
!!v   other than 'CH', It should be renamed once used for other families.
!!
!--------------------------------------------------------------------------
  function burp_chem_read_all(obsfam,stnid,varno,nlev,ndim,bkstp,block_type,match_nlev,  &
                              codtyplist_opt) result(burp_out)

    implicit none

    character(len=9), intent(in) :: stnid
    character(len=4), intent(in) :: block_type
    integer, intent(in)          :: ndim,varno,nlev,bkstp
    logical, intent(in)          :: match_nlev
    integer, intent(in), optional :: codtyplist_opt(:)
    character(len=*), intent(in) :: obsfam
    type(struct_oss_obsdata) :: burp_out

    character(len=500) :: filename
    logical :: found

    filename = burp_get_filename(obsfam,found)

    if (found) then
       burp_out = burp_chem_read(filename,stnid,varno,nlev,ndim,bkstp,block_type,match_nlev, &
                                 codtyplist_opt=codtyplist_opt)
    else
       if (burp_split_L) then
          ! Must allocate burp_out so that it is available from ALL processors when
          ! requiring of rpn_comm_allgather via oss_obsdata_MPIallgather.
          write(*,*) "burp_chem_read_all: Could not find/open BURP file: ",trim(filename)
          if (ndim == 1) then
             call oss_obsdata_alloc(burp_out,1,dim1=nlev)
          else
             call oss_obsdata_alloc(burp_out,1,dim1=nlev,dim2_opt=nlev)
          end if
          burp_out%nrep=0
          write(*,*) "burp_chem_read_all: Number of reports set to ",burp_out%nrep
       else
          call utl_abort('burp_chem_read_all: Could not find BURP file: ' // trim(filename))
       end if
    end if

    if (burp_split_L) call oss_obsdata_MPIallgather(burp_out)

  end function burp_chem_read_all

!--------------------------------------------------------------------------
!! *Purpose*: Retrieve information from observation BURP file. Can retrieve
!!            either 1D or 2D data from a report. Currently only used for
!!            chemical constituent but can potentially be used for any family.
!!
!! @author M. Sitwell, ARQI/AQRD, March 2015
!!
!! Revision:
!!v    M. Sitwell, ARQI/AQRD, May 2015
!!v      - Modified to read both 1D and 2D data from a report
!!v    Y. Rochon, ARQI/AQRD, May 2016
!!v      - Updated to increment 'icount' only if varno is also found in addition
!!v        stnid, nlev, block_type and bkstp
!!v    M. Sitwell, ARQI/AQRD, June 2016
!!v      - Added match_nlev input argument
!!v    Y. Rochon, ARQI/AQRD, Oct 2016
!!v      - Added optional codtyplist argument and option of input varno<=0.
!!
!! Input:
!!v           filename      BURP file name
!!v           stnid         station ID of observation
!!v           varno         BUFR code
!!v                         If <=0, search through all codes to obtain first
!!v                         between 10000 and 16000.
!!v           nlev          number of levels in the observation
!!v           ndim          number of dimensions for the retrieved data in
!!v                         each report (e.g. ndim=1 for std, ndim=2 for
!!v                         averagine kernels) 
!!v           bkstp         bkstp number of requested block
!!v           block_type    block type indicated by the two rightmost bits
!!v                         of bknat. Valid values are 'DATA', 'INFO', '3-D',
!!v                         and 'MRQR'.
!!v           match_nlev    determines if the report matching criteria includes checking
!!v                         if the report number of levels is the same as the input
!!v                         argument nlev
!!v           codtyplist    optional CODTYP list for search
!!
!! Output: 
!!v           burp_out      struct_oss_obsdata object
!!
!! Comments:
!!
!!v   - BUFR power 10 exponent element (i.e. data with BUFR number BUFR_SCALE_EXPONENT) 
!!v     will only be applied to 1D data if present.
!!v   - As burp_out is for a specific input stnid, burp_out%code only contains the (lat/long and
!!v     time coord.) with 22 characters.
!!v   - Exponent BUFR data (i.e. data with BUFR number BUFR_SCALE_EXPONENT) will only 
!!v     be applied to 1D data.
!!v   - This routine is general enough to be used by observation families
!!v     other than 'CH', It should be renamed once used for other families.
!!
!--------------------------------------------------------------------------
  function burp_chem_read(filename,stnid,varno,nlev,ndim,bkstp,block_type,match_nlev,  &
                          codtyplist_opt) result(burp_out)
    
    implicit none

    character(len=*), intent(in) :: filename
    character(len=9), intent(in) :: stnid
    character(len=4), intent(in) :: block_type
    integer, intent(in)          :: ndim,varno,nlev,bkstp
    logical, intent(in)          :: match_nlev
    integer, intent(in), optional :: codtyplist_opt(:)
    type(struct_oss_obsdata) :: burp_out

    character(len=9)  :: rep_stnid
    type(burp_file)   :: brp
    type(burp_rpt)    :: rep
    type(burp_block)  :: blk
    integer           :: error,ref_rpt,nrep,ref_blk,varno_ivar
    integer           :: ref_bkstp,nval,ivar,iexp,ilev,icount,icodtyp
    integer           :: date,time,ilat,ilon,iele,nele,icol
    real(8)           :: val,exponent

    ! initialize burp file, report, and block
    call BURP_Init(brp, iostat=error)
    call BURP_Init(rep, iostat=error)
    call BURP_Init(blk, iostat=error)

    ! open the burp file
    call BURP_New(brp, FILENAME=filename, MODE=FILE_ACC_READ, IOSTAT=error)
    
    if (error == 0) then
       write(*,*) "burp_chem_read: Reading file " // trim(filename)
       write(*,*) "burp_chem_read: Selecting STNID = ",stnid," BUFR = ",varno," block type = ",block_type
       write(*,*) "burp_chem_read:           bkstp = ",bkstp," nlev = ",nlev," match_nlev = ",match_nlev
       if (present(codtyplist_opt)) write(*,*) "burp_chem_read: CodeTypeList: ",codtyplist_opt(:)
    else
       call utl_abort('burp_chem_read: Could not find/open BURP file: ' // trim(filename))
    end if

    ! get number of reports in file
    call BURP_Get_Property(brp, NRPTS=nrep)

    ! allocate memory
    if (ndim == 1) then
       call oss_obsdata_alloc(burp_out,nrep,dim1=nlev)
    else
       call oss_obsdata_alloc(burp_out,nrep,dim1=nlev,dim2_opt=nlev)
    end if
    
    icount = 0  ! counter of reports with same stnid, number of levels, and varno as input 
    ref_rpt = 0
    
    ! loop through reports    
    REPORTS: do

       ref_rpt = BURP_Find_Report(brp, REPORT=rep, SEARCH_FROM=ref_rpt, IOSTAT=error)

       if (ref_rpt<0) exit REPORTS
       
       call BURP_Get_Property(rep, STNID=rep_stnid, DATE=date, TEMPS=time, LATI=ilat, LONG=ilon, IDTYP=icodtyp) 

       if (present(codtyplist_opt)) then
          if (.not.any(codtyplist_opt(:) == icodtyp)) cycle REPORTS
       end if

       if (.not.utl_stnid_equal(stnid,rep_stnid)) cycle REPORTS

       ! loop through blocks
       ref_blk = 0
       BLOCKS: do
          
          ref_blk = BURP_Find_Block(rep, BLOCK=blk, SEARCH_FROM=ref_blk, IOSTAT=error)          
          if (ref_blk<0) exit BLOCKS
          
          call BURP_Get_Property(blk, NELE=nele, NVAL=nval, BKSTP=ref_bkstp, IOSTAT=error)

          if (.not.IS_Burp_Btyp(trim(block_type),BLOCK=blk) .or. bkstp /= ref_bkstp .or. (match_nlev.and.nval /= nlev)) cycle BLOCKS

          if (varno > 0) then
             ivar = BURP_Find_Element(blk, ELEMENT=varno, IOSTAT=error)
             if (ivar < 0) cycle BLOCKS
          else 
             ! Search for first data element within elements 10000 and 16000.
             varno_ivar=-1
             do ivar=1,nele
                varno_ivar=BURP_Get_Element(blk, INDEX=ivar, IOSTAT=error)
                if (varno_ivar >= 10000.and.varno_ivar < 16000) exit
             end do
             if (varno_ivar < 10000.or.varno_ivar >= 16000) call utl_abort('burp_chem_read: No valid element found for STNID ' // rep_stnid )
          end if

          ! required block found if code reaches this point, retrieve data and store in burp_out
          
          if (nval > nlev) call utl_abort('burp_chem_read: number of levels in the report (' // trim(utl_str(nval)) // &
                                         ') exceeds the specified maximum number of levels (' // trim(utl_str(nlev)) // &
                                         ') for STNID ' // rep_stnid )

          icount=icount+1
          burp_out%code(icount) = oss_obsdata_get_header_code(ilon,ilat,date,time,rep_stnid)  ! this code is a unique identifier for this report

          if (ndim == 1) then
             ! retrieve 1D data

             iexp = BURP_Find_Element(blk, ELEMENT=BUFR_SCALE_EXPONENT, IOSTAT=error)
                
             if (iexp < 0) then
                ! No exponent found in block
                do ilev=1,nval                   
                   burp_out%data1d(ilev,icount) = BURP_Get_Rval(blk, NELE_IND=ivar, NVAL_IND=ilev, NT_IND=1, IOSTAT=error)                
                end do
             else
                ! Apply exponent
                do ilev=1,nval                   
                   val = BURP_Get_Rval(blk, NELE_IND=ivar, NVAL_IND=ilev, NT_IND=1, IOSTAT=error)                
                   exponent = BURP_Get_Rval(blk, NELE_IND=iexp, NVAL_IND=ilev, NT_IND=1, IOSTAT=error)
                   burp_out%data1d(ilev,icount) = val * 10**exponent
                end do
             end if
                   
          else if (ndim == 2) then
             ! retrieve 2D data

             icol = 0
             do iele=1,nele
                ivar = BURP_Get_Element(blk, INDEX=iele, IOSTAT=error)
                if (ivar == varno) then
                   icol = icol+1
                   do ilev=1,nval                  
                      burp_out%data2d(ilev,icol,icount) = BURP_Get_Rval(blk, NELE_IND=iele, NVAL_IND=ilev, NT_IND=1, IOSTAT=error)
                   end do
                end if
             end do

          end if

          exit BLOCKS
          
       end do BLOCKS
      
    end do REPORTS

    ! resize first dimension of data arrays from length of nrep to icount
    call utl_resize(burp_out%code,icount)
    if (ndim == 1) then
       call utl_resize(burp_out%data1d,nlev,icount)
    else if (ndim == 2) then
       call utl_resize(burp_out%data2d,nlev,nlev,icount)
    end if

    burp_out%nrep = icount

    write(*,*) "burp_chem_read: Reading of file complete. Number of reports found: ",burp_out%nrep
    
    ! deallocate
    Call BURP_Free(brp,iostat=error)
    Call BURP_Free(rep,iostat=error)
    Call BURP_Free(blk,iostat=error)
    
  end function burp_chem_read

!--------------------------------------------------------------------------
!! *Purpose*: Add or modify information from all BURP files. Currently only used for
!!            chemical constituent but can potentially be used for any family.
!!            If BURP files are not split, a single files is updated. If BURP files
!!            are split, then each split file is updated. See the function burp_chem_update
!!            for more details.
!!
!! @author M. Sitwell, ARQI/AQRD, Sept 2016
!!
!! Input:
!!v           obsfam        observation family name
!!v           varno         BUFR descriptors. Number of elements must be 
!!v                         max(1,obsdata%dim2)
!!v           bkstp         bkstp number of requested block
!!v           block_type    block type indicated by the two rightmost bits
!!v                         of bknat. Valid values are 'DATA', 'INFO', '3-D',
!!v                         and 'MRQR'.
!!v           obsdata       Input struct_oss_obsdata object for varno.
!!v           multi         Indicates if intended report are for 'UNI' or 'MULTI' level data (optional)
!!v           fname_prefix  filename prefix, e.g. brpua_ (optional)
!!v
!!
!! Output: 
!!v           nrep_modified Number of modified reports
!
!! Comment:
!!
!!v   This routine is general enough to be used by observation families
!!v   other than 'CH', It should be renamed once used for other families.
!!
!--------------------------------------------------------------------------
  function burp_chem_update_all(obsfam,varno,bkstp,block_type,obsdata,multi_opt) result(nrep_modified)

    implicit none

    character(len=4), intent(in) :: block_type
    character(len=*), intent(in) :: obsfam
    type(struct_oss_obsdata), intent(inout) :: obsdata
    integer, intent(in) :: varno(:),bkstp    
    character(len=*), intent(in), optional :: multi_opt
    integer :: nrep_modified

    integer :: ierr,nrep_modified_global

    if (burp_split_L .or. mpi_myid == 0) then
       nrep_modified = burp_chem_update(burp_get_filename(obsfam),varno,bkstp,block_type,obsdata,multi_opt=multi_opt)
    end if

    if (burp_split_L) then
       call rpn_comm_allreduce(nrep_modified,nrep_modified_global,1,"MPI_INTEGER","MPI_SUM","GRID",ierr)
       nrep_modified = nrep_modified_global
    end if

  end function burp_chem_update_all

!--------------------------------------------------------------------------
!! *Purpose*: Add or modify information from BURP file in existing block and
!!            for specified BUFR descriptor varno(s). Provided data can be either
!!            1D or 2D data. Currently only used for chemical constituent but can
!!            potentially be used for any family.
!!
!! @author Y. Rochon, ARQI/AQRD, June 2016 (partly based on burp_chem_read by M. Sitwell)
!!
!! Revisions:
!!v        M. Sitwell, ARQI/AQRD, Aug 2016
!!v          - Modified to preserve order of reports.
!!v        Y. Rochon, ARQI/AQRD, Jan 2017
!!v          - Added optional check for multi using 'DATA' when block_type='INFO'.
!!
!! Input:
!!v           filename      BURP file name
!!v           varno         BUFR descriptors. Number of elements must be 
!!v                         max(1,obsdata%dim2)
!!v           bkstp         bkstp number of requested block
!!v           block_type    block type indicated by the two rightmost bits
!!v                         of bknat. Valid values are 'DATA', 'INFO', '3-D',
!!v                         and 'MRQR'.
!!v           obsdata       Input struct_oss_obsdata object for varno.
!!v           multi         Indicates if intended report are for 'UNI' or 'MULTI' level data (optional)
!!
!! Output:
!!v           nrep_modified Number of modified reports
!!
!! Comments:
!!v  - Currently assumes that all elements of varno(:) are distinct from each other.
!!v  - In blocks with new data to be added/modified, if the varno already exists in the block, the
!!v    new data will overwrite the existing varno data, otherwise will append the new data
!!v    to the block.
!!v  - The settings for BURP_Write_Block should have ENCODE_BLOCK and CONVERT_BLOCK set to
!!v    .true. in all cases, including when the block has not been modified, due to problems
!!v    that can occur when writing blocks containing negative integers with datyp=4.
!!v  - This routine is general enough to be used by observation families
!!v    other than 'CH', It should be renamed once used for other families.
!!
!--------------------------------------------------------------------------
  function burp_chem_update(filename,varno,bkstp,block_type,obsdata,multi_opt) result(nrep_modified)

    implicit none

    character(len=*), intent(in) :: filename
    character(len=4), intent(in) :: block_type
    type(struct_oss_obsdata), intent(inout) :: obsdata
    integer, intent(in) :: varno(:),bkstp
    
    character(len=*), intent(in), optional :: multi_opt

    integer :: nrep_modified,ncount
    logical :: blk_found
    integer, parameter :: LNMX=100000, code_len=90

    character(len=9)  :: stnid
    character(len=code_len) :: code    ! Must be at least as large as burp_code_len
    type(burp_file)   :: brp
    type(burp_rpt)    :: rep,rep_new
    type(burp_block)  :: blk
    integer           :: error,ref_rpt,nrep,ref_blk,ndim,dim1,dim2
    integer           :: ref_bkstp,nval,ivar,ilev,istat
    integer           :: date,time,ilat,ilon,iele,nele,k
    integer, allocatable :: address(:)
    real(4), allocatable :: new_vals(:,:,:)
    logical, allocatable :: modify(:)
    
    ! Check presence of data to update
    if (obsdata%nrep <= 0) then
       write(*,*) 'burp_chem_update: Skipped due to absence of data to update.'
       return
    end if
    
    ! Identify dimensions for the input data    
    ndim=obsdata%ndim
    dim1=obsdata%dim1
    if (ndim == 1) then
       dim2=1
    else
       dim2=obsdata%dim2
    end if
    
    if (size(varno) < dim2) call utl_abort('burp_chem_update: Number of BUFR elements not sufficient. ' // &
                                          trim(utl_str(size(varno))) // ' vs ' // trim(utl_str(dim2)))

    if (code_len < oss_obsdata_code_len()) call utl_abort('burp_chem_update: Length of code string' &
                                          // ' needs to be increased to ' // trim(utl_str(oss_obsdata_code_len())))
     
    ! initialize burp file, report, and block system resources
    call BURP_Init(brp, iostat=error)
    call BURP_Init(rep, R2=rep_new, iostat=error)
    call BURP_Init(blk, iostat=error)

    ! open the burp file in append mode (to replace or add data in a block)
    call BURP_New(brp, FILENAME=filename, MODE=FILE_ACC_APPEND, IOSTAT=error)
    if (error /= 0) call utl_abort('burp_chem_update: Could not open BURP file: ' // trim(filename))

    ! get number of reports in file
    call BURP_Get_Property(brp, NRPTS=nrep)

    allocate(address(nrep),modify(nrep),new_vals(dim1,dim2,nrep))
    address(:)=0
    modify(:)=.false.
    new_vals(:,:,:)=0.

    ! First loop through reports to identify addresses of original file as well as identify if new
    ! information should be included to that report.
    ! NOTE: The addresses of all reports have to be saved in their original order to ensure the
    !       order of the reports in the file is unchanged.
    ref_rpt=0
    ncount=0
    obsdata%irep=1
    REPORTS1: do

       ref_rpt = BURP_Find_Report(brp, REPORT=rep, SEARCH_FROM=ref_rpt, IOSTAT=error)
       if (ref_rpt<0) exit REPORTS1

       ncount=ncount+1
       address(ncount)=ref_rpt

       call BURP_Get_Property(rep, STNID=stnid, DATE=date, TEMPS=time, LATI=ilat, LONG=ilon)

       if (stnid(1:2) == '>>') cycle REPORTS1

       ! Get unique identifier for search from input data
       code = oss_obsdata_get_header_code(ilon,ilat,date,time,stnid)
       
       ! Determine if replacement/additional data likely present for this report
       if (dim1 == 1.and.dim2 == 1) then
          new_vals(1,1,ncount)=real(oss_obsdata_get_element(obsdata,trim(code),1,stat_opt=istat))
       else if (dim2 == 1) then
          new_vals(:,1,ncount)=real(oss_obsdata_get_array1d(obsdata,trim(code),stat_opt=istat))
       else 
          new_vals(:,:,ncount)=real(oss_obsdata_get_array2d(obsdata,trim(code),stat_opt=istat))
       end if

       if (istat == 0) then
          if (present(multi_opt)) then
             ! loop through blocks to find first data block
             ref_blk = 0
             BLOCKS1: do
                ref_blk = BURP_Find_Block(rep, BLOCK=blk, SEARCH_FROM=ref_blk, IOSTAT=error)          
                if (ref_blk<0) exit BLOCKS1
                if (IS_Burp_Btyp('DATA',BLOCK=blk)) then
                   if (IS_Burp_Btyp(trim(multi_opt),BLOCK=blk)) modify(ncount) = .true.
                   exit BLOCKS1
                end if
             end do BLOCKS1
          else
             modify(ncount) = .true.
          end if
       end if

    end do REPORTS1
    
    nrep_modified = count(modify)   ! number of reports with same code and, possibly, same number of obs data levels

    ! Generate new report
    Call BURP_New(rep_new, ALLOC_SPACE=10*LNMX, IOSTAT=error)

    ! second loop through reports to include the new information to the file    
    REPORTS2: do k=1,ncount
    
       call BURP_Get_Report(brp, REPORT=rep, REF=address(k), IOSTAT=error)
       
       ! Copy report header
       Call BURP_Copy_Header(TO=rep_new,FROM=rep)
       Call BURP_Init_Report_Write(brp,rep_new,IOSTAT=error)

       ! loop through blocks
       ref_blk = 0
       BLOCKS: do
          
          ref_blk = BURP_Find_Block(rep, BLOCK=blk, SEARCH_FROM=ref_blk, IOSTAT=error)          
          if (ref_blk<0) exit BLOCKS
          
          if (modify(k)) then

             call BURP_Get_Property(blk, NELE=nele, NVAL=nval, BKSTP=ref_bkstp, IOSTAT=error)

             blk_found = IS_Burp_Btyp(trim(block_type),BLOCK=blk) .and. bkstp == ref_bkstp .and. dim1 == nval
         
             if (blk_found) then
                ! Block to be modified has been found, add new data to block.
                ! If the varno is already in the block, the new data will overwrite the
                ! existing data, otherwise will append the new data to the block.

                do iele=1,dim2
                   ivar = BURP_Find_Element(blk, ELEMENT=varno(iele), IOSTAT=error)           
                   if (ivar < 0) then
                      ivar=nele+1
                      call BURP_Resize_Block(blk,ADD_NELE=1,IOSTAT=error)
                      call BURP_Set_Element(blk,NELE_IND=ivar,ELEMENT=varno(iele),IOSTAT=error)
                   end if
                
                   do ilev=1,nval 
                      call BURP_Set_Rval(blk,NELE_IND=ivar,NVAL_IND=ilev,NT_IND=1,RVAL=new_vals(ilev,iele,k),IOSTAT=error)                 
                   end do
                end do
        
             end if
          end if
             
          call BURP_Write_Block(rep_new, BLOCK=blk, ENCODE_BLOCK=.true., CONVERT_BLOCK=.true., IOSTAT=error)
         
       end do BLOCKS
       
       call BURP_Delete_Report(brp,rep,IOSTAT=error)
       call BURP_Write_Report(brp,rep_new,IOSTAT=error) 
  
    end do REPORTS2
        
    ! deallocate
    deallocate(address,modify,new_vals)
    Call BURP_Free(brp,iostat=error)
    Call BURP_Free(rep,R2=rep_new,iostat=error)
    Call BURP_Free(blk,iostat=error)
    
  end function burp_chem_update

!--------------------------------------------------------------------------
!! *Purpose*: Returns the BURP file name assigned to the calling processor. If the
!!            input family has more than one file, the first file found will be
!!            returned. File names are assigned in the module burpFiles_mod.
!!
!! @author M. Sitwell  Sept 2016
!!
!! Input:
!!v    obsfam            observation family name
!!
!! Output:
!!v    burp_filename  file name of associated BURP file
!!v    found          logical indicating if the BURP file could be found (optional)
!--------------------------------------------------------------------------
  function burp_get_filename(obsfam,found_opt) result(burp_filename)

    implicit none

    character(len=2), intent(in) :: obsfam
    logical, intent(out), optional :: found_opt
    character(len=128) :: burp_filename
    
    logical :: file_found
    integer :: ifile

    burp_filename = ""
    file_found = .false.
       
    do ifile=1,burp_nfiles
       if (obsfam == burp_cfamtyp(ifile)) then
          burp_filename = burp_cfilnam(ifile)
          inquire(file=trim(burp_filename), exist=file_found)
          exit
       end if
    end do

    if (.not.file_found) write(*,*) "burp_get_filename: File not found for observation family " // trim(obsfam)

    if (present(found_opt)) found_opt = file_found

  end function burp_get_filename

end module burpFiles_mod
