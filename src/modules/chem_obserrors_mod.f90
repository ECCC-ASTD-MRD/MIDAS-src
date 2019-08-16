!-------------------------------------- LICENCE BEGIN ------------------------------------
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

module chem_obserrors_mod
  ! MODULE chem_obserrors_mod (prefix='chm' category='5. B and R matrices')
  !
  ! :Purpose: Holds observation standard-deviation information for chemical
  !           species
  !

  use utilities_mod
  use obsSubSpaceData_mod
  use burpFiles_mod    
  use obsFiles_mod
  use MathPhysConstants_mod
  use chem_setup_mod, only: chm_setup_get_str

  implicit none
  private
  
! public procedures
! -----------------

  public :: chm_dealloc_obs_err_stddev,chm_read_obs_err_stddev, &
            chm_get_obs_err_stddev

! module structures
! -----------------

  type :: struct_chm_std
     !
     ! Structure containing information retrieved from auxiliary file for holding 
     ! observation std dev information
     !
     !  Variable               Description
     !  --------               -----------
     !  n_stnid                Number of sub-families (identified via STNIDs)
     !  stnids                 Sub-families (STNIDs; * are wild cards)
     !  element                BUFR element in data block 
     !  source                 0: Set entirely from the auxiliary file being read. No 
     !                            initial values read from observation files
     !                         1: Initial values in observation files 
     !                            (may be adjusted after input)
     !                         2: Initial values in observation files for variable number
     !                            of vertical levels (for error std deviations only)
     !  std_type               Index of setup approach (used in combination with source)
     !                         For source value 0 or 1, 
     !                         0: std1 or observation file values (sigma)
     !                         1: max(std3,std2*ZVAL)  if source=0
     !                            max(std3,std2*sigma) otherwise
     !                         2: sqrt(std3**2+(std2*ZVAL)**2))  if source=0
     !                            sqrt(std3**2+(std2*sigma)**2)) otherwise
     !                         3: min(std3,max(std2,std1_chm*ZVAL)) if source=0
     !                            min(std3,max(std2,std1_chm*sigma))  otherwise
     !                         4: sqrt(std2**2+(std1*ZVAL)**2))  if source=0 
     !                            sqrt(std2**2+(std1*sigma)**2)) otherwise
     !  ibegin                 Position index of start of data for given
     !                         sub-family in the arrays std1,levels,lat
     !  n_lvl                  Number of vertical levels (max number when source=2)
     !  levels                 Vertical levels (in coordinate of sub-family data)
     !  n_lat                  Number of latitudes
     !  lat                    Latitudes (degrees; ordered in increasing size)
     !  std1                   See std_type for usage (dependent on vertical level)
     !  std2                   See std_type for usage (independent of vertical level)
     !  std3                   See std_type for usage (independent of vertical level)

     integer ::  n_stnid
     character(len=12), allocatable :: stnids(:)
     integer, allocatable :: element(:),std_type(:),n_lat(:)
     integer, allocatable :: source(:),ibegin(:),n_lvl(:)
     real(8), allocatable :: std1(:),std2(:),std3(:)
     real(8), allocatable :: levels(:),lat(:)

  end type struct_chm_std

  type(struct_chm_std)  :: chm_std
 
! Array of pointers to hold std's read from observation files.
! Note: Ideally, these should be an element in the 'struct_chm_std' derived 
! types, but currently this result in an internal compiler error.

  type(struct_oss_obsdata), allocatable :: chm_obsSub_std(:)

contains

  subroutine chm_read_obs_err_stddev
    !
    !:Purpose: To read observation errors from auxiliary file or observation
    !          file.
    !
    implicit none

    integer, parameter :: ndim=1
    integer :: istnid

    ! read the error std. dev. information from the auxiliary file
    call chm_read_obs_err_stddev_auxfile

    ! set size of observation sub-space std array
    allocate(chm_obsSub_std(chm_std%n_stnid))

    ! read from the observation file if requested
    do istnid=1,chm_std%n_stnid
       if (chm_std%source(istnid).ge.1) then
          
          ! retrieve data from stats blocks (with bkstp=14 and block_type='DATA')
          chm_obsSub_std(istnid) = obsf_obsSub_read('CH',chm_std%stnids(istnid),chm_std%element(istnid), &
                                 chm_std%n_lvl(istnid),ndim,bkstp_opt=14,block_opt='DATA', &
                                 match_nlev_opt=chm_std%source(istnid).eq.1 )

       end if
    end do

  end subroutine chm_read_obs_err_stddev

  subroutine chm_read_obs_err_stddev_auxfile
    !
    !:Purpose:  To read and store observation error std. dev. as needed for CH
    !           family obs.
    !
  implicit none

  integer :: FNOM, FCLOS
  integer :: IERR, JLEV, JELM, nulstat, ios, isize, icount
  logical :: LnewExists
  
  character (len=128) :: ligne

  EXTERNAL FNOM,FCLOS

! Initialization

  chm_std%n_stnid=0
!
! CHECK THE EXISTENCE OF THE NEW FILE WITH STATISTICS
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  INQUIRE(FILE=trim(chm_setup_get_str('auxfile')),EXIST=LnewExists)
  IF (.not.LnewExists) then
    WRITE(*,*) '---------------------------------------------------------------'
    WRITE(*,*) 'WARNING! chm_read_obs_err_stddev: auxiliary file ' // trim(chm_setup_get_str('auxfile'))
    WRITE(*,*) 'WARNING! not available. Default CH family stddev to be applied if needed.'
    WRITE(*,*) '---------------------------------------------------------------'
    return
  ENDIF
!
! Read observation error std dev. from file auxiliary file for constituent data
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  NULSTAT=0
  IERR=FNOM(NULSTAT,trim(chm_setup_get_str('auxfile')),'SEQ',0)
  IF ( IERR .EQ. 0 ) THEN
    open(unit=nulstat, file=trim(chm_setup_get_str('auxfile')), status='OLD')
  ELSE
    CALL utl_abort('chm_read_obs_err_stddev_auxfile: COULD NOT OPEN AUXILIARY FILE ' // trim(chm_setup_get_str('auxfile')))
  ENDIF

! Read error standard deviations for constituents if available.
! (CH family; ozone and others)
  
  ios=0
  read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne
  do while (trim(adjustl(ligne(1:12))).ne.'SECTION I:') 
      read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne
  end do    
  
! Read number of observation set sub-families (STNIDs and ...) and allocate space

  read(nulstat,*,iostat=ios,err=10,end=10) chm_std%n_stnid
  read(nulstat,*,iostat=ios,err=10,end=10) isize
  
  allocate(chm_std%stnids(chm_std%n_stnid))
  allocate(chm_std%std_type(chm_std%n_stnid),chm_std%n_lat(chm_std%n_stnid))
  allocate(chm_std%source(chm_std%n_stnid),chm_std%ibegin(chm_std%n_stnid))
  allocate(chm_std%element(chm_std%n_stnid),chm_std%n_lvl(chm_std%n_stnid))
  allocate(chm_std%std1(isize),chm_std%std2(chm_std%n_stnid),chm_std%std3(chm_std%n_stnid))
  allocate(chm_std%levels(isize),chm_std%lat(isize))
 
  chm_std%element(:)=0
  chm_std%source(:)=0
  chm_std%std_type(:)=0
  chm_std%n_lvl(:)=1
  chm_std%n_lat(:)=1

! Begin reading for each sub-family
! Important: Combination of STNID, BUFR element and number of vertical levels
!            to determine association to the observations.

  icount=0
  STNIDLOOP: do jelm=1,chm_std%n_stnid
    chm_std%ibegin(jelm)=icount+1

    ! disregard line of dashes
    read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne

    ! Read STNID (* as wildcard)    
    read(nulstat,'(2X,A9)',iostat=ios,err=10,end=10) chm_std%stnids(jelm) 
    
!   Read (1) BUFR element,
!        (2) Flag indication if EOR provided from this auxiliary file or
!            to be read from the observation file,
!        (3) Index specifying OER setup method,
!        (4) Number of vertical levels
!        (5) Number of latitudes
!
!   Important: Combination of STNID, BUFR element and number of vertical levels
!              to determine association to the observations.
!
    read(nulstat,*,iostat=ios,err=10,end=10) chm_std%element(jelm),chm_std%source(jelm),  &
       chm_std%std_type(jelm), chm_std%n_lvl(jelm), chm_std%n_lat(jelm),  &
       chm_std%std2(jelm), chm_std%std3(jelm)

    if (chm_std%n_lvl(jelm).lt.1) chm_std%n_lvl(jelm)=1
    if (chm_std%n_lat(jelm).lt.1) chm_std%n_lat(jelm)=1
    
    if (icount+chm_std%n_lvl(jelm)*chm_std%n_lat(jelm).gt.isize) then
       write(*,'(10X,"Max array size exceeded: ",I6)') isize
       CALL utl_abort('chm_read_obs_err_stddev_auxfile: PROBLEM READING OBSERR STD DEV.')    
    end if

    ! disregard line of dashes
    read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne
    
    ! disregard data section if not needed
    if (chm_std%std_type(jelm).eq.1.or.chm_std%std_type(jelm).eq.2.or.(chm_std%source(jelm).ge.1.and.chm_std%std_type(jelm).eq.0)) &
         cycle STNIDLOOP 

    if (chm_std%n_lvl(jelm).eq.1.and.chm_std%n_lat(jelm).eq.1) then
    
!      Read one value only (independent of level and latitude)
       
       icount=icount+1
       read(nulstat,*,iostat=ios,err=10,end=10) chm_std%std1(icount)

    else if (chm_std%n_lvl(jelm).eq.1.and.chm_std%n_lat(jelm).gt.1) then
    
!      Value dependent on latitude only
       
!      Read reference latitudes (must be in order of increasing size)
       
       read(nulstat,*,iostat=ios,err=10,end=10)                      &
              chm_std%lat(icount+1:icount+chm_std%n_lat(jelm))
      
!      Read OER-related values
  
       read(nulstat,*,iostat=ios,err=10,end=10)                 &
                 chm_std%std1(icount+1:icount+chm_std%n_lat(jelm))

       icount=icount+chm_std%n_lat(jelm)

    else if (chm_std%n_lvl(jelm).gt.1.and.chm_std%n_lat(jelm).eq.1) then
    
!      Value dependent on vertical level only
      
       do jlev=1,chm_std%n_lvl(jelm)
          icount=icount+1
          
!         Read vertical level and OER-related value.
          
          read(nulstat,*,iostat=ios,err=10,end=10)                 &
                 chm_std%levels(icount),chm_std%std1(icount)

       end do
   
    else if (chm_std%n_lvl(jelm).gt.1.and.chm_std%n_lat(jelm).gt.1) then
    
!      Value dependent on vertical level and latitude 
       
!      Read reference latitudes (must be in order of increasing size)
       read(nulstat,*,iostat=ios,err=10,end=10)                      &
              chm_std%lat(icount+1:icount+chm_std%n_lat(jelm))
!       write(*, '(10X,20F9.3)') chm_std%lat(icount+1:icount+chm_std%n_lat(jelm))
      
       do jlev=1,chm_std%n_lvl(jelm)
          
!         Read vertical level and OER-related lat-dependent values.
          
          read(nulstat,*,iostat=ios,err=10,end=10)                   &
                 chm_std%levels(icount+jlev),                           &
                 chm_std%std1(icount+(jlev-1)*chm_std%n_lat(jelm)+1:icount+jlev*chm_std%n_lat(jelm))

       end do
       icount=icount+chm_std%n_lat(jelm)*chm_std%n_lvl(jelm)
    end if
 end do STNIDLOOP
   
 10 if (ios.gt.0) then
       WRITE(*,*) 'File read error message number: ',ios
       CALL utl_abort('chm_read_obs_err_stddev_auxfile: PROBLEM READING OBSERR STD DEV.')    
    end if
 
 11 CLOSE(UNIT=NULSTAT)
    IERR=FCLOS(NULSTAT)    

  end subroutine chm_read_obs_err_stddev_auxfile

  subroutine chm_obs_err_stddev_index(CSTNID,NLEV,VARNO,ZLAT,ISTNID,JINT)
    ! 
    !:Purpose: To return the station ID and latitude indices corresponding to a
    !          measurement.
    !
    implicit none

    character(len=*), intent(in)  :: CSTNID
    integer, intent(in)           :: NLEV,VARNO
    real(8), intent(in)           :: ZLAT
    integer, intent(out)          :: ISTNID,JINT
    integer                       :: JN,ibegin
    real(8)                       :: lat

 !  Important: Combination of STNID, BUFR element and number of vertical levels
 !             to determine association to the observations.

 !             Find stnid with same number of vertical levels and same BUFR element.
 !             Note: * in chm_std%stnids stands for a wildcard
     
    ISTNID=0
    DO JN=1,chm_std%n_stnid

       ! First compare STNID values allowing for * and blanks in 
       ! chm_std%stnids(JN) as wildcards

       IF (utl_stnid_equal(chm_std%stnids(JN),CSTNID)) THEN
          IF ( (NLEV.EQ.chm_std%n_lvl(JN) .OR. chm_std%source(JN).eq.2) .AND. VARNO.EQ.chm_std%element(JN) ) THEN
             ISTNID=JN
             exit
          END IF
       END IF

    END DO

    IF (ISTNID.EQ.0) THEN
       write(*,*) 'chm_obs_err_stddev_index: Error std. dev. is unavailable for STNID ' // trim(CSTNID) // & 
                  ' and NLEV = ' // trim(utl_str(NLEV)) // ' and VARNO = ' // trim(utl_str(VARNO))
       write(*,*)
       write(*,*) ' Contents of chm_std (n_stnid = ' // trim(utl_str(chm_std%n_stnid)) // '):'
       if (chm_std%n_stnid.gt.0) then
          write(*,'(A)') '  STNID                 BUFR     NLEVELS'
          do jn=1,chm_std%n_stnid
             write(*,'(2X,A,2X,I12,2X,I10)') chm_std%stnids(jn),chm_std%element(jn),chm_std%n_lvl(jn)
          end do
       end if
       write(*,*)
       call utl_abort('chm_obs_err_stddev_index: Check section I of the auxiliary file.')
    ELSE

       IF (chm_std%n_lat(ISTNID) .GT. 1) THEN

          ! Find latitude index for interpolation.
          ! Assuming increasing latitudes in chm_std%lat

          lat = zlat / MPC_RADIANS_PER_DEGREE_R8  ! radians to degrees

          ibegin=chm_std%ibegin(ISTNID)-1
          IF (lat .GE. chm_std%lat(ibegin+chm_std%n_lat(ISTNID))) THEN
             JINT=chm_std%n_lat(ISTNID)+1
          ELSE
             DO JINT=1,chm_std%n_lat(ISTNID)
                IF (lat .LE. chm_std%lat(ibegin+JINT) ) exit
             END DO
          END IF
                                           
       END IF       
    END IF         

  end subroutine chm_obs_err_stddev_index

  function chm_get_obs_err_stddev(cstnid,nlev,varno,zlat,zlon,idate,itime,zval,&
                                  zlev,ilev,ifirst) result(obs_err_stddev) 
    ! 
    !:Purpose: To return the observational error std dev for a CH family
    !          measurement
    implicit none
    real(8)  :: obs_err_stddev 
   
    ! Arguments:
    character(len=*), intent(in) :: CSTNID ! station ID
    integer, intent(in) :: NLEV ! number of levels
    integer, intent(in) :: VARNO ! BUFR number
    real(8), intent(in) :: ZLAT ! latitude (radians)
    real(8), intent(in) :: ZLON ! longitude (radians)
    integer, intent(in) :: IDATE ! date in YYYYMMDD format
    integer, intent(in) :: ITIME ! time in HHMM format
    real(8), intent(in) :: ZVAL  ! observation values
    real(8), intent(in) :: ZLEV  ! vertical coordinate value
    integer, intent(in) :: ILEV  ! observation number in the profile
    logical, intent(in) :: IFIRST! true:  first call for a profile

    ! Locals:
    real(8) :: wgt,zwb,sigma
    integer :: ibegin,JLEV,JN,stat

    integer, save :: ISTNID,JINT
    
    ! If this call is for the first level for this measurement, get
    ! the station ID and latitude indices corresponding to this measurement
    if (ifirst) call chm_obs_err_stddev_index(CSTNID,NLEV,VARNO,ZLAT,ISTNID,JINT)                  
            
    ! Get weighting of error std. dev. if required

    if (chm_std%std_type(ISTNID).gt.2 .or. &
       (chm_std%source(ISTNID).eq.0 .and. chm_std%std_type(ISTNID).eq.0) ) then

       IF (chm_std%n_lvl(ISTNID) .GT. 1) THEN
                 
          ! Find nearest vertical level (no interpolation)
                 
          zwb=1.E10
          ibegin=chm_std%ibegin(ISTNID)-1
          DO JN=1,chm_std%n_lvl(ISTNID)
             IF (zwb .GT. abs(ZLEV-chm_std%levels(ibegin+JN)) ) THEN
                JLEV=JN
                zwb=abs(ZLEV-chm_std%levels(ibegin+JN))
             END IF
          END DO
          JLEV=ibegin+(JLEV-1)*chm_std%n_lat(ISTNID)+1
       ELSE
          JLEV=chm_std%ibegin(ISTNID)
       END IF

       IF (chm_std%n_lat(ISTNID) .GT. 1) THEN
                
          ! Apply interpolation

          JLEV=JLEV+JINT-1
          ibegin=chm_std%ibegin(ISTNID)-1
          IF (JINT.EQ.1.OR.JINT.GT.chm_std%n_lat(ISTNID)) THEN
             wgt=chm_std%std1(JLEV)
          ELSE
             wgt=(chm_std%std1(JLEV-1)*(chm_std%lat(ibegin+JINT)-ZLAT)+ &
                  chm_std%std1(JLEV)*(ZLAT-chm_std%lat(ibegin+JINT-1)))/ &
                  (chm_std%lat(ibegin+JINT)-chm_std%lat(ibegin+JINT-1))
          END IF
       ELSE
          wgt=chm_std%std1(JLEV)             
       END IF
         
    end if
             
    ! Set the error std. dev.
                   
    IF (chm_std%source(ISTNID).EQ.0) THEN
               
       ! Set error standard deviations from scratch using content of
       ! previously read content of the auxiliary file.
                
       select case(chm_std%std_type(ISTNID))
       case(0)
          obs_err_stddev = wgt
       case(1)
          obs_err_stddev = max(chm_std%std3(ISTNID),chm_std%std2(ISTNID)*ZVAL)
       case(2)
          obs_err_stddev = sqrt(chm_std%std3(ISTNID)**2+(chm_std%std2(ISTNID)*ZVAL)**2)
       case(3)
          obs_err_stddev = min(chm_std%std3(ISTNID),max(chm_std%std2(ISTNID),wgt*ZVAL))
       case(4)
          obs_err_stddev = sqrt(chm_std%std2(ISTNID)**2+(wgt*ZVAL)**2)
       case default
          call utl_abort('chm_get_obs_err_stddev: std_type = ' // trim(utl_str(chm_std%std_type(ISTNID))) // &
               ' for STNID = ' // trim(CSTNID) // ' is not recognized.')
       end select

    ELSE
          
       ! Adjust error standard deviations read from observation file if requested.
       
       sigma = oss_obsdata_get_element(chm_obsSub_std(istnid), oss_obsdata_get_header_code(zlon,zlat,idate,itime,cstnid), ilev, stat_opt=stat)

       select case(stat)
       case(1)
          call utl_abort("chm_get_obs_err_stddev: No reports available for STNID = " // trim(cstnid) // &
                       ", nlev = " // trim(utl_str(nlev)) // ", varno = " // trim(utl_str(varno)) )
       case(2)
          call utl_abort("chm_get_obs_err_stddev: Report not found for STNID = " // trim(cstnid) // &
                       ", nlev = " // trim(utl_str(nlev)) // ", varno = " // trim(utl_str(varno)) )
       end select

       select case(chm_std%std_type(ISTNID))
       case(0)
          obs_err_stddev = sigma
       case(1)
          obs_err_stddev = max(chm_std%std3(ISTNID),chm_std%std2(ISTNID)*sigma)
       case(2)
          obs_err_stddev = sqrt(chm_std%std3(ISTNID)**2+(chm_std%std2(ISTNID)*sigma)**2)
       case(3)
          obs_err_stddev = min(chm_std%std3(ISTNID),max(chm_std%std2(ISTNID),wgt*sigma))
       case(4)
          obs_err_stddev = sqrt(chm_std%std2(ISTNID)**2+(wgt*sigma)**2)
       case default
          call utl_abort('chm_get_obs_err_stddev: std_type = ' // trim(utl_str(chm_std%std_type(ISTNID))) // &
               ' for STNID = ' // trim(CSTNID) // ' is not recognized.')
       end select
       
    END IF
    
  end function chm_get_obs_err_stddev

  subroutine chm_dealloc_obs_err_stddev
    ! 
    !:Purpose: To deallocate temporary storage space used for observation errors
    !          for the CH family.
    !
    implicit none

    integer :: istnid

    if (chm_std%n_stnid.eq.0) return
    
    if (allocated(chm_obsSub_std)) then
       do istnid=1,chm_std%n_stnid
          if (chm_std%source(istnid).ge.1) call oss_obsdata_dealloc(chm_obsSub_std(istnid))
       end do
       deallocate(chm_obsSub_std)       
    end if

    if (allocated(chm_std%stnids))   deallocate(chm_std%stnids)
    if (allocated(chm_std%n_lvl))    deallocate(chm_std%n_lvl)
    if (allocated(chm_std%std_type)) deallocate(chm_std%std_type)
    if (allocated(chm_std%ibegin))   deallocate(chm_std%ibegin)
    if (allocated(chm_std%element))  deallocate(chm_std%element)
    if (allocated(chm_std%source))   deallocate(chm_std%source)
    if (allocated(chm_std%n_lat))    deallocate(chm_std%n_lat)
    if (allocated(chm_std%std1))     deallocate(chm_std%std1)
    if (allocated(chm_std%std2))     deallocate(chm_std%std2)
    if (allocated(chm_std%std3))     deallocate(chm_std%std3)
    if (allocated(chm_std%levels))   deallocate(chm_std%levels)
    if (allocated(chm_std%lat))      deallocate(chm_std%lat)

  end subroutine chm_dealloc_obs_err_stddev

end module chem_obserrors_mod
