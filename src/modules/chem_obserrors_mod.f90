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

!--------------------------------------------------------------------------
!! MODULE chem_obserrors_mod
!!
!! *Purpose*: Provides routines regarding observation error standard deviations.
!!
!! @author Mike Sitwell and Yves Rochon (ARQI/AQRD)
!!
!! Public routines:
!!v       - "chm_read_obs_err_stddev,chm_get_obs_err_stddev,chm_dealloc_obs_err_stddev":
!!v         Routines and strucure for setting of obs error std. dev. used in
!!v        'obserrors_mod.ftn90'. 
!!
!--------------------------------------------------------------------------
module chem_obserrors_mod

  use utilities_mod
  use obsSubSpaceData_mod
  use burpFiles_mod    
  use MathPhysConstants_mod
  
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
     ! Structure containing information retrieved from file obsinfo_chm for holding 
     ! observation std dev information
     !
     !  Variable               Description
     !  --------               -----------
     !  n_stnid                Number of sub-families (identified via STNIDs)
     !  stnids                 Sub-families (STNIDs; * are wild cards)
     !  bfr                    BUFR/BURP element in data block 
     !  brp                    0: Set entirely from the ascii file being read. No 
     !                            initial values read from BURP files
     !                         1: Initial values in obs BURP files 
     !                            (may be adjusted after input)
     !                         2: Initial values in obs BURP files for variable number
     !                            of vertical levels (for error std deviations only)
     !  std_type               Index of setup approach (used in combination with
     !                            nstd_brp_chm)
     !                         For nstd_brp_chm value 0 or 1, 
     !                         0: std1 or BURP file vales (sigma)
     !                         1: max(std3,std2*ZVAL)  if brp=0
     !                            max(std3,std2*sigma) otherwise
     !                         2: sqrt(std3**2+(std2*ZVAL)**2))  if brp=0
     !                            sqrt(std3**2+(std2*sigma)**2)) otherwise
     !                         3: min(std3,max(std2,std1_chm*ZVAL)) if brp=0
     !                            min(std3,max(std2,std1_chm*sigma))  otherwise
     !                         4: sqrt(std2**2+(std1*ZVAL)**2))  if brp=0 
     !                            sqrt(std2**2+(std1*sigma)**2)) otherwise
     !  ibegin                 Position index of start of data for given
     !                         sub-family in the arrays std1,levels,lat
     !  n_lvl                  Number of vertical levels (max number when brp=2)
     !  levels                 Vertical levels (in coordinate of sub-family data)
     !  n_lat                  Number of latitudes
     !  lat                    Latitudes (degrees; ordered in increasing size)
     !  std1                   See std_type for usage (dependent on vertical level)
     !  std2                   See std_type for usage (independent of vertical level)
     !  std3                   See std_type for usage (independent of vertical level)

     integer ::  n_stnid
     character(len=12), allocatable :: stnids(:)
     integer, allocatable :: bfr(:),std_type(:),n_lat(:)
     integer, allocatable :: brp(:),ibegin(:),n_lvl(:)
     real(8), allocatable :: std1(:),std2(:),std3(:)
     real(8), allocatable :: levels(:),lat(:)

  end type struct_chm_std

  type(struct_chm_std)  :: chm_std
 
! Array of pointers to hold std's read from BURP file.
! Note: Ideally, these should be an element in the 'struct_chm_std' derived 
! types, but currently this result in an internal compiler error.

  type(struct_oss_obsdata), allocatable :: chm_burp_std(:)

contains

!-------------------------------------------------------------------------
!
! SUBROUTINE chem_read_obs_err_stddev
!!
!! *Purpose*: Read observation errors from ascii file and BURP file if specified in ascii file.
!!
!! @author M. Sitwell, March 2015
!!
!--------------------------------------------------------------------------
  subroutine chm_read_obs_err_stddev

    implicit none

    integer :: istnid

    ! read the error std. dev. information from the ascii file
    call chm_read_obs_err_stddev_ascii

    ! set size of BURP file std array
    allocate(chm_burp_std(chm_std%n_stnid))

    ! read from BURP file
    do istnid=1,chm_std%n_stnid
       if (chm_std%brp(istnid).ge.1) then
          
          ! retrieve data from stats blocks (with bkstp=14 and block_type='DATA')
          chm_burp_std(istnid) = burp_chem_read_all('CH',chm_std%stnids(istnid), &
               chm_std%bfr(istnid), chm_std%n_lvl(istnid), 1, 14, 'DATA',       &
               match_nlev=chm_std%brp(istnid).eq.1 )

       end if
    end do

  end subroutine chm_read_obs_err_stddev

!-------------------------------------------------------------------------
!
! SUBROUTINE chm_read_obs_err_stddev_ascii
!!
!! *Purpose*:  Read and store observation error std. dev. as needed for CH family obs.
!!
!! @author Y. Rochon, Nov 2014
!!
!--------------------------------------------------------------------------
  subroutine chm_read_obs_err_stddev_ascii

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
  INQUIRE(FILE='obsinfo_chm',EXIST=LnewExists)
  IF (.not.LnewExists) then
    WRITE(*,*) '---------------------------------------------------------------'
    WRITE(*,*) 'WARNING! chm_read_obs_err_stddev: obsinfo_chm not available.   '
    WRITE(*,*) 'WARNING! Default CH family stddev to be applied if needed.     '
    WRITE(*,*) '---------------------------------------------------------------'
    return
  ENDIF
!
! Read observation error std dev. from file obsinfo_chm for constituent data
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  NULSTAT=0
  IERR=FNOM(NULSTAT,'obsinfo_chm','SEQ',0)
  IF ( IERR .EQ. 0 ) THEN
    open(unit=nulstat, file='obsinfo_chm', status='OLD')
  ELSE
    CALL utl_abort('chm_read_obs_err_stddev_ascii: COULD NOT OPEN FILE obsinfo_chm!')
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
  allocate(chm_std%brp(chm_std%n_stnid),chm_std%ibegin(chm_std%n_stnid))
  allocate(chm_std%bfr(chm_std%n_stnid),chm_std%n_lvl(chm_std%n_stnid))
  allocate(chm_std%std1(isize),chm_std%std2(chm_std%n_stnid),chm_std%std3(chm_std%n_stnid))
  allocate(chm_std%levels(isize),chm_std%lat(isize))
 
  chm_std%bfr(:)=0
  chm_std%brp(:)=0
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
!        (2) Flag indication if EOR provided from this ascii file or
!            to be read from the BURP file,
!        (3) Index specifying OER setup method,
!        (4) Number of vertical levels
!        (5) Number of latitudes
!
!   Important: Combination of STNID, BUFR element and number of vertical levels
!              to determine association to the observations.
!
    read(nulstat,*,iostat=ios,err=10,end=10) chm_std%bfr(jelm),chm_std%brp(jelm),  &
       chm_std%std_type(jelm), chm_std%n_lvl(jelm), chm_std%n_lat(jelm),  &
       chm_std%std2(jelm), chm_std%std3(jelm)

    if (chm_std%n_lvl(jelm).lt.1) chm_std%n_lvl(jelm)=1
    if (chm_std%n_lat(jelm).lt.1) chm_std%n_lat(jelm)=1
    
    if (icount+chm_std%n_lvl(jelm)*chm_std%n_lat(jelm).gt.isize) then
       write(*,'(10X,"Max array size exceeded: ",I6)') isize
       CALL utl_abort('chm_read_obs_err_stddev_ascii: PROBLEM READING OBSERR STD DEV.')    
    end if

    ! disregard line of dashes
    read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne
    
    ! disregard data section if not needed
    if (chm_std%std_type(jelm).eq.1.or.chm_std%std_type(jelm).eq.2.or.(chm_std%brp(jelm).ge.1.and.chm_std%std_type(jelm).eq.0)) &
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
       CALL utl_abort('chm_read_obs_err_stddev_ascii: PROBLEM READING OBSERR STD DEV.')    
    end if
 
 11 CLOSE(UNIT=NULSTAT)
    IERR=FCLOS(NULSTAT)    

  end subroutine chm_read_obs_err_stddev_ascii

!-------------------------------------------------------------------------
!
! SUBROUTINE chm_obs_err_stddev_index(CSTNID,NLEV,VARNO,ZLAT,ISTNID,JINT)
!! 
!! *Purpose*: Returns the station ID and latitude indices corresponding to a measurement.
!!
!! @author M. Sitwell Feb 2015 
!!
!--------------------------------------------------------------------------
  subroutine chm_obs_err_stddev_index(CSTNID,NLEV,VARNO,ZLAT,ISTNID,JINT)

    implicit none

    character(len=*), intent(in)  :: CSTNID
    integer, intent(in)           :: NLEV,VARNO
    real(8), intent(in)           :: ZLAT
    integer, intent(out)          :: ISTNID,JINT
    integer                       :: JN,ilen1,ilen2,ji,ibegin
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
          IF ( (NLEV.EQ.chm_std%n_lvl(JN) .OR. chm_std%brp(JN).eq.2) .AND. VARNO.EQ.chm_std%bfr(JN) ) THEN
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
             write(*,'(2X,A,2X,I12,2X,I10)') chm_std%stnids(jn),chm_std%bfr(jn),chm_std%n_lvl(jn)
          end do
       end if
       write(*,*)
       call utl_abort('chm_obs_err_stddev_index: Check section I of file obsinfo_chm.')
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

!-------------------------------------------------------------------------
!
! FUNCTION chm_get_obs_err_stddev(CSTNID,NLEV,VARNO,ZLAT,ZLON,IDATE,ITIME,ZVAL,ZLEV,ILEV,IFIRST) result(obs_err_stddev)
!! 
!! *Purpose*: Returns the observational error std for a CH family measurement
!!
!! @author Y. Rochon and M. Sitwell, Feb 2015
!!
!! Input
!!v    cstnid        station ID
!!v    nlev          number of levels
!!v    varno         BUFR number
!!v    zlat          latitude (radians)
!!v    zlon          longitude (radians)
!!v    idate         date in YYYYMMDD format
!!v    itime         time in HHMM format
!!v    zval          observation values
!!v    zlev          vertical coordinate value
!!v    ilev          observation number in the profile
!!v    ifirst        logical indicating if first call for a profile
!!
!--------------------------------------------------------------------------
  function chm_get_obs_err_stddev(CSTNID,NLEV,VARNO,ZLAT,ZLON,IDATE,ITIME,ZVAL,ZLEV,ILEV,IFIRST) result(obs_err_stddev) 

    implicit none
   
    character(len=*), intent(in) :: CSTNID
    real(8), intent(in) :: ZVAL,ZLEV,ZLAT,ZLON
    integer, intent(in) :: NLEV,VARNO,IDATE,ITIME,ILEV
    logical, intent(in) :: IFIRST
    real(8)  :: obs_err_stddev 

    real(8) :: wgt,zwb,sigma
    integer :: ibegin,JLEV,JN,stat

    integer, save :: ISTNID,JINT
    
    ! If this call is for the first level for this measurement, get
    ! the station ID and latitude indices corresponding to this measurement
    if (ifirst) call chm_obs_err_stddev_index(CSTNID,NLEV,VARNO,ZLAT,ISTNID,JINT)                  
            
    ! Get weighting of error std. dev. if required

    if (chm_std%std_type(ISTNID).gt.2 .or. &
       (chm_std%brp(ISTNID).eq.0 .and. chm_std%std_type(ISTNID).eq.0) ) then

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
                   
    IF (chm_std%brp(ISTNID).EQ.0) THEN
               
       ! Set error standard deviations from scratch using content of
       ! previously read content of the "obsinfo_chm" file.
                
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
          
       ! Adjust error standard deviations read from BURP file if requested.
       
       sigma = oss_obsdata_get_element(chm_burp_std(istnid), oss_obsdata_get_header_code(zlon,zlat,idate,itime,cstnid), ilev, stat=stat)

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

!-------------------------------------------------------------------------
!
!! SUBROUTINE chm_dealloc_obs_err_stddev
!! 
!! *Purpose*: Deallocate temporary storage space used for observation errors for the CH family.
!!
!! @author Y. Rochon, Nov 2014
!!
!--------------------------------------------------------------------------
  subroutine chm_dealloc_obs_err_stddev
  
    implicit none

    integer :: istnid

    if (chm_std%n_stnid.eq.0) return
    
    if (allocated(chm_burp_std)) then
       do istnid=1,chm_std%n_stnid
          if (chm_std%brp(istnid).ge.1) call oss_obsdata_dealloc(chm_burp_std(istnid))
       end do
       deallocate(chm_burp_std)       
    end if

    if (allocated(chm_std%stnids))   deallocate(chm_std%stnids)
    if (allocated(chm_std%n_lvl))    deallocate(chm_std%n_lvl)
    if (allocated(chm_std%std_type)) deallocate(chm_std%std_type)
    if (allocated(chm_std%ibegin))   deallocate(chm_std%ibegin)
    if (allocated(chm_std%bfr))      deallocate(chm_std%bfr)
    if (allocated(chm_std%brp))      deallocate(chm_std%brp)
    if (allocated(chm_std%n_lat))    deallocate(chm_std%n_lat)
    if (allocated(chm_std%std1))     deallocate(chm_std%std1)
    if (allocated(chm_std%std2))     deallocate(chm_std%std2)
    if (allocated(chm_std%std3))     deallocate(chm_std%std3)
    if (allocated(chm_std%levels))   deallocate(chm_std%levels)
    if (allocated(chm_std%lat))      deallocate(chm_std%lat)

  end subroutine chm_dealloc_obs_err_stddev

end module chem_obserrors_mod