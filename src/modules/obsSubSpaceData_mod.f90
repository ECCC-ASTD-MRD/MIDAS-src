
module obsSubSpaceData_mod
  ! MODULE obsSubSpaceData_mod (prefix='oss' category='6. High-level data objects')
  !
  !:Purpose: Repository of obs space structures, arrays, and routines specific
  !          to obs data pertinent to subspaces of the overall ObsSpaceData.
  !          Most tools are independent of ObsSpaceData and can be used by
  !          themselves for the users' application(s).
  !


  ! Public routines:
  !
  !       - "oss_obsdata_get_*" to get arrays or elements from input 
  !         'struc_oss_obsdata' structure.
  !         Requires companion availabity of "oss_obsdata_alloc" and "oss_obsdata_dealloc"           
  ! 
  !       - "oss_obsdata_add_data1d" to add data value(s) to obsdata%data1d 
  !         with associated identifier code.
  !
  !       - "oss_obsdata_MPIallgather" gathers previously saved obsdata 
  !         from all processors.
  !
  !       - "oss_obsdata_code_len" to pass on oss_code_len value
  !
  !       - "oss_comboIdList" to provide list of fixed or accumulated stnid, (stnid,varno) 
  !          or (stnid,varno,multi/uni-level) combinations to be used in searches.
  !
  !       - "oss_get_comboIdList" uses the subroutine oss_comboIdlist to compile a unique list of stnid,  
  !          (stnid,varno) or (stnid,varno,multi/uni-level) combinations to be used in searches.
  !

  use codePrecision_mod
  use utilities_mod    
  use MathPhysConstants_mod
  use midasMpi_mod
  use bufr_mod
  use obsSpaceData_mod  ! for use in oss_get_comboIdList 
   
  implicit none
  private

! public procedures
! -----------------

  public :: oss_obsdata_get_data1d,oss_obsdata_add_data1d,oss_obsdata_alloc,oss_obsdata_dealloc
  public :: oss_obsdata_get_element,oss_obsdata_get_array1d,oss_obsdata_get_array2d
  public :: oss_obsdata_get_header_code,oss_obsdata_MPIallgather
  public :: oss_obsdata_code_len, oss_comboIdList, oss_get_comboIdList
  
! public types
! ------------

  public :: struct_oss_obsdata
  
! module constants
! -----------------

  integer, parameter :: oss_code_len=40            ! Max string size for code in struct_oss_obsdata.
                                                   ! Minimum required size:
                                                   ! 22 (lat/long and time coord) + 9 (stnid) = 31
  integer, parameter :: oss_code_sublen=22         ! Length of lat/long and time coord
  integer, parameter :: oss_code_latlen=5          ! Length of lat segment 

! interface for generating obsdata BURP header codes from (lat,long,date,hhmm,stnid)
  interface oss_obsdata_get_header_code
     module procedure obsdata_get_header_code_i
     module procedure obsdata_get_header_code_r
  end interface oss_obsdata_get_header_code

! module structures
! -----------------

  type :: struct_oss_obsdata

     !  Structure storing information associated to observations such
     !  as BURP file reports (irep) either retrieved from BURP files 
     !  themselves or from other sources.  
     !     
     !  Variable               Description
     !  --------               -----------
     !  ndim                   number of dimensions of the data arrays
     !                         (i.e. ndim=1 for std and ndim=2 for averaging kernels)
     !  nrep                   number of reports or observations
     !  dim1                   first array dimension for each  report/observation
     !  dim2                   second array dimension for each  report/observation (only relevant for ndim=2)
     !  data1d                 obs data of dimension (dim1,nrep)
     !  data2d                 obs data of dimension (dim1,dim2,nrep)
     !  irep                   current report position
     !  code                   unique character string for identifying the report/observation
     !
     !  Follow-up to be made: modify dim1 and dim2 as pointer arrays dependent on irep
     !
     
     real(8), pointer :: data1d(:,:),data2d(:,:,:)
     character(len=oss_code_len), pointer :: code(:)
     integer :: ndim,nrep,dim1,dim2,irep
  
  end type struct_oss_obsdata

contains


  subroutine oss_obsdata_alloc(obsdata,nrep,dim1,dim2_opt)
    !
    ! :Purpose: Allocates memory for structure struct_oss_obsdata to hold obs data file information.
    !           If dim2 is specified, then the data array associated with each observation/report
    !           will be 2D array. will be a 1D array if dim2 is not specified.
    ! 
    ! :Arguments:
    !      :obsdata: data structure to allocate
    !      :nrep: max number of associated observations/reports for data 
    !      :dim1: first dimension length of the array associated to each observation/report
    !      :dim2: second dimension length of the array associated to each observation/report (optional)
    !
    implicit none

    ! Arguments:
    type(struct_oss_obsdata), intent(inout) :: obsdata
    integer                 , intent(in)    :: nrep
    integer                 , intent(in)    :: dim1 
    integer       , optional, intent(in)    :: dim2_opt
    
    obsdata%nrep = nrep
    obsdata%dim1 = dim1

    if (present(dim2_opt)) then
       obsdata%dim2 = dim2_opt
       obsdata%ndim = 2
       allocate(obsdata%data2d(dim1,dim2_opt,nrep))
       obsdata%data2d(:,:,:) = 0.0D0
       obsdata%data1d => null()
    else
       obsdata%dim2 = 0
       obsdata%ndim = 1
       allocate(obsdata%data1d(dim1,nrep))
       obsdata%data1d(:,:) = 0.0D0
       obsdata%data2d => null()
    end if

    ! code is a character string assigned to each observation/report to uniquely identify it
    allocate(obsdata%code(nrep))
    
    ! obsdata%irep is a counter used to keep tract of position in the
    ! data arrays when extracting values via oss_obsdata_get_array*.
    ! The value is initialized to one, pointing to the very first element.
    obsdata%irep = 1

  end subroutine oss_obsdata_alloc


  subroutine oss_obsdata_dealloc(obsdata)
    !
    ! :Purpose: Deallocates memory for structure struct_oss_obsdata
    !
    implicit none

    ! Arguments:
    type(struct_oss_obsdata), intent(inout) :: obsdata
    
    if (associated(obsdata%data1d)) deallocate(obsdata%data1d)
    if (associated(obsdata%data2d)) deallocate(obsdata%data2d)
    if (associated(obsdata%code))   deallocate(obsdata%code)

  end subroutine oss_obsdata_dealloc


  function oss_obsdata_get_element( obsdata, code, idim1, stat_opt ) result(element)
    ! 
    ! :Purpose: Returns element of array from obsdata 1D data array. The returned element is the
    !           one with the specified identifying code.
    !
    ! :Arguments: 
    !           :obsdata: struct_oss_obsdata instance
    !           :code:    unique identifying code
    !           :idim1:   position of element in the dim1 axis
    !           :element: retrieved element from obsdata%data1d
    !           :stat:    status code
    !
    implicit none

    ! Arguments:
    type(struct_oss_obsdata), intent(inout) :: obsdata
    character(len=*)        , intent(in)    :: code
    integer                 , intent(in)    :: idim1
    integer       , optional, intent(out)   :: stat_opt
    ! Result:
    real(8) :: element
    
    ! find obsdata%irep for current observation
    call obsdata_set_index(obsdata,code,stat_opt=stat_opt)
    if (present(stat_opt)) then
       if (stat_opt.ne.0) then
          element = 0.
          return
       end if
    end if

    ! get element from data array at current position
    element = obsdata%data1d(idim1,obsdata%irep)
    
    ! increment position in data array
    if (obsdata%irep.eq.obsdata%nrep) then
       obsdata%irep = 1
    else
       obsdata%irep=obsdata%irep+1
    end if

  end function oss_obsdata_get_element


  function oss_obsdata_get_array1d(obsdata,code,stat_opt) result(array)
    ! 
    ! :Purpose: Returns 1D data array from obsdata. The returned array is the one with the specified
    !           identifying code.
    ! :Arguments:
    !           :obsdata: struct_oss_obsdata instance
    !           :code:    unique identifying code
    !           :array:   retrieved array from obsdata%data1d of dimension obsdata%dim1
    !           :stat:    search success (0 - found; 1 = no data; 2 = not found)
    !
    implicit none

    ! Arguments:
    type(struct_oss_obsdata), intent(inout) :: obsdata
    character(len=*)        , intent(in)    :: code
    integer       , optional, intent(out)   :: stat_opt
    ! Result:
    real(8) :: array(obsdata%dim1)
    
    ! find obsdata%irep for current observation
    call obsdata_set_index(obsdata,code,stat_opt=stat_opt)
    if ( present( stat_opt ) ) then
      if ( stat_opt /= 0 ) then
        array(:) = 0.
        return
      end if
    end if
    
    ! get element from data array at current position
    array = obsdata%data1d(:,obsdata%irep)
    
    ! increment position in data array
    if (obsdata%irep.eq.obsdata%nrep) then
       obsdata%irep = 1
    else
       obsdata%irep=obsdata%irep+1
    end if

  end function oss_obsdata_get_array1d


  function oss_obsdata_get_data1d( obsdata, lon, lat, date, time, stnid, stat_opt ) result(array)
    !
    ! :Purpose: Extract 1D data array from structure according to input (lat,long,date,time,stnid)
    !
    ! :Arguments:
    !           :obsdata: Structure from which data elements are to be found  
    !           :lon:     longitude real (radians)
    !           :lat:     latitude real (radians)
    !           :date:    date in YYYYMMDD
    !           :time:    time in HHMM
    !           :stnid:   station ID
    !           :array:   Identified 1D array
    !           :stat:    search success (0 - found; 1 = no data; 2 = not found)
    !    
    implicit none

    ! Arguments:
    type(struct_oss_obsdata), intent(inout) :: obsdata
    real(8)                 , intent(in)    :: lon
    real(8)                 , intent(in)    :: lat
    integer                 , intent(in)    :: date
    integer                 , intent(in)    :: time
    character(len=*)        , intent(in)    :: stnid
    integer       , optional, intent(out)   :: stat_opt
    ! Result:
    real(8) :: array(obsdata%dim1)

    ! Locals:
    character(len=oss_code_len) :: code
    
    ! Set desired identifier code
    code=oss_obsdata_get_header_code(lon,lat,date,time,stnid) 
  
    ! Get array corresponding to code.
    array=oss_obsdata_get_array1d(obsdata,code,stat_opt=stat_opt)
    
  end function oss_obsdata_get_data1d


  function oss_obsdata_get_array2d( obsdata, code, stat_opt ) result(array)
    ! 
    ! :Purpose: Returns 2D data array from obsdata. The returned array is the one with the specified
    !           identifying code.
    !
    ! :Arguments: 
    !           :obsdata: struct_oss_obsdata instance
    !           :code:    unique identifying code
    !           :array:   retrieved array from obsdata%data2d of dimension (obsdata%dim1,obsdata%dim2)
    !           :stat:    search success (0 - found; 1 = no data; 2 = not found)
    !
    implicit none

    ! Arguments:
    type(struct_oss_obsdata), intent(inout) :: obsdata
    character(len=*)        , intent(in)    :: code
    integer       , optional, intent(out)   :: stat_opt
    ! Result:
    real(8) :: array(obsdata%dim1,obsdata%dim2)
    
    ! find obsdata%irep for current observation
    call obsdata_set_index(obsdata,code,stat_opt=stat_opt)
    if (present(stat_opt)) then
       if (stat_opt.ne.0) then
          array(:,:) = 0.
          return
       end if
    end if
    
    ! get elements from data array at current position
    array = obsdata%data2d(:,:,obsdata%irep)
    
    ! increment position in data array
    if (obsdata%irep.eq.obsdata%nrep) then
       obsdata%irep = 1
    else
       obsdata%irep=obsdata%irep+1
    end if

  end function oss_obsdata_get_array2d


  subroutine obsdata_set_index( obsdata, code, stat_opt )
    ! 
    ! :Purpose: Sets the position variable (irep) in struct_oss_obsdata to reference the record
    !           that matches the input identifying code.
    !
    ! :Arguments: 
    !           :obsdata:      struct_oss_obsdata instance
    !           :index:        obs index
    !           :code:         code for comparison to those in obsdata
    !           :obsdata%irep: current index position for the observation/report
    !           :stat_opt:     status of call (optional)
    !                    :0: no errors
    !                    :1: no reports available
    !                    :2: report not found
    !
    ! :Comments: If the optional argument stat_opt is provided and an error occurs, the error code will
    !            be returned and the abort will not be called to allow for error handling.
    !
    implicit none

    ! Arguments:
    type(struct_oss_obsdata), intent(inout) :: obsdata
    character(len=*)        , intent(in)    :: code
    integer       , optional, intent(out)   :: stat_opt

    ! Locals:
    integer :: i
    integer :: ref_lat
    
    if ( obsdata%nrep <= 0 ) then
      if ( present( stat_opt ) ) then
        stat_opt = 1
        return
      else
        call utl_abort("obsdata_set_index: No reports available. Check for consistency " // &
                       "between input BURP file and input NAMBURP_FILTER_*  namelist " // &
                       "(and input auxiliary file if part of CH family).")
      end if
    end if

    i=0
    
    ! Search for matching identifier code
    do while (trim(obsdata%code(obsdata%irep)) /= trim(code))
       obsdata%irep=obsdata%irep+1
       if (obsdata%irep > obsdata%nrep) obsdata%irep=1
       if (i > obsdata%nrep) then
          if (len_trim(code) >= oss_code_sublen) then
       
             ! Assumes codes of the form "LAT--LON--YYYYMMDDHHMM*" when len(code)>=oss_code_sublen.
	  
             ! Upper loop search did not find a match. For valid data near the poles over
             ! the global analysis grid, the lat could have been moved to the nearest analysis grid
             ! latitude. The following is to account for this, assuming this is the only exception
             ! for points near the poles. 
             !
             ! This is done as a second search step so as not to slow down the normally sufficient 
             ! search performed above.

             i=0
             read(code(1:oss_code_latlen),*) ref_lat
             
             ! Search for matching identifier code
             do while (.not.obsdata_extra_code_test(trim(obsdata%code(obsdata%irep)),code,ref_lat))
                obsdata%irep=obsdata%irep+1
                if (obsdata%irep > obsdata%nrep) obsdata%irep=1
                if (i > obsdata%nrep) exit
                i=i+1       
             end do
             if (i > obsdata%nrep) then
                if (present(stat_opt)) then
                   stat_opt = 2
                   return
                else
                   call utl_abort("obsdata_set_index: Obs index not found for nrep = " // trim(utl_str(obsdata%nrep)) // " and code = '" // code // "'")
                end if
             end if
          else
             if (present(stat_opt)) then
                stat_opt = 2
                return
             else
                call utl_abort("obsdata_set_index: Obs index not found for nrep = " // trim(utl_str(obsdata%nrep)) // " and code = '" // code // "'")
             end if
          end if
          exit
       end if 
       i=i+1       
    end do
         
    if (present(stat_opt)) stat_opt = 0

  end subroutine obsdata_set_index
             

  function obsdata_extra_code_test( test_code, ref_code, ref_lat ) result(found)
    ! 
    ! :Purpose: Test matching of code values accounting for rare differences
    !           in stored lat (and lon) value(s) when codes are stored as strings in the form  
    !           LAT--LON--YYYYMMDDHHMM* (ie. with >= oss_code_sublen characters).
    !
    ! :Caveat: The current version assumes the only source of difference would stem from
    !          a shift to the nearest latitude of the analysis grid from points near the pole.
    !          (this source of difference identified by M. Sitwell)
    !          Also currently assumes that most poleward analysis grid latitudes are within 1 degree
    !          away from a pole.
    !
    ! :Arguments:
    !           :test_code: code for comparison to ref_code  
    !           :ref_lat:   latitude  (x100) part of reference code   
    !           :ref_code:  reference code
    !           :found:     logical indicating if a match has been found.  
    ! 
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: test_code
    character(len=*), intent(in) :: ref_code
    integer         , intent(in) :: ref_lat
    ! Result:
    logical :: found

    ! Locals:
    integer, parameter :: lat_lim1=-8900    ! Lat*100
    integer, parameter :: lat_lim2=8900
    integer :: lat
          
    if (test_code(6:len_trim(test_code)) /= ref_code(oss_code_latlen+1:len_trim(ref_code))) then
      found=.false.
      return
    else 
      read(test_code(1:oss_code_latlen),*) lat
      if ((lat < lat_lim1 .and. ref_lat < lat_lim1 .and. lat < ref_lat ).or. &
          (lat > lat_lim2 .and. ref_lat > lat_lim2 .and. lat > ref_lat ) ) then
        found=.true.	     
        write(*,*) 'obsdata_extra_code_test: Accounted for lat. mismatch in codes near poles: ', &
              lat,ref_lat
      else
        found=.false.
      end if
    end if

  end function obsdata_extra_code_test
    

  function obsdata_get_header_code_i( ilon, ilat, date, time, stnid ) result(code)
    ! 
    ! :Purpose: Generates a string code to identify an obervation by the header information in
    !           a BURP report. The BURP header information is saved as a string in the form  
    !           LAT--LON--YYYYMMDDHHMMSTNID----. Intention of this function is to be used for
    !           setting the unique identifier 'code' in struct_oss_obsdata. Can be called under
    !           the interface oss_obsdata_get_header_code.
    !
    ! :Arguments:
    !           :ilon:  longitude integer
    !           :ilat:  latitude integer
    !           :date:  date in YYYYMMDD
    !           :time:  time in HHMM
    !           :stnid: station ID
    !           :code:  unique code
    ! 
    implicit none

    ! Arguments:
    integer,          intent(in) :: ilon
    integer,          intent(in) :: ilat
    integer,          intent(in) :: date
    integer,          intent(in) :: time
    character(len=*), intent(in) :: stnid
    ! Result:
    character(len=oss_code_len) :: code

    write(code(1:5),'(I5.5)') ilat

    if (ilon > 0) then
      write(code(6:10),'(I5.5)') ilon
    else
      write(code(6:10),'(I5.5)') 36000 + ilon
    end if
    
    write(code(11:18),'(I8.8)') date
    write(code(19:22),'(I4.4)') time

    write(code(23:oss_code_len),'(A)') stnid(1:min(len_trim(stnid),oss_code_len-oss_code_sublen))

  end function obsdata_get_header_code_i


  function obsdata_get_header_code_r( lon, lat, date, time, stnid ) result(code)
    ! 
    ! :Purpose: Generates a string code to identify an obervation by the header information in
    !           a BURP report. The BURP header information is saved as a string in the form  
    !           LAT--LON--YYYYMMDDHHMMSTNID----. Intention of this function is to be used for
    !           setting the unique identifier 'code' in struct_oss_obsdata. Can be called under
    !           the interface oss_obsdata_get_header_code.
    !
    ! :Arguments:
    !           :lon:   longitude real (radians)
    !           :lat:   latitude real (radians)
    !           :date:  date in YYYYMMDD
    !           :time:  time in HHMM
    !           :stnid: station ID
    !           :code:  unique code
    ! 
    implicit none

    ! Arguments:
    real(8)         , intent(in) :: lon
    real(8)         , intent(in) :: lat
    integer         , intent(in) :: date
    integer         , intent(in) :: time
    character(len=*), intent(in) :: stnid
    ! Result:
    character(len=oss_code_len)  :: code

    ! Locals:
    integer :: ilon, ilat
    
    ilon = nint(100*(lon/MPC_RADIANS_PER_DEGREE_R8))
    ilat = nint(100*(90. + lat/MPC_RADIANS_PER_DEGREE_R8))

    code = obsdata_get_header_code_i(ilon,ilat,date,time,stnid)

  end function obsdata_get_header_code_r

!----------------------------------------------------------------------------------------

  subroutine oss_obsdata_add_data1d(obsdata,val,code,maxsize,dim1_opt)
    !
    ! :Purpose: Add data value(s) to obsdata%data1d with associated identifier code.
    !
    ! :Arguments:
    !           obsdata       struct_oss_obsdata instance
    !           val           data array to store in obsdata%data1d
    !           code          identifying code based on (lat,long,date,hhmm) if not also stnid
    !           maxsize       max allowed size for obsdata
    !           dim1          value() dimension (optional)
    !           obsdata       Updated obsdata 
    ! 
    ! :Comments:
    !
    !          - Retrieval of values from obsdata%data1d to be done via oss_obsdata_get_element (or oss_obsdata_get_array1d).
    !          - If obsdata allocation is required for all processors (such as for use later with obsdata_MPIGather), 
    !            allocation and/or initialization of arrays needs to be done at a corresponding appropriate location 
    !            outside the obs operations such as in oss_setup to ensure allocation is done for all processors, 
    !            including those without associated data. This is to ensure that rpn_comm_allgather will work 
    !            in routine obsdata_MPIGather.
    !
    implicit none
    
    ! Arguments:
    type(struct_oss_obsdata), intent(inout) :: obsdata
    real(8)                 , intent(in)    :: val(:)
    integer                 , intent(in)    :: maxsize
    character(len=*)        , intent(in)    :: code
    integer       , optional, intent(in)    :: dim1_opt

    if (.not.associated(obsdata%data1d)) then
      if (present(dim1_opt)) then 
         call oss_obsdata_alloc(obsdata,maxsize,dim1=dim1_opt)
      else
         call oss_obsdata_alloc(obsdata,maxsize,dim1=1)
      end if
      obsdata%nrep=0
    end if

    if (obsdata%dim1 > size(val)) &
         call utl_abort('obsdata_add_data1d: Insufficient data values provided. ' // trim(utl_str(size(val))) )
     
    ! nrep counts the number data values/profiles in the data arrays
    obsdata%nrep = obsdata%nrep+1 

    if (obsdata%nrep > maxsize) &
         call utl_abort('obsdata_add_data1d: Reach max size of array ' // trim(utl_str(maxsize)) )
  
    ! Save unique code
    obsdata%code(obsdata%nrep)=trim(code)
    
    ! Save value(s)
    obsdata%data1d(1:obsdata%dim1,obsdata%nrep) = val(1:obsdata%dim1)
    
  end subroutine oss_obsdata_add_data1d

    
  subroutine oss_obsdata_MPIallgather(obsdata) 
    !
    ! :Purpose: Gathers previously saved obsdata from all processors.
    !
    ! :Arguments:
    !           :obsdata: Local struct_oss_obsdata to become global
    ! 
    ! :Comments:
    !
    !           - Assumes obsdata%dim1 (and obsdata%dim2) the same over all processors.
    !
    implicit none

    ! Arguments:
    type(struct_oss_obsdata), intent(inout) :: obsdata

    ! Locals:
    integer, allocatable :: nrep(:)
    character(len=oss_code_len), allocatable :: code_local(:),code_global(:,:)
    real(8), allocatable :: data1d_local(:,:),data1d_global(:,:,:),data2d_local(:,:,:),data2d_global(:,:,:,:)
    integer :: i,ierr,nproc,nrep_total,nrep_max,irep,array_size

    write(*,*) 'Begin oss_obsdata_MPIallgather'
      
    ! Identify number of processors.
    call rpn_comm_size("GRID",nproc,ierr)
      
    ! Get number of reports on each processor

    allocate(nrep(nproc))
    nrep(:)=0

    call rpn_comm_allgather(obsdata%nrep,1,"MPI_INTEGER",nrep,1,"MPI_INTEGER","GRID",ierr)

    nrep_total=sum(nrep)

    if (nrep_total == 0) then
       deallocate(nrep)
       write(*,*) 'Exit oss_obsdata_MPIallgather: no reports'
       return
    end if
        
    nrep_max = maxval(nrep)

    ! Get values from all processors into global arrays

    allocate(code_local(nrep_max))
    allocate(code_global(nrep_max,nproc))

    code_local(:)=''
    if (obsdata%nrep > 0) code_local(1:obsdata%nrep) = obsdata%code(1:obsdata%nrep)

    call mmpi_allgather_string(code_local,code_global,nrep_max,oss_code_len,nproc,"GRID",ierr)

    if (obsdata%ndim == 1) then

       allocate(data1d_local(obsdata%dim1,nrep_max))
       allocate(data1d_global(obsdata%dim1,nrep_max,nproc))

       data1d_local(:,:)=0.0D0
       if (obsdata%nrep > 0) data1d_local(:,1:obsdata%nrep) = obsdata%data1d(:,1:obsdata%nrep)
       
       array_size = nrep_max*obsdata%dim1

       call rpn_comm_allgather(data1d_local,array_size,pre_obsMpiReal,data1d_global,array_size,pre_obsMpiReal,"GRID",ierr)
    else 
       
       allocate(data2d_local(obsdata%dim1,obsdata%dim2,nrep_max))
       allocate(data2d_global(obsdata%dim1,obsdata%dim2,nrep_max,nproc))
       
       data2d_local(:,:,:)=0.0D0
       if (obsdata%nrep > 0) data2d_local(:,:,1:obsdata%nrep) = obsdata%data2d(:,:,1:obsdata%nrep)

       array_size = nrep_max*obsdata%dim1*obsdata%dim2

       call rpn_comm_allgather(data2d_local,array_size,pre_obsMpiReal,data2d_global,array_size,pre_obsMpiReal,"GRID",ierr)

    end if
  
    deallocate(code_local)
    if (obsdata%ndim == 1) then 
        deallocate(data1d_local)
    else
        deallocate(data2d_local)
    end if 

    ! Concatenate values from all processors
         
    irep = 0
    call oss_obsdata_dealloc(obsdata)    

    if (obsdata%ndim == 1) then      
       
       call oss_obsdata_alloc(obsdata,nrep_total,obsdata%dim1)

       do i=1,nproc
          if (nrep(i) > 0) then
             obsdata%code(irep+1:irep+nrep(i)) = code_global(1:nrep(i),i)
             obsdata%data1d(1:obsdata%dim1,irep+1:irep+nrep(i)) = data1d_global(1:obsdata%dim1,1:nrep(i),i)
             irep = irep+nrep(i)
          end if
       end do

    else

       call oss_obsdata_alloc(obsdata,nrep_total,obsdata%dim1,obsdata%dim2)

       do i=1,nproc
          if (nrep(i) > 0) then
             obsdata%code(irep+1:irep+nrep(i)) = code_global(1:nrep(i),i)
             obsdata%data2d(1:obsdata%dim1,1:obsdata%dim2,irep+1:irep+nrep(i)) = &
                         data2d_global(1:obsdata%dim1,1:obsdata%dim2,1:nrep(i),i)  
             irep = irep+nrep(i)
          end if
       end do
 
    end if  
     
    deallocate(nrep,code_global)
    if (obsdata%ndim == 1) then 
        deallocate(data1d_global)
    else
        deallocate(data2d_global)
    end if 
    
    write(*,*) 'Exit oss_obsdata_MPIallgather'
 
  end subroutine oss_obsdata_MPIallgather


  subroutine oss_get_comboIdlist( obsSpaceData, stnid_list, varno_list, unilev_list, num_elements, nset )
    ! 
    ! :Purpose: Uses the subroutine oss_comboIdlist to compile a unique list of stnid,  
    !           (stnid,varno) or (stnid,varno,multi/uni-level) combinations to be used in searches.
    !
    ! :Arguments:
    !           :obsSpaceData: Observation space data
    !           :stnid_list:   List of unique stnids
    !           :varno_list:   List of unique varno
    !           :unilev_list:  List of unique uni/multi-level identifications
    !           :num_elements: Number of unique elements in *_list arrrays
    !           :nset:         Integer indicating grouping, with values indicating
    !
    !                    - 1: group by stnid
    !                    - 2: group by (stnid,bufr)
    !                    - 3: group by (stnid,bufr,multi/uni-level)
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData

    ! Locals:
    integer, parameter :: nmax=100
    integer, intent(out) :: num_elements
    integer, intent(out) :: nset
    integer, intent(out) :: varno_list(nmax)
    character(len=9), intent(out) :: stnid_list(nmax)
    logical, intent(out) :: unilev_list(nmax)
    integer :: headerIndex,bodyIndex,vco,nlev_obs,varno
    logical :: all_combos

    call oss_comboIdlist(all_combos_opt=all_combos)
    
    if (all_combos) then
    
       ! Loop over obs to find all (stnid,varno) pairs to form a common sequence of search pairs
       ! over all processors. The prescribed starting stnid selections can have wild card 
       ! characters (via *, see routine oss_comboIdlist).

       call obs_set_current_header_list(obsSpaceData,'CH')
       HEADER: do
          headerIndex = obs_getHeaderIndex(obsSpaceData)
          if (headerIndex < 0) exit HEADER
         
          ! Body info that we only need for first point in the profile
          bodyIndex = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)     
          vco = obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex)

          if (vco.ne.1 .and. vco.ne.2 .and. vco.ne.4 .and. vco.ne.5) then
             ! Vertical coordinate not handled
             write(*,*) 'oss_get_comboIdlist: Currently unaccounted VCO = ',vco
             cycle HEADER
          end if
          
          ! Identify varno and nlev_obs (exclude BUFR_SCALE_EXPONENT elements)
          nlev_obs = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex)
          call obs_set_current_body_list(obsSpaceData,headerIndex)
          do
             bodyIndex = obs_getBodyIndex(obsSpaceData)
             if (bodyIndex < 0) exit
             if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex).eq.BUFR_SCALE_EXPONENT) then
                nlev_obs = nlev_obs-1
             else
                varno=obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
             end if
          end do

          ! Adds to running list of unique pairs if unique
          call oss_comboIdlist(stnid_add_opt=obs_elem_c(obsSpaceData,'STID',headerIndex), varno_add_opt=varno, unilev_add_opt=(nlev_obs.eq.1.and.vco.ge.4))
                       
       end do HEADER
       
       ! Get a common sequence of search pairs over all processors. 
       
       call oss_comboIdlist(gather_mpi_opt=.true.)
       
    end if
    
    ! Get list of unique pairs
    call oss_comboIdlist(stnid_list_opt=stnid_list, varno_list_opt=varno_list, unilev_list_opt=unilev_list, num_elements_opt=num_elements, nset_opt=nset)

  end subroutine oss_get_comboIdlist

!-----------------------------------------------------------------------------------
  
  subroutine oss_comboIdList( stnid_add_opt, varno_add_opt, unilev_add_opt, stnid_list_opt, &
                              varno_list_opt, unilev_list_opt, &
                              num_elements_opt, initialize_opt, nset_opt, gather_mpi_opt, all_combos_opt )
    ! 
    ! :Purpose: Provide list of fixed or accumulated stnid, (stnid,varno) or 
    !           (stnid,varno,multi/uni-level) combinations to be used in searches.
    !
    !           Can be used for both single processor and  MPI mode, where 'gather_mpi' must be set
    !           to .true. at some point for use with MPI.
    !            Called from osd_chem_diagnmostics in file obspacediag_mod.ftn90.
    !  
    ! :Arguments:
    !           :stnid_add_opt:    stnid to add to stnid_list if part of unique set
    !           :varno_add_opt:    varno to add to varno_list if part of unique set
    !           :unilev_add_opt:   unilev logical to add to unilev_list if part of unique set
    !           :initialize_opt:   Initialize internal arrays and counters
    !           :gather_mpi_opt:   Gather all local MPI process and recompile unique lists  
    !           :nset_opt:         Integer indicating grouping of diagnostics. Input variable
    !                              if initialize=.true., output variable otherwise.
    !                              Values indicate
    !
    !                              - 1: group by stnid
    !                              - 2: group by (stnid,bufr)
    !                              - 3: group by (stnid,bufr,multi/uni-level)
    !           :all combos_opt:   Indicates if all combinations specified by nset are to
    !                              be use, or only those specified in the namelist NAMCHEM
    !                              Input variable if initialize=.true., output variable otherwise.
    !           :stnid_list_opt:   List of unique stnids
    !           :varno_list_opt:   List of unique varno
    !           :unilev_list_opt:  List of unique uni/multi-level identifications
    !           :num_elements_opt: Number of unique elements in *_list arrrays
    !
    implicit none
    
    integer, parameter :: nmax=100, stnid_len=9

    ! Arguments:
    logical                 , intent(in)   , optional :: initialize_opt
    logical                 , intent(in)   , optional :: gather_mpi_opt
    logical                 , intent(in)   , optional :: unilev_add_opt
    integer                 , intent(in)   , optional :: varno_add_opt
    character(len=stnid_len), intent(in)   , optional :: stnid_add_opt
    integer                 , intent(inout), optional :: nset_opt
    logical                 , intent(inout), optional :: all_combos_opt
    integer                 , intent(out)  , optional :: varno_list_opt(nmax)
    integer                 , intent(out)  , optional :: num_elements_opt
    character(len=stnid_len), intent(out)  , optional :: stnid_list_opt(nmax)
    logical                 , intent(out)  , optional :: unilev_list_opt(nmax)

    ! Locals:
    integer, save :: varno_unique(nmax)
    character(len=stnid_len), save :: stnid_unique(nmax)
    logical, save :: unilev_unique(nmax)
    integer, allocatable :: num_unique_all(:),varno_unique_all(:,:)
    character(len=stnid_len),allocatable :: stnid_unique_all(:,:)
    logical, allocatable :: unilev_unique_all(:,:)
    integer, save :: num_unique ! running count of number of unique elements
    integer, save :: iset=2
    logical, save :: lall_combos=.true.
    logical :: same,init
    integer :: i,j,nproc,iproc,ierr

    init=.false.
    if (present(initialize_opt)) init = initialize_opt

    
    ! Initialize internal arrays and counters
    if (init) then
       stnid_unique(:) = ''
       varno_unique(:) = 0
       unilev_unique(:) = .false.
       num_unique = 0
       if (present(nset_opt)) iset = nset_opt
       if (present(all_combos_opt)) lall_combos = all_combos_opt
    end if      


    ! Add new elements to internal arrays if not there already
    if (present(stnid_add_opt)) then
      
       if (iset >= 2 .and. (.not. present(varno_add_opt))) call utl_abort('oss_comboIdlist: varno_add must be present to add element for nset>=2.')
       if (iset >= 3 .and. (.not. present(unilev_add_opt))) call utl_abort('oss_comboIdlist: unilev_add must be present to add element for nset>=3.')

       same = .false.

       do i=1,num_unique
          same = utl_stnid_equal(stnid_add_opt,stnid_unique(i))
          if (iset >= 2) same = same .and. varno_add_opt == varno_unique(i)
          if (iset >= 3) same = same .and. unilev_add_opt.eqv.unilev_unique(i)
          if (same) exit
       end do

       if (.not.same) then
          num_unique=num_unique+1
          if (num_unique > nmax) call utl_abort("oss_comboIDlist: Max allowed dimension exceeded.")
          stnid_unique(num_unique) = stnid_add_opt
          if (iset >= 2) varno_unique(num_unique) = varno_add_opt
          if (iset >= 3) unilev_unique(num_unique) = unilev_add_opt
       end if
    end if        

    ! Gather unique arrays from each local mpi process and compile global unique arrays
    if (present(gather_mpi_opt)) then
       if (gather_mpi_opt) then

          call rpn_comm_size("GRID",nproc,ierr)

          allocate(num_unique_all(nproc))
          allocate(stnid_unique_all(nmax,nproc))
          if (iset >= 2) allocate(varno_unique_all(nmax,nproc))
          if (iset >= 3) allocate(unilev_unique_all(nmax,nproc))
          
          num_unique_all(:) = 0
          stnid_unique_all(:,:) = ''
          if (iset >= 2) varno_unique_all(:,:) = 0
          if (iset >= 3) unilev_unique_all(:,:) = .false.
          
          if(mmpi_doBarrier) call rpn_comm_barrier("GRID",ierr)

          call rpn_comm_allgather(num_unique,1,"MPI_INTEGER",num_unique_all,1,"MPI_INTEGER","GRID",ierr)
          call mmpi_allgather_string(stnid_unique,stnid_unique_all,nmax,stnid_len,nproc,"GRID",ierr)
          if (iset >= 2) call rpn_comm_allgather(varno_unique,nmax,"MPI_INTEGER",varno_unique_all,nmax,"MPI_INTEGER","GRID",ierr)
          if (iset >= 3) call rpn_comm_allgather(unilev_unique,nmax,"MPI_LOGICAL",unilev_unique_all,nmax,"MPI_LOGICAL","GRID",ierr)
          
          stnid_unique(:) = ''
          if (iset >= 2) varno_unique(:) = 0
          if (iset >= 3) unilev_unique(:) = .false.
          num_unique = 0
          
          ! Amalgamate unique lists
          do iproc=1,nproc
             do j=1,num_unique_all(iproc)

                same = .false.

                do i=1,num_unique
                   same = utl_stnid_equal(stnid_unique_all(j,iproc),stnid_unique(i))
                   if (iset >= 2) same = same .and. varno_unique_all(j,iproc) == varno_unique(i)
                   if (iset >= 3) same = same .and. unilev_unique_all(j,iproc).eqv.unilev_unique(i)
                   if (same) exit
                end do
              
                if (.not.same) then
                   num_unique=num_unique+1
                   if (num_unique > nmax) call utl_abort("oss_comboIDlist: Max allowed dimension exceeded.")
                   stnid_unique(num_unique) = stnid_unique_all(j,iproc)
                   if (iset >= 2) varno_unique(num_unique) = varno_unique_all(j,iproc)
                   if (iset >= 3) unilev_unique(num_unique) = unilev_unique_all(j,iproc)
                end if
             end do
          end do

          deallocate(num_unique_all,stnid_unique_all)
          if (allocated(varno_unique_all)) deallocate(varno_unique_all)
          if (allocated(unilev_unique_all)) deallocate(unilev_unique_all)

       end if
    end if


    ! Return internal arrays (and other info) if requested
    if (present(varno_list_opt)) varno_list_opt = varno_unique
    if (present(stnid_list_opt)) stnid_list_opt = stnid_unique
    if (present(unilev_list_opt)) unilev_list_opt = unilev_unique
    if (present(num_elements_opt)) num_elements_opt = num_unique
    if (.not. init) then
       if (present(nset_opt)) nset_opt = iset
       if (present(all_combos_opt)) all_combos_opt = lall_combos
    end if
    
    
  end subroutine oss_comboIdList
  

  integer function oss_obsdata_code_len()  
    !
    ! :Purpose: Pass on oss_code_len value.
    !
    implicit none
    
    oss_obsdata_code_len=oss_code_len
    
   end function oss_obsdata_code_len

end module obsSubSpaceData_mod
