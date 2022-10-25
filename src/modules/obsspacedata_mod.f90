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


! ObsSpaceData_mod:  the module, ObsSpaceData_mod, follows IndexListDepot_mod

! NOTE:  Throughout this file:
!             column_index   - is not (in general) indexed from one.  Each column
!                              index has an equivalent name, OBS_*.
!             active_index   - is indexed from one by definition (a column index)
!             row_index      - is indexed from one.  It has no equivalent name.
!             bodyIndex, etc.- necessarily a row index
!             HeaderIndex,etc.-necessarily a row index

module IndexListDepot_mod
  ! MODULE indexListDepot_mod (prefix='ild' category='7. Low-level data objects and utilities')
  !
  ! :Purpose: The raison d'etre of this module is to support ObsSpaceData_mod in
  !           facilitating the traversal of a selection of the rows in its table.
  !           The selection of rows could be from either the header table or the
  !           body table. ObsSpaceData_mod currently populates the list with one
  !           of:
  !                 all header members of a given family
  !                 all   body members of a given family
  !                 all   body members of a given header row index
  !
  ! :Usage:
  !    An ObsSpaceData_mod client must first call either
  !    obs_set_current_body_list or obs_set_current_header_list, specifying
  !    either the family of interest or the header row index of interest. 
  !    This does not return the list directly to the caller, but rather writes
  !    the list, as a struct_index_list, to the private contents of the obs
  !    oject that is returned to the caller but which cannot be examined by the
  !    caller.  Two lists can be active simultaneously:  one header list and
  !    one body list.
  !
  !    In order to access the indices that are in the list, the
  !    ObsSpaceData_mod client must call either obs_getHeaderIndex or
  !    obs_getBodyIndex, giving the ObsSpaceData object as the only argument. 
  !    On each call, one index is returned.  On calls after the last index in
  !    the list, a value of -1 is returned.
  !
  !    This is not a fully fledged module.  It is better described as a
  !    structure definition with a couple of helpful methods.  It is intended
  !    that the client, ObsSpaceData_mod, read/write directly from/to instances
  !    of these structures.



  ! STRUCT_INDEX_LIST:
  !    A struct_index_list contains the identity of the family or header that
  !    was used to create the list, the actual list of indices, and the number
  !    of the list element that was last returned to the user.
  !
  ! STRUCT_INDEX_LIST_DEPOT:
  !    Because it is typical for a client to traverse a small group of lists
  !    several times each, performance is improved by retaining recent lists,
  !    thus avoiding having to regenerate them on each request.  Recent lists
  !    are stored in a struct_index_list_depot.  ObsSpaceData_mod contains one
  !    struct_index_list_depot for header lists and another for body lists.
  !    The struct_index_list_depot structure contains current_list, a pointer
  !    to the list in the depot that was last requested by the ObsSpaceData_mod
  !    client.
  !
  ! OMP:
  !    ObsSpaceData_mod has been designed so that it may be called from within
  !    an OMP section of code.  If there are n OMP threads, it is possible that
  !    there be as many as n lists in use simultaneously.  The parameter, 
  !    NUMBER_OF_LISTS, has been set to n to accommodate that many threads.  In
  !    this case, because the current_list of the depot is not OMP-private, it
  !    cannot be asked to remember the current list for each of the OMP
  !    threads.  Therefore, obs_set_current_header/body_list returns to the
  !    client an additional pointer to the list itself.  This pointer must then
  !    be passed as an optional parameter to obs_getHeader/BodyIndex.  When a
  !    new list is requested by an OMP thread, the same physical memory is
  !    re-used for the new list.
  !
  !    In order to ensure that the same physical memory is not initially
  !    distributed to more than one OMP thread, the two small sections of
  !    IndexListDepot_mod that accomplish this are marked omp-critical.
  !

   implicit none
   save
   public

   ! methods
   public :: ild_initialize             ! initialize a list depot
   public :: ild_finalize               ! finalize a list depot
   public :: ild_get_empty_index_list   ! return a pointer to an initialized list
   public :: ild_get_next_index         ! return the next element in the list

   interface ild_get_next_index
      module procedure ild_get_next_index_depot
      module procedure ild_get_next_index_private
   end interface ild_get_next_index

                                        ! This dimension must accommodate the
                                        ! maximum number of OMP processors 
   integer, parameter :: NUMBER_OF_LISTS = 32

   type struct_index_list
      ! a list of integers, not to say indices into a struct_obs
      character(len=2) :: family        ! current_element's belong to this family
                                        ! Used only for a body list:
      integer :: header                 ! current_element's belong to this header
      integer :: current_element        ! the element that has just been returned

                                        ! the actual list of integers
                                        ! N.B.:  These are the indices that are
                                        !        being held for the client.
                                        !        In the context of this module,
                                        !        these are elements; i.e. row
                                        !        indices.  Current_element is the
                                        !        last index that was returned
                                        !        to the client.
      integer, dimension(:), allocatable :: indices
   end type struct_index_list

   type struct_index_list_depot
      ! A collection of lists, either empty or populated
                                        ! the collection of lists
      type(struct_index_list), dimension(NUMBER_OF_LISTS) :: index_lists
      integer :: list_last_attributed   ! list that was populated most recently
                                        ! list that was   used    most recently
      type(struct_index_list), pointer :: current_list
   end type struct_index_list_depot


contains

   subroutine ild_finalize( depot )
      !
      ! :Purpose: Finalize the indicated list depot
      !
      implicit none

      type(struct_index_list_depot), intent(inout) :: depot ! the depot to be finalized

      integer :: list                   ! an index

      ! Deallocate each list
      do list = 1,NUMBER_OF_LISTS
         depot%index_lists(list)%family = 'xx'
         depot%index_lists(list)%header = 0
         depot%index_lists(list)%current_element = 0
         deallocate(depot%index_lists(list)%indices)
      end do
      depot%list_last_attributed = 0
   end subroutine ild_finalize


   function ild_get_empty_index_list( depot, private_list ) &
                                                         result(empty_index_list)
      !
      ! :Purpose: From the given depot, return an index-list structure that contains
      !           no data, as a pointer.
      !           In other words, clear data from the (cyclicly) next (i.e. oldest)
      !           list and return a pointer to it.
      !
      !           If the list is being used within an OMP block, the ObsSpaceData
      !           client is responsible for holding a pointer to its own list.  This
      !           is supplied as the parameter, private_list.
      !
      implicit none

      type(struct_index_list), pointer                     :: empty_index_list  ! the returned list
      type(struct_index_list_depot), intent(inout), target :: depot             ! the depot containing the list
      type(struct_index_list), pointer, intent(inout), optional :: private_list ! used instead of depot (for OMP blocks)

      nullify(empty_index_list)

      if(present(private_list)) then
         ! This is an OMP thread
         if(associated(private_list)) then
            ! Memory has already been assigned for that thread.  Re-use it.
            empty_index_list => private_list ! Set the return pointer
         end if
      end if

      if(.not. associated(empty_index_list)) then
!$omp critical
         ! Increment (cyclicly) the index to the next list
         depot%list_last_attributed = depot%list_last_attributed + 1
         if (depot%list_last_attributed > NUMBER_OF_LISTS) &
            depot%list_last_attributed = 1

         ! Set the return pointer
         empty_index_list => depot%index_lists(depot%list_last_attributed)
!$omp end critical
      end if

      ! Initialize some values in the list
      ! empty_index_list%indices(:) = -1 --> No, the array is too big.
      empty_index_list%family = '  '
      empty_index_list%header = -1
      empty_index_list%current_element = 0

      return
   end function ild_get_empty_index_list


   function ild_get_next_index_depot(depot, no_advance) result(next_index)
      !
      ! :Purpose: From the given depot, increment the index to the current element,
      !           and return the element itself, the new current element.
      !
      implicit none

      integer :: next_index                                                ! the returned index
      type(struct_index_list_depot), intent(inout), target :: depot        ! the depot containing the list
                                                                           ! if present, do not increment
      logical, intent(in), optional                        :: no_advance   ! current_element, just return next one
      type(struct_index_list), pointer                     :: current_list ! current list of the depot
      integer                                              :: next_element ! next element of the current list

      current_list => depot%current_list
!$omp critical
                                        ! Obtain the next element from the list
      next_element = current_list%current_element + 1
      next_index = current_list%indices(next_element)
      if(.not. present(no_advance) .and. next_index /= -1) then
                                        ! Increment the current element
         current_list%current_element = next_element
      end if
!$omp end critical
   end function ild_get_next_index_depot


   function ild_get_next_index_private(private_list, no_advance) &
                                                               result(next_index)
      !
      ! :Purpose:  From the given list, increment the index to the current element, and
      !            return the element itself, the new current element.
      !
      implicit none
      integer :: next_index                                           ! the returned index
      type(struct_index_list), pointer, intent(inout) :: private_list ! the list of interest
                                                                      ! if present, do not increment
      logical, intent(in), optional                   :: no_advance   ! current_element, just return next one
      integer                                         :: next_element ! next element of the list

                                        ! Obtain the next element from the list
      next_element = private_list%current_element + 1
      next_index = private_list%indices(next_element)

      if(.not. present(no_advance) .and. next_index /= -1) then
                                        ! Increment the current element
         private_list%current_element = next_element
      end if
   end function ild_get_next_index_private


   subroutine ild_initialize(depot, numHeaderBody_max)
      !
      ! :Purpose: Initialize the indicated list depot
      !   
      ! :Note: indices is allocated with 2 extra elements to make room for
      !        the end-of-list flag that is set in
      !        obs_set_current_header/body_list
      !
      implicit none

      type(struct_index_list_depot), intent(inout) :: depot ! the depot to be initialized
                                                            ! max size of header or body of
                                                            ! struct_obs & hence of depot
      integer, intent(in) :: numHeaderBody_max
      integer :: list                   ! an index

      ! Allocate each list
      do list = 1,NUMBER_OF_LISTS
         depot%index_lists(list)%family = 'xx'
         depot%index_lists(list)%header = 0
         depot%index_lists(list)%current_element = 0
         allocate(depot%index_lists(list)%indices(numHeaderBody_max+2))
         depot%index_lists(list)%indices(:)=0
      end do
      depot%list_last_attributed = 0
      nullify(depot%current_list)
   end subroutine ild_initialize

end module IndexListDepot_mod



module ObsColumnNames_mod
  ! MODULE obsColumnNames_mod (prefix='obs' category='7. Low-level data objects and utilities')
  !
  ! :Purpose: This module simply groups together many fortran paramters that
  !           serve as column names for ObsSpaceData_mod.
  !
  ! :NOTE: This module is logistically a part of the ObsSpaceData_mod module.
  !        In fact, if fortran allowed it, ObsColumnNames_mod would be
  !        'contain'ed inside the ObsSpaceData_mod module. For this reason, and
  !        more importantly because these parameters constitute a part of the
  !        visible (from outside ObsSpaceData_mod) interface to
  !        ObsSpaceData_mod, the parameters defined in this module carry the
  !        prefix, OBS, and not CN.
  !

   public


   !
   ! INTEGER-HEADER COLUMN NUMBERS
   !  
   ! the first column index for integer header variables defined below
   ! (chosen such that every column number is unique, so that a mismatch between
   !  a column number and a column type (real, int, head, body) can be detected)
   integer, parameter :: NHDR_INT_BEG = 101
   integer, parameter, public :: OBS_RLN = NHDR_INT_BEG ! report location
                                             ! unique(within obsdat), possibly
   integer, parameter, public :: OBS_ONM = OBS_RLN+1  ! ordered, station id number
   integer, parameter, public :: OBS_INS = OBS_ONM+1  ! instrument ID  
   integer, parameter, public :: OBS_OTP = OBS_INS+1  ! observation Type (file index)
   integer, parameter, public :: OBS_ITY = OBS_OTP+1  ! code: instrument & retrieval type
   integer, parameter, public :: OBS_SAT = OBS_ITY+1  ! satellite code 
   integer, parameter, public :: OBS_TEC = OBS_SAT+1  ! satellite processing technique
   integer, parameter, public :: OBS_SEN = OBS_TEC+1  ! satellite sensor code
   integer, parameter, public :: OBS_DAT = OBS_SEN+1  ! observation date YYYYMMD
   integer, parameter, public :: OBS_ETM = OBS_DAT+1  ! observation time HHMM
   integer, parameter, public :: OBS_NLV = OBS_ETM+1  ! number of data at this location
   integer, parameter, public :: OBS_PAS = OBS_NLV+1  ! batch no. in sequential analysis
   integer, parameter, public :: OBS_REG = OBS_PAS+1  ! region number in the batch
   integer, parameter, public :: OBS_IP  = OBS_REG+1  ! number of mpi processors
   integer, parameter, public :: OBS_IPF = OBS_IP+1   ! mpi task id for file
   integer, parameter, public :: OBS_IPC = OBS_IPF+1  ! mpi task id for column/obsspacedata
   integer, parameter, public :: OBS_IPT = OBS_IPC+1  ! mpi task id for latlontile
   integer, parameter, public :: OBS_ST1 = OBS_IPT+1  ! header level status/rejection flag
   integer, parameter, public :: OBS_IDF = OBS_ST1+1  ! id. no. of observation-source file

   integer, parameter, public :: OBS_GQF = OBS_IDF+1  ! iasi GQISFLAGQUAL: general quality flag which indicates,
                                                      ! when it is set to TRUE (=1), that some anomaly has been
                                                      ! detected at some step in the Level 0 or Level 1 IASI Processing.
                                                      ! The products should not be used if this occurs.

   integer, parameter, public :: OBS_GQL = OBS_GQF+1  ! iasi GQISQUALINDEXLOC: It is defined as the uncertainty of the
                                                      ! coregistration between IIS (IASI internal imager) and AVHRR.

   integer, parameter, public :: OBS_NCO2= OBS_GQL+1  ! NCO2: number of valid CO2 slicing estimates (AIRS,IASI,CrIS)
   integer, parameter, public :: OBS_STYP= OBS_NCO2+1 ! surface type in obs file (0,1,2)
   integer, parameter, public :: OBS_ROQF= OBS_STYP+1 ! QUALITY FLAGS FOR RADIO OCCULTATION DATA
   integer, parameter, public :: OBS_SWQ1= OBS_ROQF+1 ! QUALITY VALUES FOR SW (AMV) DATA
   integer, parameter, public :: OBS_SWQ2= OBS_SWQ1+1 ! QUALITY VALUES FOR SW (AMV) DATA
   integer, parameter, public :: OBS_SWMT= OBS_SWQ2+1 ! QUALITY VALUES FOR SW (AMV) DATA
   integer, parameter, public :: OBS_SWLS= OBS_SWMT+1 ! QUALITY VALUES FOR SW (AMV) DATA
   integer, parameter, public :: OBS_SWGA= OBS_SWLS+1 ! QUALITY VALUES FOR SW (AMV) DATA
   integer, parameter, public :: OBS_SWHA= OBS_SWGA+1 ! QUALITY VALUES FOR SW (AMV) DATA
   integer, parameter, public :: OBS_CHM = OBS_SWHA+1 ! BUFR code (table 08046) of constituent type (for the CH family)
   integer, parameter, public :: OBS_FOV = OBS_CHM+1  ! field of view 
   integer, parameter, public :: OBS_PRFL= OBS_FOV+1  ! profile id. number
   integer, parameter, public :: OBS_PHAS= OBS_PRFL+1 ! phase of flight
   integer, parameter, public :: OBS_ORI = OBS_PHAS+1 ! originating centre code
   integer, parameter, public :: OBS_LCH = OBS_ORI+1  ! launch time (hhmm)
   integer, parameter, public :: OBS_RTP = OBS_LCH+1  ! radiosonde type code
   integer, parameter, public :: OBS_HDD = OBS_RTP+1  ! date in burp header
   integer, parameter, public :: OBS_HDT = OBS_HDD+1  ! time in burp header
   integer, parameter, public :: OBS_TFLG= OBS_HDT+1  ! flag for hi-res time element
   integer, parameter, public :: OBS_LFLG= OBS_TFLG+1 ! flag for hi-res lat element
   integer, parameter, public :: OBS_ORBI= OBS_LFLG+1! satellite orbit index
   integer, parameter, public :: OBS_AQF1= OBS_ORBI+1! ATMS Geolocalisation Quality Control Flag 
   integer, parameter, public :: OBS_AQF2= OBS_AQF1+1! ATMS Granule Level Quality Control Flag 
   integer, parameter, public :: OBS_AQF3= OBS_AQF2+1! ATMS Scan Level Quality Control Flag 
   integer, parameter, public :: OBS_TTYP= OBS_AQF3+1! TERRAIN TYP INDICE for TOVS QC
   integer, parameter, public :: OBS_INFG= OBS_TTYP+1! SATPLOT INFO FLAG for TOVS
   integer, parameter, public :: OBS_RAIN= OBS_INFG+1! UKMO rain flag for ssmis obs
   integer, parameter, public :: OBS_CHID= OBS_RAIN+1! id. no. of ice chart


   ! the last column index for integer header variables defined just above
   integer, parameter :: NHDR_INT_END = OBS_CHID

   integer, parameter :: NHDR_INT_SIZE = NHDR_INT_END - NHDR_INT_BEG + 1

   !
   ! INTEGER-HEADER COLUMN NAMES
   !
   character(len=4), target :: ocn_ColumnNameList_IH(NHDR_INT_BEG:NHDR_INT_END) = &
      (/ 'RLN ','ONM ','INS ','OTP ','ITY ','SAT ','TEC ','SEN ','DAT ','ETM ', &  
         'NLV ','PAS ','REG ','IP  ','IPF ','IPC ','IPT ', &  
         'ST1 ','IDF ','GQF ','GQL ','NCO2','STYP','ROQF', &
         'SWQ1','SWQ2','SWMT','SWLS','SWGA','SWHA','CHM ','FOV ', &
         'PRFL','PHAS','ORI ','LCH ','RTP ','HDD ','HDT ','TFLG',&
         'LFLG','ORBI','AQF1','AQF2','AQF3','TTYP','INFG','RAIN','CHID'/)

   !
   ! REAL-HEADER COLUMN NUMBERS
   !
   ! the first column index for real header variables defined below
   ! (chosen such that every column number is unique, so that a mismatch between
   !  a column number and a column type (real, int, head, body) can be detected)
   integer, parameter :: NHDR_REAL_BEG = 201
   integer, parameter, public :: OBS_LAT = NHDR_REAL_BEG  ! latitude  in radians (N positive)
   integer, parameter, public :: OBS_LON = OBS_LAT+1 ! longitude in radians (E positive)
   integer, parameter, public :: OBS_ALT = OBS_LON+1 ! station altitude
   integer, parameter, public :: OBS_BX  = OBS_ALT+1 ! x-coordinate of block in R3
   integer, parameter, public :: OBS_BY  = OBS_BX +1 ! y-coordinate of block in R3
   integer, parameter, public :: OBS_BZ  = OBS_BY +1 ! z-coordinate of block in R3

   integer, parameter, public :: OBS_M1C1  = OBS_BZ +1   ! mean for class 1 AVHRR channel 1
   integer, parameter, public :: OBS_M1C2  = OBS_M1C1 +1 ! mean for class 1 AVHRR channel 2
   integer, parameter, public :: OBS_M1C3  = OBS_M1C2 +1 ! mean for class 1 AVHRR channel 3
   integer, parameter, public :: OBS_M1C4  = OBS_M1C3 +1 ! mean for class 1 AVHRR channel 4
   integer, parameter, public :: OBS_M1C5  = OBS_M1C4 +1 ! mean for class 1 AVHRR channel 5
   integer, parameter, public :: OBS_M1C6  = OBS_M1C5 +1 ! mean for class 1 AVHRR channel 6

   integer, parameter, public :: OBS_M2C1  = OBS_M1C6 +1 ! mean for class 2 AVHRR channel 1
   integer, parameter, public :: OBS_M2C2  = OBS_M2C1 +1 ! mean for class 2 AVHRR channel 2
   integer, parameter, public :: OBS_M2C3  = OBS_M2C2 +1 ! mean for class 2 AVHRR channel 3
   integer, parameter, public :: OBS_M2C4  = OBS_M2C3 +1 ! mean for class 2 AVHRR channel 4
   integer, parameter, public :: OBS_M2C5  = OBS_M2C4 +1 ! mean for class 2 AVHRR channel 5
   integer, parameter, public :: OBS_M2C6  = OBS_M2C5 +1 ! mean for class 2 AVHRR channel 6

   integer, parameter, public :: OBS_M3C1  = OBS_M2C6 +1 ! mean for class 3 AVHRR channel 1
   integer, parameter, public :: OBS_M3C2  = OBS_M3C1 +1 ! mean for class 3 AVHRR channel 2
   integer, parameter, public :: OBS_M3C3  = OBS_M3C2 +1 ! mean for class 3 AVHRR channel 3
   integer, parameter, public :: OBS_M3C4  = OBS_M3C3 +1 ! mean for class 3 AVHRR channel 4
   integer, parameter, public :: OBS_M3C5  = OBS_M3C4 +1 ! mean for class 3 AVHRR channel 5
   integer, parameter, public :: OBS_M3C6  = OBS_M3C5 +1 ! mean for class 3 AVHRR channel 6

   integer, parameter, public :: OBS_M4C1  = OBS_M3C6 +1 ! mean for class 4 AVHRR channel 1
   integer, parameter, public :: OBS_M4C2  = OBS_M4C1 +1 ! mean for class 4 AVHRR channel 2
   integer, parameter, public :: OBS_M4C3  = OBS_M4C2 +1 ! mean for class 4 AVHRR channel 3
   integer, parameter, public :: OBS_M4C4  = OBS_M4C3 +1 ! mean for class 4 AVHRR channel 4
   integer, parameter, public :: OBS_M4C5  = OBS_M4C4 +1 ! mean for class 4 AVHRR channel 5
   integer, parameter, public :: OBS_M4C6  = OBS_M4C5 +1 ! mean for class 4 AVHRR channel 6

   integer, parameter, public :: OBS_M5C1  = OBS_M4C6 +1 ! mean for class 5 AVHRR channel 1
   integer, parameter, public :: OBS_M5C2  = OBS_M5C1 +1 ! mean for class 5 AVHRR channel 2
   integer, parameter, public :: OBS_M5C3  = OBS_M5C2 +1 ! mean for class 5 AVHRR channel 3
   integer, parameter, public :: OBS_M5C4  = OBS_M5C3 +1 ! mean for class 5 AVHRR channel 4
   integer, parameter, public :: OBS_M5C5  = OBS_M5C4 +1 ! mean for class 5 AVHRR channel 5
   integer, parameter, public :: OBS_M5C6  = OBS_M5C5 +1 ! mean for class 5 AVHRR channel 6

   integer, parameter, public :: OBS_M6C1  = OBS_M5C6 +1 ! mean for class 6 AVHRR channel 1
   integer, parameter, public :: OBS_M6C2  = OBS_M6C1 +1 ! mean for class 6 AVHRR channel 2
   integer, parameter, public :: OBS_M6C3  = OBS_M6C2 +1 ! mean for class 6 AVHRR channel 3
   integer, parameter, public :: OBS_M6C4  = OBS_M6C3 +1 ! mean for class 6 AVHRR channel 4
   integer, parameter, public :: OBS_M6C5  = OBS_M6C4 +1 ! mean for class 6 AVHRR channel 5
   integer, parameter, public :: OBS_M6C6  = OBS_M6C5 +1 ! mean for class 6 AVHRR channel 6

   integer, parameter, public :: OBS_M7C1  = OBS_M6C6 +1 ! mean for class 7 AVHRR channel 1
   integer, parameter, public :: OBS_M7C2  = OBS_M7C1 +1 ! mean for class 7 AVHRR channel 2
   integer, parameter, public :: OBS_M7C3  = OBS_M7C2 +1 ! mean for class 7 AVHRR channel 3
   integer, parameter, public :: OBS_M7C4  = OBS_M7C3 +1 ! mean for class 7 AVHRR channel 4
   integer, parameter, public :: OBS_M7C5  = OBS_M7C4 +1 ! mean for class 7 AVHRR channel 5
   integer, parameter, public :: OBS_M7C6  = OBS_M7C5 +1 ! mean for class 7 AVHRR channel 6


   integer, parameter, public :: OBS_S1C1  = OBS_M7C6 +1 ! stdev for class 1 AVHRR channel 1
   integer, parameter, public :: OBS_S1C2  = OBS_S1C1 +1 ! stdev for class 1 AVHRR channel 2
   integer, parameter, public :: OBS_S1C3  = OBS_S1C2 +1 ! stdev for class 1 AVHRR channel 3
   integer, parameter, public :: OBS_S1C4  = OBS_S1C3 +1 ! stdev for class 1 AVHRR channel 4
   integer, parameter, public :: OBS_S1C5  = OBS_S1C4 +1 ! stdev for class 1 AVHRR channel 5
   integer, parameter, public :: OBS_S1C6  = OBS_S1C5 +1 ! stdev for class 1 AVHRR channel 6

   integer, parameter, public :: OBS_S2C1  = OBS_S1C6 +1 ! stdev for class 2 AVHRR channel 1
   integer, parameter, public :: OBS_S2C2  = OBS_S2C1 +1 ! stdev for class 2 AVHRR channel 2
   integer, parameter, public :: OBS_S2C3  = OBS_S2C2 +1 ! stdev for class 2 AVHRR channel 3
   integer, parameter, public :: OBS_S2C4  = OBS_S2C3 +1 ! stdev for class 2 AVHRR channel 4
   integer, parameter, public :: OBS_S2C5  = OBS_S2C4 +1 ! stdev for class 2 AVHRR channel 5
   integer, parameter, public :: OBS_S2C6  = OBS_S2C5 +1 ! stdev for class 2 AVHRR channel 6

   integer, parameter, public :: OBS_S3C1  = OBS_S2C6 +1 ! stdev for class 3 AVHRR channel 1
   integer, parameter, public :: OBS_S3C2  = OBS_S3C1 +1 ! stdev for class 3 AVHRR channel 2
   integer, parameter, public :: OBS_S3C3  = OBS_S3C2 +1 ! stdev for class 3 AVHRR channel 3
   integer, parameter, public :: OBS_S3C4  = OBS_S3C3 +1 ! stdev for class 3 AVHRR channel 4
   integer, parameter, public :: OBS_S3C5  = OBS_S3C4 +1 ! stdev for class 3 AVHRR channel 5
   integer, parameter, public :: OBS_S3C6  = OBS_S3C5 +1 ! stdev for class 3 AVHRR channel 6

   integer, parameter, public :: OBS_S4C1  = OBS_S3C6 +1 ! stdev for class 4 AVHRR channel 1
   integer, parameter, public :: OBS_S4C2  = OBS_S4C1 +1 ! stdev for class 4 AVHRR channel 2
   integer, parameter, public :: OBS_S4C3  = OBS_S4C2 +1 ! stdev for class 4 AVHRR channel 3
   integer, parameter, public :: OBS_S4C4  = OBS_S4C3 +1 ! stdev for class 4 AVHRR channel 4
   integer, parameter, public :: OBS_S4C5  = OBS_S4C4 +1 ! stdev for class 4 AVHRR channel 5
   integer, parameter, public :: OBS_S4C6  = OBS_S4C5 +1 ! stdev for class 4 AVHRR channel 6

   integer, parameter, public :: OBS_S5C1  = OBS_S4C6 +1 ! stdev for class 5 AVHRR channel 1
   integer, parameter, public :: OBS_S5C2  = OBS_S5C1 +1 ! stdev for class 5 AVHRR channel 2
   integer, parameter, public :: OBS_S5C3  = OBS_S5C2 +1 ! stdev for class 5 AVHRR channel 3
   integer, parameter, public :: OBS_S5C4  = OBS_S5C3 +1 ! stdev for class 5 AVHRR channel 4
   integer, parameter, public :: OBS_S5C5  = OBS_S5C4 +1 ! stdev for class 5 AVHRR channel 5
   integer, parameter, public :: OBS_S5C6  = OBS_S5C5 +1 ! stdev for class 5 AVHRR channel 6

   integer, parameter, public :: OBS_S6C1  = OBS_S5C6 +1 ! stdev for class 6 AVHRR channel 1
   integer, parameter, public :: OBS_S6C2  = OBS_S6C1 +1 ! stdev for class 6 AVHRR channel 2
   integer, parameter, public :: OBS_S6C3  = OBS_S6C2 +1 ! stdev for class 6 AVHRR channel 3
   integer, parameter, public :: OBS_S6C4  = OBS_S6C3 +1 ! stdev for class 6 AVHRR channel 4
   integer, parameter, public :: OBS_S6C5  = OBS_S6C4 +1 ! stdev for class 6 AVHRR channel 5
   integer, parameter, public :: OBS_S6C6  = OBS_S6C5 +1 ! stdev for class 6 AVHRR channel 6

   integer, parameter, public :: OBS_S7C1  = OBS_S6C6 +1 ! stdev for class 7 AVHRR channel 1
   integer, parameter, public :: OBS_S7C2  = OBS_S7C1 +1 ! stdev for class 7 AVHRR channel 2
   integer, parameter, public :: OBS_S7C3  = OBS_S7C2 +1 ! stdev for class 7 AVHRR channel 3
   integer, parameter, public :: OBS_S7C4  = OBS_S7C3 +1 ! stdev for class 7 AVHRR channel 4
   integer, parameter, public :: OBS_S7C5  = OBS_S7C4 +1 ! stdev for class 7 AVHRR channel 5
   integer, parameter, public :: OBS_S7C6  = OBS_S7C5 +1 ! stdev for class 7 AVHRR channel 6
   integer, parameter, public :: OBS_CF1   = OBS_S7C6 +1 ! AVHRR fraction of class 1
   integer, parameter, public :: OBS_CF2   = OBS_CF1  +1 ! AVHRR fraction of class 2
   integer, parameter, public :: OBS_CF3   = OBS_CF2  +1 ! AVHRR fraction of class 3
   integer, parameter, public :: OBS_CF4   = OBS_CF3  +1 ! AVHRR fraction of class 4
   integer, parameter, public :: OBS_CF5   = OBS_CF4  +1 ! AVHRR fraction of class 5
   integer, parameter, public :: OBS_CF6   = OBS_CF5  +1 ! AVHRR fraction of class 6
   integer, parameter, public :: OBS_CF7   = OBS_CF6  +1 ! AVHRR fraction of class 7
   integer, parameter, public :: OBS_ETOP  = OBS_CF7  +1 ! CO2 slicing consensus (median) cloud top pressure
   integer, parameter, public :: OBS_VTOP  = OBS_ETOP +1 ! estimated error on CO2 slicing cloud top pressure
   integer, parameter, public :: OBS_ECF   = OBS_VTOP +1 ! CO2 slicing effective cloud fraction 
   integer, parameter, public :: OBS_VCF   = OBS_ECF  +1 ! estimated error on CO2 CO2 slicing cloud  fraction 
   integer, parameter, public :: OBS_HE    = OBS_VCF  +1 ! cloud effective height (one channel)
   integer, parameter, public :: OBS_ZTSR  = OBS_HE   +1 ! retrieved skin temperature from window channel in K
   integer, parameter, public :: OBS_ZTM   = OBS_ZTSR +1 ! model temperature, eta=1, in K (should not be there)
   integer, parameter, public :: OBS_ZTGM  = OBS_ZTM  +1 ! surface model temperature (skin) in K
   integer, parameter, public :: OBS_ZLQM  = OBS_ZTGM +1 ! specific humidity at surface (2m) in kg/kg
   integer, parameter, public :: OBS_ZPS   = OBS_ZLQM +1 ! surface model pressure in Pa
   integer, parameter, public :: OBS_TRAD  = OBS_ZPS  +1 ! Local EARTH Radius Metres
   integer, parameter, public :: OBS_GEOI  = OBS_TRAD +1 ! Geoid Undulation  Metres
   integer, parameter, public :: OBS_CLF   = OBS_GEOI +1 ! cloud fraction
   integer, parameter, public :: OBS_SUN   = OBS_CLF  +1 ! sun zenith angle
   integer, parameter, public :: OBS_SZA   = OBS_SUN  +1 ! satellite zenith angle
   integer, parameter, public :: OBS_AZA   = OBS_SZA  +1 ! satellite azimuthal angle
   integer, parameter, public :: OBS_SAZ   = OBS_AZA  +1 ! sun azimuth angle
   integer, parameter, public :: OBS_CLWO  = OBS_SAZ  +1 ! cloud liquid water retrieved from observation
   integer, parameter, public :: OBS_CLWB  = OBS_CLWO +1 ! cloud liquid water retrieved from background
   integer, parameter, public :: OBS_MWS   = OBS_CLWB +1 ! model wind speed (in ASCAT data)
   integer, parameter, public :: OBS_SCAT  = OBS_MWS  +1 ! atmospheric scatering index
   integer, parameter, public :: OBS_IWV   = OBS_SCAT +1 ! atmospheric integrated water vapor for ssmis
   integer, parameter, public :: OBS_RZAM  = OBS_IWV  +1 ! Azimuth of the Radar beam   (radians)
   integer, parameter, public :: OBS_RELE  = OBS_RZAM +1 ! Elevation of the Radar beam (radians)
   integer, parameter, public :: OBS_RANS  = OBS_RELE +1 ! Initial range of the Radar beam 
   integer, parameter, public :: OBS_RANE  = OBS_RANS +1 ! Final range of the Radar beam 

   ! the last column index for real header variables defined just above
   integer, parameter :: NHDR_REAL_END = OBS_RANE 
   integer, parameter :: NHDR_REAL_SIZE = NHDR_REAL_END - NHDR_REAL_BEG + 1

   !
   ! REAL-HEADER COLUMN NAMES
   !
   character(len=4), target :: ocn_ColumnNameList_RH(NHDR_REAL_BEG:NHDR_REAL_END) =  &
      (/'LAT ','LON ','ALT ','BX  ','BY  ','BZ  ', &
        'M1C1','M1C2','M1C3','M1C4','M1C5','M1C6', &
        'M2C1','M2C2','M2C3','M2C4','M2C5','M2C6', &
        'M3C1','M3C2','M3C3','M3C4','M3C5','M3C6', &
        'M4C1','M4C2','M4C3','M4C4','M4C5','M4C6', &
        'M5C1','M5C2','M5C3','M5C4','M5C5','M5C6', &
        'M6C1','M6C2','M6C3','M6C4','M6C5','M6C6', &
        'M7C1','M7C2','M7C3','M7C4','M7C5','M7C6', &
        'S1C1','S1C2','S1C3','S1C4','S1C5','S1C6', &
        'S2C1','S2C2','S2C3','S2C4','S2C5','S2C6', &
        'S3C1','S3C2','S3C3','S3C4','S3C5','S3C6', &
        'S4C1','S4C2','S4C3','S4C4','S4C5','S4C6', &
        'S5C1','S5C2','S5C3','S5C4','S5C5','S5C6', &
        'S6C1','S6C2','S6C3','S6C4','S6C5','S6C6', &
        'S7C1','S7C2','S7C3','S7C4','S7C5','S7C6', &
        'CF1 ','CF2 ','CF3 ','CF4 ','CF5 ','CF6 ', &
        'CF7 ','ETOP','VTOP','ECF ','VCF ','HE  ', &
        'ZTSR','ZTM ','ZTGM','ZLQM','ZPS ','TRAD', &
        'GEOI','CLF ','SUN ','SZA ','AZA ','SAZ ', &
        'CLW1','CLW2','MWS ','SCAT','IWV','RZAM',   &
        'RELE','RANS','RANE'/)
   !
   ! INTEGER-BODY COLUMN NUMBERS
   !
   ! the first column index for integer body variables defined below
   ! (chosen such that every column number is unique, so that a mismatch between
   !  a column number and a column type (real, int, head, body) can be detected)
   integer, parameter :: NBDY_INT_BEG = 401
   integer, parameter, public :: OBS_VNM = NBDY_INT_BEG ! variable number
   integer, parameter, public :: OBS_FLG = OBS_VNM+1  ! flags
   integer, parameter, public :: OBS_KFA = OBS_FLG+1  ! marker for forward interp problems
   integer, parameter, public :: OBS_ASS = OBS_KFA+1  ! flag to indicate if assimilated
   integer, parameter, public :: OBS_HIND= OBS_ASS+1  ! corresponding header row index
   integer, parameter, public :: OBS_VCO = OBS_HIND+1 ! type of vertical coordinate
   integer, parameter, public :: OBS_LYR = OBS_VCO+1  ! Index of anal level above observ'n
                                             ! Flag: extrapolation necessary of
   integer, parameter, public :: OBS_XTR = OBS_LYR+1  ! anal variables to obs'n location
   integer, parameter, public :: OBS_QCF2= OBS_XTR+1  ! TOVs Data Level Qc Flag
                                                      ! flag values for btyp=9248 block ele 033081
   integer, parameter, public :: OBS_CLA= OBS_QCF2+1  ! csr cloud amount btyp 9248 ele 020081

   ! the last column index for integer body variables defined just above
   integer, parameter :: NBDY_INT_END = OBS_CLA
   integer, parameter :: NBDY_INT_SIZE = NBDY_INT_END - NBDY_INT_BEG + 1

   !
   ! INTEGER-BODY COLUMN NAMES
   !
   character(len=4), target :: ocn_ColumnNameList_IB(NBDY_INT_BEG:NBDY_INT_END) = &
      (/ 'VNM ','FLG ','KFA ','ASS ','HIND','VCO ','LYR ','XTR ', 'QCFL ', &
          'CLA '/)  

   !
   ! REAL-BODY COLUMN NUMBERS
   !
   ! the first column index for real body variables defined below
   ! (chosen such that every column number is unique, so that a mismatch between
   !  a column number and a column type (real, int, head, body) can be detected)
   integer, parameter :: NBDY_REAL_BEG   = 501
   integer, parameter, public :: OBS_PPP = NBDY_REAL_BEG ! vertical coordinate
   integer, parameter, public :: OBS_SEM = OBS_PPP +1 ! surface emissivity
   integer, parameter, public :: OBS_VAR = OBS_SEM +1 ! value of the observation
   integer, parameter, public :: OBS_OMP = OBS_VAR +1 ! obs - H (trial field)
   integer, parameter, public :: OBS_OMA = OBS_OMP +1 ! obs - H (analysis)
   integer, parameter, public :: OBS_OER = OBS_OMA +1 ! sigma(obs)
   integer, parameter, public :: OBS_HPHT= OBS_OER +1 ! root of (hpht with hx scalar)
   integer, parameter, public :: OBS_HAHT= OBS_HPHT +1 ! root of (hp_{a}ht with hx scalar)
   integer, parameter, public :: OBS_ZHA = OBS_HAHT+1 ! vert coordinate for Schur product
   integer, parameter, public :: OBS_OMP6= OBS_ZHA +1 ! obs - H (6-h trial field)
   integer, parameter, public :: OBS_OMA0= OBS_OMP6+1 ! obs - H (analysis at central time)
   integer, parameter, public :: OBS_SIGI= OBS_OMA0+1 ! ensemble-based estimate of the innov std dev
   integer, parameter, public :: OBS_SIGO= OBS_SIGI+1 ! ensemble-based estimate of obs std dev
   integer, parameter, public :: OBS_POB = OBS_SIGO+1 ! initial value of "gamma" for variational QC 
   integer, parameter, public :: OBS_WORK= OBS_POB +1 ! temporary values
   integer, parameter, public :: OBS_PRM = OBS_WORK+1 ! (adjusted) observed value for tovs in variational assimilation
   integer, parameter, public :: OBS_JOBS= OBS_PRM +1 ! contribution to obs cost function
   integer, parameter, public :: OBS_QCV = OBS_JOBS+1 ! weight-reduction factor for var QC
   integer, parameter, public :: OBS_FSO = OBS_QCV+1  ! forecast sensitivity to observations
   integer, parameter, public :: OBS_CRPS= OBS_FSO+1  ! Continuous Ranked Probability Score
   integer, parameter, public :: OBS_BCOR= OBS_CRPS+1 ! observation bias correction
   integer, parameter, public :: OBS_OMPE= OBS_BCOR+1 ! error standard deviation of [obs - H (trial field)]
   integer, parameter, public :: OBS_LATD= OBS_OMPE+1 ! obs LATitude  in Data table (radians)
   integer, parameter, public :: OBS_LOND= OBS_LATD+1 ! obs LONgitude in Data table (radians)
   integer, parameter, public :: OBS_BTCL= OBS_LOND+1 ! clear-sky simulated observation
   integer, parameter, public :: OBS_LOCI= OBS_BTCL+1 ! LOCation Information for observation (e.g. range along radar beam)

   ! the number of real body variables defined just above
   integer, parameter :: NBDY_REAL_END = OBS_LOCI
   integer, parameter :: NBDY_REAL_SIZE = NBDY_REAL_END - NBDY_REAL_BEG + 1

   !
   ! REAL-BODY COLUMN NAMES
   !
   character(len=4), target :: ocn_ColumnNameList_RB(NBDY_REAL_BEG:NBDY_REAL_END) = &
      (/ 'PPP ','SEM ','VAR ','OMP ','OMA ','OER ','HPHT','HAHT','ZHA ','OMP6','OMA0',     &
         'SIGI','SIGO','POB ','WORK','PRM ','JOBS','QCV ','FSO ','CRPS','BCOR','OMPE',     &
         'ROLA','ROLO','VAR2','LOCI' /)
end module ObsColumnNames_mod



module ObsDataColumn_mod
  !
  ! MODULE obsDataColumn_mod (prefix='odc' category='7. Low-level data objects and utilities')
  !
  ! :Purpose:
  !    This module is used exclusively by the obsSpaceData module which follows
  !    in this file. The derived type is used to represent a "column" of
  !    observation data in an instance of the struct_obs defined in obsSpaceData.
  !    It contains a pointer for each possible type of data stored in a column,
  !    but only one should be allocated at any time.
  !
   use codePrecision_mod
   use ObsColumnNames_mod
   implicit none
   save
   private

   ! CLASS-CONSTANT:
   ! CLASS-CONSTANT:
   ! CLASS-CONSTANT:

   ! This type gathers together into one structure the various CLASS-CONSTANT
   ! characteristics of a data column.  Four instances (flavours) of this derived
   ! type are defined below.
   type, public :: struct_odc_flavour
                                        ! These 2 values are informational only.
                                        ! They are used in error messages.
      character(len=4) :: dataType      ! REAL or INT
      character(len=4) :: headOrBody    ! HEAD or BODY

      integer :: ncol_beg, ncol_end
      logical         , dimension(:), pointer :: columnActive   ! indexed from 1
      character(len=4), dimension(:), pointer :: columnNameList

      integer,dimension(:), pointer ::activeIndexFromColumnIndex
      logical :: activeIndexFromColumnIndex_defined

      integer,dimension(:), pointer ::columnIndexFromActiveIndex
      logical :: columnIndexFromActiveIndex_defined
   end type struct_odc_flavour

   ! The four CLASS-CONSTANT flavours of data columns:
   !    Integer / Real  +  Body / Header
   type(struct_odc_flavour), public, target :: odc_flavour_IB, &
                                               odc_flavour_IH, &
                                               odc_flavour_RB, &
                                               odc_flavour_RH
   ! end of CLASS-CONSTANT objects
   ! end of CLASS-CONSTANT objects
   ! end of CLASS-CONSTANT objects


   ! methods
   public :: odc_allocate, odc_deallocate, odc_class_initialize
   public :: odc_activateColumn, odc_numActiveColumn
   public :: odc_columnElem, odc_columnSet
   public :: odc_columnIndexFromActiveIndex, odc_activeIndexFromColumnIndex

   ! This type allows a single derived type to contain either a real or an int
   type, public :: struct_obsDataColumn
      logical          :: allocated = .false.
                                        ! For these arrays:
                                        !   1st dim'n:  row index (element index)
      integer,          pointer :: value_i(:)     => NULL()
      real(pre_obsReal),pointer :: value_r(:)     => NULL()
      character(len=4) :: dataType
   end type struct_obsDataColumn


   ! This type contains one array of data columns.  Four of these are necessary
   ! to constitute a complete set of observation data (struct_obs).
   type, public :: struct_obsDataColumn_Array
                                        ! CLASS-CONSTANT values (1 of 4 flavours)
      type(struct_odc_flavour), pointer :: odc_flavour => NULL()
                                        ! object-specific values
                                        !   1st dim'n: column index (column name)
      type(struct_obsDataColumn), dimension(:), pointer :: columns => NULL()
   end type struct_obsDataColumn_Array
  

   ! These arrays store the status of the columns.  An active column is 
   ! allocated (and can therefore be used) in any object that is instantiated.
   logical, target :: columnActive_IH(NHDR_INT_BEG:NHDR_INT_END ) = .false.
   logical, target :: columnActive_RH(NHDR_REAL_BEG:NHDR_REAL_END) = .false.
   logical, target :: columnActive_IB(NBDY_INT_BEG:NBDY_INT_END ) = .false.
   logical, target :: columnActive_RB(NBDY_REAL_BEG:NBDY_REAL_END) = .false.

   integer, target :: activeIndexFromColumnIndex_IB(NBDY_INT_BEG:NBDY_INT_END)
   integer, target :: activeIndexFromColumnIndex_IH(NHDR_INT_BEG:NHDR_INT_END)
   integer, target :: activeIndexFromColumnIndex_RB(NBDY_REAL_BEG:NBDY_REAL_END)
   integer, target :: activeIndexFromColumnIndex_RH(NHDR_REAL_BEG:NHDR_REAL_END)

   integer, target :: columnIndexFromActiveIndex_IB(NBDY_INT_SIZE)
   integer, target :: columnIndexFromActiveIndex_IH(NHDR_INT_SIZE)
   integer, target :: columnIndexFromActiveIndex_RB(NBDY_REAL_SIZE)
   integer, target :: columnIndexFromActiveIndex_RH(NHDR_REAL_SIZE)

   integer, public, parameter :: odc_ENKF_bdy_int_column_list(8) = &
      (/OBS_VNM, OBS_FLG, OBS_ASS, OBS_HIND, OBS_VCO, OBS_LYR, &
        OBS_QCF2, OBS_CLA /)
   integer, public, parameter :: odc_ENKF_bdy_real_column_list(15) = &
      (/OBS_PPP, OBS_SEM, OBS_VAR,  OBS_OMP,  OBS_OMA,  OBS_OER,  OBS_HPHT,&
        OBS_HAHT,OBS_ZHA, OBS_OMP6, OBS_OMA0, OBS_SIGI, OBS_SIGO, OBS_LATD,&
        OBS_LOND /)

contains

   subroutine odc_abort(cdmessage)
     !
     ! :Purpose: Abort a job on error (same as OBS_ABORT)
     !
     ! :Arguments:
     !           :cdmessage: message to be printed
     !
     ! :Note: For debugging (i.e. UNIT_TESTING is defined), obs_abort should
     !        generally be followed by a 'return' in the calling routine.

!#if defined(UNIT_TESTING)
!      use pFUnit
!#endif
      implicit none
      character(len=*), intent(in) :: cdmessage

      write(*,'(//,4X,"ABORTING IN ObsDataColumn_mod:-------",/,8X,A)')cdmessage
      call flush(6)

!#if defined(UNIT_TESTING)
!      call throw(Exception('exiting in odc_abort:' // cdmessage))
!#else
      call qqexit(1)
      stop
!#endif
   end subroutine odc_abort


   function odc_activeIndexFromColumnIndex(odc_flavour,column_index_in, &
                                           recompute) result(active_index_out)
      !
      ! :Purpose: The list of active columns is only a subset of all possible
      !           columns.  Return the index into the list of active columns, given
      !           the index into the list of all columns.
      !
      implicit none
      type(struct_odc_flavour), intent(inout) :: odc_flavour
      integer                 , intent(in)    :: column_index_in
      logical, optional       , intent(in)    :: recompute
      integer                                 :: active_index_out

      integer :: active_index, &
                 column_index
      character(len=100) :: message

      if(present(recompute)) then
         if(recompute) odc_flavour%activeIndexFromColumnIndex_defined=.false.
      endif

      if(.not. odc_flavour%activeIndexFromColumnIndex_defined) then
         odc_flavour%activeIndexFromColumnIndex_defined=.true.
         active_index=0
         odc_flavour%activeIndexFromColumnIndex(:)=-1

         do column_index = odc_flavour%ncol_beg, odc_flavour%ncol_end
            if(odc_flavour%columnActive(column_index)) then
               active_index=active_index+1
               odc_flavour%activeIndexFromColumnIndex(column_index) =active_index
            endif
         enddo
      endif

      active_index_out=odc_flavour%activeIndexFromColumnIndex(column_index_in)

      if(active_index_out == -1) then
         write(message,*)'ODC_activeIndexFromColumnIndex: requested column is ',&
                         'not active!  Column name is ', &
                         odc_flavour%columnNameList(column_index_in)
         call odc_abort(message)
      end if
   end function odc_activeIndexFromColumnIndex


   subroutine odc_allocate( odc, numRows, name, dataType, scratchReal, scratchInt )
     !
     ! :Purpose: Allocate a single column of obs data according to 
     !           specifications in input arguments
     !
     ! :Arguments:
     !           :odc:        instance of the obsDataColumn type
     !           :numRows:    number of column rows to allocate
     !           :name:       character string name of column
     !           :dataType:   character string type of column data: REAL or INT
     !           :headOrBody: character string indicating HEAD or BODY
     !
     !
      implicit none
      type(struct_obsDataColumn), intent(inout) :: odc
      integer, intent(in) :: numRows
      character(len=*), intent(in) :: name,dataType
      real(pre_obsReal), pointer, intent(in) :: scratchReal(:)
      integer          , pointer, intent(in) :: scratchInt(:)

      if(odc%allocated) then
         call odc_abort('ODC_ALLOCATE: column is already allocated. name=' &
                        // name)
         return
      endif

      odc%allocated=.true.
      odc%dataType = dataType

      select case (trim(dataType))
      case ('INT')
         allocate(odc%value_i(numRows))
         odc%value_i(:)=0
         odc%value_r   => scratchReal

      case ('REAL')
         allocate(odc%value_r(numRows))
         odc%value_r(:)=real(0.0D0, pre_obsReal)
         odc%value_i   => scratchInt

      case default
         call odc_abort('ODC_ALLOCATE: unknown data type. type=' // dataType)
      end select

   end subroutine odc_allocate


   subroutine odc_activateColumn(odc_flavour, column_index)
      !
      ! :Purpose: Set the 'active' flag for the indicated column.  This enables memory
      !           allocation for this column without actually allocating the memory.
      !
      implicit none
      type(struct_odc_flavour), intent(inout)  :: odc_flavour
      integer, intent(in) :: column_index

      integer :: active_index, dummy_index

      if(.not.odc_flavour%columnActive(column_index)) then
         odc_flavour%columnActive(column_index) = .true.
      endif

      ! force the recalculation of indices to go between activeColumnIndex and
      ! columnIndex
      active_index=odc_activeIndexFromColumnIndex(odc_flavour,column_index, &
                                                  recompute=.true.)
      dummy_index =odc_columnIndexFromActiveIndex(odc_flavour,active_index, &
                                                  recompute=.true.)
   end subroutine odc_activateColumn


   subroutine odc_initColumnFlavour(odc_flavour, dataType_in, headOrBody_in)
      !
      ! :Purpose: Set pointers according to the four column flavours (header /
      !           body, integer / real).
      !      
      type(struct_odc_flavour), intent(inout) :: odc_flavour
      character(len=*)        , intent(in)    :: dataType_in    ! REAL or INT
      character(len=*)        , intent(in)    :: headOrBody_in  ! HEAD or BODY

      odc_flavour%dataType   = trim(dataType_in)
      odc_flavour%headOrBody = trim(headOrBody_in)

      select case (trim(dataType_in))
      case ('REAL')
         select case (trim(headOrBody_in))
         case ('HEAD')
            odc_flavour%ncol_beg = NHDR_REAL_BEG
            odc_flavour%ncol_end = NHDR_REAL_END
            odc_flavour%columnActive   => columnActive_RH
            odc_flavour%columnNameList => ocn_ColumnNameList_RH
            odc_flavour%activeIndexFromColumnIndex=>activeIndexFromColumnIndex_RH
            odc_flavour%activeIndexFromColumnIndex_defined = .false.
            odc_flavour%columnIndexFromActiveIndex=>columnIndexFromActiveIndex_RH
            odc_flavour%columnIndexFromActiveIndex_defined = .false.
         case ('BODY')
            odc_flavour%ncol_beg = NBDY_REAL_BEG
            odc_flavour%ncol_end = NBDY_REAL_END
            odc_flavour%columnActive   => columnActive_RB
            odc_flavour%columnNameList => ocn_ColumnNameList_RB
            odc_flavour%activeIndexFromColumnIndex=>activeIndexFromColumnIndex_RB
            odc_flavour%activeIndexFromColumnIndex_defined = .false.
            odc_flavour%columnIndexFromActiveIndex=>columnIndexFromActiveIndex_RB
            odc_flavour%columnIndexFromActiveIndex_defined = .false.
         end select

      case ('INT')
         select case (trim(headOrBody_in))
         case ('HEAD')
            odc_flavour%ncol_beg = NHDR_INT_BEG
            odc_flavour%ncol_end = NHDR_INT_END
            odc_flavour%columnActive   => columnActive_IH
            odc_flavour%columnNameList => ocn_ColumnNameList_IH
            odc_flavour%activeIndexFromColumnIndex=>activeIndexFromColumnIndex_IH
            odc_flavour%activeIndexFromColumnIndex_defined = .false.
            odc_flavour%columnIndexFromActiveIndex=>columnIndexFromActiveIndex_IH
            odc_flavour%columnIndexFromActiveIndex_defined = .false. 
         case ('BODY') 
            odc_flavour%ncol_beg = NBDY_INT_BEG
            odc_flavour%ncol_end = NBDY_INT_END
            odc_flavour%columnActive   => columnActive_IB
            odc_flavour%columnNameList => ocn_ColumnNameList_IB
            odc_flavour%activeIndexFromColumnIndex=>activeIndexFromColumnIndex_IB
            odc_flavour%activeIndexFromColumnIndex_defined = .false.
            odc_flavour%columnIndexFromActiveIndex=>columnIndexFromActiveIndex_IB
            odc_flavour%columnIndexFromActiveIndex_defined = .false.
         end select
      end select

   end subroutine odc_initColumnFlavour


   subroutine odc_class_initialize(obsColumnMode)
      !
      ! :Purpose: Set variables that use the same values for all instances of the class.
      !
      implicit none
      ! mode controlling the subset of columns that are activated in all objects
      character(len=*), intent(in) :: obsColumnMode

      integer :: column_index, list_index, ii
      integer, parameter :: COLUMN_LIST_SIZE = 100
      integer, dimension(COLUMN_LIST_SIZE) :: hdr_int_column_list, &
                  hdr_real_column_list, bdy_int_column_list, bdy_real_column_list

      ! Initialize the four column flavours:
      call odc_initColumnFlavour(odc_flavour_IB, 'INT',  'BODY')
      call odc_initColumnFlavour(odc_flavour_IH, 'INT',  'HEAD')
      call odc_initColumnFlavour(odc_flavour_RB, 'REAL', 'BODY')
      call odc_initColumnFlavour(odc_flavour_RH, 'REAL', 'HEAD')

      COLUMN_MODE:if(trim(obsColumnMode) == 'ALL') then
         do column_index=NHDR_INT_BEG,NHDR_INT_END
            call odc_activateColumn(odc_flavour_IH, column_index)
         enddo
         do column_index=NHDR_REAL_BEG,NHDR_REAL_END
            call odc_activateColumn(odc_flavour_RH, column_index)
         enddo
         do column_index=NBDY_INT_BEG,NBDY_INT_END
            call odc_activateColumn(odc_flavour_IB, column_index)
         enddo
         do column_index=NBDY_REAL_BEG,NBDY_REAL_END
            call odc_activateColumn(odc_flavour_RB, column_index)
         enddo

      elseif(trim(obsColumnMode) == 'ENKF') then COLUMN_MODE

         hdr_int_column_list= &
            (/OBS_RLN, OBS_ONM, OBS_INS, OBS_OTP, OBS_ITY, OBS_SAT, OBS_TEC, &
              OBS_DAT, OBS_ETM, OBS_NLV, OBS_STYP, OBS_PAS, OBS_REG, OBS_IP,  &
              OBS_ST1, OBS_IDF, OBS_SEN, &
              OBS_GQF, OBS_GQL, OBS_ROQF, (0,ii=21,100) /)

         hdr_real_column_list= &
            (/OBS_LAT, OBS_LON, OBS_ALT, OBS_BX,  OBS_BY,  OBS_BZ, OBS_TRAD, &
              OBS_GEOI,OBS_CLF,OBS_SUN,OBS_SZA,OBS_AZA,OBS_SAZ, &
              OBS_RZAM, OBS_RELE, OBS_RANS, OBS_RANE, (0,ii=18,100)/)

         bdy_int_column_list(:)    = 0
         bdy_int_column_list(1:size(odc_ENKF_bdy_int_column_list)) = &
            odc_ENKF_bdy_int_column_list(:)

         bdy_real_column_list(:)   = 0
         bdy_real_column_list(1:size(odc_ENKF_bdy_real_column_list)) = &
            odc_ENKF_bdy_real_column_list(:)

         do list_index=1,COLUMN_LIST_SIZE
            column_index = hdr_int_column_list(list_index)
            if(column_index == 0) exit
            call odc_activateColumn(odc_flavour_IH, column_index)
         end do

         do list_index=1,COLUMN_LIST_SIZE
            column_index = hdr_real_column_list(list_index)
            if(column_index == 0) exit
            call odc_activateColumn(odc_flavour_RH, column_index)
         end do

         do list_index=1,COLUMN_LIST_SIZE
            column_index = bdy_int_column_list(list_index)
            if(column_index == 0) exit
            call odc_activateColumn(odc_flavour_IB, column_index)
         end do

         do list_index=1,COLUMN_LIST_SIZE
            column_index = bdy_real_column_list(list_index)
            if(column_index == 0) exit
            call odc_activateColumn(odc_flavour_RB, column_index)
         end do

      elseif(trim(obsColumnMode) == 'ENKFMIDAS') then COLUMN_MODE

         hdr_int_column_list= &
            (/OBS_RLN, OBS_ONM, OBS_INS, OBS_OTP, OBS_ITY, OBS_SAT, OBS_TEC, &
              OBS_DAT, OBS_ETM, OBS_NLV, OBS_STYP,OBS_PAS, OBS_REG, OBS_IP,  &
              OBS_ST1, OBS_IDF, OBS_SEN, &
              OBS_SWQ1,OBS_SWQ2,OBS_SWMT,OBS_SWLS,OBS_SWGA,OBS_SWHA, &
              OBS_GQF, OBS_GQL, OBS_ROQF, OBS_TTYP, (0,ii=28,100) /)

         hdr_real_column_list= &
            (/OBS_LAT, OBS_LON, OBS_ALT, OBS_BX,  OBS_BY,  OBS_BZ, OBS_TRAD, &
              OBS_GEOI,OBS_CLF, OBS_SUN, OBS_SZA, OBS_AZA, OBS_SAZ, &
              OBS_RZAM, OBS_RELE, OBS_RANS, OBS_RANE, (0,ii=18,100)/)

         bdy_int_column_list= &
            (/OBS_VNM, OBS_FLG, OBS_ASS, OBS_HIND,OBS_VCO, OBS_LYR, &
              OBS_XTR, OBS_QCF2, OBS_CLA, (0,ii=10,100) /)

         bdy_real_column_list= &
            (/OBS_PPP, OBS_SEM, OBS_VAR, OBS_OMP, OBS_OMA, OBS_OER, OBS_HPHT,&
              OBS_HAHT,OBS_ZHA, OBS_OMP6, OBS_OMA0, OBS_SIGI, OBS_SIGO, OBS_WORK, OBS_PRM,&
              OBS_JOBS, OBS_BCOR,OBS_LOND,OBS_LATD, (0,ii=20,100) /)

         do list_index=1,COLUMN_LIST_SIZE
            column_index = hdr_int_column_list(list_index)
            if(column_index == 0) exit
            call odc_activateColumn(odc_flavour_IH, column_index)
         end do

         do list_index=1,COLUMN_LIST_SIZE
            column_index = hdr_real_column_list(list_index)
            if(column_index == 0) exit
            call odc_activateColumn(odc_flavour_RH, column_index)
         end do

         do list_index=1,COLUMN_LIST_SIZE
            column_index = bdy_int_column_list(list_index)
            if(column_index == 0) exit
            call odc_activateColumn(odc_flavour_IB, column_index)
         end do

         do list_index=1,COLUMN_LIST_SIZE
            column_index = bdy_real_column_list(list_index)
            if(column_index == 0) exit
            call odc_activateColumn(odc_flavour_RB, column_index)
         end do

      elseif(trim(obsColumnMode) == 'VAR') then COLUMN_MODE

         do column_index=NHDR_INT_BEG,NHDR_INT_END
            if( column_index < OBS_GQF .or. column_index > OBS_NCO2  &
              ) call odc_activateColumn(odc_flavour_IH, column_index)
         enddo
         do column_index=NHDR_REAL_BEG,NHDR_REAL_END
            if(     column_index /= OBS_BX &
               .and.column_index /= OBS_BY &
               .and.column_index /= OBS_BZ &
               .and. (column_index < OBS_M1C1 .or. &
                      column_index > OBS_ZPS) &
              ) call odc_activateColumn(odc_flavour_RH, column_index)
         enddo

         do column_index=NBDY_INT_BEG,NBDY_INT_END
            if( column_index /= OBS_KFA ) call odc_activateColumn(odc_flavour_IB, column_index)
         enddo
         do column_index = NBDY_REAL_BEG, NBDY_REAL_END
            if(      column_index /= OBS_OMP6 &
               .and. column_index /= OBS_OMA0 &
               .and. column_index /= OBS_HAHT &
               .and. column_index /= OBS_SIGI &
               .and. column_index /= OBS_SIGO &
              )call odc_activateColumn(odc_flavour_RB, column_index)
         enddo

      endif COLUMN_MODE
   end subroutine odc_class_initialize


   subroutine odc_columnElem(odc_array, column_index, row_index, value_i,value_r)
     !
     ! :Purpose: Returns the value of the row_index'th element in the column array
     !           with the indicated column_index.
     !
     !           The column array can be of any one of the four possible column-array
     !           flavours.  The flavour is selected by one of four wrappers to this
     !           method.
     !
     implicit none
     type(struct_obsDataColumn_Array), intent(in)  :: odc_array
     integer                         , intent(in)  :: column_index
     integer                         , intent(in)  :: row_index
     integer                         , intent(out) :: value_i
     real(pre_obsReal)               , intent(out) :: value_r

     character(len=100) :: message
      
      if(      column_index >= odc_array%odc_flavour%ncol_beg &
         .and. column_index <= odc_array%odc_flavour%ncol_end) then
         if(odc_array%odc_flavour%columnActive(column_index)) then
            ! Return the value (return int AND real, just to make it simple)
            value_i = odc_array%columns(column_index)%value_i(row_index)
            value_r = odc_array%columns(column_index)%value_r(row_index)
         else
            write(message,*)'abort in odc_columnElem [' &
                          // odc_array%odc_flavour%dataType //',' &
                          // odc_array%odc_flavour%headOrBody // &
                          ']: column not active: ', &
                          odc_array%odc_flavour%columnNameList(column_index)
            call odc_abort(message)
         endif
      else
         write(message,*) 'abort in odc_columnElem [' &
                          // odc_array%odc_flavour%dataType //',' &
                          // odc_array%odc_flavour%headOrBody // &
                          ']: column index out of range: ', column_index
         call odc_abort(message); return
      endif
   end subroutine odc_columnElem


   function odc_columnIndexFromActiveIndex(odc_flavour,active_index_in, &
                                           recompute) result(column_index_out)
      !
      ! :Purpose: The list of active columns is only a subset of all possible
      !           columns.  Return the index into the list of all columns, given
      !           the index into the list of active columns, and given the column
      !           flavour.
      !
      implicit none
      type(struct_odc_flavour), intent(inout) :: odc_flavour
      integer                 , intent(in)    :: active_index_in
      logical, optional       , intent(in)    :: recompute
      integer                                 :: column_index_out

      integer :: active_index, &
                 column_index

      if(present(recompute)) then
         if(recompute) odc_flavour%columnIndexFromActiveIndex_defined=.false.
      endif

      if(.not. odc_flavour%columnIndexFromActiveIndex_defined) then
         odc_flavour%columnIndexFromActiveIndex_defined=.true.
         active_index=0

         do column_index = odc_flavour%ncol_beg, odc_flavour%ncol_end
            if(odc_flavour%columnActive(column_index)) then
               active_index=active_index+1
               odc_flavour%columnIndexFromActiveIndex(active_index) =column_index
            endif
         enddo
      endif

      column_index_out=odc_flavour%columnIndexFromActiveIndex(active_index_in)
   end function odc_columnIndexFromActiveIndex


   subroutine odc_columnSet(odc_array, column_index, row_index, &
                            value_i, value_r, numElements, numElements_max)
      !
      ! :Purpose: Sets the value of the row_index'th element in the column array with
      !           the indicated column_index.
      !
      !           The column array can be of any one of the four possible column-array
      !           flavours.  The flavour is selected by one of four wrappers to this
      !           method.
      !
      implicit none
      type(struct_obsDataColumn_Array), intent(inout) :: odc_array
      integer                         , intent(in)    :: column_index
      integer                         , intent(in)    :: row_index
      integer                         , intent(in)    :: value_i
      real(pre_obsReal)               , intent(in)    :: value_r
      integer                         , intent(inout) :: numElements
      integer                         , intent(in)    :: numElements_max

      character(len=100) :: message

      ! Validate the requested row_index, and
      ! Increment the number of elements, if necessary
      if(row_index > numElements_max) then
         write(message,*)'The requested ', &
                         trim(odc_array%odc_flavour%headOrBody), ' row_index, ',&
                         row_index,', is greater than the maximum, ', &
                         numElements_max
         call odc_abort(message)

      else if(row_index > numElements+1) then
         write(message,*)'The requested ', &
                         trim(odc_array%odc_flavour%headOrBody), &
                         ' row_index, ', row_index, &
                         ', is beyond the next available index, ', &
                         numElements+1
         call odc_abort(message)

      else if(row_index == numElements+1) then
         numElements = numElements+1
      end if


      ! Validate the requested column index, and
      ! Record the value
      if(      column_index >= odc_array%odc_flavour%ncol_beg &
         .and. column_index <= odc_array%odc_flavour%ncol_end) then
         if(odc_array%odc_flavour%columnActive(column_index)) then
            ! Record the value (record int AND real, just to make it simple)
            odc_array%columns(column_index)%value_i(row_index) = value_i
            odc_array%columns(column_index)%value_r(row_index) = value_r
         else
            write(message,*) 'abort in odc_columnSet [' &
                          // odc_array%odc_flavour%dataType //',' &
                          // odc_array%odc_flavour%headOrBody // &
                          ']: column not active: ', &
                          odc_array%odc_flavour%columnNameList(column_index)
            call odc_abort(message)
         endif
      else
         write(message,*) 'abort in odc_columnSet [' &
                          // odc_array%odc_flavour%dataType //',' &
                          // odc_array%odc_flavour%headOrBody // &
                          ']: column index out of range: ', column_index
         call odc_abort(message); return
      endif

   end subroutine odc_columnSet


   subroutine odc_deallocate(odc)
      !
      ! :Purpose: Deallocate a single column of obs data
      !
      ! :Arguments:
      !           :odc: instance of the obsDataColumn type
      !
      implicit none
      type(struct_obsDataColumn), intent(inout) :: odc

      if(.not.odc%allocated) then
         call odc_abort('ODC_DEALLOCATE: column is not already allocated.')
      endif

      odc%allocated=.false.
      if(trim(odc%dataType) == 'INT' .and. associated(odc%value_i)) then
         deallocate(odc%value_i)
         nullify(odc%value_i)
         nullify(odc%value_r)    ! Dont deallocate:  this is not the only pointer
      end if

      if(trim(odc%dataType) == 'REAL' .and. associated(odc%value_r)) then
         deallocate(odc%value_r)
         nullify(odc%value_r)
         nullify(odc%value_i)    ! Dont deallocate:  this is not the only pointer
      endif

   end subroutine odc_deallocate


   function odc_numActiveColumn(odc_array) result(numActiveColumn)
      !
      ! :Purpose: Return the number of active columns that are contained in the given
      !           column array.
      !
      !           The column array can be of any one of the four possible column-array
      !           flavours.
      !
      implicit none
      type(struct_obsDataColumn_Array), intent(in) :: odc_array
      integer :: numActiveColumn

      integer :: column_index

      numActiveColumn=0
      do column_index = odc_array%odc_flavour%ncol_beg, &
                        odc_array%odc_flavour%ncol_end
         if(odc_array%odc_flavour%columnActive(column_index))  &
            numActiveColumn=numActiveColumn+1
      enddo
   end function odc_numActiveColumn

end module ObsDataColumn_mod


module ObsSpaceData_mod
  ! MODULE obsSpaceData_mod (prefix='obs' category='2. High-level data objects')
  !
  ! :Purpose: This module contains a structure definition, struct_obs, and
  !           methods for accessing its data.  An instance of struct_obs contains
  !           all information pertaining to a set of observation-space data.
  !           (This had evolved from the CMA structure, originated in work by
  !           D. Vasiljevic at ECMWF.)
  !
   use codePrecision_mod
   use ObsColumnNames_mod
   use ObsDataColumn_mod
   use IndexListDepot_mod
   use mathPhysConstants_mod
   implicit none
   save
   private



   ! CLASS-CONSTANT:
   ! CLASS-CONSTANT:
   ! CLASS-CONSTANT:

   logical, save :: obs_class_initialized = .false.

   ! end of CLASS-CONSTANT variables.
   ! end of CLASS-CONSTANT variables
   ! end of CLASS-CONSTANT variables.


   ! PUBLIC METHODS:
   public obs_append     ! append an obsdat object to another obsdat object
   public obs_bodyElem_i ! obtain an integer body element from observation object
   public obs_bodyElem_r ! obtain a real body element from the observation object
   public obs_bodyPrimaryKey ! obtain the body primary key value
   public obs_bodyIndex_mpiglobal ! obtain mpiglobal body row index
   public obs_bodySet_i  ! set an integer body value in the observation object
   public obs_bodySet_r  ! set a real body value in the observation object
   public obs_class_initialize ! initialize class variables: column mode
   public obs_clean      ! remove from obs data those that not to be assimilated
   public obs_clean2     ! modified version of obs_clean used by MIDAS
   public obs_columnActive_IB ! return the active status for a column (T/F)
   public obs_columnActive_IH !                 "
   public obs_columnActive_RB !                 "
   public obs_columnActive_RH !                 "
   public obs_isColumnNameValid      ! Check if column name is a valid obsSpaceData name
   public obs_columnIndexFromName    ! get the index from the name (any flavour)
   public obs_columnIndexFromName_IB ! get the index from the name
   public obs_columnIndexFromName_IH !         "
   public obs_columnIndexFromName_RB !         "
   public obs_columnIndexFromName_RH !         "
   public obs_columnDataType ! tell user if column index has real or integer data
   public obs_copy       ! copy an obsdat object
   public obs_count_headers ! count the stations and observations in the object
   public obs_elem_c     ! obtain character element from the observation object
   public obs_enkf_prntbdy! print all data records associated with an observation
   public obs_enkf_prnthdr! print the header of an observation record
   public obs_expandToMpiGlobal ! restore data for the mpi-global context
   public obs_extractObsRealBodyColumn    ! return entire selected column (real/body)
   public obs_extractObsRealBodyColumn_r4 ! return entire selected column (real/body), real4 version
   public obs_extractObsIntBodyColumn     ! return entire selected column (int/body)
   public obs_extractObsRealHeaderColumn  ! return entire selected column (real/header)
   public obs_extractObsIntHeaderColumn   ! return entire selected column (int/header)
   public obs_finalize   ! object clean-up
   public obs_generate_header ! fill in observation-data header, from burp files
   public obs_getBodyIndex ! obtain an element from the current body list
   public obs_getFamily  ! return the family of a datum
   public obs_getHeaderIndex ! obtain an element from the current header list
   public obs_getNchanAvhrr  ! to get the number of AVHRR channels
   public obs_getNclassAvhrr ! to get the number of AVHRR radiance classes
   public obs_headElem_i ! obtain an integer header element from the obs'n object
   public obs_headElem_r ! obtain real header element from the observation object
   public obs_headPrimaryKey ! obtain the header primary key value
   public obs_headerIndex_mpiglobal ! obtain mpiglobal header row index
   public obs_headSet_i  ! set an integer header value in the observation object
   public obs_headSet_r  ! set a real header value in the observation object
   public obs_initialize ! variable initialization
   public obs_mpiLocal   ! obtain the current mpi state of the observation object
   public obs_MpiRedistribute ! do a general redistribution of obs over mpi tasks
   public obs_numBody    ! returns the number of bodies recorded
   public obs_numBody_max! returns the dimensioned number of bodies
   public obs_numBody_mpiglobal ! returns mpi-global number of bodies recorded
   public obs_numHeader  ! returns the number of headers recorded
   public obs_numHeader_max ! returns the dimensioned number of headers
   public obs_numHeader_mpiglobal ! returns mpi-global number of headers recorded
   public obs_order      ! put obs data in the order required for assimilation
   public obs_print      ! obs_enkf_prnthdr & obs_enkf_prntbdy for each station
   public obs_prnt_csv   ! call obs_tosqlhdr and obs_tosqlbdy for each station
   public obs_prntbdy    ! print the body data for one header
   public obs_prnthdr    ! print the data contained in one header
   public obs_read       ! read the observation data from binary files
   public obs_creatSubCMA ! create a sub-CMA from the global CMA.  
   public obs_reduceToMpiLocal ! retain only data pertinent to the mpi-local PE
   public obs_squeeze    ! reallocate objects arrays to reduce memory use
   public obs_select     ! select observations in a vertical range
   public obs_setBodyPrimaryKey ! set the value of the body primary key
   public obs_setHeadPrimaryKey ! set the value of the body primary key
   public obs_set_c      ! set a character value in the observation object
   public obs_set_current_body_list   ! set a body list for a family as current
   public obs_set_current_header_list ! set a header list for a family as current
   public obs_setFamily  ! set the family of a datum
   public obs_sethind    ! set the header index in body table
   public obs_famExist   ! checks if a family is present in the observation set
   public obs_status     ! returns the values of the object's status variables
   public obs_write      ! write the observation data to binary files
                         ! (calls obs_write_hdr, obs_write_bdy, obs_write_hx
                         !  for each station)
   public obs_write_hx   ! write to binary files a station's interpolated values

   interface obs_getBodyIndex
      module procedure obs_getBodyIndex_depot
      module procedure obs_getBodyIndex_private
   end interface obs_getBodyIndex

   interface obs_set_current_body_list
      module procedure obs_set_current_body_list_from_family
      module procedure obs_set_current_body_list_from_header
      module procedure obs_set_current_body_list_all
   end interface obs_set_current_body_list

   interface obs_set_current_header_list
      module procedure obs_set_current_header_list_from_family
      module procedure obs_set_current_header_list_all
   end interface

   interface obs_bodySet_r
      module procedure obs_bodySet_r4
      module procedure obs_bodySet_r8
   end interface obs_bodySet_r

   interface obs_headSet_r
      module procedure obs_headSet_r4
      module procedure obs_headSet_r8
   end interface obs_headSet_r

   ! PRIVATE METHODS:
   private obs_abort     ! abort a job on error
   private obs_allocate  ! array allocation
   private obs_columnIndexFromNameForFlavour ! get the index from the name
   private obs_deallocate! array de-allocation
   private obs_mpiDistributeIndices ! distribute header & body indices for mpi parallelization
   private obs_tosqlbdy  ! write the observation data in comma-separated format
   private obs_tosqlhdr  ! write the observation header in comma-separated format
   private obs_write_bdy ! write the observation data to binary files
   private obs_write_hdr ! write the observation header to binary files


   ! PARAMETERS INHERITED FROM ObsColumnNames_mod (make them public)
   !    integer-header column numbers
   public :: OBS_RLN, OBS_ONM, OBS_INS, OBS_OTP, OBS_ITY, OBS_SAT, OBS_TEC
   public :: OBS_DAT, OBS_ETM, OBS_NLV, OBS_PAS, OBS_REG, OBS_IP , OBS_SEN
   public :: OBS_IPF, OBS_IPC, OBS_IPT
   public :: OBS_ST1, OBS_IDF
   public :: OBS_GQF, OBS_GQL
   public :: OBS_NCO2,OBS_STYP,OBS_ROQF
   public :: OBS_SWQ1,OBS_SWQ2,OBS_SWMT,OBS_SWLS,OBS_SWGA,OBS_SWHA
   public :: OBS_CHM, OBS_FOV, OBS_PRFL, OBS_PHAS, OBS_ORI
   public :: OBS_LCH, OBS_RTP, OBS_HDD, OBS_HDT, OBS_TFLG, OBS_LFLG
   public :: OBS_ORBI,OBS_AQF1, OBS_AQF2, OBS_AQF3, OBS_TTYP, OBS_INFG
   public :: OBS_RAIN, OBS_CHID

   !    real-header column numbers
   public :: OBS_LAT, OBS_LON, OBS_ALT, OBS_BX,  OBS_BY,  OBS_BZ
   public :: OBS_M1C1, OBS_M1C2, OBS_M1C3, OBS_M1C4, OBS_M1C5, OBS_M1C6
   public :: OBS_M2C1, OBS_M2C2, OBS_M2C3, OBS_M2C4, OBS_M2C5, OBS_M2C6
   public :: OBS_M3C1, OBS_M3C2, OBS_M3C3, OBS_M3C4, OBS_M3C5, OBS_M3C6
   public :: OBS_M4C1, OBS_M4C2, OBS_M4C3, OBS_M4C4, OBS_M4C5, OBS_M4C6
   public :: OBS_M5C1, OBS_M5C2, OBS_M5C3, OBS_M5C4, OBS_M5C5, OBS_M5C6
   public :: OBS_M6C1, OBS_M6C2, OBS_M6C3, OBS_M6C4, OBS_M6C5, OBS_M6C6
   public :: OBS_M7C1, OBS_M7C2, OBS_M7C3, OBS_M7C4, OBS_M7C5, OBS_M7C6
   public :: OBS_S1C1, OBS_S1C2, OBS_S1C3, OBS_S1C4, OBS_S1C5, OBS_S1C6
   public :: OBS_S2C1, OBS_S2C2, OBS_S2C3, OBS_S2C4, OBS_S2C5, OBS_S2C6
   public :: OBS_S3C1, OBS_S3C2, OBS_S3C3, OBS_S3C4, OBS_S3C5, OBS_S3C6
   public :: OBS_S4C1, OBS_S4C2, OBS_S4C3, OBS_S4C4, OBS_S4C5, OBS_S4C6
   public :: OBS_S5C1, OBS_S5C2, OBS_S5C3, OBS_S5C4, OBS_S5C5, OBS_S5C6
   public :: OBS_S6C1, OBS_S6C2, OBS_S6C3, OBS_S6C4, OBS_S6C5, OBS_S6C6
   public :: OBS_S7C1, OBS_S7C2, OBS_S7C3, OBS_S7C4, OBS_S7C5, OBS_S7C6
   public :: OBS_CF1,  OBS_CF2,  OBS_CF3,  OBS_CF4,  OBS_CF5,  OBS_CF6, OBS_CF7
   public :: OBS_ETOP, OBS_VTOP, OBS_ECF,  OBS_VCF , OBS_HE  , OBS_ZTSR
   public :: OBS_ZTM , OBS_ZTGM, OBS_ZLQM, OBS_ZPS , OBS_TRAD, OBS_GEOI
   public :: OBS_CLF , OBS_SUN,  OBS_SZA,  OBS_AZA , OBS_SAZ , OBS_CLWO, OBS_CLWB, OBS_MWS
   public :: OBS_SCAT, OBS_IWV,  OBS_RZAM, OBS_RELE, OBS_RANS, OBS_RANE
   !    integer-body column numbers
   public :: OBS_VNM, OBS_FLG, OBS_KFA, OBS_ASS, OBS_HIND,OBS_VCO, OBS_LYR
   public :: OBS_XTR, OBS_QCF2, OBS_CLA

   !    real-body column numbers
   public :: OBS_PPP, OBS_SEM, OBS_VAR, OBS_OMP, OBS_OMA, OBS_OER, OBS_HPHT
   public :: OBS_HAHT,OBS_ZHA, OBS_OMP6,OBS_OMA0,OBS_SIGI,OBS_SIGO,OBS_POB
   public :: OBS_WORK,OBS_PRM, OBS_JOBS,OBS_QCV, OBS_FSO, OBS_CRPS,OBS_BCOR
   public :: OBS_OMPE,OBS_LATD,OBS_LOND,OBS_BTCL,OBS_LOCI

   ! OBSERVATION-SPACE FUNDAMENTAL PARAMETERS
   integer, public, parameter :: obs_assimilated    = 1 ! OBS_ASS value for assimilated obs
   integer, public, parameter :: obs_notAssimilated = 0 ! OBS_ASS value for non assimilated obs

   real(pre_obsReal), public, parameter :: obs_missingValue_R = real(MPC_missingValue_R8, pre_obsReal) ! Missing value

   ! DERIVED TYPE AND MODULE VARIABLE DECLARATIONS

                                        ! It is intended that these null values
                                        ! be used with scratchRealHeader, etc.
   real(pre_obsReal), parameter :: NULL_COLUMN_VALUE_R = real(-9.99D9, pre_obsReal)
   integer          , parameter :: NULL_COLUMN_VALUE_I = -9.99

   ! This type is the goal of the ObsSpaceData and supporting modules.  An
   ! instance of this derived type contains all information pertaining to a set
   ! of observation-space data.
   type, public :: struct_obs
      private
      logical          :: allocated = .false.
                                        ! Two internal structures for
                                        ! facilitating multiple traversals of
                                        ! subsets of the observation-space data.
      type(struct_index_list_depot) :: header_index_list_depot
      type(struct_index_list_depot) :: body_index_list_depot
                                        ! For the cstnid and cfamily arrays:
                                        !   1st dim'n:  row index
      character(len=12),  pointer :: cstnid(:)
      character(len=2),   pointer :: cfamily(:)
                                        ! The primary keys for both tables
      integer(8), pointer :: headerPrimaryKey(:), bodyPrimaryKey(:)
                                        ! The four arrays of data columns
      type(struct_obsDataColumn_Array) :: &
                             realHeaders, & ! real header columns
                             intHeaders,  & ! integer header columns
                             realBodies,  & ! real body columns
                             intBodies      ! integer body columns

      ! These scratch arrays are the basis for allowing the manipulation of both
      ! (real and integer) values of any obsDataColumn, without having to decide
      ! which one (real or integer) contains the sought value.  This ability
      ! simplifies the code.
      real(pre_obsReal), pointer :: scratchRealHeader(:), &
                                    scratchRealBody  (:)
      integer          , pointer :: scratchIntHeader (:), &
                                    scratchIntBody   (:)

      integer :: numHeader              ! Actual number of headers on record
      integer :: numHeader_max          ! maximum number of headers(i.e.stations)
      integer :: numBody                ! Actual total number of bodies on record
      integer :: numBody_max            ! maximum number of bodies (i.e. data)

                                        ! row indices of mpiglobal data, useful
                                        ! for transforming from mpiLocal back to
                                        ! mpiGlobal
                                        !    1st dim'n: row index, only in
                                        !               mpiLocal context
      integer, pointer, dimension(:) :: headerIndex_mpiglobal  => NULL()
      integer, pointer, dimension(:) :: bodyIndex_mpiglobal    => NULL()

      logical :: mpi_local       ! T: keep only data needed by this PE (mpilocal)
   end type struct_obs


contains


   subroutine obs_abort(cdmessage)
      !
      ! :Purpose: To stop a job when an error occurred
      !
      ! :Arguments:
      !           :cdmessage: message to be printed
      !

!#if defined(UNIT_TESTING)
!      use pFUnit
!#endif
      implicit none
      character(len=*), intent(in) :: cdmessage

      write(*,'(//,4X,"ABORTING IN ObsSpaceData_mod:-------",/,8X,A)')cdmessage
      call flush(6)

!#if defined(UNIT_TESTING)
!      call throw(Exception('exiting in obs_abort'))
!#else
      call qqexit(13)
      stop
!#endif
   end subroutine obs_abort


   subroutine obs_allocate(obsdat, numHeader_max, numBody_max, silent)
      !
      ! :Purpose: Allocate arrays according to the parameters, numHeader_max and
      !           numBody_max.  This is a private method.
      !
      implicit none

      type(struct_obs), intent(inout) :: obsdat
      integer,          intent(in)    :: numHeader_max
      integer,          intent(in)    :: numBody_max
      logical, optional,intent(in)    :: silent

      logical :: silent_
      integer :: column_index

      if(obsdat%allocated) then
         call obs_abort('OBS_ALLOCATE:  a second allocation of ObsSpaceData has been attempted.')
         return
      endif
      obsdat%allocated=.true.

      if(present(silent)) then
         silent_  = silent
      else
         silent_  = .false.
      end if

      if(.not. silent_) then
         write(*,*) ' DIMENSIONS OF OBSERVATION ARRAYS:'
         write(*,*) ' numHeader_max = ',numHeader_max,'  numBody_max = ', &
                    numBody_max
      end if
      obsdat%numHeader_max=numHeader_max
      obsdat%numBody_max=numBody_max

      !
      ! ALLOCATE THE ARRAYS and initialize the contents
      !
      allocate(obsdat%realHeaders%columns(NHDR_REAL_BEG:NHDR_REAL_END))
      obsdat%realHeaders%odc_flavour => odc_flavour_RH

      allocate(obsdat%intHeaders%columns(NHDR_INT_BEG:NHDR_INT_END))
      obsdat%intHeaders%odc_flavour => odc_flavour_IH

      HEADER:if(numHeader_max > 0) then
         allocate(obsdat%headerPrimaryKey(numHeader_max))
         obsdat%headerPrimaryKey(:) = 0

         allocate(obsdat%cfamily(numHeader_max))
         obsdat%cfamily(:)='XX'

         allocate(obsdat%cstnid(numHeader_max))
         obsdat%cstnid(:)='XXXXXXXXXXXX'

         allocate(obsdat%scratchRealHeader(numHeader_max))
         allocate(obsdat%scratchIntHeader (numHeader_max))
         obsdat%scratchRealHeader(:) = NULL_COLUMN_VALUE_R
         obsdat%scratchIntHeader (:) = NULL_COLUMN_VALUE_I

         do column_index=NHDR_REAL_BEG,NHDR_REAL_END
            if(obsdat%realHeaders%odc_flavour%columnActive(column_index))  &
               call odc_allocate(obsdat%realHeaders%columns(column_index), &
                                 numHeader_max, &
                                 ocn_ColumnNameList_RH(column_index),'REAL',&
                                 obsdat%scratchRealHeader, &
                                 obsdat%scratchIntHeader)
         enddo

         do column_index=NHDR_INT_BEG,NHDR_INT_END
            if(obsdat%intHeaders%odc_flavour%columnActive(column_index))  &
               call odc_allocate(obsdat%intHeaders%columns(column_index), &
                                 numHeader_max, &
                                 ocn_ColumnNameList_IH(column_index),'INT ', &
                                 obsdat%scratchRealHeader, &
                                 obsdat%scratchIntHeader)
         enddo
      endif HEADER

      allocate(obsdat%realBodies%columns(NBDY_REAL_BEG:NBDY_REAL_END))
      obsdat%realBodies%odc_flavour => odc_flavour_RB

      allocate(obsdat%intBodies%columns(NBDY_INT_BEG:NBDY_INT_END))
      obsdat%intBodies%odc_flavour => odc_flavour_IB

      BODY:if(numBody_max > 0) then
         allocate(obsdat%bodyPrimaryKey(numBody_max))
         obsdat%bodyPrimaryKey(:)   = 0

         allocate(obsdat%scratchRealBody(numBody_max))
         allocate(obsdat%scratchIntBody (numBody_max))
         obsdat%scratchRealBody(:) = NULL_COLUMN_VALUE_R
         obsdat%scratchIntBody (:) = NULL_COLUMN_VALUE_I

         do column_index=NBDY_REAL_BEG,NBDY_REAL_END
            if(obsdat%realBodies%odc_flavour%columnActive(column_index))  &
               call odc_allocate(obsdat%realBodies%columns(column_index), &
                                 numBody_max, &
                                 ocn_ColumnNameList_RB(column_index),'REAL', &
                                 obsdat%scratchRealBody, obsdat%scratchIntBody)
         enddo

         do column_index=NBDY_INT_BEG,NBDY_INT_END
            if(obsdat%intBodies%odc_flavour%columnActive(column_index))  &
               call odc_allocate(obsdat%intBodies%columns(column_index), &
                                 numBody_max, &
                                 ocn_ColumnNameList_IB(column_index),'INT ', &
                                 obsdat%scratchRealBody, obsdat%scratchIntBody)
         enddo
      endif BODY

      call ild_initialize(obsdat%header_index_list_depot, obsdat%numHeader_max)
      call ild_initialize(obsdat%body_index_list_depot,   obsdat%numBody_max)
   end subroutine obs_allocate


   subroutine obs_append( obsdat, hx, obs_out, hx_out )
      !
      ! :Purpose: with a call of type obs_append(obs_1,obs_2) append obs_1 to obs_2
      !

      type (struct_obs), intent(in)    :: obsdat
      type (struct_obs), intent(inout) :: obs_out
      real(8),           intent(in)    :: hx(:,:)
      real(8),           intent(inout) :: hx_out(:,:)

      integer :: i_data_read,i_data_read_first,i_data_read_last,i_data_write
      integer :: i_last,istation,i_station_write,nens,i_write_first
      integer :: loc,loc_last
      integer :: column_index, pass_offset

      if (obsdat%numHeader == 0) then
         write(*,*) 'odd input for routine obs_append'
         write(*,*) 'no stations need to be added to the obsdat.'
         return
      endif

      nens=size(hx,1)

      ! Locate the first available locations in the output
      if (obs_out%numHeader >= 1) then
                                        ! Memorize pass_offset, since numHeader
                                        ! will change
         pass_offset=obs_headElem_i(obs_out,OBS_PAS,obs_out%numHeader) 

         i_last=1
         loc_last=obs_headElem_i(obs_out, OBS_RLN, 1)
         storedlast: do istation=2,obs_out%numHeader
            loc=obs_headElem_i(obs_out, OBS_RLN, istation)
            if (loc > loc_last) then
               i_last=istation
               loc_last=loc
            endif
         enddo storedlast
         ! The first available locations in the output:
         i_station_write=obs_out%numHeader+1
         i_data_write=obs_headElem_i(obs_out, OBS_RLN, i_last) &
                     +obs_headElem_i(obs_out, OBS_NLV, i_last)
      else
         pass_offset=0
         i_station_write=1
         i_data_write=1
      endif

      stations: do istation=1,obsdat%numHeader
         i_data_read_first=obs_headElem_i(obsdat, OBS_RLN, istation)
         i_data_read_last=i_data_read_first &
                         +obs_headElem_i(obsdat, OBS_NLV, istation) -1
         i_write_first=i_data_write

         observations: do i_data_read=i_data_read_first,i_data_read_last
            do column_index=NBDY_INT_BEG,NBDY_INT_END
               if(obsdat%intBodies%odc_flavour%columnActive(column_index)) &
                  call obs_bodySet_i(obs_out, column_index, i_data_write, &
                       obs_bodyElem_i(obsdat, column_index, i_data_read))
            enddo
            do column_index=NBDY_REAL_BEG,NBDY_REAL_END
               if(obsdat%realBodies%odc_flavour%columnActive(column_index)) &
                  call obs_bodySet_r(obs_out, column_index, i_data_write, &
                       obs_bodyElem_r(obsdat, column_index, i_data_read))
            enddo
            hx_out(1:nens,i_data_write)=hx(:,i_data_read)
                                        ! Make HIND point to new header row_index
            call obs_bodySet_i(obs_out, OBS_HIND, i_data_write, i_station_write)
            i_data_write=i_data_write+1
         enddo observations

         do column_index=NHDR_INT_BEG,NHDR_INT_END
            if(obsdat%intHeaders%odc_flavour%columnActive(column_index)) &
               call obs_headSet_i(obs_out, column_index, i_station_write, &
                    obs_headElem_i(obsdat, column_index, istation))
         enddo
         call obs_headSet_i(obs_out, OBS_ONM, i_station_write, i_station_write)
         call obs_headSet_i(obs_out, OBS_RLN, i_station_write, i_write_first)

         if (obs_out%numHeader > 0) then
            call obs_headSet_i(obs_out, OBS_PAS, i_station_write, &
                               obs_headElem_i(obs_out,OBS_PAS,i_station_write) &
                               +pass_offset&
                              )
         end if

         do column_index=NHDR_REAL_BEG,NHDR_REAL_END
            if(obsdat%realHeaders%odc_flavour%columnActive(column_index)) &
               call obs_headSet_r(obs_out, column_index, i_station_write, &
                    obs_headElem_r(obsdat, column_index, istation))
         enddo
         obs_out%headerPrimaryKey(i_station_write) = obsdat%headerPrimaryKey(istation)

         obs_out%cstnid(  i_station_write)=obsdat%cstnid ( istation)
         obs_out%cfamily( i_station_write)=obsdat%cfamily( istation)

         i_station_write=i_station_write+1
      enddo stations
   end subroutine obs_append


   function obs_bodyElem_i( obsdat, column_index, row_index) result(value_i)
      !
      ! :Purpose: To control access to the observation object.  Returns the (integer)
      !           value of the row_index'th ObsData element with the indicated column
      !           index from the "body".
      !
      implicit none
      integer                       :: value_i
      type(struct_obs), intent(in)  :: obsdat
      integer         , intent(in)  :: column_index
      integer         , intent(in)  :: row_index

      real(pre_obsReal) :: value_r         ! not used

      call odc_columnElem(obsdat%intBodies, column_index, row_index, &
                          value_i, value_r)
   end function obs_bodyElem_i


   function obs_bodyElem_r(obsdat,column_index,row_index) result(value_r)
      !
      ! :Purpose: Get a real-valued body observation-data element
      !      To control access to the observation object.  Returns the (real)
      !      value of the row_index'th ObsData element with the indicated column
      !      index from the "body".
      !
      implicit none
      real(pre_obsReal)             :: value_r
      type(struct_obs), intent(in)  :: obsdat
      integer         , intent(in)  :: column_index
      integer         , intent(in)  :: row_index

      integer :: value_i                ! not used

      call odc_columnElem(obsdat%realBodies, column_index, row_index, &
                          value_i, value_r)
   end function obs_bodyElem_r


   function obs_bodyPrimaryKey(obsdat,row_index) result(primaryKey)
      !
      ! :Purpose: Get the body primary key value.
      !
      implicit none
      integer(8)                    :: primaryKey
      type(struct_obs), intent(in)  :: obsdat
      integer         , intent(in)  :: row_index

      primaryKey = obsdat%bodyPrimaryKey(row_index)

   end function obs_bodyPrimaryKey


   function obs_bodyIndex_mpiglobal(obsdat,row_index) result(value)
      !
      ! :Purpose: Get the mpiglobal body row_index.
      !      To control access to the mpiglobal row_index into the "body".
      !
      implicit none
      integer value
      type(struct_obs), intent(in)  :: obsdat
      integer         , intent(in)  :: row_index

      value=obsdat%bodyIndex_mpiglobal(row_index)
   end function obs_bodyIndex_mpiglobal


   subroutine obs_bodySet_i(obsdat, column_index, row_index, value_i)
      !
      ! :Purpose: Set an integer-valued body observation-data element.
      !      To control access to the observation object.  Sets the (integer)
      !      value of the row_index'th ObsData element with the indicated column
      !      index from the "body".
      !
      implicit none
      type(struct_obs), intent(inout)  :: obsdat
      integer         , intent(in)     :: column_index
      integer         , intent(in)     :: row_index
      integer         , intent(in)     :: value_i

      call odc_columnSet(obsdat%intBodies, column_index, row_index, &
                         value_i, NULL_COLUMN_VALUE_R, &
                         obsdat%numBody, obsdat%numBody_max)
   end subroutine obs_bodySet_i


   subroutine obs_bodySet_r4(obsdat, column_index, row_index, value_r4)
      !
      ! :Purpose: Set a real-valued body observation-data element.
      !      To control access to the observation object.  Sets the (real)
      !      value of the row_index'th ObsData element with the indicated column
      !      index from the "body".
      !
      implicit none
      type(struct_obs), intent(inout)  :: obsdat
      integer         , intent(in)     :: column_index
      integer         , intent(in)     :: row_index
      real(4)         , intent(in)     :: value_r4

      real(pre_obsReal)                :: value_r

      value_r = value_r4

      call odc_columnSet(obsdat%realBodies, column_index, row_index, &
                         NULL_COLUMN_VALUE_I, value_r, &
                         obsdat%numBody, obsdat%numBody_max)
   end subroutine obs_bodySet_r4


   subroutine obs_bodySet_r8(obsdat, column_index, row_index, value_r8)
      !
      ! :Purpose: Set a real-valued body observation-data element.
      !      To control access to the observation object.  Sets the (real)
      !      value of the row_index'th ObsData element with the indicated column
      !      index from the "body".
      !
      implicit none
      type(struct_obs), intent(inout)  :: obsdat
      integer         , intent(in)     :: column_index
      integer         , intent(in)     :: row_index
      real(8)         , intent(in)     :: value_r8

      real(pre_obsReal)                :: value_r

      value_r = value_r8

      call odc_columnSet(obsdat%realBodies, column_index, row_index, &
                         NULL_COLUMN_VALUE_I, value_r, &
                         obsdat%numBody, obsdat%numBody_max)
   end subroutine obs_bodySet_r8


   subroutine obs_class_initialize(obsColumnMode_in, myip)
      !
      ! :Purpose: Set observation-data class variables.
      !      Set variables that take the same value for all instances of the
      !      class.
      !
      implicit none
      ! mode controlling the subset of columns that are activated in all objects
      character(len=*), intent(in), optional :: obsColumnMode_in
      integer, intent(in), optional :: myip

      integer :: myip_

      ! The 'save' makes this a CLASS-CONSTANT variable
      character(len=12), save :: obsColumnMode_class = '            '

      INITIALIZED: if(.not. obs_class_initialized) then
         obs_class_initialized = .true.

         if(present(myip)) then
            myip_ = myip
         else
            ! Make an assumption:  cause every processor to write much to stdout
            myip_ = 0
         end if

         ! Determine which columns will be initially active
         if(present(obsColumnMode_in)) then
            obsColumnMode_class=trim(obsColumnMode_in)
         else
            obsColumnMode_class='ALL '
         endif

         write(*,*)'OBS_CLASS_INITIALIZE: obsColumnMode=', &
                   trim(obsColumnMode_class)

         ! DELEGATE THE REAL WORK TO THE ODC CLASS
         ! DELEGATE THE REAL WORK TO THE ODC CLASS
         call odc_class_initialize(obsColumnMode_class)

      else INITIALIZED
         write(*,*) 'obs_class_initialize: !!! WARNING WARNING WARNING!!!'
         write(*,*) 'obs_class_initialize: already called before, not ' &
                    // 're-activating columns'
         write(*,*) 'obs_class_initialize: !!! WARNING WARNING WARNING!!!'
         if(present(obsColumnMode_in)) then
            if(trim(obsColumnMode_in) /= trim(obsColumnMode_class)) then
               call obs_abort('obs_class_initialize: called with different '&
                            //'value of obsColumnMode than during first call: '&
                            // trim(obsColumnMode_class) // ' /= ' &
                            // trim(obsColumnMode_in))
               return
            endif
         endif
      endif INITIALIZED
   end subroutine obs_class_initialize


   subroutine obs_clean(obsdat,hx,nens,nobsout,qcvar,checkZha_opt)
      !
      ! :Purpose: remove all observations from the obsdat  
      !           that will not be assimilated. 
      !
      ! :Arguments:
      !           :nobsout: unit number for the ASCII output
      !           :qcvar:   input logical indicating if the input obsdat 
      !                     data have benefited from a qc-var procedure
      ! :The logic applied:
      !     A body (and its associated header)
      !     will be retained if these three conditions are all met:
      !
      !         - 1) either of:
      !
      !             - 1a) btest(obsdat%intBodies%columns(OBS_FLG,jdata),12)
      !             - 1b) .not. qcvar (the 5th parameter of obs_clean)
      !         - 2) obsdat% intBodies%columns(OBS_ASS,jdata) == 1
      !         - 3) obsdat%realBodies%columns(OBS_ZHA,jdata) >= 0.0
      !
      implicit none

      type (struct_obs), intent(inout) :: obsdat

      real(8), intent(inout) :: hx(:,:)
      integer, intent(in)    :: nens
      integer, intent(in)    :: nobsout
      logical, intent(in)    :: qcvar
      logical, intent(in), optional :: checkZha_opt

      integer :: iaccept,idata,ipnt,iwrite
      integer :: jdata,kobs,var3d,kobsout
      integer :: column_index
      integer :: active_index
      logical :: checkZha

      if (nobsout > 0) write(nobsout,'(1x,A,I7)') 'stations prior to cleanup: ', obsdat%numHeader
      write(*,*) 'enter obs_clean'

      ! User can choose if check on negativity of OBS_ZHA is done (is done by default)
      if (present(checkZha_opt)) then
        checkZha = checkZha_opt
      else
        checkZha = .true.
      end if

      kobsout=0 
      iwrite=0
      stations: do kobs=1,obsdat%numHeader
         ipnt  = obs_headElem_i(obsdat, OBS_RLN, kobs)
         idata = obs_headElem_i(obsdat, OBS_NLV, kobs)
         iaccept=0
         observations: do jdata = ipnt, ipnt + idata - 1
            if (     btest(obs_bodyElem_i(obsdat, OBS_FLG, jdata),12) &
                .or. .not. qcvar) then 
               ! data will be accepted if they went through the variational
               ! system  including the qcvar. They will also be accepted if the
               ! qcvar procedure was not applied (i.e. when backalt files are
               ! used as input).
               var3d=1
            else
               var3d=0
            endif

            ! To remove observations for which the height in the atmosphere has
            ! not been assigned (for instance because they are above the model
            ! top for the EnKF system)
            if (checkZha .and. obs_bodyElem_r(obsdat, OBS_ZHA, jdata) < 0.) then
               call obs_bodySet_i(obsdat, OBS_ASS, jdata, -1)
            endif

            if (      (obs_bodyElem_i(obsdat, OBS_ASS, jdata) == obs_assimilated) &
                 .and.(var3d == 1)) then 
               ! the observation will be used in the analysis
               iaccept=iaccept+1
               iwrite=iwrite+1
               do active_index=1,odc_numActiveColumn(obsdat%intBodies)
                  column_index=odc_columnIndexFromActiveIndex( &
                                       obsdat%intBodies%odc_flavour,active_index)
                     call obs_bodySet_i (obsdat, column_index, iwrite, &
                          obs_bodyElem_i(obsdat, column_index, jdata))
               enddo
               do active_index=1,odc_numActiveColumn(obsdat%realBodies)
                  column_index=odc_columnIndexFromActiveIndex( &
                                      obsdat%realBodies%odc_flavour,active_index)
                     call obs_bodySet_r (obsdat, column_index, iwrite, &
                          obs_bodyElem_r(obsdat, column_index, jdata))
               enddo
               hx(1:nens,iwrite)=hx(1:nens,jdata)
               ! Revise the header row cross-index to match the revised headers
               call obs_bodySet_i (obsdat, OBS_HIND, iwrite, kobsout+1)
            endif
         enddo observations

         ! adjust obsdat%realHeaders%columns 
         if (iaccept > 0) then
            kobsout=kobsout+1
            do active_index=1,odc_numActiveColumn(obsdat%intHeaders)
               column_index = odc_columnIndexFromActiveIndex( &
                                      obsdat%intHeaders%odc_flavour,active_index)
               call obs_headSet_i (obsdat, column_index, kobsout, &
                    obs_headElem_i(obsdat, column_index, kobs))
            enddo
            do active_index=1,odc_numActiveColumn(obsdat%realHeaders)
               column_index = odc_columnIndexFromActiveIndex( &
                                     obsdat%realHeaders%odc_flavour,active_index)
               call obs_headSet_r (obsdat, column_index, kobsout, &
                    obs_headElem_r(obsdat, column_index, kobs))
            enddo
            obsdat%headerPrimaryKey(kobsout) = obsdat%headerPrimaryKey(kobs)
            obsdat%cstnid(kobsout)=obsdat%cstnid(kobs)
            obsdat%cfamily(kobsout)=obsdat%cfamily(kobs)
            ! Revise the body cross-indices to match the revised bodies
            call obs_headSet_i(obsdat, OBS_NLV, kobsout, iaccept)
            call obs_headSet_i(obsdat, OBS_RLN, kobsout, iwrite-iaccept+1)
         endif
      enddo stations
      obsdat%numHeader=kobsout
      obsdat%numBody = iwrite

      if (nobsout > 0) then
        write(nobsout,'(1x,A)') 'after cleanup of the cma: '
        write(nobsout,'(1x,A,I7)') &
             'number of stations containing valid data   ',obsdat%numHeader
        write(nobsout,'(1x,A,I7)') & 
             'number of observations now in the cma file ',obsdat%numBody
      end if

   end subroutine obs_clean


   subroutine obs_clean2(obsdat)
     !
     ! :Purpose: remove all observations from the obsdat  
     !           that will not be assimilated. 
     !           modified version of obs_clean, used by MIDAS
     !
     implicit none

     type (struct_obs), intent(inout) :: obsdat

     integer :: numBodyAccept,bodyIndexEnd,bodyIndexBeg,numBodyCleaned
     integer :: bodyIndex,headerIndex,numHeaderCleaned
     integer :: column_index
     integer :: active_index

     write(*,*) 'before cleanup of obsSpaceData:'
     write(*,*) 'number of headers ', obsdat%numHeader
     write(*,*) 'number of bodies  ', obsdat%numBody

     numHeaderCleaned = 0 
     numBodyCleaned = 0
     stations: do headerIndex = 1, obsdat%numHeader
       bodyIndexBeg = obs_headElem_i(obsdat, OBS_RLN, headerIndex)
       bodyIndexEnd = obs_headElem_i(obsdat, OBS_NLV, headerIndex) + bodyIndexBeg - 1
       numBodyAccept = 0
       observations: do bodyIndex = bodyIndexBeg, bodyIndexEnd
         if ( obs_bodyElem_i(obsdat,OBS_ASS,bodyIndex) == obs_assimilated ) then
           ! the observation will be used in the analysis
           numBodyAccept = numBodyAccept + 1
           numBodyCleaned = numBodyCleaned + 1
           do active_index = 1, odc_numActiveColumn(obsdat%intBodies)
             column_index = odc_columnIndexFromActiveIndex( &
                              obsdat%intBodies%odc_flavour,active_index)
             call obs_bodySet_i (obsdat, column_index, numBodyCleaned, &
                  obs_bodyElem_i(obsdat, column_index, bodyIndex))
           enddo
           do active_index = 1, odc_numActiveColumn(obsdat%realBodies)
             column_index = odc_columnIndexFromActiveIndex( &
                                 obsdat%realBodies%odc_flavour,active_index)
             call obs_bodySet_r (obsdat, column_index, numBodyCleaned, &
                                 obs_bodyElem_r(obsdat, column_index, bodyIndex))
           enddo
           ! Revise the header row cross-index to match the revised headers
           call obs_bodySet_i (obsdat, OBS_HIND, numBodyCleaned, numHeaderCleaned+1)
         endif
       enddo observations

       ! adjust obsdat%realHeaders%columns 
       if (numBodyAccept > 0) then
         numHeaderCleaned = numHeaderCleaned + 1
         do active_index = 1, odc_numActiveColumn(obsdat%intHeaders)
           column_index = odc_columnIndexFromActiveIndex( &
                                  obsdat%intHeaders%odc_flavour,active_index)
           call obs_headSet_i (obsdat, column_index, numHeaderCleaned, &
                obs_headElem_i(obsdat, column_index, headerIndex))
         enddo
         do active_index = 1, odc_numActiveColumn(obsdat%realHeaders)
           column_index = odc_columnIndexFromActiveIndex( &
                                 obsdat%realHeaders%odc_flavour,active_index)
           call obs_headSet_r (obsdat, column_index, numHeaderCleaned, &
                obs_headElem_r(obsdat, column_index, headerIndex))
         enddo
         obsdat%headerPrimaryKey(numHeaderCleaned) = obsdat%headerPrimaryKey(headerIndex)
         obsdat%cstnid(numHeaderCleaned) = obsdat%cstnid(headerIndex)
         obsdat%cfamily(numHeaderCleaned) = obsdat%cfamily(headerIndex)
         ! Revise the body cross-indices to match the revised bodies
         call obs_headSet_i(obsdat, OBS_NLV, numHeaderCleaned, numBodyAccept)
         call obs_headSet_i(obsdat, OBS_RLN, numHeaderCleaned, numBodyCleaned-numBodyAccept+1)
       endif
     enddo stations
     obsdat%numHeader = numHeaderCleaned
     obsdat%numBody = numBodyCleaned

     ! reallocate obsdat to free up memory
     call obs_squeeze(obsdat)

     write(*,*) 'after cleanup of obsSpaceData: '
     write(*,*) 'number of headers with valid data   ',obsdat%numHeader
     write(*,*) 'number of bodies kept               ',obsdat%numBody

   end subroutine obs_clean2

   function obs_columnActive_IB(obsdat,column_index) result(columnActive)
      !
      ! :Purpose:
      !      Return the active status for a column
      !
      implicit none
      type (struct_obs), intent(inout) :: obsdat
      integer :: column_index
      logical :: columnActive

      columnActive = obsdat%intBodies%odc_flavour%columnActive(column_index)

   end function obs_columnActive_IB


   function obs_columnActive_IH(obsdat,column_index) result(columnActive)
      !
      ! :Purpose:
      !      Return the active status for a column
      !
      implicit none
      type (struct_obs), intent(inout) :: obsdat
      integer :: column_index
      logical :: columnActive

      columnActive = obsdat%intHeaders%odc_flavour%columnActive(column_index)

   end function obs_columnActive_IH


   function obs_columnActive_RB(obsdat,column_index) result(columnActive)
      !
      ! :Purpose:
      !      Return the active status for a column
      !
      implicit none
      type (struct_obs), intent(inout) :: obsdat
      integer :: column_index
      logical :: columnActive

      columnActive = obsdat%realBodies%odc_flavour%columnActive(column_index)

   end function obs_columnActive_RB


   function obs_columnActive_RH(obsdat,column_index) result(columnActive)
      !
      ! :Purpose:
      !      Return the active status for a column
      !
      implicit none
      type (struct_obs), intent(inout) :: obsdat
      integer :: column_index
      logical :: columnActive

      columnActive = obsdat%realHeaders%odc_flavour%columnActive(column_index)

   end function obs_columnActive_RH


   function obs_columnIndexFromName(column_name) result(column_index_out)
      !
      ! :Purpose:
      !      Situations do occur where the client knows only the name of a
      !      column, but needs to know its index. This method supplies the index.
      !
      implicit none
      character(len=*)        , intent(in) :: column_name
      integer                              :: column_index_out

      integer            :: column_index
      logical            :: lfound
      character(len=100) :: message

      lfound=.false.

      ! check integer-header
      do column_index=odc_flavour_IH%ncol_beg, odc_flavour_IH%ncol_end
         if(trim(column_name) == &
            trim(odc_flavour_IH%columnNameList(column_index)))then
            lfound=.true.
            column_index_out = column_index
            exit
         endif
      enddo
      if (lfound) return

      ! check real-header
      do column_index=odc_flavour_RH%ncol_beg, odc_flavour_RH%ncol_end
         if(trim(column_name) == &
            trim(odc_flavour_RH%columnNameList(column_index)))then
            lfound=.true.
            column_index_out = column_index
            exit
         endif
      enddo
      if (lfound) return

      ! check integer-body
      do column_index=odc_flavour_IB%ncol_beg, odc_flavour_IB%ncol_end
         if(trim(column_name) == &
            trim(odc_flavour_IB%columnNameList(column_index)))then
            lfound=.true.
            column_index_out = column_index
            exit
         endif
      enddo
      if (lfound) return

      ! check real-body
      do column_index=odc_flavour_RB%ncol_beg, odc_flavour_RB%ncol_end
         if(trim(column_name) == &
            trim(odc_flavour_RB%columnNameList(column_index)))then
            lfound=.true.
            column_index_out = column_index
            exit
         endif
      enddo
      if (lfound) return

      write(message,*)'abort in obs_columnIndexFromName: ' // &
                      'name not found='// trim(column_name)
      call obs_abort(message); return

   end function obs_columnIndexFromName

  function obs_isColumnNameValid(column_name) result(isValid)  
    !
    ! :Purpose: Check if the obsSpaceData column name is valid.
    !      
    implicit none

    ! arguements:
    character(len=*), intent(in)  :: column_name
    logical                       :: isValid

    ! locals: 
    integer                       :: column_index

    isValid = .false.

    ! check "STID" 
    if (trim(column_name) == "STID") isValid = .true.
 
    if (isValid) return

    ! check integer-header
    do column_index = odc_flavour_IH%ncol_beg, odc_flavour_IH%ncol_end
      if (trim(column_name) == &
           trim(odc_flavour_IH%columnNameList(column_index))) then
        isValid = .true.
        exit
      end if
    end do
    if (isValid) return

    ! check real-header
    do column_index = odc_flavour_RH%ncol_beg, odc_flavour_RH%ncol_end
      if (trim(column_name) == &
           trim(odc_flavour_RH%columnNameList(column_index))) then
        isValid = .true.
        exit
      end if
    end do
    if (isValid) return
  
    ! check integer-body
    do column_index = odc_flavour_IB%ncol_beg, odc_flavour_IB%ncol_end
      if (trim(column_name) == &
           trim(odc_flavour_IB%columnNameList(column_index))) then
        isValid = .true.
        exit
      end if
    end do
    if (isValid) return

    ! check real-body
    do column_index = odc_flavour_RB%ncol_beg, odc_flavour_RB%ncol_end
      if (trim(column_name) == &
           trim(odc_flavour_RB%columnNameList(column_index))) then
        isValid = .true.
        exit
      end if
    end do
    if (isValid) return
   
  end function obs_isColumnNameValid


   function obs_columnIndexFromNameForFlavour(odc_flavour, column_name) &
                                                         result(column_index_out)
      !
      ! :Purpose:
      !      Situations do occur where the client knows only the name of a
      !      column, but needs to know its index. This method supplies the index.
      !
      implicit none
      type(struct_odc_flavour), intent(in) :: odc_flavour
      character(len=*)        , intent(in) :: column_name
      integer                              :: column_index_out

      integer            :: column_index
      logical            :: lfound
      character(len=100) :: message

      lfound=.false.
      do column_index=odc_flavour%ncol_beg, odc_flavour%ncol_end
         if(trim(column_name) == &
            trim(odc_flavour%columnNameList(column_index)))then
            lfound=.true.
            column_index_out = column_index
            exit
         endif
      enddo

      if(.not.lfound) then
         write(message,*)'abort in obs_columnIndexFromNameForFlavour [' &
                       // odc_flavour%dataType //','// odc_flavour%headOrBody //&
                       ']: name not found='// column_name
         call obs_abort(message); return
      end if
   end function obs_columnIndexFromNameForFlavour


   function obs_columnIndexFromName_IB(column_name) result(column_index)
      !
      ! :Purpose:
      !      This wrapper around obs_columnIndexFromName selects the data-column
      !      flavour.
      !
      implicit none
      character(len=*), intent(in) :: column_name
      integer                      :: column_index

      column_index = obs_columnIndexFromNameForFlavour(odc_flavour_IB, column_name)
   end function obs_columnIndexFromName_IB


   function obs_columnIndexFromName_IH(column_name) result(column_index)
      !
      ! :Purpose:
      !      This wrapper around obs_columnIndexFromName selects the data-column
      !      flavour.
      !
      implicit none
      character(len=*), intent(in) :: column_name
      integer                      :: column_index

      column_index = obs_columnIndexFromNameForFlavour(odc_flavour_IH, column_name)
   end function obs_columnIndexFromName_IH


   function obs_columnIndexFromName_RB(column_name) result(column_index)
      !
      ! :Purpose:
      !      This wrapper around obs_columnIndexFromName selects the data-column
      !      flavour.
      !
      implicit none
      character(len=*), intent(in) :: column_name
      integer                      :: column_index

      column_index = obs_columnIndexFromNameForFlavour(odc_flavour_RB, column_name)
   end function obs_columnIndexFromName_RB


   function obs_columnIndexFromName_RH(column_name) result(column_index)
      !
      ! :Purpose:
      !      This wrapper around obs_columnIndexFromName selects the data-column
      !      flavour.
      !
      implicit none
      character(len=*), intent(in) :: column_name
      integer                      :: column_index

      column_index = obs_columnIndexFromNameForFlavour(odc_flavour_RH, column_name)
   end function obs_columnIndexFromName_RH


   function obs_columnDataType(columnIndex) result(dataType)
     !
     ! :Purpose: return the data type of column, either 'real' or 'integer'
     !
     implicit none

     integer, intent(in) :: columnIndex
     character(len=7)    :: dataType

     if (columnIndex >= NHDR_INT_BEG .and. columnIndex <= NHDR_INT_END) then
       dataType = 'integer'
     else if (columnIndex >= NHDR_REAL_BEG .and. columnIndex <= NHDR_REAL_END) then
       dataType = 'real'
     else if (columnIndex >= NBDY_INT_BEG .and. columnIndex <= NBDY_INT_END) then
       dataType = 'integer'
     else if (columnIndex >= NBDY_REAL_BEG .and. columnIndex <= NBDY_REAL_END) then
       dataType = 'real'
     else
       call obs_abort('obs_columnDataType: invalid value of column index')
     end if

   end function obs_columnDataType


   subroutine obs_copy( obs_a, obs_b )
      !
      ! :Purpose: copy an obsdat object
      !
      ! :Arguments:
      !           :obs_a: input object
      !           :obs_b: a copy of obs_a
      !
      ! :Note: this method assumes that obs_b has already been initialized
      !
      implicit none

      type(struct_obs), intent(in)    :: obs_a
      type(struct_obs), intent(inout) :: obs_b

      integer :: column_index

      ! check if object to be copied is empty and react appropriately
      if (obs_a%numHeader_max == 0 .or. obs_a%numBody_max == 0) then
         obs_b%numHeader     = obs_a%numHeader
         obs_b%numHeader_max = obs_a%numHeader_max
         obs_b%numBody       = obs_a%numBody
         obs_b%numBody_max   = obs_a%numBody_max
         return
      endif

      !** Commented out by M. Buehner to allow use in EnVar (also added copy of 
      !** headerIndex_mpiglobal and bodyIndex_mpiglobal, if they exist)
      !if(obs_a%mpi_local)then
      !   call obs_abort( &
      !         'obs_copy() is not equipped to handle the case, mpi_local=.true.')
      !   return
      !end if

      do column_index=NHDR_REAL_BEG,NHDR_REAL_END
         if(obs_a%realHeaders%odc_flavour%columnActive(column_index))  &
              obs_b%realHeaders%columns(column_index)%value_r(:) &
            = obs_a%realHeaders%columns(column_index)%value_r(:)
      enddo
      do column_index=NHDR_INT_BEG,NHDR_INT_END
         if(obs_a%intHeaders%odc_flavour%columnActive(column_index))  &
              obs_b%intHeaders%columns(column_index)%value_i(:) &
            = obs_a%intHeaders%columns(column_index)%value_i(:)
      enddo
      do column_index=NBDY_REAL_BEG,NBDY_REAL_END
         if(obs_a%realBodies%odc_flavour%columnActive(column_index))  &
              obs_b%realBodies%columns(column_index)%value_r(:) &
            = obs_a%realBodies%columns(column_index)%value_r(:)
      enddo
      do column_index=NBDY_INT_BEG,NBDY_INT_END
         if(obs_a%intBodies%odc_flavour%columnActive(column_index))  &
              obs_b%intBodies%columns(column_index)%value_i(:) &
            = obs_a%intBodies%columns(column_index)%value_i(:)
      enddo
      obs_b%headerPrimaryKey(:) = obs_a%headerPrimaryKey(:)
      obs_b%cstnid(:)     = obs_a%cstnid(:)
      obs_b%cfamily(:)    = obs_a%cfamily(:)

      obs_b%numHeader     = obs_a%numHeader
      obs_b%numHeader_max = obs_a%numHeader_max
      obs_b%numBody       = obs_a%numBody
      obs_b%numBody_max   = obs_a%numBody_max

      if(associated(obs_a%headerIndex_mpiglobal)) then
        write(*,*) 'obs_copy: copying headerIndex_mpiglobal'
        if(associated(obs_b%headerIndex_mpiglobal)) then
          deallocate(obs_b%headerIndex_mpiglobal)
        endif
        allocate(obs_b%headerIndex_mpiglobal(size(obs_a%headerIndex_mpiglobal)))
        obs_b%headerIndex_mpiglobal(:) = obs_a%headerIndex_mpiglobal(:)
      endif

      if(associated(obs_a%bodyIndex_mpiglobal)) then
        write(*,*) 'obs_copy: copying bodyIndex_mpiglobal'
        if(associated(obs_b%bodyIndex_mpiglobal)) then
          deallocate(obs_b%bodyIndex_mpiglobal)
        endif
        allocate(obs_b%bodyIndex_mpiglobal(size(obs_a%bodyIndex_mpiglobal)))
        obs_b%bodyIndex_mpiglobal(:) = obs_a%bodyIndex_mpiglobal(:)
      endif

      obs_b%mpi_local     = obs_a%mpi_local
   end subroutine obs_copy


   subroutine obs_count_headers( obsdat, kulout )
     !
     ! :Purpose: count the number of stations and
     !           observations that are in the obsdat.
     !
     ! :Arguments:
     !           :obsdat: obsSpaceData object
     !           :kulout: unit number for ASCII error messages and 
     !                    observation counts.
     !
      implicit none

     ! Arguments: 
      type (struct_obs), intent(in) :: obsdat
      integer          , intent(in) :: kulout

      integer, parameter :: MAXID = 256
      integer  :: allstn,allobs,id,idata,kobs
      integer, dimension(MAXID)   :: numobs,numstn
      character(len=100) :: message

      ! initialize totals to zero
      numstn(:)=0
      numobs(:)=0
      allstn=0
      allobs=0

      do kobs=1,obsdat%numHeader
         id=obs_headElem_i(obsdat, OBS_ITY, kobs)
         if(id > MAXID) then
            id=mod(id,1000)
         endif
         if ((id < 1).or.(id > MAXID)) then 
            write(message,*)'OBS_COUNT_HEADERS: ITY (instrument and ' // &
                            'retrieval type) out of range: ', id
            write(kulout,*)message
            call obs_abort(message); return
         endif
         numstn(id)=numstn(id)+1
         ! idata: number of obs for this station
         idata = obs_headElem_i(obsdat, OBS_NLV, kobs)
         numobs(id)=numobs(id)+idata
      enddo

      write(kulout,*) 'number of stations and observations'
      write(kulout,*) ' idtype #stations #observations '
      do id=1,MAXID
         if (numstn(id) > 0) then
            write(kulout,'(i3,3x,i7,2x,i8)') id,numstn(id),numobs(id)
         endif
         allstn=allstn+numstn(id)
         allobs=allobs+numobs(id)
      enddo
      write(kulout,'(1x,A,I7)') 'total number of stations:     ',allstn
      write(kulout,'(1x,A,I7)') 'total number of observations: ',allobs

      return
   end subroutine obs_count_headers


   subroutine obs_deallocate(obsdat)
      !
      ! :Purpose:
      !      De-allocate arrays.  This is a private method.
      !
      implicit none

      type(struct_obs), intent(inout) :: obsdat

      integer :: ierr
      integer :: column_index

      if(.not.obsdat%allocated) then
         ! The object content has already been de-allocated.  Don't repeat it.
         return
      endif

      obsdat%allocated=.false.
      HEADER:if(obsdat%numHeader_max > 0) then
         if (associated(obsdat%headerPrimaryKey)) then
            deallocate(obsdat%headerPrimaryKey,STAT=ierr)
            nullify(obsdat%headerPrimaryKey)
            if(ierr /= 0)write(*,*) 'Problem detected with headerPrimaryKey. IERR =',ierr
         end if
          
         if (associated(obsdat%cfamily)) then
            deallocate(obsdat%cfamily,STAT=ierr)
            nullify(obsdat%cfamily)
            if(ierr /= 0)write(*,*) 'Problem detected with CFAMILY. IERR =',ierr
         end if

         if (associated(obsdat%cstnid))then
            deallocate(obsdat%cstnid,STAT=ierr)
            nullify(obsdat%cstnid)
            if(ierr /= 0)write(*,*) 'Problem detected with CSTNID. IERR =',ierr
         end if

         if (associated(obsdat%realHeaders%columns))then
            do column_index=NHDR_REAL_BEG,NHDR_REAL_END
               if(obsdat%realHeaders%odc_flavour%columnActive(column_index)) &
                  call odc_deallocate(obsdat%realHeaders%columns(column_index))
            enddo
         end if

         if (associated(obsdat%intHeaders%columns))then
            do column_index=NHDR_INT_BEG,NHDR_INT_END
               if(obsdat%intHeaders%odc_flavour%columnActive(column_index)) &
                  call odc_deallocate(obsdat%intHeaders%columns(column_index))
            enddo
         end if

         deallocate(obsdat%scratchRealHeader)
         nullify   (obsdat%scratchRealHeader)

         deallocate(obsdat%scratchIntHeader)
         nullify   (obsdat%scratchIntHeader)
      endif HEADER

      BODY:if(obsdat%numBody_max > 0) then
         if (associated(obsdat%bodyPrimaryKey)) then
            deallocate(obsdat%bodyPrimaryKey,STAT=ierr)
            nullify(obsdat%bodyPrimaryKey)
            if(ierr /= 0)write(*,*) 'Problem detected with bodyPrimaryKey. IERR =',ierr
         end if

         if (associated(obsdat%realBodies%columns))then
            do column_index=NBDY_REAL_BEG,NBDY_REAL_END
               if(obsdat%realBodies%odc_flavour%columnActive(column_index)) &
                  call odc_deallocate(obsdat%realBodies%columns(column_index))
            enddo
         end if

         if (associated(obsdat%intBodies%columns))then
            do column_index=NBDY_INT_BEG,NBDY_INT_END
               if(obsdat%intBodies%odc_flavour%columnActive(column_index)) &
                  call odc_deallocate(obsdat%intBodies%columns(column_index))
            enddo
         end if

         deallocate(obsdat%scratchRealBody)
         nullify   (obsdat%scratchRealBody)

         deallocate(obsdat%scratchIntBody)
         nullify   (obsdat%scratchIntBody)
      endif BODY

      obsdat%numHeader_max=0
      obsdat%numBody_max=0

      call ild_finalize(obsdat%header_index_list_depot)
      call ild_finalize(obsdat%body_index_list_depot)

   end subroutine obs_deallocate


   function obs_elem_c(obsdat,name,row_index) result(value)
      !
      ! :Purpose: Get a character-string-valued observation-data element.
      !      To control access to the observation object.  Returns the
      !      (character) value of the ObsData element with the indicated name
      !      and row_index.
      !
      implicit none

      character(len=12) :: value
      type(struct_obs), intent(in)  :: obsdat
      character(len=*), intent(in)  :: name
      integer         , intent(in)  :: row_index

      select case (trim(name))
      case ('STID'); value=obsdat%cstnid(row_index)

      case default
         call obs_abort('ERROR:  ' // name(1:4) // &
                        ' is not a character(len=*) observation.');return
      end select
   end function obs_elem_c


   subroutine obs_enkf_prntbdy(obsdat,kstn,kulout)
      !
      ! :Purpose:  print all data records associated with an observation
      !
      ! :Arguments:
      !     :kstn: no. of station 
      !     :kulout: unit used for printing
      !
      implicit none

      type(struct_obs), intent(in) :: obsdat
      integer,      intent(in)  :: kstn, kulout

      integer :: ipnt, idata, idata2, jdata, var3d

      ! general information

      ipnt  = obs_headElem_i(obsdat, OBS_RLN, kstn)
      idata = obs_headElem_i(obsdat, OBS_NLV, kstn)

      if(idata == 1) then
         write(kulout,fmt=9101)idata,kstn
      else
         write(kulout,fmt=9100)idata,kstn
      end if
9100  format(2x,'there are ',i6,1x,'data in record no.',1x,i6)
9101  format(2x,'there is ', i6,1x,'data in record no.',1x,i6)

      ! print all data records

      write(kulout,'(a,a,a,a)') '   no.   var.  press. ass observ. ', &
         '  o minus p  o minus p6  o minus a   o min a0     obserr. root(hpht) ', &
         'root(haht) innov std dev obs std dev   acc ', &
         'zhad  vco'
      do jdata = ipnt, ipnt + idata - 1
         idata2 = jdata -ipnt + 1
         if (btest(obsdat%intBodies%columns(OBS_FLG)%value_i(jdata),12)) then 
            var3d=1
         else
            var3d=0
         endif

         write(kulout,fmt=9201) idata2,obsdat%intBodies%columns(OBS_VNM)%value_i(jdata), &
            obsdat%realBodies%columns(OBS_PPP )%value_r(jdata), &
            obsdat%intBodies%columns(OBS_ASS )%value_i(jdata), &
            obsdat%realBodies%columns(OBS_VAR )%value_r(jdata), &
            obsdat%realBodies%columns(OBS_OMP )%value_r(jdata), &
            obsdat%realBodies%columns(OBS_OMP6)%value_r(jdata), &
            obsdat%realBodies%columns(OBS_OMA )%value_r(jdata), &
            obsdat%realBodies%columns(OBS_OMA0)%value_r(jdata), &
            obsdat%realBodies%columns(OBS_OER )%value_r(jdata), &
            obsdat%realBodies%columns(OBS_HPHT)%value_r(jdata), &
            obsdat%realBodies%columns(OBS_HAHT)%value_r(jdata), &
            obsdat%realBodies%columns(OBS_SIGI)%value_r(jdata), &
            obsdat%realBodies%columns(OBS_SIGO)%value_r(jdata), &
            var3d, &
            obsdat%realBodies%columns(OBS_ZHA )%value_r(jdata), &
            obsdat%intBodies%columns(OBS_VCO )%value_i(jdata)
      enddo

9201  format(1x,i6,1x,i6,1x,f7.0,1x,i3,10(1x,f10.3),1x,i2, &
           1x,f10.3,1x,i2)

      return

   end subroutine obs_enkf_prntbdy


   subroutine obs_enkf_prnthdr(obsdata,kobs,kulout)
      !
      ! :Purpose: printing of the header of an observation record
      !
      ! :Arguments:
      !           :kobs: no. of observation
      !           :kulout: unit used for optional printing
      !
      implicit none

      type(struct_obs), intent(in) :: obsdata
      integer,      intent(in)  :: kobs, kulout

      integer :: obsRLN, obsONM, obsDAT, obsETM, obsINS, obsOTP, obsITY
      integer :: obsNLV, obsPAS, obsREG, obsIP
      real(pre_obsReal) ::obsLAT, obsLON, obsALT, obsBX, obsBY, obsBZ, obsAZA

      ! general information

      write(kulout,fmt=9100)kobs,obsdata%cstnid(KOBS)
9100  format(//,2x,'-- observation record no.' &
         ,1x,i6,3x,' station id:',A12)

      ! print header's content

9202  format(2x,'position within realBodies:',i6)
      obsRLN = obs_headElem_i(obsdata, OBS_RLN, kobs)
      obsONM = obs_headElem_i(obsdata, OBS_ONM, kobs)
      obsDAT = obs_headElem_i(obsdata, OBS_DAT, kobs)
      obsETM = obs_headElem_i(obsdata, OBS_ETM, kobs)
      obsINS = obs_headElem_i(obsdata, OBS_INS, kobs)
      obsOTP = obs_headElem_i(obsdata, OBS_OTP, kobs)
      obsITY = obs_headElem_i(obsdata, OBS_ITY, kobs)
      obsLAT = obs_headElem_r(obsdata, OBS_LAT, kobs)
      obsLON = obs_headElem_r(obsdata, OBS_LON, kobs)
      obsALT = obs_headElem_r(obsdata, OBS_ALT, kobs)
      obsBX  = obs_headElem_r(obsdata, OBS_BX , kobs)
      obsBY  = obs_headElem_r(obsdata, OBS_BY , kobs)
      obsBZ  = obs_headElem_r(obsdata, OBS_BZ , kobs)
      obsAZA = obs_headElem_r(obsdata, OBS_AZA, kobs)
      write(kulout,fmt=9200) obsRLN, obsONM, obsDAT, obsETM, &
           obsINS, obsOTP, obsITY, obsLAT, obsLON, obsALT, &
           obsBX , obsBY , obsBZ , obsAZA

      obsNLV = obs_headElem_i(obsdata, OBS_NLV, kobs)
      obsPAS = obs_headElem_i(obsdata, OBS_PAS, kobs)
      obsREG = obs_headElem_i(obsdata, OBS_REG, kobs)
      obsIP  = obs_headElem_i(obsdata, OBS_IP , kobs)
      write(kulout,fmt=9201) obsNLV, obsPAS, obsREG, obsIP
         

9200  format(2x,'position within realBodies:',i8,1x,'stn. number:',i6,1x,/, &
         '  date: ',i10,1x,' time: ',i8,/, &
         '  model box:',i12,1x,'instrument: ',i6,1x, &
         'obs. type:',i8,1x,/, &
         '  (lat,lon):',f12.6,1x,f12.6,1x, &
         'stations altitude:',f12.6,1x,/,2x, &
         'block location: ',3(f12.6,1x),/,2x, &
         'azimuth angle:',f12.6)
9201  format('  number of data:',i6,1x, &
         ' pass: ',i6,' region: ',i6,/,2x,'processor: ',i6)

      return
   end subroutine obs_enkf_prnthdr


   subroutine obs_expandToMpiGlobal(obsdat)
      !
      ! :Purpose: restore Global array realBodies and intBodies.
      !      To reconstitute the mpi-global observation object by gathering the
      !      necessary data from all processors (to all processors).
      !
      ! :Note: for the character data cstnid(:), this is converted to integers
      !       with IACHAR and back to characters with ACHAR, to facilitate this
      !       gather through rpn_comm_allreduce
      !
      implicit none

      type(struct_obs), intent(inout) :: obsdat

      integer, allocatable :: headerIndex_mpiglobal(:),all_headerIndex_mpiglobal(:,:)
      integer, allocatable :: bodyIndex_mpiglobal(:),all_bodyIndex_mpiglobal(:,:)
      integer(8), allocatable :: headerPrimaryKey_mpilocal(:), all_headerPrimaryKey_mpilocal(:,:)
      integer(8), allocatable :: bodyPrimaryKey_mpilocal(:), all_bodyPrimaryKey_mpilocal(:,:)
      integer, allocatable :: intHeaders_mpilocal(:,:),all_intHeaders_mpilocal(:,:,:)
      real(pre_obsReal), allocatable :: realHeaders_mpilocal(:,:),all_realHeaders_mpilocal(:,:,:)
      integer, allocatable :: intStnid_mpilocal(:,:),all_intStnid_mpilocal(:,:,:)
      integer, allocatable :: intFamily_mpilocal(:,:),all_intFamily_mpilocal(:,:,:)
      integer, allocatable :: intBodies_mpilocal(:,:),all_intBodies_mpilocal(:,:,:)
      integer, allocatable :: all_numHeader_mpilocal(:), all_numBody_mpilocal(:)
      real(pre_obsReal), allocatable :: realBodies_mpilocal(:,:),all_realBodies_mpilocal(:,:,:)

      integer :: ierr
      integer :: get_max_rss
      integer :: numHeader_mpilocalmax,numBody_mpilocalmax
      integer :: numHeader_mpiGlobal,numBody_mpiGlobal
      integer :: numHeader_mpilocal, numBody_mpilocal
      integer :: bodyIndex_mpilocal,bodyIndex
      integer :: headerIndex_mpilocal,headerIndex
      integer :: headerIndexOffset,bodyIndexOffset
      integer :: nsize,sourcePE,nprocs_mpi,myid_mpi,procIndex
      integer :: charIndex,activeIndex,columnIndex
      !!---------------------------------------------------------------

      write(*,*) 'Entering obs_expandToMpiGlobal'
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      if(.not. obsdat%mpi_local)then
         call obs_abort('OBS_EXPANDTOMPIGLOBAL has been called, but the ' &
                        // 'obsSpaceData object is already in mpi-global state')
         return
      endif

      ! determine rank and number of mpi tasks
      call rpn_comm_size("GRID",nprocs_mpi,ierr)
      call rpn_comm_rank("GRID",myid_mpi,ierr)

      ! determine number of rows in mpiglobal arrays
      numHeader_mpiGlobal = obs_numHeader_mpiglobal(obsdat)
      numBody_mpiGlobal   = obs_numBody_mpiglobal  (obsdat)

      ! determine number of rows in mpilocal arrays
      numHeader_mpilocal = obs_numHeader(obsdat)
      numBody_mpilocal   = obs_numBody(obsdat)

      ! construct the header/body mpiglobal indices if they do not exist
      if ( .not. associated(obsdat%headerIndex_mpiglobal) ) then
         write(*,*) 'obs_expandToMpiGlobal: This object was never mpiglobal. ' // &
                    'Just gather the body and header rows in a simple way.'

         if(numHeader_mpilocal > 0) then
            ! Allocate the list of global header indices
            allocate(obsdat%headerIndex_mpiglobal(numHeader_mpilocal))
            obsdat%headerIndex_mpiglobal(:)=0
         else
            nullify(obsdat%headerIndex_mpiglobal)
            write(*,*) 'obs_expandToMpiGlobal: This mpi processor has zero headers.'
         end if

         if(numBody_mpilocal > 0) then
            ! Allocate the list of global body indices
            allocate(obsdat%bodyIndex_mpiglobal(numBody_mpilocal))
            obsdat%bodyIndex_mpiglobal(:)=0
         else
            nullify(obsdat%bodyIndex_mpiglobal)
            write(*,*) 'obs_expandToMpiGlobal: This mpi processor has zero bodies.'
         endif

         allocate( all_numHeader_mpilocal(nprocs_mpi) )
         call rpn_comm_allgather(numHeader_mpilocal,     1, "mpi_integer", &
                                 all_numHeader_mpilocal, 1, "mpi_integer", &
                                 "GRID",ierr)
         allocate( all_numBody_mpilocal(nprocs_mpi) )
         call rpn_comm_allgather(numBody_mpilocal,     1, "mpi_integer", &
                                 all_numBody_mpilocal, 1, "mpi_integer", &
                                 "GRID",ierr)

         headerIndexOffset = 0
         bodyIndexOffset = 0
         do procIndex = 1, myid_mpi
            headerIndexOffset = headerIndexOffset + all_numHeader_mpilocal(procIndex)
            bodyIndexOffset   = bodyIndexOffset   + all_numBody_mpilocal(procIndex)
         end do

         do headerIndex = 1, numHeader_mpilocal
            obsdat%headerIndex_mpiglobal(headerIndex) = headerIndex + headerIndexOffset
         end do

         do bodyIndex = 1, numBody_mpilocal
            obsdat%bodyIndex_mpiglobal(bodyIndex) = bodyIndex + bodyIndexOffset
         end do

      end if ! construct mpiglobal indices

      ! first set the mpiglobal header index value stored in the body table
      do bodyIndex_mpilocal=1,obsdat%numBody
         headerIndex_mpilocal=obs_bodyElem_i(obsdat,OBS_HIND,bodyIndex_mpilocal)
         call obs_bodySet_i(obsdat,OBS_HIND,bodyIndex_mpilocal, &
            obsdat%headerIndex_mpiglobal(headerIndex_mpilocal))
      enddo

      ! gather the lists of mpiglobal header indices on proc 0 to know where everything goes
      call rpn_comm_allreduce(obsdat%numHeader,numHeader_mpilocalmax,1,"mpi_integer","mpi_max","GRID",ierr)
      allocate(headerIndex_mpiglobal(numHeader_mpilocalmax))
      headerIndex_mpiglobal(:)=0
      do headerIndex_mpilocal=1,obsdat%numHeader
         headerIndex_mpiglobal(headerIndex_mpilocal)=obsdat%headerIndex_mpiglobal(headerIndex_mpilocal)
      enddo

      if(myid_mpi == 0) then
         allocate(all_headerIndex_mpiglobal(numHeader_mpilocalmax,0:nprocs_mpi-1))
      else
         allocate(all_headerIndex_mpiglobal(1,1))
      end if

      call rpn_comm_gather(headerIndex_mpiglobal    ,numHeader_mpilocalmax,"mpi_integer", &
                           all_headerIndex_mpiglobal,numHeader_mpilocalmax,"mpi_integer", &
                           0,"GRID",ierr)
      deallocate(headerIndex_mpiglobal)
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      ! make the header primary key mpiglobal
      allocate(headerPrimaryKey_mpilocal(numHeader_mpilocalmax))
      headerPrimaryKey_mpilocal(:) = 0
      do headerIndex_mpilocal=1,obsdat%numHeader
         headerPrimaryKey_mpilocal(headerIndex_mpilocal)=  &
            obsdat%headerPrimaryKey(headerIndex_mpilocal)
      enddo
      if(myid_mpi == 0) then
         allocate(all_headerPrimaryKey_mpilocal(numHeader_mpilocalmax,0:nprocs_mpi-1))
      else
         allocate(all_headerPrimaryKey_mpilocal(1,1))
      end if
      nsize=size(headerPrimaryKey_mpilocal)
      call rpn_comm_gather(headerPrimaryKey_mpilocal    ,nsize,"mpi_integer8", &
                           all_headerPrimaryKey_mpilocal,nsize,"mpi_integer8", &
                           0,"GRID",ierr)
      deallocate(headerPrimaryKey_mpilocal)
      
      ! make header-level integer data mpiglobal
      allocate(intHeaders_mpilocal(odc_numActiveColumn(obsdat%intHeaders),numHeader_mpilocalmax))
      intHeaders_mpilocal(:,:)=0
      do headerIndex_mpilocal=1,obsdat%numHeader
         do activeIndex=1,odc_numActiveColumn(obsdat%intHeaders)
            columnIndex=odc_columnIndexFromActiveIndex( &
                                     obsdat%intHeaders%odc_flavour, activeIndex)
            intHeaders_mpilocal(activeIndex,headerIndex_mpilocal)=  &
               obs_headElem_i(obsdat, columnIndex, headerIndex_mpilocal)
         enddo
      enddo

      if(myid_mpi == 0) then
         allocate(all_intHeaders_mpilocal(odc_numActiveColumn(obsdat%intHeaders),numHeader_mpilocalmax,0:nprocs_mpi-1))
      else
         allocate(all_intHeaders_mpilocal(1,1,1))
      end if

      nsize=size(intHeaders_mpilocal)
      call rpn_comm_gather(intHeaders_mpilocal    ,nsize,"mpi_integer", &
                           all_intHeaders_mpilocal,nsize,"mpi_integer", &
                           0,"GRID",ierr)
      deallocate(intHeaders_mpilocal)
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      ! make header-level real data mpiglobal
      allocate(realHeaders_mpilocal(odc_numActiveColumn(obsdat%realHeaders),numHeader_mpilocalmax))
      realHeaders_mpilocal(:,:)=real(0.0d0,pre_obsReal)
      do headerIndex_mpilocal=1,obsdat%numHeader
         do activeIndex=1,odc_numActiveColumn(obsdat%realHeaders)
            columnIndex=odc_columnIndexFromActiveIndex( &
                                     obsdat%realHeaders%odc_flavour, activeIndex)
            realHeaders_mpilocal(activeIndex,headerIndex_mpilocal)=  &
               obs_headElem_r(obsdat, columnIndex, headerIndex_mpilocal)
         enddo
      enddo
      
      if(myid_mpi == 0) then
         allocate(all_realHeaders_mpilocal(odc_numActiveColumn(obsdat%realHeaders),numHeader_mpilocalmax,0:nprocs_mpi-1))
      else
         allocate(all_realHeaders_mpilocal(1,1,1))
      end if

      nsize=size(realHeaders_mpilocal)
      call rpn_comm_gather(realHeaders_mpilocal    ,nsize,pre_obsMpiReal, &
                           all_realHeaders_mpilocal,nsize,pre_obsMpiReal, &
                           0,"GRID",ierr)
      deallocate(realHeaders_mpilocal)
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      ! make station-id data mpiglobal
      allocate(intStnid_mpilocal(len(obsdat%cstnid(1)),numHeader_mpilocalmax))
      intStnid_mpilocal(:,:)=0
      do headerIndex_mpilocal=1,obsdat%numHeader
         do charIndex=1,len(obsdat%cstnid(1))
            intStnid_mpilocal(charIndex,headerIndex_mpilocal)=  &
               iachar(obsdat%cstnid(headerIndex_mpilocal)(charIndex:charIndex))
         enddo
      enddo

      if(myid_mpi == 0) then
         allocate(all_intStnid_mpilocal(len(obsdat%cstnid(1)),numHeader_mpilocalmax,0:nprocs_mpi-1))
      else
         allocate(all_intStnid_mpilocal(1,1,1))
      end if

      nsize=size(intStnid_mpilocal)
      call rpn_comm_gather(intStnid_mpilocal    ,nsize,"mpi_integer", &
                           all_intStnid_mpilocal,nsize,"mpi_integer", &
                           0,"GRID",ierr)
      deallocate(intStnid_mpilocal)
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      ! make obs family data mpiglobal
      allocate(intFamily_mpilocal(len(obsdat%cfamily(1)),numHeader_mpilocalmax))
      intFamily_mpilocal(:,:)=0
      do headerIndex_mpilocal=1,obsdat%numHeader
         do charIndex=1,len(obsdat%cfamily(1))
            intFamily_mpilocal(charIndex,headerIndex_mpilocal)=  &
               iachar(obsdat%cfamily(headerIndex_mpilocal)(charIndex:charIndex))
         enddo
      enddo

      if(myid_mpi == 0) then
         allocate(all_intFamily_mpilocal(len(obsdat%cfamily(1)),numHeader_mpilocalmax,0:nprocs_mpi-1))
      else
         allocate(all_intFamily_mpilocal(1,1,1))
      end if

      nsize=size(intFamily_mpilocal)
      call rpn_comm_gather(intFamily_mpilocal    ,nsize,"mpi_integer", &
                           all_intFamily_mpilocal,nsize,"mpi_integer", &
                           0,"GRID",ierr)
      deallocate(intFamily_mpilocal)
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      ! gather the lists of mpiglobal body indices on proc 0 to know where everything goes
      call rpn_comm_allreduce(obsdat%numBody,numBody_mpilocalmax,1,"mpi_integer","mpi_max","GRID",ierr)
      allocate(bodyIndex_mpiglobal(numBody_mpilocalmax))
      bodyIndex_mpiglobal(:)=0
      do bodyIndex_mpilocal=1,obsdat%numBody
         bodyIndex_mpiglobal(bodyIndex_mpilocal)=obsdat%bodyIndex_mpiglobal(bodyIndex_mpilocal)
      enddo

      if(myid_mpi == 0) then
         allocate(all_bodyIndex_mpiglobal(numBody_mpilocalmax,0:nprocs_mpi-1))
      else
         allocate(all_bodyIndex_mpiglobal(1,1))
      end if

      call rpn_comm_gather(bodyIndex_mpiglobal    ,numBody_mpilocalmax,"mpi_integer", &
                           all_BodyIndex_mpiglobal,numBody_mpilocalmax,"mpi_integer", &
                           0,"GRID",ierr)
      deallocate(bodyIndex_mpiglobal)
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      ! make body primary key mpiglobal
      allocate(bodyPrimaryKey_mpilocal(numBody_mpilocalmax))
      bodyPrimaryKey_mpilocal(:)=0
      do bodyIndex_mpilocal=1,obsdat%numBody
         bodyPrimaryKey_mpilocal(bodyIndex_mpilocal)=  &
            obsdat%bodyPrimaryKey(bodyIndex_mpilocal)
      enddo
      if(myid_mpi == 0) then
         allocate(all_bodyPrimaryKey_mpilocal(numBody_mpilocalmax,0:nprocs_mpi-1))
      else
         allocate(all_bodyPrimaryKey_mpilocal(1,1))
      end if
      nsize=size(bodyPrimaryKey_mpilocal)
      call rpn_comm_gather(bodyPrimaryKey_mpilocal    ,nsize,"mpi_integer8", &
                           all_bodyPrimaryKey_mpilocal,nsize,"mpi_integer8", &
                           0,"GRID",ierr)
      deallocate(bodyPrimaryKey_mpilocal)

      ! make body-level integer data mpiglobal
      allocate(intBodies_mpilocal(odc_numActiveColumn(obsdat%intBodies),numBody_mpilocalmax))
      intBodies_mpilocal(:,:)=0
      do bodyIndex_mpilocal=1,obsdat%numBody
         do activeIndex=1,odc_numActiveColumn(obsdat%intBodies)
            columnIndex=odc_columnIndexFromActiveIndex( &
                                     obsdat%intBodies%odc_flavour, activeIndex)
            intBodies_mpilocal(activeIndex,bodyIndex_mpilocal)=  &
               obs_bodyElem_i(obsdat, columnIndex, bodyIndex_mpilocal)
         enddo
      enddo

      if(myid_mpi == 0) then
         allocate(all_intBodies_mpilocal(odc_numActiveColumn(obsdat%intBodies),numBody_mpilocalmax,0:nprocs_mpi-1))
      else
         allocate(all_intBodies_mpilocal(1,1,1))
      end if

      nsize=size(intBodies_mpilocal)
      call rpn_comm_gather(intBodies_mpilocal    ,nsize,"mpi_integer", &
                           all_intBodies_mpilocal,nsize,"mpi_integer", &
                           0,"GRID",ierr)
      deallocate(intBodies_mpilocal)
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      ! make body-level real data mpiglobal
      allocate(realBodies_mpilocal(odc_numActiveColumn(obsdat%realBodies),numBody_mpilocalmax))
      realBodies_mpilocal(:,:)=real(0.0d0,pre_obsReal)
      do bodyIndex_mpilocal=1,obsdat%numBody
         do activeIndex=1,odc_numActiveColumn(obsdat%realBodies)
            columnIndex=odc_columnIndexFromActiveIndex( &
                                     obsdat%realBodies%odc_flavour, activeIndex)
            realBodies_mpilocal(activeIndex,bodyIndex_mpilocal)=  &
               obs_bodyElem_r(obsdat, columnIndex, bodyIndex_mpilocal)
         enddo
      enddo

      if(myid_mpi == 0) then
         allocate(all_realBodies_mpilocal(odc_numActiveColumn(obsdat%realBodies),numBody_mpilocalmax,0:nprocs_mpi-1))
      else
         allocate(all_realBodies_mpilocal(1,1,1))
      end if

      nsize=size(realBodies_mpilocal)
      call rpn_comm_gather(realBodies_mpilocal    ,nsize,pre_obsMpiReal, &
                           all_realBodies_mpilocal,nsize,pre_obsMpiReal, &
                           0,"GRID",ierr)
      deallocate(realBodies_mpilocal)
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      ! destroy object's mpilocal data and allocate mpiglobal data
      call obs_deallocate(obsdat)

      ! Only processor 0 does any work hereafter
      if(myid_mpi == 0) then
          call obs_allocate(obsdat,numHeader_mpiGlobal,numBody_mpiGlobal)
      else
          call obs_allocate(obsdat,0,0)
      endif

      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      if(myid_mpi == 0) then

         do sourcePE=0,nprocs_mpi-1
            do headerIndex_mpilocal=1,numHeader_mpilocalmax
               ! grab the mpiglobal header index
               headerIndex=all_headerIndex_mpiglobal(headerIndex_mpilocal,sourcePE)
               if(headerIndex > 0) then

                  obsdat%headerPrimaryKey(headerIndex)= &
                         all_headerPrimaryKey_mpilocal(headerIndex_mpilocal,sourcePE)

                  do activeIndex=1,odc_numActiveColumn(obsdat%realHeaders)
                     columnIndex=odc_columnIndexFromActiveIndex(obsdat%realHeaders%odc_flavour,activeIndex)
                     obsdat%realHeaders%columns(columnIndex)%value_r(headerIndex)= &
                            all_realHeaders_mpilocal(activeIndex,headerIndex_mpilocal,sourcePE)
                  enddo

                  do activeIndex=1,odc_numActiveColumn(obsdat%intHeaders)
                     columnIndex=odc_columnIndexFromActiveIndex(obsdat%intHeaders%odc_flavour,activeIndex)
                     obsdat%intHeaders%columns(columnIndex)%value_i(headerIndex)= &
                            all_intHeaders_mpilocal(activeIndex,headerIndex_mpilocal,sourcePE)
                  enddo

                  do charIndex=1,len(obsdat%cstnid(1))
                     obsdat%cstnid(headerIndex)(charIndex:charIndex) = &
                       achar(all_intStnid_mpilocal(charIndex,headerIndex_mpilocal,sourcePE))
                  enddo

                  do charIndex=1,len(obsdat%cfamily(1))
                     obsdat%cfamily(headerIndex)(charIndex:charIndex) = &
                        achar(all_intFamily_mpilocal(charIndex,headerIndex_mpilocal,sourcePE))
                  enddo

               endif
            enddo
         enddo

         do sourcePE=0,nprocs_mpi-1
            do bodyIndex_mpilocal=1,numBody_mpilocalmax
               bodyIndex=all_bodyIndex_mpiglobal(bodyIndex_mpilocal,sourcePE)
               if(bodyIndex > 0) then

                  obsdat%bodyPrimaryKey(bodyIndex)=  &
                     all_bodyPrimaryKey_mpilocal(bodyIndex_mpilocal,sourcePE)

                  do activeIndex=1,odc_numActiveColumn(obsdat%realBodies)
                     columnIndex=odc_columnIndexFromActiveIndex(obsdat%realBodies%odc_flavour,activeIndex)
                     obsdat%realBodies%columns(columnIndex)%value_r(bodyIndex)= &
                        all_realBodies_mpilocal(activeIndex,bodyIndex_mpilocal,sourcePE)
                  enddo

                  do activeIndex=1,odc_numActiveColumn(obsdat%intBodies)
                     columnIndex=odc_columnIndexFromActiveIndex(obsdat%intBodies%odc_flavour,activeIndex)
                     obsdat%intBodies%columns(columnIndex)%value_i(bodyIndex)=  &
                        all_intBodies_mpilocal(activeIndex,bodyIndex_mpilocal,sourcePE)
                  enddo

               endif
            enddo
         enddo
 

         ! Make RLN point to global data
         do headerIndex=1,numHeader_mpiGlobal
            if(headerIndex == 1) then
               obsdat%intHeaders%columns(OBS_RLN)%value_i(headerIndex) = 1
            else
               obsdat%intHeaders%columns(OBS_RLN)%value_i(headerIndex) = &
                   obsdat%intHeaders%columns(OBS_RLN)%value_i(headerIndex-1) + &
                   obsdat%intHeaders%columns(OBS_NLV)%value_i(headerIndex-1)
            endif
         enddo

         obsdat%numBody     =  numBody_mpiGlobal
         obsdat%numHeader   =numHeader_mpiGlobal

      else

         obsdat%numBody     = 0
         obsdat%numHeader   = 0

      endif ! myid_mpi == 0

      ! deallocate the complete temporary arrays
      deallocate(all_headerIndex_mpiglobal)
      deallocate(all_bodyIndex_mpiglobal)
      deallocate(all_headerPrimaryKey_mpilocal)
      deallocate(all_bodyPrimaryKey_mpilocal)
      deallocate(all_intStnid_mpilocal)
      deallocate(all_intFamily_mpilocal)
      deallocate(all_intHeaders_mpilocal)
      deallocate(all_realHeaders_mpilocal)
      deallocate(all_intBodies_mpilocal)
      deallocate(all_realBodies_mpilocal)

      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
      obsdat%mpi_local = .false.
      write(*,*) 'Leaving obs_expandToMpiGlobal'
      return
   end subroutine obs_expandToMpiGlobal


   subroutine obs_extractObsRealBodyColumn(realBodyColumn, obsSpaceData, obsColumnIndex)
     ! :Purpose: Extract contents of a real body column into a vector.
     !
     implicit none

     ! arguments
     real(pre_obsReal), intent(out) :: realBodyColumn(:)
     type(struct_obs)     :: obsSpaceData
     integer              :: obsColumnIndex

     ! locals
     integer :: bodyIndex

     do bodyIndex = 1, obs_numBody(obsSpaceData)
       realBodyColumn(bodyIndex) = obs_bodyElem_r(obsSpaceData,obsColumnIndex,bodyIndex)
     end do

   end subroutine obs_extractObsRealBodyColumn


   subroutine obs_extractObsRealBodyColumn_r4(realBodyColumn, obsSpaceData, obsColumnIndex)
     ! :Purpose: Extract contents of a real body column into a vector.
     !
     implicit none

     ! arguments
     real(4), intent(out) :: realBodyColumn(:)
     type(struct_obs)     :: obsSpaceData
     integer              :: obsColumnIndex

     ! locals
     integer :: bodyIndex

     do bodyIndex = 1, obs_numBody(obsSpaceData)
       realBodyColumn(bodyIndex) = obs_bodyElem_r(obsSpaceData,obsColumnIndex,bodyIndex)
     end do

   end subroutine obs_extractObsRealBodyColumn_r4


   subroutine obs_extractObsIntBodyColumn(intBodyColumn, obsSpaceData, obsColumnIndex)
     ! :Purpose: Extract contents of an integer body column into a vector.
     !
     implicit none

     ! arguments
     integer, intent(out) :: intBodyColumn(:)
     type(struct_obs)     :: obsSpaceData
     integer              :: obsColumnIndex

     ! locals
     integer :: bodyIndex

     do bodyIndex = 1, obs_numBody(obsSpaceData)
       intBodyColumn(bodyIndex) = obs_bodyElem_i(obsSpaceData,obsColumnIndex,bodyIndex)
     end do

   end subroutine obs_extractObsIntBodyColumn


   subroutine obs_extractObsRealHeaderColumn(realHeaderColumn, obsSpaceData, obsColumnIndex)
     ! :Purpose: Extract contents of a real header column into a vector. Note
     !           that the output can be either in the form of a vector with
     !           length equal to the number of rows in the body OR header table.
     !
     implicit none

     ! arguments
     real(pre_obsReal), intent(out) :: realHeaderColumn(:)
     type(struct_obs)     :: obsSpaceData
     integer              :: obsColumnIndex

     ! locals
     integer :: bodyIndex, headerIndex

     if (size(realHeaderColumn) == obs_numBody(obsSpaceData)) then
       do bodyIndex = 1, obs_numBody(obsSpaceData)
         headerIndex = obs_bodyElem_i(obsSpaceData,OBS_HIND,bodyIndex)
         realHeaderColumn(bodyIndex) = obs_headElem_r(obsSpaceData,obsColumnIndex,headerIndex)
       end do
     else if(size(realHeaderColumn) == obs_numHeader(obsSpaceData)) then
       do headerIndex = 1, obs_numHeader(obsSpaceData)
         realHeaderColumn(headerIndex) = obs_headElem_r(obsSpaceData,obsColumnIndex,headerIndex)
       end do
     else
       call obs_abort('extractObsRealHeaderColumn: size of output vector invalid')
     end if

   end subroutine obs_extractObsRealHeaderColumn


   subroutine obs_extractObsIntHeaderColumn(intHeaderColumn, obsSpaceData, obsColumnIndex)
     ! :Purpose: Extract contents of an integer header column into a vector. Note
     !           that the output can be either in the form of a vector with
     !           length equal to the number of rows in the body OR header table.
     !
     implicit none

     ! arguments
     integer, intent(out) :: intHeaderColumn(:)
     type(struct_obs)     :: obsSpaceData
     integer              :: obsColumnIndex

     ! locals
     integer :: bodyIndex, headerIndex

     if (size(intHeaderColumn) == obs_numBody(obsSpaceData)) then
       do bodyIndex = 1, obs_numBody(obsSpaceData)
         headerIndex = obs_bodyElem_i(obsSpaceData,OBS_HIND,bodyIndex)
         intHeaderColumn(bodyIndex) = obs_headElem_i(obsSpaceData,obsColumnIndex,headerIndex)
       end do
     else if(size(intHeaderColumn) == obs_numHeader(obsSpaceData)) then
       do headerIndex = 1, obs_numHeader(obsSpaceData)
         intHeaderColumn(headerIndex) = obs_headElem_i(obsSpaceData,obsColumnIndex,headerIndex)
       end do
     else
       call obs_abort('extractObsIntHeaderColumn: size of output vector invalid')
     end if

   end subroutine obs_extractObsIntHeaderColumn


   subroutine obs_finalize(obsdat)
      !
      ! :Purpose: De-allocate memory and clean up the object.
      !      De-allocate object arrays, and perform any other clean-up that is
      !      necessary before object deletion.
      !
      implicit none
      type(struct_obs), intent(inout) :: obsdat

      call obs_deallocate(obsdat)
   end subroutine obs_finalize


   subroutine obs_generate_header(obsdat, ilat, ilon, ialt, inbon, instrum, &
                                  satzen, isat, itech, nvtyp, ity, idate, &
                                  itime, clstnid, imask, satazim, sunza, clfr)
      !
      ! :Note: Output:
      !
      !           - obsdat%realHeaders%columns(OBS_LON,) - in degrees
      !           - obsdat%realHeaders%columns(OBS_LAT,) - in degrees, equator at 0 degrees
      !           - obsdat%realHeaders%columns(OBS_ALT,) - in metres, with no offset
      !
      use MathPhysConstants_mod
      implicit none

      type(struct_obs), intent(inout) :: obsdat
      integer, intent(in) :: ilat
      integer, intent(in) :: ilon
      integer, intent(in) :: ialt
      integer, intent(in) :: inbon
      integer, intent(in) :: instrum
      integer, intent(in) :: isat
      integer, intent(in) :: itech
      integer, intent(in) :: nvtyp
      integer, intent(in) :: ity
      integer, intent(in) :: idate
      integer, intent(in) :: itime
      integer, intent(in) :: imask
      real(kind=pre_obsReal) :: clfr
      real(kind=pre_obsReal) :: sunza
      real(kind=pre_obsReal) :: satzen
      real(kind=pre_obsReal) :: satazim
      character(len=9), intent(in) :: clstnid
      real(kind=8) :: torad

      torad=MPC_RADIANS_PER_DEGREE_R8

      !
      !     IF VALID DATA WERE FOUND GENERATE THE OBSDAT HEADER
      !      AND INCREMENT OBSDAT%numHeader
      !
      ! PLH          if  ( obsdat%numHeader < nmxobs) then
      if ( obsdat%numHeader < obsdat%numHeader_max) then
         obsdat%numHeader=obsdat%numHeader + 1

         obsdat%realHeaders%columns(OBS_LON)%value_r(obsdat%numHeader) = real(ilon) *0.01
         obsdat%realHeaders%columns(OBS_LAT)%value_r(obsdat%numHeader) = real(ilat) *0.01-90.0
         ! PLH ADDED OBS_BX OBS_BY OBS_BZ
         obsdat%realHeaders%columns(OBS_BX)%value_r(obsdat%numHeader)=0.0
         obsdat%realHeaders%columns(OBS_BY)%value_r(obsdat%numHeader)=0.0
         obsdat%realHeaders%columns(OBS_BZ)%value_r(obsdat%numHeader)=0.0

         obsdat%realHeaders%columns(OBS_ALT)%value_r(obsdat%numHeader) = real(ialt)
         ! PLH       obsdat%realHeaders%columns(ncmtlo)%value_r(obsdat%numHeader) = (real(ilon)*0.01)*ztorad
         ! PLH       obsdat%realHeaders%columns(ncmtla)%value_r(obsdat%numHeader) = (real(ilat)*0.01-90.)*ztorad
         call obs_headSet_i(obsdat, OBS_NLV, obsdat%numHeader, inbon)
         !       print*,'NOBTOTAL=',obsdat%numHeader

         if ( obsdat%numHeader == 1) then
            ! This is the first entry into the obsdat
            call obs_headSet_i(obsdat, OBS_RLN, 1, 1)
         else
            call obs_headSet_i(obsdat, OBS_RLN, obsdat%numHeader, &
                obs_headElem_i(obsdat, OBS_RLN, obsdat%numHeader-1) &
              + obs_headElem_i(obsdat, OBS_NLV, obsdat%numHeader-1))
         endif
         !
         !          REMAINDER OF HEADER
         !
         call obs_headSet_i(obsdat, OBS_ONM, obsdat%numHeader, obsdat%numHeader)
         call obs_headSet_i(obsdat, OBS_INS, obsdat%numHeader, instrum )
         call obs_headSet_i(obsdat, OBS_SAT, obsdat%numHeader, isat)
         call obs_headSet_i(obsdat, OBS_TEC, obsdat%numHeader, itech)
         call obs_headSet_i(obsdat, OBS_OTP, obsdat%numHeader, nvtyp)
         call obs_headSet_i(obsdat, OBS_ITY, obsdat%numHeader, ity)
         call obs_headSet_i(obsdat, OBS_DAT, obsdat%numHeader, idate)
         call obs_headSet_i(obsdat, OBS_ETM, obsdat%numHeader, itime)
         obsdat%cstnid(obsdat%numHeader)         = clstnid
         ! PLH       call obs_headSet_i(obsdat, ncmoec, obsdat%numHeader, 999)
         call obs_headSet_i(obsdat, OBS_STYP,obsdat%numHeader, imask)
         call obs_headSet_r(obsdat, OBS_AZA, obsdat%numHeader, satazim)
         call obs_headSet_r(obsdat, OBS_SUN, obsdat%numHeader, sunza)
         call obs_headSet_r(obsdat, OBS_SZA, obsdat%numHeader, satzen)
         call obs_headSet_r(obsdat, OBS_CLF, obsdat%numHeader, clfr)
         ! PLH       call obs_headSet_i(obsdat, ncmst1, obsdat%numHeader, iflgs)

      else

         call obs_abort('obs_generate_header: numHeader reached numHeader_max')

      endif

   end subroutine obs_generate_header


   function obs_getBodyIndex_depot(obsdat) result(row_index)
      !
      ! :Purpose:
      !      Return the next element from the current body list
      !
      implicit none
      integer :: row_index
      type(struct_obs), intent(inout) :: obsdat

      row_index = ild_get_next_index(obsdat%body_index_list_depot)
   end function obs_getBodyIndex_depot


   function obs_getBodyIndex_private(private_list) result(row_index)
      !
      ! :Purpose:
      !      Return the next element from the supplied private body list
      !
      implicit none
      integer :: row_index
      type(struct_index_list), pointer, intent(inout) :: private_list

      row_index = ild_get_next_index(private_list)
   end function obs_getBodyIndex_private


   function obs_getFamily(obsdat,headerIndex_in,bodyIndex)
      !
      ! :Purpose:
      !      Return the family for the indicated header, or else for the
      !      indicated body.
      !
      implicit none

      character(len=2)             :: obs_getFamily
      type(struct_obs), intent(in) :: obsdat
      integer,optional, intent(in) :: headerIndex_in,bodyIndex

      integer          :: headerIndex

      if(present(headerIndex_in)) then
         headerIndex=headerIndex_in
      elseif(present(bodyIndex)) then
         headerIndex=obs_bodyElem_i(obsdat,OBS_HIND,bodyIndex)
      else
         call obs_abort('OBS_GETFAMILY: Header or Body index must be specified!')
         return
      endif

      obs_getFamily=obsdat%cfamily(headerIndex)

   end function obs_getFamily


   function obs_getNclassAvhrr() 
      !
      ! :Purpose:
      !      to get the number of  AVHRR radiance classes
      !
     implicit none

     integer :: obs_getNclassAvhrr

     obs_getNclassAvhrr = ( OBS_CF7 - OBS_CF1 + 1 )

   end function obs_getNclassAvhrr

   function obs_getNchanAvhrr() 
      !
      ! :Purpose:
      !      to get the number of AVHRR channels
      !
     implicit none

     integer :: obs_getNchanAvhrr

     obs_getNchanAvhrr = ( OBS_M1C6 - OBS_M1C1 + 1 )

   end function obs_getNchanAvhrr

   function obs_getHeaderIndex(obsdat) result(row_index)
      !
      ! :Purpose:
      !      Return the next element from the current header list.
      !
      implicit none
      integer :: row_index
      type(struct_obs), intent(inout) :: obsdat

      row_index = ild_get_next_index(obsdat%header_index_list_depot)
   end function obs_getHeaderIndex


   function obs_headElem_i(obsdat,column_index,row_index) result(value_i)
      !
      ! :Purpose: Get an integer-valued header observation-data element.
      !      To control access to the observation object.  Returns the (integer)
      !      value of the row_index'th ObsData element with the indicated column
      !      index from the "header".
      !
      implicit none
      integer                       :: value_i
      type(struct_obs), intent(in)  :: obsdat
      integer         , intent(in)  :: column_index
      integer         , intent(in)  :: row_index

      real(pre_obsReal) :: value_r         ! not used

      call odc_columnElem(obsdat%intHeaders, column_index, row_index, &
                          value_i, value_r)
   end function obs_headElem_i


   function obs_headElem_r(obsdat,column_index,row_index) result(value_r)
      !
      ! :Purpose: Get a real-valued header observation-data element.
      !      To control access to the observation object.  Returns the (real)
      !      value of the row_index'th ObsData element with the indicated column
      !      index from the "header".
      !
      implicit none
      real(pre_obsReal)             :: value_r
      type(struct_obs), intent(in)  :: obsdat
      integer         , intent(in)  :: column_index
      integer         , intent(in)  :: row_index

      integer :: value_i                ! unused

      call odc_columnElem(obsdat%realHeaders, column_index, row_index, &
                          value_i, value_r)
   end function obs_headElem_r


   function obs_headPrimaryKey(obsdat,row_index) result(primaryKey)
      !
      ! :Purpose: Get the header primary key value.
      !
      implicit none
      integer(8)                    :: primaryKey
      type(struct_obs), intent(in)  :: obsdat
      integer         , intent(in)  :: row_index

      primaryKey = obsdat%headerPrimaryKey(row_index)

   end function obs_headPrimaryKey


   function obs_headerIndex_mpiglobal(obsdat,row_index) result(value)
      !
      ! :Purpose: Get the mpiglobal header row index.
      !      To control access to the mpiglobal row_index into the "header".
      !
      implicit none
      integer value
      type(struct_obs), intent(in)  :: obsdat
      integer         , intent(in)  :: row_index

      value=obsdat%headerIndex_mpiglobal(row_index)
   end function obs_headerIndex_mpiglobal


   subroutine obs_headSet_i(obsdat, column_index, row_index, value_i)
      !
      ! :Purpose: set an integer-valued header observation-data element.
      !      To control access to the observation object.  Sets the (integer)
      !      value of the row_index'th ObsData element with the indicated column
      !      index from the "header".
      !
      implicit none
      type(struct_obs), intent(inout)  :: obsdat
      integer         , intent(in)     :: column_index
      integer         , intent(in)     :: row_index
      integer         , intent(in)     :: value_i

      call odc_columnSet(obsdat%intHeaders, column_index, row_index, &
                         value_i, NULL_COLUMN_VALUE_R, &
                         obsdat%numHeader, obsdat%numHeader_max)
   end subroutine obs_headSet_i


   subroutine obs_headSet_r4(obsdat, column_index, row_index, value_r4)
      !
      ! :Purpose: set a real header value in the observation object.
      !      To control access to the observation object.  Sets the (real)
      !      value of the row_index'th ObsData element with the indicated column
      !      index from the "header".
      !
      implicit none
      type(struct_obs), intent(inout)  :: obsdat
      integer         , intent(in)     :: column_index
      integer         , intent(in)     :: row_index
      real(4)         , intent(in)     :: value_r4

      real(pre_obsReal)                :: value_r

      value_r = value_r4

      call odc_columnSet(obsdat%realHeaders, column_index, row_index, &
                         NULL_COLUMN_VALUE_I, value_r, &
                         obsdat%numHeader, obsdat%numHeader_max)
   end subroutine obs_headSet_r4


   subroutine obs_headSet_r8(obsdat, column_index, row_index, value_r8)
      !
      ! :Purpose: set a real header value in the observation object.
      !      To control access to the observation object.  Sets the (real)
      !      value of the row_index'th ObsData element with the indicated column
      !      index from the "header".
      !
      implicit none
      type(struct_obs), intent(inout)  :: obsdat
      integer         , intent(in)     :: column_index
      integer         , intent(in)     :: row_index
      real(8)         , intent(in)     :: value_r8

      real(pre_obsReal)                :: value_r

      value_r = value_r8

      call odc_columnSet(obsdat%realHeaders, column_index, row_index, &
                         NULL_COLUMN_VALUE_I, value_r, &
                         obsdat%numHeader, obsdat%numHeader_max)
   end subroutine obs_headSet_r8


   subroutine obs_initialize(obsdat, numHeader_max, numBody_max, mpi_local, &
                             silent)
      !
      ! :Purpose: Set an observation-data module to a known state.
      !      Initialize object variables, and allocate arrays according to the
      !      parameters, header_max and body_max.
      !
      implicit none
                                        ! instance of obsSpaceData
      type(struct_obs),intent(inout):: obsdat !inout allows detection of 2nd call
                                        ! number of header elements allocated
      integer, optional, intent(in) :: numHeader_max
                                        ! total no. of body elements allocated
      integer, optional, intent(in) :: numBody_max
      logical, optional, intent(in) :: mpi_local
      logical, optional, intent(in) :: silent

      logical :: silent_
      integer :: nulnam,fnom,fclos,ierr
      integer :: nmxobs,ndatamx
      namelist /namdimo/nmxobs,ndatamx
      character(len=120) :: message

      if(.not. obs_class_initialized) then
         call obs_abort('obs_class_initialize must be called before ' // &
                        'obs_initialize')
         return
      end if

      !
      ! INITIALIZE ALL OBJECT VARIABLES
      !
      nullify(obsdat%headerPrimaryKey)
      nullify(obsdat%bodyPrimaryKey)
      nullify(obsdat%cstnid)
      nullify(obsdat%cfamily)

      obsdat%numHeader     = 0
      obsdat%numHeader_max = 0
      obsdat%numBody       = 0
      obsdat%numBody_max   = 0

      if(present(mpi_local)) then
         obsdat%mpi_local  = mpi_local
      else
         obsdat%mpi_local  = .false.
      end if

      if(present(silent)) then
         silent_  = silent
      else
         silent_  = .false.
      end if

      nullify(obsdat%headerIndex_mpiglobal)
      nullify(obsdat%bodyIndex_mpiglobal)

      !
      ! DETERMINE THE ARRAY DIMENSIONS
      !
      if(present(numHeader_max)) then
         ! numBody_max is necessarily also present
         nmxobs = numHeader_max
         ndatamx = numBody_max

      else
         ! Initialize with bad values
         nmxobs=0
         ndatamx=0

         ! Open the file, flnml
         nulnam=0
         ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
         if(ierr < 0) then
            write(message,*)'Failed to open flnml to obtain nmxobs and ndatamx:'&
                            // '  ierr=', ierr
            call obs_abort(message); return
         end if

         ! Read the dimensions from a namelist
         read(nulnam,nml=namdimo,iostat=ierr)
         if(ierr /= 0) call obs_abort('obs_initialize: Error reading namelist')
         write(*,nml=namdimo)
         ierr=fclos(nulnam)
         ! Verify that the namelist contained values
         if(nmxobs <= 0 .or. ndatamx <= 0) then
            write(message,*)'From file, flnml, positive values were not ' &
                            // 'obtained for nmxobs or ndatamx:  ',nmxobs,ndatamx
            call obs_abort(message); return
         end if
      end if ! present(numHeader_max)

      !
      ! ALLOCATE
      !
      call obs_allocate(obsdat, nmxobs, ndatamx, silent_)

      return
   end subroutine obs_initialize


   subroutine obs_mpiDistributeIndices(obsdat)
      !
      ! :Purpose:
      !  Compute headerIndex_mpiglobal and bodyIndex_mpiglobal: 
      !  this determines how obs are distributed over MPI processes
      !  and is needed for converting from mpiglobal to mpilocal and vice versa.
      !  The header indices are distributed following the chosen strategy,
      !  currently either "round robin" or by latitude bands.
      !
      ! :Note: this subroutine is called before converting from mpiglobal to 
      !        mpilocal
      !
      ! :Comments:  In principle this method could have obtained
      !           my_mpi_id by use'ing the module, mpi.  However, it queries
      !           rpn_comm for itself because the mpi module belongs to the
      !           3dvar code, whereas the present module is shared code.
      !
      implicit none
      type(struct_obs), intent(inout) :: obsdat

      integer :: headerIndex_mpiglobal,headerIndex_mpilocal
      integer ::   bodyIndex_mpiglobal,  bodyIndex_mpilocal
      integer :: numHeader_mpiLocal,numBody_mpiLocal,idata,idataend
      integer :: my_mpi_id, my_mpi_idx_dummy, my_mpi_idy_dummy

      write(*,*) '-------- Start obs_mpiDistributeIndices ---------'

      if(obsdat%mpi_local) then
         call obs_abort( &
                      'obs_mpiDistributeIndices: data already mpi-local, Abort!')
         return
      end if

      call rpn_comm_mype(my_mpi_id, my_mpi_idx_dummy, my_mpi_idy_dummy)

      ! Count number of headers and bodies for each processor
      numHeader_mpiLocal=0
      numBody_mpiLocal=0
      do headerIndex_mpiglobal=1,obsdat%numHeader
         if ( my_mpi_id == obs_headElem_i( obsdat, OBS_IP, headerIndex_mpiglobal )) then
            numHeader_mpiLocal = numHeader_mpiLocal + 1
            numBody_mpilocal = numBody_mpilocal &
                            +obs_headElem_i( obsdat, OBS_NLV, headerIndex_mpiglobal )
         end if
      enddo
      write(*,*) 'obs_mpidistributeindices: numHeader_mpiLocal,global=',  &
                 numHeader_mpiLocal,obsdat%numHeader
      write(*,*) 'obs_mpidistributeindices: numBody_mpiLocal,global=',  &
                 numBody_mpiLocal,obsdat%numBody

      if ( numHeader_mpilocal > 0 ) then
         ! Allocate the list of global header indices
         allocate(obsdat%headerIndex_mpiglobal(numHeader_mpilocal))
         obsdat%headerIndex_mpiglobal(:)=0
      else
         nullify(obsdat%headerIndex_mpiglobal)
         write(*,*) 'This mpi processor has zero headers.'
      end if

      if ( numBody_mpilocal > 0 ) then
         ! Allocate the list of global body indices
         allocate(obsdat%bodyIndex_mpiglobal(numBody_mpilocal))
         obsdat%bodyIndex_mpiglobal(:)=0
      else
         nullify(obsdat%bodyIndex_mpiglobal)
         write(*,*) 'This mpi processor has zero bodies to treat'
      endif

      ! determine the list of header indices
      headerIndex_mpilocal = 0
      do headerIndex_mpiglobal=1,obsdat%numHeader
         if ( my_mpi_id == obs_headElem_i( obsdat, OBS_IP, headerIndex_mpiglobal )) then
            headerIndex_mpilocal=headerIndex_mpilocal+1
            obsdat%headerIndex_mpiglobal(headerIndex_mpilocal) &
                                                           =headerIndex_mpiglobal
         endif
      enddo

      ! determine the corresponding list of body indices
      bodyIndex_mpilocal=0
      do headerIndex_mpilocal=1,numHeader_mpilocal
         headerIndex_mpiglobal=obsdat%headerIndex_mpiglobal(headerIndex_mpilocal)
         idata= obs_headElem_i(obsdat, OBS_RLN, headerIndex_mpiglobal)
         idataend = obs_headElem_i(obsdat, OBS_NLV, headerIndex_mpiglobal) &
                  + idata -1 
         do bodyIndex_mpiglobal=idata,idataend 
            bodyIndex_mpilocal=bodyIndex_mpilocal+1 
            obsdat%bodyIndex_mpiglobal(bodyIndex_mpilocal) = bodyIndex_mpiglobal
         enddo
      enddo

      write(*,*) '-------- END OF obs_mpiDistributeIndices ---------'
      write(*,*) ' '
   end subroutine obs_mpiDistributeIndices


   logical function obs_mpiLocal(obsdat)
      !
      ! :Purpose: returns true if the object contains only data that are
      !          needed by the current mpi PE; false if it contains all
      !          data.
      !      To provide the state of the internal variable, mpiLocal.  This
      !      method exists primarily to facilitate unit tests on this module.
      !
      implicit none
      type(struct_obs) , intent(in)  :: obsdat
      obs_mpiLocal=obsdat%mpi_local
   end function obs_mpiLocal


   integer function obs_numBody(obsdat)
      !
      ! :Purpose: returns the number of mpi-local bodies recorded.
      !      To provide the number of bodies that are currently recorded in the
      !      mpi-local observation-data object.
      !
      implicit none
      type(struct_obs) , intent(in)  :: obsdat
      obs_numBody=obsdat%numBody
   end function obs_numBody


   integer function obs_numBody_max(obsdat)
      !
      ! :Purpose: returns the dimensioned mpi-local number of bodies.
      !      To provide the dimension for the number of bodies in the mpi-local
      !      observation-data object.
      !
      implicit none

      type(struct_obs) , intent(in)  :: obsdat

      obs_numBody_max=obsdat%numBody_max
   end function obs_numBody_max


   integer function obs_numBody_mpiglobal(obsdat)
      !
      ! :Purpose: returns the number of bodies recorded in the
      !           entire mpi-global obs object.
      !      To provide the number of bodies that are currently recorded in the
      !      entire mpi-global observation-data object.
      !
      implicit none

      type(struct_obs), intent(in)  :: obsdat

      integer :: numBody_mpiGlobal, sizedata, ierr

      if(obsdat%mpi_local)then
         sizedata=1
         call rpn_comm_allreduce(obsdat%numBody,numBody_mpiGlobal,sizedata, &
                                 "mpi_integer","mpi_sum","GRID",ierr)
         obs_numBody_mpiglobal = numBody_mpiGlobal
      else
         obs_numBody_mpiglobal = obsdat%numBody
      end if

   end function obs_numBody_mpiglobal


   integer function obs_numHeader(obsdat)
      !
      ! :Purpose: returns the number of mpi-local headers recorded.
      !      To provide the number of headers that are currently recorded in the
      !      observation-data object.
      !
      implicit none

      type(struct_obs) , intent(in)  :: obsdat

      obs_numHeader=obsdat%numHeader
   end function obs_numHeader


   integer function obs_numHeader_max(obsdat)
      !
      ! :Purpose: returns the dimensioned mpi-local number of headers.
      !      To provide the dimension for the number of headers in the mpi-local
      !      observation-data object.
      !
      implicit none
      type(struct_obs) , intent(in)  :: obsdat
      obs_numHeader_max=obsdat%numHeader_max
   end function obs_numHeader_max


   integer function obs_numHeader_mpiglobal(obsdat)
      !
      ! :Purpose: returns the number of headers recorded in
      !           the entire mpi-global obs object.
      !      To provide the number of headers that are currently recorded in the
      !      entire mpi-global observation-data object.
      !
      implicit none

      type(struct_obs) , intent(in)  :: obsdat

      integer :: numHeader_mpiGlobal, sizedata, ierr

      if(obsdat%mpi_local)then
         sizedata=1
         call rpn_comm_allreduce(obsdat%numHeader,numHeader_mpiGlobal,sizedata, &
                                 "mpi_integer","mpi_sum","GRID",ierr)
         obs_numHeader_mpiglobal = numHeader_mpiGlobal
      else
         obs_numHeader_mpiglobal = obsdat%numHeader
      end if
   end function obs_numHeader_mpiglobal


   subroutine obs_order(obsdat)
      !
      ! :Purpose:
      !      Put an obsdat file into the order required for the sequential
      !      assimilation. Note that it is known, as a by-product of the
      !      algorithm that was used  to determine the pass and the region for
      !      each station, at what exact location (information  in OBS_ONM) each
      !      station has to be. The algorithm requires the exchange of at most
      !      mxstn headers.  A faster algorithm likely exists.
      !
      implicit none

      type (struct_obs), intent(inout) :: obsdat

      integer :: hdr,jk
      logical :: sorted
      integer :: bodyElem, first, last

      do hdr = 1, obsdat%numHeader
         sorted = .false.
         do while(.not.sorted)
            jk = obs_headElem_i(obsdat, OBS_ONM, hdr)
            if ( jk == hdr ) then
               sorted=.true.
            else
               call obs_exchange_stations( obsdat, jk, hdr )
            endif
         end do
      enddo

      do hdr = 1, obsdat%numHeader
         ! Make the body members of header(hdr) point to the new row_index of hdr.
         first=        obs_headElem_i(obsdat, OBS_RLN, hdr)
         last =first + obs_headElem_i(obsdat, OBS_NLV, hdr) -1
         do bodyElem=first,last
            call obs_bodySet_i(obsdat, OBS_HIND, bodyElem, hdr)
         end do
      enddo

      return

   contains


      subroutine obs_exchange_stations( obsdat, j, k ) 
         !
         ! :Purpose: exchange the headers of stations j and k 
         !
         implicit none

         type (struct_obs), intent(inout) :: obsdat
         integer          , intent(in)    :: j
         integer          , intent(in)    :: k

         real(pre_obsReal):: rdum 
         integer          :: idum
         integer          :: column_index
         character(12)    :: cdum

         do column_index=NHDR_INT_BEG, NHDR_INT_END
            if(obsdat%intHeaders%odc_flavour%columnActive(column_index)) then
               idum=obs_headElem_i(obsdat, column_index, j)
               call obs_headSet_i(obsdat, column_index, j, &
                   obs_headElem_i(obsdat, column_index, k))
               call obs_headSet_i(obsdat, column_index, k, idum)
            endif
         enddo

         do column_index=NHDR_REAL_BEG,NHDR_REAL_END
            if(obsdat%realHeaders%odc_flavour%columnActive(column_index)) then
               rdum=obs_headElem_r(obsdat, column_index, j)
               call obs_headSet_r(obsdat, column_index, j, &
                   obs_headElem_r(obsdat, column_index, k))
               call obs_headSet_r(obsdat, column_index, k, rdum)
            endif
         enddo

         cdum=obsdat%cstnid(j)
         obsdat%cstnid(j)=obsdat%cstnid(k)
         obsdat%cstnid(k)=cdum

         return
      end subroutine obs_exchange_stations
   end subroutine obs_order


   subroutine obs_print( obsdat, nobsout )
      !
      ! :Purpose: print the contents of the obsdat to an ASCII file
      !
      ! :Arguments:
      !           :obsdat:  obsSpaceData object
      !           :nobsout: unit used for printing
      !
      implicit none

      type(struct_obs), intent(inout) :: obsdat
      integer,         intent(in)    :: nobsout

      integer :: jo

      do jo=1,obsdat%numHeader
         call obs_enkf_prnthdr(obsdat,jo,nobsout)
         call obs_enkf_prntbdy(obsdat,jo,nobsout)
      enddo

      return
   end subroutine obs_print


   subroutine obs_prnt_csv( obsdat, nhdrsql, nbdysql )
      !
      ! :Purpose: print the contents of the obsdat to csv (comma separated
      !           values) files
      !
      ! :Arguments:
      !           :obsdat:  obsSpaceData object
      !           :nhdrsql: unit used for printing header
      !           :nbdysql: unit used for printing body
      !
      implicit none

      type (struct_obs), intent(inout) :: obsdat
      integer          , intent(in)    :: nhdrsql
      integer          , intent(in)    :: nbdysql

      integer :: jo

      do jo=1,obsdat%numHeader
         call obs_tosqlhdr(obsdat,jo,nhdrsql)
         call obs_tosqlbdy(obsdat,jo,nbdysql)
      enddo

      return
   end subroutine obs_prnt_csv


   subroutine obs_prntbdy( obsdat , index_header, unitout_opt )
      !
      ! :Purpose: Print all data records associated with an observation
      !
      ! :Arguments:
      !           :obsdat:  obsSpaceData object
      !           :index_header: index of the group of observations to be printed
      !           :unitout: unit number on which to print
      !
      implicit none

      type(struct_obs), intent(in)           :: obsdat
      integer         , intent(in)           :: index_header
      integer         , intent(in), optional :: unitout_opt ! variable output unit facilitates unit testing

      integer :: unitout_

      integer :: ipnt, idata, idata2, jdata, ivco
      character(len=13) :: ccordtyp(4)

      integer :: obsVNM, obsFLG, obsASS
      real(pre_obsReal) :: obsPPP, obsVAR, obsOMP, obsOMA, obsOER, obsHPHT, obsOMPE

      if ( present( unitout_opt ) ) then
        unitout_ = unitout_opt
      else
        unitout_ = 6
      end if

      ccordtyp(1) = 'HEIGHT      :'
      ccordtyp(2) = 'PRESSURE    :'
      ccordtyp(3) = 'CHANNEL NUM :'
      ccordtyp(4) = 'VCO UNDEFINED'
      !
      ! 1. General information
      !
      ipnt  = obs_headElem_i(obsdat,OBS_RLN,index_header)
      idata = obs_headElem_i(obsdat,OBS_NLV,index_header)

      if ( idata == 1) then
        write ( unitout_, fmt=9101 ) idata, index_header, NBDY_INT_SIZE+NBDY_REAL_SIZE
      else
        write ( unitout_, fmt=9100 ) idata, index_header, NBDY_INT_SIZE+NBDY_REAL_SIZE
      end if

9100  format(4x,'THERE ARE ', &
         i3,1x,'DATA IN OBSERVATION RECORD NO.' &
         ,1x,i6,4x,'DATA RECORD''S LENGTH:',i6)
9101  format(4x,'THERE IS ', &
         i3,1x,'DATUM IN OBSERVATION RECORD NO.' &
         ,1x,i6,4x,'DATA RECORD''S LENGTH:',i6)
      !
      ! 2. Print all data records
      !
      do jdata = ipnt, ipnt + idata - 1
        idata2 = jdata -ipnt + 1
        if ( obs_bodyElem_i( obsdat, OBS_ASS, jdata ) >= 0) then
          ivco = obs_bodyElem_i( obsdat, OBS_VCO, jdata )
          if ( ivco < 1 .or. ivco > 3 ) ivco = 4
          obsVNM = obs_bodyElem_i( obsdat, OBS_VNM , jdata )
          obsPPP = obs_bodyElem_r( obsdat, OBS_PPP , jdata )
          obsVAR = obs_bodyElem_r( obsdat, OBS_VAR , jdata )
          obsOMP = obs_bodyElem_r( obsdat, OBS_OMP , jdata )
          obsOMA = obs_bodyElem_r( obsdat, OBS_OMA , jdata )
          obsOER = obs_bodyElem_r( obsdat, OBS_OER , jdata )
          obsHPHT= obs_bodyElem_r( obsdat, OBS_HPHT, jdata )
          obsOMPE= obs_bodyElem_r( obsdat, OBS_OMPE, jdata )
          obsFLG = obs_bodyElem_i( obsdat, OBS_FLG , jdata )
          obsASS = obs_bodyElem_i( obsdat, OBS_ASS , jdata )
          write ( unitout_, fmt=9201 ) idata2, &
               obsVNM, ccordtyp(ivco), &
               obsPPP,obsVAR,obsOMP,obsOMA, &
               obsOER,obsHPHT,obsOMPE,obsFLG,obsASS
        end if
      end do

9201  format(4x,'DATA NO.',i6,/,10x &
         ,'VARIABLE NO.:',i6,4x,a13,g13.6,4x &
         ,/,10x &
         ,'OBSERVED VALUE:',g23.16,5x,'OBSERVED - BACKGROUND VALUES:' &
         ,g23.16,4x &
         ,/,10x &
         ,'OBSERVED - ANALYZED VALUES:',g13.6,4x &
         ,/,10x &
         ,'ERROR STANDARD DEVIATIONS FOR' &
         ,/,20x &
         ,'OBSERVATION:',g13.6,4x &
         ,/,20x &
         ,'FIRST-GUESS:',g13.6,4x &
         ,/,20x &
         ,'OBS minus F-G:',g13.6,4x &
         ,/,10x &
         ,'BURP FLAGS:',i6,4x,'OBS. ASSIMILATED (1-->YES;0-->NO):',i3)

      return

   end subroutine obs_prntbdy


   subroutine obs_prnthdr( obsdata, index_hd, unitout_opt )
      !
      ! :Purpose: Printing of the header of an observation record
      ! :Arguments:
      !           :obsdata:  obsSpaceData object
      !           :index_hd: index of the header to be printed
      !           :unitout_opt: unit number on which to print
      !

      implicit none


      type(struct_obs), intent(in)           :: obsdata
      integer         , intent(in)           :: index_hd
      integer         , intent(in), optional :: unitout_opt ! variable output unit facilitates unit testing

      integer :: unitout_
      integer :: obsRLN, obsONM, obsINS, obsOTP, obsITY
      integer :: obsDAT, obsETM, obsNLV, obsST1, obsSTYP, obsTTYP
      real(pre_obsReal) :: obsLON, obsLAT, obsALT
      character(len=12) :: stnid

      if ( present( unitout_opt ) ) then
         unitout_ = unitout_opt
      else
         unitout_ = 6
      end if

      !
      ! 1. General information
      !
      write(unitout_,fmt=9100)index_hd, NHDR_INT_SIZE + NHDR_REAL_SIZE
9100  format(//,10x,'-- OBSERVATION RECORD NO.' &
         ,1x,i6,3x,'HEADER''S LENGTH:',i6)
      !
      ! 2. PRINT HEADER'S CONTENT
      !
      obsRLN = obs_headElem_i(obsdata,OBS_RLN,index_hd)
      obsONM = obs_headElem_i(obsdata,OBS_ONM,index_hd)
      obsINS = obs_headElem_i(obsdata,OBS_INS,index_hd)
      obsOTP = obs_headElem_i(obsdata,OBS_OTP,index_hd)
      obsITY = obs_headElem_i(obsdata,OBS_ITY,index_hd)
      obsLAT = obs_headElem_r(obsdata,OBS_LAT,index_hd)
      obsLON = obs_headElem_r(obsdata,OBS_LON,index_hd)
      obsDAT = obs_headElem_i(obsdata,OBS_DAT,index_hd)
      obsETM = obs_headElem_i(obsdata,OBS_ETM,index_hd)
      obsALT = obs_headElem_r(obsdata,OBS_ALT,index_hd)
      obsNLV = obs_headElem_i(obsdata,OBS_NLV,index_hd)
      obsSTYP= obs_headElem_i(obsdata,OBS_STYP,index_hd)
      obsTTYP= obs_headElem_i(obsdata,OBS_TTYP,index_hd)
      obsST1 = obs_headElem_i(obsdata,OBS_ST1,index_hd)
      stnid  = obs_elem_c(obsdata,'STID',index_hd)
      write(unitout_,fmt=9200) obsRLN, obsONM, obsINS, obsOTP, obsITY, &
           obsLAT, obsLON, obsDAT, obsETM, stnid, obsALT, obsNLV, &
           obsSTYP, obsTTYP, obsST1

9200  format(6x,'Position within realBodies:',i6,1x,'OBS. NUMBER:',i6,1x &
         ,'INSTR. ID:',i6,1x,'OBS. TYPE:',i6,1x &
         ,'INSTR./RETR. TYPE:',i6,1x &
         ,/,6x &
         ,'OBSERVATION LOCATION. (LAT,LON):',2(f10.4,1x) &
         ,'DATE:',i12,1x,'EXACT TIME: ',i6,1x &
         ,/,6x &
         ,'STATION ID:',a9,1x &
         ,'STATION''S ALTITUDE:',g13.6,1x &
         ,'NUMBER OF DATA:',i6,1x &
         ,/,6x &
         ,'SURFACE TYPE:',i4,1x &
         ,'TERRAIN TYPE:',i4,1x &
         ,'REPORT STATUS 2:',i10,1x &
         ,/,6x &
         )

      return
   end subroutine obs_prnthdr


   subroutine obs_read( obsdat, hx, nobshdr, nobsbdy, nobshx )
      !
      ! :Purpose: read the obsdat structure with observational information from
      !           unformatted files.  The files have been written by obs_write().
      !
      ! :Arguments:
      !           :obsdat:  obsSpaceData object
      !           :hx:      observation data
      !           :nobshdr: unit number of the file with obsdat header info.
      !           :nobsbdy: unit number of the file with obsdat body info.
      !           :nobshx:  unit number of the file with hx (-1 if not used)
      !
      ! :Note: It is assumed that the obsdat arrays have already been allocated
      !        with dimensions that will exactly hold the number of data to be read
      !
      implicit none

      type (struct_obs), intent(inout) :: obsdat
      integer          , intent(in)    :: nobshdr
      integer          , intent(in)    :: nobsbdy
      integer          , intent(in)    :: nobshx
      real(8)          , intent(out)   :: hx(:,:)

      integer  :: i,ifirst,ilast,iobscur,istn,j,k,nens
      integer  :: column_index
      character(len=100) :: message

      obsdat%mpi_local = .false.

      if (nobshx == -1) then
         nens=0
      else
         nens=size(hx,1)
      endif
                                        ! get index of 1st active RB column
      column_index=odc_columnIndexFromActiveIndex( &
                                                 obsdat%realBodies%odc_flavour,1)
      obsdat%numBody=size(obsdat%realBodies%columns(column_index)%value_r,1)
                                        ! get index of 1st active RH column
      column_index=odc_columnIndexFromActiveIndex( &
                                                obsdat%realHeaders%odc_flavour,1)
      obsdat%numHeader=size(obsdat%realHeaders%columns(column_index)%value_r,1)

      iobscur=0

      ! read stations

      readstn: do istn=1,obsdat%numHeader
         read(nobshdr,end=288,err=288) &
            (obsdat%intHeaders%columns(odc_columnIndexFromActiveIndex( &
                                               obsdat%intHeaders%odc_flavour,i) &
                                      )%value_i(istn),&
                 i=1,odc_numActiveColumn(obsdat%intHeaders)),&
            (obsdat%realHeaders%columns(odc_columnIndexFromActiveIndex( &
                                              obsdat%realHeaders%odc_flavour,j) &
                                       )%value_r(istn),&
                 j=1,odc_numActiveColumn(obsdat%realHeaders)),&
            obsdat%cstnid(istn), &
            obsdat%cfamily(istn)

         if (istn == 1) then 
            call obs_headSet_i(obsdat, OBS_RLN, istn, 1)
         else
            call obs_headSet_i (obsdat, OBS_RLN, istn, &
                 obs_headElem_i(obsdat, OBS_RLN, istn-1) &
                +obs_headElem_i(obsdat, OBS_NLV, istn-1))
         endif
         iobscur=iobscur+obs_headElem_i(obsdat, OBS_NLV, istn)
         ! now read the observations:
         ifirst=         obs_headElem_i(obsdat, OBS_RLN, istn)
         ilast =ifirst + obs_headElem_i(obsdat, OBS_NLV, istn) -1
         ! only read those columns specified by the user
         do i=ifirst,ilast
            read(nobsbdy) &
              (obsdat%intBodies%columns(odc_ENKF_bdy_int_column_list(j) &
                                       )%value_i(i),&
                   j=1,size(odc_ENKF_bdy_int_column_list(:))),&
              (obsdat%realBodies%columns(odc_ENKF_bdy_real_column_list(k) &
                                        )%value_r(i),&
                   k=1,size(odc_ENKF_bdy_real_column_list(:)))
         enddo
         if (nens > 0) then
            do i=ifirst,ilast
               read(nobshx)  (hx(j,i),j=1,nens)
            enddo
         endif
      enddo readstn

      if (iobscur /= obsdat%numBody) then
         write(message,*)'OBS_READ: the number of references in the header, ', &
                         iobscur, ', does not match the body size, ', &
                         obsdat%numBody
         call obs_abort(message); return
      endif
288   write(*,*) 'file is now empty'
      write(*,*) 'close nobshdr which is on unit: ',nobshdr
      close(nobshdr)
      write(*,*) 'close nobsbdy which is on unit: ',nobsbdy
      close(nobsbdy)
      if (nens > 0) then
         write(*,*) 'close nobshx which is on unit: ',nobshx
         close(nobshx)
      endif
      write(*,*) 'exit from obs_read'
      return

   end subroutine obs_read


   subroutine obs_creatSubCMA(cma,cma_sub,ipasscur,np)
   !
   ! :Purpose: create a sub-CMA from the global CMA.
   !   
   ! :Arguments: 
   !           :cma: the global CMA. 
   !           :ipasscur:  current pass number    
   !           :np:  number of processors  used .  
   !           :cma_sub:  the sub-CMA created from the global CMA.  
   !
   implicit none

   integer         , intent(in)    :: ipasscur
   integer         , intent(in)    :: np
   type(struct_obs), intent(in)    :: cma
   type(struct_obs), intent(inout) :: cma_sub

! local variables ..
   integer :: istn,ipass,ireg ,idata
   integer :: mxstn ,ifirst,ilast, i,j
   integer :: isize,rsize,active_index ,column_index
   integer :: ii0,iifirst,ii
   integer :: numheader, numBody
   integer :: intHeader
   integer :: iregcur
   real(pre_obsReal) :: realHeader
   integer,save :: istart=1
!
   iregcur=1
   numheader=0
   numBody=0
   mxstn=obs_numheader(cma)
   isize=odc_numActiveColumn(cma%intHeaders)
   rsize=odc_numActiveColumn(cma%realHeaders)
   do istn=istart,mxstn

      ipass=obs_headElem_i(cma,OBS_PAS,istn)
      !ireg=obs_headElem_i(cma,OBS_REG,istn)

      if (ipass /= ipasscur) then
         istart=istn
         !cycle
         exit
      endif
      ! default assignment of all input variables to cma_ipasscur .
      numheader=numheader+1
      cma_sub%numHeader=numheader
      idata=numheader
      do active_index=1,isize
            column_index=odc_columnIndexFromActiveIndex( &
                                     cma%intHeaders%odc_flavour, active_index)
            intHeader=cma%intHeaders%columns(column_index)%value_i(istn)
            if(cma%intHeaders%odc_flavour%columnActive(column_index))  &
               call obs_headSet_i(cma_sub, column_index, idata, &
                                  intHeader)
      enddo
!
      do active_index=1,rsize
            column_index=odc_columnIndexFromActiveIndex( &
                                    cma%realHeaders%odc_flavour, active_index)
            realHeader=cma%realHeaders%columns(column_index)%value_r(istn)
            if(cma%realHeaders%odc_flavour%columnActive(column_index))  &
               call obs_headSet_r(cma_sub, column_index, idata, &
                                  realHeader)
       enddo
       cma_sub%cstnid(idata)  =cma%cstnid(istn)
       cma_sub%cfamily(idata) =cma%cfamily(istn)

       ! determine which process will handle this station.
       ! the corresponding scatter operation is in program scattercma.
       !! call obs_headSet_i(cma_sub,OBS_IP,idata,mod(ipasscur,np))
       ! Use Round-robin obs. distribution.
         call obs_headSet_i(cma_sub,OBS_IP,idata, mod(obs_headElem_i(cma,OBS_ONM,istn),np))

         ireg=obs_headElem_i(cma_sub, OBS_REG, idata)
         if (iregcur /= ireg) then
            iregcur=ireg
         endif

         if (idata == 1) then
            call obs_headSet_i(cma_sub, OBS_RLN, idata, 1)
          else
            call obs_headSet_i(cma_sub, OBS_RLN, idata, &
                obs_headElem_i(cma_sub, OBS_RLN, idata-1) &
               +obs_headElem_i(cma_sub, OBS_NLV, idata-1))
         endif

         ! now read the bodies:
         ifirst=         obs_headElem_i(cma_sub, OBS_RLN, idata)
         ilast =ifirst + obs_headElem_i(cma_sub, OBS_NLV, idata)-1
         iifirst=obs_headElem_i(cma, OBS_RLN, istn)
         ii0=0
         do i=ifirst,ilast
            ii=iifirst+ii0
            do j=1,odc_numActiveColumn(cma%intBodies)
               column_index=odc_columnIndexFromActiveIndex( &
                            cma%intBodies%odc_flavour,j)
               cma_sub%intBodies%columns(column_index)%value_i(i)=  &
                           cma%intBodies%columns(column_index)%value_i(ii)
           enddo !===> j
!
           do j=1,odc_numActiveColumn(cma%realBodies)
               column_index=odc_columnIndexFromActiveIndex(&
                             cma%realBodies%odc_flavour,j)
              cma_sub%realBodies%columns(column_index)%value_r(i)= &
                  cma%realBodies%columns(column_index)%value_r(ii)
           enddo ! ===> j
              ! Make HIND point to new header row_index
            call obs_bodySet_i(cma_sub, OBS_HIND, i, idata)
            ii0=ii0+1
         enddo
         numBody=numBody+ii0


         ! go back to read the next station
      enddo !  =========>.
      cma_sub%numBody=numBody
      cma_sub%mpi_local=.false.

   end subroutine obs_creatSubCMA



   subroutine obs_reduceToMpiLocal(obsdat)
      !
      ! :Purpose: re-construct observation data object by
      !          giving local Obs TAG. 
      !      To retain in the observation object only those data that are
      !      pertinent to the present mpi processor, i.e. convert from mpiglobal
      !      to mpilocal.
      !
      implicit none

      type(struct_obs), intent(inout) :: obsdat

      ! Declare Local Variables
      integer(8),        allocatable,dimension(:)   :: headerPrimaryKey_tmp
      integer(8),        allocatable,dimension(:)   :: bodyPrimaryKey_tmp
      character(len=12), allocatable,dimension(:)   :: cstnid_tmp
      character(len=2),  allocatable,dimension(:)   :: cfamily_tmp
      real(pre_obsReal), allocatable,dimension(:,:) :: realHeaders_tmp
      real(pre_obsReal), allocatable,dimension(:,:) :: realBodies_tmp

      integer,allocatable,dimension(:,:) :: intHeaders_tmp,intBodies_tmp

      integer :: numHeader_mpilocal,numHeader_mpiglobal
      integer ::   numBody_mpilocal,  numBody_mpiglobal
      integer :: bodyIndex_mpilocal,bodyIndex_mpiglobal
      integer :: headerIndex_mpilocal,headerIndex_mpiglobal
      integer :: idataend,idata,active_index
      integer :: column_index

      !!---------------------------------------------------------------
      WRITE(*,*) '============= Enter obs_reduceToMpiLocal =============='

      if(obsdat%mpi_local)then
         call obs_abort('OBS_REDUCETOMPILOCAL() has been called, but the ' &
                        // 'obsSpaceData object is already in mpi-local state')
         return
      end if

      if(obsdat%numHeader_max == 0)then
         call obs_abort('OBS_REDUCETOMPILOCAL() has been called when there are '&
                        // 'no data.  Obs_reduceToMpiLocal cannot be called ' &
                        // 'after a call to obs_expandToMpiGlobal.')
         return
      end if

      ! compute the mpilocal lists of indices into the mpiglobal data
      call obs_mpiDistributeIndices(obsdat)

      ! calculate the size of the local obs data  
      if(associated(obsdat%headerIndex_mpiglobal)) then
         numHeader_mpilocal=size(obsdat%headerIndex_mpiglobal) 
      else
         numHeader_mpilocal=0
      endif
      numBody_mpiLocal=0
      do headerIndex_mpilocal=1,numHeader_mpilocal
         headerIndex_mpiglobal=obsdat%headerIndex_mpiglobal(headerIndex_mpilocal)
         idata=obs_headElem_i(obsdat, OBS_NLV, headerIndex_mpiglobal)
         numBody_mpiLocal = numBody_mpiLocal + idata
      enddo

      numHeader_mpiGlobal = obs_numHeader(obsdat)
      numBody_mpiGlobal   = obs_numBody(obsdat)

      ! allocate temporary arrays to hold mpilocal data
      if(numHeader_mpiLocal > 0) then
         allocate(headerPrimaryKey_tmp(numHeader_mpiLocal)) 
         allocate(cfamily_tmp(    numHeader_mpiLocal)) 
         allocate( cstnid_tmp(    numHeader_mpiLocal)) 
         allocate(realHeaders_tmp(odc_numActiveColumn(obsdat%realHeaders), &
                                  numHeader_mpilocal))
         allocate( intHeaders_tmp(odc_numActiveColumn(obsdat%intHeaders), &
                                  numHeader_mpilocal))
      endif

      if(numBody_mpiLocal > 0) then
         allocate(bodyPrimaryKey_tmp(numBody_mpiLocal)) 
         allocate( realBodies_tmp(odc_numActiveColumn(obsdat%realBodies), &
                                  numBody_mpilocal))
         allocate(  intBodies_tmp(odc_numActiveColumn(obsdat%intBodies), &
                                  numBody_mpilocal))
      endif

      ! copy the mpilocal data to temporary arrays: header-level data
      do headerIndex_mpilocal=1,numHeader_mpilocal 
         headerIndex_mpiglobal=obsdat%headerIndex_mpiglobal(headerIndex_mpilocal)

         do active_index=1,odc_numActiveColumn(obsdat%realHeaders)
            column_index=odc_columnIndexFromActiveIndex( &
                                    obsdat%realHeaders%odc_flavour, active_index)
            realHeaders_tmp(active_index,headerIndex_mpilocal)= &
                      obs_headElem_r(obsdat, column_index, headerIndex_mpiglobal)
         enddo

         do active_index=1,odc_numActiveColumn(obsdat%intHeaders)
            column_index=odc_columnIndexFromActiveIndex( &
                                     obsdat%intHeaders%odc_flavour, active_index)
            intHeaders_tmp(active_index,headerIndex_mpilocal)= &
                      obs_headElem_i(obsdat, column_index, headerIndex_mpiglobal)
         enddo

         headerPrimaryKey_tmp(headerIndex_mpilocal) = &
                obsdat%headerPrimaryKey(headerIndex_mpiglobal)

         cstnid_tmp (headerIndex_mpilocal) =obsdat%cstnid (headerIndex_mpiglobal)
         cfamily_tmp(headerIndex_mpilocal) =obsdat%cfamily(headerIndex_mpiglobal)

                                        ! Make RLN point to local data
         if(headerIndex_mpilocal == 1) then
            intHeaders_tmp &
                (odc_activeIndexFromColumnIndex( &
                                        obsdat%intHeaders%odc_flavour,OBS_RLN), &
                 1 &
                ) = 1
         else
            intHeaders_tmp &
                (odc_activeIndexFromColumnIndex( &
                                        obsdat%intHeaders%odc_flavour,OBS_RLN), &
                 headerIndex_mpilocal &
                ) = intHeaders_tmp(odc_activeIndexFromColumnIndex( &
                                        obsdat%intHeaders%odc_flavour,OBS_RLN), &
                                   headerIndex_mpilocal-1 &
                                  ) &
                  + intHeaders_tmp(odc_activeIndexFromColumnIndex( &
                                        obsdat%intHeaders%odc_flavour,OBS_NLV), &
                                   headerIndex_mpilocal-1 &
                                  ) 
         endif
      enddo

      ! copy the mpilocal data to temporary arrays: body-level data
      bodyIndex_mpilocal=0 
      do headerIndex_mpilocal=1,numHeader_mpilocal
         headerIndex_mpiglobal=obsdat%headerIndex_mpiglobal(headerIndex_mpilocal)
 
                                        ! Make HIND point to local header
         idata    = obs_headElem_i(obsdat, OBS_RLN,headerIndex_mpiglobal)
         idataend = obs_headElem_i(obsdat, OBS_NLV,headerIndex_mpiglobal)+idata-1
         do bodyIndex_mpiglobal=idata,idataend 
            bodyIndex_mpilocal=bodyIndex_mpilocal+1 

            bodyPrimaryKey_tmp(bodyIndex_mpilocal)= &
                        obsdat%bodyPrimaryKey(bodyIndex_mpiglobal)

            do active_index=1,odc_numActiveColumn(obsdat%realBodies)
               column_index=odc_columnIndexFromActiveIndex( &
                                      obsdat%realBodies%odc_flavour,active_index)
               realBodies_tmp(active_index,bodyIndex_mpilocal)= &
                        obs_bodyElem_r(obsdat, column_index, bodyIndex_mpiglobal)
            enddo

            do active_index=1,odc_numActiveColumn(obsdat%intBodies)
               column_index=odc_columnIndexFromActiveIndex( &
                                       obsdat%intBodies%odc_flavour,active_index)
               intBodies_tmp(active_index,bodyIndex_mpilocal)= &
                        obs_bodyElem_i(obsdat, column_index, bodyIndex_mpiglobal)
            enddo

            intBodies_tmp(odc_activeIndexFromColumnIndex( &
                                       obsdat%intBodies%odc_flavour, OBS_HIND), &
                          bodyIndex_mpilocal &
                         ) = headerIndex_mpilocal
         enddo
      enddo

      ! destroy object's mpiglobal data and allocate mpilocal data
      obsdat%numHeader=numHeader_mpiLocal
      obsdat%numBody=numBody_mpiLocal
      call obs_deallocate(obsdat)
      call obs_allocate(obsdat,obsdat%numHeader,obsdat%numBody)

      ! copy all data from temporary arrays to object's arrays
      HEADER:if(numHeader_mpiLocal > 0) then
         obsdat%headerPrimaryKey(:)=headerPrimaryKey_tmp(:  ) 
         obsdat%cfamily(:  )=cfamily_tmp(:  ) 
         obsdat%cstnid (:  )= cstnid_tmp(:  )
         do active_index=1,odc_numActiveColumn(obsdat%realHeaders)
            column_index=odc_columnIndexFromActiveIndex( &
                                    obsdat%realHeaders%odc_flavour, active_index)
            obsdat%realHeaders%columns(column_index)%value_r(:) &
                                                 =realHeaders_tmp(active_index,:)
         enddo
         do active_index=1,odc_numActiveColumn(obsdat%intHeaders)
            column_index=odc_columnIndexFromActiveIndex( &
                                     obsdat%intHeaders%odc_flavour, active_index)
            obsdat%intHeaders%columns(column_index)%value_i(:) &
                                                  =intHeaders_tmp(active_index,:)
         enddo

         ! deallocate temporary arrays
         deallocate(headerPrimaryKey_tmp)
         deallocate(cfamily_tmp)
         deallocate(cstnid_tmp)
         deallocate(realHeaders_tmp)
         deallocate(intHeaders_tmp)
      endif HEADER

      BODY:if(numBody_mpiLocal > 0) then
         obsdat%bodyPrimaryKey(:) = bodyPrimaryKey_tmp(:)
         do active_index=1,odc_numActiveColumn(obsdat%realBodies)
            column_index=odc_columnIndexFromActiveIndex( &
                                     obsdat%realBodies%odc_flavour, active_index)
            obsdat%realBodies%columns(column_index)%value_r(:) &
                                                  =realBodies_tmp(active_index,:)
         enddo
         do active_index=1,odc_numActiveColumn(obsdat%intBodies)
            column_index=odc_columnIndexFromActiveIndex( &
                                      obsdat%intBodies%odc_flavour, active_index)
            obsdat%intBodies%columns(column_index)%value_i(:) &
                                                  =intBodies_tmp(active_index,:) 
         enddo

         ! deallocate temporary arrays
         deallocate(bodyPrimaryKey_tmp)
         deallocate(realBodies_tmp)
         deallocate(intBodies_tmp)
      endif BODY


      obsdat%mpi_local = .true.

      write(*,*) '============= Leave obs_reduceToMpiLocal =============='

      return
   end subroutine obs_reduceToMpiLocal


   subroutine obs_squeeze( obsdat )
      !
      ! :Purpose: re-construct observation data object to save memory
      !
      implicit none

      type(struct_obs), intent(inout) :: obsdat

      ! Declare Local Variables
      integer(8),        allocatable :: headerPrimaryKey_tmp(:)
      integer(8),        allocatable :: bodyPrimaryKey_tmp(:)
      character(len=12), allocatable :: cstnid_tmp(:)
      character(len=2),  allocatable :: cfamily_tmp(:)
      real(pre_obsReal), allocatable :: realHeaders_tmp(:,:), realBodies_tmp(:,:)
      integer,           allocatable :: intHeaders_tmp(:,:),intBodies_tmp(:,:)
      integer :: bodyIndex, headerIndex, active_index, column_index
      integer :: idataend, idata

      write(*,*) '============= Enter obs_squeeze =============='

      if(obsdat%numHeader_max == obsdat%numHeader .and. &
         obsdat%numBody_max == obsdat%numBody)then
         write(*,*) 'obs_squeeze: obsdata instance is already sized correctly, do nothing'
         return
      end if

      ! allocate temporary arrays to hold data
      if(obsdat%numHeader > 0) then
         allocate(headerPrimaryKey_tmp(obsdat%numHeader))
         allocate(cfamily_tmp(    obsdat%numHeader)) 
         allocate( cstnid_tmp(    obsdat%numHeader)) 
         allocate(realHeaders_tmp(odc_numActiveColumn(obsdat%realHeaders), &
                                  obsdat%numHeader))
         allocate( intHeaders_tmp(odc_numActiveColumn(obsdat%intHeaders), &
                                  obsdat%numHeader))
      endif

      if(obsdat%numBody > 0) then
         allocate(bodyPrimaryKey_tmp(obsdat%numBody))
         allocate( realBodies_tmp(odc_numActiveColumn(obsdat%realBodies), &
                                  obsdat%numBody))
         allocate(  intBodies_tmp(odc_numActiveColumn(obsdat%intBodies), &
                                  obsdat%numBody))
      endif

      ! copy the data to temporary arrays: header-level data
      do headerIndex=1,obsdat%numHeader 
         headerPrimaryKey_tmp(headerIndex)= &
                      obsdat%headerPrimaryKey(headerIndex)
        
         do active_index=1,odc_numActiveColumn(obsdat%realHeaders)
            column_index=odc_columnIndexFromActiveIndex( &
                                    obsdat%realHeaders%odc_flavour, active_index)
            realHeaders_tmp(active_index,headerIndex)= &
                      obs_headElem_r(obsdat, column_index, headerIndex)
         enddo

         do active_index=1,odc_numActiveColumn(obsdat%intHeaders)
            column_index=odc_columnIndexFromActiveIndex( &
                                     obsdat%intHeaders%odc_flavour, active_index)
            intHeaders_tmp(active_index,headerIndex)= &
                      obs_headElem_i(obsdat, column_index, headerIndex)
         enddo

         cstnid_tmp (headerIndex) =obsdat%cstnid (headerIndex)
         cfamily_tmp(headerIndex) =obsdat%cfamily(headerIndex)

      enddo

      ! copy the data to temporary arrays: body-level data
      do headerIndex=1,obsdat%numHeader
                                        ! Make HIND point to local header
         idata    = obs_headElem_i(obsdat, OBS_RLN,headerIndex)
         idataend = obs_headElem_i(obsdat, OBS_NLV,headerIndex)+idata-1
         do bodyIndex=idata,idataend 
            bodyPrimaryKey_tmp(bodyIndex)= &
                        obsdat%bodyPrimaryKey(bodyIndex)
           
            do active_index=1,odc_numActiveColumn(obsdat%realBodies)
               column_index=odc_columnIndexFromActiveIndex( &
                                      obsdat%realBodies%odc_flavour,active_index)
               realBodies_tmp(active_index,bodyIndex)= &
                        obs_bodyElem_r(obsdat, column_index, bodyIndex)
            enddo

            do active_index=1,odc_numActiveColumn(obsdat%intBodies)
               column_index=odc_columnIndexFromActiveIndex( &
                                       obsdat%intBodies%odc_flavour,active_index)
               intBodies_tmp(active_index,bodyIndex)= &
                        obs_bodyElem_i(obsdat, column_index, bodyIndex)
            enddo

            intBodies_tmp(odc_activeIndexFromColumnIndex(obsdat%intBodies%odc_flavour, OBS_HIND), &
                          bodyIndex) = headerIndex
         enddo
      enddo

      ! destroy object's data and allocate at correct size data
      call obs_deallocate(obsdat)
      call obs_allocate(obsdat,obsdat%numHeader,obsdat%numBody)

      ! copy all data from temporary arrays to object's arrays
      HEADER:if(obsdat%numHeader > 0) then
         obsdat%headerPrimaryKey(:)=headerPrimaryKey_tmp(:)
         obsdat%cfamily(:)=cfamily_tmp(:) 
         obsdat%cstnid (:)= cstnid_tmp(:)
         do active_index=1,odc_numActiveColumn(obsdat%realHeaders)
            column_index=odc_columnIndexFromActiveIndex( &
                                    obsdat%realHeaders%odc_flavour, active_index)
            obsdat%realHeaders%columns(column_index)%value_r(:) &
                                                 =realHeaders_tmp(active_index,:)
         enddo
         do active_index=1,odc_numActiveColumn(obsdat%intHeaders)
            column_index=odc_columnIndexFromActiveIndex( &
                                     obsdat%intHeaders%odc_flavour, active_index)
            obsdat%intHeaders%columns(column_index)%value_i(:) &
                                                  =intHeaders_tmp(active_index,:)
         enddo

         ! deallocate temporary arrays
         deallocate(headerPrimaryKey_tmp)
         deallocate(cfamily_tmp)
         deallocate(cstnid_tmp)
         deallocate(realHeaders_tmp)
         deallocate(intHeaders_tmp)
      endif HEADER

      BODY:if(obsdat%numBody > 0) then
         obsdat%bodyPrimaryKey(:) = bodyPrimaryKey_tmp(:)
         do active_index=1,odc_numActiveColumn(obsdat%realBodies)
            column_index=odc_columnIndexFromActiveIndex( &
                                     obsdat%realBodies%odc_flavour, active_index)
            obsdat%realBodies%columns(column_index)%value_r(:) &
                                                  =realBodies_tmp(active_index,:)
         enddo
         do active_index=1,odc_numActiveColumn(obsdat%intBodies)
            column_index=odc_columnIndexFromActiveIndex( &
                                      obsdat%intBodies%odc_flavour, active_index)
            obsdat%intBodies%columns(column_index)%value_i(:) &
                                                  =intBodies_tmp(active_index,:) 
         enddo

         ! deallocate temporary arrays
         deallocate(bodyPrimaryKey_tmp)
         deallocate(realBodies_tmp)
         deallocate(intBodies_tmp)
      endif BODY

      write(*,*) '============= Leave obs_squeeze =============='

   end subroutine obs_squeeze


   subroutine obs_MpiRedistribute( obsdat_inout, target_ip_index )
      !
      ! :Purpose: Redistribute obs over mpi tasks according to mpi task id stored 
      !                            in the integer header column "target_ip_index"
      !
      implicit none

      type(struct_obs), intent(inout) :: obsdat_inout
      integer :: target_ip_index

      ! Declare Local Variables
      type(struct_obs) :: obsdat_tmp
      integer, allocatable :: numHeaderPE_mpilocal(:), numHeaderPE_mpiglobal(:)
      integer, allocatable :: numBodyPE_mpilocal(:), numBodyPE_mpiglobal(:)
      integer,           allocatable :: intcstnid_send(:,:,:), intcstnid_recv(:,:,:)
      integer,           allocatable :: intcfamily_send(:,:,:), intcfamily_recv(:,:,:)
      integer(8),        allocatable :: primaryKey_send(:,:), primaryKey_recv(:,:)
      real(pre_obsReal), allocatable :: real_send(:,:,:), real_recv(:,:,:)
      integer,           allocatable :: int_send(:,:,:), int_recv(:,:,:), message_onm(:,:)
      real(pre_obsReal), allocatable :: real_send_2d(:,:), real_recv_2d(:,:)
      integer,           allocatable :: int_send_2d(:,:), int_recv_2d(:,:)
      integer :: numHeader_in, numBody_in
      integer :: numHeader_out, numBody_out
      integer :: numHeader_mpimessage, numBody_mpimessage
      integer :: bodyIndex, headerIndex, headerIndex_out, bodyIndex_out, columnIndex, activeIndex, procIndex, charIndex
      integer :: bodyIndexBeg, bodyIndexEnd
      integer :: nprocs_mpi, myid_mpi, ierr, nsize, target_ip
      logical :: needToRedistribute, needToRedistribute_mpiglobal
      integer :: get_max_rss

      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
      write(*,*) '============= Enter obs_MpiRedistribute =============='
      write(*,*) 'redistribute data according to mpi task ID stored in column :', &
                 ocn_ColumnNameList_IH(target_ip_index)

      ! determine rank and number of mpi tasks
      call rpn_comm_size("GRID",nprocs_mpi,ierr)
      call rpn_comm_rank("GRID",myid_mpi,ierr)


      ! Number of headers and bodies per task before redistribution
      numHeader_in = obs_numHeader(obsdat_inout)
      numBody_in   = obs_numBody(obsdat_inout)

      ! check if redistribution is really needed
      needToRedistribute = .false.
      do headerIndex = 1, numHeader_in
         target_ip = obs_headElem_i(obsdat_inout,target_ip_index,headerIndex)
         if (target_ip /= myid_mpi) needToRedistribute = .true.
      enddo
      call rpn_comm_allreduce(needToRedistribute,needToRedistribute_mpiglobal,1,  &
                              "MPI_LOGICAL","MPI_LOR","world",ierr)
      if(.not.needToRedistribute_mpiglobal) then
         write(*,*) 'obs_MpiRedistribute: do not need to redistribute, returning'
         return
      endif

      ! allocate arrays used for counting on each mpi task
      allocate(numHeaderPE_mpilocal(nprocs_mpi))
      allocate(numHeaderPE_mpiglobal(nprocs_mpi))
      allocate(numBodyPE_mpilocal(nprocs_mpi))
      allocate(numBodyPE_mpiglobal(nprocs_mpi))

      ! Compute number of headers and bodies per task after redistribution
      numHeaderPE_mpilocal(:) = 0
      numBodyPE_mpilocal(:) = 0
      do headerIndex = 1, numHeader_in
         target_ip = obs_headElem_i(obsdat_inout,target_ip_index,headerIndex)
         numHeaderPE_mpilocal(1+target_ip) = numHeaderPE_mpilocal(1+target_ip) + 1
         numBodyPE_mpilocal(1+target_ip)   = numBodyPE_mpilocal(1+target_ip)   + obs_headElem_i(obsdat_inout,OBS_NLV,headerIndex)
      enddo
      call rpn_comm_allreduce(numHeaderPE_mpilocal,numHeaderPE_mpiglobal,nprocs_mpi,  &
                              "MPI_INTEGER","MPI_SUM","world",ierr)
      call rpn_comm_allreduce(numBodyPE_mpilocal,numBodyPE_mpiglobal,nprocs_mpi,  &
                              "MPI_INTEGER","MPI_SUM","world",ierr)
      numHeader_out = numHeaderPE_mpiglobal(myid_mpi+1)
      numBody_out   = numBodyPE_mpiglobal(myid_mpi+1)
      write(*,*) 'obs_MpiRedistribute: num mpi header and body before redistribution =', numHeader_in, numBody_in
      write(*,*) 'obs_MpiRedistribute: num mpi header and body after redistribution  =', numHeader_out, numBody_out

      ! Compute the max number of headers and bodies in each mpi message sent/received in the transpose
      call rpn_comm_allreduce(numHeaderPE_mpilocal,numHeaderPE_mpiglobal,nprocs_mpi,  &
                              "MPI_INTEGER","MPI_MAX","GRID",ierr)
      call rpn_comm_allreduce(numBodyPE_mpilocal,numBodyPE_mpiglobal,nprocs_mpi,  &
                              "MPI_INTEGER","MPI_MAX","GRID",ierr)
      if(myid_mpi == 0) write(*,*) 'obs_MpiRedistribute: num mpi header messages =', numHeaderPE_mpilocal
      if(myid_mpi == 0) write(*,*) 'obs_MpiRedistribute: num mpi body messages =', numBodyPE_mpilocal
      if(myid_mpi == 0) write(*,*) 'obs_MpiRedistribute: num mpi header messages (max) =', numHeaderPE_mpiglobal
      if(myid_mpi == 0) write(*,*) 'obs_MpiRedistribute: num mpi body messages (max) =', numBodyPE_mpiglobal
      numHeader_mpimessage = maxval(numHeaderPE_mpiglobal(:))
      numBody_mpimessage   = maxval(numBodyPE_mpiglobal(:))

      call obs_initialize(obsdat_tmp,numHeader_out,  &
                          numBody_out,mpi_local=.true.)
      obsdat_tmp%numHeader = numHeader_out
      obsdat_tmp%numBody   = numBody_out

      ! allocate temporary arrays to hold header-level data for mpi communication
      allocate(intcfamily_send(len(obsdat_inout%cfamily(1)),numHeader_mpimessage,nprocs_mpi)) 
      allocate(intcfamily_recv(len(obsdat_inout%cfamily(1)),numHeader_mpimessage,nprocs_mpi)) 

      allocate(intcstnid_send(len(obsdat_inout%cstnid(1)),numHeader_mpimessage,nprocs_mpi)) 
      allocate(intcstnid_recv(len(obsdat_inout%cstnid(1)),numHeader_mpimessage,nprocs_mpi)) 

      allocate(real_send(odc_numActiveColumn(obsdat_inout%realHeaders), &
                         numHeader_mpimessage,nprocs_mpi))
      allocate(real_recv(odc_numActiveColumn(obsdat_inout%realHeaders), &
                         numHeader_mpimessage,nprocs_mpi))
      allocate(int_send(odc_numActiveColumn(obsdat_inout%intHeaders), &
                        numHeader_mpimessage,nprocs_mpi))
      allocate(int_recv(odc_numActiveColumn(obsdat_inout%intHeaders), &
                        numHeader_mpimessage,nprocs_mpi))

      allocate(primaryKey_send(numHeader_mpimessage,nprocs_mpi))
      allocate(primaryKey_recv(numHeader_mpimessage,nprocs_mpi))

      ! copy the data to temporary arrays: header-level data
      numHeaderPE_mpilocal(:) = 0
      int_send(:,:,:) = -99999
      do headerIndex=1,numHeader_in
         target_ip = obs_headElem_i(obsdat_inout,target_ip_index,headerIndex)
         numHeaderPE_mpilocal(1+target_ip) = numHeaderPE_mpilocal(1+target_ip) + 1

         primaryKey_send(numHeaderPE_mpilocal(1+target_ip),1+target_ip)= &
               obsdat_inout%headerPrimaryKey(headerIndex)

         do activeIndex=1,odc_numActiveColumn(obsdat_inout%realHeaders)
            columnIndex=odc_columnIndexFromActiveIndex(obsdat_inout%realHeaders%odc_flavour, &
                                                       activeIndex)
            real_send(activeIndex,numHeaderPE_mpilocal(1+target_ip),1+target_ip)= &
               obs_headElem_r(obsdat_inout, columnIndex, headerIndex)
         enddo

         do activeIndex=1,odc_numActiveColumn(obsdat_inout%intHeaders)
            columnIndex=odc_columnIndexFromActiveIndex(obsdat_inout%intHeaders%odc_flavour, &
                                                       activeIndex)
            int_send(activeIndex,numHeaderPE_mpilocal(1+target_ip),1+target_ip)= &
               obs_headElem_i(obsdat_inout, columnIndex, headerIndex)
         enddo

         do charIndex=1, len(obsdat_inout%cstnid(1))
            intcstnid_send(charIndex,numHeaderPE_mpilocal(1+target_ip),1+target_ip)= &
               iachar(obsdat_inout%cstnid(headerIndex)(charIndex:charIndex))
         enddo

         do charIndex=1,len(obsdat_inout%cfamily(1))
            intcfamily_send(charIndex,numHeaderPE_mpilocal(1+target_ip),1+target_ip)=  &
               iachar(obsdat_inout%cfamily(headerIndex)(charIndex:charIndex))
         enddo

      enddo

      ! do mpi communication: header-level data
      if(nprocs_mpi > 1) then
        nsize = numHeader_mpimessage
        call rpn_comm_alltoall(primaryKey_send,nsize,"mpi_integer8",  &
                               primaryKey_recv,nsize,"mpi_integer8","GRID",ierr)

        nsize = numHeader_mpimessage*odc_numActiveColumn(obsdat_inout%realHeaders)
        call rpn_comm_alltoall(real_send,nsize,"mpi_double_precision",  &
                               real_recv,nsize,"mpi_double_precision","GRID",ierr)

        nsize = numHeader_mpimessage*odc_numActiveColumn(obsdat_inout%intHeaders)
        call rpn_comm_alltoall(int_send,nsize,"mpi_integer",  &
                               int_recv,nsize,"mpi_integer","GRID",ierr)

        nsize = numHeader_mpimessage*len(obsdat_inout%cstnid(1))
        call rpn_comm_alltoall(intcstnid_send,nsize,"mpi_integer",  &
                               intcstnid_recv,nsize,"mpi_integer","GRID",ierr)

        nsize = numHeader_mpimessage*len(obsdat_inout%cfamily(1))
        call rpn_comm_alltoall(intcfamily_send,nsize,"mpi_integer",  &
                               intcfamily_recv,nsize,"mpi_integer","GRID",ierr)
      else
        primaryKey_recv(:,1)   = primaryKey_send(:,1)
        real_recv(:,:,1)       = real_send(:,:,1)
        int_recv(:,:,1)        = int_send(:,:,1)
        intcstnid_recv(:,:,1)  = intcstnid_send(:,:,1)
        intcfamily_recv(:,:,1) = intcfamily_send(:,:,1)
      endif

      allocate(message_onm(numHeader_mpimessage,nprocs_mpi))
      activeIndex = odc_activeIndexFromColumnIndex(obsdat_inout%intHeaders%odc_flavour,OBS_ONM)
      message_onm(:,:) = int_recv(activeIndex,:,:)

      ! copy the data from temporary arrays: header-level data
      headerIndex_out = 0
      do procIndex = 1, nprocs_mpi
         do headerIndex=1,numHeader_mpimessage
            if(int_recv(1,headerIndex,procIndex) /= -99999) then
               if(target_ip_index == OBS_IPF) then
                 ! put headers back in original order as in the files
                 headerIndex_out = message_onm(headerIndex,procIndex)
               else
                 headerIndex_out = headerIndex_out + 1
               endif

               obsdat_tmp%headerPrimaryKey(headerIndex_out)= &
                     primaryKey_recv(headerIndex,procIndex)

               do activeIndex=1,odc_numActiveColumn(obsdat_inout%realHeaders)
                  columnIndex=odc_columnIndexFromActiveIndex(obsdat_inout%realHeaders%odc_flavour, &
                                                             activeIndex)
                  obsdat_tmp%realHeaders%columns(columnIndex)%value_r(headerIndex_out)= &
                     real_recv(activeIndex,headerIndex,procIndex)
               enddo

               do activeIndex=1,odc_numActiveColumn(obsdat_inout%intHeaders)
                  columnIndex=odc_columnIndexFromActiveIndex(obsdat_inout%intHeaders%odc_flavour, &
                                                             activeIndex)
                  obsdat_tmp%intHeaders%columns(columnIndex)%value_i(headerIndex_out)= &
                     int_recv(activeIndex,headerIndex,procIndex)
               enddo

               do charIndex=1, len(obsdat_inout%cstnid(1))
                  obsdat_tmp%cstnid(headerIndex_out)(charIndex:charIndex) = &
                     achar(intcstnid_recv(charIndex,headerIndex,procIndex))
               enddo

               do charIndex=1, len(obsdat_inout%cfamily(1))
                  obsdat_tmp%cfamily(headerIndex_out)(charIndex:charIndex) = &
                     achar(intcfamily_recv(charIndex,headerIndex,procIndex))
               enddo

            endif
         enddo
      enddo

      ! recompute RLN to point to correct body data
      do headerIndex=1,numHeader_out
         if(headerIndex == 1) then
            obsdat_tmp%intHeaders%columns(OBS_RLN)%value_i(headerIndex) = 1
         else
            obsdat_tmp%intHeaders%columns(OBS_RLN)%value_i(headerIndex) = &
                obsdat_tmp%intHeaders%columns(OBS_RLN)%value_i(headerIndex-1) + &
                obsdat_tmp%intHeaders%columns(OBS_NLV)%value_i(headerIndex-1)
         endif
      enddo

      deallocate(primaryKey_send)
      deallocate(primaryKey_recv)
      deallocate(intcfamily_send)
      deallocate(intcfamily_recv)
      deallocate(intcstnid_send)
      deallocate(intcstnid_recv)
      deallocate(real_send)
      deallocate(real_recv)
      deallocate(int_send)
      deallocate(int_recv)

      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      ! Do communication for bodyPrimaryKey

      allocate(primaryKey_send(numBody_mpimessage,nprocs_mpi))
      allocate(primaryKey_recv(numBody_mpimessage,nprocs_mpi))

      numBodyPE_mpilocal(:) = 0
      primaryKey_send(:,:) = -99999
      do bodyIndex=1,numBody_in
         headerIndex = obs_bodyElem_i(obsdat_inout,OBS_HIND,bodyIndex)
         target_ip = obs_headElem_i(obsdat_inout,target_ip_index,headerIndex)
         numBodyPE_mpilocal(1+target_ip) = numBodyPE_mpilocal(1+target_ip) + 1
         primaryKey_send(numBodyPE_mpilocal(1+target_ip),1+target_ip)= &
              obsdat_inout%bodyPrimaryKey(bodyIndex)
      enddo
      if(nprocs_mpi > 1) then
         nsize = numBody_mpimessage
         call rpn_comm_alltoall(primaryKey_send,nsize,"mpi_integer8",  &
                                primaryKey_recv,nsize,"mpi_integer8","GRID",ierr)
      else
         primaryKey_recv(:,1) = primaryKey_send(:,1)
      endif
      if(target_ip_index == OBS_IPF) then
         ! copy the data in the same order as in the original files
         do procIndex = 1, nprocs_mpi
            bodyIndex = 0
            do headerIndex=1,numHeader_mpimessage
               headerIndex_out = message_onm(headerIndex,procIndex)
               if(headerIndex_out /= -99999) then
                  bodyIndexBeg = obs_headElem_i(obsdat_tmp,OBS_RLN,headerIndex_out)
                  bodyIndexEnd = obs_headElem_i(obsdat_tmp,OBS_NLV,headerIndex_out) + bodyIndexBeg - 1
                  do bodyIndex_out = bodyIndexBeg, bodyIndexEnd
                     bodyIndex = bodyIndex + 1
                     obsdat_tmp%bodyPrimaryKey(bodyIndex_out)= &
                          primaryKey_recv(bodyIndex,procIndex)
                  enddo
               endif
            enddo
         enddo
         ! copy the data in sequential order
         bodyIndex_out = 0
         do procIndex = 1, nprocs_mpi
            do bodyIndex=1,numBody_mpimessage
               if(primaryKey_recv(bodyIndex,procIndex) /= -99999) then
                  bodyIndex_out = bodyIndex_out + 1
                  obsdat_tmp%bodyPrimaryKey(bodyIndex_out)= &
                       primaryKey_recv(bodyIndex,procIndex)
               endif
            enddo
         enddo
      endif

      deallocate(primaryKey_send)
      deallocate(primaryKey_recv)

      ! First do REAL body columns

      allocate(real_send_2d(numBody_mpimessage,nprocs_mpi))
      allocate(real_recv_2d(numBody_mpimessage,nprocs_mpi))

      do activeIndex=1,odc_numActiveColumn(obsdat_inout%realBodies)
         columnIndex=odc_columnIndexFromActiveIndex(obsdat_inout%realBodies%odc_flavour, &
                                                    activeIndex)

         ! copy the data to temporary arrays: body-level data
         numBodyPE_mpilocal(:) = 0
         real_send_2d(:,:) = -99999.0d0
         do bodyIndex=1,numBody_in
            headerIndex = obs_bodyElem_i(obsdat_inout,OBS_HIND,bodyIndex)
            target_ip = obs_headElem_i(obsdat_inout,target_ip_index,headerIndex)
            numBodyPE_mpilocal(1+target_ip) = numBodyPE_mpilocal(1+target_ip) + 1
            real_send_2d(numBodyPE_mpilocal(1+target_ip),1+target_ip)= &
               obs_bodyElem_r(obsdat_inout, columnIndex, bodyIndex)
         enddo

         ! do mpi communication: body-level data
         if(nprocs_mpi > 1) then
           nsize = numBody_mpimessage
           call rpn_comm_alltoall(real_send_2d,nsize,"mpi_double_precision",  &
                                  real_recv_2d,nsize,"mpi_double_precision","GRID",ierr)
         else
           real_recv_2d(:,1)       = real_send_2d(:,1)
         endif

         ! copy the data from temporary arrays: body-level data
         if(target_ip_index == OBS_IPF) then

            ! copy the data in the same order as in the original files
            do procIndex = 1, nprocs_mpi
               bodyIndex = 0
               do headerIndex=1,numHeader_mpimessage
                  headerIndex_out = message_onm(headerIndex,procIndex)
                  if(headerIndex_out /= -99999) then
                     bodyIndexBeg = obs_headElem_i(obsdat_tmp,OBS_RLN,headerIndex_out)
                     bodyIndexEnd = obs_headElem_i(obsdat_tmp,OBS_NLV,headerIndex_out) + bodyIndexBeg - 1
                     do bodyIndex_out = bodyIndexBeg, bodyIndexEnd
                        bodyIndex = bodyIndex + 1
                        obsdat_tmp%realBodies%columns(columnIndex)%value_r(bodyIndex_out)= &
                           real_recv_2d(bodyIndex,procIndex)
                     enddo
                  endif
               enddo
            enddo

         else

            ! copy the data in sequential order
            bodyIndex_out = 0
            do procIndex = 1, nprocs_mpi
               do bodyIndex=1,numBody_mpimessage
                  if(real_recv_2d(bodyIndex,procIndex) /= -99999.0d0) then
                     bodyIndex_out = bodyIndex_out + 1
                     obsdat_tmp%realBodies%columns(columnIndex)%value_r(bodyIndex_out)= &
                        real_recv_2d(bodyIndex,procIndex)
                  endif
               enddo
            enddo

         endif

      enddo ! activeIndex

      deallocate(real_send_2d)
      deallocate(real_recv_2d)

      ! Now do INTEGER body columns

      allocate(int_send_2d(numBody_mpimessage,nprocs_mpi))
      allocate(int_recv_2d(numBody_mpimessage,nprocs_mpi))

      do activeIndex=1,odc_numActiveColumn(obsdat_inout%intBodies)
         columnIndex=odc_columnIndexFromActiveIndex(obsdat_inout%intBodies%odc_flavour, &
                                                    activeIndex)

         ! copy the data to temporary arrays: body-level data
         numBodyPE_mpilocal(:) = 0
         int_send_2d(:,:) = -99999
         do bodyIndex=1,numBody_in
            headerIndex = obs_bodyElem_i(obsdat_inout,OBS_HIND,bodyIndex)
            target_ip = obs_headElem_i(obsdat_inout,target_ip_index,headerIndex)
            numBodyPE_mpilocal(1+target_ip) = numBodyPE_mpilocal(1+target_ip) + 1
            int_send_2d(numBodyPE_mpilocal(1+target_ip),1+target_ip)= &
               obs_bodyElem_i(obsdat_inout, columnIndex, bodyIndex)
         enddo

         ! do mpi communication: body-level data
         if(nprocs_mpi > 1) then
           nsize = numBody_mpimessage
           call rpn_comm_alltoall(int_send_2d,nsize,"mpi_integer",  &
                                  int_recv_2d,nsize,"mpi_integer","GRID",ierr)
         else
           int_recv_2d(:,1)        = int_send_2d(:,1)
         endif

         ! copy the data from temporary arrays: body-level data
         if(target_ip_index == OBS_IPF) then

            ! copy the data in the same order as in the original files
            do procIndex = 1, nprocs_mpi
               bodyIndex = 0
               do headerIndex=1,numHeader_mpimessage
                  headerIndex_out = message_onm(headerIndex,procIndex)
                  if(headerIndex_out /= -99999) then
                     bodyIndexBeg = obs_headElem_i(obsdat_tmp,OBS_RLN,headerIndex_out)
                     bodyIndexEnd = obs_headElem_i(obsdat_tmp,OBS_NLV,headerIndex_out) + bodyIndexBeg - 1
                     do bodyIndex_out = bodyIndexBeg, bodyIndexEnd
                        bodyIndex = bodyIndex + 1
                        obsdat_tmp%intBodies%columns(columnIndex)%value_i(bodyIndex_out)= &
                           int_recv_2d(bodyIndex,procIndex)
                     enddo
                  endif
               enddo
            enddo

         else

            ! copy the data in sequential order
            bodyIndex_out = 0
            do procIndex = 1, nprocs_mpi
               do bodyIndex=1,numBody_mpimessage
                  if(int_recv_2d(bodyIndex,procIndex) /= -99999) then
                     bodyIndex_out = bodyIndex_out + 1
                     obsdat_tmp%intBodies%columns(columnIndex)%value_i(bodyIndex_out)= &
                        int_recv_2d(bodyIndex,procIndex)
                  endif
               enddo
            enddo

         endif

      enddo ! activeIndex

      deallocate(int_send_2d)
      deallocate(int_recv_2d)

      ! reset the headerIndex within the body
      do headerIndex = 1, obs_numheader(obsdat_tmp)
        bodyIndexBeg = obs_headElem_i(obsdat_tmp,OBS_RLN,headerIndex)
        bodyIndexEnd = obs_headElem_i(obsdat_tmp,OBS_NLV,headerIndex) + bodyIndexBeg - 1
        do bodyIndex = bodyIndexBeg, bodyIndexEnd
          call obs_bodySet_i(obsdat_tmp,OBS_HIND,bodyIndex, headerIndex)
        enddo
      enddo

      deallocate(message_onm)

      deallocate(numHeaderPE_mpilocal)
      deallocate(numHeaderPE_mpiglobal)
      deallocate(numBodyPE_mpilocal)
      deallocate(numBodyPE_mpiglobal)

      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      ! reallocate obsdat_inout to the new size and copy over redistributed data
      call obs_deallocate(obsdat_inout)
      call obs_allocate(obsdat_inout,obsdat_tmp%numHeader,obsdat_tmp%numBody)
      call obs_copy(obsdat_tmp,obsdat_inout)
      call obs_finalize(obsdat_tmp)

      write(*,*) '============= Leave obs_MpiRedistribute =============='
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      return
   end subroutine obs_MpiRedistribute


   subroutine obs_select( obsdat ,hx, obs_sel, hx_sel, zhamin, zhamax, nens, nobsout )
      !
      ! :Purpose: select only the observations with zhamin < lop(P) <= zhamax.   
      !
      ! :Arguments:
      !           :obsdat:  input obsSpaceData object 
      !           :hx:      interpolated values of obsdat
      !           :obs_sel: selected obsdat 
      !           :hx_sel:  interpolated values of obs_sel
      !           :zhamin,zhamax:  range of zha values to be selected.
      !           :nens: number of ensemble members
      !           :nobsout: unit number for the ASCII output
      !
      implicit none

      type (struct_obs), intent(in) :: obsdat
      type (struct_obs), intent(inout) :: obs_sel
      real(8),           intent(in) :: hx(:,:)
      real(pre_obsReal), intent(in) :: zhamin
      real(pre_obsReal), intent(in) :: zhamax
      real(8),           intent(out):: hx_sel(:,:)
      integer, intent(in)    :: nens
      integer, intent(in)    :: nobsout

      integer :: iaccept,idata,ipnt,iwrite
      integer :: jdata,kobs,kobsout
      integer :: column_index
      integer :: active_index

      if(obsdat%mpi_local)then
         call obs_abort('obs_select() is not equipped to handle the case, ' // &
                        'mpi_local=.true.')
         return
      end if

      write(nobsout,'(1x,A,I7)')'stations prior to selection: ', obsdat%numHeader
      write(*,*) 'enter obs_select'

      kobsout=0 
      iwrite=0
      stations: do kobs=1,obsdat%numHeader
         ipnt  = obs_headElem_i(obsdat, OBS_RLN, kobs)
         idata = obs_headElem_i(obsdat, OBS_NLV, kobs)
         iaccept=0
         observations: do jdata = ipnt, ipnt + idata - 1
            ! To remove observations that are not in the desired vertical layer 
            if ((obs_bodyElem_r(obsdat, OBS_ZHA, jdata) > zhamin).and. &
                (obs_bodyElem_r(obsdat, OBS_ZHA, jdata) <= zhamax)) then
               iaccept=iaccept+1
               iwrite=iwrite+1
               do active_index=1,odc_numActiveColumn(obsdat%intBodies)
                  column_index=odc_columnIndexFromActiveIndex( &
                                      obsdat%intBodies%odc_flavour, active_index)
                  call obs_bodySet_i(obs_sel, column_index, iwrite, &
                       obs_bodyElem_i(obsdat, column_index, jdata))
               enddo
               do active_index=1,odc_numActiveColumn(obsdat%realBodies)
                  column_index=odc_columnIndexFromActiveIndex( &
                                     obsdat%realBodies%odc_flavour, active_index)
                  call obs_bodySet_r(obs_sel, column_index, iwrite, &
                       obs_bodyElem_r(obsdat, column_index, jdata))
               enddo
               hx_sel(1:nens,iwrite)=hx(1:nens,jdata)
            endif
         enddo observations

         ! adjust obs_sel%*Headers%columns
         if (iaccept > 0) then
            kobsout=kobsout+1
            do active_index=1,odc_numActiveColumn(obsdat%intHeaders)
               column_index=odc_columnIndexFromActiveIndex( &
                                     obsdat%intHeaders%odc_flavour, active_index)
               call obs_headSet_i(obs_sel, column_index, kobsout, &
                    obs_headElem_i(obsdat, column_index, kobs))
            enddo
            do active_index=1,odc_numActiveColumn(obsdat%realHeaders)
               column_index=odc_columnIndexFromActiveIndex( &
                                    obsdat%realHeaders%odc_flavour, active_index)
               call obs_headSet_r(obs_sel, column_index, kobsout, &
                    obs_headElem_r(obsdat, column_index, kobs))
            enddo
            obs_sel%cstnid(kobsout)=obsdat%cstnid(kobs)
            obs_sel%cfamily(kobsout)=obsdat%cfamily(kobs)
            call obs_headSet_i(obs_sel, OBS_NLV, kobsout, iaccept)
            call obs_headSet_i(obs_sel, OBS_RLN, kobsout, iwrite-iaccept+1)
         endif
      enddo stations
      obs_sel%mpi_local = .false.

      write(nobsout, '(1x, A, 2F11.6)') &
                 'after selection of observations in the range: ', zhamin, zhamax
      write(nobsout,'(1x,A,I7)') &
         'number of stations containing valid data      ',obs_sel%numHeader
      write(nobsout,'(1x,A,I7)') & 
         'number of observations now in the obsdat file ',obs_sel%numBody

   end subroutine obs_select


   subroutine obs_set_c(obsdat, name, row_index, value)
      !
      ! :Purpose: set a character(len=9) in the observation object
      !      To control access to the observation object.
      !
      implicit none
      type(struct_obs), intent(inout)  :: obsdat
      character(len=*), intent(in)     :: name
      integer         , intent(in)     :: row_index
      character(len=*), intent(in)     :: value

      select case (trim(name))
      case ('STID'); obsdat%cstnid (row_index) = value
         if(row_index == (obsdat%numHeader+1)) obsdat%numHeader = obsdat%numHeader+1

      case default
         call obs_abort('ERROR writing:  ' // trim(name) // &
                        ' is not a character(len=9) observation.')
         return
      end select
   end subroutine obs_set_c


   subroutine obs_set_current_body_list_from_family(obsdat, family, &
      list_is_empty, current_list)
      !
      ! :Purpose:
      !      Create a row_index list from the indicated family and place it in
      !      the body depot.
      !
      implicit none
      type(struct_obs), intent(inout), target :: obsdat
      character(len=*), intent(in) :: family
      logical, intent(out), optional :: list_is_empty
      type(struct_index_list), pointer, intent(out), optional :: current_list

      type(struct_index_list_depot), pointer :: depot
      type(struct_index_list), pointer :: index_list
      integer :: index_header, list, list_index, row_index
      integer :: first, last

      nullify(index_list)
      depot => obsdat%body_index_list_depot

      ! Search for an existing list
      if(present(current_list)) then
         if(associated(current_list)) then
            if (current_list%family == family) then
               index_list => current_list
            end if ! family matches
         end if ! associated

      else ! not present(current_list)
         do list = 1, NUMBER_OF_LISTS
            if (depot%index_lists(list)%family == family) then
               index_list => depot%index_lists(list)
               exit                     ! Don't look any further
            end if
         end do
      end if

      ! If the list does not already exist
      if (.not. associated(index_list)) then

         ! Acquire memory for the list
         if(present(current_list)) then
            ! This is an OMP thread. Re-use the same physical memory for the list
            index_list => ild_get_empty_index_list(depot, current_list)
         else
            index_list => ild_get_empty_index_list(depot)
         end if

         ! Initialize the new list
         index_list%family = family
         index_list%header = -99        ! not used

         !
         ! Populate the list
         !

         ! Loop over all header indices of the family
         list_index = 0
         call obs_set_current_header_list(obsdat, family)
         HEADER: do
            index_header = obs_getHeaderIndex(obsdat)
            if (index_header < 0) exit HEADER
            first= obs_headElem_i(obsdat,OBS_RLN,index_header)
            last = obs_headElem_i(obsdat,OBS_NLV,index_header) + first - 1
            do row_index=first,last     ! For each item indicated in the header
               ! Add the row_index to the list
               list_index = list_index + 1
               index_list%indices(list_index) = row_index
            end do
         end do HEADER
         index_list%indices(list_index+1)= -1 ! Flag the end of the list ...
         index_list%indices(list_index+2)= -1 ! ... clearly
      end if ! list does not already exist

      index_list%current_element = 0    ! Set pointer to the start of the list
      depot%current_list => index_list  ! Note the current list

      if(present(list_is_empty)) then
         ! Return whether the list is empty
         list_is_empty = (ild_get_next_index(depot, no_advance=.true.) < 0)
      end if

      if(present(current_list)) then
         ! Return a pointer to the current list
         current_list => index_list
      end if
   end subroutine obs_set_current_body_list_from_family


   subroutine obs_set_current_body_list_from_header(obsdat, header, &
      list_is_empty, current_list)
      !
      ! :Purpose:
      !      Create a row_index list from the indicated header and place it in
      !      the body depot.
      !
      implicit none
      type(struct_obs), intent(inout), target :: obsdat
      integer, intent(in) :: header
      logical, intent(out), optional :: list_is_empty
      type(struct_index_list), pointer, intent(out), optional :: current_list

      type(struct_index_list_depot), pointer :: depot
      type(struct_index_list), pointer :: index_list
      integer :: list, list_index, row_index
      integer :: first, last

      nullify(index_list)
      depot => obsdat%body_index_list_depot

      ! Search for an existing list
      if(present(current_list)) then
         if(associated(current_list)) then
            if (current_list%header == header) then
               index_list => current_list
            end if ! header matches
         end if ! associated

      else ! not present(current_list)
         do list = 1, NUMBER_OF_LISTS
            if (depot%index_lists(list)%header == header) then
               index_list => depot%index_lists(list)
               exit                     ! Don't look any further
            end if
         end do
      end if

      ! If the list does not already exist
      if (.not. associated(index_list)) then

         ! Acquire memory for the list
         if(present(current_list)) then
            ! This is an OMP thread. Re-use the same physical memory for the list
            index_list => ild_get_empty_index_list(depot, current_list)
         else
            index_list => ild_get_empty_index_list(depot)
         end if

         ! Initialize the new list
         index_list%family = 'xx'       ! not used
         index_list%header = header

         ! Populate the list
         first= obs_headElem_i(obsdat,OBS_RLN,header)
         last = obs_headElem_i(obsdat,OBS_NLV,header) + first - 1
         list_index = 0
         do row_index=first,last        ! For each item indicated in the header
            ! Add the row_index to the list
            list_index = list_index + 1
            index_list%indices(list_index) = row_index
         end do
         index_list%indices(list_index+1)= -1 ! Flag the end of the list ...
         index_list%indices(list_index+2)= -1 ! ... clearly
      end if ! list does not already exist

      index_list%current_element = 0    ! Set pointer to the start of the list
      depot%current_list => index_list  ! Note the current list

      if(present(list_is_empty)) then
         ! Return whether the list is empty
         list_is_empty = (ild_get_next_index(depot, no_advance=.true.) < 0)
      end if

      if(present(current_list)) then
         ! Return a pointer to the current list
         current_list => index_list
      end if
   end subroutine obs_set_current_body_list_from_header


   subroutine obs_set_current_body_list_all(obsdat, list_is_empty, current_list)
      !
      ! :Purpose:
      !      Create a row_index list containing all bodies and place it in the
      !      body depot.
      !
      implicit none
      type(struct_obs), intent(inout), target :: obsdat
      logical, intent(out), optional :: list_is_empty
      type(struct_index_list), pointer, intent(out), optional :: current_list

      type(struct_index_list_depot), pointer :: depot
      type(struct_index_list), pointer :: index_list
      integer :: list, list_index, row_index, index_header
      integer :: first, last

      nullify(index_list)
      depot => obsdat%body_index_list_depot

      ! Search for an existing list
      if(present(current_list)) then
         if(associated(current_list)) then
            if (      current_list%header == -1 &
                .and. current_list%family == '  ') then
               index_list => current_list
            end if ! null header and family
         end if ! associated

      else ! not present(current_list)
         do list = 1, NUMBER_OF_LISTS
            if (      depot%index_lists(list)%header == -1 &
                .and. depot%index_lists(list)%family == '  ') then
               index_list => depot%index_lists(list)
               exit                     ! Don't look any further
            end if
         end do
      end if

      ! If the list does not already exist
      if (.not. associated(index_list)) then

         ! Acquire memory for the list
         if(present(current_list)) then
            ! This is an OMP thread. Re-use the same physical memory for the list
            index_list => ild_get_empty_index_list(depot, current_list)
         else
            index_list => ild_get_empty_index_list(depot)
         end if

         ! Initialize the new list
         index_list%family = '  '       ! null
         index_list%header = -1         ! null

         !
         ! Populate the list
         !

         ! Loop over all header indices
         list_index = 0
         call obs_set_current_header_list(obsdat)
         HEADER: do
            index_header = obs_getHeaderIndex(obsdat)
            if (index_header < 0) exit HEADER
            first= obs_headElem_i(obsdat,OBS_RLN,index_header)
            last = obs_headElem_i(obsdat,OBS_NLV,index_header) + first - 1
            do row_index=first,last     ! For each item indicated in the header
               ! Add the row_index to the list
               list_index = list_index + 1
               index_list%indices(list_index) = row_index
            end do
         end do HEADER
         index_list%indices(list_index+1)= -1 ! Flag the end of the list ...
         index_list%indices(list_index+2)= -1 ! ... clearly
      end if ! list does not already exist

      index_list%current_element = 0    ! Set pointer to the start of the list
      depot%current_list => index_list  ! Note the current list

      if(present(list_is_empty)) then
         ! Return whether the list is empty
         list_is_empty = (ild_get_next_index(depot, no_advance=.true.) < 0)
      end if

      if(present(current_list)) then
         ! Return a pointer to the current list
         current_list => index_list
      end if
   end subroutine obs_set_current_body_list_all


   subroutine obs_set_current_header_list_from_family(obsdat, family)
      !
      ! :Purpose:
      !      Find or create a row_index list for the indicated family and place
      !      it in the header depot.
      !
      implicit none
      type(struct_obs), intent(inout), target :: obsdat
      character(len=*), intent(in) :: family

      type(struct_index_list_depot), pointer :: depot
      type(struct_index_list), pointer :: index_list
      integer :: list, list_index, row_index

      nullify(index_list)
      depot => obsdat%header_index_list_depot

      ! Search for an existing list
      do list = 1, NUMBER_OF_LISTS
         if (depot%index_lists(list)%family == family) then
            index_list => depot%index_lists(list)
            index_list%current_element=0! Start at the beginning of the list
            exit                        ! Don't look any further
         end if
      end do

      ! If the list does not already exist
      if (.not. associated(index_list)) then
         ! Create a new list
         index_list => ild_get_empty_index_list(depot)
         index_list%family = family
         index_list%header = -1

         ! Populate the list
         list_index = 0
         do row_index = 1, obsdat%numHeader
            ! If the station is of the right family
            if(obsdat%cfamily(row_index) == family) then
               ! Add the row_index to the list
               list_index = list_index + 1
               index_list%indices(list_index) = row_index
            end if
         end do
         index_list%indices(list_index+1)= -1 ! Flag the end of the list ...
         index_list%indices(list_index+2)= -1 ! ... clearly
      end if ! list does not already exist

      index_list%current_element = 0    ! Set pointer to the start of the list
      depot%current_list => index_list  ! Note the current list
   end subroutine obs_set_current_header_list_from_family


   subroutine obs_set_current_header_list_all(obsdat)
      !
      ! :Purpose:
      !      Find or create a row_index list for all headers and place it in the
      !      header depot.
      !
      implicit none
      type(struct_obs), intent(inout), target :: obsdat

      type(struct_index_list_depot), pointer :: depot
      type(struct_index_list), pointer :: index_list
      integer :: list, list_index, row_index

      nullify(index_list)
      depot => obsdat%header_index_list_depot

      ! Search for an existing list
      do list = 1, NUMBER_OF_LISTS
         if (depot%index_lists(list)%family == '  ') then
            index_list => depot%index_lists(list)
            index_list%current_element=0! Start at the beginning of the list
            exit                        ! Don't look any further
         end if
      end do

      ! If the list does not already exist
      if (.not. associated(index_list)) then
         ! Create a new list
         index_list => ild_get_empty_index_list(depot)
         index_list%family = '  '
         index_list%header = -1

         ! Populate the list
         list_index = 0
         do row_index = 1, obsdat%numHeader
            ! Add the row_index to the list
            list_index = list_index + 1
            index_list%indices(list_index) = row_index
         end do
         index_list%indices(list_index+1)= -1 ! Flag the end of the list ...
         index_list%indices(list_index+2)= -1 ! ... clearly
      end if ! list does not already exist

      index_list%current_element = 0    ! Set pointer to the start of the list
      depot%current_list => index_list  ! Note the current list
   end subroutine obs_set_current_header_list_all


   subroutine obs_setBodyPrimaryKey(obsdat,bodyIndex,primaryKey)
      !
      ! :Purpose:
      !      Set to the indicated value the body primary key for the indicated body.
      !
      implicit none
      type(struct_obs), intent(inout) :: obsdat
      integer(8),       intent(in)    :: primaryKey
      integer,          intent(in)    :: bodyIndex

      obsdat%bodyPrimaryKey(bodyIndex) = primaryKey

      if(bodyIndex == (obsdat%numBody+1)) then
         obsdat%numBody=obsdat%numBody+1
      endif

   end subroutine obs_setBodyPrimaryKey


   subroutine obs_setHeadPrimaryKey(obsdat,headerIndex,primaryKey)
      !
      ! :Purpose:
      !      Set to the indicated value the header primary key for the indicated header.
      !
      implicit none
      type(struct_obs), intent(inout) :: obsdat
      integer(8),       intent(in)    :: primaryKey
      integer,optional, intent(in)    :: headerIndex

      obsdat%headerPrimaryKey(headerIndex) = primaryKey
      if(headerIndex == (obsdat%numHeader+1)) then
         obsdat%numHeader=obsdat%numHeader+1
      endif

   end subroutine obs_setHeadPrimaryKey


   subroutine obs_setFamily(obsdat,Family_in,headerIndex_in,bodyIndex)
      !
      ! :Purpose:
      !      Set to the indicated value the family for the indicated header, or
      !      else for the indicated body.
      !
      implicit none
      type(struct_obs), intent(inout) :: obsdat
      character(len=*), intent(in)    :: Family_in
      integer,optional, intent(in)    :: headerIndex_in,bodyIndex

      integer          :: headerIndex

      if(present(headerIndex_in)) then
         headerIndex=headerIndex_in
      elseif(present(bodyIndex)) then
         headerIndex=obs_bodyElem_i(obsdat,OBS_HIND,bodyIndex)
      else
         call obs_abort('OBS_SETFAMILY: Header or Body index must be specified!')
         return
      endif

      obsdat%cfamily(headerIndex)=Family_in
      if(headerIndex == (obsdat%numHeader+1)) then
         obsdat%numHeader=obsdat%numHeader+1
      endif

   end subroutine obs_setFamily


   subroutine obs_sethind(obsSpaceData)
     !
     ! :PURPOSE: Set the header index in the body table
     !
     implicit none
     type(struct_obs) :: obsSpaceData
     integer :: ij,idata,idatend,bodyIndex,headerIndex
     !
     ! Set the header index in the body of obsSpaceData
     !
     ij=0
     do headerIndex = 1, obs_numheader(obsSpaceData)
        idata   = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
        idatend = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + idata - 1
        do bodyIndex= idata, idatend
           ij   = ij+1
           call obs_bodySet_i(obsSpaceData,OBS_HIND,IJ, headerIndex)
        end do
     end do
   end subroutine obs_sethind


   subroutine obs_status(obsdat, obs_full, numstns_out, numobs_out, kulout)
      !
      ! :Purpose: obtain basic status of the observation object
      !      Return the values of the object's status variables.
      !
      type (struct_obs), intent(in) :: obsdat
      logical, intent(out) :: obs_full
      integer, intent(out) :: numstns_out
      integer, intent(out) :: numobs_out
      integer, intent(in)  :: kulout


      ! PLH   if ( obsdat%numHeader >= nmxobs ) then
      ! PLH     if ( obsdat%numBody >= ndatamx .or. obsdat%numHeader >= nmxobs ) then
      if (     obsdat%numBody   >= obsdat%numBody_max &
          .or. obsdat%numHeader >= obsdat%numHeader_max) then
         write(kulout,*) ' OBSDAT FILE FULL'
         obs_full = .true.

      else
         obs_full = .false.
      end if

      numstns_out = obsdat%numHeader
      numobs_out  = obsdat%numBody
   end subroutine obs_status


   subroutine obs_tosqlbdy(obsdat,kobs,kulout)
      !
      ! :Purpose: print all data records associated with a station
      !
      ! :Arguments:
      !           :obsdat: obsSpaceData object
      !           :kobs: no. of observation
      !           :kulout: unit used for printing
      !
      implicit none

      type (struct_obs), intent(inout) :: obsdat
      integer,      intent(in) :: kobs
      integer,      intent(in) :: kulout

      integer  :: idata,idata2,ihaht,ihpht,ioer,ioma,ioma0,iomp,iomp6,ipnt
      integer  :: ippp, ivnm,ivnmc,istat,isigi, isigo,ivar,jdata,jtrans,var3d
      integer  :: obsVNM, obsASS, obsZHA, obsVCO, obsFLG
      integer  :: mrbcol,mrbcvt
      real     :: rppp
      character(len=100) :: message
      external :: mrbcol,mrbcvt

      ipnt  = obs_headElem_i(obsdat, OBS_RLN, kobs)
      idata = obs_headElem_i(obsdat, OBS_NLV, kobs)

      do jdata = ipnt, ipnt + idata - 1
         idata2 = jdata -ipnt + 1
         if (btest(obs_bodyElem_i(obsdat, OBS_FLG, jdata),12)) then
            var3d=1
         else
            var3d=0
         endif

         ippp=obs_bodyElem_r(obsdat, OBS_PPP, jdata)
         rppp=float(ippp)

         ivnm=obs_bodyElem_i(obsdat, OBS_VNM, jdata)
         istat=mrbcol(ivnm,ivnmc,1) 
         istat=mrbcvt(ivnmc,ivar ,obs_bodyElem_r(obsdat, OBS_VAR, jdata),1,1,1,1)
         istat=mrbcvt(ivnmc,iomp ,obs_bodyElem_r(obsdat, OBS_OMP, jdata),1,1,1,1)
         istat=mrbcvt(ivnmc,iomp6,obs_bodyElem_r(obsdat, OBS_OMP6,jdata),1,1,1,1)
         istat=mrbcvt(ivnmc,ioma ,obs_bodyElem_r(obsdat, OBS_OMA, jdata),1,1,1,1)
         istat=mrbcvt(ivnmc,ioma0,obs_bodyElem_r(obsdat, OBS_OMA0,jdata),1,1,1,1)
         istat=mrbcvt(ivnmc,ioer ,obs_bodyElem_r(obsdat, OBS_OER, jdata),1,1,1,1)
         istat=mrbcvt(ivnmc,ihpht,obs_bodyElem_r(obsdat, OBS_HPHT,jdata),1,1,1,1)
         istat=mrbcvt(ivnmc,ihaht,obs_bodyElem_r(obsdat, OBS_HAHT,jdata),1,1,1,1)
         istat=mrbcvt(ivnmc,isigi,obs_bodyElem_r(obsdat, OBS_SIGI,jdata),1,1,1,1)
         istat=mrbcvt(ivnmc,isigo,obs_bodyElem_r(obsdat, OBS_SIGO,jdata),1,1,1,1)
         jtrans=obs_bodyElem_i(obsdat, OBS_VCO, jdata)
         if (jtrans == 1) then
            istat=mrbcol(7001,ivnmc,1)
            istat=mrbcvt(ivnmc,ippp,rppp,1,1,1,1)
         elseif (jtrans == 2) then
            istat=mrbcol(7004,ivnmc,1)
            istat=mrbcvt(ivnmc,ippp,rppp,1,1,1,1)
         elseif (jtrans == 3) then
            istat=mrbcol(2150,ivnmc,1)
            istat=mrbcvt(ivnmc,ippp,rppp,1,1,1,1)
         else
            write(message,*) &
               'OBS_TOSQLBDY: attention, mauvaise coordonnee verticale, ', jtrans
            call obs_abort(message)
            return
         endif

         obsVNM = obs_bodyElem_i(obsdat, OBS_VNM, jdata)
         obsASS = obs_bodyElem_i(obsdat, OBS_ASS, jdata)
         obsZHA = obs_bodyElem_r(obsdat, OBS_ZHA, jdata)
         obsVCO = obs_bodyElem_i(obsdat, OBS_VCO, jdata)
         obsFLG = obs_bodyElem_i(obsdat, OBS_FLG, jdata)
         write(kulout,fmt=9201) kobs,idata2, obsVNM, ippp, obsASS, &
            ivar,iomp,iomp6,ioma,ioma0,ioer,ihpht,ihaht,isigi,isigo,var3d,  &
            obsZHA, obsVCO, obsFLG
      enddo

9201  format(1x,i9,',',i3,2(',',i6),',',i3,10(',',i8), &
         ',',i2,',',f10.3,',',i2,',',i12)

      return
   end subroutine obs_tosqlbdy


   subroutine obs_tosqlhdr(obsdata,kobs,kulout)
      !
      ! :Purpose: printing of the header of a station record for sql
      !
      !
      ! :Arguments:
      !          :obsdata: obsSpaceData object
      !          :kobs: no. of observation
      !          :kulout: unit used for output
      !
      implicit none

      type (struct_obs), intent(inout) :: obsdata
      integer, intent(in) :: kobs
      integer, intent(in) :: kulout

      integer :: ialt,idburp,ilon,ilat,iout
      integer :: obsDAT, obsETM, obsRLN, obsONM, obsINS, obsOTP
      integer :: obsNLV, obsPAS, obsREG, obsIP
      character(len=12) :: ccstnid
      real(8) :: torad 

      torad=4.d0*atan(1.d0)/180.d0

      ccstnid=obsdata%cstnid(kobs)

      ! Replace occasional appearance of "," by "b" in CCSTNID to avoid problem
      ! when converting this output to sqlite. - Xingxiu Deng, March 2009
      do
         iout=index(ccstnid,',')
         if (iout > 0 ) then
            ccstnid(iout:iout)='b'
         else
            exit
         endif
      enddo

      ialt=obs_headElem_r(obsdata, OBS_ALT, kobs)+400
      ilon=nint((obs_headElem_r(obsdata, OBS_LON, kobs)/torad)*100.0)
      ilat=nint((obs_headElem_r(obsdata, OBS_LAT, kobs)/torad+90.0)*100.0)

      idburp=mod(obs_headElem_i(obsdata, OBS_ITY, kobs),1000)
      obsDAT = obs_headElem_i(obsdata, OBS_DAT, kobs)
      obsETM = obs_headElem_i(obsdata, OBS_ETM, kobs)
      obsRLN = obs_headElem_i(obsdata, OBS_RLN, kobs)
      obsONM = obs_headElem_i(obsdata, OBS_ONM, kobs)
      obsINS = obs_headElem_i(obsdata, OBS_INS, kobs)
      obsOTP = obs_headElem_i(obsdata, OBS_OTP, kobs)
      obsNLV = obs_headElem_i(obsdata, OBS_NLV, kobs)
      obsPAS = obs_headElem_i(obsdata, OBS_PAS, kobs)
      obsREG = obs_headElem_i(obsdata, OBS_REG, kobs)
      obsIP  = obs_headElem_i(obsdata, OBS_IP , kobs)
      write(kulout,fmt=9200) kobs, CCSTNID, obsDAT, obsETM, &
           obsRLN, obsONM, obsINS, obsOTP, &
           idburp,ilat,ilon,ialt,              &
           obsNLV, obsPAS, obsREG, obsIP

9200  format(2x,i9,',',a9,',',i10,',',i8,',',i6,',',i6, &
         ',',i12,',',i6,4(',',i8),4(',',i6))

      return

   end subroutine obs_tosqlhdr


   subroutine obs_write(obsdat,hx, &
      nens,nobshdrout,nobsbdyout,nobshxout,nobsdimout)
      ! 
      ! :Purpose: 
      !      Write the obsdat info to unformatted files.
      !
      ! :Note: the body information is written in the order that it will
      !      be used by sekfeta.f
      !
      implicit none
      type(struct_obs), intent(in) :: obsdat
      real(8),      intent(in), dimension(:,:) :: hx
      integer,      intent(in) :: nens,nobshdrout,nobsbdyout, &
         nobshxout,nobsdimout

      integer :: irealBodies,jo,nrealBodies

      irealBodies=1
      do jo=1,obsdat%numHeader
         call obs_write_hdr(obsdat,jo,nobshdrout,irealBodies,nrealBodies)
         call obs_write_bdy(obsdat,jo,nobsbdyout)
         if (nens > 0) then
            call obs_write_hx(obsdat,hx,jo,nobshxout)
         endif
         irealBodies=irealBodies+nrealBodies
      enddo
      write(nobsdimout,*) obsdat%numHeader
      write(nobsdimout,*) irealBodies-1
      write(nobsdimout,*) nens

      return

   end subroutine obs_write


   subroutine obs_write_bdy(obsdat,kobs,kulout)
      !
      ! :Purpose: write the data records associated with a
      !                 station in unformatted form.
      !
      ! :Arguments:
      !           :obsdat: obsSpaceData object
      !           :kobs: no. of observation
      !           :kulout: unit used for writing 
      !
      implicit none 

      type(struct_obs), intent(in) :: obsdat
      integer, intent(in) ::  kobs,kulout

      integer :: ipnt,idata,j,jdata,k


      ipnt  = obs_headElem_i(obsdat, OBS_RLN, kobs) 
      idata = obs_headElem_i(obsdat, OBS_NLV, kobs)

      ! write the data records
      do jdata=ipnt,ipnt+idata-1
         write(kulout) &
           (obsdat%intBodies%columns(odc_ENKF_bdy_int_column_list(k) &
                                    )%value_i(jdata), &
                 k=1,size(odc_ENKF_bdy_int_column_list(:))),&
           (obsdat%realBodies%columns(odc_ENKF_bdy_real_column_list(j) &
                                     )%value_r(jdata), &
                 j=1,size(odc_ENKF_bdy_real_column_list(:)))
      enddo

      return

   end subroutine obs_write_bdy


   subroutine obs_write_hdr( obsdat, kobs, kulout, irealBodies, nrealBodies )
      !
      ! :Purpose: writing of the header of a station record
      !
      ! :Arguments:
      !           :obsdat: obsSpaceData object
      !           :kobs: no. of observation
      !           :kulout: unit used for output 
      !           :irealBodies: location in the sorted realBodies
      !           :nrealBodies: number of observations for this header
      !
      implicit none

      type(struct_obs), intent(in)  :: obsdat
      integer         , intent(in)  :: kobs
      integer         , intent(in)  :: kulout
      integer         , intent(in)  :: irealBodies
      integer         , intent(out) :: nrealBodies

      integer :: i,j,nprocs_mpi,ierr


      ! (note that as a part of the writing, the body is being sorted
      !  so that the order of the observations in the body array 
      !  corresponds with the order of the headers in the header array).
      call rpn_comm_size("GRID",nprocs_mpi,ierr)

      if(obsdat%mpi_local .and. nprocs_mpi>1) then
         call obs_abort('obs_write_hdr() is not equipped to handle the ' // &
                        'case, mpi_local=.true.')
         return
      end if

      nrealBodies=obs_headElem_i(obsdat, OBS_NLV, kobs)
      ! write the header's content 
      write(kulout) irealBodies, &
            (obs_headElem_i &
                  (obsdat, &
                   odc_columnIndexFromActiveIndex(obsdat%intHeaders%odc_flavour,&
                                                  i), &
                   kobs &
                  ), &
             i=2,odc_numActiveColumn(obsdat%intHeaders) &
            ), &

            (obs_headElem_r &
                 (obsdat, &
                  odc_columnIndexFromActiveIndex(obsdat%realHeaders%odc_flavour,&
                                                 j), &
                  kobs &
                 ), &
             j=1,odc_numActiveColumn(obsdat%realHeaders) &
            ), &

            obsdat%cstnid(kobs), &
            obsdat%cfamily(kobs)

      return

   end subroutine obs_write_hdr


   subroutine obs_write_hx(obsdat,hx,kobs,kulout)
      !
      ! :Purpose: write the interpolated values associated with a
      !                 station in unformatted form.
      !
      ! :Arguments:
      !           :obsdat: obsSpaceData object
      !           :hx: interpolated values    
      !           :kobs: no. of station 
      !           :kulout: unit used for writing 
      !
      implicit none 

      type(struct_obs), intent(in) :: obsdat
      real(8), intent(in), dimension(:,:) :: hx
      integer, intent(in) :: kobs
      integer, intent(in) :: kulout

      integer :: ipnt,idata,iens,jdata,nens

      nens = size(hx,1)

      ipnt  = obs_headElem_i(obsdat, OBS_RLN, kobs) 
      idata = obs_headElem_i(obsdat, OBS_NLV, kobs)

      ! write the data records
      do jdata=ipnt,ipnt+idata-1
         write(kulout) (hx(iens,jdata),iens=1,nens)
      enddo

      return

   end subroutine obs_write_hx


   function obs_famExist( obsdat, family, localMPI_opt )
     !
     ! :Purpose: Check if any observations of a particular family are present
     !          in obsdat. Returned result will be the MPI local value if
     !          the optional argument local_mpi is set to .true., will be
     !          the MPI global value otherwise.
     !
     ! :Arguments:
     !           :obsdat: ObsSpaceData structure
     !           :family: Obs family
     !           :localMPI_opt: return MPI local result; optional, default is .false.
     !
     implicit none 

     logical                                :: obs_famExist  ! Logical indicating if 'family' is part of the list
     type(struct_obs), intent(in)           :: obsdat
     character(len=2), intent(in)           :: family
     logical         , intent(in), optional :: localMPI_opt
     
     integer :: index_header,ierr
     logical :: famExist,local

     if (present(localMPI_opt)) then
        local = localMPI_opt
     else
        local = .false.
     end if
     
     famExist = .false.

     ! check if family exists in local MPI process
     do index_header = 1,obs_numheader(obsdat)
        if (obs_getFamily(obsdat,index_header) == family) then
           famExist = .true.
           exit
        end if
     end do

     if (local) then
        ! return MPI local value
        obs_famExist = famExist
     else
        ! return MPI global value
        call rpn_comm_allreduce(famExist,obs_famExist,1,"MPI_LOGICAL","MPI_LOR","GRID",ierr)
     end if

   end function obs_famExist

end module ObsSpaceData_mod
