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

module obsOperatorsChem_mod
  ! MODULE obsOperatorsChem_mod (prefix='oopc' category='5. Observation operators')
  !
  ! :Purpose: Observation operators for CH obs family, including nonlinear, tangent-linear
  !           and adjoint versions, and related setup and input routines.
  !
  use earthConstants_mod
  use mathPhysConstants_mod
  use obsSpaceData_mod
  use columnData_mod 
  use bufr_mod
  use physicsFunctions_mod
  use mpi_mod
  use utilities_mod
  use varNameList_mod
  use obsSubSpaceData_mod
  use obsFiles_mod
  use codtyp_mod
  use stateToColumn_mod   ! for use of s2c_column_hbilin
  use bmatrixchem_mod   

  implicit none
  save
  private

  ! public procedures
  public :: oopc_CHobsoperators, oopc_diagn_only, oopc_add_efftemp_obsfile

  !-------------------------------------------------------------------------
  ! Various structures and parameters for the CH family

  ! module structures
  ! -----------------

  type :: struct_chm_obsoperators  
  
     !  Structure holding work variables for observation operators
     !     
     !  Variable               Description
     !  --------               -----------
     !  nobslev                Number of observations in the profile
     !  nmodlev                Number of model levels in the column
     !  varno                  BUFR descriptor element for obs units
     !  constituent_id         BUFR code element of local GRIB Table 08046 identifying the constituent
     !                         (similar to BUFR Table 08043)
     !  modelIndex             Obs operator index
     !                         0 - vertical interpolator
     !                         1 - layer averaging
     !                         2 - layer integration
     !  layer_identified       .true. if a layer (with identified layer boundaries)
     !                         .false if layer boundaries are not available.
     !  vmodpress              Model layer boundaries taken as middle between model level
     !  vlayertop              Layer top (final work values in Pa)
     !  vlayerbottom           Layer bottom (final work values in Pa)
     !  vweights               Second order Lagrangian interp integration weights
     !  zh                     Initial innovation model array (other than conversion constants)
     !  zhp                    Part of innovation operator not related to resolution.
     !  imodlev_top            Top level of non-zero values in zh
     !  imodlev_bot            Bottom level of non-zero values in zh
     !  trial                  Trial (background) profile at observation location
     !  tt                     Temperature profile on model levels (Kelvin)
     !  hu                     Specific humidity 
     !  height                 Height on model levels (m)
     !  pp                     Pressure on model levels (Pa)
     !  lat                    Latitude of observation (radians)
     !  lon                    Longitude of observation (radians)
     !  obslev                 Observation profile level values (OBS_PPP)
     !  varName                Variable/obs nomvar
     !  stnid                  Observation station ID
     !  date                   YYYYMMDD (date of obs)
     !  hhmm                   HHMM (time of obs)
     !  obs_index              Observation index
     !                         Note: Depending on the data of interest, the index of a required array element or 
     !                               profile associated to an observation can be identified from (lat,long,date,hhmm,
     !                               stnid,optional task-dependent identifier if needed) or obs_index. 
     !                               The latter is for associations of data identified within
     !                               processing of individual CPUs. Each of the two index identifiers is represented 
     !                               by the unique character string identifier 'code' of struct_oss_obsdata
     !                               (e.g. see obsdata_get_header_code_r for use of (lat,long,date,hhmm,stnid)).
     !  vco                    Index of vertical coord type for obs
     !                           1 - Altitudes (m)
     !                           2 - Pressure (Pa)
     !                           3 - Channel index
     !                           4 - not provided with obs. Obs is for total column values.
     !                           5 - not provided with obs. Obs is a surface point value.
     !  iavgkern               Integer indicating if averaging kernels are to be applied. Value
     !                         of zero indicates no averaging kernel to be applied. Non-zero value
     !                         indicates index in chm_avgkern%obsSubSpace arrays.
     !  apply_genoper          Indicates if the generalized observation operator should be applied
     !  column_bound           Boudary imporsed on a column measurement
     !  dtransform             Derivative for any transform that needs to be applied to a profile
     
     integer :: nobslev,nmodlev,modelIndex,constituent_id,vco,varno,date,hhmm,iavgkern,obs_index
     logical :: layer_identified,apply_genoper
     real(8) :: lat,lon,column_bound
     character(len=12) :: stnid
     character(len=4)  :: varName
     real(8), allocatable :: vlayertop(:),vlayerbottom(:),vmodpress(:),tt(:),height(:),pp(:)
     real(8), allocatable :: zh(:,:),zhp(:,:),vweights(:,:),obslev(:),dtransform(:),hu(:)
     real(8), pointer     :: trial(:)
     integer, allocatable :: imodlev_top(:),imodlev_bot(:)

  end type struct_chm_obsoperators

  type :: struct_chm_info
     !  Information arrays retrieved from auxiliary file regarding vertical levels 
     !  or averaging kernels
     !
     !  Variable               Description
     !  --------               -----------
     !  n_stnid                Number of sub-families (identified via STNIDs)
     !  stnids                 Sub-families (STNIDs; * are wild cards)
     !  element                BUFR element in data block
     !  source                 0: Set entirely from the auxiliary file being read. No 
     !                            initial values read from observation files
     !                         1: Initial values in observation files for constant number
     !                            of vertical levels (may be adjusted after input)
     !  vco                    Vertical coordinate type (1, 2, or 3, see bufr_read_mod)
     !
     !  ibegin                 Position index of start of data for given
     !                         sub-family.
     !  n_lvl                  Number of vertical levels 
     !  n_lat                  Number of latitudes
     !  lat                    Latitudes (degrees; ordered in increasing size)
     !
     !  vlayertop              Layer top 
     !  vlayerbottom           Layer bottom
     !  rak                    Averaging kernel matrices
     
     integer ::  n_stnid
     character(len=12), allocatable :: stnids(:)
     integer, allocatable :: element(:),source(:)
     integer, allocatable :: vco(:),n_lat(:)
     integer, allocatable :: ibegin(:),n_lvl(:)
     real(8), allocatable :: rak(:),vlayertop(:),vlayerbottom(:)
     real(8), allocatable :: lat(:)
 
     type(struct_oss_obsdata), allocatable :: obsSubSpace(:)
 
  end type struct_chm_info

  type :: struct_chm_griddata

     !  Structure storing gridded fields 
     !     
     !  Variable               Description
     !  --------               -----------
     !  field2d                Gridded 2 field
     !  field3d                Gridded 3 field (lon,lat,vlev)
     !  nlat                   number of latitudes
     !  nlon                   number of longitudes
     !  nlev                   number of vertical levels
     !  lat,lon                grid lat,lon in radians
     !  vlev                   vertical levels
     !  ivkind                 Index of vertical coordinate type. Defintion may vary according to source.
     !                         For fields read for RPN files and use of convip:
     !                             0: P is in height [m] (metres) with respect to sea level 
     !                             1: P is in sigma [sg] (0.0 -> 1.0) 
     !                             2: P is in pressure [mb] (millibars) 
     !                             3: P is in an arbitrary code 
     !                             4: P is in height [M] (metres) with respect to ground level 
     !                             5: P is in hybrid coordinates [hy] 
     !                             6: P is in theta [th] 
     !                         For use with obs                      
     
     real(8), pointer :: field2d(:,:),field3d(:,:,:),lat(:),lon(:),vlev(:)
     integer :: nlev,nlon,nlat,ivkind
  
  end type struct_chm_griddata

  type(struct_chm_info) :: chm_layers
  type(struct_chm_info) :: chm_avgkern

  ! Arrays for integration upper boundary of retrieved total column measurements 
  type(struct_oss_obsdata) :: chm_column_boundary

  ! Arrays for background error std. dev.
  type(struct_oss_obsdata) :: chm_sigma_trial

  ! File name of auxiliary text file constaining supplemental observation information
  character(len=50), parameter :: chm_aux_filename="obsinfo_chm" 

  ! Max nummber of constituents (max size of related arrays)
  integer, parameter :: chm_constituents_size=30   ! = max allowed value of "iconstituent_id" for Table 08046.
                                                   ! Value to be increased as needed up to a max of 6999 as values
                                                   ! > 7000 (and less 0) are assumed assigned to non-constituent fields  

  ! Arrays containing input reference fields and fields interpolated 
  ! to obs locations

  type(struct_oss_obsdata)  :: chm_ref_trial
  type(struct_chm_griddata) :: chm_ref_fields(0:chm_constituents_size,2)

  ! Arrays to contain the calculated concentration-weighted effective temperature
  ! associated to total column data. It will be stored in the observation file.

  type(struct_oss_obsdata) :: chm_efftemp

  ! General config/setup information parameters 
  ! See description list of NAMCHEM namelist parameters in routine chm_setup

  integer, parameter :: assim_maxfamnum=1           ! Could be used for other families as well with >1
  integer, parameter :: assim_maxsize=100           ! max size of assim_* arrays  

  integer :: assim_famNum
  logical :: assim_all(assim_maxfamnum)
  integer :: assim_num(assim_maxfamnum),assim_varno(assim_maxfamnum,assim_maxsize),assim_nlev(assim_maxfamnum,assim_maxsize)
  integer :: assim_exclude_nflag(assim_maxfamnum),assim_exclude_flag(assim_maxfamnum,assim_maxsize)
  character(len=9) :: assim_stnid(assim_maxfamnum,assim_maxsize)
  character(len=2) :: assim_fam(assim_maxfamnum)
  
  integer :: chm_generalized_operator(0:chm_constituents_size)   ! Same as genoper in NAMCHEM
  integer :: chm_tropo_mode(0:chm_constituents_size),chm_tropo_bound(0:chm_constituents_size)
  integer :: chm_obsdata_maxsize
  real(8) :: chm_tropo_column_top(0:chm_constituents_size)
  ! Identification of the model 
  character(len=10) :: modelName = 'GEM-MACH' 

  ! Setup initialization key
  logical :: initializedChem = .false.   
  
  
  !--------------------------------------------------------------------------
  
contains

  !============================== CH obs family =============================
   
  !--------------------------------------------------------------------------
  ! oopc_setupCH
  !--------------------------------------------------------------------------
  subroutine oopc_setupCH(datestamp_opt)
    !
    !:Purpose: To set up additional information required by constituent obs and
    !          not provided in obsSpaceData.  Also to assign observation layer
    !          top and bottom levels (and averaging kernel matrices).
    !          See 'oopc_CHobsoperators'. 
    !
    implicit none

    ! Arguments:
    integer, intent(in), optional :: datestamp_opt

    write(*,*) 'Begin oopc_setupCH'
    !allocate(chm_layers%chm_obsSubx(2))

    ! Read NAMCHEM namelist and set related parameters

    call chm_read_namchem
      
    ! Read top and bottom layer boundaries of partial (or total) column meausurements
  
    call chm_read_layers
      
    ! To deallocate space if required elsewhere, one should use
    ! call chm_dealloc_layers
   
    ! Read averaging kernel matrices
  
    call chm_read_avgkern
  
    ! To deallocate space if required elsewhere, one should use
    ! call chm_dealloc_avgkern
  
    ! Read reference (e.g. climatological) fields
  
    call chm_read_ref_fields(datestamp_opt=datestamp_opt)

    ! Allocation of chm_efftemp done in chm_setup instead of obsdata_add_data1d
    ! to ensure allocation is done for all processors, including those without associated data.
    ! This is to ensure that rpn_comm_allgather will work in routine obsdata_MPIGather.

    if (.not.associated(chm_efftemp%data1d)) then
      call oss_obsdata_alloc(chm_efftemp,chm_obsdata_maxsize,dim1=1)
      chm_efftemp%nrep=0
    end if

    write(*,*) 'Completed oopc_setupCH'
  
  end subroutine oopc_setupCH

  !--------------------------------------------------------------------------
  ! chm_read_namchem
  !--------------------------------------------------------------------------
  subroutine chm_read_namchem
    !:Purpose: Read and store miscellaneous flags and constants.
    !
    !:Comment: assim_* arrays could instead be made available to all families
    !          by moving them to a different input namelist (and changing its
    !          dimensions settings).
    !
    !:Output:
    !
    !  :Read from NAMCHEM namelist:
    ! 
    !     :genoper:
    !                           Whether generalized observation operator should
    !                           be used and selection of approach
    !                             ===  =======================================
    !                             <=0  not applied
    !                              1   use trial field xb for mass weighted
    !                                  increment distribution
    !                              2   use a combination of the difference of
    !                                  an external reference xc and the trial
    !                                  field xb, i.e. mass weighted increment
    !                                  distribution as a(xc-xb) + b*xc where a
    !                                  and b depend on the size of
    !                                  sum[(xc-xb)/sig(xb)]^2 over the profile
    !                             ===  =======================================
    !
    !     :assim_fam:           List of families to which filt_diagn_only is to
    !                           apply.
    ! 
    !     :assim_exclude_flag:  Array specifying bits for identifying
    !                           diagnostic-only observations for observations
    !                           that would otherwise be assimilated according to
    !                           the other assim_* arrays
    ! 
    !     :assim_exclude_nflag: Number of bit flags to specify in
    !                           assim_exclude_flag array
    ! 
    !     :assim_all:           Logical indicating if all assimilatable obs of
    !                           the specified family will be assimilated
    !                           (default is .true.)
    ! 
    !     :assim_num:           Number combinations (stnid, bufr element,
    !                           multi/uni-level) identified for assimilation.
    !                           All others will not be assimilated. OmP and OmA
    !                           diagnostics and output will still be produced
    !                           for non-assimilated datasets.
    !                             ===  =======================================
    !                              0   none are to be assimilated
    !                             >0   sets of (stnid, bufr varno,
    !                                  multi/uni-levels) to be assimilated
    !                             ===  =======================================
    ! 
    !     :assim_varno:         Bufr elements of obs sets for assimilation. A
    !                           value of 0 implies that all are to be used.
    !
    !     :assim_stnid:         Stnids of obs sets for assimilation. '*' denote
    !                           wild cards
    ! 
    !     :assim_nlev:            ===  =========================
    !                              0   multi-level and uni-level
    !                              1   uni_level
    !                             >1   multi-level 
    !                             ===  =========================
    ! 
    !     :tropo_mode:          Integer indicating if special treatment is to be
    !                           given to the troposphere when assimilating total
    !                           column measurements. Values indicate
    !                             ===  =======================================
    !                              0   No special treatment given (default)
    !                              1   Values of the adjoint model above
    !                                  obsoper%column_bound set to zero. If
    !                                  specified, generalized innovation
    !                                  operator only applied below
    !                                  obsoper%column_bound in the tangent
    !                                  linear model.
    !                              2   Values of tangent linear model and
    !                                  adjoint model above obsoper%column_bound
    !                                  set to zero.
    !                             ===  =======================================
    !                           Array index refers to BUFR code element of Table
    !                           08046 (iconstituent_id) identifying the
    !                           constituent. Relevant for total column
    !                           measurements only.
    ! 
    !     :tropo_bound:         Integer indicating which column top value to use
    !                           if tropo_mode is non-zero.
    !                             ===  =======================================
    !                              0   Use fixed value of tropo_column_top
    !                              1   Use model determination of tropopause
    !                              2   Use model determination of PBL
    !                             ===  =======================================
    !                           Options 1 and 2 will default to the value set
    !                           in tropo_column_top if the model derived column
    !                           top could not be determined. Relevant for total
    !                           column measurements only.
    !                          
    !     :tropo_column_top:    Default value to use for the column boundary
    !                           (in Pa). Array index refers to BUFR code element
    !                           of Table 08046 (iconstituent_id) identifying
    !                           the constituent. Relevant for total column
    !                           measurements only.
    ! 
    !     :obsdata_maxsize:     Max allowed size of work arrays (in terms of
    !                           number of obs) associated to ordered observation
    !                           indices
    ! 
    !      :modelName:          Nome/identifiier of forecast model
    !                           Default: 'GEM-MACH'
    !                           Set to 'GEM' for varNames of 'O3L', 'CH4L', and 'N2OL'
    !
    implicit none

    ! Locals:
    integer :: FNOM, FCLOS
    integer :: IERR, ios, nulnam, i

    integer :: genoper(0:chm_constituents_size)
    integer :: tropo_mode(0:chm_constituents_size),tropo_bound(0:chm_constituents_size)
    integer :: obsdata_maxsize
    real(8) :: tropo_column_top(0:chm_constituents_size)
    
    character(len=10)  :: namfile 

    EXTERNAL FNOM,FCLOS

    namelist /namchem/ assim_fam,assim_all,assim_num,assim_stnid,assim_varno,    &
                     assim_nlev, assim_exclude_nflag,assim_exclude_flag,       &
                     tropo_mode,tropo_bound,tropo_column_top,obsdata_maxsize,  &
                     genoper,modelName
  
    ! Default NAMCHEM values

    genoper(:)=0
    chm_obsdata_maxsize=90000
 
    assim_fam(:)=''
    assim_fam(1)='CH'
    assim_famNum=1
    assim_all(:)=.true.
    assim_num(:)=0  
    assim_stnid(:,:)='*********'
    assim_varno(:,:)=0
    assim_nlev(:,:)=0
 
    assim_exclude_nflag(:)=1
    assim_exclude_flag(:,1)=6  ! this is for the 'in reserve' bit for a BURP marker
    assim_exclude_flag(:,2:)=0

    tropo_mode(:) = 0
    tropo_bound(:) = 0
    tropo_column_top(:) = 0.0

    ! Read from namelist file NAMCHEM

    namfile=trim("flnml")
    nulnam=0
    ierr=FNOM(nulnam,namfile,'R/O',0)

    read(nulnam,nml=namchem,iostat=ios)
    if (ios.lt.-4.or.ios.gt.0) then 
      call utl_abort('chm_read_namchem: Error in reading NAMCHEM namelist. iostat = ' // trim(utl_str(ios)) )   
    else if (mpi_myid.eq.0) then
      write(*,nml=namchem)   
    end if
  
    ierr=FCLOS(nulnam)      

    do i=size(assim_fam),1,-1
      if (assim_fam(i).ne.'') exit
    end do
    assim_famnum=i

    chm_generalized_operator(:) = genoper(:)
    chm_tropo_mode(:) = tropo_mode(:)
    chm_tropo_bound(:) = tropo_bound(:)
    chm_tropo_column_top(:) = tropo_column_top(:)
  
  end subroutine chm_read_namchem

  !-------------------------------------------------------------------------------  
  !--------------------- Routines related to layer top & bottom levels------------

  !--------------------------------------------------------------------------
  ! chm_read_layers
  !--------------------------------------------------------------------------
  subroutine chm_read_layers
    !
    !:Purpose: To read and to store top- and bottom-layer boundaries for CH
    !          sub-families
    !
    !:Comments:
    !
    !  A) The option of reading from observation files is TBD. This will change
    !     the approach in allocating the arrays size as the sizes will become
    !     dependent on the number of related obs for which the observation files
    !     will need to be read.
    !
    implicit none

    ! Locals:
    integer :: FNOM, FCLOS
    integer :: IERR, JLEV, JELM, nulstat, ios, isize, icount
    logical :: LnewExists
  
    character (len=128) :: ligne

    EXTERNAL FNOM,FCLOS
  
    ! Initialization

    chm_layers%n_stnid=0

    INQUIRE(FILE=trim(chm_aux_filename),EXIST=LnewExists)
    IF (.not.LnewExists )then
      WRITE(*,*)   '----------------------------------------------'
      WRITE(*,*)   'chm_read_layers: COULD NOT FIND AUXILIARY FILE ' // trim(chm_aux_filename)
      WRITE(*,*)   '----------------------------------------------'
      return
    ENDIF

    ! Check for available layer info.

    NULSTAT=0
    IERR=FNOM(NULSTAT,trim(chm_aux_filename),'SEQ',0)
    IF ( IERR .EQ. 0 ) THEN
      open(unit=nulstat, file=trim(chm_aux_filename), status='OLD')
    ELSE
      CALL utl_abort('chm_read_layers: COULD NOT OPEN AUXILIARY FILE ' // trim(chm_aux_filename))
    ENDIF

    ios=0
    read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne
    do while (trim(adjustl(ligne(1:13))).ne.'SECTION II:') 
      read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne
    end do    
  
    ! Read number of observation set sub-families (STNIDs and ...) and allocate space
   
    read(nulstat,*,iostat=ios,err=10,end=10) chm_layers%n_stnid
    read(nulstat,*,iostat=ios,err=10,end=10) isize

    allocate(chm_layers%stnids(chm_layers%n_stnid))
    allocate(chm_layers%vco(chm_layers%n_stnid))
    allocate(chm_layers%source(chm_layers%n_stnid),chm_layers%ibegin(chm_layers%n_stnid))
    allocate(chm_layers%element(chm_layers%n_stnid),chm_layers%n_lvl(chm_layers%n_stnid))
    allocate(chm_layers%vlayertop(isize),chm_layers%vlayerbottom(isize))
 
    chm_layers%element(:)=0
    chm_layers%vco(:)=0
    chm_layers%source(:)=0
    chm_layers%n_lvl(:)=1

    ! Begin reading for each sub-family
    ! Important: Combination of STNID, BUFR element and number of vertical levels
    !            to determine association to the observations.

    icount=0
    do jelm=1,chm_layers%n_stnid
      chm_layers%ibegin(jelm)=icount+1

      read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne

      ! Read STNID (* is a wildcard)
    
      read(nulstat,'(2X,A9)',iostat=ios,err=10,end=10) chm_layers%stnids(jelm) 

      ! Read (1) Obs BUFR element.
      !      (2) Vertical coord type (1, 2, or 3)
      !     (3) Flag indication if EOR provided from this auxiliary file or
      !         to be read from the observation file,
      !     (4) Number of vertical levels 
      !
      ! Important: Combination of STNID, BUFR element and number of vertical levels  
      !            to determine association to the observations.

      read(nulstat,*,iostat=ios,err=10,end=10) chm_layers%element(jelm),chm_layers%vco(jelm),  &
        chm_layers%source(jelm),chm_layers%n_lvl(jelm)  
    
      if (icount+chm_layers%n_lvl(jelm).gt.isize) CALL utl_abort('chm_read_layers: READING PROBLEM. Max array size exceeded: ' &
                                                                // trim(utl_str(isize)))    


      read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne
    
      if (chm_layers%n_lvl(jelm).ge.1) then   
        do jlev=1,chm_layers%n_lvl(jelm)
          icount=icount+1
          
          ! Read top and bottom levels
          
          read(nulstat,*,iostat=ios,err=10,end=10)                 &
                 chm_layers%vlayertop(icount),chm_layers%vlayerbottom(icount)
        end do
      end if

      ! if (chm_layers%source(jelm).eq.1) then
      !    
      !   Read from observation files
      ! 
      !      .....
      !
      ! end if

    end do
   
 10 if (ios.gt.0) CALL utl_abort('chm_read_layers: READING PROBLEM. File read error message number: ' // trim(utl_str(ios)))    
 
 11 CLOSE(UNIT=NULSTAT)
    IERR=FCLOS(NULSTAT)    

  end subroutine chm_read_layers

  !--------------------------------------------------------------------------
  ! chm_get_layer_boundaries
  !--------------------------------------------------------------------------
  subroutine chm_get_layer_boundaries(cstnid,varno,ivco,nlev,default_top, &
                                      default_bottom,lfound,layertop, &
                                      layerbottom)
    !
    !:Purpose: To return layer boundaries for an observation. Combination of
    !          STNID, element variable number and number of vertical levels to
    !          determine association to the observations. Default values for top
    !          and bottom layers for total column measurements are to be
    !          provided.
    !
    !:Arguments:        
    !    - ivco            type of vertical coordinate (see burpread_mod.ftn90
    !                      or routine chm_obsoperators for definitions)
    implicit none

    ! Arguments:
    character(len=12), intent(in) :: cstnid ! station id
    integer, intent(in)           :: varno  ! BUFR element
    integer, intent(in)           :: ivco
    integer, intent(in)           :: nlev   ! number of levels in the observation
    real(8), intent(in)           :: default_top ! default value for top layer for total column measurement
    real(8), intent(in)           :: default_bottom ! default value for bottom layer for total column measurement
    logical, intent(inout)        :: lfound ! .true. if layer boundaries found
    real(8), intent(out)          :: layertop(nlev) ! top layer values
    real(8), intent(out)          :: layerbottom(nlev) ! bottom layer values

    ! Locals:
    integer                       :: ISTNID,JN,start_index
    logical                       :: iset

    ! Find stnid with same number of vertical levels, and same BUFR element
          
    ISTNID=0
    lfound=.false.

    DO JN=1,chm_layers%n_stnid

       ! First compare STNID values allowing for * and blanks in 
       ! chm_layers%stnids(JN) as wildcards
       iset = utl_stnid_equal(chm_layers%stnids(JN),CSTNID)

       ! Check if number of levels, code, and vertical coordinate type are equal.
       ! If number of levels is one and no vertical coordinate provided for total column measurement (i.e. IVCO.eq.4),
       ! then check of vertical coordinate type is disregarded afterwards.
       IF (iset) THEN
          IF ( varno.EQ.chm_layers%element(JN) .AND. NLEV.EQ.chm_layers%n_lvl(JN) .AND. &
              (IVCO.EQ.chm_layers%vco(JN).OR.IVCO.EQ.4) ) THEN
             ISTNID=JN
             exit
          END IF
       END IF
       
    END DO

    IF (ISTNID.EQ.0) THEN
       ! If integrated layer information not found, if a total column measurement
       ! set to defaults, else do nothing

       if (bufr_IsIntegral(varno) .and. nlev.eq.1) then          
          lfound=.true.
          layertop(1) = default_top
          layerbottom(1) = default_bottom
       end if

    ELSE
       ! layer information has been found in auxiliary file
       lfound=.true.
       start_index = chm_layers%ibegin(ISTNID)
       layertop(:) = chm_layers%vlayertop(start_index:start_index+nlev-1)
       layerbottom(:) = chm_layers%vlayerbottom(start_index:start_index+nlev-1)  
    END IF

  end subroutine chm_get_layer_boundaries

  subroutine chm_dealloc_layers
    !
    !:Purpose: To deallocate temporary storage space used for layer info
    !

    implicit none

    if (chm_layers%n_stnid.eq.0) return

    call chm_dealloc_info(chm_layers)
 
  end subroutine chm_dealloc_layers

  !-------------------------------------------------------------------------------
  !------------------- Routines related to averaging kernel matrices -------------

  !--------------------------------------------------------------------------
  ! chm_read_avgkern
  !--------------------------------------------------------------------------
  subroutine chm_read_avgkern
    !
    !:Purpose: To read averaging kernels from auxiliary file or observation file
    !

    implicit none

    ! Locals:
    integer, parameter :: ndim=2

    integer :: istnid

    ! read the averaging kernel information from the auxiliary file
    call chm_read_avgkern_auxfile

    ! set size of observation file array
    allocate(chm_avgkern%obsSubSpace(chm_avgkern%n_stnid))

    ! read from observation file
    do istnid=1,chm_avgkern%n_stnid
       if (chm_avgkern%source(istnid) == 1) then
          
          ! retrieve data from stats blocks (with bkstp=14 and block_type='DATA')
          chm_avgkern%obsSubSpace(istnid) = obsf_obsSub_read('CH',chm_avgkern%stnids(istnid),bufr_avgkern, &
                                     chm_avgkern%n_lvl(istnid), ndim, bkstp_opt=14, &
                                     block_opt='DATA', match_nlev_opt=.true.)
          
       end if
    end do

  end subroutine chm_read_avgkern

  !--------------------------------------------------------------------------
  ! chm_read_avgkern_auxfile
  !--------------------------------------------------------------------------
  subroutine chm_read_avgkern_auxfile
    !
    !:Purpose: To read and to store averaging kernel matricesfor CH sub-families
    !
    !:Comments:
    !      - Currently implemented for only one latitude band
    !

    implicit none

    !Locals:
    integer :: FNOM, FCLOS
    integer :: IERR, JLEV, JELM, nulstat, ios, isize, icount, iend  
    logical :: LnewExists
  
    character (len=128) :: ligne

    EXTERNAL FNOM,FCLOS

    ! Initialization

    chm_avgkern%n_stnid=0

    INQUIRE(FILE=trim(chm_aux_filename),EXIST=LnewExists)
    IF (.not.LnewExists )then
      WRITE(*,*)   '--------------------------------------------------------'
      WRITE(*,*)   'chm_read_avgkern_auxfile: COULD NOT FIND AUXILIARY FILE ' // trim(chm_aux_filename)
      WRITE(*,*)   '--------------------------------------------------------'
      return
    ENDIF

    ! Check for available layer info.

    NULSTAT=0
    IERR=FNOM(NULSTAT,trim(chm_aux_filename),'SEQ',0)
    IF ( IERR .EQ. 0 ) THEN
      open(unit=nulstat, file=trim(chm_aux_filename), status='OLD')
    ELSE
      CALL utl_abort('chm_read_avgkern_auxfile: COULD NOT OPEN AUXILIARY FILE ' // trim(chm_aux_filename))
    ENDIF

    ios=0
    read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne
    do while (trim(adjustl(ligne(1:14))).ne.'SECTION III:') 
      read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne
    end do    
  
    ! Read number of observation set sub-families (STNIDs and ...) and allocate space
   
    read(nulstat,*,iostat=ios,err=10,end=10) chm_avgkern%n_stnid
    read(nulstat,*,iostat=ios,err=10,end=10) isize

    allocate(chm_avgkern%stnids(chm_avgkern%n_stnid))
    allocate(chm_avgkern%source(chm_avgkern%n_stnid),chm_avgkern%ibegin(chm_avgkern%n_stnid))
    allocate(chm_avgkern%element(chm_avgkern%n_stnid),chm_avgkern%n_lvl(chm_avgkern%n_stnid))
    allocate(chm_avgkern%rak(isize))
 
    chm_avgkern%element(:)=0
    chm_avgkern%source(:)=0
    chm_avgkern%n_lvl(:)=1

    ! Begin reading for each sub-family
    ! Important: Combination of STNID, BUFR element and number of vertical levels 
    !            to determine association to the observations.

    icount=1
    STNIDLOOP: do jelm=1,chm_avgkern%n_stnid
      chm_avgkern%ibegin(jelm)=icount

      read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne

      ! Read STNID (* is a wildcard)
    
      read(nulstat,'(2X,A9)',iostat=ios,err=10,end=10) chm_avgkern%stnids(jelm) 

      ! Read (1) Obs BUFR element.
      !      (2) Flag indication if avgkern provided from this auxiliary file or
      !          to be read from an observation file,
      !      (3) Number of vertical levels
      !
      ! Important: Combination of STNID, BUFR element and number of vertical levels
      !            to determine association to the observations.

      read(nulstat,*,iostat=ios,err=10,end=10) chm_avgkern%element(jelm),  &
        chm_avgkern%source(jelm),chm_avgkern%n_lvl(jelm)  
    
      if (icount+chm_avgkern%n_lvl(jelm).gt.isize) &
        CALL utl_abort('chm_read_avgkern_auxfile: READING PROBLEM.Max array size exceeded:' // trim(utl_str(isize)))    


      read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne

      ! disregard data section if values to be specified in BUFR file
      if (chm_avgkern%source(jelm).eq.1) cycle STNIDLOOP
    
      if (chm_avgkern%n_lvl(jelm).gt.1) then   
        do jlev=1,chm_avgkern%n_lvl(jelm)

          iend=icount+chm_avgkern%n_lvl(jelm)-1

          ! Read averaging kernel matrix   
          read(nulstat,*,iostat=ios,err=10,end=10) chm_avgkern%rak(icount:iend)

          icount=iend+1

        end do
      end if

    end do STNIDLOOP
   
 10 if (ios.gt.0) CALL utl_abort('chm_read_avgkern_auxfile: READING PROBLEM. File read error message number: ' // trim(utl_str(ios)))    
 
 11 CLOSE(UNIT=NULSTAT)
    IERR=FCLOS(NULSTAT)    

  end subroutine chm_read_avgkern_auxfile

  !--------------------------------------------------------------------------
  ! chm_dealloc_avgkern
  !--------------------------------------------------------------------------
  subroutine chm_dealloc_avgkern
    !
    !:Purpose: To deallocate temporary storage space used for averaging kernels
    !
    implicit none

    ! Locals:
    integer :: istnid

    if (chm_avgkern%n_stnid.eq.0) return

    if (allocated(chm_avgkern%obsSubSpace)) then
       do istnid=1,chm_avgkern%n_stnid
          if (chm_avgkern%source(istnid).eq.1) call oss_obsdata_dealloc(chm_avgkern%obsSubSpace(istnid))
       end do
       deallocate(chm_avgkern%obsSubSpace)
    end if

    call chm_dealloc_info(chm_avgkern)
  
  end subroutine chm_dealloc_avgkern

  !--------------------------------------------------------------------------
  ! chm_find_avgkern
  !--------------------------------------------------------------------------
  function chm_find_avgkern(cstnid,varno,nlev) result(ISTNID)
    !
    !:Purpose: To find the averaging kernel for an observation if one is
    !          specified. Returns 0 if either not found or not specified.
    !          Combination of STNID, BUFR element and number of vertical levels
    !          to determine association to the observations.
    !
    implicit none
    integer :: ISTNID ! Index of averaging kernel in chm_avgkern if found. Zero indicates  averaging kernel not found.

    ! Arguments:
    character(len=12), intent(in) :: cstnid ! station id
    integer, intent(in) :: varno ! BUFR descriptor element
    integer, intent(in) :: nlev  ! number of levels in the observation

    ! Locals:
    integer :: JN
    logical :: iset

    ! Find stnid with same number of vertical levels, and same BUFR element
          
    ISTNID=0

    DO JN=1,chm_avgkern%n_stnid

       ! First compare STNID values allowing for * and blanks in 
       ! chm_avgkern%stnids(JN) as wildcards
       iset = utl_stnid_equal(chm_avgkern%stnids(JN),CSTNID)

       ! Check if number of levels and BUFR code are equal.
       IF (iset) THEN
          IF ( varno.EQ.chm_avgkern%element(JN) .AND. NLEV.EQ.chm_avgkern%n_lvl(JN) ) THEN
             ISTNID=JN
             exit
          END IF
       END IF
       
    END DO

  end function chm_find_avgkern

  !--------------------------------------------------------------------------
  ! chm_get_avgkern
  !--------------------------------------------------------------------------
  subroutine chm_get_avgkern(istnid,stnid,nlev,zlat,zlon,idate,itime,avg_kern)
    !
    !:Purpose: To return averaging kernel for an observation.
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: istnid ! index of averaging kernel in chm_avgkern
    character(len=*), intent(in) :: stnid
    integer, intent(in)  :: nlev   ! number of observation levels
    real(8), intent(in)  :: zlat,zlon
    integer, intent(in)  :: idate  ! YYYYMMDD
    integer, intent(in)  :: itime  ! HHMM
    real(8), intent(out) :: avg_kern(nlev,nlev) ! the averaging kernel

    ! Locals:
    integer :: start_index,end_index

    if (istnid.gt.0 .and. istnid.le.chm_avgkern%n_stnid) then
       
       if (chm_avgkern%source(istnid).eq.0) then
          ! get averaging kernel from auxiliary file
          start_index = chm_avgkern%ibegin(ISTNID)
          end_index = nlev*(start_index+nlev-1)
          avg_kern = RESHAPE(chm_avgkern%rak(start_index:end_index),(/nlev,nlev/),ORDER =(/2,1/))
       else
          ! get averaging kernel from observation file
          avg_kern = oss_obsdata_get_array2d(chm_avgkern%obsSubSpace(istnid), oss_obsdata_get_header_code(zlon,zlat,idate,itime,stnid))
       end if

    else
       call utl_abort("chm_get_avgkern: Invalid station ID index.")
    end if

  end subroutine chm_get_avgkern

  !-------------------------------------------------------------------------------
  !------------------------- Routines related to gridded reference fields  -------

  !--------------------------------------------------------------------------
  ! chm_read_ref_fields
  !--------------------------------------------------------------------------
  subroutine chm_read_ref_fields(datestamp_opt)
    !
    !:Purpose:  To read reference fields as directed by the content of the
    !           auxiliary file.
    !
    ! Comments:
    !      - ****** NOT TESTED *********
    !      - Fields are provided in RPN/fst files specified in the auxiliary
    !        file (with path and filename)
    !      - Reference fields can be in a separate RPN file with name provided
    !        by the auxiliary or in monthly static background stats file
    !        (glbchemcov or bgcov; see 'isrc' below).
    !      - Fields assumed to be of the same units as those of the
    !        corresponding input trial fields
    !
    implicit none

    ! Arguments:
    integer, intent(in), optional :: datestamp_opt

    ! Locals:
    character(len=128) :: fname
    character(len=4) :: varName
    character(len=12) :: etiket
    integer :: i,id,nd,j,ndim,ijour,imonth,iday,itime,isrc
    real(8) :: day
    integer, external :: newdate
   
    integer, external :: FNOM, FCLOS
    integer :: IERR, nulstat, ios
    logical :: LExists
    
    logical, parameter :: linterp=.true.
    
    integer :: ni, nj, nkeys, kind
    real(8), allocatable :: array1(:,:,:),array2(:,:,:),lvls(:),xlat(:),xlong(:) ! Allocated in chm_fst_read
  
    character (len=128) :: ligne

    ! Initialize dimensions to zero

    chm_ref_fields(:,:)%nlon=0
    chm_ref_fields(:,:)%nlat=0
    chm_ref_fields(:,:)%nlev=1
    
    inquire(FILE=trim(chm_aux_filename),EXIST=LExists)
    IF (.not.LExists )then
      WRITE(*,*)   '---------------------------------------------------'
      WRITE(*,*)   'chm_read_ref_fields: COULD NOT FIND AUXILIARY FILE ' // trim(chm_aux_filename)
      WRITE(*,*)   '---------------------------------------------------'
      return
    ENDIF

    ! Check for file names containing ref fields

    NULSTAT=0
    IERR=FNOM(NULSTAT,trim(chm_aux_filename),'SEQ',0)
    IF ( IERR .EQ. 0 ) THEN
       open(unit=nulstat, file=trim(chm_aux_filename), status='OLD')
    ELSE
       CALL utl_abort('chm_read_ref_fields: COULD NOT OPEN AUXILIARY FILE ' // trim(chm_aux_filename))
    ENDIF

    ios=0
    read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne
    do while (trim(adjustl(ligne(1:14))).ne.'SECTION IV:') 
       read(nulstat,'(A)',iostat=ios,err=10,end=11) ligne
    end do    
    
    ! Read number of constituents with associated input file(s)
   
    read(nulstat,*,iostat=ios,err=10,end=10) ndim
    if (ndim.le.0) go to 10
    
    ! Initialization

    if (linterp.and.present(datestamp_opt)) then
       ierr = newdate(datestamp_opt,ijour,itime,-3)
       if (ierr<0) then
          Write(*,*) "Invalid datestamp ",datestamp_opt,ijour,itime,ierr
          call utl_abort('chm_read_ref_fields')
       endif
       imonth = MOD(ijour/100,100)
       iday = MOD(ijour,100)
       day=iday+itime*1.0D-8
       if (day.gt.15.) then
          day=day-15.0
       else
          day=day+15.0
       end if
    endif
    
    ! Get needed fields for each file(s)

    do i=1,ndim 

       ! Read id,nd,isrc. id: constituent code; nd: number of sets; 1 or 2;
       ! isrc: 1 for fname being in the auxiliary file, 0 for glbchemcov (or bgcov if glbchemcov not present) 
       
       read(nulstat,*,iostat=ios,err=10,end=10)
       read(nulstat,*,iostat=ios,err=10,end=10) id,nd,isrc    
       varName=vnl_varnameFromVarnum(0,id)

       if (isrc.eq.1) then
          read(nulstat,*,iostat=ios,err=10,end=10) fname
       else
         inquire(file='./glbchemcov',exist=LExists)
         if (LExists) then
            fname='./glbchemcov'
         else
            inquire(file='./bgcov',exist=LExists)
            if(LExists) then
              fname='./bgcov'
            else               
               call utl_abort('chm_read_ref_fields: did not find file.')
            end if
          end if
       end if
       
       do j=1,nd
          read(nulstat,*,iostat=ios,err=10,end=10) etiket             
                           
          call utl_readFstField(trim(fname),varName,-1,imonth,-1,etiket,ni,nj,nkeys,array1,xlat_opt=xlat,xlong_opt=xlong,lvls_opt=lvls,kind_opt=kind)

          if (j.eq.1) then
              chm_ref_fields(id,1)%nlon=ni
              chm_ref_fields(id,1)%nlat=nj
              chm_ref_fields(id,1)%nlev=nkeys
              chm_ref_fields(id,1)%ivkind=kind   
                         
              allocate(chm_ref_fields(id,1)%field3d(ni,nj,nkeys))
              allocate(chm_ref_fields(id,1)%vlev(nkeys),chm_ref_fields(id,1)%lon(ni),chm_ref_fields(id,1)%lat(nj))
              
              chm_ref_fields(id,1)%lat(1:nj)=xlat(1:nj)*MPC_RADIANS_PER_DEGREE_R8
              chm_ref_fields(id,1)%lon(1:ni)=xlong(1:ni)*MPC_RADIANS_PER_DEGREE_R8
              where (chm_ref_fields(id,1)%lon(1:ni).lt.0.0) chm_ref_fields(id,1)%lon(1:ni)=2.0*MPC_PI_R8 + &
                                                            chm_ref_fields(id,1)%lon(1:ni)
              chm_ref_fields(id,1)%vlev(1:nkeys)=lvls(1:nkeys)              
          else
              chm_ref_fields(id,2)%nlon=ni
              chm_ref_fields(id,2)%nlat=nj
              chm_ref_fields(id,2)%nlev=nkeys
              chm_ref_fields(id,2)%ivkind=kind

              allocate(chm_ref_fields(id,2)%field3d(ni,nj,nkeys))
              allocate(chm_ref_fields(id,2)%vlev(nkeys),chm_ref_fields(id,2)%lon(ni),chm_ref_fields(id,2)%lat(nj))
              
              chm_ref_fields(id,2)%lat(1:nj)=xlat(1:nj)*MPC_RADIANS_PER_DEGREE_R8
              chm_ref_fields(id,2)%lon(1:ni)=xlong(1:ni)*MPC_RADIANS_PER_DEGREE_R8
              where (chm_ref_fields(id,2)%lon(1:ni).lt.0.0) chm_ref_fields(id,2)%lon(1:ni)=2.0*MPC_PI_R8 + &
                                                            chm_ref_fields(id,2)%lon(1:ni)
              chm_ref_fields(id,2)%vlev(1:nkeys)=lvls(1:nkeys)
          end if

          if (.not.linterp .or. (.not.present(datestamp_opt))) then

             if (j.eq.1) then 
                 chm_ref_fields(id,1)%field3d(:,:,:) = array1(:,:,:)
             else
                 chm_ref_fields(id,2)%field3d(:,:,:) = array1(:,:,:)
             end if

          else

             ! Following for interpolation as a function of days from mid-months.
             
             if (iday.gt.15) then
                 if (imonth.eq.12) then
                    call utl_readFstField(trim(fname),varName,-1,1,-1,etiket,ni,nj,nkeys,array2,lvls_opt=lvls,kind_opt=kind)
                else
                   call utl_readFstField(trim(fname),varName,-1,imonth+1,-1,etiket,ni,nj,nkeys,array2,lvls_opt=lvls,kind_opt=kind)
                end if
          
                ! Linearly interpolate in time (approximately - assumes 30 day months)

                if (j.eq.1) then              
                   chm_ref_fields(id,1)%field3d(:,:,:) = (array1(:,:,:)*(30.0-day)+array2(:,:,:)*day)/30.0
                else
                   chm_ref_fields(id,2)%field3d(:,:,:) = (array1(:,:,:)*(30.0-day)+array2(:,:,:)*day)/30.0
                end if
             
             else if (iday.le.15) then
                if (imonth.eq.1) then
                   call utl_readFstField(trim(fname),varName,-1,12,-1,etiket,ni,nj,nkeys,array2,lvls_opt=lvls,kind_opt=kind)
                else
                   call utl_readFstField(trim(fname),varName,-1,imonth-1,-1,etiket,ni,nj,nkeys,array2,lvls_opt=lvls,kind_opt=kind)
                end if

                ! Linearly interpolate in time (approximately - assumes 30 day months)

                if (j.eq.1) then
                   chm_ref_fields(id,1)%field3d(:,:,:) = (array2(:,:,:)*(30.0-day)+array1(:,:,:)*day)/30.0
                else
                   chm_ref_fields(id,2)%field3d(:,:,:) = (array2(:,:,:)*(30.0-day)+array1(:,:,:)*day)/30.0
                end if
             
             end if
          
          end if
 
          if (allocated(array1)) deallocate(array1)
          if (allocated(array2)) deallocate(array2)   
                 
       end do
    end do 
     
 10 if (ios.gt.0) CALL utl_abort('chm_read_ref_fields: READING PROBLEM. File read error message number: ' // trim(utl_str(ios)))    
 
 11 CLOSE(UNIT=NULSTAT)
    IERR=FCLOS(NULSTAT)    
    
  end subroutine chm_read_ref_fields

  !--------------------------------------------------------------------------
  ! chm_set_reference_obsdata
  !--------------------------------------------------------------------------
  subroutine chm_set_reference_obsdata(obsoper)
    !
    !:Purpose: To determine and to store reference profile at obs location if
    !          needed by the observation operators.
    !
    !:Input:
    !
    !    :obsoper%constituent_id: Constituent id
    !    :obsoper%nmodlev:        Number of model levels for variables other
    !                             than uu and vv
    !    :obsoper%pressmod:       Model pressure array
    !    :obsoper%tt:             Model temperature (Kelvin)
    !    :obsoper%height:         Model height (m)
    !    :obsoper%hu:             Specific humidity 
    !    :obsoper%lat:            Latitude (rad)
    !    :obsoper%lon:            Longitude (rad)
    !
    !:Output:
    ! 
    !    :chm_ref_trial:          Reference profile object
    !
    !:Comments:
    !      - ***** NOT TESTED *****
    implicit none

    ! Arguments
    type(struct_chm_obsoperators), intent(inout) :: obsoper

    ! Locals
    integer :: i,istart,id
    real(8) :: tropo_press, refprof(obsoper%nmodlev),refprof2(obsoper%nmodlev),dt
    real(8), allocatable :: pressrefin(:)
    logical, allocatable :: lsuccess(:)
    
    if (obsoper%constituent_id < 0 .or. obsoper%constituent_id > chm_constituents_size) return

    id=obsoper%constituent_id
    
    if (chm_generalized_operator(id).le.1) return           
    if (chm_ref_fields(id,1)%nlat.eq.0) return
    
    ! Set vertical levels of reference.
    ! Convert to pressure coordinate if needed.
    
    if (allocated(pressrefin)) deallocate(pressrefin)
    allocate(pressrefin(chm_ref_fields(id,1)%nlev))
    pressrefin(:)=chm_ref_fields(id,1)%vlev(1:chm_ref_fields(id,1)%nlev)

    if (allocated(lsuccess)) deallocate(lsuccess)
    allocate(lsuccess(chm_ref_fields(id,1)%nlev))
    lsuccess(:)=.true.
    
    if (chm_ref_fields(id,1)%ivkind.eq.2) then
        pressrefin(:)=pressrefin(:)*100. ! Conversion from hPa to Pa.
    else if (chm_ref_fields(id,1)%ivkind.eq.0) then
        where (pressrefin.lt.obsoper%height(obsoper%nmodlev)) pressrefin=obsoper%height(obsoper%nmodlev)
        pressrefin(:) = phf_convert_z_to_pressure(pressrefin,obsoper%height,obsoper%pp, &
                        chm_ref_fields(id,1)%nlev,obsoper%nmodlev,obsoper%lat,lsuccess)
    else if (chm_ref_fields(id,1)%ivkind.eq.4) then
        pressrefin(:)=pressrefin(:) + obsoper%height(obsoper%nmodlev)
        pressrefin(:) = phf_convert_z_to_pressure(pressrefin,obsoper%height,obsoper%pp, &
                        chm_ref_fields(id,1)%nlev,obsoper%nmodlev,obsoper%lat,lsuccess)
    else if (chm_ref_fields(id,1)%ivkind.eq.1) then
        pressrefin(:)=pressrefin(:)*obsoper%pp(obsoper%nmodlev) ! Convert from sigma to Pa   
    else
       call utl_abort('chm_get_reference_obsdata: Cannot handle vertical coordinate of kind ' // trim(utl_str(chm_ref_fields(id,1)%ivkind)))
    end if
    
    ! Interpolate to obs lat/long location and model level

    call s2c_column_hbilin(chm_ref_fields(id,1)%field3d,pressrefin, &
                    chm_ref_fields(id,1)%nlon,chm_ref_fields(id,1)%nlat,chm_ref_fields(id,1)%nlev, &
                    chm_ref_fields(id,1)%lon,chm_ref_fields(id,1)%lat,obsoper%lon,obsoper%lat, &
                    refprof,obsoper%pp,obsoper%nmodlev)
    
    if (chm_ref_fields(id,2)%nlat.gt.0.and.chm_ref_fields(id,2)%nlon.gt.0.and.chm_ref_fields(id,2)%nlev.gt.0) then
        
        if (any(obsoper%tt.le.0.0)) call utl_abort('chm_get_reference_obsdata: Missing TT for determining tropopause pressure')
        
        ! Get second reference field (for troposphere)
        
        tropo_press=-1.0
        
        if (all(obsoper%hu.ge.0.0D0)) then
           tropo_press=phf_get_tropopause(obsoper%nmodlev,obsoper%pp,obsoper%tt,obsoper%height,hu_opt=obsoper%hu)
        else
           tropo_press=phf_get_tropopause(obsoper%nmodlev,obsoper%pp,obsoper%tt,obsoper%height)
         end if

        if (tropo_press.gt.0) then
            
           ! Set vertical levels of reference.
           ! Convert to pressure coordinate if needed
 
           if (allocated(pressrefin)) deallocate(pressrefin)
           allocate(pressrefin(chm_ref_fields(id,2)%nlev))    
           pressrefin(:)=chm_ref_fields(id,2)%vlev(1:chm_ref_fields(id,2)%nlev)

           if (allocated(lsuccess)) deallocate(lsuccess)
           allocate(lsuccess(chm_ref_fields(id,2)%nlev))
           lsuccess(:)=.true.

           if (chm_ref_fields(id,2)%ivkind.eq.2) then
               pressrefin(:)=pressrefin(:)*100. ! Conversion from hPa to Pa.
           else if (chm_ref_fields(id,2)%ivkind.eq.0) then
               where (pressrefin.lt.obsoper%height(obsoper%nmodlev)) pressrefin=obsoper%height(obsoper%nmodlev)
               pressrefin(:) = phf_convert_z_to_pressure(pressrefin,obsoper%height,obsoper%pp, &
                               chm_ref_fields(id,2)%nlev,obsoper%nmodlev,obsoper%lat,lsuccess)
           else if (chm_ref_fields(id,2)%ivkind.eq.4) then
               pressrefin(:)=pressrefin(:) + obsoper%height(obsoper%nmodlev)
               pressrefin(:) = phf_convert_z_to_pressure(pressrefin,obsoper%height,obsoper%pp, &
                               chm_ref_fields(id,2)%nlev,obsoper%nmodlev,obsoper%lat,lsuccess)
           else if (chm_ref_fields(id,2)%ivkind.eq.1) then
               pressrefin(:)=pressrefin(:)*obsoper%pp(obsoper%nmodlev) ! Convert from sigma to Pa        
           else 
               call utl_abort('chm_get_reference_obsdata: Cannot handle vertical coordinate of kind ' // trim(utl_str(chm_ref_fields(id,2)%ivkind)))
           end if
      
           ! Interpolate to obs lat/long and model levels

           call s2c_column_hbilin(chm_ref_fields(id,2)%field3d,pressrefin, &
                    chm_ref_fields(id,2)%nlon,chm_ref_fields(id,2)%nlat, &
                    chm_ref_fields(id,2)%nlev,chm_ref_fields(id,2)%lon, &
                    chm_ref_fields(id,2)%lat,obsoper%lon,obsoper%lat,refprof2,obsoper%pp,obsoper%nmodlev)
    
        end if

        ! Combine with upper level profile
       
        do i=obsoper%nmodlev,3,-1
           if (obsoper%pp(i).lt.tropo_press) exit
           refprof(i)=refprof2(i)            
        end do
        istart=i
             
        ! Apply linear combination of four levels just above the tropopause
        
        do i=istart,max(2,istart-3),-1
            dt=(istart+1.0-i)/5.0
            refprof(i)=dt*refprof2(i) + (1.0-dt)*refprof(i)
        end do
                    
    end if 

    if (allocated(pressrefin)) deallocate(pressrefin)
    if (allocated(lsuccess)) deallocate(lsuccess) 

    ! ------- Save in chm_ref_trial ---------
       
    if (.not.associated(chm_ref_trial%data1d)) then
       call oss_obsdata_alloc(chm_ref_trial, chm_obsdata_maxsize, dim1=obsoper%nmodlev)
       chm_ref_trial%nrep = 0
    end if

    ! Here, nrep will count the number of filled elements in the data arrays
    chm_ref_trial%nrep = chm_ref_trial%nrep+1 

    if (chm_ref_trial%nrep.gt.chm_obsdata_maxsize) &
         call utl_abort('chm_get_ref_obsdata: Reach max size of array ' // trim(utl_str(chm_obsdata_maxsize)) )
  
    ! The obsoper%obs_index (header index) serves as the unique locator code 
    write(chm_ref_trial%code(chm_ref_trial%nrep),'(I22)') obsoper%obs_index

    ! Save profile in chm_ref_trial
    
    chm_ref_trial%data1d(:,chm_ref_trial%nrep) = refprof(:)

  end subroutine chm_set_reference_obsdata
 
  !--------------------------------------------------------------------------
  ! chm_get_ref_column
  !--------------------------------------------------------------------------
  function chm_get_ref_column(code) result(array)
    !
    !:Purpose: To extract and to provide column from chm_ref_field associated to
    !          code.     
    !  
    implicit none
    real(8) :: array(chm_ref_trial%dim1) ! retrieved array from obsdata%data1d of dimension obsdata%dim1

    ! Arguments
    character(len=*), intent(in) :: code ! unique identifying code

    ! Locals:
    integer :: stat ! search success (0 - found; 1 = no data; 2 = not found)

    array = oss_obsdata_get_array1d(chm_ref_trial,code,stat)
    if (stat.gt.0) call utl_abort("chm_get_ref_column: Code not found - " // &
                                  trim(code))
    
  end function chm_get_ref_column

  !-----------------------------------------Misc ---------------------------------
  !-------------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! chm_dealoc_info
  !--------------------------------------------------------------------------
  subroutine chm_dealloc_info(info)
    !
    !:Purpose: To deallocate struct_chm_info instance
    !
    
    implicit none

    type(struct_chm_info), intent(inout) :: info

    if (allocated(info%stnids))       deallocate(info%stnids)
    if (allocated(info%element))      deallocate(info%element)
    if (allocated(info%source))       deallocate(info%source)
    if (allocated(info%vco))          deallocate(info%vco)
    if (allocated(info%n_lat))        deallocate(info%n_lat)
    if (allocated(info%ibegin))       deallocate(info%ibegin)
    if (allocated(info%n_lvl))        deallocate(info%n_lvl)
    if (allocated(info%rak))          deallocate(info%rak)
    if (allocated(info%vlayertop))    deallocate(info%vlayertop)
    if (allocated(info%vlayerbottom)) deallocate(info%vlayerbottom)
    if (allocated(info%lat))          deallocate(info%lat)

  end subroutine chm_dealloc_info

  !--------------------------------------------------------------------------
  ! oopc_diag_only
  !--------------------------------------------------------------------------
  logical function oopc_diagn_only(cfamName,cstnid,varno,nobslev,flag)
    ! 
    !:Purpose: To identify whether or not the obs set identified by the
    !          combination of (cstnidin,bufrin,nlevs) will be assimilated or
    !          else used for independent verifications after
    !          assimilation/minimization
    ! 
    implicit none

    ! Arguments
    character(len=*), intent(in) :: cfamName ! Family name
    character(len=*), intent(in) :: cstnid   ! Input station id
    integer, intent(in) :: varno   ! Obs BUFR number
    integer, intent(in) :: nobslev ! Number of levels
    integer, intent(in) :: flag    ! observation integer flag

    ! Locals:
    integer :: i,elemId,ifam
    
    ifam=0
    if (assim_famNum.gt.0) then
       do i=1,assim_famNum
          if (assim_fam(i).eq.cfamName) then
             ifam=i
             exit
          end if
       end do
    end if

    if (ifam.eq.0) then
       ! assimilate all observations
       oopc_diagn_only = .false.
       return    
    end if
    
    if (assim_all(ifam)) then
       ! assimilate all observations
       oopc_diagn_only = .false.
    else if (assim_num(ifam).le.0) then
       ! assimilate no observations
       oopc_diagn_only = .true.
    else if (assim_num(ifam).gt.0) then
       ! check if this observation is listed in the assim_* arrays
       elemId=0
       do i=1,assim_num(ifam)
          if (utl_stnid_equal(trim(assim_stnid(ifam,i)),trim(cstnid))) then
             if (assim_varno(ifam,i).eq.0.or.assim_varno(ifam,i).eq.varno) then
                if (assim_nlev(ifam,i).eq.0.or.(nobslev.eq.1.and.assim_nlev(ifam,i).eq.1).or. &
                     (nobslev.gt.1.and.assim_nlev(ifam,i).gt.1)) then
                   elemId=i
                   exit
                end if
             end if
          end if
       end do
       oopc_diagn_only = elemId.eq.0
    end if
    
    if (oopc_diagn_only) return

    ! check if the observation integer flag has a bit marked by assim_exclude_flag (same flagging as in filt_suprep)
    if (assim_exclude_nflag(ifam).gt.0) then
       do i=1,assim_exclude_nflag(ifam)
          if (btest(flag, 13 - assim_exclude_flag(ifam,i) )) then
             oopc_diagn_only = .true.
             return
          end if
       end do
    end if

  end function oopc_diagn_only

  !-------------------------------------------------------------------------------
  !------------------ Routines associated to chm_efftemp -------------------------

  !--------------------------------------------------------------------------
  ! oopc_add_efftemp_obsdata
  !--------------------------------------------------------------------------
  subroutine oopc_add_efftemp_obsdata(code,temp_eff)
    !
    !:Purpose: To add effective temperature value to its obsdata object
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: code ! unique identifying code
    real(8), intent(in) :: temp_eff(:)   ! effective temperature

    call oss_obsdata_add_data1d(chm_efftemp,temp_eff,code,chm_obsdata_maxsize)
    
  end subroutine oopc_add_efftemp_obsdata

  !--------------------------------------------------------------------------
  ! oopc_add_efftemp_obsfile
  !--------------------------------------------------------------------------
  subroutine oopc_add_efftemp_obsfile()
    !          
    !:Purpose: To add effective temperatures in obs file.
    !
    implicit none

    ! Locals:
    integer :: nrep_modified,varno(1)

    ! If needed, add effective temperature values in obs file for total column measurements

    !if ( .not.obsf_fileTypeIsBurp() ) call utl_abort('oopc_add_efftemp_obsfile: only compatible with BURP files')

    call oss_obsdata_MPIallgather(chm_efftemp)
    
    if (chm_efftemp%nrep.gt.0) then
        varno(1)=12001
        nrep_modified = obsf_obsSub_update(chm_efftemp,'CH',varno(1:max(1,chm_efftemp%dim2)),bkstp_opt=0, &
                                        block_opt='INFO',multi_opt='UNI') 
        write(*,*) 'oopc_add_efftemp_obsfile: Added ',nrep_modified,' effective temperature values in the obs file.'
    end if 

  end subroutine oopc_add_efftemp_obsfile

  !------------------ CONTROL ROUTINES FOR OBSERVATION OPERATORS -----------------
  !-------------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! oopc_CHobsoperators
  !--------------------------------------------------------------------------
  subroutine oopc_CHobsoperators(columnTrl, obsSpaceData, kmode, &
                                 columnAnlInc_opt, jobs_opt, destObsColumn_opt)
    !
    !:Purpose: To apply the observation operators for chemical constituents.
    !          Mode of operator set by kmode.
    !
    !:Comments:
    !      - See type struct_chm_obsoperators for description of obsoper elements.
    !      - Currently can only handle the case when nlev_bkgrnd == nlev_inc
    !
    !:Arguments:
    !   :columnTrl:  Column of x_background interpolated to observation
    !                location. Can have the same vertical levels as the
    !                trial field (columnTrlOnTrlLev) or as the increment field
    !                (columnTrlOnAnlIncLev)
    !   :kmode:
    !        +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
    !        |kmode|       Mode of         |             Results               |
    !        |     | Observation Operator  |                                   |
    !        +=====+=======================+===================================+
    !        |  0  |for general simulation |OmP and total Jo(x_background)     |
    !        |     |operator               |for CH. OmP saved in OBS_OMP of    |
    !        |     |(non-linear and linear)|obsSpaceData                       |
    !        +-----+-----------------------+-----------------------------------+
    !        |  1  |for identification or  |background error standard dev.     |
    !        |     |determination of       |in observation space saved in      |
    !        |     |of sqrt(diag(H*B*H^T). |in OBS_HPHT of obsSpaceData        |
    !        |     |Depends on the presence|if OmP error std dev not initially
    !        |     |of OmP error std dev   |available in OBS_OMPE              |
   !        +-----+-----------------------+-----------------------------------+
    !        |  2  |for tangent linear     |Hdx saved in OBS_WORK of           |
    !        |     |operator               |obsSpaceData                       |
    !        +-----+-----------------------+-----------------------------------+
    !        |  3  |for adjoint of tangent |H^T * R^-1 (OmP-Hdx) in            |
    !        |     |linear operator        |columnAnlInc_opt                   |
    !        +-----+-----------------------+-----------------------------------+    
    !
    !   :columnAnlInc_opt: Optional argument for input/output of column of
    !                      increment (column). For kmode=2, used as input for
    !                      increment H_horiz dx interpolated to observation
    !                      location. For kmode=3, used as output for H^T * R^-1
    !                      (OmP-Hdx). Required for kmode=2,3.
    !
    !
    !   :jobs_opt:       Optional output of total Jo(x_background) for chemical
    !                    constituents. Required for kmode=0 and not provided
    !                    otherwise.
    !
    !


    ! More Comments (not rendered in Sphinx):
    !      - Two equivalent methods for looping over a report body.
    !
    !        Method 1:
    !
    !             call obs_set_current_body_list(obsSpaceData,headerIndex)
    !             BODY: do
    !
    !                bodyIndex = obs_getBodyIndex(obsSpaceData)
    !                if (bodyIndex < 0) exit BODY1
    !
    !                ... obs_bodyElem_r(obsSpaceData, ... ,bodyIndex)
    !  
    !             enddo BODY
    !
    !        Method 2:
    !
    !             bodyIndex_start = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
    !             bodyIndex_end = bodyIndex_start+obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex)-1
    !             do  bodyIndex=bodyIndex_start,bodyIndex_end
    !                ... obs_bodyElem_r(obsSpaceData, ... ,bodyIndex)   
    !             end do
    !
    implicit none
    
    ! Arguments:
    type(struct_columnData), intent(inout) :: columnTrl
    type(struct_obs),intent(inout)::obsSpaceData ! Observation-space data object
    integer, intent(in) :: kmode
    type(struct_columnData), intent(inout), optional :: columnAnlInc_opt
    real(8), intent(out), optional :: jobs_opt
    integer, intent(in), optional :: destObsColumn_opt

    ! Local variables
    real(8) :: zomp,zinc,zoer,zhbht
    integer, external :: fclos

    ! Obs space local variables

    integer :: headerIndex,bodyIndex,bodyIndex_start,bodyIndex_end
    integer :: icodtyp,iobslev,nobslev,varno
    integer :: destObsColumn
    character(len=12) :: stnid

    integer, allocatable :: ixtr(:),iass(:),flag(:)
    logical, allocatable :: success(:),process_obs(:)

    ! Model space profile local variables

    real(8), allocatable :: obs_col(:)
    real(8), pointer :: col(:),model_col(:)
    integer :: nlev_bkgrnd,nlev_inc,imodlev
    character(len=4), parameter :: varLevel = 'TH'
    
    type(struct_chm_obsoperators) :: obsoper

    ! Apply setup on first call
    if (.not.initializedChem) then
       call oopc_setupCH
       initializedChem = .true.
    end if
    
    if ((kmode.eq.2.or.kmode.eq.3) .and. (.not.present(columnAnlInc_opt))) &
       call utl_abort("oopc_CHobsoperators: columnAnlInc_opt must be specified for kmode = " // utl_str(kmode))
    
    ! Initializations
    if ( present(destObsColumn_opt) ) then
      destObsColumn = destObsColumn_opt
    else
      destObsColumn = obs_omp
    end if

    if (present(jobs_opt)) jobs_opt = 0.d0

    nlev_bkgrnd = col_getNumLev(columnTrl,varLevel)
    
    ! Allocate memory for model_col. Not necessary for kmode=0 since model_col points to obsoper%trial.
    select case(kmode)
    case(2)
       nlev_inc = col_getNumLev(columnAnlInc_opt,varLevel)
       allocate(model_col(nlev_inc))
    case(1,3)
       allocate(model_col(nlev_bkgrnd))
    end select

    ! Allocations outside chm_obsoper_init since this can be done outside the HEADER loop.
    ! See chm_obsoper_init for assignment of array content.
    
    ! Model obs background, height, TT, and HU profiles.
    allocate(obsoper%trial(nlev_bkgrnd),obsoper%height(nlev_bkgrnd),obsoper%tt(nlev_bkgrnd),obsoper%hu(nlev_bkgrnd))
    ! Model PP and pressure model layer boundaries taken as the middle between model levels.
    allocate(obsoper%pp(nlev_bkgrnd),obsoper%vmodpress(nlev_bkgrnd+1))
    ! Work array: derivative for any variable transform to be applied to a profile
    allocate(obsoper%dtransform(nlev_bkgrnd))
    ! Work array: integration weights array associated to the second order Lagrangian interpolation - see routines
    ! chm_vertintg* 
    allocate(obsoper%vweights(nlev_bkgrnd,nlev_bkgrnd))
    
    ! Loop over all header indices of the 'CH' family:
    
    call obs_set_current_header_list(obsSpaceData,'CH')
    HEADER: do

      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER
  
      icodtyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
      if (icodtyp.ne.codtyp_get_codtyp('CHEMREMOTE').and.icodtyp.ne.codtyp_get_codtyp('CHEMINSITU')) cycle HEADER
      
      stnid = obs_elem_c(obsSpaceData,'STID',headerIndex)
      bodyIndex_start = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
      bodyIndex_end = bodyIndex_start+obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex)-1

      ! Set number of obs profile elements by removing count of BUFR_SCALE_EXPONENT elements
      nobslev = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex)  
      do bodyIndex=bodyIndex_start,bodyIndex_end
         if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex).eq.BUFR_SCALE_EXPONENT) nobslev = nobslev-1
      end do

      ! varno is expected to be the same for all profile points where OBS_VNM value .ne. BUFR_SCALE_EXPONENT
      do bodyIndex=bodyIndex_start,bodyIndex_end
         if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex).ne.BUFR_SCALE_EXPONENT) then
            varno = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
            exit
         end if
      end do

      ! Allocate memory for remaining profile data not in obsoper
      allocate(obs_col(nobslev),success(nobslev),ixtr(nobslev),iass(nobslev),process_obs(nobslev),flag(nobslev))

      ! Check to see if background error variances available
      if (kmode.eq.1) process_obs(:) = bchm_StatsExistForVarName( &
         vnl_varnameFromVarnum(varno,obs_headElem_i(obsSpaceData,OBS_CHM,headerIndex),modelName))
 
      ! Prepare for checking if any processing is needed according to initial flag values     
      iobslev=0
 
      do bodyIndex=bodyIndex_start,bodyIndex_end
         if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex).ne.BUFR_SCALE_EXPONENT) then

            iobslev=iobslev+1

            ixtr(iobslev) = obs_bodyElem_i(obsSpaceData,OBS_XTR,bodyIndex) ! indicates if obs extends outside model profile vertical range
            iass(iobslev) = obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) ! indicates if obs is to be assimilated
            flag(iobslev) = obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex) ! observation integer flag
               
            ! Indicates if this obs should be processed by chm_obsoperators
            if (kmode.eq.1) then
               process_obs(iobslev) = ixtr(iobslev).eq.0.and.iass(iobslev).eq.obs_assimilated.and.process_obs(iobslev)
            else
               process_obs(iobslev) = ixtr(iobslev).eq.0.and.iass(iobslev).eq.obs_assimilated
            end if

         end if
      end do

      ! Initialize processing success flag
      success(1:nobslev) = process_obs(1:nobslev)
      
      if (all(.not.process_obs)) then

         ! All observations in the profile flagged so can skip obs operator for current measurement

         if (kmode.eq.3) then
            model_col(:) = 0.0D0
            obsoper%varName = vnl_varnameFromVarnum(varno,obs_headElem_i(obsSpaceData,OBS_CHM,headerIndex),modelName)
         end if

      else  

         if ( kmode == 1 ) then
           ! Check if "OmP error std dev" is already available
            call obs_set_current_body_list(obsSpaceData,headerIndex)
            iobslev=0
            BODYINDEX1: do bodyIndex=bodyIndex_start,bodyIndex_end
              if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) == BUFR_SCALE_EXPONENT) cycle
              iobslev = iobslev + 1 
              if ( process_obs(iobslev) ) then
                if (obs_bodyElem_r(obsSpaceData,OBS_OMPE,bodyIndex) > 0.0d0 ) then
                  ! "OmP error std dev" is already available for this measurement. Go to the next measurement.
                  
                  ! write(*,*) 'OMPE output ',stnid,obs_bodyElem_r(obsSpaceData,OBS_OMPE,bodyIndex)
                  ! TEMPORARY: First, estimate OBS_HPHT for storage in output files (in the event it is needed externally for
                  ! other purposes (e.g. total column ozone bias correction and corresponding re-doing for marker settings)               
                  if ( obs_bodyElem_r(obsSpaceData,OBS_OMPE,bodyIndex) > 1.1d0*obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex) ) then
                    zhbht = sqrt(obs_bodyElem_r(obsSpaceData,OBS_OMPE,bodyIndex)**2-obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex)**2)
                    call obs_bodySet_r(obsSpaceData,OBS_HPHT,bodyIndex,zhbht)
                  else
                    call obs_bodySet_r(obsSpaceData,OBS_HPHT,bodyIndex,0.5d0*obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex))
                  end if
     
                  deallocate(process_obs,success,ixtr,iass,obs_col,flag)
                  
                  cycle HEADER
                else
                  ! Proceed with the calc of sqrt(diag(HBHT))
                  exit BODYINDEX1
                end if
              end if
            end do BODYINDEX1
         end if

         ! Initialize obsoper variables and allocate arrays
         call chm_obsoper_init(obsoper,obsSpaceData,headerIndex,columnTrl,nlev_bkgrnd,nobslev,kmode,varno,stnid)
 
         ! Initialize model_col, dependent on kmode. Used for input for kmode=0,2, output for kmode=3.
         ! model_col represents for kmode 0) the horizontally interpolated background H_horiz(x_b)
         !                                1) not used
         !                                2) the analysis increment H_horiz dx
         !                                3) the result of applying the adjoint of H_vert 

         select case(kmode)
         case(0)
            model_col => obsoper%trial
         case(2)
            do imodlev=1,nlev_inc
               model_col(imodlev) = col_getElem(columnAnlInc_opt,imodlev,headerIndex,obsoper%varName)
            end do
         case(1,3)
            model_col(:) = 0.0D0
         end select
               
         ! Loop over all body indices (profile elements) to aquire remaining data
      
         iobslev=0
         call obs_set_current_body_list(obsSpaceData,headerIndex)
         BODY1: do

            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if (bodyIndex < 0) exit BODY1
            
            ! Get position in profile and skip over BUFR_SCALE_EXPONENT elements

            if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex).ne.BUFR_SCALE_EXPONENT) then
               iobslev=iobslev+1
            else
               cycle BODY1
            end if

            ! Get vertical coordinate data. Valid for point data values in profiles.
            ! For layer data values, vertical coordinate data will instead be assigned within chm_obsoperators.

            obsoper%obslev(iobslev) = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)

            ! Get normalized increment
            if (kmode.eq.3) then
               if (iass(iobslev).eq.1) then
                  obs_col(iobslev) = obs_bodyElem_r(obsSpaceData,OBS_WORK,bodyIndex)
               else
                  obs_col(iobslev) = 0.0D0
               end if
            end if

         enddo BODY1
      
         ! Apply observation operator. model_col,obs_col are the inputs/outputs in model,observation
         ! space, respectively. Other required inputs are in obsoper. Input/output is as follows:
         !
         !    kmode      model_col      obs_col
         !    -----      ---------      -------
         !      0           in            out
         !      1         not used        out
         !      2           in            out
         !      3           out           in

         call chm_obsoperators(obsoper,model_col,obs_col,kmode,ixtr,success)
 
      end if

      ! Output results
      
      if (kmode.eq.3) then
         ! Store H^T * R^-1 (OmP-Hdx) in columnAnlInc
                     
         col => col_getColumn(columnAnlInc_opt,headerIndex,obsoper%varName)
         col(1:nlev_bkgrnd) = model_col(1:nlev_bkgrnd)

      else
         ! Store results in obsSpaceData

         iobslev=0
         call obs_set_current_body_list(obsSpaceData,headerIndex)
         BODY2: do
            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if (bodyIndex < 0) exit BODY2

            ! Get position in profile and skip over BUFR_SCALE_EXPONENT elements
            if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex).ne.BUFR_SCALE_EXPONENT) then
               iobslev=iobslev+1
            else
               cycle BODY2
            end if

            ! Check for success in calculations
            if (process_obs(iobslev).and..not.success(iobslev)) then
               ! Observation was flagged within this call of oopc_CHobsoperators
               call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
               call obs_bodySet_r(obsSpaceData,OBS_OMP,bodyIndex,0.0D0)
               call obs_bodySet_r(obsSpaceData,OBS_OMA,bodyIndex,0.0D0)
               call obs_bodySet_r(obsSpaceData,OBS_HPHT,bodyIndex,0.0D0)
               call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex,0.0D0)
               call obs_bodySet_i(obsSpaceData,OBS_FLG,bodyIndex, ibset(obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex),9) )
               cycle BODY2
            else if (iass(iobslev).eq.0) then
               ! Observation was flagged previous to this call of oopc_CHobsoperators
               call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
               call obs_bodySet_r(obsSpaceData,OBS_OMP,bodyIndex,0.0D0)
               call obs_bodySet_r(obsSpaceData,OBS_OMA,bodyIndex,0.0D0)
               call obs_bodySet_r(obsSpaceData,OBS_HPHT,bodyIndex,0.0D0)
               call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex,0.0D0)
               cycle BODY2
            else if (.not.process_obs(iobslev) .and. kmode == 1) then
               call obs_bodySet_r(obsSpaceData,OBS_HPHT,bodyIndex,0.0D0)
               cycle BODY2
            end if

            ! Store result in appropriate location in obsSpaceData
            select case(kmode)
            case(0)
            
               ! Store OmP in OBS_OMP and add to Jo(x_background) of CH.   
                           
               zomp = obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex) - obs_col(iobslev)
               call obs_bodySet_r(obsSpaceData,destObsColumn,bodyIndex,zomp)

               if (oopc_diagn_only('CH',stnid,varno,nobslev,flag(iobslev))) then
                  ! Observation is for diagnostics and is not to be assimilated
                  call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
               else if (present(jobs_opt).and.iass(iobslev).eq.1) then
                  ! Add to Jo contribution (factor of 0.5 to be applied outside report loop)
                  zinc = zomp/obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex)
                  jobs_opt = jobs_opt + zinc**2
               end if

            case(1)
            
               ! Background error standard deviations in
               ! observation space, sqrt(diag(H*B_static*H^T)), 
               ! saved in OBS_HPHT of obsSpaceData.
               ! Resulting OmP error std dev estimate saved in OBS_OMPE

               call obs_bodySet_r(obsSpaceData,OBS_HPHT,bodyIndex,obs_col(iobslev))
               
               zoer = obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex)
               zomp = sqrt(obs_col(iobslev)*obs_col(iobslev) + zoer*zoer)
               call obs_bodySet_r(obsSpaceData,OBS_OMPE,bodyIndex,zomp)

            case(2)
            
               !   Store Hdx in OBS_WORK of obsSpaceData               
               call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex,obs_col(iobslev))
            end select

         enddo BODY2

      end if

      ! Deallocate profile data
      deallocate(process_obs,success,ixtr,iass,obs_col,flag)
      call chm_obsoper_dealloc(obsoper)
  
    enddo HEADER

    deallocate(obsoper%trial,obsoper%pp,obsoper%tt,obsoper%dtransform)
    deallocate(obsoper%height,obsoper%vmodpress,obsoper%vweights,obsoper%hu)
    if (kmode.ne.0) deallocate(model_col)
  
    if (present(jobs_opt)) jobs_opt = 0.5d0*jobs_opt
    
  end subroutine oopc_CHobsoperators

  !--------------------------------------------------------------------------
  ! chm_obsoper_init
  !--------------------------------------------------------------------------
  subroutine chm_obsoper_init(obsoper,obsSpaceData,headerIndex,columnTrl, &
                              nmodlev,nobslev,kmode,varno,stnid)
    !
    !:Purpose: To initialize struct_chm_obsoperators variables and to allocate
    !          arrays.
    !
    !:Comments: 
    !           - Allocation of arrays that are dependent on only nlev_bkgrd
    !             (nmodlev) have been moved outside this subroutine so that they
    !             are allocated only once.
    !
    !:Arguments:
    !
    !     :columnTrl:  Column of x_background interpolated to observation
    !                  location. Can have the same vertical levels as the
    !                  trial field (columnTrlOnTrlLev) or as the increment field
    !                  (columnTrlOnAnlIncLev)
    !     :kmode:       
    !   
    !                - 0 for non-linear/linear model in assimilation (all models
    !                  included are currently linear)
    !                - 1 for determination of sqrt(diag(H*B*H^T))
    !                - 2 for tangent linear model
    !                - 3 for adjoint model 
    implicit none

    ! Arguments:
    type(struct_chm_obsoperators), intent(inout) :: obsoper ! Object for constituents associated to obs
    type(struct_obs), intent(inout) :: obsSpaceData ! Obs-Space Data object
    integer, intent(in) :: headerIndex ! Measurement index in obsSpaceData
    type(struct_columnData), intent(inout) :: columnTrl
    integer, intent(in) :: nmodlev ! Number of background field (model) levels
    integer, intent(in) :: nobslev ! Number of obs elements (see chm_obsoper_proceed)
    integer, intent(in) :: kmode
    integer, intent(in) :: varno ! BUFR number
    character(len=12), intent(in) :: stnid ! Station ID

    ! Locals:
    integer :: bodyIndex ! Measurement element index in obsSpaceDate (see chm_obsoper_proceed)
    integer :: jl,nmodlev_uv
    real(8), pointer :: col_height_ptr(:)
    real(8), allocatable :: uu(:),vv(:)
    character(len=4), parameter :: varLevel = 'TH'
    real(8) :: checkID

    obsoper%nmodlev = nmodlev
    obsoper%nobslev = nobslev
    obsoper%obs_index = headerIndex
    obsoper%varno = varno
    obsoper%stnid = stnid

    ! Get obs space info that are part of the profile header
    obsoper%date  = obs_headElem_i(obsSpaceData,OBS_DAT,headerIndex)
    obsoper%hhmm  = obs_headElem_i(obsSpaceData,OBS_ETM,headerIndex)
    ! Constituent identifyer following local version of WMO GRIB Table 08046 (similar to BUFR Table 08043)
    obsoper%constituent_id = obs_headElem_i(obsSpaceData,OBS_CHM,headerIndex)
    
    ! Check if constituent id is recognized (function will abort if not recognized)
    if ( obsoper%constituent_id >= 0 .and. obsoper%constituent_id < 7000) checkID = vnl_varMassFromVarnum(obsoper%constituent_id)

    obsoper%lat = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex) 
    obsoper%lon = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex) 

    ! Body info that we only need for first point in the profile
    bodyIndex       = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex) 
    ! Index of vertical coordinate type    
    obsoper%vco     = obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex)  
    ! Model field name (NOMVAR value)
    obsoper%varName = vnl_varnameFromVarnum(obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex),obsoper%constituent_id,modelName)

    ! Allocate arrays

    allocate(obsoper%obslev(nobslev))       ! Reference vertical levels
    allocate(obsoper%vlayertop(nobslev))    ! Layer tops for layer measurements
    allocate(obsoper%vlayerbottom(nobslev)) ! Layer bottoms for layer measurements
    allocate(obsoper%imodlev_top(nobslev))  ! Index of highest model level (lowest index) involved with obs element
    allocate(obsoper%imodlev_bot(nobslev))  ! Index of lowest model level (highest index) involved with obs element

    allocate(obsoper%zh(nobslev,nmodlev))   ! Local model operator H (excluding conversion constants and horizontal interpolation)
    allocate(obsoper%zhp(nobslev,nmodlev))  ! Part of zh that excludes aspects related to vertical resolition

    obsoper%vweights(:,:)=0.0D0            
    obsoper%zh(:,:)=0.0D0
    obsoper%zhp(:,:)=0.0D0
    obsoper%imodlev_top(:)=1
    obsoper%imodlev_bot(:)=nmodlev
    obsoper%dtransform(:)=1.0D0

    if (.not.col_varExist(columnTrl,'TT')) then
       if (chm_required_field('TT',obsoper%varno)) then
          call utl_abort("chem_opsoper_init: TT required for BUFR code " // trim(utl_str(obsoper%varno)))
       end if
    end if

    ! Get background profiles at observation location
    do jl=1,nmodlev
       obsoper%pp(jl) = col_getPressure(columnTrl,jl,headerIndex,varLevel)
       obsoper%trial(jl) = col_getElem(columnTrl,jl,headerIndex,obsoper%varName)
       obsoper%tt(jl) = col_getElem(columnTrl,jl,headerIndex,'TT')
    enddo

    if (col_varExist(columnTrl,'TT').and.col_varExist(columnTrl,'HU').and.col_varExist(columnTrl,'P0')) then     
       ! Height would have been generated in the call to sugomobs. 
       ! Convert from geopotential to geopotential height.
       col_height_ptr => col_getColumn(columnTrl,headerIndex,'Z_T')
       obsoper%height(1:nmodlev) = col_height_ptr(1:nmodlev)
    else
       obsoper%height(:) = -1.
    end if

    ! Get specific humidity if available
    if (col_varExist(columnTrl,'HU')) then
       do jl=1,nmodlev
         obsoper%hu(jl) = col_getElem(columnTrl,jl,headerIndex,'HU')       ! lnq was replaced by q
       enddo
    else
       obsoper%hu(:)=-1
    end if
    
    ! If applicable, get column upper boundaries for use with total column measurements when the related increment profile
    ! is to be restricted to the lower atmosphere (e.g. troposphere or PBL; 
    ! when tropo_bound>0 )
    if (obsoper%vco.eq.4.and.nobslev.eq.1.and.kmode.ne.1) then
       if (kmode.eq.0) then
          
          if (col_varExist(columnTrl,'HU')) then
             if (col_varExist(columnTrl,'UU').and.col_varExist(columnTrl,'VV')) then
                nmodlev_uv=col_getNumLev(columnTrl,'MM')
                allocate(uu(nmodlev_uv),vv(nmodlev_uv))
                do jl=1,nmodlev_uv
                   uu(jl) = col_getElem(columnTrl,jl,headerIndex,'UU')
                   vv(jl) = col_getElem(columnTrl,jl,headerIndex,'VV')
                enddo
                obsoper%column_bound = chm_get_col_boundary(obsoper%constituent_id,nmodlev,obsoper%pp,obsoper%tt,obsoper%height,hu_opt=obsoper%hu,uu_opt=uu,vv_opt=vv)
                deallocate(uu,vv)
             else 
                obsoper%column_bound = chm_get_col_boundary(obsoper%constituent_id,nmodlev,obsoper%pp,obsoper%tt,obsoper%height,hu_opt=obsoper%hu)   
             end if

          else
              obsoper%column_bound = chm_get_col_boundary(obsoper%constituent_id,nmodlev,obsoper%pp,obsoper%tt,obsoper%height)
          end if

          call chm_add_col_boundary(headerIndex,obsoper%column_bound)  ! save boundary for kmode>0 calls using headerIndex

       else
          obsoper%column_bound = chm_retrieve_col_boundary(headerIndex)
       end if
    else
       obsoper%column_bound = -1.
    end if

  end subroutine chm_obsoper_init

  !--------------------------------------------------------------------------
  ! chm_obsoper_dealloc
  !--------------------------------------------------------------------------
  subroutine chm_obsoper_dealloc(obsoper)
    !
    !:Purpose: To deallocate arrays for struct_chm_obsoperators.
    !
    !:Comments:
    !           - Deallocation of arrays that are dependent on only nmodlev have
    !             been moved outside this subroutine so that they are
    !             deallocated only once.
    !
    implicit none

    ! Arguments:
    type(struct_chm_obsoperators), intent(inout) :: obsoper

    if (allocated(obsoper%obslev))       deallocate(obsoper%obslev)
    if (allocated(obsoper%vlayertop))    deallocate(obsoper%vlayertop)
    if (allocated(obsoper%vlayerbottom)) deallocate(obsoper%vlayerbottom)
    if (allocated(obsoper%zh))           deallocate(obsoper%zh)
    if (allocated(obsoper%zhp))          deallocate(obsoper%zhp)
    if (allocated(obsoper%imodlev_top))  deallocate(obsoper%imodlev_top)
    if (allocated(obsoper%imodlev_bot))  deallocate(obsoper%imodlev_bot)

  end subroutine chm_obsoper_dealloc
    
  !--------------------------- Routines for observation operators ----------------
  !-------------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! chm_obsoperators
  !--------------------------------------------------------------------------
  subroutine chm_obsoperators(obsoper,model_col,obs_col,kmode,ixtr,success)
    !
    !:Purpose: To apply observation operator for indicated observation data and
    !          condition.
    !
    !:Arguments:
    !     :kmode: mode of observation observational operator
    !               0) general (potentially non-linear) simulation operator
    !               1) determination of sqrt(diag(H*B*H^T))
    !               2) tangent linear operator
    !               3) linear adjoint operator
    !
    !     :obs_col:   Observation space input/output profile
    !
    !                  +------+------+------+------+------+------+------+
    !                  |kmode |    input/output    | profile            |
    !                  +======+====================+====================+
    !                  |  0   |        out         | H(xb)              |
    !                  +------+--------------------+--------------------+
    !                  |  1   |        out         | sqrt(diag(H*B*H^T))|
    !                  +------+--------------------+--------------------+
    !                  |  2   |        out         | H*dx               |
    !                  +------+--------------------+--------------------+
    !                  |  3   |        in          | R**-1 (Hdx-d)      |
    !                  +------+--------------------+--------------------+
    !
    !     :model_col: Model space input/output profile
    !
    !                  +------+------+------+------+------+------+------+
    !                  |kmode |   input/output     | profile            |
    !                  +======+====================+====================+
    !                  |  0   |       in           | xb                 |
    !                  +------+--------------------+--------------------+
    !                  |  1   |       not used     | not used           |
    !                  +------+--------------------+--------------------+
    !                  |  2   |       in           | dx at obs location |
    !                  +------+--------------------+--------------------+
    !                  |  3   |       out          | adjoint product    |
    !                  |      |                    | H^T(...)           |
    !                  +------+--------------------+--------------------+
    !
    !     :ixtr:      Flag indicating whether obs within the model vertical
    !                 coord range (.ne.0 for no). Can be modified internally
    !                 - hence intent(inout) - even though these changes will not
    !                 be needed outside this routine.
    !
    ! Further changes required for generalization:
    !
    !  1) Add layer average operators.
    !  2) Add AOD operators (summation over model layers).
    !  3) Add option to include use of obs error correlation matrix for kmode=2,3 
    !     (This may/will need to be done in oop_Hchm and oop_HTchm where the division 
    !      by sigma_obs is applied. A new routine will be needed for this
    !      operation - and others for reading the correlation matrices similarly to the
    !      averaging kernels.)
    !
    ! Comments:
    !      - When kmode=0, call from chm_observation_operators passes model_col as a pointer to
    !        obsoper%trial.
    !      - Does not yet account for potential future applications of obs 
    !        vertical correlation matrices.
    !      - Potential specification of background error std. dev. (sigma_trial(:,2)) and correlation matrices 
    !        for the ensemble-based and lam cases to be done when stats for these become in use with constituents.
    !
  implicit none

  ! Arguments
  type(struct_chm_obsoperators), intent(inout) :: obsoper ! information needed by the observation operator

  ! I/O arguments: obs space variables
    
  integer, intent(in) :: kmode 
  integer, intent(inout) :: ixtr(obsoper%nobslev)
  logical, intent(inout) :: success(obsoper%nobslev) ! Indicates whether the observation was successfully assimilated

  ! I/O arguments: model space profile data and others

  real(8), intent(inout) :: model_col(obsoper%nmodlev), obs_col(obsoper%nobslev)

  ! Local variables
  
  real(8) :: press_obs(obsoper%nobslev),zwork(obsoper%nmodlev),unit_conversion(obsoper%nmodlev)
  real(8) :: sigma_trial(obsoper%nmodlev,2),temp_eff(1),zhmin
  integer :: iobslev,imodlev,jvar
  integer, parameter :: code_len=90
  character(len=code_len) :: code    ! Must be at least as large as chm_code_len
  real(8), allocatable :: avg_kern(:,:)

  if (code_len.lt.oss_obsdata_code_len()) call utl_abort('chm_obsoperators: Length of code string' &
                                          // ' needs to be increased to ' // trim(utl_str(oss_obsdata_code_len())))
 
  ! Determine if layer boundaries are assigned to this data source.
  ! If so, obtain them for use in this routine. 
  ! Routine provides obsoper%layer_identified,
  !                  obsoper%vlayertop(nobslev), 
  !                  obsoper%vlayerbottom(nobslev) 
 
  call chm_get_layer_boundaries(obsoper%stnid,obsoper%varno,obsoper%vco, &
       obsoper%nobslev,obsoper%pp(1),obsoper%pp(obsoper%nmodlev), &
       obsoper%layer_identified,obsoper%vlayertop,obsoper%vlayerbottom)
         
  ! Identify observation operator based on observation units and presence or
  ! not of layer boundaries

  if (bufr_IsIntegral(obsoper%varno)) then
     if (.not.obsoper%layer_identified) then
        write(*,*)   '----------------------------------------------------------'
        write(*,*)   'STNID, BUFR index, nobslev: ',obsoper%stnid,' ',obsoper%varno,obsoper%nobslev
        call utl_abort('chm_obsoperators: Required layer boundaries not available!')
     else
        ! Vertical integration operator
        obsoper%modelIndex=3
     end if
  else if (obsoper%layer_identified) then
      ! Layer averaging operator
     obsoper%modelIndex=4
  else if (obsoper%vco.eq.5.and.obsoper%nobslev.eq.1) then  
      ! Surface point (in-situ) measurement
      obsoper%modelIndex=2
  else    
      ! Vertical interpolation operator
     obsoper%modelIndex=1
  end if  

  ! Indicate if the generalized innovation operator is to be applied.

  obsoper%apply_genoper=.false.
  if (obsoper%constituent_id.ge.0) then
    if (chm_generalized_operator(obsoper%constituent_id) > 0 .and. kmode /= 1 .and. &
        (obsoper%modelIndex == 3 .or. obsoper%modelIndex == 4)) then
       
      if (kmode.eq.0) then
          ! Set reference profile for use with generalized innovation operator when kmode>=2
          call chm_set_reference_obsdata(obsoper)
  
          ! Get background error std dev profile at obs locations
          sigma_trial(:,:)=0.D0
          call bchm_getsigma(obsoper%varName,obsoper%nmodlev,obsoper%lat,obsoper%lon,sigma_trial(:,1),vlev_opt=obsoper%pp) 
!!          call blamchm_getsigma(obsoper%varName,obsoper%nmodlev,obsoper%lat,obsoper%lon,sigma_trial(:,1)) 
!!          call benschm_getsigma(obsoper%varName,obsoper%nmodlev,obsoper%lat,obsoper%lon,sigma_trial(:,2)) 

          call chm_add_sigma_trial(obsoper%obs_index,sigma_trial)

      else if (kmode.ge.2) then
          obsoper%apply_genoper = .true.
      end if
    end if
  end if

  ! Apply unit conversion (apply unit conversion later for kmode=3)

  select case(kmode)
  case(0)
     ! Perform transformation and unit conversion on obsoper%trial.
     ! Note that model_col => obsoper%trial for kmode=0.
     call chm_convert_units(obsoper,obsoper%trial)
  case(1)
     ! Save the conversion factor in <unit_conversion>
     unit_conversion(:) = 1.0
     call chm_convert_units(obsoper,unit_conversion,incr_opt=.true.)
  case(2)
     call chm_convert_units(obsoper,model_col,incr_opt=.true.)
  end select

  if (obsoper%apply_genoper) then
     ! Perform unit conversion on obsoper%trial when applying the generalized
     ! obs operator for kmode=2,3. Keep obsoper%trial in ug/kg in this case.
     call chm_convert_units(obsoper,obsoper%trial,ppb_opt=.true.)
  end if

  ! Convert observation vertical coordinate value(s) to pressure if needed

  select case(obsoper%vco)
  case(1)
     ! Convert altitude to pressure
     select case(obsoper%modelIndex)
     case(1)
        press_obs = phf_convert_z_to_pressure(obsoper%obslev,obsoper%height,obsoper%pp,obsoper%nobslev,obsoper%nmodlev,obsoper%lat,success)
     case(2,3)
        obsoper%vlayertop = phf_convert_z_to_pressure(obsoper%vlayertop,obsoper%height,obsoper%pp,obsoper%nobslev,obsoper%nmodlev,obsoper%lat,success)
        obsoper%vlayerbottom = phf_convert_z_to_pressure(obsoper%vlayerbottom,obsoper%height,obsoper%pp,obsoper%nobslev,obsoper%nmodlev,obsoper%lat,success)
     end select
  case(2)
     ! Pressure, no conversion needed
     if (obsoper%modelIndex.eq.1) press_obs = obsoper%obslev
  case(4,5)
     ! No actions taken
  case default
     call utl_abort("chm_obsoperators: vertical coordinate type vco = " // trim(utl_str(obsoper%vco)) // " not available for this operator.")
  end select

  ! Determine if averaging kernel is to be applied

  obsoper%iavgkern = chm_find_avgkern(obsoper%stnid,obsoper%varno,obsoper%nobslev)

  ! Apply appropriate core observation operator
   
  select case(obsoper%modelIndex)
  case(1)

     ! Vertical interpolation operator

     call chm_vert_interp_operator(obsoper,press_obs,ixtr,success)

  case(2)

     ! Surface point measurement

     ! Set weight to unity for lowest model level.
     obsoper%zh(1,1:obsoper%nmodlev-1)=0.0
     obsoper%zh(1,obsoper%nmodlev) = 1.0

     ! Set range of elements for model vertical levels
     obsoper%imodlev_top(1) = obsoper%nmodlev
     obsoper%imodlev_bot(1) = obsoper%nmodlev 
    
  case(3)
  
     ! Layer integration operator

     call chm_layer_integ_operator(obsoper,ixtr,success,kmode)

!  case(4)

!    Layer averaging operator

!    call chm_layer_avg_operator  ! see 3dvar_chem routine ch_vavg

  end select
  
  ! Apply averaging kernels if requested

  if (obsoper%iavgkern.gt.0) then

     allocate(avg_kern(obsoper%nobslev,obsoper%nobslev))
     
     call chm_get_avgkern(obsoper%iavgkern,obsoper%stnid,obsoper%nobslev,obsoper%lat,obsoper%lon,obsoper%date,obsoper%hhmm,avg_kern)
     do iobslev=1,obsoper%nobslev
        if (success(iobslev)) then

           ! Apply averaging kernels to observation operator(s)
           obsoper%zh(iobslev,:)= matmul(avg_kern(iobslev,:),obsoper%zh(:,:))
           if (obsoper%apply_genoper) obsoper%zhp(iobslev,:) = matmul(avg_kern(iobslev,:),obsoper%zhp(:,:))

           ! Extend vertical range of obs operator according to the influence of
           ! the averaging kernel. Either extend to the entire model vertical range
           ! (commented out below) or to the vertical range with non-negligable values.

!          obsoper%imodlev_top(iobslev) = 1
!          obsoper%imodlev_bot(iobslev) = obsoper%nmodlev

           zhmin=1.0D-10*maxval(abs(obsoper%zh(iobslev,:)))
           do imodlev=1,obsoper%imodlev_top(iobslev)
              if (abs(obsoper%zh(iobslev,imodlev)).gt.zhmin) exit
           end do
           if (imodlev.gt.obsoper%imodlev_top(iobslev)) imodlev=obsoper%imodlev_top(iobslev)
           obsoper%imodlev_top(iobslev) = imodlev
           do imodlev=obsoper%nmodlev,obsoper%imodlev_bot(iobslev),-1
              if (abs(obsoper%zh(iobslev,imodlev)).gt.zhmin) exit
           end do
           if (imodlev.lt.obsoper%imodlev_bot(iobslev)) imodlev=obsoper%imodlev_bot(iobslev)
           obsoper%imodlev_bot(iobslev) = imodlev

        end if
     end do

     deallocate(avg_kern)

  end if

  ! Apply generalized innovation operator if requested

  if (obsoper%apply_genoper) call chm_genoper(obsoper,kmode,success)

  ! Finalize required quantities depending on kmode
   
  select case(kmode)

  case(0,2)

      ! Finalize non-linear/linear operator step

      do iobslev=1,obsoper%nobslev
        if (success(iobslev)) then
           obs_col(iobslev)=dot_product(obsoper%zh(iobslev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)), &
                                 model_col(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)))
        end if
      end do

      ! Calculate concentration-weighted effective temperature (for output purpose)

      if (kmode.eq.0) then 

         if ((obsoper%constituent_id.ge.0.and.obsoper%constituent_id.lt.7000).and. &
            obsoper%modelIndex.eq.3.and.obsoper%nobslev.eq.1.and.obsoper%vco.eq.4.and.success(1)) then
            if (all(obsoper%tt.gt.0.0).and.obs_col(1).gt.0.0) then
                temp_eff(1)=dot_product(obsoper%zh(1,obsoper%imodlev_top(1):obsoper%imodlev_bot(1)), &
                        obsoper%tt(obsoper%imodlev_top(1):obsoper%imodlev_bot(1))*model_col(obsoper%imodlev_top(1):obsoper%imodlev_bot(1))) &
                        /obs_col(1)
                code=oss_obsdata_get_header_code(obsoper%lon,obsoper%lat,obsoper%date,obsoper%hhmm,obsoper%stnid)
                call oopc_add_efftemp_obsdata(code,temp_eff)
             end if
          end if
       end if
                       
  case(1)

     ! Compute sqrt(diag(H*B*H^T))

     ! Apply unit conversion to observation operator
     do iobslev=1,obsoper%nobslev
       if (success(iobslev)) then
          obsoper%zh(iobslev,:) = unit_conversion * obsoper%zh(iobslev,:)
       end if
     end do

     ! Get background error std dev profile at obs location
     sigma_trial(:,:)=0
     call bchm_getsigma(obsoper%varName,obsoper%nmodlev,obsoper%lat,obsoper%lon,sigma_trial(:,1)) 
!!      call blamchm_getsigma(obsoper%varName,obsoper%nmodlev,obsoper%lat,obsoper%lon,sigma_trial(:,1)) 
!!      call benschm_getsigma(obsoper%varName,obsoper%nmodlev,obsoper%lat,obsoper%lon,sigma_trial(:,2)) 

     ! Identify variable position index in background error correlation matrices

     jvar=1
     do while (trim(bchm_varnamelist(jvar)).ne.'') 
        if (trim(bchm_varnamelist(jvar)).eq.trim(obsoper%varName)) exit
        jvar=jvar+1
     end do
     if (trim(bchm_varnamelist(jvar)).eq.'') call utl_abort('chm_genoper: Correlation matrix not found for ' // trim(obsoper%varName) )

     do iobslev=1,obsoper%nobslev
        if (success(iobslev)) then
!!           call chm_corvert_mult(obsoper%varName,obsoper%zh(iobslev,:),obs_col(iobslev), &
!!                obsoper%imodlev_top(iobslev),obsoper%imodlev_bot(iobslev),1,obsoper%nmodlev, &
!!                1,.true.,3,sigma_trial) ! get h*B*h^T

           do imodlev=obsoper%imodlev_top(iobslev),obsoper%imodlev_bot(iobslev)
               zwork(imodlev)=sum(obsoper%zh(iobslev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)) &
                         *bchm_corvert(imodlev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev),jvar) &
                         *sigma_trial(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev),1))*sigma_trial(imodlev,1)
           end do
           obs_col(iobslev)=sum(obsoper%zh(iobslev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)) &
                *zwork(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)))

           obs_col(iobslev) = sqrt(obs_col(iobslev))  ! save as sqrt(h*B*h^T)
        else
           obs_col(iobslev) = 0.0
        end if
     end do

  case(3)

      ! H^T*grad contribution from adjoint of tangent linear model.

      model_col(:) = 0.0

      do iobslev=1,obsoper%nobslev
        if (success(iobslev)) then
           zwork(:)=0.0D0
            zwork(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)) = &
               obs_col(iobslev)*obsoper%zh(iobslev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev))
                
            call chm_convert_units(obsoper,zwork,incr_opt=.true.)
           
            model_col(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)) = &
               model_col(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)) + &
               zwork(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev))
        end if
      end do

  end select

  end subroutine chm_obsoperators

  !--------------------------------------------------------------------------
  ! chm_convert_units
  !--------------------------------------------------------------------------
  subroutine chm_convert_units(obsoper,model_col,ppb_opt,incr_opt)
    !
    !:Purpose: To set unit-conversion factor for consistency of Hx units with
    !          obs units.
    !
    !:Arguments:
    !     :ppb_opt:        indicates whether model_col should be kept in
    !                      ug/kg instead of the units dictated by the BUFR
    !                      number (optional, .false. by default)
    !     :incr_opt:       indicates if model_col is actually an increment 
    !                      (optional, .false. by default). Needed for non-linear
    !                      transformations (i.e. for 'HU')
    !     :model_col       Array to be converted. Either trial or increment-related.
    !
    !:Comments:
    !    
    !      A. Standard model/analysis species field provided as mass mixing 
    !         ratio in ug/kg (ppb). Conversion to ppb is applied when this is 
    !         not the case except for AOD and surface emissions.
    !         As this is hard-coded, any changes in analysis variable must
    !         be reflected by correspondingly modifying this module.
    !
    !      B. Unit-conversion factor is calculated in chm_convert_units
    !         from the following factors:
    !           (1) physical constants
    !           (2) parameters related to a particular species such as molecular
    !               mass
    !           (3) variables such as T and P from background field at each
    !               iteration
    !
    !      C. The baseline integral observation operator can be interpreted as
    !         being integrals of the gas partial pressure, giving products in
    !         kg/m^2, e.g. with sample discretized layer integrals::
    !                (mass density) * dz = - d(gas partial pressure)/g 
    !                                    = - [rho(gas)/rho(air)]*dP/g
    !                                    = - 1E-9 * [mass mixing ratio in parts per billion (ppb)]*dP/g 
    ! 
    !         The actual integration in pressure (in Pascal) is performed
    !         outside this routine. For integral products in kg/m^2, the output
    !         of this routine is to be in mmr/(m/s^2)  (mmr=mass mixing ratio),
    !         which is equivalent to (1E-9 ug/kg)/(m/s^2) and kg/(m^2*Pa). 
    !         Therefore, the input value in ug/kg has to be multiplied by 1E-9/g
    !         (g=RG below).
    ! 
    !         For integral products in other units, additional conversion 
    !         factors are also to be applied.
    !
    !      D. List should be revised following changes to the 'tableburp' file.
    !
    !
    !      E. Coefficients related to unit conversion
    !
    !         rho_stp=1.293                     Air density at STP (1.293 kg/m^3)
    !         RG=9.807 (=g)                     Acceleration due to gravity (m/s^2)
    !         MPC_AVOGADRO_R8 = Na              Avogadro's number. 6.023E23 molecules/mole
    !         MPC_MOLAR_MASS_DRY_AIR_R8 (m_air) Dry air molecular mass. 28.9644 g/mole
    !         MPC_RGAS_IDEAL_R8 = R             Ideal gas constant. 8.341 J/mole/K  (J=kg m^2/s^2)
    !
    !                                             PV = nRT (n=number of moles)
    !
    !         MPC_RGAS_DRY_AIR_R8  = Rd         Dry air constant. 287.1 J/kg/K  (J=kg m^2/s^2)
    !                                           = MPC_RGAS_IDEAL_R8 * 1000 g/km / MPC_MOLAR_MASS_DRY_AIR_R8
    !
    !                                             P=rho*Rd*T = [n*m_air*0.001 kg/g]*Rd*T
    !                                                        = n*[m_air*0.001*Rd]*T
    !                                                        = nRT
    ! Further changes required for generalization 
    !
    !  1) Conversion for surface emissions not included as yet (if any is
    !     needed)
    !
    implicit none

    ! Arguments:
    type(struct_chm_obsoperators), intent(inout) :: obsoper ! basic information related to the observation operator
    real(8), intent(inout) :: model_col(obsoper%nmodlev) ! Model-space profile to have its units changed
    logical, intent(in), optional :: ppb_opt, incr_opt
    
    ! Locals:
    real(8) :: zcoef
    integer :: exp_P,exp_T
    logical :: ppb_out, incr_out
    real(8), parameter :: rho_stp=1.293  ! kg/m^3
    
    if (obsoper%constituent_id.lt.0) return
    
    ! No conversion necessary for these BUFR numbers
    if (any( obsoper%varno .eq. (/ BUFR_UNIT_OptDepth,BUFR_UNIT_OptDepth2, &
            BUFR_UNIT_OptDepth3, BUFR_UNIT_MR_NVaerosol, BUFR_NETT /)  )) return
    
    if (present(ppb_opt)) then
       ppb_out = ppb_opt
    else
       ppb_out = .false.
    end if
      
    if (present(incr_opt)) then
       incr_out = incr_opt
    else
       incr_out = .false.
    end if
    
    zcoef = 1.
    exp_T = 0   ! exponent of multiplicative factor obsoper%tt
    exp_P = 0   ! exponent of multiplicative factor obsoper%pp
    
    ! Convert to ug/kg if not already in those units

    if (obsoper%varName(1:2).eq.'AF' .or. obsoper%varName(1:2).eq.'AC') then

      ! PM2.5 or PM10
       
      if (any(obsoper%varno .eq. (/ BUFR_UNIT_VMR, BUFR_UNIT_VMR2,  &
                    BUFR_UNIT_MolePerMole, BUFR_UNIT_MolePerMole2,   &
                    BUFR_UNIT_NumberDensity, BUFR_UNIT_MolarDensity, & 
                    BUFR_UNIT_PartPress, BUFR_UNIT_PartPress2, &
                    BUFR_UNIT_DU, BUFR_UNIT_DU2, BUFR_UNIT_DU3, BUFR_UNIT_DU4, &
                    BUFR_UNIT_IntegND, BUFR_UNIT_IntegND2 /) )) &
          call utl_abort("chm_convert_units: BUFR # " // trim(utl_str(obsoper%varno)) // " is not valid for PM" )
       
      ! Conversion from ug/m^3 to ug/kg  (scaling by Rd*T/P)
      zcoef = zcoef * MPC_RGAS_DRY_AIR_R8
      exp_T = exp_T+1 ! multiply by T
      exp_P = exp_P-1 ! divide by P
   
    elseif (obsoper%varName(1:2).eq.'HU') then
       
      if (.not.incr_out) then
        ! Converts specific humidity (q) to mass mixing ratio mmr = q/(1-q)
        model_col =  model_col / (1.0d0 - model_col)
      else
        ! For conversion of q increment (dq) to mass mixing ratio increment (dmmr)
        ! dmmr = dq/(1-q)^2         
        model_col = model_col/(1.0d0 - obsoper%trial)**2
      end if  
      ! Conversion factor for kg/kg to ug/kg
      zcoef = zcoef * 1.0d9      
        
    end if
  
    ! Convert from ug/kg to desired unit if ppb_out = .false.

    if (.not.ppb_out) then
       select case (obsoper%varno)
       
       ! The first four cases below are for integral observations which
       ! require a conversion, in this routine, from ug/kg to the units of 
       ! the integrand values for integrals in pressure (Pascal). Comment C above.
       ! Note: 1 ug/kg = 1 ppb = 1E9 mmr (mass mixing ratio)
       
       case(BUFR_UNIT_IntegDens, BUFR_UNIT_IntegDens2, BUFR_UNIT_IntegDens3) 
        
          ! For conversion from ug/kg to integrand values in kg/(m^2*Pa) 
          ! Note: 1 kg/(m^2**Pa) = = 1 mmr / RG = 1E-9 ug/kg / RG 
          !
          zcoef = zcoef * 1.0d-9 / ec_rg
          
       case(BUFR_UNIT_IntegMolarDens) 
        
          ! For conversion from ug/kg to integrand values in moles/(m^2*Pa)
          ! Note: 1 moles/(m^2*Pa) = 1E-9 ug/kg / [RG * (1E-3 kg/g * m_gas)]
          !
          ! To convert from kg/m^2 for the gas to moles/m^2, one must 
          ! divide by the molar mass of the gas in kg/mole.
          !
          ! Note: One u or Da (unified atomic mass unit or dalton) is numerically equivalent to 1 g/mole.
          ! So 1 kg is equivalent to (1E3/(atomic mass)) moles
          
          zcoef = zcoef * 1.0d-6 / (ec_rg*vnl_varMassFromVarNum(obsoper%constituent_id))
          
       case(BUFR_UNIT_IntegND, BUFR_UNIT_IntegND2)
          
          ! For conversion from ug/kg to integrand values in molecules/(m^2*Pa)
          ! Note: 1 molecule/(m^2*Pa) = 1E-9 ug/kg * Na / [RG * (1E-3 kg/g * m_gas)]
          ! 
          ! To convert from kg/m^2 for the gas to molecules/m^2, one must 
          ! divide by the gas molar mass (kg/mole) and multiply by the Avogrado number          

          zcoef = zcoef * 1.0d-6 * MPC_AVOGADRO_R8 &
               / (ec_rg*vnl_varMassFromVarNum(obsoper%constituent_id))

       case(BUFR_UNIT_DU, BUFR_UNIT_DU2, BUFR_UNIT_DU3, BUFR_UNIT_DU4) 
      
          ! For conversion from ug/kg to integrand values in DU/Pa
          ! Note: 1 DU/Pa =  1E-9 ug/kg / [RG * m_gas * rho_stp/m_air * 1E-5 m] 
          !
          ! 1 DU = 0.01 mm of gas at STP = 1E-5 m of gas at STP
          !      = integral of gas number density at STP over 1E-5 m.
          !      = Na*P/(RT) * 1E-5   at STP  (Na=Avogadro's number)
          !      = Na*(molar density) * 1E-5   at STP
          !      = Na*rho(STP)/(molar mass) * 1E-5 m
          !      = integral of air number density at STP over 1E-5 m.
          !      = Na*rho(air,STP)/m_air * 1.E-5 m
          !               
          ! Hence 1 DU equivalent to Na*rho_stp/m_air * 1E-5 m (= 2.69E20 molecules/m^2)
          !
          ! To convert from kg/m^2 for the gas to molecules/m^2, one must 
          ! divide by the gas molar mass (kg/mole) and multiply by the Avogrado number 
          !
          ! To convert from molecules/m^2 to DU, one must divide by 2.69E20 or (Na*rho_stp/m_air * 1E-5).
          ! So for conversion from kg/m^2 to DU, must divide by (m_gas*rho_stp/m_air * 1E-5)       
          
          zcoef = zcoef * 1.0d-4 * MPC_MOLAR_MASS_DRY_AIR_R8 &
                /(vnl_varMassFromVarNum(obsoper%constituent_id)*ec_rg*rho_stp)
        
       case(BUFR_UNIT_Density, BUFR_UNIT_Density2, &
            BUFR_UNIT_AirDensity, BUFR_UNIT_PMDensity)

          ! For conversion from ug/kg to kg/m^3
          !
          ! rho(gas) = mass mixing ratio * rho(air) = mass mixing ratio * P/Rd/T
          
          zcoef = zcoef * 1.0d-9 / MPC_RGAS_DRY_AIR_R8
          exp_T = exp_T-1 ! divide by T
          exp_P = exp_P+1 ! multiply by P
          
       case(BUFR_UNIT_MMR, BUFR_UNIT_MMR2) 

          ! For conversion from ug/kg to kg/kg
     
          zcoef = zcoef * 1.0d-9 
     
       case(BUFR_UNIT_PartPress, BUFR_UNIT_PartPress2) 
     
          ! For conversion from ug/kg to partial pressure (PA)
          !
          ! parial pressure = P * vmr
          !                 = P * m_air/m_gas * mass mixing ratio
          
          zcoef = zcoef * 1.0d-9 * MPC_MOLAR_MASS_DRY_AIR_R8 &
               /vnl_varMassFromVarNum(obsoper%constituent_id)
          exp_P = exp_P+1 ! multiply by P

       case(BUFR_UNIT_NumberDensity)
          
          ! For conversion from ug/kg to molecules/m^3
          !
          ! Number density of gas = Na*rho(gas)/m_gas = Na*rho(air) * mass mixing ratio /m_gas
          !                       = Na * P/Rd/T * mass mixing ratio /m_gas
        
          zcoef = zcoef * 1.0d-6 * MPC_AVOGADRO_R8/MPC_RGAS_DRY_AIR_R8 &
               /vnl_varMassFromVarNum(obsoper%constituent_id)
          exp_T = exp_T-1 ! divide by T
          exp_P = exp_P+1 ! multiply by P

       case(BUFR_UNIT_MolarDensity)
          
          ! For conversion from ug/kg to moles/m^3
          !
          ! Mole density of gas = rho(gas)/m_gas = rho(air) * mass mixing ratio /m_gas
          !                       = P/Rd/T * mass mixing ratio /m_gas

          zcoef = zcoef * 1.0d-6 /MPC_RGAS_DRY_AIR_R8 &
               /vnl_varMassFromVarNum(obsoper%constituent_id)
          exp_T = exp_T-1 ! divide by T
          exp_P = exp_P+1 ! multiply by P

       case(BUFR_UNIT_VMR, BUFR_UNIT_VMR2, BUFR_UNIT_MolePerMole, BUFR_UNIT_MolePerMole2)
          
          ! For conversion from ug/kg to vmr (or moles/mole)
          
          zcoef = zcoef * 1.0d-9 * MPC_MOLAR_MASS_DRY_AIR_R8 &
               /vnl_varMassFromVarNum(obsoper%constituent_id)

       !case(15192,15011) 
           
          ! Code to be revised when actually applied for the first time
          ! according to model field units.
           
       case default 
        
          call utl_abort('CHM_CONVERT_UNITS: Unknown obs units for varno = ' // trim(utl_str(obsoper%varno)) )
         
       end select
    end if
  
    ! Apply constant scaling
    model_col = model_col * zcoef
    
    if (exp_T.ne.0) then
       if (any(obsoper%tt.le.0.)) call utl_abort("CHM_CONVERT_UNITS: Missing valid temperature for conversion.")
       model_col = model_col * obsoper%tt**exp_T
    end if
    
    if (exp_P.ne.0) model_col = model_col * obsoper%pp**exp_P
    
  end subroutine chm_convert_units

  !--------------------------------------------------------------------------
  ! chm_required_field
  !--------------------------------------------------------------------------
  function chm_required_field(varName,varno) result(needed)
    !
    !:Purpose: To determine whether the specifed field name is required
    !          somewhere in the observation operators for a particular
    !          observation type.
    !
    implicit none
    logical :: needed

    ! Arguments:
    character(len=*), intent(in) :: varName ! Name of field
    integer,          intent(in) :: varno   ! BUFR descriptor element

    select case(trim(varName))
    case('TT')
 
       select case (varno)
       case(BUFR_UNIT_Density,BUFR_UNIT_Density2,BUFR_UNIT_AirDensity, &
            BUFR_UNIT_PMDensity,BUFR_UNIT_NumberDensity,BUFR_UNIT_MolarDensity)
          needed = .true.
       case default
          needed = .false.
       end select
         
    case default
       needed = .false.
    end select

  end function chm_required_field

  !--------------------------------------------------------------------------
  ! chm_layer_integ_operator
  !--------------------------------------------------------------------------
  subroutine chm_layer_integ_operator(obsoper,ixtr,success,kmode)
    !
    !:Purpose: To perform layer integration and calculations
    !
    !:Arguments:
    !     :kmode:       Observation model stage used to allow option of tropo
    !                   increment determination from total column data when
    !                   obsoper%column_bound > obsoper%vlayertop for
    !                   nobslev=1. For kmode=3, calc only for region between
    !                   obsoper%column_bound and surface. For kmode=2, split
    !                   calc for model top to obsoper%column_bound and 
    !                   obsoper%column_bound and surface.
    !
    implicit none

    ! Arguments:
    type(struct_chm_obsoperators), intent(inout) :: obsoper ! observation operator object
    integer, intent(inout) :: ixtr(obsoper%nobslev) ! Flag indicating if obs outside model vertical range: 0 for no.
    logical, intent(inout) :: success(obsoper%nobslev) ! success of integration
    integer, intent(in) :: kmode

    ! Locals:
    integer :: iobslev,tropo_mode
    real(8) :: vlayertop_ref,vlayerbottom_ref,imodlev_bot_ref,checkID
    
    ! Conduct initial setup for vertical integration components

    call chm_vertintg_setup(obsoper)

    ! Ensure that each layer is within model vertical range.
    
    do iobslev=1,obsoper%nobslev
     
       if (obsoper%vlayerbottom(iobslev).lt.obsoper%vlayertop(iobslev)) then
          success(1:obsoper%nobslev)=.false.
          write(*,*) 'chm_layer_integ_operator: WARNING. Layer top/bot value problem.', &
               obsoper%vlayertop(iobslev), obsoper%vlayerbottom(iobslev), &
               '. Entire profile skipped over.'
          return
       else if (obsoper%vlayerbottom(iobslev).lt.obsoper%pp(1)*1.01 .or. &  
            obsoper%vlayertop(iobslev).gt.obsoper%pp(obsoper%nmodlev)*0.99) then
          success(iobslev)=.false.
          if (obsoper%vlayerbottom(iobslev).lt.obsoper%pp(1)*1.01) then
             ixtr(iobslev)=1
          else
             ixtr(iobslev)=2
          end if
          write(*,*) 'chm_layer_integ_operator: WARNING. Layer top/bot value problem.', &
               obsoper%vlayertop(iobslev), obsoper%vlayerbottom(iobslev)
          cycle
       end if
       if (obsoper%vlayerbottom(iobslev).gt.obsoper%pp(obsoper%nmodlev)*0.999) obsoper%vlayerbottom(iobslev)=obsoper%pp(obsoper%nmodlev)*0.999
       if (obsoper%vlayertop(iobslev).lt.obsoper%pp(1)*1.001) obsoper%vlayertop(iobslev)=obsoper%pp(1)*1.001
      
    end do

    ! Check for special treatment if tropo_mode>=1, kmode=2,3, and nobslev=1 for
    ! column observations (obsoper%vco=4)  that extend to the surface.

    tropo_mode = chm_tropo_mode(obsoper%constituent_id)
    
    if (obsoper%nobslev.eq.1.and.kmode.ge.2.and.obsoper%vlayerbottom(1).gt.obsoper%pp(obsoper%nmodlev)*0.99.and. &
         obsoper%constituent_id.ge.0.and.obsoper%vco.eq.4) then
       
       ! Check if constituent id is recognized (function will abort if not recognized)
       if ( obsoper%constituent_id >= 0 ) checkID = vnl_varMassFromVarnum(obsoper%constituent_id)
       
       vlayerbottom_ref=obsoper%vlayerbottom(1)
       
       if ( tropo_mode.ge.1 .and. obsoper%column_bound.gt.obsoper%vlayertop(1) ) then
          
          if (obsoper%iavgkern.ne.0) &
               call utl_abort("chm_layer_integ_operator: Use of averaging kernels not possible with reduced range of increment profile.")
          
          if (kmode.eq.2.and.tropo_mode.eq.1) then
             
             ! When kmode=2, split calc in two. This is done due to difference in 
             ! calc at the interface region when producing zh and zhp. The tangent linear
             ! model in the lower region for kmode=2 must be consistent with 
             ! that associated to kmode=3.
             
             ! Start with bottom region in order to use correct zhp with chm_genoper
             ! when use of this operator is requested.
             
             vlayertop_ref=obsoper%vlayertop(1)
             
             obsoper%vlayertop(1) = obsoper%column_bound
             
             call chm_vertintg(obsoper,ixtr,success)
             
             ! Apply generalized innovation operator if requested
             
             if (obsoper%apply_genoper) call chm_genoper(obsoper,kmode,success)
             
             obsoper%apply_genoper=.false.
             
             imodlev_bot_ref=obsoper%imodlev_bot(1)
             
             ! Reset top and bottom values for integration of the remaining region.
             ! The second integration provides the change in upper level contributions to the
             ! total column from the assimilation of other observations.
             
             obsoper%vlayertop(1)=vlayertop_ref
             obsoper%vlayerbottom(1)=obsoper%column_bound
          else
             ! Reset top new value. Restricts adjoint/tangent linear calcs to this reduced region.            
             obsoper%vlayertop(1)=obsoper%column_bound
          end if
          
       end if
    else
       vlayerbottom_ref=0.0
    end if
    
    ! Calculate vertical integration components for specified obs layer.
    
    call chm_vertintg(obsoper,ixtr,success)
    
    ! If tropo_mode=1, reset original vertical range 
    ! for the tangent linear operator 
    if (obsoper%nobslev.eq.1.and.kmode.eq.2.and.obsoper%constituent_id.ge.0.and.obsoper%vco.eq.4) then
       if (tropo_mode.eq.1.and.  &
            obsoper%column_bound.gt.obsoper%vlayertop(1).and. &
            vlayerbottom_ref.gt.obsoper%pp(obsoper%nmodlev)*0.99) then    
          
          obsoper%vlayerbottom(1)=vlayerbottom_ref
          obsoper%imodlev_bot(1)=imodlev_bot_ref
       end if
    end if
   
  end subroutine chm_layer_integ_operator

  !--------------------------------------------------------------------------
  ! chm_vertintg_setup
  !--------------------------------------------------------------------------
  subroutine chm_vertintg_setup(obsoper)
    !
    !:Purpose: Preliminary calculations for producing components required for
    !          vertical  integration w.r.t. pressure to calculate partial (or
    !          total) column value of model state profile or used for adjoint
    !          calc.
    !
    !:Output:
    !         :obsoper%vmodpress(nmodlev+1):       Model layer boundaries given
    !                                              that pressmod are taken at
    !                                              mid-layer values.
    !         :obsoper%vweights(nmodlev,nmodlev):  Second order Lagrangian
    !                                              interp integration weights
    !
    !:Comments:
    !         - This subroutine does the following:
    !
    !           - Setting of model layer boundaries
    !           - Determining integration weights associated to second order
    !             Lagrangian interpolation.
    !
    !         - Layer boundaries are taken as mid-point between eta levels in
    !           lnP coordinate. Layer values are set to be the values
    !           interpolated to the mid-point in P within the various layers.
    !           Interpolation in P is done quadratically. 
    !
    implicit none

    ! Arguments:
    type(struct_chm_obsoperators), intent(inout) :: obsoper

    ! Locals:

    integer   :: jk
    real(8)   :: zp, zp1, zp2, zp3, zr1, zr2, zr3

    ! Determine P boundaries of analysis layers and save weights for
    ! use in setting innovation operator array.
    ! N.B.: Boundaries of layers set to mid-point of model levels
      
    ! Calculate layer boundaries

    obsoper%vmodpress(1)=obsoper%pp(1)
    obsoper%vmodpress(obsoper%nmodlev+1)= obsoper%pp(obsoper%nmodlev)

    DO JK = 2, obsoper%nmodlev
       obsoper%vmodpress(jk)=sqrt(obsoper%pp(jk-1)*obsoper%pp(jk))
    END DO

    ! Interpolation to mid-layer level in P using
    ! second degree Lagrangian interpolator.
    ! N.B.: Integration is w.r.t. P
    
    ! Calculating for jk=1
    
    zp1= obsoper%pp(1)
    zp2= obsoper%pp(2)
    zp3= obsoper%pp(3)
    zp = (obsoper%vmodpress(2)+obsoper%vmodpress(1))/2.0
    zr1=(zp-zp2)*(zp-zp3)/(zp1-zp2)/(zp1-zp3)
    zr2=(zp-zp1)*(zp-zp3)/(zp2-zp1)/(zp2-zp3)
    zr3=(zp-zp2)*(zp-zp1)/(zp3-zp2)/(zp3-zp1)
    obsoper%vweights(1,1)=zr1
    obsoper%vweights(2,1)=zr2
    obsoper%vweights(3,1)=zr3

    DO JK=2,obsoper%nmodlev-1
       zp1=obsoper%pp(jk-1)
       zp2=obsoper%pp(jk)
       zp3=obsoper%pp(jk+1)
       zp=(obsoper%vmodpress(jk+1)+obsoper%vmodpress(jk))/2.0
       zr1=(zp-zp2)*(zp-zp3)/(zp1-zp2)/(zp1-zp3)
       zr2=(zp-zp1)*(zp-zp3)/(zp2-zp1)/(zp2-zp3)
       zr3=(zp-zp2)*(zp-zp1)/(zp3-zp2)/(zp3-zp1)
       obsoper%vweights(jk-1,jk)=zr1
       obsoper%vweights(jk,jk)=zr2
       obsoper%vweights(jk+1,jk)=zr3
    ENDDO
    
    ! Calculating  for jk=obsoper%nmodlev
    
    zp1= obsoper%pp(obsoper%nmodlev-2)
    zp2= obsoper%pp(obsoper%nmodlev-1)
    zp3= obsoper%pp(obsoper%nmodlev)
    zp = (obsoper%vmodpress(obsoper%nmodlev+1)+obsoper%vmodpress(obsoper%nmodlev))/2.0
    zr1=(zp-zp2)*(zp-zp3)/(zp1-zp2)/(zp1-zp3)
    zr2=(zp-zp1)*(zp-zp3)/(zp2-zp1)/(zp2-zp3)
    zr3=(zp-zp2)*(zp-zp1)/(zp3-zp2)/(zp3-zp1)
    obsoper%vweights(obsoper%nmodlev-2,obsoper%nmodlev)=zr1
    obsoper%vweights(obsoper%nmodlev-1,obsoper%nmodlev)=zr2
    obsoper%vweights(obsoper%nmodlev,obsoper%nmodlev)=zr3

  end subroutine chm_vertintg_setup

  !--------------------------------------------------------------------------
  ! chm_vertintg
  !--------------------------------------------------------------------------
  subroutine chm_vertintg(obsoper,ixtr,success)
    !
    !:Purpose: To calculate components required for vertical integration w.r.t.
    !          pressure to calculate partial (or total) column value of model
    !          state profile or used for adjoint calc.
    !
    !:Input:
    !         :obsoper%nmodlev:      # of model vertical levels
    !         :obsoper%nobslev:      # of obs vertical levels
    !         :obsoper%vweights:     See routine chm_vertintg_setup
    !         :obsoper%vmodpress:
    !
    !:Output:
    !         :obsoper%zh(obsoper%nobslev,obsoper%nmodlev):
    !                      Initial innovation model array (other than
    !                      conversion constants)
    !         :obsoper%zhp(obsoper%nobslev,obsoper%nmodlev): 
    !                      Part of innovation operator not related to resolution
    !
    implicit none

    ! Arguments:
    type(struct_chm_obsoperators), intent(inout) :: obsoper
    integer, intent(in) ::ixtr(obsoper%nobslev) ! Flag indicating if obs outside model vertical range (0 for no)
    logical, intent(in) :: success(obsoper%nobslev) ! indicating if calc are to be performed

    ! Locals:
    integer, parameter :: ivweights=2  ! Order of Lagrangian interpolation.

    ! Declaration of local variables
    
    integer   :: J,JK,ILMAX2,ILMIN2
    integer   :: ILMIN, ILMAX, iobslev
    real(8)   :: zp, zp1, zp2, zp3, zr1, zr2, zr3, ptop, pbtm

    do iobslev=1,obsoper%nobslev

       if (success(iobslev).or.(ixtr(iobslev).eq.0.and.obsoper%iavgkern.ne.0)) then

          ptop = obsoper%vlayertop(iobslev)
          pbtm = obsoper%vlayerbottom(iobslev)
         
          ! Find the range of vertical levels over which to perform the integration
          ! and set innovation operator ZH over this range.
          
          ilmin=1
          ilmax=obsoper%nmodlev
          if (ptop.le.obsoper%vmodpress(1)*1.01.and.pbtm.ge.obsoper%vmodpress(obsoper%nmodlev+1)*0.99) then

             ! Total column integration part
             
             do jk = 1,obsoper%nmodlev
                do j=max(1,jk-ivweights),min(obsoper%nmodlev,jk+ivweights)
                   obsoper%zh(iobslev,jk)=obsoper%zh(iobslev,jk)+(obsoper%vmodpress(j+1) &
                        -obsoper%vmodpress(j))*obsoper%vweights(jk,j)
                   obsoper%zhp(iobslev,jk)=obsoper%zhp(iobslev,jk)+obsoper%vweights(jk,j)
                end do
             end do
             
          else

             ! Partial column integration part (special treatment at boundaries)
     
             ! Identify analysis layer boundaries just within obs layer.
             
             ilmin = chm_igetmodlev(ptop, obsoper%vmodpress, 'top', obsoper%nmodlev+1)
             ilmax = chm_igetmodlev(pbtm, obsoper%vmodpress, 'btm', obsoper%nmodlev+1)
               
             if (ilmin.eq.ilmax+1) then

                ! Entire obs layer within one analysis layer
                
                j=ilmin
                if (j.lt.3) j=3
                if (j.gt.obsoper%nmodlev) j=obsoper%nmodlev
                zp1=obsoper%vmodpress(j-2)
                zp2=obsoper%vmodpress(j-1)
                zp3=obsoper%vmodpress(j)
                zp=(ptop+pbtm)/2.0
                zr1=(zp-zp2)*(zp-zp3)/(zp1-zp2)/(zp1-zp3)
                zr2=(zp-zp1)*(zp-zp3)/(zp2-zp1)/(zp2-zp3)
                zr3=(zp-zp2)*(zp-zp1)/(zp3-zp2)/(zp3-zp1)
                
                obsoper%zh(iobslev,j-2)=(pbtm-ptop)*zr1
                obsoper%zh(iobslev,j-1)=(pbtm-ptop)*zr2
                obsoper%zh(iobslev,j)=(pbtm-ptop)*zr3
                obsoper%zhp(iobslev,j-2)=zr1
                obsoper%zhp(iobslev,j-1)=zr2
                obsoper%zhp(iobslev,j)=zr3
                ilmin=j-2
                ilmax=j
                  
             else
                
                ! Determine terms from the inner layers (excluding the lower and upper
                ! boundary layers when these layers not covering entire analyses layers)
                
                if (pbtm.ge.obsoper%vmodpress(obsoper%nmodlev)*0.99) then
                   ilmax2=obsoper%nmodlev
                else
                   ilmax2=ilmax-1
                end if
                if (ptop.le.obsoper%vmodpress(1)*1.01) then
                   ilmin=1
                   ilmin2=ilmin
                else
                   ilmin2=ilmin
                end if
                if (ilmin2.le.ilmax2) then
                   do jk = ilmin2,ilmax2
                      do j=max(1,jk-ivweights),min(obsoper%nmodlev,jk+ivweights)
                         obsoper%zh(iobslev,jk)=obsoper%zh(iobslev,jk)+(obsoper%vmodpress(j+1) &
                               -obsoper%vmodpress(j))*obsoper%vweights(jk,j)
                         obsoper%zhp(iobslev,jk)=obsoper%zhp(iobslev,jk)+obsoper%vweights(jk,j)
                      end do
                   end do
                end if
                
                ! Determine terms from the lower and upper boundary layers
                ! when these layers do not cover entire analyses layers.
                
                if (pbtm.lt.obsoper%vmodpress(obsoper%nmodlev)*0.99) then
                     
                   j=ilmax+1
                   if (j.gt.obsoper%nmodlev) j=obsoper%nmodlev
                   if (j.lt.3) j=3
                   zp1=obsoper%vmodpress(j-2)
                   zp2=obsoper%vmodpress(j-1)
                   zp3=obsoper%vmodpress(j)
                   zp=(obsoper%vmodpress(ilmax)+pbtm)/2.0
                   zr1=(zp-zp2)*(zp-zp3)/(zp1-zp2)/(zp1-zp3)
                   zr2=(zp-zp1)*(zp-zp3)/(zp2-zp1)/(zp2-zp3)
                   zr3=(zp-zp2)*(zp-zp1)/(zp3-zp2)/(zp3-zp1)
                   
                   obsoper%zh(iobslev,j-2)=obsoper%zh(iobslev,j-2)+(pbtm - obsoper%vmodpress(ilmax))*zr1
                   obsoper%zh(iobslev,j-1)=obsoper%zh(iobslev,j-1)+(pbtm - obsoper%vmodpress(ilmax))*zr2
                   obsoper%zh(iobslev,j)=obsoper%zh(iobslev,j)+(pbtm - obsoper%vmodpress(ilmax))*zr3
                   obsoper%zhp(iobslev,j-2)=obsoper%zhp(iobslev,j-2)+zr1
                   obsoper%zhp(iobslev,j-1)=obsoper%zhp(iobslev,j-1)+zr2
                   obsoper%zhp(iobslev,j)=obsoper%zhp(iobslev,j)+zr3
                   ilmax=j
                  
                end if
                  
                if (ptop.gt.obsoper%vmodpress(1)*1.01) then
                     
                   j=ilmin-1
                   if (j.lt.1) j=1
                   if (j.gt.obsoper%nmodlev-2) j=obsoper%nmodlev-2
                   zp1= obsoper%vmodpress(j)
                   zp2= obsoper%vmodpress(j+1)
                   zp3= obsoper%vmodpress(j+2)
                   zp = (obsoper%vmodpress(ilmin)+ptop)/2.0
                   zr1=(zp-zp2)*(zp-zp3)/(zp1-zp2)/(zp1-zp3)
                   zr2=(zp-zp1)*(zp-zp3)/(zp2-zp1)/(zp2-zp3)
                   zr3=(zp-zp2)*(zp-zp1)/(zp3-zp2)/(zp3-zp1)
                   
                   obsoper%zh(iobslev,j)=obsoper%zh(iobslev,j)+(obsoper%vmodpress(ilmin)-ptop)*zr1
                   obsoper%zh(iobslev,j+1)=obsoper%zh(iobslev,j+1)+(obsoper%vmodpress(ilmin)-ptop)*zr2
                   obsoper%zh(iobslev,j+2)=obsoper%zh(iobslev,j+2)+(obsoper%vmodpress(ilmin)-ptop)*zr3
                   obsoper%zhp(iobslev,j)=obsoper%zhp(iobslev,j)+zr1
                   obsoper%zhp(iobslev,j+1)=obsoper%zhp(iobslev,j+1)+zr2
                   obsoper%zhp(iobslev,j+2)=obsoper%zhp(iobslev,j+2)+zr3
                   ilmin=j
                   if (ilmax.lt.j+2) ilmax=j+2
                   
                end if
                if (ilmin.gt.ilmax-2) ilmin=ilmax-2
             end if
          end if

          obsoper%imodlev_top(iobslev)=ilmin
          obsoper%imodlev_bot(iobslev)=ilmax
            
       else
          obsoper%zh(iobslev,:) = 0.0D0
          obsoper%zhp(iobslev,:) = 0.0D0
          
          obsoper%imodlev_top(iobslev)=1
          obsoper%imodlev_bot(iobslev)=1
       end if
       
    end do

  end subroutine chm_vertintg

  !--------------------------------------------------------------------------
  ! chm_vert_interp_operator
  !--------------------------------------------------------------------------
  subroutine chm_vert_interp_operator(obsoper,pres_obs,ixtr,success)
    !
    !:Purpose: Interpolation to point in profile. Uses piecewise linear vertical
    !          interpolation in log(Pressure).
    !
    !:Input:
    !       :obsoper%pp:     pressure on model levels, assumed to be in
    !                        ascending order
    !
    !:Output:
    !
    !       :obsoper%zh:     interpolation coefficients
    !
    !:Comments:
    !       - Current implementation searches for index of nearest model level.
    !         This step is redundant since this information is already saved in
    !         obsSpaceData in OBS_LYR. This step is repeated so that the
    !         routines used under oopc_CHobsoperators are more independent of the
    !         rest of the EnVar code. If it is desired to skip this redundant
    !         step, the content of the OBS_LYR column could be passed to
    !         chm_obsoperators and subsequently to this subroutine.
    !
    implicit none

    ! Arguments:
    type(struct_chm_obsoperators), intent(inout) :: obsoper
    integer, intent(inout) :: ixtr(obsoper%nobslev)  ! Flag indicating if obs outside model vertical range (0 for no)
    real(8), intent(in) :: pres_obs(obsoper%nobslev) ! pressure on observation levels
    logical, intent(inout) :: success(obsoper%nobslev)! success of interpolation

    ! Locals:
    integer :: iobslev,jmodlev

    do iobslev=1,obsoper%nobslev

       ! check if obs is above or below model boundaries
       if ( pres_obs(iobslev).lt.obsoper%pp(1) .or. &
            pres_obs(iobslev).gt.obsoper%pp(obsoper%nmodlev) ) then
          success(iobslev)=.false.
          if (pres_obs(iobslev).lt.obsoper%pp(1)) then
              ixtr(iobslev)=1
          else
              ixtr(iobslev)=2
          end if 
       end if

       if (success(iobslev).or.(ixtr(iobslev).eq.0.and.obsoper%iavgkern.ne.0)) then

          ! Find model layers directly above and below obs.
          ! After exit of loop, the obs will be between model
          ! levels jmodlev and jmodlev+1.
          do jmodlev=1,obsoper%nmodlev-1
             if ( pres_obs(iobslev).ge.obsoper%pp(jmodlev) .and. &
                  pres_obs(iobslev).lt.obsoper%pp(jmodlev+1) ) then
                exit
             end if
          end do

          ! Set interpolation weights
          obsoper%zh(iobslev,jmodlev+1) = LOG(pres_obs(iobslev)/obsoper%pp(jmodlev)) &
                                        / LOG(obsoper%pp(jmodlev+1)/obsoper%pp(jmodlev))
          obsoper%zh(iobslev,jmodlev) = 1.0D0 - obsoper%zh(iobslev,jmodlev+1)

          ! set range of nonzero elements for model vertical levels
          obsoper%imodlev_top(iobslev) = jmodlev
          obsoper%imodlev_bot(iobslev) = jmodlev+1

       else
          obsoper%imodlev_top(iobslev) = 1
          obsoper%imodlev_bot(iobslev) = 1
       end if

    end do

  end subroutine chm_vert_interp_operator

  !--------------------------------------------------------------------------
  ! chm_igetmodlev
  !--------------------------------------------------------------------------
  integer function chm_igetmodlev(rpress, rppobs, topbtm, ntotlev)
    !
    !:Purpose: To get the vertical level index for the pressure in rppobs
    !          within obs layer and nearest specified obs layer boundary.  
    implicit none

    ! Arguments:
    real(8), intent(in) :: rpress          ! pressure value in Pascal
    real(8), intent(in) :: rppobs(ntotlev) ! profile of pressure at obs location
    character(len=*), intent(in) :: topbtm ! indicating whether we are looking for top or bottom pressure
    integer, intent(in) :: ntotlev         ! total number of levels of rppobs

    ! Locals
    integer     :: ilev1, ilev2
    integer     :: jk

    ! Find the model levels adjacent to pressure level rpress

    ! Default values
    
    if (rpress .lt. 0.) then
       if ((topbtm .eq. 'btm') .or. (topbtm .eq. 'BTM')) then
          chm_igetmodlev = ntotlev
       endif
       if ((topbtm .eq. 'top') .or. (topbtm .eq. 'TOP')) then
          chm_igetmodlev = 1
       endif                                                  
    endif
      
    ilev1=0
    ilev2=1
    do jk=1,ntotlev
       if (rpress.gt.rppobs(jk)) then
          ilev1=jk
          ilev2=jk+1
       else
          exit
       endif
    enddo

    ! Find the model level index

    ! If we are looking for top level, the index is the level immediately 
    ! below. if looking for bottom level, the index is the one immediately 
    ! above.
    
    if ((topbtm .eq. 'btm') .or. (topbtm .eq. 'BTM')) then
       chm_igetmodlev=ilev1
    else if ((topbtm .eq. 'top') .or. (topbtm .eq. 'TOP')) then
       chm_igetmodlev=ilev2
    endif

    if (chm_igetmodlev .lt. 1) chm_igetmodlev=1
    if (chm_igetmodlev .gt. ntotlev) chm_igetmodlev=ntotlev
  
  end function chm_igetmodlev

  !--------------------------------------------------------------------------
  ! chm_genoper
  !--------------------------------------------------------------------------
  subroutine chm_genoper(obsoper,kmode,success)
    !
    !:Purpose: Set generalized innovation operator for integral or layer avg
    !          obs. Relevant only for 3D incremental fields. This version is
    !          intended to vertically distribute the obs increments
    !          proportionally to the background state.
    !
    !:Input:
    !       :obsoper%zh,zhp:       see routine chm_vertintg_setup
    !
    !:Output:
    !       :obsoper%zh(obsoper%nmodlev):     a*w: Final innovation model array
    !                                         (other than conversion constants)
    !       :obsoper%zhp(obsoper%nmodlev):    w (see comments section)
    !
    ! Comments:
    !
    !      (1) This routine prepares an alternative innovation operator g,
    !          called the generalized innovation operator, to take the place of
    !          the innovation (TLM) operator h (row of zh). The operator g is
    !          specified as:
    ! 
    !             g = a*w
    ! 
    !          where the modified innovation operator 'w' can be set as:
    ! 
    !             w = P[ (h'x)^T ] *  B^{-1}     
    ! 
    !          with  h' is the part of h which excludes resolution dependence
    !          (only/mostly contains the physics part of h; zhp),
    !
    !            - x is the state profile rval
    !            - P is a window cutoff operator (sets small values to zero)
    !            - B is the original/initial total "vertical" covariance matrix
    !              (in 2D)
    !          
    ! 
    !          and 'a' is a proportionality constant ensuring that the
    !          innovation increment remains unchanged for the 1D case in the
    !          absence of other obs., i.e.,
    ! 
    !                  a^2 = (h*B*h^T)(w*B*w^T)^{-1},
    ! 
    !          Application of the state profile x (rval) is to make the
    !          increment profile be more proportional to the state profile.
    ! 
    !          The presence of B^{-1} is to negate the weight re-distribution
    !          from the later application of B in grad(Jo).
    ! 
    !          While dx is provided to the obs operator, the minimization is
    !          done for dx/sigma where sigma is the background error std. dev.
    !          in B and so C (correlation matrix) is used instead of B in the
    !          minimization. Moreover, the transformation from dx/sigma to dx is
    !          done outside the forward model operators (at the spectral to
    !          physical space  transformation step). For this reason, the
    !          expression for w should technically be replaced by
    !  
    !             w = P[ (h'x/sigma^2)^T ] *  C^{-1}            (Option 1 below)
    ! 
    !          still with
    ! 
    !             a^2 = (h*B*h^T)(w*B*w^T)^{-1}
    !  
    !          The presence of C^{-1} does/can give difficulty to the iterative
    !          variational  minimization. It can result in oscillations in the
    !          increment profile depending on where the iterations are stopped.
    !          Moreover, if the spectral space C matrix is for  non-separable
    !          vertical and horizontal correlations, there will be oscillations 
    !          due to inconstencies between the total inverse vertical
    !          correlation C^{-1} in physical space  and the inverse vertical
    !          correlation matrix for each spectral wavenumber.
    ! 
    !          As alternatives, one can completely omit the role of  C^{-1} from
    !          w, i.e.
    ! 
    !             w = P[ (h'x/sigma^2)^T ]                        (Option 3)
    ! 
    !          or use the following substitute to approximate the role of C^{-1}
    !          in approximately negating the weight re-distribution from later
    !          application of C in grad(Jo), i.e.
    ! 
    !             w(i) = P[ (h'x/sigma^2)^T ]_i / sum(C(:,i))  (Option 2 - preferred)
    ! 
    !      (2) The matrices B and B^{-1} are the total error covariance matrix
    !          (in physical space)
    !          and its inverse with the related error correlation matrices
    !          *corvert* and *corverti*  provided from *bmatrixchem_mod.ftn90*.
    ! 
    !      (3) In the presence of both (a) neighbouring measurements and (b)
    !          horizontal background error  correlation lengths that vary in the
    !          vertical, the increments will also be subject to the latter,
    !          displaying larger increments in vertical regions with larger
    !          horizontal error correlation lengths - this distorting the
    !          vertical increment distribution stemming from chm_genoper alone.
    !          This stems from w accounting for vertical correlations via
    !          C^{^-1} and not the horizontal correlations in B.
    ! 
    !      (4) NOTE: Cases with ensemble-based and or lam-based background
    !          covariances are not taken into account in this version.
    !
    implicit none

    ! Arguments:
    type(struct_chm_obsoperators), intent(inout) :: obsoper
    integer, intent(in) :: kmode ! index specifying if content to be applied (i.e. if kmode>1)
    logical, intent(in) :: success(obsoper%nobslev) ! indicating if calc are to be performed

    ! Locals:
    real(8), parameter :: pwin=0.01
    integer  :: iobslev,imodlev,irmse,jvar
    logical  :: lrgsig
    real(8)  :: zwbw(1),zhbh(1),za,work(obsoper%nmodlev),sigma_trial(obsoper%nmodlev,4)
    real(8), parameter :: threshold=1.D-20
    real(8)  :: zmin,rvalw(obsoper%nmodlev),rvalr(obsoper%nmodlev)
    real(8)  :: rvalc(obsoper%nmodlev),rmse,w1 
    character(len=22) :: code
   
    if (kmode.le.1) return

    ! Retrieve from stored background error std dev [elemements (:,1-2)] at obs location [and inverses at elements (:,3-4)]    
    sigma_trial=chm_retrieve_sigma_trial(obsoper%obs_index)               

    ! Identify variable position index in background error correlation matrices
    
    jvar=1
    do while (trim(bchm_varnamelist(jvar)).ne.'') 
       if (trim(bchm_varnamelist(jvar)).eq.trim(obsoper%varName)) exit
       jvar=jvar+1
    end do
    if (trim(bchm_varnamelist(jvar)).eq.'') call utl_abort('chm_genoper: Correlation matrix not found for ' // trim(obsoper%varName) )
    
    ! Initialize reference mass (mixing ratio) weighting profile as profile from trial field
     
    rvalr(1:obsoper%nmodlev)=obsoper%trial(1:obsoper%nmodlev)
    
    ! Check on rvalr
            
    if (any(abs(rvalr(:)).le.threshold)) then
       ! write(*,*) 'chm_genoper: Trial field profile segment is ~zero. Could cause an abort in this routine.'
       ! write(*,*) 'chm_genoper: To prevent abort and still be effective, corresponding elements set to a non-zero constant.'
       ! write(*,*) 'chm_genoper: ',varName,rlat,rlong,nobslev
       if (all(abs(rvalr(:)).le.threshold)) then
          rvalr(:)=1.0 
       else  
          zmin=minval(abs(rvalr(:)), &
               mask=abs(rvalr(:)).gt.threshold)
          where (abs(rvalr(:)).le.threshold) &
               rvalr(:)=zmin  
       end if
    end if
    
    ! Loop over obs elements
    
    lrgsig=.true. 
    write(code,'(I22)') obsoper%obs_index
    
    do iobslev=1,obsoper%nobslev
       
       if (.not.success(iobslev)) cycle   
       
       ! Set reference mass (mixing ratio) weighting profile
       
       rvalw(1:obsoper%nmodlev)=rvalr(1:obsoper%nmodlev)
       if (chm_generalized_operator(obsoper%constituent_id) > 1) then

          ! Set reference mass weighting profile according to the difference between 
          ! an external reference (such as a climatology) and trial field profiles.
          !
          ! This is a mechanism to force the solution profile shape somewhat towards that 
          ! of the external reference when the reliability of the vertical structure of the 
          ! trial field is not high.
          !
          ! This can be used at the beginning of long assimilation periods if the initial
          ! trial field is not as realistic as may be desired or somewhat mitigate
          ! gradual biasing of vertical structures that might otherwise occur from assimilation 
          ! of integrated quantities (when there is insufficient data from other observation types).
          !
          ! The larger the rms difference of xc (external reference) with xb (trial field profile), 
          ! the greater is the influence of this difference in the weighting. If there is 
          ! little difference, then either xb or xc can be directly used as the weighting 
          ! profile; xb is used below.
          
          rmse=0.0
          irmse=0
          rvalc(1:obsoper%nmodlev)=chm_get_ref_column(code)
            
          zmin=pwin*maxval(abs(obsoper%zhp(iobslev,1:obsoper%nmodlev)))
          do imodlev=1,obsoper%nmodlev
             if (obsoper%zhp(iobslev,imodlev).gt.zmin) then
                irmse=irmse+1
                rmse=rmse+((rvalc(imodlev)-obsoper%trial(imodlev))*sigma_trial(imodlev,3))**2
             end if
          end do
          rmse=rmse*2./irmse
          if (rmse.lt.1.0) then
             rvalw(1:obsoper%nmodlev)=rvalc(1:obsoper%nmodlev)
          else 
             rvalw(1:obsoper%nmodlev)=rvalc(1:obsoper%nmodlev)-obsoper%trial(1:obsoper%nmodlev)
             if (rmse.lt.5.0) then
                w1=(5.0-rmse)/(5.0-1.0)
                rvalw(1:obsoper%nmodlev)=(1.0-w1)*rvalw(1:obsoper%nmodlev)+w1*obsoper%trial(1:obsoper%nmodlev)
             end if
             where (abs(rvalw(1:obsoper%nmodlev)).lt.0.01*rvalc(1:obsoper%nmodlev)) &
                  rvalw(1:obsoper%nmodlev)=sign(0.01*rvalc(1:obsoper%nmodlev),rvalw(1:obsoper%nmodlev))
          end if
       end if
       
       ! Begin preparation of the new innovation operator w (=new zhp)
       
       if (obsoper%nobslev.eq.1.and.obsoper%imodlev_top(iobslev).eq.1.and. &
            obsoper%imodlev_bot(iobslev).eq.obsoper%nmodlev) then
          
          ! Treat as total column obs. Here, zhp would be approx. equal
          ! to 1 except for the near-end points of the model vertical domain,
          ! the latter due to the discretized domain. Not using zhp avoids this
          ! discretization issue from weakly affecting results at the boundaries.
          
          work(1:obsoper%nmodlev)=obsoper%trial(1:obsoper%nmodlev)*sigma_trial(1:obsoper%nmodlev,3) 
       else 

          ! Account for localized obs function (e.g. partial columns, Jacobians. For Jacobians,
          ! zhp must also be independent of the model layer thicknesses.)
          work(1:obsoper%nmodlev)=obsoper%zhp(iobslev,1:obsoper%nmodlev)* &
               obsoper%trial(1:obsoper%nmodlev)*sigma_trial(1:obsoper%nmodlev,3) 
       end if

       ! Apply cutoff (apply to zhp*obsoper%trial/sigma_trial(:,1) instead of the resultant zh)
       ! Resultant outside cutoff region should be zhp*obsoper%trial/sigma_trial(:,1)/sigma_trial(:,1).
       ! Note: sigma_trial(:,3)=1.0/sigma_trial(:,1)
 
       zmin=pwin*maxval(abs(work(1:obsoper%nmodlev)))
       where (abs(work(1:obsoper%nmodlev)).lt.zmin) 
          work(1:obsoper%nmodlev)=0.0D0 
       elsewhere        
          work(1:obsoper%nmodlev)=work(1:obsoper%nmodlev)*sigma_trial(1:obsoper%nmodlev,3)
       endwhere

       ! Application of C^{-1} or substitute (for negating the weight impact of the
       ! later application of C in finalizing grad(Jo). Option 2 is favoured.
       !  
       ! Option 1: Application of C^{-1}
       ! call chm_corvert_mult(obsoper%varName,work(1:obsoper%nmodlev),obsoper%zhp(iobslev,1:obsoper%nmodlev), &
       !                       obsoper%imodlev_top(iobslev),obsoper%imodlev_bot(iobslev),1,obsoper%nmodlev,obsoper%nmodlev, &
       !                       .false.,-1)
       !
       ! Option 2: Application of 1/sum(C(:,j)) to approximately negate the weight re-distribution from C
       !           in the calc of grad(Jo).
       !
       ! call chm_corvert_mult(obsoper%varName,work(1:obsoper%nmodlev),obsoper%zhp(iobslev,1:obsoper%nmodlev), &
       !                      obsoper%imodlev_top(iobslev),obsoper%imodlev_bot(iobslev),1,obsoper%nmodlev,obsoper%nmodlev, &
       !                      .false.,0)
       obsoper%zhp(iobslev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev))= &
            obsoper%zhp(iobslev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)) &
            *bchm_invsum(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev),jvar)

       ! Option 3: Just skip over consideration attempt at negating the weight re-distribution from C.
       ! chm_obsoper%zhp(iobslev,1:obsoper%nmodlev)=work(1:obsoper%nmodlev)
       !
       ! Determine proportionality factor 'a' = (h*B*h^T)(w*B*w^T)^{-1}
       !
       ! Determine/estimate w*B*w^T (zwbw(1))
       !
       ! call chm_corvert_mult(obsoper%varName,obsoper%zhp(iobslev,1:obsoper%nmodlev), &
       !                       zwbw,obsoper%imodlev_top(iobslev),obsoper%imodlev_bot(iobslev),1,obsoper%nmodlev, &
       !                       1,lrgsig,3,sigma_trial))
       do imodlev=obsoper%imodlev_top(iobslev),obsoper%imodlev_bot(iobslev)
          work(imodlev)=sum(obsoper%zhp(iobslev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)) &
               *bchm_corvert(imodlev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev),jvar) &
               *sigma_trial(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev),1))*sigma_trial(imodlev,1)
       end do
       zwbw(1)=sum(obsoper%zhp(iobslev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)) &
            *work(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)))
       
       ! Determine/estimate h*B*h^T (zhbh(1))
       !
       !  call chm_corvert_mult(obsoper%varName,obsoper%zh(iobslev,1:obsoper%nmodlev), &
       !                        zhbh,obsoper%imodlev_top(iobslev),obsoper%imodlev_bot(iobslev),1,obsoper%nmodlev, &
       !                        1,lrgsig,3,sigma_trial)
       do imodlev=obsoper%imodlev_top(iobslev),obsoper%imodlev_bot(iobslev)
          work(imodlev)=sum(obsoper%zh(iobslev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)) &
                       *bchm_corvert(imodlev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev),jvar) &
                       *sigma_trial(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev),1))*sigma_trial(imodlev,1)
       end do
       zhbh(1)=sum(obsoper%zh(iobslev,obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)) &
            *work(obsoper%imodlev_top(iobslev):obsoper%imodlev_bot(iobslev)))

       ! Set proportionality factor 'a'
         
       za=sqrt(zhbh(1)/zwbw(1))
       ! if (abs(obsoper%lat*180./3.1415-78.).lt.2.0.and.abs(rlon*180./3.1415-185.).lt.2.0) then
       !     write(6,*) 'ZA  ',obsoper%lat*180.0/3.1415,rlon*180.0/3.1415,za,zhbh(1),zwbw(1)
       !     write(6,*) 'obsoper%trial',obsoper%trial(1:obsoper%nmodlev)
       !     write(6,*) 'sigma_trial',sigma_trial(1:obsoper%nmodlev,1)
       !     write(6,*) 'ZH  ',chm_obsoper%zh(iobslev,1:obsoper%nmodlev)
       !     write(6,*) 'ZHP ',chm_obsoper%zhp(iobslev,1:obsoper%nmodlev)*za
       ! end if
         
       ! Set final innovation operator
         
       obsoper%zh(iobslev,1:obsoper%nmodlev)=obsoper%zhp(iobslev,1:obsoper%nmodlev)*za
         
    end do

      
  end subroutine chm_genoper

  !--------------------------------------------------------------------------
  ! chm_corvert_mult
  !--------------------------------------------------------------------------
  subroutine chm_corvert_mult(varName,rmat_in,rmat_out,imodlev_top,imodlev_bot,&
                              ndim1,ndim2,ndim3,lrgsig,itype,rsig_opt)
    !
    !:Purpose: To multiply a matrix by the covariance matrix or its inverse
    !          (or combinations of these).
    !
    !:Arguments:
    !    :ndim3:         output dimension expected by the calling routine
    !
    !                        - for itype=+/-3  ndim3=ndim1
    !                        - otherwise       ndim3=ndim2 
    !
    !    :itype:         type of multiplication, for an input matrix A the
    !                    output matrix is
    !                      === ===========================
    !                        0 D(i,j)=A(i,j)/sum(C(1:n,i))
    !                        1 A*C
    !                        2 C*A
    !                        3 A*C*A^T
    !                       -1 A*CI
    !                       -2 CI*A
    !                       -3 A*CI*A^T
    !                      === ===========================
    !                    where C is the covariance matrix and CI is its inverse
    !
    !:Comments:
    !
    !    - If rmat_in is a 1-D vector, then
    !
    !         - for cases +/- 2, one should have set ndim2=1 and
    !           ndim1=vector-length.
    !         - for cases +/- 1,3, one should have set ndim1=1 and
    !           ndim2=vector-length.
    ! 
    !    - Revisions required whem LAM and ensembles cases become available.
    !
    implicit none

    ! Arguments:
      character(len=*), intent(in) :: varName ! variable name
      real(8), intent(in)    :: rmat_in(ndim1,ndim2) ! input matrix/vector A (see comments section)
      real(8), intent(out)   :: rmat_out(ndim1,ndim3) ! output matrix/vector
      integer, intent(in)    :: imodlev_top(ndim1) ! top level of non-zero values in rmat_in
      integer, intent(in)    :: imodlev_bot(ndim1) ! bottom level of non-zero values in rmat_in
      integer, intent(in)    :: ndim1 ! matrix dimension (one 1D input vector)
      integer, intent(in)    :: ndim2 ! matrix dimension
      integer, intent(in)    :: ndim3
      logical, intent(in)    :: lrgsig ! indicates whether rgsig to be included as part of C or CI below.
      integer, intent(in)    :: itype
      real(8), intent(in), optional :: rsig_opt(ndim2,2) ! background error std. dev. at obs locations (must be provided when lrgsig=.true.)

    ! Locals:
      real(8) :: rsig(ndim2,2)
      
      rmat_out(:,:)=0.0D0
      rsig=0.0
      if (present(rsig_opt)) rsig=rsig_opt
        
      ! Apply operation related to static background error covariance/correlation matrix.
      ! Applicability tests within b*chm_corvert_mult.

      call bchm_corvert_mult(varName,rmat_in,rmat_out,imodlev_top,imodlev_bot,ndim1,ndim2,ndim3, &
                             lrgsig,itype,rsig(:,1))
!      call blamchm_corvert_mult(varName,rmat_in,rmat_out,imodlev_top,imodlev_bot,ndim1,ndim2,ndim3, &
!                             lrgsig,itype,rsig(:,1))

      ! Apply operation related to ensemble-based background error covariance/correlation matrix.
      ! Applicability test within benschm_corvert_mult.

!      call benschm_corvert_mult(varName,rmat_in,rmat_out,imodlev_top,imodlev_bot,ndim1,ndim2,ndim3, &
!                                lrgsig,itype,rsig(:,2))

  end subroutine chm_corvert_mult

  function chm_get_col_boundary(iconstituent_id,nmodlev,pressmod,tt,height,hu_opt, &
                                uu_opt,vv_opt) result(bound_press)
    !
    !:Purpose: To determine and to store the boundary (e.g. tropopause or PBL)
    !          pressure levels if needed by the observation operators.    
    implicit none
    real(8) :: bound_press ! pressure level of boundary to be imposed

    ! Arguments:
    integer, intent(in) :: iconstituent_id  ! BUFR code element of Table 08046 identifying the constituent
    integer, intent(in) :: nmodlev          ! number of model levels for variables other than uu and vv
    real(8), intent(in) :: pressmod(nmodlev)! model pressure array
    real(8), intent(in) :: tt(nmodlev)      ! model temperature (Kelvin)
    real(8), intent(in) :: height(nmodlev)  ! height (meters)
    real(8), optional, intent(in) :: hu_opt(nmodlev) ! specific humidity
    real(8), optional, intent(in) :: uu_opt(:) ! model zonal wind component(m/s)
    real(8), optional, intent(in) :: vv_opt(:) ! model meridional wind component (m/s)
   
    ! Locals:
    integer :: tropo_bound
    
    bound_press = -1.
    
    if (chm_tropo_mode(iconstituent_id) == 0) return
   
    tropo_bound = chm_tropo_bound(iconstituent_id)
    
    if (tropo_bound.gt.0) then
       if (.not.all(tt.lt.0.) .and. .not.all(height.lt.0.) ) then
    
          select case(tropo_bound)
          case(1)
    
             ! Get tropopause pressure level
      
             bound_press = phf_get_tropopause(nmodlev,pressmod,tt,height,hu_opt=hu_opt)
    
          case(2)
 
             ! Get PBL pressure level
      
             bound_press = phf_get_pbl(nmodlev,pressmod,tt,height,hu_opt=hu_opt,uu_opt=uu_opt,vv_opt=vv_opt) 
      
          case default
             call utl_abort("chm_get_col_boundary: Unrecognized value for tropo_bound of " &
                  // trim(utl_str(tropo_bound)) )
          end select
                                     
      end if     
    end if
      
    ! Use tropo_column_top value if tropo_bound=0 or model derived boundary was unsuccessful
    if (bound_press.lt.0.0) &
         bound_press = chm_tropo_column_top(iconstituent_id)
      
  end function chm_get_col_boundary
       
  !--------------------------------------------------------------------------
  ! chm_add_col_boundary
  !--------------------------------------------------------------------------
  subroutine chm_add_col_boundary(headerIndex,bound_press)
    !
    !:Purpose: To add column boundary data to chm_column_boundary which can be
    !          retrieved later using a header index.
    !
    implicit none 

    ! Arguments:
    integer, intent(in) :: headerIndex
    real(8), intent(in) :: bound_press

    ! Locals:
    integer obsdata_maxsize

    obsdata_maxsize = chm_obsdata_maxsize
     
    if (.not.associated(chm_column_boundary%data1d)) then
       call oss_obsdata_alloc(chm_column_boundary,obsdata_maxsize, dim1=1)
       chm_column_boundary%nrep = 0
    end if

    ! In this case nrep will count the number of filled reps in the data arrays
    chm_column_boundary%nrep = chm_column_boundary%nrep+1 

    if (chm_column_boundary%nrep.gt.obsdata_maxsize) &
         call utl_abort('chm_add_col_boundary: Reach max size of array ' // trim(utl_str(obsdata_maxsize)) )
  
    ! Use the header number as the unique code for this obs data
    write(chm_column_boundary%code(chm_column_boundary%nrep),'(I22)') headerIndex

    chm_column_boundary%data1d(1,chm_column_boundary%nrep) = bound_press

  end subroutine chm_add_col_boundary

  !--------------------------------------------------------------------------
  ! chm_retrieve_col_boundary
  !--------------------------------------------------------------------------
  function chm_retrieve_col_boundary(headerIndex) result(bound_press)
    !
    !:Purpose: To retrieve previously saved column boundary data in
    !          chm_column_boundary from the header index.
    !
    implicit none
    real(8) :: bound_press

    ! Arguments:
    integer, intent(in) :: headerIndex

    ! Locals:
    character(len=22) :: code

    write(code,'(I22)') headerIndex
    
    bound_press = oss_obsdata_get_element(chm_column_boundary,code,1)

  end function chm_retrieve_col_boundary

  !--------------------------------------------------------------------------
  ! chm_add_sigma_trial
  !--------------------------------------------------------------------------
  subroutine chm_add_sigma_trial(headerIndex,sigma)
    !
    !:Purpose: To add background sigma profiles (and inverse) to chm_sigma_trial
    !          which can be retrieved later using a header index.
    !
    implicit none 

    ! Arguments:
    integer, intent(in) :: headerIndex
    real(8), intent(in) :: sigma(:,:)

    ! Locals:
    integer :: obsdata_maxsize
    
    obsdata_maxsize =  chm_obsdata_maxsize
     
    if (.not.associated(chm_sigma_trial%data2d)) then
       call oss_obsdata_alloc(chm_sigma_trial, obsdata_maxsize, dim1=size(sigma,dim=1), dim2_opt=max(4,size(sigma,dim=2)))
       chm_sigma_trial%nrep = 0
    end if

    ! In this case nrep will count the number of filled reps in the data arrays
    chm_sigma_trial%nrep = chm_sigma_trial%nrep+1 

    if (chm_sigma_trial%nrep.gt.obsdata_maxsize) &
         call utl_abort('chm_sigma_trial: Reached max size of array ' // trim(utl_str(obsdata_maxsize)) )
  
    ! Use the header number as the unique code for this obs data
    write(chm_sigma_trial%code(chm_sigma_trial%nrep),'(I22)') headerIndex

    chm_sigma_trial%data2d(:,1:2,chm_sigma_trial%nrep) = sigma(:,1:2)

    where (sigma(:,1).gt.0.0D0)
       chm_sigma_trial%data2d(:,3,chm_sigma_trial%nrep) = 1.0D0/sigma(:,1)
    elsewhere
       chm_sigma_trial%data2d(:,3,chm_sigma_trial%nrep) = 0.0D0
    end where

    where (sigma(:,2).gt.0.0D0)
       chm_sigma_trial%data2d(:,4,chm_sigma_trial%nrep) = 1.0D0/sigma(:,2)
    elsewhere
       chm_sigma_trial%data2d(:,4,chm_sigma_trial%nrep) = 0.0D0
    end where

  end subroutine chm_add_sigma_trial

  !--------------------------------------------------------------------------
  ! chm_retrieve_sigma_triol
  !--------------------------------------------------------------------------
  function chm_retrieve_sigma_trial(headerIndex) result(sigma)
    !
    !:Purpose: To retrieve previously saved background sigma profiles
    !          chm_sigma_trial from the header index.
    implicit none
    real(8) :: sigma(chm_sigma_trial%dim1,chm_sigma_trial%dim2)

    ! Arguments:
    integer, intent(in) :: headerIndex

    ! Locals:
    character(len=22) :: code

    write(code,'(I22)') headerIndex
    
    sigma = oss_obsdata_get_array2d(chm_sigma_trial,code)

  end function chm_retrieve_sigma_trial

end module obsOperatorsChem_mod
