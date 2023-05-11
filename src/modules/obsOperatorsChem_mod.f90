
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
  use midasMpi_mod
  use utilities_mod
  use varNameList_mod
  use obsSubSpaceData_mod
  use obsFiles_mod
  use codtyp_mod
  use ozoneClim_mod
  use presProfileOperators_mod
  use timeCoord_mod
  use bCovarSetupChem_mod
  
  implicit none
  save
  private

  ! public procedures
  public :: oopc_CHobsoperators, oopc_diagnOnly, oopc_addEfftempObsfile

  !-------------------------------------------------------------------------
  ! Various structures and parameters for the CH family

  ! module structures
  ! -----------------

  type :: struct_oopc_obsoperators  
  
    !  Structure holding work variables for individual observation operator
    !  applications
    !     
    !  Variable               Description
    !  --------               -----------
    !  nobslev                Number of observations in the profile
    !  nmodlev                Number of model levels in the column
    !  varno                  BUFR descriptor element for obs units
    !  constituentId          BUFR code element of local GRIB Table 08046 identifying the constituent
    !                         (similar to BUFR Table 08043)
    !  operatorCategory       CH obs operator category. One of the following.
    !                           'Interp': interpolators
    !                           'Surface': surface point operators
    !                           'Integ': integration operators
    !                           'LayerAvg': layer averaging operators
    !                         Automatically assigned based on obs BUFR element
    !                         and, when relevant, assocciated obsinfo_chm content.
    !  layerIdentified        .true. if a layer (with identified layer boundaries)
    !                         .false if layer boundaries are not available.
    !  valertop               Layer top (final work values in Pa)
    !  vlayerbottom           Layer bottom (final work values in Pa)
    !  zh                     Initial innovation model array (other than conversion constants)
    !  zhp                    Part of innovation operator not related to resolution.
    !  modlevindexTop         Top level of non-zero values in zh
    !  modlevindexBot         Bottom level of non-zero values in zh
    !  trial                  Trial (background) profile on model levels at observation location
    !  tt                     Temperature profile on model levels (Kelvin)
    !  hu                     Specific humidity  on model levels
    !  height                 Height on model levels (m)
    !  pp                     Pressure on model levels (Pa)
    !  lat                    Latitude of observation (radians)
    !  lon                    Longitude of observation (radians)
    !  obslev                 Observation profile level values (OBS_PPP)
    !  obsSpaceTrial          obs estimate on obs levels based on trial/background profile          
    !  varName                Variable/obs nomvar
    !  stnid                  Observation station ID
    !  date                   YYYYMMDD (date of obs)
    !  hhmm                   HHMM (time of obs)
    !  obs_index              Observation index
    !                         Note: Depending on the data of interest, the index of a required array element or 
    !                               profile associated to an observation can be identified from (lat,long,date,hhmm,
    !                               stnid, optional task-dependent identifier if needed) or "obs_index". 
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
    !                         indicates index in oopc_avgkern%obsSubSpace arrays.
    !  applyGenOper           Indicates if the generalized observation operator should be applied
    !  columnBound            Boundary imposed on a column measurement
    !  ixtr                   Indicates when outside model vertical range (when >0)
    !  success                Indicates if the observation is to be assimilated
    !                         or was successfully assimilated 
    
    integer :: nobslev,nmodlev,constituentId,vco,varno,date,hhmm,iavgkern,obs_index
    logical :: layerIdentified,applyGenOper
    real(8) :: lat,lon,columnBound
    character(len=12) :: stnid
    character(len=10) :: operatorCategory
    character(len=4)  :: varName
    real(8), allocatable :: vlayertop(:),vlayerbottom(:),tt(:),height(:),pp(:)
    real(8), allocatable :: zh(:,:),zhp(:,:),obslev(:),hu(:)
    integer, allocatable :: ixtr(:)
    logical, allocatable :: success(:)
    real(8), pointer     :: trial(:)
    integer, allocatable :: modlevindexTop(:),modlevindexBot(:)
    real(8), allocatable :: obsSpaceTrial(:)

  end type struct_oopc_obsoperators
 
  type(struct_oopc_obsoperators) :: obsoper
  
  type :: struct_oopc_operatorsDepot 
  
    !  Structure holding saved arrays for observation operators
    !  A subset of 'struct_oopc_obsoperators'
    !     
    !  Variable               Description
    !  --------               -----------
    !  nobslev                Number of observations in the profile
    !  valertop               Layer top (final work values in Pa)
    !  vlayerbottom           Layer bottom (final work values in Pa)
    !  zh                     Initial innovation model array (other than conversion constants)
    !  zhp                    Part of innovation operator not related to resolution.
    !  modlevindexTop         Top level of non-zero values in zh
    !  modlevindexBot         Bottom level of non-zero values in zh
    !  trial                  Trial (background) profile on model levels at observation location
    !  iavgkern               Integer indicating if averaging kernels are to be applied. Value
    !                         of zero indicates no averaging kernel to be applied. Non-zero value
    !                         indicates index in oopc_avgkern%obsSubSpace arrays.
    !  applyGenOper           Indicates if the generalized observation operator should be applied
    !  ixtr                   Indicates when outside model vertical range (when >0)
    !  success                Indicates if the observation is to be assimilated
    !                         or was successfully assimilated 
    
    integer :: nobslev,iavgkern
    logical :: applyGenOper
    real(8), allocatable :: vlayertop(:),vlayerbottom(:)
    real(8), allocatable :: zh(:,:),zhp(:,:)
    integer, allocatable :: ixtr(:)
    logical, allocatable :: success(:)
    real(8), pointer     :: trial(:)
    integer, allocatable :: modlevindexTop(:),modlevindexBot(:)

  end type struct_oopc_operatorsDepot

  type :: struct_oopc_info
  
    !  Information arrays retrieved from auxiliary file regarding vertical levels 
    !  or averaging kernels
    !
    !  Variable               Description
    !  --------               -----------
    !  n_stnid                Number of sub-families (identified via STNIDs)
    !  stnids                 Sub-families (STNIDs; * are wild cards)
    !  element                BUFR element in data block
    !  profElement            BUFR element of profile required by obs operator if
    !                         if it differs from 'element' (e.g. see 'n_col' and 'type' below).
    !  source                 0: Set entirely from the auxiliary file being read. No 
    !                            initial values read from observation files
    !                         1: Initial values in observation files for constant number
    !                            of vertical levels (may be adjusted after input)
    !  vco                    Index of vertical coord type for obs
    !                           1 - Altitudes (m)
    !                           2 - Pressure (Pa)
    !                           3 - Channel index
    !                           4 - not provided with obs. Obs is for total column values.
    !                           5 - not provided with obs. Obs is a surface point value.
    !                           6 - sigma vertical coordinate (*Psfc to give levels in Pa); in this module only.
    !  ibegin                 Position index of start of data for given
    !                         sub-family.
    !  n_lvl                  Lengh of corresponding obs profile (number of rows)
    !  n_col                  Length of array for each obs element (number of columns; optional) 
    !                         - Usually n_lvl will be the number of vertical levels of profiles 
    !                           and n_col is not needed.
    !                         - For averaging kernels:
    !                            - When the obs are single level data, n_lvl=1 or 2. If it requires more levels
    !                              or data for observation operators, n_col = number of additional data 
    !                            - The a-priori contribution (I-A)xa can be added as the element on each row.
    !                              When done and n_lvl>1, n_col would be expected to be n_lvl+1
    !                            - If an integeration over the n_lvl values (>1) needs to be performed,
    !                              the column n_col = n_lvl+2 provides the integration weights for a summations
    !                              that gives the integral.
    !                         - For vertical levels:
    !                            - n_col will be set to 2 as default for layer boundaries
    !                            - n_col needs to be set to 1 to use this simply as a single vertical level profile
    !  n_lat                  Number of latitudes
    !  lat                    Latitudes (degrees; ordered in increasing size)
    !
    !  vlayertop              Layer top 
    !  vlayerbottom           Layer bottom
    !  rak                    Averaging kernel matrices
    !  type                   Operation subtypes.
    !
    !                         Currently relevant only for averaging kernels:
    !
    !                           Type name          Description
    !                           ============================================================
    !                           'default'          No special treatment (default)
    !                           'log'              Application of log-space averaging kernel
    !                           ============================================================
    !
    integer ::  n_stnid
    character(len=12), allocatable :: stnids(:)
    character(len=15), allocatable :: type(:)  
    integer, allocatable :: element(:),source(:),profElement(:)
    integer, allocatable :: vco(:),n_lat(:)
    integer, allocatable :: ibegin(:),n_lvl(:),n_col(:)
    real(8), allocatable :: rak(:),vlayertop(:),vlayerbottom(:)
    real(8), allocatable :: lat(:)
 
    type(struct_oss_obsdata), allocatable :: obsSubSpace(:)
 
  end type struct_oopc_info

  type(struct_oopc_info) :: oopc_levels
  type(struct_oopc_info) :: oopc_avgkern

  ! Arrays for integration upper boundary of retrieved total column measurements 
  type(struct_oss_obsdata) :: oopc_columnBoundary

  ! File name of auxiliary text file constaining supplemental observation information
  character(len=50), parameter :: oopc_aux_filename="obsinfo_chm" 

  ! Max nummber of constituents (max size of related arrays)
  integer, parameter :: oopc_constituentsSize=30  ! = max allowed value of "iconstituentId" for Table 08046.
                                                  ! Value to be increased as needed up to a max of 6999 as values
                                                  ! > 7000 (and less 0) are assumed assigned to non-constituent fields  

  ! Arrays containing input reference fields and fields interpolated 
  ! to obs locations
  
  type :: struct_oopc_field

    !  Structure for storing reference (climatological) fields needed for
    !  operatorSubType(2,:) == 'genOper' with 
    !  genOperConstraintType == 'Diff' (see below).
    !     
    !  Variable               Description
    !  --------               -----------
    !  field                  Gridded 3D field (lon,lat,vlev) or 2D field (1,lat,vlev)
    !  nlat                   number of latitudes
    !  nlon                   number of longitudes
    !  nlev                   number of vertical levels
    !  lat,lon                lat,lon grid in radians
    !  vlev                   vertical levels
    !  ivkind                 Index of vertical coordinate type. Defintion may vary according to source.
    !                         For fields read from RPN files and use of convip:
    !                             0: P is in height [m] (metres) with respect to sea level 
    !                             1: P is in stddev [sg] (0.0 -> 1.0) 
    !                             2: P is in pressure [mb] (millibars) 
    !                             3: P is in an arbitrary code 
    !                             4: P is in height [M] (metres) with respect to ground level 
    !                             5: P is in hybrid coordinates [hy] 
    !                             6: P is in theta [th] 
    !                         For use with obs                      
    
    real(8), allocatable :: field(:,:,:),lat(:),lon(:),vlev(:)
    integer :: nlev,nlon,nlat,ivkind
  
  end type struct_oopc_field

  type(struct_oss_obsdata)  :: oopc_bgRef
  type(struct_oopc_field)   :: oopc_climFields(0:oopc_constituentsSize,2)

  ! Arrays to contain the calculated concentration-weighted effective temperature
  ! associated to total column data. It will be stored in the observation file.

  type(struct_oss_obsdata) :: oopc_efftemp

  type(struct_bcsc_bgStats) :: bgStats ! Background covariances

  ! Setup initialization key
  logical :: initializedChem = .false.    
  
  ! Operator storage initialization key
  logical :: initializedOperators = .false.    
 
  ! General config/setup information parameters 
  ! See description list of NAMCHEM namelist parameters and others 
  ! in routine oopc_readNamchem.
  ! Following variables/parameters could be placed in a data structure/type
  ! (e.g. struct_oopc_nmlparm)

  character(len=5) :: oopc_genOperConstraintType(0:oopc_constituentsSize) 
  real(8) :: oopc_genOperHCorrlenExpnt(0:oopc_constituentsSize)  
  real(8) :: oopc_genOperOmAStatsFactor(0:oopc_constituentsSize) 
  integer :: oopc_tropo_mode(0:oopc_constituentsSize),oopc_tropo_bound(0:oopc_constituentsSize)
  integer :: oopc_obsdata_maxsize
  real(8) :: oopc_tropo_column_top(0:oopc_constituentsSize)
  logical :: oopc_storeOperators
  
  integer, parameter :: assim_maxfamnum=1           ! Could be used for other families as well with >1
  integer, parameter :: assim_maxsize=100           ! max size of assim_* arrays  

  integer :: assim_famNum
  character(len=2) :: assim_fam(assim_maxfamnum) ! List of families to which filt_diagnOnly is to apply
  logical :: assim_all(assim_maxfamnum) ! Choose to assimilate all obs of specified family
  integer :: assim_num(assim_maxfamnum) ! Number of combinations identified for assimilation
  integer :: assim_varno(assim_maxfamnum,assim_maxsize) ! List of bufr elements to assimilate (0 means all)
  integer :: assim_nlev(assim_maxfamnum,assim_maxsize) ! 0: multi- and uni-lev; 1: uni-lev; >1 multi-lev
  integer :: assim_exclude_nflag(assim_maxfamnum) ! List of bits for excluding obs from assimilation
  integer :: assim_exclude_flag(assim_maxfamnum,assim_maxsize) ! Number of bits for excluding obs
  character(len=9) :: assim_stnid(assim_maxfamnum,assim_maxsize) ! List of stnid to assimilation '*' for wild card
  character(len=20) :: operatorSubType(2,assim_maxsize) ! Operator sub-type name
  character(len=10) :: modelName = 'GEM-MACH' ! Identification of the model

  !**************************************************************************
  
  contains

  !--------------------------------------------------------------------------
  ! oopc_setupCH
  !--------------------------------------------------------------------------
  subroutine oopc_setupCH(kmode)
    !
    !:Purpose: To set up additional information required by constituent obs and
    !          not provided in obsSpaceData.  Also to assign observation layer
    !          top and bottom levels (and averaging kernel matrices).
    !          See 'oopc_CHobsoperators'. 
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: kmode ! Mode of observation operator

    ! Locals:
    logical :: success
    
    write(*,*) 'Begin oopc_setupCH'

    ! Read NAMCHEM namelist and set related parameters

    call oopc_readNamchem
    write(*,*) 'oopc_setupCH: Completed oopc_readNamchem'
      
    ! Read reference vertical levels (where needed) or 
    ! top and bottom layer boundaries of partial (or total) column meausurements
  
    call oopc_readLevels
    write(*,*) 'oopc_setupCH: Completed oopc_readLevels'
      
    ! To deallocate space if required elsewhere, one should use
    ! call oopc_deallocLevels
   
    ! Read averaging kernel matrices
  
    call oopc_readAvgkern
    write(*,*) 'oopc_setupCH: Completed oopc_readAvgkern'
  
    ! To deallocate space if required elsewhere, one should use
    ! call oopc_deallocAvgkern
  
    ! Read reference (e.g. climatological) fields
  
    call oopc_readFields(oopc_climFields,oopc_aux_filename,'CH', &
                         oopc_constituentsSize,2,oopc_genOperConstraintType, &
			 success,filetype_opt='TXT')
    if ( .not. success ) then
      call utl_abort('oopc_setupCH: Failed in oopc_readFields')
    end if
    write(*,*) 'oopc_setupCH: Completed oopc_readFields'
    
    ! Allocation of oopc_efftemp done in oopc_setupCH instead of obsdata_add_data1d
    ! to ensure allocation is done for all processors, including those without associated data.
    ! This is to ensure that rpn_comm_allgather will work in routine obsdata_MPIGather.

    if (.not.associated(oopc_efftemp%data1d)) then
      call oss_obsdata_alloc(oopc_efftemp,oopc_obsdata_maxsize,dim1=1)
      oopc_efftemp%nrep=0
    end if

    ! Get background error std dev stats when they may be required
    
    bgStats%initialized = .false.
    
    if (kmode == 1 .or. any(operatorSubType(2,:) == 'genOper') ) then
      call bcsc_getCovarCH(bgStats)
    end if

    write(*,*) 'Completed oopc_setupCH'
  
  end subroutine oopc_setupCH

  !--------------------------------------------------------------------------
  ! oopc_readNamchem
  !--------------------------------------------------------------------------
  subroutine oopc_readNamchem
    !:Purpose: Read and store miscellaneous flags and constants.
    !
    !:Comment: assim_* arrays could instead be made available to all families
    !          by moving them to a different input namelist (and changing its
    !          dimensions settings).
    !
    ! 
    !     :genOperConstraintType:
    !                           Reference profile type for weighted integration
    !                           or layer averaging (generalized observation) 
    !                           operator.
    !                           Relevant for operatorSubType(2,i)='genOper'.
    !                           ================================================
    !                           'Trial'  use trial field xb for mass weighted
    !                                    increment distribution
    !                           'Diff'   use a combination of the difference of
    !                                    an external reference xc and the trial
    !                                    field xb, i.e. mass weighted increment
    !                                    distribution as a(xc-xb) + b*xc where a
    !                                    and b depend on the size of
    !                                    sum[(xc-xb)/sig(xb)]^2 over the profile
    !                             ==============================================
    !
    !     :genOperHCorrlenExpnt: Used with operatorSubType(2,i) ='genOper'
    !                           Exponent for partially mitigating the effect of 
    !                           the influence of neighbouring column amonunt obs
    !                           from background error correlations.
    !                           Emperically obtained exponent value.
    !                           Not optimal for all possible local horizontal data densities.
    !
    !     :genOperOmAStatsFactor:  OmA RMS (or std dev) conservation factor for 
    !                           operatorSubType(2,i) ='genOper'.
    !
    !
    !     :assim_fam:           List of families to which filt_diagnOnly is to
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
    !                           When assim_all is .false., account for the setttings 
    !                           of assim_num, assim_varno, assim_stnid, assim_nlev.
    !
    !     :assim_num:           Relevant when assim_all = ,false.
    !                           Number combinations (stnid, bufr element,
    !                           multi/uni-level) identified for assimilation.
    !                           All others will not be assimilated. OmP and OmA
    !                           diagnostics and output will still be produced
    !                           for non-assimilated datasets.
    ! 
    !                             ===  =======================================
    !                              0   none are to be assimilated when assim_all
    !                                  is .false. (default)
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
    !                                  obsoper%columnBound set to zero. If
    !                                  specified, generalized innovation
    !                                  operator only applied below
    !                                  obsoper%columnBound in the tangent
    !                                  linear model.
    !                              2   Values of tangent linear model and
    !                                  adjoint model above obsoper%columnBound
    !                                  set to zero.
    !                             ===  =======================================
    !                           Array index refers to BUFR code element of Table
    !                           08046 (iconstituentId) identifying the
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
    !                           of Table 08046 (iconstituentId) identifying
    !                           the constituent. Relevant for total column
    !                           measurements only.
    ! 
    !     :obsdata_maxsize:     Max allowed size of work arrays (in terms of
    !                           number of obs) associated to ordered observation
    !                           indices
    ! 
    ! 
    !     :modelName:           Identifier of forecast model
    !                           Default: 'GEM-MACH'
    !                           Set to 'GEM' for varNames of 'O3L', 'CH4L', and 'N2OL'
    !
    !     :operatorSubType(2,i):Operator sub-type name.
    !                           Index (2,i) for sub-type to apply for stnid in element (1,i)
    !                           See related "obsoper@operatorCategory" automatically assigned
    !                           based on obs BUFR element and obsinfo_chm content.
    !
    !                           Operator        Sub-type name    Description
    !                           Category 
    !                           =============  =====================================================================
    !                           'Interp'        'default'        Piecewise linear interpolation (default)
    !                                           'wgtAvg'         Piecewise weighted averaging  interpolator
    !                           'Surface'       'default'        No special treatment (default)
    !                           'Integ'         'default'        Simple/basic vertical integration (default)
    !                                           'genOper'        Weighted vertical integration - see 'genOper*' parameters
    !                           'LayerAvg'      'default'        Simple layer averaging (default)
    !                                           'genOper'        Weighted vertical layer averaging - see 'genOper*' parameters
    !                           ====================================================================================
    !
    !                           Notes:
    !
    !                           - 'genOper' requires NAMBCHM namelist parameter settings
    !                             getPhysSpaceStats=.true. and 
    !                             getPhysSpaceHCorrel=.true.
    !                           - Application of averaging kernels is directed only
    !                             by the content of the obsinfo_chm file 'SECTION III'
    !                                             
    !     :storeOperators:      Logical indicating if linear operators are stored for re-use in TL and AD calc.
    !                           If so, the linear operators will not be re-calculated at different iterations.
    !                           Not used when tropo_mode>=1
    !
    
    implicit none

    ! Locals:
    integer :: fnom, fclos
    integer :: ierr, ios, nulnam, i
    character(len=10)  :: namfile 

    ! Namelist variables (local)
    integer :: tropo_mode(0:oopc_constituentsSize) ! Special treatment for troposphere of total column obs
    integer :: tropo_bound(0:oopc_constituentsSize) ! Indicate which column top value used for special treatment
    real(8) :: tropo_column_top(0:oopc_constituentsSize) ! Default for column boundary (in Pa) of total column obs
    logical :: storeOperators ! Choose to store linear operators for re-use in TL/AD
    character(len=5) :: genOperConstraintType(0:oopc_constituentsSize) ! Strong constraint for generalized obs operator (see oopc_genOper)
    real(8) :: genOperHCorrlenExpnt(0:oopc_constituentsSize)  ! Exponent for horiz. correl. length weighting in oopc_genOper
    real(8) :: genOperOmAStatsFactor(0:oopc_constituentsSize) ! Additional OmAStats normalization factor for oopc_genOper
    integer :: obsdata_maxsize ! Max number of obs associated with ordered obs indices

    external fnom,fclos

    namelist /namchem/ assim_fam,assim_all,assim_num,assim_stnid,assim_varno,    &
                       assim_nlev, assim_exclude_nflag,assim_exclude_flag,       &
                       tropo_mode,tropo_bound,tropo_column_top,obsdata_maxsize,  &
                       modelName,operatorSubType,storeOperators, &
		       genOperConstraintType,genOperHCorrlenExpnt,genOperOmAStatsFactor
  
    ! Default NAMCHEM values

    genOperConstraintType(:)='Trial'
    genOperHCorrlenExpnt(:)=1.7
    genOperOmAStatsFactor(:)=1.0
    
    obsdata_maxsize=90000
    
    storeOperators=.false.
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
    
    operatorSubType(1,:) = ''
    operatorSubType(2,:) = 'default'

    ! Read from namelist file NAMCHEM

    namfile=trim("flnml")
    nulnam=0
    ierr=fnom(nulnam,namfile,'R/O',0)

    read(nulnam,nml=namchem,iostat=ios)
    if (ios < -4 .or. ios > 0) then 
      call utl_abort('oopc_readNamchem: Error in reading NAMCHEM namelist. iostat = ' // trim(utl_str(ios)) )   
    else if (mmpi_myid == 0) then
      write(*,nml=namchem)   
    end if
  
    ierr=fclos(nulnam)      

    do i=size(assim_fam),1,-1
      if (assim_fam(i) /= '') exit
    end do
    assim_famnum=i

    oopc_genOperConstraintType(:) = genOperConstraintType(:)
    oopc_genOperHCorrlenExpnt(:) =  genOperHCorrlenExpnt(:)
    oopc_genOperOmAStatsFactor(:) = genOperOmAStatsFactor(:)
    oopc_tropo_mode(:) = tropo_mode(:)
    oopc_tropo_bound(:) = tropo_bound(:)
    oopc_tropo_column_top(:) = tropo_column_top(:)
    oopc_storeOperators = storeOperators
    oopc_obsdata_maxsize=obsdata_maxsize
  
  end subroutine oopc_readNamchem

  !--------------------------------------------------------------------------  
  !------- Routines related to reference or layer top & bottom levels -------

  !--------------------------------------------------------------------------
  ! oopc_readLevels
  !--------------------------------------------------------------------------
  subroutine oopc_readLevels
    !
    !:Purpose: To read and to store reference levels (where needed) or top- and 
    !          bottom-layer boundaries for CH sub-families.
    !
    implicit none

    ! Locals:
    integer :: fnom, fclos
    integer :: ierr, jlev, jelm, nulstat, ios, isize, icount
    logical :: LnewExists,newread
  
    character (len=128) :: ligne

    external fnom,fclos
  
    ! Initialization

    oopc_levels%n_stnid=0

    inquire(file=trim(oopc_aux_filename),EXIST=LnewExists)
    if (.not.LnewExists )then
      write(*,*)   '---------------------------------------------------------------'
      write(*,*)   'oopc_readLevels: COULD NOT FIND AUXILIARY file ' // trim(oopc_aux_filename)
      write(*,*)   '---------------------------------------------------------------'
      return
    end if

    ! Check for available layer info.

    nulstat=0
    ierr=fnom(nulstat,trim(oopc_aux_filename),'SEQ',0)
    if ( ierr == 0 ) then
      open(unit=nulstat, file=trim(oopc_aux_filename), status='OLD')
    else
      call utl_abort('oopc_readLevels: COULD NOT OPEN AUXILIARY file ' // trim(oopc_aux_filename))
    end if

    ios=0
    read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne
    do while (trim(adjustl(ligne(1:13))) /= 'SECTION II:') 
      read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne
    end do    
  
    ! Read number of observation set sub-families (STNIDs and ...) and allocate space
   
    read(nulstat,*,iostat=ios,err=10,end=10) oopc_levels%n_stnid
    read(nulstat,*,iostat=ios,err=10,end=10) isize

    allocate(oopc_levels%stnids(oopc_levels%n_stnid))
    allocate(oopc_levels%vco(oopc_levels%n_stnid))
    allocate(oopc_levels%source(oopc_levels%n_stnid),oopc_levels%ibegin(oopc_levels%n_stnid))
    allocate(oopc_levels%element(oopc_levels%n_stnid),oopc_levels%n_lvl(oopc_levels%n_stnid))
    allocate(oopc_levels%vlayertop(isize),oopc_levels%vlayerbottom(isize))
    allocate(oopc_levels%n_col(oopc_levels%n_stnid))
 
    oopc_levels%element(:)=0
    oopc_levels%vco(:)=0
    oopc_levels%source(:)=0
    oopc_levels%n_lvl(:)=1
    oopc_levels%n_col(:)=2

    ! Begin reading for each sub-family
    ! Important: Combination of STNID, BUFR element and number of vertical levels
    !            to determine association to the observations.

    icount=0
    do jelm=1,oopc_levels%n_stnid
      oopc_levels%ibegin(jelm)=icount+1

      read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne

      ! Read STNID (* is a wildcard)
    
      read(nulstat,'(2X,A9)',iostat=ios,err=10,end=10) oopc_levels%stnids(jelm) 

      ! Read (1) Obs BUFR element.
      !      (2) Vertical coord type (1, 2, or 3)
      !      (3) Flag indication if EOR provided from this auxiliary file or
      !          to be read from the observation file,
      !      (4) Number of vertical levels (rows)
      !      (5) Number of vertical level profiles/columns (optional; default is 2)
      !
      ! Important: Combination of STNID, BUFR element and number of vertical levels  
      !            to determine association to the observations.

      newread=.true.
      read(nulstat,*,iostat=ios,err=20,end=20) oopc_levels%element(jelm),oopc_levels%vco(jelm),  &
        oopc_levels%source(jelm),oopc_levels%n_lvl(jelm),oopc_levels%n_col(jelm) 
      newread=.false.
 20   if (newread) then
        backspace(nulstat,iostat=ios,err=10)
        backspace(nulstat,iostat=ios,err=10)       
        read(nulstat,*,iostat=ios,err=10,end=10) oopc_levels%element(jelm),oopc_levels%vco(jelm),  &
          oopc_levels%source(jelm),oopc_levels%n_lvl(jelm) 
      end if

      if (icount+oopc_levels%n_lvl(jelm) > isize) then
        call utl_abort('oopc_readLevels: READING PROBLEM. ' // &
	                'Max array size exceeded: ' // &
                        trim(utl_str(isize)))
      end if    

      read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne
          
      if ( oopc_levels%n_lvl(jelm) > 0 ) then
        if (oopc_levels%n_col(jelm) > 1) then
      
          do jlev=1,oopc_levels%n_lvl(jelm)
            icount=icount+1
          
            ! Read top and bottom levels
           
            read(nulstat,*,iostat=ios,err=10,end=10)                 &
                   oopc_levels%vlayertop(icount),oopc_levels%vlayerbottom(icount)
          end do
          
        else
        
          do jlev=1,oopc_levels%n_lvl(jelm)
            icount=icount+1
          
            ! Read single level
           
            read(nulstat,*,iostat=ios,err=10,end=10)                 &
                   oopc_levels%vlayertop(icount)
            
          end do
          oopc_levels%vlayerbottom(:)= oopc_levels%vlayertop(:)

        end if              
      end if

    end do
   
 10 if (ios > 0) then
      call utl_abort('oopc_readLevels: READING PROBLEM. ' // &
           'File read error message number: ' // trim(utl_str(ios)))
    end if    
 
 11 CLOSE(UNIT=nulstat)
    ierr=fclos(nulstat)    

  end subroutine oopc_readLevels

  !--------------------------------------------------------------------------
  ! oopc_getLevels
  !--------------------------------------------------------------------------
  subroutine oopc_getLevels
    !
    !:Purpose: To return reference model levels or layer boundaries for an observation. 
    !          Combination of STNID, element variable number and number of vertical 
    !          levels to determine association to the observations. Default values
    !          for top and bottom layers for total column measurements are to be
    !          provided.
    !
    
    implicit none

    ! Arguments:

    ! Locals:
    integer :: nlev              ! number of levels in the observation
    integer :: istnid,stnidIndex,startIndex,levelIndex,level
    logical :: iset
    real(8), allocatable :: levelsTop(:),levelsBot(:)
    
    nlev=obsoper%nobslev

    ! Find stnid with same number of vertical levels, and same BUFR element
          
    istnid=0
    obsoper%layerIdentified=.false.

    do stnidIndex=1,oopc_levels%n_stnid

       ! First compare STNID values allowing for * and blanks in 
       ! oopc_levels%stnids(JN) as wildcards
       iset = utl_stnid_equal(oopc_levels%stnids(stnidIndex),obsoper%stnid)

       ! Check if number of levels, code, and vertical coordinate type are equal.
       ! If number of levels is one and no vertical coordinate provided for total column measurement (i.e. obsoper%vco == 4),
       ! then check of vertical coordinate type is disregarded afterwards.       
       if (iset) then	      
          if ( obsoper%varno == oopc_levels%element(stnidIndex) .and. &
             (nlev == oopc_levels%n_lvl(stnidIndex) .or. &
	     (nlev == 1 .and. oopc_levels%n_lvl(stnidIndex) > 1) ) .and. &
             (obsoper%vco == oopc_levels%vco(stnidIndex) .or. &
	      obsoper%vco == 4) ) then
	    
            istnid=stnidIndex
            exit
	      
          end if
       end if
       
    end do
    
    if (istnid == 0) then
       ! If integrated layer information not found, if a total column measurement
       ! set to defaults, else do nothing

       if (bufr_IsIntegral(obsoper%varno) .and. nlev == 1) then          
          obsoper%layerIdentified=.true.
          obsoper%vlayertop(1) = obsoper%pp(1)
          obsoper%vlayerbottom(1) = obsoper%pp(obsoper%nmodlev)
       end if

    else
       ! level or layer information has been found in auxiliary file
       
       if (nlev == oopc_levels%n_lvl(stnidIndex) .and. nlev > 1) then
         ! layer information has been found in auxiliary file
         obsoper%layerIdentified=.true.
         startIndex = oopc_levels%ibegin(stnidIndex)
         obsoper%vlayertop(:) = oopc_levels%vlayertop(startIndex:startIndex+nlev-1)
         obsoper%vlayerbottom(:) = oopc_levels%vlayerbottom(startIndex:startIndex+nlev-1) 
       else 
                
         ! Alternative level or layer information has been found
         ! Must be pressure or sigma
         
         if (oopc_levels%vco(stnidIndex) /= 2 &
	   .and. oopc_levels%vco(stnidIndex) /= 6) then
	   
           call utl_abort('oopc_getLevels: Cannot handle this vertical coordinate type')
	 end if 
          
         obsoper%nobslev=oopc_levels%n_lvl(stnidIndex)
         
         if (oopc_levels%n_col(stnidIndex) > 1) obsoper%layerIdentified=.true.
        
         startIndex = oopc_levels%ibegin(stnidIndex)

         allocate(levelsTop(obsoper%nobslev))   
         allocate(levelsBot(obsoper%nobslev))
         levelsTop(:) = oopc_levels%vlayertop(startIndex:startIndex+obsoper%nobslev-1)
         levelsBot(:) = oopc_levels%vlayerbottom(startIndex:startIndex+obsoper%nobslev-1)
         if (oopc_levels%vco(stnidIndex) == 6) then
           levelsTop(:) = levelsTop(:)*obsoper%pp(obsoper%nmodlev)
           levelsBot(:) = levelsBot(:)*obsoper%pp(obsoper%nmodlev)
         end if
          
         ! Ensure levels do not go below the surface

         do levelIndex=1,obsoper%nobslev
           if (levelsBot(levelIndex) >= obsoper%pp(obsoper%nmodlev)) then  
             obsoper%nobslev = levelIndex
             levelsBot(levelIndex) = obsoper%pp(obsoper%nmodlev)
             if (levelsTop(levelIndex) >= obsoper%pp(obsoper%nmodlev)) then
               levelsTop(levelIndex) = obsoper%pp(obsoper%nmodlev)
	     end if
             exit
           end if
         end do
           
         ! Reset work array sizes and settings
         deallocate(obsoper%vlayertop,obsoper%vlayerbottom,obsoper%obslev)
         deallocate(obsoper%modlevindexTop,obsoper%modlevindexBot,obsoper%zh,obsoper%zhp)

         allocate(obsoper%obslev(obsoper%nobslev))       ! Reference vertical levels
         allocate(obsoper%vlayertop(obsoper%nobslev))    ! Layer tops for layer measurements
         allocate(obsoper%vlayerbottom(obsoper%nobslev)) ! Layer bottoms for layer measurements
         allocate(obsoper%modlevindexTop(obsoper%nobslev))  ! Index of highest model level (lowest index) involved with obs element
         allocate(obsoper%modlevindexBot(obsoper%nobslev))  ! Index of lowest model level (highest index) involved with obs element
         allocate(obsoper%zh(obsoper%nobslev,obsoper%nmodlev))   ! Local model operator H (excluding conversion constants and horizontal interpolation)
         allocate(obsoper%zhp(obsoper%nobslev,obsoper%nmodlev))  ! Part of zh that excludes aspects related to vertical resolition

         obsoper%vlayertop(:) = levelsTop(1:obsoper%nobslev)
         obsoper%vlayerbottom(:) = levelsBot(1:obsoper%nobslev)
         obsoper%obslev(:) = levelsTop(1:obsoper%nobslev)
         obsoper%zh(:,:)=0.0d0
         obsoper%zhp(:,:)=0.0d0

         deallocate(levelsTop,levelsBot)

       end if

       obsoper%modlevindexTop(:)=1
       startIndex=2
       do level=1,obsoper%nobslev
         if ( obsoper%vlayertop(level) < obsoper%pp(2)) then
	   obsoper%modlevindexTop(level)=1
	 else 
           do levelIndex=startIndex,obsoper%nmodlev
             if (obsoper%vlayertop(level) == obsoper%pp(levelIndex)) then
                obsoper%modlevindexTop(level)=levelIndex
		exit
             else if (obsoper%vlayertop(level) < obsoper%pp(levelIndex)) then
                obsoper%modlevindexTop(level)=levelIndex-1
                exit
             end if
           end do
           if (level > 1) then
             if (obsoper%modlevindexTop(level) <  obsoper%modlevindexTop(level-1)) &
               obsoper%modlevindexTop(level)=obsoper%modlevindexTop(level-1)
           end if
           startIndex=levelIndex
         end if
       end do
       
                                   
       obsoper%modlevindexBot(:)=obsoper%modlevindexTop(:)+1
       startIndex=obsoper%nmodlev-1
       do level=obsoper%nobslev,1,-1
         if ( obsoper%vlayertop(level) > obsoper%pp(obsoper%nmodlev-1)) then
	   obsoper%modlevindexTop(level)=obsoper%nmodlev
	 else 
           do levelIndex=startIndex,1,-1
             if (obsoper%vlayerbottom(level) == obsoper%pp(levelIndex)) then
                obsoper%modlevindexBot(level)=levelIndex
	        exit
             else if (obsoper%vlayerbottom(level) > obsoper%pp(levelIndex)) then
                obsoper%modlevindexBot(level)=levelIndex+1
                exit
             end if
           end do
           if (level < obsoper%nobslev) then
             if (obsoper%modlevindexBot(level) >  obsoper%modlevindexBot(level+1)) &
               obsoper%modlevindexBot(level)=obsoper%modlevindexBot(level+1)
           end if
	   startIndex=levelIndex
	 end if
       end do
       
    end if

  end subroutine oopc_getLevels

  !--------------------------------------------------------------------------
  ! oopc_deallocLevels
  !--------------------------------------------------------------------------
  subroutine oopc_deallocLevels
    !
    !:Purpose: To deallocate temporary storage space used for layer info
    !

    implicit none

    if (oopc_levels%n_stnid == 0) return

    call oopc_deallocInfo(oopc_levels)
 
  end subroutine oopc_deallocLevels

  !--------------------------------------------------------------------------
  !-------------- Routines related to averaging kernel matrices -------------

  !--------------------------------------------------------------------------
  ! oopc_readAvgkern
  !--------------------------------------------------------------------------
  subroutine oopc_readAvgkern
    !
    !:Purpose: To read averaging kernels from auxiliary file or observation file
    !              
    implicit none

    ! Locals:
    integer, parameter :: ndim=2

    integer :: istnid

    ! read the averaging kernel information from the auxiliary file
    call oopc_readAvgkernAuxfile

    ! set size of observation file array
    allocate(oopc_avgkern%obsSubSpace(oopc_avgkern%n_stnid))

    ! read from observation file
    do istnid=1,oopc_avgkern%n_stnid
      if (oopc_avgkern%source(istnid) == 1) then
        
        ! retrieve data from stats blocks (with bkstp=14 and block_type='DATA')
        oopc_avgkern%obsSubSpace(istnid) = obsf_obsSub_read('CH', &
	  oopc_avgkern%stnids(istnid),bufr_avgkern,oopc_avgkern%n_lvl(istnid), &
          ndim, numColumns_opt=oopc_avgkern%n_col(istnid), &
	  bkstp_opt=14, block_opt='DATA', match_nlev_opt=.true.)
     
      end if
    end do

  end subroutine oopc_readAvgkern

  !--------------------------------------------------------------------------
  ! oopc_readAvgkernAuxfile
  !--------------------------------------------------------------------------
  subroutine oopc_readAvgkernAuxfile
    !
    !:Purpose: To read and to store averaging kernel matricesfor CH sub-families
    !
    !:Comments:
    !      - Currently implemented for only one latitude band
    !

    implicit none

    !Locals:
    integer :: fnom, fclos
    integer :: ierr, jlev, jelm, nulstat, ios, isize, icount, iend
    logical :: LnewExists,newread
  
    character (len=128) :: ligne

    external fnom,fclos

    ! Initialization

    oopc_avgkern%n_stnid=0

    inquire(file=trim(oopc_aux_filename),EXIST=LnewExists)
    if (.not.LnewExists )then
      write(*,*)   '-------------------------------------------------------'
      write(*,*)   'oopc_readAvgkernAuxfile: COULD NOT FIND AUXILIARY file ' // trim(oopc_aux_filename)
      write(*,*)   '-------------------------------------------------------'
      return
    end if

    ! Check for available layer info.

    nulstat=0
    ierr=fnom(nulstat,trim(oopc_aux_filename),'SEQ',0)
    if ( ierr == 0 ) then
      open(unit=nulstat, file=trim(oopc_aux_filename), status='OLD')
    else
      call utl_abort('oopc_readAvgkernAuxfile: COULD NOT OPEN AUXILIARY file ' // trim(oopc_aux_filename))
    end if

    ios=0
    read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne
    do while (trim(adjustl(ligne(1:14))) /= 'SECTION III:') 
      read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne
    end do    
  
    ! Read number of observation set sub-families (STNIDs and ...) and allocate space
   
    read(nulstat,*,iostat=ios,err=10,end=10) oopc_avgkern%n_stnid
    read(nulstat,*,iostat=ios,err=10,end=10) isize

    allocate(oopc_avgkern%stnids(oopc_avgkern%n_stnid))
    allocate(oopc_avgkern%source(oopc_avgkern%n_stnid))
    allocate(oopc_avgkern%ibegin(oopc_avgkern%n_stnid))
    allocate(oopc_avgkern%element(oopc_avgkern%n_stnid))
    allocate(oopc_avgkern%n_lvl(oopc_avgkern%n_stnid))
    allocate(oopc_avgkern%n_col(oopc_avgkern%n_stnid))
    allocate(oopc_avgkern%type(oopc_avgkern%n_stnid))
    allocate(oopc_avgkern%profElement(oopc_avgkern%n_stnid))
    allocate(oopc_avgkern%rak(isize))
 
    oopc_avgkern%element(:)=0
    oopc_avgkern%source(:)=0
    oopc_avgkern%n_lvl(:)=1
    oopc_avgkern%n_col(:)=1
    oopc_avgkern%type(:) = 'default'

    ! Begin reading for each sub-family
    ! Important: Combination of STNID, BUFR element and number of vertical levels 
    !            to determine association to the observations.

    icount=1
    STNIDLOOP: do jelm=1,oopc_avgkern%n_stnid
      oopc_avgkern%ibegin(jelm)=icount

      read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne

      ! Read STNID (* is a wildcard)
    
      read(nulstat,'(2X,A9,1X,A15)',iostat=ios,err=10,end=10) oopc_avgkern%stnids(jelm),oopc_avgkern%type(jelm)
      if (trim(oopc_avgkern%type(jelm)) /= 'log') then
        oopc_avgkern%type(jelm) = 'default'
      end if

      ! Read (1) Obs BUFR element.
      !      (2) Flag indication if avgkern provided from this auxiliary file or
      !          to be read from an observation file,
      !      (3) Number of obs in profile (number of rows)
      !      (4) Number of columns (optional read; number of obs operator vertical levels with or without the
      !          resultant a priori contribution (I-A)xa).
      !          If number of rows (n_lvl) equals 1, the actual levels for the columns need to be specified from 
      !          section II (see routine oopc_readLevels)
      !
      ! Important: Combination of STNID, BUFR element and number of vertical levels
      !            to determine association to the observations.

      newread=.true.
      read(nulstat,*,iostat=ios,err=20,end=20) oopc_avgkern%element(jelm),  &
        oopc_avgkern%source(jelm),oopc_avgkern%n_lvl(jelm),oopc_avgkern%n_col(jelm),oopc_avgkern%profElement(jelm)
      newread=.false.
 20   if (newread) then
        backspace(nulstat,iostat=ios,err=10)
        backspace(nulstat,iostat=ios,err=10)
        read(nulstat,*,iostat=ios,err=10,end=10) oopc_avgkern%element(jelm),  &
          oopc_avgkern%source(jelm),oopc_avgkern%n_lvl(jelm)
        oopc_avgkern%n_col(jelm)=oopc_avgkern%n_lvl(jelm) 
        oopc_avgkern%profElement(jelm)= oopc_avgkern%element(jelm)
      end if
      
      if (icount+oopc_avgkern%n_lvl(jelm) > isize) then
        call utl_abort('oopc_readAvgkernAuxfile: READING PROBLEM. ' // &
	               'Max array size exceeded:' // trim(utl_str(isize)))    
      end if
      read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne

      ! disregard data section if values to be specified in BUFR file
      if (oopc_avgkern%source(jelm) == 1) cycle STNIDLOOP
    
      if (oopc_avgkern%n_lvl(jelm) > 0 .and. oopc_avgkern%n_col(jelm) > 1) then   
        do jlev=1,oopc_avgkern%n_lvl(jelm)

          iend=icount+oopc_avgkern%n_col(jelm)-1
       
          ! Read averaging kernel matrix   
          read(nulstat,*,iostat=ios,err=10,end=10) oopc_avgkern%rak(icount:iend)

          icount=iend+1

        end do
      end if

    end do STNIDLOOP
   
 10 if (ios > 0) then
      call utl_abort('oopc_readAvgkernAuxfile: READING PROBLEM. ' // &
                     'File read error message number: ' // trim(utl_str(ios)))
    end if   
 
 11 CLOSE(UNIT=nulstat)
    ierr=fclos(nulstat)    

  end subroutine oopc_readAvgkernAuxfile

  !--------------------------------------------------------------------------
  ! oopc_deallocAvgkern
  !--------------------------------------------------------------------------
  subroutine oopc_deallocAvgkern
    !
    !:Purpose: To deallocate temporary storage space used for averaging kernels
    !
    implicit none

    ! Locals:
    integer :: istnid

    if (oopc_avgkern%n_stnid == 0) return

    if (allocated(oopc_avgkern%obsSubSpace)) then
      do istnid=1,oopc_avgkern%n_stnid
        if (oopc_avgkern%source(istnid) == 1) then
	  call oss_obsdata_dealloc(oopc_avgkern%obsSubSpace(istnid))
        end if
      end do
      deallocate(oopc_avgkern%obsSubSpace)
    end if

    call oopc_deallocInfo(oopc_avgkern)
  
  end subroutine oopc_deallocAvgkern

  !--------------------------------------------------------------------------
  ! oopc_findAvgkern
  !--------------------------------------------------------------------------
  function oopc_findAvgkern(cstnid,varno,nlev) result(istnid)
    !
    !:Purpose: To find the averaging kernel for an observation if one is
    !          specified. Returns 0 if either not found or not specified.
    !          Combination of STNID, BUFR element and number of vertical levels
    !          to determine association to the observations.
    !
    implicit none
    integer :: istnid ! Index of averaging kernel in oopc_avgkern if found. Zero indicates  averaging kernel not found.

    ! Arguments:
    character(len=12), intent(in) :: cstnid ! station id
    integer, intent(in) :: varno ! BUFR descriptor element
    integer, intent(in) :: nlev  ! number of levels in the observation

    ! Locals:
    integer :: stnidIndex
    logical :: iset

    ! Find stnid with same number of vertical levels, and same BUFR element
          
    istnid=0

    do stnidIndex=1,oopc_avgkern%n_stnid

      ! First compare STNID values allowing for * and blanks in 
      ! oopc_avgkern%stnids(stnidIndex) as wildcards
      iset = utl_stnid_equal(oopc_avgkern%stnids(stnidIndex),CSTNID)

      ! Check if number of levels and BUFR code are equal.
      if (iset) then
        if ( varno == oopc_avgkern%element(stnidIndex) ) then
          if (nlev < 1 .or. ( nlev == oopc_avgkern%n_lvl(stnidIndex) .or. &
	      nlev == oopc_avgkern%n_col(stnidIndex) .or. &
              nlev == oopc_avgkern%n_col(stnidIndex)-1 .or. &
              nlev == oopc_avgkern%n_col(stnidIndex)-2) ) then
	      
            istnid=stnidIndex
            exit
          end if
        end if
      end if
              
    end do

  end function oopc_findAvgkern

  !--------------------------------------------------------------------------
  ! oopc_getAvgkern
  !--------------------------------------------------------------------------
  subroutine oopc_getAvgkern(istnid,nlev,ncol,code,avg_kern)
    !
    !:Purpose: To return averaging kernel for an observation.
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: istnid       ! index of averaging kernel in oopc_avgkern
    character(len=*), intent(in) :: code ! measurement identifier
    integer, intent(in)  :: nlev         ! number of observation levels
    integer, intent(in)  :: ncol         ! number of columns for avg kernel info (without the a priori contribution)
    real(8), intent(out) :: avg_kern(nlev,ncol+2) ! the averaging kernel, plus the possible a priori contribution
                                                  ! and integration weights.

    ! Locals:
    integer :: startIndex,endIndex

    if (istnid > 0 .and. istnid <= oopc_avgkern%n_stnid) then
       
      if (oopc_avgkern%source(istnid) == 0) then
        ! Check number of columns
        if (ncol < oopc_avgkern%n_col(istnid) .or. &
	    ncol+2 > oopc_avgkern%n_col(istnid) ) then
	    
          call utl_abort('oopc_getAvgkern: Inconsistency ' // &
	                 'in avg kern size for ' // oopc_avgkern%stnids(istnid) )
        end if
	      
        ! get averaging kernel from auxiliary file
        startIndex = oopc_avgkern%ibegin(istnid)
        endIndex = nlev*(startIndex+oopc_avgkern%n_col(istnid)-1)
        avg_kern = RESHAPE(oopc_avgkern%rak(startIndex:endIndex), &
	                 (/nlev,oopc_avgkern%n_col(istnid)/),ORDER =(/2,1/))
      else
        ! get averaging kernel from observation file
        avg_kern(1:nlev,1:oopc_avgkern%n_col(istnid)) = &
	  oss_obsdata_get_array2d(oopc_avgkern%obsSubSpace(istnid),code)
      end if

    else
      call utl_abort("oopc_getAvgkern: Invalid station ID index.")
    end if

  end subroutine oopc_getAvgkern

  !----------------------------------- Misc ---------------------------------
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! oopc_dealocInfo
  !--------------------------------------------------------------------------
  subroutine oopc_deallocInfo(info)
    !
    !:Purpose: To deallocate struct_oopc_info instance
    !
    
    implicit none

    type(struct_oopc_info), intent(inout) :: info

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

  end subroutine oopc_deallocInfo

  !--------------------------------------------------------------------------
  ! oopc_diagnOnly
  !--------------------------------------------------------------------------
  logical function oopc_diagnOnly(cfamName,cstnid,varno,nobslev,flag)
    ! 
    !:Purpose: To identify whether or not the obs set identified by the
    !          combination of (cstnid,varno,nobslev) will be assimilated or
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
    if (assim_famNum > 0) then
      do i=1,assim_famNum
        if (assim_fam(i) == cfamName) then
          ifam=i
          exit
        end if
      end do
    end if

    if (ifam == 0) then
      ! assimilate all observations
      oopc_diagnOnly = .false.
      return    
    end if
    
    if (assim_all(ifam)) then
      ! assimilate all observations
      oopc_diagnOnly = .false.
    else if (assim_num(ifam) <= 0) then
      ! assimilate no observations
      oopc_diagnOnly = .true.
    else if (assim_num(ifam) > 0) then
      ! check if this observation is listed in the assim_* arrays
      elemId=0
      do i=1,assim_num(ifam)
        if (utl_stnid_equal(trim(assim_stnid(ifam,i)),trim(cstnid))) then
          if (assim_varno(ifam,i) == 0 .or. assim_varno(ifam,i) == varno) then
            if (assim_nlev(ifam,i) == 0 .or. (nobslev == 1 .and. &
	        assim_nlev(ifam,i) == 1) .or.  &
                (nobslev > 1 .and. assim_nlev(ifam,i) > 1)) then
		    
              elemId=i
              exit
            end if
          end if
        end if
      end do
      oopc_diagnOnly = elemId == 0
    end if
    
    if (oopc_diagnOnly) return

    ! check if the observation integer flag has a bit marked by assim_exclude_flag (same flagging as in filt_suprep)
    if (assim_exclude_nflag(ifam) > 0) then
      do i=1,assim_exclude_nflag(ifam)
        if (btest(flag, 13 - assim_exclude_flag(ifam,i) )) then
          oopc_diagnOnly = .true.
          return
        end if
      end do
    end if

  end function oopc_diagnOnly

  !--------------------------------------------------------------------------
  ! oopc_checkType
  !--------------------------------------------------------------------------
  function oopc_checkType(StnidSet,TypeSet,stnid,type) result(sameType)
    !
    !:Purpose: To determine if specified combination of (stnid,type) found
    !          in (StnidSet,TypeSet).
    !
    implicit none

    logical :: sametype

    ! Arguments:
    character(len=*), intent(in) :: StnidSet(:),TypeSet(:)
    character(len=*), intent(in) :: stnid,type
    
    ! Locals:
    
    integer :: stnidIndex
          
    sameType = .false.
    
    do stnidIndex=1,size(StnidSet)
    
      if ( trim(TypeSet(stnidIndex)) == '' .or. &
           trim(StnidSet(stnidIndex)) == '' ) exit

      if (utl_stnid_equal(StnidSet(stnidIndex),stnid)) then
        if ( trim(StnidSet(stnidIndex)) == trim(stnid) ) then
          if ( trim(TypeSet(stnidIndex)) == trim(type) ) then
            sameType=.true.
          else
            sameType=.false.
          end if
          exit
        else
          if ( trim(TypeSet(stnidIndex)) == trim(type) ) sameType=.true.           
        end if
      end if
       
    end do

    if ( .not.sameType .and. trim(type) == 'default' ) sameType = .true.
    
  end function oopc_checkType

  !--------------------------------------------------------------------------
  ! oopc_getType
  !--------------------------------------------------------------------------
  function oopc_getType(StnidSet,TypeSet,stnid) result(type)
    !
    !:Purpose: To determine "type" for specified "stnid" found in "(StnidSet,Typesef)".
    !
    implicit none

    character(len=15) :: type

    ! Arguments:
    character(len=*), intent(in)  :: StnidSet(:),TypeSet(:)
    character(len=*), intent(in)  :: stnid
    
    ! Locals:
    integer :: stnidIndex
          
    type = 'default'
    
    do stnidIndex=1,size(StnidSet)
    
      if ( trim(TypeSet(stnidIndex)) == '' .or. &
           trim(StnidSet(stnidIndex)) == '' ) exit

      if (utl_stnid_equal(StnidSet(stnidIndex),stnid)) then
        type=trim(TypeSet(stnidIndex))
        if ( trim(StnidSet(stnidIndex)) == trim(stnid) ) exit
      end if
       
    end do
    
  end function oopc_getType

  !--------------------------------------------------------------------------
  !------------------ Routines associated to oopc_efftemp --------------------

  !--------------------------------------------------------------------------
  ! oopc_addEfftempObsdata
  !--------------------------------------------------------------------------
  subroutine oopc_addEfftempObsdata(code,temp_eff)
    !
    !:Purpose: To add effective temperature value to its obsdata object
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: code ! unique identifying code
    real(8), intent(in) :: temp_eff(:)   ! effective temperature

    call oss_obsdata_add_data1d(oopc_efftemp,temp_eff,code,oopc_obsdata_maxsize)
    
  end subroutine oopc_addEfftempObsdata

  !--------------------------------------------------------------------------
  ! oopc_addEfftempObsfile
  !--------------------------------------------------------------------------
  subroutine oopc_addEfftempObsfile()
    !          
    !:Purpose: To add effective temperatures in obs file.
    !
    implicit none

    ! Locals:
    integer :: nrep_modified,varno(1)

    ! If needed, add effective temperature values in obs file for total column measurements

    call oss_obsdata_MPIallgather(oopc_efftemp)
    
    if (oopc_efftemp%nrep > 0) then
      varno(1)=12001
      nrep_modified = obsf_obsSub_update(oopc_efftemp,'CH', &
                      varno(1:max(1,oopc_efftemp%dim2)),bkstp_opt=0, &
                      block_opt='INFO',multi_opt='UNI') 
      write(*,*) 'oopc_addEfftempObsfile: Added ',nrep_modified, &
                 ' effective temperature values in the obs file.'
    end if 

  end subroutine oopc_addEfftempObsfile

  !-------------- CONTROL ROUTINES FOR OBSERVATION OPERATORS ----------------
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! oopc_CHobsoperators
  !--------------------------------------------------------------------------
  subroutine oopc_CHobsoperators(columnTrl, obsSpaceData, mode, &
                                 columnAnlInc_opt, jobs_opt, destObsColumn_opt)
    !
    !:Purpose: To apply the observation operators for chemical constituents.
    !          Mode of operator set by mode (see also kmode below).
    !
    !:Comments:
    !      - See type struct_oopc_obsoperators for description of obsoper elements.
    !      - Currently can only handle the case when nlev_bkgrnd == nlev_inc
    !
    !:Arguments:
    !   :columnTrl:      Column of x_background interpolated to observation
    !                    location. Can have the same vertical levels as the
    !                    trial field (columnhr) or as the increment field
    !                    (columng)
    !   :mode: (kmode retained internally following the switch to mode as input)
    !    +------+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
    !    | mode  |kmode|       Mode of         |             Results               |
    !    |       |     | Observation Operator  |                                   |
    !    +=======+=====+=======================+===================================+
    !    |  nl   |  0  |for general simulation |OmP and total Jo(x_background)     |
    !    |       |     |operator               |for CH. OmP saved in OBS_OMP of    |
    !    |       |     |(non-linear and linear)|obsSpaceData                       |
    !    +-------+-----+-----------------------+-----------------------------------+
    !    | HBHT  |  1  |for identification or  |background error standard dev.     |
    !    |       |     |determination of       |in observation space saved in      |
    !    |       |     |of sqrt(diag(H*B*H^T). |in OBS_HPHT of obsSpaceData        |
    !    |       |     |Depends on the presence|if OmP error std dev not initially
    !    |       |     |of OmP error std dev   |available in OBS_OMPE              |
    !    +-------+-----+-----------------------+-----------------------------------+
    !    |  tl   |  2  |for tangent linear     |Hdx saved in OBS_WORK of           |
    !    |       |     |operator               |obsSpaceData                       |
    !    +------+-----+-----------------------+-----------------------------------+
    !    |adjoint|  3  |for adjoint of tangent |H^T * R^-1 (OmP-Hdx) in            |
    !    |       |     |linear operator        |columnAnlInc_opt                      |
    !    +-------+-----+-----------------------+-----------------------------------+    
    !
    !   :columnAnlInc_opt:  Optional argument for input/output of column of
    !                    increment (column). For kmode=2, used as input for
    !                    increment H_horiz dx interpolated to observation
    !                    location. For kmode=3, used as output for H^T * R^-1
    !                    (OmP-Hdx). Required for kmode=2,3.
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
    !             end do BODY
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
    character(len=*), intent(in) :: mode
    type(struct_columnData), intent(inout), optional :: columnAnlInc_opt
    real(8), intent(out), optional :: jobs_opt
    integer, intent(in), optional :: destObsColumn_opt

    ! Local variables
    real(8) :: zomp,zinc,zoer,zhbht
    integer :: kmode ! Mode of observation operator
    integer, external :: fclos

    ! Obs space local variables

    integer :: headerIndex,bodyIndex,bodyIndex_start,bodyIndex_end
    integer :: icodtyp,obslevIndex,nobslev,varno,maxnumHeaders,headerCount
    integer :: destObsColumn
    character(len=12) :: stnid

    integer, allocatable :: iass(:),flag(:)
    logical, allocatable :: process_obs(:)

    ! Model space profile local variables

    real(8), allocatable :: obs_col(:)
    real(8), pointer :: col(:),model_col(:)
    integer :: nlev_bkgrnd,nlev_inc,modlevIndex
    character(len=2), parameter :: varLevel = 'TH'

    if ( mode == 'nl' ) then
      kmode = 0
    else if ( mode == 'HBHT' ) then
      kmode = 1
    else if ( mode == 'tl' ) then
      kmode = 2
    else if ( mode == 'adjoint' ) then
      kmode = 3
    end if
    
    ! Apply setup on first call
    if (.not.initializedChem) then
      call oopc_setupCH(kmode)
      initializedChem = .true.
    else if ( mode == 'HBHT' .and. .not.bgStats%initialized ) then
      call bcsc_getCovarCH(bgStats)
    end if
    
    if ((kmode == 2 .or. kmode == 3) .and. (.not.present(columnAnlInc_opt))) then
      call utl_abort("oopc_CHobsoperators: columnAnlInc_opt must " // &
                     "be specified for kmode = " // utl_str(kmode))
    end if
     
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
      if (nlev_inc /= nlev_bkgrnd) then
        write(*,*) "oopc_CHobsoperators: nlev_inc ", &
                   "and nlev_bkgrnd not the same: ",nlev_inc, nlev_bkgrnd
      end if
    case(1,3)
      allocate(model_col(nlev_bkgrnd))
    end select

    ! Allocations outside oopc_obsoperInit since this can be done outside the HEADER loop.
    ! See oopc_obsoperInit for assignment of array content.
    
    ! Model obs background, height, TT, and HU profiles.
    allocate(obsoper%trial(nlev_bkgrnd),obsoper%height(nlev_bkgrnd),obsoper%tt(nlev_bkgrnd),obsoper%hu(nlev_bkgrnd))
    ! Model PP and pressure model layer boundaries taken as the middle between model levels.
    allocate(obsoper%pp(nlev_bkgrnd))
    
    ! Determine number of obs
    maxnumHeaders=0
    call obs_set_current_header_list(obsSpaceData,'CH')
    do
      maxnumHeaders=maxnumHeaders+1
      if (obs_getHeaderIndex(obsSpaceData) < 0) exit
    end do

    ! Loop over all header indices of the 'CH' family:
    
    headerCount=0
    call obs_set_current_header_list(obsSpaceData,'CH')
    HEADER: do

      headerCount = headerCount+1
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER
  
      icodtyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
      if (icodtyp /= codtyp_get_codtyp('CHEMREMOTE') .and. &
        icodtyp /= codtyp_get_codtyp('CHEMINSITU')) cycle HEADER
      
      stnid = obs_elem_c(obsSpaceData,'STID',headerIndex)
      bodyIndex_start = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
      bodyIndex_end = bodyIndex_start+obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex)-1

      ! Set number of obs profile elements by removing count of BUFR_SCALE_EXPONENT elements
      nobslev = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex)  
      do bodyIndex=bodyIndex_start,bodyIndex_end
        if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) == &
          BUFR_SCALE_EXPONENT) then
	  nobslev = nobslev-1
	end if
      end do

      ! varno is expected to be the same for all profile points where OBS_VNM value /= BUFR_SCALE_EXPONENT
      do bodyIndex=bodyIndex_start,bodyIndex_end
        if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) &
	    /= BUFR_SCALE_EXPONENT) then
	     
          varno = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
          exit
        end if
      end do

      ! Allocate memory for remaining profile data not in obsoper
      allocate(obs_col(nobslev),iass(nobslev),process_obs(nobslev),flag(nobslev))

      if (allocated(obsoper%success)) deallocate(obsoper%success,obsoper%ixtr)
      allocate(obsoper%success(nobslev),obsoper%ixtr(nobslev))

      ! Check to see if background error variances available
      if (kmode == 1) then
        process_obs(:) = bcsc_StatsExistForVarName(vnl_varnameFromVarnum(varno, &
			 obs_headElem_i(obsSpaceData,OBS_CHM,headerIndex), &
			 modelName))
      end if
 
      ! Prepare for checking if any processing is needed according to initial flag values     
      obslevIndex=0
 
      do bodyIndex=bodyIndex_start,bodyIndex_end
        if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) &
	     /= BUFR_SCALE_EXPONENT) then

          obslevIndex=obslevIndex+1

          obsoper%ixtr(obslevIndex) = &
	    obs_bodyElem_i(obsSpaceData,OBS_XTR,bodyIndex) ! indicates if obs extends outside model profile vertical range
          iass(obslevIndex) = obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) ! indicates if obs is to be assimilated
          flag(obslevIndex) = obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex) ! observation integer flag
               
          ! Indicates if this obs should be processed by oopc_obsoperators
          if (kmode == 1) then
            process_obs(obslevIndex) = obsoper%ixtr(obslevIndex) == 0 &
	             .and. iass(obslevIndex) == obs_assimilated .and. &
	              process_obs(obslevIndex)
          else
            process_obs(obslevIndex) = obsoper%ixtr(obslevIndex) == 0 &
	                    .and. iass(obslevIndex) == obs_assimilated
          end if

        end if
      end do

      ! Initialize processing success flag
      
      obsoper%success(1:nobslev) = process_obs(1:nobslev)
      
      if (all(.not.process_obs)) then

        ! All observations in the profile flagged so can skip obs operator for current measurement

        if (kmode == 3) then
          model_col(:) = 0.0D0
          obsoper%varName = vnl_varnameFromVarnum(varno, &
	    obs_headElem_i(obsSpaceData,OBS_CHM,headerIndex),modelName)
        end if

      else  

        if ( kmode == 1 ) then
          ! Check if "OmP error std dev" is already available
          call obs_set_current_body_list(obsSpaceData,headerIndex)
          obslevIndex=0
          BODYINDEX1: do bodyIndex=bodyIndex_start,bodyIndex_end
	  
            if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) == &
	        BUFR_SCALE_EXPONENT) cycle
		
            obslevIndex = obslevIndex + 1 
            if ( process_obs(obslevIndex) ) then
              if (obs_bodyElem_r(obsSpaceData,OBS_OMPE,bodyIndex) > 0.0d0 ) then
                ! "OmP error std dev" is already available for this measurement.
		! Go to the next measurement.
                  
                ! TEMPORARY: First, estimate OBS_HPHT for storage in output files (in the event it is needed externally for
                ! other purposes (e.g. total column ozone bias correction and corresponding re-doing for marker settings)               
                if ( obs_bodyElem_r(obsSpaceData,OBS_OMPE,bodyIndex) > &
		     1.1d0*obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex) ) then

                  zhbht = sqrt(obs_bodyElem_r(obsSpaceData,OBS_OMPE,bodyIndex)**2 &
		          -obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex)**2)
                  call obs_bodySet_r(obsSpaceData,OBS_HPHT,bodyIndex,zhbht)
                else
                  call obs_bodySet_r(obsSpaceData,OBS_HPHT,bodyIndex, &
		       0.5d0*obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex))
                end if
     
                deallocate(process_obs,obsoper%success,obsoper%ixtr,iass,obs_col,flag)
                  
                cycle HEADER
              else
                ! Proceed with the calc of sqrt(diag(HBHT))
                exit BODYINDEX1
              end if
            end if
          end do BODYINDEX1
        end if

        ! Initialize obsoper variables and allocate arrays        
	call oopc_obsoperInit(obsSpaceData,headerIndex,columnTrl,nlev_bkgrnd,nobslev,kmode,varno,stnid)
 
        ! Initialize model_col, dependent on kmode. Used for input for kmode=0,2, output for kmode=3.
        ! model_col represents for kmode 0) the horizontally interpolated background H_horiz(x_b)
        !                                1) not used
        !                                2) the analysis increment H_horiz dx
        !                                3) the result of applying the adjoint of H_vert 

        select case(kmode)
        case(0)
          model_col => obsoper%trial
        case(2)
          do modlevIndex=1,nlev_inc
            model_col(modlevIndex) = col_getElem(columnAnlInc_opt,modlevIndex,headerIndex,obsoper%varName)
          end do
        case(1,3)
          model_col(:) = 0.0D0
        end select
               
        ! Loop over all body indices (profile elements) to acquire remaining data
      
        obslevIndex=0
        call obs_set_current_body_list(obsSpaceData,headerIndex)
        BODY1: do

          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY1
            
          ! Get position in profile and skip over BUFR_SCALE_EXPONENT elements

          if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) /= &
	      BUFR_SCALE_EXPONENT) then
	      
            obslevIndex=obslevIndex+1
          else
            cycle BODY1
          end if

          ! Get vertical coordinate data. Valid for point data values in profiles.
          ! For layer data values, vertical coordinate data will instead be assigned within oopc_obsoperators.

          obsoper%obslev(obslevIndex) = &
	    obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)

          ! Get normalized increment
          if (kmode == 3) then
            if (iass(obslevIndex) == 1) then
              obs_col(obslevIndex) = obs_bodyElem_r(obsSpaceData,OBS_WORK,bodyIndex)
            else
              obs_col(obslevIndex) = 0.0D0
            end if
          end if

          ! Get obs estimate from trial field (if kmode>1)
            
          if (kmode <= 1) then
            obsoper%obsSpaceTrial(obslevIndex) = 0.0
          else
            ! Store for use by TL and AD of non-linear operators
            obsoper%obsSpaceTrial(obslevIndex) = &
	      obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex) &
              - obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)
          end if
            
        end do BODY1
      
        ! Apply observation operator. model_col,obs_col are the inputs/outputs in model,observation
        ! space, respectively. Other required inputs are in obsoper. Input/output is as follows:
        !
        !    kmode      model_col      obs_col
        !    -----      ---------      -------
        !      0           in            out
        !      1         not used        out
        !      2           in            out
        !      3           out           in

        call oopc_obsoperators(model_col,obs_col,maxnumHeaders,headerCount,kmode)

      end if

      ! Output results
      
      if (kmode == 3) then
        ! Store H^T * R^-1 (OmP-Hdx) in columnInc
                     
        col => col_getColumn(columnAnlInc_opt,headerIndex,obsoper%varName)
        col(1:nlev_bkgrnd) = model_col(1:nlev_bkgrnd)

      else
        ! Store results in obsSpaceData

        obslevIndex=0
        call obs_set_current_body_list(obsSpaceData,headerIndex)
        BODY2: do
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY2

          ! Get position in profile and skip over BUFR_SCALE_EXPONENT elements
          if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) /= &
	      BUFR_SCALE_EXPONENT) then
	      
            obslevIndex=obslevIndex+1
          else
            cycle BODY2
          end if

          ! Check for success in calculations
          if (process_obs(obslevIndex) .and. .not.obsoper%success(obslevIndex)) then
            ! Observation was flagged within this call of oopc_CHobsoperators
            call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
            call obs_bodySet_r(obsSpaceData,OBS_OMP,bodyIndex,0.0D0)
            call obs_bodySet_r(obsSpaceData,OBS_OMA,bodyIndex,0.0D0)
            call obs_bodySet_r(obsSpaceData,OBS_HPHT,bodyIndex,0.0D0)
            call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex,0.0D0)
            call obs_bodySet_i(obsSpaceData,OBS_FLG,bodyIndex, &
	         ibset(obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex),9) )
            cycle BODY2
          else if (iass(obslevIndex) == 0) then
            ! Observation was flagged previous to this call of oopc_CHobsoperators
            call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
            call obs_bodySet_r(obsSpaceData,OBS_OMP,bodyIndex,0.0D0)
            call obs_bodySet_r(obsSpaceData,OBS_OMA,bodyIndex,0.0D0)
            call obs_bodySet_r(obsSpaceData,OBS_HPHT,bodyIndex,0.0D0)
            call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex,0.0D0)
            cycle BODY2
          else if (.not.process_obs(obslevIndex) .and. kmode == 1) then
            call obs_bodySet_r(obsSpaceData,OBS_HPHT,bodyIndex,0.0D0)
            cycle BODY2
          end if

          ! Store result in appropriate location in obsSpaceData
          select case(kmode)
          case(0)
            
            ! Store OmP in OBS_OMP and add to Jo(x_background) of CH.   
                           
            zomp = obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex) - obs_col(obslevIndex)
            call obs_bodySet_r(obsSpaceData,destObsColumn,bodyIndex,zomp)

            if (oopc_diagnOnly('CH',stnid,varno,nobslev,flag(obslevIndex))) then
              ! Observation is for diagnostics and is not to be assimilated
              call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
            else if (present(jobs_opt) .and. iass(obslevIndex) == 1) then
              ! Add to Jo contribution (factor of 0.5 to be applied outside report loop)
              zinc = zomp/obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex)
              jobs_opt = jobs_opt + zinc**2
            end if

          case(1)
            
            ! Background error standard deviations in
            ! observation space, sqrt(diag(H*B_static*H^T)), 
            ! saved in OBS_HPHT of obsSpaceData.
            ! Resulting OmP error std dev estimate saved in OBS_OMPE

            call obs_bodySet_r(obsSpaceData,OBS_HPHT,bodyIndex,obs_col(obslevIndex))
               
            zoer = obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex)
            zomp = sqrt(obs_col(obslevIndex)*obs_col(obslevIndex) + zoer*zoer)
            call obs_bodySet_r(obsSpaceData,OBS_OMPE,bodyIndex,zomp)

          case(2)
            
            !   Store Hdx in OBS_WORK of obsSpaceData               
            call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex,obs_col(obslevIndex))
          end select

        end do BODY2

      end if

      ! Deallocate profile data
      deallocate(process_obs,iass,obs_col,flag)
      call oopc_obsoper_dealloc
  
    end do HEADER

    deallocate(obsoper%trial,obsoper%pp,obsoper%tt)
    deallocate(obsoper%height,obsoper%hu)
    if (kmode.ne.0) deallocate(model_col)
  
    if (present(jobs_opt)) jobs_opt = 0.5d0*jobs_opt
    
  end subroutine oopc_CHobsoperators

  !--------------------------------------------------------------------------
  ! oopc_obsoperInit
  !--------------------------------------------------------------------------
  subroutine oopc_obsoperInit(obsSpaceData,headerIndex,columnTrl, &
                              nmodlev,nobslev,kmode,varno,stnid)
    !
    !:Purpose: To initialize struct_oopc_obsoperators variables and to allocate
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
    !                  trial field (columnTrl) or as the increment field
    !                  (columnAnlInc)
    !     :kmode:      Mode of observation operator  
    !                  - 0 for non-linear/linear model in assimilation (all models
    !                    included are currently linear)
    !                  - 1 for determination of sqrt(diag(H*B*H^T))
    !                  - 2 for tangent linear model
    !                  - 3 for adjoint model 
    !    
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData ! Obs-Space Data object
    integer, intent(in) :: headerIndex ! Measurement index in obsSpaceData
    type(struct_columnData), intent(inout) :: columnTrl
    integer, intent(in) :: nmodlev ! Number of background field (model) levels
    integer, intent(in) :: nobslev ! Number of obs elements (see oopc_obsoper_proceed)
    integer, intent(in) :: kmode   ! Mode of observation operator
    integer, intent(in) :: varno   ! obs unit BUFR number
    character(len=12), intent(in) :: stnid ! Station ID

    ! Locals:
    integer :: bodyIndex ! Measurement element index in obsSpaceDate (see oopc_obsoper_proceed)
    integer :: jl,nmodlev_uv
    real(8), pointer :: col_height_ptr(:)
    real(8), allocatable :: uu(:),vv(:)
    character(len=2), parameter :: varLevel = 'TH'
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
    obsoper%constituentId = obs_headElem_i(obsSpaceData,OBS_CHM,headerIndex)
    
    ! Check if constituent id is recognized (function will abort if not recognized)
    if ( obsoper%constituentId >= 0 .and. obsoper%constituentId < 7000) then
      checkID = vnl_varMassFromVarnum(obsoper%constituentId)
    end if

    obsoper%lat = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex) 
    obsoper%lon = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex) 

    ! Body info that we only need for first point in the profile
    bodyIndex       = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex) 
    ! Index of vertical coordinate type    
    obsoper%vco     = obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex)  
    ! Model field name (NOMVAR value)
    obsoper%varName = vnl_varnameFromVarnum(obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex),obsoper%constituentId,modelName)

    ! Allocate arrays

    allocate(obsoper%obslev(nobslev))       ! Reference vertical levels
    allocate(obsoper%vlayertop(nobslev))    ! Layer tops for layer measurements
    allocate(obsoper%vlayerbottom(nobslev)) ! Layer bottoms for layer measurements
    allocate(obsoper%modlevindexTop(nobslev))  ! Index of highest model level (lowest index) involved with obs element
    allocate(obsoper%modlevindexBot(nobslev))  ! Index of lowest model level (highest index) involved with obs element
    allocate(obsoper%obsSpaceTrial(nobslev))    ! Obs estimate from trial field 

    allocate(obsoper%zh(nobslev,nmodlev))   ! Local model operator H (excluding conversion constants and horizontal interpolation)
    allocate(obsoper%zhp(nobslev,nmodlev))  ! Part of zh that excludes aspects related to vertical resolution

    obsoper%zh(:,:)=0.0D0
    obsoper%zhp(:,:)=0.0D0
    obsoper%modlevindexTop(:)=1
    obsoper%modlevindexBot(:)=nmodlev
    
    if (.not.col_varExist(columnTrl,'TT')) then
      if (oopc_required_field('TT',obsoper%varno)) then
        call utl_abort("chem_opsoper_init: TT required for BUFR code " // trim(utl_str(obsoper%varno)))
      end if
    end if

    ! Get background profiles at observation location
    do jl=1,nmodlev
      obsoper%pp(jl) = col_getPressure(columnTrl,jl,headerIndex,varLevel)
      obsoper%trial(jl) = col_getElem(columnTrl,jl,headerIndex,obsoper%varName)
      obsoper%tt(jl) = col_getElem(columnTrl,jl,headerIndex,'TT')
    end do

    if (col_varExist(columnTrl,'TT').and.col_varExist(columnTrl,'HU') &
        .and.col_varExist(columnTrl,'P0')) then   
        
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
        obsoper%hu(jl) = col_getElem(columnTrl,jl,headerIndex,'HU')  ! lnq was replaced by q
      end do
    else
      obsoper%hu(:)=-1
    end if
    
    ! If applicable, get column upper boundaries for use with total column 
    ! measurements when the related increment profile is to be restricted to 
    ! the lower atmosphere (e.g. troposphere or PBL; when tropo_bound>0 )
    if (obsoper%vco == 4 .and. nobslev == 1 .and. kmode /= 1) then
      if (kmode == 0) then
          
        if (col_varExist(columnTrl,'HU')) then
          if (col_varExist(columnTrl,'UU').and.col_varExist(columnTrl,'VV')) then
            nmodlev_uv=col_getNumLev(columnTrl,'MM')
            allocate(uu(nmodlev_uv),vv(nmodlev_uv))
            do jl=1,nmodlev_uv
              uu(jl) = col_getElem(columnTrl,jl,headerIndex,'UU')
              vv(jl) = col_getElem(columnTrl,jl,headerIndex,'VV')
            end do
            obsoper%columnBound = oopc_getColBoundary(obsoper%constituentId,&
	      nmodlev,obsoper%pp,obsoper%tt,obsoper%height, &
	      hu_opt=obsoper%hu,uu_opt=uu,vv_opt=vv)
            deallocate(uu,vv)
          else 
            obsoper%columnBound = oopc_getColBoundary(obsoper%constituentId,nmodlev,obsoper%pp,obsoper%tt,obsoper%height,hu_opt=obsoper%hu)   
          end if

        else
          obsoper%columnBound = oopc_getColBoundary(obsoper%constituentId,nmodlev,obsoper%pp,obsoper%tt,obsoper%height)
        end if

        call oopc_addColBoundary(headerIndex,obsoper%columnBound)  ! save boundary for kmode>0 calls using headerIndex

      else
        obsoper%columnBound = oopc_retrieveColBoundary(headerIndex)
      end if
    else
      obsoper%columnBound = -1.
    end if

  end subroutine oopc_obsoperInit

  !--------------------------------------------------------------------------
  ! oopc_obsoper_dealloc
  !--------------------------------------------------------------------------
  subroutine oopc_obsoper_dealloc
    !
    !:Purpose: To deallocate arrays for struct_oopc_obsoperators.
    !
    !:Comments:
    !           - Deallocation of arrays that are dependent on only nmodlev have
    !             been moved outside this subroutine so that they are
    !             deallocated only once.
    !
    implicit none

    ! Arguments:
    ! type(struct_oopc_obsoperators), intent(inout) :: obsoper ! (see module oopc_obsoper) 

    if (allocated(obsoper%obslev))       deallocate(obsoper%obslev)
    if (allocated(obsoper%vlayertop))    deallocate(obsoper%vlayertop)
    if (allocated(obsoper%vlayerbottom)) deallocate(obsoper%vlayerbottom)
    if (allocated(obsoper%zh))           deallocate(obsoper%zh)
    if (allocated(obsoper%zhp))          deallocate(obsoper%zhp)
    if (allocated(obsoper%modlevindexTop))  deallocate(obsoper%modlevindexTop)
    if (allocated(obsoper%modlevindexBot))  deallocate(obsoper%modlevindexBot)
    if (allocated(obsoper%obsSpaceTrial))  deallocate(obsoper%obsSpaceTrial)
    if (allocated(obsoper%success))      deallocate(obsoper%success)
    if (allocated(obsoper%ixtr))         deallocate(obsoper%ixtr)

  end subroutine oopc_obsoper_dealloc
    
  !---------------------- Routines for observation operators ----------------
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! oopc_obsoperators
  !--------------------------------------------------------------------------
  subroutine oopc_obsoperators(model_col,obs_col,maxnumHeaders,headerCount,kmode)
    !
    !:Purpose: To apply observation operator for indicated observation data and
    !          condition.
    !
    !:Arguments:
    !     :kmode: mode of observation operator
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
    !     :maxnumHeaders:  Total number of CH obs for this CPU
    !     :headerCount:    Obs number up to maxnumHeaders
    !
    ! Further changes required for generalization:
    !
    !  1) Add layer average operators.
    !  2) Add AOD operators (summation over model layers).
    !  3) Add option to include use of obs error correlation matrix for kmode=2,3 
    !     (This may/will need to be done in oop_Hchm and oop_HTchm where the division 
    !      by stddev_obs is applied. A new routine will be needed for this
    !      operation - and others for reading the correlation matrices similarly to the
    !      averaging kernels.)
    !
    ! Comments:
    !      - When kmode=0, call from oopc_CHobsoperators passes model_col as a pointer to
    !        obsoper%trial.
    !      - Does not yet account for potential future applications of obs 
    !        vertical correlation matrices.
    !      - Potential specification of background error std. dev. (fdeStddev(:,1:2)) and correlation matrices 
    !        for the ensemble-based and lam cases to/could be done when stats for these become in use with constituents.
    !
    implicit none

    ! Arguments:
    ! type(struct_oopc_obsoperators), intent(inout) :: obsoper ! (see module oopc_obsoper) 

    ! I/O arguments: obs space variables
    
    integer, intent(in) :: kmode,maxnumHeaders,headerCount

    ! I/O arguments: model space profile data and others

    real(8), intent(inout) :: model_col(obsoper%nmodlev), obs_col(obsoper%nobslev)

    ! Locals: 
    logical :: successLocal(obsoper%nobslev) 
    real(8) :: zwork(obsoper%nmodlev),unit_conversion(obsoper%nmodlev)
    real(8) :: fdeStddev(obsoper%nmodlev,2),temp_eff(1)
    real(8), allocatable :: avg_kern(:,:)
    integer :: obslevIndex,modlevIndex,varIndex,nobslevOriginal,varnoOriginal
    integer, parameter :: code_len=90
    character(len=code_len) :: code    ! Must be at least as large as oopc_code_len
     
    if (code_len < oss_obsdata_code_len()) then
      call utl_abort('oopc_obsoperators: Length of code string' // &
                     ' needs to be increased to ' // &
	             trim(utl_str(oss_obsdata_code_len())))
    end if
 
    ! Determine if layer boundaries are assigned to this data source OR, for obsoper%nobslev = 1, 
    ! if the number of observation operator levels differ (i.e. use of 1D averaging kernels) 
    ! If so, obtain them for use in this routine. 
    ! Routine provides obsoper%layerIdentified,
    !                  obsoper%vlayertop(nobslev), 
    !                  obsoper%vlayerbottom(nobslev)
    ! and resets obsoper%nobslev, obsoper%obslev for alternative work vertical levels. 
    ! The original obsoper%nobslev is saved as nobslevOriginal
 
    nobslevOriginal = obsoper%nobslev
    varnoOriginal = obsoper%varno
    call oopc_getLevels

    ! Prepare observation operator
  
    call oopc_prepareOperator(model_col,avg_kern,nobslevOriginal,maxnumHeaders, &
                              headerCount,kmode,unit_conversion,successLocal)

    ! Finalize required quantities depending on kmode
   
    select case(kmode)

    case(0)

      ! Finalize non-linear/linear operator step
     
      if (obsoper%iavgkern > 0) then
        if (oopc_checkType(oopc_avgkern%stnids, &
	                   oopc_avgkern%type,obsoper%stnid,'log')) then
          
          zwork(1:obsoper%nobslev) = 0.0d0
          do obslevIndex=1,obsoper%nobslev
            if (successLocal(obslevIndex)) &
              zwork(obslevIndex)=log(dot_product(obsoper%zh(obslevIndex, &
	        obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex)), &
                model_col(obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex)) ))
          end do

          do obslevIndex=1,oopc_avgkern%n_lvl(obsoper%iavgkern)
            if (obsoper%success(obslevIndex)) then
              obs_col(obslevIndex)= dot_product(avg_kern(obslevIndex, &
	                            1:obsoper%nobslev),zwork(1:obsoper%nobslev))
            end if
	    
            ! Add a priori contribution when provided
            if (oopc_avgkern%n_col(obsoper%iavgkern) >= obsoper%nobslev+1) then
	      obs_col(obslevIndex) = obs_col(obslevIndex) + &
	                             avg_kern(obslevIndex,obsoper%nobslev+1)
            end if
	    
            obs_col(obslevIndex) = exp(obs_col(obslevIndex))
          end do
         
          if ( oopc_avgkern%n_col(obsoper%iavgkern) == obsoper%nobslev+2 &
	       .and. nobslevOriginal == 1 ) then
	    
            obs_col(1) = sum( avg_kern(1:oopc_avgkern%n_lvl(obsoper%iavgkern), &
	      obsoper%nobslev+2)*obs_col(1:oopc_avgkern%n_lvl(obsoper%iavgkern)) )
	  end if
        else 

          ! Standard treatment for linear model
            
          do obslevIndex=1,oopc_avgkern%n_lvl(obsoper%iavgkern) 
            if (obsoper%success(obslevIndex)) then
              obs_col(obslevIndex)=dot_product(obsoper%zh(obslevIndex,&
	        obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex)), &
                model_col(obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex)))
            end if
            ! Add a priori contribution when provided
            if ( oopc_avgkern%n_col(obsoper%iavgkern) >= obsoper%nobslev+1) then
              obs_col(obslevIndex) = obs_col(obslevIndex) +  &
	      avg_kern(obslevIndex,obsoper%nobslev+1)
            end if
          end do
         
          ! Account for integration via weighted summation
          if ( oopc_avgkern%n_col(obsoper%iavgkern) == obsoper%nobslev+2 &
	       .and. nobslevOriginal == 1 .and. &
	       oopc_avgkern%n_lvl(obsoper%iavgkern) > 1 ) then
	    
            obs_col(1) = sum( avg_kern(1:oopc_avgkern%n_lvl(obsoper%iavgkern), &
	      obsoper%nobslev+2)*obs_col(1:oopc_avgkern%n_lvl(obsoper%iavgkern)) )
          end if
        end if

      else

        ! Standard treatment for linear model
            
        do obslevIndex=1,nobslevOriginal
          if (obsoper%success(obslevIndex)) then
            obs_col(obslevIndex)=dot_product(obsoper%zh(obslevIndex, &
	      obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex)), &
              model_col(obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex)))
	  end if
        end do
       
      end if

      ! Calculate concentration-weighted effective temperature (for output purpose)

      if ((obsoper%constituentId >= 0 .and. obsoper%constituentId < 7000) .and.  &
        trim(obsoper%operatorCategory) == 'Integ' .and. obsoper%nobslev == 1 .and. &
	obsoper%vco == 4 .and. obsoper%success(1)) then
	
        if (all(obsoper%tt > 0.0) .and. obs_col(1) > 0.0) then
          temp_eff(1)=dot_product(obsoper%zh(1, &
	    obsoper%modlevindexTop(1):obsoper%modlevindexBot(1)), &
            obsoper%tt(obsoper%modlevindexTop(1):obsoper%modlevindexBot(1))* &
	    model_col(obsoper%modlevindexTop(1):obsoper%modlevindexBot(1))) &
            /obs_col(1)
          code=oss_obsdata_get_header_code(obsoper%lon,obsoper%lat,&
	    obsoper%date,obsoper%hhmm,obsoper%stnid)
          call oopc_addEfftempObsdata(code,temp_eff)
        end if
      end if
                       
    case(1)

      ! Compute sqrt(diag(H*B*H^T))

      ! Apply unit conversion to observation operator
      do obslevIndex=1,nobslevOriginal
        if (obsoper%success(obslevIndex)) then
          obsoper%zh(obslevIndex,:) = unit_conversion * obsoper%zh(obslevIndex,:)
        end if
      end do

      ! Get background error std dev profile at obs location
      fdeStddev(:,:)=0
      call bcsc_getBgStddev(obsoper%varName,obsoper%nmodlev, &
                            obsoper%lat,obsoper%lon,fdeStddev(:,1)) 

      ! Identify variable position index in background error correlation matrices

      varIndex=1
      do while (trim(bgStats%varNameList(varIndex)) /= '') 
        if (trim(bgStats%varNameList(varIndex)) == trim(obsoper%varName)) exit
        varIndex=varIndex+1
      end do
      if (trim(bgStats%varNameList(varIndex)) == '') then
        call utl_abort('oopc_obsoperators: Correlation matrix not found for ' &
	               // trim(obsoper%varName) )
      end if

      do obslevIndex=1,nobslevOriginal
        if (obsoper%success(obslevIndex)) then
          do modlevIndex=obsoper%modlevindexTop(obslevIndex), &
	                  obsoper%modlevindexBot(obslevIndex)
		  
            zwork(modlevIndex)=sum(obsoper%zh(obslevIndex, &
	      obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex)) &
              *bgStats%corvert(modlevIndex,obsoper%modlevindexTop(obslevIndex): &
	        obsoper%modlevindexBot(obslevIndex),varIndex) &
              *fdeStddev(obsoper%modlevindexTop(obslevIndex): &
	        obsoper%modlevindexBot(obslevIndex),1)) &
              *fdeStddev(modlevIndex,1)
          end do
          obs_col(obslevIndex)=sum(obsoper%zh(obslevIndex, &
	    obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex)) &
            *zwork(obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex)))
    
          obs_col(obslevIndex) = sqrt(obs_col(obslevIndex))  ! save as sqrt(h*B*h^T)
        else
          obs_col(obslevIndex) = 0.0
        end if
      end do

    case(2)

      ! Finalize TL operator application

      ! Standard treatment for linear and linearized model
            
      if (obsoper%iavgkern > 0) then

        do obslevIndex=1,oopc_avgkern%n_lvl(obsoper%iavgkern) 
          if (obsoper%success(obslevIndex)) then
            obs_col(obslevIndex)=dot_product(obsoper%zh(obslevIndex, &
	      obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex)), &
              model_col(obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex)))
	  end if 
        end do
         
        ! Account for integration via weighted summation
        if ( oopc_avgkern%n_col(obsoper%iavgkern) == obsoper%nobslev+2 &
	     .and. nobslevOriginal == 1 .and. &
	     oopc_avgkern%n_lvl(obsoper%iavgkern) > 1 ) then
	  
          obs_col(1) = sum( avg_kern(1:oopc_avgkern%n_lvl(obsoper%iavgkern), &
	    obsoper%nobslev+2)*obs_col(1:oopc_avgkern%n_lvl(obsoper%iavgkern)) )

        end if
      else
     
        do obslevIndex=1,nobslevOriginal
          if (obsoper%success(obslevIndex)) then
            obs_col(obslevIndex)=dot_product(obsoper%zh(obslevIndex, &
	      obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex)), &
              model_col(obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex)))
          else
            obs_col(obslevIndex)=0.0d0
          end if
        end do
     
      end if
     
    case(3)

      ! H^T*grad contribution from adjoint of tangent linear model.

      model_col(:) = 0.0
       
      ! Standard treatment for linear and linearized model
       
      if (obsoper%iavgkern > 0) then

        if ( oopc_avgkern%n_col(obsoper%iavgkern) == obsoper%nobslev+2 &
	     .and. nobslevOriginal == 1 .and. &
	     oopc_avgkern%n_lvl(obsoper%iavgkern) > 1 ) then
       
          ! Account for integration via weighted summation
          zwork(:)=0.0d0         
          do obslevIndex=obsoper%modlevindexTop(obslevIndex),obsoper%modlevindexBot(obslevIndex)
            zwork(obslevIndex)=  &
	      sum(obsoper%zh(1:oopc_avgkern%n_lvl(obsoper%iavgkern),obslevIndex)* &
              avg_kern(1:oopc_avgkern%n_lvl(obsoper%iavgkern),obsoper%nobslev+2) )
          end do
          obsoper%zh(1,:) = 0.0d0
          obsoper%zh(1,obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex)) = &
            zwork(obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex))
        end if
      end if          

      do obslevIndex=1,nobslevOriginal
        if (obsoper%success(obslevIndex)) then
          zwork(:)=0.0D0
          zwork(obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex)) = &
             obs_col(obslevIndex)*obsoper%zh(obslevIndex, &
	     obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex))
                
          call oopc_convertUnits(zwork,incr_opt=.true.)
          
          model_col(obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex)) = &
             model_col(obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex)) + &
             zwork(obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex))
        end if
      end do
     
    end select

    if (obsoper%iavgkern > 0) deallocate(avg_kern)
     
  end subroutine oopc_obsoperators

  !--------------------------------------------------------------------------
  ! oopc_prepareOperator
  !--------------------------------------------------------------------------
  subroutine oopc_prepareOperator(model_col,avg_kern,nobslevOriginal, &
                                  maxnumHeaders,headerCount,kmode,unit_conversion,successLocal)
    !
    !:Purpose: To prepare observation operator
    !
    
    implicit none
    
    ! Arguments:
    integer, intent(in) :: kmode ! Mode of observation operator
    integer, intent(in) :: nobslevOriginal ! Number of actual obs elements in measurement
    integer, intent(in) :: maxnumHeaders   ! Total number of CH obs for this CPU
    integer, intent(in) :: headerCount     ! Obs counter/index
    real(8), intent(inout) :: model_col(obsoper%nmodlev)
    real(8), intent(out)   :: unit_conversion(obsoper%nmodlev)
    logical, intent(out)   :: successLocal(obsoper%nobslev)
    real(8), intent(out), allocatable :: avg_kern(:,:)

    ! Locals:  
    real(8), allocatable :: press_obs(:)
    integer, allocatable :: ixtrLocal(:)
    real(8) :: fdeStddev(obsoper%nmodlev,2),zhmin
    integer :: obslevIndex,modlevIndex
    integer, parameter :: code_len=90
    character(len=code_len) :: code    ! Must be at least as large as oopc_code_len
    character(len=20) :: message
    logical :: successDepot

    ! Check if obs BUFR element needs to be changed with use of averaging kernels
    
    if ( obsoper%iavgkern > 0 ) then
      if ( oopc_avgkern%n_col(obsoper%iavgkern) == obsoper%nobslev+2 &
           .and. nobslevOriginal == 1 ) then
	   
        obsoper%varno=oopc_avgkern%ProfElement(obsoper%iavgkern)
      end if
    end if

    ! Identify observation operator based on observation units and presence or
    ! not of layer boundaries

    if (bufr_IsIntegral(obsoper%varno)) then
      if (.not.obsoper%layerIdentified) then
        write(*,*)   '----------------------------------------------------------'
        write(*,*)   'STNID, BUFR index, nobslev: ',obsoper%stnid,' ', &
	             obsoper%varno,obsoper%nobslev
        call utl_abort('oopc_obsoperators: Required layer boundaries not available!')
      else      
        ! Vertical integration operator
        obsoper%operatorCategory='Integ'
      end if
    else if (obsoper%layerIdentified) then
      ! Layer averaging operator
      obsoper%operatorCategory='LayerAvg'
    else if (obsoper%vco == 5 .and. obsoper%nobslev == 1) then  
      ! Surface point (in-situ) measurement
      obsoper%operatorCategory='Surface'
    else    
      ! Vertical interpolation operator
      obsoper%operatorCategory='Interp'
    end if  

    ! Look for pre-calculated operator

    successDepot = oopc_operatorDepot(headerCount,maxnumHeaders,kmode,'get')

    if ( successDepot ) then

      ! Get averaging kernels if requested

      if (obsoper%iavgkern > 0) then

        if (allocated(avg_kern)) deallocate(avg_kern)
        allocate(avg_kern(oopc_avgkern%n_lvl(obsoper%iavgkern), &
	                  oopc_avgkern%n_col(obsoper%iavgkern)))
    
        code=oss_obsdata_get_header_code(obsoper%lon,obsoper%lat,obsoper%date, &
	  obsoper%hhmm,obsoper%stnid)
        call oopc_getAvgkern(obsoper%iavgkern,oopc_avgkern%n_lvl(obsoper%iavgkern), &
                             oopc_avgkern%n_col(obsoper%iavgkern),code,avg_kern)
      end if

      ! Apply unit conversion (apply later when kmode=3)
      if (kmode == 2) call oopc_convertUnits(model_col,incr_opt=.true.)

      return
    end if

    ! Indicate if the generalized innovation operator is to be applied.

    obsoper%applyGenOper=.false.
    if (obsoper%constituentId >= 0 .and. &
      oopc_checkType(operatorSubType(1,:),operatorSubType(2,:), &
      obsoper%stnid,'genOper')) then 
      
      if (kmode /= 1 .and. (trim(obsoper%operatorCategory) == 'Integ' .or. &
	                    trim(obsoper%operatorCategory) == 'LayerAvg'   )) then
       
        if ( kmode == 2) then
          ! Set reference profiles for use with generalized innovation operator
	  ! when kmode>=2
          call oopc_addToProfileSet(oopc_climFields,oopc_bgRef, &
	         oopc_constituentsSize,2,obsoper%nmodlev,obsoper%pp, &
		 obsoper%height,obsoper%lat,obsoper%lon, obsoper%obs_index, &
		 oopc_obsdata_maxsize,varKind_opt='CH', &
                 varNumber_opt=obsoper%constituentId, &
		 tt_opt=obsoper%tt,hu_opt=obsoper%hu)
		 
          ! Get background error std dev profile at obs locations
          fdeStddev(:,:)=0.D0
	  
          call bcsc_getBgStddev(obsoper%varName,obsoper%nmodlev,obsoper%lat, &
	    obsoper%lon,fdeStddev(:,1),vlev_opt=obsoper%pp) 

          call bcsc_addBgStddev(obsoper%obs_index,fdeStddev, &
	    oopc_obsdata_maxsize)
  
        end if 	
        if (kmode >= 2) obsoper%applyGenOper = .true.

      end if
    end if
    
    ! Apply unit conversion (apply unit conversion later for kmode=3)

    select case(kmode)
    case(0)
      ! Perform transformation and unit conversion on obsoper%trial.
      ! Note that model_col => obsoper%trial for kmode=0.
      call oopc_convertUnits(obsoper%trial)
    case(1)
      ! Save the conversion factor in <unit_conversion>
      unit_conversion(:) = 1.0
      call oopc_convertUnits(unit_conversion,incr_opt=.true.)
    case(2)
      call oopc_convertUnits(model_col,incr_opt=.true.)
    end select

    if (obsoper%applyGenOper) then
      ! Perform unit conversion on obsoper%trial when applying the generalized
      ! obs operator for kmode=2,3. Keep obsoper%trial in ug/kg in this case.
      call oopc_convertUnits(obsoper%trial,ppb_opt=.true.)
    end if

    if (allocated(press_obs)) deallocate(press_obs,ixtrLocal)
    allocate(press_obs(obsoper%nobslev),ixtrLocal(obsoper%nobslev))
    if (nobslevOriginal == obsoper%nobslev) then
      ixtrLocal(:) = obsoper%ixtr(:)
      successLocal(:) = obsoper%success(:)
    else
      if (any(obsoper%ixtr(:) == 0)) then
        ixtrLocal(:)=0
      else
        ixtrLocal(:)=1
      end if
      if (any(obsoper%success)) then
        successLocal(:) = .true.
      else
        successLocal(:) = .false.
      end if
    end if 

    ! Convert observation vertical coordinate value(s) to pressure if needed
    
    select case(obsoper%vco)
    case(1)
      ! Convert altitude to pressure
      if (trim(obsoper%operatorCategory) == 'Interp') then
        press_obs = phf_convert_z_to_pressure(obsoper%obslev,obsoper%height, &
	            obsoper%pp, &
                    obsoper%nobslev,obsoper%nmodlev,obsoper%lat,successLocal)
        ! Allows for obs levels below the lowest TH level and above the surface
        where (obsoper%obslev(1:obsoper%nobslev) < &
	  obsoper%height(obsoper%nmodlev)) &
	  press_obs(1:obsoper%nobslev)= obsoper%pp(obsoper%nmodlev)
	  
      else if (trim(obsoper%operatorCategory) == 'Integ' .or. &
               trim(obsoper%operatorCategory) == 'LayerAvg') then
        obsoper%vlayertop = phf_convert_z_to_pressure(obsoper%vlayertop, &
	                    obsoper%height, &
                            obsoper%pp,obsoper%nobslev,obsoper%nmodlev, &
			    obsoper%lat,successLocal)
        obsoper%vlayerbottom = phf_convert_z_to_pressure(obsoper%vlayerbottom, &
	                       obsoper%height, &
                               obsoper%pp,obsoper%nobslev,obsoper%nmodlev, &
			       obsoper%lat,successLocal)
      end if
    case(2)
      ! Pressure, no conversion needed
      if (trim(obsoper%operatorCategory) == 'Interp') press_obs = obsoper%obslev
    case(4,5)
      ! No actions taken
    case default
      call utl_abort("oopc_obsoperators: vertical coordinate type vco = " &
                      // trim(utl_str(obsoper%vco)) //  &
		      " not available for this operator.")
    end select

    ! Determine if averaging kernel is to be applied

    if (nobslevOriginal == obsoper%nobslev) then
      obsoper%iavgkern = oopc_findAvgkern(obsoper%stnid,obsoper%varno, &
                         obsoper%nobslev)
    else
      obsoper%iavgkern = oopc_findAvgkern(obsoper%stnid,obsoper%varno,0)
    end if
  
    ! Apply appropriate core observation operator
   
    select case(trim(obsoper%operatorCategory))
    case('Interp')

      ! Vertical interpolation operator
     
      if ( obsoper%iavgkern /= 0 ) then
        message = 'doAll&noExtrap'
      else
        message = 'noExtrap'
      end if
      call ppo_vertInterpWgts(obsoper%pp,press_obs,obsoper%nmodlev, &
           obsoper%nobslev,obsoper%zh,obsoper%modlevindexTop, &
	   obsoper%modlevindexBot,method_opt=oopc_getType(operatorSubType(1,:), &
	   operatorSubType(2,:),obsoper%stnid),skipType_opt=message, &
	   outbound_opt=ixtrLocal,success_opt=successLocal)
 
    case('Surface')

      ! Surface point measurement

      ! Set weight to unity for lowest model level.
      obsoper%zh(1,1:obsoper%nmodlev-1)=0.0
      obsoper%zh(1,obsoper%nmodlev) = 1.0

      ! Set range of elements for model vertical levels
      obsoper%modlevindexTop(1) = obsoper%nmodlev
      obsoper%modlevindexBot(1) = obsoper%nmodlev 
    
    case('Integ')
  
      ! Layer integration operator

      if ( obsoper%iavgkern /= 0 ) then
        message = 'doAll&noExtrap'
      else
        message = 'default'
      end if
      
      call oopc_vertObsLayersWgts('integ',ixtrLocal,successLocal,kmode,message)

    case('LayerAvg')
    
      ! Layer averaging operator

      if ( obsoper%iavgkern /= 0 ) then
        message = 'doAll&noExtrap'
      else
        message = 'default'
      end if

      call oopc_vertObsLayersWgts('avg',ixtrLocal,successLocal,kmode,message)

    end select

    if (nobslevOriginal == obsoper%nobslev) then
      obsoper%ixtr(:) = ixtrLocal(:)
      obsoper%success(:) = successLocal(:)
    else
      if (all(ixtrLocal(:) == 0)) then
        obsoper%ixtr(:)=0
      else
        obsoper%ixtr(:)=1
      end if
      if (all(successLocal)) then
        obsoper%success(:) = .true.
      else
        obsoper%success(:) = .false.
      end if
    end if 
    deallocate(press_obs,ixtrLocal)
  
    ! Apply averaging kernels if requested

    if (obsoper%iavgkern > 0) then

      allocate(avg_kern(oopc_avgkern%n_lvl(obsoper%iavgkern), &
               oopc_avgkern%n_col(obsoper%iavgkern)))
       
      code=oss_obsdata_get_header_code(obsoper%lon,obsoper%lat, &
                                       obsoper%date,obsoper%hhmm,obsoper%stnid)
      call oopc_getAvgkern(obsoper%iavgkern, &
                           oopc_avgkern%n_lvl(obsoper%iavgkern), &
                           oopc_avgkern%n_col(obsoper%iavgkern),code,avg_kern)
     
      if (oopc_checkType(oopc_avgkern%stnids,oopc_avgkern%type,obsoper%stnid, &
                         'default')) then
        do obslevIndex=1,oopc_avgkern%n_lvl(obsoper%iavgkern) 
          if (obsoper%success(obslevIndex)) then

            ! Apply averaging kernels to observation operator(s)
            
            obsoper%zh(obslevIndex,:) = matmul(avg_kern(obslevIndex, &
	                                1:obsoper%nobslev),obsoper%zh(:,:))
            if (obsoper%applyGenOper) then
	      obsoper%zhp(obslevIndex,:) = &
	        matmul(avg_kern(obslevIndex,1:obsoper%nobslev),obsoper%zhp(:,:))
	    end if
           
            ! Extend vertical range of obs operator according to the influence of
            ! the averaging kernel. Either extend to the entire model vertical range
            ! (commented out below) or to the vertical range with non-negligable values.

            ! obsoper%modlevindexTop(obslevIndex) = 1
            ! obsoper%modlevindexBot(obslevIndex) = obsoper%nmodlev

            zhmin=1.0D-10*maxval(abs(obsoper%zh(obslevIndex,:)))
            do modlevIndex=1,obsoper%modlevindexTop(obslevIndex)
              if (abs(obsoper%zh(obslevIndex,modlevIndex)) > zhmin) exit
            end do
            if (modlevIndex > obsoper%modlevindexTop(obslevIndex)) then
	      modlevIndex=obsoper%modlevindexTop(obslevIndex)
	    end if
            obsoper%modlevindexTop(obslevIndex) = modlevIndex
            do modlevIndex=obsoper%nmodlev,obsoper%modlevindexBot(obslevIndex),-1
              if (abs(obsoper%zh(obslevIndex,modlevIndex)) > zhmin) exit
            end do
            if (modlevIndex.lt.obsoper%modlevindexBot(obslevIndex)) then
	      modlevIndex=obsoper%modlevindexBot(obslevIndex)
	    end if
            obsoper%modlevindexBot(obslevIndex) = modlevIndex
           
          end if
        end do
 
      else if (oopc_checkType(oopc_avgkern%stnids,oopc_avgkern%type, &
               obsoper%stnid,'log')) then
              
        if (kmode == 0) then
        
          ! Apply log-space averaging kernel below - no transformation needed here
          ! Do not merge the averaging kernel (avgkern) and the vertical interpolator (zh) 
          
        else 
              
          ! Apply linearization of operator involving the log-space averaging kernel 
                      
          do obslevIndex=1,obsoper%nobslev
            avg_kern(1:oopc_avgkern%n_lvl(obsoper%iavgkern),obslevIndex) = &
	      avg_kern(1:oopc_avgkern%n_lvl(obsoper%iavgkern),obslevIndex) /   &
              dot_product( obsoper%zh(obslevIndex, &
	      obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex)), &
              obsoper%trial(obsoper%modlevindexTop(obslevIndex): &
	                    obsoper%modlevindexBot(obslevIndex)) )
          end do                                 

          do obslevIndex=1,oopc_avgkern%n_lvl(obsoper%iavgkern) 
            if (obsoper%success(obslevIndex)) then
              avg_kern(obslevIndex,1:obsoper%nobslev) =  &
	        obsoper%obsSpaceTrial(obslevIndex)*avg_kern(obslevIndex, &
		1:obsoper%nobslev)

              ! Merge the averaging kernel matrix (avgkern) and the 
	      ! vertical interpolator (initial zh) 
              obsoper%zh(obslevIndex,:) = matmul(avg_kern(obslevIndex, &
	        1:obsoper%nobslev),obsoper%zh(:,:))
             
            end if
          end do
           
        end if  
       
        ! Note that obsoper%applyGenOper=.true. is not set up for this case
        if (obsoper%applyGenOper) then
	  call utl_abort('prepareOperator: Log ' // &
	                 'space averaging kernels not currently usable with ' // &
	                 'obsoper%applyGenOper=.true.')
	end if
      else
        call utl_abort('oopc_prepareOperator: This averaging kernel ' // &
	               'application not yet available')
      end if
    else
      if ( nobslevOriginal /= obsoper%nobslev ) then
        call utl_abort('oops_prepareOperator: Case of differing obs ' // &
	               'and calculation levels not recognized.')
      end if
    end if

    ! Apply generalized innovation operator if requested

    if (obsoper%applyGenOper) call oopc_genOper(kmode)

    ! Save operator if needed
    
    successDepot = oopc_operatorDepot(headerCount,maxnumHeaders,kmode,'save')
    
  end subroutine oopc_prepareOperator

  !--------------------------------------------------------------------------
  ! oopc_operatorDepot
  !--------------------------------------------------------------------------
  function oopc_operatorDepot(headerCount,maxnumHeaders,kmode,action) result(success)
    !
    !:Purpose: To save or get previously calculated linear observation operator
    !          for use by TL and AD operators
    !
    implicit none
    
    ! Arguments:
    integer, intent(in) :: kmode           ! Mode of observation operator
    integer, intent(in) :: headerCount     ! Obs counter/index
    integer, intent(in) :: maxnumHeaders   ! Total number of CH obs for this CPU
    character(len=*), intent(in) :: action ! 'save' or 'get' operator for current obs
    logical :: success                     ! Succes of action succeeded  

    ! Locals:  
    integer, save :: maxnumOperators
    type(struct_oopc_operatorsDepot), allocatable, save :: operators(:)
    
    success = .false.
        
    if ( .not.oopc_storeOperators .or. kmode < 2 .or. headerCount <= 0 .or. &
         ( any(oopc_tropo_mode(:) >= 1) .and. &
	   trim(obsoper%operatorCategory) == 'Integ') ) return
    
    if ( trim(action) == 'get' .and. initializedOperators .and. kmode >=2 ) then
    
      ! Check for available pre-calculated operator

      if (operators(headerCount)%nobslev > 0) then

        ! Update arrays/values that might have changed or would otherwise be set 
        ! below when not stored earlier.

        obsoper%trial(:) = operators(headerCount)%trial(:)
        obsoper%vlayertop(:) = operators(headerCount)%vlayertop(:)
        obsoper%vlayerbottom(:) = operators(headerCount)%vlayerbottom(:)
        obsoper%modlevindexTop(:) = operators(headerCount)%modlevindexTop(:)
        obsoper%modlevindexBot(:) = operators(headerCount)%modlevindexBot(:)
        obsoper%zh(:,:) = operators(headerCount)%zh(:,:)
        obsoper%ixtr(:) = operators(headerCount)%ixtr(:)
        obsoper%success(:) = operators(headerCount)%success(:)
        obsoper%iavgkern = operators(headerCount)%iavgkern
        obsoper%applyGenOper = operators(headerCount)%applyGenOper

        if ( obsoper%applyGenOper ) then
	  obsoper%zhp(:,:) = operators(headerCount)%zhp(:,:) 
	end if

        success = .true. 
      else
        ! Calculations will be performed in th calling routine to prepare 
	! the operator fields    
      end if  
         
    else if ( trim(action) == 'save' .and. kmode == 2 ) then 
    
      ! Save operator for later use
      
      if ( .not.initializedOperators ) then
        if (allocated(operators)) deallocate(operators)
        maxnumOperators = maxnumHeaders
        allocate(operators(maxnumOperators))
        initializedOperators = .true.
         
        operators(:)%nobslev = 0

        write(*,*) 'oopc_operatorDepot: Max number of operators to save: ', &
	  maxnumOperators
         
      end if

      if (allocated(operators(headerCount)%zh)) then 
        deallocate(operators(headerCount)%trial)
        deallocate(operators(headerCount)%vlayertop)
        deallocate(operators(headerCount)%vlayertop)
        deallocate(operators(headerCount)%modlevindexTop)
        deallocate(operators(headerCount)%modlevindexBot)
        deallocate(operators(headerCount)%zh)
        deallocate(operators(headerCount)%ixtr)
        deallocate(operators(headerCount)%success)
        if ( obsoper%applyGenOper ) deallocate(operators(headerCount)%zhp)
      end if
      
      allocate(operators(headerCount)%trial(obsoper%nmodlev))    
      allocate(operators(headerCount)%vlayertop(obsoper%nobslev))    
      allocate(operators(headerCount)%vlayerbottom(obsoper%nobslev))
      allocate(operators(headerCount)%modlevindexTop(obsoper%nobslev)) 
      allocate(operators(headerCount)%modlevindexBot(obsoper%nobslev)) 
      allocate(operators(headerCount)%zh(obsoper%nobslev,obsoper%nmodlev)) 
      allocate(operators(headerCount)%success(obsoper%nobslev))    
      allocate(operators(headerCount)%ixtr(obsoper%nobslev))    
      
      operators(headerCount)%nobslev=obsoper%nobslev
      operators(headerCount)%iavgkern=obsoper%iavgkern
      operators(headerCount)%applyGenOper=obsoper%applyGenOper
      operators(headerCount)%trial(:) = obsoper%trial(:)
      operators(headerCount)%vlayertop(:) = obsoper%vlayertop(:)
      operators(headerCount)%vlayerbottom(:) = obsoper%vlayerbottom(:)
      operators(headerCount)%modlevindexTop(:) = obsoper%modlevindexTop(:)
      operators(headerCount)%modlevindexBot(:) = obsoper%modlevindexBot(:)
      operators(headerCount)%zh(:,:) = obsoper%zh(:,:)
      operators(headerCount)%ixtr(:) = obsoper%ixtr(:)
      operators(headerCount)%success(:) = obsoper%success(:)
 
      if ( obsoper%applyGenOper ) then
        allocate(operators(headerCount)%zhp(obsoper%nobslev,obsoper%nmodlev)) 
        operators(headerCount)%zhp(:,:) = obsoper%zhp(:,:)
      end if
      
      success = .true. 

    end if
    
  end function oopc_operatorDepot
  
  !--------------------------------------------------------------------------
  ! oopc_convertUnits
  !--------------------------------------------------------------------------
  subroutine oopc_convertUnits(model_col,ppb_opt,incr_opt)
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
    !      B. Unit-conversion factor is calculated in oopc_convertUnits
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
    real(8), intent(inout) :: model_col(obsoper%nmodlev) ! Model-space profile to have its units changed
    logical, intent(in), optional :: ppb_opt, incr_opt
    
    ! Locals:
    real(8) :: zcoef
    integer :: exp_P,exp_T
    logical :: ppb_out, incr_out
    real(8), parameter :: rho_stp=1.293  ! kg/m^3
    
    if (obsoper%constituentId < 0) return
    
    ! No conversion necessary for these BUFR numbers
    if (any( obsoper%varno == (/ BUFR_UNIT_OptDepth,BUFR_UNIT_OptDepth2, &
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

    if (obsoper%varName(1:2) == 'AF' .or. obsoper%varName(1:2) == 'AC') then

      ! PM2.5 or PM10
       
      if (any(obsoper%varno == (/ BUFR_UNIT_VMR, BUFR_UNIT_VMR2,  &
              BUFR_UNIT_MolePerMole, BUFR_UNIT_MolePerMole2,   &
              BUFR_UNIT_NumberDensity, BUFR_UNIT_MolarDensity, & 
              BUFR_UNIT_PartPress, BUFR_UNIT_PartPress2, &
              BUFR_UNIT_DU, BUFR_UNIT_DU2, BUFR_UNIT_DU3, BUFR_UNIT_DU4, &
              BUFR_UNIT_IntegND, BUFR_UNIT_IntegND2 /) )) then
		    
        call utl_abort("oopc_convertUnits: BUFR # " // trim(utl_str(obsoper%varno)) // " is not valid for PM" )
      end if
        
      ! Conversion from ug/m^3 to ug/kg  (scaling by Rd*T/P)
      zcoef = zcoef * MPC_RGAS_DRY_AIR_R8
      exp_T = exp_T+1 ! multiply by T
      exp_P = exp_P-1 ! divide by P
   
    else if (obsoper%varName(1:2) == 'HU') then
       
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
         
        zcoef = zcoef * 1.0d-6 / (ec_rg*vnl_varMassFromVarNum(obsoper%constituentId))
         
      case(BUFR_UNIT_IntegND, BUFR_UNIT_IntegND2)
         
        ! For conversion from ug/kg to integrand values in molecules/(m^2*Pa)
        ! Note: 1 molecule/(m^2*Pa) = 1E-9 ug/kg * Na / [RG * (1E-3 kg/g * m_gas)]
        ! 
        ! To convert from kg/m^2 for the gas to molecules/m^2, one must 
        ! divide by the gas molar mass (kg/mole) and multiply by the Avogrado number          

        zcoef = zcoef * 1.0d-6 * MPC_AVOGADRO_R8 &
                / (ec_rg*vnl_varMassFromVarNum(obsoper%constituentId))
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
                /(vnl_varMassFromVarNum(obsoper%constituentId)*ec_rg*rho_stp)
        
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
                /vnl_varMassFromVarNum(obsoper%constituentId)
        exp_P = exp_P+1 ! multiply by P

      case(BUFR_UNIT_NumberDensity)
          
        ! For conversion from ug/kg to molecules/m^3
        !
        ! Number density of gas = Na*rho(gas)/m_gas = Na*rho(air) * mass mixing ratio /m_gas
        !                       = Na * P/Rd/T * mass mixing ratio /m_gas
        
        zcoef = zcoef * 1.0d-6 * MPC_AVOGADRO_R8/MPC_RGAS_DRY_AIR_R8 &
                /vnl_varMassFromVarNum(obsoper%constituentId)
        exp_T = exp_T-1 ! divide by T
        exp_P = exp_P+1 ! multiply by P

      case(BUFR_UNIT_MolarDensity)
          
        ! For conversion from ug/kg to moles/m^3
        !
        ! Mole density of gas = rho(gas)/m_gas = rho(air) * mass mixing ratio /m_gas
        !                       = P/Rd/T * mass mixing ratio /m_gas

        zcoef = zcoef * 1.0d-6 /MPC_RGAS_DRY_AIR_R8 &
                /vnl_varMassFromVarNum(obsoper%constituentId)
        exp_T = exp_T-1 ! divide by T
        exp_P = exp_P+1 ! multiply by P

      case(BUFR_UNIT_VMR, BUFR_UNIT_VMR2, BUFR_UNIT_MolePerMole, &
           BUFR_UNIT_MolePerMole2)
          
        ! For conversion from ug/kg to vmr (or moles/mole)
          
        zcoef = zcoef * 1.0d-9 * MPC_MOLAR_MASS_DRY_AIR_R8 &
                /vnl_varMassFromVarNum(obsoper%constituentId)

      case default 
        
        call utl_abort('oopc_convertUnits: Unknown obs units ' // &
	               'for varno = ' //  trim(utl_str(obsoper%varno)) )
         
      end select
    end if
  
    ! Apply constant scaling
    model_col = model_col * zcoef
    
    if (exp_T /= 0) then
      if (any(obsoper%tt <= 0.)) then
        call utl_abort("oopc_convertUnits: " // &
                       "Missing valid temperature for conversion.")
      end if	  
      model_col = model_col * obsoper%tt**exp_T
    end if
    
    if (exp_P /= 0) model_col = model_col * obsoper%pp**exp_P
    
  end subroutine oopc_convertUnits

  !--------------------------------------------------------------------------
  ! oopc_required_field
  !--------------------------------------------------------------------------
  function oopc_required_field(varName,varno) result(needed)
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

  end function oopc_required_field

  !--------------------------------------------------------------------------
  ! oopc_vertObsLayersWgts
  !--------------------------------------------------------------------------
  subroutine oopc_vertObsLayersWgts(operator,ixtr,success,kmode,skipType)
    !
    !:Purpose: To calculate integration (or averaging) weights "wgts" required for vertical 
    !          integration (or averaging) for the full vertical range or a set of target layers. 
    !          Given the calculated weights and a user input array vector X, the integral 
    !          (or average) for a given layer i would be given by sum(wgts(i,:)*X(:))
    !
    !:Arguments:
    !     :operator:    Operator type: 'integ' or 'avg'
    !     :ixtr:        Flag indicating if obs outside model vertical range: 0 for no.
    !     :successs:    Success of integration (or averaging)
    !     :kmode:       Observation model stage used to allow option of tropo
    !                   increment determination from total (or avg) column data when
    !                   obsoper%columnBound > obsoper%vlayertop for
    !                   nobslev=1. For kmode=3, calc only for region between
    !                   obsoper%columnBound and surface. For kmode=2, split
    !                   calc for model top to obsoper%columnBound and 
    !                   obsoper%columnBound and surface.
    !
    !     :skipType:    Skipping processing of specific target layers depending on case:
    !                   'default' - skipping application via input success_opt only
    !                   'doAll&noExtrap' - application of both success_opt and outbound_opt
    !
    
    implicit none

    ! Arguments:
    integer, intent(inout) :: ixtr(obsoper%nobslev) ! Flag indicating if obs outside model vertical range: 0 for no.
    logical, intent(inout) :: success(obsoper%nobslev) ! success of integration (or averaging)
    integer, intent(in) :: kmode ! Mode of observation operator 
    character(len=*), intent(in) :: operator ! Operator type: 'integ' or 'avg'
    character(len=*), intent(in) :: skipType ! Skipping processing of specific target layers depending on case

    ! Locals:
    integer :: obslevIndex,tropo_mode
    real(8) :: vlayertop_ref,vlayerbottom_ref,modlevindexBot_ref,checkID
    
    ! Conduct initial setup for vertical integration (or avegaging) components

    call ppo_vertLayersSetup(operator,obsoper%pp,obsoper%nmodlev)

    ! Ensure that each layer is within model vertical range.
    
    do obslevIndex=1,obsoper%nobslev
    
      if (obsoper%vlayerbottom(obslevIndex) < obsoper%vlayertop(obslevIndex)) then
        success(1:obsoper%nobslev)=.false.
        write(*,*) 'oopc_vertObsLayersWgts: WARNING. ' // &
	          'Layer top/bot value problem.', &
                   obsoper%vlayertop(obslevIndex), &
		   obsoper%vlayerbottom(obslevIndex), &
                   '. Entire profile skipped over.'
        return
      else if (obsoper%vlayerbottom(obslevIndex) < obsoper%pp(1)*1.01 .or. &  
               obsoper%vlayertop(obslevIndex) > obsoper%pp(obsoper%nmodlev)*0.99) then

        success(obslevIndex)=.false.
        if (obsoper%vlayerbottom(obslevIndex) < obsoper%pp(1)*1.01) then
           ixtr(obslevIndex)=1
        else
           ixtr(obslevIndex)=2
        end if
        write(*,*) 'oopc_vertObsLayersWgts: WARNING. Layer top/bot value problem.', &
             obsoper%vlayertop(obslevIndex), obsoper%vlayerbottom(obslevIndex)
        cycle
      end if
      if (obsoper%vlayerbottom(obslevIndex) > &
          obsoper%pp(obsoper%nmodlev)*0.999) then
	   
	 obsoper%vlayerbottom(obslevIndex)=obsoper%pp(obsoper%nmodlev)*0.999

      end if
      if (obsoper%vlayertop(obslevIndex) < obsoper%pp(1)*1.001) then
        obsoper%vlayertop(obslevIndex)=obsoper%pp(1)*1.001
      end if
      
    end do

    tropo_mode=0
    
    if (obsoper%nobslev == 1 .and. kmode >= 2 .and. trim(operator) == 'integ' &
        .and. obsoper%vlayerbottom(1) > obsoper%pp(obsoper%nmodlev)*0.99 .and.  &
        obsoper%constituentId >= 0 .and. obsoper%vco == 4) then
       
      ! Check if constituent id is recognized (function will abort if not recognized)
      if ( obsoper%constituentId >= 0 ) then
        checkID = vnl_varMassFromVarnum(obsoper%constituentId)
      end if
       
      vlayerbottom_ref=obsoper%vlayerbottom(1)

      ! Check for special treatment if tropo_mode>=1, kmode=2,3, and nobslev=1 for
      ! column observations (obsoper%vco=4)  that extend to the surface.

      tropo_mode = oopc_tropo_mode(obsoper%constituentId)
       
      if ( tropo_mode >= 1 .and. obsoper%columnBound > obsoper%vlayertop(1) ) then
          
        if (obsoper%iavgkern /= 0) then
          call utl_abort("oopc_vertObsLayersWgts: Use of averaging ' // & 
	       'kernels not possible with reduced range of increment profile.")
	end if
          
        if (kmode == 2 .and. tropo_mode == 1) then
             
          ! When kmode=2, split calc in two. This is done due to difference in 
          ! calc at the interface region when producing zh and zhp. The tangent linear
          ! model in the lower region for kmode=2 must be consistent with 
          ! that associated with kmode=3.
             
          ! Start with bottom region in order to use correct zhp with oopc_genOper
          ! when use of this operator is requested.
           
          vlayertop_ref=obsoper%vlayertop(1)
             
          obsoper%vlayertop(1) = obsoper%columnBound
             
          call ppo_vertIntegWgts(obsoper%vlayertop,obsoper%vlayerbottom, &
	       obsoper%nmodlev,obsoper%nobslev,obsoper%modlevindexTop, &
               obsoper%modlevindexBot,obsoper%zh,wgts_opt=obsoper%zhp, &
               skipType_opt=skipType,outbound_opt=ixtr,success_opt=success)
	       	 
          ! Apply generalized innovation operator if requested             
          if (obsoper%applyGenOper) call oopc_genoper(kmode)
            
          obsoper%applyGenOper=.false.
             
          modlevindexBot_ref=obsoper%modlevindexBot(1)
             
          ! Reset top and bottom values for integration (or averaging) of the remaining region.
          ! The second integration (or averaging) provides the change in upper level contributions to the
          ! total column (or average) from the assimilation of other observations.
             
          obsoper%vlayertop(1)=vlayertop_ref
          obsoper%vlayerbottom(1)=obsoper%columnBound
        else
          ! Reset top new value. Restricts adjoint/tangent linear calcs to this reduced region.            
          obsoper%vlayertop(1)=obsoper%columnBound
        end if
          
      end if
    else
      vlayerbottom_ref=0.0
    end if
    
    ! Calculate vertical integration (or averaging) components for specified obs layer(s).

    if ( trim(operator) == 'integ' ) then 
           
      call ppo_vertIntegWgts(obsoper%vlayertop,obsoper%vlayerbottom, &
           obsoper%nmodlev,obsoper%nobslev,obsoper%modlevindexTop, &
	   obsoper%modlevindexBot,obsoper%zh,wgts_opt=obsoper%zhp, &
           skipType_opt=skipType,outbound_opt=ixtr,success_opt=success, &
	   dealloc_opt=.true.)
			      
      ! If tropo_mode=1, reset original vertical range 
      ! for the tangent linear operator 
      if (obsoper%nobslev == 1 .and. kmode == 2 .and. &
        obsoper%constituentId >= 0 .and. obsoper%vco == 4) then
	
        if (tropo_mode == 1 .and.   &
          obsoper%columnBound > obsoper%vlayertop(1) .and.  &
          vlayerbottom_ref > obsoper%pp(obsoper%nmodlev)*0.99) then    
          
          obsoper%vlayerbottom(1)=vlayerbottom_ref
          obsoper%modlevindexBot(1)=modlevindexBot_ref
        end if
      end if
    else
    
      call ppo_vertAvgWgts(obsoper%vlayertop,obsoper%vlayerbottom, &
           obsoper%nmodlev,obsoper%nobslev,obsoper%modlevindexTop, &
	   obsoper%modlevindexBot,obsoper%zh,wgts_opt=obsoper%zhp, &
           skipType_opt=skipType,outbound_opt=ixtr,success_opt=success, &
	   dealloc_opt=.true.)
			   
    end if
        
  end subroutine oopc_vertObsLayersWgts

  !--------------------------------------------------------------------------
  ! oopc_genOper
  !--------------------------------------------------------------------------
  subroutine oopc_genOper(kmode)
    !
    !:Purpose: Set generalized innovation operator for integral or layer avg
    !          obs. Relevant only for incremental fields. This version is
    !          intended to vertically distribute (approximately) the obs increments
    !          proportionally to the (a) background state or (b) the differences
    !          between a reference/climatological state and the background state.
    !          Option (b) is recommended for a spinup period when the shape of the
    !          background state is not physically representative. This is to force
    !          the analysis shape toward the reference/climatology shape in the absence
    !          of local profiles observations.
    !
    !:Input:
    !       :kmode:              index indicating if the operator is to be applied             
    !       :obsoper%zh,zhp:     see routine ppo_vertIntegSetup
    !       :oopc_genOperConstraintType: index specifying the reference state
    !       :oopc_genOperhCorrlenExpnt:  Exponent for horiz. correl. length weighting
    !       :oopc_genOperOmaStatsFactor: Additional OmAStats normalization factor
    !       :bgStats:            structure containing the background stats data
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
    !          where the innovation weight operator 'w' can be set as, 
    !          for the 1D case,
    ! 
    !             w = P[ (h'x)^T ] *  B^{-1}     
    ! 
    !          with  h' is the part of h which excludes resolution dependence
    !          (only/mostly contains the physics part of h; zhp),
    !
    !            - x is the state profile rval
    !            - P is a window cutoff operator (sets small values to zero)
    !            - B is the original/initial total "vertical" covariance matrix
    !              (2D)
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
    !          from the later application of B in grad(Jo). Its specification
    !          is approximate to reduce the effect of numerical error in the 
    !          generating the inverse of B.
    !
    !          A scaling by the 1/bgStats%hcorrlen^q profile is addded to 
    !          partially mitigate effect of the variation in the vertical 
    !          of the influence of neighbouring obs via background horizontal 
    !          error correlations. 
    !   
    !          The specification of the q value is to roughly give (a) a more 
    !          constant-like increment in the absence of the trial profile form 
    !          constraint and (b), in the presence of the form constraint, a 
    !          better match of the max increment level to the max concentration
    !          level. 
    !
    !          As the horizontal density of obs is somewhat variable
    !          and not known locally, the specified q value is not optimal for each
    !          individual case but better than not applying any related mitigation.
    !
    !      (2) The matrix B is the total error covariance matrix (in physical space)
    !          with the related error correlation matrix
    !          *corvert* and std dev provided from *bCovarSetupChem_mod.ftn90*.
    ! 
    !      (3) NOTE: Cases with ensemble-based and or lam-based background
    !          covariances are not taken into account in this version.
    !
    
    implicit none

    ! Arguments:
    integer, intent(in) :: kmode ! Index specifying if content to be applied (i.e. if kmode>1)

    ! Locals:
    real(8), parameter :: pwin=0.01
    integer  :: obslevIndex,modlevIndex,irmse,varIndex
    real(8)  :: zwbw,zhbh,za,work(obsoper%nmodlev),fdeStddev(obsoper%nmodlev,2)
    real(8), parameter :: threshold=1.D-20
    real(8)  :: zmin,rvalw(obsoper%nmodlev),rvalr(obsoper%nmodlev)
    real(8)  :: rvalc(obsoper%nmodlev),rmse
    character(len=22) :: code
   
    if (kmode <= 1) return

    ! Retrieve from stored background error std dev [elemements (:,1-2)] at obs location [and inverses at elements (:,3-4)]    
    fdeStddev = bcsc_retrieveBgStddev(obsoper%nmodlev,2,obsoper%obs_index) 
    if (obsoper%constituentId == 0) then
      ! To mitigate combined effects of approximations in expressions below and 
      ! the forced strong reduction of fdeStddev(:,1) in the top level(s) of the
      ! model. Bkgd ozone error std dev was set to near zero at the lid. The
      ! absence of a resettting such as below would otherwise
      ! lead to oscillatory numerical error effects in the first few levels.
      ! This change will still force a near-zero (small) increment at the lid.
      fdeStddev(1,2) = fdeStddev(2,2)
    end if
    
    ! Identify variable position index in background error correlation matrices
   
    varIndex=1
    do while (trim(bgStats%varNameList(varIndex)) /= '') 
      if (trim(bgStats%varNameList(varIndex)) == trim(obsoper%varName)) exit
      varIndex=varIndex+1
    end do

    if (trim(bgStats%varNameList(varIndex)) == '') then
      call utl_abort('oopc_genOper: Background stats not found for ' // &
                     trim(obsoper%varName) )
    end if
    
    ! Initialize reference mass (mixing ratio) weighting profile as profile from trial field
     
    rvalr(1:obsoper%nmodlev)=obsoper%trial(1:obsoper%nmodlev)
    
    ! Check on rvalr
            
    if (any(abs(rvalr(:)) <= threshold)) then
      if (all(abs(rvalr(:)) <= threshold)) then
        rvalr(:)=1.0 
      else  
        zmin=minval(abs(rvalr(:)), &
             mask=abs(rvalr(:)) > threshold)
        where (abs(rvalr(:)) <= threshold) rvalr(:)=zmin  
      end if
    end if
    
    ! Loop over obs elements
    
    write(code,'(I22)') obsoper%obs_index

    do obslevIndex=1,obsoper%nobslev
       
      if (.not.obsoper%success(obslevIndex)) cycle   
       
      ! Set reference mass (mixing ratio) weighting profile
       
      rvalw(1:obsoper%nmodlev)=rvalr(1:obsoper%nmodlev)
      if (trim(oopc_genOperConstraintType(obsoper%constituentId)) == 'Diff') then

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
        
        rmse=0.0d0
        irmse=0
        rvalc(1:obsoper%nmodlev) = oopc_getProfile(oopc_bgRef,code)
            
        zmin=pwin*maxval(abs(obsoper%zhp(obslevIndex,1:obsoper%nmodlev)))
        do modlevIndex=1,obsoper%nmodlev
          if (obsoper%zhp(obslevIndex,modlevIndex) > zmin) then
            irmse=irmse+1
            rmse=rmse+((rvalc(modlevIndex)-obsoper%trial(modlevIndex))* &
	         fdeStddev(modlevIndex,2))**2
          end if
        end do
        if (irmse > 0) rmse=rmse*2.0d0/irmse
        if (rmse < 1.0d0) then
          rvalw(1:obsoper%nmodlev)=rvalc(1:obsoper%nmodlev)
        else 
          rvalw(1:obsoper%nmodlev)=rvalc(1:obsoper%nmodlev)- &
                                  obsoper%trial(1:obsoper%nmodlev)
          where (abs(rvalw(1:obsoper%nmodlev)) < 0.01d0* &
	         rvalc(1:obsoper%nmodlev)) 
            rvalw(1:obsoper%nmodlev)=sign(0.01d0*rvalc(1:obsoper%nmodlev), &
                                   rvalw(1:obsoper%nmodlev))
	  end where
        end if
        if (all(abs(rvalw(:)) <= threshold)) then
          rvalw(1:obsoper%nmodlev)=rvalc(1:obsoper%nmodlev)
        else if (any(abs(rvalw(:)) <= threshold)) then 
          zmin=minval(abs(rvalw(:)),mask=abs(rvalw(:)) > threshold)
          where (abs(rvalw(:)) <= threshold) rvalw(:)= zmin
	end if  
      end if

      ! Begin preparation of the new innovation operator w (=new zhp)
       
      if (obsoper%nobslev == 1 .and. obsoper%modlevindexTop(obslevIndex) == 1 .and.  &
          obsoper%modlevindexBot(obslevIndex) == obsoper%nmodlev) then
          
        ! Treat as total column obs. Here, zhp would otherwise be approx. equal
        ! to 1 except for the near-end points of the model vertical domain,
        ! the latter due to the discretized domain. Not using zhp avoids this
        ! discretization issue from weakly affecting results at the boundaries.

	work(1:obsoper%nmodlev)=1.0

      else 

        ! Account for localized obs function (e.g. partial columns, Jacobians. For Jacobians,
        ! zhp must also be independent of the model layer thicknesses.)
	
	! Apply cutoff       
        zmin=pwin*maxval(abs(obsoper%zhp(obslevIndex,1:obsoper%nmodlev)))					     
        where (abs(obsoper%zhp(obslevIndex,1:obsoper%nmodlev)) < zmin) 
          work(1:obsoper%nmodlev)=0.0d0 
        elsewhere
          work(1:obsoper%nmodlev)=obsoper%zhp(obslevIndex,1:obsoper%nmodlev)
	endwhere 
      end if

      ! Application of 1D space vertical covariance inverse B^{-1} for partially
      ! mitigate the weight impact of the later application of B in finalizing 
      ! grad(Jo). Note: fdeStddev(:,2)=1.0/fdeStddev(:,1)

      ! Also impose 
      ! - approximate mitigation of horizontal correlation effect
      ! - shape forcing profile via rvalw   

      work(1:obsoper%nmodlev)=work(1:obsoper%nmodlev)*rvalw(1:obsoper%nmodlev) &
        *fdeStddev(1:obsoper%nmodlev,2) 
	
      obsoper%zhp(obslevIndex,:)=0.0D0
      !$OMP PARALLEL DO PRIVATE(modlevIndex)       
      do modlevIndex=obsoper%modlevindexTop(obslevIndex), &
                     obsoper%modlevindexBot(obslevIndex)
        obsoper%zhp(obslevIndex,modlevIndex) = fdeStddev(modlevIndex,2)    &
	     /bgStats%hcorrlen(modlevIndex,varIndex)**oopc_genOperHCorrlenExpnt(varIndex) &
	     *sum(work(1:obsoper%nmodlev)                                  &
	         *bgStats%corverti(1:obsoper%nmodlev,modlevIndex,varIndex))
      end do
      !$OMP END PARALLEL DO

      ! Determine proportionality factor 'a' = (h*B*h^T)(w*B*w^T)^{-1}
       
      ! First determine/estimate w*B*w^T (zwbw)
       
      !$OMP PARALLEL DO PRIVATE(modlevIndex)       
      do modlevIndex=obsoper%modlevindexTop(obslevIndex), &
                     obsoper%modlevindexBot(obslevIndex)
	 
        work(modlevIndex)=sum(obsoper%zhp(obslevIndex, &
	  obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex)) &
           *bgStats%corvert(modlevIndex,obsoper%modlevindexTop(obslevIndex): &
	     obsoper%modlevindexBot(obslevIndex),varIndex) &
           *fdeStddev(obsoper%modlevindexTop(obslevIndex): &
	     obsoper%modlevindexBot(obslevIndex),1))*fdeStddev(modlevIndex,1)
      end do
      !$OMP END PARALLEL DO
       
      zwbw=sum(obsoper%zhp(obslevIndex,obsoper%modlevindexTop(obslevIndex): &
                           obsoper%modlevindexBot(obslevIndex)) &
               *work(obsoper%modlevindexTop(obslevIndex):           &
	             obsoper%modlevindexBot(obslevIndex)))
       
      ! Determine/estimate h*B*h^T (zhbh)

      !$OMP PARALLEL DO PRIVATE(modlevIndex)       
      do modlevIndex=obsoper%modlevindexTop(obslevIndex), &
                     obsoper%modlevindexBot(obslevIndex)
	 
        work(modlevIndex)=sum(obsoper%zh(obslevIndex, &
	  obsoper%modlevindexTop(obslevIndex):obsoper%modlevindexBot(obslevIndex)) &
            *bgStats%corvert(modlevIndex,obsoper%modlevindexTop(obslevIndex): &
	      obsoper%modlevindexBot(obslevIndex),varIndex) &
            *fdeStddev(obsoper%modlevindexTop(obslevIndex): &
	      obsoper%modlevindexBot(obslevIndex),1))*fdeStddev(modlevIndex,1)
      end do
      !$OMP END PARALLEL DO
       
      zhbh=sum(obsoper%zh(obslevIndex,obsoper%modlevindexTop(obslevIndex): &
                          obsoper%modlevindexBot(obslevIndex)) &
               *work(obsoper%modlevindexTop(obslevIndex):          &
	             obsoper%modlevindexBot(obslevIndex)))

      ! Set proportionality factor 'a'

      za=sqrt(zhbh/zwbw)*oopc_genOperOmAStatsFactor(varIndex)
         
      ! Set final innovation operator
      
      obsoper%zh(obslevIndex,1:obsoper%nmodlev)= &
                                 obsoper%zhp(obslevIndex,1:obsoper%nmodlev)*za
         
    end do

  end subroutine oopc_genOper

  !--------------------------------------------------------------------------
  ! oopc_getColBoundary
  !--------------------------------------------------------------------------
  function oopc_getColBoundary(iconstituentId,nmodlev,pressmod,tt,height,hu_opt, &
                               uu_opt,vv_opt) result(boundPress)
    !
    !:Purpose: To determine and to store the boundary (e.g. tropopause or PBL)
    !          pressure levels if needed by the observation operators.    
    implicit none
    real(8) :: boundPress ! pressure level of boundary to be imposed

    ! Arguments:
    integer, intent(in) :: iconstituentId  ! BUFR code element of Table 08046 identifying the constituent
    integer, intent(in) :: nmodlev          ! number of model levels for variables other than uu and vv
    real(8), intent(in) :: pressmod(nmodlev)! model pressure array
    real(8), intent(in) :: tt(nmodlev)      ! model temperature (Kelvin)
    real(8), intent(in) :: height(nmodlev)  ! height (meters)
    real(8), optional, intent(in) :: hu_opt(nmodlev) ! specific humidity
    real(8), optional, intent(in) :: uu_opt(:) ! model zonal wind component(m/s)
    real(8), optional, intent(in) :: vv_opt(:) ! model meridional wind component (m/s)
   
    ! Locals:
    integer :: tropo_bound
    
    boundPress = -1.0d0
    
    if (oopc_tropo_mode(iconstituentId) == 0) return
   
    tropo_bound = oopc_tropo_bound(iconstituentId)
    
    if (tropo_bound > 0) then
      if (.not.all(tt < 0.) .and. .not.all(height < 0.) ) then
    
        select case(tropo_bound)
        case(1)
    
          ! Get tropopause pressure level
      
          boundPress = phf_get_tropopause(nmodlev,pressmod,tt,height,hu_opt=hu_opt)
    
        case(2)
 
          ! Get PBL pressure level
      
          boundPress = phf_get_pbl(nmodlev,pressmod,tt,height,hu_opt=hu_opt, &
	               uu_opt=uu_opt,vv_opt=vv_opt) 
      
        case default
           call utl_abort("oopc_getColBoundary: Unrecognized value for tropo_bound of " &
                          // trim(utl_str(tropo_bound)) )
        end select
                                     
      end if     
    end if
      
    ! Use tropo_column_top value if tropo_bound=0 or model derived boundary was unsuccessful
    if (boundPress < 0.0d0) boundPress = oopc_tropo_column_top(iconstituentId)
      
  end function oopc_getColBoundary
       
  !--------------------------------------------------------------------------
  ! oopc_addColBoundary
  !--------------------------------------------------------------------------
  subroutine oopc_addColBoundary(headerIndex,boundPress)
    !
    !:Purpose: To add column boundary data to oopc_columnBoundary which can be
    !          retrieved later using a header index.
    !
    implicit none 

    ! Arguments:
    integer, intent(in) :: headerIndex ! Header index
    real(8), intent(in) :: boundPress  ! Pressure boundary
     
    if (.not.associated(oopc_columnBoundary%data1d)) then
      call oss_obsdata_alloc(oopc_columnBoundary,oopc_obsdata_maxsize, dim1=1)
      oopc_columnBoundary%nrep = 0
    end if

    ! In this case nrep will count the number of filled reps in the data arrays
    oopc_columnBoundary%nrep = oopc_columnBoundary%nrep+1 

    if (oopc_columnBoundary%nrep > oopc_obsdata_maxsize) then
      call utl_abort('oopc_addColBoundary: Reach max size of array ' // &
	             trim(utl_str(oopc_obsdata_maxsize)) )
    end if
  
    ! Use the header number as the unique code for this obs data
    write(oopc_columnBoundary%code(oopc_columnBoundary%nrep),'(I22)') headerIndex

    oopc_columnBoundary%data1d(1,oopc_columnBoundary%nrep) = boundPress

  end subroutine oopc_addColBoundary

  !--------------------------------------------------------------------------
  ! oopc_retrieveColBoundary
  !--------------------------------------------------------------------------
  function oopc_retrieveColBoundary(headerIndex) result(boundPress)
    !
    !:Purpose: To retrieve previously saved column boundary data in
    !          oopc_columnBoundary from the header index.
    !
    implicit none
    real(8) :: boundPress

    ! Arguments:
    integer, intent(in) :: headerIndex

    ! Locals:
    character(len=22) :: code

    write(code,'(I22)') headerIndex
    
    boundPress = oss_obsdata_get_element(oopc_columnBoundary,code,1)

  end function oopc_retrieveColBoundary

  !==========================================================================
  !-------- Stand-alone routines to read and extract climatology fields -----

  ! Could be placed in a category 4 module.
  ! public :: oopc_readFields, oopc_addToProfileSet, oopc_getProfile

  !--------------------------------------------------------------------------
  ! oopc_readFields
  !--------------------------------------------------------------------------
  subroutine oopc_readFields(climatFields,filename,variable,    &
                             maxNumFields,maxNumTypes,         &
                             fieldRequired,success,filetype_opt)
    !
    !:Purpose:  To read climatrology (reference) fields as directed by input
    !
    ! Comments:
    !      - Fields are provided in RPN/fst files 
    !      - Reference fields can be in a separate RPN file with name provided
    !        within 'filename' if filetype='TXT' or provided as 'filename' if it
    !        refers to an RPN standaard file.
    !      - Fields assumed to be of the same units as those of the
    !        corresponding input trial fields
    !
    implicit none

    ! Arguments:    
    type(struct_oopc_field), intent(out) :: climatFields(0:maxNumFields,maxNumTypes)
    character(len=*), intent(in) :: filename
    integer, intent(in):: maxNumFields,maxNumTypes
    logical, intent(out) :: success
    character(len=*), intent(in) :: variable
    character(len=*), intent(in) :: fieldRequired(0:maxNumFields) 
    character(len=*), intent(in), optional :: filetype_opt

    ! Locals:
    character(len=3) :: filetype
    character(len=256) :: fname
    character(len=4) :: varName
    character(len=12) :: etiket
    integer :: varIndex,id,nd,j,numvar,ijour,imonth,iday,itime,latIndex
    real(8) :: day
    integer :: datestamp
    integer, external :: newdate
   
    integer, external :: fnom, fclos
    integer :: ierr, nulun, ios
    logical :: fileExists
    
    logical :: timeInterp
    
    integer :: ni, nj, nkeys, kind
    real(8), allocatable :: array1(:,:,:),array2(:,:,:),lvls(:),xlat(:),xlong(:) 
    real(8), allocatable :: pressclim(:),ozoneclim(:,:)
    
    character (len=128) :: ligne

    ! Initialize dimensions to zero
    
    climatFields(:,:)%nlon=0
    climatFields(:,:)%nlat=0
    climatFields(:,:)%nlev=1
 
    if ( trim(variable) == 'CH' ) then
      if ( all(fieldRequired(:) == 'Trial') ) then
        ! Not needed
        success=.true.
	return
      end if
    end if
   
    inquire(file=trim(filename),exist=fileExists)
    if ( .not.fileExists ) then
      write(*,*)  '----------------------------------------------------'
      write(*,*)  'oopc_readFields: COULD NOT FIND file ' // trim(filename)
      write(*,*)  '----------------------------------------------------'
      success = .false.
      return
    else
      success = .true.
    end if

    ! Check for file names containing climatological fields or input directives

    if ( present(filetype_opt) ) then
      filetype = trim(fileType_opt) 
    else
      filetype = 'RPN'      
    end if
    
    nulun=0
    ierr=0
    if ( filetype == 'TXT' ) then
      ierr=fnom(nulun,trim(filename),'SEQ',0)
      if ( ierr == 0 ) then
        open(unit=nulun, file=trim(filename), status='OLD')
        ios=0
        
        if ( trim(variable) == 'CH' ) then 
        
          ! CH variable kind (for constituent fields)
          
          read(nulun,'(A)',iostat=ios,err=10,end=10) ligne
          do while (trim(adjustl(ligne(1:14))) /= 'SECTION IV:') 
            read(nulun,'(A)',iostat=ios,err=10,end=11) ligne
          end do    
          
          ! Read number of constituents with associated input file(s)
   
          read(nulun,*,iostat=ios,err=10,end=10) numvar
          if (numvar <= 0) go to 10
        else
          numvar=1
          nd=1
          call utl_abort('oopc_readFields: Variable kind or name ' // &
	                 trim(variable) // ' not taken into account')  
        end if
      end if
    else if ( filetype == 'RPN' ) then
      numvar=1
      nd=1
    else if ( filetype /= 'RPN' ) then
      call utl_abort('oopc_readFields: File type ' // trim(filetype) // &
                     ' not recognized') 
    end if

    if ( ierr /= 0 ) then
      call utl_abort('oopc_readFields: COULD NOT OPEN file ' // trim(filename))
    end if
       
    ! Initialization

    timeInterp = .true.
    datestamp=tim_getDateStamp()
    ierr = newdate(datestamp,ijour,itime,-3)
    if ( ierr < 0 ) then
      call utl_abort('oopc_readFields: Invalid datestamp ' // &
                     trim(utl_str(datestamp)) )
    end if
    imonth = MOD(ijour/100,100)
    iday = MOD(ijour,100)
    day=iday+itime*1.0D-8
    if (day > 15.) then
      day=day-15.0
    else
      day=day+15.0
    end if
    
    ! Get needed fields for each file/varIndex

    do varIndex=1,numvar

      if ( trim(variable) == 'CH' ) then 
       
        ! Read id,nd
        ! id: constituent identifier code; (0 for ozone, ...)
        ! nd: number of sets; 1 or 2 (nd=2 required when different profile 
        !       sets need to be merged according to the tropopause height 
	!       when the first set referring to strato files and teh second 
	!       to tropo fields)
       
        read(nulun,*,iostat=ios,err=10,end=10)
        read(nulun,*,iostat=ios,err=10,end=10) id,nd   
        varName=vnl_varnameFromVarnum(0,id)

        read(nulun,*,iostat=ios,err=10,end=10) fname
        inquire(file=trim(fname),exist=fileExists)
        if ( .not. fileExists ) then
          call utl_abort('oopc_readFields: Did not find file ' // trim(fname))
        end if
      else
        id=varIndex
        ! Currently assumes nunmar = 1 and fname = filename. Could be extended
        fname = filename
        varname = trim(variable)
      end if
              
      do j=1,nd
       
        if ( trim(fname) ==  'ozoneclim98' ) then
          timeInterp = .false.
          call ozo_read_Climatology(datestamp,nlat_opt=nj,nlev_opt=nkeys, &
	                            press_opt=pressclim,ozone_opt=ozoneclim) 
          id=0
          ni=1
          allocate(array1(1,nj,nkeys),lvls(nkeys),xlat(nj),xlong(1))
          ! Convert from ppmv to microgram/kg
          array1(1,1:nj,1:nkeys) =  ozoneclim(1:nj,1:nkeys) * &
	         MPC_MOLAR_MASS_O3_R8 / (1.0d-3 * MPC_MOLAR_MASS_DRY_AIR_R8)
          lvls(1:nkeys) = pressclim(1:nkeys)
	  deallocate(ozoneclim,pressclim)
          kind = 2
          xlong(1)=0.0d0
          do latIndex = 1, nj
            xlat(latIndex) = (latIndex-1)*180.0d0/(nj-1) - 90.0d0
          end do
          etiket = '            '       
        else   
          if ( nd == 2 ) then
            read(nulun,*,iostat=ios,err=10,end=10) etiket    
          else
            etiket = '            '       
          end if                   
          call utl_readFstField(trim(fname),varName,-1,imonth,-1,etiket, &
	       ni,nj,nkeys,array1,xlat_opt=xlat,xlong_opt=xlong,         &
               lvls_opt=lvls,kind_opt=kind)
        end if      

        climatFields(id,j)%nlon=ni
        climatFields(id,j)%nlat=nj
        climatFields(id,j)%nlev=nkeys
        climatFields(id,j)%ivkind=kind   
                         
        allocate(climatFields(id,j)%field(ni,nj,nkeys))
        allocate(climatFields(id,j)%vlev(nkeys),climatFields(id,j)%lon(ni))
        allocate(climatFields(id,j)%lat(nj))
              
        climatFields(id,j)%lat(1:nj)=xlat(1:nj)*MPC_RADIANS_PER_DEGREE_R8
        climatFields(id,j)%lon(1:ni)=xlong(1:ni)*MPC_RADIANS_PER_DEGREE_R8
        where (climatFields(id,j)%lon(1:ni) < 0.0) 
	  climatFields(id,j)%lon(1:ni)=2.0*MPC_PI_R8 + climatFields(id,j)%lon(1:ni)
	end where
        climatFields(id,j)%vlev(1:nkeys)=lvls(1:nkeys)              

        if (.not.timeInterp) then

          climatFields(id,j)%field(:,:,:) = array1(:,:,:)

        else

          ! Following for interpolation as a function of days from mid-months.
             
          if (iday > 15) then
            if (imonth == 12) then
              call utl_readFstField(trim(fname),varName,-1,1,-1,etiket, &
                   ni,nj,nkeys,array2,lvls_opt=lvls,kind_opt=kind)
            else
              call utl_readFstField(trim(fname),varName,-1,imonth+1,-1, &
	           etiket,ni,nj,nkeys,array2,lvls_opt=lvls,kind_opt=kind)
            end if
          
            ! Linearly interpolate in time 
	    ! (approximately - assumes 30 day months)

            climatFields(id,j)%field(:,:,:) = (array1(:,:,:)*(30.0-day)+array2(:,:,:)*day)/30.0
             
          else if (iday <= 15) then
            if (imonth == 1) then
              call utl_readFstField(trim(fname),varName,-1,12,-1,etiket, &
                   ni,nj,nkeys,array2,lvls_opt=lvls,kind_opt=kind)
            else
              call utl_readFstField(trim(fname),varName,-1,imonth-1,-1, &
                   etiket,ni,nj,nkeys,array2,lvls_opt=lvls,kind_opt=kind)
            end if

            ! Linearly interpolate in time 
            ! (approximately - assumes 30 day months)

            climatFields(id,j)%field(:,:,:) = (array2(:,:,:)* &
		                             (30.0-day)+array1(:,:,:)*day)/30.0             
          end if          
        end if
 
        if (allocated(array1)) deallocate(array1,lvls,xlat,xlong)
        if (allocated(array2)) deallocate(array2)   
                 
      end do
    end do 
     
 10 if (ios > 0) then
      call utl_abort('oopc_readFields: READING PROBLEM.' // &
                     ' File read error message number: ' // trim(utl_str(ios))) 
    end if   
    close(unit=nulun)
    ierr = fclos(nulun)
    if ( any(fieldRequired(:) == 'Diff') .and. trim(variable) == 'CH' ) then
      do j=0,maxNumFields
        if ( climatFields(j,1)%nlon == 0 .and. trim(fieldRequired(j)) == 'Diff' ) then
          call utl_abort('oopc_readFields: READING PROBLEM. Did not' // &
	                 ' find SECTION IV required for constituent ID ' // &
			 trim(utl_str(j)))
	end if
      end do
    end if 
	
    return    

 11 close(unit=nulun)
    ierr = fclos(nulun)
    if ( any(fieldRequired(:) == 'Diff') .and. trim(variable) == 'CH' ) then
      call utl_abort('oopc_readFields: READING PROBLEM. Did not find ' // &
                      'SECTION IV.') 
    end if 
	     
  end subroutine oopc_readFields

  !--------------------------------------------------------------------------
  ! oopc_addToProfileSet
  !--------------------------------------------------------------------------
  subroutine oopc_addToProfileSet(climatFields,climatProfileSet,maxNumFields,maxNumTypes, &
                                  numModelLevs,modelPressLevs,modelHeightLevs,obsLat, &
                                  obsLong,obsIndex,maxsize,varKind_opt,varNumber_opt,tt_opt,hu_opt)
    
    !:Purpose: To determine and to store a profile at obs location as part of a cumulative
    !          profile set for a specific variable
    !
    !:Input:
    !
    !    :climatFields:           Input fields from which interpolations are done
    !    :climatProfileSet:       Input profile set
    !    :maxNumFields:           Size of first dimension for climatFields
    !    :maxNumTypes:            Size of second dimension for climatFields
    !    :numModelLevs:           Number of model levels
    !    :modelPressLevs          Model pressure array (Pa)
    !    :modelHeightLevs:        Model height (m)
    !    :obsLat:                 Latitude (rad)
    !    :obsLong:                Longitude (rad)
    !    :obsIndex:               Unique measurement identifier    
    !    :varKind_opt:            variable kind (currently only relevant for 'CH')
    !    :varNumber_opt:          Constituent id
    !    :tt_opt:                 Model temperature (Kelvin)
    !    :hu_opt:                 Specific humidity 
    !    :maxsize:                Max number of obs for which climatProfileSet will be used
    !
    !:Output:
    ! 
    !    :climatProfileSet:       Updated profile set (with one profile added for (obs_long,obs_lat))
    !
    !:Comments:
    !
    implicit none

    ! Arguments
    type(struct_oopc_field), intent(in) :: climatFields(0:maxNumFields,maxNumTypes)
    type(struct_oss_obsdata), intent(inout)  :: climatProfileSet
    integer, intent(in):: maxNumFields,maxNumTypes
    integer, intent(in) :: obsIndex,numModelLevs,maxsize
    real(8), intent(in) :: modelPressLevs(numModelLevs),modelHeightLevs(numModelLevs)
    real(8), intent(in) :: obsLat,obsLong
    
    integer, intent(in), optional :: varNumber_opt
    real(8), intent(in), optional :: tt_opt(:),hu_opt(:)
    character(len=*), optional :: varKind_opt
    
    ! Locals
    integer :: level,start,id
    real(8) :: tropo_press, refprof(numModelLevs),refprof2(numModelLevs),dt
    real(8), allocatable :: pressrefin(:)
    logical, allocatable :: success(:)

    if ( present(varNumber_opt) ) then
      if ( varNumber_opt  < 0 ) return
      id = varNumber_opt
    else
      id = 0
    end if

    if ( present(varKind_opt) ) then
      ! Not currently used
    end if
    
    if (climatFields(id,1)%nlat == 0) return
    
    ! Set vertical levels of reference.
    ! Convert to pressure coordinate if needed.
    
    if (allocated(pressrefin)) deallocate(pressrefin)
    allocate(pressrefin(climatFields(id,1)%nlev))
    pressrefin(:) = climatFields(id,1)%vlev(1:climatFields(id,1)%nlev)

    if (allocated(success)) deallocate(success)
    allocate(success(climatFields(id,1)%nlev))
    success(:)=.true.
    
    if (climatFields(id,1)%ivkind == 2) then
      pressrefin(:)=pressrefin(:)*100. ! Conversion from hPa to Pa.
    else if (climatFields(id,1)%ivkind == 0) then
      where (pressrefin < modelHeightLevs(numModelLevs))
        pressrefin=modelHeightLevs(numModelLevs)
      end where
      pressrefin(:) = phf_convert_z_to_pressure(pressrefin,modelHeightLevs,modelPressLevs, &
                      climatFields(id,1)%nlev,numModelLevs,obsLat,success)
    else if (climatFields(id,1)%ivkind == 4) then
      pressrefin(:)=pressrefin(:) + modelHeightLevs(numModelLevs)
      pressrefin(:) = phf_convert_z_to_pressure(pressrefin,modelHeightLevs,modelPressLevs, &
                      climatFields(id,1)%nlev,numModelLevs,obsLat,success)
    else if (climatFields(id,1)%ivkind == 1) then
      pressrefin(:)=pressrefin(:)*modelPressLevs(numModelLevs) ! Convert from sigma to Pa   
    else
       call utl_abort('oopc_addToProfileSet: Cannot handle vertical coordinate of kind ' // trim(utl_str(climatFields(id,1)%ivkind)))
    end if
    
    ! Interpolate to obs lat/long (or lat) location and model level

    call oopc_column_hbilin(climatFields(id,1)%field,pressrefin, &
                  climatFields(id,1)%nlon,climatFields(id,1)%nlat,climatFields(id,1)%nlev, &
                  climatFields(id,1)%lon,climatFields(id,1)%lat,obsLong,obsLat, &
                  refprof,modelPressLevs,numModelLevs)

    if (climatFields(id,2)%nlat > 0 .and. climatFields(id,2)%nlon > 0 &
        .and. climatFields(id,2)%nlev > 0) then
        
      if ( .not. present(tt_opt) ) then
        call utl_abort('oopc_addToProfileSet: Missing TT for determining ' // &
	               'tropopause pressure')
      end if
      if ( any(tt_opt <= 0.0d0) ) then
        call utl_abort('oopc_addToProfileSet: Invalid TT for determining ' // &
	               'tropopause pressure')
      end if
        
      ! Get second reference field (for troposphere)
        
      tropo_press=-1.0
        
      if ( present(hu_opt) ) then
        if (all(hu_opt >= 0.0D0)) then
          tropo_press=phf_get_tropopause(numModelLevs,modelPressLevs, &
	              tt_opt,modelHeightLevs,hu_opt=hu_opt)
        else
          tropo_press=phf_get_tropopause(numModelLevs,modelPressLevs, &
	              tt_opt,modelHeightLevs)
        end if
      else
        tropo_press=phf_get_tropopause(numModelLevs,modelPressLevs,tt_opt,modelHeightLevs)
      end if
	
      if (tropo_press > 0) then
          
        ! Set vertical levels of reference.
        ! Convert to pressure coordinate if needed
 
        if (allocated(pressrefin)) deallocate(pressrefin)
        allocate(pressrefin(climatFields(id,2)%nlev))    
        pressrefin(:)= climatFields(id,2)%vlev(1:climatFields(id,2)%nlev)

        if (allocated(success)) deallocate(success)
        allocate(success(climatFields(id,2)%nlev))
        success(:)=.true.

        if (climatFields(id,2)%ivkind == 2) then
          pressrefin(:)=pressrefin(:)*100. ! Conversion from hPa to Pa.
        else if (climatFields(id,2)%ivkind == 0) then
          where (pressrefin < modelHeightLevs(numModelLevs)) 
	    pressrefin=modelHeightLevs(numModelLevs)
	  end where 
          pressrefin(:) = phf_convert_z_to_pressure(pressrefin, &
	                  modelHeightLevs,modelPressLevs, &
                          climatFields(id,2)%nlev,numModelLevs, &
	                  obsLat,success)
        else if (climatFields(id,2)%ivkind == 4) then
          pressrefin(:)=pressrefin(:) + modelHeightLevs(numModelLevs)
          pressrefin(:) = phf_convert_z_to_pressure(pressrefin, &
	                   modelHeightLevs,modelPressLevs, &
                           climatFields(id,2)%nlev,numModelLevs,obsLat, &
	                   success)
        else if (climatFields(id,2)%ivkind == 1) then
          pressrefin(:)=pressrefin(:)*modelPressLevs(numModelLevs) ! Convert from sigma to Pa   
        else
          call utl_abort('oopc_addToProfileSet: Cannot handle vertical ' // &
	      'coordinate of kind ' // trim(utl_str(climatFields(id,2)%ivkind)))
        end if
            
        ! Interpolate to obs lat/long (or lat) and model levels
            
        call oopc_column_hbilin(climatFields(id,2)%field,pressrefin, &
             climatFields(id,2)%nlon,climatFields(id,2)%nlat,climatFields(id,2)%nlev, &
             climatFields(id,2)%lon,climatFields(id,2)%lat,obsLong,obsLat, &
             refprof2,modelPressLevs,numModelLevs)
    
      end if

       ! Combine with upper level profile
       
       do level=numModelLevs,3,-1
         if (modelPressLevs(level) < tropo_press) exit
         refprof(level)=refprof2(level)            
       end do
       start=level
            
       ! Apply linear combination of four levels just above the tropopause
        
       do level=start,max(2,start-3),-1
         dt=(start+1.0-level)/5.0
         refprof(level)=dt*refprof2(level) + (1.0-dt)*refprof(level)
      end do
                    
    end if 

    if (allocated(pressrefin)) deallocate(pressrefin)
    if (allocated(success)) deallocate(success) 

    ! ------- Save in climatProfileSet ---------
       
    if (.not.associated(climatProfileSet%data1d)) then
      call oss_obsdata_alloc(climatProfileSet, maxsize, dim1=numModelLevs)
      climatProfileSet%nrep = 0
    end if

    ! Here, nrep will count the number of filled elements in the data arrays
    climatProfileSet%nrep = climatProfileSet%nrep+1 

    if (climatProfileSet%nrep > maxsize) then
      call utl_abort('oopc_addToProfilesSet: Reach max size of array ' // &
	             trim(utl_str(maxsize)) )
    end if
    
    ! obsIndex serves as the unique locator code 
    write(climatProfileSet%code(climatProfileSet%nrep),'(I22)') obsIndex
    
    ! Save profile in climatProfileSet
    
    climatProfileSet%data1d(:,climatProfileSet%nrep) = refprof(:)

  end subroutine oopc_addToProfileSet
  
  !--------------------------------------------------------------------------
  ! oopc_getProfile
  !--------------------------------------------------------------------------
  function oopc_getProfile(climatProfileSet,code) result(profile)
    !
    !:Purpose: To extract and provide profile from climatProfileSet according to 
    !          code value.     
    !  
    implicit none
  
    ! Arguments
    type(struct_oss_obsdata), intent(inout)  :: climatProfileSet  ! Profile set
    character(len=*), intent(in) :: code      ! unique obs identifying code    
    real(8) :: profile(climatProfileSet%dim1) ! retrieved array from obsdata%data1d of dimension obsdata%dim1

    ! Locals:
    integer :: status ! search success (0 = found; 1 = no data; 2 = not found)

    profile = oss_obsdata_get_array1d(climatProfileSet,code,status)
    if (status > 0) then
      call utl_abort("oopc_getProfile: Code not found - " // trim(code))
    end if
    
  end function oopc_getProfile

  !--------------------------------------------------------------------------
  ! oopc_column_hbilin
  !--------------------------------------------------------------------------
  subroutine oopc_column_hbilin(field,vlev,nlong,nlat,nlev,xlong,xlat, &
                               plong,plat,vprof,vlevout,nlevout)
    !
    ! :Purpose: Horizontal bilinear interpolation from a 3D field to a profile at (plong,plat).
    !           Assumes vertical interpolation not needed or already done.
    !
    !           This version can be used with fields that are not part of the background state,
    !           such as climatologies.
    !
    !           This version does not depend in column_data and gridstatevector modules.
    !
    implicit none

    ! arguments:
    integer, intent(in) :: nlong            ! number or longitudes
    integer, intent(in) :: nlat             ! number or latitudes
    integer, intent(in) :: nlev             ! number of vertical levels
    integer, intent(in) :: nlevout          ! number of target vertical levels
    real(8), intent(in) :: field(nlong,nlat,nlev) ! 3D field
    real(8), intent(in) :: vlev(nlev)       ! vertical levels of input field (in pressure)
    real(8), intent(in) :: xlong(nlong)     ! longitudes (radians)
    real(8), intent(in) :: xlat(nlat)       ! latitudes (radians)
    real(8), intent(in) :: plong            ! target longitude (radians)
    real(8), intent(in) :: plat             ! target latitude (radian)
    real(8), intent(in) :: vlevout(nlevout) ! target vertical levels (in pressure)
    real(8), intent(out) :: vprof(nlevout)  ! profile at (plong,plat)
    
    ! locals:
    real(8) :: lnvlev(nlev),lnvlevout(nlevout),plong2
    integer :: ilev,lonIndex,latIndex,i,j

    real(8) :: DLDX, DLDY, DLDP, DLW1, DLW2, DLW3, DLW4

    call utl_tmg_start(30,'--StateToColumn')

    ! Find near lat/long grid points

    if ( nlong > 1 ) then
      plong2 = plong
      if (plong2 < 0.0) plong2 = 2.D0*MPC_PI_R8 + plong2
      do lonIndex = 2, nlong
        if  (xlong(lonIndex-1) < xlong(lonIndex)) then
          if (plong2 >= xlong(lonIndex-1) .and. plong2 <= xlong(lonIndex)) exit
        else 
          ! Assumes this is a transition between 360 to 0 (if it exists). Skip over.
        end if
      end do
      lonIndex = lonIndex-1
    else
      lonIndex=0
    end if

    do latIndex = 2, nlat
      if (plat <= xlat(latIndex)) exit
    end do
    latIndex = latIndex-1

    if ( lonIndex == 0 ) then
    
      ! Set lat interpolation weights

      DLDY = (plat - xlat(latIndex))/(xlat(latIndex+1)-xlat(latIndex))

      DLW1 = (1.d0-DLDY)
      DLW2 = DLDY

      ! Set vertical interpolation weights (assumes pressure vertical coordinate)

      lnvlevout(:) = log(vlevout(:))    
      lnvlev(:) = log(vlev(:))    

      ilev = 1
      do i = 1, nlevout
        do j = ilev, nlev          
          if (lnvlevout(i) < lnvlev(j)) exit ! assumes lnvlevout and lnvlev increase with index
        end do
        ilev = j-1
        if (ilev < 1) then
          ilev = 1
        else if (ilev >= nlev) then
          ilev = nlev-1
        end if

        DLDP = (lnvlev(ilev+1)-lnvlevout(i))/(lnvlev(ilev+1)-lnvlev(ilev))
          
        vprof(i) = DLDP* (DLW1 * field(lonIndex,latIndex,ilev)      &
                        + DLW2 * field(lonIndex,latIndex+1,ilev))   & 
          + (1.d0-DLDP)* (DLW1 * field(lonIndex,latIndex,ilev+1)    &
                        + DLW2 * field(lonIndex,latIndex+1,ilev+1))  
      end do
      
    else
    
      ! Set lat/long interpolation weights

      DLDX = (plong - xlong(lonIndex))/(xlong(lonIndex+1)-xlong(lonIndex))
      DLDY = (plat - xlat(latIndex))/(xlat(latIndex+1)-xlat(latIndex))

      DLW1 = (1.d0-DLDX) * (1.d0-DLDY)
      DLW2 =       DLDX  * (1.d0-DLDY)
      DLW3 = (1.d0-DLDX) *       DLDY
      DLW4 =       DLDX  *       DLDY

      ! Set vertical interpolation weights (assumes pressure vertical coordinate)

      lnvlevout(:) = log(vlevout(:))    
      lnvlev(:) = log(vlev(:))    

      ilev = 1
      do i = 1, nlevout
        do j = ilev, nlev          
          if (lnvlevout(i) < lnvlev(j)) exit ! assumes lnvlevout and lnvlev increase with index
        end do
        ilev = j-1
        if (ilev < 1) then
          ilev = 1
        else if (ilev >= nlev) then
          ilev = nlev-1
        end if

        DLDP = (lnvlev(ilev+1)-lnvlevout(i))/(lnvlev(ilev+1)-lnvlev(ilev))
          
        vprof(i) = DLDP* (DLW1 * field(lonIndex,latIndex,ilev)      &
                        + DLW2 * field(lonIndex+1,latIndex,ilev)    &
                        + DLW3 * field(lonIndex,latIndex+1,ilev)    &
                        + DLW4 * field(lonIndex+1,latIndex+1,ilev)) &
          + (1.d0-DLDP)* (DLW1 * field(lonIndex,latIndex,ilev+1)    &
                        + DLW2 * field(lonIndex+1,latIndex,ilev+1)  &
                        + DLW3 * field(lonIndex,latIndex+1,ilev+1)  &
                        + DLW4 * field(lonIndex+1,latIndex+1,ilev+1))                               
      end do
    end if

    call utl_tmg_stop(30)

  end subroutine oopc_column_hbilin

end module obsOperatorsChem_mod
