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

!-----------------------------------------------------------------------------
!! MODULE chem_setup_mod (prefix='chm')
!!         
!! *Purpose*: Provides pre-processing (setup) routines and 
!!            related tools reqarding observation averaging kernels, 
!!            and std. (fst) gridded format files.
!!  
!!
!! @author Mike Sitwell and Yves Rochon (ARQI/AQRD)
!!
!! Revisions:  
!!   
!! Public routines:
!!v
!!v       - "chm_setup": Routines and structure for setting and assignment 
!!v         of observation layer top and bottom levels (and averaging kernel matrices - tbc). 
!!v         See 'preproc.ftn90' and 'chm_obsoperators'. 
!!v
!!v       - "chm_apply_2dfieldr4_transform" applied to 2D field variable transformation if requested.
!!v
!!v       - "chm_find_avgkern" to find averaging kernel matrix associated to an obs
!!v
!!v       - "chm_get_layer_boundaries" to return layer boundaries for an observation. 
!!v
!!v       - "chm_get_avgkern" to return averaging kernel for an observation.
!!v
!!v       - "chm_var_maxnumber" to pass on value of chm_constituents_size
!!v
!!v       - "chm_add_efftemp_*" adds effective temperatures in obsfile and or 
!!v          obsdata subspace
!!v
!!v        - "chm_set_reference_obsdata" determines and stores reference 
!!v          profile at obs location if needed by the observation operators.
!!v
!!v       - "chm_get_ref_column" to extract and provide column from chm_ref_field 
!!v          associated to code.
!!v
!!v       - "chm_setup_get_*" allows access to desired private parameters, most
!!v          of them provided via the NAMCHEM namelist
!!v
!!v       - "chm_setup_set_int" allows setting of some integer element(s) 
!!v         specific to the CH family outside this module. 
!!v
!!v       - "chm_diagn_only": Identify if data is to be assimilated or only used
!!v         for independent verifications.
!!v
!!
!! Comment:
!!
!! See !x regarding commented output in *2dfieldr4* to ascii file as sample
!! usage of message output file. Commented out as the same file cannot be
!! simultaneoulsly opened for wriring by different processors. Other similar
!! usage in code has been removed.
!!
!----------------------------------------------------------------------------------------
module chem_setup_mod
 
  use mpi_mod, only: mpi_myid, mpi_allgather_string
  use utilities_mod
  use obsSubSpaceData_mod
  use bufr_mod
  use varNameList_mod
  use MathPhysConstants_mod
  use burpFiles_mod
  use obsFiles_mod
  use physicsfunctions_mod
  use stateToColumn_mod, only: s2c_column_hbilin
  
  implicit none
  private
  
! public procedures
! -----------------

  public :: chm_setup, chm_apply_2dfieldr4_transform
  public :: chm_find_avgkern, chm_get_layer_boundaries, chm_get_avgkern
  public :: chm_var_maxnumber, chm_diagn_only
  public :: chm_add_efftemp_obsfile, chm_add_efftemp_obsdata
  public :: chm_get_ref_column, chm_set_reference_obsdata
  public :: chm_setup_get_float, chm_setup_get_str, chm_setup_get_int
  public :: chm_setup_set_int
  
! public types (for use in chm_obsoperators_mod.ftn90)
! ----------------------------------------------------
  public :: struct_chm_obsoperators

! module constants
! -----------------
  integer, parameter :: chm_cfg_size=100           ! max size of config arrays
  integer, parameter :: chm_constituents_size=30   ! max size of constituents arrays
                                                   ! = max allowed value of "iconstituent_id" for Table 08046.
                                                   ! Value to be increased as needed up to a max of 6999 as values
                                                   ! > 7000 (and less 0) are assumed assigned to non-constituent fields  
  integer, parameter :: chm_filename_size=50       ! Max size of *_filename 
  character(len=chm_filename_size), parameter :: chm_aux_filename="obsinfo_chm"  ! auxiliary file name with supplimental
                                                                                 ! observational information

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
     !  gz                     Geopotential height on model levels (m)
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
     !                         indicates index in chm_avgkern,chm_obsSub_avgkern arrays.
     !  apply_genoper          Indicates if the generalized observation operator should be applied
     !  column_bound           Boudary imporsed on a column measurement
     !  dtransform             Derivative for any transform that needs to be applied to a profile
     
     integer :: nobslev,nmodlev,modelIndex,constituent_id,vco,varno,date,hhmm,iavgkern,obs_index
     logical :: layer_identified,apply_genoper
     real(8) :: lat,lon,column_bound
     character(len=12) :: stnid
     character(len=4)  :: varName
     real(8), allocatable :: vlayertop(:),vlayerbottom(:),vmodpress(:),tt(:),gz(:),pp(:)
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

! Array of pointers to averaging kernels read from observation files.
! Note: Ideally, these should be an element in the 'struct_chm_info' derived 
! types, but currently this result in an internal compiler error.

  type(struct_oss_obsdata), allocatable :: chm_obsSub_avgkern(:)

! Arrays containing input reference fields and fields interpolated 
! to obs locations

  type(struct_oss_obsdata)  :: chm_ref_trial
  type(struct_chm_griddata) :: chm_ref_fields(0:chm_constituents_size,2)

! Arrays to contain the calculated concentration-weighted effective temperature
! associated to total column data. It will be stored in the observation file.

  type(struct_oss_obsdata) :: chm_efftemp

! General config/setup information parameters 
! See description list of NAMCHEM namelist parameters in routine chm_setup

  integer, parameter :: chm_famnum=1
  integer :: assim_famNum
  logical :: assim_all(chm_famnum)
  integer :: assim_num(chm_famnum),assim_varno(chm_famnum,chm_cfg_size),assim_nlev(chm_famnum,chm_cfg_size)
  integer :: assim_exclude_nflag(chm_famnum),assim_exclude_flag(chm_famnum,chm_cfg_size)
  character(len=9) :: assim_stnid(chm_famnum,chm_cfg_size)
  character(len=2) :: assim_fam(chm_famnum)
  
  integer :: generalized_operator(0:chm_constituents_size)   ! Same as genoper in NAMCHEM
  integer :: tropo_mode(0:chm_constituents_size),tropo_bound(0:chm_constituents_size)
  integer :: obsdata_maxsize
  integer :: transform(0:chm_constituents_size),message_fileunit
  real(8) :: amu(0:chm_constituents_size),tropo_column_top(0:chm_constituents_size)
  real(8) :: low_cutoff(0:chm_constituents_size),high_cutoff(0:chm_constituents_size)
  real(8) :: sigma_cutoff(0:chm_constituents_size)
  character(len=chm_filename_size) :: message_filename 

contains

!----------------------------------- Setup -----------------------------------
!---------------------------------------------------------------------------
!
! SUBROUTINE chm_setup
!!
!! *Purpose*: Setup additional information required by constituent obs and 
!!            not provided in obsSpaceData.
!!
!! @author Y. Rochon, Dec 2014 
!! 
!! Revisions: 
!!v            M. Sitwell, Feb 2015
!!v            - Removed references to earlier structure
!!v            Y. Rochon, Dec 2015
!!v            - Added call to chm_read_ref_fields and related optional input
!!v              'datestamp'
!!v          
!----------------------------------------------------------------------------------------
  subroutine chm_setup(datestamp_opt)
  
  implicit none

  integer, intent(in), optional :: datestamp_opt

  write(*,*) 'Begin chm_setup'

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
      call oss_obsdata_alloc(chm_efftemp,obsdata_maxsize,dim1=1)
      chm_efftemp%nrep=0
  end if

  write(*,*) 'Completed chm_setup'
  
  end subroutine chm_setup

!-----------------------------------------------------------------------------------
!
! SUBROUTINE chm_read_namchem
!!
!! *Purpose*: Read and store miscellaneous flags and constants.
!!
!! @author Y. Rochon, ARQI/AQRD, Apr 2015 
!!
!! Revisions: 
!!            Y. Rochon, ARQI/AQRD, Sept 2015
!!            - Added reading of variables related to options for attributing only 
!!              tropospheric portion of the increment to total column observations
!!            M. Sitwell, ARQI/AQRD, Oct 2015
!!            - Moved amu from auxiliary file to namelist namchem    
!!
!! Output:
!!
!!v   Read from NAMCHEM namelist
!!v
!!v        genoper             If generalized observation operator should be used and selection of approach
!!v                            <=0: not applied
!!v                            1: use trial field xb for mass weighted increment distribution
!!v                            2: use a combination of the difference of an external reference xc 
!!v                               and the trial field xb, i.e. mass weighted increment distribution 
!!v                               as a(xc-xb) + b*xc where a and b depend on the size of 
!!v                               sum[(xc-xb)/sig(xb)]^2 over the profile.
!!v
!!v        assim_fam           List of families to which filt_diagn_only is to apply.
!!v
!!v        assim_exclude_flag  Array specifying bits for identifying diagnostic-only observations for
!!v                            observations that would otherwise be assimilated according to the other
!!v                            assim_* arrays
!!v
!!v        assim_exclude_nflag Number of bit flags to specify in assim_exclude_flag array
!!v
!!v        assim_all           Logical indicating if all assimilatable obs of the specified family
!!v                            will be assimilated (default is .true.)
!!v
!!v        assim_num           Number combinations (stnid, bufr element, multi/uni-level) 
!!v                            identified for assimilation. All others will not
!!v                            be assimilated. OmP and OmA diagnostics and output
!!v                            will still be produced for non-assimilated datasets.
!!v                         
!!v                             0:  none are to be assimilated
!!v                            >0:  sets of (stnid, bufr varno, multi/uni-levels) to be assimilated
!!v
!!v        assim_varno         Bufr elements of obs sets for assimilation. A value of
!!v                            0 implies that all are to be used.
!/v
!!v        assim_stnid         Stnids of obs sets for assimilation. '*' denote wild cards
!!v
!!v        assim_nlev           0:  multi-level and uni-level
!!v                             1:  uni_level
!!v                            >1: multi-level 
!!v
!!v        tropo_mode          Integer indicating if special treatment is to be given to the troposphere
!!v                            when assimilating total column measurements. Values indicate
!!v                             0:  No special treatment given (default)
!!v                             1:  Values of the adjoint model above obsoper%column_bound set to zero.
!!v                                 If specified, generalized innovation operator only applied below
!!v                                 obsoper%column_bound in the tangent linear model.
!!v                             2:  Values of tangent linear model and adjoint model above
!!v                                 obsoper%column_bound set to zero.
!!v                            Array index refers to BUFR code element of Table 08046 (iconstituent_id)
!!v                            identifying the constituent. Relevant for total column measurements only.
!!v
!!v        tropo_bound         Integer indicating which column top value to use if tropo_mode is non-zero.
!!v                              0: Use fixed value of tropo_column_top
!!v                              1: Use model determination of tropopause
!!v                              2: Use model determination of PBL
!!v                            Options 1 and 2 will default to the value set in tropo_column_top if the model
!!v                            derived column top could not be determined. Relevant for total column measurements only.
!!v                         
!!v        tropo_column_top    Default value to use for the column boundary (in Pa). Array index refers to BUFR code
!!v                            element of Table 08046 (iconstituent_id) identifying the constituent.
!!v                            Relevant for total column measurements only.
!!v
!!v        amu                 Molecular mass  of constituents in g/moles (needed for unit conversions)
!!v                            Array index refers to BUFR code element of Table 08046 (iconstituent_id)
!!v                            identifying the constituent.
!!v
!!v       obsdata_maxsize      Max allowed size of work arrays (in terms of number of obs) associated to
!!v                            ordered observation indices
!!v
!!v       sigma_cutoff         fractional value of background error std. dev. to apply to increment
!!v                            sizes for setting them to zero prior to storage. Not applied if 
!!v                            cutoff value is <= 0.0 or >0.1
!!v
!!v       low_cutoff           min value allowed for increments prior to storage in rebm expressed as a fraction
!!v                            of the background field (generally < 1)
!!v
!!v       high_cutoff          max value allowed for increments prior to storage in rebm expressed as a multiple
!!v                            of the background field (generally > 1)
!!v
!!v       transform            Index specifying form of analysis increment (and related adjoint operation)
!!v                               -1: no transformation (dx given input trial field denoted as x)
!!v                              >=0: Ensure positive values. Provide warning if
!!v                                   non-positive values encountered.
!!v                                1: dlnx
!!v
!!v       message_filename     File name for file containing various messages and warnings related to chemical
!!v                            constituents that are not included in the listing file.
!!v
!!v                            Not currently used in code.
!!v                            Writing in MPI would need to use different files for each processor 
!!v                            or output via a single processor.
!!v
!---------------------------------------------------------------------------------------- 
  subroutine chm_read_namchem

  implicit none

  integer :: FNOM, FCLOS
  integer :: IERR, JELM, nulstat, ios, isize, nulnam, i
  integer :: genoper(0:chm_constituents_size)
  
  character(len=128) :: ligne  
  character(len=10)  :: namfile 

  EXTERNAL FNOM,FCLOS

  namelist /namchem/ assim_fam,assim_all,assim_num,assim_stnid,assim_varno,assim_nlev, &
                     assim_exclude_nflag,assim_exclude_flag,                   &
                     tropo_mode,tropo_bound,tropo_column_top,amu,sigma_cutoff, &
                     low_cutoff,high_cutoff,transform,obsdata_maxsize,         &
                     message_filename,genoper

  
  ! Default NAMCHEM values
  
  genoper(:)=0
  obsdata_maxsize=90000
  sigma_cutoff(:)=0.01
  low_cutoff(:)=0.1
  high_cutoff(:)=10.0
  transform(:)=0   ! At least ensure positive values and provide related warning as needed.
  
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

  message_filename = 'chem_message_'   ! Not used
  
  amu(:) = -1.0
  amu(0) = 48.0    ! Molecular mass in g/mole for O3
  amu(1) = 18.02   ! H2O
  amu(2) = 16.04   ! CH4
  amu(3) = 44.01   ! CO2
  amu(4) = 28.01   ! CO
  amu(5) = 46.01   ! NO2
  amu(6) = 44.01   ! N2O
  amu(7) = 30.03   ! HCHO=Formaldehyde
  amu(8) = 64.06   ! SO2
  amu(9) = 17.03   ! NH3
  amu(11) = 30.0   ! NO
  amu(26) = 1.0    ! PM2.5 - Not applicable
  amu(27) = 1.0    ! PM10  - Not applicable

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

  generalized_operator(:) = genoper(:)
  do i=chm_famNum,1,-1
     if (assim_fam(i).ne.'') exit
  end do
  assim_famnum=i
  
  end subroutine chm_read_namchem

!-----------------------------------------------------------------------------------------
! -------------------- Routines related to layer top & bottom levels----------------------
!
! SUBROUTINE chm_read_layers
!!
!! *Purpose*: Read and store top and bottom layer boundaries for CH sub-families
!!
!!
!! @author  Y. Rochon, ARQI/AQRD, Dec 2014 
!!v         - Partially based on oer_read_obs_erreurs_conv.
!!
!! Revisions: 
!!v            M. Sitwell, ARQI/AQRD, Feb 2015
!!v            - Renaming of routine and removal of lines no longer required.
!!          
!! Comments:
!!
!!v A) The option of reading from observation files is TBD. This will change the approach in allocating
!!v    the arrays size as the sizes will become dependent on the number of related obs for
!!v    which the observation files will need to be read.
!!
!----------------------------------------------------------------------------------------
  subroutine chm_read_layers

  implicit none

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
!
! Check for available layer info.
!
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

!   Read STNID (* is a wildcard)
    
    read(nulstat,'(2X,A9)',iostat=ios,err=10,end=10) chm_layers%stnids(jelm) 

!   Read (1) Obs BUFR element.
!        (2) Vertical coord type (1, 2, or 3)
!        (3) Flag indication if EOR provided from this auxiliary file or
!            to be read from the observation file,
!        (4) Number of vertical levels
!
!   Important: Combination of STNID, BUFR element and number of vertical levels
!              to determine association to the observations.
!
    read(nulstat,*,iostat=ios,err=10,end=10) chm_layers%element(jelm),chm_layers%vco(jelm),  &
       chm_layers%source(jelm),chm_layers%n_lvl(jelm)  
    
    if (icount+chm_layers%n_lvl(jelm).gt.isize) then
       write(*,'(10X,"Max array size exceeded: ",I6)') isize
       CALL utl_abort('chm_read_layers: READING PROBLEM.')    
    end if

    read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne
    
    if (chm_layers%n_lvl(jelm).ge.1) then   
       do jlev=1,chm_layers%n_lvl(jelm)
          icount=icount+1
          
          ! Read top and bottom levels
          
          read(nulstat,*,iostat=ios,err=10,end=10)                 &
                 chm_layers%vlayertop(icount),chm_layers%vlayerbottom(icount)
       end do
    end if

!    if (chm_layers%source(jelm).eq.1) then
!    
!      Read from observation files
!
!      .....
!
!    end if

  end do
   
 10 if (ios.gt.0) then
       WRITE(*,*) 'File read error message number: ',ios
       CALL utl_abort('chm_read_layers: READING PROBLEM')    
    end if
 
 11 CLOSE(UNIT=NULSTAT)
    IERR=FCLOS(NULSTAT)    

  end subroutine chm_read_layers

!----------------------------------------------------------------------------------------
!
! SUBROUTINE chm_get_layer_boundaries
!!
!! *Purpose*: Return layer boundaries for an observation. Combination of STNID, element
!!            variable number and number of vertical levels to determine association to the
!!            observations. Default values for top and bottom layers for total column measurements
!!            are to be provided.
!!
!! @author M. Sitwell, Y. Rochon, Feb 2015
!!
!! Revision: 
!!
!! Inputs:
!!v   - cstnid          station id
!!v   - varno           BUFR element
!!v   - ivco            type of vertical coordinate (see burpread_mod.ftn90 or
!!v                     routine chm_obsoperators for definitions)
!!v   - nlev            number of levels in the observation
!!v   - default_top     default value for top layer for total column measurement
!!v   - default_bottom  default value for bottom layer for total column measurement
!!
!! Outputs:
!!v   - lfound          Logical being .true. if layer boundaries found.
!!v   - layertop        top layer values
!!v   - layerbottom     bottom layer values
!!
! ---------------------------------------------------------------------------------------
  subroutine chm_get_layer_boundaries(cstnid,varno,ivco,nlev,default_top,default_bottom, &
                                      lfound,layertop,layerbottom)

    implicit none

    character(len=12), intent(in) :: cstnid
    integer, intent(in)           :: varno,ivco,nlev
    real(8), intent(in)           :: default_top,default_bottom
    real(8), intent(out)          :: layertop(nlev),layerbottom(nlev)
    logical, intent(inout)        :: lfound
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

!----------------------------------------------------------------------------------------
!
! SUBROUTINE chm_dealloc_layers
!!
!! *Purpose*: Deallocate temporary storage space used for layer info
!!
!! @author Y. Rochon, Dec 2014
!!
!! Revision: 
!!
! ---------------------------------------------------------------------------------------
  subroutine chm_dealloc_layers

    implicit none

    if (chm_layers%n_stnid.eq.0) return

    call chm_dealloc_info(chm_layers)
 
  end subroutine chm_dealloc_layers

!-----------------------------------------------------------------------------------------
!------------------- Routines related to averaging kernel matrices ----------------------
!  
! SUBROUTINE chm_read_avgkern
!!
!! *Purpose*:  Read averaging kernels from auxiliary file or observation file.
!!
!! @author M. Sitwell, ARQI/AQRD, March 2015
!!
!! Revision: 
!!
!-----------------------------------------------------------------------------------------
  subroutine chm_read_avgkern

    implicit none

    integer, parameter :: ndim=2

    integer :: istnid

    ! read the averaging kernel information from the auxiliary file
    call chm_read_avgkern_auxfile

    ! set size of observation file array
    allocate(chm_obsSub_avgkern(chm_avgkern%n_stnid))

    ! read from observation file
    do istnid=1,chm_avgkern%n_stnid
       if (chm_avgkern%source(istnid) == 1) then
          
          ! retrieve data from stats blocks (with bkstp=14 and block_type='DATA')
          chm_obsSub_avgkern(istnid) = obsf_obsSub_read('CH',chm_avgkern%stnids(istnid),bufr_avgkern, &
                                     chm_avgkern%n_lvl(istnid), ndim, bkstp_opt=14, &
                                     block_opt='DATA', match_nlev_opt=.true.)
          
       end if
    end do

  end subroutine chm_read_avgkern

!---------------------------------------------------------------------------------------
!
! SUBROUTINE chm_read_avgkern_auxfile
!!
!! *Purpose*: Read and store averaging kernel matricesfor CH sub-families
!!
!! @author Y. Rochon, M. Sitwell, ARQI/AQRD, Feb 2015 
!!v            - Currently implemented for only one latitude band
!!
!! Revisions: 
!!          
!----------------------------------------------------------------------------------------
  subroutine chm_read_avgkern_auxfile

  implicit none

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
!
! Check for available layer info.
!
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

!   Read STNID (* is a wildcard)
    
    read(nulstat,'(2X,A9)',iostat=ios,err=10,end=10) chm_avgkern%stnids(jelm) 

!   Read (1) Obs BUFR element.
!        (2) Flag indication if avgkern provided from this auxiliary file or
!            to be read from an observation file,
!        (3) Number of vertical levels
!
!   Important: Combination of STNID, BUFR element and number of vertical levels
!              to determine association to the observations.
!
    read(nulstat,*,iostat=ios,err=10,end=10) chm_avgkern%element(jelm),  &
       chm_avgkern%source(jelm),chm_avgkern%n_lvl(jelm)  
    
    if (icount+chm_avgkern%n_lvl(jelm).gt.isize) then
       write(*,'(10X,"Max array size exceeded: ",I6)') isize
       CALL utl_abort('chm_read_avgkern_auxfile: READING PROBLEM.')    
    end if

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
   
 10 if (ios.gt.0) then
       WRITE(*,*) 'File read error message number: ',ios
       CALL utl_abort('chm_read_avgkern_auxfile: READING PROBLEM')    
    end if
 
 11 CLOSE(UNIT=NULSTAT)
    IERR=FCLOS(NULSTAT)    

  end subroutine chm_read_avgkern_auxfile

!----------------------------------------------------------------------------------------
!
! SUBROUTINE chm_dealloc_avgkern
!!
!! *Purpose*: Deallocate temporary storage space used for averaging kernels
!!
!! @author Y. Rochon, Dec 2014
!!
!! Revision: 
!!
!---------------------------------------------------------------------------------------
  subroutine chm_dealloc_avgkern

    implicit none

    integer :: istnid

    if (chm_avgkern%n_stnid.eq.0) return

    if (allocated(chm_obsSub_avgkern)) then
       do istnid=1,chm_avgkern%n_stnid
          if (chm_avgkern%source(istnid).eq.1) call oss_obsdata_dealloc(chm_obsSub_avgkern(istnid))
       end do
       deallocate(chm_obsSub_avgkern)
    end if

    call chm_dealloc_info(chm_avgkern)
  
  end subroutine chm_dealloc_avgkern

!----------------------------------------------------------------------------------------
!
! FUNCTION chm_find_avgkern(cstnid,varno,nlev) result(ISTNID)
!!
!! *Purpose*: Finds the averaging kernel for an observation if one is specified. Returns 0 if
!!            either not found or not specified. Combination of STNID, BUFR element and number
!!            of vertical levels to determine association to the observations.
!!
!! @author M. Sitwell, March 2015
!!
!! Revision: 
!!          
!! Input:
!!v   - cstnid          station id
!!v   - varno           BUFR descriptor element
!!v   - nlev            number of levels in the observation
!!
!! Output:
!!v   - ISTNID          Index of averaging kernel in chm_avgkern if found. Zero indicates
!!v                     averaging kernel not found.
!!v
! ---------------------------------------------------------------------------------------
  function chm_find_avgkern(cstnid,varno,nlev) result(ISTNID)
    
    implicit none

    character(len=12), intent(in) :: cstnid
    integer, intent(in) :: varno,nlev
    integer :: ISTNID,JN
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

!----------------------------------------------------------------------------------------
!
! SUBROUTINE chm_get_avgkern(istnid,stnid,nlev,zlat,zlon,idate,itime,avg_kern)
!!
!! *Purpose*: Return averaging kernel for an observation.
!!
!! @author M. Sitwell, March 2015
!!
!! Revisions: 
!!
!! Input:
!!v   istnid          index of averaging kernel in chm_avgkern
!!v   nlev            number of observation levels
!!v   idate           YYYYMMDD
!!v   itime           HHMM
!!
!! Output:
!!v   avg_kern        the averaging kernel
!!v
!---------------------------------------------------------------------------------------
  subroutine chm_get_avgkern(istnid,stnid,nlev,zlat,zlon,idate,itime,avg_kern)

    implicit none

    integer, intent(in)  :: istnid,nlev,idate,itime
    real(8), intent(in)  :: zlat,zlon
    character(len=*), intent(in) :: stnid
    real(8), intent(out) :: avg_kern(nlev,nlev)
    integer :: start_index,end_index

    if (istnid.gt.0 .and. istnid.le.chm_avgkern%n_stnid) then
       
       if (chm_avgkern%source(istnid).eq.0) then
          ! get averaging kernel from auxiliary file
          start_index = chm_avgkern%ibegin(ISTNID)
          end_index = nlev*(start_index+nlev-1)
          avg_kern = RESHAPE(chm_avgkern%rak(start_index:end_index),(/nlev,nlev/),ORDER =(/2,1/))
       else
          ! get averaging kernel from observation file
          avg_kern = oss_obsdata_get_array2d(chm_obsSub_avgkern(istnid), oss_obsdata_get_header_code(zlon,zlat,idate,itime,stnid))
       end if

    else
       call utl_abort("chm_get_avgkern: Invalid station ID index.")
    end if

  end subroutine chm_get_avgkern

!-----------------------------------------------------------------------------------------
!------------------------- Routines related to gridded reference fields  ----------------
!  
! SUBROUTINE chm_read_ref_fields(datestamp)
!!
!! *Purpose*:  Read reference fields as directed by the content of the auxiliary file.
!!             Fields are provided in RPN/fst files specified in the auxiliary file (with path
!!             and filename)
!!
!!             Reference fields can be in a separate RPN file with name provided by the auxiliary
!!             or in monthly static background stats file (glbchemcov or bgcov; see 'isrc' below).
!!
!! ****** NOT TESTED *********
!!
!! @author Y. Rochon, ARQI/AQRD, Feb 2015
!!
!! Revision: 
!!
!! Comments:
!!
!!v     - Fields assumed to be of the same units as those of the corresponding input trial fields
!!
!-----------------------------------------------------------------------------------------
  subroutine chm_read_ref_fields(datestamp_opt)
    
    implicit none

    integer, intent(in), optional :: datestamp_opt
    
    character(len=128) :: fname
    character(len=4) :: varName
    character(len=12) :: etiket
    integer :: i,id,nd,j,ndim,ijour,imonth,iday,itime,isrc
    real(8) :: day
    integer, external :: newdate
   
    integer, external :: FNOM, FCLOS
    integer :: IERR, JLEV, JELM, nulstat, ios, isize, icount, iend
    logical :: LExists
    
    logical, parameter :: linterp=.true.
    
    integer :: ni, nj, nkeys, kind
    real(8), allocatable :: array1(:,:,:),array2(:,:,:),lvls(:),xlat(:),xlong(:) ! Allocated in chm_fst_read
  
    character (len=128) :: ligne

!   Initialize dimensions to zero

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

!   Check for file names containing ref fields

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
    
!   Read number of constituents with associated input file(s)
   
    read(nulstat,*,iostat=ios,err=10,end=10) ndim
    if (ndim.le.0) go to 10
    
!   Initialization

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
    
!   Get needed fields for each file(s)

    do i=1,ndim 

!      Read id,nd,isrc. id: constituent code; nd: number of sets; 1 or 2;
!      isrc: 1 for fname being in the auxiliary file, 0 for glbchemcov (or bgcov if glbchemcov not present) 
       
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

!            Following for interpolation as a function of days from mid-months.
             
             if (iday.gt.15) then
                 if (imonth.eq.12) then
                    call utl_readFstField(trim(fname),varName,-1,1,-1,etiket,ni,nj,nkeys,array2,lvls_opt=lvls,kind_opt=kind)
                else
                   call utl_readFstField(trim(fname),varName,-1,imonth+1,-1,etiket,ni,nj,nkeys,array2,lvls_opt=lvls,kind_opt=kind)
                end if
          
!               Linearly interpolate in time (approximately - assumes 30 day months)

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

!               Linearly interpolate in time (approximately - assumes 30 day months)

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
     
 10 if (ios.gt.0) then
       WRITE(*,*) 'File read error message number: ',ios
       CALL utl_abort('chm_read_ref_fields: READING PROBLEM')    
    end if
 
 11 CLOSE(UNIT=NULSTAT)
    IERR=FCLOS(NULSTAT)    
    
  end subroutine chm_read_ref_fields
   
!---------------------------------------------------------------------------------------
!
! SUBROUTINE chm_set_reference_obsdata(obsoper)
!!
!! *Purpose*: Determine and store reference profile at obs location if needed by
!!          the observation operators.
!!
!! ***** NOT TESTED *****
!!
!! @author Y. Rochon, ARQI/AQRD May 2016
!!
!! Revisions:
!!
!! Input:
!!v
!!v      obsoper%constituent_id  Constituent id
!!v      obsoper%nmodlev         Number of model levels for variables other than uu and vv
!!v      obsoper%pressmod        Model pressure array
!!v      obsoper%tt              Model temperature (Kelvin)
!!v      obsoper%gz              Model geopotential height (m)
!!v      obsoper%hu              Specific humidity 
!!v      obsoper%lat             Latitude (rad)
!!v      obsoper%lon             Longitude (rad)
!!v
!! Output:
!!v
!!v      chm_ref_trial           Reference profile object
!!v
!!v Comments
!! 
!----------------------------------------------------------------------------------------
  subroutine chm_set_reference_obsdata(obsoper)
    
    implicit none

    type(struct_chm_obsoperators), intent(inout) :: obsoper
       
    integer :: i,istart,id
    real(8) :: tropo_press, refprof(obsoper%nmodlev),refprof2(obsoper%nmodlev),dt
    real(8), allocatable :: pressrefin(:)
    logical, allocatable :: lsuccess(:)
    
    if (obsoper%constituent_id.lt.0.or.obsoper%constituent_id.gt.chm_var_maxnumber()) return

    id=obsoper%constituent_id
    
    if (generalized_operator(id).le.1) return           
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
        where (pressrefin.lt.obsoper%gz(obsoper%nmodlev)) pressrefin=obsoper%gz(obsoper%nmodlev)
        pressrefin(:) = phf_convert_z_to_pressure(pressrefin,obsoper%gz,obsoper%pp, &
                        chm_ref_fields(id,1)%nlev,obsoper%nmodlev,obsoper%lat,lsuccess)
    else if (chm_ref_fields(id,1)%ivkind.eq.4) then
        pressrefin(:)=pressrefin(:) + obsoper%gz(obsoper%nmodlev)
        pressrefin(:) = phf_convert_z_to_pressure(pressrefin,obsoper%gz,obsoper%pp, &
                        chm_ref_fields(id,1)%nlev,obsoper%nmodlev,obsoper%lat,lsuccess)
    else if (chm_ref_fields(id,1)%ivkind.eq.1) then
        pressrefin(:)=pressrefin(:)*obsoper%pp(obsoper%nmodlev) ! Convert from sigma to Pa   
    else
       write(*,*) 'Vertical coordinate kind is ',chm_ref_fields(id,1)%ivkind
       call utl_abort('chm_get_reference_obsdata: Cannot handle vertical coordinate of this kind.')
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
           tropo_press=phf_get_tropopause(obsoper%nmodlev,obsoper%pp,obsoper%tt,obsoper%gz,hu_opt=obsoper%hu)
        else
           tropo_press=phf_get_tropopause(obsoper%nmodlev,obsoper%pp,obsoper%tt,obsoper%gz)
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
               where (pressrefin.lt.obsoper%gz(obsoper%nmodlev)) pressrefin=obsoper%gz(obsoper%nmodlev)
               pressrefin(:) = phf_convert_z_to_pressure(pressrefin,obsoper%gz,obsoper%pp, &
                               chm_ref_fields(id,2)%nlev,obsoper%nmodlev,obsoper%lat,lsuccess)
           else if (chm_ref_fields(id,2)%ivkind.eq.4) then
               pressrefin(:)=pressrefin(:) + obsoper%gz(obsoper%nmodlev)
               pressrefin(:) = phf_convert_z_to_pressure(pressrefin,obsoper%gz,obsoper%pp, &
                               chm_ref_fields(id,2)%nlev,obsoper%nmodlev,obsoper%lat,lsuccess)
           else if (chm_ref_fields(id,2)%ivkind.eq.1) then
               pressrefin(:)=pressrefin(:)*obsoper%pp(obsoper%nmodlev) ! Convert from sigma to Pa        
           else 
               write(*,*) 'Vertical coordinate kind is ',chm_ref_fields(id,2)%ivkind
               call utl_abort('chm_get_reference_obsdata: Cannot handle vertical coordinate of this kind.')
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
       call oss_obsdata_alloc(chm_ref_trial, obsdata_maxsize, dim1=obsoper%nmodlev)
       chm_ref_trial%nrep = 0
    end if

    ! Here, nrep will count the number of filled elements in the data arrays
    chm_ref_trial%nrep = chm_ref_trial%nrep+1 

    if (chm_ref_trial%nrep.gt.obsdata_maxsize) &
         call utl_abort('chm_get_ref_obsdata: Reach max size of array ' // trim(utl_str(obsdata_maxsize)) )
  
    ! The obsoper%obs_index (header index) serves as the unique locator code 
    write(chm_ref_trial%code(chm_ref_trial%nrep),'(I22)') obsoper%obs_index

    ! Save profile in chm_ref_trial
    
    chm_ref_trial%data1d(:,chm_ref_trial%nrep) = refprof(:)

  end subroutine chm_set_reference_obsdata
 
!---------------------------------------------------------------------------------------
!
! FUNCTION chm_get_ref_column(code) result(array)
!!
!! *Purpose*: Extract and provide column from chm_ref_field associated to code.
!!
!! @author Y. Rochon, ARQI/AQRD Feb 2017
!!
!! Revisions:
!!
!! Input:
!!v
!!v      code         unique identifying code
!!v
!! Output:
!!v
!!v      array        retrieved array from obsdata%data1d of dimension obsdata%dim1
!!v      stat         search success (0 - found; 1 = no data; 2 = not found)
!!v 
!----------------------------------------------------------------------------------------    
  function chm_get_ref_column(code) result(array)

    implicit none
        
    character(len=*), intent(in) :: code
    real(8) :: array(chm_ref_trial%dim1)

    integer :: stat

    array = oss_obsdata_get_array1d(chm_ref_trial,code,stat)
    if (stat.gt.0) call utl_abort("chm_get_ref_column: Code not found - " // &
                                  trim(code))
    
  end function chm_get_ref_column

!-----------------------------------------Misc ----------------------------------------
!--------------------------------------------------------------------------------------
!  
! SUBROUTINE chm_dealloc_info(info)
!!
!! *Purpose*: Deallocates struct_chm_info instance
!!
!! @author M. Sitwell  May 2015
!!
!---------------------------------------------------------------------------------------
  subroutine chm_dealloc_info(info)
    
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

! ------------------------------------------------------------------------------------------
!
! LOGICAL FUNCTION chm_diagn_only(cfamName,cstnid,varno,nobslev,flag)
!! 
!! *Purpose*: Identify whether or not the obs set identified by 
!!            the combination of (cstnidin, bufrin,nlevs) will
!!            be assimilated or used for independent verifications.
!!
!! @author Y.J. Rochon, ARQI/AQRD, July 2015
!!    
!! Revisions:
!!
!! Input:
!!v
!!v       cfamName         Family name
!!v       cstnid           Input station id
!!v       varno            Obs BUFR number
!!v       nobslev          Number of levels
!!v       flag             observation integer flag
!!v
!! Output:
!!v
!!v       chm_diagn_only    Indicating if assimilation to be skipped but data
!!v                         is to be used for independent verifications after 
!!v                         assimilation/minimization.
!!v
!-------------------------------------------------------------------------------------------
 logical function chm_diagn_only(cfamName,cstnid,varno,nobslev,flag)
 
    implicit none

    integer, intent(in) :: varno,nobslev,flag
    character(len=*), intent(in) :: cstnid,cfamName
  
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
       chm_diagn_only = .false.
       return    
    end if
    
    if (assim_all(ifam)) then
       ! assimilate all observations
       chm_diagn_only = .false.
    else if (assim_num(ifam).le.0) then
       ! assimilate no observations
       chm_diagn_only = .true.
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
       chm_diagn_only = elemId.eq.0
    end if
    
    if (chm_diagn_only) return

    ! check if the observation integer flag has a bit marked by assim_exclude_flag (same flagging as in filt_suprep)
    if (assim_exclude_nflag(ifam).gt.0) then
       do i=1,assim_exclude_nflag(ifam)
          if (btest(flag, 13 - assim_exclude_flag(ifam,i) )) then
             chm_diagn_only = .true.
             return
          end if
       end do
    end if

  end function chm_diagn_only

!------------------------------------------------------------------------------------------
!
! SUBROUTINE chm_apply_2dfieldr4_transform(iconstituent_id,varName,jlev,jstep,field,l_reverse)
!!
!! *Purpose*: Apply transform (or its inverse) of 2D field.
!!            Called by routine readTrialField in file innovation_mod.ftn90. 
!!
!! @author Y. Rochon, Feb 2016
!!
!! Revisions:
!!
!! Input
!!
!!v     iconstituent_id BUFR code element of Table 08046 identifying the constituent.
!!v     varName         Field name (nomvar)
!!v     l_reverse       Reverse/inverse transformation if present (default value .false.).
!!v     jlev            vertical level index
!!v     jstep           Time step index
!!v
!! InOut
!!v
!!v     field           2D field
!!v
!! Comments:
!!v
!!v 1. The EnVar assumes that the input background error covariances are provided
!!v    for the transformed field if a variable transformation is requested!
!!v
!-----------------------------------------------------------------------------------
  subroutine chm_apply_2dfieldr4_transform(iconstituent_id,varName,jlev,jstep,field,l_reverse_opt)

    implicit none
    
    integer, intent(in) :: iconstituent_id,jlev,jstep
    character(len=*), intent(in) :: varName
    logical, intent(in), optional :: l_reverse_opt

    real(4), intent(inout) :: field(:,:)            

    integer :: i,j,ier,unit,icount
    real(4) :: valmin
    real(4), parameter :: valmin_ref=1.0E-20
    integer, external :: fclos
    logical :: lrev
    
    if (iconstituent_id.lt.0.or.iconstituent_id.gt.chm_constituents_size) return
    
    if (transform(iconstituent_id).lt.0) return
    
    if (present(l_reverse_opt)) then
       lrev = l_reverse_opt
    else
       lrev = .false.
    end if

    ! Check for non-positive values if forward transformation 

    if (.not.lrev) then

       !x Following open cannot be used at the same time by more than one processor 
       !x unless opening file for reading only.
       !x call utl_open_asciifile(message_filename,unit)

       icount=0
       LOOP1: do j=size(field,2),1,-1
          valmin=minval(field(:,j),mask=field(:,j).gt.0.0)
          if (valmin.gt.1.E30) valmin=valmin_ref
          icount=icount+count(field(:,j).lt.valmin)
          where (field(:,j).lt.valmin) field(:,j)=valmin
          !x do i=1,size(field,1) 
          !x   if (field(i,j).lt.valmin) then
          !x      write(unit,'(A,G9.2,A,A)') "Unexpected/undesired negative value of ",field(i,j), " for input field ",trim(varName)
          !x      write(unit,'(A,3I4,A,I3,A,G9.2)')  "at location (",i,j,jlev,") and step ",jstep,". Value replaced by ",valmin,"."                    
          !x      field(i,j)=valmi
          !x      icount=icount+1
          !x   end if
          !x   if (icount.gt.50) then
          !x      write(unit,*) "Stopped counting... WARNING: Large number of zero or negative values for ",trim(varName),"trial field."
          !x      icount=icount+count(field.lt.valmin)
          !x      where (field.lt.valmin) field=valmin
          !x	     exit LOOP1
          !x   end if
          !x end do
       end do LOOP1
       write(*,*) "chm_apply_2dfield4_transform: WARNING - ",icount," zero and or negative values found. Reset to small positive value for ", &
                  trim(varName)," trial field."
       
       !x ier=fclos(unit)
       
    end if

    ! Apply transformation
    
    select case(transform(iconstituent_id))
    case(0)
       ! No tranformation
    case(1)

       ! Transform lnx to/from x or dlnx to/from dx
           
       if (.not.lrev) then
       
          ! Forward transformation
          
          field=log(field)
     
       else
    
          ! Reverse transformation          
       
          field=exp(field)
             
       end if
   
    case default
       call utl_abort('chm_apply_2dfieldr4_transform: Transformation #' // trim(utl_str(iconstituent_id)) // &
            ' for constituent ' // trim(utl_str(iconstituent_id)) // " and variable name " // trim(varname) // &
            ' is not defined.')
    end select
      
  end subroutine chm_apply_2dfieldr4_transform

!---------------------------------------------------------------------------------------
!
! INTEGER FUNCTION chm_var_maxnumber()  
!!
!! *Purpose*: Pass on chm_constituents_size
!!
!! @author Y. Rochon, ARQI/AQRD, Feb 2017
!!
!! Revisions:
!!
!----------------------------------------------------------------------------------------
  integer function chm_var_maxnumber()  

    implicit none
    
    chm_var_maxnumber=chm_constituents_size
    
   end function chm_var_maxnumber

!---------------------------------------------------------------------------------------
!
! FUNCTION chm_setup_get_str(stype) result(val) 
!!
!! *Purpose*: Pass on associated 'stype'_filename
!!
!! @author Y. Rochon, ARQI/AQRD, Feb 2017
!!
!! Revisions:
!!
!! In
!!
!!v   stype     reference string for filename identifier
!!
!----------------------------------------------------------------------------------------
  function chm_setup_get_str(stype) result(val) 

    implicit none
    character(len=*), intent(in) :: stype
    character(len=chm_filename_size) :: val
    
    select case(trim(stype))    
    case('message')
       val=message_filename
    case('auxfile')
       val=chm_aux_filename
    case default
       call utl_abort('chm_setup_get_str: Selection not found ' // trim(stype))
    end select
    
  end function chm_setup_get_str
  
!---------------------------------------------------------------------------------------
!   
! FUNCTION chm_setup_get_int(stype,index)  result(val)
!!
!! *Purpose*: Pass on associated *'stype'*
!!
!! @author Y. Rochon, ARQI/AQRD, Feb 2017
!!
!! Revisions:
!!
!! In
!!
!!v   stype     reference string for integer-based identifier/value
!!v   index     array index (optional)
!!v
!----------------------------------------------------------------------------------------
  function chm_setup_get_int(stype,index_opt)  result(val)

    implicit none
    character(len=*), intent(in) :: stype
    integer, intent(in), optional :: index_opt
    integer :: val
    
    select case(trim(stype))    
    case('message_unit')
       val=message_fileunit
    case('obsdata_maxsize')
       val=obsdata_maxsize
    case('genoper')
       if (present(index_opt)) then
          val=generalized_operator(index_opt)
       else
          call utl_abort('chm_setup_get_int: Missing index for ' // trim(stype))
       end if
    case('tropo_mode')
       if (present(index_opt)) then
          val=tropo_mode(index_opt)
       else
          call utl_abort('chm_setup_get_int: Missing index for ' // trim(stype))
       end if
    case('tropo_bound')
       if (present(index_opt)) then
          val=tropo_bound(index_opt)
       else
          call utl_abort('chm_setup_get_int: Missing index for ' // trim(stype))
       end if
    case('transform')
       if (present(index_opt)) then
          val=transform(index_opt)
       else
          call utl_abort('chm_setup_get_int: Missing index for ' // trim(stype))
       end if
    case default
       call utl_abort('chm_setup_get_int: Selection not found ' // trim(stype))
    end select
    
  end function chm_setup_get_int

!---------------------------------------------------------------------------------------
!   
! FUNCTION chm_setup_set_int(stype,val)  result(ier)
!!
!! *Purpose*: Assign value to *'stype'*
!!
!! @author Y. Rochon, ARQI/AQRD, Feb 2017
!!
!! Revisions:
!!
!! In
!!
!!v   stype     reference string for integer-based identifier/value
!!v   val       integer value to assign
!!
!----------------------------------------------------------------------------------------
  function chm_setup_set_int(stype,val)  result(ier)

    implicit none
    character(len=*), intent(in) :: stype
    integer, intent(in) :: val
    integer :: ier
    
    ier=0
    select case(trim(stype))    
    case('message_unit')
       message_fileunit=val
    case default
       ier=-1
       call utl_abort('chm_setup_set_int: Selection not found ' // trim(stype))
    end select
    
  end function chm_setup_set_int

!---------------------------------------------------------------------------------------
!    
! FUNCTION chm_setup_get_float(stype,index)  result(val)
!!
!! *Purpose*: Pass on associated *'stype'*
!!
!! @author Y. Rochon, ARQI/AQRD, Feb 2017
!!
!! Revisions:
!!
!! In
!!
!!v   stype     reference string for float-based identifier/value
!!v   index     array index (optional)
!!
!----------------------------------------------------------------------------------------
  function chm_setup_get_float(stype,index_opt)  result(val)

    implicit none
    character(len=*), intent(in) :: stype
    integer, intent(in), optional :: index_opt
    real(8) :: val
    
    select case(trim(stype))    
    case('amu')
       if (present(index_opt)) then
          val=amu(index_opt)
       else
          call utl_abort('chm_setup_get_float: Missing index for ' // trim(stype))
       end if
    case('tropo_column_top')
       if (present(index_opt)) then
          val=tropo_column_top(index_opt)
       else
          call utl_abort('chm_setup_get_float: Missing index for ' // trim(stype))
       end if
    case('low_cutoff')
       if (present(index_opt)) then
          val=low_cutoff(index_opt)
       else
          call utl_abort('chm_setup_get_float: Missing index for ' // trim(stype))
       end if
    case('high_cutoff')
       if (present(index_opt)) then
          val=high_cutoff(index_opt)
       else
          call utl_abort('chm_setup_get_float: Missing index for ' // trim(stype))
       end if
    case('sigma_cutoff')
       if (present(index_opt)) then
          val=sigma_cutoff(index_opt)
       else
          call utl_abort('chm_setup_get_float: Missing index for ' // trim(stype))
       end if
    case default
       call utl_abort('chm_setup_get_float: Selection not found ' // trim(stype))
    end select
    
  end function chm_setup_get_float

!-----------------------------------------------------------------------------------
!------------------ Routines associated to chm_efftemp -----------------------------
!
! SUBROUTINE chm_add_efftemp_obsdata(code,temp_eff)
!!
!! *Purpose*: Add effective temperature value to its obsdata object
!!
!! @author Y. Rochon, ARQI/AQRD Feb 2017
!!
!! Revisions:
!!
!! Input
!!
!!v      code         unique identifying code
!!v      temp_eff     effective temperature
!! 
!----------------------------------------------------------------------------------------    
  subroutine chm_add_efftemp_obsdata(code,temp_eff)

    implicit none
        
    character(len=*), intent(in) :: code
    real(8), intent(in) :: temp_eff(:)

    call oss_obsdata_add_data1d(chm_efftemp,temp_eff,code,obsdata_maxsize)
    
  end subroutine chm_add_efftemp_obsdata
   
!---------------------------------------------------------------------------------------
!
! SUBROUTINE chm_add_efftemp_obsfile()
!!          
!! *Purpose*: Add effective temperatures in obs file.
!!
!! @author Y. Rochon, Feb 2017 
!! 
!! Revisions:
!!
!----------------------------------------------------------------------------------------
  subroutine chm_add_efftemp_obsfile()

    implicit none
    
    integer :: ierr,nrep_modified,varno(1)

!   If needed, add effective temperature values in obs file for total column measurements

    !if ( .not.obsf_fileTypeIsBurp() ) call utl_abort('chm_add_efftemp_obsfile: only compatible with BURP files')

    call oss_obsdata_MPIallgather(chm_efftemp)
    
    if (chm_efftemp%nrep.gt.0) then
        varno(1)=12001
        nrep_modified = obsf_obsSub_update(chm_efftemp,'CH',varno(1:max(1,chm_efftemp%dim2)),bkstp_opt=0, &
                                        block_opt='INFO',multi_opt='UNI') 
        write(*,*) 'chm_add_efftemp_obsfile: Added ',nrep_modified,' effective temperature values in the obs file.'
    end if 

  end subroutine chm_add_efftemp_obsfile

!-----------------------------------------------------------------------------------

end module chem_setup_mod
