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

module BCovarSetupChem_mod 
  ! MODULE BCovarSetupChem_mod (prefix='bcsc' category='6. High-level data objects')
  !
  ! :Purpose: Contains routines for the reading and preparation of
  !           background-error covariance elements. Correlation matrices 
  !           based on horizontally homogeneous/isotropic correlations.  
  !
  ! :Comments:
  !
  !   1. Covariances uncoupled from weather variable.
  !
  !   2. Handles univariate and multivariate covariances.
  !      See routines bcsc_readcorns2 and bcsc_sucorns2. 
  ! 
  !   3. For multiple univariate variables (or univarite blocks of one to
  !      multiple variables), one can alternatively have multiple sets of
  !      covariance matrices within this module instead of a single covariance
  !      matrix setup (similarly to what was done for corvert*).
  !
  ! Public Subroutines:
  
  !    bcsc_setupCH:    Must be called first. Sets of background covariance
  !                     matrix (and balance operators if any are eventually
  !                     added)
  !    bcsc_getCovarCH: Provides covariances and related
  !    bcsc_getScaleFactor : Provides std dev scaling factor
  !                     elements for background check and obs operators. 
  !    bcsc_finalize    Deallocate internal module arrays.
  !    bcsc_StatsExistForVarName: Determine is Covar available for specified variable
  !    bcsc_getBgStddev: Interpolation background error std dev to obs location
  !    bcsc_retrieveBgStddev: Retrieve previously saved background stddev
  !                      profiles in bgStddev from the header index.
  !    bcsc_addBgStddev: Add background stddev profiles (and inverse) to 
  !                      bgStddev which can be retrieved later using a header index.

  use midasMpi_mod
  use MathPhysConstants_mod
  use earthConstants_mod
  use obsSubSpaceData_mod
  use gridStateVector_mod
  use globalSpectralTransform_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use varNameList_mod
  use utilities_mod
  use calcHeightAndPressure_mod

  implicit none
  save
  private

  ! public procedures
  ! ------------------
    
  public :: bcsc_setupCH,bcsc_finalize
  public :: bcsc_getCovarCH
  public :: bcsc_getScaleFactor
  public :: bcsc_StatsExistForVarName
  public :: bcsc_retrieveBgStddev,bcsc_addBgStddev,bcsc_getBgStddev
  
  ! public types
  ! ------------
    
  public :: struct_bcsc_bgStats

  ! module shared variables
  ! ----------------------- 
  
  integer             :: nlev_M,nlev_T
  integer             :: gstID         
  real(8), allocatable   :: rlatr(:),rlongr(:)     
  type(struct_vco),pointer :: vco_anl
  
  character(len=15) :: bcsc_mode
                            
  integer, external   :: get_max_rss
  integer             :: nulbgst=0

  ! Bacgkround error covariance matrix elements.
  ! One could add an additional dimension to corns  
  ! for separate block-univariate correlation matrices.
  ! This would also permit merging of bmatrixhi_mod and bmatrixchem_mod
  ! into one module.

  real(8),allocatable :: stddev(:,:,:)
  real(8),allocatable :: rstddev(:,:)

  ! Parameters of the NAMBCHM namelist

  integer             :: ntrunc
  real(8)             :: rpor(vnl_numvarmax)
  real(8)             :: rvloc(vnl_numvarmax)
  integer,parameter   :: maxNumLevels=200
  real(8)             :: scaleFactor(vnl_numvarmax,maxNumLevels)
  integer             :: numModeZero  ! number of eigenmodes to set to zero
  logical             :: ReadWrite_sqrt

  ! Indicate of physical spaces covariances (stddev, corverti and corverti) are
  ! to be calculated/saved
  logical             :: getPhysSpaceStats

  ! Indicate if physical space correlation lengths are calculated from spectral
  ! space covariances. Needed for some CH family obs operator settings.
  logical             :: getPhysSpaceHCorrel

  ! Indicates if physial space stats are output in file 'bCovarSetupChem_out.fst'.
  ! Provided as extra utility - not used my MIDAS main programs.
  logical             :: WritePhysSpaceStats
  character(len=4)    :: stddevMode
  character(len=4)    :: IncludeAnlVarKindCH(vnl_numvarmax)
  character(len=4)    :: CrossCornsVarKindCH(vnl_numvarmax)
  character(len=20)   :: TransformVarKindCH

  ! Square root of scaleFactor   
  real(8) :: scaleFactor_stddev(vnl_numvarmax,maxNumLevels) 
      
  real(8), parameter :: zps = 101000.D0 ! Reference surface pressure

  ! module structures
  ! -----------------
  
  type :: struct_bcsc_bgStats

     !  Structure storing background error stats 
     !  and related elements

     ! logical indicating if stats are available     
     logical              :: initialized=.false.
     
     integer              :: numvar3d   ! number of 3D fields
     integer              :: numvar2d   ! number of 2D fields
     integer              :: ni         ! number of latitudes
     integer              :: nj         ! number of longitudes
     integer              :: nlev       ! number of vertical levels
     
     ! Total number of elements over all verticallevels and 
     ! variables (varNAmeList)
     integer              :: nkgdim
     
     integer              :: ntrunc     ! spectral dimension
     character(len=4), allocatable :: varNameList(:) ! list of variable names
                      
     integer, allocatable :: nsposit(:)  ! start positions of fields
     real(8), allocatable :: lat(:)      ! grid lat in radians
     real(8), allocatable :: lon(:)      ! grid lon in radians
     real(8), allocatable :: vlev(:)     ! vertical levels
     real(8), allocatable :: corns(:,:,:)    ! spectral space correlation matrix

     ! Phys. Space vertical correlation matrix
     real(8), allocatable :: corvert(:,:,:)
     
     real(8), allocatable :: corverti(:,:,:) ! Inverse of 'corvert'
     
     ! 1 / sum of vertical correlations in each row
     real(8), allocatable :: invsum(:,:)
     
     real(8), allocatable :: stddev(:,:,:)   ! error std dev
     real(8), allocatable :: hcorrlen(:,:)   ! horizontal correlation lengths
  
  end type struct_bcsc_bgStats

  ! Assigned type variables
  ! -----------------------
    
  type(struct_bcsc_bgStats)  :: bgStats  
  type(struct_oss_obsdata)   :: bgStddev ! Arrays for background error 
                                           ! std dev in obs space

  !*************************************************************************
    
  contains

  !--------------------------------------------------------------------------
  ! bcsc_setupCH
  !--------------------------------------------------------------------------
  subroutine bcsc_setupCH(hco_in,vco_in,covarNeeded,mode)
    !
    ! :Purpose: Set up for constituents static background error covariances.
    !
    implicit none

    !Arguments
    type(struct_hco), intent(in), pointer :: hco_in
    type(struct_vco), intent(in), pointer :: vco_in
    logical, intent(out)     :: covarNeeded
    character(len=*), intent(in) :: mode ! 'Analysis' or 'BackgroundCheck'

    !Locals
    integer :: nulnam, ierr, fnom, fclos, status
    integer :: varIndex,nChmVars,varIndex2
    character(len=4) :: BchmVars(vnl_numvarmax)
    real(8), pointer    :: pressureProfile_T(:)
        
    NAMELIST /NAMBCHM/ntrunc,rpor,rvloc,scaleFactor,numModeZero,ReadWrite_sqrt, &
                      stddevMode,IncludeAnlVarKindCH,getPhysSpaceHCorrel, &
		      CrossCornsVarKindCH,WritePhysSpaceStats, &
                      TransformVarKindCH,getPhysSpaceStats

    write(*,*) 'Started bcsc_setupCH'  
     
    ! First check if there are any CH fields 
    
    covarNeeded = .true.
    varIndex2=0
    do varIndex = 1, vnl_numvarmax
      if (gsv_varExist(varName = vnl_varNameList(varIndex))) then
        if (vnl_varKindFromVarname(vnl_varNameList(varIndex)) == 'CH') then
          varIndex2 = 1
          exit
        end if 
      end if      
    end do
    if (varIndex2 == 0) then
      ! Assume there is no need for Bchm
      covarNeeded = .false.
      return
    end if

    bgStats%numvar3d = 0
    bgStats%numvar2d = 0
 
    allocate(bgStats%varNameList(vnl_numvarmax))
    bgStats%varNameList(:) = ''
    allocate(bgStats%nsposit(vnl_numvarmax+1))
    bgStats%nsposit(1) = 1

    if ( trim(mode) == 'Analysis' .or. trim(mode) == 'BackgroundCheck') then
      bcsc_mode = trim(mode)
      if(mmpi_myid == 0) write(*,*)
      if(mmpi_myid == 0) write(*,*) 'bcsc_setupCH: Mode activated = ', &
        trim(bcsc_mode)
    else
      write(*,*)
      write(*,*) 'mode = ', trim(mode)
      call utl_abort('bcsc_setupCH: unknown mode')
    end if

    ! Initialization of namelist NAMBCHM parameters
    
    ntrunc=108
    rpor(:)=3000.D3
    rvloc(:)=4.0D0
    scaleFactor(:,:) = 0.0d0
    numModeZero = 0
    ReadWrite_sqrt = .false.
    WritePhysSpaceStats = .false.
    getPhysSpaceHCorrel = .false.
    getPhysSpaceStats = .false.
    stddevMode = 'GD3D'    
    IncludeAnlVarKindCH(:) = ''
    CrossCornsVarKindCH(:) = ''
    TransformVarKindCH = ''
            
    ! Read namelist input
    
    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=NAMBCHM,iostat=ierr)
    if(ierr /= 0) call utl_abort('bcsc_setupCH: Error reading namelist')
    if(mmpi_myid == 0) write(*,nml=NAMBCHM)
    ierr = fclos(nulnam)

    ! Set BchmVars
    nChmVars=0
    BChmVars(:)=''
    if (trim(IncludeAnlVarKindCH(1)) == '') then
      do varIndex = 1, vnl_numvarmax
        if (.not. gsv_varExist(varName = vnl_varNameList(varIndex))) cycle
        if (vnl_varKindFromVarname(vnl_varNameList(varIndex)) /= 'CH') cycle
        nChmVars = nChmVars+1
        BchmVars(nChmVars) = trim(vnl_varNameList(varIndex))
      end do
    else
      do varIndex = 1, vnl_numvarmax
        if (.not. gsv_varExist(varName = vnl_varNameList(varIndex))) cycle
        if (vnl_varKindFromVarname(vnl_varNameList(varIndex)) /= 'CH') cycle
        do varIndex2 = 1, vnl_numvarmax
          if (trim(IncludeAnlVarKindCH(varIndex2)) == '') exit
          if (trim(vnl_varNameList(varIndex)) == &
	    trim(IncludeAnlVarKindCH(varIndex2))) then
            nChmVars = nChmVars+1
            BchmVars(nChmVars)= trim(vnl_varNameList(varIndex))
            exit
          end if
        end do
      end do
    end if

    if (nChmVars == 0) then 
      if(mmpi_myid == 0) then
        write(*,*) 'Size of BchmVars is zero. B matrix for CH family ', &
	           'not produced.'
        write(*,*) 'No chemical assimilation to be performed.'
        write(*,*) 'Completed bcsc_setupCH'
      end if
      covarNeeded = .false.
      return
    end if
    
    ! Set vertical dimensions

    vco_anl => vco_in
    nLev_M = vco_anl%nlev_M
    nLev_T = vco_anl%nlev_T

    ! Find the 3D variables (within NAMSTATE namelist)

    do varIndex = 1, vnl_numvarmax3D    
      if (gsv_varExist(varName=vnl_varNameList3D(varIndex)) .and. &
          any(trim(vnl_varNameList3D(varIndex))==BchmVars(1:nChmVars)) ) then

        if (vnl_varKindFromVarname(vnl_varNameList3D(varIndex)) /= 'CH') cycle
	  
        bgStats%numvar3d = bgStats%numvar3d + 1
        bgStats%nsposit(bgStats%numvar3d+1)= &
	  bgStats%nsposit(bgStats%numvar3d)+nLev_T
        bgStats%varNameList(bgStats%numvar3d)= &
	  vnl_varNameList3D(varIndex)
      end if
    end do
 
    ! Find the 2D variables (within NAMSTATE namelist)

    do varIndex = 1, vnl_numvarmax2D
      if (gsv_varExist(varName=vnl_varNameList2D(varIndex)) .and. &
          any(trim(vnl_varNameList2D(varIndex)) == BchmVars(1:nChmVars)) ) then

        if (vnl_varKindFromVarname(vnl_varNameList2D(varIndex)) /= 'CH') cycle
        bgStats%numvar2d = bgStats%numvar2d + 1
        bgStats%nsposit(bgStats%numvar3d+bgStats%numvar2d+1) = &
	  bgStats%nsposit(bgStats%numvar3d+bgStats%numvar2d)+1
        bgStats%varNameList(bgStats%numvar2d) = vnl_varNameList2D(varIndex)
      end if       
    end do

    if (bgStats%numvar3d+bgStats%numvar2d == 0) then    
      if (mmpi_myid == 0) then
        write(*,*) 'B matrix for CH family not produced.'
        write(*,*) 'No chemical assimilation to be performed.'
        write(*,*) 'Completed bcsc_setupCH'
      end if
      covarNeeded = .false.
      return
    else if (mmpi_myid == 0) then
      if (bgStats%numvar3d > 0) then
        write(*,*) 'bcsc_setupCH: Number of 3D variables', &
	  bgStats%numvar3d,bgStats%varNameList(1:bgStats%numvar3d)
      end if
      if (bgStats%numvar2d > 0) then
        write(*,*) 'bcsc_setupCH: Number of 2D variables', &
	  bgStats%numvar2d,bgStats%varNameList(bgStats%numvar3d+1: &
	    bgStats%numvar3d+bgStats%numvar2d)
      end if
    end if
    
    bgStats%nkgdim = &
      bgStats%nsposit(bgStats%numvar3d+bgStats%numvar2d+1)-1

    ! Scalefactors must be > 0 until ensembles for constituents can be used.
    
    if (bgStats%numvar3d > 0) then    
      if (any(scaleFactor(1:bgStats%numvar3d,1:nLev_T) <= 0.0D0)) then
        write(*,*) 'Scalefactors: ',scaleFactor(1:bgStats%numvar3d,1:nLev_T) 
        call utl_abort('bcsc_setupCH: Scalefactors values must be > 0 for now.')
      end if
    end if
    if (bgStats%numvar2d > 0) then     
      if (any(scaleFactor(1:bgStats%numvar2d,1) <= 0.0D0)) then 
        write(*,*) 'Scalefactors: ',scaleFactor(1:bgStats%numvar3d,1:nLev_T) 
        call utl_abort('bcsc_setupCH: Scalefactors values must be > 0 for now.') 
      end if
    end if
    
    ! Set scalefactor_stddev

    scaleFactor_stddev(1:bgStats%numvar3d+bgStats%numvar2d, &
      1:max(nLev_M,nLev_T))=0.0d0
      
    where (scaleFactor(1:bgStats%numvar3d+bgStats%numvar2d, &
                       1:max(nLev_M,nLev_T)) > 0.0d0) 
      
      scaleFactor_stddev(1:bgStats%numvar3d+bgStats%numvar2d, &
        1:max(nLev_M,nLev_T)) = sqrt(scaleFactor(1:bgStats%numvar3d &
	                             + bgStats%numvar2d,1:max(nLev_M,nLev_T)))
	  
    end where
    
    bgStats%ni = hco_in%ni
    bgStats%nj = hco_in%nj
    bgStats%ntrunc = ntrunc
    gstID  = gst_setup(bgStats%ni,bgStats%nj,bgStats%ntrunc, &
      bgStats%nkgdim)

    if (allocated(bgStats%lat)) deallocate(bgStats%lat)
    if (allocated(bgStats%lon)) deallocate(bgStats%lon)    
    allocate(bgStats%lat(bgStats%nj))  
    allocate(bgStats%lon(bgStats%ni+1))
    
    bgStats%lat(1:bgStats%nj) = hco_in%lat(1:bgStats%nj)
    bgStats%lon(1:bgStats%ni) = hco_in%lon(1:bgStats%ni)
    bgStats%lon(bgStats%ni+1) = 2*MPC_PI_R8
    
    ! Assign sizes and transfer some fields
    
    if ( bgStats%numvar3d > 0 ) then
      bgStats%nlev = nlev_T
    else if (bgStats%numvar2d > 0 ) then
      bgStats%nlev = 1
    else
      bgStats%nlev = 0
    end if
    
    allocate(stddev(bgStats%ni+1, bgStats%nj, bgStats%nkgdim))
    allocate(rstddev(bgStats%nkgdim, 0:bgStats%ntrunc))    
      
    allocate(bgStats%corns(bgStats%nkgdim,bgStats%nkgdim, &
      0:bgStats%ntrunc))  
    allocate(bgStats%stddev(bgStats%ni+1,bgStats%nj, &
      bgStats%nkgdim ))
    if ( bgStats%nlev > 1 ) then 
      allocate(bgStats%vlev(bgStats%nlev))
      if (getPhysSpaceStats .or. trim(bcsc_mode) == 'BackgroundCheck') then
        allocate(bgStats%corvert(bgStats%nlev,bgStats%nlev, &
          bgStats%numvar3d))
      end if
      if (getPhysSpaceStats) then
        allocate(bgStats%corverti(bgStats%nlev,bgStats%nlev,&
          bgStats%numvar3d))
        allocate(bgStats%invsum(bgStats%nlev,bgStats%numvar3d))
      end if
    end if
    
    ! Get vertical levels
    
    if (bgStats%nlev > 1) then
      call czp_fetchProfile(vco_anl, zps, profT_opt=pressureProfile_T)
      bgStats%vlev(1:bgStats%nlev) = pressureProfile_T(1:bgStats%nlev) 
    else if (bgStats%nlev == 1) then
      bgStats%vlev(1)=zps     
    end if
      
    ! Read covar stats, scale standard deviations,  and apply localization 
    ! to vertical correlation matrices in horizontal spectral space
    
    call bcsc_rdstats
    bgStats%stddev(:,:,:)=stddev(:,:,:)
    if (allocated(stddev)) deallocate(stddev)
    
    ! Generate or read correlation matrix square roots
    
    call bcsc_sucorns2
    
    if(mmpi_myid == 0) write(*,*) 'Completed bcsc_setupCH'
    
    bgStats%initialized = .true.
    if (allocated(rstddev)) deallocate(rstddev)

  end subroutine bcsc_setupCH

  !--------------------------------------------------------------------------
  ! bcsc_getScaleFactor
  !--------------------------------------------------------------------------
  subroutine bcsc_getScaleFactor(scaleFactorOut)
    !
    ! :Purpose: To set scaling factors for background error std. dev.
    !
    implicit none

    !Arguments
    real(8), intent(out) :: scaleFactorOut(:,:) ! Error std. dev. scale factor

    !Locals
    integer :: levelIndex,varIndex

    do varIndex = 1, bgStats%numvar3d+bgStats%numvar2d
    do levelIndex = 1, bgStats%nsposit(varIndex+1)-bgStats%nsposit(varIndex)
      scaleFactorOut(varIndex,levelIndex) = scaleFactor_stddev(varIndex,levelIndex)
    end do
    end do

  end subroutine bcsc_getScaleFactor

  !--------------------------------------------------------------------------
  ! bcsc_rdstats
  !--------------------------------------------------------------------------
  subroutine bcsc_rdstats
    !
    ! :Purpose: To read chemical constituents background stats file.
    !
    implicit none

    !Locals
    integer :: ierr, fnom, fstouv, fstfrm, fclos
    logical :: lExists
    character(len=12) :: bFileName = './bgchemcov'
    type(struct_vco),pointer :: vco_file => null()

    inquire(file=bFileName,exist=lExists)
    if ( lexists ) then
      ierr = fnom(nulbgst,bFileName,'RND+OLD+R/O',0)
      if ( ierr == 0 ) then
        ierr =  fstouv(nulbgst,'RND+OLD')
      else
        call utl_abort('bcsc_rdstats: Problem in opening the background ' // &
	               'chemical constituent stat file')
      end if
    else
      call utl_abort('bcsc_rdstats: Background chemical constituent stat ' // &
                     'file is missing')
    end if

    ! check if analysisgrid and covariance file have the same vertical levels
    call vco_SetupFromFile( vco_file,  & ! OUT
                            bFileName )  ! IN
    if (.not. vco_equal(vco_anl,vco_file)) then
      call utl_abort('bcsc_rstats: vco from analysisgrid and chem cov file ' // &
                     'do not match')
    end if

    ! Read spectral space correlations
    
    call bcsc_readcorns2
    
    ! Read error standard deviations
    
    call bcsc_rdstddev 

    ! Scale error standard deviations (and save to file if requested)
    
    call bcsc_scalestd
    
    ierr = fstfrm(nulbgst)
    ierr = fclos(nulbgst)

  end subroutine bcsc_rdstats

  !--------------------------------------------------------------------------
  ! bcsc_scalestd
  !--------------------------------------------------------------------------
  subroutine bcsc_scalestd
    !
    ! :Purpose: To scale error standard-deviation values.
    !
    implicit none

    !Locals
    integer :: lonIndex, latIndex, varIndex, levelIndex, nlev, nulsig
    integer :: ierr, fnom, fstouv, fstfrm, fclos
  
    if (WritePhysSpaceStats .and. mmpi_myid == 0) then
      nulsig = 0
      ierr = fnom(nulsig,'bCovarSetupChem_out.fst','STD+RND',0)
      ierr = fstouv(nulsig,'RND')
      ierr = utl_fstecr(bgStats%vlev,-32,nulsig,0,0,0,1,1,bgStats%nlev, &
                        0,0,0,'X','PX','Pressure','X',0,0,0,0,5,.true.)
      ierr = utl_fstecr(bgStats%lat,-32,nulsig,0,0,0,1,bgStats%nj,1, &
                        0,0,0,'X','^^','latitude','X',0,0,0,0,5,.true.)
      ierr = utl_fstecr(bgStats%lon,-32,nulsig,0,0,0,bgStats%ni+1,1,1, &
                        0,0,0,'X','>>','longitude','X',0,0,0,0,5,.true.)
    end if
    
    do varIndex = 1,bgStats%numvar3d+bgStats%numvar2d
      nlev=bgStats%nsposit(varIndex+1)-bgStats%nsposit(varIndex)
      do lonIndex = 1, bgStats%ni+1
      do latIndex = 1, bgStats%nj
        stddev(lonIndex,latIndex,bgStats%nsposit(varIndex): &
	       bgStats%nsposit(varIndex+1)-1) = &
               scaleFactor_stddev(varIndex,1:nlev)* &
               stddev(lonIndex,latIndex,bgStats%nsposit(varIndex): &
	       bgStats%nsposit(varIndex+1)-1)
      end do
      end do

      if (WritePhysSpaceStats .and. mmpi_myid == 0) then
        do levelIndex=1,nlev
          ierr = utl_fstecr(stddev(1:bgStats%ni+1,1:bgStats%nj, &
	         bgStats%nsposit(varIndex)-1+levelIndex),-32,nulsig, &
		 0,0,0,bgStats%ni+1,bgStats%nj,1,levelIndex,0,nlev, &
                 'X',bgStats%varNameList(varIndex),'STDDEV','X',0,0,0,0, &
		 5,.true.)
        end do
      end if

    end do
    
    if (WritePhysSpaceStats .and. mmpi_myid == 0) then
      ierr = fstfrm(nulsig)  
      ierr = fclos(nulsig)
    end if

  end subroutine bcsc_scalestd
  
  !--------------------------------------------------------------------------
  ! bcsc_readCorns2
  !--------------------------------------------------------------------------
  subroutine bcsc_readCorns2
    !
    ! :Purpose: To read correlation information and to form the correlation
    !          matrix.
    !
    !:Notes: Can read distinct block diagonal matrices with or without
    !        cross-correlations.
    !

    ! Based on bhi_readcorns2.
    implicit none

    !Locals
    integer :: jn,ierr,varIndex,varIndex2
    integer :: jcol,jrow,jstart,jnum,jstart2,jnum2
    real(8), allocatable, dimension(:) :: zstdsrc
    real(8), allocatable, dimension(:,:) :: zcornssrc

    ! Standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3,idateo
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=4), allocatable :: clnomvarCrosscorns(:)
    character(len=12) :: cletiket
    integer :: fstinf

    rstddev(:,:) = 0.0d0
    bgStats%corns(:,:,:) = 0.0d0
    if (any(CrossCornsVarKindCH(:) /= '')) then
      allocate(clnomvarCrosscorns(bgStats%numvar3d+bgStats%numvar2d))
      clnomvarCrosscorns(:)=''
    end if
    
    ! Read auto-correlations matrices
    
    do varIndex = 1, bgStats%numvar3d+bgStats%numvar2d
   
      clnomvar = bgStats%varNameList(varIndex)
      jstart = bgStats%nsposit(varIndex)
      jnum = bgStats%nsposit(varIndex+1)-bgStats%nsposit(varIndex)
      allocate(zcornssrc(jnum,jnum))
      allocate(zstdsrc(jnum))

      do jn = 0, bgStats%ntrunc

        ! Looking for FST record parameters.
      
        idateo = -1
        cletiket = 'RSTDDEV'
        ip1 = -1
        ip2 = jn
        ip3 = -1
        cltypvar = 'X'
          
        ierr = utl_fstlir(ZSTDSRC,nulbgst,INI,INJ,INK,idateo,cletiket,ip1, &
	                  ip2,ip3,cltypvar,clnomvar)
          
        if (ierr < 0 .and. ip2 < 10 .and. all(CrossCornsVarKindCH(:) == '')) then
          write(*,*) 'bcsc_readcorns2: RSTDDEV ',ip2,jnum,clnomvar
          call utl_abort('bcsc_readcorns2: Problem with constituents ' // &
	                 'background stat file')
        else if (ierr < 0 .and. ip2 == 0 .and. any(CrossCornsVarKindCH(:) /= '')) then
          write(*,*) 'bcsc_readcorns2: Assumes content from cross-corrns ' // &
	             'input for ',clnomvar
          clnomvarCrosscorns(varIndex)=clnomvar
          exit
        end if
        if (ini /= jnum)  then
	  call utl_abort('bcsc_readcorns2: Constituents ' // &
	                 'background stat levels inconsistencies')
        end if

        ! Looking for FST record parameters..

        if (ierr >= 0) then
          idateo = -1
          cletiket = 'CORRNS'
          ip1 = -1
          ip2 = jn
          ip3 = -1
          cltypvar = 'X'
          ierr = utl_fstlir(ZCORNSSRC,nulbgst,INI,INJ,INK,idateo,cletiket, &
	                    ip1,ip2,ip3,cltypvar,clnomvar)

          if (ierr < 0) then
            write(*,*) 'bcsc_readcorns2: CORRNS ',ip2,jnum,clnomvar
            call utl_abort('bcsc_readcorns2: Problem with constituents ' // &
	                   'background stat file')
          end if
          if (ini /= jnum .and. inj /= jnum) then
	    call utl_abort('bcsc_readcorns2: Constituents BG stat levels ' // &
	                   'inconsistencies')
	  end if
        else
          write(*,*) 'WARNING from bcsc_readcorns2: Setting RSDTDDEV to ' // &
	             '1.D-15 for NOMVAR and JN: ',clnomvar,' ',jn
          zstdsrc(1:jnum) = 1.D-15
          zcornssrc(1:jnum,1:jnum) = 0.0D0
          do jrow = 1, jnum
            zcornssrc(jrow,jrow) = 1.0D0
          end do
        end if
          
        rstddev(jstart:jstart+jnum-1,jn) = zstdsrc(1:jnum)
        bgStats%corns(jstart:jstart+jnum-1,jstart:jstart+jnum-1,jn)= &
	  zcornssrc(1:jnum,1:jnum)
          
      end do
       
      deallocate(zcornssrc)
      deallocate(zstdsrc)

    end do

    ! Read cross-correlation matrices
    
    if (any(CrossCornsVarKindCH(:) /= '')) then  
     
      do varIndex = 1, bgStats%numvar3d+bgStats%numvar2d
   
        if (all(CrossCornsVarKindCH(:) /= bgStats%varNameList(varIndex))) cycle

        clnomvar = bgStats%varNameList(varIndex)
        jstart = bgStats%nsposit(varIndex)
        jnum = bgStats%nsposit(varIndex+1)-bgStats%nsposit(varIndex)

        if (clnomvarCrosscorns(varIndex) == clnomvar) then 
	  clnomvarCrosscorns(varIndex)=''
	end if
        
        do varIndex2 = 1, bgStats%numvar3d+bgStats%numvar2d
        
          if (varIndex == varIndex2) cycle
          
          cletiket='CORRNS '//bgStats%varNameList(varIndex2)
          ierr = fstinf(nulbgst,INI,INJ,INK,-1,cletiket,-1,-1,-1,'X',clnomvar)
          if (ierr < 0 ) cycle
          
          jstart2 = bgStats%nsposit(varIndex2)
          jnum2 =  bgStats%nsposit(varIndex2+1)-bgStats%nsposit(varIndex2)
          allocate(zcornssrc(jnum,jnum2))
          
          if (clnomvarCrosscorns(varIndex2) == bgStats%varNameList(varIndex2)) then
	    clnomvarCrosscorns(varIndex2)=''
	  end if
          
          do jn = 0, bgStats%ntrunc
 
            ierr = utl_fstlir(ZCORNSSRC,nulbgst,INI,INJ,INK,-1,cletiket,-1, &
	                      jn,-1,'X',clnomvar)

            if (ierr < 0) then
              if (jn < 10) then
                write(*,*) 'bcsc_readcorns2: CORRNS ',ip2,jnum,clnomvar, &
                bgStats%varNameList(varIndex2)
                call utl_abort('bcsc_readcorns2: Problem with constituents ' // &
                               'background stat file')
              else
                exit
              end if
            end if
            if (ini /= jnum .and. inj /= jnum2) then
	      call utl_abort('bcsc_readcorns2: Constituents BG2 stat ' // &
	                     'levels inconsistencies')
            end if
	     
            bgStats%corns(jstart:jstart+jnum-1,jstart2:jstart2+jnum2-1,jn)= &
	      zcornssrc(1:jnum,1:jnum2)
            bgStats%corns(jstart2:jstart2+jnum2-1,jstart:jstart+jnum-1,jn)= &
	      transpose(zcornssrc(1:jnum,1:jnum2))
          
          end do
          deallocate(zcornssrc)
        end do
      end do   
    end if
    
    if (any(CrossCornsVarKindCH(:) /= '')) then
      if (any(clnomvarCrosscorns(:) /= '')) then
        write(*,*) 'bcsc_readcorns2: Missing matrix for ', &
	           clnomvarCrosscorns(1:bgStats%numvar3d+bgStats%numvar2d)
        call utl_abort('bcsc_readcorns2: Missing correlations matrix')      
      end if
      deallocate(clnomvarCrosscorns)
    end if
     
    ! Apply convolution to RSTDDEV correlation

    call bcsc_convol
    
    do jn = 0, bgStats%ntrunc

      ! Re-build correlation matrix: factorization of corns with convoluted RSTDDEV
      do jcol = 1, bgStats%nkgdim
        bgStats%corns(1:bgStats%nkgdim,jcol,jn) = &
	  rstddev(1:bgStats%nkgdim,jn) * &
	  bgStats%corns(1:bgStats%nkgdim,jcol,jn) * rstddev(jcol,jn)
      end do

    end do

  end subroutine bcsc_readCorns2

  !--------------------------------------------------------------------------
  ! bcsc_convol
  !--------------------------------------------------------------------------
  subroutine bcsc_convol
    implicit none

    !Locals
    real(8) :: dlfact2,dlc,dsummed
    real(8) :: dtlen,zr,dlfact
    integer :: jn,latIndex,jk,varIndex,levelIndex,nsize,ierr
    real(8) :: zleg(0:bgStats%ntrunc,bgStats%nj)
    real(8) :: zsp(0:bgStats%ntrunc,bgStats%nkgdim)
    real(8) :: zgr(bgStats%nj,bgStats%nkgdim)
    real(8) :: zrmu(bgStats%nj)

    integer :: nlev_MT,ini,inj,ink,nulcorns
    real(8), allocatable :: wtemp(:,:,:)    
    real(8), allocatable :: hcorrel(:,:,:),hdist(:)
    logical :: lfound
    integer :: fnom, fstouv, fstfrm, fclos

    do latIndex = 1, bgStats%nj
      zrmu(latIndex)  = gst_getrmu(latIndex,gstID)
    end do

    ! CONVERT THE CORRELATIONS IN SPECTRAL SPACE INTO SPECTRAL
    ! COEFFICIENTS OF THE CORRELATION FUNCTION AND FUNCTION TO BE
    ! SELF-CONVOLVED

    do jn = 0, bgStats%ntrunc
      dlfact = ((2.0d0*jn+1.0d0)/2.0d0)**(0.25d0)
      dlfact2 = ((2.0d0*jn +1.0d0)/2.0d0)**(0.25d0)
      do jk = 1, bgStats%nkgdim
        zsp(jn,jk) = rstddev(jk,jn)*dlfact*dlfact2
      end do
    end do

    ! Transform to physical space
    call gst_zleginv(gstID,zgr,zsp,bgStats%nkgdim)
    
    ! Truncate in horizontal extent with Gaussian window
    
    varIndex=1
    do jk = 1, bgStats%nkgdim
      if (jk == bgStats%nsposit(varIndex)) then
        dtlen = rpor(varIndex)
        varIndex=varIndex+1 
      end if
      if (dtlen > 0.0d0) then
        dlc = 1.d0/dble(dtlen)
        dlc = 0.5d0*dlc*dlc
        do latIndex = 1, bgStats%nj
          zr = ec_ra * acos(zrmu(latIndex))
          dlfact = dexp(-(zr**2)*dlc)
          zgr(latIndex,jk) = dlfact*zgr(latIndex,jk)
        end do
      end if

      !write(*,*) 'zeroing length (km)=',jk,dtlen/1000.0
    end do

    ! Transform back to spectral space
    call gst_zlegdir(gstID,zgr,zsp,bgStats%nkgdim)

    ! Convert back to correlations
    do jk = 1, bgStats%nkgdim
      do jn = 0, bgStats%ntrunc
        zsp(jn,jk) = zsp(jn,jk)*(2.0d0/(2.0d0*jn+1.0d0))**(0.25d0)
      end do
    end do

    ! PUT BACK INTO RSTDDEV
    do jn = 0, bgStats%ntrunc
      do jk = 1, bgStats%nkgdim
        rstddev(jk,jn) = zsp(jn,jk)
      end do
    end do

    ! Re-normalize to ensure correlations
    do jk = 1, bgStats%nkgdim
      dsummed = 0.d0
      do jn = 0, bgStats%ntrunc
        dsummed = dsummed+ dble(rstddev(jk,jn)**2)*sqrt(((2.d0*jn)+1.d0)/2.d0)
      end do
      dsummed = sqrt(dsummed)
      do jn = 0, bgStats%ntrunc
        if(dsummed > 1.d-30) rstddev(jk,jn) = rstddev(jk,jn)/dsummed
      end do
    end do

    ! CONVERT THE SPECTRAL COEFFICIENTS OF THE CORRELATION FUNCTION
    ! BACK INTO CORRELATIONS OF SPECTRAL COMPONENTS
    do jn = 0, bgStats%ntrunc
      dlfact = sqrt(0.5d0)*(1.0d0/((2.0d0*jn+1)/2.0d0))**0.25d0
      do jk = 1, bgStats%nkgdim
        rstddev(jk,jn) = rstddev(jk,jn)*dlfact
      end do
    end do

    if ( .not.getPhysSpaceHCorrel .or. .not.getPhysSpaceStats ) return

    ! Compute resultant physical space horizontal correlations and
    ! 1/e correlation length from correlation array if not available
    
    if (allocated(hcorrel)) deallocate(hcorrel)
    allocate(hcorrel(bgStats%nj,bgStats%nlev, &
      bgStats%numvar3d+bgStats%numvar2d))
    if (allocated(wtemp)) deallocate(wtemp)
    allocate(wtemp(0:bgStats%ntrunc,bgStats%nj,1))
    if (allocated(bgStats%hcorrlen)) deallocate(bgStats%hcorrlen)
    allocate(bgStats%hcorrlen(bgStats%nlev,  &
      bgStats%numvar3d+bgStats%numvar2d))
    if (allocated(hdist)) deallocate(hdist)
    allocate(hdist(bgStats%nj))

    lfound=.true.
    do varIndex = 1, bgStats%numvar3d+bgStats%numvar2d
      nlev_MT = bgStats%nsposit(varIndex+1)-bgStats%nsposit(varIndex)
      do levelIndex = 1, nlev_MT
        ierr = utl_fstlir(hcorrel(:,levelIndex,varIndex),nulbgst,INI,INJ, &
	                  INK,-1,'HCORREL',levelIndex,-1,-1,'X', &
	                  bgStats%varNameList(varIndex))
        if (ierr < 0 ) then
          lfound=.false.
          exit
        end if
      end do
      ierr = utl_fstlir(bgStats%hcorrlen(1:nlev_MT,varIndex),nulbgst,INI, &
             INJ,INK,-1,'HCORRLEN',-1,-1,-1,'X',bgStats%varNameList(varIndex))
      if (ierr < 0 ) then
        lfound=.false.
        exit
      end if
    end do 
    
    if (lfound) return

    do latIndex = 1, bgStats%nj
      hdist(latIndex)=ec_ra*acos(zrmu(latIndex))
    end do

    zleg(:,:)=0.0d0
    wtemp(:,:,:)=0.0d0
    hcorrel(:,:,:)=0.0d0
    bgStats%hcorrlen(:,:)=0.0
    
    do latIndex = mmpi_myid+1, bgStats%nJ, mmpi_nprocs
      do jn = 0, bgStats%ntrunc
        wtemp(jn,latIndex,1) = gst_getzleg(jn,latIndex,gstID)
      end do
    end do
    
    nsize=bgStats%nJ*(bgStats%ntrunc+1)    
    call rpn_comm_allreduce(wtemp(0:bgStats%ntrunc,1:bgStats%ni,1), &
         zleg(0:bgStats%ntrunc,1:bgStats%nJ),nsize,"mpi_double_precision", &
         "mpi_sum","GRID",ierr)

    deallocate(wtemp)
    allocate(wtemp(bgStats%nj, bgStats%nlev, &
      bgStats%numvar3d+bgStats%numvar2d))
    wtemp(:,:,:)=0.0
    
    varIndex = 1
    levelIndex = 1
    do jk = 1, bgStats%nkgdim
      if (jk == bgStats%nsposit(varIndex+1)) then
        varIndex = varIndex+1 
        levelIndex = 1
      end if

      do latIndex = mmpi_myid+1, bgStats%nj, mmpi_nprocs
        do jn = 0, bgStats%ntrunc
          wtemp(latIndex,levelIndex,varIndex) = wtemp(latIndex,levelIndex, &
	        varIndex)+rstddev(jk,jn)*rstddev(jk,jn)*  &
                sqrt(2.0)*sqrt(2.0*jn+1.0)*zleg(jn,latIndex)
        end do       
      end do
      levelIndex = levelIndex+1
    end do
    
    nsize=bgStats%nj*bgStats%nkgdim   
    call rpn_comm_allreduce(wtemp,hcorrel,nsize,"mpi_double_precision", &
                            "mpi_sum","GRID",ierr)
    deallocate(wtemp)
    
    if ( mmpi_myid == 0 ) then

      varIndex = 1
      levelIndex = 1
      do jk = 1, bgStats%nkgdim
        if (jk == bgStats%nsposit(varIndex+1)) then
          varIndex = varIndex+1 
          levelIndex = 1
        end if
        do latIndex=bgStats%nj-1,2,-1
          if (hcorrel(latIndex,levelIndex,varIndex) <= 0.368) then  ! 1/e ~ 0.368
            bgStats%hcorrlen(levelIndex,varIndex) = (hdist(latIndex)* &
	      (hcorrel(latIndex+1,levelIndex,varIndex)-0.368) &
              + hdist(latIndex+1)*(0.368-hcorrel(latIndex,    &
	      levelIndex,varIndex))) &
              /(hcorrel(latIndex+1,levelIndex,varIndex)- &
	      hcorrel(latIndex,levelIndex,varIndex))
            exit
          end if
        end do
        levelIndex = levelIndex+1
      end do  
    
      if (WritePhysSpaceStats) then 
        nulcorns = 0
        ierr = fnom(nulcorns,'bCovarSetupChem_out.fst','STD+RND',0)
        ierr = fstouv(nulcorns,'RND')

        do varIndex = 1, bgStats%numvar3d+bgStats%numvar2d
          nlev_MT = bgStats%nsposit(varIndex+1)-bgStats%nsposit(varIndex)
          do levelIndex = 1, nlev_MT
            ierr = utl_fstecr(hcorrel(:,levelIndex,varIndex),-32,nulcorns,0,0,0, &
	                      1,bgStats%nj,1,levelIndex,0,nlev_MT,'X', &
		              bgStats%varNameList(varIndex), &
                              'HCORREL','X',0,0,0,0,5,.true.)
          end do
          ierr = utl_fstecr(bgStats%hcorrlen(1:nlev_MT,varIndex),-32,nulcorns, &
	         0,0,0,1,1,nlev_MT,0,0,0,'X',bgStats%varNameList(varIndex), &
                 'HCORRLEN','X',0,0,0,0,5,.true.)
          ierr = utl_fstecr(hdist(1:bgStats%nj),-32,nulcorns,0,0,0,1, &
	         bgStats%nj,1,0,0,0,'X',bgStats%varNameList(varIndex), &
                 'HDIST','X',0,0,0,0,5,.true.)
        end do
      
        ierr = fstfrm(nulcorns)  
        ierr = fclos(nulcorns)
	
      end if
      
      write(*,*)
      write(*,*) 'bcsc_convol: Horizontal correlations'
      write(*,*)
      write(*,*) 'Separation distances (km)'
      write(*,*) bgStats%nj-bgStats%nj*4/5+1
      write(*,*) hdist(bgStats%nj*4/5:bgStats%nj)/1000.00
      write(*,*)
      do varIndex = 1, bgStats%numvar3d+bgStats%numvar2d
        if (varIndex <= bgStats%numvar3d) then
          write(*,*) bgStats%varNameList(varIndex), bgStats%nlev
          do levelIndex = 1, bgStats%nlev
            write(*,'(i3,f10.2,3000f6.2)') levelIndex, &
	                  bgStats%hcorrlen(levelIndex,varIndex)/1000.00, &
	                  hcorrel(bgStats%nj*4/5:bgStats%nj,levelIndex,varIndex)
          end do
        else
          write(*,*) bgStats%varNameList(varIndex), 1
          write(*,'(i3,f10.2,3000f6.2)') 1, &
	                           bgStats%hcorrlen(1,varIndex)/1000.00, &
	                           hcorrel(bgStats%nj*4/5:bgStats%nj,1,varIndex)
        end if
        write(*,*)
      end do
    end if

    if (allocated(hcorrel)) deallocate(hcorrel)
    if (allocated(hdist)) deallocate(hdist)

  end subroutine bcsc_convol

  !--------------------------------------------------------------------------
  ! bcsc_rdstddev
  !--------------------------------------------------------------------------
  subroutine bcsc_rdstddev
    !
    ! :Purpose: To read stddev and to set as 3D fields.
    !
    implicit none

    !Locals
    integer :: ikey
    real(8) :: rcoord(10000)
    
    ! standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3
    integer :: idate(100)

    ! Reading the data

    idate(1) = -1
    ip2      = -1
    ip3      = -1

    ! Get latitudes and longitudes if available

    ip1=-1
    ikey = utl_fstlir(rcoord,nulbgst,ini,inj,ink, &
                      idate(1),' ',ip1,ip2,ip3,' ','^^')
                         
    if (ikey >= 0) then
      if (allocated(rlatr)) deallocate(rlatr)
      allocate(rlatr(inj))
      rlatr(1:inj) = rcoord(1:inj) 
      rlatr(1:inj) = rlatr(1:inj)*MPC_RADIANS_PER_DEGREE_R8
    else 
      ! Assume same as bgStats%lat
      if (allocated(rlatr)) deallocate(rlatr)
      allocate(rlatr(bgStats%nj))
      inj = bgStats%nj
      rlatr(1:inj) = bgStats%lat(1:inj) 
    end if    

    ikey = utl_fstlir(rcoord,nulbgst,ini,inj,ink, &
                      idate(1),' ',ip1,ip2,ip3,' ','>>')
                         
    if (ikey >= 0) then
      if (allocated(rlongr)) deallocate(rlongr)
      allocate(rlongr(ini+1))
      rlongr(1:ini) = rcoord(1:ini) 
      rlongr(1:ini) = rlongr(1:ini)*MPC_RADIANS_PER_DEGREE_R8
    else if (stddevMode /= 'SP2D') then
      ! Assume same as bgStats%lon
      if (allocated(rlongr)) deallocate(rlongr)
      allocate(rlongr(bgStats%ni+1))
      ini = bgStats%ni
      rlongr(1:ini) = bgStats%lon(1:ini) 
    end if 
    rlongr(ini+1) = 360.*MPC_RADIANS_PER_DEGREE_R8
     
    ! Read specified input type for error std. dev.
    
    if(stddevMode == 'GD3D') then
      call bcsc_rdstd3D
    elseif(stddevMode == 'GD2D') then
      call bcsc_rdstd
    elseif(stddevMode == 'SP2D') then
      call bcsc_rdspstd
    else
      call utl_abort('bcsc_rdstddev: unknown stddevMode')
    end if
    
  end subroutine bcsc_rdstddev

  !--------------------------------------------------------------------------
  ! bcsc_rdspstd
  !--------------------------------------------------------------------------
  subroutine bcsc_rdspstd
    implicit none

    !Locals
    integer :: varIndex,jn,inix,injx,inkx
    integer :: ikey, levelIndexo, firstn,lastn
    real(8) :: zsp(0:bgStats%ntrunc,max(nlev_M,nlev_T))
    real(8) :: zspbuf(max(nlev_M,nlev_T))
    real(8) :: zgr(bgStats%nj,max(nlev_M,nlev_T))
   
    ! standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3
    integer :: idate(100),nlev_MT
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fstinf

    stddev(:,:,:) = 0.0d0
    
    ! Reading the Legendre poly representation of the 2D background 
    ! error std. dev. field

    idate(1) = -1
    ip1      = -1
    ip2      = -1
    ip3      = -1

    cletiket = 'SPSTDDEV'
    cltypvar = 'X'
    
    do varIndex = 1,bgStats%numvar3d+bgStats%numvar2d
      clnomvar = bgStats%varNameList(varIndex)
      nlev_MT = bgStats%nsposit(varIndex+1)-bgStats%nsposit(varIndex)
      
      firstn = -1
      do jn = 0, bgStats%ntrunc
        ip2 = jn
        ikey = fstinf(nulbgst,inix,injx,inkx,idate(1),cletiket,ip1,ip2,ip3, &
	       cltypvar,clnomvar)

        if(ikey >= 0 ) then
          ikey = utl_fstlir(zspbuf(1:nlev_MT),nulbgst,ini,inj,ink,idate(1), &
	                    cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
        else
          if(firstn == -1) firstn = jn
          lastn = jn
          zspbuf(:) = 0.0d0
        end if

        if (ini /= nlev_MT) then
          if (ini == nlev_MT-1) then
            zspbuf(nlev_MT) = zspbuf(ini)
            write(*,*) 'WARNING in CHM_RDSPSTD: ini and nlev_MT not ', &
	               'same size - ',ini,nlev_MT
          else
            write(*,*) 'JN, INI, nlev_MT, NOMVAR: ',jn,ini,nlev_MT,' ',clnomvar
            call utl_abort('CHM_RDSPSTD: Constituents background stats ' // &
	                   'levels inconsistency')
          end if    
        end if
     
        do levelIndexo = 1, nlev_MT
          zsp(jn,levelIndexo) = zspbuf(levelIndexo)
        end do
      end do
      
      if (mmpi_myid == 0 .and. firstn /= -1) then
        write(*,*) 'WARNING: CANNOT FIND SPSTD FOR ',clnomvar, &
                   ' AT N BETWEEN ',firstn,' AND ',lastn,', SETTING TO ZERO!!!'
      end if

      call gst_zleginv(gstID,zgr(:,1:nlev_MT),zsp(:,1:nlev_MT),nlev_MT)

      do jn = 1, bgStats%ni+1
        stddev(jn,1:bgStats%nj, &
	  bgStats%nsposit(varIndex):bgStats%nsposit(varIndex+1)-1) = &
	  zgr(1:bgStats%nj,1:nlev_MT)
      end do

    end do 

  end subroutine bcsc_rdspstd

  !--------------------------------------------------------------------------
  ! bcsc_rdspstd_newfmt
  !--------------------------------------------------------------------------
  subroutine bcsc_rdspstd_newfmt
    implicit none

    !Locals
    integer :: varIndex,jn,inix,injx,inkx,ntrunc_file
    integer :: ikey,levelIndexo
    real(8) :: zsp(0:bgStats%ntrunc,max(nlev_M,nlev_T))
    real(8), allocatable :: zspbuf(:)
    real(8) :: zgr(bgStats%nj,max(nlev_M,nlev_T))

    ! standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3
    integer :: idate(100),nlev_MT
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fstinf

    stddev(:,:,:) = 0.0d0

    ! Reading the data

    idate(1) = -1
    ip2      = -1
    ip3      = -1

    cletiket = 'SPSTDDEV'
    cltypvar = 'X'

    ! Check if file is old format
    
    ip1 = -1
    clnomvar = bgStats%varNameList(1) 
    ikey = fstinf(nulbgst,inix,injx,inkx,idate(1),cletiket,ip1,ip2,ip3, &
                  cltypvar,clnomvar)
    write(*,*) 'ini,inj,ink=',inix,injx,inkx
    if(inix > 1) then
      write(*,*) 'bcsc_rdspstd_newfmt: ini>1, SPSTDDEV is in old format, ', &
                 ' calling bcsc_RDSPSTD...'
      call bcsc_rdspstd
      return
    end if

    ! write(*,*) 'Reading for 3D and 2D variables'
    
    do varIndex = 1,bgStats%numvar3d+bgStats%numvar2d
      clnomvar = bgStats%varNameList(varIndex)
      nlev_MT = bgStats%nsposit(varIndex+1)-bgStats%nsposit(varIndex)

      !write(*,*) 'Reading ',clnomvar

      do levelIndexo = 1, nlev_MT
        if (nlev_MT == 1) then
           ip1 = -1
        else if(vnl_varLevelFromVarName(clnomvar) == 'MM') then
          ip1 = vco_anl%ip1_M(levelIndexo)
        else
          ip1 = vco_anl%ip1_T(levelIndexo)
        end if
        
        ikey = fstinf(nulbgst,inix,ntrunc_file,inkx,idate(1),cletiket, &
	       ip1,ip2,ip3,cltypvar,clnomvar)
        ntrunc_file = ntrunc_file-1

        allocate(zspbuf(0:ntrunc_file))
        
        if(ikey >= 0 ) then
          ikey = utl_fstlir(zspbuf(0:ntrunc_file),nulbgst,ini,inj,ink, &
	                    idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
        else
          write(*,*) 'bcsc_rdspstd_newfmt: ',varIndex,clnomvar,nlev_MT, &
	             levelIndexo,ikey,bgStats%ntrunc,ntrunc_file
          call utl_abort('bcsc_rdspstd_newfmt: SPSTDDEV record not found')
        end if

        zsp(:,levelIndexo) = 0.0d0
        do jn = 0, min(bgStats%ntrunc,ntrunc_file)
          zsp(jn,levelIndexo) = zspbuf(jn)
        end do
        deallocate(zspbuf)
      end do

      call gst_zleginv(gstID,zgr(:,1:nlev_MT),zsp(:,1:nlev_MT),nlev_MT)

      do jn = 1, bgStats%ni+1
        stddev(jn,1:bgStats%nj,bgStats%nsposit(varIndex): &
	  bgStats%nsposit(varIndex+1)-1) = zgr(1:bgStats%nj,1:nlev_MT)
      end do

    end do

  end subroutine bcsc_rdspstd_newfmt

  !--------------------------------------------------------------------------
  ! bcsc_rdstd
  !--------------------------------------------------------------------------
  subroutine bcsc_rdstd
    !
    ! :Purpose: To read 2D stddev and to store as 3D
    !
    implicit none

    !Locals
    integer :: varIndex,in
    integer :: ikey
    real(8), allocatable :: stddev3d(:,:,:)
    
    ! standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3
    integer :: idate(100),nlev_MT
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fstinf
    
    real(8), allocatable  :: vlev(:),vlevout(:)

    stddev(:,:,:) = 0.0d0

    ! Reading the data

    idate(1) = -1
    ip1      = -1
    ip2      = -1
    ip3      = -1
    
    cletiket = 'STDDEV'
    cltypvar=' '

    ! Reading for 3D and 2D variables
    
    do varIndex = 1,bgStats%numvar3d+bgStats%numvar2d
      clnomvar = bgStats%varNameList(varIndex)
      nlev_MT = bgStats%nsposit(varIndex+1)-bgStats%nsposit(varIndex)
 
      ikey = fstinf(nulbgst,ini,inj,ink,idate(1),cletiket,ip1,ip2,ip3, &
                    cltypvar,clnomvar)
      if (ikey < 0 .or. ini > 1 .or. ink /= nlev_MT) then
        write(*,*) 'bcsc_rdstd: ',varIndex,clnomvar,ikey,ini,ink,nlev_MT
        call utl_abort(': bcsc_rdstd record not found or incorrect')          
      end if
      
      allocate(stddev3d(1,inj,ink))
      stddev3d(1,:,:) = 0.0D0
      allocate(vlev(ink),vlevout(nlev_MT))
      vlev(:) = 1.0D0
      vlevout(:) = 1.0D0
      
      ikey = utl_fstlir(stddev3d(1,:,:), nulbgst, ini, inj, ink, &
                        idate(1), cletiket, ip1, ip2, ip3, cltypvar, clnomvar)

      if (ikey < 0) then
        write(*,*) 'bcsc_rdstd: ',varIndex,clnomvar,nlev_MT,ikey
        call utl_abort(': bcsc_rdstd record not found')
      end if
        
      ! Extend to 3D
      if (inj == bgStats%nj) then
        do in = 1, bgStats%ni+1
          stddev(in,:,bgStats%nsposit(varIndex): &
	    bgStats%nsposit(varIndex+1)-1) = stddev3d(1,:,:) 
        end do
      else
        ! Interpolate in lat
        call gsv_field3d_hbilin(stddev3d, 1, inj, ink, rlongr, rlatr, vlev, &
             stddev(:,:,bgStats%nsposit(varIndex): &
	     bgStats%nsposit(varIndex+1)-1), bgStats%ni+1, bgStats%ni, &
	     nlev_MT, bgStats%lon, bgStats%lat, vlevout)
      end if
      deallocate(stddev3d, vlev, vlevout)
       
    end do

  end subroutine bcsc_rdstd

  !--------------------------------------------------------------------------
  ! bcsc_rdstd3d
  !--------------------------------------------------------------------------
  subroutine bcsc_rdstd3d
    !
    ! :Purpose: To read 3D stddev.
    !
    ! Originally based on bcsc_rdspstd_newfmt
    !
    implicit none

    !Locals
    integer :: varIndex
    integer :: ikey,levelIndexo
    real(8), allocatable :: stddev3d(:,:,:)
    
    ! standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3
    integer :: idate(100),nlev_MT
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fstinf

    real(8) :: vlev(1),vlevout(1)

    stddev(:,:,:) = 0.0d0
    vlev(:) = 1.0D0
    vlevout(:) = 1.0D0

    ! Reading the data

    idate(1) = -1
    ip2      = -1
    ip3      = -1
    
    cletiket = 'STDDEV3D'
    cltypvar=' '

    ! Reading for 3D and 2D variables
    
    do varIndex = 1,bgStats%numvar3d+bgStats%numvar2d
      clnomvar = bgStats%varNameList(varIndex)
      nlev_MT = bgStats%nsposit(varIndex+1)-bgStats%nsposit(varIndex)

      !write(*,*) 'Reading ',clnomvar

      do levelIndexo = 1, nlev_MT
        if (nlev_MT == 1) then
          ip1 = -1
        else if(vnl_varLevelFromVarName(clnomvar) == 'MM') then
          ip1 = vco_anl%ip1_M(levelIndexo)
        else
          ip1 = vco_anl%ip1_T(levelIndexo)
        end if
        
        ikey = fstinf(nulbgst,ini,inj,ink,idate(1),cletiket,ip1,ip2,ip3, &
	       cltypvar,clnomvar)
        if (ikey < 0) then
            write(*,*) 'bcsc_rdstd3d: ',varIndex,clnomvar,ikey,levelIndexo
            call utl_abort(': bcsc_RDSTD record not foun0d')          
        end if
      
        allocate(stddev3d(ini+1, inj, 1))
        stddev3d(:,:,1) = 0.0D0
        
        ikey = utl_fstlir(stddev3d(1:ini,:,1), nulbgst, ini, inj, ink, &
                          idate(1), cletiket, ip1, ip2, ip3, cltypvar, clnomvar)
        if (ikey < 0) then
          write(*,*) 'bcsc_rdstd3d: ',varIndex,clnomvar,nlev_MT,levelIndexo, &
	              ikey,ip1
          call utl_abort(': bcsc_RDSTD3D record not found')
        end if
        
        ! Move to stddev
        if (inj == bgStats%nj .and. ini == bgStats%ni) then
          stddev(1:bgStats%ni,:,bgStats%nsposit(varIndex)+(levelIndexo-1)) = &
	    stddev3d(1:bgStats%ni,:,1) 
          stddev(bgStats%ni+1,:,bgStats%nsposit(varIndex)+(levelIndexo-1)) = &
	    stddev3d(1,:,1)
        else
          ! Interpolate in lat and long
          stddev3d(ini+1,:,1) = stddev3d(1,:,1)
          call gsv_field3d_hbilin(stddev3d, ini+1, inj, 1, rlongr, rlatr, vlev, &
               stddev(:,:,bgStats%nsposit(varIndex)+(levelIndexo-1)), &
	       bgStats%ni+1, bgStats%nj, 1, &
               bgStats%lon, bgStats%lat, vlevout)
       end if
       deallocate(stddev3d)

      end do
    end do

  end subroutine bcsc_rdstd3d

  !--------------------------------------------------------------------------
  ! bcsc_sucorns2
  !--------------------------------------------------------------------------
  subroutine bcsc_sucorns2
    implicit none

    !Locals
    real(8) :: eigenval(bgStats%nkgdim)
    real(8) :: eigenvalsqrt(bgStats%nkgdim)
    real(8), allocatable :: eigenvec(:,:),result(:,:)

    integer :: jn,jk1,jk2,varIndex,ierr
    integer :: ilwork,info,jnum,jstart,nsize

    real(8) :: zwork(2*4*bgStats%nkgdim)
    real(8) :: ztlen,zcorr,zr,zpres1,zpres2,eigenvalmax
    real(8), allocatable :: corns_temp(:,:,:)
    logical, allocatable :: lfound_sqrt(:)
    
    ! Apply vertical localization to correlations of 3D fields.
    ! Currently assumes no-cross correlations for variables (block diagonal matrix)
  
    if ( bgStats%numvar3d > 0 ) then   
      do varIndex = 1, bgStats%numvar3d
        ztlen = rvloc(varIndex) ! specify length scale (in units of ln(Pressure))
        
        jstart = bgStats%nsposit(varIndex)
        jnum = bgStats%nsposit(varIndex+1)-bgStats%nsposit(varIndex)
       
        if(ztlen > 0.0d0) then
          ! calculate 5'th order function (from Gaspari and Cohn)
          do jk1 = 1, jnum
            zpres1 = log(bgStats%vlev(jk1))
            do jk2 = 1, jnum
              zpres2 = log(bgStats%vlev(jk2))
              zr = abs(zpres2 - zpres1)
              zcorr = gasparicohn(ztlen,zr)
              bgStats%corns(jstart-1+jk1,jstart-1+jk2,0:bgStats%ntrunc) = &
                bgStats%corns(jstart-1+jk1,jstart-1+jk2, &
	        0:bgStats%ntrunc)*zcorr  
            end do
          end do
        end if
      end do
    end if
    
    ! Compute total vertical correlations and its inverse 
    ! (currently for each block matrix).
    
    if (getPhysSpaceStats .or. trim(bcsc_mode) == 'BackgroundCheck') then
      call bcsc_corvertSetup
    end if
    
    if (trim(bcsc_mode) == 'BackgroundCheck') return

    allocate(lfound_sqrt(bgStats%numvar3d+bgStats%numvar2d))
    if (ReadWrite_sqrt) then
      ! if desired, read precomputed sqrt of corns
      call readcorns(lfound_sqrt,bgStats%ntrunc,'CORNS_SQRT')
    else
      lfound_sqrt(:)=.false.
    end if

    do varIndex = 1, bgStats%numvar3d+bgStats%numvar2d

      if (any(CrossCornsVarKindCH(:) /= '')) then
        jstart=1
        jnum=bgStats%nkgdim
      else
        jstart = bgStats%nsposit(varIndex)
        jnum = bgStats%nsposit(varIndex+1)-bgStats%nsposit(varIndex)
      end if
      
      if (.not.lfound_sqrt(varIndex)) then
        
        ! compute square-root of corns for each total wavenumber
    
        allocate(corns_temp(jnum,jnum,0:bgStats%ntrunc))
        allocate(eigenvec(jnum,jnum),result(jnum,jnum))
        corns_temp(:,:,:) = 0.0d0
        do jn = mmpi_myid, bgStats%ntrunc, mmpi_nprocs

          eigenvec(1:jnum,1:jnum) = bgStats%corns(jstart:jstart-1+jnum, &
	                            jstart:jstart-1+jnum,jn)

          ! CALCULATE EIGENVALUES AND EIGENVECTORS.
          ilwork = 4*jnum*2
          call dsyev('V','U',jnum,eigenvec,jnum,eigenval,zwork,ilwork,info)
          if(info /= 0) then
            write(*,*) 'bcsc_sucorns2: non-zero value of info =',info, &
	               ' returned by dsyev for wavenumber ',jn,varIndex
            call utl_abort('bcsc_sucorns2')
          end if
     
          ! set selected number of eigenmodes to zero
          if(numModeZero > 0) then
            do jk1 = 1, numModeZero
              eigenval(jk1) = 0.0d0
            end do
           end if

           eigenvalmax = maxval(eigenval(1:jnum))
           do jk1 = 1, jnum
             ! if(eigenval(jk1) < 1.0d-15) then
             if(eigenval(jk1) < 1.0d-8*eigenvalmax) then
               eigenvalsqrt(jk1) = 0.0d0
             else
               eigenvalsqrt(jk1) = sqrt(eigenval(jk1))
             end if
           end do

           ! E * lambda^1/2
           do jk1 = 1, jnum
             result(1:jnum,jk1) = eigenvec(1:jnum,jk1)*eigenvalsqrt(jk1)
           end do
  
           ! (E * lambda^1/2) * E^T
           do jk1 = 1, jnum
           do jk2 = 1, jnum
             corns_temp(1:jnum,jk1,jn) = corns_temp(1:jnum,jk1,jn) &
	                               + result(1:jnum,jk2)*eigenvec(jk1,jk2)
           end do
           end do

         end do ! jn
  
         nsize = jnum*jnum*(bgStats%ntrunc+1)
         call rpn_comm_allreduce(corns_temp, &
	   bgStats%corns(jstart:jstart+jnum-1,jstart:jstart+jnum-1,:), &
	   nsize,"mpi_double_precision","mpi_sum","GRID",ierr)
         deallocate(corns_temp,eigenvec,result)

      end if
      
      if (any(CrossCornsVarKindCH(:) /= '')) exit

    end do
   
    deallocate(lfound_sqrt)
    if(ReadWrite_sqrt) then
      ! Write computed sqrt to a separate file.
      call writecorns(bgStats%ntrunc,'CORNS_SQRT',-1)
    end if
    
  end subroutine bcsc_sucorns2

  !--------------------------------------------------------------------------
  ! bcsc_corvertSetup
  !--------------------------------------------------------------------------
  subroutine bcsc_corvertSetup
    !
    ! :Purpose: To compute total vertical correlations (bcsc_corvert) and its
    !           inverse (bgStats%corverti; currently for each block matrix).
    !
    !
    ! :Note: Currently assumes no (or neglects) cross-correlations 
    !
    implicit none

    !Locals
    real(8) :: eigenval(bgStats%nkgdim)
    real(8), allocatable :: eigenvec(:,:),result(:,:)

    integer :: jn,jk1,jk2,varIndex,numvartot,jstart
    integer :: ilwork,info,jnum,nvlev,ierr

    real(8) :: zwork(2*4*bgStats%nkgdim)
    real(8) :: eigenvalmax
    integer iulcorvert
    
    ! Standard file variables
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fnom,fstouv,fstfrm,fclos
    logical, allocatable :: lfound(:)
   
    ! Compute total vertical correlations and its inverse 
    ! (currently for each block matrix).

    if (mmpi_myid == 0 .and. WritePhysSpaceStats) then
      iulcorvert = 0
      ierr = fnom(iulcorvert,'bCovarSetupChem_out.fst','STD+RND',0)
      ierr = fstouv(iulcorvert,'RND')
    end if 
    
    nvlev=-1
    numvartot=bgStats%numvar3d
    do varIndex = 1, numvartot
   
      jstart = bgStats%nsposit(varIndex)
      jnum = bgStats%nsposit(varIndex+1)-bgStats%nsposit(varIndex)
       
      bgStats%corvert(1:jnum,1:jnum,varIndex) = 0.0d0
      do jn = 0, bgStats%ntrunc
        bgStats%corvert(1:jnum,1:jnum,varIndex) = &
	  bgStats%corvert(1:jnum,1:jnum,varIndex)+ ((2*jn+1)* &
          bgStats%corns(bgStats%nsposit(varIndex): &
	  bgStats%nsposit(varIndex+1)-1,&
	  bgStats%nsposit(varIndex):bgStats%nsposit(varIndex+1)-1,jn))
      end do
       
      ! Inverse (and possible vertical interpolation to model levels) 
      ! not needed if not in analysis mode
      
      if (trim(bcsc_mode) == 'BackgroundCheck') cycle

      ! Generate 1/sum(covert(:,i,VarIndex)
        
      do jn = 1, jnum 
        bgStats%invsum(jn,varIndex) = 1.0d0/sum(bgStats%corvert(1:jnum,jn,varIndex)) 
      end do
      
      if (varIndex == 1) then
        allocate(lfound(numvartot))
        call readcorns(lfound,0,'CORVERTI') 
      end if
      
      if (.not.lfound(varIndex)) then
      
        write(*,*) 'bcsc_corvertSetup: Calculating CORVERTI ', &
	           'for varIndex =',varIndex

        allocate(eigenvec(jnum,jnum),result(jnum,jnum))
        eigenvec(1:jnum,1:jnum)=bgStats%corvert(1:jnum,1:jnum,varIndex)       

        ! CALCULATE EIGENVALUES AND EIGENVECTORS.
        ilwork = 4*jnum*2
        call dsyev('V','U',jnum,eigenvec,jnum,eigenval,zwork,ilwork,info)
        if (info /= 0) then
          write(*,*) 'bcsc_corvertSetup: non-zero value of info =',info, &
	             ' returned by dsyev for wavenumber ',jn
          call utl_abort('bcsc_corvertSetup')
        end if

        ! Set selected number of eigenmodes to zero
        if (numModeZero > 0) then
          do jk1 = 1, numModeZero
            eigenval(jk1) = 0.0d0
          end do
          write(*,*) 'bcsc_corvertSetup: modified eigenvalues=',eigenval(:)
        end if

        ! E * lambda^{-1}
        eigenvalmax=maxval(eigenval(1:jnum))
        do jk1 = 1, jnum
          if (eigenval(jk1) > 1.0d-8*eigenvalmax) then
            result(1:jnum,jk1) = eigenvec(1:jnum,jk1)/eigenval(jk1)
          else
            result(1:jnum,jk1) = 0.0D0
          end if
        end do

        ! E * lambda^{-1} * E^T
        bgStats%corverti(1:jnum,1:jnum,varIndex)=0.0D0
        do jk1 = 1, jnum
          do jk2 = 1, jnum
            bgStats%corverti(1:jnum,jk1,varIndex) = &
	      bgStats%corverti(1:jnum,jk1,varIndex) + result(1:jnum,jk2)* &
	      eigenvec(jk1,jk2)
          end do
        end do

        ! zr=maxval(abs(bgStats%corverti(1:jnum,1:jnum,varIndex)))
        ! where (abs(bgStats%corverti(1:jnum,1:jnum,varIndex)) < 1.E-5*zr) &
        !   bgStats%corverti(1:jnum,1:jnum,varIndex)=0.0D0

        ! Check inverse (for output when mmpi_myid is 0)
        result(1:jnum,1:jnum)=0.0D0
        do jk1 = 1, jnum
          do jk2 = 1, jnum
            result(1:jnum,jk1) = result(1:jnum,jk1) + &
              bgStats%corvert(1:jnum,jk2,varIndex)* &
	      bgStats%corverti(jk2,jk1,varIndex)
          end do
        end do

        cletiket = 'C*C^-1'
        clnomvar = bgStats%varNameList(varIndex)

        if (mmpi_myid == 0 .and. writePhysSpaceStats) then          
          ierr = utl_fstecr(result(1:jnum,1:jnum),-32,iulcorvert,0,0,0,jnum,jnum, &
                 1,0,0,bgStats%ntrunc,'X',clnomvar,cletiket,'X',0,0,0,0,5,.true.)
        end if
	 
        deallocate(eigenvec,result)
      end if
      
    end do
    
    if (mmpi_myid /= 0 .or. .not.WritePhysSpaceStats) return
    
    ierr = fstfrm(iulcorvert)  
    ierr = fclos(iulcorvert)

    jnum=nvlev
    call writecorns(0,'CORVERT',nvlev)
    if ( allocated(lfound) ) then
      if (any(.not.lfound(1:numvartot))) call writecorns(0,'CORVERTI',-1)
      deallocate(lfound)
    end if

  end subroutine bcsc_corvertSetup

  !--------------------------------------------------------------------------
  ! writecorns
  !--------------------------------------------------------------------------
  subroutine writecorns(nmat,cletiket,nlev)
  
    implicit none

    !Arguments
    character(len=*), intent(in) :: cletiket
    integer, intent(in) :: nmat,nlev

    !Locals
    integer :: jn, nulcorns,ierr,varIndex,jstart,jnum,numvartot

    ! standard file variables
    integer :: ip1,ip2,ip3
    integer :: idateo, ipak, idatyp
    character(len=4)  :: clnomvar
    integer :: fnom, fstouv, fstfrm, fclos
    
    if(mmpi_myid==0) then

      write(*,*) 'WRITECORNS: ', trim(cletiket), ' being written to file ', &
                 'bCovaarSetupChem_out.fst for number of matrices - 1 =', nmat

      nulcorns = 0
      ierr = fnom(nulcorns,'bCovarSetupChem_out.fst','STD+RND',0)
      ierr = fstouv(nulcorns,'RND')

      ipak = -32
      idatyp = 5
      ip1 = 0
      ip2 = 0
      ip3 = bgStats%ntrunc
      idateo = 0

      if (nmat == 0) then
        numvartot=bgStats%numvar3d
      else
        numvartot=bgStats%numvar3d+bgStats%numvar2d
      end if
      
      do varIndex = 1, numvartot
   
        if (any(CrossCornsVarKindCH(:) /= '') .and. nmat > 0) then
          jstart=1
          jnum=bgStats%nkgdim
          clnomvar = 'ZZ'
        else
          jstart = bgStats%nsposit(varIndex)
          jnum = bgStats%nsposit(varIndex+1)-bgStats%nsposit(varIndex)
          if (nlev > 0) jnum=nlev
          clnomvar = bgStats%varNameList(varIndex)
        end if
        
        if (nmat > 0 ) then
          do jn = 0, nmat
            ip2 = jn
            ierr = utl_fstecr(bgStats%corns(jstart:jstart+jnum-1, &
	           jstart:jstart+jnum-1,jn),ipak,nulcorns,idateo,0,0,jnum,  &
                   jnum,1,ip1,ip2,ip3,'X',clnomvar,cletiket,'X',0,0,0,0, &
		   idatyp,.true.)
          end do
        else
          if (trim(cletiket) == 'CORVERT') then
            ierr = utl_fstecr(bgStats%corvert(1:jnum,1:jnum,varIndex), &
	           ipak,nulcorns,idateo,0,0,jnum,jnum,1,  &
                   ip1,ip2,ip3,'X',clnomvar,cletiket,'X',0,0,0,0,idatyp,.true.)
          else
            ierr = utl_fstecr(bgStats%corverti(1:jnum,1:jnum,varIndex), &
	           ipak,nulcorns,idateo,0,0,jnum,jnum,1,  &
                   ip1,ip2,ip3,'X',clnomvar,cletiket,'X',0,0,0,0,idatyp,.true.)
          end if
        end if
        
        if (any(CrossCornsVarKindCH(:) /= '') .and. nmat > 0) exit
        
      end do
      
      ierr = fstfrm(nulcorns)  
      ierr = fclos(nulcorns)
    end if

  end subroutine writecorns

  !--------------------------------------------------------------------------
  ! readcorns
  !--------------------------------------------------------------------------
  subroutine readcorns(lfound,nmat,cletiket)

    implicit none

    !Arguments
    logical, intent(out) :: lfound(:)
    integer, intent(in) :: nmat
    character(len=*), intent(in) :: cletiket

    !Locals
    integer :: jn, icornskey,varIndex,jnum,jstart,numvartot
    real(8), allocatable :: zcornssrc(:,:)

    ! standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3
    integer :: idateo
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar

    if (nmat == 0) then
      numvartot=bgStats%numvar3d
      if (trim(cletiket) == 'CORVERT') then
        bgStats%corvert(:,:,:)=0.0d0
      else if (trim(cletiket) == 'CORVERTI') then
        bgStats%corverti(:,:,:)=0.0d0
      end if
    else if (nmat == bgStats%ntrunc) then
      numvartot=bgStats%numvar3d+bgStats%numvar2d
    end if
    
    write(*,*) 'readcorns: ', trim(cletiket), &
               ' being searched for number of matrices -1 =',nmat

    idateo = -1
    ip1 = -1
    ip3 = bgStats%ntrunc
    cltypvar = 'X'

    lfound(:)=.false.    
    VARCYCLE: do varIndex = 1, numvartot
   
      if (any(CrossCornsVarKindCH(:) /= '') .and. nmat > 0) then
        jstart=1
        jnum=bgStats%nkgdim
        clnomvar = 'ZZ'
      else
        jstart = bgStats%nsposit(varIndex)
        jnum = bgStats%nsposit(varIndex+1)-bgStats%nsposit(varIndex) 
        clnomvar = bgStats%varNameList(varIndex)
      end if
      allocate(zcornssrc(jnum,jnum))

      do jn = 0, nmat
        ip2 = jn

        ! Looking for FST record parameters..
  
        icornskey = utl_fstlir(ZCORNSSRC,nulbgst,INI,INJ,INK,idateo,cletiket, &
	                       ip1,ip2,ip3,cltypvar,clnomvar)

        if(icornskey  < 0 ) then
          write(*,*) 'readcorns: matrix not found in stats file for variable ', &
	              clnomvar
          deallocate(zcornssrc)
          cycle VARCYCLE
        end if

        if (ini /= jnum .or. inj /= jnum) then
	  call utl_abort('readcorns: BG stat levels inconsitencies')
	end if

        if (nmat > 0) then
          if (jn == 0) then
	    bgStats%corns(jstart:jstart+jnum-1,jstart:jstart+jnum-1,:)=0.0d0
	  end if
          bgStats%corns(jstart:jstart+jnum-1,jstart:jstart+jnum-1,jn)= &
	    zcornssrc(1:jnum,1:jnum)
        else if (trim(cletiket) == 'CORVERT') then
          bgStats%corvert(1:jnum,1:jnum,varIndex)=zcornssrc(1:jnum,1:jnum)
        else
          bgStats%corverti(1:jnum,1:jnum,varIndex)=zcornssrc(1:jnum,1:jnum)        
        end if
      end do
      lfound(varIndex)=.true.

      deallocate(zcornssrc)

      if (any(CrossCornsVarKindCH(:) /= '') .and. nmat > 0) exit

    end do VARCYCLE
    
  end subroutine readcorns

  !--------------------------------------------------------------------------
  ! gaspariCohn
  !--------------------------------------------------------------------------
  function gaspariCohn(ztlen,zr)

    !Arguments
    real(8) :: gasparicohn
    real(8), intent(in) :: ztlen,zr
    
    !Locals
    real(8)  :: zlc

    zlc = ztlen/2.0d0
    if(zr <= zlc) then
      gasparicohn = -0.250d0*(zr/zlc)**5+0.5d0*(zr/zlc)**4             &
                  +0.625d0*(zr/zlc)**3-(5.0d0/3.0d0)*(zr/zlc)**2+1.0d0
    elseif(zr <= (2.0d0*zlc)) then
      gasparicohn = (1.0d0/12.0d0)*(zr/zlc)**5-0.5d0*(zr/zlc)**4         &
                  +0.625d0*(zr/zlc)**3+(5.0d0/3.0d0)*(zr/zlc)**2       &
                  -5.0d0*(zr/zlc)+4.0d0-(2.0d0/3.0d0)*(zlc/zr)
    else
      gasparicohn = 0.0d0
    end if
    if(gasparicohn < 0.0d0) gasparicohn = 0.0d0

  end function gaspariCohn

  !--------------------------------------------------------------------------
  ! bcsc_finalize
  !--------------------------------------------------------------------------
  subroutine bcsc_finalize()
    implicit none

    if(bgStats%initialized) then
      if (allocated(bgStats%nsposit))  deallocate(bgStats%nsposit)
      if (allocated(bgStats%vlev))     deallocate(bgStats%vlev)
      if (allocated(bgStats%lat))      deallocate(bgStats%lat)
      if (allocated(bgStats%lon))      deallocate(bgStats%lon)       
      if (allocated(bgStats%stddev))   deallocate(bgStats%stddev)
      if (allocated(bgStats%corns))    deallocate(bgStats%corns)
      if (allocated(bgStats%hcorrlen)) deallocate(bgStats%hcorrlen)
      if (allocated(bgStats%corvert))  deallocate(bgStats%corvert)
      if (allocated(bgStats%corverti)) deallocate(bgStats%corverti)
      if (allocated(bgStats%invsum))   deallocate(bgStats%invsum)         
    end if

  end subroutine bcsc_finalize

  !--------------------------------------------------------------------------
  ! bcsc_getCovarCH
  !--------------------------------------------------------------------------
  subroutine bcsc_getCovarCH(bgStatsOut,transformVarKind_opt)
    !
    ! :Purpose: Pass on covariances in bgStats
    !
    implicit none

    !Arguments
    type(struct_bcsc_bgStats), intent(out) :: bgStatsOut  ! Structure with covariance elements

    ! Analysis variable transform prefix (*) of transforms *CH_tlm, *CH_ad,
    ! *CH_ens applicable via routine gvt_transform
    character(len=20), intent(out), optional :: TransformVarKind_opt  
			 
    if (present(TransformVarKind_opt)) TransformVarKind_opt=TransformVarKindCH

    if(.not.bgStats%initialized) then
      bgStatsOut%initialized=.false.
      call utl_abort('bcsc_getCovarCH: Covariances not set up.')
    end if
        
    bgStatsOut=bgStats
    
  end subroutine bcsc_getCovarCH

  !--------------------------------------------------------------------------
  ! bcsc_resetCorvert
  !--------------------------------------------------------------------------
  subroutine bcsc_resetCorvert(nlev_T,vlev_T)
    !
    ! :Purpose: Vertically interpolate error correlation matrix fields to 
    !           generate approximate matrices/vectors in trial field vertical
    !           levels via interpolation. No need to make matrix positive  
    !           definite for this approximation.
    !
    implicit none
    
    !Arguments
    integer, intent(in) :: nlev_T    ! Number of vertical levels for trial fields
    real(8), intent(in), pointer :: vlev_T(:) ! Trial field vertical levels
      
    integer :: ilev1,ilev2,j,d1,d2
    real(8) :: dz

    real(8), allocatable :: wtemp1(:,:), wtemp2(:,:,:)
        
    d2=0 
    d1=nlev_T-bgStats%nlev
    if (d1 < 0) then
      d1=0
      d2=d1
    end if
    
    allocate(wtemp2(nlev_T, nlev_T,bgStats%numvar3d))
    allocate(wtemp1(nlev_T,bgStats%numvar3d))
    wtemp2(:,:,:)=0.0d0
    wtemp1(:,:)=0.0d0
    
    ! Apply interpolation

    do ilev1 = 1, nlev_T
      if (bgStats%vlev(1) >= vlev_T(ilev1)) then
        wtemp1(ilev1,1:bgStats%numvar3d)= &
	      bgStats%invsum(1,1:bgStats%numvar3d)

        wtemp2(ilev1,1:bgStats%nlev+d2,1:bgStats%numvar3d)= &
	      bgStats%corvert(1,1:bgStats%nlev+d2,1:bgStats%numvar3d)
        wtemp2(1:bgStats%nlev+d2,ilev1,1:bgStats%numvar3d)= &
	      bgStats%corvert(1:bgStats%nlev+d2,1,1:bgStats%numvar3d)          
      else if (bgStats%vlev(bgStats%nlev) <= vlev_T(ilev1)) then
        wtemp1(ilev1,1:bgStats%numvar3d)= &
	      bgStats%invsum(bgStats%nlev,1:bgStats%numvar3d)

        wtemp2(ilev1,1+ilev1-bgStats%nlev-d2:ilev1,1:bgStats%numvar3d)= &
	      bgStats%corvert(bgStats%nlev,1-d2:bgStats%nlev, &
	                      1:bgStats%numvar3d)
        wtemp2(1+ilev1-bgStats%nlev-d2:ilev1,ilev1,1:bgStats%numvar3d)= &
	      bgStats%corvert(1-d2:bgStats%nlev,bgStats%nlev, &
	                      1:bgStats%numvar3d)          
      else 
        do ilev2=1,bgStats%nlev-1
          if (vlev_T(ilev1) >= bgStats%vlev(ilev2) .and. &
              vlev_T(ilev1) < bgStats%vlev(ilev2+1)) exit
        end do

        dz=log(vlev_T(ilev1)/bgStats%vlev(ilev2)) &
            /log(bgStats%vlev(ilev2+1)/bgStats%vlev(ilev2)) 
                  
        wtemp1(ilev1,1:bgStats%numvar3d)=bgStats%invsum(ilev2, &
	  1:bgStats%numvar3d)*(1.0-dz) &
          +bgStats%invsum(ilev2+1,1:bgStats%numvar3d)*dz

        do j=1-max(ilev1-d1,1),nlev_T-min(d1+ilev1,nlev_T)
          wtemp2(ilev1,ilev1+j,1:bgStats%numvar3d)= &
	        (bgStats%corvert(ilev2,ilev2+j,1:bgStats%numvar3d) &
                *(1.0-dz) + bgStats%corvert(ilev2+1,ilev2+1+j, &
		                            1:bgStats%numvar3d)*dz)
          wtemp2(ilev1+j,ilev1,1:bgStats%numvar3d)= &
	    wtemp2(ilev1,ilev1+j,1:bgStats%numvar3d)
        end do               
      end if
    end do
    
    bgStats%invsum(1:nlev_T,:) = wtemp1(:,:)
    bgStats%corvert(1:nlev_T,1:nlev_T,:) = wtemp2(:,:,:)
    
    deallocate(wtemp1,wtemp2)
    
  end subroutine bcsc_resetCorvert

  !--------------------------------------------------------------------------
  ! bcsc_resetHcorrlen
  !--------------------------------------------------------------------------
  subroutine bcsc_resetHcorrlen(nlev_T,vlev_T)
    !
    ! :Purpose: To interpolate horizontal correlation length 
    !           to trial field vertical levels.
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: nlev_T      ! Number of target vertical levels
    real(8), intent(in)  :: vlev_T(:)   ! Target vertical levels
    
    ! Locals:
    real(8) :: hcorrlen(nlev_T,bgStats%numvar3d) ! correlation lengths
    integer :: nlev,ilev1,ilev2
    real(8) :: dz

    nlev = bgStats%nlev
    
    ! Apply interpolation

    do ilev1=1, nlev_T
      if (bgStats%vlev(1) >= vlev_T(ilev1)) then
        hcorrlen(ilev1,:)=bgStats%hcorrlen(1,1:bgStats%numvar3d)
      else if (bgStats%vlev(nlev) <= vlev_T(ilev1)) then
        hcorrlen(ilev1,:)=bgStats%hcorrlen(nlev,1:bgStats%numvar3d)
      else  
        do ilev2=1,nlev-1
          if (vlev_T(ilev1) >= bgStats%vlev(ilev2) .and. &
              vlev_T(ilev1) <  bgStats%vlev(ilev2+1)) exit
        end do
        dz=log(vlev_T(ilev1)/ bgStats%vlev(ilev2)) &
          /log( bgStats%vlev(ilev2+1)/ bgStats%vlev(ilev2))
        hcorrlen(ilev1,:)=bgStats%hcorrlen(ilev2,1:bgStats%numvar3d)* &
	                  (1.0-dz)+bgStats%hcorrlen(ilev2+1, &
			  1:bgStats%numvar3d)*dz
      end if
    end do
    
    bgStats%hcorrlen(1:nlev_T,1:bgStats%numvar3d) = hcorrlen(:,:)
    
  end subroutine bcsc_resetHcorrlen

  !--------------------------------------------------------------------------
  ! bcsc_StatsExistForVarName
  !--------------------------------------------------------------------------
  logical function bcsc_StatsExistForVarName(varName)
    !
    ! :Purpose: To check whether covariances have been made available for the
    !           specified variable
    !
    implicit none

    !Arguments
    character(len=4), intent(in) :: varName

    if (allocated(bgStats%varNameList)) then
      if (any(bgStats%varNameList(:) == varName)) then
        bcsc_StatsExistForVarName = .true.
      else
        bcsc_StatsExistForVarName = .false.
      end if
    else
      bcsc_StatsExistForVarName = .false.
    end if
    
  end function bcsc_StatsExistForVarName

  !--------------------------------------------------------------------------
  ! bcsc_getBgStddev
  !--------------------------------------------------------------------------
  subroutine bcsc_getBgStddev(varName,maxsize,xlat,xlong,stddevOut,vlev_opt)
    !
    ! :Purpose: To interpolate error std. dev. to obs location.
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: varName ! Variable name
    integer, intent(in) :: maxsize  ! Max array size
    real(8), intent(in) :: xlat     ! Target latitude
    real(8), intent(in) :: xlong    ! Target longitude
    real(8), intent(out) :: stddevOut(:) ! Error std. dev.
    real(8), intent(in), optional :: vlev_opt(:) ! Target vertical levels
    
    ! Locals:
    integer :: varIndex,latIndex,lonIndex,levIndex,nlev,startPosition
    real(8) :: work(maxsize,2),zc1,zc2,rlat1,rlat2,rlong1,rlong2,zd1,zd2
    real(8) :: stddev_max
    integer :: ilev1,ilev2
    real(8) :: dz
      
    integer, parameter :: itype=0

    if (.not.bgStats%initialized) return
    
    ! Determine location and size of background error std. dev.

    do varIndex = 1, bgStats%numvar3d+bgStats%numvar2d
      if (trim(varName) == trim(bgStats%varNameList(varIndex))) then
        if (varIndex <= bgStats%numvar3d ) then
          nlev = bgStats%nlev
        else
          nlev = 1
        end if
        startPosition = bgStats%nsposit(varIndex)	
        exit
      end if
    end do
    if  (varIndex > bgStats%numvar3d+bgStats%numvar2d) &
      call utl_abort('bcsc_getbgStddev: Variable not found')
    
    if (.not.present(vlev_opt) .and. nlev /= maxsize ) then
      write(*,*) 'nlev, maxsize: ',nlev,maxsize
      call utl_abort('bcsc_getbgStddev: Inconsistent size')
    end if
 
    ! Determine reference longitude index 

    lonIndex = 2    
    do while (xlong > bgStats%lon(lonIndex) .and. lonIndex < bgStats%ni+1) 
      lonIndex = lonIndex+1
    end do

    ! Set interpolation weights

    rlong2 = bgStats%lon(lonIndex)
    rlong1 = bgStats%lon(lonIndex-1)
    
    zd2 = (xlong-rlong1)/(rlong2-rlong1)
    zd1 = 1.0-zd2
     
    ! Determine reference latitude index

    latIndex = 2 
    do while (xlat > bgStats%lat(latIndex) .and. latIndex < bgStats%nj) 
      latIndex = latIndex+1
    end do

    ! Set interpolation weights

    rlat2 = bgStats%lat(latIndex)
    rlat1 = bgStats%lat(latIndex-1)
    
    zc2 = (xlat-rlat1)/(rlat2-rlat1)
    zc1 = 1.0-zc2

    if (xlat <= rlat1) then
      zc1 = 1.0
      zc2 = 0.0
    else if (xlat >= rlat2) then
      zc1 = 0.0
      zc2 = 1.0
    end if
        
    ! Apply interpolation

    if (itype == 0) then
    
      ! Interpolation of variances and take square root

      work(1:nlev,1) = zc1*bgStats%stddev(lonIndex-1,latIndex-1, &
                                     startPosition:startPosition-1+nlev)**2 + &
                       zc2*bgStats%stddev(lonIndex-1,latIndex, &
		                     startPosition:startPosition-1+nlev)**2
    
      work(1:nlev,2) = zc1*bgStats%stddev(lonIndex,latIndex-1, &
                                     startPosition:startPosition-1+nlev)**2 + &
                      zc2*bgStats%stddev(lonIndex,latIndex, &
                                    startPosition:startPosition-1+nlev)**2 
       
      work(1:nlev,1) = zd1*work(1:nlev,1) + zd2*work(1:nlev,2)
      
      if (nlev /= maxsize) then
        do ilev1=1, maxsize
          if (bgStats%vlev(1) >= vlev_opt(ilev1)) then
            stddevOut(ilev1)=sqrt(work(1,1))
          else if (bgStats%vlev(nlev) <= vlev_opt(ilev1)) then
            stddevOut(ilev1)=sqrt(work(nlev,1))
          else  
            do ilev2=1,nlev-1
              if (vlev_opt(ilev1) >= bgStats%vlev(ilev2) .and. &
                  vlev_opt(ilev1) <  bgStats%vlev(ilev2+1)) exit
            end do
            dz=log(vlev_opt(ilev1)/ bgStats%vlev(ilev2)) &
               /log( bgStats%vlev(ilev2+1)/ bgStats%vlev(ilev2))
            stddevOut(ilev1)=sqrt(work(ilev2,1)*(1.0-dz)+work(ilev2+1,1)*dz)
          end if
        end do
      else
        stddevOut(1:nlev) = sqrt(work(1:nlev,1))
      end if 
    
    else 

      ! Interpolation of std. dev. (to reduce execution time)

      work(1:nlev,1) = zc1*bgStats%stddev(lonIndex-1,latIndex-1, &
                                     startPosition:startPosition-1+nlev)**2 + &      
                       zc2*bgStats%stddev(lonIndex-1,latIndex, &
                                     startPosition:startPosition-1+nlev)**2 

      work(1:nlev,2) = zc1*bgStats%stddev(lonIndex,latIndex-1, &
                                     startPosition:startPosition-1+nlev)**2 + &
                       zc2*bgStats%stddev(lonIndex,latIndex, &
                                     startPosition:startPosition-1+nlev)**2

      work(1:nlev,1) = zd1*work(1:nlev,1) + zd2*work(1:nlev,2)
    
      if (nlev /= maxsize) then
        do ilev1=1, maxsize
          if (bgStats%vlev(1) >= vlev_opt(ilev1)) then
            stddevOut(ilev1)=work(1,1)
          else if (bgStats%vlev(nlev) <= vlev_opt(ilev1)) then
            stddevOut(ilev1)=work(nlev,1)
          else  
            do ilev2=1,nlev-1
              if (vlev_opt(ilev1) >= bgStats%vlev(ilev2) .and. &
                  vlev_opt(ilev1) < bgStats%vlev(ilev2+1)) exit
            end do
            dz=log(vlev_opt(ilev1)/bgStats%vlev(ilev2)) &
               /log(bgStats%vlev(ilev2+1)/bgStats%vlev(ilev2))
            stddevOut(ilev1)=work(ilev2,1)*(1.0-dz)+work(ilev2+1,1)*dz
          end if
        end do
      else
        stddevOut(1:nlev) = work(1:nlev,1)
      end if 
      
    end if
    
    stddev_max = maxval(bgStats%stddev(lonIndex-1:lonIndex, &
                        latIndex-1:latIndex,  &                 
			startPosition:startPosition-1+nlev)) 
    do levIndex = 1, nlev
      if (stddevOut(levIndex) < 0.0 .or. stddevOut(levIndex) > 1.1*stddev_max) then
        write(*,*) 'bcsc_getbgStddev: Interpolated std dev incorrect:'
        write(*,*) 'bcsc_getbgStddev:   zc1,zc2,zd1,zd2 = ',zc1,zc2,zd1,zd2
        write(*,*) 'bcsc_getbgStddev:   stddevOut,stddev_max = ', &
	           stddevOut(levIndex),stddev_max
        write(*,*) 'bcsc_getbgStddev:   lonIndex,latIndex,j,xlong,xlat = ', &
	           lonIndex,latIndex,levIndex,xlong,xlat
        write(*,*) 'bcsc_getbgStddev:   rlong2,rlong1,rlat1,rlat2 = ', &
	           rlong2,rlong1,rlat1,rlat2
        call utl_abort('bcsc_getbgStddev')
      end if
    end do

  end subroutine bcsc_getBgStddev

  !--------------------------------------------------------------------------
  ! bcsc_addBgStddev
  !--------------------------------------------------------------------------
  subroutine bcsc_addBgStddev(headerIndex,stddevIn,obsdataMaxsize)
    !
    ! :Purpose: To add background stddev profiles (and inverse) to bgStddev
    !          which can be retrieved later using a header index.
    !
    implicit none 

    ! Arguments:
    integer, intent(in) :: headerIndex
    real(8), intent(in) :: stddevIn(:,:)
    integer, intent(in) :: obsdataMaxsize
     
    if (.not.associated(bgStddev%data2d)) then
      call oss_obsdata_alloc(bgStddev, obsdataMaxsize, dim1= &
           size(stddevIn,dim=1), dim2_opt=size(stddevIn,dim=2))
      bgStddev%nrep = 0
    end if

    ! In this case nrep will count the number of filled reps in the data arrays
    bgStddev%nrep = bgStddev%nrep+1 

    if (bgStddev%nrep > obsdataMaxsize) then
      call utl_abort('bcsc_addbgStddev: Reached max size of array ' // &
	             trim(utl_str(obsdataMaxsize)) )
    end if
    
    ! Use the header number as the unique code for this obs data
    write(bgStddev%code(bgStddev%nrep),'(I22)') headerIndex

    bgStddev%data2d(:,1,bgStddev%nrep) = stddevIn(:,1)

    where (stddevIn(:,1) > 0.0D0)
      bgStddev%data2d(:,2,bgStddev%nrep) = 1.0D0/stddevIn(:,1)
    elsewhere
      bgStddev%data2d(:,2,bgStddev%nrep) = 0.0D0
    end where

  end subroutine bcsc_addBgStddev

  !--------------------------------------------------------------------------
  ! bcsc_retrievebgStddev
  !--------------------------------------------------------------------------
  function bcsc_retrieveBgStddev(dim1,dim2,headerIndex) result(stddevOut)
    !
    ! :Purpose: To retrieve previously saved background stddev profiles
    !           in bgStddev from the header index of the observation.
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: dim1, dim2  ! Dimensions of output array
    integer, intent(in) :: headerIndex ! Index of observation

    real(8) :: stddevOut(dim1,dim2)

    ! Locals:
    character(len=22) :: code
    
    if (bgStddev%dim1 /= dim1 .or. bgStddev%dim2 /= dim2) then
      call utl_abort("bcsc_retrievebgStddev: Inconsitent dimensions")
    end if

    write(code,'(I22)') headerIndex
	  
    stddevOut = oss_obsdata_get_array2d(bgStddev,code)
    
  end function bcsc_retrieveBgStddev

end module BCovarSetupChem_mod
