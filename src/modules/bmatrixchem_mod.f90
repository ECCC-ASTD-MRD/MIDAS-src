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

module BmatrixChem_mod 
  ! MODULE BmatrixChem_mod (prefix='bchm' category='2. B and R matrices')
  !
  ! :Purpose: Contains routines involving the preparation and application of
  !           background-error covariance matrix(ces). Matrix based on
  !           horizontally homogeneous/isotropic correlations.  Specifically,
  !           the methods can:
  !
  !           1. read and prepare the static background error covariance matrix
  !              (if reading of balance operators if any are eventually added)
  !
  !           2. transform from control vector (spectral space) to analysis
  !              increments (and its adjoint) using the background error
  !              covariance matrix (and the balance operators when present).
  !
  !          Based on bmatrixHI_mod.ftn90
  !
  ! :Comments and questions:
  !
  !   1. Covariances uncoupled from weather variable.
  !
  !   2.  Handles univariate and multivariate covariances.
  !       See routines bchm_readcorns2 and bchm_sucorns2. 
  !
  !   3. One could potentially make public the functions/routines which are
  !      identical to those in bmatrixhi_mod.ftn90 (except possibly in name) so
  !      that one copy is present in the code.
  !
  !   4. For multiple univariate variables (or univarite blocks of one to
  !      multiple variables), one can alternatively have multiple sets of
  !      covariance matrices within this module instead of a single covariance
  !      matrix setup (similarly to what was done for bchm_corvert*).
  !
  !
  ! Public Subroutines (which call other internal routines/functions):
  !    bchm_setup:      Must be called first. Sets of background covariance
  !                     matrix (and balance operators if any are eventually
  !                     added)
  !
  !    bchm_BSqrt:      Transformations from control vector to analysis
  !                     increments in the minimization process.
  !
  !    bchm_BSqrtAd:    Adjoint of bchm_BSqrt.
  !    bchm_Finalize    Deallocate internal module arrays.
  !    bchm_corvert_mult: Multiple an input matrix/array with 'bchm_corvert' or
  !                     'bchm_corverti'
  !    bchm_getsigma:   Obtain background error std. dev. profile at
  !                     obs/specified location. 
  !    bchm_getBgSigma: Obtain background error std. dev. at specified grid
  !                     point for specified field.
  !    bchm_is_initialized: checks if B_chm has been intialized.
  !    bchm_StatsExistForVarname: Checfs if covariances available for specified
  !                     variable.

  use mpi_mod
  use mpivar_mod
  use MathPhysConstants_mod
  use earthConstants_mod
  use gridStateVector_mod
  use gridVariableTransforms_mod
  use globalSpectralTransform_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use varNameList_mod
  use obsSpaceData_mod
  use utilities_mod

  implicit none
  save
  private

  ! public procedures
  public :: bchm_Setup,bchm_BSqrt,bchm_BSqrtAd,bchm_Finalize
  public :: bchm_expandToMPIglobal,bchm_expandToMPIglobal_r4,bchm_reduceToMPIlocal,bchm_reduceToMPIlocal_r4,bchm_getScaleFactor
  public :: bchm_corvert_mult, bchm_getsigma, bchm_is_initialized
  public :: bchm_truncateCV,bchm_getBgSigma, bchm_StatsExistForVarname
  public :: bchm_resetCorvert

  ! public arrays
  public :: bchm_corvert, bchm_corverti, bchm_invsum, bchm_varnamelist 
  
  logical             :: initialized = .false.                    
  integer             :: nj_l,ni_l                    
  integer             :: nlev_M,nlev_T,nlev_T_even,nkgdim
  integer             :: nla_mpiglobal,nla_mpilocal           
  integer             :: cvDim_mpilocal,cvDim_mpiglobal           
  integer             :: gstID, gstID2          
  integer             :: nlev_bdl    
  real(8), allocatable   :: rlat(:),rlong(:),rlatr(:),rlongr(:)     
  type(struct_vco),pointer :: vco_anl

  real(8), pointer    :: pressureProfile_M(:),pressureProfile_T(:)

  integer             :: mymBeg,mymEnd,mymSkip,mymCount
  integer             :: mynBeg,mynEnd,mynSkip,mynCount
  integer             :: maxMyNla
  integer             :: myLatBeg,myLatEnd
  integer             :: myLonBeg,myLonEnd
  integer, pointer    :: ilaList_mpiglobal(:)
  integer, pointer    :: ilaList_mpilocal(:)
  
  character(len=15) :: bchm_mode
                            
  integer, external   :: get_max_rss
  integer             :: nulbgst=0

  ! Bacgkround error covariance matrix elements.
  ! One could add an additional dimension to corns  
  ! for separate block-univariate correlation matrices.
  ! This would also permit merging of bmatrixhi_mod and bmatrixchem_mod
  ! into one module.

  real(8),allocatable :: rgsig(:,:,:)
  real(8),allocatable :: corns(:,:,:)
  real(8),allocatable :: rstddev(:,:)
  real(8),allocatable :: bchm_hcorrel(:,:,:),hdist(:),bchm_hcorrlen(:,:)

  ! Physical space (total) vertical correlation matrices and its inverse.
  
  real(8), allocatable, dimension(:,:,:) :: bchm_corvert,bchm_corverti
  real(8), allocatable :: bchm_invsum(:,:)

  ! Parameters of the NAMBCHM namelist

  integer             :: ntrunc
  real(8)             :: rpor(vnl_numvarmax)
  real(8)             :: rvloc(vnl_numvarmax)
  integer,parameter   :: maxNumLevels=200
  real(8)             :: scaleFactor(maxNumLevels,vnl_numvarmax)
  real(8)             :: scaleFactor_sigma(maxNumLevels,vnl_numvarmax)
  integer             :: numModeZero  ! number of eigenmodes to set to zero
  logical             :: ReadWrite_sqrt,ReadWrite_PSStats
  logical             :: getPhysSpace_hcorrel
  character(len=4)    :: stddevMode
  character(len=4)    :: IncludeAnlVarKindCH(vnl_numvarmax)
  character(len=4)    :: CrossCornsVarKindCH(vnl_numvarmax)
  character(len=20)   :: TransformVarKindCH
 
  ! Number of incremental variables/fields
  integer             :: numvar3d,numvar2d
  ! Start position of each field in composite arrays
  integer, allocatable :: nsposit(:)
  ! Name list of incremental variables/fields
  character(len=4),allocatable    :: bchm_varNameList(:)
  
  ! Indicates if vertical levels of bchm_corvert and bchm_invsum
  ! have been reset for consistency with trial field vertical coordinate
  logical :: lresetCorvert= .false.
  
  ! Reference surface pressure
  real(8), parameter :: zps = 101000.D0

  !*************************************************************************
    
  contains

  !--------------------------------------------------------------------------
  ! bchm_setup
  !--------------------------------------------------------------------------
  subroutine bchm_setup(hco_in,vco_in,CVDIM_OUT,mode_opt)
    !
    !:Purpose: To set up for constituents static background error covariances.
    !
    implicit none

    !Arguments
    type(struct_hco),pointer :: hco_in
    type(struct_vco),pointer :: vco_in
    integer                  :: cvDim_out
    character(len=*), intent(in), optional :: mode_opt

    !Locals
    integer :: nulnam, ierr, fnom, fclos, jm, jn, status
    integer :: latPerPE, latPerPEmax, lonPerPE, lonPerPEmax

    integer :: VarIndex,nChmVars,VarIndex2
    character(len=4) :: BchmVars(vnl_numvarmax)
    
    NAMELIST /NAMBCHM/ntrunc,rpor,rvloc,scaleFactor,numModeZero,ReadWrite_sqrt,stddevMode, &
                      IncludeAnlVarKindCH,getPhysSpace_hcorrel,CrossCornsVarKindCH,ReadWrite_PSStats, &
                      TransformVarKindCH
   
   ! First check if there are any CH fields 
    
    VarIndex2=0
    do VarIndex = 1, vnl_numvarmax
      if (gsv_varExist(varName = vnl_varNameList(VarIndex))) then
        if (vnl_varKindFromVarname(vnl_varNameList(VarIndex)) == 'CH') then
          VarIndex2 = 1
          exit
        end if 
      end if      
    end do
    if (VarIndex2 == 0) then
       ! Assume there is no need for Bchm
       cvDim_out = 0
       return
    end if

    numvar3d = 0
    numvar2d = 0
 
    allocate(bchm_varNameList(vnl_numvarmax))
    bchm_varNameList(:) = ''
    allocate(nsposit(vnl_numvarmax+1))
    nsposit(1) = 1

    if ( present(mode_opt) ) then
       if ( trim(mode_opt) == 'Analysis' .or. trim(mode_opt) == 'BackgroundCheck') then
         bchm_mode = trim(mode_opt)
         if(mpi_myid == 0) write(*,*)
         if(mpi_myid == 0) write(*,*) 'bchm_setup: Mode activated = ', trim(bchm_mode)
       else
          write(*,*)
          write(*,*) 'mode = ', trim(mode_opt)
          call utl_abort('bchm_setup: unknown mode')
       end if
    else
       bchm_mode = 'Analysis'
       if(mpi_myid == 0) write(*,*)
       if(mpi_myid == 0) write(*,*) 'bchm_setup: Analysis mode activated (by default)'
    end if

    ! Initialization of namelist NAMBCHM parameters
    
    ntrunc=108
    rpor(:)=3000.D3
    rvloc(:)=4.0D0
    scaleFactor(:,:) = 1.0d0
    numModeZero = 0
    ReadWrite_sqrt = .false.
    ReadWrite_PSStats = .false.
    getPhysSpace_hcorrel = .false.
    stddevMode = 'GD3D'    
    IncludeAnlVarKindCH(:) = ''
    CrossCornsVarKindCH(:) = ''
    TransformVarKindCH = ''
            
    ! Read namelist input
    
    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=NAMBCHM,iostat=ierr)
    if(ierr.ne.0) call utl_abort('bchm_setup: Error reading namelist')
    if(mpi_myid == 0) write(*,nml=NAMBCHM)
    ierr = fclos(nulnam)

    ! Set BchmVars
    nChmVars=0
    BChmVars(:)=''
    if (trim(IncludeAnlVarKindCH(1)) == '') then
      do VarIndex = 1, vnl_numvarmax
        if (.not. gsv_varExist(varName = vnl_varNameList(VarIndex))) cycle
        if (vnl_varKindFromVarname(vnl_varNameList(VarIndex)) /= 'CH') cycle
        nChmVars = nChmVars+1
        BchmVars(nChmVars) = trim(vnl_varNameList(VarIndex))
      end do
    else
      do VarIndex = 1, vnl_numvarmax
        if (.not. gsv_varExist(varName = vnl_varNameList(VarIndex))) cycle
        if (vnl_varKindFromVarname(vnl_varNameList(VarIndex)) /= 'CH') cycle
        do VarIndex2 = 1, vnl_numvarmax
          if (trim(IncludeAnlVarKindCH(VarIndex2)) == '') exit
          if (trim(vnl_varNameList(VarIndex)) == trim(IncludeAnlVarKindCH(VarIndex2))) then
            nChmVars = nChmVars+1
            BchmVars(nChmVars)= trim(vnl_varNameList(VarIndex))
            exit
          end if
        end do
      end do
    end if

    if (nChmVars.eq.0) then 
       if(mpi_myid == 0) then
           write(*,*) 'Size of BchmVars is zero. Bhi matrix for CH family not produced.'
           write(*,*) 'No chemical assimilation to be performed.'
           write(*,*) 'END OF BCHM_SETUP'
       end if
       cvdim_out = 0
       return
    end if
    
    ! Set vertical dimensions

    vco_anl => vco_in
    nLev_M = vco_anl%nlev_M
    nLev_T = vco_anl%nlev_T
    if(mod(nLev_T,2).ne.0) then
      nLev_T_even = nLev_T+1
    else
      nLev_T_even = nLev_T
    endif
    if(mpi_myid == 0) write(*,*) 'bchm_setup: nLev_T, nLev_T_even=',nLev_T, nLev_T_even

    !   Find the 3D variables (within NAMSTATE namelist)

    do VarIndex = 1, vnl_numvarmax3D    
       if (gsv_varExist(varName=vnl_varNameList3D(VarIndex)) .and. &
           any(trim(vnl_varNameList3D(VarIndex))==BchmVars(1:nChmVars)) ) then

           if (vnl_varKindFromVarname(vnl_varNameList3D(VarIndex)) /= 'CH') cycle
           numvar3d = numvar3d + 1
           nsposit(numvar3d+1)=nsposit(numvar3d)+nLev_T
           bchm_varNameList(numvar3d)=vnl_varNameList3D(VarIndex)
       end if
    end do
 
    !   Find the 2D variables (within NAMSTATE namelist)

    do VarIndex = 1, vnl_numvarmax2D
      if (gsv_varExist(varName=vnl_varNameList2D(VarIndex)) .and. &
          any(trim(vnl_varNameList2D(VarIndex)) == BchmVars(1:nChmVars)) ) then

        if (vnl_varKindFromVarname(vnl_varNameList2D(VarIndex)) /= 'CH') cycle
        numvar2d = numvar2d + 1
        nsposit(numvar3d+numvar2d+1) = nsposit(numvar3d+numvar2d)+1
        bchm_varNameList(numvar2d) = vnl_varNameList2D(VarIndex)
      end if       
    end do
    
    if (numvar3d+numvar2d == 0) then    
      if (mpi_myid == 0) then
        write(*,*) 'Bhi matrix for CH family not produced.'
        write(*,*) 'No chemical assimilation to be performed.'
        write(*,*) 'END OF BCHM_SETUP'
      end if
      cvDim_out = 0
      return
    else if (mpi_myid == 0) then
      if (numvar3d > 0) &
        write(*,*) 'bchm_setup: Number of 3D variables', numvar3d,bchm_varNameList(1:numvar3d)
      if (numvar2d > 0) &
        write(*,*) 'bchm_setup: Number of 2D variables', numvar2d,bchm_varNameList(numvar3d+1:numvar3d+numvar2d)
    end if
    
    nkgdim =  nsposit(numvar3d+numvar2d+1)-1

    ! Assumes the input 'scalefactor' is a scaling factor of the variances.
    
    where (scaleFactor(1:max(nLev_M,nLev_T),1:numvar3d+numvar2d) < 0.0) &
           scaleFactor(1:max(nLev_M,nLev_T),1:numvar3d+numvar2d) = 1.0d0
    scaleFactor_sigma(1:max(nLev_M,nLev_T),1:numvar3d+numvar2d) =   &
       sqrt(scaleFactor(1:max(nLev_M,nLev_T),1:numvar3d+numvar2d))

    nla_mpiglobal = (ntrunc+1)*(ntrunc+2)/2
    ni_l = hco_in%ni
    nj_l = hco_in%nj
    if (allocated(rlat)) deallocate(rlat)
    if (allocated(rlong)) deallocate(rlong)
    allocate(rlat(nj_l))
    allocate(rlong(ni_l+1))
    rlat(1:nj_l) = hco_in%lat(1:nj_l)
    rlong(1:ni_l) = hco_in%lon(1:ni_l)
    rlong(ni_l+1) = 2*MPC_PI_R8
    ! rlong(1) = 0.0
    ! do jn = 1, ni_l+1
    !   rlong(jn) = 2*MPC_PI_R8/(ni_l)*(jn-1)
    ! end do
    
    ! Begin calcs.

    ! Need an even number of levels for spectral transform (gstID2)
    
    gstID  = gst_setup(ni_l,nj_l,ntrunc,nkgdim)
    gstID2 = gst_setup(ni_l,nj_l,ntrunc,nlev_T_even)
    if(mpi_myid == 0) write(*,*) 'BCHM:returned value of gstID =',gstID
    if(mpi_myid == 0) write(*,*) 'BCHM:returned value of gstID2=',gstID2

    call mpivar_setup_latbands(nj_l, latPerPE, latPerPEmax, myLatBeg, myLatEnd)
    call mpivar_setup_lonbands(ni_l, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd)

    call mpivar_setup_m(ntrunc,mymBeg,mymEnd,mymSkip,mymCount)
    call mpivar_setup_n(ntrunc,mynBeg,mynEnd,mynSkip,mynCount)

    call gst_ilaList_mpiglobal(ilaList_mpiglobal,nla_mpilocal,maxMyNla,gstID,mymBeg,mymEnd,mymSkip,mynBeg,mynEnd,mynSkip)
    call gst_ilaList_mpilocal(ilaList_mpilocal,gstID,mymBeg,mymEnd,mymSkip,mynBeg,mynEnd,mynSkip)

    ! compute mpilocal control vector size
    do jm = mymBeg, mymEnd, mymSkip
      do jn = mynBeg, mynEnd, mynSkip
        if(jm.le.jn) then
          if(jm == 0) then
            ! only real component for jm=0
            cvDim_mpilocal = cvDim_mpilocal + 1*nkgdim
          else
            ! both real and imaginary components for jm>0
            cvDim_mpilocal = cvDim_mpilocal + 2*nkgdim
          endif
        endif
      enddo
    enddo
    cvDim_out = cvDim_mpilocal

    ! also compute mpiglobal control vector dimension
    call rpn_comm_allreduce(cvDim_mpilocal, cvDim_mpiglobal, 1, &
                            "mpi_integer", "mpi_sum", "GRID", ierr)

    allocate(rgsig(ni_l+1, nj_l, nkgdim))
    allocate(corns(nkgdim, nkgdim, 0:ntrunc))  
    allocate(rstddev(nkgdim, 0:ntrunc))  
    allocate(bchm_corvert(nlev_T, nlev_T, numvar3d+numvar2d))
    allocate(bchm_corverti(nlev_T, nlev_T, numvar3d+numvar2d))
    allocate(bchm_invsum(nlev_T, numvar3d+numvar2d))

    status = vgd_levels( vco_anl%vgrid, ip1_list=vco_anl%ip1_M, levels=pressureProfile_M, &
                         sfc_field=zps, in_log=.false.)
    status = vgd_levels( vco_anl%vgrid, ip1_list=vco_anl%ip1_T, levels=pressureProfile_T, &
                         sfc_field=zps, in_log=.false.)
    
    ! Read and apply localization to vertical correlation matrices in horizontal 
    ! spectral space and read and scale standard deviations.
    
    call bchm_rdstats

    ! Generate or read correlation matrix square roots
    
    call bchm_sucorns2
    
    if(mpi_myid == 0) write(*,*) 'END OF BCHM_SETUP'
    
    initialized = .true.

  end subroutine bchm_setup

  !--------------------------------------------------------------------------
  ! bchm_StatsExistForVarName
  !--------------------------------------------------------------------------
  logical function bchm_StatsExistForVarName(VarName)
    !
    !:Purpose: To check whether covariances have been made available for the
    !          specified variable
    !
    implicit none

    !Arguments
    character(len=4), intent(in) :: VarName

    if (allocated(bchm_varNameList)) then
       if (any(bchm_varNameList(:) == VarName)) then
          bchm_StatsExistForVarName = .true.
       else
          bchm_StatsExistForVarName = .false.
       end if
    else
       bchm_StatsExistForVarName = .false.
    end if
    
  end function bchm_StatsExistForVarName

  !--------------------------------------------------------------------------
  ! bchm_is_initialized
  !--------------------------------------------------------------------------
  logical function bchm_is_initialized()
    !
    !:Purpose: To check whether B_chm has been initialized
    !
    implicit none

    bchm_is_initialized = initialized

  end function bchm_is_initialized

  !--------------------------------------------------------------------------
  ! bchm_getScaleFactor
  !--------------------------------------------------------------------------
  subroutine bchm_getScaleFactor(scaleFactor_out)
    !
    !:Purpose: To set scaling factors for background error std. dev.
    !
    implicit none

    !Arguments
    real(8) :: scaleFactor_out(:,:)

    !Locals
    integer :: levelIndex,VarIndex

    do VarIndex = 1, numvar3d+numvar2d
    do levelIndex = 1, nsposit(VarIndex+1)-nsposit(VarIndex)
      scaleFactor_out(levelIndex,VarIndex) = scaleFactor(levelIndex,VarIndex)
    end do
    end do

  end subroutine bchm_getScaleFactor

  !--------------------------------------------------------------------------
  ! bchm_rdstats
  !--------------------------------------------------------------------------
  subroutine bchm_rdstats
    !
    !:Purpose: To read chemical constituents background stats file.
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
        call utl_abort('bchm_RDSTATS: Problem in opening the background chemical constituent stat file')
      endif
    else
      call utl_abort('bchm_RDSTATS: Background chemical constituent stat file is missing')
    endif

    ! check if analysisgrid and covariance file have the same vertical levels
    call vco_SetupFromFile( vco_file,  & ! OUT
                            bFileName )  ! IN
    if (.not. vco_equal(vco_anl,vco_file)) then
      call utl_abort('bmatrixchem: vco from analysisgrid and chem cov file do not match')
    end if

    ! Read spectral space correlations
    
    call bchm_readcorns2
    
    ! Read error standard deviations
    
    call bchm_rdstddev 

    ! Scale error standard deviations
    
    call bchm_scalestd
    
    ierr = fstfrm(nulbgst)
    ierr = fclos(nulbgst)

  end subroutine bchm_rdstats

  !--------------------------------------------------------------------------
  ! bchm_scalestd
  !--------------------------------------------------------------------------
  subroutine bchm_scalestd
    !
    !:Purpose: To scale error standard-deviation values.
    !
    implicit none

    !Locals
    integer :: lonIndex, latIndex, VarIndex, levelIndex, nlev, nulsig
    integer :: ierr, fnom, fstouv, fstfrm, fclos
  
    if (ReadWrite_PSStats .and. mpi_myid == 0) then
       nulsig = 0
       ierr = fnom(nulsig,'bmatrixchem_out.fst','STD+RND',0)
       ierr = fstouv(nulsig,'RND')
       ierr = utl_fstecr(pressureProfile_T,-32,nulsig,0,0,0,1,1,nlev_T,0,0,0, &
                   'X','PX','Pressure','X',0,0,0,0,5,.true.)
    end if
    
    do VarIndex = 1,numvar3d+numvar2d
      nlev=nsposit(VarIndex+1)-nsposit(VarIndex)
      do lonIndex = 1, ni_l+1
      do latIndex = 1, nj_l
!        rgsig(latIndex,nsposit(VarIndex):nsposit(VarIndex+1)-1) = &
!               scaleFactor_sigma(1:nlev,VarIndex)* &
!               rgsig(latIndex,nsposit(VarIndex):nsposit(VarIndex+1)-1)
        rgsig(lonIndex,latIndex,nsposit(VarIndex):nsposit(VarIndex+1)-1) = &
               scaleFactor_sigma(1:nlev,VarIndex)* &
               rgsig(lonIndex,latIndex,nsposit(VarIndex):nsposit(VarIndex+1)-1)
      enddo
      enddo

      if (ReadWrite_PSStats .and. mpi_myid == 0) then
         do levelIndex=1,nlev
            ierr = utl_fstecr(rgsig(1:ni_l+1,1:nj_l,nsposit(VarIndex)-1+levelIndex),-32,nulsig,0,0,0,ni_l+1,nj_l,1,levelIndex,0,nlev, &
                   'X',bchm_varNameList(VarIndex),'RGSIG','X',0,0,0,0,5,.true.)
         end do
      end if

    enddo
    
    if (ReadWrite_PSStats .and. mpi_myid == 0) then
      ierr = fstfrm(nulsig)  
      ierr = fclos(nulsig)
    end if

  end subroutine bchm_scalestd
  !--------------------------------------------------------------------------
  ! bchm_readCorns2
  !--------------------------------------------------------------------------
  subroutine bchm_readCorns2
    !
    !:Purpose: To read correlation information and to form the correlation
    !          matrix.
    !
    !:Notes: Can read distinct block diagonal matrices with or without
    !        cross-correlations.
    !

    ! Based on bhi_readcorns2.
    implicit none

    !Locals
    integer :: jn,ierr,VarIndex,VarIndex2
    integer :: jcol,jrow,jstart,jnum,jstart2,jnum2
    real(8), allocatable, dimension(:) :: zstdsrc
    real(8), allocatable, dimension(:,:) :: zcornssrc

    ! Standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3,idateo
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=4), allocatable :: clnomvar_crosscorns(:)
    character(len=12) :: cletiket
    integer :: fstinf

    rstddev(:,:) = 0.0d0
    corns(:,:,:) = 0.0d0
    if (any(CrossCornsVarKindCH(:) /= '')) then
      allocate(clnomvar_crosscorns(numvar3d+numvar2d))
      clnomvar_crosscorns(:)=''
    end if
    
    ! Read auto-correlations matrices
    
    do VarIndex = 1, numvar3d+numvar2d
   
      clnomvar = bchm_varNameList(VarIndex)
      jstart = nsposit(VarIndex)
      jnum = nsposit(VarIndex+1)-nsposit(VarIndex)
      allocate(zcornssrc(jnum,jnum))
      allocate(zstdsrc(jnum))

      do jn = 0, ntrunc

        ! Looking for FST record parameters..
      
        idateo = -1
        cletiket = 'RSTDDEV'
        ip1 = -1
        ip2 = jn
        ip3 = -1
        cltypvar = 'X'
          
        ierr = utl_fstlir(ZSTDSRC,nulbgst,INI,INJ,INK,idateo,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
          
        if (ierr < 0 .and. ip2 < 10 .and. all(CrossCornsVarKindCH(:) == '')) then
          write(*,*) 'bchm_READCORNS2: RSTDDEV ',ip2,jnum,clnomvar
          call utl_abort('bchm_READCORNS2: Problem with constituents background stat file')
        else if (ierr < 0 .and. ip2 == 0 .and. any(CrossCornsVarKindCH(:) /= '')) then
          write(*,*) 'bchm_READCORNS2: Assumes content from cross-corrns input for ',clnomvar
          clnomvar_crosscorns(VarIndex)=clnomvar
          exit
        end if
        if (ini /= jnum)  call utl_abort('bchm_READCORNS2: Constituents background stat levels inconsistencies')

        ! Looking for FST record parameters..

        if (ierr >= 0) then
          idateo = -1
          cletiket = 'CORRNS'
          ip1 = -1
          ip2 = jn
          ip3 = -1
          cltypvar = 'X'
          ierr = utl_fstlir(ZCORNSSRC,nulbgst,INI,INJ,INK,idateo,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)

          if (ierr < 0) then
            write(*,*) 'bchm_READCORNS2: CORRNS ',ip2,jnum,clnomvar
            call utl_abort('bchm_READCORNS2: Problem with constituents background stat file')
          end if
          if (ini /= jnum .and. inj /= jnum) call utl_abort('bchm_READCORNS2: Constituents BG stat levels inconsistencies')
        else
          write(*,*) 'WARNING from bchm_READCORNS2: Setting RSDTDDEV to 1.D-15 for NOMVAR and JN: ',clnomvar,' ',jn
          zstdsrc(1:jnum) = 1.D-15
          zcornssrc(1:jnum,1:jnum) = 0.0D0
          do jrow = 1, jnum
            zcornssrc(jrow,jrow) = 1.0D0
          end do
        end if
          
        rstddev(jstart:jstart+jnum-1,jn) = zstdsrc(1:jnum)
        corns(jstart:jstart+jnum-1,jstart:jstart+jnum-1,jn)=zcornssrc(1:jnum,1:jnum)
          
      enddo
       
      deallocate(zcornssrc)
      deallocate(zstdsrc)

    enddo

    ! Read cross-correlation matrices
    
    if (any(CrossCornsVarKindCH(:) /= '')) then  
     
      do VarIndex = 1, numvar3d+numvar2d
   
        if ( all(CrossCornsVarKindCH(:) /= bchm_varNameList(VarIndex)) ) cycle

        clnomvar = bchm_varNameList(VarIndex)
        jstart = nsposit(VarIndex)
        jnum = nsposit(VarIndex+1)-nsposit(VarIndex)

        if (clnomvar_crosscorns(VarIndex) == clnomvar) clnomvar_crosscorns(VarIndex)=''
        
        do VarIndex2 = 1, numvar3d+numvar2d
        
          if (VarIndex == VarIndex2) cycle
          
          cletiket='CORRNS '//bchm_varNameList(VarIndex2)
          ierr = fstinf(nulbgst,INI,INJ,INK,-1,cletiket,-1,-1,-1,'X',clnomvar)
          if (ierr < 0 ) cycle
          
          jstart2 = nsposit(VarIndex2)
          jnum2 =  nsposit(VarIndex2+1)-nsposit(VarIndex2)
          allocate(zcornssrc(jnum,jnum2))
          
          if (clnomvar_crosscorns(VarIndex2) == bchm_varNameList(VarIndex2)) clnomvar_crosscorns(VarIndex2)=''
          
          do jn = 0, ntrunc
 
             ierr = utl_fstlir(ZCORNSSRC,nulbgst,INI,INJ,INK,-1,cletiket,-1,jn,-1,'X',clnomvar)

             if (ierr < 0) then
               if (jn < 10) then
                 write(*,*) 'bchm_READCORNS2: CORRNS ',ip2,jnum,clnomvar,bchm_varNameList(VarIndex2)
                 call utl_abort('bchm_READCORNS2: Problem with constituents background stat file')
               else
                 exit
               end if
             end if
             if (ini /= jnum .and. inj /= jnum2) call utl_abort('bchm_READCORNS2: Constituents BG2 stat levels inconsistencies')
          
             corns(jstart:jstart+jnum-1,jstart2:jstart2+jnum2-1,jn)=zcornssrc(1:jnum,1:jnum2)
             corns(jstart2:jstart2+jnum2-1,jstart:jstart+jnum-1,jn)=transpose(zcornssrc(1:jnum,1:jnum2))
          
          end do
          deallocate(zcornssrc)
        end do
      end do   
    end if
    
    if (any(CrossCornsVarKindCH(:) /= '')) then
      if (any(clnomvar_crosscorns(:) /= '')) then
         write(*,*) 'bchm_READCORNS2: Missing matrix for ',clnomvar_crosscorns(1:numvar3d+numvar2d)
         call utl_abort('bchm_READCORNS2: Missing correlations matrix')      
      end if
      deallocate(clnomvar_crosscorns)
    end if
     
    ! Apply convolution to RSTDDEV correlation

    call bchm_convol
    
    do jn = 0, ntrunc

      ! Re-build correlation matrix: factorization of corns with convoluted RSTDDEV
      do jcol = 1, nkgdim
        corns(1:nkgdim,jcol,jn) = rstddev(1:nkgdim,jn) * corns(1:nkgdim,jcol,jn) * rstddev(jcol,jn)
      enddo

    enddo

    !write(*,*) 'Done in bchm_READCORNS2'
  end subroutine bchm_readCorns2

  !--------------------------------------------------------------------------
  ! bchm_convol
  !--------------------------------------------------------------------------
  subroutine bchm_convol
    implicit none

    !Locals
    real(8) :: dlfact2,dlc,dsummed
    real(8) :: dtlen,zr,dlfact
    integer :: jn,latIndex,jk,VarIndex,levelIndex,nsize,ierr
    real(8) :: zleg(0:ntrunc,nj_l),zsp(0:ntrunc,nkgdim),zgr(nj_l,nkgdim)
    real(8) :: dlwti(nj_l),zrmu(nj_l)

    integer :: nlev_MT,ini,inj,ink,nulcorns
    real(8), allocatable :: wtemp(:,:,:)    
    logical :: lfound
    integer :: fnom, fstouv, fstfrm, fclos

    do latIndex = 1, nj_l
      dlwti(latIndex) = gst_getrwt(latIndex,gstID)
      zrmu(latIndex)  = gst_getrmu(latIndex,gstID)
    end do

    ! CONVERT THE CORRELATIONS IN SPECTRAL SPACE INTO SPECTRAL
    ! COEFFICIENTS OF THE CORRELATION FUNCTION AND FUNCTION TO BE
    ! SELF-CONVOLVED

    do jn = 0, ntrunc
      dlfact = ((2.0d0*jn+1.0d0)/2.0d0)**(0.25d0)
      dlfact2 = ((2.0d0*jn +1.0d0)/2.0d0)**(0.25d0)
      do jk = 1, nkgdim
        zsp(jn,jk) = rstddev(jk,jn)*dlfact*dlfact2
      enddo
    enddo

    ! Transform to physical space
    call gst_zleginv(gstID,zgr,zsp,nkgdim)
    
    ! Truncate in horizontal extent with Gaussian window
    
    VarIndex=1
    do jk = 1, nkgdim
      if (jk == nsposit(VarIndex)) then
        dtlen = rpor(VarIndex)
        VarIndex=VarIndex+1 
      endif
      if(dtlen.gt.0.0d0) then
        dlc = 1.d0/dble(dtlen)
        dlc = 0.5d0*dlc*dlc
        do latIndex = 1, nj_l
          zr = ec_ra * acos(zrmu(latIndex))
          dlfact = dexp(-(zr**2)*dlc)
          zgr(latIndex,jk) = dlfact*zgr(latIndex,jk)
        enddo
      endif

      !write(*,*) 'zeroing length (km)=',jk,dtlen/1000.0
    enddo

    ! Transform back to spectral space
    call gst_zlegdir(gstID,zgr,zsp,nkgdim)

    ! Convert back to correlations
    do jk = 1, nkgdim
      do jn = 0, ntrunc
         zsp(jn,jk) = zsp(jn,jk)*(2.0d0/(2.0d0*jn+1.0d0))**(0.25d0)
      enddo
    enddo

    ! PUT BACK INTO RSTDDEV
    do jn = 0, ntrunc
      do jk = 1, nkgdim
         rstddev(jk,jn) = zsp(jn,jk)
      enddo
    enddo

    ! Re-normalize to ensure correlations
    do jk = 1, nkgdim
      dsummed = 0.d0
      do jn = 0, ntrunc
        dsummed = dsummed+ dble(rstddev(jk,jn)**2)*sqrt(((2.d0*jn)+1.d0)/2.d0)
      enddo
      dsummed = sqrt(dsummed)
      do jn = 0, ntrunc
        if(dsummed.gt.1.d-30) rstddev(jk,jn) = rstddev(jk,jn)/dsummed
      enddo
    enddo

    ! CONVERT THE SPECTRAL COEFFICIENTS OF THE CORRELATION FUNCTION
    ! BACK INTO CORRELATIONS OF SPECTRAL COMPONENTS
    do jn = 0, ntrunc
      dlfact = sqrt(0.5d0)*(1.0d0/((2.0d0*jn+1)/2.0d0))**0.25d0
      do jk = 1, nkgdim
        rstddev(jk,jn) = rstddev(jk,jn)*dlfact
      enddo
    enddo

    if ( .not.getPhysSpace_hcorrel .or. .not.ReadWrite_PSStats ) return

    ! Compute resultant physical space horizontal correlations and
    ! 1/e correlation length from correlation array if not available
    
    if (allocated(bchm_hcorrel)) deallocate(bchm_hcorrel)
    allocate(bchm_hcorrel(nj_l, nlev_T, numvar3d+numvar2d))
    if (allocated(wtemp)) deallocate(wtemp)
    allocate(wtemp(0:ntrunc,nj_l,1))
    if (allocated(bchm_hcorrlen)) deallocate(bchm_hcorrlen)
    allocate(bchm_hcorrlen(nlev_T, numvar3d+numvar2d))
    if (allocated(hdist)) deallocate(hdist)
    allocate(hdist(nj_l))

    lfound=.true.
    do VarIndex = 1, numvar3d+numvar2d
      nlev_MT = nsposit(VarIndex+1)-nsposit(VarIndex)
      do levelIndex = 1, nlev_MT
         ierr = utl_fstlir(bchm_hcorrel(:,levelIndex,VarIndex),nulbgst,INI,INJ,INK,-1,'HCORREL',levelIndex,-1,-1,'X',bchm_varNameList(VarIndex))
         if (ierr < 0 ) then
           lfound=.false.
           exit
         end if
      end do
      ierr = utl_fstlir(bchm_hcorrlen(1:nlev_MT,VarIndex),nulbgst,INI,INJ,INK,-1,'HCORRLEN',-1,-1,-1,'X',bchm_varNameList(VarIndex))
      if (ierr < 0 ) then
        lfound=.false.
        exit
      end if     
    end do
    
    if (lfound) return

    do latIndex = 1, nj_l
      hdist(latIndex)=ec_ra*acos(zrmu(latIndex))
    end do

    zleg(:,:)=0.0d0
    wtemp(:,:,:)=0.0d0
    bchm_hcorrel(:,:,:)=0.0d0
    bchm_hcorrlen(:,:)=0.0
    
    do latIndex = mpi_myid+1, nj_l, mpi_nprocs
       do jn = 0, ntrunc
         wtemp(jn,latIndex,1) = gst_getzleg(jn,latIndex,gstID)
       end do
    end do
    
    nsize=nj_l*(ntrunc+1)    
    call rpn_comm_allreduce(wtemp(0:ntrunc,1:nj_l,1),zleg(0:ntrunc,1:nj_l),nsize,"mpi_double_precision","mpi_sum","GRID",ierr)

    deallocate(wtemp)
    allocate(wtemp(nj_l, nlev_T, numvar3d+numvar2d))
    wtemp(:,:,:)=0.0
    
    VarIndex = 1
    levelIndex = 1
    do jk = 1, nkgdim
      if (jk == nsposit(VarIndex+1)) then
         VarIndex = VarIndex+1 
         levelIndex = 1
      end if

      do latIndex = mpi_myid+1, nj_l, mpi_nprocs
         do jn = 0, ntrunc
            wtemp(latIndex,levelIndex,VarIndex) = wtemp(latIndex,levelIndex,VarIndex)+rstddev(jk,jn)*rstddev(jk,jn)*  &
                                     sqrt(2.0)*sqrt(2.0*jn+1.0)*zleg(jn,latIndex)
         end do       
      end do
      levelIndex = levelIndex+1
    end do
    
    nsize=nj_l*nkgdim   
    call rpn_comm_allreduce(wtemp,bchm_hcorrel,nsize,"mpi_double_precision","mpi_sum","GRID",ierr)
    deallocate(wtemp)
    
    if ( mpi_myid == 0 ) then

      VarIndex = 1
      levelIndex = 1
      do jk = 1, nkgdim
        if (jk == nsposit(VarIndex+1)) then
          VarIndex = VarIndex+1 
          levelIndex = 1
        end if
        do latIndex=nj_l-1,2,-1
          if (bchm_hcorrel(latIndex,levelIndex,VarIndex) <= 0.368) then     ! 1/e ~ 0.368
            bchm_hcorrlen(levelIndex,VarIndex) = (hdist(latIndex)*(bchm_hcorrel(latIndex+1,levelIndex,VarIndex)-0.368) &
               + hdist(latIndex+1)*(0.368-bchm_hcorrel(latIndex,levelIndex,VarIndex))) &
               /(bchm_hcorrel(latIndex+1,levelIndex,VarIndex)-bchm_hcorrel(latIndex,levelIndex,VarIndex))
            exit
          end if
        end do
        levelIndex = levelIndex+1
      end do  
    
      nulcorns = 0
      ierr = fnom(nulcorns,'bmatrixchem_out.fst','STD+RND',0)
      ierr = fstouv(nulcorns,'RND')

      do VarIndex = 1, numvar3d+numvar2d
        nlev_MT = nsposit(VarIndex+1)-nsposit(VarIndex)
        do levelIndex = 1, nlev_MT
          ierr = utl_fstecr(bchm_hcorrel(:,levelIndex,VarIndex),-32,nulcorns,0,0,0,1,nj_l,1,levelIndex,0,nlev_MT,'X',bchm_varNameList(VarIndex), &
                 'HCORREL','X',0,0,0,0,5,.true.)
        end do
        ierr = utl_fstecr(bchm_hcorrlen(1:nlev_MT,VarIndex),-32,nulcorns,0,0,0,1,1,nlev_MT,0,0,0,'X',bchm_varNameList(VarIndex), &
                'HCORRLEN','X',0,0,0,0,5,.true.)
        ierr = utl_fstecr(hdist(1:nj_l),-32,nulcorns,0,0,0,1,nj_l,1,0,0,0,'X',bchm_varNameList(VarIndex), &
                'HDIST','X',0,0,0,0,5,.true.)
      end do
      
      ierr = fstfrm(nulcorns)  
      ierr = fclos(nulcorns)

      write(*,*)
      write(*,*) 'bchm_convol: Horizontal correlations'
      write(*,*)
      write(*,*) 'Separation distances (km)'
      write(*,*) nj_l-nj_l*4/5+1
      write(*,*) hdist(nj_l*4/5:nj_l)/1000.00
      write(*,*)
      do VarIndex = 1, numvar3d+numvar2d
        if (VarIndex <= numvar3d) then
          write(*,*) bchm_varNameList(VarIndex), nlev_T
          do levelIndex = 1, nlev_T
            write(*,'(i3,f10.2,3000f6.2)') levelIndex,bchm_hcorrlen(levelIndex,VarIndex)/1000.00,bchm_hcorrel(nj_l*4/5:nj_l,levelIndex,VarIndex)
          end do
        else
          write(*,*) bchm_varNameList(VarIndex), 1
          write(*,'(i3,f10.2,3000f6.2)') 1,bchm_hcorrlen(1,VarIndex)/1000.00,bchm_hcorrel(nj_l*4/5:nj_l,1,VarIndex)
        end if
        write(*,*)
      end do
    endif

  end subroutine bchm_convol

  !--------------------------------------------------------------------------
  ! bchm_rdstddev
  !--------------------------------------------------------------------------
  subroutine bchm_rdstddev
    !
    !:Purpose: To read stddev and to set as 3D fields.
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
      ! Assume same as rlat 
      if (allocated(rlatr)) deallocate(rlatr)
      allocate(rlatr(nj_l))
      inj = nj_l
      rlatr(1:inj) = rlat(1:inj) 
    end if    

    ikey = utl_fstlir(rcoord,nulbgst,ini,inj,ink, &
                         idate(1),' ',ip1,ip2,ip3,' ','>>')
                         
    if (ikey >= 0) then
      if (allocated(rlongr)) deallocate(rlongr)
      allocate(rlongr(ini+1))
      rlongr(1:ini) = rcoord(1:ini) 
      rlongr(1:ini) = rlongr(1:ini)*MPC_RADIANS_PER_DEGREE_R8
    else if (stddevMode /= 'SP2D') then
      ! Assume same as rlong
      if (allocated(rlongr)) deallocate(rlongr)
      allocate(rlongr(ni_l+1))
      ini = ni_l
      rlongr(1:ini) = rlong(1:ini) 
    end if 
    rlongr(ini+1) = 360.*MPC_RADIANS_PER_DEGREE_R8
     
    ! Read specified input type for error std. dev.
    
    if(stddevMode == 'GD3D') then
      call bchm_rdstd3D
    elseif(stddevMode == 'GD2D') then
      call bchm_rdstd
    elseif(stddevMode == 'SP2D') then
      call bchm_rdspstd
    else
      call utl_abort('bchm_RDSTDDEV: unknown stddevMode')
    endif
    
  end subroutine bchm_rdstddev

  !--------------------------------------------------------------------------
  ! bchm_rdspstd
  !--------------------------------------------------------------------------
  subroutine bchm_rdspstd
    implicit none

    !Locals
    integer :: VarIndex,jn,inix,injx,inkx
    integer :: ikey, levelIndexo, firstn,lastn
    real(8) :: zsp(0:ntrunc,max(nlev_M,nlev_T)),zspbuf(max(nlev_M,nlev_T))
    real(8) :: zgr(nj_l,max(nlev_M,nlev_T))
    
    ! standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3
    integer :: idate(100),nlev_MT
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fstinf

    rgsig(:,:,:) = 0.0d0
    
    ! Reading the Legendre poly representation of the 2D background error std. dev. field

    idate(1) = -1
    ip1      = -1
    ip2      = -1
    ip3      = -1

    cletiket = 'SPSTDDEV'
    cltypvar = 'X'
    
    do VarIndex = 1,numvar3d+numvar2d
      clnomvar = bchm_varNameList(VarIndex)
      nlev_MT = nsposit(VarIndex+1)-nsposit(VarIndex)
      
      firstn = -1
      do jn = 0, ntrunc
        ip2 = jn
        ikey = fstinf(nulbgst,inix,injx,inkx,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)

        if(ikey >= 0 ) then
          ikey = utl_fstlir(zspbuf(1:nlev_MT),nulbgst,ini,inj,ink,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
        else
          if(firstn == -1) firstn = jn
          lastn = jn
          zspbuf(:) = 0.0d0
        endif

        if (ini /= nlev_MT) then
          if (ini == nlev_MT-1) then
            zspbuf(nlev_MT) = zspbuf(ini)
            write(*,*) 'WARNING in CHM_RDSPSTD: ini and nlev_MT not same size - ',ini,nlev_MT
          else
            write(*,*) 'JN, INI, nlev_MT, NOMVAR: ',jn,ini,nlev_MT,' ',clnomvar
            call utl_abort('CHM_RDSPSTD: Constituents background stats levels inconsistency')
          end if    
        end if
     
        do levelIndexo = 1, nlev_MT
          zsp(jn,levelIndexo) = zspbuf(levelIndexo)
        enddo
      enddo
      
      if(mpi_myid == 0 .and. firstn /= -1) write(*,*) 'WARNING: CANNOT FIND SPSTD FOR ',clnomvar, &
            ' AT N BETWEEN ',firstn,' AND ',lastn,', SETTING TO ZERO!!!'

      call gst_zleginv(gstID,zgr(:,1:nlev_MT),zsp(:,1:nlev_MT),nlev_MT)

      do jn = 1, ni_l+1
         rgsig(jn,1:nj_l,nsposit(VarIndex):nsposit(VarIndex+1)-1) = zgr(1:nj_l,1:nlev_MT)
      end do

    enddo 

  end subroutine bchm_rdspstd

  !--------------------------------------------------------------------------
  ! bchm_rdspstd_newfmt
  !--------------------------------------------------------------------------
  subroutine bchm_rdspstd_newfmt
    implicit none

    !Locals
    integer :: VarIndex,jn,inix,injx,inkx,ntrunc_file
    integer :: ikey,levelIndexo
    real(8) :: zsp(0:ntrunc,max(nlev_M,nlev_T))
    real(8), allocatable :: zspbuf(:)
    real(8) :: zgr(nj_l,max(nlev_M,nlev_T))

    ! standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3
    integer :: idate(100),nlev_MT
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fstinf

    rgsig(:,:,:) = 0.0d0

    ! Reading the data

    idate(1) = -1
    ip2      = -1
    ip3      = -1

    cletiket = 'SPSTDDEV'
    cltypvar = 'X'

    ! Check if file is old format
    
    ip1 = -1
    clnomvar = bchm_varNameList(1) 
    ikey = fstinf(nulbgst,inix,injx,inkx,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
    write(*,*) 'ini,inj,ink=',inix,injx,inkx
    if(inix > 1) then
      write(*,*) 'bchm_RDSPSTD_NEWFMT: ini>1, SPSTDDEV is in old format, calling bchm_RDSPSTD...'
      call bchm_rdspstd
      return
    endif

    ! write(*,*) 'Reading for 3D and 2D variables'
    
    do VarIndex = 1,numvar3d+numvar2d
      clnomvar = bchm_varNameList(VarIndex)
      nlev_MT = nsposit(VarIndex+1)-nsposit(VarIndex)

      !write(*,*) 'Reading ',clnomvar

      do levelIndexo = 1, nlev_MT
        if (nlev_MT == 1) then
           ip1 = -1
        else if(vnl_varLevelFromVarName(clnomvar) == 'MM') then
          ip1 = vco_anl%ip1_M(levelIndexo)
        else
          ip1 = vco_anl%ip1_T(levelIndexo)
        endif
        
        ikey = fstinf(nulbgst,inix,ntrunc_file,inkx,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
        ntrunc_file = ntrunc_file-1

        allocate(zspbuf(0:ntrunc_file))
        
        if(ikey >= 0 ) then
          ikey = utl_fstlir(zspbuf(0:ntrunc_file),nulbgst,ini,inj,ink,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
        else
          write(*,*) 'bchm_RDSPSTD_NEWFMT: ',VarIndex,clnomvar,nlev_MT,levelIndexo,ikey,ntrunc,ntrunc_file
          call utl_abort('bchm_RDSPSTD_NEWFMT: SPSTDDEV record not found')
        endif

        zsp(:,levelIndexo) = 0.0d0
        do jn = 0, min(ntrunc,ntrunc_file)
          zsp(jn,levelIndexo) = zspbuf(jn)
        enddo
        deallocate(zspbuf)
      enddo

      call gst_zleginv(gstID,zgr(:,1:nlev_MT),zsp(:,1:nlev_MT),nlev_MT)

      do jn = 1, ni_l+1
         rgsig(jn,1:nj_l,nsposit(VarIndex):nsposit(VarIndex+1)-1) = zgr(1:nj_l,1:nlev_MT)
      end do

    enddo

  end subroutine bchm_rdspstd_newfmt

  !--------------------------------------------------------------------------
  ! bchm_rdstd
  !--------------------------------------------------------------------------
  subroutine bchm_rdstd
    !
    !:Purpose: To read 2D stddev and to store as 3D
    !
    implicit none

    !Locals
    integer :: VarIndex,in
    integer :: ikey
    real(8), allocatable :: rgsig3d(:,:,:)
    
    ! standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3
    integer :: idate(100),nlev_MT
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fstinf
    
    real(8), allocatable  :: vlev(:),vlevout(:)

    rgsig(:,:,:) = 0.0d0

    ! Reading the data

    idate(1) = -1
    ip1      = -1
    ip2      = -1
    ip3      = -1
    
    cletiket = 'STDDEV'
    cltypvar=' '

    ! Reading for 3D and 2D variables
    
    do VarIndex = 1,numvar3d+numvar2d
      clnomvar = bchm_varNameList(VarIndex)
      nlev_MT = nsposit(VarIndex+1)-nsposit(VarIndex)

      !write(*,*) 'Reading ',clnomvar
 
      ikey = fstinf(nulbgst,ini,inj,ink,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
      if (ikey < 0 .or. ini > 1 .or. ink /= nlev_MT) then
        write(*,*) 'bchm_RDSTD: ',VarIndex,clnomvar,ikey,ini,ink,nlev_MT
        call utl_abort(': bchm_RDSTD record not found or incorrect')          
      end if
      
      allocate(rgsig3d(1,inj,ink))
      rgsig3d(1,:,:) = 0.0D0
      allocate(vlev(ink),vlevout(nlev_MT))
      vlev(:) = 1.0D0
      vlevout(:) = 1.0D0
      
      ikey = utl_fstlir(rgsig3d(1,:,:), nulbgst, ini, inj, ink, &
                         idate(1), cletiket, ip1, ip2, ip3, cltypvar, clnomvar)

      if (ikey < 0) then
        write(*,*) 'bchm_RDSTD: ',VarIndex,clnomvar,nlev_MT,ikey
        call utl_abort(': bchm_RDSTD record not found')
      endif
        
      ! Extend to 3D
      if (inj == nj_l) then
        do in = 1, ni_l+1
          rgsig(in,:,nsposit(VarIndex):nsposit(VarIndex+1)-1) = rgsig3d(1,:,:) 
        end do
      else
         ! Interpolate in lat
         call gsv_field3d_hbilin(rgsig3d, 1, inj, ink, rlongr, rlatr, vlev, &
              rgsig(:,:,nsposit(VarIndex):nsposit(VarIndex+1)-1), ni_l+1, nj_l, nlev_MT, &
              rlong, rlat, vlevout)
      end if
      deallocate(rgsig3d, vlev, vlevout)
       
    enddo

  end subroutine bchm_rdstd

  !--------------------------------------------------------------------------
  ! bchm_rdstd3d
  !--------------------------------------------------------------------------
  subroutine bchm_rdstd3d
    !
    !:Purpose: To read 3D stddev.
    !

    ! originally based on bchm_rdspstd_newfmt
    !
    implicit none

    !Locals
    integer :: VarIndex
    integer :: ikey,levelIndexo
    real(8), allocatable :: rgsig3d(:,:,:)
    
    ! standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3
    integer :: idate(100),nlev_MT
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fstinf

    real(8) :: vlev(1),vlevout(1)

    rgsig(:,:,:) = 0.0d0
    vlev(:) = 1.0D0
    vlevout(:) = 1.0D0

    ! Reading the data

    idate(1) = -1
    ip2      = -1
    ip3      = -1
    
    cletiket = 'STDDEV3D'
    cltypvar=' '

    ! Reading for 3D and 2D variables
    
    do VarIndex = 1,numvar3d+numvar2d
      clnomvar = bchm_varNameList(VarIndex)
      nlev_MT = nsposit(VarIndex+1)-nsposit(VarIndex)

      !write(*,*) 'Reading ',clnomvar

      do levelIndexo = 1, nlev_MT
        if (nlev_MT == 1) then
          ip1 = -1
        else if(vnl_varLevelFromVarName(clnomvar) == 'MM') then
          ip1 = vco_anl%ip1_M(levelIndexo)
        else
          ip1 = vco_anl%ip1_T(levelIndexo)
        endif
        
        ikey = fstinf(nulbgst,ini,inj,ink,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
        if (ikey < 0) then
            write(*,*) 'bchm_RDSTD: ',VarIndex,clnomvar,ikey,levelIndexo
            call utl_abort(': bchm_RDSTD record not foun0d')          
        end if
      
        allocate(rgsig3d(ini+1, inj, 1))
        rgsig3d(:,:,1) = 0.0D0
        
        ikey = utl_fstlir(rgsig3d(1:ini,:,1), nulbgst, ini, inj, ink, &
                         idate(1), cletiket, ip1, ip2, ip3, cltypvar, clnomvar)
        if (ikey < 0) then
          write(*,*) 'bchm_RDSTD3D: ',VarIndex,clnomvar,nlev_MT,levelIndexo,ikey,ip1
          call utl_abort(': bchm_RDSTD3D record not found')
        endif
        
        ! Move to rgsig
        if (inj == nj_l .and. ini == ni_l) then
          rgsig(1:ni_l,:,nsposit(VarIndex)+(levelIndexo-1)) = rgsig3d(1:ni_l,:,1) 
          rgsig(ni_l+1,:,nsposit(VarIndex)+(levelIndexo-1)) = rgsig3d(1,:,1)
        else
          ! Interpolate in lat and long
          rgsig3d(ini+1,:,1) = rgsig3d(1,:,1)
          call gsv_field3d_hbilin(rgsig3d, ini+1, inj, 1, rlongr, rlatr, vlev, &
               rgsig(:,:,nsposit(VarIndex)+(levelIndexo-1)), ni_l+1, nj_l, 1, &
               rlong, rlat, vlevout)
       end if
       ! write(*,*) 'STDDDEV ',levelIndexo,rgsig(1,1,nsposit(VarIndex)+(levelIndexo-1)),rgsig(ni_l+1,1,nsposit(VarIndex)+(levelIndexo-1)),rgsig(ni_l/2,nj_l/2,nsposit(VarIndex)+(levelIndexo-1)),rgsig(ni_l+1,nj_l,nsposit(VarIndex)+(levelIndexo-1))
       deallocate(rgsig3d)

      enddo
    enddo

  end subroutine bchm_rdstd3d

  !--------------------------------------------------------------------------
  ! bchm_sucorns2
  !--------------------------------------------------------------------------
  subroutine bchm_sucorns2
    implicit none

    !Locals
    real(8) :: eigenval(nkgdim)
    real(8) :: eigenvalsqrt(nkgdim)
    real(8), allocatable :: eigenvec(:,:),result(:,:)

    integer :: jn,jk1,jk2,VarIndex,ierr
    integer :: ilwork,info,jnum,jstart,nsize

    real(8) :: zwork(2*4*nkgdim)
    real(8) :: ztlen,zcorr,zr,zpres1,zpres2,eigenvalmax
    real(8), allocatable :: corns_temp(:,:,:)
    logical, allocatable :: lfound_sqrt(:)
    
    ! Apply vertical localization to correlations of 3D fields.
    ! Currently assumes no-cross correlations for variables (block diagonal matrix)
    
    do VarIndex = 1, numvar3d
      ztlen = rvloc(VarIndex)    ! specify length scale (in units of ln(Pressure))
        
      jstart = nsposit(VarIndex)
      jnum = nsposit(VarIndex+1)-nsposit(VarIndex)
       
      if(ztlen > 0.0d0) then
        ! calculate 5'th order function (from Gaspari and Cohn)
        do jk1 = 1, jnum
          zpres1 = log(pressureProfile_T(jk1))
          do jk2 = 1, jnum
            zpres2 = log(pressureProfile_T(jk2))
            zr = abs(zpres2 - zpres1)
            zcorr = gasparicohn(ztlen,zr)
            corns(jstart-1+jk1,jstart-1+jk2,0:ntrunc) = &
                  corns(jstart-1+jk1,jstart-1+jk2,0:ntrunc)*zcorr  
          enddo
        enddo
      endif
    enddo

    ! Compute total vertical correlations and its inverse (currently for each block matrix).
    
    call bchm_corvert_setup
    
    if (trim(bchm_mode) == 'BackgroundCheck') return

    allocate(lfound_sqrt(numvar3d+numvar2d))
    if(ReadWrite_sqrt) then
      ! if desired, read precomputed sqrt of corns
      call readcorns(lfound_sqrt,ntrunc,'CORNS_SQRT')
    else
      lfound_sqrt(:)=.false.
    endif

    do VarIndex = 1, numvar3d+numvar2d

      if (any(CrossCornsVarKindCH(:) /= '')) then
         jstart=1
         jnum=nkgdim
      else
         jstart = nsposit(VarIndex)
         jnum = nsposit(VarIndex+1)-nsposit(VarIndex)
      end if
      
      if (.not.lfound_sqrt(VarIndex)) then
        
         ! compute square-root of corns for each total wavenumber
    
         allocate(corns_temp(jnum,jnum,0:ntrunc),eigenvec(jnum,jnum),result(jnum,jnum))
         corns_temp(:,:,:) = 0.0d0
         do jn = mpi_myid, ntrunc, mpi_nprocs

           eigenvec(1:jnum,1:jnum) = corns(jstart:jstart-1+jnum,jstart:jstart-1+jnum,jn)

           ! CALCULATE EIGENVALUES AND EIGENVECTORS.
           ilwork = 4*jnum*2
           call dsyev('V','U',jnum,eigenvec,jnum,eigenval,zwork,ilwork,info)
           if(info /= 0) then
              write(*,*) 'bchm_sucorns2: non-zero value of info =',info,' returned by dsyev for wavenumber ',jn,VarIndex
              call utl_abort('bchm_SUCORNS')
           endif
      
           ! set selected number of eigenmodes to zero
           if(numModeZero > 0) then
             ! write(*,*) 'bchm_sucorns2: setting ',numModeZero,' eigenvalues to zero for wavenumber n=',jn
             ! write(*,*) 'bchm_sucorns2: original eigenvalues=',eigenval(:)
             do jk1 = 1, numModeZero
               eigenval(jk1) = 0.0d0
             enddo
             ! write(*,*) 'bchm_sucorns2: modified eigenvalues=',eigenval(:)
           endif

           eigenvalmax = maxval(eigenval(1:jnum))
           do jk1 = 1, jnum
             ! if(eigenval(jk1) < 1.0d-15) then
             if(eigenval(jk1) < 1.0d-8*eigenvalmax) then
               eigenvalsqrt(jk1) = 0.0d0
             else
               eigenvalsqrt(jk1) = sqrt(eigenval(jk1))
             endif
           enddo

           ! E * lambda^1/2
           do jk1 = 1, jnum
              result(1:jnum,jk1) = eigenvec(1:jnum,jk1)*eigenvalsqrt(jk1)
           enddo
  
           ! (E * lambda^1/2) * E^T
           do jk1 = 1, jnum
           do jk2 = 1, jnum
              corns_temp(1:jnum,jk1,jn) = corns_temp(1:jnum,jk1,jn) + result(1:jnum,jk2)*eigenvec(jk1,jk2)
           enddo
           enddo

         enddo ! jn
  
         nsize = jnum*jnum*(ntrunc+1)
         call rpn_comm_allreduce(corns_temp,corns(jstart:jstart+jnum-1,jstart:jstart+jnum-1,:),nsize,"mpi_double_precision","mpi_sum","GRID",ierr)
         deallocate(corns_temp,eigenvec,result)

      end if
      
      if (any(CrossCornsVarKindCH(:) /= '')) exit

    end do
   
    deallocate(lfound_sqrt)
    if(ReadWrite_sqrt) then
      ! Write computed sqrt to a separate file.
      call writecorns(ntrunc,'CORNS_SQRT',-1)
    endif
    
  end subroutine bchm_sucorns2

  !--------------------------------------------------------------------------
  ! bchm_corvert_setup
  !--------------------------------------------------------------------------
  subroutine bchm_corvert_setup
    !
    !:Purpose: To compute total vertical correlations (bchm_corvert) and its
    !          inverse (bchm_corverti; currently for each block matrix).
    !
    !
    !:Note: Currently assumes no (or neglects) cross-correlations 
    !
    implicit none

    !Locals
    real(8) :: eigenval(nkgdim)
    real(8), allocatable :: eigenvec(:,:),result(:,:)

    integer :: jn,jk1,jk2,VarIndex,numvartot,jstart
    integer :: ilwork,info,jnum,nvlev,ierr

    real(8) :: zwork(2*4*nkgdim)
    real(8) :: eigenvalmax
    integer iulcorvert
    
    ! Standard file variables
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fnom,fstouv,fstfrm,fclos
    logical, allocatable :: lfound(:)
   
    ! Compute total vertical correlations and its inverse (currently for each block matrix).

    if (mpi_myid == 0 .and. ReadWrite_PSStats) then
       iulcorvert = 0
       ierr = fnom(iulcorvert,'bmatrixchem_out.fst','STD+RND',0)
       ierr = fstouv(iulcorvert,'RND')
    end if 
    
    nvlev=-1
    numvartot=numvar3d
    do VarIndex = 1, numvartot
   
      jstart = nsposit(VarIndex)
      jnum = nsposit(VarIndex+1)-nsposit(VarIndex)
       
      bchm_corvert(1:jnum,1:jnum,VarIndex) = 0.0d0
      do jn = 0, ntrunc
        bchm_corvert(1:jnum,1:jnum,VarIndex) = bchm_corvert(1:jnum,1:jnum,VarIndex)+ ((2*jn+1)* &
             corns(nsposit(VarIndex):nsposit(VarIndex+1)-1,nsposit(VarIndex):nsposit(VarIndex+1)-1,jn))
      enddo

      !  where (abs(bchm_corvert(1:jnum,1:jnum,VarIndex)) .lt. 1.E-3) bchm_corvert(1:jnum,1:jnum,VarIndex)=0.0D0
       
      ! Inverse (and possible vertical interpolation to model levels) not needed if not in analysis mode
      if (trim(bchm_mode) == 'BackgroundCheck') cycle
      
      if (mpi_myid /= 0 .or. .not.ReadWrite_PSStats) cycle
         
      if (VarIndex.eq.1) then
        allocate(lfound(numvartot))
        call readcorns(lfound,0,'CORVERTI') 
      end if
      
      if (.not.lfound(VarIndex)) then
      
         write(*,*) 'bchm_corvert_setup: Calculating CORVERT/CORVERTI for VarIndex =',VarIndex

         allocate(eigenvec(jnum,jnum),result(jnum,jnum))
         eigenvec(1:jnum,1:jnum)=bchm_corvert(1:jnum,1:jnum,VarIndex)       

         ! CALCULATE EIGENVALUES AND EIGENVECTORS.
         ilwork = 4*jnum*2
         call dsyev('V','U',jnum,eigenvec,jnum,eigenval,zwork,ilwork,info)
         if (info.ne.0) then
            write(*,*) 'bchm_corvert_setup: non-zero value of info =',info,' returned by dsyev for wavenumber ',jn
            call utl_abort('bchm_corvert_setup')
         endif

         ! Set selected number of eigenmodes to zero
         if (numModeZero > 0) then
           do jk1 = 1, numModeZero
             eigenval(jk1) = 0.0d0
           enddo
           write(*,*) 'bchm_corvert_setup: modified eigenvalues=',eigenval(:)
         endif

         ! E * lambda^{-1}
         eigenvalmax=maxval(eigenval(1:jnum))
         do jk1 = 1, jnum
           if (eigenval(jk1) > 1.0d-8*eigenvalmax) then
             result(1:jnum,jk1) = eigenvec(1:jnum,jk1)/eigenval(jk1)
           else
             result(1:jnum,jk1) = 0.0D0
           end if
         enddo

         ! E * lambda^{-1} * E^T
         bchm_corverti(1:jnum,1:jnum,VarIndex)=0.0D0
         do jk1 = 1, jnum
         do jk2 = 1, jnum
            bchm_corverti(1:jnum,jk1,VarIndex) = bchm_corverti(1:jnum,jk1,VarIndex) + result(1:jnum,jk2)*eigenvec(jk1,jk2)
         enddo
         enddo

         ! zr=maxval(abs(bchm_corverti(1:jnum,1:jnum,VarIndex)))
         ! where (abs(bchm_corverti(1:jnum,1:jnum,VarIndex)) .lt. 1.E-5*zr) bchm_corverti(1:jnum,1:jnum,VarIndex)=0.0D0

         ! Check inverse (for output when mpi_myid is 0)
         result(1:jnum,1:jnum)=0.0D0
         do jk1 = 1, jnum
         do jk2 = 1, jnum
            result(1:jnum,jk1) = result(1:jnum,jk1) + &
              bchm_corvert(1:jnum,jk2,VarIndex)*bchm_corverti(jk2,jk1,VarIndex)
         enddo
         enddo

         cletiket = 'C*C^-1'
         clnomvar = bchm_varNameList(VarIndex)
          
         ierr = utl_fstecr(result(1:jnum,1:jnum),-32,iulcorvert,0,0,0,jnum,jnum,1,  &
                       0,0,ntrunc,'X',clnomvar,cletiket,'X',0,0,0,0,5,.true.)
                         
         deallocate(eigenvec,result)
      end if
      
      ! Generate 1/sum(covert(:,i,VarIndex)
        
      do jk1 = 1, jnum  
        bchm_invsum(jk1,VarIndex) = 1.0D0/sum(bchm_corvert(1:jnum,jk1,VarIndex))
      end do
      
    end do
    
    if (mpi_myid /= 0 .or. .not.ReadWrite_PSStats) return
    
    deallocate(lfound)

    ! Reset bchm_corvert and bchm_invsum if vertical dimensions differ from that of trial fields
        
    nvlev=0
    if (trim(bchm_mode) /= 'BackgroundCheck') call bchm_resetCorvert(nvlev)

    ! Write bchm_invsum
        
    do VarIndex = 1, numvartot
       jnum = nsposit(VarIndex+1)-nsposit(VarIndex)
       if (nvlev.gt.0) jnum=nvlev
       ierr = utl_fstecr(bchm_invsum(1:jnum,VarIndex),-32,iulcorvert,0,0,0,1,1,jnum,  &
                      0,0,0,'X',bchm_varNameList(VarIndex),'VCOR INVSUM','X',0,0,0,0,5,.true.)
    end do
    ierr = fstfrm(iulcorvert)  
    ierr = fclos(iulcorvert)

    call writecorns(0,'CORVERT',nvlev)
    if (any(.not.lfound(1:numvartot))) call writecorns(0,'CORVERTI',-1)

  end subroutine bchm_corvert_setup

  !--------------------------------------------------------------------------
  ! bchm_writecorns
  !--------------------------------------------------------------------------
  subroutine writecorns(nmat,cletiket,nlev)
  
    implicit none

    !Arguments
    character(len=*) :: cletiket
    integer :: nmat,nlev

    !Locals
    integer :: jn, nulcorns,ierr,VarIndex,jstart,jnum,numvartot

    ! standard file variables
    integer :: ip1,ip2,ip3
    integer :: idateo, ipak, idatyp
    character(len=4)  :: clnomvar
    integer :: fnom, fstouv, fstfrm, fclos
    
    if(mpi_myid==0) then

      write(*,*) 'WRITECORNS: ', trim(cletiket), ' being written to file bmatrixchem_out.fst for number of matrices - 1 =', nmat

      nulcorns = 0
      ierr = fnom(nulcorns,'bmatrixchem_out.fst','STD+RND',0)
      ierr = fstouv(nulcorns,'RND')

      ipak = -32
      idatyp = 5
      ip1 = 0
      ip2 = 0
      ip3 = ntrunc
      idateo = 0

      if (nmat == 0) then
         numvartot=numvar3d
      else
         numvartot=numvar3d+numvar2d
      end if
      
      do VarIndex = 1, numvartot
   
        if (any(CrossCornsVarKindCH(:) /= '') .and. nmat > 0) then
           jstart=1
           jnum=nkgdim
           clnomvar = 'ZZ'
        else
           jstart = nsposit(VarIndex)
           jnum = nsposit(VarIndex+1)-nsposit(VarIndex)
           if (nlev.gt.0) jnum=nlev
           clnomvar = bchm_varNameList(VarIndex)
        end if
        
        if (nmat > 0 ) then
          do jn = 0, nmat
            ip2 = jn
            ierr = utl_fstecr(corns(jstart:jstart+jnum-1,jstart:jstart+jnum-1,jn),ipak,nulcorns,idateo,0,0,jnum,jnum,1,  &
                         ip1,ip2,ip3,'X',clnomvar,cletiket,'X',0,0,0,0,idatyp,.true.)
          enddo
        else
          if (trim(cletiket) == 'CORVERT') then
             ierr = utl_fstecr(bchm_corvert(1:jnum,1:jnum,VarIndex),ipak,nulcorns,idateo,0,0,jnum,jnum,1,  &
                          ip1,ip2,ip3,'X',clnomvar,cletiket,'X',0,0,0,0,idatyp,.true.)
          else
             ierr = utl_fstecr(bchm_corverti(1:jnum,1:jnum,VarIndex),ipak,nulcorns,idateo,0,0,jnum,jnum,1,  &
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
  ! bchm_readcorns
  !--------------------------------------------------------------------------
  subroutine readcorns(lfound,nmat,cletiket)

    implicit none

    !Arguments
    logical :: lfound(:)
    integer :: nmat
    character(len=*) :: cletiket

    !Locals
    integer :: jn, icornskey,VarIndex,jnum,jstart,numvartot
    real(8), allocatable :: zcornssrc(:,:)

    ! standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3
    integer :: idateo
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar

    if (nmat == 0) then
       numvartot=numvar3d
       if (trim(cletiket) == 'CORVERT') then
          bchm_corvert(:,:,:)=0.0d0
       else if (trim(cletiket) == 'CORVERTI') then
          bchm_corverti(:,:,:)=0.0d0
       end if
    else if (nmat == ntrunc) then
       numvartot=numvar3d+numvar2d
    end if
    
    write(*,*) 'READCORNS: ', trim(cletiket), ' being searched for number of matrices -1 =',nmat

    idateo = -1
    ip1 = -1
    ip3 = ntrunc
    cltypvar = 'X'

    lfound(:)=.false.    
    VARCYCLE: do VarIndex = 1, numvartot
   
      if (any(CrossCornsVarKindCH(:) /= '') .and. nmat > 0) then
        jstart=1
        jnum=nkgdim
        clnomvar = 'ZZ'
      else
        jstart = nsposit(VarIndex)
        jnum = nsposit(VarIndex+1)-nsposit(VarIndex) 
        clnomvar = bchm_varNameList(VarIndex)
      end if
      allocate(zcornssrc(jnum,jnum))

      do jn = 0, nmat
        ip2 = jn

        ! Looking for FST record parameters..
  
        icornskey = utl_fstlir(ZCORNSSRC,nulbgst,INI,INJ,INK,idateo,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)

        if(icornskey .lt.0 ) then
          write(*,*) 'READCORNS: matrix not found in stats file for variable ', clnomvar
          deallocate(zcornssrc)
          cycle VARCYCLE
        endif

        if (ini /= jnum .or. inj /= jnum) call utl_abort('READCORNS: BG stat levels inconsitencies')

        if (nmat > 0) then
           if (jn == 0) corns(jstart:jstart+jnum-1,jstart:jstart+jnum-1,:)=0.0d0
           corns(jstart:jstart+jnum-1,jstart:jstart+jnum-1,jn)=zcornssrc(1:jnum,1:jnum)
        else if (trim(cletiket) == 'CORVERT') then
           bchm_corvert(1:jnum,1:jnum,VarIndex)=zcornssrc(1:jnum,1:jnum)
        else
           bchm_corverti(1:jnum,1:jnum,VarIndex)=zcornssrc(1:jnum,1:jnum)        
        end if
      enddo
      lfound(VarIndex)=.true.

      deallocate(zcornssrc)

      if (any(CrossCornsVarKindCH(:) /= '') .and. nmat > 0) exit

    end do VARCYCLE
    
  end subroutine readcorns

  !--------------------------------------------------------------------------
  ! gaspariCohn
  !--------------------------------------------------------------------------
  function gaspariCohn(ztlen,zr)

    !Arguments
    real(8)  :: gasparicohn
    real(8)  :: ztlen,zr
    
    !Locals
    real(8)  :: zlc

    zlc = ztlen/2.0d0
    if(zr.le.zlc) then
      gasparicohn = -0.250d0*(zr/zlc)**5+0.5d0*(zr/zlc)**4             &
                  +0.625d0*(zr/zlc)**3-(5.0d0/3.0d0)*(zr/zlc)**2+1.0d0
    elseif(zr.le.(2.0d0*zlc)) then
      gasparicohn = (1.0d0/12.0d0)*(zr/zlc)**5-0.5d0*(zr/zlc)**4         &
                  +0.625d0*(zr/zlc)**3+(5.0d0/3.0d0)*(zr/zlc)**2       &
                  -5.0d0*(zr/zlc)+4.0d0-(2.0d0/3.0d0)*(zlc/zr)
    else
      gasparicohn = 0.0d0
    endif
    if(gasparicohn.lt.0.0d0) gasparicohn = 0.0d0

  end function gaspariCohn

  !--------------------------------------------------------------------------
  ! bchm_truncateCV
  !--------------------------------------------------------------------------
  subroutine bchm_truncateCV(controlVector_inout,ntruncCut)
    !
    !:Purpose: To set to zero all coefficients with total wavenumber > ntruncCut
    !

    ! Based on bhi_truncateCV.
    implicit none

    !Arguments
    real(8), pointer :: controlVector_inout(:)
    integer          :: ntruncCut
    
    !Locals
    integer          :: jn, jm, ila_mpiglobal, ila_mpilocal, levelIndex, jdim

    if(.not. initialized) then
      if(mpi_myid == 0) write(*,*) 'bchm_truncateCV: bMatrixChem not initialized'
      return
    endif

    if(ntruncCut.gt.ntrunc) then
      write(*,*) ntruncCut, ntrunc
      call utl_abort('bchm_truncateCV: ntruncCut is greater than ntrunc!')
    endif

    jdim = 0
    do levelIndex = 1, nkgdim
      do jm = mymBeg, mymEnd, mymSkip
        do jn = mynBeg, mynEnd, mynSkip
          if(jm.le.jn) then
            ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
            ila_mpilocal  = ilaList_mpilocal(ila_mpiglobal)
            if(jm == 0) then
              ! only real component for jm=0
              jdim = jdim + 1
              if(jn.gt.ntruncCut) controlVector_inout(jdim) = 0.0d0
            else
              ! both real and imaginary components for jm>0
              jdim = jdim + 1
              if(jn.gt.ntruncCut) controlVector_inout(jdim) = 0.0d0
              jdim = jdim + 1
              if(jn.gt.ntruncCut) controlVector_inout(jdim) = 0.0d0
            endif
          endif
        enddo
      enddo
   enddo
   
  end subroutine bchm_truncateCV

  !--------------------------------------------------------------------------
  ! bchm_bSqrt
  !--------------------------------------------------------------------------
  subroutine bchm_bSqrt(controlvector_in,statevector, stateVectorRef_opt)
    !
    !:Purpose: To apply B_CHM^1/2 to a control vector.
    !

    ! Based on bhi_bsqrt
    implicit none

    !Arguments
    real(8)   :: controlVector_in(cvDim_mpilocal)
    type(struct_gsv) :: statevector
    type(struct_gsv), optional :: stateVectorRef_opt
    
    !Locals
    real(8),allocatable :: gd_out(:,:,:)
    real(8)   :: hiControlVector(nla_mpilocal,2,nkgdim)
    character(len=30) :: transform
    integer :: varIndex
    
    if(.not. initialized) return

    if(mpi_myid == 0) write(*,*) 'bchm_bsqrt: starting'
    if(mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    allocate(gd_out(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim))

    call bchm_cain(controlVector_in,hiControlVector)

    call bchm_spa2gd(hiControlVector,gd_out)
    
    call bchm_copyToStatevector(statevector,gd_out)

    if ( trim(transformVarKindCH) /= '' ) then  
     
      transform = trim(transformVarKindCH)//'CH_tlm'
      do varIndex= 1, numvar3d+numvar2d
   
        if (vnl_varKindFromVarname(vnl_varNameList(varIndex)) /= 'CH') cycle
    
        if ( present(stateVectorRef_opt) ) then
          call gvt_transform( statevector,  &                          ! INOUT
                              trim(transform), &                       ! IN
                              stateVectorRef_opt=stateVectorRef_opt, varName_opt=bchm_varNameList(varIndex) ) ! IN
        else
          call gvt_transform( statevector,  &                          ! INOUT
                              trim(transform), &                       ! IN
                              varName_opt=bchm_varNameList(varIndex) ) ! IN
        end if

      end do
    end if
    
    deallocate(gd_out)

    if(mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if(mpi_myid == 0) write(*,*) 'bchm_bsqrt: done'

  end subroutine bchm_bSqrt

  !--------------------------------------------------------------------------
  ! bchm_bSqrtAd
  !--------------------------------------------------------------------------
  subroutine bchm_bSqrtAd(statevector,controlVector_out, stateVectorRef_opt)
    !
    !:Purpose: To apply adjoint of B_CHM^1/2 to a statevector.
    !

    ! Based on bhi_bSqrtAd.
    !
    implicit none

    !Arguments
    real(8)   :: controlVector_out(cvDim_mpilocal)
    type(struct_gsv) :: statevector
    type(struct_gsv), optional :: stateVectorRef_opt
    
    !Locals
    real(8), allocatable :: gd_in(:,:,:)
    real(8)   :: hiControlVector(nla_mpilocal,2,nkgdim)
    character(len=30) :: transform
    integer :: varIndex  

    if(.not. initialized) then
      if(mpi_myid == 0) write(*,*) 'bMatrixChem not initialized'
      return
    endif

    if(mpi_myid == 0) write(*,*) 'bchm_bsqrtad: starting'
    if(mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    allocate(gd_in(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim))

    if ( trim(transformVarKindCH) /= '' ) then  
          
      transform = trim(transformVarKindCH)//'CH_ad'            
      do varIndex = 1, numvar3d+numvar2d
   
        if (vnl_varKindFromVarname(vnl_varNameList(varIndex)) /= 'CH') cycle
    
        if ( present(stateVectorRef_opt) ) then
          call gvt_transform( statevector,  &                          ! INOUT
                              trim(transform), &                       ! IN
                              stateVectorRef_opt=stateVectorRef_opt, varName_opt=bchm_varNameList(varIndex) ) ! IN
        else
          call gvt_transform( statevector,  &                          ! INOUT
                              trim(transform), &                       ! IN
                              varName_opt=bchm_varNameList(varIndex) ) ! IN
        end if

      end do
    end if

    call bchm_copyFromStatevector(statevector,gd_in)

    call bchm_spa2gdad(gd_in,hiControlVector)

    call bchm_cainad(hiControlVector,controlVector_out)

    deallocate(gd_in)

    if(mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if(mpi_myid == 0) write(*,*) 'bchm_bsqrtad: done'

  end subroutine bchm_bSqrtAd

  !--------------------------------------------------------------------------
  ! bchm_cain
  !--------------------------------------------------------------------------
  subroutine bchm_cain(controlVector_in,hiControlVector_out)
    !
    implicit none

    !Arguments
    real(8) :: controlVector_in(cvDim_mpilocal)
    real(8) :: hiControlVector_out(nla_mpilocal,2,nkgdim)

    !Locals
    integer :: jdim, levelIndex, jm, jn, ila_mpilocal, ila_mpiglobal

    jdim = 0
    hiControlVector_out(:,:,:) = 0.0d0
    do levelIndex = 1, nkgdim
      do jm = mymBeg, mymEnd, mymSkip
        do jn = mynBeg, mynEnd, mynSkip
          if(jm.le.jn) then
            ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
            ila_mpilocal  = ilaList_mpilocal(ila_mpiglobal)
            if(jm == 0) then
              ! only real component for jm=0
              jdim = jdim + 1
              hiControlVector_out(ila_mpilocal,1,levelIndex) = controlVector_in(jdim)
            else
              ! both real and imaginary components for jm>0
              jdim = jdim + 1
              hiControlVector_out(ila_mpilocal,1,levelIndex) = controlVector_in(jdim)
              jdim = jdim + 1
              hiControlVector_out(ila_mpilocal,2,levelIndex) = controlVector_in(jdim)
            endif
          endif
        enddo
      enddo
    enddo

  end subroutine bchm_cain

  !--------------------------------------------------------------------------
  ! bchm_cainAd
  !--------------------------------------------------------------------------
  subroutine bchm_cainAd(hiControlVector_in,controlVector_out)
    implicit none

    !Arguments
    real(8) :: controlVector_out(cvDim_mpilocal)
    real(8) :: hiControlVector_in(nla_mpilocal,2,nkgdim)

    !Locals
    integer :: jdim, levelIndex, jm, jn, ila_mpilocal, ila_mpiglobal

    jdim = 0
    do levelIndex = 1, nkgdim
      do jm = mymBeg, mymEnd, mymSkip
        do jn = mynBeg, mynEnd, mynSkip
          if(jm.le.jn) then
            ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
            ila_mpilocal  = ilaList_mpilocal(ila_mpiglobal)
            if(jm == 0) then
              ! only real component for jm=0
              jdim = jdim + 1
              controlVector_out(jdim) = controlVector_out(jdim) + hiControlVector_in(ila_mpilocal,1,levelIndex)
            else
              ! both real and imaginary components for jm>0
              jdim = jdim + 1
              controlVector_out(jdim) = controlVector_out(jdim) + hiControlVector_in(ila_mpilocal,1,levelIndex)*2.0d0
              jdim = jdim + 1
              controlVector_out(jdim) = controlVector_out(jdim) + hiControlVector_in(ila_mpilocal,2,levelIndex)*2.0d0
            endif
          endif
        enddo
      enddo
    enddo

  end subroutine bchm_cainAd

  !--------------------------------------------------------------------------
  ! bchm_spa2gd
  !--------------------------------------------------------------------------
  subroutine bchm_spa2gd(hiControlVector_in,gd_out)
    implicit none

    !Arguments
    real(8) :: hiControlVector_in(nla_mpilocal,2,nkgdim)
    real(8) :: gd_out(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)

    !Locals
    real(8) :: sp(nla_mpilocal,2,nkgdim)

    integer :: jn,jm,ila_mpilocal,ila_mpiglobal,icount
    real(8) :: sq2
    real(8) , allocatable :: zsp(:,:,:), zsp2(:,:,:)
    integer :: levelIndex, lonIndex, latIndex
    real(8), target  :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)

    ! maybe not needed:
    sp(:,:,:) = 0.0d0
    sq2 = sqrt(2.0d0)

    allocate(zsp(nkgdim,2,mymCount))
    allocate(zsp2(nkgdim,2,mymCount))

!$OMP PARALLEL DO PRIVATE(jn,jm,levelIndex,ila_mpiglobal,ila_mpilocal,zsp2,zsp,icount)
    do jn = mynBeg, mynEnd, mynSkip

      icount = 0
      do jm = mymBeg, mymEnd, mymSkip
        if(jm.le.jn) then
          icount = icount+1
          ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
          ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
          do levelIndex = 1, nkgdim
            zsp(levelIndex,1,icount) = hiControlVector_in(ila_mpilocal,1,levelIndex)
            zsp(levelIndex,2,icount) = hiControlVector_in(ila_mpilocal,2,levelIndex)
          enddo
        endif
      enddo
      if(icount.gt.0) then

        CALL DGEMM('N','N',nkgdim,2*icount,nkgdim,1.0d0,corns(1,1,jn),nkgdim,zsp(1,1,1),nkgdim,0.0d0,zsp2(1,1,1),nkgdim)
        ! CALL DGEMUL(corns(1,1,jn),nkgdim,'N',zsp(1,1,1),nkgdim,'N',zsp2(1,1,1),nkgdim,nkgdim,nkgdim,2*icount)

        icount = 0
        do jm = mymBeg, mymEnd, mymSkip
          if(jm.le.jn) then
            icount = icount+1
            ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
            ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
            do levelIndex = 1, nkgdim
              sp(ila_mpilocal,1,levelIndex) = zsp2(levelIndex,1,icount)
              sp(ila_mpilocal,2,levelIndex) = zsp2(levelIndex,2,icount)
            enddo
          endif
        enddo

      endif

      ! make adjustments for jm=0
      if(mymBeg == 0) then

        ila_mpiglobal = gst_getNind(0,gstID) + jn
        ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)

        do levelIndex = 1, nkgdim
          sp(ila_mpilocal,1,levelIndex) = sp(ila_mpilocal,1,levelIndex)*sq2
          sp(ila_mpilocal,2,levelIndex) = 0.0d0
        enddo

      endif

    enddo
!$OMP END PARALLEL DO
    deallocate(zsp)
    deallocate(zsp2)

!$OMP PARALLEL DO PRIVATE(latIndex,levelIndex,lonIndex)
    do levelIndex = 1, nkgdim
       do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
             gd(lonIndex,latIndex,levelIndex) = 0.0d0
          enddo
       enddo
    enddo
!$OMP END PARALLEL DO

    call gst_setID(gstID)
    call gst_speree(sp,gd)
    call gst_setID(gstID2)

!$OMP PARALLEL DO PRIVATE(latIndex,levelIndex,lonIndex)
    do levelIndex = 1, nkgdim
       do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
             gd(lonIndex,latIndex,levelIndex) = gd(lonIndex,latIndex,levelIndex)*rgsig(lonIndex,latIndex,levelIndex)
          enddo
       enddo
    enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(latIndex,levelIndex,lonIndex)
    do levelIndex = 1, nkgdim
       do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
             gd_out(lonIndex,latIndex,levelIndex) = gd(lonIndex,latIndex,levelIndex)
          enddo
       enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine bchm_spa2gd

  !--------------------------------------------------------------------------
  ! bchm_spa2gdad
  !--------------------------------------------------------------------------
  subroutine bchm_spa2gdad(gd_in,hiControlVector_out)
    implicit none

    !Arguments
    real(8) :: hiControlVector_out(nla_mpilocal,2,nkgdim)
    real(8) :: gd_in(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)

    !Locals
    real(8) :: sp(nla_mpilocal,2,nkgdim)

    integer :: jn, jm, ila_mpilocal, ila_mpiglobal, icount
    real(8) :: sq2
    real(8), allocatable :: zsp(:,:,:), zsp2(:,:,:)

    integer :: levelIndex, lonIndex, latIndex
    real(8), target :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)


!$OMP PARALLEL DO PRIVATE(latIndex,levelIndex,lonIndex)
    do levelIndex = 1, nkgdim
       do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd                                                      
             gd(lonIndex,latIndex,levelIndex) = gd_in(lonIndex,latIndex,levelIndex)
          enddo
       enddo
    enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(latIndex,levelIndex,lonIndex)
   do levelIndex = 1, nkgdim
      do latIndex = myLatBeg, myLatEnd
         do lonIndex = myLonBeg, myLonEnd
            gd(lonIndex,latIndex,levelIndex) = gd(lonIndex,latIndex,levelIndex)*rgsig(lonIndex,latIndex,levelIndex)
         enddo
      enddo
   enddo
!$OMP END PARALLEL DO

    call gst_setID(gstID)
    call gst_speree_ad(sp,gd)

    hiControlVector_out(:,:,:) = 0.0d0
    sq2 = sqrt(2.0d0)
    allocate(zsp(nkgdim,2,mymCount))
    allocate(zsp2(nkgdim,2,mymCount))
    
!$OMP PARALLEL DO PRIVATE(JN,JM,levelIndex,ILA_MPILOCAL,ILA_MPIGLOBAL,zsp,zsp2,icount)
    do jn = mynBeg, mynEnd, mynSkip

      icount = 0
      do jm = mymBeg, mymEnd, mymSkip
        if(jm.le.jn) then
          icount = icount+1
          ila_mpiglobal = gst_getNind(jm,gstID) + jn - jm
          ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
          do levelIndex = 1, nkgdim
            zsp2(levelIndex,1,icount) = sp(ila_mpilocal,1,levelIndex)
            zsp2(levelIndex,2,icount) = sp(ila_mpilocal,2,levelIndex)
          enddo
        endif
      enddo

      if(icount.gt.0) then

        CALL DGEMM('T','N',nkgdim,2*icount,nkgdim,1.0d0,corns(1,1,jn),nkgdim,zsp2(1,1,1),nkgdim,0.0d0,zsp(1,1,1),nkgdim)
        ! CALL DGEMUL(corns(1,1,jn),nkgdim,'T',zsp2(1,1,1),nkgdim,'N',zsp(1,1,1),nkgdim,nkgdim,nkgdim,2*icount)

        icount = 0
        do jm = mymBeg, jn, mymSkip
          icount=icount+1
          ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
          ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
          do levelIndex = 1, nkgdim
            hiControlVector_out(ila_mpilocal,1,levelIndex) = zsp(levelIndex,1,icount)
            hiControlVector_out(ila_mpilocal,2,levelIndex) = zsp(levelIndex,2,icount)
          enddo
        enddo

      endif

      ! make adjustments for jm=0
      if(mymBeg == 0) then

        ila_mpiglobal = gst_getNIND(0,gstID) + jn
        ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)

       do levelIndex = 1, nkgdim
          hiControlVector_out(ila_mpilocal,1,levelIndex) = hiControlVector_out(ila_mpilocal,1,levelIndex)*sq2
          hiControlVector_out(ila_mpilocal,2,levelIndex) = hiControlVector_out(ila_mpilocal,2,levelIndex)*sq2
        enddo

      endif

    enddo
!$OMP END PARALLEL DO
    deallocate(zsp)
    deallocate(zsp2)

  end subroutine bchm_spa2gdad

  !--------------------------------------------------------------------------
  ! bchm_copyToStatevector
  !--------------------------------------------------------------------------
  subroutine bchm_copyToStatevector(statevector,gd)
    implicit none
    
    !Arguments
    type(struct_gsv) :: statevector
    real(8) :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)
    
    !Locals
    integer :: lonIndex, levelIndex, levelIndex2, latIndex, VarIndex, ilev1, ilev2
    real(4), pointer :: field_r4(:,:,:)
    real(8), pointer :: field_r8(:,:,:)

    do VarIndex = 1,numvar3d+numvar2d
      ilev1 = nsposit(VarIndex)
      ilev2 = nsposit(VarIndex+1)-1
      if (gsv_getDataKind(statevector) == 4) then
        call gsv_getField(statevector,field_r4,bchm_varNameList(VarIndex))
        do levelIndex = ilev1, ilev2
          levelIndex2 = levelIndex-ilev1+1
          do latIndex = myLatBeg, myLatEnd
            do lonIndex = myLonBeg, myLonEnd
              field_r4(lonIndex,latIndex,levelIndex2) = gd(lonIndex,latIndex,levelIndex)
            enddo
          enddo
        enddo
      else
        call gsv_getField(statevector,field_r8,bchm_varNameList(VarIndex))
        do levelIndex = ilev1, ilev2
          levelIndex2 = levelIndex-ilev1+1
          do latIndex = myLatBeg, myLatEnd
            do lonIndex = myLonBeg, myLonEnd
              field_r8(lonIndex,latIndex,levelIndex2) = gd(lonIndex,latIndex,levelIndex)
            enddo
          enddo
        enddo
      end if

    enddo
  end subroutine bchm_copyToStatevector

  !--------------------------------------------------------------------------
  ! bchm_copyFromStatevector
  !--------------------------------------------------------------------------
  subroutine bchm_copyFromStatevector(statevector,gd)
    implicit none

    !Arguments
    type(struct_gsv) :: statevector
    real(8)          :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)
    
    !Locals
    integer :: lonIndex, levelIndex, levelIndex2, latIndex, VarIndex, ilev1, ilev2
    real(4), pointer :: field_r4(:,:,:)
    real(8), pointer :: field_r8(:,:,:)

    do VarIndex = 1,numvar3d+numvar2d
      ilev1 = nsposit(VarIndex)
      ilev2 = nsposit(VarIndex+1)-1 

      if (gsv_getDataKind(statevector) == 4) then
        call gsv_getField(statevector,field_r4,bchm_varNameList(VarIndex))
        do levelIndex = ilev1, ilev2
          levelIndex2 = levelIndex-ilev1+1
          do latIndex = myLatBeg, myLatEnd
            do lonIndex = myLonBeg, myLonEnd
              gd(lonIndex,latIndex,levelIndex) = field_r4(lonIndex,latIndex,levelIndex2)
            enddo
          enddo
        enddo
      else
        call gsv_getField(statevector,field_r8,bchm_varNameList(VarIndex))
        do levelIndex = ilev1, ilev2
          levelIndex2 = levelIndex-ilev1+1
          do latIndex = myLatBeg, myLatEnd
            do lonIndex = myLonBeg, myLonEnd
              gd(lonIndex,latIndex,levelIndex) = field_r8(lonIndex,latIndex,levelIndex2)
            enddo
          enddo
        enddo
      end if
    enddo

  end subroutine bchm_copyFromStatevector

  !--------------------------------------------------------------------------
  ! bchm_finalize
  !--------------------------------------------------------------------------
  subroutine bchm_finalize()
    implicit none

    if(initialized) then
       deallocate(pressureProfile_M)
       deallocate(pressureProfile_T)
       deallocate(rgsig)
       deallocate(corns)
       deallocate(rstddev)
       deallocate(bchm_corvert,bchm_corverti,bchm_invsum)
    end if

  end subroutine bchm_finalize

  !--------------------------------------------------------------------------
  ! bchm_reduceToMPILocal
  !--------------------------------------------------------------------------
  subroutine bchm_reduceToMPILocal(cv_mpilocal,cv_mpiglobal)
    implicit none

    !Arguments
    real(8), intent(out) :: cv_mpilocal(cvDim_mpilocal)
    real(8), intent(in)  :: cv_mpiglobal(:)

    !Locals
    real(8), allocatable :: cv_allMaxMpilocal(:,:)

    integer, allocatable :: cvDim_allMpilocal(:), displs(:)
    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)

    integer :: jproc,cvDim_maxmpilocal,ierr
    integer :: levelIndex,jn,jm,ila_mpiglobal,jdim_mpilocal,jdim_mpiglobal

    if (.not.initialized) return

    call rpn_comm_allreduce(cvDim_mpilocal, cvDim_maxmpilocal, &
                            1,"MPI_INTEGER","MPI_MAX","GRID",ierr)

    if(mpi_myid == 0) then
       allocate(cvDim_allMpiLocal(mpi_nprocs))
    else
       allocate(cvDim_allMpiLocal(1))
    end if

    call rpn_comm_gather(cvDim_mpiLocal   ,1,"mpi_integer",       &
                         cvDim_allMpiLocal,1,"mpi_integer",0,"GRID",ierr)

    if(mpi_myid == 0) then
       allocate(allnBeg(mpi_nprocs))
       allocate(allnEnd(mpi_nprocs))
       allocate(allnSkip(mpi_nprocs))
       allocate(allmBeg(mpi_nprocs))
       allocate(allmEnd(mpi_nprocs))
       allocate(allmSkip(mpi_nprocs))
    else
       allocate(allnBeg(1))
       allocate(allnEnd(1))
       allocate(allnSkip(1))
       allocate(allmBeg(1))
       allocate(allmEnd(1))
       allocate(allmSkip(1))
    end if

    call rpn_comm_gather(mynBeg  ,1,"mpi_integer",       &
                         allnBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynEnd  ,1,"mpi_integer",       &
                         allnEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynSkip ,1,"mpi_integer",       &
                         allnSkip,1,"mpi_integer",0,"GRID",ierr)

    call rpn_comm_gather(mymBeg  ,1,"mpi_integer",       &
                         allmBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymEnd  ,1,"mpi_integer",       &
                         allmEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymSkip ,1,"mpi_integer",       &
                         allmSkip,1,"mpi_integer",0,"GRID",ierr)

    ! Prepare to data to be distributed
    if (mpi_myid == 0) then

       allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mpi_nprocs))

!$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,levelIndex,jm,jn,ila_mpiglobal,jdim_mpiglobal)
       do jproc = 0, (mpi_nprocs-1)
          cv_allmaxmpilocal(:,jproc+1) = 0.d0
          jdim_mpilocal = 0

          do levelIndex = 1, nkgdim
             do jm = allmBeg(jproc+1), allmEnd(jproc+1), allmSkip(jproc+1)
                do jn = allnBeg(jproc+1), allnEnd(jproc+1), allnSkip(jproc+1)

                   if(jm.le.jn) then
                      
                      ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
                      
                      ! figure out index into global control vector
                      if(jm == 0) then
                         ! for jm=0 only real part
                         jdim_mpiglobal = ila_mpiglobal
                      else
                         ! for jm>0 both real and imaginary part
                         jdim_mpiglobal = 2*ila_mpiglobal-1 - (ntrunc+1)
                      endif
                      ! add offset for level
                      jdim_mpiglobal = jdim_mpiglobal + (levelIndex-1) * (ntrunc+1)*(ntrunc+1)
                      
                      ! index into local control vector computer as in cain
                      if(jm == 0) then
                         ! only real component for jm=0
                         jdim_mpilocal = jdim_mpilocal + 1
                         cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal)
                      else
                         ! both real and imaginary components for jm>0
                         jdim_mpilocal = jdim_mpilocal + 1
                         cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal)
                         jdim_mpilocal = jdim_mpilocal + 1
                         cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal+1)
                      endif
                      
                      if (jdim_mpilocal > cvDim_allMpiLocal(jproc+1)) then
                         write(*,*)
                         write(*,*) 'ERROR: jdim_mpilocal > cvDim_allMpiLocal(jproc+1)', jdim_mpilocal, cvDim_mpilocal
                         write(*,*) '       proc, levelIndex, jn, jm = ',jproc, levelIndex, jn, jm
                         call utl_abort('bhi_reduceToMPILocal')
                      end if
                      if (jdim_mpiglobal > cvDim_mpiglobal) then
                         write(*,*)
                         write(*,*) 'ERROR: jdim_mpiglobal > cvDim_mpiglobal', jdim_mpiglobal, cvDim_mpiglobal
                         write(*,*) '       proc, levelIndex, jn, jm = ',jproc, levelIndex, jn, jm
                         call utl_abort('bhi_reduceToMPILocal')
                      end if

                   endif
                enddo
             enddo
          enddo
 
       end do ! jproc
!$OMP END PARALLEL DO

    else
       allocate(cv_allmaxmpilocal(1,1))
    end if

    !- Distribute
    allocate(displs(mpi_nprocs))
    do jproc = 0, (mpi_nprocs-1)
       displs(jproc+1) = jproc*cvDim_maxMpiLocal ! displacement wrt cv_allMaxMpiLocal from which
                                                 ! to take the outgoing data to process jproc
    end do

    call rpn_comm_scatterv(cv_allMaxMpiLocal, cvDim_allMpiLocal, displs, "mpi_double_precision", &
                           cv_mpiLocal, cvDim_mpiLocal, "mpi_double_precision", &
                           0, "GRID", ierr)

    !- End
    deallocate(displs)
    deallocate(cv_allMaxMpiLocal)
    deallocate(cvDim_allMpiLocal)
    deallocate(allnBeg)
    deallocate(allnEnd)
    deallocate(allnSkip)
    deallocate(allmBeg)
    deallocate(allmEnd)
    deallocate(allmSkip)

  end subroutine bchm_reduceToMPILocal

  !--------------------------------------------------------------------------
  ! bchm_reduceToMPILocal_r4
  !--------------------------------------------------------------------------
  subroutine bchm_reduceToMPILocal_r4(cv_mpilocal,cv_mpiglobal)
    implicit none

    !Arguments
    real(4), intent(out) :: cv_mpilocal(cvDim_mpilocal)
    real(4), intent(in)  :: cv_mpiglobal(:)

    !Locals
    real(4), allocatable :: cv_allMaxMpilocal(:,:)

    integer, allocatable :: cvDim_allMpilocal(:), displs(:)
    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)

    integer :: jproc,cvDim_maxmpilocal,ierr
    integer :: levelIndex,jn,jm,ila_mpiglobal,jdim_mpilocal,jdim_mpiglobal

    if (.not.initialized) return

    call rpn_comm_allreduce(cvDim_mpilocal, cvDim_maxmpilocal, &
                            1,"MPI_INTEGER","MPI_MAX","GRID",ierr)

    if(mpi_myid == 0) then
       allocate(cvDim_allMpiLocal(mpi_nprocs))
    else
       allocate(cvDim_allMpiLocal(1))
    end if

    call rpn_comm_gather(cvDim_mpiLocal   ,1,"mpi_integer",       &
                         cvDim_allMpiLocal,1,"mpi_integer",0,"GRID",ierr)

    if(mpi_myid == 0) then
       allocate(allnBeg(mpi_nprocs))
       allocate(allnEnd(mpi_nprocs))
       allocate(allnSkip(mpi_nprocs))
       allocate(allmBeg(mpi_nprocs))
       allocate(allmEnd(mpi_nprocs))
       allocate(allmSkip(mpi_nprocs))
    else
       allocate(allnBeg(1))
       allocate(allnEnd(1))
       allocate(allnSkip(1))
       allocate(allmBeg(1))
       allocate(allmEnd(1))
       allocate(allmSkip(1))
    end if

    call rpn_comm_gather(mynBeg  ,1,"mpi_integer",       &
                         allnBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynEnd  ,1,"mpi_integer",       &
                         allnEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynSkip ,1,"mpi_integer",       &
                         allnSkip,1,"mpi_integer",0,"GRID",ierr)

    call rpn_comm_gather(mymBeg  ,1,"mpi_integer",       &
                         allmBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymEnd  ,1,"mpi_integer",       &
                         allmEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymSkip ,1,"mpi_integer",       &
                         allmSkip,1,"mpi_integer",0,"GRID",ierr)

    ! Prepare to data to be distributed
    if (mpi_myid == 0) then

       allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mpi_nprocs))

!$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,levelIndex,jm,jn,ila_mpiglobal,jdim_mpiglobal)
       do jproc = 0, (mpi_nprocs-1)
          cv_allmaxmpilocal(:,jproc+1) = 0.d0
          jdim_mpilocal = 0

          do levelIndex = 1, nkgdim
             do jm = allmBeg(jproc+1), allmEnd(jproc+1), allmSkip(jproc+1)
                do jn = allnBeg(jproc+1), allnEnd(jproc+1), allnSkip(jproc+1)

                   if(jm.le.jn) then
                      
                      ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
                      
                      ! figure out index into global control vector
                      if(jm == 0) then
                         ! for jm=0 only real part
                         jdim_mpiglobal = ila_mpiglobal
                      else
                         ! for jm>0 both real and imaginary part
                         jdim_mpiglobal = 2*ila_mpiglobal-1 - (ntrunc+1)
                      endif
                      ! add offset for level
                      jdim_mpiglobal = jdim_mpiglobal + (levelIndex-1) * (ntrunc+1)*(ntrunc+1)
                      
                      ! index into local control vector computer as in cain
                      if(jm == 0) then
                         ! only real component for jm=0
                         jdim_mpilocal = jdim_mpilocal + 1
                         cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal)
                      else
                         ! both real and imaginary components for jm>0
                         jdim_mpilocal = jdim_mpilocal + 1
                         cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal)
                         jdim_mpilocal = jdim_mpilocal + 1
                         cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal+1)
                      endif
                      
                      if (jdim_mpilocal > cvDim_allMpiLocal(jproc+1)) then
                         write(*,*)
                         write(*,*) 'ERROR: jdim_mpilocal > cvDim_allMpiLocal(jproc+1)', jdim_mpilocal, cvDim_mpilocal
                         write(*,*) '       proc, levelIndex, jn, jm = ',jproc, levelIndex, jn, jm
                         call utl_abort('bhi_reduceToMPILocal')
                      end if
                      if (jdim_mpiglobal > cvDim_mpiglobal) then
                         write(*,*)
                         write(*,*) 'ERROR: jdim_mpiglobal > cvDim_mpiglobal', jdim_mpiglobal, cvDim_mpiglobal
                         write(*,*) '       proc, levelIndex, jn, jm = ',jproc, levelIndex, jn, jm
                         call utl_abort('bhi_reduceToMPILocal')
                      end if

                   endif
                enddo
             enddo
          enddo
 
       end do ! jproc
!$OMP END PARALLEL DO

    else
       allocate(cv_allmaxmpilocal(1,1))
    end if

    !- Distribute
    allocate(displs(mpi_nprocs))
    do jproc = 0, (mpi_nprocs-1)
       displs(jproc+1) = jproc*cvDim_maxMpiLocal ! displacement wrt cv_allMaxMpiLocal from which
                                                 ! to take the outgoing data to process jproc
    end do

    call rpn_comm_scatterv(cv_allMaxMpiLocal, cvDim_allMpiLocal, displs, "mpi_real4", &
                           cv_mpiLocal, cvDim_mpiLocal, "mpi_real4", &
                           0, "GRID", ierr)

    !- End
    deallocate(displs)
    deallocate(cv_allMaxMpiLocal)
    deallocate(cvDim_allMpiLocal)
    deallocate(allnBeg)
    deallocate(allnEnd)
    deallocate(allnSkip)
    deallocate(allmBeg)
    deallocate(allmEnd)
    deallocate(allmSkip)

  end subroutine bchm_reduceToMPILocal_r4

  !--------------------------------------------------------------------------
  ! bchm_expandToMPIGlobal
  !--------------------------------------------------------------------------
  subroutine bchm_expandToMPIGlobal(cv_mpilocal,cv_mpiglobal)
    implicit none

    !Arguments
    real(8), intent(in)  :: cv_mpilocal(cvDim_mpilocal)
    real(8), intent(out) :: cv_mpiglobal(:)

    !Locals
    real(8), allocatable :: cv_maxmpilocal(:)
    real(8), pointer :: cv_allmaxmpilocal(:,:) => null()
    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)
    integer :: levelIndex, jn, jm, jproc, ila_mpiglobal, jdim_mpilocal, jdim_mpiglobal, ierr, cvDim_maxmpilocal

    if (.not.initialized) return

    !
    !- 1.  Gather all local control vectors onto mpi task 0
    !
    call rpn_comm_allreduce(cvDim_mpilocal,cvDim_maxmpilocal,1,"mpi_integer","mpi_max","GRID",ierr)

    allocate(cv_maxmpilocal(cvDim_maxmpilocal))

    nullify(cv_allmaxmpilocal)
    if(mpi_myid == 0) then
       allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mpi_nprocs))
    else
       allocate(cv_allmaxmpilocal(1,1))
    end if

    cv_maxmpilocal(:) = 0.0d0
    cv_maxmpilocal(1:cvDim_mpilocal) = cv_mpilocal(1:cvDim_mpilocal)

    call rpn_comm_gather(cv_maxmpilocal,    cvDim_maxmpilocal, "mpi_double_precision",  &
                         cv_allmaxmpilocal, cvDim_maxmpilocal, "mpi_double_precision", 0, "GRID", ierr )

    deallocate(cv_maxmpilocal)

    !
    !- 2.  Reorganize gathered mpilocal control vectors into the mpiglobal control vector
    !
    if(mpi_myid == 0) then
       allocate(allnBeg(mpi_nprocs))
       allocate(allnEnd(mpi_nprocs))
       allocate(allnSkip(mpi_nprocs))
       allocate(allmBeg(mpi_nprocs))
       allocate(allmEnd(mpi_nprocs))
       allocate(allmSkip(mpi_nprocs))
    else
       allocate(allnBeg(1))
       allocate(allnEnd(1))
       allocate(allnSkip(1))
       allocate(allmBeg(1))
       allocate(allmEnd(1))
       allocate(allmSkip(1))
    end if

    call rpn_comm_gather(mynBeg  ,1,"mpi_integer",       &
                         allnBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynEnd  ,1,"mpi_integer",       &
                         allnEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynSkip ,1,"mpi_integer",       &
                         allnSkip,1,"mpi_integer",0,"GRID",ierr)

    call rpn_comm_gather(mymBeg  ,1,"mpi_integer",       &
                         allmBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymEnd  ,1,"mpi_integer",       &
                         allmEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymSkip ,1,"mpi_integer",       &
                         allmSkip,1,"mpi_integer",0,"GRID",ierr)

    if(mpi_myid == 0) then
      cv_mpiglobal(:) = 0.0d0

!$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,levelIndex,jm,jn,ila_mpiglobal,jdim_mpiglobal)
      do jproc = 0, (mpi_nprocs-1)
        jdim_mpilocal = 0

        do levelIndex = 1, nkgdim
          do jm = allmBeg(jproc+1), allmEnd(jproc+1), allmSkip(jproc+1)
            do jn = allnBeg(jproc+1), allnEnd(jproc+1), allnSkip(jproc+1)
              if(jm.le.jn) then

                ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm

                ! figure out index into global control vector
                if(jm == 0) then
                  ! for jm=0 only real part
                  jdim_mpiglobal = ila_mpiglobal
                else
                  ! for jm>0 both real and imaginary part
                  jdim_mpiglobal = 2*ila_mpiglobal-1 - (ntrunc+1)
                endif
                ! add offset for level
                jdim_mpiglobal = jdim_mpiglobal + (levelIndex-1) * (ntrunc+1)*(ntrunc+1)

                ! index into local control vector
                if(jm == 0) then
                  ! only real component for jm=0
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_mpiglobal(jdim_mpiglobal) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                else
                  ! both real and imaginary components for jm>0
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_mpiglobal(jdim_mpiglobal) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_mpiglobal(jdim_mpiglobal+1) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                endif

                if(jdim_mpiglobal.gt.cvDim_mpiglobal)   &
                  write(*,*) 'ERROR: jdim,cvDim,mpiglobal=',jdim_mpiglobal,cvDim_mpiglobal,levelIndex,jn,jm

              endif
            enddo
          enddo
        enddo
      enddo ! jproc
!$OMP END PARALLEL DO

    endif ! myid == 0 

    deallocate(allnBeg)
    deallocate(allnEnd)
    deallocate(allnSkip)
    deallocate(allmBeg)
    deallocate(allmEnd)
    deallocate(allmSkip)
    deallocate(cv_allmaxmpilocal)

  end subroutine bchm_expandToMPIGlobal

  !--------------------------------------------------------------------------
  ! bchm_expandToMPIGlobal_r4
  !--------------------------------------------------------------------------
  subroutine bchm_expandToMPIGlobal_r4(cv_mpilocal,cv_mpiglobal)
    implicit none

    !Arguments
    real(4), intent(in)  :: cv_mpilocal(cvDim_mpilocal)
    real(4), intent(out) :: cv_mpiglobal(:)

    !Locals
    real(4), allocatable :: cv_maxmpilocal(:)
    real(4), pointer :: cv_allmaxmpilocal(:,:) => null()
    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)
    integer :: levelIndex, jn, jm, jproc, ila_mpiglobal, jdim_mpilocal, jdim_mpiglobal, ierr, cvDim_maxmpilocal

    if (.not.initialized) return

    !
    !- 1.  Gather all local control vectors onto mpi task 0
    !
    call rpn_comm_allreduce(cvDim_mpilocal,cvDim_maxmpilocal,1,"mpi_integer","mpi_max","GRID",ierr)

    allocate(cv_maxmpilocal(cvDim_maxmpilocal))

    nullify(cv_allmaxmpilocal)
    if(mpi_myid == 0) then
       allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mpi_nprocs))
    else
       allocate(cv_allmaxmpilocal(1,1))
    end if

    cv_maxmpilocal(:) = 0.0d0
    cv_maxmpilocal(1:cvDim_mpilocal) = cv_mpilocal(1:cvDim_mpilocal)

    call rpn_comm_gather(cv_maxmpilocal,    cvDim_maxmpilocal, "mpi_real4",  &
                         cv_allmaxmpilocal, cvDim_maxmpilocal, "mpi_real4", 0, "GRID", ierr )

    deallocate(cv_maxmpilocal)

    !
    !- 2.  Reorganize gathered mpilocal control vectors into the mpiglobal control vector
    !
    if(mpi_myid == 0) then
       allocate(allnBeg(mpi_nprocs))
       allocate(allnEnd(mpi_nprocs))
       allocate(allnSkip(mpi_nprocs))
       allocate(allmBeg(mpi_nprocs))
       allocate(allmEnd(mpi_nprocs))
       allocate(allmSkip(mpi_nprocs))
    else
       allocate(allnBeg(1))
       allocate(allnEnd(1))
       allocate(allnSkip(1))
       allocate(allmBeg(1))
       allocate(allmEnd(1))
       allocate(allmSkip(1))
    end if

    call rpn_comm_gather(mynBeg  ,1,"mpi_integer",       &
                         allnBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynEnd  ,1,"mpi_integer",       &
                         allnEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynSkip ,1,"mpi_integer",       &
                         allnSkip,1,"mpi_integer",0,"GRID",ierr)

    call rpn_comm_gather(mymBeg  ,1,"mpi_integer",       &
                         allmBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymEnd  ,1,"mpi_integer",       &
                         allmEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymSkip ,1,"mpi_integer",       &
                         allmSkip,1,"mpi_integer",0,"GRID",ierr)

    if(mpi_myid == 0) then
      cv_mpiglobal(:) = 0.0d0

!$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,levelIndex,jm,jn,ila_mpiglobal,jdim_mpiglobal)
      do jproc = 0, (mpi_nprocs-1)
        jdim_mpilocal = 0

        do levelIndex = 1, nkgdim
          do jm = allmBeg(jproc+1), allmEnd(jproc+1), allmSkip(jproc+1)
            do jn = allnBeg(jproc+1), allnEnd(jproc+1), allnSkip(jproc+1)
              if(jm.le.jn) then

                ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm

                ! figure out index into global control vector
                if(jm == 0) then
                  ! for jm=0 only real part
                  jdim_mpiglobal = ila_mpiglobal
                else
                  ! for jm>0 both real and imaginary part
                  jdim_mpiglobal = 2*ila_mpiglobal-1 - (ntrunc+1)
                endif
                ! add offset for level
                jdim_mpiglobal = jdim_mpiglobal + (levelIndex-1) * (ntrunc+1)*(ntrunc+1)

                ! index into local control vector
                if(jm == 0) then
                  ! only real component for jm=0
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_mpiglobal(jdim_mpiglobal) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                else
                  ! both real and imaginary components for jm>0
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_mpiglobal(jdim_mpiglobal) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_mpiglobal(jdim_mpiglobal+1) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                endif

                if(jdim_mpiglobal.gt.cvDim_mpiglobal)   &
                  write(*,*) 'ERROR: jdim,cvDim,mpiglobal=',jdim_mpiglobal,cvDim_mpiglobal,levelIndex,jn,jm

              endif
            enddo
          enddo
        enddo
      enddo ! jproc
!$OMP END PARALLEL DO

    endif ! myid == 0 

    deallocate(allnBeg)
    deallocate(allnEnd)
    deallocate(allnSkip)
    deallocate(allmBeg)
    deallocate(allmEnd)
    deallocate(allmSkip)
    deallocate(cv_allmaxmpilocal)

  end subroutine bchm_expandtompiglobal_r4

  !--------------------------------------------------------------------------
  ! bchm_corvert_mult
  !--------------------------------------------------------------------------
  subroutine bchm_corvert_mult(varName,rmat_in,rmat_out,lvl_top,lvl_bot,ndim1, &
                               ndim2,ndim3,lrgsig,itype,rsig_opt)
    !
    !:Purpose: To multiply with local matrix C=bchm_corvert (itype>0) or
    !          CI=bchm_corverti (itype<0).
    !
    !            Given that A=rmat_in (=input):
    !
    !              =====  ======
    !              itype  output
    !              =====  ======
    !                0     D(i,j)=A(i,j)/sum(C(1:n,i))
    !                1     A*C
    !                2     C*A
    !                3     A*C*A^T
    !               -1     A*CI
    !               -2     CI*A
    !               -3     A*CI*A^T
    !              =====  ======
    !
    !:Comments:
    !
    !  -   If rmat_in is a 1-D vector, then
    !
    !      - for cases +/- 2, one should have set ndim2=1 and
    !        ndim1=vector-length.
    !      - for cases +/- 1,3, one should have set ndim1=1 and
    !        ndim2=vector-length.
    !
    !  -   rmat_out assumed to be initialized prior to call to bchm_corvert_mult.
    !
    !  -   bchm_corverti is the inverse of bchm_corvert.
    !
    !:Arguments:
    !    :ndim3: Output dimension expected by the calling routine.
    !            bchm_corvert_mult insists on these values for ndim3:
    !
    !            - for itype  = +/-3, ndim3=ndim1
    !            - for itype != +/-3, ndim3=ndim2
    !
    implicit none

    ! Arguments
    character(len=*), intent(in) :: varName ! Variable name
    real(8), intent(in) :: rmat_in(ndim1,ndim2) ! Input matrix/vector A (see comments section)
    real(8), intent(inout) :: rmat_out(ndim1,ndim3) ! Output matrix/vector 
    integer, intent(in) :: lvl_top(ndim1) ! Top level of non-zero values in rmat_in
    integer, intent(in) :: lvl_bot(ndim1) ! Bottom level of non-zero values in rmat_in
    integer, intent(in) :: ndim1 ! Matrix dimension - to be one 1D input vectors
    integer, intent(in) :: ndim2 ! Matrix dimension
    integer, intent(in) :: ndim3 ! Matrix dimension
    logical, intent(in) :: lrgsig! Indicates if rsig_opt to be included as part of C or CI below.
    integer, intent(in) :: itype ! Type of operator (see above).
    real(8), intent(in), optional :: rsig_opt(ndim2) ! Input error standard deviations. Required when lrgsig=.true.
   
    !Locals
    integer :: VarIndex,jk1,jk2,jk3,nsize
    real(8) :: rmat_work(ndim2,ndim2)

    if (.not.initialized) return
 
    if (.not.present(rsig_opt).and.lrgsig) call utl_abort('bchm_corvert_mult: Missing rsig_opt')  

				! Determine location and size in bchm_corvert/bchm_corverti
    
    do VarIndex = 1, numvar3d+numvar2d
      if (trim(varName) == trim(bchm_varNameList(VarIndex))) exit
    end do
    if  (VarIndex > numvar3d+numvar2d) call utl_abort('bchm_corvert_mult: Variable not found')
    
    nsize = nsposit(VarIndex+1)-nsposit(VarIndex)
    
			! Apply matrix/vector multiplication.
    
    rmat_work(:,:) = 0.0D0

    if (itype == 3) then

      !  A*C*A^T
       
      if (ndim3 /= ndim1 .or. ndim2 /= nsize) call utl_abort('bchm_corvert_mult: Size does not match - itype=3')
      if (lrgsig) then
        do jk1 = 1, nsize
        do jk2 = 1, ndim1
          rmat_work(jk1,jk2) = sum(rmat_in(jk2,lvl_top(jk2):lvl_bot(jk2))*bchm_corvert(jk1,lvl_top(jk2):lvl_bot(jk2),VarIndex) &
                                *rsig_opt(lvl_top(jk2):lvl_bot(jk2)))*rsig_opt(jk1)
        end do
        end do
      else
        do jk1 = 1, nsize
        do jk2 = 1, ndim1
          rmat_work(jk1,jk2) = sum(rmat_in(jk2,lvl_top(jk2):lvl_bot(jk2))*bchm_corvert(jk1,lvl_top(jk2):lvl_bot(jk2),VarIndex))
        end do
        end do
      end if
      do jk1=1,ndim1
      do jk2=1,ndim1
        rmat_out(jk1,jk2) = sum(rmat_in(jk1,lvl_top(jk1):lvl_bot(jk1))*rmat_work(lvl_top(jk1):lvl_bot(jk1),jk2))
      end do
      end do
 
    else if (itype == -3) then

      ! A*CI*A^T
       
      if (ndim3 /= ndim1 .or. ndim2 /= nsize) call utl_abort('bchm_corvert_mult: Size does not match - itype=3')
      if (lrgsig) then
        do jk1 = 1, nsize
        do jk2 = 1, ndim1
          rmat_work(jk1,jk2) = sum(rmat_in(jk2,lvl_top(jk2):lvl_bot(jk2))*bchm_corverti(jk1,lvl_top(jk2):lvl_bot(jk2),VarIndex) &
                                /rsig_opt(lvl_top(jk2):lvl_bot(jk2)))/rsig_opt(jk1)
        end do
        end do
      else
        do jk1 = 1, nsize
        do jk2 = 1, ndim1
          rmat_work(jk1,jk2) = sum(rmat_in(jk2,lvl_top(jk2):lvl_bot(jk2))*bchm_corverti(jk1,lvl_top(jk2):lvl_bot(jk2),VarIndex))
        end do
        end do
      endif
       
      do jk1 = 1, ndim1
      do jk2 = 1, ndim1
        rmat_out(jk1,jk2) = sum(rmat_in(jk1,lvl_top(jk1):lvl_bot(jk1))*rmat_work(lvl_top(jk1):lvl_bot(jk1),jk2))
      end do
      end do
      
    else if (itype == 2) then
    
      ! C*A

      if (ndim3 /= ndim2 .or. ndim1 /= nsize) call utl_abort('bchm_corvert_mult: Size does not match - itype=2')
      if (lrgsig) then
        do jk1 = 1, nsize
        do jk3 = 1, nsize
          do jk2 = lvl_top(jk3), lvl_bot(jk3)     ! Instead of do jk2=1,ndim2
            rmat_out(jk1,jk2) = rmat_out(jk1,jk2)+bchm_corvert(jk1,jk3,VarIndex)*rsig_opt(jk3)*rsig_opt(jk1)*rmat_in(jk3,jk2)
          end do
        end do
        end do
      else
        do jk1 = 1, nsize
        do jk3 = 1, nsize
          do jk2 = lvl_top(jk3), lvl_bot(jk3)     ! Instead of do jk2=1,ndim2
            rmat_out(jk1,jk2) = rmat_out(jk1,jk2)+bchm_corvert(jk1,jk3,VarIndex)*rmat_in(jk3,jk2)
          end do
        end do
        end do
      endif

    else if (itype == -2) then
    
      ! CI*A

      if (ndim3 /= ndim2 .or. ndim1 /= nsize) call utl_abort('bchm_corvert_mult: Size does match - itype=-2')
      if (lrgsig) then
        do jk1 = 1, nsize
        do jk3 = 1, nsize
          do jk2 = lvl_top(jk3), lvl_bot(jk3)     ! Instead of do jk2=1,ndim2
            rmat_out(jk1,jk2) = rmat_out(jk1,jk2)+bchm_corverti(jk1,jk3,VarIndex)*rmat_in(jk3,jk2)/(rsig_opt(jk3)*rsig_opt(jk1))
          end do
        end do
        end do
      else
        do jk1 = 1, nsize
        do jk3 = 1, nsize
          do jk2 = lvl_top(jk3), lvl_bot(jk3)     ! Instead of do jk2=1,ndim2
            rmat_out(jk1,jk2) = rmat_out(jk1,jk2)+bchm_corverti(jk1,jk3,VarIndex)*rmat_in(jk3,jk2)
          end do
        end do
        end do
      endif
          
    else if (itype == 1) then
    
      ! A*C

      if (ndim3 /= ndim2 .or. ndim2 /= nsize) call utl_abort('bchm_corvert_mult: Size does not match - itype=1')
      if (lrgsig) then
        do jk1 = 1, ndim1
        do jk2 = 1, nsize
          rmat_out(jk1,jk2) = sum(rmat_in(jk1,lvl_top(jk1):lvl_bot(jk1))*bchm_corvert(lvl_top(jk1):lvl_bot(jk1),jk2,VarIndex) &
                               *rsig_opt(lvl_top(jk1):lvl_bot(jk1)))*rsig_opt(jk2)
        end do
        end do
      else
        do jk1 = 1, ndim1
        do jk2 = 1, nsize
          rmat_out(jk1,jk2) = sum(rmat_in(jk1,lvl_top(jk1):lvl_bot(jk1))*bchm_corvert(lvl_top(jk1):lvl_bot(jk1),jk2,VarIndex))
        end do
        end do
      endif
          
    else if (itype == -1) then

      ! A*CI

      if (ndim3 /= ndim2 .or. ndim2 /= nsize) call utl_abort('bchm_corvert_mult: Size does not match - itype=-1')
      if (lrgsig) then
        do jk1 = 1, ndim1
        do jk2 = 1, nsize
          rmat_out(jk1,jk2) = sum(rmat_in(jk1,lvl_top(jk1):lvl_bot(jk1))*bchm_corverti(lvl_top(jk1):lvl_bot(jk1),jk2,VarIndex) &
                                 /rsig_opt(lvl_top(jk1):lvl_bot(jk1)))/rsig_opt(jk2)
        end do
        end do
      else
        do jk1 = 1, ndim1
        do jk2 = 1, nsize
          rmat_out(jk1,jk2) = sum(rmat_in(jk1,lvl_top(jk1):lvl_bot(jk1))*bchm_corverti(lvl_top(jk1):lvl_bot(jk1),jk2,VarIndex))
        end do
        end do
      endif

    else if (itype == 0) then

      ! Instead of A*CI do D(i,j)=A(i,j)/sum(C(1:n,i))

      if (ndim3 /= ndim2 .or. ndim2 /= nsize) call utl_abort('bchm_corvert_mult: Size does not match - itype=0')
      if (lrgsig) then
        do jk1 = 1, ndim1
        do jk2 = lvl_top(jk1), lvl_bot(jk1)   ! instead of do jk2=1,nsize
          rmat_out(jk1,jk2) = rmat_in(jk1,jk2)/rsig_opt(jk2) &
              /sum(bchm_corvert(1:nsize,jk2,VarIndex)/rsig_opt(1:nsize)) 
        end do
        end do
      else
        do jk1 = 1, ndim1
          rmat_out(jk1,lvl_top(jk1):lvl_bot(jk1)) = rmat_in(jk1,lvl_top(jk1):lvl_bot(jk1))*bchm_invsum(lvl_top(jk1):lvl_bot(jk1),VarIndex)            
!          do jk2 = lvl_top(jk1), lvl_bot(jk1)   ! instead of do jk2=1,nsize
!              rmat_out(jk1,jk2) = rmat_in(jk1,jk2)/sum(bchm_corvert(1:nsize,jk2,VarIndex))
!          end do
        end do
      endif

    else
      call utl_abort('bchm_corvert_mult: Requested type not found')
    end if
   
  end subroutine bchm_corvert_mult

  !--------------------------------------------------------------------------
  ! bchm_getsigma
  !--------------------------------------------------------------------------
  subroutine bchm_getsigma(varName,ndim2,xlat,xlong,rsig,vlev_opt)
    !
    !:Purpose: To interpolate error std. dev. to obs location.
    !
    !:Comment: Does not currently account for transform_varKindCH
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: varName ! Variable name
    integer, intent(in) :: ndim2    ! Max array size
    real(8), intent(in) :: xlat     ! Target latitude
    real(8), intent(in) :: xlong    ! Target longitude
    real(8), intent(out) :: rsig(:) ! Error std. dev.
    real(8), intent(in), optional :: vlev_opt(:) ! Target vertical levels
    
    ! Locals:
    integer :: VarIndex,levelIndex1,levelIndex2,latIndex,lonIndexg,j,nsize
    real(8) :: rvar(ndim2,2),zc1,zc2,rlat1,rlat2,rlong1,rlong2,zd1,zd2,rsig_max
    integer :: ilev1,ilev2
    real(8) :: dz
      
    integer, parameter :: itype=0

    if (.not.initialized) return
    
    ! Determine location and size of background error std. dev.

    do VarIndex = 1, numvar3d+numvar2d
      if (trim(varName) == trim(bchm_varNameList(VarIndex))) exit
    end do
    if  (VarIndex > numvar3d+numvar2d) call utl_abort('bchm_getsigma: Variable not found')
    
    levelIndex1 = nsposit(VarIndex)
    levelIndex2 = nsposit(VarIndex+1)-1

    nsize = levelIndex2-levelIndex1+1
    
    if (.not.present(vlev_opt) .and. nsize /= ndim2 ) then
      write(6,*) 'NSIZE, NDIM2: ',nsize,ndim2
      call utl_abort('bchm_getsigma: Inconsistent size')
    end if
 
    ! Determine reference longitude index 

    lonIndexg = 2    
    do while (xlong > rlong(lonIndexg) .and. lonIndexg < ni_l+1) 
      lonIndexg = lonIndexg+1
    end do

    ! Set interpolation weights

    rlong2 = rlong(lonIndexg)
    rlong1 = rlong(lonIndexg-1)
    
    zd2 = (xlong-rlong1)/(rlong2-rlong1)
    zd1 = 1.0-zd2
     
    ! Determine reference latitude index (could alternatively use gst_getrlati)

    latIndex = 2 
    do while (xlat > rlat(latIndex) .and. latIndex < nj_l) 
      ! Note: gst_getrlati needs to be consistent with rlat (except start of indexing and extrapolated lats)
      latIndex = latIndex+1
    enddo

    ! Set interpolation weights

    rlat2 = rlat(latIndex)
    rlat1 = rlat(latIndex-1)
    
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

      rvar(1:nsize,1) = zc1*rgsig(lonIndexg-1,latIndex-1,levelIndex1:levelIndex2)**2 + zc2*rgsig(lonIndexg-1,latIndex,levelIndex1:levelIndex2)**2
    
      rvar(1:nsize,2) = zc1*rgsig(lonIndexg,latIndex-1,levelIndex1:levelIndex2)**2 + zc2*rgsig(lonIndexg,latIndex,levelIndex1:levelIndex2)**2
       
      rvar(1:nsize,1) = zd1*rvar(1:nsize,1) + zd2*rvar(1:nsize,2)
      
      if (nsize /= ndim2) then
         do ilev1=1, ndim2
            if (pressureProfile_T(1).ge.vlev_opt(ilev1)) then
              rsig(ilev1)=sqrt(rvar(1,1))
            else if (pressureProfile_T(nsize).le.vlev_opt(ilev1)) then
              rsig(ilev1)=sqrt(rvar(nsize,1))
            else  
              do ilev2=1,nsize-1
                 if (vlev_opt(ilev1) >= pressureProfile_T(ilev2) .and. &
                     vlev_opt(ilev1) < pressureProfile_T(ilev2+1)) exit
              end do
              dz=log(vlev_opt(ilev1)/pressureProfile_T(ilev2)) &
                 /log(pressureProfile_T(ilev2+1)/pressureProfile_T(ilev2))
              rsig(ilev1)=sqrt(rvar(ilev2,1)*(1.0-dz)+rvar(ilev2+1,1)*dz)
            end if
        end do
      else
         rsig(1:nsize) = sqrt(rvar(1:nsize,1))
      end if 
    
    else 

      ! Interpolation of std. dev. (to reduce execution time)

      rvar(1:nsize,1) = zc1*rgsig(lonIndexg-1,latIndex-1,levelIndex1:levelIndex2) + zc2*rgsig(lonIndexg-1,latIndex,levelIndex1:levelIndex2)     

      rvar(1:nsize,2) = zc1*rgsig(lonIndexg,latIndex-1,levelIndex1:levelIndex2) + zc2*rgsig(lonIndexg,latIndex,levelIndex1:levelIndex2)

      rvar(1:nsize,1) = zd1*rvar(1:nsize,1) + zd2*rvar(1:nsize,2)
    
      if (nsize /= ndim2) then
         do ilev1=1, ndim2
            if (pressureProfile_T(1).ge.vlev_opt(ilev1)) then
              rsig(ilev1)=rvar(1,1)
            else if (pressureProfile_T(nsize).le.vlev_opt(ilev1)) then
              rsig(ilev1)=rvar(nsize,1)
            else  
              do ilev2=1,nsize-1
                 if (vlev_opt(ilev1) >= pressureProfile_T(ilev2) .and. &
                     vlev_opt(ilev1) < pressureProfile_T(ilev2+1)) exit
              end do
              dz=log(vlev_opt(ilev1)/pressureProfile_T(ilev2)) &
                 /log(pressureProfile_T(ilev2+1)/pressureProfile_T(ilev2))
              rsig(ilev1)=rvar(ilev2,1)*(1.0-dz)+rvar(ilev2+1,1)*dz
            end if
         end do
      else
         rsig(1:nsize) = rvar(1:nsize,1)
      end if 
      
    end if
    
    rsig_max = maxval(rgsig(lonIndexg-1:lonIndexg,latIndex-1:latIndex,:))
    do j = 1, nsize
      if (rsig(j) < 0.0 .or. rsig(j) > 1.1*rsig_max) then
        write(*,*) 'bchm_getsigma: Interpolated sigma incorrect:'
        write(*,*) 'bchm_getsigma:   zc1,zc2,zd1,zd2 = ',zc1,zc2,zd1,zd2
        write(*,*) 'bchm_getsigma:   rsig,rsig_max = ',rsig(j),rsig_max
        write(*,*) 'bchm_getsigma:   lonIndexg,latIndex,j,xlong,xlat = ',lonIndexg,latIndex,j,xlong,xlat
        write(*,*) 'bchm_getsigma:   rlong2,rlong1,rlat1,rlat2 = ',rlong2,rlong1,rlat1,rlat2
        call utl_abort('bchm_getsigma')
      end if
    end do

  end subroutine bchm_getsigma

  !--------------------------------------------------------------------------
  ! bchm_getbgsigma
  !--------------------------------------------------------------------------
  real(8) function bchm_getbgsigma(lonIndex,latIndex,levelIndex,VarIndex)
    !
    !:Purpose: To get error std. dev. a specified grid point and for specified
    !          field
    !
    implicit none

    !Arguments
    integer :: lonIndex, latIndex, VarIndex, levelIndex

    bchm_getbgsigma = rgsig(lonIndex, latIndex, nsposit(VarIndex)-1+levelIndex)

  end function bchm_getbgsigma
  
  !--------------------------------------------------------------------------
  ! bchm_resetCorvert
  !--------------------------------------------------------------------------
  subroutine bchm_resetCorvert(nvlev)
    !
    !:Purpose: Vertically interpolate error correlation matrix fields to generate 
    !          approximate matrices/vectors in trial field vertical levels via
    !          interpolation. No need to make matrix positive definite for this approximation.
    !
    implicit none
    
    !Arguments
    integer :: nvlev
    
    !Locals
    type(struct_vco),pointer :: vco_trl => null()
    real(8), pointer :: vlev(:)
    integer :: status
       
    integer :: nsize,ilev1,ilev2,j,d1,d2
    real(8) :: dz

    real(8), allocatable :: wtemp1(:,:), wtemp2(:,:,:)
    
    if (.not.initialized .or. numvar3d == 0 .or. lresetCorvert) return

    ! Set reference pressure levels of trial for a surface pressure of zps
    call vco_SetupFromFile( vco_trl, './trlm_01' )
    status = vgd_levels( vco_trl%vgrid, ip1_list=vco_trl%ip1_T, levels=vlev, &
                         sfc_field=zps, in_log=.false.)
    nvlev = vco_trl%nlev_T
    
    !nsize = nsposit(2)-nsposit(1)+1
    nsize=nlev_T
    if (nsize == nvlev ) then
       lresetCorvert=.true.
       return
    end if
        
    d2=0 
    d1=nvlev-nsize
    if (d1 < 0) then
       d1=0
       d2=d1
    end if
    
    allocate(wtemp2(nvlev, nvlev,numvar3d))
    allocate(wtemp1(nvlev,numvar3d))
    wtemp2(:,:,:)=0.0d0
    
    ! Apply interpolation

    do ilev1 = 1, nvlev
       if (pressureProfile_T(1).ge.vlev(ilev1)) then
          wtemp1(ilev1,1:numvar3d)=bchm_invsum(1,1:numvar3d)

          wtemp2(ilev1,1:nsize+d2,1:numvar3d)=bchm_corvert(1,1:nsize+d2,1:numvar3d)
          wtemp2(1:nsize+d2,ilev1,1:numvar3d)=bchm_corvert(1:nsize+d2,1,1:numvar3d)          
       else if (pressureProfile_T(nsize).le.vlev(ilev1)) then
          wtemp1(ilev1,1:numvar3d)=bchm_invsum(nsize,1:numvar3d)
          
          wtemp2(ilev1,1+ilev1-nsize-d2:ilev1,1:numvar3d)=bchm_corvert(nsize,1-d2:nsize,1:numvar3d)
          wtemp2(1+ilev1-nsize-d2:ilev1,ilev1,1:numvar3d)=bchm_corvert(1-d2:nsize,nsize,1:numvar3d)          
       else  
          do ilev2=1,nsize-1
             if (vlev(ilev1) >= pressureProfile_T(ilev2) .and. &
                 vlev(ilev1) < pressureProfile_T(ilev2+1)) exit
          end do

          dz=log(vlev(ilev1)/pressureProfile_T(ilev2)) &
              /log(pressureProfile_T(ilev2+1)/pressureProfile_T(ilev2))     
          wtemp1(ilev1,1:numvar3d)=bchm_invsum(ilev2,1:numvar3d)*(1.0-dz) &
                                  +bchm_invsum(ilev2+1,1:numvar3d)*dz
                    
          do j=1-max(ilev1-d1,1),nvlev-min(d1+ilev1,nvlev)
            wtemp2(ilev1,ilev1+j,1:numvar3d)=(bchm_corvert(ilev2,ilev2+j,1:numvar3d) &
                      *(1.0-dz) + bchm_corvert(ilev2+1,ilev2+1+j,1:numvar3d)*dz)
            wtemp2(ilev1+j,ilev1,1:numvar3d)=wtemp2(ilev1,ilev1+j,1:numvar3d)
          end do               
       end if
    end do
    
    deallocate(bchm_invsum,bchm_corvert)
    
    allocate(bchm_corvert(nvlev,nvlev,numvar3d))
    allocate(bchm_invsum(nvlev,numvar3d))

    bchm_invsum(:,:) = wtemp1(:,:)
    bchm_corvert(:,:,:) = wtemp2(:,:,:)
    
    deallocate(wtemp1,wtemp2)
    
    lresetCorvert = .true.
    
  end subroutine bchm_resetCorvert

end module BmatrixChem_mod
