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

MODULE BmatrixChem_mod 
  ! MODULE BmatrixChem_mod (prefix='bchm' category='5. B and R matrices')
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
  !   2. Currently assumes univariate constituent assimilation (constituent
  !      variables currently uncorrelated). See routines BCHM_READCORNS2 and
  !      BCHM_SUCORNS2. These  routines would need to be modified if
  !      cross-correlations were included.
  !
  !   3. One could potentially make public the functions/routines which are
  !      identical to those in bmatrixhi_mod.ftn90 (except possibly in name) so
  !      that one copy is present in the code.
  !
  !   4. For multiple univariate variables (or univarite blocks of one to
  !      multiple variables), one could alternatively have multiple sets of
  !      covariance matrices within this moduleinstead of a single covariance
  !      matrix setup (similarly to what was done for bchm_corvert*).
  !

  ! Public Subroutines (which call other internal routines/functions):
  !    BCHM_setup:      Must be called first. Sets of background covariance
  !                     matrix (and balance operators if any are eventually
  !                     added)
  !
  !    BCHM_BSqrt:      Transformations from control vector to analysis
  !                     increments in the minimization process.
  !
  !    BCHM_BSqrtAd:    Adjoint of BCHM_BSqrt.
  !    BCHM_Finalize    Deallocate internal module arrays.
  !    BCHM_corvert_mult: Multiple an input matrix/array with 'bchm_corvert' or
  !                     'bchm_corverti'
  !    BCHM_getsigma:   Obtain background error std. dev. profile at
  !                     obs/specified location. 
  !    BCHM_getBgSigma: Obtain background error std. dev. at specified grid
  !                     point for specified field.
  !    BCHM_is_initialized: checks if B_chm has been intialized.
  !    BCHM_StatsExistForVarname: Checfs if covariances available for specified
  !                     variable.

  use mpi_mod
  use mpivar_mod
  use MathPhysConstants_mod
  use earthConstants_mod
  use gridStateVector_mod
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
  real(8),allocatable :: corrlong(:,:,:),hwhm(:,:)

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
  logical             :: ReadWrite_sqrt
  character(len=4)    :: stddevMode
  character(len=4)    :: IncludeAnlVarCH(vnl_numvarmax)
 
! Number of incremental variables/fields
  integer             :: numvar3d,numvar2d
! Start position of each field in composite arrays
  integer, allocatable :: nsposit(:)
! Name list of incremental variables/fields
  character(len=4),allocatable    :: bchm_varNameList(:)


CONTAINS

  SUBROUTINE BCHM_setup(hco_in,vco_in,CVDIM_OUT,mode_opt)
    !
    !:Purpose: To set up for constituents static background error covariances.
    !
    implicit none

    type(struct_hco),pointer :: hco_in
    type(struct_vco),pointer :: vco_in
    integer                  :: cvDim_out
    character(len=*), intent(in), optional :: mode_opt

    integer :: jlev, nulnam, ierr, fnom, fclos, jm, jn, status
    integer :: latPerPE, latPerPEmax, lonPerPE, lonPerPEmax
    real(8) :: zps

    integer :: jvar,nChmVars,jvar2
    character(len=4) :: BchmVars(vnl_numvarmax)
    
    NAMELIST /NAMBCHM/ntrunc,rpor,rvloc,scaleFactor,numModeZero,ReadWrite_sqrt,stddevMode, &
                      IncludeAnlVarCH

   ! First check if there are any CH fields 
    
    jvar2=0
    do jvar = 1, vnl_numvarmax
      if (gsv_varExist(varName = vnl_varNameList(jvar))) then
        if (vnl_varKindFromVarname(vnl_varNameList(jvar)) == 'CH') then
          jvar2 = 1
          exit
        end if 
      end if      
    end do
    if (jvar2 == 0) then
       ! Assume there is no need for Bchm
       cvDim_out = 0
       return
    end if

    call tmg_start(120,'BCHM_SETUP')

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
    stddevMode = 'GD3D'    
    IncludeAnlVarCH(:) = ''
    
    ! Read namelist input
    
    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=NAMBCHM,iostat=ierr)
    if(ierr.ne.0) call utl_abort('BCHM_setup: Error reading namelist')
    if(mpi_myid == 0) write(*,nml=NAMBCHM)
    ierr = fclos(nulnam)

    ! Set BchmVars
    nChmVars=0
    BChmVars(:)=''
    if (trim(IncludeAnlVarCH(1)) == '') then
      do jvar = 1, vnl_numvarmax
        if (.not. gsv_varExist(varName = vnl_varNameList(jvar))) cycle
        if (vnl_varKindFromVarname(vnl_varNameList(jvar)) /= 'CH') cycle
        nChmVars = nChmVars+1
        BchmVars(nChmVars) = trim(vnl_varNameList(jvar))
      end do
    else
      do jvar = 1, vnl_numvarmax
        if (.not. gsv_varExist(varName = vnl_varNameList(jvar))) cycle
        if (vnl_varKindFromVarname(vnl_varNameList(jvar)) /= 'CH') cycle
        do jvar2 = 1, vnl_numvarmax
          if (trim(IncludeAnlVarCH(jvar2)) == '') exit
          if (trim(vnl_varNameList(jvar)) == trim(IncludeAnlVarCH(jvar2))) then
            nChmVars = nChmVars+1
            BchmVars(nChmVars) = trim(vnl_varNameList(jvar))
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
       call tmg_stop(120)
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
    if(mpi_myid == 0) write(*,*) 'BCHM_setup: nLev_T, nLev_T_even=',nLev_T, nLev_T_even

    !   Find the 3D variables (within NAMSTATE namelist)

    do jvar = 1, vnl_numvarmax3D    
       if (gsv_varExist(varName=vnl_varNameList3D(jvar)) .and. &
           any(trim(vnl_varNameList3D(jvar))==BchmVars(1:nChmVars)) ) then

           if (vnl_varKindFromVarname(vnl_varNameList3D(jvar)) /= 'CH') cycle
           numvar3d = numvar3d + 1
           nsposit(numvar3d+1)=nsposit(numvar3d)+nLev_T
           bchm_varNameList(numvar3d)=vnl_varNameList3D(jvar)
       end if
    end do

!   Find the 2D variables (within NAMSTATE namelist)

    do jvar = 1, vnl_numvarmax2D
      if (gsv_varExist(varName=vnl_varNameList2D(jvar)) .and. &
          any(trim(vnl_varNameList2D(jvar)) == BchmVars(1:nChmVars)) ) then

        if (vnl_varKindFromVarname(vnl_varNameList2D(jvar)) /= 'CH') cycle
        numvar2d = numvar2d + 1
        nsposit(numvar3d+numvar2d+1) = nsposit(numvar3d+numvar2d)+1
        bchm_varNameList(numvar2d) = vnl_varNameList2D(jvar)
      end if       
    end do
    
    if (numvar3d+numvar2d == 0) then    
      if (mpi_myid == 0) then
        write(*,*) 'Bhi matrix for CH family not produced.'
        write(*,*) 'No chemical assimilation to be performed.'
        write(*,*) 'END OF BCHM_SETUP'
      end if
      call tmg_stop(120)
      cvDim_out = 0
      return
    else if (mpi_myid == 0) then
      if (numvar3d > 0) &
        write(*,*) 'BCHM_setup: Number of 3D variables', numvar3d,bchm_varNameList(1:numvar3d)
      if (numvar2d > 0) &
        write(*,*) 'BCHM_setup: Number of 2D variables', numvar2d,bchm_varNameList(numvar3d+1:numvar3d+numvar2d)
    end if

    nkgdim =  nsposit(numvar3d+numvar2d+1)-1

    ! Initialization of namelist NAMBCHM parameters
    
    ntrunc = 108
    rpor(:) = 3000.D3
    rvloc(:) = 4.0D0
    scaleFactor(:,:) = 1.0d0
    numModeZero = 0
    ReadWrite_sqrt = .false.
    stddevMode = 'GD3D'
 
    ! Read namelist input
    
    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam, nml=NAMBCHM, iostat = ierr)
    if(ierr /= 0) call utl_abort('BCHM_setup: Error reading namelist')
    if(mpi_myid == 0) write(*, nml = NAMBCHM)
    ierr = fclos(nulnam)

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

    zps = 101000.D0
    status = vgd_levels( vco_anl%vgrid, ip1_list=vco_anl%ip1_M, levels=pressureProfile_M, &
                         sfc_field=zps, in_log=.false.)
    status = vgd_levels( vco_anl%vgrid, ip1_list=vco_anl%ip1_T, levels=pressureProfile_T, &
                         sfc_field=zps, in_log=.false.)

    call BCHM_rdstats
    call BCHM_sucorns2

    if(mpi_myid == 0) write(*,*) 'END OF BCHM_SETUP'
    
    initialized = .true.

    call tmg_stop(120)

  END SUBROUTINE BCHM_setup

  logical function BCHM_StatsExistForVarName(VarName)
    !
    !:Purpose: To check whether covariances have been made available for the
    !          specified variable
    !
    implicit none
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
    
  end function BCHM_StatsExistForVarName

  logical function BCHM_is_initialized()
    !
    !:Purpose: To check whether B_chm has been initialized
    !
    implicit none

    bchm_is_initialized = initialized

  end function BCHM_is_initialized

  subroutine BCHM_getScaleFactor(scaleFactor_out)
    !
    !:Purpose: To set scaling factors for background error std. dev.
    !
    implicit none
    real(8) :: scaleFactor_out(:,:)
    integer :: jlev,jvar

    do jvar = 1, numvar3d+numvar2d
    do jlev = 1, nsposit(jvar+1)-nsposit(jvar)
      scaleFactor_out(jlev,jvar) = scaleFactor(jlev,jvar)
    end do
    end do

  end subroutine BCHM_getScaleFactor

  subroutine BCHM_rdstats
    !
    !:Purpose: To read chemical constituents background stats file.
    !
   implicit none

    integer :: ierr, fnom, fstouv, fstfrm, fclos
    logical :: lExists
    character(len=12) :: bFileName1 = './bgchemcov'
    character(len=8)  :: bFileName2 = './bgcov'

    inquire(file=bFileName1,exist=lExists)
    if ( lexists ) then
      ierr = fnom(nulbgst,bFileName1,'RND+OLD+R/O',0)
      if ( ierr == 0 ) then
        ierr =  fstouv(nulbgst,'RND+OLD')
      else
        call utl_abort('BCHM_RDSTATS: NO BACKGROUND CHEMICAL CONSTITUENT STAT FILE!!')
      endif
    else
      ! Assume chemical constituent stats in file bgcov. 
      inquire(file=bFileName2,exist=lExists)  
      if (lexists) then 
        ierr = fnom(nulbgst,bFileName2,'RND+OLD+R/O',0)
        if ( ierr == 0 ) then
          ierr =  fstouv(nulbgst,'RND+OLD')
        else
          call utl_abort('BCHM_RDSTATS: NO BACKGROUND CHEMICAL CONSTITUENT STAT FILE!!')
        endif 
      else          
        call utl_abort('BCHM_RDSTATS: NO BACKGROUND CHEMICAL CONSTITUENT STAT FILE!!')
      end if
    endif

    call BCHM_readcorns2

    call BCHM_rdstddev 
    
    call BCHM_scalestd
    
    ierr = fstfrm(nulbgst)
    ierr = fclos(nulbgst)

  end subroutine BCHM_rdstats

  subroutine BCHM_scalestd
    !
    !:Purpose: To scale error standard-deviation values.
    !
    implicit none

    integer :: jlon, jlat, jvar, nlev

    do jvar = 1,numvar3d+numvar2d
      nlev=nsposit(jvar+1)-nsposit(jvar)
      do jlon = 1, ni_l+1
      do jlat = 1, nj_l
!        rgsig(jlat,nsposit(jvar):nsposit(jvar+1)-1) = &
!               scaleFactor_sigma(1:nlev,jvar)* &
!               rgsig(jlat,nsposit(jvar):nsposit(jvar+1)-1)
        rgsig(jlon,jlat,nsposit(jvar):nsposit(jvar+1)-1) = &
               scaleFactor_sigma(1:nlev,jvar)* &
               rgsig(jlon,jlat,nsposit(jvar):nsposit(jvar+1)-1)
      enddo
      enddo
    enddo

  end subroutine BCHM_scalestd

  subroutine BCHM_truncateCV(controlVector_inout,ntruncCut)
    !
    !:Purpose: To set to zero all coefficients with total wavenumber > ntruncCut
    !

    ! Based on bhi_truncateCV.
    implicit none

    real(8), pointer :: controlVector_inout(:)
    integer          :: ntruncCut
    integer          :: jn, jm, ila_mpiglobal, ila_mpilocal, jlev, jdim

    if(.not. initialized) then
      if(mpi_myid == 0) write(*,*) 'bchm_truncateCV: bMatrixChem not initialized'
      return
    endif

    if(ntruncCut.gt.ntrunc) then
      write(*,*) ntruncCut, ntrunc
      call utl_abort('bchm_truncateCV: ntruncCut is greater than ntrunc!')
    endif

    jdim = 0
    do jlev = 1, nkgdim
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
   
 end subroutine BCHM_truncateCV

 SUBROUTINE BCHM_bSqrt(controlvector_in,statevector)
    !
    !:Purpose: To apply B_CHM^1/2 to a control vector.
    !

    ! Based on bhi_bsqrt
   implicit none

    real(8)   :: controlVector_in(cvDim_mpilocal)
    type(struct_gsv) :: statevector
    real(8),allocatable :: gd_out(:,:,:)
    real(8)   :: hiControlVector(nla_mpilocal,2,nkgdim)
    integer   :: jvar, ilev1, ilev2

    if(.not. initialized) return

    if(mpi_myid == 0) write(*,*) 'BCHM_bsqrt: starting'
    if(mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    allocate(gd_out(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim))

    call bchm_cain(controlVector_in,hiControlVector)

    call bchm_spa2gd(hiControlVector,gd_out)
    
    call bchm_copyToStatevector(statevector,gd_out)

    deallocate(gd_out)

    if(mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if(mpi_myid == 0) write(*,*) 'bchm_bsqrt: done'

  end subroutine BCHM_bSqrt

!-----------------------------------------------------------------------------------------------

  subroutine BCHM_cain(controlVector_in,hiControlVector_out)
!
    implicit none

    real(8) :: controlVector_in(cvDim_mpilocal)
    real(8) :: hiControlVector_out(nla_mpilocal,2,nkgdim)

    integer :: jdim, jlev, jm, jn, ila_mpilocal, ila_mpiglobal

    jdim = 0
    hiControlVector_out(:,:,:) = 0.0d0
    do jlev = 1, nkgdim
      do jm = mymBeg, mymEnd, mymSkip
        do jn = mynBeg, mynEnd, mynSkip
          if(jm.le.jn) then
            ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
            ila_mpilocal  = ilaList_mpilocal(ila_mpiglobal)
            if(jm == 0) then
              ! only real component for jm=0
              jdim = jdim + 1
              hiControlVector_out(ila_mpilocal,1,jlev) = controlVector_in(jdim)
            else
              ! both real and imaginary components for jm>0
              jdim = jdim + 1
              hiControlVector_out(ila_mpilocal,1,jlev) = controlVector_in(jdim)
              jdim = jdim + 1
              hiControlVector_out(ila_mpilocal,2,jlev) = controlVector_in(jdim)
            endif
          endif
        enddo
      enddo
    enddo

  end SUBROUTINE BCHM_cain

!-----------------------------------------------------------------------------------------------

  SUBROUTINE BCHM_copyToStatevector(statevector,gd)
    implicit none
    type(struct_gsv) :: statevector
    real(8) :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)
    integer :: jlon, jlev, jlev2, jlat, jvar, ilev1, ilev2
    real(8), pointer :: field(:,:,:)

    do jvar = 1,numvar3d+numvar2d
      field => gsv_getField3D_r8(statevector,bchm_varNameList(jvar))
      ilev1 = nsposit(jvar)
      ilev2 = nsposit(jvar+1)-1 
        
!!!$OMP PARALLEL DO PRIVATE(jlat,jlev,jlev2,jlon)
      do jlev = ilev1, ilev2
        jlev2 = jlev-ilev1+1
        do jlat = myLatBeg, myLatEnd
          do jlon = myLonBeg, myLonEnd
            field(jlon,jlat,jlev2) = gd(jlon,jlat,jlev)
          enddo
        enddo
      enddo
!!!$OMP END PARALLEL DO
    enddo
  end subroutine BCHM_copyToStatevector

  SUBROUTINE BCHM_bSqrtAd(statevector,controlVector_out)
    !
    !:Purpose: To apply adjoint of B_CHM^1/2 to a statevector.
    !

    ! Based on bhi_bSqrtAd.
    !
    implicit none

    real(8)   :: controlVector_out(cvDim_mpilocal)
    type(struct_gsv) :: statevector
    real(8), allocatable :: gd_in(:,:,:)
    real(8)   :: hiControlVector(nla_mpilocal,2,nkgdim)
    integer   :: jvar, ilev1, ilev2

    if(.not. initialized) then
      if(mpi_myid == 0) write(*,*) 'bMatrixChem not initialized'
      return
    endif

    if(mpi_myid == 0) write(*,*) 'bchm_bsqrtad: starting'
    if(mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    allocate(gd_in(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim))
    call bchm_copyFromStatevector(statevector,gd_in)

    call bchm_spa2gdad(gd_in,hiControlVector)

    call bchm_cainad(hiControlVector,controlVector_out)

    deallocate(gd_in)

    if(mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if(mpi_myid == 0) write(*,*) 'BCHM_bsqrtad: done'

  end subroutine BCHM_bSqrtAd

!-----------------------------------------------------------------------------------------------

  SUBROUTINE BCHM_copyFromStatevector(statevector,gd)
    implicit none
    type(struct_gsv) :: statevector
    real(8)          :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)
    integer :: jlon, jlev, jlev2, jlat, jvar, ilev1, ilev2
    real(8), pointer :: field(:,:,:)

    do jvar = 1,numvar3d+numvar2d
      field => gsv_getField3D_r8(statevector,bchm_varNameList(jvar))

      ilev1 = nsposit(jvar)
      ilev2 = nsposit(jvar+1)-1 

!!!$OMP PARALLEL DO PRIVATE(jlat,jlev,jlev2,jlon)
      do jlev = ilev1, ilev2
        jlev2 = jlev-ilev1+1
        do jlat = myLatBeg, myLatEnd
          do jlon = myLonBeg, myLonEnd
            gd(jlon,jlat,jlev) = field(jlon,jlat,jlev2)
          enddo
        enddo
      enddo
!!!$OMP END PARALLEL DO
     enddo

  end subroutine BCHM_copyFromStatevector

  subroutine BCHM_readCorns2
    !
    !:Purpose: To read correlation information and to form the correlation
    !          matrix.
    !
    !:Notes: Currently assumes distinct block diagonal matrices (no
    !        cross-correlations)
    !

    ! Based on bhi_readcorns2.
    implicit none

    integer :: jn, istdkey,icornskey,jvar
    integer :: jcol,jrow,jstart,jnum
    real(8), allocatable, dimension(:) :: zstdsrc
    real(8), allocatable, dimension(:,:) :: zcornssrc

    ! Standard file variables
    integer :: ini,inj,ink,in1,in2
    integer :: ip1,ip2,ip3,idateo
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fstinf

    rstddev(:,:) = 0.0d0
    corns(:,:,:) = 0.0d0
    
    do jvar = 1, numvar3d+numvar2d
   
      clnomvar = bchm_varNameList(jvar)
      jstart = nsposit(jvar)
      jnum = nsposit(jvar+1)-nsposit(jvar)
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
          
        istdkey = utl_fstlir(ZSTDSRC,nulbgst,INI,INJ,INK,idateo,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
          
        if (istdkey < 0 .and. ip2 < 10) then
          write(*,*) 'BCHM_READCORNS2: RSTDDEV ',ip2,jnum,clnomvar
          call utl_abort('BCHM_READCORNS2: Problem with constituents background stat file')
        end if
        if (ini /= jnum)  call utl_abort('BCHM_READCORNS2: Constituents background stat levels inconsistencies')

        ! Looking for FST record parameters..

        if (istdkey >= 0) then
          idateo = -1
          cletiket = 'CORRNS'
          ip1 = -1
          ip2 = jn
          ip3 = -1
          cltypvar = 'X'
          icornskey = utl_fstlir(ZCORNSSRC,nulbgst,INI,INJ,INK,idateo,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)

          if (icornskey < 0) then
            write(*,*) 'BCHM_READCORNS2: CORRNS ',ip2,jnum,clnomvar
            call utl_abort('BCHM_READCORNS2: Problem with constituents background stat file')
          end if
          if (ini /= jnum .and. inj /= jnum) call utl_abort('BCHM_READCORNS2: Constituents BG stat levels inconsistencies')
        else
          write(*,*) 'WARNING from BCHM_READCORNS2: Setting RSDTDDEV to 1.D-15 for NOMVAR and JN: ',clnomvar,' ',jn
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

    ! Apply convolution to RSTDDEV correlation

    call BCHM_convol
    
    do jn = 0, ntrunc

      ! Re-build correlation matrix: factorization of corns with convoluted RSTDDEV
      do jcol = 1, nkgdim
        do jrow = 1, nkgdim
          corns(jrow,jcol,jn) = rstddev(jrow,jn) * corns(jrow,jcol,jn)* rstddev(jcol,jn)
        enddo
      enddo

    enddo

    !write(*,*) 'Done in BCHM_READCORNS2'
  end subroutine bchm_readCorns2

!-----------------------------------------------------------------------------------------------

  SUBROUTINE BCHM_convol
    implicit none

    real(8) dlfact2,dlc,dsummed
    real(8) dtlen,zr,dlfact
    integer ilen,jn,jlat,jk,jvar,jlev
    real(8) zleg(0:ntrunc,nj_l),zsp(0:ntrunc,nkgdim),zgr(nj_l,nkgdim)
    real(8) dlwti(nj_l),zrmu(nj_l)

    integer inracp
    real(8) zpg(nj_l),zsia(nj_l),zrad(nj_l),zpgssin2(nj_l)
    real(8) zsinm1(nj_l),zsinm2(nj_l),zsin2(nj_l),zsinlat(nj_l)
    real(8) dlfact1, dln
    real(8) dlnorm(0:ntrunc)

    logical lldebug
    
    lldebug=.true.
    
    do jlat = 1, nj_l
      dlwti(jlat) = gst_getrwt(jlat,gstID)
      zrmu(jlat)  = gst_getrmu(jlat,gstID)
    end do

!   CONVERT THE CORRELATIONS IN SPECTRAL SPACE INTO SPECTRAL
!   COEFFICIENTS OF THE CORRELATION FUNCTION AND FUNCTION TO BE
!   SELF-CONVOLVED

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
    
    jvar=1
    do jk = 1, nkgdim
      if (jk == nsposit(jvar)) then
        dtlen = rpor(jvar)
        jvar=jvar+1 
      endif
      if(dtlen.gt.0.0d0) then
        dlc = 1.d0/dble(dtlen)
        dlc = 0.5d0*dlc*dlc
        do jlat = 1, nj_l
          zr = ra * acos(zrmu(jlat))
          dlfact = dexp(-(zr**2)*dlc)
          zgr(jlat,jk) = dlfact*zgr(jlat,jk)
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

    ! Compute resultant physical space correlations

    if (allocated(corrlong)) deallocate(corrlong)
    allocate(corrlong(nj_l, nlev_T, numvar3d+numvar2d))
    if (allocated(hwhm)) deallocate(hwhm)
    allocate(hwhm(nlev_T, numvar3d+numvar2d))
    
    do jlat= 1, nj_l
    do jn= 0, ntrunc
      zleg(jn,jlat) = gst_getzleg(jn,jlat,gstID)
    end do
    end do
    
    jvar = 1
    jlev = 1
    do jk = 1, nkgdim
      if (jk == nsposit(jvar+1)) then
        jvar = jvar+1 
        jlev = 1
      end if
      do jlat = 1, nj_l
        corrlong(jlat,jlev,jvar) = 0.0D0
        do jn = 0, ntrunc
          corrlong(jlat,jlev,jvar) = corrlong(jlat,jlev,jvar)+rstddev(jk,jn)*rstddev(jk,jn)*  &
                                     sqrt(2.0)*sqrt(2.0*jn+1.0)*zleg(jn,jlat)
         end do
      enddo
      jlev = jlev+1 
    end do
   
   ! Get approx. half-width at half-max of correlation function
   
    jvar = 1
    jlev = 1
    do jk = 1, nkgdim
      if (jk == nsposit(jvar+1)) then
        jvar = jvar+1 
        jlev = 1
      end if
      do jlat= nj_l-1, 2, -1
        if (corrlong(jlat,jlev,jvar) <= 0.5) then
          hwhm(jlev,jvar) = ra*acos(zrmu(jlat))
          exit
        end if
      end do
      jlev = jlev+1
    end do  
    if(lldebug) then
      do jk = 1, numvar3d+numvar2d
        write(701,*)
        write(701,*) 'Horizontal correlations'
        write(701,*)
        write(701,*) 'Separation distances (km)'
        write(701,*) nj_l
        write(701,*) abs(ra*acos(zrmu(1:nj_l)))
        write(701,*)
        if (jk <= numvar3d) then
          write(701,*) bchm_varNameList(jvar), nlev_T
          do jlev = 1, nlev_T
            write(701,'(i3,f10.2,3000f6.2)') jlev,hwhm(jlev,jvar),corrlong(1:nj_l,jlev,jvar)
          end do
        else
          write(701,*) bchm_varNameList(jvar), 1
          write(701,'(i3,f10.2,3000f6.2)') 1,hwhm(jlev,jvar),corrlong(1:nj_l,jlev,jvar)
        end if
      end do
    endif

  end subroutine BCHM_convol

  subroutine BCHM_rdstddev
    !
    !:Purpose: To read stddev and to set as 3D fields.
    !
    implicit none

    integer :: ikey
    real(8) :: rcoord(10000)
    
    ! standard file variables
    integer :: ini,inj,ink, inpas, inbits, idatyp, ideet
    integer :: ip1,ip2,ip3,ig1,ig2,ig3,ig4,iswa,ilength,idltf
    integer :: iubc,iextr1,iextr2,iextr3,ipak,ipas
    integer :: iliste(100),idate(100),idimax,infon,iheures,idateo,nlev_MT
    character(len=1)  :: clgrtyp
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket

!   Reading the data

    idate(1) = -1
    ip2      = -1
    ip3      = -1

!   Get latitudes and longitudes if available

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
      call BCHM_rdstd3D
    elseif(stddevMode == 'GD2D') then
      call BCHM_rdstd
    elseif(stddevMode == 'SP2D') then
      call BCHM_rdspstd
    else
      call utl_abort('BCHM_RDSTDDEV: unknown stddevMode')
    endif
    
  end subroutine BCHM_rdstddev

!-----------------------------------------------------------------------------------------------

  SUBROUTINE BCHM_RDSPSTD
    implicit none

    integer :: jvar,jn,inix,injx,inkx
    integer :: ikey, jlevo, firstn,lastn
    real(8) :: zsp(0:ntrunc,max(nlev_M,nlev_T)),zspbuf(max(nlev_M,nlev_T))
    real(8) :: zgr(nj_l,max(nlev_M,nlev_T)),zstddev(nkgdim,nj_l)
    
    ! standard file variables
    integer :: ini,inj,ink, inpas, inbits, idatyp, ideet
    integer :: ip1,ip2,ip3,ig1,ig2,ig3,ig4,iswa,ilength,idltf
    integer :: iubc,iextr1,iextr2,iextr3,ipak,ipas
    integer :: iliste(100),idate(100),idimax,infon,iheures,idateo,nlev_MT
    character(len=1)  :: clgrtyp
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fstinf

    rgsig(:,:,:) = 0.0d0
    
!   Reading the Legendre poly representation of the 2D background error std. dev. field

    idate(1) = -1
    ip1      = -1
    ip2      = -1
    ip3      = -1

    cletiket = 'SPSTDDEV'
    cltypvar = 'X'
    
    do jvar = 1,numvar3d+numvar2d
      clnomvar = bchm_varNameList(jvar)
      nlev_MT = nsposit(jvar+1)-nsposit(jvar)
      
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
        
        do jlevo = 1, nlev_MT
          zsp(jn,jlevo) = zspbuf(jlevo)
        enddo
      enddo
      
      if(mpi_myid == 0 .and. firstn /= -1) write(*,*) 'WARNING: CANNOT FIND SPSTD FOR ',clnomvar, &
            ' AT N BETWEEN ',firstn,' AND ',lastn,', SETTING TO ZERO!!!'

      call gst_zleginv(gstID,zgr(:,1:nlev_MT),zsp(:,1:nlev_MT),nlev_MT)

      do jn = 1, ni_l+1
         rgsig(jn,1:nj_l,nsposit(jvar):nsposit(jvar+1)-1) = zgr(1:nj_l,1:nlev_MT)
      end do

    enddo 

  END SUBROUTINE BCHM_RDSPSTD

!-----------------------------------------------------------------------------------------------

  SUBROUTINE BCHM_RDSPSTD_NEWFMT
    implicit none

    integer :: jvar,jn,inix,injx,inkx,ntrunc_file
    integer :: ikey,jlevo
    real(8) :: zsp(0:ntrunc,max(nlev_M,nlev_T)),work
    real(8), allocatable :: zspbuf(:)
    real(8) :: zgr(nj_l,max(nlev_M,nlev_T)),zstddev(nkgdim,nj_l)

    ! standard file variables
    integer :: ini,inj,ink, inpas, inbits, idatyp, ideet
    integer :: ip1,ip2,ip3,ig1,ig2,ig3,ig4,iswa,ilength,idltf
    integer :: iubc,iextr1,iextr2,iextr3,ipak,ipas
    integer :: iliste(100),idate(100),idimax,infon,iheures,idateo,nlev_MT
    character(len=1)  :: clgrtyp
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fstinf

    rgsig(:,:,:) = 0.0d0

!   Reading the data

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
      write(*,*) 'BCHM_RDSPSTD_NEWFMT: ini>1, SPSTDDEV is in old format, calling BCHM_RDSPSTD...'
      call BCHM_rdspstd
      return
    endif

    ! write(*,*) 'Reading for 3D and 2D variables'
    
    do jvar = 1,numvar3d+numvar2d
      clnomvar = bchm_varNameList(jvar)
      nlev_MT = nsposit(jvar+1)-nsposit(jvar)

      !write(*,*) 'Reading ',clnomvar

      do jlevo = 1, nlev_MT
        if (nlev_MT == 1) then
           ip1 = -1
        else if(vnl_varLevelFromVarName(clnomvar) == 'MM') then
          ip1 = vco_anl%ip1_M(jlevo)
        else
          ip1 = vco_anl%ip1_T(jlevo)
        endif
        
        ikey = fstinf(nulbgst,inix,ntrunc_file,inkx,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
        ntrunc_file = ntrunc_file-1

        allocate(zspbuf(0:ntrunc_file))
        
        if(ikey >= 0 ) then
          ikey = utl_fstlir(zspbuf(0:ntrunc_file),nulbgst,ini,inj,ink,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
        else
          write(*,*) 'BCHM_RDSPSTD_NEWFMT: ',jvar,clnomvar,nlev_MT,jlevo,ikey,ntrunc,ntrunc_file
          call utl_abort('BCHM_RDSPSTD_NEWFMT: SPSTDDEV record not found')
        endif

        zsp(:,jlevo) = 0.0d0
        do jn = 0, min(ntrunc,ntrunc_file)
          zsp(jn,jlevo) = zspbuf(jn)
        enddo
        deallocate(zspbuf)
      enddo

      call gst_zleginv(gstID,zgr(:,1:nlev_MT),zsp(:,1:nlev_MT),nlev_MT)

      do jn = 1, ni_l+1
         rgsig(jn,1:nj_l,nsposit(jvar):nsposit(jvar+1)-1) = zgr(1:nj_l,1:nlev_MT)
      end do

    enddo

  end subroutine BCHM_rdspstd_newfmt

  subroutine BCHM_rdstd
    !
    !:Purpose: To read 2D stddev and to store as 3D
    !
    implicit none

    integer :: jvar,in
    integer :: ikey,jlevo
    real(8), allocatable :: rgsig3d(:,:,:)
    
    ! standard file variables
    integer :: ini,inj,ink, inpas, inbits, idatyp, ideet
    integer :: ip1,ip2,ip3,ig1,ig2,ig3,ig4,iswa,ilength,idltf
    integer :: iubc,iextr1,iextr2,iextr3,ipak,ipas
    integer :: iliste(100),idate(100),idimax,infon,iheures,idateo,nlev_MT
    character(len=1)  :: clgrtyp
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fstinf
    
    real(8), allocatable  :: vlev(:),vlevout(:)

    rgsig(:,:,:) = 0.0d0

!   Reading the data

    idate(1) = -1
    ip1      = -1
    ip2      = -1
    ip3      = -1
    
    cletiket = 'STDDEV'
    cltypvar=' '

    ! Reading for 3D and 2D variables
    
    do jvar = 1,numvar3d+numvar2d
      clnomvar = bchm_varNameList(jvar)
      nlev_MT = nsposit(jvar+1)-nsposit(jvar)

      !write(*,*) 'Reading ',clnomvar
 
      ikey = fstinf(nulbgst,ini,inj,ink,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
      if (ikey < 0 .or. ini > 1 .or. ink /= nlev_MT) then
        write(*,*) 'BCHM_RDSTD: ',jvar,clnomvar,ikey,ini,ink,nlev_MT
        call utl_abort(': BCHM_RDSTD record not found or incorrect')          
      end if
      
      allocate(rgsig3d(1,inj,ink))
      rgsig3d(1,:,:) = 0.0D0
      allocate(vlev(ink),vlevout(nlev_MT))
      vlev(:) = 1.0D0
      vlevout(:) = 1.0D0
      
      ikey = utl_fstlir(rgsig3d(1,:,:), nulbgst, ini, inj, ink, &
                         idate(1), cletiket, ip1, ip2, ip3, cltypvar, clnomvar)

      if (ikey < 0) then
        write(*,*) 'BCHM_RDSTD: ',jvar,clnomvar,nlev_MT,ikey
        call utl_abort(': BCHM_RDSTD record not found')
      endif
        
      ! Extend to 3D
      if (inj == nj_l) then
        do in = 1, ni_l+1
          rgsig(in,:,nsposit(jvar):nsposit(jvar+1)-1) = rgsig3d(1,:,:) 
        end do
      else
         ! Interpolate in lat
         call gsv_field3d_hbilin(rgsig3d, 1, inj, ink, rlongr, rlatr, vlev, &
              rgsig(:,:,nsposit(jvar):nsposit(jvar+1)-1), ni_l+1, nj_l, nlev_MT, &
              rlong, rlat, vlevout)
      end if
      deallocate(rgsig3d, vlev, vlevout)
       
    enddo

  end subroutine bchm_rdstd

  subroutine bchm_rdstd3d
    !
    !:Purpose: To read 3D stddev.
    !

    ! originally based on bchm_rdspstd_newfmt
    !
    implicit none

    integer :: jvar
    integer :: ikey,jlevo
    real(8), allocatable :: rgsig3d(:,:,:)
    
    ! standard file variables
    integer :: ini,inj,ink, inpas, inbits, idatyp, ideet
    integer :: ip1,ip2,ip3,ig1,ig2,ig3,ig4,iswa,ilength,idltf
    integer :: iubc,iextr1,iextr2,iextr3,ipak,ipas
    integer :: iliste(100),idate(100),idimax,infon,iheures,idateo,nlev_MT
    character(len=1)  :: clgrtyp
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fstinf

    real(8) :: vlev(1),vlevout(1)

    rgsig(:,:,:) = 0.0d0
    vlev(:) = 1.0D0
    vlevout(:) = 1.0D0

!   Reading the data

    idate(1) = -1
    ip2      = -1
    ip3      = -1
    
    cletiket = 'STDDEV3D'
    cltypvar=' '

    ! Reading for 3D and 2D variables
    
    do jvar = 1,numvar3d+numvar2d
      clnomvar = bchm_varNameList(jvar)
      nlev_MT = nsposit(jvar+1)-nsposit(jvar)

      !write(*,*) 'Reading ',clnomvar

      do jlevo = 1, nlev_MT
        if (nlev_MT == 1) then
          ip1 = -1
        else if(vnl_varLevelFromVarName(clnomvar) == 'MM') then
          ip1 = vco_anl%ip1_M(jlevo)
        else
          ip1 = vco_anl%ip1_T(jlevo)
        endif
        
        ikey = fstinf(nulbgst,ini,inj,ink,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
        if (ikey < 0) then
            write(*,*) 'BCHM_RDSTD: ',jvar,clnomvar,ikey,jlevo
            call utl_abort(': BCHM_RDSTD record not foun0d')          
        end if
      
        allocate(rgsig3d(ini+1, inj, 1))
        rgsig3d(:,:,1) = 0.0D0
        
        ikey = utl_fstlir(rgsig3d(1:ini,:,1), nulbgst, ini, inj, ink, &
                         idate(1), cletiket, ip1, ip2, ip3, cltypvar, clnomvar)
        if (ikey < 0) then
          write(*,*) 'BCHM_RDSTD3D: ',jvar,clnomvar,nlev_MT,jlevo,ikey,ip1
          call utl_abort(': BCHM_RDSTD3D record not found')
        endif
        
        ! Move to rgsig
        if (inj == nj_l .and. ini == ni_l) then
          rgsig(1:ni_l,:,nsposit(jvar)+(jlevo-1)) = rgsig3d(:,:,1) 
          rgsig(ni_l+1,:,nsposit(jvar)+(jlevo-1)) = rgsig3d(1,:,1)
        else
          ! Interpolate in lat and long
          rgsig3d(ini+1,:,1) = rgsig3d(1,:,1)
          call gsv_field3d_hbilin(rgsig3d, ini+1, inj, 1, rlongr, rlatr, vlev, &
               rgsig(:,:,nsposit(jvar)+(jlevo-1)), ni_l+1, nj_l, 1, &
               rlong, rlat, vlevout)
       end if
       ! write(*,*) 'STDDDEV ',jlevo,rgsig(1,1,nsposit(jvar)+(jlevo-1)),rgsig(ni_l+1,1,nsposit(jvar)+(jlevo-1)),rgsig(ni_l/2,nj_l/2,nsposit(jvar)+(jlevo-1)),rgsig(ni_l+1,nj_l,nsposit(jvar)+(jlevo-1))
       deallocate(rgsig3d)

      enddo
    enddo

  end subroutine bchm_rdstd3d

!-----------------------------------------------------------------------------------------------

  SUBROUTINE BCHM_SUCORNS2
    implicit none

    real(8) :: eigenval(nkgdim), eigenvec(nkgdim,nkgdim),result(nkgdim,nkgdim)
    real(8) :: eigenvalsqrt(nkgdim)

    integer :: jlat,jn,jk1,jk2,jk3,jr,jvar
    integer :: ilwork,info,jnum
    integer :: ikey, nsize

    real(8) :: zwork(2*4*nkgdim)
    real(8) :: ztlen,zcorr,zr,zpres1,zpres2,eigenvalmax
    real(8),allocatable :: corns_temp(:,:,:)
    logical lfound_sqrt
    
    ! Standard file variables
    integer :: ini,inj,ink, inpas, inbits, idatyp, ideet
    integer :: ip1,ip2,ip3,ig1,ig2,ig3,ig4,iswa,ilength,idltf
    integer :: iubc,iextr1,iextr2,iextr3,ierr,ipak,ipas,ntrials
    integer :: iliste(100),idate(100),idimax,infon,iheures,idateo
    character(len=2)  :: cltypvar
    character(len=1)  :: clgrtyp
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket

    ! Apply vertical localization to correlations of 3D fields.
    ! Currently assumes no-cross correlations for variables (block diagonal matrix)
    
    do jvar = 1, numvar3d
      ztlen = rvloc(jvar)    ! specify length scale (in units of ln(Pressure))
        
      jnum = nsposit(jvar+1)-nsposit(jvar)
       
      if(ztlen > 0.0d0) then
        ! calculate 5'th order function (from Gaspari and Cohn)
        do jk1 = 1, jnum
          zpres1 = log(pressureProfile_T(jk1))
          do jk2 = 1, jnum
            zpres2 = log(pressureProfile_T(jk2))
            zr = abs(zpres2 - zpres1)
            zcorr = gasparicohn(ztlen,zr)
            corns(nsposit(jvar)-1+jk1,nsposit(jvar)-1+jk2,0:ntrunc) = &
                  corns(nsposit(jvar)-1+jk1,nsposit(jvar)-1+jk2,0:ntrunc)*zcorr  
          enddo
        enddo
      endif
    enddo

    ! Compute total vertical correlations and its inverse (currently for each block matrix).

    call bchm_corvert_setup
    
    if (trim(bchm_mode) == 'BackgroundCheck') return

    ! Following assumes full matrices. It does not take advantage of any block diagonal structure.
    ! Accounting for block diagonal structures would/might improve computation time.
    
    ! compute square-root of corns for each total wavenumber
    
    allocate(corns_temp(nkgdim,nkgdim,0:ntrunc))
    corns_temp(:,:,:) = 0.0d0
    do jn = mpi_myid, ntrunc, mpi_nprocs

      eigenvec(1:nkgdim,1:nkgdim) = corns(1:nkgdim,1:nkgdim,jn)

      ! CALCULATE EIGENVALUES AND EIGENVECTORS.
      ilwork = 4*nkgdim*2
      call dsyev('V','U',nkgdim,eigenvec,nkgdim,eigenval,zwork,ilwork,info)
      if(info /= 0) then
        write(*,*) 'BCHM_sucorns2: non-zero value of info =',info,' returned by dsyev for wavenumber ',jn
        call utl_abort('BCHM_SUCORNS')
      endif
      
      ! set selected number of eigenmodes to zero
      if(numModeZero > 0) then
!        write(*,*) 'bchm_sucorns2: setting ',numModeZero,' eigenvalues to zero for wavenumber n=',jn
!        write(*,*) 'bchm_sucorns2: original eigenvalues=',eigenval(:)
        do jk1 = 1, numModeZero
          eigenval(jk1) = 0.0d0
        enddo
!        write(*,*) 'bchm_sucorns2: modified eigenvalues=',eigenval(:)
      endif

      eigenvalmax = maxval(eigenval(1:jnum))
      do jk1 = 1, nkgdim
!        if(eigenval(jk1) < 1.0d-15) then
        if(eigenval(jk1) < 1.0d-8*eigenvalmax) then
          eigenvalsqrt(jk1) = 0.0d0
        else
          eigenvalsqrt(jk1) = sqrt(eigenval(jk1))
        endif
      enddo

      ! E * lambda^1/2
      do jk1 = 1, nkgdim
         result(1:nkgdim,jk1) = eigenvec(1:nkgdim,jk1)*eigenvalsqrt(jk1)
      enddo

      ! (E * lambda^1/2) * E^T
      do jk1 = 1, nkgdim
      do jk2 = 1, nkgdim
         do jk3 = 1, nkgdim
           corns_temp(jk2,jk1,jn) = corns_temp(jk2,jk1,jn) + result(jk2,jk3)*eigenvec(jk1,jk3)
         enddo
      enddo
      enddo

    enddo ! jn

    nsize = nkgdim*nkgdim*(ntrunc+1)
    call rpn_comm_allreduce(corns_temp,corns,nsize,"mpi_double_precision","mpi_sum","GRID",ierr)
    deallocate(corns_temp)

    if(ReadWrite_sqrt) then
      ! if desired, read precomputed sqrt of corns which overwrites computed matrix
      call readcorns_sqrt(lfound_sqrt)
      if(.not.lfound_sqrt) then
        ! if precomputed sqrt does not exist in stats, then write it out to a separate file
        call writecorns_sqrt
      endif
    endif

  END SUBROUTINE BCHM_SUCORNS2

!-----------------------------------------------------------------------------------------------

  SUBROUTINE WRITECORNS_SQRT
    implicit none

    integer :: jn, nulcorns_sqrt, ierr

    ! standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3
    integer :: idateo, ipak, idatyp
    character(len=2)  :: cltypvar
    character(len=1)  :: clgrtyp
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fnom, fstouv, fstfrm, fclos

    write(*,*) 'WRITECORNS_SQRT: CORNS_SQRT will be written to file corns_sqrt_chm.fst for NTRUNC =', ntrunc

    if(mpi_myid==0) then
      nulcorns_sqrt = 0
      ierr = fnom(nulcorns_sqrt,'corns_sqrt_chm.fst','RND',0)
      ierr = fstouv(nulcorns_sqrt,'RND')

      ipak = -32
      idatyp = 5
      ip1 = 0
      ip3 = ntrunc
      idateo = 0

      do jn = 0, ntrunc
        ip2 = jn
        ierr = utl_fstecr(corns(1,1,jn),ipak,nulcorns_sqrt,idateo,0,0,nkgdim,nkgdim,1,  &
                       ip1,ip2,ip3,'X','ZZ','CORNS_SQRT','X',0,0,0,0,idatyp,.true.)
      enddo

      ierr = fstfrm(nulcorns_sqrt)  
      ierr = fclos(nulcorns_sqrt)
    endif

  END SUBROUTINE WRITECORNS_SQRT

!-----------------------------------------------------------------------------------------------

  SUBROUTINE READCORNS_SQRT(lfound_sqrt)
    implicit none
    logical :: lfound_sqrt

    integer :: jn, icornskey
    integer :: jcol,jrow
    real(8), allocatable :: zcornssrc(:,:)

    ! standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3
    integer :: idateo
    character(len=2)  :: cltypvar
    character(len=1)  :: clgrtyp
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket

    allocate(zcornssrc(nkgdim,nkgdim))

    do jn = 0, ntrunc

      ! Looking for FST record parameters..

      write(*,*) 'READCORNS_SQRT: looking for CORNS_SQRT with NTRUNC =', ntrunc
      idateo = -1
      cletiket = 'CORNS_SQRT'
      ip1 = -1
      ip2 = jn
      ip3 = ntrunc
      cltypvar = 'X'
      clnomvar = 'ZZ'
      icornskey = utl_fstlir(ZCORNSSRC,nulbgst,INI,INJ,INK,idateo,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)

      if(icornskey .lt.0 ) then
        write(*,*) 'READCORNS_SQRT: CORNS_SQRT not found in stats file, use computed sqrt'
        lfound_sqrt = .false.
        return
      else
        write(*,*) 'READCORNS_SQRT: CORNS_SQRT found in stats file, will use it instead of computed sqrt'
        lfound_sqrt = .true.
      endif

      if (ini .ne. nkgdim .or. inj .ne. nkgdim) then
        call utl_abort('READCORNS_SQRT: BG stat levels inconsitencies')
      endif

      do jcol = 1, nkgdim
        do jrow = 1, nkgdim
          corns(jrow,jcol,jn) = zcornssrc(jrow,jcol)
        enddo
      enddo

    enddo

    deallocate(zcornssrc)

  END SUBROUTINE READCORNS_SQRT

!-----------------------------------------------------------------------------------------------

  FUNCTION GASPARICOHN(ztlen,zr)

    real(8)  :: gasparicohn
    real(8)  :: ztlen,zr,zlc

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

  END FUNCTION GASPARICOHN

!-----------------------------------------------------------------------------------------------

  SUBROUTINE BCHM_SPA2GD(hiControlVector_in,gd_out)
    IMPLICIT NONE

    real(8) :: hiControlVector_in(nla_mpilocal,2,nkgdim)
    real(8) :: gd_out(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)

    real(8) :: sp(nla_mpilocal,2,nkgdim)

    integer :: jn,jm,ila_mpilocal,ila_mpiglobal,icount
    real(8) :: sq2, zp
    real(8) , allocatable :: zsp(:,:,:), zsp2(:,:,:)
    integer :: ilon, jlev, jlon, jlat, jla_mpilocal
    real(8), target  :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)
    real(8) :: dla2, dl1sa2, zpsb(myLonBeg:myLonEnd,myLatBeg:myLatEnd)

    call tmg_start(82,'BCHM_SPA2GD1')

    ! maybe not needed:
    sp(:,:,:) = 0.0d0
    sq2 = sqrt(2.0d0)

    allocate(zsp(nkgdim,2,mymCount))
    allocate(zsp2(nkgdim,2,mymCount))

!$OMP PARALLEL DO PRIVATE(jn,jm,jlev,ila_mpiglobal,ila_mpilocal,zsp2,zsp,icount)
    do jn = mynBeg, mynEnd, mynSkip

      icount = 0
      do jm = mymBeg, mymEnd, mymSkip
        if(jm.le.jn) then
          icount = icount+1
          ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
          ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
          do jlev = 1, nkgdim
            zsp(jlev,1,icount) = hiControlVector_in(ila_mpilocal,1,jlev)
            zsp(jlev,2,icount) = hiControlVector_in(ila_mpilocal,2,jlev)
          enddo
        endif
      enddo
      if(icount.gt.0) then

        CALL DGEMM('N','N',nkgdim,2*icount,nkgdim,1.0d0,corns(1,1,jn),nkgdim,zsp(1,1,1),nkgdim,0.0d0,zsp2(1,1,1),nkgdim)
!        CALL DGEMUL(corns(1,1,jn),nkgdim,'N',zsp(1,1,1),nkgdim,'N',zsp2(1,1,1),nkgdim,nkgdim,nkgdim,2*icount)

        icount = 0
        do jm = mymBeg, mymEnd, mymSkip
          if(jm.le.jn) then
            icount = icount+1
            ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
            ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
            do jlev = 1, nkgdim
              sp(ila_mpilocal,1,jlev) = zsp2(jlev,1,icount)
              sp(ila_mpilocal,2,jlev) = zsp2(jlev,2,icount)
            enddo
          endif
        enddo

      endif

      ! make adjustments for jm=0
      if(mymBeg == 0) then

        ila_mpiglobal = gst_getNind(0,gstID) + jn
        ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)

        do jlev = 1, nkgdim
          sp(ila_mpilocal,1,jlev) = sp(ila_mpilocal,1,jlev)*sq2
          sp(ila_mpilocal,2,jlev) = 0.0d0
        enddo

      endif

    enddo
!$OMP END PARALLEL DO
    deallocate(zsp)
    deallocate(zsp2)
    call tmg_stop(82)


!$OMP PARALLEL DO PRIVATE(JLAT,JLEV,JLON)
    do jlev = 1, nkgdim
       do jlat = myLatBeg, myLatEnd
          do jlon = myLonBeg, myLonEnd
             gd(jlon,jlat,jlev) = 0.0d0
          enddo
       enddo
    enddo
!$OMP END PARALLEL DO

    call tmg_start(127,'BCHM_SPEREE')
    call gst_setID(gstID)
    call gst_speree(sp,gd)
    call gst_setID(gstID2)
    call tmg_stop(127)

    call tmg_start(85,'BCHM_SPA2GD2')

!$OMP PARALLEL DO PRIVATE(jlat,jlev,jlon)
    do jlev = 1, nkgdim
       do jlat = myLatBeg, myLatEnd
          do jlon = myLonBeg, myLonEnd
             gd(jlon,jlat,jlev) = gd(jlon,jlat,jlev)*rgsig(jlon,jlat,jlev)
          enddo
       enddo
    enddo
!$OMP END PARALLEL DO

    call tmg_stop(85)

!$OMP PARALLEL DO PRIVATE(JLAT,JLEV,JLON)
    do jlev = 1, nkgdim
       do jlat = myLatBeg, myLatEnd
          do jlon = myLonBeg, myLonEnd
             gd_out(jlon,jlat,jlev) = gd(jlon,jlat,jlev)
          enddo
       enddo
    enddo
!$OMP END PARALLEL DO

  END SUBROUTINE BCHM_SPA2GD

!-----------------------------------------------------------------------------------------------

  SUBROUTINE BCHM_SPA2GDAD(gd_in,hiControlVector_out)
    implicit none

    real(8) :: hiControlVector_out(nla_mpilocal,2,nkgdim)
    real(8) :: gd_in(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)

    real(8) :: sp(nla_mpilocal,2,nkgdim)

    integer :: jn, jm, ila_mpilocal, ila_mpiglobal, icount
    real(8) :: sq2, zp
    real(8), allocatable :: zsp(:,:,:), zsp2(:,:,:)

    integer :: ilon, jlev, jlon, jlat, jla_mpilocal
    real(8) :: dl1sa2, dla2, zpsb(myLonBeg:myLonEnd,myLatBeg:myLatEnd)
    real(8), target :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)


!$OMP PARALLEL DO PRIVATE(JLAT,JLEV,JLON)
    do jlev = 1, nkgdim
       do jlat = myLatBeg, myLatEnd
          do jlon = myLonBeg, myLonEnd                                                      
             gd(jlon,jlat,jlev) = gd_in(jlon,jlat,jlev)
          enddo
       enddo
    enddo
!$OMP END PARALLEL DO

   call tmg_start(85,'BCHM_SPA2GD2')

!$OMP PARALLEL DO PRIVATE(jlat,jlev,jlon)
   do jlev = 1, nkgdim
      do jlat = myLatBeg, myLatEnd
         do jlon = myLonBeg, myLonEnd
            gd(jlon,jlat,jlev) = gd(jlon,jlat,jlev)*rgsig(jlon,jlat,jlev)
         enddo
      enddo
   enddo
!$OMP END PARALLEL DO

    call tmg_stop(85)

    call tmg_start(86,'BCHM_REESPE')
    call gst_setID(gstID)
    call gst_speree_ad(sp,gd)
    call tmg_stop(86)

    call tmg_start(82,'BCHM_SPA2GD1')

    hiControlVector_out(:,:,:) = 0.0d0
    sq2 = sqrt(2.0d0)
    allocate(zsp(nkgdim,2,mymCount))
    allocate(zsp2(nkgdim,2,mymCount))
    
!$OMP PARALLEL DO PRIVATE(JN,JM,JLEV,ILA_MPILOCAL,ILA_MPIGLOBAL,zsp,zsp2,icount)
    do jn = mynBeg, mynEnd, mynSkip

      icount = 0
      do jm = mymBeg, mymEnd, mymSkip
        if(jm.le.jn) then
          icount = icount+1
          ila_mpiglobal = gst_getNind(jm,gstID) + jn - jm
          ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
          do jlev = 1, nkgdim
            zsp2(jlev,1,icount) = sp(ila_mpilocal,1,jlev)
            zsp2(jlev,2,icount) = sp(ila_mpilocal,2,jlev)
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
          do jlev = 1, nkgdim
            hiControlVector_out(ila_mpilocal,1,jlev) = zsp(jlev,1,icount)
            hiControlVector_out(ila_mpilocal,2,jlev) = zsp(jlev,2,icount)
          enddo
        enddo

      endif

      ! make adjustments for jm=0
      if(mymBeg == 0) then

        ila_mpiglobal = gst_getNIND(0,gstID) + jn
        ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)

       do jlev = 1, nkgdim
          hiControlVector_out(ila_mpilocal,1,jlev) = hiControlVector_out(ila_mpilocal,1,jlev)*sq2
          hiControlVector_out(ila_mpilocal,2,jlev) = hiControlVector_out(ila_mpilocal,2,jlev)*sq2
        enddo

      endif

    enddo
!$OMP END PARALLEL DO
    deallocate(zsp)
    deallocate(zsp2)
    call tmg_stop(82)

  END SUBROUTINE BCHM_SPA2GDAD

!----------------------------------------------------------------------------------------------- 

  SUBROUTINE BCHM_cainAd(hiControlVector_in,controlVector_out)
    IMPLICIT NONE

    real(8) :: controlVector_out(cvDim_mpilocal)
    real(8) :: hiControlVector_in(nla_mpilocal,2,nkgdim)

    integer :: jdim, jlev, jm, jn, ila_mpilocal, ila_mpiglobal

    jdim = 0
    do jlev = 1, nkgdim
      do jm = mymBeg, mymEnd, mymSkip
        do jn = mynBeg, mynEnd, mynSkip
          if(jm.le.jn) then
            ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
            ila_mpilocal  = ilaList_mpilocal(ila_mpiglobal)
            if(jm == 0) then
              ! only real component for jm=0
              jdim = jdim + 1
              controlVector_out(jdim) = controlVector_out(jdim) + hiControlVector_in(ila_mpilocal,1,jlev)
            else
              ! both real and imaginary components for jm>0
              jdim = jdim + 1
              controlVector_out(jdim) = controlVector_out(jdim) + hiControlVector_in(ila_mpilocal,1,jlev)*2.0d0
              jdim = jdim + 1
              controlVector_out(jdim) = controlVector_out(jdim) + hiControlVector_in(ila_mpilocal,2,jlev)*2.0d0
            endif
          endif
        enddo
      enddo
    enddo

  END SUBROUTINE BCHM_cainAd

!-----------------------------------------------------------------------------------------------

  SUBROUTINE BCHM_Finalize()
    implicit none

    if(initialized) then
       deallocate(pressureProfile_M)
       deallocate(pressureProfile_T)
       deallocate(rgsig)
       deallocate(corns)
       deallocate(rstddev)
       deallocate(bchm_corvert,bchm_corverti,bchm_invsum)
    end if

  END SUBROUTINE BCHM_Finalize

!-----------------------------------------------------------------------------------------------

  SUBROUTINE BCHM_reduceToMPILocal(cv_mpilocal,cv_mpiglobal)
    implicit none
    real(8), intent(out) :: cv_mpilocal(cvDim_mpilocal)
    real(8), intent(in)  :: cv_mpiglobal(:)

    real(8), allocatable :: cv_allMaxMpilocal(:,:)

    integer, allocatable :: cvDim_allMpilocal(:), displs(:)
    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)

    integer :: jproc,cvDim_maxmpilocal,ierr
    integer :: jlev,jn,jm,ila_mpilocal,ila_mpiglobal,jdim_mpilocal,jdim_mpiglobal

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

!$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,jlev,jm,jn,ila_mpiglobal,jdim_mpiglobal)
       do jproc = 0, (mpi_nprocs-1)
          cv_allmaxmpilocal(:,jproc+1) = 0.d0
          jdim_mpilocal = 0

          do jlev = 1, nkgdim
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
                      jdim_mpiglobal = jdim_mpiglobal + (jlev-1) * (ntrunc+1)*(ntrunc+1)
                      
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
                         write(*,*) '       proc, jlev, jn, jm = ',jproc, jlev, jn, jm
                         call utl_abort('bhi_reduceToMPILocal')
                      end if
                      if (jdim_mpiglobal > cvDim_mpiglobal) then
                         write(*,*)
                         write(*,*) 'ERROR: jdim_mpiglobal > cvDim_mpiglobal', jdim_mpiglobal, cvDim_mpiglobal
                         write(*,*) '       proc, jlev, jn, jm = ',jproc, jlev, jn, jm
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

  END SUBROUTINE BCHM_reduceToMPILocal

!-----------------------------------------------------------------------------------------------

  SUBROUTINE BCHM_reduceToMPILocal_r4(cv_mpilocal,cv_mpiglobal)
    implicit none
    real(4), intent(out) :: cv_mpilocal(cvDim_mpilocal)
    real(4), intent(in)  :: cv_mpiglobal(:)

    real(4), allocatable :: cv_allMaxMpilocal(:,:)

    integer, allocatable :: cvDim_allMpilocal(:), displs(:)
    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)

    integer :: jproc,cvDim_maxmpilocal,ierr
    integer :: jlev,jn,jm,ila_mpilocal,ila_mpiglobal,jdim_mpilocal,jdim_mpiglobal

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

!$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,jlev,jm,jn,ila_mpiglobal,jdim_mpiglobal)
       do jproc = 0, (mpi_nprocs-1)
          cv_allmaxmpilocal(:,jproc+1) = 0.d0
          jdim_mpilocal = 0

          do jlev = 1, nkgdim
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
                      jdim_mpiglobal = jdim_mpiglobal + (jlev-1) * (ntrunc+1)*(ntrunc+1)
                      
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
                         write(*,*) '       proc, jlev, jn, jm = ',jproc, jlev, jn, jm
                         call utl_abort('bhi_reduceToMPILocal')
                      end if
                      if (jdim_mpiglobal > cvDim_mpiglobal) then
                         write(*,*)
                         write(*,*) 'ERROR: jdim_mpiglobal > cvDim_mpiglobal', jdim_mpiglobal, cvDim_mpiglobal
                         write(*,*) '       proc, jlev, jn, jm = ',jproc, jlev, jn, jm
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

  END SUBROUTINE BCHM_reduceToMPILocal_r4

!-----------------------------------------------------------------------------------------------

  SUBROUTINE BCHM_expandToMPIGlobal(cv_mpilocal,cv_mpiglobal)
    implicit none
    real(8), intent(in)  :: cv_mpilocal(cvDim_mpilocal)
    real(8), intent(out) :: cv_mpiglobal(:)

    real(8), allocatable :: cv_maxmpilocal(:)
    real(8), pointer :: cv_allmaxmpilocal(:,:) => null()
    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)
    integer :: jlev, jn, jm, jproc, ila_mpiglobal, jdim_mpilocal, jdim_mpiglobal, ierr, cvDim_maxmpilocal

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

    call tmg_start(128,'BCHM_COMM')
    call rpn_comm_gather(cv_maxmpilocal,    cvDim_maxmpilocal, "mpi_double_precision",  &
                         cv_allmaxmpilocal, cvDim_maxmpilocal, "mpi_double_precision", 0, "GRID", ierr )
    call tmg_stop(128)

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

!$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,jlev,jm,jn,ila_mpiglobal,jdim_mpiglobal)
      do jproc = 0, (mpi_nprocs-1)
        jdim_mpilocal = 0

        do jlev = 1, nkgdim
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
                jdim_mpiglobal = jdim_mpiglobal + (jlev-1) * (ntrunc+1)*(ntrunc+1)

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
                  write(*,*) 'ERROR: jdim,cvDim,mpiglobal=',jdim_mpiglobal,cvDim_mpiglobal,jlev,jn,jm

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

  end SUBROUTINE BCHM_expandToMPIGlobal

!-----------------------------------------------------------------------------------------------

  SUBROUTINE BCHM_expandToMPIGlobal_r4(cv_mpilocal,cv_mpiglobal)
    implicit none
    real(4), intent(in)  :: cv_mpilocal(cvDim_mpilocal)
    real(4), intent(out) :: cv_mpiglobal(:)

    real(4), allocatable :: cv_maxmpilocal(:)
    real(4), pointer :: cv_allmaxmpilocal(:,:) => null()
    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)
    integer :: jlev, jn, jm, jproc, ila_mpiglobal, jdim_mpilocal, jdim_mpiglobal, ierr, cvDim_maxmpilocal

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

    call tmg_start(128,'BCHM_COMM')
    call rpn_comm_gather(cv_maxmpilocal,    cvDim_maxmpilocal, "mpi_real4",  &
                         cv_allmaxmpilocal, cvDim_maxmpilocal, "mpi_real4", 0, "GRID", ierr )
    call tmg_stop(128)

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

!$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,jlev,jm,jn,ila_mpiglobal,jdim_mpiglobal)
      do jproc = 0, (mpi_nprocs-1)
        jdim_mpilocal = 0

        do jlev = 1, nkgdim
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
                jdim_mpiglobal = jdim_mpiglobal + (jlev-1) * (ntrunc+1)*(ntrunc+1)

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
                  write(*,*) 'ERROR: jdim,cvDim,mpiglobal=',jdim_mpiglobal,cvDim_mpiglobal,jlev,jn,jm

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

  subroutine bchm_corvert_setup
    !
    !:Purpose: To compute total vertical correlations (bchm_corvert) and its
    !          inverse (bchm_corverti; currently for each block matrix).
    !
    !
    !:Note: Currently assumes no cross-correlations 
    !
    implicit none

    real(8) :: eigenval(nkgdim), eigenvec(nkgdim,nkgdim),result(nkgdim,nkgdim)

    integer :: jn,jk1,jk2,jk3,jvar
    integer :: ilwork,info,jnum
    integer :: ikey
    logical :: lldebug

    real(8) :: zwork(2*4*nkgdim)
    real(8) :: zr,eigenvalmax
    integer iulcorvert
    
    ! Standard file variables
    integer :: ini,inj,ink, inpas, inbits, idatyp, ideet
    integer :: ip1,ip2,ip3,ig1,ig2,ig3,ig4,iswa,ilength,idltf
    integer :: iubc,iextr1,iextr2,iextr3,ierr,ipak,ipas,ntrials
    integer :: iliste(100),idate(100),idimax,infon,iheures,idateo
    character(len=2)  :: cltypvar
    character(len=1)  :: clgrtyp
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fstprm,fstinf
    integer :: fnom,fstouv,fstfrm,fclos

    lldebug=.false.
    
    ! Compute total vertical correlations and its inverse (currently for each block matrix).

    if (mpi_myid == 0) then
      iulcorvert = 0
      ierr = fnom(iulcorvert,'corvert_modular_chm.fst','RND',0)
      ierr = fstouv(iulcorvert,'RND')
    end if
    
    do jvar = 1, numvar3d
   
      jnum = nsposit(jvar+1)-nsposit(jvar)
       
      bchm_corvert(1:jnum,1:jnum,jvar) = 0.0d0
      do jn = 0, ntrunc
        bchm_corvert(1:jnum,1:jnum,jvar) = bchm_corvert(1:jnum,1:jnum,jvar)+ ((2*jn+1)* &
             corns(nsposit(jvar):nsposit(jvar+1)-1,nsposit(jvar):nsposit(jvar+1)-1,jn))
      enddo

      !  where (abs(bchm_corvert(1:jnum,1:jnum,jvar)) .lt. 1.E-3) bchm_corvert(1:jnum,1:jnum,jvar)=0.0D0

      if (mpi_myid == 0) then
        ikey = fstinf(NULBGST,ini,inj,ink,-1,'CORRNS',-1,0,-1,' ',bchm_varNameList(jvar))
        ierr = fstprm(ikey,idateo,ideet,inpas,ini,inj,ink, inbits        &
             ,idatyp,ip1,ip2,ip3,cltypvar,clnomvar,cletiket,clgrtyp      &
             ,ig1,ig2,ig3,ig4,iswa,ilength,idltf,iubc,iextr1,iextr2      &
             ,iextr3)

        ini = jnum
        inj = jnum
        ink = 1
        ip1 = 0
        ip2 = ntrunc
        ip3 = 0
        clnomvar = bchm_varNameList(jvar)
        cletiket = 'CORVERT'
        idatyp = 5

        ierr = utl_fstecr(bchm_corvert(1:jnum,1:jnum,jvar) &
             , -inbits, iulcorvert, idateo    &
             , ideet,inpas, ini, inj, ink, ip1, ip2, ip3, cltypvar,     &
             clnomvar,cletiket,clgrtyp,ig1, ig2, ig3, ig4, idatyp,      &
             .true.)

        if (lldebug) then
          write(701,*)
          write(701,*) bchm_varNameList(jvar)
          write(701,*) 'Total correlation matrix'
          write(701,*) bchm_corvert(1:jnum,1:jnum,jvar)
        endif
                      
      endif
       
      ! Inverse not needed if not in analysis mode
      if (trim(bchm_mode) == 'BackgroundCheck') cycle
       
      eigenvec(1:jnum,1:jnum)=bchm_corvert(1:jnum,1:jnum,jvar)       

      ! CALCULATE EIGENVALUES AND EIGENVECTORS.
      ilwork = 4*jnum*2
      call dsyev('V','U',jnum,eigenvec(1:jnum,1:jnum),jnum,eigenval,zwork,ilwork,info)
      if (info.ne.0) then
        write(*,*) 'BCHM_corvert_setup: non-zero value of info =',info,' returned by dsyev for wavenumber ',jn
        call utl_abort('BCHM_corvert_setup')
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
      bchm_corverti(1:jnum,1:jnum,jvar)=0.0D0
      do jk1 = 1, jnum
      do jk2 = 1, jnum
        do jk3 = 1, jnum
          bchm_corverti(jk2,jk1,jvar) = bchm_corverti(jk2,jk1,jvar) + result(jk2,jk3)*eigenvec(jk1,jk3)
        enddo
      enddo
      enddo

      ! Check inverse (for output when mpi_myid is 0)
      result(1:jnum,1:jnum)=0.0D0
      do jk1 = 1, jnum
      do jk2 = 1, jnum
        do jk3 = 1, jnum
          result(jk2,jk1) = result(jk2,jk1) + &
            bchm_corvert(jk2,jk3,jvar)*bchm_corverti(jk3,jk1,jvar)
        enddo
      enddo
      enddo

      ! zr=maxval(abs(bchm_corverti(1:jnum,1:jnum,jvar)))
      ! where (abs(bchm_corverti(1:jnum,1:jnum,jvar)) .lt. 1.E-5*zr) bchm_corverti(1:jnum,1:jnum,jvar)=0.0D0
            
      if (mpi_myid == 0) then

        cletiket = 'CORVERTI'

        ierr = utl_fstecr(bchm_corverti(1:jnum,1:jnum,jvar) &
             ,  -inbits, iulcorvert, idateo    &
             , ideet,inpas, ini, inj, ink, ip1, ip2, ip3, cltypvar,     &
             clnomvar,cletiket,clgrtyp,ig1, ig2, ig3, ig4, idatyp,      &
             .true.)

        cletiket = 'C*C^-1'
          
        ierr = utl_fstecr(result(1:jnum,1:jnum) &
             ,  -inbits, iulcorvert, idateo    &
             , ideet,inpas, ini, inj, ink, ip1, ip2, ip3, cltypvar,     &
             clnomvar,cletiket,clgrtyp,ig1, ig2, ig3, ig4, idatyp,      &
             .true.)

        if (lldebug) then
           write(701,*) 'Inverse'
           write(701,*) bchm_corverti(1:jnum,1:jnum,jvar)
           write(701,*) 'Product'
           write(701,*) result(1:jnum,1:jnum)
        endif
                     
      end if

      ! Generate 1/sum(covert(:,i,jvar)
       
      do jk2 = 1, jnum  
        bchm_invsum(jk2,jvar) = 1.0D0/sum(bchm_corvert(1:jnum,jk2,jvar))
      end do

    end do
    
    if (mpi_myid == 0) then
      ierr = fstfrm(iulcorvert)
      ierr = fclos(iulcorvert)
    end if

  end subroutine bchm_corvert_setup

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

    ! Arguments:
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
   
    integer :: jvar,jk1,jk2,jk3,nsize
    real(8) :: rmat_work(ndim2,ndim2),rsum

    if (.not.initialized) return
 
    if (.not.present(rsig_opt).and.lrgsig) call utl_abort('BCHM_corvert_mult: Missing rsig_opt')  

				! Determine location and size in bchm_corvert/bchm_corverti
    
    do jvar = 1, numvar3d+numvar2d
      if (trim(varName) == trim(bchm_varNameList(jvar))) exit
    end do
    if  (jvar > numvar3d+numvar2d) call utl_abort('BCHM_corvert_mult: Variable not found')
    
    nsize = nsposit(jvar+1)-nsposit(jvar)
    
			! Apply matrix/vector multiplication.
    
    rmat_work(:,:) = 0.0D0

    if (itype == 3) then

      !  A*C*A^T
       
      if (ndim3 /= ndim1 .or. ndim2 /= nsize) call utl_abort('BCHM_corvert_mult: Size does not match - itype=3')
      if (lrgsig) then
        do jk1 = 1, nsize
        do jk2 = 1, ndim1
          rmat_work(jk1,jk2) = sum(rmat_in(jk2,lvl_top(jk2):lvl_bot(jk2))*bchm_corvert(jk1,lvl_top(jk2):lvl_bot(jk2),jvar) &
                                *rsig_opt(lvl_top(jk2):lvl_bot(jk2)))*rsig_opt(jk1)
        end do
        end do
      else
        do jk1 = 1, nsize
        do jk2 = 1, ndim1
          rmat_work(jk1,jk2) = sum(rmat_in(jk2,lvl_top(jk2):lvl_bot(jk2))*bchm_corvert(jk1,lvl_top(jk2):lvl_bot(jk2),jvar))
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
       
      if (ndim3 /= ndim1 .or. ndim2 /= nsize) call utl_abort('BCHM_corvert_mult: Size does not match - itype=3')
      if (lrgsig) then
        do jk1 = 1, nsize
        do jk2 = 1, ndim1
          rmat_work(jk1,jk2) = sum(rmat_in(jk2,lvl_top(jk2):lvl_bot(jk2))*bchm_corverti(jk1,lvl_top(jk2):lvl_bot(jk2),jvar) &
                                /rsig_opt(lvl_top(jk2):lvl_bot(jk2)))/rsig_opt(jk1)
        end do
        end do
      else
        do jk1 = 1, nsize
        do jk2 = 1, ndim1
          rmat_work(jk1,jk2) = sum(rmat_in(jk2,lvl_top(jk2):lvl_bot(jk2))*bchm_corverti(jk1,lvl_top(jk2):lvl_bot(jk2),jvar))
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

      if (ndim3 /= ndim2 .or. ndim1 /= nsize) call utl_abort('BCHM_corvert_mult: Size does not match - itype=2')
      if (lrgsig) then
        do jk1 = 1, nsize
        do jk3 = 1, nsize
          do jk2 = lvl_top(jk3), lvl_bot(jk3)     ! Instead of do jk2=1,ndim2
            rmat_out(jk1,jk2) = rmat_out(jk1,jk2)+bchm_corvert(jk1,jk3,jvar)*rsig_opt(jk3)*rsig_opt(jk1)*rmat_in(jk3,jk2)
          end do
        end do
        end do
      else
        do jk1 = 1, nsize
        do jk3 = 1, nsize
          do jk2 = lvl_top(jk3), lvl_bot(jk3)     ! Instead of do jk2=1,ndim2
            rmat_out(jk1,jk2) = rmat_out(jk1,jk2)+bchm_corvert(jk1,jk3,jvar)*rmat_in(jk3,jk2)
          end do
        end do
        end do
      endif

    else if (itype == -2) then
    
      ! CI*A

      if (ndim3 /= ndim2 .or. ndim1 /= nsize) call utl_abort('BCHM_corvert_mult: Size does match - itype=-2')
      if (lrgsig) then
        do jk1 = 1, nsize
        do jk3 = 1, nsize
          do jk2 = lvl_top(jk3), lvl_bot(jk3)     ! Instead of do jk2=1,ndim2
            rmat_out(jk1,jk2) = rmat_out(jk1,jk2)+bchm_corverti(jk1,jk3,jvar)*rmat_in(jk3,jk2)/(rsig_opt(jk3)*rsig_opt(jk1))
          end do
        end do
        end do
      else
        do jk1 = 1, nsize
        do jk3 = 1, nsize
          do jk2 = lvl_top(jk3), lvl_bot(jk3)     ! Instead of do jk2=1,ndim2
            rmat_out(jk1,jk2) = rmat_out(jk1,jk2)+bchm_corverti(jk1,jk3,jvar)*rmat_in(jk3,jk2)
          end do
        end do
        end do
      endif
          
    else if (itype == 1) then
    
      ! A*C

      if (ndim3 /= ndim2 .or. ndim2 /= nsize) call utl_abort('BCHM_corvert_mult: Size does not match - itype=1')
      if (lrgsig) then
        do jk1 = 1, ndim1
        do jk2 = 1, nsize
          rmat_out(jk1,jk2) = sum(rmat_in(jk1,lvl_top(jk1):lvl_bot(jk1))*bchm_corvert(lvl_top(jk1):lvl_bot(jk1),jk2,jvar) &
                               *rsig_opt(lvl_top(jk1):lvl_bot(jk1)))*rsig_opt(jk2)
        end do
        end do
      else
        do jk1 = 1, ndim1
        do jk2 = 1, nsize
          rmat_out(jk1,jk2) = sum(rmat_in(jk1,lvl_top(jk1):lvl_bot(jk1))*bchm_corvert(lvl_top(jk1):lvl_bot(jk1),jk2,jvar))
        end do
        end do
      endif
          
    else if (itype == -1) then

      ! A*CI

      if (ndim3 /= ndim2 .or. ndim2 /= nsize) call utl_abort('BCHM_corvert_mult: Size does not match - itype=-1')
      if (lrgsig) then
        do jk1 = 1, ndim1
        do jk2 = 1, nsize
          rmat_out(jk1,jk2) = sum(rmat_in(jk1,lvl_top(jk1):lvl_bot(jk1))*bchm_corverti(lvl_top(jk1):lvl_bot(jk1),jk2,jvar) &
                                 /rsig_opt(lvl_top(jk1):lvl_bot(jk1)))/rsig_opt(jk2)
        end do
        end do
      else
        do jk1 = 1, ndim1
        do jk2 = 1, nsize
          rmat_out(jk1,jk2) = sum(rmat_in(jk1,lvl_top(jk1):lvl_bot(jk1))*bchm_corverti(lvl_top(jk1):lvl_bot(jk1),jk2,jvar))
        end do
        end do
      endif

    else if (itype == 0) then

      ! Instead of A*CI do D(i,j)=A(i,j)/sum(C(1:n,i))

      if (ndim3 /= ndim2 .or. ndim2 /= nsize) call utl_abort('BCHM_corvert_mult: Size does not match - itype=0')
      if (lrgsig) then
        do jk1 = 1, ndim1
        do jk2 = lvl_top(jk1), lvl_bot(jk1)   ! instead of do jk2=1,nsize
          rmat_out(jk1,jk2) = rmat_in(jk1,jk2)/rsig_opt(jk2) &
              /sum(bchm_corvert(1:nsize,jk2,jvar)/rsig_opt(1:nsize)) 
        end do
        end do
      else
        do jk1 = 1, ndim1 
          rmat_out(jk1,lvl_top(jk1):lvl_bot(jk1)) = rmat_in(jk1,lvl_top(jk1):lvl_bot(jk1))*bchm_invsum(lvl_top(jk1):lvl_bot(jk1),jvar)            
!          do jk2 = lvl_top(jk1), lvl_bot(jk1)   ! instead of do jk2=1,nsize
!              rmat_out(jk1,jk2) = rmat_in(jk1,jk2)/sum(bchm_corvert(1:nsize,jk2,jvar))
!          end do
        end do
      endif

    else
      call utl_abort('BCHM_corvert_mult: Requested type not found')
    end if
   
  end subroutine bchm_corvert_mult

  subroutine bchm_getsigma(varName,ndim2,xlat,xlong,rsig)
    !
    !:Purpose: To interpolate error std. dev. to obs location.
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: varName ! Variable name
    integer, intent(in) :: ndim2    ! Max array size
    real(8), intent(in) :: xlat     ! Target latitude
    real(8), intent(in) :: xlong    ! Target longitude
    real(8), intent(out) :: rsig(:) ! Error std. dev.
   
    ! Locals:
    integer :: jvar,jlev1,jlev2,jlat,jlong,j,nsize
    real(8) :: rvar(ndim2,2),zc1,zc2,rlat1,rlat2,rlong1,rlong2,zd1,zd2,rsig_max

    integer, parameter :: itype=0

    if (.not.initialized) return
    
   ! Determine location and size of background error std. dev.

    do jvar = 1, numvar3d+numvar2d
      if (trim(varName) == trim(bchm_varNameList(jvar))) exit
    end do
    if  (jvar > numvar3d+numvar2d) call utl_abort('BCHM_getsigma: Variable not found')
    
    jlev1 = nsposit(jvar)
    jlev2 = nsposit(jvar+1)-1

    nsize = jlev2-jlev1+1
    if (nsize /= ndim2) then
      write(6,*) 'NSIZE, NDIM2: ',nsize,ndim2
      call utl_abort('BCHM_getsigma: Inconsistent size')
    end if
 
    ! Determine reference longitude index 

    jlong = 2    
    do while (xlong > rlong(jlong) .and. jlong < ni_l+1) 
      jlong = jlong+1
    end do

    ! Set interpolation weights

    rlong2 = rlong(jlong)
    rlong1 = rlong(jlong-1)
    
    zd2 = (xlong-rlong1)/(rlong2-rlong1)
    zd1 = 1.0-zd2
     
    ! Determine reference latitude index (could alternatively use gst_getrlati)

    jlat = 2 
    do while (xlat > rlat(jlat) .and. jlat < nj_l) 
      ! Note: gst_getrlati needs to be consistent with rlat (except start of indexing and extrapolated lats)
      jlat = jlat+1
    enddo

    ! Set interpolation weights

    rlat2 = rlat(jlat)
    rlat1 = rlat(jlat-1)
    
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

      rvar(1:nsize,1) = zc1*rgsig(jlong-1,jlat-1,jlev1:jlev2)**2 + zc2*rgsig(jlong-1,jlat,jlev1:jlev2)**2
    
      rvar(1:nsize,2) = zc1*rgsig(jlong,jlat-1,jlev1:jlev2)**2 + zc2*rgsig(jlong,jlat,jlev1:jlev2)**2
       
      rsig(1:nsize) = sqrt( zd1*rvar(1:nsize,1) + zd2*rvar(1:nsize,2) )
    
    else 

      ! Interpolation of std. dev. (to reduce execution time)

      rvar(1:nsize,1) = zc1*rgsig(jlong-1,jlat-1,jlev1:jlev2) + zc2*rgsig(jlong-1,jlat,jlev1:jlev2)     

      rvar(1:nsize,2) = zc1*rgsig(jlong,jlat-1,jlev1:jlev2) + zc2*rgsig(jlong,jlat,jlev1:jlev2)

      rsig(1:nsize) = zd1*rvar(1:nsize,1) + zd2*rvar(1:nsize,2)
    
    end if
    
    do j = 1, nsize
      rsig_max = maxval(rgsig(jlong-1:jlong,jlat-1:jlat,jlev1-1+j))
      if (rsig(j) < 0.0 .or. rsig(j) > 1.00001*rsig_max) then
        write(6,*) 'bchm_getsigma: Interpolated sigma incorrect:'
        write(6,*) 'bchm_getsigma:   zc1,zc2,zd1,zd2 = ',zc1,zc2,zd1,zd2
        write(6,*) 'bchm_getsigma:   rgsig = ',rgsig(jlong-1,jlat-1,jlev1-1+j),rgsig(jlong,jlat-1,jlev1-1+j),rgsig(jlong-1,jlat,jlev1-1+j),rgsig(jlong,jlat,jlev1-1+j)
        write(6,*) 'bchm_getsigma:   rsig,rsig_max = ',rsig(j),rsig_max
        write(6,*) 'bchm_getsigma:   jlong,jlat,j,xlong,xlat = ',jlong,jlat,j,xlong,xlat
        write(6,*) 'bchm_getsigma:   rlong2,rlong1,rlat1,rlat2 = ',rlong2,rlong1,rlat1,rlat2
        call utl_abort('bchm_getsigma')
      end if
    end do

  end subroutine bchm_getsigma

  real(8) function BCHM_getbgsigma(jlon,jlat,jlev,jvar)
    !
    !:Purpose: To get error std. dev. a specified grid point and for specified
    !          field
    !
    implicit none

    integer :: jlon, jlat, jvar, jlev

    bchm_getbgsigma = rgsig(jlon, jlat, nsposit(jvar)-1+jlev)

  end function BCHM_getbgsigma
!-----------------------------------------------------------------------------------------------
  
end module BmatrixChem_mod
