!--------------------------------------- LICENCE BEGIN -----------------------------------
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

MODULE BmatrixDiff_mod
  ! MODULE BmatrixDiff_mod (prefix='bdiff' category='5. B and R matrices')
  !
  ! :Purpose: Performs transformation from control vector to analysis increment 
  !           using the background-error covariance matrix based on correlations
  !           modelled using a diffusion operator.
  !
  use mpi_mod
  use mpivar_mod
  use MathPhysConstants_mod
  use earthConstants_mod
  use gridStateVector_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use varNameList_mod
  use utilities_mod
  use diffusion_mod
  implicit none
  save
  private

  ! public procedures
  public :: bdiff_Setup, bdiff_BSqrt, bdiff_BSqrtAd, bdiff_Finalize
  public :: bdiff_getScaleFactor

  logical             :: initialized = .false.
  integer             :: nj_l, ni_l
  integer             :: cvDim_mpilocal, cvDim_mpiglobal

  integer, allocatable :: diffID(:)

  ! Bacgkround-error covariance matrix elements.
  real(8), allocatable :: stddev(:,:,:)

  ! read in from the namelist:
  integer, parameter  :: maxNumVars = 200
  real(8)             :: scaleFactor(maxNumVars)
  real(8)             :: scaleFactor_sigma(maxNumVars)

  character(len=4)    :: stddevMode

  ! Homogeneous background-error standard deviation (when stddevMode == 'HOMO')
  real(8) :: homogeneous_std(maxNumVars)

  ! Number of incremental variables/fields
  integer             :: numvar2d
  ! Start position of each field in composite arrays
  integer, allocatable :: nsposit(:)
  ! Name list of incremental variables/fields
  character(len=4), allocatable :: bdiff_varNameList(:)

  integer             :: myLatBeg, myLatEnd
  integer             :: myLonBeg, myLonEnd

  integer,external    :: get_max_rss

  integer             :: nulbgst = 0

CONTAINS

  SUBROUTINE BDIFF_setup(hco_in, vco_in, CVDIM_OUT, mode_opt)
    implicit none

    type(struct_hco), pointer :: hco_in
    type(struct_vco), pointer :: vco_in
    integer, intent(out)      :: cvDim_out
    character(len=*), intent(in), optional :: mode_opt

    character(len=15) :: bdiff_mode

    integer :: nulnam, ierr, fnom, fclos
    integer :: status, latPerPE, latPerPEmax, lonPerPE, lonPerPEmax
    integer :: jvar

    type(struct_vco), pointer :: vco_anl

    ! namelist variables
    ! Horizontal correlation length scale (km)
    real :: corr_len(maxNumVars)
    ! Stability criteria (definitely < 0.5)
    real :: stab(maxNumVars)
    ! Number of samples in the estimation of the normalization factors by randomization.
    integer :: nsamp(maxNumVars)
    ! Indicate to use the implicit formulation of the diffusion operator (.true.) or
    ! the explicit version (.false.).
    logical :: limplicit(maxNumVars)

    NAMELIST /NAMBDIFF/ corr_len, stab, nsamp, limplicit, scaleFactor, stddevMode, homogeneous_std

    call tmg_start(17,'BDIFF_SETUP')
    if(mpi_myid == 0) write(*,*) 'bdiff_setup: starting'
    if(mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if ( present(mode_opt) ) then
       if ( trim(mode_opt) == 'Analysis' .or. trim(mode_opt) == 'BackgroundCheck') then
         bdiff_mode = trim(mode_opt)
         if(mpi_myid == 0) write(*,*)
         if(mpi_myid == 0) write(*,*) 'bmatrixDiff: Mode activated = ', trim(bdiff_mode)
       else
          write(*,*)
          write(*,*) 'mode = ', trim(mode_opt)
          call utl_abort('bmatrixDiff: unknown mode')
       end if
    else
       bdiff_mode = 'Analysis'
       if(mpi_myid == 0) write(*,*)
       if(mpi_myid == 0) write(*,*) 'bmatrixDiff: Analysis mode activated (by default)'
    end if

    vco_anl => vco_in
    if (vco_anl%Vcode  /= 5002 .and. vco_anl%Vcode /= 5005 .and. vco_anl%Vcode /= 0) then
      write(*,*) 'vco_anl%Vcode = ', vco_anl%Vcode
      call utl_abort('bmatrixDiff: unknown vertical coordinate type!')
    end if

    numvar2d = 0

    allocate(bdiff_varNameList(vnl_numvarmax))
    bdiff_varNameList(:)=''
    allocate(nsposit(vnl_numvarmax+1))
    nsposit(1) = 1

    ! Find the 2D variables (within NAMSTATE namelist)

    if(gsv_varExist(varName='GL  ')) then

       numvar2d = numvar2d + 1
       nsposit(numvar2d+1) = nsposit(numvar2d)+1
       bdiff_varNameList(numvar2d) = 'GL  '

    end if
    if(gsv_varExist(varName='TM  ')) then

       numvar2d = numvar2d + 1
       nsposit(numvar2d+1) = nsposit(numvar2d)+1
       bdiff_varNameList(numvar2d) = 'TM  '

    end if

    if (numvar2d == 0) then    
       if(mpi_myid == 0) then
          write(*,*) 'Bdiff matrix not produced.'
          write(*,*) 'END OF BDIFF_SETUP'
       end if
       call tmg_stop(17)
       cvdim_out = 0
       return
    else if (mpi_myid == 0) then
       write(*,*) 'BDIFF_setup: Number of 2D variables', numvar2d, bdiff_varNameList(1:numvar2d)
    end if

    ! default values for namelist variables
    corr_len(:) = 10.0
    stab(:)     = 0.2
    nsamp(:)    = 10000
    limplicit(:) = .false.
    scaleFactor(:) = 0.0d0
    stddevMode  = 'GD2D'
    homogeneous_std(:) = -1.0d0

    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=nambdiff,iostat=ierr)
    if(ierr /= 0) call utl_abort('bdiff_setup: Error reading namelist')
    if(mpi_myid == 0) write(*,nml=nambdiff)
    ierr = fclos(nulnam)

    if (sum(scaleFactor(:)) == 0.0d0 ) then
       if(mpi_myid == 0) write(*,*) 'bmatrixDiff: scaleFactor=0, skipping rest of setup'
       cvdim_out = 0
       call tmg_stop(17)
       return
    end if

    if ( trim(bdiff_mode) == 'BackgroundCheck' ) then
       cvDim_out = 9999 ! Dummy value > 0 to indicate to the background check (s/r compute_HBHT_ensemble) 
                        ! that Diff is used
       call tmg_stop(17)
       return
    end if

    ! Assumes the input 'scalefactor' is a scaling factor of the variances.

    do jvar = 1, numvar2d
       if(scaleFactor(jvar) > 0.0d0) then 
          scaleFactor_sigma(jvar) = sqrt(scaleFactor(jvar))
       else
          scaleFactor_sigma(jvar) = 0.0d0
       end if
    end do

    ni_l = hco_in%ni
    nj_l = hco_in%nj

    allocate(diffID(numvar2d))
    do jvar = 1, numvar2d
       diffID(jvar) = diff_setup(hco_in, corr_len(jvar), stab(jvar), nsamp(jvar), limplicit(jvar))
    end do

    call mpivar_setup_latbands(nj_l, latPerPE, latPerPEmax, myLatBeg, myLatEnd)
    call mpivar_setup_lonbands(ni_l, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd)

    ! compute mpilocal control vector size
    cvDim_mpilocal = ni_l*nj_l*numvar2d
    cvDim_out = cvDim_mpilocal

    ! also compute mpiglobal control vector dimension
    call rpn_comm_allreduce(cvDim_mpilocal,cvDim_mpiglobal,1,"mpi_integer","mpi_sum","GRID",ierr)

    allocate(stddev(ni_l, nj_l, numvar2d))

    call BDIFF_rdstats

    if(mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if(mpi_myid == 0) write(*,*) 'END OF BDIFF_SETUP'

    initialized = .true.

    call tmg_stop(17)

  END SUBROUTINE BDIFF_setup


  subroutine bdiff_getScaleFactor(scaleFactor_out)
    implicit none

    real(8), intent(out) :: scaleFactor_out(:)

    integer :: jvar

    do jvar = 1, numvar2d
       scaleFactor_out(jvar) = scaleFactor(jvar)
    end do

  end subroutine bdiff_getScaleFactor


  subroutine bdiff_rdstats
    !
    !:Purpose: To read background-error stats file.
    !
    implicit none

    integer :: ierr, nmax, fnom, fstouv, fstfrm, fclos
    integer :: jvar
    logical :: lExists
    character(len=12) :: bFileName1 = './bgstddev'
    character(len=8)  :: bFileName2 = './bgcov'

    if(stddevMode == 'GD2D') then

       inquire(file=bFileName1, exist=lExists)
       if ( lexists ) then
          ierr = fnom(nulbgst, bFileName1, 'RND+OLD+R/O', 0)
          if ( ierr == 0 ) then
             nmax = fstouv(nulbgst, 'RND+OLD')
          else
             call utl_abort('BDIFF_RDSTATS: ERROR OPENING FILE '//trim(bFileName1))
          end if
       else
          ! Assume background-error stats in file bgcov. 
          inquire(file=bFileName2, exist=lExists)  
          if (lexists) then 
             ierr = fnom(nulbgst, bFileName2, 'RND+OLD+R/O', 0)
             if ( ierr == 0 ) then
                nmax = fstouv(nulbgst, 'RND+OLD')
             else
                call utl_abort('BDIFF_RDSTATS: ERROR OPENING FILE '//trim(bFileName2))
             end if
          else
             call utl_abort('BDIFF_RDSTATS: NO BACKGROUND-ERROR STAT FILE!!')
          end if
       end if

       call BDIFF_rdstd

       ierr = fstfrm(nulbgst)
       ierr = fclos(nulbgst)

    elseif(stddevMode == 'HOMO') then

       do jvar = 1, numvar2d
          stddev(:,:,jvar) = homogeneous_std(jvar)
       end do
       
    else

       call utl_abort('BDIFF_RDSTATS: unknown stddevMode: '//trim(stddevMode))

    end if

    call BDIFF_scalestd

  end subroutine bdiff_rdstats

  subroutine bdiff_rdstd
    !
    !:Purpose: To read 2D stddev and store as 3D
    !
    implicit none

    integer :: jvar, in
    integer :: ikey
    real(8), allocatable :: rgsig2d(:,:)

    ! standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3
    integer :: idate(100)
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket

!   Reading the data

    idate(1) = -1
    ip1      = -1
    ip2      = -1
    ip3      = -1

    cletiket = 'STDDEV'
    cltypvar = ' '

    ! Reading for 2D variables

    do jvar = 1, numvar2d

       clnomvar = bdiff_varNameList(jvar)

       allocate(rgsig2d(ni_l,nj_l))
       rgsig2d(:,:) = 0.0D0

       ikey = utl_fstlir(rgsig2d(:,:),nulbgst,ini,inj,ink, &
                         idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)

       if (ikey < 0) then
          write(*,*) 'BDIFF_RDSTD: ',jvar, clnomvar, ikey
          call utl_abort(': BDIFF_RDSTD record not found')
       end if

       stddev(:,:,jvar) = rgsig2d(:,:)

       deallocate(rgsig2d)

    end do

  end subroutine bdiff_rdstd

  subroutine bdiff_scalestd
    !
    !:Purpose: To scale background-error standard-deviation values.
    !
    implicit none

    integer :: jlon, jlat, jvar

    do jvar = 1, numvar2d
       do jlat = 1, nj_l
          do jlon = 1, ni_l
             stddev(jlon,jlat,jvar) = scaleFactor_sigma(jvar)*       &
               stddev(jlon,jlat,jvar)
          end do
       end do
    end do

  end subroutine bdiff_scalestd


  SUBROUTINE BDIFF_bSqrt(controlVector_in, statevector)
    implicit none

    real(8),          intent(in)    :: controlVector_in(cvDim_mpilocal)
    type(struct_gsv), intent(inout) :: statevector

    real(8) :: gd_in( myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)
    real(8) :: gd_out(myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)

    integer :: jvar

    if(.not. initialized) then
      if(mpi_myid == 0) write(*,*) 'bdiff_bsqrt: bMatrixDIFF not initialized'
      return
    end if

    if(mpi_myid == 0) write(*,*) 'bdiff_bsqrt: starting'
    if(mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    call bdiff_cain(controlVector_in, gd_in)

    do jvar = 1, numvar2d

       ! Apply square root of the diffusion operator.
       call diff_Csqrt(diffID(jvar), gd_in(:,:,jvar), gd_out(:,:,jvar))

       ! Multiply by the diagonal matrix of background-error standard deviations.
       gd_out(:,:,jvar) = gd_out(:,:,jvar)*stddev(:,:,jvar)

    end do

    call copyToStatevector(statevector, gd_out)

    if(mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if(mpi_myid == 0) write(*,*) 'bdiff_bsqrt: done'

  END SUBROUTINE BDIFF_bSqrt


  SUBROUTINE BDIFF_bSqrtAd(statevector, controlVector_out)
    implicit none

    type(struct_gsv), intent(in)  :: statevector
    real(8),          intent(out) :: controlVector_out(cvDim_mpilocal)

    real(8) :: gd_in( myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)
    real(8) :: gd_out(myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)

    integer :: jvar

    if(.not. initialized) then
      if(mpi_myid == 0) write(*,*) 'bdiff_bsqrtad: bMatrixDIFF not initialized'
      return
    end if

    if(mpi_myid == 0) write(*,*) 'bdiff_bsqrtad: starting'
    if(mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    call copyFromStatevector(statevector, gd_in)

    do jvar = 1, numvar2d

       ! Multiply by the diagonal matrix of background-error standard deviations.
       gd_in(:,:,jvar) = gd_in(:,:,jvar)*stddev(:,:,jvar)

       ! Apply the adjoint of the square root of the diffusion operator.
       call diff_Csqrtadj(diffID(jvar), gd_in(:,:,jvar), gd_out(:,:,jvar))

    end do

    call bdiff_cainad(gd_out, controlVector_out)

    if(mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if(mpi_myid == 0) write(*,*) 'bdiff_bsqrtad: done'

  END SUBROUTINE BDIFF_bSqrtAd


  SUBROUTINE copyToStatevector(statevector, gd)
    implicit none

    real(8),          intent(in)    :: gd(myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)
    type(struct_gsv), intent(inout) :: statevector

    integer :: jlon, jlev, jlev2, jlat, jvar, ilev1, ilev2
    real(8), pointer :: field(:,:,:)

    do jvar = 1, numvar2d

       if(mpi_myid == 0) write(*,*) 'copyToStatevector: ',bdiff_varNameList(jvar)
       field => gsv_getField3D_r8(statevector, bdiff_varNameList(jvar))
       if(mpi_myid == 0) write(*,*) 'copyToStatevector: gsv_getField3D_r8 done.'

       ilev1 = nsposit(jvar)
       ilev2 = nsposit(jvar+1)-1 

!!!$OMP PARALLEL DO PRIVATE(jlat,jlev,jlev2,jlon)
       do jlev = ilev1, ilev2
          jlev2 = jlev-ilev1+1
          do jlat = myLatBeg, myLatEnd
             do jlon = myLonBeg, myLonEnd
                field(jlon,jlat,jlev2) = gd(jlon,jlat,jlev)
             end do
          end do
       end do
!!!$OMP END PARALLEL DO
    end do

  END SUBROUTINE copyToStatevector


  SUBROUTINE copyFromStatevector(statevector, gd)
    implicit none

    type(struct_gsv), intent(in)  :: statevector
    real(8),          intent(out) :: gd(myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)

    integer :: jlon, jlev, jlev2, jlat, jvar, ilev1, ilev2
    real(8), pointer :: field(:,:,:)

    do jvar = 1, numvar2d

       field => gsv_getField3D_r8(statevector, bdiff_varNameList(jvar))

       ilev1 = nsposit(jvar)
       ilev2 = nsposit(jvar+1)-1 

!!!$OMP PARALLEL DO PRIVATE(jlat,jlev,jlev2,jlon)
       do jlev = ilev1, ilev2
          jlev2 = jlev-ilev1+1
          do jlat = myLatBeg, myLatEnd
             do jlon = myLonBeg, myLonEnd
                gd(jlon,jlat,jlev) = field(jlon,jlat,jlev2)
             end do
          end do
       end do
!!!$OMP END PARALLEL DO
    end do

  END SUBROUTINE copyFromStatevector


  SUBROUTINE BDIFF_cain(controlVector_in, gd_out)
    implicit none

    real(8), intent(in)  :: controlVector_in(cvDim_mpilocal)
    real(8), intent(out) :: gd_out(myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)

    integer :: jn, jlev, jlon, jlat

    jn = 0
    do jlev = 1, numvar2d
      do jlat = myLatBeg, myLatEnd
        do jlon = myLonBeg, myLonEnd
           jn = jn + 1
           gd_out(jlon,jlat,jlev) = ControlVector_in(jn)
        end do
      end do
    end do

  end SUBROUTINE BDIFF_cain


  SUBROUTINE BDIFF_cainAd(gd_in, diffControlVector_out)
    implicit none

    real(8), intent(in)  :: gd_in(myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)
    real(8), intent(out) :: diffControlVector_out(cvDim_mpilocal)

    integer :: jn, jlev, jlon, jlat

    jn = 0
    do jlev = 1, numvar2d
      do jlat = myLatBeg, myLatEnd
        do jlon = myLonBeg, myLonEnd
           jn = jn + 1
           diffControlVector_out(jn) = gd_in(jlon,jlat,jlev)
        end do
      end do
    end do

  END SUBROUTINE BDIFF_cainAd


  SUBROUTINE BDIFF_Finalize()
    implicit none

    if (initialized) then
       initialized = .false.
       deallocate(stddev)
       deallocate(diffID)
       deallocate(nsposit)
       deallocate(bdiff_varNameList)
    end if

  END SUBROUTINE BDIFF_Finalize


END MODULE BmatrixDiff_mod
