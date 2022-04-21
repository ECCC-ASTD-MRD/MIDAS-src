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
  integer             :: latPerPE, latPerPEmax, lonPerPE, lonPerPEmax

  integer,external    :: get_max_rss

  integer             :: nulbgst = 0

CONTAINS

  subroutine bdiff_setup ( hco_in, vco_in, cvDim_out, mode_opt )
  
    implicit none

    ! Arguments:
    type(struct_hco), intent(inout), pointer  :: hco_in
    type(struct_vco), intent(inout), pointer  :: vco_in
    integer         , intent(out)             :: cvDim_out
    character(len=*), intent(in)   , optional :: mode_opt
    
    ! locals:
    character(len=15)         :: bdiff_mode
    type(struct_vco), pointer :: vco_anl
    integer                   :: nulnam, ierr, fnom, fclos
    integer                   :: variableIndex
    ! namelist variables
    real    :: corr_len( maxNumVars )  ! Horizontal correlation length scale (km)
    real    :: stab( maxNumVars )      ! Stability criteria (definitely < 0.5)
    integer :: nsamp(maxNumVars)       ! Number of samples in the estimation of the normalization factors by randomization.
    logical :: useImplicit(maxNumVars) ! Indicate to use the implicit formulation of the diffusion operator (.true.) or
                                       ! the explicit version (.false.).
    character(len=*), parameter :: myName = 'bdiff_setup'
    
    NAMELIST /NAMBDIFF/ corr_len, stab, nsamp, useImplicit, scaleFactor, stddevMode, homogeneous_std

    call utl_tmg_start(65,'----B_DIFF_Setup')
    if(mpi_myid == 0) write(*,*) myName//': starting'
    if(mpi_myid == 0) write(*,*) myName//': Memory Used: ',get_max_rss()/1024,'Mb'

    if ( present( mode_opt ) ) then
      if ( trim( mode_opt ) == 'Analysis' .or. trim( mode_opt ) == 'BackgroundCheck' ) then
        bdiff_mode = trim( mode_opt )
        if( mpi_myid == 0 ) write(*,*)
        if( mpi_myid == 0 ) write(*,*) myName//': Mode activated = ', trim(bdiff_mode)
      else
        write(*,*)
        write(*,*)  myName//'mode = ', trim(mode_opt)
        call utl_abort( myName//': unknown mode' )
      end if
    else
      bdiff_mode = 'Analysis'
      if( mpi_myid == 0 ) write(*,*)
      if( mpi_myid == 0 ) write(*,*) myName//': analysis mode activated (by default)'
    end if

    vco_anl => vco_in
    if (vco_anl%Vcode  /= 5002 .and. vco_anl%Vcode /= 5005 .and. vco_anl%Vcode /= 0 ) then
      write(*,*)  myName//'vco_anl%Vcode = ', vco_anl%Vcode
      call utl_abort( myName//': unknown vertical coordinate type!')
    end if

    numvar2d = 0

    allocate( bdiff_varNameList( vnl_numvarmax ) )
    bdiff_varNameList( : )=''
    allocate( nsposit( vnl_numvarmax + 1 ) )
    nsposit(1) = 1

    ! Find the 2D variables (within NAMSTATE namelist)

    if ( gsv_varExist( varName = 'GL  ' )) then

      numvar2d = numvar2d + 1
      nsposit( numvar2d + 1 ) = nsposit( numvar2d ) + 1
      bdiff_varNameList( numvar2d ) = 'GL  '

    end if
   
    if ( gsv_varExist( varName = 'TM  ')) then

      numvar2d = numvar2d + 1
      nsposit( numvar2d + 1 ) = nsposit( numvar2d ) + 1
      bdiff_varNameList( numvar2d ) = 'TM  '

    end if

    if ( numvar2d == 0) then
       
      if ( mpi_myid == 0) then
        write(*,*) myName//': Bdiff matrix not produced.'
        write(*,*) myName//': END'
      end if
      call utl_tmg_stop(65)
      cvdim_out = 0
      return
      
    else if (mpi_myid == 0) then

      write(*,*) myName//': number of 2D variables', numvar2d, bdiff_varNameList( 1 : numvar2d )

    end if

    ! default values for namelist variables
    corr_len(:) = 10.0
    stab(:)     = 0.2
    nsamp(:)    = 10000
    useImplicit(:) = .false.
    scaleFactor(:) = 0.0d0
    stddevMode  = 'GD2D'
    homogeneous_std(:) = -1.0d0

    nulnam = 0
    ierr = fnom( nulnam,'./flnml','FTN+SEQ+R/O',0)
    read( nulnam, nml = nambdiff, iostat = ierr )
    if ( ierr /= 0 ) call utl_abort( myName//': Error reading namelist')
    if ( mpi_myid == 0) write( *, nml = nambdiff )
    ierr = fclos( nulnam )

    if ( sum(scaleFactor(:) ) == 0.0d0 ) then
      if( mpi_myid == 0) write(*,*) myName//': scaleFactor=0, skipping rest of setup'
      cvdim_out = 0
      call utl_tmg_stop(65)
      return
    end if

    if ( trim(bdiff_mode) == 'BackgroundCheck' ) then
      cvDim_out = 9999 ! Dummy value > 0 to indicate to the background check (s/r ose_compute_HBHT_ensemble) 
                       ! that Diff is used
      call utl_tmg_stop(65)
      return
    end if

    ! Assumes the input 'scalefactor' is a scaling factor of the variances.

    do variableIndex = 1, numvar2d
      if ( scaleFactor( variableIndex ) > 0.0d0 ) then 
        scaleFactor_sigma( variableIndex ) = sqrt( scaleFactor( variableIndex ) )
      else
        scaleFactor_sigma( variableIndex ) = 0.0d0
      end if
    end do

    ni_l = hco_in%ni
    nj_l = hco_in%nj

    allocate( diffID( numvar2d ) )
    do variableIndex = 1, numvar2d
      write(*,*) myName//': setup the diffusion operator for the variable ', bdiff_varNameList( variableIndex ) 
      diffID( variableIndex ) = diff_setup ( variableIndex, bdiff_varNameList(1:numvar2d), hco_in, vco_in, corr_len( variableIndex ), &
                                             stab( variableIndex ), nsamp( variableIndex ), useImplicit( variableIndex ) )
    end do

    call mpivar_setup_latbands( nj_l, latPerPE, latPerPEmax, myLatBeg, myLatEnd )
    call mpivar_setup_lonbands( ni_l, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd )

    ! compute mpilocal control vector size
    cvDim_mpilocal = lonPerPE * latPerPE * numvar2d
    cvDim_out = cvDim_mpilocal

    ! also compute mpiglobal control vector dimension
    call rpn_comm_allreduce( cvDim_mpilocal, cvDim_mpiglobal, 1, "mpi_integer", "mpi_sum", "GRID", ierr )

    allocate(stddev(myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d))

    call bdiff_rdstats( hco_in, vco_in )

    if(mpi_myid == 0) write(*,*) myName//': Memory Used: ',get_max_rss()/1024,'Mb'

    if(mpi_myid == 0) write(*,*) myName//': END'

    initialized = .true.

    call utl_tmg_stop(65)

  end subroutine bdiff_setup

  
  subroutine bdiff_getScaleFactor( scaleFactor_out )
    
    implicit none

    real(8), intent(out) :: scaleFactor_out(:)

    integer :: variableIndex

    do variableIndex = 1, numvar2d
       
      scaleFactor_out( variableIndex ) = scaleFactor( variableIndex )
       
    end do

  end subroutine bdiff_getScaleFactor
  

  subroutine bdiff_rdstats( hco_in, vco_in )
    !
    !:Purpose: To read background-error stats file.
    !
    implicit none
    
    type(struct_hco), pointer :: hco_in
    type(struct_vco), pointer :: vco_in
    
    integer :: ierr, nmax, fnom, fstouv, fstfrm, fclos
    integer :: variableIndex
    logical :: lExists
    character(len=12) :: bFileName1 = './bgstddev'
    character(len=8)  :: bFileName2 = './bgcov'
    character(len=*), parameter :: myName = 'bdiff_rdstats'
    
    write(*,*) myName//': stddevMode is ', stddevMode
    write(*,*) myName//': Number of 2D variables', numvar2d, bdiff_varNameList( 1 : numvar2d )
    
    if ( stddevMode == 'GD2D' ) then

      inquire( file = bFileName1, exist = lExists )
       
      if ( lexists ) then
        ierr = fnom(nulbgst, bFileName1, 'RND+OLD+R/O', 0)
        if ( ierr == 0 ) then
          nmax = fstouv(nulbgst, 'RND+OLD')
        else
          call utl_abort( myName//': error opening file '//trim(bFileName1))
        end if
      else
        ! Assume background-error stats in file bgcov. 
        inquire( file = bFileName2, exist = lExists )  
        if ( lexists ) then 
          ierr = fnom( nulbgst, bFileName2, 'RND+OLD+R/O', 0 )
          if ( ierr == 0 ) then
            nmax = fstouv(nulbgst, 'RND+OLD')
          else
            call utl_abort( myName//': error opening file '//trim(bFileName2) )
          end if
        else
          call utl_abort( myName//': no background error statistics file found!!' )
        end if
      end if

      call bdiff_readBGstdField( hco_in, vco_in )

      ierr = fstfrm(nulbgst)
      ierr = fclos(nulbgst)

    else if ( stddevMode == 'HOMO' ) then

      do variableIndex = 1, numvar2d
        write(*,*) myName//': stdev = ', homogeneous_std( variableIndex ), ' for variable ', bdiff_varNameList( variableIndex )    
        stddev( :, :, variableIndex) = homogeneous_std( variableIndex )
      end do
       
    else

      call utl_abort( myName//': unknown stddevMode: '//trim(stddevMode) )

    end if

    call bdiff_scalestd

    write(*,*) myName//': END '
    
  end subroutine bdiff_rdstats

  subroutine bdiff_readBGstdField( hco_in, vco_in )
    !
    !:Purpose: to read 2D background error standard deviation field
    !          stored on Z, U or G grid and interpolate it to the analysis grid
    !
    implicit none
    
    ! Arguments
    type(struct_hco), pointer, intent(in) :: hco_in 
    type(struct_vco), pointer, intent(in) :: vco_in

    ! locals
    integer                     :: variableIndex, ierr
    type(struct_gsv)            :: statevector
    real(4), pointer            :: field3D_r4_ptr(:,:,:)
    real(8)                     :: minStddev, maxStddev
    character(len=*), parameter :: myName = 'bdiff_readBGstdField'

    write(*,*) myName//': Reading 2D fields from ./bgstddev...'
    write(*,*) myName//': Number of 2D variables', numvar2d, bdiff_varNameList( 1 : numvar2d )

    call gsv_allocate( statevector, 1, hco_in, vco_in, dateStamp_opt=-1, &
                       dataKind_opt=4, mpi_local_opt=.true., &
                       hInterpolateDegree_opt='LINEAR', &
                       varNames_opt=bdiff_varNameList(1:numvar2d) )
    call gsv_zero( statevector )
    call gsv_readFromFile(statevector, './bgstddev', 'STDDEV', ' ', unitConversion_opt = .false. )

    do variableIndex = 1, numvar2d
      call gsv_getField( statevector, field3D_r4_ptr, bdiff_varNameList( variableIndex ) )
      stddev( :, :, variableIndex ) = dble( field3D_r4_ptr( :, :, 1 ) )
      if (mpi_nprocs > 1) then
        call rpn_comm_allreduce(minval(stddev(:,:,variableIndex)),minStddev,1,'mpi_real8','mpi_min','GRID',ierr)
        call rpn_comm_allreduce(maxval(stddev(:,:,variableIndex)),maxStddev,1,'mpi_real8','mpi_max','GRID',ierr)
      else
        minStddev = minval(stddev(:,:,variableIndex))
        maxStddev = maxval(stddev(:,:,variableIndex))
      end if
      write(*,*) myName//': variable ', bdiff_varNameList( variableIndex ),' min/max: ', &
                 minStddev, maxStddev
    end do
   
    call gsv_deallocate( statevector )

  end subroutine bdiff_readBGstdField


  subroutine bdiff_scalestd
    !
    !:Purpose: To scale background-error standard-deviation values.
    !
    implicit none

    integer :: variableIndex
    character(len=*), parameter :: myName = 'bdiff_scalestd'
    
    do variableIndex = 1, numvar2d
       
      write(*,*) myName//': scaling ', bdiff_varNameList( variableIndex ), ' STD field with the factor ',  scaleFactor_sigma( variableIndex )
       
      stddev( :, :, variableIndex ) = scaleFactor_sigma( variableIndex ) * stddev( : , : , variableIndex )
      
    end do

  end subroutine bdiff_scalestd


  subroutine bdiff_bSqrt( controlVector_in, statevector )
    
    implicit none

    real(8),          intent(in)    :: controlVector_in( cvDim_mpilocal )
    type(struct_gsv), intent(inout) :: statevector

    real(8) :: gd_in( myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)
    real(8) :: gd_out(myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)

    integer :: variableIndex
    character(len=*), parameter :: myName = 'bdiff_bSqrt'
    
    if( .not. initialized) then
      if( mpi_myid == 0 ) write(*,*) myName//': bMatrixDIFF not initialized'
      return
    end if

    if(mpi_myid == 0) write(*,*) myName//': starting'

    call bdiff_cain( controlVector_in, gd_in )

    do variableIndex = 1, numvar2d

      ! Apply square root of the diffusion operator.
      call diff_Csqrt( diffID( variableIndex ), gd_in( :, :, variableIndex ), gd_out( :, :, variableIndex ) )

      ! Multiply by the diagonal matrix of background error standard deviations.
      gd_out( :, :, variableIndex ) = gd_out( :, :, variableIndex ) * stddev( :, :, variableIndex )

    end do

    call bdiff_copyToStatevector( statevector, gd_out )

    if(mpi_myid == 0) write(*,*) myName//': Memory Used: ',get_max_rss()/1024,'Mb'
    if(mpi_myid == 0) write(*,*) myName//': done'

  end subroutine bdiff_bSqrt


  subroutine bdiff_bSqrtAd( statevector, controlVector_out )
    
    implicit none

    type(struct_gsv), intent(in)  :: statevector
    real(8),          intent(out) :: controlVector_out(cvDim_mpilocal)

    real(8) :: gd_in( myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)
    real(8) :: gd_out(myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)

    integer :: variableIndex
    character(len=*), parameter :: myName = 'bdiff_bSqrtAd'

    if ( .not. initialized ) then
      if ( mpi_myid == 0 ) write(*,*) myName//': bMatrixDIFF not initialized'
      return
    end if

    if(mpi_myid == 0) write(*,*)  myName//': starting'

    call bdiff_copyFromStatevector( statevector, gd_in )

    do variableIndex = 1, numvar2d

      ! Multiply by the diagonal matrix of background-error standard deviations.
      gd_in( :, :, variableIndex ) = gd_in( :, :, variableIndex ) * stddev( :, :, variableIndex )

      ! Apply the adjoint of the square root of the diffusion operator.
      call diff_Csqrtadj( diffID( variableIndex ), gd_in( :, :, variableIndex ), gd_out( :, :, variableIndex) )

    end do

    call bdiff_cainad(gd_out, controlVector_out)

    if ( mpi_myid == 0) write(*,*) myName//': Memory Used: ', get_max_rss()/1024,'Mb'
    if ( mpi_myid == 0) write(*,*) myNAme//': done'

  end subroutine bdiff_bSqrtAd


  subroutine bdiff_copyToStatevector ( statevector, gd )
    
    implicit none

    real(8),          intent(in)    :: gd(myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)
    type(struct_gsv), intent(inout) :: statevector

    integer :: jlon, jlev, jlev2, jlat, variableIndex, ilev1, ilev2
    real(4), pointer :: field_r4(:,:,:)
    real(8), pointer :: field_r8(:,:,:)
    character(len=*), parameter :: myName = 'bdiff_copyToStatevector'
    
    do variableIndex = 1, numvar2d

      ilev1 = nsposit( variableIndex )
      ilev2 = nsposit( variableIndex + 1 ) - 1 

      if ( mpi_myid == 0) write(*,*) myName//': ',bdiff_varNameList( variableIndex )
      
      if (gsv_getDataKind(statevector) == 4) then
        call gsv_getField( statevector, field_r4, bdiff_varNameList( variableIndex ) )
        do jlev = ilev1, ilev2
          jlev2 = jlev-ilev1+1
          do jlat = myLatBeg, myLatEnd
            do jlon = myLonBeg, myLonEnd
              field_r4(jlon,jlat,jlev2) = gd(jlon,jlat,jlev)
            end do
          end do
        end do
      else
        call gsv_getField( statevector, field_r8, bdiff_varNameList( variableIndex ) )
        do jlev = ilev1, ilev2
          jlev2 = jlev-ilev1+1
          do jlat = myLatBeg, myLatEnd
            do jlon = myLonBeg, myLonEnd
              field_r8(jlon,jlat,jlev2) = gd(jlon,jlat,jlev)
            end do
          end do
        end do
      end if

    end do

  end subroutine bdiff_copyToStatevector


  subroutine bdiff_copyFromStatevector( statevector, gd )
    
    implicit none

    type(struct_gsv), intent(in)  :: statevector
    real(8),          intent(out) :: gd(myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)

    integer :: jlon, jlev, jlev2, jlat, variableIndex, ilev1, ilev2
    real(4), pointer :: field_r4(:,:,:)
    real(8), pointer :: field_r8(:,:,:)
    character(len=*), parameter :: myName = 'bdiff_copyFromStatevector'

    do variableIndex = 1, numvar2d

      ilev1 = nsposit( variableIndex )
      ilev2 = nsposit( variableIndex + 1 ) - 1 

      if ( mpi_myid == 0) write(*,*) myName//': ',bdiff_varNameList( variableIndex )
      if (gsv_getDataKind(statevector) == 4) then
        call gsv_getField(statevector, field_r4, bdiff_varNameList( variableIndex ))
        do jlev = ilev1, ilev2
          jlev2 = jlev-ilev1+1
          do jlat = myLatBeg, myLatEnd
            do jlon = myLonBeg, myLonEnd
              gd( jlon, jlat, jlev ) = field_r4( jlon, jlat, jlev2 )
            end do
          end do
        end do
      else
        call gsv_getField(statevector, field_r8, bdiff_varNameList( variableIndex ))
        do jlev = ilev1, ilev2
          jlev2 = jlev-ilev1+1
          do jlat = myLatBeg, myLatEnd
            do jlon = myLonBeg, myLonEnd
              gd( jlon, jlat, jlev ) = field_r8( jlon, jlat, jlev2 )
            end do
          end do
        end do
      end if
    end do

  end subroutine bdiff_copyFromStatevector


  subroutine bdiff_cain( controlVector_in, gd_out )
    
    implicit none

    real(8), intent(in)  :: controlVector_in(cvDim_mpilocal)
    real(8), intent(out) :: gd_out(myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)

    integer :: jn, jlev, jlon, jlat

    jn = 0
    do jlev = 1, numvar2d
      do jlat = myLatBeg, myLatEnd
        do jlon = myLonBeg, myLonEnd
          jn = jn + 1
          gd_out( jlon, jlat, jlev ) = ControlVector_in( jn )
        end do
      end do
    end do

  end subroutine bdiff_cain


  subroutine bdiff_cainAd( gd_in, diffControlVector_out )
    
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

  end subroutine bdiff_cainAd


  subroutine bdiff_Finalize()
    
    implicit none

    if ( initialized ) then
      initialized = .false.
      deallocate( stddev )
      deallocate( diffID )
      deallocate( nsposit )
      deallocate( bdiff_varNameList )
    end if

  end subroutine bdiff_Finalize


end module BmatrixDiff_mod
