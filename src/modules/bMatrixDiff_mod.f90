
MODULE bMatrixDiff_mod
  ! MODULE bMatrixDiff_mod (prefix='bdiff' category='2. B and R matrices')
  !
  ! :Purpose: Performs transformation from control vector to analysis increment 
  !           using the background-error covariance matrix based on correlations
  !           modelled using a diffusion operator.
  !
  use midasMpi_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use varNameList_mod
  use physicsFunctions_mod
  use utilities_mod
  use diffusion_mod
  implicit none
  save
  private

  ! public procedures
  public :: bdiff_Setup, bdiff_BSqrt, bdiff_BSqrtAd, bdiff_Finalize
  public :: bdiff_getScaleFactor, bdiff_reduceToMPILocal

  logical             :: initialized = .false.
  integer             :: nj_l, ni_l
  integer             :: cvDim_mpilocal, cvDim_mpiglobal

  integer, allocatable :: diffID(:)

  ! Background-error covariance matrix elements.
  real(8), allocatable :: stddev(:,:,:)

  integer, parameter  :: maxNumVars = 200
  real(8)             :: scaleFactor_sigma(maxNumVars)

  ! read in from the namelist:
  real(8)          :: scaleFactor(maxNumVars)     ! scale factor applied to variances
  character(len=4) :: stddevMode                  ! can be 'GD2D' or 'HOMO'
  real(8)          :: homogeneous_std(maxNumVars) ! homogeneous standard deviation (when stddevMode is 'HOMO')

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

  !--------------------------------------------------------------------------
  ! bdiff_setup
  !--------------------------------------------------------------------------
  subroutine bdiff_setup ( hco_in, vco_in, cvDim_out, mode_opt )
    !
    !:Purpose: Setup the diffusion B matrix
    !
    implicit none

    ! Arguments:
    type(struct_hco), pointer,  intent(in)  :: hco_in
    type(struct_vco), pointer,  intent(in)  :: vco_in
    integer         ,           intent(out) :: cvDim_out
    character(len=*), optional, intent(in)  :: mode_opt

    ! Locals:
    character(len=15)         :: bdiff_mode
    type(struct_vco), pointer :: vco_anl
    integer                   :: nulnam, ierr, fnom, fclos
    integer                   :: variableIndex, latIndex, latIndexIgnore
    real(8)                   :: maxDistance
    real(8), allocatable      :: distance(:)
    character(len=*), parameter :: myName = 'bdiff_setup'
    ! namelist variables
    real    :: corr_len( maxNumVars )  ! Horizontal correlation length scale (km)
    real    :: stab( maxNumVars )      ! Stability criteria (definitely < 0.5)
    integer :: nsamp(maxNumVars)       ! Number of samples in the estimation of the normalization factors by randomization
    real(8) :: latIgnoreFraction       ! Relative zonal grid spacing limit where lats near each numerical pole are ignored
    logical :: useImplicit(maxNumVars) ! choose to use implicit formulation of diffusion operator (.true.) or explicit version (.false.)
    
    NAMELIST /NAMBDIFF/ corr_len, stab, nsamp, useImplicit, scaleFactor, stddevMode, homogeneous_std, latIgnoreFraction

    call utl_tmg_start(65,'----B_DIFF_Setup')
    if(mmpi_myid == 0) write(*,*) myName//': starting'
    if(mmpi_myid == 0) write(*,*) myName//': Memory Used: ',get_max_rss()/1024,'Mb'

    ! default values for namelist variables
    corr_len(:) = 10.0
    stab(:)     = 0.2
    nsamp(:)    = 10000
    useImplicit(:) = .false.
    scaleFactor(:) = 0.0d0
    stddevMode  = 'GD2D'
    homogeneous_std(:) = -1.0d0
    latIgnoreFraction = 1.0d6 ! large value so that nothing is ignored by default

    if ( .not. utl_isNamelistPresent('NAMBDIFF','./flnml') ) then
      if ( mmpi_myid == 0 ) then
        write(*,*) 'bdiff_setup: nambdiff is missing in the namelist.'
        write(*,*) '             The default values will be taken.'
      end if
    else
      nulnam = 0
      ierr = fnom( nulnam,'./flnml','FTN+SEQ+R/O',0)
      read( nulnam, nml = nambdiff, iostat = ierr )
      if ( ierr /= 0 ) call utl_abort( myName//': Error reading namelist')
      if ( mmpi_myid == 0) write( *, nml = nambdiff )
      ierr = fclos( nulnam )
    end if

    if ( sum(scaleFactor(:) ) == 0.0d0 ) then
      if( mmpi_myid == 0) write(*,*) myName//': scaleFactor=0, skipping rest of setup'
      cvdim_out = 0
      call utl_tmg_stop(65)
      return
    end if

    if ( present( mode_opt ) ) then
      if ( trim( mode_opt ) == 'Analysis' .or. trim( mode_opt ) == 'BackgroundCheck' ) then
        bdiff_mode = trim( mode_opt )
        if( mmpi_myid == 0 ) write(*,*)
        if( mmpi_myid == 0 ) write(*,*) myName//': Mode activated = ', trim(bdiff_mode)
      else
        write(*,*)
        write(*,*)  myName//'mode = ', trim(mode_opt)
        call utl_abort( myName//': unknown mode' )
      end if
    else
      bdiff_mode = 'Analysis'
      if( mmpi_myid == 0 ) write(*,*)
      if( mmpi_myid == 0 ) write(*,*) myName//': analysis mode activated (by default)'
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
       
      if ( mmpi_myid == 0) then
        write(*,*) myName//': Bdiff matrix not produced.'
        write(*,*) myName//': END'
      end if
      call utl_tmg_stop(65)
      cvdim_out = 0
      return
      
    else if (mmpi_myid == 0) then

      write(*,*) myName//': number of 2D variables', numvar2d, bdiff_varNameList( 1 : numvar2d )

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

    ! Compute latIndexIgnore from latIgnoreFraction
    if (latIgnoreFraction < 1.0) then
      allocate(distance(nj_l))
      do latIndex = 1, nj_l
        distance(latIndex) = &
             phf_calcDistance(real(hco_in%lat2d_4(ni_l/2,latIndex),8), real(hco_in%lon2d_4((ni_l/2)+1,latIndex),8), &
                              real(hco_in%lat2d_4(ni_l/2,latIndex),8), real(hco_in%lon2d_4((ni_l/2)  ,latIndex),8))
      end do
      maxDistance = maxval(distance(:))
      if (mmpi_myid==0) write(*,*) myName//': maxDistance = ', maxDistance
      latIndexIgnore = 0
      latLoop: do latIndex = 1, nj_l
        if (distance(latIndex)/maxDistance > latIgnoreFraction) then
          if (mmpi_myid==0) then
            write(*,*) '   latIndex-1, distance, fraction = ', latIndex-1, distance(latIndex-1), distance(latIndex-1)/maxDistance 
            write(*,*) '***latIndex  , distance, fraction = ', latIndex  , distance(latIndex  ), distance(latIndex  )/maxDistance 
            write(*,*) '   latIndex+1, distance, fraction = ', latIndex+1, distance(latIndex+1), distance(latIndex+1)/maxDistance
          end if
          latIndexIgnore = latIndex
          exit latLoop
        end if
      end do latLoop
      deallocate(distance)
    else
      latIndexIgnore = 0
    end if
    if (mmpi_myid==0) write(*,*) myName//': latIndexIgnore = ', latIndexIgnore

    allocate( diffID( numvar2d ) )
    do variableIndex = 1, numvar2d
      write(*,*) myName//': setup the diffusion operator for the variable ', bdiff_varNameList( variableIndex ) 
      diffID( variableIndex ) = diff_setup ( variableIndex, bdiff_varNameList(1:numvar2d), hco_in, vco_in, corr_len( variableIndex ), &
                                             stab( variableIndex ), nsamp( variableIndex ), useImplicit( variableIndex ), &
                                             latIndexIgnore )
    end do

    call mmpi_setup_latbands( nj_l, latPerPE, latPerPEmax, myLatBeg, myLatEnd )
    call mmpi_setup_lonbands( ni_l, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd )

    ! compute mpilocal control vector size
    cvDim_mpilocal = lonPerPE * latPerPE * numvar2d
    cvDim_out = cvDim_mpilocal

    ! also compute mpiglobal control vector dimension
    call rpn_comm_allreduce( cvDim_mpilocal, cvDim_mpiglobal, 1, "mpi_integer", "mpi_sum", "GRID", ierr )

    allocate(stddev(myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d))

    call bdiff_rdstats( hco_in, vco_in )

    if(mmpi_myid == 0) write(*,*) myName//': Memory Used: ',get_max_rss()/1024,'Mb'

    if(mmpi_myid == 0) write(*,*) myName//': END'

    initialized = .true.

    call utl_tmg_stop(65)

  end subroutine bdiff_setup
  
  !--------------------------------------------------------------------------
  ! bdiff_getScaleFactor
  !--------------------------------------------------------------------------
  subroutine bdiff_getScaleFactor( scaleFactor_out )
    !
    !:Purpose: Return the specified scaleFactor.
    !    
    implicit none

    ! Arguments:
    real(8), intent(out) :: scaleFactor_out(:)

    ! Locals:
    integer :: variableIndex

    do variableIndex = 1, numvar2d
       
      scaleFactor_out( variableIndex ) = scaleFactor( variableIndex )
       
    end do

  end subroutine bdiff_getScaleFactor
  
  !--------------------------------------------------------------------------
  ! bdiff_rdstats
  !--------------------------------------------------------------------------
  subroutine bdiff_rdstats( hco_in, vco_in )
    !
    !:Purpose: To read background-error stats file.
    !
    implicit none
    
    ! Arguments:
    type(struct_hco), pointer, intent(in) :: hco_in
    type(struct_vco), pointer, intent(in) :: vco_in
    
    ! Locals:
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

  !--------------------------------------------------------------------------
  ! bdiff_readBGstdField
  !--------------------------------------------------------------------------
  subroutine bdiff_readBGstdField( hco_in, vco_in )
    !
    !:Purpose: to read 2D background error standard deviation field
    !          stored on Z, U or G grid and interpolate it to the analysis grid
    !
    implicit none
    
    ! Arguments:
    type(struct_hco), pointer, intent(in) :: hco_in 
    type(struct_vco), pointer, intent(in) :: vco_in

    ! Locals:
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
    call gio_readFromFile(statevector, './bgstddev', 'STDDEV', ' ', unitConversion_opt = .false. )

    do variableIndex = 1, numvar2d
      call gsv_getField( statevector, field3D_r4_ptr, bdiff_varNameList( variableIndex ) )
      stddev( :, :, variableIndex ) = dble( field3D_r4_ptr( :, :, 1 ) )
      if (mmpi_nprocs > 1) then
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

  !--------------------------------------------------------------------------
  ! bdiff_scalestd
  !--------------------------------------------------------------------------
  subroutine bdiff_scalestd
    !
    !:Purpose: To scale background-error standard-deviation values.
    !
    implicit none

    ! Locals:
    integer :: variableIndex
    character(len=*), parameter :: myName = 'bdiff_scalestd'
    
    do variableIndex = 1, numvar2d
       
      write(*,*) myName//': scaling ', bdiff_varNameList( variableIndex ), ' STD field with the factor ',  scaleFactor_sigma( variableIndex )
       
      stddev( :, :, variableIndex ) = scaleFactor_sigma( variableIndex ) * stddev( : , : , variableIndex )
      
    end do

  end subroutine bdiff_scalestd

  !--------------------------------------------------------------------------
  ! bdiff_bSqrt
  !--------------------------------------------------------------------------
  subroutine bdiff_bSqrt( controlVector_in, statevector )
    !
    !:Purpose: Apply the sqrt of the B matrix 
    !    
    implicit none

    ! Arguments:
    real(8),          intent(in)    :: controlVector_in( cvDim_mpilocal )
    type(struct_gsv), intent(inout) :: statevector

    ! Locals:
    real(8) :: gd_in( myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)
    real(8) :: gd_out(myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)
    integer :: variableIndex
    character(len=*), parameter :: myName = 'bdiff_bSqrt'
    
    if( .not. initialized) then
      if( mmpi_myid == 0 ) write(*,*) myName//': bMatrixDIFF not initialized'
      return
    end if

    if(mmpi_myid == 0) write(*,*) myName//': starting'

    call bdiff_cain( controlVector_in, gd_in )

    do variableIndex = 1, numvar2d

      ! Apply square root of the diffusion operator.
      call diff_Csqrt( diffID( variableIndex ), gd_in( :, :, variableIndex ), gd_out( :, :, variableIndex ) )

      ! Multiply by the diagonal matrix of background error standard deviations.
      gd_out( :, :, variableIndex ) = gd_out( :, :, variableIndex ) * stddev( :, :, variableIndex )

    end do

    call bdiff_copyToStatevector( statevector, gd_out )

    if(mmpi_myid == 0) write(*,*) myName//': Memory Used: ',get_max_rss()/1024,'Mb'
    if(mmpi_myid == 0) write(*,*) myName//': done'

  end subroutine bdiff_bSqrt

  !--------------------------------------------------------------------------
  ! bdiff_bSqrtAd
  !--------------------------------------------------------------------------
  subroutine bdiff_bSqrtAd( statevector, controlVector_out )
    !
    !:Purpose: Apply the adjoint (i.e. transpose) of the sqrt of the B matrix
    !    
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)  :: statevector
    real(8),          intent(out) :: controlVector_out(cvDim_mpilocal)

    ! Locals:
    real(8) :: gd_in( myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)
    real(8) :: gd_out(myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)
    integer :: variableIndex
    character(len=*), parameter :: myName = 'bdiff_bSqrtAd'

    if ( .not. initialized ) then
      if ( mmpi_myid == 0 ) write(*,*) myName//': bMatrixDIFF not initialized'
      return
    end if

    if(mmpi_myid == 0) write(*,*)  myName//': starting'

    call bdiff_copyFromStatevector( statevector, gd_in )

    do variableIndex = 1, numvar2d

      ! Multiply by the diagonal matrix of background-error standard deviations.
      gd_in( :, :, variableIndex ) = gd_in( :, :, variableIndex ) * stddev( :, :, variableIndex )

      ! Apply the adjoint of the square root of the diffusion operator.
      call diff_Csqrtadj( diffID( variableIndex ), gd_in( :, :, variableIndex ), gd_out( :, :, variableIndex) )

    end do

    call bdiff_cainad(gd_out, controlVector_out)

    if ( mmpi_myid == 0) write(*,*) myName//': Memory Used: ', get_max_rss()/1024,'Mb'
    if ( mmpi_myid == 0) write(*,*) myNAme//': done'

  end subroutine bdiff_bSqrtAd

  !--------------------------------------------------------------------------
  ! bdiff_copyToStatevector
  !--------------------------------------------------------------------------
  subroutine bdiff_copyToStatevector ( statevector, gd )
    !
    !:Purpose: Copy the working array to a statevector object 
    !    
    implicit none

    ! Arguments:
    real(8),          intent(in)    :: gd(myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)
    type(struct_gsv), intent(inout) :: statevector

    ! Locals:
    integer :: jlon, jlev, jlev2, jlat, variableIndex, ilev1, ilev2
    real(4), pointer :: field_r4(:,:,:)
    real(8), pointer :: field_r8(:,:,:)
    character(len=*), parameter :: myName = 'bdiff_copyToStatevector'
    
    do variableIndex = 1, numvar2d

      ilev1 = nsposit( variableIndex )
      ilev2 = nsposit( variableIndex + 1 ) - 1 

      if ( mmpi_myid == 0) write(*,*) myName//': ',bdiff_varNameList( variableIndex )
      
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

  !--------------------------------------------------------------------------
  ! bdiff_copyFromStatevector
  !--------------------------------------------------------------------------
  subroutine bdiff_copyFromStatevector( statevector, gd )
    !
    !:Purpose: Copy the contents of the statevector object to the working array
    !    
    implicit none

    ! Arguments:
    type(struct_gsv), intent(in)  :: statevector
    real(8),          intent(out) :: gd(myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)

    ! Locals:
    integer :: jlon, jlev, jlev2, jlat, variableIndex, ilev1, ilev2
    real(4), pointer :: field_r4(:,:,:)
    real(8), pointer :: field_r8(:,:,:)
    character(len=*), parameter :: myName = 'bdiff_copyFromStatevector'

    do variableIndex = 1, numvar2d

      ilev1 = nsposit( variableIndex )
      ilev2 = nsposit( variableIndex + 1 ) - 1 

      if ( mmpi_myid == 0) write(*,*) myName//': ',bdiff_varNameList( variableIndex )
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

  !--------------------------------------------------------------------------
  ! bdiff_reduceToMPILocal
  !--------------------------------------------------------------------------
  subroutine bdiff_reduceToMPILocal(cv_mpilocal,cv_mpiglobal)
    !
    !:Purpose: Extract the subset of the global control vector needed for local MPI task
    !
    implicit none

    ! Arguments:
    real(8), intent(out) :: cv_mpilocal(:)
    real(8), intent(in)  :: cv_mpiglobal(:)

    ! Locals:
    integer :: jn, jlat, jlon, jlev, ierr
    real(8), allocatable :: gd_mpiGlobal(:,:,:)

    allocate(gd_mpiGlobal(ni_l,nj_l,numvar2d))
    gd_mpiGlobal(:,:,:) = 0.0d0

    jn = 0
    if (mmpi_myid == 0) then
      do jlev = 1, numvar2d
        do jlat = 1, nj_l
          do jlon = 1, ni_l
            jn = jn + 1
            gd_mpiGlobal( jlon, jlat, jlev ) = cv_mpiglobal( jn )
          end do
        end do
      end do
    end if 
    call rpn_comm_bcast(gd_mpiGlobal, size(gd_mpiGlobal), 'MPI_REAL8', 0, 'GRID', ierr)
    
    jn = 0
    do jlev = 1, numvar2d
      do jlat = myLatBeg, myLatEnd
        do jlon = myLonBeg, myLonEnd
          jn = jn + 1
          cv_mpilocal(jn) = gd_mpiGlobal(jlon,jlat,jlev)
        end do
      end do
    end do

    deallocate(gd_mpiGlobal)

  end subroutine bdiff_reduceToMPILocal

  !--------------------------------------------------------------------------
  ! bdiff_cain
  !--------------------------------------------------------------------------
  subroutine bdiff_cain( controlVector_in, gd_out )
    !
    !:Purpose: Transform from control vector to working array
    !    
    implicit none

    ! Arguments:
    real(8), intent(in)  :: controlVector_in(cvDim_mpilocal)
    real(8), intent(out) :: gd_out(myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)

    ! Locals:
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

  !--------------------------------------------------------------------------
  ! bdiff_cainAd
  !--------------------------------------------------------------------------
  subroutine bdiff_cainAd( gd_in, diffControlVector_out )
    !
    !:Purpose: Transform from working array to control vector
    !        
    implicit none

    ! Arguments:
    real(8), intent(in)  :: gd_in(myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d)
    real(8), intent(out) :: diffControlVector_out(cvDim_mpilocal)

    ! Locals:
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

  !--------------------------------------------------------------------------
  ! bdiff_Finalize
  !--------------------------------------------------------------------------
  subroutine bdiff_Finalize()
    !
    !:Purpose: Deallocate some arrays after we don't need the B matrix anymore. 
    !        
    implicit none

    if ( initialized ) then
      initialized = .false.
      deallocate( stddev )
      deallocate( diffID )
      deallocate( nsposit )
      deallocate( bdiff_varNameList )
    end if

  end subroutine bdiff_Finalize

end module bMatrixDiff_mod
