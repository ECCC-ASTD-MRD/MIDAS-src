
MODULE bMatrixDiff_mod
  ! MODULE bMatrixDiff_mod (prefix='bdiff' category='2. B and R matrices')
  !
  !:Purpose:  Performs transformation from control vector to analysis increment
  !           (and the adjoint transformation) using the background-error
  !           covariance matrix based on correlations modelled using a
  !           diffusion operator.
  !
  use rmn_fst98
  use midasMpi_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use varNameList_mod
  use physicsFunctions_mod
  use utilities_mod
  use runtimeInfo_mod
  use diffusion_mod
  use message_mod
  use timeCoord_mod

  implicit none
  save
  private

  ! Public procedures
  public :: bdiff_Setup, bdiff_BSqrt, bdiff_BSqrtAd, bdiff_Finalize
  public :: bdiff_getScaleFactor, bdiff_reduceToMPILocal
  public :: bdiff_getSSTBGstdFromFourSeasons

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
  logical          :: fourSeasonsBgstdSST         ! Compute daily BG stddev from 4 seasonal fields valid on the 15th of Feb, May, Aug and Nov

  ! Number of incremental variables/fields
  integer             :: numvar2d
  ! Start position of each field in composite arrays
  integer, allocatable :: nsposit(:)
  ! Name list of incremental variables/fields
  character(len=4), allocatable :: bdiff_varNameList(:)

  integer             :: myLatBeg, myLatEnd
  integer             :: myLonBeg, myLonEnd
  integer             :: latPerPE, latPerPEmax, lonPerPE, lonPerPEmax

  integer             :: nulbgst = 0

CONTAINS

  !--------------------------------------------------------------------------
  ! bdiff_setup
  !--------------------------------------------------------------------------
  subroutine bdiff_setup (hco_in, vco_in, cvDim_out, mode_opt)
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
    integer                   :: ierr
    integer                   :: variableIndex, latIndex, latIndexIgnore
    real(8)                   :: maxDistance
    real(8), allocatable      :: distance(:)

    ! namelist variables
    real    :: corr_len( maxNumVars )  ! Horizontal correlation length scale (km)
    real    :: stab( maxNumVars )      ! Stability criteria (definitely < 0.5)
    integer :: nsamp(maxNumVars)       ! Number of samples in the estimation of the normalization factors by randomization
    real(8) :: latIgnoreFraction       ! Relative zonal grid spacing limit where lats near each numerical pole are ignored
    logical :: useImplicit(maxNumVars) ! choose to use implicit formulation of diffusion operator (.true.) or explicit version (.false.)

    NAMELIST /NAMBDIFF/ corr_len, stab, nsamp, useImplicit, scaleFactor, stddevMode, &
                        homogeneous_std, latIgnoreFraction, fourSeasonsBgstdSST

    call rti_tmg_start(65,'----B_DIFF_Setup')
    if(mmpi_myid == 0) call msg('bdiff_setup', 'Starting')
    call msg_memUsage('bdiff_setup', mpiAll_opt=.false.)

    ! default values for namelist variables
    corr_len(:) = 10.0
    stab(:)     = 0.2
    nsamp(:)    = 10000
    useImplicit(:) = .false.
    scaleFactor(:) = 0.0d0
    stddevMode  = 'GD2D'
    homogeneous_std(:) = -1.0d0
    latIgnoreFraction = 1.0d6 ! large value so that nothing is ignored by default
    fourSeasonsBgstdSST = .false.

    if (.not. utl_isNamelistPresent('NAMBDIFF','./flnml')) then
      if ( mmpi_myid == 0 ) then
        call msg('bdiff_setup', 'nambdiff is missing in the namelist.')
        call msg('bdiff_setup', 'The default values will be taken.')
      end if
    else
      call rti_tmg_start(181,'low-level--readNML')
      read(utl_flnml, nml = nambdiff, iostat = ierr)
      if (ierr /= 0) call rti_abort('bdiff_setup: Error reading namelist')
      if (mmpi_myid == 0) write(*, nml = nambdiff)
      call rti_tmg_stop(181)
    end if

    if ( utl_isEqual(sum(scaleFactor(:)),0.0d0) ) then
      if(mmpi_myid == 0) call msg('bdiff_setup', 'scaleFactor=0, skipping rest of setup')
      cvdim_out = 0
      call rti_tmg_stop(65)
      return
    end if

    if (present(mode_opt)) then
      if (trim(mode_opt) == 'Analysis' .or. trim( mode_opt) == 'BackgroundCheck') then
        bdiff_mode = trim( mode_opt )
        if(mmpi_myid == 0) write(*,*)
        if(mmpi_myid == 0) call msg('bdiff_setup', 'Mode activated = '//trim(bdiff_mode))
      else
        write(*,*)
        call msg('bdiff_setup', 'mode = '//trim(mode_opt))
        call rti_abort('bdiff_setup: unknown mode')
      end if
    else
      bdiff_mode = 'Analysis'
      if(mmpi_myid == 0) write(*,*)
      if(mmpi_myid == 0) call msg('bdiff_setup', 'Analysis mode activated (by default)')
    end if

    vco_anl => vco_in
    if (vco_anl%Vcode /= 5002 .and. vco_anl%Vcode /= 5005 .and. vco_anl%Vcode /= 0) then
      call msg('bdiff_setup', 'vco_anl%Vcode = '//str(vco_anl%Vcode))
      call rti_abort('bdiff_setup: unknown vertical coordinate type!')
    end if

    numvar2d = 0

    allocate(bdiff_varNameList(vnl_numvarmax))
    bdiff_varNameList(:)=''
    allocate(nsposit(vnl_numvarmax + 1))
    nsposit(1) = 1

    ! Find the 2D variables (within NAMSTATE namelist)

    if (gsv_varExist(varName = 'GL  ')) then

      numvar2d = numvar2d + 1
      nsposit(numvar2d + 1) = nsposit(numvar2d) + 1
      bdiff_varNameList(numvar2d) = 'GL  '

    end if

    if (gsv_varExist(varName = 'TM  ')) then

      numvar2d = numvar2d + 1
      nsposit( numvar2d + 1 ) = nsposit( numvar2d ) + 1
      bdiff_varNameList(numvar2d) = 'TM  '

    end if

    if (gsv_varExist(varName = 'GE  ')) then

      numvar2d = numvar2d + 1
      nsposit(numvar2d + 1) = nsposit(numvar2d) + 1
      bdiff_varNameList(numvar2d) = 'GE  '

    end if

    if (numvar2d == 0) then

      if (mmpi_myid == 0) then
        call msg('bdiff_setup', 'Bdiff matrix not produced.')
        call msg('bdiff_setup', 'END')
      end if
      call rti_tmg_stop(65)
      cvdim_out = 0
      return

    else if (mmpi_myid == 0) then
      write(*,*) 'bdiff_setup: number of 2D variables', numvar2d, bdiff_varNameList(1:numvar2d)
    end if

    if (trim(bdiff_mode) == 'BackgroundCheck') then
      cvDim_out = 9999 ! Dummy value > 0 to indicate to the background check (s/r ose_compute_HBHT_ensemble)
                       ! that Diff is used
      call rti_tmg_stop(65)
      return
    end if

    ! Assumes the input 'scalefactor' is a scaling factor of the variances.

    do variableIndex = 1, numvar2d
      if (scaleFactor( variableIndex ) > 0.0d0) then
        scaleFactor_sigma(variableIndex) = sqrt(scaleFactor(variableIndex))
      else
        scaleFactor_sigma(variableIndex) = 0.0d0
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
      if (mmpi_myid==0) write(*,*) 'bdiff_setup: maxDistance = ', maxDistance
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
    if (mmpi_myid==0) write(*,*) 'bdiff_setup: latIndexIgnore = ', latIndexIgnore

    allocate( diffID( numvar2d ) )
    do variableIndex = 1, numvar2d
      write(*,*) 'bdiff_setup: setup the diffusion operator for the variable ', &
                 bdiff_varNameList(variableIndex)
      diffID(variableIndex) = diff_setup (variableIndex, bdiff_varNameList(1:numvar2d), hco_in, vco_in, &
                                          corr_len(variableIndex), stab(variableIndex), nsamp(variableIndex), &
                                          useImplicit(variableIndex), latIndexIgnore)
    end do

    call mmpi_setup_latbands(nj_l, latPerPE, latPerPEmax, myLatBeg, myLatEnd)
    call mmpi_setup_lonbands(ni_l, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd)

    ! compute mpilocal control vector size
    cvDim_mpilocal = lonPerPE * latPerPE * numvar2d
    cvDim_out = cvDim_mpilocal

    ! also compute mpiglobal control vector dimension
    call mmpi_allReduce(cvDim_mpilocal, cvDim_mpiglobal, mmpi_sum)

    allocate(stddev(myLonBeg:myLonEnd, myLatBeg:myLatEnd, numvar2d))

    call bdiff_rdstats(hco_in, vco_in)

    call msg_memUsage('bdiff_setup', mpiAll_opt=.false.)

    if(mmpi_myid == 0) call msg('bdiff_setup', 'Completed')

    initialized = .true.

    call rti_tmg_stop(65)

  end subroutine bdiff_setup

  !--------------------------------------------------------------------------
  ! bdiff_getScaleFactor
  !--------------------------------------------------------------------------
  subroutine bdiff_getScaleFactor(scaleFactor_out)
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
  subroutine bdiff_rdstats(hco_in, vco_in)
    !
    !:Purpose: To read background-error stats file.
    !
    implicit none

    ! Arguments:
    type(struct_hco), pointer, intent(in) :: hco_in
    type(struct_vco), pointer, intent(in) :: vco_in

    ! Locals:
    integer :: ierr, nmax
    integer :: variableIndex
    logical :: lExists

    call msg('bdiff_rdstats', 'stddevMode is '//stddevMode)
    write(*,*) 'bdiff_rdstats: Number of 2D variables', numvar2d, bdiff_varNameList(1:numvar2d)

    if (stddevMode == 'GD2D') then

      inquire(file = './bgstddev', exist = lExists)

      if (lexists) then
        ierr = fnom(nulbgst, './bgstddev', 'RND+OLD+R/O', 0)
        if (ierr == 0) then
          nmax = fstouv(nulbgst, 'RND+OLD')
        else
          call rti_abort('bdiff_rdstats: Error opening file bgstddev')
        end if
      else
        ! Assume background-error stats in file bgcov.
        inquire( file = './bgcov', exist = lExists )
        if (lexists) then
          ierr = fnom(nulbgst, './bgcov', 'RND+OLD+R/O', 0)
          if (ierr == 0) then
            nmax = fstouv(nulbgst, 'RND+OLD')
          else
            call rti_abort('bdiff_rdstats: error opening file bgcov')
          end if
        else
          call rti_abort('bdiff_rdstats: No background error statistics file found.')
        end if
      end if

      if (fourSeasonsBgstdSST) then
        call bdiff_getSSTBGstdFromFourSeasons(hco_in, vco_in)
      else
        call bdiff_readBGstdField(hco_in, vco_in)
      end if

      ierr = fstfrm(nulbgst)
      ierr = fclos(nulbgst)

    else if (stddevMode == 'HOMO') then

      do variableIndex = 1, numvar2d
        write(*,*) 'bdiff_rdstats: stdev = ', homogeneous_std(variableIndex), &
                   ' for variable ', bdiff_varNameList(variableIndex)
        stddev(:,:,variableIndex) = homogeneous_std(variableIndex)
      end do

    else

      call rti_abort('bdiff_rdstats: unknown stddevMode: '//trim(stddevMode))

    end if

    call bdiff_scalestd

    call msg('bdiff_rdstats', 'Completed.')

  end subroutine bdiff_rdstats

  !--------------------------------------------------------------------------
  ! bdiff_readBGstdField
  !--------------------------------------------------------------------------
  subroutine bdiff_readBGstdField(hco, vco)
    !
    !:Purpose: to read 2D background error standard deviation field
    !          stored on Z, U or G grid and interpolate it to the analysis grid
    !
    implicit none

    ! Arguments:
    type(struct_hco), pointer, intent(in) :: hco
    type(struct_vco), pointer, intent(in) :: vco

    ! Locals:
    integer          :: variableIndex
    type(struct_gsv) :: stateVector
    real(8), pointer :: field3D_r8_ptr(:,:,:)
    real(8)          :: minStddev, maxStddev

    call msg('bdiff_readBGstdField', 'Reading 2D fields from ./bgstddev...')
    call msg('bdiff_readBGstdField', 'Number of 2D variables'//str(numvar2d)//&
                                     str(bdiff_varNameList(1:numvar2d)))

    call gsv_allocate(stateVector, 1, hco, vco, dateStamp_opt = -1, &
                      dataKind_opt = 8, mpi_local_opt = .true., &
                      hInterpolateDegree_opt = 'LINEAR', &
                      varNames_opt = bdiff_varNameList(1:numvar2d))
    call gsv_zero(stateVector)

    call gio_readFromFile(stateVector, './bgstddev', 'STDDEV', ' ', unitConversion_opt = .false.)

    do variableIndex = 1, numvar2d

      call gsv_getField(statevector, field3D_r8_ptr, bdiff_varNameList(variableIndex))
      stddev(:,:,variableIndex) = dble(field3D_r8_ptr(:,:,1))
      if (mmpi_nprocs > 1) then
        call mmpi_allreduce(minval(stddev(:,:,variableIndex)), minStddev, mmpi_min)
        call mmpi_allreduce(maxval(stddev(:,:,variableIndex)), maxStddev, mmpi_max)
      else
        minStddev = minval(stddev(:,:,variableIndex))
        maxStddev = maxval(stddev(:,:,variableIndex))
      end if

      call msg('bdiff_readBGstdField', 'Variable '//bdiff_varNameList(variableIndex)//&
                                       ' min/max: '//str(minStddev)//', '//str(maxStddev))
    end do

    call gsv_deallocate(statevector)

    call msg('bdiff_readBGstdField', 'Completed.')

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

    call msg_memUsage('bdiff_bSqrt', mpiAll_opt=.false.)
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

    call msg_memUsage('bdiff_bSqrtAd', mpiAll_opt=.false.)
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
              field_r4(jlon,jlat,jlev2) = real(gd(jlon,jlat,jlev),4)
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
    integer :: jn, jlat, jlon, jlev
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
    call mmpi_bcast(gd_mpiGlobal)

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

  !--------------------------------------------------------------------------
  ! bdiff_getBGstdFromFourSeasons
  !--------------------------------------------------------------------------
  subroutine bdiff_getSSTBGstdFromFourSeasons(hco, vco, stateVector_opt)
    !
    !:Purpose: to get SST 2D background error standard deviation field
    !          stored on Z, U or G grid and interpolate it to the analysis grid.
    !          The operationally used BG std field is updated daily using
    !          four seasonal fields. The estimation for current day is computed
    !          as a weighted average between two fields: from the left and from the right
    !          of the current date.
    !
    implicit none

    ! Arguments:
    type(struct_hco), pointer , intent(in)    :: hco             ! horizontal coordinates
    type(struct_vco), pointer , intent(in)    :: vco             ! vertical coordinates
    type(struct_gsv), optional, intent(inout) :: stateVector_opt ! state vector with resulting estimation
                                                                 ! of SST BG std

    ! Locals:
    integer          :: indexLeft, indexRight
    type(struct_gsv) :: stateVectorLeft, stateVectorRight
    real(8), pointer :: field3DLeft_r8_ptr(:,:,:), field3DRight_r8_ptr(:,:,:)
    real(8), pointer :: ptr_r8(:,:,:)
    real(8)          :: minStddev, maxStddev
    integer          :: hour, day, month, yyyy, ndays, indexMonth
    integer          :: dateStamp ! date stamp for the current day
    real(8)          :: weight
    logical          :: updatedDatestampNov
    real(8)          :: deltaDatestampLeft, deltaRightLeft, delta

    type :: bgops      ! operationally used estimations of background error standard deviation
                       ! are valid on the 15th of Feb, 15th May, 15th Aug and 15th Nov.
      character(len=3) :: validMonthName(4)   = (/    'Feb',    'May',    'Aug',    'Nov' /)
      character(len=6) :: etiket(4)           = (/ 'STDFEB', 'STDMAY', 'STDAUG', 'STDNOV' /)
      integer          :: validMonthNumber(4) = (/        2,        5,        8,       11 /)
      integer          :: validDay  = 15
      integer          :: validHour = 0
      integer          :: dataStamp(4)
    end type bgops
    type(bgops) :: bgstd

    call msg('bdiff_getSSTBGstdFromFourSeasons', 'Reading SST BG std fields from ./bgstddev...')

    dateStamp = tim_getDateStamp()
    call tim_dateStampToYYYYMMDDHH(dateStamp, hour, day, month, ndays, yyyy)
    call msg('bdiff_getSSTBGstdFromFourSeasons', &
             'Interpolating SST BG std field for day/month/year (datestamp): '//&
                                                 str(day) // ',' // str(month) //','//&
                                                 str(yyyy) //', ('//str(dateStamp)//')')

    ! compute datestamps for 4 seasonal fields in bgstddev file
    do indexMonth = 1, 4
      bgstd%dataStamp(indexMonth) = tim_yyyymmddhhToDatestamp(yyyy, bgstd%validMonthNumber(indexMonth), &
                                                              bgstd%validDay, bgstd%validHour)
    end do

    ! compute lower and upper position indexes of the seasonal BG std fields
    updatedDatestampNov = .false.
    do indexMonth = 1, 4
      call difdatr(dateStamp, bgstd%dataStamp(indexMonth), delta)
      if (delta <= 0.d0) then
        indexRight = indexMonth
        if (indexRight == 1) then
          call msg('bdiff_getSSTBGstdFromFourSeasons', 'Seasonal field from the right is in Feb.')
          indexLeft = 4
          call msg('bdiff_getSSTBGstdFromFourSeasons', &
                   'It requires an update of the corresponding bgstd%datestamp(4): '//&
                   str(bgstd%dataStamp(indexLeft))//&
                   'to Nov of the previous year: '//str(yyyy - 1))
          bgstd%dataStamp(indexLeft) = tim_yyyymmddhhToDatestamp(yyyy - 1, bgstd%validMonthNumber(indexLeft), &
                                                                 bgstd%validDay, bgstd%validHour)
          call msg('bdiff_getSSTBGstdFromFourSeasons', 'Updated bgstd%datestamp(4): '//&
                                                       str(bgstd%dataStamp(indexLeft)))
          updatedDatestampNov = .true.
        else
          indexLeft = indexRight - 1
        end if
        exit
      end if
    end do

    ! for dates in range 16 Nov - 31 Dec:
    call difdatr(dateStamp, bgstd%dataStamp(4), delta)
    if (delta > 0.d0 .and. .not. updatedDatestampNov) then
      indexLeft  = 4
      indexRight = 1
      call msg('bdiff_getSSTBGstdFromFourSeasons', &
               'Seasonal field from the right is in Feb of the next year.'//&
               'Updating bgstd%datestamp(1): '//str(bgstd%dataStamp(indexRight)))
      bgstd%dataStamp(indexRight) = tim_yyyymmddhhToDatestamp(yyyy + 1, bgstd%validMonthNumber(indexRight), &
                                                              bgstd%validDay, bgstd%validHour)
      call msg('bdiff_getSSTBGstdFromFourSeasons', &
               'Updated bgstd%datestamp(1): '//str(bgstd%dataStamp(indexRight)))
    end if

    call difdatr(dateStamp, bgstd%dataStamp(indexLeft), deltaDatestampLeft)
    call difdatr(bgstd%dataStamp(indexRight), bgstd%dataStamp(indexLeft), deltaRightLeft)

    if (deltaDatestampLeft < 0. .or. deltaRightLeft < 0.) then
      call rti_abort('bdiff_getSSTBGstdFromFourSeasons: Confusion! '//&
                     'Both distances between two dates must be positive! '//&
                     str(deltaDatestampLeft)//', '//str(deltaRightLeft))
    end if

    weight  = deltaDatestampLeft / deltaRightLeft

    call msg('bdiff_getSSTBGstdFromFourSeasons', 'BG std will be computed between: '//&
                                                 bgstd%validMonthName(indexLeft)//' and '//&
                                                 bgstd%validMonthName(indexRight))
    call msg('bdiff_getSSTBGstdFromFourSeasons', &
             'Weight for: '//bgstd%validMonthName(indexLeft) //': '//str(1. - weight))
    call msg('bdiff_getSSTBGstdFromFourSeasons', &
             'Weight for: '//bgstd%validMonthName(indexRight)//': '//str(     weight))

    ! read BG std from the left
    call msg('bdiff_getSSTBGstdFromFourSeasons', 'Reading '//bgstd%validMonthName(indexLeft)//&
                                                 ' BG std field with etiket: '//bgstd%etiket(indexLeft))
    call gsv_allocate(stateVectorLeft, 1, hco, vco, dataKind_opt = 8, &
                      datestamp_opt = -1, mpi_local_opt = .true., &
                      hInterpolateDegree_opt = 'LINEAR', varNames_opt = (/'TM'/))
    call gio_readFromFile(stateVectorLeft, './bgstddev', bgstd%etiket(indexLeft), &
                          ' ', unitConversion_opt=.false., containsFullField_opt=.true.)
    call gsv_getField(stateVectorLeft, field3DLeft_r8_ptr)

    ! read BG std from the right
    call msg('bdiff_getSSTBGstdFromFourSeasons', 'Reading '//bgstd%validMonthName(indexRight)//&
                                                 ' BG std field with etiket: '//bgstd%etiket(indexRight))
    call gsv_allocate(stateVectorRight, 1, hco, vco, dataKind_opt = 8, &
                      datestamp_opt = -1, mpi_local_opt = .true., &
                      hInterpolateDegree_opt = 'LINEAR', varNames_opt = (/'TM'/))
    call gio_readFromFile(stateVectorRight, './bgstddev', bgstd%etiket(indexRight), &
                          ' ', unitConversion_opt=.false., containsFullField_opt=.true.)
    call gsv_getField(stateVectorRight, field3DRight_r8_ptr)

    if (present(stateVector_opt)) then
      call gsv_getField(stateVector_opt, ptr_r8)
      ptr_r8(:,:,1) = (1. - weight) * field3DLeft_r8_ptr(:,:,1) + &
                            weight  * field3DRight_r8_ptr(:,:,1)
    else
      stddev(:,:,1) = (1. - weight) * field3DLeft_r8_ptr(:,:,1) + &
                            weight  * field3DRight_r8_ptr(:,:,1)
      if (mmpi_nprocs > 1) then
        call mmpi_allreduce(minval(stddev(:,:,1)), minStddev, mmpi_min)
        call mmpi_allreduce(maxval(stddev(:,:,1)), maxStddev, mmpi_max)
      else
        minStddev = minval(stddev(:,:,1))
        maxStddev = maxval(stddev(:,:,1))
      end if
      call msg('bdiff_getSSTBGstdFromFourSeasons', &
               'SST BG std min/max: '//str(minStddev)//', '//str(maxStddev))
    end if

    call gsv_deallocate(stateVectorLeft)
    call gsv_deallocate(stateVectorRight)

    call msg('bdiff_getSSTBGstdFromFourSeasons', 'Completed.')

  end subroutine bdiff_getSSTBGstdFromFourSeasons

end module bMatrixDiff_mod
