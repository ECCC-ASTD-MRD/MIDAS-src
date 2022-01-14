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
!CANADA, H9P 1J3; or send e-mail to ec.service.rpn.ec@canada.ca
!-------------------------------------- LICENCE END --------------------------------------

module diffusion_mod
  ! MODULE diffusion_mod (prefix='diff' category='3. High-level transformations')
  !
  ! :Purpose: Diffusion operator and storage of related diffusion configuration used to model
  !           background-error horizontal correlations. Both explicit and implicit
  !           formulations are included, with implicit being preferred when using large
  !           correlation length scales relative to the grid spacing. For implicit the
  !           MPI topology must have npex=1, i.e. 1xNPEYxNUMTHREADS. A 2D MPI topology
  !           can be used for the explicit formulation.
  !
  ! :Reference: Weaver, A. T., and P. Courtier, 2001: Correlation modelling on
  !             the sphere using a generalized diffusion equation.
  !             Q. J. R. Meteorol. Soc., 127, 1815-1846.
  !
  ! :Basic equations: Lcorr^2 = 2*k*dt*numt   (1)
  !                   stab    = k*dt/dx^2     (2)
  !
  use mpi_mod
  use mpivar_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use oceanMask_mod
  use earthConstants_mod
  use randomNumber_mod
  use utilities_mod
  use gridStateVector_mod

  implicit none
  save
  private

  ! Public subroutines and functions
  public :: diff_setup, diff_finalize, diff_Csqrt, diff_Csqrtadj

  type struct_diff
    ! number of grid points in x and y directions (includes perimeter of land points)
    integer :: ni, nj
    integer :: myLonBeg, myLonEnd, myLatBeg, myLatEnd
    integer :: lonPerPE, lonPerPEmax, latPerPE, latPerPEmax
    integer :: myLonBeg_transpose, myLonEnd_transpose
    integer :: lonPerPE_transpose, lonPerPEmax_transpose
    integer :: numt
    real(8) :: dt
    real(8), allocatable :: cosyhalf(:), cosyinv(:), cosyinvsq(:)
    real(8), allocatable :: mhalfx(:,:), mhalfy(:,:)
    real(8), allocatable :: khalfx(:,:), khalfy(:,:)
    real(8) :: dlon, dlat        ! grid spacing in radians
    ! These are used in subroutines diff_Csqrt and diff_Csqrtadj
    real(8), allocatable :: Winv(:,:), Wsqrt(:,:), Winvsqrt(:,:)
    real(8), allocatable :: Lambda(:,:)
    real(8), allocatable :: diff1x_ap(:,:),diff1x_bp_inv(:,:)
    real(8), allocatable :: diff1x_c(:,:)
    real(8), allocatable :: diff1y_ap(:,:),diff1y_bp_inv(:,:)
    real(8), allocatable :: diff1y_c(:,:)
    logical :: limplicit
  end type struct_diff

  integer, parameter :: nMaxDiff = 10
  integer            :: nDiffAlreadyAllocated = 0
  type(struct_diff)  :: diff(nMaxDiff)

!*************************************************************************

contains

  integer function diff_setup (variableIndex, bdiff_varNameList, hco, vco, corr_len, stab, numberSamples, limplicit)

    implicit none

    ! Arguments:
    integer,                   intent(in)    :: variableIndex        ! Variable index in bdiff_varNameList(:)
    character(len=4),          intent(in)    :: bdiff_varNameList(:) ! list of 2D analysis variables  
    type(struct_hco), pointer, intent(inout) :: hco                  ! Horizontal grid structure
    type(struct_vco), pointer, intent(in)    :: vco                  ! Vertical grid structure
    real,                      intent(in)    :: corr_len             ! Horizontal correlation length scale (km);
                                                                     ! if it is equal to -1, a 2D field of it is read from a file
    real,                      intent(in)    :: stab                 ! Stability criteria (definitely < 0.5)
    integer,                   intent(in)    :: numberSamples        ! Number of samples to estimate normalization factors by randomization.
    logical,                   intent(in)    :: limplicit            ! Indicate to use the implicit formulation

    ! Locals:    
    real(8), allocatable :: latr(:) ! latitudes on the analysis rotated grid, in radians
    real(4), allocatable :: buf2d(:,:)
    integer :: lonIndex, latIndex, sampleIndex, timeStep
    real(8) :: mindxy, maxL, currentMin, currentLatSpacing, currentLonSpacing
    real(8) :: a,b
    real(8), allocatable :: Lcorr(:,:)
    real(8), allocatable :: kappa(:,:)
    real(8), allocatable :: W(:,:)
    real(8), allocatable :: m(:,:)
    real(8), allocatable :: xin(:,:), xin_transpose(:,:)
    real(8), allocatable :: lambdaLocal(:,:) ! auxiliary variable to to MPI_ALLREDUCE of diff % Lambda

    ! diff_norm_fact is the name of the RPN format file for the normalization factors.
    character(len=*), parameter :: diff_norm_fact = './diffusmod.std'

    ! Variables and functions required to write to RPN Standard files.

    integer  :: nmax

    integer  :: std_unit, ierr, ikey, npak, datyp
    integer  :: nii, njj, nkk
    integer  :: dateo
    integer  :: deet, npas
    integer  :: ip1, ip2, ip3
    integer  :: ig1, ig2, ig3, ig4
    character(len=1) :: grtyp
    character(len=2) :: typvar
    character(len=12) :: etiket
    logical  :: rewrit, file_exist
    real     :: dumwrk(1)

    integer  :: fnom,fstouv,fstecr,fstlir,fstfrm,fclos
    external :: fnom,fstouv,fstecr,fstlir,fstfrm,fclos

    integer  :: diffID

    integer :: ni, nj
    integer :: nsize
    integer :: myLonBeg, myLonEnd, myLatBeg, myLatEnd
    integer :: seed                                   ! Seed for random number generator

    character(len=*), parameter :: myName = 'diff_setup'
    character(len=*), parameter :: correlationLengthFileName = './bgstddev'
    type(struct_gsv)            :: statevector
    type(struct_ocm)            :: oceanMask
    real(4),pointer             :: field3D_r4_ptr(:,:,:)

    if ( nDiffAlreadyAllocated == nMaxDiff ) then
      write(*,*) myName//': the maximum number of diffusion operators have already been allocated! ',nMaxDiff
      call utl_abort( myName )
    end if

    nDiffAlreadyAllocated = nDiffAlreadyAllocated + 1
    diffID = nDiffAlreadyAllocated

    ni = hco % ni
    nj = hco % nj
    diff(diffID) % dlon = hco % dlon
    diff(diffID) % dlat = hco % dlat
    diff(diffID) % limplicit = limplicit

    ! domain partionning
    diff(diffID) % ni = ni
    diff(diffID) % nj = nj
    call mpivar_setup_latbands(diff(diffID) % nj, diff(diffID) % latPerPE,  &
         diff(diffID) % latPerPEmax, diff(diffID) % myLatBeg, diff(diffID) % myLatEnd)
    call mpivar_setup_lonbands(diff(diffID) % ni, diff(diffID) % lonPerPE,  &
         diff(diffID) % lonPerPEmax, diff(diffID) % myLonBeg, diff(diffID) % myLonEnd)

    ! also, determine lonIndex begin and end for when array is transposed (for implicit diffusion only)
    call mpivar_setup_latbands(diff(diffID) % ni, diff(diffID) % lonPerPE_transpose,  &
         diff(diffID) % lonPerPEmax_transpose, diff(diffID) % myLonBeg_transpose, diff(diffID) % myLonEnd_transpose)
	 
    myLonBeg = diff(diffID) % myLonBeg
    myLonEnd = diff(diffID) % myLonEnd
    myLatBeg = diff(diffID) % myLatBeg
    myLatEnd = diff(diffID) % myLatEnd
    
    write(*,*) myName//' ***** Starting using the following parameters: *****'
    write(*,*) myName//' Variable : ', bdiff_varNameList( variableIndex ) 
    write(*,*) myName//' Horizontal correlation length scale (km): ', corr_len 
    write(*,*) myName//' Stability criteria: ', stab
    write(*,*) myName//' Indicate implicit diffusion operator (.true.) or explicit version (.false.).: ', limplicit
    write(*,*) myName//' ni/nj: ', ni, nj   
    write(*,*) myName//' ### MPI domain partitionning ###:'
    write(*,*) myName//' [ myLonBeg, myLonEnd ]: [ ', myLonBeg, ' ', myLonEnd, ' ]'  
    write(*,*) myName//' [ myLatBeg, myLatEnd ]: [ ', myLatBeg, ' ', myLatEnd, ' ]'

    ! For implicit diffusion we only allow decomposition by latitude bands
    if ( limplicit .and.  mpi_npex > 1 ) then
      call utl_abort( myName//' Error: for implicit diffusion NPEX must be 1 (i.e. 1xNPEYxNUMTHREADS)' )
    end if
    
    allocate( diff(diffID) % cosyhalf ( nj         ), diff(diffID) % cosyinv ( nj )       , diff(diffID) % cosyinvsq ( nj )    )
    allocate( diff(diffID) % Winv     ( ni    , nj ), diff(diffID) % Wsqrt   ( ni, nj )   , diff(diffID) % Winvsqrt  ( ni, nj ))
    allocate( diff(diffID) % khalfx   ( ni - 1, nj ), diff(diffID) % khalfy  ( ni, nj -1 ) )
    allocate( diff(diffID) % mhalfx   ( ni - 1, nj ), diff(diffID) % mhalfy  ( ni, nj -1 ) )
    allocate( diff(diffID) % Lambda   ( ni    , nj ) )
    allocate( lambdaLocal               ( ni    , nj ) )

    allocate( latr( nj )      )
    allocate( Lcorr( ni, nj ) )
    allocate( kappa( ni, nj ) )
    allocate(     W( ni, nj ) )
    allocate(     m( ni, nj ) )
    allocate(   xin( myLonBeg : myLonEnd, myLatBeg : myLatEnd ) )

    allocate( diff(diffID) % diff1x_ap     ( ni, nj ) )
    allocate( diff(diffID) % diff1x_bp_inv ( ni, nj ) )
    allocate( diff(diffID) % diff1x_c      ( ni, nj ) )
    allocate( diff(diffID) % diff1y_ap     ( nj, ni ) )
    allocate( diff(diffID) % diff1y_bp_inv ( nj, ni ) )
    allocate( diff(diffID) % diff1y_c      ( nj, ni ) )

    latr(:) = hco % lat(:)

    diff(diffID) % cosyinv(:)   = 1.0d0 / cos( latr(:) )
    diff(diffID) % cosyinvsq(:) = diff(diffID) % cosyinv(:) * diff(diffID) % cosyinv(:)

    ! cosinus of latitudes on staggered grid
    diff(diffID) % cosyhalf(:) = cos( latr(:) + 0.5d0 * diff(diffID) % dlat )

    ! Get mask from analysisgrid file
    call ocm_readMaskFromFile(oceanMask, hco, vco, './analysisgrid')

    ! land mask (1=water, 0=land)
    do latIndex = 1, nj
      do lonIndex = 1, ni
          
        if ( oceanMask%mask ( lonIndex, latIndex, 1 ) ) then
          m ( lonIndex, latIndex ) = 1.0d0
        else
          m ( lonIndex, latIndex ) = 0.0d0
        end if
       
      end do
    end do
    call ocm_deallocate(oceanMask)
   
    m (  :, 1  ) = 0.0d0
    m (  1, :  ) = 0.0d0
    m (  :, nj ) = 0.0d0
    m ( ni, :  ) = 0.0d0

    ! define mask on staggered grids
    do latIndex = 1, nj
      do lonIndex = 1, ni - 1
        if ( sum ( m ( lonIndex : lonIndex + 1, latIndex ) ) < 2.0d0 ) then
          diff (diffID) % mhalfx( lonIndex, latIndex ) = 0.0d0
        else
          diff (diffID) % mhalfx( lonIndex, latIndex ) = 1.0d0
        end if
      end do
    end do

    do latIndex = 1, nj - 1
      do lonIndex = 1, ni
        if ( sum( m ( lonIndex, latIndex : latIndex + 1 ) ) < 2.0d0 ) then
          diff (diffID) % mhalfy ( lonIndex, latIndex ) = 0.0d0
        else
          diff (diffID) % mhalfy ( lonIndex, latIndex ) = 1.0d0
        end if
      end do
    end do

    ! minimum grid spacing over the domain
    mindxy = 100000.0d0
    do latIndex = 1, nj - 1
      do lonIndex = 1, ni - 1

        if ( ( diff (diffID) % mhalfy ( lonIndex, latIndex ) == 1.0d0 ) .and. ( diff (diffID) % mhalfx ( lonIndex, latIndex ) == 1.0d0 ) ) then

          currentLonSpacing = cos( latr( latIndex ) ) * diff(diffID) % dlon
          currentLatSpacing =                           diff(diffID) % dlat
          currentMin = min ( currentLatSpacing, currentLonSpacing )  

          if ( currentMin < mindxy ) then
            mindxy = currentMin
          end if

        end if

      end do
    end do

    mindxy = min( mindxy, diff(diffID) % dlat )
    write(*,*) myName//': Minimim grid spacing: mindxy = ', mindxy

    if ( corr_len == -1 ) then

      write(*,*) myName//': Correlation length scale 2D field will be read from the file: ', correlationLengthFileName
      call gsv_allocate( statevector, 1, hco, vco, dateStamp_opt=-1, dataKind_opt=4, &
                         hInterpolateDegree_opt='LINEAR', varNames_opt=bdiff_varNameList, &
                         mpi_local_opt=.false. )
      call gsv_zero( statevector )
      call gsv_readFromFile(statevector, correlationLengthFileName, 'CORRLEN', ' ', unitConversion_opt = .false. )

      call gsv_getField( statevector, field3D_r4_ptr, bdiff_varNameList( variableIndex ) )
      Lcorr(:,:) = dble( field3D_r4_ptr( :, :, 1 ) )       
      write(*,*) myName//': correlation length scale 2D field for variable ', bdiff_varNameList( variableIndex ),' min/max: ', &
                 minval( Lcorr(:,:) ), maxval( Lcorr(:,:) )
      call gsv_deallocate( statevector )

    else
      
      Lcorr(:,:) = corr_len

   end if

    Lcorr(:,:) =  Lcorr(:,:) / ( ec_rayt / 1000.0 )     ! lengthscale in radians
    maxL = maxval( Lcorr( 2 : ni - 1, 2 : nj - 1 ) ) ! maximum lengthscale over domain

    ! set main parameters for diffusion operator
    kappa(:,:) = Lcorr(:,:)**2                                              ! arbitrarily set k to L^2 (in radians)
    diff(diffID) % dt = stab * ( mindxy**2 ) / ( maxL**2 )                ! determine dt from stability criteria (2)
    ! diff(diffID)%numt = 1.0d0/(2.0d0*diff(diffID)%dt)                     ! determine number of timesteps from (1)
    diff(diffID) % numt = ceiling(1.0d0/(4.0d0*diff(diffID)%dt))*2        ! make sure it is an even integer
    diff(diffID) % dt = 1.0d0 / ( 2.0d0 * dble( diff(diffID) % numt ) ) ! recompute dt

    ! interpolate diffusion coefficient onto 2 staggered lat-lon grids
    do latIndex = 1, nj
      do lonIndex = 1, ni - 1
        diff(diffID) % khalfx( lonIndex, latIndex ) = ( kappa( lonIndex, latIndex ) + kappa( lonIndex + 1, latIndex ) ) / 2.0d0
      end do
    end do
    do latIndex = 1, nj - 1
      do lonIndex = 1, ni
        diff(diffID) % khalfy( lonIndex, latIndex ) = ( kappa( lonIndex, latIndex ) + kappa( lonIndex, latIndex + 1 ) ) / 2.0d0
      end do
    end do

    ! print this stuff in listing file for user information:
    if ( .not. limplicit ) then
      write(*,*)
      write(*,*) myName//': Number of timesteps = ', diff(diffID) % numt 
      write(*,*) myName//': Stability           = ', maxval( kappa ) * diff(diffID) % dt / ( mindxy**2 )
      write(*,*)
    end if

    ! this is the matrix necessary for defining the inner product: for lat-lon grid, only cos(y)
    ! Actually, only the diagonal of the matrix is stored in the array W, since the matrix is diagonal.
    W(1,:) = cos( latr(:) )
    do lonIndex = 2, ni
      W( lonIndex, : ) = W( 1, : )
    end do
    diff(diffID) % Winv(:,:)     = 1.0d0 / W(:,:)
    diff(diffID) % Wsqrt(:,:)    = sqrt( W(:,:) )
    diff(diffID) % Winvsqrt(:,:) = 1.0d0 / diff(diffID) % Wsqrt(:,:)

    ! specify number of timesteps and timestep length for implicit 1D diffusion
    if ( limplicit ) then
      diff(diffID) % numt = 5
      diff(diffID) % dt   = 1.0d0 / ( 2.0d0 * dble( 2 * diff(diffID) % numt ) - 3.0d0 )
    end if

    ! compute the LU decomposition for the implicit 1D diffusion
    diff(diffID) % diff1x_ap(:,:) = 0.0d0
    diff(diffID) % diff1x_bp_inv(:,:) = 0.0d0      
    diff(diffID) % diff1x_c(:,:) = 0.0d0
    diff(diffID) % diff1y_ap(:,:) = 0.0d0
    diff(diffID) % diff1y_bp_inv(:,:) = 0.0d0      
    diff(diffID) % diff1y_c(:,:) = 0.0d0

    !$OMP PARALLEL DO PRIVATE(latIndex,lonIndex,a,b)
    do latIndex = 2, nj - 1
      lonIndex = 2
      diff(diffID) % diff1x_bp_inv( lonIndex, latIndex ) = 1.0d0 / ( 1.0d0 + &
          diff(diffID) % dt * diff(diffID) % cosyinvsq( latIndex ) * ( diff(diffID) % mhalfx( lonIndex, latIndex ) * diff(diffID) % khalfx( lonIndex, latIndex ) + &
          diff(diffID) % mhalfx( lonIndex - 1, latIndex ) * diff(diffID) % khalfx( lonIndex - 1, latIndex ) ) / ( diff(diffID) % dlon * diff(diffID) % dlon ))
      do lonIndex = 3, ni - 1
        ! elements of the tri-diagonal coefficient matrix
        a = - diff(diffID) % dt * diff(diffID) % cosyinvsq( latIndex ) * diff(diffID) % mhalfx( lonIndex - 1, latIndex ) * diff(diffID) % khalfx( lonIndex - 1, latIndex ) / &
            ( diff(diffID) % dlon *diff(diffID) % dlon )
        b = 1 + diff(diffID) % dt * diff(diffID) % cosyinvsq( latIndex ) * ( diff(diffID) % mhalfx( lonIndex, latIndex ) * diff(diffID) % khalfx( lonIndex, latIndex ) + &
            diff(diffID) % mhalfx( lonIndex - 1, latIndex ) * diff(diffID) % khalfx( lonIndex - 1, latIndex ) ) / ( diff(diffID) % dlon * diff(diffID) % dlon )
        diff(diffID) % diff1x_c( lonIndex, latIndex ) = - diff(diffID) % dt * diff(diffID) % cosyinvsq( latIndex ) * diff(diffID) % mhalfx( lonIndex, latIndex ) * &
            diff(diffID) % khalfx( lonIndex, latIndex ) / ( diff(diffID) % dlon * diff(diffID) % dlon )
        diff(diffID) % diff1x_ap( lonIndex, latIndex ) = a * diff(diffID) % diff1x_bp_inv( lonIndex - 1, latIndex )
        diff(diffID) % diff1x_bp_inv( lonIndex, latIndex ) = 1.0d0/(b-a*diff(diffID)%diff1x_c(lonIndex-1,latIndex) * &
            diff(diffID) % diff1x_bp_inv( lonIndex - 1, latIndex ))
      end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(latIndex,lonIndex,a,b)
    do lonIndex = 2, ni - 1
      latIndex = 2
      diff(diffID) % diff1y_bp_inv( latIndex, lonIndex ) = 1.0d0 / ( 1.0d0 + &
           diff(diffID) % dt * diff(diffID) % cosyinvsq( latIndex ) * ( diff(diffID) % cosyhalf( latIndex ) * diff(diffID) % mhalfy( lonIndex, latIndex ) * &
           diff(diffID) % khalfy( lonIndex, latIndex ) + &
           diff(diffID) % cosyhalf( latIndex - 1 ) * diff(diffID) % mhalfy( lonIndex, latIndex - 1 ) * diff(diffID) % khalfy( lonIndex, latIndex - 1 ) ) / &
           ( diff(diffID) % dlat * diff(diffID) % dlat ) )
      do latIndex = 3, nj-1
        ! elements of the tri-diagonal coefficient matrix
        a = - diff(diffID)%dt*diff(diffID)%cosyinv(latIndex)*diff(diffID)%cosyhalf(latIndex-1)* &
            diff(diffID)%mhalfy(lonIndex,latIndex-1)*diff(diffID)%khalfy(lonIndex,latIndex-1)/(diff(diffID)%dlat*diff(diffID)%dlat)
        b = 1 + diff(diffID)%dt*diff(diffID)%cosyinv(latIndex)*(diff(diffID)%cosyhalf(latIndex)*diff(diffID)%mhalfy(lonIndex,latIndex)* &
            diff(diffID)%khalfy(lonIndex,latIndex) + &
            diff(diffID) % cosyhalf( latIndex - 1 ) * diff(diffID) % mhalfy( lonIndex, latIndex - 1 ) * diff(diffID) % khalfy( lonIndex, latIndex - 1 ) ) &
            / ( diff(diffID) % dlat * diff(diffID) % dlat )
        diff(diffID) % diff1y_c( latIndex, lonIndex ) = - diff(diffID) % dt * diff(diffID) % cosyinv( latIndex ) * diff(diffID) % cosyhalf( latIndex ) * &
            diff(diffID) % mhalfy( lonIndex, latIndex ) * diff(diffID) % khalfy( lonIndex, latIndex ) / ( diff(diffID) % dlat * diff(diffID) % dlat )
        diff(diffID)%diff1y_ap(latIndex,lonIndex) = a*diff(diffID) % diff1y_bp_inv( latIndex - 1, lonIndex )
        diff(diffID) % diff1y_bp_inv( latIndex, lonIndex ) = 1.0d0 / ( b - a * diff(diffID) % diff1y_c( latIndex - 1, lonIndex ) &
            *diff(diffID) % diff1y_bp_inv( latIndex - 1, lonIndex ))
      end do
    end do
    !$OMP END PARALLEL DO

    std_unit = 0
    inquire (file=diff_norm_fact, exist=file_exist)
        
    if ( file_exist ) then

      write(*,*) myName//': file containing normalization factors exists: ', trim(diff_norm_fact)
       
      ierr = fnom(std_unit, diff_norm_fact, 'RND+R/O', 0)
      nmax = fstouv(std_unit, 'RND')

      allocate(buf2d( ni, nj ))
      ! Looking for FST record parameters..
      dateo = -1
      etiket = ''
      ip1 = -1
      ip2 = -1
      ip3 = -1
      grtyp = ' '
      ikey = utl_fstlir_r4( buf2d, std_unit, nii, njj, nkk, dateo, etiket, ip1, ip2, ip3, grtyp, 'LAMB')
      write(*,*) myName//': nii / njj: ', nii, njj
      write(*,*) myName//': ni  / nj : ', ni, nj
      if ( ikey <= 0 ) write(*,*) myName//': Attention! ikey = ', ikey

      ierr = fstfrm(std_unit)
      ierr = fclos(std_unit)

      diff(diffID) % Lambda = dble( buf2d )

      deallocate( buf2d )

    end if

    if ( nii /= ni .or. njj /= nj .or. ikey <= 0 .or. ( .not. file_exist) ) then

      if ( .not. file_exist) write(*,*) myName//': file containing normalization factors does not exist!!! ',  trim(diff_norm_fact)
      
      seed = 1 
      call rng_setup( abs( seed + mpi_myid ))

      write(*,*) myName//': Number of samples, ni * nj: ', numberSamples, ni * nj


      ! compute normalization:  Lambda = inverse stddev of (Diffuse * W^-1/2)
      write(*,*)  myName//': Estimate normalization factors for diffusion using randomization method...'
      write(*,*)  myName//': will use ', numberSamples,' samples.',' ni and nj: ', ni, nj
      call flush(6)
      diff(diffID) % Lambda = 0.0d0
      lambdaLocal             = 0.0d0

      SAMPLE: do sampleIndex = 1, numberSamples

        if ( modulo( sampleIndex, 100 ) == 0 ) write(*,*) myName//': Computing sample: ', sampleIndex

        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            xin( lonIndex, latIndex ) = diff(diffID) % Winvsqrt( lonIndex, latIndex ) * rng_gaussian()
          end do
        end do

        if ( limplicit ) then

          allocate(xin_transpose(diff(diffID)%myLonBeg_transpose:diff(diffID)%myLonEnd_transpose,diff(diffID)%nj))
          do timeStep = 1, diff(diffID) % numt
	    
            call diffusion1x_implicit( diffID, xin, xin )

            call transposeLatToLonBands( diffID, xin, xin_transpose )
            call diffusion1y_implicit( diffID, xin_transpose, xin_transpose )
            call transposeLonToLatBands( diffID, xin_transpose, xin )

          end do
          deallocate(xin_transpose)

        else

          call diffusion_explicit( diffID, xin, xin )

        end if

        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd

            diff(diffID) % Lambda( lonIndex, latIndex ) = diff(diffID) % Lambda( lonIndex, latIndex ) + xin( lonIndex, latIndex ) * xin( lonIndex, latIndex )

	  end do
	end do

      end do SAMPLE

      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd

          diff(diffID) % Lambda( lonIndex, latIndex ) = sqrt( diff(diffID) % Lambda( lonIndex, latIndex ) / dble( numberSamples - 1 ) ) ! normalization: inverse of rms of ens

          if ( diff(diffID) % Lambda( lonIndex, latIndex ) > 0.0d0 ) then
            diff(diffID) % Lambda( lonIndex, latIndex ) = 1.0d0 / diff(diffID) % Lambda( lonIndex, latIndex )
          end if

        end do
      end do

      lambdaLocal( myLonBeg : myLonEnd, myLatBeg : myLatEnd ) = diff(diffID) % Lambda( myLonBeg : myLonEnd, myLatBeg : myLatEnd )    
      nsize = ni * nj
      call rpn_comm_allreduce( lambdaLocal, diff(diffID) % Lambda, nsize, "mpi_double_precision", "mpi_sum", "GRID", ierr )

      if ( mpi_myid == 0 ) then

	write(*,*)  myName//': Save normalization coefficient on proc ', mpi_myid, ' into the file: ', trim(diff_norm_fact)
        npak = 0
        dateo = 0
        deet = 0
        npas = 0
        ip1 = 0
        ip2 = 0
        ip3 = 0
        typvar = 'X'
        grtyp = 'X'
        ig1 = 0
        ig2 = 0
        ig3 = 0
        ig4 = 0
        datyp = 1
        rewrit = .FALSE.

        if ( limplicit ) then
          write (etiket, FMT='(''KM'',i3.3,''IMPLICI'')') int(corr_len)
        else
          write (etiket, FMT='(''KM'',i3.3,''STAB'',f3.1)') int(corr_len), stab
        end if

        ierr = fnom( std_unit, diff_norm_fact, 'RND', 0 )
        nmax = fstouv( std_unit, 'RND')

        ierr = fstecr( real(diff(diffID) % Lambda ), dumwrk, npak, std_unit,          &
                       dateo, deet, npas, NI, NJ, 1, ip1, ip2, ip3,                     &
                       typvar, 'LAMB', etiket, grtyp, ig1, ig2, ig3, ig4, datyp, rewrit )

        ierr = fstfrm( std_unit )
        ierr = fclos( std_unit )

      end if 

    end if

    deallocate( xin )
    deallocate( m )
    deallocate( W )
    deallocate( kappa )
    deallocate( Lcorr )
    deallocate( latr )
    deallocate( lambdaLocal )

    diff_setup = diffID

    write(*,*) myName//' ***** END *****'

  end function diff_setup


  subroutine diff_finalize(diffID)
    !
    !:Purpose: To finalize the diffusion operator module, and to free up memory.
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: diffID

    write(*,*) 'diff_finalize: deallocating arrays fordiffID= ', diffID

    deallocate( diff(diffID) % diff1y_c      )
    deallocate( diff(diffID) % diff1y_bp_inv )
    deallocate( diff(diffID) % diff1y_ap     )
    deallocate( diff(diffID) % diff1x_c      )
    deallocate( diff(diffID) % diff1x_bp_inv )
    deallocate( diff(diffID) % diff1x_ap     )

    deallocate( diff(diffID) % Lambda )
    deallocate( diff(diffID) % mhalfy, diff(diffID) % mhalfx )
    deallocate( diff(diffID) % khalfy, diff(diffID) % khalfx )
    deallocate( diff(diffID) % Winvsqrt, diff(diffID) % Wsqrt, diff(diffID)%Winv )
    deallocate( diff(diffID) % cosyinvsq, diff(diffID) % cosyinv, diff(diffID) % cosyhalf )

  end subroutine diff_finalize


  subroutine diffusion_explicit( diffID, xin, xout )
    !
    !:Purpose: To compute Lsqrt*xin (diffusion over numt/2 timesteps), and to
    !          specify initial conditions
    implicit none

    ! Arguments:
    integer, intent(in)  :: diffID
    real(8), intent(in)  :: xin ( :, : )
    real(8), intent(out) :: xout( :, : )
    
    ! Locals:
    integer :: timeIndex, latIndex, lonIndex
    integer :: myLonBeg, myLonEnd, myLatBeg, myLatEnd
    integer :: myLonBegNoB, myLonEndNoB, myLatBegNoB, myLatEndNoB
    integer :: lonPerPE, latPerPE
    integer :: sendToPE, recvFromPE, sendTag, recvTag, status, ierr
    real(8) :: xhalo( (diff(diffID)%myLonBeg-1):(diff(diffID)%myLonEnd+1), (diff(diffID)%myLatBeg-1):(diff(diffID)%myLatEnd+1) )
    real(8) :: xlast( (diff(diffID)%myLonBeg-1):(diff(diffID)%myLonEnd+1), (diff(diffID)%myLatBeg-1):(diff(diffID)%myLatEnd+1) )
    real(8), allocatable :: sendBufLon(:), recvBufLon(:), sendBufLat(:), recvBufLat(:)

    call tmg_start(185, 'diffusion_explicit' )

    lonPerPE = diff(diffID) % lonPerPE
    latPerPE = diff(diffID) % latPerPE

    myLonBeg = diff(diffID) % myLonBeg
    myLonEnd = diff(diffID) % myLonEnd
    myLatBeg = diff(diffID) % myLatBeg
    myLatEnd = diff(diffID) % myLatEnd

    ! remove global border from range of grid points where output is calculated
    myLonBegNoB = max(myLonBeg, 2)
    myLonEndNoB = min(myLonEnd, diff(diffID)%ni-1)
    myLatBegNoB = max(myLatBeg, 2)
    myLatEndNoB = min(myLatEnd, diff(diffID)%nj-1)

    ! setup for mpi halo exchange
    allocate(sendBufLon(lonPerPE))
    allocate(recvBufLon(lonPerPE))
    allocate(sendBufLat(latPerPE))
    allocate(recvBufLat(latPerPE))

    xhalo(:,:) = 0.0d0
    xlast(:,:) = 0.0d0
    xlast( myLonBeg:myLonEnd, myLatBeg:myLatEnd ) = xin(:,:)
    
    ! iterate difference equations
    TIME: do timeIndex = 1, diff(diffID) % numt / 2

      ! exchange 4 arrays of halo information between mpi tasks

      ! halo array #1: send western edge interior, recv eastern edge halo
      if (mpi_npex > 1) then
        sendToPE   = mpi_myidx - 1
        recvFromPE = mpi_myidx + 1
        sendTag = (mpi_myidx)*10000 + (mpi_myidx - 1)
        recvTag = (mpi_myidx + 1)*10000 + (mpi_myidx)
        if (mpi_myidx == 0) then
          ! western-most tile only recv
          call rpn_comm_recv( recvBufLat, latPerPE, "mpi_real8", recvFromPE, recvTag, "EW", status, ierr )
          xlast(myLonEnd+1,myLatBeg:myLatEnd) = recvBufLat(:)
        else if (mpi_myidx == (mpi_npex-1)) then
          ! eastern-most tile only send
          sendBufLat(:) = xlast(myLonBeg,myLatBeg:myLatEnd)
          call rpn_comm_send( sendBufLat, latPerPE, "mpi_real8", sendToPE, sendTag, "EW", ierr )
        else
          ! interior tiles both send and recv
          sendBufLat(:) = xlast(myLonBeg,myLatBeg:myLatEnd)
          call rpn_comm_sendrecv( sendBufLat, latPerPE, "mpi_real8", sendToPE,   sendTag, &
                                  recvBufLat, latPerPE, "mpi_real8", recvFromPE, recvTag, &
                                  "EW", status, ierr )
          xlast(myLonEnd+1,myLatBeg:myLatEnd) = recvBufLat(:)
        end if
      end if

      ! halo array #2: send eastern edge interior, recv western edge halo
      if (mpi_npex > 1) then
        sendToPE   = mpi_myidx + 1
        recvFromPE = mpi_myidx - 1
        sendTag = (mpi_myidx)*10000 + (mpi_myidx + 1)
        recvTag = (mpi_myidx - 1)*10000 + (mpi_myidx)
        if (mpi_myidx == (mpi_npex-1)) then
          ! eastern-most tile only recv
          call rpn_comm_recv( recvBufLat, latPerPE, "mpi_real8", recvFromPE, recvTag, "EW", status, ierr )
          xlast(myLonBeg-1,myLatBeg:myLatEnd) = recvBufLat(:)
        else if (mpi_myidx == 0) then
          ! western-most tile only send
          sendBufLat(:) = xlast(myLonEnd,myLatBeg:myLatEnd)
          call rpn_comm_send( sendBufLat, latPerPE, "mpi_real8", sendToPE, sendTag, "EW", ierr )
        else
          ! interior tiles both send and recv
          sendBufLat(:) = xlast(myLonEnd,myLatBeg:myLatEnd)
          call rpn_comm_sendrecv( sendBufLat, latPerPE, "mpi_real8", sendToPE,   sendTag, &
                                  recvBufLat, latPerPE, "mpi_real8", recvFromPE, recvTag, &
                                  "EW", status, ierr )
          xlast(myLonBeg-1,myLatBeg:myLatEnd) = recvBufLat(:)
        end if
      end if

      ! halo array #3: send southern edge interior, recv northern edge halo
      if (mpi_npey > 1) then
        sendToPE   = mpi_myidy - 1
        recvFromPE = mpi_myidy + 1
        sendTag = (mpi_myidy)*10000 + (mpi_myidy - 1)
        recvTag = (mpi_myidy + 1)*10000 + (mpi_myidy)
        if (mpi_myidy == 0) then
          ! southern-most tile only recv
          call rpn_comm_recv( recvBufLon, lonPerPE, "mpi_real8", recvFromPE, recvTag, "NS", status, ierr )
          xlast(myLonBeg:myLonEnd,myLatEnd+1) = recvBufLon(:)
        else if (mpi_myidy == (mpi_npey-1)) then
          ! northern-most tile only send
          sendBufLon(:) = xlast(myLonBeg:myLonEnd,myLatBeg)
          call rpn_comm_send( sendBufLon, lonPerPE, "mpi_real8", sendToPE, sendTag, "NS", ierr )
        else
          ! interior tiles both send and recv
          sendBufLon(:) = xlast(myLonBeg:myLonEnd,myLatBeg)
          call rpn_comm_sendrecv( sendBufLon, lonPerPE, "mpi_real8", sendToPE,   sendTag, &
                                  recvBufLon, lonPerPE, "mpi_real8", recvFromPE, recvTag, &
                                  "NS", status, ierr )
          xlast(myLonBeg:myLonEnd,myLatEnd+1) = recvBufLon(:)
        end if
      end if

      ! halo array #4: send northern edge interior, recv southern edge halo
      if (mpi_npey > 1) then
        sendToPE   = mpi_myidy + 1
        recvFromPE = mpi_myidy - 1
        sendTag = (mpi_myidy)*10000 + (mpi_myidy + 1)
        recvTag = (mpi_myidy - 1)*10000 + (mpi_myidy)
        if (mpi_myidy == (mpi_npey-1)) then
          ! northern-most tile only recv
          call rpn_comm_recv( recvBufLon, lonPerPE, "mpi_real8", recvFromPE, recvTag, "NS", status, ierr )
          xlast(myLonBeg:myLonEnd,myLatBeg-1) = recvBufLon(:)
        else if (mpi_myidy == 0) then
          ! southern-most tile only send
          sendBufLon(:) = xlast(myLonBeg:myLonEnd,myLatEnd)
          call rpn_comm_send( sendBufLon, lonPerPE, "mpi_real8", sendToPE, sendTag, "NS", ierr )
        else
          ! interior tiles both send and recv
          sendBufLon(:) = xlast(myLonBeg:myLonEnd,myLatEnd)
          call rpn_comm_sendrecv( sendBufLon, lonPerPE, "mpi_real8", sendToPE,   sendTag, &
                                  recvBufLon, lonPerPE, "mpi_real8", recvFromPE, recvTag, &
                                  "NS", status, ierr )
          xlast(myLonBeg:myLonEnd,myLatBeg-1) = recvBufLon(:)
        end if
      end if

      !$OMP PARALLEL DO PRIVATE(latIndex,lonIndex)
      do latIndex = myLatBegNoB, myLatEndNoB
        do lonIndex = myLonBegNoB, myLonEndNoB
          xhalo(lonIndex,latIndex) = xlast(lonIndex,latIndex) + diff(diffID)%dt *                                                             &
                ( diff(diffID)%cosyinvsq(latIndex) *                                                                                          &
                  ( diff(diffID)%mhalfx(lonIndex,latIndex) * diff(diffID)%khalfx(lonIndex,latIndex) *                                         &
                    ( xlast(lonIndex+1,latIndex) - xlast(lonIndex,latIndex) ) / diff(diffID)%dlon -                                           &
                    diff(diffID)%mhalfx(lonIndex-1,latIndex) * diff(diffID)%khalfx(lonIndex-1,latIndex) *                                     &
                    ( xlast(lonIndex,latIndex) - xlast(lonIndex-1,latIndex) ) / diff(diffID)%dlon                                             &
                  ) / diff(diffID)%dlon +                                                                                                     &
                  diff(diffID)%cosyinv(latIndex) *                                                                                            &
                  ( diff(diffID)%cosyhalf(latIndex) * diff(diffID)%mhalfy(lonIndex,latIndex) * diff(diffID)%khalfy(lonIndex,latIndex) *       &
                    ( xlast(lonIndex,latIndex+1) - xlast(lonIndex,latIndex) ) / diff(diffID)%dlat -                                           &
                    diff(diffID)%cosyhalf(latIndex-1) * diff(diffID)%mhalfy(lonIndex,latIndex-1) * diff(diffID)%khalfy(lonIndex,latIndex-1) * &
                    ( xlast(lonIndex,latIndex) - xlast(lonIndex,latIndex-1) ) / diff(diffID)%dlat                                             &
                  ) / diff(diffID)%dlat                                                                                                       &
                )
        end do
      end do
      !$OMP END PARALLEL DO

      xlast(:,:) = xhalo(:,:)

    end do TIME

    xout(:,:) = xlast( myLonBeg:myLonEnd, myLatBeg:myLatEnd )

    deallocate(sendBufLon)
    deallocate(recvBufLon)
    deallocate(sendBufLat)
    deallocate(recvBufLat)

    call tmg_stop(185)

  end subroutine diffusion_explicit


  subroutine diff_Csqrt( diffID, xin, xout )

    implicit none

    ! Arguments
    integer, intent(in)  :: diffID
    real(8), intent(in)  :: xin ( diff(diffID)%myLonBeg:diff(diffID)%myLonEnd, diff(diffID)%myLatBeg:diff(diffID)%myLatEnd )
    real(8), intent(out) :: xout( diff(diffID)%myLonBeg:diff(diffID)%myLonEnd, diff(diffID)%myLatBeg:diff(diffID)%myLatEnd )

    ! Locals:
    integer :: tIndex
    real(8), allocatable :: xout_transpose(:,:)

    ! compute Csqrt

    ! this is the C^1/2 required for the forward model: Csqrt = Lambda * Diffuse * W^-1/2
    xout(:,:) = xin(:,:) *  &
                diff(diffID) % Winvsqrt(diff(diffID)%myLonBeg:diff(diffID)%myLonEnd, &
                                          diff(diffID)%myLatBeg:diff(diffID)%myLatEnd)
    if ( diff(diffID) % limplicit ) then
      allocate(xout_transpose(diff(diffID)%myLonBeg_transpose:diff(diffID)%myLonEnd_transpose,diff(diffID)%nj))
      do tIndex = 1, diff(diffID) % numt

        call diffusion1x_implicit( diffID, xout, xout )

        call transposeLatToLonBands( diffID, xout, xout_transpose )
        call diffusion1y_implicit( diffID, xout_transpose, xout_transpose )
        call transposeLonToLatBands( diffID, xout_transpose, xout )

      end do
      deallocate(xout_transpose)
    else
      call diffusion_explicit ( diffID, xout, xout )
    end if

    xout(:,:) = xout(:,:) *  &
                diff(diffID) % Lambda(diff(diffID)%myLonBeg:diff(diffID)%myLonEnd, &
                                      diff(diffID)%myLatBeg:diff(diffID)%myLatEnd)

  end subroutine diff_Csqrt


  subroutine diff_Csqrtadj( diffID, xin, xout )

    implicit none

    ! Arguments:
    integer, intent(in)  :: diffID
    real(8), intent(in)  :: xin ( diff(diffID)%myLonBeg:diff(diffID)%myLonEnd, diff(diffID)%myLatBeg:diff(diffID)%myLatEnd )
    real(8), intent(out) :: xout( diff(diffID)%myLonBeg:diff(diffID)%myLonEnd, diff(diffID)%myLatBeg:diff(diffID)%myLatEnd )

    ! Locals:
    integer :: tIndex
    real(8), allocatable :: xout_transpose(:,:)

    ! compute Csqrtadj

    ! this is the (C^1/2)^T required for the adjoint: Csqrt^T = W^1/2 * Diffuse * W^-1 * Lambda
    xout(:,:) = xin (:,:) * diff(diffID) % Lambda(diff(diffID)%myLonBeg:diff(diffID)%myLonEnd, &
                                                    diff(diffID)%myLatBeg:diff(diffID)%myLatEnd) 
    xout(:,:) = xout(:,:) * diff(diffID) %   Winv(diff(diffID)%myLonBeg:diff(diffID)%myLonEnd, &
                                                    diff(diffID)%myLatBeg:diff(diffID)%myLatEnd) 
    if ( diff(diffID) % limplicit ) then
      allocate(xout_transpose(diff(diffID)%myLonBeg_transpose:diff(diffID)%myLonEnd_transpose,diff(diffID)%nj))
      do tIndex = 1, diff(diffID) % numt

        call transposeLatToLonBands( diffID, xout, xout_transpose )
        call diffusion1y_implicit( diffID, xout_transpose, xout_transpose )
        call transposeLonToLatBands( diffID, xout_transpose, xout )

        call diffusion1x_implicit( diffID, xout, xout )

      end do
      deallocate(xout_transpose)
    else
      call diffusion_explicit( diffID, xout, xout )
    end if

    xout(:,:) = xout(:,:) * diff(diffID) % Wsqrt(diff(diffID)%myLonBeg:diff(diffID)%myLonEnd, &
                                                   diff(diffID)%myLatBeg:diff(diffID)%myLatEnd)

  end subroutine diff_Csqrtadj


  subroutine transposeLatToLonBands(diffID, xin, xout)
    implicit none

    ! Arguments:
    integer, intent(in)  :: diffID
    real(8), intent(in)  :: xin ( diff(diffID)%ni, diff(diffID)%myLatBeg:diff(diffID)%myLatEnd )
    real(8), intent(out) :: xout( diff(diffID)%myLonBeg_transpose:diff(diffID)%myLonEnd_transpose, diff(diffID)%nj )

    ! Locals:
    integer :: yourid, ierr, nsize
    integer :: allLonBeg(mpi_nprocs), allLonEnd(mpi_nprocs)
    integer :: allLatBeg(mpi_nprocs), allLatEnd(mpi_nprocs)
    real(8), allocatable :: xsend(:,:,:),xrecv(:,:,:)

    call tmg_start(188,'diff-transposeLatToLon')

    ! Abort if grid cannot be evenly divided over the mpi tasks
    if (diff(diffID)%lonPerPE_transpose /= diff(diffID)%lonPerPEmax_transpose .or. &
        diff(diffID)%latPerPE           /= diff(diffID)%latPerPEmax ) then
      call utl_abort('diffusion_mod-transposeLatToLonBands: grid cannot be evenly distributed over mpi tasks')
    end if

    allocate(xsend(diff(diffID)%lonPerPE_transpose,diff(diffID)%latPerPE, mpi_nprocs))
    allocate(xrecv(diff(diffID)%lonPerPE_transpose,diff(diffID)%latPerPE, mpi_nprocs))

    call rpn_comm_allgather(diff(diffID)%myLonBeg_transpose,1,'mpi_integer',       &
                            allLonBeg                      ,1,'mpi_integer','GRID',ierr)
    call rpn_comm_allgather(diff(diffID)%myLonEnd_transpose,1,'mpi_integer',       &
                            allLonEnd                      ,1,'mpi_integer','GRID',ierr)
    call rpn_comm_allgather(diff(diffID)%myLatBeg,1,'mpi_integer',       &
                            allLatBeg            ,1,'mpi_integer','GRID',ierr)
    call rpn_comm_allgather(diff(diffID)%myLatEnd,1,'mpi_integer',       &
                            allLatEnd            ,1,'mpi_integer','GRID',ierr)

    xout(:,:) = 0.0d0
!    xout(diff(diffID)%myLonBeg_transpose:diff(diffID)%myLonEnd_transpose,diff(diffID)%myLatBeg:diff(diffID)%myLatEnd) = &
!     xin(diff(diffID)%myLonBeg_transpose:diff(diffID)%myLonEnd_transpose,diff(diffID)%myLatBeg:diff(diffID)%myLatEnd)

    !$OMP PARALLEL DO PRIVATE(yourid)
    do yourid = 1, mpi_nprocs
      xsend(:,:,yourid) = xin(allLonBeg(yourid):allLonEnd(yourid),:) 
    end do
    !$OMP END PARALLEL DO

    nsize = diff(diffID)%lonPerPE_transpose * diff(diffID)%latPerPE
    if (mpi_nprocs > 1) then
      call rpn_comm_alltoall(xsend, nsize, 'mpi_real8',  &
                             xrecv, nsize, 'mpi_real8', 'grid', ierr)
    else
      xrecv(:,:,1) = xsend(:,:,1)
    end if

    !$OMP PARALLEL DO PRIVATE(yourid)
    do yourid = 1, mpi_nprocs
      xout(:,allLatBeg(yourid):allLatEnd(yourid)) = xrecv(:,:,yourid)
    end do
    !$OMP END PARALLEL DO

    deallocate(xsend)
    deallocate(xrecv)

    call tmg_stop(188)

  end subroutine transposeLatToLonBands

  
  subroutine transposeLonToLatBands(diffID, xin, xout)
    implicit none

    ! Arguments:
    integer, intent(in)  :: diffID
    real(8), intent(in)  :: xin ( diff(diffID)%myLonBeg_transpose:diff(diffID)%myLonEnd_transpose, diff(diffID)%nj )
    real(8), intent(out) :: xout( diff(diffID)%ni, diff(diffID)%myLatBeg:diff(diffID)%myLatEnd )

    ! Locals:
    integer :: yourid, ierr, nsize
    integer :: allLonBeg(mpi_nprocs), allLonEnd(mpi_nprocs)
    integer :: allLatBeg(mpi_nprocs), allLatEnd(mpi_nprocs)
    real(8), allocatable :: xsend(:,:,:),xrecv(:,:,:)

    call tmg_start(189,'diff-transposeLonToLat')

    ! Abort if grid cannot be evenly divided over the mpi tasks
    if (diff(diffID)%lonPerPE_transpose /= diff(diffID)%lonPerPEmax_transpose .or. &
        diff(diffID)%latPerPE           /= diff(diffID)%latPerPEmax ) then
      call utl_abort('diffusion_mod-transposeLonToLatBands: grid cannot be evenly distributed over mpi tasks')
    end if

    allocate(xsend(diff(diffID)%lonPerPE_transpose,diff(diffID)%latPerPE, mpi_nprocs))
    allocate(xrecv(diff(diffID)%lonPerPE_transpose,diff(diffID)%latPerPE, mpi_nprocs))

    call rpn_comm_allgather(diff(diffID)%myLonBeg_transpose,1,'mpi_integer',       &
                            allLonBeg                      ,1,'mpi_integer','GRID',ierr)
    call rpn_comm_allgather(diff(diffID)%myLonEnd_transpose,1,'mpi_integer',       &
                            allLonEnd                      ,1,'mpi_integer','GRID',ierr)
    call rpn_comm_allgather(diff(diffID)%myLatBeg,1,'mpi_integer',       &
                            allLatBeg            ,1,'mpi_integer','GRID',ierr)
    call rpn_comm_allgather(diff(diffID)%myLatEnd,1,'mpi_integer',       &
                            allLatEnd            ,1,'mpi_integer','GRID',ierr)

    xout(:,:) = 0.0d0
!    xout(diff(diffID)%myLonBeg_transpose:diff(diffID)%myLonEnd_transpose,diff(diffID)%myLatBeg:diff(diffID)%myLatEnd) = &
!     xin(diff(diffID)%myLonBeg_transpose:diff(diffID)%myLonEnd_transpose,diff(diffID)%myLatBeg:diff(diffID)%myLatEnd)

    !$OMP PARALLEL DO PRIVATE(yourid)
    do yourid = 1, mpi_nprocs
      xsend(:,:,yourid) = xin(:,allLatBeg(yourid):allLatEnd(yourid))
    end do
    !$OMP END PARALLEL DO

    nsize = diff(diffID)%lonPerPE_transpose * diff(diffID)%latPerPE
    if (mpi_nprocs > 1) then
      call rpn_comm_alltoall(xsend, nsize, 'mpi_real8',  &
                             xrecv, nsize, 'mpi_real8', 'grid', ierr)
    else
      xrecv(:,:,1) = xsend(:,:,1)
    end if

    !$OMP PARALLEL DO PRIVATE(yourid)
    do yourid = 1, mpi_nprocs
      xout(allLonBeg(yourid):allLonEnd(yourid),:) = xrecv(:,:,yourid) 
    end do
    !$OMP END PARALLEL DO

    deallocate(xsend)
    deallocate(xrecv)

    call tmg_stop(189)

  end subroutine transposeLonToLatBands

  
  subroutine diffusion1x_implicit(diffID, xin, xout)
    !
    !:Purpose: To compute Lsqrt*xin (diffusion over 1 timestep, loop over
    !          timesteps is external to the subroutine).
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: diffID
    real(8), intent(in)  :: xin ( diff(diffID)%ni, diff(diffID)%myLatBeg:diff(diffID)%myLatEnd )
    real(8), intent(out) :: xout( diff(diffID)%ni, diff(diffID)%myLatBeg:diff(diffID)%myLatEnd )

    ! Locals:
    integer :: latIndex, lonIndex
    real(8) :: xlast( diff(diffID)%ni, diff(diffID)%myLatBeg:diff(diffID)%myLatEnd )
    real(8) :: dp( diff(diffID)%ni )

    call tmg_start(186,'diffusion_implicitx')

    !$OMP PARALLEL DO PRIVATE( latIndex, lonIndex )
    do latIndex = diff(diffID)%myLatBeg, diff(diffID)%myLatEnd
      do lonIndex = 1, diff (diffID) % ni
        xlast ( lonIndex, latIndex ) = xin ( lonIndex, latIndex )
      end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(latIndex,lonIndex,dp)
    do latIndex =  max(2, diff(diffID)%myLatBeg), min(diff(diffID)%nj-1, diff(diffID)%myLatEnd)
      lonIndex = 2
      dp( lonIndex ) = xlast ( lonIndex, latIndex )
      do lonIndex = 3, diff (diffID) % ni - 1
        dp( lonIndex ) = xlast( lonIndex, latIndex ) - diff(diffID) % diff1x_ap( lonIndex, latIndex ) * dp( lonIndex - 1 )
      end do
      lonIndex = diff(diffID) % ni - 1
      xout( lonIndex, latIndex ) = dp( lonIndex ) * diff(diffID) % diff1x_bp_inv( lonIndex, latIndex )
      do lonIndex = diff(diffID) % ni - 2, 2, -1
        xout( lonIndex, latIndex ) = ( dp( lonIndex ) - diff(diffID) % diff1x_c( lonIndex, latIndex ) * xout( lonIndex + 1, latIndex ) ) &
                                   * diff(diffID) % diff1x_bp_inv( lonIndex, latIndex )
      end do
    end do
    !$OMP END PARALLEL DO

    do latIndex = diff(diffID)%myLatBeg, diff(diffID)%myLatEnd
      xout ( 1, latIndex )               = xin ( 1, latIndex )
      xout ( diff(diffID)%ni, latIndex ) = xin ( diff(diffID)%ni, latIndex )
    end do

    call tmg_stop(186)

  end subroutine diffusion1x_implicit


  subroutine diffusion1y_implicit(diffID, xin, xout)
    !
    !:Purpose: To compute Lsqrt*xin (diffusion over 1 timestep, loop over 
    !          timesteps is external to the subroutine).
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: diffID
    real(8), intent(in)  :: xin (diff(diffID)%myLonBeg_transpose:diff(diffID)%myLonEnd_transpose, diff(diffID)%nj)
    real(8), intent(out) :: xout(diff(diffID)%myLonBeg_transpose:diff(diffID)%myLonEnd_transpose, diff(diffID)%nj)

    ! Locals:
    integer :: latIndex, lonIndex
    real(8) :: xlast(diff(diffID)%nj, diff(diffID)%myLonBeg_transpose:diff(diffID)%myLonEnd_transpose)
    real(8) :: dp(diff(diffID)%nj)

    ! NOTE:for improved efficiency, the 2D fields used internally are
    !      ordered (diff(diffID)%nj,diff(diffID)%ni) and
    !      NOT (diff(diffID)%ni,diff(diffID)%nj) as in the rest of the code!

    call tmg_start(187,'diffusion_implicity')

    !$OMP PARALLEL DO PRIVATE(latIndex,lonIndex)
    do latIndex = 1, diff (diffID) % nj
      do lonIndex = diff(diffID)%myLonBeg_transpose, diff(diffID)%myLonEnd_transpose
        xlast ( latIndex, lonIndex ) = xin ( lonIndex, latIndex )
      end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(latIndex,lonIndex,dp)
    do lonIndex = max(2, diff(diffID)%myLonBeg_transpose), min(diff(diffID)%ni-1, diff(diffID)%myLonEnd_transpose)
      latIndex = 2
      dp ( latIndex ) = xlast ( latIndex, lonIndex )
      do latIndex = 3, diff (diffID) % nj - 1
        dp ( latIndex ) = xlast ( latIndex, lonIndex ) - diff (diffID) % diff1y_ap ( latIndex, lonIndex ) * dp ( latIndex - 1 )
      end do
      latIndex = diff (diffID) % nj - 1
      xout ( lonIndex, latIndex ) = dp ( latIndex ) * diff (diffID) % diff1y_bp_inv ( latIndex, lonIndex )
      do latIndex = diff(diffID) % nj - 2, 2, -1
        xout ( lonIndex, latIndex ) = ( dp ( latIndex ) - diff(diffID) % diff1y_c ( latIndex, lonIndex ) * xout ( lonIndex, latIndex + 1 ) ) &
                                    * diff (diffID) % diff1y_bp_inv ( latIndex, lonIndex )
      end do
    end do
    !$OMP END PARALLEL DO

    do lonIndex = diff(diffID)%myLonBeg_transpose, diff(diffID)%myLonEnd_transpose
      xout( lonIndex, 1 )                    = xin ( lonIndex, 1 )
      xout( lonIndex, diff (diffID) % nj ) = xin ( lonIndex, diff (diffID) % nj )
    end do

    call tmg_stop(187)

  end subroutine diffusion1y_implicit

end module diffusion_mod
