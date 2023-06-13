
MODULE localizationFunction_mod
  ! MODULE localizationFunction_mod (prefix='lfn' category='2. B and R matrices')
  !
  ! :Purpose: To store various functions for horizontal and vertical covariance
  !           localization. Length-scale estimation for function available in
  !           this module is also possible through curve fitting using
  !           pre-computed optimal separation-distance localization values.
  !
  use earthConstants_mod
  use utilities_mod
  implicit none
  save
  private

  ! public procedures
  public :: lfn_Setup, lfn_Response, lfn_lengthscale, lfn_createBiPerFunction

  logical             :: initialized = .false.

  character(len=128)  :: LocFunction

CONTAINS

  !--------------------------------------------------------------------------
  ! lfn_Setup
  !--------------------------------------------------------------------------
  subroutine lfn_Setup(LocFunctionWanted)
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: LocFunctionWanted
  
    select case(trim(LocFunctionWanted))
    case ('FifthOrder')
       write(*,*)
       write(*,*) 'lfn_setup: Using Gaspari and Cohn 5th order piecewise rationale function'
       LocFunction='FifthOrder'
    case default
       write(*,*)
       write(*,*) 'Unsupported HORIZONTAL localization function : ', trim(LocFunctionWanted)
       call utl_abort('lfn_setup')
    end select

    initialized  = .true.

  end subroutine lfn_Setup

  !--------------------------------------------------------------------------
  ! lfn_Response
  !--------------------------------------------------------------------------
  function lfn_Response(distance,lengthscale) result(correlation)
    implicit none

    ! Arguments:
    real(8) :: distance
    real(8) :: lengthscale
    ! Result:
    real(8) :: correlation

    if ( .not. initialized ) then
      write(*,*)
      write(*,*) 'The localixation function module was NOT initialized'
      call utl_abort('lfn_horizResponse')
    endif

    select case(trim(LocFunction))
    case ('FifthOrder')
       correlation=lfn_FifthOrderFunction(distance,lengthscale)
    case default
       write(*,*)
       write(*,*) 'Unsupported localization function : ', trim(LocFunction)
       call utl_abort('lfn_Response')
    end select

  end function lfn_Response

  !--------------------------------------------------------------------------
  ! lfn_gradient
  !--------------------------------------------------------------------------
  function lfn_gradient(distance,lengthscale) result(gradient)
    implicit none

    ! Arguments:
    real(8) :: distance
    real(8) :: lengthscale
    ! Result:
    real(8) :: gradient

    if ( .not. initialized ) then
      write(*,*)
      write(*,*) 'The localixation function module was NOT initialized'
      call utl_abort('lfn_gradient')
    endif

    select case(trim(LocFunction))
    case ('FifthOrder')
       gradient=lfn_FifthOrderFunctionGradient(distance,lengthscale)
    case default
       write(*,*)
       write(*,*) 'Unsupported localization function : ', trim(LocFunction)
       call utl_abort('lfn_gradient')
    end select

  end function lfn_gradient

  !--------------------------------------------------------------------------
  ! lfn_FifthOrderFunction
  !--------------------------------------------------------------------------
  function lfn_FifthOrderFunction(distance,lengthscale) result(correlation)
    implicit none

    ! Arguments:
    real(8) :: distance
    real(8) :: lengthscale
    ! Result:
    real(8) :: correlation

    ! Locals:
    real(8) :: halflength

    halflength = lengthscale/2.d0

    if ( distance <= halflength ) then
       correlation = -0.250d0*(distance/halflength)**5  &
                     +         0.5d0*(distance/halflength)**4  &
                     +       0.625d0*(distance/halflength)**3  &
                     - (5.0d0/3.0d0)*(distance/halflength)**2  &
                     + 1.0d0
    else if ( distance <= (2.0d0*halflength) ) then
       correlation =  (1.0d0/12.0d0)*(distance/halflength)**5  &
                     -         0.5d0*(distance/halflength)**4  &
                     +       0.625d0*(distance/halflength)**3  &
                     + (5.0d0/3.0d0)*(distance/halflength)**2  &
                     -         5.0d0*(distance/halflength)     &
                     + 4.0d0                                 &
                     - (2.0d0/3.0d0)*(halflength/distance) 
    else
       correlation = 0.d0
    endif

  end function lfn_FifthOrderFunction

  !--------------------------------------------------------------------------
  ! lfn_FifthOrderFunctionGradient
  !--------------------------------------------------------------------------
  function lfn_FifthOrderFunctionGradient(distance,lengthscale) result(gradient)
    implicit none

    ! Arguments:
    real(8) :: distance
    real(8) :: lengthscale
    ! Result:
    real(8) :: gradient

    ! Locals:
    real(8) :: halflength

    halflength = lengthscale/2.d0

    ! Compute df/dhalflength
    if ( distance <= halflength ) then
       gradient =     5.d0*0.250d0*(distance/halflength)**5 / halflength  &
                  -     4.d0*0.5d0*(distance/halflength)**4 / halflength  &
                  -   3.d0*0.625d0*(distance/halflength)**3 / halflength  &
                  + (10.0d0/3.0d0)*(distance/halflength)**2 / halflength
    else if ( distance <= (2.0d0*halflength) ) then
       gradient =  (-5.0d0/12.0d0)*(distance/halflength)**5 / halflength  &
                  +    4.0d0*0.5d0*(distance/halflength)**4 / halflength  &
                  -  3.0d0*0.625d0*(distance/halflength)**3 / halflength  &
                  - (10.0d0/3.0d0)*(distance/halflength)**2 / halflength  &
                  +          5.0d0*(distance/halflength)    / halflength  &
                  -   (2.0d0/3.0d0)*(1.0d0/distance) 
    else
       gradient = 0.d0
    endif

    ! Compute df/dlengthscale (= df/dhalflength * dhalflength/dlengthscale)
    gradient = gradient / 2.d0

  end function lfn_FifthOrderFunctionGradient

!--------------------------------------------------------------------------
! lfn_CreateBiPerFunction
!--------------------------------------------------------------------------
  SUBROUTINE  lfn_CreateBiPerFunction(gridpoint,CorrelLength,dlon,ni,nj,nk)
    implicit none

    ! Arguments:
    integer, intent(in)  :: ni
    integer, intent(in)  :: nj
    integer, intent(in)  :: nk
    real(8), intent(in)  :: dlon             ! in radian
    real(8), intent(in)  :: CorrelLength(nk) ! in km
    real(8), intent(out) :: gridpoint(ni,nj,nk)

    ! Locals:
    integer          :: i, j, k, iref, jref
    real(8)          :: distance, distance_ref

    gridpoint(:,:,:) = 0.d0

    distance_ref = dlon * ec_ra

    !
    !- Create a bi-periodic correlation function by centering the function in each 4 corners
    !

    !- Lower-Left Corner (reference corner)
    iref = 1
    jref = 1
    do j = 1, nj
       do i = 1, ni
          distance = distance_ref * sqrt( real((i-iref)**2 + (j-jref)**2,8) )
          do k = 1, nk
             gridpoint(i,j,k) = gridpoint(i,j,k) + lfn_response(distance,1000.d0*CorrelLength(k)) 
          enddo
       enddo
    enddo

    !- Upper-Left Corner
    iref = 1
    jref = nj+1 ! the lam spectral transform is periodic at nj+1
    do j = 1, nj
       do i = 1, ni
          distance = distance_ref * sqrt( real((i-iref)**2 + (j-jref)**2,8) )
          do k = 1, nk
             gridpoint(i,j,k) = gridpoint(i,j,k) + lfn_response(distance,1000.d0*CorrelLength(k)) 
          enddo
       enddo
    enddo

    !- Lower-Right Corner
    iref = ni+1 ! the lam spectral transform is periodic at ni+1
    jref = 1
    do j = 1, nj
       do i = 1, ni
          distance = distance_ref * sqrt( real((i-iref)**2 + (j-jref)**2,8) )
          do k = 1, nk
             gridpoint(i,j,k) = gridpoint(i,j,k) + lfn_response(distance,1000.d0*CorrelLength(k)) 
          enddo
       enddo
    enddo

    !- Upper-Right Corner
    iref = ni+1 ! the lam spectral transform is periodic at ni+1
    jref = nj+1 ! the lam spectral transform is periodic at nj+1
    do j = 1, nj
       do i = 1, ni
          distance = distance_ref * sqrt( real((i-iref)**2 + (j-jref)**2,8) )
          do k = 1, nk
             gridpoint(i,j,k) = gridpoint(i,j,k) + lfn_response(distance,1000.d0*CorrelLength(k)) 
          enddo
       enddo
    enddo

  END SUBROUTINE Lfn_CreateBiPerFunction

  !--------------------------------------------------------------------------
  ! LFN_LENGTHSCALE
  !--------------------------------------------------------------------------
  subroutine lfn_lengthscale(lengthscale, rmse, localizationValues, &
                             distance, weight, numbins)
    implicit none

    ! Arguments:
    integer, intent(in)    :: numBins
    real(8), intent(in)    :: localizationValues(numBins)
    real(8), intent(in)    :: distance(numBins)
    real(8), intent(in)    :: weight(numBins)
    real(8), intent(inout) :: lengthscale
    real(8), intent(out)   :: rmse

    ! Locals:
    integer :: ierr, ilist
    real(8) :: minv, fit

    ! 1.  Initialization
    ilist       = 1 
    ierr        = 0
    minv        = 0.01

    ! 2.  Apply Curve fitting algorithm
    call lfn_curveFit(numBins, distance,      & ! IN
                      localizationValues,     & ! IN
                      weight, numBins,        & ! IN
                      lengthscale,            & ! INOUT
                      minv, ilist,            & ! IN
                      ierr,                   & ! INOUT
                      fit)                      ! OUT

    ! 3.  root mean squared deviation output as 'degree of fit'
    rmse = sqrt(fit)

  end subroutine lfn_lengthscale

  !--------------------------------------------------------------------------
  ! LFN_CURVEFIT
  !--------------------------------------------------------------------------
  subroutine lfn_curveFit(Nn, x, y, w, nmax, param, minv, ilist, ierr, ssq)
    !
    !:Purpose: This subroutine computes the lengthscale (sl, here: param) of a
    !          FUNCTION that best fits (in a least-square perspective) a
    !          data set. The method follows routine CURVEFIT from:
    !          Heeswijk, M.V., and C.G. Fox, 1988: Iterative Method and Fortran
    !          Code for Nonlinear Curve Fitting, Computers and Geosciences, 14,
    !          4, pp. 489-503.
    !
    implicit none

    ! Arguments:
    INTEGER,            INTENT(IN)    :: Nn
    REAL(8),            INTENT(IN)    :: x(Nn)
    REAL(8),            INTENT(IN)    :: y(Nn)
    REAL(8),            INTENT(IN)    :: w(Nn)
    INTEGER,            INTENT(IN)    :: nmax
    REAL(8),            INTENT(INOUT) :: param
    REAL(8),            INTENT(IN)    :: minv
    INTEGER,            INTENT(IN)    :: ilist
    INTEGER,            INTENT(INOUT) :: ierr
    REAL(8),            INTENT(OUT)   :: ssq

    ! Locals:
    INTEGER, PARAMETER    :: icon = 100   ! max iteration
    INTEGER, PARAMETER    :: lbis = 10    ! max bisection    
    REAL(8)              :: d(Nn)
    REAL(8)              :: g(Nn)
    REAL(8)              :: m
    REAL(8)              :: gtd
    REAL(8)              :: gtg
    REAL(8)              :: gtginv
    REAL(8)              :: ssqold, rootmsq
    REAL(8)              :: w_t    
    INTEGER           :: iter,nbis,i                     ! Loop counters.
    CHARACTER(LEN=60) :: FMT1
    CHARACTER(LEN=12), PARAMETER :: RoutineName = 'lfn_curveFit'
    
    !-------------------------------------------------------------------------------
    
    !-------------------------------------------------------------------------------
    ! 1. Initialize variables
    !-------------------------------------------------------------------------------
    ierr   = 0          ! Error flags
    iter   = 0          ! Iteration counter
    nbis   = 0          ! Bisection counter
    m      = 9.9E9      ! dummy value to initialize m
    ssqold = 9.9E9      ! dummy value to initialize ssqold
    
    !-------------------------------------------------------------------------------
    ! 2. Perform curve fitting using a minimization (least-square) approach
    !-------------------------------------------------------------------------------
    !PRINT*, RoutineName, ' fitting ', trim(locFunction)
    
    DO WHILE ( (ABS(m) > minv) .AND. (iter <= icon) .AND. (nbis <= lbis) )
       
       ! 2.1  Form vector D; i.e. a data point minus the isolated Taylor Series
       !      term for that data point. In addition calculate the sum of the 
       !      square of the residuals
       ssq = 0.d0
       DO i = 1, nmax
          d(i) = y(i) - lfn_response(x(i), param)
          ssq = ssq + d(i)**2 * w(i)
       END DO
       rootmsq = sqrt(ssq)
       
       ! 2.2  If we are converging
       IF (ssq < ssqold) THEN
          
          ! The matrix formation in the following program segment follows
          ! equation 3.12 on page 41 of Geophysical Data Analysis: Discrete
          ! Inverse Theory by Menke (1984)
          
          ! 2.2.1  Form matrix G
          DO i = 1, nmax
             g(i)=lfn_gradient(x(i), param) ! curve derivative wrt to param
          END DO
          
          ! 2.2.2  Form matrix GTD
          gtd = 0.d0
          DO i = 1, nmax
             gtd = gtd + d(i) * g(i) * w(i)
          END DO
          
          ! 2.2.3  Form matrix GTG
          gtg = 0.d0
          DO i = 1, nmax
             gtg = gtg + g(i) * g(i) * w(i)
          END DO
          
          ! 2.2.4  Find the inverse of matrix GTG
          if (gtg /= 0.d0) then
             gtginv = 1.d0 / gtg
          else
             write(*,*) 'Lfn_Curvefit: gtg = 0 ! Aborting the iteration process'
             DO i = 1, nmax
                write(*,*) g(i), w(i), x(i), param 
             END DO
             gtginv = tiny(1.d0)
             !call utl_abort('Lfn_Curvefit')
          end if
          
          ! 2.2.5  Find vector M (i.e. increment of PARAM)
          m = gtginv * gtd
          
          ssqold = ssq          ! keep
          
          iter = iter + 1
          
          ! 2.3  If we are diverging
       ELSE 
          
          m = m / 2.d0           ! Bisect the correction
          nbis = nbis + 1
          
       END IF
       
       ! 2.4  Update the parameter
       param = param + m
       
       ! 2.5  Print some info
       FMT1 = "(A12,2X,A2,I6,2X,A2,I3,2X,A2,G10.4,2X,A5,G10.4,2X,A2,G10.4)"
       
       IF (ilist == 1) WRITE(*,FMT1) RoutineName,'i=',iter,'n=',nbis,'m=',m,  &
            'rmsd=',rootmsq,'s=',param
       
    END DO
    
    PRINT*
    PRINT*," Estimated lengthscale = ",param
    
    !-------------------------------------------------------------------------------
    ! 3. Compute the final mean squared deviation
    !-------------------------------------------------------------------------------
    ssq = 0.d0
    w_t = 0.d0
    DO i = 1, nmax
       d(i) = y(i) - lfn_response(x(i), param)
       w_t = w_t + w(i)
       ssq = ssq + w(i) * d(i)**2
    END DO
    ssq = ssq / w_t
    
  end SUBROUTINE Lfn_Curvefit

END MODULE localizationFunction_mod
