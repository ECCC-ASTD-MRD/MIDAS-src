
module linearAlgebra_mod
  !
  ! MODULE midasMpi_mod (prefix='linalg' category='8. Low-level utilities and constants')
  !
  !:Purpose:  Subroutines/functions and public variables related to linear algebra routines
  !
#if MKL_SUPPORT
  use mkl_service
#endif
  use utilities_mod

  implicit none
  save
  private

  ! Public variables
  !integer, public, protected :: mmpi_myid      = 0

  ! module constants
  ! global communicator
  ! Since, we are only considering a single grid, we can assume here that there is only one world
 ! type(mpi_comm), public, parameter :: mmpi_comm_GRID = MPI_COMM_WORLD

  ! Public procedures
  public :: linalg_setMKLThreads
  public :: linalg_matSqrt
  public :: linalg_fastMatMul
  public :: linalg_matInverse
  public :: linalg_eigenDecomp
  public :: linalg_fastInverse
  public :: linalg_pseudo_inverse

contains

  subroutine linalg_setMKLThreads()

#if MKL_SUPPORT
    integer :: nulnam, ierr
    ! external definitions
    integer, external :: fnom, fclos

    ! Namelist variables for 'namMKL'
    logical :: oneThreadMKL ! choose to use only 1 thread for MKL subroutines
    logical :: dynamicMKL   ! choose to use dynamic assignment of threads for MKL subroutines

    namelist /nammkl/ oneThreadMKL, dynamicMKL

    ! default values for MKL namelist
    oneThreadMKL = .false.
    dynamicMKL = .true.

    ! read the MKL namelist
    if (utl_isNamelistPresent('namMKL','./flnml')) then
       ! reading namelist variables
       nulnam = 0
       ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
       read(nulnam, nml = namMKL, iostat = ierr)
       if (ierr /= 0) call utl_abort('linalg_setMKLThreads:: Error reading namelist')
       ierr = fclos(nulnam)
    else
      write(*,*) 'linalg_setMKLThreads: namMKL is missing in the namelist.'
      write(*,*) '                      the default values will be taken.'
    endif

    ! Modify the MKL thread configuration based on namelist variables

    if (dynamicMKL) then
      call mkl_set_dynamic(1)
    else
      call mkl_set_dynamic(0)
    end if

    if (oneThreadMKL) then
      call mkl_set_num_threads(1)
      write(*,*) 'linalg_setMKLThreads: number of threads used for MKL set to one'
    else
      write(*,*) 'linalg_setMKLThreads: default number of threads used for MKL'
    end if
#else
    if (utl_isNamelistPresent('namMKL','./flnml')) then
      write(*,*) 'linalg_setMKLThreads: Ignoring ''namMKL'' namelist because there is no MKL support'
    end if
#endif

  end subroutine linalg_setMKLThreads

  subroutine linalg_matsqrt(matrix, rank, exponentSign, printInformation_opt )
    !
    !:Purpose: Calculate square root of an error covariance matrix
    !
    implicit none

    ! Arguments:
    integer,           intent(in)    :: rank
    real(8),           intent(inout) :: matrix(rank,rank)
    real(8),           intent(in)    :: exponentSign
    logical, optional, intent(in)    :: printInformation_opt ! switch to print be more verbose

    ! Locals:
    real(8), allocatable :: eigenValues(:)
    real(8), allocatable :: work(:)
    real(8), allocatable :: eigenVectors(:,:)
    integer :: sizework, info, index, index1, index2
    logical :: printInformation

    if (present(printInformation_opt)) then
       printInformation = printInformation_opt
    else
       printInformation = .false.
    end if

    if (printInformation) then
      write(*,*)
      write(*,*) 'linalg_matsqrt: Starting...'
    end if

    sizework = 64 * rank
    allocate(work(sizework))

    allocate(eigenValues (rank))
    allocate(eigenVectors(rank,rank))

    !- Calculate EigenVectors (V) and EigenValues (D) of B matrix
    eigenVectors(:,:) = matrix(:,:)

    call dsyev('V','U',rank,   & ! IN
               eigenVectors,   & ! INOUT
               rank,           & ! IN
               eigenValues,    & ! OUT
               work, sizework, & ! IN
               info )            ! OUT

    if ( info /= 0 ) then
      write(*,*)
      write(*,*) 'dsyev: ',info
      call utl_abort('linalg_matsqrt: DSYEV failed!')
    end if

    if (printInformation) then
      write(*,*)
      write(*,'(1x,"Original EIGEN VALUES: ")')
      write(*,'(1x,10f7.3)') (eigenValues(index),index=1,rank)
      if (exponentSign < 0.d0) then
        write(*,*)
        write(*,'(A,1x,e14.6)') "Condition number:", &
             maxval(eigenValues(:))/minval(eigenValues(:))
      end if
    end if

    !- Calculate Matrix^0.5 = V D^0.5 V^t
    where(eigenValues(:) < 0.d0)
      eigenValues = 0.d0
    end where

    eigenValues(:) = eigenValues(:)**(0.5d0*exponentSign)

    do index1 = 1, rank
      do index2 = 1, rank
        matrix(index1,index2) = sum ( eigenVectors (index1,1:rank)   &
                                    * eigenVectors (index2,1:rank)   &
                                    * eigenValues(1:rank) )
      end do
    end do

    deallocate(eigenVectors)
    deallocate(eigenValues)
    deallocate(work)

    if (printInformation) then
      write(*,*)
      write(*,*) 'linalg_matsqrt: Ending...'
    end if

  end subroutine linalg_matsqrt

  !--------------------------------------------------------------------------
  ! linalg_matInverse
  !--------------------------------------------------------------------------
  subroutine linalg_matInverse(matrix, rank, inverseSqrt_opt, printInformation_opt, &
                            eigenValueRelThreshold_opt)
    !
    !:Purpose: Calculate the inverse of a covariance matrix
    !          and, optionally, also the inverse square-root.
    !
    implicit none

    ! Arguments:
    integer,           intent(in)    :: rank                 ! order of the matrix
    real(8),           intent(inout) :: matrix(:,:)          ! on entry, the original matrix; on exit, the inverse
    real(8), optional, intent(inout) :: inverseSqrt_opt(:,:) ! if present, the inverse sqrt matrix on exit
    real(8), optional, intent(in)    :: eigenValueRelThreshold_opt
    logical, optional, intent(in)    :: printInformation_opt ! switch to print be more verbose

    ! Locals:
    integer :: index1, index2, info, sizework
    real(8) :: sizework_r8(1), eigenValueMin
    real(8), allocatable :: work(:), eigenVectors(:,:), eigenValues(:)
    logical :: printInformation

    call utl_tmg_start(180,'low-level--linalg_matInverse')

    if (present(printInformation_opt)) then
      printInformation = printInformation_opt
    else
      printInformation = .false.
    end if

    if (printInformation) then
      write(*,*)' linalg_matInvers: Inverse matrix of a symmetric matrix'
    end if

    ! 1. Computation of eigenvalues and eigenvectors

    allocate(eigenVectors(rank,rank))
    allocate(eigenValues(rank))

    do index2=1,rank
      do index1=1,rank
        eigenVectors(index1,index2)=matrix(index1,index2)
      end do
    end do

    ! query the size of the 'work' vector by calling 'DSYEV' with 'sizework=-1'
    sizework = -1
    info = -1
    call dsyev('V','U',rank, eigenVectors, rank, eigenValues, sizework_r8, sizework, info)

    ! compute the eigenvalues
    sizework=int(sizework_r8(1))
    allocate(work(sizework))
    call dsyev('V','U',rank, eigenVectors,rank, eigenValues, work, sizework, info)
    deallocate(work)

    if (printInformation) then
      write(*,'(1x,"Original eigen values: ")')
      write(*,'(1x,10f7.3)') (eigenValues(index1),index1=1,rank)

      if(minval(eigenValues) > 1.0d-10) then
        write(*,'(A,1x,e14.6)') "Condition number:", &
             maxval(eigenValues)/minval(eigenValues)
      end if
    end if

    ! 2.  Take inverse of eigenvalues

    if (present(eigenValueRelThreshold_opt)) then
      eigenValueMin = eigenValueRelThreshold_opt*maxval(eigenValues(1:rank))
    else
      eigenValueMin = 1.0d-10
    end if

    do index1=1,rank
      if(eigenValues(index1) > eigenValueMin) then
        eigenValues(index1)= 1.0d0/eigenValues(index1)
      else
        write(*,*) 'linalg_matInverse: WARNING eigenvalue is too small = ', index1, eigenValues(index1)
        eigenValues(index1) = 0.0d0
      end if
    end do

    if (printInformation) then
      write(*,'(1x,"Inverse of original eigen values: ")')
      write(*,'(1x,10f7.3)') (eigenValues(index1),index1=1,rank)
    end if

    ! 3.  Compute the inverse matrix

    do index2 = 1, rank
      do index1 = 1, rank
        matrix(index1,index2) = sum ( eigenVectors (index1,1:rank)   &
                                    * eigenVectors (index2,1:rank)   &
                                    * eigenValues(1:rank) )
      end do
    end do

    ! 4.  If requested, computed the inverse square-root also

    if (present(inverseSqrt_opt)) then
      do index1=1,rank
        if(eigenValues(index1) > 1.0d-10) then
          eigenValues(index1)= sqrt(eigenValues(index1))
        else
          eigenValues(index1) = 0.0d0
        end if
      end do
      do index2 = 1, rank
        do index1 = 1, rank
          inverseSqrt_opt(index1,index2) = sum ( eigenVectors (index1,1:rank)   &
                                               * eigenVectors (index2,1:rank)   &
                                               * eigenValues(1:rank) )
        end do
      end do
    end if

    ! 5. Deallocate local arrays
    deallocate(eigenVectors,eigenValues)

    if (printInformation) then
      write(*,*) 'linalg_matInverse: done'
      write(*,*) ' '
    end if

    call utl_tmg_stop(180)

  end subroutine linalg_matInverse

  !--------------------------------------------------------------------------
  ! linalg_eigenDecomp
  !--------------------------------------------------------------------------
  subroutine linalg_eigenDecomp(matrix, eigenValues, eigenVectors, tolerance, numReturned, printInformation_opt)
    !
    !:Purpose: Calculate eigenValues/Vectors and return only those with eigenValues
    !          whose magnitude is greater than the specified tolerance.
    !
    implicit none

    ! Arguments:
    real(8),           intent(in)  :: matrix(:,:)          ! on entry, the original matrix; on exit, the inverse
    real(8),           intent(out) :: eigenValues(:)       ! computed eigenValues
    real(8),           intent(out) :: eigenVectors(:,:)    ! computed eigenVectors
    real(8),           intent(in)  :: tolerance            ! threshold for eigenValue magnitude to be returned
    integer,           intent(out) :: numReturned          ! number of eigenValues/Vectors returned
    logical, optional, intent(in)  :: printInformation_opt ! switch to print be more verbose

    ! Locals:
    integer :: rank, index1, index2, info, sizework
    real(8) :: sizework_r8(1)
    real(8), allocatable :: work(:), eigenVectorsOrig(:,:), eigenValuesOrig(:)
    logical :: printInformation

    if (present(printInformation_opt)) then
      printInformation = printInformation_opt
    else
      printInformation = .false.
    end if

    if (printInformation) then
      write(*,*)' linalg_eigenDecomp: Eigen decomposition of a symmetric matrix'
    end if

    !     1. Computation of eigenvalues and eigenvectors

    rank = size(matrix,1)
    allocate(eigenVectorsOrig(rank,rank))
    allocate(eigenValuesOrig(rank))

    do index2 = 1, rank
      do index1 = 1, rank
        eigenVectorsOrig(index1,index2)=matrix(index1,index2)
      end do
    end do

    ! Query the size of the 'work' vector by calling 'DSYEV' with 'sizework=-1'
    sizework = -1
    info = -1
    call dsyev('V', 'U', rank, eigenVectorsOrig, rank, eigenValuesOrig,  &
               sizework_r8, sizework, info)

    ! Compute the eigenvalues/vectors
    sizework = int(sizework_r8(1))
    allocate(work(sizework))
    call dsyev('V', 'U', rank, eigenVectorsOrig, rank, eigenValuesOrig,  &
               work, sizework, info)
    deallocate(work)

    if (printInformation) then
      write(*,'(1x,"Original eigen values: ")')
      write(*,'(1x,10f7.3)') (eigenValuesOrig(index1),index1=1,rank)

      if(minval(eigenValuesOrig) > tolerance) then
        write(*,'(A,1x,e14.6)') "Condition number:", &
             maxval(eigenValuesOrig)/minval(eigenValuesOrig)
      end if
    end if

    !     2.  Determine which eigen values/vectors to return

    numReturned = 0
    do index1 = rank, 1, -1
      if (eigenValuesOrig(index1) > tolerance) then
        numReturned = numReturned + 1
      else
        exit
      end if
    end do

    if (printInformation) then
      write(*,*) 'Number of eigen values returned =', numReturned, ' out of', rank
    end if

    !     3.  Copy eigenValues/Vectors into output arrays with reversed order

    do index1 = 1, numReturned
      ! And set negative values to zero
      eigenValues(index1) = max(0.0D0,eigenValuesOrig(rank-index1+1))
    end do
    do index1 = numReturned+1, rank
      eigenValues(index1) = 0.0D0
    end do

    do index2 = 1, numReturned
      do index1 = 1, rank
        eigenVectors(index1,index2) = eigenVectorsOrig(index1,rank-index2+1)
      end do
    end do
    do index2 = numReturned+1, rank
      do index1 = 1, rank
        eigenVectors(index1,index2) = 0.0D0
      end do
    end do

    !     4. Deallocate local arrays
    deallocate(eigenVectorsOrig,eigenValuesOrig)

    if (printInformation) then
      write(*,*) 'linalg_eigenDecomp: done'
      write(*,*) ' '
    end if

  end subroutine linalg_eigenDecomp

  !-----------------------------------------
  ! linalg_fastInverse
  !-----------------------------------------
  subroutine linalg_fastInverse(inputMatrix, inverse)
    !
    !:Purpose: for fast computation of matrix inverse using LU decomposition
    !
    implicit none

    ! Arguments:
    real(8),           intent(in)  :: inputMatrix(:,:)   ! Input Matrix
    real(8),           intent(out) :: inverse(:,:)       ! its inverse

    ! Locals:
    integer :: info, lwork, lineDim, columnDim
    real(8), allocatable :: work(:)
    integer, allocatable :: pivot(:)

    lineDim = size(inputMatrix, dim=1)
    columnDim = size(inputMatrix, dim=2)
    if (lineDim /= columnDim) then
      call utl_abort('linalg_fastInverse: the input matrix should be square !')
    end if

    inverse(1:lineDim,1:columnDim) = inputMatrix(:,:)

    allocate(pivot(lineDim))
    call dgetrf(lineDim, lineDim, inverse, columnDim, pivot, info)
    if (info < 0) then
      call utl_abort('linalg_fastInverse: invalid value for parameter ' // utl_str(-info) // ' in lapack subroutine dgetrf')
    else if (info > 0) then
      call utl_abort('linalg_fastInverse: in dgetrf the U matrix is exactly singular ' // utl_str(info) )
    end if

    lwork = -1
    allocate(work(1))
    !first call to query work array work size
    call dgetri(columnDim, inverse, columnDim, pivot, work, lwork, info)
    lwork = int(work(1))
    call utl_reallocate(work, lwork)
    call dgetri(columnDim, inverse, columnDim, pivot, work, lwork, info)
    if (info < 0) then
      call utl_abort('linalg_fastInverse: invalid value for parameter ' // utl_str(-info) // ' in lapack subroutine dgetri')
    else if (info > 0) then
      call utl_abort('linalg_fastInverse: in dgetri singular matrix' // utl_str(info) )
    end if

    deallocate(work)
    deallocate(pivot)

  end subroutine linalg_fastInverse

  !-----------------------------------------
  ! linalg_pseudo_inverse
  !-----------------------------------------
  subroutine linalg_pseudo_inverse(inputMatrix, pseudoInverse, threshold_opt)
    !
    !:Purpose: to calculate the More-Penrose pseudo inverse of the matrix inputMatrix
    !
    implicit none

    ! Arguments:
    real(8),           intent(in)  :: inputMatrix(:,:)   ! Input Matrix
    real(8),           intent(out) :: pseudoInverse(:,:) ! its Moore Penrose Pseudo-Inverse
    real(8), optional, intent(in)  :: threshold_opt

    ! Locals:
    real(8), allocatable :: copyMatrix(:,:), leftSingularVector(:,:), rightSingularVectorT(:,:)
    real(8), allocatable :: singularValues(:)
    integer :: info, lwork, lineIndex
    integer :: lineDim, columnDim, minDim
    real(8) :: thresh
    real(8), allocatable :: work(:)
    character(len=80) :: errorMessage

    lineDim = size(inputMatrix, dim=1)
    columnDim = size(inputMatrix, dim=2)
    minDim = min(lineDim, columnDim)

    allocate( copyMatrix(lineDim,columnDim), leftSingularVector(lineDim,lineDim) )
    allocate( rightSingularVectorT(columnDim,columnDim), singularValues(minDim) )

    copyMatrix(:,:) = inputMatrix(:,:) ! Work with a copy because copyMatrix will be overwriten
    lwork = max(10000, max(1, 3 * min(lineDim,columnDim) + max(lineDim,columnDim), 5 * minDim ))
    allocate(work(lwork))
    call dgesvd("A", "A", lineDim, columnDim, copyMatrix, lineDim, singularValues, &
         leftSingularVector, lineDim, rightSingularVectorT, columnDim, work, lwork, info )

    if (info /= 0) then
      write(errorMessage,*) "linalg_pseudo_inverse: Problem in DGESVD ! ",info
      call utl_abort(errorMessage)
    end if

    deallocate(work)

    if (present(threshold_opt)) then
      thresh = threshold_opt
    else
      !according to wikipedia... as in matlab or numpy
      thresh = epsilon(thresh) * max(lineDim, columnDim) * maxval(singularValues)
    end if
    print *,"linalg_pseudo_inverse: threshold= ",thresh

    pseudoInverse(:,:)=0.d0
    do lineIndex = 1, minDim
      If (singularValues(lineIndex) > thresh) then
        pseudoInverse(lineIndex,:) = ( 1.d0 / singularValues(lineIndex) ) * leftSingularVector(:,lineIndex)
      end if
    end do

    pseudoInverse = matmul( transpose(rightSingularVectorT), pseudoInverse)

    deallocate( singularValues, rightSingularVectorT )
    deallocate( leftSingularVector, copyMatrix )

  end subroutine linalg_pseudo_inverse

  subroutine linalg_fastMatMul(AmatrixIn, BmatrixIn, CmatrixOut, &
                            isATransposed_opt, isBTransposed_opt, isCSymmetric_opt, &
                            firstDim_opt, lastDim_opt, summationDim_opt)
    !
    !:Purpose: Compute matrix multiplication CmatrixOut=AmatrixIn*BmatrixIn
    !             AmatrixIn  is a matrix MxK
    !             BmatrixIn  is a matrix KxN
    !             CmatrixOut is a matrix MxN
    !
    !          The matrix dimensions (M, N and K) are usually inferred
    !          from argument array allocated dimensions but the arrays
    !          'AmatrixIn', 'BmatrixIn' and 'CmatrixOut' can be larger
    !          allocated arrays than 'M', 'N' and 'K' in which case
    !          you have to specify those numbers with arguments:
    !               firstDim_opt     for M
    !               lastDim_opt      for N
    !               summationDim_opt for K
    !
    !          If the result matrix is expected to be symmetric, you
    !          can avoid to do the full matrix multiplication with
    !          'isCSymmeric_opt = .true.' and compute only the upper
    !          half of the result matrix.  The lower half will be
    !          copied from the just computed upper half.
    implicit none

    ! Arguments:
    real(8),           intent(in)  :: AmatrixIn(:,:)  ! Input  matrix
    real(8),           intent(in)  :: BmatrixIn(:,:)  ! Input  matrix
    real(8),           intent(out) :: CmatrixOut(:,:) ! Output matrix
    logical, optional, intent(in)  :: isATransposed_opt ! Should the matrix 'AmatrixIn' be transposed before multiplication?
    logical, optional, intent(in)  :: isBTransposed_opt ! Should the matrix 'BmatrixIn' be transposed before multiplication?
    logical, optional, intent(in)  :: isCSymmetric_opt  ! Is the result matrix 'CmatrixOut' expected to be symmetric?
    integer, optional, intent(in)  :: firstDim_opt     ! First  dimension of 'AmatrixIn' or 'CmatrixOut' (defaults to first  dimension of 'CmatrixOut')
    integer, optional, intent(in)  :: lastDim_opt      ! Second dimension of 'BmatrixIn' or 'CmatrixOut' (defaults to second dimension of 'CmatrixOut')
    integer, optional, intent(in)  :: summationDim_opt ! Second dimension of 'AmatrixIn' or first dimension of 'BmatrixIn' (defaults to second dimension of 'AmatrixIn')

    ! Locals:
    integer :: firstDim, lastDim, summationDim
    integer :: firstDimA, firstDimB, firstDimC
    integer :: iIndex, jIndex
    character :: transposeA, transposeB
    logical :: isCSymmetric

    firstDimA = size(AmatrixIn, 1)
    firstDimB = size(BmatrixIn, 1)
    firstDimC = size(CmatrixOut,1)

    if (present(firstDim_opt)) then
      firstDim = firstDim_opt
    else
      firstDim = size(CmatrixOut,1)
    end if

    if (present(lastDim_opt)) then
      lastDim = lastDim_opt
    else
      lastDim = size(CmatrixOut,2)
    end if

    ! default value
    transposeA = 'N'
    if (present(isATransposed_opt)) then
      if (isATransposed_opt) then
        transposeA = 'T'
      end if
    end if

    ! default value
    transposeB = 'N'
    if (present(isBTransposed_opt)) then
      if (isBTransposed_opt) then
        transposeB = 'T'
      end if
    end if

    if (present(summationDim_opt)) then
      summationDim = summationDim_opt
    else
      if (transposeA == 'N') then
         summationDim = size(AmatrixIn,2)
      else
         summationDim = size(AmatrixIn,1)
      end if
    end if

    if (present(isCSymmetric_opt)) then
      isCSymmetric = isCSymmetric_opt
    else
      isCSymmetric = .false.
    end if

    call utl_tmg_start(184,'low-level--linalg_fastMatMul')

    if (isCSymmetric) then
      ! https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2025-0/gemmt.html
      call dgemmt('U', transposeA, transposeB, &
                   lastDim, summationDim,      & ! N, K
                   1.0d0,                      & ! alpha
                   AmatrixIn,  firstDimA,      & ! A
                   BmatrixIn,  firstDimB,      & ! B
                   0.0d0,                      & ! beta
                   CmatrixOut, firstDimC)        ! C

      ! Copy upper triangle to lower triangle (symmetric matrix)
      !$OMP PARALLEL DO PRIVATE (iIndex,jIndex)
      do jIndex = 1, lastDim
        do iIndex = jIndex+1, lastDim
          CmatrixOut(iIndex,jIndex) = CmatrixOut(jIndex,iIndex)
        end do
      end do
      !$OMP END PARALLEL DO
    else
      ! https://www.netlib.org/lapack/explore-html/dd/d09/group__gemm_ga1e899f8453bcbfde78e91a86a2dab984.html#ga1e899f8453bcbfde78e91a86a2dab984
      call dgemm(transposeA, transposeB, &
                 firstDim, lastDim, summationDim, & ! M, N, K
                 1.0d0,                  & ! alpha
                 AmatrixIn,  firstDimA,  & ! A
                 BmatrixIn,  firstDimB,  & ! B
                 0.0d0,                  & ! beta
                 CmatrixOut, firstDimC)    ! C
    end if

    call utl_tmg_stop(184)

  end subroutine linalg_fastMatMul

end module linearAlgebra_mod
