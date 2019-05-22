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

module rMatrix_mod
  ! MODULE rMatrix_mod (prefix='rmat' category='5. B and R matrices')
  !
  ! :Purpose: Module to handle non-diagonal observation-error covariance
  !           matrices for assimilation of radiances
  !
  use rttov_interfaces_mod
  use mpi_mod
  use mpivar_mod
  use rttov_const, only  : errorstatus_success
  use utilities_mod
  implicit none
  private
  save

  ! public variables
  public :: rmat_lnondiagr
  ! public subroutines
  public :: rmat_init,rmat_cleanup,rmat_readCMatrix,rmat_setFullRMatrix,rmat_sqrtRm1

 type rmat_matrix
    real(8) ,pointer     :: Rmat(:,:)=>null()
    integer ,pointer     :: listChans(:)=>null()
    integer              :: nchans=0
 END type rmat_matrix

  type(rmat_matrix),target,allocatable  :: R_inst(:) ! non diagonal Covariance matrices (R) for each instrument
  type(rmat_matrix),target,allocatable  :: C_inst(:) ! non diagonal Correlation matrices for each instrument
  type(rmat_matrix),target,allocatable  :: R_tovs(:) ! non diagonal R matrices used for the assimilation of all radiances
 
  logical :: rmat_lnondiagr

  contains

  subroutine rmat_init(nsensors,nobtovs)
   
    implicit none

    ! Arguments:
    integer, intent (in) :: nsensors
    integer, intent (in) :: nobtovs

    ! locals
    integer :: nulnam,ierr
    integer ,external:: fnom,fclos
    namelist /NAMRMAT/rmat_lnonDiagR

    ! Default value for parameter rmat_lnondiagr, don't use interchannel correlation by default
    rmat_lnonDiagR = .false.

    ! Read the parameters from NAMRMAT
    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=namrmat, iostat=ierr)
    if (ierr /= 0) call utl_abort('rmat_init: Error reading namelist')
    if (mpi_myid == 0) write(*,nml=namrmat)
    ierr = fclos(nulnam)
    if (rmat_lnonDiagR) then
      allocate(R_inst(nsensors))
      allocate(C_inst(nsensors))
      allocate(R_tovs(nobtovs))
    end if

  end subroutine rmat_init

  subroutine rmat_cleanup()
    implicit none
    if (rmat_lnondiagr) then
      deallocate(R_inst)
      deallocate(C_inst)
      deallocate(R_tovs)
    end if
  end subroutine rmat_cleanup

  subroutine rmat_readCMatrix(instrument, sensor_id, ichan )
   
    implicit none
    
    ! Arguments
    integer ,intent (in) :: instrument(3)
    integer ,intent (in) :: sensor_id
    integer ,intent (in) :: ichan(:)

    ! locals
    character (len=64) :: filename
    integer :: err

    call rttov_coeffname (err, instrument, coeffname=filename, filetype="Cmat")
    
    if (err == errorstatus_success) then
       filename = trim(filename) // ".dat"
       call rmat_readCMatrixByFileName(filename,C_inst(sensor_id), ichan )
    else
      write(*,*) "Unknown instrument ",instrument(:)
      call utl_abort("rmat_read_C_matrix")
    end if
    
  end subroutine rmat_readCMatrix

  subroutine rmat_setFullRMatrix(sigma,sensor_id,offset)
    implicit none
   
    ! Arguments:
    integer, intent (in) :: sensor_id
    integer, intent (in) :: offset
    real(8), intent (in) :: sigma(:)

    ! locals
    integer :: i,nchn,j,ii,jj,nsigma
    real (8) :: product

    write(*,*) "rmat_setFullRMatrix: "

    write(*,*) "sensor_id:",sensor_id

    nsigma = size( sigma )
    if (nsigma < 1) then
      write(*,*) "rmat_setFullRMatrix: Strange sigma array size",nsigma
      write(*,*) "Please check !"
      return
    end if

    R_inst(sensor_id) % nchans  =  C_inst(sensor_id) % nchans
    if ( R_inst(sensor_id) % nchans == 0) return  
    allocate( R_inst(sensor_id)%Rmat(R_inst(sensor_id) % nchans,R_inst(sensor_id) % nchans) )
    allocate( R_inst(sensor_id)%listChans(R_inst(sensor_id) % nchans) )
    
    R_inst(sensor_id)%listChans(:) = C_inst(sensor_id)%listChans(:)

    do i=1,C_inst(sensor_id)%nchans
      ii = R_inst(sensor_id)%listChans(i) + offset
      do j=1,C_inst(sensor_id)%nchans
        jj = R_inst(sensor_id)%listChans(j) + offset
        product = sigma(ii) * sigma(jj)
        if (product <= 0.) then
          write(*,*) "Invalid variance: missing channel in stat_tovs !"
          write(*,*) ii,jj,offset,sigma(ii),sigma(jj)
          call utl_abort('rmat_setFullRMatrix')
        end if
        R_inst(sensor_id)%Rmat(i,j) = product * C_inst(sensor_id)%Rmat(i,j)
      end do
    end do

  end subroutine rmat_setFullRMatrix


  subroutine rmat_readCMatrixByFileName(infile,C,chanList_opt)
    implicit none

    ! Arguments:
    character (len=*),intent(in) :: infile ! name of input file
    type(rmat_matrix),intent(inout) :: C    ! correlation matrix structure
    integer ,intent(in),optional :: chanList_opt(:) ! list of requested channels (if missing will read all file content)

    ! locals
    integer :: i,j,iu,ierr,count,ich,nchn,nch
    integer ,external :: fnom,fclos
    real(8) :: x
    integer ,allocatable :: index(:)

    nchn = -1
    if (present(chanList_opt)) then
      nchn=size(chanList_opt)
    end if

    iu = 0
    ierr = fnom(iu,trim(infile),'FTN+SEQ+R/O',0)
    if (ierr /= 0) then
      write(*,*) "Cannot open "//trim(infile)
      call utl_abort("rmat_readCMatrixByFileName")
    end if

    write(*,*) "rmat_readCMatrixByFileName: Reading "//trim(infile)
    
    read(iu,*) nch
    if (nchn == -1) then
      nchn = nch
    else
      if(nchn > nch) then
        write(*,*) "Not enough channels in "//trim(infile)
        write(*,*) nchn,nch
        call utl_abort("rmat_readCMatrixByFileName")
      end if
    end if
    allocate(index(nch))
    
    C%nchans = nchn
    allocate(C%Rmat(nchn,nchn))
    allocate(C%listChans(nchn))
    C%Rmat = 0.d0
    do i = 1,nchn
      C%Rmat(i,i) = 1.d0
    end do
    count = 0
    index = -1
    do i=1,nch
      read(iu,*) ich
      if (present(chanList_opt)) then
        bj:do j=1,nchn
          if (ich == chanList_opt(j)) then
            count = count + 1
            index(i) = j
            C%listChans(count) = ich
            exit bj
          end if
        end do bj
      else
        index(i) = i
        C%listChans(i) = ich
        count = count + 1
      end if
    end do
    if (count /= nchn) then
      write(*,*) "Warning: Missing information in "//trim(infile)
      do j=1,nchn
        write(*,*) j, chanList_opt(j) 
      end do
      write(*,*) "Not important if there is no observation of this family"
    end if

    do
      read(iu,*,iostat=ierr) i,j,x
      if (ierr /= 0) exit
      if (index(i) /= -1 .and. index(j) /= -1) then
        C%Rmat(index(i),index(j)) = x
        C%Rmat(index(j),index(i)) = x
      end if
    end do

    ierr= fclos(iu)
    deallocate(index)

  end subroutine rmat_readCMatrixByFileName


  subroutine rmat_sqrtRm1(sensor_id,nsubset,x,y,list_sub,indexTovs)
    !
    ! :Purpose: Apply the operator R**-1/2 to x
    !           result in y for the subset of channels specified
    !           in list_sub
    !
    implicit none

    ! Arguments:
    integer , intent (in) :: sensor_id
    integer , intent (in) :: nsubset
    integer , intent(in)  :: list_sub(nsubset)
    real(8) , intent(in)  :: x(nsubset)
    real(8) , intent(out) :: y(nsubset)
    integer , intent(in)  :: indexTovs

    ! locals
    real (8) :: Rsub(nsubset,nsubset),alpha,beta
    integer :: index(nsubset)
    integer :: i,j
    type(rmat_matrix),pointer :: R

    if (R_tovs(indexTovs)%nchans == 0) then

      if (sensor_id > 0 .and. sensor_id <= size(R_inst)) then
        R => R_inst(sensor_id)
      else 
        write(*,*) "invalid sensor_id",sensor_id,size(R_inst)
        call utl_abort('rmat_sqrtRm1')
      end if

      index = -1
      do i=1,nsubset
        bj: do j=1,R%nchans
          if (list_sub(i) == R%listChans(j)) then
            index(i) = j
            exit bj 
          end if
        end do bj
      end do
      if (any(index == -1)) then
        write(*,*) "Missing information for some channel !"
        write(*,*) list_sub(:)
        write(*,*) index(:)
        call utl_abort('rmat_sqrtRm1')
      end if
      R_tovs(indexTovs)%nchans = nsubset
      allocate(R_tovs(indexTovs)%listChans(nsubset))
      R_tovs(indexTovs)%listChans(1:nsubset) = list_sub(1:nsubset)
      do j=1,nsubset
        do i=1,nsubset
          Rsub(i,j) = R%Rmat(index(i),index(j))
        end do
      end do
      ! Calculation of R**-1/2
      call tmg_start(95,'RMAT_MATSQRT')
      call utl_matSqrt(Rsub,nsubset,-1.d0,.false.)
      call tmg_stop(95)
      allocate(R_tovs(indexTovs)%Rmat(nsubset,nsubset))
      do j=1,nsubset
        do i=1,nsubset
          R_tovs(indexTovs)%Rmat(i,j) = Rsub(i,j)
        end do
      end do
    end if

    call tmg_start(96,'RMAT_MATMUL')
    alpha = 1.d0
    beta = 0.d0
    y = 0.d0
    ! Optimized symetric matrix vector product from Lapack
    call dsymv("L", nsubset, alpha, R_tovs(indexTovs)%Rmat, nsubset,x, 1, beta, y, 1)
    call tmg_stop(96)

  end subroutine rmat_sqrtRm1

end module rMatrix_mod
