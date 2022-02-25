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
  use obsSpaceData_mod
  use tovs_nl_mod
  implicit none
  private
  save

  ! public variables
  public :: rmat_lnondiagr
  ! public subroutines
  public :: rmat_init,rmat_cleanup,rmat_readCMatrix,rmat_RsqrtInverseOneObs, rmat_RsqrtInverseAllObs

 type rmat_matrix
    real(8) ,pointer     :: Rmat(:,:)=>null()
    integer ,pointer     :: listChans(:)=>null()
    integer              :: nchans=0
 END type rmat_matrix

  type(rmat_matrix),target,allocatable  :: Rcorr_inst(:) ! non diagonal Correlation matrices for each instrument
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
      allocate(Rcorr_inst(nsensors))
      allocate(R_tovs(nobtovs))
    end if

  end subroutine rmat_init

  subroutine rmat_cleanup()
    implicit none
    if (rmat_lnondiagr) then
      deallocate(Rcorr_inst)
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
       call rmat_readCMatrixByFileName(filename,Rcorr_inst(sensor_id), ichan )
    else
      write(*,*) "Unknown instrument ",instrument(:)
      call utl_abort("rmat_read_C_matrix")
    end if
    
  end subroutine rmat_readCMatrix


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


  subroutine rmat_RsqrtInverseOneObs(sensor_id,nsubset,obsIn,obsOut,list_sub,list_oer,indexTovs)
    !
    ! :Purpose: Apply the operator R**-1/2 to obsIn
    !           result in obsOut for the subset of channels specified
    !           in list_sub
    !
    implicit none

    ! Arguments:
    integer , intent (in) :: sensor_id
    integer , intent (in) :: nsubset
    integer , intent(in)  :: list_sub(nsubset)
    real(8) , intent(in)  :: list_oer(nsubset)
    real(8) , intent(in)  :: obsIn(nsubset)
    real(8) , intent(out) :: obsOut(nsubset)
    integer , intent(in)  :: indexTovs

    ! locals
    real (8) :: Rsub(nsubset,nsubset), alpha, beta, product 
    integer :: index(nsubset)
    integer :: i,j

    if (R_tovs(indexTovs)%nchans == 0) then

      if (sensor_id <= 0 .or. sensor_id > size(Rcorr_inst)) then
        write(*,*) "invalid sensor_id",sensor_id,size(Rcorr_inst)
        call utl_abort('rmat_RsqrtInverseOneObs')
      end if

      index = -1
      do i=1,nsubset
        bj: do j = 1, Rcorr_inst(sensor_id)%nchans
          if (list_sub(i) == Rcorr_inst(sensor_id)%listChans(j)) then
            index(i) = j
            exit bj 
          end if
        end do bj
      end do
      if (any(index == -1)) then
        write(*,*) "Missing information for some channel !"
        write(*,*) list_sub(:)
        write(*,*) index(:)
        call utl_abort('rmat_RsqrtInverseOneObs')
      end if
      R_tovs(indexTovs)%nchans = nsubset
      allocate(R_tovs(indexTovs)%listChans(nsubset))
      R_tovs(indexTovs)%listChans(1:nsubset) = list_sub(1:nsubset)
      do j=1,nsubset
        do i=1,nsubset
          product = list_oer(i) * list_oer(j)
          Rsub(i,j) = product * Rcorr_inst(sensor_id)%Rmat(index(i),index(j))
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
    obsOut = 0.d0
    ! Optimized symetric matrix vector product from Lapack
    call dsymv("L", nsubset, alpha, R_tovs(indexTovs)%Rmat, nsubset,obsIn, 1, beta, obsOut, 1)
    call tmg_stop(96)

  end subroutine rmat_RsqrtInverseOneObs


  !--------------------------------------------------------------------------
  ! rmat_RsqrtInverseAllObs
  !--------------------------------------------------------------------------
  subroutine rmat_RsqrtInverseAllObs( obsSpaceData, elem_dest_i, elem_src_i )
    !
    !:Purpose: To apply observation-error variances to ROBDATA8(k_src,*) and to
    !          store it in the elem_src_s of obsspacedata
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsspacedata
    integer, intent(in)  :: elem_dest_i ! destination index
    integer, intent(in)  :: elem_src_i  ! source index

    ! Locals:
    integer :: bodyIndex, headerIndex
    integer :: idata, idatend, idatyp, count, channelNumber, channelIndex
    real(8) :: obsIn( tvs_maxChannelNumber ), obsOut( tvs_maxChannelNumber )
    integer :: list_chan( tvs_maxChannelNumber )
    real(8) :: list_OER( tvs_maxChannelNumber )

    ! NOTE I tried using openMP on this loop, but it increased the cost from 4sec to 80sec!!!
    do headerIndex =1, obs_numHeader(obsspacedata)

      idata   = obs_headElem_i( obsspacedata, OBS_RLN, headerIndex )
      idatend = obs_headElem_i( obsspacedata, OBS_NLV, headerIndex ) + idata - 1
      idatyp  = obs_headElem_i( obsspacedata, OBS_ITY, headerIndex )

      if ( tvs_isIdBurpTovs(idatyp) .and. rmat_lnondiagr ) then

        count = 0

        do bodyIndex = idata, idatend

          if (obs_bodyElem_i( obsspacedata, OBS_ASS, bodyIndex ) == obs_assimilated ) then
            call tvs_getChannelNumIndexFromPPP( obsSpaceData, headerIndex, bodyIndex, &
                                                channelNumber, channelIndex )
            count = count + 1
            list_chan( count ) = channelNumber
            list_OER( count ) = obs_bodyElem_r( obsspacedata, OBS_OER, bodyIndex )
            obsIn( count ) = obs_bodyElem_r( obsspacedata, elem_src_i, bodyIndex )
          end if
        
        end do

        if ( count > 0 .and. tvs_tovsIndex( headerIndex ) > 0 ) then

          call rmat_RsqrtInverseOneObs( tvs_lsensor( tvs_tovsIndex( headerIndex )), count, obsIn(1:count), obsOut(1:count), list_chan(1:count), list_OER(1:count), tvs_tovsIndex(headerIndex) )

          count = 0
          do bodyIndex = idata, idatend
            if ( obs_bodyElem_i( obsspacedata, OBS_ASS, bodyIndex ) == obs_assimilated) then
              count = count + 1
              call obs_bodySet_r(obsspacedata, elem_dest_i, bodyIndex,obsOut(count))
            end if
          end do

        else

          do bodyIndex = idata, idatend
            call obs_bodySet_r(obsspacedata, elem_dest_i, bodyIndex, 0.d0)
          end do

        end if 

      else
        if ( obs_getFamily(obsSpaceData, headerIndex_in=headerIndex) == 'TO' .and. rmat_lnondiagr) then
          call utl_abort('rmat_RsqrtInverseAllObs: should not happen check NAMTOVSINST for missing sensor')
        end if
        do bodyIndex = idata, idatend
          if (obs_bodyElem_i( obsspacedata, OBS_ASS, bodyIndex ) == obs_assimilated) then
            call obs_bodySet_r( obsspacedata, elem_dest_i, bodyIndex, &
                 obs_bodyElem_r( obsspacedata, elem_src_i, bodyIndex) / obs_bodyElem_r( obsspacedata, OBS_OER, bodyIndex ))
          end if
        end do

      end if ! is it a radiance in non diagonal R mode ?

    end do !loop on header
 
  end subroutine rmat_RsqrtInverseAllObs

end module rMatrix_mod
