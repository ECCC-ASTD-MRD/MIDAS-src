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

module mpivar_mod
  ! MODULE mpivar_mod (prefix='mpivar' category='8. Low-level utilities and constants')
  !
  ! :Purpose: Subroutine and public variables related to the mpi decomposition
  !           specific to the MIDAS code. Depends on the more general mpi_mod
  !           module.
  !
  use mpi_mod
  use utilities_mod
  implicit none
  save
  private

  ! public procedures
  public :: mpivar_setup_latbands, mpivar_setup_lonbands
  public :: mpivar_setup_m, mpivar_setup_n
  public :: mpivar_setup_levels
  public :: mpivar_setup_varslevels
  public :: mpivar_myidXfromLon, mpivar_myidYfromLat

  contains


  subroutine mpivar_setup_latbands(nj, latPerPE, latPerPEmax, myLatBeg, myLatEnd,  &
                                   myLatHalfBeg_opt, myLatHalfEnd_opt, divisible_opt)
    ! :Purpose: compute parameters that define the mpi distribution of
    !          latitudes over tasks in Y direction (npey)
    implicit none
    integer           :: nj
    integer           :: latPerPE
    integer           :: latPerPEmax
    integer           :: myLatBeg
    integer           :: myLatEnd
    integer           :: njlath
    integer, optional :: myLatHalfBeg_opt
    integer, optional :: myLatHalfEnd_opt
    logical, optional :: divisible_opt

    integer :: latPerPEmin, ierr
    logical, save :: firstCall = .true.

    latPerPEmin = floor(real(nj) / real(mpi_npey))
    myLatBeg = 1 + (mpi_myidy * latPerPEmin)
    if( mpi_myidy < (mpi_npey-1) ) then
      myLatEnd = (1 + mpi_myidy) * latPerPEmin
    else
      myLatEnd = nj
    end if
    latPerPE = myLatEnd - myLatBeg + 1
    call rpn_comm_allreduce(latPerPE,latPerPEmax,1,'MPI_INTEGER','MPI_MAX','NS',ierr)

    if( firstCall ) then
      write(*,'(a,4i8)') 'mpivar_setup_latbands: latPerPE, latPerPEmax, myLatBeg, myLatEnd = ',  &
           latPerPE, latPerPEmax, myLatBeg, myLatEnd
      firstCall = .false.
    end if

    if (present(myLatHalfBeg_opt).and.present(myLatHalfEnd_opt)) then
      njlath = (nj + 1) / 2
      if (myLatBeg <= njlath .and. myLatEnd <= njlath) then
        myLatHalfBeg_opt = myLatBeg
        myLatHalfEnd_opt = myLatEnd
      elseif (myLatBeg >= njlath .and. myLatEnd >= njlath) then
        myLatHalfBeg_opt = 1 + nj - myLatEnd
        myLatHalfEnd_opt = 1 + nj - myLatBeg
      else
        myLatHalfBeg_opt = min(myLatBeg, 1 + nj - myLatEnd)
        myLatHalfEnd_opt = njlath
      end if
    end if

    if( present(divisible_opt) ) then
      divisible_opt = (latPerPEmin * mpi_npey == nj)
    end if

  end subroutine mpivar_setup_latbands


  function mpivar_myidYfromLat(latIndex, nj) result(IP_y)
    ! :Purpose: use same logic as setup_latbands to compute myidy
    !          corresponding to a latitude grid index
    implicit none
    integer :: latIndex
    integer :: nj
    integer :: IP_y

    IP_y = (latIndex-1) / floor( real(nj) / real(mpi_npey) )
    IP_y = min( mpi_npey-1, IP_y )

  end function mpivar_myidYfromLat


  subroutine mpivar_setup_lonbands(ni, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd, divisible_opt)
    ! :Purpose: compute parameters that define the mpi distribution of
    !          longitudes over tasks in X direction (npex)
    implicit none

    integer          :: ni
    integer          :: lonPerPE
    integer          :: lonPerPEmax
    integer          :: myLonBeg
    integer          :: myLonEnd
    logical, optional :: divisible_opt

    integer :: lonPerPEmin, ierr
    logical, save :: firstCall = .true.

    lonPerPEmin = floor(real(ni) / real(mpi_npex))
    myLonBeg = 1 + (mpi_myidx * lonPerPEmin)
    if( mpi_myidx < (mpi_npex-1) ) then
      myLonEnd = (1 + mpi_myidx) * lonPerPEmin
    else
      myLonEnd = ni
    end if
    lonPerPE = myLonEnd - myLonBeg + 1
    call rpn_comm_allreduce(lonPerPE,lonPerPEmax,1,'MPI_INTEGER','MPI_MAX','EW',ierr)

    if( firstCall ) then
      write(*,'(a,4i8)') 'mpivar_setup_lonbands: lonPerPE, lonPerPEmax, myLonBeg, myLonEnd = ', &
           lonPerPE, lonPerPEmax, myLonBeg, myLonEnd
      firstCall = .false.
    end if

    if( present(divisible_opt) ) then
      divisible_opt = (lonPerPEmin * mpi_npex == ni)
    end if

  end subroutine mpivar_setup_lonbands


  function mpivar_myidXfromLon(lonIndex, ni) result(IP_x)
    ! :Purpose: use same logic as setup_lonbands to compute myidx
    !          corresponding to a longitude grid index
    implicit none

    integer :: lonIndex
    integer :: ni
    integer :: IP_x

    IP_x = (lonIndex-1) / floor( real(ni) / real(mpi_npex) )
    IP_x = min( mpi_npex-1, IP_x )

  end function mpivar_myidXfromLon


  subroutine mpivar_setup_m(ntrunc, mymBeg, mymEnd, mymSkip, mymCount)
    ! :Purpose: compute parameters that define the mpi distribution of
    !          wavenumber m over tasks in Y direction (npey)
    implicit none
    integer :: ntrunc
    integer :: mymBeg
    integer :: mymEnd
    integer :: mymSkip
    integer :: mymCount
    integer :: jm

    mymBeg = mpi_myidy
    mymEnd = ntrunc
    mymSkip = mpi_npey
    mymCount = 0
    do jm = mymBeg, mymEnd, mymSkip
      mymCount = mymCount + 1
    end do

    write(*,'(a,4i8)') 'mpivar_setup_m: mymBeg, mymEnd, mymSkip, mymCount = ', mymBeg, mymEnd, mymSkip, mymCount

  end subroutine mpivar_setup_m

 
  subroutine mpivar_setup_n(ntrunc, mynBeg, mynEnd, mynSkip, mynCount)
    ! :Purpose: compute parameters that define the mpi distribution of
    !          wavenumber n over tasks in X direction (npex)
    implicit none

    integer :: ntrunc
    integer :: mynBeg
    integer :: mynEnd
    integer :: mynSkip
    integer :: mynCount
    integer :: jn

    mynBeg = mpi_myidx
    mynEnd = ntrunc
    mynSkip = mpi_npex
    mynCount = 0
    do jn = mynBeg, mynEnd, mynSkip
      mynCount = mynCount + 1
    end do

    write(*,'(a,4i8)') 'mpivar_setup_n: mynBeg, mynEnd, mynSkip, mynCount = ', mynBeg, mynEnd, mynSkip, mynCount

  end subroutine mpivar_setup_n


  subroutine mpivar_setup_levels(numlevels, myLevBeg, myLevEnd, myLevCount)
    ! :Purpose: compute parameters that define the mpi distribution of
    !          levels over tasks in X direction (npex)
    implicit none

    integer :: numlevels
    integer :: myLevBeg
    integer :: myLevEnd
    integer :: myLevCount
    integer :: jlev
    integer :: jproc
    integer :: factor
    integer :: myLevCounts(mpi_npex)

    ! when possible, always divide into even number of levels per MPI task
    if(mod(numlevels, 2) /= 0) then
      write(*,*) 'mpivar_setup_levels: total number of levels is not even, now=', numlevels
      write(*,*) '                     therefore, if global grid, may not be able to do '
      write(*,*) '                     transforms of vor/div <-> u/v'
      factor = 1
    else
      factor = 2
    end if

    myLevCounts(:) = 0
    do jproc = 1, mpi_npex
      do jlev = jproc, (numlevels / factor), mpi_npex
        myLevCounts(jproc) = myLevCounts(jproc) + 1
      end do
    end do
    do jproc = 1, mpi_npex
      myLevCounts(jproc) = myLevCounts(jproc) * factor
    end do

    myLevCount = myLevCounts(mpi_myidx + 1)

    if (myLevCount > 0) then
      myLevBeg = 1
      do jproc = 1, mpi_myidx
        myLevBeg = myLevBeg + myLevCounts(jproc)
      end do
      myLevEnd = myLevBeg + myLevCount - 1
    else
      myLevBeg = 1
      myLevEnd = 0
    end if

    write(*,'(a,3i8)') 'mpivar_setup_levels: myLevBeg, myLevEnd, myLevCount = ',  &
         myLevBeg, myLevEnd, myLevCount

  end subroutine mpivar_setup_levels


  subroutine mpivar_setup_varslevels(numk, mykBeg, mykEnd, mykCount)
    ! :Purpose: compute parameters that define the mpi distribution of
    !          variables/levels (i.e. 1->nk) over all tasks (nprocs)
    implicit none

    integer :: numk
    integer :: mykBeg
    integer :: mykEnd
    integer :: mykCount
    integer :: jk
    integer :: jproc
    integer :: mykCounts(mpi_nprocs)

    mykCounts(:) = 0
    do jproc = 1, mpi_nprocs
      do jk = jproc, numk, mpi_nprocs
        mykCounts(jproc) = mykCounts(jproc) + 1
      end do
    end do

    mykCount = mykCounts(mpi_myid + 1)

    mykBeg = 1
    do jproc = 1, mpi_myid
      mykBeg = mykBeg + mykCounts(jproc)
    end do
    mykEnd = mykBeg + mykCount - 1

  end subroutine mpivar_setup_varslevels

end module mpivar_mod
