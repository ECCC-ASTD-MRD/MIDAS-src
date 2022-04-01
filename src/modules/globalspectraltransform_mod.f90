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

module globalSpectralTransform_mod
  ! MODULE globalSpectralTransform_mod (prefix='gst' category='3. High-level transformations')
  !
  ! :Purpose: To perform global spectral transform (spherical harmonic transform
  !           with grid-point field on a standard global Gaussian grid). 
  !
  use codePrecision_mod
  use mpi
  use mpi_mod
  use mpivar_mod
  use MathPhysConstants_mod
  use earthConstants_mod
  use utilities_mod
  implicit none
  save
  private

  ! public subroutines
  public :: gst_setup,  &
            gst_speree, gst_speree_ad, gst_speree_kij, gst_speree_kij_ad, gst_reespe, gst_reespe_kij, &
            gst_spgd, gst_spgda, gst_gdsp, gst_zleginv, gst_zlegdir,  &
            gst_setID, gst_setDefaultID, gst_setToDefaultID,  &
            gst_ilaList_mpilocal, gst_ilaList_mpiglobal
  ! public functions
  public :: gst_getRmu, gst_getRwt, gst_getnind, gst_getrlati, gst_getr1qm2, gst_getrsqm2, &
            gst_getrnnp1, gst_getr1snp1, gst_getzleg, gst_getNla


  type  :: T_gst
    real(8),allocatable   :: rmu(:)
    real(8),allocatable   :: rwt(:)
    real(8),allocatable   :: rsqm2(:)
    real(8),allocatable   :: r1qm2(:)
    real(8),allocatable   :: rlati(:)
    integer,allocatable   :: nind(:)
    integer,allocatable   :: nindrh(:)
    real(8),allocatable   :: dalp(:,:)
    real(8),allocatable   :: dealp(:,:)
    real(8),allocatable   :: rwocs(:)
    real(8),allocatable   :: r1mu2(:)
    real(8),allocatable   :: rcolat(:)
    real(8),allocatable   :: r1mui(:)
    real(8),allocatable   :: r1mua(:)
    integer,allocatable   :: nclm(:)
    real(8),allocatable   :: r1snp1(:)
    real(8),allocatable   :: rnnp1(:)
    real(8),allocatable   :: zleg(:,:)
    integer               :: ntrunc
    integer               :: nla
    integer               :: nlarh
    integer               :: njlath
    integer               :: ni
    integer               :: nj
    integer               :: nk
    integer               :: myLatBeg, myLatEnd, latPerPE, latPerPEmax
    integer,allocatable   :: allLatBeg(:), allLatEnd(:), allLatPerPE(:)
    integer               :: myLonBeg, myLonEnd, lonPerPE, lonPerPEmax
    integer,allocatable   :: allLonBeg(:), allLonEnd(:), allLonPerPE(:)
    integer               :: myLatHalfBeg, myLatHalfEnd
    integer               :: mymBeg, mymEnd, mymSkip, mymCount, maxmCount
    integer               :: mynBeg, mynEnd, mynSkip, mynCount
    integer               :: myNla, maxMyNla
    integer,allocatable   :: allNla(:)
    integer,allocatable   :: mymIndex(:)
    integer,pointer       :: ilaList(:), allIlaList(:,:)
    integer,allocatable   :: allmBeg(:), allmEnd(:), allmSkip(:)
    integer,allocatable   :: allnBeg(:), allnEnd(:), allnSkip(:)
    integer               :: myLevBeg, myLevEnd, myLevCount, maxMyLevCount
    integer,allocatable   :: allLevBeg(:), allLevEnd(:)
    integer               :: sendType_LevToLon, recvType_LevToLon
    integer               :: sendType_LonToLev, recvType_LonToLev
    logical               :: lonLatDivisible
  end type T_gst

  integer,parameter :: nMaxGst = 20
  integer      :: gstIdDefault = -1
  integer      :: nGstAlreadyAllocated = 0
  integer      :: gstID = 0
  type(T_gst)  :: gst(nmaxgst)

  integer :: nlatbd = 8

contains

  subroutine GST_setID(gstID_in)
    implicit none

    ! Arguments:
    integer :: gstID_in
  
    gstID = gstID_in

  end subroutine GST_setID


  subroutine GST_setDefaultID(gstID_in)
    implicit none

    ! Arguments:
    integer :: gstID_in
  
    gstIDDefault = gstID_in

  end subroutine GST_setDefaultID


  subroutine GST_setToDefaultID
    implicit none

    gstID = gstIdDefault

  end subroutine GST_setToDefaultID


  integer function GST_getNla(gstID_opt)
    implicit none

    ! Arguments:
    integer,optional :: gstID_opt

    ! Locals:
    integer :: gstID_l

    if(present(gstID_opt)) then
      gstID_l = gstID_opt
    else
      gstID_l = gstIdDefault
    endif

    gst_getNla = gst(gstID_l)%myNla

  end function GST_getNla


  real(8) function GST_getRmu(latIndex,gstID_opt)
    implicit none

    ! Arguments:
    integer :: latIndex
    integer,optional :: gstID_opt

    ! Locals:
    integer :: gstID_l
    integer :: latIndex2

    if(present(gstID_opt)) then
      gstID_l = gstID_opt
    else
      gstID_l = gstIdDefault
    endif

    ! use a flipped latitude (south to north)
    latIndex2 = gst(gstID_l)%nj - latIndex + 1
    gst_getRmu = gst(gstID_l)%rmu(latIndex2)

  end function GST_getRmu


  real(8) function GST_getRnnp1(ilaIndex,gstID_opt)
    implicit none

    ! Arguments:
    integer :: ilaIndex
    integer,optional :: gstID_opt

    ! Locals:
    integer :: gstID_l

    if(present(gstID_opt)) then
      gstID_l = gstID_opt
    else
      gstID_l = gstIdDefault
    endif

    gst_getRnnp1 = gst(gstID_l)%rnnp1(ilaIndex)

  end function GST_getRnnp1


  real(8) function GST_getR1snp1(ilaIndex,gstID_opt)
    implicit none

    ! Arguments:
    integer :: ilaIndex
    integer,optional :: gstID_opt

    ! Locals:
    integer :: gstID_l

    if(present(gstID_opt)) then
      gstID_l = gstID_opt
    else
      gstID_l = gstIdDefault
    endif

    gst_getR1snp1 = gst(gstID_l)%r1snp1(ilaIndex)

  end function GST_getR1snp1


  real(8) function GST_getRwt(latIndex,gstID_opt)
    implicit none

    ! Arguments:
    integer :: latIndex
    integer,optional :: gstID_opt

    ! Locals:
    integer :: gstID_l
    integer :: latIndex2
  
    if(present(gstID_opt)) then
      gstID_l = gstID_opt
    else
      gstID_l = gstIdDefault
    endif

    ! use a flipped latitude (south to north)
    latIndex2 = gst(gstID_l)%nj - latIndex + 1
    gst_getRwt = gst(gstID_l)%rwt(latIndex2)

  end function GST_getRwt


  integer function GST_getNind(mIndex,gstID_opt)
    implicit none

    ! Arguments:
    integer :: mIndex
    integer,optional :: gstID_opt

    ! Locals:
    integer :: gstID_l
  
    if(present(gstID_opt)) then
      gstID_l = gstID_opt
    else
      gstID_l = gstIdDefault
    endif

    gst_getNind = gst(gstID_l)%nind(mIndex)

  end function GST_getNind


  real(8) function GST_getRlati(latIndex,gstID_opt)
    implicit none

    ! Arguments:
    integer :: latIndex
    integer,optional :: gstID_opt
    integer :: latIndex2

    ! Locals:
    integer :: gstID_l
 
    if(present(gstID_opt)) then
      gstID_l = gstID_opt
    else
      gstID_l = gstIdDefault
    endif

    ! use a flipped latitude (south to north)
    latIndex2 = gst(gstID_l)%nj - latIndex + 1
    gst_getRlati = gst(gstID_l)%rlati(latIndex2)

  end function GST_getRlati


  real(8) function GST_getR1qm2(latIndex,gstID_opt)
    implicit none

    ! Arguments:
    integer :: latIndex
    integer,optional :: gstID_opt

    ! Locals:
    integer :: gstID_l
    integer :: latIndex2
  
    if(present(gstID_opt)) then
      gstID_l = gstID_opt
    else
      gstID_l = gstIdDefault
    endif

    ! use a flipped latitude (south to north)
    latIndex2 = gst(gstID_l)%nj - latIndex + 1
    gst_getR1qm2 = gst(gstID_l)%r1qm2(latIndex2)

  end function GST_getR1qm2


  real(8) function GST_getRsqm2(latIndex,gstID_opt)
    implicit none

    ! Arguments:
    integer :: latIndex
    integer,optional :: gstID_opt

    ! Locals:
    integer :: gstID_l
    integer :: latIndex2

    if(present(gstID_opt)) then
      gstID_l = gstID_opt
    else
      gstID_l = gstIdDefault
    endif

    ! use a flipped latitude (south to north)
    latIndex2 = gst(gstID_l)%nj - latIndex + 1
    gst_getRsqm2 = gst(gstID_l)%rsqm2(latIndex2)

  end function GST_getRsqm2


  real(8) function GST_getzleg(legendreIndex,latIndex,gstID_in)
    ! 
    !:Purpose: To pass on Legendre polynomial element 
    !
    implicit none

    ! Arguments:
    integer :: legendreIndex,latIndex,gstID_in

    gst_getzleg=gst(gstID_in)%zleg(legendreIndex,latIndex)

  end function GST_getzleg


  subroutine GST_ilaList_mpiglobal(ilaList,myNla,maxMyNla,gstID_in,mymBeg,&
                                   mymEnd,mymSkip,mynBeg,mynEnd,mynSkip)
    !
    !:Purpose: To produce an array to convert an mpilocal "ila" into an
    !          mpiglobal "ila"
    implicit none

    ! Arguments:
    integer, pointer :: ilaList(:)
    integer          :: myNla, maxMyNla
    integer          :: gstID_in, mymBeg, mymEnd, mymSkip
    integer          :: mynBeg, mynEnd, mynSkip

    ! Locals:
    integer          :: jm, jn, ierr

    ! compute mpilocal value of nla
    myNla = 0
    do jm = mymBeg, mymEnd, mymSkip
      do jn = mynBeg, mynEnd, mynSkip
        if(jm.le.jn) then
          myNla = myNla + 1
        endif
      enddo
    enddo

    ! determine maximum value of myNla over all processors (used for dimensioning)
    call rpn_comm_allreduce(myNla,maxMyNla,1,'MPI_INTEGER','MPI_MAX','GRID',ierr)

    allocate(ilaList(maxMyNla))
    ilaList(:) = 0

    myNla = 0
    do jm = mymBeg, mymEnd, mymSkip
      do jn = mynBeg, mynEnd, mynSkip
        if(jm.le.jn) then
          myNla = myNla + 1
          ilaList(myNla) = gst_getNind(jm,gstID_in) + jn - jm
        endif
      enddo
    enddo

  end subroutine GST_ilaList_mpiglobal


  subroutine GST_ilaList_mpilocal(ilaList,gstID_in,mymBeg,mymEnd,mymSkip,&
                                  mynBeg,mynEnd,mynSkip)
    !
    !:Purpose: To produce an array to convert an mpiglobal "ila" into an
    !          mpilocal "ila"
    implicit none

    ! Arguments:
    integer, pointer :: ilaList(:)
    integer          :: gstID_in, mymBeg, mymEnd, mymSkip
    integer          :: mynBeg, mynEnd, mynSkip

    ! Locals:
    integer          :: jm, jn, myNla

    ! assume mpiglobal value of nla already set in gst structure
    allocate(ilaList(gst(gstID_in)%nla))
    ilaList(:) = 0

    myNla = 0
    do jm = mymBeg, mymEnd, mymSkip
      do jn = mynBeg, mynEnd, mynSkip
        if(jm.le.jn) then
          myNla = myNla + 1
          ilaList(gst_getNind(jm,gstID_in) + jn - jm) = myNla
        endif
      enddo
    enddo

  end subroutine GST_ilaList_mpilocal


  integer function gst_setup(ni_in,nj_in,ntrunc_in,maxlevels_in)
    implicit none

    ! Arguments:
    integer  :: ni_in, nj_in, ntrunc_in
    integer  :: maxlevels_in

    ! Locals:
    integer  :: jn, jm, ila, ierr
    integer  :: latPerPE, latPerPEmax, myLatBeg, myLatEnd, myLatHalfBeg, myLatHalfEnd
    integer  :: lonPerPE, lonPerPEmax, myLonBeg, myLonEnd
    integer  :: myLevBeg, myLevEnd, myLevCount
    integer  :: mymBeg, mymEnd, mymSkip, mymCount
    integer  :: mynBeg, mynEnd, mynSkip, mynCount
    real(8)  :: znnp1, z1snp1
    integer(kind=MPI_ADDRESS_KIND) :: lowerBound, extent
    integer :: realSize, sendType, recvType
    logical :: divisibleLon, divisibleLat

    if(nGstAlreadyAllocated.eq.nMaxGst) then
      if(mpi_myid.eq.0) write(*,*) 'gst_setup: The maxmimum number of spectral transform have already been allocated! ', nMaxGst
      call utl_abort('gst_setup')
    endif

    nGstAlreadyAllocated = nGstAlreadyAllocated+1
    gstID = nGstAlreadyAllocated
    call gst_setDefaultID(gstID)
    if(mpi_myid.eq.0) write(*,*) 'gst_setup: Now setting up spectral transform #', gstID
    if(mpi_myid.eq.0) write(*,*) 'gst_setup: following *kind* used for internal reals : ', pre_specTransReal

    gst(gstID)%ni = ni_in
    gst(gstID)%nj = nj_in
    gst(gstID)%njlath = (gst(gstID)%nj + 1)/2

    if (ntrunc_in == -1) then ! no truncation case
      gst(gstID)%ntrunc = gst(gstID)%ni / 2
    else
      gst(gstID)%ntrunc = ntrunc_in
    end if
    gst(gstID)%nla = (gst(gstID)%ntrunc + 1) * (gst(gstID)%ntrunc +2)/2
    gst(gstID)%nlarh = (gst(gstID)%ntrunc+1) * (gst(gstID)%ntrunc+1)
    if(mpi_myid.eq.0) write(*,*) 'gst_setup: ntrunc=', gst(gstID)%ntrunc

    call mpivar_setup_latbands(gst(gstID)%nj,  &
         latPerPE, latPerPEmax, myLatBeg, myLatEnd, myLatHalfBeg, myLatHalfEnd, divisible_opt=divisibleLat)
    call mpivar_setup_lonbands(gst(gstID)%ni,  &
         lonPerPE, lonPerPEmax, myLonBeg, myLonEnd, divisible_opt= divisibleLon)

    gst(gstID)%lonLatDivisible = (divisibleLon .and. divisibleLat)
    if( mpi_myid == 0 ) write(*,*) 'gst_setup: lonLatDivisible = ', gst(gstID)%lonLatDivisible

    gst(gstID)%nk = maxlevels_in
    ! 2D MPI decomposition: split levels across npex
    call mpivar_setup_levels(maxlevels_in, myLevBeg, myLevEnd, myLevCount)
    write(*,*) 'gst_setup: myLevBeg,End,Count=', myLevBeg, myLevEnd, myLevCount

    !!! Distribution of lon/lat tiles (gridpoint space) and n/m (spectral space)
    ! range of lons handled by this processor
    gst(gstID)%myLonBeg = myLonBeg
    gst(gstID)%myLonEnd = myLonEnd
    gst(gstID)%lonPerPE = lonPerPE 
    gst(gstID)%lonPerPEmax = lonPerPEmax
    ! range of lats handled by this processor
    gst(gstID)%myLatBeg = myLatBeg
    gst(gstID)%myLatEnd = myLatEnd
    gst(gstID)%myLatHalfBeg = myLatHalfBeg
    gst(gstID)%myLatHalfEnd = myLatHalfEnd
    gst(gstID)%latPerPE = latPerPE 
    gst(gstID)%latPerPEmax = latPerPEmax
    ! range of n handled by this processor
    call mpivar_setup_n(gst(gstID)%ntrunc,mynBeg,mynEnd,mynSkip,mynCount)
    gst(gstID)%mynBeg = mynBeg
    gst(gstID)%mynEnd = mynEnd
    gst(gstID)%mynSkip = mynSkip
    gst(gstID)%mynCount = mynCount
    ! range of m handled by this processor
    call mpivar_setup_m(gst(gstID)%ntrunc,mymBeg,mymEnd,mymSkip,mymCount)
    gst(gstID)%mymBeg = mymBeg
    gst(gstID)%mymEnd = mymEnd
    gst(gstID)%mymSkip = mymSkip
    gst(gstID)%mymCount = mymCount
    call rpn_comm_allreduce(gst(gstID)%mymCount,gst(gstID)%maxmCount, &
                            1,'MPI_INTEGER','MPI_MAX','GRID',ierr)
    ! range of levels handled by this processor when in spectral space
    gst(gstID)%myLevBeg = myLevBeg
    gst(gstID)%myLevEnd = myLevEnd      
    gst(gstID)%myLevCount = myLevCount
    call rpn_comm_allreduce(gst(gstID)%myLevCount,gst(gstID)%maxMyLevCount, &
                            1,'MPI_INTEGER','MPI_MAX','GRID',ierr)

    if(mpi_myid.eq.0) write(*,*) 'gst_setup: allocating comleg...'
    call allocate_comleg
    if(mpi_myid.eq.0) write(*,*) 'gst_setup: calling suleg...'
    call suleg(.false.)
    if(mpi_myid.eq.0) write(*,*) 'gst_setup: calling sualp...'
    call sualp

    ! compute the list of spectral coefficient indices when distributed by n and m
    call gst_ilaList_mpiglobal(gst(gstID)%ilaList,gst(gstID)%myNla,gst(gstID)%maxMyNla,gstID,  &
                               gst(gstID)%mymBeg,gst(gstID)%mymEnd,gst(gstID)%mymSkip,  &
                               gst(gstID)%mynBeg,gst(gstID)%mynEnd,gst(gstID)%mynSkip)
    allocate(gst(gstID)%allNla(mpi_npex))
    call rpn_comm_allgather(gst(gstID)%myNla,1,'mpi_integer',  &
                            gst(gstID)%allNla,1,'mpi_integer','EW',ierr)
    allocate(gst(gstID)%allIlaList(gst(gstID)%maxMyNla,mpi_npex))
    call rpn_comm_allgather(gst(gstID)%ilaList,gst(gstID)%maxMyNla,'mpi_integer',  &
                            gst(gstID)%allIlaList,gst(gstID)%maxMyNla,'mpi_integer','EW',ierr)

    allocate(gst(gstID)%allLonBeg(mpi_npex))
    CALL rpn_comm_allgather(gst(gstID)%myLonBeg,1,'mpi_integer',       &
                            gst(gstID)%allLonBeg,1,'mpi_integer','EW',ierr)
    allocate(gst(gstID)%allLonEnd(mpi_npex))
    CALL rpn_comm_allgather(gst(gstID)%myLonEnd,1,'mpi_integer',       &
                            gst(gstID)%allLonEnd,1,'mpi_integer','EW',ierr)
    allocate(gst(gstID)%allLonPerPE(mpi_npex))
    CALL rpn_comm_allgather(gst(gstID)%lonPerPE,1,'mpi_integer',       &
                            gst(gstID)%allLonPerPE,1,'mpi_integer','EW',ierr)

    allocate(gst(gstID)%allLatBeg(mpi_npey))
    CALL rpn_comm_allgather(gst(gstID)%myLatBeg,1,'mpi_integer',       &
                            gst(gstID)%allLatBeg,1,'mpi_integer','NS',ierr)
    allocate(gst(gstID)%allLatEnd(mpi_npey))
    CALL rpn_comm_allgather(gst(gstID)%myLatEnd,1,'mpi_integer',       &
                            gst(gstID)%allLatEnd,1,'mpi_integer','NS',ierr)
    allocate(gst(gstID)%allLatPerPE(mpi_npey))
    CALL rpn_comm_allgather(gst(gstID)%latPerPE,1,'mpi_integer',       &
                            gst(gstID)%allLatPerPE,1,'mpi_integer','NS',ierr)

    allocate(gst(gstID)%mymIndex(gst(gstID)%mymBeg:gst(gstID)%mymEnd))
    gst(gstID)%mymIndex(:) = 0
    do jm = gst(gstID)%mymBeg,gst(gstID)%mymEnd,gst(gstID)%mymSkip
      if(jm.eq.gst(gstID)%mymBeg) then
        gst(gstID)%mymIndex(jm) = 1
      else
        gst(gstID)%mymIndex(jm) = gst(gstID)%mymIndex(jm-gst(gstID)%mymSkip) + 1
      endif
      write(*,*) 'gst_setup: mymIndex(',jm,')=',gst(gstID)%mymIndex(jm)
    enddo

    allocate(gst(gstID)%allnBeg(mpi_npex))
    CALL rpn_comm_allgather(gst(gstID)%mynBeg,1,'mpi_integer',       &
                            gst(gstID)%allnBeg,1,'mpi_integer','EW',ierr)
    allocate(gst(gstID)%allnEnd(mpi_npex))
    CALL rpn_comm_allgather(gst(gstID)%mynEnd,1,'mpi_integer',       &
                            gst(gstID)%allnEnd,1,'mpi_integer','EW',ierr)
    allocate(gst(gstID)%allnSkip(mpi_npex))
    CALL rpn_comm_allgather(gst(gstID)%mynSkip,1,'mpi_integer',       &
                            gst(gstID)%allnSkip,1,'mpi_integer','EW',ierr)

    allocate(gst(gstID)%allmBeg(mpi_npey))
    CALL rpn_comm_allgather(gst(gstID)%mymBeg,1,'mpi_integer',       &
                            gst(gstID)%allmBeg,1,'mpi_integer','NS',ierr)
    allocate(gst(gstID)%allmEnd(mpi_npey))
    CALL rpn_comm_allgather(gst(gstID)%mymEnd,1,'mpi_integer',       &
                            gst(gstID)%allmEnd,1,'mpi_integer','NS',ierr)
    allocate(gst(gstID)%allmSkip(mpi_npey))
    CALL rpn_comm_allgather(gst(gstID)%mymSkip,1,'mpi_integer',       &
                            gst(gstID)%allmSkip,1,'mpi_integer','NS',ierr)

    allocate(gst(gstID)%allLevBeg(mpi_npex))
    CALL rpn_comm_allgather(gst(gstID)%myLevBeg,1,'mpi_integer',       &
                            gst(gstID)%allLevBeg,1,'mpi_integer','EW',ierr)
    allocate(gst(gstID)%allLevEnd(mpi_npex))
    CALL rpn_comm_allgather(gst(gstID)%myLevEnd,1,'mpi_integer',       &
                            gst(gstID)%allLevEnd,1,'mpi_integer','EW',ierr)

    if(mpi_myid.eq.0) write(*,*) 'gst_setup: allLonBeg=',gst(gstID)%allLonBeg
    if(mpi_myid.eq.0) write(*,*) 'gst_setup: allLonEnd=',gst(gstID)%allLonEnd
    if(mpi_myid.eq.0) write(*,*) 'gst_setup: allLatBeg=',gst(gstID)%allLatBeg
    if(mpi_myid.eq.0) write(*,*) 'gst_setup: allLatEnd=',gst(gstID)%allLatEnd
    if(mpi_myid.eq.0) write(*,*) 'gst_setup: allnBeg=',gst(gstID)%allnBeg
    if(mpi_myid.eq.0) write(*,*) 'gst_setup: allnEnd=',gst(gstID)%allnEnd
    if(mpi_myid.eq.0) write(*,*) 'gst_setup: allnSkip=',gst(gstID)%allnSkip
    if(mpi_myid.eq.0) write(*,*) 'gst_setup: allmBeg=',gst(gstID)%allmBeg
    if(mpi_myid.eq.0) write(*,*) 'gst_setup: allmEnd=',gst(gstID)%allmEnd
    if(mpi_myid.eq.0) write(*,*) 'gst_setup: allmSkip=',gst(gstID)%allmSkip
    if(mpi_myid.eq.0) write(*,*) 'gst_setup: allLevBeg=',gst(gstID)%allLevBeg
    if(mpi_myid.eq.0) write(*,*) 'gst_setup: allLevEnd=',gst(gstID)%allLevEnd
    write(*,*) 'gst_setup: mymCount=',gst(gstID)%mymCount
    write(*,*) 'gst_setup: maxmCount=',gst(gstID)%maxmCount
    write(*,*) 'gst_setup: myNla=',gst(gstID)%myNla
    write(*,*) 'gst_setup: maxMyNla=',gst(gstID)%maxMyNla

    allocate(gst(gstID)%r1snp1(gst(gstID)%nla))
    allocate(gst(gstID)%rnnp1(gst(gstID)%nla))
    gst(gstID)%r1snp1(1) = 0.d0
    gst(gstID)%rnnp1(1) = 0.d0
    do jn = 1, gst(gstID)%ntrunc
       znnp1  = -real(jn,8)*real(jn+1,8)
       z1snp1 =  1.d0/znnp1
       do jm = 0, jn
          ila = gst(gstID)%nind(jm) + jn - jm
          gst(gstID)%r1snp1(ila) = z1snp1
          gst(gstID)%rnnp1(ila) = znnp1
       enddo
    enddo

    call gst_zlegpol(gstID)

    ! setup mpi derived types used in transposes (only used when grid is divisible)
    ! ... mpi_type_vector(count, blocklength, stride, ...)
    ! ... mpi_type_create_resized(oldtype, lowerbound, extent(in bytes), newtype, ierr)

    call mpi_type_size(pre_specTransMpiType, realSize, ierr)
    lowerBound = 0

    ! create the send type for LevToLon
    extent = gst(gstID)%maxMyLevCount * gst(gstID)%lonPerPE * realSize
    call mpi_type_vector(gst(gstID)%latPerPE, gst(gstID)%maxMyLevCount * gst(gstID)%lonPerPE,  &
                         gst(gstID)%maxMyLevCount * gst(gstID)%ni, pre_specTransMpiType, sendtype, ierr)
    call mpi_type_create_resized(sendtype, lowerBound , extent, gst(gstID)%sendType_LevToLon, ierr);
    call mpi_type_commit(gst(gstID)%sendType_LevToLon,ierr)

    ! create the receive type for LevToLon
    extent = gst(gstID)%maxMyLevCount * realSize
    call mpi_type_vector(gst(gstID)%lonPerPE * gst(gstID)%latPerPE , gst(gstID)%maxMyLevCount,  &
                         gst(gstID)%nk, pre_specTransMpiType, recvtype, ierr);
    call mpi_type_create_resized(recvtype, lowerBound, extent, gst(gstID)%recvType_LevToLon, ierr);
    call mpi_type_commit(gst(gstID)%recvType_LevToLon, ierr);

    ! create the send type for LonToLev
    extent = gst(gstID)%maxMyLevCount * realSize
    call mpi_type_vector(gst(gstID)%lonPerPE * gst(gstID)%latPerPE , gst(gstID)%maxMyLevCount,  &
                         gst(gstID)%nk, pre_specTransMpiType, sendtype, ierr);
    call mpi_type_create_resized(sendtype, lowerBound, extent, gst(gstID)%sendType_LonToLev, ierr);
    call mpi_type_commit(gst(gstID)%sendType_LonToLev, ierr);

    ! create the recv type for LonToLev
    extent = gst(gstID)%maxMyLevCount * gst(gstID)%lonPerPE * realSize
    call mpi_type_vector(gst(gstID)%latPerPE, gst(gstID)%maxMyLevCount * gst(gstID)%lonPerPE,  &
                         gst(gstID)%maxMyLevCount * gst(gstID)%ni, pre_specTransMpiType, recvtype, ierr)
    call mpi_type_create_resized(recvtype, lowerBound , extent, gst(gstID)%recvType_LonToLev, ierr);
    call mpi_type_commit(gst(gstID)%recvType_LonToLev,ierr)

    gst_setup = gstID
  
  end function GST_SETUP

!-------------------------------------------------------------------------------
! Data transposes with respect to 1 of 2 dimensions (i.e. NS or EW)
!-------------------------------------------------------------------------------

  subroutine transpose2d_NtoLev(psp_in,psp_out)
    implicit none

    ! Arguments:
    real(8) :: psp_in (gst(gstID)%myNla, 2, gst(gstID)%nk)
    real(8) :: psp_out(gst(gstID)%nla, 2, &
                       gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)

    ! Locals:
    real(pre_specTransReal) :: sp_send(gst(gstID)%maxMyNla, 2, &
                                       gst(gstID)%maxMyLevCount, mpi_npex)
    real(pre_specTransReal) :: sp_recv(gst(gstID)%maxMyNla, 2, &
                                       gst(gstID)%maxMyLevCount, mpi_npex)
    integer :: yourid,ila,icount,nsize,ierr,jlev,jlev2

    call utl_tmg_start(152,'low-level--gst_barr')
    if(mpi_doBarrier) call rpn_comm_barrier('EW',ierr)
    call tmg_stop(152)

    call utl_tmg_start(153,'low-level--gst_transpose_NtoLEV')

    !$OMP PARALLEL DO PRIVATE(yourid,jlev,jlev2,icount)
    do yourid = 0, (mpi_npex-1)
      do jlev = gst(gstID)%allLevBeg(yourid+1), gst(gstID)%allLevEnd(yourid+1)
        jlev2 = jlev - gst(gstID)%allLevBeg(yourid+1) + 1
        do icount = 1, gst(gstID)%myNla
          sp_send(icount,1,jlev2,yourid+1) = psp_in(icount,1,jlev)
          sp_send(icount,2,jlev2,yourid+1) = psp_in(icount,2,jlev)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    nsize = gst(gstID)%maxMyNla * 2 * gst(gstID)%maxMyLevCount
    if(mpi_npex.gt.1) then
      call rpn_comm_alltoall(sp_send, nsize, pre_specTransMpiReal,  &
                             sp_recv, nsize, pre_specTransMpiReal, 'EW', ierr)
    else
      sp_recv(:,:,:,1) = sp_send(:,:,:,1)
    endif

    !$OMP PARALLEL DO PRIVATE(yourid,jlev,jlev2,icount,ila)
    do yourid = 0, (mpi_npex-1)
      do jlev = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
        jlev2 = jlev - gst(gstID)%myLevBeg + 1
        do icount = 1, gst(gstID)%allNla(yourid+1)
          ila = gst(gstID)%allIlaList(icount,yourid+1)
          psp_out(ila,1,jlev) = sp_recv(icount,1,jlev2,yourid+1)
          psp_out(ila,2,jlev) = sp_recv(icount,2,jlev2,yourid+1)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    call tmg_stop(153)

  end subroutine transpose2d_NtoLev


  subroutine transpose2d_LevtoN(psp_in,psp_out)
    implicit none

    ! Arguments:
    real(8) :: psp_in (gst(gstID)%nla, 2, &
                       gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: psp_out(gst(gstID)%myNla, 2, gst(gstID)%nk)

    ! Locals:
    real(pre_specTransReal) :: sp_send(gst(gstID)%maxMyNla, 2, &
                                       gst(gstID)%maxMyLevCount, mpi_npex)
    real(pre_specTransReal) :: sp_recv(gst(gstID)%maxMyNla, 2, &
                                       gst(gstID)%maxMyLevCount, mpi_npex)
    integer :: yourid,ila,icount,nsize,ierr,jlev,jlev2

    call utl_tmg_start(152,'low-level--gst_barr')
    if(mpi_doBarrier) call rpn_comm_barrier('EW',ierr)
    call tmg_stop(152)

    call utl_tmg_start(153,'low-level--gst_transpose_NtoLEV')

    !$OMP PARALLEL DO PRIVATE(yourid,jlev,jlev2,icount,ila)
    do yourid = 0, (mpi_npex-1)
      do jlev = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
        jlev2 = jlev - gst(gstID)%myLevBeg + 1
        sp_send(:,:,jlev2,yourid+1) = 0.0d0
        do icount = 1, gst(gstID)%allNla(yourid+1)
          ila = gst(gstID)%allIlaList(icount,yourid+1)
          sp_send(icount,1,jlev2,yourid+1) = psp_in(ila,1,jlev)
          sp_send(icount,2,jlev2,yourid+1) = psp_in(ila,2,jlev)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    nsize = gst(gstID)%maxMyNla * 2 * gst(gstID)%maxMyLevCount
    if(mpi_npex.gt.1) then
      call rpn_comm_alltoall(sp_send, nsize, pre_specTransMpiReal,  &
                             sp_recv, nsize, pre_specTransMpiReal, 'EW', ierr)
    else
      sp_recv(:,:,:,1) = sp_send(:,:,:,1)
    endif

    !$OMP PARALLEL DO PRIVATE(yourid,jlev,jlev2,icount)
    do yourid = 0, (mpi_npex-1)
      do jlev = gst(gstID)%allLevBeg(yourid+1), gst(gstID)%allLevEnd(yourid+1)
        jlev2 = jlev - gst(gstID)%allLevBeg(yourid+1) + 1
        do icount = 1, gst(gstID)%myNla
          psp_out(icount,1,jlev) = sp_recv(icount,1,jlev2,yourid+1)
          psp_out(icount,2,jlev) = sp_recv(icount,2,jlev2,yourid+1)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    call tmg_stop(153)

  end subroutine transpose2d_LevtoN


  subroutine transpose2d_MtoLat(pgd_in,pgd_out)
    implicit none

    ! Arguments:
    real(8) :: pgd_in(2*gst(gstID)%maxmCount, gst(gstID)%nj, &
                        gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: pgd_out(gst(gstID)%ni, gst(gstID)%myLatBeg:gst(gstID)%myLatEnd, &
                       gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)

    ! Locals:
    real(pre_specTransReal) :: gd_send(gst(gstID)%maxmCount, 2, gst(gstID)%latPerPEmax, &
                                       gst(gstID)%maxMyLevCount, mpi_npey)
    real(pre_specTransReal) :: gd_recv(gst(gstID)%maxmCount, 2, gst(gstID)%latPerPEmax, &
                                       gst(gstID)%maxMyLevCount, mpi_npey)
    integer :: yourid,jm,jm2,icount,nsize,ierr,jlev,jlev2,jlat,jlat2

    call utl_tmg_start(152,'low-level--gst_barr')
    if(mpi_doBarrier) call rpn_comm_barrier('NS',ierr)
    call tmg_stop(152)

    call utl_tmg_start(154,'low-level--gst_transpose_MtoLAT')

    !$OMP PARALLEL DO PRIVATE(yourid,jlat,jlat2,jlev,jlev2,icount,jm,jm2)
    do yourid = 0, (mpi_npey-1)
      do jlev = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
        jlev2 = jlev - gst(gstID)%myLevBeg + 1
        gd_send(:, :, :, jlev2, yourid+1) = 0.0d0
        do jlat = gst(gstID)%allLatBeg(yourid+1), gst(gstID)%allLatEnd(yourid+1)
          jlat2 = jlat - gst(gstID)%allLatBeg(yourid+1) + 1
          icount = 0
          do jm = gst(gstID)%mymBeg, gst(gstID)%mymEnd, gst(gstID)%mymSkip
            jm2 = 2*gst(gstID)%mymIndex(jm)
            icount = icount + 1
            gd_send(icount, 1, jlat2, jlev2, yourid+1) = pgd_in(jm2-1, jlat, jlev)
            gd_send(icount, 2, jlat2, jlev2, yourid+1) = pgd_in(jm2  , jlat, jlev)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    nsize = gst(gstID)%maxmCount * 2 * gst(gstID)%maxMyLevCount * gst(gstID)%latPerPEmax
    if(mpi_npey.gt.1) then
      call rpn_comm_alltoall(gd_send, nsize, pre_specTransMpiReal,  &
                             gd_recv, nsize, pre_specTransMpiReal, 'NS', ierr)
    else
      gd_recv(:,:,:,:,1) = gd_send(:,:,:,:,1)
    endif

    !$OMP PARALLEL DO PRIVATE(yourid,jlat,jlat2,jlev,jlev2,icount,jm,jm2)
    do yourid = 0, (mpi_npey-1)
      do jlev = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
        jlev2 = jlev - gst(gstID)%myLevBeg + 1
        do jlat = gst(gstID)%myLatBeg, gst(gstID)%myLatEnd
          jlat2 = jlat - gst(gstID)%myLatBeg + 1
          icount = 0
          do jm = gst(gstID)%allmBeg(yourid+1), gst(gstID)%allmEnd(yourid+1), gst(gstID)%allmSkip(yourid+1)
            jm2 = 2*jm+1
            icount = icount + 1
            pgd_out(jm2  , jlat, jlev) = gd_recv(icount, 1, jlat2, jlev2, yourid+1)
            pgd_out(jm2+1, jlat, jlev) = gd_recv(icount, 2, jlat2, jlev2, yourid+1)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    call tmg_stop(154)

  end subroutine transpose2d_MtoLat


  subroutine transpose2d_MtoLat_kij(pgd_in,pgd_out)
    implicit none

    ! Arguments:
    real(8) :: pgd_in(gst(gstID)%maxMyLevCount, 2*gst(gstID)%maxmCount, &
                      gst(gstID)%nj)
    real(8) :: pgd_out(gst(gstID)%maxMyLevCount, gst(gstID)%ni, &
                       gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)

    ! Locals:
    real(pre_specTransReal), allocatable, save :: gd_send(:,:,:,:,:)
    real(pre_specTransReal), allocatable, save :: gd_recv(:,:,:,:,:)
    integer :: yourid,jm,jm2,icount,nsize,ierr,jlev,jlat,jlat2

    call utl_reAllocate(gd_send, gst(gstID)%maxMyLevCount, gst(gstID)%maxmCount, 2, gst(gstID)%latPerPEmax, mpi_npey)
    call utl_reAllocate(gd_recv, gst(gstID)%maxMyLevCount, gst(gstID)%maxmCount, 2, gst(gstID)%latPerPEmax, mpi_npey)

    call utl_tmg_start(152,'low-level--gst_barr')
    if(mpi_doBarrier) call rpn_comm_barrier('NS',ierr)
    call tmg_stop(152)

    call utl_tmg_start(154,'low-level--gst_transpose_MtoLAT')

    !$OMP PARALLEL DO PRIVATE(yourid,jlat,jlat2,jlev,icount,jm,jm2)
    do yourid = 0, (mpi_npey-1)
      gd_send(:, :, :, :, yourid+1) = 0.0d0
      do jlat = gst(gstID)%allLatBeg(yourid+1), gst(gstID)%allLatEnd(yourid+1)
        jlat2 = jlat - gst(gstID)%allLatBeg(yourid+1) + 1
        icount = 0
        do jm = gst(gstID)%mymBeg, gst(gstID)%mymEnd, gst(gstID)%mymSkip
          jm2 = 2*gst(gstID)%mymIndex(jm)
          icount = icount + 1
          do jlev = 1, gst(gstID)%maxMyLevCount
            gd_send(jlev, icount, 1, jlat2, yourid+1) = pgd_in(jlev, jm2-1, jlat)
            gd_send(jlev, icount, 2, jlat2, yourid+1) = pgd_in(jlev, jm2  , jlat)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    nsize = gst(gstID)%maxmCount * 2 * gst(gstID)%maxMyLevCount * gst(gstID)%latPerPEmax
    if(mpi_npey.gt.1) then
      call rpn_comm_alltoall(gd_send, nsize, pre_specTransMpiReal,  &
                             gd_recv, nsize, pre_specTransMpiReal, 'NS', ierr)
    else
      gd_recv(:,:,:,:,1) = gd_send(:,:,:,:,1)
    endif

    !$OMP PARALLEL DO PRIVATE(yourid,jlat,jlat2,jlev,icount,jm,jm2)
    do yourid = 0, (mpi_npey-1)
      do jlat = gst(gstID)%myLatBeg, gst(gstID)%myLatEnd
        jlat2 = jlat - gst(gstID)%myLatBeg + 1
        icount = 0
        do jm = gst(gstID)%allmBeg(yourid+1), gst(gstID)%allmEnd(yourid+1), gst(gstID)%allmSkip(yourid+1)
          jm2 = 2*jm+1
          icount = icount + 1
          do jlev = 1, gst(gstID)%maxMyLevCount
            pgd_out(jlev, jm2  , jlat) = gd_recv(jlev, icount, 1, jlat2, yourid+1)
            pgd_out(jlev, jm2+1, jlat) = gd_recv(jlev, icount, 2, jlat2, yourid+1)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    call tmg_stop(154)

  end subroutine transpose2d_MtoLat_kij


  subroutine transpose2d_LattoM(pgd_in,pgd_out)
    implicit none

    ! Arguments:
    real(8) :: pgd_in(gst(gstID)%ni, gst(gstID)%myLatBeg:gst(gstID)%myLatEnd, &
                      gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: pgd_out(2*gst(gstID)%maxmCount, gst(gstID)%nj,  &
                         gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)

    ! Locals:
    real(pre_specTransReal) :: gd_send(gst(gstID)%maxmCount, 2, gst(gstID)%latPerPEmax, &
                                       gst(gstID)%maxMyLevCount, mpi_npey)
    real(pre_specTransReal) :: gd_recv(gst(gstID)%maxmCount, 2, gst(gstID)%latPerPEmax, &
                                       gst(gstID)%maxMyLevCount, mpi_npey)
    integer :: yourid,jm,jm2,icount,nsize,ierr,jlev,jlev2,jlat,jlat2

    call utl_tmg_start(152,'low-level--gst_barr')
    if(mpi_doBarrier) call rpn_comm_barrier('NS',ierr)
    call tmg_stop(152)

    call utl_tmg_start(154,'low-level--gst_transpose_MtoLAT')

    !$OMP PARALLEL DO PRIVATE(yourid,jlat,jlat2,jlev,jlev2,icount,jm,jm2)
    do yourid = 0, (mpi_npey-1)
      do jlev = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
        jlev2 = jlev - gst(gstID)%myLevBeg + 1
        gd_send(:, :, :, jlev2, yourid+1) = 0.0d0
        do jlat = gst(gstID)%myLatBeg, gst(gstID)%myLatEnd
          jlat2 = jlat - gst(gstID)%myLatBeg + 1
          icount = 0
          do jm = gst(gstID)%allmBeg(yourid+1), gst(gstID)%allmEnd(yourid+1), gst(gstID)%allmSkip(yourid+1)
            jm2 = 2*jm+1
            icount = icount + 1
            gd_send(icount, 1, jlat2, jlev2, yourid+1) = pgd_in(jm2  , jlat, jlev)
            gd_send(icount, 2, jlat2, jlev2, yourid+1) = pgd_in(jm2+1, jlat, jlev)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    nsize = gst(gstID)%maxmCount * 2 * gst(gstID)%maxMyLevCount * gst(gstID)%latPerPEmax
    if(mpi_npey.gt.1) then
      call rpn_comm_alltoall(gd_send, nsize, pre_specTransMpiReal,  &
                             gd_recv, nsize, pre_specTransMpiReal, 'NS', ierr)
    else
      gd_recv(:,:,:,:,1) = gd_send(:,:,:,:,1)
    endif

    !$OMP PARALLEL DO PRIVATE(yourid,jlat,jlat2,jlev,jlev2,icount,jm,jm2)
    do yourid = 0, (mpi_npey-1)
      do jlev = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
        jlev2 = jlev - gst(gstID)%myLevBeg + 1
        do jlat = gst(gstID)%allLatBeg(yourid+1), gst(gstID)%allLatEnd(yourid+1)
          jlat2 = jlat - gst(gstID)%allLatBeg(yourid+1) + 1
          icount = 0
          do jm = gst(gstID)%mymBeg, gst(gstID)%mymEnd, gst(gstID)%mymSkip
            jm2 = 2*gst(gstID)%mymIndex(jm)
            icount = icount + 1
            pgd_out(jm2-1,jlat,jlev) = gd_recv(icount,1,jlat2,jlev2,yourid+1)
            pgd_out(jm2  ,jlat,jlev) = gd_recv(icount,2,jlat2,jlev2,yourid+1)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    call tmg_stop(154)

  end subroutine transpose2d_LattoM


  subroutine transpose2d_LattoM_kij(pgd_in,pgd_out)
    implicit none

    ! Arguments:
    real(8) :: pgd_in(gst(gstID)%maxMyLevCount, gst(gstID)%ni, &
                          gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)
    real(8) :: pgd_out(gst(gstID)%maxMyLevCount, 2*gst(gstID)%maxmCount, &
                       gst(gstID)%nj)

    ! Locals:
    real(pre_specTransReal), allocatable, save :: gd_send(:,:,:,:,:)
    real(pre_specTransReal), allocatable, save :: gd_recv(:,:,:,:,:)
    integer :: yourid,jm,jm2,icount,nsize,ierr,jlev,jlat,jlat2

    call utl_reAllocate(gd_send, gst(gstID)%maxMyLevCount, gst(gstID)%maxmCount, 2, gst(gstID)%latPerPEmax, mpi_npey)
    call utl_reAllocate(gd_recv, gst(gstID)%maxMyLevCount, gst(gstID)%maxmCount, 2, gst(gstID)%latPerPEmax, mpi_npey)

    call utl_tmg_start(152,'low-level--gst_barr')
    if(mpi_doBarrier) call rpn_comm_barrier('NS',ierr)
    call tmg_stop(152)

    call utl_tmg_start(154,'low-level--gst_transpose_MtoLAT')

    !$OMP PARALLEL DO PRIVATE(yourid,jlat,jlat2,jlev,icount,jm,jm2)
    do yourid = 0, (mpi_npey-1)
      gd_send(:, :, :, :, yourid+1) = 0.0d0
      do jlat = gst(gstID)%myLatBeg, gst(gstID)%myLatEnd
        jlat2 = jlat - gst(gstID)%myLatBeg + 1
        icount = 0
        do jm = gst(gstID)%allmBeg(yourid+1), gst(gstID)%allmEnd(yourid+1), gst(gstID)%allmSkip(yourid+1)
          jm2 = 2*jm+1
          icount = icount + 1
          do jlev = 1, gst(gstID)%maxMyLevCount
            gd_send(jlev, icount, 1, jlat2, yourid+1) = pgd_in(jlev, jm2  , jlat)
            gd_send(jlev, icount, 2, jlat2, yourid+1) = pgd_in(jlev, jm2+1, jlat)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    nsize = gst(gstID)%maxmCount * 2 * gst(gstID)%maxMyLevCount * gst(gstID)%latPerPEmax
    if(mpi_npey.gt.1) then
      call rpn_comm_alltoall(gd_send, nsize, pre_specTransMpiReal,  &
                             gd_recv, nsize, pre_specTransMpiReal, 'NS', ierr)
    else
      gd_recv(:,:,:,:,1) = gd_send(:,:,:,:,1)
    endif

    !$OMP PARALLEL DO PRIVATE(yourid,jlat,jlat2,jlev,icount,jm,jm2)
    do yourid = 0, (mpi_npey-1)
      do jlat = gst(gstID)%allLatBeg(yourid+1), gst(gstID)%allLatEnd(yourid+1)
        jlat2 = jlat - gst(gstID)%allLatBeg(yourid+1) + 1
        icount = 0
        do jm = gst(gstID)%mymBeg, gst(gstID)%mymEnd, gst(gstID)%mymSkip
          jm2 = 2*gst(gstID)%mymIndex(jm)
          icount = icount + 1
          do jlev = 1, gst(gstID)%maxMyLevCount
            pgd_out(jlev, jm2-1, jlat) = gd_recv(jlev, icount, 1, jlat2, yourid+1)
            pgd_out(jlev, jm2  , jlat) = gd_recv(jlev, icount, 2, jlat2, yourid+1)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    call tmg_stop(154)

  end subroutine transpose2d_LattoM_kij


  subroutine transpose2d_LevtoLon(pgd_in,pgd_out)
    implicit none

    ! Arguments:
    real(8) :: pgd_in(gst(gstID)%ni, gst(gstID)%myLatBeg:gst(gstID)%myLatEnd, &
                      gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: pgd_out(gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                       gst(gstID)%myLatBeg:gst(gstID)%myLatEnd, gst(gstID)%nk)

    ! Locals:
    real(pre_specTransReal) :: gd_send(gst(gstID)%lonPerPEmax, gst(gstID)%latPerPEmax, &
                                       gst(gstID)%maxMyLevCount, mpi_npex)
    real(pre_specTransReal) :: gd_recv(gst(gstID)%lonPerPEmax, gst(gstID)%latPerPEmax, &
                                       gst(gstID)%maxMyLevCount, mpi_npex)
    integer :: youridP1,nsize,ierr,jlev,jlev2

    call utl_tmg_start(152,'low-level--gst_barr')
    if(mpi_doBarrier) call rpn_comm_barrier('EW',ierr)
    call tmg_stop(152)

    call utl_tmg_start(155,'low-level--gst_transpose_LEVtoLON')

    !$OMP PARALLEL DO PRIVATE(youridP1,jlev,jlev2)
    do jlev = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
      jlev2 = jlev - gst(gstID)%myLevBeg + 1
      do youridP1 = 1, mpi_npex
        gd_send(:, :, jlev2, youridP1) =  0.0d0
        gd_send(1:gst(gstID)%allLonPerPE(youridP1), 1:gst(gstID)%latPerPE, jlev2, youridP1) =  &
          pgd_in(gst(gstID)%allLonBeg(youridP1):gst(gstID)%allLonEnd(youridP1), :, jlev)
      enddo
    enddo
    !$OMP END PARALLEL DO

    nsize = gst(gstID)%lonPerPEmax * gst(gstID)%maxMyLevCount * gst(gstID)%latPerPEmax
    if(mpi_npex.gt.1) then
      call rpn_comm_alltoall(gd_send,nsize,pre_specTransMpiReal,  &
                             gd_recv,nsize,pre_specTransMpiReal,'EW',ierr)
    else
      gd_recv(:,:,:,1) = gd_send(:,:,:,1)
    endif

    !$OMP PARALLEL DO PRIVATE(youridP1,jlev,jlev2)
    do youridP1 = 1, mpi_npex
      do jlev=gst(gstID)%allLevBeg(youridP1),gst(gstID)%allLevEnd(youridP1)
        jlev2=jlev-gst(gstID)%allLevBeg(youridP1)+1
        pgd_out(:, :, jlev) = gd_recv(1:gst(gstID)%lonPerPE, 1:gst(gstID)%latPerPE, jlev2, youridP1)
      enddo
    enddo
    !$OMP END PARALLEL DO

    call tmg_stop(155)

  end subroutine transpose2d_LevtoLon


  subroutine transpose2d_LevtoLon_kij_mpitypes8(pgd_in,pgd_out)
    implicit none

    ! Arguments:
    real(8) :: pgd_in(gst(gstID)%maxMyLevCount, gst(gstID)%ni, &
                      gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)
    real(8) :: pgd_out(gst(gstID)%nk, gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                       gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)

    ! Locals:
    integer :: nsize,ierr

    call utl_tmg_start(152,'low-level--gst_barr')
    if(mpi_doBarrier) call rpn_comm_barrier('EW',ierr)
    call tmg_stop(152)

    call utl_tmg_start(155,'low-level--gst_transpose_LEVtoLON')

    nsize = gst(gstID)%lonPerPE * gst(gstID)%maxMyLevCount * gst(gstID)%latPerPE
    if(mpi_npex.gt.1) then
      call mpi_alltoall(pgd_in,  1, gst(gstID)%sendType_LevToLon,  &
                        pgd_out, 1, gst(gstID)%recvType_LevToLon, mpi_comm_EW, ierr)
    else
      pgd_out(:,:,:) = pgd_in(:,:,:)
    endif

    call tmg_stop(155)

  end subroutine transpose2d_LevtoLon_kij_mpitypes8


  subroutine transpose2d_LevtoLon_kij_mpitypes4(pgd_in,pgd_out)
    implicit none

    ! Arguments:
    real(8) :: pgd_in(gst(gstID)%maxMyLevCount, gst(gstID)%ni, &
                      gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)
    real(8) :: pgd_out(gst(gstID)%nk, gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                       gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)

    ! Locals:
    integer :: nsize,ierr
    real(4) :: pgd_in_r4(gst(gstID)%maxMyLevCount, gst(gstID)%ni, &
                         gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)
    real(4) :: pgd_out_r4(gst(gstID)%nk, gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                          gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)

    call utl_tmg_start(152,'low-level--gst_barr')
    if(mpi_doBarrier) call rpn_comm_barrier('EW',ierr)
    call tmg_stop(152)

    call utl_tmg_start(155,'low-level--gst_transpose_LEVtoLON')

    nsize = gst(gstID)%lonPerPE * gst(gstID)%maxMyLevCount * gst(gstID)%latPerPE
    if(mpi_npex.gt.1) then
      pgd_in_r4(:,:,:) = pgd_in(:,:,:)
      call mpi_alltoall(pgd_in_r4,  1, gst(gstID)%sendType_LevToLon,  &
                        pgd_out_r4, 1, gst(gstID)%recvType_LevToLon, mpi_comm_EW, ierr)
      pgd_out(:,:,:) = pgd_out_r4(:,:,:)
    else
      pgd_out(:,:,:) = pgd_in(:,:,:)
    endif

    call tmg_stop(155)

  end subroutine transpose2d_LevtoLon_kij_mpitypes4


  subroutine transpose2d_LevtoLon_kij(pgd_in,pgd_out)
    implicit none

    ! Arguments:
    real(8) :: pgd_in(gst(gstID)%maxMyLevCount, gst(gstID)%ni, &
                      gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)
    real(8) :: pgd_out(gst(gstID)%nk, gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                       gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)

    ! Locals:
    real(pre_specTransReal) :: gd_send(gst(gstID)%maxMyLevCount, gst(gstID)%lonPerPEmax,&
                                       gst(gstID)%latPerPEmax, mpi_npex)
    real(pre_specTransReal) :: gd_recv(gst(gstID)%maxMyLevCount, gst(gstID)%lonPerPEmax,&
                                       gst(gstID)%latPerPEmax, mpi_npex)
    integer :: youridP1, nsize, ierr, yourNumLev

    call utl_tmg_start(152,'low-level--gst_barr')
    if(mpi_doBarrier) call rpn_comm_barrier('EW',ierr)
    call tmg_stop(152)

    call utl_tmg_start(155,'low-level--gst_transpose_LEVtoLON')

    !$OMP PARALLEL DO PRIVATE(youridP1)
    do youridP1 = 1, mpi_npex
      gd_send(:, :, :, youridP1) = 0.0d0
      gd_send(:, 1:gst(gstID)%allLonPerPE(youridP1), 1:gst(gstID)%latPerPE, youridP1) =  &
        pgd_in(:, gst(gstID)%allLonBeg(youridP1):gst(gstID)%allLonEnd(youridP1), :)
    enddo
    !$OMP END PARALLEL DO

    nsize = gst(gstID)%lonPerPEmax * gst(gstID)%maxMyLevCount * gst(gstID)%latPerPEmax
    if(mpi_npex.gt.1) then
      call rpn_comm_alltoall(gd_send, nsize, pre_specTransMpiReal,   &
                             gd_recv, nsize, pre_specTransMpiReal, 'EW', ierr)
    else
      gd_recv(:,:,:,1) = gd_send(:,:,:,1)
    endif

    !$OMP PARALLEL DO PRIVATE(youridP1,yourNumLev)
    do youridP1 = 1, mpi_npex
      yourNumLev = gst(gstID)%allLevEnd(youridP1) - gst(gstID)%allLevBeg(youridP1) + 1
      pgd_out(gst(gstID)%allLevBeg(youridP1):gst(gstID)%allLevEnd(youridP1), :, :) =  &
           gd_recv(1:yourNumLev, 1:gst(gstID)%lonPerPE, 1:gst(gstID)%latPerPE, youridP1)
    enddo
    !$OMP END PARALLEL DO

    call tmg_stop(155)

  end subroutine transpose2d_LevtoLon_kij


  subroutine transpose2d_LontoLev(pgd_in,pgd_out)
    implicit none

    ! Arguments:
    real(8) :: pgd_in(gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                      gst(gstID)%myLatBeg:gst(gstID)%myLatEnd, gst(gstID)%nk)
    real(8) :: pgd_out(gst(gstID)%ni, gst(gstID)%myLatBeg:gst(gstID)%myLatEnd, &
                       gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)

    ! Locals:
    real(pre_specTransReal) :: gd_send(gst(gstID)%lonPerPEmax, gst(gstID)%latPerPEmax, &
                                       gst(gstID)%maxMyLevCount, mpi_npex)
    real(pre_specTransReal) :: gd_recv(gst(gstID)%lonPerPEmax, gst(gstID)%latPerPEmax, &
                                       gst(gstID)%maxMyLevCount, mpi_npex)
    integer :: youridP1,nsize,ierr,jlev,jlev2

    call utl_tmg_start(152,'low-level--gst_barr')
    if(mpi_doBarrier) call rpn_comm_barrier('EW',ierr)
    call tmg_stop(152)

    call utl_tmg_start(155,'low-level--gst_transpose_LEVtoLON')

    !$OMP PARALLEL DO PRIVATE(youridP1,jlev,jlev2)
    do youridP1 = 1, mpi_npex
      do jlev=gst(gstID)%allLevBeg(youridP1),gst(gstID)%allLevEnd(youridP1)
        jlev2=jlev-gst(gstID)%allLevBeg(youridP1)+1
        gd_send(:, :, jlev2, youridP1) = 0.0d0
        gd_send(1:gst(gstID)%lonPerPE, 1:gst(gstID)%latPerPE, jlev2, youridP1) =  &
             pgd_in(:,:,jlev)
      enddo
    enddo
    !$OMP END PARALLEL DO

    nsize = gst(gstID)%lonPerPEmax * gst(gstID)%maxMyLevCount * gst(gstID)%latPerPEmax
    if(mpi_npex.gt.1) then
      call rpn_comm_alltoall(gd_send, nsize, pre_specTransMpiReal,  &
                             gd_recv, nsize, pre_specTransMpiReal, 'EW', ierr)
    else
      gd_recv(:,:,:,1) = gd_send(:,:,:,1)
    endif

    !$OMP PARALLEL DO PRIVATE(youridP1,jlev,jlev2)
    do jlev = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
      jlev2 = jlev - gst(gstID)%myLevBeg + 1
      do youridP1 = 1, mpi_npex
        pgd_out(gst(gstID)%allLonBeg(youridP1):gst(gstID)%allLonEnd(youridP1), :, jlev) =  &
          gd_recv(1:gst(gstID)%allLonPerPE(youridP1),1:gst(gstID)%latPerPE,jlev2,youridP1)
      enddo
    enddo
    !$OMP END PARALLEL DO

    call tmg_stop(155)

  end subroutine transpose2d_LontoLev


  subroutine transpose2d_LontoLev_kij_mpitypes8(pgd_in,pgd_out)
    implicit none

    ! Arguments:
    real(8) :: pgd_in(gst(gstID)%nk, gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                      gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)
    real(8) :: pgd_out(gst(gstID)%maxMyLevCount, gst(gstID)%ni, &
                       gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)

    ! Locals:
    integer :: nsize, ierr

    call utl_tmg_start(152,'low-level--gst_barr')
    if(mpi_doBarrier) call rpn_comm_barrier('EW',ierr)
    call tmg_stop(152)

    call utl_tmg_start(155,'low-level--gst_transpose_LEVtoLON')

    nsize = gst(gstID)%lonPerPE * gst(gstID)%maxMyLevCount * gst(gstID)%latPerPE
    if(mpi_npex.gt.1) then
      call mpi_alltoall(pgd_in,  1, gst(gstID)%sendType_LonToLev,  &
                        pgd_out, 1, gst(gstID)%recvType_LonToLev, mpi_comm_EW, ierr)
    else
      pgd_out(:,:,:) = pgd_in(:,:,:)
    endif

    call tmg_stop(155)

  end subroutine transpose2d_LontoLev_kij_mpitypes8


  subroutine transpose2d_LontoLev_kij_mpitypes4(pgd_in,pgd_out)
    implicit none

    ! Arguments:
    real(8) :: pgd_in(gst(gstID)%nk, gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                      gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)
    real(8) :: pgd_out(gst(gstID)%maxMyLevCount, gst(gstID)%ni, &
                       gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)

    ! Locals:
    integer :: nsize, ierr
    real(4) :: pgd_in_r4(gst(gstID)%nk, gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                         gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)
    real(4) :: pgd_out_r4(gst(gstID)%maxMyLevCount, gst(gstID)%ni, &
                          gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)

    call utl_tmg_start(152,'low-level--gst_barr')
    if(mpi_doBarrier) call rpn_comm_barrier('EW',ierr)
    call tmg_stop(152)

    call utl_tmg_start(155,'low-level--gst_transpose_LEVtoLON')

    nsize = gst(gstID)%lonPerPE * gst(gstID)%maxMyLevCount * gst(gstID)%latPerPE
    if(mpi_npex.gt.1) then
      pgd_in_r4(:,:,:) = pgd_in(:,:,:)
      call mpi_alltoall(pgd_in_r4,  1, gst(gstID)%sendType_LonToLev,  &
                        pgd_out_r4, 1, gst(gstID)%recvType_LonToLev, mpi_comm_EW, ierr)
      pgd_out(:,:,:) = pgd_out_r4(:,:,:)
    else
      pgd_out(:,:,:) = pgd_in(:,:,:)
    endif

    call tmg_stop(155)

  end subroutine transpose2d_LontoLev_kij_mpitypes4


  subroutine transpose2d_LontoLev_kij(pgd_in,pgd_out)
    implicit none

    ! Arguments:
    real(8) :: pgd_in(gst(gstID)%nk, gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                      gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)
    real(8) :: pgd_out(gst(gstID)%maxMyLevCount, gst(gstID)%ni, &
                       gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)

    ! Locals:
    real(pre_specTransReal) :: gd_send(gst(gstID)%maxMyLevCount, gst(gstID)%lonPerPEmax,&
                                       gst(gstID)%latPerPEmax, mpi_npex)
    real(pre_specTransReal) :: gd_recv(gst(gstID)%maxMyLevCount, gst(gstID)%lonPerPEmax,&
                                       gst(gstID)%latPerPEmax, mpi_npex)
    integer :: youridP1,nsize,ierr,jlev,jlev2,yourNumLev

    call utl_tmg_start(152,'low-level--gst_barr')
    if(mpi_doBarrier) call rpn_comm_barrier('EW',ierr)
    call tmg_stop(152)

    call utl_tmg_start(155,'low-level--gst_transpose_LEVtoLON')

    !$OMP PARALLEL DO PRIVATE(youridP1,yourNumLev)
    do youridP1 = 1, mpi_npex
      yourNumLev = gst(gstID)%allLevEnd(youridP1) - gst(gstID)%allLevBeg(youridP1) + 1
      gd_send(:, :, :, youridP1) = 0.0d0
      gd_send(1:yourNumLev, 1:gst(gstID)%lonPerPE, 1:gst(gstID)%latPerPE, youridP1) =  &
           pgd_in(gst(gstID)%allLevBeg(youridP1):gst(gstID)%allLevEnd(youridP1),:,:)
    enddo
    !$OMP END PARALLEL DO

    nsize = gst(gstID)%lonPerPEmax * gst(gstID)%maxMyLevCount * gst(gstID)%latPerPEmax
    if(mpi_npex.gt.1) then
      call rpn_comm_alltoall(gd_send, nsize, pre_specTransMpiReal,  &
                             gd_recv, nsize, pre_specTransMpiReal, 'EW', ierr)
    else
      gd_recv(:,:,:,1) = gd_send(:,:,:,1)
    endif

    !$OMP PARALLEL DO PRIVATE(youridP1,jlev,jlev2)
    do youridP1 = 1, mpi_npex
      pgd_out(:, gst(gstID)%allLonBeg(youridP1):gst(gstID)%allLonEnd(youridP1), :) =  &
        gd_recv(:, 1:gst(gstID)%allLonPerPE(youridP1), 1:gst(gstID)%latPerPE, youridP1)
    enddo
    !$OMP END PARALLEL DO

    call tmg_stop(155)

  end subroutine transpose2d_LontoLev_kij

!-------------------------------------------------------------------------------
! Subroutines to re-order the u and v wind components for mpi version of spgd
! and spgda
!-------------------------------------------------------------------------------

  subroutine interleaveWinds_sp(psp,nflev)
    implicit none

    ! Arguments:
    integer :: nflev
    real(8) :: psp(gst(gstID)%myNla,2,gst(gstID)%nk)

    ! Locals:
    real(8) :: tempvalues(2,nflev*2)
    integer :: jk, ila

    !$OMP PARALLEL DO PRIVATE (ILA,JK,TEMPVALUES)
    do ila = 1, gst(gstID)%myNla
       do jk = 1, nflev
          ! place u in new position in temporary array
          tempvalues(:,(jk*2)-1) = psp(ila,:,jk)
          ! place v in new position in temporary array
          tempvalues(:,jk*2)     = psp(ila,:,jk+nflev)
       enddo
       ! move contents of temporary array back to original array
       psp(ila,:,1:2*nflev) = tempvalues(:,1:2*nflev)
    enddo
    !$OMP END PARALLEL DO

  end subroutine interleaveWinds_sp


  subroutine unInterleaveWinds_sp(psp,nflev)
    implicit none

    ! Arguments:
    integer :: nflev
    real(8) :: psp(gst(gstID)%myNla,2,gst(gstID)%nk)

    ! Locals:
    real(8) :: tempvalues(2,nflev*2)
    integer :: jk, ila

    !$OMP PARALLEL DO PRIVATE (ILA,JK,TEMPVALUES)
    do ila = 1, gst(gstID)%myNla
       do jk = 1, nflev
          ! place u in original position in temporary array
          tempvalues(:,jk)       = psp(ila,:,(jk*2)-1)
          ! place v in original position in temporary array
          tempvalues(:,jk+nflev) = psp(ila,:,jk*2)
       enddo
       ! move contents of temporary array back to original array
       psp(ila,:,1:2*nflev) = tempvalues(:,1:2*nflev)
    enddo
    !$OMP END PARALLEL DO

  end subroutine unInterleaveWinds_sp


  subroutine interleaveWinds_gd(pgd,nflev)
    implicit none

    ! Arguments:
    integer :: nflev
    real(8) :: pgd(gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                   gst(gstID)%myLatBeg:gst(gstID)%myLatEnd, gst(gstID)%nk)

    ! Locals:
    real(8) :: tempvalues(nflev*2)
    integer :: jlat, jk, jlon

    !$OMP PARALLEL DO PRIVATE (JLAT,JLON,JK,TEMPVALUES)
    do jlat = gst(gstID)%myLatBeg, gst(gstID)%myLatEnd
       do jlon = gst(gstID)%myLonBeg, gst(gstID)%myLonEnd
          do jk = 1, nflev
             ! place u in original position in temporary array
             tempvalues((jk*2)-1) = pgd(jlon,jlat,jk)
             ! place v in original position in temporary array
             tempvalues(jk*2)     = pgd(jlon,jlat,jk+nflev)
          enddo
          ! move contents of temporary array back to original array
          pgd(jlon,jlat,1:2*nflev) = tempvalues(1:2*nflev)
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine interleaveWinds_gd


  subroutine unInterleaveWinds_gd(pgd,nflev)
    implicit none

    ! Arguments:
    integer :: nflev
    real(8) :: pgd(gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                   gst(gstID)%myLatBeg:gst(gstID)%myLatEnd, gst(gstID)%nk)

    ! Locals:
    real(8) :: tempvalues(nflev*2)
    integer :: jlat, jk, jlon

    !$OMP PARALLEL DO PRIVATE (JLAT,JLON,JK,TEMPVALUES)
    do jlat = gst(gstID)%myLatBeg, gst(gstID)%myLatEnd
       do jlon = gst(gstID)%myLonBeg, gst(gstID)%myLonEnd
          do jk = 1, nflev
             ! place u in original position in temporary array
             tempvalues(jk)       = pgd(jlon,jlat,(jk*2)-1)
             ! place v in original position in temporary array
             tempvalues(jk+nflev) = pgd(jlon,jlat,jk*2)
          enddo
          ! move contents of temporary array back to original array
          pgd(jlon,jlat,1:2*nflev) = tempvalues(1:2*nflev)
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine unInterleaveWinds_gd

!-------------------------------------------------------------------------------
! Main spectral transform subroutines
!-------------------------------------------------------------------------------

  subroutine gst_spgd(psp,pgd,nflev)
    implicit none

    ! Arguments:
    integer :: nflev
    real(8) :: psp(gst(gstID)%myNla,2,gst(gstID)%nk)
    real(8) :: pgd(gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                   gst(gstID)%myLatBeg:gst(gstID)%myLatEnd, gst(gstID)%nk)

    ! Locals:
    real(8), allocatable :: psp2(:,:,:),pgd2(:,:,:),pgd3(:,:,:)
    integer :: jlat, jk, jlon

    ! check if this mpi task will deal with winds during Legendre transform
    if(gst(gstID)%myLevBeg.le.2*nflev) then
      ! ensure that the number of levels on this mpi task is even to allow interleaving of u and v
      ! only necessary when number of levels on an mpi task is less than all wind levels (2nd condition)
      if( (mod(gst(gstID)%myLevCount,2).ne.0) .and. (gst(gstID)%myLevCount < 2*nflev) ) then
        write(*,*) 'GST_SPGD: myLevCount = ',gst(gstID)%myLevCount
        call utl_abort('GST_SPGD: Number of levels on this mpi task must be even!')
      endif
    endif

    allocate(psp2(gst(gstID)%nla, 2, gst(gstID)%myLevBeg:gst(gstID)%myLevEnd))
    allocate(pgd2(2*gst(gstID)%maxmcount, gst(gstID)%nj,  gst(gstID)%myLevBeg:gst(gstID)%myLevEnd))
    allocate(pgd3(gst(gstID)%ni, gst(gstID)%myLatBeg:gst(gstID)%myLatEnd, gst(gstID)%myLevBeg:gst(gstID)%myLevEnd))

    ! 1.0 First reorder wind components to have u and v for same level on same mpi task
    call interleaveWinds_sp(psp,nflev)

    ! 1.1 Transpose data along npex from N to Levels
    call transpose2d_NtoLev(psp,psp2)

    ! 1.2 Inverse Legendre transform
    call utl_tmg_start(150,'low-level--gst_lt')
    call spgdpar(psp2,pgd2,nflev)
    call tmg_stop(150)
    deallocate(psp2)

    ! 1.3 Transpose data along npey from M to Latitudes
    call transpose2d_MtoLat(pgd2,pgd3)
    deallocate(pgd2)

    !$OMP PARALLEL DO PRIVATE (JLAT,JLON,JK)
    ! 2.1 Reset to zero the modes that are not part of the truncation
    do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
      do jlat = gst(gstID)%myLatBeg, gst(gstID)%myLatEnd
        do jlon = 2*(gst(gstID)%ntrunc+1)+1, gst(gstID)%ni
          pgd3(jlon,jlat,jk) = 0.d0
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    ! 2.2 Apply the FFT 
    call utl_tmg_start(151,'low-level--gst_fft')
    call fft3dvar(pgd3,+1)
    call tmg_stop(151)

    ! 2.3 Transpose data along npex from Levels to Longitudes
    call transpose2d_LevtoLon(pgd3,pgd)
    deallocate(pgd3)

    ! 2.4 Now undo reordering of wind components 
    call unInterleaveWinds_gd(pgd,nflev)

  end subroutine gst_spgd


  ! FIRST ATTEMPT AT MODIFICATIONS FOR MPI
  subroutine gst_gdsp(psp,pgd,nflev)
    implicit none

    ! Arguments:
    integer :: nflev
    real(8) :: psp(gst(gstID)%myNla,2,gst(gstID)%nk)
    real(8) :: pgd(gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                   gst(gstID)%myLatBeg:gst(gstID)%myLatEnd,gst(gstID)%nk)

    ! Locals:
    real(8), allocatable :: psp2(:,:,:),pgd2(:,:,:),pgd3(:,:,:)

    ! check if this mpi task will deal with winds during Legendre transform
    if(gst(gstID)%myLevBeg.le.2*nflev) then
      ! ensure that the number of levels on this mpi task is even to allow interleaving of u and v
      ! only necessary when number of levels on an mpi task is less than all wind levels (2nd condition)
      if( (mod(gst(gstID)%myLevCount,2).ne.0) .and. (gst(gstID)%myLevCount < 2*nflev) ) then
        write(*,*) 'GST_GDSP: myLevCount = ',gst(gstID)%myLevCount
        call utl_abort('GST_GDSP: Number of levels on this mpi task must be even!')
      endif
    endif

    allocate(psp2(gst(gstID)%nla, 2, gst(gstID)%myLevBeg:gst(gstID)%myLevEnd))
    allocate(pgd2(2*gst(gstID)%maxmcount, gst(gstID)%nj, gst(gstID)%myLevBeg:gst(gstID)%myLevEnd))
    allocate(pgd3(gst(gstID)%ni, gst(gstID)%myLatBeg:gst(gstID)%myLatEnd, gst(gstID)%myLevBeg:gst(gstID)%myLevEnd))

    ! 1.0 First reorder wind components to have u and v for same level on same mpi task
    call interleaveWinds_gd(pgd,nflev)

    ! Transpose data along npex from Longitudes to Levels
    call transpose2d_LontoLev(pgd,pgd3)

    ! 1. Fourier transform all fields for all latitudes
    call utl_tmg_start(151,'low-level--gst_fft')
    call fft3dvar(pgd3,-1)
    call tmg_stop(151)

    ! Transpose data along npey from Latitudes to M
    call transpose2d_LattoM(pgd3,pgd2)
    deallocate(pgd3)

    ! 2. Direct Legendre transform including wind transformations
    call utl_tmg_start(150,'low-level--gst_lt')
    call gdsppar(psp2,pgd2,nflev)
    call tmg_stop(150)
    deallocate(pgd2)

    ! Transpose data along npex from Levels to N
    call transpose2d_LevtoN(psp2,psp)
    deallocate(psp2)

    ! 2.4 Now undo reordering of wind components 
    call unInterleaveWinds_sp(psp,nflev)

  end subroutine gst_gdsp


  subroutine spgdpar(psp,pgd2,nflev)
    !
    !:Purpose: Inverse spectral transform(PARALLEL LOOP)
    implicit none

    ! Arguments:
    integer :: nflev
    real(8) :: psp(gst(gstID)%nla,2,gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: pgd2(2*gst(gstID)%maxmcount,gst(gstID)%nj,&
                      gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)

    ! Locals
    integer :: jj, jj2, jm, jn, ilonr, iloni, jk, jk2, ila, inm

    real(8) :: zjm, factor
    real(8) :: dlalp(0:gst(gstID)%ntrunc,gst(gstID)%njlath)
    real(8) :: dldalp(0:gst(gstID)%ntrunc,gst(gstID)%njlath)
    real(8) :: zfms(gst(gstID)%njlath+1,2,gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: zfma(gst(gstID)%njlath+1,2,gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: dlsp(0:gst(gstID)%ntrunc,2,gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)

    ! Inverse Legendre transform

    !$OMP PARALLEL DO PRIVATE(DLALP,DLDALP,DLSP,ZFMS,ZFMA, &
    !$OMP INM,ILA,JM,JN,JK,JK2,JJ,JJ2,ZJM,ILONR,ILONI,FACTOR)
    do jm = gst(gstID)%mymBeg, gst(gstID)%mymEnd, gst(gstID)%mymSkip

       ilonr = 2*gst(gstID)%mymIndex(jm)-1
       iloni = 2*gst(gstID)%mymIndex(jm)
       zjm = real(jm,8)

       ! 2.1 Copy global spectral state into local spectral state
       do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
          do jn = jm, gst(gstID)%ntrunc
             ila = gst(gstID)%nind(jm) + jn - jm
             inm = jn - jm
             if(jk > 2*nflev) then
                ! Scalar fields
                dlsp(inm,1,jk) = psp(ila,1,jk)
                dlsp(inm,2,jk) = psp(ila,2,jk)
             else
                ! Vector fields                
                dlsp(inm,1,jk) = psp(ila,1,jk)*gst(gstID)%r1snp1(ila)
                dlsp(inm,2,jk) = psp(ila,2,jk)*gst(gstID)%r1snp1(ila)
             endif
          enddo
       enddo

       ! 2.2  Get Legendre polynomial (and its derivative) for all latitudes
       !      but for the chosen value of "m" from the global array
       call getalp (dlalp,dldalp,gst(gstID)%njlath,gst(gstID)%ntrunc,gst(gstID)%ntrunc,jm)

       ! 2.3  Perform the inverse Legendre transform for all fields
       call leginv2d(jm,zfms,zfma,dlsp,dlalp,gst(gstID)%njlath,gst(gstID)%ntrunc,gst(gstID)%ntrunc)

       ! 2.4 Passage to Fourier space

       ! Scalar fields
       do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
          if(jk > 2*nflev) then
             do jj = 1, gst(gstID)%nj
                jj2 = gst(gstID)%nj - jj + 1
                if(jj.le.gst(gstID)%njlath) then
                   pgd2(ilonr,jj2,jk) = zfms(jj,1,jk) + zfma(jj,1,jk)
                   pgd2(iloni,jj2,jk) = zfms(jj,2,jk) + zfma(jj,2,jk)
                else
                   pgd2(ilonr,jj2,jk) = zfms(jj2,1,jk) - zfma(jj2,1,jk)
                   pgd2(iloni,jj2,jk) = zfms(jj2,2,jk) - zfma(jj2,2,jk)
                endif
             enddo
          endif
       enddo

       ! Vector fields: Note that u and v are interleaved in mode 5!
       do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd, 2
          if(jk <= 2*nflev) then
             jk2 = jk + 1  ! jk is u, jk2 is v
             do jj = 1, gst(gstID)%nj
                jj2 = gst(gstID)%nj - jj + 1
                if(jj.le.gst(gstID)%njlath) then
                   pgd2(ilonr,jj2,jk ) = -zjm*(zfms(jj,2,jk2) + zfma(jj,2,jk2))
                   pgd2(iloni,jj2,jk ) =  zjm*(zfms(jj,1,jk2) + zfma(jj,1,jk2))
                   pgd2(ilonr,jj2,jk2) = -zjm*(zfms(jj,2,jk)  + zfma(jj,2,jk))
                   pgd2(iloni,jj2,jk2) =  zjm*(zfms(jj,1,jk)  + zfma(jj,1,jk))
                else
                   pgd2(ilonr,jj2,jk)  = -zjm*(zfms(jj2,2,jk2) - zfma(jj2,2,jk2))
                   pgd2(iloni,jj2,jk)  =  zjm*(zfms(jj2,1,jk2) - zfma(jj2,1,jk2))
                   pgd2(ilonr,jj2,jk2) = -zjm*(zfms(jj2,2,jk)  - zfma(jj2,2,jk))
                   pgd2(iloni,jj2,jk2) =  zjm*(zfms(jj2,1,jk)  - zfma(jj2,1,jk))
                endif
             enddo
          endif
       enddo

       ! 2.5 Completion of the computation of the winds in Fourier space
       call leginv2d(jm,zfma,zfms,dlsp,dldalp,gst(gstID)%njlath,gst(gstID)%ntrunc,gst(gstID)%ntrunc)

       do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd, 2
          if(jk <= 2*nflev) then ! only for winds
             jk2 = jk + 1
             do jj = 1, gst(gstID)%nj
                jj2 = gst(gstID)%nj - jj + 1
                factor = ec_ra / cos(gst(gstID)%rlati(jj))
                if(jj.ne.jj2) then
                   ! For latitudes not exactly at equator
                   if(jj.le.gst(gstID)%njlath) then
                      ! northern latitudes
                      ! u component
                      pgd2(ilonr,jj2,jk)  = factor * ( pgd2(ilonr,jj2,jk)  - (zfms(jj,1,jk)  + zfma(jj,1,jk)) )
                      pgd2(iloni,jj2,jk)  = factor * ( pgd2(iloni,jj2,jk)  - (zfms(jj,2,jk)  + zfma(jj,2,jk)) )
                      ! v component
                      pgd2(ilonr,jj2,jk2) = factor * ( pgd2(ilonr,jj2,jk2) + (zfms(jj,1,jk2) + zfma(jj,1,jk2)) )
                      pgd2(iloni,jj2,jk2) = factor * ( pgd2(iloni,jj2,jk2) + (zfms(jj,2,jk2) + zfma(jj,2,jk2)) )
                   else
                      ! southern latitudes
                      ! u component
                      pgd2(ilonr,jj2,jk ) = factor * ( pgd2(ilonr,jj2,jk ) -(zfms(jj2,1,jk ) - zfma(jj2,1,jk )) )
                      pgd2(iloni,jj2,jk ) = factor * ( pgd2(iloni,jj2,jk ) -(zfms(jj2,2,jk ) - zfma(jj2,2,jk )) )
                      ! v component
                      pgd2(ilonr,jj2,jk2) = factor * ( pgd2(ilonr,jj2,jk2) +(zfms(jj2,1,jk2) - zfma(jj2,1,jk2)) )
                      pgd2(iloni,jj2,jk2) = factor * ( pgd2(iloni,jj2,jk2) +(zfms(jj2,2,jk2) - zfma(jj2,2,jk2)) )
                   endif
                else
                   ! Special case for the equator (jj.eq.jj2)
                   write(*,*) 'SPGDPAR: special case of jj.eq.jj2!!!'
                   ! u component
                   pgd2(ilonr,jj2,jk ) = factor * ( pgd2(ilonr,jj2,jk ) -(zfms(jj,1,jk ) + zfma(jj,1,jk )) )
                   pgd2(iloni,jj2,jk ) = factor * ( pgd2(iloni,jj2,jk ) -(zfms(jj,2,jk ) + zfma(jj,2,jk )) )
                   ! v component
                   pgd2(ilonr,jj2,jk2) = factor * ( pgd2(ilonr,jj2,jk2) +(zfms(jj,1,jk2) + zfma(jj,1,jk2)) )
                   pgd2(iloni,jj2,jk2) = factor * ( pgd2(iloni,jj2,jk2) +(zfms(jj,2,jk2) + zfma(jj,2,jk2)) )
                endif
             enddo ! jj
          endif
       enddo ! jk
    enddo  ! end loop on m
    !$OMP END PARALLEL DO

  end subroutine spgdpar


  subroutine gdsppar(psp,pgd2,nflev)
    implicit none

    ! Arguments:
    integer :: nflev
    real(8) :: psp(gst(gstID)%nla,2,gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: pgd2(2*gst(gstID)%maxmcount,gst(gstID)%nj, &
                      gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)

    ! Locals:
    integer :: jj, jj2, jk, jk2, ilonr, iloni, jm ,ila, inm, jn
    real(8) :: zjm, dlrwt(gst(gstID)%nj), dlrwocs(gst(gstID)%nj)
    real(8) :: dlalp(0:gst(gstID)%ntrunc,gst(gstID)%njlath)
    real(8) :: dldalp(0:gst(gstID)%ntrunc,gst(gstID)%njlath)
    real(8) :: dlsp(0:gst(gstID)%ntrunc,2, &
                    gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: dlsp2(0:gst(gstID)%ntrunc,2, &
                     gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: zfms(gst(gstID)%njlath+1,2, &
                    gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: zfma(gst(gstID)%njlath+1,2, &
                    gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)

    ! 1. Adjustment needed when an odd number of latitudes is considered
    dlrwt(:)   = gst(gstID)%rwt(:)
    !dlrwocs(:) = gst(gstID)%rwocs(:)
    do jj = 1, gst(gstID)%nj
      dlrwocs(jj) = gst(gstID)%rwt(jj)/(ec_ra*cos(gst(gstID)%rlati(jj)))
    enddo
    if (mod(gst(gstID)%nj,2).ne.0) then
       dlrwt(gst(gstID)%njlath)   = dlrwt(gst(gstID)%njlath)/2.d0
       dlrwocs(gst(gstID)%njlath) = dlrwocs(gst(gstID)%njlath)/2.d0
    end if

    !$OMP PARALLEL DO PRIVATE(DLALP,DLDALP,DLSP,DLSP2,ZFMS,ZFMA, &
    !$OMP INM,ILA,JM,JN,JK,JK2,JJ,JJ2,ZJM,ILONR,ILONI)
    do jm = gst(gstID)%mymBeg, gst(gstID)%mymEnd, gst(gstID)%mymSkip

       ilonr = 2*gst(gstID)%mymIndex(jm)-1
       iloni = 2*gst(gstID)%mymIndex(jm)
       zjm   = real(jm,8)

       ! 2.1 Fetch the Legendre functions and their derivatives for this choice of "m"
       call getalp(dlalp,dldalp,gst(gstID)%njlath,gst(gstID)%ntrunc,gst(gstID)%ntrunc,jm)

       do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
          do jj = 1, gst(gstID)%njlath
             zfms(jj,1,jk) = 0.0d0
             zfms(jj,2,jk) = 0.0d0
             zfma(jj,1,jk) = 0.0d0
             zfma(jj,2,jk) = 0.0d0
          enddo
       enddo

       ! 2.2  Build the symmetric and anti-symmetric Fourier coefficients including
       !      the appropriate quadrature weights (see scientific notes)
       do jj = 1, gst(gstID)%njlath
          jj2 = gst(gstID)%nj - jj + 1

          ! 2.2.1  Coefficients for scalar fields
          do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
             if(jk .gt. 2*nflev) then
                ! symmetric coefficients
                zfms(jj,1,jk) = dlrwt(jj)*(pgd2(ilonr,jj2,jk) + pgd2(ilonr,jj,jk))
                zfms(jj,2,jk) = dlrwt(jj)*(pgd2(iloni,jj2,jk) + pgd2(iloni,jj,jk))
                ! antisymmetric coefficients
                zfma(jj,1,jk) = dlrwt(jj)*(pgd2(ilonr,jj2,jk) - pgd2(ilonr,jj,jk))
                zfma(jj,2,jk) = dlrwt(jj)*(pgd2(iloni,jj2,jk) - pgd2(iloni,jj,jk))
             endif
          enddo

          ! 2.2.2 Coefficients associated with the wind fields
          do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd, 2
             if(jk .le. 2*nflev) then
                jk2 = jk + 1  ! jk is u, jk2 is v
                ! vorticity: symmetric coefficients
                zfms(jj,1,jk) = -zjm*dlrwocs(jj)*(pgd2(iloni,jj2,jk2) + pgd2(iloni,jj,jk2))
                zfms(jj,2,jk) =  zjm*dlrwocs(jj)*(pgd2(ilonr,jj2,jk2) + pgd2(ilonr,jj,jk2))
                ! vorticity: antisymmetric coefficients
                zfma(jj,1,jk) = -zjm*dlrwocs(jj)*(pgd2(iloni,jj2,jk2) - pgd2(iloni,jj,jk2))
                zfma(jj,2,jk) =  zjm*dlrwocs(jj)*(pgd2(ilonr,jj2,jk2) - pgd2(ilonr,jj,jk2))
                ! divergence: symmetric coefficients
                zfms(jj,1,jk2) = -zjm*dlrwocs(jj)*(pgd2(iloni,jj2,jk) + pgd2(iloni,jj,jk))
                zfms(jj,2,jk2) =  zjm*dlrwocs(jj)*(pgd2(ilonr,jj2,jk) + pgd2(ilonr,jj,jk))
                ! divergence: antisymmetric coefficients
                zfma(jj,1,jk2) = -zjm*dlrwocs(jj)*(pgd2(iloni,jj2,jk) - pgd2(iloni,jj,jk))
                zfma(jj,2,jk2) =  zjm*dlrwocs(jj)*(pgd2(ilonr,jj2,jk) - pgd2(ilonr,jj,jk))
             endif
          enddo
       enddo

       ! 2.3 First one with ALP for all scalar fields and for half the terms
       !     required to define the divergence and vorticity
       call legdir2d(jm,zfms,zfma,dlsp,dlalp,gst(gstID)%njlath,gst(gstID)%ntrunc,gst(gstID)%ntrunc)
  
       ! 2.4  Second transform with DALP to complete the construction of the
       !      vorticity and divergence fields
       do jj = 1, gst(gstID)%njlath
          jj2 = gst(gstID)%nj - jj + 1
          do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd, 2
             if(jk <= 2*nflev) then ! only for winds
                jk2 = jk + 1

                ! symmetric coefficients for zonal wind
                zfms(jj,1,jk) = dlrwocs(jj)*(pgd2(ilonr,jj2,jk) + pgd2(ilonr,jj,jk))
                zfms(jj,2,jk) = dlrwocs(jj)*(pgd2(iloni,jj2,jk) + pgd2(iloni,jj,jk))
                ! antisymmetric coefficients for zonal wind
                zfma(jj,1,jk) = dlrwocs(jj)*(pgd2(ilonr,jj2,jk) - pgd2(ilonr,jj,jk))
                zfma(jj,2,jk) = dlrwocs(jj)*(pgd2(iloni,jj2,jk) - pgd2(iloni,jj,jk))

                ! symmetric coefficients for zonal wind
                zfms(jj,1,jk2) = -dlrwocs(jj)*(pgd2(ilonr,jj2,jk2) + pgd2(ilonr,jj,jk2))
                zfms(jj,2,jk2) = -dlrwocs(jj)*(pgd2(iloni,jj2,jk2) + pgd2(iloni,jj,jk2))
                ! antisymmetric coefficients for zonal wind
                zfma(jj,1,jk2) = -dlrwocs(jj)*(pgd2(ilonr,jj2,jk2) - pgd2(ilonr,jj,jk2))
                zfma(jj,2,jk2) = -dlrwocs(jj)*(pgd2(iloni,jj2,jk2) - pgd2(iloni,jj,jk2))

             endif
          enddo
       enddo

       call legdir2d(jm,zfma,zfms,dlsp2,dldalp,gst(gstID)%njlath,gst(gstID)%ntrunc,gst(gstID)%ntrunc)

       ! 2.5  Transfer the result in the global state
       do jn = jm, gst(gstID)%ntrunc
          ila = gst(gstid)%nind(jm) + jn - jm
          inm = jn - jm
          do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
             if(jk <= 2*nflev) then ! for winds
                psp(ila,1,jk) = dlsp(inm,1,jk) + dlsp2(inm,1,jk)
                psp(ila,2,jk) = dlsp(inm,2,jk) + dlsp2(inm,2,jk)
             else ! for scalar fields
                psp(ila,1,jk) = dlsp(inm,1,jk)
                psp(ila,2,jk) = dlsp(inm,2,jk)
             endif
          enddo
       enddo

    enddo ! jm
    !$OMP END PARALLEL DO
  end subroutine gdsppar


  subroutine gst_spgda(psp,pgd,nflev)
    implicit none

    ! Arguments:
    integer :: nflev
    real(8) :: psp(gst(gstID)%myNla,2,gst(gstID)%nk)
    real(8) :: pgd(gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                   gst(gstID)%myLatBeg:gst(gstID)%myLatEnd, gst(gstID)%nk)

    ! Locals:
    real(8), allocatable :: psp2(:,:,:),pgd2(:,:,:),pgd3(:,:,:)

    ! check if this mpi task will deal with winds during Legendre transform
    if(gst(gstID)%myLevBeg.le.2*nflev) then
      ! ensure that the number of levels on this mpi task is even to allow interleaving of u and v
      ! only necessary when number of levels on an mpi task is less than all wind levels (2nd condition)
      if( (mod(gst(gstID)%myLevCount,2).ne.0) .and. (gst(gstID)%myLevCount < 2*nflev) ) then
        write(*,*) 'GST_SPGDA: myLevCount = ',gst(gstID)%myLevCount
        call utl_abort('GST_SPGDA: Number of levels on this mpi task must be even!')
      endif
    endif

    allocate(psp2(gst(gstID)%nla, 2, gst(gstID)%myLevBeg:gst(gstID)%myLevEnd))
    allocate(pgd2(2*gst(gstID)%maxmcount, gst(gstID)%nj,  gst(gstID)%myLevBeg:gst(gstID)%myLevEnd))
    allocate(pgd3(gst(gstID)%ni, gst(gstID)%myLatBeg:gst(gstID)%myLatEnd, gst(gstID)%myLevBeg:gst(gstID)%myLevEnd))

    call adjnorm(pgd)

    ! First reorder wind components to have u and v for same level on same mpi task
    call interleaveWinds_gd(pgd,nflev)

    ! Transpose data along npex from Longitudes to Levels
    call transpose2d_LontoLev(pgd,pgd3)

    ! Fourier transform all fields for all latitudes
    call utl_tmg_start(151,'low-level--gst_fft')
    call fft3dvar(pgd3,-1)
    call tmg_stop(151)

    ! Transpose data along npey from Latitudes to M
    call transpose2d_LattoM(pgd3,pgd2)
    deallocate(pgd3)

    ! Direct Legendre transform including wind transformations
    call utl_tmg_start(150,'low-level--gst_lt')
    call spgdapar(psp2,pgd2,nflev)
    call tmg_stop(150)
    deallocate(pgd2)

    ! Transpose data along npex from Levels to N
    call transpose2d_LevtoN(psp2,psp)
    deallocate(psp2)

    ! Now undo reordering of wind components 
    call unInterleaveWinds_sp(psp,nflev)

  end subroutine gst_spgda


  subroutine spgdapar(psp,pgd2,nflev)
    implicit none

    ! Arguments:
    integer :: nflev
    real(8) :: psp(gst(gstID)%nla,2,gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: pgd2(2*gst(gstID)%maxmCount,gst(gstID)%nj, &
                      gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)

    ! Locals:
    integer :: jj, jj2, jk, jk2, ilonr, iloni, jm ,ila, inm, jn
    real(8) :: zjm,factor,dlrwt(gst(gstID)%nj)
    real(8) :: dlalp(0:gst(gstID)%ntrunc,gst(gstID)%njlath)
    real(8) :: dldalp(0:gst(gstID)%ntrunc,gst(gstID)%njlath)
    real(8) :: dlsp(0:gst(gstID)%ntrunc,2, &
                    gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: dlsp2(0:gst(gstID)%ntrunc,2, &
                     gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: zfms(gst(gstID)%njlath+1,2, &
                    gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: zfma(gst(gstID)%njlath+1,2, &
                    gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)

    !    1. Set up according to the desired grid
    !       ---------------------
    dlrwt(:) = gst(gstID)%rwt(:)
    if (mod(gst(gstID)%nj,2).ne.0) then
       dlrwt(gst(gstID)%njlath) = dlrwt(gst(gstID)%njlath)/2.d0
    end if

    ! 2. Fourier transform all fields for all latitudes
    !$OMP PARALLEL DO PRIVATE(DLALP,DLDALP,DLSP,DLSP2,ZFMS,ZFMA, &
    !$OMP INM,ILA,JM,JN,JK,JK2,JJ,JJ2,ZJM,ILONR,ILONI,FACTOR)
    do jm = gst(gstID)%mymBeg, gst(gstID)%mymEnd, gst(gstID)%mymSkip

          ilonr = 2*gst(gstID)%mymIndex(jm)-1
          iloni = 2*gst(gstID)%mymIndex(jm)
          zjm   = real(jm,8)

          ! 2.1 Fetch the Legendre functions and their derivatives for this choice of "m"
          call getalp(dlalp,dldalp,gst(gstID)%njlath,gst(gstID)%ntrunc,gst(gstID)%ntrunc,jm)

          ! 2.2  Build the symmetric and anti-symmetric Fourier coefficients including
          !      the appropriate quadrature weights (see scientific notes)
          do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
             do jj = 1, gst(gstID)%njlath
                zfms(jj,1,jk) = 0.0d0
                zfms(jj,2,jk) = 0.0d0
                zfma(jj,1,jk) = 0.0d0
                zfma(jj,2,jk) = 0.0d0
             enddo
          enddo

          ! 2.2.1 Coefficients associated with the scalar fields
          do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
             if(jk .gt. 2*nflev) then
                do jj = 1, gst(gstID)%nj
                   jj2 = gst(gstID)%nj-jj+1
                   if(jj.le.gst(gstID)%njlath) then
                      ! Northern hemisphere
                      zfms(jj,1,jk) = pgd2(ilonr,jj2,jk)
                      zfms(jj,2,jk) = pgd2(iloni,jj2,jk)
                      zfma(jj,1,jk) = pgd2(ilonr,jj2,jk)
                      zfma(jj,2,jk) = pgd2(iloni,jj2,jk)
                   else
                      ! Southern hemisphere
                      zfms(jj2,1,jk) = zfms(jj2,1,jk) + pgd2(ilonr,jj2,jk)
                      zfms(jj2,2,jk) = zfms(jj2,2,jk) + pgd2(iloni,jj2,jk)
                      zfma(jj2,1,jk) = zfma(jj2,1,jk) - pgd2(ilonr,jj2,jk)
                      zfma(jj2,2,jk) = zfma(jj2,2,jk) - pgd2(iloni,jj2,jk)
                   endif
                enddo
             endif
          enddo

          ! 2.2.2 Coefficients associated with the wind fields
          do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd, 2
             if(jk .le. 2*nflev) then
                jk2 = jk + 1
                do jj = 1, gst(gstID)%nj
                   jj2 = gst(gstID)%nj-jj+1
                   if(jj.le.gst(gstID)%njlath) then
                      ! Northern hemisphere
                      ! vorticity: symmetric coefficients
                      zfms(jj,1,jk ) = -pgd2(iloni,jj2,jk2)
                      zfms(jj,2,jk ) =  pgd2(ilonr,jj2,jk2)
                      ! vorticity: antisymmetric coefficients
                      zfma(jj,1,jk ) = -pgd2(iloni,jj2,jk2)
                      zfma(jj,2,jk ) =  pgd2(ilonr,jj2,jk2)
                      ! divergence: symmetric coefficients
                      zfms(jj,1,jk2) = -pgd2(iloni,jj2,jk )
                      zfms(jj,2,jk2) =  pgd2(ilonr,jj2,jk )
                      ! divergence: antisymmetric coefficients
                      zfma(jj,1,jk2) = -pgd2(iloni,jj2,jk )
                      zfma(jj,2,jk2) =  pgd2(ilonr,jj2,jk )
                   else
                      ! Southern hemisphere
                      ! vorticity: symmetric coefficients
                      zfms(jj2,1,jk ) = zfms(jj2,1,jk ) - pgd2(iloni,jj2,jk2)
                      zfms(jj2,2,jk ) = zfms(jj2,2,jk ) + pgd2(ilonr,jj2,jk2)
                      ! vorticity: antisymmetric coefficients
                      zfma(jj2,1,jk ) = zfma(jj2,1,jk ) + pgd2(iloni,jj2,jk2)
                      zfma(jj2,2,jk ) = zfma(jj2,2,jk ) - pgd2(ilonr,jj2,jk2)
                      ! divergence: symmetric coefficients
                      zfms(jj2,1,jk2) = zfms(jj2,1,jk2) - pgd2(iloni,jj2,jk )
                      zfms(jj2,2,jk2) = zfms(jj2,2,jk2) + pgd2(ilonr,jj2,jk )
                      ! divergence: antisymmetric coefficients
                      zfma(jj2,1,jk2) = zfma(jj2,1,jk2) + pgd2(iloni,jj2,jk )
                      zfma(jj2,2,jk2) = zfma(jj2,2,jk2) - pgd2(ilonr,jj2,jk )
                   endif
                enddo ! jj
             endif
          enddo ! jk

          do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
             do jj = 1, gst(gstID)%njlath
                zfms(jj,1,jk) = dlrwt(jj)*zfms(jj,1,jk)
                zfms(jj,2,jk) = dlrwt(jj)*zfms(jj,2,jk)
                zfma(jj,1,jk) = dlrwt(jj)*zfma(jj,1,jk)
                zfma(jj,2,jk) = dlrwt(jj)*zfma(jj,2,jk)
             enddo
             if(jk .le. 2*nflev) then
                do jj = 1, gst(gstID)%njlath
                   factor = ec_ra / cos(gst(gstID)%rlati(jj))
                   zfms(jj,1,jk) = factor*zjm*zfms(jj,1,jk)
                   zfms(jj,2,jk) = factor*zjm*zfms(jj,2,jk)
                   zfma(jj,1,jk) = factor*zjm*zfma(jj,1,jk)
                   zfma(jj,2,jk) = factor*zjm*zfma(jj,2,jk)
                enddo
             endif
          enddo ! jk

          ! 2.3 First one with ALP for all scalar fields and for half the terms
          !     required to define the divergence and vorticity
          call legdir2d(jm,zfms,zfma,dlsp,dlalp,gst(gstID)%njlath,gst(gstID)%ntrunc,gst(gstID)%ntrunc)
                                
          ! 2.4  Second transform with DALP to complete the construction of the
          !      vorticity and divergence fields
          do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
             if(jk .le. 2*nflev) then
                do jj = 1, gst(gstID)%njlath
                   zfms(jj,1,jk) = 0.0d0
                   zfms(jj,2,jk) = 0.0d0
                   zfma(jj,1,jk) = 0.0d0
                   zfma(jj,2,jk) = 0.0d0
                enddo
             endif
          enddo

          do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd, 2
             if(jk .le. 2*nflev) then
                jk2 = jk + 1
                do jj = 1, gst(gstID)%nj
                   jj2 = gst(gstID)%nj-jj+1
                   if(jj.le.gst(gstID)%njlath) then
                      ! Northern hemisphere
                      ! symmetric coefficients for zonal wind
                      zfms(jj,1,jk ) =  pgd2(ilonr,jj2,jk )
                      zfms(jj,2,jk ) =  pgd2(iloni,jj2,jk )
                      ! antisymmetric coefficients for zonal wind
                      zfma(jj,1,jk ) =  pgd2(ilonr,jj2,jk )
                      zfma(jj,2,jk ) =  pgd2(iloni,jj2,jk )
                      ! symmetric coefficients for meridional wind
                      zfms(jj,1,jk2) = -pgd2(ilonr,jj2,jk2)
                      zfms(jj,2,jk2) = -pgd2(iloni,jj2,jk2)
                      ! antisymmetric coefficients for meridional wind
                      zfma(jj,1,jk2) = -pgd2(ilonr,jj2,jk2)
                      zfma(jj,2,jk2) = -pgd2(iloni,jj2,jk2)
                   else
                      ! Southern hemisphere
                      ! symmetric coefficients for zonal wind
                      zfms(jj2,1,jk ) = zfms(jj2,1,jk ) + pgd2(ilonr,jj2,jk )
                      zfms(jj2,2,jk ) = zfms(jj2,2,jk ) + pgd2(iloni,jj2,jk )
                      ! antisymmetric coefficients for zonal wind
                      zfma(jj2,1,jk ) = zfma(jj2,1,jk ) - pgd2(ilonr,jj2,jk )
                      zfma(jj2,2,jk ) = zfma(jj2,2,jk ) - pgd2(iloni,jj2,jk )
                      ! symmetric coefficients for meridional wind
                      zfms(jj2,1,jk2) = zfms(jj2,1,jk2) - pgd2(ilonr,jj2,jk2)
                      zfms(jj2,2,jk2) = zfms(jj2,2,jk2) - pgd2(iloni,jj2,jk2)
                      ! antisymmetric coefficients for meridional wind
                      zfma(jj2,1,jk2) = zfma(jj2,1,jk2) + pgd2(ilonr,jj2,jk2)
                      zfma(jj2,2,jk2) = zfma(jj2,2,jk2) + pgd2(iloni,jj2,jk2)
                   endif
                enddo
             endif
          enddo

          do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
             if(jk .le. 2*nflev) then
                do jj = 1, gst(gstID)%njlath
                   factor = ec_ra / cos(gst(gstID)%rlati(jj))
                   zfms(jj,1,jk) = factor*dlrwt(jj)*zfms(jj,1,jk)
                   zfms(jj,2,jk) = factor*dlrwt(jj)*zfms(jj,2,jk)
                   zfma(jj,1,jk) = factor*dlrwt(jj)*zfma(jj,1,jk)
                   zfma(jj,2,jk) = factor*dlrwt(jj)*zfma(jj,2,jk)
                enddo
             endif
          enddo

          call legdir2d(jm,zfma,zfms,dlsp2,dldalp,gst(gstID)%njlath,gst(gstID)%ntrunc,gst(gstID)%ntrunc)

          ! 2.5  Transfer the result in the global state
          do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
             do jn = jm, gst(gstID)%ntrunc
                ila = gst(gstID)%nind(jm) + jn - jm
                inm = jn - jm
                if(jk .le. 2*nflev) then
                   psp(ila,1,jk) = -gst(gstID)%r1snp1(ila)*(dlsp(inm,1,jk) + dlsp2(inm,1,jk))
                   psp(ila,2,jk) = -gst(gstID)%r1snp1(ila)*(dlsp(inm,2,jk) + dlsp2(inm,2,jk))
                else
                   psp(ila,1,jk) = dlsp(inm,1,jk)
                   psp(ila,2,jk) = dlsp(inm,2,jk)
                endif
             enddo
          enddo

    ! End of loop on zonal wavenumbers
    enddo
    !$OMP END PARALLEL DO

  end subroutine spgdapar


  subroutine gst_speree(psp,pgd)
    implicit none

    ! Arguments:
    real(8) :: psp(gst(gstID)%myNla, 2, gst(gstID)%nk)
    real(8) :: pgd(gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                   gst(gstID)%myLatBeg:gst(gstID)%myLatEnd, gst(gstID)%nk)

    ! Locals:
    real(8), allocatable :: psp2(:,:,:),pgd2(:,:,:),pgd3(:,:,:)
    integer :: jlat, jk, jlon

    allocate(psp2(gst(gstID)%nla, 2, gst(gstID)%myLevBeg:gst(gstID)%myLevEnd))
    allocate(pgd2(2*gst(gstID)%maxmcount, gst(gstID)%nj, gst(gstID)%myLevBeg:gst(gstID)%myLevEnd))
    allocate(pgd3(gst(gstID)%ni, gst(gstID)%myLatBeg:gst(gstID)%myLatEnd, gst(gstID)%myLevBeg:gst(gstID)%myLevEnd))

    ! 1.0 Transpose data along npex from N to Levels
    call transpose2d_NtoLev(psp,psp2)

    ! 1.1 Inverse Legendre transform (lon -> m)
    call utl_tmg_start(150,'low-level--gst_lt')
    call spereepar(psp2,pgd2)
    call tmg_stop(150)
    deallocate(psp2)

    ! 1.2 Transpose data along npey from M to Latitudes
    call transpose2d_MtoLat(pgd2,pgd3)
    deallocate(pgd2)

    ! 2.1 Reset to zero the modes that are not part of the truncation
    !$OMP PARALLEL DO PRIVATE (JLAT,JLON,JK)
    do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
      do jlat = gst(gstID)%myLatBeg, gst(gstID)%myLatEnd
        do jlon = 2*(gst(gstID)%ntrunc+1)+1, gst(gstID)%ni
          pgd3(jlon,jlat,jk) = 0.d0
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    ! 2.2 Apply the inverse FFT 
    call utl_tmg_start(151,'low-level--gst_fft')
    call fft3dvar(pgd3,+1)
    call tmg_stop(151)

    ! 2.3 Transpose data along npex from Levels to Longitudes
    call transpose2d_LevtoLon(pgd3,pgd)
    deallocate(pgd3)

  end subroutine gst_speree


  subroutine gst_speree_kij(psp,pgd)
    implicit none

    ! Arguments:
    real(8) :: psp(gst(gstID)%myNla, 2, gst(gstID)%nk)
    real(8) :: pgd(gst(gstID)%nk, gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                   gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)

    ! Locals:
    real(8), allocatable :: psp2(:,:,:),pgd2(:,:,:),pgd3(:,:,:)
    integer :: jlat, jk, jlon, ierr

    call utl_tmg_start(152,'low-level--gst_barr')
    if(mpi_doBarrier) call rpn_comm_barrier('GRID',ierr)
    call tmg_stop(152)

    allocate(psp2(gst(gstID)%nla, 2, gst(gstID)%myLevBeg:gst(gstID)%myLevEnd))
    allocate(pgd2(gst(gstID)%maxMyLevCount, 2*gst(gstID)%maxmcount, gst(gstID)%nj))
    allocate(pgd3(gst(gstID)%maxMyLevCount, gst(gstID)%ni, gst(gstID)%myLatBeg:gst(gstID)%myLatEnd))

    pgd2( (gst(gstID)%myLevCount+1):gst(gstID)%maxMyLevCount, :, :) = 0.0d0

    ! 1.0 Transpose data along npex from N to Levels
    call transpose2d_NtoLev(psp,psp2)

    ! 1.1 Inverse Legendre transform (lon -> m)
    call utl_tmg_start(150,'low-level--gst_lt')
    call spereepar_kij(psp2,pgd2)
    call tmg_stop(150)
    deallocate(psp2)

    ! 1.2 Transpose data along npey from M to Latitudes
    call transpose2d_MtoLat_kij(pgd2,pgd3)
    deallocate(pgd2)

    ! 2.1 Reset to zero the modes that are not part of the truncation
    !$OMP PARALLEL DO PRIVATE (JLAT,JLON,JK)
    do jlat = gst(gstID)%myLatBeg, gst(gstID)%myLatEnd
      do jlon = 2*(gst(gstID)%ntrunc+1)+1, gst(gstID)%ni
        do jk = 1, gst(gstID)%maxMyLevCount
          pgd3(jk,jlon,jlat) = 0.d0
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    ! 2.2 Apply the inverse FFT 
    call utl_tmg_start(151,'low-level--gst_fft')
    call fft3dvar_kij(pgd3,+1)
    call tmg_stop(151)

    ! 2.3 Transpose data along npex from Levels to Longitudes
    if( gst(gstID)%lonLatDivisible ) then
      if (pre_specTransReal == 4) then
        call transpose2d_LevtoLon_kij_mpitypes4(pgd3,pgd)
      else
        call transpose2d_LevtoLon_kij_mpitypes8(pgd3,pgd)
      end if
    else
      call transpose2d_LevtoLon_kij(pgd3,pgd)
    end if
    deallocate(pgd3)

    call utl_tmg_start(152,'low-level--gst_barr')
    if(mpi_doBarrier) call rpn_comm_barrier('GRID',ierr)
    call tmg_stop(152)

  end subroutine gst_speree_kij

  subroutine gst_speree_ad(psp,pgd)
    implicit none

    ! Arguments:
    real(8) :: psp(gst(gstID)%myNla, 2, gst(gstID)%nk)
    real(8) :: pgd(gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                   gst(gstID)%myLatBeg:gst(gstID)%myLatEnd, gst(gstID)%nk)

    call adjnorm(pgd)

    call gst_reespe(psp,pgd)

  end subroutine gst_speree_ad


  subroutine gst_reespe(psp,pgd)
    implicit none

    ! Arguments:
    real(8) :: psp(gst(gstID)%myNla, 2, gst(gstID)%nk)
    real(8) :: pgd(gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                   gst(gstID)%myLatBeg:gst(gstID)%myLatEnd, gst(gstID)%nk)

    ! Locals:
    real(8), allocatable :: psp2(:,:,:),pgd2(:,:,:),pgd3(:,:,:)

    allocate(psp2(gst(gstID)%nla, 2, gst(gstID)%myLevBeg:gst(gstID)%myLevEnd))
    allocate(pgd2(2*gst(gstID)%maxmcount, gst(gstID)%nj, gst(gstID)%myLevBeg:gst(gstID)%myLevEnd))
    allocate(pgd3(gst(gstID)%ni, gst(gstID)%myLatBeg:gst(gstID)%myLatEnd, gst(gstID)%myLevBeg:gst(gstID)%myLevEnd))

    ! Transpose data along npex from Longitudes to Levels
    call transpose2d_LontoLev(pgd,pgd3)

    ! 1. Apply the FFT
    call utl_tmg_start(151,'low-level--gst_fft')
    call fft3dvar(pgd3,-1)
    call tmg_stop(151)

    ! Transpose data along npey from Latitudes to M
    call transpose2d_LattoM(pgd3,pgd2)
    deallocate(pgd3)

    ! 2. Direct Legendre transform
    call utl_tmg_start(150,'low-level--gst_lt')
    call reespepar(pgd2,psp2)
    call tmg_stop(150)
    deallocate(pgd2)

    ! Transpose data along npex from Levels to N
    call transpose2d_LevtoN(psp2,psp)
    deallocate(psp2)

  end subroutine gst_reespe


  subroutine gst_reespe_kij(psp,pgd)
    implicit none

    ! Arguments:
    real(8) :: psp(gst(gstID)%myNla, 2, gst(gstID)%nk)
    real(8) :: pgd(gst(gstID)%nk, gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                   gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)

    ! Locals:
    real(8), allocatable :: psp2(:,:,:),pgd2(:,:,:),pgd3(:,:,:)
    integer :: ierr

    call utl_tmg_start(152,'low-level--gst_barr')
    if(mpi_doBarrier) call rpn_comm_barrier('GRID',ierr)
    call tmg_stop(152)

    allocate(psp2(gst(gstID)%nla, 2, gst(gstID)%myLevBeg:gst(gstID)%myLevEnd))
    allocate(pgd2(gst(gstID)%maxMyLevCount, 2*gst(gstID)%maxmcount, gst(gstID)%nj))
    allocate(pgd3(gst(gstID)%maxMyLevCount, gst(gstID)%ni, gst(gstID)%myLatBeg:gst(gstID)%myLatEnd))

    ! Transpose data along npex from Longitudes to Levels
    if( gst(gstID)%lonLatDivisible ) then
      if (pre_specTransReal == 4) then
        call transpose2d_LontoLev_kij_mpitypes4(pgd,pgd3)
      else
        call transpose2d_LontoLev_kij_mpitypes8(pgd,pgd3)
      end if
    else
      call transpose2d_LontoLev_kij(pgd,pgd3)
    end if

    ! 1. Apply the FFT
    call utl_tmg_start(151,'low-level--gst_fft')
    call fft3dvar_kij(pgd3,-1)
    call tmg_stop(151)

    ! Transpose data along npey from Latitudes to M
    call transpose2d_LattoM_kij(pgd3,pgd2)
    deallocate(pgd3)

    ! 2. Direct Legendre transform
    call utl_tmg_start(150,'low-level--gst_lt')
    call reespepar_kij(pgd2,psp2)
    call tmg_stop(150)
    deallocate(pgd2)

    ! Transpose data along npex from Levels to N
    call transpose2d_LevtoN(psp2,psp)
    deallocate(psp2)

    call utl_tmg_start(152,'low-level--gst_barr')
    if(mpi_doBarrier) call rpn_comm_barrier('GRID',ierr)
    call tmg_stop(152)

  end subroutine gst_reespe_kij

  subroutine gst_speree_kij_ad(psp,pgd)
    implicit none

    ! Arguments:
    real(8) :: psp(gst(gstID)%myNla, 2, gst(gstID)%nk)
    real(8) :: pgd(gst(gstID)%nk, gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                   gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)

    call adjnorm_kij(pgd)

    call gst_reespe_kij(psp,pgd)

  end subroutine gst_speree_kij_ad

  subroutine adjnorm(pgd)
    implicit none

    ! Arguments:
    real(8) :: pgd(gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                   gst(gstID)%myLatBeg:gst(gstID)%myLatEnd, gst(gstID)%nk)

    ! Locals:
    integer :: jk, jlon, jlat
    integer :: lat1, lat2, lon1, lon2
    real(8) :: rwtinv(gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)

    lat1 = gst(gstID)%myLatBeg
    lat2 = gst(gstID)%myLatEnd
    lon1 = gst(gstID)%myLonBeg
    lon2 = gst(gstID)%myLonEnd

    do jlat = lat1, lat2
      rwtinv(jlat) = real(gst(gstID)%ni,8) / gst_getRWT(jlat, gstID)
    enddo

    !$OMP PARALLEL DO PRIVATE (jk,jlat,jlon)
    do jk = 1, gst(gstID)%nk
        do jlat = lat1, lat2
          do jlon = lon1, lon2
            pgd(jlon,jlat,jk) = rwtinv(jlat) * pgd(jlon,jlat,jk)
          end do
        end do
    end do
    !$OMP END PARALLEL DO

  end subroutine adjnorm

  subroutine adjnorm_kij(pgd)
    implicit none

    ! Arguments:
    real(8) :: pgd(gst(gstID)%nk, gst(gstID)%myLonBeg:gst(gstID)%myLonEnd, &
                   gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)

    ! Locals:
    integer :: jk, jlon, jlat
    integer :: lat1, lat2, lon1, lon2
    real(8) :: rwtinv(gst(gstID)%myLatBeg:gst(gstID)%myLatEnd)

    lat1 = gst(gstID)%myLatBeg
    lat2 = gst(gstID)%myLatEnd
    lon1 = gst(gstID)%myLonBeg
    lon2 = gst(gstID)%myLonEnd

    do jlat = lat1, lat2
      rwtinv(jlat) = real(gst(gstID)%ni,8) / gst_getRWT(jlat, gstID)
    enddo

    !$OMP PARALLEL DO PRIVATE (jlat,jlon,jk)
    do jlat = lat1, lat2
       do jlon = lon1, lon2
         do jk = 1, gst(gstID)%nk
            pgd(jk,jlon,jlat) = rwtinv(jlat) * pgd(jk,jlon,jlat)
         end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine adjnorm_kij

  subroutine spereepar(psp,pgd2)
    !
    !:Purpose: Inverse spectral transform(MPI PARALLEL LOOP)
    implicit none

    ! Arguments:
    real(8) :: psp(gst(gstID)%nla,2,gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: pgd2(2*gst(gstID)%maxmcount,gst(gstID)%nj, &
                      gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)

    ! Locals
    integer :: jj, jj2, jm, jn, ilonr, iloni, jk, ila, inm

    real(8) :: dlalp(0:gst(gstID)%ntrunc,gst(gstID)%njlath)
    real(8) :: dldalp(0:gst(gstID)%ntrunc,gst(gstID)%njlath)
    real(8) :: zfms(gst(gstID)%njlath+1,2, &
                        gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: zfma(gst(gstID)%njlath+1,2, &
                        gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: dlsp(0:gst(gstID)%ntrunc,2, &
                      gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)

    ! Inverse Legendre transform

    !$OMP PARALLEL DO PRIVATE(DLALP,DLDALP,DLSP,ZFMS,ZFMA, &
    !$OMP INM,ILA,JM,JN,JK,JJ,JJ2,ILONR,ILONI)
    do jm = gst(gstID)%mymBeg, gst(gstID)%mymEnd, gst(gstID)%mymSkip

       ilonr = 2*gst(gstID)%mymIndex(jm)-1
       iloni = 2*gst(gstID)%mymIndex(jm)

       ! 2.1 Copy global spectral state into local spectral state
       do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
          do jn = jm, gst(gstID)%ntrunc
             ila = gst(gstID)%nind(jm) + jn - jm
             inm = jn - jm
             dlsp(inm,1,jk) = psp(ila,1,jk)
             dlsp(inm,2,jk) = psp(ila,2,jk)
          enddo
       enddo

       ! 2.2  Get Legendre polynomial (and its derivative) for all latitudes
       !      but for the chosen value of "m" from the global array
       call getalp(dlalp,dldalp,gst(gstID)%njlath,gst(gstID)%ntrunc,gst(gstID)%ntrunc,jm)

       ! 2.3  Perform the inverse Legendre transform for all fields
       call leginv2d(jm,zfms,zfma,dlsp,dlalp,gst(gstID)%njlath,gst(gstID)%ntrunc,gst(gstID)%ntrunc)

       ! 2.4 Passage to Fourier space

       do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
          do jj = 1, gst(gstID)%nj
             jj2 = gst(gstID)%nj - jj + 1
             if(jj.le.gst(gstID)%njlath) then
                pgd2(ilonr,jj2,jk) = zfms(jj,1,jk) + zfma(jj,1,jk)
                pgd2(iloni,jj2,jk) = zfms(jj,2,jk) + zfma(jj,2,jk)
             else
                pgd2(ilonr,jj2,jk) = zfms(jj2,1,jk) - zfma(jj2,1,jk)
                pgd2(iloni,jj2,jk) = zfms(jj2,2,jk) - zfma(jj2,2,jk)
             endif
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine spereepar


  subroutine spereepar_kij(psp,pgd2)
    !
    !:Purpose:  Inverse spectral transform(MPI PARALLEL LOOP)
    implicit none

    ! Arguments:
    real(8) :: psp(gst(gstID)%nla,2,gst(gstID)%myLevCount)
    real(8) :: pgd2(gst(gstID)%maxMyLevCount,2*gst(gstID)%maxmcount, &
                    gst(gstID)%nj)

    ! Locals
    integer :: jj, jj2, jm, jn, ilonr, iloni, jk, ila, inm

    real(8) :: dlalp(0:gst(gstID)%ntrunc,gst(gstID)%njlath)
    real(8) :: dldalp(0:gst(gstID)%ntrunc,gst(gstID)%njlath)
    real(8) :: zfms(gst(gstID)%myLevCount,gst(gstID)%njlath+1,2)
    real(8) :: zfma(gst(gstID)%myLevCount,gst(gstID)%njlath+1,2)
    real(8) :: dlsp(gst(gstID)%myLevCount,0:gst(gstID)%ntrunc,2)

    ! Inverse Legendre transform

    !$OMP PARALLEL DO PRIVATE(DLALP,DLDALP,DLSP,ZFMS,ZFMA, &
    !$OMP INM,ILA,JM,JN,JK,JJ,JJ2,ILONR,ILONI)
    do jm = gst(gstID)%mymBeg, gst(gstID)%mymEnd, gst(gstID)%mymSkip

       ilonr = 2*gst(gstID)%mymIndex(jm)-1
       iloni = 2*gst(gstID)%mymIndex(jm)

       ! 2.1 Copy global spectral state into local spectral state
       do jn = jm, gst(gstID)%ntrunc
          ila = gst(gstID)%nind(jm) + jn - jm
          inm = jn - jm
          do jk = 1, gst(gstID)%myLevCount
             dlsp(jk,inm,1) = psp(ila,1,jk)
             dlsp(jk,inm,2) = psp(ila,2,jk)
          enddo
       enddo

       ! 2.2  Get Legendre polynomial (and its derivative) for all latitudes
       !      but for the chosen value of "m" from the global array
       call getalp(dlalp,dldalp,gst(gstID)%njlath,gst(gstID)%ntrunc,gst(gstID)%ntrunc,jm)

       ! 2.3  Perform the inverse Legendre transform for all fields
       call leginv2d_kij(jm,zfms,zfma,dlsp,dlalp,gst(gstID)%njlath,gst(gstID)%ntrunc,gst(gstID)%ntrunc)

       ! 2.4 Passage to Fourier space

       do jj = 1, gst(gstID)%nj
          jj2 = gst(gstID)%nj - jj + 1
          if(jj.le.gst(gstID)%njlath) then
             do jk = 1, gst(gstID)%myLevCount
                pgd2(jk,ilonr,jj2) = zfms(jk,jj,1) + zfma(jk,jj,1)
                pgd2(jk,iloni,jj2) = zfms(jk,jj,2) + zfma(jk,jj,2)
             enddo
          else
             do jk = 1, gst(gstID)%myLevCount
                pgd2(jk,ilonr,jj2) = zfms(jk,jj2,1) - zfma(jk,jj2,1)
                pgd2(jk,iloni,jj2) = zfms(jk,jj2,2) - zfma(jk,jj2,2)
             enddo
          endif
       enddo

    enddo
    !$OMP END PARALLEL DO

  end subroutine spereepar_kij


  subroutine reespepar(pgd2,psp)
    implicit none

    ! Arguments:
    real(8) :: psp(gst(gstID)%nla,2,gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: pgd2(2*gst(gstID)%maxmcount,gst(gstID)%nj, &
                      gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)

    ! Locals:
    integer :: jj,jj2,jk,ilonr, iloni
    integer :: jm, ila, inm, jn

    real(8) :: dlalp(0:gst(gstID)%ntrunc,gst(gstID)%njlath)
    real(8) :: dldalp(0:gst(gstID)%ntrunc, gst(gstID)%njlath)
    real(8) :: dlsp(0:gst(gstID)%ntrunc,2, &
                    gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: zfms(gst(gstID)%njlath+1,2, &
                    gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: zfma(gst(gstID)%njlath+1,2, &
                    gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: dlrwt(gst(gstID)%nj)

    ! 1. Adjustment needed when an odd number of latitudes is considered
    dlrwt(:) = gst(gstID)%rwt(:)
    if (mod(gst(gstID)%nj,2).ne.0) then
       dlrwt(gst(gstID)%njlath) = dlrwt(gst(gstID)%njlath)/2.d0
    end if

    !$OMP PARALLEL DO PRIVATE(DLALP,DLDALP,DLSP,ZFMS,ZFMA, &
    !$OMP INM,ILA,JM,JN,JK,JJ,JJ2,ILONR,ILONI)
    do jm = gst(gstID)%mymBeg, gst(gstID)%mymEnd, gst(gstID)%mymSkip

      ilonr = 2*gst(gstID)%mymIndex(jm)-1
      iloni = 2*gst(gstID)%mymIndex(jm)

      ! 2.1 Fetch the Legendre functions and their derivatives for this choice of "m"
      call getalp(dlalp,dldalp,gst(gstID)%njlath,gst(gstID)%ntrunc,gst(gstID)%ntrunc,jm)

      do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
          ! 2.2  Build the symmetric and anti-symmetric Fourier coefficients including
          !      the appropriate quadrature weights (see scientific notes)
          do jj = 1, gst(gstID)%njlath
             zfms(jj,1,jk) = 0.0d0
             zfms(jj,2,jk) = 0.0d0
             zfma(jj,1,jk) = 0.0d0
             zfma(jj,2,jk) = 0.0d0
          enddo

          do jj = 1, gst(gstID)%nj
             jj2 = gst(gstID)%nj-jj+1
             if(jj.le.gst(gstID)%njlath) then
                ! Northern hemisphere
                zfms(jj,1,jk) = pgd2(ilonr,jj2,jk)
                zfms(jj,2,jk) = pgd2(iloni,jj2,jk)
                zfma(jj,1,jk) = pgd2(ilonr,jj2,jk)
                zfma(jj,2,jk) = pgd2(iloni,jj2,jk)
             else
                ! Southern hemisphere
                zfms(jj2,1,jk) = zfms(jj2,1,jk) + pgd2(ilonr,jj2,jk)
                zfms(jj2,2,jk) = zfms(jj2,2,jk) + pgd2(iloni,jj2,jk)
                zfma(jj2,1,jk) = zfma(jj2,1,jk) - pgd2(ilonr,jj2,jk)
                zfma(jj2,2,jk) = zfma(jj2,2,jk) - pgd2(iloni,jj2,jk)
             endif
          enddo

          do jj = 1,gst(gstID)%njlath
             zfms(jj,1,jk) = dlrwt(jj)*zfms(jj,1,jk)
             zfms(jj,2,jk) = dlrwt(jj)*zfms(jj,2,jk)
             zfma(jj,1,jk) = dlrwt(jj)*zfma(jj,1,jk)
             zfma(jj,2,jk) = dlrwt(jj)*zfma(jj,2,jk)
          enddo

      enddo ! jk

      ! 2.3 First one with ALP for all scalar fields and for half the terms
      !     required to define the divergence and vorticity
      call legdir2d(jm,zfms,zfma,dlsp,dlalp,gst(gstID)%njlath,gst(gstID)%ntrunc,gst(gstID)%ntrunc)

      ! 2.4 Transfer the result in the global state
      do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
          do jn = jm, gst(gstID)%ntrunc
            ila = gst(gstID)%nind(jm) + jn - jm
            inm = jn - jm
            psp(ila,1,jk) = dlsp(inm,1,jk)
            psp(ila,2,jk) = dlsp(inm,2,jk)
          enddo
      enddo ! jk

    ! End of loop on zonal wavenumbers
    enddo
    !$OMP END PARALLEL DO

  end subroutine reespepar


  subroutine reespepar_kij(pgd2,psp)
    implicit none

    ! Arguments:
    real(8) :: psp(gst(gstID)%nla,2,gst(gstID)%myLevCount)
    real(8) :: pgd2(gst(gstID)%maxMyLevCount,2*gst(gstID)%maxmcount, &
                    gst(gstID)%nj)

    ! Locals:
    integer :: jj,jj2,jk,ilonr, iloni
    integer :: jm, ila, inm, jn

    real(8) :: dlalp(0:gst(gstID)%ntrunc,gst(gstID)%njlath)
    real(8) :: dldalp(0:gst(gstID)%ntrunc, gst(gstID)%njlath)
    real(8) :: dlsp(gst(gstID)%myLevCount,0:gst(gstID)%ntrunc,2)
    real(8) :: zfms(gst(gstID)%myLevCount,gst(gstID)%njlath+1,2)
    real(8) :: zfma(gst(gstID)%myLevCount,gst(gstID)%njlath+1,2)
    real(8) :: dlrwt(gst(gstID)%nj)

    ! 1. Adjustment needed when an odd number of latitudes is considered
    dlrwt(:) = gst(gstID)%rwt(:)
    if (mod(gst(gstID)%nj,2).ne.0) then
       dlrwt(gst(gstID)%njlath) = dlrwt(gst(gstID)%njlath)/2.d0
    end if

    !$OMP PARALLEL DO PRIVATE(DLALP,DLDALP,DLSP,ZFMS,ZFMA, &
    !$OMP INM,ILA,JM,JN,JK,JJ,JJ2,ILONR,ILONI)
    do jm = gst(gstID)%mymBeg, gst(gstID)%mymEnd, gst(gstID)%mymSkip

      ilonr = 2*gst(gstID)%mymIndex(jm)-1
      iloni = 2*gst(gstID)%mymIndex(jm)

      ! 2.1 Fetch the Legendre functions and their derivatives for this choice of "m"
      call getalp(dlalp,dldalp,gst(gstID)%njlath,gst(gstID)%ntrunc,gst(gstID)%ntrunc,jm)

      ! 2.2  Build the symmetric and anti-symmetric Fourier coefficients including
      !      the appropriate quadrature weights (see scientific notes)
      do jj = 1, gst(gstID)%njlath
        do jk = 1, gst(gstID)%myLevCount
          zfms(jk,jj,1) = 0.0d0
          zfms(jk,jj,2) = 0.0d0
          zfma(jk,jj,1) = 0.0d0
          zfma(jk,jj,2) = 0.0d0
        enddo
      enddo

      do jj = 1, gst(gstID)%nj
        jj2 = gst(gstID)%nj-jj+1
        if(jj.le.gst(gstID)%njlath) then
          ! Northern hemisphere
          do jk = 1, gst(gstID)%myLevCount
            zfms(jk,jj,1) = pgd2(jk,ilonr,jj2)
            zfms(jk,jj,2) = pgd2(jk,iloni,jj2)
            zfma(jk,jj,1) = pgd2(jk,ilonr,jj2)
            zfma(jk,jj,2) = pgd2(jk,iloni,jj2)
          enddo
        else
          do jk = 1, gst(gstID)%myLevCount
            ! Southern hemisphere
            zfms(jk,jj2,1) = zfms(jk,jj2,1) + pgd2(jk,ilonr,jj2)
            zfms(jk,jj2,2) = zfms(jk,jj2,2) + pgd2(jk,iloni,jj2)
            zfma(jk,jj2,1) = zfma(jk,jj2,1) - pgd2(jk,ilonr,jj2)
            zfma(jk,jj2,2) = zfma(jk,jj2,2) - pgd2(jk,iloni,jj2)
          enddo
        endif
      enddo

      do jj = 1,gst(gstID)%njlath
        do jk = 1, gst(gstID)%myLevCount
          zfms(jk,jj,1) = dlrwt(jj)*zfms(jk,jj,1)
          zfms(jk,jj,2) = dlrwt(jj)*zfms(jk,jj,2)
          zfma(jk,jj,1) = dlrwt(jj)*zfma(jk,jj,1)
          zfma(jk,jj,2) = dlrwt(jj)*zfma(jk,jj,2)
        enddo
      enddo


      ! 2.3 First one with ALP for all scalar fields and for half the terms
      !     required to define the divergence and vorticity
      call legdir2d_kij(jm,zfms,zfma,dlsp,dlalp,gst(gstID)%njlath,gst(gstID)%ntrunc,gst(gstID)%ntrunc)

      ! 2.4 Transfer the result in the global state
      do jn = jm, gst(gstID)%ntrunc
        ila = gst(gstID)%nind(jm) + jn - jm
        inm = jn - jm
        do jk = 1, gst(gstID)%myLevCount
          psp(ila,1,jk) = dlsp(jk,inm,1)
          psp(ila,2,jk) = dlsp(jk,inm,2)
        enddo
      enddo

    ! End of loop on zonal wavenumbers
    enddo
    !$OMP END PARALLEL DO

  end subroutine reespepar_kij


  subroutine legdir2d(km,pfms,pfma,ddsp,ddalp,klath,ktrunc,ktruncdim)
    implicit none

    ! Arguments:
    integer :: km, ktrunc, ktruncdim, klath
    real(8) :: pfms(gst(gstID)%njlath+1,2, &
                    gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: pfma(gst(gstID)%njlath+1,2, &
                    gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: ddalp(0:ktruncdim,klath)
    real(8) :: ddsp(0:ktruncdim,2,gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)

    ! Locals:
    integer :: jk, jlat, jn, inm, itrunc, inmp1

    itrunc = ktrunc
    if(mod(ktrunc-km+1,2).eq.1) itrunc = ktrunc-1

    if(km.ne.ktrunc)then
       ddsp(0:ktrunc,:,:) = 0.d0
       do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
          do jlat = 1, klath
             do jn = km, itrunc, 2
                inm = jn - km
                inmp1 = inm + 1
                ddsp(inm,  1,jk) = ddsp(inm,  1,jk) + ddalp(inm,  jlat)*pfms(jlat,1,jk)
                ddsp(inmp1,1,jk) = ddsp(inmp1,1,jk) + ddalp(inmp1,jlat)*pfma(jlat,1,jk)
                ddsp(inm,  2,jk) = ddsp(inm,  2,jk) + ddalp(inm,  jlat)*pfms(jlat,2,jk)
                ddsp(inmp1,2,jk) = ddsp(inmp1,2,jk) + ddalp(inmp1,jlat)*pfma(jlat,2,jk)
             enddo
          enddo
       enddo
    end if

    if(mod(ktrunc-km+1,2).eq.1) then
       jn = ktrunc
       inm = jn - km
       ddsp(inm,:,:) = 0.d0
       do jk = gst(gstID)%myLevBeg, gst(gstID)%myLevEnd
          do jlat = 1, klath
             ddsp(inm,1,jk) = ddsp(inm,1,jk) + ddalp(inm,jlat)*pfms(jlat,1,jk )
             ddsp(inm,2,jk) = ddsp(inm,2,jk) + ddalp(inm,jlat)*pfms(jlat,2,jk )
          enddo
       enddo
    end if

  end subroutine legdir2d


  subroutine legdir2d_kij(km,pfms,pfma,ddsp,ddalp,klath,ktrunc,ktruncdim)
    implicit none

    ! Arguments:
    integer :: km, ktrunc, ktruncdim, klath
    real(8) :: pfms(gst(gstID)%myLevCount,gst(gstID)%njlath+1,2)
    real(8) :: pfma(gst(gstID)%myLevCount,gst(gstID)%njlath+1,2)
    real(8) :: ddalp(0:ktruncdim,klath)
    real(8) :: ddsp(gst(gstID)%myLevCount,0:ktruncdim,2)

    ! Locals:
    integer :: jk, jlat, jn, inm, itrunc, inmp1

    itrunc = ktrunc
    if(mod(ktrunc-km+1,2).eq.1) itrunc = ktrunc-1

    if(km.ne.ktrunc)then
       ddsp(:,0:ktrunc,:) = 0.d0
       do jlat = 1, klath
          do jn = km, itrunc, 2
             inm = jn - km
             inmp1 = inm + 1
             do jk = 1, gst(gstID)%myLevCount
                ddsp(jk,inm,  1) = ddsp(jk,inm,  1) + ddalp(inm,  jlat)*pfms(jk,jlat,1)
                ddsp(jk,inmp1,1) = ddsp(jk,inmp1,1) + ddalp(inmp1,jlat)*pfma(jk,jlat,1)
                ddsp(jk,inm,  2) = ddsp(jk,inm,  2) + ddalp(inm,  jlat)*pfms(jk,jlat,2)
                ddsp(jk,inmp1,2) = ddsp(jk,inmp1,2) + ddalp(inmp1,jlat)*pfma(jk,jlat,2)
             enddo
          enddo
       enddo
    end if

    if(mod(ktrunc-km+1,2).eq.1) then
       jn = ktrunc
       inm = jn - km
       ddsp(:,inm,:) = 0.d0
       do jlat = 1, klath
          do jk = 1, gst(gstID)%myLevCount
             ddsp(jk,inm,1) = ddsp(jk,inm,1) + ddalp(inm,jlat)*pfms(jk,jlat,1 )
             ddsp(jk,inm,2) = ddsp(jk,inm,2) + ddalp(inm,jlat)*pfms(jk,jlat,2 )
          enddo
       enddo
    end if

  end subroutine legdir2d_kij


  subroutine leginv2d(km,pfms,pfma,ddsp,ddalp,klath,ktrunc,ktruncdim)
    implicit none

    ! Arguments:
    integer :: km, ktrunc, ktruncdim, klath
    real(8) :: pfms(gst(gstID)%njlath+1,2, &
                    gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: pfma(gst(gstID)%njlath+1,2, &
                    gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)
    real(8) :: ddalp(0:ktruncdim,klath)
    real(8) :: ddsp(0:ktruncdim,2,gst(gstID)%myLevBeg:gst(gstID)%myLevEnd)

    ! Locals:
    integer :: jk, jlat, jn, inm, itrunc, inmp1

    itrunc = ktrunc
    if(mod(ktrunc-km+1,2).eq.1) itrunc = ktrunc-1

    if(km .ne. ktrunc)then
       do jk = gst(gstID)%myLevBeg,gst(gstID)%myLevEnd
          pfms(:,:,jk) = 0.d0
          pfma(:,:,jk) = 0.d0
          do jn = km, itrunc, 2
             inm = jn - km
             inmp1  = inm + 1
             do jlat = 1, klath
                pfms(jlat,1,jk) = pfms(jlat,1,jk) + ddalp(inm,  jlat) * ddsp(inm,  1,jk)
                pfma(jlat,1,jk) = pfma(jlat,1,jk) + ddalp(inmp1,jlat) * ddsp(inmp1,1,jk)
                pfms(jlat,2,jk) = pfms(jlat,2,jk) + ddalp(inm,  jlat) * ddsp(inm,  2,jk)
                pfma(jlat,2,jk) = pfma(jlat,2,jk) + ddalp(inmp1,jlat) * ddsp(inmp1,2,jk)
             enddo
          enddo
       enddo
    end if

    if(mod(ktrunc-km+1,2).eq.1) then
       jn = ktrunc
       if (km .ne. ktrunc) then
          inm = jn - km
          do jk = gst(gstID)%myLevBeg,gst(gstID)%myLevEnd
             do jlat = 1, klath
                pfms(jlat,1,jk) = pfms(jlat,1,jk) + ddalp(inm,jlat) * ddsp(inm,1,jk)
                pfms(jlat,2,jk) = pfms(jlat,2,jk) + ddalp(inm,jlat) * ddsp(inm,2,jk)
             enddo
          enddo
       else
          inm = jn - km
          do jk = gst(gstID)%myLevBeg,gst(gstID)%myLevEnd
             pfms(:,:,jk) = 0.d0
             pfma(:,:,jk) = 0.d0
             do jlat = 1, klath
                pfms(jlat,1,jk) = ddalp(inm,jlat) * ddsp(inm,1,jk)
                pfms(jlat,2,jk) = ddalp(inm,jlat) * ddsp(inm,2,jk)
             enddo
          enddo
       end if
    end if

  end subroutine leginv2d


  subroutine leginv2d_kij(km,pfms,pfma,ddsp,ddalp,klath,ktrunc,ktruncdim)
    implicit none

    ! Arguments:
    integer :: km, ktrunc, ktruncdim, klath
    real(8) :: pfms(gst(gstID)%myLevCount,gst(gstID)%njlath+1,2)
    real(8) :: pfma(gst(gstID)%myLevCount,gst(gstID)%njlath+1,2)
    real(8) :: ddalp(0:ktruncdim,klath)
    real(8) :: ddsp(gst(gstID)%myLevCount,0:ktruncdim,2)

    ! Locals:
    integer :: jk, jlat, jn, inm, itrunc, inmp1

    itrunc = ktrunc
    if(mod(ktrunc-km+1,2).eq.1) itrunc = ktrunc-1

    if(km .ne. ktrunc)then
       pfms(:,:,:) = 0.d0
       pfma(:,:,:) = 0.d0
       do jn = km, itrunc, 2
          inm = jn - km
          inmp1  = inm + 1
          do jlat = 1, klath
             do jk = 1, gst(gstID)%myLevCount
                pfms(jk,jlat,1) = pfms(jk,jlat,1) + ddalp(inm,  jlat) * ddsp(jk,inm,  1)
                pfma(jk,jlat,1) = pfma(jk,jlat,1) + ddalp(inmp1,jlat) * ddsp(jk,inmp1,1)
                pfms(jk,jlat,2) = pfms(jk,jlat,2) + ddalp(inm,  jlat) * ddsp(jk,inm,  2)
                pfma(jk,jlat,2) = pfma(jk,jlat,2) + ddalp(inmp1,jlat) * ddsp(jk,inmp1,2)
             enddo
          enddo
       enddo
    end if

    if(mod(ktrunc-km+1,2).eq.1) then
       jn = ktrunc
       if (km .ne. ktrunc) then
          inm = jn - km
          do jlat = 1, klath
             do jk = 1, gst(gstID)%myLevCount
                pfms(jk,jlat,1) = pfms(jk,jlat,1) + ddalp(inm,jlat) * ddsp(jk,inm,1)
                pfms(jk,jlat,2) = pfms(jk,jlat,2) + ddalp(inm,jlat) * ddsp(jk,inm,2)
             enddo
          enddo
       else
          inm = jn - km
          pfms(:,:,:) = 0.d0
          pfma(:,:,:) = 0.d0
          do jlat = 1, klath
             do jk = 1, gst(gstID)%myLevCount
                pfms(jk,jlat,1) = ddalp(inm,jlat) * ddsp(jk,inm,1)
                pfms(jk,jlat,2) = ddalp(inm,jlat) * ddsp(jk,inm,2)
             enddo
          enddo
       end if
    end if

  end subroutine leginv2d_kij


  subroutine allocate_comleg
    !
    !:Purpose: Subroutine for initializing the Legendre transform
    implicit none

    allocate(gst(gstID)%rmu(gst(gstID)%nj))  
    allocate(gst(gstID)%rwt(gst(gstID)%nj))
    allocate(gst(gstID)%rwocs(gst(gstID)%nj))
    allocate(gst(gstID)%r1mu2(gst(gstID)%nj))
    allocate(gst(gstID)%rsqm2(gst(gstID)%nj))
    allocate(gst(gstID)%rcolat(gst(gstID)%nj))
    allocate(gst(gstID)%r1qm2(gst(gstID)%nj))
    allocate(gst(gstID)%r1mui(gst(gstID)%nj))
    allocate(gst(gstID)%r1mua(gst(gstID)%nj))
    allocate(gst(gstID)%rlati((-1):(gst(gstID)%nj+2)))
    allocate(gst(gstID)%nind(0:gst(gstID)%ntrunc))
    allocate(gst(gstID)%nindrh(0:gst(gstID)%ntrunc))
    allocate(gst(gstID)%nclm(0:gst(gstID)%ntrunc))

  end subroutine allocate_comleg


  subroutine suleg(lverbose_opt)
    !
    !:Purpose: To initializethe Gaussian latitudes, weights and related
    !          quantities
    implicit none

    ! Arguments:
    logical, optional :: lverbose_opt

    ! Locals:
    logical :: lverbose
    integer :: jlat, jm
    real(8) :: zpisu2
    external gauss

    ! Some explanation:
    ! rmu = sin(latitude)
    ! colat = pi/2 - latitude
    ! rwt = gauss weight (complicated)
    ! rwocs = rwt / cos(latitude)^2

    if(present(lverbose_opt)) then
      lverbose = lverbose_opt
    else
      lverbose = .false.
    endif

    if(mpi_myid.eq.0) write(*,fmt='(//,6(" ***********"))')
    if(mpi_myid.eq.0) write(*,*)'     SULEG: initialisation of Gaussian', &
         ' latitudes, weights, etc...'
    if(mpi_myid.eq.0) write(*,fmt='(6(" ***********"))')

    !     1. GAUSSIAN LATITUDES AND WEIGHTS OVER AN HEMISPHERE
    !     -------------------------------------------------

    call gauss8(gst(gstID)%njlath,gst(gstID)%rmu(1),gst(gstID)%rwt(1), &
                gst(gstID)%rsqm2(1),gst(gstID)%rcolat(1),gst(gstID)%rwocs(1), &
                gst(gstID)%r1qm2(1),gst(gstID)%r1mui(1),gst(gstID)%r1mu2(1))

    do jlat = 1, gst(gstID)%njlath
       gst(gstID)%rlati(jlat) = asin(gst(gstID)%rmu(jlat))
       gst(gstID)%r1mua(jlat) = ec_r1sa*gst(gstID)%r1mui(jlat)
    enddo

    !     2. COMPLETION FOR THE SOUTHERN HEMISPHERE
    !     --------------------------------------

    do jlat = gst(gstID)%njlath +1, gst(gstID)%nj
       gst(gstID)%rmu(jlat)   =  -gst(gstID)%rmu(2*gst(gstID)%njlath +1 - jlat)
       gst(gstID)%rwocs(jlat) =   gst(gstID)%rwocs(2*gst(gstID)%njlath +1 - jlat)
       gst(gstID)%r1mu2(jlat) =   gst(gstID)%r1mu2(2*gst(gstID)%njlath +1 - jlat)
       gst(gstID)%rsqm2(jlat) =   gst(gstID)%rsqm2(2*gst(gstID)%njlath +1 - jlat)
       gst(gstID)%r1qm2(jlat) =   gst(gstID)%r1qm2(2*gst(gstID)%njlath +1 - jlat)
       gst(gstID)%r1mui(jlat) =   gst(gstID)%r1mui(2*gst(gstID)%njlath +1 - jlat)
       gst(gstID)%r1mua(jlat) =   gst(gstID)%r1mua(2*gst(gstID)%njlath +1 - jlat)
       gst(gstID)%rwt(jlat)   =   gst(gstID)%rwt(2*gst(gstID)%njlath +1 - jlat)
       gst(gstID)%rlati(jlat) = - gst(gstID)%rlati (2*gst(gstID)%njlath +1 - jlat)
    enddo

    zpisu2 = MPC_PI_R8/2.d0
    do jlat = 1, gst(gstID)%nj
       gst(gstID)%rcolat(jlat) = zpisu2 - gst(gstID)%rlati(jlat)
    enddo

    !*    3. Overdimensioning for interpolation

    gst(gstID)%rlati(-1) =   MPC_PI_R8-gst(gstID)%rlati(1)
    gst(gstID)%rlati(0) =   MPC_PI_R8*.5d0
    gst(gstID)%rlati(gst(gstID)%nj+1) =-MPC_PI_R8*.5d0
    gst(gstID)%rlati(gst(gstID)%nj+2) =-MPC_PI_R8-gst(gstID)%rlati(gst(gstID)%nj)

    !*    4. Print the content of GAUS

    if(lverbose.and.mpi_myid.eq.0) write(*,fmt='(" JLAT:",4X," RLATI",8X,"RCOLAT",8X,"RMU",10X ,"RWT",12X,"RW0CS")')
    do jlat = 1, gst(gstID)%nj
       if(lverbose.and.mpi_myid.eq.0) write(*,fmt='(2X,I4,5(2X,G23.16))')  &
            jlat,gst(gstID)%rlati(jlat),gst(gstID)%rcolat(jlat), gst(gstID)%rmu(jlat)  &
            ,gst(gstID)%rwt(jlat),gst(gstID)%rwocs(jlat)
    enddo

    if(lverbose.and.mpi_myid.eq.0) write(*,fmt='(//," JLAT:",4X,"R1MU2",8X,"RSQM2",9X,"R1QM2",10X,"R1MUI",10X,"R1MUA")')

    do jlat = 1, gst(gstID)%nj
       if(lverbose.and.mpi_myid.eq.0)  &
         write(*,fmt='(2X,I4,5(2X,G23.16))') jlat,gst(gstID)%r1mu2(jlat),gst(gstID)%rsqm2(jlat),gst(gstID)%r1qm2(jlat)  &
              ,gst(gstID)%r1mui(jlat),gst(gstID)%r1mua(jlat)
    enddo

    !*    5.  Positioning within spectral arrays

    do jm = 0, gst(gstID)%ntrunc
       gst(gstID)%nind(jm)   = jm*(gst(gstID)%ntrunc+1) - (jm*(jm-1))/2 + 1
       gst(gstID)%nindrh(jm) = jm*(gst(gstID)%ntrunc+1) + 1
       gst(gstID)%nclm(jm)   = gst(gstID)%ntrunc - jm + 1
    enddo

    if(lverbose.and.mpi_myid.eq.0) write(*,fmt='(/," NIND(0:NTRUNC):",/,10(2X,I8))')  &
         (gst(gstID)%nind(jm),jm=0,gst(gstID)%ntrunc)
    if(lverbose.and.mpi_myid.eq.0) write(*,fmt='(" NINDRH(0:NTRUNC):",/,10(2X,I8))')  &
         (gst(gstID)%nindrh(jm),jm=0,gst(gstID)%ntrunc)
    if(lverbose.and.mpi_myid.eq.0) write(*,fmt='("   NCLM(0:NTRUNC):",/,10(2X,I8))')  &
         (gst(gstID)%nclm(jm),jm=0,gst(gstID)%ntrunc)

  end subroutine suleg


  subroutine gauss8(nracp,racp,pg,sia,rad,pgssin2,sinm1,sinm2,sin2)
    implicit none

    ! Arguments:
    integer :: nracp
    real(8) :: racp(*)
    real(8) :: pg(*)
    real(8) :: sia(*)
    real(8) :: rad(*)
    real(8) :: pgssin2(*)
    real(8) :: sinm1(*)
    real(8) :: sinm2(*)
    real(8) :: sin2(*)

    ! Locals:
    real(8) :: xlim,pi,fi,fi1,fn,dot,dn,dn1,a,b,c,g,gm,gp,gt,ractemp,gtemp
    integer :: i,ir,irm,irp

    xlim = 1.d-13
    pi = 4.d0*atan(1.d0)
    ir = 2*nracp
    fi = dble(ir)
    fi1 = fi+1.d0
    fn = dble(nracp)

    do i = 1,nracp
       dot = dble(i-1)
       racp(i) = -pi*.5d0*(dot+.5d0)/fn + pi*.5d0
       racp(i) =  sin(racp(i))
    enddo

    dn = fi/sqrt(4.d0*fi*fi-1.d0)
    dn1 = fi1/sqrt(4.d0*fi1*fi1-1.d0)
    a = dn1*fi
    b = dn*fi1
    irp = ir + 1
    irm = ir -1

    do i = 1,nracp
42     call ordleg8(g,racp(i),ir)
       call ordleg8(gm,racp(i),irm)
       call ordleg8(gp,racp(i),irp)
       gt = (a*gp-b*gm)/(racp(i)*racp(i)-1.d0)
       ractemp = racp(i) - g/gt
       gtemp = racp(i) - ractemp
       racp(i) = ractemp
       if( abs(gtemp).gt.xlim) go to 42
    enddo

    do i = 1,nracp
       a = 2.d0*(1.d0-racp(i)**2)
       call ordleg8(b,racp(i),irm)
       b = b*b*fi*fi
       pg(i) = a*(fi-.5d0)/b
       rad(i) =   acos(racp(i))
       sia(i) =  sin(rad(i))
       c = (sia(i))**2
       sinm1(i) = 1.d0/sia(i)
       sinm2(i) = 1.d0/c
       pgssin2(i) = pg(i)/c
       sin2(i) = c
    enddo

  end subroutine gauss8


  subroutine ordleg8(sx,coa,ir)
    implicit none

    ! Arguments:
    real(8) :: sx
    real(8) :: coa
    integer :: ir

    ! Locals:
    integer :: n,kk,k,n1,irpp,irppm
    real(8) :: pi,sqr2,delta,sia,theta,c1,c4,s1,ang,fk,fn,fn2,fn2sq,a,b

    pi    = 4.d0*atan(1.d0)
    sqr2  = sqrt(2.d0)
    irpp  = ir   + 1
    irppm = irpp - 1
    delta = acos(coa)
    sia   = sin(delta)

    theta = delta
    c1    = sqr2

    do n = 1,irppm
       fn2   = dble(2*n)
       fn2sq = fn2*fn2
       c1    =  c1*sqrt(1.d0 - 1.d0/fn2sq)
    enddo

    n   = irppm
    fn  = dble(n)
    ang = fn*theta
    s1  = 0.d0
    c4  = 1.d0
    a   =-1.d0
    b   = 0.d0
    n1  = n+1

    do kk = 1,n1,2
       k   = kk-1
       if (k.eq.n) c4 = 0.5d0*c4
       s1  = s1+c4* cos(ang)
       a   =  a+2.d0
       b   =  b+1.d0
       fk  = dble(k)
       ang = theta*(fn-fk-2.d0)
       c4  = ( a * (fn-b+1.d0) / (b*(fn2-a)) )*c4
    enddo

    sx = s1*c1

  end subroutine ordleg8


  subroutine sualp
    implicit none

    ! Locals:
    integer :: jj,jlat,jm,jn,ilat
    integer :: ilarh,ila,ilatbd
    real(8) :: dlalp(gst(gstID)%nlarh,nlatbd)
    real(8) :: dldalp(gst(gstID)%nlarh,nlatbd)
    !     
    !     Memory allocation for Legendre polynomials
    !     
    if(mpi_myid.eq.0) write(*,*) 'allocating dalp:',gst(gstID)%nla,gst(gstID)%njlath,gst(gstID)%nla*gst(gstID)%njlath
    allocate( gst(gstID)%dalp(gst(gstID)%nla,gst(gstID)%njlath) )
    allocate( gst(gstID)%dealp(gst(gstID)%nla,gst(gstID)%njlath))
    if(mpi_myid.eq.0) write(*,*) 'succeeded'

    latitudes: do jlat = 1, gst(gstID)%njlath, nlatbd
       ilatbd = min(nlatbd,gst(gstID)%njlath - jlat + 1)

       if(ilatbd.eq.8) then
          call allp(dlalp,dldalp,gst(gstID)%rmu(jlat),gst(gstID)%nclm(0),gst(gstID)%ntrunc,ilatbd)
       else
          call allp2(dlalp,dldalp,gst(gstID)%rmu(jlat),gst(gstID)%ntrunc,ilatbd)
       endif

       do jm = 0,gst(gstID)%ntrunc
          do jn = jm,gst(gstID)%ntrunc
             ila = gst(gstID)%nind(jm) + jn -jm
             ilarh = gst(gstID)%nindrh(jm) + jn-jm
             do jj = 1,ilatbd
                ilat = jlat+jj-1
                gst(gstID)%dalp (ila,jlat+jj-1) = dlalp (ilarh,jj)
                gst(gstID)%dealp(ila,jlat+jj-1) = dldalp(ilarh,jj)
             enddo
          enddo
       enddo
    enddo latitudes

  end subroutine sualp


  subroutine getalp(ddalp,dddalp,klath,ktrunc,ktruncdim ,km)
    implicit none

    ! Arguments:
    real(8) :: ddalp(0:ktruncdim,klath)
    real(8) :: dddalp(0:ktruncdim,klath)
    integer :: klath
    integer :: ktrunc
    integer :: ktruncdim
    integer :: km

    ! Locals: 
    integer :: ila,ind
    integer :: jlat,jn, jlen

    do jlat = 1,klath
       do jlen = 0, ktrunc
          ddalp(jlen,jlat) = 0.d0
          dddalp(jlen,jlat)= 0.d0
       end do
    end do

    do jlat = 1, klath
       do jn = km, ktrunc
          ila = gst(gstID)%nind(km) + jn-km
          ind = jn-km
          ddalp(ind,jlat) =  gst(gstID)%dalp(ila,jlat)
          dddalp(ind,jlat) = gst(gstID)%dealp(ila,jlat)
       end do
    end do

  end subroutine getalp


  subroutine allp( p , g , x , lr , r , nlatp) 

    implicit none

    ! Arguments:
    integer :: r, nlatp, lr(0:r)
    real(8) :: p(0:r,0:r,nlatp) , g(0:r,0:r,nlatp) 
    real(8) :: x(nlatp) 

    ! Locals:
    real(8) :: onehalf   
    real(8) :: xp , xp2,  p0, enm, fnm
    integer :: ilat , m , l , n

    data onehalf /0.5d0/

    do ilat = 1,nlatp
       xp2 = sqrt( 1.0d0 - x(ilat) ** 2 ) 
       p(0,0,ilat) = sqrt(onehalf) 
       do m = 1,r 
          xp = real(m,8)
          p(0,m,ilat) = sqrt( (2.0d0*xp+1.0d0)/(2.0d0*xp) ) * xp2 * p(0,m-1,ilat)
       enddo
    enddo

    do ilat = 1,nlatp
       do m = 0,r 
          xp = real(m,8)
          g(0,m,ilat) = - x(ilat)*xp * p(0,m,ilat) 
       enddo
    enddo
    do n = 1,r
       do m = 0,lr(n)-1
          l =  1
          p0 = real(m+n,8)
          xp = real(m,8)
          enm = sqrt( ((p0*p0-xp*xp)*(2.0d0*p0+1.0d0))/(2.0d0*p0-1.0d0) )
          fnm = sqrt( (2.0d0*p0+1.0d0)/((p0*p0-xp*xp)*(2.0d0*p0-1.0d0)) )

          p(n,m,l) = ( x(l) * p0 * p(n-1,m,l) -  g(n-1,m,l) ) * fnm 
          g(n,m,l) = enm * p(n-1,m,l) - x(l) * p0 * p(n,m,l) 
          l = l + 1
          p(n,m,l) = ( x(l) * p0 * p(n-1,m,l) -  g(n-1,m,l) ) * fnm 
          g(n,m,l) = enm * p(n-1,m,l) - x(l) * p0 * p(n,m,l) 
          l = l + 1
          p(n,m,l) = ( x(l) * p0 * p(n-1,m,l) -  g(n-1,m,l) ) * fnm 
          g(n,m,l) = enm * p(n-1,m,l) - x(l) * p0 * p(n,m,l) 
          l = l + 1
          p(n,m,l) = ( x(l) * p0 * p(n-1,m,l) -  g(n-1,m,l) ) * fnm 
          g(n,m,l) = enm * p(n-1,m,l) - x(l) * p0 * p(n,m,l) 
          l = l + 1
          p(n,m,l) = ( x(l) * p0 * p(n-1,m,l) -  g(n-1,m,l) ) * fnm 
          g(n,m,l) = enm * p(n-1,m,l) - x(l) * p0 * p(n,m,l) 
          l = l + 1
          p(n,m,l) = ( x(l) * p0 * p(n-1,m,l) -  g(n-1,m,l) ) * fnm 
          g(n,m,l) = enm * p(n-1,m,l) - x(l) * p0 * p(n,m,l) 
          l = l + 1
          p(n,m,l) = ( x(l) * p0 * p(n-1,m,l) -  g(n-1,m,l) ) * fnm 
          g(n,m,l) = enm * p(n-1,m,l) - x(l) * p0 * p(n,m,l) 
          l = l + 1
          p(n,m,l) = ( x(l) * p0 * p(n-1,m,l) -  g(n-1,m,l) ) * fnm 
          g(n,m,l) = enm * p(n-1,m,l) - x(l) * p0 * p(n,m,l) 
          l = l + 1
       enddo
    enddo

  end subroutine allp


  subroutine allp2( p , g , x , r , nlatp)

    implicit none

    ! Arguments:
    integer :: r, nlatp
    real(8) :: p(0:r,0:r,nlatp) , g(0:r,0:r,nlatp) 
    real(8) :: x(nlatp)

    ! Locals:
    real(8) :: onehalf   
    real(8) :: xp , xp2,  p0, enm, fnm
    integer :: ilat , m , l , n, jlat

    data onehalf /0.5d0/

    do ilat = 1,nlatp
       xp2 = sqrt( 1.0d0 - x(ilat) ** 2 )
       p(0,0,ilat) = sqrt(onehalf)
       do m = 1,r
          xp = real(m,8)
          p(0,m,ilat) = sqrt( (2.0d0*xp+1.0d0)/(2.0d0*xp) ) * xp2 * p(0,m-1,ilat)
       enddo
    enddo

    do ilat = 1,nlatp
       do m = 0,r
          xp = real(m,8)
          g(0,m,ilat) = - x(ilat)*xp * p(0,m,ilat)
       enddo
    enddo

    do n = 1,r
       do m = 0, r
          p0 = real(m+n,8)
          xp = real(m,8)
          enm = sqrt( ((p0*p0-xp*xp)*(2.0d0*p0+1.0d0))/(2.0d0*p0-1.0d0) )
          fnm = sqrt( (2.0d0*p0+1.0d0)/((p0*p0-xp*xp)*(2.0d0*p0-1.0d0)) )

          do jlat = 1, nlatp
             l = jlat
             p(n,m,l) = ( x(l) * p0 * p(n-1,m,l) - g(n-1,m,l) ) * fnm
             g(n,m,l) = enm * p(n-1,m,l) - x(l) * p0 * p(n,m,l)
          enddo
       enddo
    enddo

  end subroutine allp2


  subroutine gst_zlegpol(gstID_in)
    !
    !:Purpose: To evaluate Legendre polynomials restricted to (n,m) = (n,0)
    implicit none

    ! Arguments:
    integer :: gstID_in

    ! Loclas:
    integer :: jn, jlat
    real(8) :: dlfact1, dlfact2, dln
    real(8) :: dlnorm(0:gst(gstID)%ntrunc)

    allocate(gst(gstID_in)%zleg(0:gst(gstID)%ntrunc,gst(gstID)%nj))

    do jlat = 1, gst(gstID_in)%nj
       gst(gstID_in)%zleg(0,gst(gstID_in)%nj-jlat+1) = sqrt(0.5d0)
       gst(gstID_in)%zleg(1,gst(gstID_in)%nj-jlat+1) = sqrt(1.5d0)*gst(gstID_in)%rmu(jlat)
    enddo

    do jn = 0, gst(gstID_in)%ntrunc
       dln = 1.d0*real(jn,8)
       dlnorm(jn) = dsqrt((2.d0*dln + 1.d0)/2.d0)
    enddo

    do jn = 1, gst(gstID_in)%ntrunc-1
       dln = real(jn,8)
       dlfact1 = ((2.d0*dln+1.d0)/(dln+1.d0))*(dlnorm(jn+1)/dlnorm(jn))
       dlfact2 = (dln/(dln+1.d0))*(dlnorm(jn+1)/dlnorm(jn-1))
       do jlat = 1,gst(gstID_in)%nj
          gst(gstID_in)%zleg(jn+1,gst(gstID_in)%nj-jlat+1) =   &
                  dlfact1*gst(gstID_in)%rmu(jlat)*dble(gst(gstID_in)%zleg(jn,gst(gstID_in)%nj-jlat+1))   &
                - dlfact2*dble(gst(gstID_in)%zleg(jn-1,gst(gstID_in)%nj-jlat+1))
       enddo
    enddo

  end subroutine gst_zlegpol


  subroutine gst_zlegdir(gstID_in,pf,pn,klev)
    !
    !:Purpose: Direct Legendre transform restricted to
    !
    implicit none

    ! Arguments:
    integer :: gstID_in
    real(8) :: pf(gst(gstID_in)%nj,klev) ! PF(NJ,KLEV): field in physical space
    real(8) :: pn(0:gst(gstID_in)%ntrunc,klev) ! PN(0:ntrunc, KLEV): spectral coefficients
    integer :: klev ! number of fields to transform

    ! Locals:
    integer :: jlat, jn
    real(8), allocatable :: zwork(:,:)

    allocate(zwork(0:gst(gstID_in)%ntrunc,gst(gstID_in)%nj))

    do jlat = 1, gst(gstID_in)%nj
       do jn = 0, gst(gstID_in)%ntrunc
          zwork(jn,jlat) = gst(gstID_in)%zleg(jn,jlat)*gst_getRWT(jlat,gstID_in)
       end do
    end do

    call dgemm('N','N',gst(gstID_in)%ntrunc+1, klev, gst(gstID_in)%nj, 1.0d0, zwork(0,1),  &
                       gst(gstID_in)%ntrunc+1, pf(1,1),  &
                       gst(gstID_in)%nj, 0.0d0, pn(0,1), gst(gstID_in)%ntrunc+1) 

    deallocate(zwork)

  end subroutine gst_zlegdir


  subroutine gst_zleginv(gstID_in,pf,pn,klev)
    !
    !:Purpose: Direct Legendre transform restricted to fields that vary with
    !          latitude only
    !
    implicit none

    ! Arguments:
    integer :: gstID_in
    real(8) :: pf(gst(gstID_in)%nj,klev) ! PF(KNJDIM,KLEVDIM)  : field in physical space
    real(8) :: pn(0:gst(gstID_in)%ntrunc,klev) ! PN(0:KNDIM, KLEVDIM): spectral coefficients
    integer :: klev ! number of fields to transform

    ! Locals:
    integer :: jlat, jn
    real(8), allocatable :: zwork(:,:)

    allocate(zwork(0:gst(gstID_in)%ntrunc,gst(gstID_in)%nj))

    do jlat = 1, gst(gstID_in)%nj
       do jn = 0, gst(gstID_in)%ntrunc
          zwork(jn,jlat) = gst(gstID_in)%zleg(jn,jlat)
       end do
    end do

    call dgemm('T','N',gst(gstID_in)%nj, klev, gst(gstID_in)%ntrunc+1, 1.0d0, zwork(0,1),  &
                       gst(gstID_in)%ntrunc+1, pn(0,1),   &
                       gst(gstID_in)%ntrunc+1, 0.0d0, pf(1,1), gst(gstID_in)%nj) 

    deallocate(zwork)

  end subroutine gst_zleginv


! ---------------------------------------------
! FFT subroutines
! ---------------------------------------------

  subroutine fft3dvar(pgd,kdir)
    implicit none

    ! Arguments:
    real(8) :: pgd(:,:,:)
    integer :: kdir

    ! Locals:
    integer :: kfield,ni,nj
    real(8), allocatable :: pgd2(:,:)
    integer :: ji,jj,jk,ijump,i

    ni = size(pgd,1)
    nj = size(pgd,2)
    kfield = size(pgd,3)
    ijump = ni + 2

    i = ni
    call ngfft( i )
    if ( i.ne.ni ) then
       write(*,*) 'fft3dvar: NI = ',ni,' I = ',i
       call utl_abort('fft3dvar: vector length is not compatible with FFT')
    endif
    call setfft8(ni)

    allocate(pgd2(ni+2,nj))

    ! Copy over input data into over-dimensioned array
    !$OMP PARALLEL DO PRIVATE (jj,jk,ji,pgd2)
    do jk=1,kfield
      do jj=1, nj
        do ji=1,ni
          pgd2(ji,jj)=pgd(ji,jj,jk)
        enddo
        do ji=ni+1,ni+2
          pgd2(ji,jj)=0.0d0
        enddo
      enddo

      call ffft8(pgd2(1,1),1,ijump,nj,kdir)
      !*     subroutine ffft8( a, inc, jump, lot, isign )
      !*     a      is the array containing input & output data
      !*     inc    is the increment within each data 'vector'
      !*            (e.g. inc=1 for consecutively stored data)
      !*     jump   is the increment between the start of each data vector
      !*     lot    is the number of data vectors
      !*     isign  = +1 for transform from spectral to gridpoint
      !*            = -1 for transform from gridpoint to spectral

      do jj=1,nj
        do ji=1,ni
          pgd(ji,jj,jk)=pgd2(ji,jj)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(pgd2)

  end subroutine fft3dvar

  subroutine fft3dvar_kij(pgd,kdir)
    implicit none

    ! Arguments:
    real(8) :: pgd(:,:,:)
    integer :: kdir

    ! Locals:
    integer :: kfield,ni,nj
    real(8), allocatable :: pgd2(:,:)
    integer :: ji,jj,jk,inc,ijump,i

    kfield = size(pgd,1)
    ni = size(pgd,2)
    nj = size(pgd,3)

    i = ni
    call ngfft( i )
    if ( i.ne.ni ) then
       write(*,*) 'fft3dvar: NI = ',ni,' I = ',i
       call utl_abort('fft3dvar: vector length is not compatible with FFT')
    endif
    call setfft8(ni)

    inc = kfield
    ijump = 1
    allocate(pgd2(kfield,ni+2))

    ! Copy over input data into over-dimensioned array
    !$OMP PARALLEL DO PRIVATE (jj,jk,ji,pgd2)
    do jj=1, nj
      do ji=1,ni
        do jk=1,kfield
          pgd2(jk,ji)=pgd(jk,ji,jj)
        enddo
      enddo
      do ji=ni+1,ni+2
        do jk=1,kfield
          pgd2(jk,ji)=0.0d0
        enddo
      enddo

      call ffft8(pgd2(1,1),inc,ijump,kfield,kdir)
      !*     subroutine ffft8( a, inc, jump, lot, isign )
      !*     a      is the array containing input & output data
      !*     inc    is the increment within each data 'vector'
      !*            (e.g. inc=1 for consecutively stored data)
      !*     jump   is the increment between the start of each data vector
      !*     lot    is the number of data vectors
      !*     isign  = +1 for transform from spectral to gridpoint
      !*            = -1 for transform from gridpoint to spectral

      do ji=1,ni
        do jk=1,kfield
          pgd(jk,ji,jj)=pgd2(jk,ji)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(pgd2)

  end subroutine fft3dvar_kij

  subroutine ngfft( n )
    implicit none

    ! Arguments:
    integer n

    ! Locals:
    integer l
    parameter ( l = 3 )
    integer k( l ) , m
    data m , k / 8 , 2 , 3 , 5 /

    integer i,j

    if ( n.le.m ) n = m + 1
    n = n - 1
1   n = n + 1
    i = n
2   do 3 j=1,l
       if( mod(i,k(j)) .eq. 0 ) go to 4
3   continue
    go to 1
4   i = i/k(j)
    if( i .ne. 1 ) go to 2

  end subroutine ngfft


end module globalSpectralTransform_mod
