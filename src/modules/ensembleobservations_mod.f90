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

MODULE ensembleObservations_mod
  ! MODULE ensembleObservations (prefix='eob' category='2. High-level data objects')
  !
  ! :Purpose: Store and manipulate ensemble of quanitites in observation space.
  !           This module uses the kdtree2 module for efficiently finding the
  !           nearest observations within the local volume.
  !
  use kdtree2_mod
  use columnData_mod
  use tovs_nl_mod
  use rttov_types, only: rttov_transmission
  use parkind1, only: jpim, jprb
  use ramDisk_mod
  use mpi_mod
  use mpivar_mod
  use obsSpaceData_mod
  use randomNumber_mod
  use mathPhysConstants_mod
  use physicsFunctions_mod
  use utilities_mod
  use earthconstants_mod
  use bufr_mod
  use codePrecision_mod
  implicit none
  save
  private

  ! public types
  public :: struct_eob

  ! public procedures
  public :: eob_allocate, eob_deallocate, eob_allGather, eob_getLocalBodyIndices
  public :: eob_setYb, eob_setDeterYb, eob_setLatLonObs, eob_setObsErrInv, eob_setMeanOMP
  public :: eob_setHPHT, eob_calcAndRemoveMeanYb, eob_setLogPres, eob_copy, eob_zero
  public :: eob_calcRandPert
  public :: eob_backgroundCheck, eob_huberNorm, eob_rejectRadNearSfc
  public :: eob_writeToFiles, eob_readFromFiles

  integer, parameter :: maxNumLocalObsSearch = 500000
  integer,external   :: get_max_rss

  type struct_eob
    logical                       :: allocated = .false.
    integer                       :: numMembers       ! number of ensemble members
    integer                       :: numObs           ! number of observations
    type(struct_obs), pointer     :: obsSpaceData     ! pointer to obsSpaceData object
    real(8), allocatable          :: lat(:), lon(:)   ! lat/lon of observation
    real(8), allocatable          :: logPres(:)       ! ln(pres) of obs, used for localization
    real(8), allocatable          :: obsErrInv(:)     ! inverse of obs error variances
    real(4), allocatable          :: Yb_r4(:,:)       ! background ensemble perturbation in obs space
    real(4), allocatable          :: randPert_r4(:,:) ! unbiased random perturbations with covariance equal to R
    real(8), allocatable          :: meanYb(:)        ! ensemble mean background state in obs space
    real(8), allocatable          :: deterYb(:)       ! deterministic background state in obs space
    real(8), allocatable          :: obsValue(:)      ! the observed value
    integer, allocatable          :: assFlag(:)       ! assimilation flag
  end type struct_eob

  type(kdtree2), pointer :: tree => null()

CONTAINS

  !--------------------------------------------------------------------------
  ! eob_allocate
  !--------------------------------------------------------------------------
  subroutine eob_allocate(ensObs, numMembers, numObs, obsSpaceData)
    !
    ! :Purpose: Allocate an ensObs object
    !
    implicit none

    ! arguments
    type(struct_eob)         :: ensObs
    integer                  :: numMembers
    integer                  :: numObs
    type(struct_obs), target :: obsSpaceData

    if ( ensObs%allocated ) then
      write(*,*) 'eob_allocate: this object is already allocated, deallocating first.'
      call eob_deallocate( ensObs )
    end if

    ensObs%obsSpaceData => obsSpaceData
    ensObs%numMembers = numMembers
    ensObs%numObs     = numObs

    allocate( ensObs%lat(ensObs%numObs) )
    allocate( ensObs%lon(ensObs%numObs) )
    allocate( ensObs%logPres(ensObs%numObs) )
    allocate( ensObs%obsValue(ensObs%numObs) )
    allocate( ensObs%obsErrInv(ensObs%numObs) )
    allocate( ensObs%Yb_r4(ensObs%numMembers,ensObs%numObs) )
    allocate( ensObs%randPert_r4(ensObs%numMembers,ensObs%numObs) )
    allocate( ensObs%meanYb(ensObs%numObs) )
    allocate( ensObs%deterYb(ensObs%numObs) )
    allocate( ensObs%assFlag(ensObs%numObs) )

    ensObs%allocated = .true.

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine eob_allocate

  !--------------------------------------------------------------------------
  ! eob_deallocate
  !--------------------------------------------------------------------------
  subroutine eob_deallocate( ensObs )
    implicit none

    ! arguments
    type(struct_eob) :: ensObs

    if ( .not. ensObs%allocated ) return

    deallocate( ensObs%lat )
    deallocate( ensObs%lon )
    deallocate( ensObs%logPres )
    deallocate( ensObs%obsValue )
    deallocate( ensObs%obsErrInv )
    deallocate( ensObs%Yb_r4 )
    deallocate( ensObs%randPert_r4 )
    deallocate( ensObs%meanYb )
    deallocate( ensObs%deterYb )
    deallocate( ensObs%assFlag )

    ensObs%allocated = .false.

  end subroutine eob_deallocate

  !--------------------------------------------------------------------------
  ! eob_zero
  !--------------------------------------------------------------------------
  subroutine eob_zero(ensObs)
    !
    ! :Purpose: Initialize an ensObs object to zero
    !
    implicit none

    ! arguments
    type(struct_eob)         :: ensObs

    if ( .not.ensObs%allocated ) then
      call utl_abort('eob_zero: this object is not allocated')
    end if

    ensObs%lat(:)        = 0.0d0
    ensObs%lon(:)        = 0.0d0
    ensObs%logPres(:)    = 0.0d0
    ensObs%obsValue(:)   = 0.0d0
    ensObs%obsErrInv(:)  = 0.0d0
    ensObs%Yb_r4(:,:)    = 0.0
    ensObs%randPert_r4(:,:) = 0.0
    ensObs%meanYb(:)     = 0.0d0
    ensObs%deterYb(:)    = 0.0d0
    ensObs%assFlag(:)    = 0

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine eob_zero

  !--------------------------------------------------------------------------
  ! eob_clean (private routine)
  !--------------------------------------------------------------------------
  subroutine eob_clean(ensObs)
    !
    ! :Purpose: Remove all obs from the ensObs object that are not 
    !           flagged for assimilation. All arrays will be reallocated
    !           in place to the smaller size after cleaning. 
    !
    ! :Note: After this procedure, the link with the bodyIndex of obsSpaceData will
    !        be lost!
    !
    implicit none

    ! arguments
    type(struct_eob) :: ensObs

    ! locals
    integer :: obsIndex, obsCleanIndex, numObsClean
    type(struct_eob) :: ensObsClean

    call eob_setAssFlag(ensObs)

    numObsClean = 0
    do obsIndex = 1, ensObs%numObs
      if (ensObs%assFlag(obsIndex) == 1) numObsClean = numObsClean + 1
    end do

    write(*,*) 'eob_clean: reducing numObs from ', ensObs%numObs, ' to ', numObsClean
    call eob_allocate(ensObsClean, ensObs%numMembers, numObsClean, ensObs%obsSpaceData)

    obsCleanIndex = 0
    do obsIndex = 1, ensObs%numObs
      if (ensObs%assFlag(obsIndex) == 1) then
        obsCleanIndex = obsCleanIndex + 1
        ensObsClean%lat(obsCleanIndex)        = ensObs%lat(obsIndex)
        ensObsClean%lon(obsCleanIndex)        = ensObs%lon(obsIndex)
        ensObsClean%logPres(obsCleanIndex)    = ensObs%logPres(obsIndex)
        ensObsClean%obsErrInv(obsCleanIndex)  = ensObs%obsErrInv(obsIndex)
        ensObsClean%Yb_r4(:,obsCleanIndex)    = ensObs%Yb_r4(:,obsIndex)
        ensObsClean%randPert_r4(:,obsCleanIndex) = ensObs%randPert_r4(:,obsIndex)
        ensObsClean%meanYb(obsCleanIndex)     = ensObs%meanYb(obsIndex)
        ensObsClean%deterYb(obsCleanIndex)    = ensObs%deterYb(obsIndex)
        ensObsClean%obsValue(obsCleanIndex)   = ensObs%obsValue(obsIndex)
        ensObsClean%assFlag(obsCleanIndex)    = ensObs%assFlag(obsIndex)
      end if
    end do

    ! reallocate the original object with the new size
    call eob_deallocate(ensObs)
    call eob_allocate(ensObs, ensObsClean%numMembers, numObsClean, ensObsClean%obsSpaceData)

    ! copy the cleaned object into the original object and deallocate local copy
    call eob_copy(ensObsClean, ensObs)
    call eob_deallocate(ensObsClean)

  end subroutine eob_clean

  !--------------------------------------------------------------------------
  ! eob_copy
  !--------------------------------------------------------------------------
  subroutine eob_copy(ensObsIn,ensObsOut)
    implicit none

    ! arguments
    type(struct_eob) :: ensObsIn
    type(struct_eob) :: ensObsOut

    ensObsOut%lat(:)        = ensObsIn%lat(:)
    ensObsOut%lon(:)        = ensObsIn%lon(:)
    ensObsOut%logPres(:)    = ensObsIn%logPres(:)
    ensObsOut%obsErrInv(:)  = ensObsIn%obsErrInv(:)
    ensObsOut%Yb_r4(:,:)    = ensObsIn%Yb_r4(:,:)
    ensObsOut%randPert_r4(:,:) = ensObsIn%randPert_r4(:,:)
    ensObsOut%meanYb(:)     = ensObsIn%meanYb(:)
    ensObsOut%deterYb(:)    = ensObsIn%deterYb(:)
    ensObsOut%obsValue(:)   = ensObsIn%obsValue(:)
    ensObsOut%assFlag(:)    = ensObsIn%assFlag(:)

  end subroutine eob_copy

  !--------------------------------------------------------------------------
  ! eob_allGather
  !--------------------------------------------------------------------------
  subroutine eob_allGather(ensObs,ensObs_mpiglobal)
    !
    ! :Purpose: Collect obs information distributed over all mpi tasks and
    !           make it available on all mpi tasks. The output ensObs object
    !           will be allocated within this subroutine.
    !
    implicit none

    ! arguments
    type(struct_eob) :: ensObs
    type(struct_eob) :: ensObs_mpiglobal

    ! locals
    integer :: ierr, procIndex, memberIndex, numObs_mpiglobal
    integer :: allNumObs(mpi_nprocs), displs(mpi_nprocs)

    write(*,*) 'eob_allGather: starting'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! refresh assimilation flag and then clean ensObs before communicating and writing
    call eob_setAssFlag(ensObs)
    call eob_clean(ensObs)

    call rpn_comm_allgather( ensObs%numObs, 1, 'mpi_integer',  &
                             allNumObs, 1, 'mpi_integer', &
                             'GRID', ierr )
    numObs_mpiglobal = sum(allNumObs(:))

    if (ensObs_mpiglobal%allocated) then
      call utl_abort('eob_allGather: output ensObs object must not be already allocated')
    end if
    call eob_allocate(ensObs_mpiglobal, ensObs%numMembers, numObs_mpiglobal, ensObs%obsSpaceData)

    if ( mpi_myid == 0 ) then
      displs(1) = 0
      do procIndex = 2, mpi_nprocs
        displs(procIndex) = displs(procIndex-1) + allNumObs(procIndex-1)
      end do
    else
      displs(:) = 0
    end if

    call rpn_comm_gatherv( ensObs%lat          , ensObs%numObs, 'mpi_real8', &
                           ensObs_mpiglobal%lat, allNumObs, displs, 'mpi_real8',  &
                           0, 'GRID', ierr )
    call rpn_comm_gatherv( ensObs%lon          , ensObs%numObs, 'mpi_real8', &
                           ensObs_mpiglobal%lon, allNumObs, displs, 'mpi_real8',  &
                           0, 'GRID', ierr )
    call rpn_comm_gatherv( ensObs%logPres          , ensObs%numObs, 'mpi_real8', &
                           ensObs_mpiglobal%logPres, allNumObs, displs, 'mpi_real8',  &
                           0, 'GRID', ierr )
    call rpn_comm_gatherv( ensObs%obsValue          , ensObs%numObs, 'mpi_real8', &
                           ensObs_mpiglobal%obsValue, allNumObs, displs, 'mpi_real8',  &
                           0, 'GRID', ierr )
    call rpn_comm_gatherv( ensObs%obsErrInv          , ensObs%numObs, 'mpi_real8', &
                           ensObs_mpiglobal%obsErrInv, allNumObs, displs, 'mpi_real8',  &
                           0, 'GRID', ierr )
    call rpn_comm_gatherv( ensObs%meanYb          , ensObs%numObs, 'mpi_real8', &
                           ensObs_mpiglobal%meanYb, allNumObs, displs, 'mpi_real8',  &
                           0, 'GRID', ierr )
    call rpn_comm_gatherv( ensObs%deterYb          , ensObs%numObs, 'mpi_real8', &
                           ensObs_mpiglobal%deterYb, allNumObs, displs, 'mpi_real8',  &
                           0, 'GRID', ierr )
    call rpn_comm_gatherv( ensObs%assFlag          , ensObs%numObs, 'mpi_integer', &
                           ensObs_mpiglobal%assFlag, allNumObs, displs, 'mpi_integer',  &
                           0, 'GRID', ierr )
    do memberIndex = 1, ensObs%numMembers
      call rpn_comm_gatherv( ensObs%Yb_r4(memberIndex,:)          , ensObs%numObs, 'mpi_real4', &
                             ensObs_mpiglobal%Yb_r4(memberIndex,:), allNumObs, displs, 'mpi_real4',  &
                             0, 'GRID', ierr )
    end do
    do memberIndex = 1, ensObs%numMembers
      call rpn_comm_gatherv( ensObs%randPert_r4(memberIndex,:)          , ensObs%numObs, 'mpi_real4', &
                             ensObs_mpiglobal%randPert_r4(memberIndex,:), allNumObs, displs, 'mpi_real4',  &
                             0, 'GRID', ierr )
    end do

    call rpn_comm_bcast(ensObs_mpiglobal%lat, ensObs_mpiglobal%numObs, 'mpi_real8',  &
                        0, 'GRID', ierr)
    call rpn_comm_bcast(ensObs_mpiglobal%lon, ensObs_mpiglobal%numObs, 'mpi_real8',  &
                        0, 'GRID', ierr)
    call rpn_comm_bcast(ensObs_mpiglobal%logPres, ensObs_mpiglobal%numObs, 'mpi_real8',  &
                        0, 'GRID', ierr)
    call rpn_comm_bcast(ensObs_mpiglobal%obsValue, ensObs_mpiglobal%numObs, 'mpi_real8',  &
                        0, 'GRID', ierr)
    call rpn_comm_bcast(ensObs_mpiglobal%obsErrInv, ensObs_mpiglobal%numObs, 'mpi_real8',  &
                        0, 'GRID', ierr)
    call rpn_comm_bcast(ensObs_mpiglobal%meanYb, ensObs_mpiglobal%numObs, 'mpi_real8',  &
                        0, 'GRID', ierr)
    call rpn_comm_bcast(ensObs_mpiglobal%deterYb, ensObs_mpiglobal%numObs, 'mpi_real8',  &
                        0, 'GRID', ierr)
    call rpn_comm_bcast(ensObs_mpiglobal%assFlag, ensObs_mpiglobal%numObs, 'mpi_integer',  &
                        0, 'GRID', ierr)
    do memberIndex = 1, ensObs%numMembers
      call rpn_comm_bcast(ensObs_mpiglobal%Yb_r4(memberIndex,:), ensObs_mpiglobal%numObs, 'mpi_real4',  &
                          0, 'GRID', ierr)
    end do
    do memberIndex = 1, ensObs%numMembers
      call rpn_comm_bcast(ensObs_mpiglobal%randPert_r4(memberIndex,:), ensObs_mpiglobal%numObs, 'mpi_real4',  &
                          0, 'GRID', ierr)
    end do

    write(*,*) 'eob_allGather: total number of obs to be assimilated =', sum(ensObs_mpiglobal%assFlag(:))

    write(*,*) 'eob_allGather: finished'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine eob_allGather

  !--------------------------------------------------------------------------
  ! eob_writeToFiles
  !--------------------------------------------------------------------------
  subroutine eob_writeToFiles(ensObs)
    !
    ! :Purpose: Write the contents of an ensObs object to files
    !
    implicit none

    ! arguments
    type(struct_eob)         :: ensObs

    ! locals
    integer :: unitNum, ierr, obsIndex, memberIndex
    character(len=40) :: fileName
    character(len=4)  :: memberIndexStr
    integer :: fnom, fclos

    ! only the first mpi task does writing, assuming mpi gather already done
    if (mpi_myid /= 0) return

    if ( .not.ensObs%allocated ) then
      call utl_abort('eob_writeToFiles: this object is not allocated')
    end if

    ! write the lat, lon and obs values to a file
    fileName = 'eob_Lat_Lon_ObsValue'
    write(*,*) 'eob_writeToFiles: writing ',trim(filename)
    unitNum = 0
    ierr = fnom(unitNum, fileName, 'FTN+SEQ+UNF+R/W', 0)
    write(unitNum) (ensObs%lat(obsIndex), obsIndex = 1, ensObs%numObs)
    write(unitNum) (ensObs%lon(obsIndex), obsIndex = 1, ensObs%numObs)
    write(unitNum) (ensObs%obsValue(obsIndex), obsIndex = 1, ensObs%numObs)
    ierr = fclos(unitNum)

    ! write the contents of Yb, 1 member per file
    do memberIndex = 1, ensObs%numMembers
      write(memberIndexStr,'(I0.4)') memberIndex
      fileName = 'eob_HX_' // memberIndexStr
      write(*,*) 'eob_writeToFiles: writing ',trim(filename)
      unitNum = 0
      ierr = fnom(unitNum, fileName, 'FTN+SEQ+UNF+R/W', 0)
      write(unitNum) (ensObs%Yb_r4(memberIndex,obsIndex), obsIndex = 1, ensObs%numObs)
      ierr = fclos(unitNum)
    end do

  end subroutine eob_writeToFiles

  !--------------------------------------------------------------------------
  ! eob_readFromFiles
  !--------------------------------------------------------------------------
  subroutine eob_readFromFiles(ensObs)
    !
    ! :Purpose: Read the contents of an ensObs object from files
    !
    implicit none

    ! arguments
    type(struct_eob)         :: ensObs

    call utl_abort('eob_readFromFiles: not yet implemented')

  end subroutine eob_readFromFiles

  !--------------------------------------------------------------------------
  ! eob_getLocalBodyIndices
  !--------------------------------------------------------------------------
  function eob_getLocalBodyIndices(ensObs,localBodyIndices,distances,lat,lon,logPres,  &
                                   hLocalize,vLocalize,numLocalObsFound) result(numLocalObs)
    !
    ! :Purpose: Return a list of values of bodyIndex for all observations within 
    !           the local volume around the specified lat/lon used for assimilation
    !           (as defined by h/vLocalize). The kdtree2 module is used to efficiently
    !           perform this task. The kdtree itself is constructed on the first call.
    !
    implicit none

    ! arguments
    integer          :: numLocalObs
    type(struct_eob) :: ensObs
    integer          :: localBodyIndices(:)
    real(8)          :: distances(:)
    real(8)          :: lat, lon, logPres, hLocalize, vLocalize
    integer          :: numLocalObsFound

    ! locals
    integer :: bodyIndex, numLocalObsFoundSearch, maxNumLocalObs, localObsIndex
    real(8) :: distance
    real(kdkind), allocatable         :: positionArray(:,:)
    type(kdtree2_result), allocatable :: searchResults(:)
    real(kdkind)                      :: maxRadius
    real(kdkind)                      :: refPosition(3)

    ! create the kdtree on the first call
    if (.not. associated(tree)) then
      write(*,*) 'eob_getLocalBodyIndices: start creating kdtree'
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
      allocate(positionArray(3,ensObs%numObs))
      do bodyIndex = 1, ensObs%numObs
        positionArray(1,bodyIndex) = RA * sin(ensObs%lon(bodyIndex)) * cos(ensObs%lat(bodyIndex))
        positionArray(2,bodyIndex) = RA * cos(ensObs%lon(bodyIndex)) * cos(ensObs%lat(bodyIndex))
        positionArray(3,bodyIndex) = RA *                              sin(ensObs%lat(bodyIndex))
      end do
      tree => kdtree2_create(positionArray, sort=.true., rearrange=.true.) 
      write(*,*) 'eob_getLocalBodyIndices: done creating kdtree'
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    end if

    ! do the search
    maxNumLocalObs = size(localBodyIndices)
    maxRadius = hLocalize**2
    refPosition(1) = RA * sin(lon) * cos(lat)
    refPosition(2) = RA * cos(lon) * cos(lat)
    refPosition(3) = RA *            sin(lat)
    allocate(searchResults(maxNumLocalObsSearch))
    call kdtree2_r_nearest(tp=tree, qv=refPosition, r2=maxRadius, nfound=numLocalObsFoundSearch,&
                           nalloc=maxNumLocalObsSearch, results=searchResults)
    if (numLocalObsFoundSearch > maxNumLocalObsSearch) then
      call utl_abort('eob_getLocalBodyIndices: the parameter maxNumLocalObsSearch must be increased')
    end if

    ! copy search results to output vectors, only those within vertical localization distance
    numLocalObsFound = 0
    numLocalObs = 0
    do localObsIndex=1, numLocalObsFoundSearch
      distance = abs( logPres - ensObs%logPres(searchResults(localObsIndex)%idx) )
      if (distance <= vLocalize .and. ensObs%assFlag(searchResults(localObsIndex)%idx)==1) then
        numLocalObsFound = numLocalObsFound + 1
        if (numLocalObs < maxNumLocalObs) then
          numLocalObs = numLocalObs + 1
          localBodyIndices(numLocalObs) = searchResults(localObsIndex)%idx
          distances(numLocalObs) = sqrt(searchResults(localObsIndex)%dis)
        end if
      end if
    end do
    deallocate(searchResults)

  end function eob_getLocalBodyIndices

  !--------------------------------------------------------------------------
  ! eob_setLatLonObs
  !--------------------------------------------------------------------------
  subroutine eob_setLatLonObs(ensObs)
    implicit none

    ! arguments
    type(struct_eob) :: ensObs

    call obs_extractObsRealHeaderColumn(ensObs%lat, ensObs%obsSpaceData, OBS_LAT)
    call obs_extractObsRealHeaderColumn(ensObs%lon, ensObs%obsSpaceData, OBS_LON)
    call obs_extractObsRealBodyColumn(ensObs%obsValue, ensObs%obsSpaceData, OBS_VAR)

  end subroutine eob_setLatLonObs

  !--------------------------------------------------------------------------
  ! eob_setobsErrInv
  !--------------------------------------------------------------------------
  subroutine eob_setobsErrInv(ensObs)
    implicit none

    ! arguments
    type(struct_eob) :: ensObs

    ! locals
    integer :: obsIndex

    call obs_extractObsRealBodyColumn(ensObs%obsErrInv, ensObs%obsSpaceData, OBS_OER)
    do obsIndex = 1, ensObs%numObs
      if(ensObs%obsErrInv(obsIndex) > 0.0d0) then
        ensObs%obsErrInv(obsIndex) = 1.0d0/(ensObs%obsErrInv(obsIndex)**2)
      else
        ensObs%obsErrInv(obsIndex) = 0.0d0
      end if
    end do

  end subroutine eob_setobsErrInv

  !--------------------------------------------------------------------------
  ! eob_setLogPres
  !--------------------------------------------------------------------------
  subroutine eob_setLogPres(ensObs, columnMeanTrl)
    !
    ! :Purpose: Set the ln(pressure) value for each observation that 
    !           will be used when doing vertical localization. For
    !           radiance observations, the level of the maximum value
    !           of the derivative of transmission is used.
    !
    implicit none

    ! arguments
    type(struct_eob)        :: ensObs
    type(struct_columnData) :: columnMeanTrl

    ! locals
    integer          :: obsIndex, headerIndex, channelIndex, tovsIndex, numTovsLevels, nosensor
    integer          :: levIndex, levIndexBelow, levIndexAbove, nLev_M
    integer          :: varNumber(ensObs%numObs), obsVcoCode(ensObs%numObs), codType(ensObs%numObs)
    real(8)          :: obsHeight, interpFactor, obsPPP(ensObs%numObs)
    real(8), pointer :: sfcPres_ptr(:,:), presM_ptr(:,:), heightM_ptr(:,:)

    call eob_setAssFlag(ensObs)

    presM_ptr   => col_getAllColumns(columnMeanTrl,'P_M')
    heightM_ptr => col_getAllColumns(columnMeanTrl,'Z_M')
    sfcPres_ptr => col_getAllColumns(columnMeanTrl,'P0')
    nLev_M = col_getNumLev(columnMeanTrl,'MM')

    call obs_extractObsRealBodyColumn(obsPPP, ensObs%obsSpaceData, OBS_PPP) ! this needs to work also for non-pressure level obs!!!
    call obs_extractObsIntBodyColumn(varNumber, ensObs%obsSpaceData, OBS_VNM)
    call obs_extractObsIntBodyColumn(obsVcoCode, ensObs%obsSpaceData, OBS_VCO)
    call obs_extractObsIntHeaderColumn(codType, ensObs%obsSpaceData, OBS_ITY)
    do obsIndex = 1, ensObs%numObs
      headerIndex = obs_bodyElem_i(ensObs%obsSpaceData,OBS_HIND,obsIndex)

      if( varNumber(obsIndex) == BUFR_NETS .or. varNumber(obsIndex) == BUFR_NEPS   .or.  &
          varNumber(obsIndex) == BUFR_NEUS .or. varNumber(obsIndex) == BUFR_NEVS   .or.  &
          varNumber(obsIndex) == BUFR_NESS .or. varNumber(obsIndex) == BUFR_NEPN   .or.  &
          varNumber(obsIndex) == BUFR_VIS  .or. varNumber(obsIndex) == BUFR_LOGVIS .or.  &
          varNumber(obsIndex) == BUFR_GUST .or.  &
          varNumber(obsIndex) == BUFR_radarPrecip .or. varNumber(obsIndex) == BUFR_logRadarPrecip ) then

        ! all surface observations
        ensObs%logPres(obsIndex) = log(sfcPres_ptr(1,headerIndex))

      else if (varNumber(obsIndex) == BUFR_NEZD) then

        ! ZTD observation, try 0.7*Psfc (i.e. ~700hPa when Psfc=1000hPa)
        ensObs%logPres(obsIndex) = log(0.7D0 * sfcPres_ptr(1,headerIndex))

      else if (obsPPP(obsIndex) > 0.0d0 .and. obsVcoCode(obsIndex)==2) then

        ! all pressure level observations
        ensObs%logPres(obsIndex) = log(obsPPP(obsIndex))

      else if(obsVcoCode(obsIndex) == 1) then

        ! all height level observations (not including surface obs)
        obsHeight = obsPPP(obsIndex)

        ! find level just below the observation
        levIndexBelow = 0
        LEV_LOOP: do levIndex = 1, nLev_M
          if (obsHeight > heightM_ptr(levIndex,headerIndex)) then
            levIndexBelow = levIndex
            exit LEV_LOOP
          end if
        end do LEV_LOOP

        ! set the log pressure for observation
        if (levIndexBelow == 1) then
          ! above top level, use top level pressure
          ensObs%logPres(obsIndex) = log(presM_ptr(1,headerIndex))
        else if (levIndexBelow == 0) then
          ! below bottom level, use surface pressure
          ensObs%logPres(obsIndex) = log(sfcPres_ptr(1,headerIndex))
        else
          ! interpolate
          levIndexAbove = levIndexBelow - 1
          interpFactor = ( obsHeight                              - heightM_ptr(levIndexBelow,headerIndex) ) /  &
                         ( heightM_ptr(levIndexAbove,headerIndex) - heightM_ptr(levIndexBelow,headerIndex) )
          ensObs%logPres(obsIndex) = interpFactor           * log(presM_ptr(levIndexAbove,headerIndex)) +  &
                                     (1.0D0 - interpFactor) * log(presM_ptr(levIndexBelow,headerIndex))
        end if

      else if(tvs_isIdBurpTovs(codType(obsIndex))) then

        tovsIndex = tvs_tovsIndex(headerIndex)
        nosensor = tvs_lsensor(tovsIndex)
        numTovsLevels   = size(tvs_transmission(tovsIndex)%tau_levels,1)
        channelIndex = nint(obsPPP(obsIndex))
        channelIndex = max(0,min(channelIndex,tvs_maxChannelNumber+1))
        channelIndex = channelIndex - tvs_channelOffset(nosensor)
        channelIndex = utl_findArrayIndex(tvs_ichan(:,nosensor), tvs_nchan(nosensor), channelIndex)
        if (channelIndex > 0 .and. ensObs%assFlag(obsIndex)==1) then
          call max_transmission(tvs_transmission(tovsIndex), numTovsLevels, &
                                channelIndex, tvs_profiles(tovsIndex)%p, ensObs%logPres(obsIndex))
          if(mpi_myid == 0) then
            write(*,*) 'eob_setLogPres for tovs: ', codType(obsIndex), obsPPP(obsIndex), 0.01*exp(ensObs%logPres(obsIndex))
          end if
        else
          ensObs%logPres(obsIndex) = log(500.0D2)
        end if

      else if(ensObs%assFlag(obsIndex)==1) then

        write(*,*) 'eob_setLatLonPresObs: ERROR! cannot compute pressure for this observation: ',  &
                   obsPPP(obsIndex), varNumber(obsIndex), obsVcoCode(obsIndex)
        call utl_abort('eob_setLogPres')

      end if
    end do

  end subroutine eob_setLogPres

  !--------------------------------------------------------------------------
  ! eob_setAssFlag
  !--------------------------------------------------------------------------
  subroutine eob_setAssFlag(ensObs)
    implicit none

    type(struct_eob) :: ensObs

    call obs_extractObsIntBodyColumn(ensObs%assFlag, ensObs%obsSpaceData, OBS_ASS)

  end subroutine eob_setAssFlag

  !--------------------------------------------------------------------------
  ! eob_setYb
  !--------------------------------------------------------------------------
  subroutine eob_setYb(ensObs, memberIndex)
    implicit none

    type(struct_eob) :: ensObs
    integer          :: memberIndex

    ! get the Y-HX value from obsSpaceData
    call obs_extractObsRealBodyColumn_r4(ensObs%Yb_r4(memberIndex,:), ensObs%obsSpaceData, OBS_OMP)

    ! now compute HX = Y - (Y-HX)
    ensObs%Yb_r4(memberIndex,:) = ensObs%obsValue(:) - ensObs%Yb_r4(memberIndex,:)

  end subroutine eob_setYb

  !--------------------------------------------------------------------------
  ! eob_setDeterYb
  !--------------------------------------------------------------------------
  subroutine eob_setDeterYb(ensObs)
    implicit none

    type(struct_eob) :: ensObs

    ! get the Y-HX value from obsSpaceData
    call obs_extractObsRealBodyColumn(ensObs%DeterYb(:), ensObs%obsSpaceData, OBS_OMP)

    ! now compute HX = Y - (Y-HX)
    ensObs%DeterYb(:) = ensObs%obsValue(:) - ensObs%DeterYb(:)

  end subroutine eob_setDeterYb

  !--------------------------------------------------------------------------
  ! eob_calcAndRemoveMeanYb
  !--------------------------------------------------------------------------
  subroutine eob_calcAndRemoveMeanYb(ensObs)
    implicit none

    ! arguments
    type(struct_eob) :: ensObs

    ! locals
    integer :: obsIndex

    do obsIndex = 1, ensObs%numObs
      ensObs%meanYb(obsIndex) = sum(ensObs%Yb_r4(:,obsIndex)) / ensObs%numMembers
      ensObs%Yb_r4(:,obsIndex) = ensObs%Yb_r4(:,obsIndex) - ensObs%meanYb(obsIndex)
    end do

  end subroutine eob_calcAndRemoveMeanYb

  !--------------------------------------------------------------------------
  ! eob_calcRandPert
  !--------------------------------------------------------------------------
  subroutine eob_calcRandPert(ensObs, randomSeed)
    implicit none

    ! arguments
    type(struct_eob) :: ensObs
    integer          :: randomSeed

    ! locals
    integer :: obsIndex, memberIndex
    real(4) :: meanRandPert, sigObs

    call rng_setup(abs(randomSeed))

    do obsIndex = 1, ensObs%numObs
      sigObs = obs_bodyElem_r(ensObs%obsSpaceData, OBS_OER, obsIndex)
      do memberIndex = 1, ensObs%numMembers
        ensObs%randPert_r4(memberIndex,obsIndex) = sigObs * rng_gaussian()
      end do

      meanRandPert = sum(ensObs%randPert_r4(:,obsIndex)) / real(ensObs%numMembers,4)
      ensObs%randPert_r4(:,obsIndex) = ensObs%randPert_r4(:,obsIndex) - meanRandPert
    end do

  end subroutine eob_calcRandPert

  !--------------------------------------------------------------------------
  ! eob_setMeanOMP
  !--------------------------------------------------------------------------
  subroutine eob_setMeanOMP(ensObs)
    implicit none

    ! arguments
    type(struct_eob) :: ensObs

    ! locals
    integer :: obsIndex

    do obsIndex = 1, ensObs%numObs
      call obs_bodySet_r(ensObs%obsSpaceData, OBS_OMP, obsIndex,  &
                         ensObs%obsValue(obsIndex)-ensObs%meanYb(obsIndex))
    end do

  end subroutine eob_setMeanOMP

  !--------------------------------------------------------------------------
  ! eob_setHPHT
  !--------------------------------------------------------------------------
  subroutine eob_setHPHT(ensObs)
    implicit none

    ! arguments
    type(struct_eob) :: ensObs

    ! locals
    integer :: obsIndex, memberIndex
    real(8) :: hpht, sigo, sigi

    do obsIndex = 1, ensObs%numObs
      hpht = 0.0d0
      do memberIndex = 1, ensObs%numMembers
        hpht = hpht + ensObs%Yb_r4(memberIndex,obsIndex)**2 / ensObs%numMembers
      end do
      if (hpht > 0.0D0) then 
        hpht = sqrt(hpht)
      else
        hpht = 0.0D0
      end if
      call obs_bodySet_r(ensObs%obsSpaceData, OBS_HPHT, obsIndex, hpht)

      ! also set 'sigi' and 'sigo'
      sigo = obs_bodyElem_r(ensObs%obsSpaceData, OBS_OER, obsIndex)
      sigi = (sigo**2 + hpht**2)**0.5
      call obs_bodySet_r(ensObs%obsSpaceData, OBS_SIGI, obsIndex, sigi)
      call obs_bodySet_r(ensObs%obsSpaceData, OBS_SIGO, obsIndex, sigo)
    end do

  end subroutine eob_setHPHT

  !--------------------------------------------------------------------------
  ! eob_backgroundCheck
  !--------------------------------------------------------------------------
  subroutine eob_backgroundCheck(ensObs)
    !
    ! :Purpose: Apply additional background using the ensemble spread.
    !
    implicit none

    ! arguments:
    type(struct_eob) :: ensObs

    ! locals:
    integer :: bodyIndexBeg, bodyIndexEnd, headerIndex, ivar, bodyIndex
    integer :: numRejected, numRejectedMpiGlobal, ierr, windCount
    real    :: sigo, sigb, omp, sig, reject_limit
    logical :: reject_wind

    numRejected = 0
    reject_limit = 5.0
    do headerIndex = 1, obs_numheader(ensObs%obsSpaceData)
      bodyIndexBeg = obs_headElem_i(ensObs%obsSpaceData, OBS_RLN, headerIndex)
      bodyIndexEnd = obs_headElem_i(ensObs%obsSpaceData, OBS_NLV, headerIndex) + bodyIndexBeg - 1
      reject_wind = .false.

      do bodyIndex = bodyIndexBeg, bodyIndexEnd
        if (obs_bodyElem_i(ensObs%obsSpaceData, OBS_ASS, bodyIndex) == obs_notAssimilated) cycle
        sigo = obs_bodyElem_r(ensObs%obsSpaceData, OBS_OER, bodyIndex)
        sigb = obs_bodyElem_r(ensObs%obsSpaceData, OBS_HPHT, bodyIndex)
        ! cut off at reject_limit standard deviations
        sig = reject_limit*(sigo**2 + sigb**2)**0.5
        omp = abs(obs_bodyElem_r(ensObs%obsSpaceData, OBS_OMP, bodyIndex))
        if (omp > sig) then
          call obs_bodySet_i(ensObs%obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
          call obs_bodySet_i(ensObs%obsSpaceData, OBS_FLG, bodyIndex,  &
                             IBSET(obs_bodyElem_i(ensObs%obsSpaceData, OBS_FLG, bodyIndex),9))
          ivar = obs_bodyElem_i(ensObs%obsSpaceData, OBS_VNM, bodyIndex)
          reject_wind = bufr_isWindComponent(ivar)
          write(*,*) 'eob_backgroundCheck: headerIndex, omp, sig, ivar, reject_wind', &
                                           headerIndex, omp, sig, ivar, reject_wind
          numRejected = numRejected + 1
        end if
      end do

      ! take the same reject decision for both components of the
      ! wind vector (and any other winds for that station).
      ! N.B.: This seems to assume only one level per header, this is generally 
      ! the case for wind observations currently, since radiosondes are 4D
      if (reject_wind) then 
        ! first count how many wind observations we have for this station
        windCount = 0
        do bodyIndex = bodyIndexBeg, bodyIndexEnd
          if (bufr_isWindComponent(obs_bodyElem_i(ensObs%obsSpaceData, OBS_VNM, bodyIndex))) then
            windCount = windCount + 1
          end if
        end do
        if (windCount > 2) then
          write(*,*) 'eob_backgroundCheck: WARNING' 
          write(*,*) 'Station ',headerIndex,' has ',windCount,' wind observations '
          write(*,*) 'Perhaps old radiosonde format - std dev not changed for other wind component'
        else
          do bodyIndex = bodyIndexBeg, bodyIndexEnd
            if (obs_bodyElem_i(ensObs%obsSpaceData, OBS_ASS, bodyIndex) == obs_notassimilated) cycle
            ivar = obs_bodyElem_i(ensObs%obsSpaceData, OBS_VNM, bodyIndex)
            if (bufr_isWindComponent(ivar) .and.  &
                obs_bodyElem_i(ensObs%obsSpaceData, OBS_ASS, bodyIndex) == obs_assimilated) then
              call obs_bodySet_i(ensObs%obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
              call obs_bodySet_i(ensObs%obsSpaceData, OBS_FLG, bodyIndex,  &
                                 IBSET(obs_bodyElem_i(ensObs%obsSpaceData, OBS_FLG, bodyIndex),9))
              write(*,*) 'eob_backgroundCheck: other wind component headerIndex, ivar', &
                                               headerIndex, ivar
              numRejected = numRejected + 1
            end if
          end do
        end if
      end if
    end do

    call rpn_comm_allreduce(numRejected,numRejectedMpiGlobal,1,'mpi_integer','mpi_sum','GRID',ierr)
    write(*,*)
    write(*,*) 'eob_backgroundCheck: number of observations rejected (local) =', numRejected
    write(*,*) 'eob_backgroundCheck: number of observations rejected (global)=', numRejectedMpiGlobal

  end subroutine eob_backgroundCheck

  !--------------------------------------------------------------------------
  ! eob_huberNorm
  !--------------------------------------------------------------------------
  subroutine eob_huberNorm(ensObs)
    !
    ! :Purpose: Apply huber norm quality control procedure.
    !
    implicit none

    ! arguments
    type(struct_eob) :: ensObs

    ! locals
    integer        :: huberCount, huberCountMpiGlobal, ivar, windCount, ierr
    integer        :: bodyIndex, bodyIndexBeg, bodyIndexEnd, headerIndex
    real(obs_real) :: c_limit, sig, sigo, sigb, omp, sigo_hub, sigo_hub_wind
    logical        :: reject_wind

    c_limit = 2.0
    huberCount = 0
    do headerIndex = 1,obs_numheader(ensObs%obsSpaceData)
      bodyIndexBeg = obs_headElem_i(ensObs%obsSpaceData, OBS_RLN, headerIndex)
      bodyIndexEnd = obs_headElem_i(ensObs%obsSpaceData, OBS_NLV, headerIndex) + bodyIndexBeg - 1
      reject_wind = .false.
      sigo_hub_wind = 0.0
      do bodyIndex = bodyIndexBeg, bodyIndexEnd
        sigo = obs_bodyElem_r(ensObs%obsSpaceData, OBS_OER, bodyIndex)
        sigb = obs_bodyElem_r(ensObs%obsSpaceData, OBS_HPHT, bodyIndex)
        ! cut off at reject_limit standard deviations
        sig = c_limit*(sigo**2 + sigb**2)**0.5
        omp = abs(obs_bodyElem_r(ensObs%obsSpaceData, OBS_OMP, bodyIndex))
        if (omp > sig) then
          ! redefining the observational error such that the innovation
          ! is at exactly c_limit standard deviations.
          huberCount = huberCount + 1
          sigo_hub = ((omp/c_limit)**2-sigb**2)**0.5
          ivar = obs_bodyElem_i(ensObs%obsSpaceData, OBS_VNM, bodyIndex)
          if (bufr_isWindComponent(ivar)) then
            ! the other wind components will be changed at the same time (next if)
            reject_wind = .true.
            call obs_bodySet_r(ensObs%obsSpaceData, OBS_OER, bodyIndex, sigo_hub)
            ! this is for the special case both components of the wind innovation
            ! suggest using the Huber norm. 
            if (sigo_hub > sigo_hub_wind) then
              sigo_hub_wind = sigo_hub
            end if
          else
            call obs_bodySet_r(ensObs%obsSpaceData, OBS_OER, bodyIndex, sigo_hub)
          end if
        end if
      end do

      ! Give the same inflated observation error to both wind components.
      ! this part assumes that modern sondes are used (at most two
      ! wind observations per station.
      if (reject_wind) then
        ! first count how many wind observations we have for this station
        windCount = 0
        do bodyIndex = bodyIndexBeg, bodyIndexEnd
          if (bufr_isWindComponent(obs_bodyElem_i(ensObs%obsSpaceData, OBS_VNM, bodyIndex))) then
            windCount = windCount + 1
          end if
        end do
        if (windCount > 2) then
          write(*,*) 'Warning Hubernorm' 
          write(*,*) 'Station ',headerIndex,' has ',windCount,' wind observations '
          write(*,*) 'Perhaps old radiosonde format - std dev not changed for other wind component'
        else
          huberCount = huberCount + 1
          do bodyIndex = bodyIndexBeg, bodyIndexEnd
            if (bufr_isWindComponent(obs_bodyElem_i(ensObs%obsSpaceData, OBS_VNM, bodyIndex))) then
              call obs_bodySet_r(ensObs%obsSpaceData, OBS_OER, bodyIndex, sigo_hub_wind)
            end if
          end do
        end if
      end if 
    end do

    call rpn_comm_allreduce(huberCount, huberCountMpiGlobal, 1, 'mpi_integer', 'mpi_sum', 'GRID', ierr)
    write(*,*)
    write(*,*) 'eob_huberNorm: number of obs with increased error stddev (local) = ', huberCount
    write(*,*) 'eob_huberNorm: number of obs with increased error stddev (global)= ', huberCountMpiGlobal

  end subroutine eob_huberNorm

  !--------------------------------------------------------------------------
  ! eob_rejectRadNearSfc
  !--------------------------------------------------------------------------
  subroutine eob_rejectRadNearSfc(ensObs)
    !
    ! :Purpose: Reject all radiance observations with peak sensitivity
    !           too close to the surface.
    !
    implicit none

    ! arguments
    type(struct_eob) :: ensObs

    ! locals
    integer :: acceptCount, rejectCount, acceptCountMpiGlobal, rejectCountMpiGlobal
    integer :: varNumber, bodyIndex, ierr
    ! reject lower than 975 hPa
    real, parameter    :: logPresRadianceLimit = 11.4876E0

    acceptCount = 0
    rejectCount = 0
    do bodyIndex = 1, obs_numbody(ensObs%obsSpaceData)
      varNumber = obs_bodyElem_i(ensObs%obsSpaceData, OBS_VNM, bodyIndex)
      if (varNumber == bufr_nbt1 .or.  &
          varNumber == bufr_nbt2 .or.  &
          varNumber == bufr_nbt3) then
        if (ensObs%logPres(bodyIndex) < logPresRadianceLimit) then
          acceptCount = acceptCount + 1
        else
          rejectCount = rejectCount + 1
          call obs_bodySet_i(ensObs%obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
          call obs_bodySet_i(ensObs%obsSpaceData, OBS_FLG, bodyIndex,  &
                             IBSET(obs_bodyElem_i(ensObs%obsSpaceData, OBS_FLG, bodyIndex),9))
        end if 
      end if
    end do

    call rpn_comm_allreduce(acceptCount, acceptCountMpiGlobal, 1, 'mpi_integer', 'mpi_sum', 'GRID', ierr)
    call rpn_comm_allreduce(rejectCount, rejectCountMpiGlobal, 1, 'mpi_integer', 'mpi_sum', 'GRID', ierr)
    write(*,*)
    write(*,*) 'eob_rejectRadNearSfc: Number of accepted, rejected observations (local) : ',  &
               acceptCount, rejectCount
    write(*,*) 'eob_rejectRadNearSfc: Number of accepted, rejected observations (global): ',  &
               acceptCountMpiGlobal, rejectCountMpiGlobal

  end subroutine eob_rejectRadNearSfc

  !--------------------------------------------------------------------------
  ! max_transmission (private routine)
  !--------------------------------------------------------------------------
  subroutine max_transmission(transmission, numLevels, transIndex, rttovPres, maxLnP)
    !
    ! :Purpose: Determine the height in log pressure where we find the maximum 
    !           value of the first derivative of transmission with respect to 
    !           log pressure
    !
    implicit none

    ! arguments
    type(rttov_transmission), intent(in) :: transmission ! transmission (rttov type)
    integer(kind=jpim), intent(in)       :: numLevels    ! number of RTTOV levels
    integer, intent(in)                  :: transIndex   ! index of transmission%tau_levels
    real(kind=jprb), pointer, intent(in) :: rttovPres(:) ! pressure of RTTOV levels
    real(8), intent(out)                 :: maxLnP       ! log pressure of maximum

    ! locals
    integer :: levIndex
    real(8) :: lnPres(numLevels), avgPres(numLevels-1)
    real(8) :: diffTau, derivTau(numLevels), maxDeriv
    integer :: nAvgLev, low, high, maxIndex

    nAvgLev = numLevels - 1
    lnPres(:) = log(rttovPres(:)*MPC_PA_PER_MBAR_R8)
    ! calculate the first derivative of transmission with respect to log pressure
    ! and find the level index for its maximum
    maxDeriv = -0.1d0
    derivTau(1) = 0.0d0
    maxIndex = numLevels
    do levIndex = 2, numLevels
      avgPres(levIndex-1) = 0.5d0*(lnPres(levIndex)+lnPres(levIndex-1))
      diffTau = transmission%tau_levels(levIndex-1,transIndex) - transmission%tau_levels(levIndex,transIndex)
      derivTau(levIndex) = diffTau / (lnPres(levIndex)-lnPres(levIndex-1))
      if (derivTau(levIndex)>maxDeriv) then
        maxDeriv = derivTau(levIndex)
        maxIndex = levIndex
      end if
    end do

    ! get the height in log pressure for the level index (maxIndex) found above
    if (maxIndex==1) maxIndex = maxIndex + 1
    if ((maxIndex==2).or.(maxIndex==numLevels)) then
      maxLnP = avgPres(maxIndex-1)
    else
      call get_peak(maxIndex,nAvgLev,avgPres,derivTau,maxLnP)
    end if

  end subroutine max_transmission

  !--------------------------------------------------------------------------
  ! get_peak (private routine)
  !--------------------------------------------------------------------------
  subroutine get_peak(maxIndex,nlev,lnp,deriv,maxLnP)
    !
    ! :Purpose: Do quadratic interpolation to find pressure of peak transmission.
    !
    implicit none

    ! arguments
    integer, intent(in)    :: maxIndex, nlev
    real(8), intent(in)    :: lnp(nlev), deriv(nlev+1)
    real(8), intent(inout) :: maxLnP

    ! locals
    external :: dgesv
    integer, parameter :: N=3
    integer :: info
    integer, parameter :: lda=N, ldb=N, nrhs=1
    integer :: ipiv(N)
    real(8) :: A(lda,N),B(ldb,nrhs)
    integer :: index1, index2

    index2 = 0
    do index1=maxIndex-1,maxIndex+1
      index2 = index2 + 1
      A(index2,1) = lnp(index1-1)*lnp(index1-1)
      A(index2,2) = lnp(index1-1)
      A(index2,3) = 1.0d0
      B(index2,1) = deriv(index1)
    end do

    call dgesv(N,nrhs,A,lda,ipiv,B,ldb,info)

    if (info==0) then
      maxLnP = -0.5*(B(2,1)/B(1,1))
    else
      maxLnP = lnp(maxIndex-1)
    end if

  end subroutine get_peak

end module ensembleObservations_mod
