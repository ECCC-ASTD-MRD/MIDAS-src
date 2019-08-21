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
  use mathPhysConstants_mod
  use physicsFunctions_mod
  use utilities_mod
  use earthconstants_mod
  use bufr_mod
  implicit none
  save
  private

  ! public types
  public :: struct_eob

  ! public procedures
  public :: eob_allocate, eob_deallocate, eob_allGather, eob_getLocalBodyIndices
  public :: eob_setYb, eob_setDeterYb, eob_setLatLonObs, eob_setMeanOMP, eob_setAssFlag
  public :: eob_setMeanHPHT, eob_calcRemoveMeanYb, eob_setPres, eob_clean, eob_copy

  integer, parameter :: maxNumLocalObsSearch = 500000
  integer,external   :: get_max_rss

  type struct_eob
    logical                       :: allocated = .false.
    integer                       :: numMembers     ! number of ensemble members
    integer                       :: numObs         ! number of observations
    type(struct_obs), pointer     :: obsSpaceData   ! pointer to obsSpaceData object
    real(8), allocatable          :: lat(:), lon(:) ! lat/lon of observation
    real(8), allocatable          :: logPres(:)     ! ln(pres) of obs, used for localization
    real(8), allocatable          :: varObsInv(:)   ! inverse of obs error variances
    real(4), allocatable          :: Yb_r4(:,:)     ! background ensemble perturbation in obs space
    real(8), allocatable          :: meanYb(:)      ! ensemble mean background state in obs space
    real(8), allocatable          :: deterYb(:)     ! deterministic background state in obs space
    real(8), allocatable          :: obsValue(:)    ! the observed value
    integer, allocatable          :: assFlag(:)     ! assimilation flag
  end type struct_eob

  type(kdtree2), pointer :: tree => null()

CONTAINS

  !--------------------------------------------------------------------------
  ! eob_allocate
  !--------------------------------------------------------------------------
  subroutine eob_allocate(ensObs, numMembers, numObs, obsSpaceData)
    ! :Purpose: Allocate an ensObs object
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
    allocate( ensObs%varObsInv(ensObs%numObs) )
    allocate( ensObs%Yb_r4(ensObs%numMembers,ensObs%numObs) )
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
    deallocate( ensObs%varObsInv )
    deallocate( ensObs%Yb_r4 )
    deallocate( ensObs%meanYb )
    deallocate( ensObs%deterYb )
    deallocate( ensObs%assFlag )

    ensObs%allocated = .false.

  end subroutine eob_deallocate

  !--------------------------------------------------------------------------
  ! eob_clean
  !--------------------------------------------------------------------------
  subroutine eob_clean(ensObs)
    ! :Purpose: Remove all obs from the ensObs object that are not 
    !           flagged for assimilation. All arrays will be reallocated
    !           in place to the smaller size after cleaning.
    implicit none

    ! arguments
    type(struct_eob) :: ensObs

    ! locals
    integer :: obsIndex, obsCleanIndex, numObsClean
    type(struct_eob) :: ensObsClean

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
        ensObsClean%lat(obsCleanIndex)       = ensObs%lat(obsIndex)
        ensObsClean%lon(obsCleanIndex)       = ensObs%lon(obsIndex)
        ensObsClean%logPres(obsCleanIndex)   = ensObs%logPres(obsIndex)
        ensObsClean%varObsInv(obsCleanIndex) = ensObs%varObsInv(obsIndex)
        ensObsClean%Yb_r4(:,obsCleanIndex)   = ensObs%Yb_r4(:,obsIndex)
        ensObsClean%meanYb(obsCleanIndex)    = ensObs%meanYb(obsIndex)
        ensObsClean%deterYb(obsCleanIndex)   = ensObs%deterYb(obsIndex)
        ensObsClean%obsValue(obsCleanIndex)  = ensObs%obsValue(obsIndex)
        ensObsClean%assFlag(obsCleanIndex)   = ensObs%assFlag(obsIndex)
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

    ensObsOut%lat(:)       = ensObsIn%lat(:)
    ensObsOut%lon(:)       = ensObsIn%lon(:)
    ensObsOut%logPres(:)   = ensObsIn%logPres(:)
    ensObsOut%varObsInv(:) = ensObsIn%varObsInv(:)
    ensObsOut%Yb_r4(:,:)   = ensObsIn%Yb_r4(:,:)
    ensObsOut%meanYb(:)    = ensObsIn%meanYb(:)
    ensObsOut%deterYb(:)   = ensObsIn%deterYb(:)
    ensObsOut%obsValue(:)  = ensObsIn%obsValue(:)
    ensObsOut%assFlag(:)   = ensObsIn%assFlag(:)

  end subroutine eob_copy

  !--------------------------------------------------------------------------
  ! eob_allGather
  !--------------------------------------------------------------------------
  subroutine eob_allGather(ensObs,ensObs_mpiglobal)
    ! :Purpose: Collect obs information distributed over all mpi tasks and
    !           make it available on all mpi tasks. The output ensObs object
    !           will be allocated within this subroutine.
    implicit none

    ! arguments
    type(struct_eob) :: ensObs
    type(struct_eob) :: ensObs_mpiglobal

    ! locals
    integer :: ierr, procIndex, memberIndex, numObs_mpiglobal
    integer :: allNumObs(mpi_nprocs), displs(mpi_nprocs)

    write(*,*) 'eob_allGather: starting'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

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
    call rpn_comm_gatherv( ensObs%varObsInv          , ensObs%numObs, 'mpi_real8', &
                           ensObs_mpiglobal%varObsInv, allNumObs, displs, 'mpi_real8',  &
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

    call rpn_comm_bcast(ensObs_mpiglobal%lat, ensObs_mpiglobal%numObs, 'mpi_real8',  &
                        0, 'GRID', ierr)
    call rpn_comm_bcast(ensObs_mpiglobal%lon, ensObs_mpiglobal%numObs, 'mpi_real8',  &
                        0, 'GRID', ierr)
    call rpn_comm_bcast(ensObs_mpiglobal%logPres, ensObs_mpiglobal%numObs, 'mpi_real8',  &
                        0, 'GRID', ierr)
    call rpn_comm_bcast(ensObs_mpiglobal%obsValue, ensObs_mpiglobal%numObs, 'mpi_real8',  &
                        0, 'GRID', ierr)
    call rpn_comm_bcast(ensObs_mpiglobal%varObsInv, ensObs_mpiglobal%numObs, 'mpi_real8',  &
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

    write(*,*) 'eob_allGather: total number of obs to be assimilated =', sum(ensObs_mpiglobal%assFlag(:))

    write(*,*) 'eob_allGather: finished'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine eob_allGather

  !--------------------------------------------------------------------------
  ! eob_getLocalBodyIndices
  !--------------------------------------------------------------------------
  function eob_getLocalBodyIndices(ensObs,localBodyIndices,distances,lat,lon,logPres,  &
                                   hLocalize,vLocalize,numLocalObsFound) result(numLocalObs)
    ! :Purpose: Return a list of values of bodyIndex for all observations within 
    !           the local volume around the specified lat/lon used for assimilation
    !           (as defined by h/vLocalize). The kdtree2 module is used to efficiently
    !           perform this task. The kdtree itself is constructed on the first call.
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

    ! locals
    integer :: obsIndex

    call obs_extractObsRealHeaderColumn(ensObs%lat, ensObs%obsSpaceData, OBS_LAT)
    call obs_extractObsRealHeaderColumn(ensObs%lon, ensObs%obsSpaceData, OBS_LON)
    call obs_extractObsRealBodyColumn(ensObs%obsValue, ensObs%obsSpaceData, OBS_VAR)
    call obs_extractObsRealBodyColumn(ensObs%varObsInv, ensObs%obsSpaceData, OBS_OER)
    do obsIndex = 1, ensObs%numObs
      if(ensObs%varObsInv(obsIndex) > 0.0d0) then
        ensObs%varObsInv(obsIndex) = 1.0d0/(ensObs%varObsInv(obsIndex)**2)
      else
        ensObs%varObsInv(obsIndex) = 0.0d0
      end if
    end do

  end subroutine eob_setLatLonObs

  !--------------------------------------------------------------------------
  ! eob_setPres
  !--------------------------------------------------------------------------
  subroutine eob_setPres(ensObs, columnMeanTrl)
    ! :Purpose: Set the ln(pressure) value for each observation that 
    !           will be used when doing vertical localization. For
    !           radiance observations, the level of the maximum value
    !           of the derivative of transmission is used.
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

      if( varNumber(obsIndex) == BUFR_NETS .or. varNumber(obsIndex) == BUFR_NEPS .or.  &
          varNumber(obsIndex) == BUFR_NEUS .or. varNumber(obsIndex) == BUFR_NEVS .or.  &
          varNumber(obsIndex) == BUFR_NESS .or. varNumber(obsIndex) == BUFR_NEPN .or. &
          varNumber(obsIndex) == BUFR_VIS  .or. varNumber(obsIndex) == BUFR_GUST .or. &
          varNumber(obsIndex) == BUFR_radarPrecip ) then

        ! all surface observations
        ensObs%logPres(obsIndex) = log(sfcPres_ptr(1,headerIndex))

      else if (varNumber(obsIndex) == BUFR_NEZD) then

        ! ZTD observation, try 0.7*Psfc (i.e. ~700hPa when Psfc=1000hPa)
        ensObs%logPres(obsIndex) = log(0.7D0 * sfcPres_ptr(1,headerIndex))

      else if (obsPPP(obsIndex) > 0.0d0 .and. obsVcoCode(obsIndex)==2) then

        ! all pressure level observations
        ensObs%logPres(obsIndex) = log(obsPPP(obsIndex))

      else if(obsVcoCode(obsIndex)==1) then

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
            write(*,*) 'eob_setPres for tovs: ', codType(obsIndex), obsPPP(obsIndex), 0.01*exp(ensObs%logPres(obsIndex))
          end if
        else
          ensObs%logPres(obsIndex) = log(500.0D2)
        end if

      else if(ensObs%assFlag(obsIndex)==1) then

        write(*,*) 'eob_setLatLonPresObs: ERROR! cannot compute pressure for this observation: ',  &
                   obsPPP(obsIndex), varNumber(obsIndex), obsVcoCode(obsIndex)
        call utl_abort('eob_setPres')

      end if
    end do

  end subroutine eob_setPres

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
  ! eob_calcRemoveMeanYb
  !--------------------------------------------------------------------------
  subroutine eob_calcRemoveMeanYb(ensObs)
    implicit none

    ! arguments
    type(struct_eob) :: ensObs

    ! locals
    integer :: obsIndex

    do obsIndex = 1, ensObs%numObs
      ensObs%meanYb(obsIndex) = sum(ensObs%Yb_r4(:,obsIndex)) / ensObs%numMembers
      ensObs%Yb_r4(:,obsIndex) = ensObs%Yb_r4(:,obsIndex) - ensObs%meanYb(obsIndex)
    end do

  end subroutine eob_calcRemoveMeanYb

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
  ! eob_setMeanHPHT
  !--------------------------------------------------------------------------
  subroutine eob_setMeanHPHT(ensObs)
    implicit none

    ! arguments
    type(struct_eob) :: ensObs

    ! locals
    integer :: obsIndex, memberIndex
    real(8) :: hpht

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
    end do

  end subroutine eob_setMeanHPHT

  !--------------------------------------------------------------------------
  ! max_transmission
  !--------------------------------------------------------------------------
  subroutine max_transmission(transmission, numLevels, transIndex, rttovPres, maxLnP)
    ! :Purpose: Determine the height in log pressure where we find the maximum 
    !           value of the first derivative of transmission with respect to 
    !           log pressure
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
  ! get_peak
  !--------------------------------------------------------------------------
  subroutine get_peak(maxIndex,nlev,lnp,deriv,maxLnP)
    ! :Purpose: Do quadratic interpolation to find pressure of peak transmission.
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
