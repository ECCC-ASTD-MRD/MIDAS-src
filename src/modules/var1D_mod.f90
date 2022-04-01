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

module var1D_mod
  ! MODULE var1D_mod (prefix='var1D' category='3. High-level transformations')
  !
  ! :Purpose: contains all 1Dvar-related methods.
  !
  use columnData_mod
  use columnVariableTransforms_mod
  use controlVector_mod
  use gridStatevector_mod
  use horizontalCoord_mod
  use mpi_mod 
  use obsSpaceData_mod
  use timeCoord_mod
  use utilities_mod
  use verticalCoord_mod
  use codeprecision_mod
  use tovs_nl_mod
  use mathphysconstants_mod

  implicit none
  save
  private

  ! public procedures
  public :: var1D_Setup,var1D_Finalize
  public :: var1D_bsetup
  public :: var1D_transferColumnToYGrid, var1D_get1DVarIncrement
  public :: var1D_sqrtB, var1D_sqrtBT

  type(struct_hco), pointer :: hco_yGrid
  logical             :: initialized = .false.
  integer             :: nlev_M, nlev_T, nkgdim
  integer             :: cvDim_mpilocal
  type(struct_vco), pointer :: vco_anl
  integer, parameter   :: maxNumLevels=200
  integer              :: nulbgst=0
  real(8), allocatable :: bMatrix(:,:)
  real(8), allocatable :: bSqrtLand(:,:,:), bSqrtSea(:,:,:)
  real(4), allocatable :: latLand(:), lonLand(:), latSea(:), lonSea(:)
  integer              :: nLonLatPosLand, nLonLatPosSea, var1D_varCount
  character(len=4), allocatable :: var1D_varList(:)
  integer, external    :: get_max_rss
  integer, allocatable :: var1D_validHeaderIndex(:)    ! pointeur vers les colonnes assimilables pour minimiser la taille du vecteur de controle
  integer              :: var1D_validHeaderCount !taille effective de  var1D_validHeaderIndex
  type(struct_vco), target  :: vco_1Dvar
  integer,          parameter :: numMasterBmat = 1
  character(len=4), parameter :: masterBmatTypeList (numMasterBmat) = (/'HI'   /)
  character(len=8), parameter :: masterBmatLabelList(numMasterBmat) = (/'B_HI' /)
  logical,          parameter :: masterbmatIs3dList (numMasterBmat) = (/.true. /) 
  integer            :: numBmat
  integer, parameter :: numBmatMax = 10
  character(len=4) :: bmatTypeList  (numBmatMax)
  character(len=9) :: bmatLabelList (numBmatMax)
  integer          :: bmatInstanceID(numBmatMax)
  logical          :: bmatIs3dList  (numBmatMax)
  logical          :: bmatActive    (numBmatMax)
  !Namelist variables
  real(8)             :: scaleFactor(maxNumLevels)    ! scaling factors for variances (no
  real(8)             :: scaleFactorLQ(maxNumLevels)  ! scaling factors for LQ variances


contains

  !--------------------------------------------------------------------------
  !  var1D_setup
  !--------------------------------------------------------------------------
  subroutine var1D_setup(vco_in, obsdat, CVDIM_OUT)
    !
    ! :Purpose: to setup var1D module
    !
    implicit none
    ! arguments:
    type(struct_vco), pointer, intent(in):: vco_in
    type (struct_obs), intent(in)        :: obsdat
    integer, intent(out)                 :: cvDim_out
    ! locals:
    integer :: levelIndex, nulnam, ierr
    integer, external ::  fnom, fclos
    integer :: status, Vcode_anl
    logical :: fileExists
    type(struct_vco), pointer :: vco_file => null()
    character(len=24) :: oneDBmatLand = './Bmatrix_land.bin'
    character(len=24) :: oneDBmatSea = './Bmatrix_sea.bin'
    integer :: extractDate, locationIndex, countGood, headerIndex
    integer :: bodyStart, bodyEnd, bodyIndex
 
    NAMELIST /NAMVAR1D/ scaleFactor, scaleFactorLQ

    call utl_tmg_start(50,'--Bmatrix')
    call utl_tmg_start(51,'----B_HI_Setup')
    if(mpi_myid == 0) write(*,*) 'var1D_setup: Starting'
    if(mpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

    ! default values for namelist variables
    scaleFactor(:) = 1.0d0
    scaleFactorLQ(:) = 1.0d0
    nulnam = 0
    ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    read(nulnam, nml=namvar1D, iostat=ierr)
    if ( ierr /= 0 ) call utl_abort( 'var1D_setup: Error reading namelist' )
    if ( mpi_myid == 0 ) write( *, nml = namvar1D )
    ierr = fclos( nulnam )

    do levelIndex = 1, maxNumLevels
      if( scaleFactor( levelIndex ) > 0.0d0 ) then 
        scaleFactor( levelIndex ) = sqrt( scaleFactor( levelIndex ))
      else
        scaleFactor( levelIndex ) = 0.0d0
      end if
    end do
    if ( sum( scaleFactor( 1 : maxNumLevels ) ) == 0.0d0 ) then
      if ( mpi_myid == 0 ) write(*,*) 'var1D_setup: scaleFactor=0, skipping rest of setup'
      cvdim_out = 0
      call tmg_stop(51)
      call tmg_stop(50)
      return
    end if
    do levelIndex = 1, maxNumLevels
      if(scaleFactorLQ(levelIndex) > 0.0d0) then 
        scaleFactorLQ(levelIndex) = sqrt(scaleFactorLQ(levelIndex))
      else
        scaleFactorLQ(levelIndex) = 0.0d0
      end if
    end do
    vco_anl => vco_in
    nLev_M = vco_anl%nlev_M
    nLev_T = vco_anl%nlev_T
    if (mpi_myid == 0) write(*,*) 'var1D_setup: nLev_M, nLev_T=',nLev_M, nLev_T

    if (mpi_myid == 0) write(*,*) 'var1D_setup: Read 1DVar background statistics'
    inquire(file=trim(oneDBmatLand), exist=fileExists)
    if ( fileExists ) then
      ierr = fnom(nulbgst, trim(oneDBmatLand), 'FTN+SEQ+UNF+OLD+R/O', 0)
    else
      call utl_abort('var1D_setup: No 1DVar BACKGROUND STAT FILE ' // trim(oneDBmatLand))
    end if
    read(nulbgst) extractDate, vco_1Dvar%nLev_T, vco_1Dvar%nLev_M, vco_1Dvar%Vcode, &
         vco_1Dvar%ip1_sfc, vco_1Dvar%ip1_T_2m, vco_1Dvar%ip1_M_10m, var1D_varCount, nkgdim, nLonLatPosLand
    allocate( vco_1Dvar%ip1_T(nLev_T), vco_1Dvar%ip1_M(nLev_M) )
    allocate( var1D_varList(var1D_varCount) )
    allocate( bMatrix(nkgdim,nkgdim) )
    allocate( latLand(nLonLatPosLand), lonLand(nLonLatPosLand))       
    allocate( bSqrtLand(nLonLatPosLand, nkgdim, nkgdim) )
    read(nulbgst) vco_1Dvar%ip1_T(:), vco_1Dvar%ip1_M(:), var1D_varList(:)
    do locationIndex = 1, nLonLatPosLand
      read(nulbgst) latLand(locationIndex), lonLand(locationIndex), bMatrix(:,:)      
      bSqrtLand(locationIndex, :, :) = bMatrix(:, :)
      call utl_matsqrt(bSqrtLand(locationIndex, :, :), nkgdim, 1.d0, printInformation_opt=.false. )
    end do
    ierr = fclos(nulbgst)

    inquire(file=trim(oneDBmatSea), exist=fileExists)
    if ( fileExists ) then
      ierr = fnom(nulbgst, trim(oneDBmatSea), 'FTN+SEQ+UNF+OLD+R/O', 0)
    else
      call utl_abort('var1D_setup: No 1DVar BACKGROUND STAT FILE ' // trim(oneDBmatSea))
    end if
    read(nulbgst) extractDate, vco_1Dvar%nLev_T, vco_1Dvar%nLev_M, vco_1Dvar%Vcode, &
         vco_1Dvar%ip1_sfc, vco_1Dvar%ip1_T_2m, vco_1Dvar%ip1_M_10m, var1D_varCount, nkgdim, nLonLatPosSea
    allocate( bSqrtSea(nLonLatPosSea, nkgdim, nkgdim) )
    allocate( latSea(nLonLatPosSea), lonSea(nLonLatPosSea))
    read(nulbgst) vco_1Dvar%ip1_T(:), vco_1Dvar%ip1_M(:), var1D_varList(:)
    do locationIndex = 1, nLonLatPosSea
      read(nulbgst) latSea(locationIndex), lonSea(locationIndex), bMatrix(:,:)
      bSqrtSea(locationIndex, :, :) = bMatrix(:, :)
      call utl_matsqrt(bSqrtSea(locationIndex,:,:), nkgdim, 1.d0, printInformation_opt=.false. )
    end do
    ierr = fclos(nulbgst)

    vco_1Dvar%initialized = .true.
    vco_1Dvar%vGridPresent = .false.
    vco_file => vco_1Dvar
    if (.not. vco_equal(vco_anl,vco_file)) then
      call utl_abort('var1D_setup: vco from analysisgrid and cov file do not match')
    end if
    status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)
    if(Vcode_anl /= 5002 .and. Vcode_anl /= 5005) then
      write(*,*) 'Vcode_anl = ',Vcode_anl
      call utl_abort('var1D_setup: unknown vertical coordinate type!')
    end if
    if (.not. (gsv_varExist(varName='TT') .and.  &
               gsv_varExist(varName='UU') .and.  &
               gsv_varExist(varName='VV') .and.  &
               (gsv_varExist(varName='HU').or.gsv_varExist(varName='LQ')) .and.  &
               gsv_varExist(varName='P0')) ) then
      call utl_abort('var1D_setup: Some or all weather fields are missing. If it is desired to deactivate&
           & the weather assimilation, then all entries of the array SCALEFACTOR in the namelist NAMVAR1D&
           & should be set to zero.')
    end if
    if (.not. gsv_varExist(varName='TG')) then
      write(*,*) 'var1D_setup: WARNING: The TG field is missing. This must be present when assimilating'
      write(*,*) 'radiance observations.'
    end if

    !we want to count how many obs are really assimilable to minimize controlvector size
    var1D_validHeaderCount = 0
    allocate(  var1D_validHeaderIndex(obs_numHeader(obsdat)) )
    do headerIndex = 1, obs_numHeader(obsdat)
      bodyStart = obs_headElem_i(obsdat, OBS_RLN, headerIndex)
      bodyEnd = obs_headElem_i(obsdat, OBS_NLV, headerIndex) + bodyStart - 1
      countGood = 0
      do bodyIndex = bodyStart, bodyEnd
        if (obs_bodyElem_i(obsdat, OBS_ASS, bodyIndex) == obs_assimilated) countGood = countGood + 1
      end do
      if (countGood > 0)  then
        var1D_validHeaderCount = var1D_validHeaderCount + 1
        var1D_validHeaderIndex(var1D_validHeaderCount) = headerIndex
        if (var1D_validHeaderCount == 1) write(*,*) 'first OBS', headerIndex
      end if
    end do
    write(*,*) 'var1D_setup: var1D_validHeaderCount, obs_numHeader(obsdat)', var1D_validHeaderCount, obs_numHeader(obsdat)
    cvDim_out = nkgdim * var1D_validHeaderCount
    cvDim_mpilocal = cvDim_out

    initialized = .true.

    call tmg_stop(51)
    call tmg_stop(50)

  end subroutine var1D_setup
  
  !--------------------------------------------------------------------------
  ! var1D_bSqrtHi
  !--------------------------------------------------------------------------
  subroutine var1D_bSqrtHi(controlVector_in, column, dataObs)
    !
    ! :Purpose: HI component of B square root in 1DVar mode
    !
    implicit none
    ! arguments:
    real(8), intent(in)                    :: controlVector_in(cvDim_mpilocal)
    type(struct_columnData), intent(inout) :: column
    type(struct_obs), intent(in)           :: dataObs
    ! locals:
    integer :: headerIndex, latitudeBandIndex(1), varIndex, columnIndex
    real(8), pointer :: currentColumn(:)
    real(8), allocatable ::  oneDProfile(:)
    real(8) :: latitude
    integer :: surfaceType, offset

    if (mpi_myid == 0) write(*,*) 'var1D_bsqrtHi: starting'
    if (mpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

    if (.not. initialized) then
      if (mpi_myid == 0) write(*,*) 'var1D_bsqrtHi: 1Dvar B matrix not initialized'
      return
    end if
    allocate(oneDProfile(nkgdim))
    do columnIndex = 1, var1D_validHeaderCount 
      headerIndex = var1D_validHeaderIndex(columnIndex)
      latitude = obs_headElem_r(dataObs, OBS_LAT, headerIndex) !radian 
      surfaceType = tvs_ChangedStypValue(dataObs, headerIndex)
      if (surfaceType == 1) then !Sea
        latitudeBandIndex = minloc( abs( latitude - latSea(:)) )
        oneDProfile(:) = matmul(bSqrtSea(latitudeBandIndex(1), :, :), controlVector_in(1+(columnIndex-1)*nkgdim:columnIndex*nkgdim))
      else ! Land or Sea Ice
        latitudeBandIndex = minloc( abs( latitude - latLand(:)) )
        oneDProfile(:) = matmul(bSqrtLand(latitudeBandIndex(1), :, :), controlVector_in(1+(columnIndex-1)*nkgdim:columnIndex*nkgdim))
      end if
      offset = 0
      do varIndex = 1, var1D_varCount
        currentColumn => col_getColumn(column, headerIndex, varName_opt=var1D_varList(varIndex))
        currentColumn(:) = 0.d0
        currentColumn(:) = oneDProfile(offset+1:offset+size(currentColumn))
        offset = offset + size(currentColumn)
      end do
      if (offset /= nkgdim) then
        write(*,*) 'var1D_bsqrtHi: offset, nkgdim', offset, nkgdim
        call utl_abort('var1D_bSqrtHi: inconsistency between Bmatrix and statevector size')
      end if
    end do
    deallocate(oneDProfile)
    if (mpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
    if (mpi_myid == 0) write(*,*) 'var1D_bSqrtHi: done'

  end subroutine var1D_bSqrtHi

  !--------------------------------------------------------------------------
  ! var1D_bSqrtHiAd
  !--------------------------------------------------------------------------
  subroutine var1D_bSqrtHiAd(controlVector_in, column, dataObs)
    !
    ! :Purpose: HI component of B square root adjoint in 1DVar mode
    !
    implicit none
    ! arguments:
    real(8), intent(inout)                 :: controlVector_in(cvDim_mpilocal)
    type(struct_columnData), intent(inout) :: column
    type (struct_obs), intent(in)          :: dataObs
    ! locals:
    integer :: headerIndex, latitudeBandIndex(1), varIndex, columnIndex
    real(8), pointer :: currentColumn(:)
    real(8), allocatable ::  oneDProfile(:)
    real(8) :: latitude
    integer :: surfaceType, offset

    if (mpi_myid == 0) write(*,*) 'var1D_bSqrtHiAd: starting'
    if (mpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
    if (.not. initialized) then
      if (mpi_myid == 0) write(*,*) 'var1D_bSqrtHiAd: 1dvar Bmatrix not initialized'
      return
    end if
    allocate(oneDProfile(nkgdim))
    do columnIndex = 1, var1D_validHeaderCount
      headerIndex = var1D_validHeaderIndex(columnIndex)
      offset = 0
      do varIndex = 1, var1D_varCount
        currentColumn => col_getColumn(column, headerIndex, varName_opt=var1D_varList(varIndex))
        oneDProfile(offset+1:offset+size(currentColumn)) = currentColumn(:)
        offset = offset + size(currentColumn)
      end do
      if (offset /= nkgdim) then
        write(*,*) 'var1D_bSqrtHiAd: offset, nkgdim', offset, nkgdim
        call utl_abort('var1D_bSqrtHiAd: inconsistency between Bmatrix and statevector size')
      end if
      latitude = obs_headElem_r(dataObs, OBS_LAT, headerIndex) !radian
      surfaceType =  tvs_ChangedStypValue(dataObs, headerIndex)
      if (surfaceType == 1) then !Sea
        latitudeBandIndex = minloc( abs( latitude - latSea(:)) )
        controlVector_in(1+(columnIndex-1)*nkgdim:columnIndex*nkgdim) =  &
             controlVector_in(1+(columnIndex-1)*nkgdim:columnIndex*nkgdim) + &
             matmul(bSqrtSea(latitudeBandIndex(1),:,:), oneDProfile)
      else ! Land or Sea Ice
        latitudeBandIndex = minloc( abs( latitude - latLand(:)) )
        controlVector_in(1+(columnIndex-1)*nkgdim:columnIndex*nkgdim) =  &
             controlVector_in(1+(columnIndex-1)*nkgdim:columnIndex*nkgdim) + &
             matmul(bSqrtLand(latitudeBandIndex(1),:,:), oneDProfile)
      end if
    end do
    deallocate(oneDProfile)
    if (mpi_myid == 0) write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
    if (mpi_myid == 0) write(*,*) 'var1D_bSqrtHiAd: done'

  end subroutine var1D_bSqrtHiAd

  !--------------------------------------------------------------------------
  ! var1D_Finalize
  !--------------------------------------------------------------------------
  subroutine var1D_Finalize()
    !
    ! :Purpose: to deallocate memory used by internal module strucures
    !
    implicit none
    if (initialized) then
       deallocate( bMatrix)
       deallocate( bSqrtLand )
       deallocate( bSqrtSea )
       deallocate( latLand, lonLand, latSea, lonSea )
       deallocate( var1D_validHeaderIndex )
    end if

  end subroutine var1D_Finalize

  !--------------------------------------------------------------------------
  ! var1D_transferColumnToYGrid
  !--------------------------------------------------------------------------
  subroutine var1D_transferColumnToYGrid( stateVector, obsSpaceData, column)
    !
    ! :Purpose: to transfer content of a columndata object to a statevector object
    !           without interpolation (to be used in 1DVar mode to write increments on Y grid).
    !
    implicit none
    ! arguments:
    type(struct_gsv), intent(in)           :: stateVector
    type(struct_obs), intent(in)           :: obsSpaceData
    type(struct_columnData), intent(inout) :: column
    ! locals:
    integer :: varIndex, globalObsIndex, obsIndex, taskIndex, headerIndex
    integer, allocatable :: var1D_validHeaderCountAllTasks(:), obsOffset(:)
    real(8), pointer :: myColumn(:), myField(:,:,:)
    real(8), allocatable, target :: dummy(:)
    integer :: var1D_validHeaderCountMpiGlobal, var1D_validHeaderCountMax, ierr, status
    real(8) :: lat, lon
    integer :: varDim, tag

    call rpn_comm_barrier("GRID",ierr)
    allocate( obsOffset(0:mpi_nprocs-1) )
    if (mpi_myid ==0) then
      allocate( var1D_validHeaderCountAllTasks(mpi_nprocs) )
    else
      allocate(var1D_validHeaderCountAllTasks(1))
    end if

    call rpn_comm_gather(var1D_validHeaderCount  , 1, 'MPI_INTEGER', var1D_validHeaderCountAllTasks, 1,'MPI_INTEGER', 0, "GRID", ierr )
    if (mpi_myId ==0) then
      var1D_validHeaderCountMpiGlobal = sum( var1D_validHeaderCountAllTasks(:) )
      var1D_validHeaderCountMax = maxval( var1D_validHeaderCountAllTasks(:) )
      obsOffset(0) = 0
      do taskIndex = 1, mpi_nprocs - 1
        obsOffset(taskIndex) = obsOffset(taskIndex - 1) + var1D_validHeaderCountAllTasks(taskIndex)
      end do
      write(*,*) 'var1D_transferColumnToYGrid: obsOffset: ', obsOffset(:)
    end if
    call rpn_comm_bcast( obsOffset, mpi_nprocs, 'MPI_INTEGER', 0,  "GRID",ierr )
    call rpn_comm_bcast( var1D_validHeaderCountMax, 1, 'MPI_INTEGER', 0,  "GRID",ierr )

    call hco_setupYgrid(hco_Ygrid, 1, var1D_validHeaderCountMpiGlobal)
    if (mpi_myId ==0) then
      call gsv_allocate(stateVector, numstep=tim_nstepobsinc, hco_ptr=hco_Ygrid, vco_ptr=column%vco, &
           datestamp_opt=tim_getDatestamp(), mpi_local_opt=.false., &
           dataKind_opt=pre_incrReal, allocHeight_opt=.false., allocPressure_opt=.false., &
           besilent_opt=.false.)
      write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'	
    end if

    write(*,*) 'var1D_transferColumnToYGrid: start of lat-lon dissemination'
    do obsIndex = 1, var1D_validHeaderCountMax
      if (obsIndex <= var1D_validHeaderCount ) then
        headerIndex = var1D_validHeaderIndex(obsIndex)      
        lat = obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex)
        lon = obs_headElem_r(obsSpaceData, OBS_LON, headerIndex)
      else
        lat = MPC_missingValue_R8
        lon = MPC_missingValue_R8
      end if
      if (mpi_myId == 0) then
        if ( obsIndex <= var1D_validHeaderCount ) then
          hco_yGrid%lat2d_4(1, obsIndex) = lat
          hco_yGrid%lon2d_4(1, obsIndex) = lon
        end if
      else
        tag = 2 * mpi_myID
        call rpn_comm_send( lat, 1, 'mpi_real8', 0, tag,     'GRID', ierr )
        call rpn_comm_send( lon, 1, 'mpi_real8', 0, tag + 1, 'GRID', ierr )
      end if

      if (mpi_myId == 0) then
        do taskIndex = 1,  mpi_nprocs - 1
          tag = 2 * taskIndex
          call rpn_comm_recv( lat, 1, 'mpi_real8', taskIndex, tag, 'GRID', status, ierr )
          call rpn_comm_recv( lon, 1, 'mpi_real8', taskIndex, tag+1, 'GRID', status, ierr )
          if (lat /= MPC_missingValue_R8 .and. lon /= MPC_missingValue_R8) then 
            globalObsIndex = obsIndex + obsOffset(taskIndex)
            hco_yGrid%lat2d_4(1, globalObsIndex) = lat
            hco_yGrid%lon2d_4(1, globalObsIndex) = lon
          end if
        end do
      end if
      call rpn_comm_barrier("GRID",ierr)
    end do

    call rpn_comm_barrier("GRID",ierr)
    write(*,*) 'var1D_transferColumnToYGrid: end of lat-lon dissemination'
    
    do varIndex = 1, var1D_varCount
      write(*,*) 'var1D_transferColumnToYGrid: start of dissemination for ', var1D_varList(varIndex)
      if (mpi_myId == 0 ) then
        call gsv_getField(stateVector, myField, varName_opt=var1D_varList(varIndex), stepIndex_opt=1)
        varDim = gsv_getNumLevFromVarName(stateVector, var1D_varList(varIndex))
      end if
      call rpn_comm_bcast(varDim, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
      allocate(dummy(varDim))
      dummy(:) = MPC_missingValue_R8
      do obsIndex = 1, var1D_validHeaderCountMax
        if (obsIndex <= var1D_validHeaderCount ) then
          headerIndex = var1D_validHeaderIndex(obsIndex) 
          myColumn => col_getColumn(column, headerIndex, varName_opt=var1D_varList(varIndex))
        else
          myColumn => dummy
        end if
        if (mpi_myId == 0) then
          if ( obsIndex <= var1D_validHeaderCount ) then
            myField(1, obsIndex, :) = myColumn(:)
          end if
        else
          tag = mpi_myId
          call rpn_comm_send(myColumn , varDim, 'mpi_real8', 0, tag, 'GRID', ierr )
        end if

        if (mpi_myId == 0) then
          do taskIndex = 1,  mpi_nprocs - 1
            tag = taskIndex
            call rpn_comm_recv(myColumn,  varDim, 'mpi_real8', taskIndex, tag, 'GRID', status, ierr )
            if (all( myColumn /=  MPC_missingValue_R8)) then
              globalObsIndex = obsIndex + obsOffset(taskIndex)
              myField(1, globalObsIndex, :) = myColumn(:)
            end if
          end do
        end if
      end do

      write(*,*) 'var1D_transferColumnToYGrid: end of dissemination for ', var1D_varList(varIndex)
      deallocate(dummy)

    end do

    call rpn_comm_barrier("GRID",ierr)
    deallocate( obsOffset )
    deallocate( var1D_validHeaderCountAllTasks )
    

  end subroutine var1D_transferColumnToYGrid

  !--------------------------------------------------------------------------
  ! var1D_get1DVarIncrement
  !--------------------------------------------------------------------------
  subroutine var1D_get1DVarIncrement(incr_cv, column, columnTrlOnAnlIncLev, &
                                     obsSpaceData, nvadim_mpilocal)
    !
    ! :Purpose: to compute 1Dvar increment from control vector
    !
    implicit none
    ! arguments:
    real(8), intent(in)                    :: incr_cv(:)
    type(struct_columnData), intent(inout) :: column
    type(struct_columnData), intent(in)    :: columnTrlOnAnlIncLev
    type(struct_obs),        intent(in)    :: obsSpaceData
    integer,                 intent(in)    :: nvadim_mpilocal
    ! compute increment from control vector (multiply by B^1/2)
    call var1D_sqrtB(incr_cv, nvadim_mpilocal, column, obsSpaceData)
    call cvt_transform(column, 'ZandP_tl', columnTrlOnAnlIncLev)

  end subroutine var1D_get1DVarIncrement

  !--------------------------------------------------------------------------
  ! var1D_bsetup
  !--------------------------------------------------------------------------
  subroutine var1D_bsetup(vco_anl, obsdat)
    !
    !:Purpose: To initialize the 1Dvar analysis Background term.
    !
    implicit none
    ! arguments:
    type(struct_vco), pointer, intent(in) :: vco_anl
    type (struct_obs), intent(in)         :: obsdat
    ! locals:
    integer, allocatable :: cvDimPerInstance(:)
    integer :: cvdim
    integer :: masterBmatIndex, bMatInstanceIndex, nBmatInstance, bmatIndex
    character(len=2) :: bMatInstanceIndexString
    character(len=3) :: bMatExtraLabel
    logical :: active
    !
    !- 1.  Setup the B matrices
    !
    numBmat = 0
    do masterBmatIndex = 1, numMasterBmat
      select case( trim(masterBmatTypeList(masterBmatIndex)) )
      case ('HI')
        !- 1.1 Time-Mean Homogeneous and Isotropic...
        nBmatInstance = 1 ! hardwired
        allocate(cvdimPerInstance(nBmatInstance))
        write(*,*) 'var1D_bsetup: Setting up the modular GLOBAL HI 1D covariances...'
        call var1D_Setup(vco_anl, obsdat, cvdim)
        write(*,*) " var1D_bsetup: cvdim= ", cvdim
        cvdimPerInstance(1) = cvdim
      case default
        call utl_abort( 'var1D_bSetup: requested bmatrix type does not exist ' // trim(masterBmatTypeList(masterBmatIndex)) )
      end select
      !- 1.2 Append the info to the B matrix info arrays and setup the proper control sub-vectors
      do bMatInstanceIndex = 1, nBmatInstance
        numBmat = numBmat + 1
        if (nBmatInstance == 1) then
          bMatExtraLabel= ""
        else
          write(bMatInstanceIndexString,'(I2.2)') bMatInstanceIndex 
          bMatExtraLabel="_"//trim(bMatInstanceIndexString)
        end if
        bmatLabelList (numBmat) = trim(masterbmatLabelList(masterBmatIndex))//trim(bMatExtraLabel)
        bmatTypeList  (numBmat) = masterBmatTypeList(masterBmatIndex)
        bmatIs3dList  (numBmat) = masterbmatIs3dList(masterBmatIndex)
        bmatInstanceID(numBmat) = bMatInstanceIndex
        call cvm_setupSubVector(bmatLabelList(numBmat), bmatTypeList(numBmat), cvdimPerInstance(bMatInstanceIndex))
      end do
      deallocate(cvdimPerInstance)
    end do
    !
    !- 2. Print a summary and set the active B matrices array
    !
    write(*,*)
    write(*,*) "var1D_setup SUMMARY, number of B matrices found = ", numBmat
    do bmatIndex = 1, numBmat
      write(*,*) "  B matrix #", bmatIndex
      active = cvm_subVectorExists(bmatLabelList(bmatIndex))
      if (active) then
        write(*,*) "   ACTIVE"
      else
        write(*,*) "   NOT USED"
      end if
      write(*,*) "     -> label       = ", bmatLabelList (bmatIndex)
      write(*,*) "     -> type        = ", bmatTypeList  (bmatIndex)
      if (active) then
        write(*,*) "     -> is 3D       = ", bmatIs3dList  (bmatIndex)
        write(*,*) "     -> instance ID = ", bmatInstanceID(bmatIndex)
      end if
      bmatActive(bmatIndex) = active
    end do

  end subroutine var1D_bsetup

  !--------------------------------------------------------------------------
  ! var1D_sqrtB
  !-------------------------------------------------------------------------- 
  subroutine var1D_sqrtB(controlVector, cvdim, column, obsSpaceData)
    !
    !:Purpose: To transform model state from control-vector space to grid-point
    !          space.    
    !
    implicit none
    ! arguments:
    integer, intent(in)                    :: cvdim
    real(8), intent(in)                    :: controlVector(cvdim)
    type(struct_columnData), intent(inout) :: column
    type(struct_obs), intent(in)           :: obsSpaceData
    ! locals:
    integer :: bmatIndex
    real(8), pointer :: subVector(:)

    call utl_tmg_start(50,'--Bmatrix')

    !
    !- 1.  Compute the analysis increment
    !
    bmat_loop: do bmatIndex = 1, numBmat
      if ( .not. bmatActive(bmatIndex) ) cycle bmat_loop
      subVector => cvm_getSubVector( controlVector, bmatLabelList(bmatIndex) )
      select case( trim(bmatTypeList(bmatIndex)) )
      case ('HI')
        !- 1.1 Time-Mean Homogeneous and Isotropic...
        call utl_tmg_start(52,'----B_HI_TL')
        call var1D_bsqrtHi( subVector,   &  ! IN
                            column,      &  ! OUT
                            obsspacedata )  ! IN
        call tmg_stop(52)
      case default
        call utl_abort( 'var1D_sqrtB: requested bmatrix type does not exist ' // trim(bmatTypeList(bmatIndex)) )
      end select
    end do bmat_loop

    call tmg_stop(50)

  end subroutine var1D_sqrtB

  !--------------------------------------------------------------------------
  ! var1D_sqrtBT
  !--------------------------------------------------------------------------
  subroutine var1D_sqrtBT(controlVector, cvdim, column, obsData)
    !
    !:Purpose: To transform model state from grid-point space to
    !          error-covariance space.
    !
    implicit none
    ! arguments:
    integer, intent(in)                    :: cvdim
    real(8), intent(in)                    :: controlVector(cvdim)
    type(struct_columnData), intent(inout) :: column
    type(struct_obs), intent(in)           :: obsData
    ! locals:
    integer :: bmatIndex
    real(8), pointer :: subVector(:)

    call utl_tmg_start(50,'--Bmatrix')

    ! Process components in opposite order as forward calculation
    bmat_loop: do bmatIndex = numBmat, 1, -1
      if ( .not. bmatActive(bmatIndex) ) cycle bmat_loop
      subVector => cvm_getSubVector( controlVector, bmatLabelList(bmatIndex) )
      select case( trim(bmatTypeList(bmatIndex)) )
      case ('HI')
        !- 2.5 Time-Mean Homogeneous and Isotropic...
        call utl_tmg_start(53,'----B_HI_AD')
        call var1D_bsqrtHiAd( subvector, &  ! IN
                              column,    &  ! OUT
                              obsData )     ! IN
        call tmg_stop(53)
      case default
        call utl_abort( 'var1D_sqrtBT: requested bmatrix type does not exist ' // trim(bmatTypeList(bmatIndex)) )
      end select
    end do bmat_loop

    call tmg_stop(50)

  end subroutine var1D_sqrtBT

end module var1D_mod
