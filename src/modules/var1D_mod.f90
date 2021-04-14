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
  ! MODULE BmatrixHI_mod (prefix='var1D' category='')
  !
  ! :Purpose: contains all 1dvar-related methods.
  !
  use columnData_mod
  use MathPhysConstants_mod
  use ObsSpaceData_mod
  use analysisgrid_mod
  use codePrecision_mod
  use columnData_mod
  use columnVariableTransforms_mod
  use controlVector_mod
  use earthConstants_mod
  use getGridPosition_mod
  use gridStatevector_mod
  use gridVariableTransforms_mod
  use gridstatevector_mod
  use horizontalCoord_mod
  use mathPhysConstants_mod
  use mpi, only : mpi_status_size ! this is the mpi library module
  use mpi_mod
  use mpivar_mod
  use mpivar_mod 
  use obsSpaceData_mod
  use obsTimeInterp_mod
  use physicsFunctions_mod
  use slantprofilelatlon_mod
  use timeCoord_mod
  use tt2phi_mod
  use utilities_mod
  use varNameList_mod
  use verticalCoord_mod



  implicit none
  save
  private

  ! public procedures
  public :: var1D_Setup,var1D_Finalize
  public :: var1D_bsetup
  public :: var1D_transferColumnToYGrid, var1D_get1DVarIncrement
  public :: var1D_sqrtB, var1D_sqrtBT

  public :: var1D_obsCount, var1D_obsPointer, var1D_varCount, var1D_varList


  type(struct_hco), target        :: hco_yGrid
  logical             :: initialized = .false.
  integer             :: nj_l,ni_l
  integer, save       :: nlev_M, nlev_T, nkgdim
  integer             :: cvDim_mpilocal
  type(struct_vco),pointer :: vco_anl
 
  ! originally from common blocks and possibly from the namelist:
  integer,parameter   :: maxNumLevels=200
  real(8)             :: scaleFactor(maxNumLevels)
  real(8)             :: scaleFactorLQ(maxNumLevels)
  real(8)             :: scaleFactorCC(maxNumLevels)
  logical             :: scaleTG
  logical             :: TweakTG
  integer             :: nulbgst=0
  

  ! this should come from state vector object
  integer             :: nspositUU 
  integer             :: nspositVV 
  integer             :: nspositTT 
  integer             :: nspositQ
  integer             :: nspositPS 
  integer             :: nspositTG

  !variables added for oneDvar in Midas
  real(8),allocatable :: bMatrix(:,:)!, bMatrixSea(:,:,:)
  real(8),allocatable :: bHalfLand(:,:,:), bHalfSea(:,:,:)
  real(4),allocatable :: latLand(:), lonLand(:), latSea(:), lonSea(:)
  integer              :: nLonLatPosLand, nLonLatPosSea, var1D_varCount
  character(len=4),allocatable :: var1D_varList(:)
  integer,external    :: get_max_rss
  integer,save,allocatable :: var1D_obsPointer(:) ! pointeur vers les colonnes assimilables pour minimiser la taille du vecteur de controle
  integer, save :: var1D_obsCount !taille effective de  var1D_obsPointer
  type(struct_vco),target  :: vco_1Dvar
 integer,          parameter :: numMasterBmat = 1
  character(len=4), parameter :: masterBmatTypeList (numMasterBmat) = (/'HI'   /) !, 'LATB'  , 'ENS'  , 'CHM'  , 'DIFF'  /)
  character(len=8), parameter :: masterBmatLabelList(numMasterBmat) = (/'B_HI' /) !, 'B_LATB', 'B_ENS', 'B_CHM', 'B_DIFF'/)
  logical,          parameter :: masterbmatIs3dList (numMasterBmat) = (/.true. /) !, .true.  , .false., .true. , .true.  /)

  integer            :: numBmat
  integer, parameter :: numBmatMax = 10

  character(len=4) :: bmatTypeList  (numBmatMax)
  character(len=9) :: bmatLabelList (numBmatMax)
  integer          :: bmatInstanceID(numBmatMax)
  logical          :: bmatIs3dList  (numBmatMax)
  logical          :: bmatActive    (numBmatMax)

contains

  subroutine var1D_setup(vco_in, obsdat, CVDIM_OUT)
    implicit none
    type(struct_vco),pointer :: vco_in
    type (struct_obs)        :: obsdat
    integer                  :: cvDim_out

    integer :: jlev, nulnam, ierr
    integer, external ::  fnom, fclos
    integer :: status, Vcode_anl
    logical :: fileExists
    type(struct_vco),pointer :: vco_file => null()
   
    character(len=24) :: oneDBmatLand = './Bmatrix_land.bin'
    character(len=24) :: oneDBmatSea = './Bmatrix_sea.bin'
    integer :: extractDate, locationIndex, countBad, headerIndex
    integer :: bodyStart, bodyEnd, bodyIndex
 
    NAMELIST /NAMVAR1D/scaleFactor,scaleFactorLQ,scaleFactorCC,scaleTG, &
         TweakTG

    call tmg_start(15,'BHI1D_SETUP')
    if(mpi_myid == 0) write(*,*) 'var1D_setup: starting'
    if(mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! default values for namelist variables
    scaleFactor(:) = 1.0d0
    scaleFactorLQ(:) = 1.0d0
    scaleFactorCC(:) = 1.0d0
    scaleTG = .true.
    TweakTG = .false.

    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=namvar1D,iostat=ierr)
    if ( ierr /= 0 ) call utl_abort( 'var1D_setup: Error reading namelist' )
    if ( mpi_myid == 0 ) write( *, nml = namvar1D )
    ierr = fclos( nulnam )

    do jlev = 1, maxNumLevels
      if( scaleFactor( jlev ) > 0.0d0 ) then 
        scaleFactor( jlev ) = sqrt( scaleFactor( jlev ))
      else
        scaleFactor( jlev ) = 0.0d0
      end if
    end do

    if ( sum( scaleFactor( 1 : maxNumLevels ) ) == 0.0d0 ) then
      if ( mpi_myid == 0 ) write(*,*) 'bmatrixHI: scaleFactor=0, skipping rest of setup'
      cvdim_out = 0
      call tmg_stop(15)
      return
    end if

    vco_anl => vco_in
    nLev_M = vco_anl%nlev_M
    nLev_T = vco_anl%nlev_T

    if (mpi_myid == 0) write(*,*) 'BHI1D_setup: nLev_M, nLev_T=',nLev_M, nLev_T
    
    if (mpi_myid == 0) write(*,*) 'BHI1D_setup: read 1DVar background statistics'
    inquire(file=trim(oneDBmatLand), exist=fileExists)
    if ( fileExists ) then
      ierr = fnom(nulbgst, trim(oneDBmatLand), 'FTN+SEQ+UNF+OLD+R/O', 0)
    else
      call utl_abort('BHI1D_setup:No 1DVar BACKGROUND STAT FILE ' // trim(oneDBmatLand))
    end if
    
    read(nulbgst) extractDate, vco_1Dvar % nLev_T,  vco_1Dvar %nLev_M,  vco_1Dvar %Vcode, &
         vco_1Dvar %ip1_sfc, vco_1Dvar %ip1_T_2m, vco_1Dvar %ip1_M_10m, var1D_varCount, nkgdim, nLonLatPosLand
    
    allocate( vco_1Dvar %ip1_T(nLev_T), vco_1Dvar %ip1_M(nLev_M) )
    allocate( var1D_varList(var1D_varCount) )
    allocate( bMatrix(nkgdim,nkgdim) )
    allocate( latLand(nLonLatPosLand), lonLand(nLonLatPosLand))
    read(nulbgst) vco_1Dvar %ip1_T(:), vco_1Dvar %ip1_M(:), var1D_varList(:)
    do locationIndex = 1, nLonLatPosLand
      read(nulbgst) latLand(locationIndex), lonLand(locationIndex),  bMatrix(:,:)
      if (locationIndex == 1) then
        
        nkgdim = nkgdim
       
        allocate( bHalfLand(nLonLatPosLand,nkgdim,nkgdim) )

      end if
      
      bHalfLand(locationIndex,:,:) = bMatrix(1:nkgdim,1:nkgdim)
         
      !print *,"avant utl_matsqrt", locationIndex, latLand(locationIndex), lonLand(locationIndex)
      call utl_matsqrt(bHalfLand(locationIndex,:,:), nkgdim, 1.d0, printInformation_opt=.false. )
    end do
    ierr = fclos(nulbgst)

    inquire(file=trim(oneDBmatSea), exist=fileExists)
    if ( fileExists ) then
      ierr = fnom(nulbgst, trim(oneDBmatSea), 'FTN+SEQ+UNF+OLD+R/O', 0)
    else
      call utl_abort('BHI1D_setup:No 1DVar BACKGROUND STAT FILE ' // trim(oneDBmatSea))
    end if
    read(nulbgst) extractDate, vco_1Dvar %nLev_T,  vco_1Dvar %nLev_M,  vco_1Dvar %Vcode, &
         vco_1Dvar %ip1_sfc, vco_1Dvar %ip1_T_2m, vco_1Dvar %ip1_M_10m, var1D_varCount, nkgdim, nLonLatPosSea
    
    allocate( bHalfSea(nLonLatPosSea,       nkgdim,nkgdim) )
    allocate( latSea(nLonLatPosSea), lonSea(nLonLatPosSea))
    read(nulbgst) vco_1Dvar %ip1_T(:), vco_1Dvar %ip1_M(:), var1D_varList(:)
    
    do locationIndex = 1, nLonLatPosSea
      read(nulbgst) latSea(locationIndex), lonSea(locationIndex),  bMatrix(:,:)
      bHalfSea(locationIndex,:,:) = bMatrix(1:nkgdim, 1:nkgdim)
      call utl_matsqrt(bHalfSea(locationIndex,:,:), nkgdim, 1.d0, printInformation_opt=.false. )
    end do
    ierr = fclos(nulbgst)
    vco_1Dvar % initialized = .true.
    vco_1Dvar % vGridPresent = .false.
    vco_file => vco_1Dvar
    
    
    if (.not. vco_equal(vco_anl,vco_file)) then
      call utl_abort('bmatrixHI: vco from analysisgrid and cov file do not match')
    end if
    
    status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)
    if(Vcode_anl .ne. 5002 .and. Vcode_anl .ne. 5005) then
      write(*,*) 'Vcode_anl = ',Vcode_anl
      call utl_abort('bmatrixHI: unknown vertical coordinate type!')
    end if
    
    if (.not. (gsv_varExist(varName='TT') .and.  &
               gsv_varExist(varName='UU') .and.  &
               gsv_varExist(varName='VV') .and.  &
               (gsv_varExist(varName='HU').or.gsv_varExist(varName='LQ')) .and.  &
               gsv_varExist(varName='P0')) ) then
      call utl_abort('bmatrixHI: Some or all weather fields are missing. If it is desired to deactivate the weather assimilation, then all entries of the array SCALEFACTOR in the namelist NAMBHI should be set to zero.')
    end if
    if (.not. gsv_varExist(varName='TG')) then
      write(*,*) 'bmatrixHI: WARNING: The TG field is missing. This must be present when assimilating'
      write(*,*) '                    radiance observations'
    end if

    do jlev = 1, max(nLev_M,nLev_T)
      if(scaleFactorLQ(jlev).gt.0.0d0) then 
        scaleFactorLQ(jlev) = sqrt(scaleFactorLQ(jlev))
      else
        scaleFactorLQ(jlev) = 0.0d0
      end if
    end do
      
    do jlev = 1, max(nLev_M,nLev_T)
      if(scaleFactorCC(jlev).gt.0.0d0) then 
        scaleFactorCC(jlev) = sqrt(scaleFactorCC(jlev))
      else
        scaleFactorCC(jlev) = 0.0d0
      end if
    end do

    !we want to count how many obs are really assimilable to minimize controlvector size
    var1D_obsCount = 0
    allocate(  var1D_obsPointer(obs_numHeader(obsdat)) )
    do headerIndex = 1, obs_numHeader(obsdat)
      bodyStart = obs_headElem_i(obsdat, OBS_RLN,headerIndex)
      bodyEnd= obs_headElem_i(obsdat, OBS_NLV,headerIndex) + bodyStart - 1
      countBad = 0
      do bodyIndex = bodyStart, bodyEnd
        if (obs_bodyElem_i(obsdat, OBS_ASS, bodyIndex) == obs_notAssimilated) countBad = countBad + 1
      end do
      if (countBad <  obs_headElem_i(obsdat, OBS_NLV,headerIndex) )  then
        var1D_obsCount = var1D_obsCount + 1
        var1D_obsPointer(var1D_obsCount) = headerIndex
        if (var1D_obsCount == 1) print *,'first OBS',headerIndex
      end if
    end do
    print *,'var1D_obsCount',var1D_obsCount, obs_numHeader(obsdat)

    cvDim_out = nkgdim * var1D_obsCount
    cvDim_mpilocal = cvDim_out
    initialized = .true.

  end subroutine var1D_setup


!  subroutine var1D_getScaleFactor(scaleFactor_out)
!    implicit none
!    real(8) :: scaleFactor_out(:)
!    integer :: jlev
!
!    do jlev = 1, max(nLev_M,nLev_T)
!      scaleFactor_out(jlev) = scaleFactor(jlev)
!    end do
!
!  end subroutine var1D_getScaleFactor


  
  subroutine bSqrt(controlVector_in, column, dataObs)
    implicit none

    ! Arguments
    real(8), intent(in)        :: controlVector_in(cvDim_mpilocal)
    type(struct_columnData)    :: column
    type (struct_obs):: dataObs

    ! Locals
    integer :: headerIndex, latitudeBandIndex(1), varIndex, obsIndex
    real(8), pointer :: myColumn(:)
    real(8), allocatable ::  myVector(:)
    real(8) :: latitude
    integer :: surfaceType, myOffset

    if (mpi_myid == 0) write(*,*) 'bsqrt: starting'
    if (mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if (.not. initialized) then
      if (mpi_myid == 0) write(*,*) '1Dvar B matrix not initialized'
      return
    end if
    allocate(myVector(nkgdim))
    do obsIndex = 1, var1D_obsCount 
      headerIndex = var1D_obsPointer(obsIndex)
      latitude = obs_headElem_r(dataObs, OBS_LAT, headerIndex) !radian 
      surfaceType = obs_headElem_i(dataObs, OBS_STYP, headerIndex)
      !    extract land/sea/sea-ice flag (0=land, 1=sea, 2=sea-ice)
      !       profiles(tovsIndex) % skin % surftype = tvs_ChangedStypValue(obsSpaceData,headerIndex)
      ! tvs_ChangedStypValue could be moved to another module (no real dependency with TOVS)
      if (surfaceType == 1) then !Sea
        latitudeBandIndex = minloc( abs( latitude - latSea(:)) )
        !Could be replaced later with a call to an optimized Matrix-Vector Product routine from Lapack for efficiency
        myVector(:) = matmul(bHalfSea(latitudeBandIndex(1),:,:), controlVector_in(1+(obsIndex-1)*nkgdim:obsIndex*nkgdim))
      else ! Land or Sea Ice
        latitudeBandIndex = minloc( abs( latitude - latLand(:)) )
        !Could be replaced later with a call to an optimized Matrix-Vector Product routine from Lapack for efficiency
        myVector(:) = matmul(bHalfLand(latitudeBandIndex(1),:,:), controlVector_in(1+(obsIndex-1)*nkgdim:obsIndex*nkgdim))
      end if

      myOffset = 0
      do varIndex = 1, var1D_varCount
        myColumn => col_getColumn(column, headerIndex, varName_opt=var1D_varList(varIndex))
        myColumn(:) = 0.d0
        myColumn(:) = myVector(myOffset+1:myOffset+size(myColumn))
        myOffset = myOffset + size(myColumn)
      end do
      if (myOffset /= nkgdim) then
        print *,'YY', myOffset, nkgdim
        call utl_abort('size problem tl')
      end if
    end do

    deallocate(myVector)

    if (mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if (mpi_myid == 0) write(*,*) 'var1D_bsqrt: done'

  end subroutine bSqrt

  subroutine bSqrtAd(controlVector_in, column, dataObs)
    implicit none

    ! Arguments
    real(8), intent(inout)     :: controlVector_in(cvDim_mpilocal)
    type(struct_columnData)    :: column
    type (struct_obs)          :: dataObs

    ! Locals
    integer :: headerIndex, latitudeBandIndex(1), varIndex, obsIndex
    real(8), pointer :: myColumn(:)
    real(8), allocatable ::  myVector(:)
    real(8) :: latitude
    integer :: surfaceType, myOffset

    if (mpi_myid == 0) write(*,*) 'bsqrtAd: starting'
    if (mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if (.not. initialized) then
      if (mpi_myid == 0) write(*,*) '1dvar Bmatrix not initialized'
      return
    end if

    allocate(myVector(nkgdim))
    do obsIndex = 1, var1D_obsCount
      headerIndex = var1D_obsPointer(obsIndex)
      myOffset = 0
      do varIndex = 1, var1D_varCount
        myColumn => col_getColumn(column, headerIndex, varName_opt=var1D_varList(varIndex))
        myVector(myOffset+1:myOffset+size(myColumn)) = myColumn(:)
        myOffset = myOffset + size(myColumn)
      end do
      if (myOffset /= nkgdim) then
        print *,'XX', myOffset, nkgdim
        call utl_abort('size problem ad')
      end if
      latitude = obs_headElem_r(dataObs, OBS_LAT, headerIndex) !radian
      surfaceType = obs_headElem_i(dataObs, OBS_STYP, headerIndex)
      !    extract land/sea/sea-ice flag (0=land, 1=sea, 2=sea-ice)
      !       profiles(tovsIndex) % skin % surftype = tvs_ChangedStypValue(obsSpaceData,headerIndex)
      ! tvs_ChangedStypValue could be moved to another module (no real dependency with TOVS)
      if (surfaceType == 1) then !Sea
        latitudeBandIndex = minloc( abs( latitude - latSea(:)) )
        !Could be replaced later with a call to an optimized Matrix-Vector Product routine from Blas like dsymv for efficiency
        controlVector_in(1+(obsIndex-1)*nkgdim:obsIndex*nkgdim) =  &
             controlVector_in(1+(obsIndex-1)*nkgdim:obsIndex*nkgdim) + &
             matmul(bHalfSea(latitudeBandIndex(1),:,:), myVector)
      else ! Land or Sea Ice
        latitudeBandIndex = minloc( abs( latitude - latLand(:)) )
        !Could be replaced later with a call to an optimized Matrix-Vector Product routine from Blas for efficiency
        controlVector_in(1+(obsIndex-1)*nkgdim:obsIndex*nkgdim) =  &
             controlVector_in(1+(obsIndex-1)*nkgdim:obsIndex*nkgdim) + &
             matmul(bHalfLand(latitudeBandIndex(1),:,:), myVector)
      end if

    end do
    deallocate(myVector)
    if (mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if (mpi_myid == 0) write(*,*) 'bsqrtAd: done'

  end subroutine bSqrtAd

  subroutine var1D_Finalize()
    implicit none

    if (initialized) then
       deallocate( bMatrix)
       deallocate( bHalfLand )
       deallocate( bHalfSea )
       deallocate(latLand, lonLand, latSea, lonSea )
       deallocate( var1D_obsPointer )
    end if

  end subroutine var1D_Finalize

  subroutine var1D_transferColumnToYGrid( stateVector, obsSpaceData, column)
    ! :Purpose: to transfer content of a columndata object to a statevector object
    !           without interpolation (to be used in 1DVar mode to write increments on Y grid).
    !
    implicit none
    ! arguments
    type(struct_gsv)           :: stateVector
    type(struct_obs)           :: obsSpaceData
    type(struct_columnData)    :: column
    ! locals
    integer :: obsIndex, varIndex, searchIndex, i
    integer :: startIndex
    integer, allocatable :: targetIndex(:)
    real(8), pointer :: myColumn(:), myField(:,:,:)
   
    type(struct_hco), pointer       :: hco_myGrid
    integer :: nObs1DVarTotal,  nObs1DVarMax, ierr
    integer, allocatable :: tempo(:), obsPointerMpiGlobal(:)
    real(8) :: lat, lon, latMpiGlobal, lonMpiGlobal
    real(8),allocatable,target :: dummy(:)
    real(8),allocatable :: myColumnMpiGlobal(:)
    integer :: varDim, numStep
    integer, allocatable :: dateStampList(:)


    call rpn_comm_barrier("GRID",ierr)

    call rpn_comm_allreduce(var1D_obsCount, nobs1DVarTotal, 1, "mpi_integer", "mpi_sum", "GRID", ierr)
    call rpn_comm_allreduce(var1D_obsCount, nobs1DVarMax, 1, "mpi_integer", "mpi_max", "GRID", ierr)
    if (mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    allocate( obsPointerMpiGlobal(nobs1DVarMax * mpi_nprocs),stat=ierr )
    if (mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    allocate( tempo(nobs1DVarMax),stat=ierr )
    if (mpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    tempo(:) = 0
    do obsIndex = 1, var1D_obsCount
      tempo(obsIndex) = obs_headPrimaryKey(obsSpaceData, var1D_obsPointer(obsIndex) )
    end do

    call rpn_comm_gather(tempo, nobs1DVarMax, 'MPI_INTEGER', obsPointerMpiGlobal, &
         nobs1DVarMax, 'MPI_INTEGER', 0, 'GRID', ierr)

    if ( mpi_myid == 0 ) then
      call isort(obsPointerMpiGlobal, mpi_nprocs*nobs1DVarMax)
      do i=1, mpi_nprocs * nobs1DvarMax 
        if (obsPointerMpiGlobal(i) > 0) then
          startIndex = i
          exit
        end if
      end do
    end if
       
    call rpn_comm_bcast(startIndex, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(obsPointerMpiGlobal, size(obsPointerMpiGlobal), 'MPI_INTEGER', 0, 'GRID', ierr)
    print *,"startindex",startIndex
   
    if (mpi_myId == 0 ) then
      hco_yGrid %initialized = .true.
      hco_yGrid %ni = 1
      hco_yGrid %nj = nObs1DVarTotal
      hco_yGrid %grtyp='Y'
      hco_yGrid %grtypTicTac='L'
      if (allocated(hco_yGrid %lat2d_4) ) then
        deallocate( hco_yGrid %lat2d_4)
        deallocate( hco_yGrid %lon2d_4) 
      end if
      allocate( hco_yGrid %lat2d_4(1,nObs1DvarTotal))
      allocate( hco_yGrid %lon2d_4(1,nObs1DVarTotal)) 
      hco_yGrid %xlat1 = 0.d0
      hco_yGrid %xlon1 = 0.d0
      hco_yGrid %xlat2 = 1.d0 
      hco_yGrid %xlon2 = 1.d0

      hco_myGrid => hco_ygrid
      !numStep = tim_nstepobsinc
      !allocate(dateStampList(numStep))
      numStep = 1
      !call tim_getstamplist(dateStampList,numStep,tim_getDatestamp())
      call gsv_allocate(stateVector, numstep=tim_nstepobsinc, hco_ptr=hco_myGrid, vco_ptr=column%vco, &
           datestamp_opt=tim_getDatestamp(), mpi_local_opt=.false., &
           dataKind_opt=pre_incrReal, allocHeight_opt=.false., allocPressure_opt=.false., &
           besilent_opt=.false.) !, dateStampList_opt=dateStampList )
      !deallocate(dateStampList)
      
    end if

    allocate( targetIndex(size(obsPointerMpiGlobal)) )
    targetIndex(:) = -1
    do i = startIndex, size(obsPointerMpiGlobal)
      search:do searchIndex = 1, var1D_obsCount
        if (obsPointerMpiGlobal(i) == tempo(searchIndex) ) then
          targetIndex(i) = searchIndex
          exit search
        end if
      end do search
      lat = huge(lat)
      lon = huge(lon)
      if (targetIndex(i) > 0) then
        lat = obs_headElem_r(obsSpaceData, OBS_LAT, targetIndex(i))
        lon = obs_headElem_r(obsSpaceData, OBS_LON, targetIndex(i))
      end if
      call rpn_comm_allreduce(lat, latMpiGLobal, 1, "mpi_real8", "mpi_min", "GRID", ierr)
      call rpn_comm_allreduce(lon, lonMpiGlobal, 1, "mpi_real8", "mpi_min", "GRID", ierr)
      if (mpi_myId == 0 ) then
        hco_yGrid %lat2d_4(1,i - startIndex + 1) = latMpiGLobal
        hco_yGrid %lon2d_4(1,i - startIndex + 1) = lonMpiGlobal
      end if
    end do

    do varIndex = 1, var1D_varCount      
      if (mpi_myId == 0 ) then
        call gsv_getField(stateVector, myField, varName_opt=var1D_varList(varIndex), stepIndex_opt=1)
        varDim = gsv_getNumLevFromVarName(stateVector, var1D_varList(varIndex))
      end if
      call rpn_comm_bcast(varDim, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
      allocate(dummy(varDim))
      allocate(myColumnMpiGLobal(varDim))
      dummy(:) = huge( dummy(1) )
      do i = startIndex, size(obsPointerMpiGlobal)
        if (targetIndex(i) >0 ) then
          myColumn => col_getColumn(column, targetIndex(i), varName_opt=var1D_varList(varIndex)) 
        else
          myColumn => dummy
        end if
        call rpn_comm_allreduce(myColumn, myColumnMpiGLobal, varDim, "mpi_real8", "mpi_min", "GRID", ierr)
        if (mpi_myId == 0 ) then
          myField(1, i-startIndex+1, :) = myColumnMpiGLobal(:)
        end if
      end do
      deallocate(dummy)
      deallocate(myColumnMpiGLobal)
    end do

    deallocate(obsPointerMpiGlobal)
    deallocate(targetIndex)
    deallocate(tempo)
     
  end subroutine var1D_transferColumnToYGrid

  !--------------------------------------------------------------------------
  ! var1D_get1DVarIncrement
  !--------------------------------------------------------------------------
  subroutine var1D_get1DVarIncrement(incr_cv,column,columng,obsSpaceData,nvadim_mpilocal)

    implicit none

    ! arguments
    real(8) :: incr_cv(:)
    type(struct_columnData), intent(inout)  :: column
    type(struct_columnData), intent(in)     :: columng
    type(struct_obs),        intent(in)     :: obsSpaceData
    integer :: nvadim_mpilocal

    ! compute increment from control vector (multiply by B^1/2)
    call var1D_sqrtB(incr_cv, nvadim_mpilocal, column, obsSpaceData)
    call cvt_transform(column, columng, 'PsfcToP_tl')

  end subroutine var1D_get1DVarIncrement


  !--------------------------------------------------------------------------
  ! bmat_setup
  !--------------------------------------------------------------------------
  subroutine var1D_bsetup(vco_anl, obsdat)
    !
    !:Purpose: To initialize the 1Dvar analysis Background term.
    !
    implicit none

    ! Arguments:
    type(struct_vco), pointer :: vco_anl
    type (struct_obs)         :: obsdat

    ! Locals:
    integer, allocatable :: cvDimPerInstance(:)
    integer :: cvdim
    integer :: masterBmatIndex, bMatInstanceIndex, nBmatInstance, bmatIndex

    character(len=2) :: bMatInstanceIndexString
    character(len=3) :: bMatExtraLabel

    logical :: active

    !
    !- 2.  Setup the B matrices
    !
    numBmat = 0

    do masterBmatIndex = 1, numMasterBmat

      select case( trim(masterBmatTypeList(masterBmatIndex)) )
      case ('HI')

        !- 2.1 Time-Mean Homogeneous and Isotropic...
        nBmatInstance = 1 ! hardwired
        allocate(cvdimPerInstance(nBmatInstance))

        
        write(*,*) 'Setting up the modular GLOBAL HI 1D covariances...'
        call var1D_Setup(vco_anl, obsdat, cvdim)
        
        print *,"cvdim= ",cvdim
        cvdimPerInstance(1) = cvdim


      case default

        call utl_abort( 'bmat_setup: requested bmatrix type does not exist ' // trim(masterBmatTypeList(masterBmatIndex)) )

      end select

      !- 2.6 Append the info to the B matrix info arrays and setup the proper control sub-vectors
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
    !- 3. Print a summary and set the active B matrices array
    !
    write(*,*)
    write(*,*) " var1D_setup SUMMARY, number of B matrices found = ", numBmat
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
  ! bmat_sqrtB
  !-------------------------------------------------------------------------- 
  subroutine var1D_sqrtB(controlVector, cvdim, column, obsSpaceData)
    !
    !:Purpose: To transform model state from control-vector space to grid-point
    !          space.
    !          
    !
    implicit none

    ! arguments
    integer                    :: cvdim
    real(8)                    :: controlVector(cvdim)
    type(struct_columnData)    :: column
    type(struct_obs)           :: obsSpaceData

    ! locals
    integer :: bmatIndex
    real(8), pointer :: subVector(:)
   

    !
    !- 1.  Set analysis increment to zero and allocate a temporary statevector
    !


    !
    !- 2.  Compute the analysis increment
    !
    bmat_loop: do bmatIndex = 1, numBmat

      if ( .not. bmatActive(bmatIndex) ) cycle bmat_loop

      subVector => cvm_getSubVector( controlVector, bmatLabelList(bmatIndex) )
    

      select case( trim(bmatTypeList(bmatIndex)) )
      case ('HI')

        !- 2.1 Time-Mean Homogeneous and Isotropic...
        call tmg_start(50,'B_HI')
        call bsqrt( subVector,   &  ! IN
                    column,      &  ! OUT
                    obsspacedata )  ! IN
        call tmg_stop(50)

      end select

      ! Make latest increment contribution 4D, if necessary
      !if ( bmatIs3dList(bmatIndex) ) call gsv_3dto4d( statevector_temp )

      ! Add latest contribution to total increment in statevector
      !call gsv_add( statevector_temp, statevector )

    end do bmat_loop

    !call gsv_deallocate( statevector_temp )

  end subroutine var1D_sqrtB

  !--------------------------------------------------------------------------
  ! bmat_sqrtBT
  !--------------------------------------------------------------------------
  subroutine var1D_sqrtBT(controlVector, cvdim, column, obsData)
    !
    !:Purpose: To transform model state from grid-point space to
    !          error-covariance space.
    !
    implicit none

    ! Arguments
    integer :: cvdim
    real(8) :: controlVector(cvdim)
    type(struct_columnData)    :: column
    type(struct_obs)           :: obsData

    ! Locals
    integer :: bmatIndex
    real(8),pointer :: subVector(:)


    ! Process components in opposite order as forward calculation
    bmat_loop: do bmatIndex = numBmat, 1, -1

      if ( .not. bmatActive(bmatIndex) ) cycle bmat_loop

      subVector => cvm_getSubVector( controlVector, bmatLabelList(bmatIndex) )
     

      select case( trim(bmatTypeList(bmatIndex)) )

      case ('HI')

        !- 2.5 Time-Mean Homogeneous and Isotropic...
        call tmg_start(51,'B_HI_T')
    
        call bsqrtad( subvector, &  ! IN
                      column,    &  ! OUT
                      obsData )     ! IN

        call tmg_stop(51)

      end select

    end do bmat_loop

  end subroutine var1D_sqrtBT

end module var1D_mod
