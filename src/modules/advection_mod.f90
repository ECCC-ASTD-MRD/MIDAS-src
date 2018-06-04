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

!--------------------------------------------------------------------------
!! MODULE advection (prefix="adv")
!!
!! *Purpose*: To perform forward and/or backward advection (based on 
!!            semi-lagrangian trajectories) for both gridStateVector and
!!            ensemble of gridStateVectors
!!
!--------------------------------------------------------------------------
MODULE advection_mod
  use ramDisk_mod
  use mpivar_mod
  use mathPhysConstants_mod
  use EarthConstants_mod
  use timeCoord_mod
  use ensembleStateVector_mod
  use gridStateVector_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use utilities_mod
  implicit none
  save
  private

  ! public procedures
  public :: struct_adv, adv_Setup
  public :: adv_ensemble_tl, adv_ensemble_ad
  public :: adv_statevector_tl, adv_statevector_ad

  type :: struct_adv_LevType
    integer, allocatable :: lonIndex       (:,:,:) ! lon, lat, lev
    integer, allocatable :: latIndex       (:,:,:)
    real(8), allocatable :: interpWeight_BL(:,:,:)
    real(8), allocatable :: interpWeight_BR(:,:,:)
    real(8), allocatable :: interpWeight_TL(:,:,:)
    real(8), allocatable :: interpWeight_TR(:,:,:)
  end type struct_adv_LevType

  type :: struct_adv_timeStep
    type(struct_adv_LevType), allocatable :: levType(:) 
  end type struct_adv_timeStep

  type :: struct_adv
    private
    integer :: nLev_M
    integer :: nLev_T
    integer :: ni, nj
    integer :: lonPerPE, lonPerPEmax, myLonBeg, myLonEnd
    integer :: latPerPE, latPerPEmax, myLatBeg, myLatEnd
    integer, allocatable :: allLonBeg(:), allLatBeg(:)
    integer :: nTimeStep
    integer :: timeStepIndexMainSource
    integer, allocatable :: timeStepIndexSource(:)
    logical :: singleTimeStepIndexSource
    type(struct_adv_timeStep), allocatable :: timeStep(:)
  end type struct_adv

  integer, external  :: get_max_rss

  integer, parameter :: MMindex = 1
  integer, parameter :: THindex = 2
  integer, parameter :: SFindex = 3

  real(8), parameter :: numGridPts    =  1.0d0 ! used to compute numSubStep
  real(8), parameter :: latitudePatch = 80.0d0 ! this defines latitude where rotated grid used

  integer :: numStepReferenceFlow
  integer :: lonPerPE, latPerPE
  integer, allocatable :: allLonBeg(:), allLatBeg(:)
  integer, allocatable :: numSubStep(:)

  real(8), pointer     :: uu_referenceFlow_ptr4d(:,:,:,:)
  real(8), pointer     :: vv_referenceFlow_ptr4d(:,:,:,:)
  real(8), allocatable :: uu_referenceFlow_mpiGlobal(:,:,:)
  real(8), allocatable :: vv_referenceFlow_mpiGlobal(:,:,:)
  real(8), allocatable :: uu_referenceFlow_ThermoLevel(:,:)
  real(8), allocatable :: vv_referenceFlow_ThermoLevel(:,:)

  real(8) :: delT_sec, referenceFlowFactor

  type(struct_hco), pointer :: hco

  ! Control parameter for the level of listing output
  logical, parameter :: verbose = .false.

CONTAINS

  !--------------------------------------------------------------------------
  ! adv_Setup
  !--------------------------------------------------------------------------
  SUBROUTINE adv_setup(adv, mode, hco_in, vco_in, numStepAdvectedField, &
                       dateStampListAdvectedField, numStepReferenceFlow_in, delT_hour, &
                       referenceFlowFactor_in, levTypeList, referenceFlowFilename_opt, &
                       statevector_referenceFlow_opt)
    implicit none

    type(struct_adv) :: adv
    type(struct_hco), pointer :: hco_in
    type(struct_vco), pointer :: vco_in

    character(len=*), intent(in) :: mode
    character(len=*), intent(in) :: levTypeList
    character(len=*), optional, intent(in) :: referenceFlowFilename_opt
    integer, intent(in) :: numStepAdvectedField, numStepReferenceFlow_in
    integer, intent(in) :: dateStampListAdvectedField(numStepAdvectedField)
    real(8), intent(in) :: referenceFlowFactor_in, delT_hour

    type(struct_gsv), optional :: statevector_referenceFlow_opt

    integer :: latIndex0, lonIndex0, latIndex, lonIndex, levIndex, jsubStep, stepIndexRF, stepIndexAF, ierr, gdxyfll
    integer :: nsize, latIndex_mpiglobal, lonIndex_mpiglobal
    integer :: alfa, nLevType

    integer, allocatable :: dateStampListReferenceFlow(:)
    integer, allocatable :: advectedFieldAssociatedStepIndexRF(:)
    integer, allocatable :: advectionReferenceFlowStartingStepIndex(:)
    integer, allocatable :: advectionReferenceFlowEndingStepIndex  (:)

    real(8) :: uu, vv, subDelT, lonAdvect, latAdvect, delx, dely, sumWeight
    real(8) :: uu_p, vv_p, lonAdvect_p, latAdvect_p, Gcoef, Scoef
    real(8) :: interpWeight_BL, interpWeight_BR, interpWeight_TL, interpWeight_TR
    real(8), allocatable :: uu_referenceFlow_mpiGlobalTiles(:,:,:,:)
    real(8), allocatable :: vv_referenceFlow_mpiGlobalTiles(:,:,:,:)

    real(4) :: lonAdvect_deg_r4, latAdvect_deg_r4, xpos_r4, ypos_r4

    character(len=64) :: filename
    character(len=3)  :: filenumber

    type(struct_gsv) :: statevector_referenceFlow

    logical :: AdvectFileExists

    integer :: nLev, kIndex, levTypeIndex, stepIndexRF_start, stepIndexRF_end 
    integer :: myLonBeg, myLonEnd
    integer :: myLatBeg, myLatEnd

    !
    !- 1.  Set low-level variables
    !
    numStepReferenceFlow = numStepReferenceFlow_in
    referenceFlowFactor  = referenceFlowFactor_in
    adv%nTimeStep        = numStepAdvectedField
    
    allocate(adv%timeStepIndexSource(numStepAdvectedField))

    if (vco_in%Vcode /= 5002 ) then
      call utl_abort('adv_setup: Only vCode=5002 is currently supported!')
    end if

    select case(trim(levTypeList))
    case ('MMLevsOnly')
      nLevType = 1
    case ('allLevs')
      nLevType = 3
    case default
      write(*,*)
      write(*,*) 'Unsupported levTypeList: ', trim(levTypeList)
      call utl_abort('adv_setup')
    end select
    
    !- 1.1 Mode
    select case(trim(mode))
    case ('fromFirstTimeIndex')
      adv%timeStepIndexMainSource  = 1
      adv%timeStepIndexSource(:)   = adv%timeStepIndexMainSource
      adv%singleTimeStepIndexSource= .true.
    case ('fromMiddleTimeIndex')
      if (mod(numStepAdvectedField,2) == 0) then
        call utl_abort('adv_setup: numStepAdvectedField cannot be even with direction=fromMiddleTimeIndex') 
      end if
      adv%timeStepIndexMainSource  = (numStepAdvectedField+1)/2
      adv%timeStepIndexSource(:)   = adv%timeStepIndexMainSource
      adv%singleTimeStepIndexSource= .true.
    case ('towardFirstTimeIndex','towardFirstTimeIndexInverse')
      adv%timeStepIndexMainSource = 1
      do stepIndexAF = 1, numStepAdvectedField
        adv%timeStepIndexSource(stepIndexAF) = stepIndexAF
      end do
      adv%singleTimeStepIndexSource = .false.
    case('towardMiddleTimeIndex','towardMiddleTimeIndexInverse')
      if (mod(numStepAdvectedField,2) == 0) then
        call utl_abort('adv_setup: numStepAdvectedField cannot be even with direction=fromMiddleTimeIndex') 
      end if
      adv%timeStepIndexMainSource = (numStepAdvectedField+1)/2
      do stepIndexAF = 1, numStepAdvectedField
        adv%timeStepIndexSource(stepIndexAF) = stepIndexAF
      end do
      adv%singleTimeStepIndexSource = .false.
    case default
      write(*,*)
      write(*,*) 'Unsupported mode : ', trim(mode)
      call utl_abort('adv_setup')
    end select

    !- Set some important values
    delT_sec = delT_hour*3600.0D0

    !- 1.2 Grid Size
    hco => hco_in
    adv%ni = hco%ni
    adv%nj = hco%nj
    adv%nLev_M = vco_in%nLev_M
    adv%nLev_T = vco_in%nLev_T

    call mpivar_setup_latbands(adv%nj, adv%latPerPE, adv%latPerPEmax, adv%myLatBeg, adv%myLatEnd)
    call mpivar_setup_lonbands(adv%ni, adv%lonPerPE, adv%lonPerPEmax, adv%myLonBeg, adv%myLonEnd)
    allocate(adv%allLonBeg(mpi_npex))
    call rpn_comm_allgather(adv%myLonBeg,1,"mpi_integer",       &
         adv%allLonBeg,1,"mpi_integer","EW",ierr)
    allocate(adv%allLatBeg(mpi_npey))
    call rpn_comm_allgather(adv%myLatBeg,1,"mpi_integer",       &
         adv%allLatBeg,1,"mpi_integer","NS",ierr)

    lonPerPE = adv%lonPerPE 
    latPerPE = adv%latPerPE
    allocate(allLonBeg(mpi_npex))
    allocate(allLatBeg(mpi_npey))
    allLonBeg(:) = adv%allLonBeg(:)
    allLatBeg(:) = adv%allLatBeg(:)

    !- 1.3 Memory allocation
    myLonBeg = adv%myLonBeg
    myLonEnd = adv%myLonEnd
    myLatBeg = adv%myLatBeg
    myLatEnd = adv%myLatEnd

    nLev = adv%nLev_M

    allocate(adv%timeStep(numStepAdvectedField))

    do stepIndexAF = 1, numStepAdvectedField

      if (stepIndexAF == adv%timeStepIndexMainSource) cycle ! no interpolation needed for this time step

      allocate(adv%timeStep(stepIndexAF)%levType(nLevType))
      do levTypeIndex = 1,  nLevType ! ( 1=MM, 2=TH, 3=SF )
        if      (levTypeIndex == MMindex ) then
          nLev = adv%nLev_M
        else if (levTypeIndex == THindex) then
          nLev = adv%nLev_T
        else
          nLev = 1
        end if
        allocate(adv%timeStep(stepIndexAF)%levType(levTypeIndex)%lonIndex       (myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLev))
        allocate(adv%timeStep(stepIndexAF)%levType(levTypeIndex)%latIndex       (myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLev))
        allocate(adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_BL(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLev))
        allocate(adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_BR(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLev))
        allocate(adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_TL(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLev))
        allocate(adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_TR(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLev))
        adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_BL(:,:,:) = 1.0d0
        adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_BR(:,:,:) = 0.0d0
        adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_TL(:,:,:) = 0.0d0
        adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_TR(:,:,:) = 0.0d0
      end do
    end do

    !
    !- 2.  Read in the reference forecasts (winds) to use for advection
    !
    allocate(dateStampListReferenceFlow(numStepReferenceFlow))

    if (present(referenceFlowFilename_opt) ) then

      if (mpi_myid == 0)  then
        write(*,*)
        write(*,*) 'referenceFlow source taken from input file = ', trim(referenceFlowFilename_opt) 
      end if
      !- Read in the forecasts (winds) to use for advection
      do stepIndexRF = 1, numStepReferenceFlow
        call incdatr(dateStampListReferenceFlow(stepIndexRF), dateStampListAdvectedField(1), &
                     real(stepIndexRF-1,8)*delT_hour)
      end do

      call gsv_allocate(statevector_referenceFlow,numStepReferenceFlow, hco, vco_in, &
                        dateStampList_opt=dateStampListReferenceFlow, &
                        varNames_opt=(/'UU','VV','P0'/), mpi_local_opt=.true.)
      
      fileName = ram_fullWorkingPath(trim(referenceFlowFilename_opt))
      inquire(file=trim(fileName),exist=AdvectFileExists)
      write(*,*) 'AdvectFileExists', AdvectFileExists
      do stepIndexRF = 1, numStepReferenceFlow
        call gsv_readFromFile(statevector_referenceFlow,fileName,' ',' ',stepIndex_opt=stepIndexRF)
      end do

      uu_referenceFlow_ptr4d => gsv_getField_r8(statevector_referenceFlow, 'UU')
      vv_referenceFlow_ptr4d => gsv_getField_r8(statevector_referenceFlow, 'VV')

    else if (present(statevector_referenceFlow_opt)) then

      if (mpi_myid == 0)  then
        write(*,*)
        write(*,*) 'referenceFlow source = input gridStateVector'
        write(*,*) numStepReferenceFlow
        write(*,*) statevector_referenceFlow_opt%dateStampList(:)
      end if
      dateStampListReferenceFlow(:) = statevector_referenceFlow_opt%dateStampList(:)
      uu_referenceFlow_ptr4d => gsv_getField_r8(statevector_referenceFlow_opt, 'UU')
      vv_referenceFlow_ptr4d => gsv_getField_r8(statevector_referenceFlow_opt, 'VV')

    else
      call utl_abort('adv_setup: referenceFlow source was not provided!')
    end if

    !
    !-  3.  Advection setup
    !

    !-  3.1 Match the stepIndex between the reference flow and the fields to be advected
    allocate(advectedFieldAssociatedStepIndexRF(numStepAdvectedField))
    advectedFieldAssociatedStepIndexRF(:) = -1
    do stepIndexAF = 1, numStepAdvectedField
      do stepIndexRF = 1, numStepReferenceFlow
        if ( dateStampListAdvectedField(stepIndexAF) == dateStampListReferenceFlow(stepIndexRF) ) then
          advectedFieldAssociatedStepIndexRF(stepIndexAF) = stepIndexRF
          if (mpi_myid == 0)  then
            write(*,*)
            write(*,*) 'stepIndex Match', stepIndexAF, stepIndexRF
          end if
          exit 
        end if
      end do
      if ( advectedFieldAssociatedStepIndexRF(stepIndexAF) == -1 ) then
        call utl_abort('adv_setup: no match between dateStampListAdvectedField and dateStampListReferenceFlow')
      end if
    end do

    !-  3.2 Set starting, ending and direction parameters
    allocate(advectionReferenceFlowStartingStepIndex(numStepAdvectedField))
    allocate(advectionReferenceFlowEndingStepIndex  (numStepAdvectedField))

    select case(trim(mode))
    case ('fromFirstTimeIndex','fromMiddleTimeIndex','towardFirstTimeIndexInverse','towardMiddleTimeIndexInverse')
      do stepIndexAF = 1, numStepAdvectedField
        advectionReferenceFlowStartingStepIndex(stepIndexAF) = advectedFieldAssociatedStepIndexRF(adv%timeStepIndexMainSource)
        advectionReferenceFlowEndingStepIndex  (stepIndexAF) = advectedFieldAssociatedStepIndexRF(stepIndexAF)
      end do
    case ('towardFirstTimeIndex','towardMiddleTimeIndex')
      do stepIndexAF = 1, numStepAdvectedField
        advectionReferenceFlowStartingStepIndex(stepIndexAF) = advectedFieldAssociatedStepIndexRF(stepIndexAF)
        advectionReferenceFlowEndingStepIndex  (stepIndexAF) = advectedFieldAssociatedStepIndexRF(adv%timeStepIndexMainSource)
      end do
    case default
      write(*,*)
      write(*,*) 'Oops! This should never happen. Check the code...'
      call utl_abort('adv_setup')
    end select

    !
    !- 4.  Perform the advection (backward and/or forward) 
    !
    if (mpi_myid == 0) write(*,*) 'setupAdvectAmplitude: starting'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    allocate(numSubStep(adv%nj))
    allocate(uu_referenceFlow_mpiGlobal(numStepReferenceFlow, adv%ni, adv%nj))
    allocate(vv_referenceFlow_mpiGlobal(numStepReferenceFlow, adv%ni, adv%nj))

    allocate(uu_referenceFlow_mpiGlobalTiles(numStepReferenceFlow, adv%lonPerPE, adv%latPerPE, mpi_nprocs))
    allocate(vv_referenceFlow_mpiGlobalTiles(numStepReferenceFlow, adv%lonPerPE, adv%latPerPE, mpi_nprocs))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    do levTypeIndex = 1, min(nLevType,2) ! make sure that SF level is skipped
      if      (levTypeIndex == MMindex ) then
        nLev = adv%nLev_M
      else if (levTypeIndex == THindex) then
        nLev = adv%nLev_T
        allocate(uu_referenceFlow_ThermoLevel(myLonBeg:myLonEnd,myLatBeg:myLatEnd))
        allocate(vv_referenceFlow_ThermoLevel(myLonBeg:myLonEnd,myLatBeg:myLatEnd))
      else
        call utl_abort('adv_setup: unknown levTypeIndex')
      end if

      if (mpi_myid == 0) write(*,*)
      if (mpi_myid == 0) write(*,*) 'setupAdvectAmplitude: levTypeIndex = ', levTypeIndex
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      do levIndex = 1, nLev ! loop over levels

        if (mpi_myid == 0) write(*,*) 'setupAdvectAmplitude: levIndex = ', levIndex

        call processReferenceFlow(levTypeIndex, levIndex,                                           & ! IN
                                  uu_referenceFlow_mpiGlobalTiles, vv_referenceFlow_mpiGlobalTiles, & ! IN
                                  adv%nLev_M, adv%nLev_T, myLatBeg, myLatEnd)                         ! IN

        do stepIndexAF = 1, numStepAdvectedField
          if (stepIndexAF == adv%timeStepIndexMainSource) cycle ! no interpolation needed for this time step

          stepIndexRF_start = advectionReferenceFlowStartingStepIndex(stepIndexAF)
          stepIndexRF_end   = advectionReferenceFlowEndingStepIndex  (stepIndexAF) 

          ! loop over all initial grid points within tile for determining trajectories
          do latIndex0 = adv%myLatBeg, adv%myLatEnd
            do lonIndex0 = adv%myLonBeg, adv%myLonEnd

              call calcTrajAndWeights(lonIndex, latIndex, interpWeight_BL, interpWeight_BR,   & ! OUT
                                      interpWeight_TL, interpWeight_TR,                       & ! OUT
                                      latIndex0, lonIndex0, stepIndexRF_start, stepIndexRF_end) ! IN

              ! store the final position of the trajectory and interp weights
              adv%timeStep(stepIndexAF)%levType(levTypeIndex)%lonIndex       (lonIndex0,latIndex0,levIndex) = lonIndex
              adv%timeStep(stepIndexAF)%levType(levTypeIndex)%latIndex       (lonIndex0,latIndex0,levIndex) = latIndex
              adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_BL(lonIndex0,latIndex0,levIndex) = interpWeight_BL
              adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_BR(lonIndex0,latIndex0,levIndex) = interpWeight_BR
              adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_TL(lonIndex0,latIndex0,levIndex) = interpWeight_TL
              adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_TR(lonIndex0,latIndex0,levIndex) = interpWeight_TR

            end do ! lonIndex0
          end do ! latIndex0

        end do ! stepIndexAF

      end do ! levIndex

      if (levTypeIndex == THindex) then
        deallocate(vv_referenceFlow_ThermoLevel)
        deallocate(uu_referenceFlow_ThermoLevel)
      end if

    end do ! levTypeIndex

    ! Surface Level: Use info from lower momentum levels
    if ( nLevType == 3 ) then
      do stepIndexAF = 1, numStepAdvectedField
        if (stepIndexAF == adv%timeStepIndexMainSource) cycle ! no interpolation needed for this time step
        adv%timeStep(stepIndexAF)%levType(SFindex)%lonIndex       (:,:,1) = &
             adv%timeStep(stepIndexAF)%levType(MMindex)%lonIndex       (:,:,adv%nLev_M)
        adv%timeStep(stepIndexAF)%levType(SFindex)%latIndex       (:,:,1) = &
             adv%timeStep(stepIndexAF)%levType(MMindex)%latIndex       (:,:,adv%nLev_M)
        adv%timeStep(stepIndexAF)%levType(SFindex)%interpWeight_BL(:,:,1) = &
             adv%timeStep(stepIndexAF)%levType(MMindex)%interpWeight_BL(:,:,adv%nLev_M)
        adv%timeStep(stepIndexAF)%levType(SFindex)%interpWeight_BR(:,:,1) = &
             adv%timeStep(stepIndexAF)%levType(MMindex)%interpWeight_BR(:,:,adv%nLev_M)
        adv%timeStep(stepIndexAF)%levType(SFindex)%interpWeight_TL(:,:,1) = &
             adv%timeStep(stepIndexAF)%levType(MMindex)%interpWeight_TL(:,:,adv%nLev_M)
        adv%timeStep(stepIndexAF)%levType(SFindex)%interpWeight_TR(:,:,1) = &
             adv%timeStep(stepIndexAF)%levType(MMindex)%interpWeight_TR(:,:,adv%nLev_M)
      end do
    end if

    deallocate(numSubStep)
    deallocate(allLonBeg)
    deallocate(allLatBeg)
    deallocate(uu_referenceFlow_mpiGlobalTiles)
    deallocate(vv_referenceFlow_mpiGlobalTiles)
    deallocate(uu_referenceFlow_mpiGlobal)
    deallocate(vv_referenceFlow_mpiGlobal)
    if (present(referenceFlowFilename_opt) ) then
      call gsv_deallocate(statevector_referenceFlow)
    end if

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if (mpi_myid == 0) write(*,*) 'adv_setup: done'

  end SUBROUTINE adv_setup

  !--------------------------------------------------------------------------
  ! processReferenceFlow
  !--------------------------------------------------------------------------
  SUBROUTINE processReferenceFlow (levTypeIndex, levIndex, &
                                   uu_referenceFlow_mpiGlobalTiles, vv_referenceFlow_mpiGlobalTiles, &
                                   nLev_M, nLev_T, myLatBeg, myLatEnd)
    implicit none
    integer, intent(in) :: levTypeIndex, levIndex, nLev_M, nLev_T, myLatBeg, myLatEnd
    real(8) :: uu_referenceFlow_mpiGlobalTiles(:,:,:,:)
    real(8) :: vv_referenceFlow_mpiGlobalTiles(:,:,:,:)

    integer :: stepIndexRF, nsize, ierr 
    integer :: procID, procIDx, procIDy, lonIndex, latIndex
    integer :: lonIndex_mpiglobal, latIndex_mpiglobal

    real(8) :: uu, vv, latAdvect

    nsize = lonPerPE*latPerPE

    if ( levTypeIndex == MMindex ) then
      ! No vertical interpolation is needed
      do stepIndexRF = 1, numStepReferenceFlow
        ! gather the winds for this level
        call rpn_comm_allgather(uu_referenceFlow_ptr4d(:,:,levIndex,stepIndexRF)  , nsize, "mpi_double_precision", &
                                uu_referenceFlow_mpiGlobalTiles(stepIndexRF,:,:,:), nsize, "mpi_double_precision", &
                                "GRID", ierr )
        call rpn_comm_allgather(vv_referenceFlow_ptr4d(:,:,levIndex,stepIndexRF)  , nsize, "mpi_double_precision", &
                                vv_referenceFlow_mpiGlobalTiles(stepIndexRF,:,:,:), nsize, "mpi_double_precision", &
                                "GRID", ierr )
      end do
    else if (levTypeIndex == THindex ) then
      ! Vertical interpolation is needed...
      ! The adopted approach follows the vertical interpolation for the amplitude fields in bMatrixEnsemble_mod
      do stepIndexRF = 1, numStepReferenceFlow
        if (levIndex == 1) then
          ! use top momentum level amplitudes for top thermo level
          uu_referenceFlow_ThermoLevel(:,:) = uu_referenceFlow_ptr4d(:,:,1,stepIndexRF)
          vv_referenceFlow_ThermoLevel(:,:) = vv_referenceFlow_ptr4d(:,:,1,stepIndexRF)
        else if (levIndex == nLev_T) then
          ! use surface momentum level amplitudes for surface thermo level
          uu_referenceFlow_ThermoLevel(:,:) = uu_referenceFlow_ptr4d(:,:,nLev_M,stepIndexRF)
          vv_referenceFlow_ThermoLevel(:,:) = vv_referenceFlow_ptr4d(:,:,nLev_M,stepIndexRF)
        else
          ! for other levels, interpolate momentum winds to get thermo winds
!$OMP PARALLEL DO PRIVATE (latIndex)
          do latIndex = myLatBeg, myLatEnd
            uu_referenceFlow_ThermoLevel(:,latIndex) = 0.5d0*( uu_referenceFlow_ptr4d(:,latIndex,levIndex-1,stepIndexRF) + &
                 uu_referenceFlow_ptr4d(:,latIndex,levIndex,stepIndexRF) )
            vv_referenceFlow_ThermoLevel(:,latIndex) = 0.5d0*( vv_referenceFlow_ptr4d(:,latIndex,levIndex-1,stepIndexRF) + &
                 vv_referenceFlow_ptr4d(:,latIndex,levIndex,stepIndexRF) )
          end do
!$OMP END PARALLEL DO
        end if

        ! gather the INTERPOLATED winds for this level
        call rpn_comm_allgather(uu_referenceFlow_ThermoLevel                      , nsize, "mpi_double_precision",  &
                                uu_referenceFlow_mpiGlobalTiles(stepIndexRF,:,:,:), nsize, "mpi_double_precision",  &
                                "GRID", ierr )
        call rpn_comm_allgather(uu_referenceFlow_ThermoLevel                      , nsize, "mpi_double_precision",  &
                                vv_referenceFlow_mpiGlobalTiles(stepIndexRF,:,:,:), nsize, "mpi_double_precision",  &
                                "GRID", ierr )

      end do

    else
      call utl_abort('processReferenceFlow: invalid levTypeIndex')
    end if

    ! rearrange gathered winds for convenience
    do procIDy = 0, (mpi_npey-1)
      do procIDx = 0, (mpi_npex-1)
        procID = procIDx + procIDy*mpi_npex
        do latIndex = 1, latPerPE
          latIndex_mpiglobal = latIndex + allLatBeg(procIDy+1) - 1
          do lonIndex = 1, lonPerPE
            lonIndex_mpiglobal = lonIndex + allLonBeg(procIDx+1) - 1
            uu_referenceFlow_mpiGlobal(:, lonIndex_mpiglobal, latIndex_mpiglobal) = uu_referenceFlow_mpiGlobalTiles(:, lonIndex, latIndex, procID+1)
            vv_referenceFlow_mpiGlobal(:, lonIndex_mpiglobal, latIndex_mpiglobal) = vv_referenceFlow_mpiGlobalTiles(:, lonIndex, latIndex, procID+1)
          end do
        end do
      end do
    end do

    ! determine the number of time steps required as a function of latitude
    do latIndex = 1, hco%nj
      latAdvect = hco%lat(latIndex)
      if (abs(latAdvect) < latitudePatch*MPC_RADIANS_PER_DEGREE_R8) then
        uu = maxval(abs(uu_referenceFlow_mpiGlobal(:,:,latIndex) /(RA*cos(latAdvect)))) ! in rad/s
        vv = maxval(abs(vv_referenceFlow_mpiGlobal(:,:,latIndex) / RA)) ! in rad/s
      else
        uu = maxval(abs(uu_referenceFlow_mpiGlobal(:,:,latIndex) / RA)) ! in rad/s
        vv = maxval(abs(vv_referenceFlow_mpiGlobal(:,:,latIndex) / RA)) ! in rad/s
      end if
      numSubStep(latIndex) = max( 1,  &
           nint( (delT_sec * referenceFlowFactor * uu) / (numGridPts*(hco%lon(2)-hco%lon(1))) ),  &
           nint( (delT_sec * referenceFlowFactor * vv) / (numGridPts*(hco%lat(2)-hco%lat(1))) ) )
    end do
    if (mpi_myid == 0) write(*,*) 'min and max of numSubStep',minval(numSubStep(:)),maxval(numSubStep(:))

  end SUBROUTINE processReferenceFlow

  !--------------------------------------------------------------------------
  ! calcTrajAndWeights
  !--------------------------------------------------------------------------
  SUBROUTINE calcTrajAndWeights(lonIndex, latIndex, interpWeight_BL, interpWeight_BR, &
                                interpWeight_TL, interpWeight_TR, latIndex0, lonIndex0, &
                                stepIndexRF_start, stepIndexRF_end)
    implicit none

    integer, intent(out) :: lonIndex, latIndex
    real(8), intent(out) :: interpWeight_BL, interpWeight_BR, interpWeight_TL, interpWeight_TR
    integer, intent(in)  :: latIndex0, lonIndex0, stepIndexRF_start, stepIndexRF_end

    integer :: subStepIndex, stepIndexRF, ierr, gdxyfll
    integer :: alfa, ni, nj,  stepIndex_direction, stepIndex_first, stepIndex_last

    real(8) :: uu, vv, subDelT, lonAdvect, latAdvect, delx, dely, sumWeight
    real(8) :: uu_p, vv_p, lonAdvect_p, latAdvect_p, Gcoef, Scoef
    real(4) :: lonAdvect_deg_r4, latAdvect_deg_r4, xpos_r4, ypos_r4

    ni = hco%ni
    nj = hco%nj

    subDelT = delT_sec/real(numSubStep(latIndex0),8)  ! in seconds

    ! position at the initial time of back trajectory
    lonAdvect = hco%lon(lonIndex0)  ! in radians
    latAdvect = hco%lat(latIndex0)
    lonIndex = lonIndex0  ! index
    latIndex = latIndex0
    xpos_r4 = real(lonIndex,4)
    ypos_r4 = real(latIndex,4)

    ! initial positions in rotated coordinate system
    lonAdvect_p = 0.0d0
    latAdvect_p = 0.0d0

    if (verbose) then
      write(*,*) 'final lonAdvect,latAdvect=', &
           lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
           latAdvect*MPC_DEGREES_PER_RADIAN_R8
      write(*,*) 'numSubStep=', numSubStep(latIndex0)
    end if

    ! time stepping strategy
    if      (stepIndexRF_end > stepIndexRF_start) then
      ! back trajectory   , stepping backwards
      stepIndex_first     = stepIndexRF_end-1
      stepIndex_last      = stepIndexRF_start
      stepIndex_direction = -1
    else if (stepIndexRF_end < stepIndexRF_start) then
      ! forward trajectory, stepping forward
      stepIndex_first     = stepIndexRF_end
      stepIndex_last      = stepIndexRF_start-1
      stepIndex_direction = 1
    else
      call utl_abort('calcTrajAndWeights: fatal error with stepIndexRF')
    end if

    do stepIndexRF = stepIndex_first, stepIndex_last, stepIndex_direction

      if (verbose) write(*,*) 'stepIndexRF,lonIndex,latIndex=',stepIndexRF,lonIndex,latIndex

      do subStepIndex = 1, numSubStep(latIndex0)

        alfa = (subStepIndex-1)/numSubStep(latIndex0)
        ! perform one timestep of back trajectory
        if (abs(hco%lat(latIndex0)) < latitudePatch*MPC_RADIANS_PER_DEGREE_R8) then
          ! points away from pole, handled normally
          ! determine wind at current location (now at BL point)
          uu = (  alfa *uu_referenceFlow_mpiGlobal(stepIndexRF  ,lonIndex,latIndex) + &
               (1-alfa)*uu_referenceFlow_mpiGlobal(stepIndexRF+1,lonIndex,latIndex) ) &
               /(RA*cos(hco%lat(latIndex))) ! in rad/s
          vv = (   alfa*vv_referenceFlow_mpiGlobal(stepIndexRF  ,lonIndex,latIndex) + &
               (1-alfa)*vv_referenceFlow_mpiGlobal(stepIndexRF+1,lonIndex,latIndex) ) &
               /RA
          ! apply user-specified scale factor to advecting winds
          uu = referenceFlowFactor * uu
          vv = referenceFlowFactor * vv

          ! compute next position
          lonAdvect = lonAdvect + real(stepIndex_direction,8)*subDelT*uu  ! in radians
          latAdvect = latAdvect + real(stepIndex_direction,8)*subDelT*vv

          if (verbose) then
            write(*,*) 'not near pole, lonAdvect,latAdvect,uu,vv=', &
                 lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                 latAdvect*MPC_DEGREES_PER_RADIAN_R8,uu,vv
          end if
        else
          ! points near pole, handled in a special way
          ! determine wind at current location (now at BL point)
          uu =    alfa *uu_referenceFlow_mpiGlobal(stepIndexRF  ,lonIndex,latIndex) + &
               (1-alfa)*uu_referenceFlow_mpiGlobal(stepIndexRF+1,lonIndex,latIndex)  ! in m/s
          vv =    alfa *vv_referenceFlow_mpiGlobal(stepIndexRF  ,lonIndex,latIndex) + &
               (1-alfa)*vv_referenceFlow_mpiGlobal(stepIndexRF+1,lonIndex,latIndex)
          ! transform wind vector into rotated coordinate system
          Gcoef = ( cos(latAdvect)*cos(hco%lat(latIndex0)) + &
               sin(latAdvect)*sin(hco%lat(latIndex0))*cos(lonAdvect-hco%lon(lonIndex0)) ) / &
               cos(latAdvect_p)
          Scoef = ( sin(hco%lat(latIndex0))*sin(lonAdvect-hco%lon(lonIndex0)) ) / &
               cos(latAdvect_p)
          uu_p = Gcoef * uu - Scoef * vv ! in m/s
          vv_p = Scoef * uu + Gcoef * vv 

          ! apply user-specified scale factor to advecting winds
          uu_p = referenceFlowFactor * uu_p ! in m/s
          vv_p = referenceFlowFactor * vv_p

          ! compute next position (in rotated coord system)
          lonAdvect_p = lonAdvect_p + real(stepIndex_direction,8)*subDelT*uu_p/(RA*cos(latAdvect_p))  ! in radians
          latAdvect_p = latAdvect_p + real(stepIndex_direction,8)*subDelT*vv_p/RA

          if (verbose) then
            write(*,*) '    near pole, uu_p,vv_p,Gcoef,Scoef=', &
                 uu_p, vv_p, Gcoef, Scoef
            write(*,*) '    near pole, lonAdvect_p,latAdvect_p=', &
                 lonAdvect_p*MPC_DEGREES_PER_RADIAN_R8, &
                 latAdvect_p*MPC_DEGREES_PER_RADIAN_R8
          end if

          ! compute lon/lat in original coordinate system
          lonAdvect = hco%lon(lonIndex0) +                                                  &
               atan2( cos(latAdvect_p)*sin(lonAdvect_p) ,                            &
               ( cos(latAdvect_p)*cos(lonAdvect_p)*cos(hco%lat(latIndex0)) -  &
               sin(latAdvect_p)*sin(hco%lat(latIndex0)) ) )
          latAdvect = asin( cos(latAdvect_p)*cos(lonAdvect_p)*sin(hco%lat(latIndex0)) + &
               sin(latAdvect_p)*cos(hco%lat(latIndex0)) )

          if (verbose) then
            write(*,*) '    near pole, lonAdvect,latAdvect=', &
                 lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                 latAdvect*MPC_DEGREES_PER_RADIAN_R8
          end if
        end if

        ! convert lon/lat position into index
        lonAdvect_deg_r4 = real(lonAdvect,4)* MPC_DEGREES_PER_RADIAN_R4
        latAdvect_deg_r4 = real(latAdvect,4)* MPC_DEGREES_PER_RADIAN_R4
        ierr = gdxyfll(hco%EZscintID, xpos_r4, ypos_r4, &
             latAdvect_deg_r4, lonAdvect_deg_r4, 1)

        ! determine the bottom-left grid point
        lonIndex = floor(xpos_r4)
        latIndex = floor(ypos_r4)

        ! check if position is east of the grid
        if (floor(xpos_r4) > ni) then
          if (verbose) then
            write(*,*) 'lonIndex0,latIndex0,lon,lat,stepIndexRF,x,y xpos_r4 > ni :', &
                 lonIndex0,latIndex0,lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                 latAdvect*MPC_DEGREES_PER_RADIAN_R8,stepIndexRF,xpos_r4,ypos_r4
          end if
          ! add 10*epsilon(real*4) to ensure do not go too far due to limited precision
          lonAdvect = lonAdvect - 2.0D0*MPC_PI_R8 + 10.0D0*real(epsilon(1.0),8)
          lonAdvect_deg_r4 = real(lonAdvect,4)* MPC_DEGREES_PER_RADIAN_R4
          latAdvect_deg_r4 = real(latAdvect,4)* MPC_DEGREES_PER_RADIAN_R4
          ierr = gdxyfll(hco%EZscintID, xpos_r4, ypos_r4, &
               latAdvect_deg_r4, lonAdvect_deg_r4, 1)
          if (verbose) then
            write(*,*) 'new                            xpos_r4 > ni :', &
                 lonIndex0,latIndex0,lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                 latAdvect*MPC_DEGREES_PER_RADIAN_R8,stepIndexRF,xpos_r4,ypos_r4
          end if
        end if

        ! check if position is west of the grid
        if (floor(xpos_r4) < 1) then
          if (verbose) then
            write(*,*) 'lonIndex0,latIndex0,lon,lat,stepIndexRF,x,y xpos_r4 <  1 :', &
                 lonIndex0,latIndex0,lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                 latAdvect*MPC_DEGREES_PER_RADIAN_R8,stepIndexRF,xpos_r4,ypos_r4
          end if
          ! subtract 10*epsilon(real*4) to ensure do not go too far due to limited precision
          lonAdvect = lonAdvect + 2.0D0*MPC_PI_R8 - 10.0D0*real(epsilon(1.0),8)
          lonAdvect_deg_r4 = real(lonAdvect,4)* MPC_DEGREES_PER_RADIAN_R4
          latAdvect_deg_r4 = real(latAdvect,4)* MPC_DEGREES_PER_RADIAN_R4
          ierr = gdxyfll(hco%EZscintID, xpos_r4, ypos_r4, &
               latAdvect_deg_r4, lonAdvect_deg_r4, 1)
          if (verbose) then
            write(*,*) 'new                            xpos_r4 <  1 :', &
                 lonIndex0,latIndex0,lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                 latAdvect*MPC_DEGREES_PER_RADIAN_R8,stepIndexRF,xpos_r4,ypos_r4
          end if
        end if

        ! longitude is still outside grid - should not happen!
        if (floor(xpos_r4) > ni) then 
          write(*,*) '***still outside lonIndex > ni: stepIndexRF,subStepIndex,lonIndex0,latIndex0,x,y,uu=', &
               stepIndexRF,subStepIndex,lonIndex0,latIndex0,xpos_r4,ypos_r4,uu
          xpos_r4 = real(ni)
          lonAdvect = hco%lon(ni)
        end if
        if (floor(xpos_r4) <  1) then 
          write(*,*) '***still outside lonIndex < 1 : stepIndexRF,subStepIndex,lonIndex0,latIndex0,x,y,uu=', &
               stepIndexRF,subStepIndex,lonIndex0,latIndex0,xpos_r4,ypos_r4,uu
          xpos_r4 = 1.0
          lonAdvect = hco%lon(1)
        end if

        ! if position is poleward of last lat circle, ensure valid lat index
        if (latIndex > nj) then
          if (verbose) then
            write(*,*) 'lonIndex0,latIndex0,lon,lat,stepIndexRF,x,y ypos_r4 > nj :', &
                 lonIndex0,latIndex0,lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                 latAdvect*MPC_DEGREES_PER_RADIAN_R8,stepIndexRF,xpos_r4,ypos_r4
          end if
          ypos_r4 = real(nj)
          latAdvect = hco%lat(nj)
        end if

        ! if position is poleward of first lat circle, ensure valid lat index
        if (latIndex < 1) then
          if (verbose) then
            write(*,*) 'lonIndex0,latIndex0,lon,lat,stepIndexRF,x,y ypos_r4 <  1 :', &
                 lonIndex0,latIndex0,lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                 latAdvect*MPC_DEGREES_PER_RADIAN_R8,stepIndexRF,xpos_r4,ypos_r4
          end if
          ypos_r4 = 1.0
          latAdvect = hco%lat(1)
        end if

        ! determine bottom left grid point again after possible adjustments
        lonIndex = floor(xpos_r4)
        latIndex = floor(ypos_r4)

      end do ! subStepIndex

    end do ! stepIndexRF

    if (verbose) write(*,*) 'final, initial lonIndex,latIndex', lonIndex0,latIndex0,lonIndex,latIndex

    delx = real(xpos_r4,8) - real(lonIndex,8)
    dely = real(ypos_r4,8) - real(latIndex,8)

    interpWeight_BL = min(max( (1.d0-delx) * (1.d0-dely), 0.0d0), 1.0d0)
    interpWeight_BR = min(max(       delx  * (1.d0-dely), 0.0d0), 1.0d0)
    interpWeight_TL = min(max( (1.d0-delx) *       dely , 0.0d0), 1.0d0)
    interpWeight_TR = min(max(       delx  *       dely , 0.0d0), 1.0d0)

    sumWeight = interpWeight_BL + interpWeight_BR + interpWeight_TL +  interpWeight_TR

    if (sumWeight > 1.1d0) then
      write(*,*) 'sumWeight > 1.1 : ', sumWeight
      write(*,*) '          BL, BR, TL, TR=', &
           interpWeight_BL, interpWeight_BR, interpWeight_TL, interpWeight_TR
      write(*,*) '          lonIndex0, latIndex0, lonIndex, latIndex, delx, dely  =', &
           lonIndex0, latIndex0, lonIndex, latIndex, delx, dely
      interpWeight_BL = 0.25d0
      interpWeight_BR = 0.25d0
      interpWeight_TL = 0.25d0
      interpWeight_TR = 0.25d0
    end if

  end SUBROUTINE calcTrajAndWeights
  
  !--------------------------------------------------------------------------
  ! adv_ensemble_tl
  !--------------------------------------------------------------------------
  SUBROUTINE adv_ensemble_tl( ens, adv, nEns )
    implicit none
    type(struct_ens)    :: ens
    type(struct_adv)    :: adv
    integer, intent(in) :: nEns

    if ( adv%nLev_M /= ens_getnLev_M(ens) .or. adv%nLev_T /= ens_getnLev_T(ens) ) then
      call utl_abort('adv_ensemble_tl: vertical levels are not compatible')
    end if

    if      ( ens_getDataKind(ens) == 8 ) then
      call adv_ensemble_tl_r8( ens, adv, nEns )
    else if ( ens_getDataKind(ens) == 4 ) then
      call adv_ensemble_tl_r4( ens, adv, nEns )
    else
      call utl_abort('adv_ensemble_tl: ens%dataKind not valid')
    end if

  END SUBROUTINE adv_ensemble_tl

  !--------------------------------------------------------------------------
  ! adv_ensemble_tl_r8
  !--------------------------------------------------------------------------
  SUBROUTINE adv_ensemble_tl_r8( ens, adv, nEns )
    implicit none
    type(struct_ens)    :: ens
    type(struct_adv)    :: adv
    integer, intent(in) :: nEns

    real(8), pointer     :: ens_oneLev(:,:,:,:)
    real(8), allocatable :: ens1_mpiglobal_tiles(:,:,:,:)
    real(8), allocatable :: ens1_mpiglobal(:,:,:)

    integer :: memberIndex, stepIndex, levIndex, lonIndex, latIndex, kIndex
    integer :: lonIndex2, latIndex2, lonIndex2_p1, latIndex2_p1, nsize, ierr
    integer :: procID, procIDx, procIDy, lonIndex_mpiglobal, latIndex_mpiglobal
    integer :: levTypeIndex, stepIndexAF

    logical :: gatheringDone

    character(len=4) :: varName

    allocate(ens1_mpiglobal_tiles(nEns,adv%lonPerPE,adv%latPerPE,mpi_nprocs))
    allocate(ens1_mpiglobal(nEns,adv%ni,adv%nj))

    do kIndex = 1, ens_getNumK(ens)

      levIndex = ens_getLevFromK    (ens,kIndex)
      varName  = ens_getVarNameFromK(ens,kIndex)
      if      (vnl_varLevelFromVarname(varName) == 'MM') then
        levTypeIndex = MMindex
      else if (vnl_varLevelFromVarname(varName) == 'TH') then
        levTypeIndex = THindex
      else if (vnl_varLevelFromVarname(varName) == 'SF') then
        levTypeIndex = SFindex
      end if

      ens_oneLev => ens_getOneLev_r8(ens,kIndex)

      gatheringDone = .false.

      do stepIndexAF = 1, adv%nTimeStep

        if (stepIndexAF == adv%timeStepIndexMainSource) cycle ! no interpolation needed for this time step

        if (.not. gatheringDone ) then

          ! gather the global field to be interpolated on all tasks
          nsize = nEns*adv%lonPerPE*adv%latPerPE
          call rpn_comm_allgather(ens_oneLev(1:nEns,adv%timeStepIndexSource(stepIndexAF),:,:), nsize, "mpi_double_precision",  &
                                  ens1_mpiglobal_tiles(:,:,:,:), nsize, "mpi_double_precision",  &
                                  "GRID", ierr )

          ! rearrange gathered winds for convenience
!$OMP PARALLEL DO PRIVATE (procIDy,procIDx,procID,latIndex,lonIndex,latIndex_mpiglobal,lonIndex_mpiglobal,memberIndex)
          do procIDy = 0, (mpi_npey-1)
            do procIDx = 0, (mpi_npex-1)
              procID = procIDx + procIDy*mpi_npex
              do latIndex = 1, adv%latPerPE
                latIndex_mpiglobal = latIndex + adv%allLatBeg(procIDy+1) - 1
                do lonIndex = 1, adv%lonPerPE
                  lonIndex_mpiglobal = lonIndex + adv%allLonBeg(procIDx+1) - 1
                  do memberIndex = 1, nEns
                    ens1_mpiglobal(memberIndex,lonIndex_mpiglobal, latIndex_mpiglobal) = &
                         ens1_mpiglobal_tiles(memberIndex, lonIndex, latIndex, procID+1)
                  end do ! memberIndex
                end do ! lonIndex
              end do ! latIndex
            end do ! procIDx
          end do ! procIDy
!$OMP END PARALLEL DO

          if (adv%singleTimeStepIndexSource) gatheringDone = .true.

        end if

!$OMP PARALLEL DO PRIVATE (latIndex,lonIndex,lonIndex2,latIndex2,lonIndex2_p1,latIndex2_p1,memberIndex)
        do latIndex = adv%myLatBeg, adv%myLatEnd
          do lonIndex = adv%myLonBeg, adv%myLonEnd
            lonIndex2 = adv%timeStep(stepIndexAF)%levType(levTypeIndex)%lonIndex(lonIndex,latIndex,levIndex)
            latIndex2 = adv%timeStep(stepIndexAF)%levType(levTypeIndex)%latIndex(lonIndex,latIndex,levIndex)
            lonIndex2_p1 = mod(lonIndex2,adv%ni)+1 ! assume periodic
            latIndex2_p1 = min(latIndex2+1,adv%nj)
            do memberIndex = 1, nEns
              ens_oneLev(memberIndex,stepIndexAF,lonIndex,latIndex) =                                           &
                   adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_BL(lonIndex,latIndex,levIndex)* &
                   ens1_mpiglobal(memberIndex, lonIndex2   ,latIndex2   ) +                                 &
                   adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_BR(lonIndex,latIndex,levIndex)* &
                   ens1_mpiglobal(memberIndex, lonIndex2_p1,latIndex2   ) +                                 &
                   adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_TL(lonIndex,latIndex,levIndex)* &
                   ens1_mpiglobal(memberIndex, lonIndex2   ,latIndex2_p1) +                                 &
                   adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_TR(lonIndex,latIndex,levIndex)* &
                   ens1_mpiglobal(memberIndex, lonIndex2_p1,latIndex2_p1)
            end do ! memberIndex
          end do ! lonIndex
        end do ! latIndex
!$OMP END PARALLEL DO

      end do ! stepIndexAF
      
    end do ! kIndex

    deallocate(ens1_mpiglobal_tiles)
    deallocate(ens1_mpiglobal)

  END SUBROUTINE adv_ensemble_tl_r8

  !--------------------------------------------------------------------------
  ! adv_ensemble_tl_r8
  !--------------------------------------------------------------------------
  SUBROUTINE adv_ensemble_tl_r4( ens, adv, nEns )
    implicit none
    type(struct_ens)    :: ens
    type(struct_adv)    :: adv
    integer, intent(in) :: nEns

    real(4), pointer     :: ens_oneLev(:,:,:,:)
    real(4), allocatable :: ens1_mpiglobal_tiles(:,:,:,:)
    real(4), allocatable :: ens1_mpiglobal(:,:,:)

    integer :: memberIndex, stepIndex, levIndex, lonIndex, latIndex, kIndex
    integer :: lonIndex2, latIndex2, lonIndex2_p1, latIndex2_p1, nsize, ierr
    integer :: procID, procIDx, procIDy, lonIndex_mpiglobal, latIndex_mpiglobal
    integer :: levTypeIndex, stepIndexAF

    logical :: gatheringDone

    character(len=4) :: varName

    allocate(ens1_mpiglobal_tiles(nEns,adv%lonPerPE,adv%latPerPE,mpi_nprocs))
    allocate(ens1_mpiglobal(nEns,adv%ni,adv%nj))

    do kIndex = 1, ens_getNumK(ens)

      levIndex = ens_getLevFromK    (ens,kIndex)
      varName  = ens_getVarNameFromK(ens,kIndex)
      if      (vnl_varLevelFromVarname(varName) == 'MM') then
        levTypeIndex = MMindex
      else if (vnl_varLevelFromVarname(varName) == 'TH') then
        levTypeIndex = THindex
      else if (vnl_varLevelFromVarname(varName) == 'SF') then
        levTypeIndex = SFindex
      end if

      ens_oneLev => ens_getOneLev_r4(ens,kIndex)

      gatheringDone = .false.

      do stepIndexAF = 1, adv%nTimeStep

        if (stepIndexAF == adv%timeStepIndexMainSource) cycle ! no interpolation needed for this time step

        if (.not. gatheringDone ) then

          ! gather the global field to be interpolated on all tasks
          nsize = nEns*adv%lonPerPE*adv%latPerPE
          call rpn_comm_allgather(ens_oneLev(1:nEns,adv%timeStepIndexSource(stepIndexAF),:,:), nsize, "mpi_real4",  &
                                  ens1_mpiglobal_tiles(:,:,:,:), nsize, "mpi_real4",  &
                                  "GRID", ierr )

          ! rearrange gathered winds for convenience
!$OMP PARALLEL DO PRIVATE (procIDy,procIDx,procID,latIndex,lonIndex,latIndex_mpiglobal,lonIndex_mpiglobal,memberIndex)
          do procIDy = 0, (mpi_npey-1)
            do procIDx = 0, (mpi_npex-1)
              procID = procIDx + procIDy*mpi_npex
              do latIndex = 1, adv%latPerPE
                latIndex_mpiglobal = latIndex + adv%allLatBeg(procIDy+1) - 1
                do lonIndex = 1, adv%lonPerPE
                  lonIndex_mpiglobal = lonIndex + adv%allLonBeg(procIDx+1) - 1
                  do memberIndex = 1, nEns
                    ens1_mpiglobal(memberIndex,lonIndex_mpiglobal, latIndex_mpiglobal) = &
                         ens1_mpiglobal_tiles(memberIndex, lonIndex, latIndex, procID+1)
                  end do ! memberIndex
                end do ! lonIndex
              end do ! latIndex
            end do ! procIDx
          end do ! procIDy
!$OMP END PARALLEL DO

          if (adv%singleTimeStepIndexSource) gatheringDone = .true. 

        end if

!$OMP PARALLEL DO PRIVATE (latIndex,lonIndex,lonIndex2,latIndex2,lonIndex2_p1,latIndex2_p1,memberIndex)
        do latIndex = adv%myLatBeg, adv%myLatEnd
          do lonIndex = adv%myLonBeg, adv%myLonEnd
            lonIndex2 = adv%timeStep(stepIndexAF)%levType(levTypeIndex)%lonIndex(lonIndex,latIndex,levIndex)
            latIndex2 = adv%timeStep(stepIndexAF)%levType(levTypeIndex)%latIndex(lonIndex,latIndex,levIndex)
            lonIndex2_p1 = mod(lonIndex2,adv%ni)+1 ! assume periodic
            latIndex2_p1 = min(latIndex2+1,adv%nj)
            do memberIndex = 1, nEns
              ens_oneLev(memberIndex,stepIndexAF,lonIndex,latIndex) =                                           &
                   real(adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_BL(lonIndex,latIndex,levIndex),4)* &
                   ens1_mpiglobal(memberIndex, lonIndex2   ,latIndex2   ) +                                 &
                   real(adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_BR(lonIndex,latIndex,levIndex),4)* &
                   ens1_mpiglobal(memberIndex, lonIndex2_p1,latIndex2   ) +                                 &
                   real(adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_TL(lonIndex,latIndex,levIndex),4)* &
                   ens1_mpiglobal(memberIndex, lonIndex2   ,latIndex2_p1) +                                 &
                   real(adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_TR(lonIndex,latIndex,levIndex),4)* &
                   ens1_mpiglobal(memberIndex, lonIndex2_p1,latIndex2_p1)
            end do ! memberIndex
          end do ! lonIndex
        end do ! latIndex
!$OMP END PARALLEL DO

      end do ! stepIndexAF
      
    end do ! kIndex

    deallocate(ens1_mpiglobal_tiles)
    deallocate(ens1_mpiglobal)

  END SUBROUTINE adv_ensemble_tl_r4

  !--------------------------------------------------------------------------
  ! adv_ensemble_ad
  !--------------------------------------------------------------------------
  SUBROUTINE adv_ensemble_ad( ens, adv, nEns )
    implicit none
    type(struct_ens)    :: ens
    type(struct_adv)    :: adv
    integer, intent(in) :: nEns

    real(8), pointer     :: ens_oneLev(:,:,:,:)
    real(8), allocatable :: ens1_mpiglobal(:,:,:)
    real(8), allocatable :: ens1_mpiglobal_tiles(:,:,:,:)
    real(8), allocatable :: ens1_mpiglobal_tiles2(:,:,:,:)

    integer :: memberIndex, stepIndex, levIndex, lonIndex, latIndex, kIndex
    integer :: lonIndex2, latIndex2, lonIndex2_p1, latIndex2_p1, nsize, ierr
    integer :: procID, procIDx, procIDy, lonIndex_mpiglobal, latIndex_mpiglobal
    integer :: levTypeIndex, stepIndexAF

    character(len=4) :: varName

    if ( .not. adv%singleTimeStepIndexSource ) then
      call utl_abort('adv_ensemble_ad cannot deal with multiple timeStep index source')
    end if
    if ( ens_getDataKind(ens) /= 8 ) then
      call utl_abort('adv_ensemble_ad can only deal with double precision (real8) ensembleStateVector')
    end if
    if ( adv%nLev_M /= ens_getnLev_M(ens) .or. adv%nLev_T /= ens_getnLev_T(ens) ) then
      call utl_abort('adv_ensemble_ad: vertical levels are not compatible')
    end if

    allocate(ens1_mpiglobal(nEns,adv%ni,adv%nj))
    allocate(ens1_mpiglobal_tiles (nEns,adv%lonPerPE,adv%latPerPE,mpi_nprocs))
    allocate(ens1_mpiglobal_tiles2(nEns,adv%lonPerPE,adv%latPerPE,mpi_nprocs))

    do kIndex = 1, ens_getNumK(ens)
            
      levIndex = ens_getLevFromK    (ens,kIndex)
      varName  = ens_getVarNameFromK(ens,kIndex)
      if      (vnl_varLevelFromVarname(varName) == 'MM') then
        levTypeIndex = MMindex
      else if (vnl_varLevelFromVarname(varName) == 'TH') then
        levTypeIndex = THindex
      else if (vnl_varLevelFromVarname(varName) == 'SF') then
        levTypeIndex = SFindex
      end if

      ens1_mpiglobal(:,:,:) = 0.0d0
      ens_oneLev => ens_getOneLev_r8(ens,kIndex)

      do latIndex = adv%myLatBeg, adv%myLatEnd
        do lonIndex = adv%myLonBeg, adv%myLonEnd
          do stepIndexAF = 1, adv%nTimeStep
            if (stepIndexAF == adv%timeStepIndexMainSource) cycle ! no interpolation needed for this time step
            ! this is the bottom-left grid point
            lonIndex2 = adv%timeStep(stepIndexAF)%levType(levTypeIndex)%lonIndex(lonIndex,latIndex,levIndex)
            latIndex2 = adv%timeStep(stepIndexAF)%levType(levTypeIndex)%latIndex(lonIndex,latIndex,levIndex)
            lonIndex2_p1 = mod(lonIndex2,adv%ni)+1 ! assume periodic
            latIndex2_p1 = min(latIndex2+1,adv%nj)
!$OMP PARALLEL DO PRIVATE(memberIndex)
            do memberIndex = 1, nEns
              ens1_mpiglobal(memberIndex, lonIndex2   ,latIndex2) = &
                   ens1_mpiglobal(memberIndex, lonIndex2   ,latIndex2) +  &
                   adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_BL(lonIndex,latIndex,levIndex)* &
                      ens_oneLev(memberIndex,stepIndexAF,lonIndex,latIndex)
              ens1_mpiglobal(memberIndex, lonIndex2_p1,latIndex2) = &
                   ens1_mpiglobal(memberIndex, lonIndex2_p1,latIndex2) +  &
                   adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_BR(lonIndex,latIndex,levIndex)* &
                      ens_oneLev(memberIndex,stepIndexAF,lonIndex,latIndex)
              ens1_mpiglobal(memberIndex, lonIndex2   ,latIndex2_p1) = &
                   ens1_mpiglobal(memberIndex, lonIndex2   ,latIndex2_p1) +  &
                   adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_TL(lonIndex,latIndex,levIndex)* &
                      ens_oneLev(memberIndex,stepIndexAF,lonIndex,latIndex)
              ens1_mpiglobal(memberIndex, lonIndex2_p1,latIndex2_p1) = &
                   ens1_mpiglobal(memberIndex, lonIndex2_p1,latIndex2_p1) +  &
                   adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_TR(lonIndex,latIndex,levIndex)* &
                      ens_oneLev(memberIndex,stepIndexAF,lonIndex,latIndex)
            end do ! memberIndex
!$OMP END PARALLEL DO
          end do ! stepIndexAF
        end do ! lonIndex
      end do ! latIndex

      ! redistribute the global initial time field across mpi tasks by tiles
!$OMP PARALLEL DO PRIVATE(procIDy,procIDx,procID,latIndex,latIndex_mpiglobal,lonIndex,lonIndex_mpiglobal,memberIndex)
      do procIDy = 0, (mpi_npey-1)
        do procIDx = 0, (mpi_npex-1)
          procID = procIDx + procIDy*mpi_npex
          do latIndex = 1, adv%latPerPE
            latIndex_mpiglobal = latIndex + adv%allLatBeg(procIDy+1) - 1
            do lonIndex = 1, adv%lonPerPE
              lonIndex_mpiglobal = lonIndex + adv%allLonBeg(procIDx+1) - 1
              do memberIndex = 1, nEns
                ens1_mpiglobal_tiles(memberIndex, lonIndex, latIndex, procID+1) =  &
                     ens1_mpiglobal(memberIndex, lonIndex_mpiglobal, latIndex_mpiglobal)
              end do ! memberIndex
            end do ! lonIndex
          end do ! latIndex
        end do ! procIDx
      end do ! procIDy
!$OMP END PARALLEL DO

      nsize = nEns*adv%lonPerPE*adv%latPerPE
      if (mpi_nprocs > 1) then
        call rpn_comm_alltoall(ens1_mpiglobal_tiles, nsize,"mpi_double_precision",  &
                               ens1_mpiglobal_tiles2,nsize,"mpi_double_precision","GRID",ierr)
      else
        ens1_mpiglobal_tiles2(:,:,:,1) = ens1_mpiglobal_tiles(:,:,:,1)
      end if

      do procID = 0, (mpi_nprocs-1)
!$OMP PARALLEL DO PRIVATE(latIndex,latIndex2,lonIndex,lonIndex2,memberIndex)
        do latIndex = 1, adv%latPerPE
          latIndex2= latIndex + adv%myLatBeg - 1
          do lonIndex = 1, adv%lonPerPE
            lonIndex2 = lonIndex + adv%myLonBeg - 1
            do memberIndex = 1, nEns
              ens_oneLev(memberIndex, adv%timeStepIndexMainSource, lonIndex2, latIndex2) = &
                   ens_oneLev(memberIndex, adv%timeStepIndexMainSource, lonIndex2, latIndex2) +  &
                   ens1_mpiglobal_tiles2(memberIndex, lonIndex, latIndex, procID+1)
            end do ! memberIndex
          end do ! lonIndex
        end do ! latIndex
!$OMP END PARALLEL DO
      end do ! procID

    end do ! levIndex

    deallocate(ens1_mpiglobal)
    deallocate(ens1_mpiglobal_tiles)
    deallocate(ens1_mpiglobal_tiles2)

  END SUBROUTINE adv_ensemble_ad

  !--------------------------------------------------------------------------
  ! adv_statevector_tl
  !--------------------------------------------------------------------------
  SUBROUTINE adv_statevector_tl( statevector, adv)
    implicit none
    type(struct_gsv)    :: statevector
    type(struct_adv)    :: adv

    real(8), pointer     :: field4D(:,:,:,:)
    real(8), allocatable :: field2D_mpiglobal_tiles(:,:,:)
    real(8), allocatable :: field2D_mpiglobal(:,:)

    integer :: stepIndex, levIndex, lonIndex, latIndex, kIndex
    integer :: lonIndex2, latIndex2, lonIndex2_p1, latIndex2_p1, nsize, ierr
    integer :: procID, procIDx, procIDy, lonIndex_mpiglobal, latIndex_mpiglobal
    integer :: levTypeIndex, stepIndexAF

    logical :: gatheringDone

    character(len=4) :: varName

    if ( gsv_getDataKind(statevector) /= 8 ) then
      call utl_abort('adv_statevector_tl can only deal with double precision (real8) gridStateVector')
    end if
    if ( adv%nLev_M /= statevector%vco%nLev_M .or. adv%nLev_T /= statevector%vco%nLev_T ) then
      call utl_abort('adv_statevector_tl: vertical levels are not compatible')
    end if

    allocate(field2D_mpiglobal_tiles(adv%lonPerPE,adv%latPerPE,mpi_nprocs))
    allocate(field2D_mpiglobal(adv%ni,adv%nj))

    do kIndex = 1, gsv_getNumK(statevector)

      levIndex = gsv_getLevFromK    (statevector,kIndex)
      varName  = gsv_getVarNameFromK(statevector,kIndex)
      if      (vnl_varLevelFromVarname(varName) == 'MM') then
        levTypeIndex = MMindex
      else if (vnl_varLevelFromVarname(varName) == 'TH') then
        levTypeIndex = THindex
      else if (vnl_varLevelFromVarname(varName) == 'SF') then
        levTypeIndex = SFindex
      end if

      field4D => gsv_getField_r8(statevector, varName)

      gatheringDone = .false.

      do stepIndexAF = 1, adv%nTimeStep

        if (stepIndexAF == adv%timeStepIndexMainSource) cycle ! no interpolation needed for this time step

        if (.not. gatheringDone ) then

          ! gather the global field to be interpolated on all tasks
          nsize = adv%lonPerPE*adv%latPerPE
          call rpn_comm_allgather(field4D(:,:,levIndex,adv%timeStepIndexSource(stepIndexAF)), nsize, "mpi_double_precision",  &
                                  field2D_mpiglobal_tiles(:,:,:), nsize, "mpi_double_precision",  &
                                  "GRID", ierr )

          ! rearrange gathered winds for convenience
!$OMP PARALLEL DO PRIVATE (procIDy,procIDx,procID,latIndex,lonIndex,latIndex_mpiglobal,lonIndex_mpiglobal)
          do procIDy = 0, (mpi_npey-1)
            do procIDx = 0, (mpi_npex-1)
              procID = procIDx + procIDy*mpi_npex
              do latIndex = 1, adv%latPerPE
                latIndex_mpiglobal = latIndex + adv%allLatBeg(procIDy+1) - 1
                do lonIndex = 1, adv%lonPerPE
                  lonIndex_mpiglobal = lonIndex + adv%allLonBeg(procIDx+1) - 1
                    field2D_mpiglobal(lonIndex_mpiglobal, latIndex_mpiglobal) = &
                         field2D_mpiglobal_tiles(lonIndex, latIndex, procID+1)
                end do ! lonIndex
              end do ! latIndex
            end do ! procIDx
          end do ! procIDy
!$OMP END PARALLEL DO

          if (adv%singleTimeStepIndexSource) gatheringDone = .true. 

        end if

!$OMP PARALLEL DO PRIVATE (latIndex,lonIndex,lonIndex2,latIndex2,lonIndex2_p1,latIndex2_p1)
        do latIndex = adv%myLatBeg, adv%myLatEnd
          do lonIndex = adv%myLonBeg, adv%myLonEnd
            lonIndex2 = adv%timeStep(stepIndexAF)%levType(levTypeIndex)%lonIndex(lonIndex,latIndex,levIndex)
            latIndex2 = adv%timeStep(stepIndexAF)%levType(levTypeIndex)%latIndex(lonIndex,latIndex,levIndex)
            lonIndex2_p1 = mod(lonIndex2,adv%ni)+1 ! assume periodic
            latIndex2_p1 = min(latIndex2+1,adv%nj)
            field4D(lonIndex,latIndex,levIndex,stepIndexAF) =                                                 &
                 adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_BL(lonIndex,latIndex,levIndex)* &
                 field2D_mpiglobal(lonIndex2   ,latIndex2   ) +                                               &
                 adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_BR(lonIndex,latIndex,levIndex)* &
                 field2D_mpiglobal(lonIndex2_p1,latIndex2   ) +                                               &
                 adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_TL(lonIndex,latIndex,levIndex)* &
                 field2D_mpiglobal(lonIndex2   ,latIndex2_p1) +                                               &
                 adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_TR(lonIndex,latIndex,levIndex)* &
                 field2D_mpiglobal(lonIndex2_p1,latIndex2_p1)
          end do ! lonIndex
        end do ! latIndex
!$OMP END PARALLEL DO

      end do ! stepIndexAF

    end do ! kIndex

    deallocate(field2D_mpiglobal_tiles)
    deallocate(field2D_mpiglobal)

  END SUBROUTINE adv_statevector_tl

  !--------------------------------------------------------------------------
  ! adv_statevector_ad
  !--------------------------------------------------------------------------
  SUBROUTINE adv_statevector_ad( statevector, adv)
    implicit none
    type(struct_gsv)    :: statevector
    type(struct_adv)    :: adv

    real(8), pointer     :: field4D(:,:,:,:)
    real(8), allocatable :: field2D_mpiglobal(:,:)
    real(8), allocatable :: field2D_mpiglobal_tiles (:,:,:)
    real(8), allocatable :: field2D_mpiglobal_tiles2(:,:,:)

    integer :: stepIndex, levIndex, lonIndex, latIndex, kIndex
    integer :: lonIndex2, latIndex2, lonIndex2_p1, latIndex2_p1, nsize, ierr
    integer :: procID, procIDx, procIDy, lonIndex_mpiglobal, latIndex_mpiglobal
    integer :: levTypeIndex, stepIndexAF

    character(len=4) :: varName

    if ( adv%singleTimeStepIndexSource ) then
      call utl_abort('adv_statevector_ad cannot work for singleTimeStepIndexSource')
    end if
    if ( gsv_getDataKind(statevector) /= 8 ) then
      call utl_abort('adv_statevector_ad can only deal with double precision (real8) ensembleStateVector')
    end if
    if ( adv%nLev_M /= statevector%vco%nLev_M .or. adv%nLev_T /= statevector%vco%nLev_T ) then
      call utl_abort('adv_statevector_ad: vertical levels are not compatible')
    end if

    allocate(field2D_mpiglobal(adv%ni,adv%nj))
    allocate(field2D_mpiglobal_tiles (adv%lonPerPE,adv%latPerPE,mpi_nprocs))
    allocate(field2D_mpiglobal_tiles2(adv%lonPerPE,adv%latPerPE,mpi_nprocs))

    do kIndex = 1, gsv_getNumK(statevector)
            
      levIndex = gsv_getLevFromK    (statevector,kIndex)
      varName  = gsv_getVarNameFromK(statevector,kIndex)
      if      (vnl_varLevelFromVarname(varName) == 'MM') then
        levTypeIndex = MMindex
      else if (vnl_varLevelFromVarname(varName) == 'TH') then
        levTypeIndex = THindex
      else if (vnl_varLevelFromVarname(varName) == 'SF') then
        levTypeIndex = SFindex
      end if

      field4D => gsv_getField_r8(statevector, varName)

      do stepIndexAF = 1, adv%nTimeStep

        if (stepIndexAF == adv%timeStepIndexMainSource) cycle ! no interpolation needed for this time step

        field2D_mpiglobal(:,:)                      = 0.0d0

        do latIndex = adv%myLatBeg, adv%myLatEnd
          do lonIndex = adv%myLonBeg, adv%myLonEnd
            ! this is the bottom-left grid point
            lonIndex2 = adv%timeStep(stepIndexAF)%levType(levTypeIndex)%lonIndex(lonIndex,latIndex,levIndex)
            latIndex2 = adv%timeStep(stepIndexAF)%levType(levTypeIndex)%latIndex(lonIndex,latIndex,levIndex)
            lonIndex2_p1 = mod(lonIndex2,adv%ni)+1 ! assume periodic
            latIndex2_p1 = min(latIndex2+1,adv%nj)
            field2D_mpiglobal(lonIndex2   ,latIndex2) = &
            field2D_mpiglobal(lonIndex2   ,latIndex2) +  &
                 adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_BL(lonIndex,latIndex,levIndex)* &
                 field4D(lonIndex,latIndex,levIndex,stepIndexAF)
            field2D_mpiglobal(lonIndex2_p1,latIndex2) = &
            field2D_mpiglobal(lonIndex2_p1,latIndex2) +  &
                 adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_BR(lonIndex,latIndex,levIndex)* &
                 field4D(lonIndex,latIndex,levIndex,stepIndexAF)
            field2D_mpiglobal(lonIndex2   ,latIndex2_p1) = &
            field2D_mpiglobal(lonIndex2   ,latIndex2_p1) +  &
                 adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_TL(lonIndex,latIndex,levIndex)* &
                 field4D(lonIndex,latIndex,levIndex,stepIndexAF)
            field2D_mpiglobal(lonIndex2_p1,latIndex2_p1) = &
            field2D_mpiglobal(lonIndex2_p1,latIndex2_p1) +  &
                 adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_TR(lonIndex,latIndex,levIndex)* &
                 field4D(lonIndex,latIndex,levIndex,stepIndexAF)
          end do ! lonIndex
        end do ! latIndex

        ! redistribute the global initial time field across mpi tasks by tiles
!$OMP PARALLEL DO PRIVATE(procIDy,procIDx,procID,latIndex,latIndex_mpiglobal,lonIndex,lonIndex_mpiglobal)
        do procIDy = 0, (mpi_npey-1)
          do procIDx = 0, (mpi_npex-1)
            procID = procIDx + procIDy*mpi_npex
            do latIndex = 1, adv%latPerPE
              latIndex_mpiglobal = latIndex + adv%allLatBeg(procIDy+1) - 1
              do lonIndex = 1, adv%lonPerPE
                lonIndex_mpiglobal = lonIndex + adv%allLonBeg(procIDx+1) - 1
                field2D_mpiglobal_tiles(lonIndex, latIndex, procID+1) =  &
                     field2D_mpiglobal(lonIndex_mpiglobal, latIndex_mpiglobal)
              end do ! lonIndex
            end do ! latIndex
          end do ! procIDx
        end do ! procIDy
!$OMP END PARALLEL DO


        nsize = adv%lonPerPE*adv%latPerPE
        if (mpi_nprocs > 1) then
          call rpn_comm_alltoall(field2D_mpiglobal_tiles, nsize,"mpi_double_precision",  &
                                 field2D_mpiglobal_tiles2,nsize,"mpi_double_precision","GRID",ierr)
        else
          field2D_mpiglobal_tiles2(:,:,1) = field2D_mpiglobal_tiles(:,:,1)
        end if

        field4D(:, :, levIndex, adv%timeStepIndexSource(stepIndexAF)) = 0.d0

        do procID = 0, (mpi_nprocs-1)
!$OMP PARALLEL DO PRIVATE(latIndex,latIndex2,lonIndex,lonIndex2)
          do latIndex = 1, adv%latPerPE
            latIndex2= latIndex + adv%myLatBeg - 1
            do lonIndex = 1, adv%lonPerPE
              lonIndex2 = lonIndex + adv%myLonBeg - 1
              field4D(lonIndex2, latIndex2, levIndex, adv%timeStepIndexSource(stepIndexAF)) = &
              field4D(lonIndex2, latIndex2, levIndex, adv%timeStepIndexSource(stepIndexAF)) + &
                   field2D_mpiglobal_tiles2(lonIndex, latIndex, procID+1)
            end do ! lonIndex
          end do ! latIndex
!$OMP END PARALLEL DO
        end do ! procID

      end do ! stepIndexAF

    end do ! levIndex

    deallocate(field2D_mpiglobal)
    deallocate(field2D_mpiglobal_tiles)
    deallocate(field2D_mpiglobal_tiles2)

  END SUBROUTINE adv_statevector_ad

END MODULE advection_mod
