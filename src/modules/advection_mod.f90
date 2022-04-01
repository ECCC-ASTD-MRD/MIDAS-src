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

MODULE advection_mod
  ! MODULE advection_mod (prefix="adv" category='3. High-level transformations')
  !
  ! :Purpose: To perform forward and/or backward advection (based on 
  !           semi-lagrangian trajectories) for both gridStateVector and
  !           ensemble of gridStateVectors
  !
  use ramDisk_mod
  use mpi_mod
  use mpivar_mod
  use mathPhysConstants_mod
  use earthConstants_mod
  use timeCoord_mod
  use ensembleStateVector_mod
  use gridStateVector_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use utilities_mod
  use varNameList_mod
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

  integer :: numStepSteeringFlow
  integer :: lonPerPE, latPerPE
  integer, allocatable :: allLonBeg(:), allLatBeg(:)
  integer, allocatable :: numSubStep(:)

  real(8), pointer     :: uu_steeringFlow_ptr4d(:,:,:,:)
  real(8), pointer     :: vv_steeringFlow_ptr4d(:,:,:,:)
  real(8), allocatable :: uu_steeringFlow_mpiGlobal(:,:,:)
  real(8), allocatable :: vv_steeringFlow_mpiGlobal(:,:,:)
  real(8), allocatable :: uu_steeringFlow_ThermoLevel(:,:)
  real(8), allocatable :: vv_steeringFlow_ThermoLevel(:,:)

  real(8), allocatable :: steeringFlowFactor(:)

  real(8) :: steeringFlowDelTsec

  type(struct_hco), pointer :: hco
  type(struct_vco), pointer :: vco

  logical :: nlat_equalAcrossMpiTasks, nlon_equalAcrossMpiTasks

  ! Control parameter for the level of listing output
  logical, parameter :: verbose = .false.

CONTAINS

  !--------------------------------------------------------------------------
  ! adv_Setup
  !--------------------------------------------------------------------------
  SUBROUTINE adv_setup(adv, mode, hco_in, vco_in, numStepAdvectedField, &
                       dateStampListAdvectedField, numStepSteeringFlow_in, steeringFlowDelThour, &
                       steeringFlowFactor_in, levTypeList, steeringFlowFilename_opt, &
                       statevector_steeringFlow_opt)
    implicit none

    type(struct_adv) :: adv
    type(struct_hco), pointer :: hco_in
    type(struct_vco), pointer :: vco_in

    character(len=*), intent(in) :: mode
    character(len=*), intent(in) :: levTypeList
    character(len=*), optional, intent(in) :: steeringFlowFilename_opt
    integer, intent(in) :: numStepAdvectedField, numStepSteeringFlow_in
    integer, intent(in) :: dateStampListAdvectedField(numStepAdvectedField)
    real(8), intent(in) :: steeringFlowFactor_in(vco_in%nLev_M)
    real(8), intent(in) :: steeringFlowDelThour

    type(struct_gsv), optional :: statevector_steeringFlow_opt

    integer :: latIndex0, lonIndex0, latIndex, lonIndex, levIndex, stepIndexSF, stepIndexAF, ierr
    integer :: levIndexBelow, levIndexAbove
    integer :: gdxyfll, gdllfxy
    integer :: nLevType

    integer, allocatable :: dateStampListSteeringFlow(:)
    integer, allocatable :: advectedFieldAssociatedStepIndexSF(:)
    integer, allocatable :: advectionSteeringFlowStartingStepIndex(:)
    integer, allocatable :: advectionSteeringFlowEndingStepIndex  (:)

    real(8) :: interpWeight_BL, interpWeight_BR, interpWeight_TL, interpWeight_TR
    real(8), allocatable :: uu_steeringFlow_mpiGlobalTiles(:,:,:,:)
    real(8), allocatable :: vv_steeringFlow_mpiGlobalTiles(:,:,:,:)

    real(4) :: xpos_r4, ypos_r4, xposTH_r4, yposTH_r4
    real(4) :: lonMMbelow_deg_r4, lonMMabove_deg_r4, latMMbelow_deg_r4, latMMabove_deg_r4, lonTH_deg_r4, latTH_deg_r4 
    real(4), allocatable :: xposMM_r4(:,:,:,:), yposMM_r4(:,:,:,:)

    character(len=64) :: filename

    type(struct_gsv) :: statevector_steeringFlow

    logical :: AdvectFileExists

    integer :: nLev, levTypeIndex, stepIndexSF_start, stepIndexSF_end 
    integer :: myLonBeg, myLonEnd
    integer :: myLatBeg, myLatEnd

    !
    !- 1.  Set low-level variables
    !
    numStepSteeringFlow   = numStepSteeringFlow_in
    adv%nTimeStep         = numStepAdvectedField

    allocate(steeringFlowFactor(vco_in%nLev_M))
    do levIndex = 1, vco_in%nLev_M
      steeringFlowFactor(levIndex) = steeringFlowFactor_in(levIndex)
      write(*,*) 'adv_setup: steeringFlowFactor = ', levIndex, steeringFlowFactor(levIndex)
    end do

    allocate(adv%timeStepIndexSource(numStepAdvectedField))

    if (vco_in%Vcode /= 5002 .and. vco_in%Vcode /= 5005 ) then
      call utl_abort('adv_setup: Only vCode 5002 and 5005 are currently supported!')
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
    case ('fromLastTimeIndex')
      adv%timeStepIndexMainSource  = numStepAdvectedField
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
        call utl_abort('adv_setup: numStepAdvectedField cannot be even with direction=towardMiddleTimeIndex') 
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
    steeringFlowDelTsec = steeringFlowDelThour*3600.0D0

    !- 1.2 Grid Size
    hco => hco_in
    adv%ni = hco%ni
    adv%nj = hco%nj
    vco => vco_in
    adv%nLev_M = vco%nLev_M
    adv%nLev_T = vco%nLev_T

    call mpivar_setup_latbands(adv%nj, adv%latPerPE, adv%latPerPEmax, adv%myLatBeg, adv%myLatEnd, & 
         divisible_opt=nlat_equalAcrossMpiTasks)
    call mpivar_setup_lonbands(adv%ni, adv%lonPerPE, adv%lonPerPEmax, adv%myLonBeg, adv%myLonEnd, &
         divisible_opt=nlon_equalAcrossMpiTasks)
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
    !- 2.  Read in the wind data to use for advection
    !
    allocate(dateStampListSteeringFlow(numStepSteeringFlow))

    if (present(steeringFlowFilename_opt) ) then

      if (mpi_myid == 0)  then
        write(*,*)
        write(*,*) 'steeringFlow source taken from input file = ', trim(steeringFlowFilename_opt) 
      end if

      do stepIndexSF = 1, numStepSteeringFlow
        call incdatr(dateStampListSteeringFlow(stepIndexSF), dateStampListAdvectedField(1), &
                     real(stepIndexSF-1,8)*steeringFlowDelThour)
      end do

      call gsv_allocate(statevector_steeringFlow,numStepSteeringFlow, hco, vco, &
                        dateStampList_opt=dateStampListSteeringFlow, &
                        varNames_opt=(/'UU','VV','P0'/), mpi_local_opt=.true., &
                        hInterpolateDegree_opt='LINEAR')
      
      fileName = ram_fullWorkingPath(trim(steeringFlowFilename_opt))
      inquire(file=trim(fileName),exist=AdvectFileExists)
      write(*,*) 'AdvectFileExists', AdvectFileExists
      do stepIndexSF = 1, numStepSteeringFlow
        call gsv_readFromFile( statevector_steeringFlow, fileName, ' ', ' ', stepIndex_opt=stepIndexSF, &
                               containsFullField_opt=.true.)
      end do

      call gsv_getField(statevector_steeringFlow, uu_steeringFlow_ptr4d, 'UU')
      call gsv_getField(statevector_steeringFlow, vv_steeringFlow_ptr4d, 'VV')

    else if (present(statevector_steeringFlow_opt)) then

      if (mpi_myid == 0)  then
        write(*,*)
        write(*,*) 'steeringFlow source = input gridStateVector'
        write(*,*) numStepSteeringFlow
        write(*,*) statevector_steeringFlow_opt%dateStampList(:)
      end if
      dateStampListSteeringFlow(:) = statevector_steeringFlow_opt%dateStampList(:)
      call gsv_getField(statevector_steeringFlow_opt, uu_steeringFlow_ptr4d, 'UU')
      call gsv_getField(statevector_steeringFlow_opt, vv_steeringFlow_ptr4d, 'VV')

    else
      call utl_abort('adv_setup: steeringFlow source was not provided!')
    end if

    !
    !-  3.  Advection setup
    !

    !-  3.1 Match the stepIndex between the reference flow and the fields to be advected
    allocate(advectedFieldAssociatedStepIndexSF(numStepAdvectedField))
    advectedFieldAssociatedStepIndexSF(:) = -1
    do stepIndexAF = 1, numStepAdvectedField
      do stepIndexSF = 1, numStepSteeringFlow
        if ( dateStampListAdvectedField(stepIndexAF) == dateStampListSteeringFlow(stepIndexSF) ) then
          advectedFieldAssociatedStepIndexSF(stepIndexAF) = stepIndexSF
          if (mpi_myid == 0)  then
            write(*,*)
            write(*,*) 'stepIndex Match', stepIndexAF, stepIndexSF
          end if
          exit 
        end if
      end do
      if ( advectedFieldAssociatedStepIndexSF(stepIndexAF) == -1 ) then
        call utl_abort('adv_setup: no match between dateStampListAdvectedField and dateStampListSteeringFlow')
      end if
    end do

    !-  3.2 Set starting, ending and direction parameters
    allocate(advectionSteeringFlowStartingStepIndex(numStepAdvectedField))
    allocate(advectionSteeringFlowEndingStepIndex  (numStepAdvectedField))

    select case(trim(mode))
    case ('fromFirstTimeIndex','fromMiddleTimeIndex','fromLastTimeIndex', &
          'towardFirstTimeIndexInverse','towardMiddleTimeIndexInverse')
      do stepIndexAF = 1, numStepAdvectedField
        advectionSteeringFlowStartingStepIndex(stepIndexAF) = advectedFieldAssociatedStepIndexSF(adv%timeStepIndexMainSource)
        advectionSteeringFlowEndingStepIndex  (stepIndexAF) = advectedFieldAssociatedStepIndexSF(stepIndexAF)
      end do
    case ('towardFirstTimeIndex','towardMiddleTimeIndex')
      do stepIndexAF = 1, numStepAdvectedField
        advectionSteeringFlowStartingStepIndex(stepIndexAF) = advectedFieldAssociatedStepIndexSF(stepIndexAF)
        advectionSteeringFlowEndingStepIndex  (stepIndexAF) = advectedFieldAssociatedStepIndexSF(adv%timeStepIndexMainSource)
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
    allocate(uu_steeringFlow_mpiGlobal(numStepSteeringFlow, adv%ni, adv%nj))
    allocate(vv_steeringFlow_mpiGlobal(numStepSteeringFlow, adv%ni, adv%nj))

    allocate(uu_steeringFlow_mpiGlobalTiles(numStepSteeringFlow, adv%lonPerPE, adv%latPerPE, mpi_nprocs))
    allocate(vv_steeringFlow_mpiGlobalTiles(numStepSteeringFlow, adv%lonPerPE, adv%latPerPE, mpi_nprocs))
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if ( nLevType >= THindex ) then
      allocate(xposMM_r4(numStepSteeringFlow,myLonBeg:myLonEnd,myLatBeg:myLatEnd,adv%nLev_M))
      allocate(yposMM_r4(numStepSteeringFlow,myLonBeg:myLonEnd,myLatBeg:myLatEnd,adv%nLev_M))
    end if

    !- 4.1 Compute the trajectories on the momentum levels
    levTypeIndex = MMindex
    nLev         = adv%nLev_M

    do levIndex = 1, nLev ! loop over levels
      
      if (mpi_myid == 0) write(*,*) 'setupAdvectAmplitude: levIndex = ', levIndex

      call processSteeringFlow(levTypeIndex, levIndex,                                         & ! IN
                               uu_steeringFlow_mpiGlobalTiles, vv_steeringFlow_mpiGlobalTiles, & ! OUT
                               adv%nLev_M, adv%nLev_T, myLatBeg, myLatEnd)                       ! IN

      do stepIndexAF = 1, numStepAdvectedField
        if (stepIndexAF == adv%timeStepIndexMainSource) cycle ! no interpolation needed for this time step

        stepIndexSF_start = advectionSteeringFlowStartingStepIndex(stepIndexAF)
        stepIndexSF_end   = advectionSteeringFlowEndingStepIndex  (stepIndexAF) 

        ! loop over all initial grid points within tile for determining trajectories
        do latIndex0 = adv%myLatBeg, adv%myLatEnd
          do lonIndex0 = adv%myLonBeg, adv%myLonEnd

            call calcTrajectory(xpos_r4, ypos_r4,                                                 & ! OUT
                                latIndex0, lonIndex0, levIndex, stepIndexSF_start, stepIndexSF_end) ! IN

            if ( nLevType >= THindex ) then
              xposMM_r4(stepIndexAF,lonIndex0,latIndex0,levIndex) = xpos_r4
              yposMM_r4(stepIndexAF,lonIndex0,latIndex0,levIndex) = ypos_r4
            end if
            
            call calcWeights(lonIndex, latIndex, interpWeight_BL, interpWeight_BR,   & ! OUT
                             interpWeight_TL, interpWeight_TR,                       & ! OUT
                             xpos_r4, ypos_r4)                                         ! IN

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

    !- 4.2 Thermodynamic levels: Interpolate vertically the positions found on the momentum levels
    if ( nLevType >= THindex ) then
      levTypeIndex = THindex
      nLev         = adv%nLev_T

      do levIndex = 1, nLev ! loop over levels
        do stepIndexAF = 1, numStepAdvectedField
          if (stepIndexAF == adv%timeStepIndexMainSource) cycle ! no interpolation needed for this time step

          do latIndex0 = adv%myLatBeg, adv%myLatEnd
            do lonIndex0 = adv%myLonBeg, adv%myLonEnd

              if (levIndex == 1 .and. vco%Vcode == 5002) then
                ! use top momentum level amplitudes for top thermo level
                xposTH_r4 = xposMM_r4(stepIndexAF,lonIndex0,latIndex0,1)
                yposTH_r4 = yposMM_r4(stepIndexAF,lonIndex0,latIndex0,1)
              else if (levIndex == nLev) then
                ! use surface momentum level amplitudes for surface thermo level
                xposTH_r4 = xposMM_r4(stepIndexAF,lonIndex0,latIndex0,adv%nLev_M)
                yposTH_r4 = yposMM_r4(stepIndexAF,lonIndex0,latIndex0,adv%nLev_M)
              else
                ! for other levels, interpolate momentum positions to get thermo positions (as in GEM)
                if (vco%Vcode == 5002) then
                  levIndexBelow = levIndex
                  levIndexAbove = levIndex-1
                else
                  levIndexBelow = levIndex+1
                  levIndexAbove = levIndex
                end if

                ierr = gdllfxy(hco%EZscintID, latMMbelow_deg_r4, lonMMbelow_deg_r4, &
                               xposMM_r4(stepIndexAF,lonIndex0,latIndex0,levIndexBelow), &
                               yposMM_r4(stepIndexAF,lonIndex0,latIndex0,levIndexBelow), 1) 

                ierr = gdllfxy(hco%EZscintID, latMMabove_deg_r4, lonMMabove_deg_r4, &
                               xposMM_r4(stepIndexAF,lonIndex0,latIndex0,levIndexAbove), &
                               yposMM_r4(stepIndexAF,lonIndex0,latIndex0,levIndexAbove), 1)

                if (lonMMbelow_deg_r4 < 0.0) lonMMbelow_deg_r4 = lonMMbelow_deg_r4 + 360.0
                if (lonMMabove_deg_r4 < 0.0) lonMMabove_deg_r4 = lonMMabove_deg_r4 + 360.0

                if ( abs(lonMMbelow_deg_r4 - lonMMabove_deg_r4) > 180.0 ) then
                  if (lonMMbelow_deg_r4 > 180.0 ) then
                    lonMMbelow_deg_r4 = lonMMbelow_deg_r4 - 360.0
                  else
                    lonMMabove_deg_r4 = lonMMabove_deg_r4 - 360.0
                  end if
                end if
                lonTH_deg_r4 = 0.5 * (lonMMbelow_deg_r4 + lonMMabove_deg_r4)
                if (lonTH_deg_r4 < 0.0) lonTH_deg_r4 = lonTH_deg_r4 + 360.0
                latTH_deg_r4 = 0.5 * (latMMbelow_deg_r4 + latMMabove_deg_r4)

                ierr = gdxyfll(hco%EZscintID, xposTH_r4, yposTH_r4, &
                               latTH_deg_r4, lonTH_deg_r4, 1)
              end if

              ! Compute weights
              call calcWeights(lonIndex, latIndex, interpWeight_BL, interpWeight_BR,   & ! OUT
                               interpWeight_TL, interpWeight_TR,                       & ! OUT
                               xposTH_r4, yposTH_r4)                                     ! IN

              ! Store the final position of the trajectory and interp weights
              adv%timeStep(stepIndexAF)%levType(levTypeIndex)%lonIndex       (lonIndex0,latIndex0,levIndex) = lonIndex
              adv%timeStep(stepIndexAF)%levType(levTypeIndex)%latIndex       (lonIndex0,latIndex0,levIndex) = latIndex
              adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_BL(lonIndex0,latIndex0,levIndex) = interpWeight_BL
              adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_BR(lonIndex0,latIndex0,levIndex) = interpWeight_BR
              adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_TL(lonIndex0,latIndex0,levIndex) = interpWeight_TL
              adv%timeStep(stepIndexAF)%levType(levTypeIndex)%interpWeight_TR(lonIndex0,latIndex0,levIndex) = interpWeight_TR
            end do
          end do
        end do
      end do

      deallocate(xposMM_r4)
      deallocate(yposMM_r4)

    end if

    !- 4.3 Surface level: Use the positions from the lowest momentum levels
    if ( nLevType == SFindex ) then
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
    deallocate(uu_steeringFlow_mpiGlobalTiles)
    deallocate(vv_steeringFlow_mpiGlobalTiles)
    deallocate(uu_steeringFlow_mpiGlobal)
    deallocate(vv_steeringFlow_mpiGlobal)
    if (present(steeringFlowFilename_opt) ) then
      call gsv_deallocate(statevector_steeringFlow)
    end if

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if (mpi_myid == 0) write(*,*) 'adv_setup: done'

  end SUBROUTINE adv_setup

  !--------------------------------------------------------------------------
  ! processSteeringFlow
  !--------------------------------------------------------------------------
  SUBROUTINE processSteeringFlow (levTypeIndex, levIndex, &
                                  uu_steeringFlow_mpiGlobalTiles, vv_steeringFlow_mpiGlobalTiles, &
                                  nLev_M, nLev_T, myLatBeg, myLatEnd)
    implicit none
    integer, intent(in) :: levTypeIndex, levIndex, nLev_M, nLev_T, myLatBeg, myLatEnd
    real(8) :: uu_steeringFlow_mpiGlobalTiles(:,:,:,:)
    real(8) :: vv_steeringFlow_mpiGlobalTiles(:,:,:,:)

    integer :: stepIndexSF, nsize, ierr 
    integer :: procID, procIDx, procIDy, lonIndex, latIndex
    integer :: lonIndex_mpiglobal, latIndex_mpiglobal
    integer :: levIndexBelow, levIndexAbove

    real(8) :: uu, vv, latAdvect

    nsize = lonPerPE*latPerPE

    if ( levTypeIndex == MMindex ) then
      ! No vertical interpolation is needed
      do stepIndexSF = 1, numStepSteeringFlow
        ! gather the winds for this level
        call rpn_comm_allgather(uu_steeringFlow_ptr4d(:,:,levIndex,stepIndexSF)  , nsize, "mpi_double_precision", &
                                uu_steeringFlow_mpiGlobalTiles(stepIndexSF,:,:,:), nsize, "mpi_double_precision", &
                                "GRID", ierr )
        call rpn_comm_allgather(vv_steeringFlow_ptr4d(:,:,levIndex,stepIndexSF)  , nsize, "mpi_double_precision", &
                                vv_steeringFlow_mpiGlobalTiles(stepIndexSF,:,:,:), nsize, "mpi_double_precision", &
                                "GRID", ierr )
      end do
    else if (levTypeIndex == THindex ) then
      ! Vertical interpolation is needed...
      ! The adopted approach follows the vertical interpolation for the amplitude fields in bMatrixEnsemble_mod
      do stepIndexSF = 1, numStepSteeringFlow
        if (levIndex == 1 .and. vco%Vcode == 5002) then
          ! use top momentum level amplitudes for top thermo level
          uu_steeringFlow_ThermoLevel(:,:) = uu_steeringFlow_ptr4d(:,:,1,stepIndexSF)
          vv_steeringFlow_ThermoLevel(:,:) = vv_steeringFlow_ptr4d(:,:,1,stepIndexSF)
        else if (levIndex == nLev_T) then
          ! use surface momentum level amplitudes for surface thermo level
          uu_steeringFlow_ThermoLevel(:,:) = uu_steeringFlow_ptr4d(:,:,nLev_M,stepIndexSF)
          vv_steeringFlow_ThermoLevel(:,:) = vv_steeringFlow_ptr4d(:,:,nLev_M,stepIndexSF)
        else
          ! for other levels, interpolate momentum winds to get thermo winds
          if (vco%Vcode == 5002) then
            levIndexBelow = levIndex
            levIndexAbove = levIndex-1
          else
            levIndexBelow = levIndex+1
            levIndexAbove = levIndex
          end if
          !$OMP PARALLEL DO PRIVATE (latIndex)
          do latIndex = myLatBeg, myLatEnd
            uu_steeringFlow_ThermoLevel(:,latIndex) = 0.5d0*( uu_steeringFlow_ptr4d(:,latIndex,levIndexAbove,stepIndexSF) + &
                 uu_steeringFlow_ptr4d(:,latIndex,levIndexBelow,stepIndexSF) )
            vv_steeringFlow_ThermoLevel(:,latIndex) = 0.5d0*( vv_steeringFlow_ptr4d(:,latIndex,levIndexAbove,stepIndexSF) + &
                 vv_steeringFlow_ptr4d(:,latIndex,levIndexBelow,stepIndexSF) )
          end do
          !$OMP END PARALLEL DO
        end if

        ! gather the INTERPOLATED winds for this level
        call rpn_comm_allgather(uu_steeringFlow_ThermoLevel                      , nsize, "mpi_double_precision",  &
                                uu_steeringFlow_mpiGlobalTiles(stepIndexSF,:,:,:), nsize, "mpi_double_precision",  &
                                "GRID", ierr )
        call rpn_comm_allgather(uu_steeringFlow_ThermoLevel                      , nsize, "mpi_double_precision",  &
                                vv_steeringFlow_mpiGlobalTiles(stepIndexSF,:,:,:), nsize, "mpi_double_precision",  &
                                "GRID", ierr )

      end do

    else
      call utl_abort('processSteeringFlow: invalid levTypeIndex')
    end if

    ! rearrange gathered winds for convenience
    do procIDy = 0, (mpi_npey-1)
      do procIDx = 0, (mpi_npex-1)
        procID = procIDx + procIDy*mpi_npex
        do latIndex = 1, latPerPE
          latIndex_mpiglobal = latIndex + allLatBeg(procIDy+1) - 1
          do lonIndex = 1, lonPerPE
            lonIndex_mpiglobal = lonIndex + allLonBeg(procIDx+1) - 1
            uu_steeringFlow_mpiGlobal(:, lonIndex_mpiglobal, latIndex_mpiglobal) = uu_steeringFlow_mpiGlobalTiles(:, lonIndex, latIndex, procID+1)
            vv_steeringFlow_mpiGlobal(:, lonIndex_mpiglobal, latIndex_mpiglobal) = vv_steeringFlow_mpiGlobalTiles(:, lonIndex, latIndex, procID+1)
          end do
        end do
      end do
    end do

    ! determine the number of time steps required as a function of latitude
    do latIndex = 1, hco%nj
      latAdvect = hco%lat(latIndex)
      if (abs(latAdvect) < latitudePatch*MPC_RADIANS_PER_DEGREE_R8) then
        uu = maxval(abs(uu_steeringFlow_mpiGlobal(:,:,latIndex) /(ec_ra*cos(latAdvect)))) ! in rad/s
        vv = maxval(abs(vv_steeringFlow_mpiGlobal(:,:,latIndex) / ec_ra)) ! in rad/s
      else
        uu = maxval(abs(uu_steeringFlow_mpiGlobal(:,:,latIndex) / ec_ra)) ! in rad/s
        vv = maxval(abs(vv_steeringFlow_mpiGlobal(:,:,latIndex) / ec_ra)) ! in rad/s
      end if
      numSubStep(latIndex) = max( 1,  &
           nint( (steeringFlowDelTsec * steeringFlowFactor(levIndex) * uu) / (numGridPts*(hco%lon(2)-hco%lon(1))) ),  &
           nint( (steeringFlowDelTsec * steeringFlowFactor(levIndex) * vv) / (numGridPts*(hco%lat(2)-hco%lat(1))) ) )
    end do
    if (mpi_myid == 0) write(*,*) 'min and max of numSubStep',minval(numSubStep(:)),maxval(numSubStep(:))

  end SUBROUTINE processSteeringFlow

  !--------------------------------------------------------------------------
  ! calcTrajectory
  !--------------------------------------------------------------------------
  SUBROUTINE calcTrajectory(xpos_r4, ypos_r4, latIndex0, lonIndex0, &
                            levIndex, stepIndexSF_start, stepIndexSF_end)
    implicit none

    real(4), intent(out) :: xpos_r4, ypos_r4
    integer, intent(in)  :: latIndex0, lonIndex0, stepIndexSF_start, stepIndexSF_end
    integer, intent(in)  :: levIndex

    integer :: subStepIndex, stepIndexSF, ierr, gdxyfll, latIndex, lonIndex
    integer :: alfa, ni, nj,  stepIndex_direction, stepIndex_first, stepIndex_last

    real(8) :: uu, vv, subDelT, lonAdvect, latAdvect
    real(8) :: uu_p, vv_p, lonAdvect_p, latAdvect_p, Gcoef, Scoef
    real(4) :: lonAdvect_deg_r4, latAdvect_deg_r4

    ni = hco%ni
    nj = hco%nj

    subDelT = steeringFlowDelTsec/real(numSubStep(latIndex0),8)  ! in seconds

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
    if      (stepIndexSF_end > stepIndexSF_start) then
      ! back trajectory   , stepping backwards
      stepIndex_first     = stepIndexSF_end-1
      stepIndex_last      = stepIndexSF_start
      stepIndex_direction = -1
    else if (stepIndexSF_end < stepIndexSF_start) then
      ! forward trajectory, stepping forward
      stepIndex_first     = stepIndexSF_end
      stepIndex_last      = stepIndexSF_start-1
      stepIndex_direction = 1
    else
      call utl_abort('calcTrajAndWeights: fatal error with stepIndexSF')
    end if

    do stepIndexSF = stepIndex_first, stepIndex_last, stepIndex_direction

      if (verbose) write(*,*) 'stepIndexSF,lonIndex,latIndex=',stepIndexSF,lonIndex,latIndex

      do subStepIndex = 1, numSubStep(latIndex0)

        alfa = (subStepIndex-1)/numSubStep(latIndex0)
        ! perform one timestep
        if (abs(hco%lat(latIndex0)) < latitudePatch*MPC_RADIANS_PER_DEGREE_R8) then
          ! points away from pole, handled normally
          ! determine wind at current location (now at BL point)
          uu = (  alfa *uu_steeringFlow_mpiGlobal(stepIndexSF  ,lonIndex,latIndex) + &
               (1-alfa)*uu_steeringFlow_mpiGlobal(stepIndexSF+1,lonIndex,latIndex) ) &
               /(ec_ra*cos(hco%lat(latIndex))) ! in rad/s
          vv = (   alfa*vv_steeringFlow_mpiGlobal(stepIndexSF  ,lonIndex,latIndex) + &
               (1-alfa)*vv_steeringFlow_mpiGlobal(stepIndexSF+1,lonIndex,latIndex) ) &
               /ec_ra
          ! apply user-specified scale factor to advecting winds
          uu = steeringFlowFactor(levIndex) * uu
          vv = steeringFlowFactor(levIndex) * vv

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
          uu =    alfa *uu_steeringFlow_mpiGlobal(stepIndexSF  ,lonIndex,latIndex) + &
               (1-alfa)*uu_steeringFlow_mpiGlobal(stepIndexSF+1,lonIndex,latIndex)  ! in m/s
          vv =    alfa *vv_steeringFlow_mpiGlobal(stepIndexSF  ,lonIndex,latIndex) + &
               (1-alfa)*vv_steeringFlow_mpiGlobal(stepIndexSF+1,lonIndex,latIndex)
          ! transform wind vector into rotated coordinate system
          Gcoef = ( cos(latAdvect)*cos(hco%lat(latIndex0)) + &
               sin(latAdvect)*sin(hco%lat(latIndex0))*cos(lonAdvect-hco%lon(lonIndex0)) ) / &
               cos(latAdvect_p)
          Scoef = ( sin(hco%lat(latIndex0))*sin(lonAdvect-hco%lon(lonIndex0)) ) / &
               cos(latAdvect_p)
          uu_p = Gcoef * uu - Scoef * vv ! in m/s
          vv_p = Scoef * uu + Gcoef * vv 

          ! apply user-specified scale factor to advecting winds
          uu_p = steeringFlowFactor(levIndex) * uu_p ! in m/s
          vv_p = steeringFlowFactor(levIndex) * vv_p

          ! compute next position (in rotated coord system)
          lonAdvect_p = lonAdvect_p + real(stepIndex_direction,8)*subDelT*uu_p/(ec_ra*cos(latAdvect_p))  ! in radians
          latAdvect_p = latAdvect_p + real(stepIndex_direction,8)*subDelT*vv_p/ec_ra

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
            write(*,*) 'lonIndex0,latIndex0,lon,lat,stepIndexSF,x,y xpos_r4 > ni :', &
                 lonIndex0,latIndex0,lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                 latAdvect*MPC_DEGREES_PER_RADIAN_R8,stepIndexSF,xpos_r4,ypos_r4
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
                 latAdvect*MPC_DEGREES_PER_RADIAN_R8,stepIndexSF,xpos_r4,ypos_r4
          end if
        end if

        ! check if position is west of the grid
        if (floor(xpos_r4) < 1) then
          if (verbose) then
            write(*,*) 'lonIndex0,latIndex0,lon,lat,stepIndexSF,x,y xpos_r4 <  1 :', &
                 lonIndex0,latIndex0,lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                 latAdvect*MPC_DEGREES_PER_RADIAN_R8,stepIndexSF,xpos_r4,ypos_r4
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
                 latAdvect*MPC_DEGREES_PER_RADIAN_R8,stepIndexSF,xpos_r4,ypos_r4
          end if
        end if

        ! longitude is still outside grid - should not happen!
        if (floor(xpos_r4) > ni) then 
          write(*,*) '***still outside lonIndex > ni: stepIndexSF,subStepIndex,lonIndex0,latIndex0,x,y,uu=', &
               stepIndexSF,subStepIndex,lonIndex0,latIndex0,xpos_r4,ypos_r4,uu
          xpos_r4 = real(ni)
          lonAdvect = hco%lon(ni)
        end if
        if (floor(xpos_r4) <  1) then 
          write(*,*) '***still outside lonIndex < 1 : stepIndexSF,subStepIndex,lonIndex0,latIndex0,x,y,uu=', &
               stepIndexSF,subStepIndex,lonIndex0,latIndex0,xpos_r4,ypos_r4,uu
          xpos_r4 = 1.0
          lonAdvect = hco%lon(1)
        end if

        ! if position is poleward of last lat circle, ensure valid lat index
        if (latIndex > nj) then
          if (verbose) then
            write(*,*) 'lonIndex0,latIndex0,lon,lat,stepIndexSF,x,y ypos_r4 > nj :', &
                 lonIndex0,latIndex0,lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                 latAdvect*MPC_DEGREES_PER_RADIAN_R8,stepIndexSF,xpos_r4,ypos_r4
          end if
          ypos_r4 = real(nj)
          latAdvect = hco%lat(nj)
        end if

        ! if position is poleward of first lat circle, ensure valid lat index
        if (latIndex < 1) then
          if (verbose) then
            write(*,*) 'lonIndex0,latIndex0,lon,lat,stepIndexSF,x,y ypos_r4 <  1 :', &
                 lonIndex0,latIndex0,lonAdvect*MPC_DEGREES_PER_RADIAN_R8, &
                 latAdvect*MPC_DEGREES_PER_RADIAN_R8,stepIndexSF,xpos_r4,ypos_r4
          end if
          ypos_r4 = 1.0
          latAdvect = hco%lat(1)
        end if

        ! determine bottom left grid point again after possible adjustments
        lonIndex = floor(xpos_r4)
        latIndex = floor(ypos_r4)

      end do ! subStepIndex

    end do ! stepIndexSF

    if (verbose) write(*,*) 'final, initial xpos,ypos', lonIndex0,latIndex0,xpos_r4, ypos_r4

  end SUBROUTINE calcTrajectory
  
  !--------------------------------------------------------------------------
  ! calcWeights
  !--------------------------------------------------------------------------
  SUBROUTINE calcWeights(lonIndex, latIndex, interpWeight_BL, interpWeight_BR, &
                         interpWeight_TL, interpWeight_TR, xpos_r4, ypos_r4)
    implicit none

    integer, intent(out) :: lonIndex, latIndex
    real(8), intent(out) :: interpWeight_BL, interpWeight_BR, interpWeight_TL, interpWeight_TR
    real(4), intent(in)  :: xpos_r4, ypos_r4

    real(8) :: delx, dely, sumWeight

    ! Determine bottom left grid point
    lonIndex = floor(xpos_r4)
    latIndex = floor(ypos_r4)

    if (lonIndex < 1 .or. lonIndex > hco%ni .or. &
        latIndex < 1 .or. latIndex > hco%nj ) then
      write(*,*)
      write(*,*) 'calcWeights: the input positions are wrong'
      write(*,*) '             xpos_r4, ypos_r4, lonIndex, latIndex = ', xpos_r4, ypos_r4, lonIndex, latIndex
      call utl_abort('calcWeights')
    end if

    ! Compute the surrounding four gridpoint interpolation weights
    delx = real(xpos_r4,8) - real(lonIndex,8)
    dely = real(ypos_r4,8) - real(latIndex,8)

    interpWeight_BL = min(max( (1.d0-delx) * (1.d0-dely), 0.0d0), 1.0d0)
    interpWeight_BR = min(max(       delx  * (1.d0-dely), 0.0d0), 1.0d0)
    interpWeight_TL = min(max( (1.d0-delx) *       dely , 0.0d0), 1.0d0)
    interpWeight_TR = min(max(       delx  *       dely , 0.0d0), 1.0d0)

    ! Verification
    sumWeight = interpWeight_BL + interpWeight_BR + interpWeight_TL +  interpWeight_TR

    if (sumWeight > 1.0d0) then
      write(*,*) 'sumWeight > 1.0 : ', sumWeight
      write(*,*) '          BL, BR, TL, TR=', &
           interpWeight_BL, interpWeight_BR, interpWeight_TL, interpWeight_TR
      write(*,*) '          xpos_r4, ypos_r4, lonIndex, latIndex, delx, dely  =', &
           xpos_r4, ypos_r4, lonIndex, latIndex, delx, dely
      call utl_abort('calcWeights')
    end if

  end SUBROUTINE calcWeights

  !--------------------------------------------------------------------------
  ! adv_ensemble_tl
  !--------------------------------------------------------------------------
  SUBROUTINE adv_ensemble_tl( ens, adv, nEns )
    implicit none
    type(struct_ens)    :: ens
    type(struct_adv)    :: adv
    integer, intent(in) :: nEns

    if ( adv%nLev_M /= ens_getNumLev(ens,'MM') .or. adv%nLev_T /= ens_getNumLev(ens,'TH') ) then
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

    integer :: memberIndex, levIndex, lonIndex, latIndex, kIndex
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

          ! rearrange gathered fields for convenience
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

    integer :: memberIndex, levIndex, lonIndex, latIndex, kIndex
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

          ! rearrange gathered fields for convenience
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

    integer :: memberIndex, levIndex, lonIndex, latIndex, kIndex
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
    if ( adv%nLev_M /= ens_getNumLev(ens,'MM') .or. adv%nLev_T /= ens_getNumLev(ens,'TH') ) then
      call utl_abort('adv_ensemble_ad: vertical levels are not compatible')
    end if
    if ( .not. nlat_equalAcrossMpiTasks .or. .not. nlon_equalAcrossMpiTasks) then
      call utl_abort('adv_ensemble_ad can only deal with even nlon and nlat across all MPI tasks')
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

    integer :: levIndex, lonIndex, latIndex, kIndex
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

    call utl_tmg_start(140,'--ADV_GSV')

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

      call gsv_getField(statevector, field4D, varName)

      gatheringDone = .false.

      do stepIndexAF = 1, adv%nTimeStep

        if (stepIndexAF == adv%timeStepIndexMainSource) cycle ! no interpolation needed for this time step

        if (.not. gatheringDone ) then

          ! gather the global field to be interpolated on all tasks
          call rpn_comm_barrier('GRID',ierr)
          call utl_tmg_start(141,'----ADV_GSV_Comm')
          nsize = adv%lonPerPE*adv%latPerPE
          call rpn_comm_allgather(field4D(:,:,levIndex,adv%timeStepIndexSource(stepIndexAF)), nsize, "mpi_double_precision",  &
                                  field2D_mpiglobal_tiles(:,:,:), nsize, "mpi_double_precision",  &
                                  "GRID", ierr )
          call tmg_stop(141)

          ! rearrange gathered fields for convenience
          call utl_tmg_start(142,'----ADV_GSV_Shuffling')
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
          call tmg_stop(142)

          if (adv%singleTimeStepIndexSource) gatheringDone = .true. 

        end if

        call utl_tmg_start(143,'----ADV_GSV_Calc')

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

        call tmg_stop(143)

      end do ! stepIndexAF

    end do ! kIndex

    deallocate(field2D_mpiglobal_tiles)
    deallocate(field2D_mpiglobal)

    call tmg_stop(140)

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

    integer :: levIndex, lonIndex, latIndex, kIndex
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
    if ( .not. nlat_equalAcrossMpiTasks .or. .not. nlon_equalAcrossMpiTasks) then
      call utl_abort('adv_ensemble_ad can only deal with even nlon and nlat across all MPI tasks')
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

      call gsv_getField(statevector, field4D, varName)

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
