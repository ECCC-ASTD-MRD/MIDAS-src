!-------------------------------------- LICENCE BEGIN ------------------------------------
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

program midas_advector
  !
  ! :Purpose: Main program for the propagation of fields based on 
  !           Lagrangian advection
  !
  use version_mod
  use ramDisk_mod
  use utilities_mod
  use mpi_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use advection_mod
  use ensembleStateVector_mod
  implicit none

  type(struct_hco), pointer :: hco => null()
  type(struct_vco), pointer :: vco => null()

  type(struct_gsv) :: statevector

  type(struct_adv) :: adv_forward
  type(struct_adv) :: adv_backward

  integer :: fnom, fclos, ierr, nulnam, stepIndex, mode
  integer :: get_max_rss, newdate, dateStamp, yyyymmdd, hhmmsshh

  integer,allocatable :: dateStampListAdvectedFields(:)

  integer,parameter   :: maxNumLevels=200

  real(8), allocatable :: advectFactor_M(:)

  ! Namelist variables
  character(len=256)  :: steeringFlowFile
  character(len=256)  :: fileToAdvec
  character(len=256)  :: direction
  integer             :: nEns
  integer             :: dateStart
  integer             :: advectedFieldNumStep, steeringFlowNumStep
  real(8)             :: advectedFieldDelThour, steeringFlowDelThour
  real(8)             :: advectFactor(maxNumLevels)

  namelist /namadvector/ fileToAdvec, steeringFlowFile, nEns, &
                         advectedFieldNumStep, advectedFieldDelThour, &
                         steeringFlowNumStep , steeringFlowDelThour, advectFactor, &
                         dateStart, direction

  call ver_printNameAndVersion('advector','Propagation based on advection')

  !
  !- 1.  Settings and module initializations
  !
  write(*,*)
  write(*,*) '> midas-advector: setup - START'

  !- 1.1 mpi
  call mpi_initialize

  !- 1.2 timings
  call tmg_init(mpi_myid, 'TMG_INFO')
  call utl_tmg_start(0,'Main')

  !- 1.3 Read namelist options
  fileToAdvec           = 'missing'
  steeringFlowFile      = 'missing'
  direction             = 'backward-forward'
  advectFactor(:)       = 0.75d0
  nEns                  = 1
  advectedFieldNumStep  = 2
  steeringFlowNumStep   = 2
  advectedFieldDelThour = 1.d0
  steeringFlowDelThour  = 1.d0
  dateStart             = 1912013100

  nulnam = 0
  ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read(nulnam,nml=namadvector,iostat=ierr)
  if (ierr /=0) call utl_abort('midas-advector: Error reading namelist')
  if (mpi_myid == 0) write(*,nml=namadvector)
  ierr = fclos(nulnam)

  !- 1.4 RAM disk usage
  call ram_setup

  !- 1.5 Temporal grid
  yyyymmdd = dateStart/100
  hhmmsshh = (dateStart-yyyymmdd*100) * 1000000
  mode = 3 ! printable to stamp
  ierr = newdate(dateStamp, yyyymmdd, hhmmsshh, mode)
  write(*,*) ''
  write(*,*) 'midas-advector: starting date '
  write(*,*) '                yyyymmdd  = ',yyyymmdd
  write(*,*) '                hhmmsshh  = ',hhmmsshh
  write(*,*) '                datestamp = ',dateStamp

  select case(trim(direction))
  case('forward')
    write(*,*) '   ... going FORWARD'
  case('forward-backward')
    write(*,*) '   ... going FORWARD and BACKWARD'
  case('backward')
    write(*,*) '   ... going BACKWARD'
  case('backward-forward')
    write(*,*) '   ... going BACKWARD and FORWARD'
  case default
    call utl_abort('midas-advector: unknown direction')
  end select

  allocate(dateStampListAdvectedFields(advectedFieldNumStep))
  if (trim(direction) == 'forward'          .or. &
      trim(direction) == 'forward-backward' ) then
    do stepIndex = 1, advectedFieldNumStep
      call incdatr(dateStampListAdvectedFields(stepIndex),dateStamp, &
                   real(stepIndex-1,8)*advectedFieldDelThour)
    end do
  else
    do stepIndex = advectedFieldNumStep, 1, -1
      call incdatr(dateStampListAdvectedFields(stepIndex),dateStamp, &
                   real(stepIndex-advectedFieldNumStep,8)*advectedFieldDelThour)
    end do
  end if

  !- 1.6 Variables of the model states
  call gsv_setup

  !- 1.7 Set the horizontal domain
  call hco_SetupFromFile(hco, fileToAdvec, ' ') ! IN

  !- 1.8 Initialize the vertical coordinate
  call vco_SetupFromFile( vco,        & ! OUT
                          fileToAdvec ) ! IN

  write(*,*)
  write(*,*) '> midas-advector: setup - END'
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- 2.  Setup the advection
  !
  allocate(advectFactor_M(vco%nLev_M))
  advectFactor_M(:) = advectFactor(1:vco%nLev_T)

  if (trim(direction) == 'forward'          .or. &
      trim(direction) == 'backward-forward' .or. &
      trim(direction) == 'forward-backward' ) then
    call adv_setup( adv_forward,                                               & ! OUT
                    'fromFirstTimeIndex', hco, vco,                            & ! IN
                    advectedFieldNumStep, dateStampListAdvectedFields,         & ! IN
                    steeringFlowNumStep, steeringFlowDelThour, advectFactor_M, & ! IN
                    'allLevs', steeringFlowFilename_opt=steeringFlowFile )       ! IN
  end if

  if (trim(direction) == 'backward'         .or. &
      trim(direction) == 'backward-forward' .or. &
      trim(direction) == 'forward-backward' ) then
    call adv_setup( adv_backward,                                              & ! OUT
                    'fromLastTimeIndex', hco, vco,                             & ! IN
                    advectedFieldNumStep, dateStampListAdvectedFields,         & ! IN
                    steeringFlowNumStep, steeringFlowDelThour, advectFactor_M, & ! IN
                    'allLevs', steeringFlowFilename_opt=steeringFlowFile )       ! IN
  end if

  deallocate(advectFactor_M)

  !
  !- 3.  Read the input file and advect it
  !
  call gsv_allocate(statevector,advectedFieldNumStep, hco, vco, &
                    dateStampList_opt=dateStampListAdvectedFields, &
                    mpi_local_opt=.true., &
                    allocHeight_opt=.false., allocPressure_opt=.false.)

  !- 3.1 Forward advection
  if (trim(direction) == 'forward'          .or. &
      trim(direction) == 'forward-backward' ) then
    call gio_readFromFile(statevector,fileToAdvec,' ',' ',stepIndex_opt=1, &
                          unitConversion_opt=.false.)

    call adv_statevector_tl( statevector,  & ! INOUT
                             adv_forward )   ! IN

    do stepIndex = 1, advectedFieldNumStep
      call gio_writeToFile(statevector,'./advectedFields.fst','FORWARD',     & ! IN
                           stepIndex_opt=stepIndex, unitConversion_opt=.false.)! IN
    end do
  end if

  !- Backward advection of the forwardly advected fields
  if (trim(direction) == 'forward-backward' ) then
    call adv_statevector_tl( statevector,  & ! INOUT
                             adv_backward )  ! IN

    do stepIndex = 1, advectedFieldNumStep
      call gio_writeToFile(statevector,'./advectedFields.fst','FORW_BACK',    & ! IN
                           stepIndex_opt=stepIndex, unitConversion_opt=.false.) ! IN
    end do
  end if

  !- 3.2 Backward advection
  if (trim(direction) == 'backward'          .or. &
      trim(direction) == 'backward-forward' ) then
    call gio_readFromFile(statevector,fileToAdvec,' ',' ',stepIndex_opt=advectedFieldNumStep, &
                          unitConversion_opt=.false.)

    call adv_statevector_tl( statevector,  & ! INOUT
                             adv_backward )  ! IN

    do stepIndex = 1, advectedFieldNumStep
      call gio_writeToFile(statevector,'./advectedFields.fst','BACKWARD',     & ! IN
                           stepIndex_opt=stepIndex, unitConversion_opt=.false.) ! IN
    end do
  end if

  !- Forward advection of the backwardly advected fields
  if (trim(direction) == 'backward-forward' ) then
    call adv_statevector_tl( statevector,  & ! INOUT
                             adv_forward )   ! IN

    do stepIndex = 1, advectedFieldNumStep
      call gio_writeToFile(statevector,'./advectedFields.fst','BACK_FORW',    & ! IN
                           stepIndex_opt=stepIndex, unitConversion_opt=.false.) ! IN
    end do
  end if

  call gsv_deallocate(statevector)

  !
  !- 4.  Ending
  !
  write(*,*)
  write(*,*) '> midas-advector: Ending'
  call utl_tmg_stop(0)
  call tmg_terminate(mpi_myid, 'TMG_INFO')

  call rpn_comm_finalize(ierr) 

end program midas_advector
