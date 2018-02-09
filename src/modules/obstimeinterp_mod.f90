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
!! MODULE obsTimeInterp (prefix="oti")
!!
!! *Purpose*: To store public variables and procedures related to the time coordinate.
!!
!--------------------------------------------------------------------------
module obsTimeInterp_mod
  use mpivar_mod
  use utilities_mod
  use timecoord_mod
  use obsSpaceData_mod
  implicit none
  save
  private
  
  ! public procedures
  public :: oti_setup, oti_setupMpiGlobal, oti_initialized
  public :: oti_timeBinning
  public :: oti_setTimeInterpWeight, oti_getTimeInterpWeight, oti_getTimeInterpWeightMpiGlobal
  public :: oti_timeInterpWeightAllZero

  real(8), pointer :: timeInterpWeight(:,:) => NULL() ! weights for linear temporal interpolation to obs times
  real(8), pointer :: timeInterpWeightMpiGlobal(:,:,:) => NULL() ! mpi global version of weights
  logical          :: initialized = .false.

  integer, external :: get_max_rss

contains

  function oti_initialized() result(initialized_out)
    implicit none
    logical initialized_out

    initialized_out = initialized

  end function oti_initialized

  subroutine oti_timeBinning(obsSpaceData,nstepobs)
    implicit none

    ! arguments
    type(struct_obs) :: obsSpaceData
    integer :: nstepobs

    ! locals
    integer :: stepIndex, headerIndex, familyIndex
    integer :: bodyIndex, bodyIndexBeg, bodyIndexEnd, nsize, ierr
    integer, allocatable :: idataass(:,:), inumheader(:,:)
    integer, allocatable :: my_idataass(:,:), my_inumheader(:,:)
    integer, parameter   :: numFamily = 10
    character(len=2)     :: familylist(numFamily)
    character(len=256)   :: formatspec, formatspec2
    real(8)              :: stepObsIndex

    if ( .not.tim_initialized() ) call utl_abort('oti_timeBinning: timeCoord module not initialized')

    allocate(idataass(numFamily,nStepObs+1))
    allocate(my_idataass(numFamily,nStepObs+1))
    my_idataass(:,:) = 0
    allocate(inumheader(numFamily,nStepObs+1))
    allocate(my_inumheader(numFamily,nStepObs+1))
    my_inumheader(:,:) = 0

    familylist(1) = 'UA'
    familylist(2) = 'AI'
    familylist(3) = 'SF'
    familylist(4) = 'TO'
    familylist(5) = 'SW'
    familylist(6) = 'SC'
    familylist(7) = 'PR'
    familylist(8) = 'RO'
    familylist(9) = 'GP'
    familylist(10) = 'CH'

    do headerIndex = 1, obs_numheader(obsSpaceData)
      call tim_getStepObsIndex(stepObsIndex,tim_getDatestamp(), &
           obs_headElem_i(obsSpaceData,OBS_DAT,headerIndex), &
           obs_headElem_i(obsSpaceData,OBS_ETM,headerIndex),nstepobs)
      if (stepObsIndex > 0.0d0) then
        stepIndex = nint(stepObsIndex)
        bodyIndexBeg = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
        bodyIndexEnd = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + bodyIndexBeg - 1          
        do familyIndex = 1, numFamily
          if (obs_getfamily(obsSpaceData,headerIndex) == familylist(familyIndex)) then
            my_inumheader(familyIndex,stepIndex) = my_inumheader(familyIndex,stepIndex)+1
            my_inumheader(familyIndex,nStepObs+1) = my_inumheader(familyIndex,nStepObs+1)+1
            do bodyIndex = bodyIndexBeg, bodyIndexEnd
              if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == 1) then
                my_idataass(familyIndex,stepIndex) = my_idataass(familyIndex,stepIndex) + 1
                my_idataass(familyIndex,nStepObs+1) = &
                     my_idataass(familyIndex,nStepObs+1) + 1
              end if
            end do
          end if
        end do
      else
        write(*,*) 'oti_timeBinning: observation outside time window:',headerIndex,stepObsIndex
      end if
    end do

    formatspec ='(1X,A6,":"'
    do stepIndex = 1,nStepObs
      formatspec = trim(formatspec)//',1X,I7' ! this is for each time bin
    end do
    formatspec = trim(formatspec)//',1X,I9' ! this is for the total
    formatspec = trim(formatspec)//')'

    formatspec2 = '(1X,A6,":"'
    do stepIndex = 1,nStepObs
      formatspec2 = trim(formatspec2)//',1X,I7'
    end do
    formatspec2 = trim(formatspec2)//',1X,A9)'

    write(*,*) '-----------------------------------------------------------------'
    write(*,*) 'Distribution of number of headers over stepobs ON LOCAL PROCESSOR'
    write(*,trim(formatspec2)) 'Bin#',(stepIndex, stepIndex = 1, nStepObs),'Total'
    do familyIndex = 1, numFamily
      write(*,trim(formatspec)) familylist(familyIndex),(my_inumheader(familyIndex,stepIndex), &
            stepIndex = 1, nStepObs+1)
    end do
    write(*,trim(formatspec)) 'ALL',(sum(my_inumheader(:,stepIndex)), stepIndex = 1, nStepObs+1)
    write(*,*) '----------------------------------------------------------------'
    write(*,*) 'Distribution of assimilated data over stepobs ON LOCAL PROCESSOR'
    write(*,trim(formatspec2)) 'Bin#', (stepIndex, stepIndex = 1, nStepObs),'Total'
    do familyIndex = 1, numFamily
      write(*,trim(formatspec)) familylist(familyIndex), (my_idataass(familyIndex,stepIndex), &
            stepIndex = 1, nStepObs+1)
    end do
    write(*,trim(formatspec)) 'ALL',(sum(my_idataass(:,stepIndex)),stepIndex=1,nStepObs+1)
    write(*,*) '----------------------------------------------------------------'

    nsize = size(inumheader)
    call rpn_comm_allreduce(my_inumheader, inumheader, nsize, &
         "mpi_integer", "mpi_sum", "GRID", ierr)
    deallocate(my_inumheader) 
    nsize = size(idataass)
    call rpn_comm_allreduce(my_idataass, idataass, nsize, &
         "mpi_integer", "mpi_sum", "GRID", ierr)
    deallocate(my_idataass) 
    if (mpi_myid == 0) then
      write(*,*) '----------------------------------------------------------------'
      write(*,*) 'Distribution of number of headers over stepobs ON ALL PROCESSORS'
      write(*,trim(formatspec2)) 'Bin#', (stepIndex, stepIndex = 1, nStepObs), 'Total'
      do familyIndex = 1, numFamily
        write(*,trim(formatspec)) familylist(familyIndex), (inumheader(familyIndex,stepIndex), &
              stepIndex = 1, nStepObs+1)
      end do
      write(*,trim(formatspec)) 'ALL', (sum(inumheader(:,stepIndex)), stepIndex = 1, nStepObs+1)
      write(*,*) '---------------------------------------------------------------'
      write(*,*) 'Distribution of assimilated data over stepobs ON ALL PROCESSORS'
      write(*,trim(formatspec2)) 'Bin#', (stepIndex, stepIndex = 1, nStepObs), 'Total'
      do familyIndex = 1, numFamily
        write(*,trim(formatspec)) familylist(familyIndex), (idataass(familyIndex,stepIndex), &
              stepIndex = 1, nStepObs+1)
      end do
      write(*,trim(formatspec)) 'ALL', (sum(idataass(:,stepIndex)), stepIndex = 1, nStepObs+1)
      write(*,*) '---------------------------------------------------------------'
    end if

    deallocate(idataass)
    deallocate(inumheader)

  end subroutine oti_timeBinning


  subroutine oti_setup(obsSpaceData, numStep)
    implicit none

    ! arguments
    type(struct_obs) :: obsSpaceData
    integer          :: numStep

    ! locals
    integer :: headerIndex
    real(8) :: stepObsIndex

    if ( .not.tim_initialized() ) call utl_abort('oti_setup: timeCoord module not initialized')

    if (mpi_myid == 0) write(*,*) ' '
    if (mpi_myid == 0) write(*,*) '-------- Entering oti_setup ---------'
    if (mpi_myid == 0) write(*,*) ' '

    if (mpi_myid == 0) write(*,*) 'oti_setup: Number of step obs for time interpolation : ', numStep

    if (associated(timeInterpWeight)) deallocate(timeInterpWeight)
    allocate(timeInterpWeight(obs_numHeader(obsSpaceData),numStep))
    timeInterpWeight(:,:) = 0.0d0

    do headerIndex = 1, obs_numHeader(obsSpaceData)

      if (numStep == 1) then
        call oti_setTimeInterpWeight(1.0d0,headerIndex,1)
      else
        ! building floating point step index
        call tim_getStepObsIndex(stepObsIndex,tim_getDatestamp(),  &
                             obs_headElem_i(obsSpaceData,OBS_DAT,headerIndex),  &
                             obs_headElem_i(obsSpaceData,OBS_ETM,headerIndex), numStep)
        ! leave all weights zero if obs time is out of range, otherwise set weights
        if (floor(stepObsIndex) >= numStep) then
          write(*,*) 'oti_setup: stepObsIndex too big=', headerIndex, stepObsIndex
        else if (floor(stepObsIndex) <= 0) then
          write(*,*) 'oti_setup: stepObsIndex too small=',headerIndex, stepObsIndex
        else
          call oti_setTimeInterpWeight(1.0d0-(stepObsIndex-floor(stepObsIndex)), headerIndex, floor(stepObsIndex))
          call oti_setTimeInterpWeight(stepObsIndex-floor(stepObsIndex), headerIndex, floor(stepObsIndex)+1)
        end if
      end if

    end do

    initialized = .true.

    if (mpi_myid == 0) write(*,*) ' '
    if (mpi_myid == 0) write(*,*) '-------- End of oti_setup ---------'
    if (mpi_myid == 0) write(*,*) ' '

  end subroutine oti_setup


  subroutine oti_setupMpiGlobal()
    implicit none

    ! locals
    integer              :: numHeader, numHeaderMax, numStep, nsize, ierr
    real(8), allocatable :: timeInterpWeightMax(:,:)

    if ( .not.associated(timeInterpWeight) ) then
      call utl_abort('oti_setupMpiGlobal: oti_setup must first be called')
    end if

    numHeader = size(timeInterpWeight,1)
    numStep = size(timeInterpWeight,2)
    write(*,*) 'oti_setupMpiGlobal: before allreduce ', numHeader ; call flush(6)
    call rpn_comm_allreduce(numHeader, numHeaderMax, 1,  &
                            'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)

    write(*,*) 'oti_setupMpiGlobal: allocating array of dimension ', &
               numHeaderMax, numStep, mpi_nprocs 
    allocate(timeInterpWeightMpiGlobal(numHeaderMax,numStep,mpi_nprocs))

    ! copy over timeInterpWeight into a local array with same size for all mpi tasks
    allocate(timeInterpWeightMax(numHeaderMax,numStep))
    timeInterpWeightMax(:,:) = 0.0d0
    timeInterpWeightMax(1:numHeader,1:numStep) = timeInterpWeight(1:numHeader,1:numStep)

    nsize = numHeaderMax * numStep 
    call rpn_comm_allgather(timeInterpWeightMax,       nsize, 'MPI_REAL8',  &
                            timeInterpWeightMpiGlobal, nsize, 'MPI_REAL8',  &
                            'GRID', ierr)

    deallocate(timeInterpWeightMax)

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine oti_setupMpiGlobal


  subroutine oti_setTimeInterpWeight(weight_in,headerIndex,stepObs)
    implicit none
    integer, intent(in) :: headerIndex, stepObs
    real(8), intent(in) :: weight_in

    timeInterpWeight(headerIndex, stepObs) = weight_in

  end SUBROUTINE oti_setTimeInterpWeight


  function oti_getTimeInterpWeight(headerIndex,stepObs) result(weight_out)
    implicit none
    real(8)             :: weight_out
    integer, intent(in) :: headerIndex, stepObs

    weight_out = timeInterpWeight(headerIndex, stepObs)

  end function oti_getTimeInterpWeight


  function oti_getTimeInterpWeightMpiGlobal(headerIndex,stepObs,procIndex) result(weight_out)
    implicit none
    real(8)             :: weight_out
    integer, intent(in) :: headerIndex, stepObs, procIndex

    weight_out = timeInterpWeightMpiGlobal(headerIndex, stepObs, procIndex)

  end function oti_getTimeInterpWeightMpiGlobal


  function oti_timeInterpWeightAllZero(headerIndex) result(allZero)
    implicit none
    logical             :: allZero
    integer, intent(in) :: headerIndex

    allZero = all(timeInterpWeight(headerIndex, :) == 0.0d0)

  end function oti_timeInterpWeightAllZero


end module obsTimeInterp_mod
