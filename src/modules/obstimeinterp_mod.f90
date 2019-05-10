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

module obsTimeInterp_mod
  ! MODULE obsTimeInterp_mod (prefix='oti' category='7. Low-level data objects and utilities')
  !
  ! :Purpose: To store public variables and procedures related to the time
  !           coordinate.
  !
  use mpi_mod
  use mpivar_mod
  use utilities_mod
  use timecoord_mod
  use obsSpaceData_mod
  implicit none
  save
  private
  
  ! public derived type
  public :: struct_oti

  ! public procedures
  public :: oti_setup, oti_initialized, oti_deallocate
  public :: oti_timeBinning
  public :: oti_setTimeInterpWeight, oti_getTimeInterpWeight, oti_getTimeInterpWeightMpiGlobal
  public :: oti_timeInterpWeightAllZero

  type struct_oti
    real(8), pointer :: timeInterpWeight(:,:) => NULL() ! weights for temporal interpolation to obs times
    real(8), pointer :: timeInterpWeightMpiGlobal(:,:,:) => NULL() ! mpi global version of weights
    logical          :: initialized = .false.
  end type struct_oti

  integer, external :: get_max_rss

contains

  function oti_initialized(oti) result(initialized_out)
    implicit none
    type(struct_oti), pointer :: oti
    logical                   :: initialized_out

    initialized_out = oti%initialized

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
    integer, parameter   :: numFamily = 13
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

    familylist( 1) = 'UA'
    familylist( 2) = 'AI'
    familylist( 3) = 'SF'
    familylist( 4) = 'TO'
    familylist( 5) = 'SW'
    familylist( 6) = 'SC'
    familylist( 7) = 'PR'
    familylist( 8) = 'RO'
    familylist( 9) = 'GP'
    familylist(10) = 'CH'
    familylist(11) = 'TM'
    familylist(12) = 'AL'
    familylist(13) = 'GL'

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
              if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) then
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


  subroutine oti_setup(oti, obsSpaceData, numStep, interpType, flagObsOutside_opt)
    implicit none

    ! arguments
    type(struct_oti), pointer  :: oti
    type(struct_obs)           :: obsSpaceData
    integer                    :: numStep
    character(len=*)           :: interpType
    logical, optional          :: flagObsOutside_opt

    ! locals
    integer           :: headerIndex
    real(8)           :: stepObsIndex

    if ( associated(oti) ) then
      call utl_abort('oti_setup: the supplied oti pointer is not null!')
    endif

    allocate(oti)

    if ( .not.tim_initialized() ) call utl_abort('oti_setup: timeCoord module not initialized')

    if (mpi_myid == 0) write(*,*) ' '
    if (mpi_myid == 0) write(*,*) '-------- Entering oti_setup ---------'
    if (mpi_myid == 0) write(*,*) ' '

    if (mpi_myid == 0) write(*,*) 'oti_setup: Number of step obs for time interpolation : ', numStep

    if (associated(oti%timeInterpWeight)) deallocate(oti%timeInterpWeight)
    allocate(oti%timeInterpWeight(obs_numHeader(obsSpaceData),numStep))
    oti%timeInterpWeight(:,:) = 0.0d0

    do headerIndex = 1, obs_numHeader(obsSpaceData)

      if (numStep == 1) then
        call oti_setTimeInterpWeight(oti, 1.0d0, headerIndex, 1)
      else
        ! building floating point step index
        call tim_getStepObsIndex(stepObsIndex,tim_getDatestamp(),  &
                             obs_headElem_i(obsSpaceData,OBS_DAT,headerIndex),  &
                             obs_headElem_i(obsSpaceData,OBS_ETM,headerIndex), numStep)
        ! leave all weights zero if obs time is out of range, otherwise set weights
        if (floor(stepObsIndex) > numStep) then
          write(*,*) 'oti_setup: stepObsIndex too big=', headerIndex, stepObsIndex
          !call utl_abort('oti_setup: this case should not occur')
        else if (floor(stepObsIndex) < 1) then
          write(*,*) 'oti_setup: stepObsIndex too small=',headerIndex, stepObsIndex
          !call utl_abort('oti_setup: this case should not occur')
        else
          if ( trim(interpType) == 'LINEAR' ) then
            if ( stepObsIndex >= real(numStep,8) ) then
              ! special case not handled by general approach
              call oti_setTimeInterpWeight(oti, 1.0d0, headerIndex, numStep)
            else
              ! general approach for most observations
              call oti_setTimeInterpWeight(oti, 1.0d0-(stepObsIndex-floor(stepObsIndex)), headerIndex, floor(stepObsIndex))
              call oti_setTimeInterpWeight(oti, stepObsIndex-floor(stepObsIndex), headerIndex, floor(stepObsIndex)+1)
            end if
          else if ( trim(interpType) == 'NEAREST' ) then
            call oti_setTimeInterpWeight(oti, 1.0d0, headerIndex, nint(stepObsIndex))
          else
            call utl_abort('oti_setup: unknown interpolation type : ' // trim(interpType))
          end if
        end if
      end if

    end do

    oti%initialized = .true.

    if ( present(flagObsOutside_opt) ) then
      if (flagObsOutside_opt) call oti_flagObsOutsideWindow(oti, obsSpaceData)
    end if

    ! also setup MPI global version of weights, needed for s2c_nl
    call oti_setupMpiGlobal(oti)

    if (mpi_myid == 0) write(*,*) ' '
    if (mpi_myid == 0) write(*,*) '-------- End of oti_setup ---------'
    if (mpi_myid == 0) write(*,*) ' '

  end subroutine oti_setup


  subroutine oti_deallocate(oti)
    implicit none

    ! arguments
    type(struct_oti), pointer :: oti

    if (associated(oti%timeInterpWeight)) deallocate(oti%timeInterpWeight)
    if (associated(oti%timeInterpWeightMpiGlobal)) deallocate(oti%timeInterpWeightMpiGlobal)
    oti%initialized = .false.
    nullify(oti)

  end subroutine oti_deallocate


  subroutine oti_setupMpiGlobal(oti)
    implicit none

    ! arguments
    type(struct_oti), pointer :: oti

    ! locals
    integer              :: numHeader, numHeaderMax, numStep, nsize, ierr
    real(8), allocatable :: timeInterpWeightMax(:,:)

    if ( .not.associated(oti%timeInterpWeight) ) then
      call utl_abort('oti_setupMpiGlobal: oti_setup must first be called')
    end if

    numHeader = size(oti%timeInterpWeight,1)
    numStep = size(oti%timeInterpWeight,2)
    write(*,*) 'oti_setupMpiGlobal: before allreduce ', numHeader ; call flush(6)
    call rpn_comm_allreduce(numHeader, numHeaderMax, 1,  &
                            'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)

    write(*,*) 'oti_setupMpiGlobal: allocating array of dimension ', &
               numHeaderMax, numStep, mpi_nprocs 
    allocate(oti%timeInterpWeightMpiGlobal(numHeaderMax,numStep,mpi_nprocs))

    ! copy over timeInterpWeight into a local array with same size for all mpi tasks
    allocate(timeInterpWeightMax(numHeaderMax,numStep))
    timeInterpWeightMax(:,:) = 0.0d0
    timeInterpWeightMax(1:numHeader,1:numStep) = oti%timeInterpWeight(1:numHeader,1:numStep)

    nsize = numHeaderMax * numStep 
    call rpn_comm_allgather(timeInterpWeightMax,           nsize, 'MPI_REAL8',  &
                            oti%timeInterpWeightMpiGlobal, nsize, 'MPI_REAL8',  &
                            'GRID', ierr)

    deallocate(timeInterpWeightMax)

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine oti_setupMpiGlobal


  subroutine oti_setTimeInterpWeight(oti, weight_in, headerIndex, stepObs)
    implicit none
    type(struct_oti), pointer :: oti
    integer, intent(in)       :: headerIndex, stepObs
    real(8), intent(in)       :: weight_in

    oti%timeInterpWeight(headerIndex, stepObs) = weight_in

  end SUBROUTINE oti_setTimeInterpWeight


  function oti_getTimeInterpWeight(oti, headerIndex, stepObs) result(weight_out)
    implicit none
    type(struct_oti), pointer :: oti
    real(8)                   :: weight_out
    integer, intent(in)       :: headerIndex, stepObs

    weight_out = oti%timeInterpWeight(headerIndex, stepObs)

  end function oti_getTimeInterpWeight


  function oti_getTimeInterpWeightMpiGlobal(oti, headerIndex, stepObs, procIndex) result(weight_out)
    implicit none
    type(struct_oti), pointer :: oti
    real(8)                   :: weight_out
    integer, intent(in)       :: headerIndex, stepObs, procIndex

    weight_out = oti%timeInterpWeightMpiGlobal(headerIndex, stepObs, procIndex)

  end function oti_getTimeInterpWeightMpiGlobal


  function oti_timeInterpWeightAllZero(oti, headerIndex) result(allZero)
    implicit none
    type(struct_oti), pointer :: oti
    logical                   :: allZero
    integer, intent(in)       :: headerIndex

    if ( .not.associated(oti%timeInterpWeight) ) then
      call utl_abort('oti_timeInterpWeightAllZero: oti_setup must first be called')
    end if

    allZero = all(oti%timeInterpWeight(headerIndex, :) == 0.0d0)

  end function oti_timeInterpWeightAllZero


  subroutine oti_flagObsOutsideWindow(oti, obsSpaceData)
    implicit none

    ! arguments
    type(struct_oti), pointer :: oti
    type(struct_obs)          :: obsSpaceData

    ! locals
    integer :: headerIndex, bodyIndex, bodyIndexBeg, bodyIndexEnd

    if ( .not.associated(oti%timeInterpWeight) ) then
      call utl_abort('oti_flagObsOutsideWindow: oti_setup must first be called')
    end if

    do headerIndex = 1, obs_numheader(obsSpaceData)

      if ( oti_timeInterpWeightAllZero(oti, headerIndex) ) then
        ! obs is outside of assimilation window

        write(*,*) 'oti_flagObsOutsideWindow: Observation time outside assimilation window: ',  &
             obs_headElem_i(obsSpaceData,OBS_DAT,headerIndex),obs_headElem_i(obsSpaceData,OBS_ETM,headerIndex)

        ! flag these observations as out of time domain and turn off its assimilation flag
        bodyIndexBeg = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
        bodyIndexEnd = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) + bodyIndexBeg -1
        do bodyIndex = bodyIndexBeg, bodyIndexEnd
          call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
        end do
        call obs_headSet_i(obsSpaceData, OBS_ST1, headerIndex,  &
             ibset( obs_headElem_i(obsSpaceData, OBS_ST1, headerIndex), 05))

        ! set the weight to 1 for time step 1 so that the column will be interpolated
        call oti_setTimeInterpWeight(oti, 1.0d0, headerIndex, 1)
      end if

    end do

  end subroutine oti_flagObsOutsideWindow

end module obsTimeInterp_mod
