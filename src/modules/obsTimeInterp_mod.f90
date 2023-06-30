
module obsTimeInterp_mod
  ! MODULE obsTimeInterp_mod (prefix='oti' category='4. Data Object transformations')
  !
  !:Purpose:  To store public variables and procedures related to the time
  !           coordinate.
  !
  use midasMpi_mod
  use utilities_mod
  use timecoord_mod
  use obsSpaceData_mod
  use obsFamilyList_mod
  
  implicit none
  save
  private
  
  ! public derived type
  public :: struct_oti

  ! public procedures
  public :: oti_setup, oti_deallocate
  public :: oti_timeBinning
  public :: oti_setTimeInterpWeight, oti_getTimeInterpWeight, oti_getTimeInterpWeightMpiGlobal
  public :: oti_timeInterpWeightAllZero

  type struct_oti
    real(8), pointer :: timeInterpWeight(:,:) => NULL() ! weights for temporal interpolation to obs times
    real(8), pointer :: timeInterpWeightMpiGlobal(:,:,:) => NULL() ! mpi global version of weights
  end type struct_oti

  integer, parameter :: maxNumWrites = 50

  integer, external :: get_max_rss

contains

  subroutine oti_timeBinning( obsSpaceData, nstepobs )
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(in) :: obsSpaceData
    integer         , intent(in) :: nstepobs

    ! Locals:
    integer :: stepIndex, headerIndex, familyIndex
    integer :: bodyIndex, bodyIndexBeg, bodyIndexEnd, nsize, ierr
    integer, allocatable :: idataass(:,:), inumheader(:,:)
    integer, allocatable :: my_idataass(:,:), my_inumheader(:,:)
    character(len=256)   :: formatspec, formatspec2
    real(8)              :: stepObsIndex
    integer, save :: numWrites = 0

    if ( .not.tim_initialized() ) call utl_abort('oti_timeBinning: timeCoord module not initialized')

    allocate(idataass(ofl_numFamily,nStepObs+1))
    allocate(my_idataass(ofl_numFamily,nStepObs+1))
    my_idataass(:,:) = 0
    allocate(inumheader(ofl_numFamily,nStepObs+1))
    allocate(my_inumheader(ofl_numFamily,nStepObs+1))
    my_inumheader(:,:) = 0

    do headerIndex = 1, obs_numheader(obsSpaceData)
      call tim_getStepObsIndex(stepObsIndex,tim_getDatestamp(), &
           obs_headElem_i(obsSpaceData,OBS_DAT,headerIndex), &
           obs_headElem_i(obsSpaceData,OBS_ETM,headerIndex),nstepobs)
      if (stepObsIndex > 0.0d0) then
        stepIndex = nint(stepObsIndex)
        bodyIndexBeg = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
        bodyIndexEnd = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + bodyIndexBeg - 1          
        do familyIndex = 1, ofl_numFamily
          if (obs_getfamily(obsSpaceData,headerIndex_opt=headerIndex) == ofl_familyList(familyIndex)) then
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
        numWrites = numWrites + 1
        if (numWrites < maxNumWrites) then
          write(*,*) 'oti_timeBinning: observation outside time window:',headerIndex,stepObsIndex
        else if (numWrites == maxNumWrites) then
          write(*,*) 'oti_timeBinning: more observations outside time window, but reached'
          write(*,*) '                 maximum number of writes to the listing.'
        end if
      end if
    end do

    formatspec ='(1x,a6,":"'
    do stepIndex = 1,nStepObs
      formatspec = trim(formatspec)//',1x,i9' ! this is for each time bin
    end do
    formatspec = trim(formatspec)//',1x,i9' ! this is for the total
    formatspec = trim(formatspec)//')'

    formatspec2 = '(1x,a6,":"'
    do stepIndex = 1,nStepObs
      formatspec2 = trim(formatspec2)//',1x,i9'
    end do
    formatspec2 = trim(formatspec2)//',1x,a9)'

    write(*,*) '-----------------------------------------------------------------'
    write(*,*) 'Distribution of number of headers over stepobs ON LOCAL PROCESSOR'
    write(*,trim(formatspec2)) 'Bin#',(stepIndex, stepIndex = 1, nStepObs),'Total'
    do familyIndex = 1, ofl_numFamily
      write(*,trim(formatspec)) ofl_familyList(familyIndex),(my_inumheader(familyIndex,stepIndex), &
            stepIndex = 1, nStepObs+1)
    end do
    write(*,trim(formatspec)) 'ALL',(sum(my_inumheader(:,stepIndex)), stepIndex = 1, nStepObs+1)
    write(*,*) '----------------------------------------------------------------'
    write(*,*) 'Distribution of assimilated data over stepobs ON LOCAL PROCESSOR'
    write(*,trim(formatspec2)) 'Bin#', (stepIndex, stepIndex = 1, nStepObs),'Total'
    do familyIndex = 1, ofl_numFamily
      write(*,trim(formatspec)) ofl_familyList(familyIndex), (my_idataass(familyIndex,stepIndex), &
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
    if (mmpi_myid == 0) then
      write(*,*) '----------------------------------------------------------------'
      write(*,*) 'Distribution of number of headers over stepobs ON ALL PROCESSORS'
      write(*,trim(formatspec2)) 'Bin#', (stepIndex, stepIndex = 1, nStepObs), 'Total'
      do familyIndex = 1, ofl_numFamily
        write(*,trim(formatspec)) ofl_familyList(familyIndex), (inumheader(familyIndex,stepIndex), &
              stepIndex = 1, nStepObs+1)
      end do
      write(*,trim(formatspec)) 'ALL', (sum(inumheader(:,stepIndex)), stepIndex = 1, nStepObs+1)
      write(*,*) '---------------------------------------------------------------'
      write(*,*) 'Distribution of assimilated data over stepobs ON ALL PROCESSORS'
      write(*,trim(formatspec2)) 'Bin#', (stepIndex, stepIndex = 1, nStepObs), 'Total'
      do familyIndex = 1, ofl_numFamily
        write(*,trim(formatspec)) ofl_familyList(familyIndex), (idataass(familyIndex,stepIndex), &
              stepIndex = 1, nStepObs+1)
      end do
      write(*,trim(formatspec)) 'ALL', (sum(idataass(:,stepIndex)), stepIndex = 1, nStepObs+1)
      write(*,*) '---------------------------------------------------------------'
    end if

    deallocate(idataass)
    deallocate(inumheader)

  end subroutine oti_timeBinning


  subroutine oti_setup(oti, obsSpaceData, numStep, headerIndexBeg, headerIndexEnd, &
                       interpType_opt, flagObsOutside_opt)
    !
    implicit none

    ! Arguments:
    type(struct_oti), pointer , intent(out)   :: oti
    type(struct_obs)          , intent(inout) :: obsSpaceData
    integer                   , intent(in)    :: numStep
    integer                   , intent(in)    :: headerIndexBeg
    integer                   , intent(in)    :: headerIndexEnd
    character(len=*), optional, intent(in)    :: interpType_opt
    logical,          optional, intent(in)    :: flagObsOutside_opt

    ! Locals:
    integer             :: headerIndex
    real(8)             :: stepObsIndex
    integer, save       :: numWrites = 0
    
    if ( associated(oti) ) then
      call utl_abort('oti_setup: the supplied oti pointer is not null!')
    endif

    allocate(oti)

    if ( .not.tim_initialized() ) call utl_abort('oti_setup: timeCoord module not initialized')

    if ( numStep > 1 .and. .not. present(interpType_opt) ) then
      call utl_abort('oti_setup: interpType_opt must be specified when numStep > 1')
    end if

    if ( trim(interpType_opt) == 'LINEAR' .and. tim_fullyUseExtremeTimeBins) then
      call utl_abort('oti_setup: LINEAR time interpolation is not compatible with ' // &
                     'tim_fullyUseExtremeTimeBins==.true.')
    end if

    if (mmpi_myid == 0) write(*,*) ' '
    if (mmpi_myid == 0) write(*,*) '-------- Entering oti_setup ---------'
    if (mmpi_myid == 0) write(*,*) ' '

    if (mmpi_myid == 0) write(*,*) 'oti_setup: Number of step obs for time interpolation : ', numStep

    allocate(oti%timeInterpWeight(headerIndexBeg:headerIndexEnd,numStep))
    oti%timeInterpWeight(:,:) = 0.0d0

    do headerIndex = headerIndexBeg, headerIndexEnd

      ! building floating point step index
      call tim_getStepObsIndex(stepObsIndex,tim_getDatestamp(),  &
                               obs_headElem_i(obsSpaceData,OBS_DAT,headerIndex),  &
                               obs_headElem_i(obsSpaceData,OBS_ETM,headerIndex), numStep)
      
      ! leave all weights zero if obs time is out of range, otherwise set weights
      if (.not.tim_fullyUseExtremeTimeBins .and. (ceiling(stepObsIndex) > numStep .or. floor(stepObsIndex) < 1)) then
        numWrites = numWrites + 1
        if (numWrites < maxNumWrites) then
          write(*,'(a,i10,f8.2,i7,2a,i10,i6)') 'oti_setup: observation outside time window, headerIndex =', &
                                               headerIndex, stepObsIndex, &
                                               obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex), ' ', &
                                               obs_elem_c    (obsSpaceData, 'STID' , headerIndex), &
                                               obs_headElem_i(obsSpaceData,OBS_DAT,headerIndex), &
                                               obs_headElem_i(obsSpaceData,OBS_ETM,headerIndex)  
        else if (numWrites == maxNumWrites) then
          write(*,*) 'oti_setup: More obs outside time window, but reached maximum number of writes to the listing.'
        end if
      else if (tim_fullyUseExtremeTimeBins .and. (nint(stepObsIndex) > numStep .or. nint(stepObsIndex) < 1)) then
        numWrites = numWrites + 1
        if (numWrites < maxNumWrites) then
          write(*,'(a,2i,2a,i10,i6)') 'oti_setup: observation outside time window, headerIndex =',headerIndex, &
                                      obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex), ' ', &
                                      obs_elem_c    (obsSpaceData, 'STID' , headerIndex), &
                                      obs_headElem_i(obsSpaceData,OBS_DAT,headerIndex), &
                                      obs_headElem_i(obsSpaceData,OBS_ETM,headerIndex)  
        else if (numWrites == maxNumWrites) then
          write(*,*) 'oti_setup: More obs outside time window, but reached maximum number of writes to the listing.'
        end if
      else
        if (numStep == 1) then
          call oti_setTimeInterpWeight(oti, 1.0d0, headerIndex, 1)
        else
          if (trim(interpType_opt) == 'LINEAR') then
            if (stepObsIndex >= real(numStep,8)) then
              ! special case not handled by general approach
              call oti_setTimeInterpWeight(oti, 1.0d0, headerIndex, numStep)
            else
              ! general approach for most observations
              call oti_setTimeInterpWeight(oti, 1.0d0-(stepObsIndex-floor(stepObsIndex)), headerIndex, floor(stepObsIndex))
              call oti_setTimeInterpWeight(oti, stepObsIndex-floor(stepObsIndex), headerIndex, floor(stepObsIndex)+1)
            end if
          else if ( trim(interpType_opt) == 'NEAREST' ) then
            if (nint(stepObsIndex) > numStep) then
              write(*,*) 'stepObsIndex = ', stepObsIndex
              call utl_abort('oti_setup: stepObsIndex is too large!')
            end if
            call oti_setTimeInterpWeight(oti, 1.0d0, headerIndex, nint(stepObsIndex))
          else
            call utl_abort('oti_setup: unknown interpolation type : ' // trim(interpType_opt))
          end if
        end if
      end if

    end do

    if ( present(flagObsOutside_opt) ) then
      if (flagObsOutside_opt) call oti_flagObsOutsideWindow(oti, obsSpaceData, headerIndexBeg, headerIndexEnd)
    end if

    ! also setup MPI global version of weights, needed for s2c_nl
    call oti_setupMpiGlobal(oti)

    if (mmpi_myid == 0) write(*,*) ' '
    if (mmpi_myid == 0) write(*,*) '-------- End of oti_setup ---------'
    if (mmpi_myid == 0) write(*,*) ' '

  end subroutine oti_setup


  subroutine oti_deallocate(oti)
    implicit none

    ! Arguments:
    type(struct_oti), pointer, intent(inout) :: oti

    if (associated(oti%timeInterpWeight)) deallocate(oti%timeInterpWeight)
    if (associated(oti%timeInterpWeightMpiGlobal)) deallocate(oti%timeInterpWeightMpiGlobal)
    if (associated(oti)) deallocate(oti)
    nullify(oti)

  end subroutine oti_deallocate


  subroutine oti_setupMpiGlobal( oti )
    !
    implicit none

    ! Arguments:
    type(struct_oti), pointer, intent(inout) :: oti

    ! Locals:
    integer              :: numHeader, numHeaderMax, numStep, nsize, ierr
    real(8), allocatable :: timeInterpWeightMax(:,:)

    if ( .not.associated(oti%timeInterpWeight) ) then
      call utl_abort('oti_setupMpiGlobal: oti_setup must first be called')
    end if

    numHeader = size(oti%timeInterpWeight,1)
    numStep = size(oti%timeInterpWeight,2)
    write(*,*) 'oti_setupMpiGlobal: before allreduce ', numHeader, numStep
    call rpn_comm_allreduce(numHeader, numHeaderMax, 1,  &
                            'MPI_INTEGER', 'MPI_MAX', 'GRID', ierr)

    write(*,*) 'oti_setupMpiGlobal: allocating array of dimension ', &
               numHeaderMax, numStep, mmpi_nprocs 
    allocate(oti%timeInterpWeightMpiGlobal(numHeaderMax,numStep,mmpi_nprocs))

    ! copy over timeInterpWeight into a local array with same size for all mpi tasks
    allocate(timeInterpWeightMax(numHeaderMax,numStep))
    timeInterpWeightMax(:,:) = 0.0d0
    timeInterpWeightMax(1:numHeader,1:numStep) = oti%timeInterpWeight(:,1:numStep)

    nsize = numHeaderMax * numStep 
    call rpn_comm_allgather(timeInterpWeightMax,           nsize, 'MPI_REAL8',  &
                            oti%timeInterpWeightMpiGlobal, nsize, 'MPI_REAL8',  &
                            'GRID', ierr)

    deallocate(timeInterpWeightMax)

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine oti_setupMpiGlobal


  subroutine oti_setTimeInterpWeight( oti, weight_in, headerIndex, stepObs )
    !
    implicit none

    ! Arguments:
    type(struct_oti), pointer, intent(inout) :: oti
    integer,                   intent(in)    :: headerIndex
    integer,                   intent(in)    :: stepObs
    real(8),                   intent(in)    :: weight_in

    oti%timeInterpWeight(headerIndex, stepObs) = weight_in

  end subroutine oti_setTimeInterpWeight


  function oti_getTimeInterpWeight( oti, headerIndex, stepObs ) result(weight_out)
    !
    implicit none

    ! Arguments:
    type(struct_oti), pointer, intent(inout) :: oti
    integer,                   intent(in)    :: headerIndex
    integer,                   intent(in)    :: stepObs
    ! Result:
    real(8)                   :: weight_out

    weight_out = oti%timeInterpWeight(headerIndex, stepObs)

  end function oti_getTimeInterpWeight


  function oti_getTimeInterpWeightMpiGlobal( oti, headerIndex, stepObs, procIndex ) result(weight_out)
    !
    implicit none
  
    ! Arguments:
    type(struct_oti), pointer, intent(inout) :: oti
    integer,                   intent(in)    :: headerIndex
    integer,                   intent(in)    :: stepObs
    integer,                   intent(in)    :: procIndex
    ! Result:
    real(8)                   :: weight_out

    weight_out = oti%timeInterpWeightMpiGlobal(headerIndex, stepObs, procIndex)

  end function oti_getTimeInterpWeightMpiGlobal


  function oti_timeInterpWeightAllZero( oti, headerIndex ) result(allZero)
    !
    implicit none

    ! Arguments:
    type(struct_oti), pointer, intent(inout) :: oti
    integer,                   intent(in)    :: headerIndex
    ! Result:
    logical                   :: allZero

    if ( .not.associated(oti%timeInterpWeight) ) then
      call utl_abort('oti_timeInterpWeightAllZero: oti_setup must first be called')
    end if

    allZero = all(oti%timeInterpWeight(headerIndex, :) == 0.0d0)

  end function oti_timeInterpWeightAllZero


  subroutine oti_flagObsOutsideWindow( oti, obsSpaceData, headerIndexBeg, headerIndexEnd )
    !
    implicit none

    ! Arguments:
    type(struct_oti), pointer, intent(inout) :: oti
    type(struct_obs)         , intent(inout) :: obsSpaceData
    integer                  , intent(in)    :: headerIndexBeg
    integer                  , intent(in)    :: headerIndexEnd

    ! Locals:
    integer :: headerIndex, bodyIndex, bodyIndexBeg, bodyIndexEnd
    integer :: obsDAT, obsETM
    integer, save :: numWrites = 0

    if ( .not.associated(oti%timeInterpWeight) ) then
      call utl_abort('oti_flagObsOutsideWindow: oti_setup must first be called')
    end if

    do headerIndex = headerIndexBeg, headerIndexEnd

      if ( oti_timeInterpWeightAllZero(oti, headerIndex) ) then
        ! obs is outside of assimilation window

        numWrites = numWrites + 1
        if (numWrites < maxNumWrites) then
          obsDAT = obs_headElem_i(obsSpaceData,OBS_DAT,headerIndex)
          obsETM = obs_headElem_i(obsSpaceData,OBS_ETM,headerIndex)
          write(*,*) 'oti_flagObsOutsideWindow: Observation time outside assimilation window: ',  &
               obsDAT, obsETM
        else if (numWrites == maxNumWrites) then
          write(*,*) 'oti_flagObsOutsideWindow: More rejects, but reached maximum number of writes to the listing.'
        end if

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
