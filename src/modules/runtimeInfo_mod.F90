
module runtimeInfo_mod
  !
  ! MODULE runtimeInfo_mod (prefix='rti' category='8. Low-level utilities and constants')
  !
  !:Purpose: Store and print the MIDAS version number for display in the listing.
  !
  use mpi_f08
  use omp_lib
#ifdef __INTEL_LLVM_COMPILER
  use ifcore, only: tracebackqq
#endif
  use rmn_fnom

  implicit none

  save
  private

  ! Public routines
  public :: rti_abort, rti_printTime
  public :: rti_tmg_start, rti_tmg_stop
  public :: rti_writeStatus

contains

  subroutine rti_abort(message)
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: message

    ! Locals:
    integer :: ierr

    write(6,9000) message
9000 format(//,4X,"!!!---ABORT---!!!",/,8X,"MIDAS stopped in ",A)
    flush(6)

#ifdef __INTEL_LLVM_COMPILER
    ! We must provide a 'user_exit_code'.  With only 'status', the
    ! code aborts rightaway with the chance to run 'mpi_abort' after.
    call tracebackqq(user_exit_code=-1, status=ierr)
#endif

    call mpi_abort(mpi_comm_world, 1, ierr)

  end subroutine rti_abort

  !--------------------------------------------------------------------------
  ! rti_printTime
  !--------------------------------------------------------------------------
  subroutine rti_printTime(reset_opt)
    !
    !:Purpose: Print the elapsed time in the listing. Use of the optional
    !          argument `reset_opt=.true.` resets the accumulator to zero.
    !
    implicit none

    ! Arguments:
    logical, optional, intent(in) :: reset_opt ! Allow user to reset the accumulator

    ! Locals:
    real(8), save :: startTime = -1.0d0
    real(8), save :: accumulatedStart = -1.0d0
    real(8), save :: previousTime = -1.0d0
    real(8)       :: currentTime
    logical, save :: firstCall = .true.
    logical       :: reset
    character(len=8)  :: dateString
    character(len=10) :: timeString

    if (present(reset_opt)) then
      reset = reset_opt
    else
      reset = .false.
    end if

    currentTime = omp_get_wtime()

    if (startTime < 0.0d0) then
      startTime = currentTime
    end if

    if (previousTime < 0.0d0) then
      previousTime = currentTime
    end if

    if (accumulatedStart < 0.0d0 .or. reset) then
      accumulatedStart = currentTime
    end if

    ! Also get the actual date and time
    call date_and_time(dateString, timeString)

    if (firstCall) then
      write(*,'(A,A)') &
           ' rti_printTime: '//dateString//' '//timeString(1:2)//'h '//timeString(3:4)//'m '//timeString(5:10)//'s, ', &
           'First call, counters initialized'
    end if

    if (reset .and. .not.firstCall) then
      write(*,'(A,A)') &
           ' rti_printTime: '//dateString//' '//timeString(1:2)//'h '//timeString(3:4)//'m '//timeString(5:10)//'s, ', &
           'Accumulator reset'
    end if

    if (.not. firstCall) then
      write(*,'(A,A,f10.4,A,A,f10.4,A,A,f10.4,A)') &
           ' rti_printTime: '//dateString//' '//timeString(1:2)//'h '//timeString(3:4)//'m '//timeString(5:10)//'s, ', &
                          'deltaT = ', (currentTime - previousTime), ' s, ', &
                          'accumT = ', (currentTime - accumulatedStart), ' s, ', &
                          'totalT = ', (currentTime - startTime), ' s'
    end if

    previousTime = currentTime
    firstCall = .false.

  end subroutine rti_printTime

  !--------------------------------------------------------------------------
  ! rti_tmg_start
  !--------------------------------------------------------------------------
  subroutine rti_tmg_start(blockIndex, blockLabel)
    !
    !:Purpose: Wrapper for 'hpcoperf' subroutine 'tmg_start'
    !
    implicit none

    ! Arguments:
    integer,          intent(in) :: blockIndex
    character(len=*), intent(in) :: blockLabel

    ! Locals:
    integer            :: labelLength
    integer, parameter :: labelPaddedLength = 40
    character(len=labelPaddedLength) :: blockLabelPadded

    ! only the first thread does the timing
    if (omp_get_thread_num() > 0) return

    blockLabelPadded = '........................................'
    labelLength = min(len_trim(blockLabel), labelPaddedLength)
    blockLabelPadded(1:labelLength) = blockLabel(1:labelLength)

    call tmg_start(blockIndex, blockLabelPadded)

  end subroutine rti_tmg_start

  !--------------------------------------------------------------------------
  ! rti_tmg_stop
  !--------------------------------------------------------------------------
  subroutine rti_tmg_stop(blockIndex)
    !
    !:Purpose: Wrapper for 'hpcoperf' subroutine 'tmg_stop'
    !
    implicit none

    ! Arguments:
    integer,          intent(in) :: blockIndex

    ! only the first thread does the timing
    if (omp_get_thread_num() > 0) return

    call tmg_stop(blockIndex)

  end subroutine rti_tmg_stop

  !--------------------------------------------------------------------------
  ! rti_writeStatus
  !--------------------------------------------------------------------------
  subroutine rti_writeStatus(cmsg)
    !
    !:Purpose: Create the file 'VAR3D_STATUS.dot' and write content to it
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: cmsg

    ! Locals:
    integer :: iulstatus, ierr
    character(len=22):: clmsg

    clmsg = 'VAR3D_STATUS='//cmsg
    iulstatus = 0
    ierr = fnom(iulstatus, 'VAR3D_STATUS.dot', 'SEQ+FMT', 0)
    rewind (iulstatus)
    WRITE(iulstatus,'(a22)') clmsg
    ierr = fclos(iulstatus)

  end subroutine rti_writeStatus

end module runtimeInfo_mod
