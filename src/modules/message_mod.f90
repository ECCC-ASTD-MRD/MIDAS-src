module message_mod
  ! MODULE message_mod (prefix='msg' category='8. Low-level utilities and constants')
  !
  ! :Purpose: Ouput message interface with configurable verbosity
  !
  use mpi_mod
  use utilities_mod
  implicit none
  save
  private

  ! public procedures
  public :: msg_message, msg_memUsage
  
  ! private variables ! DBGmad move to namelist?
  integer, parameter    :: msg_maxOriginLen = 30
  integer, parameter    :: msg_lineLen = 70

  ! Namelist variables
  !-------------------
  integer :: verbosityThreshold

  contains
  
  !--------------------------------------------------------------------------
  ! msg_readNml (private)
  !--------------------------------------------------------------------------
  subroutine msg_readNml()
    !
    ! :Purpose: Reads the module configuration namelist
    !
    ! :Namelist parameters:
    !       :verbosityThreshold:  Define until which verbosity level messages
    !                             are outputed
    !
    implicit none

    ! Locals:
    logical, save :: alreadyRead = .false.
    integer :: nulnam, ierr, fnom, fclos
    namelist /NAMMSG/verbosityThreshold
  
    if (alreadyRead) then
      return
    else
      alreadyRead = .true.
    end if
  
    ! default namelist value
    verbosityThreshold = 1
  
    if ( .not. utl_isNamelistPresent('NAMMSG','./flnml') ) then
      if ( mpi_myid == 0 ) then
        write(*,*) 'msg_message: NAMMSG is missing in the namelist.'
        write(*,*) '             The default values will be taken.'
      end if
    else
      nulnam=0
      ierr=fnom(nulnam, 'flnml','FTN+SEQ+R/O',0)
      read(nulnam,nml=nammsg,iostat=ierr)
      if (ierr /= 0) call utl_abort('msg_message: Error reading namelist NAMMSG')
      if (mpi_myid == 0) write(*,nml=nammsg)
      ierr = fclos(nulnam)
    end if

  end subroutine msg_readNml

  !--------------------------------------------------------------------------
  ! msg_message
  !--------------------------------------------------------------------------
  subroutine msg_message(origin, message, verbosityLevel_opt, mpiAll_opt)
    !
    ! :Purpose: output message if its verbosity level is greater or equal than
    !           the user provided verbosity threshold
    !
    implicit none
  
    ! Arguments:
    character(len=*),  intent(in) :: origin             ! originating subroutine/function
    character(len=*),  intent(in) :: message            ! message to be printed
    integer, optional, intent(in) :: verbosityLevel_opt ! minimal verbosity level to print the message
    logical, optional, intent(in) :: mpiAll_opt         ! if .true. prints to all MPI tasks, otherwise only to tile 0
    
    ! Locals:
    logical :: mpiAll = .false.
    integer :: verbLevel = 1
  
    call msg_readNml()

    if (present(verbosityLevel_opt)) verbLevel = verbosityLevel_opt
    if (present(mpiAll_opt)) mpiAll = mpiAll_opt
  
    if (verbLevel >= verbosityThreshold) then
      if (mpiAll) then
        call msg_write(origin, message)
      else
        if (mpi_myid == 0) call msg_write(origin, message)
      end if
    end if
    end subroutine msg_message

  !--------------------------------------------------------------------------
  ! msg_write (private)
  !--------------------------------------------------------------------------
  subroutine msg_write(origin, message)
    !
    ! :Purpose: Format and write message to default output
    !
    implicit none
  
    ! Arguments:
    character(len=*), intent(in) :: origin     ! originating subroutine/function
    character(len=*), intent(in) :: message    ! message to be printed
  
    ! Locals:
    integer :: originLen, oneLineMsgLen, i
    character(len=15) :: firstLineFormat, otherLineFormat
    character(len=msg_lineLen)  :: msgLine

    if (len(origin) > msg_maxOriginLen) then
      call utl_abort('DBGmad : fix that!')
    end if
    originLen = len(origin)
    oneLineMsgLen = msg_lineLen - originLen - 2
  
    i = 0 
    if (len(message) > oneLineMsgLen) then
      ! format: "origin: message on the first line..........."
      !         "        second line........................."
      !         "        last line"
      write(firstLineFormat,'(A,I2,A,I2,A)') '(A',originLen,',A2,A',oneLineMsgLen,')'

      write(*,firstLineFormat) origin, ': ', message(1:oneLineMsgLen)
      do
        if ( (i+1)*oneLineMsgLen > len(message) ) then
          ! message printed
          exit
        else if ( (i+2)*oneLineMsgLen > len(message) ) then
          ! last line
          write(otherLineFormat,'(A,I2,A,I2,A)') '(A',originLen+2,',A', &
               len(message)-(i+1)*oneLineMsgLen,')'
        else
          ! neither first nor last
          write(otherLineFormat,'(A,I2,A,I2,A)') '(A',originLen+2,',A',oneLineMsgLen,')'
        end if
        i = i + 1
        msgLine = message(i*oneLineMsgLen+1:(i+1)*oneLineMsgLen+1)
        write(*,otherLineFormat) repeat(' ',originLen+2),msgLine
      end do
    else
      write(firstLineFormat,'(A,I2,A,I2,A)') '(A',originLen,',A2,A',len(message),')'

      write(*,firstLineFormat) origin, ': ', message
    end if
    
  end subroutine msg_write

  !--------------------------------------------------------------------------
  ! msg_memUsage
  !--------------------------------------------------------------------------
  subroutine msg_memUsage()
    !
    ! :Purpose: 
    !

    !call msg_message()

  end subroutine msg_memUsage
end module
