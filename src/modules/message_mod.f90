module message_mod
  ! MODULE message_mod (prefix='msg' category='8. Low-level utilities and constants')
  !
  ! :Purpose: Output message interface with configurable verbosity.
  !           Also provides string representation for some intrinsic types.
  !
  use midasMpi_mod
  use utilities_mod
  implicit none
  save
  private

  integer, external :: get_max_rss

  ! public procedures
  public :: msg, msg_memUsage, msg_setVerbThreshold

  ! intrinsic type string representations
  public :: str
  interface str
    module procedure msg_log2str
    module procedure msg_int2str
    module procedure msg_real42str
    module procedure msg_real82str
    module procedure msg_charArray2str
    module procedure msg_logArray2str
    module procedure msg_intArray2str
    module procedure msg_real4Array2str
    module procedure msg_real8Array2str
  end interface

  ! private module variables
  integer, parameter    :: msg_maxIndent = 10
  integer, parameter    :: msg_lineLen = 70
  integer, parameter    :: msg_num2strBufferLen = 200

  ! Namelist variables
  !-------------------
  integer :: verbosityThreshold

  contains
  
  !--------------------------------------------------------------------------
  ! msg
  !--------------------------------------------------------------------------
  subroutine msg(origin, message, verb_opt, mpiAll_opt)
    !
    ! :Purpose: Output message if its verbosity level is greater or equal than
    !           the user provided verbosity threshold (see `msg_readNml()`).
    !           The verbosity levels are:
    !
    !                             * 0 : critical, always printed
    !                             * 1 : default priority; printed in operational context
    !                             * 2 : detailed ouptut, provides extra information
    !                             * 3 : intended for developpers, printed for debugging or specific diagnostcs
    !
    implicit none

    ! Arguments:
    character(len=*),  intent(in) :: origin     ! originating subroutine, function or program
    character(len=*),  intent(in) :: message    ! message to be printed
    integer, optional, intent(in) :: verb_opt   ! minimal verbosity level to print the message
    logical, optional, intent(in) :: mpiAll_opt ! if `.true.` prints to all MPI tasks, otherwise only to tile 0

    ! Locals:
    logical :: mpiAll
    integer :: verbLevel

    call msg_readNml()

    verbLevel = 1
    if (present(verb_opt)) verbLevel = verb_opt
    mpiAll = .true.
    if (present(mpiAll_opt)) mpiAll = mpiAll_opt

    if (verbLevel <= verbosityThreshold) then
      if (mpiAll) then
        call msg_write(origin, message)
      else
        if (mmpi_myid == 0) call msg_write(origin, message)
      end if
    end if
    end subroutine msg

  !--------------------------------------------------------------------------
  ! msg_memUsage
  !--------------------------------------------------------------------------
  subroutine msg_memUsage(origin, verb_opt, mpiAll_opt)
    !
    ! :Purpose: Report memory usage
    !
    implicit none

    ! Arguments:
    character(len=*),  intent(in) :: origin     ! originating subroutine, function or program
    integer, optional, intent(in) :: verb_opt   ! minimal verbosity level to print the message
    logical, optional, intent(in) :: mpiAll_opt ! if `.true.` prints to all MPI tasks, otherwise only to tile 0

    ! Locals:
    integer :: usageMb

    usageMb = get_max_rss()/1024
    call msg( origin, "Memory Used: "//str(usageMb)//" Mb", verb_opt, mpiAll_opt)

  end subroutine msg_memUsage

  !--------------------------------------------------------------------------
  ! msg_setVerbThreshold
  !--------------------------------------------------------------------------
  subroutine msg_setVerbThreshold(threshold_opt)
    !
    ! :Purpose: Sets the verbosity level at runtime, overrides the namelist
    !           value if set.
    !
    implicit none

    ! Arguments:
    integer, optional, intent(in) :: threshold_opt

    ! Locals:
    integer :: threshold

    if (present(threshold_opt)) then
      threshold = threshold_opt
    else
      threshold = 1
    end if
    call msg( 'msg_setVerbThreshold', 'Setting verbosity threshold to '&
              //str(threshold), verb_opt=3)
    verbosityThreshold = threshold

  end subroutine

  !--------------------------------------------------------------------------
  ! msg_readNml (private)
  !--------------------------------------------------------------------------
  subroutine msg_readNml()
    !
    ! :Purpose: Reads the module configuration namelist
    !
    ! :Namelist parameters:
    !       :verbosityThreshold:  Each call to `msg()` specifies a verbosity
    !                             level; this threshold configures until which level
    !                             messages will be outputed.
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
        write(*,*) 'msg: NAMMSG is missing in the namelist.'
        write(*,*) '             The default values will be taken.'
      end if
    else
      nulnam=0
      ierr=fnom(nulnam, 'flnml','FTN+SEQ+R/O',0)
      read(nulnam,nml=nammsg,iostat=ierr)
      if (ierr /= 0) call utl_abort('msg: Error reading namelist NAMMSG')
      if (mpi_myid == 0) write(*,nml=nammsg)
      ierr = fclos(nulnam)
    end if

  end subroutine msg_readNml

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
    integer :: originLen, oneLineMsgLen, indentLen, i
    character(len=15) :: firstLineFormat, otherLineFormat
    character(len=msg_lineLen)  :: msgLine
    character(len=msg_lineLen)  :: readLine
    character(len=:), allocatable :: originTrunc

    if (len(origin) > msg_lineLen) then
      originTrunc = origin(1:msg_lineLen)
    else
      originTrunc = origin
    end if
    originLen = len(originTrunc)

    indentLen = originLen
    if (originLen > msg_maxIndent) then
      indentLen = msg_maxIndent
    end if
    oneLineMsgLen = msg_lineLen - originLen - 2
  
    if (len(message) > oneLineMsgLen) then
      ! Multiple lines message
      ! format: "origin: message on the first line..........."
      !         "        second line........................."
      !         "        last line"
      i = 0
      readLine = message(1:oneLineMsgLen+1)
      msgLine = msg_breakOnSpace(readLine)
      i = i + len(trim(msgLine)) +1
      write(firstLineFormat,'(A,I2,A,I2,A)') '(A',originLen,',A2,A', &
                                              len(trim(msgLine)),')'
      write(*,firstLineFormat) originTrunc, ': ', message(1:oneLineMsgLen)
      oneLineMsgLen = msg_lineLen - indentLen - 2
      do
        if ( i >= len(message) ) then
          ! message printed
          return
        else if ( i + oneLineMsgLen > len(trim(message)) ) then
          ! last line
          msgLine = message(i+1:len(message))
        else
          ! neither first nor last
          readLine = message(i+1:i+oneLineMsgLen+1)
          msgLine = msg_breakOnSpace(readLine)
        end if
        i = i + len(trim(msgLine)) +1
        write(otherLineFormat,'(A,I2,A,I2,A)') '(A',indentLen+2,',A', &
                                                len(trim(msgLine)),')'
        write(*,otherLineFormat) repeat(' ',indentLen+2),trim(msgLine)
      end do
    else
      ! Single lines message
      ! format: "origin: short message"
      write(firstLineFormat,'(A,I2,A,I2,A)') '(A',originLen,',A2,A',len(message),')'
      write(*,firstLineFormat) originTrunc, ': ', message
    end if
    contains
      !----------------------------------------------------------------------
      ! msg_breakOnSpace (private)
      !----------------------------------------------------------------------
      function msg_breakOnSpace(line) result(shorterLine)
        !
        ! :Purpose: Breaks line on last full word
        !
        implicit none
    
        ! Arguments:
        character(len=msg_lineLen), intent(in)  :: line
        character(len=msg_lineLen) :: shorterLine
        integer :: i
    
        i = index(trim(line),' ',back=.true.)
        if (i == 0 .or. i == len(trim(line)) ) then
          shorterLine = trim(line)
          return
        else
          shorterLine = line(1:i-1)
          return
        end if
    
      end function msg_breakOnSpace
  end subroutine msg_write

  !--------------------------------------------------------------------------
  ! msg_log2str (private)
  !--------------------------------------------------------------------------
  function msg_log2str(num) result(string)
    !
    ! :Purpose: Returns string representation of `logical`
    !
    implicit none

    ! Arguments:
    logical, intent(in)           :: num
    character(len=:), allocatable :: string

    ! Locals:
    character(len=msg_num2strBufferLen) :: buffer

    write(buffer,*) num
    string = trim(adjustl(buffer))

  end function msg_log2str

  !--------------------------------------------------------------------------
  ! msg_int2str (private)
  !--------------------------------------------------------------------------
  function msg_int2str(num) result(string)
    !
    ! :Purpose: Returns string representation of `integer`
    !
    implicit none

    ! Arguments:
    integer, intent(in)           :: num
    character(len=:), allocatable :: string

    ! Locals:
    character(len=msg_num2strBufferLen) :: buffer

    write(buffer,*) num
    string = trim(adjustl(buffer))

  end function msg_int2str

  !--------------------------------------------------------------------------
  ! msg_real42str (private)
  !--------------------------------------------------------------------------
  function msg_real42str(num, digits_opt) result(string)
    !
    ! :Purpose: Returns string representation of `real(4)`
    !
    implicit none

    ! Arguments:
    real(4), intent(in)           :: num
    character(len=:), allocatable :: string
    integer, optional             :: digits_opt

    ! Locals:
    character(len=20)                   :: readFmt, digitBuffer
    character(len=msg_num2strBufferLen) :: buffer

    if (present(digits_opt)) then
      write(digitBuffer, *) digits_opt
      write(readFmt,*) '(F20.'//trim(adjustl(digitBuffer))//')' 
      write(buffer, readFmt) num
    else
      write(buffer,*) num
    end if
    string = trim(adjustl(buffer))

  end function msg_real42str

  !--------------------------------------------------------------------------
  ! msg_real82str (private)
  !--------------------------------------------------------------------------
  function msg_real82str(num, digits_opt) result(string)
    !
    ! :Purpose: Returns string representation of `real(8)` 
    !
    implicit none

    ! Arguments:
    real(8), intent(in)           :: num
    character(len=:), allocatable :: string
    integer, optional             :: digits_opt

    ! Locals:
    character(len=20)                   :: readFmt, digitBuffer
    character(len=msg_num2strBufferLen) :: buffer

    if (present(digits_opt)) then
      write(digitBuffer, *) digits_opt
      write(readFmt,*) '(F20.'//trim(adjustl(digitBuffer))//')' 
      write(buffer, readFmt) num
    else
      write(buffer,*) num
    end if
    string = trim(adjustl(buffer))

  end function msg_real82str

  !--------------------------------------------------------------------------
  ! msg_charArray2str (private)
  !--------------------------------------------------------------------------
  function msg_charArray2str(array, vertical_opt) result(string)
    !
    ! :Purpose: Returns string representation of `character(len=*), dimension(:)` 
    !
    implicit none

    ! Arguments
    character(len=*), dimension(:), intent(in) :: array
    character(len=:), allocatable     :: string
    logical, optional                 :: vertical_opt

    ! Locals:
    integer           :: i
    logical           :: vertical
    character(len=2)  :: sep

    vertical=.false.
    if (present(vertical_opt)) vertical = vertical_opt
    if (vertical) then
      sep=new_line('a')
      string='(/'//sep
    else
      sep=', '
      string='(/ '
    end if

    do i=1,size(array)
      string = string//array(i)
      if (i /= size(array)) string = string//sep
    end do
    string = string//' /)'
    
  end function msg_charArray2str

  !--------------------------------------------------------------------------
  ! msg_intArray2str (private)
  !--------------------------------------------------------------------------
  function msg_logArray2str(array, vertical_opt) result(string)
    !
    ! :Purpose: Returns string representation of `logical, dimension(:)` 
    !
    implicit none

    ! Arguments
    logical, dimension(:), intent(in) :: array
    character(len=:), allocatable     :: string
    logical, optional                 :: vertical_opt

    ! Locals:
    integer           :: i
    logical           :: vertical
    character(len=2)  :: sep

    vertical=.false.
    if (present(vertical_opt)) vertical = vertical_opt
    if (vertical) then
      sep=new_line('a')
      string='(/'//sep
    else
      sep=', '
      string='(/ '
    end if

    do i=1,size(array)
      string = string//msg_log2str(array(i))
      if (i /= size(array)) string = string//sep
    end do
    string = string//' /)'
    
  end function msg_logArray2str

  !--------------------------------------------------------------------------
  ! msg_intArray2str (private)
  !--------------------------------------------------------------------------
  function msg_intArray2str(array, vertical_opt) result(string)
    !
    ! :Purpose: Returns string representation of `integer, dimension(:)` 
    !
    implicit none

    ! Arguments
    integer, dimension(:), intent(in) :: array
    character(len=:), allocatable     :: string
    logical, optional                 :: vertical_opt

    ! Locals:
    integer           :: i
    logical           :: vertical
    character(len=2)  :: sep

    vertical=.false.
    if (present(vertical_opt)) vertical = vertical_opt
    if (vertical) then
      sep=new_line('a')
      string='(/'//sep
    else
      sep=', '
      string='(/ '
    end if

    do i=1,size(array)
      string = string//msg_int2str(array(i))
      if (i /= size(array)) string = string//sep
    end do
    string = string//' /)'
    
  end function msg_intArray2str

  !--------------------------------------------------------------------------
  ! msg_real4Array2str (private)
  !--------------------------------------------------------------------------
  function msg_real4Array2str(array, digits_opt, vertical_opt) result(string)
    !
    ! :Purpose: Returns string representation of `real(4), dimension(:)` 
    !
    implicit none

    ! Arguments
    real(4), dimension(:), intent(in) :: array
    character(len=:), allocatable     :: string
    integer, optional                 :: digits_opt
    logical, optional                 :: vertical_opt

    ! Locals:
    integer           :: i
    logical           :: vertical
    character(len=2)  :: sep

    vertical=.false.
    if (present(vertical_opt)) vertical = vertical_opt
    if (vertical) then
      sep=new_line('')
      string='(/'//sep
    else
      sep=', '
      string='(/ '
    end if

    do i=1,size(array)
      string = string//msg_real42str(array(i), digits_opt=digits_opt)
      if (i /= size(array)) string = string//sep
    end do
    string = string//' /)'

    
  end function msg_real4Array2str

  !--------------------------------------------------------------------------
  ! msg_real8Array2str (private)
  !--------------------------------------------------------------------------
  function msg_real8Array2str(array, digits_opt, vertical_opt) result(string)
    !
    ! :Purpose: Returns string representation of `real(8), dimension(:)` 
    !
    implicit none

    ! Arguments
    real(8), dimension(:), intent(in) :: array
    character(len=:), allocatable     :: string
    integer, optional                 :: digits_opt
    logical, optional                 :: vertical_opt

    ! Locals:
    integer           :: i
    logical           :: vertical
    character(len=2)  :: sep

    vertical=.false.
    if (present(vertical_opt)) vertical = vertical_opt
    if (vertical) then
      sep=new_line('')
      string='(/'//sep
    else
      sep=', '
      string='(/ '
    end if

    do i=1,size(array)
      string = string//msg_real82str(array(i), digits_opt=digits_opt)
      if (i /= size(array)) string = string//sep
    end do
    string = string//' /)'
    
  end function msg_real8Array2str
end module
