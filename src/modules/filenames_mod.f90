
module fileNames_mod
  ! MODULE fileNames_mod (prefix='fln' category='7. Low-level data objects')
  !
  ! :Purpose: Routines related to file names
  !
  use utilities_mod
  use clib_interfaces_mod
  use ramDisk_mod
  use timeCoord_mod
  implicit none
  save
  private

  ! public procedures
  public :: fln_ensFileName, fln_ensAnlFileName, fln_ensTrlFileName

contains

  !--------------------------------------------------------------------------
  ! fln_ensFileName
  !--------------------------------------------------------------------------
  subroutine fln_ensFileName(ensFileName, ensPathName, memberIndex_opt, ensFileNamePrefix_opt,  &
                             ensFileBaseName_opt, shouldExist_opt, ensembleFileExtLength_opt, &
                             copyToRamDisk_opt, resetFileInfo_opt, fileMemberIndex1_opt)
    ! :Purpose: Return the filename of an ensemble member. Will also call routine in 
    !           ramdisk_mod module that will copy the file (if shouldExist_opt is true)
    !           to the ram disk. If the memberIndex_opt is not specified, the filename
    !           is returned without the member index extension (used to read deterministic
    !           background state that is stored in the ensemble directory for LETKF).
    !
    implicit none

    ! Arguments:
    character(len=*)  :: ensFileName
    character(len=*)  :: ensPathName
    integer, optional :: memberIndex_opt
    character(len=*),optional  :: ensFileBaseName_opt, ensFileNamePrefix_opt
    logical, optional :: shouldExist_opt
    integer, optional :: ensembleFileExtLength_opt
    logical, optional :: copyToRamDisk_opt
    logical, optional :: resetFileInfo_opt
    integer, optional :: fileMemberIndex1_opt

    ! Locals:
    integer          :: numFiles, returnCode, totalLength, ensembleBaseFileNameLength
    character(len=10):: ensNumber  !! this is sufficient until we reach 10^10 members
    logical          :: shouldExist
    character(len=2000) :: fileList(10), fileNamePattern
    character        :: ensembleFileExtLengthStr
    logical, save    :: firstTime = .true.
    integer, save    :: ensembleFileExtLength = 4
    integer, save    :: fileMemberIndex1 = 1
    character(len=4),   save :: fileMemberIndex1Str
    character(len=200), save :: ensFileBaseName

    if ( present(resetFileInfo_opt) ) then
      if (resetFileInfo_opt) firstTime = .true.
      return
    end if

    if ( present(shouldExist_opt) ) then
      shouldExist = shouldExist_opt
    else
      shouldExist = .true.
    end if

    ! Do this step only once in the program since this should not change during the program is running.
    if (firstTime) then
      if (present(fileMemberIndex1_opt)) fileMemberIndex1 = fileMemberIndex1_opt
      if (fileMemberIndex1 > 999) then
        write(fileMemberIndex1Str,'(i4.4)') fileMemberIndex1
      else
        write(fileMemberIndex1Str,'(i3.3)') fileMemberIndex1
      end if

      write(*,*) 'fln_ensFileName: looking for ./' // trim(enspathname) // '/' // '*_*' // trim(fileMemberIndex1Str)
      fileNamePattern = './' // trim(enspathname) // '/' // '*_*' // trim(fileMemberIndex1Str)
      returnCode = clib_glob(fileList,numFiles,trim(fileNamePattern),10)
      if (returnCode /= 1) then
        if (shouldExist) then
          call utl_abort('fln_ensFileName: reached maximum number of files or no file is available')
        else
          write(*,*) 'fln_ensFileName: reached maximum number of files or no file is available'
          ensFileName='missing'
          return
        end if
      end if

      ensFileName = trim(fileList(1))
      totalLength = len_trim(ensFileName)
      if ( totalLength == 0 ) then
        call utl_abort('fln_ensFileName: ensFileName seems empty: ''ensFileName=' // trim(ensFileName) // '''')
      end if

      ! find number of digits used to identify ensemble member index
      ensembleFileExtLength = 0
      do
        if ( ensFileName((totalLength-ensembleFileExtLength):(totalLength-ensembleFileExtLength)) == '_' ) exit
        ensembleFileExtLength = ensembleFileExtLength + 1
      end do

      if (ensembleFileExtLength == 0) then
        call utl_abort('fln_ensFileName: Cannot determine the ensemble file extention length with ' // trim(ensFileName))
      end if

      ! find the last '/' in the file name to get the basename of the file
      ensembleBaseFileNameLength = 0
      do
        if ( totalLength == ensembleBaseFileNameLength ) exit
        if ( ensFileName((totalLength-ensembleBaseFileNameLength):(totalLength-ensembleBaseFileNameLength)) == '/' ) exit
        ensembleBaseFileNameLength = ensembleBaseFileNameLength + 1
      end do

      ! if 'ensFileName = ./abc/def/ghi/123_456_001' then
      !    totalLength = 25
      !    ensembleFileExtLength = 3
      !    ensembleBaseFileNameLength = 11
      !    ensFileBaseName = '123_456'
      ! if 'ensFileName = 123_456_001' then
      !    totalLength = 11
      !    ensembleFileExtLength = 3
      !    ensembleBaseFileNameLength = 11
      !    ensFileBaseName = '123_456'
      ensFileBasename = ensFileName((totalLength-ensembleBaseFileNameLength+1):(totalLength-ensembleFileExtLength-1))

      firstTime = .false.
    end if

    if (present(memberIndex_opt)) then
      write(ensembleFileExtLengthStr,'(i1.1)') ensembleFileExtLength
      write(ensNumber,'(i' // ensembleFileExtLengthStr // '.' // ensembleFileExtLengthStr // ')') &
          memberIndex_opt + fileMemberIndex1 - 1
    end if

    if (present(memberIndex_opt)) then
      if (present(ensFileNamePrefix_opt)) then
        ensFileName = trim(enspathname) // '/' // trim(ensFileNamePrefix_opt) //  &
                      trim(ensFileBaseName) // '_' // trim(ensNumber)
      else
        ensFileName = trim(enspathname) // '/' // trim(ensFileBaseName) // '_' // trim(ensNumber)
      end if
    else
      if (present(ensFileNamePrefix_opt)) then
        ensFileName = trim(enspathname) // '/' // trim(ensFileNamePrefix_opt) // trim(ensFileBaseName)
      else
        ensFileName = trim(enspathname) // '/' // trim(ensFileBaseName)
      end if
    end if

    write(*,*) 'fln_ensFileName: ensFileName = ', trim(ensFileName)

    if ( shouldExist ) ensFileName = ram_fullWorkingPath(ensFileName, copyToRamDisk_opt=copyToRamDisk_opt)

    if (present(ensFileBaseName_opt)) ensFileBaseName_opt = trim(ensFileBaseName)
    if (present(ensembleFileExtLength_opt)) ensembleFileExtLength_opt = ensembleFileExtLength

  end subroutine fln_ensFileName

  !--------------------------------------------------------------------------
  ! fln_ensAnlFileName
  !--------------------------------------------------------------------------
  subroutine fln_ensAnlFileName( ensFileName, ensPathName, dateStamp,  &
                                 memberIndex_opt, ensFileNamePrefix_opt, ensFileNameSuffix_opt )
    ! :Purpose: Return the filename for an analysis state, including for
    !           ensemble members (by specifying memberIndex_opt). The member
    !           index extension is assumed to be 4 digits.
    !
    implicit none

    ! arguments
    character(len=*)  :: ensFileName
    character(len=*)  :: ensPathName
    integer           :: dateStamp
    integer, optional :: memberIndex_opt
    character(len=*), optional :: ensFileNamePrefix_opt
    character(len=*), optional :: ensFileNameSuffix_opt

    ! locals
    integer :: imode, ierr, hours, prntdate, prnttime, newdate
    character(len=4)  :: ensNumber
    character(len=10) :: dateStrAnl

    ! Set the printable date for analysis related file names
    imode = -3 ! stamp to printable date and time: YYYYMMDD, HHMMSShh
    ierr = newdate(dateStamp, prntdate, prnttime, imode)
    hours = prnttime/1000000
    write(dateStrAnl,'(i10.10)') prntdate*100 + hours

    if (present(memberIndex_opt)) then
      write(ensNumber,'(i4.4)') memberIndex_opt
    end if

    if (present(memberIndex_opt)) then
      if (present(ensFileNamePrefix_opt)) then
        if (present(ensFileNameSuffix_opt)) then
          ensFileName = trim(ensPathName) // '/' // trim(ensFileNamePrefix_opt) //  &
                        dateStrAnl // '_000_' // trim(ensFileNameSuffix_opt) // '_' // trim(ensNumber)
        else
          ensFileName = trim(ensPathName) // '/' // trim(ensFileNamePrefix_opt) //  &
                        dateStrAnl // '_000_' // trim(ensNumber)
        end if
      else
        if (present(ensFileNameSuffix_opt)) then
          ensFileName = trim(ensPathName) // '/' // dateStrAnl // '_000_' //  &
                        trim(ensFileNameSuffix_opt) // '_' // trim(ensNumber)
        else
          ensFileName = trim(ensPathName) // '/' // dateStrAnl // '_000_' // trim(ensNumber)
        end if
      end if
    else
      if (present(ensFileNamePrefix_opt)) then
        if (present(ensFileNameSuffix_opt)) then
          ensFileName = trim(ensPathName) // '/' // trim(ensFileNamePrefix_opt) //  &
                        dateStrAnl // '_000_' // trim(ensFileNameSuffix_opt)
        else
          ensFileName = trim(ensPathName) // '/' // trim(ensFileNamePrefix_opt) //  &
                        dateStrAnl // '_000'
        end if
      else
        if (present(ensFileNameSuffix_opt)) then
          ensFileName = trim(ensPathName) // '/' // dateStrAnl // '_000_' //  &
                        trim(ensFileNameSuffix_opt)
        else
          ensFileName = trim(ensPathName) // '/' // dateStrAnl // '_000'
        end if
      end if
    end if

    write(*,*) 'fln_ensAnlFileName: ensFileName = ', trim(ensFileName)

  end subroutine fln_ensAnlFileName

  !--------------------------------------------------------------------------
  ! fln_ensTrlFileName
  !--------------------------------------------------------------------------
  subroutine fln_ensTrlFileName( ensFileName, ensPathName, dateStamp,  &
                                 memberIndex_opt, ensFileNamePrefix_opt, ensFileNameSuffix_opt )
    ! :Purpose: Return the filename for a trial state, including for
    !           ensemble members (by specifying memberIndex_opt). The member
    !           index extension is assumed to be 4 digits.
    !
    implicit none

    ! arguments
    character(len=*)  :: ensFileName
    character(len=*)  :: ensPathName
    integer           :: dateStamp
    integer, optional :: memberIndex_opt
    character(len=*), optional :: ensFileNamePrefix_opt
    character(len=*), optional :: ensFileNameSuffix_opt

    ! locals
    integer :: imode, ierr, hours, prntdate, prnttime, newdate, dateStampTrl
    character(len=4)  :: ensNumber
    character(len=3)  :: leadTimeStr
    character(len=10) :: dateStrTrl

    ! Compute the datestamp for the origin time of the trial forecasts
    call incdatr(dateStampTrl, dateStamp, -tim_windowsize)
    
    ! Set the printable date for trial related file names
    imode = -3 ! stamp to printable date and time: YYYYMMDD, HHMMSShh
    ierr = newdate(dateStampTrl, prntdate, prnttime, imode)
    hours = prnttime/1000000
    write(dateStrTrl,'(i10.10)') prntdate*100 + hours

    if (present(memberIndex_opt)) then
      write(ensNumber,'(i4.4)') memberIndex_opt
    end if
    write(leadTimeStr,'(i3.3)') nint(tim_windowsize)

    if (present(memberIndex_opt)) then
      if (present(ensFileNamePrefix_opt)) then
        if (present(ensFileNameSuffix_opt)) then
          ensFileName = trim(ensPathName) // '/' // trim(ensFileNamePrefix_opt) //  &
                        dateStrTrl // '_' // leadTimeStr // '_' // trim(ensFileNameSuffix_opt) // '_' // trim(ensNumber)
        else
          ensFileName = trim(ensPathName) // '/' // trim(ensFileNamePrefix_opt) //  &
                        dateStrTrl // '_' // leadTimeStr // '_' // trim(ensNumber)
        end if
      else
        if (present(ensFileNameSuffix_opt)) then
          ensFileName = trim(ensPathName) // '/' // dateStrTrl // '_' // leadTimeStr // '_' //  &
                        trim(ensFileNameSuffix_opt) // '_' // trim(ensNumber)
        else
          ensFileName = trim(ensPathName) // '/' // dateStrTrl // '_' // leadTimeStr // '_' // trim(ensNumber)
        end if
      end if
    else
      if (present(ensFileNamePrefix_opt)) then
        if (present(ensFileNameSuffix_opt)) then
          ensFileName = trim(ensPathName) // '/' // trim(ensFileNamePrefix_opt) //  &
                        dateStrTrl // '_' // leadTimeStr // '_' // trim(ensFileNameSuffix_opt)
        else
          ensFileName = trim(ensPathName) // '/' // trim(ensFileNamePrefix_opt) //  &
                        dateStrTrl // '_' // leadTimeStr
        end if
      else
        if (present(ensFileNameSuffix_opt)) then
          ensFileName = trim(ensPathName) // '/' // dateStrTrl // '_' // leadTimeStr // '_' //  &
                        trim(ensFileNameSuffix_opt)
        else
          ensFileName = trim(ensPathName) // '/' // dateStrTrl // '_' // leadTimeStr
        end if
      end if
    end if

    write(*,*) 'fln_ensTrlFileName: ensFileName = ', trim(ensFileName)

  end subroutine fln_ensTrlFileName

end module fileNames_mod
