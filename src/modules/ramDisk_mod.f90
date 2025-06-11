
module ramDisk_mod
  ! MODULE ramDisk_mod (prefix='ram' category='8. Low-level utilities and constants')
  !
  !:Purpose: To control the file manipulations/enquiries on the RAM disk.
  !
  use utilities_mod
  use clibInterfaces_mod

  implicit none
  save
  private

  ! Public procedures
  public :: ram_setup, ram_fullWorkingPath, ram_remove, ram_getRamDiskDir
  public :: ram_removeRamDiskFromName, ram_fileIsOnRamDisk

  ! Public constants:
  integer, public, parameter :: ram_fileNameLength = 256
  integer, public, parameter :: ram_fullPathLength = 512

  character(len=ram_fileNameLength) :: ram_disk_dir

  logical :: ram_disk_dir_exists
  logical :: initialized = .false.

contains

  !--------------------------------------------------------------------------
  ! ram_setup
  !--------------------------------------------------------------------------
  subroutine ram_setup()
    implicit none

    ! Locals:
    integer :: length_ram_disk_dir, status

    !
    !- Determine the ramdisk directory, if available
    !
    status = 0
    call get_environment_variable('MIDAS_RAMDISKDIR',ram_disk_dir,length_ram_disk_dir,status,.true.)

    if (status.gt.1) then
      write(*,*) 'ram_setup: Problem when getting the environment variable MIDAS_RAMDISKDIR'
    end if
    if (status.eq.1) then
      write(*,*) 'ram_setup: The environment variable MIDAS_RAMDISKDIR has not been detected!'
      write(*,*) '           Assume all files in current working directory'
      ram_disk_dir_exists = .false.
      ram_disk_dir = 'DOES_NOT_EXIST'
    else
      write(*,*)
      write(*,*) 'ram_setup: The environment variable MIDAS_RAMDISKDIR has correctly been detected'
      write(*,*) 'ram_setup: Files will first be opened from directory: ', trim(ram_disk_dir)
      ram_disk_dir_exists = .true.
    end if

    initialized = .true.

  end subroutine ram_setup

  !--------------------------------------------------------------------------
  ! ram_fullWorkingPath
  !--------------------------------------------------------------------------
  function ram_fullWorkingPath(fileName, noAbort_opt, copyToRamDisk_opt) result(fullWorkingPath)
    !
    !:Purpose: Given a filename, return the full path by either adding the
    !          current working directory or the ram disk directory. By default,
    !          will copy the file to the ram disk directory, if it exists.
    !
    implicit none

    ! Arguments:
    logical, optional, intent(in) :: noAbort_opt
    logical, optional, intent(in) :: copyToRamDisk_opt
    character(len=*) , intent(in) :: fileName
    ! Result:
    character(len=ram_fullPathLength) :: fullWorkingPath

    ! Locals:
    character(len=ram_fileNameLength) :: fileName2, subDirectory
    logical            :: fileExists, noAbort, copyToRamDisk
    integer            :: status

    if ( .not. initialized ) then
      call ram_setup()
    end if

    if ( present(noAbort_opt) ) then
      noAbort = noAbort_opt
    else
      noAbort = .false.
    end if

    if ( present(copyToRamDisk_opt) ) then
      copyToRamDisk = copyToRamDisk_opt
    else
      copyToRamDisk = .true.
    end if

    ! this should make it safe for calls where input and output are the same variable
    fileName2 = trim(fileName)

    ! first look for file in the ram disk directory
    if ( ram_disk_dir_exists ) then
      fullWorkingPath = trim(ram_disk_dir) // '/' // trim(fileName2)
    else
      fullWorkingPath = trim(fileName2)
    end if
    inquire(file=trim(fullWorkingPath),exist=fileExists)

    ! treat case when no ramdisk exists
    if ( .not. ram_disk_dir_exists ) then
      if ( fileExists ) then
        return
      else
        if ( noAbort ) then
          fullWorkingPath = ' '
          return
        else
          write(*,*) 'ram_fullWorkingPath: file name = ', trim(fileName2)
          call utl_abort('ram_fullWorkingPath: this file cannot be found on disk.')
        end if
      end if
    end if

    ! treat the case when ramdisk DOES exists
    if ( fileExists ) then
      write(*,*) 'ram_fullWorkingPath: this file found on ram disk: ', trim(fileName2)
    else
      ! not found on ram disk, so now look in working directory
      inquire(file='./' // trim(fileName2),exist=fileExists)

      if ( .not. fileExists ) then
        if ( noAbort ) then
          fullWorkingPath = ' '
        else
          write(*,*) 'ram_fullWorkingPath: file name          = ', trim(fileName2)
          write(*,*) 'ram_fullWorkingPath: ram disk directory = ', trim(ram_disk_dir)
          call utl_abort('ram_fullWorkingPath: this file cannot be found.')
        end if
      else

        if ( copyToRamDisk ) then
          ! copy the file from disk to the ramdisk
          if ( index(trim(filename2),'/') /= 0 ) then
            status = clib_dirname(trim(filename2),subDirectory)
            status = clib_mkdir_r(trim(ram_disk_dir) // '/' // trim(subDirectory))
            if ( status /= clib_ok ) then
              call utl_abort('ram_fullWorkingPath: problem with mkdir')
            end if
            status = clib_isdir(trim(ram_disk_dir) // '/' // trim(subDirectory))
            if ( status /= clib_ok ) then
              call utl_abort('ram_fullWorkingPath: problem with checking existence of directory')
            end if
          end if

          ! copy the file from disk to the ramdisk
          status = utl_copyFile(trim(fileName2), trim(ram_disk_dir) // '/' // trim(fileName2))

          fullWorkingPath = trim(ram_disk_dir) // '/' // trim(fileName2)
          write(*,*) 'ram_fullWorkingPath: file copied to ramdisk: ', trim(fullWorkingPath)
        else
          fullWorkingPath = trim(fileName2)
          write(*,*) 'ram_fullWorkingPath: file left on disk, as requested: ', trim(fullWorkingPath)
        end if

      end if

    end if

  end function ram_fullWorkingPath

  !--------------------------------------------------------------------------
  ! ram_remove
  !--------------------------------------------------------------------------
  function ram_remove(fullWorkingPath) result(returnCode)
    !
    !:Purpose:  Given the full path+filename, remove the file only if
    !           it is located on the ram disk (to free up memory)
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: fullWorkingPath
    ! Result:
    integer          :: returnCode

    ! Locals:
    logical          :: fileExists

    if ( .not. initialized ) then
      call ram_setup()
    end if

    inquire(file=trim(fullWorkingPath),exist=fileExists)
    if ( .not. fileExists) then
      write(*,*) 'ram_Remove: file does not exist: ',trim(fullWorkingPath)
      returnCode = 0
      return
    end if

    if ( .not. ram_disk_dir_exists ) then
      write(*,*) 'ram_remove: no ram disk in use.'
      returnCode = 0
      return
    end if

    if ( index(trim(fullWorkingPath), trim(ram_disk_dir)) == 1 ) then
      write(*,*) 'ram_remove: removing file that is on the ram disk: ', trim(fullWorkingPath)
      returnCode = clib_remove(fullWorkingPath)
    else
      write(*,*) 'ram_remove: this file not on ram disk: ', trim(fullWorkingPath)
      returnCode = 0
    end if

  end function ram_remove

  !--------------------------------------------------------------------------
  ! ram_fileIsOnRamDisk
  !--------------------------------------------------------------------------
  function ram_fileIsOnRamDisk(fileName) result(fileIsOnRamDisk)
    !
    !:Purpose:  Given the full path+filename, determine if this file is on
    !           the ram disk by simply checking if its path contains the ram Disk path
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: fileName
    ! Result:
    logical                      :: fileIsOnRamDisk

    if ( .not. initialized ) then
      call ram_setup()
    end if

    ! treat case when no ramdisk exists
    if (.not. ram_disk_dir_exists) then
      fileIsOnRamDisk = .false.
      return
    end if

    ! check if path contains ram disk path
    if (index(trim(fileName), trim(ram_disk_dir)) == 1) then
      write(*,*) 'ram_fileIsOnRamDisk: file path contains ram disk path'
      fileIsOnRamDisk = .true.
    else
      fileIsOnRamDisk = .false.
    end if

  end function ram_fileIsOnRamDisk

  !--------------------------------------------------------------------------
  ! ram_removeRamDiskFromName
  !--------------------------------------------------------------------------
  function ram_removeRamDiskFromName(fileName) result(fileNameWithoutRamDisk)
    !
    !:Purpose:  Given the full path+filename, remove the ram disk path from the name.
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: fileName
    ! Result:
    character(len=ram_fileNameLength) :: fileNameWithoutRamDisk

    if ( .not. initialized ) then
      call ram_setup()
    end if

    ! treat case when no ramdisk exists
    if (.not. ram_disk_dir_exists) then
      fileNameWithoutRamDisk = trim(fileName)
      return
    end if

    ! check if path contains ram disk path
    if (index(trim(fileName), trim(ram_disk_dir)) == 1) then
      write(*,*) 'ram_removeRamDiskFromName: file path contains ram disk path'
      fileNameWithoutRamDisk = fileName(len_trim(ram_disk_dir)+2 : len_trim(fileName))
      write(*,*) 'ram_removeRamDiskFromName: Name with and without ram disk path = '
      write(*,*) '                           With:    ', trim(fileName)
      write(*,*) '                           Without: ', trim(fileNameWithoutRamDisk)
    else
      fileNameWithoutRamDisk = trim(fileName)
    end if

  end function ram_removeRamDiskFromName

  !--------------------------------------------------------------------------
  ! ram_getRamDiskDir
  !--------------------------------------------------------------------------
  function ram_getRamDiskDir() result(fullWorkingPath)
    !
    !:Purpose:  Return the ram disk path.
    !
    implicit none

    ! Result:
    character(len=ram_fullPathLength) :: fullWorkingPath ! Ram disk path returned

    if ( initialized ) then
      if ( ram_disk_dir_exists ) then
        fullWorkingPath = trim(ram_disk_dir) // '/'
      else
        fullWorkingPath = ' '
      end if
    else
      fullWorkingPath = ' '
    end if

  end function ram_getRamDiskDir

end module ramDisk_mod
