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
!! MODULE ramDisk_mod (prefix="ram")
!!
!! *Purpose*: Control the file manipulations/enquiries on the RAM disk
!!
!--------------------------------------------------------------------------
module ramDisk_mod
  !
  ! MODULE ramDisk_mod (prefix="ram")
  !
  ! *Purpose*: Control the file manipulations/enquiries on the RAM disk
  !
  !
  use utilities_mod
  use clib_interfaces_mod
  implicit none
  save
  private

  ! public procedures
  public :: ram_setup, ram_fullWorkingPath, ram_remove

  character(len=256) :: ram_disk_dir

  logical :: ram_disk_dir_exists
  logical :: initialized = .false.

contains

!--------------------------------------------------------------------------
! ram_setup
!--------------------------------------------------------------------------
  subroutine ram_setup
    implicit none

    integer :: length_ram_disk_dir, status

    !
    !- Determine the ramdisk directory, if available
    !
    status = 0
    call get_environment_variable('OAVAR_RAM_DISK_DIR',ram_disk_dir,length_ram_disk_dir,status,.true.)

    if (status.gt.1) then
      write(*,*) 'ram_setup: Problem when getting the environment variable OAVAR_RAM_DISK_DIR'
    end if
    if (status.eq.1) then
      write(*,*) 'ram_setup: The environment variable OAVAR_RAM_DISK_DIR has not been detected!'
      write(*,*) '           Assume all files in current working directory'
      ram_disk_dir_exists = .false.  
      ram_disk_dir = 'DOES_NOT_EXIST'
    else
      write(*,*)
      write(*,*) 'ram_setup: The environment variable OAVAR_RAM_DISK_DIR has correctly been detected'
      write(*,*) 'ram_setup: Files will first be opened from directory: ', trim(ram_disk_dir)
      ram_disk_dir_exists = .true.
    end if

    initialized = .true.

  end subroutine ram_setup

!--------------------------------------------------------------------------
! ram_fullWorkingPath - given a filename, return the full path by either
!                       adding the current working directory or the ram disk
!                       directory
!--------------------------------------------------------------------------
  function ram_fullWorkingPath(fileName, noAbort_opt) result(fullWorkingPath)
    implicit none
    character(len=512) :: fullWorkingPath
    logical, optional  :: noAbort_opt
    character(len=*)   :: fileName

    logical            :: fileExists, noAbort
    character(len=256) :: fileName2, subDirectory
    integer            :: status, fileSize

    if ( .not. initialized ) then
      call utl_abort('ram_fullWorkingPath: ramDisk module has not been initialized.')
    end if

    if ( present(noAbort_opt) ) then
      noAbort = noAbort_opt
    else
      noAbort = .false.
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
      ! now look in working directory
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
        status = copyFile(trim(fileName2), trim(ram_disk_dir) // '/' // trim(fileName2))

        fullWorkingPath = trim(ram_disk_dir) // '/' // trim(fileName2)
        fileSize = clib_size(trim(fullWorkingPath))
        write(*,*) 'ram_fullWorkingPath: file copied to ramdisk: ', trim(fullWorkingPath)
        write(*,*) 'ram_fullWorkingPath: size of copied file = ', fileSize

      end if

    end if

  end function ram_fullWorkingPath

!--------------------------------------------------------------------------
! ram_remove - given the full path+filename, remove the file only if 
!                     it is located on the ram disk (to free up memory)
!--------------------------------------------------------------------------
  function ram_remove(fullWorkingPath) result(returnCode)
    implicit none
    character(len=*) :: fullWorkingPath
    integer          :: returnCode
    logical          :: fileExists

    if ( .not. initialized ) then
      call utl_abort('ram_remove: ramDisk module has not been initialized.')
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
! copyFile - copy the specified file to the new location and/or name
!            This function is very general, but was initially written to
!            copy files from the disk to the ram disk
!--------------------------------------------------------------------------
  function copyFile(filein, fileout) result(status)
    implicit none
    character(len=*) :: filein
    character(len=*) :: fileout
    integer :: status

    integer :: ierr, unitin, unitout, numBuffer, numChar
    integer, parameter :: bufferSize = 1024
    character :: buffer(bufferSize), onechar

    call tmg_start(100,'CopyFile')

    write(*,*) 'copyFile: copy from ', trim(filein), ' to ', trim(fileout)

    unitin=10
    open(unit=unitin, file=trim(filein), status='OLD', form='UNFORMATTED', action='READ', access='STREAM')
    unitout=11
    open(unit=unitout, file=trim(fileout), status='NEW', form='UNFORMATTED', action='WRITE', access='STREAM')

    numBuffer = 0
    do 
      read(unitin,iostat=ierr) buffer
      if (ierr < 0) exit
      numBuffer = numBuffer + 1
      write(unitout) buffer
    end do
    numChar = numBuffer*bufferSize

    do 
      read(unitin,iostat=ierr,pos=numChar+1) oneChar
      if (ierr < 0) exit
      numChar = numChar + 1
      write(unitout) oneChar
    end do
    write(*,*) 'copyFile: copied ', numChar, ' bytes'

    close(unit=unitin)
    close(unit=unitout)

    if (numChar > 0) then
      status = 0
    else
      status = -1
      call utl_abort('ramdisk_mod copyFile: ERROR, zero bytes copied')
    end if

    call tmg_stop(100)

  end function copyFile

end module ramDisk_mod
