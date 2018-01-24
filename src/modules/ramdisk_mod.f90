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
  use utilities_mod
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
  function ram_fullWorkingPath(fileName, noAbort) result(fullWorkingPath)
    implicit none
    character(len=512) :: fullWorkingPath
    logical, optional  :: noAbort
    character(len=*)   :: fileName

    logical            :: fileExists, noAbort2
    character(len=256) :: fileName2

    if ( .not. initialized ) then
          call utl_abort('ram_fullWorkingPath: ramDisk module has not been initialized.')
    end if

    if ( present(noAbort) ) then
      noAbort2 = noAbort
    else
      noAbort2 = .false.
    end if

    ! this should make it safe for calls where input and output are the same variable
    fileName2 = trim(fileName)

    ! first look for file in the ram disk directory
    fullWorkingPath = trim(ram_disk_dir) // '/' // trim(fileName2)
    inquire(file=trim(fullWorkingPath),exist=fileExists)

    if ( fileExists ) then
      write(*,*) 'ram_fullWorkingPath: this file found on ram disk: ', trim(fileName2)
    else
      ! now look in working directory
      fullWorkingPath = './' // trim(fileName2)
      inquire(file=trim(fullWorkingPath),exist=fileExists)

      if ( .not. fileExists ) then
        if ( noAbort2 ) then
          fullWorkingPath = ' '
        else
          write(*,*) 'ram_fullWorkingPath: file name          = ', trim(fileName2)
          write(*,*) 'ram_fullWorkingPath: ram disk directory = ', trim(ram_disk_dir)
          call utl_abort('ram_fullWorkingPath: this file cannot be found.')
        end if
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

    ! The following interface was extracted from #include <clib_interface.cdk>
    interface clib_remove
      integer function clib_remove_schhide(path)
        implicit none
        character(len=*),intent(IN) :: path
      end function
    end interface

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

end module ramDisk_mod