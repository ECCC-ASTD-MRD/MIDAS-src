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
!! MODULE obsFiles_mod (prefix="obsf")
!!
!! *Purpose*: High-level module to handle reading/writing of observations that
!!            can be stored in one of several different formats. Currently, the
!!            only supported formats are:
!!
!!            1. BURP
!!            2. CMA (binary format of obsSpaceData contents)
!!
!--------------------------------------------------------------------------
module obsFiles_mod
  use mpi_mod
  use utilities_mod
  use obsSpaceData_mod
  use burpFiles_mod
  use cmaFiles_mod
  implicit none
  save
  private

  ! public procedures
  public :: obsf_setup, obsf_filesSplit, obsf_getFileType, obsf_fileTypeIsBurp
  public :: obsf_readFiles, obsf_writeFiles

  character(len=10) :: obsFileType
  logical           :: obsFilesSplit
  logical           :: initialized = .false.

contains

  subroutine obsf_setup(dateStamp_out,obsFileMode)
    implicit none

    ! arguments
    integer :: dateStamp_out
    character(len=*) :: obsFileMode

    ! locals
    integer :: status, length_obsFileType

    status = 0
    call get_environment_variable('MIDAS_OBS_FILE_TYPE',obsFileType,length_obsFileType,status,.true.)

    if (status.gt.1) then
      write(*,*) 'obsf_setup: Problem when getting the environment variable MIDAS_OBS_FILE_TYPE'
    end if
    if (status.eq.1) then
      write(*,*) 'obsf_setup: The environment variable MIDAS_OBS_FILE_TYPE has not been detected!'
      write(*,*) '            Assume the obs file type is: BURP'
      obsFileType = 'BURP'
    else
      write(*,*)
      write(*,*) 'obsf_setup: The environment variable MIDAS_OBS_FILE_TYPE has correctly been detected'
      write(*,*) 'obsf_setup: The obs file type is : ', trim(obsFileType)
    end if

    !
    ! Determine if obsFiles are split
    !
    if ( obsFileType == 'BURP' ) then
      obsFilesSplit = .true.
    else if ( obsFileType == 'CMA' ) then
      obsFilesSplit = .false.
    else
      call utl_abort('obsf_setup: invalid observation file type: ' // trim(obsFileType))
    end if

    !
    ! Do some setup of observation files
    !
    if ( obsFileType == 'BURP' ) then
      call burp_setupFiles(dateStamp_out,obsFileMode)
    else
      dateStamp_out = -1
    end if

    initialized = .true.

  end subroutine obsf_setup


  function obsf_filesSplit() result(obsFilesSplit_out)
    implicit none

    ! arguments
    logical :: obsFilesSplit_out

    if ( .not.initialized ) call utl_abort('obsf_filesSplit: obsFiles_mod not initialized!')

    obsFilesSplit_out = obsFilesSplit
        
  end function obsf_filesSplit


  subroutine obsf_getFileType(obsFileType_out)
    implicit none

    ! arguments
    character(len=*) :: obsFileType_out

    if ( .not.initialized ) call utl_abort('obsf_getFileType: obsFiles_mod not initialized!')

    obsFileType_out = obsFileType
        
  end subroutine obsf_getFileType


  function obsf_fileTypeIsBurp() result(fileTypeIsBurp)
    implicit none
 
    ! arguments
    logical :: fileTypeIsBurp

    fileTypeIsBurp = ( trim(obsFileType) == 'BURP' )

  end function obsf_fileTypeIsBurp


  subroutine obsf_readFiles(obsSpaceData)
    implicit none

    ! arguments
    type(struct_obs) :: obsSpaceData

    if ( .not.initialized ) call utl_abort('obsf_readFiles: obsFiles_mod not initialized!')

    if ( obsFileType == 'BURP' ) then
      call burp_readFiles(obsSpaceData)
    else if ( obsFileType == 'CMA' ) then
      ! read same global CMA file on all mpi tasks
      call cma_readFiles(obsSpaceData)
    end if

  end subroutine obsf_readFiles


  subroutine obsf_writeFiles(obsSpaceData,HXensT_mpiglobal_opt,asciDumpObs_opt)
    implicit none

    ! arguments
    type(struct_obs)           :: obsSpaceData
    real(8), pointer, optional :: HXensT_mpiglobal_opt(:,:)
    logical, optional          :: asciDumpObs_opt

    if ( .not.initialized ) call utl_abort('obsf_writeFiles: obsFiles_mod not initialized!')

    if ( (present(HXensT_mpiglobal_opt) .or. present(asciDumpObs_opt)) .and. obsFileType .ne. 'CMA' ) then
      call utl_abort('obsf_writeFiles: writing of ensemble Hx and asciDump ' //  &
                     'are only compatible with CMA obs file type')
    end if

    if ( obsFileType == 'BURP' ) then
      call burp_updateFiles(obsSpaceData)
    else if ( obsFileType == 'CMA' ) then
      ! only 1 mpi task should do the writing
      call cma_writeFiles(obsSpaceData,HXensT_mpiglobal_opt,asciDumpObs_opt)
    end if

  end subroutine obsf_writeFiles

end module obsFiles_mod
