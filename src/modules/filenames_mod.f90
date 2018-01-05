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
!! MODULE fileNames_mod (prefix="fln")
!!
!! *Purpose*: Routines related to file names
!!
!--------------------------------------------------------------------------
module fileNames_mod
  use utilities_mod
  use ramDisk_mod
  implicit none
  save
  private

  ! public procedures
  public :: fln_ensFileName

contains

 !--------------------------------------------------------------------------
 ! fln_ensFileName
 !--------------------------------------------------------------------------
  subroutine fln_ensFileName(ensFileName, ensPathName, memberIndex, ensFileNamePrefix_opt, ensFileBaseName_opt, shouldExist_opt)
    implicit none

    ! arguments
    character(len=*)  :: ensFileName
    character(len=*)  :: ensPathName
    integer           :: memberIndex
    character(len=*),optional  :: ensFileBaseName_opt, ensFileNamePrefix_opt
    logical, optional :: shouldExist_opt

    ! locals
    integer          :: numFiles, returnCode, totalLength, ensembleBaseFileNameLength
    character(len=4) :: ensNumber
    logical          :: shouldExist
    character(len=2000) :: fileList(10), fileNamePattern
    character        :: ensembleFileExtLengthStr
    logical, save    :: firstTime = .true.
    integer, save    :: ensembleFileExtLength = 4
    character(len=200), save :: ensFileBaseName

    ! The following interface was extracted from #include <clib_interface.cdk>
    interface clib_glob
      integer function clib_glob_schhide(filelist,nfiles,pattern,maxnfiles)
        implicit none
        integer,intent(IN)  :: maxnfiles
        character(len=*),intent(IN) :: pattern
        integer,intent(OUT) :: nfiles
        character(len=*),dimension(maxnfiles),intent(OUT):: filelist
      end function clib_glob_schhide
    end interface clib_glob

    if ( present(shouldExist_opt) ) then
      shouldExist = shouldExist_opt
    else
      shouldExist = .true.
    end if

    if (firstTime) then
      write(*,*) 'ens_fileName: trying to get ensemble file name'
      fileNamePattern = './' // trim(enspathname) // '/' // '*[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9]_*001'
      returnCode = clib_glob(fileList,numFiles,trim(fileNamePattern),10)
      if (returnCode /= 1) then
        call utl_abort('ens_fileName: reached maximum number of files or no files are available')
      end if
      write(*,*) 'ens_fileName: fileList = ',trim(fileList(1))
      write(*,*) 'ens_fileName: numFiles, returnCode = ', numFiles,returnCode

      ensFileName = trim(fileList(1))
      totalLength = len_trim(ensFileName)

      ! find number of digits used to identify ensemble member index
      ensembleFileExtLength = 0
      do
        if ( ensFileName((totalLength-ensembleFileExtLength):(totalLength-ensembleFileExtLength)) == '_' ) exit
        ensembleFileExtLength = ensembleFileExtLength + 1
      end do

      ! find number of digits used to identify the basename of the file
      ensembleBaseFileNameLength = 0
      do
        if ( ensFileName((totalLength-ensembleBaseFileNameLength):(totalLength-ensembleBaseFileNameLength)) == '/' ) exit
        ensembleBaseFileNameLength = ensembleBaseFileNameLength + 1
      end do

      if (ensembleFileExtLength == 0) then
        call utl_abort('ens_fileName: Cannot determine the ensemble file extention length with ' // trim(ensFileName))
      end if

      ensFileBasename = ensFileName((totalLength-ensembleBaseFileNameLength+1):(totalLength-ensembleFileExtLength-1))
      
      firstTime = .false.
    end if

    write(ensembleFileExtLengthStr,'(i1.1)') ensembleFileExtLength
    write(ensNumber,'(i' // ensembleFileExtLengthStr // '.' // ensembleFileExtLengthStr // ')') memberIndex
    ensFileName = trim(enspathname) // '/' // trim(ensFileBaseName) // '_' // trim(ensNumber)

    if ( shouldExist ) ensFileName = ram_fullWorkingPath(ensFileName)

    if (present(ensFileBaseName_opt)) ensFileBaseName_opt = trim(ensFileBaseName)

  end subroutine fln_ensFileName

end module fileNames_mod
