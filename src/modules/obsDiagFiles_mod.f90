
module obsDiagFiles_mod
  ! MODULE obsDiagFiles_mod (prefix='diaf' category='3. Observation input/output')
  !
  !:Purpose:  To write the "diag" format SQLITE observation files. Data is stored in
  !           obsSpaceData object.
  !
  use midasMpi_mod
  use codePrecision_mod
  use mathPhysConstants_mod
  use obsSpaceData_mod
  use utilities_mod
  use ramDisk_mod
  use codtyp_mod
  use ensembleObservations_mod

  implicit none
  save
  private

  ! Public procedures
  public :: diaf_writeAllSqlDiagFiles

  contains

  !--------------------------------------------------------------------------
  ! diaf_writeAllSqlDiagFiles
  !--------------------------------------------------------------------------
  subroutine diaf_writeAllSqlDiagFiles(obsdat, sfFileName, onlyAssimObs, addFSOdiag, ensObs_opt)
    !
    ! :Purpose: To prepare the writing of obsSpaceData content into SQLite format files
    !
    implicit none

    ! Arguments:
    type(struct_obs),           intent(inout) :: obsdat         ! obsSpaceData object
    character(len=*),           intent(in)    :: sfFileName     ! fileName acronym used for surface obs file
    logical         ,           intent(in)    :: onlyAssimObs   ! only write assimilated obs
    logical         ,           intent(in)    :: addFSOdiag     ! include FSO column in body table
    type(struct_eob), optional, intent(in)    :: ensObs_opt     ! ensObs object

    ! Locals:
    integer                :: familyIndex
    character(len=2)       :: obsFamilyList(50)
    integer                :: obsFamilyListSize
    character(len=codtyp_name_length) :: fileName

    ! ensure all mpi tasks have same list of common obs family names
    call diaf_getObsFamilyListMpiGlobal(obsdat, obsFamilyListSize, obsFamilyList)

    ! get list of all possible tovs codetype values and unique list of corresponding filenames

    do familyIndex = 1, obsFamilyListSize

      write(*,*) 'diaf_writeAllSqlDiagFiles: Family = ', familyIndex, obsFamilyList(familyIndex)

      fileName = diaf_getObsFileName(obsFamilyList(familyIndex), sfFileName_opt=sfFileName)
      call diaf_writeSqlDiagFile(obsdat, obsFamilyList(familyIndex), &
                                 onlyAssimObs, addFSOdiag, fileName, &
                                 ensObs_opt=ensObs_opt )

    end do

  end subroutine diaf_writeAllSqlDiagFiles

  !--------------------------------------------------------------------------
  ! diaf_writeSqlDiagFile
  !--------------------------------------------------------------------------
  subroutine diaf_writeSqlDiagFile(obsdat, obsFamily, onlyAssimObs, addFSOdiag, instrumentFileName, codeTypeList_opt, ensObs_opt)
    !
    ! :Purpose: To write the obsSpaceData content into SQLite format files
    !
    implicit none

    ! Arguments:
    type(struct_obs)           , intent(inout) :: obsdat
    character(len=*)           , intent(in)    :: obsFamily
    logical                    , intent(in)    :: onlyAssimObs
    logical                    , intent(in)    :: addFSOdiag
    character(len=*)           , intent(in)    :: instrumentFileName
    integer          , optional, intent(in)    :: codeTypeList_opt(:)
    type(struct_eob) , optional, intent(in)    :: ensObs_opt

    call utl_abort('diaf_writeSqlDiagFile: The code related to SQL has been removed!')
  end subroutine diaf_writeSqlDiagFile

  !--------------------------------------------------------------------------
  ! diaf_getObsFileName
  !--------------------------------------------------------------------------
  function diaf_getObsFileName(obsFamily, sfFileName_opt, codetype_opt) result(fileName)
    !
    ! :Purpose: Return the part of the observation file name associated
    !           with the type of observation it contains.
    !
    implicit none

    ! Arguments:
    character(len=*)          , intent(in) :: obsFamily
    character(len=*), optional, intent(in) :: sfFileName_opt ! fileName acronym used for surface obs file
    integer,          optional, intent(in) :: codetype_opt
    ! Result:
    character(len=codtyp_name_length)      :: fileName

    if (obsFamily == 'TO') then
      if (.not. present(codetype_opt)) then
        call utl_abort('diaf_getObsFileName: codetype_opt must be specified for TO family')
      end if

      if (codtyp_get_name(codeType_opt) == 'radianceclear') then
        fileName  = 'csr'
      else if (codtyp_get_name(codeType_opt) == 'mhs' .or. codtyp_get_name(codeType_opt) == 'amsub') then
        fileName = 'to_amsub'
      else if (codtyp_get_name(codeType_opt) == 'amsua') then
        fileName = 'to_amsua'
      else if (codtyp_get_name(codeType_opt) == 'ssmi') then
        fileName = 'ssmis'
      else if (codtyp_get_name(codeType_opt) == 'crisfsr') then
        fileName = 'cris'
      else if (codtyp_get_name(codeType_opt) == 'atms') then
        fileName = 'atms'
      else
        fileName = codtyp_get_name(codeType_opt)
      end if
    else
      if (.not. present(sfFileName_opt)) then
        call utl_abort('diaf_getObsFileName: sfFileName_opt must be specified')
      end if
      call up2low(obsFamily, fileName)
      if (fileName == 'ra') fileName = 'radar'
      if (fileName == 'sf') then
        ! use either 'sf' or 'sfc' for filename with surface obs
        fileName = sfFileName_opt
      end if
    end if

  end function diaf_getObsFileName

  !--------------------------------------------------------------------------
  ! diaf_getObsFamilyListMpiGlobal
  !--------------------------------------------------------------------------
  subroutine diaf_getObsFamilyListMpiGlobal(obsdat, obsFamilyListSizeCommon,  &
                                            obsFamilyListCommon)
    !
    ! :Purpose: Obtain a common set of obs family names over all mpi tasks
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsdat
    integer         , intent(out)   :: obsFamilyListSizeCommon
    character(len=*), intent(out)   :: obsFamilyListCommon(:)

    ! Locals:
    integer                       :: headerIndex, familyIndex, charIndex, procIndex
    integer                       :: obsFamilyListSizeMpiLocal, obsFamilyListSizeMaxMpiLocal, obsFamilyListSizeMax
    character(len=2), allocatable :: obsFamilyListMpiLocal(:)
    character(len=2), allocatable :: obsFamilyListMpiGlobal(:,:)
    character(len=2)              :: currentObsFamily
    integer, allocatable          :: intObsFamilyListMpiLocal(:,:)
    integer, allocatable          :: intObsFamilyListMpiGlobal(:,:,:)
    integer, allocatable          :: allObsFamilyListSizeMpiLocal(:)

    obsFamilyListSizeMax = size(obsFamilyListCommon)
    write(*,*) 'obsFamilyListSizeMax =', obsFamilyListSizeMax

    ! get family list for this mpi task
    obsFamilyListSizeMpiLocal = 0
    allocate(obsFamilyListMpiLocal(obsFamilyListSizeMax))
    obsFamilyListMpiLocal(:) = 'XX'
    HEADER: do headerIndex = 1, obs_numHeader(obsdat)
      currentObsFamily = obs_getFamily(obsdat, headerIndex_opt=headerIndex)
      if (any(obsFamilyListMpiLocal(:) == currentObsFamily)) cycle HEADER
      obsFamilyListSizeMpiLocal = obsFamilyListSizeMpiLocal + 1
      obsFamilyListMpiLocal(obsFamilyListSizeMpiLocal) = currentObsFamily
      write(*,*) 'add the family: ', currentObsFamily
    end do HEADER
    write(*,*) 'obsFamilyListSizeMpiLocal =', obsFamilyListSizeMpiLocal
    write(*,*) 'obsFamilyListMpiLocal = ', obsFamilyListMpiLocal(1:obsFamilyListSizeMpiLocal)

    allocate(allObsFamilyListSizeMpiLocal(mmpi_nprocs))
    call mmpi_allGather(obsFamilyListSizeMpiLocal, allObsFamilyListSizeMpiLocal)
    call mmpi_allReduce(obsFamilyListSizeMpiLocal, obsFamilyListSizeMaxMpiLocal, mmpi_max)

    ! convert local family list from characters to integers
    allocate(intObsFamilyListMpiLocal(len(currentObsFamily),obsFamilyListSizeMaxMpiLocal))
    intObsFamilyListMpiLocal(:,:)=0
    do familyIndex = 1, obsFamilyListSizeMpiLocal
      do charIndex = 1, len(currentObsFamily)
        intObsFamilyListMpiLocal(charIndex,familyIndex) =  &
             iachar(obsFamilyListMpiLocal(familyIndex)(charIndex:charIndex))
      end do
    end do

    ! communicate obs family list to all mpi tasks as integers
    allocate(intObsFamilyListMpiGlobal(len(currentObsFamily),obsFamilyListSizeMaxMpiLocal,mmpi_nprocs))
    call mmpi_allGather(intObsFamilyListMpiLocal, intObsFamilyListMpiGlobal)

    ! convert global family lists from integers to characters
    allocate(obsFamilyListMpiGlobal(obsFamilyListSizeMaxMpiLocal,mmpi_nprocs))
    obsFamilyListMpiGlobal(:,:) = 'XX'
    do procIndex = 1, mmpi_nprocs
      do familyIndex = 1, allObsFamilyListSizeMpiLocal(procIndex)
        do charIndex=1,len(currentObsFamily)
          obsFamilyListMpiGlobal(familyIndex,procIndex)(charIndex:charIndex) =  &
               achar(intObsFamilyListMpiGlobal(charIndex,familyIndex,procIndex))
        end do
      end do
      write(*,*) 'obsFamilyListMpiGlobal = ', procIndex,  &
           obsFamilyListMpiGlobal(1:allObsFamilyListSizeMpiLocal(procIndex),procIndex)
    end do

    ! construct single common list of families to be used for all mpi tasks
    obsFamilyListCommon(:) = 'YY'
    obsFamilyListSizeCommon = 0
    do procIndex = 1, mmpi_nprocs
      FAMILY: do familyIndex = 1, obsFamilyListSizeMaxMpiLocal
        if (obsFamilyListMpiGlobal(familyIndex,procIndex) == 'XX') cycle FAMILY
        if (any(obsFamilyListCommon(:) == obsFamilyListMpiGlobal(familyIndex,procIndex))) cycle FAMILY
        obsFamilyListSizeCommon = obsFamilyListSizeCommon + 1
        obsFamilyListCommon(obsFamilyListSizeCommon) = obsFamilyListMpiGlobal(familyIndex,procIndex)
      end do FAMILY
    end do
    write(*,*) 'obsFamilyListSizeCommon = ', obsFamilyListSizeCommon
    write(*,*) 'obsFamilyListCommon = ', obsFamilyListCommon(1:obsFamilyListSizeCommon)

    deallocate(allObsFamilyListSizeMpiLocal)
    deallocate(obsFamilyListMpiGlobal)
    deallocate(intObsFamilyListMpiGlobal)
    deallocate(intObsFamilyListMpiLocal)
    deallocate(obsFamilyListMpiLocal)

  end subroutine diaf_getObsFamilyListMpiGlobal

end module obsDiagFiles_mod
