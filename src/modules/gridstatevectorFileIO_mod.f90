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

module gridStateVectorFileIO_mod
  ! MODULE gridStateVectorFile_mod (prefix='gio' category='1. High-level functionality')
  !
  ! :Purpose: The grid-point state vector I/O methods.
  !
  use mpi
  use mpi_mod
  use gridStateVector_mod
  use utilities_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use oceanMask_mod
  use varNameList_mod
  use ramDisk_mod
  use timeCoord_mod
  use mathPhysConstants_mod
  use codePrecision_mod

  ! public subroutines and functions
  public :: gio_readFromFile, gio_readTrials, gio_readFile
  public :: gio_readMaskFromFile
  public :: gio_fileUnitsToStateUnits
  public :: gio_getMaskLAM

  integer, external :: get_max_rss

  contains
  !--------------------------------------------------------------------------
  ! gio_readFromFile
  !--------------------------------------------------------------------------
  subroutine gio_readFromFile(statevector_out, fileName, etiket_in, typvar_in, stepIndex_opt,  &
                              unitConversion_opt, PsfcReference_opt, readHeightSfc_opt,   &
                              containsFullField_opt, vcoFileIn_opt)
    implicit none

    ! arguments
    type(struct_gsv)              :: statevector_out
    character(len=*), intent(in)  :: fileName
    character(len=*), intent(in)  :: etiket_in
    character(len=*), intent(in)  :: typvar_in
    integer, optional             :: stepIndex_opt
    logical, optional             :: unitConversion_opt
    logical, optional             :: readHeightSfc_opt
    logical, optional,intent(in)  :: containsFullField_opt
    real(8), optional             :: PsfcReference_opt(:,:)
    type(struct_vco), optional, pointer, intent(in)  :: vcoFileIn_opt

    ! locals
    integer :: stepIndex, varIndex
    character(len=4) :: varName
    logical :: doHorizInterp, doVertInterp, unitConversion
    logical :: readHeightSfc, containsFullField
    type(struct_vco), pointer :: vco_file
    type(struct_hco), pointer :: hco_file
    logical :: foundVarNameInFile 

    nullify(vco_file, hco_file)

    write(*,*) ''
    write(*,*) 'gio_readFromFile: START'
    call tmg_start(158,'gio_readFromFile')

    if ( present(stepIndex_opt) ) then
      stepIndex = stepIndex_opt
    else
      stepIndex = statevector_out%anltime
    end if
    if ( stepIndex > stateVector_out%numStep .or. stepIndex < 1 ) then
      write(*,*) 'stepIndex = ', stepIndex
      call utl_abort('gio_readFromFile: invalid value for stepIndex')
    end if

    if ( present(unitConversion_opt) ) then
      unitConversion = unitConversion_opt
    else
      unitConversion = .true.
    end if

    if (present(containsFullField_opt)) then
      containsFullField = containsFullField_opt
    else
      containsFullField = .true.
    end if
    write(*,*) 'gio_readFromFile: containsFullField = ', containsFullField

    if ( present(readHeightSfc_opt) ) then
      readHeightSfc = readHeightSfc_opt
    else
      readHeightSfc = .false.
    end if

    ! set up vertical and horizontal coordinate for input file
    if ( present(vcoFileIn_opt)) then
      vco_file => vcoFileIn_opt
      write(*,*)
      write(*,*) 'gio_readFromFile: vertical levels defined in user-supplied vco object will be read'
    else
      call vco_setupFromFile(vco_file,trim(fileName),beSilent_opt=.true.)
      write(*,*)
      write(*,*) 'gio_readFromFile: all the vertical levels will be read from ', trim(fileName) 
    end if

    foundVarNameInFile = .false.

    do varIndex = 1, vnl_numvarmax
      varName = vnl_varNameList(varIndex)

      if ( .not. gsv_varExist(statevector_out,varName) ) cycle

      ! make sure variable is in the file
      if ( .not. utl_varNamePresentInFile(varName,fileName_opt=trim(fileName)) ) cycle

      ! adopt a variable on the full/dynamic LAM grid
      if ( .not. statevector_out%hco%global .and. (trim(varName) == 'TM' .or. trim(varName) == 'MG')) cycle

      foundVarNameInFile = .true.

      exit

    end do

    ! special case when only TM (Surface Temperature) is in the file:
    if ( .not. foundVarNameInFile ) then
      varname = 'TM'
      if ( gsv_varExist( statevector_out, varname ) .and. &
           utl_varNamePresentInFile( varname, fileName_opt = trim( fileName ))) &
        foundVarNameInFile = .true.
    end if   

    ! to be safe for situations where, e.g. someone wants to only read MG from a file
    if ( .not. foundVarNameInFile ) then
      varname = 'P0'
      if ( utl_varNamePresentInFile( varname, fileName_opt = trim( fileName ))) &
        foundVarNameInFile = .true.
    end if   

    if ( .not. foundVarNameInFile) call utl_abort('gio_readFromFile: NO variables found in the file!!!')

    write(*,*) 'gio_readFromFile: defining hco by varname= ', varName

    call hco_setupFromFile( hco_file, trim(fileName), ' ', gridName_opt='FILEGRID', varName_opt = varName )

    ! test if horizontal and/or vertical interpolation needed for statevector grid
    doVertInterp = .not.vco_equal(vco_file,statevector_out%vco)
    doHorizInterp = .not.hco_equal(hco_file,statevector_out%hco)
    write(*,*) 'gio_readFromFile: doVertInterp = ', doVertInterp, ', doHorizInterp = ', doHorizInterp

    ! call appropriate subroutine to do actual work
    if ( (doVertInterp .or. doHorizInterp) .and. statevector_out%mpi_distribution=='Tiles' ) then
      call readFromFileAndInterpToTiles(statevector_out, fileName,  &
             vco_file, hco_file, etiket_in, typvar_in, stepIndex, unitConversion,  &
             readHeightSfc, containsFullField, PsfcReference_opt)
    else if ( (doVertInterp .or. doHorizInterp) .and. .not.stateVector_out%mpi_local ) then
      call readFromFileAndInterp1Proc(statevector_out, fileName,  &
             vco_file, hco_file, etiket_in, typvar_in, stepIndex, unitConversion,  &
             readHeightSfc, containsFullField)
    else if ( .not.(doVertInterp .or. doHorizInterp) .and. stateVector_out%mpi_local ) then
      call readFromFileAndTransposeToTiles(statevector_out, fileName,  &
             etiket_in, typvar_in, stepIndex, unitConversion,  &
             readHeightSfc, containsFullField)
    else
      call readFromFileOnly(statevector_out, fileName,  &
                                etiket_in, typvar_in, stepIndex, unitConversion,  &
                                readHeightSfc, containsFullField)
    end if

    call tmg_stop(158)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'gio_readFromFile: END'

  end subroutine gio_readFromFile

  !--------------------------------------------------------------------------
  ! readFromFileAndInterpToTiles
  !--------------------------------------------------------------------------
  subroutine readFromFileAndInterpToTiles(statevector_out, fileName,  &
             vco_file, hco_file, etiket_in, typvar_in, stepIndex, unitConversion,  &
             readHeightSfc, containsFullField, PsfcReference_opt)
    implicit none
    ! Note this routine currently only works correctly for reading FULL FIELDS,
    ! not increments or perturbations... because of the HU -> LQ conversion

    ! arguments
    type(struct_gsv)              :: statevector_out
    character(len=*), intent(in)  :: fileName
    type(struct_vco), pointer     :: vco_file
    type(struct_hco), pointer     :: hco_file
    character(len=*), intent(in)  :: etiket_in
    character(len=*), intent(in)  :: typvar_in
    integer                       :: stepIndex
    logical                       :: unitConversion
    logical                       :: readHeightSfc
    logical                       :: containsFullField
    real(8), optional             :: PsfcReference_opt(:,:)

    ! locals
    type(struct_gsv) :: statevector_file_r4, statevector_tiles, statevector_hinterp_r4, statevector_vinterp

    real(4), pointer     :: field3d_r4_ptr(:,:,:)
    real(8), pointer     :: field_in_ptr(:,:,:,:), field_out_ptr(:,:,:,:)
    real(8), allocatable :: PsfcReference3D(:,:,:)

    character(len=4), pointer :: varNamesToRead(:)

    nullify(field3d_r4_ptr, field_in_ptr, field_out_ptr)

    write(*,*) ''
    write(*,*) 'readFromFileAndInterpToTiles: START'

    nullify(varNamesToRead)
    call gsv_varNamesList(varNamesToRead, statevector_out)

    !-- 1.0 Read the file, distributed over mpi task with respect to variables/levels

    ! initialize single precision 3D working copy of statevector for reading file
    call gsv_allocate(statevector_file_r4, 1, hco_file, vco_file,                &
                      dateStamp_opt=statevector_out%datestamplist(stepIndex),    &
                      mpi_local_opt=.true., mpi_distribution_opt='VarsLevs',     &
                      dataKind_opt=4, allocHeightSfc_opt=readHeightSfc,          &
                      varNames_opt=varNamesToRead,                               &
                      hInterpolateDegree_opt=statevector_out%hInterpolateDegree, &
                      hExtrapolateDegree_opt=statevector_out%hExtrapolateDegree )

    call gio_readFile(statevector_file_r4, filename, etiket_in, typvar_in,  &
                      containsFullField, readHeightSfc_opt=readHeightSfc)

    !-- 2.0 Horizontal Interpolation

    ! initialize single precision 3D working copy of statevector for horizontal interpolation result
    call gsv_allocate(statevector_hinterp_r4, 1, statevector_out%hco, vco_file,  &
                      dateStamp_opt=statevector_out%datestamplist(stepIndex),    &
                      mpi_local_opt=.true., mpi_distribution_opt='VarsLevs',     &
                      dataKind_opt=4, allocHeightSfc_opt=readHeightSfc,          &
                      varNames_opt=varNamesToRead,                               &
                      hInterpolateDegree_opt=statevector_out%hInterpolateDegree, &
                      hExtrapolateDegree_opt=statevector_out%hExtrapolateDegree )

    write(*,*) 'DBGmad calling gsv_hInterpolate_r4'
    call gsv_hInterpolate_r4(statevector_file_r4, statevector_hinterp_r4)

    call gsv_deallocate(statevector_file_r4)

    !-- 3.0 Unit conversion

    if ( unitConversion ) then
      write(*,*) 'DBGmad calling gsv_fileUnitsToStateUnits'
      call gio_fileUnitsToStateUnits( statevector_hinterp_r4, containsFullField )
    end if

    !-- 4.0 MPI communication from vars/levels to lat/lon tiles

    ! initialize double precision 3D working copy of statevector for mpi communication result
    call gsv_allocate(statevector_tiles, 1, statevector_out%hco, vco_file,    &
                      dateStamp_opt=statevector_out%datestamplist(stepIndex), &
                      mpi_local_opt=.true., mpi_distribution_opt='Tiles',     &
                      dataKind_opt=8, allocHeightSfc_opt=readHeightSfc,               &
                      varNames_opt=varNamesToRead )

    write(*,*) 'DBGmad calling gsv_transposeVarsLevsToTiles'
    call gsv_transposeVarsLevsToTiles(statevector_hinterp_r4, statevector_tiles)

    call gsv_deallocate(statevector_hinterp_r4)

    !-- 5.0 Vertical interpolation

    ! initialize double precision 3D working copy of statevector for mpi communication result
    call gsv_allocate(statevector_vinterp, 1, statevector_out%hco, statevector_out%vco, &
                      dateStamp_opt=statevector_out%datestamplist(stepIndex), &
                      mpi_local_opt=.true., mpi_distribution_opt='Tiles', dataKind_opt=8, &
                      allocHeightSfc_opt=readHeightSfc, varNames_opt=varNamesToRead )

    if (present(PsfcReference_opt) ) then
      write(*,*) 'DBGmad PsfcReference_opt present, caling gsv_vInterpolate'
      allocate(PsfcReference3D(statevector_tiles%myLonBeg:statevector_tiles%myLonEnd, &
                               statevector_tiles%myLatBeg:statevector_tiles%myLatEnd,1))
      PsfcReference3D(:,:,1) = PsfcReference_opt(:,:)
      call gsv_vInterpolate(statevector_tiles,statevector_vinterp,PsfcReference_opt=PsfcReference3D)
      deallocate(PsfcReference3D)
    else
      write(*,*) 'DBGmad PsfcReference_opt not present, caling gsv_vInterpolate'
      call gsv_vInterpolate(statevector_tiles,statevector_vinterp)
    end if

    call gsv_deallocate(statevector_tiles)

    !-- 6.0 Copy result to output statevector
    call gsv_copy(statevector_vinterp, statevector_out, stepIndexOut_opt=stepIndex)

    call gsv_deallocate(statevector_vinterp)
    deallocate(varNamesToRead)

    write(*,*) 'readFromFileAndInterpToTiles: END'

  end subroutine readFromFileAndInterpToTiles

  !--------------------------------------------------------------------------
  ! readFromFileAndTransposeToTiles
  !--------------------------------------------------------------------------
  subroutine readFromFileAndTransposeToTiles(statevector_out, fileName,  &
             etiket_in, typvar_in, stepIndex, unitConversion,  &
             readHeightSfc, containsFullField)
    implicit none
    ! Note this routine currently only works correctly for reading FULL FIELDS,
    ! not increments or perturbations... because of the HU -> LQ conversion

    ! arguments
    type(struct_gsv)              :: statevector_out
    character(len=*), intent(in)  :: fileName
    character(len=*), intent(in)  :: etiket_in
    character(len=*), intent(in)  :: typvar_in
    integer                       :: stepIndex
    logical                       :: unitConversion
    logical                       :: readHeightSfc
    logical                       :: containsFullField

    ! locals
    type(struct_gsv) :: statevector_file_r4, statevector_tiles

    real(4), pointer     :: field3d_r4_ptr(:,:,:)
    real(8), pointer     :: field_in_ptr(:,:,:,:), field_out_ptr(:,:,:,:)

    character(len=4), pointer :: varNamesToRead(:)

    nullify(field3d_r4_ptr, field_in_ptr, field_out_ptr)

    write(*,*) ''
    write(*,*) 'readFromFileAndTransposeToTiles: START'

    nullify(varNamesToRead)
    call gsv_varNamesList(varNamesToRead, statevector_out)

    !-- 1.0 Read the file, distributed over mpi task with respect to variables/levels

    ! initialize single precision 3D working copy of statevector for reading file
    call gsv_allocate(statevector_file_r4, 1, statevector_out%hco, statevector_out%vco, &
                      dateStamp_opt=statevector_out%datestamplist(stepIndex),           &
                      mpi_local_opt=.true., mpi_distribution_opt='VarsLevs',            &
                      dataKind_opt=4, allocHeightSfc_opt=readHeightSfc,                         &
                      varNames_opt=varNamesToRead )

    call gio_readFile(statevector_file_r4, filename, etiket_in, typvar_in,  &
                      containsFullField, readHeightSfc_opt=readHeightSfc)

    !-- 2.0 Unit conversion
    if ( unitConversion ) then
      call gio_fileUnitsToStateUnits( statevector_file_r4, containsFullField )
    end if

    !-- 3.0 MPI communication from vars/levels to lat/lon tiles
    call gsv_allocate(statevector_tiles, 1, statevector_out%hco, statevector_out%vco, &
                      dateStamp_opt=statevector_out%datestamplist(stepIndex), &
                      mpi_local_opt=.true., mpi_distribution_opt='Tiles',     &
                      dataKind_opt=8, allocHeightSfc_opt=readHeightSfc,               &
                      varNames_opt=varNamesToRead)

    call gsv_transposeVarsLevsToTiles(statevector_file_r4, statevector_tiles)

    !-- 4.0 Copy result to output statevector
    call gsv_copy(statevector_tiles, statevector_out, stepIndexOut_opt=stepIndex)

    call gsv_deallocate(statevector_file_r4)
    deallocate(varNamesToRead)

    write(*,*) 'readFromFileAndTransposeToTiles: END'

  end subroutine readFromFileAndTransposeToTiles

  !--------------------------------------------------------------------------
  ! readFromFileAndInterp1Proc
  !--------------------------------------------------------------------------
  subroutine readFromFileAndInterp1Proc(statevector_out_r4, fileName,  &
             vco_file, hco_file, etiket_in, typvar_in, stepIndex, unitConversion,  &
             readHeightSfc, containsFullField)
    implicit none

    ! arguments
    type(struct_gsv)              :: statevector_out_r4
    character(len=*), intent(in)  :: fileName
    type(struct_vco), pointer     :: vco_file
    type(struct_hco), pointer     :: hco_file
    character(len=*), intent(in)  :: etiket_in
    character(len=*), intent(in)  :: typvar_in
    integer                       :: stepIndex
    logical                       :: unitConversion
    logical                       :: readHeightSfc
    logical                       :: containsFullField

    ! locals
    type(struct_gsv) :: statevector_file_r4, statevector_hinterp_r4, statevector_vinterp_r4

    real(4), pointer     :: field_in_ptr(:,:,:,:), field_out_ptr(:,:,:,:)

    character(len=4), pointer :: varNamesToRead(:)

    nullify(field_in_ptr, field_out_ptr)

    write(*,*) ''
    write(*,*) 'readFromFileAndInterp1Proc: START'

    nullify(varNamesToRead)
    call gsv_varNamesList(varNamesToRead, statevector_out_r4)

    !-- 1.0 Read the file

    ! initialize single precision 3D working copy of statevector for reading file
    call gsv_allocate(statevector_file_r4, 1, hco_file, vco_file,                    &
                      dateStamp_opt=statevector_out_r4%datestamplist(stepIndex),     &
                      mpi_local_opt=.false., dataKind_opt=4,                         &
                      allocHeightSfc_opt=readHeightSfc, varNames_opt=varNamesToRead, &
                      hInterpolateDegree_opt=statevector_out_r4%hInterpolateDegree,  &
                      hExtrapolateDegree_opt=statevector_out_r4%hExtrapolateDegree)

    call gio_readFile(statevector_file_r4, filename, etiket_in, typvar_in,  &
                      containsFullField, readHeightSfc_opt=readHeightSfc)

    !-- 2.0 Horizontal Interpolation

    ! initialize single precision 3D working copy of statevector for horizontal interpolation result
    call gsv_allocate(statevector_hinterp_r4, 1, statevector_out_r4%hco, vco_file,   &
                      dateStamp_opt=statevector_out_r4%datestamplist(stepIndex),     &
                      mpi_local_opt=.false., dataKind_opt=4,                         &
                      allocHeightSfc_opt=readHeightSfc, varNames_opt=varNamesToRead, &
                      hInterpolateDegree_opt=statevector_out_r4%hInterpolateDegree,  &
                      hExtrapolateDegree_opt=statevector_out_r4%hExtrapolateDegree)

    call gsv_hInterpolate_r4(statevector_file_r4, statevector_hinterp_r4)

    call gsv_deallocate(statevector_file_r4)

    !-- 3.0 Unit conversion (must come before vertical interp to get Psfc in Pascals)

    if ( unitConversion ) then
      call gio_fileUnitsToStateUnits( statevector_hinterp_r4, containsFullField )
    end if

    !-- 4.0 Vertical interpolation

    ! initialize double precision 3D working copy of statevector for mpi communication result
    call gsv_allocate(statevector_vinterp_r4, 1, statevector_out_r4%hco, statevector_out_r4%vco, &
                      dateStamp_opt=statevector_out_r4%datestamplist(stepIndex),                 &
                      mpi_local_opt=.false., dataKind_opt=4,                                     &
                      allocHeightSfc_opt=readHeightSfc, varNames_opt=varNamesToRead)

    call gsv_vInterpolate_r4(statevector_hinterp_r4,statevector_vinterp_r4)

    call gsv_deallocate(statevector_hinterp_r4)

    !-- 5.0 Copy result to output statevector
    call gsv_copy(statevector_vinterp_r4, statevector_out_r4, stepIndexOut_opt=stepIndex)

    call gsv_deallocate(statevector_vinterp_r4)
    deallocate(varNamesToRead)

    write(*,*) 'readFromFileAndInterp1Proc: END'

  end subroutine readFromFileAndInterp1Proc

  !--------------------------------------------------------------------------
  ! readFromFileOnly
  !--------------------------------------------------------------------------
  subroutine readFromFileOnly(statevector_out, fileName,  &
             etiket_in, typvar_in, stepIndex, unitConversion,  &
             readHeightSfc, containsFullField)
    implicit none

    ! arguments
    type(struct_gsv)              :: statevector_out
    character(len=*), intent(in)  :: fileName
    character(len=*), intent(in)  :: etiket_in
    character(len=*), intent(in)  :: typvar_in
    integer                       :: stepIndex
    logical                       :: unitConversion
    logical                       :: readHeightSfc
    logical                       :: containsFullField

    write(*,*) ''
    write(*,*) 'readFromFileOnly: Do simple reading with no interpolation and no mpi redistribution'

    if ( statevector_out%dataKind /= 4) then
      call utl_abort('readFromFileOnly: Only compatible with dataKind=4')
    end if

    call gio_readFile( statevector_out, filename, etiket_in, typvar_in,  &
                       containsFullField, readHeightSfc_opt=readHeightSfc, stepIndex_opt=stepIndex )

    if ( unitConversion ) then
      call gio_fileUnitsToStateUnits( statevector_out, containsFullField, stepIndex_opt=stepIndex )
    end if

    write(*,*) 'readFromFileOnly: END'

  end subroutine readFromFileOnly

  !--------------------------------------------------------------------------
  ! gio_readFile
  !--------------------------------------------------------------------------
  subroutine gio_readFile(statevector, filename, etiket_in, typvar_in, &
                          containsFullField, readHeightSfc_opt, stepIndex_opt, &
                          ignoreDate_opt)
    !
    ! :Purpose: Read an RPN standard file and put the contents into a
    !           stateVector object.
    !
    implicit none

    ! arguments
    type(struct_gsv),  intent(inout) :: statevector
    character(len=*),  intent(in)    :: fileName
    character(len=*),  intent(in)    :: etiket_in
    character(len=*),  intent(in)    :: typvar_in
    logical,           intent(in)    :: containsFullField
    logical, optional, intent(in)    :: readHeightSfc_opt
    integer, optional, intent(in)    :: stepIndex_opt
    logical, optional, intent(in)    :: ignoreDate_opt

    ! locals
    type(struct_hco), pointer :: hco_physics

    integer :: nulfile, ierr, ip1, ni_file, nj_file, nk_file, kIndex, stepIndex, ikey, levIndex
    integer :: stepIndexBeg, stepIndexEnd, ni_var, nj_var, nk_var
    integer :: fnom, fstouv, fclos, fstfrm, fstlir, fstinf
    integer :: fstprm, EZscintID_var, ezdefset, ezqkdef

    integer :: dateo_var, deet_var, npas_var, nbits_var, datyp_var
    integer :: ip1_var, ip2_var, ip3_var, swa_var, lng_var, dltf_var, ubc_var
    integer :: extra1_var, extra2_var, extra3_var
    integer :: ig1_var, ig2_var, ig3_var, ig4_var
    integer :: varIndex, dateStampList(statevector%numStep)

    character(len=4 ) :: nomvar_var
    character(len=2 ) :: typvar_var
    character(len=1 ) :: grtyp_var
    character(len=12) :: etiket_var

    real(4), pointer :: field_r4_ptr(:,:,:,:)
    real(4), pointer :: gd2d_file_r4(:,:)
    real(4), allocatable :: gd2d_var_r4(:,:)
    integer, allocatable :: mask(:,:)

    character(len=4)  :: varName, varNameToRead
    character(len=4)  :: varLevel

    type(struct_vco), pointer :: vco_file
    type(struct_hco), pointer :: hco_file
    logical :: foundVarNameInFile, ignoreDate

    vco_file => gsv_getVco(statevector)

    if ( statevector%mpi_distribution /= 'VarsLevs' .and. &
         statevector%mpi_local ) then
      call utl_abort('gio_readFile: statevector must have ' //   &
                     'complete horizontal fields on each mpi task.')
    end if

    if ( present(stepIndex_opt) ) then
      stepIndexBeg = stepIndex_opt
      stepIndexEnd = stepIndex_opt
    else
      stepIndexBeg = 1
      stepIndexEnd = statevector%numStep
    end if

    if ( present(ignoreDate_opt) ) then
      ignoreDate = ignoreDate_opt
    else
      ignoreDate = .false.
    end if

    if ( .not. associated( statevector%dateStampList )) then
      call utl_abort('gio_readFile: dateStampList of statevector is not associated with a target!')
    else
      dateStampList(:) = statevector%dateStampList(:)
      if (ignoreDate) then
        write(*,*) 'gio_readFile: as requested, ignoring the date when reading fields'
        dateStampList(:) = -1
      end if
    end if

    !- Open input field
    nulfile = 0
    write(*,*) 'gio_readFile: file name = ',trim(fileName)
    ierr = fnom(nulfile,trim(fileName),'RND+OLD+R/O',0)

    if ( ierr >= 0 ) then
      ierr  =  fstouv(nulfile,'RND+OLD')
    else
      call utl_abort('gio_readFile: problem opening input file')
    end if

    if (nulfile == 0 ) then
      call utl_abort('gio_readFile: unit number for input file not valid')
    end if

    ! Read surface height if requested
    if ( present(readHeightSfc_opt) ) then
      if ( readHeightSfc_opt .and. gsv_isAssocHeightSfc(statevector) ) then
        write(*,*) 'gio_readFile: reading the surface height'
        varName = 'GZ'
        ip1 = gsv_getVco(statevector)%ip1_sfc
        typvar_var = typvar_in
        ikey = fstinf(nulfile, ni_file, nj_file, nk_file,  &
                      -1, etiket_in, &
                      -1, -1, -1, typvar_var, varName)

        if ( ikey < 0 ) then
          if ( trim(typvar_in) /= "" ) then
            typvar_var(2:2) = '@'
            ikey = fstinf(nulfile, ni_file, nj_file, nk_file,  &
                 -1, etiket_in, &
                 -1, -1, -1, typvar_var, varName)
          end if
          if ( ikey < 0 ) then
            write(*,*) 'gio_readFile: etiket_in = ',etiket_in
            write(*,*) 'gio_readFile: typvar_in = ',typvar_in
            call utl_abort('gio_readFile: Problem with reading surface height from file')
          end if
        end if

        if ( ni_file /= gsv_getHco(statevector)%ni .or. nj_file /= gsv_getHco(statevector)%nj ) then
          write(*,*) 'ni, nj in file        = ', ni_file, nj_file
          write(*,*) 'ni, nj in statevector = ', gsv_getHco(statevector)%ni, gsv_getHco(statevector)%nj
          call utl_abort('gio_readFile: Dimensions of surface height not consistent')
        end if

        allocate(gd2d_file_r4(ni_file,nj_file))
        gd2d_file_r4(:,:) = 0.0d0
        ierr = fstlir(gd2d_file_r4(:,:),nulfile,ni_file, nj_file, nk_file,  &
                      -1,etiket_in,ip1,-1,-1,  &
                      typvar_var,varName)
        if ( ierr < 0 ) then
          write(*,*) 'ip1 = ',ip1
          write(*,*) 'etiket_in = ',etiket_in
          write(*,*) 'typvar_var = ',typvar_var
          call utl_abort('gio_readFile: Problem with reading surface height from file')
        end if
        call gsv_setHeightSfc(statevector, &
             real(gd2d_file_r4(1:gsv_getHco(statevector)%ni,1:gsv_getHco(statevector)%nj),8)*10.0d0)
        deallocate(gd2d_file_r4)
      end if
    end if

    nullify(hco_file)
    nullify(gd2d_file_r4)
    if ( statevector%mykCount > 0 ) then
      if (gsv_getHco(statevector)%global) then

        foundVarNameInFile = .false.
        do varIndex = 1, vnl_numvarmax
          varName = vnl_varNameList(varIndex)

          if (.not. gsv_varExist(statevector,varName)) cycle

          ! make sure variable is in the file
          if ( .not. utl_varNamePresentInFile(varName,fileName_opt=trim(fileName)) ) cycle

          ! adopt a variable on the full/dynamic LAM grid
          if ( (trim(varName) == 'TM'   .or. trim(varName) == 'MG' ) ) cycle

          foundVarNameInFile = .true.

          exit
        end do

        ! special case when only TM (Surface Temperature) is in the file:
        if ( .not. foundVarNameInFile ) then
          varname = 'TM'
          if ( gsv_varExist( statevector, varname ) .and. &
               utl_varNamePresentInFile( varname, fileName_opt = trim( fileName ))) &
            foundVarNameInFile = .true.
        end if   

        ! to be safe for situations where, e.g. someone wants to only read MG from a file
        if ( .not. foundVarNameInFile ) then
          varname = 'P0'
          if ( utl_varNamePresentInFile( varname, fileName_opt = trim( fileName ))) &
            foundVarNameInFile = .true.
        end if

        if ( .not. foundVarNameInFile) call utl_abort('gio_readFile: NO variable is in the file')

        call hco_setupFromFile(hco_file, filename, ' ', 'INPUTFILE', varName_opt=varName)

      else
        ! In LAM mode, force the input file dimensions to be always identical to the input statevector dimensions
        hco_file => gsv_getHco(statevector)

        ! Also attempt to set up the physics grid
        if (gsv_interpToPhysicsGrid) then
          var_loop: do varIndex = 1, vnl_numvarmax
            varName = vnl_varNameList(varIndex)
            hco_physics = gsv_getHco_physics(statevector)
            if ( .not. gsv_varExist(statevector,varName)) cycle var_loop
            if ( .not. vnl_isPhysicsVar(varName) ) cycle var_loop
            if ( utl_varNamePresentInFile(varName, fileName_opt=filename) .and. &
               .not. associated(hco_physics) ) then
              write(*,*) 'gio_readFile: set up physics grid using the variable:', varName
              call hco_SetupFromFile(hco_physics, filename, ' ', 'INPUTFILE', varName_opt=varName)
              exit var_loop
            end if
          end do var_loop
        end if

      end if
      allocate(gd2d_file_r4(hco_file%ni,hco_file%nj))
      gd2d_file_r4(:,:) = 0.0
    end if

    ! Read all other fields needed for this MPI task
    call gsv_getField(statevector,field_r4_ptr)
    do stepIndex = stepIndexBeg, stepIndexEnd
      k_loop: do kIndex = statevector%mykBeg, statevector%mykEnd
        varName = gsv_getVarNameFromK(statevector,kIndex)
        levIndex = gsv_getLevFromK(statevector,kIndex)

        if (.not.gsv_varExist(statevector,varName)) cycle k_loop

        ! Check that the wanted field is present in the file
        if (utl_varNamePresentInFile(varName,fileUnit_opt=nulfile)) then
          varNameToRead = varName
        else
          select case (trim(varName))
          case ('LVIS')
            varNameToRead = 'VIS'
          case ('Z_T','Z_M','P_T','P_M')
            cycle k_loop
          case ('LPR')
            varNameToRead = 'PR'
          case default
            call utl_abort('gio_readFile: variable '//trim(varName)//' was not found in '//trim(fileName))
          end select
        end if

        varLevel = vnl_varLevelFromVarname(varNameToRead)
        if (varLevel == 'MM') then
          ip1 = vco_file%ip1_M(levIndex)
        else if (varLevel == 'TH') then
          ip1 = vco_file%ip1_T(levIndex)
        else if (varLevel == 'SF') then
          ip1 = -1
        else if (varLevel == 'SFTH') then
          ip1 = vco_file%ip1_T_2m
        else if (varLevel == 'SFMM') then
          ip1 = vco_file%ip1_M_10m
        else if (varLevel == 'OT') then
          ip1 = vco_ip1_other(levIndex)
        else if (varLevel == 'DP') then
          ip1 = vco_file%ip1_depth(levIndex)
        else if (varLevel == 'SFDP') then
          ip1 = -1
        else
          write(*,*) 'varLevel =', varLevel
          call utl_abort('gio_readFile: unknown varLevel')
        end if

        typvar_var = typvar_in

        ! Make sure that the input variable has the same grid size than hco_file
        ikey = fstinf(nulfile, ni_var, nj_var, nk_var,         &
                      datestamplist(stepIndex), etiket_in, &
                      -1, -1, -1, typvar_var, varNameToRead)

        if ( ikey < 0 ) then
          if ( trim(typvar_in) /= "" ) then
            typvar_var(2:2) = '@'
            ikey = fstinf(nulfile, ni_var, nj_var, nk_var,         &
                          datestamplist(stepIndex), etiket_in, &
                          -1, -1, -1, typvar_var, varNameToRead)
          end if
          if (ikey < 0) then
            write(*,*) 'gio_readFile: looking for datestamp = ', datestamplist(stepIndex)
            write(*,*) 'gio_readFile: etiket_in = ',etiket_in
            write(*,*) 'gio_readFile: typvar_in = ',typvar_in
            call utl_abort('gio_readFile: cannot find field ' // trim(varNameToRead) // ' in file ' // trim(fileName))
          end if
        end if

        ierr = fstprm( ikey,                                                               & ! IN
                       dateo_var, deet_var, npas_var, ni_var, nj_var, nk_var, nbits_var,   & ! OUT
                       datyp_var, ip1_var, ip2_var, ip3_var, typvar_var, nomvar_var,       & ! OUT
                       etiket_var, grtyp_var, ig1_var, ig2_var, ig3_var, ig4_var, swa_var, & ! OUT
                       lng_var, dltf_var, ubc_var, extra1_var, extra2_var, extra3_var )      ! OUT
        statevector%deet                      = deet_var
        statevector%ip2List(stepIndex)        = ip2_var
        statevector%npasList(stepIndex)       = npas_var
        statevector%dateOriginList(stepIndex) = dateo_var
        statevector%etiket                    = etiket_var

        ! Check if we found a mask field by mistake - if yes, need to fix the code!
        if (typvar_var == '@@') then
          call utl_abort('gio_readFile: read a mask file by mistake - need to modify file or fix the code')
        end if

        if ( ni_var == hco_file%ni .and. nj_var == hco_file%nj ) then
          ierr = fstlir(gd2d_file_r4(:,:),nulfile,ni_file, nj_file, nk_file,  &
                        datestamplist(stepIndex),etiket_in,ip1,-1,-1,  &
                        typvar_var,varNameToRead)
        else
          ! Special cases for variables that are on a different horizontal grid in LAM (e.g. TG)
          write(*,*)
          write(*,*) 'gio_readFile: variable on a different horizontal grid = ',trim(varNameToRead)
          write(*,*) ni_var, hco_file%ni, nj_var, hco_file%nj
          if ( gsv_interpToPhysicsGrid ) then
            hco_physics = gsv_getHco_physics(statevector)
            if ( associated(hco_physics) ) then
              if ( ni_var == hco_physics%ni .and. &
                   nj_var == hco_physics%nj ) then
                write(*,*) 'gio_readFile: this variable on same grid as other physics variables'
                statevector%onPhysicsGrid(vnl_varListIndex(varName)) = .true.
              else
                call utl_abort('gio_readFile: this variable not on same grid as other physics variables')
              end if
            else
              call utl_abort('gio_readFile: physics grid has not been set up')
            end if
          end if

          if (gsv_getHco(statevector)%global) then
            call utl_abort('gio_readFile: This is not allowed in global mode!')
          end if

          EZscintID_var  = ezqkdef( ni_var, nj_var, grtyp_var, ig1_var, ig2_var, ig3_var, ig4_var, nulfile ) ! IN

          allocate(gd2d_var_r4(ni_var,nj_var))
          gd2d_var_r4(:,:) = 0.0

          ierr = fstlir(gd2d_var_r4(:,:),nulfile,ni_var, nj_var, nk_var,  &
                        datestamplist(stepIndex),etiket_in,ip1,-1,-1,  &
                        typvar_in,varNameToRead)

          ierr = ezdefset(hco_file%EZscintID,EZscintID_var)
          ierr = utl_ezsint( gd2d_file_r4, gd2d_var_r4, interpDegree='NEAREST', extrapDegree_opt='NEUTRAL' )

          ! read the corresponding mask if it exists
          if (typvar_var(2:2) == '@') then
            write(*,*) 'gio_readFile: read mask that needs interpolation for variable name: ', nomvar_var
            call utl_abort('gio_readFile: not implemented yet')
          end if

          deallocate(gd2d_var_r4)
        end if

        if (varNameToRead == varName .or. .not. containsFullField) then
          field_r4_ptr(:,:,kIndex,stepIndex) = gd2d_file_r4(1:gsv_getHco(statevector)%ni,1:gsv_getHco(statevector)%nj)
        else
          select case (trim(varName))
          case ('LVIS')
            field_r4_ptr(:,:,kIndex,stepIndex) = &
                 log(max(min(gd2d_file_r4(1:gsv_getHco(statevector)%ni,1:gsv_getHco(statevector)%nj),mpc_maximum_vis_r4),mpc_minimum_vis_r4))
          case ('LPR')
            field_r4_ptr(:,:,kIndex,stepIndex) = &
                 log(mpc_minimum_pr_r4 + max(0.0,gd2d_file_r4(1:gsv_getHco(statevector)%ni,1:gsv_getHco(statevector)%nj)))
          case default
            call utl_abort('gio_readFile: Oups! This should not happen... Check the code.')
          end select
        endif

        if (ierr.lt.0)then
          write(*,*) varNameToRead,ip1,datestamplist(stepIndex)
          call utl_abort('gio_readFile: Problem with reading file')
        end if

        ! When mpi distribution could put UU on a different mpi task than VV
        ! or only one wind component present in statevector
        ! then we re-read the corresponding UV component and store it
        if ( statevector%extraUVallocated ) then
          if (varName == 'UU') then
            ierr = fstlir(gd2d_file_r4(:,:),nulfile, ni_file, nj_file, nk_file,  &
                          datestamplist(stepIndex),etiket_in,ip1,-1,-1,  &
                          typvar_in,'VV')
            call gsv_setFieldUV(statevector, &
                   gd2d_file_r4(1:gsv_getHco(statevector)%ni,1:gsv_getHco(statevector)%nj), &
                   kIndex, stepIndex)
          else if (varName == 'VV') then
            ierr = fstlir(gd2d_file_r4(:,:),nulfile, ni_file, nj_file, nk_file,  &
                          datestamplist(stepIndex),etiket_in,ip1,-1,-1,  &
                          typvar_in,'UU')
            call gsv_setFieldUV(statevector, &
                   gd2d_file_r4(1:gsv_getHco(statevector)%ni,1:gsv_getHco(statevector)%nj), &
                   kIndex, stepIndex)
          end if
        end if

      end do k_loop
    end do

    if (gsv_getHco(statevector)%global .and. statevector%mykCount > 0) call hco_deallocate(hco_file)
    if (allocated(mask)) deallocate(mask)

    ierr = fstfrm(nulfile)
    ierr = fclos(nulfile)        
    if ( associated(gd2d_file_r4) ) deallocate(gd2d_file_r4)

    ! Read in an oceanMask if it is present in the file
    call gio_readMaskFromFile(statevector, trim(filename))

  end subroutine gio_readFile

  !--------------------------------------------------------------------------
  ! gio_readMaskFromFile
  !--------------------------------------------------------------------------
  subroutine gio_readMaskFromFile(stateVector, filename)
    !
    ! :Purpose: Check if any ocean mask fields exist. If so, read for the surface
    !           or all ocean depth levels.
    !
    implicit none

    ! arguments
    type(struct_gsv)              :: stateVector
    character(len=*), intent(in)  :: fileName

    call ocm_readMaskFromFile(stateVector%oceanMask,gsv_getHco(statevector), gsv_getVco(statevector), filename)

  end subroutine gio_readMaskFromFile

  !--------------------------------------------------------------------------
  ! gio_getMaskLAM
  !--------------------------------------------------------------------------
  subroutine gio_getMaskLAM(statevector_mask, hco_ptr, vco_ptr, hInterpolateDegree_opt)
    !:Purpose: To read a LAM mask from a file (./analinc_mask by default).
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(inout) :: statevector_mask
    type(struct_hco), intent(in), pointer :: hco_ptr
    type(struct_vco), intent(in), pointer :: vco_ptr
    character(len=*), intent(in), optional :: hInterpolateDegree_opt

    ! Locals
    character(len=12) :: hInterpolationDegree

    if (present(hInterpolateDegree_opt)) then
      hInterpolationDegree = hInterpolateDegree_opt
    else
      hInterpolationDegree = 'LINEAR'
    end if

    call gsv_allocate(statevector_mask, 1, hco_ptr, vco_ptr, dateStamp_opt=-1, &
                      dataKind_opt=pre_incrReal, &
                      mpi_local_opt=.true., varNames_opt=(/'MSKC'/),           &
                      hInterpolateDegree_opt=hInterpolationDegree)
    call gio_readFromFile(statevector_mask, './analinc_mask', ' ', ' ', unitConversion_opt=.false., &
                          vcoFileIn_opt=vco_ptr)

  end subroutine gio_getMaskLAM

  !--------------------------------------------------------------------------
  ! gio_fileUnitsToStateUnits
  !--------------------------------------------------------------------------
  subroutine gio_fileUnitsToStateUnits(statevector, containsFullField, stepIndex_opt)
    !
    !:Purpose: Unit conversion needed after reading rpn standard file
    !
    implicit none

    ! Arguments:
    type(struct_gsv)            :: statevector
    logical                     :: containsFullField
    integer, optional           :: stepIndex_opt

    ! Locals:
    real(4), pointer :: field_r4_ptr(:,:,:,:)
    real(4), pointer :: gdUV_r4(:,:,:)
    real(8) :: multFactor
    integer :: stepIndex, stepIndexBeg, stepIndexEnd, kIndex
    character(len=4) :: varName
    

    if ( present(stepIndex_opt) ) then
      stepIndexBeg = stepIndex_opt
      stepIndexEnd = stepIndex_opt
    else
      stepIndexBeg = 1
      stepIndexEnd = statevector%numStep
    end if

    call gsv_getField(statevector,field_r4_ptr)

    step_loop: do stepIndex = stepIndexBeg, stepIndexEnd

      ! Do unit conversion for all variables
      do kIndex = statevector%mykBeg, statevector%mykEnd
        varName = gsv_getVarNameFromK(statevector,kIndex)

        if ( trim(varName) == 'UU' .or. trim(varName) == 'VV') then
          multFactor = mpc_m_per_s_per_knot_r8 ! knots -> m/s
        else if ( trim(varName) == 'P0' ) then
          multFactor = mpc_pa_per_mbar_r8 ! hPa -> Pa
        else if ( vnl_varKindFromVarname(trim(varName)) == 'CH' ) then 
          if ( gsv_conversionVarKindCHtoMicrograms ) then
            if ( trim(varName) == 'TO3' .or. trim(varName) == 'O3L' ) then
              ! Convert from volume mixing ratio to micrograms/kg
              ! Standard ozone input would not require this conversion as it is already in micrograms/kg
              multFactor = 1.0d9*vnl_varMassFromVarName(trim(varName)) &
                           /mpc_molar_mass_dry_air_r8 ! vmr -> micrograms/kg
            else
              multFactor = 1.0d0 ! no conversion
            end if
          else
            multFactor = 1.0d0 ! no conversion
          end if
        else
          multFactor = 1.0d0 ! no conversion
        end if
        write(*,*) 'DBGmad conversion ', varName, multFactor

        if ( multFactor /= 1.0d0 ) then
          field_r4_ptr(:,:,kIndex,stepIndex) = real( multFactor * field_r4_ptr(:,:,kIndex,stepIndex), 4 )
        end if

        if ( trim(varName) == 'TT' .and. containsFullField ) then
          field_r4_ptr(:,:,kIndex,stepIndex) = real( field_r4_ptr(:,:,kIndex,stepIndex) +  &
                                                     mpc_k_c_degree_offset_r8, 4 )
        end if

        if ( trim(varName) == 'TM' .and. containsFullField ) then
          if (maxval(field_r4_ptr(:,:,kIndex,stepIndex)) < 50.0) then
            field_r4_ptr(:,:,kIndex,stepIndex) = real( field_r4_ptr(:,:,kIndex,stepIndex) + &
                                                       mpc_k_c_degree_offset_r8, 4 )
          end if
        end if

        if ( trim(varName) == 'VIS' .and. containsFullField ) then
          field_r4_ptr(:,:,kIndex,stepIndex) = min(field_r4_ptr(:,:,kIndex,stepIndex),mpc_maximum_vis_r4)
        end if

        if ( vnl_varKindFromVarname(trim(varName)) == 'CH' .and. containsFullField ) then 
          if ( gsv_minValVarKindCH(vnl_varListIndex(varName)) > 1.01*mpc_missingValue_r8 ) &
            field_r4_ptr(:,:,kIndex,stepIndex) = max( field_r4_ptr(:,:,kIndex,stepIndex), &
              real(gsv_minValVarKindCH(vnl_varListIndex(trim(varName)))) )
        end if

        if ( trim(varName) == 'PR' .and. containsFullField ) then
          field_r4_ptr(:,:,kIndex,stepIndex) = max(field_r4_ptr(:,:,kIndex,stepIndex),0.0)
        end if
      end do

      ! Do unit conversion for extra copy of winds, if present
      if ( statevector%extraUVallocated ) then
        multFactor = mpc_m_per_s_per_knot_r8 ! knots -> m/s

        !$OMP PARALLEL DO PRIVATE (kIndex)
        do kIndex = statevector%myUVkBeg, statevector%myUVkEnd
          call gsv_getFieldUV(statevector, gdUV_r4, kIndex)
          call gsv_setFieldUV(statevector, &
                 real( multFactor * gdUV_r4(:,:,stepIndex), 4 ), kIndex, stepIndex)
        end do
        !$OMP END PARALLEL DO

      end if

    end do step_loop

  end subroutine gio_fileUnitsToStateUnits

  !--------------------------------------------------------------------------
  ! gio_readTrials
  !--------------------------------------------------------------------------
  subroutine gio_readTrials(stateVectorTrialIn)
    !
    !:Purpose: Reading trials
    !
    implicit none

    ! Arguments
    type(struct_gsv), target, intent(inout) :: stateVectorTrialIn

    ! Locals
    type(struct_gsv),  target :: stateVectorTrial
    type(struct_gsv), pointer :: stateVectorTrial_ptr 
    type(struct_gsv)     :: stateVector_1step_r4
    integer              :: fnom, fstouv, fclos, fstfrm, fstinf
    integer              :: ierr, ikey, stepIndex, stepIndexToRead, trialIndex, nulTrial
    integer, parameter   :: maxNumTrials = 100
    integer              :: ni_file, nj_file, nk_file, dateStamp, varNameIndex
    integer              :: procToRead, numBatch, batchIndex, stepIndexBeg, stepIndexEnd
    character(len=2)     :: fileNumber
    character(len=512)   :: fileName
    logical              :: fileExists, allocHeightSfc
    logical              :: useInputStateVectorTrial 
    character(len=4), pointer :: varNamesToRead(:)
    character(len=4)     :: varNameForDateStampSearch

    call tmg_start(150,'gio_readTrials')

    if ( mpi_myid == 0 ) then
      write(*,*) ''
      write(*,*) 'gio_readTrials: STARTING'
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    end if

    if ( gsv_varExist(stateVectorTrialIn,'Z_T') .or. &
         gsv_varExist(stateVectorTrialIn,'Z_M') .or. &
         gsv_varExist(stateVectorTrialIn,'P_T') .or. &
         gsv_varExist(stateVectorTrialIn,'P_M') ) then

      useInputStateVectorTrial = .false.

      allocHeightSfc = ( stateVectorTrialIn%vco%Vcode /= 0 )

      ! Allocate single-precision statevector without Z/P to read trials
      call gsv_allocate( stateVectorTrial, stateVectorTrialIn%numStep, &
                         stateVectorTrialIn%hco, stateVectorTrialIn%vco, &
                         dateStamp_opt=tim_getDateStamp(), &
                         mpi_local_opt=stateVectorTrialIn%mpi_local, &
                         mpi_distribution_opt='Tiles', dataKind_opt=4,  &
                         allocHeightSfc_opt=allocHeightSfc, &
                         hInterpolateDegree_opt=stateVectorTrialIn%hInterpolateDegree, &
                         allocHeight_opt=.false., allocPressure_opt=.false., &
                         beSilent_opt=.false. )
      call gsv_zero( stateVectorTrial )
      stateVectorTrial_ptr => stateVectorTrial
    else
      useInputStateVectorTrial = .true.

      stateVectorTrial_ptr => stateVectorTrialIn
    end if

    nullify(varNamesToRead)
    call gsv_varNamesList(varNamesToRead, stateVectorTrial_ptr)

    varNameForDateStampSearch = ' '
    do varNameIndex = 1, size(varNamesToRead)
      select case (trim(varNamesToRead(varNameIndex)))
      case ('Z_T','Z_M','P_T','P_M')
        cycle
      case default
        varNameForDateStampSearch = varNamesToRead(varNameIndex)
        exit
      end select
    end do

    ! warn if not enough mpi tasks
    if ( mpi_nprocs < stateVectorTrial_ptr%numStep ) then
      write(*,*) 'gio_readTrials: number of trial time steps, mpi tasks = ', stateVectorTrial_ptr%numStep, mpi_nprocs
      write(*,*) 'gio_readTrials: for better efficiency, the number of mpi tasks should '
      write(*,*) '                be at least as large as number of trial time steps'
    end if

    allocHeightSfc = stateVectorTrial_ptr%heightSfcPresent

    ! figure out number of batches of time steps for reading
    numBatch = ceiling(real(stateVectorTrial_ptr%numStep) / real(mpi_nprocs))
    write(*,*) 'gio_readTrials: reading will be done by number of batches = ', numBatch

    BATCH: do batchIndex = 1, numBatch

      stepIndexBeg = 1 + (batchIndex - 1) * mpi_nprocs
      stepIndexEnd = min(stateVectorTrial_ptr%numStep, stepIndexBeg + mpi_nprocs - 1)
      write(*,*) 'gio_readTrials: batchIndex, stepIndexBeg/End = ', batchIndex, stepIndexBeg, stepIndexEnd

      ! figure out which time step I will read, if any (-1 if none)
      stepIndexToRead = -1
      do stepIndex = stepIndexBeg, stepIndexEnd
        procToRead = nint( real(stepIndex - stepIndexBeg) * real(mpi_nprocs) / real(stepIndexEnd - stepIndexBeg + 1) )
        if ( procToRead == mpi_myid ) stepIndexToRead = stepIndex
        if ( mpi_myid == 0 ) write(*,*) 'gio_readTrials: stepIndex, procToRead = ', stepIndex, procToRead
      end do

      ! loop over all times for which stateVector is allocated
      if ( stepIndexToRead /= -1 ) then
        dateStamp = stateVectorTrial_ptr%dateStampList(stepIndexToRead)
        write(*,*) 'gio_readTrials: reading background for time step: ',stepIndexToRead, dateStamp
        write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

        ! identify which trial file corresponds with current datestamp
        ikey = 0
        do trialIndex = 1, maxNumTrials
          write(fileNumber,'(I2.2)') trialIndex
          fileName = 'trlm_' // trim(fileNumber)
          inquire(file=trim(fileName),exist=fileExists)
          if ( .not. fileExists ) exit
          nulTrial = 0
          ierr = fnom(nulTrial,trim(fileName),'RND+OLD+R/O',0)
          ierr = fstouv(nulTrial,'RND+OLD')
          ikey = fstinf(nulTrial, ni_file, nj_file, nk_file,  &
                        dateStamp, ' ', -1, -1, -1, ' ', varNameForDateStampSearch)
          ierr = fstfrm(nulTrial)
          ierr = fclos(nulTrial)
          if ( ikey > 0 ) exit
        end do

        if ( ikey <= 0 .or. .not.fileExists ) then 
          write(*,*) 'stepIndexToRead, dateStamp = ', stepIndexToRead, dateStamp
          call utl_abort('gio_readTrials: trial file not found for this increment timestep')
        end if

        ! allocate stateVector for storing just 1 time step
        if ( batchIndex == 1 ) then
          call gsv_allocate( stateVector_1step_r4, 1, stateVectorTrial_ptr%hco, stateVectorTrial_ptr%vco, &
                             dateStamp_opt=dateStamp, mpi_local_opt=.false., dataKind_opt=4,        &
                             allocHeightSfc_opt=allocHeightSfc, varNames_opt=varNamesToRead,        &
                             hInterpolateDegree_opt=stateVectorTrial_ptr%hInterpolateDegree,           &
                             hExtrapolateDegree_opt=stateVectorTrial_ptr%hExtrapolateDegree)
          call gsv_zero( stateVector_1step_r4 )
        else
          call gsv_modifyDate( stateVector_1step_r4, dateStamp )
        end if

        ! read the trial file for this timestep
        fileName = ram_fullWorkingPath(fileName)
        call gio_readFromFile(stateVector_1step_r4, fileName, ' ', 'P',  &
                              readHeightSfc_opt=allocHeightSfc)

        ! remove from ram disk to save some space
        ierr = ram_remove(fileName)
      else

        if ( gsv_isAllocated(stateVector_1step_r4) ) call gsv_deallocate(stateVector_1step_r4)

      end if ! I read a time step

      if ( stateVectorTrial_ptr%mpi_distribution == 'VarsLevs' ) then
        call gsv_transposeStepToVarsLevs(stateVector_1step_r4, stateVectorTrial_ptr, stepIndexBeg)
      else if ( stateVectorTrial_ptr%mpi_distribution == 'Tiles' ) then
        call gsv_transposeStepToTiles(stateVector_1step_r4, stateVectorTrial_ptr, stepIndexBeg)
      else
        call utl_abort( 'gio_readTrials: not compatible with mpi_distribution = ' // &
                        trim(stateVectorTrial_ptr%mpi_distribution) )
      end if

      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      if ( gsv_isAllocated(stateVector_1step_r4) .and. batchIndex == numBatch ) then
        call gsv_deallocate(stateVector_1step_r4)
      end if

    end do BATCH

    if ( .not. useInputStateVectorTrial ) then
      call gsv_copy( stateVectorTrial_ptr, stateVectorTrialIn, allowVarMismatch_opt=.true. )
      call gsv_deallocate( stateVectorTrial )
    end if

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'gio_readTrials: FINISHED'
    write(*,*) ''

    call tmg_stop(150)

  end subroutine gio_readTrials

end module gridStateVectorFileIO_mod
