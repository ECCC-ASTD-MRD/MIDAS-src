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
  ! MODULE gridStateVectorFile_mod (prefix='gio' category='4. Data Object transformations')
  !
  ! :Purpose: The grid-point state vector I/O methods.
  !
  use mpi
  use midasMpi_mod
  use gridStateVector_mod
  use interpolation_mod
  use utilities_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use oceanMask_mod
  use varNameList_mod
  use ramDisk_mod
  use timeCoord_mod
  use mathPhysConstants_mod
  use codePrecision_mod
  use Vgrid_Descriptors
  implicit none
  save
  private

  ! public subroutines and functions
  public :: gio_readFromFile, gio_readTrials, gio_readFile
  public :: gio_readMaskFromFile
  public :: gio_getMaskLAM
  public :: gio_writeToFile
  public :: gio_fileUnitsToStateUnits

  integer, external :: get_max_rss

  contains
  !--------------------------------------------------------------------------
  ! gio_readFromFile
  !--------------------------------------------------------------------------
  subroutine gio_readFromFile(statevector_out, fileName, etiket_in, typvar_in, &
                              stepIndex_opt, unitConversion_opt, &
                              statevectorRef_opt, readHeightSfc_opt, &
                              containsFullField_opt, vcoFileIn_opt)
    !
    ! :Purpose: Read an RPN standard file and put the contents into a
    !           stateVector object. Main high level wrapper subroutine.
    !
    implicit none

    ! arguments
    type(struct_gsv),                    intent(inout) :: statevector_out
    character(len=*),                    intent(in)    :: fileName
    character(len=*),                    intent(in)    :: etiket_in
    character(len=*),                    intent(in)    :: typvar_in
    integer,          optional,          intent(in)    :: stepIndex_opt
    logical,          optional,          intent(in)    :: unitConversion_opt
    type(struct_gsv), optional,          intent(in)    :: statevectorRef_opt ! Reference statevector providing optional fields (P0, TT, HU)
    logical,          optional,          intent(in)    :: readHeightSfc_opt
    logical,          optional,          intent(in)    :: containsFullField_opt
    type(struct_vco), optional, pointer, intent(in)    :: vcoFileIn_opt

    ! locals
    logical           :: doHorizInterp, doVertInterp, unitConversion
    logical           :: readHeightSfc, containsFullField
    logical           :: foundVarNameInFile 
    integer           :: stepIndex, varIndex
    character(len=4)  :: varName

    type(struct_vco), pointer :: vco_file
    type(struct_hco), pointer :: hco_file

    nullify(vco_file, hco_file)

    write(*,*)
    write(*,*) 'gio_readFromFile: START'
    call utl_tmg_start(160,'low-level--gsv_readFromFile')

    if ( present(stepIndex_opt) ) then
      stepIndex = stepIndex_opt
    else
      stepIndex = statevector_out%anltime
    end if
    if (stepIndex > stateVector_out%numStep .or. stepIndex < 1) then
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
    if (present(vcoFileIn_opt)) then
      vco_file => vcoFileIn_opt
      write(*,*)
      write(*,*) 'gio_readFromFile: vertical levels defined in user-supplied vco object will be read'
    else
      call vco_setupFromFile(vco_file,trim(fileName), beSilent_opt=.true.)
      write(*,*)
      write(*,*) 'gio_readFromFile: all the vertical levels will be read from ', trim(fileName) 
    end if

    foundVarNameInFile = .false.

    do varIndex = 1, vnl_numvarmax
      varName = vnl_varNameList(varIndex)

      if (.not. gsv_varExist(statevector_out,varName)) cycle

      ! make sure variable is in the file
      if (.not. utl_varNamePresentInFile(varName,fileName_opt=trim(fileName))) cycle

      ! adopt a variable on the full/dynamic LAM grid
      if (.not. statevector_out%hco%global .and. (trim(varName) == 'TM' .or. trim(varName) == 'MG')) cycle

      foundVarNameInFile = .true.

      exit

    end do

    ! special case when only TM (Surface Temperature) is in the file:
    if (.not. foundVarNameInFile) then
      varname = 'TM'
      if (gsv_varExist( statevector_out, varname) .and. &
          utl_varNamePresentInFile( varname, fileName_opt = trim( fileName ))) &
        foundVarNameInFile = .true.
    end if   

    ! to be safe for situations where, e.g. someone wants to only read MG from a file
    if (.not. foundVarNameInFile) then
      varname = 'P0'
      if (utl_varNamePresentInFile( varname, fileName_opt = trim( fileName))) &
        foundVarNameInFile = .true.
    end if   

    if (.not. foundVarNameInFile) call utl_abort('gio_readFromFile: NO variables found in the file!!!')

    write(*,*) 'gio_readFromFile: defining hco by varname= ', varName

    call hco_setupFromFile(hco_file, trim(fileName), etiket_in, gridName_opt='FILEGRID', varName_opt = varName)

    ! test if horizontal and/or vertical interpolation needed for statevector grid
    doVertInterp = .not.vco_equal(vco_file,statevector_out%vco)
    doHorizInterp = .not.hco_equal(hco_file,statevector_out%hco)
    write(*,*) 'gio_readFromFile: doVertInterp = ', doVertInterp, ', doHorizInterp = ', doHorizInterp

    ! call appropriate subroutine to do actual work
    if ((doVertInterp .or. doHorizInterp) .and. statevector_out%mpi_distribution=='Tiles') then
      call readFromFileAndInterpToTiles(statevector_out, fileName,  &
             vco_file, hco_file, etiket_in, typvar_in, stepIndex, unitConversion,  &
             readHeightSfc, containsFullField, statevectorRef_opt=statevectorRef_opt)
    else if ((doVertInterp .or. doHorizInterp) .and. .not.stateVector_out%mpi_local) then
      call readFromFileAndInterp1Proc(statevector_out, fileName,  &
             vco_file, hco_file, etiket_in, typvar_in, stepIndex, unitConversion,  &
             readHeightSfc, containsFullField)
    else if (.not.(doVertInterp .or. doHorizInterp) .and. stateVector_out%mpi_local) then
      call readFromFileAndTransposeToTiles(statevector_out, fileName,  &
             etiket_in, typvar_in, stepIndex, unitConversion,  &
             readHeightSfc, containsFullField)
    else
      call readFromFileOnly(statevector_out, fileName,  &
                                etiket_in, typvar_in, stepIndex, unitConversion,  &
                                readHeightSfc, containsFullField)
    end if

    call utl_tmg_stop(160)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'gio_readFromFile: END'

  end subroutine gio_readFromFile

  !--------------------------------------------------------------------------
  ! readFromFileAndInterpToTiles
  !--------------------------------------------------------------------------
  subroutine readFromFileAndInterpToTiles(statevector_out, fileName,  &
             vco_file, hco_file, etiket_in, typvar_in, stepIndex, unitConversion,  &
             readHeightSfc, containsFullField, statevectorRef_opt)
    !
    ! :Purpose: Read an RPN standard file and put the contents into a
    !           stateVector object. Wrapper subroutine that also proceed with
    !           distributed interpolation on MPI tiles.
    !
    ! :Note: this routine currently only works correctly for reading FULL FIELDS,
    !        not increments or perturbations... because of the HU -> LQ conversion
    implicit none

    ! arguments
    type(struct_gsv),           intent(inout) :: statevector_out
    character(len=*),           intent(in)    :: fileName
    type(struct_vco), pointer,  intent(in)    :: vco_file
    type(struct_hco), pointer,  intent(in)    :: hco_file
    character(len=*),           intent(in)    :: etiket_in
    character(len=*),           intent(in)    :: typvar_in
    integer,                    intent(in)    :: stepIndex
    logical,                    intent(in)    :: unitConversion
    logical,                    intent(in)    :: readHeightSfc
    logical,                    intent(in)    :: containsFullField
    type(struct_gsv), optional, intent(in)    :: statevectorRef_opt ! Reference statevector providing optional fields (P0, TT, HU)

    ! locals
    real(4), pointer     :: field3d_r4_ptr(:,:,:)
    real(8), pointer     :: field_in_ptr(:,:,:,:), field_out_ptr(:,:,:,:)

    character(len=4), pointer :: varNamesToRead(:)

    type(struct_gsv) :: statevector_file_r4, statevector_tiles 
    type(struct_gsv) :: statevector_hinterp_r4, statevector_vinterp

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

    call int_hInterp_gsv(statevector_file_r4, statevector_hinterp_r4)

    call gsv_deallocate(statevector_file_r4)

    !-- 3.0 Unit conversion

    if ( unitConversion ) then
      call gio_fileUnitsToStateUnits( statevector_hinterp_r4, containsFullField )
    end if

    !-- 4.0 MPI communication from vars/levels to lat/lon tiles

    ! initialize double precision 3D working copy of statevector for mpi communication result
    call gsv_allocate(statevector_tiles, 1, statevector_out%hco, vco_file,    &
                      dateStamp_opt=statevector_out%datestamplist(stepIndex), &
                      mpi_local_opt=.true., mpi_distribution_opt='Tiles',     &
                      dataKind_opt=8, allocHeightSfc_opt=readHeightSfc,               &
                      varNames_opt=varNamesToRead )

    call gsv_transposeVarsLevsToTiles(statevector_hinterp_r4, statevector_tiles)

    call gsv_deallocate(statevector_hinterp_r4)

    !-- 5.0 Vertical interpolation

    ! initialize double precision 3D working copy of statevector for mpi communication result
    call gsv_allocate(statevector_vinterp, 1, statevector_out%hco, statevector_out%vco, &
                      dateStamp_opt=statevector_out%datestamplist(stepIndex), &
                      mpi_local_opt=.true., mpi_distribution_opt='Tiles', dataKind_opt=8, &
                      allocHeightSfc_opt=readHeightSfc, varNames_opt=varNamesToRead )

    call int_vInterp_gsv( statevector_tiles, statevector_vinterp, & 
                          statevectorRef_opt=statevectorRef_opt)

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
    !
    ! :Purpose: Read an RPN standard file and put the contents into a
    !           stateVector object. Wrapper subroutine that also proceed with
    !           distributed transposition on MPI tiles.
    !
    ! :Note: this routine currently only works correctly for reading FULL FIELDS,
    !        not increments or perturbations... because of the HU -> LQ conversion
    implicit none

    ! arguments
    type(struct_gsv), intent(inout) :: statevector_out
    character(len=*), intent(in)    :: fileName
    character(len=*), intent(in)    :: etiket_in
    character(len=*), intent(in)    :: typvar_in
    integer,          intent(in)    :: stepIndex
    logical,          intent(in)    :: unitConversion
    logical,          intent(in)    :: readHeightSfc
    logical,          intent(in)    :: containsFullField

    ! locals
    real(4), pointer     :: field3d_r4_ptr(:,:,:)
    real(8), pointer     :: field_in_ptr(:,:,:,:), field_out_ptr(:,:,:,:)

    character(len=4), pointer :: varNamesToRead(:)

    type(struct_gsv) :: statevector_file_r4, statevector_tiles

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

    call gsv_deallocate(statevector_tiles)
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
    !
    ! :Purpose: Read an RPN standard file and put the contents into a
    !           stateVector object. Wrapper subroutine that also proceed with
    !           (serial) interpolation.
    !
    implicit none

    ! arguments
    type(struct_gsv),          intent(inout)  :: statevector_out_r4
    character(len=*),          intent(in)     :: fileName
    type(struct_vco), pointer, intent(in)     :: vco_file
    type(struct_hco), pointer, intent(in)     :: hco_file
    character(len=*),          intent(in)     :: etiket_in
    character(len=*),          intent(in)     :: typvar_in
    integer,                   intent(in)     :: stepIndex
    logical,                   intent(in)     :: unitConversion
    logical,                   intent(in)     :: readHeightSfc
    logical,                   intent(in)     :: containsFullField

    ! locals
    real(4), pointer :: field_in_ptr(:,:,:,:), field_out_ptr(:,:,:,:)

    character(len=4), pointer :: varNamesToRead(:)

    type(struct_gsv) :: statevector_file_r4, statevector_hinterp_r4
    type(struct_gsv) :: statevector_vinterp_r4

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

    call int_hInterp_gsv(statevector_file_r4, statevector_hinterp_r4)

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

    call int_vInterp_gsv(statevector_hinterp_r4,statevector_vinterp_r4)

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
    !
    ! :Purpose: Read an RPN standard file and put the contents into a
    !           stateVector object.  Wrapper subroutine
    !
    implicit none

    ! arguments
    type(struct_gsv), intent(inout) :: statevector_out
    character(len=*), intent(in)    :: fileName
    character(len=*), intent(in)    :: etiket_in
    character(len=*), intent(in)    :: typvar_in
    integer,          intent(in)    :: stepIndex
    logical,          intent(in)    :: unitConversion
    logical,          intent(in)    :: readHeightSfc
    logical,          intent(in)    :: containsFullField

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
    !           stateVector object.  Low level subroutine that does the actual
    !           file reading.
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
    integer :: nulfile, ierr, ip1, ni_file, nj_file, nk_file, kIndex, stepIndex
    integer :: ikey, levIndex
    integer :: stepIndexBeg, stepIndexEnd, ni_var, nj_var, nk_var
    integer :: fnom, fstouv, fclos, nulnam, fstfrm, fstlir, fstinf
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
    real(4), pointer :: gd2d_file_r4(:,:), gd2d_r4_UV_ptr(:,:,:)
    real(8), pointer :: heightSfc_ptr(:,:)
    real(4), allocatable :: gd2d_var_r4(:,:)

    character(len=4)  :: varName, varNameToRead
    character(len=4)  :: varLevel

    type(struct_vco), pointer :: vco_file
    type(struct_hco), pointer :: hco_file

    logical :: foundVarNameInFile, ignoreDate

    ! Namelist variables
    logical :: interpToPhysicsGrid

    NAMELIST /NAMSTIO/interpToPhysicsGrid

    write(*,*) 'gio_readFile: starting'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    interpToPhysicsGrid = .false.
    if ( .not. utl_isNamelistPresent('NAMSTIO','./flnml') ) then
      if ( mmpi_myid == 0 ) then
        write(*,*) 'gio_readFile: namstio is missing in the namelist. The default values will be taken.'
      end if
    else
      ! Read namelist NAMSTIO
      nulnam=0
      ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam,nml=namstio,iostat=ierr)
      if (ierr.ne.0) call utl_abort('gio_readfile: Error reading namelist')
      if (mmpi_myid.eq.0) write(*,nml=namstio)
      ierr=fclos(nulnam)
    end if

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

    if (.not. associated(statevector%dateStampList)) then
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

    if (ierr >= 0) then
      ierr  =  fstouv(nulfile,'RND+OLD')
    else
      call utl_abort('gio_readFile: problem opening input file')
    end if

    if (nulfile == 0) then
      call utl_abort('gio_readFile: unit number for input file not valid')
    end if

    ! Read surface height if requested
    if (present(readHeightSfc_opt)) then
      if (readHeightSfc_opt .and. gsv_isAssocHeightSfc(statevector)) then
        write(*,*) 'gio_readFile: reading the surface height'
        varName = 'GZ'
        ip1 = statevector%vco%ip1_sfc
        typvar_var = typvar_in
        ikey = fstinf(nulfile, ni_file, nj_file, nk_file,  &
                      -1, etiket_in, &
                      -1, -1, -1, typvar_var, varName)

        if (ikey < 0) then
          if (trim(typvar_in) /= "") then
            typvar_var(2:2) = '@'
            ikey = fstinf(nulfile, ni_file, nj_file, nk_file,  &
                 -1, etiket_in, &
                 -1, -1, -1, typvar_var, varName)
          end if
          if (ikey < 0) then
            write(*,*) 'gio_readFile: etiket_in = ', etiket_in
            write(*,*) 'gio_readFile: typvar_in = ', typvar_in
            call utl_abort('gio_readFile: Problem with reading surface height from file')
          end if
        end if

        if (ni_file /= statevector%hco%ni .or. nj_file /= statevector%hco%nj) then
          write(*,*) 'ni, nj in file        = ', ni_file, nj_file
          write(*,*) 'ni, nj in statevector = ', statevector%hco%ni, statevector%hco%nj
          call utl_abort('gio_readFile: Dimensions of surface height not consistent')
        end if

        allocate(gd2d_file_r4(ni_file,nj_file))
        gd2d_file_r4(:,:) = 0.0d0
        ierr = fstlir(gd2d_file_r4(:,:), nulfile, ni_file, nj_file, nk_file,  &
                      -1,etiket_in,ip1,-1,-1,  &
                      typvar_var,varName)
        if (ierr < 0) then
          write(*,*) 'ip1 = ', ip1
          write(*,*) 'etiket_in = ', etiket_in
          write(*,*) 'typvar_var = ', typvar_var
          call utl_abort('gio_readFile: Problem with reading surface height from file')
        end if
        heightSfc_ptr => gsv_getHeightSfc(statevector)
        heightSfc_ptr = real(gd2d_file_r4(1:gsv_getHco(statevector)%ni,1:gsv_getHco(statevector)%nj),8)*10.0d0
        deallocate(gd2d_file_r4)
      end if
    end if

    nullify(hco_file)
    nullify(gd2d_file_r4)
    if ( statevector%mykCount > 0 ) then
      if (statevector%hco%global) then

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
          if (gsv_varExist( statevector, varname ) .and. &
              utl_varNamePresentInFile( varname, fileName_opt = trim( fileName ))) &
            foundVarNameInFile = .true.
        end if   

        ! to be safe for situations where, e.g. someone wants to only read MG from a file
        if (.not. foundVarNameInFile) then
          varname = 'P0'
          if (utl_varNamePresentInFile( varname, fileName_opt = trim( fileName))) &
            foundVarNameInFile = .true.
        end if

        if (.not. foundVarNameInFile) call utl_abort('gio_readFile: NO variable is in the file')

        call hco_setupFromFile(hco_file, filename, ' ', 'INPUTFILE', varName_opt=varName)

      else
        ! In LAM mode, force the input file dimensions to be always identical to the input statevector dimensions
        hco_file => statevector%hco

        ! Also attempt to set up the physics grid
        if (interpToPhysicsGrid) then
          var_loop: do varIndex = 1, vnl_numvarmax
            varName = vnl_varNameList(varIndex)
            if (.not. gsv_varExist(statevector,varName)) cycle var_loop
            if (.not. vnl_isPhysicsVar(varName)) cycle var_loop
            if (utl_varNamePresentInFile(varName, fileName_opt=filename) .and. &
               .not. associated(statevector%hco_physics)) then
              write(*,*) 'gio_readFile: set up physics grid using the variable:', varName
              call hco_SetupFromFile(statevector%hco_physics, filename, ' ', &
                                     'INPUTFILE', varName_opt=varName)
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
        else if (varLevel == 'SS') then
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

        if (ikey < 0) then
          if (trim(typvar_in) /= '') then
            typvar_var(2:2) = '@'
            ikey = fstinf(nulfile, ni_var, nj_var, nk_var, &
                          datestamplist(stepIndex), etiket_in, &
                          -1, -1, -1, typvar_var, varNameToRead)
          end if
          if (ikey < 0) then
            write(*,*) 'gio_readFile: looking for datestamp = ', datestamplist(stepIndex)
            write(*,*) 'gio_readFile: etiket_in = ', etiket_in
            write(*,*) 'gio_readFile: typvar_in = ', typvar_in
            call utl_abort('gio_readFile: cannot find field ' // trim(varNameToRead) // ' in file ' // trim(fileName))
          end if
        end if

        ierr = fstprm(ikey,                                                               & ! IN
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

        if (ni_var == hco_file%ni .and. nj_var == hco_file%nj) then
          ierr = fstlir(gd2d_file_r4(:,:),nulfile,ni_file, nj_file, nk_file,  &
                        datestamplist(stepIndex),etiket_in,ip1,-1,-1,  &
                        typvar_var,varNameToRead)
        else
          ! Special cases for variables that are on a different horizontal grid in LAM (e.g. TG)
          write(*,*)
          write(*,*) 'gio_readFile: variable on a different horizontal grid = ',trim(varNameToRead)
          write(*,*) ni_var, hco_file%ni, nj_var, hco_file%nj
          if (interpToPhysicsGrid) then
            if (associated(statevector%hco_physics)) then
              if (ni_var == statevector%hco_physics%ni .and. &
                  nj_var == statevector%hco_physics%nj) then
                write(*,*) 'gio_readFile: this variable on same grid as other physics variables'
                statevector%onPhysicsGrid(vnl_varListIndex(varName)) = .true.
              else
                call utl_abort('gio_readFile: this variable not on same grid as other physics variables')
              end if
            else
              call utl_abort('gio_readFile: physics grid has not been set up')
            end if
          end if

          if (statevector%hco%global) then
            call utl_abort('gio_readFile: This is not allowed in global mode!')
          end if

          EZscintID_var  = ezqkdef(ni_var, nj_var, grtyp_var, ig1_var, ig2_var, ig3_var, ig4_var, nulfile) ! IN

          allocate(gd2d_var_r4(ni_var,nj_var))
          gd2d_var_r4(:,:) = 0.0

          ierr = fstlir(gd2d_var_r4(:,:),nulfile,ni_var, nj_var, nk_var,  &
                        datestamplist(stepIndex),etiket_in,ip1,-1,-1,  &
                        typvar_in,varNameToRead)

          ierr = ezdefset(hco_file%EZscintID,EZscintID_var)
          ierr = int_hInterpScalar( gd2d_file_r4, gd2d_var_r4, &
                                    interpDegree='NEAREST', extrapDegree_opt='NEUTRAL' )

          ! read the corresponding mask if it exists
          if (typvar_var(2:2) == '@') then
            write(*,*) 'gio_readFile: read mask that needs interpolation for variable name: ', nomvar_var
            call utl_abort('gio_readFile: not implemented yet')
          end if

          deallocate(gd2d_var_r4)
        end if

        if (varNameToRead == varName .or. .not. containsFullField) then
          field_r4_ptr(:,:,kIndex,stepIndex) = gd2d_file_r4(1:statevector%hco%ni,1:statevector%hco%nj)
        else
          select case (trim(varName))
          case ('LVIS')
            field_r4_ptr(:,:,kIndex,stepIndex) = &
                 log(max(min(gd2d_file_r4(1:statevector%hco%ni,1:statevector%hco%nj),mpc_maximum_vis_r4),mpc_minimum_vis_r4))
          case ('LPR')
            field_r4_ptr(:,:,kIndex,stepIndex) = &
                 log(mpc_minimum_pr_r4 + max(0.0,gd2d_file_r4(1:statevector%hco%ni,1:statevector%hco%nj)))
          case default
            call utl_abort('gio_readFile: Oups! This should not happen... Check the code.')
          end select
        endif

        if (ierr < 0)then
          write(*,*) varNameToRead, ip1, datestamplist(stepIndex)
          call utl_abort('gio_readFile: Problem with reading file')
        end if

        ! When mpi distribution could put UU on a different mpi task than VV
        ! or only one wind component present in statevector
        ! then we re-read the corresponding UV component and store it
        if (statevector%extraUVallocated) then
          if (varName == 'UU') then
            ierr = fstlir(gd2d_file_r4(:,:),nulfile, ni_file, nj_file, nk_file,  &
                          datestamplist(stepIndex),etiket_in,ip1,-1,-1,  &
                          typvar_in,'VV')
            call gsv_getFieldUV(statevector, gd2d_r4_UV_ptr, kIndex)
            gd2d_r4_UV_ptr(1:gsv_getHco(statevector)%ni,1:gsv_getHco(statevector)%nj, stepIndex) &
                 = gd2d_file_r4(1:gsv_getHco(statevector)%ni,1:gsv_getHco(statevector)%nj)
          else if (varName == 'VV') then
            ierr = fstlir(gd2d_file_r4(:,:),nulfile, ni_file, nj_file, nk_file,  &
                          datestamplist(stepIndex),etiket_in,ip1,-1,-1,  &
                          typvar_in,'UU')
            call gsv_getFieldUV(statevector, gd2d_r4_UV_ptr, kIndex)
            gd2d_r4_UV_ptr(1:gsv_getHco(statevector)%ni,1:gsv_getHco(statevector)%nj, stepIndex) &
                 = gd2d_file_r4(1:gsv_getHco(statevector)%ni,1:gsv_getHco(statevector)%nj)
          end if
        end if

      end do k_loop
    end do

    if (statevector%hco%global .and. statevector%mykCount > 0) then
      write(*,*) 'deallocating hco_file'
      call hco_deallocate(hco_file)
    end if

    ierr = fstfrm(nulfile)
    ierr = fclos(nulfile)        
    if ( associated(gd2d_file_r4) ) deallocate(gd2d_file_r4)

    ! Read in an oceanMask if it is present in the file
    call gio_readMaskFromFile(statevector, trim(filename))

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'gio_readFile: finished'

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
    type(struct_gsv), intent(inout) :: stateVector
    character(len=*), intent(in)    :: fileName

    call ocm_readMaskFromFile(stateVector%oceanMask,gsv_getHco(statevector), &
                              gsv_getVco(statevector), filename)

  end subroutine gio_readMaskFromFile

  !--------------------------------------------------------------------------
  ! gio_getMaskLAM
  !--------------------------------------------------------------------------
  subroutine gio_getMaskLAM(statevector_mask, hco_ptr, vco_ptr, hInterpolateDegree_opt)
    !
    ! :Purpose: To read a LAM mask from a file (./analinc_mask by default).
    !
    implicit none

    ! Arguments
    type(struct_gsv),           intent(inout) :: statevector_mask
    type(struct_hco), pointer,  intent(in)    :: hco_ptr
    type(struct_vco), pointer,  intent(in)    :: vco_ptr
    character(len=*), optional, intent(in)    :: hInterpolateDegree_opt

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
  ! gio_readTrials
  !--------------------------------------------------------------------------
  subroutine gio_readTrials(stateVectorTrialIn)
    !
    ! :Purpose: Reading trials
    !
    implicit none

    ! Arguments
    type(struct_gsv), target, intent(inout) :: stateVectorTrialIn

    ! Locals
    logical             :: fileExists, allocHeightSfc
    logical             :: useInputStateVectorTrial 
    integer, parameter  :: maxNumTrials = 100
    integer             :: fnom, fstouv, fclos, fstfrm, fstinf
    integer             :: ierr, ikey, stepIndex, stepIndexToRead, trialIndex, nulTrial
    integer             :: ni_file, nj_file, nk_file, dateStamp, varNameIndex
    integer             :: procToRead, numBatch, batchIndex, stepIndexBeg, stepIndexEnd

    character(len=2)          :: fileNumber
    character(len=512)        :: fileName
    character(len=4)          :: varNameForDateStampSearch
    character(len=4), pointer :: varNamesToRead(:)

    type(struct_gsv), target  :: stateVectorTrial
    type(struct_gsv), pointer :: stateVectorTrial_ptr 
    type(struct_gsv)          :: stateVector_1step_r4

    call utl_tmg_start(1,'--ReadTrials')

    if ( mmpi_myid == 0 ) then
      write(*,*) ''
      write(*,*) 'gio_readTrials: START'
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
    if ( mmpi_nprocs < stateVectorTrial_ptr%numStep ) then
      write(*,*) 'gio_readTrials: number of trial time steps, mpi tasks = ', stateVectorTrial_ptr%numStep, mmpi_nprocs
      write(*,*) 'gio_readTrials: for better efficiency, the number of mpi tasks should '
      write(*,*) '                be at least as large as number of trial time steps'
    end if

    allocHeightSfc = stateVectorTrial_ptr%heightSfcPresent

    ! figure out number of batches of time steps for reading
    numBatch = ceiling(real(stateVectorTrial_ptr%numStep) / real(mmpi_nprocs))
    write(*,*) 'gio_readTrials: reading will be done by number of batches = ', numBatch

    BATCH: do batchIndex = 1, numBatch

      stepIndexBeg = 1 + (batchIndex - 1) * mmpi_nprocs
      stepIndexEnd = min(stateVectorTrial_ptr%numStep, stepIndexBeg + mmpi_nprocs - 1)
      write(*,*) 'gio_readTrials: batchIndex, stepIndexBeg/End = ', batchIndex, stepIndexBeg, stepIndexEnd

      ! figure out which time step I will read, if any (-1 if none)
      stepIndexToRead = -1
      do stepIndex = stepIndexBeg, stepIndexEnd
        procToRead = nint( real(stepIndex - stepIndexBeg) * real(mmpi_nprocs) / real(stepIndexEnd - stepIndexBeg + 1) )
        if ( procToRead == mmpi_myid ) stepIndexToRead = stepIndex
        if ( mmpi_myid == 0 ) write(*,*) 'gio_readTrials: stepIndex, procToRead = ', stepIndex, procToRead
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

    call utl_tmg_stop(1)

  end subroutine gio_readTrials

  !--------------------------------------------------------------------------
  ! gio_writeToFile
  !--------------------------------------------------------------------------
  subroutine gio_writeToFile(statevector_in, fileName, etiket_in, &
                             scaleFactor_opt, ip3_opt, stepIndex_opt, typvar_opt,&
                             HUcontainsLQ_opt, unitConversion_opt, &
                             writeHeightSfc_opt, numBits_opt, containsFullField_opt)
    !
    ! :Purpose: Write a statevector object to an RPN standard file.
    !
    implicit none

    ! arguments
    type(struct_gsv), target,   intent(in) :: statevector_in
    character(len=*),           intent(in) :: fileName
    character(len=*),           intent(in) :: etiket_in
    real(8),          optional, intent(in) :: scaleFactor_opt
    integer,          optional, intent(in) :: ip3_opt, stepIndex_opt
    character(len=*), optional, intent(in) :: typvar_opt
    logical,          optional, intent(in) :: HUcontainsLQ_opt
    logical,          optional, intent(in) :: unitConversion_opt
    logical,          optional, intent(in) :: writeHeightSfc_opt
    integer,          optional, intent(in) :: numBits_opt
    logical,          optional, intent(in) :: containsFullField_opt

    ! locals
    logical :: iDoWriting, unitConversion, containsFullField

    integer :: fclos, fnom, fstouv, fstfrm, nulnam
    integer :: nulfile, stepIndex
    integer :: ierr, fstecr, ezdefset
    integer :: ni, nj, nk
    integer :: dateo, npak, levIndex, nlev, varIndex, maskLevIndex
    integer :: ip1, ip2, ip3, deet, npas, datyp
    integer :: ig1 ,ig2 ,ig3 ,ig4
    integer :: yourid, nsize, youridy, youridx

    real(4) :: factor_r4, work_r4

    character(len=1)          :: grtyp
    character(len=4)          :: nomvar
    character(len=2)          :: typvar
    character(len=12)         :: etiket
    character(len=4)          :: varLevel
    character(len=4), pointer :: varNamesToRead(:)

    integer, allocatable :: mask(:,:)
    real(4), allocatable :: work2d_r4(:,:), work2dFile_r4(:,:), gd_send_r4(:,:), gd_recv_r4(:,:,:)
    real(8), pointer :: field_r8(:,:,:,:), heightSfc_ptr(:,:)
    real(4), pointer :: field_r4(:,:,:,:)

    type(struct_gsv), pointer :: statevector
    type(struct_gsv), target  :: statevector_tiles

    logical :: interpToPhysicsGrid
    NAMELIST /NAMSTIO/interpToPhysicsGrid

    write(*,*) 'gio_writeToFile: START'

    call utl_tmg_start(161,'low-level--gsv_writeToFile')

    interpToPhysicsGrid = .false.
    if ( .not. utl_isNamelistPresent('NAMSTIO','./flnml') ) then
      if ( mmpi_myid == 0 ) then
        write(*,*) 'gio_writeToFile: namstio is missing in the namelist. The default values will be taken.'
      end if
    else
      ! Read namelist NAMSTIO
      nulnam=0
      ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam,nml=namstio,iostat=ierr)
      if (ierr.ne.0) call utl_abort('gio_writeToFile: Error reading namelist')
      if (mmpi_myid.eq.0) write(*,nml=namstio)
      ierr=fclos(nulnam)
    end if

    !
    !- 1.  Since this routine can only work with 'Tiles' distribution when mpi_local = .true., 
    !      transpose a statevector using 'VarsLevs' distribution
    !
    if ( stateVector_in%mpi_distribution == 'VarsLevs' .and. &
         stateVector_in%mpi_local ) then
      nullify(varNamesToRead)
      call gsv_varNamesList(varNamesToRead,statevector_in) 
      call gsv_allocate(statevector_tiles, statevector_in%numStep, statevector_in%hco, &
                        statevector_in%vco, dataKind_opt=statevector_in%dataKind,      &
                        mpi_local_opt=.true., mpi_distribution_opt='Tiles',            &
                        dateStampList_opt=statevector_in%dateStampList,                &
                        varNames_opt=varNamesToRead)
      call gsv_transposeVarsLevsToTiles(statevector_in, statevector_tiles)
      statevector => stateVector_tiles
      deallocate(varNamesToRead)
    else
      statevector => stateVector_in
    end if

    !
    !- 2.  Set some variables
    !
    if ( .not. statevector%mpi_local ) then
      write(*,*) 'gio_writeToFile: writing statevector that is already mpiglobal!'
    end if

    if (present(ip3_opt)) then
      ip3 = ip3_opt
    else
      ip3 = 0
    end if

    if (present(unitConversion_opt)) then
      unitConversion = unitConversion_opt
    else
      unitConversion = .true.
    end if

    if (present(containsFullField_opt)) then
      containsFullField = containsFullField_opt
    else
      containsFullField = .false.
    end if
    write(*,*) 'gio_writeToFile: containsFullField = ', containsFullField

    ! if step index not specified, choose anltime (usually center of window)
    if (present(stepIndex_opt)) then
      stepIndex = stepIndex_opt
    else
      stepIndex = statevector%anltime
    end if

    if ( present(numBits_opt) ) then
      npak = -numBits_opt
    else
      npak = -32
    end if

    ! initialization of parameters for writing to file
    if (statevector%dateOriginList(stepIndex) /= mpc_missingValue_int) then
      dateo = statevector%dateOriginList(stepIndex)
    else
      dateo  = statevector%dateStampList(stepIndex)
    end if

    if (statevector%deet /= mpc_missingValue_int) then
      deet = statevector%deet
    else
      deet = 0
    end if

    if (statevector%npasList(stepIndex) /= mpc_missingValue_int) then
      npas = statevector%npasList(stepIndex)
    else
      npas = 0
    end if

    if (statevector%ip2List(stepIndex) /= mpc_missingValue_int) then
      ip2 = statevector%ip2List(stepIndex)
    else
      ip2 = 0
    end if

    if (statevector%etiket /= 'UNDEFINED') then
      etiket = statevector%etiket
    else
      etiket = trim(etiket_in)
    end if

    ni     = statevector%ni
    nj     = statevector%nj
    nk     = 1
    if ( present(typvar_opt) ) then
      typvar = trim(typvar_opt)
    else
      typvar = 'R'
    end if
    if ( statevector%oceanMask%maskPresent ) then
      typvar(2:2) = '@'
    end if
    grtyp  = statevector%hco%grtyp
    ig1    = statevector%hco%ig1
    ig2    = statevector%hco%ig2
    ig3    = statevector%hco%ig3
    ig4    = statevector%hco%ig4
    datyp  = 134

    ! only proc 0 does writing or each proc when data is global 
    ! (assuming only called for proc with global data)
    iDoWriting = (mmpi_myid == 0) .or. (.not. statevector%mpi_local)

    !
    !- 3.  Write the global StateVector
    !
    if (iDoWriting) then

      !- Open output field
      nulfile = 0
      write(*,*) 'gio_writeToFile: file name = ',trim(fileName)
      ierr = fnom(nulfile,trim(fileName),'RND+APPEND',0)

      if ( ierr >= 0 ) then
        ierr  =  fstouv(nulfile,'RND')
      else
        call utl_abort('gio_writeToFile: problem opening output file')
      end if

      if (nulfile == 0 ) then
        call utl_abort('gio_writeToFile: unit number for output file not valid')
      end if

      !- Write TicTacToc
      if ( (mmpi_myid == 0 .and. statevector%mpi_local) .or. .not.statevector%mpi_local ) then
        call writeTicTacToc(statevector,nulfile,etiket) ! IN
      endif

    end if

    allocate(gd_send_r4(statevector%lonPerPEmax,statevector%latPerPEmax))
    if ( mmpi_myid == 0 .or. (.not. statevector%mpi_local) ) then
      allocate(work2d_r4(statevector%ni,statevector%nj))
      if (statevector%mpi_local) then
        ! Receive tile data from all mpi tasks
        allocate(gd_recv_r4(statevector%lonPerPEmax,statevector%latPerPEmax,mmpi_nprocs))
      else
        ! Already have entire domain on mpi task (lat/lonPerPEmax == nj/ni)
        allocate(gd_recv_r4(statevector%lonPerPEmax,statevector%latPerPEmax,1))
      end if
    else
      allocate(gd_recv_r4(1,1,1))
      allocate(work2d_r4(1,1))
    end if

    ! Write surface height, if requested
    if ( present(writeHeightSfc_opt) ) then
      if ( writeHeightSfc_opt .and. gsv_isAssocHeightSfc(statevector) ) then
        write(*,*) 'gio_writeToFile: writing surface height'
        heightSfc_ptr => gsv_getHeightSfc(statevector)
        ! MPI communication
        gd_send_r4(1:statevector%lonPerPE, &
                   1:statevector%latPerPE) =  &
             real(heightSfc_ptr(statevector%myLonBeg:statevector%myLonEnd, &
                                    statevector%myLatBeg:statevector%myLatEnd),4)
        if ( (mmpi_nprocs > 1) .and. statevector%mpi_local ) then
          nsize = statevector%lonPerPEmax * statevector%latPerPEmax
          call rpn_comm_gather(gd_send_r4, nsize, 'mpi_real4',  &
                               gd_recv_r4, nsize, 'mpi_real4', 0, 'grid', ierr )
        else
          ! just copy when either nprocs is 1 or data is global
          gd_recv_r4(:,:,1) = gd_send_r4(:,:)
        end if
        if ( mmpi_myid == 0 .and. statevector%mpi_local ) then
          !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
          do youridy = 0, (mmpi_npey-1)
            do youridx = 0, (mmpi_npex-1)
              yourid = youridx + youridy*mmpi_npex
                work2d_r4(statevector%allLonBeg(youridx+1):statevector%allLonEnd(youridx+1),  &
                          statevector%allLatBeg(youridy+1):statevector%allLatEnd(youridy+1)) = &
                  gd_recv_r4(1:statevector%allLonPerPE(youridx+1),  &
                             1:statevector%allLatPerPE(youridy+1),yourid+1)
            end do
          end do
          !$OMP END PARALLEL DO
        else if ( .not. statevector%mpi_local ) then
          work2d_r4(:,:) = gd_recv_r4(:,:,1)
        end if

        ! now do writing
        if (iDoWriting) then
          ip1 = statevector%vco%ip1_sfc
          nomvar = 'GZ'

          !- Scale
          factor_r4 = real(1.0d0/10.0d0,4)
          work2d_r4(:,:) = factor_r4 * work2d_r4(:,:)

          !- Writing to file
          ierr = fstecr(work2d_r4, work_r4, npak, nulfile, dateo, deet, npas, ni, nj, &
                        nk, ip1, ip2, ip3, typvar, nomvar, etiket, grtyp,      &
                        ig1, ig2, ig3, ig4, datyp, .false.)
        end if ! iDoWriting

      end if
    end if

    do varIndex = 1, vnl_numvarmax 
 
      if (gsv_varExist(statevector,vnl_varNameList(varIndex)) ) then

        nlev = statevector%varNumLev(varIndex)

        do levIndex = 1, nlev

          if ( statevector%dataKind == 8 ) then
            call gsv_getField(statevector,field_r8,vnl_varNameList(varIndex))
            gd_send_r4(1:statevector%lonPerPE,  &
                       1:statevector%latPerPE) =  &
                real(field_r8(statevector%myLonBeg:statevector%myLonEnd, &
                              statevector%myLatBeg:statevector%myLatEnd,levIndex,stepIndex),4)
          else
            call gsv_getField(statevector,field_r4,vnl_varNameList(varIndex))
            gd_send_r4(1:statevector%lonPerPE,  &
                       1:statevector%latPerPE) =  &
                field_r4(statevector%myLonBeg:statevector%myLonEnd, &
                         statevector%myLatBeg:statevector%myLatEnd,levIndex,stepIndex)
          end if

          nsize = statevector%lonPerPEmax*statevector%latPerPEmax
          if ( (mmpi_nprocs > 1) .and. (statevector%mpi_local) ) then
            call rpn_comm_gather(gd_send_r4, nsize, 'mpi_real4',  &
                                 gd_recv_r4, nsize, 'mpi_real4', 0, 'grid', ierr )
          else
            ! just copy when either nprocs is 1 or data is global
            gd_recv_r4(:,:,1) = gd_send_r4(:,:)
          end if

          if ( mmpi_myid == 0 .and. statevector%mpi_local ) then
            !$OMP PARALLEL DO PRIVATE(youridy,youridx,yourid)
            do youridy = 0, (mmpi_npey-1)
              do youridx = 0, (mmpi_npex-1)
                yourid = youridx + youridy*mmpi_npex
                work2d_r4(statevector%allLonBeg(youridx+1):statevector%allLonEnd(youridx+1),  &
                            statevector%allLatBeg(youridy+1):statevector%allLatEnd(youridy+1)) = &
                    gd_recv_r4(1:statevector%allLonPerPE(youridx+1),  &
                               1:statevector%allLatPerPE(youridy+1), yourid+1)
              end do
            end do
            !$OMP END PARALLEL DO
          else if ( .not. statevector%mpi_local ) then
            work2d_r4(:,:) = gd_recv_r4(:,:,1)
          end if

          ! now do writing
          if (iDoWriting) then

            ! Set the ip1 value
            if (vnl_varLevelFromVarname(vnl_varNameList(varIndex)) == 'MM') then
              ip1 = statevector%vco%ip1_M(levIndex)
            else if (vnl_varLevelFromVarname(vnl_varNameList(varIndex)) == 'TH') then
              ip1 = statevector%vco%ip1_T(levIndex)
            else if (vnl_varLevelFromVarname(vnl_varNameList(varIndex)) == 'SF') then
              ip1 = 0
            else if (vnl_varLevelFromVarname(vnl_varNameList(varIndex)) == 'SFTH') then
              ip1 = statevector%vco%ip1_T_2m
            else if (vnl_varLevelFromVarname(vnl_varNameList(varIndex)) == 'SFMM') then
              ip1 = statevector%vco%ip1_M_10m
            else if (vnl_varLevelFromVarname(vnl_varNameList(varIndex)) == 'OT') then
              ip1 = vco_ip1_other(levIndex)
            else if (vnl_varLevelFromVarname(vnl_varNameList(varIndex)) == 'DP') then
              ip1 = statevector%vco%ip1_depth(levIndex)
            else if (vnl_varLevelFromVarname(vnl_varNameList(varIndex)) == 'SS') then
              ip1 = statevector%vco%ip1_seaLevel
            else
              varLevel = vnl_varLevelFromVarname(vnl_varNameList(varIndex))
              write(*,*) 'gio_writeToFile: unknown type of vertical level: ', varLevel
              call utl_abort('gio_writeToFile')
            end if

            ! Set the level index for the mask (if present)
            if (vnl_varLevelFromVarname(vnl_varNameList(varIndex)) == 'DP') then
              maskLevIndex = levIndex
            else
              maskLevIndex = 1
            end if
             
            ! Set the output variable name
            nomvar = trim(vnl_varNameList(varIndex))
            if ( trim(nomvar) == 'HU' .and. present(HUcontainsLQ_opt) ) then
               if ( HUcontainsLQ_opt ) nomvar = 'LQ'
            end if

            if ( vnl_varKindFromVarname(trim(nomvar)) == 'CH' .and. containsFullField ) then 
              ! Impose lower limits
              if ( gsv_minValVarKindCH(vnl_varListIndex(nomvar)) > 1.01*mpc_missingValue_r8 ) &
                work2d_r4(:,:) = max( work2d_r4(:,:), real(gsv_minValVarKindCH(vnl_varListIndex(trim(nomvar)))) )
            end if
 
            ! Set the conversion factor
            if ( unitConversion ) then

              if ( trim(nomvar) == 'UU' .or. trim(nomvar) == 'VV') then
                factor_r4 = mpc_knots_per_m_per_s_r4 ! m/s -> knots
              else if ( trim(nomvar) == 'P0' .or. trim(nomvar) == 'UP' .or.  &
                        trim(nomvar) == 'PB' ) then
                factor_r4 = 0.01 ! Pa -> hPa
              else if ( vnl_varKindFromVarname(trim(nomvar)) == 'CH' ) then 
                if ( gsv_conversionVarKindCHtoMicrograms ) then
                  ! Apply inverse transform of unit conversion
                  if ( trim(nomvar) == 'TO3' .or. trim(nomvar) == 'O3L' ) then
                    factor_r4 = 1.0E-9*mpc_molar_mass_dry_air_r4 &
                              /vnl_varMassFromVarName(trim(nomvar)) ! micrograms/kg -> vmr
                  else
                    factor_r4 = 1.0d0 ! no conversion
                  end if
                else
                  factor_r4 = 1.0d0 ! no conversion
                end if
              else
                factor_r4 = 1.0d0 ! no conversion
              end if
            else
              factor_r4 = 1.0
            end if

            if (present(scaleFactor_opt)) factor_r4 = factor_r4 * real(scaleFactor_opt,4)

            !- Scale
            work2d_r4(:,:) = factor_r4 * work2d_r4(:,:)

            !- Convert Kelvin to Celcius only if full field
            if (containsFullField .and. (trim(nomvar) == 'TT' .or. trim(nomvar) == 'TM')  ) then
              where (work2d_r4(:,:) > 100.0)
                work2d_r4(:,:) = work2d_r4(:,:) - mpc_k_c_degree_offset_r4
              end where
            end if

            !- Do interpolation back to physics grid, if needed
            if ( interpToPhysicsGrid .and. statevector%onPhysicsGrid(varIndex) ) then
              write(*,*) 'writeToFile: interpolate this variable back to physics grid: ', &
                         nomvar, associated(statevector%hco_physics)
              allocate(work2dFile_r4(statevector%hco_physics%ni,statevector%hco_physics%nj))
              work2dFile_r4(:,:) = 0.0
              ierr = ezdefset( statevector%hco_physics%EZscintID, statevector%hco%EZscintID )
              ierr = int_hInterpScalar( work2dFile_r4, work2d_r4, &
                                        interpDegree='NEAREST', extrapDegree_opt='NEUTRAL' )

              !- Writing to file
              ierr = fstecr(work2dFile_r4, work_r4, npak, nulfile, dateo, deet, npas, &
                            statevector%hco_physics%ni, statevector%hco_physics%nj, &
                            nk, ip1, ip2, ip3, typvar, nomvar, etiket, grtyp,      &
                            statevector%hco_physics%ig1, statevector%hco_physics%ig2, &
                            statevector%hco_physics%ig3, statevector%hco_physics%ig4, &
                            datyp, .false.)
              deallocate(work2dFile_r4)

            else

              !- Writing to file
              ierr = fstecr(work2d_r4, work_r4, npak, nulfile, dateo, deet, npas, ni, nj, &
                            nk, ip1, ip2, ip3, typvar, nomvar, etiket, grtyp,      &
                            ig1, ig2, ig3, ig4, datyp, .false.)

            end if

            if ( statevector%oceanMask%maskPresent ) then
              if (.not.allocated(mask)) allocate(mask(ni,nj))
              call ocm_copyToInt(statevector%oceanMask,mask,maskLevIndex)
              ierr = fstecr(mask, work_r4, -1, nulfile, dateo, deet, npas, ni, nj, &
                            nk, ip1, ip2, ip3, '@@', nomvar, etiket, grtyp,      &
                            ig1, ig2, ig3, ig4, 2, .false.)
            end if

          end if ! iDoWriting

        end do ! levIndex

      end if ! varExist

    end do ! varIndex

    deallocate(work2d_r4)
    deallocate(gd_send_r4)
    deallocate(gd_recv_r4)
    if (allocated(mask)) deallocate(mask)

    if (iDoWriting) then
      ierr = fstfrm(nulfile)
      ierr = fclos(nulfile)        
    end if

    !
    !- 4.  Ending
    !
    if ( stateVector_in%mpi_distribution == 'VarsLevs' .and. &
         stateVector_in%mpi_local ) then
      call gsv_deallocate(statevector_tiles)
    end if

    call utl_tmg_stop(161)
    write(*,*) 'gio_writeToFile: END'

  end subroutine gio_writeToFile

  !--------------------------------------------------------------------------
  ! writeTicTacToc
  !--------------------------------------------------------------------------
  subroutine writeTicTacToc(statevector,iun,etiket)
    !
    ! :Purpose: Write a statevector object grid descriptors to an RPN standard
    !           file.
    !
    implicit none

    ! Arguments
    type(struct_gsv), intent(in) :: statevector
    integer,          intent(in) :: iun
    character(len=*), intent(in) :: etiket

    ! Locals
    integer :: ier
    integer :: dateo, npak, status, fstecr
    integer :: ip1,ip2,ip3,deet,npas,datyp,ig1,ig2,ig3,ig4
    integer :: ig1_tictac,ig2_tictac,ig3_tictac,ig4_tictac

    character(len=1)  :: grtyp
    character(len=2)  :: typvar

    !
    !- 1.  Writing Tic-Tac
    !
    if ( statevector % hco % grtyp == 'Z' ) then
      npak     = -32
      deet     =  0
      ip1      =  statevector%hco%ig1
      ip2      =  statevector%hco%ig2
      ip3      =  statevector%hco%ig3
      npas     =  0
      datyp    =  1
      grtyp    =  statevector%hco%grtypTicTac
      typvar   = 'X'
      dateo =  0

      call cxgaig ( grtyp,                                          & ! IN
                    ig1_tictac, ig2_tictac, ig3_tictac, ig4_tictac, & ! OUT
                    real(statevector%hco%xlat1), real(statevector%hco%xlon1),   & ! IN
                    real(statevector%hco%xlat2), real(statevector%hco%xlon2)  )   ! IN

      ig1      =  ig1_tictac
      ig2      =  ig2_tictac
      ig3      =  ig3_tictac
      ig4      =  ig4_tictac

      ier = utl_fstecr(statevector%hco%lon*mpc_degrees_per_radian_r8, npak, &
                       iun, dateo, deet, npas, statevector%ni, 1, 1, ip1,    &
                       ip2, ip3, typvar, '>>', etiket, grtyp, ig1,          &
                       ig2, ig3, ig4, datyp, .true.)

      ier = utl_fstecr(statevector%hco%lat*mpc_degrees_per_radian_r8, npak, &
                       iun, dateo, deet, npas, 1, statevector%nj, 1, ip1,    &
                       ip2, ip3, typvar, '^^', etiket, grtyp, ig1,          &
                       ig2, ig3, ig4, datyp, .true.)

      ! Also write the tic tac for the physics grid
      if ( any(statevector%onPhysicsGrid(:)) ) then

        ip1      =  statevector%hco_physics%ig1
        ip2      =  statevector%hco_physics%ig2
        ip3      =  statevector%hco_physics%ig3
        grtyp    =  statevector%hco_physics%grtypTicTac
        
        call cxgaig ( grtyp,                                          & ! IN
                      ig1_tictac, ig2_tictac, ig3_tictac, ig4_tictac, & ! OUT
                      real(statevector%hco_physics%xlat1), real(statevector%hco_physics%xlon1),   & ! IN
                      real(statevector%hco_physics%xlat2), real(statevector%hco_physics%xlon2)  )   ! IN

        ig1      =  ig1_tictac
        ig2      =  ig2_tictac
        ig3      =  ig3_tictac
        ig4      =  ig4_tictac

        ier = utl_fstecr(statevector%hco_physics%lon*mpc_degrees_per_radian_r8, npak, &
                         iun, dateo, deet, npas, statevector%hco_physics%ni, 1, 1, ip1,    &
                         ip2, ip3, typvar, '>>', etiket, grtyp, ig1,          &
                         ig2, ig3, ig4, datyp, .true.)

        ier = utl_fstecr(statevector%hco_physics%lat*mpc_degrees_per_radian_r8, npak, &
                         iun, dateo, deet, npas, 1, statevector%hco_physics%nj, 1, ip1,    &
                         ip2, ip3, typvar, '^^', etiket, grtyp, ig1,          &
                         ig2, ig3, ig4, datyp, .true.)
      end if

    else if ( statevector % hco % grtyp == 'U' ) then
      npak     = -32
      ier = fstecr(statevector%hco%tictacU, statevector%hco%tictacU, npak, iun, 0, 0, 0, size(statevector%hco%tictacU), 1, 1  , &
                   statevector%hco%ig1, statevector%hco%ig2,  statevector%hco%ig3, 'X', '^>', etiket, &
                   'F', 1, 0, 0, 0, 5, .false.)

    else if ( statevector % hco % grtyp == 'Y' ) then
      npak     = -32
      deet     =  0
      ip1      =  statevector%hco%ig1
      ip2      =  statevector%hco%ig2
      ip3      =  statevector%hco%ig3
      npas     =  0
      datyp    =  1
      grtyp    =  statevector%hco%grtypTicTac
      typvar   = 'X'
      dateo =  0

      call cxgaig ( grtyp,                                          & ! IN
                    ig1_tictac, ig2_tictac, ig3_tictac, ig4_tictac, & ! OUT
                    real(statevector%hco%xlat1), real(statevector%hco%xlon1),   & ! IN
                    real(statevector%hco%xlat2), real(statevector%hco%xlon2)  )   ! IN

      ig1      =  ig1_tictac
      ig2      =  ig2_tictac
      ig3      =  ig3_tictac
      ig4      =  ig4_tictac

      ier = utl_fstecr(statevector%hco%lon2d_4*mpc_degrees_per_radian_r8, npak, &
                       iun, dateo, deet, npas, statevector%ni, statevector%nj, 1,    &
                       ip1, ip2, ip3, typvar, '>>', etiket, grtyp,          &
                       ig1, ig2, ig3, ig4, datyp, .true.)

      ier = utl_fstecr(statevector%hco%lat2d_4*mpc_degrees_per_radian_r8, npak, &
                       iun, dateo, deet, npas, statevector%ni, statevector%nj, 1,    &
                       ip1, ip2, ip3, typvar, '^^', etiket, grtyp,          &
                       ig1, ig2, ig3, ig4, datyp, .true.)


    end if

    !
    !- Writing Toc-Toc
    !
    if ( statevector%vco%vgridPresent ) then
      status = vgd_write(statevector%vco%vgrid,iun,'fst')
      if ( status /= VGD_OK ) then
        call utl_abort('writeTicTacToc: ERROR with vgd_write')
      end if
    end if

  end subroutine writeTicTacToc

  !--------------------------------------------------------------------------
  ! gio_fileUnitsToStateUnits
  !--------------------------------------------------------------------------
  subroutine gio_fileUnitsToStateUnits(statevector, containsFullField, stepIndex_opt)
    !
    ! :Purpose: Unit conversion needed after reading RPN standard file
    !
    implicit none

    ! Arguments:
    type(struct_gsv),  intent(inout)  :: statevector
    logical,           intent(in)     :: containsFullField
    integer, optional, intent(in)     :: stepIndex_opt

    ! Locals:
    integer :: stepIndex, stepIndexBeg, stepIndexEnd, kIndex

    real(4), pointer :: field_r4_ptr(:,:,:,:), fieldUV_r4_ptr(:,:,:)
    real(8)          :: multFactor

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

        if ( multFactor /= 1.0d0 ) then
          field_r4_ptr(:,:,kIndex,stepIndex) = real( multFactor * field_r4_ptr(:,:,kIndex,stepIndex), 4 )
        end if

        if ( trim(varName) == 'TT' .and. containsFullField ) then
          field_r4_ptr(:,:,kIndex,stepIndex) = real( field_r4_ptr(:,:,kIndex,stepIndex) +  &
                                                     mpc_k_c_degree_offset_r8, 4 )
        end if

        if ( trim(varName) == 'TM' .and. containsFullField ) then
          where (field_r4_ptr(:,:,kIndex,stepIndex) < 100.0)
            field_r4_ptr(:,:,kIndex,stepIndex) = real( field_r4_ptr(:,:,kIndex,stepIndex) + &
                                                       mpc_k_c_degree_offset_r8, 4 )
          end where
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

        !$OMP PARALLEL DO PRIVATE (kIndex, fieldUV_r4_ptr)
        do kIndex = statevector%myUVkBeg, statevector%myUVkEnd
          nullify(fieldUV_r4_ptr)
          call gsv_getFieldUV(statevector,fieldUV_r4_ptr,kIndex)
          fieldUV_r4_ptr(:,:,stepIndex) =  &
               real( multFactor * fieldUV_r4_ptr(:,:,stepIndex), 4)
        end do
        !$OMP END PARALLEL DO

      end if

    end do step_loop

  end subroutine gio_fileUnitsToStateUnits

end module gridStateVectorFileIO_mod
