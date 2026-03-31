program midas_createLamTemplateGrids
  !
  !:Purpose: Program for computing the LAM analysis grids (core and extended)
  !
  !          ---
  !
  !:Algorithm: Mimics was is done in calcStats for Bhi computation for LAM
  !
  !           --
  !
  !:File I/O: The required input files and produced output files are listed as follows.
  !
  !           --
  !
  !============================================== ====================================================================
  ! Input and Output Files                         Description of file
  !============================================== ====================================================================
  ! ``inputFile``                                  In  - File containing P0 on the desired core analysis grid
  ! ``analysisGrid``                               Out - File containing the needed LAM analysis grid info
  !============================================== ====================================================================
  !
  !           --
  !
  !:Input file format: FST file provide through the -inputFile calling argument
  !
  !:Output file format: FST file named analysisGrid
  !
  !           --
  !
  !:Synopsis: Below is a summary of the ``createLamTemplateGrids`` program calling sequence:
  !
  !           - **Initial setups:**
  !
  !             - Read calling arguements (see below)
  !
  !             - Setup horizontal and vertical grid objects based on the input file.
  !
  !           - **Computation:**
  !
  !             - Adjust the demanded extension zone to ensure that the exteded grid is compatible with fast FFT
  !
  !             - Define and write the analysis grids
  !
  !           --
  !
  !:Options: No namelist is needed for this program. The needed info is provided by 3 program calling arguments
  !          * -inputFile : The FST file containing the P0 on the desired core analysis grid
  !          * -grd_ext_x : The minimum extension in x direction 
  !          * -grd_ext_y : The minimum extension in y direction
  !
  !
 
  use version_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use lamAnalysisGridTransforms_mod
  use lamSpectralTransform_mod
  implicit none

  type(struct_vco), pointer :: vco_inputFile => null()
  type(struct_hco), pointer :: hco_inputFile => null()

  integer            :: fstopc
  integer            :: ierr
  integer            :: ni_ext, nj_ext
  integer            :: ni_ext_adjusted, nj_ext_adjusted
  integer            :: grd_ext_x, grd_ext_y
  integer            :: grd_ext_x_adjusted, grd_ext_y_adjusted

  character(len=256) :: inputFileName

  call ver_printNameAndVersion('createLamTemplateGrids','Create LAM analysis template grids')

  !
  !- 1.  Initilization
  !
  ierr = fstopc('MSGLVL','ERRORS',0)

  ! Read calling arguments
  call get_inputArguments(inputFileName, grd_ext_x, grd_ext_y)

  ! Initialize the horizontal grid
  call hco_SetupFromFile(hco_inputFile, inputFileName, ' ', ' ', varName_opt='P0')

  ! Initialize the vertical grid
  call vco_SetupFromFile(vco_inputFile, inputFileName, ' ')

  !
  !- 2. Create LAM template grids
  !
  ni_ext = hco_inputFile%ni + grd_ext_x
  ni_ext_adjusted = ni_ext
  call lst_nFFT(ni_ext_adjusted)
  nj_ext = hco_inputFile%nj + grd_ext_y
  nj_ext_adjusted = nj_ext
  call lst_nFFT(nj_ext_adjusted)

  grd_ext_x_adjusted = ni_ext_adjusted - hco_inputFile%ni
  grd_ext_y_adjusted = nj_ext_adjusted - hco_inputFile%nj

  write(*,*)
  write(*,'(A,I4,I4)') '           Requested extension zone ', grd_ext_x, grd_ext_y
  write(*,'(A,I4,I4)') ' Fast FFT compatible extension zone ', grd_ext_x_adjusted, grd_ext_y_adjusted
  
  call lgt_createLamTemplateGrids('./analysisgrid', hco_inputFile, vco_inputFile, & ! IN
                                  grd_ext_x_adjusted, grd_ext_y_adjusted)           ! IN

  !
  !- 3.  Ending...
  !
  write(*,*)
  write(*,*) '--------------------------------'
  write(*,*) '> ENDING createLamTemplateGrids '
  write(*,*) '--------------------------------'
  
end program midas_createLamTemplateGrids

!----------------------------------------------------------------------
!- get_inputArguments
!----------------------------------------------------------------------
subroutine get_inputArguments(inputFileName, grd_ext_x, grd_ext_y)
  !
  ! :Purpose: To read program calling arguments
  !
  use utilities_mod
  
  implicit none

  ! Arguments:
  character(len=256), intent(out) :: inputFileName
  integer, intent(out) :: grd_ext_x
  integer, intent(out) :: grd_ext_y

  ! Locals:
  integer, parameter :: nArguments = 3

  character(len=24)  :: arg(nArguments)
  character(len=256) :: argDefaultValue1(nArguments)
  character(len=256) :: argDefaultValue2(nArguments) 

  !- Program argument...
  !  names
  data arg /'input.','grd_ext_x','grd_ext_y'/
  !  values if argument not added in program call
  data argDefaultValue1 /'none','-1','-1'/
  !  values if program called with argument but without assigned value
  data argDefaultValue2 /'none','-1','-1'/

  !- Linkages
  call ccard(arg, argDefaultValue2, argDefaultValue1, nArguments, -1)

  read (argDefaultValue1(1),fmt='(a256)') inputFileName
  read (argDefaultValue1(2),*) grd_ext_x
  read (argDefaultValue1(3),*) grd_ext_y

  !- Error traps
  if (inputFileName == 'none') then
    call utl_abort('get_inputArguments: -input missing or empty')
  end if
  if (grd_ext_x < 0) then
    call utl_abort('get_inputArguments: -grd_ext_x missing, empty or negative')
  end if
  if (grd_ext_y < 0) then
    call utl_abort('get_inputArguments: -grd_ext_y missing, empty or negative')
  end if

end subroutine get_inputArguments
