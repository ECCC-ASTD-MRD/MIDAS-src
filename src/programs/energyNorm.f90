program midas_energyNorm
  !
  !:Purpose: Main program for computing the energyNorm of an atmospheric state
  !
  !          ---
  !
  !:Algorithm: The energy norm definition is defined in the paper `HFSOI approach <https://doi.org/10.1175/MWR-D-17-0252.1>`_
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
  ! ``flnml``                                      In  - Main namelist file with parameters user may modify
  ! ``inputFiles``                                 In  - File containing a list of files which are the input atmopheric state
  ! ``energyNorm_ascii``                           Out - File containing the energy norms for each input file
  !============================================== ====================================================================
  !
  !
  !           --
  !
  !:Input file format: Here is a description of the input file format ``inputFiles``
  !
  !                    Each line is a path to an input atmospheric
  !                    state.
  !
  !                    The first file specified is considered as the
  !                    reference state
  !
  !                    The energy norm of each other files will be
  !                    computed using that reference state.
  !
  !                    A line starting with '#' or '!' will be ignored
  !                    and considered as a comment.
  !
  !:Output file format: Here is a description of the output file format ``energyNorm_ascii``
  !
  !                     The first line is the path to the reference
  !                     state.
  !
  !                     The second line is the multiplicative factor
  !                     applied to the energy norms.
  !
  !                     The third line is an empty line.
  !
  !                     The fourth line is the legend of the following
  !                     lines.
  !
  !                     Each other line starts with the file name and
  !                     the energy norm total and followed with the
  !                     energy of each component.
  !
  !           --
  !
  !:Synopsis: Below is a summary of the ``energyNorm`` program calling sequence:
  !
  !           - **Initial setups:**
  !
  !             - Parse ``inputFiles`` to get the reference state and
  !                  the other states to compute the energy norm using
  !                  the reference state.
  !
  !             - Setup horizontal and vertical grid objects the reference state.
  !
  !             - Allocate the ``stateVector`` objects on the input grid
  !
  !             - Read all the input files with ``gio_readFromFile``
  !
  !           - **Computation:**
  !
  !             - For each atmospheric state:
  !
  !                - Compute the difference between the atmospheric state and the reference state with ``gsv_add``
  !
  !                - Call ``gvt_energyNorm``: Compute the energy norm for that atmospheric state
  !
  !                - Write the energy norm in the file ``energyNorm_ascii``
  !
  !           - **Final steps:**
  !
  !             - Deallocate the arrays and state vectors
  !
  !           --
  !
  !:Options: `List of namelist blocks <../namelists_in_each_program.html#energynorm>`_
  !          that can affect the ``energyNorm`` program.
  !
  !          * The namelist blocks used to configure the program are listed in the following table:
  !
  !
  !========================= ====================== =============================================================
  ! Module                   Namelist               Description of what is controlled
  !========================= ====================== =============================================================
  ! ``midas_energy``         ``NAMENERGYNORM``      The variable ``fullStates`` if input states are full or not.
  !                                                 You can apply a scaling to energy norm values with ``multiplicativeFactor``.
  ! ``timeCoord_mod``        ``NAMTIME``            assimilation time window length, temporal resolution of
  !                                                 the background state and increment
  ! ``gridStateVector_mod``  ``NAMSTATE``           set the variables to read
  !========================= ====================== =============================================================
  !
  !
  use rmn_fnom
  use midasMpi_mod
  use version_mod
  use utilities_mod
  use runtimeInfo_mod
  use message_mod
  use timeCoord_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use gridVariableTransforms_mod

  implicit none

  character(len=256), parameter :: inputFileName  = 'inputFiles'
  character(len=256), parameter :: outputFileName = 'energyNorm_ascii'
  integer :: ierr,nulFileOutput
  integer :: fileIndex, numberOfFiles, maxFileLength
  character(len=1024) :: referenceFileName, fileNameHeader, fileNameFormat
  character(len=1024), allocatable :: fileNames(:)
  type(struct_gsv) :: stateVector, stateVectorReference
  type(struct_vco), pointer :: vco => null()
  type(struct_hco), pointer :: hco => null()

  ! Namelist variables
  logical :: fullStates ! If '.true.', then the files will be considered as full states and the energy norm will be compute with the difference of the state and the reference state (default is ``.true``).
  real(8) :: multiplicativeFactor ! Multiplicative factor to apply to energy norm outputs (default is 1)

  namelist /namEnergyNorm/ fullStates, multiplicativeFactor

  call ver_printNameAndVersion('energyNorm','Compute the energy norm of an atmospheric state')

  ! MPI initilization
  call mmpi_initialize

  call tmg_init(mmpi_myid, 'TMG_INFO')

  call rti_tmg_start(0,'Main')
  call rti_printTime()

  ! Read the namelists
  call utl_readNml()

  !- Read the namelist for energyNorm program (if it exists)
  ! set default values for the namelist variables
  fullStates = .true.
  multiplicativeFactor = 1.0D0
  if (utl_isNamelistPresent('namEnergyNorm', './flnml')) then
    call rti_tmg_start(181,'low-level--readNML')
    read(utl_flnml, nml = namEnergyNorm, iostat = ierr)
    if (ierr /= 0) call rti_abort('midas-energyNorm: Error reading namelist namEnergyNorm')
    if (mmpi_myid == 0) write(*,nml = namEnergyNorm)
    call rti_tmg_stop(181)
  else
    write(*,*)
    write(*,*) 'midas-energyNorm: Namelist block namEnergyNorm is missing in the namelist.'
    write(*,*) '                    The default value will be taken.'
    if (mmpi_myid == 0) write(*, nml = namEnergyNorm)
  end if

  !- Initialize the Temporal grid and set dateStamp from env variable
  call tim_setup()

  ! Setup the format of the output RPN standard files to 'XDF' or 'RSF'
  call gio_setup

  !- Initialize variables of the model states
  call gsv_setup
  call msg_memUsage('midas-energyNorm')

  ! parse the input file to get the files names
  call findFileNames(inputFileName, numberOfFiles, referenceFileName, fileNames, maxFileLength)

  write(*,*) 'midas-energyNorm: Opening ascii output file: ', trim(outputFileName)
  nulFileOutput = 0
  if ( mmpi_myid == 0 ) then
    ierr = fnom(nulFileOutput, trim(outputFileName), 'SEQ+R/W', 0)
    if (ierr /= 0) then
      call rti_abort('midas-energyNorm: Cannot open ascii output file')
    end if
  end if

  call initializeReferenceState(referenceFileName, stateVectorReference, hco, vco)

  call mmpi_barrier

  if ( mmpi_myid == 0 ) then
    ! write the reference file in the output
    write(nulFileOutput,'(a,a)') 'The reference file is ', trim(referenceFileName)
    write(nulFileOutput,'(a,ES14.4)') 'multiplicative factor = ', multiplicativeFactor
    write(nulFileOutput,*)

    ! Here, we add 4 to the length for the string to be written.  This
    ! number comes from the fact that the numeric format in
    ! 'computeEnergyNorm' is 'ES14.4' which is leaving 4 blank spaces.
    fileNameFormat = 'a' // str(maxFileLength+4)

    ! write header in file
    fileNameHeader = 'fileName' ! We put this string in a character array to have 'fileName' left justified
    write(nulFileOutput,'(' // trim(fileNameFormat) // ',6a14)') fileNameHeader,    &
                                                                  'total         ', &
                                                                  'tt            ', &
                                                                  'uu            ', &
                                                                  'vv            ', &
                                                                  'hu            ', &
                                                                  'p0            '
  end if

  call gsv_allocate(stateVector, tim_nstepobs, hco, vco,  &
                    dateStamp_opt=tim_getDateStamp(),     &
                    mpi_local_opt=.true., dataKind_opt=8, &
                    hInterpolateDegree_opt='LINEAR',      &
                    beSilent_opt=.false. )

  fileNameFormat = 'a' // str(maxFileLength)
  do fileIndex = 1, numberOfFiles
    call computeEnergyNorm(stateVector, stateVectorReference, fileNames(fileIndex),       &
                           fullStates, multiplicativeFactor, nulFileOutput, fileNameFormat)
    call mmpi_barrier
  end do ! fileIndex

  ! closing 'outputFileName'
  if ( mmpi_myid == 0 ) then
    ierr = fclos(nulFileOutput)
  end if

  call gsv_deallocate(stateVector)
  call gsv_deallocate(stateVectorReference)
  deallocate(fileNames)

  !
  !- 3. Job termination
  !

  call rti_printTime()
  call rti_tmg_stop(0)

  call tmg_terminate(mmpi_myid, 'TMG_INFO')

  call mmpi_finalize

contains

  !--------------------------------------------------------------------------
  ! findFileNames
  !--------------------------------------------------------------------------
  subroutine findFileNames(inputFileName, numberOfFilesToProcess, referenceFileName, &
                           fileNames, maxFileLength)
    !
    ! :Purpose: Helper function which parses the input file to extract
    !           the input files names
    !
    implicit none

    ! Arguments:
    character(len=*),              intent(in)  :: inputFileName ! input file name in which all the atmospheric state vectors are defined
    integer,                       intent(out) :: numberOfFilesToProcess ! number of files described in 'inputFiles' (excluding the reference file)
    character(len=*),              intent(out) :: referenceFileName ! reference file name with which the energy norm will be computed
    character(len=*), allocatable, intent(out) :: fileNames(:) ! list of file names for which the energy norms will be computed
    integer,                       intent(out) :: maxFileLength ! maximum length of the input file names

    ! Locals:
    integer :: ierr, nulFileInput

    write(*,*) 'midas-energyNorm: findFileNames: Opening file: ', trim(inputFileName)
    nulFileInput = 0
    ierr = fnom(nulFileInput, trim(inputFileName), 'SEQ+R/O', 0)
    if (ierr /= 0) then
      call rti_abort('midas-energyNorm: findFileNames: Cannot open ascii output file')
    end if

    ! Read a first time the file 'inputFileName' to find the number of files to
    ! process ('numberOfFilesToProcess')
    call parseInputFiles(nulFileInput, numberOfFilesToProcess)

    ! Allocate the array 'fileNames'
    allocate(fileNames(numberOfFilesToProcess))

    ! We will read again each lines so we rewind the file unit
    ! position to the beginning of the file.
    rewind(nulFileInput)

    ! Read a second time the file to set 'referenceFileName' and fill
    ! the array 'fileNames' with the file names to process
    call parseInputFiles(nulFileInput, referenceFileName_opt = referenceFileName, &
                         fileNames_opt = fileNames, maxFileLength_opt = maxFileLength)

    ! closing 'inputFileName'
    ierr = fclos(nulFileInput)

  end subroutine findFileNames

  !--------------------------------------------------------------------------
  ! parseInputFiles
  !--------------------------------------------------------------------------
  subroutine parseInputFiles(nulFile, numberOfFilesToProcess_opt, &
                             referenceFileName_opt, fileNames_opt, maxFileLength_opt)
    !
    ! :Purpose: Helper function which parses the input file line by line to extract
    !           the number of files or the input files names themselves
    !
    implicit none

    ! Arguments:
    integer,                       intent(in)              :: nulFile ! file unit of the 'inputFiles'
    integer,                       intent(out),   optional :: numberOfFilesToProcess_opt ! number of files described in 'inputFiles' (excluding the reference file)
    character(len=*),              intent(out),   optional :: referenceFileName_opt ! reference file name with which the energy norm will be computed
    integer,                       intent(out),   optional :: maxFileLength_opt ! maximum length of the input file names
    character(len=*), allocatable, intent(inout), optional :: fileNames_opt(:) ! list of file names for which the energy norms will be computed

    ! Locals:
    integer :: readStatus, lineNumber, charIndex
    character(len=1024) :: line, trimmedLine
    integer :: numberOfInputFiles
    logical :: initializeAllFileNames, anyFileNamesOptArgIsGiven

    initializeAllFileNames = present(referenceFileName_opt) .and. present(maxFileLength_opt) .and. present(fileNames_opt)
    anyFileNamesOptArgIsGiven = present(referenceFileName_opt) .or. present(maxFileLength_opt) .or. present(fileNames_opt)

    if (.not. initializeAllFileNames .and. anyFileNamesOptArgIsGiven ) then
      call rti_abort('midas-energyNorm: parseInputFiles has been called with one or two of ''referenceFileName_opt'', ''fileNames_opt'', ''maxFileLength_opt''.  All must be specified or none.')
    end if

    if ( initializeAllFileNames ) then
      maxFileLength_opt = 0
      if (.not.allocated(fileNames_opt)) then
        call rti_abort('midas-energyNorm: parseInputFiles has been called without a proper allocated ''fileNames_opt'' argument.')
      end if
    end if

    lineNumber = 0
    numberOfInputFiles = 0
    readLoopLineByLine: do
      lineNumber = lineNumber + 1
      read(nulFile, '(a)', iostat=readStatus) line

      ! If we reached the end of file, exit the loop
      if ( readStatus < 0 ) exit readLoopLineByLine
      ! We encountered an error while reading the file
      if ( readStatus > 0 ) then
        call rti_abort('midas-energyNorm: parseInputFiles: Problem reading line ' // str(lineNumber) // ' of file ' // trim(inputFileName))
      end if

      ! build 'trimmedLine' by removing leading spaces in 'line'
      ! the function 'trim' is only removing the trailing spaces in 'line'
      trimmedLine = ''
      trimLoop: do charIndex = 1, len_trim(line)
        if ( line(charIndex:charIndex) /= ' ' ) then
          trimmedLine = line(charIndex:len_trim(line))
          exit trimLoop
        end if
      end do trimLoop

      ! If the line starts with '#' or '!' or is empty, we ignore the line.
      if ( trimmedLine(1:1) == '#' .or. trimmedLine(1:1) == '!' .or. len_trim(trimmedLine) == 0 ) then
        cycle readLoopLineByLine
      end if

      numberOfInputFiles = numberOfInputFiles + 1

      if ( initializeAllFileNames ) then
        if ( numberOfInputFiles == 1 ) then
          referenceFileName_opt = trimmedLine
        else
          maxFileLength_opt = max(len_trim(trimmedLine), maxFileLength_opt)
          fileNames_opt(numberOfInputFiles-1) = trimmedLine
        end if
      end if
    end do readLoopLineByLine

    if ( numberOfInputFiles == 0 ) then
      call rti_abort('midas-energyNorm: parseInputFiles: No input state has been given in the ''' // trim(inputFileName) // '''')
    else if ( numberOfInputFiles == 1 ) then
      call rti_abort('midas-energyNorm: parseInputFiles: No state has been given in the ''' // trim(inputFileName) // ''' other than the reference state')
    end if

    if (present(numberOfFilesToProcess_opt)) then
      ! numberOfInputFiles includes the reference file which
      ! numberOfFilesToProcess does not include
      numberOfFilesToProcess_opt = numberOfInputFiles-1
    end if

  end subroutine parseInputFiles

  !--------------------------------------------------------------------------
  ! initializeReferenceState
  !--------------------------------------------------------------------------
  subroutine initializeReferenceState(inputFileName, stateVector, hco, vco)
    !
    ! :Purpose: Helper function which initializes the reference state
    !           vector and horizonal and vertical definitions from
    !           'inputileName'
    !
    implicit none

    ! Arguments:
    character(len=*),          intent(in)  :: inputFileName ! The path to the file read the reference state
    type(struct_gsv),          intent(out) :: stateVector ! gridStateVector structure to store the reference state
    type(struct_hco), pointer, intent(out) :: hco ! 'horizontalCoord' which will be used for all gridStateVector reads after
    type(struct_vco), pointer, intent(out) :: vco ! 'verticalCoord' which will be used for all gridStateVector reads after

    call hco_SetupFromFile(hco, inputFileName, etiketName=' ')
    call vco_SetupFromFile(vco, inputFileName)

    call gsv_allocate(stateVector, tim_nstepobs, hco, vco, &
                      dateStamp_opt=tim_getDateStamp(),             &
                      mpi_local_opt=.true., dataKind_opt=8,         &
                      hInterpolateDegree_opt='LINEAR',              &
                      beSilent_opt=.false.)

    call rti_tmg_start(1,'--ReadingStateVectorRef')
    call gio_readFromFile(stateVector, inputFileName, etiket_in=' ', typvar_in=' ')
    call rti_tmg_stop(1)
    call msg_memUsage('midas-energyNorm')

  end subroutine initializeReferenceState

  !--------------------------------------------------------------------------
  ! computeEnergyNorm
  !--------------------------------------------------------------------------
  subroutine computeEnergyNorm(stateVector, stateVectorReference, fileName, &
                               fullState, multiplicativeFactor, nulFile,    &
                               fileNameFormat)
    !
    ! :Purpose: Helper function which computes the energy norm for an
    !           input file with respect to a reference. It then prints
    !           the result in the outputFile identified by 'nulFile'.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: stateVector ! state vector structure filled by the content of 'fileName'
    type(struct_gsv), intent(in)    :: stateVectorReference ! state vector representing the reference state
    character(len=*), intent(in)    :: fileName ! input file to read the state vector on which the energy norm will be computed
    logical,          intent(in)    :: fullState ! Does the state vector in 'fileName' representing a full state or not
    real(8),          intent(in)    :: multiplicativeFactor ! factor applied to energy norm values to scale them
    integer,          intent(in)    :: nulFile ! output file unit in which the energy norm will be written
    character(len=*), intent(in)    :: fileNameFormat ! format when writing the file name in the output file

    ! Constants:
    ! The energy norm will be computed on the globe
    real(8), parameter :: latMin = -95.0d0
    real(8), parameter :: latMax =  95.0d0
    real(8), parameter :: lonMin = -185.0d0
    real(8), parameter :: lonMax =  365.0d0
    ! The energy norm will include all those components
    logical, parameter :: includeUVNorm = .true.
    logical, parameter :: includeTTNorm = .true.
    logical, parameter :: includeP0Norm = .true.
    logical, parameter :: includeHUNorm = .true.
    logical, parameter :: includeTGNorm = .false.
    ! If 'straNorm' is .true. then the error norm is from 100hPa to 1hPa.
    ! If .false., it is from surface to 100hPa.
    ! The energy norm is computed for the troposphere.
    logical, parameter :: straNorm = .false.

    ! Locals:
    type(struct_gvt_energyNorm) :: energyNorm

    call msg_memUsage('midas-energyNorm')

    call rti_tmg_start(2,'--ReadingStateVector')
    call gio_readFromFile(stateVector, fileName, etiket_in=' ', typvar_in=' ')
    call rti_tmg_stop(2)
    call msg_memUsage('midas-energyNorm')

    ! If 'fullState' is true then we must compute the difference
    ! between the reference state and the state itself.
    ! If 'fullState' is false, then we don't need it because it is
    ! already a difference (for example, a standard deviation or RMS).
    if ( fullState ) then
      ! compute the difference between the state vector and the reference
      ! stateVector = stateVector - stateVectorReference
      call rti_tmg_start(3,'--computeStateVectorDifference')
      call gsv_add(stateVectorReference, stateVector, -1.0d0)
      call rti_tmg_stop(3)
      call msg_memUsage('midas-energyNorm')
    end if

    call rti_tmg_start(4,'--computeEnergyNorm')
    energyNorm = gvt_energyNorm(stateVector, stateVectorReference, &
                                latMin=latMin, latMax=latMax, lonMin=lonMin, lonMax=lonMax, &
                                uvNorm=includeUVNorm, ttNorm=includeTTNorm, p0Norm=includeP0Norm, &
                                huNorm=includeHUNorm, tgNorm=includeTGNorm, straNorm=straNorm)
    call rti_tmg_stop(4)

    call msg_memUsage('midas-energyNorm')

    !! Write the result in the file
    if ( mmpi_myid == 0 ) then
      write(nulFile,'(' // trim(fileNameFormat) // ',6ES14.4)') fileName, multiplicativeFactor*energyNorm%total, &
                                                                          multiplicativeFactor*energyNorm%tt,    &
                                                                          multiplicativeFactor*energyNorm%uu,    &
                                                                          multiplicativeFactor*energyNorm%vv,    &
                                                                          multiplicativeFactor*energyNorm%hu,    &
                                                                          multiplicativeFactor*energyNorm%p0
    end if

  end subroutine computeEnergyNorm

end program midas_energyNorm
