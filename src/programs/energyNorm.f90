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
  !                     The second line is the legend of the following
  !                     lines
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
  !             - Allocate the stateVector objects on the input grid
  !
  !             - Read all the input files with 'gio_readFromFile'
  !
  !           - **Computation:**
  !
  !             - For each atmospheric state:
  !
  !                - Compute the difference between the atmospheric state and the reference state with 'gsv_add'
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
  ! ``timeCoord_mod``        ``NAMTIME``            assimilation time window length, temporal resolution of
  !                                                 the background state and increment
  ! ``gridStateVector_mod``  ``NAMSTATE``           set the variables to read
  !========================= ====================== =============================================================
  !
  !
  use version_mod
  use utilities_mod
  use midasMpi_mod
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
  integer :: istamp,exdb,exfin,fnom,fclos,ierr,nulFileOutput
  integer :: fileIndex, numberOfFiles
  character(len=1024) :: referenceFileName
  character(len=1024), allocatable :: fileNames(:)
  type(struct_gsv), target  :: stateVector, stateVectorReference
  type(struct_vco), pointer :: vco => null()
  type(struct_hco), pointer :: hco => null()

  ! Namelist variables
  logical :: fullStates

  namelist /namEnergyNorm/ fullStates

  istamp = exdb('ENERGYNORM','DEBUT','NON')

  call ver_printNameAndVersion('energyNorm','Compute the energy norm of an atmospheric state')

  ! MPI initilization
  call mmpi_initialize

  call tmg_init(mmpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')
  call utl_printTime()

  ! Read the namelists
  call utl_readNml()

  !- Read the namelist for energyNorm program (if it exists)
  ! set default values for the namelist variable
  fullStates = .true.
  if (utl_isNamelistPresent('namEnergyNorm', './flnml')) then
    call utl_tmg_start(181,'low-level--readNML')
    read(utl_flnml, nml = namEnergyNorm, iostat = ierr)
    if (ierr /= 0) call utl_abort('midas-energyNorm: Error reading namelist namEnergyNorm')
    if (mmpi_myid == 0) write(*,nml = namEnergyNorm)
    call utl_tmg_stop(181)
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
  call parseInputFiles(inputFileName,referenceFileName,fileNames,numberOfFiles)

  write(*,*) 'midas-energyNorm: Opening ascii output file: ', trim(outputFileName)
  nulFileOutput = 0
  if ( mmpi_myid == 0 ) then
    ierr = fnom(nulFileOutput, trim(outputFileName), 'SEQ+R/W', 0)
    if (ierr /= 0) then
      call utl_abort('midas-energyNorm: Cannot open ascii output file')
    end if
  end if

  call initializeReferenceState(referenceFileName, stateVectorReference, hco, vco)

  call rpn_comm_barrier('GRID',ierr)

  if ( mmpi_myid == 0 ) then
    ! write the reference file in the output
    write(nulFileOutput,'(a,a)') 'The reference file is ', trim(referenceFileName)
    write(nulFileOutput,*)
    ! write header in file
    write(nulFileOutput,'("fileName",6a14)') 'total', 'tt', 'uu', 'vv', 'hu', 'p0'
  end if

  call gsv_allocate(stateVector, tim_nstepobs, hco, vco,  &
                    dateStamp_opt=tim_getDateStamp(),     &
                    mpi_local_opt=.true., dataKind_opt=8, &
                    hInterpolateDegree_opt='LINEAR',      &
                    beSilent_opt=.false. )

  do fileIndex = 1, numberOfFiles
    call compute_energyNorm(stateVectorReference, fileNames(fileIndex), stateVector, &
                            fullStates, nulFileOutput)
    call rpn_comm_barrier('GRID',ierr)
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
  istamp = exfin('ENERGYNORM','FIN','NON')

  call utl_printTime()
  call utl_tmg_stop(0)

  call tmg_terminate(mmpi_myid, 'TMG_INFO')

  call rpn_comm_finalize(ierr)

contains

  !--------------------------------------------------------------------------
  ! parseInputFiles
  !--------------------------------------------------------------------------
  subroutine parseInputFiles(inputFileName, referenceFileName, fileNames, &
                             numberOfFilesToComputeTheEnergyNormAgainstTheReference)
    !
    ! :Purpose: Helper function which parses the input file to extract
    !           the input files names
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in)  :: inputFileName
    character(len=*), intent(out) :: referenceFileName
    character(len=*), allocatable, intent(out) :: fileNames(:)
    integer,          intent(out) :: numberOfFilesToComputeTheEnergyNormAgainstTheReference

    ! Locals:
    integer :: ierr,readStatus, lineNumber, charIndex
    character(len=1024) :: line, trimmedLine
    integer :: nulFileInput, numberOfInputFiles
    logical :: readFileIsFirstPass

    write(*,*) 'midas-energyNorm: Opening file: ', trim(inputFileName)
    nulFileInput = 0
    ierr = fnom(nulFileInput, trim(inputFileName), 'SEQ+R/O', 0)
    if (ierr /= 0) then
      call utl_abort('midas-energyNorm: Cannot open ascii output file')
    end if

    ! Read a first time the file to find the number of files: readFileIsFirstPass = .true.
    !     This allows to allocate the array 'fileNames'
    ! Read a second time the file to set 'referenceFileName' and fill the array 'fileNames' : readFileIsFirstPass = .false.
    readFileIsFirstPass = .true.
    readFileIterationLoop: do
      lineNumber = 0
      numberOfInputFiles = 0
      readLoopLineByLine: do
        lineNumber = lineNumber + 1
        read(nulFileInput, '(a)', iostat=readStatus) line

        ! If we reached the end of file, exit the loop
        if ( readStatus < 0 ) exit readLoopLineByLine
        ! We encountered an error while reading the file
        if ( readStatus > 0 ) then
          call utl_abort('midas-energyNorm: Problem reading line ' // str(lineNumber) // ' of file ' // trim(inputFileName))
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
          ! write(*,*) 'The line ' // str(lineNumber) // ' is a comment or is empty'
          cycle readLoopLineByLine
        end if

        numberOfInputFiles = numberOfInputFiles + 1

        if ( .not. readFileIsFirstPass ) then
          if ( numberOfInputFiles == 1 ) then
            referenceFileName = trimmedLine
          else
            fileNames(numberOfInputFiles-1) = trimmedLine
          end if
        end if
      end do readLoopLineByLine

      if ( readFileIsFirstPass ) then
        if ( numberOfInputFiles == 0 ) then
          call utl_abort('No input state has been given in the ''' // trim(inputFileName) // '''')
        else if ( numberOfInputFiles == 1 ) then
          call utl_abort('No state has been given in the ''' // trim(inputFileName) // ''' other than the reference state')
        end if

        numberOfFilesToComputeTheEnergyNormAgainstTheReference = numberOfInputFiles-1
        allocate(fileNames(numberOfFilesToComputeTheEnergyNormAgainstTheReference))

        rewind(nulFileInput)
        readFileIsFirstPass = .false.
      else
        exit readFileIterationLoop
      end if
    end do readFileIterationLoop

    ! closing 'inputFileName'
    ierr = fclos(nulFileInput)

  end subroutine parseInputFiles

  !--------------------------------------------------------------------------
  ! initializeReferenceState
  !--------------------------------------------------------------------------
  subroutine initializeReferenceState(inputFileName, stateVectorReference, hco, vco)
    !
    ! :Purpose: Helper function which initializes the reference state
    !           vector and horizonal and vertical definitions from
    !           'inputileName'
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: inputFileName
    type(struct_gsv), pointer, intent(in)  :: stateVectorReference
    type(struct_vco), pointer, intent(out) :: vco
    type(struct_hco), pointer, intent(out) :: hco

    call hco_SetupFromFile(hco, inputFileName, etiketName=' ')
    call vco_SetupFromFile(vco, inputFileName)

    call gsv_allocate(stateVectorReference, tim_nstepobs, hco, vco, &
                      dateStamp_opt=tim_getDateStamp(),             &
                      mpi_local_opt=.true., dataKind_opt=8,         &
                      hInterpolateDegree_opt='LINEAR',              &
                      beSilent_opt=.false.)

    call utl_tmg_start(1,'--ReadingStateVectorRef')
    call gio_readFromFile(stateVectorReference, inputFileName, etiket_in=' ', typvar_in=' ')
    call utl_tmg_stop(1)
    call msg_memUsage('midas-energyNorm')

  end subroutine initializeReferenceState

  !--------------------------------------------------------------------------
  ! compute_energyNorm
  !--------------------------------------------------------------------------
  subroutine compute_energyNorm(stateVectorReference, fileName, stateVector, fullState, nulFile)
    !
    ! :Purpose: Helper function which computes the energy norm for an
    !           input file with respect to a reference. It then prints
    !           the result in the outputFile identified by 'nulFile'.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), pointer, intent(in)  :: stateVectorReference
    character(len=*), intent(in)  :: fileName
    type(struct_gsv), pointer, intent(in) :: stateVector
    logical :: fullState
    integer :: nulFile

    ! Constants:
    real(8), parameter :: latMin = -95.0d0
    real(8), parameter :: latMax =  95.0d0
    real(8), parameter :: lonMin = -185.0d0
    real(8), parameter :: lonMax =  365.0d0
    logical, parameter :: includeUVNorm = .true.
    logical, parameter :: includeTTNorm = .true.
    logical, parameter :: includeP0Norm = .true.
    logical, parameter :: includeHUNorm = .true.
    logical, parameter :: includeTGNorm = .false.
    ! if 'straNorm' is .true. then the error norm is from 100hPa to 1hPa,
    ! if .false., it is from surface to 100hPa
    logical, parameter :: straNorm = .false.

    ! Locals:
    type(struct_gvt_energyNorm) :: energyNorm

    call msg_memUsage('midas-energyNorm')

    call utl_tmg_start(2,'--ReadingStateVector')
    call gio_readFromFile(stateVector, fileName, etiket_in=' ', typvar_in=' ')
    call utl_tmg_stop(2)
    call msg_memUsage('midas-energyNorm')

    !! If 'fullState' is true then we must compute the difference
    !! between the reference state and the state itself.
    !! If 'fullState' is false, then we don't need it because it is
    !! already a difference (for example, a standard deviation or RMS).
    if ( fullState ) then
      ! compute the difference between the state vector and the reference
      ! stateVector = stateVector - stateVectorReference
      call utl_tmg_start(3,'--computeStateVectorDifference')
      call gsv_add(stateVectorReference, stateVector, -1.0d0)
      call utl_tmg_stop(3)
      call msg_memUsage('midas-energyNorm')
    end if

    call utl_tmg_start(4,'--computeEnergyNorm')
    energyNorm = gvt_energyNorm(stateVector, stateVectorReference, &
                                latMin=latMin, latMax=latMax, lonMin=lonMin, lonMax=lonMax, &
                                uvNorm=includeUVNorm, ttNorm=includeTTNorm, p0Norm=includeP0Norm, &
                                huNorm=includeHUNorm, tgNorm=includeTGNorm, straNorm=straNorm)
    call utl_tmg_stop(4)

    call msg_memUsage('midas-energyNorm')

    !! Write the result in the file
    if ( mmpi_myid == 0 ) then
      write(nulFile,'(a,6ES14.4)') trim(fileName), energyNorm%total, energyNorm%tt, energyNorm%uu, energyNorm%vv, energyNorm%hu, energyNorm%p0
    end if

  end subroutine compute_energyNorm

end program midas_energyNorm
