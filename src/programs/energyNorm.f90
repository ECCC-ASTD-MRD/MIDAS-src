program midas_energyNorm
  !
  !:Purpose: Main program for computing the energyNorm of an atmospheric state
  !
  !          ---
  !
  !:Algorithm: The energy norm definition is defined in the paper `HFSOI approach <https://doi.org/10.1175/MWR-D-17-0252.1>`_
  !
  !            --
  !
  !:File I/O: The required input files and produced output files are listed as follows.
  !
  !           --
  !
  !============================================== ====================================================================
  ! Input and Output Files                         Description of file
  !============================================== ====================================================================
  ! ``flnml``                                      In - Main namelist file with parameters user may modify
  ! ``inputFile``                                  In - File containing the atmospheric state with which the energy norm will be computed
  ! ``inputFile_ref``                              In - File containing the reference state
  !============================================== ====================================================================
  !
  !           --
  !
  !:Synopsis: Below is a summary of the ``energyNorm`` program calling sequence:
  !
  !           - **Initial setups:**
  !
  !             - Setup horizontal and vertical grid objects from ``inputFile``.
  !
  !             - Allocate the stateVector objects on the input grid
  !
  !             - Read the files 'inputFile_ref' and 'inputFile' with 'gio_readFromFile'
  !
  !           - **Computation:**
  !
  !             - Compute the difference between the two objects with 'gsv_add'
  !
  !             - ``gvt_energyNorm``: Compute the energy norm
  !  !
  !           - **Final steps:**
  !
  !             - Write the energy norm value in the file ``energyNorm_ascii``
  !
  !             - Deallocate the state vectors
  !
  !           --
  !
  !:Options: `List of namelist blocks <../namelists_in_each_program.html#energyNorm>`_
  !          that can affect the ``energyNorm`` program.
  !
  !          * The namelist blocks used to configure FSOI are listed in the following table:
  !
  !
  !========================= ====================== =============================================================
  ! Module                   Namelist               Description of what is controlled
  !========================= ====================== =============================================================
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
  integer :: istamp,exdb,exfin,fnom,fclos,ierr,nulFileInput,nulFileOutput
  integer :: readStatus, lineNumber, charIndex
  type(struct_gsv), target  :: stateVector, stateVectorReference
  type(struct_vco), pointer :: vco => null()
  type(struct_hco), pointer :: hco => null()
  character(len=1024) :: line, trimmedLine
  logical :: isReferenceStateInitialized = .false.
  logical :: isAtLeastOneEnergyNormHasBeenComputed = .false.

  istamp = exdb('ENERGYNORM','DEBUT','NON')

  call ver_printNameAndVersion('energyNorm','Compute the energy norm of an atmospheric state')

  ! MPI initilization
  call mmpi_initialize

  call tmg_init(mmpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')
  call utl_printTime()

  ! Read the namelists
  call utl_readNml()

  !- Initialize the Temporal grid and set dateStamp from env variable
  call tim_setup()

  ! Setup the format of the output RPN standard files to 'XDF' or 'RSF'
  call gio_setup

  !- Initialize variables of the model states
  call gsv_setup
  call msg_memUsage('midas-energyNorm')

  write(*,*) 'midas-energyNorm: Opening file: ', trim(inputFileName)
  nulFileInput = 0
  ierr = fnom(nulFileInput, trim(inputFileName), 'SEQ+R/O', 0)
  if (ierr /= 0) then
    call utl_abort('midas-energyNorm: Cannot open ascii output file')
  end if

  write(*,*) 'midas-energyNorm: Opening ascii output file: ', trim(outputFileName)
  nulFileOutput = 0
  if ( mmpi_myid == 0 ) then
    ierr = fnom(nulFileOutput, trim(outputFileName), 'SEQ+R/W', 0)
    if (ierr /= 0) then
      call utl_abort('midas-energyNorm: Cannot open ascii output file')
    end if
  end if

  lineNumber = 0
  readLoop: do
    lineNumber = lineNumber + 1
    read(nulFileInput, '(a)', iostat=readStatus) line

    ! If we reached the end of file, exit the loop
    if ( readStatus < 0 ) exit readLoop
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
      cycle readLoop
    end if

    if ( .not. isReferenceStateInitialized ) then
      call initializeReferenceState(trimmedLine, stateVectorReference, hco, vco)

      if ( mmpi_myid == 0 ) then
        ! write the reference file in the output
        write(nulFileOutput,'(a,a)') 'The reference file is ', trim(trimmedLine)
        write(nulFileOutput,*)
        ! write header in file
        write(nulFileOutput,'("fileName",6a14)') 'total', 'tt', 'uu', 'vv', 'hu', 'p0'
      end if

      call gsv_allocate(stateVector, tim_nstepobs, hco, vco,  &
                        dateStamp_opt=tim_getDateStamp(),     &
                        mpi_local_opt=.true., dataKind_opt=8, &
                        hInterpolateDegree_opt='LINEAR',      &
                        beSilent_opt=.false. )

      isReferenceStateInitialized = .true.
      call msg_memUsage('midas-energyNorm')
    else
      ! Compute the energy norm with state vector and its reference
      ! and print in the output file given by 'line'
      call compute_energyNorm(stateVectorReference, trimmedLine, stateVector, nulFileOutput)
      isAtLeastOneEnergyNormHasBeenComputed = .true.
    endif ! .not. isReferenceStateInitialized

    call utl_printTime()
  end do readLoop

  call rpn_comm_barrier('GRID',ierr)

  if ( .not. isReferenceStateInitialized ) then
    write(nulFileOutput,'(a)') 'No input state has been given in the ''' // trim(inputFileName) // ''''
  else if ( .not. isAtLeastOneEnergyNormHasBeenComputed ) then
    write(*,*) 'No state has been given in the ''' // trim(inputFileName) // ''' other than the reference state'
    if ( mmpi_myid == 0 ) then
      write(nulFileOutput,*) 'No state has been given in the ''' // trim(inputFileName) // ''' other than the reference state'
    end if
  end if

  ! closing 'inputFileName'
  ierr = fclos(nulFileInput)
  ! closing 'outputFileName'
  if ( mmpi_myid == 0 ) then
    ierr = fclos(nulFileOutput)
  end if

  if ( isReferenceStateInitialized ) then
    call gsv_deallocate(stateVector)
    call gsv_deallocate(stateVectorReference)
  end if

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
  ! compute_energyNorm
  !--------------------------------------------------------------------------
  subroutine initializeReferenceState(inputFileName, stateVectorReference, hco, vco)
    !
    ! :Purpose: Helper function initializes the reference state vector
    !           and horizonal and vertical definitions from
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
  ! energyNorm
  !--------------------------------------------------------------------------
  subroutine compute_energyNorm(stateVectorReference, fileName, stateVector, nulFile)
    !
    ! :Purpose: Helper function which computes the energy norm for an
    !           input file with respect to a reference. If then prints
    !           the result in the outputFile identified by 'nulFile'.
    !
    implicit none

    ! Arguments:
    type(struct_gsv), pointer, intent(in)  :: stateVectorReference
    character(len=*), intent(in)  :: fileName
    type(struct_gsv), pointer, intent(in) :: stateVector
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

    ! compute the difference between the state vector and the reference
    ! stateVector = stateVector - stateVectorReference
    call utl_tmg_start(3,'--computeStateVectorDifference')
    call gsv_add(stateVectorReference, stateVector, -1.0d0)
    call utl_tmg_stop(3)

    call msg_memUsage('midas-energyNorm')

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
