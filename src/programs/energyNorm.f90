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

  character(len=256) :: outFileName
  integer            :: istamp,exdb,exfin,fnom,fclos,ierr,nulFile,energyNormIndex
  type(struct_gsv)   :: stateVector, stateVectorReference
  type(struct_vco), pointer :: vco => null()
  type(struct_hco), pointer :: hco => null()
  type(struct_gvt_energyNorm), allocatable :: energyNorms(:)
  integer, parameter :: numberOfInputFiles = 1

  istamp = exdb('ENERGYNORM','DEBUT','NON')

  call ver_printNameAndVersion('energyNorm','Compute the energy norm of an atmospheric state')

  ! MPI initilization
  call mmpi_initialize

  call tmg_init(mmpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')
  call utl_printTime()

  ! Read the namelists
  call utl_readNml()

  !
  !- Initialize the Temporal grid and set dateStamp from env variable
  !
  call tim_setup()

  ! Setup the format of the output RPN standard files to 'XDF' or 'RSF'
  call gio_setup

  !
  !- Initialize variables of the model states
  !
  call gsv_setup
  call msg_memUsage('midas-energyNorm')

  !
  !- Initialize the Analysis grid
  !
  call hco_SetupFromFile(hco, 'inputFile', etiketName=' ')

  !
  !- Initialisation of the analysis grid vertical coordinate from analysisgrid file
  call vco_SetupFromFile(vco, 'inputFile')

  allocate(energyNorms(numberOfInputFiles))

  call gsv_allocate(stateVector, tim_nstepobs, hco, vco,  &
                    dateStamp_opt=tim_getDateStamp(), &
                    mpi_local_opt=.true., dataKind_opt=8,  &
                    hInterpolateDegree_opt='LINEAR', &
                    beSilent_opt=.false. )

  call gsv_allocate(stateVectorReference, tim_nstepobs, hco, vco,  &
                    dateStamp_opt=tim_getDateStamp(), &
                    mpi_local_opt=.true., dataKind_opt=8,  &
                    hInterpolateDegree_opt='LINEAR', &
                    beSilent_opt=.false. )

  call utl_tmg_start(1,'--ReadingStateVector')
  call gio_readFromFile(stateVector, 'inputFile', etiket_in=' ', typvar_in=' ')
  call utl_tmg_stop(1)

  call gsv_getInfo(stateVector, 'Reading from ''inputFile''')

  call utl_tmg_start(2,'--ReadingStateVectorRef')
  call gio_readFromFile(stateVectorReference, 'inputFile_ref', etiket_in=' ', typvar_in=' ')
  call utl_tmg_stop(2)
  call msg_memUsage('midas-energyNorm')

  call gsv_getInfo(stateVectorReference, 'Reading from ''inputFile_ref''')

  ! Compute the energy norm with state vector and its reference
  energyNorms(1) = compute_energyNorm(stateVector, stateVectorReference)

  call utl_tmg_start(5,'--WriteEnergyNormsToAscii')

  ! outFileName = trim(outFileName) // '_ascii'
  outFileName = 'energyNorm_ascii'
  write(*,*) 'epp_printRmsStats: Opening ascii output file: ', trim(outFileName)
  nulFile = 0
  ierr = fnom (nulFile, outFileName, 'SEQ+R/W', 0)
  if (ierr /= 0) then
    call utl_abort('epp_printRmsStats: Cannot open ascii output file')
  end if

  !! Write the energy norm value in the file 'outFileName'
  do energyNormIndex = 1, numberOfInputFiles
    write(nulFile,*) 'energy norm = ', energyNorms(energyNormIndex)
  end do

  ierr = fclos(nulFile)

  call utl_tmg_stop(5)

  call gsv_deallocate(stateVector)
  call gsv_deallocate(stateVectorReference)
  deallocate(energyNorms)

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
  function compute_energyNorm(stateVector, stateVectorReference) result(energyNorm)
    !
    ! :Purpose: Helper function which computes the energy norms by
    !           taking the difference with the reference and calling
    !           'gvt_energyNorm' on a stateVector
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: stateVector
    type(struct_gsv), intent(in)    :: stateVectorReference
    ! Result:
    type(struct_gvt_energyNorm) :: energyNorm

    ! Constants
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

  end function compute_energyNorm

end program midas_energyNorm
