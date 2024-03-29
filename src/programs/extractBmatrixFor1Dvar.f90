program midas_extractBmatrixFor1Dvar
  !
  !:Purpose: Main program to extract B matrix to binary file Bmatrix.bin suitable for 1Dvar applications
  !          The B matrix is defined at a set of locations specified by the variable lonlatExtract and date 
  !          extractDate
  !
  !          ---
  !
  !:Algorithm:  The B matrix is computed column by column by application of operators bmat_sqrtBT and bmat_sqrtB 
  !
  !            --
  !
  !:File I/O: The required input files and produced output files are listed as follows.
  !
  !============================================== ==============================================================
  ! Input and Output Files                          Description of file
  !============================================== ==============================================================
  ! ``flnml``                                      In - Main namelist file with parameters user may modify
  ! ``flnml_static``                               In - The "static" namelist that should not be modified
  ! ``bgcov``                                      In - The B NMC matrix file
  ! ``analysisgrid``                               In - analysis grid 
  ! ``Bmatrix.bin``                                Out - The Bmatrix binary file for 1DVar
  !============================================== ==============================================================
  !
  !           --
  !
  !:Synopsis: Below is a summary of the ``extractBmatrixFor1Dvar`` program calling sequence:
  !    
  !             - **Initial setups:**
  !
  !               - Read parameters from the program namelist section NAMEXTRACT
  !
  !               - Setup temporal grid
  !
  !               - Get vertical and horizontal grid information from analysisgrid file
  !
  !               - Allocate a gridStatevector object
  !
  !               - Initialize the B matrix
  !
  !               - Initialize the gridded variable transform module 
  !
  !             - **Computation:**
  !
  !               - For each location specified in the namelist: 
  !
  !                   - for each element of the stateVector
  !
  !                      - set this elemnt to one and the other to zero
  !
  !                      - apply bmat_sqrtBT
  !
  !                      - apply bmat_sqrtB
  !
  !                      - get the corresponding column of the B matrix
  !
  !                   - write the resulting B matrix to file
  !           --
  !:Options: `List of namelist blocks <../namelists_in_each_program.html#extractbmatrixfor1dvar>`_
  !          that can affect the ``extractBmatrixFor1Dvar`` program.
  !
  !          - The use of ``extractBmatrixFor1Dvar`` program is controlled by the namelist block
  !            ``&namextract`` read by the ``extractBmatrixFor1Dvar`` program.
  !              *  ``extractdate`` date (YYYYMMDDHH format) for the B matrix extracted
  !              *  ``lonlatExtract`` (longitudes, latitudes) pairs definining the locations where the B matrix is to be extracted
  !              *  ``varNameExtract`` name of the variable to extract or ``all`` to extract everything in namstate
  !              *  ``stepBinExtract`` should be one of ``first``, ``middle`` or ``last`` to define when in the assimilation window the B matrix is valid
  !
  use version_mod
  use midasMpi_mod
  use mathPhysConstants_mod
  use controlVector_mod
  use gridVariableTransforms_mod
  use varNameList_mod
  use gridStateVector_mod
  use bMatrix_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use utilities_mod
  use ramDisk_mod
  implicit none

  type(struct_gsv) :: statevector
  type(struct_hco), pointer :: hco_anl  => null()
  type(struct_hco), pointer :: hco_core => null()
  type(struct_vco), pointer :: vco_anl  => null()
  real(4) :: factor1, factor2
  real(8), pointer :: field4d(:,:,:,:)
  real(8), allocatable :: controlVector(:)
  integer, parameter :: nmaxLevs = 100
  real(4) :: latitude, longitude
  integer, external :: fclos, fnom, fstopc, newdate, get_max_rss
  integer :: ierr
  integer :: varIndex, nkgdim, levIndex1, lonIndex, latIndex, levIndex2
  integer :: kIndex1, kIndex2, columnProcIdLocal, columnProcIdGlobal, nulmat, varCount
  integer :: idate, itime, nulnam, dateStamp
  integer :: stepBinExtractIndex
  integer :: nLonLatPos, lonLatPosIndex

  integer, parameter :: lonColumn = 1
  integer, parameter :: latColumn = 2

  character(len=10)  :: datestr
  character(len=4)   :: varName1, varName2
  character(len=4),allocatable :: varList(:)
  real(8), allocatable :: Bmatrix(:,:)

  ! namelist variables
  integer :: extractdate               ! date for the B matrix extracted
  integer :: lonlatExtract(nmaxLevs,2) ! lon lat pairs definining the locations where the B matrix is to be extracted
  character(len=128) :: stepBinExtract ! number of step bins to extract (1 typically for B NMC)
  character(len=4)   :: varNameExtract ! variable to extract (all to extract everything in namstate)
  namelist /namextract/ extractdate, lonlatExtract, &
                    varNameExtract, stepBinExtract

  call ver_printNameAndVersion('extractBmatrix', 'Extract 1Dvar B matrix')

  ! MPI, tmg initialization
  call mmpi_initialize
  call tmg_init(mmpi_myid, 'TMG_INFO')
  call utl_tmg_start(0, 'Main')
  ierr = fstopc('MSGLVL', 'ERRORS', 0)

  ! Set default values for namelist NAMEXTRACT parameters
  extractdate        =  2011020100
  varNameExtract     = 'all'
  stepBinExtract     = 'middle'
  lonlatExtract(:,:) = -1

  ! Read the parameters from NAMEXTRACT
  nulnam = 0
  ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
  read(nulnam, nml=namextract, iostat=ierr)
  if (ierr /= 0) call utl_abort('midas-extractBmatrix: Error reading namelist')
  write(*, nml=namextract)
  ierr = fclos(nulnam)

  nLonLatPos = 0
  do lonlatPosIndex = 1, size(lonlatExtract(:,lonColumn))
    if (lonlatExtract(lonlatPosIndex,lonColumn) >= 1 .and. lonlatExtract(lonlatPosIndex,latColumn) >= 1) nLonLatPos = nLonLatPos + 1  
  end do

  if ( nLonLatPos == 0 ) then
    call utl_abort('midas-extractBmatrix: You should specify at least one lonlat position !')
  end if

  ! Decompose extractdate(yyyymmddhh) into idate(YYYYMMDD) itime(HHMMSShh)
  ! and calculate date-time stamp
  idate = extractdate/100
  itime = (extractdate-idate*100)*1000000
  ierr = newdate(dateStamp, idate, itime, 3)
  write(datestr,'(i10.10)') extractdate
  write(*,*)' idate= ',idate,' time= ',itime
  write(*,*)' date= ',extractdate,' stamp= ',dateStamp

  ! Top Level Control setup
  call ram_setup
  !- Initialize the Temporal grid
  call tim_setup
  if (tim_getDateStamp() == 0) then
    call tim_setDatestamp(dateStamp)
  else
    call utl_abort('midas-extractBmatrix: Code changes needed to use dateStamp from env variable')
  end if
  ! Initialize variables of the model states
  call gsv_setup
  ! Initialize the Analysis horizontal grid
  call hco_SetupFromFile( hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN
  if ( hco_anl % global ) then
    hco_core => hco_anl
  else
    !- Iniatilized the core (Non-Extended) analysis grid
    call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore' ) ! IN
  end if
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  ! Initialize the vertical coordinate from the statistics file
  call vco_SetupFromFile( vco_anl,        & ! OUT
                          './analysisgrid') ! IN
  ! Allocate the statevector
  call gsv_allocate(statevector, tim_nstepobsinc, hco_anl, vco_anl, &
                    datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true., &
                    allocHeight_opt=.false., allocPressure_opt=.false.)
  call gsv_zero(statevector)
  nkgdim = statevector%nk
  allocate( Bmatrix(nkgdim, nkgdim) )
  ! Setup the B matrix
  call bmat_setup(hco_anl, hco_core, vco_anl)
  !- Initialize the gridded variable transform module
  call gvt_setup(hco_anl, hco_core, vco_anl)
  call gvt_setupRefFromTrialFiles('HU')
  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  !
  !==============================================
  !- Compute columns of B matrix
  !==============================================
  !
  select case(trim(stepBinExtract))
  case ('first')
   stepBinExtractIndex = 1
  case ('middle')
    if (mod(tim_nstepobsinc,2) == 0) then
      write(*,*)
      write(*,*) 'odd number of nstepobsinc a required for obs place in the middle of the analysis window '
      write(*,*) 'tim_nstepobsinc = ', tim_nstepobsinc
      call utl_abort('midas-extractBmatrix')
    end if
   stepBinExtractIndex = (tim_nstepobsinc+1)/2
  case ('last')
   stepBinExtractIndex = tim_nstepobsinc
  case default
    write(*,*)
    write(*,*) 'Unsupported stepBinExtract : ', trim(stepBinExtract)
    call utl_abort('midas-extractBmatrix')
  end select
  
  allocate(controlVector(cvm_nvadim))

  write(*,*) '************************************************************'
  write(*,*) 'midas-extractBmatrix: Compute columns of B matrix for 1Dvar'
  write(*,*) '************************************************************'
  write(*,*)
  write(*,*) ' temporal location           = ', trim(stepBinExtract), stepBinExtractIndex
  write(*,*) ' number of lon-lat positions = ', nLonLatPos
  
  if (mmpi_myId == 0) then
    varCount = 0
    do varIndex=1, vnl_numvarmax
      if ( .not. gsv_varExist(varName=vnl_varNameList(varIndex)) ) cycle
      if ( trim(varNameExtract)  /= 'all' .and. (trim(varNameExtract) /= trim(vnl_varNameList(varIndex))) ) cycle
      varCount = varCount + 1
    end do
    allocate(varList(varCount))
    varCount = 0
    do varIndex=1, vnl_numvarmax
      if ( .not. gsv_varExist(varName=vnl_varNameList(varIndex)) ) cycle
      if ( trim(varNameExtract) /= 'all' .and. (trim(varNameExtract) /= trim(vnl_varNameList(varIndex))) ) cycle
      varCount = varCount + 1
      varList(varCount) = vnl_varNameList(varIndex)
    end do
    nulmat = 0
    ierr = fnom(nulmat, './Bmatrix.bin', 'FTN+SEQ+UNF', 0)
    write(nulmat) extractDate, vco_anl % nlev_T, vco_anl % nlev_M, vco_anl % Vcode, &
         vco_anl % ip1_sfc, vco_anl % ip1_T_2m, vco_anl % ip1_M_10m, varCount, nkgdim, nLonLatPos
    write(nulmat) vco_anl % ip1_T(:), vco_anl % ip1_M(:), varList(:)
  end if
  
  locationLoop:do lonLatPosIndex = 1, nLonLatPos

    Bmatrix(:,:) = MPC_missingValue_R8

    latIndex = lonlatExtract(lonLatPosIndex,latColumn)
    lonIndex = lonlatExtract(lonLatPosIndex,lonColumn)
    
    latitude = hco_anl%lat2d_4(lonIndex, latIndex)
    longitude = hco_anl%lon2d_4(lonIndex, latIndex)

    variableLoop1:do kIndex1 = 1, nkgdim
      varName1 = gsv_getVarNameFromK(statevector, kIndex1)
      if ( .not. gsv_varExist(varName=varName1) ) cycle
      if ( trim(varNameExtract) /= 'all' .and. trim(varNameExtract) /= trim(varName1) ) cycle
      
      write(*,*)
      write(*,*) 'midas-extractBmatrix: simulating a pseudo-observation of ', trim(varName1)
      
      factor1 = getConversionFactor( varName1 )
      levIndex1 = gsv_getLevFromK(statevector, kIndex1)
      call gsv_zero(statevector)
      call gsv_getField(statevector,field4d, varName1)
      if ( latIndex >= statevector%myLatBeg .and. latIndex <= statevector%myLatEnd .and. &
           lonIndex >= statevector%myLonBeg .and. lonIndex <= statevector%myLonEnd ) then
        if (vnl_varLevelFromVarname(varName1) == 'SF') then
          field4d(lonIndex, latIndex, 1, stepBinExtractIndex) = 1.0D0
        else
          field4d(lonIndex, latIndex, levIndex1, stepBinExtractIndex) = 1.0D0
        end if
      end if
      !ici field4d est initialise a zero sauf au point qui nous interesse ou on a 1
      !et donc aussi statevector (puisque field4V est un pointeur vers la partie correspondante de statevector)
      controlVector(:) = 0.0d0
      call bmat_sqrtBT(controlVector, cvm_nvadim, statevector)
      call bmat_sqrtB (controlVector, cvm_nvadim, statevector)
      ! ici statevector contient la colonne correspondante de la matrice B
      write(*,*) 'midas-extractBmatrix: writing out the column of B. levIndex1,lonIndex,latIndex=', &
                 levIndex1,lonIndex,latIndex
      
      variableLoop2:do kIndex2 = 1, nkgdim
        varName2 = gsv_getVarNameFromK(statevector, kIndex2)
        if ( .not. gsv_varExist(varName= varName2) ) cycle
        if ( trim(varNameExtract) /= 'all' .and. trim(varNameExtract) /= trim(varName2) ) cycle
        columnProcIdLocal = -1
        if ( latIndex >= statevector % myLatBeg .and. latIndex <= statevector % myLatEnd .and. &
             lonIndex >= statevector % myLonBeg .and. lonIndex <= statevector % myLonEnd ) then
          columnProcIdLocal = mmpi_myId
          call gsv_getField(statevector, field4d, varName2)
          factor2 = getConversionFactor( varName2 )
          levIndex2 = gsv_getLevFromK(statevector, kIndex2)
          bmatrix(kIndex2, kIndex1) = factor1 * factor2 * field4d(lonIndex, latIndex, levIndex2, stepBinExtractIndex)
        end if
        call rpn_comm_allreduce(columnProcIdLocal, columnProcIdGlobal, 1, "mpi_integer", "mpi_max", "GRID", ierr)
      end do variableLoop2

    end do variableLoop1

    call RPN_COMM_bcast(Bmatrix, nkgdim * nkgdim, 'MPI_REAL8', columnProcIdGlobal, 'GRID', ierr )
    if (mmpi_myId ==0) then
      write(nulmat) latitude, longitude,  Bmatrix(:,:)
    end if

  end do locationLoop
    
  ierr = fclos(nulmat)

  deallocate(Bmatrix)
  deallocate(controlVector)

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

  call gsv_deallocate(statevector)

  write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
  
  ! MPI, tmg finalize
  call utl_tmg_stop(0)
  call tmg_terminate(mmpi_myid, 'TMG_INFO')
  call rpn_comm_finalize(ierr) 

  write(*,*) ' --------------------------------'
  write(*,*) ' midas-extractBmatrix ENDS'
  write(*,*) ' --------------------------------'

contains

  real(4) function getConversionFactor(varName)
    character(len=*), intent(in) :: varName
    if ( trim(varName) == 'UU' .or. trim(varName) == 'VV') then
      getConversionFactor = mpc_knots_per_m_per_s_r4 ! m/s -> knots
    else if ( trim(varName) == 'P0' ) then
      getConversionFactor = 0.01 ! Pa -> hPa
    else 
      getConversionFactor = 1.0 ! no conversion
    end if
  end function getConversionFactor

end program midas_extractBmatrixFor1Dvar
