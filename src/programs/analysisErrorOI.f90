!-------------------------------------- LICENCE BEGIN ------------------------------------
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

program midas_analysisErrorOI
  !
  !:Purpose: Calculate analysis-error standard deviation given
  !          new assimilated observations. It only works for sea ice variables and
  !          uses a simple OI approach.
  !
  !          ---
  !
  !:Algorithm: The Optimal Interpolation (OI) is a data assimilation approach
  !            where both the state and its estimated error are computed using
  !            observations while taking into account the specified 
  !            uncertainties for both the observations and the background
  !            state (i.e. the R and B covariance matrices, respectively).
  !            The computations are done independently at each analysis grid
  !            point with only local observations, those that have a 
  !            significant influence on the analysis at the grid point
  !            location. Here the code only implements the calculation to
  !            update the diagonal of the B covariance matrix.
  !
  !            --
  !
  !============================================== ==============================================================
  ! Input and Output Files                         Description of file
  !============================================== ==============================================================
  ! ``flnml``                                      In - Main namelist file with parameters user may modify
  ! ``trlm_01``                                    In - Background-error standard deviation
  ! ``analysisgrid``                               In - File defining grid for computing the analysis error
  ! ``sea_ice_obs-err``                            In - Observation error statistics
  ! ``bgstddev``                                   In - Static background-error statistics
  ! ``bgSeaIceConc``                               In - Background sea ice concentration
  ! ``obsfiles_$FAM/obs$FAM_0001_0001``            In - Observation file for each "family" (only 1 MPI task)
  ! ``anlm_000m``                                  Out - Analysis-error on the analysis grid
  ! ``obsfiles_$FAM.updated/obs$FAM_0001_0001``    Out - Updated obs file for each "family" (only 1 MPI task)
  !============================================== ==============================================================
  !
  !           --
  !
  !:Synopsis: Below is a summary of the ``analysisErrorOI`` program calling sequence:
  !
  !             - **Initial setups:**
  !
  !               - Setup horizontal and vertical grid objects for "analysis
  !                 grid" from ``trlm_01`` file.
  !
  !               - Setup ``obsSpaceData`` object and read observations from
  !                 files: ``inn_setupObs``.
  !
  !               - Setup ``columnData`` and ``gridStateVector`` modules (read
  !                 list of analysis variables from namelist) and allocate column
  !                 object for storing trial on analysis levels.
  !
  !               - Setup the observation error statistics in ``obsSpaceData``
  !                 object: ``oer_setObsErrors``.
  !
  !               - Filter out observations from satellites
  !                 not specified in the name list: ``filt_iceConcentration``.
  !
  !               - Filter scatterometer backscatter anisotropy observations
  !                 where wind speed is too small: ``filt_backScatAnisIce``.
  !
  !               - Setup observation-error for scatterometer backscatter
  !                 anisotropy observations: ``oer_setErrBackScatAnisIce``.
  !
  !             - **Main calculation:**
  !
  !               - Compute the analysis-error: ``aer_analysisError``.
  !
  !               - Update the Days Since Last Obs: ``aer_daysSinceLastObs``.
  !
  !               - Update the observation files: ``obsf_writeFiles``.
  !
  !           --
  !
  !:Options: `List of namelist blocks <../namelists_in_each_program.html#analysiserroroi>`_
  !          that can affect the ``analysisErrorOI`` program.
  !
  !          * The relevant namelist blocks used to configure the
  !            analysis-error calculation are listed in the following table:
  ! 
  !======================== ================== ==============================================================
  ! Module                   Namelist           Description of what is controlled
  !======================== ================== ==============================================================
  ! ``timeCoord_mod``       ``NAMTIME``         assimilation time window length, temporal resolution of
  !                                             the background state
  ! ``columndata_mod``      ``NAMSTATE``        name of the analysis variable (RPN nomvar, 4-character long),
  !                                             only sea ice concentration (GL) is allowed for now
  ! ``gridstatevector_mod``       "                                   "
  ! ``obsspacedata_mod``    ``NAMDIMO``         specify the maximum number of header and body elements
  ! ``obsfilter_mod``       ``NAMFILT``         list of varno to use and bit flags (13-bit#) for filtering
  ! ``obsfilter_mod``       ``namPlatformIce``  list of observation platforms to assimilate
  ! ``sqliteread_mod``      ``NAMSQLgl``        list of varno to read from SQLite files for sea ice
  ! ``sqliteread_mod``      ``namSQLUpdate``    list of elements to update in SQLite files
  ! ``sqliteread_mod``      ``namSQLInsert``    place holder, could be empty
  ! ``analysisErrorOI_mod`` ``NAMAER``          set the maximum analysis-error std dev. allowed
  !======================== ================== ==============================================================
  !
  use version_mod
  use ramDisk_mod
  use utilities_mod
  use midasMpi_mod
  use message_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use obsSpaceData_mod
  use columnData_mod
  use gridStateVector_mod
  use obsFiles_mod
  use innovation_mod
  use obsErrors_mod
  use obsFilter_mod  
  use analysisErrorOI_mod

  implicit none

  integer :: istamp, exdb, exfin
  integer :: ierr, dateStampFromObs
  integer :: get_max_rss
  character(len=48) :: obsMpiStrategy, varMode
  character(len=20) :: trlmFileName
  character(len=15), parameter :: myName = 'analysisErrorOI'

  type(struct_obs)       , target :: obsSpaceData
  type(struct_columnData), target :: trlColumnOnAnlLev
  type(struct_hco)      , pointer :: hco_anl => null()
  type(struct_vco)      , pointer :: vco_anl => null()

  istamp = exdb('ANALYSISERROROI','DEBUT','NON')

  call ver_printNameAndVersion(myName,'Program to calculate the analysis-error standard deviation for sea ice using OI.')

  ! MPI initialization
  call mmpi_initialize

  if( mmpi_nprocs > 1 ) then
    call msg(myName,'mmpi_nprocs = '//str(mmpi_nprocs))
    call utl_abort(myName//': this version of the code should only be used with one mpi task.')
  end if

  call tmg_init(mmpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')

  varMode='analysis'

  ! Setup the ram disk
  call ram_setup

  obsMpiStrategy = 'LIKESPLITFILES'

  !
  !- Initialize the Temporal grid and set dateStamp from env variable
  !
  call tim_setup()

  if (tim_nstepobs > 1 .or. tim_nstepobsinc > 1) then
    call utl_abort(myName//': The program assumes only one time step.')
  end if

  !
  !- Initialize observation file names and set datestamp if not already
  !
  call obsf_setup(dateStampFromObs, varMode)
  if (tim_getDateStamp() == 0) then
    if (dateStampFromObs > 0) then
      call tim_setDateStamp(dateStampFromObs)
    else
      call utl_abort(myName//': DateStamp was not set')
    end if
  end if

  !
  !- Initialize constants
  !
  if (mmpi_myid == 0) call mpc_printConstants(6)

  !
  !- Initialize list of analyzed variables.
  !
  call gsv_setup
  call msg_memUsage(myName)

  trlmFileName = './trlm_01'

  !
  !- Initialize the Analysis grid
  !
  call msg(myName,'Set hco parameters for analysis grid', mpiAll_opt=.false.)
  call hco_SetupFromFile(hco_anl, trlmFileName, aer_backgroundEtiket) ! IN

  !
  !- Initialisation of the analysis grid vertical coordinate from analysisgrid file
  !
  call vco_SetupFromFile( vco_anl,        & ! OUT
                          './analysisgrid') ! IN

  call col_setVco(trlColumnOnAnlLev,vco_anl)
  call msg_memUsage(myName)

  !
  !- Setup and read observations
  !
  call inn_setupObs(obsSpaceData, hco_anl, 'VAR', obsMpiStrategy, varMode) ! IN
  call msg_memUsage(myName)

  !
  !- Basic setup of columnData module
  !
  call col_setup
  call msg_memUsage(myName)

  !
  !- Memory allocation for background column data
  !
  call col_allocate(trlColumnOnAnlLev,obs_numheader(obsSpaceData),mpiLocal_opt=.true.)

  !
  !- Initialize the observation error covariances
  !
  call oer_setObsErrors(obsSpaceData, varMode) ! IN
  call msg_memUsage(myName)

  ! Sea ice concentration
  call filt_iceConcentration(obsSpaceData, beSilent=.false.)
  call filt_backScatAnisIce(obsSpaceData, beSilent=.false.)
  call oer_setErrBackScatAnisIce(obsSpaceData, beSilent=.false.)

  ! Compute the analysis-error
  call aer_analysisError(obsSpaceData, hco_anl, vco_anl, trlmFileName)

  ! Update the Days Since Last Obs
  call aer_daysSinceLastObs(obsSpaceData, hco_anl, vco_anl, trlmFileName)

  ! Now write out the observation data files
  if ( .not. obsf_filesSplit() ) then 
    call msg(myName,'reading/writing global observation files')
    call obs_expandToMpiGlobal(obsSpaceData)
    if (mmpi_myid == 0) call obsf_writeFiles(obsSpaceData)
  else
    ! redistribute obs data to how it was just after reading the files
    call obs_MpiRedistribute(obsSpaceData,OBS_IPF)
    call obsf_writeFiles(obsSpaceData)
  end if

  ! Deallocate copied obsSpaceData
  call obs_finalize(obsSpaceData)

  !
  ! 3. Job termination
  !
  istamp = exfin('ANALYSISERROROI','FIN','NON')

  call utl_tmg_stop(0)

  call tmg_terminate(mmpi_myid, 'TMG_INFO')

  call rpn_comm_finalize(ierr) 

end program midas_analysisErrorOI
