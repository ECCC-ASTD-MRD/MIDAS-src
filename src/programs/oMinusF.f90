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

program midas_ominusf
  !
  !:Purpose: Main program for Observation minus Forecast (O-F) computation.
  !
  !          ---
  !
  !:Algorithm: The non-linear observation operators map a gridded state vector into the 
  !            observation space to compute the difference between the observations 
  !            and that state in observation space. The gridded state vector can be 
  !            background state or the analysis. In case of background state, the 
  !            difference is the innovation vector: ``y-H(xb)``. If asked by 
  !            the user, the diagonal of the background errors standard deviation in 
  !            observation space, :math:`{diag(H B H^{T})}^{1/2}`, are also computed.
  !            The bias corrections are applied for satellite radiances before writing 
  !            the observation files.
  !
  !            --
  !
  !============================================== ==============================================================
  ! Input and Output Files                        Description of file
  !============================================== ==============================================================
  ! ``flnml``                                      In - Main namelist file with parameters user may modify
  ! ``flnml_static``                               In - The "static" namelist that should not be modified
  ! ``trlm_$NN`` (e.g. ``trlm_01``)                In - Background state (a.k.a. trial) files for each timestep
  ! ``analysisgrid``                               In - File defining grid for computing the analysis increment
  ! ``obsfiles_$FAM/obs$FAM_$NNNN_$NNNN``          In - Observation file for each "family" and MPI task
  ! ``obserr``                                     In - Observation error statistics
  ! ``obsfiles_$FAM.updated/obs$FAM_$NNNN_$NNNN``  Out - Updated obs file for each "family" and MPI task
  ! Remainder are files related to radiance obs:
  ! ``stats_$SENSOR_assim``                        In - Satellite radiance observation errors of different sensors
  ! ``stats_tovs``                                 In - Satellite radiance observation errors
  ! ``stats_tovs_symmetricObsErr``                 In - User-defined symmetric TOVS errors for all sky
  ! ``Cmat_$PLATFORM_$SENSOR.dat``                 In - Inter-channel observation-error correlations
  ! ``dynbcor.coeffs.$SENSOR.*.coeffs_$SENSOR``    In - Dynamic bias correction file
  ! ``ceres_global.std``                           In - Surface emmissivity and type?
  ! ``rtcoef_$PLATFORM_$SENSOR.dat``               In - RTTOV coefficient files
  ! ``rttov_h2o_limits.dat``                       In - Min/max humidity limits applied to analysis
  ! ``ozoneclim98``                                In - Ozone climatology
  !============================================== ==============================================================
  !
  !           --
  !
  !:Synopsis: Below is a summary of the ``oMinusF`` program calling sequence:
  !
  !             - **Initial setups:**
  !
  !               - Read the NAMOMF namelist and check/modify some values.
  !
  !               - Various modules are setup: ``obsFiles_mod``, ``timeCoord_mod``.
  !
  !               - Setup ``gridStateVector`` module.
  !
  !               - Setup horizontal and vertical grid objects from either
  !                 ``analysisgrid`` or  first trial file: ``trlm_01``.
  !
  !               - Setup ``obsSpaceData`` object and read observations from
  !                 files: ``inn_setupObs``.
  !
  !               - Applying optional bias corrections to some observation types.
  !
  !               - Setup ``columnData`` module (read list of analysis variables 
  !                 from namelist) and allocate column object.
  !
  !               - Allocate a stateVector object and then read the state 
  !                 (either trials or analysis): ``gio_readTrials``.
  !
  !             - **Computation**
  !
  !               - Compute interpolated column on trial level ``columnTrlOnTrlLev`` 
  !                 from the state: ``inn_setupColumnsOnTrlLev``.
  !
  !               - Compute innovation from background state: ``inn_computeInnovation``.
  !
  !               - For computing background errors in observation space, 
  !                 ``columnTrlOnTrlLev`` are interpolated from background to analysis 
  !                 levels, ``columnTrlOnAnlIncLev``, and the linearize operators are
  !                 initialized: ``inn_setupColumnsOnAnlIncLev``.
  !
  !               - Update radiance bias correction in ``obsSpaceData`` and apply 
  !                 the bias corrections to the observations and innovations for
  !                 radiances: ``bcs_calcBias``, ``bcs_applyBiasCorrection``.
  !
  !               - Write the final bias corrected results into the observation file.
  !
  !             --
  !
  !:Options: `List of namelist blocks <../namelists_in_each_program.html#ominusf>`_
  !          that can affect the ``oMinusF`` program.
  !
  !          * The use of ``oMinusF`` program is controlled by the namelist block
  !           ``&NAMOMF`` read by the ``oMinusF`` program.
  !
  !          * Some of the other relevant namelist blocks used to configure the
  !            ``oMinusF`` are listed in the following table:
  ! 
  !=========================== ========================= =========================================
  ! Module                      Namelist                  Description of what is controlled
  !=========================== ========================= =========================================
  ! ``biasCorrectionConv_mod``  ``NAMBIASCONV``           variables to perform bias correction 
  !                                                       for conventional observations.
  ! ``biasCorrectionConv_mod``  ``NAMSONDETYPES``         additional variables to perform bias 
  !                                                       correction for radiosondes conventional 
  !                                                       observations.
  ! ``biasCorrectionSat_mod``   ``NAMBIASSAT``            variables to perform bias correction 
  !                                                       for satellite radiances.
  ! ``burpread_mod``            ``NAMADDTOBURP``          element IDs to add to the BURP file
  !=========================== ========================= =========================================
  !
  use version_mod
  use oMinusF_mod
  use obsSpaceData_mod
  use columnData_mod
  use obsFiles_mod
  use utilities_mod
  use midasMpi_mod
  use biasCorrectionSat_mod
  use ensembleObservations_mod
  use fileNames_mod
  implicit none

  ! Namelist
  integer :: nEns       ! ensemble size
  logical :: addHBHT    ! choose to add the value of HBHT to obsSpaceData so it can be output
  logical :: addSigmaO  ! choose to add the value of sigma_obs to obsSpaceData so it can be output
  NAMELIST /NAMOMF/addHBHT, addSigmaO, nEns
  
  integer :: fnom, fclos, nulnam, ierr, headerIndex
  
  type(struct_columnData),target  :: columnTrlOnAnlIncLev
  type(struct_columnData),target  :: columnTrlOnTrlLev
  type(struct_obs),       target  :: obsSpaceData
  type(struct_eob),       target  :: ensObs
  
  character(len=20)  :: oMinusFmode
  character(len=256) :: ensPathName = 'ensemble'
  character(len=256) :: ensFileName
  character(len=256) :: trlFileName = './trlm_01'

  logical :: ensFileExists, trlFileExists
  
  call ver_printNameAndVersion('oMinusF','Computation of the innovation')

  !
  !- 1.  Initialization
  !
  
  !- 1.1 mpi
  call mmpi_initialize

  !- 1.2 timings
  call tmg_init(mmpi_myid, 'TMG_INFO')
  call utl_tmg_start(0,'Main')

  if ( mmpi_myid == 0 ) then
    call utl_writeStatus('VAR3D_BEG')
  endif

  !- 1.3 Namelist
  addHBHT   = .false. ! default value
  addSigmaO = .false.
  nEns      = 20

  nulnam = 0
  ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read(nulnam,nml=namomf,iostat=ierr)
  if (ierr /= 0) call utl_abort('midas-OminusF: Error reading namelist')
  if (mmpi_myid == 0) write(*,nml=namomf)
  ierr = fclos(nulnam)

  !- 1.4 Set mode
  inquire(file=trim(trlFileName),exist=trlFileExists)
  call fln_ensFileName(ensFileName, ensPathName, memberIndex_opt=1, &
                       shouldExist_opt=.false.)
  inquire(file=trim(ensFileName),exist=ensFileExists)

  if      (trlFileExists) then
    write(*,*)
    write(*,*) 'Trial/Prog file found'
    write(*,*) 'Setting mode to DETERMINISTIC'
    oMinusFmode = 'deterministic'
  else if (ensFileExists) then
    write(*,*)
    write(*,*) 'Ensemble file found'
    write(*,*) 'Setting mode to ENSEMBLE'
    oMinusFmode = 'ensemble'
  else
    write(*,*)
    write(*,*) 'trlFileName = ', trim(trlFileName)
    write(*,*) 'ensFileName = ', trim(ensFileName)
    call utl_abort('oMinusF : did not find a trial/prog or ensemble file')
  end if

  !
  !- 2.  Calculate the Observation - Forecast difference
  !
  if (trim(oMinusFmode) == 'deterministic') then

    !- 2.1 Compute O-F and store in obsSpaceDate
    call omf_oMinusF(columnTrlOnAnlIncLev, columnTrlOnTrlLev, obsSpaceData, &
                     'OminusF', addHBHT, addSigmaO)

    !- 2.2 Remove biases
    call bcs_calcBias(obsSpaceData,columnTrlOnTrlLev) ! Fill in OBS_BCOR obsSpaceData column with computed bias correction

    call bcs_applyBiasCorrection(obsSpaceData,OBS_VAR,"TO") ! Apply bias correction to OBS
    call bcs_applyBiasCorrection(obsSpaceData,OBS_OMP,"TO") ! Apply bias correction to O-F

    !- 2.3 Write the results

    !  2.3.1 Into the listings
    write(*,*)
    write(*,*) '> midas-OminusF: printing the FIRST header and body'
    do headerIndex = 1, min(1,obs_numHeader(obsSpaceData))
      call obs_prnthdr(obsSpaceData,headerIndex)
      call obs_prntbdy(obsSpaceData,headerIndex)
    end do
    ! 2.3.2 Into the observation files
    write(*,*)
    write(*,*) '> midas-OminusF: writing to file'
    call obsf_writeFiles(obsSpaceData)

  else ! ensemble

    !- 2.1 Compute O-F and store in ensObs
    call omf_oMinusFens( ensObs, obsSpaceData, nEns, ensPathName, &
                         'OminusF', addHBHT, addSigmaO)

    !- 2.3 Write the results
    call obsf_writeFiles(obsSpaceData, ensObs_opt=ensObs)
    
  end if  
    
  !
  !- 3.  Ending
  !
  write(*,*)
  write(*,*) '> midas-OminusF: Ending'

  call obs_finalize(obsSpaceData) ! deallocate obsSpaceData
  
  call utl_tmg_stop(0)
  call tmg_terminate(mmpi_myid, 'TMG_INFO')

  call rpn_comm_finalize(ierr)

  if ( mmpi_myid == 0 ) then
    call utl_writeStatus('VAR3D_END')
  endif

end program midas_ominusf
