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

program midas_var
  !
  ! :Purpose: Main program for variational minimization or background check 
  !           (depending on the mode selected in the namelist).
  !
  use ramDisk_mod
  use utilities_mod
  use mpi_mod
  use mpivar_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use obsSpaceData_mod
  use columnData_mod  
  use gridStateVector_mod
  use obsSpaceDiag_mod
  use controlVector_mod
  use obsFiles_mod
  use obsFilter_mod  
  use minimization_mod
  use innovation_mod
  use minimization_mod
  use analysisGrid_mod
  use bmatrix_mod
  use tovs_nl_mod
  use obsErrors_mod
  use gridVariableTransforms_mod
  use obsOperators_mod
  use statetocolumn_mod
  use multi_ir_bgck_mod
  use biasCorrectionSat_mod
  use increment_mod
  use residual_mod
  use stateToColumn_mod
  use backgroundCheck_mod
  use biasCorrectionConv_mod

  implicit none

  integer :: istamp,exdb,exfin
  integer :: ierr,nconf

  type(struct_obs),        target  :: obsSpaceData
  type(struct_columnData), target  :: trlColumnOnAnlLev
  type(struct_columnData), target  :: trlColumnOnTrlLev
  type(struct_gsv)                 :: stateVectorIncr
  type(struct_gsv)                 :: stateVectorTrial
  type(struct_hco), pointer        :: hco_anl => null()
  type(struct_vco), pointer        :: vco_anl => null()

  real(8), allocatable :: controlVectorIncr(:)

  character(len=9)  :: clmsg
  character(len=48) :: obsMpiStrategy, varMode

  logical :: writeAnalysis

  NAMELIST /NAMCT0/nconf,writeAnalysis

  integer :: nulnam, fnom, fclos, get_max_rss, headerIndex

  istamp = exdb('VAR','DEBUT','NON')

  write(*,'(/,' //                                                &
            '3(" *****************"),/,' //                       &
            '14x,"-- START OF MAIN PROGRAM MIDAS-VAR: --",/,' //   &
            '14x,"-- VARIATIONAL ASSIMILATION          --",/, ' //&
            '14x,"-- Revision : ",a," --",/,' //       &
            '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  ! MPI initialization
  call mpi_initialize  

  call tmg_init(mpi_myid, 'TMG_VAR' )

  call tmg_start(1,'MAIN')

  if(mpi_myid == 0) then
    clmsg = 'VAR3D_BEG'
    call utl_writeStatus(clmsg)
  endif 

  ! 1. Top level setup
  nconf             = 141
  writeAnalysis = .false.

  nulnam=0
  ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  if(ierr.ne.0) call utl_abort('midas-var: Error opening file flnml')
  read(nulnam,nml=namct0,iostat=ierr)
  if(ierr.ne.0) call utl_abort('midas-var: Error reading namelist')
  write(*,nml=namct0)
  ierr=fclos(nulnam)

  write(*,*)
  select case(nconf)
  case (141)
    write(*,*) 'midas-var: Analysis mode selected'
    varMode='analysis'
  case (111)
    write(*,*) 'midas-var: Background check for IR sat. data mode selected'
    varMode='bgckIR'
  case (101)
    write(*,*) 'midas-var: Background check for conventional obs mode selected'
    varMode='bgckConv'
  case default
    write(*,*) 'midas-var: Unknown mode ', nconf
    call utl_abort('midas-var')
  end select

  call ram_setup

  !
  !- Read variational bias correction namelist (default is to not use it)
  !
  call bcs_readConfig()

  ! 2. Decide on configuration of job

  ! ---BGCHECK (conventional obs)--- !
  if ( trim(varMode) == 'bgckConv' ) then
    if(mpi_myid == 0) write(*,*) 'MIDAS-VAR: CONVENTIONNAL BGCHECK MODE'

    ! Do initial set up
    call tmg_start(2,'PREMIN')

    obsMpiStrategy = 'LIKESPLITFILES'

    call var_setup('ALL') ! obsColumnMode   
    
    ! Apply optional bias corrections when namelist logicals aiBiasActive, gpBiasActive are TRUE
    call bcc_applyAIBcor(obsSpaceData)    
    call bcc_applyGPBcor(obsSpaceData)
    
    ! Reading, horizontal interpolation and unit conversions of the 3D trial fields
    call inn_setupBackgroundColumns( trlColumnOnTrlLev, obsSpaceData )

    ! Interpolate trial columns to analysis levels and setup for linearized H
    call inn_setupBackgroundColumnsAnl(trlColumnOnTrlLev,trlColumnOnAnlLev,obsSpaceData)

    ! Compute observation innovations and prepare obsSpaceData for minimization
    call inn_computeInnovation(trlColumnOnTrlLev,obsSpaceData)
    call tmg_stop(2)

    ! Do the background check and output the observation data files
    call bgck_bgcheck_conv(trlColumnOnAnlLev,trlColumnOnTrlLev,obsSpaceData)
 
    ! Write out contents of obsSpaceData into observation files
    call obsf_writeFiles(obsSpaceData)

    ! Print the FIRST header and body
    if (mpi_myid == 0 ) then
      do headerIndex = 1, min(1,obs_numHeader(obsSpaceData))
        call obs_prnthdr(obsSpaceData,headerIndex)
        call obs_prntbdy(obsSpaceData,headerIndex)
      end do
    end if
 
    ! deallocate obsSpaceData
    call obs_finalize(obsSpaceData)

  ! ---BGCHECK (AIRS, IASI, CrIS)--- !
  else if ( trim(varMode) == 'bgckIR' ) then
    if(mpi_myid == 0) write(*,*) 'MIDAS-VAR: HYPERSPECTRAL IR BGCHECK MODE'

    ! Do initial set up
    call tmg_start(2,'PREMIN')

    obsMpiStrategy = 'LIKESPLITFILES'

    call var_setup('ALL') ! obsColumnMode   

    ! Reading, horizontal interpolation and unit conversions of the 3D trial fields
    call inn_setupBackgroundColumns( trlColumnOnTrlLev, obsSpaceData )

    ! Compute observation innovations and prepare obsSpaceData for minimization
    call inn_computeInnovation(trlColumnOnTrlLev,obsSpaceData)
    call tmg_stop(2)

    call bcs_calcBias(obsSpaceData,trlColumnOnTrlLev) ! Fill in OBS_BCOR obsSpaceData column with computed bias correction

    call bcs_applyBiasCorrection(obsSpaceData,OBS_VAR,"TO") ! Apply bias correction to OBS

    call bcs_applyBiasCorrection(obsSpaceData,OBS_OMP,"TO") ! Apply bias correction to O-F

    ! Do the IR background check

    call irbg_bgCheckIR(trlColumnOnTrlLev,obsSpaceData)

    !  Write out contents of obsSpaceData into BURP files
    call obsf_writeFiles(obsSpaceData)

    !  Add cloud parameter data to burp files (AIRS,IASI,CrIS,...)
    call obsf_addCloudParametersAndEmissivity(obsSpaceData)

    do headerIndex =1, min(1,obs_numHeader(obsSpaceData))
      call obs_prnthdr(obsSpaceData,headerIndex)
      call obs_prntbdy(obsSpaceData,headerIndex)
    end do

    ! Deallocate obsSpaceData
    call obs_finalize(obsSpaceData)


  ! ---ANALYSIS MODE--- !
  else if ( trim(varMode) == 'analysis' ) then
    write(*,*) 'MIDAS-VAR: ANALYSIS MODE'

    ! Do initial set up
    call tmg_start(2,'PREMIN')

    obsMpiStrategy = 'LIKESPLITFILES'

    call var_setup('VAR') ! obsColumnMode
    call tmg_stop(2)

    ! Read trials and horizontally interpolate to columns
    call tmg_start(2,'PREMIN')
    if (writeAnalysis) then
      call inn_setupBackgroundColumns( trlColumnOnTrlLev, obsSpaceData,  &
                                       stateVectorTrialOut_opt=stateVectorTrial )
    else
      call inn_setupBackgroundColumns( trlColumnOnTrlLev, obsSpaceData )
    end if

    !
    !- Initialize the background-error covariance, also sets up control vector module (cvm)
    !
    call bmat_setup(hco_anl,vco_anl)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    ! - Initialize the gridded variable transform module
    !
    call gvt_setup(hco_anl,vco_anl)

    !
    !- Set up the minimization module, now that the required parameters are known
    !  NOTE: some global variables remain in minimization_mod that must be initialized before
    !        inn_setupBackgroundColumns
    !
    call min_setup( cvm_nvadim ) ! IN
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! Interpolate trial columns to analysis levels and setup for linearized H
    call inn_setupBackgroundColumnsAnl(trlColumnOnTrlLev,trlColumnOnAnlLev,obsSpaceData)

    ! Compute observation innovations and prepare obsSpaceData for minimization
    call inn_computeInnovation(trlColumnOnTrlLev,obsSpaceData)
    call tmg_stop(2)

    allocate(controlVectorIncr(cvm_nvadim),stat=ierr)
    if(ierr.ne.0) then
      write(*,*) 'var: Problem allocating memory for ''controlVectorIncr''',ierr
      call utl_abort('aborting in VAR')
    endif

    ! Do minimization of cost function
    call min_minimize(trlColumnOnAnlLev,obsSpaceData,controlVectorIncr)

    ! Compute satellite bias correction increment and write to file
    call bcs_writebias(controlVectorIncr)

    call gsv_allocate(stateVectorIncr, tim_nstepobsinc, hco_anl, vco_anl, &
         datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true., &
         allocHeight_opt=.false., allocPressure_opt=.false.)

    ! get final increment
    call inc_getIncrement(controlVectorIncr,stateVectorIncr,cvm_nvadim)

    ! output the analysis increment
    call tmg_start(6,'WRITEINCR')
    call inc_writeIncrement(stateVectorIncr) ! IN
    call tmg_stop(6)

    ! Conduct obs-space post-processing diagnostic tasks (some diagnostic
    ! computations controlled by NAMOSD namelist in flnml)
    call osd_ObsSpaceDiag(obsSpaceData,trlColumnOnAnlLev)

    ! Deallocate memory related to B matrices
    call bmat_finalize()

    ! compute and write the analysis (as well as the increment on the trial grid)
    if (writeAnalysis) then
      call tmg_start(18,'ADDINCREMENT')
      call inc_computeAndWriteAnalysis(stateVectorIncr,                    &  ! IN
                                       stateVectorTrial_opt=stateVectorTrial) ! IN
      call tmg_stop(18)
    end if

    if (mpi_myid == 0) then
      clmsg = 'REBM_DONE'
      call utl_writeStatus(clmsg)
    end if

    call gsv_deallocate(stateVectorIncr)

    ! write the Hessian
    call min_writeHessian(controlVectorIncr)
    deallocate(controlVectorIncr)

    ! Deallocate memory related to variational bias correction
    call bcs_finalize()

    ! Now write out the observation data files
    if(min_niter.gt.0) then
      if ( .not. obsf_filesSplit() ) then 
        write(*,*) 'We read/write global observation files'
        call obs_expandToMpiGlobal(obsSpaceData)
        if (mpi_myid == 0) call obsf_writeFiles(obsSpaceData)
      else
        ! redistribute obs data to how it was just after reading the files
        call obs_MpiRedistribute(obsSpaceData,OBS_IPF)
        call obsf_writeFiles(obsSpaceData)
      end if
    end if

    ! Deallocate copied obsSpaceData
    call obs_finalize(obsSpaceData)

  else

    write(*,*) ' MIDAS-VAR: ERROR, UNKNOWN NCONF SPECIFIED'

  end if

  ! 3. Job termination

  istamp = exfin('VAR','FIN','NON')

  if(mpi_myid == 0) then
    clmsg = 'VAR3D_END'
    call utl_writeStatus(clmsg)
  endif

  call tmg_stop(1)

  call tmg_terminate(mpi_myid, 'TMG_VAR' )

  call rpn_comm_finalize(ierr) 

contains
  !--------------------------------------------------------------------------
  !- var_setup
  !--------------------------------------------------------------------------
  subroutine var_setup(obsColumnMode)
    implicit none

    character (len=*) :: obsColumnMode
    integer :: datestamp
    type(struct_hco),pointer :: hco_core => null()

    integer :: get_max_rss

    write(*,*) ''
    write(*,*) '-----------------------------------'
    write(*,*) '-- Starting subroutine var_setup --'
    write(*,*) '-----------------------------------'

    !
    !- Initialize the Temporal grid
    !
    call tim_setup

    !     
    !- Initialize observation file names and set datestamp
    !
    call obsf_setup( dateStamp, varMode )
    if ( dateStamp > 0 ) then
      call tim_setDatestamp(datestamp)     ! IN
    else
      call utl_abort('var_setup: Problem getting dateStamp from observation file')
    end if

    !
    !- Initialize constants
    !
    if(mpi_myid.eq.0) call mpc_printConstants(6)

    !
    !- Initialize variables of the model states
    !
    call gsv_setup
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Initialize the Analysis grid
    !
    if(mpi_myid.eq.0) write(*,*)''
    if(mpi_myid.eq.0) write(*,*)'var_setup: Set hco parameters for analysis grid'
    call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

    if ( hco_anl % global ) then
      call agd_SetupFromHCO( hco_anl ) ! IN
    else
      !- Initialize the core (Non-Extended) analysis grid
      if(mpi_myid.eq.0) write(*,*)'var_setup: Set hco parameters for core grid'
      call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore' ) ! IN
      !- Setup the LAM analysis grid metrics
      call agd_SetupFromHCO( hco_anl, hco_core ) ! IN
    end if

    !     
    !- Initialisation of the analysis grid vertical coordinate from analysisgrid file
    !
    call vco_SetupFromFile( vco_anl,        & ! OUT
                            './analysisgrid') ! IN

    call col_setVco(trlColumnOnAnlLev,vco_anl)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Setup and read observations
    !
    call inn_setupObs(obsSpaceData, obsColumnMode, obsMpiStrategy, varMode) ! IN
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Basic setup of columnData module
    !
    call col_setup
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Memory allocation for background column data
    !
    call col_allocate(trlColumnOnAnlLev,obs_numheader(obsSpaceData),mpiLocal_opt=.true.)

    !
    !- Initialize the observation error covariances
    !
    call oer_setObsErrors(obsSpaceData, varMode) ! IN
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine var_setup

end program midas_var
