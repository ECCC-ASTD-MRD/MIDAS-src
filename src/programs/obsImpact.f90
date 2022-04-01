!------------------------------------- LICENCE BEGIN -------------------------------------
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

program midas_obsimpact
  !
  ! :Purpose: Main program for Observation Impact computation (FSOI)
  !
  use version_mod
  use codePrecision_mod
  use ramDisk_mod
  use utilities_mod
  use mpi_mod
  use mpivar_mod
  use MathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use columnData_mod
  use obsSpaceData_mod
  use controlVector_mod
  use gridStateVector_mod
  use bmatrix_mod
  use bmatrixensemble_mod
  use stateToColumn_mod
  use obsOperators_mod
  use costFunction_mod
  use quasinewton_mod
  use innovation_mod
  use obsFiles_mod
  use obsFilter_mod
  use obsErrors_mod
  use gridVariableTransforms_mod
  use rttov_const, only :inst_name, platform_name
  use tovs_nl_mod
  use rMatrix_mod

  implicit none

  integer :: istamp,exdb,exfin,ierr

  type(struct_obs),       target :: obsSpaceData
  type(struct_columnData),target :: columnTrlOnAnlIncLev
  type(struct_columnData),target :: columnTrlOnTrlLev
  type(struct_gsv)               :: stateVectorTrialHighRes

  character(len=48) :: obsMpiStrategy
  character(len=3)  :: obsColumnMode

  type(struct_obs),        pointer :: obsSpaceData_ptr
  type(struct_columnData), pointer :: columnTrlOnAnlIncLev_ptr
  type(struct_columnData), pointer :: column_ptr
  type(struct_vco),        pointer :: vco_anl => null()
  type(struct_hco),        pointer :: hco_anl => null()
  type(struct_hco),        pointer :: hco_trl => null()
  type(struct_vco),        pointer :: vco_trl => null()
  type(struct_hco),        pointer :: hco_core => null()
  integer             :: fso_nsim, nvadim_mpilocal, dateStamp
  integer,external    :: get_max_rss
  real(8),allocatable :: vhat(:)
  logical             :: allocHeightSfc

  ! namelist variables
  integer             :: nvamaj, nitermax, nsimmax
  real(8)             :: leadTime, repsg, rdf1fac,latMinNorm, latMaxNorm, lonMinNorm, lonMaxNorm
  logical             :: includeUVnorm, includeTTnorm, includeP0norm, includeHUnorm, includeTGnorm 
  character(len=256)  :: forecastPath
  character(len=4)    :: fsoMode

  NAMELIST /NAMFSO/leadTime, nvamaj, nitermax, nsimmax
  NAMELIST /NAMFSO/repsg, rdf1fac, forecastPath, fsoMode
  NAMELIST /NAMFSO/latMinNorm, latMaxNorm, lonMinNorm, lonMaxNorm
  NAMELIST /NAMFSO/includeUVnorm, includeTTnorm, includeP0norm, includeHUnorm, includeTGnorm

  istamp = exdb('OBSIMPACT','DEBUT','NON')

  call ver_printNameAndVersion('obsImpact','Calculation of observation impact')

  ! MPI initilization
  call mpi_initialize

  call tmg_init(mpi_myid, 'TMG_INFO')

  call utl_tmg_start(0,'Main')

  if (mpi_myid == 0) then
    call utl_writeStatus('VAR3D_BEG')
  end if

  call ram_setup

  !
  !- 1. Settings 
  !
  obsColumnMode  = 'VAR'
  obsMpiStrategy = 'LIKESPLITFILES'

  !- Do initial set up
  call fso_setup

  !
  !- Initialize the Temporal grid
  !
  call tim_setup
  !     
  !- Initialize burp file names and set datestamp
  !
  call obsf_setup( dateStamp, 'FSO' )
  if ( dateStamp > 0 ) then
    call tim_setDatestamp(datestamp)     ! IN
  else
    call utl_abort('obsImpact: Problem getting dateStamp from observation file')
  end if
  !
  !- Initialize constants
  !
  if ( mpi_myid == 0 ) then
    call mpc_printConstants(6)
    call pre_printPrecisions
  end if

  !
  !- Initialize variables of the model states
  !
  call gsv_setup
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
  !
  !- Initialize the Analysis grid
  !
  if (mpi_myid.eq.0) write(*,*)''
  if (mpi_myid.eq.0) write(*,*)'obsImpact: Set hco parameters for analysis grid'
  call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

  if ( hco_anl % global ) then
    hco_core => hco_anl
  else
    !- Iniatilized the core (Non-Exteded) analysis grid
    call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore' ) ! IN
  end if

  !     
  !- Initialisation of the analysis grid vertical coordinate from analysisgrid file !
  call vco_SetupFromFile( vco_anl,        & ! OUT
                          './analysisgrid') ! IN

  call col_setVco(columnTrlOnAnlIncLev,vco_anl)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- Setup and read observations
  !
  call inn_setupObs(obsSpaceData, hco_anl, obsColumnMode, obsMpiStrategy, 'FSO') ! IN
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- Basic setup of columnData module
  !
  call col_setup
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- Memory allocation for background column data
  !
  call col_allocate(columnTrlOnAnlIncLev,obs_numheader(obsSpaceData),mpiLocal_opt=.true.)

  !
  !- Initialize the observation error covariances
  !
  call oer_setObsErrors(obsSpaceData, 'FSO') ! IN
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Reading trials
  call inn_getHcoVcoFromTrlmFile( hco_trl, vco_trl )
  allocHeightSfc = ( vco_trl%Vcode /= 0 )

  call gsv_allocate( stateVectorTrialHighRes, tim_nstepobs, hco_trl, vco_trl,  &
                     dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                     mpi_distribution_opt='Tiles', dataKind_opt=4,  &
                     allocHeightSfc_opt=allocHeightSfc, hInterpolateDegree_opt='LINEAR', &
                     beSilent_opt=.false. )
  call gsv_zero( stateVectorTrialHighRes )
  call gsv_readTrials( stateVectorTrialHighRes )
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Horizontally interpolate trials to trial columns
  call inn_setupColumnsOnTrlLev( columnTrlOnTrlLev, obsSpaceData, hco_core, &
                                   stateVectorTrialHighRes )
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !- Initialize the background-error covariance, also sets up control vector module (cvm)
  !
  call bmat_setup(hco_anl,hco_core,vco_anl)
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  ! - Initialize the gridded variable transform module
  !
  call gvt_setup(hco_anl,hco_core,vco_anl)
  call gvt_setupRefFromTrialFiles('HU')
  call gvt_setupRefFromTrialFiles('height')

  !
  !- 2. Do the actual job
  !

  ! Interpolate trial columns to analysis levels and setup for linearized H
  call inn_setupColumnsOnAnlIncLev(columnTrlOnTrlLev,columnTrlOnAnlIncLev)

  ! Compute observation innovations and prepare obsSpaceData for minimization
  call inn_computeInnovation(columnTrlOnTrlLev,obsSpaceData)

  ! Perform forecast sensitivity to observation calculation using ensemble approach 
  call fso_ensemble(columnTrlOnAnlIncLev,obsSpaceData)

  ! Deallocate memory related to B matrices
  call bmat_finalize()

  ! Now write out the observation data files
  if ( .not. obsf_filesSplit() ) then
    write(*,*) 'We read/write global observation files'
    call obs_expandToMpiGlobal(obsSpaceData)
    if (mpi_myid == 0) call obsf_writeFiles(obsSpaceData)
  else
    ! redistribute obs data to how it was just after reading the files
    call obs_MpiRedistribute(obsSpaceData,OBS_IPF)
    call obsf_writeFiles(obsSpaceData)
  end if
  
  ! Deallocate copied obsSpaceData
  call obs_finalize(obsSpaceData)

  !
  !- 3. Job termination
  !
  istamp = exfin('OBSIMPACT','FIN','NON')

  if (mpi_myid == 0) then
    call utl_writeStatus('VAR3D_END')
  endif

  call tmg_stop(0)

  call tmg_terminate(mpi_myid, 'TMG_INFO')

  call rpn_comm_finalize(ierr)

contains

  subroutine fso_setup
    implicit none
  
    integer :: ierr,nulnam
    integer :: fnom,fclos

    ! set default values for namelist variables
    leadtime = 12.0d0
    nvamaj = 6
    nitermax = 100
    nsimmax  = 120
    repsg    = 1d-5
    rdf1fac  = 0.25d0
    latMinNorm = -95.0d0
    latMaxNorm = 95.0d0
    lonMinNorm = -185.0d0
    lonMaxNorm = 365.0d0
    includeUVnorm=.true.
    includeTTnorm=.true.
    includeP0norm=.true.
    includeHUnorm=.false.
    includeTGnorm=.false.
    forecastPath = './forecasts'
    fsoMode  = 'HFSO' 

    ! read in the namelist NAMFSO
    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=namfso,iostat=ierr)
    if(ierr /= 0) call utl_abort('fso_setup: Error reading namelist')
    write(*,nml=namfso)
    ierr = fclos(nulnam)
    ! convert Latmin,max and Lonmin,max from degres into RAD
    latMinNorm = latMinNorm*MPC_RADIANS_PER_DEGREE_R8 
    latMaxNorm = latMaxNorm*MPC_RADIANS_PER_DEGREE_R8 
    lonMinNorm = lonMinNorm*MPC_RADIANS_PER_DEGREE_R8 
    lonMaxNorm = lonMaxNorm*MPC_RADIANS_PER_DEGREE_R8 

    call ben_setFsoLeadTime(leadTime)
    fso_nsim = 0

  end subroutine fso_setup
 

  subroutine fso_ensemble(columnTrlOnAnlIncLev,obsSpaceData)
    implicit none

    type(struct_columnData),target  :: columnTrlOnAnlIncLev
    type(struct_obs),target         :: obsSpaceData
    type(struct_columnData),target  :: column
    type(struct_gsv)                :: statevector_FcstErr, statevector_fso
    type(struct_vco), pointer       :: vco_anl
    real(8),allocatable             :: ahat(:), zhat(:)
    integer                         :: dateStamp_fcst, dateStamp
    
    !for Observation space 
    integer                         :: headerIndex, bodyIndexBeg, bodyIndexEnd, bodyIndex
    real(8)                         :: fso_ori, fso_fin

    if (mpi_myid == 0) write(*,*) 'fso_ensemble: starting'

    vco_anl => col_getVco(columnTrlOnAnlIncLev)

    nvadim_mpilocal = cvm_nvadim

    ! initialize column object for storing "increment"
    call col_setVco(column,col_getVco(columnTrlOnAnlIncLev))
    call col_allocate(column,col_getNumCol(columnTrlOnAnlIncLev),mpiLocal_opt=.true.)

    ! compute dateStamp_fcst
    call incdatr(dateStamp_fcst, tim_getDatestamp(), leadTime)
    dateStamp = tim_getDatestamp()
    write(*,*) 'fso_ensemble: analysis datestamp = ',dateStamp
    write(*,*) 'fso_ensemble: forecast datestamp = ',dateStamp_fcst

    ! allocate control vector related arrays (these are all mpilocal)
    allocate(ahat(nvadim_mpilocal))
    allocate(vhat(nvadim_mpilocal))
    allocate(zhat(nvadim_mpilocal))

    ! initialize control vector related arrays to zero
    ahat(:) = 0.0d0
    vhat(:) = 0.0d0

    ! for statevector_FcstErr
    call gsv_allocate(statevector_FcstErr, 1, hco_anl, vco_anl, &
                      datestamp_opt=datestamp_fcst, mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)

    ! for statevector_fso 
    call gsv_allocate(statevector_fso, tim_nstepobsinc, hco_anl, vco_anl, &
                      dataKind_opt=pre_incrReal, &
                      datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true.)

    ! compute forecast error = C * (error_t^fa + error_t^fb)  
    call fso_calcFcstError(columnTrlOnAnlIncLev,statevector_FcstErr)
   

    ! compute vhat = B_t^T/2 * C * (error_t^fa + error_t^fb)  
    call bmat_sqrtBT(vhat, nvadim_mpilocal, statevector_FcstErr, useFSOFcst_opt = .true.)

    if (mpi_myid == 0) write(*,*) maxval(vhat),minval(vhat)

    if( trim(fsoMode) == 'HFSO' ) then
      call fso_minimize(nvadim_mpilocal, zhat, column, columnTrlOnAnlIncLev, obsSpaceData)
      ahat = zhat + vhat
      call bmat_sqrtB(ahat, nvadim_mpilocal, statevector_fso)
    elseif( trim(fsoMode) == 'EFSO' ) then
      call bmat_sqrtB(vhat, nvadim_mpilocal, statevector_fso)
    end if

    ! Compute yhat = [R^-1 H B^1/2 ahat], and put in OBS_FSO
    call s2c_tl(statevector_fso,column,columnTrlOnAnlIncLev,obsSpaceData)  ! put in column H_horiz B^1/2 ahat
    call oop_Htl(column,columnTrlOnAnlIncLev,obsSpaceData,1)          ! Save as OBS_WORK: H_vert H_horiz B^1/2 vhat = H B^1/2 ahat
    call rmat_RsqrtInverseAllObs(obsSpaceData,OBS_FSO,OBS_WORK) ! Save as OBS_FSO : R**-1/2 H B^1/2 ahat
    call rmat_RsqrtInverseAllObs(obsSpaceData,OBS_FSO,OBS_FSO)  ! Save as OBS_FSO : R**-1 H B^1/2 ahat\

    do headerIndex = 1, obs_numHeader(obsSpaceData)

      bodyIndexBeg = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
      bodyIndexEnd = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + bodyIndexBeg - 1

      do bodyIndex = bodyIndexBeg, bodyIndexEnd
        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == 1 ) then
          fso_ori = obs_bodyElem_r(obsSpaceData,OBS_FSO,bodyIndex)
          fso_fin = fso_ori * obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)
          call obs_bodySet_r(obsSpaceData,OBS_FSO,bodyIndex, fso_fin)
        end if
      end do

    end do

    ! print out the information of total FSO for each family
    call fso_sumFSO(obsSpaceData)

    ! deallocate the control vector related arrays
    deallocate(ahat)
    deallocate(vhat)
    deallocate(zhat)
    call col_deallocate(column)

    if (mpi_myid == 0) write(*,*) 'end of fso_ensemble'

  end subroutine fso_ensemble

  subroutine fso_calcFcstError(columnTrlOnAnlIncLev,statevector_out)
    !
    ! In this subroutine it reads the forecast from background and analysis, the verifying analysis
    ! Based on these inputs, it calculates the Forecast error
    !
    implicit none
    
    type(struct_columnData),target  :: columnTrlOnAnlIncLev
    type(struct_gsv)                :: statevector_fa, statevector_fb, statevector_a
    type(struct_gsv)                :: statevector_out
    character(len=256)              :: fileName_fa, fileName_fb, fileName_a
    logical                         :: faExists
    type(struct_gsv)                :: statevector_tempfa, statevector_tempfb

    type(struct_vco), pointer       :: vco_anl
    integer                         :: dateStamp_fcst, dateStamp
    
    vco_anl => col_getVco(columnTrlOnAnlIncLev)

    ! compute dateStamp_fcst
    call incdatr(dateStamp_fcst, tim_getDatestamp(), leadTime)
    dateStamp = tim_getDatestamp()
    write(*,*) 'fso_ensemble: analysis datestamp = ',dateStamp
    write(*,*) 'fso_ensemble: forecast datestamp = ',dateStamp_fcst

    ! read forecasts from the analysis and background state
    fileName_fa = trim(forecastPath) // '/forecast_a'
    inquire(file=trim(fileName_fa),exist=faExists)
    write(*,*) 'faExists', faExists
    call gsv_allocate(statevector_fa, 1, hco_anl, vco_anl, &
                      datestamp_opt=datestamp_fcst, mpi_local_opt=.true., &
                      hInterpolateDegree_opt='LINEAR', &
                      allocHeight_opt=.false., allocPressure_opt=.false.)
    call gsv_readFromFile(statevector_fa, fileName_fa, ' ', 'P', containsFullField_opt=.true.)

    !for statevecotr_tempfa
    call gsv_allocate(statevector_tempfa, 1, hco_anl, vco_anl, &
                      datestamp_opt=datestamp_fcst, mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)

    !for statevector_fb
    fileName_fb = trim(forecastPath) // '/forecast_b'
    call gsv_allocate(statevector_fb, 1, hco_anl, vco_anl, &
                      datestamp_opt=datestamp_fcst, mpi_local_opt=.true., &
                      hInterpolateDegree_opt='LINEAR', &
                      allocHeight_opt=.false., allocPressure_opt=.false.)
    call gsv_readFromFile(statevector_fb, fileName_fb, ' ', 'P', containsFullField_opt=.true.)

    !for statevecotr_tempfb
    call gsv_allocate(statevector_tempfb, 1, hco_anl, vco_anl, &
                      datestamp_opt=datestamp_fcst, mpi_local_opt=.true., &
                      allocHeight_opt=.false., allocPressure_opt=.false.)

    ! read verifying analysis
    fileName_a = trim(forecastPath) // '/analysis'
    call gsv_allocate(statevector_a, 1,hco_anl, vco_anl, &
                      datestamp_opt=datestamp_fcst, mpi_local_opt=.true., &
                      hInterpolateDegree_opt='LINEAR', &
                      allocHeight_opt=.false., allocPressure_opt=.false.)
    call gsv_readFromFile(statevector_a, fileName_a, ' ', 'A', containsFullField_opt=.true.)

    ! compute error of both forecasts (overwrite forecasts with error)
    call gsv_add(statevector_a, statevector_fa, -1.0d0)
    call gsv_add(statevector_a, statevector_fb, -1.0d0)

    call gsv_copy(statevector_fa,statevector_tempfa)
    call gsv_copy(statevector_fb,statevector_tempfb)
    call gsv_multEnergyNorm(statevector_tempfa, statevector_a,  &
                            latMinNorm, latMaxNorm, lonMinNorm, lonMaxNorm, & 
                            includeUVnorm, includeTTnorm, includeP0norm,  &
                            includeHUnorm, includeTGnorm) ! use analysis as reference state
    call gsv_multEnergyNorm(statevector_tempfb, statevector_a,  &
                            latMinNorm, latMaxNorm, lonMinNorm, lonMaxNorm, & 
                            includeUVnorm, includeTTnorm, includeP0norm,  &
                            includeHUnorm, includeTGnorm) ! use analysis as reference state

    ! compute error Norm =  C * (error_t^fa + error_t^fb)
    call gsv_add(statevector_fa, statevector_fb, 1.0d0)
    call gsv_multEnergyNorm(statevector_fb, statevector_a,  &
                            latMinNorm, latMaxNorm, lonMinNorm, lonMaxNorm, & 
                            includeUVnorm, includeTTnorm, includeP0norm,  &
                            includeHUnorm, includeTGnorm) ! use analysis as reference state
    call gsv_copy(statevector_fb,statevector_out)

  end subroutine fso_calcFcstError

  subroutine fso_minimize(nvadim,zhat,column,columnTrlOnAnlIncLev,obsSpaceData)
    implicit none

    type(struct_columnData),target  :: columnTrlOnAnlIncLev, column
    type(struct_obs),target         :: obsSpaceData
    real(8),dimension(nvadim)       :: zhat
    real(8),allocatable             :: gradJ(:), vatra(:)
    ! for minimization
    integer                         :: imode, itermax, isimmax, indic,nvadim, nmtra
    real(8)                         :: zjsp, zxmin, zdf1, zeps, dlgnorm, dlxnorm,zspunused(1)
    integer                         :: impres, iztrl(10), intunused(1)
    real                            :: rspunused(1)
    integer                         :: nulout = 6

   
    call utl_tmg_start(90,'--Minimization')

    if (mpi_myid == 0) write(*,*) 'fso_minimize: starting'

    nmtra = (4 + 2*nvamaj)*nvadim
    write(*,'(4X,"NVAMAJ = ",I3,/5X,"NMTRA =",I14)') nvamaj,nmtra

    columnTrlOnAnlIncLev_ptr => columnTrlOnAnlIncLev
    column_ptr  => column
    obsSpaceData_ptr => obsSpaceData
  
    allocate(gradJ(nvadim))
    allocate(vatra(nmtra))

    gradJ(:) = 0.0d0
    vatra(:) = 0.0d0
    zhat(:)  = 0.0d0
    
    ! Compute zhat by performing variational minimization
    ! Set-up for the minimization
    if (mpi_myid == 0) then
      impres = 5
    else
      impres = 0
    end if

    imode = 0
    zeps = repsg
    itermax = nitermax
    isimmax = nsimmax
    zxmin = epsilon(zxmin)
    ! initial gradient calculation
    indic = 2
    call utl_tmg_start(91,'----QuasiNewton')
    call simvar(indic,nvadim,zhat,zjsp,gradJ)
    call tmg_stop(91)
    zdf1 =  rdf1fac * abs(zjsp)

    ! print amplitude of initial gradient and cost function value
    call prscal(nvadim,gradJ,gradJ,dlgnorm)
    dlgnorm = dsqrt(dlgnorm)
    call prscal(nvadim,zhat,zhat,dlxnorm)
    dlxnorm = dsqrt(dlxnorm)
    write(*,*)' |X| = ', dlxnorm
    write(*,'(/4X,"J(X) = ",G23.16,4X,"|Grad J(X)| = ",G23.16)') zjsp, dlgnorm
    write(*,'(//,10X," Minimization QNA_N1QN3 starts ...",/  &
             10x,"DXMIN =",G23.16,2X,"DF1 =",G23.16,2X,"EPSG =",G23.16  &
             /,10X,"IMPRES =",I3,2X,"NITER = ",I3,2X,"NSIM = ",I3)') zxmin,zdf1,zeps,impres,itermax,isimmax

    ! Do the minimization
    call utl_tmg_start(91,'----QuasiNewton')
    call qna_n1qn3(simvar, dscalqn, dcanonb, dcanab, nvadim, zhat,  &
                   zjsp, gradJ, zxmin, zdf1, zeps, impres, nulout, imode,   &
                   itermax,isimmax, iztrl, vatra, nmtra, intunused, rspunused,  &
                   zspunused)
    call tmg_stop(91)
    call fool_optimizer(obsSpaceData)

    write(*,'(//,20X,20("*"),2X    &
        ,/,20X,"              Minimization ended with MODE:",I4  &
        ,/,20X,"                Total number of iterations:",I4  &
        ,/,20X,"               Total number of simulations:",I4)' ) imode,itermax,isimmax

    if (mpi_myid == 0) write(*,*) 'end of fso_minimize'

    deallocate(vatra)
    deallocate(gradJ)

    call tmg_stop(90)

  end subroutine fso_minimize

  subroutine fso_sumFSO(obsSpaceData)
    implicit none
   
    real(8)            :: pfso_1
    type(struct_obs)   :: obsSpaceData
    integer            :: bodyIndex,itvs,isens,headerIndex
    integer            :: bodyIndexBeg, bodyIndexEnd 

    integer, parameter :: numFamily = 10
    character(len=2), parameter :: familyList(numFamily) = (/'UA','AI','SF','SC','TO','SW','PR','RO','GP','CH'/)
    real(8)            :: tfso(numFamily), tfsotov_sensors(tvs_nsensors),totFSO
    integer            :: numAss_local(numFamily), numAss_global(numFamily)
    integer            :: numAss_sensors_loc(tvs_nsensors), numAss_sensors_glb(tvs_nsensors)
    integer            :: ierr, familyIndex

    if (mpi_myid == 0) write(*,*) 'sum of FSO information' 
    
    ! initialize 
    do familyIndex = 1, numFamily
      tfso(familyIndex) = 0.d0
      numAss_local(familyIndex) = 0
      numAss_global(familyIndex) = 0
    end do

    tfsotov_sensors(:) = 0.d0
    numAss_sensors_loc(:) = 0
    numAss_sensors_glb(:) = 0
    totFSO = 0.d0

    do familyIndex = 1, numFamily
      do bodyIndex = 1, obs_numbody(obsSpaceData)

        pfso_1 = obs_bodyElem_r(obsSpaceData,OBS_FSO,bodyIndex)
        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == 1 ) then
          ! FSO for each family         
          if (obs_getFamily(obsSpaceData,bodyIndex=bodyIndex) == familyList(familyIndex) ) then
            tfso(familyIndex) = tfso(familyIndex) + pfso_1
            numAss_local(familyIndex) = numAss_local(familyIndex) + 1
          end if
        end if 

      end do
    end do

    do itvs = 1, tvs_nobtov
      headerIndex  = tvs_headerIndex(itvs)
      if (headerIndex > 0 ) then
        bodyIndexBeg = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
        bodyIndexEnd = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex) + bodyIndexBeg - 1
        do bodyIndex = bodyIndexBeg, bodyIndexEnd 
          pfso_1 = obs_bodyElem_r(obsSpaceData,OBS_FSO,bodyIndex)
          if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == 1 ) then
            isens = tvs_lsensor (itvs)
            tfsotov_sensors(isens) =  tfsotov_sensors(isens) + pfso_1
            numAss_sensors_loc(isens) = numAss_sensors_loc(isens) + 1
          end if
        end do
      end if
    end do

    do familyIndex = 1, numFamily
      call mpi_allreduce_sumreal8scalar(tfso(familyIndex),"GRID")
      totFSO = totFSO + tfso(familyIndex)
      call rpn_comm_allreduce(numAss_local(familyIndex), numAss_global(familyIndex) ,1,"MPI_INTEGER","MPI_SUM","GRID",ierr)
    end do
    
    do isens = 1, tvs_nsensors
      call mpi_allreduce_sumreal8scalar(tfsotov_sensors(isens),"GRID")
      call rpn_comm_allreduce(numAss_sensors_loc(isens), numAss_sensors_glb(isens) ,1,"MPI_INTEGER","MPI_SUM","GRID",ierr)
    end do


    if (mpi_myid == 0) then

      write(*,*) ' '
      write(*,'(a15,f15.8)') 'Total FSO=', totFSO 
      write(*,*) ' '

      do familyIndex = 1, numFamily
        write(*,'(a4,a2,a2,f15.8,a16,i10)') 'FSO-', familyList(familyIndex), '=', tfso(familyIndex),'  Count Number=', numAss_global(familyIndex) 
      end do
      write(*,*) ' '

      if (tvs_nsensors > 0) then
        write(*,'(1x,a)') 'For TOVS decomposition by sensor:'
        write(*,'(1x,a)') '#  plt sat ins    FSO'
        do isens = 1, tvs_nsensors
          write(*,'(i2,1x,a,1x,a,1x,i2,1x,f15.8,i10)') isens,inst_name(tvs_instruments(isens)), &
                platform_name(tvs_platforms(isens)),tvs_satellites(isens),tfsotov_sensors(isens), numAss_sensors_glb(isens)
        end do
        write(*,*) ' '
      end if

    end if

  end subroutine fso_sumFSO

  subroutine simvar(indic,nvadim,zhat,Jtotal,gradJ)
    implicit none
    ! Argument declarations
    integer :: nvadim ! Dimension of the control vector in forecast error coraviances space
    ! Value of indic
    ! Note: 1 and 4 are reserved values for call back from m1qn3.
    !       For direct calls use other value than 1 and 4.
    ! =1 No action taken; =4 Both J(u) and its gradient are computed.
    ! =2 Same as 4 (compute J and gradJ) but do not interrupt timer of the
    !    minimizer.
    ! =3 Compute Jo and gradJo only.
    integer :: indic
    real(8)  :: Jtotal ! Cost function of the Variational algorithm
    real(8), dimension(nvadim) :: gradJ ! Gradient of the Variational Cost funtion
    real(8), dimension(nvadim) :: zhat ! Control variable in forecast error covariances space
    !
    ! Purpose: Implement the Variational solver as described in
    ! Courtier, 1997, Dual formulation of four-dimentional variational assimilation,
    ! Q.J.R., pp2449-2461.
    !
    ! Author : Simon Pellerin *ARMA/MSC October 2005
    !          (Based on previous versions of evaljo.ftn, evaljg.ftn and evaljgns.ftn).
    !
    ! Local declaration
    real(8) :: ahat_vhat(nvadim)
    real(8) :: Jb, Jobs
    type(struct_gsv) :: statevector
    type(struct_vco), pointer :: vco_anl

    call tmg_stop(91)
    call tmg_stop(90)

    if (indic /= 1) then ! No action taken if indic == 1
      fso_nsim = fso_nsim + 1

      if (mpi_myid == 0) then
        write(*,*) 'Entering simvar for simulation ',fso_nsim
        write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
        call flush(6)
      end if

      ! note: vhat = B_t^T/2 hat(del x_t)
      ahat_vhat(1:nvadim_mpilocal) = zhat(1:nvadim_mpilocal) + vhat(1:nvadim_mpilocal)

      ! Computation of background term of cost function:
      Jb = dot_product(zhat(1:nvadim_mpilocal),zhat(1:nvadim_mpilocal))/2.d0
      call mpi_allreduce_sumreal8scalar(Jb,"GRID")

      vco_anl => col_getVco(columnTrlOnAnlIncLev_ptr)
      call gsv_allocate(statevector,tim_nstepobsinc, hco_anl, vco_anl, &
                        dataKind_opt=pre_incrReal, mpi_local_opt=.true.)

      call bmat_sqrtB(ahat_vhat,nvadim_mpilocal,statevector)

      call s2c_tl(statevector,column_ptr,columnTrlOnAnlIncLev_ptr,obsSpaceData_ptr)  ! put in column H_horiz dx
      call utl_tmg_start(10,'--Observations')
      call utl_tmg_start(18,'----ObsOper_TL')
      call oop_Htl(column_ptr,columnTrlOnAnlIncLev_ptr,obsSpaceData_ptr,fso_nsim)  ! Save as OBS_WORK: H_vert H_horiz dx = Hdx

      call rmat_RsqrtInverseAllObs(obsSpaceData_ptr,OBS_WORK,OBS_WORK)  ! Save as OBS_WORK : R**-1/2 (Hdx)
      call tmg_stop(18)
      call tmg_stop(10)

      call cfn_calcJo(obsSpaceData_ptr)  ! Store J-obs in OBS_JOBS : 1/2 * R**-1 (Hdx)**2

      Jobs = 0.d0
      call utl_tmg_start(90,'--Minimization')
      call utl_tmg_start(92,'----SumCostFunction')
      call cfn_sumJo(obsSpaceData_ptr,Jobs)
      Jtotal = Jb + Jobs
      if (indic == 3) then
        Jtotal = Jobs
        if (mpi_myid == 0) write(*,FMT='(6X,"SIMVAR:  JO = ",G23.16,6X)') Jobs
      else
        Jtotal = Jb + Jobs
        if (mpi_myid == 0) write(*,FMT='(6X,"SIMVAR:  Jb = ",G23.16,6X,"JO = ",G23.16,6X,"Jt = ",G23.16)') Jb,Jobs,Jtotal
      end if
      call tmg_stop(92)
      call tmg_stop(90)

      call utl_tmg_start(10,'--Observations')
      call rmat_RsqrtInverseAllObs(obsSpaceData_ptr,OBS_WORK,OBS_WORK)  ! Modify OBS_WORK : R**-1 (Hdx)
      call tmg_stop(10)

      call col_zero(column_ptr)

      call utl_tmg_start(10,'--Observations')
      call utl_tmg_start(19,'----ObsOper_AD')
      call oop_Had(column_ptr,columnTrlOnAnlIncLev_ptr,obsSpaceData_ptr)   ! Put in column : H_vert**T R**-1 (Hdx)
      call tmg_stop(19)
      call tmg_stop(10)

      call s2c_ad(statevector,column_ptr,columnTrlOnAnlIncLev_ptr,obsSpaceData_ptr)  ! Put in statevector H_horiz**T H_vert**T R**-1 (Hdx)

      gradJ(:) = 0.d0
      call bmat_sqrtBT(gradJ,nvadim_mpilocal,statevector)
      call gsv_deallocate(statevector)

      if (indic /= 3) then
        gradJ(1:nvadim_mpilocal) = zhat(1:nvadim_mpilocal) + gradJ(1:nvadim_mpilocal)
      end if
    end if

    call utl_tmg_start(90,'--Minimization')
    call utl_tmg_start(91,'----QuasiNewton')

    if (mpi_myid == 0) write(*,*) 'end of simvar'

  end subroutine simvar

  SUBROUTINE DSCALQN(KDIM,PX,PY,DDSC)
    ! DSCALQN: inner product in canonical space
    !
    ! Purpose: interface for the inner product to be used
    ! by the minimization subroutines N1QN3.
    !
    ! Arguments
    !     i : KDIM      : dimension of the vectors
    !     i : PX, PY    : vector for which <PX,PY> is being calculated
    !     o : DDSC      : result of the inner product
    IMPLICIT NONE

    INTEGER KDIM
    REAL(8) PX(KDIM), PY(KDIM)
    REAL(8) DDSC

    CALL PRSCAL(KDIM,PX,PY,DDSC)
    RETURN
  END SUBROUTINE DSCALQN

  SUBROUTINE PRSCAL(KDIM,PX,PY,DDSC)
    ! PRSCAL: inner product in canonical space
    !
    ! Author  : P. Gauthier *ARMA/AES  January 27, 1993
    ! Purpose: evaluation of the inner product used in the minimization
    !
    ! Arguments
    !     i : KDIM     : dimension of the vectors
    !     i : PX, PY   : vector for which <PX,PY> is being calculated
    !     o : DDSC     : result of the inner product
    !
    IMPLICIT NONE

    INTEGER KDIM, J
    REAL(8) PX(KDIM), PY(KDIM)
    REAL(8) DDSC

    DDSC = 0.D0

    do j=1,nvadim_mpilocal
      DDSC = DDSC + PX(J)*PY(J)
    end do

    call mpi_allreduce_sumreal8scalar(ddsc,"GRID")

    RETURN

  END SUBROUTINE PRSCAL

  SUBROUTINE DCANAB(KDIM,PY,PX)
    ! DCANAB  - Change of variable associated with the canonical inner product
    !
    ! Author    JM Belanger CMDA/SMC   May 2001
    ! Double precision version based on single precision CTCAB.
    ! Refered to  as dummy argument DTCAB by N1QN3 minimization
    ! package.
    !
    ! Purpose: to compute PX = L^-1 * Py with L related to the inner product
    ! <PX,PY> = PX^t  L^t  L PY
    ! (see the modulopt documentation aboutn DTCAB)
    !
    IMPLICIT NONE

    INTEGER KDIM
    REAL(8) PX(KDIM), PY(KDIM)

    INTEGER JDIM

    DO JDIM = 1, KDIM
      PX(JDIM) = PY(JDIM)
    END DO

    RETURN

  END SUBROUTINE DCANAB

  SUBROUTINE DCANONB(KDIM,PX,PY)
    ! DCANONB  - Change of variable associated with the canonical inner product
    !
    ! Author    JM Belanger CMDA/SMC  May 2001
    ! Double precision version based on single precision CANONB.
    ! Refered to as dummy argument DTONB by N1QN3 minimization
    ! package.
    !
    ! Purpose: to compute PY = L * PX with L related to the inner product
    ! <PX,PY> = PX^t  L^t  L PY
    !(see the modulopt documentation about DTONB)
    !
    IMPLICIT NONE
    INTEGER KDIM
    REAL(8) PX(KDIM), PY(KDIM)

    INTEGER JDIM

    DO JDIM = 1, KDIM
      PY(JDIM) = PX(JDIM)
    END DO

    RETURN
  END SUBROUTINE DCANONB

end program midas_obsimpact 
