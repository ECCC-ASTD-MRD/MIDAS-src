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

program midas_diagHBHt
  !
  ! :Purpose: Main program for computing background error variance in observation
  !           space.
  !
  use codePrecision_mod
  use ramDisk_mod
  use utilities_mod
  use mpi_mod
  use mpivar_mod
  use MathPhysConstants_mod
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
  use randomnumber_mod
  use obsTimeInterp_mod
  use stateToColumn_mod
  use innovation_mod
  use analysisGrid_mod
  use bmatrix_mod
  use tovs_nl_mod
  use obsErrors_mod
  use gridVariableTransforms_mod
  use obsOperators_mod
  implicit none

  integer :: istamp,exdb,exfin
  integer :: ierr

  type(struct_obs),       target :: obsSpaceData
  type(struct_columnData),target :: trlColumnOnAnlLev
  type(struct_columnData),target :: trlColumnOnTrlLev

  character(len=48) :: obsMpiStrategy, varMode

  istamp = exdb('diagHBHt','DEBUT','NON')

  write(*,'(/,' //                                                &
            '3(" *****************"),/,' //                       &
            '14x,"-- START OF MAIN PROGRAM MIDAS-diagHBHt: --",/,' //   &
            '14x,"-- RANDOMIZED DIAGNOSTIC of HBHt         --",/,' //   &
            '14x,"-- VAR Revision number   ",a," --",/,' //       &
            '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  ! MPI initilization
  call mpi_initialize  

  call tmg_init(mpi_myid, 'TMG_VAR' )

  call tmg_start(1,'MAIN')

  varMode='analysis'

  call ram_setup

  ! Do initial set up
  call tmg_start(2,'PREMIN')

  obsMpiStrategy = 'LIKESPLITFILES'

  call var_setup('VAR') ! obsColumnMode
  call tmg_stop(2)

  ! Reading, horizontal interpolation and unit conversions of the 3D trial fields
  call tmg_start(2,'PREMIN')
  call inn_setupBackgroundColumns( trlColumnOnTrlLev, obsSpaceData )

  ! Interpolate trial columns to analysis levels and setup for linearized H
  call inn_setupBackgroundColumnsAnl(trlColumnOnTrlLev,trlColumnOnAnlLev)

  ! Compute observation innovations and prepare obsSpaceData for minimization
  call inn_computeInnovation(trlColumnOnTrlLev,obsSpaceData)
  call tmg_stop(2)

  ! Compute perturbed
  call diagHBHt(trlColumnOnTrlLev,trlColumnOnAnlLev,obsSpaceData)
  ! Conduct obs-space post-processing diagnostic tasks (some diagnostic 
  ! computations controlled by NAMOSD namelist in flnml)
  call osd_ObsSpaceDiag(obsSpaceData,trlColumnOnAnlLev)

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

  ! 3. Job termination

  istamp = exfin('diagHBHt','FIN','NON')

  call tmg_stop(1)

  call tmg_terminate(mpi_myid, 'TMG_VAR' )

  call rpn_comm_finalize(ierr) 

contains

  !--------------------------------------------------------------------------
  !! *Purpose*: Control of the preprocessing of the variational assimilation
  !!
  !! Revisions:
  !!           Y.J. Rochon, Jan 2016
  !!           - Addition of test on availability of input trial fields according
  !!             to related observation families.
  !--------------------------------------------------------------------------
  subroutine var_setup(obsColumnMode)
    implicit none

    character (len=*) :: obsColumnMode
    integer :: datestamp
    type(struct_vco),pointer :: vco_anl => null()
    type(struct_hco),pointer :: hco_anl => null()
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
    !- Initialize burp file names and set datestamp
    !
    call obsf_setup( dateStamp, 'analysis' )
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
    if(mpi_myid.eq.0) write(*,*)' preproc: Set hco parameters for analysis grid'
    call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

    if ( hco_anl % global ) then
      call agd_SetupFromHCO( hco_anl ) ! IN
    else
      !- Iniatilized the core (Non-Exteded) analysis grid
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

    !
    !- Initialize the background-error covariance, also sets up control vector module (cvm)
    !
    call bmat_setup(hco_anl,vco_anl)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    ! - Initialize the gridded variable transform module
    !
   
    call gvt_setup(hco_anl,vco_anl)
    
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine var_setup

  subroutine diagHBHt(columng,column,obsSpaceData)
    implicit none
    type(struct_gsv)        :: statevector

    type(struct_obs)        :: obsSpaceData ! Observation-related data
    type(struct_columnData) :: columng      ! Columns of the background interpolated 
                                          ! to analysis levels and to obs horizontal locations
    type(struct_columnData) :: column

    type(struct_hco), pointer :: hco_anl
    type(struct_vco), pointer :: vco_anl
    real(8) ,allocatable :: random_vector(:)
    real(8) ,allocatable :: local_random_vector(:)
    integer :: index_body, local_dimension, jj, ierr, dateprnt,timeprnt,nrandseed, istat
    integer ,external :: newdate,get_max_rss
    real(8) ,external :: gasdev
    !
    !- 1.  Initialization

    write(*,*)
    write(*,*) 'Computing perturbations for randomized HBHT evaluation START'

    vco_anl => col_getVco(columng)
    hco_anl => agd_getHco('ComputationalGrid')
    !- 1.3 Create a gridstatevector to store the perturbations
    call gsv_allocate(statevector,tim_nstepobsinc,hco_anl,vco_anl, &
                      dataKind_opt=pre_incrReal,mpi_local_opt=.true.)

    !- 1.4 Create column vectors to store the perturbation interpolated to obs horizontal locations
    call col_setVco(column,vco_anl)
    call col_allocate(column,col_getNumCol(columng),mpiLocal_opt=.true.)

    !- 1.6
    call oti_timeBinning(obsSpaceData,tim_nstepobsinc)

    !
    !- 2.  Compute the perturbations
    !

    !- 2.1 Random perturbations
    write(*,*)
    write(*,*) 'Generating random perturbation:'

    !- Global vector (same for each processors)
    allocate(random_vector(cvm_nvadim_mpiglobal),stat =istat )
    allocate(local_random_vector(cvm_nvadim),stat =istat )

    !- Initialize random number generator
    ierr = newdate(tim_getDatestamp(), dateprnt, timeprnt, -3)
    nrandseed=100*dateprnt + int(timeprnt/100.0) 
    write(*,*) 'diagHBHt: Random seed set to ',nrandseed ; call flush(6)
    call rng_setup(nrandseed)
    ! Generate a random vector from N(0,1)
    do jj = 1, cvm_nvadim_mpiglobal
      random_vector(jj) = rng_gaussian()
    enddo
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    !- Extract only the subvector for this processor
    call bmat_reduceToMPILocal(local_random_vector,  & ! OUT
         random_vector)      ! IN
    local_dimension = size(local_random_vector)
    !- Transform to control variables in physical space
    call bmat_sqrtB(local_random_vector,local_dimension,statevector)
    !- 2.2 Interpolation to the observation horizontal locations

    call s2c_tl( statevector,           & ! IN
                 column,                & ! OUT (H_horiz EnsPert)
                 columng, obsSpaceData )  ! IN
    !- 2.3 Interpolation to observation space
    call oop_Htl(column,columng,obsSpaceData,min_nsim=1)

  !- Copy from OBS_WORK to OBS_HPHT

    do index_body = 1, obs_numBody(obsSpaceData)
      call obs_bodySet_r(obsSpaceData,OBS_HPHT,index_body,obs_bodyElem_r(obsSpaceData,OBS_WORK,index_body) )
    end do
  !
  !- 3.  Ending/Deallocation
  !
    deallocate(random_vector, stat=istat)
    deallocate(local_random_vector, stat=istat)
    call col_deallocate(column)
    call gsv_deallocate(statevector)
    write(*,*)
    write(*,*) 'Computing perturbations for randomized HBHT evaluation END' ; call flush(6)

end subroutine diagHBHt




end program midas_diagHBHt
