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

!--------------------------------------------------------------------------
!!
!! *Purpose*: Main program for variational minimization and background check 
!!            (depending on the mode selected in the namelist).
!!
!--------------------------------------------------------------------------
program midas_var
  !
  ! **Purpose**: Main program for variational minimization and background check 
  ! (depending on the mode selected in the namelist).
  !
  use ramDisk_mod
  use utilities_mod
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
  use minimization_mod
  use innovation_mod
  use WindRotation_mod
  use minimization_mod
  use analysisGrid_mod
  use bmatrix_mod
  use tovs_nl_mod
  use obsErrors_mod
  use variableTransforms_mod
  use obsOperators_mod
  use multi_ir_bgck_mod
  use biasCorrection_mod
  use increment_mod
  use residual_mod
  use stateToColumn_mod

  implicit none

  integer :: istamp,exdb,exfin
  integer :: ierr,nconf

  type(struct_obs),       target :: obsSpaceData
  type(struct_columnData),target :: trlColumnOnAnlLev
  type(struct_columnData),target :: trlColumnOnTrlLev
  type(struct_vco),       pointer :: vco_anl
  type(struct_hco),       pointer :: hco_anl
  type(struct_gsv) :: statevector_incr

  real*8,allocatable :: controlVector_incr(:)

  character(len=9) :: clmsg
  character(len=48) :: obsMpiStrategy, varMode

  logical :: writeAnalysis
  NAMELIST /NAMCT0/NCONF,writeAnalysis

  integer :: nulnam, fnom, fclos, get_max_rss

  istamp = exdb('VAR','DEBUT','NON')

  write(*,'(/,' //                                                &
            '3(" *****************"),/,' //                       &
            '14x,"-- START OF MAIN PROGRAM MIDAS-VAR: --",/,' //   &
            '14x,"-- VARIATIONAL ASSIMILATION          --",/, ' //&
            '14x,"-- VAR Revision number   ",a," --",/,' //       &
            '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  ! MPI initilization
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

  select case(nconf)
  case (141)
    write(*,*)
    write(*,*) 'midas-var: Analysis mode selected'
    varMode='analysis'
  case (111)
    write(*,*)
    write(*,*) 'midas-var: Background check for IR sat. data mode selected'
    varMode='bgckIR'
  case (101)
    write(*,*)
    write(*,*) 'midas-var: Background check for conventional obs mode selected'
    varMode='bgckConv'
  case default
    write(*,*)
    write(*,*) 'midas-var: Unknown mode ', nconf
    call utl_abort('midas-var')
  end select

  call ram_setup

  ! 2. Decide on configuration of job

  ! ---BGCHECK (conventional obs)--- !
  if ( trim(varMode) == 'bgckConv' ) then
    if(mpi_myid == 0) write(*,*) 'MIDAS-VAR: CONVENTIONNAL BGCHECK MODE'

    ! Do initial set up
    call tmg_start(2,'PREMIN')

    obsMpiStrategy = 'LIKESPLITFILES'

    call var_setup('ALL') ! obsColumnMode   

    ! Reading, horizontal interpolation and unit conversions of the 3D trial fields
    call inn_setupBackgroundColumns(trlColumnOnTrlLev,obsSpaceData)

    ! Interpolate trial columns to analysis levels and setup for linearized H
    call inn_setupBackgroundColumnsAnl(trlColumnOnTrlLev,trlColumnOnAnlLev)

    ! Compute observation innovations and prepare obsSpaceData for minimization
    call inn_computeInnovation(trlColumnOnTrlLev,obsSpaceData)
    call tmg_stop(2)

    ! Do the background check and output the observation data files
    call bgcheck_conv(trlColumnOnAnlLev,trlColumnOnTrlLev,obsSpaceData)

  ! ---BGCHECK (AIRS, IASI, CrIS)--- !
  else if ( trim(varMode) == 'bgckIR' ) then
    if(mpi_myid == 0) write(*,*) 'MIDAS-VAR: HYPERSPECTRAL IR BGCHECK MODE'

    ! Do initial set up
    call tmg_start(2,'PREMIN')

    obsMpiStrategy = 'LIKESPLITFILES'

    call var_setup('ALL') ! obsColumnMode   

    ! Reading, horizontal interpolation and unit conversions of the 3D trial fields
    call inn_setupBackgroundColumns(trlColumnOnTrlLev,obsSpaceData)

    ! Compute observation innovations and prepare obsSpaceData for minimization
    call inn_computeInnovation(trlColumnOnTrlLev,obsSpaceData)
    call tmg_stop(2)

    ! Do the background check and output the observation data files
    call irbg_bgCheckIR(trlColumnOnTrlLev,obsSpaceData)

  ! ---ANALYSIS MODE--- !
  else if ( trim(varMode) == 'analysis' ) then
    write(*,*) 'MIDAS-VAR: ANALYSIS MODE'

    ! Do initial set up
    call tmg_start(2,'PREMIN')

    obsMpiStrategy = 'LATLONTILESBALANCED'

    call var_setup('VAR') ! obsColumnMode
    call tmg_stop(2)

    !
    !- Initialize the background-error covariance, also sets up control vector module (cvm)
    !
    call bmat_setup(hco_anl,vco_anl)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Initialize variational bias correction (default is to not use it)
    !
    call bias_setup()
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    ! - Initialize the gridded variable transform module
    !
    call vtr_setup(hco_anl,vco_anl)

    !
    !- Set up the minimization module, now that the required parameters are known
    !  NOTE: some global variables remain in minimization_mod that must be initialized before
    !        inn_setupBackgroundColumns
    !
    call min_setup( cvm_nvadim ) ! IN
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! Reading, horizontal interpolation and unit conversions of the 3D trial fields
    call tmg_start(2,'PREMIN')
    call inn_setupBackgroundColumns(trlColumnOnTrlLev,obsSpaceData)

    ! Interpolate trial columns to analysis levels and setup for linearized H
    call inn_setupBackgroundColumnsAnl(trlColumnOnTrlLev,trlColumnOnAnlLev)

    ! Compute observation innovations and prepare obsSpaceData for minimization
    call inn_computeInnovation(trlColumnOnTrlLev,obsSpaceData)
    call tmg_stop(2)

    allocate(controlVector_incr(cvm_nvadim),stat=ierr)
    if(ierr.ne.0) then
      write(*,*) 'var: Problem allocating memory for ''controlVector_incr''',ierr
      call utl_abort('aborting in VAR')
    endif

    ! Do minimization of cost function
    call min_minimize(trlColumnOnAnlLev,obsSpaceData,controlVector_incr)

    ! Compute satellite bias correction increment and write to file
    call bias_writebias(controlVector_incr,cvm_nvadim)

    call gsv_allocate(statevector_incr, tim_nstepobsinc, hco_anl, vco_anl, &
         datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true.)

    ! get final increment
    call inc_getIncrement(controlVector_incr,statevector_incr,cvm_nvadim)

    deallocate(controlVector_incr)

    ! output the analysis increment
    call tmg_start(6,'WRITEINCR')
    call inc_writeIncrement(statevector_incr) ! IN
    call tmg_stop(6)

    ! Conduct obs-space post-processing diagnostic tasks (some diagnostic
    ! computations controlled by NAMOSD namelist in flnml)
    call osd_ObsSpaceDiag(obsSpaceData,trlColumnOnAnlLev)

    ! Deallocate memory related to B matrices
    call bmat_finalize()

    ! compute and write the analysis (as well as the increment on the trial grid)
    if (writeAnalysis) then
      call tmg_start(129,'ADDINCREMENT')
      call inc_computeAndWriteAnalysis(statevector_incr) ! IN
      call tmg_stop(129)
    end if

    if (mpi_myid == 0) then
      clmsg = 'REBM_DONE'
      call utl_writeStatus(clmsg)
    end if

    ! calculate OBS_OMA for diagnostic (i.e. non-assimilated) observations
    call var_calcOmA(statevector_incr,trlColumnOnAnlLev,obsSpaceData,obsAssVal=3)

    call gsv_deallocate(statevector_incr)

    ! Deallocate memory related to variational bias correction
    call bias_finalize()

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
    type(struct_vco),pointer :: vco_trl => null()
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
    !- Set vertical coordinate parameters from !! record in trial file
    !
    if(mpi_myid.eq.0) write(*,*)''
    if(mpi_myid.eq.0) write(*,*)'var_setup: Set vcoord parameters for trial grid'
    call vco_SetupFromFile( vco_trl,     & ! OUT
                            './trlm_01')   ! IN
    call col_setVco(trlColumnOnTrlLev,vco_trl)

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
       !- Iniatilized the core (Non-Exteded) analysis grid
       call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore' ) ! IN
       !- Setup the LAM analysis grid metrics
       call agd_SetupFromHCO( hco_anl, hco_core ) ! IN
    end if

    if ( hco_anl % rotated ) then
       call uvr_Setup(hco_anl) ! IN 
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
    !- Setup observation operators
    !
    call oop_setup(varMode) ! IN

    !
    !- Basic setup of columnData module
    !
    call col_setup
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Memory allocation for background column data
    !
    call col_allocate(trlColumnOnAnlLev,obs_numheader(obsSpaceData),mpiLocal_opt=.true.)
    call col_allocate(trlColumnOnTrlLev,obs_numheader(obsSpaceData),mpiLocal_opt=.true.)

    !
    !- Initialize the observation error covariances
    !
    call oer_setObsErrors(obsSpaceData, varMode) ! IN
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  end subroutine var_setup

!--------------------------------------------------------------------------
!! *Purpose*: Calculates the OmA for diagnostic (i.e. valid but non-assimilated)
!!            observations.
!!
!!            The last argument, 'obsAssVal', contains the value of
!!            'OBS_ASS' to test against to know which observations are
!!            not assimilated.
!!
!! @author M. Sitwell Sept 2015
!--------------------------------------------------------------------------
  subroutine var_calcOmA(statevector_incr,columng,obsSpaceData,obsAssVal)

    implicit none

    type(struct_gsv), intent(inout) :: statevector_incr
    type(struct_columnData), intent(inout) :: columng
    type(struct_obs), intent(inout) :: obsSpaceData
    integer :: obsAssVal

    type(struct_columnData) :: column
    integer :: bodyIndex,headerIndex,ierr,diagnosticObsAssValue
    logical :: calc_OmA,calc_OmA_global

    ! Check for the presence of diagnostic only observations
    calc_OmA = .false.
    call obs_set_current_body_list(obsSpaceData)
    BODY: do
       bodyIndex = obs_getBodyIndex(obsSpaceData)
       if (bodyIndex < 0) exit BODY
       if (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex).eq.obsAssVal) then
          calc_OmA = .true.
          exit BODY
       end if
    end do BODY

    call rpn_comm_allreduce(calc_OmA,calc_OmA_global,1,"MPI_LOGICAL","MPI_LOR","GRID",ierr)
    if (.not.calc_OmA_global) return

    ! Initialize columnData object for increment
    call col_setVco(column,col_getVco(columng))
    call col_allocate(column,col_getNumCol(columng),mpiLocal_opt=.true.)
    call col_copyLatLon(columng,column)

    ! Put H_horiz dx in column
    call s2c_tl(statevector_incr,column,columng,obsSpaceData)

    ! Save as OBS_WORK: H_vert H_horiz dx = Hdx
    call oop_Htl(column,columng,obsSpaceData,min_nsim,obsAssVal_opt=obsAssVal)

    ! Calculate OBS_OMA from OBS_WORK : d-Hdx
    call res_compute(obsSpaceData,obsAssVal_opt=obsAssVal)

  end subroutine var_calcOmA

end program midas_var
