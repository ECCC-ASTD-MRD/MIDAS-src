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
program midas_gencoeff
  !
  ! **Purpose**: Main program for variational minimization and background check 
  ! (depending on the mode selected in the namelist).
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
  use obsFiles_mod
  use obsFilter_mod  
  use innovation_mod
  use tovs_nl_mod
  use obsErrors_mod
  use variableTransforms_mod
  use obsOperators_mod
  use statetocolumn_mod
  use biasCorrection_mod
  use increment_mod
  use stateToColumn_mod
  use backgroundCheck_mod
  use analysisGrid_mod

  implicit none

  integer,external :: exdb,exfin,fnom, fclos, get_max_rss
  integer ::  nulnam, headerIndex
  integer :: ierr,istamp

  type(struct_obs),        target  :: obsSpaceData
  type(struct_columnData), target  :: trlColumnOnTrlLev
  type(struct_columnData), target  :: trlColumnOnAnlLev
  type(struct_hco), pointer        :: hco_anl => null()
  type(struct_vco), pointer        :: vco_anl => null()


  character(len=48),parameter :: obsMpiStrategy = 'LIKESPLITFILES', &
                                 varMode        = 'analysis'



  istamp = exdb('GENCOEFF','DEBUT','NON')

  write(*,'(/,' //                                                &
            '3(" *****************"),/,' //                       &
            '14x,"-- START OF MAIN PROGRAM MIDAS-GENCOEFF: --",/,' //   &
            '14x,"-- BIAS CORRECTION COEFFICIENT COMPUTATION  --",/, ' //&
            '14x,"-- Revision : ",a," --",/,' //       &
            '3(" *****************"))') 'GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE'

  ! MPI initialization
  call mpi_initialize  

  call tmg_init(mpi_myid, 'TMG_GENCOEFF' )

  call tmg_start(1,'MAIN')

 
  ! 1. Top level setup

  nulnam=0
  ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  if(ierr.ne.0) call utl_abort('midas-var: Error opening file flnml')
  !read(nulnam,nml=namct0,iostat=ierr)
  if(ierr.ne.0) call utl_abort('midas-var: Error reading namelist')
  !write(*,nml=namct0)
  ierr=fclos(nulnam)

  call ram_setup()
 

  ! Do initial set up
  call tmg_start(2,'SETUP')

  call gencoeff_setup('VAR') ! obsColumnMode
  call tmg_stop(2)

  ! Read trials and horizontally interpolate to columns
  call tmg_start(3,'TRIALS')
  call inn_setupBackgroundColumns( trlColumnOnAnlLev, obsSpaceData )
  call tmg_stop(3)



  call bias_setup()
  !
  ! Remove bias correction if requested
  !
  call tmg_start(4,'REMOVE_BCOR')
  call bias_removeBiasCorrection(obsSpaceData,"TO")
  call tmg_stop(4)

   
  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Compute observation innovations
  call tmg_start(5,'COMP_INOV')
  call inn_computeInnovation(trlColumnOnAnlLev,obsSpaceData)
  call tmg_stop(5)

  
    !
    ! Refresh bias correction if requested
    !
  call tmg_start(6,'REFRESH_BCOR')
  call bias_refreshBiasCorrection(obsSpaceData,trlColumnOnAnlLev)
  call tmg_stop(6)

  !
  ! Filter obs if requested
  !

  call tmg_start(7,'REGRESSION')
  call bias_do_regression(trlColumnOnAnlLev,obsSpaceData)
  call tmg_stop(7)

  ! Write coefficients to file
  call tmg_start(8,'WRITECOEFFS')
  call bias_writebias()
  call tmg_stop(8)
  

  !
  ! output O-F statistics befor bias coorection
  !
  call tmg_start(9,'STATS')
  call bias_computeResidualsStatistics(obsSpaceData,"_raw")
  call tmg_stop(9)
  !
  ! fill OBS_BCOR with computed bias correction
  !
  call tmg_start(10,'COMPBIAS')
  call bias_calcBias(obsSpaceData,trlColumnOnAnlLev)
  call tmg_stop(10)

  !
  ! output  O-F statistics after bias coorection
  !
  call tmg_start(9,'STATS')
  call  bias_computeResidualsStatistics(obsSpaceData,"_corrected")
  call tmg_stop(9)

  ! Deallocate internal bias correction structures 

  call bias_finalize()

   ! Deallocate copied obsSpaceData
   call obs_finalize(obsSpaceData)

  
  ! 3. Job termination

   istamp = exfin('GENCOEFF','FIN','NON')

   call tmg_stop(1)

   call tmg_terminate(mpi_myid, 'TMG_GENCOEFF' )

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
  subroutine gencoeff_setup(obsColumnMode)
    implicit none

    character (len=*) :: obsColumnMode
    integer :: datestamp
    type(struct_hco),pointer :: hco_core => null()

    write(*,*) ''
    write(*,*) '----------------------------------------'
    write(*,*) '-- Starting subroutine gencoeff_setup --'
    write(*,*) '----------------------------------------'

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
    if(mpi_myid.eq.0) write(*,*)'gencoeff_setup: Set hco parameters for analysis grid'
    call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

    if ( hco_anl % global ) then
      call agd_SetupFromHCO( hco_anl ) ! IN
    else
      !- Initialize the core (Non-Extended) analysis grid
      if(mpi_myid.eq.0) write(*,*)'gencoeff_setup: Set hco parameters for core grid'
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
    !- Setup observation operators
    !
    call oop_setup(varMode) ! IN

    !
    !- Basic setup of columnData module
    !
    call col_setup
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Initialize the observation error covariances
    !
    call oer_setObsErrors(obsSpaceData, varMode) ! IN
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

 end subroutine gencoeff_setup



end program midas_gencoeff
