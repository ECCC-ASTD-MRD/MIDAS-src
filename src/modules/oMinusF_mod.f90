
module oMinusF_mod
  ! MODULE oMinusF_mod (prefix='omf' category='1. High-level functionality')
  !
  ! :Purpose: Module for Observation minus Forecast (O-F) computation
  !
  use codePrecision_mod
  use ramDisk_mod
  use midasMpi_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use timeCoord_mod
  use obsSpaceData_mod
  use columnData_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use obsFiles_mod
  use innovation_mod
  use obsErrors_mod
  use biasCorrectionConv_mod
  use obsSpaceErrorStdDev_mod
  use ensembleObservations_mod
  use ensembleStateVector_mod
  use fileNames_mod
  use statetocolumn_mod
  implicit none
  private

  ! public subroutines and functions
  public :: omf_oMinusF, omf_oMinusFens

  contains

    !--------------------------------------------------------------------------
    ! omf_oMinusF
    !--------------------------------------------------------------------------
    subroutine omf_oMinusF(columnTrlOnAnlIncLev, columnTrlOnTrlLev, obsSpaceData, &
                           midasMode, addHBHT, addSigmaO)
      !
      ! :Purpose: compute Observation-minus-Forecast (OmF)
      !
      implicit none

      ! Arguments:
      type(struct_columnData),target, intent(inout)  :: columnTrlOnAnlIncLev
      type(struct_columnData),target, intent(inout)  :: columnTrlOnTrlLev
      type(struct_obs),       target, intent(inout)  :: obsSpaceData
      type(struct_gsv)             :: stateVectorTrialHighRes
      character(len=*), intent(in) :: midasMode
      logical, intent(in) :: addHBHT
      logical, intent(in) :: addSigmaO

      ! Locals:
      type(struct_vco),       pointer :: vco_anl  => null()
      type(struct_hco),       pointer :: hco_anl  => null()
      type(struct_hco),       pointer :: hco_trl => null()
      type(struct_vco),       pointer :: vco_trl => null()
      type(struct_hco),       pointer :: hco_core => null()
      logical           :: allocHeightSfc
      character(len=48) :: obsMpiStrategy
      character(len=3)  :: obsColumnMode
      character(len=10) :: trialFileName
      integer :: dateStampFromObs, get_max_rss

      write(*,*) " ---------------------------------------"
      write(*,*) " ---  START OF SUBROUTINE oMinusF    ---"
      write(*,*) " ---  Computation of the innovation  ---"
      write(*,*) " ---------------------------------------"

      !
      !- 1.  Settings and module initializations
      !
      write(*,*)
      write(*,*) '> omf_oMinusF: setup - START'

      obsMpiStrategy = 'LIKESPLITFILES'
      obsColumnMode  = 'VAR'
      trialFileName  = './trlm_01'

      !- 1.3 RAM disk usage
      call ram_setup

      !- 1.4 Temporal grid and dateStamp from trial file
      call tim_setup(fileNameForDate_opt=trim(trialFileName))

      !- 1.5 Observation file names and get datestamp, but do not use it
      call obsf_setup(dateStampFromObs, trim(midasMode))

      !- 1.6 Constants
      if ( mmpi_myid == 0 ) then
        call mpc_printConstants(6)
        call pre_printPrecisions
      end if

      !- 1.7 Variables of the model states
      call gsv_setup

      !- 1.8 Set the horizontal domain
      if ( addHBHT ) then
        call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS') ! IN
        if ( hco_anl % global ) then
          hco_core => hco_anl
        else
          !- Iniatilized the core (Non-Exteded) analysis grid
          call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID') ! IN
        end if
      else
        call hco_SetupFromFile(hco_anl, trim(trialFileName), ' ') ! IN
        hco_core => hco_anl
      end if

      ! 1.9 Setup a column vector following the analysis vertical grid
      if ( addHBHT ) then
        call vco_SetupFromFile(vco_anl,        & ! OUT
                               './analysisgrid') ! IN
        call col_setVco(columnTrlOnAnlIncLev,vco_anl)
      end if

      !- 1.10 Setup and read observations
      call inn_setupObs(obsSpaceData, hco_anl, obsColumnMode, obsMpiStrategy,trim(midasMode))!IN

      ! Apply optional bias corrections when namelist logicals {fam}BiasActive are TRUE
      ! (Only reverse existing corrections when namelist logicals {fam}RevOnly are TRUE)
      if (obs_famExist(obsSpaceData,'AI')) call bcc_applyAIBcor(obsSpaceData)    
      if (obs_famExist(obsSpaceData,'GP')) call bcc_applyGPBcor(obsSpaceData)
      if (obs_famExist(obsSpaceData,'UA')) call bcc_applyUABcor(obsSpaceData)
      
      !- 1.11 Basic setup of columnData module
      call col_setup

      !- 1.12 Memory allocation for background column data
      if ( addHBHT ) then
        call col_allocate(columnTrlOnAnlIncLev, obs_numheader(obsSpaceData),mpiLocal_opt=.true.)
      end if

      if ( addSigmaO ) then
        !- 1.13 Initialize the observation error covariances
        write(*,*)
        write(*,*) '> omf_oMinusF: Adding sigma_O'
        call oer_setObsErrors(obsSpaceData, trim(midasMode))
      end if

      ! Reading trials
      call inn_getHcoVcoFromTrlmFile( hco_trl, vco_trl )
      allocHeightSfc = ( vco_trl%Vcode /= 0 )

      call gsv_allocate( stateVectorTrialHighRes, tim_nstepobs, hco_trl, vco_trl,  &
                         dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                         mpi_distribution_opt='Tiles', dataKind_opt=4,  &
                         allocHeightSfc_opt=allocHeightSfc, hInterpolateDegree_opt='LINEAR', &
                         beSilent_opt=.false. )
      call gsv_zero( stateVectorTrialHighRes )
      call gio_readTrials( stateVectorTrialHighRes )
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      ! Horizontally interpolate trials to trial columns
      call inn_setupColumnsOnTrlLev( columnTrlOnTrlLev, obsSpaceData, hco_core, &
                                       stateVectorTrialHighRes )

      write(*,*)
      write(*,*) '> omf_oMinusF: setup - END'
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      !
      !- 2.  O-F computation
      !

      !- 2.1 Compute observation innovations
      write(*,*)
      write(*,*) '> omf_oMinusF: compute innovation'
      call inn_computeInnovation(columnTrlOnTrlLev, obsSpaceData, analysisMode_opt=.false.)

      if ( addHBHT ) then
        write(*,*)
        write(*,*) '> omf_oMinusF: Adding HBH^T'
        !- 2.2 Interpolate background columns to analysis levels and setup for linearized H
        call inn_setupColumnsOnAnlIncLev( columnTrlOnTrlLev, columnTrlOnAnlIncLev )
        !- 2.3 Compute the background errors in observation space
        call ose_computeStddev(columnTrlOnAnlIncLev,hco_anl,obsSpaceData)
      end if

    end subroutine omf_oMinusF

    !--------------------------------------------------------------------------
    ! omf_oMinusFens
    !--------------------------------------------------------------------------
    subroutine omf_oMinusFens(ensObs, obsSpaceData, nEns, ensPathName, &
                              midasMode, addHBHT, addSigmaO)
      !
      ! :Purpose: compute Observation-minus-Forecast (OmF) for ensembles
      !
      implicit none

      ! Arguments:
      type(struct_eob), target, intent(inout)  :: ensObs
      type(struct_obs), target, intent(inout)  :: obsSpaceData
      integer, intent(in) :: nEns
      character(len=*), intent(in) :: ensPathName
      character(len=*), intent(in) :: midasMode
      logical, intent(in) :: addHBHT
      logical, intent(in) :: addSigmaO

      ! Locals:
      type(struct_columnData)   :: columTrlOnTrlLev
      type(struct_columnData)   :: columnTrlOnAnlIncLev
      type(struct_ens)          :: ensembleTrl4D
      type(struct_gsv)          :: stateVector4D
      type(struct_gsv)          :: stateVectorWithZandP4D
      type(struct_gsv)          :: stateVectorHeightSfc            
      type(struct_vco), pointer :: vco_anl  => null()
      type(struct_hco), pointer :: hco_anl  => null()
      type(struct_hco), pointer :: hco_ens  => null()
      type(struct_vco), pointer :: vco_ens  => null()
      type(struct_hco), pointer :: hco_core => null()
      character(len=256) :: ensFileName
      character(len=48)  :: obsMpiStrategy
      character(len=3)   :: obsColumnMode
      integer, allocatable :: dateStampList(:)      
      integer :: datestamp, get_max_rss, memberIndex

      write(*,*) " -----------------------------------------"
      write(*,*) " ---  START OF SUBROUTINE oMinusFens   ---"
      write(*,*) " ---  Computation of the innovation    ---"
      write(*,*) " -----------------------------------------"

      !
      !- 1.  Settings and module initializations
      !
      write(*,*)
      write(*,*) '> omf_oMinusFens: setup - START'

      obsMpiStrategy = 'LIKESPLITFILES'
      obsColumnMode  = 'VAR'
 
      !- 1.3 RAM disk usage
      call ram_setup

      !- 1.4 Horizontal, vertical and temporal dimensions
      call fln_ensFileName(ensFileName, ensPathName, memberIndex_opt=1)
      
      call tim_setup( fileNameForDate_opt=trim(ensFileName) )
      allocate(dateStampList(tim_nstepobs))
      call tim_getstamplist(dateStampList,tim_nstepobs,tim_getDatestamp())
      
      call hco_setupFromFile(hco_ens, ensFileName, ' ', 'ENSFILEGRID')
      call vco_setupFromFile(vco_ens, ensFileName)
      
      !- 1.5 Observation file names and get datestamp, but do not use it
      call obsf_setup(dateStamp, trim(midasMode))

      !- 1.6 Constants
      if ( mmpi_myid == 0 ) then
        call mpc_printConstants(6)
        call pre_printPrecisions
      end if

      !- 1.7 Variables of the model states
      call gsv_setup

      !- 1.8 Set the horizontal domain
      if ( addHBHT ) then
        call hco_SetupFromFile(hco_anl, './analysisgrid', 'ANALYSIS') ! IN
        if ( hco_anl % global ) then
          hco_core => hco_anl
        else
          !- Iniatilized the core (Non-Exteded) analysis grid
          call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID') ! IN
        end if
      else
        hco_core => hco_ens
      end if

      ! 1.9 Setup a column vector following the analysis vertical grid
      if ( addHBHT ) then
        call vco_SetupFromFile(vco_anl,        & ! OUT
                               './analysisgrid') ! IN
        call col_setVco(columnTrlOnAnlIncLev,vco_anl)
      end if

      !- 1.10 Setup and read observations
      call inn_setupObs(obsSpaceData, hco_core, obsColumnMode, obsMpiStrategy, trim(midasMode)) !IN

      ! Apply optional bias corrections when namelist logicals {fam}BiasActive are TRUE
      ! (Only reverse existing corrections when namelist logicals {fam}RevOnly are TRUE)
      if (obs_famExist(obsSpaceData,'AI')) call bcc_applyAIBcor(obsSpaceData)    
      if (obs_famExist(obsSpaceData,'GP')) call bcc_applyGPBcor(obsSpaceData)
      if (obs_famExist(obsSpaceData,'UA')) call bcc_applyUABcor(obsSpaceData)
      
      !- 1.11 Basic setup of columnData module
      call col_setup
      call col_setVco(columTrlOnTrlLev, vco_ens)
      call col_allocate(columTrlOnTrlLev, obs_numheader(obsSpaceData))
      write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'

      !- 1.12 Memory allocation for background column data
      if ( addHBHT ) then
        call col_allocate(columnTrlOnAnlIncLev, obs_numheader(obsSpaceData),mpiLocal_opt=.true.)
      end if

      if ( addSigmaO ) then
        !- 1.13 Initialize the observation error covariances
        write(*,*)
        write(*,*) '> omf_oMinusF: Adding sigma_O'
        call oer_setObsErrors(obsSpaceData, trim(midasMode))
      end if

      !- 1.14 Allocate and initialize eob object for storing y-HX values
      call eob_allocate(ensObs, nEns, obs_numBody(obsSpaceData), obsSpaceData)
      call eob_zero(ensObs)
      call eob_setLatLonObs(ensObs)

      !- 1.15 Allocate statevector for storing state with heights and pressures allocated (for s2c_nl)
      call gsv_allocate(stateVectorWithZandP4D, tim_nstepobs, hco_ens, vco_ens, &
                        dateStamp_opt=tim_getDateStamp(),  &
                        mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                        dataKind_opt=4, allocHeightSfc_opt=.true.)
      call gsv_zero(stateVectorWithZandP4D)

      call gsv_allocate(stateVector4D, tim_nstepobs, hco_ens, vco_ens, &
                        dateStamp_opt=tim_getDateStamp(),  &
                        mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                        dataKind_opt=4, allocHeightSfc_opt=.true., &
                        allocHeight_opt=.false., allocPressure_opt=.false.)
      call gsv_zero(stateVector4D)
      
      !- 1.16 Read the sfc height from ensemble member 1
      call gsv_allocate(stateVectorHeightSfc, 1, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                        mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                        dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0','TT'/))
      call gio_readFromFile(stateVectorHeightSfc, ensFileName, ' ', ' ',  &
                            containsFullField_opt=.true., readHeightSfc_opt=.true.)

      !- 1.17 Reading ensemble
      call ens_allocate(ensembleTrl4D, nEns, tim_nstepobs, hco_ens, vco_ens, dateStampList)
      call ens_readEnsemble(ensembleTrl4D, ensPathName, biPeriodic=.false.)
      
      write(*,*)
      write(*,*) '> omf_oMinusFens: setup - END'
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

      !
      !- 2.  O-F computation
      !
      do memberIndex = 1, nEns
        write(*,*) ''
        write(*,*) 'oMinusFens: compute O-P for ensemble member ', memberIndex
        write(*,*) 'Memory Used: ', get_max_rss()/1024, 'Mb'
        
        !- 2.1 Copy selected member to a stateVector
        call ens_copyMember(ensembleTrl4D, stateVector4D, memberIndex)
        call gsv_copy(stateVector4D, stateVectorWithZandP4D, allowVarMismatch_opt=.true., &
                      beSilent_opt=.true.)
        call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorWithZandP4D)

        !- 2.2 Compute and set Yb in ensObs
        call s2c_nl(stateVectorWithZandP4D, obsSpaceData, columTrlOnTrlLev, hco_ens, &
                   timeInterpType='LINEAR', dealloc_opt=.false., &
                   beSilent_opt=.true.)
        
        !- 2.3 Compute innovation
        call inn_computeInnovation(columTrlOnTrlLev, obsSpaceData, analysisMode_opt=.false., beSilent_opt=.true.)

        !- 2.4 Copy to ensObs
        call eob_setYb(ensObs, memberIndex)

        !- 2.5 Add background errors in observation space
        if (addHBHT .and. memberIndex == 1) then
          write(*,*)
          write(*,*) '> omf_oMinusFens: Adding HBH^T'
          !  Interpolate background columns to analysis levels and setup for linearized H
          call inn_setupColumnsOnAnlIncLev(columTrlOnTrlLev, columnTrlOnAnlIncLev)
          !  Compute the background errors in observation space
          call ose_computeStddev(columnTrlOnAnlIncLev,hco_anl,obsSpaceData)
        end if

      end do

      !
      !- 3.  Ending
      !
      call gsv_deallocate(stateVectorWithZandP4D)
      call gsv_deallocate(stateVector4D)
      call gsv_deallocate(stateVectorHeightSfc)

    end subroutine omf_oMinusFens
    
end module oMinusF_mod
