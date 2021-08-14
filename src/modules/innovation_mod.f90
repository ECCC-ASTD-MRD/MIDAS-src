!--------------------------------------- LICENCE BEGIN -----------------------------------
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

module innovation_mod
  ! MODULE innovation_mod (prefix='inn' category='1. High-level functionality')
  !
  ! :Purpose: Several high-level subroutines used to compute the innovations:
  !           that is, the observation-minus-background values. This includes
  !           the subroutine that reads in the gridded high-res background state
  !           from standard files.
  !
  use codePrecision_mod
  use mpi_mod
  use ramDisk_mod
  use obsSpaceData_mod
  use columnData_mod
  use timeCoord_mod
  use obsTimeInterp_mod
  use obsOperators_mod
  use earthConstants_mod
  use mathPhysConstants_mod
  use mpivar_mod
  use horizontalCoord_mod
  use varNameList_mod
  use verticalCoord_mod
  use gridStateVector_mod
  use tt2phi_mod
  use utilities_mod
  use obsFilter_mod  
  use gps_mod
  use tovs_nl_mod
  use tovs_lin_mod
  use multi_ir_bgck_mod
  use obsFiles_mod
  use randomNumber_mod
  use obsErrors_mod
  use bufr_mod
  use statetocolumn_mod
  use biascorrectionSat_mod
  implicit none
  save
  private

  ! public procedures
  public :: inn_setupObs, inn_computeInnovation
  public :: inn_perturbObs, inn_setupBackgroundColumns, inn_setupBackgroundColumnsAnl

  character(len=48) :: innovationMode

contains

  !--------------------------------------------------------------------------
  ! inn_setupObs
  !--------------------------------------------------------------------------
  subroutine inn_setupobs(obsSpaceData, hco_anl, obsColumnMode, obsMpiStrategy, &
       innovationMode_in, obsClean_opt )
    !
    !:Purpose: To initialize the observation parameters and constants
    implicit none

    ! Arguments:
    type(struct_obs)                        :: obsSpaceData
    type(struct_hco), pointer               :: hco_anl
    character(len=*)                        :: obsMpiStrategy
    character(len=*)                        :: obsColumnMode
    character(len=*), intent(in)            :: innovationMode_in
    logical,                       optional :: obsClean_opt

    ! Locals:
    character(len=20) :: nameDimFile
    integer :: get_max_rss, ierr, fnom, fclos, unitDimFile, mxstn, mxobs
    logical :: obsDimFileExists

    write(*,FMT=9000)
9000 FORMAT(/,1x,' INN_SETUPOBS - Initialisation of observations',/,1x,3('- -----------'))

    !
    !- Setup de the mode
    !
    innovationMode = innovationMode_in

    !
    !- Specify the active observation-array columns
    !
    call obs_class_initialize(obsColumnMode)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    !
    !- Allocate memory for observation arrays
    !
    nameDimFile = './obs/cmadim'
    inquire( file=trim(nameDimFile), exist=obsDimFileExists )
    if ( obsDimFileExists ) then
      unitDimFile = 0
      ierr = fnom( unitDimFile, trim(nameDimFile), 'FTN+SEQ+R/O', 0 )
      read(unitDimFile,*) mxstn
      read(unitDimFile,*) mxobs
      ierr=fclos(unitDimFile)  
      call obs_initialize( obsSpaceData, numHeader_max=mxstn, numBody_max=mxobs, mpi_local=obsf_filesSplit() )
    else
      call obs_initialize( obsSpaceData, mpi_local=obsf_filesSplit() )
    end if
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Set up the list of elements to be assimilated and flags for rejection
    !
    call filt_setup(innovationMode) ! IN

    !
    !- Initialize TOVS processing
    !
    call tvs_setup

    !
    !- Read the observations from files
    !
    call tmg_start(11,'READ_OBS')
    call obsf_readFiles( obsSpaceData )
    call tmg_stop(11)

    !
    !- Initialize GPS processing
    !
    if (obs_famExist(obsSpaceData,'RO')) call gps_setupro
    if (obs_famExist(obsSpaceData,'GP')) call gps_setupgb
    
    !
    !- Filter out data from the obs data base
    !
    call filt_suprep(obsSpaceData)

    ! Additional filtering for bias correction if requested 
    call bcs_setup()
    call bcs_filterObs(obsSpaceData)

    if ( present(obsClean_opt) ) then
      if ( obsClean_opt ) then
        write(*,*) ''
        write(*,*) 'inn_setupObs: !!WARNING!! Performing a cleanup of obsSpaceData'
        write(*,*) '              This may make it impossible to update burp files'
        write(*,*) ''
        call obs_clean2(obsSpaceData)
      end if
    end if

    !
    !- Set (OBS_IPC and OBS_IPT) or OBS_IP columns according to the chosen strategy
    !
    write(*,*)
    write(*,*) 'inn_setupObs - Using obsMpiStrategy = ', trim(obsMpiStrategy)
    call setObsMpiStrategy(obsSpaceData,hco_anl,obsMpiStrategy)

    !
    !- Check if burp files already split
    !
    if ( obsf_filesSplit() ) then 
      ! local observations files, so just do reallocation to reduce memory used
      call obs_squeeze(obsSpaceData)
      if ( obs_columnActive_IH(obsSpaceData,OBS_IPC) ) then
        call obs_MpiRedistribute(obsSpaceData,OBS_IPC)
      end if
      write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    else
      ! complete set of obs on each MPI process, only keep subset according to OBS_IP
      call obs_reduceToMpiLocal(obsSpaceData)
    end if

    !
    !- Initialization and memory allocation for TOVS processing
    !
    if (obs_famExist(obsSpaceData,'TO') .and. (trim(innovationMode) /= 'thinning')) then
      call tvs_setupAlloc(obsSpaceData)
      if (trim(innovationMode) == 'bgck' ) call irbg_setup()
      ! Initialize non diagonal observation error matrices
      if ( trim(innovationMode) == 'analysis' .or. trim(innovationMode) == 'FSO') call oer_setInterchanCorr()
    end if

  end subroutine inn_setupobs
  
  !--------------------------------------------------------------------------
  ! inn_setupBackgroundColumns
  !--------------------------------------------------------------------------
  subroutine inn_setupBackgroundColumns(columnhr, obsSpaceData, hco_core, stateVectorTrialOut_opt)
    !
    !:Purpose: To compute vertical (and potentially slanted) columns of trial data interpolated to obs location
    implicit none
    
    ! arguments
    type(struct_columnData)    :: columnhr
    type(struct_obs)           :: obsSpaceData
    type(struct_hco), pointer  :: hco_core
    type(struct_gsv), optional :: stateVectorTrialOut_opt

    ! locals
    type(struct_gsv)          :: stateVectorTrial
    type(struct_gsv)          :: stateVectorTrialNoZorP 
    type(struct_hco), pointer :: hco_trl => null()
    type(struct_vco), pointer :: vco_trl => null()
    integer                   :: ierr, nulnam, fnom, fclos
    logical                   :: deallocInterpInfo, allocHeightSfc
    real(8), pointer          :: onecolumn(:)

    character(len=4), pointer :: anlVar(:)

    character(len=20) :: timeInterpType_nl  ! 'NEAREST' or 'LINEAR'
    integer           :: numObsBatches      ! number of batches for calling interp setup

    NAMELIST /NAMINN/timeInterpType_nl, numObsBatches

    write(*,*)
    write(*,*) 'inn_setupBackgroundColumns: START'
    nullify(hco_trl,vco_trl)

    timeInterpType_nl = 'NEAREST'
    numObsBatches     = 20

    if (utl_isNamelistPresent('naminn','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      if (ierr /= 0) call utl_abort('inn_setupBackgroundColumns: Error opening file flnml')
      read(nulnam,nml=naminn,iostat=ierr)
      if (ierr /= 0) call utl_abort('inn_setupBackgroundColumns: Error reading namelist')
      if (mpi_myid == 0) write(*,nml=naminn)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'inn_setupBackgroundColumns: Namelist block NAMINN is missing in the namelist.'
      write(*,*) '                            The default values will be taken.'
      if (mpi_myid == 0) write(*,nml=naminn)
    end if

    call tmg_start(10,'INN_SETUPBACKGROUNDCOLUMNS')

    ! check if gsv is initialized.
    if ( .not. gsv_isInitialized() ) then
       call utl_abort('inn_setupBackgroundColumns: add call to gsv_setup in the main program.')
    end if

    nullify(anlVar)
    call gsv_varNamesList(anlVar)
    call hco_SetupFromFile(hco_trl, './trlm_01', ' ', 'Trial', varName_opt=anlVar(1))

    call vco_SetupFromFile(vco_trl, './trlm_01')

    call col_setVco(columnhr,vco_trl)
    call col_allocate(columnhr,obs_numHeader(obsSpaceData),mpiLocal_opt=.true.)

    ! copy latitude from obsSpaceData
    if ( obs_numHeader(obsSpaceData) > 0 ) then
      call obs_extractObsRealHeaderColumn(columnhr%lat(:), obsSpaceData, OBS_LAT)
    end if

    if (vco_trl%Vcode == 0) then
      allocHeightSfc = .false.
    else
      allocHeightSfc = .true.
    end if

    deallocInterpInfo = .true.

    call gsv_allocate( stateVectorTrial, tim_nstepobs, hco_trl, vco_trl,  &
                       dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                       mpi_distribution_opt='Tiles', dataKind_opt=4,  &
                       allocHeightSfc_opt=allocHeightSfc, hInterpolateDegree_opt='LINEAR', &
                       beSilent_opt=.false. )
    call gsv_zero( stateVectorTrial )

    if ( gsv_varExist(stateVectorTrial,'Z_T') .or. &
         gsv_varExist(stateVectorTrial,'Z_M') .or. &
         gsv_varExist(stateVectorTrial,'P_T') .or. &
         gsv_varExist(stateVectorTrial,'P_M') ) then

      ! Use statevector without Z/P allocated to read trials 
      call gsv_allocate( stateVectorTrialNoZorP, tim_nstepobs, hco_trl, vco_trl,  &
                         dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                         mpi_distribution_opt='Tiles', dataKind_opt=4,  &
                         allocHeightSfc_opt=allocHeightSfc, hInterpolateDegree_opt='LINEAR', &
                         allocHeight_opt=.false., allocPressure_opt=.false., &
                         beSilent_opt=.false. )
      call gsv_zero( stateVectorTrialNoZorP )
      call gsv_readTrials( stateVectorTrialNoZorP )
      call gsv_copy( stateVectorTrialNoZorP, stateVectorTrial, allowVarMismatch_opt=.true. )
      call gsv_deallocate( stateVectorTrialNoZorP )
    else
      call gsv_readTrials( stateVectorTrial )
    end if

    ! if requested, make trials available to calling routine after degrading timesteps
    if (present(stateVectorTrialOut_opt)) then
      call gsv_allocate( stateVectorTrialOut_opt, tim_nstepobsinc, hco_trl, vco_trl,  &
                         dataKind_opt=pre_incrReal, &
                         dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                         allocHeightSfc_opt=allocHeightSfc, hInterpolateDegree_opt='LINEAR', &
                         allocHeight_opt=.false., allocPressure_opt=.false. )
      call gsv_copy( stateVectorTrial, stateVectorTrialOut_opt,  &
                     allowTimeMismatch_opt=.true., allowVarMismatch_opt=.true. )
    end if

    call s2c_nl( stateVectorTrial, obsSpaceData, columnhr, hco_core, &
                 timeInterpType=timeInterpType_nl, &
                 moveObsAtPole_opt=.true., numObsBatches_opt=numObsBatches, &
                 dealloc_opt=deallocInterpInfo )
    call gsv_deallocate(stateVectorTrial)

    if ( col_getNumCol(columnhr) > 0 .and. col_varExist(columnhr,'Z_T ') ) then
      write(*,*) 'inn_setupBackgroundColumns, statevector->Column 1:'
      write(*,*) 'Z_T:'
      onecolumn => col_getColumn(columnhr,1,'Z_T ')
      write(*,*) onecolumn(:)
      write(*,*) 'Z_M:'
      onecolumn => col_getColumn(columnhr,1,'Z_M ')
      write(*,*) onecolumn(:)

      nullify(onecolumn)
    end if
    if ( col_getNumCol(columnhr) > 0 .and. col_varExist(columnhr,'P_T ') ) then
      write(*,*) 'inn_setupBackgroundColumns, statevector->Column 1:'
      write(*,*) 'P_T:'
      onecolumn => col_getColumn(columnhr,1,'P_T ')
      write(*,*) onecolumn(:)
      write(*,*) 'P_M:'
      onecolumn => col_getColumn(columnhr,1,'P_M ')
      write(*,*) onecolumn(:)

      nullify(onecolumn)
    end if

    call tmg_stop(10)

    write(*,*) 'inn_setupBackgroundColumns: END'

  end subroutine inn_setupBackgroundColumns

  !--------------------------------------------------------------------------
  ! inn_setupBackgroundColumnsAnl
  !--------------------------------------------------------------------------
  subroutine inn_setupBackgroundColumnsAnl(columnhr,columng)
    !
    !:Purpose: To create trial data columns on analysis increment levels
    implicit none

    ! arguments
    type(struct_columnData) :: columng, columnhr

    ! locals
    integer :: jvar, jlev, columnIndex
    real(8), pointer :: columng_ptr(:), columnhr_ptr(:)

    write(*,*)
    write(*,*) 'inn_setupBackgroundColumnsAnl: START'

    call tmg_start(10,'INN_SETUPBACKGROUNDCOLUMNS')

    !
    !- Data copying from columnh to columng
    !
    
    ! copy latitude
    if ( col_getNumCol(columng) > 0 ) then
      columng%lat(:) = columnhr%lat(:)
    end if

    ! copy 2D surface variables
    do jvar = 1, vnl_numvarmax2D
      if ( .not. col_varExist(columng,vnl_varNameList2D(jvar)) ) cycle
      if ( col_getNumCol(columng) > 0 ) then       
        do columnIndex = 1, col_getNumCol(columng)
          columng_ptr  => col_getColumn( columng , columnIndex, vnl_varNameList2D(jvar) )
          columnhr_ptr => col_getColumn( columnhr, columnIndex, vnl_varNameList2D(jvar) )
          columng_ptr(:) = columnhr_ptr(:)
        end do
      end if
    end do

    !
    !- Calculate pressure profiles on analysis increment levels
    !
    if ( col_getNumCol(columng) > 0 .and. col_varExist(columng,'P_T') ) then
      call col_calcPressure(columng)

      ! Print pressure on thermo levels for the first original and destination column
      if ( mpi_myid == 0 ) then
        write(*,*) 'inn_setupBackgroundColumnsAnl, before vintprof, COLUMNHR(1):'
        write(*,*) 'P_T:'
        columnhr_ptr => col_getColumn(columnhr,1,'P_T')
        write(*,*) columnhr_ptr (:)

        write(*,*) 'inn_setupBackgroundColumnsAnl, before vintprof, COLUMNG(1):'
        write(*,*) 'P_T:'
        columng_ptr => col_getColumn(columng,1,'P_T')
        write(*,*) columng_ptr (:)
        write(*,*)
      end if
      
    end if

    !
    !- Vertical interpolation of 3D variables from trials levels to analysis increment levels
    !
    do jvar = 1, vnl_numvarmax3D

      if ( .not. col_varExist(columng,vnl_varNameList3D(jvar)) ) cycle
      
      call col_vintprof( columnhr, columng, vnl_varNameList3D(jvar), useColumnPressure_opt=.false. )

      if ( vnl_varNameList3D(jvar) == 'HU  ') then
        ! Imposing a minimum value for HU
        do columnIndex = 1, col_getNumCol(columng)
          columng_ptr => col_getColumn(columng,columnIndex,'HU')
          do jlev=1,col_getNumLev(columng,'TH')
            columng_ptr(jlev) = max(columng_ptr(jlev),col_rhumin)
          end do
        end do

      ! Imposing minimum value for LWCR at surface
      else if ( vnl_varNameList3D(jvar) == 'LWCR') then
        do columnIndex = 1, col_getNumCol(columng)
          columng_ptr => col_getColumn(columng,columnIndex,'LWCR')
          jlev = col_getNumLev(columng,'TH')
          columng_ptr(jlev) = max(columng_ptr(jlev),col_minClwAtSfc)
        end do

      else if (trim(vnl_varKindFromVarname(vnl_varNameList3D(jvar))) == 'CH') then
        ! Imposing boundary values for CH kind variables. This is to prevent
        ! undesired values usually from vertical extrapolation.
        do columnIndex = 1, col_getNumCol(columng)
          columng_ptr => col_getColumn(columng,columnIndex,trim(vnl_varNameList3D(jvar)))
          if ( col_minValVarKindCH(vnl_varListIndex(vnl_varNameList3D(jvar))) > 1.01*MPC_missingValue_R8 ) then
            do jlev=1,col_getNumLev(columng,'TH')
              columng_ptr(jlev) = max(columng_ptr(jlev),col_minValVarKindCH(vnl_varListIndex(vnl_varNameList3D(jvar))))
            end do
          end if
        end do
      end if

    end do

    ! Print pressure on thermo levels for the first column
    if ( col_getNumCol(columng) > 0 .and. col_varExist(columng,'P_T') ) then
      if ( mpi_myid == 0 ) then
        write(*,*) 'inn_setupBackgroundColumnsAnl, after vintprof, COLUMNG(1):'
        write(*,*) 'P_T:'
        columng_ptr => col_getColumn(columng,1,'P_T')
        write(*,*) columng_ptr (:)
        write(*,*)
      end if
    end if

    !
    !- Height adjustments
    !

    ! Set surface height
    do columnIndex = 1, col_getNumCol(columng)
      columng%HeightSfc(1,columnIndex) = columnhr%HeightSfc(1,columnIndex)
    end do

    ! Remove the height offset for the diagnostic levels for backward compatibility only
    if ( col_varExist(columng,'Z_T') .and. .not.columng%addHeightSfcOffset ) then 
      do columnIndex = 1, col_getNumCol(columng)
        columng_ptr => col_getColumn(columng,columnIndex,'Z_T')
        columng_ptr(col_getNumLev(columng,'TH')) = columng%HeightSfc(1,columnIndex)
      end do
    end if
    if ( col_varExist(columng,'Z_M') .and. .not.columng%addHeightSfcOffset ) then
      do columnIndex = 1, col_getNumCol(columng)
        columng_ptr => col_getColumn(columng,columnIndex,'Z_M')
        columng_ptr(col_getNumLev(columng,'MM')) = columng%HeightSfc(1,columnIndex)
      end do
    end if

    ! Print height info of the first original and interpolated columns
    if (col_getNumCol(columng) > 0) then
      write(*,*)
      write(*,*) 'inn_setupBackgroundColumnsAnl, vIntProf output:'

      if ( col_getNumLev(columng,'TH') > 0 .and. col_varExist(columng,'Z_T') ) then
        write(*,*) 'Z_T (columnhr):'
        columng_ptr => col_getColumn(columnhr,1,'Z_T')
        write(*,*) columng_ptr(:)
        write(*,*) 'Z_T (columng):'
        columng_ptr => col_getColumn(columng,1,'Z_T')
        write(*,*) columng_ptr(:)
      end if
      
      if ( col_getNumLev(columng,'MM') > 0 .and. col_varExist(columng,'Z_M') ) then
        write(*,*) 'Z_M(columnhr):'
        columng_ptr => col_getColumn(columnhr,1,'Z_M')
        write(*,*) columng_ptr(:)
        write(*,*) 'Z_M(columng):'
        columng_ptr => col_getColumn(columng,1,'Z_M')
        write(*,*) columng_ptr(:)
      end if
 
      write(*,*) 'HeightSfc:', columng%HeightSfc(1,1)
    end if

    call tmg_stop(10)

    write(*,*) 'inn_setupBackgroundColumnsAnl: END'

  end subroutine inn_setupBackgroundColumnsAnl

  !--------------------------------------------------------------------------
  ! inn_computeInnovation
  !--------------------------------------------------------------------------
  subroutine inn_computeInnovation(columnhr,obsSpaceData,destObsColumn_opt,beSilent_opt)
    !
    !:Purpose: To initialize observation innovations using the nonlinear H
    implicit none

    ! Arguments:
    type(struct_columnData) :: columnhr
    type(struct_obs)        :: obsSpaceData
    integer, optional       :: destObsColumn_opt ! column where result stored, default is OBS_OMP
    logical, optional       :: beSilent_opt

    ! Locals:
    real(8) :: Jo, JoRaob, JoSatWind, JoSurfc
    real(8) :: JoSfcSF, JoSfcUA, JoTov, JoAirep, JoSfcSC, JoProf, JoAladin, JoSfcTM
    real(8) :: JoGpsRO, JoGpsGB, JoSfcGP, JoSfcRA, JoChm, JoSfcGL, JoSfchy, JoRadvel
    integer :: destObsColumn, get_max_rss
    logical :: lgpdata, beSilent

    write(*,*)
    write(*,*) '--Starting subroutine inn_computeInnovation--'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    !
    !- Subroutine options seting
    !
    if ( present(destObsColumn_opt) ) then
      destObsColumn = destObsColumn_opt
    else
      destObsColumn = obs_omp
    end if

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .false.
    end if

    if ( .not.beSilent ) write(*,*) 'oti_timeBinning: Before filtering done in inn_computeInnovation'
    if ( .not.beSilent ) call oti_timeBinning(obsSpaceData,tim_nstepobs)

    !
    !- Pre-processing
    !

    ! Reject observed elements too far below the surface
    call filt_topo(columnhr,obsSpaceData,beSilent)
    ! ( Note: Pressure values for elements slightly below the surface are replaced by the surface
    !   pressure values of the trial field. GB-GPS (met and ZTD) observations are processed in
    !   s/r filt_topoSFC in obsFilter_mod.ftn90 )

    ! Remove surface station wind observations
    if (trim(innovationMode) == 'analysis' .or. trim(innovationMode) == 'FSO') call filt_surfaceWind(obsSpaceData,beSilent)

    ! Find interpolation layer in model profiles

    if ( col_getNumLev(columnhr,'MM') > 1 ) call oop_vobslyrs(columnhr, obsSpaceData, beSilent)

    !
    !- Calculate the innovations [Y - H(Xb)] and place the result in obsSpaceData in OBS_OMP column
    !
    call tmg_start(48,'NL_OBS_OPER')
    
    ! Radiosondes
    call oop_ppp_nl(columnhr, obsSpaceData, beSilent, JoRaob, 'UA', destObsColumn)

    ! Aircrafts
    call oop_ppp_nl(columnhr, obsSpaceData, beSilent, JoAirep, 'AI', destObsColumn)

    ! SatWinds
    call oer_sw(columnhr,obsSpaceData)
    call oop_ppp_nl(columnhr, obsSpaceData, beSilent, JoSatWind, 'SW', destObsColumn)

    ! Surface (SF, UA, SC, GP and RA families)
    call oop_sfc_nl(columnhr, obsSpaceData, beSilent, JoSfcSF, 'SF', destObsColumn)
    call oop_sfc_nl(columnhr, obsSpaceData, beSilent, JoSfcUA, 'UA', destObsColumn)
    call oop_sfc_nl(columnhr, obsSpaceData, beSilent, JoSfcSC, 'SC', destObsColumn)
    call oop_sfc_nl(columnhr, obsSpaceData, beSilent, JoSfcGP, 'GP', destObsColumn)
    call oop_sfc_nl(columnhr, obsSpaceData, beSilent, JoSfcRA, 'RA', destObsColumn)

    ! RADAR Doppler velocity
    call oop_raDvel_nl(columnhr,obsSpaceData, beSilent,JoRadVel,'RA', destObsColumn)  

    ! Sea surface temperature
    call oop_sst_nl(columnhr, obsSpaceData, beSilent, JoSfcTM, 'TM', destObsColumn)

    ! Sea ice concentration
    call filt_iceConcentration(obsSpaceData, beSilent)
    call filt_backScatAnisIce(obsSpaceData, beSilent)
    call oer_setErrBackScatAnisIce(columnhr, obsSpaceData, beSilent)
    call oop_ice_nl(columnhr, obsSpaceData, beSilent, JoSfcGL, 'GL', destObsColumn)

    ! Hydrology
    call oop_hydro_nl(columnhr, obsSpaceData, beSilent, JoSFCHY, 'HY', destObsColumn)

    ! TOVS / Radiances
    if (trim(innovationMode) == 'bgck'  ) then
      call oop_tovs_nl(columnhr, obsSpaceData, tim_getDatestamp(),  &
                       beSilent, JoTov, bgckMode_opt=.true., destObs_opt=destObsColumn)
    else
      call oop_tovs_nl(columnhr, obsSpaceData, tim_getDatestamp(),  &
                       beSilent, JoTov, bgckMode_opt=.false., destObs_opt=destObsColumn)
    end if

    ! Profilers
    call oop_zzz_nl(columnhr, obsSpaceData, beSilent, JoProf, 'PR', destObsColumn)

    ! Aladin winds
    call oop_zzz_nl(columnhr, obsSpaceData, beSilent, JoAladin, 'AL', destObsColumn)

    ! GPS radio occultation
    JoGpsRO=0.0D0
    if (obs_famExist(obsSpaceData,'RO', localMPI_opt = .true. )) then
       call filt_gpsro(columnhr, obsSpaceData, beSilent)
       call oer_SETERRGPSRO(columnhr, obsSpaceData, beSilent)
       call oop_gpsro_nl(columnhr, obsSpaceData, beSilent, JoGpsRO, destObsColumn)
    end if

    ! Chemical constituents
    call oop_chm_nl(columnhr, obsSpaceData, JoChm, destObsColumn)

    ! GPS ground-based zenith delay
    JoGpsGB=0.0D0
    if (obs_famExist(obsSpaceData,'GP', localMPI_opt = .true. )) then
      if (trim(innovationMode) == 'analysis' .or. trim(innovationMode) == 'FSO') then
        call oer_SETERRGPSGB(columnhr, obsSpaceData, beSilent, lgpdata, .true.)
        if (lgpdata) call oop_gpsgb_nl(columnhr, obsSpaceData, beSilent, JoGpsGB,  &
                                       destObsColumn, analysisMode_opt=.true.)
      else
        call oer_SETERRGPSGB(columnhr, obsSpaceData, beSilent, lgpdata, .false.)
        if (lgpdata) call oop_gpsgb_nl(columnhr, obsSpaceData, beSilent, JoGpsGB,  &
                                       destObsColumn, analysisMode_opt=.false.)
      end if
    end if

    call tmg_stop(48)

    ! Combine surface contributions
    JoSurfc = JoSfcUA + JoSfcSF + JoSfcSC + JoSfcGP + JoSfcRA + JoSfcTM + JoSfcGL + JoSFCHY

    ! Compute total Jo
    Jo = JoRaob  + JoAirep + JoSatWind + JoSurfc + JoTov + JoProf + JoAladin + &
         JoGpsRO + JoGpsGB + JoChm     + JoRadVel

    ! Print Jo components
    if ( .not.beSilent ) then
      write(*,*) 'Cost function values for this MPI task:'
      write(*,'(a15,f30.16)') 'JoRaob   = ',JoRaob
      write(*,'(a15,f30.16)') 'JoAirep  = ',JoAirep
      write(*,'(a15,f30.16)') 'JoSurfc  = ',JoSurfc
      write(*,'(a15,f30.16)') 'JoSfcTM  = ',JoSfcTM
      write(*,'(a15,f30.16)') 'JoSfcGL  = ',JoSfcGL
      write(*,'(a15,f30.16)') 'JoSfcSF  = ',JoSfcSF
      write(*,'(a15,f30.16)') 'JoSfcUA  = ',JoSfcUA
      write(*,'(a15,f30.16)') 'JoSfcSC  = ',JoSfcSC
      write(*,'(a15,f30.16)') 'JoSfcGP  = ',JoSfcGP
      write(*,'(a15,f30.16)') 'JoSfcRA  = ',JoSfcRA
      write(*,'(a15,f30.16)') 'JoTov    = ',JoTov
      write(*,'(a15,f30.16)') 'JoSatWind= ',JoSatWind
      write(*,'(a15,f30.16)') 'JoProf   = ',JoProf
      write(*,'(a15,f30.16)') 'JoAladin = ',JoAladin
      write(*,'(a15,f30.16)') 'JoGpsRO  = ',JoGpsRO
      write(*,'(a15,f30.16)') 'JoGpsGB  = ',JoGpsGB
      write(*,'(a15,f30.16)') 'JoChm    = ',JoChm
      write(*,'(a15,f30.16)') 'JoRadVel = ',JoRadVel
      write(*,'(a15,f30.16)') 'Total Jo = ',Jo

      call mpi_allreduce_sumreal8scalar(JoRaob,'GRID')
      call mpi_allreduce_sumreal8scalar(JoAirep,'GRID')
      call mpi_allreduce_sumreal8scalar(JoSurfc,'GRID')
      call mpi_allreduce_sumreal8scalar(JoSfcSF,'GRID')
      call mpi_allreduce_sumreal8scalar(JoSfcUA,'GRID')
      call mpi_allreduce_sumreal8scalar(JoSfcSC,'GRID')
      call mpi_allreduce_sumreal8scalar(JoSfcGP,'GRID')
      call mpi_allreduce_sumreal8scalar(JoSfcRA,'GRID')
      call mpi_allreduce_sumreal8scalar(JoTov,'GRID')
      call mpi_allreduce_sumreal8scalar(JoSatWind,'GRID')
      call mpi_allreduce_sumreal8scalar(JoProf,'GRID')
      call mpi_allreduce_sumreal8scalar(JoAladin,'GRID')
      call mpi_allreduce_sumreal8scalar(JoGpsRO,'GRID')
      call mpi_allreduce_sumreal8scalar(JoGpsGB,'GRID')
      call mpi_allreduce_sumreal8scalar(JoChm,'GRID')
      call mpi_allreduce_sumreal8scalar(JoSfcTM,'GRID')
      call mpi_allreduce_sumreal8scalar(JoSfcGL,'GRID')
      call mpi_allreduce_sumreal8scalar(JoRadVel,'GRID')

      write(*,*) 'Cost function values summed for all MPI tasks:'
      write(*,'(a15,f30.16)') 'JoRaob   = ',JoRaob
      write(*,'(a15,f30.16)') 'JoAirep  = ',JoAirep
      write(*,'(a15,f30.16)') 'JoSurfc  = ',JoSurfc
      write(*,'(a15,f30.16)') 'JoSfcTM  = ',JoSfcTM
      write(*,'(a15,f30.16)') 'JoSfcGL  = ',JoSfcGL
      write(*,'(a15,f30.16)') 'JoSfcSF  = ',JoSfcSF
      write(*,'(a15,f30.16)') 'JoSfcUA  = ',JoSfcUA
      write(*,'(a15,f30.16)') 'JoSfcSC  = ',JoSfcSC
      write(*,'(a15,f30.16)') 'JoSfcGP  = ',JoSfcGP
      write(*,'(a15,f30.16)') 'JoSfcRA  = ',JoSfcRA
      write(*,'(a15,f30.16)') 'JoTov    = ',JoTov
      write(*,'(a15,f30.16)') 'JoSatWind= ',JoSatWind
      write(*,'(a15,f30.16)') 'JoProf   = ',JoProf
      write(*,'(a15,f30.16)') 'JoAladin = ',JoAladin
      write(*,'(a15,f30.16)') 'JoGpsRO  = ',JoGpsRO
      write(*,'(a15,f30.16)') 'JoGpsGB  = ',JoGpsGB
      write(*,'(a15,f30.16)') 'JoChm    = ',JoChm
      write(*,'(a15,f30.16)') 'JoRadVel = ',JoRadVel
    end if ! beSilent

    ! Sum Jo over all MPI tasks
    call mpi_allreduce_sumreal8scalar(Jo,'GRID')
    write(*,'(a15,f30.16)') 'Total Jo = ',Jo

    if ( .not.beSilent ) write(*,*) 'oti_timeBinning: After filtering done in inn_computeInnovation'
    if ( .not.beSilent ) call oti_timeBinning(obsSpaceData,tim_nstepobs)

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) '--Done subroutine inn_computeInnovation--'

  end subroutine inn_computeInnovation

  !--------------------------------------------------------------------------
  ! setObsMpiStrategy
  !--------------------------------------------------------------------------
  subroutine setObsMpiStrategy(obsSpaceData, hco_anl, mpiStrategy)
    !
    !:Purpose: To distribute header indices following the chosen strategy,
    !          current options: "LIKESPLITFILES", "ROUNDROBIN", "LATLONTILES".
    implicit none

    ! Arguments:
    type(struct_obs),          intent(inout) :: obsSpaceData
    type(struct_hco), pointer, intent(in)    :: hco_anl
    character(len=*),          intent(in)    :: mpiStrategy

    ! Locals:
    real(8) :: lat_r8, lon_r8
    real    :: lat_r4, lon_r4
    real    :: xpos_r4, ypos_r4

    integer :: numHeaderFile, headerIndex, latIndex, lonIndex, ierr
    integer :: IP, IP_x, IP_y
    integer :: gdxyfll
    !
    !- Determine obs_ipc (column) and obs_ipt (tile) according to distribution strategy
    !
    write(*,*)
    numHeaderFile = obs_numheader(obsSpaceData)
    write(*,*) 'setObsMpiStrategy: numHeader after reading files = ',numHeaderFile
    write(*,*) 'setObsMpiStrategy: strategy = ',trim(mpiStrategy)

    ! make sure old PE index is not used (set to -1)
    do headerIndex = 1, numHeaderFile
       call obs_headSet_i(obsSpaceData,OBS_IP,headerIndex,-1)
    end do

    ! set PE index for the file (where reading and updating was/will be done)
    if ( obs_columnActive_IH(obsSpaceData,OBS_IPF) ) then
      do headerIndex = 1, numHeaderFile
        call obs_headSet_i(obsSpaceData,OBS_IPF,headerIndex,mpi_myid)
      end do
    end if

    select case (trim(mpiStrategy))
    case ('LIKESPLITFILES')
       !- Keep distribution exactly as it is in the split files:
       do headerIndex = 1, numHeaderFile
          if ( obs_columnActive_IH(obsSpaceData,OBS_IPC) ) &
            call obs_headSet_i(obsSpaceData,OBS_IPC,headerIndex, mpi_myid)
          if ( obs_columnActive_IH(obsSpaceData,OBS_IPT) ) &
            call obs_headSet_i(obsSpaceData,OBS_IPT,headerIndex, mpi_myid)
          if ( obs_columnActive_IH(obsSpaceData,OBS_IP) ) &
            call obs_headSet_i(obsSpaceData,OBS_IP,headerIndex, mpi_myid)
       end do
    case ('ROUNDROBIN')
       !- Distribute by a round-robin strategy for both obs_ipc and obs_ipt:
       !  (Only use if files already not split by round robin)
       do headerIndex = 1, numHeaderFile
          IP = mod((headerIndex-1),mpi_nprocs)
          if ( obs_columnActive_IH(obsSpaceData,OBS_IPC) ) &
            call obs_headSet_i(obsSpaceData,OBS_IPC,headerIndex, IP)
          if ( obs_columnActive_IH(obsSpaceData,OBS_IPT) ) &
            call obs_headSet_i(obsSpaceData,OBS_IPT,headerIndex, IP)
          if ( obs_columnActive_IH(obsSpaceData,OBS_IP) ) &
            call obs_headSet_i(obsSpaceData,OBS_IP ,headerIndex, IP)
       end do
    case ('LATLONTILES')
       !- Distribute by latitude/longitude tiles for both obs_ipc and obs_ipt:
       do headerIndex = 1, numHeaderFile
          ! compute grid index for each observation header
          lat_r8 = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)
          lon_r8 = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex)
          lat_r4 = real(lat_r8) * MPC_DEGREES_PER_RADIAN_R4
          lon_r4 = real(lon_r8) * MPC_DEGREES_PER_RADIAN_R4
          ierr = gdxyfll( hco_anl%EZscintID,   & ! IN 
               xpos_r4, ypos_r4,    & ! OUT
               lat_r4, lon_r4, 1 )    ! IN

          ! compute correponding mpi task id for each observation
          latIndex = floor(ypos_r4)
          lonIndex = floor(xpos_r4)
          IP_y = mpivar_myidYfromLat(latIndex, hco_anl%nj)
          IP_x = mpivar_myidXfromLon(lonIndex, hco_anl%ni)
          IP = IP_x + IP_y*mpi_npex

          call obs_headSet_i(obsSpaceData,OBS_IPC,headerIndex, IP)
          call obs_headSet_i(obsSpaceData,OBS_IPT,headerIndex, IP)
       end do
    case ('LATLONTILESBALANCED')
       !- Distribute by latitude/longitude tiles, but with simple & cheap balancing for obs_ipc (1 send or recv):
       call utl_abort('setObsMpiStrategy: Sorry, LATLONTILESBALANCED no longer available')
    case default
       write(*,*)
       write(*,*) 'ERROR unknown mpiStrategy: ', trim(mpiStrategy)
       call utl_abort('setObsMpiStrategy')
    end select

  end subroutine setObsMpiStrategy

  !--------------------------------------------------------------------------
  ! inn_perturbObs
  !--------------------------------------------------------------------------
  subroutine inn_perturbObs(obsSpaceData,numAnalyses,indexAnalysis, &
                            indexBatch,obs_column_index_src,obs_column_index_dest)
    !
    !:Purpose: To perturb the innovation vector to simulate effect of
    !          observation uncertainty
    !
    !.. WARNING:: perturbations are not the same when MPI topology changes!!!
    !
    implicit none

    ! Arguments:
    type(struct_obs) :: obsSpaceData
    integer :: numAnalyses,indexAnalysis,indexBatch
    integer :: obs_column_index_src,obs_column_index_dest

    ! Locals:
    integer :: numPerturbations
    integer :: nrandseed,iseed,indexAnalysis2,indexBody,indexFamily
    integer, parameter :: numFamily=10
    real*8  :: zmean,originalOmp
    real*8  :: scaleFactor(numFamily)
    character(len=2) :: familyList(numFamily)
    real*8, save, pointer :: obsPerturbations(:,:) => NULL()
    logical, save :: firstTime = .true.

    familyList(1)='UA' ; scaleFactor(1)=1.00d0
    familyList(2)='AI' ; scaleFactor(2)=1.00d0
    familyList(3)='SF' ; scaleFactor(3)=1.00d0
    familyList(4)='TO' ; scaleFactor(4)=1.00d0
    familyList(5)='SW' ; scaleFactor(5)=1.00d0
    familyList(6)='SC' ; scaleFactor(6)=1.00d0
    familyList(7)='PR' ; scaleFactor(7)=1.00d0
    familyList(8)='RO' ; scaleFactor(8)=1.00d0
    familyList(9)='GP' ; scaleFactor(9)=1.00d0
    familyList(10)='RA'; scaleFactor(10)=1.00d0    

    numPerturbations = numAnalyses

    if (firstTime) then

       if (.not.associated(obsPerturbations)) then
          write(*,*) 'perturbObs: allocating space for all perturbations'
          allocate(obsPerturbations(obs_numBody(obsSpaceData),numPerturbations))
       end if

       write(*,*) 'perturbObs: computing random numbers'

       ! compute random perturbations
       do indexAnalysis2 = 1,numPerturbations
          nrandseed   = indexAnalysis2 + (indexBatch-1)*numAnalyses
          iseed     = ABS(nrandseed)
          write(*,*) 'perturbobs: indexAnalysis, iseed=',indexAnalysis2, iseed
          call rng_setup(iseed) ! JFC: should be called only once, no???
          do indexBody = 1,obs_numBody(obsSpaceData)
             obsPerturbations(indexBody,indexAnalysis2)=rng_gaussian()*obs_bodyElem_r(obsSpaceData,OBS_OER,indexBody)
          end do
       end do

       ! apply scale factor
       do indexFamily = 1,numFamily
          do indexBody = 1,obs_numBody(obsSpaceData)      
             if (obs_getFamily(obsSpaceData,bodyIndex=indexBody).eq.familyList(indexFamily)) then
                do indexAnalysis2 = 1,numPerturbations
                   obsPerturbations(indexBody,indexAnalysis2)=obsPerturbations(indexBody,indexAnalysis2)*scaleFactor(indexFamily)
                end do
             end if
          end do
       end do

       ! remove ensemble mean
       do indexBody = 1,obs_numBody(obsSpaceData)      
          zmean=0.0d0
          do indexAnalysis2 = 1,numPerturbations
             zmean=zmean+obsPerturbations(indexBody,indexAnalysis2)
          end do
          zmean=zmean/numPerturbations
          do indexAnalysis2 = 1,numPerturbations
             obsPerturbations(indexBody,indexAnalysis2)=obsPerturbations(indexBody,indexAnalysis2)-zmean
          end do
       end do

       firstTime=.false.

    end if

    ! apply perturbation for current analysis
    write(*,*) 'perturbObs: applying perturbation for analysis number: ',indexAnalysis
    do indexBody = 1,obs_numBody(obsSpaceData)      
       if (obs_bodyElem_i(obsSpaceData,OBS_ASS,indexBody) == obs_assimilated) then
          originalOmp = obs_bodyElem_r(obsSpaceData,obs_column_index_src,indexBody)
          call obs_bodySet_r(obsSpaceData,obs_column_index_dest,indexBody,originalOmp+obsPerturbations(indexBody,indexAnalysis))
       end if
    end do

  end subroutine inn_perturbObs

end module innovation_mod
