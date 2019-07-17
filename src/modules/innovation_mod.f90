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
  use analysisGrid_mod
  use verticalCoord_mod
  use gridStateVector_mod
  use tt2phi_mod
  use utilities_mod
  use obsFilter_mod  
  use gps_mod
  use tovs_nl_mod
  use tovs_lin_mod
  use multi_ir_bgck_mod
  use chem_setup_mod
  use obsFiles_mod
  use randomNumber_mod
  use obsErrors_mod
  use bufr_mod
  use statetocolumn_mod
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
  subroutine inn_setupobs(obsSpaceData, obsColumnMode, obsMpiStrategy, &
       innovationMode_in, obsClean_opt )
    !
    !:Purpose: To initialize the observation parameters and constants
    implicit none

    ! Arguments:
    type(struct_obs)                        :: obsSpaceData
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
    ! Setup de the mode
    !
    innovationMode = innovationMode_in

    !
    ! Specify the active observation-array columns
    !
    call obs_class_initialize(obsColumnMode)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    !
    ! Allocate memory for observation arrays
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
    ! Set up the list of elements to be assimilated and flags for rejection
    !
    call filt_setup(innovationMode) ! IN

    !
    ! Read the observations from files
    !
    call tmg_start(11,'READ_OBS')
    call obsf_readFiles( obsSpaceData )
    call tmg_stop(11)

    !
    ! Read the NAMELIST NAMGGPSRO
    !
    if (obs_famExist(obsSpaceData,'RO')) call gps_setupro
    !
    ! Initialize GB-GPS processing (read NAMGPSGB in namelist file)
    !
    if (obs_famExist(obsSpaceData,'GP')) call gps_setupgb
    !
    ! Initialize TOVS processing
    !
    if (obs_famExist(obsSpaceData,'TO')) call tvs_setup
    !
    ! Read the NAMELIST NAMCHEM and set up additional constituent
    ! obs related info not found in obsSpaceData.
    !
    if (obs_famExist(obsSpaceData,'CH')) call chm_setup

    !
    ! Filter out data from CMA
    !
    call tmg_start(14,'SUPREP')
    call filt_suprep(obsSpaceData)
    call tmg_stop(14)

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
    ! set (OBS_IPC and OBS_IPT) or OBS_IP columns according to the chosen strategy
    !
    write(*,*)
    write(*,*) 'inn_setupObs - Using obsMpiStrategy = ', trim(obsMpiStrategy)
    call setObsMpiStrategy(obsSpaceData,obsMpiStrategy)

    !
    ! Check if burp files already split
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
    if (obs_famExist(obsSpaceData,'TO')) then
      call tvs_setupAlloc(obsSpaceData)
      if (trim(innovationMode) == 'bgckIR' ) call irbg_setup(obsSpaceData)
      ! Initialize non diagonal observation error matrices
      if ( trim(innovationMode) == 'analysis' .or. trim(innovationMode) == 'FSO') call oer_setInterchanCorr()
    end if

  end subroutine inn_setupobs


  subroutine inn_setupBackgroundColumns(columnhr,obsSpaceData)
    implicit none

    ! arguments
    type(struct_columnData) :: columnhr
    type(struct_obs)        :: obsSpaceData

    ! locals
    type(struct_gsv)          :: stateVector_trial
    type(struct_hco), pointer :: hco_trl => null()
    type(struct_vco), pointer :: vco_trl => null()
    integer                   :: varIndex, ierr, nulnam, fnom, fclos
    character(len=4)          :: varNamesToRead(1), lastVarNameToRead
    logical                   :: deallocInterpInfo, allocHeightSfc
    real(8), pointer          :: onecolumn(:)

    character(len=4), pointer :: anlVar(:)

    character(len=20) :: timeInterpType_nl  ! 'NEAREST' or 'LINEAR'
    NAMELIST /NAMINN/timeInterpType_nl

    write(*,*) 'inn_setupBackgroundColumns: START'

    timeInterpType_nl='NEAREST'

    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    if(ierr /= 0) call utl_abort('inn_setupBackgroundColumns: Error opening file flnml')
    read(nulnam,nml=naminn,iostat=ierr)
    if(ierr /= 0) write(*,*) 'WARNING: namelist block NAMINN not found, use the default values'
    write(*,nml=naminn)
    ierr = fclos(nulnam)

    call tmg_start(10,'INN_SETUPBACKGROUNDCOLUMNS')

    nullify(anlVar)
    call gsv_varNamesList(anlVar)
    call hco_SetupFromFile(hco_trl, './trlm_01', ' ', 'Trial', varName_opt=anlVar(1))

    call vco_SetupFromFile(vco_trl, './trlm_01')

    call col_setVco(columnhr,vco_trl)
    call col_allocate(columnhr,obs_numHeader(obsSpaceData),mpiLocal_opt=.true.)

    if (vco_trl%Vcode == 0) then
      allocHeightSfc = .false.
    else
      allocHeightSfc = .true.
    end if

    deallocInterpInfo = .true.

    call gsv_allocate( stateVector_trial, tim_nstepobs, hco_trl, vco_trl,  &
                       dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                       mpi_distribution_opt='Tiles', dataKind_opt=4,  &
                       allocHeightSfc_opt=allocHeightSfc, hInterpolateDegree_opt='LINEAR', &
                       beSilent_opt=.false.)
    call gsv_zero( stateVector_trial )
    call gsv_readTrials( stateVector_trial )
    call s2c_nl( stateVector_trial, obsSpaceData, columnhr, timeInterpType=timeInterpType_nl, &
                 moveObsAtPole_opt=.true., dealloc_opt=deallocInterpInfo )
    call gsv_deallocate(stateVector_trial)

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


  subroutine inn_setupBackgroundColumnsAnl(columnhr,columng,obsSpaceData)
    implicit none

    ! arguments
    type(struct_columnData) :: columng,columnhr
    type(struct_obs)        :: obsSpaceData

    ! locals
    integer :: jvar, jlev, columnIndex
    real(8), pointer :: columng_ptr(:), columnhr_ptr(:)

    write(*,*) 'inn_setupBackgroundColumnsAnl: START'

    call tmg_start(10,'INN_SETUPBACKGROUNDCOLUMNS')

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

    ! calculate pressure profiles on analysis levels
    if ( col_getNumCol(columng) > 0 .and. col_varExist(columng,'P_T') ) then
      call col_calcPressure(columng)
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
    endif

    ! vertical interpolation of 3D variables
    do jvar = 1, vnl_numvarmax3D
      if ( .not. col_varExist(columng,vnl_varNameList3D(jvar)) ) cycle
      call col_vintprof( columnhr, columng, vnl_varNameList3D(jvar), useColumnPressure_opt=.false. )

      ! Imposing a minimum value for HU
      if ( vnl_varNameList3D(jvar) == 'HU  ') then
        do columnIndex = 1, col_getNumCol(columng)
          columng_ptr => col_getColumn(columng,columnIndex,'HU')
          do jlev=1,col_getNumLev(columng,'TH')
            columng_ptr(jlev) = max(columng_ptr(jlev),col_rhumin)
          enddo
        end do
      end if
    end do

    if ( col_getNumCol(columng) > 0 .and. col_varExist(columng,'P_T') ) then
      if ( mpi_myid == 0 ) then
        write(*,*) 'inn_setupBackgroundColumnsAnl, after vintprof, COLUMNG(1):'
        write(*,*) 'P_T:'
        columng_ptr => col_getColumn(columng,1,'P_T')
        write(*,*) columng_ptr (:)
        write(*,*)
      end if
    endif

    if ( col_varExist(columng,'TT') .and. col_varExist(columng,'HU') .and. col_varExist(columng,'P0') ) then
      !
      !- Using T, q and PS to compute height for columng
      !
      do columnIndex = 1, col_getNumCol(columng)
        columng%HeightSfc(1,columnIndex) = columnhr%HeightSfc(1,columnIndex)
      end do

      ! remove the height offset for the diagnostic levels for backward compatibility only
      if ( col_varExist(columng,'Z_T') .and. .not.columng%addHeightSfcOffset ) then
        do columnIndex = 1, col_getNumCol(columng)
          columng_ptr => col_getColumn(columng,columnIndex,'Z_T')
          columng_ptr(col_getNumLev(columng,'TH')) = columng%HeightSfc(1,columnIndex)

          columng_ptr => col_getColumn(columng,columnIndex,'Z_M')
          columng_ptr(col_getNumLev(columng,'MM')) = columng%HeightSfc(1,columnIndex)
        end do
      end if

    else
      write(*,*) 'inn_setupBackgroundColumnsAnl:  height TLM calcs not generated since TT, HU and P0 not all present'
    end if

    if ( col_getNumCol(columng) > 0 .and. col_varExist(columng,'Z_T') ) then
      write(*,*) 'inn_setupBackgroundColumnsAnl, vIntProf output:'
      write(*,*) 'Z_T (columnhr):'
      columng_ptr => col_getColumn(columnhr,1,'Z_T')
      write(*,*) columng_ptr(:)
      write(*,*) 'Z_T (columng):'
      columng_ptr => col_getColumn(columng,1,'Z_T')
      write(*,*) columng_ptr(:)

      write(*,*) 'Z_M(columnhr):'
      columng_ptr => col_getColumn(columnhr,1,'Z_M')
      write(*,*) columng_ptr(:)
      write(*,*) 'Z_M(columng):'
      columng_ptr => col_getColumn(columng,1,'Z_M')
      write(*,*) columng_ptr(:)

      write(*,*) 'HeightSfc:', columng%HeightSfc(1,1)
    end if

    call tmg_stop(10)

    write(*,*) 'inn_setupBackgroundColumnsAnl: END'

  end subroutine inn_setupBackgroundColumnsAnl


  subroutine inn_computeInnovation(columnhr,obsSpaceData,beSilent_opt)
    !
    !:Purpose: To initialize Observation Innovations using the nonlinear H
    implicit none

    ! Arguments:
    type(struct_columnData) :: columnhr
    type(struct_obs)        :: obsSpaceData
    logical, optional       :: beSilent_opt

    ! Locals:
    real(8) :: zjo,zjoraob,zjosatwind,zjosurfc
    real(8) :: zjosfcsf,zjosfcua,zjotov,zjoairep,zjosfcsc,zjoprof,zjoaladin,zjosfctm
    real(8) :: zjogpsro,zjogpsgb,zjosfcgp,zjochm,zjosfcgl,zjosfchy
    integer :: ierr, get_max_rss
    logical :: lgpdata, beSilent

    write(*,*) '--Starting subroutine inn_computeInnovation--'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if ( present(beSilent_opt) ) then
      beSilent = beSilent_opt
    else
      beSilent = .false.
    end if

    if ( .not.beSilent ) write(*,*) 'oti_timeBinning: Before filtering done in inn_computeInnovation'
    if ( .not.beSilent ) call oti_timeBinning(obsSpaceData,tim_nstepobs)
    !
    !     Reject observed elements too far below the surface. Pressure values
    !     for elements slightly below the surface are replaced by the surface
    !     pressure values of the trial field.
    !
    !     GB-GPS (met and ZTD) observations are processed in s/r filt_topoSFC (in obsFilter_mod.ftn90)
    !
    call filt_topo(columnhr,obsSpaceData,beSilent)
    !
    !     Remove surface station wind observations
    !
    if (trim(innovationMode) == 'analysis' .or. trim(innovationMode) == 'FSO') call filt_surfaceWind(obsSpaceData,beSilent)
    !
    !     Find interpolation layer in model profiles
    !
    if ( col_getNumLev(columnhr,'MM') > 1 ) call oop_vobslyrs(columnhr, obsSpaceData, beSilent)
    !
    !
    !------ Calculate the innovations Y - H(Xb) and place
    !       the result in obsSpaceData in OBS_OMP column
    !
    !        RAOBS
    !------------------------------
    !
    call tmg_start(48,'NL_OBS_OPER')
    call oop_ppp_nl(columnhr, obsSpaceData, beSilent, ZJORAOB, 'UA')
    !
    !        AIREPS
    !--------------------------------
    call oop_ppp_nl(columnhr, obsSpaceData, beSilent, ZJOAIREP, 'AI')
    !
    !        SATWINDS
    !--------------------------------
    call oer_sw(columnhr,obsSpaceData)
    call oop_ppp_nl(columnhr, obsSpaceData, beSilent, ZJOSATWIND, 'SW')
    !
    !        SURFACE (SF, UA, SC AND GP FAMILIES)
    !-------------------------------
    call oop_sfc_nl  (columnhr, obsSpaceData, beSilent, ZJOSFCSF, 'SF')
    call oop_sfc_nl  (columnhr, obsSpaceData, beSilent, ZJOSFCUA, 'UA')
    call oop_sfc_nl  (columnhr, obsSpaceData, beSilent, ZJOSFCSC, 'SC')
    call oop_sfc_nl  (columnhr, obsSpaceData, beSilent, ZJOSFCGP, 'GP')
    call oop_sst_nl  (columnhr, obsSpaceData, beSilent, ZJOSFCTM, 'TM')
    call oop_ice_nl  (columnhr, obsSpaceData, beSilent, ZJOSFCGL, 'GL')
    call oop_hydro_nl(columnhr, obsSpaceData, beSilent, ZJOSFCHY, 'HY')

    ZJOSURFC = ZJOSFCUA + ZJOSFCSF + ZJOSFCSC + ZJOSFCGP + ZJOSFCTM + ZJOSFCGL + ZJOSFCHY
    !
    !        TOVS - RADIANCE
    !-------------------------------
    if (trim(innovationMode) == 'bgckIR'  ) then
      call oop_tovs_nl(columnhr,obsSpaceData,tim_getDatestamp(),filt_rlimlvhu,beSilent,ZJOTOV,bgckMode_opt=.true.)
    else
      call oop_tovs_nl(columnhr,obsSpaceData,tim_getDatestamp(),filt_rlimlvhu,beSilent,ZJOTOV,bgckMode_opt=.false.)
    end if
    !
    !        PROFILER
    !------------------------------
    call oop_zzz_nl(columnhr, obsSpaceData, beSilent, ZJOPROF, 'PR')
    !
    !        GEOMETRIC HEIGHT - ALADIN WINDS
    !------------------------------
    call oop_zzz_nl(columnhr, obsSpaceData, beSilent, ZJOALADIN, 'AL')
    !
    !        GPS - RADIO OCCULTATION
    !-------------------------------
    ZJOGPSRO=0.0D0
    if (obs_famExist(obsSpaceData,'RO', localMPI_opt = .true. )) then
       CALL filt_gpsro(columnhr, obsSpaceData, beSilent)
       CALL oer_SETERRGPSRO(columnhr, obsSpaceData, beSilent)
       call oop_gpsro_nl(columnhr, obsSpaceData, beSilent, ZJOGPSRO)
    end if
    !
    !        CH - CHEMICAL CONSTITUENTS
    !-------------------------------
    call oop_chm_nl(columnhr, obsSpaceData, beSilent, zjochm)
    !
    !        GPS - GROUND-BASED ZENITH DELAY
    !-------------------------------
    !
    ZJOGPSGB=0.0D0
    if (obs_famExist(obsSpaceData,'GP', localMPI_opt = .true. )) then
      if (trim(innovationMode) == 'analysis' .or. trim(innovationMode) == 'FSO') then
        call oer_SETERRGPSGB(columnhr, obsSpaceData, beSilent, lgpdata, .true.)
        if (lgpdata) call oop_gpsgb_nl(columnhr, obsSpaceData, beSilent, ZJOGPSGB, .true.)
      else
        call oer_SETERRGPSGB(columnhr, obsSpaceData, beSilent, lgpdata, .false.)
        if (lgpdata) call oop_gpsgb_nl(columnhr, obsSpaceData, beSilent, ZJOGPSGB, .false.)
      end if
    end if

    call tmg_stop(48)

    if ( .not.beSilent ) then
      write(*,*) 'Cost function values for this MPI task:'
      write(*,'(a15,f30.16)') 'JORAOB   = ',ZJORAOB
      write(*,'(a15,f30.16)') 'JOAIREP  = ',ZJOAIREP
      write(*,'(a15,f30.16)') 'JOSURFC  = ',ZJOSURFC
      write(*,'(a15,f30.16)') 'JOSURTM  = ',ZJOSFCTM
      write(*,'(a15,f30.16)') 'JOSURGL  = ',ZJOSFCGL
      write(*,'(a15,f30.16)') 'JOSFCSF  = ',ZJOSFCSF
      write(*,'(a15,f30.16)') 'JOSFCUA  = ',ZJOSFCUA
      write(*,'(a15,f30.16)') 'JOSFCSC  = ',ZJOSFCSC
      write(*,'(a15,f30.16)') 'JOSFCGP  = ',ZJOSFCGP
      write(*,'(a15,f30.16)') 'JOTOV    = ',ZJOTOV
      write(*,'(a15,f30.16)') 'JOSATWIND= ',ZJOSATWIND
      write(*,'(a15,f30.16)') 'JOPROF   = ',ZJOPROF
      write(*,'(a15,f30.16)') 'JOALADIN = ',ZJOALADIN
      write(*,'(a15,f30.16)') 'JOGPSRO  = ',ZJOGPSRO
      write(*,'(a15,f30.16)') 'JOGPSGB  = ',ZJOGPSGB
      write(*,'(a15,f30.16)') 'JOCHM    = ',ZJOCHM
    end if

    !
    !=======================================================================
    ZJO =  ZJORAOB + ZJOAIREP + ZJOSATWIND + &
         ZJOSURFC + ZJOTOV + ZJOPROF + ZJOALADIN + ZJOGPSRO + ZJOGPSGB + ZJOCHM
    !=======================================================================

    if ( .not.beSilent ) then
      write(*,'(a15,f30.16)') 'Total Jo = ',ZJO

      call mpi_allreduce_sumreal8scalar(ZJORAOB,'GRID')
      call mpi_allreduce_sumreal8scalar(ZJOAIREP,'GRID')
      call mpi_allreduce_sumreal8scalar(ZJOSURFC,'GRID')
      call mpi_allreduce_sumreal8scalar(ZJOSFCSF,'GRID')
      call mpi_allreduce_sumreal8scalar(ZJOSFCUA,'GRID')
      call mpi_allreduce_sumreal8scalar(ZJOSFCSC,'GRID')
      call mpi_allreduce_sumreal8scalar(ZJOSFCGP,'GRID')
      call mpi_allreduce_sumreal8scalar(ZJOTOV,'GRID')
      call mpi_allreduce_sumreal8scalar(ZJOSATWIND,'GRID')
      call mpi_allreduce_sumreal8scalar(ZJOPROF,'GRID')
      call mpi_allreduce_sumreal8scalar(ZJOALADIN,'GRID')
      call mpi_allreduce_sumreal8scalar(ZJOGPSRO,'GRID')
      call mpi_allreduce_sumreal8scalar(ZJOGPSGB,'GRID')
      call mpi_allreduce_sumreal8scalar(ZJOCHM,'GRID')
      call mpi_allreduce_sumreal8scalar(ZJOSFCTM,'GRID')
      call mpi_allreduce_sumreal8scalar(ZJOSFCGL,'GRID')

      write(*,*) 'Cost function values summed for all MPI tasks:'
      write(*,'(a15,f30.16)') 'JORAOB   = ',ZJORAOB
      write(*,'(a15,f30.16)') 'JOAIREP  = ',ZJOAIREP
      write(*,'(a15,f30.16)') 'JOSURFC  = ',ZJOSURFC
      write(*,'(a15,f30.16)') 'JOSURTM  = ',ZJOSFCTM
      write(*,'(a15,f30.16)') 'JOSURGL  = ',ZJOSFCGL
      write(*,'(a15,f30.16)') 'JOSFCSF  = ',ZJOSFCSF
      write(*,'(a15,f30.16)') 'JOSFCUA  = ',ZJOSFCUA
      write(*,'(a15,f30.16)') 'JOSFCSC  = ',ZJOSFCSC
      write(*,'(a15,f30.16)') 'JOSFCGP  = ',ZJOSFCGP
      write(*,'(a15,f30.16)') 'JOTOV    = ',ZJOTOV
      write(*,'(a15,f30.16)') 'JOSATWIND= ',ZJOSATWIND
      write(*,'(a15,f30.16)') 'JOPROF   = ',ZJOPROF
      write(*,'(a15,f30.16)') 'JOALADIN = ',ZJOALADIN
      write(*,'(a15,f30.16)') 'JOGPSRO  = ',ZJOGPSRO
      write(*,'(a15,f30.16)') 'JOGPSGB  = ',ZJOGPSGB
      write(*,'(a15,f30.16)') 'JOCHM    = ',ZJOCHM

    end if ! beSilent

    call mpi_allreduce_sumreal8scalar(ZJO,'GRID')
    write(*,'(a15,f30.16)') 'Total Jo = ',ZJO

    if ( .not.beSilent ) write(*,*) 'oti_timeBinning: After filtering done in inn_computeInnovation'
    if ( .not.beSilent ) call oti_timeBinning(obsSpaceData,tim_nstepobs)

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) '--Done subroutine inn_computeInnovation--'

  end subroutine inn_computeInnovation


  subroutine setObsMpiStrategy(obsSpaceData, mpiStrategy)
    !
    !:Purpose: To distribute header indices following the chosen strategy,
    !          current options: "LIKESPLITFILES", "ROUNDROBIN", "LATLONTILES".
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData
    character(len=*), intent(in)    :: mpiStrategy

    ! Locals:
    type(struct_hco), pointer :: hco_anl

    real(8) :: lat_r8, lon_r8
    real    :: lat_r4, lon_r4
    real    :: xpos_r4, ypos_r4

    integer :: numHeaderFile, headerIndex, latIndex, lonIndex, ierr
    integer :: IP, IP_x, IP_y
    integer :: gdxyfll, get_max_rss
    !
    !- 1.  Get some info
    !

    !- 1.1 Get the horizontal coordinate of the analysis grid
    hco_anl => agd_getHco('ComputationalGrid')

    !
    !- 2.  Determine obs_ipc (column) and obs_ipt (tile) according to distribution strategy
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
       !- 2.1 Keep distribution exactly as it is in the split files:
       do headerIndex = 1, numHeaderFile
          if ( obs_columnActive_IH(obsSpaceData,OBS_IPC) ) &
            call obs_headSet_i(obsSpaceData,OBS_IPC,headerIndex, mpi_myid)
          if ( obs_columnActive_IH(obsSpaceData,OBS_IPT) ) &
            call obs_headSet_i(obsSpaceData,OBS_IPT,headerIndex, mpi_myid)
          if ( obs_columnActive_IH(obsSpaceData,OBS_IP) ) &
            call obs_headSet_i(obsSpaceData,OBS_IP,headerIndex, mpi_myid)
       end do
    case ('ROUNDROBIN')
       !- 2.2 Distribute by a round-robin strategy for both obs_ipc and obs_ipt:
       !      (Only use if files already not split by round robin)
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
       !- 2.3 Distribute by latitude/longitude tiles for both obs_ipc and obs_ipt:
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
       !- 2.4 Distribute by latitude/longitude tiles, but with simple & cheap balancing for obs_ipc (1 send or recv):

       call utl_abort('setObsMpiStrategy: Sorry, LATLONTILESBALANCED no longer available')

    case default
       write(*,*)
       write(*,*) 'ERROR unknown mpiStrategy: ', trim(mpiStrategy)
       call utl_abort('setObsMpiStrategy')
    end select

  end subroutine setObsMpiStrategy


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
    integer, parameter :: numFamily=9
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

    numPerturbations = numAnalyses

    if(firstTime) then

       if(.not.associated(obsPerturbations)) then
          write(*,*) 'perturbObs: allocating space for all perturbations'
          allocate(obsPerturbations(obs_numBody(obsSpaceData),numPerturbations))
       endif

       write(*,*) 'perturbObs: computing random numbers'

       ! compute random perturbations
       do indexAnalysis2 = 1,numPerturbations
          nrandseed   = indexAnalysis2 + (indexBatch-1)*numAnalyses
          iseed     = ABS(nrandseed)
          write(*,*) 'perturbobs: indexAnalysis, iseed=',indexAnalysis2, iseed
          call rng_setup(iseed) ! JFC: should be called only once, no???
          do indexBody = 1,obs_numBody(obsSpaceData)
             obsPerturbations(indexBody,indexAnalysis2)=rng_gaussian()*obs_bodyElem_r(obsSpaceData,OBS_OER,indexBody)
          enddo
       enddo

       ! apply scale factor
       do indexFamily = 1,numFamily
          do indexBody = 1,obs_numBody(obsSpaceData)      
             if(obs_getFamily(obsSpaceData,bodyIndex=indexBody).eq.familyList(indexFamily)) then
                do indexAnalysis2 = 1,numPerturbations
                   obsPerturbations(indexBody,indexAnalysis2)=obsPerturbations(indexBody,indexAnalysis2)*scaleFactor(indexFamily)
                enddo
             endif
          enddo
       enddo

       ! remove ensemble mean
       do indexBody = 1,obs_numBody(obsSpaceData)      
          zmean=0.0d0
          do indexAnalysis2 = 1,numPerturbations
             zmean=zmean+obsPerturbations(indexBody,indexAnalysis2)
          enddo
          zmean=zmean/numPerturbations
          do indexAnalysis2 = 1,numPerturbations
             obsPerturbations(indexBody,indexAnalysis2)=obsPerturbations(indexBody,indexAnalysis2)-zmean
          enddo
       enddo

       firstTime=.false.

    endif

    ! apply perturbation for current analysis
    write(*,*) 'perturbObs: applying perturbation for analysis number: ',indexAnalysis
    do indexBody = 1,obs_numBody(obsSpaceData)      
       if(obs_bodyElem_i(obsSpaceData,OBS_ASS,indexBody) == obs_assimilated) then
          originalOmp = obs_bodyElem_r(obsSpaceData,obs_column_index_src,indexBody)
          call obs_bodySet_r(obsSpaceData,obs_column_index_dest,indexBody,originalOmp+obsPerturbations(indexBody,indexAnalysis))
       endif
    enddo

  end subroutine inn_perturbObs


end module innovation_mod
