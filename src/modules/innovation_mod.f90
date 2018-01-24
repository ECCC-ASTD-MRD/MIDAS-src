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

!--------------------------------------------------------------------------
!! MODULE innovation_mod (prefix="inn")
!!
!! *Purpose*: Several high-level subroutines used to compute the innovations,
!!            that is, the observation-minus-background values. This includes
!!            the subroutine that reads in the gridded high-res background state
!!            from standard files.
!!
!--------------------------------------------------------------------------
module innovation_mod
  use mpi_mod
  use ramDisk_mod
  use obsSpaceData_mod
  use columnData_mod
  use timeCoord_mod
  use obsTimeInterp_mod
  use obsOperators_mod
  use EarthConstants_mod
  use MathPhysConstants_mod
  use mpivar_mod
  use horizontalCoord_mod
  use columnData_mod
  use varNameList_mod
  use analysisGrid_mod
  use gridStateVector_mod
  use WindRotation_mod
  use tt2phi_mod
  use utilities_mod
  use obsFilter_mod  
  use gps_mod
  use tovs_nl_mod
  use tovs_lin_mod
  use multi_ir_bgck_mod
  use chem_setup_mod, only: chm_setup, chm_apply_2dfieldr4_transform
  use obsFiles_mod
  use randomNumber_mod
  use obsErrors_mod
  use bufr_mod
  implicit none
  save
  private

  ! public procedures
  public :: inn_setupObs, inn_setupBackgroundColumns, inn_computeInnovation
  public :: inn_perturbObs, inn_setupBackgroundColumnsAnl

  character(len=48) :: innovationMode

contains

  !--------------------------------------------------------------------------
  ! inn_setupObs
  !--------------------------------------------------------------------------
  subroutine inn_setupobs(obsSpaceData, obsColumnMode, obsMpiStrategy, &
       innovationMode_in, obsClean_opt )
    !
    !**s/r INN_SETUPOBS  - Initialisation of observation parameters and constants
    !
    ! Revisions:
    !           Y. Rochon and M. Sitwell, Jan 2016
    !           - Use of obs_famExist and corresponding move of calls
    !             to gps_setupro, gps_setupgb and tovs_setup after the call to obsf_readFiles
    !             as initialization of family list in obs_famExist must follow
    !             saving of the obs in obsSpaceData.
    !           - Addition of call to chm_setup for inclusion of constituent data
    !             info not included in obsSpaceData_mod and reading of NAMCHEM.
    !
    IMPLICIT NONE

    type(struct_obs)                        :: obsSpaceData
    character(len=*)                        :: obsMpiStrategy
    character(len=*)                        :: obsColumnMode
    character(len=*), intent(in)            :: innovationMode_in
    logical,                       optional :: obsClean_opt

    character(len=20) :: nameDimFile
    integer :: get_max_rss, ierr, fnom, fclos, unitDimFile, mxstn, mxobs
    logical :: obsDimFileExists

    WRITE(*,FMT=9000)
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
    call filt_sethind(obsSpaceData)
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


  subroutine inn_setupBackgroundColumnsAnl(columnhr,columng)
    implicit none

    ! arguments
    type(struct_columnData) :: columng,columnhr

    ! locals
    integer :: jvar, jlev, columnIndex
    real(8), pointer :: columng_ptr(:), columnhr_ptr(:)

    call tmg_start(10,'INN_SETUPBACKGROUNDCOLUMNS')

    ! copy Lat/Lon and positions needed for interpolation
    call col_copyLatLon(columnhr,columng)

    ! copy 2D surface variables
    do jvar = 1, vnl_numvarmax2D
       if ( .not. col_varExist(vnl_varNameList2D(jvar)) ) cycle
       if ( col_getNumCol(columng) > 0 ) then       
          do columnIndex = 1, col_getNumCol(columng)
             columng_ptr  => col_getColumn( columng , columnIndex, vnl_varNameList2D(jvar) )
             columnhr_ptr => col_getColumn( columnhr, columnIndex, vnl_varNameList2D(jvar) )
             columng_ptr(:) = columnhr_ptr(:)
          enddo
       endif
    enddo

    ! calculate pressure profiles on analysis levels
    if (col_getNumCol(columng) > 0 .and. col_varExist('P0')) then
       call col_calcPressure(columng)
       do jlev = 1,col_getNumLev(columng,'MM')
          if (mpi_myid.eq.0) write(*,*) 'inn_setupBackgroundColumnsAnl: jlev, col_getPressure(COLUMNG,jlev,1,MM) = ',  &
               jlev,col_getPressure(columng,jlev,1,'MM')
       end do
       do jlev = 1,col_getNumLev(columng,'TH')
          if (mpi_myid.eq.0) write(*,*) 'inn_setupBackgroundColumnsAnl: jlev, col_getPressure(COLUMNG,jlev,1,TH) = ',  &
               jlev,col_getPressure(columng,jlev,1,'TH')
       end do
       do jlev = 1,col_getNumLev(columng,'MM')
         if (mpi_myid.eq.0) write(*,*) 'inn_setupBackgroundColumnsAnl: jlev, col_getPressureDeriv(COLUMNG,jlev,1,MM) = ',  &
              jlev,col_getPressureDeriv(columng,jlev,1,'MM')
       end do
       do jlev = 1,col_getNumLev(columng,'TH')
          if (mpi_myid.eq.0) write(*,*) 'inn_setupBackgroundColumnsAnl: jlev, col_getPressureDeriv(COLUMNG,jlev,1,TH) = ',  &
              jlev,col_getPressureDeriv(columng,jlev,1,'TH')
       end do
    endif

    ! vertical interpolation of 3D variables
    do jvar = 1, vnl_numvarmax3D
       if ( .not. col_varExist( vnl_varNameList3D(jvar) ) ) cycle
       if ( vnl_varNameList3D(jvar) == 'GZ  ') cycle

       ! do extra manipulation of HU to obtain results closer to previous 
       ! version of code that did vertical interpolation before converting HU->LQ
       !if ( vnl_varNameList3D(jvar).eq.'HU  ' ) then
       !   ! conversion from log(humidity) to specific humidity 
       !   do columnIndex = 1, col_getNumCol(columng)
       !      columnhr_ptr => col_getColumn(columnhr,columnIndex,'HU')
       !      do jlev=1,col_getNumLev(columnhr,'TH')
       !         columnhr_ptr(jlev)=exp(columnhr_ptr(jlev))
       !      enddo
       !   enddo
       !endif

       call col_vintprof( columnhr, columng, vnl_varNameList3D(jvar) )

       !if ( vnl_varNameList3D(jvar).eq.'HU  ' ) then
       !   ! conversion from specific humidity to log(humidity)
       !   do columnIndex = 1, col_getNumCol(columng)
       !      columng_ptr  => col_getColumn(columng ,columnIndex,'HU')
       !      columnhr_ptr => col_getColumn(columnhr,columnIndex,'HU')
       !      do jlev=1,col_getNumLev(columnhr,'TH')
       !         columnhr_ptr(jlev) = log(max(columnhr_ptr(jlev),col_rhumin))
       !      enddo
       !      do jlev=1,col_getNumLev(columng,'TH')
       !         columng_ptr(jlev) = log(max(columng_ptr(jlev),col_rhumin))
       !      enddo
       !   enddo
       !endif
    enddo

    if (col_varExist('TT') .and. col_varExist('HU') .and. col_varExist('P0')) then
       !
       !- Using T, q and PS to compute GZ for columng
       !
       do columnIndex = 1, col_getNumCol(columng)
          call col_setGZsfc(columng ,columnIndex, col_getGZsfc(columnhr, columnIndex))
       enddo
       if (col_getNumLev(columng,'MM') > 1) call tt2phi(columng)
    else
       write(*,*) 'inn_setupBackgroundColumnsAnl:  GZ TLM calcs not generated since TT, HU and P0 not all present'
    end if

    call tmg_stop(10)

  end subroutine inn_setupBackgroundColumnsAnl


  subroutine inn_setupBackgroundColumns(columnhr,obsSpaceData)
    !
    ! Purpose: Fill in COLUMNHR with trial profiles
    !
    ! Arguments: COLUMNHR, OBSSPACEDATA
    !
    ! Revisions:
    !   Ping Du, Oct. 2014
    !   - Added input of constituents trial fields and related setting of trial profiles.
    !     This was done by including the "case(default)" in the loop over 3D variables.
    !   Y. Rochon, Feb 2016
    !   - Additions introduced to allow variable transformations for constituent fields.
    !   - Addition of varnamelist_mod and valvarkindfromvarname
    !
    implicit none

    ! arguments
    type(struct_columnData) :: columnhr
    type(struct_obs) :: obsSpaceData

    ! locals
    type(struct_vco), pointer :: vco_trl
    type(struct_hco), pointer :: hco_anl
    integer, allocatable :: nultrl(:)
    integer              :: jlev, jobs, jvar, ierr, status, iset, jstep, jlatlontile
    integer              :: ni_trl, nj_trl, nk_trl, nultrl_forEZ
    integer              :: ig1obs,ig2obs,ig3obs,ig4obs
    integer              :: idata, idatend, jdata
    integer              :: ip1_pak_trl, ip1_vco_trl, nultrl2
    integer              :: nlevtrl_T, nlevtrl_M, numColumns, numColumn_maxmpiglobal
    integer, parameter   :: maxLevels = 200
    integer              :: EZscintID_trl, iip1s(maxLevels), iip2, iip3
    integer, allocatable :: idate(:), itime(:)
    integer, allocatable :: notag(:,:) ! (nobtot,nstepobs) obs tag associated to observations of each bin
    integer, allocatable :: nobs(:), nobs_maxmpiglobal(:) ! number of headers for each stepobs bin
    integer, allocatable :: nobsgid_mpiglobal(:,:), nobs_mpiglobal(:,:)
    integer, allocatable :: datestamplist(:)
    real(8), allocatable :: varInterphr_T(:,:), varInterphr_M(:,:), varInterphr_VV(:,:)
    real(4)              :: lat_r4, lon_r4, lat_deg_r4, lon_deg_r4, xpos_r4, ypos_r4
    real(4)              :: xposLowerBoundAnl_r4, xposUpperBoundAnl_r4
    real(8)              :: lat_r8, lon_r8, ypos_r8, xpos_r8, lat_rot, lon_rot
    real(8)              :: zig1,zig2,zig3,zig4,stepObsIndex
    real(8), allocatable :: dlonfld(:), dlatfld(:)
    real(8), allocatable :: dlonfld_mpiglobal(:,:), dlatfld_mpiglobal(:,:)
    real(8), pointer     :: column_ptr(:) => null()
    logical              :: trialExists, noUpperGZ
    character(len=1)     :: clgrtyp
    character(len=2)     :: cltypvar
    character(len=12)    :: cletiket
    character(len=2)     :: flnum
    character(len=128)   :: trialfile(100), trialfile_forEZ
    character(len=4)     :: varnameForEZ
    integer :: key, dateo, deet, npas, nbits, datyp
    integer :: ip1, ip2, ip3, swa, lng, dltf, ubc
    integer :: extra1, extra2, extra3
    integer :: ig1, ig2, ig3, ig4

    ! external functions
    integer :: get_max_rss, newdate
    integer :: ezgprm, ezqkdef, gdxyfll, gdllfxy
    integer :: fnom, fclos, fstouv, fstfrm, fstinf, fstprm

    write(*,*) ' '
    write(*,*) '-------- ENTERING INN_SETUPBACKGROUNDCOLUMNS ---------'
    write(*,*) ' '
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    call tmg_start(10,'INN_SETUPBACKGROUNDCOLUMNS')

    !
    !     Ensure that all trial field files exist and
    !     open all trial field files (assume 1 file per time step)
    !
    allocate(nultrl(tim_nStepObs))
    nultrl(:)=0
    trialfile(:) = 'NOT_DEFINED'
    do jstep = (1+mpi_myid), tim_nStepObs, mpi_nprocs
       write(flnum,'(I2.2)') jstep
       trialfile(jstep)='trlm_'//trim(flnum)
       trialfile(jstep) = ram_fullWorkingPath(trialfile(jstep))
       inquire(file=trim(trialfile(jstep)),exist=trialExists)
       if (.not.trialExists) then
          write(*,*) 'File missing=',trialfile(jstep)
          call utl_abort('INN_SETUPBACKGROUNDCOLUMNS:DID NOT FIND A TRIAL FIELD FILE')
       else
          ierr=fnom(nultrl(jstep),trim(trialfile(jstep)),'RND+OLD+R/O',0)
          ierr=fstouv(nultrl(jstep),'RND+OLD')
          write(*,*) 'ITRIAL - File :', trialfile(jstep)
          write(*,*) ' opened as unit file ',nultrl(jstep)
       end if
    enddo

    !
    !     Vertical coordinate parameters 
    !
    vco_trl => col_getVco(columnhr)
    nlevtrl_M = vco_getNumLev(vco_trl,'MM')
    nlevtrl_T = vco_getNumLev(vco_trl,'TH')
    if (mpi_myid.eq.0) write(*,*)'INN_SETUPBACKGROUNDCOLUMNS:niv thermo:',nlevtrl_T,' momentum',nlevtrl_M

    !
    !     Compute the maximum number of columns over all processors (lat-lon tiles)
    !
    numColumns = col_getNumCol(columnhr)
    call rpn_comm_allreduce(numColumns,numColumn_maxmpiglobal,1,  &
         'MPI_INTEGER','MPI_MAX','GRID',ierr)

    !
    !     Allocate trial field column object and other local arrays
    !
    if (numColumns.gt.0) then
       allocate(notag(numColumns,tim_nStepObs))
       allocate(varInterphr_T(nlevtrl_T,numColumns))
       allocate(varInterphr_M(nlevtrl_M,numColumns))
       allocate(varInterphr_VV(nlevtrl_M,numColumns))
       varInterphr_T(:,:)=0.0d0
       varInterphr_M(:,:)=0.0d0
       varInterphr_VV(:,:)=0.0d0
    endif

    allocate(dlonfld(numColumn_maxmpiglobal))
    allocate(dlatfld(numColumn_maxmpiglobal))
    allocate(dlonfld_mpiglobal(numColumn_maxmpiglobal,mpi_nprocs))
    allocate(dlatfld_mpiglobal(numColumn_maxmpiglobal,mpi_nprocs))
    allocate(nobs(tim_nStepObs))
    allocate(nobs_maxmpiglobal(tim_nStepObs))
    allocate(datestamplist(tim_nStepObs))
    allocate(idate(tim_nStepObs))
    allocate(itime(tim_nStepObs))
    allocate(nobsgid_mpiglobal(tim_nStepObs,mpi_nprocs))
    allocate(nobs_mpiglobal(tim_nStepObs,mpi_nprocs))

    !
    !     Computing date and time of step obs for error message
    !
    call tim_getstamplist(datestamplist,tim_nStepObs,tim_getDatestamp())
    do jstep = 1,tim_nStepObs
       ierr = newdate(datestamplist(jstep),idate(jstep),itime(jstep),-3)
       if (mpi_myid.eq.0) write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS: datestamplist=',jstep,datestamplist(jstep)
    end do

    !
    !-    Get the Analysis Grid structure
    !
    hco_anl => agd_getHco('CoreGrid')

    if ( hco_anl % global ) then
       xposLowerBoundAnl_r4 = - huge(1.0) ! no limit since grid is global (periodic)
       xposUpperBoundAnl_r4 = + huge(1.0) ! no limit since grid is global (periodic)
    else
       xposLowerBoundAnl_r4 = 1.0
       xposUpperBoundAnl_r4 = real(hco_anl % ni)
    end if

    !
    !- Get horizontal grid parameters to be used to test grid bounds
    !
    varnameForEZ='NONE'
    do jvar=1,vnl_numvarmax2D
       if (gsv_varExist(varName=vnl_varNameList2D(jvar))) then
          varnameForEZ=vnl_varNameList2D(jvar)
          exit
       end if
    end do
    if (trim(varnameForEZ) == 'NONE' ) then
       do jvar=1,vnl_numvarmax3D
          if (gsv_varExist(varName=vnl_varNameList3D(jvar))) then
             varnameForEZ=vnl_varNameList3D(jvar)
             exit
          end if
       end do
    end if
    if (mpi_myid.eq.0) write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS: varname for determining grid =',trim(varnameForEZ)

    if ( (1+mpi_myid) <= tim_nStepObs ) then
      ! already open (possibly on ram disk)
      nultrl_forEZ = nultrl(1+mpi_myid)
    else
      ! get from working directory (disk)
      trialfile_forEZ = 'trlm_01'
      trialfile_forEZ = ram_fullWorkingPath(trialfile_forEZ)
      nultrl_forEZ = 0
      ierr = fnom(nultrl_forEZ, trim(trialfile_forEZ),'RND+OLD+R/O',0)
      ierr = fstouv(nultrl_forEZ,'RND+OLD')
    endif

    dateo  = -1
    cletiket = ' '
    ip1    = -1
    ip2    = -1
    ip3    = -1
    cltypvar = ' '
    key = fstinf( nultrl_forEZ,                                          & ! IN
                  ni_trl, nj_trl, nk_trl,                                & ! OUT
                  dateo, cletiket, ip1, ip2, ip3, cltypvar, varnameForEZ ) ! IN

    if (key < 0) then
       write(*,*)
       write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS: nultrl_forEZ = ',nultrl_forEZ
       write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS: Unable to find trial field = ',varnameForEZ
       call utl_abort('INN_SETUPBACKGROUNDCOLUMNS')
    end if

    ierr = fstprm( key,                                                     & ! IN
                   dateo, deet, npas, ni_trl, nj_trl, nk_trl, nbits,        & ! OUT
                   datyp, ip1, ip2, ip3, cltypvar, varnameForEZ, cletiket,  & ! OUT
                   clgrtyp, ig1, ig2, ig3,                                  & ! OUT
                   ig4, swa, lng, dltf, ubc, extra1, extra2, extra3 )         ! OUT

    EZscintID_trl = ezqkdef(ni_trl,nj_trl,clgrtyp,ig1,ig2,ig3,ig4,nultrl_forEZ)

    nobs(:) = 0

    do jstep = 1,tim_nStepObs

       dlonfld(:)=0.0d0
       dlatfld(:)=0.0d0

       do jobs=1, obs_numheader(obsSpaceData)

          call tim_getStepObsIndex(stepObsIndex,tim_getDatestamp(),  &
               obs_headElem_i(obsSpaceData,OBS_DAT,jobs),  &
               obs_headElem_i(obsSpaceData,OBS_ETM,jobs),tim_nstepobs)

          ! check if obs is outside of assimilation window when jstep = 1
          if (jstep.eq.1 .and.  &
               (stepobsIndex.lt.1.0 .or. stepObsIndex.gt.real(tim_nstepobs,8)) ) then
             write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS: Observation time outside assimilation window: ',  &
                  obs_headElem_i(obsSpaceData,OBS_DAT,jobs),obs_headElem_i(obsSpaceData,OBS_ETM,jobs)

             ! put the obs in the first time bin (it has to go somewhere!)
             stepObsIndex=1.0d0

             ! flag it as out of time domain and turn off its assimilation flag
             idata = obs_headElem_i(obsSpaceData,OBS_RLN,jobs)
             idatend = obs_headElem_i(obsSpaceData,OBS_NLV,jobs) + idata -1
             do jdata = idata, idatend
                call obs_bodySet_i(obsSpaceData,OBS_ASS,JDATA, 0)
             end do
             call obs_headSet_i(obsSpaceData,OBS_ST1,jobs,  &
                  ibset( obs_headElem_i(obsSpaceData,OBS_ST1,jobs), 05))
          end if

          if ( nint(stepObsIndex) == jstep ) then

             nobs(jstep) = nobs(jstep) + 1
             notag(nobs(jstep),jstep) = jobs

             !- Get LatLon of observation location
             lat_r8=obs_headElem_r(obsSpaceData,OBS_LAT,jobs)
             lon_r8=obs_headElem_r(obsSpaceData,OBS_LON,jobs)
             lat_r4=real(lat_r8)
             lon_r4=real(lon_r8)
             if (lon_r4.lt.0.0         ) lon_r4 = lon_r4 + 2.0*MPC_PI_R4
             if (lon_r4.ge.2.*MPC_PI_R4) lon_r4 = lon_r4 - 2.0*MPC_PI_R4

             lat_deg_r4=lat_r4 * MPC_DEGREES_PER_RADIAN_R4 ! Radian To Degree
             lon_deg_r4=lon_r4 * MPC_DEGREES_PER_RADIAN_R4

             !
             !- Find the position in the analysis grid
             !
             ierr = gdxyfll( hco_anl % EZscintID, xpos_r4, ypos_r4, &
                  lat_deg_r4, lon_deg_r4, 1)

             !- Test if the obs is outside the analysis grid
             if ( xpos_r4 < xposLowerBoundAnl_r4  .or. &
                  xpos_r4 > xposUpperBoundAnl_r4  .or. &
                  ypos_r4 < 1.0                   .or. &
                  ypos_r4 > real(hco_anl % nj) ) then

                if ( hco_anl % global ) then
                   ! Modify latitude if we have an observation at or near the poles
                   write(*,*) ''
                   write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS: Moving OBS inside the GLOBAL ANALYSIS grid, ', jobs
                   write(*,*) '  true position : ', lat_deg_r4, lon_deg_r4, ypos_r4, xpos_r4

                   !- Move the observation to the nearest grid point
                   if ( ypos_r4 < 1.0 )                ypos_r4 = 1.0
                   if ( ypos_r4 > real(hco_anl % nj) ) ypos_r4 = real(hco_anl % nj)

                   ierr = gdllfxy( hco_anl % EZscintID, &    ! IN
                        lat_deg_r4, lon_deg_r4, & ! OUT
                        xpos_r4, ypos_r4, 1)      ! IN

                   write(*,*) '  new  position : ', lat_deg_r4, lon_deg_r4, ypos_r4, xpos_r4

                   lat_r8 = real(lat_deg_r4,8) * MPC_RADIANS_PER_DEGREE_R8
                   lon_r8 = real(lon_deg_r4,8) * MPC_RADIANS_PER_DEGREE_R8
                   call obs_headSet_r(obsSpaceData,OBS_LAT,jobs, lat_r8) ! IN
                   call obs_headSet_r(obsSpaceData,OBS_LON,jobs, lon_r8) ! IN

                else
                   ! The observation is outside the domain
                   ! In LAM Analysis mode we must discard this observation
                   write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS: Rejecting OBS outside the LAM ANALYSIS grid domain, ', jobs
                   write(*,*) '  position : ', lat_deg_r4, lon_deg_r4, ypos_r4, xpos_r4

                   idata   = obs_headElem_i(obsSpaceData,OBS_RLN,jobs)
                   idatend = obs_headElem_i(obsSpaceData,OBS_NLV,jobs) + idata -1
                   do jdata = idata, idatend
                      call obs_bodySet_i(obsSpaceData,OBS_ASS,JDATA, 0)
                   end do
                   call obs_headSet_i(obsSpaceData,OBS_ST1,jobs,  &
                        ibset( obs_headElem_i(obsSpaceData,OBS_ST1,jobs), 05))
                end if

             end if

             !- Convert to rotated grid if needed
             if (hco_anl % rotated) then
                call uvr_RotateLatLon( lat_rot, lon_rot,       & ! OUT (radians)
                                       lat_r8,                 & ! IN  (radians)
                                       lon_r8,                 & ! IN  (radians)
                                       'ToLatLonRot')            ! IN
             else
                lat_rot = lat_r8
                lon_rot = lon_r8
             end if

             !- Store the above 3 pairs of values in column structure
             ypos_r8 = real(ypos_r4,8)
             xpos_r8 = real(xpos_r4,8)
             call col_setLatLon( columnhr, jobs, lat_r8, lon_r8,   & ! IN
                  ypos_r8, xpos_r8, lat_rot, lon_rot ) ! IN

             !
             !- Find the position in the trial field grid
             !
             ierr=gdxyfll(EZscintID_trl, xpos_r4, ypos_r4, lat_deg_r4, lon_deg_r4, 1)

             if ( xpos_r4 >= 1.0 .and. xpos_r4 <= real(ni_trl) .and.  &
                  ypos_r4 >= 1.0 .and. ypos_r4 <= real(nj_trl) ) then

                dlonfld(nobs(jstep)) = lon_r8
                dlatfld(nobs(jstep)) = lat_r8
                if (dlonfld(nobs(jstep)).lt.0.0d0)  &
                     dlonfld(nobs(jstep)) = dlonfld(nobs(jstep)) +  &
                     2*MPC_PI_R8
                if (dlonfld(nobs(jstep)).ge.2.0d0*MPC_PI_R8)  &
                     dlonfld(nobs(jstep)) =dlonfld(nobs(jstep)) -  &
                     2*MPC_PI_R8
                dlonfld(nobs(jstep))=dlonfld(nobs(jstep))*MPC_DEGREES_PER_RADIAN_R8
                dlatfld(nobs(jstep))=dlatfld(nobs(jstep))*MPC_DEGREES_PER_RADIAN_R8

             else
                ! The observation is outside the domain
                ! With a LAM trial field we must discard this observation
                write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS: Rejecting OBS outside the TRIAL field domain, ', jobs
                write(*,*) '  position : ', lat_deg_r4, lon_deg_r4, ypos_r4, xpos_r4

                idata   = obs_headElem_i(obsSpaceData,OBS_RLN,jobs)
                idatend = obs_headElem_i(obsSpaceData,OBS_NLV,jobs) + idata -1
                do jdata = idata, idatend
                   call obs_bodySet_i(obsSpaceData,OBS_ASS,JDATA, 0)
                end do
                call obs_headSet_i(obsSpaceData,OBS_ST1,jobs,  &
                     ibset( obs_headElem_i(obsSpaceData,OBS_ST1,jobs), 05))

                ! However, we must assigned a realistic lat-lon to this point
                ! to avoid problem later in Hx computation.
                ierr=gdllfxy(EZscintID_trl, lat_deg_r4, lon_deg_r4, real(ni_trl)/2.0,  &
                     real(nj_trl)/2.0, 1) ! Middle of the domain
                dlonfld(nobs(jstep)) = real(lon_deg_r4,8)
                dlatfld(nobs(jstep)) = real(lat_deg_r4,8)
             end if

          end if
       end do ! jobs

       ! gather and compute the max number of obs over all processors for each timestep
       call rpn_comm_allreduce(nobs(jstep),nobs_maxmpiglobal(jstep),1,  &
            'MPI_INTEGER','MPI_MAX','GRID',ierr)
       call rpn_comm_allgather(nobs(jstep),1,'mpi_integer',       &
            nobs_mpiglobal(jstep,:),1,'mpi_integer', &
            'GRID',ierr)
       ! gather lon-lat of observations from all processors
       call rpn_comm_allgather(dlonfld,numColumn_maxmpiglobal,'mpi_double_precision',       &
            dlonfld_mpiglobal,numColumn_maxmpiglobal,'mpi_double_precision', &
            'GRID',ierr)
       call rpn_comm_allgather(dlatfld,numColumn_maxmpiglobal,'mpi_double_precision',       &
            dlatfld_mpiglobal,numColumn_maxmpiglobal,'mpi_double_precision', &
            'GRID',ierr)

       zig1 = 0.0D0
       zig2 = 0.0D0
       zig3 = 1.0D0
       zig4 = 1.0D0
       call utl_cxgaig('L',ig1obs,ig2obs,ig3obs,ig4obs,zig1,zig2,zig3,zig4)

       do jlatlontile = 1,mpi_nprocs
          if (nobs_mpiglobal(jstep,jlatlontile).gt.0) then
             nobsgid_mpiglobal(jstep,jlatlontile) = utl_ezgdef(nobs_mpiglobal(jstep,jlatlontile),  &
                  1,'Y','L',ig1obs,ig2obs,ig3obs,ig4obs,  &
                  dlonfld_mpiglobal(1:nobs_mpiglobal(jstep,jlatlontile),jlatlontile),  &
                  dlatfld_mpiglobal(1:nobs_mpiglobal(jstep,jlatlontile),jlatlontile))
          else
             !write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS: NO OBS found for this time/lat bin =',jstep,jlatlontile
             nobsgid_mpiglobal(jstep,jlatlontile) = -999
          end if
       end do

    end do ! jstep

    !
    !     reading 2D fields
    !
    do jvar=1,vnl_numvarmax2D
       if (.not.gsv_varExist(varName=vnl_varNameList2D(jvar))) cycle

       call readTrialField(varInterphr_M,varInterphr_VV,vnl_varNameList2D(jvar),'SF')

       if (numColumns.gt.0) then       
          if (vnl_varNameList2D(jvar).eq.'P0  ') then
             varInterphr_M(:,:)=varInterphr_M(:,:)*MPC_PA_PER_MBAR_R8
          endif
          call col_fillmvo(columnhr,varInterphr_M,vnl_varNameList2D(jvar))
       endif

    enddo

    !
    !     Derive the pressure fields at observation points
    !      
    if (numColumns > 0 .and. col_varExist('P0')) then
       call col_calcPressure(columnhr)

       do jlev = 1,col_getNumLev(columnhr,'MM')
          if (mpi_myid.eq.0) write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS: jlev, col_getPressure(COLUMNHR,jlev,1,MM) = ',  &
               jlev,col_getPressure(columnhr,jlev,1,'MM')
       end do
       do jlev = 1,col_getNumLev(columnhr,'TH')
          if (mpi_myid.eq.0) write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS: jlev, col_getPressure(COLUMNHR,jlev,1,TH) = ',  &
               jlev,col_getPressure(columnhr,jlev,1,'TH')
       end do
       if (mpi_myid.eq.0) write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS: surface Pressure=',col_getElem(columnhr,1,1,'P0')

    end if

    !      
    !     Variable GZ qui se trouve sur les niveaux momentum et thermodynamiques
    !
    write(*,*)' ----- Initializing GZ ----'

    !
    !     Lire les GZ des niveaux Momentum
    !
    call readTrialField(varInterphr_M,varInterphr_VV,'GZ  ','MM',noUpperGZ)

    if (numColumns.gt.0) then       
       varInterphr_M(:,:)=varInterphr_M(:,:)*10.0d0*RG
       call col_fillmvo(columnhr,varInterphr_M,'GZ  ','MM')
       if (mpi_myid.eq.0) write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS:GZ_M'
       do jlev = 1,nlevtrl_M
          if (mpi_myid.eq.0) write(*,*) 'GZ,',jlev,varInterphr_M(jlev,1)
       enddo
    endif

    !
    !     Lire les GZ des niveaux Thermodynamique
    !
    call readTrialField(varInterphr_T,varInterphr_VV,'GZ  ','TH',noUpperGZ)

    if (numColumns.gt.0) then       
       varInterphr_T(:,:)=varInterphr_T(:,:)*10.0d0*RG
       call col_fillmvo(columnhr,varInterphr_T,'GZ  ','TH')
       if (mpi_myid.eq.0) write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS:GZ_TH'
       do jlev = 1,nlevtrl_T
          if (mpi_myid.eq.0) write(*,*)'GZ,',jlev,varInterphr_T(jlev,1)
       enddo
    endif

    !
    !     Now all of the other 3D variables
    !
    do jvar=1, vnl_numvarmax3D

       if (.not.gsv_varExist(varName=vnl_varNameList3D(jvar))) cycle

       if (vnl_varNameList3D(jvar).eq.'VV') cycle  ! to avoid cycle for VV

       select case ( vnl_varNameList3D(jvar) )
          !
          !       Variables sur les niveaux momentum
          !
       case ('UU')
          write(*,*)' ----- Initializing UU and VV  ----'

          call readTrialField(varInterphr_M,varInterphr_VV,'UV  ','MM')

          if (numColumns.gt.0) then       

             if (mpi_myid.eq.0) write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS: UU ,nlev= ',nlevtrl_M
             do jlev = 1,nlevtrl_M
                if (mpi_myid.eq.0) write(*,*) 'UU',jvar,jlev,varInterphr_M(jlev,1)
             enddo
             if (mpi_myid.eq.0) write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS: VV ,nlev= ',nlevtrl_M
             do jlev = 1,nlevtrl_M
                if (mpi_myid.eq.0) write(*,*) 'VV',jvar,jlev,varInterphr_VV(jlev,1)
             enddo

             call col_fillmvo(columnhr,varInterphr_M,'UU  ')
             call col_fillmvo(columnhr,varInterphr_VV,'VV  ')

             ! conversion from knots to m/s
             do jobs=1,numColumns
                column_ptr => col_getColumn(columnhr,jobs,'UU')
                do jlev=1,col_getNumLev(columnhr,'MM')
                   column_ptr(jlev)=column_ptr(jlev)*MPC_M_PER_S_PER_KNOT_R8
                enddo
                column_ptr => col_getColumn(columnhr,jobs,'VV')
                do jlev=1,col_getNumLev(columnhr,'MM')
                   column_ptr(jlev)=column_ptr(jlev)*MPC_M_PER_S_PER_KNOT_R8
                enddo
             enddo

          endif

          !
          !       Variable sur les niveaux thermodynamiques
          !
       case ('TT','HU')
          write(*,*)' ----- Initializing ',vnl_varNameList3D(jvar),' ----'

          call readTrialField(varInterphr_T,varInterphr_VV,vnl_varNameList3D(jvar),'TH')

          if (numColumns.gt.0) then       

             if (mpi_myid.eq.0) write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS:',vnl_varNameList3D(jvar)
             do jlev = 1,nlevtrl_T
                if (mpi_myid.eq.0) write(*,*) trim(vnl_varNameList3D(jvar)),',',jlev,varInterphr_T(jlev,1)
             enddo

             call col_fillmvo(columnhr,varInterphr_T,vnl_varNameList3D(jvar))

             if (vnl_varNameList3D(jvar).eq.'TT  ') then
                ! conversion from Celcius to Kelvin
                do jobs=1,numColumns
                   column_ptr => col_getColumn(columnhr,jobs,'TT')
                   do jlev=1,col_getNumLev(columnhr,'TH')
                      column_ptr(jlev)=column_ptr(jlev)+MPC_K_C_DEGREE_OFFSET_R8
                   enddo
                enddo
             elseif (vnl_varNameList3D(jvar).eq.'HU  ') then
                !!! conversion from specific humidity to log(humidity)
               ! Imposing a minimum value for HU
                do jobs=1,numColumns
                   column_ptr => col_getColumn(columnhr,jobs,'HU')
                   do jlev=1,col_getNumLev(columnhr,'TH')
                      column_ptr(jlev)=max(column_ptr(jlev),col_rhumin)  !log(max(column_ptr(jlev),col_rhumin))
                   enddo
                enddo
             endif

          endif

       case default

          call readTrialField(varInterphr_T,varInterphr_VV,vnl_varNameList3D(jvar),vnl_varLevelFromVarname(vnl_varNameList3D(jvar)))

          if(numColumns.gt.0) then

             if(mpi_myid.eq.0) write(*,*) 'inn_setupBackgroundColumns:',vnl_varNameList3D(jvar)
             do jlev = 1,nlevtrl_T
                if(mpi_myid.eq.0) write(*,*) trim(vnl_varNameList3D(jvar)),',',jlev,varInterphr_T(jlev,1)
             enddo
             
             call col_fillmvo(columnhr,varInterphr_T,vnl_varNameList3D(jvar))

          end if

       end select
    enddo

    if (col_varExist('TT') .and. col_varExist('HU') .and. col_varExist('P0')) then
       !
       !- Using T, q and PS to compute GZ for columnhr
       !
       if (noUpperGZ .and. (col_getNumLev(columnhr,'MM') > 1) ) then
          write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS: No upper level GZ found, computing'
          call tt2phi(columnhr)
       endif
    else
       write(*,*) ':inn_setupBackgroundColumns  GZ TLM calcs not generated since TT, HU and P0 not all present'
    end if

    !
    !- Close the files
    !
    do jstep = (1+mpi_myid), tim_nStepObs, mpi_nprocs
       ierr=fstfrm(nultrl(jstep))  
       ierr=fclos(nultrl(jstep))  
       ierr = ram_remove(trialfile(jstep))
    enddo

    !
    !- Deallocate the local arrays
    !
    if (numColumns.gt.0) then       
       deallocate(notag)
       deallocate(varInterphr_T)
       deallocate(varInterphr_M)
       deallocate(varInterphr_VV)
    endif
    deallocate(datestamplist)
    deallocate(nobs,nobs_maxmpiglobal)
    deallocate(nultrl)
    deallocate(idate)
    deallocate(itime)
    deallocate(dlonfld)
    deallocate(dlatfld)
    deallocate(dlonfld_mpiglobal)
    deallocate(dlatfld_mpiglobal)
    deallocate(nobsgid_mpiglobal)
    deallocate(nobs_mpiglobal)

    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) ' '
    write(*,*) '-------- Leaving INN_SETUPBACKGROUNDCOLUMNS ---------'
    write(*,*) ' '
    call tmg_stop(10)

  contains

    subroutine readTrialField(varInterphr_MT,varInterphr_VV,varName_in,varLevel,noUpperGZ_opt)
      !
      ! s/r readTrialField
      !
      !     Author  : M. Buehner, Dec 2012
      !
      !     Purpose: Read and interpolate all levels/time steps for a single variable of trial field
      !
      !     Revisions:
      !               Y. Rochon, ARQI/AQRD, Feb 2016
      !               - Added calls for constituent variable transformations
      !                 (routine chm_apply_2dfieldr4_transform)
      !
      implicit none
      character(len=*) :: varName_in
      character(len=*) :: varLevel
      character(len=4) :: varName
      logical, optional :: noUpperGZ_opt
      real*8 :: varInterphr_MT(:,:),varInterphr_VV(:,:)
      real*4, allocatable :: varTrial_r4(:,:),varTrial_VV_r4(:,:)
      real*4, allocatable :: varTrial_zero_r4(:,:)
      real*4, allocatable :: varInterp_r4(:,:,:),varInterp_VV_r4(:,:,:)
      real*4, allocatable :: varInterp2_r4(:),varInterp2_VV_r4(:)
      real*4, allocatable :: varInterp_recv_r4(:,:),varInterp_recv_VV_r4(:,:)
      integer :: nlevel,nsize,iip1,pe_send,pe_recv,tag,tag2
      integer :: fstlir,ezdefset

      integer, parameter :: maxnumkeys = 500
      integer :: numkeys, keys(maxnumkeys)
      integer :: key, fstinf, fstprm, fstinl, EZscintID, ezqkdef
      integer :: ni, nj, nk
      integer :: dateo, deet, npas, nbits, datyp
      integer :: ip1, ip2, ip3, swa, lng, dltf, ubc
      integer :: extra1, extra2, extra3
      integer :: ig1, ig2, ig3, ig4

      if(mpi_doBarrier) call rpn_comm_barrier('GRID',ierr)
      call tmg_start(91,'READTRIALFIELD')
      !
      ! Determine the type and number of vertical levels
      !
      if (trim(varName_in).eq.'UV') then
         varName='UU  '
      else
         varName=varName_in
      endif

      nlevel = col_getNumLev(columnhr,varLevel)

      if ( (1+mpi_myid) <= tim_nStepObs ) then
        ! already open (possibly on ram disk)
        nultrl_forEZ = nultrl(1+mpi_myid)
      else
        ! get from working directory (disk)
        trialfile_forEZ = 'trlm_01'
        trialfile_forEZ = ram_fullWorkingPath(trialfile_forEZ)
        nultrl_forEZ = 0
        ierr = fnom(nultrl_forEZ, trim(trialfile_forEZ),'RND+OLD+R/O',0)
        ierr = fstouv(nultrl_forEZ,'RND+OLD')
      endif

      !
      ! Check if too few GZ levels are in file, if so then just read surface
      !
      if (trim(varName).eq.'GZ' .and. present(noUpperGZ_opt)) then
         ierr = fstinl(nultrl_forEZ, ni, nj, nk, -1, ' ', -1, -1, -1, &
              ' ','GZ',keys, numkeys, maxnumkeys)
         write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS: nlevel =  ', nlevel, ', numkeys = ', numkeys
         if (numkeys.lt.nlevel) then
            write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS: Too few GZ levels found, will only read surface'
            noUpperGZ_opt = .true.
         else
            write(*,*) 'INN_SETUPBACKGROUNDCOLUMNS: sufficient levels of GZ found, will read all'
            noUpperGZ_opt = .false.
         endif
      endif

      !
      ! Determine grid size and EZSCINT ID
      !
      dateo  = -1
      cletiket = ' '
      ip1    = -1
      ip2    = -1
      ip3    = -1
      cltypvar = ' '

      key = fstinf( nultrl_forEZ,                            & ! IN
           ni, nj, nk,                                       & ! OUT
           dateo, cletiket, ip1, ip2, ip3, cltypvar, varName ) ! IN

      if (key < 0) then
         write(6,*)
         write(6,*) 'INN_SETUPBACKGROUNDCOLUMNS: Unable to find trial field = ',varName
         call utl_abort('INN_SETUPBACKGROUNDCOLUMNS')
      end if

      ierr = fstprm( key,                                              & ! IN
           dateo, deet, npas, ni, nj, nk, nbits,              & ! OUT
           datyp, ip1, ip2, ip3, cltypvar, varName, cletiket, & ! OUT
           clgrtyp, ig1, ig2, ig3,                            & ! OUT
           ig4, swa, lng, dltf, ubc, extra1, extra2, extra3 )   ! OUT

      EZscintID  = ezqkdef( ni, nj, clgrtyp, ig1, ig2, ig3, ig4, nultrl_forEZ )   ! IN

      allocate(varTrial_r4(ni,nj))
      allocate(varTrial_VV_r4(ni,nj))
      allocate(varTrial_zero_r4(ni,nj))
      allocate(varInterp_r4(maxval(nobs_mpiglobal),tim_nStepObs,mpi_nprocs))
      allocate(varInterp_VV_r4(maxval(nobs_mpiglobal),tim_nStepObs,mpi_nprocs))
      allocate(varInterp2_r4(maxval(nobs_mpiglobal)))
      allocate(varInterp2_VV_r4(maxval(nobs_mpiglobal)))
      allocate(varInterp_recv_r4(maxval(nobs_mpiglobal),tim_nStepObs))
      allocate(varInterp_recv_VV_r4(maxval(nobs_mpiglobal),tim_nStepObs))
      varTrial_zero_r4(:,:) = 0.0

      ! in the case that not all variable have the same etiket or typvar
      ! (this is necessary for extra 3d-var done before gen_coeff)
      cletiket='            '
      cltypvar='  '
      
      LEV_LOOP: do jlev = 1, nlevel

         STEP_LOOP1: do jstep = (1+mpi_myid), tim_nStepObs, mpi_nprocs

            if (nobs_maxmpiglobal(jstep) <= 0) cycle STEP_LOOP1

            if (trim(varName).eq.'GZ'.and.jlev.eq.nlevel) then
               ! use surface level IP1 for GZ (essential for Vcode=5005)
               IIP1 = vco_trl%ip1_sfc
            elseif (varLevel.eq.'MM') then
               IIP1 = vco_trl%ip1_M(jlev)
            elseif (varLevel.eq.'TH') then
               IIP1 = vco_trl%ip1_T(jlev)
            elseif (varLevel.eq.'SF') then
               IIP1 = -1
            else
               call utl_abort('INN_SETUPBACKGROUNDCOLUMNS: unknown varLevel')
            endif

            if (trim(varName).eq.'GZ'.and.present(noUpperGZ_opt)) then
               if (noUpperGZ_opt.and.jlev.ne.nlevel) then
                  ! do not try to read GZ above surface if it does not exist
                  varTrial_r4(:,:) = 0.0
               else
                  ierr=fstlir(varTrial_r4(:,:),nultrl(jstep),ni,nj,nk,  &
                       datestamplist(jstep) ,cletiket,iip1,-1,-1,  &
                       cltypvar,varName)
               endif
            else
               ierr=fstlir(varTrial_r4(:,:),nultrl(jstep),ni,nj,nk,  &
                    datestamplist(jstep) ,cletiket,iip1,-1,-1,  &
                    cltypvar,varName)
            endif

            if (ierr.lt.0)then
               write(*,2001) varName,iip1,idate(jstep),itime(jstep)
               call utl_abort('INN_SETUPBACKGROUNDCOLUMNS: Problem with background file')
            end if

            if (vnl_varKindFromVarname(varName).eq.'CH') &
                 call chm_apply_2dfieldr4_transform(vnl_varnumFromVarName(varName),varName,jlev,jstep,varTrial_r4)

            if (varName.eq.'UU') then
               ierr=fstlir(varTrial_VV_r4(:,:),nultrl(jstep),ni,nj,nk,  &
                    datestamplist(jstep) ,cletiket,iip1,-1,-1,  &
                    cltypvar,'VV')
               if (ierr.lt.0)then
                  write(*,2001) 'VV',iip1,idate(jstep),itime(jstep)
                  call utl_abort('INN_SETUPBACKGROUNDCOLUMNS: Problem with background file')
               end if
            endif

            ! Interpolate to mpiglobal set of columns for a subset of levels
            do jlatlontile = 1,mpi_nprocs
               if (nobs_mpiglobal(jstep,jlatlontile).gt.0) then
                  iset = ezdefset(nobsgid_mpiglobal(jstep,jlatlontile),EZscintID)
                  if (trim(varName).eq.'UU') then

                     ierr = utl_ezuvint(varInterp_r4(1:nobs_mpiglobal(jstep,jlatlontile),jstep,jlatlontile),  &
                                        varInterp_VV_r4(1:nobs_mpiglobal(jstep,jlatlontile),jstep,jlatlontile),  &
                                        varTrial_r4,varTrial_zero_r4, interpDegree='LINEAR')
                     ierr = utl_ezuvint(varInterp2_r4(1:nobs_mpiglobal(jstep,jlatlontile)),  &
                                        varInterp2_VV_r4(1:nobs_mpiglobal(jstep,jlatlontile)),  &
                                        varTrial_zero_r4,varTrial_VV_r4, interpDegree='LINEAR')

                     varInterp_r4(1:nobs_mpiglobal(jstep,jlatlontile),jstep,jlatlontile) =  &
                          real(varInterp_r4(1:nobs_mpiglobal(jstep,jlatlontile),jstep,jlatlontile),8) +  &
                          real(varInterp2_r4(1:nobs_mpiglobal(jstep,jlatlontile)),8)
                     varInterp_VV_r4(1:nobs_mpiglobal(jstep,jlatlontile),jstep,jlatlontile) =  &
                          real(varInterp_VV_r4(1:nobs_mpiglobal(jstep,jlatlontile),jstep,jlatlontile),8) +  &
                          real(varInterp2_VV_r4(1:nobs_mpiglobal(jstep,jlatlontile)),8)

                     ! This could replace the code above, but results are changed, so more testing needed
                     !ierr = utl_ezuvint2(varInterp(1:nobs_mpiglobal(jstep,jlatlontile),jstep,jlatlontile),  &
                     !                    varInterp_VV(1:nobs_mpiglobal(jstep,jlatlontile),jstep,jlatlontile),  &
                     !                    varTrial_r4,varTrial_VV_r4,  &
                     !                    nobs_mpiglobal(jstep,jlatlontile),ni*nj)
                  else
                     ierr = utl_ezsint(varInterp_r4(1:nobs_mpiglobal(jstep,jlatlontile),jstep,jlatlontile),  &
                                       varTrial_r4, interpDegree='LINEAR')
                  endif
               endif
            enddo

         enddo STEP_LOOP1

         TILE_LOOP: do jlatlontile = 1,mpi_nprocs
            STEP_LOOP2: do jstep = 1, tim_nStepObs

               pe_send = mod(jstep-1,mpi_nprocs)
               pe_recv = jlatlontile-1
               tag  = pe_recv*500 + pe_send
               tag2 = pe_recv*500 + pe_send + 1000000

               if ( mpi_myid == pe_recv ) then
                  varInterp_recv_r4(:,jstep) = 0.0
                  varInterp_recv_VV_r4(:,jstep) = 0.0
               endif

               if ( nobs_mpiglobal(jstep,jlatlontile) <= 0 ) cycle STEP_LOOP2

               if (pe_send == pe_recv) then
                  if (mpi_myid == pe_send) then
                     varInterp_recv_r4(1:nobs_mpiglobal(jstep,jlatlontile),jstep) =  &
                          varInterp_r4(1:nobs_mpiglobal(jstep,jlatlontile),jstep,jlatlontile)
                     if (trim(varName) == 'UU') then
                        varInterp_recv_VV_r4(1:nobs_mpiglobal(jstep,jlatlontile),jstep) =  &
                             varInterp_VV_r4(1:nobs_mpiglobal(jstep,jlatlontile),jstep,jlatlontile)
                     endif
                  endif
               else
                  if (mpi_myid == pe_send) then
                     nsize=nobs_mpiglobal(jstep,jlatlontile)
                     call rpn_comm_send(varInterp_r4(1:nobs_mpiglobal(jstep,jlatlontile),jstep,jlatlontile), &
                          nsize,'mpi_real4',pe_recv,tag,'GRID',ierr)
                     if (trim(varName).eq.'UU') then
                        call rpn_comm_send(varInterp_VV_r4(1:nobs_mpiglobal(jstep,jlatlontile),jstep,jlatlontile), &
                             nsize,'mpi_real4',pe_recv,tag2,'GRID',ierr)
                     endif
                  endif

                  if (mpi_myid == pe_recv) then
                     nsize=nobs_mpiglobal(jstep,jlatlontile)
                     call rpn_comm_recv(varInterp_recv_r4(1:nobs_mpiglobal(jstep,jlatlontile),jstep), &
                          nsize,'mpi_real4',pe_send,tag,'GRID',status,ierr)
                     if (trim(varName) == 'UU') then
                        call rpn_comm_recv(varInterp_recv_VV_r4(1:nobs_mpiglobal(jstep,jlatlontile),jstep), &
                             nsize,'mpi_real4',pe_send,tag2,'GRID',status,ierr)
                     endif
                  endif
               endif

            enddo STEP_LOOP2
         enddo TILE_LOOP

         do jstep = 1, tim_nStepObs
            do jobs = 1, nobs(jstep)
               varInterphr_MT(jlev,notag(jobs,jstep)) = varInterp_recv_r4(jobs,jstep)
               if (trim(varName).eq.'UU') then
                  varInterphr_VV(jlev,notag(jobs,jstep)) = varInterp_recv_VV_r4(jobs,jstep)
               endif
            enddo
         enddo

      enddo LEV_LOOP

      deallocate(varTrial_r4,varTrial_VV_r4)
      deallocate(varTrial_zero_r4)
      deallocate(varInterp_r4,varInterp_VV_r4)
      deallocate(varInterp_recv_r4,varInterp_recv_VV_r4)
      deallocate(varInterp2_r4,varInterp2_VV_r4)

2001  format(1x,'INN_SETUPBACKGROUNDCOLUMNS: Problem finding variable',1x,a4,1x,'at level',  &
           i10,1x,', on',1x,i8,1x,'at',1x,i8.8,1x,'HHMMSSss')

      call tmg_stop(91)

    end subroutine readTrialField

  end subroutine inn_setupBackgroundColumns


  subroutine inn_computeInnovation(columnhr,obsSpaceData,beSilent_opt)
    !
    ! Initialise Observation Innovations using the nonlinear H
    !
    implicit none

    type(struct_columnData) :: columnhr
    type(struct_obs)        :: obsSpaceData
    logical, optional       :: beSilent_opt

    real(8) :: zjo,zjoraob,zjosatwind,zjosurfc
    real(8) :: zjosfcsf,zjosfcua,zjotov,zjoairep,zjosfcsc,zjoprof
    real(8) :: zjogpsro,zjogpsgb,zjosfcgp,zjochm
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
    if ( col_getNumLev(columnhr,'MM') > 1 ) call oop_vobslyrs(columnhr,obsSpaceData)
    !
    !
    !------ Calculate the innovations Y - H(Xb) and place
    !       the result in obsSpaceData in OBS_OMP column
    !
    !        RAOBS
    !------------------------------
    !
    call tmg_start(48,'NL_OBS_OPER')
    call oop_ppp_nl(columnhr,obsSpaceData,ZJORAOB,'UA')
    !
    !        AIREPS
    !--------------------------------
    call oop_ppp_nl(columnhr,obsSpaceData,ZJOAIREP,'AI')
    !
    !        SATWINDS
    !--------------------------------
    call oer_sw(columnhr,obsSpaceData)
    call oop_ppp_nl(columnhr,obsSpaceData,ZJOSATWIND,'SW')
    !
    !        SURFACE (SF, UA, SC AND GP FAMILIES)
    !-------------------------------
    call oop_sfc_nl(columnhr,obsSpaceData,ZJOSFCSF,'SF')
    call oop_sfc_nl(columnhr,obsSpaceData,ZJOSFCUA,'UA')
    call oop_sfc_nl(columnhr,obsSpaceData,ZJOSFCSC,'SC')
    call oop_sfc_nl(columnhr,obsSpaceData,ZJOSFCGP,'GP')
    ZJOSURFC = ZJOSFCUA + ZJOSFCSF + ZJOSFCSC + ZJOSFCGP
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
    call oop_zzz_nl(columnhr,obsSpaceData,ZJOPROF,'PR')
    !
    !        GPS - RADIO OCCULTATION
    !-------------------------------
    ZJOGPSRO=0.0D0
    if (obs_famExist(obsSpaceData,'RO',local_mpi=.true.)) then
       CALL filt_gpsro(columnhr,obsSpaceData)
       CALL oer_SETERRGPSRO(columnhr,obsSpaceData)
       call oop_gpsro_nl(columnhr,obsSpaceData,beSilent,ZJOGPSRO)
    end if
    !
    !        CH - CHEMICAL CONSTITUENTS
    !-------------------------------
    call oop_chm_nl(columnhr,obsSpaceData,zjochm)
    !
    !        GPS - GROUND-BASED ZENITH DELAY
    !-------------------------------
    !
    ZJOGPSGB=0.0D0
    if (obs_famExist(obsSpaceData,'GP',local_mpi=.true.)) then
      if (trim(innovationMode) == 'analysis' .or. trim(innovationMode) == 'FSO') then
        call oer_SETERRGPSGB(columnhr,obsSpaceData,lgpdata,.true.)
        if (lgpdata) call oop_gpsgb_nl(columnhr,obsSpaceData,beSilent,ZJOGPSGB,.true.)
      else
        call oer_SETERRGPSGB(columnhr,obsSpaceData,lgpdata,.false.)
        if (lgpdata) call oop_gpsgb_nl(columnhr,obsSpaceData,beSilent,ZJOGPSGB,.false.)
      end if
    end if

    call tmg_stop(48)
    !
    !=======================================================================
    ZJO =  ZJORAOB + ZJOAIREP + ZJOSATWIND + &
         ZJOSURFC + ZJOTOV + ZJOPROF + ZJOGPSRO + ZJOGPSGB + ZJOCHM
    !=======================================================================

    if ( .not.beSilent ) then
      write(*,*) 'Cost function values for this MPI task:'
      write(*,'(a15,f30.16)') 'JORAOB   = ',ZJORAOB
      write(*,'(a15,f30.16)') 'JOAIREP  = ',ZJOAIREP
      write(*,'(a15,f30.16)') 'JOSURFC  = ',ZJOSURFC
      write(*,'(a15,f30.16)') 'JOSFCSF  = ',ZJOSFCSF
      write(*,'(a15,f30.16)') 'JOSFCUA  = ',ZJOSFCUA
      write(*,'(a15,f30.16)') 'JOSFCSC  = ',ZJOSFCSC
      write(*,'(a15,f30.16)') 'JOSFCGP  = ',ZJOSFCGP
      write(*,'(a15,f30.16)') 'JOTOV    = ',ZJOTOV
      write(*,'(a15,f30.16)') 'JOSATWIND= ',ZJOSATWIND
      write(*,'(a15,f30.16)') 'JOPROF   = ',ZJOPROF
      write(*,'(a15,f30.16)') 'JOGPSRO  = ',ZJOGPSRO
      write(*,'(a15,f30.16)') 'JOGPSGB  = ',ZJOGPSGB
      write(*,'(a15,f30.16)') 'JOCHM    = ',ZJOCHM
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
      call mpi_allreduce_sumreal8scalar(ZJOGPSRO,'GRID')
      call mpi_allreduce_sumreal8scalar(ZJOGPSGB,'GRID')
      call mpi_allreduce_sumreal8scalar(ZJOCHM,'GRID')

      write(*,*) 'Cost function values summed for all MPI tasks:'
      write(*,'(a15,f30.16)') 'JORAOB   = ',ZJORAOB
      write(*,'(a15,f30.16)') 'JOAIREP  = ',ZJOAIREP
      write(*,'(a15,f30.16)') 'JOSURFC  = ',ZJOSURFC
      write(*,'(a15,f30.16)') 'JOSFCSF  = ',ZJOSFCSF
      write(*,'(a15,f30.16)') 'JOSFCUA  = ',ZJOSFCUA
      write(*,'(a15,f30.16)') 'JOSFCSC  = ',ZJOSFCSC
      write(*,'(a15,f30.16)') 'JOSFCGP  = ',ZJOSFCGP
      write(*,'(a15,f30.16)') 'JOTOV    = ',ZJOTOV
      write(*,'(a15,f30.16)') 'JOSATWIND= ',ZJOSATWIND
      write(*,'(a15,f30.16)') 'JOPROF   = ',ZJOPROF
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
    ! PURPOSE:
    !  Header indices are distributed following the chosen strategy,
    !  current options: "LIKESPLITFILES", "ROUNDROBIN", "LATLONTILES" or "LATLONTILESBALANCED".
    !
    implicit none

    type(struct_obs), intent(inout) :: obsSpaceData
    character(len=*), intent(in)    :: mpiStrategy

    type(struct_hco), pointer :: hco_anl

    real(8) :: lat_r8, lon_r8
    real    :: lat_r4, lon_r4
    real    :: xpos_r4, ypos_r4

    integer :: headerIndex
    integer :: latIndex, lonIndex
    integer :: ierr, nsize
    integer :: IP, IP_x, IP_y, IP2
    integer :: gdxyfll
    integer :: numHeaderFile,numHeaderFile_mpiglobal(mpi_nprocs)
    integer :: obsLoadTile_mpilocal(mpi_nprocs),obsLoadTile_mpiglobal(mpi_nprocs)
    integer :: obsLoadTile_mpiglobal_tmp(mpi_nprocs)
    integer :: obsLoadTile_sorted(mpi_nprocs), PE_sorted(mpi_nprocs), totalLoadToSend
    integer :: numHeaderTile_mpilocal(mpi_nprocs),numHeaderTile_mpiglobal(mpi_nprocs)
    integer :: PE_sender, PE_receiver, loadSent, thisLoadToSend, loadDifference
    integer, allocatable :: IPT_mpiglobal(:,:), IPT_mpilocal(:)
    integer, allocatable :: obsLoad_mpiglobal(:,:), obsLoad_mpilocal(:)
    integer :: codtyp, numtovs(mpi_nprocs), numir(mpi_nprocs), numtovs_mpiglobal(mpi_nprocs), numir_mpiglobal(mpi_nprocs)
    integer :: totalObsLoad_mpilocal(mpi_nprocs), totalObsLoad_mpiglobal(mpi_nprocs)
    integer :: get_max_rss
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

       write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
       call rpn_comm_allgather(numHeaderFile, 1, 'mpi_integer',  &
            numHeaderFile_mpiglobal, 1, 'mpi_integer', 'GRID', ierr)

       ! set PE index for the tile (where interpolation will be done)
       do headerIndex = 1, numHeaderFile
          ! compute grid index for each observation header
          lat_r8 = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)
          lon_r8 = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex)
          lat_r4 = real(lat_r8) * MPC_DEGREES_PER_RADIAN_R4
          lon_r4 = real(lon_r8) * MPC_DEGREES_PER_RADIAN_R4
          ierr = gdxyfll( hco_anl%EZscintID,   & ! IN 
               xpos_r4, ypos_r4,    & ! OUT
               lat_r4, lon_r4, 1 )    ! IN

          ! compute corresponding mpi task id for each observation
          latIndex = floor(ypos_r4)
          lonIndex = floor(xpos_r4)
          IP_y = mpivar_myidYfromLat(latIndex, hco_anl%nj)
          IP_x = mpivar_myidXfromLon(lonIndex, hco_anl%ni)
          IP = IP_x + IP_y*mpi_npex
          call obs_headSet_i(obsSpaceData,OBS_IPT,headerIndex,IP)
       end do

       ! set PE index for the column (where the obs operator will be done)

       ! initialize IPC = IPT
       do headerIndex = 1, numHeaderFile
          IP = obs_headElem_i(obsSpaceData,OBS_IPT,headerIndex)
          call obs_headSet_i(obsSpaceData,OBS_IPC,headerIndex,IP)
       end do

       ! make all of the IPT and obsLoad values visible on all processors
       nsize = maxval(numHeaderFile_mpiglobal(:))
       allocate(IPT_mpiglobal(nsize,mpi_nprocs))
       allocate(IPT_mpilocal(nsize))
       allocate(obsLoad_mpiglobal(nsize,mpi_nprocs))
       allocate(obsLoad_mpilocal(nsize))
       IPT_mpilocal(:) = -1
       obsLoad_mpilocal(:) = -1
       do headerIndex = 1, numHeaderFile
          IPT_mpilocal(headerIndex) = obs_headElem_i(obsSpaceData,OBS_IPT,headerIndex)
          obsLoad_mpilocal(headerIndex) = obsLoad(headerIndex)
       enddo
       call rpn_comm_allgather(IPT_mpilocal, nsize, 'mpi_integer',  &
            IPT_mpiglobal, nsize, 'mpi_integer', 'GRID', ierr)
       call rpn_comm_allgather(obsLoad_mpilocal, nsize, 'mpi_integer',  &
            obsLoad_mpiglobal, nsize, 'mpi_integer', 'GRID', ierr)

       write(*,*) ' '
       write(*,*) 'setObsMpiStrategy: do balancing according to approximate load'
       write(*,*) ' '

       ! compute total number of headers per PE on tiles
       numHeaderTile_mpilocal(:) = 0
       do headerIndex = 1, numHeaderFile
          IP = obs_headElem_i(obsSpaceData,OBS_IPT,headerIndex)
          numHeaderTile_mpilocal(IP+1) = numHeaderTile_mpilocal(IP+1) + 1
       end do
       call rpn_comm_allreduce(numHeaderTile_mpilocal,numHeaderTile_mpiglobal,mpi_nprocs,  &
            'MPI_INTEGER','MPI_SUM','GRID',ierr)
       if(mpi_myid == 0) write(*,*) 'setObsMpiStrategy: original numHeaderTile_mpiglobal     =',numHeaderTile_mpiglobal(:)

       ! compute total load per PE on tiles
       obsLoadTile_mpilocal(:) = 0
       do headerIndex = 1, numHeaderFile
          if(obsLoad(headerIndex).gt.0) then
             IP = obs_headElem_i(obsSpaceData,OBS_IPT,headerIndex)
             obsLoadTile_mpilocal(IP+1) = obsLoadTile_mpilocal(IP+1) + obsLoad(headerIndex)
          endif
       end do
       call rpn_comm_allreduce(obsLoadTile_mpilocal,obsLoadTile_mpiglobal,mpi_nprocs,  &
            'MPI_INTEGER','MPI_SUM','GRID',ierr)
       if(mpi_myid == 0) write(*,*) 'setObsMpiStrategy: original obsLoadTile_mpiglobal     =',obsLoadTile_mpiglobal(:)

       ! set obsLoad to zero if only 1 header on the tile to avoid sending it
       do IP = 0, mpi_nprocs-1
          if(numHeaderTile_mpiglobal(IP+1)==1) then
             obsLoadTile_mpiglobal(IP+1) = 0
             if(mpi_myid == 0) write(*,*) 'setObsMpiStrategy: setting the obsLoad to zero, ', & 
                  'since only 1 header index on the tile with myid = ',IP
          endif
       enddo

       ! sort the list of loads per tile
       obsLoadTile_mpiglobal_tmp(:) = obsLoadTile_mpiglobal(:)
       do IP = 0, mpi_nprocs-1
          obsLoadTile_sorted(IP+1) = -1
          do IP2 = 0, mpi_nprocs-1
             if(obsLoadTile_mpiglobal_tmp(IP2+1).gt.obsLoadTile_sorted(IP+1)) then
                obsLoadTile_sorted(IP+1) = obsLoadTile_mpiglobal_tmp(IP2+1)
                PE_sorted(IP+1) = IP2
             endif
          enddo
          obsLoadTile_mpiglobal_tmp(PE_sorted(IP+1)+1) = -1
       enddo
       if(mpi_myid == 0) write(*,*) 'setObsMpiStrategy: sorted obsLoadTile_mpiglobal       =',obsLoadTile_sorted(:)
       if(mpi_myid == 0) write(*,*) 'setObsMpiStrategy: corresponding sorted PE_mpiglobal  =',PE_sorted(:)

       ! modify IPC for headers that need to be sent to balance load
       IP_LOOP: do IP = 0, (mpi_nprocs/2)-1
          PE_sender=PE_sorted(IP+1)
          PE_receiver=PE_sorted(mpi_nprocs-IP)
          loadDifference = obsLoadTile_mpiglobal(PE_sender+1) -  &
                           obsLoadTile_mpiglobal(PE_receiver+1)
          totalLoadToSend = floor(0.5 * loadDifference)
          if(mpi_myid == 0) write(*,*) 'setObsMpiStrategy: PE_sender, PE_receiver, totalLoadToSend = ', &
                                       PE_sender, PE_receiver, totalLoadToSend

          ! modify IPC value from PE_sender to PE_receiver for headers whose load adds up to totalLoadToSend
          loadSent=0
          do IP2 = 0, (mpi_nprocs-1)
             if( loadSent < totalLoadToSend ) then
                do headerIndex = 1, numHeaderFile_mpiglobal(IP2+1)
                   thisLoadToSend = obsLoad_mpiglobal(headerIndex,IP2+1)
                   if( IPT_mpiglobal(headerIndex,IP2+1) == PE_sender .and. &
                       thisLoadToSend > 0 .and. &
                       loadSent < totalLoadToSend ) then
                      ! check if sending this header would mean entire load difference is sent
                      if( (loadSent + thisLoadToSend) == loadDifference ) then
                        if( mpi_myid == 0 ) write(*,*) 'setObsMpiStrategy: prevent sending all of load difference!'
                        cycle IP_LOOP
                      end if
                      loadSent = loadSent + thisLoadToSend
                      if(mpi_myid.eq.IP2) then
                         ! header to send is currently on my PE
                         !write(*,*) 'setObsMpiStrategy: changing obs_ipc from ', &
                         !           obs_headElem_i(obsSpaceData,OBS_IPC,headerIndex),' to ',PE_receiver
                         call obs_headSet_i(obsSpaceData,OBS_IPC,headerIndex,PE_receiver)
                      endif
                   endif
                enddo
             endif
          enddo
       enddo IP_LOOP

       deallocate(IPT_mpiglobal)
       deallocate(IPT_mpilocal)
       deallocate(obsLoad_mpiglobal)
       deallocate(obsLoad_mpilocal)

       ! ** The rest is just diagnostics used when trying to improve the 
       ! ** formula for estimating the load - could be removed

       if (obs_famExist(obsSpaceData,'TO')) then

       ! count the number of tovs and IR observations for each file
       numtovs(:)=0
       numir(:)=0
       totalObsLoad_mpilocal(:)=0
       do headerIndex = 1, numHeaderFile
          codtyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
          if(tvs_isIdBurpTovs(codtyp) ) numtovs(mpi_myid+1)=numtovs(mpi_myid+1)+1
          if(tvs_isIdBurpInst(codtyp,'IASI') ) numir(mpi_myid+1)=numir(mpi_myid+1)+1
          if(tvs_isIdBurpInst(codtyp,'AIRS') ) numir(mpi_myid+1)=numir(mpi_myid+1)+1
          if(tvs_isIdBurpInst(codtyp,'CRIS') ) numir(mpi_myid+1)=numir(mpi_myid+1)+1
          totalObsLoad_mpilocal(mpi_myid+1) = totalObsLoad_mpilocal(mpi_myid+1) + obsLoad(headerIndex)
       enddo
       call rpn_comm_allreduce(numtovs,numtovs_mpiglobal,mpi_nprocs,  &
            'MPI_INTEGER','MPI_SUM','GRID',ierr)
       call rpn_comm_allreduce(numir,numir_mpiglobal,mpi_nprocs,  &
            'MPI_INTEGER','MPI_SUM','GRID',ierr)
       call rpn_comm_allreduce(totalObsLoad_mpilocal,totalObsLoad_mpiglobal,mpi_nprocs,  &
            'MPI_INTEGER','MPI_SUM','GRID',ierr)
       if(mpi_myid == 0) write(*,*) 'setObsMpiStrategy: number of FILE headers for all TOVS = ',numtovs_mpiglobal(:)
       if(mpi_myid == 0) write(*,*) '                   number of FILE headers for only IR  = ',numir_mpiglobal(:)
       if(mpi_myid == 0) write(*,*) '                   estimated total obsLoad             = ',totalObsLoad_mpiglobal(:)

       ! count the number of tovs and IR observations for each tile
       numtovs(:)=0
       numir(:)=0
       totalObsLoad_mpilocal(:)=0
       do headerIndex = 1, numHeaderFile
          IP = obs_headElem_i(obsSpaceData,OBS_IPT,headerIndex)
          codtyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
          if(tvs_isIdBurpTovs(codtyp) ) numtovs(IP+1)=numtovs(IP+1)+1
          if(tvs_isIdBurpInst(codtyp,'IASI') ) numir(IP+1)=numir(IP+1)+1
          if(tvs_isIdBurpInst(codtyp,'AIRS') ) numir(IP+1)=numir(IP+1)+1
          if(tvs_isIdBurpInst(codtyp,'CRIS') ) numir(IP+1)=numir(IP+1)+1
          totalObsLoad_mpilocal(IP+1) = totalObsLoad_mpilocal(IP+1) + obsLoad(headerIndex)
       enddo
       call rpn_comm_allreduce(numtovs,numtovs_mpiglobal,mpi_nprocs,  &
            'MPI_INTEGER','MPI_SUM','GRID',ierr)
       call rpn_comm_allreduce(numir,numir_mpiglobal,mpi_nprocs,  &
            'MPI_INTEGER','MPI_SUM','GRID',ierr)
       call rpn_comm_allreduce(totalObsLoad_mpilocal,totalObsLoad_mpiglobal,mpi_nprocs,  &
            'MPI_INTEGER','MPI_SUM','GRID',ierr)
       if(mpi_myid == 0) write(*,*) 'setObsMpiStrategy: number of TILE headers for all TOVS = ',numtovs_mpiglobal(:)
       if(mpi_myid == 0) write(*,*) '                   number of TILE headers for only IR  = ',numir_mpiglobal(:)
       if(mpi_myid == 0) write(*,*) '                   estimated total obsLoad             = ',totalObsLoad_mpiglobal(:)

       ! count the number of tovs and IR observations for columns
       numtovs(:)=0
       numir(:)=0
       totalObsLoad_mpilocal(:)=0
       do headerIndex = 1, numHeaderFile
          IP = obs_headElem_i(obsSpaceData,OBS_IPC,headerIndex)
          codtyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
          if(tvs_isIdBurpTovs(codtyp) ) numtovs(IP+1)=numtovs(IP+1)+1
          if(tvs_isIdBurpInst(codtyp,'IASI') ) numir(IP+1)=numir(IP+1)+1
          if(tvs_isIdBurpInst(codtyp,'AIRS') ) numir(IP+1)=numir(IP+1)+1
          if(tvs_isIdBurpInst(codtyp,'CRIS') ) numir(IP+1)=numir(IP+1)+1
          totalObsLoad_mpilocal(IP+1) = totalObsLoad_mpilocal(IP+1) + obsLoad(headerIndex)
       enddo
       call rpn_comm_allreduce(numtovs,numtovs_mpiglobal,mpi_nprocs,  &
            'MPI_INTEGER','MPI_SUM','GRID',ierr)
       call rpn_comm_allreduce(numir,numir_mpiglobal,mpi_nprocs,  &
            'MPI_INTEGER','MPI_SUM','GRID',ierr)
       call rpn_comm_allreduce(totalObsLoad_mpilocal,totalObsLoad_mpiglobal,mpi_nprocs,  &
            'MPI_INTEGER','MPI_SUM','GRID',ierr)
       if(mpi_myid == 0) write(*,*) 'setObsMpiStrategy: number of COLS headers for all TOVS = ',numtovs_mpiglobal(:)
       if(mpi_myid == 0) write(*,*) '                   number of COLS headers for only IR  = ',numir_mpiglobal(:)
       if(mpi_myid == 0) write(*,*) '                   estimated total obsLoad             = ',totalObsLoad_mpiglobal(:)

       numtovs(:)=0
       do headerIndex = 1, numHeaderFile
          IP = obs_headElem_i(obsSpaceData,OBS_IPC,headerIndex)
          codtyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
          if(tvs_isIdBurpInst(codtyp,'AMSUA')) numtovs(IP+1)=numtovs(IP+1)+1
       enddo
       call rpn_comm_allreduce(numtovs,numtovs_mpiglobal,mpi_nprocs,  &
            'MPI_INTEGER','MPI_SUM','GRID',ierr)
       if(mpi_myid == 0) write(*,*) 'setObsMpiStrategy: number of COLS headers for AMSUA     = ',numtovs_mpiglobal(:)

       numtovs(:)=0
       do headerIndex = 1, numHeaderFile
          IP = obs_headElem_i(obsSpaceData,OBS_IPC,headerIndex)
          codtyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
          if(tvs_isIdBurpInst(codtyp,'AMSUB') .or. tvs_isIdBurpInst(codtyp,'MHS')) numtovs(IP+1)=numtovs(IP+1)+1
       enddo
       call rpn_comm_allreduce(numtovs,numtovs_mpiglobal,mpi_nprocs,  &
            'MPI_INTEGER','MPI_SUM','GRID',ierr)
       if(mpi_myid == 0) write(*,*) 'setObsMpiStrategy: number of COLS headers for AMSUB/MHS = ',numtovs_mpiglobal(:)

       numtovs(:)=0
       do headerIndex = 1, numHeaderFile
          IP = obs_headElem_i(obsSpaceData,OBS_IPC,headerIndex)
          codtyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
          if(tvs_isIdBurpInst(codtyp,'METEOSAT') .or. tvs_isIdBurpInst(codtyp,'GOES')) numtovs(IP+1)=numtovs(IP+1)+1
       enddo
       call rpn_comm_allreduce(numtovs,numtovs_mpiglobal,mpi_nprocs,  &
            'MPI_INTEGER','MPI_SUM','GRID',ierr)
       if(mpi_myid == 0) write(*,*) 'setObsMpiStrategy: number of COLS headers for GEORAD    = ',numtovs_mpiglobal(:)

       numtovs(:)=0
       do headerIndex = 1, numHeaderFile
          IP = obs_headElem_i(obsSpaceData,OBS_IPC,headerIndex)
          codtyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
          if(tvs_isIdBurpInst(codtyp,'SSMIS')) numtovs(IP+1)=numtovs(IP+1)+1
       enddo
       call rpn_comm_allreduce(numtovs,numtovs_mpiglobal,mpi_nprocs,  &
            'MPI_INTEGER','MPI_SUM','GRID',ierr)
       if(mpi_myid == 0) write(*,*) 'setObsMpiStrategy: number of COLS headers for SSMIS     = ',numtovs_mpiglobal(:)

      end if

    case default
       write(*,*)
       write(*,*) 'ERROR unknown mpiStrategy: ', trim(mpiStrategy)
       call utl_abort('setObsMpiStrategy')
    end select

  CONTAINS

    function obsLoad(headerIndex)
      implicit none
      integer :: headerIndex, codtyp, obsLoad
      integer :: bodyIndexBeg, bodyIndexEnd, bodyIndex

      ! this is a very simple recipe for estimating the computational load 
      ! based only on the codtyp - it works better than more complicated 
      ! approaches that tried to count the number of assimilated elements in
      ! each header and then multiply this by a factor based on the number 
      ! of RTTOV predictors for each type of radiance - probably should write
      ! a separate program to do a more objective evalution of this

      codtyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
      if(tvs_isIdBurpInst(codtyp,'IASI')) then
         obsLoad = 200
      elseif(tvs_isIdBurpInst(codtyp,'AIRS')) then
         obsLoad = 200
      elseif(tvs_isIdBurpInst(codtyp,'CRIS')) then
         obsLoad = 150
      elseif(tvs_isIdBurpTovs(codtyp)) then
         ! all other types of radiance obs
         obsLoad = 5
      else
         ! all non-radiance obs
         obsLoad = 1
      end if
    end function obsLoad

  end subroutine setObsMpiStrategy


  subroutine inn_perturbObs(obsSpaceData,numAnalyses,indexAnalysis,indexBatch,obs_column_index_src,obs_column_index_dest)
    !
    !Purpose:
    ! Perturb the innovation vector to simulate effect of observation uncertainty
    !
    ! WARNING: perturbations are not the same when MPI topology changes!!!
    !
    !Author  : M. Buehner, Dec, 2013
    !
    implicit none

    type(struct_obs) :: obsSpaceData
    integer :: numAnalyses,indexAnalysis,indexBatch,numPerturbations
    integer :: obs_column_index_src,obs_column_index_dest

    integer :: nrandseed,iseed,indexAnalysis2,indexBody,indexFamily,iass
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
       iass = obs_bodyElem_i(obsSpaceData,OBS_ASS,indexBody)
       if(iass.eq.1 .or. iass.eq.-1) then
          originalOmp = obs_bodyElem_r(obsSpaceData,obs_column_index_src,indexBody)
          call obs_bodySet_r(obsSpaceData,obs_column_index_dest,indexBody,originalOmp+obsPerturbations(indexBody,indexAnalysis))
       endif
    enddo

  end subroutine inn_perturbObs


end module innovation_mod
