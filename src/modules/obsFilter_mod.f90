
module obsFilter_mod
  ! MODULE obsFilter_mod (prefix='filt' category='5. Observation operators')
  !
  !:Purpose:  Various types of filters that are applied to the observations
  !           mostly to reject them so that they will not be assimilated.
  !
  use midasMpi_mod
  use codePrecision_mod
  use earthConstants_mod
  use obsSpaceData_mod
  use columnData_mod
  use bufr_mod
  use gps_mod
  use utilities_mod
  use varNameList_mod
  use physicsFunctions_mod
  use codtyp_mod
  use radialVelocity_mod
  use obsFlags_mod

  implicit none
  save
  private

  ! Public variables
  real(8), public, protected :: filt_rlimlvhu

  ! Public procedures
  public :: filt_setup, filt_topo, filt_suprep
  public :: filt_surfaceWind, filt_gpsro,  filt_backScatAnisIce, filt_iceConcentration, filt_radvel
  public :: filt_bufrCodeAssimilated, filt_getBufrCodeAssimilated, filt_nBufrCodeAssimilated
  public :: filt_getSfcBufferZoneCHheight

  integer, parameter :: nelemsMax = 30
  integer, parameter :: nflagsMax = 15
  integer :: filt_nelems, filt_nflags
  integer, target :: filt_nlist(nelemsMax)
  integer :: filt_nlistflg(nflagsMax)

  ! topographic rejection criteria
  integer, parameter :: numElem = 21
  real(8)            :: altDiffMax(numElem) =   & ! default values (in metres)
       (/     50.d0,    50.d0,     50.d0,      50.d0,     50.d0,    800.d0,    800.d0,  &
             800.d0,   800.d0,   1000.d0,      50.d0,     50.d0,     50.d0,     50.d0,  &
              50.d0,    50.d0,     50.d0,      50.d0,     50.d0,     50.d0,   1000.d0 /)
  integer, parameter :: elemList(numElem) =  &
       (/ BUFR_NEDS, BUFR_NEFS, BUFR_NEUS, BUFR_NEVS, BUFR_NESS, BUFR_NETS, BUFR_NEPS, &
          BUFR_NEPN, BUFR_NEGZ, BUFR_NEZD, BUFR_NEDD, BUFR_NEFF, BUFR_NEUU, BUFR_NEVV, &
          BUFR_NEES, BUFR_NETT, BUFR_NEAL, bufr_vis , bufr_logVis, bufr_gust, bufr_radvel /)

  integer, parameter :: nTopoFiltFam = 8
  character(len=2) :: filtTopoList(nTopoFiltFam) = '  '

  character(len=48) :: filterMode

  logical :: initialized = .false.

  ! Namelist variables:
  logical :: discardlandsfcwind          ! choose to reject surface wind obs over land
  real(8) :: surfaceBufferZone_Pres      ! height of buffer zone (in Pa) for rejecting obs near sfc
  real(8) :: surfaceBufferZone_Height    ! height of buffer zone (in m) for rejecting obs near sfc
  real(8) :: surfaceBufferZoneCH_Pres    ! height buffer zone (in Pa) for rejecting CH family obs near sfc
  real(8) :: surfaceBufferZoneCH_Height  ! height buffer zone (in m) for rejecting CH family obs near sfc
  logical :: useEnkfTopoFilt             ! choose to use simpler approach (originally in EnKF) for rejecting obs near sfc
  logical :: rejectGZforAnalysis         ! whether to reject geopotential height for analysis update

contains

  !--------------------------------------------------------------------------
  ! findElemIndex
  !-------------------------------------------------------------------------
  function findElemIndex(varNum) result(listIndex)
    implicit none

    ! Arguments:
    integer, intent(in) :: varNum
    ! Result:
    integer :: listIndex

    ! Locals:
    integer :: elemIndex

    listIndex = -1
    do elemIndex=1,numElem
       if(varNum == elemList(elemIndex)) listIndex = elemIndex
    end do

    if (listIndex == -1) then
       write(*,*) 'filterobs_mod-findElemIndex: WARNING: varNum value not found: ',varNum
    end if

  end function findElemIndex

  !--------------------------------------------------------------------------
  ! filt_setup
  !--------------------------------------------------------------------------
  subroutine filt_setup(filterMode_in)
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: filterMode_in

    ! Locals:
    integer :: ierr, elem, refValue, bitNumber, obsFamilyIndex, elemIndex
    integer :: flagIndex, elementIndex

    ! Namelist variables: (local)
    integer :: nelems                    ! MUST NOT BE INCLUDED IN NAMELIST!
    integer :: nlist(nelemsMax)          ! list of bufr element IDs to consider for assimilation
    integer :: nflags                    ! MUST NOT BE INCLUDED IN NAMELIST!
    integer :: nlistflg(nFLagsMax)       ! list of flag 'reference numbers' to use for rejecting obs
    integer :: nelems_altDiffMax         ! MUST NOT BE INCLUDED IN NAMELIST!
    integer :: list_altDiffMax(numElem)  ! list of bufr element IDs to apply maximum altitude difference
    character(len=2) :: list_topoFilt(nTopoFiltFam) ! list of obs family names for applying max altitude
    real(8) :: value_altDiffMax(numElem) ! value of maximum difference between model sfc and obs altitude
    real(8) :: rlimlvhu                  ! pressure level (in hPa) above which humidity (ES) obs are rejected

    namelist /namfilt/nelems, nlist, nflags, nlistflg, rlimlvhu, discardlandsfcwind, &
         nelems_altDiffMax, list_altDiffMax, value_altDiffMax, surfaceBufferZone_Pres, &
         surfaceBufferZone_Height, list_topoFilt, useEnkfTopoFilt, rejectGZforAnalysis, &
         surfaceBufferZoneCH_Pres,surfaceBufferZoneCH_Height

    filterMode = filterMode_in

    ! set default values for namelist variables
    nlist(:) =  MPC_missingValue_INT
    nelems = MPC_missingValue_INT
    nlistflg(:) = MPC_missingValue_INT
    nflags = MPC_missingValue_INT
    list_altDiffMax (:) = MPC_missingValue_INT
    value_altDiffMax(:) = MPC_missingValue_R8
    nelems_altDiffMax = MPC_missingValue_INT

    list_topoFilt(:) = '**'

    rlimlvhu = 300.d0
    discardlandsfcwind = .true.

    surfaceBufferZone_Pres   = 5000.0d0   ! default value in Pascals
    surfaceBufferZone_Height =  400.0d0   ! default value in Metres
    surfaceBufferZoneCH_Pres = 5000.0d0   ! default value in Pascals
    surfaceBufferZoneCH_Height = 400.0d0  ! default value in Metres

    useEnkfTopoFilt = .false.
    rejectGZforAnalysis = .true.

    call utl_tmg_start(181,'low-level--readNML')
    read(utl_flnml,nml=namfilt,iostat=ierr)
    if(ierr /= 0) call utl_abort('filt_setup: Error reading namelist! Hint: did you replace ltopofilt by list_topoFilt?')
    if(mmpi_myid == 0) write(*,nml=namfilt)
    call utl_tmg_stop(181)

    filt_rlimlvhu    = rlimlvhu

    if (nelems /= MPC_missingValue_INT) then
      call utl_abort('filt_setup: check NAMFILT namelist section; NELEMS should be removed')
    end if
    filt_nelems = 0
    do elementIndex = 1, nelemsMax
      if (nlist(elementIndex) /= MPC_missingValue_INT) then
        filt_nlist(elementIndex) = nlist(elementIndex)
        filt_nelems = filt_nelems + 1
      end if
    end do

    if (nflags /= MPC_missingValue_INT) then
      call utl_abort('filt_setup: check NAMFILT namelist section; NFLAGS should be removed')
    end if
    filt_nflags = 0
    do flagIndex = 1, nflagsMax
      if (nlistflg(flagIndex) /= MPC_missingValue_INT) then
        filt_nlistflg(flagIndex) = nlistflg(flagIndex)
        filt_nflags = filt_nflags + 1
      end if
    end do

    if(mmpi_myid == 0) then
      write(*,'(1X,"***********************************")')
      write(*,'(1X," ELEMENTS SELECTED FOR ASSIMILATION:",/)')
      write(*,'(1X,"***********************************")')
      do elem=1,filt_nelems
        write(*,'(15X,I5)') filt_nlist(elem)
      end do
      write(*,'(1X,"***********************************")')
      write(*,*) ' REJECT ELEMENTS WITH REJECT FLAG '
      do flagIndex = 1, filt_nflags
        refValue = filt_nlistflg(flagIndex)
        bitNumber = flg_refValueToBitNum(refValue)
        write(*,*) 'refValue, bitNumber, description = ', refValue, bitNumber,' ',flg_descrip(bitNumber)
      end do
      write(*,'(1X,"***********************************")')
    end if

    !
    !- Set values for altDiffMax
    !
    if (nelems_altDiffMax /= MPC_missingValue_INT) then
      call utl_abort('filt_setup: check namelist section NAMFILT: nelems_altDiffMax should be removed')
    end if
    do elem = 1, numElem
      if ( list_altDiffMax(elem) /= MPC_missingValue_INT .and. .not. utl_isEqual(value_altDiffMax(elem), MPC_missingValue_R8) ) then
        elemIndex = findElemIndex(list_altDiffMax(elem))
        if ( elemIndex >= 1 .and. elemIndex <= numElem ) then
          altDiffMax(elemIndex) = value_altDiffMax(elem)
          write(*,*) ' filt_setup: altDiffMax value for ', elemList(elemIndex), ' is set to ', altDiffMax(elemIndex)
        else
          call utl_abort('filt_setup: Error in value setting for altDiffMax')
        end if
      end if
    end do

    !
    !- Set the topographic rejection list
    !
    if (all(list_topoFilt(:) == '**')) then
      ! default list
      filtTopoList(1) = 'SF'
      filtTopoList(2) = 'UA'
      filtTopoList(3) = 'AI'
      filtTopoList(4) = 'SW'
      filtTopoList(5) = 'PR'
      filtTopoList(6) = 'AL'
      filtTopoList(7) = 'TO'
      filtTopoList(8) = 'CH'
    else
      do obsFamilyIndex = 1, nTopoFiltFam
        if (list_topoFilt(obsFamilyIndex) /= '**') then
          filtTopoList(obsFamilyIndex) = list_topoFilt(obsFamilyIndex)
        end if
      end do
    end if

    initialized = .true.

  end subroutine filt_setup

  !--------------------------------------------------------------------------
  ! filt_getSfcBufferZoneCHheight
  !--------------------------------------------------------------------------
  function filt_getSfcBufferZoneCHheight()   result(sfcBufferZoneCHheight)
    implicit none

    ! Result
    real(8) :: sfcBufferZoneCHheight  ! height buffer zone (in m) for rejecting CH family obs near sfc

    sfcBufferZoneCHheight = surfaceBufferZoneCH_Height

  end function filt_getSfcBufferZoneCHheight

  !--------------------------------------------------------------------------
  ! filt_suprep
  !--------------------------------------------------------------------------
  subroutine filt_suprep(obsSpaceData)
    !
    ! :Purpose: Select the data in the obsSpaceData which are to be assimilated
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData

    ! Locals:
    integer :: bodyIndex, headerIndex
    integer :: ipres, ivco, loopIndex
    integer :: idburp, ivnm, refValue, iknt, iknt_mpiglobal
    logical :: llok, llrej, llbogus

    call utl_tmg_start(22,'----ObsFiltSuprep')

    if(mmpi_myid == 0) write(*,*) 'starting subroutine filt_suprep'

    iknt = 0

    BODY: do bodyIndex = 1, obs_numbody( obsSpaceData )
      headerIndex = obs_bodyElem_i( obsSpaceData, OBS_HIND, bodyIndex   )
      ivnm        = obs_bodyElem_i( obsSpaceData, OBS_VNM , bodyIndex   )
      idburp      = obs_headElem_i( obsSpaceData, OBS_ITY , headerIndex )
      !
      ! Unwanted data types via types specified in NLIST
      !
      llok = .false.
      do loopIndex = 1, filt_nelems
        llok = ( ivnm == filt_nlist( loopIndex ) ) .or. llok
      end do
      if (.not.llok) then
        call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
        cycle BODY
      end if
      !
      ! Allow gz for bogus data only in analysis case
      !
      llbogus = ( idburp == 150 .or. idburp == 151 .or. idburp == 152 .or. idburp == 153 )
      if  ( (filterMode == 'analysis' .or. filterMode == 'FSO') .and. llok .and. ivnm == BUFR_NEGZ .and. .not.llbogus ) then
        if (rejectGZforAnalysis) llok=.false.
      end if
      !
      ! Ground-based GPS (GP) data (codtyp 189)
      ! LLOK = .TRUE. DY DEFAULT IF ELEMENT IS IN NLIST
      ! If gps_gb_LASSMET = .FALSE. don't want to assimilate Ps (BUFR_NEPS),
      ! Ts (BUFR_NETS), or (T-Td)s (BUFR_NESS)
      !
      if ( idburp == 189 ) then
        if (.not.gps_gb_lassmet .and. ( ivnm == BUFR_NEPS .or.  &
                                        ivnm == BUFR_NETS .or.  &
                                        ivnm == BUFR_NESS )) then
          llok = .false.
        end if
      end if
      !
      ! Exclude T-Td above level RLIMLVHU (mbs)
      !
      ivco  = obs_bodyElem_i( obsSpaceData, OBS_VCO, bodyIndex )
      ipres = nint( obs_bodyElem_r( obsSpaceData, OBS_PPP, bodyIndex))
      if ( ( ivco == obs_vcoPressure) .and. ( ivnm == BUFR_NEES ) .and.  &
           ( ipres < nint( filt_rlimlvhu *100.0d0 )) ) then
        llok=.false.
      end if
      !
      ! Bad data with quality control flags via bit list specified in NLISTFLG
      !
      llrej = .false.
      do loopIndex = 1, filt_nflags
        refValue = filt_nlistflg(loopIndex)
        llrej= flg_flagIsOnFromRefValue(obsSpaceData, bodyIndex, refValue) .or. llrej
      end do
      !
      ! SAR winds: assimilates wind speed only for SAR winds
      !
      if ( ivnm == BUFR_NEFS ) then
        if ( idburp /= 204 ) then
          llok = .false.
        end if
      end if
      if ( llok .and. .not. llrej ) then
        call obs_bodySet_i( obsSpaceData, OBS_ASS, bodyIndex, obs_assimilated )
        iknt = iknt + 1
      else
        call obs_bodySet_i( obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated )
      end if

    end do body

    call mmpi_allReduce(iknt, iknt_mpiglobal, mmpi_sum)
    if(mmpi_myid == 0) write(*,*) '  Number of data to be assimilated: ', iknt_mpiglobal

    if(mmpi_myid == 0) write(*,*) 'end of filt_suprep'

    ! abort if there is no data to be assimilated
    if (iknt_mpiglobal == 0 ) then
       call utl_abort('filt_suprep: No data to be assimilated')
    end if

    call utl_tmg_stop(22)

  end subroutine filt_suprep

  !--------------------------------------------------------------------------
  ! filt_topo
  !--------------------------------------------------------------------------
  subroutine filt_topo(columnTrlOnTrlLev, obsSpaceData, beSilent)
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)    :: columnTrlOnTrlLev
    type(struct_obs)       , intent(inout) :: obsSpaceData
    logical,                 intent(in)    :: beSilent

    if (all(filtTopoList(:) == '  ')) then

      if (mmpi_myid == 0 .and. .not.beSilent) then
        write(*,*)
        write(*,*)' --------------------------------------------------------------'
        write(*,*)' - filt_topo: NO topographic filtering                         '
        write(*,*)' --------------------------------------------------------------'
      end if

    else

      if (mmpi_myid == 0 .and. .not.beSilent) then
        write(*,*)
        write(*,*)' --------------------------------------------------------------'
        write(*,*)' - filt_topo: topographic filtering of the following obs family'
        write(*,*)' --------------------------------------------------------------'
      end if

      if (any(filtTopoList(:) == 'SF')) call filt_topoSurface   (columnTrlOnTrlLev,obsSpaceData,beSilent)
      if (any(filtTopoList(:) == 'UA')) call filt_topoRadiosonde(columnTrlOnTrlLev,obsSpaceData,beSilent)
      if (any(filtTopoList(:) == 'AI')) call filt_topoAISW      (columnTrlOnTrlLev,obsSpaceData,'AI',beSilent)
      if (any(filtTopoList(:) == 'SW')) call filt_topoAISW      (columnTrlOnTrlLev,obsSpaceData,'SW',beSilent)
      if (any(filtTopoList(:) == 'PR')) call filt_topoProfiler  (columnTrlOnTrlLev,obsSpaceData,beSilent)
      if (any(filtTopoList(:) == 'AL')) call filt_topoAladin    (columnTrlOnTrlLev,obsSpaceData,beSilent)
      if (any(filtTopoList(:) == 'TO')) call filt_topoTovs      (columnTrlOnTrlLev,obsSpaceData,beSilent)
      if (any(filtTopoList(:) == 'CH')) call filt_topoChemistry (columnTrlOnTrlLev,obsSpaceData,beSilent)

    end if

  end subroutine filt_topo

  !--------------------------------------------------------------------------
  ! filt_topoSurface
  !--------------------------------------------------------------------------
  subroutine filt_topoSurface(columnTrlOnTrlLev, obsSpaceData, beSilent)
    !
    ! :Purpose: Refuse elements which are too far away from the surface.
    !           Replace the pressure of elements which are slightly below
    !           the model surface by the pressure of the trial field.
    !
    implicit none

    ! Arguments:s
    type(struct_columnData), intent(in)    :: columnTrlOnTrlLev
    type(struct_obs)       , intent(inout) :: obsSpaceData
    logical                , intent(in)    :: beSilent

    ! Locals:
    real(8) :: altitudeDiff
    integer :: headerIndex, bodyIndex, familyIndex, elemIndex
    integer :: ivnm,countAssim
    integer :: countAcc(numElem),countRej(numElem)
    integer, parameter :: numFamily = 3
    character(len=2) :: list_family(numFamily)
    character(len=4) :: varLevel

    !
    ! reset gps_gb_dzmax for gb gps ztd to value in namelist file
    !
    altDiffMax(findElemIndex(BUFR_NEZD)) = gps_gb_dzmax

    if (  .not.beSilent ) then
      write(*,*) ' '
      write(*,*) ' filt_topoSurface: '
      write(*,*) ' '
      write(*,*) '*****************************************************'
      write(*,222) 'ELEMENTS                  ',(elemList(elemIndex),elemIndex=1,numElem)
      write(*,223) 'REJECTION BOUNDARY(METRE) ',(altDiffMax(elemIndex),elemIndex=1,numElem)
      write(*,*) '*****************************************************'
      write(*,*) ' '
223    format(2x,a29,16(2x,f10.0))
    end if

    ! Loop over the families of interest
    list_family(1) = 'SF'
    list_family(2) = 'UA'
    list_family(3) = 'GP'
    FAMILY: do familyIndex = 1, numFamily

       ! set counters to zero
       countRej(:)=0
       countAcc(:)=0

       ! loop over all header indices of each family
       call obs_set_current_header_list(obsSpaceData, &
            list_family(familyIndex))
       HEADER: do
          headerIndex = obs_getHeaderIndex(obsSpaceData)
          if (headerIndex < 0) exit HEADER

          ! loop over all body indices (still in the same family)
          call obs_set_current_body_list(obsSpaceData, headerIndex)
          BODY: do
             bodyIndex = obs_getBodyIndex(obsSpaceData)
             if (bodyIndex < 0) exit BODY

             ! skip this obs if it is not on height levels
             if (obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) /= obs_vcoHeight) cycle BODY

             ! skip this obs if already flagged to not be assimilated
             if (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_notAssimilated) cycle BODY

             ivnm   = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
             varLevel = vnl_varLevelFromVarnum(ivnm)
             altitudeDiff = abs( obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex) -  &
                  col_getHeight(columnTrlOnTrlLev,col_getNumLev(columnTrlOnTrlLev,varLevel),headerIndex,varLevel) )
             !
             ! apply filter to selected elements
             !
             elemIndex = findElemIndex(ivnm)
             if (elemIndex == -1) cycle BODY

             if (altitudeDiff <= altDiffMax(elemIndex)) then
                ! obs passes the acceptance criteria
                countAcc(elemIndex) = countAcc(elemIndex)+1
             else
                ! obs does not pass the acceptance criteria, reject it
                call flg_setFlag(obsSpaceData, bodyIndex, flg_18rejOro)
                call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
                countRej(elemIndex) = countRej(elemIndex)+1
             end if
          end do BODY
       end do HEADER

       if ( .not.beSilent ) then
         write(*,*) ' '
         write(*,*) '*****************************************************'
         write(*,*) 'FAMILY = ',list_family(familyIndex)
         write(*,222) 'ELEMENTS            ', (elemList(elemIndex),elemIndex=1,numElem)
         write(*,222) 'ACCEPTED  ',(countAcc(elemIndex),elemIndex=1,numElem)
         write(*,222) 'REJECTED  ',(countRej(elemIndex),elemIndex=1,numElem)
         write(*,*) '*****************************************************'
         write(*,*) ' '
       end if
222    format(2x,a29,16(2x,i10))

    end do FAMILY

    countAssim=0
    do bodyIndex=1,obs_numbody(obsSpaceData)
       if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) countAssim=countAssim+1
    end do
    if ( .not.beSilent ) write(*,'(1X," NUMBER OF DATA TO BE ASSIMILATED AFTER ADJUSTMENTS: ",i10)') countAssim
    if ( .not.beSilent ) write(*,* ) ' '

  end subroutine filt_topoSurface

  !--------------------------------------------------------------------------
  ! filt_topoRadiosonde
  !--------------------------------------------------------------------------
  subroutine filt_topoRadiosonde(columnTrlOnTrlLev, obsSpaceData, beSilent)
    !
    ! :Purpose: Refuse elements which are too far away from the surface of the model
    !           Refuse elements which are considered in the free atmosphere of
    !           the RAOB but fall in the surface boundary layer of the model atmosphere.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)    :: columnTrlOnTrlLev
    type(struct_obs),        intent(inout) :: obsSpaceData
    logical,                 intent(in)    :: beSilent

    ! Locals:
    integer :: headerIndex, bodyIndex, listIndex, elemIndex
    integer :: ivnm, countAssim
    integer :: itotacc(numElem), itotrej(numElem), isblrej(numElem)
    integer :: igzacc(numElem), igzrej(numElem), ibndrej(numElem)
    real(8) :: zval, obsPressure, altitudeDiff
    real(8) :: obsSfcAltitude, colSfcAltitude, colPressureBelow, colPressureAbove, zdelp
    logical :: llok
    real(8) :: geopotential(1), height(1)
    integer :: nlev_M
    real(8) :: lat

    if ( .not.beSilent ) then
      write(*,*) ' '
      write(*,*) ' filt_topoRadiosonde: '
      write(*,*) ' '
      write(*,*) '************************************************'
      write(*,222) ' ELEMENTS                  ',(elemList(elemIndex),elemIndex=1,numElem)
      write(*,223) ' REJECTION BOUNDARY(METRE) ',(altDiffMax(elemIndex),elemIndex=1,numElem)
      write(*,223) ' REJECTION SBL (PASCAL)    ',(surfaceBufferZone_Pres,elemIndex=1,numElem)
      write(*,*) '************************************************'
      write(*,*) ' '
223 format(2x,a29,16(2x,f6.0))
    end if

    ! set counters to zero
    itotrej(:)=0
    itotacc(:)=0
    isblrej(:)=0
    igzacc(:)=0
    igzrej(:)=0
    ibndrej(:)=0

    nlev_M = col_getNumLev(columnTrlOnTrlLev,'MM')

    ! loop over all header indices of the 'UA' family
    call obs_set_current_header_list(obsSpaceData, 'UA')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER

      ! HEIGHT GZ

      ! loop over all body indices (still in the 'UA' family)
      call obs_set_current_body_list(obsSpaceData, headerIndex)
      BODY: do
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if (bodyIndex < 0) exit BODY

        ! skip this obs if it is not on pressure level
        if( obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) /= obs_vcoPressure ) cycle BODY

        ! skip this obs if already flagged to not be assimilated
        if( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_notAssimilated ) cycle BODY

        ! skip this obs if it is not GZ
        ivnm=obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
        listIndex = findElemIndex(ivnm)
        llok = (ivnm == BUFR_NEGZ .and. listIndex /= -1)
        if (.not. llok ) cycle BODY

        ! convert altitude read from column to geopotential
        height(1) = col_getHeight(columnTrlOnTrlLev,0,headerIndex,'SF')
        lat = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)
        call phf_height2geopotential(height,lat,geopotential)

        zval = obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex)
        altitudeDiff = ( zval - geopotential(1) )/ec_rg
        ! obs is above surface, so it is ok, lets jump to the next obs
        if(altitudeDiff >= 0.0d0) cycle BODY

        if(altitudeDiff >= -1.0d0*altDiffMax(listIndex)) then
          ! obs is an acceptably small distance below the surface
          itotacc(listIndex) = itotacc(listIndex)+1
          igzacc(listIndex) = igzacc(listIndex)+1
        else
          ! too far below surface, reject
          call flg_setFlag(obsSpaceData, bodyIndex, flg_18rejOro)
          call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
          itotrej(listIndex) = itotrej(listIndex)+1
          igzrej(listIndex) = igzrej(listIndex)+1
        end if
      end do BODY
      !
      !   REJECT ELEMENTS OF U,V,T-TD,T BELOW THE MODEL SURFACE
      !   AND THOSE NON SURFACE ELEMENTS PRESENT IN THE SURFACE
      !   BOUNDARY LAYER OF THE RAOB OR OF THE MODEL.
      !   AT THIS POINT WE WANT TO KEEP OBSERVATIONS IN THE FREE
      !   ATMOSPHERE
      !
      !---Special case if station elevation is above model elevation
      !   we want to define colPressureAbove at a level above the station.
      !   To approximate that value, we will transform the difference
      !   between the 2 elevations into a difference in pressure using
      !   the rule of thumb (1Mb =8 metres)
      !---Even though TT(element=12001) is not assimmilated
      !   it is treated as if it were for the evaluation step.
      !   Otherwise we use observations of TT that are too far
      !   from the model topography in the verification.

      obsSfcAltitude = obs_headElem_r(obsSpaceData,OBS_ALT,headerIndex)
      colSfcAltitude = col_getHeight(columnTrlOnTrlLev,0,headerIndex,'SF')
      altitudeDiff = obsSfcAltitude - colSfcAltitude

      ! Set the body list & start at the beginning of the list
      call obs_set_current_body_list(obsSpaceData, headerIndex)
      BODY2: do
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if (bodyIndex < 0) exit BODY2

        if( obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) /= obs_vcoPressure ) cycle BODY2

        ! skip this obs if already flagged to not be assimilated
        if( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_notAssimilated ) cycle BODY2

        ivnm = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
        listIndex = findElemIndex(ivnm)
        llok = (ivnm /= BUFR_NEGZ .and. listIndex /= -1)
        if (.not. llok ) cycle BODY2 ! Proceed with the next bodyIndex

        obsPressure = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
        colPressureBelow = col_getElem(columnTrlOnTrlLev,1,headerIndex,'P0')
        colPressureAbove = colPressureBelow - surfaceBufferZone_Pres

        if (useEnkfTopoFilt) then
          ! Simpler rules used in the EnKF
          if(obsPressure >= colPressureAbove ) then
            if(abs(altitudeDiff) >= 50.0D0 .or. obsPressure >= colPressureBelow) then
              call flg_setFlag(obsSpaceData, bodyIndex, flg_18rejOro)
              call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
              itotrej(listIndex) = itotrej(listIndex) + 1
              ibndrej(listIndex) = ibndrej(listIndex) + 1
            end if
          end if
        else
          ! Original (and confusing) rules used in Var
          if (altitudeDiff > 0.0d0) then
            zdelp = altitudeDiff * 100.d0 / 8.0d0
            colPressureAbove = colPressureBelow - (zdelp + surfaceBufferZone_Pres)
          end if

          if(abs(altitudeDiff) <= altDiffMax(listIndex)) then
            !--Model surface and station altitude are very close
            !  Accept observation if obsPressure is within the domain
            !  of the trial field
            colPressureAbove = col_getPressure(columnTrlOnTrlLev,col_getNumLev(columnTrlOnTrlLev,'MM')-1,headerIndex,'MM')
          end if
          if(obsPressure > colPressureBelow ) then
            call flg_setFlag(obsSpaceData, bodyIndex, flg_18rejOro)
            call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
            itotrej(listIndex) = itotrej(listIndex) + 1
            ibndrej(listIndex) = ibndrej(listIndex) + 1
          else if(obsPressure <= colPressureBelow .and. obsPressure > colPressureAbove ) then
            call flg_setFlag(obsSpaceData, bodyIndex, flg_18rejOro)
            call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
            itotrej(listIndex) = itotrej(listIndex) + 1
            isblrej(listIndex) = isblrej(listIndex) + 1
          end if
        end if

      end do BODY2
    end do HEADER

    if ( .not.beSilent ) then
      write(*,*) ' '
      write(*,*) '***************************************'
      write(*,*) 'FAMILY = UA'
      write(*,222) ' ELEMENTS          ', (elemList(elemIndex),elemIndex=1,numElem)
      write(*,222) ' ACC GZ EXT  ',(igzacc(elemIndex),elemIndex=1,numElem)
      write(*,222) ' ACC TOTAL   ',(itotacc(elemIndex),elemIndex=1,numElem)
      write(*,*) '***************'
      write(*,222) ' REJ GZ EXT  ',(igzrej(elemIndex),elemIndex=1,numElem)
      write(*,222) ' REJ OUT BND ',(ibndrej(elemIndex),elemIndex=1,numElem)
      write(*,222) ' REJ SBL     ',(isblrej(elemIndex),elemIndex=1,numElem)
      write(*,222) ' REJ TOTAL   ',(itotrej(elemIndex),elemIndex=1,numElem)
      write(*,*) '***************************************'
      write(*,*) ' '
    end if
222 format(2x,a29,16(2x,i5))

    countAssim=0
    do bodyIndex=1,obs_numbody(obsSpaceData)
      if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) countAssim=countAssim+1
    end do
    if ( .not.beSilent ) write(*,'(1X," NUMBER OF DATA TO BE ASSIMILATED AFTER ADJUSTMENTS:",i10)') countAssim
    if ( .not.beSilent ) write(*,*) ' '

  end subroutine filt_topoRadiosonde

  !--------------------------------------------------------------------------
  ! filt_topoAISW
  !--------------------------------------------------------------------------
  subroutine filt_topoAISW(columnTrlOnTrlLev, obsSpaceData, obsFamily, beSilent)
    !
    ! :Purpose:  Refuse elements which are too close to the surface.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)    :: columnTrlOnTrlLev
    type(struct_obs),        intent(inout) :: obsSpaceData
    logical,                 intent(in)    :: beSilent
    character(len=2),        intent(in)    :: obsFamily

    ! Locals:
    integer :: headerIndex, bodyIndex, elemIndex, listIndex
    integer :: ivnm, countRej(numElem), countAssim
    real(8) :: obsPressure, pressureDiff

    if (obsFamily /= 'AI' .and. obsFamily /= 'SW') then
      call utl_abort('filt_topoAISW: only AI and SW family are handled by this routine. You ask for '//obsFamily)
    end if

    if ( .not.beSilent ) then
      write(*,*) ' '
      write(*,*) ' filt_topoAISW for obsFamily = ', obsFamily
      write(*,*) ' '
      write(*,*) '****************************************************'
      write(*,222) 'ELEMENTS                 ', (elemList(elemIndex),elemIndex=1,numElem)
      write(*,223) 'REJECTION BOUNDARY(HPA)  ', (surfaceBufferZone_Pres,elemIndex=1,numElem)
      write(*,*) '****************************************************'
      write(*,*) ' '
223 format(2x,a29,16(2x,f5.0))
    end if

    ! set counters to zero
    countRej(:)=0

    ! loop over all body indices of each family
    call obs_set_current_body_list(obsSpaceData, obsFamily)
    BODY: do
      bodyIndex = obs_getBodyIndex(obsSpaceData)
      if (bodyIndex < 0) exit BODY

      ! skip this observation if already flagged to not assimilate
      if(obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_notAssimilated) cycle BODY

      !
      ! reject data too close to the model orography, put to
      ! model orography, data which is below , but close to the surface.
      !
      obsPressure = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
      headerIndex = obs_bodyElem_i(obsSpaceData,OBS_HIND,bodyIndex)
      pressureDiff = col_getElem(columnTrlOnTrlLev,1,headerIndex,'P0') - obsPressure
      if ( pressureDiff < surfaceBufferZone_Pres ) then
        ivnm=obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
        listIndex = findElemIndex(ivnm)
        if(listIndex == -1) cycle BODY
        call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
        countRej(listIndex)=countRej(listIndex)+1
        call flg_setFlag(obsSpaceData, bodyIndex, flg_18rejOro)
      end if
    end do BODY

    if ( .not.beSilent ) then
      write(*,*) ' '
      write(*,*) '*****************************************************************'
      write(*,*) ' FAMILY = ',obsFamily
      write(*,222) 'ELEMENTS            ',(elemList(elemIndex),elemIndex=1,numElem)
      write(*,222) 'REJECTED  ',(countRej(elemIndex),elemIndex=1,numElem)
      write(*,*) '*****************************************************************'
      write(*,*) ' '
222 format(2x,a29,16(2x,i5))
    end if

    countAssim=0
    do bodyIndex=1,obs_numbody(obsSpaceData)
      if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) countAssim=countAssim+1
    end do
    if ( .not.beSilent ) write(*,'(1X," NUMBER OF DATA TO BE ASSIMILATED AFTER ADJUSTMENTS:",i10)') countAssim
    if ( .not.beSilent ) write(*,*) ' '

end subroutine filt_topoAISW

  !--------------------------------------------------------------------------
  ! filt_topoProfiler
  !--------------------------------------------------------------------------
  subroutine filt_topoProfiler(columnTrlOnTrlLev,obsSpaceData,beSilent)
    !
    ! :Purpose: Refuse elements which are too far away from the surface of the model
    !           Refuse elements which are considered in the free atmosphere of
    !           the RAOB but fall in the surface boundary layer of the model atmosphere.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)    :: columnTrlOnTrlLev
    type(struct_obs),        intent(inout) :: obsSpaceData
    logical,                 intent(in)    :: beSilent

    ! Locals:
    integer :: headerIndex, bodyIndex, listIndex, elemIndex
    integer :: ivnm, countAssim
    integer :: itotrej(numElem), isblrej(numElem), ibndrej(numElem)
    real(8) :: obsAltitude
    real(8) :: obsSfcAltitude,colSfcAltitude,colAltitudeBelow,colAltitudeAbove
    logical :: llok, list_is_empty

    if ( .not.beSilent ) then
      write(*,*) ' '
      write(*,*) ' filt_topoProfiler: '
      write(*,*) ' '
      write(*,*) '************************************************'
      write(*,222) ' ELEMENTS                  ',(elemList(elemIndex),elemIndex=1,numElem)
      write(*,223) ' REJECTION BOUNDARY(METRE) ',(altDiffMax(elemIndex),elemIndex=1,numElem)
      write(*,223) ' REJECTION SBL (METRE) ',(surfaceBufferZone_Height,elemIndex=1,numElem)
      write(*,*) '************************************************'
      write(*,*) ' '
223 format(2x,a29,16(2x,f6.0))
    end if

    ! set counters to zero
    itotrej(:)=0
    isblrej(:)=0
    ibndrej(:)=0

    ! loop over all header indices of the 'PR' family
    call obs_set_current_header_list(obsSpaceData, 'PR')
    HEADER: do
       headerIndex = obs_getHeaderIndex(obsSpaceData)
       if (headerIndex < 0) exit HEADER

       ! Set the body list & start at the beginning of the list
       call obs_set_current_body_list(obsSpaceData, headerIndex,list_is_empty)
       if (list_is_empty) cycle HEADER ! Proceed to the next HEADER

       !
       ! REJECT OBS BELOW THE MODEL SURFACE
       ! AND THOSE NON SURFACE ELEMENTS PRESENT IN THE SURFACE
       ! BOUNDARY LAYER OF THE RAOB OR OF THE MODEL.
       ! AT THIS POINT WE WANT TO KEEP OBSERVATIONS IN THE FREE
       ! ATMOSPHERE
       !
       colSfcAltitude = col_getHeight(columnTrlOnTrlLev,0,headerIndex,'SF')
       obsSfcAltitude = obs_headElem_r(obsSpaceData,OBS_ALT,headerIndex)

       ! loop over all body indices (still in the 'PR' family)
       BODY: do
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY

          ! skip this obs if already flagged to not be assimilated
          if (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_notAssimilated) cycle BODY

          ivnm = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
          listIndex = findElemIndex(ivnm)
          llok = (obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) == obs_vcoHeight  &
               .and. ivnm /= BUFR_NEGZ .and. listIndex /= -1)
          if (.not. llok ) cycle BODY ! Proceed to the next bodyIndex

          obsAltitude = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
          colAltitudeBelow = colSfcAltitude
          if (obsSfcAltitude > colSfcAltitude) then
             colAltitudeAbove = obsSfcAltitude + surfaceBufferZone_Height
          else
             colAltitudeAbove = colSfcAltitude + surfaceBufferZone_Height
          end if
          if(abs(obsSfcAltitude-colSfcAltitude) <= altDiffMax(listIndex)) then
             !----Model surface and station altitude are very close
             !    Accept observation if obsAltitude is within the domain
             !    of the trial field
             colAltitudeBelow = colSfcAltitude
             colAltitudeAbove = col_getHeight(columnTrlOnTrlLev,col_getNumLev(columnTrlOnTrlLev,'MM')-1,headerIndex,'MM')
          end if
          if(obsAltitude < colAltitudeBelow ) then
             call flg_setFlag(obsSpaceData, bodyIndex, flg_18rejOro)
             call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
             itotrej(listIndex)=itotrej(listIndex)+1
             ibndrej(listIndex)=ibndrej(listIndex)+1
          else if(obsAltitude >= colAltitudeBelow .and. obsAltitude < colAltitudeAbove ) then
             call flg_setFlag(obsSpaceData, bodyIndex, flg_18rejOro)
             call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
             itotrej(listIndex)=itotrej(listIndex)+1
             isblrej(listIndex)=isblrej(listIndex)+1
          end if
       end do BODY
    end do HEADER

    if ( .not.beSilent ) then
      write(*,*) ' '
      write(*,*) '***************************************'
      write(*,*) 'FAMILY = PR'
      write(*,222) ' ELEMENTS          ', (elemList(elemIndex),elemIndex=1,numElem)
      write(*,*) '************'
      write(*,222) ' REJ OUT BND ',(ibndrej(elemIndex),elemIndex=1,numElem)
      write(*,222) ' REJ SBL     ',(isblrej(elemIndex),elemIndex=1,numElem)
      write(*,222) ' REJ TOTAL   ',(itotrej(elemIndex),elemIndex=1,numElem)
      write(*,*) '***************************************'
      write(*,*) ' '
    end if
222 format(2x,a29,16(2x,i5))

    countAssim=0
    do bodyIndex=1,obs_numbody(obsSpaceData)
       if (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) countAssim=countAssim+1
    end do
    if ( .not.beSilent ) write(*,'(1X," NUMBER OF DATA TO BE ASSIMILATED AFTER ADJUSTMENTS:",i10)') countAssim
    if ( .not.beSilent ) write(*,*) ' '

  end subroutine filt_topoProfiler

  !--------------------------------------------------------------------------
  ! filt_topoAladin
  !--------------------------------------------------------------------------
  subroutine filt_topoAladin(columnTrlOnTrlLev, obsSpaceData, beSilent)
    !
    ! :Purpose: Refuse elements which are considered to be in the free atmosphere
    !           of the Aladin instrument but which fall in the surface boundary
    !           layer of the model atmosphere.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)    :: columnTrlOnTrlLev
    type(struct_obs),        intent(inout) :: obsSpaceData
    logical,                 intent(in)    :: beSilent

    ! Locals:
    integer :: headerIndex, bodyIndex, elemIndex
    integer :: ivnm, countAssim
    integer :: countAcc(numElem), countRej(numElem)
    real(8) :: obsAltitude      ! altitide of the observation
    real(8) :: colSfcAltitude   ! altitude of the model's lowest layer
    real(8) :: colAltitudeAbove ! top of the boundary layer
    logical :: list_is_empty

    if(.not. beSilent )then
      write(*,*) ' '
      write(*,*) ' filt_topoAladin: '
      write(*,*) ' '
      write(*,*) '************************************************'
      write(*,222) ' ELEMENTS              ',(elemList(elemIndex),elemIndex=1,numElem)
      write(*,223) ' REJECTION SBL (METRE) ',(surfaceBufferZone_Height,elemIndex=1,numElem)
      write(*,*) '************************************************'
      write(*,*) ' '
223 format(2x,a29,16(2x,f6.0))
    end if

    ! set counter to zero
    countAcc(:)=0
    countRej(:)=0

    ! loop over all header indices of the 'AL' family
    call obs_set_current_header_list(obsSpaceData, 'AL')
    HEADER: do
       headerIndex = obs_getHeaderIndex(obsSpaceData)
       if (headerIndex < 0) exit HEADER

       ! Set the body list & start at the beginning of the list
       call obs_set_current_body_list(obsSpaceData, headerIndex,list_is_empty)
       if (list_is_empty) cycle HEADER ! Proceed to the next HEADER

       !
       ! REJECT OBS IN THE SURFACE BOUNDARY LAYER OF THE MODEL.
       ! AT THIS POINT WE WANT TO KEEP OBSERVATIONS THAT ARE IN THE FREE
       ! ATMOSPHERE
       !
       colSfcAltitude = col_getHeight(columnTrlOnTrlLev,col_getNumLev(columnTrlOnTrlLev,'MM'), &
                               headerIndex,'MM')
       colAltitudeAbove = colSfcAltitude + surfaceBufferZone_Height

       ! Loop over all body indices (still in the 'AL' family)
       BODY: do
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY

          ! Skip this obs if it is not on height levels
          if (obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) /= obs_vcoHeight) cycle BODY

          ! Skip this obs if already flagged not to be assimilated
          if (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_notAssimilated) cycle BODY

          ! Skip this obs if it is not in the list of element types
          ivnm = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
          elemIndex = findElemIndex(ivnm)
          if(elemIndex == -1) cycle BODY


          !
          ! apply filter to selected elements
          !

          ! Reject this obs if it is in the boundary layer or below
          obsAltitude = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
          if(obsAltitude > colAltitudeAbove) then
             ! obs passes the acceptance criterion
             countAcc(elemIndex) = countAcc(elemIndex)+1

          else
             ! Flag rejection due to orography
             call flg_setFlag(obsSpaceData, bodyIndex, flg_18rejOro)

             ! Do not assimilate the observation
             call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
             countRej(elemIndex)=countRej(elemIndex)+1
          end if
       end do BODY
    end do HEADER

    if(.not. besilent)then
      write(*,*) ' '
      write(*,*) '***************************************'
      write(*,*) 'FAMILY = AL'
      write(*,222) ' ELEMENTS  ', (elemList(elemIndex),elemIndex=1,numElem)
      write(*,222) ' ACCEPTED  ', (countAcc(elemIndex),elemIndex=1,numElem)
      write(*,222) ' REJECTED  ', (countRej(elemIndex),elemIndex=1,numElem)
      write(*,*) '***************************************'
      write(*,*) ' '
    end if
222 format(2x,a29,16(2x,i5))

    countAssim=0
    do bodyIndex=1,obs_numbody(obsSpaceData)
       if (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) countAssim=countAssim+1
    end do

    if(.not. besilent)then
      write(*,'(1X," NUMBER OF DATA TO BE ASSIMILATED AFTER ADJUSTMENTS:",i10)') countAssim
      write(*,*) ' '
    end if

  end subroutine filt_topoAladin

  !--------------------------------------------------------------------------
  ! filt_topoTovs
  !--------------------------------------------------------------------------
  subroutine filt_topoTovs(columnTrlOnTrlLev, obsSpaceData, beSilent)
    !
    ! :Purpose:  Refuse data which are too close to the surface.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)    :: columnTrlOnTrlLev
    type(struct_obs),        intent(inout) :: obsSpaceData
    logical,                 intent(in)    :: beSilent

    ! Locals:
    integer :: headerIndex, bodyIndex
    integer :: idatyp, countAssim, countRej
    real(8), parameter :: minSfcPressure = 80000.d0

    if ( .not.beSilent ) then
      write(*,* ) ' '
      write(*,* ) ' filt_topoTovs: '
      write(*,* ) ' '
      write(*,* ) '****************************************************'
      write(*,222) 'ELEMENTS                  ', BUFR_NBT3
      write(*,223) 'MINIMUM SFC PRESSURE (PA) ', minSfcPressure
      write(*,* ) '****************************************************'
      write(*,* ) ' '
223 format(2x,a29,1(2x,f7.0))
    end if

    ! set counters to zero
    countRej=0

    ! loop over all header indices of the 'TO' family
    call obs_set_current_header_list(obsSpaceData, 'TO')
    HEADER: do
       headerIndex = obs_getHeaderIndex(obsSpaceData)
       if (headerIndex < 0) exit HEADER

       idatyp   = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
       if (idatyp /= 185) cycle HEADER ! Proceed to the next headerIndex

       ! loop over all body indices (still in the 'TO' family)
       call obs_set_current_body_list(obsSpaceData, headerIndex)
       BODY: do
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY

          if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) /= BUFR_NBT3) cycle BODY

          ! reject obs if the model surface pressure is below the minimum specified value
          if (col_getElem(columnTrlOnTrlLev,1,headerIndex,'P0') < minSfcPressure) then
             call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
             countRej=countRej+1
             call flg_setFlag(obsSpaceData, bodyIndex, flg_09rejBgck)
             call flg_setFlag(obsSpaceData, bodyIndex, flg_18rejOro)
          end if

       end do BODY
    end do HEADER

    if ( .not.beSilent ) then
      write(*,*) ' '
      write(*,*) '*****************************************************************'
      write(*,*) ' FAMILY = TO'
      write(*,222) 'ELEMENTS            ', BUFR_NBT3
      write(*,222) 'REJECTED  ',countRej
      write(*,*) '*****************************************************************'
      write(*,*) ' '
    end if
222 format(2x,a29,1(4x,i5))

    countAssim=0
    do bodyIndex=1,obs_numbody(obsSpaceData)
       if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) countAssim=countAssim+1
    end do
    if ( .not.beSilent ) write(*,'(1X," NUMBER OF DATA TO BE ASSIMILATED AFTER ADJUSTMENTS:",i10)') countAssim
    if ( .not.beSilent ) write(*,* ) ' '

  end subroutine filt_topoTovs

  !--------------------------------------------------------------------------
  ! filt_surfaceWind
  !--------------------------------------------------------------------------
  subroutine filt_surfaceWind(obsSpaceData, beSilent)
    !
    ! :Purpose: zap sfc wind components at land stations
    !
    IMPLICIT NONE

    ! Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData
    logical,          intent(in)    :: beSilent

    ! Locals:
    INTEGER, parameter :: JPINEL=2,JPIDLND=9
    INTEGER :: J,JID,JDATA
    LOGICAL :: LLPRINT
    INTEGER :: ITYP,IDBURP
    INTEGER :: ILISTEL(JPINEL), IDLND(JPIDLND)
    INTEGER :: IKOUNTREJ(JPINEL), IKOUNTT
    character(len=2), dimension(2) :: list_family
    integer :: index_family, headerIndex, bodyIndex
    character(len=obs_stnidLength) :: stnid
    real(pre_obsReal) :: obsLAT, obsLON, obsPPP

    DATA    IDLND / 12, 14, 146, 32, 35, 135, 136, 137, 138 /
    !
    if ( .not. discardlandsfcwind ) return
    !
    ILISTEL(1)=BUFR_NEUS
    ILISTEL(2)=BUFR_NEVS
    if ( .not.beSilent ) then
      WRITE(*,* ) ' '
      WRITE(*,* ) ' filt_surfaceWind:'
      WRITE(*,* ) ' '
      WRITE(*,* ) '*****************************************************'
      WRITE(*,222)'ELEMENTS REJECTED         ',(  ILISTEL(J),J=1,jpinel)
      WRITE(*,222)'LIST OF IDTYP             ',(   idlnd(J),J=1,jpidlnd)
      WRITE(*,* ) '*****************************************************'
      WRITE(*,* ) ' '
    end if
    LLPRINT = .FALSE.
    !cc      LLPRINT = .TRUE.
    !
    !     SET COUNTERS TO ZERO
    !
    DO J=1,JPINEL
       IKOUNTREJ(J)=0
    END DO

    !
    ! Loop over the families of interest
    !
    list_family(1) = 'SF'
    list_family(2) = 'UA'
    do index_family = 1,2
       if ( .not.beSilent ) WRITE(*,'(2x,A9,2x,A2)')'FAMILY = ',list_family(index_family)

       !
       ! loop over all header indices of each family
       !
       ! Set the header list
       ! (& start at the beginning of the list)
       call obs_set_current_header_list(obsSpaceData, &
            list_family(index_family))
       HEADER: do
          headerIndex = obs_getHeaderIndex(obsSpaceData)
          if (headerIndex < 0) exit HEADER

          !
          ! loop over all body indices (still in the same family)
          !
          ! Set the body list
          ! (& start at the beginning of the list)
          call obs_set_current_body_list(obsSpaceData, headerIndex)
          BODY: do
             bodyIndex = obs_getBodyIndex(obsSpaceData)
             if (bodyIndex < 0) exit BODY

             !             UNCONDITIONALLY REJECT SURFACE WINDS AT SYNOP/TEMP LAND STATIONS
             ITYP=obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
             IDBURP = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
             IF ( ITYP == BUFR_NEUS .OR. ITYP == BUFR_NEVS) THEN
                DO JID = 1, JPIDLND
                   IF(IDBURP == IDLND(JID) .AND. &
                        obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) THEN
                      call flg_setFlag(obsSpaceData, bodyIndex, flg_19rejLandSea)
                      call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
                      DO J = 1, JPINEL
                         IF(ITYP ==ILISTEL(J)) THEN
                            IKOUNTREJ(J)=IKOUNTREJ(J)+1
                         END IF
                      END DO
                      IF(LLPRINT .and. .not.beSilent ) THEN
                        stnid = obs_elem_c(obsSpaceData,'STID',headerIndex)
                        obsLAT = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)
                        obsLON = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex)
                        obsPPP = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
                        WRITE(*,225) 'Rej sfc wind lnd',headerIndex,ITYP,stnid, &
                             IDBURP, obsLAT, obsLON, obsPPP
225    FORMAT(2x,a13,2x,I6,2X,I5,1x,a9,1x,I6,1x,3(2x,f9.2))
                      END IF
                   END IF
                END DO
             END IF ! BUFR_NEUS or BUFR_NEVS
          END DO BODY
       END DO HEADER
       !
       if ( .not.beSilent ) then
         WRITE(*,* ) ' '
         WRITE(*,* ) '*****************************************************'
         WRITE(*,222 )'ELEMENTS            ', (  ILISTEL(J),J=1,JPINEL)
         WRITE(*,222)'REJECTED             ',(IKOUNTREJ(J),J=1,JPINEL)
         WRITE(*,* ) '*****************************************************'
         WRITE(*,* ) ' '
222    FORMAT(2x,a29,10(2x,i5))
       END IF
       !
    END DO ! family
    !
    IKOUNTT=0
    DO JDATA=1,obs_numbody(obsSpaceData)
       IF ( obs_bodyElem_i(obsSpaceData,OBS_ASS,JDATA) == obs_assimilated) IKOUNTT=IKOUNTT+1
    END DO
    if ( .not.beSilent ) WRITE(*, &
         '(1X," NUMBER OF DATA ASSIMILATED BY MIDAS AFTER ADJUSTMENTS: ",i10)') &
         IKOUNTT
    if ( .not.beSilent ) WRITE(*,* ) ' '

  end subroutine filt_surfaceWind

  !--------------------------------------------------------------------------
  ! filt_radvel
  !--------------------------------------------------------------------------
  subroutine filt_radvel(columnTrlOnTrlLev, obsSpaceData, beSilent)
    !
    ! :Purpose: Filter Radvel observations
    !           guarantee that altitude values are
    !           within bounds Altitude and check the horizontal distance between levels
    !           for further processing
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)    :: columnTrlOnTrlLev
    type(struct_obs)       , intent(inout) :: obsSpaceData
    logical                , intent(in)    :: beSilent

    ! Locals:
    integer :: bodyIndex, headerIndex, numLevels, bufrCode
    integer :: ierr, levelIndex
    real(8) :: obsAltitude, radarAltitude, beamElevation
    real(8) :: levelAltLow, levelAltHigh
    real(8) :: levelBracketLow, levelBracketHigh
    real(8) :: levelRangeNear, levelRangeFar
    logical, save :: firstCall = .true.

    ! Namelist variables:
    real(8), save :: maxRangeInterp ! max allowable horizontal distance between levels (in m) for radar winds

    namelist /namradvel/ maxRangeInterp

    if (.not.beSilent) then
      write(*,*)
      write(*,*) 'filt_radvel: begin'
    end if

    ! reading namelist variables
    if (firstCall) then
      ! default value
      maxRangeInterp = -1.0D0

      if ( utl_isNamelistPresent('namradvel', './flnml') ) then
        call utl_tmg_start(181,'low-level--readNML')
        read(utl_flnml, nml=namradvel, iostat=ierr)
        if ( ierr /= 0 ) call utl_abort('oop_raDvel_nl: Error reading namelist namradvel')
        if ( .not. beSilent ) write(*,nml=namradvel)
        call utl_tmg_stop(181)
      else if ( .not. beSilent ) then
        write(*,*)
        write(*,*) 'filt_radvel: namradvel is missing in the namelist. The default value will be taken.'
      end if
      firstCall = .false.
    end if
    !
    ! Loop over all header indices of the 'RA' family (Doppler Velocity)
    !
    call obs_set_current_header_list(obsSpaceData, 'RA')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if ( headerIndex < 0 ) exit HEADER
      !
      numLevels = col_getNumLev(columnTrlOnTrlLev, 'MM')
      ! Elevation beam (PPI)
      beamElevation = obs_headElem_r(obsSpaceData, OBS_RELE, headerIndex)
      ! Altitude radar
      radarAltitude = obs_headElem_r(obsSpaceData, OBS_ALT,  headerIndex)
      !
      ! Loop over all body indices of the 'RA' family (Doppler Velocity)
      !
      call obs_set_current_body_list(obsSpaceData, headerIndex)
      BODY: do
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if ( bodyIndex < 0 ) exit BODY
        ! Check that this observation has the expected bufr element ID
        bufrCode = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex)
        if ( bufrCode /= bufr_radvel ) cycle BODY
        ! only process observations flagged to be assimilated
        if ( obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) /= obs_assimilated ) cycle BODY

        ! Altitude of observation
        obsAltitude = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex)

        ! TODO we need to understand why model height at "numLevels" is not always the lowest
        ! Observations are rejected if their altitude is below the height 1 which may not be lowest model level...
        levelAltLow  = col_getHeight(columnTrlOnTrlLev, numLevels, headerIndex,'MM')
        if ( obsAltitude < levelAltLow ) then
          call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyindex, obs_notAssimilated)
          call flg_setFlag(obsSpaceData, bodyIndex, flg_11rejSelect)
          cycle BODY
        end if

        ! Observations are rejected if observation is above the height of first (highest) model level
        levelAltHigh = col_getHeight(columnTrlOnTrlLev, 1, headerIndex,'MM')
        if ( obsAltitude > levelAltHigh ) then
          call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyindex, obs_notAssimilated)
          call flg_setFlag(obsSpaceData, bodyIndex, flg_11rejSelect)
          cycle BODY
        end if

        ! Levels that bracket the observation from OBS_LYR
        !   note to self:   like in GEM, level=1 is the highest level
        levelIndex = obs_bodyElem_i(obsSpaceData, OBS_LYR, bodyIndex)
        levelBracketHigh = col_getHeight(columnTrlOnTrlLev, levelIndex,   headerIndex,'MM')
        levelBracketLow  = col_getHeight(columnTrlOnTrlLev, levelIndex+1, headerIndex,'MM')

        ! Observations are rejected if horizontal distance between levels is too large
        if ( maxRangeInterp > 0.0 ) then
          call rdv_getRangefromH(levelBracketHigh, radarAltitude, beamElevation, levelRangeFar )
          call rdv_getRangefromH(levelBracketLow,  radarAltitude, beamElevation, levelRangeNear)

          if ( abs(levelRangeFar-levelRangeNear) > maxRangeInterp ) then
            call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyindex, obs_notAssimilated)
            call flg_setFlag(obsSpaceData, bodyIndex, flg_11rejSelect)
            cycle BODY
          end if
        end if

      end do BODY
    end do HEADER
    if ( .not. beSilent ) write(*,*) 'filt_radvel: end'
  end subroutine filt_radvel

  !--------------------------------------------------------------------------
  ! filt_gpsro
  !--------------------------------------------------------------------------
  subroutine FILT_GPSRO(columnTrlOnTrlLev, obsSpaceData, beSilent)
    !
    ! :Purpose: Filter GPSRO observations
    !           Guarantee that altitude and observation values are
    !           within bounds for further processing
    !
    ! :Note: For noncompliant GPSRO observations:
    !
    !                   - Set assimilable flag to 0
    !                   - Set bit of cma flag 11 ON
    !
    use gps_mod
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)    :: columnTrlOnTrlLev
    type(struct_obs)       , intent(inout) :: obsSpaceData
    logical                , intent(in)    :: beSilent

    ! Locals:
    INTEGER :: headerIndex, IDATYP, bodyIndex, numROProfiles
    INTEGER :: JL, ISAT, IQLF, iProfile, varNum, IDSC, IBD
    REAL(8) :: ZMT, Rad, Geo, LAT, LON, AZM
    REAL(8) :: HNH1, HSF, HTP, HMIN, HMAX, ZOBS, ZREF, ZSAT
    LOGICAL :: LLEV, LOBS, LNOM, LSAT, LAZM, LALL, LDSC, LEDR
    REAL(8) :: PRad, CCoC(3)
    REAL(8) :: latrd, lonrd, dR(gps_ro_maxprfsize)

    if (.not.beSilent) then
      write(*,*)
      write(*,*) 'filt_gpsro: begin'
    end if
    !
    ! Loop over all header indices of the 'RO' family:
    !
    call obs_set_current_header_list(obsSpaceData,'RO')
    numROProfiles=0
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER
      !
      ! Process only refractivity data (codtyp 169):
      !
      IDATYP = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
      if ( IDATYP == 169 ) then
        numROProfiles=numROProfiles+1
        !
        ! Basic variables of the profile:
        !
        AZM  = obs_headElem_r(obsSpaceData,OBS_AZA ,headerIndex)
        ISAT = obs_headElem_i(obsSpaceData,OBS_SAT ,headerIndex)
        IQLF = obs_headElem_i(obsSpaceData,OBS_ROQF,headerIndex)
        Rad  = obs_headElem_r(obsSpaceData,OBS_TRAD,headerIndex)
        Geo  = obs_headElem_r(obsSpaceData,OBS_GEOI,headerIndex)
        LNOM = .NOT.BTEST(IQLF,16-1)
        LAZM = .TRUE.
        !
        ! Check if the satellite is within the accepted set:
        !
        ZSAT = ABS(gps_WGPS(ISAT,1))+ABS(gps_WGPS(ISAT,2))+ABS(gps_WGPS(ISAT,3))+ABS(gps_WGPS(ISAT,4))
        LSAT = ( ZSAT > 0.d0 )
        !
        ZMT = col_getHeight(columnTrlOnTrlLev,0,headerIndex,'SF')
        !
        ! Acceptable height limits:
        !
        JL = 1
        HTP = col_getHeight(columnTrlOnTrlLev,JL,headerIndex,'TH')
        HSF = ZMT+gps_SurfMin
        IF (HSF < gps_HsfMin) HSF=gps_HsfMin
        IF (HTP > gps_HtpMax) HTP=gps_HtpMax
        HMIN=Geo+HSF
        HMAX=Geo+HTP
        !
        ! Loop over all body indices for this headerIndex:
        ! (start at the beginning of the list)
        !
        call obs_set_current_body_list(obsSpaceData, headerIndex)
        BODY: do
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY
          !
          ! Altitude and reference order of magnitude value:
          !
          HNH1= obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
          varNum = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
          if (varNum == bufr_nebd) then
            HNH1=HNH1-Rad
            ZREF = 0.025d0*exp(-HNH1/6500.d0)
            LAZM = (-0.1d0 < AZM .AND. AZM < 360.1)
          else
            ZREF = 300.d0*exp(-HNH1/6500.d0)
            LAZM = .TRUE.
          end if
          !
          ! Observation:
          !
          ZOBS= obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex)
          !
          ! Positively verify that the altitude is within bounds:
          !
          LLEV= (HNH1 > HMIN) .AND. (HNH1 < HMAX)
          !
          ! Positively verify that the observable is within bounds:
          !
          LOBS= (ZOBS > (0.3d0*ZREF)) .AND. (ZOBS < (3.d0*ZREF))
          !
          ! Mark as not assimilable unless all conditions are satisfied:
          !
          LALL = LLEV .AND. LOBS .AND. LAZM .AND. LNOM .AND. LSAT
          if ( .NOT.LALL ) then
            call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex, obs_notAssimilated)
            call flg_setFlag(obsSpaceData, bodyIndex, flg_11rejSelect)
          end if
          ! Do not assimilate bending in mode gps_Level_RO = gps_Level_RO_Ref:
          if (varNum == bufr_nebd .and. gps_Level_RO == gps_Level_RO_Ref) then
            call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex, obs_notAssimilated)
          endif
          ! Do not assimilate refractivity in mode gps_Level_RO = gps_Level_RO_Bnd:
          if (varNum == bufr_nerf .and. gps_Level_RO == gps_Level_RO_Bnd) then
            call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex, obs_notAssimilated)
          endif
        end do BODY
      end if
    end do HEADER
    call gps_setNumROProfiles(numROProfiles)

    !
    ! List to enumerate and cross-link GPSRO headers 0 <= iProfile < gps_numROProfiles):
    !
    if (gps_numROProfiles > 0) then
      iProfile=0
      !
      ! Loop over all header indices of the 'RO' family:
      !
      call obs_set_current_header_list(obsSpaceData,'RO')
      HEADER2: do
        headerIndex = obs_getHeaderIndex(obsSpaceData)
        if (headerIndex < 0) exit HEADER2
        !
        ! Process only refractivity data (codtyp 169):
        !
        IDATYP = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
        if ( IDATYP == 169 ) then
          iProfile=iProfile+1
          LAT  = obs_headElem_r(obsSpaceData,OBS_LAT ,headerIndex)
          LON  = obs_headElem_r(obsSpaceData,OBS_LON ,headerIndex)
          AZM  = obs_headElem_r(obsSpaceData,OBS_AZA ,headerIndex)*MPC_RADIANS_PER_DEGREE_R8
          ISAT = obs_headElem_i(obsSpaceData,OBS_SAT ,headerIndex)
          IQLF = obs_headElem_i(obsSpaceData,OBS_ROQF,headerIndex)
          LDSC = .not.btest(IQLF,16-3)
          if (LDSC) then
            IDSC = 0
          else
            IDSC = 1
          end if
          varNum = -1
          call obs_set_current_body_list(obsSpaceData, headerIndex)
          !
          ! Loop over all body indices
          ! For storing varNum of each profile for the 'RO' family
          !
          LAZM = (-1.58d0 < AZM .AND. AZM < 6.29d0)
          LEDR = (LAZM .and. gps_roCurvAnisot)
          if (LEDR) then
            ! Correction for curvature anisotropy will be applied.
            ! Evaluate here the reference center of curvature and radius
            call phf_Rad_CCoC(LAT, LON, AZM, PRad, CCoC)
          end if
          dR(:) = 0.d0
          ibd = 1
          BODY2: do
            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if (bodyIndex < 0) exit BODY2
            varNum = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
            if (LEDR) then
              latrd = obs_bodyElem_r(obsSpaceData,OBS_LATD,bodyIndex)
              lonrd = obs_bodyElem_r(obsSpaceData,OBS_LOND,bodyIndex)
              if (-1.58d0 < latrd .and. latrd < 1.58d0 .and. -3.15d0 < lonrd .and. lonrd < 6.29d0) then
                ! Evaluate the offset of the ellipsoid from the reference sphere
                call phf_delR(latrd, lonrd, CCoC, PRad, dR(ibd))
              end if
            end if
            ibd = ibd+1
          end do BODY2
          ! Normal values of dR are approx -20 < dR(i) < 20, in m.
          ! To guard against gross latlon error, dR offsets are limited to 60 (m).
          ! If limit is reached, assume latlons are bad in this specific profile, and revert to
          ! no curvature anisotropy correction in it.
          if ( any(abs(dR) > 60.d0)) dR=0.d0
          ! Store profile info, including dR, in a table.
          call gps_setROIndexPrf(iProfile, headerIndex, varNum, ISAT, IDSC, dR)
          if (.not.beSilent) write(*,*)'RO Prf', gps_numROProfiles, iProfile, varNum, ISAT, IDSC, LEDR
        end if
      end do HEADER2
    end if

    if (.not.beSilent) write(*,*) 'filt_gpsro: end'
  end subroutine FILT_GPSRO

  !--------------------------------------------------------------------------
  !  filt_backScatAnisIce
  !--------------------------------------------------------------------------
  subroutine filt_backScatAnisIce(obsSpaceData, beSilent)
    !
    ! :Purpose: Filter scatterometer backscatter anisotropy observations
    !           where wind speed is too small
    !
    ! :Note: For noncompliant observations:
    !
    !                   - Set assimilable flag to 0
    !                   - Set bit of cma flag 13 ON
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData
    logical,          intent(in)    :: beSilent

    ! Locals:
    integer :: bufrCode, headerIndex, bodyIndex
    real(8) :: modelWindSpeed

    if (.not. obs_famExist(obsSpaceData,'GL')) return

    if (.not. beSilent) then
      write(*,*)
      write(*,*) ' filt_backScatAnisIce: begin'
    end if

    ! loop over all body indices
    call obs_set_current_body_list( obsSpaceData, 'GL' )

    BODY: do

      bodyIndex = obs_getBodyIndex( obsSpaceData )
      if ( bodyIndex < 0 ) exit BODY

      bufrCode = obs_bodyElem_i( obsSpaceData, OBS_VNM, bodyIndex )

      if ( bufrCode == BUFR_ICES ) then

        headerIndex = obs_bodyElem_i( obsSpaceData, OBS_HIND, bodyIndex )
        modelWindSpeed = obs_headElem_r( obsSpaceData, OBS_MWS, headerIndex )

        if ( modelWindSpeed < 4.0 ) then
          call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)
          call flg_setFlag(obsSpaceData, bodyIndex, flg_13ompLevel1)
        end if

      else

        cycle BODY

     end if

    end do BODY

    if (.not. beSilent) write(*,*) ' filt_backScatAnisIce: end'

  end subroutine  filt_backScatAnisIce

  !--------------------------------------------------------------------------
  ! filt_iceConcentration
  !--------------------------------------------------------------------------
  subroutine filt_iceConcentration(obsSpaceData, beSilent)
    !
    ! :Purpose: Filter out observations from satellites
    !           not specified in the name list
    !
    ! :Note: For noncompliant observations:
    !
    !                   - Set assimilable flag to 0
    !
    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData
    logical,          intent(in)    :: beSilent

    ! Locals:
    character(len=obs_stnidLength) :: cstnid
    integer           :: ierr
    integer           :: headerIndex, bodyIndex, codeType, platformIndex
    logical           :: inPlatformList
    logical, save     :: firstCall = .true.
    integer, parameter :: maxPlatformIce = 50

    ! Namelist variables:
    integer, save            :: nPlatformIce                    ! MUST NOT BE INCLUDED IN NAMELIST!
    character(len=12), save  :: listPlatformIce(maxPlatformIce) ! list of ice obs 'platforms' (station IDs) to assimilate

    namelist /namPlatformIce/ nPlatformIce, listPlatformIce

    if (.not. obs_famExist(obsSpaceData,'GL')) return

    if (firstCall) then
      ! set default values for namelist variables
      nPlatformIce = MPC_missingValue_INT
      listPlatformIce(:) = '1234567890ab'

      if (utl_isNamelistPresent('namPlatformIce','./flnml')) then
        call utl_tmg_start(181,'low-level--readNML')
        read (utl_flnml, nml = namPlatformIce, iostat = ierr)
        if ( ierr /= 0 ) call utl_abort('filt_iceConcentration: Error reading namelist')
        if ( mmpi_myid == 0 ) write(*,nml=namPlatformIce)
        call utl_tmg_stop(181)
        if (nPlatformIce /= MPC_missingValue_INT) then
          call utl_abort('filt_iceConcentration: check namPlatformIce namelist section: nPlatformIce should be removed')
        end if
        nPlatformIce = 0
        do platformIndex = 1, maxPlatformIce
          if (listPlatformIce(platformIndex) == '1234567890ab') exit
          nPlatformIce = nPlatformIce + 1
        end do
      else
        write(*,*)
        write(*,*) 'filt_iceConcentration: namPlatformIce is missing in the namelist. Filtering will be skipped.'
      end if
      firstCall = .false.
    end if

    if ( nPlatformIce < 1 ) return

    if ( nPlatformIce > maxPlatformIce ) then
      call utl_abort('filt_setup: too many elements for listPlatformIce')
    end if

    if (.not. beSilent) then
      write(*,*)
      write(*,*) 'filt_iceConcentration: begin'
    end if

    ! loop over all header indices of the 'GL' family
    call obs_set_current_header_list(obsSpaceData, 'GL')

    HEADER: do

      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER

      cstnid = obs_elem_c ( obsSpaceData, 'STID' , headerIndex )
      codeType = obs_headElem_i( obsSpaceData, OBS_ITY, headerIndex )

      inPlatformList = .false.

      PLATFORM: do platformIndex = 1, nPlatformIce

        if ( index(cstnid,trim(listPlatformIce(platformIndex))) > 0 .or. &
             index(codtyp_get_name(codeType),trim(listPlatformIce(platformIndex))) > 0) then

          inPlatformList = .true.
          exit PLATFORM

        end if

      end do PLATFORM

      if ( .not. inPlatformList ) then

        call obs_set_current_body_list(obsSpaceData, headerIndex)
        BODY: do
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY

          call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyIndex, obs_notAssimilated)

        end do BODY

      end if

    end do HEADER

    if (.not. beSilent) write(*,*) 'filt_iceConcentration: end'

  end subroutine filt_iceConcentration

  !--------------------------------------------------------------------------
  ! filt_topoChemistry
  !--------------------------------------------------------------------------
  subroutine filt_topoChemistry(columnTrlOnTrlLev, obsSpaceData, beSilent)
    !
    ! :Purpose: Rejects elements which are too far below the model surface
    !           or above the model top.
    !
    ! :Comments:
    !    Flagging of bit 4 in OBS_FLG done in filt_topoChemistry instead of
    !    set_scale_chm since this subroutine is called after chm_setup,
    !    allowing use of utl_open_asciifile
    !
    !    Assumes/requires obs profiles ordered from top to bottom for
    !    application of highestLvlBelowSfc
    !
    !    The vertical boundary criteria in this routine are applied only to
    !    observations as a function of altitude and pressure.
    !
    !    If the thickness of the buffer zone below the surface, i.e. values of
    !    surfaceBufferZoneCH_Pres and surfaceBufferZoneCH_Height, is set to
    !    zero, then only the obs levels at or above the model topography
    !    (and below the model top) will be accepted.
    !
    implicit none

    ! Arguments:
    type(struct_columnData), intent(in)    :: columnTrlOnTrlLev
    type(struct_obs),        intent(inout) :: obsSpaceData
    logical,                 intent(in)    :: beSilent

    ! Locals:
    integer :: headerIndex, bodyIndex, listIndex, elemIndex, listIndex_stnid
    integer :: ivnm, countAssim, jl, icount
    real(8) :: obsAltitude, obsPressure, colTopPressure, colSfcPressure
    real(8) :: colSfcAltitude, colTopAltitude,sfcAltitude
    real(8) :: previousAltitude, previousPressure
    real(8) :: stationAltitude
    logical :: list_is_empty, highestLvlBelowSfc,stationBelowSurface
    integer, parameter :: Nmax=100
    integer :: Num_stnid_chm,nobslev,Num_chm
    character(len=13) :: CstnidList_chm(Nmax)
    integer :: countAcc_stnid(Nmax),countRej_stnid(Nmax)
    integer :: countRejflg_stnid(Nmax),countRejflg(Nmax)
    integer :: countAcc(Nmax),countRej(Nmax),iConstituentList(Nmax)
    integer :: nlev_TH

    if (.not.obs_famExist(obsSpaceData,'CH')) return

    ! Set counters to zero
    countAcc_stnid(:)=0
    countRej_stnid(:)=0
    countRejflg_stnid(:)=0
    Num_stnid_chm=0

    countAcc(:)=0
    countRej(:)=0
    countRejflg(:)=0
    Num_chm=0

    ! Identify index of TH surface/bottom level
    nlev_TH = col_getNumLev(columnTrlOnTrlLev,'TH')
    ! Loop over all header indices of the 'CH' family
    call obs_set_current_header_list(obsSpaceData, 'CH')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER

      ! Set the body list & start at the beginning of the list
      call obs_set_current_body_list(obsSpaceData, headerIndex,list_is_empty)

      if (list_is_empty) cycle HEADER ! Proceed to next HEADER

      ! Identify max number of profile points in the profile (exclude BUFR_SCALE_EXPONENT elements)
      nobslev = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex)
      bodyIndex =obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
      do jl=0,obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex)-1
        if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex+jl) == BUFR_SCALE_EXPONENT) &
          nobslev = nobslev-1
      end do

      ! Identify element index of stnid list for the CH family
      call utl_get_stringId(obs_elem_c(obsSpaceData,'STID',headerIndex),&
               nobslev,CstnidList_chm,Num_stnid_chm,Nmax,listIndex_stnid)

      ! Set pressure and geopotential height vertical boundaries.
      call filt_topoChemSetBounds

      ! Initialize reference values used in the following loop
      highestLvlBelowSfc = .false. ! Where required, to identify
                                   ! if highest accepted level below surface was found.
      previousAltitude = 1.0D10 ! Set to large value in meters
      previousPressure = 0.0d0

      ! Loop over all body indices (still in the 'CH' family) to apply rejection criteria
      icount = 0
      BODY: do
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if (bodyIndex < 0) exit BODY

        ivnm = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
        if (ivnm == BUFR_SCALE_EXPONENT) cycle BODY

        !  Identify element index of constituent ID list for the CH family.
        if (icount == 0) call utl_get_Id(obs_headElem_i(obsSpaceData,OBS_CHM,headerIndex), &
                              iConstituentList,Num_chm,Nmax,listIndex)
        icount=icount+1

        ! Check for bit 4 of OBS_FLG, indicating a 'Suspicious element'
        if (flg_flagIsOn(obsSpaceData, bodyIndex, flg_04doubtful)) then
          call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
        end if

        ! Apply conditions for acceptability/rejection of the observation height

        if (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_notAssimilated) then

          ! Already rejected from input marker/flag.
          countRej(listIndex)=countRej(listIndex)+1
          countRej_stnid(listIndex_stnid)=countRej_stnid(listIndex_stnid)+1
          countRejflg(listIndex)=countRejflg(listIndex)+1
          countRejflg_stnid(listIndex_stnid)=countRejflg_stnid(listIndex_stnid)+1

        else if (obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) == obs_vcoHeight) then

          ! Check as a function of altitude.
          call filt_topoChemAltitudeCheck

        else if (obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) == obs_vcoPressure) then

          ! Check as a function of pressure.
          call filt_topoChemPressureCheck

        else
          countAcc(listIndex)=countAcc(listIndex)+1
          countAcc_stnid(listIndex_stnid)=countAcc_stnid(listIndex_stnid)+1
        end if

      end do BODY

    end do HEADER

    if (Num_stnid_chm > 0 .and. .not.beSilent) then
      write(*,*) ' '
      write(*,*) '*****************************************************************'
      write(*,*) ' filt_topoChemistry: '
      write(*,*) ' FAMILY = CH'
      write(*,221) ' REJECTION SFC BUFFER ZONE (metres)  ',surfaceBufferZoneCH_Height
      write(*,221) ' REJECTION SFC BUFFER ZONE (Pascals) ',surfaceBufferZoneCH_Pres
      write(*,*) ' '
      write(*,222) 'ELEMENTS for CH stnids',(CstnidList_chm(elemIndex),elemIndex=1,Num_stnid_chm)
      write(*,223) 'ACCEPTED for CH stnids',(countAcc_stnid(elemIndex),elemIndex=1,Num_stnid_chm)
      write(*,223) 'REJECTED for CH stnids',(countRej_stnid(elemIndex),elemIndex=1,Num_stnid_chm)
      write(*,223) 'REJECTED due to marker',(countRejflg_stnid(elemIndex),elemIndex=1,Num_stnid_chm)
      write(*,*) ' '
      write(*,224) 'ELEMENTS for CH       ',(iConstituentList(elemIndex),elemIndex=1,Num_chm)
      write(*,224) 'ACCEPTED for CH       ',(countAcc(elemIndex),elemIndex=1,Num_chm)
      write(*,224) 'REJECTED for CH       ',(countRej(elemIndex),elemIndex=1,Num_chm)
      write(*,223) 'REJECTED due to marker',(countRejflg(elemIndex),elemIndex=1,Num_stnid_chm)
      write(*,*) '*****************************************************************'
      write(*,*) ' '

      countAssim=0
      do bodyIndex=1,obs_numbody(obsSpaceData)
        if (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) countAssim=countAssim+1
      end do
      write(*,'(1X," NUMBER OF DATA TO BE ASSIMILATED AFTER ADJUSTMENTS (after filter_topoChemistry):",i10)') countAssim
      write(*,*) ' '
    end if
221 format(2x,a36,2x,f6.0)
222 format(2x,a29,100(2x,a10))
223 format(2x,a29,100(2x,i8,2x))
224 format(2x,a29,100(2x,i6))

  contains

    !--------------------------------------------------------------------------
    ! filt_topoChemSetBounds
    !--------------------------------------------------------------------------
    subroutine filt_topoChemSetBounds()
      !
      ! :Purpose: Set pressure and geopotential height vertical boundaries.
      !

      if (obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) == obs_vcoHeight) then

        ! Set  surface and lid height vertical boundaries from the trial field
        colSfcAltitude = col_getHeight(columnTrlOnTrlLev,nlev_TH,headerIndex,'TH')
        colTopAltitude = col_getHeight(columnTrlOnTrlLev,1,headerIndex,'TH')

        ! Check acceptability of surface station elevation when relevant

        ! To identify if the station elevation is below the surface where relevant.
        stationBelowSurface = .false.
        ! Identify station altitude (not provided when value = 0.0)
        stationAltitude = obs_headElem_r(obsSpaceData,OBS_ALT,headerIndex)
        ! Check if station height is far below column sfc altitude
        sfcAltitude = col_getHeight(columnTrlOnTrlLev,0,headerIndex,'SF')
        if ( .not. utl_isEqual(stationAltitude, 0.0d0) ) then
          if (stationAltitude <  sfcAltitude - surfaceBufferZoneCH_Height) then
            ! Station elevation is much lower than the surface. Provides a
            ! warning in the event the station info needs to be checked or
            ! the acceptability threshold value (surfaceBufferZoneCH_Height)
            ! would need to be revised.
            stationBelowSurface = .true.
            write(*,*) 'WARNING from filt_topoChemistry: Obs rejected as ', &
                       'station height ',int(stationAltitude) , &
                       ' m is severely below the surface at ',int(sfcAltitude), &
                       ' m for station ', obs_elem_c(obsSpaceData,'STID',headerIndex)

          else if (stationAltitude >  sfcAltitude + surfaceBufferZoneCH_Height) then
            ! Station elevation is much higher than the surface. Provides a
            ! warning in the event the station info needs to be checked.
            write(*,*) 'WARNING from filt_topoChemistry: Station height ', &
                       int(stationAltitude) , &
                       ' m is severely above the surface at ',int(sfcAltitude), &
                       ' m for station ',obs_elem_c(obsSpaceData,'STID',headerIndex)
          end if
        end if

      else if (obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) == obs_vcoPressure) then

        ! Set  surface and lid pressure vertical boundaries from the trial field
        colSfcPressure = col_getPressure(columnTrlOnTrlLev,nlev_TH,headerIndex,'TH')
        colTopPressure = col_getPressure(columnTrlOnTrlLev,1,headerIndex,'TH')
      end if

    end subroutine filt_topoChemSetBounds

    !--------------------------------------------------------------------------
    ! filt_topoChemAltitudeCheck
    !--------------------------------------------------------------------------
    subroutine filt_topoChemAltitudeCheck()
      !
      ! :Purpose: Set for rejection or acceptance of altitude relative to
      !           vertical boundaries. This includes checking for rejection
      !           or acceptance of the station elevation, when present
      !           (when /= 0.0), relative to the column surface. When the
      !           station elevation is present and accepted, also re-adjust
      !           obs altitude relative to the column surface (when within
      !           the difference threshold for checking the obs altitude.
      !

      obsAltitude = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)

      if (nobslev > 1 .and. obsAltitude > previousAltitude) then
        ! Requires obs profiles ordered from top to bottom.
        call utl_abort('filt_topoChemistry: Profile not ordered from top ' // &
                       'to bottom for ' // &
                       obs_elem_c(obsSpaceData,'STID',headerIndex))
      end if

      if (obsAltitude > colTopAltitude) then

        ! The observation altitude is beyond the accepted vertical range.
        call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
        call flg_setFlag(obsSpaceData, bodyIndex, flg_09rejBgck)
        countRej(listIndex)=countRej(listIndex)+1
        countRej_stnid(listIndex_stnid)=countRej_stnid(listIndex_stnid)+1

      else if (stationBelowSurface) then

        ! Accept only obs levels above the surface and the highest below (or at)
        ! the surface within the buffer zone. The latter is done to allow at
        ! least one accepted observation at each location, this being more
        ! relevant for cases with small datasets. Here, the obs levels are not
        ! adjusted relative to the difference between the station and surface
        ! elevations.
        if (.not.highestLvlBelowSfc) then
          if (obsAltitude < colSfcAltitude - surfaceBufferZoneCH_Height) then
            ! Reject observation level below the buffer zone
            call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
            call flg_setFlag(obsSpaceData, bodyIndex, flg_09rejBgck)
            call flg_setFlag(obsSpaceData, bodyIndex, flg_18rejOro)
            countRej(listIndex)=countRej(listIndex)+1
            countRej_stnid(listIndex_stnid)=countRej_stnid(listIndex_stnid)+1
            highestLvlBelowSfc = .true.
          else if (obsAltitude <= colSfcAltitude) then
            ! Accept topmost observation level below the surface and within
            ! the buffer zone.
            countAcc(listIndex)=countAcc(listIndex)+1
            countAcc_stnid(listIndex_stnid)=countAcc_stnid(listIndex_stnid)+1
            ! Identify as first (topmost) accepted level below the surface.
            highestLvlBelowSfc = .true.
          end if
        else
          ! Reject obs levels below the surface (except for the topmost level
          ! below the surface withing the buffer zone)
          call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
          call flg_setFlag(obsSpaceData, bodyIndex, flg_09rejBgck)
          call flg_setFlag(obsSpaceData, bodyIndex, flg_18rejOro)
          countRej(listIndex)=countRej(listIndex)+1
          countRej_stnid(listIndex_stnid)=countRej_stnid(listIndex_stnid)+1
        end if

      else

        ! Where relevant, if the station elevation is provided reset the obs
        ! altitude as if the station elevation was at the surface before
        ! proceeding with the remaining criteria application.  This is most
        ! relevant for observations on towers or profile measurements taken
        ! above any specified station elevation. This local resetting of the
        ! obs altitude via re-assigning of the station elevation is currently
        ! only done in this routine for accepting or rejecting the obs.
        ! This is also done separately in the obs operator for the local
        ! resetting of the obs levels.
        if (.not. utl_isEqual(stationAltitude, 0.0d0) .and. stationAltitude < sfcAltitude) then
          obsAltitude = obsAltitude - stationAltitude + sfcAltitude
        end if
        if (stationAltitude > 0.0d0 .and. obsAltitude < colSfcAltitude) then
          ! The station altitude was provided.
          ! The observation altitude is below the acceptable vertical range.
          call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
          call flg_setFlag(obsSpaceData, bodyIndex, flg_09rejBgck)
          call flg_setFlag(obsSpaceData, bodyIndex, flg_18rejOro)
          countRej(listIndex)=countRej(listIndex)+1
          countRej_stnid(listIndex_stnid)=countRej_stnid(listIndex_stnid)+1
        else if ( utl_isEqual(stationAltitude, 0.0d0) .and. &
                 (obsAltitude < colSfcAltitude - surfaceBufferZoneCH_Height .or. &
                 (highestLvlBelowSfc .and. obsAltitude < colSfcAltitude))) then
          ! The station altitude was not provided.
          ! The observation altitude is below the acceptable vertical range.
          ! Reject level below the topmost first accepted level within
          ! buffer zone
          call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
          call flg_setFlag(obsSpaceData, bodyIndex, flg_09rejBgck)
          call flg_setFlag(obsSpaceData, bodyIndex, flg_18rejOro)
          countRej(listIndex)=countRej(listIndex)+1
          countRej_stnid(listIndex_stnid)=countRej_stnid(listIndex_stnid)+1
        else
          ! Accept the obs above the surface (and below the lid)
          countAcc(listIndex)=countAcc(listIndex)+1
          countAcc_stnid(listIndex_stnid)=countAcc_stnid(listIndex_stnid)+1
          if ( utl_isEqual(stationAltitude, 0.0d0) .and. obsAltitude < colSfcAltitude) then
            ! Identify as first (topmost) accepted level below the surface.
            ! (not relevant when the station height is provided)
            highestLvlBelowSfc = .true.
          end if
        end if

      end if

      ! Updated to check if obs profiles are ordered from top to bottom.
      previousAltitude = obsAltitude

    end subroutine filt_topoChemAltitudeCheck

    !--------------------------------------------------------------------------
    ! filt_topoChemPressureCheck
    !--------------------------------------------------------------------------
    subroutine filt_topoChemPressureCheck()
      !
      ! :Purpose: Set for rejection or acceptance of pressure relative to
      !           vertical boundaries.
      !

      obsPressure = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
      if (nobslev > 1 .and.  obsPressure < previousPressure) then
        ! Requires obs profiles ordered from top to bottom.
        call utl_abort('filt_topoChemistry: Profile not ordered from top ' // &
                       'to bottom for ' // &
                       obs_elem_c(obsSpaceData,'STID',headerIndex))
      end if

      if (obsPressure < colTopPressure) then
        ! The observation pressure is above the accepted vertical range.
        call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
        call flg_setFlag(obsSpaceData, bodyIndex, flg_09rejBgck)
        countRej(listIndex)=countRej(listIndex)+1
        countRej_stnid(listIndex_stnid)=countRej_stnid(listIndex_stnid)+1
      else if (obsPressure > colSfcPressure + surfaceBufferZoneCH_Pres) then
        ! The observation pressure is below the accepted vertical range.
        call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
        call flg_setFlag(obsSpaceData, bodyIndex, flg_09rejBgck)
        call flg_setFlag(obsSpaceData, bodyIndex, flg_18rejOro)
        countRej(listIndex)=countRej(listIndex)+1
        countRej_stnid(listIndex_stnid)=countRej_stnid(listIndex_stnid)+1
      else if (highestLvlBelowSfc .and. obsPressure > colSfcPressure) then
        ! Reject obs level (below the first accepted obs level) within buffer zone
        call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
        call flg_setFlag(obsSpaceData, bodyIndex, flg_09rejBgck)
        call flg_setFlag(obsSpaceData, bodyIndex, flg_18rejOro)
        countRej(listIndex)=countRej(listIndex)+1
        countRej_stnid(listIndex_stnid)=countRej_stnid(listIndex_stnid)+1
      else
        ! Accept the obs above the surface and below the lid
        ! Also accept the topmost level below the surface and within the buffer zone.
        countAcc(listIndex)=countAcc(listIndex)+1
        countAcc_stnid(listIndex_stnid)=countAcc_stnid(listIndex_stnid)+1
        highestLvlBelowSfc = .true.
      end if

      ! Updated to check if obs profiles are ordered from top to bottom.
      previousPressure = obsPressure

    end subroutine filt_topoChemPressureCheck

  end subroutine filt_topoChemistry

  !--------------------------------------------------------------------------
  ! filt_bufrCodeAssimilated
  !-------------------------------------------------------------------------
  function filt_bufrCodeAssimilated(bufrCode) result(assimilated)
    !
    ! :Purpose: To test if a bufr code part of the assimilated observation list
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: bufrCode    ! The input bufr code
    ! Result:
    logical             :: assimilated ! Assimilated of not

    ! Locals:
    integer :: elemIndex

    if (.not. initialized) call filt_setup('none')

    assimilated = .false.

    do elemIndex = 1, filt_nelems
      if (filt_nlist(elemIndex) == bufrCode) then
        assimilated = .true.
        return
      end if
    end do

  end function filt_bufrCodeAssimilated

  !--------------------------------------------------------------------------
  ! filt_getBufrCodeAssimilated
  !-------------------------------------------------------------------------
  subroutine filt_getBufrCodeAssimilated(bufrCodeList)
    !
    ! :Purpose: To get the assimilated observation list
    !
    implicit none

    ! Arguments:
    integer, intent(out) :: bufrCodeList(filt_nelems) ! The list of assimilated bufr codes

    if (.not. initialized) call filt_setup('none')

    bufrCodeList(:) = filt_nlist(1:filt_nelems)

  end subroutine filt_getBufrCodeAssimilated

  !--------------------------------------------------------------------------
  ! filt_nBufrCodeAssimilated
  !-------------------------------------------------------------------------
  function filt_nBufrCodeAssimilated() result(nBufrCode)
    !
    ! :Purpose: To get the number of assimilated observations
    !
    implicit none

    ! Result:
    integer :: nBufrCode  ! The number of assimilated observations

    if (.not. initialized) call filt_setup('none')

    nBufrCode = filt_nelems

  end function filt_nBufrCodeAssimilated

end module obsFilter_mod
