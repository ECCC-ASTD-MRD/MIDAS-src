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

module obsFilter_mod
  ! MODULE obsFilter_mod (prefix='filt' category='1. High-level functionality')
  !
  ! :Purpose: Various types of filters that are applied to the observations
  !           mostly to reject them so that they will not be assimilated.
  !
  use codePrecision_mod
  use mpi_mod
  use mpivar_mod
  use earthConstants_mod
  use MathPhysConstants_mod
  use obsSpaceData_mod
  use columnData_mod
  use bufr_mod
  use tovs_nl_mod
  use gps_mod
  use utilities_mod
  use varNameList_mod
  use physicsFunctions_mod
  use codtyp_mod
  use radvel_mod
  implicit none
  save
  private

  ! public variables
  public :: filt_rlimlvhu
  ! public procedures
  public :: filt_setup, filt_topo, filt_suprep
  public :: filt_surfaceWind, filt_gpsro,  filt_backScatAnisIce, filt_iceConcentration, filt_radvel
  public :: filt_bufrCodeAssimilated, filt_getBufrCodeAssimilated, filt_nBufrCodeAssimilated

  integer :: filt_nelems, filt_nflags
  integer, target :: filt_nlist(30)
  integer :: filt_nlistflg(15)

  logical :: discardlandsfcwind

  real(8) :: filt_rlimlvhu

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

  real(8) :: surfaceBufferZone_Pres
  real(8) :: surfaceBufferZone_Height

  integer, parameter :: nTopoFiltFam = 8
  character(len=2) :: filtTopoList(nTopoFiltFam) = '  '
  logical :: useEnkfTopoFilt

  character(len=48) :: filterMode

  logical :: initialized = .false.

contains

  !--------------------------------------------------------------------------
  ! findElemIndex
  !------------------------------------------------------------------------- 
  function findElemIndex(varNum) result(listIndex)
    implicit none
    integer :: varNum, listIndex, elemIndex

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

    character(len=*), intent(in) :: filterMode_in

    integer :: nulnam, ierr, elem, jflag, ibit
    integer :: fnom, fclos
    integer :: nelems, nlist(30)
    integer :: nflags, nlistflg(15), obsFamilyIndex
    integer :: nelems_altDiffMax, list_altDiffMax(numElem), elemIndex
    
    character(len=2) :: list_topoFilt(nTopoFiltFam)

    real(8) :: value_altDiffMax(numElem)
    real(8) :: rlimlvhu

    character(len=35) :: CREASON(-8:13)
    data creason/'JACOBIAN IMPORTANT ABOVE MODEL TOP', &
         'ABS OROGRAPH-PHI                  ', &
         'MASQUE TERRE-MER                  ', &
         'OROGRAPHIE                        ', &
         'REJECTED BY QCVAR                 ', &
         'REJECTED BY BACKGROUND CHECK      ', &
         'BACKGROUND CHECK  LEVEL 3         ', &
         'BACKGROUND CHECK  LEVEL 2         ', &
         'BACKGROUND GHECK  LEVEL 1         ', &
         'RESERVED                          ', &
         'REJECTED BY SELECTION PROCESS     ', &
         'GENERATED BY OI                   ', &
         'REJECTION BY  OI                  ', &
         'ELEMENT ON BLACK LIST             ', &
         'RESERVED                          ', &
         'CORRECTED ELEMENT                 ', &
         'INTERPOLATED ELEMENT              ', &
         'DOUBTFUL ELEMENT                  ', &
         'POSSIBLY ERRONEOUS ELEMENT        ', &
         'ERRONEOUS ELEMENT                 ', &
         'ELEMENT EXCEEDS CLIMATE EXTREME   ', &
         'ELEMENT MODIFIED OR GEN BY  ADE   '/

    namelist /namfilt/nelems, nlist, nflags, nlistflg, rlimlvhu, discardlandsfcwind, &
         nelems_altDiffMax, list_altDiffMax, value_altDiffMax, surfaceBufferZone_Pres, &
         surfaceBufferZone_Height, list_topoFilt, useEnkfTopoFilt

    filterMode = filterMode_in

    ! set default values for namelist variables
    nlist(:) = 0
    nelems = 6
    nlist(1)=11003
    nlist(2)=11004
    nlist(3)=10194
    nlist(4)=12192
    nlist(5)=12062
    nlist(6)=12063

    nlistflg(:) = 0
    nflags=6
    nlistflg(1)=2
    nlistflg(2)=4
    nlistflg(3)=5
    nlistflg(4)=9
    nlistflg(5)=11
    nlistflg(6)=12

    list_altDiffMax (:) = 0
    value_altDiffMax(:) = -1.d0
    nelems_altDiffMax = 0

    list_topoFilt(:) = '**'

    rlimlvhu = 300.d0
    discardlandsfcwind = .true.

    surfaceBufferZone_Pres   = 5000.0d0 ! default value in Pascals
    surfaceBufferZone_Height =  400.0d0 ! default value in Metres

    useEnkfTopoFilt = .false.

    nulnam=0
    ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=namfilt,iostat=ierr)
    if(ierr.ne.0) call utl_abort('filt_setup: Error reading namelist! Hint: did you replaced ltopofilt by list_topoFilt?')
    if(mpi_myid == 0) write(*,nml=namfilt)
    ierr=fclos(nulnam)

    filt_rlimlvhu    = rlimlvhu
    filt_nelems      = nelems
    filt_nlist(:)    = nlist(:)
    filt_nflags      = nflags
    filt_nlistflg(:) = nlistflg(:)

    if(mpi_myid == 0) then
      write(*,'(1X,"***********************************")')
      write(*,'(1X," ELEMENTS SELECTED FOR ASSIMILATION:",/)')
      write(*,'(1X,"***********************************")')
      do elem=1,filt_nelems
        write(*,'(15X,I5)') filt_nlist(elem)
      end do
      write(*,'(1X,"***********************************")')
      write(*,*) ' REJECT ELEMENTS WITH REJECT FLAG '
      write(*,*)'           BIT :  '
      do jflag=1,filt_nflags
        ibit= filt_nlistflg(jflag)
        write(*,*) ibit,' ',creason(ibit)
      end do
      write(*,'(1X,"***********************************")')
    end if

    !
    !- Set values for altDiffMax
    !
    if ( nelems_altDiffMax > 0 ) then
      if ( nelems_altDiffMax > numElem ) then
        call utl_abort('filt_setup: You have specified too many altDiffMax elements')
      end if
      do elem = 1, nelems_altDiffMax
        elemIndex = findElemIndex(list_altDiffMax(elem))
        if ( elemIndex >= 1 .and. elemIndex <= numElem ) then
          altDiffMax(elemIndex) = value_altDiffMax(elem)
          write(*,*) ' filt_setup: altDiffMax value for ', elemList(elemIndex), ' is set to ', altDiffMax(elemIndex)
        else
          call utl_abort('filt_setup: Error in value setting for altDiffMax')
        end if
      end do
    end if

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
  ! filt_suprep
  !--------------------------------------------------------------------------
  subroutine filt_suprep(obsSpaceData)
    !
    ! :Purpose: Select the data in the obsSpaceData which are to be assimilated
    !
    implicit none
    type(struct_obs) :: obsSpaceData
    integer :: bodyIndex, headerIndex
    integer :: ipres, ivco, ierr, loopIndex
    integer :: idburp, ivnm, iflg, ibad, iknt, iknt_mpiglobal, ilansea
    logical :: llok, llrej, llbogus

    if(mpi_myid == 0) write(*,*) 'starting subroutine filt_suprep'

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
        llok=.false.
      end if
      !
      ! Ground-based GPS (GP) data (codtyp 189)
      ! LLOK = .TRUE. DY DEFAULT IF ELEMENT IS IN NLIST
      ! If LASSMET = .FALSE. don't want to assimilate Ps (BUFR_NEPS),
      ! Ts (BUFR_NETS), or (T-Td)s (BUFR_NESS)
      !
      if ( idburp == 189 ) then
        if (.not.lassmet .and. ( ivnm == BUFR_NEPS .or.  &
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
      if ( ( ivco == 2) .and. ( ivnm == BUFR_NEES ) .and.  &
           ( ipres < nint( filt_rlimlvhu *100.0d0 )) ) then
        llok=.false.
      end if
      !
      ! Bad data with quality control flags via bit list specified in NLISTFLG
      !
      iflg = obs_bodyElem_i( obsSpaceData, OBS_FLG, bodyIndex )
      llrej = .false.
      do loopIndex = 1, filt_nflags
        ibad = 13 - filt_nlistflg( loopIndex )
        llrej=( btest(iflg,ibad) ) .or. llrej
      end do
      !
      ! Filter TOVS data: check for invalid land/sea/sea-ice flag
      !
      if (ivnm == BUFR_NBT1 .or. ivnm == BUFR_NBT2 .or. ivnm == BUFR_NBT3) then
        ! Here we have to exclude the ssmis data since the burp file does not contain the ele 8021 land_sea
        if ( ( tvs_isIdBurpTovs(idburp) ) .and. & 
             ( idburp /= codtyp_get_codtyp('ssmis') ) ) then
          ilansea  = tvs_ChangedStypValue( obsSpaceData,  headerIndex )
          if (ilansea < 0 .or. ilansea > 2  ) llok = .false.
        end if
      end if
      !
      ! SAR winds: assimilates wind speed only for SAR winds
      !
      if ( ivnm == BUFR_NEFS ) then 
        if ( idburp .ne. 204 ) then 
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

    call rpn_comm_allreduce( iknt, iknt_mpiglobal, 1, "MPI_INTEGER", "MPI_SUM", "GRID", ierr )
    if(mpi_myid == 0) write(*,*) '  Number of data to be assimilated: ', iknt_mpiglobal

    if(mpi_myid == 0) write(*,*) 'end of filt_suprep'

    ! abort if there is no data to be assimilated
    if (iknt_mpiglobal == 0 ) then
       call utl_abort('SUPREP. NO DATA TO BE ASSIMILATED')
    end if

  end subroutine filt_suprep

  !--------------------------------------------------------------------------
  ! filt_topo
  !--------------------------------------------------------------------------
  subroutine filt_topo(columnTrlOnTrlLev, obsSpaceData, beSilent)
    implicit none

    type(struct_columnData) :: columnTrlOnTrlLev
    type(struct_obs)        :: obsSpaceData
    logical, intent(in)     :: beSilent

    if (all(filtTopoList(:) == '  ')) then

      if (mpi_myid == 0 .and. .not.beSilent) then
        write(*,*)
        write(*,*)' --------------------------------------------------------------'
        write(*,*)' - filt_topo: NO topographic filtering                         '
        write(*,*)' --------------------------------------------------------------'
      end if

    else

      if (mpi_myid == 0 .and. .not.beSilent) then
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

    type(struct_columnData) :: columnTrlOnTrlLev
    type(struct_obs)        :: obsSpaceData
    logical                 :: beSilent

    real(8) :: altitudeDiff
    integer :: headerIndex, bodyIndex, familyIndex, elemIndex
    integer :: ivnm,countAssim
    integer :: countAcc(numElem),countRej(numElem)
    integer, parameter :: numFamily = 3
    character(len=2) :: list_family(numFamily)
    character(len=4) :: varLevel

    !
    ! reset dzmax for gb gps ztd to value in namelist file
    !
    altDiffMax(findElemIndex(BUFR_NEZD)) = dzmax

    if (  .not.beSilent ) then
      write(*,*) ' '
      write(*,*) ' filt_topoSurface: '
      write(*,*) ' '
      write(*,*) '*****************************************************'
      write(*,222) 'ELEMENTS                  ',(elemList(elemIndex),elemIndex=1,numElem)
      write(*,223) 'REJECTION BOUNDARY(METRE) ',(altDiffMax(elemIndex),elemIndex=1,numElem)
      write(*,*) '*****************************************************'
      write(*,*) ' '
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
             if (obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex).ne.1) cycle BODY

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
                call obs_bodySet_i( obsSpaceData,OBS_FLG,bodyIndex,  &
                     ibset(obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex),18) )
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
222    format(2x,a29,16(2x,i5))
223    format(2x,a29,16(2x,f5.0))

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
   
    ! arguments:
    type(struct_columnData) :: columnTrlOnTrlLev
    type(struct_obs) :: obsSpaceData
    logical :: beSilent

    ! locals:
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
        if( obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex).ne.2 ) cycle BODY

        ! skip this obs if already flagged to not be assimilated
        if( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_notAssimilated ) cycle BODY

        ! skip this obs if it is not GZ
        ivnm=obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
        listIndex = findElemIndex(ivnm)
        llok = (ivnm == BUFR_NEGZ .and. listIndex.ne.-1)
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
          call obs_bodySet_i(obsSpaceData,OBS_FLG,bodyIndex,  &
               ibset( obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex), 18 ))
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

        if( obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex).ne.2 ) cycle BODY2

        ! skip this obs if already flagged to not be assimilated
        if( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_notAssimilated ) cycle BODY2

        ivnm = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
        listIndex = findElemIndex(ivnm)
        llok = (ivnm.ne.BUFR_NEGZ .and. listIndex.ne.-1)
        if (.not. llok ) cycle BODY2 ! Proceed with the next bodyIndex

        obsPressure = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
        colPressureBelow = col_getElem(columnTrlOnTrlLev,1,headerIndex,'P0')
        colPressureAbove = colPressureBelow - surfaceBufferZone_Pres

        if (useEnkfTopoFilt) then
          ! Simpler rules used in the EnKF
          if(obsPressure >= colPressureAbove ) then
            if(abs(altitudeDiff) >= 50.0D0 .or. obsPressure >= colPressureBelow) then
              call obs_bodySet_i(obsSpaceData,OBS_FLG,bodyIndex,  &
                   ibset(obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex), 18 ))
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
            call obs_bodySet_i(obsSpaceData,OBS_FLG,bodyIndex,  &
                 ibset(obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex), 18 ))
            call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
            itotrej(listIndex) = itotrej(listIndex) + 1
            ibndrej(listIndex) = ibndrej(listIndex) + 1
          else if(obsPressure <= colPressureBelow .and. obsPressure > colPressureAbove ) then
            call obs_bodySet_i(obsSpaceData,OBS_FLG,bodyIndex,  &
                 ibset( obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex), 18 ))
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
223 format(2x,a29,16(2x,f6.0))

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
    type(struct_columnData) :: columnTrlOnTrlLev
    type(struct_obs) :: obsSpaceData
    logical,          intent(in) :: beSilent
    character(len=2), intent(in) :: obsFamily

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
        call obs_bodySet_i(obsSpaceData,OBS_FLG,bodyIndex,  &
             ibset(obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex),18 ))
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
    end if

    countAssim=0
    do bodyIndex=1,obs_numbody(obsSpaceData)
      if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) countAssim=countAssim+1
    end do
    if ( .not.beSilent ) write(*,'(1X," NUMBER OF DATA TO BE ASSIMILATED AFTER ADJUSTMENTS:",i10)') countAssim
    if ( .not.beSilent ) write(*,*) ' '

222 format(2x,a29,16(2x,i5))
223 format(2x,a29,16(2x,f5.0))

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
    type(struct_columnData) :: columnTrlOnTrlLev
    type(struct_obs) :: obsSpaceData
    logical :: beSilent

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
          llok = (obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) == 1  &
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
             call obs_bodySet_i(obsSpaceData,OBS_FLG,bodyIndex,  &
                  ibset( obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex), 18 ))
             call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
             itotrej(listIndex)=itotrej(listIndex)+1
             ibndrej(listIndex)=ibndrej(listIndex)+1
          else if(obsAltitude >= colAltitudeBelow .and. obsAltitude < colAltitudeAbove ) then
             call obs_bodySet_i(obsSpaceData,OBS_FLG,bodyIndex,  &
                  ibset( obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex), 18 ))
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
223 format(2x,a29,16(2x,f6.0))

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
    type(struct_columnData) :: columnTrlOnTrlLev
    type(struct_obs) :: obsSpaceData

    integer :: headerIndex, bodyIndex, elemIndex
    integer :: ivnm, countAssim
    integer :: countAcc(numElem), countRej(numElem)
    real(8) :: obsAltitude      ! altitide of the observation
    real(8) :: colSfcAltitude   ! altitude of the model's lowest layer
    real(8) :: colAltitudeAbove ! top of the boundary layer
    logical :: list_is_empty, beSilent

    if(.not. beSilent )then
      write(*,*) ' '
      write(*,*) ' filt_topoAladin: '
      write(*,*) ' '
      write(*,*) '************************************************'
      write(*,222) ' ELEMENTS              ',(elemList(elemIndex),elemIndex=1,numElem)
      write(*,223) ' REJECTION SBL (METRE) ',(surfaceBufferZone_Height,elemIndex=1,numElem)
      write(*,*) '************************************************'
      write(*,*) ' '
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
          if (obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) /= 1) cycle BODY

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
             call obs_bodySet_i(obsSpaceData,OBS_FLG,bodyIndex,  &
                  ibset( obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex), 18 ))

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
223 format(2x,a29,16(2x,f6.0))

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
    type(struct_columnData) :: columnTrlOnTrlLev
    type(struct_obs) :: obsSpaceData
    logical :: beSilent

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
    end if

    ! set counters to zero
    countRej=0

    ! loop over all header indices of the 'TO' family
    call obs_set_current_header_list(obsSpaceData, 'TO')
    HEADER: do
       headerIndex = obs_getHeaderIndex(obsSpaceData)
       if (headerIndex < 0) exit HEADER

       idatyp   = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
       if (idatyp .ne. 185) cycle HEADER ! Proceed to the next headerIndex

       ! loop over all body indices (still in the 'TO' family)
       call obs_set_current_body_list(obsSpaceData, headerIndex)
       BODY: do 
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY

          if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex).ne.BUFR_NBT3) cycle BODY

          ! reject obs if the model surface pressure is below the minimum specified value
          if (col_getElem(columnTrlOnTrlLev,1,headerIndex,'P0') < minSfcPressure) then
             call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
             countRej=countRej+1
             call obs_bodySet_i(obsSpaceData,OBS_FLG,bodyIndex, &
                  ibset(obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex),9))
             call obs_bodySet_i(obsSpaceData,OBS_FLG,bodyIndex, &
                  ibset(obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex),18))
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
223 format(2x,a29,1(2x,f7.0))

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
    type(struct_obs) :: obsSpaceData
    logical :: beSilent

    INTEGER, parameter :: JPINEL=2,JPIDLND=9
    INTEGER :: J,JID,JDATA
    LOGICAL :: LLPRINT
    INTEGER :: ITYP,IDBURP
    INTEGER :: ILISTEL(JPINEL), IDLND(JPIDLND)
    INTEGER :: IKOUNTREJ(JPINEL), IKOUNTT
    character(len=2), dimension(2) :: list_family
    integer :: index_family, index_header, index_body
    character(len=12) :: stnid
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
          index_header = obs_getHeaderIndex(obsSpaceData)
          if (index_header < 0) exit HEADER

          !
          ! loop over all body indices (still in the same family)
          !
          ! Set the body list
          ! (& start at the beginning of the list)
          call obs_set_current_body_list(obsSpaceData, index_header)
          BODY: do 
             index_body = obs_getBodyIndex(obsSpaceData)
             if (index_body < 0) exit BODY

             !             UNCONDITIONALLY REJECT SURFACE WINDS AT SYNOP/TEMP LAND STATIONS
             ITYP=obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY)
             IDBURP = obs_headElem_i(obsSpaceData,OBS_ITY,INDEX_HEADER)
             IF ( ITYP == BUFR_NEUS .OR. ITYP == BUFR_NEVS) THEN
                DO JID = 1, JPIDLND
                   IF(IDBURP == IDLND(JID) .AND. &
                        obs_bodyElem_i(obsSpaceData,OBS_ASS,INDEX_BODY) == obs_assimilated) THEN
                      call obs_bodySet_i(obsSpaceData,OBS_FLG,INDEX_BODY, &
                           ibset( obs_bodyElem_i(obsSpaceData,OBS_FLG,INDEX_BODY), 19))
                      call obs_bodySet_i(obsSpaceData,OBS_ASS,INDEX_BODY,obs_notAssimilated)
                      DO J = 1, JPINEL
                         IF(ITYP ==ILISTEL(J)) THEN
                            IKOUNTREJ(J)=IKOUNTREJ(J)+1
                         END IF
                      END DO
                      IF(LLPRINT .and. .not.beSilent ) THEN
                        stnid = obs_elem_c(obsSpaceData,'STID',INDEX_HEADER)
                        obsLAT = obs_headElem_r(obsSpaceData,OBS_LAT,INDEX_HEADER)
                        obsLON = obs_headElem_r(obsSpaceData,OBS_LON,INDEX_HEADER)
                        obsPPP = obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY)
                        WRITE(*,225) 'Rej sfc wind lnd',INDEX_HEADER,ITYP,stnid, &
                             IDBURP, obsLAT, obsLON, obsPPP
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
       END IF
222    FORMAT(2x,a29,10(2x,i5))
223    FORMAT(2x,a29,10(2x,f5.0))
224    FORMAT(2x,a17,2x,I6,2X,I5,1x,a9,1x,2(2x,f9.2))
225    FORMAT(2x,a13,2x,I6,2X,I5,1x,a9,1x,I6,1x,3(2x,f9.2))
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
  ! filt_radlvel
  !--------------------------------------------------------------------------
  subroutine filt_radvel(columnTrlOnTrlLev, obsSpaceData, beSilent)
    !
    ! :Purpose: Filter Radvel observations
    !           guarantee that altitude values are
    !           within bounds Altitude and check the horizontal distance between levels
    !           for further processing
    implicit none
    
    ! arguments
    type(struct_columnData), intent(in)    :: columnTrlOnTrlLev
    type(struct_obs)       , intent(inout) :: obsSpaceData
    logical                , intent(in)    :: beSilent

    ! locals
    integer :: bodyIndex, headerIndex, numLevels, bufrCode, obsflag
    integer :: fnom, fclos, nulnam, ierr, levelIndex
    real(8) :: obsaltitude, radaraltitude, beamelevation, levelaltlow, levelAltHigh
    real(8) :: maxRangeInterp, levelRangeNear, levelRangeFar
    
    namelist /namradvel/ maxRangeInterp

    if (.not.beSilent) then
      write(*,*)
      write(*,*) 'filt_radvel: begin'
    end if
    
    ! default value 
    maxRangeInterp = -1.0D0

    ! reading namelist variables
    if ( utl_isNamelistPresent('namradvel', './flnml') ) then
      nulnam=0
      ierr=fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
      read(nulnam, nml=namradvel, iostat=ierr)
      if ( ierr /= 0 ) call utl_abort('oop_raDvel_nl: Error reading namelist namradvel')
      if ( .not. beSilent ) write(*,nml=namradvel)
      ierr = fclos(nulnam)
    else if ( .not. beSilent ) then
      write(*,*)
      write(*,*) 'filt_radvel: namradvel is missing in the namelist. The default value will be taken.'
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
      beamElevation = obs_headElem_r(obsSpaceData, OBS_RELE, headerIndex) * MPC_RADIANS_PER_DEGREE_R8
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

        ! Levels that bracket the observation from OBS_LYR
        !   note to self:   like in GEM, level=1 is the highest level
        levelIndex = obs_bodyElem_i(obsSpaceData, OBS_LYR, bodyIndex)

        levelAltHigh = col_getHeight(columnTrlOnTrlLev, levelIndex,   headerIndex,'MM')
        levelAltLow  = col_getHeight(columnTrlOnTrlLev, levelIndex+1, headerIndex,'MM')


        ! Observations are rejected if horizontal distance between levels is too large
        if ( maxRangeInterp > 0.0 ) then
          if ( levelAltLow < radarAltitude ) then 
            levelRangeNear = 0.0
          else
            call rdv_getRangefromH(levelAltLow, radarAltitude, beamElevation, levelRangeNear)
          end if
          call rdv_getRangefromH(levelAltHigh, radarAltitude, beamElevation, levelRangeFar )

          if ( abs(levelRangeFar-levelRangeNear) > maxRangeInterp ) then
            call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyindex, obs_notAssimilated)
            obsFlag = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyindex)
            call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyindex, IBSET(obsFlag, 11))
            cycle BODY
          end if
        end if

        ! Observations are rejected if observation is below the height of last model level
        levelAltLow  = col_getHeight(columnTrlOnTrlLev, numLevels, headerIndex, 'MM')
        if ( levelAltLow > obsAltitude ) then
          call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyindex, obs_notAssimilated)
          obsFlag = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyindex)
          call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyindex, IBSET(obsFlag, 11))
        end if

        ! Observations are rejected if observation is above the height of first model level
        levelAltHigh = col_getHeight(columnTrlOnTrlLev, 1, headerIndex,'MM')
        if ( levelAltHigh < obsAltitude ) then
          call obs_bodySet_i(obsSpaceData, OBS_ASS, bodyindex, obs_notAssimilated)
          obsFlag = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyindex)
          call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyindex, IBSET(obsFlag, 11))
        end if
        
      end do BODY
    end do HEADER  
    if ( .not. beSilent ) write(*,*) 'filt_radvel: end'
  end subroutine filt_radvel

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
    IMPLICIT NONE
    !
    type(struct_columnData) :: columnTrlOnTrlLev
    type(struct_obs)        :: obsSpaceData
    logical                 :: beSilent
    !
    INTEGER :: INDEX_HEADER, IDATYP, INDEX_BODY
    INTEGER :: JL, ISAT, IQLF, iProfile, IFLG, varNum, IDSC
    REAL(8) :: ZMT, Rad, Geo, AZM
    REAL(8) :: HNH1, HSF, HTP, HMIN, HMAX, ZOBS, ZREF, ZSAT
    LOGICAL :: LLEV, LOBS, LNOM, LSAT, LAZM, LALL, LDSC
    !
    if (.not.beSilent) then
      write(*,*)
      write(*,*) 'filt_gpsro: begin'
    end if
    !
    ! Loop over all header indices of the 'RO' family:
    !
    call obs_set_current_header_list(obsSpaceData,'RO')
    gps_numROProfiles=0
    HEADER: do
      index_header = obs_getHeaderIndex(obsSpaceData)
      if (index_header < 0) exit HEADER
      !
      ! Process only refractivity data (codtyp 169):
      !
      IDATYP = obs_headElem_i(obsSpaceData,OBS_ITY,INDEX_HEADER)
      if ( IDATYP == 169 ) then
        gps_numROProfiles=gps_numROProfiles+1
        !
        ! Basic variables of the profile:
        !
        AZM  = obs_headElem_r(obsSpaceData,OBS_AZA ,INDEX_HEADER)
        ISAT = obs_headElem_i(obsSpaceData,OBS_SAT ,INDEX_HEADER)
        IQLF = obs_headElem_i(obsSpaceData,OBS_ROQF,INDEX_HEADER)
        Rad  = obs_headElem_r(obsSpaceData,OBS_TRAD,INDEX_HEADER)
        Geo  = obs_headElem_r(obsSpaceData,OBS_GEOI,INDEX_HEADER)
        LNOM = .NOT.BTEST(IQLF,16-1)
        LAZM = .TRUE.
        !
        ! Check if the satellite is within the accepted set:
        !
        ZSAT = ABS(WGPS(ISAT,1))+ABS(WGPS(ISAT,2))+ABS(WGPS(ISAT,3))+ABS(WGPS(ISAT,4))
        LSAT = ( ZSAT > 0.d0 )
        !
        ZMT = col_getHeight(columnTrlOnTrlLev,0,index_header,'SF')
        !
        ! Acceptable height limits:
        !
        JL = 1
        HTP = col_getHeight(columnTrlOnTrlLev,JL,INDEX_HEADER,'TH')
        HSF = ZMT+SURFMIN
        IF (HSF < HSFMIN) HSF=HSFMIN
        IF (HTP > HTPMAX) HTP=HTPMAX
        HMIN=Geo+HSF
        HMAX=Geo+HTP
        !
        ! Loop over all body indices for this index_header:
        ! (start at the beginning of the list)
        !
        call obs_set_current_body_list(obsSpaceData, INDEX_HEADER)
        BODY: do 
          index_body = obs_getBodyIndex(obsSpaceData)
          if (index_body < 0) exit BODY
          !
          ! Altitude and reference order of magnitude value:
          !
          HNH1= obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY)
          varNum = obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY)
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
          ZOBS= obs_bodyElem_r(obsSpaceData,OBS_VAR,INDEX_BODY)
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
            call obs_bodySet_i(obsSpaceData,OBS_ASS,INDEX_BODY, obs_notAssimilated)
            IFLG = obs_bodyElem_i(obsSpaceData,OBS_FLG,INDEX_BODY)
            call obs_bodySet_i(obsSpaceData,OBS_FLG,INDEX_BODY, IBSET(IFLG,11))
          end if
          ! Do not assimilate bending in mode levelgpsro = LEVELGPSRO_REF:
          if (varNum == bufr_nebd .and. levelgpsro == LEVELGPSRO_REF) then
            call obs_bodySet_i(obsSpaceData,OBS_ASS,INDEX_BODY, obs_notAssimilated)
          endif
          ! Do not assimilate refractivity in mode levelgpsro = LEVELGPSRO_BND:
          if (varNum == bufr_nerf .and. levelgpsro == LEVELGPSRO_BND) then
            call obs_bodySet_i(obsSpaceData,OBS_ASS,INDEX_BODY, obs_notAssimilated)
          endif
        end do BODY
      end if
    end do HEADER
    !
    ! List to enumerate and cross-link GPSRO headers 0 <= iProfile < gps_numROProfiles):
    !
    if (gps_numROProfiles > 0) then
      if(.not.allocated(gps_vRO_IndexPrf)) allocate(gps_vRO_IndexPrf(gps_numROProfiles, 10))
      iProfile=0
      !
      ! Loop over all header indices of the 'RO' family:
      !
      call obs_set_current_header_list(obsSpaceData,'RO')
      HEADER2: do
        index_header = obs_getHeaderIndex(obsSpaceData)
        if (index_header < 0) exit HEADER2
        !     
        ! Process only refractivity data (codtyp 169):
        !
        IDATYP = obs_headElem_i(obsSpaceData,OBS_ITY,INDEX_HEADER)
        if ( IDATYP == 169 ) then
          iProfile=iProfile+1
          gps_vRO_IndexPrf(iProfile, 1) = INDEX_HEADER
          ISAT = obs_headElem_i(obsSpaceData,OBS_SAT ,INDEX_HEADER)
          IQLF = obs_headElem_i(obsSpaceData,OBS_ROQF,INDEX_HEADER)
          LDSC = .not.btest(IQLF,16-3)
          if (LDSC) then
            IDSC = 0
          else
            IDSC = 1
          end if
          varNum = -1
          call obs_set_current_body_list(obsSpaceData, INDEX_HEADER)
          !
          ! Loop over all body indices
          ! For storing varNum of each profile for the 'RO' family
          !
          BODY2: do 
            index_body = obs_getBodyIndex(obsSpaceData)
            if (index_body < 0) exit BODY2
            varNum = obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY)
            if (varNum > 0) exit BODY2
          end do BODY2
          gps_vRO_IndexPrf(iProfile, 2) = varNum
          gps_vRO_IndexPrf(iProfile, 3) = ISAT
          gps_vRO_IndexPrf(iProfile, 4) = IDSC
          write(*,*)'RO Prf', gps_numROProfiles, iProfile, varNum, ISAT, IDSC
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

    ! arguments
    type(struct_obs), intent(inout) :: obsSpaceData
    logical,          intent(in)    :: beSilent

    ! locals
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
          call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex, IBSET(obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex),13))
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

    ! arguments
    type(struct_obs), intent(inout) :: obsSpaceData
    logical,          intent(in)    :: beSilent

    ! locals
    character(len=12) :: cstnid
    integer           :: nulnam, ierr
    integer           :: headerIndex, bodyIndex, codeType, iplat
    integer           :: fnom, fclos
    logical           :: inPlatformList

    ! List of satellites (id_stn in SQLite files) used for sea ice concentration
    integer            :: nPlatformIce
    integer, parameter :: maxPlatformIce = 50
    character(len=12)  :: listPlatformIce(maxPlatformIce)

    namelist /namPlatformIce/ nPlatformIce, listPlatformIce

    if (.not. obs_famExist(obsSpaceData,'GL')) return

    ! set default values for namelist variables
    nPlatformIce = 0
    listPlatformIce(:) = '1234567890ab'

    if (utl_isNamelistPresent('namPlatformIce','./flnml')) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read (nulnam, nml = namPlatformIce, iostat = ierr)
      if ( ierr /= 0 ) call utl_abort('filt_iceConcentration: Error reading namelist')
      if ( mpi_myid == 0 ) write(*,nml=namPlatformIce)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'filt_iceConcentration: namPlatformIce is missing in the namelist. The default values will be taken.'
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

      PLATFORM: do iplat = 1, nPlatformIce

        if ( index(cstnid,trim(listPlatformIce(iplat))) > 0 .or. &
             index(codtyp_get_name(codeType),trim(listPlatformIce(iplat))) > 0) then

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
    !    Flagging of bit 4 in OBS_FLG done in filt_topoChemistry instead of set_scale_chm
    !    since this subroutine is called after chm_setup, allowing use of utl_open_asciifile
    !
    implicit none
    type(struct_columnData) :: columnTrlOnTrlLev
    type(struct_obs) :: obsSpaceData
    logical :: beSilent

    integer :: headerIndex, bodyIndex, listIndex, elemIndex, listIndex_stnid
    integer :: ivnm, countAssim, jl, icount
    real(8) :: obsAltitude, obsPressure, colTopPressure, colSfcPressure
    real(8) :: colAltitudeBelow, colAltitudeAbove
    logical :: list_is_empty

    integer, parameter :: Nmax=100
    integer :: Num_stnid_chm,nobslev,Num_chm
    character(len=13) :: CstnidList_chm(Nmax)
    integer :: countAcc_stnid(Nmax),countRej_stnid(Nmax)
    integer :: countRejflg_stnid(Nmax),countRejflg(Nmax)
    integer :: countAcc(Nmax),countRej(Nmax),iConstituentList(Nmax)

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

    ! Loop over all header indices of the 'CH' family
    call obs_set_current_header_list(obsSpaceData, 'CH')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER

      ! Set the body list & start at the beginning of the list
      call obs_set_current_body_list(obsSpaceData, headerIndex,list_is_empty)

      if (list_is_empty) cycle HEADER ! Proceed to next HEADER

      ! Set geopotential height and pressure boundaries.

      colAltitudeBelow = col_getHeight(columnTrlOnTrlLev,0,headerIndex,'SF')
      colAltitudeAbove = col_getHeight(columnTrlOnTrlLev,1,headerIndex,'MM')
      colSfcPressure = col_getElem(columnTrlOnTrlLev,1,headerIndex,'P0')
      colTopPressure = col_getPressure(columnTrlOnTrlLev,1,headerIndex,'MM')

      ! Identify max number of profile points in the profile (exclude BUFR_SCALE_EXPONENT elements)

      nobslev = obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex)
      bodyIndex =obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
      do jl=0,obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex)-1
          if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex+jl) == BUFR_SCALE_EXPONENT) &
             nobslev = nobslev-1
      end do

      ! Identify element index of stnid list for the CH family.

      call utl_get_stringId(obs_elem_c(obsSpaceData,'STID',headerIndex),&
               nobslev,CstnidList_chm,Num_stnid_chm,Nmax,listIndex_stnid)

      ! Loop over all body indices (still in the 'CH' family)
      icount=0
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
        if (btest(obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex),4)) &
           call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)

        if (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_notAssimilated) then

            ! Already rejected from input marker/flag.

            countRej(listIndex)=countRej(listIndex)+1
            countRej_stnid(listIndex_stnid)=countRej_stnid(listIndex_stnid)+1
            countRejflg(listIndex)=countRejflg(listIndex)+1
            countRejflg_stnid(listIndex_stnid)=countRejflg_stnid(listIndex_stnid)+1

        else if (obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) == 1) then

           ! Check as a function of altitude.

           obsAltitude = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
           !write(*,*) 'rejected zzz ',obs_elem_c(obsSpaceData,'STID',headerIndex),obsAltitude,colSfcPressure,colTopPressure
           if( obsAltitude < colAltitudeBelow .or. obsAltitude > colAltitudeAbove ) then
               call obs_bodySet_i(obsSpaceData,OBS_FLG,bodyIndex,  &
                   ibset( obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex), 18 ))
               call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
               countRej(listIndex)=countRej(listIndex)+1
               countRej_stnid(listIndex_stnid)=countRej_stnid(listIndex_stnid)+1
           else
               countAcc(listIndex)=countAcc(listIndex)+1
               countAcc_stnid(listIndex_stnid)=countAcc_stnid(listIndex_stnid)+1
           end if 

        else if (obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) == 2) then

           ! Check as a function of pressure.

           obsPressure = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)

           if ( obsPressure > colSfcPressure .or. obsPressure < colTopPressure) then
               call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
               call obs_bodySet_i(obsSpaceData,OBS_FLG,bodyIndex,  &
                 ibset(obs_bodyElem_i(obsSpaceData,OBS_FLG,bodyIndex),18 ))
               countRej(listIndex)=countRej(listIndex)+1
               countRej_stnid(listIndex_stnid)=countRej_stnid(listIndex_stnid)+1
           else
               countAcc(listIndex)=countAcc(listIndex)+1
               countAcc_stnid(listIndex_stnid)=countAcc_stnid(listIndex_stnid)+1
           end if
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
222 format(2x,a29,100(2x,a10))
223 format(2x,a29,100(2x,i8,2x))
224 format(2x,a29,100(2x,i6))

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

    ! Argument:
    integer :: bufrCodeList(filt_nelems) ! The list of assimilated bufr codes

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

    ! Argument:
    integer :: nBufrCode  ! The number of assimilated observations

    if (.not. initialized) call filt_setup('none')

    nBufrCode = filt_nelems

  end function filt_nBufrCodeAssimilated

end module obsFilter_mod
