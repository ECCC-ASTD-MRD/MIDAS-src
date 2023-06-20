
module obsOperators_mod
  ! MODULE obsOperators_mod (prefix='oop' category='5. Observation operators')
  !
  ! :Purpose: All observation operators, including nonlinear, tangent-linear
  !           and adjoint versions.
  !
  use codePrecision_mod
  use earthConstants_mod
  use mathPhysConstants_mod
  use obsSpaceData_mod
  use columnData_mod 
  use bufr_mod
  use physicsFunctions_mod
  use gps_mod
  use midasMpi_mod
  use timeCoord_mod
  use tovsNL_mod
  use utilities_mod
  use tovsLin_mod
  use verticalCoord_mod
  use varNameList_mod
  use obsOperatorsChem_mod
  use obserrors_mod
  implicit none
  save
  private

  ! public procedures
  public :: oop_ppp_nl, oop_sfc_nl, oop_zzz_nl, oop_gpsro_nl, oop_hydro_nl
  public :: oop_gpsgb_nl, oop_tovs_nl, oop_chm_nl, oop_sst_nl, oop_ice_nl, oop_raDvel_nl
  public :: oop_Htl, oop_Had, oop_vobslyrs, oop_iceScaling

  integer, external :: get_max_rss

  real(8), parameter :: temperatureLapseRate = 0.0065D0 ! K/m (i.e. 6.5 K/km)

  ! Jacobian caches
  real*8 , allocatable :: oop_vRO_Jacobian(:,:,:)
  logical, allocatable :: oop_vRO_lJac(:)
  real*8 , allocatable :: oop_vZTD_Jacobian(:,:)

contains

  !--------------------------------------------------------------------------
  ! oop_vobslyrs
  !--------------------------------------------------------------------------
  subroutine oop_vobslyrs( columnTrl, obsSpaceData, beSilent )
    ! :Purpose:
    !      Find which model levels to use for the vertical interpolation
    !      of model fields to obs data.
    implicit none
    type(struct_columnData) :: columnTrl
    type(struct_obs) :: obsSpaceData
    logical beSilent

    integer :: levIndex,JDATA,NLEV
    real(8) :: ZLEV,ZPT,ZPB
    integer :: IOBS,layerIndex,bufrCode
    character(len=4) :: varLevel
    integer :: bodyIndex

    if (.not.beSilent) write(*,*) "Entering subroutine OOP_VOBSLYRS"

    ! 2D mode patch
    if ( col_getNumLev(columnTrl,'MM') <= 1 ) then 
      do bodyIndex = 1, obs_numbody( obsSpaceData )
        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated .and. &
             obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) == 1 ) then
          call obs_bodySet_i(obsSpaceData,OBS_LYR,bodyIndex, 0) ! set OBS_LYR = 0
        end if
      end do
      return
    end if

    !
    !-----------------------------------------------------------------------
    !         --------
    !           ETA
    !         --------
    !
    !     1. Find where extrapolation is needed
    !        ----------------------------------
    !
    !     1.1 PPP Vertical coordinate
    ! 

!$OMP PARALLEL DO PRIVATE(jdata,zlev,iobs,bufrCode,varLevel,zpt,zpb)
    do JDATA= 1,obs_numbody(obsSpaceData)
       if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,JDATA) == obs_assimilated .and. &
            obs_bodyElem_i(obsSpaceData,OBS_VCO,JDATA) == 2 ) THEN
          if (obs_bodyElem_i(obsSpaceData,OBS_VNM,JDATA) /= bufr_nedz ) THEN
             ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,JDATA)
          else
             call utl_abort('oop_vobslyr: ZLEV cannot be set, bufr_nedz not supported!')
          end if
          IOBS = obs_bodyElem_i(obsSpaceData,OBS_HIND,JDATA)
          bufrCode = obs_bodyElem_i(obsSpaceData,OBS_VNM,JDATA)
          if (bufr_IsAtmosConstituent(bufrCode)) then
             varLevel = vnl_varLevelFromVarnum(bufrCode, &
                        obs_headElem_i(obsSpaceData,OBS_CHM,IOBS))
          else
             varLevel = vnl_varLevelFromVarnum(bufrCode)
          end if
          ZPT= col_getPressure(columnTrl,1,IOBS,varLevel)
          ZPB= col_getPressure(columnTrl,COL_GETNUMLEV(columnTrl,varLevel),IOBS,varLevel)
          if ( ZLEV < ZPT ) THEN
             call obs_bodySet_i(obsSpaceData,OBS_XTR,JDATA,1)
             !
             !- !!! WARNING !!! This obs is above the model lid. 
             !  We must turn off its assimilation flag  because the
             !  current obs operators cannot deal with this situation (JFC)                  
             if (varLevel /= 'SF') then
                write(*,*) 'oop_vobslyrs: Rejecting OBS above model lid, pressure = ', ZLEV,' < ',ZPT
                call obs_bodySet_i(obsSpaceData,OBS_ASS,JDATA, obs_notAssimilated)
             end if
          else if ( ZLEV > ZPB ) THEN
             call obs_bodySet_i(obsSpaceData,OBS_XTR,JDATA,2)
          else
             call obs_bodySet_i(obsSpaceData,OBS_XTR,JDATA,0)
          end if
       end if
    end do
!$OMP END PARALLEL DO
    !
    !     1.2 ZZZ Vertical coordinate
    !
!$OMP PARALLEL do PRIVATE(jdata,zlev,iobs,bufrCode,varLevel,zpt,zpb,nlev)
    do JDATA= 1,obs_numbody(obsSpaceData)
      if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,JDATA) == obs_assimilated .and. &
           obs_bodyElem_i(obsSpaceData,OBS_VCO,JDATA) == 1 ) then
        IOBS = obs_bodyElem_i(obsSpaceData,OBS_HIND,JDATA)
        bufrCode = obs_bodyElem_i(obsSpaceData,OBS_VNM,JDATA)
        if ( bufrCode /= bufr_nedz ) then
          ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,JDATA)
          if ( bufrCode == bufr_nebd ) THEN
             ZLEV = ZLEV - obs_headElem_r(obsSpaceData,OBS_TRAD,IOBS)
          endif
        else
          call utl_abort('oop_vobslyr: ZLEV cannot be set, bufr_nedz not supported!')
        end if
        if (bufr_IsAtmosConstituent(bufrCode)) then
          varLevel = vnl_varLevelFromVarnum(bufrCode, &
                     obs_headElem_i(obsSpaceData,OBS_CHM,IOBS))
        else
          varLevel = vnl_varLevelFromVarnum(bufrCode)
        end if
        if (varLevel == 'SF') then
          ZPT= col_getHeight(columnTrl,1,IOBS,'TH')
          ZPB= col_getHeight(columnTrl,0,IOBS,'SF')
        else
          nlev=col_getNumLev(columnTrl,varLevel)
          ZPT= col_getHeight(columnTrl,1,IOBS,varLevel)
          ZPB= col_getHeight(columnTrl,NLEV,IOBS,varLevel)
        end if
        if ( ZLEV > ZPT ) then
          call obs_bodySet_i(obsSpaceData,OBS_XTR,JDATA,1)
          write(*,*) 'oop_vobslyrs: Rejecting OBS above model lid, height =', ZLEV,' > ',ZPT
          call obs_bodySet_i(obsSpaceData,OBS_ASS,JDATA, obs_notAssimilated)
        else if ( ZLEV < ZPB ) then
          call obs_bodySet_i(obsSpaceData,OBS_XTR,JDATA,2)
        else
          call obs_bodySet_i(obsSpaceData,OBS_XTR,JDATA,0)
        end if
      end if
    end do
!$OMP END PARALLEL DO
    !
    !
    !     2. FInd interpolation layer
    !        ------------------------
    !        (Model levels are assumed to be in increasing order in Mbs)
    !
    !     2.1  PPP Vertical coordinate
    !
!$OMP PARALLEL DO PRIVATE(jdata,iobs,zlev,bufrCode,varLevel,layerIndex,nlev,levIndex,zpt,zpb)
    do JDATA = 1, obs_numbody(obsSpaceData)
      call obs_bodySet_i(obsSpaceData,OBS_LYR,JDATA,0)
      if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,JDATA) == obs_assimilated .and. &
           obs_bodyElem_i(obsSpaceData,OBS_VCO,JDATA) == 2 ) then
        IOBS = obs_bodyElem_i(obsSpaceData,OBS_HIND,JDATA)
        ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,JDATA)
        bufrCode = obs_bodyElem_i(obsSpaceData,OBS_VNM,JDATA)
        if (bufr_IsAtmosConstituent(bufrCode)) then
           varLevel = vnl_varLevelFromVarnum(bufrCode, &
                      obs_headElem_i(obsSpaceData,OBS_CHM,IOBS))
        else
           varLevel = vnl_varLevelFromVarnum(bufrCode)
        end if
        layerIndex = 1
        nlev=COL_GETNUMLEV(columnTrl,varLevel)
        do levIndex = 2,NLEV - 1
          ZPT = col_getPressure(columnTrl,levIndex,IOBS,varLevel)
          if ( ZLEV > ZPT ) layerIndex = levIndex
        end do
        ZPT = col_getPressure(columnTrl,layerIndex,IOBS,varLevel)
        ZPB = col_getPressure(columnTrl,layerIndex+1,IOBS,varLevel) 
        call obs_bodySet_i(obsSpaceData,OBS_LYR,JDATA, layerIndex)
      end if
    end do
!$OMP END PARALLEL DO
    !
    !     2.2  ZZZ Vertical coordinate and surface observations
    !
!$OMP PARALLEL DO PRIVATE(jdata,iobs,zlev,bufrCode,varLevel,layerIndex,nlev,levIndex,zpt)
    do JDATA = 1, obs_numbody(obsSpaceData)
      if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,JDATA) == obs_assimilated .and. &
           obs_bodyElem_i(obsSpaceData,OBS_VCO,JDATA) == 1 ) then
        IOBS = obs_bodyElem_i(obsSpaceData,OBS_HIND,JDATA)
        ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,JDATA)
        bufrCode = obs_bodyElem_i(obsSpaceData,OBS_VNM,JDATA)
        if (bufr_IsAtmosConstituent(bufrCode)) then
          varLevel = vnl_varLevelFromVarnum(bufrCode, &
                     obs_headElem_i(obsSpaceData,OBS_CHM,IOBS))
        else
          varLevel = vnl_varLevelFromVarnum(bufrCode)
        end if
        layerIndex = 1
        nlev=COL_GETNUMLEV(columnTrl,varLevel)
        do levIndex = 2, NLEV - 1
          ZPT = col_getHeight(columnTrl,levIndex,IOBS,varLevel)
          if ( ZLEV < ZPT ) layerIndex = levIndex
        end do
        if ( bufrCode == bufr_neps .or. bufrCode == bufr_nepn .or. &
             bufrCode == bufr_nezd .or. bufrCode == bufr_gust .or. &
             bufrCode == bufr_radarPrecip .or. bufrCode == bufr_logRadarPrecip ) THEN
          ! for surface observations associated with surface analysis variables
          layerIndex = 0
        else if ( bufrCode == bufr_nets .or. bufrCode == bufr_ness .or. &
                  bufrCode == bufr_neus .or. bufrCode == bufr_nevs .or. &
                  bufrCode == bufr_nehs .or. bufrCode == bufr_vis  .or.  &
                  bufrCode == bufr_logVis ) then
          ! for surface observations associated with NON-surface analysis variables
          layerIndex = nlev - 1
        end if
        call obs_bodySet_i(obsSpaceData,OBS_LYR,JDATA, layerIndex)
      end if
    end do
!$OMP END PARALLEL DO
    !
  end subroutine oop_vobslyrs

  !--------------------------------------------------------------------------
  ! oop_ppp_nl
  !--------------------------------------------------------------------------
  subroutine oop_ppp_nl( columnTrlOnTrlLev, obsSpaceData, beSilent, cdfam, destObsColumn )
    ! :Purpose: Computation of y - H(x)
    !           for pressure-level observations.
    !           Interpolate vertically columnTrlOnTrlLev to
    !           the pressure levels of the observations.
    !           A linear interpolation in ln(p) is performed.
    !
    implicit none

    ! arguments
    type(struct_columnData) :: columnTrlOnTrlLev
    type(struct_obs)        :: obsSpaceData
    logical                 :: beSilent
    character(len=*)        :: cdfam ! family of obsservation
    integer                 :: destObsColumn

    ! locals
    integer :: headerIndex,bodyIndex,ilyr
    integer :: iass,ixtr,ivco,bufrCode,nlev_T
    real(8) :: zvar,zwb,zwt,zexp
    real(8) :: zlev,zpt,zpb,zomp,ztvg
    real(8) :: trlValueBot,trlValueTop,lat
    character(len=4) :: varName
    character(len=4) :: varLevel
    real(8),pointer :: col_ptr(:),col_ptr_tt(:),col_ptr_hu(:)
    real(8), allocatable :: geopotential(:)
    real(8) :: heightSfc(1), geopotentialSfc(1)

    if (.not.beSilent) write(*,*) 'Entering subroutine oop_ppp_nl'

    zexp = MPC_RGAS_DRY_AIR_R8 * temperatureLapseRate / ec_rg

    nlev_T = col_getNumLev(columnTrlOnTrlLev,'TH')
    allocate(geopotential(nlev_T))

    call obs_set_current_body_list(obsSpaceData, cdfam)
    BODY: do
       bodyIndex = obs_getBodyIndex(obsSpaceData)
       if (bodyIndex < 0) exit BODY

       ! Only process pressure level observations flagged to be assimilated
       iass=obs_bodyElem_i (obsSpaceData,OBS_ASS,bodyIndex)
       ivco=obs_bodyElem_i (obsSpaceData,OBS_VCO,bodyIndex)
       if (iass /= 1 .or. ivco /= 2) cycle BODY

       ixtr=obs_bodyElem_i (obsSpaceData,OBS_XTR,bodyIndex)
       bufrCode=obs_bodyElem_i (obsSpaceData,OBS_VNM,bodyIndex)
       zvar=obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex)
       zlev=obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
       headerIndex=obs_bodyElem_i (obsSpaceData,OBS_HIND,bodyIndex)

       if ( ixtr == 0 ) then

         ! Process all data within the domain of the model
         ilyr  =obs_bodyElem_i (obsSpaceData,OBS_LYR,bodyIndex)
         varName = vnl_varNameFromVarnum(bufrCode)
         varLevel = vnl_varLevelFromVarnum(bufrCode)
         zpt= col_getPressure(columnTrlOnTrlLev,ilyr  ,headerIndex,varLevel)
         zpb= col_getPressure(columnTrlOnTrlLev,ilyr+1,headerIndex,varLevel)
         zwb  = log(zlev/zpt)/log(zpb/zpt)
         zwt  = 1.d0 - zwb
         lat = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)
         if (bufrCode == bufr_nees) then
           col_ptr_hu=>col_getColumn(columnTrlOnTrlLev,headerIndex,'HU')
           col_ptr_tt=>col_getColumn(columnTrlOnTrlLev,headerIndex,'TT')
           trlValueBot=hutoes(col_ptr_hu(ilyr+1),col_ptr_tt(ilyr+1),zpb)
           trlValueTop=hutoes(col_ptr_hu(ilyr  ),col_ptr_tt(ilyr  ),zpt)
         else
           if (trim(varName) == 'Z_T') then
             col_ptr=>col_getColumn(columnTrlOnTrlLev,headerIndex,'Z_T')
             call phf_height2geopotential(col_ptr,lat,geopotential)
             trlValueBot=geopotential(ilyr+1)
             trlValueTop=geopotential(ilyr  )
           else
             col_ptr=>col_getColumn(columnTrlOnTrlLev,headerIndex,varName)
             trlValueBot=col_ptr(ilyr+1)
             trlValueTop=col_ptr(ilyr  )
           end if
         end if
         zomp = zvar-(zwb*trlValueBot+zwt*trlValueTop)
         call obs_bodySet_r(obsSpaceData,destObsColumn,bodyIndex,zomp)

       else if (ixtr == 2) then

         ! Process only GZ that is data below model's orography
         if (bufrCode == bufr_negz ) then
           ! Forward nonlinear model for geopotential data below model's orography
           ztvg = (1.0d0 + MPC_DELTA_R8 * col_getElem(columnTrlOnTrlLev,col_getNumLev(columnTrlOnTrlLev,'TH'),headerIndex,'HU'))*  &
                col_getElem(columnTrlOnTrlLev,col_getNumLev(columnTrlOnTrlLev,'TH'),headerIndex,'TT')

           ! convert height of surface to geopotential
           heightSfc(1) = col_getHeight(columnTrlOnTrlLev,0,headerIndex,'SF')
           lat = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)
           call phf_height2geopotential(heightSfc,lat,geopotentialSfc)

           zomp = (  zvar - geopotentialSfc(1) -  &
                ztvg/(temperatureLapseRate/ec_rg)*(1.D0-(zlev/col_getElem(columnTrlOnTrlLev,1,headerIndex,'P0'))**zexp))
           call obs_bodySet_r(obsSpaceData,destObsColumn,bodyIndex,zomp)
         end if

       end if

    end do body

    deallocate(geopotential)

  end subroutine oop_ppp_nl

  !--------------------------------------------------------------------------
  ! oop_zzz_nl
  !--------------------------------------------------------------------------
  subroutine oop_zzz_nl( columnTrlOnTrlLev, obsSpaceData, beSilent, cdfam,  &
                         destObsColumn )
    ! :Purpose: Computation of y - H(x) for geometric-height observations
    !           Interpolate vertically columnTrlOnTrlLev to the geometric heights (in
    !           meters) of the observations.
    !           A linear interpolation in z is performed.
    !
    ! :Notes:
    !     As a first approximation, use the geopotential height.  Once this is
    !     working, this should be changed for a calculation of the geometric
    !     height.
    !
    !     Note that, in the case of an aladin HLOS wind, the correction to zvar
    !     (OBS_VAR) is not written back to obsSpaceData.  It is simply used to
    !     calculate OMP (which is written to obsSpaceData) and then is discarded.
    !     Thereafter, if one calculates OMP - O (this will be the uncorrected O),
    !     the result will be a corrected P.
    implicit none

    ! Arguments
    type(struct_columnData), intent(in)    :: columnTrlOnTrlLev
    type(struct_obs),        intent(inout) :: obsSpaceData
    logical,                 intent(in)    :: beSilent
    integer,                 intent(in)    :: destObsColumn
    character(len=*),        intent(in)    :: cdfam ! family of observation

    ! Locals
    integer :: headerIndex,bodyIndex,ilyr,bufrCode,levIndexTop,levIndexBot
    integer :: bodyIndexStart,bodyIndexEnd,bodyIndex2
    integer :: found  ! a group of bit flags
    integer :: ierr, nulnam, fnom,fclos
    real(8) :: zvar,zwb,zwt
    real(8) :: zlev,zpt,zpb,zomp
    real(8) :: trlValueBot,trlValueTop
    character(len=4) :: varLevel
    real(8) :: value   ! temporary holder
    real(8) :: azimuth ! HLOS wind direction CW from true north
    real(8) :: tempRef ! reference temperature used to calculate HLOS wind
    real(8) :: presRef ! reference pressure used to calculate HLOS wind
    real(8) :: dwdp    ! derivative of HLOS wind wrt P
    real(8) :: dwdt    ! derivative of HLOS wind wrt T
    real(8) :: uuLyr, vvLyr   ! wind on layer, OBS_LYR
    real(8) :: uuLyr1,vvLyr1  ! wind on layer plus 1
    real(8) :: ttLyr, ppLyr   ! T, P on layer, OBS_LYR
    real(8) :: ttLyr1,ppLyr1  ! T, P on layer, OBS_LYR plus 1
    real(8) :: ttbg,  ppbg    ! background T, P at the observation location
    logical :: list_is_empty

    ! namelist variables
    logical :: do_adjust_aladin ! choose to adjust obs value as if it was retrieved using our temp and pressure 

    namelist /NAMALADIN_OBS/do_adjust_aladin

    if (.not.beSilent) write(*,*) "Entering subroutine oop_zzz_nl"

    call obs_set_current_body_list(obsSpaceData, cdfam, list_is_empty)

    if (list_is_empty)then
      return
    end if

    ! Read in the namelist NAMALADIN_OBS
    do_adjust_aladin = .false.
    nulnam=0
    ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=namaladin_obs,iostat=ierr)
    if (ierr.ne.0) call utl_abort('oop_zzz_nl: Error reading namelist')
    if (.not.beSilent) write(*,nml=namaladin_obs)
    ierr=fclos(nulnam)

    BODY: do
      bodyIndex = obs_getBodyIndex(obsSpaceData)
      if (bodyIndex < 0) exit BODY

      ! Process all geometric-height data within the domain of the model
      if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated .or.  &
          obs_bodyElem_i(obsSpaceData,OBS_XTR,bodyIndex) /= 0 .or.  &
          obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) /= 1 ) &
        cycle BODY
      ! So, OBS_VCO==1 => OBS_PPP is a height in m

      bufrCode=obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
      zvar=obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex)
      zlev=obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
      headerIndex=obs_bodyElem_i(obsSpaceData,OBS_HIND,bodyIndex)

      ilyr = obs_bodyElem_i(obsSpaceData,OBS_LYR,bodyIndex)
      varLevel = vnl_varLevelFromVarnum(bufrCode)
      zpt= col_getHeight(columnTrlOnTrlLev,ilyr  ,headerIndex,varLevel)
      zpb= col_getHeight(columnTrlOnTrlLev,ilyr+1,headerIndex,varLevel)
      zwb  = (zpt-zlev)/(zpt-zpb)
      zwt  = 1.d0 - zwb

      select case (bufrCode)
      case (bufr_neal) ! Aladin HLOS wind observation
        ! Scan body indices for the needed attributes
        found = 0
        bodyIndexStart = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
        bodyIndexEnd   = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) &
                       + bodyIndexStart - 1
        tempRef = 0.0d0
        dwdt    = 0.0d0
        presRef = 0.0d0
        dwdp    = 0.0d0
        BODY_SUPP: do bodyIndex2 = bodyIndexStart, bodyIndexEnd

          value = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex2)
          select case(obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2))
          case(bufr_NEAZ)
            azimuth = value * MPC_RADIANS_PER_DEGREE_R8
            found = ibset(found,0)

          case(bufr_nett)
            tempRef = value
            found = ibset(found,1)

          case(bufr_neps)
            presRef = value
            found = ibset(found,2)

          case(bufr_NEDWDP)
            dwdp = value
            found = ibset(found,3)

          case(bufr_NEDWDT)
            dwdt = value
            found = ibset(found,4)
          end select

          if (popcnt(found) == 5) exit BODY_SUPP
        end do BODY_SUPP

        if (.not. btest(found,0))then
          ! The azimuth was not found.  The observation cannot be treated
          ! Set the assimilation flag to 0 to ignore this datum later.
          call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
          cycle BODY
        end if

        ! Obtain the needed forecast data
        uuLyr =col_getElem(columnTrlOnTrlLev,ilyr,  headerIndex,'UU')
        uuLyr1=col_getElem(columnTrlOnTrlLev,ilyr+1,headerIndex,'UU')
        vvLyr =col_getElem(columnTrlOnTrlLev,ilyr,  headerIndex,'VV')
        vvLyr1=col_getElem(columnTrlOnTrlLev,ilyr+1,headerIndex,'VV')
        ttLyr =col_getElem(columnTrlOnTrlLev,ilyr,  headerIndex,'TT')
        ttLyr1=col_getElem(columnTrlOnTrlLev,ilyr+1,headerIndex,'TT')
        ppLyr =col_getPressure(columnTrlOnTrlLev,ilyr  ,headerIndex,'MM')
        ppLyr1=col_getPressure(columnTrlOnTrlLev,ilyr+1,headerIndex,'MM')

        ! Interpolate forecast T, P to the observation location
        ttbg  = zwb*ttLyr1 + zwt*ttLyr
        ppbg  = zwb*ppLyr1 + zwt*ppLyr

        ! Adjust zvar, the HLOS wind observation, if all attributes are available
        if (do_adjust_aladin .and. (popcnt(found) == 5)) then
          ! Adjust in situ the HLOS wind data from obsSpaceData to account for
          ! the differences between our T, P forecast fields and those of the NWP
          ! site that calculated the HLOS wind values.  The goal is to produce an
          ! HLOS wind observation as if it had been calculated by us.
          zvar = zvar + (ttbg - tempRef) * dwdt &
                      + (ppbg - presRef) * dwdp
        end if

        ! Apply the nonlinear aladin observation operator
        trlValueBot= -vvLyr1*cos(azimuth) - uuLyr1*sin(azimuth)
        trlValueTop= -vvLyr *cos(azimuth) - uuLyr *sin(azimuth)

        ! For aladin data, the temperature and pressure are really *reference*
        ! values.  They must not be assimilated.  Mark them so.
      case (bufr_nett)
        call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
        cycle BODY
      case (bufr_neps)
        call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex,obs_notAssimilated)
        cycle BODY
      case (bufr_nees)
        call utl_abort('oop_zzz_nl: CANNOT ASSIMILATE ES!!!')

      case default
        ! These are the profiler observations
        levIndexTop = ilyr + col_getOffsetFromVarno(columnTrlOnTrlLev,bufrCode)
        levIndexBot = levIndexTop+1
        trlValueBot=col_getElem(columnTrlOnTrlLev,levIndexBot,headerIndex)
        trlValueTop=col_getElem(columnTrlOnTrlLev,levIndexTop,headerIndex)
      end select

      zomp = zvar-(zwb*trlValueBot+zwt*trlValueTop)
      call obs_bodySet_r(obsSpaceData,destObsColumn,bodyIndex,zomp)

    enddo BODY

  end subroutine oop_zzz_nl

  !--------------------------------------------------------------------------
  ! oop_sfc_nl
  !--------------------------------------------------------------------------
  subroutine oop_sfc_nl( columnTrlOnTrlLev, obsSpaceData, beSilent, cdfam,  &
                         destObsColumn )
    ! :Purpose:  Computation of the residuals to the observations
    !            FOR SURFACE DATA (except ground-based GPS zenith delay).
    implicit none

    ! arguments
    type(struct_columnData) :: columnTrlOnTrlLev
    type(struct_obs)        :: obsSpaceData
    logical                 :: beSilent
    character(len=*)        :: cdfam ! family of observation
    integer                 :: destObsColumn

    ! locals
    integer :: columnLevelIndex, bufrCode, headerIndex, bodyIndex
    integer :: ierr, nulnam, fnom, fclos
    real(8) :: obsValue, trlVirtTemp
    real(8) :: pCorrectionFactor, coeffA, coeffB
    real(8) :: obsHeight, deltaT, delTdelZ, trlLevelHeight
    real(8) :: trlValue
    real(8) :: trlUwind, trlVwind, squareSum, trlWindSpeed 
    character(len=4) :: varLevel
    logical, save :: firstCall = .true.

    ! namelist variables
    logical, save :: adjustTemperature ! choose to adjust near-sfc temperature using lapse rate and height difference

    namelist /namSurfaceObs/adjustTemperature

    if (.not.beSilent) write(*,*) "Entering subroutine oop_sfc_nl"

    ! Read in the namelist namSurfaceObs
    if (firstCall) then
      adjustTemperature = .true. ! default value

      if (utl_isNamelistPresent('namSurfaceObs','./flnml')) then
        nulnam=0
        ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
        read(nulnam,nml=namSurfaceObs,iostat=ierr)
        if (ierr /= 0) call utl_abort('oop_sfc_nl: Error reading namelist namSurfaceObs')
        if (.not. beSilent) write(*,nml=namSurfaceObs)
        ierr=fclos(nulnam)
      else if (.not. beSilent) then
        write(*,*)
        write(*,*) 'oop_sfc_nl: namSurfaceObs is missing in the namelist. The default value will be taken.'
      end if
      firstCall = .false.
    end if

    ! loop over all header indices of the specified family with surface obs
    call obs_set_current_header_list(obsSpaceData,cdfam)
    HEADER: do
       headerIndex = obs_getHeaderIndex(obsSpaceData)
       if (headerIndex < 0) exit HEADER
       ! loop over all body indices for this headerIndex
       call obs_set_current_body_list(obsSpaceData, headerIndex)
       BODY: do 
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY

          ! only process height level observations flagged to be assimilated
          if ( obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) /= 1 .or.  &
               obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) cycle BODY

          ! only process this set of surface observations
          bufrCode=obs_bodyElem_i (obsSpaceData,OBS_VNM,bodyIndex)
          if (bufrCode /= bufr_nets .and. bufrCode /= bufr_neps   .and.  &
              bufrCode /= bufr_neus .and. bufrCode /= bufr_nevs   .and.  &
              bufrCode /= bufr_ness .and. bufrCode /= bufr_nepn   .and.  &
              bufrCode /= bufr_vis  .and. bufrCode /= bufr_logVis .and.  &
              bufrCode /= bufr_gust .and.  &
              bufrCode /= bufr_radarPrecip .and.  &
              bufrCode /= bufr_logRadarPrecip .and. &
              bufrCode /= bufr_nefs ) cycle BODY

          obsValue = obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex)
          obsHeight= obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
          ! obsHeigth = station height + obs validity height offset as defined in obsUtil_mod->surfvcord
          varLevel = vnl_varLevelFromVarnum(bufrCode)

          if (bufrCode == bufr_nets .or. bufrCode == bufr_ness .or.  &
              bufrCode == bufr_neus .or. bufrCode == bufr_nevs .or.  &
              bufrCode == bufr_gust .or. bufrCode == bufr_vis  .or.  &
              bufrCode == bufr_logVis .or. &
              bufrCode == bufr_radarPrecip .or.  &
              bufrCode == bufr_logRadarPrecip ) then

             ! T1.5m,(T-Td)1.5m,u10m,v10m
             ! In this section we always extrapolate linearly the trial
             ! field at the model surface to the height of the
             ! surface observation whether the observation is above or
             ! below the model surface.
             ! NOTE: For (T-Td)1.5m,u10m,v10m we do a zero order extrapolation

             if (bufrCode == bufr_nets .and. adjustTemperature) then
               delTdelZ = temperatureLapseRate
             else
               delTdelZ = 0.0d0
             end if

             if (bufrCode == bufr_ness) then
                trlValue = hutoes(col_getElem(columnTrlOnTrlLev,col_getNumLev(columnTrlOnTrlLev,'TH'),headerIndex,'HU'),  &
                     col_getElem(columnTrlOnTrlLev,col_getNumLev(columnTrlOnTrlLev,'TH'),headerIndex,'TT'),  &
                     col_getPressure(columnTrlOnTrlLev,col_getNumLev(columnTrlOnTrlLev,'TH'),headerIndex,'TH'))
             else
                columnLevelIndex = col_getNumLev(columnTrlOnTrlLev,varLevel) + col_getOffsetFromVarno(columnTrlOnTrlLev,bufrCode)
                trlValue = col_getElem(columnTrlOnTrlLev,columnLevelIndex,headerIndex)
             end if
             trlLevelHeight = col_getHeight(columnTrlOnTrlLev,col_getNumLev(columnTrlOnTrlLev,varLevel),headerIndex,varLevel)

             call obs_bodySet_r(obsSpaceData,destObsColumn,bodyIndex,  &
                  (obsValue-trlValue + delTdelZ*(obsHeight-trlLevelHeight)) )

          else if ( bufrCode == bufr_neps .or. bufrCode == bufr_nepn ) then

            ! Surface (PS) & mean sea level (PN) pressure cases
            ! Background surface pressure are corrected for the height difference with the 
            ! observation. For mean sea level observation, the observation height = 0.

            ! 1) Temperature difference = lapse-rate (6.5 degree/km) * height difference (dz)
            deltaT = temperatureLapseRate*(obsHeight-col_getHeight(columnTrlOnTrlLev,0,headerIndex,'SF'))

            ! 2) Compute the 1.5m background virtual temperature: Tv = T*(1+0.608*HU)
            trlVirtTemp = (1.0d0 + MPC_DELTA_R8 *  &
                   col_getElem(columnTrlOnTrlLev,col_getNumLev(columnTrlOnTrlLev,'TH'),headerIndex,'HU')) *  &
                   col_getElem(columnTrlOnTrlLev,col_getNumLev(columnTrlOnTrlLev,'TH'),headerIndex,'TT')
            
            ! 3) Compute the correction coefficient
            ! The legacy code says...
            coeffA = (trlVirtTemp-deltaT)/trlVirtTemp
            coeffB = 1.0D0/(MPC_RGAS_DRY_AIR_R8*temperatureLapseRate/ec_rg)
            pCorrectionFactor = coeffA**coeffB
            ! However, the U.S. Standard Atmosphere (1976, U.S. Government Printing Office, Washington, D.C., 1976*)
            ! at page 12 says...
            ! coeffA = trlVirtTemp/(trlVirtTemp+deltaT)
            ! but the former was found to perform better (gives lower O-P values) than the latter by J-F Caron in 2018

            ! 4) O-P, where P = P0 * pCorrectionFactor (Eq. 33a at page 12 of the U.S. Standard Atmosphere, 1976, 
            !                                           U.S. Government Printing Office, Washington, D.C., 1976*)
            call obs_bodySet_r(obsSpaceData,destObsColumn,bodyIndex,  &
                               obsValue-(col_getElem(columnTrlOnTrlLev,1,headerIndex,'P0')*pCorrectionFactor))

            ! (*) available at https://ntrs.nasa.gov/api/citations/19770009539/downloads/19770009539.pdf

          else if (bufrCode == bufr_nefs) then 

            trlUwind = col_getElem(columnTrlOnTrlLev,col_getNumLev(columnTrlOnTrlLev,'MM'),headerIndex,'UU')
            trlVwind = col_getElem(columnTrlOnTrlLev,col_getNumLev(columnTrlOnTrlLev,'MM'),headerIndex,'VV')
            squareSum = trlUwind**2 + trlVwind**2
            if ( squareSum > 1.d-10 ) then
              trlWindSpeed = sqrt(squareSum)
            else
              trlWindSpeed = 0.0
            end if
            call obs_bodySet_r(obsSpaceData,OBS_OMP,bodyIndex,  &
                              obsValue-trlWindSpeed)  

          end if

       end do BODY

    end do HEADER

  end subroutine oop_sfc_nl

  !--------------------------------------------------------------------------
  ! oop_sst_nl
  !--------------------------------------------------------------------------
  subroutine oop_sst_nl( columnTrlOnTrlLev, obsSpaceData, beSilent, cdfam,  &
                         destObsColumn )
    ! :Purpose: Computation of Jo and the residuals to the observations
    !           for Sea Surface Temperature (SST) data.
    implicit none

    ! arguments
    type(struct_columnData) :: columnTrlOnTrlLev
    type(struct_obs)        :: obsSpaceData
    logical                 :: beSilent
    character(len=*)        :: cdfam        ! family of observation
    integer                 :: destObsColumn

    ! locals
    integer          :: bufrCode, headerIndex, bodyIndex
    real(8)          :: obsValue
    character(len=4) :: varName

    if (.not.beSilent) write(*,*) "Entering subroutine oop_sst_nl, family: ", trim(cdfam)

    ! loop over all header indices of the specified family with surface obs
    call obs_set_current_header_list( obsSpaceData, cdfam )

    HEADER: do

      headerIndex = obs_getHeaderIndex( obsSpaceData )
      if ( headerIndex < 0 ) exit HEADER

      ! loop over all body indices for this headerIndex
      call obs_set_current_body_list( obsSpaceData, headerIndex )

      BODY: do

        bodyIndex = obs_getBodyIndex( obsSpaceData )
        if ( bodyIndex < 0 ) exit BODY

        ! only process observations flagged to be assimilated
        if ( obs_bodyElem_i( obsSpaceData, OBS_ASS, bodyIndex ) /= obs_assimilated ) cycle BODY

        bufrCode = obs_bodyElem_i( obsSpaceData, OBS_VNM, bodyIndex )

        if ( bufrCode /= bufr_sst ) cycle BODY
        
        if ( col_varExist(columnTrlOnTrlLev,'TM') ) then
          varName = 'TM'
        else
          varName = 'TG'
        end if

        obsValue = obs_bodyElem_r( obsSpaceData, OBS_VAR, bodyIndex )
        call obs_bodySet_r( obsSpaceData, destObsColumn, bodyIndex, &
                            obsValue - ( col_getElem( columnTrlOnTrlLev, 1, headerIndex, varName ) ))

      end do BODY

    end do HEADER

  end subroutine oop_sst_nl

  !--------------------------------------------------------------------------
  ! oop_hydro_nl
  !--------------------------------------------------------------------------
  subroutine oop_hydro_nl(columnTrlOnTrlLev, obsSpaceData, beSilent, cdfam,  &
                          destObsColumn)
    ! :Purpose: To computate Jo and the residuals to the observations
    !           for hydrological data
    implicit none

    ! arguments
    type(struct_columnData) :: columnTrlOnTrlLev
    type(struct_obs)        :: obsSpaceData
    logical                 :: beSilent
    character(len=*)        :: cdfam        ! family of observation
    integer                 :: destObsColumn

    ! locals
    integer          :: bufrCode, headerIndex, bodyIndex
    real(8)          :: obsValue
    character(len=4) :: varName

    if (.not.beSilent) write(*,*) "Entering subroutine oop_hydro_nl, family: ", trim(cdfam)

    ! loop over all header indices of the specified family with surface obs
    call obs_set_current_header_list( obsSpaceData, cdfam )

    HEADER: do
      headerIndex = obs_getHeaderIndex( obsSpaceData )
      if ( headerIndex < 0 ) exit HEADER

      ! loop over all body indices for this headerIndex
      call obs_set_current_body_list( obsSpaceData, headerIndex )

      BODY: do
        bodyIndex = obs_getBodyIndex( obsSpaceData )
        if ( bodyIndex < 0 ) exit BODY

        ! only process observations flagged to be assimilated
        if ( obs_bodyElem_i( obsSpaceData, OBS_ASS, bodyIndex ) /= obs_assimilated ) cycle BODY

        bufrCode = obs_bodyElem_i( obsSpaceData, OBS_VNM, bodyIndex )

        if ( bufrCode /= bufr_riverFlow ) cycle BODY

        obsValue = obs_bodyElem_r( obsSpaceData, OBS_VAR, bodyIndex )
        varName = vnl_varNameFromVarNum(bufrCode)
        call obs_bodySet_r( obsSpaceData, destObsColumn, bodyIndex, &
                            obsValue - col_getElem(columnTrlOnTrlLev,1,headerIndex, varName_opt = varName) )

      end do BODY

    end do HEADER

  end subroutine oop_hydro_nl

  !--------------------------------------------------------------------------
  ! oop_ice_nl
  !--------------------------------------------------------------------------
  subroutine oop_ice_nl( columnTrlOnTrlLev, obsSpaceData, beSilent, cdfam,  &
                         destObsColumn )
    ! :Purpose: Computation of Jo and the residuals to the observations
    !           FOR SEA ICE CONCENTRATION DATA
    implicit none

    ! arguments
    type(struct_columnData), intent(in)    :: columnTrlOnTrlLev
    type(struct_obs)       , intent(inout) :: obsSpaceData
    logical                , intent(in)    :: beSilent
    character(len=*)       , intent(in)    :: cdfam        ! family of observation
    integer                , intent(in)    :: destObsColumn

    ! locals
    integer          :: bufrCode, headerIndex, bodyIndex
    integer          :: obsDate, monthIndex
    integer          :: trackCellNum
    real(8)          :: obsValue, backValue
    real(8)          :: conc
    character(len=4) :: varName
    character(len=8) :: ccyymmdd

    if (.not. beSilent) write(*,*) "Entering subroutine oop_ice_nl, family: ", trim(cdfam)

    ! loop over all body indices
    call obs_set_current_body_list( obsSpaceData, cdfam )

    BODY: do

      bodyIndex = obs_getBodyIndex( obsSpaceData )
      if ( bodyIndex < 0 ) exit BODY

      ! only process observations flagged to be assimilated
      if ( obs_bodyElem_i( obsSpaceData, OBS_ASS, bodyIndex ) /= obs_assimilated ) cycle BODY

      bufrCode = obs_bodyElem_i( obsSpaceData, OBS_VNM, bodyIndex )

      headerIndex = obs_bodyElem_i( obsSpaceData, OBS_HIND, bodyIndex )
      varName = vnl_varNameFromVarNum(bufrCode)

      select case (bufrCode)
      case(bufr_icec, bufr_icep)
        backValue = 100.0d0*col_getElem( columnTrlOnTrlLev, 1, headerIndex, varName )
      case(bufr_icev)
        backValue = 1.0d0*col_getElem( columnTrlOnTrlLev, 1, headerIndex, varName )
      case(bufr_ices)
        obsDate = obs_headElem_i( obsSpaceData, OBS_DAT, headerIndex ) 
        write(ccyymmdd, FMT='(i8.8)') obsDate
        read(ccyymmdd(5:6), FMT='(i2)') monthIndex
        conc = col_getElem( columnTrlOnTrlLev, 1, headerIndex, varName)
        trackCellNum = obs_headElem_i( obsSpaceData, OBS_FOV, headerIndex )
        backValue = (1.0d0-conc)*oer_ascatAnisOpenWater(trackCellNum,monthIndex) + &
                         conc*oer_ascatAnisIce(trackCellNum,monthIndex)
      case default
        cycle BODY
      end select

      obsValue = obs_bodyElem_r( obsSpaceData, OBS_VAR, bodyIndex )
      call obs_bodySet_r( obsSpaceData, destObsColumn, bodyIndex, obsValue - backValue )

    end do BODY

  end subroutine oop_ice_nl

  !--------------------------------------------------------------------------
  ! oop_raDvel_nl
  !--------------------------------------------------------------------------
  subroutine oop_raDvel_nl(columnTrlOnTrlLev, obsSpaceData, beSilent, cdfam, &
                           destObsColumn)
    ! :Purpose: Computation of Jo and OMP for Radar Doppler velocity observations 
    implicit none

    ! arguments
    type(struct_columnData), intent(in)    :: columnTrlOnTrlLev
    type(struct_obs)       , intent(inout) :: obsSpaceData
    logical                , intent(in)    :: beSilent
    character(len=*)       , intent(in)    :: cdfam        ! family of observation
    integer                , intent(in)    :: destObsColumn

    ! locals
    integer :: bodyIndex, headerIndex, levelIndex, bufrCode
    real(8) :: observedDoppler, simulatedDoppler
    real(8) :: levelAltLow, levelAltHigh
    real(8) :: radarAltitude, beamAzimuth, beamElevation, obsAltitude
    real(8) :: uuLow, uuHigh, vvLow, vvHigh, uuInterpolated, vvInterpolated
    real(8) :: interpWeight

    call obs_set_current_header_list(obsSpaceData, cdfam)
    if (.not.beSilent) write(*,*) 'Entering subroutine oop_raDvel_nl, family: ', trim(cdfam)

    
    !
    ! Loop over all header indices of the 'RA' family with schema 'radvel':
    !
    HEADER: do  

      headerIndex = obs_getHeaderIndex(obsSpaceData)  
      if (headerIndex < 0) exit HEADER
  
      radarAltitude = obs_headElem_r(obsSpaceData, OBS_ALT,  headerIndex)
      beamAzimuth   = obs_headElem_r(obsSpaceData, OBS_RZAM, headerIndex) 
      beamElevation = obs_headElem_r(obsSpaceData, OBS_RELE, headerIndex)
      call obs_set_current_body_list(obsSpaceData, headerIndex)
      !
      ! Loop over all body indices of the 'RA' family with schema 'radvel':
      !
      BODY: do
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if (bodyIndex < 0) exit BODY
        ! Check that this observation has the expected bufr element ID
        bufrCode = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex)
        if (bufrCode /= bufr_radvel) cycle BODY
        ! only process observations flagged to be assimilated
        if (obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) /= obs_assimilated) cycle BODY

        obsAltitude  = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex) 

        ! Levels that bracket the observation from OBS_LYR
        !   note to self:   like in GEM, level=1 is the highest level
        levelIndex = obs_bodyElem_i(obsSpaceData, OBS_LYR, bodyIndex)

        levelAltHigh = col_getHeight(columnTrlOnTrlLev, levelIndex,   headerIndex,'MM')
        levelAltLow  = col_getHeight(columnTrlOnTrlLev, levelIndex+1, headerIndex,'MM')

        ! Vertical interpolation of model wind at observation height
        interpWeight = (obsAltitude - levelAltLow)/(levelAltHigh - levelAltLow)
        uuHigh = col_getElem(columnTrlOnTrlLev, levelIndex,   headerIndex, 'UU')
        uuLow  = col_getElem(columnTrlOnTrlLev, levelIndex+1, headerIndex, 'UU')
        vvHigh = col_getElem(columnTrlOnTrlLev, levelIndex,   headerIndex, 'VV')
        vvLow  = col_getElem(columnTrlOnTrlLev, levelIndex+1, headerIndex, 'VV')
        uuInterpolated = uuLow + interpWeight*(uuHigh - uuLow)
        vvInterpolated = vvLow + interpWeight*(vvHigh - vvLow)

        ! Doppler velocity is the projection of wind along direction of radar beam
        ! Positive values indicates velocities "away" from the radar
        simulatedDoppler = uuInterpolated*sin(beamAzimuth) + vvInterpolated*cos(beamAzimuth)

        observedDoppler = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex)     

        call obs_bodySet_r(obsSpaceData, destObsColumn, bodyIndex, observedDoppler-simulatedDoppler)

      end do BODY
    end do HEADER
    if (.not. beSilent) write(*,*) 'Ending subroutine oop_raDvel_nl, family: ', trim(cdfam)

  end subroutine oop_raDvel_nl

  !--------------------------------------------------------------------------
  ! oop_gpsro_nl
  !--------------------------------------------------------------------------
  subroutine oop_gpsro_nl(columnTrlOnTrlLev, obsSpaceData, beSilent, destObsColumn)
    ! :Purpose: Computation of Jo and the residuals to the GPSRO observations
    !
    ! :Note: gps_struct1sw_v2 allows calculation of partial derivatives of refractivity 
    !        in gps_diff object w.r.t TT/HU/GZ/P0. The indirect dependency refractivity 
    !        to TT/HU/P0 through GZ is now attributed to direct dependency of refractivity on GZ.
    implicit none

    ! Arguments
    type(struct_columnData) :: columnTrlOnTrlLev
    type(struct_obs)        :: obsSpaceData
    type(struct_vco), pointer :: vco_hr
    logical                 :: beSilent
    integer                 :: destObsColumn

    ! Locals
    real(8) :: pjob, pjo1
    real(8) :: zlat, lat
    real(8) :: zlon, lon
    real(8) :: zazm, azm
    integer :: isat, iclf
    real(8) :: rad, geo
    real(8), allocatable :: zpp(:)
    real(8), allocatable :: ztt(:)
    real(8), allocatable :: zhu(:)
    real(8), allocatable :: zHeight(:)
    real(8), allocatable :: zuu(:)
    real(8), allocatable :: zvv(:)
    real(8) :: zp0, zmt
    real(8) :: hnh1, zobs, zmhx, zoer, zinc
    integer headerIndex, idatyp, bodyIndex
    integer jl, ngpslev, nwndlev
    logical  assim, firstheader, ldsc
    integer :: nh, nh1, iProfile, varNum
    type(gps_profile)           :: prf
    real(8)       , allocatable :: h   (:),azmv(:)
    type(gps_diff), allocatable :: rstv(:)
    integer :: Vcode

    if (.not.beSilent) write(*,*)'ENTER oop_gpsro_nl'

    vco_hr => col_getVco(columnTrlOnTrlLev)
    vcode = vco_hr%vcode

    !
    ! Initializations
    !
    ngpslev=col_getNumLev(columnTrlOnTrlLev,'TH')
    nwndlev=col_getNumLev(columnTrlOnTrlLev,'MM')
    allocate(zpp(ngpslev))
    allocate(ztt(ngpslev))
    allocate(zhu(ngpslev))
    allocate(zHeight(ngpslev))
    allocate(zuu(ngpslev))
    allocate(zvv(ngpslev))

    allocate( h    (gps_ro_maxprfsize) )
    allocate( azmv (gps_ro_maxprfsize) )
    allocate( rstv (gps_ro_maxprfsize) )

    !
    ! Loop over all header indices of the 'RO' family:
    !
    call obs_set_current_header_list(obsSpaceData,'RO')
    firstheader = .true.

    HEADER: do
       headerIndex = obs_getHeaderIndex(obsSpaceData)
       if (headerIndex < 0) exit HEADER
       !
       ! Process only refractivity data (codtyp 169)
       !
       idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
       if ( idatyp /= 169 ) cycle HEADER
       !
       ! Scan for requested data values of the profile, and count them
       !
       assim = .false.
       nh = 0
       call obs_set_current_body_list(obsSpaceData, headerIndex)
       BODY: do 
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY
          if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated ) then
             assim = .true.
             nh = nh + 1
          end if
       end do BODY
       !
       ! If no assimilations are requested, skip to next header
       !
       if (.not.assim) cycle HEADER
       iProfile = gps_iprofile_from_index(headerIndex)
       varNum   = gps_vRO_IndexPrf(iProfile, 2)
       !
       ! Basic geometric variables of the profile:
       !
       isat = obs_headElem_i(obsSpaceData,OBS_SAT,headerIndex)
       iclf = obs_headElem_i(obsSpaceData,OBS_ROQF,headerIndex)
       rad  = obs_headElem_r(obsSpaceData,OBS_TRAD,headerIndex)
       geo  = obs_headElem_r(obsSpaceData,OBS_GEOI,headerIndex)
       azm  = obs_headElem_r(obsSpaceData,OBS_AZA,headerIndex)
       zmt  = col_getHeight(columnTrlOnTrlLev,0,headerIndex,'SF')
       !
       ! Profile at the observation location:
       !
       zlat = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)
       zlon = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex)
       lat  = zlat * MPC_DEGREES_PER_RADIAN_R8
       lon  = zlon * MPC_DEGREES_PER_RADIAN_R8
       zazm = azm / MPC_DEGREES_PER_RADIAN_R8
       zp0  = col_getElem(columnTrlOnTrlLev,1,headerIndex,'P0')
       do jl = 1, ngpslev
          !
          ! Profile x
          !
          zpp(jl) = col_getPressure(columnTrlOnTrlLev,jl,headerIndex,'TH')
          ztt(jl) = col_getElem(columnTrlOnTrlLev,jl,headerIndex,'TT') - MPC_K_C_DEGREE_OFFSET_R8
          zhu(jl) = col_getElem(columnTrlOnTrlLev,jl,headerIndex,'HU')
          zHeight(jl) = col_getHeight(columnTrlOnTrlLev,jl,headerIndex,'TH')
       end do

       if (Vcode == 5002) then
          ! case with top thermo level above top momentum level (Vcode=5002)
          do jl = 1, nwndlev
             zuu(jl) = col_getElem(columnTrlOnTrlLev,jl  ,headerIndex,'UU')
             zvv(jl) = col_getElem(columnTrlOnTrlLev,jl  ,headerIndex,'VV')
          end do
       elseif (Vcode == 5005 .or. Vcode == 5100 .or. Vcode == 21001) then
          ! case without top thermo above top momentum level or unstaggered (Vcode=5001/4/5)
          do jl = 1, nwndlev-1
             zuu(jl) = col_getElem(columnTrlOnTrlLev,jl+1,headerIndex,'UU')
             zvv(jl) = col_getElem(columnTrlOnTrlLev,jl+1,headerIndex,'VV')
          end do
          zuu(nwndlev) = zuu(nwndlev-1)
          zvv(nwndlev) = zuu(nwndlev-1)
       else 
          ! case with Vcode /= 5002 and Vcode /= 5005 and Vcode /= 5100
          call utl_abort('oop_gpsro_nl: invalid vertical coord!')
       end if
       zuu(ngpslev) = zuu(nwndlev)
       zvv(ngpslev) = zuu(nwndlev)
       !
       ! GPS profile structure:
       !
       call gps_struct1sw_v2(ngpslev,zLat,zLon,zAzm,zMT,Rad,geo,zP0,zPP,zTT,zHU,zUU,zVV,zHeight,prf)
       ldsc=.not.btest(iclf,16-3)
       !
       ! Prepare the vector of all the observations:
       !
       nh1 = 0
       !
       ! Loop over all body indices for this headerIndex:
       ! (start at the beginning of the list)
       !
       call obs_set_current_body_list(obsSpaceData, headerIndex)
       BODY_2: do
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY_2
          if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated ) then
             nh1      = nh1 + 1
             h(nh1)   = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
             azmv(nh1)= zazm
          end if
       end do BODY_2
       !
       ! Apply the observation operator:
       !
       ! varNum = bufr_nebd (15037) or varNum = bufr_nerf (15036) for GPS-RO
       iProfile = gps_iprofile_from_index(headerIndex)
       if (varNum == bufr_nebd) then
          call gps_bndopv1(h, azmv, nh, prf, rstv)
       else
          call gps_refopv (h,       nh, prf, rstv)
       end if
       !
       ! Perform the (H(x)-Y)/S operation:
       !
       nh1 = 0
       pjob = 0.d0
       !
       ! Loop over all body indices for this headerIndex:
       ! (start at the beginning of the list)
       !
       call obs_set_current_body_list(obsSpaceData, headerIndex)
       BODY_3: do 
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY_3
          if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated ) then
             nh1 = nh1 + 1
             !
             ! Altitude:
             !
             hnh1= obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
             if (varNum == bufr_nebd) hnh1=hnh1-rad
             !
             ! Observation operator H(x)
             !
             zmhx = rstv(nh1)%var
             !
             ! Observation value    Y
             !
             zobs = obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex)
             !
             ! Observation error    S
             !
             zoer = obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex)
             !
             ! Normalized increment
             !
             zinc = (zmhx - zobs) / zoer
             !                           
             ! Datum contribution to Jo:
             !
             pjo1 = 0.5d0 * zinc * zinc
             !
             ! Per profile (PJOB) cumulatives:
             !
             pjob= pjob + pjo1
             !
             if (firstheader .and. .not.beSilent) then
                write(*,  &
                     '(A9,i10,3f7.2,f11.1,4f12.6,15f12.4)') 'DOBSGPSRO',  &
                     headerIndex,lat,lon,azm,hnh1,zobs,zoer,  &
                     zmhx,zinc,pjob,prf%gst(ngpslev)%var  
             end if
             call obs_bodySet_r(obsSpaceData,destObsColumn,bodyIndex, zobs - zmhx)
          end if
       end do BODY_3

       if ( .not.beSilent ) write(*,'(A9,i10,2f7.2,f18.10,f12.4,2I6)')  &
            'GPSRO_JO',headerIndex,lat,lon,pjob,zmt,isat,ldsc
       firstheader = .false.
    end do HEADER

    deallocate( rstv )
    deallocate( azmv )
    deallocate( h    )

    deallocate(zvv)
    deallocate(zuu)
    deallocate(zhu)
    deallocate(zHeight)
    deallocate(ztt)
    deallocate(zpp)

    if (.not.beSilent) write(*,*)'EXIT oop_gpsro_nl'

  end subroutine oop_gpsro_nl

  !--------------------------------------------------------------------------
  ! oop_gpsgb_nl
  !--------------------------------------------------------------------------
  subroutine oop_gpsgb_nl( columnTrlOnTrlLev, obsSpaceData, beSilent, &
                           destObsColumn, analysisMode_opt )
    ! :Purpose: Computation of the residuals to the GB-GPS ZTD observations
    implicit none

    ! Arguments
    type(struct_columnData) :: columnTrlOnTrlLev
    type(struct_obs) :: obsSpaceData
    logical           :: beSilent
    integer           :: destObsColumn
    logical, optional :: analysisMode_opt

    ! Locals
    real(8), allocatable :: zpp (:)
    real(8), allocatable :: ztt (:)
    real(8), allocatable :: zhu (:)
    real(8), allocatable :: zHeight (:)
    real(8) :: zlat, lat, zlon, lon
    real(8) :: zp0, zmt, zdzmin
    real(8) :: zobs, zoer, zinc, zhx, zlev
    real(8) :: zdz, zpsobs, zpsmod, zpwmod, zpomp, zpomps
    real(8) :: ztdomp(gps_gb_maxdata)
    real(8) :: bias, std
    integer :: headerIndex, bodyIndex, ioneobs, idatyp, bufrCode, index_ztd, iztd
    integer :: jl, nlev_T, nobs2p
    integer :: icount1, icount2, icount3, icount, icountp
    logical  :: assim, llrej, analysisMode, lfsl
    character(len=12) :: cstnid
    type(gps_profilezd)    :: prf
    type(gps_diff)         :: ztdopv
    !
    !     PW lower limit (mm) and Ps O-P upper limit (Pa) for ZTD assimilation
    !       Note:  1 mb = 100 Pa --> 2.2 mm ZTD
    !
    real(8) :: zpwmin, zpompmx
    data zpwmin    /   2.0d0 /
    data zpompmx   / 200.0d0 /
    !
    !     Criteria to select single observation (1-OBS mode)
    !
    !     Minimum value for ZTD O-P (m)
    real(8) :: xompmin
    data xompmin   / 0.015d0 /
    ! Minimum value for background (trial) PW (mm)
    real(8) :: xpwmin
    data xpwmin    / 20.0d0  /
    ! Maximum height difference between observation and background surface (m)
    real(8) :: xdzmax
    data xdzmax    / 400.0d0 /

    if (.not.beSilent) write(*,*)'ENTER oop_gpsgb_nl'

    if (present(analysisMode_opt)) then
       analysisMode = analysisMode_opt
    else
       analysisMode = .true.
    end if

    zpomps = 0.0d0

    ! Ensure Jacobian-related arrays are not allocated to force them to be recalculated in oop_H
    if (allocated(oop_vZTD_Jacobian)) deallocate(oop_vZTD_Jacobian)

    zdzmin = gps_gb_dzmin      
    nobs2p = 50

    nlev_T = col_getNumLev(columnTrlOnTrlLev,'TH')
    if (gps_gb_ltestop .and. .not.beSilent) write(*,*) '  col_getNumLev[columnTrlOnTrlLev,TH] = ',nlev_T

    !
    ! Initializations
    !
    allocate(ztt(nlev_T))
    allocate(zhu(nlev_T))
    allocate(zHeight(nlev_T))
    allocate(zpp(nlev_T))

    if ( .not.beSilent ) then
      write(*, *) ' '
      write(*, *) ' '
      write(*,'(A11,A9,3A8,A9,4A8,2A9,A7,A10,A11)')  &
           'OOP_GPSGB_NL','CSTNID','ZLAT','ZLON','ZLEV','ZDZ','ZOBS','ZOER','ZHX','O-P',  &
           'ZPOMPS','ZPOMP','ZPWMOD','ZINC2'
    end if

    icount  = 0
    icount1 = 0
    icount2 = 0
    icount3 = 0
    icountp = 0
    ioneobs = -1

    ! loop over all header indices of the 'GP' family (all obs locations/times)
    call obs_set_current_header_list(obsSpaceData,'GP')

    HEADER: do
       headerIndex = obs_getHeaderIndex(obsSpaceData)
       if (headerIndex < 0) exit HEADER

       ! Process only GP data (codtyp 189)
       idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
       if ( idatyp /= 189 ) cycle HEADER

       assim = .false.
       zpsobs = -100.0d0
       lfsl = .false.
       cstnid = obs_elem_c(obsSpaceData,'STID',headerIndex)
       if (index(cstnid,'FSL_') > 0 .or. index(cstnid,'-NOAA') > 0) lfsl = .true.

       ! Scan for requested ZTD assimilation.
       ! Get GPS antenna height ZLEV and Ps(ZLEV) (ZPSOBS)
       !
       ! loop over all body indices for this headerIndex (observations at location/time)
       call obs_set_current_body_list(obsSpaceData, headerIndex)
       BODY: do 
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY
          bufrCode = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
          if ( (bufrCode == bufr_nezd) .and. &
               (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) ) then
             zlev = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
             assim = .true.
             ! Index in body of ZTD datum (assume at most 1 per header)
             index_ztd = bodyIndex
             icount = icount + 1
          end if
          if ( bufrCode == bufr_neps ) then
             if ( (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) .or. gps_gb_llblmet ) then
                zpsobs = obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex)
                zpomps = obs_bodyElem_r(obsSpaceData,destObsColumn,bodyIndex)
             end if
          end if
       end do BODY

       ! If no ZTD assimilation requested, jump to next header
       if (.not.assim) cycle HEADER

       ! Profile at the observation location:
       lat  = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)
       lon  = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex)
       zlat = lat * MPC_DEGREES_PER_RADIAN_R8
       zlon = lon * MPC_DEGREES_PER_RADIAN_R8
       zmt  = col_getHeight(columnTrlOnTrlLev,0,headerIndex,'SF')
       zp0  = col_getElem(columnTrlOnTrlLev,1,headerIndex,'P0')
       do jl = 1, nlev_T
          zpp(jl) = col_getPressure(columnTrlOnTrlLev,jl,headerIndex,'TH')
          ztt(jl) = col_getElem(columnTrlOnTrlLev,jl,headerIndex,'TT')-MPC_K_C_DEGREE_OFFSET_R8
          zhu(jl) = col_getElem(columnTrlOnTrlLev,jl,headerIndex,'HU')
          zHeight(jl) = col_getHeight(columnTrlOnTrlLev,jl,headerIndex,'TH')
       end do
       zdz = zlev - zmt

       ! Fill GPS ZTD profile structure (PRF):
       call gps_structztd_v2(nlev_T,lat,lon,zmt,zp0,zpp,ztt,zhu,zHeight,gps_gb_lbevis,gps_gb_irefopt,prf)

       ! Apply the GPS ZTD observation operator
       ! --> output is model ZTD (type gps_diff) and P at obs height ZLEV
       call gps_ztdopv(zlev,prf,gps_gb_lbevis,zdzmin,ztdopv,zpsmod,gps_gb_iztdop)

       ! Get model profile PW
       call gps_pw(prf,zpwmod)
       ! ZTD (m)
       zhx    = ztdopv%var

       ! If analysis mode, reject ZTD data for any of the following conditions:
       !    (1) the trial PW is too low (extremely dry) 
       !    and if gps_gb_LASSMET=true and for NOAA/FSL sites only:
       !      (2) Ps observation is missing or out of normal range
       !      (3) the ABS(Ps(obs)-Ps(mod)) difference is too large
       llrej = .false.
       zpomp = -9999.0D0
       if ( analysisMode ) then
          llrej = ( zpwmod < zpwmin )
          if ( gps_gb_lassmet .and. lfsl ) then
             if ( .not. llrej ) then
                if ( zpsobs > 40000.0d0 .and. zpsobs <= 110000.0d0 ) then
                   zpomp = zpsobs - zpsmod
                   llrej = ( abs(zpomp) > zpompmx )
                   if ( llrej ) icount3 = icount3 + 1
                else
                   llrej = .true.
                   icount2 = icount2 + 1
                end if
             else
                icount1 = icount1 + 1
             end if
          end if
       end if

       if ( llrej ) then
          call obs_bodySet_i(obsSpaceData,OBS_ASS,index_ztd, obs_notAssimilated)
          if ( .not. gps_gb_lassmet ) icount1 = icount1 + 1
       end if

       ! Perform the (H(x)-Y)/SDERR operation
       !
       ! loop over all body indices for this headerIndex
       call obs_set_current_body_list(obsSpaceData, headerIndex)
       BODY_2: do 
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY_2
          bufrCode = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
          if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated .and.  &
               bufrCode == bufr_nezd ) then
             icountp = icountp + 1
             !
             ! Observation value    Y
             !
             zobs = obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex)
             !
             ! Observation error    SDERR
             !
             zoer = obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex)
             if ( zoer <= 0.0d0 ) then
                write(*,*) ' Problem with ZTD observation error!'
                write(*,*) ' Station =',cstnid
                write(*,*) ' Error =', zoer
                call utl_abort('OOP_GPSGB_NL: ABORT! BAD ZTD OBSERR') 
             end if

             ! Observation height (m)
             !
             zlev = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
             !
             ! Normalized increment ZINC
             !
             ztdomp(icountp) = zobs - zhx
             zinc  = (zhx - zobs) / zoer
             call obs_bodySet_r(obsSpaceData,destObsColumn,bodyIndex, zobs - zhx)

             !
             ! Apply data selection criteria for 1-OBS Mode
             !
             if ( gps_gb_l1obs .and. ioneobs == -1 ) then
                if ( (zobs-zhx) > xompmin .and. zpwmod > xpwmin .and.  &
                     abs(zdz) < xdzmax ) then
                   ioneobs = headerIndex
                   write(*,*) 'SINGLE OBS SITE = ',cstnid
                end if
             end if
             !
             ! Print data for first NOBS2P observations
             !
             if ( icountp <= nobs2p .and. .not.beSilent ) then
                write(*,  &
                     '(A14,A9,3(1x,f7.2),1x,f8.2,4(1x,f8.5),2(1x,f8.4),2x,f5.2,1x,f9.2,1x,f10.5)')  &
                     'OOP_GPSGB_NL: ',cstnid,zlat,zlon,zlev,zdz,zobs,zoer/gps_gb_yzderrwgt,zhx,-zinc*zoer,  &
                     zpomps/100.d0,zpomp/100.d0,zpwmod,zinc/zoer
             end if

          end if

       end do BODY_2

    end do HEADER

    deallocate(ztt)
    deallocate(zhu)
    deallocate(zHeight)
    deallocate(zpp)

    if ( .not. beSilent ) then
      write(*,*) ' '
      write(*,*) 'NUMBER OF GPS ZTD DATA FLAGGED FOR ASSIMILATION = ', icountp
      if ( icountp > 0 ) then
         bias = sum(ztdomp(1:icountp))/real(icountp,8)
         std = 0.d0
         do jl = 1, icountp
            std = std + (ztdomp(jl)-bias)**2
         end do
         write(*, *) '     MEAN O-P (BIAS) [mm] = ', bias*1000.d0
         if (icountp > 1) then
            std = sqrt(std/(real(icountp,8)-1.d0))
            write(*, *) '     STD  O-P        [mm] = ', std*1000.d0
         else
            write(*, *) '     STD  O-P        Uncomputable since number of GPS ZTD observations is 1'
         end if
         write(*, *) ' '
      end if
    end if

    if ( gps_gb_l1obs .and. analysisMode ) then
       ! Set assim flag to 0 for all observations except for selected record (site/time)
       if ( ioneobs /= -1 ) then
          call obs_set_current_header_list(obsSpaceData,'GP')
          icountp = 1
          HEADER_1: do
             headerIndex = obs_getHeaderIndex(obsSpaceData)
             if (headerIndex < 0) exit HEADER_1
             if (headerIndex /= ioneobs ) then
                call obs_set_current_body_list(obsSpaceData, headerIndex)
                BODY_1: do 
                   bodyIndex = obs_getBodyIndex(obsSpaceData)
                   if (bodyIndex < 0) exit BODY_1
                   call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex, obs_notAssimilated)
                end do BODY_1
             end if
          end do HEADER_1
       else
          call utl_abort('ERROR: FAILED TO SELECT SINGLE OBSERVATION!')
       end if
    end if

    gps_gb_numztd = icountp

    if ( analysisMode .and. icount > 0 .and. .not.gps_gb_l1obs .and. .not.beSilent ) then
       write(*,*) ' '
       write(*,*) '-----------------------------------------'
       write(*,*) ' SUMMARY OF ZTD REJECTIONS IN OOP_GPSGB_NL '
       write(*,*) '-----------------------------------------'
       write(*,*) ' TOTAL NUMBER OF ZTD DATA ORIGINALLY FLAGGED FOR ASSMILATION = ', icount
       write(*,*) '       NUMBER OF ZTD DATA       REJECTED DUE TO LOW TRIAL PW = ', icount1
       write(*,*) '       NUMBER OF ZTD DATA       REJECETD DUE TO    NO PS OBS = ', icount2
       write(*,*) '       NUMBER OF ZTD DATA       REJECETD DUE TO LARGE PS O-P = ', icount3
       write(*,*) ' TOTAL NUMBER OF REJECTED ZTD DATA                           = ', icount1+icount2+icount3
       write(*,*) '       PERCENT   REJECTED                                    = ',   &
            (real(icount1+icount2+icount3,8) / real(icount,8))*100.0d0
       write(*, *) ' TOTAL NUMBER OF ASSIMILATED ZTD DATA                        = ', icountp
       write(*,*) ' '
    end if

    if ( icount > 0 .and. gps_gb_numztd > 0) then

       if ( analysisMode ) then
          if ( .not.beSilent ) write(*,*) ' Number of GPS ZTD data to be assimilated (gps_gb_numZTD) = ', gps_gb_numztd
       else
          if ( .not.beSilent ) write(*,*) ' Number of GPS ZTD data for background check (gps_gb_numZTD) = ', gps_gb_numztd
       end if

       if ( .not.beSilent ) write(*,*) ' Allocating and setting gps_ZTD_Index(gps_gb_numZTD)...'
       if (allocated(gps_ZTD_Index)) deallocate(gps_ZTD_Index)
       allocate(gps_ZTD_Index(gps_gb_numztd))
       iztd = 0
       call obs_set_current_header_list(obsSpaceData,'GP')
       HEADER_2: do
          headerIndex = obs_getHeaderIndex(obsSpaceData)
          if (headerIndex < 0) exit HEADER_2
          idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
          if ( idatyp == 189 ) then
             call obs_set_current_body_list(obsSpaceData, headerIndex)
             BODY_3: do 
                bodyIndex = obs_getBodyIndex(obsSpaceData)
                if (bodyIndex < 0) exit BODY_3
                bufrCode = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
                if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated .and.  &
                     bufrCode == bufr_nezd ) then  
                   iztd = iztd + 1
                   gps_ZTD_Index(iztd) = headerIndex
                end if
             end do BODY_3
          end if
       end do HEADER_2

       if ( iztd /= gps_gb_numztd ) then
          call utl_abort('ERROR: gps_ZTD_Index init: iztd /= gps_gb_numZTD!')
       end if

    end if

    if (.not.beSilent) write(*,*)'EXIT oop_gpsgb_nl'

  end subroutine oop_gpsgb_nl

  !--------------------------------------------------------------------------
  ! oop_tovs_nl
  !--------------------------------------------------------------------------
  subroutine oop_tovs_nl( columnTrl, obsSpaceData, datestamp, beSilent,  &
                          bgckMode_opt, option_opt, sourceObs_opt, destObs_opt )
    ! :Purpose: Computation of the residuals to the tovs observations
    !           option_opt: defines input state:
    !              - 'HR': High Resolution background state,
    !              - 'LR': Low  Resolution background state, (CURRENTLY NOT SUPPORTED)
    !              - 'MO': Model state. (CURRENTLY NOT SUPPORTED)
    implicit none

    ! Arguments
    type(struct_columnData) :: columnTrl
    type(struct_obs) :: obsSpaceData
    integer :: datestamp
    logical :: beSilent
    logical, optional :: bgckMode_opt
    character(len=*), optional :: option_opt       ! only valid value is HR
    integer, optional, intent(in) :: sourceObs_opt ! usually set to OBS_VAR
    integer, optional, intent(in) :: destObs_opt   ! usually set to OBS_OMP

    ! locals
    integer :: jdata, sourceObs, destObs
    logical :: llprint,bgckMode
    character(len=2) :: option

    integer :: channelIndex, tovsIndex
    real(pre_obsReal) :: zdtb, obsPRM
    integer :: idatyp, channelNumber
    integer :: headerIndex, bodyIndex

    if (.not.obs_famExist(obsSpaceData,'TO', localMPI_opt = .true. )) return

    ! 0. set default values if bgckMode, option and source/dest columns not specified
    !

    if (.not.beSilent) write(*,*) "Entering subroutine oop_tovs_nl"

    if (present(bgckMode_opt)) then
       bgckMode = bgckMode_opt
    else
       bgckMode = .false.
    end if

    if (present(option_opt)) then
       option = option_opt(1:2)
    else
       option = 'HR'
    end if
    if ( option /= 'HR' ) call utl_abort('oop_tovs_nl: Invalid option for input state')

    if (present(sourceObs_opt)) then
       sourceObs = sourceObs_opt
    else
       sourceObs = OBS_VAR
    end if

    if (present(destObs_opt)) then
       destObs = destObs_opt
    else
       destObs = OBS_OMP
    end if

    ! 1.   Prepare atmospheric profiles for all tovs observation points for use in rttov
    ! .    -----------------------------------------------------------------------------
    call tvs_fillProfiles(columnTrl,obsSpaceData,datestamp,"nl",beSilent)

    if ( .not.beSilent ) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! 2.   Compute radiance
    ! .    ----------------
    call tvs_rttov(obsSpaceData,bgckMode,beSilent)
    if ( .not.beSilent ) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! 3.   Compute the residuals
    ! .    ----------------------------
    if ( option == 'HR' .or. option == 'LR' ) then
       do jdata=1,obs_numbody(obsSpaceData)
          call obs_bodySet_r(obsSpaceData,OBS_PRM,jdata, obs_bodyElem_r(obsSpaceData,sourceObs,jdata))
       end do
    end if

    ! loop over all header indices of the 'TO' family
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER

      ! Extract general information for this observation point
      !      ------------------------------------------------------

      ! process only radiance data to be assimilated?
      idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
      if ( .not. tvs_isIdBurpTovs(idatyp) ) then
        write(*,*) 'oop_tovs_nl: warning unknown radiance codtyp present check NAMTOVSINST', idatyp
        cycle HEADER
      end if
      tovsIndex = tvs_tovsIndex(headerIndex)
      if ( tovsIndex == -1 ) cycle HEADER

      ! Set the body list
      ! (& start at the beginning of the list)
      call obs_set_current_body_list(obsSpaceData, headerIndex)
      BODY: do
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if ( bodyIndex < 0 ) exit BODY

        ! Only consider if flagged for assimilation
        if ( obs_bodyElem_i(obsSpaceData,obs_ASS,bodyIndex) /= obs_assimilated ) cycle BODY

        call tvs_getChannelNumIndexFromPPP( obsSpaceData, headerIndex, bodyIndex, &
                                            channelNumber, channelIndex )

        if ( channelIndex == 0 ) call utl_abort('oop_tovs_nl: error with channel number')

        zdtb = obs_bodyElem_r(obsSpaceData,OBS_PRM,bodyIndex) - &
             tvs_radiance (tovsIndex) % bt(channelIndex)
        if ( tvs_debug ) then
          obsPRM = obs_bodyElem_r(obsSpaceData,OBS_PRM,bodyIndex)
          write(*,'(a,i4,2f8.2,f6.2)') ' channelNumber,sim,obs,diff= ', &
               channelNumber,  tvs_radiance (tovsIndex) % bt(channelIndex), &
               obsPRM, -zdtb
        end if
        call obs_bodySet_r(obsSpaceData,destObs,bodyIndex, zdtb)

        ! inflate OBS_OER for all-sky assimilation
        call oer_inflateErrAllsky(obsSpaceData, bodyIndex, destObs, beSilent_opt=.true.)

      end do BODY

    end do HEADER

    if (option == 'HR') then
      llprint = .true.
    else
      llprint = .false.
    end if
    if ( beSilent ) llprint = .false.
    if ( llprint ) call tvs_printDetailledOmfStatistics(obsSpaceData)

  end subroutine oop_tovs_nl

  !--------------------------------------------------------------------------
  ! oop_chm_nl
  !--------------------------------------------------------------------------
  subroutine oop_chm_nl( columnTrlOnTrlLev, obsSpaceData, destObsColumn )
    ! :Purpose: Computation of the residuals to the observations
    !           for all observations of the CH (chemical constituents) family.
    !           The array columnTrlOnTrlLev contains the input model array.
    !           Stores OmP in OBS_OMP in obsSpaceData.
    implicit none

    ! Arguments
    type(struct_columnData) :: columnTrlOnTrlLev
    type(struct_obs)        :: obsSpaceData
    integer                 :: destObsColumn
    
    if (.not.obs_famExist(obsSpaceData,'CH', localMPI_opt = .true. )) return

    if ( destObsColumn /= obs_omp ) then
      write(*,*) 'oop_chm_nl: WARNING: Storing results in an obs column other than OBS_OMP. Not fully implemented.'
    end if

    call oopc_CHobsoperators( columnTrlOnTrlLev,obsSpaceData,'nl', & ! 'nl' for non-linear operator
                              destObsColumn_opt=destObsColumn )

  end subroutine oop_chm_nl

  !--------------------------------------------------------------------------
  ! oop_Htl
  !--------------------------------------------------------------------------
  subroutine oop_Htl( columnAnlInc, columnTrlOnAnlIncLev, obsSpaceData, &
                      min_nsim, initializeLinearization_opt )
    !
    ! :Purpose: Compute simulated observations from profiled model increments.
    !           It returns Hdx in OBS_WORK. Calls the several linear observation operators.
    implicit none

    type(struct_columnData)   :: columnAnlInc, columnTrlOnAnlIncLev
    type(struct_obs)          :: obsSpaceData
    type(struct_vco), pointer :: vco_anl
    integer, intent(in)       :: min_nsim
    logical, optional         :: initializeLinearization_opt 

    logical, save :: initializeLinearization = .true.

    if ( mmpi_myid == 0 ) then
      write(*,*) 'OOP_Htl - Linearized observation operators'
    end if

    vco_anl => col_getVco(columnTrlOnAnlIncLev)

    ! Re-linearlize if it is asked for
    if ( present(initializeLinearization_opt) ) then
      initializeLinearization = initializeLinearization_opt
    end if

    if ( initializeLinearization ) then
      ! Find interpolation layer in model profiles (used by several operators)
      if ( col_getNumLev(columnTrlOnAnlIncLev,'MM') > 1 ) call oop_vobslyrs(columnTrlOnAnlIncLev, obsSpaceData, beSilent=.false.)

      ! Initialize some operators needed by linearized H
      call subasic_obs(columnTrlOnAnlIncLev)

      initializeLinearization = .false.
    end if

    call oop_Hpp()           ! fill in OBS_WORK : Hdx

    call oop_Hsf()           ! fill in OBS_WORK : Hdx

    call oop_Hto()           ! fill in OBS_WORK : Hdx

    call oop_Hro( initializeLinearization_opt=initializeLinearization_opt )


    if ( gps_gb_numZTD > 0 ) then 
      call oop_Hgp( initializeLinearization_opt=initializeLinearization_opt )
    end if

    call oop_Hchm()          ! fill in OBS_WORK : Hdx

    call oop_Hsst()          ! fill in OBS_WORK : Hdx

    call oop_Hice()          ! fill in OBS_WORK : Hdx

    call oop_Hhydro()        ! fill in OBS_WORK : Hdx

    call oop_HheightCoordObs()      ! fill in OBS_WORK : Hdx
 
  contains
    
    !--------------------------------------------------------------------------

    subroutine subasic_obs( columnTrlOnAnlIncLev )
      ! :Purpose:   Initialise background state dependant factors
      !             and vectors for use in TLM and adjoint of
      !             non-linear operator
      implicit none

      type(struct_columnData) :: columnTrlOnAnlIncLev

      ! locals
      type(struct_vco), pointer :: vco_anl
      integer :: jlev,columnIndex,nlev_T,vcode_anl,status
      real(8) :: zhu,one

      if ( .not.col_varExist(columnTrlOnAnlIncLev,'TT') .or. .not.col_varExist(columnTrlOnAnlIncLev,'HU') ) return

      write(*,*) 'subasic_obs: setting up linearized Tv operator'

      vco_anl => col_getVco(columnTrlOnAnlIncLev)
      one=1.0D0
      nlev_T = col_getNumLev(columnTrlOnAnlIncLev,'TH')
      Vcode_anl = vco_anl%vCode

      if ( Vcode_anl /= 5002 .and. Vcode_anl /= 5005 ) then
         call utl_abort('subasic_obs: invalid vertical coord!')
      end if

      ! initialize virtual temperature operator

!$OMP PARALLEL DO PRIVATE(jlev,columnIndex,zhu)
      do jlev = 1, nlev_T
         do columnIndex=1,col_getNumCol(columnTrlOnAnlIncLev)
            zhu=col_getElem(columnTrlOnAnlIncLev,jlev,columnIndex,'HU')
            columnTrlOnAnlIncLev%oltv(1,jlev,columnIndex) = phf_fottva(zhu,one)
            columnTrlOnAnlIncLev%oltv(2,jlev,columnIndex) = phf_folnqva(zhu,col_getElem(columnTrlOnAnlIncLev,  &
                 jlev,columnIndex,'TT'),one)
         end do
      end do
!$OMP END PARALLEL DO

    end subroutine subasic_obs

    !--------------------------------------------------------------------------

    subroutine oop_Hpp()
      ! : Purpose: Compute simulated Upper Air observations from profiled model
      !            increments.
      !            It returns Hdx in OBS_WORK
      !            Interpolate vertically the contents of commvo to
      !            the pressure levels of the observations.
      !            A linear interpolation in ln(p) is performed.
      implicit none

      integer levIndexBot,levIndexTop
      integer headerIndex,INDEX_FAMILY,layerIndex
      integer J,bodyIndex,bufrCode,nlev_T
      real(8) ZDADPS
      real(8) ZWB,ZWT,ZLTV
      real(8) ZLEV,ZPT,ZPB
      real(8) delPT,delPB
      real(8) anlIncValueBot,anlIncValueTop,trlValueBot,trlValueTop
      integer, parameter :: numFamily=3
      character(len=2) :: list_family(numFamily)
      character(len=4) :: varLevel

      list_family(1) = 'UA'
      list_family(2) = 'AI'
      list_family(3) = 'SW'

      FAMILY: do index_family=1,numFamily

         call obs_set_current_body_list(obsSpaceData,list_family(index_family))
         BODY: do 
            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if (bodyIndex < 0) exit BODY

            if (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated .and. &
                obs_bodyElem_i(obsSpaceData,OBS_XTR,bodyIndex) == 0               .and. &
                obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) == 2 ) then
               headerIndex = obs_bodyElem_i(obsSpaceData,OBS_HIND,bodyIndex)
               ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
               bufrCode = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
               varLevel = vnl_varLevelFromVarnum(bufrCode)
               layerIndex   = obs_bodyElem_i(obsSpaceData,OBS_LYR,bodyIndex)
               levIndexTop  = layerIndex + col_getOffsetFromVarno(columnTrlOnAnlIncLev,bufrCode)
               levIndexBot  = levIndexTop+1
               ZPT    = col_getPressure(columnTrlOnAnlIncLev,layerIndex  ,headerIndex,varLevel)
               ZPB    = col_getPressure(columnTrlOnAnlIncLev,layerIndex+1,headerIndex,varLevel)
               delPT  = col_getPressure(columnAnlInc,layerIndex  ,headerIndex,varLevel)
               delPB  = col_getPressure(columnAnlInc,layerIndex+1,headerIndex,varLevel)
               ZWB  = LOG(ZLEV/ZPT)/LOG(ZPB/ZPT)
               ZWT  = 1.0D0 - ZWB

               ZDADPS   = ( LOG(ZLEV/ZPB)*delPT/ZPT -   &
                    LOG(ZLEV/ZPT)*delPB/ZPB )/  &
                    LOG(ZPB/ZPT)**2

               if ( bufrCode == bufr_nees ) then
                  anlIncValueBot=hutoes_tl(col_getElem(columnAnlInc,layerIndex+1,headerIndex,'HU'), &
                       col_getElem(columnAnlInc,layerIndex+1,headerIndex,'TT'), &
                       delPB, &
                       col_getElem(columnTrlOnAnlIncLev,layerIndex+1,headerIndex,'HU'), &
                       col_getPressure(columnTrlOnAnlIncLev,layerIndex+1,headerIndex,'TH'))
                  anlIncValueTop=hutoes_tl(col_getElem(columnAnlInc,layerIndex  ,headerIndex,'HU'), &
                       col_getElem(columnAnlInc,layerIndex  ,headerIndex,'TT'), &
                       delPT, &
                       col_getElem(columnTrlOnAnlIncLev,layerIndex  ,headerIndex,'HU'), &
                       col_getPressure(columnTrlOnAnlIncLev,layerIndex  ,headerIndex,'TH'))
                  trlValueBot=hutoes(col_getElem(columnTrlOnAnlIncLev,layerIndex+1,headerIndex,'HU'), &
                       col_getElem(columnTrlOnAnlIncLev,layerIndex+1,headerIndex,'TT'), &
                       col_getPressure(columnTrlOnAnlIncLev,layerIndex+1,headerIndex,'TH'))
                  trlValueTop=hutoes(col_getElem(columnTrlOnAnlIncLev,layerIndex  ,headerIndex,'HU'), &
                       col_getElem(columnTrlOnAnlIncLev,layerIndex  ,headerIndex,'TT'), &
                       col_getPressure(columnTrlOnAnlIncLev,layerIndex  ,headerIndex,'TH'))
               else
                  anlIncValueBot=col_getElem(columnAnlInc,levIndexBot,headerIndex)
                  anlIncValueTop=col_getElem(columnAnlInc,levIndexTop,headerIndex)
                  trlValueBot=col_getElem(columnTrlOnAnlIncLev,levIndexBot,headerIndex)
                  trlValueTop=col_getElem(columnTrlOnAnlIncLev,levIndexTop,headerIndex)
               end if
               call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex,   &
                    ZWB*anlIncValueBot + ZWT*anlIncValueTop+  &
                    (trlValueBot - trlValueTop)*  &
                    ZDADPS)

            end if

         end do BODY

      end do FAMILY

    end subroutine oop_Hpp

    !--------------------------------------------------------------------------

    subroutine oop_Hsf()
      ! :Purpose: Compute simulated surface observations from profiled model
      !           increments.
      !           It returns Hdx in OBS_WORK
      !
      implicit none
      
      integer :: levIndexBot,levIndexTop
      integer :: headerIndex,layerIndex
      integer :: J, bodyIndex, bufrCode, INDEX_FAMILY, nlev, nLev_T
      real(8) :: coeffA, coeffB
      real(8) :: ZWB, ZWT, ZLTV, trlVirtTemp
      real(8) :: ZPT, ZPB, ZDELPS, ZDELTV, deltaT, obsHeight
      real(8) :: anlIncValue
      real(8) :: delP
      real(8) :: anlIncUwind, anlIncVwind, trlUwind, trlVwind, squareSum, trlWindSpeed, anlIncWindSpeed
      integer, parameter :: numFamily=5
      character(len=2) :: list_family(numFamily)
      character(len=4) :: varLevel
      !
      !     Temperature lapse rate for extrapolation of height below model surface
      !
      coeffB   = 1.0d0/(MPC_RGAS_DRY_AIR_R8*temperatureLapseRate/ec_rg)
      !
      !
      list_family(1) = 'UA'
      list_family(2) = 'SF'
      list_family(3) = 'SC'
      list_family(4) = 'GP'
      list_family(5) = 'RA'

      FAMILY: do index_family=1,numFamily

         call obs_set_current_body_list(obsSpaceData, list_family(index_family))
         BODY: do
            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if (bodyIndex < 0) exit BODY

            ! Process all data within the domain of the model
            bufrCode = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
            if ( bufrCode == bufr_nezd .or. bufrCode == bufr_radvel ) cycle BODY
            if (    (obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) == 1) &
                 .and. (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) &
                 .and. (bufrCode == bufr_nets .or. bufrCode == bufr_neps  &
                 .or. bufrCode == bufr_nepn .or. bufrCode == bufr_ness  &
                 .or. bufrCode == bufr_neus .or. bufrCode == bufr_nevs  &
                 .or. bufrCode == bufr_vis  .or. bufrCode == bufr_logVis  &
                 .or. bufrCode == bufr_gust .or. bufrCode == bufr_nefs &
                 .or. bufrCode == bufr_radarPrecip .or. bufrCode == bufr_logRadarPrecip  &
                 .or. obs_bodyElem_i(obsSpaceData,OBS_XTR,bodyIndex) == 0) ) then

              varLevel    = vnl_varLevelFromVarnum(bufrCode)
              nlev        = col_getNumLev(columnAnlInc,varLevel)
              nlev_T      = col_getNumLev(columnAnlInc,'TH')
              headerIndex = obs_bodyElem_i(obsSpaceData,OBS_HIND,bodyIndex)
              bufrCode    = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
              layerIndex  = obs_bodyElem_i(obsSpaceData,OBS_LYR,bodyIndex)
              obsHeight   = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
              levIndexTop = nlev - 1 + col_getOffsetFromVarno(columnTrlOnAnlIncLev,bufrCode)
              levIndexBot = levIndexTop + 1

              if (bufrCode == bufr_nets .OR. bufrCode == bufr_ness .OR.  &
                  bufrCode == bufr_neus .OR. bufrCode == bufr_nevs .OR. &
                  bufrCode == bufr_vis  .or. bufrCode == bufr_logVis  .or.  &
                  bufrCode == bufr_gust .or.  &
                  bufrCode == bufr_radarPrecip .or. bufrCode == bufr_logRadarPrecip) THEN
                if (bufrCode == bufr_ness ) THEN
                  delP = col_getPressure(columnAnlInc,nlev_T,headerIndex,'TH')
                  anlIncValue = hutoes_tl(col_getElem(columnAnlInc,nlev_T,headerIndex,'HU'), &
                                       col_getElem(columnAnlInc,nlev_T,headerIndex,'TT'), &
                                       delP, &
                                       col_getElem(columnTrlOnAnlIncLev,nlev_T,headerIndex,'HU'), &
                                       col_getPressure(columnTrlOnAnlIncLev,nlev,headerIndex,varLevel))
                else
                  anlIncValue=col_getElem(columnAnlInc,levIndexBot,headerIndex)
                end if
                call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex,anlIncValue)
              else if (bufrCode == bufr_nefs) then
                 anlIncUwind=col_getElem(columnAnlInc,col_getNumLev(columnAnlInc,varLevel),headerIndex,'UU') 
                 anlIncVwind=col_getElem(columnAnlInc,col_getNumLev(columnAnlInc,varLevel),headerIndex,'VV') 
                 trlUwind=col_getElem(columnTrlOnAnlIncLev,col_getNumLev(columnTrlOnAnlIncLev,varLevel),headerIndex,'UU') 
                 trlVwind=col_getElem(columnTrlOnAnlIncLev,col_getNumLev(columnTrlOnAnlIncLev,varLevel),headerIndex,'VV') 
                 squareSum=trlUwind**2+trlVwind**2
                 if ( squareSum .gt. 1.d-10 ) then 
                   trlWindSpeed=sqrt(squareSum)
                   anlIncWindSpeed=( trlUwind * anlIncUwind + trlVwind * anlIncVwind ) / trlWindSpeed
                 else
                   trlWindSpeed=0.0
                   anlIncWindSpeed=0.0 
                 end if
                 call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex,anlIncWindSpeed)  
              else if (bufrCode == bufr_neps .OR. bufrCode == bufr_nepn) THEN
                ZLTV  = columnTrlOnAnlIncLev%OLTV(1,nlev_T,headerIndex)*col_getElem(columnAnlInc,nlev_T,headerIndex,'TT')  & 
                      + columnTrlOnAnlIncLev%OLTV(2,nlev_T,headerIndex)*col_getElem(columnAnlInc,nlev_T,headerIndex,'HU')
                trlVirtTemp  = columnTrlOnAnlIncLev%OLTV(1,nlev_T,headerIndex)*col_getElem(columnTrlOnAnlIncLev,nlev_T,headerIndex,'TT')
                deltaT= temperatureLapseRate*(obsHeight-col_getHeight(columnTrlOnAnlIncLev,0,headerIndex,'SF'))
                coeffA  = ((trlVirtTemp-deltaT)/trlVirtTemp)
                ZDELPS= (col_getElem(columnAnlInc,1,headerIndex,'P0')*coeffA**coeffB)
                ZDELTV= ((col_getElem(columnTrlOnAnlIncLev,1,headerIndex,'P0')*coeffB*coeffA**(coeffB-1))  &
                     *(deltaT/(trlVirtTemp*trlVirtTemp)*ZLTV))
                call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex, ZDELPS+ZDELTV)
              else
                
                call utl_abort('oop_Hsf: You have entered the twilight zone!')

              end if

            end if

         end do BODY

      end do FAMILY

    end subroutine oop_Hsf

    !--------------------------------------------------------------------------

    subroutine oop_Hsst()
      ! :Purpose: Compute simulated sea surface temperature observations 
      !           from profiled model increments.
      !           It returns Hdx in OBS_WORK
      implicit none

      integer :: headerIndex, bodyIndex, bufrCode
      real(8) :: anlIncValueBot
      character(len=4) :: varName

      call obs_set_current_body_list( obsSpaceData, 'TM' )

      !$OMP PARALLEL DO PRIVATE(bodyIndex, bufrCode, varName, headerIndex, anlIncValueBot)
      BODY: do bodyIndex = 1, obs_numBody( obsSpaceData )
        if ( obs_getFamily(obsSpaceData, bodyIndex_opt=bodyIndex) /= 'TM' ) cycle BODY

        bufrCode = obs_bodyElem_i( obsSpaceData, OBS_VNM, bodyIndex )
        if ( bufrCode /= bufr_sst ) cycle BODY

        if ( col_varExist(columnAnlInc,'TM') ) then
          varName = 'TM'
        else
          varName = 'TG'
        end if

        if ( obs_bodyElem_i( obsSpaceData, OBS_ASS, bodyIndex ) == obs_assimilated ) then

          headerIndex = obs_bodyElem_i( obsSpaceData, OBS_HIND, bodyIndex )
          anlIncValueBot = col_getElem(columnAnlInc, 1, headerIndex, varName_opt = varName )
          call obs_bodySet_r( obsSpaceData, OBS_WORK, bodyIndex, anlIncValueBot )
        end if

      end do BODY
      !$OMP END PARALLEL DO

    end subroutine oop_Hsst

    !--------------------------------------------------------------------------

    subroutine oop_Hhydro()
      ! :Purpose: Compute simulated hydrological observations 
      !           from profiled model increments.
      !           It returns Hdx in OBS_WORK
      implicit none

      integer :: headerIndex, bodyIndex, bufrCode
      real(8) :: anlIncValueBot
      character(len=4) :: varName

      call obs_set_current_body_list( obsSpaceData, 'HY' )

      BODY: do
        bodyIndex = obs_getBodyIndex( obsSpaceData )
        if (bodyIndex < 0) exit BODY

        ! Process all data within the domain of the model
        bufrCode = obs_bodyElem_i( obsSpaceData, OBS_VNM, bodyIndex )

        if ( bufrCode /= bufr_riverFlow ) cycle BODY

        if ( obs_bodyElem_i( obsSpaceData, OBS_ASS, bodyIndex ) == obs_assimilated ) then
          headerIndex = obs_bodyElem_i( obsSpaceData, OBS_HIND, bodyIndex )
          varName     = vnl_varNameFromVarNum(bufrCode)
          anlIncValueBot  = col_getElem(columnAnlInc, 1, headerIndex,varName_opt = varName)
          call obs_bodySet_r( obsSpaceData, OBS_WORK, bodyIndex, anlIncValueBot )
        end if

      end do BODY

    end subroutine oop_Hhydro

    !--------------------------------------------------------------------------

    subroutine oop_Hice()
      ! :Purpose: Compute simulated sea ice concentration observations 
      !           from profiled model increments.
      !           It returns Hdx in OBS_WORK
      implicit none

      integer          :: headerIndex, bodyIndex, bufrCode
      integer          :: idate, imonth
      integer          :: trackCellNum
      real(8)          :: anlIncValue, scaling
      character(len=4) :: varName

      call obs_set_current_body_list( obsSpaceData, 'GL' )

      BODY: do

        bodyIndex = obs_getBodyIndex( obsSpaceData )
        if (bodyIndex < 0) exit BODY

        ! Process all data within the domain of the model
        scaling = oop_iceScaling(obsSpaceData, bodyIndex)

        if (scaling == 0.0d0) cycle BODY

        if ( obs_bodyElem_i( obsSpaceData, OBS_ASS, bodyIndex ) == obs_assimilated ) then
          bufrCode = obs_bodyElem_i( obsSpaceData, OBS_VNM, bodyIndex )
          varName = vnl_varNameFromVarNum(bufrCode)
          headerIndex = obs_bodyElem_i( obsSpaceData, OBS_HIND, bodyIndex )
          anlIncValue = scaling*col_getElem( columnAnlInc, 1, headerIndex, varName_opt = varName )
          call obs_bodySet_r( obsSpaceData, OBS_WORK, bodyIndex, anlIncValue )
        end if

      end do BODY

    end subroutine oop_Hice

    !--------------------------------------------------------------------------

    subroutine oop_Hto()
      ! :Purpose: Compute simulated radiances observations from profiled model
      !          increments.
      !          It returns Hdx in OBS_WORK
      implicit none

      integer :: datestamp

      if (.not.obs_famExist(obsSpaceData,'TO', localMPI_opt = .true. )) return

      !     1.   Prepare atmospheric profiles for all tovs observation points for use in rttov
      !     .    -----------------------------------------------------------------------------
      !
      if (min_nsim == 1) then
        datestamp = tim_getDatestamp()
        call tvs_fillProfiles(columnTrlOnAnlIncLev, obsSpaceData, datestamp, "tlad", .false.)
      end if


      !     2.   Compute radiance
      !     .    ----------------
      !
      call tvslin_rttov_tl(columnAnlInc, columnTrlOnAnlIncLev, obsSpaceData)


    end subroutine oop_Hto


    subroutine oop_HheightCoordObs()
      !
      ! :Purpose: Compute simulated geometric-height based observations
      !           such as Radar Doppler velocity.
      !           It returns Hdx in OBS_WORK
      implicit none

      ! locals
      integer :: bodyIndex, headerIndex, levelIndex, bufrCode, layerIndex, familyIndex
      integer :: bodyIndexStart, bodyIndexEnd, bodyIndex1
      real(8) :: levelAltLow, levelAltHigh, HDX, Azimuth, obsAltitude
      real(8) :: uuLow, uuHigh, vvLow, vvHigh, vInterpWeightHigh, vInterpWeightLow
      real(8) :: anlIncValueLow, anlIncValueHigh, trlValueLow, trlValueHigh
      real(8), pointer :: du_column(:), dv_column(:), height_column(:)
      integer, parameter :: NUMFAMILY=3
      character(len=2) :: listFamily(NUMFAMILY), cfam

      listFamily(1) = 'RA' ! Doppler velocity (Radial Wind) burf_radvel
      listFamily(2) = 'PR' ! Dew point difference           burf_nees
      listFamily(3) = 'AL' ! Aladin HLOS wind               burf_neal
      ! Loop over all family
      FAMILY: do familyIndex=1, NUMFAMILY

        ! Loop over all header indices 
        call obs_set_current_header_list(obsSpaceData, listFamily(familyIndex))
        HEADER: do  
          headerIndex = obs_getHeaderIndex(obsSpaceData)
          if ( headerIndex < 0 ) exit HEADER

          if ( listFamily(familyIndex) == 'RA' ) then

            ! Azimuth of the radar beam
            azimuth   = obs_headElem_r(obsSpaceData, OBS_RZAM, headerIndex)

          end if

          ! Local vector state
          du_column  => col_getColumn(columnAnlInc, headerIndex, 'UU') 
          dv_column  => col_getColumn(columnAnlInc, headerIndex, 'VV')

          ! Loop over all body indices 
          call obs_set_current_body_list(obsSpaceData, headerIndex)
          BODY: do
            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if ( bodyIndex < 0 ) exit BODY

            ! only process observations flagged to be assimilated
            if ( obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) /= obs_assimilated ) cycle BODY

            ! Check that this observation has the expected bufr element ID
            bufrCode = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex)
            if ( bufrCode == BUFR_logRadarPrecip ) cycle BODY

            if  ( bufrCode == bufr_radvel ) then
              ! Operator Hx for tangential wind is linear in x
              ! that is:
              !    h(x) = Hx
              ! and
              ! h(xb) + Hdx  =  Hxb + Hdx 
              !
              ! H includes vertical interpolation  
              !   and projection of U and V wind components along the direction of the beam
              !
              ! In matrix form for one Doppler velocity observation:
              !   VDoppler = Hx = [ iwHigh*sin(az) iwHigh*cos(az) iwLow*sin(az) iwLow*cos(az) ][ uuHigh ]
              !                                                                                [ vvHigh ]
              !                                                                                [ uuLow  ] 
              !                                                                                [ vvLow  ] 
              ! such that
              !   VDoppler =  iwHigh*sin(az)*uuHigh 
              !             + iwHigh*cos(az)*vvHigh
              !             + iwLow *sin(az)*uuLow 
              !             + iwLow *cos(az)*vvLow
              !
              ! With
              !   az     = beam azimuth (radians, met convention) 
              !   iwHigh = Interpolation Weight for model level just above observation
              !   iwLow  =                 "                         below     "
              !   uuHigh and vvHigh = wind components on model level just above observation
              !   uuLow  and vvLow  =                 "                   below     "
              !
              ! The dependence of model levels on surface pressure is neglected here

              ! OBS_LYR returns the index of the model level just above the observation
              !   level=1 is the highest level such that level+1 is lower 
              levelIndex = obs_bodyElem_i(obsSpaceData, OBS_LYR, bodyIndex)
              levelAltHigh = col_getHeight(columnTrlOnAnlIncLev, levelIndex,   headerIndex, 'MM')
              levelAltLow  = col_getHeight(columnTrlOnAnlIncLev, levelIndex+1, headerIndex, 'MM')
              obsAltitude = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex) 
           
              !vertical interpolation weights 
              vInterpWeightHigh = (obsAltitude - levelAltLow)/(levelAltHigh - levelAltLow)
              vInterpWeightLow = 1.0D0 - vInterpWeightHigh

              HDX =  vInterpWeightHigh*sin(azimuth)*du_column(levelIndex)   &
                   + vInterpWeightHigh*cos(azimuth)*dv_column(levelIndex)   &
                   + vInterpWeightLow *sin(azimuth)*du_column(levelIndex+1) &
                   + vInterpWeightLow *cos(azimuth)*dv_column(levelIndex+1)
                
              ! Store HDX in OBS_WORK  
              call obs_bodySet_r(obsSpaceData, OBS_WORK, bodyIndex, HDX)
              
            else 

              cfam = obs_getfamily(obsSpaceData,headerIndex_opt=headerIndex)
              write(*,*) 'CANNOT ASSIMILATE OBSERVATION!!!', &
                         'bufrCode =', bufrCode, 'cfam =',  cfam
              call utl_abort('oop_HheightCoordObs')

            end if
          end do BODY
        end do HEADER
      end do FAMILY

    end subroutine oop_HheightCoordObs

    !--------------------------------------------------------------------------

    subroutine oop_Hro( initializeLinearization_opt )
      !
      ! :Purpose: Compute the tangent operator for GPSRO observations.
      !
      implicit none
      logical, optional :: initializeLinearization_opt 

      real(8) :: ZMHXL
      real(8) :: DX(gps_ncvmx)

      integer :: IDATYP
      integer :: JL, JV, NGPSLEV, NWNDLEV
      integer :: headerIndex, bodyIndex, iProfile

      logical :: ASSIM

      integer :: NH, NH1

      ! Initializations
      NGPSLEV = col_getNumLev(columnAnlInc,'TH')
      NWNDLEV = col_getNumLev(columnAnlInc,'MM')

      ! call to calculate the GPSRO Jacobians
      call oop_calcGPSROJacobian( columnTrlOnAnlIncLev, obsSpaceData, &
                                  initializeLinearization_opt=initializeLinearization_opt )

      ! Loop over all header indices of the 'RO' family (Radio Occultation)
      ! Set the header list (start at the beginning of the list)
      call obs_set_current_header_list(obsSpaceData,'RO')
      HEADER: do
         headerIndex = obs_getHeaderIndex(obsSpaceData)
         if (headerIndex < 0) exit HEADER

         ! Process only refractivity data (codtyp 169)
         IDATYP = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
         DATYP: if ( IDATYP == 169 ) THEN

            ! Scan for requested data values of the profile, and count them
            ASSIM = .FALSE.
            NH = 0

            ! Loop over all body indices for this headerIndex:
            ! (start at the beginning of the list)
            call obs_set_current_body_list(obsSpaceData, headerIndex)
            BODY: do 
               bodyIndex = obs_getBodyIndex(obsSpaceData)
               if (bodyIndex < 0) exit BODY
               if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated ) THEN
                  ASSIM = .TRUE.
                  NH = NH + 1
               end if
            end do BODY

            ! If assimilations are requested, prepare and apply the observation operator
            ASSIMILATE: if (ASSIM) THEN
               iProfile=gps_iprofile_from_index(headerIndex)

               ! Local vector state
               do JL = 1, NGPSLEV
                  DX (        JL) = col_getElem(columnAnlInc,JL,headerIndex,'TT')
                  DX (NGPSLEV+JL) = col_getElem(columnAnlInc,JL,headerIndex,'HU')
               end do
               DX (2*NGPSLEV+1:3*NGPSLEV) = col_getColumn(columnAnlInc,headerIndex,'Z_T')
               DX (3*NGPSLEV+1:4*NGPSLEV) = col_getColumn(columnAnlInc,headerIndex,'P_T')

               ! Perform the (H(xb)DX-Y') operation
               ! Loop over all body indices for this headerIndex:
               NH1 = 0
               call obs_set_current_body_list(obsSpaceData, headerIndex)
               BODY_3: do 
                  bodyIndex = obs_getBodyIndex(obsSpaceData)
                  if (bodyIndex < 0) exit BODY_3
                  if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated ) THEN
                     NH1 = NH1 + 1

                     ! Evaluate H(xb)DX
                     ZMHXL = 0.d0
                     do JV = 1, 4*NGPSLEV
                        ZMHXL = ZMHXL + oop_vRO_Jacobian(iProfile,NH1,JV) * DX(JV)
                     end do

                     ! Store in CMA
                     call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex, ZMHXL)
                  end if
               end do BODY_3
            end if ASSIMILATE
         end if DATYP
      end do HEADER

      
    end subroutine oop_Hro

    !--------------------------------------------------------------------------

    subroutine oop_Hgp( initializeLinearization_opt )
      !
      ! :Purpose: Compute H'dx for all GPS ZTD observations
      !           oop_Hgp TL of DOBSGPSGB (Jo for GB-GPS ZTD observations)
      implicit none
      logical, optional :: initializeLinearization_opt 

      real(8) :: ZHX
      real(8) :: DX(gps_ncvmx)

      integer :: headerIndex, bodyIndex
      integer :: JL, NFLEV, status, iztd, icount

      logical :: ASSIM
      character(len=12) :: cstnid

      NFLEV  = col_getNumLev(columnTrlOnAnlIncLev,'TH')

      icount = 0

      ! call to calculate the GPSGB Jacobians
      call oop_calcGPSGBJacobian( columnTrlOnAnlIncLev, obsSpaceData, &
                                  initializeLinearization_opt=initializeLinearization_opt )

      ! loop over all header indices of the 'GP' family (GPS observations)
      call obs_set_current_header_list(obsSpaceData,'GP')
      HEADER: do
         headerIndex = obs_getHeaderIndex(obsSpaceData)
         if (headerIndex < 0) exit HEADER

         ! Scan for ZTD assimilation at this location
         ASSIM = .FALSE.
         ! loop over all body indices for this headerIndex
         call obs_set_current_body_list(obsSpaceData, headerIndex)
         BODY: do 
            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if (bodyIndex < 0) exit BODY

            if (   (obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex) == 189) &
                 .and. (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) == bufr_nezd) &
                 .and. (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) ) then
               ASSIM = .TRUE.
            end if
         end do BODY

         ! If ZTD assimilation, apply the TL observation operator
         if ( ASSIM ) THEN
            iztd = gps_iztd_from_index(headerIndex)
            if ( iztd < 1 .or. iztd > gps_gb_numZTD ) then
               call utl_abort('oop_Hgp: ERROR: index from gps_iztd_from_index() is out of range!')
            end if

            ! Local vector state (analysis increments)
            do JL = 1, NFLEV
               DX (JL)        = col_getElem(columnAnlInc,JL,headerIndex,'TT')
               DX (NFLEV+JL)  = col_getElem(columnAnlInc,JL,headerIndex,'HU')
            end do
            DX (2*NFLEV+1:3*NFLEV) = col_getColumn(columnAnlInc,headerIndex,'Z_T')
            DX (3*NFLEV+1:4*NFLEV) = col_getColumn(columnAnlInc,headerIndex,'P_T')

            ! Evaluate H'(xb)*dX
            ZHX = 0.D0
            do JL = 1, 4*NFLEV
               ZHX = ZHX + oop_vZTD_Jacobian(iztd,JL)*DX(JL)
            end do

            ! Store ZHX = H'dx in OBS_WORK
            ! loop over all body indices for this headerIndex
            call obs_set_current_body_list(obsSpaceData, headerIndex)
            BODY_2: do 
               bodyIndex = obs_getBodyIndex(obsSpaceData)
               if (bodyIndex < 0) exit BODY_2

               if (   (obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex) == 189) &
                    .and. (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) == bufr_nezd) &
                    .and. (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) ) then
                  call obs_bodySet_r(obsSpaceData,OBS_WORK,bodyIndex, ZHX)
                  icount = icount + 1
                  if ( icount <= 3 .and. gps_gb_LTESTOP ) then
                    cstnid = obs_elem_c(obsSpaceData,'STID',headerIndex)
                    write(*,*) iztd, cstnid
                    write(*,*) 'JAC(ncv) = ', (oop_vZTD_Jacobian(iztd,JL),JL=1,4*NFLEV)
                    write(*,*) 'DTT(JL)  = ', (DX(JL),JL=1,NFLEV)
                    write(*,*) 'DHU(JL)  = ', (DX(JL),JL=NFLEV+1,2*NFLEV)
                    write(*,*) 'DAL(JL)  = ', (DX(JL),JL=2*NFLEV+1,3*NFLEV)
                    write(*,*) 'DP (JL)  = ', (DX(JL),JL=3*NFLEV+1,4*NFLEV)
                    write(*,*) 'ZHX (mm) = ', ZHX*1000.D0
                  end if
               end if
            end do BODY_2

         end if ! ASSIM

      end do HEADER

      !      WRITE(*,*) 'oop_Hgp: Number of ZTD data locations with obs_bodySet_r(OBS_OMA) = ', icount

      !      WRITE(*,*)'EXIT oop_Hgp'

      
    end subroutine oop_Hgp

    !--------------------------------------------------------------------------

    subroutine oop_Hchm()
      ! :Purpose: Compute simulated chemical constituents observations from profiled model
      !           increments, and returns Hdx in OBS_WORK
      implicit none

      if (.not.obs_famExist(obsSpaceData,'CH', localMPI_opt = .true. )) return
      
      call oopc_CHobsoperators(columnTrlOnAnlIncLev,obsSpaceData,'tl',columnAnlInc_opt=columnAnlInc) ! 'tl' for tangent linear operator

    end subroutine oop_Hchm

  end subroutine oop_Htl

  !--------------------------------------------------------------------------
  ! oop_Had
  !--------------------------------------------------------------------------
  subroutine oop_Had(columnAnlInc, columnTrlOnAnlIncLev, obsSpaceData, initializeLinearization_opt)
    !
    ! :Purpose: Call the several adjoint of observation operators
    implicit none

    type(struct_columnData) :: columnAnlInc
    type(struct_columnData) :: columnTrlOnAnlIncLev
    type(struct_obs)        :: obsSpaceData
    logical, optional       :: initializeLinearization_opt 

    type(struct_vco), pointer :: vco_anl
    logical, save :: initializeLinearization = .true.

    if ( mmpi_myid == 0 ) then
      write(*,*)'OOP_HT- Adjoint of linearized observation operators'
    end if

    vco_anl => col_getVco(columnTrlOnAnlIncLev)

    ! Re-linearlize if it is asked for
    if ( present(initializeLinearization_opt) ) then
      initializeLinearization = initializeLinearization_opt
    end if

    !     Find interpolation layer in model profiles (used by several operators)
    if ( initializeLinearization ) then
      if ( col_getNumLev(columnTrlOnAnlIncLev,'MM') > 1 ) call oop_vobslyrs(columnTrlOnAnlIncLev, obsSpaceData, beSilent=.false.)
      initializeLinearization = .false.
    end if

    call oop_HTchm

    if ( gps_gb_numZTD > 0 ) then
      call oop_HTgp( initializeLinearization_opt=initializeLinearization_opt )
    end if

    call oop_HTheighCoordObs()

    call oop_HTro( initializeLinearization_opt=initializeLinearization_opt )

    call oop_HTto

    call oop_HTsf

    call oop_HTpp

    call oop_HTsst

    call oop_HTice

    call oop_HThydro()


  CONTAINS

    !--------------------------------------------------------------------------

    subroutine oop_HTpp
      ! :Purpose: Adjoint of the "vertical" interpolation, based on vint3d,
      !           for "UPPER AIR" data files.
      implicit none

      integer levIndexBot,levIndexTop,bufrCode
      real(8) ZRES
      real(8) ZWB,ZWT
      real(8) ZLEV,ZPT,ZPB,ZDADPS,ZPRESBPB,ZPRESBPT
      integer headerIndex,layerIndex,nlev_T
      integer bodyIndex,INDEX_FAMILY
      real(8) trlValueTop,trlValueBot
      real(8), pointer :: all_column(:),tt_column(:),hu_column(:),p_column(:)
      real(8) :: delPT,delPB
      integer, parameter :: numFamily=3
      character(len=2) :: list_family(numFamily)
      character(len=4) :: varLevel
      !
      list_family(1) = 'UA'
      list_family(2) = 'AI'
      list_family(3) = 'SW'

      FAMILY: do index_family=1,numFamily

         call obs_set_current_body_list(obsSpaceData,list_family(index_family))
         BODY: do 
            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if (bodyIndex < 0) exit BODY

             ! Process all data within the domain of the model
            if (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated .and. &
                obs_bodyElem_i(obsSpaceData,OBS_XTR,bodyIndex) == 0               .and. &
                obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) == 2 ) then
 
              headerIndex = obs_bodyElem_i(obsSpaceData,OBS_HIND,bodyIndex)
               ZRES = obs_bodyElem_r(obsSpaceData,OBS_WORK,bodyIndex)
               ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
               bufrCode = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
               varLevel = vnl_varLevelFromVarnum(bufrCode)
               layerIndex   = obs_bodyElem_i(obsSpaceData,OBS_LYR,bodyIndex)
               levIndexTop  = layerIndex  + col_getOffsetFromVarno(columnTrlOnAnlIncLev,bufrCode)
               levIndexBot  = levIndexTop+1
               ZPT  = col_getPressure(columnTrlOnAnlIncLev,layerIndex,headerIndex,varLevel)
               ZPB  = col_getPressure(columnTrlOnAnlIncLev,layerIndex+1,headerIndex,varLevel)
               delPT= col_getPressure(columnAnlInc,layerIndex  ,headerIndex,varLevel)
               delPB= col_getPressure(columnAnlInc,layerIndex+1,headerIndex,varLevel)
               ZWB  = LOG(ZLEV/ZPT)/LOG(ZPB/ZPT)
               ZWT  = 1.0D0 - ZWB

               ZDADPS   = ( LOG(ZLEV/ZPB)*delPT/ZPT -   &
                    LOG(ZLEV/ZPT)*delPB/ZPB )/  &
                    LOG(ZPB/ZPT)**2

               all_column => col_getColumn(columnAnlInc,headerIndex)
               tt_column  => col_getColumn(columnAnlInc,headerIndex,'TT')
               hu_column  => col_getColumn(columnAnlInc,headerIndex,'HU')
               if ( varLevel == 'TH' ) then
                 p_column => col_getColumn(columnAnlInc,headerIndex,'P_T')
               else if ( varLevel == 'MM' ) then
                 p_column => col_getColumn(columnAnlInc,headerIndex,'P_M')
               end if

               if (bufrCode.eq.bufr_nees) then
                  call hutoes_ad(hu_column(layerIndex+1),  &
                       tt_column(layerIndex+1),  &
                       p_column(layerIndex+1),   &
                       ZWB*ZRES,         &
                       col_getElem(columnTrlOnAnlIncLev,layerIndex+1,headerIndex,'HU'),      &
                       col_getPressure(columnTrlOnAnlIncLev,layerIndex+1,headerIndex,'TH'))
                  call hutoes_ad(hu_column(layerIndex  ),  &
                       tt_column(layerIndex  ),  &
                       p_column(layerIndex),     &
                       ZWT*ZRES,         &
                       col_getElem(columnTrlOnAnlIncLev,layerIndex  ,headerIndex,'HU'),      &
                       col_getPressure(columnTrlOnAnlIncLev,layerIndex  ,headerIndex,'TH'))
                  trlValueBot=hutoes(col_getElem(columnTrlOnAnlIncLev,layerIndex+1,headerIndex,'HU'),  &
                       col_getElem(columnTrlOnAnlIncLev,layerIndex+1,headerIndex,'TT'),  &
                       col_getPressure(columnTrlOnAnlIncLev,layerIndex+1,headerIndex,'TH'))
                  trlValueTop=hutoes(col_getElem(columnTrlOnAnlIncLev,layerIndex  ,headerIndex,'HU'),  &
                       col_getElem(columnTrlOnAnlIncLev,layerIndex  ,headerIndex,'TT'),  &
                       col_getPressure(columnTrlOnAnlIncLev,layerIndex  ,headerIndex,'TH'))
               else
                  all_column(levIndexBot) = all_column(levIndexBot) + ZWB*ZRES
                  all_column(levIndexTop) = all_column(levIndexTop) + ZWT*ZRES
                  trlValueBot=col_getElem(columnTrlOnAnlIncLev,levIndexBot,headerIndex)
                  trlValueTop=col_getElem(columnTrlOnAnlIncLev,levIndexTop,headerIndex)
               end if
               p_column(layerIndex  ) = p_column(layerIndex  ) + (trlValueBot - trlValueTop) * &
                                (LOG(ZLEV/ZPB)/ZPT/LOG(ZPB/ZPT)**2) * ZRES
               p_column(layerIndex+1) = p_column(layerIndex+1) - (trlValueBot - trlValueTop) * &
                                (LOG(ZLEV/ZPT)/ZPB/LOG(ZPB/ZPT)**2) * ZRES
            end if

         end do BODY

      end do FAMILY

    end subroutine oop_HTpp

    !--------------------------------------------------------------------------

    subroutine oop_HTsf
      ! :Purpose: based on surfc1dz to build the adjoint of the
      !          vertical interpolation for SURFACE data files.
      implicit none

      integer :: levIndexBot,levIndexTop
      real(8) :: ZRES
      real(8) :: ZWB, ZWT,coeffA,coeffB,ZATV,trlVirtTemp
      real(8) :: obsHeight, ZPT, ZPB, ZDADPS, ZDELPS, ZDELTV, deltaT
      real(8) :: trlValueBot
      real(8) :: trlUwind, trlVwind, sumSquare, trlWindSpeed, anlIncWindSpeed 
      integer :: headerIndex, layerIndex, nlev, nlev_T
      integer :: bodyIndex, bufrCode, INDEX_FAMILY
      real(8), pointer :: all_column(:), tt_column(:), hu_column(:), ps_column(:), p_column(:)
      real(8), pointer :: du_column(:), dv_column(:)
      real(8) :: dPdPsfc
      integer, parameter :: numFamily=5
      character(len=2) :: list_family(numFamily)
      character(len=4) :: varLevel

      !- Temperature lapse rate for extrapolation of height below model surface
      coeffB   = 1.0d0/(MPC_RGAS_DRY_AIR_R8*temperatureLapseRate/ec_rg)

      !
      !-   1. Fill in COMMVO by using the adjoint of the "vertical" interpolation
      !
      list_family(1) = 'UA'
      list_family(2) = 'SF'
      list_family(3) = 'SC'
      list_family(4) = 'GP'
      list_family(5) = 'RA'

      FAMILY: do index_family=1,numFamily

         call obs_set_current_body_list(obsSpaceData, list_family(index_family))
         BODY: do
            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if (bodyIndex < 0) exit BODY

            ! Process all data within the domain of the model
            bufrCode = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
            if ( bufrCode == bufr_nezd .or. bufrCode == bufr_radvel ) cycle BODY
            if (    (obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) == 1) &
                 .and. (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) &
                 .and. (bufrCode == bufr_nets .or. bufrCode == bufr_neps  &
                 .or. bufrCode == bufr_nepn .or. bufrCode == bufr_ness  &
                 .or. bufrCode == bufr_neus .or. bufrCode == bufr_nevs  &
                 .or. bufrCode == bufr_vis  .or. bufrCode == bufr_logVis  &
                 .or. bufrCode == bufr_gust .or. bufrCode == bufr_nefs &
                 .or. bufrCode == bufr_radarPrecip .or. bufrCode == bufr_logRadarPrecip  &
                 .or. obs_bodyElem_i(obsSpaceData,OBS_XTR,bodyIndex) == 0) ) then

               varLevel    = vnl_varLevelFromVarnum(bufrCode)
               nlev        = col_getNumLev(columnAnlInc,varLevel)
               nlev_T      = col_getNumLev(columnAnlInc,'TH')
               headerIndex = obs_bodyElem_i(obsSpaceData,OBS_HIND,bodyIndex)
               bufrCode    = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
               layerIndex  = obs_bodyElem_i(obsSpaceData,OBS_LYR,bodyIndex)
               obsHeight   = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
               levIndexTop = nlev - 1 + col_getOffsetFromVarno(columnTrlOnAnlIncLev,bufrCode)
               levIndexBot = levIndexTop+1
               ZRES        = obs_bodyElem_r(obsSpaceData,OBS_WORK,bodyIndex)
               
               if (bufrCode == bufr_nets .or. bufrCode == bufr_ness .or.  &
                   bufrCode == bufr_neus .or. bufrCode == bufr_nevs .or. & 
                   bufrCode == bufr_vis  .or. bufrCode == bufr_logVis  .or.  &
                   bufrCode == bufr_gust .or.  &
                   bufrCode == bufr_radarPrecip .or. bufrCode == bufr_logRadarPrecip) then
                 if ( bufrCode == bufr_ness ) then
                   tt_column  => col_getColumn(columnAnlInc,headerIndex,'TT')
                   hu_column  => col_getColumn(columnAnlInc,headerIndex,'HU')
                   p_column   => col_getColumn(columnAnlInc,headerIndex,'P_T')
                   call hutoes_ad(hu_column(nlev_T),  &
                        tt_column(nlev_T),  &
                        p_column(nlev_T),   &
                        ZRES,               &
                        col_getElem(columnTrlOnAnlIncLev,nlev_T,headerIndex,'HU'),  &
                        col_getPressure(columnTrlOnAnlIncLev,nlev_T,headerIndex,'TH'))
                 else
                   all_column => col_getColumn(columnAnlInc,headerIndex) 
                   all_column(levIndexBot) = all_column(levIndexBot) + ZRES
                 end if
               else if ( bufrCode == bufr_nefs ) then
                   du_column  => col_getColumn(columnAnlInc,headerIndex,'UU') 
                   dv_column  => col_getColumn(columnAnlInc,headerIndex,'VV')
                   anlIncWindSpeed=obs_bodyElem_r(obsSpaceData,OBS_WORK,bodyIndex) 
                   sumSquare = col_getElem(columnTrlOnAnlIncLev,nlev,headerIndex,'UU')**2 + &
                               col_getElem(columnTrlOnAnlIncLev,nlev,headerIndex,'VV')**2 
                   if ( sumSquare > 1.0d-10 ) then
                     trlWindSpeed=sqrt(sumSquare)
                     du_column(nlev) = du_column(nlev) + &
                          anlIncWindSpeed*col_getElem(columnTrlOnAnlIncLev,nlev,headerIndex,'UU')/trlWindSpeed
                     dv_column(nlev) = dv_column(nlev) + &
                          anlIncWindSpeed*col_getElem(columnTrlOnAnlIncLev,nlev,headerIndex,'VV')/trlWindSpeed
                   else
                     trlWindSpeed=0.
                   end if
               else if ( bufrCode == bufr_neps .or. bufrCode == bufr_nepn ) then
                 tt_column  => col_getColumn(columnAnlInc,headerIndex,'TT')
                 hu_column  => col_getColumn(columnAnlInc,headerIndex,'HU')
                 ps_column  => col_getColumn(columnAnlInc,headerIndex,'P0')
                 trlVirtTemp = columnTrlOnAnlIncLev%OLTV(1,nlev_T,headerIndex)*col_getElem(columnTrlOnAnlIncLev,nlev_T,headerIndex,'TT')
                 deltaT = temperatureLapseRate*(obsHeight-col_getHeight(columnTrlOnAnlIncLev,0,headerIndex,'SF'))
                 coeffA = ((trlVirtTemp-deltaT)/trlVirtTemp)
                 ZDELTV = (col_getElem(columnTrlOnAnlIncLev,1,headerIndex,'P0')*coeffB*coeffA**(coeffB-1))  &
                      *(deltaT/(trlVirtTemp*trlVirtTemp))
                 ZDELPS = coeffA**coeffB
                 ZATV   = ZDELTV*ZRES
                 ps_column(1)= ps_column(1) + ZDELPS*ZRES
                 tt_column(nlev_T) = tt_column(nlev_T) + &
                                     columnTrlOnAnlIncLev%OLTV(1,nlev_T,headerIndex)*ZATV
                 hu_column(nlev_T) = hu_column(nlev_T) + &
                                     columnTrlOnAnlIncLev%OLTV(2,nlev_T,headerIndex)*ZATV
               else

                 call utl_abort('oop_HTsf: You have entered the twilight zone!')

               end if

             end if

         end do BODY

      end do FAMILY

    end subroutine oop_HTsf

    !--------------------------------------------------------------------------

    subroutine oop_HTsst
      ! :Purpose: Adjoint of the "vertical" interpolation for SST data
      implicit none

      real(8) :: residual
      integer :: headerIndex, bodyIndex, bufrCode 
      real(8), pointer :: columnTG(:)
      character(len=4) :: varName

      call obs_set_current_body_list( obsSpaceData, 'TM' )

      !$OMP PARALLEL DO PRIVATE(bodyIndex, bufrCode, varName, headerIndex, residual, columnTG)
      BODY: do bodyIndex = 1, obs_numBody( obsSpaceData )
        if ( obs_getFamily(obsSpaceData, bodyIndex_opt=bodyIndex) /= 'TM' ) cycle BODY

        bufrCode = obs_bodyElem_i( obsSpaceData, OBS_VNM, bodyIndex )
        if ( bufrCode /= bufr_sst ) cycle BODY

        if ( col_varExist(columnAnlInc,'TM') ) then
          varName = 'TM'
        else
          varName = 'TG'
        end if

        if ( obs_bodyElem_i( obsSpaceData, OBS_ASS, bodyIndex ) == obs_assimilated ) then

          headerIndex = obs_bodyElem_i( obsSpaceData, OBS_HIND, bodyIndex )
          residual = obs_bodyElem_r( obsSpaceData, OBS_WORK, bodyIndex )
          columnTG => col_getColumn(columnAnlInc, headerIndex, varName_opt = varName ) 
          columnTG(1) = columnTG(1) + residual

        end if

      end do BODY
      !$OMP END PARALLEL DO

    end subroutine oop_HTsst

    !--------------------------------------------------------------------------

    subroutine oop_HThydro
      ! :Purpose: Adjoint of the "vertical" interpolation for Hydrological data
      implicit none

      real(8) :: residual
      integer :: headerIndex, bodyIndex, bufrCode 
      real(8), pointer :: columnHY(:)
      character(len=4) :: varName

      call obs_set_current_body_list( obsSpaceData, 'HY' )

      BODY: do

        bodyIndex = obs_getBodyIndex( obsSpaceData )
        if (bodyIndex < 0) exit BODY

        ! Process all data within the domain of the model
        bufrCode = obs_bodyElem_i( obsSpaceData, OBS_VNM, bodyIndex )

        if ( bufrCode /= bufr_riverFlow ) cycle BODY

        if ( obs_bodyElem_i( obsSpaceData, OBS_ASS, bodyIndex ) == obs_assimilated ) then

          headerIndex = obs_bodyElem_i( obsSpaceData, OBS_HIND, bodyIndex )
          residual = obs_bodyElem_r( obsSpaceData, OBS_WORK, bodyIndex )
          varName = vnl_varNameFromVarNum(bufrCode)
          columnHY => col_getColumn(columnAnlInc, headerIndex, varName_opt = varName ) 
          columnHY(1) = columnHY(1) + residual
        end if

      end do BODY

    end subroutine oop_HThydro

    !--------------------------------------------------------------------------

    subroutine oop_HTice
      ! :Purpose: Adjoint of the "vertical" interpolation for ICE data
      implicit none

      real(8)          :: residual, scaling
      integer          :: headerIndex, bodyIndex, bufrCode
      real(8), pointer :: columnGL(:)
      character(len=4) :: varName

      call obs_set_current_body_list( obsSpaceData, 'GL' )

      BODY: do

        bodyIndex = obs_getBodyIndex( obsSpaceData )
        if (bodyIndex < 0) exit BODY

        ! Process all data within the domain of the model

        scaling = oop_iceScaling(obsSpaceData, bodyIndex)

        if (scaling == 0.0d0) cycle BODY

        if ( obs_bodyElem_i( obsSpaceData, OBS_ASS, bodyIndex ) == obs_assimilated ) then
          residual = scaling*obs_bodyElem_r( obsSpaceData, OBS_WORK, bodyIndex )
          bufrCode = obs_bodyElem_i( obsSpaceData, OBS_VNM, bodyIndex )
          varName = vnl_varNameFromVarNum(bufrCode)
          headerIndex = obs_bodyElem_i( obsSpaceData, OBS_HIND, bodyIndex )
          columnGL => col_getColumn( columnAnlInc, headerIndex, varName_opt = varName )
          columnGL(1) = columnGL(1) + residual
        end if

      end do BODY

    end subroutine oop_HTice

    !--------------------------------------------------------------------------

    subroutine oop_HTto
      ! :Purpose: Adjoint of computation of residuals to the tovs observations
      implicit none

      integer :: datestamp

      if (.not.obs_famExist(obsSpaceData,'TO', localMPI_opt = .true. )) return

      ! Adjoint of computing radiance
      datestamp = tim_getDatestamp()
      call tvs_fillProfiles(columnTrlOnAnlIncLev,obsSpaceData,datestamp,"tlad",.false.)

      call tvslin_rttov_ad(columnAnlInc,columnTrlOnAnlIncLev,obsSpaceData)

    end subroutine oop_HTto

    !--------------------------------------------------------------------------

    subroutine oop_HTro ( initializeLinearization_opt ) 
      !
      ! :Purpose: Compute the adjoint operator for GPSRO observations.
      implicit none
      logical, optional :: initializeLinearization_opt 

      real(8) :: DPJO0(gps_ncvmx)
      real(8) :: DPJO1(gps_ncvmx)

      real(8) :: ZINC

      real(8), pointer :: tt_column(:),hu_column(:),height_column(:),p_column(:)
      integer :: IDATYP
      integer :: JL, NGPSLEV
      integer :: headerIndex, bodyIndex, iProfile

      logical :: ASSIM, LUSE

      integer :: NH, NH1

      ! Initializations
      NGPSLEV=col_getNumLev(columnAnlInc,'TH')

      ! call to calculate the GPSRO Jacobians
      call oop_calcGPSROJacobian( columnTrlOnAnlIncLev, obsSpaceData, &
                                  initializeLinearization_opt=initializeLinearization_opt )

      ! Loop over all header indices of the 'RO' family (Radio Occultation)
      ! Set the header list (start at the beginning of the list)
      call obs_set_current_header_list(obsSpaceData,'RO')
      HEADER: do
         headerIndex = obs_getHeaderIndex(obsSpaceData)
         if (headerIndex < 0) exit HEADER

         DPJO0 = 0.d0

         ! Process only refractivity data (codtyp 169)
         IDATYP = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
         DATYP: if ( IDATYP == 169 ) THEN

            ! Scan for requested data values of the profile, and count them
            ASSIM = .FALSE.
            NH = 0

            ! Loop over all body indices for this headerIndex:
            ! (start at the beginning of the list)
            call obs_set_current_body_list(obsSpaceData, headerIndex)
            BODY: do 
               bodyIndex = obs_getBodyIndex(obsSpaceData)
               if (bodyIndex < 0) exit BODY

               LUSE=( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated )
               if ( LUSE ) THEN
                  ASSIM = .TRUE.
                  NH = NH + 1
               end if
            end do BODY

            ! If assimilations are requested, prepare and apply the observation operator
            ASSIMILATE: if (ASSIM) THEN
               iProfile=gps_iprofile_from_index(headerIndex)

               ! Perform the (H(xb)DX-Y')/S operation
               NH1 = 0

               ! Loop over all body indices for this headerIndex:
               ! (start at the beginning of the list)
               call obs_set_current_body_list(obsSpaceData, headerIndex)
               BODY_3: do 
                  bodyIndex = obs_getBodyIndex(obsSpaceData)
                  if (bodyIndex < 0) exit BODY_3

                  LUSE=( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated )
                  if ( LUSE ) THEN
                     NH1 = NH1 + 1

                     ! Normalized increment
                     ZINC = obs_bodyElem_r(obsSpaceData,OBS_WORK,bodyIndex)

                     ! O-F Tested criteria:
                     DPJO1(1:4*NGPSLEV) = ZINC * oop_vRO_Jacobian(iProfile,NH1,1:4*NGPSLEV)

                     ! Accumulate the gradient of the observation cost function:
                     DPJO0(1:4*NGPSLEV) = DPJO0(1:4*NGPSLEV) + DPJO1(1:4*NGPSLEV)
                  end if
               end do BODY_3
            end if ASSIMILATE
         end if DATYP

         ! Store H* (HX - Z)/SIGMA in COMMVO
         tt_column => col_getColumn(columnAnlInc,headerIndex,'TT')
         hu_column => col_getColumn(columnAnlInc,headerIndex,'HU')
         height_column => col_getColumn(columnAnlInc,headerIndex,'Z_T')
         p_column => col_getColumn(columnAnlInc,headerIndex,'P_T')
         do JL = 1, NGPSLEV
            tt_column(JL) = DPJO0(JL)
            hu_column(JL) = DPJO0(JL+NGPSLEV)
            height_column(JL) = DPJO0(JL+2*NGPSLEV)
            p_column(JL)  = DPJO0(JL+3*NGPSLEV)
         end do
      end do HEADER
      
    end subroutine oop_HTro


    !--------------------------------------------------------------------------

    subroutine oop_HTgp( initializeLinearization_opt ) 
      !
      ! :Purpose: Compute Ht*grad(Jo) for all GPS ZTD observations
      !
      ! :Note:  ZTD Jacobians are computed and stored in oop_Hgp (first iter.)
      implicit none
      logical, optional :: initializeLinearization_opt 

      real(8) :: DPJO0(gps_ncvmx)
      real(8) :: JAC(gps_ncvmx)
      ! 
      real(8) :: ZINC
      integer :: JL, NFLEV, iztd
      integer :: headerIndex, bodyIndex, icount
      logical :: ASSIM

      real(8), pointer :: tt_column(:),hu_column(:),height_column(:),p_column(:)

      !      WRITE(*,*)'ENTER oop_HTgp'

      NFLEV  = col_getNumLev(columnTrlOnAnlIncLev,'TH')

      ! call to calculate the GPSGB Jacobians
      call oop_calcGPSGBJacobian( columnTrlOnAnlIncLev, obsSpaceData, &
                                  initializeLinearization_opt=initializeLinearization_opt )

      ! loop over all header indices of the 'GP' family (GPS observations)
      ! Set the header list & start at the beginning of the list
      call obs_set_current_header_list(obsSpaceData,'GP')

      icount = 0

      HEADER: do
         headerIndex = obs_getHeaderIndex(obsSpaceData)
         if (headerIndex < 0) exit HEADER

         DPJO0(:) = 0.0D0
         JAC(:)   = 0.0D0

         ! Scan for requested ZTD assimilation
         ASSIM = .FALSE.

         ! loop over all body indices (still in the 'GP' family)
         ! Set the body list & start at the beginning of the list)
         call obs_set_current_body_list(obsSpaceData, headerIndex)
         BODY: do 
            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if (bodyIndex < 0) exit BODY

            if (   (obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex) == 189) &
                 .and. (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) == bufr_nezd) &
                 .and. (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) ) then
               ASSIM = .TRUE.
            end if
         end do BODY
         !
         if (ASSIM) THEN

            icount = icount + 1
            iztd = gps_iztd_from_index(headerIndex)
            if ( iztd < 1 .or. iztd > gps_gb_numZTD ) then
               call utl_abort('oop_HTgp: ERROR: index from gps_iztd_from_index() is out of range!')
            end if

            do JL = 1, 4*NFLEV
               JAC(JL) = oop_vZTD_Jacobian(iztd,JL)
            end do

            ! Get Ht*grad(HeaderIndex) = Ht*(H'dx - d)/sigma_o^2
            ! loop over all body indices (still in the 'GP' family)
            ! Start at the beginning of the list)
            call obs_set_current_body_list(obsSpaceData, headerIndex)
            BODY_2: do 
               bodyIndex = obs_getBodyIndex(obsSpaceData)
               if (bodyIndex < 0) exit BODY_2
               if (   (obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex ) == 189) &
                    .and. (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) == bufr_nezd) &
                    .and. (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) ) then
                  ZINC = obs_bodyElem_r(obsSpaceData,OBS_WORK,bodyIndex)

                  ! Accumulate the gradient of the observation cost function
                  DPJO0(1:4*NFLEV) = ZINC * oop_vZTD_Jacobian(iztd,:)
               end if
            end do BODY_2

            ! Store Ht*grad(HeaderIndex) in COMMVO
            tt_column     => col_getColumn(columnAnlInc,headerIndex,'TT')
            hu_column     => col_getColumn(columnAnlInc,headerIndex,'HU')
            height_column => col_getColumn(columnAnlInc,headerIndex,'Z_T')
            p_column      => col_getColumn(columnAnlInc,headerIndex,'P_T')
            do JL = 1, NFLEV
               tt_column(JL) = DPJO0(JL)
               hu_column(JL) = DPJO0(JL+NFLEV)
               height_column(JL) = DPJO0(JL+2*NFLEV)
               p_column (JL) = DPJO0(JL+3*NFLEV)
            end do

         end if ! ASSIM

      end do HEADER

      !      WRITE(*,*) 'oop_HTgp: Number of ZTD data locations processed = ', icount

      !      WRITE(*,*)'EXIT oop_HTgp'

      
    end subroutine oop_HTgp

    ! --------------------------------------------------------------------------

    subroutine oop_HTheighCoordObs()
      ! :Purpose: Compute the adjoint operator of the
      !           vertical interpolation of geometric-height based data.
      !           Including Radar Doppler velocity data.
      !
      !  delta x is updated by this routine 
      implicit none

      ! locals
      integer :: bodyIndex, headerIndex, levelIndex, bufrCode, familyIndex
      integer :: bodyIndexStart, bodyIndexEnd, bodyIndex1
      real(8) :: levelAltLow, levelAltHigh, azimuth, obsAltitude, HDX 
      real(8) :: anlIncValueLow, anlIncValueHigh, trlValueLow, trlValueHigh
      real(8) :: vInterpWeightLow, vInterpWeightHigh
      real(8), pointer :: du_column(:), dv_column(:), height_column(:)
      integer, parameter :: NUMFAMILY=3
      character(len=2) :: listFamily(NUMFAMILY), cfam
     
      listFamily(1) = 'RA' ! Doppler velocity (Radial Wind) burf_radvel
      listFamily(2) = 'PR' ! Dew point difference           burf_nees
      listFamily(3) = 'AL' ! Aladin HLOS wind               burf_neal
      
      ! Loop over all families
      FAMILY: do familyIndex = 1, NUMFAMILY

        ! Loop over header indices 
        call obs_set_current_header_list(obsSpaceData, listFamily(familyIndex))
        HEADER: do  
          headerIndex = obs_getHeaderIndex(obsSpaceData)  
          if ( headerIndex < 0 ) exit HEADER

          if  ( listFamily(familyIndex) == 'RA' ) then
            ! Azimuth of the radar beam
            azimuth   = obs_headElem_r(obsSpaceData, OBS_RZAM, headerIndex)
          end if

          ! Local vector state
          du_column  => col_getColumn(columnAnlInc, headerIndex, 'UU') 
          dv_column  => col_getColumn(columnAnlInc, headerIndex, 'VV')

          ! Loop over body indices 
          call obs_set_current_body_list(obsSpaceData, headerIndex)
          BODY: do
            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if ( bodyIndex < 0 ) exit BODY
            ! only process observations flagged to be assimilated
            if ( obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) /= obs_assimilated ) cycle BODY

            ! Check that this observation has the expected bufr element ID
            bufrCode = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex)
            if ( bufrCode == BUFR_logRadarPrecip ) cycle BODY

            if ( bufrCode == bufr_radvel ) then
              ! OBS_LYR returns the index of the model level just above the observation
              !   level=1 is the highest level such that level+1 is lower 
              levelIndex = obs_bodyElem_i(obsSpaceData, OBS_LYR, bodyIndex)
              levelAltHigh = col_getHeight(columnTrlOnAnlIncLev, levelIndex,   headerIndex, 'MM')
              levelAltLow  = col_getHeight(columnTrlOnAnlIncLev, levelIndex+1, headerIndex, 'MM')
              obsAltitude = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex) 

              ! HDX from store
              HDX = obs_bodyElem_r(obsSpaceData, OBS_WORK, bodyIndex)
              vInterpWeightHigh = (obsAltitude - levelAltLow)/(levelAltHigh - levelAltLow)
              vInterpWeightLow = 1.0D0 - vInterpWeightHigh

              ! see oop_HheightCoordObs for the description of the H operator for radar Doppler velocity
              ! 
              ! In matrix form, adjoint of H operator for one observation:
              !
              !  delta x         = Ht TangentialWindResidual
              !
              !  [ deltaUUHigh ] = [ iwHigh*sin(az) ][ residual ]  
              !  [ deltaVVigh  ]   [ iwHigh*cos(az) ]
              !  [ deltaUULow  ]   [ iwLow *sin(az) ]
              !  [ deltaVVLow  ]   [ iwLow *cos(az) ]
              !
              ! With
              !   az     = beam azimuth (radians, met convention) 
              !   iwHigh = Interpolation Weight for model level just above observation
              !   iwLow  =                 "                         below     "
              !   residual = OmP - Hdx
              !   deltaUUHigh and deltaVVHigh = wind components on model level just above observation
              !   deltaUULow  and deltaVVLow  =                 "                   below     "

              !delta x is updated
              du_column(levelIndex  ) = du_column(levelIndex  ) + vInterpWeightHigh*HDX*sin(azimuth)
              dv_column(levelIndex  ) = dv_column(levelIndex  ) + vInterpWeightHigh*HDX*cos(azimuth)
              du_column(levelIndex+1) = du_column(levelIndex+1) + vInterpWeightLow *HDX*sin(azimuth)
              dv_column(levelIndex+1) = dv_column(levelIndex+1) + vInterpWeightLow *HDX*cos(azimuth)

            else
              
              cfam = obs_getfamily(obsSpaceData,headerIndex_opt=headerIndex)
              write(*,*) 'CANNOT ASSIMILATE OBSERVATION!!!', &
                         'bufrCode =', bufrCode, 'cfam =',  cfam
              call utl_abort('oop_HTheighCoordObs')

            end if
          end do BODY
        end do HEADER
      end do FAMILY

    end subroutine oop_HTheighCoordObs

    !--------------------------------------------------------------------------

    subroutine oop_HTchm
      ! :Purpose: Compute H^T * R^-1 (OmP-Hdx) for all CH observations  
      implicit none
      
      if (.not.obs_famExist(obsSpaceData,'CH', localMPI_opt = .true. )) return
      
      call oopc_CHobsoperators(columnTrlOnAnlIncLev,obsSpaceData,'adjoint', & 
                               columnAnlInc_opt=columnAnlInc) ! 'adjoint' for adjoint of the tangent linear operator

    end subroutine oop_HTchm

  end subroutine oop_Had

  !--------------------------------------------------------------------------
  ! HUtoES
  !--------------------------------------------------------------------------
  function HUtoES( hu, tt, pressure ) result( es )
    ! :Purpose:
    !          to calculate the dew point depression from specific
    !          humidity, temperature and pressure.  No ice phase
    !          is permitted and the pressure vector is given.
    implicit none

    real(8), intent(in) :: hu
    real(8), intent(in) :: tt
    real(8), intent(in) :: pressure

    real(8)             :: es

    real(8) :: husat, td

    ! get the saturated vapor pressure from specific humidity
    husat = phf_foefq8(hu,pressure)

    ! now the dewpoint temperature
    td = phf_fotw8(husat)

    ! finally the dewpoint depression
    es = min(tt-td,MPC_MAXIMUM_ES_R8)

  end function HUtoES

  !--------------------------------------------------------------------------
  ! HUtoES_tl
  !--------------------------------------------------------------------------
  function HUtoES_tl(HU_inc,TT_inc,P_inc,HU_trl,PRES_trl) result(ES_inc)
    ! :Purpose: TLM VERSION
    !           to calculate the dew point depression from specific
    !           humidity, temperature and pressure.  No ice phase
    !           is permitted and the pressure vector is given.
    implicit none

    real(8)             :: ES_inc

    ! Arguments:
    real(8), intent(in) :: HU_inc
    real(8), intent(in) :: TT_inc
    real(8), intent(in) :: P_inc
    real(8), intent(in) :: HU_trl
    real(8), intent(in) :: PRES_trl

    ! Locals:
    real(8) :: ZE, ZTD, dTDdE, ZQBRANCH
    real(8) :: dESdLQ, dESdTT, dESdP

    dESdTT = 1.0d0

    !- Forward calculations of saturation vapour pressure and dewpoint temperature
    !  and adjoint of vapour pressure from adjoint of dewpoint temperature
    ZE   = phf_FOEFQ8(HU_trl, PRES_trl)
    ZTD  = phf_FOTW8 (ZE)
    dTDdE= phf_FODTW8(ZTD,ZE)

    !- adjoint of temp. specific humidity and surface pressure due to changes in vapour pressure
    ZQBRANCH = phf_FQBRANCH(HU_trl)

    dESdLQ = - ZQBRANCH*phf_FOEFQA(1.0d0,dTDdE,HU_trl,PRES_trl)

    dESdP  = - ZQBRANCH*phf_FOEFQPSA(1.0d0,dTDdE,HU_trl,1.0d0)-  &
               (1.D0-ZQBRANCH)*(dTDdE*1.0d0)

    ES_inc =  dESdLQ*HU_inc/HU_trl + dESdP*P_inc + dESdTT*TT_inc

  end function HUtoES_tl

  !--------------------------------------------------------------------------
  ! HUtoES_ad
  !--------------------------------------------------------------------------
  subroutine HUtoES_ad(HU_inc,TT_inc,P_inc,ES_inc,HU_trl,PRES_trl)
    ! Purpose: ADJOINT VERSION
    !          to calculate the dew point depression from specific
    !          humidity, temperature and pressure.  No ice phase
    !          is permitted and the pressure vector is given.
    implicit none

    real(8), intent(inout) :: HU_inc,TT_inc,P_inc
    real(8), intent(in)  :: ES_inc,HU_trl,PRES_trl
    real(8) :: ZE,ZTD,dTDdE,ZQBRANCH
    real(8) :: dESdLQ,dESdTT,dESdP

    dESdTT = 1.0d0
   
    !- Forward calculations of saturation vapour pressure and dewpoint temperature
    !  and adjoint of vapour pressure from adjoint of dewpoint temperature
    ZE = phf_FOEFQ8(HU_trl, PRES_trl)

    ZTD=phf_FOTW8(ZE)
    dTDdE=phf_FODTW8(ZTD,ZE)

    !- adjoint of temp. specific humidity and surface pressure due to changes in vapour pressure
    ZQBRANCH = phf_FQBRANCH(HU_trl)
    dESdLQ = - ZQBRANCH*phf_FOEFQA(1.0d0,dTDdE,HU_trl,PRES_trl)

    dESdP  = - ZQBRANCH*phf_FOEFQPSA(1.0d0,dTDdE,HU_trl,1.0d0)-  &
               (1.D0-ZQBRANCH)*(dTDdE*1.0d0)

    ! TLM: ES_inc =  dESdLQ*HU_inc/HU_trl + dESdP*P_inc + dESdTT*TT_inc
    ! ADJOINT:
    HU_inc = HU_inc + dESdLQ*ES_inc/HU_trl
    P_inc  = P_inc  + dESdP *ES_inc
    TT_inc = TT_inc + dESdTT*ES_inc

  end subroutine HUtoES_ad

  !--------------------------------------------------------------------------
  ! oop_calcGPSROJacobian
  !--------------------------------------------------------------------------
  subroutine oop_calcGPSROJacobian( columnTrlOnAnlIncLev, obsSpaceData, initializeLinearization_opt )
    !
    ! :Purpose: Calculating the Jacobians of refractivity for oop_Hro/oop_HTro
    implicit none

    type(struct_columnData) :: columnTrlOnAnlIncLev
    type(struct_obs)        :: obsSpaceData
    logical, optional       :: initializeLinearization_opt 

    real(8) :: zlat, lat
    real(8) :: zlon, lon
    real(8) :: zazm, azm
    integer :: isat
    real(8) :: rad, geo, zp0
    real(8), allocatable :: zpp(:), ztt(:), zhu(:), zHeight(:), zuu(:), zvv(:)
    real(8) :: zmt
    integer :: IDATYP, varNum
    integer :: jl, ngpslev, nwndlev
    integer :: headerIndex, bodyIndex, iProfile
    logical :: ASSIM
    logical, save :: initializeLinearization = .true.
    integer :: nh, nh1
    type(gps_profile)           :: prf
    real(8)       , allocatable :: h   (:),azmv(:)
    type(gps_diff), allocatable :: rstv(:)

    ! Re-compute the Jacobian for re-linearized state
    if ( present(initializeLinearization_opt) ) then
      initializeLinearization = initializeLinearization_opt
    end if

    if ( .not. initializeLinearization ) return

    write(*,*) 'ENTER oop_calcGPSROJacobian'

    initializeLinearization = .false.

    ! Initializations
    ngpslev = col_getNumLev(columnTrlOnAnlIncLev,'TH')
    nwndlev = col_getNumLev(columnTrlOnAnlIncLev,'MM')

    allocate(zpp (ngpslev))
    allocate(ztt (ngpslev))
    allocate(zhu (ngpslev))
    allocate(zHeight (ngpslev))
    allocate(zuu (ngpslev))
    allocate(zvv (ngpslev))

    if ( .not. allocated(oop_vRO_Jacobian) ) then
      allocate( oop_vRO_Jacobian(gps_numroprofiles,gps_ro_maxprfsize,4*ngpslev) )
      allocate( oop_vRO_lJac    (gps_numROProfiles) )
      oop_vRO_Jacobian = 0.d0
      oop_vRO_lJac = .False.
    end if

    allocate( h    (gps_ro_maxprfsize) )
    allocate( azmv (gps_ro_maxprfsize) )
    allocate( rstv (gps_ro_maxprfsize) )

    ! Loop over all header indices of the 'RO' family (Radio Occultation)
    ! Set the header list (start at the beginning of the list)
    call obs_set_current_header_list(obsSpaceData,'RO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER

      ! Process only refractivity data (codtyp 169)
      IDATYP = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
      DATYP: if ( IDATYP == 169 ) then
        ! Scan for requested data values of the profile, and count them
        assim = .false.
        nh = 0

        ! Loop over all body indices for this headerIndex:
        ! (start at the beginning of the list)
        call obs_set_current_body_list(obsSpaceData, headerIndex)
        BODY: do 
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY
          if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated ) then
            assim = .true.
            nh = nh + 1
          endif
        enddo BODY

        ! If assimilations are requested, prepare and apply the observation operator
        ASSIMILATE: if (assim) then
          iProfile = gps_iprofile_from_index(headerIndex)
          if (oop_vRO_lJac(iProfile)) cycle                  ! If already done, end this HEADER
          varNum = gps_vRO_IndexPrf(iProfile, 2)

          ! Profile at the observation location:
          ! Basic geometric variables of the profile:
          zlat = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)
          zlon = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex)
          isat = obs_headElem_i(obsSpaceData,OBS_SAT,headerIndex)
          rad  = obs_headElem_r(obsSpaceData,OBS_TRAD,headerIndex)
          geo  = obs_headElem_r(obsSpaceData,OBS_GEOI,headerIndex)
          azm = obs_headElem_r(obsSpaceData,OBS_AZA,headerIndex)
          zmt  = col_getHeight(columnTrlOnAnlIncLev,0,headerIndex,'SF')
          lat  = zlat * MPC_DEGREES_PER_RADIAN_R8
          lon  = zlon * MPC_DEGREES_PER_RADIAN_R8
          zazm = azm / MPC_DEGREES_PER_RADIAN_R8
          zp0  = col_getElem(columnTrlOnAnlIncLev,1,headerIndex,'P0')
          do jl = 1, ngpslev
            ! Profile x_b
            zpp(jl) = col_getPressure(columnTrlOnAnlIncLev,jl,headerIndex,'TH')
            ztt(jl) = col_getElem(columnTrlOnAnlIncLev,jl,headerIndex,'TT') - MPC_K_C_DEGREE_OFFSET_R8
            zhu(jl) = col_getElem(columnTrlOnAnlIncLev,jl,headerIndex,'HU')
            zHeight(jl) = col_getHeight(columnTrlOnAnlIncLev,jl,headerIndex,'TH')
          enddo

          if ((col_getPressure(columnTrlOnAnlIncLev,1,headerIndex,'TH') + 1.0d-4) <  &
               col_getPressure(columnTrlOnAnlIncLev,1,headerIndex,'MM')) then
            ! case with top thermo level above top momentum level (Vcode=5002)
            do jl = 1, nwndlev
              zuu(jl) = col_getElem(columnTrlOnAnlIncLev,jl,headerIndex,'UU')
              zvv(jl) = col_getElem(columnTrlOnAnlIncLev,jl,headerIndex,'VV')
            enddo
          else
            ! case without top thermo above top momentum level or unstaggered (Vcode=5001/4/5)
            do jl = 1, nwndlev-1
              zuu(jl) = col_getElem(columnTrlOnAnlIncLev,jl+1,headerIndex,'UU')
              zvv(jl) = col_getElem(columnTrlOnAnlIncLev,jl+1,headerIndex,'VV')
            enddo
            zuu(nwndlev) = zuu(nwndlev-1)
            zvv(nwndlev) = zuu(nwndlev-1)
          endif
          zuu(ngpslev) = zuu(nwndlev)
          zvv(ngpslev) = zuu(nwndlev)

          ! GPS profile structure:
          call gps_struct1sw_v2(ngpslev,zlat,zlon,zazm,zmt,rad,geo,zp0,zpp,ztt,zhu,zuu,zvv,zHeight,prf)

          ! Prepare the vector of all the observations:
          nh1 = 0
          call obs_set_current_body_list(obsSpaceData, headerIndex)
          BODY_2: do 
            bodyIndex = obs_getBodyIndex(obsSpaceData)
            if (bodyIndex < 0) exit BODY_2
            if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated ) then
              nh1      = nh1 + 1
              h(nh1)   = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
              azmv(nh1)= zazm
            endif
          enddo BODY_2

          ! Apply the observation operator:
          ! varNum = bufr_nebd (15037) or varNum = bufr_nerf (15036) for GPS-RO
          if (varNum == bufr_nebd) then
            call gps_bndopv1(h, azmv, nh, prf, rstv)
          else
            call gps_refopv (h, nh, prf, rstv)
          end if
          do nh1 = 1, nh
            oop_vRO_Jacobian(iprofile,nh1,1:4*ngpslev)= rstv(nh1)%dvar(1:4*ngpslev)
          end do
          oop_vRO_lJac(iProfile) = .True.
        endif ASSIMILATE
      endif DATYP
    enddo HEADER

    deallocate( rstv )
    deallocate( azmv )
    deallocate( h    )

    deallocate(zvv)
    deallocate(zuu)
    deallocate(zHeight)
    deallocate(zhu)
    deallocate(ztt)
    deallocate(zpp)

    write(*,*) 'EXIT oop_calcGPSROJacobian'
    

  end subroutine oop_calcGPSROJacobian

  !--------------------------------------------------------------------------
  ! oop_calcGPSGBJacobian
  !--------------------------------------------------------------------------
  subroutine oop_calcGPSGBJacobian( columnTrlOnAnlIncLev, obsSpaceData, initializeLinearization_opt )
    !
    ! :Purpose: Calculating the Jacobians of ZTD for oop_Hgp/oop_HTgp
    implicit none

    type(struct_columnData) :: columnTrlOnAnlIncLev
    type(struct_obs)        :: obsSpaceData
    logical, optional       :: initializeLinearization_opt 

    real(8) :: ZLAT, Lat
    real(8) :: ZLON, Lon
    real(8), allocatable :: ZTTB(:)
    real(8), allocatable :: ZHUB(:)
    real(8), allocatable :: zHeight(:)
    real(8), allocatable :: ZPPB(:)

    real(8) :: ZP0B, ZPSMOD, ZPWMOD, ZPWMOD2, dZTD
    real(8) :: ZMT
    real(8) :: sfcfield
    real(8) :: dxq1, dxq2, dxq3

    real(8) :: ZLEV, ZDZMIN
    real(8) :: JAC(gps_ncvmx)

    integer :: headerIndex, bodyIndex
    integer :: JL, NFLEV, iztd, icount, vcode

    logical :: ASSIM

    type(gps_profilezd) :: PRF, PRF2
    type(gps_diff)      :: ZTDOPV, ZTDOPV2

    type(struct_vco), pointer :: vco_anl
    character(len=12) :: cstnid

    logical, save :: initializeLinearization = .true.

    ! Re-compute the Jacobian for re-linearized state
    if ( present(initializeLinearization_opt) ) then
      initializeLinearization = initializeLinearization_opt
    end if

    if ( .not. initializeLinearization ) return

    write(*,*) 'ENTER oop_calcGPSGBJacobian'

    initializeLinearization = .FALSE.

    vco_anl => col_getVco(columnTrlOnAnlIncLev)
    vcode = vco_anl%vCode

    ZDZMIN = gps_gb_DZMIN                     ! from modgpsztd_mod

    NFLEV  = col_getNumLev(columnTrlOnAnlIncLev,'TH')

    ! Initializations
    if ( .not. allocated(gps_ZTD_Index) ) call utl_abort('oop_calcGPSGBJacobian: ERROR:  gps_ZTD_Index not allocated!')
    if ( allocated(oop_vZTD_Jacobian) ) write(*,*) 'oop_calcGPSGBJacobian: WARNING, oop_vZTD_Jacobian is already allocated. Re-allocating'
    call utl_reallocate(oop_vZTD_Jacobian,gps_gb_numZTD,4*NFLEV)
    oop_vZTD_Jacobian(:,:) = 0.0d0

    allocate(ZTTB(NFLEV))
    allocate(ZHUB(NFLEV))
    allocate(zHeight(NFLEV))
    allocate(ZPPB(NFLEV))

    write(*,*) 'oop_calcGPSGBJacobian: Storing Jacobians for GPS ZTD data ...'
    write(*,*) '   INFO: Analysis grid iversion = ', vcode
    write(*,*) '         col_getNumLev[columnTrlOnAnlIncLev,TH] = ', NFLEV
    write(*,*) '         gps_gb_numZTD = ', gps_gb_numZTD

    icount = 0

    ! loop over all header indices of the 'GP' family (GPS observations)
    call obs_set_current_header_list(obsSpaceData,'GP')
    HEADER_0: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER_0

      ! Scan for ZTD assimilation at this location
      ASSIM = .FALSE.

      ! loop over all body indices for this headerIndex
      call obs_set_current_body_list(obsSpaceData, headerIndex)
      BODY_0: do 
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if (bodyIndex < 0) exit BODY_0
        if (   (obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex) == 189) &
             .and. (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) == bufr_nezd) &
             .and. (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated) ) then
          ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
          ASSIM = .TRUE.
        end if
      end do BODY_0

      if ( ASSIM ) THEN
        ! LR background profile at the observation location x :
        icount = icount + 1
        Lat  = obs_headElem_r(obsSpaceData,OBS_LAT,headerIndex)
        ZLAT = Lat * MPC_DEGREES_PER_RADIAN_R8
        Lon  = obs_headElem_r(obsSpaceData,OBS_LON,headerIndex)
        ZLON = Lon * MPC_DEGREES_PER_RADIAN_R8
        ZP0B = col_getElem(columnTrlOnAnlIncLev,1,headerIndex,'P0')
        do JL = 1, NFLEV
          ZTTB(JL) = col_getElem(columnTrlOnAnlIncLev,JL,headerIndex,'TT') - MPC_K_C_DEGREE_OFFSET_R8
          ZHUB(JL) = col_getElem(columnTrlOnAnlIncLev,JL,headerIndex,'HU')
          ZPPB(JL) = col_getPressure(columnTrlOnAnlIncLev,JL,headerIndex,'TH')
          zHeight(JL) = col_getHeight(columnTrlOnAnlIncLev,JL,headerIndex,'TH')
        end do
        if ( abs(ZPPB(NFLEV)-ZP0B) > 0.1 ) then
          write(*,*) ' oop_calcGPSGBJacobian: ERROR: |ZPPB(NFLEV)-ZP0B| > 0.1'
          write(*,*) '          ZPPB(NFLEV), ZP0B =', ZPPB(NFLEV), ZP0B
        end if
        ZMT = col_getHeight(columnTrlOnAnlIncLev,0,headerIndex,'SF')

        call gps_structztd_v2(NFLEV,Lat,Lon,ZMT,ZP0B,ZPPB,ZTTB,ZHUB,zHeight,gps_gb_LBEVIS,gps_gb_IREFOPT,PRF)
        call gps_ztdopv(ZLEV,PRF,gps_gb_LBEVIS,ZDZMIN,ZTDopv,ZPSMOD,gps_gb_IZTDOP)

        ! Observation Jacobian H'(xb)            
        JAC = ZTDopv%DVar
        iztd = gps_iztd_from_index(headerIndex)
        do JL = 1, 4*NFLEV
           oop_vZTD_Jacobian(iztd,JL) = JAC(JL)
        end do

        if ( icount <= 3 .and. gps_gb_LTESTOP ) then
          write(*,*) '--------------------------------------------------------- '
          cstnid = obs_elem_c(obsSpaceData,'STID',headerIndex)
          write(*,*) iztd, cstnid, 'ZTDopv (m) = ', ZTDopv%Var
          call gps_pw(PRF,ZPWMOD)

          sfcfield = ZP0B + 50.0d0
          call gps_structztd_v2(NFLEV,Lat,Lon,ZMT,sfcfield,ZPPB,ZTTB,ZHUB,zHeight,gps_gb_LBEVIS,gps_gb_IREFOPT,PRF2)
          call gps_ztdopv(ZLEV,PRF2,gps_gb_LBEVIS,ZDZMIN,ZTDopv2,ZPSMOD,gps_gb_IZTDOP)
          write(*,*) ' ZTD Operator Test:  dP0 = +50 Pa'
          write(*,*) ' dZTD NL     = ', ZTDopv2%Var - ZTDopv%Var
          write(*,*) ' dZTD Linear = ', oop_vZTD_Jacobian(iztd,4*NFLEV)*50.0d0
          write(*,*) ' '

          ! q dx 
          dxq1 = 0.44D-01*ZHUB(64)
          dxq2 = 0.44D-01*ZHUB(65)
          dxq3 = 0.44D-01*ZHUB(66)
          ZHUB(64) = ZHUB(64) - dxq1
          ZHUB(65) = ZHUB(65) - dxq2
          ZHUB(66) = ZHUB(66) - dxq3
          call gps_structztd_v2(NFLEV,Lat,Lon,ZMT,ZP0B,ZPPB,ZTTB,ZHUB,zHeight,gps_gb_LBEVIS,gps_gb_IREFOPT,PRF2)
          call gps_ztdopv(ZLEV,PRF2,gps_gb_LBEVIS,ZDZMIN,ZTDopv2,ZPSMOD,gps_gb_IZTDOP)
          call gps_pw(PRF2,ZPWMOD2)
          write(*,*) ' ZTD Operator Test:  dQ = -0.44E-01*Q JL = 64,65,66'
          write(*,*) ' dPW (mm)    = ', ZPWMOD2 - ZPWMOD
          write(*,*) ' dZTD NL     = ', ZTDopv2%Var - ZTDopv%Var
          dZTD = oop_vZTD_Jacobian(iztd,64+NFLEV)*(-dxq1) + oop_vZTD_Jacobian(iztd,65+NFLEV)*(-dxq2) + &
               oop_vZTD_Jacobian(iztd,66+NFLEV)*(-dxq3)
          write(*,*) ' dZTD Linear = ', dZTD
          write(*,*) ' '
          ZHUB(64) = ZHUB(64) + dxq1
          ZHUB(65) = ZHUB(65) + dxq2
          ZHUB(66) = ZHUB(66) + dxq3

          ! temperature dx                   
          ZTTB(64) = ZTTB(64) + 2.0d0
          ZTTB(65) = ZTTB(65) + 2.0d0
          ZTTB(66) = ZTTB(66) + 2.0d0
          call gps_structztd_v2(NFLEV,Lat,Lon,ZMT,ZP0B,ZPPB,ZTTB,ZHUB,zHeight,gps_gb_LBEVIS,gps_gb_IREFOPT,PRF2)
          call gps_ztdopv(ZLEV,PRF2,gps_gb_LBEVIS,ZDZMIN,ZTDopv2,ZPSMOD,gps_gb_IZTDOP)
          write(*,*) ' ZTD Operator Test:  dTT = +2.0K JL = 64,65,66'
          write(*,*) ' dZTD NL     = ', ZTDopv2%Var - ZTDopv%Var
          dZTD = oop_vZTD_Jacobian(iztd,64)*2.0d0 + oop_vZTD_Jacobian(iztd,65)*2.0d0 + &
               oop_vZTD_Jacobian(iztd,66)*2.0d0
          write(*,*) ' dZTD Linear = ', dZTD
          write(*,*) '--------------------------------------------------------- '              
        end if

      end if

    end do HEADER_0

    deallocate(ZTTB)
    deallocate(ZHUB)
    deallocate(zHeight)
    deallocate(ZPPB)

    write(*,*) 'oop_calcGPSGBJacobian:   Number of ZTD data (icount) = ', icount
    write(*,*) '           Expected number (gps_gb_numZTD) = ', gps_gb_numZTD
    write(*,*) '           Last iztd                       = ', iztd
    write(*,*) '           gps_ZTD_Index(1)                = ', gps_ZTD_Index(1)
    write(*,*) '           gps_ZTD_Index(iztd)             = ', gps_ZTD_Index(iztd)

    if ( icount /= gps_gb_numZTD ) then
      call utl_abort('oop_calcGPSGBJacobian: ERROR: icount /= gps_gb_numZTD!')
    end if
    if ( icount /= iztd ) then
      call utl_abort('oop_calcGPSGBJacobian: ERROR: icount /= iztd!')
    end if
    if ( gps_gb_numZTD /= iztd ) then
      call utl_abort('oop_calcGPSGBJacobian: ERROR: gps_gb_numZTD /= iztd!')
    end if

    write(*,*) 'EXIT oop_calcGPSGBJacobian'
    
  end subroutine oop_calcGPSGBJacobian

  function oop_iceScaling(obsSpaceData, bodyIndex) result(scaling)
    !
    ! :Purpose: Calculate the scaling factor for
    !           ice related observations to convert from
    !           model space to observation space, i.e.
    !           H(iceConc) = scaling*iceConc + constant
    !
    implicit none
    real(8) :: scaling

    ! Arguments:
    type(struct_obs), intent(in) :: obsSpaceData
    integer,          intent(in) :: bodyIndex

    ! Locals:
    integer          :: bufrCode
    integer          :: headerIndex, obsDate, monthIndex
    integer          :: trackCellNum
    character(len=8) :: ccyymmdd

    bufrCode = obs_bodyElem_i( obsSpaceData, OBS_VNM, bodyIndex )

    select case (bufrCode)
    case(BUFR_ICEC, BUFR_ICEP)
       scaling = 100.0d0
    case(BUFR_ICEV)
       scaling = 1.0d0
    case(BUFR_ICES)
       headerIndex = obs_bodyElem_i( obsSpaceData, OBS_HIND, bodyIndex )
       obsDate = obs_headElem_i( obsSpaceData, OBS_DAT, headerIndex ) 
       write(ccyymmdd, FMT='(i8.8)') obsDate
       read(ccyymmdd(5:6), FMT='(i2)') monthIndex
       trackCellNum = obs_headElem_i( obsSpaceData, OBS_FOV, headerIndex )
       scaling = oer_ascatAnisIce(trackCellNum,monthIndex) - &
            oer_ascatAnisOpenWater(trackCellNum,monthIndex)
    case default
       scaling = 0.0d0
    end select

  end function oop_iceScaling

end module obsOperators_mod
