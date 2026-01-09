
module backgroundCheck_mod
  ! MODULE backgroundCheck_mod (prefix='bgck' category='1. High-level functionality')
  !
  !:Purpose: Performs the background check on all conventional observations
  !
  use mathPhysConstants_mod
  use bufr_mod
  use obsSpaceData_mod
  use gps_mod
  use utilities_mod
  use columnData_mod
  use obsSpaceDiag_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use obsSpaceErrorStdDev_mod
  use obsFamilyList_mod
  use codtyp_mod
  
  implicit none
  private

  ! Public procedures
  public :: bgck_bgcheck_conv

  integer,          parameter :: numFamilyToProcess = 11
  character(len=2), parameter :: familyListToProcess(numFamilyToProcess)= (/ &
                                 'UA','AI','SF','SC','SW','PR','RO','GP','RA', &
                                 'CH','AL' /)

  real(8) :: uvCritDropSonde(3), swCritDropSonde(3), psCritDropSonde(3)
  
  contains

  !--------------------------------------------------------------------------
  ! bgck_bgcheck_conv
  !--------------------------------------------------------------------------
  subroutine bgck_bgcheck_conv(columnTrlOnAnlIncLev, columnTrlOnTrlLev, hco_anl, obsSpaceData)
    !
    !:Purpose: Do background check on all conventional observations
    !
    implicit none

    ! Arguments:
    type(struct_obs)         , intent(inout) :: obsSpaceData         ! Observation-related data
    type(struct_columnData)  , intent(inout) :: columnTrlOnAnlIncLev ! column data on analysis levels
    type(struct_columnData)  , intent(inout) :: columnTrlOnTrlLev    ! column data on trial levels
    type(struct_hco), pointer, intent(in)    :: hco_anl              ! horizontal grid structure

    ! Locals:
    integer                     :: familyIndex
    integer                     :: ier
    logical                     :: noObsToProcess
    character(len=*), parameter :: myName = 'bgck_bgcheck_conv'

    ! namelist variables
    logical                     :: new_bgck_sw ! choose to use the 'new' background check for SW obs
    namelist /NAMBGCKCONV/ new_bgck_sw, uvCritDropSonde, swCritDropSonde, psCritDropSonde

    write(*,*) myName//' begin conventional data background check'

    !- Check if any observations are present for conventional background check
    noObsToProcess = .true.
    do familyIndex = 1, numFamilyToProcess
      if (obs_famExist(obsSpaceData,familyListToProcess(familyIndex))) then
        noObsToProcess = .false.
      end if
    end do
    if (noObsToProcess) then
      write(*,*) myName//': No observations to process'
      return
    end if
     
    call utl_tmg_start(117,'--BgckConventional')

    !- Read namelist variables
    new_bgck_sw = .false.
    uvCritDropSonde(1) = 10.D0
    uvCritDropSonde(2) = 25.D0
    uvCritDropSonde(3) = 30.D0
    swCritDropSonde(1) = 10.D0
    swCritDropSonde(2) = 20.D0
    swCritDropSonde(3) = 30.D0
    psCritDropSonde(1) =  9.D0
    psCritDropSonde(2) = 16.D0
    psCritDropSonde(3) = 25.D0
    
    call utl_tmg_start(181,'low-level--readNML')
    read(utl_flnml, nml = NAMBGCKCONV, IOSTAT = ier)
    if (ier /= 0) write(*,*) myName//': no valid namelist NAMBGCKCONV found, default values will be taken:'
    write(*,*) myName//': new_bgck_sw = ',new_bgck_sw
    call utl_tmg_stop(181)

    !- Obtain or calc OmP-error std dev when requested and possible.
    !  Otherwise calc HBHT contribution (sigma_B in observation space)  
    call ose_computeStddev(columnTrlOnAnlIncLev, hco_anl, & ! IN
                           obsSpaceData)                    ! INOUT

    !- Do a background check of all the observations
    do familyIndex = 1, ofl_numFamily
      ! For SW only, old and new background check schemes controlled by "new_bgck_sw"
      if (obs_famExist(obsSpaceData, ofl_familyList(familyIndex))) then
        call bgck_data(ofl_familyList(familyIndex), obsSpaceData, new_bgck_sw)
      end if
    end do

    if (obs_famExist(obsSpaceData,'RO')) call bgck_gpsro(columnTrlOnTrlLev, obsSpaceData)

    ! Conduct obs-space post-processing diagnostic tasks
    ! (some diagnostic computations controlled by NAMOSD namelist in flnml)
    call osd_ObsSpaceDiag(obsSpaceData, columnTrlOnAnlIncLev, hco_anl, analysisMode_opt=.false.)

    call utl_tmg_stop(117)

  end subroutine bgck_bgcheck_conv

  !--------------------------------------------------------------------------
  ! bgck_data
  !--------------------------------------------------------------------------
  subroutine bgck_data(obsFamily, obsData, new_bgck_sw)
    !
    !:Purpose:  Calculate a background check for a data family and set the
    !           appropriate quality-control flags in obsSpaceData
    implicit none

    ! Arguments:
    character(len=2), intent(in)    :: obsFamily   ! current observation family
    type(struct_obs), intent(inout) :: obsData     ! obsSpaceData
    logical         , intent(in)    :: new_bgck_sw

    ! Locals:
    real(8) :: averageFGE, averageOER
    integer :: obsFlag, burfCode, headerIndex, codeType, bodyIndex, obsCount
    integer :: numberObs, numberObsRejected, INZOBS, INZREJ
    integer :: INPOBS, INTOBS, INDOBS, INPREJ, INTREJ, INDREJ
    real(8) :: oer,omp,fge,omper,bgchk,var,lev,lat,lon,sop
    logical :: LLOK, LLZD
    character(len=12) :: stnid
    integer :: i_ass, i_vco, i_oth, bodyIndex_u, bodyIndex_v, bodyIndex_start
    real(8) :: uu_d, uu_r, uu_f, vv_d, vv_r, vv_f, duv2, duv2_lim, zslev
    logical :: found_u, found_v
    character(len=*), parameter :: myName = 'bgck_data'

    write(*,*)' '
    write(*,*) ' ------------------------------'
    write(*,*) myName//'  background check for', obsFamily, ' data'
    write(*,*) ' ------------------------------'
    write(*,*) ' '
    write(*,'(a55,a74)')'  STNID     LATITU LONGITU  ID Elem        Level        ',  &
    ' Value        Sigmao       Sigmap         O-P       SigmaOP         qcflag  '
    write(*,'(a55,a74)')'  -----     ------ -------  -- ----        -----        ',  &
    ' -----        ------       ------         ---       -------         ------  '

    obsCount = 0
    averageFGE = 0.d0
    averageOER = 0.d0

    numberObs = 0
    numberObsRejected = 0

    ! Initialize counters for GP family observations
    if (obsFamily == 'GP') then
      INZOBS=0
      INPOBS=0
      INTOBS=0
      INDOBS=0
      INZREJ=0
      INPREJ=0
      INTREJ=0
      INDREJ=0
    end if

    if (.not. new_bgck_sw .or. obsFamily /= 'SW') then

      ! loop over all header indices of the current observation family
      call obs_set_current_header_list(obsData, obsFamily)
      HEADER: do
        headerIndex = obs_getHeaderIndex(obsData)
        if (headerIndex < 0) exit HEADER

        stnid = obs_elem_c(obsData, 'STID', headerIndex)

        !- Computation of (y-Hx)**2/(sigmaO**2+sigmaB**2)
          
        ! loop over all body indices for the current headerIndex
        call obs_set_current_body_list(obsData, headerIndex)
        BODY: do 
          bodyIndex = obs_getBodyIndex(obsData)
          if (bodyIndex < 0) exit BODY

          burfCode = obs_bodyElem_i(obsData, OBS_VNM, bodyIndex)
          LLOK = (obs_bodyElem_i(obsData, OBS_ASS, bodyIndex) == obs_assimilated)
          if (LLOK) then
            numberObs = numberObs + 1
            if (obsFamily == 'GP') then
              if (burfCode == BUFR_NEZD) INZOBS = INZOBS+1
              if (burfCode == BUFR_NEPS) INPOBS = INPOBS+1
              if (burfCode == BUFR_NETS) INTOBS = INTOBS+1
              if (burfCode == BUFR_NESS) INDOBS = INDOBS+1
            end if
            var = obs_bodyElem_r(obsData, OBS_VAR, bodyIndex)
            lev = obs_bodyElem_r(obsData, OBS_PPP, bodyIndex)              
            lat = obs_headElem_r(obsData, OBS_LAT, headerIndex) * MPC_DEGREES_PER_RADIAN_R8
            lon = obs_headElem_r(obsData, OBS_LON, headerIndex) * MPC_DEGREES_PER_RADIAN_R8

            ! BACKGROUND CHECK
            omp  = obs_bodyElem_r(obsData, OBS_OMP, bodyIndex)

            ! Get error std dev
              
            ! std(O-P) (valid/available if > 0.0)
            omper = obs_bodyElem_r(obsData, OBS_OMPE, bodyIndex)
                                         
            ! std(O)
            oer = obs_bodyElem_r(obsData, OBS_OER, bodyIndex)
            ! std(P)
            fge = obs_bodyElem_r(obsData, OBS_HPHT, bodyIndex)
            ! NOTE: For GB-GPS ZTD observations (GP family), fge is not the FGE but
            !       rather Std(O-P) set in s/r SETERRGPSGB.
            if (obsFamily == 'GP' .and. burfCode == BUFR_NEZD) then
              if (omper <= 0.0d0 .and. fge > 0.0d0) omper = fge
            end if
              
            ! Observation and background error adjustments
            if (obsFamily == 'GP') then
              ! Protect against error std dev values that are too small
              if (burfCode == BUFR_NEZD) then
                oer = oer / gps_gb_YZDERRWGT
              else
                oer = oer / gps_gb_YSFERRWGT
              end if
              if (oer < 1.0d-3 .and. burfCode /= BUFR_NEZD) then
                write(*,*)' Problem for GP STNID oer= ' , stnid, oer
                call utl_abort(myName//': PROBLEM WITH OER.')
              end if
              IF (fge < 1.0d-3) then
                write(*,*)' Problem for GP STNID FGE= ', stnid, fge
                fge=1.0d-3
              end if
              if (omper > 0.0d0 .and. omper < 1.0d-3) omper = 1.0d-3
            else if (obsFamily == 'CH') then
              if (fge**2 + oer**2 < 1.0d-60) then
                write(*,*) ' Problem for STNID FGE oer=',  stnid, fge, oer
                fge=1.0d30
                oer=1.0d-30
              end if
              if (omper > 0.0 .and. omper < 1.0d-30) omper = 1.0d-30
            else ! default
              if (fge < 0.0d0) fge = 0.0d0
              if (fge**2 + oer**2 < 1.0d-5)then
                write(*,*) ' Problem for STNID FGE oer=', stnid, fge, oer
                fge=1.0d-5
                oer=1.0d-5
              end if
              if (omper > 0.0 .and. omper < 1.0D-5) omper = 1.0d-5
            end if

            ! Calculate bgchk
            if (omper > 0.0d0) then
              bgchk = omp**2 / omper**2 
              sop   = omper
            else
              bgchk = omp**2 / (fge**2+oer**2)
              sop   = sqrt(fge**2+oer**2)
            end if

            !- UPDATE QUALITY CONTROL FLAGS, based on zbgchk
            ! (ELEMENT FLAGS + GLOBAL HEADER FLAGS)
            burfCode = obs_bodyElem_i(obsData, OBS_VNM, bodyIndex)
            if(burfCode == bufr_nees .and. obs_bodyElem_i(obsData, OBS_XTR, bodyIndex) == 0) then
              averageFGE = averageFGE + fge * fge
              averageOER = averageOER + oer * oer
              obsCount   = obsCount + 1
            end if

            ! Set flag level
            codeType = obs_headElem_i(obsData, OBS_ITY, headerIndex)
            obsFlag = setFlag(obsFamily,codeType,burfCode,bgchk)

            ! CONVERT ZTD VALUES FROM M TO MM FOR PRINTOUT
            LLZD = .FALSE.
            if (obsFamily == 'GP' .and. burfCode == BUFR_NEZD) then
              var = var * 1000.0D0
              oer = oer * 1000.0D0
              fge = fge * 1000.0D0
              omp = omp * 1000.0D0
              sop = sop * 1000.0D0
              LLZD = .TRUE.
            end if

            stnid = obs_elem_c(obsData, 'STID', headerIndex)
            if (obsFlag >= 2 .or. (LLZD .and. gps_gb_LTESTOP)) then
              if (obsFamily /= 'CH') then 
                write(*,122)  &
                   stnid, lat, lon, codeType, burfCode, lev, var, oer,  &
                   fge, omp, sop, bgchk, obsFlag
              else 
                write(*,124)  &
                   stnid, lat, lon, codeType, burfCode, lev, var, oer,  &
                   fge, omp, sop, bgchk, obsFlag
              end if
              if (obsFlag >= 2) numberObsRejected = numberObsRejected + 1
              if (obsFlag >= 2 .and. obsFamily == 'GP') then
                if (burfCode == BUFR_NEZD) INZREJ = INZREJ+1
                if (burfCode == BUFR_NEPS) INPREJ = INPREJ+1
                if (burfCode == BUFR_NETS) INTREJ = INTREJ+1
                if (burfCode == BUFR_NESS) INDREJ = INDREJ+1
              end if
            end if

            if (obsFlag == 1) then
              call obs_bodySet_i(obsData, OBS_FLG, bodyIndex, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex), 13))
            else if (obsFlag == 2) then
              call obs_bodySet_i(obsData, OBS_FLG, bodyIndex, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex), 14))
              call obs_bodySet_i(obsData, OBS_FLG, bodyIndex, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex), 16))
              call obs_bodySet_i(obsData, OBS_FLG, bodyIndex, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex), 09))
              call obs_headSet_i(obsData, OBS_ST1, headerIndex, ibset(obs_headElem_i(obsData, OBS_ST1, headerIndex), 06))
              
            else if (obsFlag == 3) then
              call obs_bodySet_i(obsData, OBS_FLG, bodyIndex, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex), 15))
              call obs_bodySet_i(obsData, OBS_FLG, bodyIndex, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex), 16))
              call obs_bodySet_i(obsData, OBS_FLG, bodyIndex, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex), 09))
              call obs_headSet_i(obsData, OBS_ST1, headerIndex, ibset(obs_headElem_i(obsData, OBS_ST1, headerIndex), 06))
            end if
          end if
        end do BODY

124     FORMAT(2x,a9,1x,f6.2,1x,f7.2,1x,I3,1x,I5,7(2x,G11.2),I3)
122     FORMAT(2x,a9,1x,f6.2,1x,f7.2,1x,I3,1x,I5,7(2x,F11.2),I3)

      end do HEADER

    else ! if (.not. new_bgck_sw .or. obsFamily .ne. 'SW') then

      call obs_set_current_body_list(obsData, 'SW')
      bodyuv: do
        bodyIndex = obs_getBodyIndex(obsData)
        if (bodyIndex < 0) exit bodyuv

        ! Only process pressure level observations flagged to be assimilated
        i_ass = obs_bodyElem_i (obsData, OBS_ASS, bodyIndex)
        i_vco = obs_bodyElem_i (obsData, OBS_VCO, bodyIndex)

        if(i_ass /= obs_assimilated .or. i_vco /= 2) cycle bodyuv
 
        headerIndex     = obs_bodyElem_i(obsData, OBS_HIND, bodyIndex)
        bodyIndex_start = obs_headElem_i(obsData, OBS_RLN, headerIndex)

        burfCode = obs_bodyElem_i(obsData, OBS_VNM, bodyIndex)
        lev  = obs_bodyElem_r(obsData, OBS_PPP, bodyIndex)

        found_u = .false.
        found_v = .false.

        if (burfCode == BUFR_NEUU) then
      
          bodyIndex_u = bodyIndex
          found_u = .true.

          do i_oth = bodyIndex_start, bodyIndex
            burfCode = obs_bodyElem_i(obsData, OBS_VNM, i_oth)
            zslev = obs_bodyElem_r(obsData, OBS_PPP, i_oth)
            if (burfCode == BUFR_NEVV .and. zslev == lev) then
              bodyIndex_v = i_oth
              found_v = .true.
            end if
          end do

        else

          bodyIndex_v = bodyIndex
          found_v = .true.
 
          do i_oth = bodyIndex_start, bodyIndex
            burfCode = obs_bodyElem_i(obsData, OBS_VNM, i_oth)
            zslev = obs_bodyElem_r(obsData, OBS_PPP, i_oth)
            if (burfCode == BUFR_NEUU .and. zslev == lev) then
              bodyIndex_u = i_oth
              found_u = .true.
            end if
          end do

        end if

        if (found_u .and. found_v) then

          uu_d = obs_bodyElem_r(obsData, OBS_OMP, bodyIndex_u)
          vv_d = obs_bodyElem_r(obsData, OBS_OMP, bodyIndex_v)
          uu_r = obs_bodyElem_r(obsData, OBS_OER, bodyIndex_u)
          vv_r = obs_bodyElem_r(obsData, OBS_OER, bodyIndex_v)
          uu_f = obs_bodyElem_r(obsData, OBS_HPHT, bodyIndex_u)
          vv_f = obs_bodyElem_r(obsData, OBS_HPHT, bodyIndex_v)

          duv2 = uu_d**2 + vv_d**2

          obsFlag = 0
          duv2_lim = (uu_r**2 + uu_f**2 + vv_r**2 + vv_f**2)*1
          if (duv2 > duv2_lim) then
            obsFlag = 1
          end if
          duv2_lim = (uu_r**2 + uu_f**2 + vv_r**2 + vv_f**2)*2
          if(duv2 > duv2_lim) then
            obsFlag = 2
          end if
          duv2_lim = (uu_r**2 + uu_f**2 + vv_r**2 + vv_f**2)*4
          if(duv2 > duv2_lim) then
            obsFlag = 3
          end if

          numberObs = numberObs + 2
          if (obsFlag >= 2) numberObsRejected = numberObsRejected + 2

          if (obsFlag == 1) then
            call obs_bodySet_i(obsData,OBS_FLG,bodyIndex_u, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex_u), 13))
            call obs_bodySet_i(obsData,OBS_FLG,bodyIndex_v, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex_v), 13))

          else if (obsFlag == 2) then
            call obs_bodySet_i(obsData, OBS_FLG, bodyIndex_u, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex_u), 14))
            call obs_bodySet_i(obsData, OBS_FLG, bodyIndex_u, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex_u), 16))
            call obs_bodySet_i(obsData, OBS_FLG, bodyIndex_u, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex_u), 09))
            call obs_bodySet_i(obsData, OBS_FLG, bodyIndex_v, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex_v), 14))
            call obs_bodySet_i(obsData, OBS_FLG, bodyIndex_v, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex_v), 16))
            call obs_bodySet_i(obsData, OBS_FLG, bodyIndex_v, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex_v), 09))
            call obs_headSet_i(obsData, OBS_ST1, headerIndex, ibset(obs_headElem_i(obsData, OBS_ST1, headerIndex), 06))

          else if (obsFlag == 3) then
            call obs_bodySet_i(obsData, OBS_FLG, bodyIndex_u, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex_u), 15))
            call obs_bodySet_i(obsData, OBS_FLG, bodyIndex_u, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex_u), 16))
            call obs_bodySet_i(obsData, OBS_FLG, bodyIndex_u, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex_u), 09))
            call obs_bodySet_i(obsData, OBS_FLG, bodyIndex_v, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex_v), 15))
            call obs_bodySet_i(obsData, OBS_FLG, bodyIndex_v, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex_v), 16))
            call obs_bodySet_i(obsData, OBS_FLG, bodyIndex_v, ibset(obs_bodyElem_i(obsData, OBS_FLG, bodyIndex_v), 09))
            call obs_headSet_i(obsData, OBS_ST1, headerIndex, ibset(obs_headElem_i(obsData, OBS_ST1, headerIndex), 06))
          end if

        end if !if (found_u .and. found_v) 

      end do bodyuv

    end if !if (.not. new_bgck_sw .or. obsFamily .ne. 'SW')

    if (numberObs > 0) then
      write(*,*)' '
      write(*,*)   myName//': FINISHED BGCHECK OF ', obsFamily, ' DATA'
      write(*,123) myName//':   ', numberObsRejected, ' OBSERVATIONS REJECTED OUT OF ', numberObs
      write(*,*)' '
    end if

    if (numberObs > 0 .and. obsFamily == 'GP') then
      write(*,*)' '
      write(*,*) '  BGCDATA:    REPORT FOR GP FAMILY OF OBSERVATIONS'
      write(*,123) 'BGCDATA:   ',INZREJ, ' ZTD  OBSERVATIONS REJECTED OUT OF ', INZOBS
      write(*,123) 'BGCDATA:   ',INPREJ, ' PSFC OBSERVATIONS REJECTED OUT OF ', INPOBS
      write(*,123) 'BGCDATA:   ',INTREJ, ' TSFC OBSERVATIONS REJECTED OUT OF ', INTOBS
      write(*,123) 'BGCDATA:   ',INDREJ, ' DPDS OBSERVATIONS REJECTED OUT OF ', INDOBS
      write(*,*)' '
    end if

123   FORMAT(2X,A,I0,A,I0)

    write(*,*)' '
    write(*,*)' ---------------------------'
    write(*,*) myName//'           DONE     '
    write(*,*)' ---------------------------'
    write(*,*)' '

    if (obsCount > 0) then
      write(*,*) myName//': obsCount: ', obsCount,'; mean(FGE): ', averageFGE / obsCount
      write(*,*) myName//': obsCount: ', obsCount,'; mean(OER): ', averageOER / obsCount
    end if

  end subroutine bgck_data

  !--------------------------------------------------------------------------
  ! bgck_gpsro
  !--------------------------------------------------------------------------
  subroutine bgck_gpsro(columnTrlOnTrlLev,obsData)
      !
      !:Purpose: Set background-check flag on GPSRO data if ABS(O-P)/P is too
      !          large
      !
      implicit none

      ! Arguments:
      type(struct_columnData), intent(in)    :: columnTrlOnTrlLev
      type(struct_obs),        intent(inout) :: obsData

      ! Locals:
      type(struct_vco), pointer :: vco_trl
      real(8) :: HNH1, ZOBS, ZMHX, ZOMF, ZREF, oer, Rad
      integer :: headerIndex
      integer :: IDATYP, iProfile, varNum
      integer :: IDATA   , IDATEND, bodyIndex
      integer :: NGPSLEV
      integer :: iversion

      vco_trl => col_getVco(columnTrlOnTrlLev)
      iversion = vco_trl%vCode
      
      write(*,*)'ENTER BGCSGPSRO'
      
      ! 1.  Initializations
      !     ---------------

      NGPSLEV=col_getNumLev(columnTrlOnTrlLev,'TH')

      ! Loop over all files

      ! loop over all header indices of the 'RO' family
      call obs_set_current_header_list(obsData, 'RO')
      HEADER: do
         headerIndex = obs_getHeaderIndex(obsData)
         if (headerIndex < 0) exit HEADER
     
         ! Process only refractivity data (codtyp 169)
         IDATYP = obs_headElem_i(obsData,OBS_ITY,headerIndex)
         IF (IDATYP == 169) THEN
            iProfile = gps_iprofile_from_index(headerIndex)
            varNum = gps_vRO_IndexPrf(iProfile, 2)

            ! Basic geometric variables of the profile
            Rad  = obs_headElem_r(obsData,OBS_TRAD,headerIndex)

            ! Loops over data in the observation
            IDATA   = obs_headElem_i(obsData,OBS_RLN,headerIndex)
            IDATEND = obs_headElem_i(obsData,OBS_NLV,headerIndex) + IDATA - 1

            ! Scan for requested assimilations, and count them
            do bodyIndex= IDATA, IDATEND
               if (obs_bodyElem_i(obsData, OBS_ASS, bodyIndex) == obs_assimilated) then
                  HNH1 = obs_bodyElem_r(obsData, OBS_PPP, bodyIndex)
                  if (varNum == bufr_nebd) HNH1 = HNH1-Rad

                  ! Increment OMF = Y - H(x)
                  ZOMF = obs_bodyElem_r(obsData, OBS_OMP, bodyIndex)

                  ! Observation value    Y
                  ZOBS = obs_bodyElem_r(obsData, OBS_VAR, bodyIndex)
                  oer = obs_bodyElem_r(obsData, OBS_OER, bodyIndex)
                  ZMHX = ZOBS-ZOMF

                  ! Reference order of magnitude value:
                  if (varNum == bufr_nebd) then
                     ZREF = 0.025d0*exp(-HNH1/6500.d0)
                  else
                     if (.not. gps_roBNorm) then
                       ZREF = 300.d0*exp(-HNH1/6500.d0)
                     else
                       ZREF = ZMHX
                     end if
                  end if
                           
                  ! OMF Tested criteria:

                  ! Reject bending whose OMB is too large (>0.1 rad)
                  if (varNum == bufr_nebd) then
                    if (DABS(ZOMF) > 0.1d0) then
                      call obs_bodySet_r(obsData,OBS_OMP,bodyIndex,0.d0)
                      call obs_bodySet_i(obsData,OBS_FLG,bodyIndex,ibset(obs_bodyElem_i(obsData,OBS_FLG,bodyIndex),16))
                      call obs_bodySet_i(obsData,OBS_FLG,bodyIndex,ibset(obs_bodyElem_i(obsData,OBS_FLG,bodyIndex),9))
                      write(*,'(A40,F10.0,3F12.4)') '0 REJECT BGCSGPSRO H  O  P (O-P/ZREF) =',HNH1,ZOBS,ZMHX,(ZOMF)/ZREF
                    end if
                  end if

                  ! Reject data outside a given absolute band, or a given relative band (n sigma)
                  if (.not. gps_roBNorm) then
                    if (DABS(ZOMF)/ZREF > gps_BgckBand .or. DABS(ZOMF)/oer > 3.d0) then
                      call obs_bodySet_i(obsData,OBS_FLG,bodyIndex,ibset(obs_bodyElem_i(obsData,OBS_FLG,bodyIndex),16))
                      call obs_bodySet_i(obsData,OBS_FLG,bodyIndex,ibset(obs_bodyElem_i(obsData,OBS_FLG,bodyIndex),9))
                      write(*,'(A40,F10.0,3F12.4)') '1 REJECT BGCSGPSRO H  O  P (O-P/ZREF) =',HNH1,ZOBS,ZMHX,(ZOMF)/ZREF
                    end if
                  else
                    if (DABS(ZOMF)/ZREF > gps_BgckBand .or. DABS(ZOMF)/oer > gps_roNsigma) then
                      call obs_bodySet_i(obsData,OBS_FLG,bodyIndex,ibset(obs_bodyElem_i(obsData,OBS_FLG,bodyIndex),16))
                      call obs_bodySet_i(obsData,OBS_FLG,bodyIndex,ibset(obs_bodyElem_i(obsData,OBS_FLG,bodyIndex),9))
                      write(*,'(A40,F10.0,3F12.4)') '2 REJECT BGCSGPSRO H  O  P (O-P/ZREF) =',HNH1,ZOBS,ZMHX,(ZOMF)/ZREF
                    end if                  
                  end if
                  
               end if
            end do
         end if

      end do HEADER

      write(*,*)'EXIT BGCSGPSRO'
      RETURN
      
  end subroutine bgck_gpsro

  !--------------------------------------------------------------------------
  ! setFlag
  !--------------------------------------------------------------------------
  function setFlag(obsFamily, codeType, bufrCode, normDeparture) result(flag)
    !
    !:Purpose: Set background-check flags according to values set in a table.
    !          Original values in table come from ECMWF.
    !
    implicit none

    ! Arguments:
    character(len=2), intent(in) :: obsFamily     ! obs family name (UA,AI,...)
    integer,          intent(in) :: codeType      ! burp code type
    integer,          intent(in) :: bufrCode      ! burp variable name
    real(8),          intent(in) :: normDeparture ! normalized background departure ((y-Hx)**2/(sigmaO**2 +sigmaB**2))
    ! Result:
    integer :: flag

    ! Locals:
    real(8), target :: saCrit(3) = (/ 10.00D0, 20.00D0, 30.00D0 /)
    real(8), target :: ttCrit(3) = (/  9.00D0, 16.00D0, 25.00D0 /)
    real(8), target :: alCrit(3) = (/ 10.00D0, 20.00D0, 30.00D0 /)
    real(8), target :: esCrit(3) = (/ 10.00D0, 20.00D0, 30.00D0 /)
    real(8), target :: pnCrit(3) = (/ 10.00D0, 20.00D0, 30.00D0 /)
    real(8), target :: tsCrit(3) = (/  5.00D0, 25.00D0, 30.00D0 /)
    real(8), target :: dzCrit(3) = (/  2.25D0,  5.06D0,  7.56D0 /)
    real(8), target :: gzCrit(3) = (/ 12.25D0, 25.00D0, 36.00D0 /)
    real(8), target :: zdCrit(3) = (/  9.00D0, 16.00D0, 25.00D0 /)
    real(8), target :: chCrit(3) = (/  9.00D0, 16.00D0, 25.00D0 /)
    real(8), target :: logVisCrit(3) = (/ 10.00D0, 20.00D0, 30.00D0 /)
    real(8), target :: radvelCrit(3) = (/ 8.00D0, 20.00D0, 30.00D0 /)
    real(8), target :: uvCrit(3), swCrit(3), psCrit(3)
    real(8), pointer :: threshold_ptr(:)

    !- Dropsonde-dependent thresholds
    if (codeType == 37) then
      ! dropsondes
      uvCrit(1) = uvCritDropSonde(1)
      uvCrit(2) = uvCritDropSonde(2)
      uvCrit(3) = uvCritDropSonde(3)
      swCrit(1) = swCritDropSonde(1)
      swCrit(2) = swCritDropSonde(2)
      swCrit(3) = swCritDropSonde(3)
      psCrit(1) = psCritDropSonde(1)
      psCrit(2) = psCritDropSonde(2)
      psCrit(3) = psCritDropSonde(3)
    else
      uvCrit(1) = 10.00D0
      uvCrit(2) = 20.00D0
      uvCrit(3) = 30.00D0
      swCrit(1) = 10.00D0
      swCrit(2) = 20.00D0
      swCrit(3) = 30.00D0
      psCrit(1) =  9.00D0
      psCrit(2) = 16.00D0
      psCrit(3) = 25.00D0
    end if

    !- Associate thresholds
    if (bufrCode == bufr_negz) then     
      ! height
      threshold_ptr => gzCrit
    else if (bufrCode == bufr_nett) then
      ! temperature
      threshold_ptr => ttCrit
    else if (bufrCode == bufr_nedz) then     
      ! satems
      threshold_ptr => dzCrit
    else if (bufrCode == bufr_nefs) then
      ! wind speed
      threshold_ptr => saCrit
    else if (bufrCode == bufr_neuu .or. bufrCode == bufr_nevv) then
      ! wind components
      threshold_ptr => uvCrit
    else if (bufrCode == bufr_neal) then
      ! aladin hlos wind observations
      threshold_ptr => alCrit
    else if (bufrCode == bufr_neus .or. bufrCode == bufr_nevs) then
      ! 10m wind components
      threshold_ptr => swCrit
    else if (bufrCode == bufr_gust) then
      ! 10m wind gust
      threshold_ptr => swCrit
    else if (bufrCode == bufr_nees) then
      ! dewpoint depression
      threshold_ptr => esCrit
    else if (bufrCode == bufr_neps) then
      ! surface pressure
      threshold_ptr => psCrit
    else if (bufrCode == bufr_nepn) then
      ! mean sea level pressures
      threshold_ptr => pnCrit
    else if (bufrCode == bufr_nets) then
      ! 1.5m temperature
      threshold_ptr => tsCrit
    else if (bufrCode == bufr_ness) then
      ! 1.5m dew point depression
      threshold_ptr => esCrit
    else if (bufrCode == bufr_vis) then
      ! visibility
      threshold_ptr => logVisCrit
    else if (bufrCode == bufr_nezd) then
      ! gb-gps zenith delay
      threshold_ptr => zdCrit
    else if (bufrCode == bufr_radvel) then    
      ! doppler velocity
      threshold_ptr => radvelCrit
    else if (obsFamily == 'CH') then
      ! chemical constituents
      threshold_ptr => chCrit
    else
      write(*,*) 'setFlag: non-defined obs = ',obsFamily, codeType, bufrCode
      flag=0
      return
    end if 

    !- Set flag
    flag=0
    if      (normDeparture >= threshold_ptr(1) .and. normDeparture < threshold_ptr(2)) then
      flag=1
    else if (normDeparture >= threshold_ptr(2) .and. normDeparture < threshold_ptr(3)) then
      flag=2
    else if (normDeparture >= threshold_ptr(3)) then
      flag=3
    end if

  end function setFlag

end module backgroundCheck_mod
