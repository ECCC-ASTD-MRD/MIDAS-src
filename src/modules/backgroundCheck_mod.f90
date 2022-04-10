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

module backgroundCheck_mod
  ! MODULE backgroundCheck_mod (prefix='bgck' category='1. High-level functionality')
  !
  ! :Purpose: Performs the background check on all conventional observations
  !
  use mathPhysConstants_mod
  use bufr_mod
  use obsSpaceData_mod
  use gps_mod
  use utilities_mod
  use mpi_mod
  use mpivar_mod
  use columnData_mod
  use obsSpaceDiag_mod
  use earthConstants_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use obsSpaceErrorStdDev_mod
  use obsFamilyList_mod
  use codtyp_mod
  
  implicit none
  private

  ! public procedures
  public :: bgck_bgcheck_conv

  integer,          parameter :: numFamilyToProcess = 11
  character(len=2), parameter :: familyListToProcess(numFamilyToProcess)= (/ &
                                 'UA','AI','SF','SC','SW','PR','RO','GP','RA', &
                                 'CH','AL' /)

  contains

  !--------------------------------------------------------------------------
  ! bgck_bgcheck_conv
  !--------------------------------------------------------------------------
  subroutine bgck_bgcheck_conv( columnTrlOnAnlIncLev, columnTrlOnTrlLev, hco_anl, obsSpaceData )
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
    integer                     :: nulNam, ier, fnom, fclos
    logical                     :: new_bgck_sw, noObsToProcess
    character(len=*), parameter :: myName = 'bgck_bgcheck_conv'

    write(*,*) myName//' begin conventional data background check'

    ! Check if any observations are present for conventional background check
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
     
    call tmg_start(3, myName )

    new_bgck_sw = .false.

    namelist /NAMBGCKCONV/ new_bgck_sw
    nulNam=0
    ier = fnom( nulNam, 'flnml', 'r/o', 0 )
    read( nulNam, nml = NAMBGCKCONV, IOSTAT = ier )
    if ( ier /= 0 ) write(*,*) myName//': no valid namelist NAMBGCKCONV found, default values will be taken:'
    write(*,*) myName//': new_bgck_sw = ',new_bgck_sw
    ier = fclos(nulNam)
     
 
    ! Obtain or calc OmP-error std dev when requested and possible.
    ! Otherwise calc HBHT contribution (sigma_B in observation space)  
    ! -------------------------------------------------------------------- 
    call ose_computeStddev( columnTrlOnAnlIncLev, hco_anl, & ! IN
                            obsSpaceData )                   ! INOUT

    ! DO A BACKGROUND CHECK ON ALL THE OBSERVATIONS
    ! ----------------------------------------------

    do familyIndex = 1, ofl_numFamily
      ! For SW only, old and new background check schemes controlled by "new_bgck_sw"
      if ( obs_famExist( obsSpaceData, ofl_familyList( familyIndex )) ) then
        call bgck_data( ofl_familyList(familyIndex), obsSpaceData, new_bgck_sw )
      end if
    end do

    if (obs_famExist( obsSpaceData, 'RO' )) call bgck_gpsro( columnTrlOnTrlLev , obsSpaceData )

    ! Conduct obs-space post-processing diagnostic tasks (some diagnostic 
    ! computations controlled by NAMOSD namelist in flnml)

    call osd_ObsSpaceDiag( obsSpaceData, columnTrlOnAnlIncLev, hco_anl, analysisMode_opt = .false. )

    call tmg_stop(3)

  end subroutine bgck_bgcheck_conv

  !--------------------------------------------------------------------------
  ! bgck_data
  !--------------------------------------------------------------------------
  subroutine bgck_data( obsFamily, obsData, new_bgck_sw )
    !
    !:Purpose:  Calculate a background check for a data family and set the
    !           appropriate quality-control flags in obsSpaceData
   
    implicit none

    ! NOTE 1: YSFERRWGT IN MODGPSZTD_MOD (FROM NML FILE) IS USED HERE FOR ERROR WEIGHTING
    !         OF TIME SERIES (FGAT) GPS MET OBSERVATIONS PS, TS, DPDS. IT IS APPLIED
    !         (ALONG WITH YZDERRWGT FOR ZTD) IN S/R SETERR AS A MULT. FACTOR TO ERRORS.
    !         YZDERRWGT AND YSFERRWGT = 1 FOR NORMAL 3D-VAR (3D-THINNING).

    ! Arguments
    character(len=2), intent(in)    :: obsFamily   ! current observation family
    type(struct_obs), intent(inout) :: obsData     ! obsSpaceData
    logical         , intent(in)    :: new_bgck_sw
    
    real(8) :: averageFGE, averageOER
    integer :: obsFlag, obsVarno, headerIndex, codeType, bodyIndex, obsCount
    integer :: numberObs, numberObsRejected, INZOBS, INZREJ
    integer :: INPOBS, INTOBS, INDOBS, INPREJ, INTREJ, INDREJ
    real(8) :: ZOER,ZOMP,ZFGE,ZOMPER,ZBGCHK,ZVAR,ZLEV,ZLAT,ZLON,ZSOP
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
    averageOER=0.d0

    numberObs = 0
    numberObsRejected = 0

    ! Initialize counters for GP family observations

    if ( obsFamily == 'GP' ) then
      INZOBS=0
      INPOBS=0
      INTOBS=0
      INDOBS=0
      INZREJ=0
      INPREJ=0
      INTREJ=0
      INDREJ=0
    end if

    if ( .not. new_bgck_sw .or. obsFamily /= 'SW' ) then

      ! loop over all header indices of the current observation family
      call obs_set_current_header_list( obsData, obsFamily )
      HEADER: do
        headerIndex = obs_getHeaderIndex( obsData )
        if ( headerIndex < 0 ) exit HEADER

        stnid = obs_elem_c( obsData, 'STID', headerIndex )

        ! 1. Computation of (HX - Z)**2/(SIGMAo**2 +SIGMAp**2)
        ! ----------------------------------------------------
          
        ! loop over all body indices for the current headerIndex
        call obs_set_current_body_list( obsData, headerIndex )
        BODY: do 
          bodyIndex = obs_getBodyIndex( obsData )
          if (bodyIndex < 0 ) exit BODY

          obsVarno = obs_bodyElem_i( obsData, OBS_VNM, bodyIndex )
          LLOK = ( obs_bodyElem_i( obsData, OBS_ASS, bodyIndex ) == obs_assimilated )
          if ( LLOK ) then
            numberObs = numberObs + 1
            if ( obsFamily == 'GP' ) then
              if ( obsVarno == BUFR_NEZD ) INZOBS = INZOBS+1
              if ( obsVarno == BUFR_NEPS ) INPOBS = INPOBS+1
              if ( obsVarno == BUFR_NETS ) INTOBS = INTOBS+1
              if ( obsVarno == BUFR_NESS ) INDOBS = INDOBS+1
            end if
            ZVAR = obs_bodyElem_r( obsData, OBS_VAR, bodyIndex )
            ZLEV = obs_bodyElem_r( obsData, OBS_PPP, bodyIndex )              
            ZLAT  = obs_headElem_r( obsData, OBS_LAT, headerIndex ) * MPC_DEGREES_PER_RADIAN_R8
            ZLON  = obs_headElem_r( obsData, OBS_LON, headerIndex ) * MPC_DEGREES_PER_RADIAN_R8

            ! BACKGROUND CHECK

            ZOMP  = obs_bodyElem_r( obsData, OBS_OMP, bodyIndex )

            ! Get error std dev
              
            ! std(O-P) (valid/available if > 0.0)
            zomper = obs_bodyElem_r( obsData, OBS_OMPE, bodyIndex )
                                         
            ! std(O)
            zoer = obs_bodyElem_r( obsData, OBS_OER, bodyIndex )
            ! std(P)
            zfge  = obs_bodyElem_r( obsData, OBS_HPHT, bodyIndex )
            ! NOTE: For GB-GPS ZTD observations (GP family), ZFGE is not the FGE but
            !       rather Std(O-P) set in s/r SETERRGPSGB.
            if ( obsFamily == 'GP' .and. obsVarno == BUFR_NEZD ) then
              if ( zomper <= 0.0d0 .and. zfge > 0.0d0 ) zomper = zfge
            end if
              
            ! Protect against error std dev values that are too small
            if ( obsFamily == 'GP' ) then
              if ( obsVarno == BUFR_NEZD ) then
                ZOER = ZOER / YZDERRWGT
              else
                ZOER = ZOER / YSFERRWGT
              end if
              if ( ZOER < 1.0d-3 .and. obsVarno /= BUFR_NEZD ) then
                write(*,*)' Problem for GP STNID ZOER= ' , stnid, ZOER
                call utl_abort( myName//': PROBLEM WITH OER.')
              end if
              IF ( ZFGE < 1.0d-3 ) then
                write(*,*)' Problem for GP STNID FGE= ', stnid, ZFGE
                ZFGE=1.0d-3
              end if
              if ( zomper > 0.0d0 .and. zomper < 1.0d-3 ) zomper = 1.0d-3
            else if ( obsFamily == 'CH' ) then
              if ( ZFGE**2 + ZOER**2 < 1.0d-60) then
                write(*,*) ' Problem for STNID FGE ZOER=',  stnid, ZFGE, ZOER
                ZFGE=1.0d30
                ZOER=1.0d-30
              end if
              if ( zomper > 0.0 .and. zomper < 1.0d-30 ) zomper = 1.0d-30
            else
              if ( zfge < 0.0d0 ) zfge = 0.0d0
              if ( ZFGE**2 + ZOER**2 < 1.0d-5)then
                write(*,*) ' Problem for STNID FGE ZOER=', stnid, ZFGE, ZOER
                ZFGE=1.0d-5
                ZOER=1.0d-5
              end if
              if ( zomper > 0.0 .and. zomper < 1.0D-5 ) zomper = 1.0d-5
            end if
              
            ! Calculate zbgchk
              
            if ( zomper > 0.0d0 ) then
              zbgchk = (zomp**2)/(zomper**2) 
              zsop = zomper
            else
              ZBGCHK=(ZOMP)**2/(ZFGE**2 + ZOER**2)
              ZSOP = SQRT(ZFGE**2 + ZOER**2)
            end if
              
            ! UPDATE QUALITY CONTROL FLAGS, based on zbgchk
            ! ( ELEMENT FLAGS + GLOBAL HEADER FLAGS)
            obsVarno = obs_bodyElem_i( obsData, OBS_VNM, bodyIndex )
            if( obsVarno == bufr_nees .and. obs_bodyElem_i( obsData, OBS_XTR, bodyIndex ) == 0 ) then
              averageFGE = averageFGE + zfge * zfge
              averageOER = averageOER + zoer * zoer
              obsCount   = obsCount + 1
            end if
            codeType = obs_headElem_i( obsData, OBS_ITY, headerIndex )
            obsFlag = ISETFLAG( obsFamily, codeType, obsVarno, ZBGCHK )

            ! CONVERT ZTD VALUES FROM M TO MM FOR PRINTOUT

            LLZD = .FALSE.
            if ( obsFamily == 'GP' .and. obsVarno == BUFR_NEZD ) then
              ZVAR = ZVAR * 1000.0D0
              ZOER = ZOER * 1000.0D0
              ZFGE = ZFGE * 1000.0D0
              ZOMP = ZOMP * 1000.0D0
              ZSOP = ZSOP * 1000.0D0
              LLZD = .TRUE.
            end if

            stnid = obs_elem_c( obsData, 'STID', headerIndex )
            if ( obsFlag >= 2 .or. ( LLZD .and. LTESTOP )) then
              if ( obsFamily /= 'CH' ) then 
                write(*,122)  &
                   stnid, zlat, zlon, codeType, obsVarno, ZLEV, ZVAR, ZOER,  &
                   ZFGE, ZOMP, ZSOP, ZBGCHK, obsFlag
              else 
                write(*,124)  &
                   stnid, zlat, zlon, codeType, obsVarno, ZLEV, ZVAR, ZOER,  &
                   ZFGE, ZOMP, ZSOP, ZBGCHK, obsFlag
              end if
              if ( obsFlag >= 2 ) numberObsRejected = numberObsRejected + 1
              if ( obsFlag >= 2 .and. obsFamily == 'GP' ) then
              if ( obsVarno == BUFR_NEZD ) INZREJ = INZREJ+1
              if ( obsVarno == BUFR_NEPS ) INPREJ = INPREJ+1
              if ( obsVarno == BUFR_NETS ) INTREJ = INTREJ+1
              if ( obsVarno == BUFR_NESS ) INDREJ = INDREJ+1
            end if
          end if

          if ( obsFlag == 1 ) then
            call obs_bodySet_i( obsData, OBS_FLG, bodyIndex, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex ), 13 ))
          else if ( obsFlag == 2 ) then
            call obs_bodySet_i( obsData, OBS_FLG, bodyIndex, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex ), 14 ))
            call obs_bodySet_i( obsData, OBS_FLG, bodyIndex, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex ), 16 ))
            call obs_bodySet_i( obsData, OBS_FLG, bodyIndex, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex ), 09 ))
            call obs_headSet_i( obsData, OBS_ST1, headerIndex, ibset( obs_headElem_i( obsData, OBS_ST1, headerIndex ), 06 ))

          else if ( obsFlag == 3 ) then
            call obs_bodySet_i( obsData, OBS_FLG, bodyIndex, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex), 15 ))
            call obs_bodySet_i( obsData, OBS_FLG, bodyIndex, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex), 16 ))
            call obs_bodySet_i( obsData, OBS_FLG, bodyIndex, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex), 09 ))
            call obs_headSet_i( obsData, OBS_ST1, headerIndex, ibset( obs_headElem_i( obsData, OBS_ST1, headerIndex), 06 ))
          end if
        end if
      end do BODY

124       FORMAT(2x,a9,1x,f6.2,1x,f7.2,1x,I3,1x,I5,7(2x,G11.2),I3 )
122       FORMAT(2x,a9,1x,f6.2,1x,f7.2,1x,I3,1x,I5,7(2x,F11.2),I3 )

    end do HEADER

    else !IF (.not. new_bgck_sw .or. obsFamily .ne. 'SW') then

      call obs_set_current_body_list( obsData, 'SW' )
      bodyuv: do
        bodyIndex = obs_getBodyIndex(obsData)
        if ( bodyIndex < 0 ) exit bodyuv

        ! Only process pressure level observations flagged to be assimilated
        i_ass = obs_bodyElem_i ( obsData, OBS_ASS, bodyIndex )
        i_vco = obs_bodyElem_i ( obsData, OBS_VCO, bodyIndex )

        if( i_ass /= obs_assimilated .or. i_vco /= 2 ) cycle bodyuv
 
        headerIndex     = obs_bodyElem_i( obsData, OBS_HIND, bodyIndex )
        bodyIndex_start = obs_headElem_i( obsData, OBS_RLN, headerIndex )

        obsVarno = obs_bodyElem_i( obsData, OBS_VNM, bodyIndex )
        zlev  = obs_bodyElem_r( obsData, OBS_PPP, bodyIndex )

        found_u = .false.
        found_v = .false.

        if ( obsVarno == BUFR_NEUU ) then
      
          bodyIndex_u = bodyIndex
          found_u = .true.

          do i_oth = bodyIndex_start, bodyIndex
            obsVarno = obs_bodyElem_i( obsData, OBS_VNM, i_oth )
            zslev = obs_bodyElem_r( obsData, OBS_PPP, i_oth )
            if ( obsVarno == BUFR_NEVV .and. zslev == zlev ) then
              bodyIndex_v = i_oth
              found_v = .true.
            end if
          end do

        else

          bodyIndex_v = bodyIndex
          found_v = .true.
 
          do i_oth = bodyIndex_start, bodyIndex
            obsVarno = obs_bodyElem_i( obsData, OBS_VNM, i_oth )
            zslev = obs_bodyElem_r( obsData, OBS_PPP, i_oth )
            if ( obsVarno == BUFR_NEUU .and. zslev == zlev ) then
              bodyIndex_u = i_oth
              found_u = .true.
            end if
          end do

        end if

        if ( found_u .and. found_v ) then

          uu_d = obs_bodyElem_r( obsData, OBS_OMP, bodyIndex_u )
          vv_d = obs_bodyElem_r( obsData, OBS_OMP, bodyIndex_v )
          uu_r = obs_bodyElem_r( obsData, OBS_OER, bodyIndex_u )
          vv_r = obs_bodyElem_r( obsData, OBS_OER, bodyIndex_v )
          uu_f = obs_bodyElem_r( obsData, OBS_HPHT, bodyIndex_u )
          vv_f = obs_bodyElem_r( obsData, OBS_HPHT, bodyIndex_v )

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
          if ( obsFlag >= 2 ) numberObsRejected = numberObsRejected + 2

          if ( obsFlag == 1 ) then
            call obs_bodySet_i( obsData,OBS_FLG,bodyIndex_u, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex_u ), 13 ))
            call obs_bodySet_i( obsData,OBS_FLG,bodyIndex_v, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex_v ), 13 ))

          else if ( obsFlag == 2 ) then
            call obs_bodySet_i( obsData, OBS_FLG, bodyIndex_u, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex_u ), 14 ))
            call obs_bodySet_i( obsData, OBS_FLG, bodyIndex_u, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex_u ), 16 ))
            call obs_bodySet_i( obsData, OBS_FLG, bodyIndex_u, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex_u ), 09 ))
            call obs_bodySet_i( obsData, OBS_FLG, bodyIndex_v, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex_v ), 14 ))
            call obs_bodySet_i( obsData, OBS_FLG, bodyIndex_v, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex_v ), 16 ))
            call obs_bodySet_i( obsData, OBS_FLG, bodyIndex_v, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex_v ), 09 ))
            call obs_headSet_i( obsData, OBS_ST1, headerIndex, ibset( obs_headElem_i( obsData, OBS_ST1, headerIndex ), 06 ))

          else if ( obsFlag == 3 ) then
            call obs_bodySet_i( obsData, OBS_FLG, bodyIndex_u, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex_u ), 15 ))
            call obs_bodySet_i( obsData, OBS_FLG, bodyIndex_u, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex_u ), 16 ))
            call obs_bodySet_i( obsData, OBS_FLG, bodyIndex_u, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex_u ), 09 ))
            call obs_bodySet_i( obsData, OBS_FLG, bodyIndex_v, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex_v ), 15 ))
            call obs_bodySet_i( obsData, OBS_FLG, bodyIndex_v, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex_v ), 16 ))
            call obs_bodySet_i( obsData, OBS_FLG, bodyIndex_v, ibset( obs_bodyElem_i( obsData, OBS_FLG, bodyIndex_v ), 09 ))
            call obs_headSet_i( obsData, OBS_ST1, headerIndex, ibset( obs_headElem_i( obsData, OBS_ST1, headerIndex ), 06 ))
          end if

        end if !if(found_u .and. found_v) 

      end do bodyuv

    end if !IF (.not. new_bgck_sw .or. obsFamily .ne. 'SW')

    if ( numberObs > 0 ) then
      write(*,*)' '
      write(*,*)   myName//': FINISHED BGCHECK OF ', obsFamily, ' DATA'
      write(*,123) myName//':   ', numberObsRejected, ' OBSERVATIONS REJECTED OUT OF ', numberObs
      write(*,*)' '
    end if

    if ( numberObs > 0 .and. obsFamily == 'GP' ) then
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

    if ( obsCount > 0) then
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

      type(struct_columnData) :: columnTrlOnTrlLev
      type(struct_obs) :: obsData
      type(struct_vco), pointer :: vco_trl
      real(8) :: HNH1, ZOBS, ZMHX, ZOMF, ZREF, ZOER, Rad

      integer :: headerIndex
      integer :: IDATYP, iProfile, varNum
      integer :: IDATA   , IDATEND, bodyIndex
      integer :: NGPSLEV
      integer :: stat, iversion

      vco_trl => col_getVco(columnTrlOnTrlLev)
      stat = vgd_get(vco_trl%vgrid,key='ig_1 - vertical coord code',value=iversion)
      
      write(*,*)'ENTER BGCSGPSRO'
      
      ! 1.  Initializations
      !     ---------------

      NGPSLEV=col_getNumLev(columnTrlOnTrlLev,'TH')

      ! Loop over all files

      ! loop over all header indices of the 'RO' family
      call obs_set_current_header_list( obsData, 'RO' )
      HEADER: do
         headerIndex = obs_getHeaderIndex(obsData)
         if (headerIndex < 0) exit HEADER
     
         ! Process only refractivity data (codtyp 169)

         IDATYP = obs_headElem_i(obsData,OBS_ITY,headerIndex)
         IF ( IDATYP == 169 ) THEN
            iProfile = gps_iprofile_from_index(headerIndex)
            varNum = gps_vRO_IndexPrf(iProfile, 2)

            ! Basic geometric variables of the profile
            
            Rad  = obs_headElem_r(obsData,OBS_TRAD,headerIndex)

            ! Loops over data in the observation
            
            IDATA   = obs_headElem_i(obsData,OBS_RLN,headerIndex)
            IDATEND = obs_headElem_i(obsData,OBS_NLV,headerIndex) + IDATA - 1

            ! Scan for requested assimilations, and count them

            do bodyIndex= IDATA, IDATEND
               if ( obs_bodyElem_i( obsData, OBS_ASS, bodyIndex ) == obs_assimilated ) then
                  HNH1 = obs_bodyElem_r( obsData, OBS_PPP, bodyIndex )
                  if ( varNum == bufr_nebd ) HNH1 = HNH1-Rad

                  ! Increment OMF = Y - H(x)

                  ZOMF = obs_bodyElem_r( obsData, OBS_OMP, bodyIndex )

                  ! Observation value    Y

                  ZOBS = obs_bodyElem_r( obsData, OBS_VAR, bodyIndex )
                  ZOER = obs_bodyElem_r( obsData, OBS_OER, bodyIndex )
                  ZMHX = ZOBS-ZOMF

                  ! Reference order of magnitude value:

                  if ( varNum == bufr_nebd ) then
                     ZREF = 0.025d0*exp(-HNH1/6500.d0)
                  else
                     if ( .not. gpsroBNorm ) then
                       ZREF = 300.d0*exp(-HNH1/6500.d0)
                     else
                       ZREF = ZMHX
                     end if
                  end if
                           
                  ! OMF Tested criteria:

                  if ( .not. gpsroBNorm ) then
                    if (DABS(ZOMF)/ZREF > BGCKBAND .or. DABS(ZOMF)/ZOER > 3.d0) then
                      call obs_bodySet_i(obsData,OBS_FLG,bodyIndex,ibset(obs_bodyElem_i(obsData,OBS_FLG,bodyIndex),16))
                      call obs_bodySet_i(obsData,OBS_FLG,bodyIndex,ibset(obs_bodyElem_i(obsData,OBS_FLG,bodyIndex),9))
                    end if
                  else
                    if ( DABS(ZOMF)/ZREF > BGCKBAND ) then
                      call obs_bodySet_i(obsData,OBS_FLG,bodyIndex,ibset(obs_bodyElem_i(obsData,OBS_FLG,bodyIndex),16))
                      call obs_bodySet_i(obsData,OBS_FLG,bodyIndex,ibset(obs_bodyElem_i(obsData,OBS_FLG,bodyIndex),9))
                      write(*,'(A40,F10.0,3F12.4)') ' REJECT BGCSGPSRO H  O  P (O-P/ZREF) =',HNH1,ZOBS,ZMHX,(ZOMF)/ZREF
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
  ! isetflag
  !--------------------------------------------------------------------------
  function isetflag( obsFamily, kodtyp, kvnam, zbgchk )
    !
    !:Purpose: Set BACKGROUND-CHECK FLAGS According to values set in a table.
    !          Original values in table come from ecmwf.
    !

    implicit none
    integer isetflag

    ! Arguments:
    character(len=2) :: obsFamily ! FAMILY  NAME ( 'UA' , 'AI'   ...etc.. )
    integer :: kodtyp ! BURP CODE TYPE
    integer :: kvnam  ! VARIABLE NAME ( BURP )
    real(8) :: zbgchk ! NORMALIZED BACKGROUND DEPARTURE

    ! Locals:      
    real(8), parameter :: zsacrit(3) = (/ 10.00D0, 20.00D0, 30.00D0 /)
    real(8), parameter :: zttcrit(3) = (/  9.00D0, 16.00D0, 25.00D0 /)
    real(8), parameter :: zalcrit(3) = (/ 10.00D0, 20.00D0, 30.00D0 /)
    real(8), parameter :: zescrit(3) = (/ 10.00D0, 20.00D0, 30.00D0 /)
    real(8), parameter :: zpscrit(3) = (/  9.00D0, 16.00D0, 25.00D0 /)
    real(8), parameter :: zpncrit(3) = (/ 10.00D0, 20.00D0, 30.00D0 /)
    real(8), parameter :: zswcrit(3) = (/ 10.00D0, 20.00D0, 30.00D0 /)
    real(8), parameter :: ztscrit(3) = (/  5.00D0, 25.00D0, 30.00D0 /)
    real(8), parameter :: zdzcrit(3) = (/  2.25D0,  5.06D0,  7.56D0 /)
    real(8), parameter :: zgzcrit(3) = (/ 12.25D0, 25.00D0, 36.00D0 /)
    real(8), parameter :: zzdcrit(3) = (/  9.00D0, 16.00D0, 25.00D0 /)
    real(8), parameter :: zchcrit(3) = (/  9.00D0, 16.00D0, 25.00D0 /)
    real(8), parameter :: zLogViscrit(3) = (/ 10.00D0, 20.00D0, 30.00D0 /)

    !real(8) :: zuvcrit(3) = (/ 10.00D0, 20.00D0, 30.00D0 /)
    real(8) :: zuvcrit(3)
      
    zuvcrit(1) = 10.00D0
    zuvcrit(3) = 30.00D0
    if ( kodtyp == 37 ) then
      zuvcrit(2)=25.00D0
    else
      zuvcrit(2)=20.00D0
    end if
      
    isetflag=0

    if ( kvnam == BUFR_NEGZ ) then
      
      ! SET FLAG FOR HEIGHTS
         
      if ( zbgchk >= zgzcrit(1) .and. zbgchk < zgzcrit(2) ) then
        isetflag=1
      else if ( zbgchk >= zgzcrit(2) .and. zbgchk < zgzcrit(3) ) then
        isetflag=2
      else if ( zbgchk >= zgzcrit(3) )then
        isetflag =3
      end if

    else if ( kvnam == BUFR_NETT ) then
      
      ! SET FLAG FOR TEMPERATURE

      if ( zbgchk >= zttcrit(1) .and. zbgchk < zttcrit(2) ) then
        isetflag=1
      else if ( zbgchk >= zttcrit(2) .and. zbgchk < zttcrit(3) ) then
        isetflag=2
      else if ( zbgchk >= zttcrit(3) )then
        isetflag =3
      end if

    else if ( kvnam == BUFR_NEDZ ) then
     
      ! SET FLAG FOR SATEMS
         
      if ( zbgchk >= zdzcrit(1) .and. zbgchk < zdzcrit(2) ) then
        isetflag=1
      else if ( zbgchk >= zdzcrit(2) .and. zbgchk < zdzcrit(3) ) then
        isetflag=2
      else if ( zbgchk >= zdzcrit(3) )then
        isetflag =3
      end if

    else if ( kvnam == BUFR_NEFS ) then

      ! SET FLAG FOR WIND SPEED

      if ( zbgchk >= zsacrit(1) .and. zbgchk < zsacrit(2) ) then
        isetflag=1
      else if ( zbgchk >= zsacrit(2) .and. zbgchk < zsacrit(3) ) then
        isetflag=2
      else if ( zbgchk >= zsacrit(3) )then
        isetflag =3
      end if

    else if ( kvnam == BUFR_NEUU .or. kvnam == BUFR_NEVV ) then

      ! SET FLAG FOR WIND COMPONENTS

      if ( zbgchk >= zuvcrit(1) .and. zbgchk < zuvcrit(2) ) then
        isetflag=1
      else if ( zbgchk >= zuvcrit(2) .and. zbgchk < zuvcrit(3) ) then
        isetflag=2
      else if ( zbgchk >= zuvcrit(3) )then
        isetflag =3
      end if
         
    else if ( kvnam == BUFR_NEAL ) then

      ! SET FLAG FOR ALADIN HLOS WIND OBSERVATIONS
 
      if ( zbgchk >= zalcrit(1) .and. zbgchk < zalcrit(2) ) then
        isetflag=1
      else if ( zbgchk >= zalcrit(2) .and. zbgchk < zalcrit(3) ) then
        isetflag=2
      else if ( zbgchk >= zalcrit(3) )then
        isetflag=3
      end if

    else if ( kvnam == BUFR_NEUS .or. kvnam == BUFR_NEVS ) then

      ! SET FLAG FOR SURFACE WIND COMPONENTS

      if ( zbgchk >= zswcrit(1) .and. zbgchk < zswcrit(2) ) then
        isetflag=1
      else if ( zbgchk >= zswcrit(2) .and. zbgchk < zswcrit(3) ) then
        isetflag=2
      else if ( zbgchk >= zswcrit(3) )then
        isetflag =3
      end if

    else if ( kvnam == bufr_gust ) then
      
      ! SET FLAG FOR SURFACE WIND GUST

      if ( zbgchk >= zswcrit(1) .and. zbgchk < zswcrit(2) ) then
        isetflag=1
      else if ( zbgchk >= zswcrit(2) .and. zbgchk < zswcrit(3) ) then
        isetflag=2
      else if ( zbgchk >= zswcrit(3) )then
        isetflag =3
      end if

    else if ( kvnam == bufr_nees ) then

      ! SET FLAG FOR DEW POINT DEPRESSION

      if ( zbgchk >= zescrit(1) .and. zbgchk < zescrit(2) ) then
        isetflag=1
      else if ( zbgchk >= zescrit(2) .and. zbgchk < zescrit(3) ) then
        isetflag=2
      else if ( zbgchk >= zescrit(3) )then
        isetflag =3
      end if

    else if ( kvnam == bufr_neps ) then

      ! SET FLAG FOR SURFACE PRESSURE
         
      if ( zbgchk >= zpscrit(1) .and. zbgchk < zpscrit(2) ) then
        isetflag=1
      else if ( zbgchk >= zpscrit(2) .and. zbgchk < zpscrit(3) ) then
        isetflag=2
      else if ( zbgchk >= zpscrit(3) )then
        isetflag =3
      end if

    else if ( kvnam == bufr_nepn ) then

      ! SET FLAG FOR MEAN SEA LEVEL PRESSURE

      if ( zbgchk >= zpncrit(1) .and. zbgchk < zpncrit(2) ) then
        isetflag=1
      else if ( zbgchk >= zpncrit(2) .and. zbgchk < zpncrit(3) ) then
        isetflag=2
      else if ( zbgchk >= zpncrit(3) ) then
        isetflag =3
      end if

    else if ( kvnam == bufr_nets ) then

      ! SET FLAG FOR SURFACE TEMPERATURE

      if ( zbgchk >= ztscrit(1) .and. zbgchk < ztscrit(2) ) then
        isetflag=1
      else if ( zbgchk >= ztscrit(2) .and. zbgchk < ztscrit(3) ) then
        isetflag=2
      else if ( zbgchk >= ztscrit(3) ) then
        isetflag =3
      end if

   else if ( kvnam == BUFR_NESS ) then
      
      ! SET FLAG FOR SURFACE DEW POINT DEPRESSION

      if ( zbgchk >= zescrit(1) .and. zbgchk < zescrit(2) ) then
        isetflag=1
      else if ( zbgchk >= zescrit(2) .and. zbgchk < zescrit(3) ) then
        isetflag=2
      else if ( zbgchk >= zescrit(3) ) then
        isetflag =3
      end if

   else if ( kvnam == bufr_vis ) then

      ! SET FLAG FOR VISIBILITY

      if ( zbgchk >= zLogViscrit(1) .and. zbgchk < zLogViscrit(2) ) then
        isetflag=1
      else if ( zbgchk >= zLogViscrit(2) .and. zbgchk < zLogViscrit(3) ) then
        isetflag=2
      else if ( zbgchk >= zLogViscrit(3) ) then
        isetflag =3
      end if

   else if ( kvnam == bufr_nezd ) then

      ! SET FLAG FOR GB-GPS ZENITH DELAY

      if ( zbgchk >= zzdcrit(1) .and. zbgchk < zzdcrit(2) ) then
        isetflag=1
      else if ( zbgchk >= zzdcrit(2) .and. zbgchk < zzdcrit(3) ) then
        isetflag=2
      else if ( zbgchk >= zzdcrit(3) ) then
        isetflag =3
      end if

    else if ( obsFamily == 'CH' ) then
    
      ! SET FLAG FOR CHEMICAL CONSTITUENTS
      if ( zbgchk >= zchcrit(1) .and. zbgchk < zchcrit(2) ) then
        isetflag=1
      else if ( zbgchk >= zchcrit(2) .and. zbgchk < zchcrit(3) ) then
        isetflag=2
      else if ( zbgchk >= zchcrit(3) ) then
        isetflag =3
      end if
         
    end if

    return
  
  end function isetflag
      
end module backgroundCheck_mod
