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
!! MODULE obsOperators (prefix="oop")
!!
!! *Purpose*: All observation operators, including nonlinear, tangent-linear
!!            and adjoint versions.
!!
!--------------------------------------------------------------------------
module obsOperators_mod
  use earthConstants_mod
  use mathPhysConstants_mod
  use obsSpaceData_mod
  use columnData_mod 
  use bufr_mod
  use physicsFunctions_mod
  use gps_mod
  use mpivar_mod
  use timeCoord_mod
  use obsFilter_mod
  use tovs_nl_mod
  use utilities_mod
  use tovs_lin_mod
  use chem_obsoperators_mod
  
  implicit none
  save
  private

  ! public procedures
  public :: oop_setup
  public :: oop_ppp_nl, oop_sfc_nl, oop_zzz_nl, oop_gpsro_nl, oop_gpsgb_nl, oop_tovs_nl, oop_chm_nl
  public :: oop_Htl, oop_Had, oop_vobslyrs

  character(len=48) :: obsoperMode

  integer, external :: get_max_rss

contains

  subroutine oop_setup(obsoperMode_in)
    character(len=*), intent(in) :: obsoperMode_in

    obsoperMode = obsoperMode_in

  end subroutine oop_setup

  subroutine oop_vobslyrs(columnghr,obsSpaceData)

    !*    Purpose:
    !      Find which model levels to use for the vertical interpolation
    !      of model fields to CMA data.
    !
    IMPLICIT NONE
    type(struct_columnData) :: columnghr
    type(struct_obs) :: obsSpaceData

    INTEGER :: JK,JDATA,NLEV
    REAL(8) :: ZLEV,ZPT,ZPB
    INTEGER :: IOBS,IK,ITYP
    LOGICAL :: LLOK
    CHARACTER(len=2) :: varLevel
    integer :: index_header, index_body
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
    Write(*,*) "Entering subroutine OOP_VOBSLYRS"

!$OMP PARALLEL DO PRIVATE(jdata,llok,zlev,iobs,ityp,varLevel,zpt,zpb)
    DO JDATA= 1,obs_numbody(obsSpaceData)
       LLOK = ( (obs_bodyElem_i(obsSpaceData,OBS_ASS,JDATA) == 1     .OR. &
            obs_bodyElem_i(obsSpaceData,OBS_ASS,JDATA) == -1) .AND. &
            obs_bodyElem_i(obsSpaceData,OBS_VCO,JDATA) == 2 )
       IF ( LLOK ) THEN
          IF(obs_bodyElem_i(obsSpaceData,OBS_VNM,JDATA) /= BUFR_NEDZ ) THEN
             ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,JDATA)
          ELSE
             call utl_abort('oop_vobslyr: ZLEV cannot be set, BUFR_NEDZ not supported!')
          END IF
          IOBS = obs_bodyElem_i(obsSpaceData,OBS_HIND,JDATA)
          ITYP = obs_bodyElem_i(obsSpaceData,OBS_VNM,JDATA)
          if (bufr_IsAtmosConstituent(ITYP)) then
             varLevel = vnl_varLevelFromVarnum(ITYP, &
                        obs_headElem_i(obsSpaceData,OBS_CHM,IOBS))
          else
             varLevel = vnl_varLevelFromVarnum(ITYP)
          end if
          ZPT= col_getPressure(COLUMNGHR,1,IOBS,varLevel)
          ZPB= col_getPressure(COLUMNGHR,COL_GETNUMLEV(COLUMNGHR,varLevel),IOBS,varLevel)
          IF ( ZLEV < ZPT ) THEN
             call obs_bodySet_i(obsSpaceData,OBS_XTR,JDATA,1)
             !
             !- !!! WARNING !!! This obs is above the model lid. 
             !  We must turn off its assimilation flag  because the
             !  current obs operators cannot deal with this situation (JFC)                  
             if(varLevel /= 'SF') then
                write(*,*) 'oop_vobslyrs: Rejecting OBS above model lid, pressure = ', ZLEV,' < ',ZPT
                call obs_bodySet_i(obsSpaceData,OBS_ASS,JDATA, 0)
             end if
          ELSE IF ( ZLEV > ZPB ) THEN
             call obs_bodySet_i(obsSpaceData,OBS_XTR,JDATA,2)
          ELSE
             call obs_bodySet_i(obsSpaceData,OBS_XTR,JDATA,0)
          END IF
       END IF
    END DO
!$OMP END PARALLEL DO
    !
    !     1.2 ZZZ Vertical coordinate
    !
!$OMP PARALLEL DO PRIVATE(jdata,llok,zlev,iobs,ityp,varLevel,zpt,zpb,nlev)
    DO JDATA= 1,obs_numbody(obsSpaceData)
       LLOK = (obs_bodyElem_i(obsSpaceData,OBS_ASS,JDATA) == 1 .AND. &
            obs_bodyElem_i(obsSpaceData,OBS_VCO,JDATA) == 1 )
       IF ( LLOK ) THEN
          IF(obs_bodyElem_i(obsSpaceData,OBS_VNM,JDATA) /= BUFR_NEDZ ) THEN
             ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,JDATA)
          ELSE
             call utl_abort('oop_vobslyr: ZLEV cannot be set, BUFR_NEDZ not supported!')
          END IF
          IOBS = obs_bodyElem_i(obsSpaceData,OBS_HIND,JDATA)
          ITYP = obs_bodyElem_i(obsSpaceData,OBS_VNM,JDATA)
          if (bufr_IsAtmosConstituent(ITYP)) then
             varLevel = vnl_varLevelFromVarnum(ITYP, &
                        obs_headElem_i(obsSpaceData,OBS_CHM,IOBS))
          else
             varLevel = vnl_varLevelFromVarnum(ITYP)
          end if
          if(varLevel == 'SF') then
             ZPT= col_getHeight(columnghr,1,IOBS,'TH')/RG
             ZPB= col_getGZsfc(columnghr,IOBS)/RG                 
          else
             nlev=col_getNumLev(columnghr,varLevel)
             ZPT= col_getHeight(columnghr,1,IOBS,varLevel)/RG
             ZPB= col_getHeight(columnghr,NLEV,IOBS,varLevel)/RG
          end if
          IF ( ZLEV > ZPT ) THEN
             call obs_bodySet_i(obsSpaceData,OBS_XTR,JDATA,1)
             write(*,*) 'oop_vobslyrs: Rejecting OBS above model lid, height =', ZLEV,' > ',ZPT
             call obs_bodySet_i(obsSpaceData,OBS_ASS,JDATA, 0)
          ELSE IF ( ZLEV < ZPB ) THEN
             call obs_bodySet_i(obsSpaceData,OBS_XTR,JDATA,2)
          ELSE
             call obs_bodySet_i(obsSpaceData,OBS_XTR,JDATA,0)
          END IF
       END IF
    END DO
!$OMP END PARALLEL DO
    !
    !
    !     2. FInd interpolation layer
    !        ------------------------
    !        (Model levels are assumed to be in increasing order in Mbs)
    !
    !     2.1  PPP Vertical coordinate
    !
!$OMP PARALLEL DO PRIVATE(jdata,llok,iobs,zlev,ityp,varLevel,ik,nlev,jk,zpt,zpb)
    DO JDATA=1,obs_numbody(obsSpaceData)
       call obs_bodySet_i(obsSpaceData,OBS_LYR,JDATA,0)
       LLOK = ( (obs_bodyElem_i(obsSpaceData,OBS_ASS,JDATA) == 1     .OR. &
            obs_bodyElem_i(obsSpaceData,OBS_ASS,JDATA) == -1) .AND. &
            obs_bodyElem_i(obsSpaceData,OBS_VCO,JDATA) == 2 )
       IF ( LLOK ) THEN
          IOBS = obs_bodyElem_i(obsSpaceData,OBS_HIND,JDATA)
          ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,JDATA)
          ITYP = obs_bodyElem_i(obsSpaceData,OBS_VNM,JDATA)
          if (bufr_IsAtmosConstituent(ITYP)) then
             varLevel = vnl_varLevelFromVarnum(ITYP, &
                        obs_headElem_i(obsSpaceData,OBS_CHM,IOBS))
          else
             varLevel = vnl_varLevelFromVarnum(ITYP)
          end if
          IK = 1
          nlev=COL_GETNUMLEV(COLUMNGHR,varLevel)
          DO JK = 2,NLEV - 1
             ZPT = col_getPressure(COLUMNGHR,JK,IOBS,varLevel)
             IF( ZLEV > ZPT ) IK = JK
          END DO
          ZPT = col_getPressure(COLUMNGHR,IK,IOBS,varLevel)
          ZPB = col_getPressure(COLUMNGHR,IK+1,IOBS,varLevel) 
          call obs_bodySet_i(obsSpaceData,OBS_LYR,JDATA, IK)
       END IF
    END DO
!$OMP END PARALLEL DO
    !
    !     2.2  ZZZ Vertical coordinate and surface observations
    !
!$OMP PARALLEL DO PRIVATE(jdata,llok,iobs,zlev,ityp,varLevel,ik,nlev,jk,zpt)
    DO JDATA= 1,obs_numbody(obsSpaceData)
       LLOK = ( (obs_bodyElem_i(obsSpaceData,OBS_ASS,JDATA) == 1     .OR. &
            obs_bodyElem_i(obsSpaceData,OBS_ASS,JDATA) == -1) .AND. &
            obs_bodyElem_i(obsSpaceData,OBS_VCO,JDATA) == 1 )
       IF ( LLOK ) THEN
          IOBS = obs_bodyElem_i(obsSpaceData,OBS_HIND,JDATA)
          ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,JDATA)
          ITYP = obs_bodyElem_i(obsSpaceData,OBS_VNM,JDATA)
          if (bufr_IsAtmosConstituent(ITYP)) then
             varLevel = vnl_varLevelFromVarnum(ITYP, &
                        obs_headElem_i(obsSpaceData,OBS_CHM,IOBS))
          else
             varLevel = vnl_varLevelFromVarnum(ITYP)
          end if
          IK = 1
          nlev=COL_GETNUMLEV(COLUMNGHR,varLevel)
          DO JK = 2,NLEV - 1
             ZPT = col_getHeight(columnghr,JK,IOBS,varLevel)/RG
             IF( ZLEV < ZPT ) IK = JK
          END DO
          IF ( ITYP == BUFR_NEPS .or. ITYP == BUFR_NEPN .or. &
               ITYP == BUFR_NEZD ) THEN
             ! for surface observations associated with surface analysis variables
             IK = 0
          ELSEIF ( ITYP == BUFR_NETS .or. ityp == BUFR_NESS .OR. &
               ITYP == BUFR_NEUS .or. ityp == BUFR_NEVS .OR. &
               ITYP == BUFR_NEHS) THEN
             ! for surface observations associated with NON-surface analysis variables
             IK = nlev - 1
          END IF
          call obs_bodySet_i(obsSpaceData,OBS_LYR,JDATA, IK)
       END IF
    END DO
!$OMP END PARALLEL DO
    !
  end subroutine oop_vobslyrs


  subroutine oop_ppp_nl(columnhr,obsSpaceData,jobs,cdfam)
    !
    !**s/r oop_ppp_nl - Computation of Jobs and y - H(x)
    !                 for pressure-level observations
    !
    !*    Purpose:  -Interpolate vertically columnhr to
    !                the pressure levels of the observations. Then compute Jobs.
    !                A linear interpolation in ln(p) is performed.
    !
    !Arguments
    !     jobs:  contribution to Jobs
    !     cdfam: family of obsservation
    !
    implicit none
    type(struct_columnData) :: columnhr
    type(struct_obs) :: obsSpaceData
    real(8) :: jobs
    character(len=*) :: cdfam

    integer :: headerIndex,bodyIndex,ilyr
    integer :: iass,ixtr,ivco,ivnm
    real(8) :: zvar,zoer
    real(8) :: zwb,zwt,zexp,zgamma,ztvg
    real(8) :: zlev,zpt,zpb,zomp
    real(8) :: columnVarB,columnVarT
    character(len=4) :: varName
    character(len=2) :: varLevel
    real(8),pointer :: col_ptr(:),col_ptr_tt(:),col_ptr_hu(:)
    !
    ! Temperature lapse rate for extrapolation of gz below model surface
    !
    Write(*,*) "Entering subroutine oop_ppp_nl"

    zgamma = 0.0065D0 / GRAV
    zexp = MPC_RGAS_DRY_AIR_R8*zgamma

    jobs = 0.d0

    call obs_set_current_body_list(obsSpaceData, cdfam)
    BODY: do
       bodyIndex = obs_getBodyIndex(obsSpaceData)
       if (bodyIndex < 0) exit BODY

       ! Only process pressure level observations flagged to be assimilated
       iass=obs_bodyElem_i (obsSpaceData,OBS_ASS,bodyIndex)
       ivco=obs_bodyElem_i (obsSpaceData,OBS_VCO,bodyIndex)
       if(iass /= 1 .or. ivco /= 2) cycle BODY

       ixtr=obs_bodyElem_i (obsSpaceData,OBS_XTR,bodyIndex)
       ivnm=obs_bodyElem_i (obsSpaceData,OBS_VNM,bodyIndex)
       zvar=obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex)
       zlev=obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
       zoer=obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex)
       headerIndex=obs_bodyElem_i (obsSpaceData,OBS_HIND,bodyIndex)

       if ( ixtr == 0 ) then

          ! Process all data within the domain of the model
          ilyr  =obs_bodyElem_i (obsSpaceData,OBS_LYR,bodyIndex)
          varName = vnl_varnameFromVarnum(ivnm)
          varLevel = vnl_varLevelFromVarnum(ivnm)
          zpt= col_getPressure(columnhr,ilyr  ,headerIndex,varLevel)
          zpb= col_getPressure(columnhr,ilyr+1,headerIndex,varLevel)
          zwb  = log(zlev/zpt)/log(zpb/zpt)
          zwt  = 1.d0 - zwb
          if(ivnm == bufr_nees) then
             col_ptr_hu=>col_getColumn(columnhr,headerIndex,'HU')
             col_ptr_tt=>col_getColumn(columnhr,headerIndex,'TT')
             !columnVarB=lqtoes(col_ptr_hu(ilyr+1),col_ptr_tt(ilyr+1),zpb)
             !columnVarT=lqtoes(col_ptr_hu(ilyr  ),col_ptr_tt(ilyr  ),zpt)
             columnVarB=hutoes(col_ptr_hu(ilyr+1),col_ptr_tt(ilyr+1),zpb)
             columnVarT=hutoes(col_ptr_hu(ilyr  ),col_ptr_tt(ilyr  ),zpt)
          else
             if(trim(varName) == 'GZ') then
                col_ptr=>col_getColumn(columnhr,headerIndex,varName,'TH')
             else
                col_ptr=>col_getColumn(columnhr,headerIndex,varName)
             end if
             columnVarB=col_ptr(ilyr+1)
             columnVarT=col_ptr(ilyr  )
          end if
          zomp = zvar-(zwb*columnVarB+zwt*columnVarT)
          jobs = jobs + zomp*zomp/(zoer*zoer)
          call obs_bodySet_r(obsSpaceData,OBS_OMP,bodyIndex,zomp)

       elseif (ixtr == 2) then

          ! Process only GZ that is data below model's orography
          if(ivnm == BUFR_NEGZ ) then
             !
             ! Forward nonlinear model for geopotential data below model's orography
             !
             !ztvg = (1.0d0 + MPC_DELTA_R8 * exp(col_getElem(columnhr,col_getNumLev(columnhr,'TH'),headerIndex,'HU')))*  &
             !     col_getElem(columnhr,col_getNumLev(columnhr,'TH'),headerIndex,'TT')
            ztvg = (1.0d0 + MPC_DELTA_R8 * col_getElem(columnhr,col_getNumLev(columnhr,'TH'),headerIndex,'HU'))*  &
                 col_getElem(columnhr,col_getNumLev(columnhr,'TH'),headerIndex,'TT')
             zomp = (  zvar - col_getGZsfc(columnhr,headerIndex) -  &
                  ztvg/zgamma*(1.D0-(zlev/col_getElem(columnhr,1,headerIndex,'P0'))**zexp))
             jobs = jobs + zomp*zomp/(zoer*zoer)
             call obs_bodySet_r(obsSpaceData,OBS_OMP,bodyIndex,zomp)
          end if

       end if

    end do body

    jobs = 0.5d0 * jobs

  end subroutine oop_ppp_nl


  subroutine oop_sfc_nl(columnhr,obsSpaceData,jobs,cdfam)
    !
    !**s/r oop_sfc_nl - Computation of Jo and the residuals to the observations
    !                 FOR SURFACE DATA (except ground-based GPS zenith delay)
    !
    !*    Purpose:  -Interpolate vertically the contents of columnhr to
    !                the pressure levels of the observations. Then
    !                compute Jo.
    !                A linear interpolation in ln(p) is performed.
    !
    !Arguments
    !     jobs   : contribution to Jo
    !     cdfam  : family of observation
    !
    implicit none
    type(struct_columnData) :: columnhr
    type(struct_obs) :: obsSpaceData
    real(8)          :: jobs
    character(len=*) :: cdfam

    integer :: ipb,ipt,ivnm,headerIndex,bodyIndex
    real(8) :: zvar,zcon,zexp,zgamma,ztvg
    real(8) :: zlev,zhhh,zgamaz,zslope,gzhr
    real(8) :: columnVarB
    character(len=2) :: varLevel
    !
    ! Temperature lapse rate for extrapolation of gz below model surface
    !

    Write(*,*) "Entering subroutine oop_sfc_nl"

    zgamma = 0.0065d0 / GRAV
    zexp = 1.0D0/(MPC_RGAS_DRY_AIR_R8*zgamma)

    jobs = 0.d0

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
          if(obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) /= 1 .or.  &
               obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= 1) cycle BODY

          ! only process this set of surface observations
          ivnm=obs_bodyElem_i (obsSpaceData,OBS_VNM,bodyIndex)
          if( ivnm /= BUFR_NETS .and. ivnm /= BUFR_NEPS .and.  &
               ivnm /= BUFR_NEUS .and. ivnm /= BUFR_NEVS .and.  &
               ivnm /= BUFR_NESS .and. ivnm /= BUFR_NEPN ) cycle BODY

          zvar = obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex)
          zlev = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
          zhhh = zlev * grav
          varLevel = vnl_varLevelFromVarnum(ivnm)

          if(ivnm == BUFR_NETS .or. ivnm == BUFR_NESS .or.  &
               ivnm == BUFR_NEUS .or. ivnm == BUFR_NEVS) then
             ! T2m,(T-TD)2m,US,VS
             ! In this section we always extrapolate linearly the trial
             ! field at the model surface to the height of the
             ! surface observation whether the observation is above or
             ! below the model surface.
             ! NOTE: For (T-TD)2m,US,VS we do a zero order extrapolation

             if(ivnm == BUFR_NETS) then
                zslope = zgamma
             else
                zslope = 0.0d0
             end if

             ipt  = col_getNumLev(COLUMNHR,varLevel)-1 + col_getOffsetFromVarno(columnhr,ivnm)
             ipb  = ipt + 1

             if(ivnm.eq.bufr_ness) then
                columnVarB=hutoes(col_getElem(columnhr,col_getNumLev(COLUMNHR,'TH'),headerIndex,'HU'),  &
                     col_getElem(columnhr,col_getNumLev(COLUMNHR,'TH'),headerIndex,'TT'),  &
                     col_getPressure(columnhr,col_getNumLev(COLUMNHR,'TH'),headerIndex,'TH'))
             else
                columnVarB=col_getElem(columnhr,ipb,headerIndex)
             end if
             gzhr=col_getHeight(columnhr,col_getNumLev(columnhr,varLevel),headerIndex,varLevel)
             call obs_bodySet_r(obsSpaceData,OBS_OMP,bodyIndex,  &
                  (zvar-columnVarB + zslope*(zhhh-gzhr)) )

          else if ( ivnm == BUFR_NEPS .or. ivnm == BUFR_NEPN ) then
            ! Surface (PS) & mean sea level (PN) pressure cases
            ! Background surface pressure are corrected for the height difference with the 
            ! observation. For mean sea level observation, the observation height = 0.

            ! 1) Temperature difference = lapse-rate (6.5 degree/km) * height difference (dz)
            zgamaz = zgamma*(zhhh-col_getGZsfc(columnhr,headerIndex))

            ! 2) Compute the 2m background virtual temperature: Tv = T*(1+0.608*HU)
            ztvg = (1.0d0 + MPC_DELTA_R8 *  &
                   exp(col_getElem(columnhr,col_getNumLev(columnhr,'TH'),headerIndex,'HU'))) *  &
                   col_getElem(columnhr,col_getNumLev(columnhr,'TH'),headerIndex,'TT')
            
            ! 3) Compute the temperature ratio 
            ! The legacy code says...
            zcon = ((ztvg-zgamaz)/ztvg)
            ! However, the U.S. Standard Atmosphere (1976, U.S. Government Printing Office, Washington, D.C., 1976*)
            ! at page 12 says...
            ! zcon = (ztvg/(ztvg+zgamaz))
            ! but the former was found to perform better (gives lower O-P values) than the latter by J-F Caron in 2018 

            ! 4) O-P, where P = P0 * zcon ** zexp (page 12 of the U.S. Standard Atmosphere, 1976, 
            !                                      U.S. Government Printing Office, Washington, D.C., 1976*)
            call obs_bodySet_r(obsSpaceData,OBS_OMP,bodyIndex,  &
                               zvar-(col_getElem(columnhr,1,headerIndex,'P0')*zcon**zexp))

            ! (*) available at https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770009539_1977009539.pdf
          end if

          ! contribution to jobs
          jobs = jobs +(obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)*   &
               obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)) / &
               (obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex)*   &
               obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex))

       end do BODY

    end do HEADER

    jobs = 0.5d0 * jobs

  end subroutine oop_sfc_nl


  subroutine oop_zzz_nl(columnhr,obsSpaceData,jobs,cdfam)
    !
    !**s/r oop_zzz_nl - Computation of Jo and the residuals to the observations
    !                 FOR UPPER AIR DATAFILES
    !
    !Author  :  J. St-James, CMDA/SMC July 2003
    !
    !Revision :
    !
    !     Purpose:  - Interpolate vertically the contents of commvo
    !                 onto the heights (in meters) of the observations.
    !                 Compute Jo.
    !                 A linear interpolation in z is performed.
    !
    !Arguments
    !     jobs:  CONTRIBUTION to Jo
    !     cdfam: FAMILY OF OBSSERVATION
    !
    implicit none
    type(struct_columnData) :: columnhr
    type(struct_obs) :: obsSpaceData
    real(8) :: jobs
    character(len=*) :: cdfam

    integer :: ipb,ipt,ivnm,ik,headerIndex,bodyIndex
    real(8) :: zvar,zwb,zwt,zlev,zpt,zpb
    character(len=2) :: varLevel, obsfam

    Write(*,*) "Entering subroutine oop_zzz_nl"

    jobs = 0.d0

    call obs_set_current_body_list(obsSpaceData, cdfam)
    BODY: do
       bodyIndex = obs_getBodyIndex(obsSpaceData)
       if (bodyIndex < 0) exit BODY

       ! Process all height-level data within the domain of the model
       if( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= 1 .or.  &
            obs_bodyElem_i(obsSpaceData,OBS_XTR,bodyIndex) /= 0 .or.  &
            obs_bodyElem_i(obsSpaceData,OBS_VCO,bodyIndex) /= 1 ) cycle BODY

       ! In case not specified, make sure only PR family is processed
       obsfam = obs_getFamily(obsSpaceData,bodyIndex=bodyIndex)
       if( obsfam /= 'PR' ) cycle BODY

       headerIndex = obs_bodyElem_i(obsSpaceData,OBS_HIND,bodyIndex)
       zvar = obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex)
       zlev = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
       ik   = obs_bodyElem_i(obsSpaceData,OBS_LYR,bodyIndex)
       ivnm = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
       ipt = ik + col_getOffsetFromVarno(columnhr,ivnm)
       ipb = ipt+1
       varLevel = vnl_varLevelFromVarnum(ivnm)
       zpt= col_getHeight(columnhr,ik  ,headerIndex,varLevel)/RG
       zpb= col_getHeight(columnhr,ik+1,headerIndex,varLevel)/RG
       zwb  = (zpt-zlev)/(zpt-zpb)
       zwt  = 1.d0 - zwb
       if(ivnm == bufr_nees) call utl_abort('oop_zzz_nl: CANNOT ASSIMILATE ES!!!')
       call obs_bodySet_r(obsSpaceData,OBS_OMP,bodyIndex,  &
            zvar-zwb*col_getElem(columnhr,ipb,headerIndex) &
            - zwt*col_getElem(columnhr,ipt,headerIndex))

       ! contribution to jobs
       jobs = jobs + obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)*   &
            obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex) /  &
            (obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex)*   &
            obs_bodyElem_r(obsSpaceData,OBS_OER,bodyIndex))

    end do BODY

    jobs = 0.5d0 * jobs

  end subroutine oop_zzz_nl


  subroutine oop_gpsro_nl(columnhr,obsSpaceData,beSilent,jobs)
    !
    !**s/r oop_gpsro_nl - Computation of Jo and the residuals to the GPSRO observations
    !
    !
    !Author  : J. M. Aparicio Jan 2004
    !          Adapted Nov 2012 for both refractivity and bending angle data
    !    -------------------
    !*    Purpose:
    !
    !Arguments
    !     jobs: total value of Jobs for GPSRO
    !
    implicit none

    type(struct_columnData) :: columnhr
    type(struct_obs)        :: obsSpaceData
    logical :: beSilent
    real(8) :: jobs

    real(8) :: pjob, pjo1
    real(8) :: zlat, lat, slat
    real(8) :: zlon, lon
    real(8) :: zazm, azm
    integer :: iazm, isat, iclf, jj
    real(8) :: rad, geo, rad1, wfgps
    real(8), allocatable :: zpp(:)
    real(8), allocatable :: zdp(:)
    real(8), allocatable :: ztt(:)
    real(8), allocatable :: zhu(:)
    real(8), allocatable :: zuu(:)
    real(8), allocatable :: zvv(:)
    real(8) :: zp0, zmt
    real(8) :: hnh1, zobs, zmhx, zoer, zinc
    integer index_header, idatyp, index_body
    integer jl, ngpslev, nwndlev, stat
    logical  assim, firstheader, ldsc
    integer :: nh, nh1
    type(gps_profile)           :: prf
    real(8)      , allocatable :: h   (:),azmv(:)
    type(gps_diff), allocatable :: rstv(:),rstvp(:),rstvm(:)

    write(*,*)'ENTER oop_gpsro_nl'
    !
    ! Initializations
    !
    ngpslev=col_getNumLev(columnhr,'TH')
    nwndlev=col_getNumLev(columnhr,'MM')
    allocate(zpp(ngpslev))
    allocate(zdp(ngpslev))
    allocate(ztt(ngpslev))
    allocate(zhu(ngpslev))
    allocate(zuu(ngpslev))
    allocate(zvv(ngpslev))

    allocate( h    (gpsro_maxprfsize) )
    allocate( azmv (gpsro_maxprfsize) )
    allocate( rstv (gpsro_maxprfsize) )
    !if (levelgpsro == 1) then
    !  allocate( rstvp(gpsro_maxprfsize) )
    !  allocate( rstvm(gpsro_maxprfsize) )
    !end if

    jobs = 0.0d0

    !
    ! Loop over all header indices of the 'RO' family:
    !
    call obs_set_current_header_list(obsSpaceData,'RO')
    firstheader = .true.

    HEADER: do
       index_header = obs_getHeaderIndex(obsSpaceData)
       if (index_header < 0) exit HEADER
       !
       ! Process only refractivity data (codtyp 169)
       !
       idatyp = obs_headElem_i(obsSpaceData,OBS_ITY,index_header)
       if ( idatyp /= 169 ) cycle HEADER
       !
       ! Scan for requested data values of the profile, and count them
       !
       assim = .false.
       nh = 0
       call obs_set_current_body_list(obsSpaceData, index_header)
       BODY: do 
          index_body = obs_getBodyIndex(obsSpaceData)
          if (index_body < 0) exit BODY
          if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,index_body) == 1 ) then
             assim = .true.
             nh = nh + 1
          end if
       end do BODY
       !
       ! If no assimilations are requested, skip to next header
       !
       if (.not.assim) cycle HEADER
       !
       ! Basic geometric variables of the profile:
       !
       iazm = obs_headElem_i(obsSpaceData,OBS_AZA,index_header)
       isat = obs_headElem_i(obsSpaceData,OBS_SAT,index_header)
       iclf = obs_headElem_i(obsSpaceData,OBS_ROQF,index_header)
       rad  = obs_headElem_r(obsSpaceData,OBS_TRAD,index_header)
       geo  = obs_headElem_r(obsSpaceData,OBS_GEOI,index_header)
       zazm = 0.01d0*iazm / MPC_DEGREES_PER_RADIAN_R8
       zmt  = col_getGZsfc(columnhr,index_header)/RG
       wfgps=0.d0
       do jj=1,numgpssats
          if (isat == igpssat(jj)) wfgps=wgps(jj)
       end do
       !
       ! Profile at the observation location:
       !
       zlat = obs_headElem_r(obsSpaceData,OBS_LAT,index_header)
       zlon = obs_headElem_r(obsSpaceData,OBS_LON,index_header)
       lat  = zlat * MPC_DEGREES_PER_RADIAN_R8
       lon  = zlon * MPC_DEGREES_PER_RADIAN_R8
       azm  = zazm * MPC_DEGREES_PER_RADIAN_R8
       slat = sin(zlat)
       zmt  = zmt * RG / gpsgravitysrf(slat)
       zp0  = col_getElem(columnhr,1,index_header,'P0')
       do jl = 1, ngpslev
          !
          ! Profile x
          !
          zpp(jl) = col_getPressure(columnhr,jl,index_header,'TH')
          zdp(jl) = col_getPressureDeriv(columnhr,jl,index_header,'TH')
          ztt(jl) = col_getElem(columnhr,jl,index_header,'TT') - p_tc
          zhu(jl) = col_getElem(columnhr,jl,index_header,'HU')
       end do

       if((col_getPressure(columnhr,1,index_header,'TH') + 1.0d-4) < &
            col_getPressure(columnhr,1,index_header,'MM')) then
          ! case with top thermo level above top momentum level (Vcode=5002)
          do jl = 1, nwndlev
             zuu(jl) = col_getElem(columnhr,jl  ,index_header,'UU') * p_knot
             zvv(jl) = col_getElem(columnhr,jl  ,index_header,'VV') * p_knot
          end do
       else
          ! case without top thermo above top momentum level or unstaggered (Vcode=5001/4/5)
          do jl = 1, nwndlev-1
             zuu(jl) = col_getElem(columnhr,jl+1,index_header,'UU') * p_knot
             zvv(jl) = col_getElem(columnhr,jl+1,index_header,'VV') * p_knot
          end do
          zuu(nwndlev) = zuu(nwndlev-1)
          zvv(nwndlev) = zuu(nwndlev-1)
       end if
       zuu(ngpslev) = zuu(nwndlev)
       zvv(ngpslev) = zuu(nwndlev)
       !     
       ! GPS profile structure:
       !
       call gps_struct1sw(ngpslev,zLat,zLon,zAzm,zMT,Rad,geo,zP0,zPP,zDP,zTT,zHU,zUU,zVV,prf)
       ldsc=.not.btest(iclf,16-3)
       !
       ! Prepare the vector of all the observations:
       !
       nh1 = 0
       !
       ! Loop over all body indices for this index_header:
       ! (start at the beginning of the list)
       !
       call obs_set_current_body_list(obsSpaceData, index_header)
       BODY_2: do 
          index_body = obs_getBodyIndex(obsSpaceData)
          if (index_body < 0) exit BODY_2
          IF ( obs_bodyElem_i(obsSpaceData,OBS_ASS,index_body) == 1 ) then
             nh1      = nh1 + 1
             h(nh1)   = obs_bodyElem_r(obsSpaceData,OBS_PPP,index_body)
             azmv(nh1)= zazm
          end if
       end do BODY_2
       !
       ! Apply the observation operator:
       !
       if (levelgpsro == 1) then
          call gps_bndopv1(h      , azmv, nh, prf, rstv)
          !call gps_bndopv1(h+wfgps, azmv, nh, prf, rstvp)
          !call gps_bndopv1(h-wfgps, azmv, nh, prf, rstvm)
          !do nh1 = 1, nh
          !  rstv(nh1)=(rstvp(nh1)+rstv(nh1)+rstvm(nh1))/3.d0
          !end do
       else
          call gps_refopv (h,       nh, prf, rstv)
       end if
       !
       ! Perform the (H(x)-Y)/S operation:
       !
       nh1 = 0
       pjob = 0.d0
       !
       ! Loop over all body indices for this index_header:
       ! (start at the beginning of the list)
       !
       call obs_set_current_body_list(obsSpaceData, index_header)
       BODY_3: do 
          index_body = obs_getBodyIndex(obsSpaceData)
          if (index_body < 0) exit BODY_3
          IF ( obs_bodyElem_i(obsSpaceData,OBS_ASS,index_body) == 1 ) then
             nh1 = nh1 + 1
             !
             ! Altitude:
             !
             hnh1= obs_bodyElem_r(obsSpaceData,OBS_PPP,index_body)
             if (levelgpsro == 1) hnh1=hnh1-rad
             !
             ! Observation operator H(x)
             !
             zmhx = rstv(nh1)%var
             !
             ! Observation value    Y
             !
             zobs = obs_bodyElem_r(obsSpaceData,OBS_VAR,index_body)
             !
             ! Observation error    S
             !
             zoer = obs_bodyElem_r(obsSpaceData,OBS_OER,index_body)
             !
             ! Normalized increment
             !
             zinc = (zmhx - zobs) / zoer
             !                           
             ! Datum contribution to Jo:
             !
             pjo1 = 0.5d0 * zinc * zinc
             !
             ! Total (PJO) and per profile (PJOB) cumulatives:
             !
             jobs = jobs + pjo1
             pjob= pjob+ pjo1
             !
             if (firstheader .and. .not.beSilent) then
                write(*,  &
                     '(A9,i10,3f7.2,f11.1,4f12.6,15f12.4)') 'DOBSGPSRO',  &
                     index_header,lat,lon,azm,hnh1,zobs,zoer,  &
                     zmhx,zinc,pjob,prf%gst(ngpslev)%var  
             end if
             call obs_bodySet_r(obsSpaceData,OBS_OMP,index_body, zobs - zmhx)
          end if
       end do BODY_3

       if ( .not.beSilent ) write(*,'(A9,i10,2f7.2,f18.10,f12.4,2I6)')  &
            'GPSRO_JO',index_header,lat,lon,pjob,zmt,isat,ldsc
       firstheader = .false.
    end do HEADER

    !if (levelgpsro == 1) then
    !  deallocate( rstvm )
    !  deallocate( rstvp )
    !end if
    deallocate( rstv )
    deallocate( azmv )
    deallocate( h    )

    deallocate(zvv)
    deallocate(zuu)
    deallocate(zhu)
    deallocate(ztt)
    deallocate(zdp)
    deallocate(zpp)

    write(*,*)'EXIT oop_gpsro_nl'

  end subroutine oop_gpsro_nl


  subroutine oop_gpsgb_nl(columnhr,obsSpaceData,beSilent,jobs,analysisMode_opt)
    !
    !**s/r oop_gpsgb_nl - Computation of Jo and the residuals to the GB-GPS ZTD observations
    !
    !
    !Author  : S. Macpherson  ARMA/MRD
    !Revisions:
    !          S. Macpherson Oct 2012
    !           -- conversion of 3dvar v11.2.2 version to Rev189 modular form.
    !           -- uses new (modified) GPS-RO modgps*.f90 for ZTD observation operator
    !           -- option to use old NL operator removed
    !           -- ZTD operator gps_ztdopv is found in MODIF modgps08refop.f90
    !           -- Uses columnData_mod.
    !
    !          S. Macpherson Dec 2012 - Jan 2013
    !           -- update from Rev189 to Rev213
    !           -- new namelist parameters in modgpsztd_mod
    !           -- ZTD operator gps_ztdopv is found in NEW modgps08ztdop.cdk90
    !           -- ZETA (eta/hybrid values) and ZGZ profiles no longer needed.
    !           -- add filter for 1-OBS option (L1OBS=.true. in namelist)
    !           -- Set vGPSZTD_Index(numGPSZTD) for Jacobian storage
    !
    !          S. Macpherson Jun 2013
    !           -- Use true implementation of ZDP (dP/dP0), although not needed here
    !
    !          S. Macpherson Nov 2014
    !           -- modifications for case where P(nlev) is not equal to P0
    !
    !          S. Macpherson Jan 2015
    !           -- adadpt for E-GVAP data (assimilate ZTD without surface met data, i.e. Psfc)
    !
    !Arguments (out)
    !     jobs: total value of Jo for all GB-GPS (ZTD) observations
    !
    implicit none

    type(struct_columnData) :: columnhr
    type(struct_obs) :: obsSpaceData
    logical           :: beSilent
    real(8)           :: jobs
    logical, optional :: analysisMode_opt

    real(8), allocatable :: zpp (:)
    real(8), allocatable :: zdp (:)
    real(8), allocatable :: ztt (:)
    real(8), allocatable :: zhu (:)
    real(8) :: zlat, lat, zlon, lon
    real(8) :: zp0, zmt, zdzmin
    real(8) :: zobs, zoer, zinc, zhx, zlev
    real(8) :: zdz, zpsobs, zpsmod, zpwmod, zpomp, zpomps
    real(8) :: ztdomp(max_gps_data)
    real(8) :: bias, std
    integer :: headerIndex, bodyIndex, ioneobs, idatyp, ityp, index_ztd, iztd
    integer :: jl, nlev_T, nobs2p, stat
    integer :: icount1, icount2, icount3, icount, icountp
    logical  :: assim, llrej, analysisMode, lfsl
    character(9) :: cstnid
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

    write(*,*)'ENTER oop_gpsgb_nl'

    if(present(analysisMode_opt)) then
       analysisMode = analysisMode_opt
    else
       analysisMode = .true.
    end if

    zpomps = 0.0d0

    ! Ensure Jacobian-related arrays are not allocated to force them to be recalculated in oop_H
    if(allocated(vGPSZTD_Jacobian)) deallocate(vGPSZTD_Jacobian)
    if(allocated(vGPSZTD_lJac)) deallocate(vGPSZTD_lJac)

    zdzmin = dzmin      
    nobs2p = 50
    jobs = 0.d0

    nlev_T = col_getNumLev(columnhr,'TH')
    if (ltestop .and. .not.beSilent) write(*,*) '  col_getNumLev(columnhr,TH) = ',nlev_T

    !
    ! Initializations
    !
    allocate(ztt(nlev_T))
    allocate(zhu(nlev_T))
    allocate(zdp(nlev_T))
    allocate(zpp(nlev_T))

    if ( .not.beSilent ) then
      write(*, *) ' '
      write(*, *) ' '
      write(*,'(A11,A9,3A8,A9,4A8,2A9,A7,A10,A11)')  &
           'OOP_GPSGB_NL','CSTNID','ZLAT','ZLON','ZLEV','ZDZ','ZOBS','ZOER','ZHX','O-P',  &
           'ZPOMPS','ZPOMP','ZPWMOD','Jobs','ZINC2'
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
          ityp = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
          if ( (ityp == BUFR_NEZD) .and. &
               (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == 1) ) then
             zlev = obs_bodyElem_r(obsSpaceData,OBS_PPP,bodyIndex)
             assim = .true.
             ! Index in body of ZTD datum (assume at most 1 per header)
             index_ztd = bodyIndex
             icount = icount + 1
          end if
          if ( ityp == bufr_neps ) then
             if ( (obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == 1) .or. llblmet ) then
                zpsobs = obs_bodyElem_r(obsSpaceData,OBS_VAR,bodyIndex)
                zpomps = obs_bodyElem_r(obsSpaceData,OBS_OMP,bodyIndex)
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
       zmt  = col_getGZsfc(columnhr,headerIndex)/RG
       zp0  = col_getElem(columnhr,1,headerIndex,'P0')
       do jl = 1, nlev_T
          zpp(jl) = col_getPressure(columnhr,jl,headerIndex,'TH')
          ! True implementation of ZDP (dP/dP0)
          zdp(jl) = col_getPressureDeriv(columnhr,jl,headerIndex,'TH')
          ztt(jl) = col_getElem(columnhr,jl,headerIndex,'TT')-MPC_K_C_DEGREE_OFFSET_R8
          zhu(jl) = col_getElem(columnhr,jl,headerIndex,'HU')
       end do
       zdz = zlev - zmt

       ! Fill GPS ZTD profile structure (PRF):
       call gps_structztd(nlev_T,lat,lon,zmt,zp0,zpp,zdp,ztt,zhu,lbevis,irefopt,prf)

       ! Apply the GPS ZTD observation operator
       ! --> output is model ZTD (type gps_diff) and P at obs height ZLEV
       call gps_ztdopv(zlev,prf,lbevis,zdzmin,ztdopv,zpsmod,iztdop)

       ! Get model profile PW
       call gps_pw(prf,zpwmod)
       ! ZTD (m)
       zhx    = ztdopv%var

       ! If analysis mode, reject ZTD data for any of the following conditions:
       !    (1) the trial PW is too low (extremely dry) 
       !    and if LASSMET=true and for NOAA/FSL sites only:
       !      (2) Ps observation is missing or out of normal range
       !      (3) the ABS(Ps(obs)-Ps(mod)) difference is too large
       llrej = .false.
       zpomp = -9999.0D0
       if ( analysisMode ) then
          llrej = ( zpwmod < zpwmin )
          if ( lassmet .and. lfsl ) then
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
          call obs_bodySet_i(obsSpaceData,OBS_ASS,index_ztd, 0)
          if ( .not. lassmet ) icount1 = icount1 + 1
       end if

       ! Perform the (H(x)-Y)/SDERR operation
       !
       ! loop over all body indices for this headerIndex
       call obs_set_current_body_list(obsSpaceData, headerIndex)
       BODY_2: do 
          bodyIndex = obs_getBodyIndex(obsSpaceData)
          if (bodyIndex < 0) exit BODY_2
          ityp = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
          if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == 1 .and.  &
               ityp == BUFR_NEZD ) then
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
             call obs_bodySet_r(obsSpaceData,OBS_OMP,bodyIndex, zobs - zhx)

             jobs = jobs + 0.5d0 * zinc * zinc
             !
             ! Apply data selection criteria for 1-OBS Mode
             !
             if ( l1obs .and. ioneobs == -1 ) then
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
                     'OOP_GPSGB_NL: ',cstnid,zlat,zlon,zlev,zdz,zobs,zoer/yzderrwgt,zhx,-zinc*zoer,  &
                     zpomps/100.d0,zpomp/100.d0,zpwmod,jobs,zinc/zoer
             end if

          end if

       end do BODY_2

    end do HEADER

    deallocate(ztt)
    deallocate(zhu)
    deallocate(zdp)
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

    if ( l1obs .and. analysisMode ) then
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
                   call obs_bodySet_i(obsSpaceData,OBS_ASS,bodyIndex, 0)
                end do BODY_1
             end if
          end do HEADER_1
       else
          call utl_abort('ERROR: FAILED TO SELECT SINGLE OBSERVATION!')
       end if
    end if

    numgpsztd = icountp

    if ( analysisMode .and. icount > 0 .and. .not.l1obs .and. .not.beSilent ) then
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
       if ( icountp > 0 ) then
          write(*, *) 'MEAN Jo = (jobs/numGPSZTD)*YZDERRWGT**2 = ',(jobs/real(icountp,8))*yzderrwgt**2
       end if
       write(*,*) ' '
    end if

    if ( icount > 0 .and. numgpsztd > 0) then

       if ( analysisMode ) then
          if ( .not.beSilent ) write(*,*) ' Number of GPS ZTD data to be assimilated (numGPSZTD) = ', numgpsztd
       else
          if ( .not.beSilent ) write(*,*) ' Number of GPS ZTD data for background check (numGPSZTD) = ', numgpsztd
       end if

       if ( .not.beSilent ) write(*,*) ' Allocating and setting vGPSZTD_Index(numGPSZTD)...'
       if(allocated(vgpsztd_index)) deallocate(vgpsztd_index)
       allocate(vgpsztd_index(numgpsztd))
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
                ityp = obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex)
                if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == 1 .and.  &
                     ityp == BUFR_NEZD ) then  
                   iztd = iztd + 1
                   vgpsztd_index(iztd) = headerIndex
                end if
             end do BODY_3
          end if
       end do HEADER_2

       if ( iztd /= numgpsztd ) then
          call utl_abort('ERROR: vGPSZTD_Index init: iztd /= numGPSZTD!')
       end if

    end if

    write(*,*)'EXIT oop_gpsgb_nl'

  end subroutine oop_gpsgb_nl


  subroutine oop_tovs_nl(columnghr,obsSpaceData,datestamp,limlvhu,beSilent,  &
                         jobs,bgckMode_opt,option_opt,sourceObs_opt,destObs_opt)
    !
    !**s/r oop_tovs_nl  - Computation of jobs and the residuals to the tovs observations
    !
    !
    !author        : j. halle *cmda/aes  april 8, 2005
    !
    !arguments
    !     option_opt: defines input state:
    !               'HR': High Resolution background state,
    !               'LR': Low  Resolution background state, (CURRENTLY NOT SUPPORTED)
    !               'MO': Model state. (CURRENTLY NOT SUPPORTED)
    !     jobs: total value of jobs for tovs
    !
    implicit none

    type(struct_columnData) :: columnghr
    type(struct_obs) :: obsSpaceData
    integer :: datestamp
    real(8) :: limlvhu
    logical :: beSilent
    real(8) :: jobs
    logical, optional :: bgckMode_opt
    character(len=*), optional :: option_opt        ! only valid value is HR
    integer, optional, intent(in) :: sourceObs_opt ! usually set to OBS_VAR
    integer, optional, intent(in) :: destObs_opt   ! usually set to OBS_OMP

    integer :: jdata, sourceObs, destObs
    logical :: llprint,bgckMode
    character(len=2) :: option

    if (.not.obs_famExist(obsSpaceData,'TO',local_mpi=.true.)) then
       jobs=0.0d0
       return
    end if

    ! 0. set default values if bgckMode, option and source/dest columns not specified
    !

    Write(*,*) "Entering subroutine oop_tovs_nl"

    if(present(bgckMode_opt)) then
       bgckMode = bgckMode_opt
    else
       bgckMode = .false.
    end if

    if(present(option_opt)) then
       option = option_opt(1:2)
    else
       option = 'HR'
    end if
    if ( option /= 'HR' ) call utl_abort('oop_tovs_nl: Invalid option for input state')

    if(present(sourceObs_opt)) then
       sourceObs = sourceObs_opt
    else
       sourceObs = OBS_VAR
    end if

    if(present(destObs_opt)) then
       destObs = destObs_opt
    else
       destObs = OBS_OMP
    end if

    ! 1.   Prepare atmospheric profiles for all tovs observation points for use in rttov
    ! .    -----------------------------------------------------------------------------
    call tvs_fillProfiles(columnghr,obsSpaceData,datestamp,limlvhu,bgckMode,beSilent)
    if ( .not.beSilent ) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! 2.   Compute radiance
    ! .    ----------------
    call tvs_rttov(columnghr,obsSpaceData,bgckMode,beSilent)
    if ( .not.beSilent ) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! 3.   Compute Jobs and the residuals
    ! .    ----------------------------
    if ( option == 'HR' .or. option == 'LR' ) then
       do jdata=1,obs_numbody(obsSpaceData)
          call obs_bodySet_r(obsSpaceData,OBS_PRM,jdata, obs_bodyElem_r(obsSpaceData,sourceObs,jdata))
       end do
    end if

    if(option == 'HR') then
       llprint = .true.
    else
       llprint = .false.
    end if
    jobs = 0.0d0

    if ( beSilent ) llprint = .false.
    call tvs_calc_jo(jobs,llprint,obsSpaceData,destObs)

  end subroutine oop_tovs_nl


  subroutine oop_chm_nl(columnhr,obsSpaceData,jobs)
    !
    !**s/r oop_chm_nl - Computation of Jo and the residuals to the observations
    !                 for all observations of the CH (chemical constituents) family.
    !                 Stores OmP in OBS_OMP in obsSpaceData.
    !
    ! Author: M. Sitwell, ARQI/AQRD, Aug 2015
    !         - Reduced to a call of chm_observation_operators
    !           which contains much of the original version by Ping Du (CMDA/MSC) and
    !           Y. Rochon (ARQI/AQRD), Jan 2015 (partially based on corresponding 
    !           pre-EnVar routine by Y.J. Rochon and Y. Yang, July 2005 to Feb 2013),
    !           with further original changes by Y. Rochon and M. Sitwell, 2015.
    !
    ! Revision:
    ! 
    ! Purpose:  Computation of Jo and the residuals to the observations
    !           for all observations of the CH (chemical constituents) family.
    !           The array columnhr contains the input model array.
    !
    !
    ! Arguments
    !
    !   Input
    !
    !     columnhr      : Columnar arrays from model fields
    !     obsSpaceData  : Obs space data structure
    !
    !   Output
    !
    !     jobs  : contribution to Jo
    !
    ! Comments:
    !
    !-----------------------------------------------------------------------------------

    implicit none
    
    type(struct_columnData) :: columnhr
    type(struct_obs) :: obsSpaceData
    real(8) :: jobs
    
    if (.not.obs_famExist(obsSpaceData,'CH',local_mpi=.true.)) then
       jobs = 0.0d0
       return
    end if

    call chm_observation_operators(columnhr,obsSpaceData,kmode=0,jobs_opt=jobs) ! kmode=0 for general operator

  end subroutine oop_chm_nl

!--------------------------------------------------------------------------
!! *Purpose*: Compute simulated observations from profiled model increments.
!!            It returns Hdx in OBS_WORK. Calls the several linear observation operators.
!!
!! Input
!!
!!v     obsAssVal_opt  value of OBS_ASS in obsSpaceData that Htl will be calculated for (optional, default = 1)
!!
!! Revisions
!!v      M. Sitwell, Feb 2017
!!v        - Added optional obsAssVal_opt input argument
!--------------------------------------------------------------------------
  subroutine oop_Htl(column,columng,obsSpaceData,min_nsim,obsAssVal_opt)
    implicit none

    type(struct_columnData) :: column,columng
    type(struct_obs) :: obsSpaceData
    type(struct_vco), pointer :: vco_anl
    integer, intent(in) :: min_nsim
    integer, intent(in), optional :: obsAssVal_opt
    integer :: obsAssVal
    logical, save :: firstTime = .true.

    IF(mpi_myid == 0) THEN
       write(*,*)'OOP_Htl - Linearized observation operators'
    end if

    if (present(obsAssVal_opt)) then
       obsAssVal = obsAssVal_opt
    else
       obsAssVal = 1
    end if

    vco_anl => col_getVco(columng)

    if ( firstTime ) then
      !     Find interpolation layer in model profiles (used by several operators)
      if ( col_getNumLev(columng,'MM') > 1 ) call oop_vobslyrs(columng,obsSpaceData)

      !     Initialize some operators needed by linearized H
      call subasic_obs(columng)

      firstTime = .false.
    end if

    call tmg_start(42,'OBS_PPP_TLAD')
    call oop_Hpp(obsAssVal)           ! fill in OBS_WORK : Hdx
    call tmg_stop(42)

    call tmg_start(43,'OBS_SFC_TLAD')
    call oop_Hsf(obsAssVal)           ! fill in OBS_WORK : Hdx
    call tmg_stop (43)

    call tmg_start(44,'OBS_TOV_TLAD')
    call oop_Hto(obsAssVal)           ! fill in OBS_WORK : Hdx
    call tmg_stop (44)

    call tmg_start(45,'OBS_GPSRO_TLAD')
    call oop_Hro(obsAssVal)
    call tmg_stop (45)

    call tmg_start(46,'OBS_ZZZ_TLAD')
    call oop_Hzp(obsAssVal)
    call tmg_stop (46)

    call tmg_start(47,'OBS_GPSGB_TLAD')
    if (numGPSZTD > 0)  call oop_Hgp(obsAssVal)
    call tmg_stop (47)

    call tmg_start(126,'OBS_CHM_TL')
    call oop_Hchm(obsAssVal)          ! fill in OBS_WORK : Hdx
    call tmg_stop (126)


  CONTAINS

    subroutine subasic_obs(columng)
      implicit none

      ! s/r SUBASIC_OBS
      !     OBJECT: Initialise background state dependant factors
      !             and vectors for use in TLM and adjoint of
      !             non-linear operator
      !
      !     Author  : S. Pellerin *ARMA/AES Sept. 98
      !
      !
      type(struct_columnData) :: columng
      type(struct_vco), pointer :: vco_anl
      integer :: jlev,columnIndex,nlev_T,vcode_anl,status
      real(8) :: zhu,one

      if ( .not.col_varExist('TT') .or. .not.col_varExist('HU') ) return

      write(*,*) 'subasic_obs: setting up linearized Tv operator'

      vco_anl => col_getVco(columng)
      one=1.0D0
      nlev_T = col_getNumLev(columng,'TH')
      status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)

      if( Vcode_anl /= 5002 .and. Vcode_anl /= 5005 ) then
         call utl_abort('subasic_obs: invalid vertical coord!')
      end if

      ! initialize virtual temperature operator

!$OMP PARALLEL DO PRIVATE(jlev,columnIndex,zhu)
      do jlev = 1, nlev_T
         do columnIndex=1,col_getNumCol(columng)
            zhu=col_getElem(columng,jlev,jobs,'HU')
            columng%oltv(1,jlev,jobs) = fottva(zhu,one)
            columng%oltv(2,jlev,jobs) = folnqva(zhu,col_getElem(columng,  &
                 jlev,jobs,'TT'),one)
         end do
      end do
!$OMP END PARALLEL DO

    end subroutine subasic_obs


    SUBROUTINE oop_Hpp(obsAssVal)
      !*
      !* Purpose: Compute simulated Upper Air observations from profiled model
      !*          increments.
      !*          It returns Hdx in OBS_WORK
      !*          Interpolate vertically the contents of commvo to
      !*          the pressure levels of the observations.
      !*          A linear interpolation in ln(p) is performed.
      !*
      !*implicits

      implicit none

      integer, intent(in) :: obsAssVal

      INTEGER IPB,IPT
      INTEGER INDEX_HEADER,INDEX_FAMILY,IK
      INTEGER J,INDEX_BODY,ITYP,nlev_T
      REAL*8 ZDADPS,ZCON
      REAL*8 ZWB,ZWT,ZLTV,ZTVG
      REAL*8 ZLEV,ZPT,ZPB
      REAL*8 dPdPsT,dPdPsB
      REAL*8 columnVarB,columnVarT,columngVarB,columngVarT
      LOGICAL LLASSIM,LLDIAG
      INTEGER, PARAMETER :: numFamily=3
      CHARACTER(len=2) :: list_family(numFamily),varLevel

      list_family(1) = 'UA'
      list_family(2) = 'AI'
      list_family(3) = 'SW'

      FAMILY: do index_family=1,numFamily

         call obs_set_current_body_list(obsSpaceData,list_family(index_family))
         BODY: DO 
            index_body = obs_getBodyIndex(obsSpaceData)
            if (index_body < 0) exit BODY

            llassim= (obs_bodyElem_i(obsSpaceData,OBS_ASS,index_body) == obsAssVal) &
                 .AND. (obs_bodyElem_i(obsSpaceData,OBS_XTR,index_body) == 0) &
                 .AND. (obs_bodyElem_i(obsSpaceData,OBS_VCO,index_body) == 2)
            lldiag = (obs_bodyElem_i(obsSpaceData,OBS_ASS,index_body) == -1) &
                 .AND. (obs_bodyElem_i(obsSpaceData,OBS_VCO,index_body) == 2)
            IF (llassim .or. lldiag) THEN
               index_header = obs_bodyElem_i(obsSpaceData,OBS_HIND,INDEX_BODY)
               ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY)
               ITYP = obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY)
               varLevel = vnl_varLevelFromVarnum(ityp)
               IK   = obs_bodyElem_i(obsSpaceData,OBS_LYR,INDEX_BODY)
               IPT  = IK + col_getOffsetFromVarno(columng,ityp)
               IPB  = IPT+1
               ZPT    = col_getPressure(COLUMNG,IK  ,INDEX_HEADER,varLevel)
               ZPB    = col_getPressure(COLUMNG,IK+1,INDEX_HEADER,varLevel)
               dPdPsT = col_getPressureDeriv(COLUMNG,IK  ,INDEX_HEADER,varLevel)
               dPdPsB = col_getPressureDeriv(COLUMNG,IK+1,INDEX_HEADER,varLevel)
               ZWB  = LOG(ZLEV/ZPT)/LOG(ZPB/ZPT)
               ZWT  = 1.0D0 - ZWB

               ZDADPS   = ( LOG(ZLEV/ZPB)*dPdPsT/ZPT -   &
                    LOG(ZLEV/ZPT)*dPdPsB/ZPB )/  &
                    LOG(ZPB/ZPT)**2

               if(ityp.eq.bufr_nees) then
                  columnVarB=hutoes_tl(col_getElem(column,IK+1,INDEX_HEADER,'HU'), &
                       col_getElem(column,IK+1,INDEX_HEADER,'TT'), &
                       col_getElem(column,1,INDEX_HEADER,'P0'), &
                       col_getElem(columng,IK+1,INDEX_HEADER,'HU'), &
                       col_getPressure(columng,IK+1,INDEX_HEADER,'TH'), &
                       dPdPsB)
                  columnVarT=hutoes_tl(col_getElem(column,IK  ,INDEX_HEADER,'HU'), &
                       col_getElem(column,IK  ,INDEX_HEADER,'TT'), &
                       col_getElem(column,1,INDEX_HEADER,'P0'), &
                       col_getElem(columng,IK  ,INDEX_HEADER,'HU'), &
                       col_getPressure(columng,IK  ,INDEX_HEADER,'TH'),  &
                       dPdPsT)
                  columngVarB=hutoes(col_getElem(columng,IK+1,INDEX_HEADER,'HU'), &
                       col_getElem(columng,IK+1,INDEX_HEADER,'TT'), &
                       col_getPressure(columng,IK+1,INDEX_HEADER,'TH'))
                  columngVarT=hutoes(col_getElem(columng,IK  ,INDEX_HEADER,'HU'), &
                       col_getElem(columng,IK  ,INDEX_HEADER,'TT'), &
                       col_getPressure(columng,IK  ,INDEX_HEADER,'TH'))
               else
                  columnVarB=col_getElem(column,IPB,INDEX_HEADER)
                  columnVarT=col_getElem(column,IPT,INDEX_HEADER)
                  columngVarB=col_getElem(columng,IPB,INDEX_HEADER)
                  columngVarT=col_getElem(columng,IPT,INDEX_HEADER)
               end if
               call obs_bodySet_r(obsSpaceData,OBS_WORK,INDEX_BODY,   &
                    ZWB*columnVarB + ZWT*columnVarT+  &
                    (columngVarB - columngVarT)*  &
                    ZDADPS*col_getElem(COLUMN,1,INDEX_HEADER,'P0'))

            end if

         end do BODY

      end do FAMILY

    end subroutine oop_Hpp


    SUBROUTINE oop_Hsf(obsAssVal)
      !*
      !* Purpose: Compute simulated surface observations from profiled model
      !*          increments.
      !*          It returns Hdx in OBS_WORK
      !*
      IMPLICIT NONE

      integer, intent(in) :: obsAssVal

      INTEGER IPB,IPT,IXTR
      INTEGER INDEX_HEADER,IK
      INTEGER J,INDEX_BODY,ITYP,INDEX_FAMILY,nlev
      REAL*8 ZCON
      REAL*8 ZWB,ZWT, ZEXP,ZGAMMA,ZLTV,ZTVG
      REAL*8 ZLEV,ZPT,ZPB,ZDELPS,ZDELTV,ZGAMAZ,ZHHH
      REAL*8 columnVarB
      REAL*8 dPdPsfc
      INTEGER, PARAMETER :: numFamily=4
      CHARACTER(len=2) :: list_family(numFamily),varLevel
      !C
      !C     Temperature lapse rate for extrapolation of gz below model surface
      !C
      zgamma = 0.0065d0 / GRAV
      zexp   = 1.0d0/(MPC_RGAS_DRY_AIR_R8*zgamma)
      !C
      !C
      list_family(1) = 'UA'
      list_family(2) = 'SF'
      list_family(3) = 'SC'
      list_family(4) = 'GP'

      FAMILY: do index_family=1,numFamily

         call obs_set_current_body_list(obsSpaceData, list_family(index_family))
         BODY: do
            index_body = obs_getBodyIndex(obsSpaceData)
            if (index_body < 0) exit BODY

            ! Process all data within the domain of the model
            ityp = obs_bodyElem_i(obsSpaceData,OBS_VNM,index_body)
            if ( ityp == bufr_nezd ) cycle BODY
            if(    (obs_bodyElem_i(obsSpaceData,OBS_VCO,index_body) == 1) &
                 .and. (obs_bodyElem_i(obsSpaceData,OBS_ASS,index_body) == obsAssVal) &
                 .and. (ityp == bufr_nets .or. ityp == bufr_neps  &
                 .or. ityp == bufr_nepn .or. ityp == bufr_ness  &
                 .or. ityp == bufr_neus .or. ityp == bufr_nevs  &
                 .or. obs_bodyElem_i(obsSpaceData,OBS_XTR,index_body) == 0) ) then

               if( ityp == bufr_neus .or. ityp == bufr_nevs ) then
                  varLevel = 'MM'
               else
                  varLevel = 'TH'
               end if
               nlev = col_getNumLev(column,varLevel)
               INDEX_HEADER = obs_bodyElem_i(obsSpaceData,OBS_HIND,INDEX_BODY)
               ITYP = obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY)
               IXTR = obs_bodyElem_i(obsSpaceData,OBS_XTR,INDEX_BODY)
               IK   = obs_bodyElem_i(obsSpaceData,OBS_LYR,INDEX_BODY)
               ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY)
               ZHHH = ZLEV * GRAV
               IPT  = nlev - 1 + col_getOffsetFromVarno(columng,ityp)
               IPB  = IPT+1

               IF (ITYP == BUFR_NETS .OR. ITYP == BUFR_NESS .OR.  &
                    ITYP == BUFR_NEUS .OR. ITYP == BUFR_NEVS) THEN
                  if(ITYP == BUFR_NESS) THEN
                     dPdPsfc = col_getPressureDeriv(columng,nlev,index_header,'TH')
                     columnVarB=hutoes_tl(col_getElem(column,nlev,INDEX_HEADER,'HU'), &
                          col_getElem(column,nlev,INDEX_HEADER,'TT'), &
                          col_getElem(column,1,INDEX_HEADER,'P0'), &
                          col_getElem(columng,nlev,INDEX_HEADER,'HU'), &
                          col_getPressure(columng,nlev,INDEX_HEADER,varLevel),  &
                          dPdPsfc)
                  else
                     columnVarB=col_getElem(COLUMN,IPB,INDEX_HEADER)
                  end if
                  call obs_bodySet_r(obsSpaceData,OBS_WORK,INDEX_BODY,columnVarB)
               ELSEIF (ITYP == BUFR_NEPS .OR. ITYP == BUFR_NEPN) THEN
                  ZLTV  = columng%OLTV(1,nlev,INDEX_HEADER)*col_getElem(COLUMN,nlev,INDEX_HEADER,'TT')  & 
                       + columng%OLTV(2,nlev,INDEX_HEADER)*col_getElem(COLUMN,nlev,INDEX_HEADER,'HU')
                  ZTVG  = columng%OLTV(1,nlev,INDEX_HEADER)*col_getElem(columng,nlev,INDEX_HEADER,'TT')
                  ZGAMAZ= ZGAMMA*(ZHHH-col_getGZsfc(columng,INDEX_HEADER))
                  ZCON  = ((ZTVG-ZGAMAZ)/ZTVG)
                  ZDELPS= (col_getElem(COLUMN,1,INDEX_HEADER,'P0')*ZCON**ZEXP)
                  ZDELTV= ((col_getElem(columng,1,INDEX_HEADER,'P0')*ZEXP*ZCON**(ZEXP-1))  &
                       *(ZGAMAZ/(ZTVG*ZTVG)*ZLTV))
                  call obs_bodySet_r(obsSpaceData,OBS_WORK,INDEX_BODY, ZDELPS+ZDELTV)
               ELSE
                  ! not sure what this block of code is for, not present in nonlinear version (Buehner)
                  IPT  = IK + col_getOffsetFromVarno(columng,ityp)
                  IPB  = IPT+1
                  ZPT  = col_getHeight(columng,IK,INDEX_HEADER,varLevel)
                  ZPB  = col_getHeight(columng,IK+1,INDEX_HEADER,varLevel)
                  ZWB  = (ZPT-ZHHH)/(ZPT-ZPB)
                  ZWT  = 1.d0 - ZWB
                  call obs_bodySet_r(obsSpaceData,OBS_WORK,INDEX_BODY,  &
                       ZWB*col_getElem(COLUMN,IPB,INDEX_HEADER) + ZWT*col_getElem(COLUMN,IPT,INDEX_HEADER)+  &
                       (col_getElem(columng,IPB,INDEX_HEADER)-col_getElem(columng,IPT,INDEX_HEADER)))
               END IF

            end if

         end do BODY

      end do FAMILY

    END subroutine oop_Hsf

    subroutine oop_Hto(obsAssVal)
      !
      ! Purpose: Compute simulated radiances observations from profiled model
      !          increments.
      !          It returns Hdx in OBS_WORK
      !
      !author        : j. halle *cmda/aes  april 8, 2005
      !
      !revision 001  : a. beaulne *cmda/smc  july 2006
      !                    -addition of geopotential field in call to
      !                     tovs_fill_profiles
      !                S. Pellerin, ARMA, August 2008
      !                    - Avoid multiple (iterative) interpolation to 43 levels
      !                      background variable profiles
      !                S. Pellerin, ARMA, January 2009
      !                    - call to oop_storeHdx_radiances instead computing Jo
      !
      implicit none

      integer, intent(in) :: obsAssVal

      integer :: datestamp

      if (.not.obs_famExist(obsSpaceData,'TO',local_mpi=.true.)) return

      !     1.   Prepare atmospheric profiles for all tovs observation points for use in rttov
      !     .    -----------------------------------------------------------------------------
      !
      if (min_nsim == 1) then
         datestamp = tim_getDatestamp()
         if ( trim(obsoperMode) == 'bgckIR') then
           call tvs_fillProfiles(columng,obsSpaceData,datestamp,filt_rlimlvhu,.true.,.false.)
         else
           call tvs_fillProfiles(columng,obsSpaceData,datestamp,filt_rlimlvhu,.false.,.false.)
         end if
      end if


      !     2.   Compute radiance
      !     .    ----------------
      !
      call tvslin_rttov_tl(column, columng, obsSpaceData, obsAssVal)


    end subroutine oop_Hto


    SUBROUTINE oop_Hro(obsAssVal)
      !*
      !* Purpose: Compute the tangent operator for GPSRO observations.
      !*
      !*Author  : J. M. Aparicio Jan 2004
      !*Modified: J. M. Aparicio Dec 2012 adapt to accept bending angle data
      !*    -------------------
      use IndexListDepot_mod, only : struct_index_list
      implicit none

      integer, intent(in) :: obsAssVal

      REAL*8 zLat, Lat, sLat
      REAL*8 zLon, Lon
      REAL*8 zAzm, Azm
      INTEGER IAZM, ISAT
      REAL*8 Rad, Geo, WFGPS, zP0
      REAL*8, allocatable :: zPP(:)
      REAL*8, allocatable :: zDP(:)
      REAL*8, allocatable :: zTT(:)
      REAL*8, allocatable :: zHU(:)
      REAL*8, allocatable :: zUU(:)
      REAL*8, allocatable :: zVV(:)
      REAL*8 zMT,radw

      REAL*8 ZMHXL
      REAL*8 DX (ngpscvmx)

      INTEGER IDATYP
      INTEGER JL, JV, NGPSLEV, NWNDLEV, stat1, JJ
      integer :: index_header, index_body, iProfile
      type(struct_index_list), pointer :: local_current_list

      LOGICAL  ASSIM, LFIRST

      INTEGER NH, NH1
      TYPE(GPS_PROFILE)           :: PRF
      REAL*8       , ALLOCATABLE :: H   (:),AZMV(:)
      TYPE(GPS_DIFF), ALLOCATABLE :: RSTV(:),RSTVP(:),RSTVM(:)
      !C      WRITE(*,*)'ENTER oop_Hro'
      !C
      !C     * 1.  Initializations
      !C     *     ---------------
      !C
      NGPSLEV=col_getNumLev(column,'TH')
      NWNDLEV=col_getNumLev(column,'MM')
      LFIRST=.FALSE.
      if ( .NOT.allocated(gps_vRO_Jacobian) ) then
         LFIRST = .TRUE.
         allocate(zPP (NGPSLEV))
         allocate(zDP (NGPSLEV))
         allocate(zTT (NGPSLEV))
         allocate(zHU (NGPSLEV))
         allocate(zUU (NGPSLEV))
         allocate(zVV (NGPSLEV))

         allocate(gps_vRO_Jacobian(gps_numROProfiles,GPSRO_MAXPRFSIZE,2*NGPSLEV+1))
         allocate(gps_vRO_lJac    (gps_numROProfiles))
         gps_vRO_lJac=.false.

         allocate( H    (GPSRO_MAXPRFSIZE) )
         allocate( AZMV (GPSRO_MAXPRFSIZE) )
         allocate( RSTV (GPSRO_MAXPRFSIZE) )
         !C         IF (LEVELGPSRO == 1) THEN
         !C            allocate( RSTVP(GPSRO_MAXPRFSIZE) )
         !C            allocate( RSTVM(GPSRO_MAXPRFSIZE) )
         !C         END IF
      end if
      !C
      !C    Loop over all header indices of the 'RO' family (Radio Occultation)
      !C
      ! Set the header list (start at the beginning of the list)
      call obs_set_current_header_list(obsSpaceData,'RO')
      !##$omp parallel default(shared) &
      !##$omp private(index_header,idatyp,assim,nh,local_current_list,index_body) &
      !##$omp private(iProfile,irad,igeo,iazm,isat,rad,geo,zazm,zmt,wfgps,jj) &
      !##$omp private(zlat,zlon,lat,lon,azm,slat) &
      !##$omp private(stat1,jl,zpp,zdp,ztt,zhu,zuu,zvv,prf,dx) &
      !##$omp private(h,azmv,rstv,rstvp,rstvm,nh1,zmhxl,jv)
      nullify(local_current_list)
      HEADER: do
         INDEX_HEADER = obs_getHeaderIndex(obsSpaceData)
         if (INDEX_HEADER < 0) exit HEADER
         !C
         !C     * Process only refractivity data (codtyp 169)
         !C
         IDATYP = obs_headElem_i(obsSpaceData,OBS_ITY,INDEX_HEADER)
         DATYP: IF ( IDATYP == 169 ) THEN
            !C
            !C     *    Scan for requested data values of the profile, and count them
            !C
            ASSIM = .FALSE.
            NH = 0
            !C
            !C     *    Loop over all body indices for this index_header:
            !C     *    (start at the beginning of the list)
            !C
            call obs_set_current_body_list(obsSpaceData, INDEX_HEADER, &
                 current_list=local_current_list)
            BODY: do 
               index_body = obs_getBodyIndex(local_current_list)
               if (index_body < 0) exit BODY
               IF ( obs_bodyElem_i(obsSpaceData,OBS_ASS,INDEX_BODY) == obsAssVal ) THEN
                  ASSIM = .TRUE.
                  NH = NH + 1
               END IF
            END DO BODY
            !C
            !C     *    If assimilations are requested, prepare and apply the observation operator
            !C
            ASSIMILATE: IF (ASSIM) THEN
               iProfile=gps_iprofile_from_index(INDEX_HEADER)
               !C
               !C     *       Profile at the observation location:
               !C
               if (.not.gps_vRO_lJac(iProfile)) then
                  !C
                  !C     *          Basic geometric variables of the profile:
                  !C
                  zLat = obs_headElem_r(obsSpaceData,OBS_LAT,INDEX_HEADER)
                  zLon = obs_headElem_r(obsSpaceData,OBS_LON,INDEX_HEADER)
                  IAZM = obs_headElem_i(obsSpaceData,OBS_AZA,INDEX_HEADER)
                  ISAT = obs_headElem_i(obsSpaceData,OBS_SAT,INDEX_HEADER)
                  Rad  = obs_headElem_r(obsSpaceData,OBS_TRAD,INDEX_HEADER)
                  Geo  = obs_headElem_r(obsSpaceData,OBS_GEOI,INDEX_HEADER)
                  zAzm = 0.01d0*IAZM / MPC_DEGREES_PER_RADIAN_R8
                  zMT  = col_getGZsfc(columng,INDEX_HEADER)/RG
                  WFGPS= 0.d0
                  DO JJ=1,NUMGPSSATS
                     IF (ISAT == IGPSSAT(JJ)) WFGPS=WGPS(JJ)
                  END DO
                  Lat  = zLat * MPC_DEGREES_PER_RADIAN_R8
                  Lon  = zLon * MPC_DEGREES_PER_RADIAN_R8
                  Azm  = zAzm * MPC_DEGREES_PER_RADIAN_R8
                  sLat = sin(zLat)
                  zMT  = zMT * RG / gpsgravitysrf(sLat)
                  zP0  = col_getElem(columng,1,index_header,'P0')
                  DO JL = 1, NGPSLEV
                     !C
                     !C     *             Profile x_b
                     !C
                     zPP(JL) = col_getPressure(columng,JL,INDEX_HEADER,'TH')
                     !C     *             True implementation of zDP (dP/dP0)
                     zDP(JL) = col_getPressureDeriv(columng,JL,INDEX_HEADER,'TH')
                     zTT(JL) = col_getElem(columng,JL,INDEX_HEADER,'TT') - p_TC
                     zHU(JL) = col_getElem(columng,JL,INDEX_HEADER,'HU')
                  END DO

                  if((col_getPressure(columng,1,index_header,'TH') + 1.0d-4) <  &
                       col_getPressure(columng,1,index_header,'MM')) then
                     ! case with top thermo level above top momentum level (Vcode=5002)
                     do jl = 1, nwndlev
                        zuu(jl) = col_getElem(columng,jl  ,index_header,'UU') * p_knot
                        zvv(jl) = col_getElem(columng,jl  ,index_header,'VV') * p_knot
                     end do
                  else
                     ! case without top thermo above top momentum level or unstaggered (Vcode=5001/4/5)
                     do jl = 1, nwndlev-1
                        zuu(jl) = col_getElem(columng,jl+1,index_header,'UU') * p_knot
                        zvv(jl) = col_getElem(columng,jl+1,index_header,'VV') * p_knot
                     end do
                     zuu(nwndlev) = zuu(nwndlev-1)
                     zvv(nwndlev) = zuu(nwndlev-1)
                  end if
                  zuu(ngpslev) = zuu(nwndlev)
                  zvv(ngpslev) = zuu(nwndlev)
                  !C     
                  !C     *          GPS profile structure:
                  !C
                  call gps_struct1sw(ngpslev,zLat,zLon,zAzm,zMT,Rad,geo,zP0,zPP,zDP,zTT,zHU,zUU,zVV,prf)
                  !C
                  !C     *          Prepare the vector of all the observations:
                  !C
                  NH1 = 0
                  call obs_set_current_body_list(obsSpaceData, index_header, &
                       current_list=local_current_list)
                  BODY_2: do 
                     INDEX_BODY = obs_getBodyIndex(local_current_list)
                     if (INDEX_BODY < 0) exit BODY_2
                     IF ( obs_bodyElem_i(obsSpaceData,OBS_ASS,INDEX_BODY) == obsAssVal ) THEN
                        NH1      = NH1 + 1
                        H(NH1)   = obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY)
                        AZMV(NH1)= zAzm
                     END IF
                  END DO BODY_2
                  !C
                  !C     *          Apply the observation operator:
                  !C
                  IF (LEVELGPSRO == 1) THEN
                     CALL GPS_BNDOPV1(H      , AZMV, NH, PRF, RSTV)
                     !C                     CALL GPS_BNDOPV1(H+WFGPS, AZMV, NH, PRF, RSTVP)
                     !C                     CALL GPS_BNDOPV1(H-WFGPS, AZMV, NH, PRF, RSTVM)
                     !C                     do nh1 = 1, nh
                     !C                        RSTV(nh1)=(RSTVP(nh1)+RSTV(nh1)+RSTVM(nh1))/3.d0
                     !C                     end do
                  ELSE
                     CALL GPS_REFOPV (H,       NH, PRF, RSTV)
                  END IF
                  DO NH1=1,NH
                     gps_vRO_Jacobian(iProfile,NH1,:)= RSTV(NH1)%DVAR(1:2*NGPSLEV+1)
                  END DO
                  gps_vRO_lJac(iProfile)=.true.
               end if
               !C
               !C     *       Local vector state
               !C
               DO JL = 1, NGPSLEV
                  DX (        JL) = col_getElem(COLUMN,JL,index_header,'TT')
                  DX (NGPSLEV+JL) = col_getElem(COLUMN,JL,index_header,'HU')
               END DO
               DX (2*NGPSLEV+1)   = col_getElem(COLUMN,1 ,index_header,'P0')
               !C
               !C     *       Perform the (H(xb)DX-Y') operation
               !C     *       Loop over all body indices for this index_header:
               !C
               NH1 = 0
               call obs_set_current_body_list(obsSpaceData, index_header, &
                    current_list=local_current_list)
               BODY_3: do 
                  INDEX_BODY = obs_getBodyIndex(local_current_list)
                  if (INDEX_BODY < 0) exit BODY_3
                  IF ( obs_bodyElem_i(obsSpaceData,OBS_ASS,INDEX_BODY) == obsAssVal ) THEN
                     NH1 = NH1 + 1
                     !C
                     !C     *             Evaluate H(xb)DX
                     !C
                     ZMHXL = 0.d0
                     DO JV = 1, 2*NGPSLEV+1
                        ZMHXL = ZMHXL + gps_vRO_Jacobian(iProfile,NH1,JV) * DX(JV)
                     END DO
                     !C
                     !C     *             Store in CMA
                     !C
                     call obs_bodySet_r(obsSpaceData,OBS_WORK,INDEX_BODY, ZMHXL)
                  END IF
               END DO BODY_3
            END IF ASSIMILATE
         END IF DATYP
      END DO HEADER
      !##$omp end parallel

      IF (LFIRST) THEN
         !C         IF (LEVELGPSRO == 1) THEN
         !C            deallocate( RSTVM )
         !C            deallocate( RSTVP )
         !C         END IF
         deallocate( RSTV )
         deallocate( AZMV )
         deallocate( H    )

         deallocate(zVV)
         deallocate(zUU)
         deallocate(zHU)
         deallocate(zTT)
         deallocate(zDP)
         deallocate(zPP)
      END IF

      !C      WRITE(*,*)'EXIT oop_Hro'
      RETURN
    END subroutine oop_Hro


    SUBROUTINE oop_Hzp(obsAssVal)
      !*
      !* Purpose: Compute simulated profiler observations from profiled model
      !*          increments.
      !*          It returns Hdx in OBS_WORK
      !*          Interpolate vertically the contents of commvo to heights
      !*          (in meters) of the observations.
      !*          A linear interpolation in z is performed.
      !*
      !*Author  :  J. St-James, CMDA/SMC July 2003

      implicit none

      integer, intent(in) :: obsAssVal

      INTEGER IPB,IPT
      INTEGER INDEX_HEADER,IK
      INTEGER J,INDEX_BODY,ITYP
      REAL*8 ZVAR,ZDA1,ZDA2
      REAL*8 ZWB,ZWT
      REAL*8 ZLEV,ZPT,ZPB,ZDENO
      character(len=2) :: varLevel

      call obs_set_current_body_list(obsSpaceData, 'PR')
      BODY: do
         index_body = obs_getBodyIndex(obsSpaceData)
         if (index_body < 0) exit BODY

         IF (   (obs_bodyElem_i(obsSpaceData,OBS_ASS,INDEX_BODY) == obsAssVal) &
              .AND. (obs_bodyElem_i(obsSpaceData,OBS_XTR,INDEX_BODY) == 0) &
              .AND. (obs_bodyElem_i(obsSpaceData,OBS_VCO,INDEX_BODY) == 1)  ) THEN
            INDEX_HEADER = obs_bodyElem_i(obsSpaceData,OBS_HIND,INDEX_BODY)
            ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY)
            IK   = obs_bodyElem_i(obsSpaceData,OBS_LYR,INDEX_BODY)
            ITYP = obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY)
            varLevel = vnl_varLevelFromVarnum(ityp)
            IPT  = IK + col_getOffsetFromVarno(columng,ityp)
            IPB  = IPT+1
            ZPT  = col_getHeight(columng,IK  ,INDEX_HEADER,varLevel)/RG
            ZPB  = col_getHeight(columng,IK+1,INDEX_HEADER,varLevel)/RG
            ZDENO= ZPT-ZPB
            ZWB  = (ZPT-ZLEV)/ZDENO
            ZWT  = 1.0D0 - ZWB

            ZDA1= (ZLEV-ZPB)/(ZDENO**2)
            ZDA2= (ZPT-ZLEV)/(ZDENO**2)

            if(ITYP == BUFR_NEES) then
               write(*,*) 'ABORTING IN OOP_HZP: CANNOT ASSIMILATE ES!!!',ityp,obs_getfamily(obsSpaceData,index_header),index_header,index_body
               call utl_abort('Aborting in oop_H')
            end if
            call obs_bodySet_r(obsSpaceData,OBS_WORK,INDEX_BODY,  &
                 ZWB*col_getElem(COLUMN,IPB,INDEX_HEADER) + ZWT*col_getElem(COLUMN,IPT,INDEX_HEADER) +  &
                 (col_getElem(columng,IPB,INDEX_HEADER) - col_getElem(columng,IPT,INDEX_HEADER))*  &
                 (ZDA1*col_getHeight(COLUMN,IK,INDEX_HEADER,varLevel)/RG + ZDA2*col_getHeight(COLUMN,IK+1,INDEX_HEADER,varLevel)/RG))
         END IF
      END DO BODY
      RETURN
    END subroutine oop_Hzp


    SUBROUTINE oop_Hgp(obsAssVal)
      !*
      !***s/r  -oop_Hgp TL of DOBSGPSGB (Jo for GB-GPS ZTD observations)
      !*
      !*
      !*Author  : S. Macpherson *ARMA October 2012
      !*    -------------------
      !**    Purpose: Compute H'dx for all GPS ZTD observations
      !*
      implicit none

      integer, intent(in) :: obsAssVal

      REAL*8 ZLAT, Lat
      REAL*8 ZLON, Lon
      REAL*8, allocatable :: ZTTB(:)
      REAL*8, allocatable :: ZHUB(:)
      REAL*8, allocatable :: ZPPB(:)
      REAL*8, allocatable :: ZDP(:)
      REAL*8 ZP0B, ZPSMOD, ZPWMOD, ZPWMOD2, dZTD
      REAL*8 ZMT
      real*8 sfcfield

      REAL*8 ZHX, ZLEV, ZDZMIN
      REAL*8 JAC(ngpscvmx)
      REAL*8 DX (ngpscvmx)

      INTEGER INDEX_HEADER, INDEX_BODY
      INTEGER JL, NFLEV, status, iztd, icount, stat, vcode

      LOGICAL      ASSIM

      real*8, dimension(:), pointer :: dpdp0 => null()

      TYPE(gps_profilezd)   :: PRF, PRF2
      TYPE(gps_diff)        :: ZTDOPV, ZTDOPV2

      !      WRITE(*,*)'ENTER oop_Hgp'

      stat = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=vcode)

      ZDZMIN = DZMIN                     ! from modgpsztd_mod

      NFLEV  = col_getNumLev(columng,'TH')

      !C
      !C     * 1.  Initializations
      !C     *     ---------------
      !C
      !     NOTE:  vGPSZTD_Index(numGPSZTD) is initialized in s/r dobsgpsgb
      !
      if (.not.allocated(vGPSZTD_Index)) then
         call utl_abort('oop_Hgp: ERROR:  vGPSZTD_Index not allocated!')
      elseif (.not.allocated(vGPSZTD_Jacobian)) then
         write(*,*) ' Allocate vGPSZTD_Jacobian(numGPSZTD,2*NFLEV+1)'
         allocate(vGPSZTD_Jacobian(numGPSZTD,2*NFLEV+1))
         allocate(vGPSZTD_lJac(numGPSZTD))
         vGPSZTD_lJac = .false.
         vGPSZTD_Jacobian = 0.0d0
      end if

      !   If first time (iteration), store the Jacobians for all ZTD data to be assimilated

      !-----------------------------------------------------------------------------------------
      INIT: IF ( .not.vGPSZTD_lJac(1) ) THEN

         allocate(ZTTB(NFLEV))
         allocate(ZHUB(NFLEV))
         allocate(ZPPB(NFLEV))
         allocate(ZDP(NFLEV))

         write(*,*) 'oop_Hgp: Storing Jacobians for GPS ZTD data ...'
         write(*,*) '   INFO: Analysis grid iversion = ', vcode
         write(*,*) '         col_getNumLev(columng,TH) = ', NFLEV
         write(*,*) '         numGPSZTD = ', numGPSZTD

         icount = 0

         ! loop over all header indices of the 'GP' family (GPS observations)
         call obs_set_current_header_list(obsSpaceData,'GP')
         HEADER_0: do
            index_header = obs_getHeaderIndex(obsSpaceData)
            if (index_header < 0) exit HEADER_0
            !C
            !C     *     Scan for ZTD assimilation at this location
            !C
            ASSIM = .FALSE.
            ! loop over all body indices for this index_header
            call obs_set_current_body_list(obsSpaceData, index_header)
            BODY_0: DO 
               index_body = obs_getBodyIndex(obsSpaceData)
               if (index_body < 0) exit BODY_0
               if (   (obs_headElem_i(obsSpaceData,OBS_ITY,INDEX_HEADER) == 189) &
                    .and. (obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY) == BUFR_NEZD) &
                    .and. (obs_bodyElem_i(obsSpaceData,OBS_ASS,INDEX_BODY) == obsAssVal) ) then
                  ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY)
                  ASSIM = .TRUE.
               end if
            END DO BODY_0

            IF ( ASSIM ) THEN
               !C
               !C     *    LR background profile at the observation location x :
               !C
               icount = icount + 1
               Lat  = obs_headElem_r(obsSpaceData,OBS_LAT,INDEX_HEADER)
               ZLAT = Lat * MPC_DEGREES_PER_RADIAN_R8
               Lon  = obs_headElem_r(obsSpaceData,OBS_LON,INDEX_HEADER)
               ZLON = Lon * MPC_DEGREES_PER_RADIAN_R8
               ZP0B = col_getElem(columng,1,INDEX_HEADER,'P0')
               DO JL = 1, NFLEV
                  ZTTB(JL) = col_getElem(columng,JL,INDEX_HEADER,'TT') - p_TC
                  ZHUB(JL) = col_getElem(columng,JL,INDEX_HEADER,'HU')
                  ZPPB(JL) = col_getPressure(columng,JL,INDEX_HEADER,'TH')
                  !C            Get ZDP = dP/dP0
                  ZDP(JL)  = col_getPressureDeriv(columng,JL,INDEX_HEADER,'TH')
               END DO
               if ( ZPPB(NFLEV) /= ZP0B ) then
                  write(*,*) ' oop_Hgp: ERROR: ZPPB(NFLEV) /= ZP0B'
                  write(*,*) '          ZPPB(NFLEV), ZP0B =', ZPPB(NFLEV), ZP0B
                  !BUE              call utl_abort('oop_Hgp')
               end if
               ZMT  = col_getGZsfc(columng,INDEX_HEADER)/RG
               if ( icount == 1 .and. LTESTOP ) write(*,*) 'ZDP (dpdp0) = ', (ZDP(JL),JL= 1,NFLEV)
               !c
               CALL gps_structztd(NFLEV,Lat,Lon,ZMT,ZP0B,ZPPB,ZDP,ZTTB,ZHUB,LBEVIS,IREFOPT,PRF)
               CALL gps_ztdopv(ZLEV,PRF,LBEVIS,ZDZMIN,ZTDopv,ZPSMOD,IZTDOP)
               !C          Observation Jacobian H'(xb)            
               JAC = ZTDopv%DVar
               iztd = gps_i_from_index(INDEX_HEADER)
               DO JL = 1, 2*NFLEV+1
                  vGPSZTD_Jacobian(iztd,JL) = JAC(JL)
               END DO
               vGPSZTD_lJac(iztd) = .true.
               !            
               if ( icount <= 3 .and. LTESTOP ) then
                  write(*,*) '--------------------------------------------------------- '
                  write(*,*) iztd, obs_elem_c(obsSpaceData,'STID',INDEX_HEADER),'ZTDopv (m) = ', ZTDopv%Var
                  CALL gps_pw(PRF,ZPWMOD)
                  !           sfc pressure dx               
                  !               ZPPB(NFLEV) = ZPPB(NFLEV) + 50.0d0
                  nullify(dpdp0)
                  sfcfield = ZP0B + 50.0d0
                  status = vgd_dpidpis(vco_anl%vgrid,vco_anl%ip1_T,dpdp0,sfcfield)
                  ZDP = dpdp0(1:NFLEV)
                  CALL gps_structztd(NFLEV,Lat,Lon,ZMT,sfcfield,ZPPB,ZDP,ZTTB,ZHUB,LBEVIS,IREFOPT,PRF2)
                  CALL gps_ztdopv(ZLEV,PRF2,LBEVIS,ZDZMIN,ZTDopv2,ZPSMOD,IZTDOP)
                  write(*,*) ' ZTD Operator Test:  dP0 = +50 Pa'
                  write(*,*) ' dZTD NL     = ', ZTDopv2%Var - ZTDopv%Var
                  write(*,*) ' dZTD Linear = ', vGPSZTD_Jacobian(iztd,2*NFLEV+1)*50.0d0
                  write(*,*) ' '
                  !               ZPPB(NFLEV) = ZPPB(NFLEV) - 50.0d0
                  !           log(q) dx 
                  ZHUB(64) = ZHUB(64) - 0.44D-01
                  ZHUB(65) = ZHUB(65) - 0.44D-01
                  ZHUB(66) = ZHUB(66) - 0.44D-01
                  CALL gps_structztd(NFLEV,Lat,Lon,ZMT,ZP0B,ZPPB,ZDP,ZTTB,ZHUB,LBEVIS,IREFOPT,PRF2)
                  CALL gps_ztdopv(ZLEV,PRF2,LBEVIS,ZDZMIN,ZTDopv2,ZPSMOD,IZTDOP)
                  CALL gps_pw(PRF2,ZPWMOD2)
                  write(*,*) ' ZTD Operator Test:  dLQ = -0.44E-01 JL = 64,65,66'
                  write(*,*) ' dPW (mm)    = ', ZPWMOD2 - ZPWMOD
                  write(*,*) ' dZTD NL     = ', ZTDopv2%Var - ZTDopv%Var
                  dZTD = vGPSZTD_Jacobian(iztd,64+NFLEV)*(-0.44D-01) + vGPSZTD_Jacobian(iztd,65+NFLEV)*(-0.44D-01) + &
                       vGPSZTD_Jacobian(iztd,66+NFLEV)*(-0.44D-01)
                  write(*,*) ' dZTD Linear = ', dZTD
                  write(*,*) ' '
                  ZHUB(64) = ZHUB(64) + 0.44D-01
                  ZHUB(65) = ZHUB(65) + 0.44D-01
                  ZHUB(66) = ZHUB(66) + 0.44D-01
                  !           temperature dx                   
                  ZTTB(64) = ZTTB(64) + 2.0d0
                  ZTTB(65) = ZTTB(65) + 2.0d0
                  ZTTB(66) = ZTTB(66) + 2.0d0
                  CALL gps_structztd(NFLEV,Lat,Lon,ZMT,ZP0B,ZPPB,ZDP,ZTTB,ZHUB,LBEVIS,IREFOPT,PRF2)
                  CALL gps_ztdopv(ZLEV,PRF2,LBEVIS,ZDZMIN,ZTDopv2,ZPSMOD,IZTDOP)
                  write(*,*) ' ZTD Operator Test:  dTT = +2.0K JL = 64,65,66'
                  write(*,*) ' dZTD NL     = ', ZTDopv2%Var - ZTDopv%Var
                  dZTD = vGPSZTD_Jacobian(iztd,64)*2.0d0 + vGPSZTD_Jacobian(iztd,65)*2.0d0 + &
                       vGPSZTD_Jacobian(iztd,66)*2.0d0
                  write(*,*) ' dZTD Linear = ', dZTD
                  write(*,*) '--------------------------------------------------------- '              
               end if

            END IF

         END DO HEADER_0

         deallocate(ZTTB)
         deallocate(ZHUB)
         deallocate(ZPPB)
         deallocate(ZDP)

         write(*,*) 'oop_Hgp:   Number of ZTD data (icount) = ', icount
         write(*,*) '           Expected number (numGPSZTD) = ', numGPSZTD
         write(*,*) '           Last iztd                   = ', iztd
         write(*,*) '           vGPSZTD_Index(1)            = ', vGPSZTD_Index(1)
         write(*,*) '           vGPSZTD_Index(iztd)         = ', vGPSZTD_Index(iztd)

         if ( icount /= numGPSZTD ) then
            call utl_abort('oop_Hgp: ERROR: icount /= numGPSZTD!')
         end if
         if ( icount /= iztd ) then
            call utl_abort('oop_Hgp: ERROR: icount /= iztd!')
         end if
         if ( numGPSZTD /= iztd ) then
            call utl_abort('oop_Hgp: ERROR: numGPSZTD /= iztd!')
         end if


      END IF INIT
      !-----------------------------------------------------------------------------------------

      icount = 0

      ! loop over all header indices of the 'GP' family (GPS observations)
      call obs_set_current_header_list(obsSpaceData,'GP')
      HEADER: do
         index_header = obs_getHeaderIndex(obsSpaceData)
         if (index_header < 0) exit HEADER
         !C
         !C     *     Scan for ZTD assimilation at this location
         !C
         ASSIM = .FALSE.
         ! loop over all body indices for this index_header
         call obs_set_current_body_list(obsSpaceData, index_header)
         BODY: DO 
            index_body = obs_getBodyIndex(obsSpaceData)
            if (index_body < 0) exit BODY

            if (   (obs_headElem_i(obsSpaceData,OBS_ITY,INDEX_HEADER) == 189) &
                 .and. (obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY) == BUFR_NEZD) &
                 .and. (obs_bodyElem_i(obsSpaceData,OBS_ASS,INDEX_BODY) == obsAssVal) ) then
               ASSIM = .TRUE.
            END IF
         END DO BODY
         !C
         !C     * If ZTD assimilation, apply the TL observation operator
         !C
         IF ( ASSIM ) THEN
            iztd = gps_i_from_index(INDEX_HEADER)
            if ( iztd < 1 .or. iztd > numGPSZTD ) then
               call utl_abort('oop_Hgp: ERROR: index from gps_i_from_index() is out of range!')
            end if
            !C
            !C     *    Local vector state (analysis increments)
            !C
            DO JL = 1, NFLEV
               DX (JL)        = col_getElem(COLUMN,JL,index_header,'TT')
               DX (NFLEV+JL)  = col_getElem(COLUMN,JL,index_header,'HU')
            END DO
            DX (2*NFLEV+1)    = col_getElem(COLUMN,1 ,index_header,'P0')

            !C     *    Evaluate H'(xb)*dX
            !C
            ZHX = 0.D0
            DO JL = 1, 2*NFLEV+1
               ZHX = ZHX + vGPSZTD_Jacobian(iztd,JL)*DX(JL)
            END DO
            !C
            !C     *    Store ZHX = H'dx in OBS_WORK
            !C
            ! loop over all body indices for this index_header
            call obs_set_current_body_list(obsSpaceData, index_header)
            BODY_2: DO 
               index_body = obs_getBodyIndex(obsSpaceData)
               if (index_body < 0) exit BODY_2

               IF (   (obs_headElem_i(obsSpaceData,OBS_ITY,INDEX_HEADER) == 189) &
                    .and. (obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY) == BUFR_NEZD) &
                    .and. (obs_bodyElem_i(obsSpaceData,OBS_ASS,INDEX_BODY) == obsAssVal) ) then
                  call obs_bodySet_r(obsSpaceData,OBS_WORK,INDEX_BODY, ZHX)
                  icount = icount + 1
                  if ( icount <= 3 .and. LTESTOP ) then
                     write(*,*) iztd, obs_elem_c(obsSpaceData,'STID',INDEX_HEADER)
                     write(*,*) 'JAC(ncv) = ', (vGPSZTD_Jacobian(iztd,JL),JL=1,2*NFLEV+1)
                     write(*,*) 'DTT(JL)  = ', (DX(JL),JL=1,NFLEV)
                     write(*,*) 'DHU(JL)  = ', (DX(JL),JL=NFLEV+1,2*NFLEV)
                     write(*,*) 'DP0(JL)  = ', DX(2*NFLEV+1)
                     write(*,*) 'ZHX (mm) = ', ZHX*1000.D0
                  end if
               END IF
            END DO BODY_2

         END IF ! ASSIM

      END DO HEADER

      !      WRITE(*,*) 'oop_Hgp: Number of ZTD data locations with obs_bodySet_r(OBS_OMA) = ', icount

      !      WRITE(*,*)'EXIT oop_Hgp'

      RETURN
    END subroutine oop_Hgp


    subroutine oop_Hchm(obsAssVal)
      !
      !**s/r oop_Hchm- Compute simulated chemical constituents observations from profiled model
      !                increments, and returns Hdx in OBS_WORK
      !
      !*Author: M. Sitwell, ARQI/AQRD, Aug 2015
      !         - Reduced to a call of chm_observation_operators
      !           which contains much of the original version by Ping Du (CMDA/MSC) and
      !           Y. Rochon (ARQI/AQRD), Jan 2015 (partially based on corresponding
      !           pre-EnVar routine by Y.J. Rochon and Y. Yang, July 2005 to Feb 2013),
      !           with further original changes by Y. Rochon and M. Sitwell, 2015. 
      !
      !**   Purpose: Compute simulated chemical constituents observations from profiled model    
      !              increments, and returns Hdx in OBS_WORK
      !
      !-----------------------------------------------------------------------------------

      implicit none

      integer, intent(in) :: obsAssVal

      if (.not.obs_famExist(obsSpaceData,'CH',local_mpi=.true.)) return
      
      call chm_observation_operators(columng,obsSpaceData,kmode=2,columnInc_opt=column,obsAssVal_opt=obsAssVal) ! kmode=2 for tangent linear operator

    end subroutine oop_Hchm


  end subroutine oop_Htl


  subroutine oop_Had(column,columng,obsSpaceData)
    implicit none
    !
    !Call the several adjoint of observation operators
    !
    type(struct_columnData) :: column,columng
    type(struct_obs) :: obsSpaceData
    type(struct_vco), pointer :: vco_anl
    logical, save :: firstTime = .true.

    IF(mpi_myid == 0) THEN
       write(*,*)'OOP_HT- Adjoint of linearized observation operators'
    end if

    vco_anl => col_getVco(columng)

    !     Find interpolation layer in model profiles (used by several operators)
    if ( firstTime ) then
      if ( col_getNumLev(columng,'MM') > 1 ) call oop_vobslyrs(columng,obsSpaceData)
      firstTime = .false.
    end if

    call tmg_start(125,'OBS_CHM_TLAD') !
    call oop_HTchm
    call tmg_stop (125)

    call tmg_start(47,'OBS_GPSGB_TLAD') !
    if (numGPSZTD > 0) call oop_HTgp
    call tmg_stop (47)        !

    call tmg_start(46,'OBS_ZZZ_TLAD') !
    call oop_HTzp
    call tmg_stop (46)

    call tmg_start(45,'OBS_GPSRO_TLAD') !
    call oop_HTro
    call tmg_stop (45)      !     !

    call tmg_start(44,'OBS_TOV_TLAD') !
    call oop_HTto
    call tmg_stop (44)      !

    call tmg_start(43,'OBS_SFC_TLAD')
    call oop_HTsf
    call tmg_stop (43)      !

    call tmg_start(42,'OBS_PPP_TLAD') !
    call oop_HTpp
    call tmg_stop (42)

  CONTAINS

    SUBROUTINE oop_HTpp
      !
      !**s/r   - Adjoint of the "vertical" interpolation
      !          for "UPPER AIR" data files.
      !
      !
      !
      !Author  : P. Koclas *CMC/AES  April 1996
      !
      !     Purpose: based on vint3d to build the adjoint of the
      !              vertical interpolation for UPPER-AIR data files.
      !
      implicit none
      INTEGER IPB,IPT,ITYP
      REAL*8 ZRES
      REAL*8 ZWB,ZWT
      REAL*8 ZLEV,ZPT,ZPB,ZDADPS,ZPRESBPB,ZPRESBPT
      INTEGER INDEX_HEADER,IK,nlev_T
      INTEGER INDEX_BODY,INDEX_FAMILY
      REAL*8 columngVarT,columngVarB
      real*8, pointer :: all_column(:),tt_column(:),hu_column(:),ps_column(:)
      REAL*8 :: dPdPsT,dPdPsB
      logical :: llassim
      INTEGER, PARAMETER :: numFamily=3
      CHARACTER(len=2) :: list_family(numFamily),varLevel
      !
      list_family(1) = 'UA'
      list_family(2) = 'AI'
      list_family(3) = 'SW'

      FAMILY: do index_family=1,numFamily

         call obs_set_current_body_list(obsSpaceData,list_family(index_family))
         BODY: DO 
            index_body = obs_getBodyIndex(obsSpaceData)
            if (index_body < 0) exit BODY

            llassim= (obs_bodyElem_i(obsSpaceData,OBS_ASS,index_body) == 1) &
                 .AND. (obs_bodyElem_i(obsSpaceData,OBS_XTR,index_body) == 0) &
                 .AND. (obs_bodyElem_i(obsSpaceData,OBS_VCO,index_body) == 2)

            ! Process all data within the domain of the model
            IF (llassim) THEN
               index_header = obs_bodyElem_i(obsSpaceData,OBS_HIND,INDEX_BODY)
               ZRES = obs_bodyElem_r(obsSpaceData,OBS_WORK,INDEX_BODY)
               ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY)
               ITYP = obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY)
               varLevel = vnl_varLevelFromVarnum(ityp)
               IK   = obs_bodyElem_i(obsSpaceData,OBS_LYR,INDEX_BODY)
               IPT  = IK  + col_getOffsetFromVarno(columng,ityp)
               IPB  = IPT+1
               ZPT  = col_getPressure(COLUMNG,IK,INDEX_HEADER,varLevel)
               ZPB  = col_getPressure(COLUMNG,IK+1,INDEX_HEADER,varLevel)
               dPdPsT = col_getPressureDeriv(COLUMNG,IK  ,INDEX_HEADER,varLevel)
               dPdPsB = col_getPressureDeriv(COLUMNG,IK+1,INDEX_HEADER,varLevel)
               ZWB  = LOG(ZLEV/ZPT)/LOG(ZPB/ZPT)
               ZWT  = 1.0D0 - ZWB

               ZDADPS   = ( LOG(ZLEV/ZPB)*dPdPsT/ZPT -   &
                    LOG(ZLEV/ZPT)*dPdPsB/ZPB )/  &
                    LOG(ZPB/ZPT)**2

               all_column => col_getColumn(column,INDEX_HEADER)
               tt_column  => col_getColumn(column,INDEX_HEADER,'TT')
               hu_column  => col_getColumn(column,INDEX_HEADER,'HU')
               ps_column  => col_getColumn(column,INDEX_HEADER,'P0')

               if(ITYP.eq.BUFR_NEES) then
                  call hutoes_ad(hu_column(IK+1),  &
                       tt_column(IK+1),  &
                       ps_column(1),     &
                       ZWB*ZRES,         &
                       col_getElem(columng,IK+1,INDEX_HEADER,'HU'),      &
                       col_getPressure(columng,IK+1,INDEX_HEADER,'TH'),  &
                       dPdPsB)
                  call hutoes_ad(hu_column(IK  ),  &
                       tt_column(IK  ),  &
                       ps_column(1),     &
                       ZWT*ZRES,         &
                       col_getElem(columng,IK  ,INDEX_HEADER,'HU'),      &
                       col_getPressure(columng,IK  ,INDEX_HEADER,'TH'),  &
                       dPdPsT)
                  columngVarB=hutoes(col_getElem(columng,IK+1,INDEX_HEADER,'HU'),  &
                       col_getElem(columng,IK+1,INDEX_HEADER,'TT'),  &
                       col_getPressure(columng,IK+1,INDEX_HEADER,'TH'))
                  columngVarT=hutoes(col_getElem(columng,IK  ,INDEX_HEADER,'HU'),  &
                       col_getElem(columng,IK  ,INDEX_HEADER,'TT'),  &
                       col_getPressure(columng,IK  ,INDEX_HEADER,'TH'))
               else
                  all_column(IPB) = all_column(IPB) + ZWB*ZRES
                  all_column(IPT) = all_column(IPT) + ZWT*ZRES
                  columngVarB=col_getElem(columng,IPB,INDEX_HEADER)
                  columngVarT=col_getElem(columng,IPT,INDEX_HEADER)
               end if
               ps_column(1)    = ps_column(1)    +         &
                    (columngVarB - columngVarT)  &
                    *ZDADPS*ZRES

            end if

         end do BODY

      end do FAMILY

    end subroutine oop_HTpp


    SUBROUTINE oop_HTsf
      !*
      !***s/r AOBSSFC  - Adjoint of the "vertical" interpolation
      !*                  for "SURFACE" data files.
      !*
      !*Author  : P. Koclas *CMC/AES  April 1996
      !*    -------------------
      !*
      !*     Purpose: based on surfc1dz to build the adjoint of the
      !*              vertical interpolation for SURFACE data files.
      !*
      implicit none
      INTEGER IPB,IPT
      REAL*8 ZRES
      REAL*8 ZWB,ZWT,zcon,zexp,zgamma,ZATV,ZTVG
      REAL*8 ZLEV,ZPT,ZPB,ZDADPS,ZDELPS,ZDELTV,ZGAMAZ,ZHHH
      INTEGER INDEX_HEADER,IK,nlev
      INTEGER INDEX_BODY,ITYP,INDEX_FAMILY
      real*8, pointer :: all_column(:),tt_column(:),hu_column(:),ps_column(:)
      REAL*8 :: dPdPsfc
      INTEGER, PARAMETER :: numFamily=4
      CHARACTER(len=2) :: list_family(numFamily),varLevel
      !C
      !C     Temperature lapse rate for extrapolation of gz below model surface
      !C
      zgamma = 0.0065d0 / GRAV
      zexp   = 1.0d0/(MPC_RGAS_DRY_AIR_R8*zgamma)
      !C
      !C*    1. Fill in COMMVO by using the adjoint of the "vertical" interpolation
      !C     .  ---------------------------------------------------------------
      list_family(1) = 'UA'
      list_family(2) = 'SF'
      list_family(3) = 'SC'
      list_family(4) = 'GP'

      FAMILY: do index_family=1,numFamily

         call obs_set_current_body_list(obsSpaceData, list_family(index_family))
         BODY: do
            index_body = obs_getBodyIndex(obsSpaceData)
            if (index_body < 0) exit BODY

            ! Process all data within the domain of the model
            ityp = obs_bodyElem_i(obsSpaceData,OBS_VNM,index_body)
            if ( ityp == bufr_nezd ) cycle BODY
            if(    (obs_bodyElem_i(obsSpaceData,OBS_VCO,index_body) == 1) &
                 .and. (obs_bodyElem_i(obsSpaceData,OBS_ASS,index_body) == 1) &
                 .and. (ityp == bufr_nets .or. ityp == bufr_neps  &
                 .or. ityp == bufr_nepn .or. ityp == bufr_ness  &
                 .or. ityp == bufr_neus .or. ityp == bufr_nevs  &
                 .or. obs_bodyElem_i(obsSpaceData,OBS_XTR,index_body) == 0) ) then

               if( ityp == bufr_neus .or. ityp == bufr_nevs ) then
                  varLevel = 'MM'
               else
                  varLevel = 'TH'
               end if
               nlev = col_getNumLev(column,varLevel)

               INDEX_HEADER = obs_bodyElem_i(obsSpaceData,OBS_HIND,INDEX_BODY)
               ITYP = obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY)
               IK   = obs_bodyElem_i(obsSpaceData,OBS_LYR,INDEX_BODY)
               ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY)
               ZHHH = ZLEV * GRAV
               IPT  = nlev - 1 + col_getOffsetFromVarno(columng,ityp)
               IPB  = IPT+1
               ZRES = obs_bodyElem_r(obsSpaceData,OBS_WORK,INDEX_BODY)
               IF (ITYP == BUFR_NETS .OR. ITYP == BUFR_NESS .OR.  &
                    ITYP == BUFR_NEUS .OR. ITYP == BUFR_NEVS ) THEN
                  if(ityp == bufr_ness) then
                     dPdPsfc = col_getPressureDeriv(columng,nlev,index_header,'TH')
                     tt_column  => col_getColumn(column,INDEX_HEADER,'TT')
                     hu_column  => col_getColumn(column,INDEX_HEADER,'HU')
                     ps_column  => col_getColumn(column,INDEX_HEADER,'P0')
                     call hutoes_ad(hu_column(nlev),  &
                          tt_column(nlev),  &
                          ps_column(1),     &
                          ZRES,             &
                          col_getElem(columng,nlev,INDEX_HEADER,'HU'),      &
                          col_getPressure(columng,nlev,INDEX_HEADER,'TH'),  &
                          dPdPsfc)
                  else
                     all_column => col_getColumn(column,INDEX_HEADER) 
                     all_column(IPB) = all_column(IPB) + ZRES
                  end if
               ELSEIF (ITYP == BUFR_NEPS .OR. ITYP == BUFR_NEPN) THEN
                  tt_column  => col_getColumn(column,INDEX_HEADER,'TT')
                  hu_column  => col_getColumn(column,INDEX_HEADER,'HU')
                  ps_column  => col_getColumn(column,INDEX_HEADER,'P0')
                  ZTVG  = columng%OLTV(1,nlev,INDEX_HEADER)*col_getElem(columng,nlev,INDEX_HEADER,'TT')
                  ZGAMAZ= ZGAMMA*(ZHHH-col_getGZsfc(columng,INDEX_HEADER))
                  ZCON  = ((ZTVG-ZGAMAZ)/ZTVG)
                  ZDELTV= (col_getElem(columng,1,INDEX_HEADER,'P0')*ZEXP*ZCON**(ZEXP-1))  &
                       *(ZGAMAZ/(ZTVG*ZTVG))
                  ZDELPS= ZCON**ZEXP
                  ZATV  = ZDELTV*ZRES
                  ps_column(1)    = ps_column(1) + ZDELPS*ZRES
                  tt_column(nlev) = tt_column(nlev)  &
                       + columng%OLTV(1,nlev,INDEX_HEADER)*ZATV
                  hu_column(nlev)= hu_column(nlev)   &
                       + columng%OLTV(2,nlev,INDEX_HEADER)*ZATV
               ELSE
                  ! not sure what this block of code is for, not present in nonlinear version (Buehner)
                  all_column => col_getColumn(column,INDEX_HEADER)
                  ps_column  => col_getColumn(column,INDEX_HEADER,'P0')
                  IPT  = IK + col_getOffsetFromVarno(columng,ityp)
                  IPB  = IPT+1
                  ZPT  = col_getHeight(columng,IK  ,INDEX_HEADER,varLevel)
                  ZPB  = col_getHeight(columng,IK+1,INDEX_HEADER,varLevel)
                  ZWB  = (ZPT-ZHHH)/(ZPT-ZPB)
                  ZWT  = 1.d0 - ZWB
                  !ccc ATTN ATTN ZDADPS EST A DEFINIR POUR UNE COORDONNEE Z
                  ZDADPS= 0.d0
                  all_column(IPB) = all_column(IPB) + ZWB*ZRES
                  all_column(IPT) = all_column(IPT) + ZWT*ZRES
                  ps_column(1)    = ps_column(1)    +         &
                       (col_getElem(columng,IPB,INDEX_HEADER) - col_getElem(columng,IPT,INDEX_HEADER))  &
                       *ZDADPS*ZRES
               END IF

            end if

         END DO BODY

      END DO FAMILY

      RETURN
    END subroutine oop_HTsf


    subroutine oop_HTto
      !
      !**s/r tovs_obs_ad  - Adjoint of computation of residuals to the tovs observations
      !
      !
      !author        : j. halle *cmda/aes  april 19, 2005
      !
      !revision 001  :
      !                S. Pellerin - ARMA, jan. 2009
      !                - call  to oop_get_radiance_ad
      !
      !    -------------------
      !     purpose:
      !
      implicit none

      if (.not.obs_famExist(obsSpaceData,'TO',local_mpi=.true.)) return

      !     1.   Getting the adjoint of the residuals
      !     .    ----------------------------------

      !     2.   Adjoint of computing radiance
      !     .    -----------------------------
      !
      call tvslin_rttov_ad(column,columng,obsSpaceData)


    end subroutine oop_HTto


    SUBROUTINE oop_HTro
      !*
      !* Purpose: Compute the adjoint operator for GPSRO observations.
      !*
      !*Author  : J. M. Aparicio Jan 2004
      !*Modified: J. M. Aparicio Dec 2012 adapt to accept bending angle data
      !*    -------------------
      use IndexListDepot_mod, only : struct_index_list
      implicit none

      REAL*8 DPJO0(ngpscvmx)
      REAL*8 DPJO1(ngpscvmx)

      REAL*8 zLat, Lat
      REAL*8 zAzm, Azm
      INTEGER IAZM, ISAT
      REAL*8 Rad, Geo, HNH1
      REAL*8 zP0, zMT

      REAL*8 ZINC, ZOER

      real*8, pointer :: tt_column(:),hu_column(:),ps_column(:)
      INTEGER IDATYP
      INTEGER JL, NGPSLEV
      integer :: index_header, index_body, iProfile
      type(struct_index_list), pointer :: local_current_list

      LOGICAL  ASSIM, LUSE

      INTEGER NH, NH1
      !C      WRITE(*,*)'ENTER oop_HTro'
      !C
      !C     * 1.  Initializations
      !C     *     ---------------
      !C
      NGPSLEV=col_getNumLev(column,'TH')
      !C
      !C    Loop over all header indices of the 'RO' family (Radio Occultation)
      !C
      ! Set the header list (start at the beginning of the list)
      call obs_set_current_header_list(obsSpaceData,'RO')
      !##$omp parallel default(shared) &
      !##$omp private(index_header,dpjo0,idatyp,assim,nh,local_current_list,index_body,luse) &
      !##$omp private(iProfile,zlat,irad,igeo,iazm,isat,rad,geo,zazm,zmt,lat) &
      !##$omp private(nh1,zinc,zoer,dpjo1) &
      !##$omp private(tt_column,hu_column,ps_column)
      nullify(local_current_list)
      HEADER: do
         INDEX_HEADER = obs_getHeaderIndex(obsSpaceData)
         if (INDEX_HEADER < 0) exit HEADER

         DPJO0 = 0.d0
         !C
         !C     * Process only refractivity data (codtyp 169)
         !C
         IDATYP = obs_headElem_i(obsSpaceData,OBS_ITY,INDEX_HEADER)
         DATYP: IF ( IDATYP == 169 ) THEN
            !C
            !C     *    Scan for requested data values of the profile, and count them
            !C
            ASSIM = .FALSE.
            NH = 0
            !C
            !C     *    Loop over all body indices for this index_header:
            !C     *    (start at the beginning of the list)
            !C
            call obs_set_current_body_list(obsSpaceData, INDEX_HEADER, &
                 current_list=local_current_list)
            BODY: do 
               index_body = obs_getBodyIndex(local_current_list)
               if (index_body < 0) exit BODY

               LUSE=( obs_bodyElem_i(obsSpaceData,OBS_ASS,INDEX_BODY) == 1 )
               IF ( LUSE ) THEN
                  ASSIM = .TRUE.
                  NH = NH + 1
               END IF
            END DO BODY
            !C
            !C     *    If assimilations are requested, prepare and apply the observation operator
            !C
            ASSIMILATE: IF (ASSIM) THEN
               iProfile=gps_iprofile_from_index(INDEX_HEADER)
               !C
               !C     *    Basic geometric variables of the profile:
               !C
               zLat = obs_headElem_r(obsSpaceData,OBS_LAT,INDEX_HEADER)
               IAZM = obs_headElem_i(obsSpaceData,OBS_AZA,INDEX_HEADER)
               ISAT = obs_headElem_i(obsSpaceData,OBS_SAT,INDEX_HEADER)
               Rad  = obs_headElem_r(obsSpaceData,OBS_TRAD,INDEX_HEADER)
               Geo  = obs_headElem_r(obsSpaceData,OBS_GEOI,INDEX_HEADER)
               zAzm = 0.01d0*IAZM / MPC_DEGREES_PER_RADIAN_R8
               zMT  = col_getGZsfc(columng,INDEX_HEADER)/RG
               Lat  = zLat * MPC_DEGREES_PER_RADIAN_R8
               !C
               !C     *       Perform the (H(xb)DX-Y')/S operation
               !C
               NH1 = 0
               !C
               !C     *       Loop over all body indices for this index_header:
               !C     *       (start at the beginning of the list)
               !C
               call obs_set_current_body_list(obsSpaceData, index_header, &
                    current_list=local_current_list)
               BODY_3: do 
                  index_body = obs_getBodyIndex(local_current_list)
                  if (index_body < 0) exit BODY_3

                  LUSE=( obs_bodyElem_i(obsSpaceData,OBS_ASS,INDEX_BODY) == 1 )
                  IF ( LUSE ) THEN
                     NH1 = NH1 + 1
                     !C
                     !C     *             Normalized increment
                     !C
                     ZINC = obs_bodyElem_r(obsSpaceData,OBS_WORK,INDEX_BODY)
                     !                     ZOER = obs_bodyElem_r(obsSpaceData,OBS_OER,INDEX_BODY)
                     !C
                     !C     *             O-F Tested criteria:
                     !C
                     DPJO1(1:2*NGPSLEV+1) = ZINC * gps_vRO_Jacobian(iProfile,NH1,:)
                     !C
                     !C     *             Accumulate the gradient of the observation cost function:
                     !C
                     DPJO0(1:2*NGPSLEV+1) = DPJO0(1:2*NGPSLEV+1) + DPJO1(1:2*NGPSLEV+1)
                  END IF
               END DO BODY_3
            END IF ASSIMILATE
         END IF DATYP
         !C
         !C     * Store H* (HX - Z)/SIGMA in COMMVO
         !C
         tt_column => col_getColumn(column,index_header,'TT')
         hu_column => col_getColumn(column,index_header,'HU')
         ps_column => col_getColumn(column,index_header,'P0')
         DO JL = 1, NGPSLEV
            tt_column(JL) = DPJO0(JL)
            hu_column(JL) = DPJO0(JL+NGPSLEV)
         END DO
         ps_column(1) = DPJO0(1+2*NGPSLEV)
      END DO HEADER
      !##$omp end parallel

      !C      WRITE(*,*)'EXIT oop_HTro'
      RETURN
    END subroutine oop_HTro


    SUBROUTINE oop_HTzp
      !*
      !***s/r AOBSZZZ  - Adjoint of the "vertical" interpolation in z
      !*                 for profiler data.
      !*
      !*Author  : J. St-James *CMDA/SMC  July 2003
      !*Revision :
      !*    -------------------
      !*
      !*     Purpose: based on vint3d to build the adjoint of the
      !*              vertical interpolation for profiler data.
      !*
      implicit none
      INTEGER IPB,IPT
      REAL*8 ZRES,ZDA1,ZDA2,ZDENO
      REAL*8 ZWB,ZWT
      REAL*8 ZLEV,ZPT,ZPB
      INTEGER INDEX_HEADER,IK,ITYP
      INTEGER INDEX_BODY
      real*8, pointer :: gz_column(:),all_column(:)
      character(len=2) :: varLevel
      !C
      !C     Process all data within the domain of the model
      !C
      call obs_set_current_body_list(obsSpaceData, 'PR')
      BODY: do
         index_body = obs_getBodyIndex(obsSpaceData)
         if (index_body < 0) exit BODY

         IF (   (obs_bodyElem_i(obsSpaceData,OBS_ASS,INDEX_BODY) == 1) &
              .AND. (obs_bodyElem_i(obsSpaceData,OBS_XTR,INDEX_BODY) == 0) &
              .AND. (obs_bodyElem_i(obsSpaceData,OBS_VCO,INDEX_BODY) == 1)  ) THEN
            INDEX_HEADER = obs_bodyElem_i(obsSpaceData,OBS_HIND,INDEX_BODY)
            ITYP = obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY)
            varLevel = vnl_varLevelFromVarnum(ityp)
            gz_column  => col_getColumn(column,INDEX_HEADER,'GZ',varLevel)
            all_column => col_getColumn(column,INDEX_HEADER)
            ZRES = obs_bodyElem_r(obsSpaceData,OBS_WORK,INDEX_BODY)
            ZLEV = obs_bodyElem_r(obsSpaceData,OBS_PPP,INDEX_BODY)
            IK   = obs_bodyElem_i(obsSpaceData,OBS_LYR,INDEX_BODY)
            IPT  = IK  + col_getOffsetFromVarno(columng,ityp)
            IPB  = IPT+1
            ZPT  = col_getHeight(columng,IK,INDEX_HEADER,varLevel)/RG
            ZPB  = col_getHeight(columng,IK+1,INDEX_HEADER,varLevel)/RG
            ZDENO= ZPT-ZPB
            ZWB  = (ZPT-ZLEV)/ZDENO
            ZWT  = 1.0D0 - ZWB

            ZDA1= (ZLEV-ZPB)/(ZDENO**2)
            ZDA2= (ZPT-ZLEV)/(ZDENO**2)
            !C
            gz_column(IK+1) = gz_column(IK+1) +    &
                 (col_getElem(columng,IPB,INDEX_HEADER)-col_getElem(columng,IPT,INDEX_HEADER))*ZDA2*ZRES/RG
            gz_column(IK) = gz_column(IK) +        &
                 (col_getElem(columng,IPB,INDEX_HEADER)-col_getElem(columng,IPT,INDEX_HEADER))*ZDA1*ZRES/RG
            all_column(IPB) = all_column(IPB) + ZWB*ZRES
            all_column(IPT) = all_column(IPT) + ZWT*ZRES

         END IF
      END DO BODY
      RETURN
    END subroutine oop_HTzp


    SUBROUTINE oop_HTgp
      !*
      !***s/r  -oop_HTgp Adjoint of TL routine oop_Hgp
      !*
      !*
      !*Author  : S. Macpherson *ARMA October 2012
      !
      !*Revisions:
      !
      ! S. Macpherson ARMA  14 Jan 2013
      !            - like oop_HTro, use OpenMP and Jacobian storage to speed up

      !*    -------------------
      !**    Purpose: Compute Ht*grad(Jo) for all GPS ZTD observations
      !
      !  NOTE:  ZTD Jacobians are computed and stored in oop_Hgp (first iter.)
      !
      !*
      implicit none

      REAL*8 DPJO0(ngpscvmx)
      REAL*8 JAC(ngpscvmx)
      ! 
      REAL*8 ZINC
      INTEGER JL, NFLEV, iztd
      integer :: index_header, index_body, icount
      LOGICAL ASSIM

      real*8, pointer :: tt_column(:),hu_column(:),ps_column(:)

      !      WRITE(*,*)'ENTER oop_HTgp'

      NFLEV  = col_getNumLev(columng,'TH')

      IF ( .not.vGPSZTD_lJac(1) ) THEN
         call utl_abort('oop_HTgp:ERROR: ZTD Jacobians not stored!')
      END IF

      ! loop over all header indices of the 'GP' family (GPS observations)
      ! Set the header list & start at the beginning of the list
      call obs_set_current_header_list(obsSpaceData,'GP')

      icount = 0

      HEADER: do
         index_header = obs_getHeaderIndex(obsSpaceData)
         if (index_header < 0) exit HEADER

         DPJO0(:) = 0.0D0
         JAC(:)   = 0.0D0
         !C
         !C       Scan for requested ZTD assimilation
         !C
         ASSIM = .FALSE.
         ! loop over all body indices (still in the 'GP' family)
         ! Set the body list & start at the beginning of the list)
         call obs_set_current_body_list(obsSpaceData, index_header)
         BODY: DO 
            index_body = obs_getBodyIndex(obsSpaceData)
            if (index_body < 0) exit BODY

            IF (   (obs_headElem_i(obsSpaceData,OBS_ITY,INDEX_HEADER) == 189) &
                 .and. (obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY) == BUFR_NEZD) &
                 .and. (obs_bodyElem_i(obsSpaceData,OBS_ASS,INDEX_BODY) == 1) ) then
               ASSIM = .TRUE.
            END IF
         END DO BODY
         !C
         IF (ASSIM) THEN

            icount = icount + 1
            iztd = gps_i_from_index(INDEX_HEADER)
            if ( iztd < 1 .or. iztd > numGPSZTD ) then
               call utl_abort('oop_HTgp: ERROR: index from gps_i_from_index() is out of range!')
            end if

            DO JL = 1, 2*NFLEV+1
               JAC(JL) = vGPSZTD_Jacobian(iztd,JL)
            END DO
            !C
            !C          Get Ht*grad(Index_header) = Ht*(H'dx - d)/sigma_o^2
            !C
            ! loop over all body indices (still in the 'GP' family)
            ! Start at the beginning of the list)
            call obs_set_current_body_list(obsSpaceData, index_header)
            BODY_2: do 
               index_body = obs_getBodyIndex(obsSpaceData)
               if (index_body < 0) exit BODY_2
               if (   (obs_headElem_i(obsSpaceData,OBS_ITY,INDEX_HEADER) == 189) &
                    .and. (obs_bodyElem_i(obsSpaceData,OBS_VNM,INDEX_BODY) == BUFR_NEZD) &
                    .and. (obs_bodyElem_i(obsSpaceData,OBS_ASS,INDEX_BODY) == 1) ) then
                  ZINC = obs_bodyElem_r(obsSpaceData,OBS_WORK,INDEX_BODY)
                  !C     *       Accumulate the gradient of the observation cost function
                  DPJO0(1:2*NFLEV+1) = ZINC * vGPSZTD_Jacobian(iztd,:)
               end if
            END DO BODY_2
            !c
            !C      *   Store Ht*grad(Index_header) in COMMVO
            !c
            tt_column => col_getColumn(column,index_header,'TT')
            hu_column => col_getColumn(column,index_header,'HU')
            ps_column => col_getColumn(column,index_header,'P0')
            DO JL = 1, NFLEV
               tt_column(JL) = DPJO0(JL)
               hu_column(JL) = DPJO0(JL+NFLEV)
            END DO
            ps_column(1) = DPJO0(2*NFLEV+1)

         END IF ! ASSIM


      END DO HEADER

      !      WRITE(*,*) 'oop_HTgp: Number of ZTD data locations processed = ', icount

      !      WRITE(*,*)'EXIT oop_HTgp'

      RETURN
    END subroutine oop_HTgp


    subroutine oop_HTchm
      !*
      !***s/r  - oop_HTchm: Adjoint of TL routine oda_Hchm
      !*
      !*
      !*Author: M. Sitwell, ARQI/AQRD, Aug 2015
      !         - Reduced to a call of chm_observation_operators
      !           which contains much of the original version by Ping Du (CMDA/MSC) and
      !           Y. Rochon (ARQI/AQRD), Jan 2015 (partially based on corresponding
      !           pre-EnVar routine by Y.J. Rochon and Y. Yang, July 2005 to Feb 2013),
      !           with further original changes by Y. Rochon and M. Sitwell, 2015.
      !
      !**   Purpose: Compute H^T * R^-1 (OmP-Hdx) for all CH observations  
      !
      !-----------------------------------------------------------------------------------

      implicit none
      
      if (.not.obs_famExist(obsSpaceData,'CH',local_mpi=.true.)) return
      
      call chm_observation_operators(columng,obsSpaceData,kmode=3,columnInc_opt=column) ! kmode=3 for adjoint of the tangent linear operator

    end subroutine oop_HTchm

  end subroutine oop_Had

  function HUtoES(hu,tt,pressure) result(es)
    !
    ! Purpose:
    !          to calculate the dew point depression from specific
    !          humidity, temperature and pressure.  No ice phase
    !          is permitted and the pressure vector is given.
    !
    implicit none
    real(8), intent(in) :: hu, tt, pressure
    real(8) :: husat, td, es

    ! get the saturated vapor pressure from lq (log of specific humidity)
    husat = foefq8(max(hu,1.d-12),pressure)

    ! now the dewpoint temperature
    td = fotw8(husat)

    ! finally the dewpoint depression
    es = min(tt-td,MPC_MAXIMUM_ES_R8)

  end function HUtoES

  function HUtoES_tl(HU_inc,TT_inc,P0_inc,HU_trial,PRES_trial,dPdPsfc) result(ES_inc)
    !
    ! Purpose: TLM VERSION
    !          to calculate the dew point depression from specific
    !          humidity, temperature and pressure.  No ice phase
    !          is permitted and the pressure vector is given.
    !
    implicit none
    REAL(8), intent(in) :: HU_inc, TT_inc, P0_inc, HU_trial, PRES_trial, dPdPsfc
    REAL(8) :: ZE, ZTD, dTDdE, ZQBRANCH, ES_inc
    REAL(8) :: dESdLQ, dESdTT, dESdP0

    dESdTT = 1.0d0

    !- Forward calculations of saturation vapour pressure and dewpoint temperature
    !  and adjoint of vapour pressure from adjoint of dewpoint temperature
    ZE   = FOEFQ8(HU_trial, PRES_trial)
    ZTD  = FOTW8 (ZE)
    dTDdE= FODTW8(ZTD,ZE)

    !- adjoint of temp. specific humidity and surface pressure due to changes in vapour pressure
    ZQBRANCH = FQBRANCH(HU_trial)

    dESdLQ = - ZQBRANCH*FOEFQA(1.0d0,dTDdE,HU_trial,PRES_trial)

    dESdP0 = - ZQBRANCH*FOEFQPSA(1.0d0,dTDdE,HU_trial,dPdPsfc)-  &
               (1.D0-ZQBRANCH)*(dTDdE*dPdPsfc)

    ES_inc =  dESdLQ*HU_inc/max(HU_trial,1.d-12) + dESdP0*P0_inc + dESdTT*TT_inc

  end function HUtoES_tl

  subroutine HUtoES_ad(HU_inc,TT_inc,P0_inc,ES_inc,HU_trial,PRES_trial,dPdPsfc)
    !
    ! Purpose: ADJOINT VERSION
    !          to calculate the dew point depression from specific
    !          humidity, temperature and pressure.  No ice phase
    !          is permitted and the pressure vector is given.
    !
    implicit none
    REAL(8), intent(inout) :: HU_inc,TT_inc,P0_inc
    REAL(8), intent(in)  :: ES_inc,HU_trial,PRES_trial,dPdPsfc
    REAL(8) :: ZE,ZTD,dTDdE,ZQBRANCH
    REAL(8) :: dESdLQ,dESdTT,dESdP0

    dESdTT = 1.0d0
   
    !- Forward calculations of saturation vapour pressure and dewpoint temperature
    !  and adjoint of vapour pressure from adjoint of dewpoint temperature
    ZE = FOEFQ8(HU_trial, PRES_trial)

    ZTD=FOTW8(ZE)
    dTDdE=FODTW8(ZTD,ZE)

    !- adjoint of temp. specific humidity and surface pressure due to changes in vapour pressure
    ZQBRANCH = FQBRANCH(HU_trial)
    dESdLQ = - ZQBRANCH*FOEFQA(1.0d0,dTDdE,HU_trial,PRES_trial)

    dESdP0 = - ZQBRANCH*FOEFQPSA(1.0d0,dTDdE,HU_trial,dPdPsfc)-  &
               (1.D0-ZQBRANCH)*(dTDdE*dPdPsfc)

    ! TLM: ES_inc =  dESdLQ*HU_inc/HU_trial + dESdP0*P0_inc + dESdTT*TT_inc
    ! ADJOINT:
    HU_inc = HU_inc + dESdLQ*ES_inc/max(HU_trial,1.d-12)
    P0_inc = P0_inc + dESdP0*ES_inc
    TT_inc = TT_inc + dESdTT*ES_inc

  end subroutine HUtoES_ad

end module obsOperators_mod
