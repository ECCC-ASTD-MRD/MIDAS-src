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
!-------------------------------------- LICENCE end --------------------------------------

module bgckmicrowave_mod
  ! MODULE bgckmicrowave_mod (prefix='mwbg' category='1. High-level functionality')
  !
  ! :Purpose: Variables for microwave background check and quality control.
  !
  use burp_module
  use MathPhysConstants_mod
  use utilities_mod

  implicit none
  save
  private

  ! Public functions/subroutines

  ! common functions/subroutines
  public :: mwbg_getData
  public :: mwbg_updateBurp 
  public :: mwbg_readStatTovs
  public :: mwbg_qcStats
  public :: mwbg_readGeophysicFieldsAndInterpolate
  public :: mwbg_findSatelliteIndex
  public :: mwbg_allocate1DCharacterArray
  public :: mwbg_allocate1DIntegerArray
  public :: mwbg_allocate2DIntegerArray
  public :: mwbg_allocate3DIntegerArray
  public :: mwbg_allocate1DLogicalArray
  public :: mwbg_allocate2DLogicalArray
  public :: mwbg_allocate1DRealArray
  public :: mwbg_allocate2DRealArray

  ! ATMS specific functions/subroutines
  public :: mwbg_landIceMaskAtms
  public :: mwbg_firstQcCheckAtms
  public :: mwbg_grossValueCheck
  public :: mwbg_nrlFilterAtms
  public :: mwbg_flagDataUsingNrlCriteria
  public :: mwbg_reviewAllcriteriaforFinalFlags
  public :: mwbg_tovCheckAtms
  ! AMSUA specific functions/subroutines
  public :: mwbg_tovCheckAmsua
  ! public variables
  public :: mwbg_debug
  public :: mwbg_clwQcThreshold
  public :: mwbg_useUnbiasedObsForClw 
  public :: mwbg_allowStateDepSigmaObs
  public :: mwbg_maxNumChan
  !
  public :: mwbg_maxNumSat
  public :: mwbg_maxNumTest
  public :: mwbg_maxNumObs
  public :: mwbg_maxScanAngle
  public :: mwbg_atmsMaxNumChan
  public :: mwbg_realMisg
  public :: mwbg_intMisg

  real    :: mwbg_clwQcThreshold
  logical :: mwbg_debug
  logical :: mwbg_useUnbiasedObsForClw 
  logical :: mwbg_allowStateDepSigmaObs

  integer, parameter :: mwbg_maxNumChan = 50
  integer, parameter :: mwbg_maxNumSat = 9
  integer, parameter :: mwbg_maxNumTest = 15
  integer, PARAMETER :: mwbg_maxNumObs = 3000
  integer, parameter :: mwbg_atmsMaxNumChan=22
  integer, parameter :: mwbg_maxScanAngle=96
  real, parameter    :: mwbg_realMisg=9.9e09
  integer, parameter :: mwbg_intMisg = -1

  ! Module variable

  integer, parameter :: mwbg_atmsNumSfcSensitiveChannel = 6

  ! Upper limit for CLW (kg/m**2) for Tb rejection over water
  real,   parameter :: clw_atms_nrl_LTrej=0.175      ! lower trop chans 1-6, 16-20
  real,   parameter :: clw_atms_nrl_UTrej=0.2        ! upper trop chans 7-9, 21-22
  ! Other NRL thresholds
  real,   parameter :: scatec_atms_nrl_LTrej=9.0     ! lower trop chans 1-6, 16-22
  real,   parameter :: scatec_atms_nrl_UTrej=18.0    ! upper trop chans 7-9
  real,   parameter :: scatbg_atms_nrl_LTrej=10.0    ! lower trop chans 1-6
  real,   parameter :: scatbg_atms_nrl_UTrej=15.0    ! upper trop chans 7-9
  real,   parameter :: mean_Tb_183Ghz_min=240.0      ! min. value for Mean(Tb) chans. 18-22 
  
contains

  !--------------------------------------------------------------------------
  !  ISRCHEQR function
  !--------------------------------------------------------------------------
  integer function ISRCHEQR (KLIST, KLEN, KENTRY)
    ! OBJET          Rechercher un element dans une liste (valeurs reelles).
    ! APPEL          INDX = ISRCHEQR (KLIST, KLEN, KENTRY)
    ! ARGUMENTS      - indx    - output -  position de l'element recherche:
    !                                   =0, element introuvable,
    !                                   >0, position de l'element trouve,
    !                - klist   - input  -  la liste
    !                - klen    - input  -  longueur de la liste
    !                - kentry  - input  -  l'element recherche

    implicit none

    integer  KLEN, JI

    real  KLIST(KLEN)
    real  KENTRY

    ISRCHEQR = 0
    do JI=1,KLEN
       if ( NINT(KLIST(JI)) .EQ. NINT(KENTRY) ) then
          ISRCHEQR = JI
          return
       end if
    end do

  end function ISRCHEQR

  !--------------------------------------------------------------------------
  !  ISRCHEQR function
  !--------------------------------------------------------------------------
  function ISRCHEQI (KLIST, KLEN, KENTRY) result(ISRCHEQI_out)
    !OBJET          Rechercher un element dans une liste (valeurs entieres).
    !ARGUMENTS      - indx    - output -  position de l'element recherche:
    !                                   =0, element introuvable,
    !                                   >0, position de l'element trouve,
    !               - klist   - input  -  la liste
    !               - klen    - input  -  longueur de la liste
    !               - kentry  - input  -  l'element recherche

    implicit none

    integer :: ISRCHEQI_out
    integer  KLEN, JI

    integer  KLIST(KLEN)
    integer  KENTRY

    ISRCHEQI_out = 0
    do JI=1,KLEN
       if ( KLIST(JI) .EQ. KENTRY ) then
          ISRCHEQI_out = JI
          return
       end if
    end do

  end function ISRCHEQI

  !--------------------------------------------------------------------------
  !  extractParamForGrodyRun
  !--------------------------------------------------------------------------  
  subroutine extractParamForGrodyRun(KCANO, ptbo, ptbomp, ptbcor, KNT, KNO, pmisg, &
                                     tb23,   tb31,   tb50,   tb53,   tb89, &
                                     tb23_P, tb31_P, tb50_P, tb53_P, tb89_P)

    !:Purpose: Compute  Grody parameters by extracting tb for required channels:
    !          - 23 Ghz = AMSU-A 1 = channel #28
    !          - 31 Ghz = AMSU-A 2 = channel #29
    !          - 50 Ghz = AMSU-A 3 = channel #30
    !          - 53 Ghz = AMSU-A 5 = channel #32
    !          - 89 Ghz = AMSU-A15 = channel #42
    ! Arguments
    integer,     intent(in)               :: KCANO(KNO,KNT)         ! observations channels
    real,        intent(in)               :: ptbo(KNO,KNT)          ! radiances
    real,        intent(in)               :: ptbomp(KNO,KNT)        ! radiances o-p
    real,        intent(in)               :: ptbcor(KNO,KNT)        ! correction aux radiances
    integer,     intent(in)               :: KNO                    ! nombre de canaux des observations 
    integer,     intent(in)               :: KNT                    ! nombre de tovs
    real,        intent(in)               :: pmisg                  ! value for missing data
    real,        intent(out)              :: tb23(KNT)              ! radiance frequence 23 Ghz   
    real,        intent(out)              :: tb31(KNT)              ! radiance frequence 31 Ghz
    real,        intent(out)              :: tb50(KNT)              ! radiance frequence 50 Ghz  
    real,        intent(out)              :: tb53(KNT)              ! radiance frequence 53 Ghz  
    real,        intent(out)              :: tb89(KNT)              ! radiance frequence 89 Ghz  
    real,        intent(out)              :: tb23_P(KNT)            ! radiance frequence 23 Ghz   
    real,        intent(out)              :: tb31_P(KNT)            ! radiance frequence 31 Ghz
    real,        intent(out)              :: tb50_P(KNT)            ! radiance frequence 50 Ghz  
    real,        intent(out)              :: tb53_P(KNT)            ! radiance frequence 53 Ghz  
    real,        intent(out)              :: tb89_P(KNT)            ! radiance frequence 89 Ghz        

    ! Locals
    integer                               :: channelval
    integer                               :: nDataIndex
    integer                               :: nChannelIndex

    do nDataIndex=1,KNT
      do nChannelIndex=1,KNO
        channelval = KCANO(nChannelIndex,nDataIndex)
        if ( ptbo(nChannelIndex,nDataIndex) .ne. pmisg ) then
          if ( ptbcor(nChannelIndex,nDataIndex) .ne. pmisg ) then
            if ( channelval .eq. 28 ) tb23(nDataIndex) = ptbo(nChannelIndex,nDataIndex) &
                 - ptbcor(nChannelIndex,nDataIndex)
            if ( channelval .eq. 29 ) tb31(nDataIndex) = ptbo(nChannelIndex,nDataIndex) &
                 - ptbcor(nChannelIndex,nDataIndex)
            if ( channelval .eq. 30 ) tb50(nDataIndex) = ptbo(nChannelIndex,nDataIndex) &
                 - ptbcor(nChannelIndex,nDataIndex)
            if ( channelval .eq. 32 ) tb53(nDataIndex) = ptbo(nChannelIndex,nDataIndex) &
                 - ptbcor(nChannelIndex,nDataIndex)
            if ( channelval .eq. 42 ) tb89(nDataIndex) = ptbo(nChannelIndex,nDataIndex) &
                 - ptbcor(nChannelIndex,nDataIndex)
          else
            if ( channelval .eq. 28 ) tb23(nDataIndex) = ptbo(nChannelIndex,nDataIndex)
            if ( channelval .eq. 29 ) tb31(nDataIndex) = ptbo(nChannelIndex,nDataIndex)
            if ( channelval .eq. 30 ) tb50(nDataIndex) = ptbo(nChannelIndex,nDataIndex)
            if ( channelval .eq. 32 ) tb53(nDataIndex) = ptbo(nChannelIndex,nDataIndex)
            if ( channelval .eq. 42 ) tb89(nDataIndex) = ptbo(nChannelIndex,nDataIndex)
          end if

          if ( channelval .eq. 28 ) tb23_P(nDataIndex) = ptbo(nChannelIndex,nDataIndex) - &
                                                         ptbomp(nChannelIndex,nDataIndex)
          if ( channelval .eq. 29 ) tb31_P(nDataIndex) = ptbo(nChannelIndex,nDataIndex) - &
                                                         ptbomp(nChannelIndex,nDataIndex)
          if ( channelval .eq. 30 ) tb50_P(nDataIndex) = ptbo(nChannelIndex,nDataIndex) - &
                                                         ptbomp(nChannelIndex,nDataIndex)
          if ( channelval .eq. 32 ) tb53_P(nDataIndex) = ptbo(nChannelIndex,nDataIndex) - &
                                                         ptbomp(nChannelIndex,nDataIndex)
          if ( channelval .eq. 42 ) tb89_P(nDataIndex) = ptbo(nChannelIndex,nDataIndex) - &
                                                         ptbomp(nChannelIndex,nDataIndex)
        else
          if ( channelval .eq. 28 ) tb23(nDataIndex) = 0.
          if ( channelval .eq. 29 ) tb31(nDataIndex) = 0.
          if ( channelval .eq. 30 ) tb50(nDataIndex) = 0.
          if ( channelval .eq. 32 ) tb53(nDataIndex) = 0.
          if ( channelval .eq. 42 ) tb89(nDataIndex) = 0.

          if ( channelval .eq. 28 ) tb23_P(nDataIndex) = 0.  
          if ( channelval .eq. 29 ) tb31_P(nDataIndex) = 0. 
          if ( channelval .eq. 30 ) tb50_P(nDataIndex) = 0. 
          if ( channelval .eq. 32 ) tb53_P(nDataIndex) = 0. 
          if ( channelval .eq. 42 ) tb89_P(nDataIndex) = 0. 
        end if
      end do
    end do

    end subroutine extractParamForGrodyRun

    !--------------------------------------------------------------------------
    !  amsuaTest10RttovRejectCheck
    !--------------------------------------------------------------------------
    subroutine amsuaTest10RttovRejectCheck(KCANO, KNOSAT, KNO, KNT, RESETQC, STNID, KMARQ, ICHECK, rejectionCodArray)

    !:Purpose:               10) test 10: RTTOV reject check (single)
    !                        Rejected datum flag has bit #9 on.

    ! Arguments
    integer,     intent(in)                :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)                :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)                :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)                :: KNT                            ! nombre de tovs
    logical,     intent(in)                :: RESETQC                        ! yes or not reset QC flag
    character *9, intent(in)               :: STNID                          ! identificateur du satellite
    integer,     intent(out)               :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(out)               :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(out)               :: rejectionCodArray(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 
    ! Locals
    integer                                :: channelval
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex
    integer                                :: IBIT

    if (.NOT.RESETQC) then
      testIndex = 10
      do nDataIndex=1,KNT
        do nChannelIndex=1,KNO
          if ( KCANO(nChannelIndex,nDataIndex) .NE. 20 ) then
            IBIT = AND(KMARQ(nChannelIndex,nDataIndex), 2**9)
            if ( IBIT .NE. 0  ) then
              ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**7)
              rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                   rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
              if ( mwbg_DEBUG ) then
                write(*,*)STNID(2:9),' RTTOV REJECT.', &
                          'CHANNEL=', KCANO(nChannelIndex,nDataIndex), &
                          ' IMARQ= ',KMARQ(nChannelIndex,nDataIndex)
              end if
            end if
          end if
        end do
      end do
    end if

  end subroutine amsuaTest10RttovRejectCheck

  !--------------------------------------------------------------------------
  !  amsuaTest1TopographyCheck
  !--------------------------------------------------------------------------
  subroutine amsuaTest1TopographyCheck (KCANO, KNOSAT, KNO, KNT, STNID, MTINTRP, KMARQ, ICHECK, rejectionCodArray)

    !:Purpose:               1) test 1: Topography check (partial)
    !                        Channel 6 is rejected for topography >  250m.
    !                        Channel 7 is rejected for topography > 2000m.

    ! Arguments
    integer,     intent(in)                :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)                :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)                :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)                :: KNT                            ! nombre de tovs
    character *9,intent(in)                :: STNID                          ! identificateur du satellite
    real,        intent(in)                :: MTINTRP(KNT)                   ! topo aux point d'obs
    integer,     intent(out)               :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(out)               :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(out)               :: rejectionCodArray(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 
    ! Locals
    integer                                :: channelval
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex

    testIndex = 1
    do nDataIndex=1,KNT
      do nChannelIndex=1,KNO
        if ( KCANO(nChannelIndex,nDataIndex) .EQ. 33 ) then
          if ( MTINTRP(nDataIndex) .GE. 250.  ) then
            ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**18)
            rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                 rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
            if ( mwbg_DEBUG ) then
              write(*,*)STNID(2:9),' TOPOGRAPHY REJECT.', &
                        'CHANNEL=', KCANO(nChannelIndex,nDataIndex), &
                        ' TOPO= ',MTINTRP(nDataIndex)
            end if
          end if
        else if ( KCANO(nChannelIndex,nDataIndex) .EQ. 34 ) then
          if ( MTINTRP(nDataIndex) .GE. 2000.  ) then
            ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**18)
            rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                 rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
            if ( mwbg_DEBUG ) then
              write(*,*)STNID(2:9),' TOPOGRAPHY REJECT.', &
                        'CHANNEL=', KCANO(nChannelIndex,nDataIndex), &
                        ' TOPO= ',MTINTRP(nDataIndex)
            end if
          end if
        end if
      end do
    end do

  end subroutine amsuaTest1TopographyCheck

  !--------------------------------------------------------------------------
  ! amsuaTest2LandSeaQualifierCheck 
  !--------------------------------------------------------------------------
  subroutine amsuaTest2LandSeaQualifierCheck (KCANO, KNOSAT, KNO, KNT, STNID, KTERMER, KMARQ, ICHECK, rejectionCodArray)

    !:Purpose:                      2) test 2: "Land/sea qualifier" code check (full)
    !                                  allowed values are: 0, land,
    !                                  1, sea,
    !                                  2, coast.

    ! Arguments
    integer,     intent(in)                :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)                :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)                :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)                :: KNT                            ! nombre de tovs
    character *9,intent(in)                :: STNID                          ! identificateur du satellite
    integer,     intent(in)                :: KTERMER(KNT)                   ! land sea qualifier
    integer,     intent(out)               :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(out)               :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(out)               :: rejectionCodArray(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 
    ! Locals
    integer                                :: channelval
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex

  
    testIndex = 2
    do nDataIndex=1,KNT
      if ( KTERMER(nDataIndex) .LT.  0  .OR. &
          KTERMER(nDataIndex) .GT.  2        ) then
        do nChannelIndex=1,KNO
          ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
          KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
          KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**7)
          rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
            rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) + 1
        end do
        if ( mwbg_DEBUG ) then
          write(*,*) STNID(2:9),'LAND/SEA QUALifIER CODE', &
                   ' REJECT. KTERMER=', KTERMER(nDataIndex)
        end if
      end if
    end do

  end subroutine amsuaTest2LandSeaQualifierCheck

  !--------------------------------------------------------------------------
  !  amsuaTest3TerrainTypeCheck
  !--------------------------------------------------------------------------
  subroutine amsuaTest3TerrainTypeCheck (KCANO, KNOSAT, KNO, KNT, STNID, ITERRAIN, mwbg_intMisg, KMARQ, ICHECK, rejectionCodArray)

    !:Purpose:                     3) test 3: "Terrain type" code check (full)
    !                                 allowed values are: -1, missing,
    !                                 0, sea-ice,
    !                                 1, snow on land.

    ! Arguments
    integer,     intent(in)               :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)               :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)               :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)               :: KNT                            ! nombre de tovs
    integer,     intent(in)               :: mwbg_intMisg                        ! terrain type missing value
    character *9,intent(in)               :: STNID                          ! identificateur du satellite
    integer,     intent(in)               :: ITERRAIN(KNT)                  ! terrain type
    integer,     intent(out)              :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(out)              :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(out)              :: rejectionCodArray(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 
    ! Locals
    integer                               :: channelval
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex

 
    testIndex = 3
    do nDataIndex=1,KNT
      if ( ITERRAIN(nDataIndex) .NE.  mwbg_intMisg ) then
        if ( ITERRAIN(nDataIndex) .LT.  0  .OR. &
           ITERRAIN(nDataIndex) .GT.  1        ) then
          do nChannelIndex=1,KNO
            ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**7)
            rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
              rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) + 1
          end do
          if ( mwbg_debug ) then
            write(*,*)STNID(2:9),'TERRAIN type CODE', &
                     ' REJECT. TERRAIN=', ITERRAIN(nDataIndex)
          end if
        end if
      end if
    end do

  end subroutine amsuaTest3TerrainTypeCheck

  !--------------------------------------------------------------------------
  ! amsuaTest4FieldOfViewCheck 
  !--------------------------------------------------------------------------
  subroutine amsuaTest4FieldOfViewCheck (KCANO, KNOSAT, KNO, KNT, STNID, ISCNPOS, mwbg_maxScanAngleAMSU, KMARQ, ICHECK, rejectionCodArray)

    !:Purpose:                          4) test 4: Field of view number check (full)
    !                                      Field of view acceptable range is [1,mwbg_maxScanAngleAMSU]  for AMSU footprints.
    ! Arguments
    integer,     intent(in)               :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)               :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)               :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)               :: KNT                            ! nombre de tovs
    character *9,intent(in)               :: STNID                          ! identificateur du satellite
    integer,     intent(in)               :: ISCNPOS(KNT)                   ! position sur le "scan" 
    integer,     intent(in)               :: mwbg_maxScanAngleAMSU                     ! max scan angle 
    integer,     intent(out)              :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(out)              :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(out)              :: rejectionCodArray(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 
    ! Locals
    integer                               :: channelval
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex

    testIndex = 4
    do nDataIndex=1,KNT
      do nChannelIndex=1,KNO
        if ( ISCNPOS(nDataIndex) .LT. 1 .OR. &
            ISCNPOS(nDataIndex) .GT. mwbg_maxScanAngleAMSU ) then
          ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
          KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
          KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**7)
          rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
            rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) + 1
          if ( mwbg_debug ) then
            write(*,*)STNID(2:9),'FIELD OF VIEW NUMBER', &
                      ' REJECT. FIELD OF VIEW= ', ISCNPOS(nDataIndex)
          end if
        end if
      end do
    end do
  
  end subroutine amsuaTest4FieldOfViewCheck 
  
  !--------------------------------------------------------------------------
  ! amsuaTest5ZenithAngleCheck 
  !--------------------------------------------------------------------------
  subroutine amsuaTest5ZenithAngleCheck (KCANO, KNOSAT, KNO, KNT, STNID, SATZEN, PMISG, KMARQ, ICHECK, rejectionCodArray)

    !:Purpose:                   5) test 5: Satellite zenith angle check (full)
    !                               Satellite zenith angle acceptable range is [0.,60.].
    ! Arguments
    integer,     intent(in)               :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)               :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)               :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)               :: KNT                            ! nombre de tovs
    character *9,intent(in)               :: STNID                          ! identificateur du satellite
    real,        intent(in)               :: SATZEN(KNT)                    ! satellite zenith angle 
    real,        intent(in)               :: PMISG                          ! missing value 
    integer,     intent(out)              :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(out)              :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(out)              :: rejectionCodArray(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 
    ! Locals
    integer                               :: channelval
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex


    testIndex = 5
    do nDataIndex=1,KNT
      if ( SATZEN(nDataIndex) .NE.  PMISG ) then
        if ( SATZEN(nDataIndex) .LT.  0.  .OR. &
           SATZEN(nDataIndex) .GT. 60.       ) then
          do nChannelIndex=1,KNO
            ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**7)
            rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
               rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) + 1
          end do
          if ( mwbg_debug ) then
            write(*,*)STNID(2:9),' SATELLITE ZENITH ANGLE', &
                      ' REJECT. SATZEN= ', &
                      SATZEN(nDataIndex)
          end if
        end if
      end if
    end do

  end subroutine amsuaTest5ZenithAngleCheck 

  !--------------------------------------------------------------------------
  ! amsuaTest6ZenAngleAndFovConsistencyCheck 
  !--------------------------------------------------------------------------
  subroutine amsuaTest6ZenAngleAndFovConsistencyCheck (KCANO, KNOSAT, KNO, KNT, STNID, SATZEN,  PMISG, &
                                                  ISCNPOS, mwbg_intMisg, mwbg_maxScanAngleAMSU, KMARQ, ICHECK, rejectionCodArray)

    !:Purpose:                            6) test 6: "Sat. zenith angle"/"field of view" consistency check.  (full)
    !                                        Acceptable difference between "Satellite zenith angle"  and
    !                                       "approximate angle computed from field of view number" is 1.8 degrees.

    ! Arguments
    integer,     intent(in)               :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)               :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)               :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)               :: KNT                            ! nombre de tovs
    character *9,intent(in)               :: STNID                          ! identificateur du satellite
    real,        intent(in)               :: SATZEN(KNT)                    ! satellite zenith angle 
    real,        intent(in)               :: PMISG                          ! missing value 
    integer,     intent(in)               :: ISCNPOS(KNT)                   ! position sur le "scan" 
    integer,     intent(in)               :: mwbg_intMisg                        ! missing value 
    integer,     intent(in)               :: mwbg_maxScanAngleAMSU                     ! max scan angle 
    integer,     intent(out)              :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(out)              :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(out)              :: rejectionCodArray(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 
    ! Locals
    integer                               :: channelval
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex
    real                                  :: ZANGL
    real                                  :: APPROXIM 
    real                                  :: ANGDif 

    ZANGL   = 3.92
    testIndex = 6
    do nDataIndex=1,KNT
      if ( SATZEN (nDataIndex) .NE.  PMISG   .AND. &
         ISCNPOS(nDataIndex) .NE.  mwbg_intMisg       ) then
        APPROXIM = ABS((ISCNPOS(nDataIndex)-mwbg_maxScanAngleAMSU/2.-0.5)*ZANGL)
        ANGDif = ABS(SATZEN (nDataIndex)-APPROXIM)
        if ( ANGDif .GT. 1.8 ) then 
          do nChannelIndex=1,KNO
            ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**7)
            rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
               rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) + 1
          end do
          if ( mwbg_debug ) then
            write(*,*)STNID(2:9),' ANGLE/FIELD OF VIEW', &
                      ' INCONSISTENCY REJECT. SATZEN= ', &
                      SATZEN(nDataIndex), ' FIELD OF VIEW= ',ISCNPOS(nDataIndex), &
                      ' ANGDif= ',ANGDif  
          end if
        end if
      end if
    end do

  end subroutine amsuaTest6ZenAngleAndFovConsistencyCheck

  !--------------------------------------------------------------------------
  !  amsuaTest7landSeaQualifyerAndModelLandSeaConsistencyCheck
  !--------------------------------------------------------------------------
  subroutine amsuaTest7landSeaQualifyerAndModelLandSeaConsistencyCheck (KCANO, KNOSAT, KNO, KNT, STNID, MGINTRP, KTERMER, &
                                                                        KMARQ, ICHECK, rejectionCodArray)

    !:Purpose:                    7) test 7: "Land/sea qual."/"model land/sea" consistency check.    (full)
    !                                Acceptable conditions are:
    !                                - both over ocean (ktermer=1; mg<0.01), new threshold 0.20, jh dec 2000,
    !                                - both over land  (ktermer=0; mg>0.80), new threshold 0.50, jh dec 2000.
    !                                - Other conditions are unacceptable.


    ! Arguments
    integer,     intent(in)                :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)                :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)                :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)                :: KNT                            ! nombre de tovs
    character *9,intent(in)                :: STNID                          ! identificateur du satellite
    real,        intent(in)                :: MGINTRP(KNT)                   ! glace mer 
    integer,     intent(in)                :: KTERMER(KNT)                   ! land sea qualifyer 
    integer,     intent(out)               :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(out)               :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(out)               :: rejectionCodArray(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 
    ! Locals
    integer                                :: channelval
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex

    testIndex = 7
    do nDataIndex=1,KNT
      if ( KTERMER (nDataIndex) .NE.  mwbg_intMisg  ) then
        if     ( KTERMER(nDataIndex) .EQ. 1       .AND. &
                MGINTRP(nDataIndex) .LT. 0.20          ) then
        elseif ( KTERMER(nDataIndex) .EQ. 0       .AND. &
                MGINTRP(nDataIndex) .GT. 0.50          ) then
        else
          do nChannelIndex=1,KNO
            ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**7)
            rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
               rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) + 1
          end do
          if ( mwbg_debug ) then
            write(*,*)STNID(2:9),' LAND/SEA QUALifIER', &
                      ' INCONSISTENCY REJECT. KTERMER= ', &
                      KTERMER(nDataIndex), ' MODEL MASK= ',MGINTRP(nDataIndex)
          end if
        end if
      end if
    end do

  end subroutine amsuaTest7landSeaQualifyerAndModelLandSeaConsistencyCheck 

  !--------------------------------------------------------------------------
  !  amsuaTest9UncorrectedTbCheck
  !--------------------------------------------------------------------------
  subroutine amsuaTest9UncorrectedTbCheck (KCANO, KNOSAT, KNO, KNT, STNID, RESETQC, KMARQ, ICHECK, rejectionCodArray)

    !:Purpose:                  9) test 9: Uncorrected Tb check (single)
    !                              Uncorrected datum (flag bit #6 off). In this case switch bit 11 ON.

    ! Arguments
    integer,     intent(in)               :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)               :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)               :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)               :: KNT                            ! nombre de tovs
    character *9,intent(in)               :: STNID                          ! identificateur du satellite
    logical,     intent(in)               :: RESETQC                        ! yes or not reset QC flag
    integer,     intent(out)              :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(out)              :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(out)              :: rejectionCodArray(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 
    ! Locals
    integer                               :: channelval
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex
    integer                               :: IBIT


    if (.NOT.RESETQC) then
      testIndex = 9
      do nDataIndex=1,KNT
        do nChannelIndex=1,KNO
          if ( KCANO(nChannelIndex,nDataIndex) .NE. 20 ) then
            IBIT = AND(KMARQ(nChannelIndex,nDataIndex), 2**6)
            if ( IBIT .EQ. 0  ) then
              ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**11)
              rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                    rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
              if ( mwbg_debug ) then
                write(*,*)STNID(2:9),' UNCORRECTED TB REJECT.', &
                           'CHANNEL=', KCANO(nChannelIndex,nDataIndex), ' IMARQ= ',KMARQ(nChannelIndex,nDataIndex)
              end if
            end if
          end if
        end do
      end do
    end if

  end subroutine amsuaTest9UncorrectedTbCheck
 
  !--------------------------------------------------------------------------
  !  amsuaTest11RadianceGrossValueCheck
  !--------------------------------------------------------------------------
  subroutine amsuaTest11RadianceGrossValueCheck (KCANO, KNOSAT, KNO, KNT, STNID, RESETQC, PTBO, PMISG, GROSSMIN, &
                                      GROSSMAX, KMARQ, ICHECK, rejectionCodArray)

    !:Purpose:                     11) test 11: Radiance observation "Gross" check (single) 
    !                                               Change this test from full to single. jh nov 2000.

    ! Arguments
    integer,     intent(in)               :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)               :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)               :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)               :: KNT                            ! nombre de tovs
    character *9,intent(in)               :: STNID                          ! identificateur du satellite
    logical,     intent(in)               :: RESETQC                        ! yes or not reset QC flag
    real,        intent(in)               :: PTBO(KNO,KNT)                      ! radiances 
    real,        intent(in)               :: PMISG                          ! missing value 
    real,        intent(in)               :: GROSSMIN(mwbg_maxNumChan)                ! Gross val min 
    real,        intent(in)               :: GROSSMAX(mwbg_maxNumChan)                ! Gross val max 
    integer,     intent(out)              :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(out)              :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(out)              :: rejectionCodArray(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 
    ! Locals
    integer                               :: channelval
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex
    logical                               :: GROSSERROR

    testIndex = 11
    do nDataIndex=1,KNT
      GROSSERROR = .FALSE.
      do nChannelIndex=1,KNO
        if ( KCANO(nChannelIndex,nDataIndex) .NE. 20     .AND. &
            KCANO(nChannelIndex,nDataIndex) .GE.  1     .AND. &
            KCANO(nChannelIndex,nDataIndex) .LE.  mwbg_maxNumChan       ) then  
          if ( PTBO(nChannelIndex,nDataIndex) .NE. PMISG .AND. &
             ( PTBO(nChannelIndex,nDataIndex).LT.GROSSMIN(KCANO(nChannelIndex,nDataIndex)).OR. &
               PTBO(nChannelIndex,nDataIndex).GT.GROSSMAX(KCANO(nChannelIndex,nDataIndex))     ) ) then
            GROSSERROR = .TRUE.
            ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**7)
            rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                   rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
            if ( mwbg_debug ) then
              write(*,*)STNID(2:9),' GROSS CHECK REJECT.', &
                        'CHANNEL=', KCANO(nChannelIndex,nDataIndex), ' TB= ',PTBO(nChannelIndex,nDataIndex)
            end if
          end if
        end if
      end do
    end do

  end subroutine amsuaTest11RadianceGrossValueCheck 
  
  !--------------------------------------------------------------------------
  !  amsuaTest12GrodyClwCheck
  !--------------------------------------------------------------------------
  subroutine amsuaTest12GrodyClwCheck (KCANO, KNOSAT, KNO, KNT, STNID, RESETQC, clw, MISGRODY, MXCLWREJ, &
                                  ICLWREJ, cloudyClwThreshold, KMARQ, ICHECK, rejectionCodArray)

    !:Purpose:                    12) test 12: Grody cloud liquid water check (partial)
    !                                 For Cloud Liquid Water > clwQcThreshold, reject AMSUA-A channels 1-5 and 15.

    ! Arguments
    integer,     intent(in)               :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)               :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)               :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)               :: KNT                            ! nombre de tovs
    character *9,intent(in)               :: STNID                          ! identificateur du satellite
    logical,     intent(in)               :: RESETQC                        ! yes or not reset QC flag
    real,        intent(in)               :: CLW(KNT)                       ! cloud liquid water 
    real,        intent(in)               :: MISGRODY                       ! MISGRODY
    integer,     intent(in)               :: MXCLWREJ                       ! cst 
    integer,     intent(in)               :: ICLWREJ(MXCLWREJ)              !
    real,        intent(in)               :: cloudyClwThreshold             !
    integer,     intent(out)              :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(out)              :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(out)              :: rejectionCodArray(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 
    ! Locals
    integer                               :: channelval
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex
    integer                               :: INDXCAN 

    testIndex = 12
    do nDataIndex=1,KNT
      if ( CLW(nDataIndex) .NE.  MISGRODY  ) then
        if ( CLW(nDataIndex) .GT. mwbg_clwQcThreshold   ) then
          do nChannelIndex=1,KNO
            INDXCAN = ISRCHEQI (ICLWREJ,MXCLWREJ,KCANO(nChannelIndex,nDataIndex))
            if ( INDXCAN.NE.0 )  then
              ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**7)
              rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                       rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
            end if
          end do
          if ( mwbg_debug ) then
            write(*,*)STNID(2:9),'Grody cloud liquid water check', &
                      ' REJECT. CLW= ',CLW(nDataIndex), ' SEUIL= ',mwbg_clwQcThreshold
          end if
        end if

        ! trun on bit=23 for cloud-affected radiances (to be used in gen_bias_corr)
        if ( CLW(nDataIndex) > cloudyClwThreshold ) then
          do nChannelIndex = 1,KNO
            INDXCAN = ISRCHEQI(ICLWREJ,MXCLWREJ,KCANO(nChannelIndex,nDataIndex))
            if ( INDXCAN /= 0 ) KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**23)
          end do
          if ( mwbg_debug ) then
            write(*,*) STNID(2:9),' Grody cloud liquid water check', &
                      ' cloud-affected obs. CLW= ',CLW(nDataIndex), ', threshold= ',cloudyClwThreshold 
          end if
        end if

      end if
    end do

  end subroutine amsuaTest12GrodyClwCheck 

  !--------------------------------------------------------------------------
  !  amsuaTest13GrodyScatteringIndexCheck
  !--------------------------------------------------------------------------
  subroutine amsuaTest13GrodyScatteringIndexCheck (KCANO, KNOSAT, KNO, KNT, STNID, RESETQC, scatw, KTERMER, ITERRAIN, &
                                              MISGRODY, MXSCATREJ, ISCATREJ, KMARQ, ICHECK, rejectionCodArray)

    !:Purpose:                  13) test 13: Grody scattering index check (partial)
    !                               For Scattering Index > 9, reject AMSUA-A channels 1-6 and 15.

    ! Arguments
    integer,     intent(in)                :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)                :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)                :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)                :: KNT                            ! nombre de tovs
    character *9,intent(in)                :: STNID                          ! identificateur du satellite
    logical,     intent(in)                :: RESETQC                        ! yes or not reset QC flag
    real,        intent(in)                :: scatw(KNT)                     ! scattering index 
    integer,     intent(in)                :: KTERMER(KNT)                   ! land sea qualifyer 
    integer,     intent(in)                :: ITERRAIN(KNT)                  ! terrain type 
    real,        intent(in)                :: MISGRODY                       ! MISGRODY
    integer,     intent(in)                :: MXSCATREJ                       ! cst 
    integer,     intent(in)                :: ISCATREJ(MXSCATREJ)              !
    integer,     intent(out)               :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(out)               :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(out)               :: rejectionCodArray(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 
    ! Locals
    integer                                :: channelval
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex
    integer                                :: INDXCAN 
    real                                   :: ZSEUILSCAT

    testIndex = 13
    ZSEUILSCAT = 9.0
    do nDataIndex=1,KNT
      if ( SCATW(nDataIndex) .NE.  MISGRODY  ) then
        if (  KTERMER (nDataIndex) .EQ.  1 .AND. &
             ITERRAIN(nDataIndex) .NE.  0 .AND. &   
             SCATW   (nDataIndex) .GT. ZSEUILSCAT   ) then
          do nChannelIndex=1,KNO
            INDXCAN = ISRCHEQI (ISCATREJ,MXSCATREJ,KCANO(nChannelIndex,nDataIndex))
            if ( INDXCAN.NE.0 )  then
              ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**7)
              rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                       rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
            end if
          end do
          if ( mwbg_debug ) then
            write(*,*)STNID(2:9),'Grody scattering index check', &
                       ' REJECT. SCATW= ',SCATW(nDataIndex), ' SEUIL= ',ZSEUILSCAT
          end if
        end if
      end if
    end do

  end subroutine amsuaTest13GrodyScatteringIndexCheck
 
  !--------------------------------------------------------------------------
  !  amsuaTest14RogueCheck
  !--------------------------------------------------------------------------
  subroutine amsuaTest14RogueCheck (KCANO, KNOSAT, KNO, KNT, STNID, ROGUEFAC, TOVERRST, PTBOMP, &
                                              PMISG, MXSFCREJ, ISFCREJ, KMARQ, ICHECK, rejectionCodArray)

    !:Purpose:                     14) test 14: "Rogue check" for (O-P) Tb residuals out of range.
    !                                  (single/full). Les observations, dont le residu (O-P) 
    !                                  depasse par un facteur (roguefac) l'erreur totale des TOVS.
    !                                  N.B.: a reject by any of the 3 surface channels produces the 
    !                                  rejection of AMSUA-A channels 1-5 and 15.

    ! Arguments
    integer,     intent(in)                :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)                :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)                :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)                :: KNT                            ! nombre de tovs
    character *9,intent(in)                :: STNID                          ! identificateur du satellite
    real,        intent(in)                :: ROGUEFAC(mwbg_maxNumChan)                  ! rogue factor 
    real,        intent(in)                :: TOVERRST(mwbg_maxNumChan,mwbg_maxNumSat)          !  erreur totale TOVs
    real,        intent(in)                :: PTBOMP(KNO,KNT)              ! radiance o-p 
    real,        intent(in)                :: PMISG                          ! missing value
    integer,     intent(in)                :: MXSFCREJ                       ! cst 
    integer,     intent(in)                :: ISFCREJ(MXSFCREJ)              !
    integer,     intent(out)               :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(out)               :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(out)               :: rejectionCodArray(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 
    ! Locals
    integer                                :: channelval
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex
    integer                                :: INDXCAN 
    real                                   :: XCHECKVAL
    logical                                :: SFCREJCT


    testIndex = 14
    do nDataIndex=1,KNT

      SFCREJCT = .FALSE.
      do nChannelIndex=1,KNO
        channelval = KCANO(nChannelIndex,nDataIndex)
        if ( channelval .NE. 20 ) then
          XCHECKVAL = ROGUEFAC(channelval)*TOVERRST(channelval,KNOSAT) 
          if ( PTBOMP(nChannelIndex,nDataIndex)      .NE. PMISG    .AND. &
              ABS(PTBOMP(nChannelIndex,nDataIndex)) .GE. XCHECKVAL     ) then
            ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**16)
            rejectionCodArray(testIndex,channelval,KNOSAT) = &
                rejectionCodArray(testIndex,channelval,KNOSAT) + 1 
            if ( mwbg_debug ) then
              write(*,*)STNID(2:9),'ROGUE CHECK REJECT.NO.', &
                      ' OBS = ',nDataIndex, &
                      ' CHANNEL= ',channelval, &
                      ' CHECK VALUE= ',XCHECKVAL, &
                      ' TBOMP= ',PTBOMP(nChannelIndex,nDataIndex)
            end if
            if ( channelval .EQ. 28 .OR. &
                 channelval .EQ. 29 .OR. &
                 channelval .EQ. 30      ) then
              SFCREJCT = .TRUE.
            end if
          end if
        end if
      end do

      if ( SFCREJCT ) then
        do nChannelIndex=1,KNO
          INDXCAN = ISRCHEQI (ISFCREJ,MXSFCREJ,KCANO(nChannelIndex,nDataIndex))
          if ( INDXCAN .NE. 0 )  then
            if ( ICHECK(nChannelIndex,nDataIndex) .NE. testIndex ) then
               ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
               KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
               KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**16)
               rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                         rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
            end if
          end if
        end do
      end if

    end do
  end subroutine amsuaTest14RogueCheck

  !--------------------------------------------------------------------------
  !  amsuaTest15ChannelSelectionWithIutilst
  !--------------------------------------------------------------------------
  subroutine amsuaTest15ChannelSelectionWithIutilst (KCANO, KNOSAT, KNO, KNT, STNID, KTERMER, ITERRAIN, &
                                              GLINTRP, IUTILST, MXSFCREJ2, ISFCREJ2, KMARQ, ICHECK, rejectionCodArray)

    !:Purpose:                         ! 15) test 15: Channel Selection using array IUTILST(chan,sat)
    !                                    IUTILST = 0 (blacklisted)
    !                                    - 1 (assmilate)
    !                                    - 2 (assimilate over open water only)
    !
    !                                    - We also set QC flag bits 7 and 9 ON for channels with IUTILST=2 
    !                                    over land or sea-ice and we set QC flag bits 7 and 9 ON for channels
    !                                    1-3,15 over land or sea-ice REGARDLESS of IUTILST value 
    !                                    (but IUTILST=0 always for these unassimilated channels).


    ! Arguments
    integer,     intent(in)                :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)                :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)                :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)                :: KNT                            ! nombre de tovs
    character *9,intent(in)                :: STNID                          ! identificateur du satellite
    integer,     intent(in)                :: KTERMER(KNT)                   ! land sea qualifyer 
    integer,     intent(in)                :: ITERRAIN(KNT)                  ! terrain type
    real  ,      intent(in)                :: GLINTRP(KNT)                   ! gl
    integer,     intent(in)                :: IUTILST(mwbg_maxNumChan,mwbg_maxNumSat)          !  channel selection
    integer,     intent(in)                :: MXSFCREJ2                       ! cst 
    integer,     intent(in)                :: ISFCREJ2(MXSFCREJ2)              !
    integer,     intent(out)               :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(out)               :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(out)               :: rejectionCodArray(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 
    ! Locals
    integer                                :: channelval
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex
    integer                                :: ITRN 
    integer                                :: INDXCAN 
    real                                   :: XCHECKVAL
    logical                                :: SFCREJCT

    testIndex = 15

    do nDataIndex=1,KNT
      ITRN = ITERRAIN(nDataIndex)
      if ( KTERMER (nDataIndex) .EQ. 1    .AND. &
           ITERRAIN(nDataIndex) .EQ. -1   .AND. &
           GLINTRP (nDataIndex) .GE. 0.01       ) then
         ITRN = 0
      end if        
      do nChannelIndex=1,KNO
          channelval = KCANO(nChannelIndex,nDataIndex)
          INDXCAN = ISRCHEQI (ISFCREJ2,MXSFCREJ2,channelval)
          if ( INDXCAN .NE. 0 )  then
            if ( KTERMER (nDataIndex) .EQ. 0 .OR. ITRN .EQ. 0 )  then
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**7)
            end if
          end if
          if ( IUTILST(channelval,KNOSAT) .NE. 1 ) then
            SFCREJCT = .FALSE.
            if ( IUTILST(channelval,KNOSAT) .EQ. 0 ) then
              SFCREJCT = .TRUE.
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**11)
            else 
              if ( KTERMER (nDataIndex) .EQ. 0 .OR. ITRN .EQ. 0 )  then
                SFCREJCT = .TRUE.
                KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
                KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**7)
              end if
            end if
            if ( SFCREJCT ) then
              ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
              rejectionCodArray(testIndex,channelval,KNOSAT) = & 
                 rejectionCodArray(testIndex,channelval,KNOSAT) + 1 
              if ( mwbg_debug ) then
                 write(*,*)STNID(2:9),'CHANNEL REJECT: ', &
                        ' OBS = ',nDataIndex, &
                        ' CHANNEL= ',channelval
              end if
            end if
          end if
        end do
    end do

    if ( mwbg_debug ) then
       write(*,*)'ICHECK = ',((ICHECK(nChannelIndex,nDataIndex),nChannelIndex=1,KNO),nDataIndex=1,KNT)
    end if

  end subroutine amsuaTest15ChannelSelectionWithIutilst

  !--------------------------------------------------------------------------
  !  mwbg_allocate1DCharacterArray
  !--------------------------------------------------------------------------
  subroutine mwbg_allocate1DCharacterArray(arrayToAllocate, firstDim)
    !:Purpose: Allocate an array on 1 dimension
    ! Arguments
    character(len=9), intent(out), allocatable :: arrayToAllocate(:)
    integer, intent(in)                 :: firstDim
    
    !locals
    integer                             :: allocStatus

    ! Allocation
    allocStatus = 0
    if (allocated(arrayToAllocate)) deallocate(arrayToAllocate)
    allocate(arrayToAllocate(firstDim), stat = allocStatus)
    if (allocStatus /= 0) then
      call utl_abort('bgckMicrowave_mod: Allocation Error in sub. mwbg_allocate1DCharacterArray')
    end if 
  end subroutine mwbg_allocate1DCharacterArray

  !--------------------------------------------------------------------------
  !  mwbg_allocate1DIntegerArray
  !--------------------------------------------------------------------------
  subroutine mwbg_allocate1DIntegerArray(arrayToAllocate, firstDim)
    !:Purpose: Allocate an array on 1 dimension
    ! Arguments
    integer, intent(inout), allocatable :: arrayToAllocate(:)
    integer, intent(in)                 :: firstDim
    
    !locals
    integer                             :: allocStatus

    ! Allocation
    allocStatus = 0
    if (allocated(arrayToAllocate)) deallocate(arrayToAllocate)
    allocate(arrayToAllocate(firstDim), stat = allocStatus)
    if (allocStatus /= 0) then
      call utl_abort('bgckMicrowave_mod: Allocation Error in sub. mwbg_allocate1DIntegerArray')
    end if 
  end subroutine mwbg_allocate1DIntegerArray

  !--------------------------------------------------------------------------
  !  mwbg_allocate1DLogicalArray
  !--------------------------------------------------------------------------
  subroutine mwbg_allocate1DLogicalArray(arrayToAllocate, firstDim)
    !:Purpose: Allocate an array on 1 dimension
    ! Arguments
    logical, intent(inout), allocatable :: arrayToAllocate(:)
    integer, intent(in)                 :: firstDim
    
    !locals
    integer                             :: allocStatus

    ! Allocation
    allocStatus = 0
    if (allocated(arrayToAllocate)) deallocate(arrayToAllocate)
    allocate(arrayToAllocate(firstDim), stat = allocStatus)
    if (allocStatus /= 0) then
      call utl_abort('bgckMicrowave_mod: Allocation Error in sub. mwbg_allocate1DLogicalArray ')
    end if 
  end subroutine mwbg_allocate1DLogicalArray

  !--------------------------------------------------------------------------
  !  mwbg_allocate2DLogicalArray
  !--------------------------------------------------------------------------
  subroutine mwbg_allocate2DLogicalArray(arrayToAllocate, firstDim, secondDim)
    !:Purpose: Allocate an array on 2 dimension
    ! Arguments
    logical, intent(inout), allocatable :: arrayToAllocate(:,:)
    integer, intent(in)                 :: firstDim
    integer, intent(in)                 :: secondDim
    
    !locals
    integer                             :: allocStatus

    ! Allocation
    allocStatus = 0
    if (allocated(arrayToAllocate)) deallocate(arrayToAllocate)
    allocate(arrayToAllocate(firstDim, secondDim), stat = allocStatus)
    if (allocStatus /= 0) then
      call utl_abort('bgckMicrowave_mod: Allocation Error in sub. mwbg_allocate2DLogicalArray ')
    end if 

  end subroutine mwbg_allocate2DLogicalArray

  !--------------------------------------------------------------------------
  !  allocate3DLogicalArray
  !--------------------------------------------------------------------------
  subroutine allocate3DLogicalArray(arrayToAllocate, firstDim, secondDim, thirdDim)
    !:Purpose: Allocate an array on 3 dimension
    ! Arguments
    logical, intent(inout), allocatable :: arrayToAllocate(:,:,:)
    integer, intent(in)                 :: firstDim
    integer, intent(in)                 :: secondDim
    integer, intent(in)                 :: thirdDim
    
    !locals
    integer                             :: allocStatus

    ! Allocation
    allocStatus = 0
    if (allocated(arrayToAllocate)) deallocate(arrayToAllocate)
    allocate(arrayToAllocate(firstDim, secondDim, thirdDim), stat = allocStatus)
    if (allocStatus /= 0) then
      call utl_abort('bgckMicrowave_mod: Allocation Error in sub. allocate3DLogicalArray ')
    end if 

  end subroutine allocate3DLogicalArray

  !--------------------------------------------------------------------------
  !  mwbg_allocate2DIntegerArray
  !--------------------------------------------------------------------------
  subroutine mwbg_allocate2DIntegerArray(arrayToAllocate, firstDim, secondDim)
    !:Purpose: Allocate an array on 2 dimension
    ! Arguments
    integer, intent(inout), allocatable :: arrayToAllocate(:,:)
    integer, intent(in)                 :: firstDim
    integer, intent(in)                 :: secondDim
    
    !locals
    integer                             :: allocStatus

    ! Allocation
    allocStatus = 0
    if (allocated(arrayToAllocate)) deallocate(arrayToAllocate)
    allocate(arrayToAllocate(firstDim, secondDim), stat = allocStatus)
    if (allocStatus /= 0) then
      call utl_abort('bgckMicrowave_mod: Allocation Error in sub. mwbg_allocate2DIntegerArray ')
    end if 

  end subroutine mwbg_allocate2DIntegerArray

  !--------------------------------------------------------------------------
  !  mwbg_allocate3DIntegerArray
  !--------------------------------------------------------------------------
  subroutine mwbg_allocate3DIntegerArray(arrayToAllocate, firstDim, secondDim, thirdDim)
    !:Purpose: Allocate an array on 3 dimension
    ! Arguments
    integer, intent(inout), allocatable :: arrayToAllocate(:,:,:)
    integer, intent(in)                 :: firstDim
    integer, intent(in)                 :: secondDim
    integer, intent(in)                 :: thirdDim
    
    !locals
    integer                             :: allocStatus

    ! Allocation
    allocStatus = 0
    if (allocated(arrayToAllocate)) deallocate(arrayToAllocate)
    allocate(arrayToAllocate(firstDim, secondDim, thirdDim), stat = allocStatus)
    if (allocStatus /= 0) then
      call utl_abort('bgckMicrowave_mod: Allocation Error in sub. mwbg_allocate3DIntegerArray ')
    end if 

  end subroutine mwbg_allocate3DIntegerArray

  !--------------------------------------------------------------------------
  !  mwbg_allocate1DRealArray
  !--------------------------------------------------------------------------
  subroutine mwbg_allocate1DRealArray(arrayToAllocate, firstDim)
    !:Purpose: Allocate an Real array on 1 dimension
    ! Arguments
    real,    intent(inout), allocatable :: arrayToAllocate(:)
    integer, intent(in)                 :: firstDim
    
    !locals
    integer                             :: allocStatus

    ! Allocation
    allocStatus = 0
    if (allocated(arrayToAllocate)) deallocate(arrayToAllocate)
    allocate(arrayToAllocate(firstDim), stat = allocStatus)
    if (allocStatus /= 0) then
      call utl_abort('bgckMicrowave_mod: Allocation Error in sub. mwbg_allocate1DRealArray ')
    end if 
  end subroutine mwbg_allocate1DRealArray

  !--------------------------------------------------------------------------
  !  mwbg_allocate2DRealArray
  !--------------------------------------------------------------------------
  subroutine mwbg_allocate2DRealArray(arrayToAllocate, firstDim, secondDim)
    !:Purpose: Allocate an Real array on 2 dimension
    ! Arguments
    real,    intent(inout), allocatable :: arrayToAllocate(:,:)
    integer, intent(in)                 :: firstDim
    integer, intent(in)                 :: secondDim
    
    !locals
    integer                             :: allocStatus

    ! Allocation
    allocStatus = 0
    if (allocated(arrayToAllocate)) deallocate(arrayToAllocate)
    allocate(arrayToAllocate(firstDim, secondDim), stat = allocStatus)
    if (allocStatus /= 0) then
      call utl_abort('bgckMicrowave_mod: Allocation Error in sub. mwbg_allocate2DRealArray ')
    end if 

  end subroutine mwbg_allocate2DRealArray

  !--------------------------------------------------------------------------
  !  allocate3DRealArray
  !--------------------------------------------------------------------------
  subroutine allocate3DRealArray(arrayToAllocate, firstDim, secondDim, thirdDim)
    !:Purpose: Allocate an array on 3 dimension
    ! Arguments
    real,    intent(inout), allocatable :: arrayToAllocate(:,:,:)
    integer, intent(in)                 :: firstDim
    integer, intent(in)                 :: secondDim
    integer, intent(in)                 :: thirdDim
    
    !locals
    integer                             :: allocStatus

    ! Allocation
    allocStatus = 0
    if (allocated(arrayToAllocate)) deallocate(arrayToAllocate)
    allocate(arrayToAllocate(firstDim, secondDim, thirdDim), stat = allocStatus)
    if (allocStatus /= 0) then
      call utl_abort('bgckMicrowave_mod: Allocation Error in sub. allocate3DRealArray ')
    end if 

  end subroutine allocate3DRealArray

  !--------------------------------------------------------------------------
  !  copy1Dimto2DimRealArray
  !--------------------------------------------------------------------------
  subroutine copy1Dimto2DimRealArray(oneDimArray, firstDim, secondDim, twoDimArray)
    !:Purpose: copy 1 dim Real Array into 2D real array given dim1 and dim2 
    ! Arguments
    integer, intent(in)                 :: firstDim
    integer, intent(in)                 :: secondDim
    real,    intent(in)                 :: oneDimArray(firstDim*secondDim)
    real,    intent(inout)              :: twoDimArray(firstDim, secondDim)
    
    !locals
    integer                             :: firstDimIndex 
    integer                             :: secondDimIndex 
    integer                             :: productDim 
    integer                             :: productDimIndex 

    ! copy the original input 1D array to 2D array. The 2D arrays are used in this s/r.
    do secondDimIndex=1,secondDim
      do firstDimIndex=1,firstDim
        productDimIndex = (secondDimIndex-1)*firstDim + firstDimIndex 
        twoDimArray(firstDimIndex,secondDimIndex) = oneDimArray(productDimIndex)
      end do
    end do

  end subroutine copy1Dimto2DimRealArray

  !--------------------------------------------------------------------------
  !  copy1Dimto2DimIntegerArray
  !--------------------------------------------------------------------------
  subroutine copy1Dimto2DimIntegerArray(oneDimArray, firstDim, secondDim, twoDimArray)
    !:Purpose: copy 1 dim Integer Array into 2D Integer array given dim1 and dim2 
    ! Arguments
    integer, intent(in)                 :: firstDim
    integer, intent(in)                 :: secondDim
    integer, intent(in)                 :: oneDimArray(firstDim*secondDim)
    integer, intent(inout)              :: twoDimArray(firstDim, secondDim)
    
    !locals
    integer                             :: firstDimIndex 
    integer                             :: secondDimIndex 
    integer                             :: productDim 
    integer                             :: productDimIndex 

    ! copy the original input 1D array to 2D array. The 2D arrays are used in this s/r.
    do secondDimIndex=1,secondDim
      do firstDimIndex=1,firstDim
        productDimIndex = (secondDimIndex-1)*firstDim + firstDimIndex 
        twoDimArray(firstDimIndex,secondDimIndex) = oneDimArray(productDimIndex)
      end do
    end do

  end subroutine copy1Dimto2DimIntegerArray

  !--------------------------------------------------------------------------
  !  mwbg_tovCheckAmsua
  !--------------------------------------------------------------------------
  subroutine mwbg_tovCheckAmsua(TOVERRST, IUTILST, KSAT,  KTERMER, KORBIT, ICANO, ICANOMP, ZO, ZCOR, &
                                ZOMP, ICHECK, KNO, KNT, PMISG, KNOSAT, KCHKPRF, &
                                ISCNPOS, MGINTRP, MTINTRP, GLINTRP, ITERRAIN, SATZEN, &
                                IMARQ, ident, clw_avg, scatw, rejectionCodArray, STNID, RESETQC, ZLAT)

  
    !:Purpose:          Effectuer le controle de qualite des radiances tovs.
    !
    !NOTES  
    !               Quinze tests sont effectues menant aux erreurs suivantes:
    !                  - 1) topography reject,
    !                  - 2) invalid land/sea qualifier,
    !                  - 3) invalid terrain type,
    !                  - 4) invalid field of view number,
    !                  - 5) satellite zenith angle out of range,
    !                  - 6) inconsistent field of view and sat. zenith angle,
    !                  - 7) inconsistent land/sea qualifier and model mask,
    !                  - 8) inconsistent terrain type and model ice, (NOT USED)
    !                  - 9) uncorrected radiance,
    !                  - 10) rejected by RTTOV,
    !                  - 11) radiance gross check failure,
    !                  - 12) cloud liquid water reject,
    !                  - 13) scattering index reject,
    !                  - 14) radiance residual rogue check failure,
    !                  - 15) channel reject (channel selection).
    !                  - **) set terrain type to sea ice given certain conditions
    implicit none 
    !Arguments:
    integer, intent(in)                    :: IUTILST(mwbg_maxNumChan,mwbg_maxNumSat) !channel Selection using array IUTILST(chan,sat)
    !                                                               IUTILST = 0 (blacklisted)
    !                                                               1 (assmilate)
    !                                                               2 (assimilate over open water only)
    real, intent(in)                       :: TOVERRST(mwbg_maxNumChan,mwbg_maxNumSat)! l'erreur totale des TOVS
    integer, intent(in)                    :: KSAT(:)            ! numero d'identificateur du satellite
    integer, intent(in)                    :: KTERMER(:)         ! indicateur terre/mer
    integer, intent(in)                    :: ISCNPOS(:)         ! position sur le "scan"
    integer, intent(in)                    :: KORBIT(:)          ! numero d'orbite
    integer, intent(in)                    :: ICANO(:)       ! canaux des observations
    integer, intent(inout)                 :: ITERRAIN(:)        ! indicateur du type de terrain
    integer, intent(in)                    :: ICANOMP(:)     ! canaux des residus (o-p)
    integer, intent(in)                    :: KNO                  ! nombre de canaux des observations 
    integer, intent(in)                    :: KNT                  ! nombre de tovs
    integer, intent(in)                    :: KNOSAT               ! numero de satellite (i.e. indice)
    integer, intent(inout)                 :: IMARQ(:)       ! marqueurs des radiances
    real, intent(in)                       :: ZO(:)          ! radiances
    real, intent(in)                       :: ZCOR(:)        ! correction aux radiances
    real, intent(in)                       :: ZOMP(:)        ! residus (o-p)
    real, intent(in)                       :: MGINTRP(:)         ! masque terre/mer du modele
    real, intent(in)                       :: MTINTRP(:)         ! topographie du modele
    real, intent(in)                       :: GLINTRP(:)         ! etendue de glace du modele
    real, intent(in)                       :: SATZEN(:)          ! angle zenith du satellite (deg.)
    real, intent(in)                       :: ZLAT(:)            ! latitude
    real, intent(in)                       :: PMISG                ! missing value
    character *9, intent(in)               :: STNID                ! identificateur du satellite
    logical, intent(in)                    :: RESETQC              ! reset du controle de qualite?
    integer,allocatable, intent(out)       :: ICHECK(:,:)          ! indicateur controle de qualite tovs par canal 
    !                                                                 =0, ok,
    !                                                                 >0, rejet,
    real,allocatable,  intent(out)         :: clw_avg(:)            ! Averaged (Background + obs)retrieved cloud liquid water, 
    !                                                                 from observation and background
    real,allocatable, intent(out)          :: scatw(:)              ! scattering index over water

    integer,allocatable, intent(out)       :: KCHKPRF(:)            ! indicateur global controle de qualite tovs. Code:
    !                                                                 =0, ok,
    !                                                                 >0, rejet d'au moins un canal
    integer,allocatable, intent(out)       :: ident(:)              !ATMS Information flag (ident) values (new BURP element 025174 in header)
    !                                                               FOR AMSUA just fill with zeros

    integer, intent(inout)                   :: rejectionCodArray(mwbg_maxNumTest,mwbg_maxNumChan,mwbg_maxNumSat)  ! cumul du nombre de rejet 
    !locals
    integer, parameter                     :: mwbg_maxScanAngleHIRS= 56 
    integer, parameter                     :: mwbg_maxScanAngleAMSU= 30 
    integer, parameter                     :: MXCLWREJ  =  6 
    integer, parameter                     :: MXSFCREJ  =  6 
    integer, parameter                     :: MXSFCREJ2 =  4 
    integer, parameter                     :: MXSCATREJ =  7 
    integer, parameter                     :: MXCANPRED =  9 
    integer, parameter                     :: JPMXSFC = 2
    real, parameter                        :: cloudyClwThreshold = 0.3
    
    integer                                :: KMARQ   (KNO,KNT)
    integer                                :: KCANO   (KNO,KNT)
    integer                                :: KCANOMP (KNO,KNT)
    real                                   :: PTBO    (KNO,KNT)
    real                                   :: PTBCOR  (KNO,KNT)
    real                                   :: PTBOMP  (KNO,KNT)
    real,allocatable                       :: clw(:)                ! obs retrieved cloud liquid water
    integer                                :: MAXVAL
    integer                                :: JI
    integer                                :: JJ
    integer                                :: INDX8
    integer                                :: INDX12
    integer                                :: INO
    integer                                :: ICHN
    integer                                :: JK
    integer                                :: IBIT
    integer                                :: JC
    integer                                :: INDX
    integer                                :: INDXCAN
    integer                                :: ITRN
    integer                                :: alloc_status 
    integer                                :: ICLWREJ (MXCLWREJ)
    integer                                :: ISFCREJ (MXSFCREJ)
    integer                                :: ISFCREJ2(MXSFCREJ2)
    integer                                :: ISCATREJ(MXSCATREJ)
    real                                   :: EPSILON
    real                                   :: ZANGL
    real                                   :: MISGRODY
    real                                   :: ZSEUILSCAT
    real                                   :: APPROXIM
    real                                   :: ANGDif
    real                                   :: XCHECKVAL
    real                                   :: GROSSMIN(mwbg_maxNumChan)
    real                                   :: GROSSMAX(mwbg_maxNumChan) 
    real                                   :: ROGUEFAC(mwbg_maxNumChan)
    real                                   :: tb23 (KNT)
    real                                   :: tb31 (KNT)
    real                                   :: tb50 (KNT)
    real                                   :: tb53 (KNT)
    real                                   :: tb89 (KNT)
    real                                   :: tb23_P (KNT)
    real                                   :: tb31_P (KNT)
    real                                   :: tb50_P (KNT)
    real                                   :: tb53_P (KNT)
    real                                   :: tb89_P (KNT)
    real                                   :: ice  (KNT)
    real                                   :: tpw  (KNT)
    real                                   :: scatl(KNT)
    integer                                :: err (KNT)
    integer                                :: rain(KNT)
    integer                                :: snow(KNT)
    logical                                :: GROSSERROR
    logical                                :: FULLREJCT
    logical                                :: SFCREJCT
    logical, save                          :: LLFIRST


    data  LLFIRST / .TRUE. /
    data  EPSILON / 0.01   /
    data  MISGRODY / -99.     /
    !      data  ROGUEFAC/ 3.0    / changed, jh, from 3 to 4, jan 2001
    !      data  ROGUEFAC/ 4.0    /
    data  ROGUEFAC / 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
                     4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
                     4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 2.0, 2.0, 2.0, &
                     3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
                     4.0, 2.0/
    data  ICLWREJ  / 28, 29, 30, 31, 32, 42 /
    data  ISFCREJ  / 28, 29, 30, 31, 32, 42 /
    data  ISCATREJ / 28, 29, 30, 31, 32, 33, 42 /
    data  ISFCREJ2 / 28, 29, 30, 42 /
                   
    data GROSSMIN / 200., 190., 190., 180., 180., 180., 170., &
                    170., 180., 170., 170., 170., 180., 180., &
                    180., 180., 170., 180., 180., 000., 120., &
                    190., 180., 180., 180., 190., 200., 120., &
                    120., 160., 190., 190., 200., 190., 180., &
                    180., 180., 180., 190., 190., 200., 130./

    data GROSSMAX / 270., 250., 250., 250., 260., 280., 290., &
                    320., 300., 320., 300., 280., 320., 300., &
                    290., 280., 330., 350., 350., 000., 310., &
                    300., 250., 250., 270., 280., 290., 310., &
                    310., 310., 300., 300., 260., 250., 250., &
                    250., 260., 260., 270., 280., 290., 330./  

    ! Allocation
    call mwbg_allocate1DRealArray(clw_avg, KNT)
    call mwbg_allocate1DRealArray(clw,     KNT)
    call mwbg_allocate1DRealArray(scatw,   KNT)

    call mwbg_allocate1DIntegerArray(kchkprf, KNT)
    call mwbg_allocate1DIntegerArray(ident, KNT)
    call mwbg_allocate2DIntegerArray(icheck, KNO, KNT)

    ! copy the original input 1D array to 2D array. The 2D arrays are used in this s/r.
    call copy1Dimto2DimIntegerArray(ICANO, KNO, KNT, KCANO)
    call copy1Dimto2DimIntegerArray(IMARQ, KNO, KNT, KMARQ)
    call copy1Dimto2DimRealArray(ZCOR, KNO, KNT, PTBCOR)
    call copy1Dimto2DimRealArray(ZO, KNO, KNT, PTBO)
    ! small check
    do INDX = 1, KNO*KNT
      if ( ICANO(INDX) /= ICANOMP(INDX) ) then
        call utl_abort('bgckMicrowave_mod: ERROR IN DIMENSIONS OF TOVS data')
      end if
    end do
    call copy1Dimto2DimIntegerArray(ICANOMP, KNO, KNT, KCANOMP)
    call copy1Dimto2DimRealArray(ZOMP, KNO, KNT, PTBOMP)

    ! Initialisation, la premiere fois seulement!
    if (LLFIRST) then
       rejectionCodArray(:,:,:) = 0
       LLFIRST = .FALSE.
    end if
    ! fill ident with zeros ONLY for consistency with ATMS
    ident(:) = 0
    ICHECK(:,:) = 0
     if ( RESETQC ) KMARQ(:,:) = 0

    do JJ=1,KNT
      do JI=1,KNO
        if ( KCANO(JI,JJ) .NE. KCANOMP(JI,JJ) ) then
          call utl_abort('bgckMicrowave_mod: INCONSISTENT CHANNEL LISTS FOR TOVS data')
        end if
      end do
    end do

    !     Grody parameters are   extract required channels:
    call extractParamForGrodyRun (KCANO, ptbo, ptbomp, ptbcor, KNT, KNO, pmisg, &
                                     tb23,   tb31,   tb50,   tb53,   tb89, &
                                     tb23_P, tb31_P, tb50_P, tb53_P, tb89_P)
    
    !  Run Grody AMSU-A algorithms.

    call grody (err, knt, tb23, tb31, tb50, tb53, tb89, &
                tb23_P, tb31_P, tb50_P, tb53_P, tb89_P, &
                satzen, zlat, ktermer, ice, tpw, clw, clw_avg, &
                rain, snow, scatl, scatw)   

    ! 10) test 10: RTTOV reject check (single)
    ! Rejected datum flag has bit #9 on.
    call amsuaTest10RttovRejectCheck (KCANO, KNOSAT, KNO, KNT, RESETQC, STNID, KMARQ, ICHECK, rejectionCodArray)

    ! 1) test 1: Topography check (partial)
    ! Channel 6 is rejected for topography >  250m.
    ! Channel 7 is rejected for topography > 2000m.
    call amsuaTest1TopographyCheck (KCANO, KNOSAT, KNO, KNT, STNID, MTINTRP, KMARQ, ICHECK, rejectionCodArray)
 
    ! 2) test 2: "Land/sea qualifier" code check (full)
    ! allowed values are: 0, land,
    !                       1, sea,
    !                       2, coast.
    call amsuaTest2LandSeaQualifierCheck (KCANO, KNOSAT, KNO, KNT, STNID, KTERMER, KMARQ, ICHECK, rejectionCodArray)

    ! 3) test 3: "Terrain type" code check (full)
    !   allowed values are: -1, missing,
    !                        0, sea-ice,
    !                        1, snow on land.
    call amsuaTest3TerrainTypeCheck (KCANO, KNOSAT, KNO, KNT, STNID, ITERRAIN, mwbg_intMisg, KMARQ, ICHECK, rejectionCodArray)
 
    ! 4) test 4: Field of view number check (full)
    !
    ! Field of view acceptable range is [1,mwbg_maxScanAngleAMSU]  for AMSU footprints.
    call amsuaTest4FieldOfViewCheck (KCANO, KNOSAT, KNO, KNT, STNID, ISCNPOS, mwbg_maxScanAngleAMSU, KMARQ, ICHECK, rejectionCodArray)

    ! 5) test 5: Satellite zenith angle check (full)
    ! Satellite zenith angle acceptable range is [0.,60.].
    call amsuaTest5ZenithAngleCheck (KCANO, KNOSAT, KNO, KNT, STNID, SATZEN, PMISG, KMARQ, ICHECK, rejectionCodArray)
    ! 6) test 6: "Sat. zenith angle"/"field of view" consistency check.  (full)
    ! Acceptable difference between "Satellite zenith angle"  and
    ! "approximate angle computed from field of view number" is 1.8 degrees.
    call amsuaTest6ZenAngleAndFovConsistencyCheck (KCANO, KNOSAT, KNO, KNT, STNID, SATZEN,  PMISG, &
                                                  ISCNPOS, mwbg_intMisg, mwbg_maxScanAngleAMSU, KMARQ, ICHECK, rejectionCodArray) 
    ! 7) test 7: "Land/sea qual."/"model land/sea" consistency check.    (full)
    ! Acceptable conditions are:
    !       a) both over ocean (ktermer=1; mg<0.01), new threshold 0.20, jh dec 2000,
    !       b) both over land  (ktermer=0; mg>0.80), new threshold 0.50, jh dec 2000.
    ! Other conditions are unacceptable.
    call amsuaTest7landSeaQualifyerAndModelLandSeaConsistencyCheck (KCANO, KNOSAT, KNO, KNT, STNID, MGINTRP, KTERMER, &
                                                               KMARQ, ICHECK, rejectionCodArray)

    ! 8) test 8: "Terrain type"/"Land/sea qual."/"model ice" consistency check. (full)
    ! Unacceptable conditions are:
    !        a) terrain is sea-ice and model has no ice(iterrain=0; gl<0.01).
    !        b) terrain is sea-ice and land/sea qualifier is land (iterrain=0; ktermer=0).
    !        c) terrain is snow on land and land/sea qualifier is sea (iterrain=1; ktermer=1).
    !        d) terrain is missing, land/sea qualifier is sea and model has ice(iterrain=-1; ktermer=1; gl>0.01). (enleve jh, jan 2001)
    ! NOT doNE ANYMORE 
    
    ! 9) test 9: Uncorrected Tb check (single)
    ! Uncorrected datum (flag bit #6 off). In this case switch bit 11 ON.
    call amsuaTest9UncorrectedTbCheck (KCANO, KNOSAT, KNO, KNT, STNID, RESETQC, KMARQ, ICHECK, rejectionCodArray) 
    ! 11) test 11: Radiance observation "Gross" check (single) 
    !  Change this test from full to single. jh nov 2000.
    call amsuaTest11RadianceGrossValueCheck (KCANO, KNOSAT, KNO, KNT, STNID, RESETQC, PTBO, PMISG, GROSSMIN, &
                                      GROSSMAX, KMARQ, ICHECK, rejectionCodArray)
    ! 12) test 12: Grody cloud liquid water check (partial)
    ! For Cloud Liquid Water > clwQcThreshold, reject AMSUA-A channels 1-5 and 15.
    call amsuaTest12GrodyClwCheck (KCANO, KNOSAT, KNO, KNT, STNID, RESETQC, clw, MISGRODY, MXCLWREJ, &
                                  ICLWREJ, cloudyClwThreshold, KMARQ, ICHECK, rejectionCodArray)
    ! 13) test 13: Grody scattering index check (partial)
    ! For Scattering Index > 9, reject AMSUA-A channels 1-6 and 15.
    call amsuaTest13GrodyScatteringIndexCheck (KCANO, KNOSAT, KNO, KNT, STNID, RESETQC, scatw, KTERMER, ITERRAIN, &
                                              MISGRODY, MXSCATREJ, ISCATREJ, KMARQ, ICHECK, rejectionCodArray)

    ! 14) test 14: "Rogue check" for (O-P) Tb residuals out of range. (single/full)
    ! Les observations, dont le residu (O-P) depasse par un facteur (roguefac) l'erreur totale des TOVS.
    ! N.B.: a reject by any of the 3 surface channels produces the rejection of AMSUA-A channels 1-5 and 15. 
    call amsuaTest14RogueCheck (KCANO, KNOSAT, KNO, KNT, STNID, ROGUEFAC, TOVERRST, PTBOMP, &
                                              PMISG, MXSFCREJ, ISFCREJ, KMARQ, ICHECK, rejectionCodArray)

    ! 15) test 15: Channel Selection using array IUTILST(chan,sat)
    !  IUTILST = 0 (blacklisted)
    !            1 (assmilate)
    !            2 (assimilate over open water only)
    !
    !  We also set QC flag bits 7 and 9 ON for channels with IUTILST=2 
    !  over land or sea-ice
    !    and 
    !  we set QC flag bits 7 and 9 ON for channels 1-3,15 over land
    !  or sea-ice REGARDLESS of IUTILST value (but IUTILST=0 always for
    !  these unassimilated channels).
    call amsuaTest15ChannelSelectionWithIutilst (KCANO, KNOSAT, KNO, KNT, STNID, KTERMER, ITERRAIN, GLINTRP, IUTILST, &
                                              MXSFCREJ2, ISFCREJ2, KMARQ, ICHECK, rejectionCodArray)
    !  Synthese de la controle de qualite au niveau de chaque point
    !  d'observation. Code:
    !            =0, aucun rejet,
    !            >0, au moins un canal rejete.

    do JJ=1,KNT
      KCHKPRF(JJ) = 0
      do JI=1,KNO 
        KCHKPRF(JJ) = MAX(KCHKPRF(JJ),ICHECK(JI,JJ))
      end do
    end do

    if ( mwbg_debug ) then
      write(*,*)'KCHKPRF = ',(KCHKPRF(JJ),JJ=1,KNT)
    end if 

    ! Copy the modified FLAG to the 1D array, used outside this s/r. 
    do JJ=1,KNT
      do JI=1,KNO
        INDX = (JJ-1)*KNO + JI 
        IMARQ(INDX) = KMARQ(JI,JJ)
      end do
    end do

    !###############################################################################
    ! FINAL STEP: set terrain type to sea ice given certain conditions
    !###############################################################################
    write(*,*) ' ==> setTerrainTypeToSeaIce : '
    call setTerrainTypeToSeaIce(GLINTRP, KTERMER, ITERRAIN)
    
  end subroutine mwbg_tovCheckAmsua

  !--------------------------------------------------------------------------
  !  mwbg_qcStats
  !--------------------------------------------------------------------------
  subroutine mwbg_qcStats(instName, ICHECK, ICAN, KNOSAT, &
                              KNO, KNT, satelliteId, LDprint, rejectionCodArray, rejectionCodArray2)

    !:Purpose:          Cumuler ou imprimer des statistiques decriptives
    !                   des rejets tovs.
    implicit none 
    !Arguments:
    character(*), intent(in)               :: instName                           ! Instrument Name
    integer, intent(in)                    :: rejectionCodArray(mwbg_maxNumTest,mwbg_maxNumChan,mwbg_maxNumSat)       ! Nombre d'element rejetes par type de test
    integer, intent(in), optional          :: rejectionCodArray2(mwbg_maxNumTest,mwbg_maxNumChan,mwbg_maxNumSat)      ! Nombre d'element rejetes par type de test
    integer, intent(in)                    :: ICHECK(:,:)                        ! indicateur controle de qualite tovs par canal 
    !                                                                              =0, ok,
    !                                                                              >0, rejet,
    integer, intent(in)                    :: ICAN(KNO*KNT)                      ! canaux des observations
    integer, intent(in)                    :: KNOSAT                             ! numero d'identificateur du satellite
    integer, intent(in)                    :: KNO                                ! nombre de canaux des observations 
    integer, intent(in)                    :: KNT                                ! nombre de tovs
    character *9, intent(in)               :: satelliteId(:)                     ! identificateur du satellite
    logical, intent(in)                    :: LDprint                            ! mode: imprimer ou cumuler?
    !Locals
    integer                                :: numSats
    integer                                :: JI
    integer                                :: JJ
    integer                                :: JK
    integer                                :: INTOTOBS
    integer                                :: INTOTACC
    integer, save                          :: INTOT(mwbg_maxNumSat)
    integer, save                          :: INTOTRJF(mwbg_maxNumSat)
    integer, save                          :: INTOTRJP(mwbg_maxNumSat)
    integer                                :: KCANO(KNO,KNT)                      ! canaux des observations


    logical, save                          :: LLFIRST = .True.
    logical                                :: FULLREJCT
    logical                                :: FULLACCPT

    ! Initialize
    if ( LLFIRST ) then
      do JJ = 1, mwbg_maxNumSat
        INTOTRJF(JJ) = 0
        INTOTRJP(JJ) = 0
        INTOT(JJ)  = 0
      end do
      LLFIRST = .false.
    end if

    if (.NOT. LDprint ) then

      ! copy the original input 1D array to 2D array. The 2D arrays are used in this s/r.
      call copy1Dimto2DimIntegerArray(ICAN, KNO, KNT, KCANO)
      ! Accumulate statistics on rejects
      do JJ = 1, KNT

        INTOT(KNOSAT) = INTOT(KNOSAT) + 1
        ! Fully accepted, fully rejected or partially rejected?
        FULLREJCT = .true.
        FULLACCPT = .true.
        if (instName == "AMSUA") then 
          do JI = 1, KNO
            if ( KCANO(JI,JJ) .NE. 20 ) then
              if ( ICHECK(JI,JJ) .NE. 0 ) then
                FULLACCPT = .false.
              else
                FULLREJCT = .false.
              end if
            end if
          end do
          if ( FULLREJCT ) then
            INTOTRJF(KNOSAT) = INTOTRJF(KNOSAT) + 1
          end if
          if ( .NOT.FULLREJCT .AND. .NOT.FULLACCPT ) then
            INTOTRJP(KNOSAT) = INTOTRJP(KNOSAT) + 1
          end if
        else if  (instName == "ATMS") then 
          do JI = 1, KNO
            if ( ICHECK(JI,JJ) .NE. 0 ) then
              FULLACCPT = .false.
            else
              FULLREJCT = .false.
            end if
          end do
          if ( FULLREJCT ) then
            INTOTRJF(KNOSAT) = INTOTRJF(KNOSAT) + 1
          end if
          if ( .NOT.FULLREJCT .AND. .NOT.FULLACCPT ) then
            INTOTRJP(KNOSAT) = INTOTRJP(KNOSAT) + 1
          end if
        end if
      end do
    else

      numSats = size(satelliteId)
      ! Print statistics
      do JK = 1, numSats

        INTOTOBS = INTOT(JK)
        INTOTACC = INTOTOBS - INTOTRJF(JK) - INTOTRJP(JK)
          write(*,'(/////50("*"))')
          write(*,'(     50("*")/)')
          write(*,'(T5,"SUMMARY OF QUALITY CONTROL FOR ", &
           A8)') satelliteId(JK) 
          write(*,'(T5,"------------------------------------- ",/)')
          write(*,'( &
           " - TOTAL NUMBER OF OBS.    = ",I10,/ &
           " - TOTAL FULL REJECTS      = ",I10,/ &
           " - TOTAL PARTIAL REJECTS   = ",I10,/ &
           "   ------------------------------------",/ &
           "   TOTAL FULLY ACCEPTED    = ",I10,/)') &
            INTOTOBS, INTOTRJF(JK), INTOTRJP(JK), INTOTACC

        if (instName == "AMSUA") then         
          write(*,'(//,1x,114("-"))')
           write(*,'(t10,"|",t47,"REJECTION CATEGORIES")')
          write(*,'(" CHANNEL",t10,"|",105("-"))')
          write(*,'(t10,"|",15i7)') (JI,JI=1,mwbg_maxNumTest)
          write(*,'(1x,"--------|",105("-"))')
          do JJ = 1, mwbg_maxNumChan
             write(*,'(3X,I2,t10,"|",15I7)') JJ,(rejectionCodArray(JI,JJ,JK), &
                                      JI=1,mwbg_maxNumTest)
          end do
          print *, ' '
          print *, ' '
          print *, ' -----------------------------------------------------'
          print *, ' Definition of rejection categories:'
          print *, ' -----------------------------------------------------'
          print *, '  1 - topography reject'
          print *, '  2 - invalid land/sea qualifier'
          print *, '  3 - invalid terrain type'
          print *, '  4 - invalid field of view number'
          print *, '  5 - satellite zenith angle out of range '
          print *, '  6 - inconsistent field of view and sat. zenith angle'
          print *, '  7 - inconsistent land/sea qualifier and model mask'
          print *, '  8 - inconsistent terrain type and land/sea', &
                   ' qualifier/model ice (NOT doNE)'
          print *, '  9 - uncorrected radiance'
          print *, ' 10 - rejected by RTTOV'
          print *, ' 11 - radiance gross check failure'
          print *, ' 12 - cloud liquid water reject'
          print *, ' 13 - scattering index reject'
          print *, ' 14 - radiance residual rogue check failure'
          print *, ' 15 - rejection by channel selection'
          print *, ' -----------------------------------------------------'
          print *, ' ' 

        else if (instName == "ATMS") then 
          write(*,'(//,1x,59("-"))')
          write(*,'(t10,"|",t19,"1. REJECTION CATEGORIES")')
          write(*,'(" CHANNEL",t10,"|",50("-"))')
          write(*,'(t10,"|",5i7)') (JI,JI=1,mwbg_maxNumTest)
          write(*,'(1x,"--------|",50("-"))')
          do JJ = 1, mwbg_maxNumChan
            write(*,'(3X,I2,t10,"|",5I7)') JJ,(rejectionCodArray(JI,JJ,JK), &
                                        JI=1,mwbg_maxNumTest)
          end do
          write(*,'(1x,59("-"))')
          write(*,'(//,1x,59("-"))')
          write(*,'(t10,"|",t19,"2. QC2 REJECT CATEGORIES")')
          write(*,'(" CHANNEL",t10,"|",50("-"))') 
          write(*,'(t10,"|",5i7)') (JI,JI=1,mwbg_maxNumTest)
          write(*,'(1x,"--------|",50("-"))')
          do JJ = 1, mwbg_maxNumChan
            write(*,'(3X,I2,t10,"|",5I7)') JJ,(rejectionCodArray2(JI,JJ,JK), &
                                        JI=1,mwbg_maxNumTest)          
          end do
          print *, ' '
          print *, ' '
          print *, ' -----------------------------------------------------'
          print *, ' Definition of rejection categories: '
          print *, ' -----------------------------------------------------'
          print *, '  1 - first bgckAtms program reject [bit 7]'
          print *, '  2 - topography reject'
          print *, '  3 - uncorrected radiance'
          print *, '  4 - innovation (O-P) based reject'
          print *, '  5 - rejection by channel selection'
          print *, ' -----------------------------------------------------'
          print *, ' '
          print *, ' QC2 REJECT numbers in Table 2 are for data that '
          print *, ' passed test 1 (data with QC flag bit 7 OFF)'
          print *, ' '
        end if 
      end do
    end if
    print *,'END OF mwbd_qcStats'
  end subroutine mwbg_qcStats

  !--------------------------------------------------------------------------
  !  resetQcC
  !--------------------------------------------------------------------------
  subroutine resetQcCases(RESETQC, KCHKPRF, globMarq)
    !:Purpose:        allumer la bit (6) indiquant que l'observation a un element
    !                 rejete par le controle de qualite de l'AO.
    !                 N.B.: si on est en mode resetqc, on remet le marqueur global a
    !                 sa valeur de defaut, soit 1024,  avant de faire la mise a jour.
    implicit none
    !Arguments:
    logical,              intent(in)     :: RESETQC       !reset the quality control flags before adding the new ones ?
    integer,              intent(in)     :: KCHKPRF(:)    !indicateur global controle de qualite tovs. Code:
    integer,              intent(inout)  :: globMarq(:)   !Marqueurs globaux  
    !Locals
    
    integer                              :: dataNum 
    integer                              :: dataIndex
    logical                              :: debug

    debug = mwbg_debug
    dataNum = size(globMarq)
    do dataIndex = 1, dataNum
      if (RESETQC) then
        globMarq(dataIndex) = 1024  
      end if
      if ( KCHKPRF(dataIndex).NE.0  ) then
        globMarq(dataIndex) = OR (globMarq(dataIndex),2**6)
      end if
    end do
    if (debug) then
      write(*,*) ' KCHKPRF   = ', (KCHKPRF(dataIndex),dataIndex=1,dataNum)
      write(*,*) ' NEW FLAGS = ', (globMarq(dataIndex),dataIndex=1,dataNum)
    end if

  end  subroutine resetQcCases

  !--------------------------------------------------------------------------
  !  setTerrainTypeToSeaIce
  !--------------------------------------------------------------------------
  subroutine setTerrainTypeToSeaIce(GLINTRP, KTERMER, ITERRAIN)

    !:Purpose:       Dans les conditions suivantes:
    !                1) l'indicateur terre/mer indique l'ocean (ktermer=1),
    !                2) le "terrain type" est manquant (iterrain=-1),
    !                3) le modele indique de la glace (gl >= 0.01),
    !                on specifie "sea ice" pour le "terrain type" (iterrain=0).
    implicit none 
    !Arguments:
    real,                 intent(in)     :: GLINTRP(:)    !sea ice
    integer,              intent(in)     :: KTERMER(:)    !land sea qualifier
    integer,              intent(inout)  :: ITERRAIN(:)   !terrain type
    !Locals
    
    integer                              :: dataNum 
    integer                              :: dataIndex
    logical                              :: debug

    debug = mwbg_debug
    dataNum = size(ITERRAIN)

    if ( mwbg_debug ) then
      write(*,*) ' OLD TERRAIN type = ', (ITERRAIN(dataIndex),dataIndex=1,dataNum)
      write(*,*) ' KTERMER = ', (KTERMER(dataIndex),dataIndex=1,dataNum)
      write(*,*) ' GLINTRP = ', (GLINTRP(dataIndex),dataIndex=1,dataNum)
    end if
    do dataIndex = 1, dataNum
      if ( KTERMER (dataIndex) == 1 .and. ITERRAIN(dataIndex) == -1 .and. GLINTRP (dataIndex) >= 0.01 ) &
           ITERRAIN(dataIndex) = 0
    end do
    if ( mwbg_debug ) then
      write(*,*) ' NEW TERRAIN type = ', (ITERRAIN(dataIndex),dataIndex=1,dataNum)
    end if
    
  end  subroutine setTerrainTypeToSeaIce

  !--------------------------------------------------------------------------
  !  addIntegerElementBurpBlock
  !--------------------------------------------------------------------------
  subroutine addIntegerElementBurpBlock(burpBlock, burpElement, burpArr, burpChannelNum, burpLocationNum, burpNele)

    !:Purpose: This subroutine takes a block, the btyp and the Intger element 
    !          and that we want to be added to the given block
    !          The element can be an array of element

    implicit none 
    ! Arguments:
    type(BURP_BLOCK),           intent(inout)       :: burpBlock               ! burp report
    integer,                    intent(in)          :: burpElement             ! burp element 
    integer,                    intent(in)          :: burpArr(:)              ! burp integer array to add
    integer,                    intent(in)          :: burpLocationNum         ! nlocs
    integer,                    intent(in)          :: burpChannelNum          ! nchannels
    integer,                    intent(in)          :: burpNele 
    
    ! Locals
    integer                                         :: ref_blk 
    logical                                         :: debug
    integer                                         :: burpLocationIndex
    integer                                         :: dataIndex
    integer                                         :: burpChannelIndex
    integer                                         :: error                   ! error status

    debug = mwbg_debug
    Call BURP_Resize_Block(burpBlock, ADD_NELE = 1, iostat = error)
    if (error /= burp_noerr) call utl_abort('bgckMicrowave_mod: BURP_Resize_Block Error')

    ! Add new elements CLW and scatw
    Call BURP_Set_Element(burpBlock,NELE_IND=burpNele,ELEMENT=burpElement,iostat=error)
    call BURP_Encode_Block(burpBlock)   ! encode the element numbers in the block
        
    dataIndex = 1 
    do burpChannelIndex = 1,burpChannelNum  
      do burpLocationIndex  =1, burpLocationNum
        Call BURP_Set_Tblval(burpBlock,NELE_IND=burpNele,NVAL_IND=burpChannelIndex,NT_IND=burpLocationIndex,TBLVAL=burpArr(dataIndex),iostat=error)
        dataIndex = dataIndex + 1        
      end do
    end do 
    !Call BURP_Convert_Block(burpBlock)

  end subroutine addIntegerElementBurpBlock

  !--------------------------------------------------------------------------
  !  addRealElementBurpBlock
  !--------------------------------------------------------------------------
  subroutine addRealElementBurpBlock(burpBlock, burpElement, burpArr, burpChannelNum, burpLocationNum, burpNele)

    !:Purpose: This subroutine takes a block, the btyp and the Real element 
    !          and that we want to be added to the given block
    !          The element can be an array of element

    implicit none 
    ! Arguments:
    type(BURP_BLOCK),           intent(inout)       :: burpBlock               ! burp report
    integer,                    intent(in)          :: burpElement             ! burp element 
    real,                       intent(in)          :: burpArr(:)              ! burp Real array to add
    integer,                    intent(in)          :: burpLocationNum         ! nlocs
    integer,                    intent(in)          :: burpChannelNum          ! nchannels
    integer,                    intent(in)          :: burpNele 
    
    ! Locals
    integer                                         :: ref_blk 
    logical                                         :: debug
    integer                                         :: burpLocationIndex
    integer                                         :: dataIndex
    integer                                         :: burpChannelIndex
    integer                                         :: error                   ! error status

    debug = mwbg_debug
    Call BURP_Resize_Block(burpBlock, ADD_NELE = 1, iostat = error)
    if (error /= burp_noerr) call utl_abort('bgckMicrowave_mod: BURP_Resize_Block Error')
    ! Add new elements CLW and scatw
    Call BURP_Set_Element(burpBlock,NELE_IND=burpNele,ELEMENT=burpElement,iostat=error)
    call BURP_Encode_Block(burpBlock)   ! encode the element numbers in the block
    dataIndex = 1 
    do burpChannelIndex = 1,burpChannelNum  
      do burpLocationIndex  =1, burpLocationNum
        Call BURP_Set_Rval(burpBlock,NELE_IND=burpNele,NVAL_IND=burpChannelIndex,NT_IND=burpLocationIndex,RVAL=burpArr(dataIndex),iostat=error)
        dataIndex = dataIndex + 1        
      end do
    end do 
    Call BURP_Convert_Block(burpBlock)

  end subroutine addRealElementBurpBlock

  !--------------------------------------------------------------------------
  ! copyIntegerElementToBurpBlock 
  !--------------------------------------------------------------------------
  subroutine copyIntegerElementToBurpBlock(burpBlock, burpElement, burpArr, burpArrName, burpChannelNum, burpLocationNum)

    !:Purpose: This subroutine takes a block, the btyp and the Real element 
    !          and that we want to be added to the given block
    !          The element can be an array of element

    implicit none 
    ! Arguments:
    type(BURP_BLOCK),           intent(inout)       :: burpBlock               ! burp report
    integer,                    intent(in)          :: burpElement             ! burp element 
    integer,                    intent(in)          :: burpArr(:)              ! burp Real array to add
    character (*),              intent(in)          :: burpArrName             ! burp array name
    integer,                    intent(in)          :: burpLocationNum         ! nlocs
    integer,                    intent(in)          :: burpChannelNum          ! nchannels
    
    ! Locals
    integer                                         :: ref_blk 
    integer                                         :: indice 
    logical                                         :: debug
    integer                                         :: burpLocationIndex
    integer                                         :: burpChannelIndex
    integer                                         :: error                  
    integer                                         :: idata 

    debug = mwbg_debug

    indice = BURP_Find_Element(burpBlock,burpElement,iostat=error)
    if (error /= burp_noerr) call utl_abort('bgckMicrowave_mod: BURP_Find_Element Error')
    if ( indice > 0 ) then
      do burpChannelIndex = 1, burpChannelNum
        do burpLocationIndex = 1, burpLocationNum
            idata = burpArr(burpLocationIndex)
            Call BURP_Set_Tblval(burpBlock,indice,burpChannelIndex,burpLocationIndex,idata)
        end do
      end do
    else
      call utl_abort('bgckMicrowave_mod: 3D block (btyp=5120) missing element '//trim(burpArrName))
    end if

  end subroutine copyIntegerElementToBurpBlock

  !--------------------------------------------------------------------------
  ! copyRealElementToBurpBlock 
  !--------------------------------------------------------------------------
  subroutine copyRealElementToBurpBlock(burpBlock, burpElement, burpArr, burpArrName, burpChannelNum, burpLocationNum)

    !:Purpose: This subroutine takes a block, the btyp and the Real element 
    !          and that we want to be added to the given block
    !          The element can be an array of element

    implicit none 
    ! Arguments:
    type(BURP_BLOCK),           intent(inout)       :: burpBlock               ! burp report
    integer,                    intent(in)          :: burpElement             ! burp element 
    real,                       intent(in)          :: burpArr(:)              ! burp Real array to add
    character (*),              intent(in)          :: burpArrName             ! burp array name
    integer,                    intent(in)          :: burpLocationNum         ! nlocs
    integer,                    intent(in)          :: burpChannelNum          ! nchannels
    
    ! Locals
    integer                                         :: ref_blk 
    integer                                         :: indice 
    logical                                         :: debug
    integer                                         :: burpLocationIndex
    integer                                         :: burpChannelIndex
    integer                                         :: error                  
    real                                            :: idata 

    debug = mwbg_debug

    indice = BURP_Find_Element(burpBlock,burpElement,iostat=error)
    if (error /= burp_noerr) call utl_abort('bgckMicrowave_mod: BURP_Find_Element Error')
    if ( indice > 0 ) then
      do burpChannelIndex = 1, burpChannelNum
        do burpLocationIndex = 1, burpLocationNum
            idata = burpArr(burpLocationIndex)
            Call BURP_Set_Rval(burpBlock,indice,burpChannelIndex,burpLocationIndex,idata)
        end do
      end do
    else
      call utl_abort('bgckMicrowave_mod: 3D block (btyp=5120) missing element '//trim(burpArrName))
    end if

  end subroutine copyRealElementToBurpBlock

  !--------------------------------------------------------------------------
  !  modifyBurpBktypAndWriteReport
  !--------------------------------------------------------------------------
  subroutine modifyBurpBktypAndWriteReport(burpReport, burpBlock, newBktyp, convertBlock)
    !
    !:Purpose:This subroutine takes a block, and its current bktyp 
    !          and then modify the bktyp then write the block to the 
    !          report

    implicit none 
    ! Arguments:
    type(BURP_BLOCK),           intent(inout)       :: burpBlock               ! burp block
    type(BURP_RPT),             intent(inout)       :: burpReport              ! burp report
    integer,                    intent(in)          :: newBktyp                ! burp element 
    logical,                    intent(in)          :: convertBlock            ! convert block yes or no
    
    ! Locals
    integer                                         :: error                  
    type(BURP_BLOCK)                                :: burpBlock_copy          ! burp block

    ! Remplacer le bloc et modifier le bktyp pour signifier "vu par AO".
    Call BURP_Set_Property(burpBlock,BKTYP=newBktyp,iostat=error)
    burpBlock_copy = burpBlock
    Call BURP_Write_Block(burpReport,BLOCK = burpBlock_copy ,CONVERT_BLOCK=convertBlock,iostat=error)
    if (error /= burp_noerr) call utl_abort('bgckMicrowave_mod: BURP_Write_Block Error')
    

  end subroutine modifyBurpBktypAndWriteReport

  !--------------------------------------------------------------------------
  !  mwbg_updateBurp
  !--------------------------------------------------------------------------
  subroutine mwbg_updateBurp(burpFileNameIn, ReportIndex, ETIKRESU, ztb, clw, scatw, ident, &
                               globMarq, RESETQC, KCHKPRF, KTERMER, ITERRAIN, IMARQ, writeTbValuesToFile, & 
                               writeModelLsqTT, writeEle25174, burpFileNameout)
    !:Purpose:      Pour AMSUA et ATMS: 
    !               Allumer les bits des marqueurs pour les tovs rejetes.
    !               Mettre a jour l'indicateur terre/mer qui a
    !               possiblement ete corrige pour la glace marine.
    !               Modifier le bktyp des donnees, marqueurs et (O-P) pourt
    !               signifier "vu par AO". 
    implicit none
    !Arguments:
    !
    character(len=90),    intent(in)     :: burpFileNameIn     !burp input Obs File
    integer,              intent(in)     :: ReportIndex        !report eportIndex
    character(len=9),     intent(in)     :: ETIKRESU           !resume report etiket
    real,                 intent(in)     :: ztb(:)             ! tempertature de brillance
    real,                 intent(in)     :: clw(:)             !cloud liquid water path 
    real,                 intent(in)     :: scatw(:)           !scattering index 
    integer,              intent(inout)  :: globMarq(:)        !Marqueurs globaux  
    integer, intent(in)                  :: ident(:)           !Information flag (ident) values (new BURP element 025174 in header)
    integer,              intent(in)     :: KTERMER(:)         !indicateur terre/mer 
    integer,              intent(inout)  :: ITERRAIN(:)        !indicateur du type de terrain 
                                                               !niveau de chaque canal  
    logical,              intent(in)     :: RESETQC            !reset the quality control flags before adding the new ones ?
    integer,              intent(in)     :: KCHKPRF(:)         !indicateur global controle de qualite tovs. Code:
                                                               !=0, ok,
                                                               !>0, rejet, 
    logical, intent(in)                  :: writeTbValuesToFile! if replace missing tb by missing value
    logical, intent(in)                  :: writeModelLsqTT    ! if replace lsq and terrain typ var by model values
    logical, intent(in)                  :: writeEle25174      ! if copie the new element 25174 in the output burp file
    integer,              intent(in)     :: IMARQ(:)           !modified flag values from mwbg_tovCheck  
    character(len=90),    intent(inout)  :: burpFileNameout    ! burp output Obs File

    !Locals:
    type(BURP_FILE), save                :: File_in
    type(BURP_FILE), save                :: File_out
    type(BURP_RPT)                       :: reportIn
    type(BURP_RPT)                       :: reportOut

    type(BURP_BLOCK)                     :: blk 
    type(BURP_BLOCK)                     :: blk_copy 

    logical, save                        :: ifFirstCall=.True.
    logical                              :: resumeReport 
    character(len=9)                     :: idStn
    integer, allocatable, save           :: adresses(:)        ! First Loop over all reports to get their adresses and save them
    integer, save                        :: nb_rpts
    integer, save                        :: nsize

    integer                              :: error
    integer                              :: ref_blk
    integer                              :: my_nt
    integer                              :: my_nval
    integer                              :: my_nele
    integer                              :: my_idtyp
    integer                              :: idtyp
    integer                              :: nlocs
    integer                              :: blat
    integer                              :: blon
    integer                              :: nblocs
    integer                              :: handle
    integer                              :: iun_burpin
    integer                              :: my_btyp
    integer                              :: my_bktyp 
    integer                              :: my_bfam
    integer                              :: new_bktyp
    integer                              :: indice 
    integer                              :: kk
    integer                              :: jj
    integer                              :: JI
    integer                              :: ipos
    integer                              :: idata 
    integer                              :: j


    if (ifFirstCall) then
      write(*,*) 'mwbg_updateBurp First Call : Initialisation'
      ! LOOP OVER ALL REPORTS OF THE INPUT FILE, APPLY PROCESSING, AND write TO OUTPUT FILE.
      call getBurpReportAdresses (burpFileNameIn, adresses)
      ! initialisation
      Call BURP_Init(File_in, F2= File_out, iostat=error)
      ! ouverture du fichier burp d'entree et de sortie
      Call BURP_Init(reportIn, R2=reportOut, iostat=error)
      ! Number of reports and maximum report size from input BURP file
      Call BURP_New(File_in,  FILENAME = burpFileNameIn,  MODE= FILE_ACC_READ,   iostat= error)
      Call BURP_Get_Property(File_in, NRPTS=nb_rpts, IO_UNIT= iun_burpin)
      nsize = MRFMXL(iun_burpin)
      ! Add nsize to report size to accomodate modified (larger) data blocks
      nsize = nsize*3
      Call BURP_New(File_out, FILENAME= burpFileNameOut, MODE= FILE_ACC_CREATE, iostat= error)
      ifFirstCall = .False.
    end if

    Call BURP_Get_Report(File_in, REPORT= reportIn, REF= adresses(ReportIndex), iostat= error) 
    if (error /= burp_noerr) call burpErrorHistory(file_in, reportIn) 
    Call BURP_Get_Property(reportIn,STNID=idStn,IDTYP=idtyp,ELEV=nlocs,LATI=blat,LONG=blon,NBLK=nblocs,HANDLE=handle)

    if ( idStn(1:2) .eq. ">>" ) then
      resumeReport = .True.
      ! change the header
      Call BURP_Set_Property(reportIn,STNID=ETIKRESU)  
      Call BURP_Write_Report(File_out,reportIn,iostat=error)
      if (error /= burp_noerr) call burpErrorHistory(File_in, reportIn, File_out, reportOut)

      if (reportIndex == nb_rpts) then
        Call BURP_Free(File_in, F2 = File_out,iostat=error)
        Call BURP_Free(reportIn, R2 = reportOut,iostat=error)
      end if
      
      return 
    end if  
    ! Create new report (reportOut) to contain modified blocks from reportIn
    Call BURP_New(reportOut, Alloc_Space = nsize,  iostat=error)
    if (error /= burp_noerr) call burpErrorHistory(file_in, reportIn, file_out, reportOut)
    
    ! initiliser pour ecriture a File_out
    Call BURP_INIT_Report_Write(File_out,reportOut,iostat=error)
    if (error /= burp_noerr) call burpErrorHistory(File_in, reportIn, File_out, reportOut)
    !  copier le header du rapport 
    Call BURP_Copy_Header(TO= reportOut, FROM= reportIn)

    Call BURP_Init(blk, B2=blk_copy, iostat=error)

    ! Read and modify the blocks in rpt and add them to reportOut
    ref_blk = 0
    BLOCKS: do

      ref_blk = BURP_Find_Block(reportIn,BLOCK= blk,SEARCH_FROM= ref_blk,iostat= error)
      if (error /= burp_noerr) call utl_abort('bgckMicrowave_mod: BURP_Find_Block Error')
      
      if (ref_blk < 0) Exit BLOCKS

      Call BURP_Get_Property(blk, &
                  NELE   = my_nele, &
                  NVAL   = my_nval, &       ! 1 or number of channels (obs per location) if Tb data/flag block
                  NT     = my_nt, &         ! 1 or number of locations in block
                  BTYP   = my_btyp, &
                  BKTYP  = my_bktyp, & 
                  BFAM   = my_bfam, &
                  iostat = error)
      if (error /= burp_noerr) call utl_abort('bgckMicrowave_mod: BURP_Get_Porperty Error')

 
      ! 1) Bloc info 3d: bloc 5120.
      !    Modifier les marqueurs globaux de 24bits pour les donnees rejetees.
      if (my_btyp == 5120) then     

        if (mwbg_debug) then
          write(*,*) ' OLD FLAGS = ', (globMarq(JJ),JJ=1,my_nt)
        end if

        call resetQcCases(RESETQC, KCHKPRF, globMarq)
        ! Remplacer les nouveaux marqueurs dans le tableau.
        call copyIntegerElementToBurpBlock(blk, 55200, globMarq, "MARQUEURS_GLOBAUX", my_nval, my_nt)
        call modifyBurpBktypAndWriteReport(reportOut, blk, my_bktyp,.False.)        
      ! 2) Bloc info (general): bloc 3072
      !    Modifier les indicateurs terre/mer possiblement corriges pour la glace
      !    marine.
      !    Add new elements CLW, scatw and ident

      else if (my_btyp == 3072) then
        if (writeEle25174) then 
          ! Add new elements IDENT
          call addIntegerElementBurpBlock(blk, 25174, ident, my_nval, my_nt, my_nele+1)
          ! Add new elements CLW
          call addRealElementBurpBlock(blk, 13209, clw, my_nval, my_nt, my_nele+2)
          ! Add new elements scatw
          call addRealElementBurpBlock(blk, 13208, scatw, my_nval, my_nt, my_nele+3)
        else
          call addRealElementBurpBlock(blk, 13209, clw, my_nval, my_nt, my_nele+1)
          ! Add new elements scatw
          call addRealElementBurpBlock(blk, 13208, scatw, my_nval, my_nt, my_nele+2)          
        end if 
        ! Remplacer les nouveaux indicateurs terre/mer dans le tableau.
        if (writeModelLsqTT) then        
          call copyIntegerElementToBurpBlock(blk, 8012, KTERMER, "Land_Sea_INDICATOR", my_nval, my_nt)
          ! Update terrain type dans le tableau.
          call copyIntegerElementToBurpBlock(blk, 13039, ITERRAIN, "TERRAIN_type", my_nval, my_nt)
          ! finally copy blk
        end if
        call modifyBurpBktypAndWriteReport(reportOut, blk, my_bktyp,.False.)        

      ! 3) Bloc multi niveaux de radiances: bloc 9218, 9248, 9264.
      !    Modifier le bktyp pour signifier "vu par AO".


      else if ( (my_btyp == 9218 .or. my_btyp == 9248 .or. my_btyp ==9264) .and. &
                my_bfam == 0 ) then 
        ! Modify Tb data if any data (ztb) were set to mwbg_realMisg

        if (writeTbValuesToFile) then
          call copyRealElementToBurpBlock(blk, 12163, ztb, "Tb_data", my_nval, my_nt)
          call modifyBurpBktypAndWriteReport(reportOut, blk, my_bktyp+4, .False.)        
        else
          call modifyBurpBktypAndWriteReport(reportOut, blk, my_bktyp+4, .True.)        
        end if



      ! 4) Bloc marqueurs multi niveaux de radiances: bloc 15362, 15392, 15408.
      !    Modifier les marqueurs de 13bits associes a chaque radiance.
      !    Modifier le bktyp pour signifier "vu par AO".
      else if ( (my_btyp == 15362 .or. my_btyp == 15392 .or. my_btyp == 15408) .and. &
                my_bfam == 0 ) then 

        ! extraire les marqueurs de 13bits des radiances; element 212163 (LEVEL 1B)
        indice = BURP_Find_Element(blk,212163,iostat=error)
        ! Remplacer les nouveaux marqueurs dans le tableau.
        ipos = 0
        do kk =1, my_nt
          do j = 1, my_nval
            ipos = ipos + 1
            Call BURP_Set_Tblval(blk,indice,j,kk,IMARQ(ipos))
          end do
        end do
        ! Remplacer le bloc et modifier le bktyp pour signifier "vu par AO".
        call modifyBurpBktypAndWriteReport(reportOut, blk, my_bktyp+4, .True.)        


      ! 5) Bloc multi niveaux de residus de radiances (O-P): bloc 9322, 9226, 9258, 9274, bfam 14
      !    Modifier le bktyp pour signifier "vu par AO".
      else if ( (my_btyp == 9322 .or. my_btyp == 9226 .or. my_btyp == 9258 .or. &
                my_btyp == 9274) .and. my_bfam == 14 ) then 
        ! Remplacer le bloc et modifier le bktyp pour signifier "vu par AO".
        call modifyBurpBktypAndWriteReport(reportOut, blk, my_bktyp+4, .True.)
      ! OTHER BLOCK 
      else
        call modifyBurpBktypAndWriteReport(reportOut, blk, my_bktyp, .True.)

      end if

    end do BLOCKS
    ! Write the modified report to the output file
    Call BURP_Write_Report(File_out,reportOut,iostat=error)
    if (error /= burp_noerr) call burpErrorHistory(File_in, reportIn, File_out, reportOut)

    if (reportIndex == nb_rpts) then
      write(*,*) "LAST REPORT : Burp Udate is closing File"
      Call BURP_Free(File_in, F2 = File_out,iostat=error)
      Call BURP_Free(reportIn, R2 = reportOut,iostat=error)
    end if 
  end subroutine mwbg_updateBurp

  !--------------------------------------------------------------------------
  !  GRODY
  !--------------------------------------------------------------------------
  subroutine GRODY (ier, ni, tb23, tb31, tb50, tb53, tb89, &
                   tb23_P, tb31_P, tb50_P, tb53_P, tb89_P, &
                   pangl, plat, ilansea, ice, tpw, clw, clw_avg, &
                   rain, snow, scatl, scatw)
    !OBJET          Compute the following parameters using 5 AMSU-A
    !               channels:
    !                  - sea ice, 
    !                  - total precipitable water, 
    !                  - cloud liquid water, 
    !                  - ocean/land rain, 
    !                  - snow cover/glacial ice,
    !                  - scattering index (sur la terre et sur l'eau).
    !               The four channels used are: 23Ghz, 31Ghz, 50Ghz and 89Ghz.
    !
    !REGERENCES     N. Grody, NOAA/NESDIS, ....
    !
    !APPEL          call   GRODY (ier, ni, tb23, tb31, tb50, tb53, tb89, pangl, plat,
    !                             ilansea, ice, tpw, clw, rain, snow, scatl, scatw) 
    !
    !ARGUMENTS      ier     - output - error return code:
    !                                  0, ok,  
    !                                  1, input parameter out of range. 
    !               - ni      - input  -  number of points to process
    !               - tb23    - input  -  23Ghz brightness temperature (K)
    !               - tb31    - input  -  31Ghz brightness temperature (K)
    !               - tb50    - input  -  50Ghz brightness temperature (K)
    !               - tb53    - input  -  53Ghz brightness temperature (K)
    !               - tb89    - input  -  89Ghz brightness temperature (K)
    !               - tb23_P  - input  -  23Ghz brightness temperature from background (K)
    !               - tb31_P  - input  -  31Ghz brightness temperature from background (K)
    !               - tb50_P  - input  -  50Ghz brightness temperature from background (K)
    !               - tb53_P  - input  -  53Ghz brightness temperature from background (K)
    !               - tb89_P  - input  -  89Ghz brightness temperature from background (K)
    !               - pangl   - input  -  satellite zenith angle (deg.)
    !               - plat    - input  -  lalitude (deg.)
    !               - ilansea - input  -  land/sea indicator (0=land;1=ocean)
    !               - ice     - output -  sea ice concentration (0-100%)
    !               - tpw     - output -  total precipitable water (0-70mm)
    !               - clw     - output -  cloud liquid water (0-3mm)
    !               - clw_avg - output -  averaged cloud liquid water from obs and 
    !                                   background (0-3mm)
    !               - rain    - output -  rain identification (0=no rain; 1=rain)
    !               - snow    - output -  snow cover and glacial ice identification: 
    !                                   (0=no snow; 1=snow; 2=glacial ice)
    !               - scatl   - output -  scattering index over land
    !               - scatw   - output -  scattering index over water
    !
    ! Notes: In the case where an output parameter cannot be calculated, the
    !        value of this parameter is to to the missing value, i.e. -99.

    implicit none

    integer ni, i

    integer ier    (:)
    integer ilansea(:)
    integer rain   (:)
    integer snow   (:)

    real mwbg_realMisgLocal, siw, sil, df1, df2, df3, a, b, c, d, e23
    real ei, cosz, tt, scat, sc31, abslat, t23, t31, t50, t89
    real sc50, par, t53
    real dif285t23, dif285t31, epsilon
    real dif285t23_P, dif285t31_P

    real tb23  (:)
    real tb31  (:)
    real tb50  (:)
    real tb53  (:)
    real tb89  (:)
    real tb23_P(:)
    real tb31_P(:)
    real tb50_P(:)
    real tb53_P(:)
    real tb89_P(:)
    real pangl (:)
    real plat  (:)
    real ice   (:)
    real tpw   (:)
    real clw   (:)
    real clw_P
    real clw_avg(:)
    real scatl (:)
    real scatw (:)

    data mwbg_realMisgLocal   / -99.     /
    data epsilon /   1.E-30 /

    ! 1) Initialise output parameters:
    do i = 1, ni
      ice  (i) = mwbg_realMisgLocal
      tpw  (i) = mwbg_realMisgLocal
      clw  (i) = mwbg_realMisgLocal
      clw_avg(i) = mwbg_realMisgLocal
      scatl(i) = mwbg_realMisgLocal
      scatw(i) = mwbg_realMisgLocal
      rain (i) = nint(mwbg_realMisgLocal)
      snow (i) = nint(mwbg_realMisgLocal)
    end do

    ! 2) Validate input parameters:
    do i = 1, ni
      if ( tb23(i)    .lt. 120.  .or. &
           tb23(i)    .gt. 350.  .or. &
           tb31(i)    .lt. 120.  .or. &
           tb31(i)    .gt. 350.  .or. &
           tb50(i)    .lt. 120.  .or. &
           tb50(i)    .gt. 350.  .or. &
           tb53(i)    .lt. 120.  .or. &
           tb53(i)    .gt. 350.  .or. &
           tb89(i)    .lt. 120.  .or. &
           tb89(i)    .gt. 350.  .or. &
           pangl(i)   .lt. -90.  .or. &
           pangl(i)   .gt.  90.  .or. &
           plat(i)    .lt. -90.  .or. &
           plat(i)    .gt.  90.  .or. &
           ilansea(i) .lt.   0   .or. &
           ilansea(i) .gt.   1        ) then
        ier(i) = 1
      else
        ier(i) = 0
      end if

    end do

    !3) Compute parameters:
    do i = 1, ni
      if ( ier(i) .eq. 0 ) then
        abslat = abs(plat(i))
        cosz   = cosd(pangl(i))
        t23 = tb23(i)
        t31 = tb31(i)
        t50 = tb50(i)
        t53 = tb53(i)
        t89 = tb89(i)
        dif285t23  =max(285.-t23,epsilon)
        dif285t23_P=max(285.-tb23_P(i),epsilon)
        dif285t31  =max(285.-t31,epsilon)
        dif285t31_P=max(285.-tb31_P(i),epsilon)

        ! scattering indices:
        siw = -113.2 + (2.41 - 0.0049*t23)*t23 &
                     + 0.454*t31 &
                     - t89 
        sil = t23 - t89 

        scatl (i) = sil
        scatw (i) = siw

        ! discriminate functions:
        df1 =  2.85 + 0.020*t23 - 0.028*t50 ! used to identify (also remove) sea ice
        df2 =  5.10 + 0.078*t23 - 0.096*t50 ! used to identify (also remove) warm deserts
        df3 = 10.20 + 0.036*t23 - 0.074*t50 ! used to identify (also remove) cold deserts

        if ( ilansea(i) .eq. 1 ) then

          ! Ocean Parameters

          !3.1) Sea ice:
          if ( abslat .lt. 50. ) then
            ice(i) = 0.0
          else
            if ( df1 .lt. 0.45 ) then
               ice(i) = 0.0
            else
              a =  1.7340  - 0.6236*cosz
              b =  0.0070  + 0.0025*cosz
              c = -0.00106 
              d = -0.00909
              e23 = a + b*t31 + c*t23 + d*t50 ! theoretical 23Ghz sfc emissivity (0.3-1.)
              if ( (t23-t31) .ge. 5. ) then   ! fov contains multiyear or new ice/water
                ei = 0.88
              else
                ei = 0.95
              end if
              ice(i) = 100*(e23-0.45)/(ei-0.45) ! sea-ice concentration within fov (0-100%) 
              ice(i) = min(100.,max(0.,ice(i)))/100.   !jh (0.-1.)
            end if
          end if

          ! 3.2) Total precipitable water:
          ! identify and remove sea ice
          if ( abslat .gt. 50.  .and. &
              df1    .gt.  0.2        ) then  
            tpw(i) = mwbg_realMisgLocal
          else
            a =  247.920  - (69.235 - 44.177*cosz)*cosz
            b = -116.270
            c =   73.409
            tpw(i) = a + b*log(dif285t23) & 
                      + c*log(dif285t31)
            tpw(i) = tpw(i)*cosz           ! theoretical total precipitable water (0-70mm)
            tpw(i) = 0.942*tpw(i) - 2.17   ! corrected   total precipitable water 
            tpw(i) = min(70.,max(0.,tpw(i)))   ! jh     
          end if

          !3.3) Cloud liquid water from obs (clw) and background state (clw_P):
          ! identify and remove sea ice
          if ( abslat .gt. 50.  .and. &
              df1    .gt.  0.0        ) then  
            clw(i) = mwbg_realMisgLocal
            clw_avg(i) = mwbg_realMisgLocal
          else
            a =  8.240 - (2.622 - 1.846*cosz)*cosz
            b =  0.754
            c = -2.265
            clw(i) = a + b*log(dif285t23) & 
                      + c*log(dif285t31)
            clw(i) = clw(i)*cosz           ! theoretical cloud liquid water (0-3mm)
            clw(i) = clw(i) - 0.03         ! corrected   cloud liquid water 
            clw(i) = min(3.,max(0.,clw(i)))   ! jh       

            clw_P = a + b*log(dif285t23_P) & 
                      + c*log(dif285t31_P)
            clw_P = clw_P*cosz           ! theoretical cloud liquid water (0-3mm)
            clw_P = clw_P - 0.03         ! corrected   cloud liquid water 
            clw_P = min(3.,max(0.,clw_P))   ! jh       

            ! averaged CLW from observation and background
            clw_avg(i) = 0.5 * (clw(i) + clw_P)
          end if

          !3.4) Ocean rain: 0=no rain; 1=rain.
          ! identify and remove sea ice
          if ( abslat .gt. 50.  .and. &
              df1    .gt.  0.0        ) then  
            rain(i) = nint(mwbg_realMisgLocal)
          else                                   ! remove non-precipitating clouds
            if ( clw(i) .gt. 0.3 .or. &
                siw    .gt. 9.0      ) then 
              rain(i) = 1
            else
              rain(i) = 0
            end if
          end if

        else

          ! Land Parameters

          ! 3.5) Rain  over land: 0=no rain; 1=rain.
          tt = 168. + 0.49*t89
          if ( sil .ge. 3. ) then
            rain(i) = 1
          else 
            rain(i) = 0
          end if
          
          ! remove snow cover
          if ( t23 .le. 261. .and. &
              t23 .lt. tt         ) then
            rain(i) = 0
          end if

          ! remove warm deserts
          if ( t89 .gt. 273.  .or. &
              df2 .lt.   0.6      ) then
            rain(i) = 0
          end if

          ! 3.6) Snow cover and glacial ice: 0=no snow; 1=snow; 2=glacial ice.
          tt = 168. + 0.49*t89
          scat = t23 - t89
          sc31 = t23 - t31
          sc50 = t31 - t50
          par  = t50 - t53

          ! re-frozen snow
          if ( t89  .lt. 255.  .and. &
              scat .lt. sc31        ) then
            scat = sc31
          end if

          ! identify glacial ice
          if ( scat .lt.   3.  .and. &
              t23  .lt. 215.        ) then
            snow(i) = 2
          end if
          if ( scat .ge. 3.       ) then
            snow(i) = 1
          else
            snow(i) = 0
          end if

          ! remove precipitation
          if ( t23 .ge. 262.  .or. &
              t23 .ge. tt           ) then
            snow(i) = 0
          end if
          if ( df3 .le. 0.35         ) then    ! remove deserts
            snow(i) = 0
          end if

          ! high elevation deserts
          if ( scat .lt.  15.  .and. &
              sc31 .lt.   3.  .and. &
              par  .gt.   2.        ) then
            snow(i) = 0
          end if

          ! remove frozen ground
          if ( scat .lt.   9.  .and. &
              sc31 .lt.   3.  .and. &
              sc50 .lt.   0.        ) then
            snow(i) = 0
          end if

        end if

      end if

      if ( mwbg_DEBUG .and. i <= 100 ) then
        print *, 'GRODY: i,tb23(i),tb31(i),tb50(i),tb89(i),pangl(i),plat(i), &
                  ilansea(i) = ', &
                  i,tb23(i),tb31(i),tb50(i),tb89(i),pangl(i),plat(i), &
                  ilansea(i)
        print *, 'GRODY: ier(i),ice(i),tpw(i),clw(i),clw_avg(i),rain(i),snow(i)=', &
                  ier(i),ice(i),tpw(i),clw(i),clw_avg(i),rain(i),snow(i)
      end if

    end do

  end subroutine GRODY

  !--------------------------------------------------------------------------
  !  mwbg_readStatTovs
  !--------------------------------------------------------------------------
  subroutine mwbg_readStatTovs(statFileName, instName, satelliteId, IUTILST, TOVERRST)
    !:Purpose:       Lire les statistiques de l'erreur totale pour les TOVS.
    !
    implicit none
    !Arguments:
    character(*), intent(in) :: statFileName  ! fichier stats des TOVS
    character(*), intent(in) :: instName      ! Instrument Name
    character(len = 9), allocatable, intent(out) :: satelliteId(:)       ! Identificateur de satellite 
    integer,                         intent(out) :: IUTILST(mwbg_maxNumChan,mwbg_maxNumSat) ! channel Selection using array IUTILST(chan,sat)
    !                                                                      IUTILST = 0 (blacklisted)
    !                                                                      1 (assmilate)
    !                                                                      2 (assimilate over open water only)
    real,                         intent(out) :: TOVERRST(mwbg_maxNumChan,mwbg_maxNumSat) ! l'erreur totale des TOVS
 
    !Locals

    integer :: iunStat
    integer :: ier 
    integer :: ISTAT 
    integer :: fnom 
    integer :: ICHN 
    integer :: satNumber 
    integer :: jpnsatIndex
    integer :: jpnchanIndex
    integer, parameter :: casesNum = 2
    integer :: satIndex 
    integer :: casesIndex 
    integer :: INDX, IPOS
    integer :: alloc_status 
    integer :: NCHNA(mwbg_maxNumSat)
    integer :: MLISCHNA(mwbg_maxNumChan,mwbg_maxNumSat)
    real*8  :: tovErrorIn(mwbg_maxNumChan,2,mwbg_maxNumSat)
    real    :: tovErrorStatus(mwbg_maxNumChan,mwbg_maxNumSat)
    integer :: channelInNum(mwbg_maxNumSat)
    integer :: identifSatId(mwbg_maxNumSat)
    real*8 :: ZDUM
    
    character(132) :: CLDUM
    character(17) :: CSATSTR

    character(12)  :: satelliteStatsType(casesNum)

    DATA satelliteStatsType     / 'Monitoring',  'Assimilation'  /  
    
    iunStat = 0
    ier = fnom(iunStat,statFileName,'SEQ+FMT',0)
    if (ier .LT. 0) then
      call utl_abort ('bgckMicrowave_mod: Problem opening TOVS total error statistics file '//trim(statFileName))
    end if



    write(*,FMT=9000)
9000 FORMAT(//,10x,"-mwbg_readStatTovs: reading total error statistics" &
          ," required for TOVS processing")

    ! 1. Initialize
100  CONTINUE
    do jpnsatIndex = 1, mwbg_maxNumSat
      NCHNA(jpnsatIndex) = 0
      channelInNum(jpnsatIndex) = 0
      identifSatId(jpnsatIndex) = 0
      do jpnchanIndex = 1, mwbg_maxNumChan
        tovErrorIn(jpnchanIndex,1,jpnsatIndex) = 0.0
        tovErrorIn(jpnchanIndex,2,jpnsatIndex) = 0.0
        IUTILST (jpnchanIndex,jpnsatIndex) = 0
        tovErrorStatus(jpnchanIndex,jpnsatIndex) = 0.0
      end do
    end do

    ! 2. Open the file
200  CONTINUE
    ! .... not done here anymore, jh, august 2000

    ! 3. Print the file contents
300  CONTINUE

    write(*,'(20X,"ASCII dump of stats_tovs file: "//)')
    do jpnchanIndex = 1, 9999999
      read (iunStat,'(A)',ERR=900,end=400) CLDUM
      write(*,'(A)')   CLDUM
    end do

    ! 4. Read number of satellites
400  CONTINUE

    rewind(iunStat)
    read (iunStat,*,ERR=900)
    read (iunStat,*,ERR=900) satNumber
    read (iunStat,*,ERR=900)
    alloc_status = 0
    if (allocated(satelliteId)) deallocate(satelliteId)
    allocate(satelliteId(satNumber), stat=alloc_status)
    if (alloc_status /= 0) then
      call utl_abort('bgckMicrowave_mod: ERROR - allocate(satelliteId) in mwbg_readStatTovs')
    end if

    ! 5. Read the satellite identification, the number of channels,
    ! the observation errors and the utilization flags
500  CONTINUE

    do jpnsatIndex = 1, satNumber
      read (iunStat,*,ERR=900)
      read (iunStat,'(A)',ERR=900) CLDUM
      CSATSTR = TRIM(ADJUSTL(CLDUM))

      ! Get satellite (e.g. NOAA18) from satellite/instrument (e.g. NOAA18 AMSUA)
      INDX = INDEX(CSATSTR, instName)
      if ( INDX .GT. 3 ) then
        if ( (instName == 'AMSUA') .and. ( INDEX(CSATSTR,'EOS-2') .GT. 0 ) ) then
          satelliteId(jpnsatIndex) = 'AQUA'
        else
          IPOS = INDX-2
          satelliteId(jpnsatIndex) = CSATSTR(1:IPOS)
        end if
      else
        write ( *,*) '(" mwbg_readStatTovs: Non- "', instName, &
                   ' instrument found in stats file!")'
        write ( *,'(A)' ) CLDUM
        call utl_abort ('bgckMicrowave_mod:  mwbg_readStatTovs No instrument found in stats file for '//trim(instName))
      end if
      read (iunStat,*,ERR=900)
      read (iunStat,*,ERR=900) identifSatId(jpnsatIndex), channelInNum(jpnsatIndex)
      do jpnchanIndex = 1, 3
        read (iunStat,*,ERR=900)
      end do

      ! Set errors to ERRBGCK column values
      do jpnchanIndex = 1, channelInNum(jpnsatIndex)
        read (iunStat,*,ERR=900) ICHN, &
                  tovErrorIn(ICHN,1,jpnsatIndex), &
                  tovErrorIn(ICHN,2,jpnsatIndex), &
                  IUTILST (ICHN,jpnsatIndex), ZDUM
        TOVERRST(ICHN,jpnsatIndex) = tovErrorIn(ICHN,1,jpnsatIndex)
      end do
      read (iunStat,*,ERR=900)
    end do

510  CONTINUE
    ! 6. Print error stats for assimilated channels

600  CONTINUE

    write(*,'(//5X,"Total errors for TOVS data"/)') 
    do jpnsatIndex = 1, satNumber
      INDX = 0
      do jpnchanIndex = 1, mwbg_maxNumChan
        if ( IUTILST(jpnchanIndex,jpnsatIndex) .NE. 0 ) then
          NCHNA(jpnsatIndex) = NCHNA(jpnsatIndex) + 1
          INDX = INDX + 1
          MLISCHNA(INDX,jpnsatIndex) = jpnchanIndex
        end if
      end do
      do casesIndex = 1, CasesNum
        write(*,'(/7X,"Satellite: ",A,5X,A)') &
          satelliteId(jpnsatIndex), satelliteStatsType(casesIndex)
        write(*,'(7X,"Channels   : ",30(T22,27I4/))') & 
         (MLISCHNA(jpnchanIndex,jpnsatIndex),jpnchanIndex=1,NCHNA(jpnsatIndex))
        write(*,'(7X,"Total errors: ",30(T22,27f4.1/))') &
         (tovErrorIn(MLISCHNA(jpnchanIndex,jpnsatIndex),casesIndex,jpnsatIndex), &
          jpnchanIndex=1,NCHNA(jpnsatIndex))
      end do
    end do
     
    ! 7. Close the file
700  CONTINUE
    ISTAT = FCLOS (iunStat)
    ! .... not done here anymore, jh, august 2000
    write(*,*) " SATID's = "
    do jpnsatIndex = 1, satNumber
      write(*,*) '  ', satelliteId(jpnsatIndex)
    end do

    return

    ! Read error
900   write ( *, '(" mwbg_readStatTovs: Problem reading ", &
                     "TOVS total error stats file ")' ) 
    ISTAT = FCLOS (iunStat)
    call utl_abort ('bgckMicrowave_mod:  mwbg_readStatTovs Problem reading TOVS total error stats file')

  end subroutine mwbg_readStatTovs

  !--------------------------------------------------------------------------
  !  atmsTest1Flagbit7Check
  !--------------------------------------------------------------------------
  subroutine atmsTest1Flagbit7Check (itest, KCANO, KMARQ, KNOSAT, ICHECK, KNO, KNT, STNID,  B7CHCK, rejectionCodArray)

    !:Purpose:               1) test 1: Check flag bit 7 on from the first bgckAtms program
    !                        Includes observations flagged for cloud liquid water, scattering index,
    !                        dryness index plus failure of several QC checks.


    ! Arguments
    integer,     intent(in)                                   :: itest(:)                 ! test number
    integer,     intent(in)                                   :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(inout)                                :: KMARQ(KNO,KNT)           ! observations channels
    integer,     intent(in)                                   :: KNOSAT                   ! numero de satellite (i.e. indice) 
    integer,     intent(in)                                   :: KNO                      ! nombre de canaux des observations 
    integer,     intent(in)                                   :: KNT                      ! nombre de tovs
    character *9,intent(in)                                   :: STNID                    ! identificateur du satellite
    integer,     intent(inout)                                :: B7CHCK(KNO,KNT)          ! 
    integer,     intent(inout)                                :: ICHECK(KNO,KNT)          ! indicateur du QC par canal
    integer,     intent(inout)                                :: rejectionCodArray(mwbg_maxNumTest,KNO,KNT) ! cumul of reject element 
    ! Locals
    integer                                                   :: channelval
    integer                                                   :: nDataIndex
    integer                                                   :: nChannelIndex
    integer                                                   :: testIndex
    integer                                                   :: IBIT 

    testIndex = 1
    if ( itest(testIndex) .eq. 1 ) then
      do nDataIndex=1,KNT
        do nChannelIndex=1,KNO
          IBIT = AND(KMARQ(nChannelIndex,nDataIndex), 2**7)
          if ( IBIT .NE. 0  ) then
            ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
            B7CHCK(nChannelIndex,nDataIndex) = 1
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
            rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                 rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
            if ( mwbg_debug ) then
              write(*,*)STNID(2:9),' first bgckAtms program REJECT.', &
                        'CHANNEL=', KCANO(nChannelIndex,nDataIndex), &
                        ' IMARQ= ',KMARQ(nChannelIndex,nDataIndex)
            end if
          end if
        end do
      end do
    end if

  end subroutine atmsTest1Flagbit7Check

  !--------------------------------------------------------------------------
  !  atmsTest2TopographyCheck
  !--------------------------------------------------------------------------
  subroutine atmsTest2TopographyCheck (itest, KCANO, KNOSAT, KNO, KNT, STNID, MTINTRP, &
                                   KMARQ, ICHTOPO, MXTOPO, ZCRIT, B7CHCK, & 
                                   ICHECK, rejectionCodArray, rejectionCodArray2)

    !:Purpose:               1) test 2: Topography check (partial)

    ! Arguments
    integer,     intent(in)                                   :: itest(:)                 ! test number
    integer,     intent(in)                                   :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)                                   :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)                                   :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)                                   :: KNT                            ! nombre de tovs
    character *9,intent(in)                                   :: STNID                          ! identificateur du satellite
    real,        intent(in)                                   :: MTINTRP(KNT)                   ! topo aux point d'obs
    integer,     intent(inout)                                :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(in)                                   :: MXTOPO 
    integer,     intent(in)                                   :: ICHTOPO(MXTOPO) 
    real ,       intent(in)                                   :: ZCRIT(MXTOPO)
    integer,     intent(inout)                                :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(inout)                                :: B7CHCK(KNO,KNT)          ! 
    integer,     intent(inout)                                :: rejectionCodArray(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 
    integer,     intent(inout)                                :: rejectionCodArray2(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 
    ! Locals
    integer                                                   :: channelval
    integer                                                   :: nDataIndex
    integer                                                   :: nChannelIndex
    integer                                                   :: testIndex
    integer                                                   :: INDXTOPO 

    testIndex = 2
    if ( itest(testIndex) .eq. 1 ) then
      do nDataIndex=1,KNT
        do nChannelIndex=1,KNO
          INDXTOPO = ISRCHEQI(ICHTOPO, MXTOPO, KCANO(nChannelIndex,nDataIndex))
          if ( INDXTOPO .GT. 0 ) then
            if ( MTINTRP(nDataIndex) .GE. ZCRIT(INDXTOPO) ) then
              ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**18)
              rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                   rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
              if ( B7CHCK(nChannelIndex,nDataIndex) .EQ. 0 ) then
                rejectionCodArray2(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                   rejectionCodArray2(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1                 
              end if
              if ( mwbg_debug ) then
                write(*,*)STNID(2:9),' TOPOGRAPHY REJECT.', &
                          'CHANNEL=', KCANO(nChannelIndex,nDataIndex), &
                          ' TOPO= ',MTINTRP(nDataIndex)
              end if
            end if
          end if
        end do
      end do
    end if

  end subroutine atmsTest2TopographyCheck

  !--------------------------------------------------------------------------
  !  atmsTest3UncorrectedTbCheck
  !--------------------------------------------------------------------------
  subroutine atmsTest3UncorrectedTbCheck (itest, KCANO, KNOSAT, KNO, KNT, STNID, RESETQC, &
                                          KMARQ, B7CHCK, ICHECK, rejectionCodArray, rejectionCodArray2)

    !:Purpose:                       Test 3: Uncorrected Tb check (single)
    !                                Uncorrected datum (flag bit #6 off). 
    !                                In this case switch bit 11 ON.

    ! Arguments
    integer,     intent(in)                                   :: itest(:)                 ! test number
    integer,     intent(in)                                   :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)                                   :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)                                   :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)                                   :: KNT                            ! nombre de tovs
    character *9,intent(in)                                   :: STNID                          ! identificateur du satellite
    logical,     intent(in)                                   :: RESETQC                        ! resetqc logical
    integer,     intent(inout)                                :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)                                :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(inout)                                :: B7CHCK(KNO,KNT)          ! 
    integer,     intent(inout)                                :: rejectionCodArray(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 
    integer,     intent(inout)                                :: rejectionCodArray2(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 
    ! Locals
    integer                                                   :: channelval
    integer                                                   :: nDataIndex
    integer                                                   :: nChannelIndex
    integer                                                   :: testIndex
    integer                                                   :: max 
    integer                                                   :: IBIT 

    
    if (.NOT.RESETQC) then
      testIndex = 3
      if ( itest(testIndex) .eq. 1 ) then
        do nDataIndex=1,KNT
          do nChannelIndex=1,KNO
            IBIT = AND(KMARQ(nChannelIndex,nDataIndex), 2**6)
            if ( IBIT .EQ. 0  ) then
              ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**11)
              rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                 rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
              if ( B7CHCK(nChannelIndex,nDataIndex) .EQ. 0 ) then
                rejectionCodArray2(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                    rejectionCodArray2(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1                 
              end if
              if ( mwbg_debug ) then
                write(*,*)STNID(2:9),' UNCORRECTED TB REJECT.', &
                          'CHANNEL=', KCANO(nChannelIndex,nDataIndex), &
                          ' IMARQ= ',KMARQ(nChannelIndex,nDataIndex)
              end if
            end if
          end do
        end do
      end if
    end if

  end subroutine atmsTest3UncorrectedTbCheck

  !--------------------------------------------------------------------------
  !  atmsTest4RogueCheck
  !--------------------------------------------------------------------------
  subroutine atmsTest4RogueCheck (itest, KCANO, KNOSAT, KNO, KNT, STNID, ROGUEFAC, TOVERRST, PTBOMP, &
                                  IDENTF, PMISG, MXSFCREJ, ISFCREJ, ICH2OMPREJ, MXCH2OMPREJ, & 
                                  KMARQ, B7CHCK, ICHECK, rejectionCodArray, rejectionCodArray2)

    !:Purpose:                         test 4: "Rogue check" for (O-P) Tb residuals out of range.
    !                                  (single/full).
    !                                  Also, over WATER remove CH.17-22 if CH.17 |O-P|>5K (partial) 
    !                                  Les observations, dont le residu (O-P) 
    !                                  depasse par un facteur (roguefac) l'erreur totale des TOVS.
    !                                  N.B.: a reject by any of the 3 amsua surface channels 1-3 produces the 
    !                                  rejection of ATMS sfc/tropospheric channels 1-6 and 16-17.
    !                                  OVER OPEN WATER
    !                                  ch. 17 Abs(O-P) > 5K produces rejection of all ATMS amsub channels 17-22.

    ! Arguments
    integer,     intent(in)              :: itest(:)                 ! test number
    integer,     intent(in)              :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)              :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)              :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)              :: KNT                            ! nombre de tovs
    character *9,intent(in)              :: STNID                          ! identificateur du satellite
    real,        intent(in)              :: ROGUEFAC(mwbg_maxNumChan)                  ! rogue factor 
    real,        intent(in)              :: TOVERRST(mwbg_maxNumChan,mwbg_maxNumSat)          !  erreur totale TOVs
    real,        intent(in)              :: PTBOMP(KNO,KNT)                ! radiance o-p 
    integer,     intent(in)              :: IDENTF(KNT)                    ! data flag ident  
    real,        intent(in)              :: PMISG                          ! missing value
    integer,     intent(in)              :: MXSFCREJ                       ! cst 
    integer,     intent(in)              :: ISFCREJ(MXSFCREJ)              !
    integer,     intent(in)              :: MXCH2OMPREJ                       ! cst 
    integer,     intent(in)              :: ICH2OMPREJ(MXCH2OMPREJ)              !
    integer,     intent(inout)           :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)           :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(inout)           :: B7CHCK(KNO,KNT)          ! 
    integer,     intent(inout)           :: rejectionCodArray(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 
    integer,     intent(inout)           :: rejectionCodArray2(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 

    ! Locals
    integer                              :: channelval
    integer                              :: nDataIndex
    integer                              :: nChannelIndex
    integer                              :: testIndex
    integer                              :: INDXCAN 
    real                                 :: XCHECKVAL
    logical                              :: SFCREJCT
    logical                              :: CH2OMPREJCT
    logical                              :: IBIT 


    testIndex = 4
    if ( itest(testIndex) .eq. 1 ) then
      do nDataIndex=1,KNT
        SFCREJCT = .FALSE.
        CH2OMPREJCT = .FALSE.
        do nChannelIndex=1,KNO
          channelval = KCANO(nChannelIndex,nDataIndex)
          XCHECKVAL = ROGUEFAC(channelval) * &
                     TOVERRST(channelval,KNOSAT) 
          if ( PTBOMP(nChannelIndex,nDataIndex)      .NE. PMISG    .AND. &
              ABS(PTBOMP(nChannelIndex,nDataIndex)) .GE. XCHECKVAL     ) then
            ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**16)
            rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) =  &
               rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) + 1
            if ( B7CHCK(nChannelIndex,nDataIndex) .EQ. 0 ) then
              rejectionCodArray2(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                 rejectionCodArray2(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1                 
            end if
            if ( mwbg_debug ) then
              write(*,*)STNID(2:9),'ROGUE CHECK REJECT.NO.', &
                     ' OBS = ',nDataIndex, &
                     ' CHANNEL= ',KCANO(nChannelIndex,nDataIndex), &
                     ' CHECK VALUE= ',XCHECKVAL, &
                     ' TBOMP= ',PTBOMP(nChannelIndex,nDataIndex)
            end if
            if ( channelval .EQ. 1 .OR. &
                channelval .EQ. 2 .OR. &
                channelval .EQ. 3    ) then
              SFCREJCT = .TRUE.
            end if
          end if
          if ( channelval .EQ. 17 .AND. PTBOMP(nChannelIndex,nDataIndex) .NE. PMISG .AND. &
              ABS(PTBOMP(nChannelIndex,nDataIndex)) .GT. 5.0 ) then
            CH2OMPREJCT = .TRUE.
          end if
        end do

        if ( SFCREJCT ) then
          do nChannelIndex=1,KNO
            INDXCAN = ISRCHEQI (ISFCREJ,MXSFCREJ,KCANO(nChannelIndex,nDataIndex))
            if ( INDXCAN .NE. 0 )  then
              if ( ICHECK(nChannelIndex,nDataIndex) .NE. testIndex ) then
                ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
                KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
                KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**16)
                rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                        rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
                if ( B7CHCK(nChannelIndex,nDataIndex) .EQ. 0 ) then
                  rejectionCodArray2(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                     rejectionCodArray2(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1                 
                end if
              end if
            end if
          end do
        end if

        !  amsub channels 17-22 obs are rejected if, for ch17 ABS(O-P) > 5K
        !    Apply over open water only (bit 0 ON in QC integer identf).
        !    Only apply if obs not rejected in this test already.
        IBIT = AND(IDENTF(nDataIndex), 2**0)
        if ( CH2OMPREJCT .AND. (IBIT .NE. 0) ) then
          do nChannelIndex=1,KNO
            INDXCAN = ISRCHEQI (ICH2OMPREJ,MXCH2OMPREJ,KCANO(nChannelIndex,nDataIndex))
            if ( INDXCAN .NE. 0 )  then
              if ( ICHECK(nChannelIndex,nDataIndex) .NE. testIndex ) then
                ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
                KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
                KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**16)
                rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                        rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
                if ( B7CHCK(nChannelIndex,nDataIndex) .EQ. 0 ) then
                  rejectionCodArray2(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                     rejectionCodArray2(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1                 
                end if
              end if
            end if
          end do
        end if

      end do
    end if

  end subroutine atmsTest4RogueCheck 

  !--------------------------------------------------------------------------
  !  atmsTest5ChannelSelectionUsingIutilst
  !--------------------------------------------------------------------------
  subroutine atmsTest5ChannelSelectionUsingIutilst(itest, KCANO, KNOSAT, KNO, KNT, STNID, &
                                                   IUTILST, KMARQ, ICHECK, rejectionCodArray)

    !:Purpose:                         test 5: Channel selection using array IUTILST(chan,sat)
    !                                  IUTILST = 0 (blacklisted)
    !                                         1 (assmilate) 

    ! Arguments
    integer,     intent(in)               :: itest(:)                 ! test number
    integer,     intent(in)               :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)               :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)               :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)               :: KNT                            ! nombre de tovs
    character *9,intent(in)               :: STNID                          ! identificateur du satellite
    integer,     intent(in)               :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(in)               :: IUTILST(mwbg_maxNumChan,mwbg_maxNumSat)          !  channsl selection
    integer,     intent(inout)            :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)            :: rejectionCodArray(mwbg_maxNumTest,KNO,KNT)       ! cumul of reject element 

    ! Locals
    integer                               :: channelval
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex
    integer                               :: INDXCAN 
    real                                  :: XCHECKVAL
    logical                               :: SFCREJCT
    logical                               :: CH2OMPREJCT

    testIndex = 5
    if ( itest(testIndex) .eq. 1 ) then
      do nDataIndex=1,KNT
        do nChannelIndex=1,KNO
           channelval = KCANO(nChannelIndex,nDataIndex)
           if ( IUTILST(channelval,KNOSAT) .EQ. 0 ) then
             KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**8)
             rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                   rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
             if ( mwbg_debug ) then
               write(*,*)STNID(2:9),'CHANNEL REJECT: ', &
                      ' OBS = ',nDataIndex, &
                      ' CHANNEL= ',channelval                  
             end if
           end if
        end do
      end do      
    end if

    if ( mwbg_debug ) then
      write(*,*) 'ICHECK = ',((ICHECK(nChannelIndex,nDataIndex),nChannelIndex=1,KNO),nDataIndex=1,KNT)
    end if

  end subroutine atmsTest5ChannelSelectionUsingIutilst

  !--------------------------------------------------------------------------
  ! mwbg_tovCheckAtms 
  !--------------------------------------------------------------------------
  subroutine mwbg_tovCheckAtms(TOVERRST, IUTILST, mglg_file, zlat, zlon, ilq, itt, &
                               zenith, qcflag2, qcflag1, KSAT, KORBIT, ICANO, ICANOMP, &
                               ztb, biasCorr, ZOMP, ICHECK, KNO, KNOMP, KNT, PMISG, KNOSAT, IDENT, &
                               KCHKPRF, ISCNPOS, MTINTRP, globMarq, IMARQ, rclw, riwv, rejectionCodArray, &
                               rejectionCodArray2, STNID, RESETQC)
                               


    !:Purpose:                   Effectuer le controle de qualite des radiances tovs.
    !
 
    implicit none 
    !Arguments
    integer, intent(in)              :: IUTILST(mwbg_maxNumChan,mwbg_maxNumSat) !channel Selection using array IUTILST(chan,sat)
    !                                                         IUTILST = 0 (blacklisted)
    !                                                         1 (assmilate)
    !                                                         2 (assimilate over open water only)
    real, intent(in)                 :: TOVERRST(mwbg_maxNumChan,mwbg_maxNumSat)! l'erreur totale des TOVS
    character(len=128), intent(in)   :: mglg_file
    integer, intent(in)              :: KNO                  ! nombre de canaux des observations 
    integer, intent(in)              :: KNOMP                  ! nombre de canaux des observations 
    integer, intent(in)              :: KNT                  ! nombre de tovs
    real,    intent(in)              :: zlat(:)
    real,    intent(in)              :: zlon(:)
    integer, intent(in)              :: ilq(:) 
    integer, intent(in)              :: itt(:)
    real,    intent(inout)           :: zenith(:)
    integer, intent(in)              :: qcflag2(:)
    integer, intent(in)              :: qcflag1(:,:)
    integer, intent(inout)           :: globMarq(:)
    integer, intent(in)              :: KSAT(KNT)            ! numero d'identificateur du satellite
    integer, intent(in)              :: ISCNPOS(KNT)         ! position sur le "scan"
    integer, intent(in)              :: KORBIT(KNT)          ! numero d'orbite
    integer, intent(in)              :: ICANO(:)       ! canaux des observations
    integer, intent(in)              :: ICANOMP(:)     ! canaux des residus (o-p)
    integer, intent(in)              :: KNOSAT               ! numero de satellite (i.e. indice)
    real, intent(inout)              :: ztb(:)        ! radiances
    real, intent(in)                 :: biasCorr(:)      ! correction aux radiances
    real, intent(in)                 :: zomp(:)      ! residus (o-p)
    real, intent(in)                 :: MTINTRP(KNT)         ! topographie du modele
    integer, allocatable, intent(out):: IDENT(:)           ! flag to identify all obs pts in report
   !                                                                 as being over land/ice, cloudy, bad IWV
    real, intent(in)                 :: PMISG                ! missing value
    character *9, intent(in)         :: STNID                ! identificateur du satellite
    logical, intent(in)              :: RESETQC              ! reset du controle de qualite?
    integer,allocatable, intent(out) :: ICHECK(:,:)          ! indicateur controle de qualite tovs par canal 
    !                                                                 =0, ok,
    !                                                                 >0, rejet,
    integer,allocatable, intent(out) :: KCHKPRF(:)            ! indicateur global controle de qualite tovs. Code:
    !                                                            =0, ok,
    !                                                            >0, rejet d'au moins un canal
    integer, intent(inout)           :: IMARQ(:)       ! marqueurs des radiances
    integer, intent(inout)           :: rejectionCodArray(mwbg_maxNumTest,mwbg_maxNumChan,mwbg_maxNumSat)           ! cumul du nombre de rejet par satellite, critere et par canal
    integer, intent(inout)           :: rejectionCodArray2(mwbg_maxNumTest,mwbg_maxNumChan,mwbg_maxNumSat)          ! cumul du nombre de rejet (chech n2) par satellite, critere et par canal
    real, allocatable, intent(out)   :: rclw (:)
    real, allocatable, intent(out)   :: riwv(:)

    !locals
    real                             :: PTBOMP(KNO,KNT)      ! residus (o-p)     2D
    integer                          :: KCANO(KNO,KNT)       ! canaux des observations 2D
    integer                          :: KCANOMP(KNO,KNT)     ! canaux des residus (o-p) 2D
    integer                          :: KMARQ(KNO,KNT)       ! marqueurs des radiances 2D
    integer, allocatable             :: lsq(:)
    integer, allocatable             :: trn(:) 
    logical, allocatable             :: waterobs(:)
    logical, allocatable             :: grossrej(:)
    logical                          :: reportHasMissingTb   ! true if Tb(ztb) are set to missing_value
    logical, allocatable             :: lqc(:,:)             ! dim(nt,mwbg_atmsMaxNumChan), lqc = .false. on input
    logical, allocatable             :: cloudobs(:)
    logical, allocatable             :: iwvreject(:)
    logical, allocatable             :: precipobs(:)
    real, allocatable                :: zdi(:)
    real, allocatable                :: scatec(:) 
    real, allocatable                :: scatbg(:) 
    real, allocatable                :: SeaIce(:) 

    integer, parameter               :: mwbg_maxScanAngleAMSU = 96 
    integer, parameter               :: MXSFCREJ   = 8 
    integer, parameter               :: MXCH2OMPREJ= 6  
    integer, parameter               :: MXTOPO     = 5 
    integer                          :: maxval 
    integer                          :: iRej 
    integer                          :: iNumSeaIce 
    integer                          :: JI
    integer                          :: JJ
    integer                          :: kk
    integer                          :: INDX8
    integer                          :: INDX12
    integer                          :: INO
    integer                          :: ICHN
    integer                          :: JK
    integer                          :: IBIT
    integer                          :: JC
    integer                          :: INDX
    integer                          :: INDXCAN
    integer                          :: INDXTOPO
    integer                          :: ISFCREJ(MXSFCREJ)
    integer                          :: ICH2OMPREJ(MXCH2OMPREJ)
    integer                          :: B7CHCK(KNO,KNT)
    real                             :: APPROXIM
    real                             :: XCHECKVAL
    real                             :: ROGUEFAC(mwbg_maxNumChan)
    real                             :: ZCRIT(MXTOPO)
    integer                          :: ITEST(mwbg_maxNumTest) 
    integer                          :: ICHTOPO(MXTOPO) 
    logical                          :: GROSSERROR
    logical                          :: FULLREJCT
    logical                          :: SFCREJCT
    logical                          :: CH2OMPREJCT
    logical, save                    :: LLFIRST
    integer, save                    :: numReportWithMissigTb           ! Number of BURP file reports where Tb set to mwbg_realMisg
    integer, save                    :: drycnt                          ! Number of pts flagged for AMSU-B Dryness Index             
    integer, save                    :: landcnt                         ! Number of obs pts found over land/ice 
    integer, save                    :: rejcnt                          ! Number of problem obs pts (Tb err, QCfail) 
    integer, save                    :: iwvcnt                          ! Number of pts with Mean 183 Ghz Tb < 240K 
    integer, save                    :: pcpcnt                          ! Number of scatter/precip obs
    integer, save                    :: cldcnt                          ! Number of water point covered by cloud 
    integer, save                    :: flgcnt                          ! Total number of filtered obs pts
    integer, save                    :: seaIcePointNum                  ! Number of waterobs points converted to sea ice points 
    integer, save                    :: clwMissingPointNum              ! Number of points where cloudLiquidWaterPath/SI missing 
    !                                                                     over water due bad data 

    DATA  LLFIRST / .TRUE. /

    DATA  ROGUEFAC / 2.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 4.0, &
                    4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 2.0, &
                    2.0, 4.0, 4.0, 4.0, 4.0, 4.0 /

    ! Channel sets for rejection in test 9 
    ! These LT channels are rejected if O-P fails rogue check for window ch. 1, 2, or 3
    DATA  ISFCREJ    / 1, 2, 3, 4, 5, 6, 16, 17 /
    !   These AMSU-B channels are rejected if ch. 17 O-P fails rogue check over OPEN WATER only    
    DATA  ICH2OMPREJ / 17, 18, 19, 20, 21, 22 /

    !  Data for TOPOGRAPHY CHECK
    !   Channel AMSUA-6 (atms ch 7) is rejected for topography  >  250m.
    !   Channel AMSUA-7 (atms ch 8) is rejected for topography  > 2000m.
    !   Channel AMSUB-3 (atms ch 22) is rejected for topography > 2500m.
    !                    atms ch 21  is rejected for topography > 2250m.
    !   Channel AMSUB-4 (atms ch 20) is rejected for topography > 2000m.
    DATA ICHTOPO  /     7,     8,    20,    21,    22  /
    DATA ZCRIT    /   250., 2000., 2000., 2250., 2500. /

    !  Test selection (0=skip test, 1=do test)
    !             1  2  3  4  5 
    DATA ITEST  / 1, 1, 1, 1, 1 /
       
    ! Initialisation, la premiere fois seulement!
    if (LLFIRST) then
      numReportWithMissigTb = 0
      flgcnt = 0
      landcnt = 0
      rejcnt = 0
      cldcnt = 0
      iwvcnt = 0
      pcpcnt = 0
      drycnt = 0
      seaIcePointNum = 0
      clwMissingPointNum = 0
      rejectionCodArray(:,:,:)  = 0
      rejectionCodArray2(:,:,:) = 0
      LLFIRST = .FALSE.
    end if

    ! PART 1 TESTS:

    !###############################################################################
    ! STEP 1 ) Determine which obs pts are over open water (i.e NOT near coasts or
    !          over/near land/ice) using model MG and LG fields from glbhyb2 ANAL    
    !###############################################################################
    write(*,*) ' ==> mwbg_landIceMaskAtms: '
    call mwbg_landIceMaskAtms(mglg_file, KNT, zlat, zlon, ilq, itt, &
                              lsq, trn, waterobs)

    !###############################################################################
    ! STEP 2 ) Check for values of TB that are missing or outside physical limits.    
    !###############################################################################

    write(*,*) ' ==> mwbg_grossValueCheck: '
    call mwbg_grossValueCheck(KNT,ztb,grossrej)
     
    if ( ANY(grossrej) ) then
      write(*,*) ' mwbg_grossValueCheck has detected bad Tb data. Number of affected &
           locations = ', COUNT(grossrej)
    end if

    !###############################################################################
    ! STEP 3 ) Preliminary QC checks --> set lqc(nt,mwbg_atmsMaxNumChan)=.true. 
    !          for data that fail QC     
    !###############################################################################

    write(*,*) ' ==> mwbg_firstQcCheckAtms: '
    call mwbg_firstQcCheckAtms(zenith, ilq, itt, zlat, zlon, ztb, ISCNPOS, stnid, &
                               KNO, KNT, lqc, grossrej, lsq, trn, qcflag1, &
                               qcflag2, ICANO, reportHasMissingTb)

    if ( reportHasMissingTb ) numReportWithMissigTb = numReportWithMissigTb + 1
    !  Exclude problem points from further calculations
    do kk = 1,KNT
      if ( COUNT(lqc(kk,:)) == mwbg_atmsMaxNumChan ) grossrej(kk) = .true.
    end do

    !###############################################################################
    ! STEP 4 ) mwbg_nrlFilterAtms returns rclw, scatec, scatbg and also does sea-ice 
    !          detection missing value for  rclw, scatec, scatbg  is -99.0 (e.g. over
    !          land or sea-ice).Sets trn=0 (sea ice) for points where retrieved SeaIce
    !          >=0.55. Does nothing if trn=0 (sea ice) and retrieved SeaIce<0.55.
    !###############################################################################
 
    !  
    write(*,*) ' ==>NNN mwbg_nrlFilterAtms: '
    call mwbg_nrlFilterAtms(KNT, ztb, biasCorr, zenith, zlat, lsq, trn, waterobs, &
                            grossrej, rclw, scatec, scatbg, iNumSeaIce, iRej, SeaIce)
    seaIcePointNum = seaIcePointNum + iNumSeaIce
    clwMissingPointNum = clwMissingPointNum + iRej
      
    !###############################################################################
    ! STEP 5 ) Apply NRL cloud filter, scattering index and sea-ice detection algorithms
    !          to OPEN WATER (waterobs=true) points.
    ! Points with SeaIce>0.55 are set to sea-ice points (waterobs --> false)
    !###############################################################################
    
    write(*,*) ' ==> mwbg_flagDataUsingNrlCriteria: '
    call mwbg_flagDataUsingNrlCriteria(KNT, ztb, biasCorr, rclw, scatec, scatbg, &
                                       SeaIce, grossrej, waterobs, mwbg_useUnbiasedObsForClw, &
                                       mwbg_realMisg, iwvreject, cloudobs, &
                                       precipobs, cldcnt , ident, riwv, zdi)

    !###############################################################################
    ! STEP 5 ) ! Review all the checks previously made to determine which obs are to be 
    !            accepted for assimilation and which are to be flagged for exclusion 
    !            (IMARQ). 
    !            grossrej()  = .true. if any channel had a gross error at the point
    !            cloudobs()  = .true. if CLW > clw_atms_nrl_LTrej (0.175) or precipobs
    !            precipobs() = .true. if precip. detected through NRL scattering indices
    !            waterobs()  = .true. if open water point
    !            iwvreject() = .true. if Mean 183 Ghz [ch. 18-22] Tb < 240K (too dry 
    !            for ch.20-22 over land)
    !###############################################################################

    write(*,*) ' ==> mwbg_reviewAllcriteriaforFinalFlags: '
    call mwbg_reviewAllcriteriaforFinalFlags(KNT,KNO, lqc, grossrej, waterobs, &
                                             precipobs, rclw, scatec, scatbg, iwvreject, riwv, &
                                             IMARQ, globMarq, zdi, &
                                             ident, mwbg_realMisg, drycnt, landcnt, rejcnt, &
                                             iwvcnt, pcpcnt, flgcnt)

    !###############################################################################
    ! PART 2 TESTS:
    !###############################################################################



    ! copy the original input 1D array to 2D array. The 2D arrays are used in this s/r.
    call copy1Dimto2DimRealArray(ZOMP, KNO, KNT, PTBOMP)
    call copy1Dimto2DimIntegerArray(ICANO, KNO, KNT, KCANO)
    call copy1Dimto2DimIntegerArray(ICANOMP, KNO, KNT, KCANOMP)
    call copy1Dimto2DimIntegerArray(IMARQ, KNO, KNT, KMARQ)
    ! allocations
    call mwbg_allocate1DIntegerArray(kchkprf, KNT)
    call mwbg_allocate2DIntegerArray(icheck, KNO, KNT)
    !  Initialisations
    ICHECK(:,:) = 0
    B7CHCK(:,:) = 0
   
    if ( RESETQC ) KMARQ(:,:) = 0

    !  Verification de l'integrite des donnees, c'est-a-dire que:
    !         i)  les dimensions des divers blocs de donnees concordent,
    !         ii) les listes de canaux concordent.
    if ( KNO .NE. KNOMP     ) then
      call utl_abort ('bgckMicrowave_mod:  ERROR IN DIMENSIONS OF TOVS DATA')
    end if

    do JJ=1,KNT
      do JI=1,KNO
        if ( KCANO(JI,JJ) .NE. KCANOMP(JI,JJ) ) then
          call utl_abort ('bgckMicrowave_mod: INCONSISTENT CHANNEL LISTS FOR TOVS DATA')
        end if
      end do
    end do

     
    ! 1) test 1: Check flag bit 7 on from the first bgckAtms program
    !  Includes observations flagged for cloud liquid water, scattering index,
    !  dryness index plus failure of several QC checks.
    call atmsTest1Flagbit7Check (itest, KCANO, KMARQ, KNOSAT, ICHECK, KNO, KNT, &
                                 STNID,  B7CHCK, rejectionCodArray)    
    ! 2) test 2: Topography check (partial)
    call atmsTest2TopographyCheck (itest, KCANO, KNOSAT, KNO, KNT, STNID, MTINTRP, &
                                   KMARQ, ICHTOPO, MXTOPO, ZCRIT, B7CHCK, & 
                                   ICHECK, rejectionCodArray, rejectionCodArray2) 
    ! 3) test 3: Uncorrected Tb check (single)
    !  Uncorrected datum (flag bit #6 off). In this case switch bit 11 ON.
    call atmsTest3UncorrectedTbCheck (itest, KCANO, KNOSAT, KNO, KNT, STNID, RESETQC, &
                                          KMARQ, B7CHCK, ICHECK, rejectionCodArray, rejectionCodArray2) 
    ! 4) test 4: "Rogue check" for (O-P) Tb residuals out of range. (single/full)
    !             Also, over WATER remove CH.17-22 if CH.17 |O-P|>5K (partial)
    !  Les observations, dont le residu (O-P) depasse par un facteur (roguefac) 
    !   l'erreur totale des TOVS.
    !  N.B.: a reject by any of the 3 amsua surface channels 1-3 produces the 
    !           rejection of ATMS sfc/tropospheric channels 1-6 and 16-17.
    !  OVER OPEN WATER
    !    ch. 17 Abs(O-P) > 5K produces rejection of all ATMS amsub channels 17-22.
    call atmsTest4RogueCheck (itest, KCANO, KNOSAT, KNO, KNT, STNID, ROGUEFAC, TOVERRST, PTBOMP, &
                                  IDENT, PMISG, MXSFCREJ, ISFCREJ, ICH2OMPREJ, MXCH2OMPREJ, & 
                                  KMARQ,  B7CHCK, ICHECK, rejectionCodArray, rejectionCodArray2)
 

    ! 5) test 5: Channel selection using array IUTILST(chan,sat)
    !  IUTILST = 0 (blacklisted)
    !            1 (assmilate)
    call atmsTest5ChannelSelectionUsingIutilst(itest, KCANO, KNOSAT, KNO, KNT, STNID, &
                                                   IUTILST, KMARQ, ICHECK, rejectionCodArray)

    !  Synthese de la controle de qualite au niveau de chaque point
    !  d'observation. Code:
    !            =0, aucun rejet,
    !            >0, au moins un canal rejete.
    do JJ=1,KNT
      KCHKPRF(JJ) = 0
      do JI=1,KNO 
        KCHKPRF(JJ) = MAX(KCHKPRF(JJ),ICHECK(JI,JJ))
      end do
    end do

    if ( mwbg_debug ) then
      write(*,*)'KCHKPRF = ',(KCHKPRF(JJ),JJ=1,KNT)
    end if 

    ! Copy the modified FLAG to the 1D array, used outside this s/r. 
    do JJ=1,KNT
      do JI=1,KNO
        INDX = (JJ-1)*KNO + JI 
        IMARQ(INDX) = KMARQ(JI,JJ)
      end do
    end do

    if(mwbg_debug) then
      write(*,*) ' --------------------------------------------------------------- '
      write(*,*) ' Number of BURP file reports where Tb set to mwbg_realMisg  = ', numReportWithMissigTb
      write(*,*) ' --------------------------------------------------------------- '
      write(*,*) ' 1. Number of obs pts found over land/ice           = ', landcnt
      write(*,*) ' 2. Number of problem obs pts (Tb err, QCfail)      = ', rejcnt
      write(*,*) ' 3. Number of cloudy obs  (CLW > clw_min)           = ', cldcnt
      write(*,*) ' 4. Number of scatter/precip obs                    = ', pcpcnt
      write(*,*) ' 5. Number of pts with Mean 183 Ghz Tb < 240K       = ', iwvcnt
      write(*,*) ' 6. Number of pts flagged for AMSU-B Dryness Index  = ', drycnt
      write(*,*) ' --------------------------------------------------------------- '
      write(*,*) ' Total number of filtered obs pts                   = ', flgcnt
      write(*,*) ' ----------------------------------------------------------------'
      write(*,*) ' '
      write(*,*) ' Number of waterobs points converted to sea ice points         = ', seaIcePointNum
      write(*,*) ' Number of points where CLW/SI missing over water due bad data = ', clwMissingPointNum
      write(*,*) ' --------------------------------------------------------------- '

      write(*,*) '   Meaning of IDENT flag bits: '
      write(*,*) ' '
      write(*,*) '      BIT    Meaning'
      write(*,*) '       0     off=land or sea-ice, on=open water away from coast'
      write(*,*) '       1     Mean 183 Ghz [ch. 18-22] is missing'
      write(*,*) '       2     NRL CLW is missing (over water)'
      write(*,*) '       3     NRL > clw_atms_nrl_LTrej (0.175 kg/m2) (cloudobs)'
      write(*,*) '       4     scatec/scatbg > Lower Troposphere limit 9/10 (precipobs)'
      write(*,*) '       5     Mean 183 Ghz [ch. 18-22] Tb < 240K'
      write(*,*) '       6     CLW > clw_atms_nrl_UTrej (0.200 kg/m2)'
      write(*,*) '       7     Dryness Index rejection (for ch. 22)'
      write(*,*) '       8     scatec/scatbg > Upper Troposphere limit 18/15'
      write(*,*) '       9     Dryness Index rejection (for ch. 21)'
      write(*,*) '      10     Sea ice > 0.55 detected'
      write(*,*) '      11     Gross error in Tb (any chan.) or other QC problem (all channels rejected)'
      write(*,*) ' '
      write(*,*) '   New Element 13209 in BURP file = CLW (kg/m2)'
      write(*,*) '   New Element 13208 in BURP file = ECMWF Scattering Index'
      write(*,*) '   New Element 25174 in BURP file = IDENT flag'
      write(*,*) ' '
    end if

  end subroutine mwbg_tovCheckAtms

  !--------------------------------------------------------------------------
  !  getBurpReportAdresses
  !--------------------------------------------------------------------------
  subroutine getBurpReportAdresses (burpFileNameIn, adresses)

    !:Purpose: Initial scan of file to get number of reports and number of data locations.
    !          Store address of each report in array adresses(nb_rpts) for main REPORTS loop
 
    implicit none 
    !Arguments:
    character(len=90),              intent(in)      :: burpFileNameIn        ! burp file
    integer,          allocatable,  intent(out)     :: adresses(:)   ! Location of the reports

    !integer,           intent(out)     :: nsize      ! 3* Number of the largest report in file

    !Locals:
    type(BURP_FILE)                    :: File_in  ! observation burp file
    type(BURP_RPT)                     :: rpt_in    ! observation burp file
    integer                            :: nobs_tot      ! Number of obs. in burp file
    integer                            :: nb_rpts      ! Number of report in burp file
    integer                            :: ref_rpt
    integer                            :: ibrptime
    integer                            :: irun
    integer                            :: iun_burpin
    integer                            :: nlocs 
    integer                            :: nsize 
    integer                            :: error 
    integer                            :: alloc_status

    integer, parameter                 :: mxnt = 3000
    character(len=9)                   :: DataIdType
    character(len=20)                  :: opt_missing

    !nsize = 0
   
    ! initialisation
    Call BURP_Init(File_in, iostat=error)
    Call BURP_Init(Rpt_in,  iostat=error)
    ! ouverture du fichier burp d'entree et de sortie
    Call BURP_New(File_in, FILENAME= burpFileNameIn, MODE= FILE_ACC_READ, iostat= error)
    ! Number of reports and maximum report size from input BURP file
    Call BURP_Get_Property(File_in, NRPTS=nb_rpts, IO_UNIT= iun_burpin)
    if (nb_rpts.le.1) then
      write(*,*) 'The input BURP file ', trim(burpFileNameIn), ' is empty!'
      stop
    end if
    nsize = MRFMXL(iun_burpin)
    write(*,*)
    write(*,*) 'Number of reports containing observations = ', nb_rpts-1
    write(*,*) 'Size of largest report = ', nsize
    write(*,*)

    ! Add nsize to report size to accomodate modified (larger) data blocks
    !nsize = nsize*3
    if(allocated(adresses)) deallocate(adresses)
    allocate(adresses(nb_rpts), stat=alloc_status)
  
    if (alloc_status /= 0) then
      call utl_abort ('bgckMicrowave_mod: ERROR - allocate(adresses(nb_rpts))')
    end if
  
    adresses(:) = 0
    ref_rpt = 0
    nb_rpts = 0
    nobs_tot = 0
    do
      ref_rpt = BURP_Find_Report(File_in, REPORT= Rpt_in, SEARCH_FROM= ref_rpt, iostat= error)
      if (error /= burp_noerr) call burpErrorHistory(file_in, rpt_in)
      if (ref_rpt < 0) Exit
      Call BURP_Get_Property(Rpt_in,TEMPS=ibrptime,ELEV=nlocs,STNID=DataIdType,RUNN=irun)  
      ! ELEV= the number of locations in the data box (for grouped data) ==> nt in each block
      if ( DataIdType(1:2) .eq. ">>" ) then
        write(*,*) 'Type de fichier a l_entree = ',DataIdType 
        if (DataIdType .ne. ">>DERIALT") then
          write(*,*) 'WARNING - le type de fichier devrait etre >>DERIALT'
        end if
      else if (DataIdType(1:1) .eq. "^" ) then
        if ( nlocs > mxnt ) then
          write(*,*) 'ERROR: Number of locations (nlocs) in report ',nb_rpts+1, ' exceeds limit (mxnt)!'
          write(*,*) '       nlocs = ', nlocs
          write(*,*) '       mxnt  = ', mxnt
          call  burpErrorHistory(File_in, Rpt_in)
        end if
        nobs_tot = nobs_tot + nlocs
      end if
      nb_rpts = nb_rpts+1
      adresses(nb_rpts) = ref_rpt
    
    end do
    write(*,*) ' Scan 1: Number of reports in input BURP file (nb_rpts) = ', nb_rpts
    write(*,*) '         Number of data locations (nobs_tot)             = ', nobs_tot
    ! if no reports abort
    if ( nb_rpts == 0 ) call burpErrorHistory(File_in, Rpt_in)
    ! if no observations STOP
    if ( nobs_tot == 0 ) call burpErrorHistory(File_in, Rpt_in)

  end subroutine  getBurpReportAdresses

  !--------------------------------------------------------------------------
  !  burpErrorHistory
  !--------------------------------------------------------------------------
  subroutine burpErrorHistory (burpFile, burpRpt, burpFile_opt, burpRpt_opt)

    ! :Purpose: Send and abort signal when there is a fail in reading burp file
    !          The history of all the actions in this attempt is listed
    !          Free the memory of the FIle and report instance

    implicit none
    ! Arguments 
    type(BURP_FILE),           intent(inout)                    :: burpFile  ! observation burp file
    type(BURP_RPT),            intent(inout)                    :: burpRpt   ! observation burp file
    type(BURP_FILE), optional, intent(inout)                    :: burpFile_opt  ! observation burp file
    type(BURP_RPT),  optional, intent(inout)                    :: burpRpt_opt   ! observation burp file

    write(*,*) BURP_STR_ERROR()
    write(*,*) "history"
    Call BURP_STR_ERROR_HISTORY()
    Call BURP_Free(burpFile)
    Call BURP_Free(burpRpt)
    if (present(burpFile_opt) .and. present(burpRpt_opt)) then
      Call BURP_Free(burpFile_opt)
      Call BURP_Free(burpRpt_opt)
    end if
    call utl_abort ('bgckMicrowave_mod: burpErrorHistory')
    

  end subroutine burpErrorHistory

  !--------------------------------------------------------------------------
  ! mwbg_readGeophysicFieldsAndInterpolate 
  !--------------------------------------------------------------------------
  subroutine mwbg_readGeophysicFieldsAndInterpolate(instName, zlat, zlon, MTINTRP, MGINTRP, GLINTRP)

    implicit none

    !:Purpose: Reads Modele Geophysical variables and save for the first time
    !         TOPOGRAPHIE (MF ou MX):
    !             MF est la topographie filtree avec unites en metres (filtered ME).
    !             MX est la topographie filtree avec unites en m2/s2  (geopotential topography).
    !         Glace de Mer (GL)
    !         Masque Terre-Mer (MG)
    !         Then Interpolate Those variables to observation location
    !Arguments: 
    character(*),       intent(in)   :: instName       ! Instrument Name
    real,               intent(in)   :: zlat(:)        ! Obseravtion Lats
    real,               intent(in)   :: zlon(:)        ! Observation Lons
    real, allocatable,  intent(out)  :: MGINTRP(:)     ! Glace de mer interpolees au pt d'obs.
    real, allocatable,  intent(out)  :: MTINTRP(:)     ! topographie filtree (en metres) et interpolees
    real ,allocatable,  intent(out)  :: GLINTRP(:)     ! Glace de mer interpolees au pt d'obs.
  
    ! Locals:
    real, allocatable, save  :: GL(:)                  ! Modele Glace de Mer (GL)
    real, allocatable, save  :: MG(:)                  ! Modele Masque Terre-Mer (MG)
    real, allocatable, save  :: MT(:)                  ! Modele Topographie (MT)
    real,              save  :: TOPOFACT               ! Facteur x topo pour avoir des unites en metre
    logical,           save  :: ifFirstCall = .True.   ! If .True. we read GL, MT and MG
    integer,           save  ::  gdmt                  ! topo interpolation param
    integer,           save  ::  gdmg                  ! mask terre-mer interpolation param
    integer,           save  ::  gdgl                  ! glace interpolation param
    integer                  ::  gdllsval          
    integer                  :: IUNGEO
    logical                  :: readGlaceMask 
    logical                  :: debug 
    integer                  :: ier, irec 
    integer                  :: ezqkdef, ezsetopt
    integer                  :: FSTINF,FSTPRM,FCLOS
    integer                  :: FSTLIR,FSTFRM, FNOM, FSTOUV
    integer                  :: NI, NJ, NK, IG1, IG2, IG3, IG4
    integer                  :: IDUM,IDUM1,IDUM2,IDUM3,IDUM4
    integer                  :: IDUM5,IDUM6,IDUM7,IDUM8
    integer                  :: IDUM9,IDUM10,IDUM11,IDUM12,IDUM13
    integer                  :: IDUM14,IDUM15,IDUM16,IDUM17,IDUM18
    character(len=12)        :: ETIKXX
    character(len=4)         :: CLNOMVAR
    character(len=4)         :: NOMVXX
    character(len=2)         :: TYPXX 
    character(len=1)         :: GRTYP
    character(len=128)       :: fileGeoName 
    integer                  :: NLAT
    integer                  :: NLON
    integer, PARAMETER       :: MXLON = 5
    integer, PARAMETER       :: MXLAT = 5
    integer, PARAMETER       :: MXELM = 40
    real,    PARAMETER       :: DLAT = 0.4
    real,    PARAMETER       :: DLON = 0.6
    real                     :: XLAT
    real                     :: XLON
    real, allocatable        :: ZLATBOX (:,:)
    real, allocatable        :: ZLONBOX (:,:)
    real, allocatable        :: MGINTBOX(:,:)
    real, allocatable        :: MTINTBOX(:,:)
    real, allocatable        :: GLINTBOX(:,:)
    integer                  :: dataIndex
    integer                  :: boxPointIndex
    integer                  :: latIndex
    integer                  :: lonIndex
    integer                  :: zlatNum
    integer                  :: zlonNum
    integer                  :: dataNum
    integer                  :: boxPointNum

    ! STEP 0: CHECK if ZLAT AND ZLON ARE SAME DIMENSION
    zlatNum = size(zlat)
    zlonNum = size(zlon)
    if (zlatNum .ne. zlonNum) then
      call utl_abort ('bgckMicrowave_mod: ERREUR: OBSERVATION ZLAT and ZLON should have SAME LENGTH')
    else 
      dataNum = zlatNum
    end if

    ! STEP 1: READ MT, GL and MG from the FST FILE 
    debug = mwbg_debug 
    if (instName == 'ATMS') then
      fileGeoName = './fstgzmx'
      readGlaceMask = .False.
    else if (instName == 'AMSUA') then
      fileGeoName = './fstglmg'
      readGlaceMask = .True.
    end if 
    if(ifFirstCall) then
      IUNGEO = 0 
      IER = FNOM(IUNGEO,fileGeoName,'STD+RND+R/O',0)

      ! 3) Lecture des champs geophysiques (MF/MX) du modele
      IER = FSTOUV(IUNGEO,'RND')

      ! TOPOGRAPHIE (MF ou MX).
      !     MF est la topographie filtree avec unites en metres (filtered ME).
      !     MX est la topographie filtree avec unites en m2/s2  (geopotential topography).
      TOPOFACT = 1.0
      IREC = FSTINF(IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1,' ','MF')
      CLNOMVAR = 'MF'
      if (IREC .LT. 0) then
        TOPOFACT = 9.80616
        IREC = FSTINF(IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1,' ','MX')
        CLNOMVAR = 'MX'
      end if
      if (IREC .LT. 0) then
        call utl_abort ('bgckMicrowave_mod: ERREUR: LA TOPOGRAPHIE (MF or MX) EST INEXISTANTE')
      else
        if(allocated(MT)) deallocate(MT)
        allocate ( MT(NI*NJ), STAT=ier)
        IER = FSTLIR(MT,IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1, &
           ' ',CLNOMVAR)
      end if
      
      IER = FSTPRM ( IREC, IDUM1, IDUM2, IDUM3, IDUM4, & 
          IDUM5, IDUM6, IDUM7, IDUM8, IDUM9, IDUM10,  &
          IDUM11, TYPXX, NOMVXX, ETIKXX, GRTYP, IG1, &
          IG2, IG3, IG4, IDUM12, IDUM13, IDUM14,  &
          IDUM15, IDUM16, IDUM17, IDUM18 )
       write (*,*) ' GRILLE MT : ',grtyp,ni,nj, &
                ig1,ig2,ig3,ig4
      ier  = ezsetopt('INTERP_DEGREE','LINEAR')  
      gdmt = ezqkdef(ni,nj,grtyp,ig1,ig2,ig3,ig4,iungeo)

      if (readGlaceMask) then 
        ! MG
        IREC = FSTINF(IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1,' ','MG')
        if (IREC .LT. 0) then
          call utl_abort ('bgckMicrowave_mod: ERREUR: LE MASQUE TERRE-MER EST INEXISTANT')
        end if

        if(allocated(MG)) deallocate(MG)
        allocate ( MG(NI*NJ), STAT=ier)
        IER = FSTLIR(MG,IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1,&
                 ' ','MG')

        IER = FSTPRM ( IREC, IDUM1, IDUM2, IDUM3, IDUM4, &
             IDUM5, IDUM6, IDUM7, IDUM8, IDUM9, IDUM10, &
             IDUM11, TYPXX, NOMVXX, ETIKXX, GRTYP, IG1,&
             IG2, IG3, IG4, IDUM12, IDUM13, IDUM14, &
             IDUM15, IDUM16, IDUM17, IDUM18 )
        write (*,*) ' GRILLE MG : ',grtyp,ni,nj, &
                ig1,ig2,ig3,ig4
        ier  = ezsetopt('INTERP_DEGREE','LINEAR')  
        gdmg = ezqkdef(ni,nj,grtyp,ig1,ig2,ig3,ig4,iungeo)
        ! GL
        IREC = FSTINF(IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1,' ','GL')
        if (IREC .LT. 0) then
          call utl_abort ('bgckMicrowave_mod: ERREUR: LE CHAMP GLACE DE MER EST INEXISTANT')
        end if

        if(allocated(GL)) deallocate(GL)
        allocate ( GL(NI*NJ), STAT=ier)
        IER = FSTLIR(GL,IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1, &
                 ' ','GL')

        IER = FSTPRM ( IREC, IDUM1, IDUM2, IDUM3, IDUM4, & 
             IDUM5, IDUM6, IDUM7, IDUM8, IDUM9, IDUM10, &
             IDUM11, TYPXX, NOMVXX, ETIKXX, GRTYP, IG1, &
             IG2, IG3, IG4, IDUM12, IDUM13, IDUM14, &
             IDUM15, IDUM16, IDUM17, IDUM18 )
        write (*,*) ' GRILLE GL : ',grtyp,ni,nj, &
                ig1,ig2,ig3,ig4
        ier  = ezsetopt('INTERP_DEGREE','LINEAR')  
        gdgl = ezqkdef(ni,nj,grtyp,ig1,ig2,ig3,ig4,iungeo)
      else 
        gdgl = -1
        gdmg = -1
      end if 
      IER = FSTFRM(IUNGEO)
      IER = FCLOS(IUNGEO)
      ifFirstCall = .False. 
    end if

    ! STEP 3:  Interpolation de la glace et le champ terre/mer du modele aux pts TOVS.
    ! N.B.: on examine ces champs sur une boite centree sur chaque obs.
    boxPointNum = MXLAT*MXLON
    if(allocated(ZLATBOX)) deallocate(ZLATBOX)
    allocate (ZLATBOX(boxPointNum, dataNum) , STAT=ier) 
    if(allocated(ZLONBOX)) deallocate(ZLONBOX)
    allocate (ZLONBOX(boxPointNum, dataNum) , STAT=ier) 
    if(allocated(MTINTBOX)) deallocate(MTINTBOX)
    allocate (MTINTBOX(boxPointNum, dataNum) , STAT=ier) 
    if(allocated(GLINTBOX)) deallocate(GLINTBOX)
    allocate (GLINTBOX(boxPointNum, dataNum) , STAT=ier) 
    if(allocated(MGINTBOX)) deallocate(MGINTBOX)
    allocate (MGINTBOX(boxPointNum, dataNum) , STAT=ier) 
    NLAT = (MXLAT-1)/2
    NLON = (MXLON-1)/2
    do dataIndex = 1, dataNum
      boxPointIndex = 0
      do latIndex = -NLAT, NLAT
        XLAT = ZLAT(dataIndex) +latIndex*DLAT
        XLAT = MAX(-90.0,MIN(90.0,XLAT))
        do lonIndex = -NLON, NLON
          boxPointIndex = boxPointIndex + 1
          XLON = ZLON(dataIndex) +lonIndex*DLON
          if ( XLON .LT. -180. ) XLON = XLON + 360.
          if ( XLON .GT.  180. ) XLON = XLON - 360.
          if ( XLON .lt.    0. ) XLON = XLON + 360.
           ZLATBOX(boxPointIndex,dataIndex) = XLAT
           ZLONBOX(boxPointIndex,dataIndex) = XLON
         end do
      end do
    end do
    ier = ezsetopt('INTERP_DEGREE','LINEAR')
    ier = gdllsval(gdmt,mtintbox,mt,ZLATBOX,ZLONBOX,boxPointNum*dataNum)
    if (ier .lt. 0) then
      call utl_abort ('bgckMicrowave_mod: ERROR in the interpolation of MT')
    end if
    if(readGlaceMask) then   
      ier = gdllsval(gdmg,mgintbox,mg,ZLATBOX,ZLONBOX,boxPointNum*dataNum)
      if (ier .lt. 0) then
        call utl_abort ('bgckMicrowave_mod: ERROR in the interpolation of MG')
      end if
      ier = gdllsval(gdgl,glintbox,gl,ZLATBOX,ZLONBOX,boxPointNum*dataNum)
      if (ier .lt. 0) then
        call utl_abort ('bgckMicrowave_mod: ERROR in the interpolation of GL')
      end if
    end if 

    if(allocated(MTINTRP)) deallocate(MTINTRP)
    allocate (MTINTRP(dataNum) , STAT=ier) 
    if(allocated(MGINTRP)) deallocate(MGINTRP)
    allocate (MGINTRP(dataNum) , STAT=ier) 
    if(allocated(GLINTRP)) deallocate(GLINTRP)
    allocate (GLINTRP(dataNum) , STAT=ier) 
    do dataIndex = 1, dataNum
      if (DEBUG) then
        print *, ' ------------------  '
        print *, ' dataIndex = ', dataIndex
        print *, '   '
        print *, ' zlat,zlon = ', zlat(dataIndex), zlon(dataIndex)
        print *, '   '
        print *, ' ZLATBOX = '
        print *,  (ZLATBOX(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
        print *, ' ZLONBOX = '
        print *,  (ZLONBOX(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
        print *, ' MGINTBOX = '
        print *,  (MGINTBOX(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
        print *, ' MTINTBOX = '
        print *,  (MTINTBOX(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
        print *, ' GLINTBOX = '
        print *,  (GLINTBOX(boxPointIndex,dataIndex),boxPointIndex=1,boxPointNum)
      end if
      MGINTRP(dataIndex) = 0.0
      MTINTRP(dataIndex) = 0.0
      GLINTRP(dataIndex) = 0.0
      do boxPointIndex=1,MXLAT*MXLON
        MTINTRP(dataIndex) = MAX(MTINTRP(dataIndex),MTINTBOX(boxPointIndex,dataIndex)/TOPOFACT)
        if(readGlaceMask) then      
          MGINTRP(dataIndex) = MAX(MGINTRP(dataIndex),MGINTBOX(boxPointIndex,dataIndex))
          GLINTRP(dataIndex) = MAX(GLINTRP(dataIndex),GLINTBOX(boxPointIndex,dataIndex))
        end if
      end do
      if (DEBUG) then
        print *, ' MGINTRP = ', MGINTRP(dataIndex)
        print *, ' MTINTRP = ', MTINTRP(dataIndex)
        print *, ' GLINTRP = ', GLINTRP(dataIndex)
      end if
    end do
  end subroutine mwbg_readGeophysicFieldsAndInterpolate

  !--------------------------------------------------------------------------
  !   readBurpInteger
  !--------------------------------------------------------------------------
  subroutine readBurpInteger (repIndex, burpRpt, burpBlkTypList, burpFam, burpEle, error, burpArr, &
                              burpArrName, burpLocationNum, burpChannelNum, abortIfMissing)
    !:Purpose: This subroutine takes the report, the family and the element 
    !          and read (out) the corresponding integer Array from the current block
    !          In some cases, if the array does not exist, it will be filled with MSING

    implicit none 
    ! Arguments:
    integer,                    intent(in)          :: repIndex                ! report index
    type(BURP_RPT),             intent(in)          :: burpRpt                 ! burp report
    integer,                    intent(in)          :: burpBlkTypList(:)       ! burp block type
    integer,                    intent(in)          :: burpFam                 ! burp family
    integer,                    intent(in)          :: burpEle                 ! burp Element num
    character (*),              intent(in)          :: burpArrName             ! burp array name
    logical      ,              intent(in)          :: abortIfMissing          ! abort if the array is missing
    integer, allocatable,       intent(out)         :: burpArr(:)              ! burp integer array read
    integer,                    intent(out)         :: error                   ! error status
    integer,                    intent(out)         :: burpLocationNum          ! my_nt value to be returned
    integer,                    intent(out)         :: burpChannelNum         ! my_nval value to be returned
    
    ! Locals
    type(BURP_BLOCK)                                :: burpBlk
    integer                                         :: positionIndex 
    integer                                         :: burpLocationIndex
    integer                                         :: burpChannelIndex
    integer                                         :: burpBlkTypListIndex
    integer                                         :: burpReadIndice
    integer                                         :: burpNele 
    integer                                         :: ref_blk 
    integer                                         :: burpBlkTypListNum
    integer                                         :: alloc_status 
    logical                                         :: debug

    debug = mwbg_debug

    burpBlkTypListNum = size(burpBlkTypList)
    burpBlkTypListIndex = 1
    ref_blk = 0
    do while ((ref_blk <= 0 ) .and. (burpBlkTypListIndex <= burpBlkTypListNum))
      ref_blk = 0
      write(*,*) 'ref_blk = ', ref_blk
      write(*,*) 'burpBlkTypListIndex = ', burpBlkTypListIndex
      ref_blk = BURP_Find_Block(burpRpt, &
                    BLOCK       = burpBlk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = burpFam, &
                    BTYP        = burpBlkTypList(burpBlkTypListIndex), &
                    iostat      = error)
      if (error /= burp_noerr) call utl_abort('bgckMicrowave_mod: BURP_Find_Block Error')
      burpBlkTypListIndex = burpBlkTypListIndex + 1
    end do
    if (ref_blk < 0) then
      if ( abortIfMissing ) then
        write(*,*) ' ERREUR - Elements ', burpArrName, ' BLOCK (', burpEle ,') are missing in Report = ', repIndex
        call utl_abort('bgckMicrowave_mod: Missing element  ' //trim(burpArrName))
      else
        burpArr(:) = MPC_missingValue_R4
        return
      end if
    end if
    call BURP_Get_Property(burpBlk, &
                     NELE = burpNele, &
                     NT   = burpLocationNum, & 
                     NVAL = burpChannelNum, iostat=error)
    if (error /= burp_noerr)  call utl_abort('bgckMicrowave_mod: BURP_Get_Property Error')
    alloc_status = 0
    if(allocated(burpArr)) deallocate(burpArr)
    if(debug) then
      write(*,*)' VAR = ', burpArrName
      write(*,*)' DIMS: NLOCs = ', burpLocationNum, '  NVALs = ',  burpChannelNum
    end if 
    allocate(burpArr(burpLocationNum*burpChannelNum), stat=alloc_status)
    if( alloc_status /= 0)  then
      call utl_abort('bgckMicrowave_mod: Memory allocation error  in readBurpReal Subroutine')
    end if
    burpReadIndice = BURP_Find_Element(burpBlk,burpEle,iostat=error)
    if ( burpReadIndice > 0 ) then
      positionIndex = 0
      do burpLocationIndex = 1, burpLocationNum
        do burpChannelIndex = 1, burpChannelNum
          positionIndex = positionIndex + 1
          burpArr(positionIndex) = BURP_Get_Tblval(burpBlk,burpReadIndice,burpChannelIndex,burpLocationIndex,error)
        end do
      end do
    else if ( burpReadIndice < 0 ) then  
      if ( abortIfMissing ) then
        write(*,*) ' ERREUR - Elements ', burpArrName, ' data (', burpEle ,') are missing in DATA block! Report = ', repIndex
        call utl_abort('bgckMicrowave_mod: Missing element  ' //trim(burpArrName))
      else
        burpArr(:) = MPC_missingValue_INT
      end if
    end if
  end subroutine readBurpInteger

  !--------------------------------------------------------------------------
  ! readBurpReal 
  !--------------------------------------------------------------------------
  subroutine readBurpReal (repIndex, burpRpt, burpBlkTypList, burpFam, burpEle, error, burpArr, &
                           burpArrName, burpLocationNum, burpChannelNum, abortIfMissing)
    !:Purpose: This subroutine takes the report, the family and the element 
    !          and read (out) the corresponding real Array from the current block
    !          In some cases, if the array does not exist, it will be filled with MSING

    implicit none 
    ! Arguments:
    integer,                    intent(in)          :: repIndex                ! report index
    type(BURP_RPT),             intent(in)          :: burpRpt                 ! burp report
    integer,                    intent(in)          :: burpBlkTypList(:)       ! burp block type
    integer,                    intent(in)          :: burpFam                 ! burp family
    integer,                    intent(in)          :: burpEle                 ! burp Element num
    character (*),              intent(in)          :: burpArrName             ! burp array name
    logical      ,              intent(in)          :: abortIfMissing          ! abort if the array is missing
    real, allocatable,          intent(out)         :: burpArr(:)              ! burp real array read
    integer,                    intent(out)         :: error                   ! error status
    integer,                    intent(out)         :: burpLocationNum          ! my_nt value to be returned
    integer,                    intent(out)         :: burpChannelNum         ! my_nval value to be returned
    
    ! Locals
    type(BURP_BLOCK)                                :: burpBlk
    integer                                         :: positionIndex 
    integer                                         :: burpLocationIndex
    integer                                         :: burpChannelIndex
    integer                                         :: burpBlkTypListIndex
    integer                                         :: burpReadIndice
    integer                                         :: burpNele 
    integer                                         :: ref_blk 
    integer                                         :: burpBlkTypListNum
    integer                                         :: alloc_status 
    logical                                         :: debug

    debug = mwbg_debug

    burpBlkTypListNum = size(burpBlkTypList)
    burpBlkTypListIndex = 1
    ref_blk = 0
    do while ((ref_blk <= 0 ) .and. (burpBlkTypListIndex <= burpBlkTypListNum))
      ref_blk = 0
      write(*,*) 'ref_blk = ', ref_blk
      write(*,*) 'burpBlkTypListIndex = ', burpBlkTypListIndex
      ref_blk = BURP_Find_Block(burpRpt, &
                    BLOCK       = burpBlk, &
                    SEARCH_FROM = ref_blk, &
                    BFAM        = burpFam, &
                    BTYP        = burpBlkTypList(burpBlkTypListIndex), &
                    iostat      = error)
      if (error /= burp_noerr) call utl_abort('bgckMicrowave_mod: BURP_Find_Block Error')
      burpBlkTypListIndex = burpBlkTypListIndex + 1
    end do
    if (ref_blk < 0) then
      if ( abortIfMissing ) then
        write(*,*) ' ERREUR - Elements ', burpArrName, ' BLOCK (', burpEle ,') are missing in Report = ', repIndex
        call utl_abort('bgckMicrowave_mod: Missing element  ' //trim(burpArrName))
      else
        burpArr(:) = MPC_missingValue_R4
        return
      end if
    end if
    call BURP_Get_Property(burpBlk, &
                     NELE = burpNele, &
                     NT   = burpLocationNum, & 
                     NVAL = burpChannelNum, iostat=error)
    if (error /= burp_noerr)  call utl_abort('bgckMicrowave_mod: BURP_Get_Property Error')
    alloc_status = 0
    if(allocated(burpArr)) deallocate(burpArr)
    if(debug) then
      write(*,*)' VAR = ', burpArrName
      write(*,*)' DIMS: NLOCs = ', burpLocationNum, '  NVALs = ',  burpChannelNum
    end if 
    allocate(burpArr(burpLocationNum*burpChannelNum), stat=alloc_status)
    if( alloc_status /= 0)  then
      call utl_abort('bgckMicrowave_mod: Memory allocation error  in readBurpReal Subroutine')
    end if
    burpReadIndice = BURP_Find_Element(burpBlk,burpEle,iostat=error)
    if ( burpReadIndice > 0 ) then
      positionIndex = 0
      do burpLocationIndex = 1, burpLocationNum
        do burpChannelIndex = 1, burpChannelNum
          positionIndex = positionIndex + 1
          burpArr(positionIndex) = BURP_Get_Rval(burpBlk,burpReadIndice,burpChannelIndex,burpLocationIndex,error)
        end do
      end do
    else if ( burpReadIndice < 0 ) then  
      if ( abortIfMissing ) then
        write(*,*) ' ERREUR - Elements ', burpArrName, ' data (', burpEle ,') are missing in DATA block! Report = ', repIndex
        call utl_abort('bgckMicrowave_mod: Missing element  ' //trim(burpArrName))
      else
        burpArr(:) = MPC_missingValue_R4
      end if
    end if
  end subroutine readBurpReal

  !--------------------------------------------------------------------------
  !   mwbg_getData
  !--------------------------------------------------------------------------
   subroutine mwbg_getData(burpFileNameIn, reportIndex, ISAT, zenith, ilq, itt, &
                          zlat, zlon, ztb, biasCorr, ZOMP, scanpos, nvalOut, &
                          ntOut, qcflag1, qcflag2, ican, icanomp, IMARQ, IORBIT, &
                          globMarq, resumeReport, iFlastReport, InstName, STNID)
    !:Purpose:   This routine extracts the needed data from the blocks in the report:
    !             rpt              = report

    ! NOTE:  reportIndex = report number (from MAIN program) **** do NOT MODifY ****
    !       kk = variable for loops over locations (nt)
    !        j = variable for loops over nval (nval = 1 or mwbg_atmsMaxNumChan)
    !Arguments:
    !
    character(len=90),    intent(in)     :: burpFileNameIn
    integer,              intent(in)     :: reportIndex    ! report index
    !type(BURP_RPT),       intent(in)     :: rpt           ! report
    character(*),         intent(in)     :: InstName       ! Instrument Name
    integer, allocatable, intent(out)    :: ISAT(:)        ! satellite identifier
    real   , allocatable, intent(out)    :: zenith(:)      ! satellite zenith angle (btyp=3072,ele=7024) 
    integer, allocatable, intent(out)    :: ilq(:)         ! land/sea qualifier     (btyp=3072,ele=8012)
    integer, allocatable, intent(out)    :: itt(:)         ! terrain-type (ice)     (btyp=3072,ele=13039)
    real   , allocatable, intent(out)    :: zlat(:)        ! latitude values (btyp=5120,ele=5002)
    real   , allocatable, intent(out)    :: zlon(:)        ! longitude values (btyp=5120,ele=6002)
    real   , allocatable, intent(out)    :: ztb(:)         ! brightness temperature (btyp=9248/9264,ele=12163) 
    real   , allocatable, intent(out)    :: biasCorr(:)    ! bias correction 
    real   , allocatable, intent(out)    :: ZOMP(:)        ! OMP values
    integer, allocatable, intent(out)    :: scanpos(:)     ! scan position (fov)    (btyp=3072,ele=5043)
    integer,              intent(out)    :: nvalOut        ! number of channels     (btyp=9248/9264)
    integer,              intent(out)    :: ntOut          ! number of locations    (btyp=5120,etc.)
    integer, allocatable, intent(out)    :: qcflag1(:,:)   ! flag values for btyp=3072 block ele 033078, 033079, 033080
    integer, allocatable, intent(out)    :: qcflag2(:)     ! flag values for btyp=9248 block ele 033081      
    integer, allocatable, intent(out)    :: ican(:)        ! channel numbers btyp=9248 block ele 5042 (= 1-22)
    integer, allocatable, intent(out)    :: icanomp(:)     ! omp channel numbers btyp= block ele  (= 1-22)
    integer, allocatable, intent(out)    :: IMARQ(:)       ! data flags
    integer, allocatable, intent(out)    :: IORBIT(:)      ! orbit number
    integer, allocatable, intent(out)    :: globMarq(:)    ! global Marqueur Data
    logical,              intent(out)    :: resumeReport   ! True if resume Report is read
    logical,              intent(out)    :: ifLastReport   ! True if last Report is read
    character(*),         intent(out)    :: STNID          ! Platform Name
    ! Locals
    type(BURP_FILE), save                :: File_in
    type(BURP_RPT)                       :: reportIn
    integer, allocatable                 :: qcflag1FirstColomn(:)
    integer, allocatable                 :: qcflag1SecondColomn(:)
    integer, allocatable                 :: qcflag1ThirdColomn(:)
    integer                              :: allocNtOut
    integer                              :: allocNvalOut
    integer                              :: burpLocationNum
    integer                              :: burpChannelNum
    integer                              :: eleChannel     
    integer                              :: eleDataQcFlag   
    integer                              :: eleQcFlag1(3) 
    integer                              :: idtyp
    integer                              :: nlocs 
    integer                              :: blat 
    integer                              :: blon 
    integer                              :: error 
    integer                              :: iun_burpin
    integer                              :: nblocs 
    integer                              :: handle 
    character(len=9)                     :: idStn0
    character(len=20)                    :: opt_missing
    logical, save                        :: ifFirstCall = .True.  ! If .True., will open file_in and do some work then set to .False.
    integer, allocatable, save           :: adresses(:)              ! First Loop over all reports to get their adresses and save them
    integer, save                        :: nb_rpts
    
    ! 0) Logical block to assign some attributs depending on InstName
    if (InstName == 'ATMS') then
      eleChannel       = 5042
      eleDataQcFlag    = 33081
      eleQcFlag1(1)    = 33078 
      eleQcFlag1(2)    = 33079 
      eleQcFlag1(3)    = 33080 
    else if (InstName == 'AMSUA') then
      eleChannel       = 2150
      eleDataQcFlag    = 33032
      eleQcFlag1(:)    = -1 
    else
      call utl_abort('bgckMicrowave_mod: ERREUR - Instrument Name not Recognized ')

    end if


    if (ifFirstCall) then
      write(*,*) 'mwbg_getData First Call : Initialisation'
      ! Set BURP "missing value" for reals
  ! Set BURP "missing value" for reals
      opt_missing = 'MISSING'
      Call BURP_Set_Options(real_OPTNAME=opt_missing,real_OPTNAME_VALUE=mwbg_realMisg)

      !opt_missing = 'MISSING'
      !Call BURP_Set_Options(real_OPTNAME=opt_missing,real_OPTNAME_VALUE=MPC_missingValue_R4)
      ! LOOP OVER ALL REPORTS OF THE INPUT FILE, APPLY PROCESSING, AND write TO OUTPUT FILE.
      call getBurpReportAdresses (burpFileNameIn, adresses)
      ! initialisation
      Call BURP_Init(File_in, iostat=error)
      ! ouverture du fichier burp d'entree et de sortie
      Call BURP_Init(reportIn,  iostat=error)
      ! Number of reports and maximum report size from input BURP file
      Call BURP_New(File_in,  FILENAME = burpFileNameIn,  MODE= FILE_ACC_READ,   iostat= error)
      Call BURP_Get_Property(File_in, NRPTS=nb_rpts, IO_UNIT= iun_burpin)
      write(*,*)
      write(*,*) 'Number of reports containing observations = ', nb_rpts-1
      write(*,*)
      ifFirstCall = .False.
    end if
    if (nb_rpts.lt.1) then
      write(*,*) 'The input BURP file ', burpFileNameIn, ' is empty!'
      call utl_abort('bgckMicrowave_mod: Empty Input File ')
    end if    
    ! last report check
    ifLastReport = reportIndex == nb_rpts
    
    Call BURP_Get_Report(File_in, REPORT= reportIn, REF= adresses(reportIndex), iostat= error) 
    if (error /= burp_noerr) call burpErrorHistory(file_in, reportIn) 
    Call BURP_Get_Property(reportIn,STNID=idStn0,IDTYP=idtyp,ELEV=nlocs,LATI=blat,LONG=blon,NBLK=nblocs,HANDLE=handle)
    if ( idStn0(1:2) .eq. ">>" ) then
      resumeReport = .True.
      return 
    else 
      resumeReport = .False.
      STNID = idStn0 
    end if
    
    !  Get OMP data from the DATA block     BTYP =  9322 or 9226 or 9258 or 9274 and bfma = 14
    call readBurpReal (reportIndex, reportIn, (/9322,9226,9258,9274/), 14, 12163, error, ZOMP, 'Omp_Data', &
                   burpLocationNum, burpChannelNum, abortIfMissing = .FALSE.) 
    
    if ( ALL(ZOMP(:) == MPC_missingValue_R4 )) then
      return
    end if

    call readBurpInteger (reportIndex, reportIn, (/9322,9226,9258,9274/), 14, eleChannel, error, ICANOMP, 'OMP_Channels', &
                          burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 

    !  Get the lat,lon from time/location block    BTYP = 5120  (also get nt)
    call readBurpReal(reportIndex, reportIn, (/5120/), 0, 5002, error, zlat, 'LAT', &
                      burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 

    call readBurpReal(reportIndex, reportIn, (/5120/), 0, 6002, error, zlon, 'LON', &
                      burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 

    !  Get info elements from the INFO block   BTYP = 3072
    call readBurpInteger(reportIndex, reportIn, (/3072/), 0, 1007, error, ISAT, 'Sat_Identifier', &
                         burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
    call readBurpInteger(reportIndex, reportIn, (/3072/), 0, 5040, error, IORBIT, 'Orbit_Number', &
                         burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
    call readBurpReal(reportIndex, reportIn, (/3072/), 0, 7024, error, ZENITH, 'Zenith_Angle', &
                      burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
    call readBurpInteger(reportIndex, reportIn, (/3072/), 0, 8012, error, ILQ, 'LandSea_Qualifier', &
                         burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
    call readBurpInteger(reportIndex, reportIn, (/3072/), 0, 13039, error, ITT, 'Terrain_Type', &
                         burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
    call readBurpInteger(reportIndex, reportIn, (/3072/), 0, 5043, error, SCANPOS, 'Scan_Position', &
                         burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 

    !  Get info elements 33078 33079 33080 if needed
    if (ALL(eleQcFlag1 /= -1)) then
      call readBurpInteger(reportIndex, reportIn, (/3072/), 0, eleQcFlag1(1), error, qcflag1FirstColomn, 'Geoloc_Quality_QcFlag', &
                           burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
      call readBurpInteger(reportIndex, reportIn, (/3072/), 0, eleQcFlag1(2), error, qcflag1SecondColomn, 'Granule_Level_QcFlag', &
                           burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
      call readBurpInteger(reportIndex, reportIn, (/3072/), 0, eleQcFlag1(3), error, qcflag1ThirdColomn, 'Scan_Level_QcFlag', &
                           burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
      if(allocated(qcflag1)) deallocate(qcflag1)
      allocate(qcflag1(size(qcflag1FirstColomn),3))
      qcflag1(:,1) = qcflag1FirstColomn
      qcflag1(:,2) = qcflag1SecondColomn
      qcflag1(:,3) = qcflag1ThirdColomn
    end if

    !  Get data from the DATA block     BTYP = 9248 or 9264    (also get nval = mwbg_atmsMaxNumChan)
    call readBurpReal (reportIndex, reportIn, (/9248,9264/), 0, 12163, error, ztb, 'Tb_data', &
                       burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
    nvalOut = burpChannelNum    ! set nvalOut (#channels) for MAIN program
    ntOut = burpLocationNum     ! set ntOut (#locations) for MAIN program
    call readBurpReal (reportIndex, reportIn, (/9248,9264/), 0, 12233, error, biasCorr, 'Bias_Corr_data', &
                       burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
    call readBurpInteger (reportIndex, reportIn, (/9248,9264/), 0, eleChannel, error, ICAN, 'Channel_Numbers', &
                          burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 
    call readBurpInteger (reportIndex, reportIn, (/9248,9264/), 0, eleDataQcFlag, error, qcflag2, 'Data_level_Qc_Flag', &
                          burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.) 

    !  Bloc marqueurs multi niveaux de radiances: bloc 15362, 15392, 15408.
    call readBurpInteger (reportIndex, reportIn, (/15362,15392,15408/), 0, 212163, error, IMARQ, 'IMARQ', &
                          burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.)

    !  Bloc info 3d: bloc 5120.
    call readBurpInteger (reportIndex, reportIn, (/5120/), 0, 55200, error, globMarq, 'Marqueurs_Gobaux', &
                          burpLocationNum, burpChannelNum, abortIfMissing = .TRUE.)

  end subroutine mwbg_getData

  !--------------------------------------------------------------------------
  !   mwbg_findSatelliteIndex
  !--------------------------------------------------------------------------
  subroutine mwbg_findSatelliteIndex(STNID, satelliteId, INOSAT)
    !:Purpose: Find the STNID Index in the satelliteId array

    ! Arguments:
    character(len=9), intent(in)  :: STNID
    character(len=9), intent(in)  :: satelliteId(:)
    integer,          intent(out) :: INOSAT 
    ! Locals
    integer                       :: satIndex   
    integer                       :: numSats 
    ! trouver l'indice du satellite
    numSats = size(satelliteId)
    INOSAT = 0
    do satIndex = 1, numSats
      if ( trim(STNID) .EQ. '^'//satelliteId(satIndex) ) then
        INOSAT = satIndex
        write(*,*)' SATELLITE = ', STNID
        write(*,*)'    INOSAT = ', INOSAT
      end if
    end do 
    if ( INOSAT .EQ. 0 ) then
      call utl_abort('bgckMicrowave_mod:  IN STATS FILE No satellite '//trim(STNID))
    end if

  end subroutine mwbg_findSatelliteIndex

  !--------------------------------------------------------------------------
  !  mwbg_landIceMaskAtms
  !--------------------------------------------------------------------------
  subroutine mwbg_landIceMaskAtms(mglg_file,npts,zlat,zlon,ilq,itt, zlq,ztt,waterobs)
    ! Adapted from: land_ice_mask_ssmis.ftn90 of mwbg_ssmis (D. Anselmo, S. Macpherson)
    !
    ! Object:   This routine sets waterobs array by performing a land/ice proximity check using
    !           using analysis MG and LG (or GL) fields used by the model which produces the trial field.
    !           The purpose of this check is to remove obs that reside close to coasts or ice,
    !           and so whose TBs may be contaminated.
    !           The GEM Global (glbhyb2) analysis contains MG and LG fields (on different grids).
    !
    !           NOTE: The 0.1 deg binary ice field check from land_ice_mask_ssmis.ftn90
    !           was removed. The land/sea qualifier (zlq) and terrain type (ztt) are modified
    !           to indicate proximity to land and sea-ice but are NOT changed in output BURP file.
    !
    !           In the application of this check, a 5x5 mesh, with spacing defined by rlat_km and
    !           rlon_km, is positioned with its center over an obs pt (2 grid pts on either side
    !           of the obs pt; size of mesh is equal to 4*rlat_km x 4*rlon_km). The values of MG
    !           and LG are evaluated at the grid points of this mesh. The maximum value of each
    !           determines whether the obs pt is too close to ice or land to be retained.
    !           **NOTE: the threshold value for MG has a very strong effect on the distance
    !                   from land that is permitted for an obs to be retained
    !
    !
    !      Maximum FOV             x---x---x---x---x     ^
    !         = 75km x 75km        |   |   |   |   |     |
    !         for Meso-sphere CHs  x---x---x---x---x     |
    !         = 74km x 47km        |   |   |   |   |     |
    !         for 19 GHz           x---x---o---x---x     | = 4*rlat_km
    !                              |   |   |   |   |     | = 4*40 km
    !                           ^  x---x---x---x---x     | = 160 km = 80 km north & south
    !                   rlat_km |  |   |   |   |   |     |
    !                           v  x---x---x---x---x     v
    !                                          <--->
    !                                         rlon_km
    !
    !                              <--------------->
    !                                 = 4*rlon_km
    !                                 = 4*40 km
    !                                 = 160 km = 80 km east & west
    !
    !
    !               MG value = 1.0  ==>  LAND       MG value = 0.0  ==>  OCEAN
    !               LG value = 1.0  ==>  ICE        LG value = 0.0  ==>  NO ICE
    !
    !
    ! Version:      Date:      Comment:
    ! --------      -----      --------
    !   0.1       16/08/12     Original adapted code.      S. Macpherson  
    !   0.2       01/03/14     Open mglg_file in R/O mode  S. Macpherson
    !
    !--------------------------------------------------------------------
    !  Variable Definitions
    !  --------------------
    ! - mglg_file  : input  -  name of file holding model MG and LG (or GL) fields
    ! - npts       : input  -  number of input obs pts in report
    ! - zlat       : input  -  array holding lat values for all obs pts in report
    ! - zlon       : input  -  array holding lon values for all obs pts in report
    ! - zlq        : in/out -  array holding land/sea qualifier values for all obs
    !                        pts of report (0 = land, 1 = sea)
    ! - ztt        : in/out -  array holding terrain-type values for all obs pts
    !                        of current report (-1 land/open water, 0 = ice)
    ! - waterobs   : output -  logical array identifying for each obs in current report
    !                        whether it is over open water, far from coast/ice
    ! - mxlat      : internal-  number of grid pts in lat. direction for mesh
    ! - mxlon      : internal-  number of grid pts in lon. direction for mesh
    ! - rlat_km    : internal-  spacing desired between mesh grid points in km
    !                        along lat. direction
    ! - rlon_km    : internal-  spacing desired between mesh grid points in km
    !                        along lon. direction
    ! - dlat       : internal-  spacing between mesh grid points along lon. direction
    !                        in degrees computed from rlat_km
    ! - dlon       : internal-  spacing between mesh grid points along lon. direction
    !                        in degrees computed from rlon_km
    ! - rkm_per_deg : internal- distance in km per degree
    !                           = Earth radius * PI/180.0
    !                           = 6371.01 km * PI/180.0
    !                           = 111.195 km
    ! - nlat,nlon  : internal-  used to define the lat/lon of the grid pts of mesh
    ! - zlatbox    : internal-  lat values at all grid pts of mesh for all obs pts
    ! - zlonbox    : internal-  lon values at all grid pts of mesh for all obs pts
    ! - latmesh    : internal-  lat values at all grid pts of mesh for 1 obs pt
    ! - lonmesh    : internal-  lon values at all grid pts of mesh for 1 obs pt
    ! - mgintob    : internal-  interpolated MG values at all grid pts of mesh for 1 obs pt
    ! - lgintob    : internal-  interpolated LG values at all grid pts of mesh for 1 obs pt
    ! - mgintrp    : internal-  max. interpolated MG value on mesh for all obs pts
    ! - lgintrp    : internal-  max. interpolated LG value on mesh for all obs pts
    ! - MGthresh   : internal-  maximum allowable land fraction for obs to be kept
    ! - LGthresh   : internal-  maximum allowable ice  fraction for obs to be kept
    implicit none

    ! Arguments:
    character(len=128), intent(in) :: mglg_file

    integer, intent(in)                   :: npts
    real,    intent(in)                   :: zlat(:)
    real,    intent(in)                   :: zlon(:)
    integer, intent(in)                   :: ilq(:) 
    integer, intent(in)                   :: itt(:)
    integer, intent(out), allocatable     :: zlq(:)
    integer, intent(out), allocatable     :: ztt(:)
    logical, intent(out), allocatable     :: waterobs(:)

    ! Locals:
    integer, parameter :: mxlat=5,mxlon=5
    integer, parameter :: iungeo=50

    integer :: ier,key,istat
    integer :: ni,nj,nk,nilg,njlg
    integer :: ig1,ig2,ig3,ig4,ig1lg,ig2lg,ig3lg,ig4lg
    integer :: idum4,idum5,idum6,idum7,idum8,idum9,idum10,idum11
    integer :: idum12,idum13,idum14,idum15,idum16,idum17,idum18

    integer :: indx,ii,jj,kk
    integer :: nlat,nlon

    integer, dimension(15) :: alloc_status = 0
  
    real, parameter :: pi=3.141592654
    real, parameter :: MGthresh=0.01,LGthresh=0.01
    real, parameter :: rlat_km=40.0,rlon_km=40.0
    real, parameter :: rkm_per_deg=111.195

    real :: xlat,xlatrad,xlon,rii,rjj
    real :: dlat,dlon

    character(len=12) :: etikxx
    character(len=4)  :: nomvxx
    character(len=2)  :: typxx
    character(len=1)  :: grtyp,grtyplg
  
    logical  :: llg

    ! F90 allocatable arrays:
    real, allocatable, dimension(:)   :: mg,lg
    real, allocatable, dimension(:)   :: latmesh,lonmesh
    real, allocatable, dimension(:)   :: mgintob,lgintob
    real, allocatable, dimension(:,:) :: zlatbox,zlonbox
    real, allocatable, dimension(:)   :: mgintrp,lgintrp
  
    ! RMNLIB interpolating functions:
    integer :: ezsetopt,ezqkdef
    integer :: gdllsval,gdid,gdidlg

    ! Define FORTRAN FST functions:
    integer, external :: fstinf,fstprm,fstlir
    integer, external :: fstouv,fstfrm,fstinl,fstvoi

    integer :: idum1,idum2,idum3

    ! Allocate space for arrays holding values on mesh grid pts.
    call mwbg_allocate1DRealArray(latmesh, mxlat*mxlon)
    call mwbg_allocate1DRealArray(lonmesh, mxlat*mxlon)
    call mwbg_allocate1DRealArray(mgintob, mxlat*mxlon)
    call mwbg_allocate1DRealArray(lgintob, mxlat*mxlon)
    call mwbg_allocate2DRealArray(zlatbox, mxlat*mxlon, npts)
    call mwbg_allocate2DRealArray(zlonbox, mxlat*mxlon, npts)
    call mwbg_allocate1DIntegerArray(zlq, npts)
    call mwbg_allocate1DIntegerArray(ztt, npts)
    call mwbg_allocate1DLogicalArray(waterobs, npts)
    zlq(:) = ilq(1:npts)  ! land/sea qualifier
    ztt(:) = itt(1:npts)  ! terrain type (sea-ice)
    
    ! Open FST file.
    ier = fnom( iungeo,mglg_file,'STD+RND+R/O',0 )
    ier = fstouv( iungeo,'RND' )

    ! Read MG field.
    key = fstinf(iungeo,ni,nj,nk,-1,' ',-1,-1,-1,' ' ,'MG')
    if ( key <  0 ) then
      call utl_abort('bgckMicrowave_mod:  mwbg_landIceMaskAtms The MG field is MISSING')
    end if

    call mwbg_allocate1DRealArray(mg, ni*nj)

    ier = fstlir(mg,iungeo,ni,nj,nk,-1,' ',-1,-1,-1,' ','MG')

    ier = fstprm(key,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,    &
                idum9,idum10,idum11,typxx,nomvxx,etikxx,grtyp,ig1,ig2,  &
                ig3,ig4,idum12,idum13,idum14,idum15,idum16,idum17,      &
                idum18)


    ! Read LG field. Use GL field as backup.
    ! **CAUTION**: Discontinuities in GL field may cause interpolation problems! LG field is preferable.
    llg=.false.
    key = fstinf(iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ' ,'LG')
    if ( key <  0 ) then
      key = fstinf(iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ' ,'GL')
      if ( key <  0 ) then
        call utl_abort('bgckMicrowave_mod:  mwbg_landIceMaskAtms The LG or GL field is MISSING')
      else
        !write(*,*) 'mwbg_landIceMaskAtms: The GL field was found and will be used.'
      end if
    else
      llg=.true.
    end if

    call mwbg_allocate1DRealArray(lg, nilg*njlg)

    if ( llg ) then
      ier = fstlir(lg,iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ','LG')
    else
      ier = fstlir(lg,iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ','GL')
    end if

    ier = fstprm(key,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,          &
                idum9,idum10,idum11,typxx,nomvxx,etikxx,grtyplg,ig1lg,ig2lg,  &
                ig3lg,ig4lg,idum12,idum13,idum14,idum15,idum16,idum17,        &
                idum18)

    ! For each obs pt, define a grid of artificial pts surrounding it.
    nlat = ( mxlat - 1 ) / 2
    nlon = ( mxlon - 1 ) / 2

    dlat = rlat_km / rkm_per_deg
    do kk = 1, npts
      indx = 0

      do ii = -nlat, nlat
        rii = float(ii)
        xlat = zlat(kk) + rii*dlat
        xlat = max( -90.0, min(90.0,xlat) )
        xlatrad = xlat*pi/180.0

        do jj = -nlon, nlon
          dlon = rlon_km / ( rkm_per_deg*cos(xlatrad) )
          rjj = float(jj)
          indx = indx + 1
          xlon = zlon(kk) + rjj*dlon
          if ( xlon < -180. ) xlon = xlon + 360.
          if ( xlon >  180. ) xlon = xlon - 360.
          if ( xlon <    0. ) xlon = xlon + 360.
          zlatbox(indx,kk) = xlat
          zlonbox(indx,kk) = xlon
        end do

      end do
    end do


    ! Interpolate values from MG and LG field to grid pts of mesh centred over each obs pt.
    ! Determine for each obs pt, the max interpolated MG and LG value within the box
    ! surrounding it.
    ier    = ezsetopt('INTERP_DEGREE','LINEAR')
    gdid   = ezqkdef(ni,nj,grtyp,ig1,ig2,ig3,ig4,iungeo)
    gdidlg = ezqkdef(nilg,njlg,grtyplg,ig1lg,ig2lg,ig3lg,ig4lg,iungeo)
    
    call mwbg_allocate1DRealArray(mgintrp, npts)
    call mwbg_allocate1DRealArray(lgintrp, npts)

    mgintrp(:) = 0.0
    lgintrp(:) = 0.0
    do kk = 1, npts

      latmesh = zlatbox(:,kk)
      lonmesh = zlonbox(:,kk)

      ier  = gdllsval(gdid,mgintob,mg,latmesh,lonmesh,mxlat*mxlon)
      ier  = gdllsval(gdidlg,lgintob,lg,latmesh,lonmesh,mxlat*mxlon)

      mgintrp(kk) = maxval(mgintob(:))
      lgintrp(kk) = maxval(lgintob(:))

    end do

    !  Initialize all obs as being over land and free of ice or snow.
    !  Determine which obs are over open water.
    waterobs(:) = .false.   ! not over open water
    ztt(:) = -1             ! no ice (reset terain type)
    zlq(:) = 0              ! land   (reset land/sea qualifier)

    do kk = 1, npts
      if ( mgintrp(kk) < MGthresh ) zlq(kk) = 1  ! ocean point away from coast
      if ( lgintrp(kk) >= LGthresh .and. zlq(kk) == 1 ) ztt(kk) = 0  ! sea-ice affected point
      if ( lgintrp(kk)  < LGthresh .and. zlq(kk) == 1 ) then
        waterobs(kk) = .true.  ! water point not in close proximity to land or sea-ice
      end if
    end do

    ier = fstfrm(iungeo)
    ier = fclos(iungeo)

  end subroutine mwbg_landIceMaskAtms

  !--------------------------------------------------------------------------
  ! mwbg_grossValueCheck  
  !--------------------------------------------------------------------------
  subroutine mwbg_grossValueCheck(npts,ztb,grossrej)

    !:Purpose: Check Tbs for values that are missing or outside physical limits.
    !          **NOTE: REJECT ALL CHANNELS OF ONE IS FOUND TO BE BAD.
    implicit none

    ! Arguments
    integer, intent(in)               :: npts             ! number of obs pts to process
    real,    intent(in)               :: ztb(:)           ! bs from input BURP file
    logical, intent(out), allocatable :: grossrej(:)      ! ogical array defining which obs are to be rejected

    ! Locals
    integer :: ii, indx1, indx2

    call mwbg_allocate1DLogicalArray(grossrej, npts)
    
    grossrej(1:npts) = .true.
    indx1 = 1
    do ii = 1, npts

      indx2 = ii*mwbg_atmsMaxNumChan
      if ( all( ztb(indx1:indx2) > 50.0 ) .and. all( ztb(indx1:indx2) < 380.0 ) ) then
        grossrej(ii) = .false.
      end if
      indx1 = indx2 + 1

    end do

  end subroutine mwbg_grossValueCheck

  !--------------------------------------------------------------------------
  !  mwbg_firstQcCheckAtms
  !--------------------------------------------------------------------------
  subroutine mwbg_firstQcCheckAtms(zenith, ilq, itt, zlat, zlon, ztb, scanpos, stnid,&
                                   nval, nt, lqc, grossrej, lsq, trn, qcflag1, qcflag2, &
                                   ican, reportHasMissingTb)
    !  This routine performs basic quality control checks on the data. It sets array
    !  lqc(nt,mwbg_atmsMaxNumChan) elements to .true. to flag data with failed checks.
    !
    !  The 7 QC checks are:
    !                 - 1) Invalid land/sea qualifier or terrain type,
    !                 - 2) Invalid field of view number,
    !                 - 3) Satellite zenith angle missing or out of range, (> 75 deg),
    !                 - 4) lat,lon check (lat,lon = O(-90.), 0(-180.))
    !                 - 5) Change in (computed) lsq,trn from (input) ilq,itt (from MG,LG fields)
    !                      ilq= 0,1 (from hi-res land/sea mask interpolated to obs point [CMDA])
    !                      itt=-1,0 (from hi-res ice analysis  interpolated to obs point [CMDA])
    !                      lsq= 0,1 (from max interp MG (0.0 to 1.0) in box surrounding obs point)
    !                      trn=-1,0 (from max interp LG (0.0 to 1.0) in box surrounding obs point)
    !                 - 6) ATMS quality flag check (qual. flag elements 33078,33079,33080,33081)

    !
    !  In most cases, lqc(ii,mwbg_atmsMaxNumChan) is set to .true. for all channels at point ii
    !  if the check detects a problem. In addition, Tb (ztb) is set to missing_value 
    !  for checks 3 and 4 fails.
    implicit none

    ! Arguments
    integer,              intent(in)                :: ilq(:)
    integer,              intent(in)                :: itt(:)
    integer,              intent(in)                :: scanpos(:)
    integer,              intent(in)                :: ican(:)
    character(len=9),     intent(in)                :: stnid
    integer,              intent(in)                :: qcflag2(:)
    integer,              intent(in)                :: qcflag1(:,:)
    integer,              intent(in)                :: nt
    integer,              intent(in)                :: nval
    integer,              intent(in)                :: lsq(:)
    integer,              intent(in)                :: trn(:)
    logical,              intent(in)                :: grossrej(:)     ! dim(nt), true if 1 or more Tb fail gross error check
    real,                 intent(in)                :: zlat(:)
    real,                 intent(in)                :: zlon(:)
    real,                 intent(inout)             :: ztb(:)
    real,                 intent(inout)             :: zenith(:)
    logical,              intent(out)               :: reportHasMissingTb ! true if Tb(ztb) are set to missing_value
    logical, allocatable, intent(out)               :: lqc(:,:)        ! dim(nt,mwbg_atmsMaxNumChan), lqc = .false. on input

    ! Locals
    integer :: ii, jj, indx1, icount
    logical :: fail, fail1, fail2

    reportHasMissingTb = .false.
    call mwbg_allocate2dLogicalArray(lqc, nt, nval)
    lqc(:,:) = .false.  ! Flag for preliminary QC checks
    ! Global rejection checks

    ! Check if number of channels is correct
    if ( nval /= mwbg_atmsMaxNumChan ) then
      write(*,*) 'WARNING: Number of channels (',nval, ') is not equal to mwbg_atmsMaxNumChan (', mwbg_atmsMaxNumChan,')'
      write(*,*) '         All data flagged as bad and returning to calling routine!'
      lqc(:,:) = .true.  ! flag all data in report as bad
      return
    end if

    ! Check for errors in channel numbers (should be 1-22 for each location ii)
    indx1 = 1
    fail = .false.
    do ii = 1,nt
      do jj = 1,mwbg_atmsMaxNumChan
        if ( ican(indx1+jj-1) /= jj ) fail = .true.
      end do
      indx1 = indx1 + mwbg_atmsMaxNumChan
    end do
    if ( fail ) then
      write(*,*) 'WARNING: Bad channel number(s) detected!'
      write(*,*) '         All data flagged as bad and returning to calling routine!'
      write(*,*) '  ican(nt*mwbg_atmsMaxNumChan) array = ', ican(:)
      lqc(:,:) = .true.  ! flag all data in report as bad
      return
    end if

    ! 1) invalid land/sea qualifier or terrain type
    !  ilq = 0 (land),     1 (sea)
    !  itt = 0 (sea-ice), -1 otherwise
    !  lsq = 1 (sea, away from land/coast [MG]),      0 otherwise
    !  trn = 0 (over or near analyzed sea-ice [LG]), -1 otherwise
    do ii = 1,nt
      fail = .false.
      if ( ilq(ii) < 0  .or. ilq(ii) > 2 ) fail = .true.
      if ( itt(ii) < -1 .or. itt(ii) > 1 ) fail = .true.
      if ( fail ) then
        write(*,*) 'WARNING: Invalid land/sea qualifier or terrain type!'
        write(*,*) '  ilq, itt, (lat, lon) = ', ilq(ii), itt(ii), '(',zlat(ii), zlon(ii),')'
      end if
      if ( ilq(ii) == 0 .and. itt(ii) == 0 ) then
        fail = .true.
        write(*,*) 'WARNING: Sea ice point (itt=0) at land point (ilq=0)!'
        write(*,*) ' lat, lon =  ', zlat(ii), zlon(ii)
      end if
      if ( fail ) lqc(ii,:) = .true.
    end do

    do ii = 1,nt
      fail = .false.
      if ( lsq(ii) < 0  .or. lsq(ii) > 2 ) fail = .true.
      if ( trn(ii) < -1 .or. trn(ii) > 1 ) fail = .true.
      if ( fail ) then
        write(*,*) 'WARNING: Invalid model-based (MG/LG) land/sea qualifier or terrain type!'
        write(*,*) '  lsq, trn, (lat, lon) = ', lsq(ii), trn(ii), '(',zlat(ii), zlon(ii),')'
      end if
      if ( fail ) lqc(ii,:) = .true.
    end do
 
    ! 2) invalid field of view number
    do ii = 1,nt
      fail = .false.
      if ( scanpos(ii) < 1  .or. scanpos(ii) > mwbg_maxScanAngle ) then
        fail = .true.
        write(*,*) 'WARNING: Invalid field of view! scanpos, lat, lon = ', scanpos(ii), zlat(ii), zlon(ii)
      end if
      if ( fail ) lqc(ii,:) = .true.
    end do

    ! 3) satellite zenith angle missing or out of range (> 75 deg)
    !  If bad zenith, then set Tb (and zenith) = missing value
    indx1 = 1
    do ii = 1,nt
      fail = .false.
      if ( zenith(ii) > 75.0 .or. zenith(ii) < 0. ) then
        fail = .true.
        write(*,*) 'WARNING: Bad or missing zenith angle! zenith, lat, lon = ', zenith(ii), zlat(ii), zlon(ii)
        zenith(ii) = mwbg_realMisg
        reportHasMissingTb = .true.
      end if
      do jj = 1,mwbg_atmsMaxNumChan
        if ( fail ) then
          lqc(ii,jj) = .true.
          ztb(indx1+jj-1) = mwbg_realMisg
        end if
      end do
      indx1 = indx1 + mwbg_atmsMaxNumChan
    end do

    ! 4) Lat,lon check
    ! Check for undecoded BURP file integer values of lat,lon = 0,0
    ! (usually associated with missing zenith angle and erroneous Tb=330K)

    icount = 0
    indx1 = 1
    do ii = 1,nt
      fail = .false.
      if ( zlat(ii) == -90.0  .and. zlon(ii) == -180.0 ) then
        fail = .true.
        icount =  icount + 1
        reportHasMissingTb = .true.
      end if
      do jj = 1,mwbg_atmsMaxNumChan
        if ( fail ) then
          lqc(ii,jj) = .true.
          ztb(indx1+jj-1) = mwbg_realMisg
        end if
      end do
      indx1 = indx1 + mwbg_atmsMaxNumChan
    end do
    if ( icount > 0 ) write(*,*) 'WARNING: Bad lat,lon pair(s) detected. Number of locations = ', icount

    icount = 0
    indx1 = 1
    do ii = 1,nt
      fail = .false.
      if ( abs(zlat(ii)) > 90.0  .or. abs(zlon(ii)) > 180.0 ) then
        fail = .true.
        icount =  icount + 1
        reportHasMissingTb = .true.
      end if
      do jj = 1,mwbg_atmsMaxNumChan
        if ( fail ) then
          lqc(ii,jj) = .true.
          ztb(indx1+jj-1) = mwbg_realMisg
        end if
      end do
      indx1 = indx1 + mwbg_atmsMaxNumChan
    end do
    if ( icount > 0 ) write(*,*) 'WARNING: Lat or lon out of range! Number of locations = ', icount

    !  5) Change in land/sea qualifier or terrain-type based on MG,LG fields
    icount = 0
    do ii = 1,nt
      fail = .false.
      if ( (ilq(ii) /= lsq(ii)) .or. (itt(ii) /= trn(ii)) ) fail = .true.
      if ( fail ) then
        icount =  icount + 1
      end if
    end do
    if ( icount > 0 ) write(*,*) 'INFO: Num. pts with land/sea qualifier or terrain type changed (MG,LG) = ', icount

    ! 6) ATMS quality flag check (qual. flag elements 33078,33079,33080,33081)

    !  33078 Geolocation quality code     qcflag1(ii,1)  code value = 0-15 (0= OK, 15=misg)
    !  33079 Granule level quality flags  qcflag1(ii,2)  16 bit flag  (start bit 6(2^5)=32) (misg=2^16-1 = 65535)
    !  33080 Scan level quality flags     qcflag1(ii,3)  20 bit flag  (start bit 7(2^6)=64) (misg=2^20-1) 
    !  33081 Channel data quality flags   qcflag2        12 bit flag  (start bit 3(2^2)=4)  (misg=2^12-1)
    !
    !  See http://www.wmo.int/pages/prog/www/WMOCodes/WMO306_vI2/2010edition/BUFRver16/BUFR_16_0_0_TableD.pdf

    indx1 = 1
    do ii = 1,nt 
      fail1 = .false.
      fail = .false.
      if ( (qcflag1(ii,1) > 0) .or. (qcflag1(ii,2) >= 32) .or. (qcflag1(ii,3) >= 64) ) then
        write(*,*) 'WARNING: INFO BLOCK QC flag(s) indicate problem with data'
        write(*,*) ' ele33078 = ',qcflag1(ii,1),' ele33079 = ',qcflag1(ii,2),' ele33080 = ', qcflag1(ii,3)
        write(*,*) ' lat, lon = ', zlat(ii), zlon(ii)
        fail1 = .true.
        if ( grossrej(ii) ) write(*,*) ' NOTE: grossrej is also true for this point!'
      end if
      do jj = 1,mwbg_atmsMaxNumChan
        fail2 = .false.
        if ( qcflag2(indx1+jj-1) >= 4 ) then
          !write(*,*) 'WARNING: DATA BLOCK QC flag ele33081 = ', qcflag2(indx1+jj-1)
          !write(*,*) '    Lat, lon, channel = ', zlat(ii), zlon(ii), ican(indx1+jj-1)
          fail2 = .true.
          fail = .true.
          !if ( (.not. fail1) .and. grossrej(ii) ) write(*,*) ' NOTE: grossrej is also true for this point!'
        end if
        if ( fail2 .or. fail1 ) lqc(ii,jj) = .true.
      end do
      if ( fail ) write(*,*) 'WARNING: DATA BLOCK QC flag ele33081 >= 4 for one or more channels! lat, lon = ', zlat(ii), zlon(ii)
      indx1 = indx1 + mwbg_atmsMaxNumChan
    end do
     
    write(*,*) 'mwbg_firstQcCheckAtms: Total number of data processed in this box = ', nt*mwbg_atmsMaxNumChan
    write(*,*) '         Total number of data flagged in this box   = ', COUNT(lqc)
    write(*,*) ' '

    return
  end subroutine mwbg_firstQcCheckAtms

  !--------------------------------------------------------------------------
  !  mwbg_nrlFilterAtms
  !--------------------------------------------------------------------------
  subroutine mwbg_nrlFilterAtms(ni, ztbcor, biasCorr, pangl, plat, ilansea, iglace, waterobs, &
                                grossrej, clw, si_ecmwf, si_bg, iNumSeaIce, iRej,SeaIce)
    !OBJET          Compute the following parameters using 5 ATMS channels:
    !                  - sea ice, 
    !                  - cloud liquid water (clw), 
    !                  - 2 scattering indices (si) (ECMWF, Bennartz-Grody)
    !               The five channels used are: 23Ghz, 31Ghz, 50Ghz, 89Ghz, and 165Ghz.
    !
    !NOTES*
    !                - open water points are converted to sea-ice points if sea ice concentration >= 0.55
    !                   and iglace (itt or terrain type) is changed accordingly
    !                - clw are missing when out-of-range parameters/Tb detected or grossrej = .true.
    !                - clw and si only computed over open water away from coasts and sea-ice
    !                - clw and si = -99.0 where value cannot be computed.
    !
    !REFERENCES     Ben Ruston, NRL Monterey
    !                  JCSDA Seminar 12/12/12: Impact of NPP Satellite Assimilation in the U.S. Navy Global Modeling System
    !
    !
    !ARGUMENTS      - ier         - output - error return code for each location:
    !                                        0, ok,  
    !                                        1, input parameter out of range or grossrej=.true. 
    !               - ni          - input  -  number of points to process (= NT)
    !               - tb23        - input  -  23Ghz brightness temperature (K) -- ch. 1
    !               - tb31        - input  -  31Ghz brightness temperature (K) -- ch. 2
    !               - tb50        - input  -  50Ghz brightness temperature (K) -- ch. 3
    !               - tb89        - input  -  89Ghz brightness temperature (K) -- ch. 16
    !               - tb165       - input  -  165Ghz brightness temperature (K) -- ch. 17
    !               - pangl       - input  -  satellite zenith angle (deg.)
    !               - plat        - input  -  latitude (deg.)
    !               - ilansea     - input  -  land/sea indicator (0=land, 1=ocean)
    !               - iglace      - in/out -  terrain type (0=ice, -1 otherwise)
    !               - waterobs    - in/out -  .true. if open water point (away from coasts and sea-ice)
    !               - grossrej    - input  -  .true. if any channel had a gross error from mwbg_grossValueCheck
    !               - clw         - output -  cloud liquid water (kg/m**2) from tb23 & tb31
    !               - si_ecmwf    - output -  ECMWF scattering index from tb89 & tb165
    !               - si_bg       - output -  Bennartz-Grody scattering index from tb89 & tb165
    !               - iNumSeaIce  - in/out -  running counter for number of open water points
    !                                       with sea-ice detected (from algorithm)
    !               - iRej        - in/out -  running counter for number of locations with bad
    !                                       pangl, plat, ilansea, or with grossrej=true
    !               - SeaIce      - output -  computed sea-ice fraction from tb23 & tb50 
    !
    !               - ice         - internal -  sea ice
    !             
    !
    ! Notes: In the case where an output parameter cannot be calculated, the
    !        value of this parameter is set to the missing value, i.e. -99.
    !
    implicit none

    integer    ::  i

    integer, intent(in)                   ::  ni
    integer, intent(out)                  ::  iNumSeaIce
    integer, intent(in)                   ::  ilansea(:)
    integer, intent(inout)                ::  iglace(:)
    integer, intent(out)                  ::  iRej
   
    
    logical, intent(in)                   ::  grossrej(:)
    logical, intent(inout)                ::  waterobs(:)

    real, intent(in)                      ::  ztbcor(:)
    real, intent(in)                      ::  biasCorr(:)
    real, intent(in)                      ::  pangl(:)
    real, intent(in)                      ::  plat(:)
    real, allocatable, intent(out)        ::  clw (:)
    real, allocatable, intent(out)        ::  si_ecmwf(:) 
    real, allocatable, intent(out)        ::  si_bg(:) 
    real, allocatable, intent(out)        ::  SeaIce(:) 

    ! Locals
    integer                               :: ier(ni)
    real                                  ::  ice(ni)
    real                                  :: tb23(ni)
    real                                  :: tb31(ni)
    real                                  :: tb50(ni)
    real                                  :: tb89(ni)
    real                                  :: tb165(ni)
    real                                  :: bcor23(ni)
    real                                  :: bcor31(ni)
    real                                  :: bcor50(ni)
    real                                  :: bcor89(ni)
    real                                  :: bcor165(ni)
    integer                               :: indx1
    integer                               :: indx2
    integer                               :: ii
    real                                  :: aa
    real                                  :: deltb
    real                                  :: abslat
    real                                  :: cosz
    real                                  :: t23
    real                                  :: t31
    real                                  :: t50
    real                                  :: t89
    real                                  :: t165
    real, parameter                       ::  rmisg = -99.


    ! Allocation
    call mwbg_allocate1DRealArray(clw,ni)
    call mwbg_allocate1DRealArray(si_ecmwf,ni)
    call mwbg_allocate1DRealArray(si_bg,ni)
    call mwbg_allocate1DRealArray(SeaIce,ni)
    ! extract required channels:
    !  23 Ghz = AMSU-A 1 = ATMS channel 1 
    !  31 Ghz = AMSU-A 2 = ATMS channel 2
    !  50 Ghz = AMSU-A 3 = ATMS channel 3
    !  53 Ghz = AMSU-A 5 = ATMS channel 6
    !  89 Ghz = AMSU-A15 = ATMS channel 16
    ! 150 Ghz = AMSU-B 2 = ATMS channel 17
    !

    indx1 = 1
    do ii = 1, ni
      indx2 = ii*mwbg_atmsMaxNumChan
      tb23(ii)      = ztbcor(indx1)
      bcor23(ii)    = biasCorr(indx1)
      tb31(ii)      = ztbcor(indx1+1)
      bcor31(ii)    = biasCorr(indx1+1)
      tb50(ii)      = ztbcor(indx1+2)
      bcor50(ii)    = biasCorr(indx1+2)
      tb89(ii)      = ztbcor(indx1+15)
      bcor89(ii)    = biasCorr(indx1+15)
      tb165(ii)    = ztbcor(indx1+16)
      bcor165(ii)    = biasCorr(indx1+16)
      indx1 = indx2 + 1
    end do
    
    ier = 0
    ! 0) Allocation
    ! a prevoir pour clw, si_ecmwf, si_bg, SeaIce


    ! 1) Initialise parameters:
    do i = 1, ni
      ice(i)      = rmisg
      clw(i)      = rmisg
      si_ecmwf(i) = rmisg
      si_bg(i)    = rmisg
      SeaIce(i)   = 0.0
    end do

    ! 2) Validate input parameters:
    do i = 1, ni
      if ( pangl(i)   .lt.   0.  .or. &
           pangl(i)   .gt.  70.  .or. &
           plat(i)    .lt. -90.  .or. & 
           plat(i)    .gt.  90.  .or. &  
           ilansea(i) .lt.   0   .or. & 
           ilansea(i) .gt.   1        ) then
         ier(i) = 1
      end if

      ! Skip computations for points where all data are rejected  (bad Tb ANY channel)       
      if ( grossrej(i) ) ier(i) = 1 

    end do

    ! 3) Compute parameters:
    do i = 1, ni

      if ( ier(i) .eq. 0 ) then

        abslat = abs(plat(i))
        cosz   = cosd(pangl(i))

        if ( mwbg_useUnbiasedObsForClw ) then
          t23 = tb23(i)
          t31 = tb31(i)
          t50 = tb50(i)
          t89 = tb89(i)
          t165 = tb165(i)
        else
          t23 = tb23(i) - bcor23(i)
          t31 = tb31(i) - bcor31(i)
          t50 = tb50(i) - bcor50(i)
          t89 = tb89(i) - bcor89(i)
          t165 = tb165(i) - bcor165(i)
        end if
        deltb = t89 - t165

        ! Check for sea-ice over water points. Set terrain type to 0 if ice>=0.55 detected.
        if ( ilansea(i) .eq. 1 ) then  ! water point

          if ( abslat .lt. 50. ) then
            ice(i) = 0.0
          else
            ice(i) = 2.85 + 0.020*t23 - 0.028*t50
          end if
          
          SeaIce(i) = ice(i)
          
          if ( ice(i) .ge. 0.55 .and. waterobs(i) ) then
            iNumSeaIce = iNumSeaIce + 1
            waterobs(i) = .false.
            iglace(i) = 0
          end if
          
        end if

        ! Compute CLW and Scattering Indices (over open water only)
        if ( waterobs(i) ) then
          if ( t23 .lt. 284. .and. t31 .lt. 284. ) then
            aa = 8.24 - (2.622 - 1.846*cosz)*cosz
            clw(i) = aa + 0.754*alog(285.0-t23) - 2.265*alog(285.0-t31)
            clw(i) = clw(i)*cosz
            if ( clw(i) .lt. 0.0 ) clw(i) = 0.0
          end if
          si_ecmwf(i) = deltb - (-46.94 + 0.248*pangl(i))
          si_bg(i)    = deltb - (-39.201 + 0.1104*pangl(i))
        end if

      else  ! ier(i) .eq. 1 case
         iRej = iRej + 1

      end if ! if ( ier(i) .eq. 0 )

      if ( mwbg_debug .and. (i .le. 100) ) then
        write(*,*) ' '
        write(*,*) ' i,tb23(i),tb31(i),tb50(i),tb89(i),tb165(i),pangl(i),plat(i), ilansea(i) = ', &
     &             i,tb23(i),tb31(i),tb50(i),tb89(i),tb165(i),pangl(i),plat(i), ilansea(i)
        write(*,*) ' ier(i),ice(i),clw(i),si_ecmwf(i),si_bg(i),iglace(i),waterobs(i) =',ier(i),ice(i),&
     &             clw(i),si_ecmwf(i),si_bg(i),iglace(i),waterobs(i)
      end if

    end do   ! i loop over ni points

  end subroutine mwbg_nrlFilterAtms

  !--------------------------------------------------------------------------
  ! mwbg_flagDataUsingNrlCriteria 
  !--------------------------------------------------------------------------
  subroutine mwbg_flagDataUsingNrlCriteria(nt,ztbcor, biasCorr, rclw, scatec, scatbg, SeaIce, grossrej, waterobs, &
                                          useUnbiasedObsForClw, mwbg_realMisg, iwvreject, & 
                                          cloudobs, precipobs,  cldcnt, ident, riwv, zdi)

    !:Purpose:                       Set the  Information flag (ident) values (new BURP element 025174 in header)
    !                                BIT    Meaning
    !                                - 0     off=land or sea-ice, on=open water away from coast
    !                                - 1     Mean 183 Ghz [ch. 18-22] is missing
    !                                - 2     CLW is missing (over water)
    !                                - 3     CLW > clw_atms_nrl_LTrej (0.175 kg/m2) (cloudobs)
    !                                - 4     scatec/scatbg > Lower Troposphere limit 9/10 (precipobs)
    !                                - 5     Mean 183 Ghz [ch. 18-22] Tb < 240K
    !                                - 6     CLW > clw_atms_nrl_UTrej (0.200 kg/m2)
    !                                - 7     Dryness Index rejection (for ch. 22)
    !                                - 8     scatec/scatbg > Upper Troposphere limit 18/15
    !                                - 9     Dryness Index rejection (for ch. 21)
    !                               - 10     Sea ice > 0.55 detected
    !                               - 11     Gross error in Tb (any chan.)  (all channels rejected)

    ! Arguments
    integer, intent(in)                        :: nt  
    real, intent(in)                           :: ztbcor(:)
    real, intent(in)                           :: biasCorr(:)
    real, intent(in)                           :: rclw (:)
    real, intent(in)                           :: scatec(:) 
    real, intent(in)                           :: scatbg(:) 
    real, intent(in)                           :: SeaIce (:)

    logical, intent(in)                        :: useUnbiasedObsForClw
    real, intent(in)                           :: mwbg_realMisg 
    logical, intent(in)                        :: grossrej(:)
    logical, intent(in)                        :: waterobs(:)
    integer, intent(inout)                     :: cldcnt 
    logical, allocatable, intent(out)          :: cloudobs(:)
    logical, allocatable, intent(out)          :: iwvreject(:)
    logical, allocatable, intent(out)          :: precipobs(:)
    integer, allocatable, intent(out)          :: ident(:)
    real, allocatable, intent(out)             :: zdi(:)
    real, allocatable, intent(out)             :: riwv(:)

    ! Locals
    integer                                    :: indx1
    integer                                    :: indx2
    integer                                    :: ii
    integer                                    :: n_cld
    real, allocatable                          :: ztb_amsub3(:)  
    real, allocatable                          :: bcor_amsub3(:)  
    real, allocatable                          :: ztb_amsub5(:) 
    real, allocatable                          :: bcor_amsub5(:) 
    real                                       ::  ztb183(5)


    call mwbg_allocate1DLogicalArray(cloudobs, nt)
    call mwbg_allocate1DLogicalArray(iwvreject, nt)
    call mwbg_allocate1DIntegerArray(ident, nt)
    call mwbg_allocate1DLogicalArray(precipobs, nt)
    call mwbg_allocate1DRealArray(riwv, nt)
    call mwbg_allocate1DRealArray(ztb_amsub3, nt)
    call mwbg_allocate1DRealArray(bcor_amsub3, nt)
    call mwbg_allocate1DRealArray(ztb_amsub5, nt)
    call mwbg_allocate1DRealArray(bcor_amsub5, nt)

    ! To begin, assume that all obs are good.
    ident(:) = 0
    cloudobs(:)  = .false.
    iwvreject(:) = .false.
    precipobs(:) = .false.

    ! Extract Tb for channels 16 (AMSU-B 1) and 17 (AMSU-B 2) for Bennartz SI
    ! Extract Tb for channels 22 (AMSU-B 3) and 18 (AMSU-B 5) for Dryness Index (DI)

    indx1 = 1
    do ii = 1, nt
      indx2 = ii*mwbg_atmsMaxNumChan
      ztb_amsub3(ii) = ztbcor(indx1+21)
      bcor_amsub3(ii) = biasCorr(indx1+21)
      ztb_amsub5(ii) = ztbcor(indx1+17)
      bcor_amsub5(ii) = biasCorr(indx1+17)
      indx1 = indx2 + 1
    end do


    ! Flag data using NRL criteria

    ! Compute Mean 183 Ghz [ch. 18-22] Tb (riwv)
    riwv = -99.0
    indx1 = 1
    do ii = 1, nt
      indx2 = ii*mwbg_atmsMaxNumChan
      if (.not.grossrej(ii)) then
        if ( useUnbiasedObsForClw ) then
          ztb183(1) = ztbcor(indx1+17)
          ztb183(2) = ztbcor(indx1+18)
          ztb183(3) = ztbcor(indx1+19)
          ztb183(4) = ztbcor(indx1+20)
          ztb183(5) = ztbcor(indx1+21)
        else
          ztb183(1) = ztbcor(indx1+17) - biasCorr(indx1+17)
          ztb183(2) = ztbcor(indx1+18) - biasCorr(indx1+18)
          ztb183(3) = ztbcor(indx1+19) - biasCorr(indx1+19)
          ztb183(4) = ztbcor(indx1+20) - biasCorr(indx1+20)
          ztb183(5) = ztbcor(indx1+21) - biasCorr(indx1+21)
        end if
        riwv(ii)  = sum(ztb183)/5.0
        if ( riwv(ii) < mean_Tb_183Ghz_min ) iwvreject(ii) = .true.
      else
        iwvreject(ii) = .true.
      end if
      indx1 = indx2 + 1
    end do

    !  Set bits in ident flag to identify where various data selection criteria are met
    !     precipobs = .true  where ECMWF or BG scattering index > min_threshold (LT)
    !     cloudobs  = .true. where CLW > min_threshold (LT) or if precipobs = .true

    where ( grossrej ) ident = IBSET(ident,11)
    where ( scatec .gt. scatec_atms_nrl_LTrej .or. scatbg .gt. scatbg_atms_nrl_LTrej ) precipobs = .true.
    n_cld = count(rclw .gt. clw_atms_nrl_LTrej)
    cldcnt  = cldcnt  + n_cld
    where ( (rclw .gt. clw_atms_nrl_LTrej) .or. precipobs ) cloudobs = .true.
    where ( waterobs )  ident = IBSET(ident,0)
    where ( iwvreject ) ident = IBSET(ident,5)
    where ( precipobs ) ident = IBSET(ident,4)
    where ( rclw .gt. clw_atms_nrl_LTrej) ident = IBSET(ident,3)
    where ( rclw .gt. clw_atms_nrl_UTrej) ident = IBSET(ident,6)
    where ( scatec .gt. scatec_atms_nrl_UTrej .or. scatbg .gt. scatbg_atms_nrl_UTrej ) ident = IBSET(ident,8)
    where ( SeaIce .ge. 0.55 ) ident = IBSET(ident,10)
      
    where ( waterobs .and. (rclw == -99.) ) ident = IBSET(ident,2)
    where ( riwv == -99.)                   ident = IBSET(ident,1)

    ! Compute the simple AMSU-B Dryness Index zdi for all points = Tb(ch.3)-Tb(ch.5)
    if ( useUnbiasedObsForClw ) then
      where ( .not.grossrej )
        zdi = ztb_amsub3 - ztb_amsub5
      elsewhere
        zdi = mwbg_realMisg
      end where
    else
      where ( .not.grossrej )
        zdi = (ztb_amsub3 - bcor_amsub3) - (ztb_amsub5 - bcor_amsub5)
      elsewhere
        zdi = mwbg_realMisg
      end where
    end if

  end subroutine mwbg_flagDataUsingNrlCriteria

  !--------------------------------------------------------------------------
  !  mwbg_reviewAllcriteriaforFinalFlags
  !--------------------------------------------------------------------------
  subroutine mwbg_reviewAllcriteriaforFinalFlags(nt,nval, lqc, grossrej, waterobs, precipobs, rclw, scatec, &
                                                 scatbg, iwvreject, riwv, IMARQ, globMarq, zdi, ident, mwbg_realMisg, &
                                                 drycnt, landcnt, rejcnt, iwvcnt, pcpcnt, flgcnt)

    !:Purpose:                   Review all the checks previously made to determine which obs are to be accepted
    !                            for assimilation and which are to be flagged for exclusion (lflagchn). 
    !                            - grossrej()  = .true. if any channel had a gross error at the point
    !                            - cloudobs()  = .true. if CLW > clw_atms_nrl_LTrej (0.175) or precipobs
    !                            - precipobs() = .true. if precip. detected through NRL scattering indices
    !                            - waterobs()  = .true. if open water point
    !                            - iwvreject() = .true. if Mean 183 Ghz [ch. 18-22] Tb < 240K (too dry for ch.20-22 over land)
    ! Arguments
    integer, intent(in)                        :: nt  
    integer, intent(in)                        :: nval 
    logical, intent(in)                        :: lqc(:,:)
    real, intent(inout)                        :: rclw (:)
    real, intent(in)                           :: scatec(:) 
    real, intent(in)                           :: scatbg(:) 
    logical, intent(in)                        :: grossrej(:)
    logical, intent(in)                        :: waterobs(:)
    logical, intent(in)                        :: iwvreject(:)
    logical, intent(in)                        :: precipobs(:)
    integer, intent(inout)                     :: ident(:)
    real, intent(in)                           :: zdi(:)
    real, intent(inout)                        :: riwv(:)
    integer, intent(inout)                     :: IMARQ(:)
    integer, intent(inout)                     :: globMarq(:)
    real, intent(in)                           :: mwbg_realMisg 
    integer, intent(inout)                     :: drycnt 
    integer, intent(inout)                     :: landcnt 
    integer, intent(inout)                     :: rejcnt 
    integer, intent(inout)                     :: iwvcnt 
    integer, intent(inout)                     :: pcpcnt 
    integer, intent(inout)                     :: flgcnt 

    ! Locals
    logical, allocatable                       :: lflagchn(:,:)
    integer                                    :: indx1
    integer                                    :: indx2
    integer                                    :: ii, kk, j, ipos, iidata
    integer                                    :: n_cld
    real, allocatable                          :: ztb_amsub3(:)  
    real, allocatable                          :: bcor_amsub3(:)  
    real, allocatable                          :: ztb_amsub5(:) 
    real, allocatable                          :: bcor_amsub5(:) 
    real                                       :: ztb183(5)


    ! Allocation
    call mwbg_allocate2DLogicalArray(lflagchn,nt, nval)

    lflagchn(:,:) = lqc(:,:)  ! initialize with flags set in mwbg_firstQcCheckAtms
    do kk = 1, nt
      ! Reject all channels if gross Tb error detected in any channel or other problems 
      if ( grossrej(kk) ) then
        lflagchn(kk,:) = .true.
      else
        ! OVER LAND OR SEA-ICE,
        !    -- CLW/SI not determined over land
        !    -- surface emissivity effects lower tropospheric and window channels     
        !    -- reject window & lower tropospheric channels 1-6, 16-19
        !    -- reject ch. 20-22 if iwvreject = .true.  [ Mean 183 Ghz [ch. 18-22] Tb < 240K ]
        !    -- check DI for AMSU-B like channels
        if  ( .not. waterobs(kk) ) then
          lflagchn(kk,1:mwbg_atmsNumSfcSensitiveChannel)     = .true.      ! AMSU-A 1-6
          lflagchn(kk,16:19)     = .true.      ! AMSU-B (like 1,2,5)
          if ( iwvreject(kk) ) lflagchn(kk,20:22) = .true.  ! AMSU-B (like 4,3)

          ! Dryness index (for AMSU-B channels 19-22 assimilated over land/sea-ice)
          ! Channel AMSUB-3 (ATMS channel 22) is rejected for a dryness index >    0.
          !                (ATMS channel 21) is rejected for a dryness index >   -5.
          ! Channel AMSUB-4 (ATMS channel 20) is rejected for a dryness index >   -8.
          if ( zdi(kk) > 0.0 ) then
            lflagchn(kk,22) = .true.
            ident(kk) = IBSET(ident(kk),7)
          end if
          if ( zdi(kk) > -5.0 ) then
            lflagchn(kk,21) = .true.
            ident(kk) = IBSET(ident(kk),9)
            drycnt = drycnt + 1
          end if
          if ( zdi(kk) > -8.0 ) then
            lflagchn(kk,20) = .true.
          end if
        end if
        ! OVER WATER,
        !    -- reject ch. 1-6, 16-20 if CLW > clw_atms_nrl_LTrej or CLW = -99.0
        !    -- reject ch. 7-9, 21-22 if CLW > clw_atms_nrl_UTrej or CLW = -99.0
        !    -- reject ch. 1-6, 16-22 if scatec > 9  or scatec = -99.0
        !    -- reject ch. 7-9        if scatec > 18 or scatec = -99.0
        !    -- reject ch. 1-6        if scatbg > 10 or scatbg = -99.0
        !    -- reject ch. 7-9        if scatbg > 15 or scatbg = -99.0
        !    -- reject ch. 16-22      if iwvreject = .true.   [ Mean 183 Ghz [ch. 18-22] Tb < 240K ]
        if  ( waterobs(kk) ) then
          if ( rclw(kk)   >  clw_atms_nrl_LTrej )  then
            lflagchn(kk,1:mwbg_atmsNumSfcSensitiveChannel) = .true.
            lflagchn(kk,16:20) = .true. 
          end if
          if ( rclw(kk)   >  clw_atms_nrl_UTrej )  then
            lflagchn(kk,7:9)   = .true.
            lflagchn(kk,21:22) = .true. 
          end if
          if ( scatec(kk) >  scatec_atms_nrl_LTrej ) then
            lflagchn(kk,1:mwbg_atmsNumSfcSensitiveChannel) = .true.
            lflagchn(kk,16:22) = .true.
          end if
          if ( scatec(kk) >  scatec_atms_nrl_UTrej ) lflagchn(kk,7:9) = .true.
          if ( scatbg(kk) >  scatbg_atms_nrl_LTrej ) lflagchn(kk,1:mwbg_atmsNumSfcSensitiveChannel) = .true.
          if ( scatbg(kk) >  scatbg_atms_nrl_UTrej ) lflagchn(kk,7:9) = .true.
          if ( iwvreject(kk) ) lflagchn(kk,16:22) = .true.
          if ( rclw(kk) == -99. ) then
            ident(kk) = IBSET(ident(kk),2)
            lflagchn(kk,1:9)   = .true.
            lflagchn(kk,16:22) = .true.
          end if
          if ( riwv(kk) == -99. ) then     ! riwv = mean_Tb_183Ghz
            ident(kk) = IBSET(ident(kk),1)
            lflagchn(kk,16:22) = .true.
          end if           
        end if
    
      end if

      if ( .not. waterobs(kk) ) landcnt  = landcnt  + 1
      if ( grossrej(kk) )  rejcnt = rejcnt + 1
      if ( iwvreject(kk))  iwvcnt = iwvcnt + 1
      if ( precipobs(kk) .and. waterobs(kk) ) then
        pcpcnt = pcpcnt + 1
      end if
        
      if ( ANY(lflagchn(kk,:)) ) flgcnt = flgcnt + 1
    end do

    ! RESET riwv array to ECMWF scattering index for output to BURP file
    riwv(:) = scatec(:)
     ! Set missing rclw and riwv to BURP missing value (mwbg_realMisg)
    where (rclw == -99. ) rclw = mwbg_realMisg
    where (riwv == -99. ) riwv = mwbg_realMisg

    ! Modify data flag values (set bit 7) for rejected data  
    ipos=0
    do kk =1, nt
      do j = 1, nval
        ipos = ipos + 1
        if (lflagchn(kk,j)) then
          IMARQ(ipos) = IBSET(IMARQ(ipos),7)
        end if
      end do
    end do

    ! Set bit 6 in 24-bit global flags if any data rejected  
    do kk =1, nt
      if ( ANY(lflagchn(kk,:)) ) globMarq(kk) = IBSET(globMarq(kk),6)
    end do
  
  end subroutine mwbg_reviewAllcriteriaforFinalFlags


  function calcStateDepObsErr_r4(clwThresh1,clwThresh2,sigmaThresh1,sigmaThresh2,clw_avg) result(sigmaObsErrUsed)
    implicit none
    real :: clwThresh1
    real :: clwThresh2
    real :: sigmaThresh1
    real :: sigmaThresh2
    real :: clw_avg
    real :: sigmaObsErrUsed

    if ( clw_avg <= clwThresh1 ) then
      sigmaObsErrUsed = sigmaThresh1
    else if ( clw_avg >  clwThresh1 .and. & 
                  clw_avg <= clwThresh2 ) then
      sigmaObsErrUsed = sigmaThresh1 + &
                      (sigmaThresh2 - sigmaThresh1) / &
                      (clwThresh2 - clwThresh1) * &
                      (clw_avg - clwThresh1) 
    else
      sigmaObsErrUsed = sigmaThresh2
    end if

  end function calcStateDepObsErr_r4

end module bgckmicrowave_mod
