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
  use midasMpi_mod
  use MathPhysConstants_mod
  use utilities_mod
  use obsSpaceData_mod
  use tovs_nl_mod
  use obsErrors_mod
  use codtyp_mod

  implicit none
  save
  private

  ! Public functions/subroutines
  public :: mwbg_bgCheckMW
  public :: mwbg_computeMwhs2SurfaceType

  real    :: mwbg_clwQcThreshold
  real    :: mwbg_cloudyClwThresholdBcorr
  real    :: mwbg_siQcOverWaterThreshold  ! for AMSUB/MHS
  real    :: mwbg_cloudySiThresholdBcorr  ! for AMSUB/MHS
  logical :: mwbg_debug
  logical :: mwbg_useUnbiasedObsForClw 

  integer, parameter :: mwbg_maxScanAngle=98
  real,    parameter :: mwbg_realMissing=-99. 
  integer, parameter :: mwbg_intMissing=-1

  ! Module variable

  integer, parameter :: mwbg_atmsNumSfcSensitiveChannel = 6
  character(len=128), parameter :: fileMgLg='fstglmg'  ! glace de mer file
  ! Upper limit for CLW (kg/m**2) for Tb rejection over water
  real,   parameter :: clw_atms_nrl_LTrej=0.175      ! lower trop chans 1-6, 16-20
  real,   parameter :: clw_atms_nrl_UTrej=0.2        ! upper trop chans 7-9, 21-22
  real,   parameter :: clw_mwhs2_nrl_LTrej=0.175
  real,   parameter :: clw_mwhs2_nrl_UTrej=0.2
  ! Other NRL thresholds
  real,   parameter :: scatec_atms_nrl_LTrej=9.0     ! lower trop chans 1-6, 16-22
  real,   parameter :: scatec_atms_nrl_UTrej=18.0    ! upper trop chans 7-9
  real,   parameter :: scatbg_atms_nrl_LTrej=10.0    ! lower trop chans 1-6
  real,   parameter :: scatbg_atms_nrl_UTrej=15.0    ! upper trop chans 7-9
  real,   parameter :: scatec_mwhs2_nrl_LTrej=9.0    ! all MWHS-2 channels (over water)
  real,   parameter :: scatbg_mwhs2_cmc_LANDrej=0.0  ! all MWHS-2 channels (all surfaces)
  real,   parameter :: scatbg_mwhs2_cmc_ICErej=40.0
  real,   parameter :: scatbg_mwhs2_cmc_SEA=15.0
  real,   parameter :: mean_Tb_183Ghz_min=240.0      ! min. value for Mean(Tb) chans. 18-22 

  integer, parameter :: mwbg_maxNumSat  = 13
  integer, parameter :: mwbg_maxNumChan = 100
  integer, parameter :: mwbg_maxNumTest = 16

  ! OPTION for values of MG (land/sea mask) and LG (ice) at each observation
  ! point using values on 5x5 mesh centered at each point.
  ! ilsmOpt = 1 --> use MAX value from all 25 mesh points
  ! ilsmOpt = 2 --> use value at central mesh point (obs location)
  ! ilsmOpt = 3 --> use AVG value from all 25 mesh points
  integer, parameter :: ilsmOpt = 2

  integer            :: rejectionCodArray (mwbg_maxNumTest, &
                                           mwbg_maxNumChan, &
                                           mwbg_maxNumSat) ! number of rejection 
  !                                                          per sat. per channl per test
  integer            :: rejectionCodArray2 (mwbg_maxNumTest, &
                                           mwbg_maxNumChan, &
                                           mwbg_maxNumSat) ! number of rejection per channl per test
  !                                                          for ATMS 2nd category of tests

  ! namelist variables
  character(len=9)              :: instName                      ! instrument name
  real                          :: clwQcThreshold                ! 
  real                          :: cloudyClwThresholdBcorr       !
  real                          :: siQcOverWaterThreshold        ! scattering index over water for AMSUB/MHS
  real                          :: cloudySiThresholdBcorr        !
  logical                       :: useUnbiasedObsForClw          !
  logical                       :: RESETQC                       ! reset Qc flags option
  logical                       :: modLSQ                        !
  logical                       :: debug                         ! debug mode


  namelist /nambgck/instName, clwQcThreshold, &
                    useUnbiasedObsForClw, debug, RESETQC,  &
                    cloudyClwThresholdBcorr, modLSQ, &
                    siQcOverWaterThreshold, cloudySiThresholdBcorr
                    

contains

  subroutine mwbg_init()
    !
    !:Purpose: This subroutine reads the namelist section NAMBGCK
    !          for the module.
    implicit none

    ! Locals:
    integer :: nulnam, ierr
    integer, external :: fnom, fclos

    ! Default values for namelist variables
    debug                   = .false.
    clwQcThreshold          = 0.3 
    useUnbiasedObsForClw    = .false.
    cloudyClwThresholdBcorr = 0.05
    siQcOverWaterThreshold  = 15.0
    cloudySiThresholdBcorr  = 5
    RESETQC                 = .false.
    modLSQ                  = .false.

    nulnam = 0
    ierr = fnom(nulnam, './flnml','FTN+SEQ+R/O', 0)
    read(nulnam, nml=nambgck, iostat=ierr)
    if (ierr /= 0) call utl_abort('mwbg_init: Error reading namelist')
    if (mmpi_myid == 0) write(*, nml=nambgck)
    ierr = fclos(nulnam)

    mwbg_debug = debug
    mwbg_clwQcThreshold = clwQcThreshold
    mwbg_useUnbiasedObsForClw = useUnbiasedObsForClw
    mwbg_cloudyClwThresholdBcorr = cloudyClwThresholdBcorr
    mwbg_siQcOverWaterThreshold = siQcOverWaterThreshold
    mwbg_cloudySiThresholdBcorr = cloudySiThresholdBcorr

  end subroutine mwbg_init 

  !--------------------------------------------------------------------------
  ! ISRCHEQI function
  !--------------------------------------------------------------------------
  function ISRCHEQI(KLIST, KLEN, KENTRY) result(ISRCHEQI_out)
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
  ! extractParamForGrodyRun
  !--------------------------------------------------------------------------  
  subroutine extractParamForGrodyRun(KCANO, ptbo, ptbomp, ptbcor, KNT, KNO, &
                                     tb23,   tb31,   tb50,   tb53,   tb89, &
                                     tb23FG, tb31FG, tb50FG, tb53FG, tb89FG)

    !:Purpose: Compute  Grody parameters by extracting tb for required channels:
    !          - 23 Ghz = AMSU-A 1 = channel #28
    !          - 31 Ghz = AMSU-A 2 = channel #29
    !          - 50 Ghz = AMSU-A 3 = channel #30
    !          - 53 Ghz = AMSU-A 5 = channel #32
    !          - 89 Ghz = AMSU-A15 = channel #42
    implicit none
    ! Arguments
    integer,     intent(in)               :: KCANO(KNO,KNT)         ! observations channels
    real,        intent(in)               :: ptbo(KNO,KNT)          ! radiances
    real,        intent(in)               :: ptbomp(KNO,KNT)        ! radiances o-p
    real,        intent(in)               :: ptbcor(KNO,KNT)        ! correction aux radiances
    integer,     intent(in)               :: KNO                    ! nombre de canaux des observations 
    integer,     intent(in)               :: KNT                    ! nombre de tovs
    real,        intent(out)              :: tb23(KNT)              ! radiance frequence 23 Ghz   
    real,        intent(out)              :: tb31(KNT)              ! radiance frequence 31 Ghz
    real,        intent(out)              :: tb50(KNT)              ! radiance frequence 50 Ghz  
    real,        intent(out)              :: tb53(KNT)              ! radiance frequence 53 Ghz  
    real,        intent(out)              :: tb89(KNT)              ! radiance frequence 89 Ghz  
    real,        intent(out)              :: tb23FG(KNT)            ! radiance frequence 23 Ghz   
    real,        intent(out)              :: tb31FG(KNT)            ! radiance frequence 31 Ghz
    real,        intent(out)              :: tb50FG(KNT)            ! radiance frequence 50 Ghz  
    real,        intent(out)              :: tb53FG(KNT)            ! radiance frequence 53 Ghz  
    real,        intent(out)              :: tb89FG(KNT)            ! radiance frequence 89 Ghz        

    ! Locals
    integer                               :: channelval
    integer                               :: nDataIndex
    integer                               :: nChannelIndex

    do nDataIndex=1,KNT
      do nChannelIndex=1,KNO
        channelval = KCANO(nChannelIndex,nDataIndex)
        if ( ptbo(nChannelIndex,nDataIndex) .ne. mwbg_realMissing ) then
          if ( ptbcor(nChannelIndex,nDataIndex) .ne. mwbg_realMissing ) then
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

          if ( channelval .eq. 28 ) tb23FG(nDataIndex) = ptbo(nChannelIndex,nDataIndex) - &
                                                         ptbomp(nChannelIndex,nDataIndex)
          if ( channelval .eq. 29 ) tb31FG(nDataIndex) = ptbo(nChannelIndex,nDataIndex) - &
                                                         ptbomp(nChannelIndex,nDataIndex)
          if ( channelval .eq. 30 ) tb50FG(nDataIndex) = ptbo(nChannelIndex,nDataIndex) - &
                                                         ptbomp(nChannelIndex,nDataIndex)
          if ( channelval .eq. 32 ) tb53FG(nDataIndex) = ptbo(nChannelIndex,nDataIndex) - &
                                                         ptbomp(nChannelIndex,nDataIndex)
          if ( channelval .eq. 42 ) tb89FG(nDataIndex) = ptbo(nChannelIndex,nDataIndex) - &
                                                         ptbomp(nChannelIndex,nDataIndex)
        else
          if ( channelval .eq. 28 ) tb23(nDataIndex) = 0.
          if ( channelval .eq. 29 ) tb31(nDataIndex) = 0.
          if ( channelval .eq. 30 ) tb50(nDataIndex) = 0.
          if ( channelval .eq. 32 ) tb53(nDataIndex) = 0.
          if ( channelval .eq. 42 ) tb89(nDataIndex) = 0.

          if ( channelval .eq. 28 ) tb23FG(nDataIndex) = 0.  
          if ( channelval .eq. 29 ) tb31FG(nDataIndex) = 0. 
          if ( channelval .eq. 30 ) tb50FG(nDataIndex) = 0. 
          if ( channelval .eq. 32 ) tb53FG(nDataIndex) = 0. 
          if ( channelval .eq. 42 ) tb89FG(nDataIndex) = 0. 
        end if
      end do
    end do

  end subroutine extractParamForGrodyRun

  !--------------------------------------------------------------------------
  ! extractParamForBennartzRun
  !--------------------------------------------------------------------------  
  subroutine extractParamForBennartzRun(KCANO, ptbo, ptbomp, ptbcor, KNT, KNO, &
                                        tb89, tb150, tb1831, tb1832, tb1833, &
                                        tb89FG, tb150FG)

    !:Purpose: Extract Parameters required to run bennaertz for required channels:
    !          extract required channels:        
    !          89 Ghz = AMSU-B 1 = channel #43
    !         150 Ghz = AMSU-B 2 = channel #44

    implicit none
    ! Arguments
    integer,     intent(in)               :: KCANO(KNO,KNT)         ! observations channels
    real,        intent(in)               :: ptbo(KNO,KNT)          ! radiances
    real,        intent(in)               :: ptbomp(KNO,KNT)        ! radiances o-p
    real,        intent(in)               :: ptbcor(KNO,KNT)        ! correction aux radiances
    integer,     intent(in)               :: KNO                    ! nombre de canaux des observations 
    integer,     intent(in)               :: KNT                    ! nombre de tovs
    real,        intent(out)              :: tb89(KNT)              ! radiance frequence 89 Ghz  
    real,        intent(out)              :: tb150(KNT)             ! radiance frequence 150 Ghz  
    real,        intent(out)              :: tb1831(KNT)            ! radiance frequence ? Ghz  
    real,        intent(out)              :: tb1832(KNT)            ! radiance frequence ? Ghz  
    real,        intent(out)              :: tb1833(KNT)            ! radiance frequence ? Ghz  
    real,        intent(out)              :: tb89FG(KNT)            ! radiance frequence 89 Ghz  
    real,        intent(out)              :: tb150FG(KNT)           ! radiance frequence 150 Ghz  

    ! Locals
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: channelval

    do nDataIndex=1,KNT
      do nChannelIndex=1,KNO
        channelval = KCANO(nChannelIndex,nDataIndex)
        if ( ptbo(nChannelIndex,nDataIndex) .ne. mwbg_realMissing ) then
          if ( ptbcor(nChannelIndex,nDataIndex) .ne. mwbg_realMissing ) then
            if ( channelval .eq. 43 ) tb89(nDataIndex) = ptbo(nChannelIndex,nDataIndex) &
                 - ptbcor(nChannelIndex,nDataIndex)
            if ( channelval .eq. 44 ) tb150(nDataIndex) = ptbo(nChannelIndex,nDataIndex) &
                 - ptbcor(nChannelIndex,nDataIndex)
            if ( channelval .eq. 45 ) tb1831(nDataIndex) = ptbo(nChannelIndex,nDataIndex) &
                 - ptbcor(nChannelIndex,nDataIndex)
            if ( channelval .eq. 46 ) tb1832(nDataIndex) = ptbo(nChannelIndex,nDataIndex) &
                 - ptbcor(nChannelIndex,nDataIndex)
            if ( channelval .eq. 47 ) tb1833(nDataIndex) = ptbo(nChannelIndex,nDataIndex) &
                 - ptbcor(nChannelIndex,nDataIndex)
          else
            if ( channelval .eq. 43 ) tb89(nDataIndex) = ptbo(nChannelIndex,nDataIndex)
            if ( channelval .eq. 44 ) tb150(nDataIndex) = ptbo(nChannelIndex,nDataIndex)
            if ( channelval .eq. 45 ) tb1831(nDataIndex) = ptbo(nChannelIndex,nDataIndex)
            if ( channelval .eq. 46 ) tb1832(nDataIndex) = ptbo(nChannelIndex,nDataIndex)
            if ( channelval .eq. 47 ) tb1833(nDataIndex) = ptbo(nChannelIndex,nDataIndex)
          end if

          if ( channelval .eq. 43 ) tb89FG(nDataIndex)  = ptbo(nChannelIndex,nDataIndex) - &
                                                          ptbomp(nChannelIndex,nDataIndex)
          if ( channelval .eq. 44 ) tb150FG(nDataIndex) = ptbo(nChannelIndex,nDataIndex) - &
                                                          ptbomp(nChannelIndex,nDataIndex)
          
        else
          if ( channelval .eq. 43 ) tb89(nDataIndex) = 0.
          if ( channelval .eq. 44 ) tb150(nDataIndex) = 0.
          if ( channelval .eq. 45 ) tb1831(nDataIndex) = 0.
          if ( channelval .eq. 46 ) tb1832(nDataIndex) = 0.
          if ( channelval .eq. 47 ) tb1833(nDataIndex) = 0.

          if ( channelval .eq. 43 ) tb89FG(nDataIndex) = 0.
          if ( channelval .eq. 44 ) tb150FG(nDataIndex) = 0.
        end if
      end do
    end do

  end subroutine extractParamForBennartzRun

  !--------------------------------------------------------------------------
  ! amsuABTest10RttovRejectCheck
  !--------------------------------------------------------------------------
  subroutine amsuABTest10RttovRejectCheck(KCANO, KNOSAT, KNO, KNT, RESETQC, &
                                          STNID, KMARQ, ICHECK)

    !:Purpose:               10) test 10: RTTOV reject check (single)
    !                        Rejected datum flag has bit #9 on.

    implicit none
    ! Arguments
    integer,     intent(in)                :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)                :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)                :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)                :: KNT                            ! nombre de tovs
    logical,     intent(in)                :: RESETQC                        ! yes or not reset QC flag
    character *9, intent(in)               :: STNID                          ! identificateur du satellite
    integer,     intent(inout)             :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)             :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    ! Locals
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

  end subroutine amsuABTest10RttovRejectCheck

  !--------------------------------------------------------------------------
  ! amsuABTest1TopographyCheck
  !--------------------------------------------------------------------------
  subroutine amsuABTest1TopographyCheck(KCANO, KNOSAT, KNO, KNT, STNID, MTINTRP, &
                                        channelForTopoFilter, altitudeForTopoFilter, KMARQ, &
                                        ICHECK)

    !:Purpose:               1) test 1: Topography check (partial)
    !                        Channel 6 is rejected for topography >  250m.
    !                        Channel 7 is rejected for topography > 2000m.

    implicit none
    ! Arguments
    integer,     intent(in)                :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)                :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)                :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)                :: KNT                            ! nombre de tovs
    character *9,intent(in)                :: STNID                          ! identificateur du satellite
    integer,     intent(in)                :: channelForTopoFilter(:)        ! channel list for filter
    real,        intent(in)                :: altitudeForTopoFilter(:)       ! altitude threshold
    real,        intent(in)                :: MTINTRP(KNT)                   ! topo aux point d'obs
    integer,     intent(inout)             :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)             :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    ! Locals
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: numFilteringTest
    integer                                :: indexFilteringTest
    integer                                :: testIndex

    testIndex = 1

    !check consistency between channelForTopoFilter and altitudeForTopoFilter
    if ( size(altitudeForTopoFilter) /= size(channelForTopoFilter) ) then 
      call utl_abort('ABORT: amsuABTest1TopographyCheck, no consistency between channel List and altitude list ')
    end if 
   
    numFilteringTest =  size(altitudeForTopoFilter) 
    indexFilteringTest = 1

    do while ( indexFilteringTest <= numFilteringTest )
      do nDataIndex=1,KNT
        do nChannelIndex=1,KNO
          if ( KCANO(nChannelIndex,nDataIndex) == channelForTopoFilter(indexFilteringTest) ) then
            if ( MTINTRP(nDataIndex) >= altitudeForTopoFilter(indexFilteringTest)  ) then
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
      indexFilteringTest = indexFilteringTest + 1
    end do

  end subroutine amsuABTest1TopographyCheck

  !--------------------------------------------------------------------------
  ! amsuABTest2LandSeaQualifierCheck 
  !--------------------------------------------------------------------------
  subroutine amsuABTest2LandSeaQualifierCheck(KCANO, KNOSAT, KNO, KNT, STNID, &
                                              KTERMER, KMARQ, ICHECK)

    !:Purpose:                      2) test 2: "Land/sea qualifier" code check (full)
    !                                  allowed values are: 0, land,
    !                                  1, sea,
    !                                  2, coast.

    implicit none
    ! Arguments
    integer,     intent(in)                :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)                :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)                :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)                :: KNT                            ! nombre de tovs
    character *9,intent(in)                :: STNID                          ! identificateur du satellite
    integer,     intent(in)                :: KTERMER(KNT)                   ! land sea qualifier
    integer,     intent(inout)             :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)             :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    ! Locals
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex

  
    testIndex = 2
    do nDataIndex=1,KNT
      if ( KTERMER(nDataIndex) <  0  .or. &
          KTERMER(nDataIndex) >  2        ) then
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

  end subroutine amsuABTest2LandSeaQualifierCheck

  !--------------------------------------------------------------------------
  !  amsuABTest3TerrainTypeCheck
  !--------------------------------------------------------------------------
  subroutine amsuABTest3TerrainTypeCheck(KCANO, KNOSAT, KNO, KNT, STNID, &
                                         ITERRAIN, KMARQ, ICHECK)

    !:Purpose:                     3) test 3: "Terrain type" code check (full)
    !                                 allowed values are: -1, missing,
    !                                 0, sea-ice,
    !                                 1, snow on land.

    implicit none
    ! Arguments
    integer,     intent(in)               :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)               :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)               :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)               :: KNT                            ! nombre de tovs
    character *9,intent(in)               :: STNID                          ! identificateur du satellite
    integer,     intent(in)               :: ITERRAIN(KNT)                  ! terrain type
    integer,     intent(inout)            :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)            :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    ! Locals
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex

 
    testIndex = 3
    do nDataIndex=1,KNT
      if ( ITERRAIN(nDataIndex) /= mwbg_intMissing ) then
        if ( ITERRAIN(nDataIndex) <  0  .or. &
           ITERRAIN(nDataIndex) >  1        ) then
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

  end subroutine amsuABTest3TerrainTypeCheck

  !--------------------------------------------------------------------------
  ! amsuABTest4FieldOfViewCheck 
  !--------------------------------------------------------------------------
  subroutine amsuABTest4FieldOfViewCheck(KCANO, KNOSAT, KNO, KNT, STNID, ISCNPOS, &
                                         maxScanAngleAMSU, KMARQ, ICHECK)

    !:Purpose:                          4) test 4: Field of view number check (full)
    !                                      Field of view acceptable range is [1,maxScanAngleAMSU]  for AMSU footprints.
    implicit none
    ! Arguments
    integer,     intent(in)               :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)               :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)               :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)               :: KNT                            ! nombre de tovs
    character *9,intent(in)               :: STNID                          ! identificateur du satellite
    integer,     intent(in)               :: ISCNPOS(KNT)                   ! position sur le "scan" 
    integer,     intent(in)               :: maxScanAngleAMSU               ! max scan angle 
    integer,     intent(inout)            :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)            :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    ! Locals
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex

    testIndex = 4
    do nDataIndex=1,KNT
      do nChannelIndex=1,KNO
        if ( ISCNPOS(nDataIndex) < 1 .or. &
            ISCNPOS(nDataIndex) > maxScanAngleAMSU ) then
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
  
  end subroutine amsuABTest4FieldOfViewCheck 
  
  !--------------------------------------------------------------------------
  ! amsuABTest5ZenithAngleCheck 
  !--------------------------------------------------------------------------
  subroutine amsuABTest5ZenithAngleCheck(KCANO, KNOSAT, KNO, KNT, STNID, &
                                         SATZEN, KMARQ, ICHECK)

    !:Purpose:                   5) test 5: Satellite zenith angle check (full)
    !                               Satellite zenith angle acceptable range is [0.,60.].
    implicit none
    ! Arguments
    integer,     intent(in)               :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)               :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)               :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)               :: KNT                            ! nombre de tovs
    character *9,intent(in)               :: STNID                          ! identificateur du satellite
    real,        intent(in)               :: SATZEN(KNT)                    ! satellite zenith angle 
    integer,     intent(inout)            :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)            :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    ! Locals
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex


    testIndex = 5
    do nDataIndex=1,KNT
      if ( SATZEN(nDataIndex) /= mwbg_realMissing ) then
        if ( SATZEN(nDataIndex) <  0.  .or. &
           SATZEN(nDataIndex) > 60.       ) then
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

  end subroutine amsuABTest5ZenithAngleCheck 

  !--------------------------------------------------------------------------
  ! amsuABTest6ZenAngleAndFovConsistencyCheck 
  !--------------------------------------------------------------------------
  subroutine amsuABTest6ZenAngleAndFovConsistencyCheck(KCANO, KNOSAT, KNO, KNT, STNID, SATZEN,  ZANGL, &
                                                  ISCNPOS, maxScanAngleAMSU, KMARQ, ICHECK)

    !:Purpose:                            6) test 6: "Sat. zenith angle"/"field of view" consistency check.  (full)
    !                                        Acceptable difference between "Satellite zenith angle"  and
    !                                       "approximate angle computed from field of view number" is 1.8 degrees.

    implicit none
    ! Arguments
    integer,     intent(in)               :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)               :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)               :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)               :: KNT                            ! nombre de tovs
    character *9,intent(in)               :: STNID                          ! identificateur du satellite
    real,        intent(in)               :: SATZEN(KNT)                    ! satellite zenith angle 
    real,        intent(in)               :: ZANGL                          ! satellite constant param
    integer,     intent(in)               :: ISCNPOS(KNT)                   ! position sur le "scan" 
    integer,     intent(in)               :: maxScanAngleAMSU               ! max scan angle 
    integer,     intent(inout)            :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)            :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    ! Locals
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex
    real                                  :: APPROXIM 
    real                                  :: ANGDif 

    testIndex = 6
    do nDataIndex=1,KNT
      if ( SATZEN (nDataIndex) /=  mwbg_realMissing   .and. &
           ISCNPOS(nDataIndex) /=  mwbg_intMissing  ) then
        APPROXIM = ABS((ISCNPOS(nDataIndex)-maxScanAngleAMSU/2.-0.5)*ZANGL)
        ANGDif = ABS(SATZEN (nDataIndex)-APPROXIM)
        if ( ANGDif > 1.8 ) then 
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

  end subroutine amsuABTest6ZenAngleAndFovConsistencyCheck

  !--------------------------------------------------------------------------
  !  amsuABTest7landSeaQualifyerAndModelLandSeaConsistencyCheck
  !--------------------------------------------------------------------------
  subroutine amsuABTest7landSeaQualifyerAndModelLandSeaConsistencyCheck(KCANO, KNOSAT, KNO, KNT, STNID, MGINTRP, KTERMER, &
                                                                        KMARQ, ICHECK)

    !:Purpose:                    7) test 7: "Land/sea qual."/"model land/sea" consistency check.    (full)
    !                                Acceptable conditions are:
    !                                - both over ocean (ktermer=1; mg<0.01), new threshold 0.20, jh dec 2000,
    !                                - both over land  (ktermer=0; mg>0.80), new threshold 0.50, jh dec 2000.
    !                                - Other conditions are unacceptable.


    implicit none
    ! Arguments
    integer,     intent(in)                :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)                :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)                :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)                :: KNT                            ! nombre de tovs
    character *9,intent(in)                :: STNID                          ! identificateur du satellite
    real,        intent(in)                :: MGINTRP(KNT)                   ! glace mer 
    integer,     intent(in)                :: KTERMER(KNT)                   ! land sea qualifyer 
    integer,     intent(inout)             :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)             :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    ! Locals
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex

    testIndex = 7
    do nDataIndex=1,KNT
      if ( KTERMER (nDataIndex) /= mwbg_intMissing  ) then
        if     ( KTERMER(nDataIndex) == 1       .and. &
                MGINTRP(nDataIndex) < 0.20          ) then
        elseif ( KTERMER(nDataIndex) == 0       .and. &
                MGINTRP(nDataIndex) > 0.50          ) then
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

  end subroutine amsuABTest7landSeaQualifyerAndModelLandSeaConsistencyCheck 

  !--------------------------------------------------------------------------
  !  amsuABTest9UncorrectedTbCheck
  !--------------------------------------------------------------------------
  subroutine amsuABTest9UncorrectedTbCheck(KCANO, KNOSAT, KNO, KNT, STNID, &
                                           RESETQC, KMARQ, ICHECK)

    !:Purpose:                  9) test 9: Uncorrected Tb check (single)
    !                              Uncorrected datum (flag bit #6 off). In this case switch bit 11 ON.

    implicit none
    ! Arguments
    integer,     intent(in)               :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)               :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)               :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)               :: KNT                            ! nombre de tovs
    character *9,intent(in)               :: STNID                          ! identificateur du satellite
    logical,     intent(in)               :: RESETQC                        ! yes or not reset QC flag
    integer,     intent(inout)            :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)            :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    ! Locals
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex
    integer                               :: IBIT


    if (.not. RESETQC) then
      testIndex = 9
      do nDataIndex=1,KNT
        do nChannelIndex=1,KNO
          if ( KCANO(nChannelIndex,nDataIndex) /= 20 ) then
            IBIT = AND(KMARQ(nChannelIndex,nDataIndex), 2**6)
            if ( IBIT == 0  ) then
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

  end subroutine amsuABTest9UncorrectedTbCheck
 
  !--------------------------------------------------------------------------
  ! amsuABTest11RadianceGrossValueCheck
  !--------------------------------------------------------------------------
  subroutine amsuABTest11RadianceGrossValueCheck(KCANO, KNOSAT, KNO, KNT, STNID, PTBO, GROSSMIN, &
                                                 GROSSMAX, KMARQ, ICHECK)

    !:Purpose:                     11) test 11: Radiance observation "Gross" check (single) 
    !                                               Change this test from full to single. jh nov 2000.

    implicit none
    ! Arguments
    integer,     intent(in)               :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)               :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)               :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)               :: KNT                            ! nombre de tovs
    character *9,intent(in)               :: STNID                          ! identificateur du satellite
    real,        intent(in)               :: PTBO(KNO,KNT)                  ! radiances 
    real,        intent(in)               :: GROSSMIN(:)                  ! Gross val min 
    real,        intent(in)               :: GROSSMAX(:)                  ! Gross val max 
    integer,     intent(inout)            :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)            :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    ! Locals
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex
    logical                               :: GROSSERROR

    testIndex = 11
    do nDataIndex=1,KNT
      GROSSERROR = .FALSE.
      do nChannelIndex=1,KNO
        if ( KCANO(nChannelIndex,nDataIndex) /= 20     .and. &
            KCANO(nChannelIndex,nDataIndex) >=  1     .and. &
            KCANO(nChannelIndex,nDataIndex) <=  KNO       ) then  
          if ( PTBO(nChannelIndex,nDataIndex) /= mwbg_realMissing .and. &
             ( PTBO(nChannelIndex,nDataIndex) < GROSSMIN(KCANO(nChannelIndex,nDataIndex)).or. &
               PTBO(nChannelIndex,nDataIndex) > GROSSMAX(KCANO(nChannelIndex,nDataIndex))     ) ) then
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

  end subroutine amsuABTest11RadianceGrossValueCheck 
  
  !--------------------------------------------------------------------------
  ! amsuaTest12GrodyClwCheck
  !--------------------------------------------------------------------------
  subroutine amsuaTest12GrodyClwCheck(KCANO, KNOSAT, KNO, KNT, STNID, clwObs, &
                                      clwFG, useStateDepSigmaObs, ktermer, MISGRODY, MXCLWREJ, &
                                      ICLWREJ, KMARQ, ICHECK)

    !:Purpose:                    12) test 12: Grody cloud liquid water check (partial)
    !                                 For Cloud Liquid Water > clwQcThreshold, reject AMSUA-A channels 1-5 and 15.

    implicit none
    ! Arguments
    integer,     intent(in)               :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)               :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)               :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)               :: KNT                            ! nombre de tovs
    character *9,intent(in)               :: STNID                          ! identificateur du satellite
    real,        intent(in)               :: clwObs(KNT)                    ! retrieved cloud liquid water from observation
    real,        intent(in)               :: clwFG(KNT)                     ! retrieved cloud liquid water from background
    logical,     intent(in)               :: useStateDepSigmaObs(:,:)       ! if using state dependent obs error
    integer,     intent(in)               :: KTERMER(KNT)                   ! land sea qualifyer 
    real,        intent(in)               :: MISGRODY                       ! MISGRODY
    integer,     intent(in)               :: MXCLWREJ                       ! cst 
    integer,     intent(in)               :: ICLWREJ(MXCLWREJ)              !
    integer,     intent(inout)            :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)            :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    ! Locals
    integer                               :: channelval
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex
    integer                               :: INDXCAN 
    real                                  :: clwUsedForQC
    real                                  :: clwObsFGaveraged
    logical                               :: surfTypeIsWater 

    testIndex = 12
    do nDataIndex=1,KNT
      if ( tvs_mwAllskyAssim ) then
        clwObsFGaveraged = 0.5 * (clwObs(nDataIndex) + clwFG(nDataIndex))
        clwUsedForQC = clwObsFGaveraged
      else
        clwUsedForQC = clwObs(nDataIndex)
      end if

      surfTypeIsWater = ( ktermer(nDataIndex) ==  1 )

      if ( clwUsedForQC /=  MISGRODY  ) then
        if ( clwUsedForQC > mwbg_clwQcThreshold ) then
          do nChannelIndex=1,KNO
            INDXCAN = ISRCHEQI (ICLWREJ,MXCLWREJ,KCANO(nChannelIndex,nDataIndex))
            if ( INDXCAN /= 0 )  then
              ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**7)
              rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                       rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
            end if
          end do
          if ( mwbg_debug ) then
            write(*,*)STNID(2:9),'Grody cloud liquid water check', &
                      ' REJECT. CLW= ',clwUsedForQC, ' SEUIL= ',mwbg_clwQcThreshold
          end if
        end if

        ! In all-sky mode, turn on bit=23 for cloud-affected radiances 
        ! when there is mismatch between clwObs and clwFG
        ! (to be used in gen_bias_corr)
        clwObsFGaveraged = 0.5 * (clwObs(nDataIndex) + clwFG(nDataIndex))
        IF ( tvs_mwAllskyAssim .and. &
            (clwObsFGaveraged > mwbg_cloudyClwThresholdBcorr .or. &
            clwObsFGaveraged == MISGRODY) ) then
          do nChannelIndex = 1,KNO
            INDXCAN = ISRCHEQI(ICLWREJ,MXCLWREJ,KCANO(nChannelIndex,nDataIndex))
            if ( INDXCAN /= 0 ) KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**23)
          end do
          if ( mwbg_debug ) then
            write(*,*) STNID(2:9),' Grody cloud liquid water check', &
                      ' cloud-affected obs. CLW= ',clwUsedForQC, ', threshold= ', &
                      mwbg_cloudyClwThresholdBcorr
          end if
        end if

      ! Reject surface sensitive observations over water, in all-sky mode, 
      ! if CLW is not retrieved, and is needed to define obs error.
      else if ( tvs_mwAllskyAssim .and. surfTypeIsWater .and. &
                clwUsedForQC == MISGRODY ) then

        loopChannel: do nChannelIndex = 1, KNO
          channelval = KCANO(nChannelIndex,nDataIndex)
          INDXCAN = ISRCHEQI(ICLWREJ,MXCLWREJ,channelval)
          if ( INDXCAN /= 0 .and. useStateDepSigmaObs(channelval,knosat) ) then
            ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**7)
            rejectionCodArray(testIndex,channelval,KNOSAT) = &
                     rejectionCodArray(testIndex,channelval,KNOSAT)+ 1
          end if
        end do loopChannel

      end if
    end do

  end subroutine amsuaTest12GrodyClwCheck 

  !-------------------------------------------------------------------------
  ! amsubTest12DrynessIndexCheck
  !-------------------------------------------------------------------------
  subroutine amsubTest12DrynessIndexCheck(KCANO, KNOSAT, KNO, KNT, STNID, tb1831, tb1833, &
                                          ktermer, glintrp, KMARQ, ICHECK)

    !:Purpose:  12) test 12: Dryness index check
    !           The difference between channels AMSUB-3 and AMSUB-5 is used as an indicator
    !           of "dryness" of the atmosphere. In extreme dry conditions, channels AMSUB-3 4 and 5
    !           are sensitive to the surface.
    !           Therefore, various thresholds are used to reject channels AMSUB-3 4 and 5 over land and ice

    implicit none
    ! Arguments
    integer,     intent(in)               :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)               :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)               :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)               :: KNT                            ! nombre de tovs
    character *9,intent(in)               :: STNID                          ! identificateur du satellite
    real,        intent(in)               :: tb1831(KNT)                    ! tb for channel  
    real,        intent(in)               :: tb1833(KNT)                    ! tb for channel  
    integer,     intent(in)               :: ktermer(:)                     ! mask terre-mer
    real,        intent(in)               :: glintrp(:)                     ! topo interpolated to obs point
    integer,     intent(inout)            :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)            :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    ! Locals
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex
    real                                  :: drynessIndex

    testIndex = 12
    do nDataIndex = 1, KNT
      drynessIndex = tb1831(nDataIndex) - tb1833(nDataIndex)
      do nChannelIndex = 1, KNO
        if ( .not. ((ktermer (nDataIndex) == 1) .and. &
                    (glintrp (nDataIndex) < 0.01)) ) then
          if ( KCANO(nChannelIndex,nDataIndex) == 45  .and. &
               drynessIndex > 0.) then
            ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**7)
            rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
            rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
            if (mwbg_debug) then
              WRITE(6,*)STNID(2:9),' DRYNESS INDEX REJECT.',        &
                       'CHANNEL=', KCANO(nChannelIndex,nDataIndex), &
                       'INDEX= ',drynessIndex
            end if
          else if ( (KCANO(nChannelIndex,nDataIndex) == 46 )  .and. &
                    drynessIndex > -10.) then
            ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**7)
            rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) =  &
            rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
            if (mwbg_debug) then
              WRITE(6,*)STNID(2:9),' DRYNESS INDEX REJECT.',       &
                      'CHANNEL=', KCANO(nChannelIndex,nDataIndex),&
                      'INDEX= ',drynessIndex
            end if
          else if ( (KCANO(nChannelIndex,nDataIndex) == 47)  .and. &
                    drynessIndex > -20.) then
            ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**7)
            rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
            rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
            if (mwbg_debug) then
              WRITE(6,*)STNID(2:9),' DRYNESS INDEX REJECT.',       &
                       'CHANNEL=', KCANO(nChannelIndex,nDataIndex),&
                       'INDEX= ',drynessIndex
            end if
          end if
        end if
      end do
    end do

  end subroutine amsubTest12DrynessIndexCheck

  !--------------------------------------------------------------------------
  ! amsuaTest13GrodyScatteringIndexCheck
  !--------------------------------------------------------------------------
  subroutine amsuaTest13GrodyScatteringIndexCheck(KCANO, KNOSAT, KNO, KNT, STNID, scatw, KTERMER, ITERRAIN, &
                                                  MISGRODY, MXSCATREJ, ISCATREJ, KMARQ, ICHECK)

    !:Purpose:                  13) test 13: Grody scattering index check (partial)
    !                               For Scattering Index > 9, reject AMSUA-A channels 1-6 and 15.

    implicit none
    ! Arguments
    integer,     intent(in)                :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)                :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)                :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)                :: KNT                            ! nombre de tovs
    character *9,intent(in)                :: STNID                          ! identificateur du satellite
    real,        intent(in)                :: scatw(KNT)                     ! scattering index 
    integer,     intent(in)                :: KTERMER(KNT)                   ! land sea qualifyer 
    integer,     intent(in)                :: ITERRAIN(KNT)                  ! terrain type 
    real,        intent(in)                :: MISGRODY                       ! MISGRODY
    integer,     intent(in)                :: MXSCATREJ                       ! cst 
    integer,     intent(in)                :: ISCATREJ(MXSCATREJ)              !
    integer,     intent(inout)             :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)             :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    ! Locals
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex
    integer                                :: INDXCAN 
    real                                   :: ZSEUILSCAT

    testIndex = 13
    ZSEUILSCAT = 9.0
    do nDataIndex=1,KNT
      if ( SCATW(nDataIndex) /=  MISGRODY  ) then
        if (  KTERMER (nDataIndex) ==  1 .and. &
             ITERRAIN(nDataIndex) /=  0 .and. &   
             SCATW   (nDataIndex) > ZSEUILSCAT   ) then
          do nChannelIndex=1,KNO
            INDXCAN = ISRCHEQI (ISCATREJ,MXSCATREJ,KCANO(nChannelIndex,nDataIndex))
            if ( INDXCAN /= 0 )  then
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
  ! amsubTest13BennartzScatteringIndexCheck
  !--------------------------------------------------------------------------
  subroutine amsubTest13BennartzScatteringIndexCheck(KCANO, KNOSAT, KNT, KNO, STNID, scatwObs, scatwFG, scatl, &
                                                     useStateDepSigmaObs, KTERMER, GLINTRP, KMARQ, ICHECK)

    !:Purpose:                  13) test 13: Bennartz scattering index check (full)
    !                               For Scattering Index > 40 sea ice
    !                                                    > 15 sea
    !                                                    > 0 land reject all AMSUB Channels
    !

    implicit none
    ! Arguments
    integer,     intent(in)                :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)                :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)                :: KNT                            ! nombre de tovs
    integer,     intent(in)                :: KNO                            ! nombre canaux de tovs
    character *9,intent(in)                :: STNID                          ! identificateur du satellite
    real,        intent(in)                :: scatwObs(KNT)                  ! scattering index over water from observation
    real,        intent(in)                :: scatwFG(KNT)                   ! scattering index over water from background
    real,        intent(in)                :: scatl(KNT)                     ! scattering index over land
    logical,     intent(in)                :: useStateDepSigmaObs(:,:)       ! if using state dependent obs error
    integer,     intent(in)                :: KTERMER(KNT)                   ! land sea qualifyer 
    real,        intent(in)                :: GLINTRP(KNT)                   ! glace de mer
    integer,     intent(inout)             :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)             :: ICHECK(KNO,KNT)                ! indicateur du QC par canal

    ! Locals
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex
    real                                   :: ZSEUILSCATICE
    real                                   :: ZSEUILSCATL
    real                                   :: ZSEUILSCATW
    real                                   :: siUsedForQcThreshold
    integer                                :: channelval
    real                                   :: siObsFGaveraged
    real                                   :: siUsedForQC
    logical                                :: FULLREJCT
    logical                                :: surfTypeIsSea 

    testIndex = 13

    ZSEUILSCATICE = 40.0
    ZSEUILSCATW   = 15.0
    ZSEUILSCATL   =  0.0
    do nDataIndex=1,KNT
      FULLREJCT = .FALSE.
      surfTypeIsSea = .false.

      if (  KTERMER (nDataIndex) == 1  ) then
        if ( GLINTRP (nDataIndex) > 0.01 ) then ! sea ice 
          if (  scatwObs(nDataIndex) /= mwbg_realMissing    .and. &
                scatwObs(nDataIndex) > ZSEUILSCATICE  ) then
            FULLREJCT = .TRUE.
          end if

        else                                    ! sea 
          surfTypeIsSea = .true.

          if ( tvs_mwAllskyAssim ) then
            siObsFGaveraged = 0.5 * (scatwObs(nDataIndex) + scatwFG(nDataIndex))
            siUsedForQC = siObsFGaveraged
            siUsedForQcThreshold = mwbg_siQcOverWaterThreshold
          else
            siUsedForQC = scatwObs(nDataIndex)
            siUsedForQcThreshold = ZSEUILSCATW
          end if

          if (  siUsedForQC /= mwbg_realMissing    .and. &
                siUsedForQC > siUsedForQcThreshold  ) then
            FULLREJCT = .TRUE.
          end if
        end if

      else                                      ! land
        if (  SCATL(nDataIndex) /= mwbg_realMissing    .and. &
              SCATL(nDataIndex) > ZSEUILSCATL  ) then
          FULLREJCT = .TRUE.
        end if
      end if
      if ( FULLREJCT )  then
        do nChannelIndex=1,KNO
          ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
          KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
          KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**7)
          rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
          rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) + 1
        end do
        if (mwbg_debug) then
          write(*,*)  STNID(2:9), ' BENNARTZ scattering index check REJECT, scatwObs=', &
                      scatwObs(nDataIndex), ', scatwFG=', scatwFG(nDataIndex), &
                     ', SCATL= ',SCATL(nDataIndex)
        end if
      end if


      if (tvs_mwAllskyAssim .and. surfTypeIsSea) then
        siObsFGaveraged = 0.5 * (scatwObs(nDataIndex) + scatwFG(nDataIndex))

        ! In all-sky mode, turn on bit=23 for cloud-affected radiances over sea
        ! when there is mismatch between scatwObs and scatwFG
        ! (to be used in gen_bias_corr)
        if (siObsFGaveraged > mwbg_cloudySiThresholdBcorr .or. &
            siObsFGaveraged == mwbg_realMissing) then
          do nChannelIndex = 1,KNO
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**23)
          end do
          if ( mwbg_debug ) then
            write(*,*) STNID(2:9),' BENNARTZ scattering index check', &
                      ' cloud-affected obs. siObsFGaveraged= ', siObsFGaveraged, ', threshold= ', &
                      mwbg_cloudySiThresholdBcorr
          end if
        end if
      end if
      
      if (tvs_mwAllskyAssim .and. surfTypeIsSea) then
        siObsFGaveraged = 0.5 * (scatwObs(nDataIndex) + scatwFG(nDataIndex))

        ! In all-sky mode, reject observations over sea if siObsFGaveraged can not be computed.
        ! siObsFGaveraged is needed to define obs error.
        if (siObsFGaveraged == mwbg_realMissing) then

          loopChannel3: do nChannelIndex = 1, KNO
            channelval = KCANO(nChannelIndex,nDataIndex)
            if (useStateDepSigmaObs(channelval,knosat)) then
              ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**7)
              rejectionCodArray(testIndex,channelval,KNOSAT) = &
                      rejectionCodArray(testIndex,channelval,KNOSAT)+ 1
            end if
          end do loopChannel3
        end if ! if (siObsFGaveraged == mwbg_realMissing)

      end if ! if (tvs_mwAllskyAssim .and. surfTypeIsSea)
      
    end do !do nDataIndex=1,KNT

  end subroutine amsubTest13BennartzScatteringIndexCheck

  !--------------------------------------------------------------------------
  ! amsuaTest14RogueCheck
  !--------------------------------------------------------------------------
  subroutine amsuaTest14RogueCheck(KCANO, KNOSAT, KNO, KNT, STNID, ROGUEFAC, TOVERRST, clwThreshArr, &
                                   useStateDepSigmaObs, sigmaObsErr, ktermer, PTBOMP, clwObs, clwFG, &
                                   MISGRODY, MXSFCREJ, ISFCREJ, KMARQ, ICHECK)

    !:Purpose:                     14) test 14: "Rogue check" for (O-P) Tb residuals out of range.
    !                                  (single/full). Les observations, dont le residu (O-P) 
    !                                  depasse par un facteur (roguefac) l'erreur totale des TOVS.
    !                                  N.B.: a reject by any of the 3 surface channels produces the 
    !                                  rejection of AMSUA-A channels 1-5 and 15.

    implicit none
    ! Arguments
    integer,     intent(in)                :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)                :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)                :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)                :: KNT                            ! nombre de tovs
    character *9,intent(in)                :: STNID                          ! identificateur du satellite
    real,        intent(in)                :: ROGUEFAC(:)                    ! rogue factor 
    real(8),     intent(in)                :: TOVERRST(:,:)                  !  erreur totale TOVs
    logical,     intent(in)                :: useStateDepSigmaObs(:,:)       ! if using state dependent obs error
    real(8),     intent(in)                :: clwThreshArr(:,:,:)            ! cloud threshold array
    real(8),     intent(in)                :: sigmaObsErr(:,:,:)             ! sigma obs error
    integer,     intent(in)                :: ktermer(KNT)                   ! land/sea identifier
    real,        intent(in)                :: clwObs(KNT)                    ! retrieved cloud liquid water from observation
    real,        intent(in)                :: clwFG(KNT)                     ! retrieved cloud liquid water from background
    real,        intent(in)                :: PTBOMP(KNO,KNT)                ! radiance o-p 
    real,        intent(in)                :: MISGRODY                       ! MISGRODY
    integer,     intent(in)                :: MXSFCREJ                       ! cst 
    integer,     intent(in)                :: ISFCREJ(MXSFCREJ)              !
    integer,     intent(inout)             :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)             :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    ! Locals
    integer                                :: channelval
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex
    integer                                :: INDXCAN 
    real                                   :: XCHECKVAL
    real                                   :: clwThresh1 
    real                                   :: clwThresh2
    real                                   :: sigmaThresh1 
    real                                   :: sigmaThresh2
    real                                   :: sigmaObsErrUsed  
    logical                                :: SFCREJCT
    logical                                :: surfTypeIsWater
    real                                   :: clwObsFGaveraged 

    testIndex = 14
    do nDataIndex=1,KNT
      surfTypeIsWater = ( ktermer(nDataIndex) ==  1 )
      SFCREJCT = .FALSE.
      do nChannelIndex=1,KNO
        channelval = KCANO(nChannelIndex,nDataIndex)
        if ( channelval /= 20 ) then
          ! using state-dependent obs error only over water.
          ! obs over sea-ice will be rejected in test 15.
          if ( tvs_mwAllskyAssim .and. useStateDepSigmaObs(channelval,KNOSAT) &
                .and. surfTypeIsWater ) then
            clwThresh1 = clwThreshArr(channelval,KNOSAT,1)
            clwThresh2 = clwThreshArr(channelval,KNOSAT,2)
            sigmaThresh1 = sigmaObsErr(channelval,KNOSAT,1)
            sigmaThresh2 = sigmaObsErr(channelval,KNOSAT,2)
            clwObsFGaveraged = 0.5 * (clwObs(nDataIndex) + clwFG(nDataIndex))
            if ( clwObsFGaveraged == MISGRODY ) then
              sigmaObsErrUsed = MPC_missingValue_R4
            else
              sigmaObsErrUsed = calcStateDepObsErr_r4(clwThresh1,clwThresh2,sigmaThresh1, &
                                                        sigmaThresh2,clwObsFGaveraged)
            end if
          else
            sigmaObsErrUsed = TOVERRST(channelval,KNOSAT)
          end if
          ! For sigmaObsErrUsed=MPC_missingValue_R4 (clwObsFGaveraged=MISGRODY
          ! in all-sky mode), the observation is already rejected in test 12.
          XCHECKVAL = ROGUEFAC(channelval) * sigmaObsErrUsed
          if ( PTBOMP(nChannelIndex,nDataIndex) /= mwbg_realMissing .and. &
              ABS(PTBOMP(nChannelIndex,nDataIndex)) >= XCHECKVAL .and. &
              sigmaObsErrUsed /= MPC_missingValue_R4 ) then
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
            if ( channelval == 28 .or. &
                 channelval == 29 .or. &
                 channelval == 30      ) then
              SFCREJCT = .TRUE.
            end if
          end if
        end if
      end do

      if ( SFCREJCT ) then
        do nChannelIndex=1,KNO
          INDXCAN = ISRCHEQI (ISFCREJ,MXSFCREJ,KCANO(nChannelIndex,nDataIndex))
          if ( INDXCAN /= 0 )  then
            if ( ICHECK(nChannelIndex,nDataIndex) /= testIndex ) then
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
  ! amsubTest14RogueCheck
  !--------------------------------------------------------------------------
  subroutine amsubTest14RogueCheck(KCANO, KNOSAT, KNO, KNT, STNID, ROGUEFAC, TOVERRST, &
                                   siThreshArr, sigmaObsErr, useStateDepSigmaObs, &
                                   iterrain, ktermer, PTBOMP, ICH2OMPREJ, MXCH2OMPREJ, KMARQ, ICHECK)

    !:Purpose:                     14) test 14: "Rogue check" for (O-P) Tb residuals out of range.
    !                                  (single)
    !                                 Also, remove CH2,3,4,5 if CH2 |O-P|>5K           (partial) 
    !                                  

    implicit none
    ! Arguments
    integer,     intent(in)                :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)                :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)                :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)                :: KNT                            ! nombre de tovs
    character *9,intent(in)                :: STNID                          ! identificateur du satellite
    real,        intent(in)                :: ROGUEFAC(:)                    ! rogue factor 
    real(8),     intent(in)                :: TOVERRST(:,:)                  ! erreur totale TOVs
    real(8),     intent(in)                :: siThreshArr(:,:,:)             ! SI threshold array
    real(8),     intent(in)                :: sigmaObsErr(:,:,:)             ! sigma obs error    
    logical,     intent(in)                :: useStateDepSigmaObs(:,:)       ! if using state dependent obs error
    integer,     intent(in)                :: ktermer(KNT)                   !
    integer,     intent(in)                :: iterrain(KNT)                  !
    real,        intent(in)                :: PTBOMP(KNO,KNT)                ! radiance o-p 
    integer,     intent(in)                :: MXCH2OMPREJ                !
    integer,     intent(in)                :: ICH2OMPREJ(MXCH2OMPREJ)                 !
    integer,     intent(inout)             :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)             :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    ! Locals
    integer                                :: channelval
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex
    integer                                :: INDXCAN 
    real                                   :: XCHECKVAL
    logical                                :: SFCREJCT
    logical                                :: CH2OMPREJCT

    testIndex = 14
    do nDataIndex=1,KNT
      CH2OMPREJCT = .FALSE.
      SFCREJCT = .FALSE.
      do nChannelIndex=1,KNO
        channelval = KCANO(nChannelIndex,nDataIndex)
        if ( channelval /= 20 ) then
          XCHECKVAL = ROGUEFAC(channelval) * TOVERRST(channelval,KNOSAT) 
          if ( PTBOMP(nChannelIndex,nDataIndex)     /=  mwbg_realMissing    .and. &
              ABS(PTBOMP(nChannelIndex,nDataIndex)) >=  XCHECKVAL     ) then
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
          end if
          if ( channelval == 44                                      .and. &
               PTBOMP(nChannelIndex,nDataIndex)  /= mwbg_realMissing .and. &
               ABS(PTBOMP(nChannelIndex,nDataIndex)) >= 5.0     ) then
            CH2OMPREJCT = .TRUE.
          end if
        end if
      end do

      if ( (CH2OMPREJCT) .and. (ktermer(nDataIndex) == 1) .and. (iterrain(nDataIndex) /= 0) ) then
        do nChannelIndex=1,KNO
          INDXCAN = ISRCHEQI (ICH2OMPREJ,MXCH2OMPREJ,KCANO(nChannelIndex,nDataIndex))
          if ( INDXCAN /= 0 )  then
            if ( ICHECK(nChannelIndex,nDataIndex) /= testIndex ) then
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
  end subroutine amsubTest14RogueCheck


  !--------------------------------------------------------------------------
  ! amsuABTest15ChannelSelectionWithIutilst
  !--------------------------------------------------------------------------
  subroutine amsuABTest15ChannelSelectionWithIutilst(KCANO, KNOSAT, KNO, KNT, STNID, KTERMER, ITERRAIN, &
                                              GLINTRP, IUTILST, MXSFCREJ2, ISFCREJ2, KMARQ, ICHECK)

    !:Purpose:                         ! 15) test 15: Channel Selection using array IUTILST(chan,sat)
    !                                    IUTILST = 0 (blacklisted)
    !                                    - 1 (assmilate)
    !                                    - 2 (assimilate over open water only)
    !
    !                                    - We also set QC flag bits 7 and 9 ON for channels with IUTILST=2 
    !                                    over land or sea-ice and we set QC flag bits 7 and 9 ON for channels
    !                                    1-3,15 over land or sea-ice REGARDLESS of IUTILST value 
    !                                    (but IUTILST=0 always for these unassimilated channels).


    implicit none
    ! Arguments
    integer,     intent(in)                :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)                :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)                :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)                :: KNT                            ! nombre de tovs
    character *9,intent(in)                :: STNID                          ! identificateur du satellite
    integer,     intent(in)                :: KTERMER(KNT)                   ! land sea qualifyer 
    integer,     intent(in)                :: ITERRAIN(KNT)                  ! terrain type
    real  ,      intent(in)                :: GLINTRP(KNT)                   ! gl
    integer,     intent(in)                :: IUTILST(:,:)          !  channel selection
    integer,     intent(in)                :: MXSFCREJ2                       ! cst 
    integer,     intent(in)                :: ISFCREJ2(MXSFCREJ2)              !
    integer,     intent(inout)             :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)             :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    ! Locals
    integer                                :: channelval
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex
    integer                                :: ITRN 
    integer                                :: INDXCAN
    logical                                :: SFCREJCT

    testIndex = 15

    do nDataIndex=1,KNT
      ITRN = ITERRAIN(nDataIndex)
      if ( KTERMER (nDataIndex) == 1    .and. &
           ITERRAIN(nDataIndex) == -1   .and. &
           GLINTRP (nDataIndex) >= 0.01       ) then
         ITRN = 0
      end if        
      do nChannelIndex=1,KNO
          channelval = KCANO(nChannelIndex,nDataIndex)
          INDXCAN = ISRCHEQI (ISFCREJ2,MXSFCREJ2,channelval)
          if ( INDXCAN /= 0 )  then
            if ( KTERMER (nDataIndex) == 0 .or. ITRN == 0 )  then
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**7)
            end if
          end if
          if ( IUTILST(channelval,KNOSAT) /= 1 ) then
            SFCREJCT = .FALSE.
            if ( IUTILST(channelval,KNOSAT) == 0 ) then
              SFCREJCT = .TRUE.
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**11)
            else 
              if ( KTERMER (nDataIndex) == 0 .or. ITRN == 0 )  then
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

  end subroutine amsuABTest15ChannelSelectionWithIutilst

  !--------------------------------------------------------------------------
  ! amsuaTest16ExcludeExtremeScattering
  !--------------------------------------------------------------------------
  subroutine amsuaTest16ExcludeExtremeScattering(KCANO, KNOSAT, KNO, KNT, STNID, KTERMER, &
                                                 PTBO, btClear2D, PTBOMP, KMARQ, ICHECK)
    !:Purpose: Exclude radiances affected extreme scattering in deep convective region.
    !          For channel 5, if BT_cld-BT_clr < -0.5 OR O-BT_clr < -0.5, reject channels 4-5.

    ! Arguments
    integer,     intent(in)                :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)                :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)                :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)                :: KNT                            ! nombre de tovs
    character *9,intent(in)                :: STNID                          ! identificateur du satellite
    integer,     intent(in)                :: KTERMER(KNT)                   ! land sea qualifyer 
    real,        intent(in)                :: PTBO(KNO,KNT)                  ! radiance o
    real,        intent(in)                :: btClear2D(KNO,KNT)             ! clear-radiance o
    real,        intent(in)                :: PTBOMP(KNO,KNT)                ! radiance o-p 
    integer,     intent(inout)             :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)             :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    ! Locals
    integer :: channelval
    integer :: nDataIndex
    integer :: nChannelIndex
    integer :: testIndex
    integer :: INDXCAN 
    real :: BTcloudy
    real :: simulatedCloudEffect
    real :: observedCloudEffect
    logical :: surfTypeIsWater 
    logical :: rejectLowPeakingChannels

    integer, dimension(2), parameter :: lowPeakingChannelsList = (/ 31, 32 /)

    testIndex = 16
    if ( .not. tvs_mwAllskyAssim ) return 

    loopObs: do nDataIndex = 1, KNT
      surfTypeIsWater = ( ktermer(nDataIndex) ==  1 )
      if ( .not. surfTypeIsWater ) cycle loopObs

      rejectLowPeakingChannels = .false.
      loopChannel2: do nChannelIndex = 1, KNO
        channelval = KCANO(nChannelIndex,nDataIndex)
        if ( channelval /= 32 ) cycle loopChannel2

        BTcloudy = PTBO(nChannelIndex,nDataIndex) - PTBOMP(nChannelIndex,nDataIndex)
        simulatedCloudEffect = BTcloudy - btClear2D(nChannelIndex,nDataIndex)
        observedCloudEffect = PTBO(nChannelIndex,nDataIndex) - btClear2D(nChannelIndex,nDataIndex)
        if ( simulatedCloudEffect < -0.5 .or. observedCloudEffect < -0.5 ) then
          rejectLowPeakingChannels = .true.
        end if

        exit loopChannel2
      end do loopChannel2

      ! reject channel 4-5
      if ( rejectLowPeakingChannels ) then
        do nChannelIndex = 1, KNO
          channelval = KCANO(nChannelIndex,nDataIndex)
          INDXCAN = ISRCHEQI(lowPeakingChannelsList,size(lowPeakingChannelsList),channelval)
          if ( INDXCAN /= 0 )  then
            ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**16)
            rejectionCodArray(testIndex,channelval,KNOSAT) = &
                rejectionCodArray(testIndex,channelval,KNOSAT) + 1 
          end if

          if ( mwbg_debug ) then
            write(*,*) STNID(2:9),' extreme scattering check reject: ', &
                    ' obs location index = ', nChannelIndex, &
                    ' channel = 1-5'
          endif
        end do
      end if

    end do loopObs

  end subroutine amsuaTest16ExcludeExtremeScattering

  !--------------------------------------------------------------------------
  ! copy1Dimto2DimRealArray
  !--------------------------------------------------------------------------
  subroutine copy1Dimto2DimRealArray(oneDimArray, firstDim, secondDim, twoDimArray)
    !:Purpose: copy 1 dim Real Array into 2D real array given dim1 and dim2 
    implicit none
    ! Arguments
    integer, intent(in)                 :: firstDim
    integer, intent(in)                 :: secondDim
    real,    intent(in)                 :: oneDimArray(firstDim*secondDim)
    real,    intent(inout)              :: twoDimArray(firstDim, secondDim)
    
    !locals
    integer                             :: firstDimIndex 
    integer                             :: secondDimIndex 
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
  ! copy1Dimto2DimIntegerArray
  !--------------------------------------------------------------------------
  subroutine copy1Dimto2DimIntegerArray(oneDimArray, firstDim, secondDim, twoDimArray)
    !:Purpose: copy 1 dim Integer Array into 2D Integer array given dim1 and dim2 
    implicit none
    ! Arguments
    integer, intent(in)                 :: firstDim
    integer, intent(in)                 :: secondDim
    integer, intent(in)                 :: oneDimArray(firstDim*secondDim)
    integer, intent(inout)              :: twoDimArray(firstDim, secondDim)
    
    !locals
    integer                             :: firstDimIndex 
    integer                             :: secondDimIndex 
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
  ! mwbg_tovCheckAmsua
  !--------------------------------------------------------------------------
  subroutine mwbg_tovCheckAmsua(TOVERRST,  clwThreshArr, sigmaObsErr, useStateDepSigmaObs, &
                                IUTILST, KTERMER, ICANO, ZO, btClear, ZCOR, &
                                ZOMP, ICHECK, KNO, KNT, KNOSAT, ISCNPOS, MGINTRP, MTINTRP, GLINTRP, ITERRAIN, SATZEN, &
                                globMarq, IMARQ, ident, clwObs, clwFG, scatw, STNID, RESETQC, ZLAT)

  
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
    integer, intent(in)                    :: IUTILST(:,:)             ! channel Selection using array IUTILST(chan,sat)
    !                                                                    IUTILST = 0 (blacklisted)
    !                                                                    1 (assmilate)
    !                                                                    2 (assimilate over open water only)
    real(8), intent(in)                    :: TOVERRST(:,:)            ! l'erreur totale des TOVS
    real(8), intent(in)                    :: clwThreshArr(:,:,:)      ! 
    real(8), intent(in)                    :: sigmaObsErr(:,:,:)       ! 
    logical, intent(in)                    :: useStateDepSigmaObs(:,:) ! if using state dependent obs error
    integer, allocatable, intent(inout)    :: globMarq(:)              ! Marqueurs globaux  
    integer, intent(in)                    :: KTERMER(:)               ! indicateur terre/mer
    integer, intent(in)                    :: ISCNPOS(:)               ! position sur le "scan"
    integer, intent(in)                    :: ICANO(:)                 ! canaux des observations
    integer, intent(inout)                 :: ITERRAIN(:)              ! indicateur du type de terrain
    integer, intent(in)                    :: KNO                      ! nombre de canaux des observations 
    integer, intent(in)                    :: KNT                      ! nombre de tovs
    integer, intent(in)                    :: KNOSAT                   ! numero de satellite (i.e. indice)
    integer, intent(inout)                 :: IMARQ(:)                 ! marqueurs des radiances
    real, intent(in)                       :: ZO(:)                    ! radiances
    real, intent(in)                       :: btClear(:)               ! clear-sky radiances
    real, intent(in)                       :: ZCOR(:)                  ! correction aux radiances
    real, intent(in)                       :: ZOMP(:)                  ! residus (o-p)
    real, intent(in)                       :: MGINTRP(:)               ! masque terre/mer du modele
    real, intent(in)                       :: MTINTRP(:)               ! topographie du modele
    real, intent(in)                       :: GLINTRP(:)               ! etendue de glace du modele
    real, intent(in)                       :: SATZEN(:)                ! angle zenith du satellite (deg.)
    real, intent(in)                       :: ZLAT(:)                  ! latitude
    character *9, intent(in)               :: STNID                    ! identificateur du satellite
    logical, intent(in)                    :: RESETQC                  ! reset du controle de qualite?
    integer, allocatable, intent(out)      :: ICHECK(:,:)              ! indicateur controle de qualite tovs par canal 
    !                                                                    =0, ok,
    !                                                                     >0, rejet,
    real, allocatable, intent(out)         :: clwObs(:)                ! retrieved cloud liquid water from observation 
    real, allocatable, intent(out)         :: clwFG(:)                 ! retrieved cloud liquid water from background 
    real, allocatable, intent(out)         :: scatw(:)                 ! scattering index over water

    integer, allocatable, intent(out)      :: ident(:)                 !ATMS Information flag (ident) values (new BURP element 025174 in header)
    !                                                                   FOR AMSUA just fill with zeros
    !locals
    integer, parameter                     :: mwbg_maxScanAngleHIRS= 56 
    integer, parameter                     :: maxScanAngleAMSU= 30 
    integer, parameter                     :: MXCLWREJ  =  6 
    integer, parameter                     :: MXSFCREJ  =  6 
    integer, parameter                     :: MXSFCREJ2 =  4 
    integer, parameter                     :: MXSCATREJ =  7 
    integer, parameter                     :: MXCANPRED =  9 
    integer, parameter                     :: JPMXSFC = 2
    real, parameter                        :: cloudyClwThreshold = 0.3
    real, parameter                        :: ZANGL = 117.6/maxScanAngleAMSU
    
    integer                                :: KMARQ   (KNO,KNT)
    integer                                :: KCANO   (KNO,KNT)
    real                                   :: PTBO    (KNO,KNT)
    real                                   :: btClear2D(KNO,KNT)
    real                                   :: PTBCOR  (KNO,KNT)
    real                                   :: PTBOMP  (KNO,KNT)
    integer, allocatable                   :: KCHKPRF(:)
    integer                                :: JI
    integer                                :: JJ
    integer                                :: INDX
    integer                                :: ICLWREJ (MXCLWREJ)
    integer                                :: ISFCREJ (MXSFCREJ)
    integer                                :: ISFCREJ2(MXSFCREJ2)
    integer                                :: ISCATREJ(MXSCATREJ)
    real                                   :: EPSILON
    real                                   :: MISGRODY
    real, allocatable                      :: GROSSMIN(:)
    real, allocatable                      :: GROSSMAX(:) 
    real, allocatable                      :: ROGUEFAC(:)
    real                                   :: tb23 (KNT)
    real                                   :: tb31 (KNT)
    real                                   :: tb50 (KNT)
    real                                   :: tb53 (KNT)
    real                                   :: tb89 (KNT)
    real                                   :: tb23FG (KNT)
    real                                   :: tb31FG (KNT)
    real                                   :: tb50FG (KNT)
    real                                   :: tb53FG (KNT)
    real                                   :: tb89FG (KNT)
    real                                   :: ice  (KNT)
    real                                   :: tpw  (KNT)
    real                                   :: scatl(KNT)
    integer                                :: err (KNT)
    integer                                :: rain(KNT)
    integer                                :: snow(KNT)
    integer                                :: channelForTopoFilter(2)
    real                                   :: altitudeForTopoFilter(2)
    logical, save                          :: LLFIRST = .true.

    EPSILON = 0.01
    MISGRODY = -99.

    call utl_reAllocate(ROGUEFAC, KNO+tvs_channelOffset(KNOSAT))
    ROGUEFAC(:) =(/ 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
                     4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
                     4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 2.0, 2.0, 2.0, &
                     3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
                     4.0, 2.0/)
    ICLWREJ(:) = (/ 28, 29, 30, 31, 32, 42 /)
    ISFCREJ(:) = (/ 28, 29, 30, 31, 32, 42 /)
    ISCATREJ(:) = (/ 28, 29, 30, 31, 32, 33, 42 /)
    ISFCREJ2(:) = (/ 28, 29, 30, 42 /)

    call utl_reAllocate(GROSSMIN, KNO+tvs_channelOffset(KNOSAT))
    GROSSMIN(:) = (/ 200., 190., 190., 180., 180., 180., 170., &
                    170., 180., 170., 170., 170., 180., 180., &
                    180., 180., 170., 180., 180., 000., 120., &
                    190., 180., 180., 180., 190., 200., 120., &
                    120., 160., 190., 190., 200., 190., 180., &
                    180., 180., 180., 190., 190., 200., 130./)

    call utl_reAllocate(GROSSMAX, KNO+tvs_channelOffset(KNOSAT))
    GROSSMAX(:) = (/ 270., 250., 250., 250., 260., 280., 290., &
                    320., 300., 320., 300., 280., 320., 300., &
                    290., 280., 330., 350., 350., 000., 310., &
                    300., 250., 250., 270., 280., 290., 310., &
                    310., 310., 300., 300., 260., 250., 250., &
                    250., 260., 260., 270., 280., 290., 330./)  
    channelForTopoFilter(:) = (/ 33, 34 /)
    altitudeForTopoFilter(:) = (/ 250., 2000./)
    ! Allocation
    call utl_reAllocate(clwObs,  KNT)
    call utl_reAllocate(clwFG,   KNT)
    call utl_reAllocate(scatw,   KNT)

    call utl_reAllocate(kchkprf, KNT)
    call utl_reAllocate(ident, KNT)
    call utl_reAllocate(icheck, KNO, KNT)

    ! copy the original input 1D array to 2D array. The 2D arrays are used in this s/r.
    call copy1Dimto2DimIntegerArray(ICANO, KNO, KNT, KCANO)
    call copy1Dimto2DimIntegerArray(IMARQ, KNO, KNT, KMARQ)
    call copy1Dimto2DimRealArray(ZCOR, KNO, KNT, PTBCOR)
    call copy1Dimto2DimRealArray(ZO, KNO, KNT, PTBO)
    call copy1Dimto2DimRealArray(btClear, KNO, KNT, btClear2D)
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

    !     Grody parameters are   extract required channels:
    call extractParamForGrodyRun (KCANO, ptbo, ptbomp, ptbcor, KNT, KNO, &
                                     tb23,   tb31,   tb50,   tb53,   tb89, &
                                     tb23FG, tb31FG, tb50FG, tb53FG, tb89FG)
    
    !  Run Grody AMSU-A algorithms.

    call grody (err, knt, tb23, tb31, tb50, tb53, tb89, tb23FG, tb31FG, &
                satzen, zlat, ktermer, ice, tpw, clwObs, clwFG, &
                rain, snow, scatl, scatw)   

    ! 10) test 10: RTTOV reject check (single)
    ! Rejected datum flag has bit #9 on.
    call amsuABTest10RttovRejectCheck (KCANO, KNOSAT, KNO, KNT, RESETQC, STNID, KMARQ, ICHECK)

    ! 1) test 1: Topography check (partial)
    ! Channel 6 is rejected for topography >  250m.
    ! Channel 7 is rejected for topography > 2000m.
    call amsuABTest1TopographyCheck (KCANO, KNOSAT, KNO, KNT, STNID, MTINTRP, &
                                        channelForTopoFilter, altitudeForTopoFilter, KMARQ, &
                                        ICHECK)
 
    ! 2) test 2: "Land/sea qualifier" code check (full)
    ! allowed values are: 0, land,
    !                       1, sea,
    !                       2, coast.
    call amsuABTest2LandSeaQualifierCheck (KCANO, KNOSAT, KNO, KNT, STNID, KTERMER, KMARQ, ICHECK)

    ! 3) test 3: "Terrain type" code check (full)
    !   allowed values are: -1, missing,
    !                        0, sea-ice,
    !                        1, snow on land.
    call amsuABTest3TerrainTypeCheck (KCANO, KNOSAT, KNO, KNT, STNID, ITERRAIN, KMARQ, ICHECK)
 
    ! 4) test 4: Field of view number check (full)
    !
    ! Field of view acceptable range is [1,maxScanAngleAMSU]  for AMSU footprints.
    call amsuABTest4FieldOfViewCheck (KCANO, KNOSAT, KNO, KNT, STNID, ISCNPOS, maxScanAngleAMSU, KMARQ, ICHECK)

    ! 5) test 5: Satellite zenith angle check (full)
    ! Satellite zenith angle acceptable range is [0.,60.].
    call amsuABTest5ZenithAngleCheck (KCANO, KNOSAT, KNO, KNT, STNID, SATZEN, KMARQ, ICHECK)
    ! 6) test 6: "Sat. zenith angle"/"field of view" consistency check.  (full)
    ! Acceptable difference between "Satellite zenith angle"  and
    ! "approximate angle computed from field of view number" is 1.8 degrees.
    call amsuABTest6ZenAngleAndFovConsistencyCheck (KCANO, KNOSAT, KNO, KNT, STNID, SATZEN,  ZANGL, &
                                                  ISCNPOS, maxScanAngleAMSU, KMARQ, ICHECK) 
    ! 7) test 7: "Land/sea qual."/"model land/sea" consistency check.    (full)
    ! Acceptable conditions are:
    !       a) both over ocean (ktermer=1; mg<0.01), new threshold 0.20, jh dec 2000,
    !       b) both over land  (ktermer=0; mg>0.80), new threshold 0.50, jh dec 2000.
    ! Other conditions are unacceptable.
    call amsuABTest7landSeaQualifyerAndModelLandSeaConsistencyCheck (KCANO, KNOSAT, KNO, KNT, STNID, MGINTRP, KTERMER, &
                                                               KMARQ, ICHECK)

    ! 8) test 8: "Terrain type"/"Land/sea qual."/"model ice" consistency check. (full)
    ! Unacceptable conditions are:
    !        a) terrain is sea-ice and model has no ice(iterrain=0; gl<0.01).
    !        b) terrain is sea-ice and land/sea qualifier is land (iterrain=0; ktermer=0).
    !        c) terrain is snow on land and land/sea qualifier is sea (iterrain=1; ktermer=1).
    !        d) terrain is missing, land/sea qualifier is sea and model has ice(iterrain=-1; ktermer=1; gl>0.01). (enleve jh, jan 2001)
    ! NOT doNE ANYMORE 
    
    ! 9) test 9: Uncorrected Tb check (single)
    ! Uncorrected datum (flag bit #6 off). In this case switch bit 11 ON.
    call amsuABTest9UncorrectedTbCheck (KCANO, KNOSAT, KNO, KNT, STNID, RESETQC, KMARQ, ICHECK) 
    ! 11) test 11: Radiance observation "Gross" check (single) 
    !  Change this test from full to single. jh nov 2000.
    call amsuABTest11RadianceGrossValueCheck (KCANO, KNOSAT, KNO, KNT, STNID, PTBO, GROSSMIN, &
                                      GROSSMAX, KMARQ, ICHECK)
    ! 12) test 12: Grody cloud liquid water check (partial)
    ! For Cloud Liquid Water > clwQcThreshold, reject AMSUA-A channels 1-5 and 15.
    call amsuaTest12GrodyClwCheck (KCANO, KNOSAT, KNO, KNT, STNID, clwObs, &
                                clwFG, useStateDepSigmaObs, ktermer, MISGRODY, MXCLWREJ, &
                                ICLWREJ, KMARQ, ICHECK)
    ! 13) test 13: Grody scattering index check (partial)
    ! For Scattering Index > 9, reject AMSUA-A channels 1-6 and 15.
    call amsuaTest13GrodyScatteringIndexCheck (KCANO, KNOSAT, KNO, KNT, STNID, scatw, KTERMER, ITERRAIN, &
                                              MISGRODY, MXSCATREJ, ISCATREJ, KMARQ, ICHECK)

    ! 14) test 14: "Rogue check" for (O-P) Tb residuals out of range. (single/full)
    ! Les observations, dont le residu (O-P) depasse par un facteur (roguefac) l'erreur totale des TOVS.
    ! N.B.: a reject by any of the 3 surface channels produces the rejection of AMSUA-A channels 1-5 and 15. 
    call amsuaTest14RogueCheck (KCANO, KNOSAT, KNO, KNT, STNID, ROGUEFAC, TOVERRST, clwThreshArr, &
                                    useStateDepSigmaObs, sigmaObsErr, ktermer, PTBOMP, clwObs, clwFG, &
                                    MISGRODY, MXSFCREJ, ISFCREJ, KMARQ, ICHECK)

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
    call amsuABTest15ChannelSelectionWithIutilst (KCANO, KNOSAT, KNO, KNT, STNID, KTERMER, ITERRAIN, GLINTRP, IUTILST, &
                                              MXSFCREJ2, ISFCREJ2, KMARQ, ICHECK)

    ! 16) test 16: exclude radiances affected by extreme scattering in deep convective region in all-sky mode.
    call amsuaTest16ExcludeExtremeScattering(KCANO, KNOSAT, KNO, KNT, STNID, KTERMER, &
                                        PTBO, btClear2D, PTBOMP,  KMARQ, ICHECK) 

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

    ! reset global marker flag (55200) and mark it if observtions are rejected
    call resetQcCases(RESETQC, KCHKPRF, globMarq)

    !###############################################################################
    ! FINAL STEP: set terrain type to sea ice given certain conditions
    !###############################################################################
    call setTerrainTypeToSeaIce(GLINTRP, KTERMER, ITERRAIN)

  end subroutine mwbg_tovCheckAmsua

  !--------------------------------------------------------------------------
  ! mwbg_tovCheckAmsub
  !--------------------------------------------------------------------------
  subroutine mwbg_tovCheckAmsub(TOVERRST, siThreshArr, sigmaObsErr, useStateDepSigmaObs, &
                                IUTILST,  KTERMER, ICANO, ZO, ZCOR, &
                                ZOMP, ICHECK, KNO, KNT, KNOSAT, ISCNPOS, MGINTRP, MTINTRP, GLINTRP, ITERRAIN, SATZEN, &
                                globMarq, IMARQ, ident, clwOBS, clwFG, scatwObs, scatwFG, STNID, RESETQC)

  
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
    !                  - 8) inconsistent terrain type and model ice
    !                  - 9) uncorrected radiance,
    !                  - 10) rejected by RTTOV,
    !                  - 11) radiance gross check failure,
    !                  - 12) drynes index reject
    !                  - 13) scattering index reject,
    !                  - 14) radiance residual rogue check failure,
    !                  - 15) channel reject (channel selection).
    !                  - **) set terrain type to sea ice given certain conditions
    implicit none 
    !Arguments:
    integer, intent(in)                    :: IUTILST(:,:)       !channel Selection using array IUTILST(chan,sat)
    !                                                               IUTILST = 0 (blacklisted)
    !                                                               1 (assmilate)
    !                                                               2 (assimilate over open water only)
    real(8), intent(in)                    :: TOVERRST(:,:)      ! l'erreur totale des TOVS
    real(8), intent(in)                    :: siThreshArr(:,:,:)  ! SI thresholds for state-dep obs err
    real(8), intent(in)                    :: sigmaObsErr(:,:,:)  ! sigmaObs limits for state-dep obs err    
    logical, intent(in)                    :: useStateDepSigmaObs(:,:) ! if using state dependent obs error
    integer, allocatable, intent(inout)    :: globMarq(:)        !Marqueurs globaux  
    integer, intent(in)                    :: KTERMER(:)         ! indicateur terre/mer
    integer, intent(in)                    :: ISCNPOS(:)         ! position sur le "scan"
    integer, intent(in)                    :: ICANO(:)           ! canaux des observations
    integer, intent(inout)                 :: ITERRAIN(:)        ! indicateur du type de terrain
    integer, intent(in)                    :: KNO                ! nombre de canaux des observations 
    integer, intent(in)                    :: KNT                ! nombre de tovs
    integer, intent(in)                    :: KNOSAT             ! numero de satellite (i.e. indice)
    integer, intent(inout)                 :: IMARQ(:)           ! marqueurs des radiances
    real, intent(in)                       :: ZO(:)              ! radiances
    real, intent(in)                       :: ZCOR(:)            ! correction aux radiances
    real, intent(in)                       :: ZOMP(:)            ! residus (o-p)
    real, intent(in)                       :: MGINTRP(:)         ! masque terre/mer du modele
    real, intent(in)                       :: MTINTRP(:)         ! topographie du modele
    real, intent(in)                       :: GLINTRP(:)         ! etendue de glace du modele
    real, intent(in)                       :: SATZEN(:)          ! angle zenith du satellite (deg.)
    character *9, intent(in)               :: STNID              ! identificateur du satellite
    logical, intent(in)                    :: RESETQC            ! reset du controle de qualite?
    integer, allocatable, intent(out)      :: ICHECK(:,:)        ! indicateur controle de qualite tovs par canal 
    !                                                              =0, ok,
    !                                                              >0, rejet,
    real, allocatable, intent(out)         :: clwObs(:)          ! retrieved cloud liquid water from observation 
    real, allocatable, intent(out)         :: clwFG(:)           ! retrieved cloud liquid water from background 
    real, allocatable, intent(out)         :: scatwObs(:)        ! scattering index over water from observation
    real, allocatable, intent(out)         :: scatwFG(:)         ! scattering index over water from background

    integer, allocatable, intent(out)       :: ident(:)          !ATMS Information flag (ident) values (new BURP element 025174 in header)
    !                                                               FOR AMSUA just fill with zeros

    !locals
    integer, parameter                     :: mwbg_maxScanAngleHIRS= 56 
    integer, parameter                     :: maxScanAngleAMSU= 90 
    integer, parameter                     :: MXCLWREJ  =  6 
    integer, parameter                     :: MXSFCREJ  =  2 
    integer, parameter                     :: MXSFCREJ2 =  1 
    integer, parameter                     :: MXCANPRED =  9 
    integer, parameter                     :: MXCH2OMPREJ= 4
    integer, parameter                     :: JPMXSFC = 2
    real, parameter                        :: cloudyClwThreshold = 0.3
    real, parameter                        :: ZANGL =  117.6/maxScanAngleAMSU
    
    integer                                :: KMARQ   (KNO,KNT)
    integer                                :: KCANO   (KNO,KNT)
    real                                   :: PTBO    (KNO,KNT)
    real                                   :: PTBCOR  (KNO,KNT)
    real                                   :: PTBOMP  (KNO,KNT)
    integer, allocatable                   :: KCHKPRF(:)          
    integer                                :: JI
    integer                                :: JJ
    integer                                :: INDX
    integer                                :: ICLWREJ (MXCLWREJ)
    integer                                :: ISFCREJ (MXSFCREJ)
    integer                                :: ICH2OMPREJ(MXCH2OMPREJ)
    integer                                :: ISFCREJ2(MXSFCREJ2)
    real                                   :: EPSILON
    real                                   :: MISGRODY
    real, allocatable                      :: GROSSMIN(:)
    real, allocatable                      :: GROSSMAX(:) 
    real, allocatable                      :: ROGUEFAC(:)
    real                                   :: tb89 (KNT)
    real                                   :: tb150 (KNT)
    real                                   :: tb1831 (KNT)
    real                                   :: tb1832 (KNT)
    real                                   :: tb1833 (KNT)
    real                                   :: tb89FG (KNT)
    real                                   :: tb150FG(KNT)
    real                                   :: scatl(KNT)
    integer                                :: err (KNT)
    integer                                :: channelForTopoFilter(3)
    real                                   :: altitudeForTopoFilter(3)
    logical, save                          :: LLFIRST = .true.

    EPSILON = 0.01
    MISGRODY = -99.

    call utl_reAllocate(ROGUEFAC, KNO+tvs_channelOffset(KNOSAT))
    ROGUEFAC(:) =(/ 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
                    4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
                    4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 2.0, 2.0, 2.0, &
                    3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
                    4.0, 2.0, 2.0, 2.0, 4.0, 4.0, 4.0/)

    ICLWREJ(:) = (/ 28, 29, 30, 31, 32, 42 /)
    ISFCREJ(:) = (/ 43, 44 /)
    ISFCREJ2(:) = (/ 43 /)
    ICH2OMPREJ(:) = (/ 44, 45, 46, 47 /)
    call utl_reAllocate(GROSSMIN, KNO+tvs_channelOffset(KNOSAT))
    GROSSMIN(:) = (/ 200., 190., 190., 180., 180., 180., 170., &
                    170., 180., 170., 170., 170., 180., 180., &
                    180., 180., 170., 180., 180., 000., 120., &
                    190., 180., 180., 180., 190., 200., 120., &
                    120., 160., 190., 190., 200., 190., 180., &
                    180., 180., 180., 190., 190., 200., 130., &
                    130., 130., 130., 130., 130./)
    call utl_reAllocate(GROSSMAX, KNO+tvs_channelOffset(KNOSAT))
    GROSSMAX(:) = (/ 270., 250., 250., 250., 260., 280., 290., &
                    320., 300., 320., 300., 280., 320., 300., &
                    290., 280., 330., 350., 350., 000., 310., &
                    300., 250., 250., 270., 280., 290., 310., &
                    310., 310., 300., 300., 260., 250., 250., &
                    250., 260., 260., 270., 280., 290., 330., &
                    330., 330., 330., 330., 330./)  

    channelForTopoFilter(:) = (/ 45, 46, 47 /)
    altitudeForTopoFilter(:) = (/ 2500., 2000., 1000./)

    ! Allocation
    call utl_reAllocate(scatwObs, KNT)
    call utl_reAllocate(scatwFG,  KNT)
    call utl_reAllocate(clwObs,   KNT)
    call utl_reAllocate(clwFG,    KNT)

    call utl_reAllocate(kchkprf, KNT)
    call utl_reAllocate(ident, KNT)
    call utl_reAllocate(icheck, KNO, KNT)

    ! copy the original input 1D array to 2D array. The 2D arrays are used in this s/r.
    call copy1Dimto2DimIntegerArray(ICANO, KNO, KNT, KCANO)
    call copy1Dimto2DimIntegerArray(IMARQ, KNO, KNT, KMARQ)
    call copy1Dimto2DimRealArray(ZCOR, KNO, KNT, PTBCOR)
    call copy1Dimto2DimRealArray(ZO, KNO, KNT, PTBO)
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

    !     Bennartz parameters are   extract required channels:
    call extractParamForBennartzRun (KCANO, ptbo, ptbomp, ptbcor, KNT, KNO, &
                                     tb89, tb150, tb1831, tb1832, tb1833, &
                                     tb89FG, tb150FG)
    
    !  Run Bennartz AMSU-B algorithms.
    call bennartz (err, knt, tb89, tb150, tb89FG, tb150FG, satzen, ktermer, &
                   scatl, scatwObs, scatwFG, clwObs, clwFG)   

    ! 10) test 10: RTTOV reject check (single)
    ! Rejected datum flag has bit #9 on.
    call amsuABTest10RttovRejectCheck (KCANO, KNOSAT, KNO, KNT, RESETQC, STNID, KMARQ, ICHECK)

    ! 1) test 1: Topography check (partial)
    ! Channel 3- 45 is rejected for topography >  2500m.
    ! Channel 4 - 46 is rejected for topography > 2000m.
    ! Channel 5 - 47 is rejected for topography > 1000m.
    call amsuABTest1TopographyCheck (KCANO, KNOSAT, KNO, KNT, STNID, MTINTRP, &
                                        channelForTopoFilter, altitudeForTopoFilter, KMARQ, &
                                        ICHECK)
 
    ! 2) test 2: "Land/sea qualifier" code check (full)
    ! allowed values are: 0, land,
    !                       1, sea,
    !                       2, coast.
    call amsuABTest2LandSeaQualifierCheck (KCANO, KNOSAT, KNO, KNT, STNID, KTERMER, KMARQ, ICHECK)

    ! 3) test 3: "Terrain type" code check (full)
    !   allowed values are: -1, missing,
    !                        0, sea-ice,
    !                        1, snow on land.
    call amsuABTest3TerrainTypeCheck (KCANO, KNOSAT, KNO, KNT, STNID, ITERRAIN, KMARQ, ICHECK)
 
    ! 4) test 4: Field of view number check (full)
    !
    ! Field of view acceptable range is [1,maxScanAngleAMSU]  for AMSU footprints.
    call amsuABTest4FieldOfViewCheck (KCANO, KNOSAT, KNO, KNT, STNID, ISCNPOS, maxScanAngleAMSU, KMARQ, ICHECK)

    ! 5) test 5: Satellite zenith angle check (full)
    ! Satellite zenith angle acceptable range is [0.,60.].
    call amsuABTest5ZenithAngleCheck (KCANO, KNOSAT, KNO, KNT, STNID, SATZEN, KMARQ, ICHECK)
    ! 6) test 6: "Sat. zenith angle"/"field of view" consistency check.  (full)
    ! Acceptable difference between "Satellite zenith angle"  and
    ! "approximate angle computed from field of view number" is 1.8 degrees.
    call amsuABTest6ZenAngleAndFovConsistencyCheck (KCANO, KNOSAT, KNO, KNT, STNID, SATZEN,  ZANGL, &
                                                  ISCNPOS, maxScanAngleAMSU, KMARQ, ICHECK) 
    ! 7) test 7: "Land/sea qual."/"model land/sea" consistency check.    (full)
    ! Acceptable conditions are:
    !       a) both over ocean (ktermer=1; mg<0.01), new threshold 0.20, jh dec 2000,
    !       b) both over land  (ktermer=0; mg>0.80), new threshold 0.50, jh dec 2000.
    ! Other conditions are unacceptable.
    call amsuABTest7landSeaQualifyerAndModelLandSeaConsistencyCheck (KCANO, KNOSAT, KNO, KNT, STNID, MGINTRP, KTERMER, &
                                                               KMARQ, ICHECK)

    ! 8) test 8: "Terrain type"/"Land/sea qual."/"model ice" consistency check. (full)
    ! Unacceptable conditions are:
    !        a) terrain is sea-ice and model has no ice(iterrain=0; gl<0.01).
    !        b) terrain is sea-ice and land/sea qualifier is land (iterrain=0; ktermer=0).
    !        c) terrain is snow on land and land/sea qualifier is sea (iterrain=1; ktermer=1).
    !        d) terrain is missing, land/sea qualifier is sea and model has ice(iterrain=-1; ktermer=1; gl>0.01). (enleve jh, jan 2001)
    ! NOT doNE ANYMORE 
    
    ! 9) test 9: Uncorrected Tb check (single) SKIP FOR NOW
    ! Uncorrected datum (flag bit #6 off). In this case switch bit 11 ON.
    ! call amsuABTest9UncorrectedTbCheck (KCANO, KNOSAT, KNO, KNT, STNID, RESETQC, KMARQ, ICHECK) 
    ! 11) test 11: Radiance observation "Gross" check (single) 
    !  Change this test from full to single. jh nov 2000.
    call amsuABTest11RadianceGrossValueCheck (KCANO, KNOSAT, KNO, KNT, STNID, PTBO, GROSSMIN, &
                                      GROSSMAX, KMARQ, ICHECK)
    ! 12) test 12:  Dryness index check 
    !The difference between channels AMSUB-3 and AMSUB-5 is used as an indicator
    !of "dryness" of the atmosphere. In extreme dry conditions, channels AMSUB-3 4 and 5
    ! are sensitive to the surface.
    ! Therefore, various thresholds are used to reject channels AMSUB-3 4 and 5
    !  over land and ice
    call amsubTest12DrynessIndexCheck (KCANO, KNOSAT, KNO, KNT, STNID, tb1831, tb1833, &
                                       ktermer, glintrp, KMARQ, ICHECK)
    ! 13) test 13: Bennartz scattering index check (full)

    call amsubTest13BennartzScatteringIndexCheck(KCANO, KNOSAT, KNT, KNO, STNID, scatwObs, scatwFG, scatl, &
                                                 useStateDepSigmaObs, KTERMER, GLINTRP, KMARQ, ICHECK)

    ! 14) test 14: "Rogue check" for (O-P) Tb residuals out of range. (single/full)
    ! Les observations, dont le residu (O-P) depasse par un facteur (roguefac) l'erreur totale des TOVS.
    ! N.B.: a reject by any of the 3 surface channels produces the rejection of AMSUA-A channels 1-5 and 15. 

    call amsubTest14RogueCheck(KCANO, KNOSAT, KNO, KNT, STNID, ROGUEFAC, TOVERRST, &
                               siThreshArr, sigmaObsErr, useStateDepSigmaObs, &
                               iterrain, ktermer, PTBOMP, ICH2OMPREJ, MXCH2OMPREJ, KMARQ, ICHECK)

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

    call amsuABTest15ChannelSelectionWithIutilst (KCANO, KNOSAT, KNO, KNT, STNID, KTERMER, ITERRAIN, GLINTRP, IUTILST, &
                                                 MXSFCREJ2, ISFCREJ2, KMARQ, ICHECK)
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

    ! reset global marker flag (55200) and mark it if observtions are rejected
    call resetQcCases(RESETQC, KCHKPRF, globMarq)

    !###############################################################################
    ! FINAL STEP: set terrain type to sea ice given certain conditions
    !###############################################################################
    call setTerrainTypeToSeaIce(GLINTRP, KTERMER, ITERRAIN)

  end subroutine mwbg_tovCheckAmsub

  !--------------------------------------------------------------------------
  ! mwbg_qcStats
  !--------------------------------------------------------------------------
  subroutine mwbg_qcStats(instName, ICHECK, ICAN, KNOSAT, &
                              KNO, KNT, satelliteId, LDprint)

    !:Purpose:          Cumuler ou imprimer des statistiques decriptives
    !                   des rejets tovs.
    implicit none 
    !Arguments:
    character(*), intent(in)               :: instName                           ! Instrument Name
    integer, intent(in)                    :: ICHECK(:,:)                        ! indicateur controle de qualite tovs par canal 
    !                                                                              =0, ok,
    !                                                                              >0, rejet,
    integer, intent(in)                    :: ICAN(KNO*KNT)                      ! canaux des observations
    integer, intent(in)                    :: KNOSAT                             ! numero d'identificateur du satellite
    integer, intent(in)                    :: KNO                                ! nombre de canaux des observations 
    integer, intent(in)                    :: KNT                                ! nombre de tovs
    character(len=15), intent(in)          :: satelliteId(:)                     ! identificateur du satellite
    logical, intent(in)                    :: LDprint                            ! mode: imprimer ou cumuler?
    !Locals
    integer                                :: numSats
    integer                                :: JI
    integer                                :: JJ
    integer                                :: JK
    integer                                :: INTOTOBS
    integer                                :: INTOTACC
    integer, allocatable, save             :: INTOT(:)!INTOT(mwbg_maxNumSat)
    integer, allocatable, save             :: INTOTRJF(:)!INTOTRJF(mwbg_maxNumSat)
    integer, allocatable, save             :: INTOTRJP(:)!INTOTRJP(mwbg_maxNumSat)
    integer                                :: KCANO(KNO,KNT)                      ! canaux des observations

    logical, save                          :: LLFIRST = .true.
    logical                                :: FULLREJCT
    logical                                :: FULLACCPT

    ! Initialize
    if ( LLFIRST ) then
      call utl_reallocate(INTOT, mwbg_maxNumSat)
      call utl_reallocate(INTOTRJF, mwbg_maxNumSat)
      call utl_reallocate(INTOTRJP, mwbg_maxNumSat)
      do JJ = 1, mwbg_maxNumSat
        INTOTRJF(JJ) = 0
        INTOTRJP(JJ) = 0
        INTOT(JJ)  = 0
      end do
      LLFIRST = .false.
    end if

    if (.not. LDprint ) then

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
            if ( KCANO(JI,JJ) /= 20 ) then
              if ( ICHECK(JI,JJ) /= 0 ) then
                FULLACCPT = .false.
              else
                FULLREJCT = .false.
              end if
            end if
          end do
          if ( FULLREJCT ) then
            INTOTRJF(KNOSAT) = INTOTRJF(KNOSAT) + 1
          end if
          if ( .not. FULLREJCT .and. .not.FULLACCPT ) then
            INTOTRJP(KNOSAT) = INTOTRJP(KNOSAT) + 1
          end if
        else if  (instName == "ATMS") then 
          do JI = 1, KNO
            if ( ICHECK(JI,JJ) /= 0 ) then
              FULLACCPT = .false.
            else
              FULLREJCT = .false.
            end if
          end do
          if ( FULLREJCT ) then
            INTOTRJF(KNOSAT) = INTOTRJF(KNOSAT) + 1
          end if
          if ( .not. FULLREJCT .and. .not.FULLACCPT ) then
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

        if (instName == "AMSUA" .or. instName == "AMSUB") then         
          write(*,'(//,1x,114("-"))')
           write(*,'(t10,"|",t47,"REJECTION CATEGORIES")')
          write(*,'(" CHANNEL",t10,"|",105("-"))')
          write(*,'(t10,"|",16i7)') (JI,JI=1,mwbg_maxNumTest)
          write(*,'(1x,"--------|",105("-"))')
          do JJ = 1, KNO 
             write(*,'(3X,I2,t10,"|",16I7)') JJ,(rejectionCodArray(JI,JJ,JK), &
                                      JI=1,mwbg_maxNumTest)
          end do
          write(*,'(1x,114("-"))')
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

        else if (instName == "ATMS" .or. instName == "MWHS2") then
          write(*,'(//,1x,59("-"))')
          write(*,'(t10,"|",t19,"1. REJECTION CATEGORIES")')
          write(*,'(" CHANNEL",t10,"|",50("-"))')
          write(*,'(t10,"|",5i7)') (JI,JI=1,mwbg_maxNumTest)
          write(*,'(1x,"--------|",50("-"))')
          do JJ = 1, KNO 
            write(*,'(3X,I2,t10,"|",5I7)') JJ,(rejectionCodArray(JI,JJ,JK), &
                                        JI=1,mwbg_maxNumTest)
          end do
          write(*,'(1x,59("-"))')
          write(*,'(//,1x,59("-"))')
          write(*,'(t10,"|",t19,"2. QC2 REJECT CATEGORIES")')
          write(*,'(" CHANNEL",t10,"|",50("-"))') 
          write(*,'(t10,"|",5i7)') (JI,JI=1,mwbg_maxNumTest)
          write(*,'(1x,"--------|",50("-"))')
          do JJ = 1, KNO
            write(*,'(3X,I2,t10,"|",5I7)') JJ,(rejectionCodArray2(JI,JJ,JK), &
                                        JI=1,mwbg_maxNumTest)          
          end do
          print *, ' '
          print *, ' '
          print *, ' -----------------------------------------------------'
          print *, ' Definition of rejection categories: '
          print *, ' -----------------------------------------------------'
          print *, '  1 - first bgckAtms/bgckMwhs2 program reject [bit 7]'
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

  end subroutine mwbg_qcStats

  !--------------------------------------------------------------------------
  ! resetQcC
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

    dataNum = size(globMarq)
    do dataIndex = 1, dataNum
      if (RESETQC) then
        globMarq(dataIndex) = 1024  
      end if
      if ( KCHKPRF(dataIndex) /= 0  ) then
        globMarq(dataIndex) = OR (globMarq(dataIndex),2**6)
      end if
    end do
    if (mwbg_debug) then
      write(*,*) ' KCHKPRF   = ', (KCHKPRF(dataIndex),dataIndex=1,dataNum)
      write(*,*) ' NEW FLAGS = ', (globMarq(dataIndex),dataIndex=1,dataNum)
    end if

  end  subroutine resetQcCases

  !--------------------------------------------------------------------------
  ! setTerrainTypeToSeaIce
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
  ! GRODY
  !--------------------------------------------------------------------------
  subroutine GRODY (ier, ni, tb23, tb31, tb50, tb53, tb89, tb23FG, tb31FG, &
                   pangl, plat, ilansea, ice, tpw, clwObs, clwFG, &
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
    !               - tb23FG  - input  -  23Ghz brightness temperature from background (K)
    !               - tb31FG  - input  -  31Ghz brightness temperature from background (K)
    !               - pangl   - input  -  satellite zenith angle (deg.)
    !               - plat    - input  -  lalitude (deg.)
    !               - ilansea - input  -  land/sea indicator (0=land;1=ocean)
    !               - ice     - output -  sea ice concentration (0-100%)
    !               - tpw     - output -  total precipitable water (0-70mm)
    !               - clwObs  - output -  retrieved cloud liquid water from observation (0-3mm)
    !               - clwFG   - output -  retrieved cloud liquid water from background (0-3mm)
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

    real zmisgLocal, siw, sil, df1, df2, df3, a, b, c, d, e23
    real ei, cosz, tt, scat, sc31, abslat, t23, t31, t50, t89
    real sc50, par, t53
    real dif285t23, dif285t31, epsilon
    real dif285t23FG, dif285t31FG

    real tb23  (:)
    real tb31  (:)
    real tb50  (:)
    real tb53  (:)
    real tb89  (:)
    real tb23FG(:)
    real tb31FG(:)
    real pangl (:)
    real plat  (:)
    real ice   (:)
    real tpw   (:)
    real clwObs(:)
    real clwFG (:)
    real scatl (:)
    real scatw (:)

    data zmisgLocal   / -99.     /
    data epsilon /   1.E-30 /

    logical skipLoopChan15Missing 

    ! 1) Initialise output parameters:
    do i = 1, ni
      ice  (i) = zmisgLocal
      tpw  (i) = zmisgLocal
      clwObs(i) = zmisgLocal
      clwFG(i) = zmisgLocal
      scatl(i) = zmisgLocal
      scatw(i) = zmisgLocal
      rain (i) = nint(zmisgLocal)
      snow (i) = nint(zmisgLocal)
    end do

    ! 2) Validate input parameters:
    do i = 1, ni
      if ( tb23(i)    < 120.  .or. &
           tb23(i)    > 350.  .or. &
           tb31(i)    < 120.  .or. &
           tb31(i)    > 350.  .or. &
           tb50(i)    < 120.  .or. &
           tb50(i)    > 350.  .or. &
           tb53(i)    < 120.  .or. &
           tb53(i)    > 350.  .or. &
           pangl(i)   < -90.  .or. &
           pangl(i)   > 90.   .or. &
           plat(i)    < -90.  .or. &
           plat(i)    > 90.   .or. &
           ilansea(i) < 0     .or. &
           ilansea(i) > 1        ) then
        ier(i) = 1
      else
        ier(i) = 0
      end if

    end do

    !3) Compute parameters:
    loopObsGrody: do i = 1, ni
      if ( ier(i) .eq. 0 ) then
        abslat = abs(plat(i))
        cosz   = cosd(pangl(i))
        t23 = tb23(i)
        t31 = tb31(i)
        t50 = tb50(i)
        t53 = tb53(i)
        dif285t23  =max(285.-t23,epsilon)
        dif285t23FG=max(285.-tb23FG(i),epsilon)
        dif285t31  =max(285.-t31,epsilon)
        dif285t31FG=max(285.-tb31FG(i),epsilon)

        skipLoopChan15Missing = .false.
        if ( tb89(i) < 120.0 .or. tb89(i) > 350.0 ) then 
          skipLoopChan15Missing = .true.
        end if

        if ( .not. skipLoopChan15Missing ) then
          t89 = tb89(i)

          ! scattering indices:
          siw = -113.2 + (2.41 - 0.0049*t23)*t23 &
                       + 0.454*t31 &
                       - t89 
          sil = t23 - t89 

          scatl (i) = sil
          scatw (i) = siw
        end if

        ! discriminate functions:
        df1 =  2.85 + 0.020*t23 - 0.028*t50 ! used to identify (also remove) sea ice
        df2 =  5.10 + 0.078*t23 - 0.096*t50 ! used to identify (also remove) warm deserts
        df3 = 10.20 + 0.036*t23 - 0.074*t50 ! used to identify (also remove) cold deserts

        if ( ilansea(i)== 1 ) then

          ! Ocean Parameters

          !3.1) Sea ice:
          if ( abslat < 50. ) then
            ice(i) = 0.0
          else
            if ( df1 < 0.45 ) then
               ice(i) = 0.0
            else
              a =  1.7340  - 0.6236*cosz
              b =  0.0070  + 0.0025*cosz
              c = -0.00106 
              d = -0.00909
              e23 = a + b*t31 + c*t23 + d*t50 ! theoretical 23Ghz sfc emissivity (0.3-1.)
              if ( (t23-t31) >= 5. ) then   ! fov contains multiyear or new ice/water
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
          if ( abslat > 50.  .and. &
              df1     >  0.2        ) then  
            tpw(i) = zmisgLocal
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

          !3.3) Cloud liquid water from obs (clwObs) and background state (clwFG):
          ! identify and remove sea ice
          if ( abslat > 50.  .and. &
              df1     >  0.0        ) then  
            clwObs(i) = zmisgLocal
            clwFG(i) = zmisgLocal
          else
            a =  8.240 - (2.622 - 1.846*cosz)*cosz
            b =  0.754
            c = -2.265
            clwObs(i) = a + b*log(dif285t23) & 
                      + c*log(dif285t31)
            clwObs(i) = clwObs(i)*cosz           ! theoretical cloud liquid water (0-3mm)
            clwObs(i) = clwObs(i) - 0.03         ! corrected   cloud liquid water 
            clwObs(i) = min(3.,max(0.,clwObs(i)))   ! jh       

            clwFG(i) = a + b*log(dif285t23FG) & 
                      + c*log(dif285t31FG)
            clwFG(i) = clwFG(i)*cosz           ! theoretical cloud liquid water (0-3mm)
            clwFG(i) = clwFG(i) - 0.03         ! corrected   cloud liquid water 
            clwFG(i) = min(3.,max(0.,clwFG(i)))   ! jh       
          endif

          if ( skipLoopChan15Missing ) cycle loopObsGrody

          !3.4) Ocean rain: 0=no rain; 1=rain.
          ! identify and remove sea ice
          if ( abslat > 50.  .and. &
              df1    >  0.0        ) then  
            rain(i) = nint(zmisgLocal)
          else                                   ! remove non-precipitating clouds
            if ( clwObs(i) > 0.3 .or. &
                siw        > 9.0      ) then 
              rain(i) = 1
            else
              rain(i) = 0
            end if
          end if

        else

          if ( skipLoopChan15Missing ) cycle loopObsGrody

          ! Land Parameters

          ! 3.5) Rain  over land: 0=no rain; 1=rain.
          tt = 168. + 0.49*t89
          if ( sil >= 3. ) then
            rain(i) = 1
          else 
            rain(i) = 0
          end if
          
          ! remove snow cover
          if ( t23 <= 261. .and. &
              t23 < tt         ) then
            rain(i) = 0
          end if

          ! remove warm deserts
          if ( t89 > 273.  .or. &
              df2  <   0.6      ) then
            rain(i) = 0
          end if

          ! 3.6) Snow cover and glacial ice: 0=no snow; 1=snow; 2=glacial ice.
          tt = 168. + 0.49*t89
          scat = t23 - t89
          sc31 = t23 - t31
          sc50 = t31 - t50
          par  = t50 - t53

          ! re-frozen snow
          if ( t89 < 255.  .and. &
              scat < sc31        ) then
            scat = sc31
          end if

          ! identify glacial ice
          if ( scat <   3.  .and. &
              t23   < 215.        ) then
            snow(i) = 2
          end if
          if ( scat >= 3.       ) then
            snow(i) = 1
          else
            snow(i) = 0
          end if

          ! remove precipitation
          if ( t23 >= 262.  .or. &
              t23  >= tt           ) then
            snow(i) = 0
          end if
          if ( df3 <= 0.35         ) then    ! remove deserts
            snow(i) = 0
          end if

          ! high elevation deserts
          if ( scat <  15.  .and. &
              sc31  <   3.  .and. &
              par   >   2.        ) then
            snow(i) = 0
          end if

          ! remove frozen ground
          if ( scat <   9.  .and. &
              sc31  <   3.  .and. &
              sc50  <   0.        ) then
            snow(i) = 0
          end if

        end if

      end if

      if ( mwbg_DEBUG .and. i <= 100 ) then
        print *, 'GRODY: i,tb23(i),tb31(i),tb50(i),tb89(i),pangl(i),plat(i), &
                  ilansea(i) = ', &
                  i,tb23(i),tb31(i),tb50(i),tb89(i),pangl(i),plat(i), &
                  ilansea(i)
        print *, 'GRODY: ier(i),ice(i),tpw(i),clwObs(i),clwFG(i),rain(i),snow(i)=', &
                  ier(i),ice(i),tpw(i),clwObs(i),clwFG(i),rain(i),snow(i)
      end if

    end do loopObsGrody

  end subroutine GRODY

  !------------------------------------------------------------------------------------
  ! bennartz
  !------------------------------------------------------------------------------------
  subroutine bennartz (ier, knt, tb89, tb150, tb89FG, tb150FG, pangl, ktermer, &
                       scatl, scatwObs, scatwFG, clwObs, clwFG)

    !:Purpose: Compute the following parameters using 2 AMSU-B channels:
    !          - scattering index (over land and ocean).*
    !          The two channels used are: 89Ghz, 150Ghz.
    !          REGERENCES     Bennartz, R., A. Thoss, A. Dybbroe and D. B. Michelson, 
    !                         1999: Precipitation Analysis from AMSU, Nowcasting SAF, 
    !                         Swedish Meteorologicali and Hydrological Institute, 
    !                         Visiting Scientist Report, November 1999.
    implicit none
    ! arguments: 
    integer, intent(out) :: ier(knt)              ! error return code:
    integer, intent(in)  :: knt              ! number of points to process
    real,    intent(in)  :: tb89(:)          ! 89Ghz AMSU-B brightness temperature (K)
    real,    intent(in)  :: tb150(:)         ! 150Ghz AMSU-B brightness temperature (K)
    real,    intent(in)  :: tb89FG(:)        ! 89Ghz AMSU-B brightness temperature from background (K)
    real,    intent(in)  :: tb150FG(:)       ! 150Ghz AMSU-B brightness temperature from background (K)
    real,    intent(in)  :: pangl(:)         !  satellite zenith angle (deg.)
    integer, intent(in)  :: ktermer(:)       ! land/sea indicator (0=land;1=ocean)
    real,    intent(out) :: scatl(:)         ! scattering index over land
    real,    intent(out) :: scatwObs(:)      ! scattering index over water from observation
    real,    intent(out) :: scatwFG(:)       ! scattering index over water from background
    real,    intent(out) :: clwObs(:)        ! obs cloud liquid water content (not computed for NOW)
    real,    intent(out) :: clwFG(:)         ! first guess cloud liquid water content (not computed for NOW)

    !     Notes: In the case where an output parameter cannot be calculated, the
    !     value of this parameter is to to the missing value, i.e. -99.
    ! Locals: 
    integer :: i 

    ! 1) Initialise output parameters
    do i = 1, knt 
      scatl(i) = mwbg_realMissing
      scatwObs(i) = mwbg_realMissing
      scatwFG(i) = mwbg_realMissing
      clwObs(i) = mwbg_realMissing
      clwFG(i) = mwbg_realMissing
    end do

    ! 2) Validate input parameters
    do i = 1, knt
      if (tb89(i)      < 120.  .or.     &
          tb89(i)     > 350.  .or.     &
          tb150(i)    < 120.  .or.     &
          tb150(i)    > 350.  .or.     & 
          pangl(i)    < -90.  .or.     &
          pangl(i)    >  90.  .or.     & 
          ktermer(i)  <   0   .or.     &
          ktermer(i)  >   1) then
          ier(i) = 1        
      else
          ier(i) = 0      
      end if 
    enddo

    ! 3) Compute parameters
    do i = 1, knt 
      if (ier(i) == 0) then
        if (ktermer(i) == 1) then
          scatwObs(i) = (tb89(i) - tb150(i)) - (-39.2010 + 0.1104 * pangl(i))
          scatwFG(i) = (tb89FG(i) - tb150FG(i)) - (-39.2010 + 0.1104 * pangl(i))
        else
          scatl(i) = (tb89(i) - tb150(i)) - (0.158 + 0.0163 * pangl(i))
        endif
      else if ((ier(i) /= 0) .and. (i <= 100 )) then 
        print *, 'bennartz: input Parameters are not all valid: '
        print *, 'bennartz: i, tb89(i), tb150(i), pangl(i), ktermer(i) = ', &
                            i, tb89(i), tb150(i), pangl(i), ktermer(i)
        print *, 'bennartz: ier(i), scatl(i), scatwObs(i), scatwFG(i)=', &
                            ier(i), scatl(i), scatwObs(i), scatwFG(i)
      endif
    end do 

  end subroutine bennartz

  !--------------------------------------------------------------------------
  ! atmsMwhs2Test1Flagbit7Check
  !--------------------------------------------------------------------------
  subroutine atmsMwhs2Test1Flagbit7Check (itest, KCANO, KMARQ, KNOSAT, ICHECK, KNO, KNT, STNID,  B7CHCK)

    !:Purpose:               1) test 1: Check flag bit 7 on from the first bgckAtms/bgckMwhs2 program
    !                        Includes observations flagged for cloud liquid water, scattering index,
    !                        dryness index plus failure of several QC checks.


    ! Arguments
    integer,     intent(in)                                   :: itest(:)                 ! test number
    integer,     intent(in)                                   :: KCANO(KNO,KNT)           ! observations channels
    integer,     intent(inout)                                :: KMARQ(KNO,KNT)           ! observations channels
    integer,     intent(in)                                   :: KNOSAT                   ! numero de satellite (i.e. indice) 
    integer,     intent(in)                                   :: KNO                      ! nombre de canaux des observations 
    integer,     intent(in)                                   :: KNT                      ! nombre de tovs
    character *9,intent(in)                                   :: STNID                    ! identificateur du satellite
    integer,     intent(inout)                                :: B7CHCK(KNO,KNT) 
    integer,     intent(inout)                                :: ICHECK(KNO,KNT)          ! indicateur du QC par canal
    ! Locals
    integer                                                   :: nDataIndex
    integer                                                   :: nChannelIndex
    integer                                                   :: testIndex
    integer                                                   :: IBIT 

    testIndex = 1
    if ( itest(testIndex) == 1 ) then
      do nDataIndex=1,KNT
        do nChannelIndex=1,KNO
          IBIT = AND(KMARQ(nChannelIndex,nDataIndex), 2**7)
          if ( IBIT /= 0  ) then
            ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
            B7CHCK(nChannelIndex,nDataIndex) = 1
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
            rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                 rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
            if ( mwbg_debug ) then
              write(*,*)STNID(2:9),' first bgckAtms/bgckMwhs2 program REJECT.', &
                        'CHANNEL=', KCANO(nChannelIndex,nDataIndex), &
                        ' IMARQ= ',KMARQ(nChannelIndex,nDataIndex)
            end if
          end if
        end do
      end do
    end if

  end subroutine atmsMwhs2Test1Flagbit7Check

  !--------------------------------------------------------------------------
  ! atmsMwhs2Test2TopographyCheck
  !--------------------------------------------------------------------------
  subroutine atmsMwhs2Test2TopographyCheck(itest, KCANO, KNOSAT, KNO, KNT, STNID, MTINTRP, &
                                           KMARQ, ICHTOPO, MXTOPO, ZCRIT, B7CHCK, ICHECK)

    !:Purpose:               1) test 2: Topography check (partial)

    ! Arguments
    integer,     intent(in)                                   :: itest(:)                       ! test number
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
    integer,     intent(inout)                                :: B7CHCK(KNO,KNT)
    ! Locals
    integer                                                   :: nDataIndex
    integer                                                   :: nChannelIndex
    integer                                                   :: testIndex
    integer                                                   :: INDXTOPO 

    testIndex = 2
    if ( itest(testIndex) == 1 ) then
      do nDataIndex=1,KNT
        do nChannelIndex=1,KNO
          INDXTOPO = ISRCHEQI(ICHTOPO, MXTOPO, KCANO(nChannelIndex,nDataIndex))
          if ( INDXTOPO > 0 ) then
            if ( MTINTRP(nDataIndex) >= ZCRIT(INDXTOPO) ) then
              ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**18)
              rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                   rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
              if ( B7CHCK(nChannelIndex,nDataIndex) == 0 ) then
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

  end subroutine atmsMwhs2Test2TopographyCheck

  !--------------------------------------------------------------------------
  ! atmsMwhs2Test3UncorrectedTbCheck
  !--------------------------------------------------------------------------
  subroutine atmsMwhs2Test3UncorrectedTbCheck(itest, KCANO, KNOSAT, KNO, KNT, STNID, RESETQC, &
                                         KMARQ, B7CHCK, ICHECK)

    !:Purpose:                       Test 3: Uncorrected Tb check (single)
    !                                Uncorrected datum (flag bit #6 off). 
    !                                In this case switch bit 11 ON.

    ! Arguments
    integer,     intent(in)                                   :: itest(:)                       ! test number
    integer,     intent(in)                                   :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)                                   :: KNOSAT                         ! numero de satellite (i.e. indice) 
    integer,     intent(in)                                   :: KNO                            ! nombre de canaux des observations 
    integer,     intent(in)                                   :: KNT                            ! nombre de tovs
    character *9,intent(in)                                   :: STNID                          ! identificateur du satellite
    logical,     intent(in)                                   :: RESETQC                        ! resetqc logical
    integer,     intent(inout)                                :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)                                :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(inout)                                :: B7CHCK(KNO,KNT)
    ! Locals
    integer                                                   :: nDataIndex
    integer                                                   :: nChannelIndex
    integer                                                   :: testIndex
    integer                                                   :: max 
    integer                                                   :: IBIT 

 
    if (.not. RESETQC) then
      testIndex = 3
      if ( itest(testIndex) == 1 ) then
        do nDataIndex=1,KNT
          do nChannelIndex=1,KNO
            IBIT = AND(KMARQ(nChannelIndex,nDataIndex), 2**6)
            if ( IBIT == 0  ) then
              ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
              KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**11)
              rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                 rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
              if ( B7CHCK(nChannelIndex,nDataIndex) == 0 ) then
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

  end subroutine atmsMwhs2Test3UncorrectedTbCheck

  !--------------------------------------------------------------------------
  ! atmsTest4RogueCheck
  !--------------------------------------------------------------------------
  subroutine atmsTest4RogueCheck(itest, KCANO, KNOSAT, KNO, KNT, STNID, ROGUEFAC, TOVERRST, &
                                 clwThreshArr, useStateDepSigmaObs, sigmaObsErr, waterobs, &
                                 PTBOMP, clwObs, clwFG, IDENTF, MXSFCREJ, ISFCREJ, ICH2OMPREJ, &
                                 MXCH2OMPREJ, KMARQ, B7CHCK, ICHECK)

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
    real,        intent(in)              :: ROGUEFAC(:)                    ! rogue factor 
    real(8),     intent(in)              :: TOVERRST(:,:)                  ! erreur totale TOVs
    logical,     intent(in)              :: useStateDepSigmaObs(:,:)       ! if using state dependent obs error
    real(8),     intent(in)              :: clwThreshArr(:,:,:)            ! cloud threshold array
    real(8),     intent(in)              :: sigmaObsErr(:,:,:)             ! sigma obs error
    logical,     intent(in)              :: waterobs(:)                    ! open water obs
    real,        intent(in)              :: PTBOMP(KNO,KNT)                ! radiance o-p 
    real,        intent(in)              :: clwObs(:)
    real,        intent(in)              :: clwFG(:)
    integer,     intent(in)              :: IDENTF(KNT)                    ! data flag ident  
    integer,     intent(in)              :: MXSFCREJ                       ! cst 
    integer,     intent(in)              :: ISFCREJ(MXSFCREJ)
    integer,     intent(in)              :: MXCH2OMPREJ                    ! cst 
    integer,     intent(in)              :: ICH2OMPREJ(MXCH2OMPREJ)
    integer,     intent(inout)           :: KMARQ(KNO,KNT)                 ! marqueur de radiance 
    integer,     intent(inout)           :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(inout)           :: B7CHCK(KNO,KNT)

    ! Locals
    integer                              :: channelval
    integer                              :: nDataIndex
    integer                              :: nChannelIndex
    integer                              :: testIndex
    integer                              :: INDXCAN 
    real                                 :: XCHECKVAL
    real                                 :: clwThresh1 
    real                                 :: clwThresh2
    real                                 :: sigmaThresh1 
    real                                 :: sigmaThresh2
    real                                 :: sigmaObsErrUsed  
    real                                 :: clwObsFGaveraged 
    logical                              :: SFCREJCT
    logical                              :: CH2OMPREJCT
    logical                              :: IBIT 

    testIndex = 4
    if ( itest(testIndex) == 1 ) then
      do nDataIndex=1,KNT
        SFCREJCT = .FALSE.
        CH2OMPREJCT = .FALSE.
        do nChannelIndex=1,KNO
          channelval = KCANO(nChannelIndex,nDataIndex)
          ! using state-dependent obs error only over water.
          ! obs over sea-ice will be rejected in test 15.
          if ( tvs_mwAllskyAssim .and. useStateDepSigmaObs(channelval,KNOSAT) &
                .and. waterobs(nDataIndex) ) then
            clwThresh1 = clwThreshArr(channelval,KNOSAT,1)
            clwThresh2 = clwThreshArr(channelval,KNOSAT,2)
            sigmaThresh1 = sigmaObsErr(channelval,KNOSAT,1)
            sigmaThresh2 = sigmaObsErr(channelval,KNOSAT,2)
            clwObsFGaveraged = 0.5 * (clwObs(nDataIndex) + clwFG(nDataIndex))
            if ( clwObs(nDataIndex) == mwbg_realMissing ) then
              sigmaObsErrUsed = MPC_missingValue_R4
            else
              sigmaObsErrUsed = calcStateDepObsErr_r4(clwThresh1,clwThresh2,sigmaThresh1, &
                                                        sigmaThresh2,clwObsFGaveraged)
            end if
          else
            sigmaObsErrUsed = TOVERRST(channelval,KNOSAT)
          end if
          ! For sigmaObsErrUsed=MPC_missingValue_R4 (clwObs=mwbg_realMissing
          ! in all-sky mode), the observation is flagged for rejection in 
          ! mwbg_reviewAllCritforFinalFlagsAtms.
          XCHECKVAL = ROGUEFAC(channelval) * sigmaObsErrUsed
          if ( PTBOMP(nChannelIndex,nDataIndex)      /= mwbg_realMissing    .and. &
              ABS(PTBOMP(nChannelIndex,nDataIndex))  >= XCHECKVAL .and. &
              sigmaObsErrUsed /= MPC_missingValue_R4 ) then
            ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**16)
            rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) =  &
               rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) + 1
            if ( B7CHCK(nChannelIndex,nDataIndex) == 0 ) then
              rejectionCodArray2(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                 rejectionCodArray2(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1                 
            end if
            if ( mwbg_debug ) then
              write(*,*)STNID(2:9),'ROGUE CHECK REJECT.NO.', &
                     ' OBS = ',nDataIndex, &
                     ' CHANNEL= ',KCANO(nChannelIndex,nDataIndex), &
                     ' CHECK VALUE= ',XCHECKVAL, &
                     ' TBOMP= ',PTBOMP(nChannelIndex,nDataIndex), &
                     ' TOVERRST= ',TOVERRST(channelval,KNOSAT)
            end if
            if ( channelval == 1 .or. &
                 channelval == 2 .or. &
                 channelval == 3    ) then
              SFCREJCT = .TRUE.
            end if
          end if
          if ( channelval == 17 .and. PTBOMP(nChannelIndex,nDataIndex) /= mwbg_realMissing .and. &
              ABS(PTBOMP(nChannelIndex,nDataIndex)) > 5.0 ) then
            CH2OMPREJCT = .TRUE.
          end if
        end do

        if ( SFCREJCT ) then
          do nChannelIndex=1,KNO
            INDXCAN = ISRCHEQI (ISFCREJ,MXSFCREJ,KCANO(nChannelIndex,nDataIndex))
            if ( INDXCAN /= 0 )  then
              if ( ICHECK(nChannelIndex,nDataIndex) /= testIndex ) then
                ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
                KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
                KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**16)
                rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                        rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
                if ( B7CHCK(nChannelIndex,nDataIndex) == 0 ) then
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
        if ( CH2OMPREJCT .and. (IBIT /= 0) ) then
          do nChannelIndex=1,KNO
            INDXCAN = ISRCHEQI (ICH2OMPREJ,MXCH2OMPREJ,KCANO(nChannelIndex,nDataIndex))
            if ( INDXCAN /= 0 )  then
              if ( ICHECK(nChannelIndex,nDataIndex) /= testIndex ) then
                ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
                KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
                KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**16)
                rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                        rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
                if ( B7CHCK(nChannelIndex,nDataIndex) == 0 ) then
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
  ! Mwhs2Test4RogueCheck
  !--------------------------------------------------------------------------
  subroutine Mwhs2Test4RogueCheck(itest, KCANO, KNOSAT, KNO, KNT, STNID, ROGUEFAC, TOVERRST, &
                                  clwThreshArr, useStateDepSigmaObs, sigmaObsErr, waterobs, &
                                  PTBOMP, clwObs, clwFG, IDENTF, ICH2OMPREJ, &
                                  MXCH2OMPREJ, KMARQ, B7CHCK, ICHECK)

    !:Purpose:                         test 4: "Rogue check" for (O-P) Tb residuals out of range.
    !                                  (single/full).
    !                                  Also, over WATER remove CH.10-15 if CH.10 |O-P|>5K (full)
    !                                  Les observations, dont le residu (O-P)
    !                                  depasse par un facteur (roguefac) l'erreur totale des TOVS.
    !                                  OVER OPEN WATER
    !                                  ch. 10 Abs(O-P) > 5K produces rejection of all ATMS amsub channels 10-15.

    ! Arguments
    integer,     intent(in)              :: itest(:)                       ! test number
    integer,     intent(in)              :: KCANO(KNO,KNT)                 ! observations channels
    integer,     intent(in)              :: KNOSAT                         ! numero de satellite (i.e. indice)
    integer,     intent(in)              :: KNO                            ! nombre de canaux des observations
    integer,     intent(in)              :: KNT                            ! nombre de tovs
    character *9,intent(in)              :: STNID                          ! identificateur du satellite
    real,        intent(in)              :: ROGUEFAC(:)                    ! rogue factor
    real(8),     intent(in)              :: TOVERRST(:,:)                  ! erreur totale TOVs
    logical,     intent(in)              :: useStateDepSigmaObs(:,:)       ! if using state dependent obs error
    real(8),     intent(in)              :: clwThreshArr(:,:,:)            ! cloud threshold array
    real(8),     intent(in)              :: sigmaObsErr(:,:,:)             ! sigma obs error
    logical,     intent(in)              :: waterobs(:)                    ! open water obs
    real,        intent(in)              :: PTBOMP(KNO,KNT)                ! radiance o-p
    real,        intent(in)              :: clwObs(:)
    real,        intent(in)              :: clwFG(:)
    integer,     intent(in)              :: IDENTF(KNT)                    ! data flag ident
    integer,     intent(in)              :: MXCH2OMPREJ                    ! cst
    integer,     intent(in)              :: ICH2OMPREJ(MXCH2OMPREJ)
    integer,     intent(inout)           :: KMARQ(KNO,KNT)                 ! marqueur de radiance
    integer,     intent(inout)           :: ICHECK(KNO,KNT)                ! indicateur du QC par canal
    integer,     intent(inout)           :: B7CHCK(KNO,KNT)

    ! Locals
    integer                              :: channelval
    integer                              :: nDataIndex
    integer                              :: nChannelIndex
    integer                              :: testIndex
    integer                              :: INDXCAN
    real                                 :: XCHECKVAL
    real                                 :: clwThresh1
    real                                 :: clwThresh2
    real                                 :: sigmaThresh1
    real                                 :: sigmaThresh2
    real                                 :: sigmaObsErrUsed
    real                                 :: clwObsFGaveraged
    logical                              :: CH2OMPREJCT
    logical                              :: IBIT

    testIndex = 4
    if ( itest(testIndex) == 1 ) then
      do nDataIndex=1,KNT
        CH2OMPREJCT = .FALSE.
        do nChannelIndex=1,KNO
          channelval = KCANO(nChannelIndex,nDataIndex)
          ! using state-dependent obs error only over water.
          ! obs over sea-ice will be rejected in test 15.
          if ( tvs_mwAllskyAssim .and. useStateDepSigmaObs(channelval,KNOSAT) &
                .and. waterobs(nDataIndex) ) then
            clwThresh1 = clwThreshArr(channelval,KNOSAT,1)
            clwThresh2 = clwThreshArr(channelval,KNOSAT,2)
            sigmaThresh1 = sigmaObsErr(channelval,KNOSAT,1)
            sigmaThresh2 = sigmaObsErr(channelval,KNOSAT,2)
            clwObsFGaveraged = 0.5 * (clwObs(nDataIndex) + clwFG(nDataIndex))
            if ( clwObs(nDataIndex) == mwbg_realMissing ) then
              sigmaObsErrUsed = MPC_missingValue_R4
            else
              sigmaObsErrUsed = calcStateDepObsErr_r4(clwThresh1,clwThresh2,sigmaThresh1, &
                                                        sigmaThresh2,clwObsFGaveraged)
            end if
          else
            sigmaObsErrUsed = TOVERRST(channelval,KNOSAT)
          end if
          ! For sigmaObsErrUsed=MPC_missingValue_R4 (clwObs=mwbg_realMissing
          ! in all-sky mode), the observation is flagged for rejection in
          ! mwbg_reviewAllCritforFinalFlagsMwhs2.
          XCHECKVAL = ROGUEFAC(channelval) * sigmaObsErrUsed
          if ( PTBOMP(nChannelIndex,nDataIndex)      /= mwbg_realMissing    .and. &
              ABS(PTBOMP(nChannelIndex,nDataIndex))  >= XCHECKVAL .and. &
              sigmaObsErrUsed /= MPC_missingValue_R4 ) then
            ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
            KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**16)
            rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) =  &
               rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) + 1
            if ( B7CHCK(nChannelIndex,nDataIndex) == 0 ) then
              rejectionCodArray2(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                 rejectionCodArray2(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
            end if
            if ( mwbg_debug ) then
              write(*,*)STNID(2:9),'ROGUE CHECK REJECT.NO.', &
                     ' OBS = ',nDataIndex, &
                     ' CHANNEL= ',KCANO(nChannelIndex,nDataIndex), &
                     ' CHECK VALUE= ',XCHECKVAL, &
                     ' TBOMP= ',PTBOMP(nChannelIndex,nDataIndex), &
                     ' TOVERRST= ',TOVERRST(channelval,KNOSAT)
            end if
          end if
          if ( channelval == 10 .and. PTBOMP(nChannelIndex,nDataIndex) /= mwbg_realMissing .and. &
              ABS(PTBOMP(nChannelIndex,nDataIndex)) > 5.0 ) then
            CH2OMPREJCT = .TRUE.
          end if
        end do

        !    Channels 10-15 are rejected if, for ch10 ABS(O-P) > 5K
        !    Apply over open water only (bit 0 ON in QC integer identf).
        !    Only apply if obs not rejected in this test already.
        IBIT = AND(IDENTF(nDataIndex), 2**0)
        if ( CH2OMPREJCT .and. (IBIT /= 0) ) then
          do nChannelIndex=1,KNO
            INDXCAN = ISRCHEQI (ICH2OMPREJ,MXCH2OMPREJ,KCANO(nChannelIndex,nDataIndex))
            if ( INDXCAN /= 0 )  then
              if ( ICHECK(nChannelIndex,nDataIndex) /= testIndex ) then
                ICHECK(nChannelIndex,nDataIndex) = MAX(ICHECK(nChannelIndex,nDataIndex),testIndex)
                KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**9)
                KMARQ(nChannelIndex,nDataIndex) = OR(KMARQ(nChannelIndex,nDataIndex),2**16)
                rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                        rejectionCodArray(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
                if ( B7CHCK(nChannelIndex,nDataIndex) == 0 ) then
                  rejectionCodArray2(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT) = &
                     rejectionCodArray2(testIndex,KCANO(nChannelIndex,nDataIndex),KNOSAT)+ 1
                end if
              end if
            end if
          end do
        end if

      end do
    end if

  end subroutine Mwhs2Test4RogueCheck

  !--------------------------------------------------------------------------
  ! atmsMwhs2Test5ChannelSelectionUsingIutilst
  !--------------------------------------------------------------------------
  subroutine atmsMwhs2Test5ChannelSelectionUsingIutilst(itest, KCANO, KNOSAT, KNO, KNT, STNID, &
                                                   IUTILST, KMARQ, ICHECK)

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
    integer,     intent(in)               :: IUTILST(:,:)          !  channsl selection
    integer,     intent(inout)            :: KMARQ(KNO,KNT)                 ! marqueur de radiance 

    ! Locals
    integer                               :: channelval
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex

    testIndex = 5
    if ( itest(testIndex) == 1 ) then
      do nDataIndex=1,KNT
        do nChannelIndex=1,KNO
           channelval = KCANO(nChannelIndex,nDataIndex)
           if ( IUTILST(channelval,KNOSAT) == 0 ) then
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

  end subroutine atmsMwhs2Test5ChannelSelectionUsingIutilst

  !--------------------------------------------------------------------------
  ! mwbg_tovCheckAtms 
  !--------------------------------------------------------------------------
  subroutine mwbg_tovCheckAtms(TOVERRST, clwThreshArr, sigmaObsErr, useStateDepSigmaObs, &
                               IUTILST, zlat, zlon, ilq, itt, zenith, qcflag2, qcflag1, &
                               ICANO, ztb, biasCorr, ZOMP, ICHECK, KNO, KNT, KNOSAT, IDENT, &
                               ISCNPOS, MTINTRP, globMarq, IMARQ, clwObs, clwFG, riwv, &
                               STNID, RESETQC)


    !:Purpose:                   Effectuer le controle de qualite des radiances tovs.
    !

    implicit none
    !Arguments
    integer, intent(in)              :: IUTILST(:,:)           ! channel Selection using array IUTILST(chan,sat)
    !                                                            IUTILST = 0 (blacklisted)
    !                                                            1 (assmilate)
    !                                                            2 (assimilate over open water only)

    real(8), intent(in)              :: TOVERRST(:,:)          ! l'erreur totale des TOVS
    real(8), intent(in)              :: clwThreshArr(:,:,:)
    real(8), intent(in)              :: sigmaObsErr(:,:,:)
    logical, intent(in)              :: useStateDepSigmaObs(:,:) ! if using state dependent obs error
    integer, intent(in)              :: KNO                    ! nombre de canaux des observations
    integer, intent(in)              :: KNT                    ! nombre de tovs
    real,    intent(in)              :: zlat(:)
    real,    intent(in)              :: zlon(:)
    integer, intent(in)              :: ilq(:)
    integer, intent(in)              :: itt(:)
    real,    intent(inout)           :: zenith(:)
    integer, intent(in)              :: qcflag2(:)
    integer, intent(in)              :: qcflag1(:,:)
    integer, intent(inout)           :: globMarq(:)
    integer, intent(in)              :: ISCNPOS(KNT)           ! position sur le "scan"
    integer, intent(in)              :: ICANO(:)               ! canaux des observations
    integer, intent(in)              :: KNOSAT                 ! numero de satellite (i.e. indice)
    real, intent(inout)              :: ztb(:)                 ! radiances
    real, intent(in)                 :: biasCorr(:)            ! correction aux radiances
    real, intent(in)                 :: zomp(:)                ! residus (o-p)
    real, intent(in)                 :: MTINTRP(KNT)           ! topographie du modele
    integer, allocatable, intent(out):: IDENT(:)               ! flag to identify all obs pts in report
    !                                                            as being over land/ice, cloudy, bad IWV
    character *9, intent(in)         :: STNID                  ! identificateur du satellite
    logical, intent(in)              :: RESETQC                ! reset du controle de qualite?
    integer,allocatable, intent(out) :: ICHECK(:,:)            ! indicateur controle de qualite tovs par canal 
    !                                                            =0, ok,
    !                                                            >0, rejet,
    integer, intent(inout)           :: IMARQ(:)               ! marqueurs des radiances
    !                                                          satellite, critere et par canal
    !                                                             (chech n2) par satellite, critere et par canal
    real, allocatable, intent(out)   :: clwObs(:)
    real, allocatable, intent(out)   :: clwFG(:)
    real, allocatable, intent(out)   :: riwv(:)

    !locals
    real                             :: PTBOMP(KNO,KNT)
    integer                          :: KCANO(KNO,KNT)
    integer                          :: KMARQ(KNO,KNT)
    integer, allocatable             :: lsq(:)
    integer, allocatable             :: trn(:)
    integer,allocatable              :: KCHKPRF(:)

    logical, allocatable             :: waterobs(:)
    logical, allocatable             :: grossrej(:)
    logical                          :: reportHasMissingTb
    logical, allocatable             :: lqc(:,:)
    logical, allocatable             :: cloudobs(:)
    logical, allocatable             :: iwvreject(:)
    logical, allocatable             :: precipobs(:)
    real, allocatable                :: zdi(:)
    real, allocatable                :: scatec(:)
    real, allocatable                :: scatbg(:)
    real, allocatable                :: SeaIce(:)

    integer, parameter               :: maxScanAngleAMSU = 96
    integer, parameter               :: MXSFCREJ   = 8
    integer, parameter               :: MXCH2OMPREJ= 6
    integer, parameter               :: MXTOPO     = 5
    integer, parameter               :: MXCLWREJ   = 6
    integer                          :: iRej
    integer                          :: iNumSeaIce
    integer                          :: JI
    integer                          :: JJ
    integer                          :: kk
    integer                          :: INDX
    integer                          :: ISFCREJ(MXSFCREJ)
    integer                          :: ICH2OMPREJ(MXCH2OMPREJ)
    integer                          :: B7CHCK(KNO,KNT)
    real, allocatable                :: ROGUEFAC(:)
    real                             :: ZCRIT(MXTOPO)
    integer                          :: ITEST(mwbg_maxNumTest)
    integer                          :: chanFlaggedForAllskyGenCoeff(MXCLWREJ)
    integer                          :: ICHTOPO(MXTOPO)
    logical, save                    :: LLFIRST = .true.
    integer, save                    :: numReportWithMissingTb
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

    call utl_reAllocate(ROGUEFAC, KNO+tvs_channelOffset(KNOSAT))
    ROGUEFAC(:) = (/2.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 4.0, &
                    4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 2.0, &
                    2.0, 4.0, 4.0, 4.0, 4.0, 4.0/)
    if ( tvs_mwAllskyAssim ) ROGUEFAC(1:3) = 3.0

    ! Channel sets for rejection in test 9 
    ! These LT channels are rejected if O-P fails rogue check for window ch. 1, 2, or 3
    ISFCREJ(:) = (/1, 2, 3, 4, 5, 6, 16, 17/)
    !   These AMSU-B channels are rejected if ch. 17 O-P fails rogue check over OPEN WATER only    
    ICH2OMPREJ(:) = (/17, 18, 19, 20, 21, 22/)

    !  Data for TOPOGRAPHY CHECK
    !   Channel AMSUA-6 (atms ch 7) is rejected for topography  >  250m.
    !   Channel AMSUA-7 (atms ch 8) is rejected for topography  > 2000m.
    !   Channel AMSUB-3 (atms ch 22) is rejected for topography > 2500m.
    !                    atms ch 21  is rejected for topography > 2250m.
    !   Channel AMSUB-4 (atms ch 20) is rejected for topography > 2000m.
    ICHTOPO(:) = (/7, 8, 20, 21, 22/)
    ZCRIT(:) = (/250., 2000., 2000., 2250., 2500./)

    !  Test selection (0=skip test, 1=do test)
    !              1  2  3  4  5
    ITEST(:)  = 0
    ITEST(1:5) = (/1, 1, 1, 1, 1/)

    ! Channels excluded from gen_bias_corr in all-sky mode
    chanFlaggedForAllskyGenCoeff(:) = (/ 1, 2, 3, 4, 5, 6/)

    ! Initialisation, la premiere fois seulement!
    if (LLFIRST) then
      numReportWithMissingTb = 0
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
    call mwbg_landIceMaskAtms(KNT, zlat, zlon, lsq, trn, waterobs)

    !###############################################################################
    ! STEP 2 ) Check for values of TB that are missing or outside physical limits.
    !###############################################################################

    call mwbg_grossValueCheck(KNT, KNO, ztb, biasCorr, 50., 380., grossrej)

    !###############################################################################
    ! STEP 3 ) Preliminary QC checks --> set lqc(KNT,KNO)=.true.
    !          for data that fail QC
    !###############################################################################

    call mwbg_firstQcCheckAtms(zenith, ilq, itt, zlat, zlon, ztb, ISCNPOS, &
                               KNO, KNT, lqc, grossrej, lsq, trn, qcflag1, &
                               qcflag2, ICANO, reportHasMissingTb)

    if ( reportHasMissingTb ) numReportWithMissingTb = numReportWithMissingTb + 1
    !  Exclude problem points from further calculations
    do kk = 1,KNT
      if ( COUNT(lqc(kk,:)) == KNO ) grossrej(kk) = .true.
    end do

    !###############################################################################
    ! STEP 4 ) mwbg_nrlFilterAtms returns clwObs, clwFG, scatec, scatbg and also does sea-ice
    !          detection missing value for  clwObs, scatec, scatbg  is -99.0 (e.g. over
    !          land or sea-ice).Sets trn=0 (sea ice) for points where retrieved SeaIce
    !          >=0.55. Does nothing if trn=0 (sea ice) and retrieved SeaIce<0.55.
    !###############################################################################
 
    call mwbg_nrlFilterAtms(KNT, KNO, ztb, zomp, biasCorr, zenith, zlat, lsq, trn, waterobs, &
                            grossrej, clwObs, clwFG, scatec, scatbg, iNumSeaIce, iRej, SeaIce)
    seaIcePointNum = seaIcePointNum + iNumSeaIce
    clwMissingPointNum = clwMissingPointNum + iRej

    !###############################################################################
    ! STEP 5 ) Apply NRL cloud filter, scattering index and sea-ice detection algorithms
    !          to OPEN WATER (waterobs=true) points.
    ! Points with SeaIce>0.55 are set to sea-ice points (waterobs --> false)
    !###############################################################################

    call mwbg_flagDataUsingNrlCritAtms(KNT, KNO, ztb, biasCorr, clwObs, scatec, scatbg, &
                                       SeaIce, grossrej, waterobs, mwbg_useUnbiasedObsForClw, &
                                       iwvreject, cloudobs, precipobs, cldcnt , ident, riwv, zdi)

    !###############################################################################
    ! STEP 6 ) ! Review all the checks previously made to determine which obs are to be
    !            accepted for assimilation and which are to be flagged for exclusion
    !            (IMARQ).
    !            grossrej()  = .true. if any channel had a gross error at the point
    !            cloudobs()  = .true. if CLW > clw_atms_nrl_LTrej (0.175) or precipobs
    !            precipobs() = .true. if precip. detected through NRL scattering indices
    !            waterobs()  = .true. if open water point
    !            iwvreject() = .true. if Mean 183 Ghz [ch. 18-22] Tb < 240K (too dry
    !            for ch.20-22 over land)
    !###############################################################################

    call mwbg_reviewAllCritforFinalFlagsAtms(KNT, KNO, lqc, grossrej, waterobs, &
                                             precipobs, clwObs, clwFG, scatec, scatbg, &
                                             iwvreject, riwv, IMARQ, globMarq, zdi, ident, &
                                             drycnt, landcnt, rejcnt, iwvcnt, pcpcnt, flgcnt, &
                                             MXCLWREJ, chanFlaggedForAllskyGenCoeff, icano)

    !###############################################################################
    ! PART 2 TESTS:
    !###############################################################################

    ! copy the original input 1D array to 2D array. The 2D arrays are used in this s/r.
    call copy1Dimto2DimRealArray(ZOMP, KNO, KNT, PTBOMP)
    call copy1Dimto2DimIntegerArray(ICANO, KNO, KNT, KCANO)
    call copy1Dimto2DimIntegerArray(IMARQ, KNO, KNT, KMARQ)
    ! allocations
    call utl_reAllocate(kchkprf, KNT)
    call utl_reAllocate(icheck, KNO, KNT)
    !  Initialisations
    ICHECK(:,:) = 0
    B7CHCK(:,:) = 0

    if ( RESETQC ) KMARQ(:,:) = 0

    ! 1) test 1: Check flag bit 7 on from the first bgckAtms program
    !  Includes observations flagged for cloud liquid water, scattering index,
    !  dryness index plus failure of several QC checks.
    call atmsMwhs2Test1Flagbit7Check (itest, KCANO, KMARQ, KNOSAT, ICHECK, KNO, KNT, &
                                      STNID,  B7CHCK)

    ! 2) test 2: Topography check (partial)
    call atmsMwhs2Test2TopographyCheck (itest, KCANO, KNOSAT, KNO, KNT, STNID, MTINTRP, &
                                        KMARQ, ICHTOPO, MXTOPO, ZCRIT, B7CHCK, ICHECK)
    ! 3) test 3: Uncorrected Tb check (single)
    !  Uncorrected datum (flag bit #6 off). In this case switch bit 11 ON.
    call atmsMwhs2Test3UncorrectedTbCheck (itest, KCANO, KNOSAT, KNO, KNT, STNID, RESETQC, &
                                           KMARQ, B7CHCK, ICHECK)
    ! 4) test 4: "Rogue check" for (O-P) Tb residuals out of range. (single/full)
    !             Also, over WATER remove CH.17-22 if CH.17 |O-P|>5K (partial)
    !  Les observations, dont le residu (O-P) depasse par un facteur (roguefac)
    !   l'erreur totale des TOVS.
    !  N.B.: a reject by any of the 3 amsua surface channels 1-3 produces the
    !           rejection of ATMS sfc/tropospheric channels 1-6 and 16-17.
    !  OVER OPEN WATER
    !    ch. 17 Abs(O-P) > 5K produces rejection of all ATMS amsub channels 17-22.
    call atmsTest4RogueCheck (itest, KCANO, KNOSAT, KNO, KNT, STNID, ROGUEFAC, TOVERRST, &
                              clwThreshArr, useStateDepSigmaObs, sigmaObsErr, waterobs, &
                              PTBOMP, clwObs, clwFG, IDENT, MXSFCREJ, ISFCREJ, ICH2OMPREJ, &
                              MXCH2OMPREJ, KMARQ, B7CHCK, ICHECK)

    ! 5) test 5: Channel selection using array IUTILST(chan,sat)
    !  IUTILST = 0 (blacklisted)
    !            1 (assmilate)
    call atmsMwhs2Test5ChannelSelectionUsingIutilst(itest, KCANO, KNOSAT, KNO, KNT, STNID, &
                                                    IUTILST, KMARQ, ICHECK)

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

    ! reset global marker flag (55200) and mark it if observtions are rejected
    call resetQcCases(RESETQC, KCHKPRF, globMarq)

    if(mwbg_debug) then
      write(*,*) ' --------------------------------------------------------------- '
      write(*,*) ' Number of BURP file reports where Tb set to mwbg_realMissing  = ', numReportWithMissingTb
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
  ! mwbg_tovCheckMwhs2
  !--------------------------------------------------------------------------
  subroutine mwbg_tovCheckMwhs2(TOVERRST, clwThreshArr, sigmaObsErr, useStateDepSigmaObs, &
                                IUTILST, zlat, zlon, ilq, itt, zenith, &
                                ICANO, ztb, biasCorr, ZOMP, ICHECK, KNO, KNT, KNOSAT, IDENT, &
                                ISCNPOS, MTINTRP, globMarq, IMARQ, clwObs, clwFG, riwv, &
                                STNID, RESETQC, modLSQ, lastHeader)


    !:Purpose:                   Effectuer le controle de qualite des radiances tovs.
    !

    implicit none
    !Arguments
    integer, intent(in)              :: IUTILST(:,:)           ! channel Selection using array IUTILST(chan,sat)
    !                                                            IUTILST = 0 (blacklisted)
    !                                                            1 (assmilate)
    !                                                            2 (assimilate over open water only)

    real(8), intent(in)              :: TOVERRST(:,:)          ! l'erreur totale des TOVS
    real(8), intent(in)              :: clwThreshArr(:,:,:)
    real(8), intent(in)              :: sigmaObsErr(:,:,:)
    logical, intent(in)              :: useStateDepSigmaObs(:,:) ! if using state dependent obs error
    integer, intent(in)              :: KNO                    ! nombre de canaux des observations
    integer, intent(in)              :: KNT                    ! nombre de tovs
    real,    intent(in)              :: zlat(:)
    real,    intent(in)              :: zlon(:)
    integer, intent(in)              :: ilq(:)
    integer, intent(in)              :: itt(:)
    real,    intent(inout)           :: zenith(:)
    integer, intent(inout)           :: globMarq(:)
    integer, intent(in)              :: ISCNPOS(KNT)           ! position sur le "scan"
    integer, intent(in)              :: ICANO(:)               ! canaux des observations
    integer, intent(in)              :: KNOSAT                 ! numero de satellite (i.e. indice)
    real, intent(inout)              :: ztb(:)                 ! radiances
    real, intent(in)                 :: biasCorr(:)            ! correction aux radiances
    real, intent(in)                 :: zomp(:)                ! residus (o-p)
    real, intent(in)                 :: MTINTRP(KNT)           ! topographie du modele
    integer, allocatable, intent(out):: IDENT(:)               ! flag to identify all obs pts in report
    !                                                            as being over land/ice, cloudy, bad IWV
    character *9, intent(in)         :: STNID                  ! identificateur du satellite
    logical, intent(in)              :: RESETQC                ! reset du controle de qualite?
    logical, intent(in)              :: modLSQ                 ! If active, recalculate values for land/sea
                                                               ! qualifier and terrain type based on LG/MG
    logical, intent(in)              :: lastHeader             ! active if last header
    integer,allocatable, intent(out) :: ICHECK(:,:)            ! indicateur controle de qualite tovs par canal
    !                                                            =0, ok,
    !                                                            >0, rejet,
    integer, intent(inout)           :: IMARQ(:)               ! marqueurs des radiances
    !                                                          satellite, critere et par canal
    !                                                             (chech n2) par satellite, critere et par canal
    real, allocatable, intent(out)   :: clwObs(:)
    real, allocatable, intent(out)   :: clwFG(:)
    real, allocatable, intent(out)   :: riwv(:)

    !locals
    real                             :: PTBOMP(KNO,KNT)
    integer                          :: KCANO(KNO,KNT)
    integer                          :: KMARQ(KNO,KNT)
    integer, allocatable             :: lsq(:)
    integer, allocatable             :: trn(:)
    integer,allocatable              :: KCHKPRF(:)

    logical, allocatable             :: waterobs(:)
    logical, allocatable             :: grossrej(:)
    logical                          :: reportHasMissingTb
    logical, allocatable             :: lqc(:,:)
    logical, allocatable             :: cloudobs(:)
    logical, allocatable             :: iwvreject(:)
    logical, allocatable             :: precipobs(:)
    real, allocatable                :: zdi(:)
    real, allocatable                :: scatec(:)
    real, allocatable                :: scatbg(:)
    real, allocatable                :: SeaIce(:)

    integer, parameter               :: maxScanAngleAMSU = 98
    integer, parameter               :: MXSFCREJ   = 8
    integer, parameter               :: MXCH2OMPREJ= 6
    integer, parameter               :: MXTOPO     = 3
    integer, parameter               :: MXCLWREJ   = 6

    integer                          :: iRej
    integer                          :: iNumSeaIce
    integer                          :: JI
    integer                          :: JJ
    integer                          :: kk
    integer                          :: INDX
    integer                          :: ICH2OMPREJ(MXCH2OMPREJ)
    integer                          :: B7CHCK(KNO,KNT)
    real, allocatable                :: ROGUEFAC(:)
    real                             :: ZCRIT(MXTOPO)
    integer                          :: ITEST(mwbg_maxNumTest)
    integer                          :: chanFlaggedForAllskyGenCoeff(MXCLWREJ)
    integer                          :: ICHTOPO(MXTOPO)
    logical, save                    :: LLFIRST = .true.
    integer, save                    :: numReportWithMissingTb
    integer, save                    :: allcnt                          ! Number of Tovs obs
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

    call utl_reAllocate(ROGUEFAC, KNO+tvs_channelOffset(KNOSAT))
    ROGUEFAC(:) = (/2.0, 9.9, 9.9, 9.9, 9.9, 9.9, 9.9, 9.9, &
                    9.9, 2.0, 4.0, 4.0, 4.0, 4.0, 4.0/)
    if ( tvs_mwAllskyAssim ) ROGUEFAC(1:3) = 9.9

    ! Channel sets for rejection in test 9
    !   These AMSU-B channels are rejected if ch. 10 O-P fails rogue check over OPEN WATER only
    ICH2OMPREJ(:) = (/10, 11, 12, 13, 14, 15/)

    !  Data for TOPOGRAPHY CHECK
    !   Channel AMSUB-3 (mwhs2 ch 11) is rejected for topography > 2500m.
    !                   (mwhs2 ch 12) is rejected for topography > 2250m.
    !   Channel AMSUB-4 (mwhs2 ch 13) is rejected for topography > 2000m.
    ICHTOPO(:) = (/11, 12, 13/)
    ZCRIT(:) = (/2500., 2250., 2000./)

    !  Test selection (0=skip test, 1=do test)
    !              1  2  3  4  5
    ITEST(:)  = 0
    ITEST(1:5) = (/1, 1, 1, 1, 1/)

    ! Channels excluded from gen_bias_corr in all-sky mode
    chanFlaggedForAllskyGenCoeff(:) = (/ 10, 11, 12, 13, 14, 15/)

    ! Initialisation, la premiere fois seulement!
    if (LLFIRST) then
      numReportWithMissingTb = 0
      allcnt = 0
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
    call mwbg_landIceMaskMwhs2(KNT, zlat, zlon, lsq, trn, waterobs, ilsmOpt)

    !###############################################################################
    ! STEP 2 ) Check for values of TB that are missing or outside physical limits.
    !###############################################################################
    call mwbg_grossValueCheck(KNT, KNO, ztb, biasCorr, 50., 380., grossrej)

    !###############################################################################
    ! STEP 3 ) Preliminary QC checks --> set lqc(KNT,KNO)=.true.
    !          for data that fail QC
    !###############################################################################
    call mwbg_firstQcCheckMwhs2(zenith, ilq, itt, zlat, zlon, ztb, ISCNPOS, &
                                KNO, KNT, lqc, lsq, trn, ICANO, reportHasMissingTb, modLSQ)

    if ( reportHasMissingTb ) numReportWithMissingTb = numReportWithMissingTb + 1
    !  Exclude problem points from further calculations
    do kk = 1,KNT
      if ( COUNT(lqc(kk,:)) == KNO ) grossrej(kk) = .true.
    end do

    !###############################################################################
    ! STEP 4 ) mwbg_nrlFilterMwhs2 returns clwObs, clwFG, scatec, scatbg and also does sea-ice
    !          detection missing value for  clwObs, scatec, scatbg  is -99.0 (e.g. over
    !          land or sea-ice).Sets trn=0 (sea ice) for points where retrieved SeaIce
    !          >=0.55. Does nothing if trn=0 (sea ice) and retrieved SeaIce<0.55.
    !###############################################################################

    call mwbg_nrlFilterMwhs2(KNT, KNO, ztb, biasCorr, zenith, zlat, lsq, trn, waterobs, &
                             grossrej, clwObs, clwFG, scatec, scatbg, iNumSeaIce, iRej, SeaIce)
    seaIcePointNum = seaIcePointNum + iNumSeaIce
    clwMissingPointNum = clwMissingPointNum + iRej

    !###############################################################################
    ! STEP 5 ) Apply NRL cloud filter, scattering index and sea-ice detection algorithms
    !          to OPEN WATER (waterobs=true) points.
    ! Points with SeaIce>0.55 are set to sea-ice points (waterobs --> false)
    !###############################################################################

    call mwbg_flagDataUsingNrlCritMwhs2(KNT, KNO, ztb, biasCorr, clwObs, scatec, &
                                        SeaIce, grossrej, waterobs, mwbg_useUnbiasedObsForClw, &
                                        iwvreject, cloudobs, precipobs, cldcnt , ident, riwv, zdi)

    !###############################################################################
    ! STEP 6 ) ! Review all the checks previously made to determine which obs are to be
    !            accepted for assimilation and which are to be flagged for exclusion
    !            (IMARQ).
    !            grossrej()  = .true. if any channel had a gross error at the point
    !            cloudobs()  = .true. if CLW > clw_mwhs2_nrl_LTrej (0.175) or precipobs
    !            precipobs() = .true. if precip. detected through NRL scattering indices
    !            waterobs()  = .true. if open water point
    !            iwvreject() = .true. if Mean 183 Ghz [ch. 11-15] Tb < 240K (too dry
    !            for ch.11-13 over land)
    !###############################################################################

    call mwbg_reviewAllCritforFinalFlagsMwhs2(KNT, KNO, lqc, grossrej, trn, waterobs, &
                                              precipobs, clwObs, clwFG, scatec, scatbg, &
                                              iwvreject, riwv, IMARQ, globMarq, zdi, ident, &
                                              allcnt, drycnt, landcnt, rejcnt, iwvcnt, pcpcnt, flgcnt, &
                                              MXCLWREJ, chanFlaggedForAllskyGenCoeff, icano)

    !###############################################################################
    ! PART 2 TESTS:
    !###############################################################################

    ! copy the original input 1D array to 2D array. The 2D arrays are used in this s/r.
    call copy1Dimto2DimRealArray(ZOMP, KNO, KNT, PTBOMP)
    call copy1Dimto2DimIntegerArray(ICANO, KNO, KNT, KCANO)
    call copy1Dimto2DimIntegerArray(IMARQ, KNO, KNT, KMARQ)
    ! allocations
    call utl_reAllocate(kchkprf, KNT)
    call utl_reAllocate(icheck, KNO, KNT)
    !  Initialisations
    ICHECK(:,:) = 0
    B7CHCK(:,:) = 0

    if ( RESETQC ) KMARQ(:,:) = 0

    ! 1) test 1: Check flag bit 7 on from the first bgckMwhs2 program
    !  Includes observations flagged for cloud liquid water, scattering index,
    !  dryness index plus failure of several QC checks.
    call atmsMwhs2Test1Flagbit7Check (itest, KCANO, KMARQ, KNOSAT, ICHECK, KNO, KNT, &
                                      STNID,  B7CHCK)

    ! 2) test 2: Topography check (partial)
    call atmsMwhs2Test2TopographyCheck (itest, KCANO, KNOSAT, KNO, KNT, STNID, MTINTRP, &
                                        KMARQ, ICHTOPO, MXTOPO, ZCRIT, B7CHCK, ICHECK)
    ! 3) test 3: Uncorrected Tb check (single)
    !  Uncorrected datum (flag bit #6 off). In this case switch bit 11 ON.
    call atmsMwhs2Test3UncorrectedTbCheck (itest, KCANO, KNOSAT, KNO, KNT, STNID, RESETQC, &
                                           KMARQ, B7CHCK, ICHECK)
    ! 4) test 4: "Rogue check" for (O-P) Tb residuals out of range. (single/full)
    !             Also, over WATER remove CH.10-15 if CH.10 |O-P|>5K (full)
    !  Les observations, dont le residu (O-P) depasse par un facteur (roguefac)
    !   l'erreur totale des TOVS.
    !  OVER OPEN WATER
    !    ch. 10 Abs(O-P) > 5K produces rejection of all ATMS amsub channels 10-15.
    call Mwhs2Test4RogueCheck (itest, KCANO, KNOSAT, KNO, KNT, STNID, ROGUEFAC, TOVERRST, &
                              clwThreshArr, useStateDepSigmaObs, sigmaObsErr, waterobs, &
                              PTBOMP, clwObs, clwFG, IDENT, ICH2OMPREJ, &
                              MXCH2OMPREJ, KMARQ, B7CHCK, ICHECK)

    ! 5) test 5: Channel selection using array IUTILST(chan,sat)
    !  IUTILST = 0 (blacklisted)
    !            1 (assmilate)
    call atmsMwhs2Test5ChannelSelectionUsingIutilst(itest, KCANO, KNOSAT, KNO, KNT, STNID, &
                                                    IUTILST, KMARQ, ICHECK)

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

    ! reset global marker flag (55200) and mark it if observtions are rejected
    call resetQcCases(RESETQC, KCHKPRF, globMarq)

    if (lastHeader) then
      write(*,*) ' --------------------------------------------------------------- '
      write(*,*) ' Number of obs pts read from BURP file                         = ', allcnt
      write(*,*) ' Number of BURP file reports where Tb set to mwbg_realMissing  = ', numReportWithMissingTb
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
      write(*,*) '       1     Mean 183 Ghz [ch. 11-15] is missing'
      write(*,*) '       2     NRL CLW is missing (over water)'
      write(*,*) '       3     NRL > clw_mwhs2_nrl_LTrej (0.175 kg/m2) (cloudobs)'
      write(*,*) '       4     scatec > Lower Troposphere limit=9 (precipobs)'
      write(*,*) '       5     Mean 183 Ghz [ch. 11-15] Tb < 240K'
      write(*,*) '       6     CLW > clw_mwhs2_nrl_UTrej (0.200 kg/m2)'
      write(*,*) '       7     Dryness Index rejection (for ch. 11)'
      write(*,*) '       8     scatbg > CMC amsu-b limit (land=0,sea=15,ice=40)'
      write(*,*) '       9     Dryness Index rejection (for ch. 12)'
      write(*,*) '      10     Sea ice > 0.55 detected'
      write(*,*) '      11     Gross error in Tb (any chan.) or other QC problem (all channels rejected)'
      write(*,*) ' '
      write(*,*) '   New Element 13209 in BURP file = CLW (kg/m2)'
      write(*,*) '   New Element 13208 in BURP file = Bennartz-Grody Scattering Index'
      write(*,*) '   New Element 25174 in BURP file = IDENT flag'
      write(*,*) ' '
    end if

  end subroutine mwbg_tovCheckMwhs2

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
    integer                  :: ier, irec
    integer                  :: ezqkdef, ezsetopt
    integer                  :: FSTINF,FSTPRM,FCLOS
    integer                  :: FSTLIR,FSTFRM, FNOM, FSTOUV
    integer                  :: NI, NJ, NK, IG1, IG2, IG3, IG4
    integer                  :: IDUM1,IDUM2,IDUM3,IDUM4
    integer                  :: IDUM5,IDUM6,IDUM7,IDUM8
    integer                  :: IDUM9,IDUM10,IDUM11,IDUM12,IDUM13
    integer                  :: IDUM14,IDUM15,IDUM16,IDUM17,IDUM18
    character(len=12)        :: ETIKXX
    character(len=4)         :: CLNOMVAR
    character(len=4)         :: NOMVXX
    character(len=2)         :: TYPXX
    character(len=1)         :: GRTYP
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
    readGlaceMask = .True.
    if (instName == 'ATMS') readGlaceMask = .False.
    if(ifFirstCall) then
      IUNGEO = 0
      IER = FNOM(IUNGEO,fileMgLg,'STD+RND+R/O',0)

      ! 3) Lecture des champs geophysiques (MF/MX) du modele
      IER = FSTOUV(IUNGEO,'RND')

      ! TOPOGRAPHIE (MF ou MX).
      !     MX est la topographie filtree avec unites en m2/s2  (geopotential topography).

      IREC = FSTINF(IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1,' ','MX')
      if (IREC .GE. 0) then
        TOPOFACT = 9.80616
        CLNOMVAR = 'MX'
        if(allocated(MT)) deallocate(MT)
        allocate ( MT(NI*NJ), STAT=ier)
        IER = FSTLIR(MT,IUNGEO,NI,NJ,NK,-1,' ',-1,-1,-1, &
           ' ',CLNOMVAR)
      else
        call utl_abort ('bgckMicrowave_mod: ERREUR: LA TOPOGRAPHIE (MF or MX) EST INEXISTANTE')
      end if

      IER = FSTPRM ( IREC, IDUM1, IDUM2, IDUM3, IDUM4, &
          IDUM5, IDUM6, IDUM7, IDUM8, IDUM9, IDUM10,  &
          IDUM11, TYPXX, NOMVXX, ETIKXX, GRTYP, IG1, &
          IG2, IG3, IG4, IDUM12, IDUM13, IDUM14,  &
          IDUM15, IDUM16, IDUM17, IDUM18 )
       write (*,*) ' GRILLE MT : ',grtyp,ni,nj, &
                ig1,ig2,ig3,ig4
      ier  = ezsetopt('INTERP_DEGREE','LINEAR')
      ier  = ezsetopt('EXTRAP_DEGREE','ABORT')
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
        ier  = ezsetopt('EXTRAP_DEGREE','ABORT')
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
        ier  = ezsetopt('EXTRAP_DEGREE','ABORT')
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
          if ( XLON < -180. ) XLON = XLON + 360.
          if ( XLON >  180. ) XLON = XLON - 360.
          if ( XLON <    0. ) XLON = XLON + 360.
           ZLATBOX(boxPointIndex,dataIndex) = XLAT
           ZLONBOX(boxPointIndex,dataIndex) = XLON
         end do
      end do
    end do
    ier = ezsetopt('INTERP_DEGREE','LINEAR')
    ier = gdllsval(gdmt,mtintbox,mt,ZLATBOX,ZLONBOX,boxPointNum*dataNum)
    if (ier < 0) then
      call utl_abort ('bgckMicrowave_mod: ERROR in the interpolation of MT')
    end if
    if(readGlaceMask) then
      ier = gdllsval(gdmg,mgintbox,mg,ZLATBOX,ZLONBOX,boxPointNum*dataNum)
      if (ier < 0) then
        call utl_abort ('bgckMicrowave_mod: ERROR in the interpolation of MG')
      end if
      ier = gdllsval(gdgl,glintbox,gl,ZLATBOX,ZLONBOX,boxPointNum*dataNum)
      if (ier < 0) then
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
      if (mwbg_debug) then
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
      if (mwbg_debug) then
        print *, ' MGINTRP = ', MGINTRP(dataIndex)
        print *, ' MTINTRP = ', MTINTRP(dataIndex)
        print *, ' GLINTRP = ', GLINTRP(dataIndex)
      end if
    end do
  end subroutine mwbg_readGeophysicFieldsAndInterpolate

  !--------------------------------------------------------------------------
  ! mwbg_landIceMaskAtms
  !--------------------------------------------------------------------------
  subroutine mwbg_landIceMaskAtms(npts,zlat,zlon,zlq,ztt,waterobs)
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
    !--------------------------------------------------------------------
    !  Variable Definitions
    !  --------------------
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
    integer, intent(in)                   :: npts
    real,    intent(in)                   :: zlat(:)
    real,    intent(in)                   :: zlon(:)
    integer, intent(out), allocatable     :: zlq(:)
    integer, intent(out), allocatable     :: ztt(:)
    logical, intent(out), allocatable     :: waterobs(:)

    ! Locals:
    logical, save :: firstCall=.true.
    integer, parameter :: mxlat=5,mxlon=5
    integer :: iungeo

    integer :: ier,key
    integer, save :: ni,nj,nk,nilg,njlg
    integer, save :: ig1,ig2,ig3,ig4,ig1lg,ig2lg,ig3lg,ig4lg
    integer :: idum4,idum5,idum6,idum7,idum8,idum9,idum10,idum11
    integer :: idum12,idum13,idum14,idum15,idum16,idum17,idum18

    integer :: indx,ii,jj,kk
    integer :: nlat,nlon

    real, parameter :: pi=3.141592654
    real, parameter :: MGthresh=0.01,LGthresh=0.01
    real, parameter :: rlat_km=40.0,rlon_km=40.0
    real, parameter :: rkm_per_deg=111.195

    real :: xlat,xlatrad,xlon,rii,rjj
    real :: dlat,dlon

    character(len=12) :: etikxx
    character(len=4)  :: nomvxx
    character(len=2)  :: typxx
    character(len=1), save :: grtyp,grtyplg

    logical  :: llg

    ! F90 allocatable arrays:
    real, allocatable, save :: mg(:),lg(:)
    real, allocatable       :: latmesh(:),lonmesh(:)
    real, allocatable       :: mgintob(:),lgintob(:)
    real, allocatable       :: zlatbox(:,:),zlonbox(:,:)
    real, allocatable       :: mgintrp(:),lgintrp(:)

    ! RMNLIB interpolating functions:
    integer :: ezsetopt,ezqkdef
    integer :: gdllsval,gdid,gdidlg

    ! Define FORTRAN FST functions:
    integer, external :: fstinf,fstprm,fstlir,fnom,fclos
    integer, external :: fstouv,fstfrm,fstinl,fstvoi

    integer :: idum1,idum2,idum3

    ! Allocate space for arrays holding values on mesh grid pts.
    call utl_reAllocate(latmesh, mxlat*mxlon)
    call utl_reAllocate(lonmesh, mxlat*mxlon)
    call utl_reAllocate(mgintob, mxlat*mxlon)
    call utl_reAllocate(lgintob, mxlat*mxlon)
    call utl_reAllocate(zlatbox, mxlat*mxlon, npts)
    call utl_reAllocate(zlonbox, mxlat*mxlon, npts)
    call utl_reAllocate(zlq, npts)
    call utl_reAllocate(ztt, npts)
    call utl_reAllocate(waterobs, npts)

    ! Open FST file.
    iungeo = 0
    ier = fnom( iungeo,fileMgLg,'STD+RND+R/O',0 )
    ier = fstouv( iungeo,'RND' )

    if (firstCall) then
      firstCall = .false.

      ! Read MG field.
      key = fstinf(iungeo,ni,nj,nk,-1,' ',-1,-1,-1,' ' ,'MG')
      if ( key <  0 ) then
        call utl_abort('mwbg_landIceMaskAtms: The MG field is MISSING')
      end if

      call utl_reAllocate(mg, ni*nj)

      ier = fstlir(mg,iungeo,ni,nj,nk,-1,' ',-1,-1,-1,' ','MG')

      ier = fstprm(key,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,    &
                   idum9,idum10,idum11,typxx,nomvxx,etikxx,grtyp,ig1,ig2,  &
                   ig3,ig4,idum12,idum13,idum14,idum15,idum16,idum17,      &
                   idum18)

      ! Read LG field. Use GL field as backup.
      ! **CAUTION**: Discontinuities in GL field may cause interpolation problems! LG field is preferable.
      llg = .false.
      key = fstinf(iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ' ,'LG')
      if ( key <  0 ) then
        key = fstinf(iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ' ,'GL')
        if ( key <  0 ) then
          call utl_abort('mwbg_landIceMaskAtms: The LG or GL field is MISSING')
        end if
      else
        llg = .true.
      end if

      call utl_reAllocate(lg, nilg*njlg)

      if ( llg ) then
        ier = fstlir(lg,iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ','LG')
      else
        ier = fstlir(lg,iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ','GL')
      end if

      ier = fstprm(key,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,          &
                   idum9,idum10,idum11,typxx,nomvxx,etikxx,grtyplg,ig1lg,ig2lg,  &
                   ig3lg,ig4lg,idum12,idum13,idum14,idum15,idum16,idum17,        &
                   idum18)

    end if ! firstCall

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

    call utl_reAllocate(mgintrp, npts)
    call utl_reAllocate(lgintrp, npts)

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
  ! mwbg_landIceMaskMwhs2
  !--------------------------------------------------------------------------
  subroutine mwbg_landIceMaskMwhs2(npts,zlat,zlon,zlq,ztt,waterobs,iopt)
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
    !           and LG are evaluated at the grid points of this mesh. For iopt=1 (or 3), The maximum
    !           (or average) determines whether the obs pt is too close to ice or land to be retained.
    !           For iopt=2, the value at the central mesh point is used.
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
    !--------------------------------------------------------------------
    !  Variable Definitions
    !  --------------------
    ! - npts       : input  -  number of input obs pts in report
    ! - zlat       : input  -  array holding lat values for all obs pts in report
    ! - zlon       : input  -  array holding lon values for all obs pts in report
    ! - zlq        : in/out -  array holding land/sea qualifier values for all obs
    !                        pts of report (0 = land, 1 = sea)
    ! - ztt        : in/out -  array holding terrain-type values for all obs pts
    !                        of current report (-1 land/open water, 0 = ice)
    ! - waterobs   : output -  logical array identifying for each obs in current report
    !                        whether it is over open water, far from coast/ice
    ! - iopt       : input  -  option for "interpolated" value of MG, LG at each location
    !                        1 = use MAX value taken from all mesh grid points
    !                        2 = use CENTRAL mesh point value (value at obs location)
    !                        3 = use AVG value of all mesh grid points
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
    integer, intent(in)                   :: npts
    integer, intent(in)                   :: iopt
    real,    intent(in)                   :: zlat(:)
    real,    intent(in)                   :: zlon(:)
    integer, intent(out), allocatable     :: zlq(:)
    integer, intent(out), allocatable     :: ztt(:)
    logical, intent(out), allocatable     :: waterobs(:)

    ! Locals:
    logical, save :: firstCall=.true.
    integer, parameter :: mxlat=5,mxlon=5
    integer :: iungeo

    integer :: ier,key
    integer, save :: ni,nj,nk,nilg,njlg
    integer, save :: ig1,ig2,ig3,ig4,ig1lg,ig2lg,ig3lg,ig4lg
    integer :: idum4,idum5,idum6,idum7,idum8,idum9,idum10,idum11
    integer :: idum12,idum13,idum14,idum15,idum16,idum17,idum18

    integer :: indx,ii,jj,kk
    integer :: nlat,nlon

    integer, parameter :: ii_obsloc=((mxlat*mxlon)/2)+1  ! 1D-index of central mesh-point (obs location)

    real, parameter :: pi=3.141592654
    real, parameter :: MGthresh=0.01,LGthresh=0.01
    real, parameter :: rlat_km=40.0,rlon_km=40.0
    real, parameter :: rkm_per_deg=111.195

    real :: xlat,xlatrad,xlon,rii,rjj
    real :: dlat,dlon

    character(len=12) :: etikxx
    character(len=4)  :: nomvxx
    character(len=2)  :: typxx
    character(len=1), save :: grtyp,grtyplg

    logical  :: llg

    ! F90 allocatable arrays:
    real, allocatable, save :: mg(:),lg(:)
    real, allocatable       :: latmesh(:),lonmesh(:)
    real, allocatable       :: mgintob(:),lgintob(:)
    real, allocatable       :: zlatbox(:,:),zlonbox(:,:)
    real, allocatable       :: mgintrp(:),lgintrp(:)

    ! RMNLIB interpolating functions:
    integer :: ezsetopt,ezqkdef
    integer :: gdllsval,gdid,gdidlg

    ! Define FORTRAN FST functions:
    integer, external :: fstinf,fstprm,fstlir,fnom,fclos
    integer, external :: fstouv,fstfrm,fstinl,fstvoi

    integer :: idum1,idum2,idum3

    ! Allocate space for arrays holding values on mesh grid pts.
    call utl_reAllocate(latmesh, mxlat*mxlon)
    call utl_reAllocate(lonmesh, mxlat*mxlon)
    call utl_reAllocate(mgintob, mxlat*mxlon)
    call utl_reAllocate(lgintob, mxlat*mxlon)
    call utl_reAllocate(zlatbox, mxlat*mxlon, npts)
    call utl_reAllocate(zlonbox, mxlat*mxlon, npts)
    call utl_reAllocate(zlq, npts)
    call utl_reAllocate(ztt, npts)
    call utl_reAllocate(waterobs, npts)

    ! Open FST file.
    iungeo = 0
    ier = fnom( iungeo,fileMgLg,'STD+RND+R/O',0 )
    ier = fstouv( iungeo,'RND' )

    if (firstCall) then
      firstCall = .false.

      ! Read MG field.
      key = fstinf(iungeo,ni,nj,nk,-1,' ',-1,-1,-1,' ' ,'MG')
      if ( key <  0 ) then
        call utl_abort('mwbg_landIceMaskMwhs2: The MG field is MISSING')
      end if

      call utl_reAllocate(mg, ni*nj)

      ier = fstlir(mg,iungeo,ni,nj,nk,-1,' ',-1,-1,-1,' ','MG')

      ier = fstprm(key,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,    &
                   idum9,idum10,idum11,typxx,nomvxx,etikxx,grtyp,ig1,ig2,  &
                   ig3,ig4,idum12,idum13,idum14,idum15,idum16,idum17,      &
                   idum18)

      ! Read LG field. Use GL field as backup.
      ! **CAUTION**: Discontinuities in GL field may cause interpolation problems! LG field is preferable.
      llg = .false.
      key = fstinf(iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ' ,'LG')
      if ( key <  0 ) then
        key = fstinf(iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ' ,'GL')
        if ( key <  0 ) then
          call utl_abort('mwbg_landIceMaskMwhs2: The LG or GL field is MISSING')
        end if
      else
        llg = .true.
      end if

      call utl_reAllocate(lg, nilg*njlg)

      if ( llg ) then
        ier = fstlir(lg,iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ','LG')
      else
        ier = fstlir(lg,iungeo,nilg,njlg,nk,-1,' ',-1,-1,-1,' ','GL')
      end if

      ier = fstprm(key,idum1,idum2,idum3,idum4,idum5,idum6,idum7,idum8,          &
                   idum9,idum10,idum11,typxx,nomvxx,etikxx,grtyplg,ig1lg,ig2lg,  &
                   ig3lg,ig4lg,idum12,idum13,idum14,idum15,idum16,idum17,        &
                   idum18)

    end if ! firstCall

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

    call utl_reAllocate(mgintrp, npts)
    call utl_reAllocate(lgintrp, npts)

    mgintrp(:) = 0.0
    lgintrp(:) = 0.0
    do kk = 1, npts

      latmesh = zlatbox(:,kk)
      lonmesh = zlonbox(:,kk)

      ier  = gdllsval(gdid,mgintob,mg,latmesh,lonmesh,mxlat*mxlon)
      ier  = gdllsval(gdidlg,lgintob,lg,latmesh,lonmesh,mxlat*mxlon)

      if ( iopt == 1 ) then
        mgintrp(kk) = maxval(mgintob(:))
        lgintrp(kk) = maxval(lgintob(:))
      elseif ( iopt == 2) then
        mgintrp(kk) = mgintob(ii_obsloc)
        lgintrp(kk) = lgintob(ii_obsloc)
      else
        mgintrp(kk) = sum(mgintob(:))/real((mxlat*mxlon))
        lgintrp(kk) = sum(lgintob(:))/real((mxlat*mxlon))
      endif

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

  end subroutine mwbg_landIceMaskMwhs2

  !--------------------------------------------------------------------------
  ! mwbg_computeMwhs2SurfaceType
  !--------------------------------------------------------------------------

  subroutine mwbg_computeMwhs2SurfaceType(obsSpaceData)
    ! :Purpose: Compute surface type element and update obsSpaceData.

    implicit none

    ! Arguments
    type(struct_obs), intent(inout) :: obsSpaceData           ! ObsSpaceData object

    integer, allocatable :: calcLandQualifierIndice(:)
    integer, allocatable :: calcTerrainTypeIndice(:)
    logical, allocatable :: waterobs(:)
    integer              :: codtyp
    integer              :: headerIndex

    logical              :: mwhs2DataPresent

    real                 :: obsLatitude(1)
    real                 :: obsLongitude(1)

    write(*,*) 'ssbg_computeMwhs2SurfaceType: Starting'

    mwhs2DataPresent = .false.
    call obs_set_current_header_list(obsSpaceData,'TO')

    HEADER0: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER0
      codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if ( tvs_isIdBurpInst(codtyp,'mwhs2') ) then
        mwhs2DataPresent = .true.
        exit HEADER0
      end if
    end do HEADER0

    if ( .not. mwhs2DataPresent ) then
      write(*,*) 'WARNING: WILL NOT RUN ssbg_computeMwhs2SurfaceType since no MWHS2 DATA is found'
      return
    end if

    call mwbg_init()

    if ( .not. modLSQ ) then
      write(*,*) 'WARNING: WILL NOT RUN ssbg_computeMwhs2SurfaceType since MODLSQ is not activated'
      return
    end if

    write(*,*) 'MWHS2 data found and modLSQ option activated!'
    write(*,*) '-->  Output file will contain recomputed values for land/sea qualifier and terrain type based on LG/MG.'

    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER1: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER1
      obsLatitude(1)  = obs_headElem_r( obsSpaceData, OBS_LAT, headerIndex )
      obsLatitude(1)  = obsLatitude(1) *MPC_DEGREES_PER_RADIAN_R8
      obsLongitude(1) = obs_headElem_r( obsSpaceData, OBS_LON, headerIndex )
      obsLongitude(1) = obsLongitude(1)*MPC_DEGREES_PER_RADIAN_R8
      if( obsLongitude(1) > 180. ) obsLongitude(1) = obsLongitude(1) - 360.
      call mwbg_landIceMaskMwhs2(1, obsLatitude, obsLongitude, calcLandQualifierIndice, &
                                 calcTerrainTypeIndice, waterobs, ilsmOpt)
      call obs_headSet_i(obsSpaceData, OBS_STYP, headerIndex, calcLandQualifierIndice(1))
      call obs_headSet_i(obsSpaceData, OBS_TTYP, headerIndex, calcTerrainTypeIndice(1))

    end do HEADER1

    write(*,*) 'ssbg_computeMwhs2SurfaceType: Finished'

  end subroutine mwbg_computeMwhs2SurfaceType

  !--------------------------------------------------------------------------
  ! mwbg_grossValueCheck  
  !--------------------------------------------------------------------------

  subroutine mwbg_grossValueCheck(npts, KNO, ztbcor, biasCorr, ztbThresholdMin, ztbThresholdMax, grossrej)

    !:Purpose: Check Tbs for values that are missing or outside physical limits.
    !          **NOTE: REJECT ALL CHANNELS OF ONE IS FOUND TO BE BAD.
    implicit none

    ! Arguments
    integer, intent(in)               :: npts             ! number of obs pts to process
    integer, intent(in)               :: KNO              ! number of ichannels
    real,    intent(in)               :: ztbcor(:)        ! bs from input BURP file
    real,    intent(in)               :: biasCorr(:)      ! Bias correction
    real,    intent(in)               :: ztbThresholdMin  ! ztb threshold for rejection
    real,    intent(in)               :: ztbThresholdMax  ! ztb threshold for rejection
    logical, intent(out), allocatable :: grossrej(:)      ! logical array defining which obs are to be rejected

    ! Locals
    integer :: pt, indx, channel

    real, allocatable                 :: ztb(:)           ! biased or unbiased radiances

    call utl_reAllocate(ztb, KNO)
    call utl_reAllocate(grossrej, npts)
    
    grossrej(1:npts) = .true.
    indx = 0
    do pt = 1, npts
      if ( mwbg_useUnbiasedObsForClw ) then
        do channel = 1, KNO
          ztb(channel) = ztbcor(indx+channel)
        end do
      else
        do channel = 1, KNO
          if (biasCorr(indx+channel) /= mwbg_realMissing) then
            ztb(channel) = ztbcor(indx+channel) - biasCorr(indx+channel)
          else
            ztb(channel) = ztbcor(indx+channel)
          end if
        end do
      end if
      if ( all( ztb > ztbThresholdMin ) .and. all( ztb < ztbThresholdMax ) ) then
        grossrej(pt) = .false.
      end if
      indx = pt*KNO

    end do

  end subroutine mwbg_grossValueCheck

  !--------------------------------------------------------------------------
  ! mwbg_firstQcCheckAtms
  !--------------------------------------------------------------------------
  subroutine mwbg_firstQcCheckAtms(zenith, ilq, itt, zlat, zlon, ztb, scanpos, &
                                   nval, nt, lqc, grossrej, lsq, trn, qcflag1, qcflag2, &
                                   ican, reportHasMissingTb)
    !  This routine performs basic quality control checks on the data. It sets array
    !  lqc(nt,nval) elements to .true. to flag data with failed checks.
    !
    !  The 6 QC checks are:
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
    !  In most cases, lqc(ii,nval) is set to .true. for all channels at point ii
    !  if the check detects a problem. In addition, Tb (ztb) is set to missing_value
    !  for checks 3 and 4 fails.
    implicit none

    ! Arguments
    integer,              intent(in)                :: ilq(:)
    integer,              intent(in)                :: itt(:)
    integer,              intent(in)                :: scanpos(:)
    integer,              intent(in)                :: ican(:)
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
    logical, allocatable, intent(out)               :: lqc(:,:)        ! dim(nt,nval), lqc = .false. on input

    ! Locals
    integer :: ii, jj, indx1, icount
    logical :: fail, fail1, fail2

    reportHasMissingTb = .false.
    call utl_reAllocate(lqc, nt, nval)
    lqc(:,:) = .false.  ! Flag for preliminary QC checks
    ! Global rejection checks

    ! Check if number of channels is correct
    !if ( nval /= mwbg_maxNumChan ) then
    !  write(*,*) 'WARNING: Number of channels (',nval, ') is not equal to mwbg_maxNumChan (', mwbg_maxNumChan,')'
    !  write(*,*) '         All data flagged as bad and returning to calling routine!'
    !  lqc(:,:) = .true.  ! flag all data in report as bad
    !  return
    !end if

    ! Check for errors in channel numbers (should be 1-22 for each location ii)
    indx1 = 1
    fail = .false.
    do ii = 1,nt
      do jj = 1, nval
        if ( ican(indx1+jj-1) /= jj ) fail = .true.
      end do
      indx1 = indx1 + nval
    end do
    if ( fail ) then
      write(*,*) 'WARNING: Bad channel number(s) detected!'
      write(*,*) '         All data flagged as bad and returning to calling routine!'
      write(*,*) '  ican(nt*nval) array = ', ican(:)
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
        zenith(ii) = mwbg_realMissing
        reportHasMissingTb = .true.
      end if
      do jj = 1,nval
        if ( fail ) then
          lqc(ii,jj) = .true.
          ztb(indx1+jj-1) = mwbg_realMissing
        end if
      end do
      indx1 = indx1 + nval
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
      do jj = 1, nval
        if ( fail ) then
          lqc(ii,jj) = .true.
          ztb(indx1+jj-1) = mwbg_realMissing
        end if
      end do
      indx1 = indx1 + nval
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
      do jj = 1, nval
        if ( fail ) then
          lqc(ii,jj) = .true.
          ztb(indx1+jj-1) = mwbg_realMissing
        end if
      end do
      indx1 = indx1 + nval
    end do
    if ( icount > 0 ) write(*,*) 'WARNING: Lat or lon out of range! Number of locations = ', icount

    !  5) Change in land/sea qualifier or terrain-type based on MG,LG fields
    icount = 0
    do ii = 1,nt
      fail = .false.
      if ( (ilq(ii) /= lsq(ii)) .or. (itt(ii) /= trn(ii)) ) then
        fail = .true.
      end if
      if ( fail ) then
        icount =  icount + 1
      end if
    end do

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
      do jj = 1, nval
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
      indx1 = indx1 + nval
    end do

    !write(*,*) 'mwbg_firstQcCheckAtms: Number of data processed and flagged = ', &
    !           nt*nval, count(lqc)

  end subroutine mwbg_firstQcCheckAtms

  !--------------------------------------------------------------------------
  ! mwbg_firstQcCheckMwhs2
  !--------------------------------------------------------------------------
  subroutine mwbg_firstQcCheckMwhs2(zenith, ilq, itt, zlat, zlon, ztb, scanpos, &
                                    nval, nt, lqc, lsq, trn, ican, reportHasMissingTb, modLSQ)
    !  This routine performs basic quality control checks on the data. It sets array
    !  lqc(nt,nval) elements to .true. to flag data with failed checks. Check 1
    !  (for ilq,itt) and check 5 are skipped if modlsqtt=.true., as the original values
    !  will be replaced in output file by lsq,trn.
    !
    !  The 5 QC checks are:
    !                 - 1) Invalid land/sea qualifier or terrain type,
    !                 - 2) Invalid field of view number,
    !                 - 3) Satellite zenith angle missing or out of range, (> 75 deg),
    !                 - 4) lat,lon check (lat,lon = O(-90.), 0(-180.))
    !                 - 5) Change in (computed) lsq,trn from (input) ilq,itt (from MG,LG fields)
    !                      ilq= 0,1 (from hi-res land/sea mask interpolated to obs point [CMDA])
    !                      itt=-1,0 (from hi-res ice analysis  interpolated to obs point [CMDA])
    !                      lsq= 0,1 (from max interp MG (0.0 to 1.0) in box surrounding obs point)
    !                      trn=-1,0 (from max interp LG (0.0 to 1.0) in box surrounding obs point)

    !
    !  In most cases, lqc(ii,nval) is set to .true. for all channels at point ii
    !  if the check detects a problem. In addition, Tb (ztb) is set to missing_value
    !  for checks 3 and 4 fails.
    implicit none

    ! Arguments
    integer,              intent(in)                :: ilq(:)
    integer,              intent(in)                :: itt(:)
    integer,              intent(in)                :: scanpos(:)
    integer,              intent(in)                :: ican(:)
    integer,              intent(in)                :: nt
    integer,              intent(in)                :: nval
    integer,              intent(in)                :: lsq(:)
    integer,              intent(in)                :: trn(:)
    real,                 intent(in)                :: zlat(:)
    real,                 intent(in)                :: zlon(:)
    real,                 intent(inout)             :: ztb(:)
    real,                 intent(inout)             :: zenith(:)
    logical,              intent(out)               :: reportHasMissingTb ! true if Tb(ztb) are set to missing_value
    logical,              intent(in)                :: modLSQ
    logical, allocatable, intent(out)               :: lqc(:,:)        ! dim(nt,nval), lqc = .false. on input

    ! Locals
    integer :: ii, jj, indx1, icount
    logical :: fail

    reportHasMissingTb = .false.
    call utl_reAllocate(lqc, nt, nval)
    lqc(:,:) = .false.  ! Flag for preliminary QC checks
    ! Global rejection checks

    ! Check if number of channels is correct
    !if ( nval /= mwbg_maxNumChan ) then
    !  write(*,*) 'WARNING: Number of channels (',nval, ') is not equal to mwbg_maxNumChan (', mwbg_maxNumChan,')'
    !  write(*,*) '         All data flagged as bad and returning to calling routine!'
    !  lqc(:,:) = .true.  ! flag all data in report as bad
    !  return
    !end if

    ! Check for errors in channel numbers (should be 1-15 for each location ii)
    indx1 = 1
    fail = .false.
    do ii = 1,nt
      do jj = 1, nval
        if ( ican(indx1+jj-1) /= jj ) fail = .true.
      end do
      indx1 = indx1 + nval
    end do
    if ( fail ) then
      write(*,*) 'WARNING: Bad channel number(s) detected!'
      write(*,*) '         All data flagged as bad and returning to calling routine!'
      write(*,*) '  ican(nt*nval) array = ', ican(:)
      lqc(:,:) = .true.  ! flag all data in report as bad
      return
    end if

    ! 1) invalid land/sea qualifier or terrain type
    !  ilq = 0 (land),     1 (sea)
    !  itt = 0 (sea-ice), -1 otherwise
    !  lsq = 1 (sea, away from land/coast [MG]),      0 otherwise
    !  trn = 0 (over or near analyzed sea-ice [LG]), -1 otherwise

    ! Checks on ilq,itt are not done if values are to be replaced in output file.

    if ( .not. modLSQ ) then
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
    end if

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
        zenith(ii) = mwbg_realMissing
        reportHasMissingTb = .true.
      end if
      do jj = 1,nval
        if ( fail ) then
          lqc(ii,jj) = .true.
          ztb(indx1+jj-1) = mwbg_realMissing
        end if
      end do
      indx1 = indx1 + nval
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
      do jj = 1, nval
        if ( fail ) then
          lqc(ii,jj) = .true.
          ztb(indx1+jj-1) = mwbg_realMissing
        end if
      end do
      indx1 = indx1 + nval
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
      do jj = 1, nval
        if ( fail ) then
          lqc(ii,jj) = .true.
          ztb(indx1+jj-1) = mwbg_realMissing
        end if
      end do
      indx1 = indx1 + nval
    end do
    if ( icount > 0 ) write(*,*) 'WARNING: Lat or lon out of range! Number of locations = ', icount

    !  5) Change in land/sea qualifier or terrain-type based on MG,LG fields
    if ( .not. modLSQ ) then
      icount = 0
      do ii = 1,nt
        fail = .false.
        if ( (ilq(ii) /= lsq(ii)) .or. (itt(ii) /= trn(ii)) ) then
          fail = .true.
        end if
        if ( fail ) then
          icount =  icount + 1
        end if
      end do
    end if

  end subroutine mwbg_firstQcCheckMwhs2

  !--------------------------------------------------------------------------
  ! mwbg_nrlFilterAtms
  !--------------------------------------------------------------------------
  subroutine mwbg_nrlFilterAtms(ni, KNO, ztbcor, zomp, biasCorr, pangl, plat, ilansea, iglace, waterobs, &
                                grossrej, clwObs, clwFG, si_ecmwf, si_bg, iNumSeaIce, iRej,SeaIce)
    !OBJET          Compute the following parameters using 5 ATMS channels:
    !                  - sea ice,
    !                  - cloud liquid water from observation (clwObs),
    !                  - cloud liquid water from first guess (clwFG),
    !                  - 2 scattering indices (si) (ECMWF, Bennartz-Grody)
    !               The five channels used are: 23Ghz, 31Ghz, 50Ghz, 89Ghz, and 165Ghz.
    !
    !NOTES*
    !                - open water points are converted to sea-ice points if sea ice concentration >= 0.55
    !                   and iglace (itt or terrain type) is changed accordingly
    !                - clwObs are missing when out-of-range parameters/Tb detected or grossrej = .true.
    !                - clwObs and si only computed over open water away from coasts and sea-ice
    !                - clwObs and si = -99.0 where value cannot be computed.
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
    !               - tb23FG      - input  -  23Ghz brightness temperature (K) from first guess -- ch. 1
    !               - tb31        - input  -  31Ghz brightness temperature (K) -- ch. 2
    !               - tb31FG      - input  -  31Ghz brightness temperature (K) from first guess -- ch. 2
    !               - tb50        - input  -  50Ghz brightness temperature (K) -- ch. 3
    !               - tb89        - input  -  89Ghz brightness temperature (K) -- ch. 16
    !               - tb165       - input  -  165Ghz brightness temperature (K) -- ch. 17
    !               - pangl       - input  -  satellite zenith angle (deg.)
    !               - plat        - input  -  latitude (deg.)
    !               - ilansea     - input  -  land/sea indicator (0=land, 1=ocean)
    !               - iglace      - in/out -  terrain type (0=ice, -1 otherwise)
    !               - waterobs    - in/out -  .true. if open water point (away from coasts and sea-ice)
    !               - grossrej    - input  -  .true. if any channel had a gross error from mwbg_grossValueCheck
    !               - clwObs      - output -  cloud liquid water from observation (kg/m**2) from tb23 & tb31
    !               - clwFG       - output -  cloud liquid water from first guess (kg/m**2) from tb23FG & tb31FG
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
    integer, intent(in)                   ::  KNO
    integer, intent(out)                  ::  iNumSeaIce
    integer, intent(in)                   ::  ilansea(:)
    integer, intent(inout)                ::  iglace(:)
    integer, intent(out)                  ::  iRej


    logical, intent(in)                   ::  grossrej(:)
    logical, intent(inout)                ::  waterobs(:)

    real, intent(in)                      ::  ztbcor(:)
    real, intent(in)                      ::  zomp(:)
    real, intent(in)                      ::  biasCorr(:)
    real, intent(in)                      ::  pangl(:)
    real, intent(in)                      ::  plat(:)
    real, allocatable, intent(out)        ::  clwObs(:)
    real, allocatable, intent(out)        ::  clwFG(:)
    real, allocatable, intent(out)        ::  si_ecmwf(:)
    real, allocatable, intent(out)        ::  si_bg(:)
    real, allocatable, intent(out)        ::  SeaIce(:)

    ! Locals
    integer                               :: ier(ni)
    real                                  ::  ice(ni)
    real                                  :: tb23(ni)
    real                                  :: tb23FG(ni)
    real                                  :: tb31(ni)
    real                                  :: tb31FG(ni)
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
    real                                  :: t23FG
    real                                  :: t31
    real                                  :: t31FG
    real                                  :: t50
    real                                  :: t89
    real                                  :: t165


    ! Allocation
    call utl_reAllocate(clwObs,ni)
    call utl_reAllocate(clwFG,ni)
    call utl_reAllocate(si_ecmwf,ni)
    call utl_reAllocate(si_bg,ni)
    call utl_reAllocate(SeaIce,ni)

    iNumSeaIce = 0

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
      indx2 = ii*KNO
      tb23(ii)      = ztbcor(indx1)
      tb23FG(ii)    = ztbcor(indx1) - zomp(indx1)
      bcor23(ii)    = biasCorr(indx1)
      tb31(ii)      = ztbcor(indx1+1)
      tb31FG(ii)    = ztbcor(indx1+1) - zomp(indx1+1)
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

    ! 1) Initialise parameters:
    do i = 1, ni
      ice(i)      = mwbg_realMissing
      clwObs(i)   = mwbg_realMissing
      clwFG(i)    = mwbg_realMissing
      si_ecmwf(i) = mwbg_realMissing
      si_bg(i)    = mwbg_realMissing
      SeaIce(i)   = 0.0
    end do

    ! 2) Validate input parameters:
    do i = 1, ni
      if ( pangl(i)   <   0.  .or. &
           pangl(i)   >  70.  .or. &
           plat(i)    < -90.  .or. &
           plat(i)    >  90.  .or. &
           ilansea(i) <   0   .or. &
           ilansea(i) >   1        ) then
         ier(i) = 1
      end if

      ! Skip computations for points where all data are rejected  (bad Tb ANY channel)
      if ( grossrej(i) ) then
        ier(i) = 1
      end if
    end do

    ! 3) Compute parameters:
    do i = 1, ni

      if ( ier(i) == 0 ) then

        abslat = abs(plat(i))
        cosz   = cosd(pangl(i))

        if (bcor23(i) == mwbg_realMissing .or. mwbg_useUnbiasedObsForClw) then
          t23 = tb23(i)
        else
          t23 = tb23(i) - bcor23(i)
        end if
        if (bcor31(i) == mwbg_realMissing .or. mwbg_useUnbiasedObsForClw) then
          t31 = tb31(i)
        else
          t31 = tb31(i) - bcor31(i)
        end if
        if (bcor50(i) == mwbg_realMissing .or. mwbg_useUnbiasedObsForClw) then
          t50 = tb50(i)
        else
          t50 = tb50(i) - bcor50(i)
        end if
        if (bcor89(i) == mwbg_realMissing .or. mwbg_useUnbiasedObsForClw) then
          t89 = tb89(i)
        else
          t89 = tb89(i) - bcor89(i)
        end if
        if (bcor165(i) == mwbg_realMissing .or. mwbg_useUnbiasedObsForClw) then
          t165 = tb165(i)
        else
          t165 = tb165(i) - bcor165(i)
        end if

        deltb = t89 - t165
        t23FG = tb23FG(i)
        t31FG = tb31FG(i)

        ! Check for sea-ice over water points. Set terrain type to 0 if ice>=0.55 detected.
        if ( ilansea(i) == 1 ) then  ! water point

          if ( abslat < 50. ) then
            ice(i) = 0.0
          else
            ice(i) = 2.85 + 0.020*t23 - 0.028*t50
          end if

          SeaIce(i) = ice(i)

          if ( ice(i) >= 0.55 .and. waterobs(i) ) then
            iNumSeaIce = iNumSeaIce + 1
            waterobs(i) = .false.
            iglace(i) = 0
          end if

        end if

        ! Compute clwObs, clwFG, and Scattering Indices (over open water only)
        if ( waterobs(i) ) then
          if ( t23 < 284. .and. t31 < 284. ) then
            aa = 8.24 - (2.622 - 1.846 * cosz) * cosz
            clwObs(i) = aa + 0.754 * alog(285.0 - t23) - 2.265 * alog(285.0 - t31)
            clwObs(i) = clwObs(i) * cosz
            if ( clwObs(i) < 0.0 ) clwObs(i) = 0.0

            clwFG(i) = aa + 0.754 * alog(285.0 - t23FG) - 2.265 * alog(285.0 - t31FG)
            clwFG(i) = clwFG(i) * cosz
            if ( clwFG(i) < 0.0 ) clwFG(i) = 0.0
          end if
          si_ecmwf(i) = deltb - (-46.94 + 0.248 * pangl(i))
          si_bg(i)    = deltb - (-39.201 + 0.1104 * pangl(i))
        end if

      else  ! ier(i) .eq. 1 case
         iRej = iRej + 1

      end if ! if ( ier(i) .eq. 0 )

      if ( mwbg_debug .and. (i <= 100) ) then
        write(*,*) ' '
        write(*,*) ' i,tb23(i),tb23FG(i),tb31(i),tb31FG(i),tb50(i),tb89(i),tb165(i),pangl(i),plat(i), ilansea(i) = ', &
     &             i,tb23(i),tb23FG(i),tb31(i),tb31FG(i),tb50(i),tb89(i),tb165(i),pangl(i),plat(i), ilansea(i)
        write(*,*) ' ier(i),ice(i),clwObs(i),clwFG(i),si_ecmwf(i),si_bg(i),iglace(i),waterobs(i) =',ier(i),ice(i),&
     &             clwObs(i),clwFG(i),si_ecmwf(i),si_bg(i),iglace(i),waterobs(i)
      end if

    end do   ! i loop over ni points

  end subroutine mwbg_nrlFilterAtms

  !--------------------------------------------------------------------------
  ! mwbg_nrlFilterMwhs2
  !--------------------------------------------------------------------------
  subroutine mwbg_nrlFilterMwhs2(ni, KNO, ztbcor, biasCorr, pangl, plat, ilansea, iglace, waterobs, &
                                 grossrej, clwObs, clwFG, si_ecmwf, si_bg, iNumSeaIce, iRej,SeaIce)
    !OBJET          Compute the following parameters using 2 MWHS2 channels:
    !                  - sea ice,
    !                  - cloud liquid water from observation (clwObs),
    !                  - cloud liquid water from first guess (clwFG),
    !                  - 2 scattering indices (si) (ECMWF, Bennartz-Grody)
    !               The two channels used are: 89Ghz, and 165Ghz.
    !
    !NOTES*
    !                - open water points are converted to sea-ice points if sea ice concentration >= 0.55
    !                   and iglace (itt or terrain type) is changed accordingly
    !                - clwObs are missing when out-of-range parameters/Tb detected or grossrej = .true.
    !                - clwObs and si_ecmwf only computed over open water away from coasts and sea-ice
    !                - si_bg is computed for all points
    !                - clwObs and si = -99.0 where value cannot be computed.
    !
    !REFERENCES     Ben Ruston, NRL Monterey
    !                  JCSDA Seminar 12/12/12: Impact of NPP Satellite Assimilation in the U.S. Navy Global Modeling System
    !
    !
    !ARGUMENTS      - ier         - output - error return code for each location:
    !                                        0, ok,
    !                                        1, input parameter out of range or grossrej=.true.
    !               - ni          - input  -  number of points to process (= NT)
    !               - tb23        - input  -  23Ghz brightness temperature (K) -- ch. 1 [missing for MWHS-2]
    !               - tb23FG      - input  -  23Ghz brightness temperature (K) from first guess -- ch. 1
    !               - tb31        - input  -  31Ghz brightness temperature (K) -- ch. 2 [missing for MWHS-2]
    !               - tb31FG      - input  -  31Ghz brightness temperature (K) from first guess -- ch. 2
    !               - tb50        - input  -  50Ghz brightness temperature (K) -- ch. 3 [missing for MWHS-2]
    !               - tb89        - input  -  89Ghz brightness temperature (K) -- ch. 16
    !               - tb165       - input  -  165Ghz brightness temperature (K) -- ch. 17
    !               - pangl       - input  -  satellite zenith angle (deg.)
    !               - plat        - input  -  latitude (deg.)
    !               - ilansea     - input  -  land/sea indicator (0=land, 1=ocean)
    !               - iglace      - in/out -  terrain type (0=ice, -1 otherwise)
    !               - waterobs    - in/out -  .true. if open water point (away from coasts and sea-ice)
    !               - grossrej    - input  -  .true. if any channel had a gross error from mwbg_grossValueCheck
    !               - clwObs      - output -  cloud liquid water from observation (kg/m**2) from tb23 & tb31
    !               - clwFG       - output -  cloud liquid water from first guess (kg/m**2) from tb23FG & tb31FG
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
    integer, intent(in)                   ::  KNO
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
    real, allocatable, intent(out)        ::  clwObs(:)
    real, allocatable, intent(out)        ::  clwFG(:)
    real, allocatable, intent(out)        ::  si_ecmwf(:)
    real, allocatable, intent(out)        ::  si_bg(:)
    real, allocatable, intent(out)        ::  SeaIce(:)

    ! Locals
    integer                               :: ier(ni)
    real                                  ::  ice(ni)
    real                                  :: tb23(ni)
    real                                  :: tb23FG(ni)
    real                                  :: tb31(ni)
    real                                  :: tb31FG(ni)
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
    real                                  :: t23FG
    real                                  :: t31
    real                                  :: t31FG
    real                                  :: t50
    real                                  :: t89
    real                                  :: t165


    ! Allocation
    call utl_reAllocate(clwObs,ni)
    call utl_reAllocate(clwFG,ni)
    call utl_reAllocate(si_ecmwf,ni)
    call utl_reAllocate(si_bg,ni)
    call utl_reAllocate(SeaIce,ni)

    iNumSeaIce = 0

    !! extract required channels:  ATMS ch.   MWHS-2 ch.
    !!   23 Ghz = AMSU-A 1 =        1          n/a
    !!   31 Ghz = AMSU-A 2 =        2          n/a
    !!   50 Ghz = AMSU-A 3 =        3          n/a
    !!   53 Ghz = AMSU-A 5 =        6          n/a
    !!   89 Ghz = AMSU-A 15/ B 1 = 16           1
    !!  150 Ghz = AMSU-B 2 =       17          10
    !!  183 Ghz = AMSU-B 3 =       22          11
    !!  183 Ghz = AMSU-B 5 =       18          15
    !
    !   Extract Tb for channels 1 (AMSU-B 1) and 10 (AMSU-B 2) for Bennartz SI
    !   Extract Tb for channels 22 (AMSU-B 3) and 18 (AMSU-B 5) for Dryness Index (DI)

    indx1 = 1
    do ii = 1, ni
      indx2 = ii*KNO
      tb23(ii)      = mwbg_realMissing
      tb23FG(ii)    = mwbg_realMissing
      bcor23(ii)    = mwbg_realMissing
      tb31(ii)      = mwbg_realMissing
      tb31FG(ii)    = mwbg_realMissing
      bcor31(ii)    = mwbg_realMissing
      tb50(ii)      = mwbg_realMissing
      bcor50(ii)    = mwbg_realMissing
      tb89(ii)      = ztbcor(indx1)
      bcor89(ii)    = biasCorr(indx1)
      tb165(ii)    = ztbcor(indx1+9)
      bcor165(ii)    = biasCorr(indx1+9)
      indx1 = indx2 + 1
    end do

    ier = 0

    ! 1) Initialise parameters:
    do i = 1, ni
      ice(i)      = mwbg_realMissing
      clwObs(i)   = mwbg_realMissing
      clwFG(i)    = mwbg_realMissing
      si_ecmwf(i) = mwbg_realMissing
      si_bg(i)    = mwbg_realMissing
      SeaIce(i)   = 0.0
    end do

    ! 2) Validate input parameters:
    do i = 1, ni
      if ( pangl(i)   <   0.  .or. &
           pangl(i)   >  70.  .or. &
           plat(i)    < -90.  .or. &
           plat(i)    >  90.  .or. &
           ilansea(i) <   0   .or. &
           ilansea(i) >   1        ) then
         ier(i) = 1
      end if

      ! Skip computations for points where all data are rejected  (bad Tb ANY channel)
      if ( grossrej(i) ) then
        ier(i) = 1
      end if
    end do

    ! 3) Compute parameters:
    do i = 1, ni

      if ( ier(i) == 0 ) then

        abslat = abs(plat(i))
        cosz   = cosd(pangl(i))

        if (bcor23(i) == mwbg_realMissing .or. mwbg_useUnbiasedObsForClw) then
          t23 = tb23(i)
        else
          t23 = tb23(i) - bcor23(i)
        end if
        if (bcor31(i) == mwbg_realMissing .or. mwbg_useUnbiasedObsForClw) then
          t31 = tb31(i)
        else
          t31 = tb31(i) - bcor31(i)
        end if
        if (bcor50(i) == mwbg_realMissing .or. mwbg_useUnbiasedObsForClw) then
          t50 = tb50(i)
        else
          t50 = tb50(i) - bcor50(i)
        end if
        if (bcor89(i) == mwbg_realMissing .or. mwbg_useUnbiasedObsForClw) then
          t89 = tb89(i)
        else
          t89 = tb89(i) - bcor89(i)
        end if
        if (bcor165(i) == mwbg_realMissing .or. mwbg_useUnbiasedObsForClw) then
          t165 = tb165(i)
        else
          t165 = tb165(i) - bcor165(i)
        end if

        deltb = t89 - t165
        t23FG = tb23FG(i)
        t31FG = tb31FG(i)

        ! Check for sea-ice over water points. Set terrain type to 0 if ice>=0.55 detected.
        if ( ilansea(i) == 1 .and. t23 .ne. mwbg_realMissing ) then  ! water point

          if ( abslat < 50. ) then
            ice(i) = 0.0
          else
            ice(i) = 2.85 + 0.020*t23 - 0.028*t50
          end if

          SeaIce(i) = ice(i)
          if ( ice(i) >= 0.55 .and. waterobs(i) ) then
            iNumSeaIce = iNumSeaIce + 1
            waterobs(i) = .false.
            iglace(i) = 0
          end if

        end if

        ! Compute clwObs, clwFG, and Scattering Indices (over open water only)
        if ( waterobs(i) ) then
          if ( t23 .ne. mwbg_realMissing ) then
            if ( t23 < 284. .and. t31 < 284. ) then
              aa = 8.24 - (2.622 - 1.846 * cosz) * cosz
              clwObs(i) = aa + 0.754 * alog(285.0 - t23) - 2.265 * alog(285.0 - t31)
              clwObs(i) = clwObs(i) * cosz
              if ( clwObs(i) < 0.0 ) clwObs(i) = 0.0

              clwFG(i) = aa + 0.754 * alog(285.0 - t23FG) - 2.265 * alog(285.0 - t31FG)
              clwFG(i) = clwFG(i) * cosz
              if ( clwFG(i) < 0.0 ) clwFG(i) = 0.0
            end if
          end if
          si_ecmwf(i) = deltb - (-46.94 + 0.248 * pangl(i))
          si_bg(i)    = deltb - (-39.201 + 0.1104 * pangl(i))
        else
          si_bg(i)    = deltb - (0.158 + 0.0163 * pangl(i))
        end if

      else  ! ier(i) .eq. 1 case
         iRej = iRej + 1

      end if ! if ( ier(i) .eq. 0 )

      if ( mwbg_debug .and. (i <= 100) ) then
        write(*,*) ' '
        write(*,*) ' i,tb23(i),tb23FG(i),tb31(i),tb31FG(i),tb50(i),tb89(i),tb165(i),pangl(i),plat(i), ilansea(i) = ', &
     &             i,tb23(i),tb23FG(i),tb31(i),tb31FG(i),tb50(i),tb89(i),tb165(i),pangl(i),plat(i), ilansea(i)
        write(*,*) ' ier(i),ice(i),clwObs(i),clwFG(i),si_ecmwf(i),si_bg(i),iglace(i),waterobs(i) =',ier(i),ice(i),&
     &             clwObs(i),clwFG(i),si_ecmwf(i),si_bg(i),iglace(i),waterobs(i)
      end if

    end do   ! i loop over ni points

  end subroutine mwbg_nrlFilterMwhs2

  !--------------------------------------------------------------------------
  ! mwbg_flagDataUsingNrlCritAtms
  !--------------------------------------------------------------------------
  subroutine mwbg_flagDataUsingNrlCritAtms(nt, nval, ztbcor, biasCorr, clwObs, scatec, scatbg, SeaIce, grossrej, waterobs, &
                                           useUnbiasedObsForClw, iwvreject, cloudobs, precipobs,  cldcnt, ident, riwv, zdi)

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
    integer, intent(in)                        :: nval
    real, intent(in)                           :: ztbcor(:)
    real, intent(in)                           :: biasCorr(:)
    real, intent(in)                           :: clwObs (:)
    real, intent(in)                           :: scatec(:)
    real, intent(in)                           :: scatbg(:)
    real, intent(in)                           :: SeaIce (:)

    logical, intent(in)                        :: useUnbiasedObsForClw
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
    integer                                    :: indx
    integer                                    :: indx1
    integer                                    :: indx2
    integer                                    :: ii
    integer                                    :: n_cld
    real, allocatable                          :: ztb_amsub3(:)
    real, allocatable                          :: bcor_amsub3(:)
    real, allocatable                          :: ztb_amsub5(:)
    real, allocatable                          :: bcor_amsub5(:)
    real                                       ::  ztb183(5)


    call utl_reAllocate(cloudobs, nt)
    call utl_reAllocate(iwvreject, nt)
    call utl_reAllocate(ident, nt)
    call utl_reAllocate(precipobs, nt)
    call utl_reAllocate(riwv, nt)
    call utl_reAllocate(ztb_amsub3, nt)
    call utl_reAllocate(bcor_amsub3, nt)
    call utl_reAllocate(ztb_amsub5, nt)
    call utl_reAllocate(bcor_amsub5, nt)

    ! To begin, assume that all obs are good.
    ident(:) = 0
    cloudobs(:)  = .false.
    iwvreject(:) = .false.
    precipobs(:) = .false.

    ! Extract Tb for channels 16 (AMSU-B 1) and 17 (AMSU-B 2) for Bennartz SI
    ! Extract Tb for channels 22 (AMSU-B 3) and 18 (AMSU-B 5) for Dryness Index (DI)

    indx1 = 1
    do ii = 1, nt
      indx2 = ii*nval
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
      indx2 = ii*nval
      if (.not.grossrej(ii)) then
        do indx = 1, 5
          if (biasCorr(indx1+indx+9) == mwbg_realMissing .or. useUnbiasedObsForClw) then
            ztb183(indx) = ztbcor(indx1+indx+16)
          else
            ztb183(indx) = ztbcor(indx1+indx+16) - biasCorr(indx1+indx+16)
          end if
        end do
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
    where ( scatec > scatec_atms_nrl_LTrej .or. scatbg > scatbg_atms_nrl_LTrej ) precipobs = .true.
    n_cld = count(clwObs > clw_atms_nrl_LTrej)
    cldcnt  = cldcnt  + n_cld
    where ( (clwObs > clw_atms_nrl_LTrej) .or. precipobs ) cloudobs = .true.
    where ( waterobs )  ident = IBSET(ident,0)
    where ( iwvreject ) ident = IBSET(ident,5)
    where ( precipobs ) ident = IBSET(ident,4)
    where ( clwObs > clw_atms_nrl_LTrej) ident = IBSET(ident,3)
    where ( clwObs > clw_atms_nrl_UTrej) ident = IBSET(ident,6)
    where ( scatec > scatec_atms_nrl_UTrej .or. scatbg > scatbg_atms_nrl_UTrej ) ident = IBSET(ident,8)
    where ( SeaIce >= 0.55 ) ident = IBSET(ident,10)

    where ( waterobs .and. (clwObs == -99.) ) ident = IBSET(ident,2)
    where ( riwv == -99.)                   ident = IBSET(ident,1)

    ! Compute the simple AMSU-B Dryness Index zdi for all points = Tb(ch.3)-Tb(ch.5)
    if ( useUnbiasedObsForClw ) then
      where ( .not.grossrej )
        zdi = ztb_amsub3 - ztb_amsub5
      elsewhere
        zdi = mwbg_realMissing
      end where
    else
      where ( .not.grossrej )
        zdi = (ztb_amsub3 - bcor_amsub3) - (ztb_amsub5 - bcor_amsub5)
      elsewhere
        zdi = mwbg_realMissing
      end where
    end if

  end subroutine mwbg_flagDataUsingNrlCritAtms

  !--------------------------------------------------------------------------
  ! mwbg_flagDataUsingNrlCritMwhs2
  !--------------------------------------------------------------------------
  subroutine mwbg_flagDataUsingNrlCritMwhs2(nt, nval, ztbcor, biasCorr, clwObs, scatec, SeaIce, grossrej, waterobs, &
                                            useUnbiasedObsForClw, iwvreject, cloudobs, precipobs,  cldcnt, ident, riwv, zdi)

    !:Purpose:                       Set the  Information flag (ident) values (new BURP element 025174 in header)
    !                                BIT    Meaning
    !                                - 0     off=land or sea-ice, on=open water away from coast
    !                                - 1     Mean 183 Ghz [ch. 18-22] is missing
    !                                - 2     CLW is missing (over water)
    !                                - 3     CLW > clw_mwhs2_nrl_LTrej (0.175 kg/m2) (cloudobs)
    !                                - 4     scatec > Lower Troposphere limit 9/10 (precipobs)
    !                                - 5     Mean 183 Ghz [ch. 18-22] Tb < 240K
    !                                - 6     CLW > clw_mwhs2_nrl_UTrej (0.200 kg/m2)
    !                               - 10     Sea ice > 0.55 detected
    !                               - 11     Gross error in Tb (any chan.)  (all channels rejected)

    ! Arguments
    integer, intent(in)                        :: nt
    integer, intent(in)                        :: nval
    real, intent(in)                           :: ztbcor(:)
    real, intent(in)                           :: biasCorr(:)
    real, intent(in)                           :: clwObs (:)
    real, intent(in)                           :: scatec(:)
    real, intent(in)                           :: SeaIce (:)

    logical, intent(in)                        :: useUnbiasedObsForClw
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
    integer                                    :: indx
    integer                                    :: indx1
    integer                                    :: indx2
    integer                                    :: ii
    integer                                    :: n_cld
    real, allocatable                          :: ztb_amsub3(:)
    real, allocatable                          :: bcor_amsub3(:)
    real, allocatable                          :: ztb_amsub5(:)
    real, allocatable                          :: bcor_amsub5(:)
    real                                       ::  ztb183(5)


    call utl_reAllocate(cloudobs, nt)
    call utl_reAllocate(iwvreject, nt)
    call utl_reAllocate(ident, nt)
    call utl_reAllocate(precipobs, nt)
    call utl_reAllocate(riwv, nt)
    call utl_reAllocate(ztb_amsub3, nt)
    call utl_reAllocate(bcor_amsub3, nt)
    call utl_reAllocate(ztb_amsub5, nt)
    call utl_reAllocate(bcor_amsub5, nt)

    ! To begin, assume that all obs are good.
    ident(:) = 0
    cloudobs(:)  = .false.
    iwvreject(:) = .false.
    precipobs(:) = .false.

    ! Extract Tb for channels 1 (AMSU-B 1) and 10 (AMSU-B 2) for Bennartz SI
    ! Extract Tb for channels 11 (AMSU-B 3) and 15 (AMSU-B 5) for Dryness Index (DI)

    indx1 = 1
    do ii = 1, nt
      indx2 = ii*nval
      ztb_amsub3(ii) = ztbcor(indx1+10)
      bcor_amsub3(ii) = biasCorr(indx1+10)
      ztb_amsub5(ii) = ztbcor(indx1+14)
      bcor_amsub5(ii) = biasCorr(indx1+14)
      indx1 = indx2 + 1
    end do


    ! Flag data using NRL criteria

    ! Compute Mean 183 Ghz [ch. 11-15] Tb (riwv)
    riwv = -99.0
    indx1 = 1
    do ii = 1, nt
      indx2 = ii*nval
      if (.not.grossrej(ii)) then
        do indx = 1, 5
          if (biasCorr(indx1+indx+9) == mwbg_realMissing .or. useUnbiasedObsForClw) then
            ztb183(indx) = ztbcor(indx1+indx+9)
          else
            ztb183(indx) = ztbcor(indx1+indx+9) - biasCorr(indx1+indx+9)
          end if
        end do
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
    where ( scatec > scatec_mwhs2_nrl_LTrej ) precipobs = .true.
    n_cld = count(clwObs > clw_mwhs2_nrl_LTrej)
    cldcnt  = cldcnt  + n_cld
    where ( (clwObs > clw_mwhs2_nrl_LTrej) .or. precipobs ) cloudobs = .true.
    where ( waterobs )  ident = IBSET(ident,0)
    where ( iwvreject ) ident = IBSET(ident,5)
    where ( precipobs ) ident = IBSET(ident,4)
    where ( clwObs > clw_mwhs2_nrl_LTrej) ident = IBSET(ident,3)
    where ( clwObs > clw_mwhs2_nrl_UTrej) ident = IBSET(ident,6)
    where ( SeaIce >= 0.55 ) ident = IBSET(ident,10)

    where ( waterobs .and. (clwObs == -99.) ) ident = IBSET(ident,2)
    where ( riwv == -99.)                   ident = IBSET(ident,1)

    ! Compute the simple AMSU-B Dryness Index zdi for all points = Tb(ch.3)-Tb(ch.5)
    if ( useUnbiasedObsForClw ) then
      where ( .not. grossrej )
        zdi = ztb_amsub3 - ztb_amsub5
      elsewhere
        zdi = mwbg_realMissing
      end where
    else
      where ( .not. grossrej )
        zdi = (ztb_amsub3 - bcor_amsub3) - (ztb_amsub5 - bcor_amsub5)
      elsewhere
        zdi = mwbg_realMissing
      end where
    end if

  end subroutine mwbg_flagDataUsingNrlCritMwhs2

  !--------------------------------------------------------------------------
  ! mwbg_reviewAllCritforFinalFlagsAtms
  !--------------------------------------------------------------------------
  subroutine mwbg_reviewAllCritforFinalFlagsAtms(nt, nval, lqc, grossrej, waterobs, &
                                                 precipobs, clwObs, clwFG, scatec, scatbg, &
                                                 iwvreject, riwv, IMARQ, globMarq, zdi, ident, &
                                                 drycnt, landcnt, rejcnt, iwvcnt, pcpcnt, flgcnt, &
                                                 MXCLWREJ, chanFlaggedForAllskyGenCoeff, icano)

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
    real, intent(inout)                        :: clwObs(:)
    real, intent(inout)                        :: clwFG(:)
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
    integer, intent(inout)                     :: drycnt
    integer, intent(inout)                     :: landcnt
    integer, intent(inout)                     :: rejcnt
    integer, intent(inout)                     :: iwvcnt
    integer, intent(inout)                     :: pcpcnt
    integer, intent(inout)                     :: flgcnt
    integer, intent(in)                        :: MXCLWREJ
    integer, intent(in)                        :: chanFlaggedForAllskyGenCoeff(:)
    integer, intent(in)                        :: icano(:)

    ! Locals
    real                                       :: clwObsFGaveraged
    logical, allocatable                       :: lflagchn(:,:)
    integer                                    :: kk, j, ipos, INDXCAN


    ! Allocation
    call utl_reAllocate(lflagchn,nt, nval)

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

        else  ! if waterobs(kk)

        ! OVER WATER,
        !    in clear-sky mode:
        !    -- reject ch. 5-6, if CLW > clw_atms_nrl_LTrej or CLW = -99.0
        !    in all-sky mode:
        !    -- reject ch. 5-6, if CLW > mwbg_clwQcThreshold or CLW = -99.0
        !
        !    -- reject ch. 1-4, if CLW > clw_atms_nrl_LTrej or CLW = -99.0
        !    -- reject ch. 16-20 if CLW > clw_atms_nrl_LTrej or CLW = -99.0
        !    -- reject ch. 7-9, 21-22 if CLW > clw_atms_nrl_UTrej or CLW = -99.0
        !    -- reject ch. 1-6, 16-22 if scatec > 9  or scatec = -99.0
        !    -- reject ch. 7-9        if scatec > 18 or scatec = -99.0
        !    -- reject ch. 1-6        if scatbg > 10 or scatbg = -99.0
        !    -- reject ch. 7-9        if scatbg > 15 or scatbg = -99.0
        !    -- reject ch. 16-22      if iwvreject = .true.   [ Mean 183 Ghz [ch. 18-22] Tb < 240K ]

          if ( clwObs(kk)   >  clw_atms_nrl_LTrej )  then
            if ( tvs_mwAllskyAssim ) then
              lflagchn(kk,1:4) = .true.
              clwObsFGaveraged = 0.5 * (clwObs(kk) + clwFG(kk))
              if ( clwObsFGaveraged > mwbg_clwQcThreshold ) lflagchn(kk,5:6) = .true.
            else
              lflagchn(kk,1:mwbg_atmsNumSfcSensitiveChannel) = .true.
            end if
            lflagchn(kk,16:20) = .true.
          end if
          if ( clwObs(kk)   >  clw_atms_nrl_UTrej )  then
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
          if ( clwObs(kk) == -99. ) then
            ident(kk) = IBSET(ident(kk),2)
            lflagchn(kk,1:9)   = .true.
            lflagchn(kk,16:22) = .true.
          end if
          if ( riwv(kk) == -99. ) then     ! riwv = mean_Tb_183Ghz
            ident(kk) = IBSET(ident(kk),1)
            lflagchn(kk,16:22) = .true.
          end if
        end if  ! if waterobs(kk)

      end if  ! if .not. grossrej(kk)

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
    ! Set missing clwObs and riwv to BURP missing value (mwbg_realMissing)
    where (clwObs == -99. ) clwObs = mwbg_realMissing
    where (clwObs == -99. ) clwFG = mwbg_realMissing
    where (riwv == -99. ) riwv = mwbg_realMissing

    ! Modify data flag values (set bit 7) for rejected data
    ! In all-sky mode, turn on bit=23 for cloud-affected radiances when
    ! there is mismatch between clwObs and clwFG (to be used in gen_bias_corr)
    ipos=0
    do kk =1, nt
      clwObsFGaveraged = 0.5 * (clwObs(kk) + clwFG(kk))

      do j = 1, nval
        ipos = ipos + 1
        if (lflagchn(kk,j)) then
          IMARQ(ipos) = IBSET(IMARQ(ipos),7)
        end if

        INDXCAN = ISRCHEQI(chanFlaggedForAllskyGenCoeff, MXCLWREJ, ICANO(ipos))
        if ( tvs_mwAllskyAssim .and. waterobs(kk) .and. INDXCAN /= 0 .and. &
             (clwObsFGaveraged > mwbg_cloudyClwThresholdBcorr .or. &
              clwObs(kk) == mwbg_realMissing) ) then
          IMARQ(ipos) = IBSET(IMARQ(ipos),23)
        end if
      end do
    end do


    ! Set bit 6 in 24-bit global flags if any data rejected
    do kk =1, nt
      if ( ANY(lflagchn(kk,:)) ) globMarq(kk) = IBSET(globMarq(kk),6)
    end do

  end subroutine mwbg_reviewAllCritforFinalFlagsAtms

  !--------------------------------------------------------------------------
  ! mwbg_reviewAllCritforFinalFlagsMwhs2
  !--------------------------------------------------------------------------
  subroutine mwbg_reviewAllCritforFinalFlagsMwhs2(nt, nval, lqc, grossrej, trn, waterobs, &
                                                  precipobs, clwObs, clwFG, scatec, scatbg, &
                                                  iwvreject, riwv, IMARQ, globMarq, zdi, ident, &
                                                  allcnt, drycnt, landcnt, rejcnt, iwvcnt, pcpcnt, flgcnt, &
                                                  MXCLWREJ, chanFlaggedForAllskyGenCoeff, icano)

    !:Purpose:                   Review all the checks previously made to determine which obs are to be accepted
    !                            for assimilation and which are to be flagged for exclusion (lflagchn).
    !                            - grossrej()  = .true. if any channel had a gross error at the point
    !                            - cloudobs()  = .true. if CLW > clw_mwhs2_nrl_LTrej (0.175) or precipobs
    !                            - precipobs() = .true. if precip. detected through NRL scattering indices
    !                            - waterobs()  = .true. if open water point
    !                            - iwvreject() = .true. if Mean 183 Ghz [ch. 18-22] Tb < 240K (too dry for ch.20-22 over land)
    ! Arguments
    integer, intent(in)                        :: nt
    integer, intent(in)                        :: nval
    logical, intent(in)                        :: lqc(:,:)
    real, intent(inout)                        :: clwObs(:)
    real, intent(inout)                        :: clwFG(:)
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
    integer, intent(inout)                     :: allcnt
    integer, intent(inout)                     :: drycnt
    integer, intent(inout)                     :: landcnt
    integer, intent(inout)                     :: rejcnt
    integer, intent(inout)                     :: iwvcnt
    integer, intent(inout)                     :: pcpcnt
    integer, intent(inout)                     :: flgcnt
    integer, intent(inout)                     :: trn(:)
    integer, intent(in)                        :: MXCLWREJ
    integer, intent(in)                        :: chanFlaggedForAllskyGenCoeff(:)
    integer, intent(in)                        :: icano(:)

    ! Locals
    real                                       :: clwObsFGaveraged
    real                                       :: scatbg_rej
    logical, allocatable                       :: lflagchn(:,:)
    integer                                    :: kk, j, ipos, INDXCAN


    ! Allocation
    call utl_reAllocate(lflagchn,nt, nval)

    lflagchn(:,:) = lqc(:,:)  ! initialize with flags set in mwbg_firstQcCheckMwhs2
    do kk = 1, nt
      allcnt = allcnt + 1  ! Counting total number of observations
      ! Reject all channels if gross Tb error detected in any channel or other problems
      if ( grossrej(kk) ) then
        lflagchn(kk,:) = .true.
      else

        ! OVER LAND OR SEA-ICE,
        !    -- CLW/SI not determined over land
        !    -- surface emissivity effects lower tropospheric and window channels
        !    -- reject window & lower tropospheric channels 1,10,14,15
        !    -- reject ch. 11-13 if iwvreject = .true.  [ Mean 183 Ghz [ch. 18-22] Tb < 240K ]
        !    -- check DI for AMSU-B like channels
        !    -- reject all channels if scatbg exceeds CMC thresholds for AMSU-B

        if  ( .not. waterobs(kk) ) then
          lflagchn(kk,(/ 1,10,14,15 /)) = .true.      ! AMSU-B (like 1,2,5)
          if ( iwvreject(kk) ) lflagchn(kk,11:13) = .true.  ! AMSU-B (like 4,3)

          ! Dryness index (for AMSU-B channels 11-14 assimilated over land/sea-ice)
          ! Channel AMSUB-3 (MWHS-2 channel 11) is rejected for a dryness index >    0.
          !                 (MWHS-2 channel 12) is rejected for a dryness index >   -5.
          ! Channel AMSUB-4 (MWHS-2 channel 13) is rejected for a dryness index >   -8.
          if ( zdi(kk) > 0.0 ) then
            lflagchn(kk,11) = .true.
            ident(kk) = IBSET(ident(kk),7)
          end if
          if ( zdi(kk) > -5.0 ) then
            lflagchn(kk,12) = .true.
            ident(kk) = IBSET(ident(kk),9)
            drycnt = drycnt + 1
          end if
          if ( zdi(kk) > -8.0 ) then
            lflagchn(kk,13) = .true.
          end if

          ! Bennartz -Grody SI check thresholds (same as for QC of AMSU-B/MHS)
          if ( trn(kk) == 0 ) then ! sea-ice
            scatbg_rej = scatbg_mwhs2_cmc_ICErej
          else                     ! land
            scatbg_rej = scatbg_mwhs2_cmc_LANDrej
          end if
          if ( scatbg(kk) > scatbg_rej ) then
            lflagchn(kk,:) = .true.
            ident(kk) = IBSET(ident(kk),8)
          end if

        else  ! if waterobs(kk)

        ! OVER WATER,
        !-----------------------------------------------------------------
        ! PLACEHOLDER VALUES FOR ALLSKY ASSIM, SINCE NOT IMPLEMENTED YET
        !    in clear-sky mode:
        !    -- reject ch. 1, if CLW > clw_mwhs2_nrl_LTrej or CLW = -99.0
        !    in all-sky mode:
        !    -- reject ch. 1, if CLW > mwbg_clwQcThreshold or CLW = -99.0
        !-----------------------------------------------------------------
        !    -- reject ch. 1, 10, 13-15 if CLW > clw_mwhs2_nrl_LTrej
        !    -- reject ch. 11-12 if CLW > clw_mwhs2_nrl_UTrej
        !    -- reject ch. 1, 10-15 if scatec > 9  or scatec = -99.0
        !    -- reject ch. 1, 10-15 if iwvreject = .true.   [ Mean 183 Ghz [ch. 11-15] Tb < 240K ]
        !    -- reject all channels if scatbg exceeds CMC SEA threshold for AMSU-B

          if ( clwObs(kk)   >  clw_mwhs2_nrl_LTrej )  then
            if ( tvs_mwAllskyAssim ) then ! NEVER TRUE SINCE NOT IMPLEMENTED YET
              clwObsFGaveraged = 0.5 * (clwObs(kk) + clwFG(kk))
              if ( clwObsFGaveraged > mwbg_clwQcThreshold ) lflagchn(kk,1) = .true.
            else
              lflagchn(kk,1) = .true.
            end if
            lflagchn(kk,(/ 10,13,14,15 /)) = .true.
          end if
          if ( clwObs(kk)   >  clw_mwhs2_nrl_UTrej )  then
            lflagchn(kk,11:12) = .true.
          end if
          if ( scatec(kk) >  scatec_mwhs2_nrl_LTrej ) then
            lflagchn(kk,1) = .true.
            lflagchn(kk,10:15) = .true.
          end if
          if ( iwvreject(kk) ) then
            lflagchn(kk,1) = .true.
            lflagchn(kk,10:15) = .true.
          end if
          if ( riwv(kk) == -99. ) then     ! riwv = mean_Tb_183Ghz
            ident(kk) = IBSET(ident(kk),1)
            lflagchn(kk,1) = .true.
            lflagchn(kk,10:15) = .true.
          end if
          ! Bennartz-Grody SI check thresholds (same as for QC of AMSU-B/MHS)
          if ( scatbg(kk) > scatbg_mwhs2_cmc_SEA ) then
            lflagchn(kk,:) = .true.
            ident(kk) = IBSET(ident(kk),8)
          end if
        end if  ! if waterobs(kk)

      end if  ! if .not. grossrej(kk)

      if ( .not. waterobs(kk) ) landcnt  = landcnt  + 1
      if ( grossrej(kk) )  rejcnt = rejcnt + 1
      if ( iwvreject(kk))  iwvcnt = iwvcnt + 1
      if ( precipobs(kk) .and. waterobs(kk) ) then
        pcpcnt = pcpcnt + 1
      end if

      if ( ANY(lflagchn(kk,:)) ) flgcnt = flgcnt + 1
    end do

    ! RESET riwv array to Bennartz-Grody scattering index for output to BURP file
    riwv(:) = scatbg(:)
    ! Set missing clwObs and riwv to BURP missing value (mwbg_realMissing)
    where (clwObs == -99. ) clwObs = mwbg_realMissing
    where (clwObs == -99. ) clwFG = mwbg_realMissing
    where (riwv == -99. ) riwv = mwbg_realMissing

    ! Modify data flag values (set bit 7) for rejected data
    ! In all-sky mode, turn on bit=23 for cloud-affected radiances when
    ! there is mismatch between clwObs and clwFG (to be used in gen_bias_corr)
    ipos=0
    do kk =1, nt
      clwObsFGaveraged = 0.5 * (clwObs(kk) + clwFG(kk))

      do j = 1, nval
        ipos = ipos + 1
        if (lflagchn(kk,j)) then
          IMARQ(ipos) = IBSET(IMARQ(ipos),7)
        end if

        INDXCAN = ISRCHEQI(chanFlaggedForAllskyGenCoeff, MXCLWREJ, ICANO(ipos))
        if ( tvs_mwAllskyAssim .and. waterobs(kk) .and. INDXCAN /= 0 .and. &
             (clwObsFGaveraged > mwbg_cloudyClwThresholdBcorr .or. &
              clwObs(kk) == mwbg_realMissing) ) then
          IMARQ(ipos) = IBSET(IMARQ(ipos),23)
        end if
      end do
    end do


    ! Set bit 6 in 24-bit global flags if any data rejected
    do kk =1, nt
      if ( ANY(lflagchn(kk,:)) ) globMarq(kk) = IBSET(globMarq(kk),6)
    end do

  end subroutine mwbg_reviewAllCritforFinalFlagsMwhs2

  function calcStateDepObsErr_r4(cloudPredictorThresh1, cloudPredictorThresh2, &
                                 sigmaThresh1, sigmaThresh2, cloudPredictorUsed) result(sigmaObsErrUsed)
    !
    ! :Purpose: Calculate single-precision state-dependent observation error.
    !                                 
    implicit none

    ! Arguments:
    real, intent(in) :: cloudPredictorThresh1
    real, intent(in) :: cloudPredictorThresh2
    real, intent(in) :: sigmaThresh1
    real, intent(in) :: sigmaThresh2
    real, intent(in) :: cloudPredictorUsed
    real :: sigmaObsErrUsed

    if (cloudPredictorUsed <= cloudPredictorThresh1) then
      sigmaObsErrUsed = sigmaThresh1
    else if (cloudPredictorUsed >  cloudPredictorThresh1 .and. & 
             cloudPredictorUsed <= cloudPredictorThresh2) then
      sigmaObsErrUsed = sigmaThresh1 + &
                        (sigmaThresh2 - sigmaThresh1) / &
                        (cloudPredictorThresh2 - cloudPredictorThresh1) * &
                        (cloudPredictorUsed - cloudPredictorThresh1) 
    else
      sigmaObsErrUsed = sigmaThresh2
    end if

  end function calcStateDepObsErr_r4

  !--------------------------------------------------------------------------
  ! mwbg_updateObsSpaceAfterQc
  !--------------------------------------------------------------------------
  subroutine mwbg_updateObsSpaceAfterQc(obsSpaceData, sensorIndex, headerIndex, obsTb, obsFlags, &
                                        cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                                        atmScatteringIndexObs, atmScatteringIndexFG, &
                                        obsGlobalMarker, newInformationFlag)

    !:Purpose:      Update obspacedata variables (obstTB and obs flags) after QC
    implicit None

    !Arguments
    type(struct_obs),     intent(inout)     :: obsSpaceData           ! obspaceData Object
    integer,              intent(in)        :: sensorIndex            ! tvs_sensorIndex
    integer,              intent(in)        :: headerIndex            ! current header index
    integer,              intent(in)        :: obsFlags(:)            ! data flags
    real,                 intent(in)        :: obsTb(:)               ! obs Tb
    real,                 intent(in)        :: cloudLiquidWaterPathObs(:)   ! obs CLW
    real,                 intent(in)        :: cloudLiquidWaterPathFG(:)    ! trial CLW
    real,                 intent(in)        :: atmScatteringIndexObs(:)  ! atmospheric scatering index from observation
    real,                 intent(in)        :: atmScatteringIndexFG(:)   ! atmospheric scatering index from background
    integer,              intent(in)        :: newInformationFlag(:)     ! information flag used with satplot
    integer,              intent(in)        :: obsGlobalMarker(:)        ! information flag used with satplot
    ! Locals
    integer                                 :: bodyIndex
    integer                                 :: obsNumCurrentLoc
    integer                                 :: bodyIndexbeg
    integer                                 :: currentChannelNumber
    integer                                 :: codtyp

    codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
    call obs_headSet_r(obsSpaceData, OBS_CLWO, headerIndex, cloudLiquidWaterPathObs(1))

    if (tvs_isInstrumAllskyTtAssim(tvs_getInstrumentId(codtyp_get_name(codtyp))) .or. &
        tvs_isInstrumAllskyTtHuAssim(tvs_getInstrumentId(codtyp_get_name(codtyp)))) then
      call obs_headSet_r(obsSpaceData, OBS_CLWB, headerIndex, cloudLiquidWaterPathFG(1))
    end if

    if (atmScatteringIndexObs(1) /= mwbg_realMissing) then
      call obs_headSet_r(obsSpaceData, OBS_SIO, headerIndex, atmScatteringIndexObs(1))
    else
      call obs_headSet_r(obsSpaceData, OBS_SIO, headerIndex, MPC_missingValue_R4)
    end if

    if (tvs_isInstrumAllskyHuAssim(tvs_getInstrumentId(codtyp_get_name(codtyp))) .or. &
        tvs_isInstrumAllskyTtHuAssim(tvs_getInstrumentId(codtyp_get_name(codtyp)))) then
      if (atmScatteringIndexFG(1) /= mwbg_realMissing) then
        call obs_headSet_r(obsSpaceData, OBS_SIB, headerIndex, atmScatteringIndexFG(1))
      else
        call obs_headSet_r(obsSpaceData, OBS_SIB, headerIndex, MPC_missingValue_R4)
      end if
    end if
    
    call obs_headSet_i(obsSpaceData, OBS_INFG, headerIndex, newInformationFlag(1))
    call obs_headSet_i(obsSpaceData, OBS_ST1, headerIndex, obsGlobalMarker(1))
    bodyIndexbeg        = obs_headElem_i( obsSpaceData, OBS_RLN, headerIndex )
    obsNumCurrentLoc    = obs_headElem_i( obsSpaceData, OBS_NLV, headerIndex )
    BODY: do bodyIndex =  bodyIndexbeg, bodyIndexbeg + obsNumCurrentLoc - 1
      currentChannelNumber=nint(obs_bodyElem_r( obsSpaceData,  OBS_PPP, bodyIndex ))-tvs_channelOffset(sensorIndex)
      call obs_bodySet_r(obsSpaceData, OBS_VAR,   bodyIndex, obsTb(currentChannelNumber))
      call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex, obsFlags(currentChannelNumber))
    end do BODY

  end subroutine mwbg_updateObsSpaceAfterQc  
   
  !--------------------------------------------------------------------------
  !  mwbg_readObsFromObsSpace
  !--------------------------------------------------------------------------
  subroutine mwbg_readObsFromObsSpace(instName, headerIndex, &
                                      satIdentifier, satZenithAngle, landQualifierIndice, &
                                      terrainTypeIndice, obsLatitude, obsLongitude, &
                                      satScanPosition, obsQcFlag1, satOrbit, & 
                                      obsGlobalMarker, burpFileSatId, obsTb, btClear, &
                                      obsTbBiasCorr, ompTb, obsQcFlag2, obsChannels, &
                                      obsFlags, sensorIndex, actualNumChannel, obsSpaceData)
    
    !:Purpose:        copy headers and bodies from obsSpaceData object to arrays

    implicit None

    !Arguments
    character(len=9),     intent(in)     :: InstName               ! Instrument Name
    integer,              intent(in)     :: headerIndex            ! current header Index 
    integer, allocatable, intent(out)    :: satIdentifier(:)       ! satellite identifier
    real   , allocatable, intent(out)    :: satZenithAngle(:)      ! satellite zenith angle (btyp=3072,ele=7024) 
    integer, allocatable, intent(out)    :: landQualifierIndice(:) ! land/sea qualifier     (btyp=3072,ele=8012)
    integer, allocatable, intent(out)    :: terrainTypeIndice(:)   ! terrain-type (ice)     (btyp=3072,ele=13039)
    real   , allocatable, intent(out)    :: obsLatitude(:)         ! latitude values (btyp=5120,ele=5002)
    real   , allocatable, intent(out)    :: obsLongitude(:)        ! longitude values (btyp=5120,ele=6002)
    integer, allocatable, intent(out)    :: satScanPosition(:)     ! scan position (fov)    (btyp=3072,ele=5043)
    integer, allocatable, intent(out)    :: obsQcFlag1(:,:)        ! flag values for btyp=3072 block ele 033078, 033079, 033080
    integer, allocatable, intent(out)    :: satOrbit(:)            ! orbit number
    integer, allocatable, intent(out)    :: obsGlobalMarker(:)     ! global Marqueur Data
    character(*),intent(out)             :: burpFileSatId          ! Platform Name
    real   , allocatable, intent(out)    :: obsTb(:)               ! brightness temperature (btyp=9248/9264,ele=12163) 
    real   , allocatable, intent(out)    :: btClear(:)             ! clear brightness temperature (btyp=9248/9264,ele=btClearElementId)
    real   , allocatable, intent(out)    :: obsTbBiasCorr(:)       ! bias correction 
    real   , allocatable, intent(out)    :: ompTb(:)               ! OMP values
    integer, allocatable, intent(out)    :: obsQcFlag2(:)          ! flag values for btyp=9248 block ele 033081      
    integer, allocatable, intent(out)    :: obsChannels(:)         ! channel numbers btyp=9248 block ele 5042 (= 1-22)
    integer, allocatable, intent(out)    :: obsFlags(:)            ! data flags
    integer,              intent(out)    :: sensorIndex            ! find tvs_sensor index corresponding to current obs
    integer,              intent(out)    :: actualNumChannel       ! actual Num channel

    type(struct_obs),     intent(inout)  :: obsSpaceData           ! obspaceData Object

    ! Locals
    integer                              :: bodyIndex
    integer                              :: obsNumCurrentLoc
    integer                              :: bodyIndexbeg
    integer                              :: headerCompt 
    integer                              :: currentChannelNumber  
    integer                              :: channelIndex
    integer                              :: numObsToProcess       
    integer                              :: iplatform
    integer                              :: instrum
    integer                              :: isat, iplatf
    integer                              :: instr
    integer                              :: codtyp
    logical                              :: sensorIndexFound

   ! find tvs_sensor index corresponding to current obs

    iplatf      = obs_headElem_i( obsSpaceData, OBS_SAT, headerIndex )
    instr       = obs_headElem_i( obsSpaceData, OBS_INS, headerIndex )

    call tvs_mapSat( iplatf, iplatform, isat )
    call tvs_mapInstrum( instr, instrum )
    
    sensorIndexFound = .false.
    do sensorIndex =1, tvs_nsensors
      if ( iplatform ==  tvs_platforms(sensorIndex)  .and. &
           isat      ==  tvs_satellites(sensorIndex) .and. &
           instrum   == tvs_instruments(sensorIndex)       ) then
          sensorIndexFound = .true. 
         exit
      end if
    end do
    if ( .not. sensorIndexFound ) call utl_abort('mwbg_readObsFromObsSpace: sensor Index not found') 

    ! find actual Number of channels
    actualNumChannel = tvs_coefs(sensorIndex)%coef%fmv_ori_nchn

    headerCompt = 1 
    numObsToProcess = 1
    ! Allocate Header elements
    call utl_reAllocate(satIdentifier, numObsToProcess)
    call utl_reAllocate(satZenithAngle, numObsToProcess)
    call utl_reAllocate(landQualifierIndice, numObsToProcess)
    call utl_reAllocate(terrainTypeIndice, numObsToProcess)
    call utl_reAllocate(obsLatitude, numObsToProcess)
    call utl_reAllocate(obsLongitude, numObsToProcess)
    call utl_reAllocate(satScanPosition, numObsToProcess)
    call utl_reAllocate(obsGlobalMarker, numObsToProcess)
    call utl_reAllocate(satOrbit, numObsToProcess)
    call utl_reAllocate(obsQcFlag1, numObsToProcess,3)
    ! Allocate Body elements
    call utl_reAllocate(obsTb, numObsToProcess*actualNumChannel)
    call utl_reAllocate(btClear, numObsToProcess*actualNumChannel)
    call utl_reAllocate(ompTb, numObsToProcess*actualNumChannel)
    call utl_reAllocate(obsTbBiasCorr, numObsToProcess*actualNumChannel)
    call utl_reAllocate(obsFlags, numObsToProcess*actualNumChannel)
    call utl_reAllocate(obsChannels, numObsToProcess*actualNumChannel)
    call utl_reAllocate(obsQcFlag2, numObsToProcess*actualNumChannel)
    !initialization
    obsTb(:) = mwbg_realMissing
    btClear(:) = mwbg_realMissing
    ompTb(:) = mwbg_realMissing
    obsTbBiasCorr(:) = mwbg_realMissing

        
    burpFileSatId                      = obs_elem_c    (obsSpaceData, 'STID' , headerIndex) 
    satIdentifier(headerCompt)         = obs_headElem_i(obsSpaceData, OBS_SAT, headerIndex) 
    satZenithAngle(headerCompt)        = obs_headElem_r(obsSpaceData, OBS_SZA, headerIndex) 
    landQualifierIndice(headerCompt)   = obs_headElem_i(obsSpaceData, OBS_STYP, headerIndex) 
    terrainTypeIndice(headerCompt)     = obs_headElem_i(obsSpaceData, OBS_TTYP, headerIndex) 
    ! If terrain type is missing, set it to -1 for the QC programs
    if (terrainTypeIndice(headerCompt) ==  99) terrainTypeIndice(headerCompt) = -1
    obsLatitude (headerCompt)          = obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex) 
    obsLongitude(headerCompt)          = obs_headElem_r(obsSpaceData, OBS_LON, headerIndex) 
    ! Convert lat/lon to degrees
    obsLongitude(headerCompt) = obsLongitude(headerCompt)*MPC_DEGREES_PER_RADIAN_R8
    if (obsLongitude(headerCompt) > 180.) obsLongitude(headerCompt) = obsLongitude(headerCompt) - 360.
    obsLatitude(headerCompt)  = obsLatitude(headerCompt) *MPC_DEGREES_PER_RADIAN_R8
    satScanPosition(headerCompt)       = obs_headElem_i(obsSpaceData, OBS_FOV , headerIndex) 
    obsGlobalMarker(headerCompt)       = obs_headElem_i(obsSpaceData, OBS_ST1, headerIndex) 
    satOrbit(headerCompt)              = obs_headElem_i(obsSpaceData, OBS_ORBI, headerIndex) 
    if (instName == 'ATMS') then  
      obsQcFlag1(headerCompt,1)        = obs_headElem_i(obsSpaceData, OBS_AQF1, headerIndex) 
      obsQcFlag1(headerCompt,2)        = obs_headElem_i(obsSpaceData, OBS_AQF2, headerIndex) 
      obsQcFlag1(headerCompt,3)        = obs_headElem_i(obsSpaceData, OBS_AQF3, headerIndex) 
    end if

    bodyIndexbeg        = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
    obsNumCurrentLoc    = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex)
    codtyp              = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)

    BODY: do bodyIndex =  bodyIndexbeg, bodyIndexbeg + obsNumCurrentLoc - 1
      currentChannelNumber = nint(obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex)) - &
                             tvs_channelOffset(sensorIndex)
      obsTb(currentChannelNumber) = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex)
      if (tvs_isInstrumAllskyTtAssim(tvs_getInstrumentId(codtyp_get_name(codtyp))) .or. &
          tvs_isInstrumAllskyHuAssim(tvs_getInstrumentId(codtyp_get_name(codtyp))) .or. &
          tvs_isInstrumAllskyTtHuAssim(tvs_getInstrumentId(codtyp_get_name(codtyp)))) then
        btClear(currentChannelNumber) = obs_bodyElem_r(obsSpaceData,  OBS_BTCL, bodyIndex)
      end if
      ompTb(currentChannelNumber)          = obs_bodyElem_r(obsSpaceData, OBS_OMP, bodyIndex)
      obsTbBiasCorr(currentChannelNumber)  = obs_bodyElem_r(obsSpaceData, OBS_BCOR,bodyIndex)
      obsFlags(currentChannelNumber)       = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex)
      obsQcFlag2(currentChannelNumber)     = obs_bodyElem_i(obsSpaceData, OBS_QCF2, bodyIndex)
      
    end do BODY
    do channelIndex=1, actualNumChannel
      obsChannels(channelIndex)    = channelIndex+tvs_channelOffset(sensorIndex)
    end do
      

  end subroutine mwbg_readObsFromObsSpace 

  !--------------------------------------------------------------------------
  ! mwbg_mwbg_bgCheckMW
  !--------------------------------------------------------------------------
  subroutine mwbg_bgCheckMW( obsSpaceData )

    !:Purpose:        Do the quality control for ATMS, AMSUA, AMSUB and MWHS2

    implicit None

    !Arguments
    type(struct_obs),     intent(inout)  :: obsSpaceData           ! obspaceData Object

    ! Locals
    integer                       :: numObsToProcess               ! number of obs in current report
    integer                       :: headerIndex                   ! header Index
    integer                       :: sensorIndex                   ! satellite index in obserror file
    integer                       :: actualNumChannel              ! iactual Number of channel for instrument
    integer                       :: codtyp                        ! codetype
    character(len=9)              :: burpFileSatId                 ! station id in burp file
    real, allocatable             :: modelInterpTerrain(:)         ! topo in standard file interpolated to obs point
    real, allocatable             :: modelInterpSeaIce(:)          ! Glace de mer " "
    real, allocatable             :: modelInterpGroundIce(:)       ! Glace de continent " "
    real,    allocatable          :: obsLatitude(:)                ! obs. point latitudes
    real,    allocatable          :: obsLongitude(:)               ! obs. point longitude
    integer, allocatable          :: satIdentifier(:)              ! Satellite identifier
    real,    allocatable          :: satZenithAngle(:)             ! sat. satZenithAngle angle
    integer, allocatable          :: landQualifierIndice(:)        ! land qualifyer
    integer, allocatable          :: terrainTypeIndice(:)          ! terrain type
    real,    allocatable          :: obsTb(:)                      ! temperature de brillance
    real,    allocatable          :: btClear(:)                    ! clear-sky BT
    real,    allocatable          :: ompTb(:)                      ! o-p temperature de "
    real,    allocatable          :: obsTbBiasCorr(:)              ! bias correction fo obsTb
    integer, allocatable          :: satScanPosition(:)            ! scan position
    integer, allocatable          :: obsQcFlag1(:,:)               ! Obs Quality flag 1
    integer, allocatable          :: obsQcFlag2(:)                 ! Obs Quality flag 2 
    integer, allocatable          :: obsChannels(:)                ! obsTb channels
    integer, allocatable          :: obsFlags(:)                   ! obs. flag
    integer, allocatable          :: satOrbit(:)                   ! orbit
    integer, allocatable          :: obsGlobalMarker(:)            ! global marker
    integer, allocatable          :: qcIndicator(:,:)              ! indicateur controle de qualite tovs par canal 
    !                                                                =0, ok,
    !                                                                >0, rejet,
    integer, allocatable          :: newInformationFlag(:)         ! ATMS Information flag (newInformationFlag) values 
    !                                                                (new BURP element  025174 in header). FOR AMSUA 
    !
    real,    allocatable          :: cloudLiquidWaterPathObs(:)    ! cloud liquid water path from observation.
    real,    allocatable          :: cloudLiquidWaterPathFG(:)     ! cloud liquid water path from background.
    real,    allocatable          :: atmScatteringIndexObs(:)      ! scattering index from observation.
    real,    allocatable          :: atmScatteringIndexFG(:)       ! scattering index from background.
    integer, external             :: exdb, exfin, fnom, fclos
    logical                       :: mwDataPresent
    logical                       :: lastHeader                    ! active while reading last report


    call utl_tmg_start(118,'--BgckMicrowave')
    mwDataPresent = .false.
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER0: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
    if (headerIndex < 0) exit HEADER0
      codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if ( tvs_isIdBurpInst(codtyp,'atms' ) .or. &
           tvs_isIdBurpInst(codtyp,'amsua') .or. &
           tvs_isIdBurpInst(codtyp,'amsub') .or. & 
           tvs_isIdBurpInst(codtyp,'mwhs2') .or. &
           tvs_isIdBurpInst(codtyp,'mhs'  ) ) then
        mwDataPresent = .true.
      end if
    end do HEADER0

    if ( .not. mwDataPresent ) then 
      write(*,*) 'WARNING: WILL NOT RUN mwbg_bgCheckMW since no ATMS or AMSUA or MWHS2'
      return
    end if 

    write(*,*) ' MWBG QC PROGRAM STARTS ....'
    ! read nambgck
    call mwbg_init()

    !Quality Control loop over all observations
    !
    ! loop over all header indices of the specified family with surface obs
    numObsToProcess = 1
    lastHeader = .false.

    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER
      if (headerIndex == obs_numHeader(obsSpaceData)) lastHeader = .true.
      codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if (instName== 'AMSUB') then
        if ( .not. tvs_isIdBurpInst(codtyp,'amsub') .and. &
             .not. tvs_isIdBurpInst(codtyp,'mhs'  )  .and. &
             .not. tvs_isIdBurpInst(codtyp,'mwhs2') ) then
          write(*,*) 'WARNING: Observation with codtyp = ', codtyp, ' is not ', instName
          cycle HEADER
        end if 
      else
        if ( .not. (tvs_isIdBurpInst(codtyp,instName)) ) then
          write(*,*) 'WARNING: Observation with codtyp = ', codtyp, ' is not ', instName
          cycle HEADER
        end if
      end if
      !###############################################################################
      ! STEP 1) read obs from obsSpacedata to start QC                               !
      !###############################################################################

      call mwbg_readObsFromObsSpace(instName, headerIndex,                            &
                                   satIdentifier, satZenithAngle,landQualifierIndice, &
                                   terrainTypeIndice, obsLatitude, obsLongitude,      &
                                   satScanPosition, obsQcFlag1, satOrbit,             &
                                   obsGlobalMarker, burpFileSatId, obsTb, btClear,    &
                                   obsTbBiasCorr, ompTb, obsQcFlag2, obsChannels,     &
                                   obsFlags, sensorIndex, actualNumChannel, obsSpaceData)

      !###############################################################################
      ! STEP 3) Interpolation de le champ MX(topogrpahy), MG et GL aux pts TOVS.
      !###############################################################################
      call mwbg_readGeophysicFieldsAndInterpolate(instName, obsLatitude, &
                                                  obsLongitude, modelInterpTerrain,     &
                                                  modelInterpGroundIce, modelInterpSeaIce)
      !###############################################################################
      ! STEP 4) Controle de qualite des TOVS. Data QC flags (obsFlags) are modified here!
      !###############################################################################

      if (instName == 'AMSUA') then
        call mwbg_tovCheckAmsua(oer_toverrst, oer_cloudPredictorThreshArr, oer_sigmaObsErr, oer_useStateDepSigmaObs, &
                                oer_tovutil, landQualifierIndice,&
                                obsChannels, obsTb, btClear, obsTbBiasCorr, &
                                ompTb, qcIndicator, actualNumChannel, numObsToProcess, sensorIndex, &
                                satScanPosition, modelInterpGroundIce, modelInterpTerrain,&
                                modelInterpSeaIce, terrainTypeIndice, satZenithAngle,     &
                                obsGlobalMarker, obsFlags, newInformationFlag, &
                                cloudLiquidWaterPathObs, cloudLiquidWaterPathFG,      &
                                atmScatteringIndexObs, burpFileSatId, RESETQC, obsLatitude)
      else if (instName == 'AMSUB') then
        call mwbg_tovCheckAmsub(oer_toverrst, oer_cloudPredictorThreshArr, oer_sigmaObsErr, oer_useStateDepSigmaObs, &
                                oer_tovutil, landQualifierIndice,&
                                obsChannels, obsTb, obsTbBiasCorr, ompTb,      & 
                                qcIndicator, actualNumChannel, numObsToProcess, sensorIndex, &
                                satScanPosition, modelInterpGroundIce, modelInterpTerrain,&
                                modelInterpSeaIce, terrainTypeIndice, satZenithAngle,     &
                                obsGlobalMarker, obsFlags, newInformationFlag,        & 
                                cloudLiquidWaterPathObs, cloudLiquidWaterPathFG,       &
                                atmScatteringIndexObs, atmScatteringIndexFG, burpFileSatId, RESETQC)
      else if (instName == 'ATMS') then
        call mwbg_tovCheckAtms(oer_toverrst, oer_cloudPredictorThreshArr, oer_sigmaObsErr, oer_useStateDepSigmaObs, &
                               oer_tovutil, obsLatitude, obsLongitude,&
                               landQualifierIndice, terrainTypeIndice, satZenithAngle,   &
                               obsQcFlag2, obsQcFlag1, &
                               obsChannels, obsTb, obsTbBiasCorr, ompTb, qcIndicator,   &
                               actualNumChannel, numObsToProcess, sensorIndex,          &
                               newInformationFlag, satScanPosition,   &
                               modelInterpTerrain, obsGlobalMarker, obsFlags,            &
                               cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                               atmScatteringIndexObs, burpFileSatId, RESETQC)
      else if (instName == 'MWHS2') then
        call mwbg_tovCheckMwhs2(oer_toverrst, oer_clwThreshArr, oer_sigmaObsErr, oer_useStateDepSigmaObs, &
                                oer_tovutil, obsLatitude, obsLongitude, &
                                landQualifierIndice, terrainTypeIndice, satZenithAngle,  &
                                obsChannels, obsTb, obsTbBiasCorr, ompTb, qcIndicator,   &
                                actualNumChannel, numObsToProcess, sensorIndex,          &
                                newInformationFlag, satScanPosition,   &
                                modelInterpTerrain, obsGlobalMarker, obsFlags,            &
                                cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                                atmScatteringIndexObs, burpFileSatId, RESETQC, modLSQ, lastHeader)
      else
        write(*,*) 'midas-bgckMW: instName = ', instName
        call utl_abort('midas-bgckMW: unknown instName')
      end if
      !###############################################################################
      ! STEP 5) Accumuler Les statistiques sur les rejets
      !###############################################################################
      call mwbg_qcStats(instName, qcIndicator, obsChannels, sensorIndex,       &
                        actualNumChannel, numObsToProcess, tvs_satelliteName(1:tvs_nsensors), &
                        .FALSE.)

      !###############################################################################
      ! STEP 6) Update Flags and obs in obsspace data
      !###############################################################################
      call mwbg_updateObsSpaceAfterQc(obsSpaceData, sensorIndex, headerIndex, obsTb, obsFlags, &
                                        cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                                        atmScatteringIndexObs, atmScatteringIndexFG, &
                                        obsGlobalMarker, newInformationFlag)

    end do HEADER
    !###############################################################################
    ! STEP 7) Print the statistics in listing file 
    !###############################################################################
    call mwbg_qcStats(instName, qcIndicator, obsChannels, sensorIndex,              &
                      actualNumChannel, numObsToProcess, tvs_satelliteName(1:tvs_nsensors), & 
                      .TRUE.)

    call utl_tmg_stop(118)

  end subroutine mwbg_bgCheckMW 

end module bgckmicrowave_mod

