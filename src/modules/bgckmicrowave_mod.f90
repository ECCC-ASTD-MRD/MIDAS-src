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
  real    :: mwbg_minSiOverWaterThreshold ! for AMSUB/MHS
  real    :: mwbg_maxSiOverWaterThreshold ! for AMSUB/MHS
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
  real                          :: minSiOverWaterThreshold       ! min scattering index over water for AMSUB/MHS
  real                          :: maxSiOverWaterThreshold       ! max scattering index over water for AMSUB/MHS
  real                          :: cloudySiThresholdBcorr        !
  logical                       :: useUnbiasedObsForClw          !
  logical                       :: RESETQC                       ! reset Qc flags option
  logical                       :: modLSQ                        !
  logical                       :: debug                         ! debug mode
  logical                       :: skipTestArr(mwbg_maxNumTest)  ! array to set to skip the test


  namelist /nambgck/instName, clwQcThreshold, &
                    useUnbiasedObsForClw, debug, RESETQC,  &
                    cloudyClwThresholdBcorr, modLSQ, &
                    minSiOverWaterThreshold, maxSiOverWaterThreshold, &
                    cloudySiThresholdBcorr, skipTestArr
                    

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
    minSiOverWaterThreshold = -10.0
    maxSiOverWaterThreshold = 30.0
    cloudySiThresholdBcorr  = 5.0
    RESETQC                 = .false.
    modLSQ                  = .false.
    skipTestArr(:)          = .false.

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
    mwbg_minSiOverWaterThreshold = minSiOverWaterThreshold
    mwbg_maxSiOverWaterThreshold = maxSiOverWaterThreshold
    mwbg_cloudySiThresholdBcorr = cloudySiThresholdBcorr

  end subroutine mwbg_init 

  !--------------------------------------------------------------------------
  ! ISRCHEQI function
  !--------------------------------------------------------------------------
  function ISRCHEQI(KLIST, KENTRY) result(ISRCHEQI_out)
    !OBJET          Rechercher un element dans une liste (valeurs entieres).
    !ARGUMENTS      - indx    - output -  position de l'element recherche:
    !                                   =0, element introuvable,
    !                                   >0, position de l'element trouve,
    !               - klist   - input  -  la liste
    !               - kentry  - input  -  l'element recherche

    implicit none
    
    ! Arguments:
    integer, intent(in) ::  KLIST(:)
    integer, intent(in) ::  KENTRY
    integer :: ISRCHEQI_out

    ! Locals:
    integer :: KLEN, JI

    ISRCHEQI_out = 0
    klen = size(KLIST)
    do JI = 1, KLEN
       if ( KLIST(JI) == KENTRY ) then
          ISRCHEQI_out = JI
          return
       end if
    end do

  end function ISRCHEQI

  !--------------------------------------------------------------------------
  ! extractParamForGrodyRun
  !--------------------------------------------------------------------------  
  subroutine extractParamForGrodyRun(obsChannels2D, obsTb2D, ompTb2D, obsTbBiasCorr2D, numObsToProcess, actualNumChannel, &
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
    integer,     intent(in)               :: actualNumChannel                                   ! nombre de canaux des observations 
    integer,     intent(in)               :: numObsToProcess                                    ! nombre de tovs
    integer,     intent(in)               :: obsChannels2D(actualNumChannel,numObsToProcess)    ! observations channels
    real,        intent(in)               :: obsTb2D(actualNumChannel,numObsToProcess)          ! radiances
    real,        intent(in)               :: ompTb2D(actualNumChannel,numObsToProcess)          ! radiances o-p
    real,        intent(in)               :: obsTbBiasCorr2D(actualNumChannel,numObsToProcess)  ! correction aux radiances
    real,        intent(out)              :: tb23(numObsToProcess)                              ! radiance frequence 23 Ghz   
    real,        intent(out)              :: tb31(numObsToProcess)                              ! radiance frequence 31 Ghz
    real,        intent(out)              :: tb50(numObsToProcess)                              ! radiance frequence 50 Ghz  
    real,        intent(out)              :: tb53(numObsToProcess)                              ! radiance frequence 53 Ghz  
    real,        intent(out)              :: tb89(numObsToProcess)                              ! radiance frequence 89 Ghz  
    real,        intent(out)              :: tb23FG(numObsToProcess)                            ! radiance frequence 23 Ghz   
    real,        intent(out)              :: tb31FG(numObsToProcess)                            ! radiance frequence 31 Ghz
    real,        intent(out)              :: tb50FG(numObsToProcess)                            ! radiance frequence 50 Ghz  
    real,        intent(out)              :: tb53FG(numObsToProcess)                            ! radiance frequence 53 Ghz  
    real,        intent(out)              :: tb89FG(numObsToProcess)                            ! radiance frequence 89 Ghz        

    ! Locals
    integer                               :: channelval
    integer                               :: nDataIndex
    integer                               :: nChannelIndex

    do nDataIndex=1,numObsToProcess
      do nChannelIndex=1,actualNumChannel
        channelval = obsChannels2D(nChannelIndex,nDataIndex)
        if ( obsTb2D(nChannelIndex,nDataIndex) /= mwbg_realMissing ) then
          if ( obsTbBiasCorr2D(nChannelIndex,nDataIndex) /= mwbg_realMissing ) then
            if ( channelval == 28 ) tb23(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex) &
                 - obsTbBiasCorr2D(nChannelIndex,nDataIndex)
            if ( channelval == 29 ) tb31(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex) &
                 - obsTbBiasCorr2D(nChannelIndex,nDataIndex)
            if ( channelval == 30 ) tb50(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex) &
                 - obsTbBiasCorr2D(nChannelIndex,nDataIndex)
            if ( channelval == 32 ) tb53(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex) &
                 - obsTbBiasCorr2D(nChannelIndex,nDataIndex)
            if ( channelval == 42 ) tb89(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex) &
                 - obsTbBiasCorr2D(nChannelIndex,nDataIndex)
          else
            if ( channelval == 28 ) tb23(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex)
            if ( channelval == 29 ) tb31(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex)
            if ( channelval == 30 ) tb50(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex)
            if ( channelval == 32 ) tb53(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex)
            if ( channelval == 42 ) tb89(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex)
          end if

          if ( channelval == 28 ) tb23FG(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex) - &
                                                         ompTb2D(nChannelIndex,nDataIndex)
          if ( channelval == 29 ) tb31FG(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex) - &
                                                         ompTb2D(nChannelIndex,nDataIndex)
          if ( channelval == 30 ) tb50FG(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex) - &
                                                         ompTb2D(nChannelIndex,nDataIndex)
          if ( channelval == 32 ) tb53FG(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex) - &
                                                         ompTb2D(nChannelIndex,nDataIndex)
          if ( channelval == 42 ) tb89FG(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex) - &
                                                         ompTb2D(nChannelIndex,nDataIndex)
        else
          if ( channelval == 28 ) tb23(nDataIndex) = 0.
          if ( channelval == 29 ) tb31(nDataIndex) = 0.
          if ( channelval == 30 ) tb50(nDataIndex) = 0.
          if ( channelval == 32 ) tb53(nDataIndex) = 0.
          if ( channelval == 42 ) tb89(nDataIndex) = 0.

          if ( channelval == 28 ) tb23FG(nDataIndex) = 0.  
          if ( channelval == 29 ) tb31FG(nDataIndex) = 0. 
          if ( channelval == 30 ) tb50FG(nDataIndex) = 0. 
          if ( channelval == 32 ) tb53FG(nDataIndex) = 0. 
          if ( channelval == 42 ) tb89FG(nDataIndex) = 0. 
        end if
      end do
    end do

  end subroutine extractParamForGrodyRun

  !--------------------------------------------------------------------------
  ! extractParamForBennartzRun
  !--------------------------------------------------------------------------  
  subroutine extractParamForBennartzRun(obsChannels2D, obsTb2D, btClear2D, ompTb2D, obsTbBiasCorr2D, numObsToProcess, actualNumChannel, &
                                        tb89, tb150, tb1831, tb1832, tb1833, &
                                        tb89FG, tb150FG, tb89FgClear, tb150FgClear)

    !:Purpose: Extract Parameters required to run bennaertz for required channels:
    !          extract required channels:        
    !          89 Ghz = AMSU-B 1 = channel #43
    !         150 Ghz = AMSU-B 2 = channel #44

    implicit none
    ! Arguments
    integer,     intent(in)               :: actualNumChannel                                   ! nombre de canaux des observations 
    integer,     intent(in)               :: numObsToProcess                                    ! nombre de tovs    
    integer,     intent(in)               :: obsChannels2D(actualNumChannel,numObsToProcess)    ! observations channels
    real,        intent(in)               :: obsTb2D(actualNumChannel,numObsToProcess)          ! radiances
    real,        intent(in)               :: btClear2D(actualNumChannel,numObsToProcess)        ! clear-sky radiances from background
    real,        intent(in)               :: ompTb2D(actualNumChannel,numObsToProcess)          ! radiances o-p
    real,        intent(in)               :: obsTbBiasCorr2D(actualNumChannel,numObsToProcess)  ! correction aux radiances
    real,        intent(out)              :: tb89(numObsToProcess)                              ! 89GHz radiance from observation
    real,        intent(out)              :: tb150(numObsToProcess)                             ! 150GHz radiance from observation
    real,        intent(out)              :: tb1831(numObsToProcess)                            ! 183GHz radiance from observation
    real,        intent(out)              :: tb1832(numObsToProcess)                            ! 183GHz radiance from observation
    real,        intent(out)              :: tb1833(numObsToProcess)                            ! 183GHz radiance from observation
    real,        intent(out)              :: tb89FG(numObsToProcess)                            ! 89GHz radiance from background
    real,        intent(out)              :: tb150FG(numObsToProcess)                           ! 150GHz radiance from background
    real,        intent(out)              :: tb89FgClear(numObsToProcess)                       ! 89GHz clear-sky radiance from background
    real,        intent(out)              :: tb150FgClear(numObsToProcess)                      ! 150GHz clear-sky radiance from background
    ! Locals
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: channelval

    do nDataIndex=1,numObsToProcess
      do nChannelIndex=1,actualNumChannel
        channelval = obsChannels2D(nChannelIndex,nDataIndex)
        if ( obsTb2D(nChannelIndex,nDataIndex) /= mwbg_realMissing ) then
          if ( obsTbBiasCorr2D(nChannelIndex,nDataIndex) /= mwbg_realMissing ) then
            if ( channelval == 43 ) tb89(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex) &
                 - obsTbBiasCorr2D(nChannelIndex,nDataIndex)
            if ( channelval == 44 ) tb150(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex) &
                 - obsTbBiasCorr2D(nChannelIndex,nDataIndex)
            if ( channelval == 45 ) tb1831(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex) &
                 - obsTbBiasCorr2D(nChannelIndex,nDataIndex)
            if ( channelval == 46 ) tb1832(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex) &
                 - obsTbBiasCorr2D(nChannelIndex,nDataIndex)
            if ( channelval == 47 ) tb1833(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex) &
                 - obsTbBiasCorr2D(nChannelIndex,nDataIndex)
          else
            if ( channelval == 43 ) tb89(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex)
            if ( channelval == 44 ) tb150(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex)
            if ( channelval == 45 ) tb1831(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex)
            if ( channelval == 46 ) tb1832(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex)
            if ( channelval == 47 ) tb1833(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex)
          end if

          if ( channelval == 43 ) tb89FG(nDataIndex)  = obsTb2D(nChannelIndex,nDataIndex) - &
                                                          ompTb2D(nChannelIndex,nDataIndex)
          if ( channelval == 44 ) tb150FG(nDataIndex) = obsTb2D(nChannelIndex,nDataIndex) - &
                                                          ompTb2D(nChannelIndex,nDataIndex)
          
        else
          if ( channelval == 43 ) tb89(nDataIndex) = 0.
          if ( channelval == 44 ) tb150(nDataIndex) = 0.
          if ( channelval == 45 ) tb1831(nDataIndex) = 0.
          if ( channelval == 46 ) tb1832(nDataIndex) = 0.
          if ( channelval == 47 ) tb1833(nDataIndex) = 0.

          if ( channelval == 43 ) tb89FG(nDataIndex) = 0.
          if ( channelval == 44 ) tb150FG(nDataIndex) = 0.
        end if

        if (btClear2D(nChannelIndex,nDataIndex) /= mwbg_realMissing) then
          if (channelval == 43) tb89FgClear(nDataIndex) = btClear2D(nChannelIndex,nDataIndex)
          if (channelval == 44) tb150FgClear(nDataIndex) = btClear2D(nChannelIndex,nDataIndex)
        else
          if (channelval == 43) tb89FgClear(nDataIndex) = 0.0
          if (channelval == 44) tb150FgClear(nDataIndex) = 0.0       
        end if

      end do
    end do

  end subroutine extractParamForBennartzRun

  !--------------------------------------------------------------------------
  ! amsuABTest10RttovRejectCheck
  !--------------------------------------------------------------------------
  subroutine amsuABTest10RttovRejectCheck(obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, RESETQC, &
                                          stnId, obsFlags2D, qcIndicator)

    !:Purpose:               10) test 10: RTTOV reject check (single)
    !                        Rejected datum flag has bit #9 on.

    implicit none
    ! Arguments
    integer,     intent(in)                :: actualNumChannel                                ! nombre de canaux des observations 
    integer,     intent(in)                :: numObsToProcess                                 ! nombre de tovs    
    integer,     intent(in)                :: obsChannels2D(actualNumChannel,numObsToProcess) ! observations channels
    integer,     intent(in)                :: sensorIndex                                     ! numero de satellite (i.e. indice) 
    logical,     intent(in)                :: RESETQC                                         ! yes or not reset QC flag
    character *9, intent(in)               :: stnId                                           ! identificateur du satellite
    integer,     intent(inout)             :: obsFlags2D(actualNumChannel,numObsToProcess)    ! marqueur de radiance 
    integer,     intent(inout)             :: qcIndicator(actualNumChannel,numObsToProcess)   ! indicateur du QC par canal
    ! Locals
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex
    integer                                :: IBIT

    if (.NOT.RESETQC) then
      testIndex = 10
      do nDataIndex=1,numObsToProcess
        do nChannelIndex=1,actualNumChannel
          if (obsChannels2D(nChannelIndex,nDataIndex) /= 20) then
            IBIT = AND(obsFlags2D(nChannelIndex,nDataIndex), 2**9)
            if (IBIT /= 0) then
              qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
              obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**7)
              rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                   rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1
              if ( mwbg_DEBUG ) then
                write(*,*)stnId(2:9),' RTTOV REJECT.', &
                          'CHANNEL=', obsChannels2D(nChannelIndex,nDataIndex), &
                          ' obsFlags2D= ',obsFlags2D(nChannelIndex,nDataIndex)
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
  subroutine amsuABTest1TopographyCheck(obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                        modelInterpTerrain, channelForTopoFilter, altitudeForTopoFilter, obsFlags2D, &
                                        qcIndicator)

    !:Purpose:               1) test 1: Topography check (partial)
    !                        Channel 6 is rejected for topography >  250m.
    !                        Channel 7 is rejected for topography > 2000m.

    implicit none
    ! Arguments
    integer,     intent(in)                :: actualNumChannel                                ! nombre de canaux des observations 
    integer,     intent(in)                :: numObsToProcess                                 ! nombre de tovs
    integer,     intent(in)                :: obsChannels2D(actualNumChannel,numObsToProcess) ! observations channels
    integer,     intent(in)                :: sensorIndex                                     ! numero de satellite (i.e. indice) 
    character *9,intent(in)                :: stnId                                           ! identificateur du satellite
    integer,     intent(in)                :: channelForTopoFilter(:)                         ! channel list for filter
    real,        intent(in)                :: altitudeForTopoFilter(:)                        ! altitude threshold
    real,        intent(in)                :: modelInterpTerrain(numObsToProcess)             ! topo aux point d'obs
    integer,     intent(inout)             :: obsFlags2D(actualNumChannel,numObsToProcess)    ! marqueur de radiance 
    integer,     intent(inout)             :: qcIndicator(actualNumChannel,numObsToProcess)   ! indicateur du QC par canal
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
      do nDataIndex=1,numObsToProcess
        do nChannelIndex=1,actualNumChannel
          if (obsChannels2D(nChannelIndex,nDataIndex) == channelForTopoFilter(indexFilteringTest)) then
            if (modelInterpTerrain(nDataIndex) >= altitudeForTopoFilter(indexFilteringTest)) then
              qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
              obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
              obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**18)
              rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                   rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1
              if ( mwbg_DEBUG ) then
                write(*,*)stnId(2:9),' TOPOGRAPHY REJECT.', &
                          'CHANNEL=', obsChannels2D(nChannelIndex,nDataIndex), &
                          ' TOPO= ',modelInterpTerrain(nDataIndex)
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
  subroutine amsuABTest2LandSeaQualifierCheck(obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                              landQualifierIndice, obsFlags2D, qcIndicator)

    !:Purpose:                      2) test 2: "Land/sea qualifier" code check (full)
    !                                  allowed values are: 0 land, 1 sea, 2 coast.

    implicit none
    ! Arguments
    integer,     intent(in)                :: actualNumChannel                                ! nombre de canaux des observations 
    integer,     intent(in)                :: numObsToProcess                                 ! nombre de tovs    
    integer,     intent(in)                :: obsChannels2D(actualNumChannel,numObsToProcess) ! observations channels
    integer,     intent(in)                :: sensorIndex                                     ! numero de satellite (i.e. indice) 
    character *9,intent(in)                :: stnId                                           ! identificateur du satellite
    integer,     intent(in)                :: landQualifierIndice(numObsToProcess)            ! land sea qualifier
    integer,     intent(inout)             :: obsFlags2D(actualNumChannel,numObsToProcess)    ! marqueur de radiance 
    integer,     intent(inout)             :: qcIndicator(actualNumChannel,numObsToProcess)   ! indicateur du QC par canal
    ! Locals
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex

  
    testIndex = 2
    do nDataIndex=1,numObsToProcess
      if (landQualifierIndice(nDataIndex) < 0  .or. &
          landQualifierIndice(nDataIndex) > 2) then
        do nChannelIndex=1,actualNumChannel
          qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
          obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
          obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**7)
          rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
            rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) + 1
        end do
        if ( mwbg_DEBUG ) then
          write(*,*) stnId(2:9),'LAND/SEA QUALifIER CODE', &
                   ' REJECT. landQualifierIndice=', landQualifierIndice(nDataIndex)
        end if
      end if
    end do

  end subroutine amsuABTest2LandSeaQualifierCheck

  !--------------------------------------------------------------------------
  !  amsuABTest3TerrainTypeCheck
  !--------------------------------------------------------------------------
  subroutine amsuABTest3TerrainTypeCheck(obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                         terrainTypeIndice, obsFlags2D, qcIndicator)

    !:Purpose:                     3) test 3: "Terrain type" code check (full)
    !                                 allowed values are: -1 missing, 0 sea-ice, 1 snow on land.

    implicit none
    ! Arguments
    integer,     intent(in)               :: actualNumChannel                                ! nombre de canaux des observations 
    integer,     intent(in)               :: numObsToProcess                                 ! nombre de tovs    
    integer,     intent(in)               :: obsChannels2D(actualNumChannel,numObsToProcess) ! observations channels
    integer,     intent(in)               :: sensorIndex                                     ! numero de satellite (i.e. indice) 
    character *9,intent(in)               :: stnId                                           ! identificateur du satellite
    integer,     intent(in)               :: terrainTypeIndice(numObsToProcess)              ! terrain type
    integer,     intent(inout)            :: obsFlags2D(actualNumChannel,numObsToProcess)    ! marqueur de radiance 
    integer,     intent(inout)            :: qcIndicator(actualNumChannel,numObsToProcess)   ! indicateur du QC par canal
    ! Locals
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex

 
    testIndex = 3
    do nDataIndex=1,numObsToProcess
      if ( terrainTypeIndice(nDataIndex) /= mwbg_intMissing ) then
        if (terrainTypeIndice(nDataIndex) <  0  .or. &
            terrainTypeIndice(nDataIndex) >  1) then
          do nChannelIndex=1,actualNumChannel
            qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**7)
            rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
              rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) + 1
          end do
          if ( mwbg_debug ) then
            write(*,*)stnId(2:9),'TERRAIN type CODE', &
                     ' REJECT. TERRAIN=', terrainTypeIndice(nDataIndex)
          end if
        end if
      end if
    end do

  end subroutine amsuABTest3TerrainTypeCheck

  !--------------------------------------------------------------------------
  ! amsuABTest4FieldOfViewCheck 
  !--------------------------------------------------------------------------
  subroutine amsuABTest4FieldOfViewCheck(obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                         satScanPosition, maxScanAngleAMSU, obsFlags2D, qcIndicator)

    !:Purpose:                          4) test 4: Field of view number check (full)
    !                                      Field of view acceptable range is [1,maxScanAngleAMSU]  for AMSU footprints.
    implicit none
    ! Arguments
    integer,     intent(in)               :: actualNumChannel                                 ! nombre de canaux des observations 
    integer,     intent(in)               :: numObsToProcess                                  ! nombre de tovs    
    integer,     intent(in)               :: obsChannels2D(actualNumChannel,numObsToProcess)  ! observations channels
    integer,     intent(in)               :: sensorIndex                                      ! numero de satellite (i.e. indice) 
    character *9,intent(in)               :: stnId                                            ! identificateur du satellite
    integer,     intent(in)               :: satScanPosition(numObsToProcess)                 ! position sur le "scan" 
    integer,     intent(in)               :: maxScanAngleAMSU                                 ! max scan angle 
    integer,     intent(inout)            :: obsFlags2D(actualNumChannel,numObsToProcess)     ! marqueur de radiance 
    integer,     intent(inout)            :: qcIndicator(actualNumChannel,numObsToProcess)    ! indicateur du QC par canal
    ! Locals
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex

    testIndex = 4
    do nDataIndex=1,numObsToProcess
      do nChannelIndex=1,actualNumChannel
        if ( satScanPosition(nDataIndex) < 1 .or. &
            satScanPosition(nDataIndex) > maxScanAngleAMSU ) then
          qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
          obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
          obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**7)
          rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
            rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) + 1
          if ( mwbg_debug ) then
            write(*,*)stnId(2:9),'FIELD OF VIEW NUMBER', &
                      ' REJECT. FIELD OF VIEW= ', satScanPosition(nDataIndex)
          end if
        end if
      end do
    end do
  
  end subroutine amsuABTest4FieldOfViewCheck 
  
  !--------------------------------------------------------------------------
  ! amsuABTest5ZenithAngleCheck 
  !--------------------------------------------------------------------------
  subroutine amsuABTest5ZenithAngleCheck(obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                         satZenithAngle, obsFlags2D, qcIndicator)

    !:Purpose:                   5) test 5: Satellite zenith angle check (full)
    !                               Satellite zenith angle acceptable range is [0.,60.].
    implicit none
    ! Arguments
    integer,     intent(in)               :: actualNumChannel                                 ! nombre de canaux des observations 
    integer,     intent(in)               :: numObsToProcess                                  ! nombre de tovs
    integer,     intent(in)               :: obsChannels2D(actualNumChannel,numObsToProcess)  ! observations channels
    integer,     intent(in)               :: sensorIndex                                      ! numero de satellite (i.e. indice) 
    character *9,intent(in)               :: stnId                                            ! identificateur du satellite
    real,        intent(in)               :: satZenithAngle(numObsToProcess)                  ! satellite zenith angle 
    integer,     intent(inout)            :: obsFlags2D(actualNumChannel,numObsToProcess)     ! marqueur de radiance 
    integer,     intent(inout)            :: qcIndicator(actualNumChannel,numObsToProcess)    ! indicateur du QC par canal
    ! Locals
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex


    testIndex = 5
    do nDataIndex=1,numObsToProcess
      if ( satZenithAngle(nDataIndex) /= mwbg_realMissing ) then
        if (satZenithAngle(nDataIndex) <  0.  .or. &
            satZenithAngle(nDataIndex) > 60.) then
          do nChannelIndex=1,actualNumChannel
            qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**7)
            rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
               rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) + 1
          end do
          if ( mwbg_debug ) then
            write(*,*)stnId(2:9),' SATELLITE ZENITH ANGLE', &
                      ' REJECT. satZenithAngle= ', &
                      satZenithAngle(nDataIndex)
          end if
        end if
      end if
    end do

  end subroutine amsuABTest5ZenithAngleCheck 

  !--------------------------------------------------------------------------
  ! amsuABTest6ZenAngleAndFovConsistencyCheck 
  !--------------------------------------------------------------------------
  subroutine amsuABTest6ZenAngleAndFovConsistencyCheck(obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, &
                                                       stnId, satZenithAngle, ZANGL, satScanPosition, &
                                                       maxScanAngleAMSU, obsFlags2D, qcIndicator)

    !:Purpose:                            6) test 6: "Sat. zenith angle"/"field of view" consistency check.  (full)
    !                                        Acceptable difference between "Satellite zenith angle"  and
    !                                       "approximate angle computed from field of view number" is 1.8 degrees.

    implicit none
    ! Arguments
    integer,     intent(in)               :: actualNumChannel                                 ! nombre de canaux des observations 
    integer,     intent(in)               :: numObsToProcess                                  ! nombre de tovs    
    integer,     intent(in)               :: obsChannels2D(actualNumChannel,numObsToProcess)  ! observations channels
    integer,     intent(in)               :: sensorIndex                                      ! numero de satellite (i.e. indice) 
    character *9,intent(in)               :: stnId                                            ! identificateur du satellite
    real,        intent(in)               :: satZenithAngle(numObsToProcess)                  ! satellite zenith angle 
    real,        intent(in)               :: ZANGL                                            ! satellite constant param
    integer,     intent(in)               :: satScanPosition(numObsToProcess)                 ! position sur le "scan" 
    integer,     intent(in)               :: maxScanAngleAMSU                                 ! max scan angle 
    integer,     intent(inout)            :: obsFlags2D(actualNumChannel,numObsToProcess)     ! marqueur de radiance 
    integer,     intent(inout)            :: qcIndicator(actualNumChannel,numObsToProcess)    ! indicateur du QC par canal
    ! Locals
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex
    real                                  :: APPROXIM 
    real                                  :: ANGDif 

    testIndex = 6
    do nDataIndex=1,numObsToProcess
      if ( satZenithAngle (nDataIndex) /=  mwbg_realMissing   .and. &
           satScanPosition(nDataIndex) /=  mwbg_intMissing  ) then
        APPROXIM = ABS((satScanPosition(nDataIndex)-maxScanAngleAMSU/2.-0.5)*ZANGL)
        ANGDif = ABS(satZenithAngle(nDataIndex)-APPROXIM)
        if ( ANGDif > 1.8 ) then 
          do nChannelIndex=1,actualNumChannel
            qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**7)
            rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
               rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) + 1
          end do
          if ( mwbg_debug ) then
            write(*,*)stnId(2:9),' ANGLE/FIELD OF VIEW', &
                      ' INCONSISTENCY REJECT. satZenithAngle= ', &
                      satZenithAngle(nDataIndex), ' FIELD OF VIEW= ',satScanPosition(nDataIndex), &
                      ' ANGDif= ',ANGDif  
          end if
        end if
      end if
    end do

  end subroutine amsuABTest6ZenAngleAndFovConsistencyCheck

  !--------------------------------------------------------------------------
  !  amsuABTest7landSeaQualifyerAndModelLandSeaConsistencyCheck
  !--------------------------------------------------------------------------
  subroutine amsuABTest7landSeaQualifyerAndModelLandSeaConsistencyCheck(obsChannels2D, sensorIndex, actualNumChannel, &
                                                                        numObsToProcess, stnId, modelInterpGroundIce, &
                                                                        landQualifierIndice, obsFlags2D, qcIndicator)

    !:Purpose: test 7: "Land/sea qual."/"model land/sea" consistency check (full)
    !          Acceptable conditions are:
    !            - both over ocean (landQualifierIndice=1; mg<0.01), new threshold 0.20, jh dec 2000,
    !            - both over land  (landQualifierIndice=0; mg>0.80), new threshold 0.50, jh dec 2000.
    !            - Other conditions are unacceptable.


    implicit none
    ! Arguments
    integer,     intent(in)                :: actualNumChannel                                ! nombre de canaux des observations 
    integer,     intent(in)                :: numObsToProcess                                 ! nombre de tovs
    integer,     intent(in)                :: obsChannels2D(actualNumChannel,numObsToProcess) ! observations channels
    integer,     intent(in)                :: sensorIndex                                     ! numero de satellite (i.e. indice) 
    character *9,intent(in)                :: stnId                                           ! identificateur du satellite
    real,        intent(in)                :: modelInterpGroundIce(numObsToProcess)           ! glace mer 
    integer,     intent(in)                :: landQualifierIndice(numObsToProcess)            ! land sea qualifyer 
    integer,     intent(inout)             :: obsFlags2D(actualNumChannel,numObsToProcess)    ! marqueur de radiance 
    integer,     intent(inout)             :: qcIndicator(actualNumChannel,numObsToProcess)   ! indicateur du QC par canal
    ! Locals
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex

    testIndex = 7
    do nDataIndex=1,numObsToProcess
      if (landQualifierIndice (nDataIndex) /= mwbg_intMissing) then
        if (landQualifierIndice(nDataIndex) == 1 .and. &
            modelInterpGroundIce(nDataIndex) < 0.20) then
        else if (landQualifierIndice(nDataIndex) == 0 .and. &
                 modelInterpGroundIce(nDataIndex) > 0.50) then
        else
          do nChannelIndex=1,actualNumChannel
            qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**7)
            rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
               rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) + 1
          end do
          if ( mwbg_debug ) then
            write(*,*)stnId(2:9),' LAND/SEA QUALifIER', &
                      ' INCONSISTENCY REJECT. landQualifierIndice= ', &
                      landQualifierIndice(nDataIndex), ' MODEL MASK= ',modelInterpGroundIce(nDataIndex)
          end if
        end if
      end if
    end do

  end subroutine amsuABTest7landSeaQualifyerAndModelLandSeaConsistencyCheck 

  !--------------------------------------------------------------------------
  !  amsuABTest9UncorrectedTbCheck
  !--------------------------------------------------------------------------
  subroutine amsuABTest9UncorrectedTbCheck(obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                           RESETQC, obsFlags2D, qcIndicator)

    !:Purpose:                  9) test 9: Uncorrected Tb check (single)
    !                              Uncorrected datum (flag bit #6 off). In this case switch bit 11 ON.

    implicit none
    ! Arguments
    integer,     intent(in)               :: actualNumChannel                                 ! nombre de canaux des observations 
    integer,     intent(in)               :: numObsToProcess                                  ! nombre de tovs    
    integer,     intent(in)               :: obsChannels2D(actualNumChannel,numObsToProcess)  ! observations channels
    integer,     intent(in)               :: sensorIndex                                      ! numero de satellite (i.e. indice) 
    character *9,intent(in)               :: stnId                                            ! identificateur du satellite
    logical,     intent(in)               :: RESETQC                                          ! yes or not reset QC flag
    integer,     intent(inout)            :: obsFlags2D(actualNumChannel,numObsToProcess)     ! marqueur de radiance 
    integer,     intent(inout)            :: qcIndicator(actualNumChannel,numObsToProcess)    ! indicateur du QC par canal
    ! Locals
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex
    integer                               :: IBIT


    if (.not. RESETQC) then
      testIndex = 9
      do nDataIndex=1,numObsToProcess
        do nChannelIndex=1,actualNumChannel
          if (obsChannels2D(nChannelIndex,nDataIndex) /= 20) then
            IBIT = AND(obsFlags2D(nChannelIndex,nDataIndex), 2**6)
            if ( IBIT == 0  ) then
              qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
              obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**11)
              rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                    rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1
              if ( mwbg_debug ) then
                write(*,*)stnId(2:9),' UNCORRECTED TB REJECT.', &
                           'CHANNEL=', obsChannels2D(nChannelIndex,nDataIndex), ' obsFlags2D= ',obsFlags2D(nChannelIndex,nDataIndex)
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
  subroutine amsuABTest11RadianceGrossValueCheck(obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                                 obsTb2D, GROSSMIN, GROSSMAX, obsFlags2D, qcIndicator)

    !:Purpose:                     11) test 11: Radiance observation "Gross" check (single) 
    !                                               Change this test from full to single. jh nov 2000.

    implicit none
    ! Arguments
    integer,     intent(in)               :: actualNumChannel                                 ! nombre de canaux des observations 
    integer,     intent(in)               :: numObsToProcess                                  ! nombre de tovs    
    integer,     intent(in)               :: obsChannels2D(actualNumChannel,numObsToProcess)  ! observations channels
    integer,     intent(in)               :: sensorIndex                                      ! numero de satellite (i.e. indice) 
    character *9,intent(in)               :: stnId                                            ! identificateur du satellite
    real,        intent(in)               :: obsTb2D(actualNumChannel,numObsToProcess)        ! radiances 
    real,        intent(in)               :: GROSSMIN(:)                                      ! Gross val min 
    real,        intent(in)               :: GROSSMAX(:)                                      ! Gross val max 
    integer,     intent(inout)            :: obsFlags2D(actualNumChannel,numObsToProcess)     ! marqueur de radiance 
    integer,     intent(inout)            :: qcIndicator(actualNumChannel,numObsToProcess)    ! indicateur du QC par canal
    ! Locals
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex
    logical                               :: GROSSERROR

    testIndex = 11
    do nDataIndex=1,numObsToProcess
      GROSSERROR = .FALSE.
      do nChannelIndex=1,actualNumChannel
        if (obsChannels2D(nChannelIndex,nDataIndex) /= 20 .and. &
            obsChannels2D(nChannelIndex,nDataIndex) >=  1 .and. &
            obsChannels2D(nChannelIndex,nDataIndex) <=  actualNumChannel) then  
          if (obsTb2D(nChannelIndex,nDataIndex) /= mwbg_realMissing .and. &
              (obsTb2D(nChannelIndex,nDataIndex) < GROSSMIN(obsChannels2D(nChannelIndex,nDataIndex)).or. &
               obsTb2D(nChannelIndex,nDataIndex) > GROSSMAX(obsChannels2D(nChannelIndex,nDataIndex)))) then
            GROSSERROR = .TRUE.
            qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**7)
            rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                   rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1
            if ( mwbg_debug ) then
              write(*,*)stnId(2:9),' GROSS CHECK REJECT.', &
                        'CHANNEL=', obsChannels2D(nChannelIndex,nDataIndex), ' TB= ',obsTb2D(nChannelIndex,nDataIndex)
            end if
          end if
        end if
      end do
    end do

  end subroutine amsuABTest11RadianceGrossValueCheck 
  
  !--------------------------------------------------------------------------
  ! amsuaTest12GrodyClwCheck
  !--------------------------------------------------------------------------
  subroutine amsuaTest12GrodyClwCheck(obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                      cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, landQualifierIndice, &
                                      ICLWREJ, obsFlags2D, qcIndicator)

    !:Purpose:                    12) test 12: Grody cloud liquid water check (partial)
    !                                 For Cloud Liquid Water > clwQcThreshold, reject AMSUA-A channels 1-5 and 15.

    implicit none
    ! Arguments
    integer,     intent(in)               :: actualNumChannel                                 ! nombre de canaux des observations 
    integer,     intent(in)               :: numObsToProcess                                  ! nombre de tovs
    integer,     intent(in)               :: obsChannels2D(actualNumChannel,numObsToProcess)  ! observations channels
    integer,     intent(in)               :: sensorIndex                                      ! numero de satellite (i.e. indice) 
    character *9,intent(in)               :: stnId                                            ! identificateur du satellite
    real,        intent(in)               :: cloudLiquidWaterPathObs(numObsToProcess)         ! retrieved cloud liquid water from observation
    real,        intent(in)               :: cloudLiquidWaterPathFG(numObsToProcess)          ! retrieved cloud liquid water from background
    integer,     intent(in)               :: landQualifierIndice(numObsToProcess)             ! land sea qualifyer 
    integer,     intent(in)               :: ICLWREJ(:)
    integer,     intent(inout)            :: obsFlags2D(actualNumChannel,numObsToProcess)     ! marqueur de radiance 
    integer,     intent(inout)            :: qcIndicator(actualNumChannel,numObsToProcess)    ! indicateur du QC par canal
    ! Locals
    integer                               :: channelval
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex
    integer                               :: INDXCAN 
    real                                  :: clwUsedForQC
    real                                  :: clwObsFGaveraged
    logical                               :: surfTypeIsWater 
    logical                               :: cldPredMissing

    testIndex = 12
    do nDataIndex=1,numObsToProcess
      if ( tvs_mwAllskyAssim ) then
        clwObsFGaveraged = 0.5 * (cloudLiquidWaterPathObs(nDataIndex) + cloudLiquidWaterPathFG(nDataIndex))
        clwUsedForQC = clwObsFGaveraged
        cldPredMissing = (cloudLiquidWaterPathObs(nDataIndex) == mwbg_realMissing .or. cloudLiquidWaterPathFG(nDataIndex) == mwbg_realMissing)
      else
        clwUsedForQC = cloudLiquidWaterPathObs(nDataIndex)
        cldPredMissing = (cloudLiquidWaterPathObs(nDataIndex) == mwbg_realMissing)
      end if

      surfTypeIsWater = (landQualifierIndice(nDataIndex) ==  1)

      if (.not. cldPredMissing) then
        if (clwUsedForQC > mwbg_clwQcThreshold) then
          do nChannelIndex=1,actualNumChannel
            INDXCAN = ISRCHEQI(ICLWREJ,obsChannels2D(nChannelIndex,nDataIndex))
            if ( INDXCAN /= 0 )  then
              qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
              obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
              obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**7)
              rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                       rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1
            end if
          end do
          if ( mwbg_debug ) then
            write(*,*)stnId(2:9),'Grody cloud liquid water check', &
                      ' REJECT. CLW= ',clwUsedForQC, ' SEUIL= ',mwbg_clwQcThreshold
          end if
        end if

        ! In all-sky mode, turn on bit=23 for channels in ICLWREJ(:) as 
        ! cloud-affected radiances over sea when there is mismatch between 
        ! cloudLiquidWaterPathObs and cloudLiquidWaterPathFG (to be used in gen_bias_corr)
        clwObsFGaveraged = 0.5 * (cloudLiquidWaterPathObs(nDataIndex) + cloudLiquidWaterPathFG(nDataIndex))
        IF (tvs_mwAllskyAssim .and. clwObsFGaveraged > mwbg_cloudyClwThresholdBcorr) then
          do nChannelIndex = 1,actualNumChannel
            INDXCAN = ISRCHEQI(ICLWREJ,obsChannels2D(nChannelIndex,nDataIndex))
            if ( INDXCAN /= 0 ) obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**23)
          end do
          if ( mwbg_debug ) then
            write(*,*) stnId(2:9),' Grody cloud liquid water check', &
                      ' cloud-affected obs. CLW= ',clwUsedForQC, ', threshold= ', &
                      mwbg_cloudyClwThresholdBcorr
          end if
        end if

      ! Reject surface sensitive observations over water, in all-sky mode, 
      ! if CLW is not retrieved, and is needed to define obs error.
      else if (tvs_mwAllskyAssim .and. surfTypeIsWater .and. cldPredMissing) then
        loopChannel: do nChannelIndex = 1, actualNumChannel
          channelval = obsChannels2D(nChannelIndex,nDataIndex)
          INDXCAN = ISRCHEQI(ICLWREJ,channelval)
          if ( INDXCAN /= 0 .and. oer_useStateDepSigmaObs(channelval,sensorIndex) ) then
            qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**7)
            rejectionCodArray(testIndex,channelval,sensorIndex) = &
                     rejectionCodArray(testIndex,channelval,sensorIndex)+ 1
          end if
        end do loopChannel

      end if
    end do

  end subroutine amsuaTest12GrodyClwCheck 

  !-------------------------------------------------------------------------
  ! amsubTest12DrynessIndexCheck
  !-------------------------------------------------------------------------
  subroutine amsubTest12DrynessIndexCheck(obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                          tb1831, tb1833, landQualifierIndice, modelInterpSeaIce, obsFlags2D, qcIndicator, &
                                          skipTestArr_opt)

    !:Purpose:  12) test 12: Dryness index check
    !           The difference between channels AMSUB-3 and AMSUB-5 is used as an indicator
    !           of "dryness" of the atmosphere. In extreme dry conditions, channels AMSUB-3 4 and 5
    !           are sensitive to the surface.
    !           Therefore, various thresholds are used to reject channels AMSUB-3 4 and 5 over land and ice

    implicit none

    ! Arguments:
    integer,     intent(in)               :: actualNumChannel                                 ! nombre de canaux des observations 
    integer,     intent(in)               :: numObsToProcess                                  ! nombre de tovs    
    integer,     intent(in)               :: obsChannels2D(actualNumChannel,numObsToProcess)  ! observations channels
    integer,     intent(in)               :: sensorIndex                                      ! numero de satellite (i.e. indice) 
    character *9,intent(in)               :: stnId                                            ! identificateur du satellite
    real,        intent(in)               :: tb1831(numObsToProcess)                          ! tb for channel  
    real,        intent(in)               :: tb1833(numObsToProcess)                          ! tb for channel  
    integer,     intent(in)               :: landQualifierIndice(:)                           ! mask terre-mer
    real,        intent(in)               :: modelInterpSeaIce(:)                             ! topo interpolated to obs point
    integer,     intent(inout)            :: obsFlags2D(actualNumChannel,numObsToProcess)     ! marqueur de radiance 
    integer,     intent(inout)            :: qcIndicator(actualNumChannel,numObsToProcess)    ! indicateur du QC par canal
    logical,     intent(in), optional     :: skipTestArr_opt(:)                               ! array to set to skip the test

    ! Locals:
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex
    real                                  :: drynessIndex
    logical, save                         :: firstCall = .true.

    testIndex = 12
    if (present(skipTestArr_opt)) then
      if (skipTestArr_opt(testIndex)) then
        if (firstCall) then
          firstCall = .false.
          write(*,*) 'amsubTest12DrynessIndexCheck: skipping this test'
        end if
        return
      end if
    end if

    do nDataIndex = 1, numObsToProcess
      drynessIndex = tb1831(nDataIndex) - tb1833(nDataIndex)
      do nChannelIndex = 1, actualNumChannel
        if ( .not. ((landQualifierIndice (nDataIndex) == 1) .and. &
                    (modelInterpSeaIce (nDataIndex) < 0.01)) ) then
          if (obsChannels2D(nChannelIndex,nDataIndex) == 45 .and. &
              drynessIndex > 0.) then
            qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**7)
            rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
              rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1
            if (mwbg_debug) then
              write(*,*)stnId(2:9),' DRYNESS INDEX REJECT.',        &
                       'CHANNEL=', obsChannels2D(nChannelIndex,nDataIndex), &
                       'INDEX= ',drynessIndex
            end if
          else if (obsChannels2D(nChannelIndex,nDataIndex) == 46 .and. &
                   drynessIndex > -10.) then
            qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**7)
            rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) =  &
            rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1
            if (mwbg_debug) then
              write(*,*)stnId(2:9),' DRYNESS INDEX REJECT.',       &
                      'CHANNEL=', obsChannels2D(nChannelIndex,nDataIndex),&
                      'INDEX= ',drynessIndex
            end if
          else if (obsChannels2D(nChannelIndex,nDataIndex) == 47 .and. &
                   drynessIndex > -20.) then
            qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**7)
            rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
            rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1
            if (mwbg_debug) then
              write(*,*)stnId(2:9),' DRYNESS INDEX REJECT.',       &
                       'CHANNEL=', obsChannels2D(nChannelIndex,nDataIndex),&
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
  subroutine amsuaTest13GrodyScatteringIndexCheck(obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                                  scatIndexOverWaterObs, landQualifierIndice, terrainTypeIndice, &
                                                  ISCATREJ, obsFlags2D, qcIndicator)

    !:Purpose:                  13) test 13: Grody scattering index check (partial)
    !                               For Scattering Index > 9, reject AMSUA-A channels 1-6 and 15.

    implicit none
    ! Arguments
    integer,     intent(in)                :: actualNumChannel                                ! nombre de canaux des observations 
    integer,     intent(in)                :: numObsToProcess                                 ! nombre de tovs    
    integer,     intent(in)                :: obsChannels2D(actualNumChannel,numObsToProcess) ! observations channels
    integer,     intent(in)                :: sensorIndex                                     ! numero de satellite (i.e. indice) 
    character *9,intent(in)                :: stnId                                           ! identificateur du satellite
    real,        intent(in)                :: scatIndexOverWaterObs(numObsToProcess)          ! scattering index 
    integer,     intent(in)                :: landQualifierIndice(numObsToProcess)            ! land sea qualifyer 
    integer,     intent(in)                :: terrainTypeIndice(numObsToProcess)              ! terrain type 
    integer,     intent(in)                :: ISCATREJ(:)
    integer,     intent(inout)             :: obsFlags2D(actualNumChannel,numObsToProcess)    ! marqueur de radiance 
    integer,     intent(inout)             :: qcIndicator(actualNumChannel,numObsToProcess)   ! indicateur du QC par canal
    ! Locals
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex
    integer                                :: INDXCAN 
    real                                   :: ZSEUILSCAT

    testIndex = 13
    ZSEUILSCAT = 9.0
    do nDataIndex=1,numObsToProcess
      if ( scatIndexOverWaterObs(nDataIndex) /=  mwbg_realMissing  ) then
        if (landQualifierIndice (nDataIndex) ==  1 .and. &
            terrainTypeIndice(nDataIndex) /=  0 .and. &   
            scatIndexOverWaterObs(nDataIndex) > ZSEUILSCAT) then
          do nChannelIndex=1,actualNumChannel
            INDXCAN = ISRCHEQI(ISCATREJ,obsChannels2D(nChannelIndex,nDataIndex))
            if ( INDXCAN /= 0 )  then
              qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
              obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
              obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**7)
              rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                       rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1
            end if
          end do
          if ( mwbg_debug ) then
            write(*,*)stnId(2:9),'Grody scattering index check', &
                       ' REJECT. scatIndexOverWaterObs= ',scatIndexOverWaterObs(nDataIndex), ' SEUIL= ',ZSEUILSCAT
          end if
        end if
      end if
    end do

  end subroutine amsuaTest13GrodyScatteringIndexCheck
 

  !--------------------------------------------------------------------------
  ! amsubTest13BennartzScatteringIndexCheck
  !--------------------------------------------------------------------------
  subroutine amsubTest13BennartzScatteringIndexCheck(obsChannels2D, sensorIndex, numObsToProcess, actualNumChannel, stnId, &
                                                     scatIndexOverWaterObs, scatIndexOverWaterFG, scatIndexOverLandObs, &
                                                     landQualifierIndice, modelInterpSeaIce, &
                                                     obsFlags2D, qcIndicator, chanIgnoreInAllskyGenCoeff, &
                                                     skipTestArr_opt)

    !:Purpose:                  13) test 13: Bennartz scattering index check (full)
    !                               For Scattering Index > 40 sea ice
    !                                                    > 15 sea
    !                                                    > 0 land reject all AMSUB Channels
    !

    implicit none
    ! Arguments
    integer,     intent(in)                :: numObsToProcess                                 ! nombre de tovs
    integer,     intent(in)                :: actualNumChannel                                ! nombre canaux de tovs    
    integer,     intent(in)                :: obsChannels2D(actualNumChannel,numObsToProcess) ! observations channels
    integer,     intent(in)                :: sensorIndex                                     ! numero de satellite (i.e. indice) 
    character *9,intent(in)                :: stnId                                           ! identificateur du satellite
    real,        intent(in)                :: scatIndexOverWaterObs(numObsToProcess)          ! scattering index over water from observation
    real,        intent(in)                :: scatIndexOverWaterFG(numObsToProcess)           ! scattering index over water from background
    real,        intent(in)                :: scatIndexOverLandObs(numObsToProcess)           ! scattering index over land
    integer,     intent(in)                :: landQualifierIndice(numObsToProcess)            ! land sea qualifyer 
    real,        intent(in)                :: modelInterpSeaIce(numObsToProcess)              ! glace de mer
    integer,     intent(inout)             :: obsFlags2D(actualNumChannel,numObsToProcess)    ! marqueur de radiance 
    integer,     intent(inout)             :: qcIndicator(actualNumChannel,numObsToProcess)   ! indicateur du QC par canal
    integer,     intent(in)                :: chanIgnoreInAllskyGenCoeff(:)                   ! channels to exclude from genCoeff
    logical,     intent(in), optional      :: skipTestArr_opt(:)                              ! array to set to skip the test


    ! Locals
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex
    integer                                :: channelval
    integer                                :: chanIndex
    real                                   :: ZSEUILSCATICE
    real                                   :: ZSEUILSCATL
    real                                   :: ZSEUILSCATW
    real                                   :: scatwUsedForQcThresh
    real                                   :: scatwObsFGaveraged
    real                                   :: scatwUsedForQC
    logical                                :: FULLREJCT
    logical                                :: surfTypeIsSea 
    logical                                :: cldPredMissing
    logical, save                          :: firstCall = .true.


    testIndex = 13
    if (present(skipTestArr_opt)) then
      if (skipTestArr_opt(testIndex)) then
        if (firstCall) then
          firstCall = .false.
          write(*,*) 'amsubTest13BennartzScatteringIndexCheck: skipping this test'
        end if
        return
      end if
    end if

    ZSEUILSCATICE = 40.0
    ZSEUILSCATW   = 15.0
    ZSEUILSCATL   =  0.0
    do nDataIndex=1,numObsToProcess
      FULLREJCT = .FALSE.
      surfTypeIsSea = .false.

      if (landQualifierIndice (nDataIndex) == 1) then
        if ( modelInterpSeaIce (nDataIndex) > 0.01 ) then ! sea ice 
          if (scatIndexOverWaterObs(nDataIndex) /= mwbg_realMissing .and. &
              scatIndexOverWaterObs(nDataIndex) > ZSEUILSCATICE) then
            FULLREJCT = .TRUE.
          end if

        else                                    ! sea 
          surfTypeIsSea = .true.

          if (tvs_mwAllskyAssim) then
            scatwObsFGaveraged = 0.5 * (scatIndexOverWaterObs(nDataIndex) + scatIndexOverWaterFG(nDataIndex))
            scatwUsedForQC = scatwObsFGaveraged
            scatwUsedForQcThresh = mwbg_maxSiOverWaterThreshold
            cldPredMissing = (scatIndexOverWaterObs(nDataIndex) == mwbg_realMissing .or. &
                              scatIndexOverWaterFG(nDataIndex) == mwbg_realMissing)
          else
            scatwUsedForQC = scatIndexOverWaterObs(nDataIndex)
            scatwUsedForQcThresh = ZSEUILSCATW
            cldPredMissing = (scatIndexOverWaterObs(nDataIndex) == mwbg_realMissing)
          end if

          if (.not. cldPredMissing .and. scatwUsedForQC > scatwUsedForQcThresh) then
            FULLREJCT = .TRUE.
          end if
        end if

      else                                      ! land
        if ( scatIndexOverLandObs(nDataIndex) /= mwbg_realMissing    .and. &
             scatIndexOverLandObs(nDataIndex) > ZSEUILSCATL ) then
          FULLREJCT = .TRUE.
        end if
      end if
      if ( FULLREJCT )  then
        do nChannelIndex=1,actualNumChannel
          qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
          obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
          obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**7)
          rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
            rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) + 1
        end do
        if (mwbg_debug) then
          write(*,*)  stnId(2:9), ' BENNARTZ scattering index check REJECT, scatIndexOverWaterObs=', &
                      scatIndexOverWaterObs(nDataIndex), ', scatIndexOverWaterFG=', scatIndexOverWaterFG(nDataIndex), &
                     ', scatIndexOverLandObs= ',scatIndexOverLandObs(nDataIndex)
        end if
      end if

      if (tvs_mwAllskyAssim .and. surfTypeIsSea) then
        scatwObsFGaveraged = 0.5 * (scatIndexOverWaterObs(nDataIndex) + scatIndexOverWaterFG(nDataIndex))

        ! In all-sky mode, turn on bit=23 for channels in chanIgnoreInAllskyGenCoeff(:)
        ! as cloud-affected radiances over sea when there is mismatch between 
        ! scatIndexOverWaterObs and scatIndexOverWaterFG (to be used in gen_bias_corr)
        if (scatwObsFGaveraged > mwbg_cloudySiThresholdBcorr .or. cldPredMissing) then
          do nChannelIndex = 1,actualNumChannel
            channelval = obsChannels2D(nChannelIndex,nDataIndex)
            chanIndex = ISRCHEQI(chanIgnoreInAllskyGenCoeff(:),channelval)
            if (chanIndex == 0) cycle
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**23)
          end do
          if ( mwbg_debug ) then
            write(*,*) stnId(2:9),' BENNARTZ scattering index check', &
                      ' cloud-affected obs. scatwObsFGaveraged= ', scatwObsFGaveraged, ', threshold= ', &
                      mwbg_cloudySiThresholdBcorr
          end if
        end if
      end if
      
      if (tvs_mwAllskyAssim .and. landQualifierIndice(nDataIndex) == 1) then
        scatwObsFGaveraged = 0.5 * (scatIndexOverWaterObs(nDataIndex) + scatIndexOverWaterFG(nDataIndex))
        cldPredMissing = (scatIndexOverWaterObs(nDataIndex) == mwbg_realMissing .or. &
                          scatIndexOverWaterFG(nDataIndex) == mwbg_realMissing)

        ! In all-sky mode, reject observations over sea if: 
        !   - scatwObsFGaveraged can not be computed.
        !   - scatwObsFGaveraged smaller than the minimum value
        !   - scatwObsFGaveraged greater than the maximum value
        ! scatwObsFGaveraged is needed to define obs error.
        if (cldPredMissing .or. scatwObsFGaveraged < mwbg_minSiOverWaterThreshold .or. &
            scatwObsFGaveraged > mwbg_maxSiOverWaterThreshold) then

          loopChannel3: do nChannelIndex = 1, actualNumChannel
            channelval = obsChannels2D(nChannelIndex,nDataIndex)
            if (oer_useStateDepSigmaObs(channelval,sensorIndex)) then
              qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
              obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
              obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**7)
              rejectionCodArray(testIndex,channelval,sensorIndex) = &
                      rejectionCodArray(testIndex,channelval,sensorIndex)+ 1
            end if
          end do loopChannel3
        end if

      end if ! if (tvs_mwAllskyAssim .and. surfTypeIsSea)
      
    end do !do nDataIndex=1,numObsToProcess

  end subroutine amsubTest13BennartzScatteringIndexCheck

  !--------------------------------------------------------------------------
  ! amsuaTest14RogueCheck
  !--------------------------------------------------------------------------
  subroutine amsuaTest14RogueCheck(obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, ROGUEFAC, &
                                   landQualifierIndice, ompTb2D, cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                                   ISFCREJ, obsFlags2D, qcIndicator)

    !:Purpose:                     14) test 14: "Rogue check" for (O-P) Tb residuals out of range.
    !                                  (single/full). Les observations, dont le residu (O-P) 
    !                                  depasse par un facteur (roguefac) l'erreur totale des TOVS.
    !                                  N.B.: a reject by any of the 3 surface channels produces the 
    !                                  rejection of AMSUA-A channels 1-5 and 15.

    implicit none
    ! Arguments
    integer,     intent(in)                :: actualNumChannel                                ! nombre de canaux des observations 
    integer,     intent(in)                :: numObsToProcess                                 ! nombre de tovs    
    integer,     intent(in)                :: obsChannels2D(actualNumChannel,numObsToProcess) ! observations channels
    integer,     intent(in)                :: sensorIndex                                     ! numero de satellite (i.e. indice) 
    character *9,intent(in)                :: stnId                                           ! identificateur du satellite
    real,        intent(in)                :: ROGUEFAC(:)                                     ! rogue factor 
    integer,     intent(in)                :: landQualifierIndice(numObsToProcess)            ! land/sea identifier
    real,        intent(in)                :: cloudLiquidWaterPathObs(numObsToProcess)        ! retrieved cloud liquid water from observation
    real,        intent(in)                :: cloudLiquidWaterPathFG(numObsToProcess)         ! retrieved cloud liquid water from background
    real,        intent(in)                :: ompTb2D(actualNumChannel,numObsToProcess)       ! radiance o-p 
    integer,     intent(in)                :: ISFCREJ(:)
    integer,     intent(inout)             :: obsFlags2D(actualNumChannel,numObsToProcess)    ! marqueur de radiance 
    integer,     intent(inout)             :: qcIndicator(actualNumChannel,numObsToProcess)   ! indicateur du QC par canal
    ! Locals
    integer                                :: channelval
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex
    integer                                :: INDXCAN 
    real                                   :: XCHECKVAL
    real                                   :: clwThresh1 
    real                                   :: clwThresh2
    real                                   :: errThresh1 
    real                                   :: errThresh2
    real                                   :: sigmaObsErrUsed  
    logical                                :: SFCREJCT
    logical                                :: surfTypeIsWater
    real                                   :: clwObsFGaveraged 

    testIndex = 14
    do nDataIndex=1,numObsToProcess
      surfTypeIsWater = (landQualifierIndice(nDataIndex) == 1)
      SFCREJCT = .FALSE.
      do nChannelIndex=1,actualNumChannel
        channelval = obsChannels2D(nChannelIndex,nDataIndex)
        if ( channelval /= 20 ) then
          ! using state-dependent obs error only over water.
          ! obs over sea-ice will be rejected in test 15.
          if ( tvs_mwAllskyAssim .and. oer_useStateDepSigmaObs(channelval,sensorIndex) &
                .and. surfTypeIsWater ) then
            clwThresh1 = oer_cldPredThresh(channelval,sensorIndex,1)
            clwThresh2 = oer_cldPredThresh(channelval,sensorIndex,2)
            errThresh1 = oer_errThreshAllsky(channelval,sensorIndex,1)
            errThresh2 = oer_errThreshAllsky(channelval,sensorIndex,2)
            clwObsFGaveraged = 0.5 * (cloudLiquidWaterPathObs(nDataIndex) + cloudLiquidWaterPathFG(nDataIndex))
            if (cloudLiquidWaterPathObs(nDataIndex) == mwbg_realMissing .or. cloudLiquidWaterPathFG(nDataIndex) == mwbg_realMissing) then
              sigmaObsErrUsed = MPC_missingValue_R4
            else
              sigmaObsErrUsed = calcStateDepObsErr_r4(clwThresh1,clwThresh2,errThresh1, &
                                                        errThresh2,clwObsFGaveraged)
            end if
          else
            sigmaObsErrUsed = oer_toverrst(channelval,sensorIndex)
          end if
          ! For sigmaObsErrUsed=MPC_missingValue_R4 (cloudLiquidWaterPathObs[FG]=mwbg_realMissing
          ! in all-sky mode), the observation is already rejected in test 12.
          XCHECKVAL = ROGUEFAC(channelval) * sigmaObsErrUsed
          if ( ompTb2D(nChannelIndex,nDataIndex) /= mwbg_realMissing .and. &
              ABS(ompTb2D(nChannelIndex,nDataIndex)) >= XCHECKVAL .and. &
              sigmaObsErrUsed /= MPC_missingValue_R4 ) then
            qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**16)
            rejectionCodArray(testIndex,channelval,sensorIndex) = &
                rejectionCodArray(testIndex,channelval,sensorIndex) + 1 
            if ( mwbg_debug ) then
              write(*,*)stnId(2:9),'ROGUE CHECK REJECT.NO.', &
                      ' OBS = ',nDataIndex, &
                      ' CHANNEL= ',channelval, &
                      ' CHECK VALUE= ',XCHECKVAL, &
                      ' TBOMP= ',ompTb2D(nChannelIndex,nDataIndex)
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
        do nChannelIndex=1,actualNumChannel
          INDXCAN = ISRCHEQI(ISFCREJ,obsChannels2D(nChannelIndex,nDataIndex))
          if ( INDXCAN /= 0 )  then
            if ( qcIndicator(nChannelIndex,nDataIndex) /= testIndex ) then
               qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
               obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
               obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**16)
               rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                         rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1
            end if
          end if
        end do
      end if

    end do
  end subroutine amsuaTest14RogueCheck

  !--------------------------------------------------------------------------
  ! amsubTest14RogueCheck
  !--------------------------------------------------------------------------
  subroutine amsubTest14RogueCheck(obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, ROGUEFAC, &
                                   terrainTypeIndice, landQualifierIndice, ompTb2D, scatIndexOverWaterObs, scatIndexOverWaterFG, &
                                   ICH2OMPREJ, obsFlags2D, qcIndicator, skipTestArr_opt)

    !:Purpose:                     14) test 14: "Rogue check" for (O-P) Tb residuals out of range.
    !                                  (single)
    !                                 Also, remove CH2,3,4,5 if CH2 |O-P|>5K           (partial) 
    !                                  

    implicit none

    ! Arguments:
    integer,     intent(in)                :: actualNumChannel                                ! nombre de canaux des observations 
    integer,     intent(in)                :: numObsToProcess                                 ! nombre de tovs    
    integer,     intent(in)                :: obsChannels2D(actualNumChannel,numObsToProcess) ! observations channels
    integer,     intent(in)                :: sensorIndex                                     ! numero de satellite (i.e. indice) 
    character *9,intent(in)                :: stnId                                           ! identificateur du satellite
    real,        intent(in)                :: ROGUEFAC(:)                                     ! rogue factor 
    integer,     intent(in)                :: landQualifierIndice(numObsToProcess)            ! land/sea identifier
    integer,     intent(in)                :: terrainTypeIndice(numObsToProcess)              ! terrain type
    real,        intent(in)                :: ompTb2D(actualNumChannel,numObsToProcess)       ! radiance o-p 
    real,        intent(in)                :: scatIndexOverWaterObs(numObsToProcess)          ! retrieved scattering-index over water from observation
    real,        intent(in)                :: scatIndexOverWaterFG(numObsToProcess)           ! retrieved scattering-index over water from background
    integer,     intent(in)                :: ICH2OMPREJ(:)
    integer,     intent(inout)             :: obsFlags2D(actualNumChannel,numObsToProcess)    ! marqueur de radiance 
    integer,     intent(inout)             :: qcIndicator(actualNumChannel,numObsToProcess)   ! indicateur du QC par canal
    logical,     intent(in), optional      :: skipTestArr_opt(:)                              ! array to set to skip the test

    ! Locals:
    integer                                :: channelval
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex
    integer                                :: INDXCAN 
    real                                   :: XCHECKVAL
    real                                   :: siThresh1 
    real                                   :: siThresh2
    real                                   :: errThresh1 
    real                                   :: errThresh2
    real                                   :: sigmaObsErrUsed
    real                                   :: scatwObsFGaveraged 
    logical                                :: CH2OMPREJCT
    logical                                :: ch2OmpRejectInAllsky
    logical                                :: channelIsAllsky
    logical                                :: surfTypeIsWater
    logical, save                          :: firstCall = .true.

    testIndex = 14
    if (present(skipTestArr_opt)) then
      if (skipTestArr_opt(testIndex)) then
        if (firstCall) then
          firstCall = .false.
          write(*,*) 'amsubTest14RogueCheck: skipping this test'
        end if
        return
      end if
    end if

    do nDataIndex=1,numObsToProcess
      surfTypeIsWater = (landQualifierIndice(nDataIndex) == 1 .and. terrainTypeIndice(nDataIndex) /= 0)
      ch2OmpRejectInAllsky = .false.
      CH2OMPREJCT = .FALSE.
      do nChannelIndex=1,actualNumChannel
        channelval = obsChannels2D(nChannelIndex,nDataIndex)
        if ( channelval /= 20 ) then
          channelIsAllsky = (tvs_mwAllskyAssim .and. &
                             oer_useStateDepSigmaObs(channelval,sensorIndex) .and. &
                             surfTypeIsWater)
          ! using state-dependent obs error only over water.
          if (channelIsAllsky) then
            siThresh1 = oer_cldPredThresh(channelval,sensorIndex,1)
            siThresh2 = oer_cldPredThresh(channelval,sensorIndex,2)
            errThresh1 = oer_errThreshAllsky(channelval,sensorIndex,1)
            errThresh2 = oer_errThreshAllsky(channelval,sensorIndex,2)
            scatwObsFGaveraged = 0.5 * (scatIndexOverWaterObs(nDataIndex) + scatIndexOverWaterFG(nDataIndex))
            if (scatIndexOverWaterObs(nDataIndex) == mwbg_realMissing .or. &
                scatIndexOverWaterFG(nDataIndex) == mwbg_realMissing) then
              sigmaObsErrUsed = MPC_missingValue_R4
            else
              sigmaObsErrUsed = calcStateDepObsErr_r4(siThresh1,siThresh2,errThresh1, &
                                                      errThresh2,scatwObsFGaveraged)
            end if
          else
            sigmaObsErrUsed = oer_toverrst(channelval,sensorIndex)
          end if
          ! For sigmaObsErrUsed=MPC_missingValue_R4 (scatIndexOverWaterObs[FG]=mwbg_realMissing
          ! in all-sky mode), the observation is already rejected in test 13.
          XCHECKVAL = ROGUEFAC(channelval) * sigmaObsErrUsed
          if (ompTb2D(nChannelIndex,nDataIndex) /= mwbg_realMissing .and. &
              abs(ompTb2D(nChannelIndex,nDataIndex)) >= XCHECKVAL .and. &
              sigmaObsErrUsed /= MPC_missingValue_R4) then
            qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**16)
            rejectionCodArray(testIndex,channelval,sensorIndex) = &
                rejectionCodArray(testIndex,channelval,sensorIndex) + 1

            ch2OmpRejectInAllSky = (channelIsAllsky .and. channelval == 44)

            if (mwbg_debug) then
              write(*,*)stnId(2:9),'ROGUE CHECK REJECT.NO.', &
                      ' OBS = ',nDataIndex, &
                      ' CHANNEL= ',channelval, &
                      ' CHECK VALUE= ',XCHECKVAL, &
                      ' TBOMP= ',ompTb2D(nChannelIndex,nDataIndex)
            end if
          end if
          if (channelval == 44 .and. &
              ompTb2D(nChannelIndex,nDataIndex) /= mwbg_realMissing) then
            if (channelIsAllsky) then
              if (ch2OmpRejectInAllSky) CH2OMPREJCT = .true.
            else
              if (abs(ompTb2D(nChannelIndex,nDataIndex)) >= 5.0) CH2OMPREJCT = .true.
            end if
          end if
        end if
      end do

      if (CH2OMPREJCT .and. landQualifierIndice(nDataIndex) == 1 .and. terrainTypeIndice(nDataIndex) /= 0) then
        do nChannelIndex=1,actualNumChannel
          INDXCAN = ISRCHEQI(ICH2OMPREJ,obsChannels2D(nChannelIndex,nDataIndex))
          if ( INDXCAN /= 0 )  then
            if ( qcIndicator(nChannelIndex,nDataIndex) /= testIndex ) then
               qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
               obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
               obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**16)
               rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                         rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1
            end if
          end if
        end do
      end if

    end do
  end subroutine amsubTest14RogueCheck


  !--------------------------------------------------------------------------
  ! amsuABTest15ChannelSelectionWithTovutil
  !--------------------------------------------------------------------------
  subroutine amsuABTest15ChannelSelectionWithTovutil(obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                                     landQualifierIndice, terrainTypeIndice, modelInterpSeaIce, &
                                                     ISFCREJ2, obsFlags2D, qcIndicator)

    !:Purpose: test 15: Channel Selection using array oer_tovutil(chan,sat)
    !          oer_tovutil = 0 (blacklisted), 1 (assmilate), 2 (assimilate over open water only)
    !          We also set QC flag bits 7 and 9 ON for channels with oer_tovutil=2
    !          over land or sea-ice and we set QC flag bits 7 and 9 ON for channels
    !          1-3,15 over land or sea-ice REGARDLESS of oer_tovutil value 
    !          (but oer_tovutil=0 always for these unassimilated channels).


    implicit none
    ! Arguments
    integer,     intent(in)                :: actualNumChannel                                ! nombre de canaux des observations 
    integer,     intent(in)                :: numObsToProcess                                 ! nombre de tovs    
    integer,     intent(in)                :: obsChannels2D(actualNumChannel,numObsToProcess) ! observations channels
    integer,     intent(in)                :: sensorIndex                                     ! numero de satellite (i.e. indice) 
    character *9,intent(in)                :: stnId                                           ! identificateur du satellite
    integer,     intent(in)                :: landQualifierIndice(numObsToProcess)            ! land sea identifier 
    integer,     intent(in)                :: terrainTypeIndice(numObsToProcess)              ! terrain type
    real  ,      intent(in)                :: modelInterpSeaIce(numObsToProcess)              ! gl
    integer,     intent(in)                :: ISFCREJ2(:)
    integer,     intent(inout)             :: obsFlags2D(actualNumChannel,numObsToProcess)    ! marqueur de radiance 
    integer,     intent(inout)             :: qcIndicator(actualNumChannel,numObsToProcess)   ! indicateur du QC par canal
    ! Locals
    integer                                :: channelval
    integer                                :: nDataIndex
    integer                                :: nChannelIndex
    integer                                :: testIndex
    integer                                :: ITRN 
    integer                                :: INDXCAN
    logical                                :: SFCREJCT

    testIndex = 15

    do nDataIndex=1,numObsToProcess
      ITRN = terrainTypeIndice(nDataIndex)
      if (landQualifierIndice (nDataIndex) == 1 .and. &
          terrainTypeIndice(nDataIndex) == -1 .and. &
          modelInterpSeaIce(nDataIndex) >= 0.01) then
        ITRN = 0
      end if        
      do nChannelIndex=1,actualNumChannel
        channelval = obsChannels2D(nChannelIndex,nDataIndex)
        INDXCAN = ISRCHEQI (ISFCREJ2,channelval)
        if ( INDXCAN /= 0 )  then
          if (landQualifierIndice (nDataIndex) == 0 .or. ITRN == 0)  then
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**7)
          end if
        end if
        if ( oer_tovutil(channelval,sensorIndex) /= 1 ) then
          SFCREJCT = .FALSE.
          if ( oer_tovutil(channelval,sensorIndex) == 0 ) then
            SFCREJCT = .TRUE.
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**11)
          else 
            if (landQualifierIndice(nDataIndex) == 0 .or. ITRN == 0)  then
              SFCREJCT = .TRUE.
              obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
              obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**7)
            end if
          end if
          if ( SFCREJCT ) then
            qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
            rejectionCodArray(testIndex,channelval,sensorIndex) = & 
                rejectionCodArray(testIndex,channelval,sensorIndex) + 1 
            if ( mwbg_debug ) then
                write(*,*)stnId(2:9),'CHANNEL REJECT: ', &
                      ' OBS = ',nDataIndex, &
                      ' CHANNEL= ',channelval
            end if
          end if
        end if
      end do
    end do

    if ( mwbg_debug ) then
       write(*,*)'qcIndicator = ',((qcIndicator(nChannelIndex,nDataIndex),nChannelIndex=1,actualNumChannel),nDataIndex=1,numObsToProcess)
    end if

  end subroutine amsuABTest15ChannelSelectionWithTovutil

  !--------------------------------------------------------------------------
  ! amsuaTest16ExcludeExtremeScattering
  !--------------------------------------------------------------------------
  subroutine amsuaTest16ExcludeExtremeScattering(obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                                 landQualifierIndice, obsTb2D, btClear2D, ompTb2D, obsFlags2D, qcIndicator)
    !:Purpose: Exclude radiances affected extreme scattering in deep convective region.
    !          For channel 5, if BT_cld-BT_clr < -0.5 OR O-BT_clr < -0.5, reject channels 4-5.

    ! Arguments
    integer,     intent(in)                :: actualNumChannel                                ! nombre de canaux des observations 
    integer,     intent(in)                :: numObsToProcess                                 ! nombre de tovs                                                 
    integer,     intent(in)                :: obsChannels2D(actualNumChannel,numObsToProcess) ! observations channels
    integer,     intent(in)                :: sensorIndex                                     ! numero de satellite (i.e. indice) 
    character *9,intent(in)                :: stnId                                           ! identificateur du satellite
    integer,     intent(in)                :: landQualifierIndice(numObsToProcess)            ! land sea identifier
    real,        intent(in)                :: obsTb2D(actualNumChannel,numObsToProcess)       ! radiance o
    real,        intent(in)                :: btClear2D(actualNumChannel,numObsToProcess)     ! clear-radiance o
    real,        intent(in)                :: ompTb2D(actualNumChannel,numObsToProcess)       ! radiance o-p 
    integer,     intent(inout)             :: obsFlags2D(actualNumChannel,numObsToProcess)    ! marqueur de radiance 
    integer,     intent(inout)             :: qcIndicator(actualNumChannel,numObsToProcess)   ! indicateur du QC par canal
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

    loopObs: do nDataIndex = 1, numObsToProcess
      surfTypeIsWater = (landQualifierIndice(nDataIndex) == 1)
      if ( .not. surfTypeIsWater ) cycle loopObs

      rejectLowPeakingChannels = .false.
      loopChannel2: do nChannelIndex = 1, actualNumChannel
        channelval = obsChannels2D(nChannelIndex,nDataIndex)
        if ( channelval /= 32 ) cycle loopChannel2

        BTcloudy = obsTb2D(nChannelIndex,nDataIndex) - ompTb2D(nChannelIndex,nDataIndex)
        simulatedCloudEffect = BTcloudy - btClear2D(nChannelIndex,nDataIndex)
        observedCloudEffect = obsTb2D(nChannelIndex,nDataIndex) - btClear2D(nChannelIndex,nDataIndex)
        if ( simulatedCloudEffect < -0.5 .or. observedCloudEffect < -0.5 ) then
          rejectLowPeakingChannels = .true.
        end if

        exit loopChannel2
      end do loopChannel2

      ! reject channel 4-5
      if ( rejectLowPeakingChannels ) then
        do nChannelIndex = 1, actualNumChannel
          channelval = obsChannels2D(nChannelIndex,nDataIndex)
          INDXCAN = ISRCHEQI(lowPeakingChannelsList,channelval)
          if ( INDXCAN /= 0 )  then
            qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**16)
            rejectionCodArray(testIndex,channelval,sensorIndex) = &
                rejectionCodArray(testIndex,channelval,sensorIndex) + 1 
          end if

          if ( mwbg_debug ) then
            write(*,*) stnId(2:9),' extreme scattering check reject: ', &
                    ' obs location index = ', nChannelIndex, &
                    ' channel = 1-5'
          end if
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
  subroutine mwbg_tovCheckAmsua(landQualifierIndice, obsChannels, obsTb, btClear, obsTbBiasCorr, &
                                ompTb, qcIndicator, actualNumChannel, numObsToProcess, sensorIndex, &
                                satScanPosition, modelInterpGroundIce, modelInterpTerrain, &
                                modelInterpSeaIce, terrainTypeIndice, satZenithAngle, &
                                obsGlobalMarker, obsFlags, newInformationFlag, &
                                cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                                scatIndexOverWaterObs, scatIndexOverWaterFG, &
                                stnId, RESETQC, obsLat)
  
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
    integer, allocatable, intent(inout)    :: obsGlobalMarker(:)          ! Marqueurs globaux  
    integer, intent(in)                    :: landQualifierIndice(:)      ! indicateur terre/mer
    integer, intent(in)                    :: satScanPosition(:)          ! position sur le "scan"
    integer, intent(in)                    :: obsChannels(:)              ! canaux des observations
    integer, intent(inout)                 :: terrainTypeIndice(:)        ! indicateur du type de terrain
    integer, intent(in)                    :: actualNumChannel            ! nombre de canaux des observations 
    integer, intent(in)                    :: numObsToProcess             ! nombre de tovs
    integer, intent(in)                    :: sensorIndex                 ! numero de satellite (i.e. indice)
    integer, intent(inout)                 :: obsFlags(:)                 ! marqueurs des radiances
    real, intent(in)                       :: obsTb(:)                    ! radiances
    real, intent(in)                       :: btClear(:)                  ! clear-sky radiances
    real, intent(in)                       :: obsTbBiasCorr(:)            ! correction aux radiances
    real, intent(in)                       :: ompTb(:)                    ! residus (o-p)
    real, intent(in)                       :: modelInterpGroundIce(:)     ! masque terre/mer du modele
    real, intent(in)                       :: modelInterpTerrain(:)       ! topographie du modele
    real, intent(in)                       :: modelInterpSeaIce(:)        ! etendue de glace du modele
    real, intent(in)                       :: satZenithAngle(:)           ! angle zenith du satellite (deg.)
    real, intent(in)                       :: obsLat(:)                   ! latitude
    character *9, intent(in)               :: stnId                       ! identificateur du satellite
    logical, intent(in)                    :: RESETQC                     ! reset du controle de qualite?
    integer, allocatable, intent(out)      :: qcIndicator(:,:)            ! indicateur controle de qualite tovs par canal 
    !                                                                       =0 ok, >0 rejet
    real, allocatable, intent(out)         :: cloudLiquidWaterPathObs(:)  ! retrieved cloud liquid water from observation 
    real, allocatable, intent(out)         :: cloudLiquidWaterPathFG(:)   ! retrieved cloud liquid water from background 
    real, allocatable, intent(out)         :: scatIndexOverWaterObs(:)    ! scattering index over water from observation
    real, allocatable, intent(out)         :: scatIndexOverWaterFG(:)     ! scattering index over water from background

    integer, allocatable, intent(out)      :: newInformationFlag(:)       !ATMS Information flag (newInformationFlag) values 
                                                                          ! (new BURP element 025174 in header)
                                                                          ! FOR AMSUA just fill with zeros
    !locals
    integer, parameter                     :: mwbg_maxScanAngleHIRS= 56 
    integer, parameter                     :: maxScanAngleAMSU= 30 
    real, parameter                        :: cloudyClwThreshold = 0.3
    real, parameter                        :: ZANGL = 117.6/maxScanAngleAMSU
    
    integer                                :: obsFlags2D(actualNumChannel,numObsToProcess)
    integer                                :: obsChannels2D(actualNumChannel,numObsToProcess)
    real                                   :: obsTb2D(actualNumChannel,numObsToProcess)
    real                                   :: btClear2D(actualNumChannel,numObsToProcess)
    real                                   :: obsTbBiasCorr2D(actualNumChannel,numObsToProcess)
    real                                   :: ompTb2D(actualNumChannel,numObsToProcess)
    integer, allocatable                   :: KCHKPRF(:)
    integer                                :: JI
    integer                                :: JJ
    integer                                :: INDX
    integer                                :: ICLWREJ (6)
    integer                                :: ISFCREJ (6)
    integer                                :: ISFCREJ2(4)
    integer                                :: ISCATREJ(7)
    real                                   :: EPSILON
    real, allocatable                      :: GROSSMIN(:)
    real, allocatable                      :: GROSSMAX(:) 
    real, allocatable                      :: ROGUEFAC(:)
    real                                   :: tb23 (numObsToProcess)
    real                                   :: tb31 (numObsToProcess)
    real                                   :: tb50 (numObsToProcess)
    real                                   :: tb53 (numObsToProcess)
    real                                   :: tb89 (numObsToProcess)
    real                                   :: tb23FG (numObsToProcess)
    real                                   :: tb31FG (numObsToProcess)
    real                                   :: tb50FG (numObsToProcess)
    real                                   :: tb53FG (numObsToProcess)
    real                                   :: tb89FG (numObsToProcess)
    real                                   :: ice  (numObsToProcess)
    real                                   :: tpw  (numObsToProcess)
    real                                   :: scatIndexOverLandObs(numObsToProcess)
    integer                                :: err (numObsToProcess)
    integer                                :: rain(numObsToProcess)
    integer                                :: snow(numObsToProcess)
    integer                                :: channelForTopoFilter(2)
    real                                   :: altitudeForTopoFilter(2)
    logical, save                          :: LLFIRST = .true.

    EPSILON = 0.01

    call utl_reAllocate(ROGUEFAC, actualNumChannel+tvs_channelOffset(sensorIndex))
    ROGUEFAC(:) =(/ 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
                     4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
                     4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 2.0, 2.0, 2.0, &
                     3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
                     4.0, 2.0/)
    ICLWREJ(:) = (/ 28, 29, 30, 31, 32, 42 /)
    ISFCREJ(:) = (/ 28, 29, 30, 31, 32, 42 /)
    ISCATREJ(:) = (/ 28, 29, 30, 31, 32, 33, 42 /)
    ISFCREJ2(:) = (/ 28, 29, 30, 42 /)

    call utl_reAllocate(GROSSMIN, actualNumChannel+tvs_channelOffset(sensorIndex))
    GROSSMIN(:) = (/ 200., 190., 190., 180., 180., 180., 170., &
                    170., 180., 170., 170., 170., 180., 180., &
                    180., 180., 170., 180., 180., 000., 120., &
                    190., 180., 180., 180., 190., 200., 120., &
                    120., 160., 190., 190., 200., 190., 180., &
                    180., 180., 180., 190., 190., 200., 130./)

    call utl_reAllocate(GROSSMAX, actualNumChannel+tvs_channelOffset(sensorIndex))
    GROSSMAX(:) = (/ 270., 250., 250., 250., 260., 280., 290., &
                    320., 300., 320., 300., 280., 320., 300., &
                    290., 280., 330., 350., 350., 000., 310., &
                    300., 250., 250., 270., 280., 290., 310., &
                    310., 310., 300., 300., 260., 250., 250., &
                    250., 260., 260., 270., 280., 290., 330./)  
    channelForTopoFilter(:) = (/ 33, 34 /)
    altitudeForTopoFilter(:) = (/ 250., 2000./)
    ! Allocation
    call utl_reAllocate(cloudLiquidWaterPathObs,  numObsToProcess)
    call utl_reAllocate(cloudLiquidWaterPathFG,   numObsToProcess)
    call utl_reAllocate(scatIndexOverWaterObs,numObsToProcess)
    call utl_reAllocate(scatIndexOverWaterFG, numObsToProcess)

    call utl_reAllocate(kchkprf, numObsToProcess)
    call utl_reAllocate(newInformationFlag, numObsToProcess)
    call utl_reAllocate(qcIndicator, actualNumChannel, numObsToProcess)

    ! copy the original input 1D array to 2D array. The 2D arrays are used in this s/r.
    call copy1Dimto2DimIntegerArray(obsChannels, actualNumChannel, numObsToProcess, obsChannels2D)
    call copy1Dimto2DimIntegerArray(obsFlags, actualNumChannel, numObsToProcess, obsFlags2D)
    call copy1Dimto2DimRealArray(obsTbBiasCorr, actualNumChannel, numObsToProcess, obsTbBiasCorr2D)
    call copy1Dimto2DimRealArray(obsTb, actualNumChannel, numObsToProcess, obsTb2D)
    call copy1Dimto2DimRealArray(btClear, actualNumChannel, numObsToProcess, btClear2D)
    call copy1Dimto2DimRealArray(ompTb, actualNumChannel, numObsToProcess, ompTb2D)

    ! Initialisation, la premiere fois seulement!
    if (LLFIRST) then
       rejectionCodArray(:,:,:) = 0
       LLFIRST = .FALSE.
    end if
    ! fill newInformationFlag with zeros ONLY for consistency with ATMS
    newInformationFlag(:) = 0
    qcIndicator(:,:) = 0
    if ( RESETQC ) obsFlags2D(:,:) = 0

    ! Grody parameters are   extract required channels:
    call extractParamForGrodyRun (obsChannels2D, obsTb2D, ompTb2D, obsTbBiasCorr2D, numObsToProcess, actualNumChannel, &
                                     tb23,   tb31,   tb50,   tb53,   tb89, &
                                     tb23FG, tb31FG, tb50FG, tb53FG, tb89FG)
    
    !  Run Grody AMSU-A algorithms.
    call grody (err, numObsToProcess, tb23, tb31, tb50, tb53, tb89, tb23FG, tb31FG, &
                satZenithAngle, obsLat, landQualifierIndice, ice, tpw, &
                cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                rain, snow, scatIndexOverLandObs, scatIndexOverWaterObs)   

    ! 10) test 10: RTTOV reject check (single)
    ! Rejected datum flag has bit #9 on.
    call amsuABTest10RttovRejectCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, RESETQC, &
                                       stnId, obsFlags2D, qcIndicator)

    ! 1) test 1: Topography check (partial)
    ! Channel 6 is rejected for topography >  250m.
    ! Channel 7 is rejected for topography > 2000m.
    call amsuABTest1TopographyCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                     modelInterpTerrain, channelForTopoFilter, altitudeForTopoFilter, obsFlags2D, &
                                     qcIndicator)
 
    ! 2) test 2: "Land/sea qualifier" code check (full)
    ! allowed values are: 0 land, 1 sea, 2 coast.
    call amsuABTest2LandSeaQualifierCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                           landQualifierIndice, obsFlags2D, qcIndicator)

    ! 3) test 3: "Terrain type" code check (full)
    ! allowed values are: -1 missing, 0 sea-ice, 1 snow on land.
    call amsuABTest3TerrainTypeCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                      terrainTypeIndice, obsFlags2D, qcIndicator)
 
    ! 4) test 4: Field of view number check (full)
    ! Field of view acceptable range is [1,maxScanAngleAMSU]  for AMSU footprints.
    call amsuABTest4FieldOfViewCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                      satScanPosition, maxScanAngleAMSU, obsFlags2D, qcIndicator)

    ! 5) test 5: Satellite zenith angle check (full)
    ! Satellite zenith angle acceptable range is [0.,60.].
    call amsuABTest5ZenithAngleCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                      satZenithAngle, obsFlags2D, qcIndicator)

    ! 6) test 6: "Sat. zenith angle"/"field of view" consistency check.  (full)
    ! Acceptable difference between "Satellite zenith angle"  and
    ! "approximate angle computed from field of view number" is 1.8 degrees.
    call amsuABTest6ZenAngleAndFovConsistencyCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, &
                                                    stnId, satZenithAngle, ZANGL, satScanPosition, &
                                                    maxScanAngleAMSU, obsFlags2D, qcIndicator) 

    ! 7) test 7: "Land/sea qual."/"model land/sea" consistency check.    (full)
    ! Acceptable conditions are:
    !       a) both over ocean (landQualifierIndice=1; mg<0.01), new threshold 0.20, jh dec 2000,
    !       b) both over land  (landQualifierIndice=0; mg>0.80), new threshold 0.50, jh dec 2000.
    ! Other conditions are unacceptable.
    call amsuABTest7landSeaQualifyerAndModelLandSeaConsistencyCheck (obsChannels2D, sensorIndex, actualNumChannel, &
                                                                     numObsToProcess, stnId, modelInterpGroundIce, &
                                                                     landQualifierIndice, obsFlags2D, qcIndicator)

    ! 8) test 8: "Terrain type"/"Land/sea qual."/"model ice" consistency check. (full)
    ! Unacceptable conditions are:
    !        a) terrain is sea-ice and model has no ice(terrainTypeIndice=0; gl<0.01).
    !        b) terrain is sea-ice and land/sea qualifier is land (terrainTypeIndice=0; landQualifierIndice=0).
    !        c) terrain is snow on land and land/sea qualifier is sea (terrainTypeIndice=1; landQualifierIndice=1).
    !        d) terrain is missing, land/sea qualifier is sea and model has ice(terrainTypeIndice=-1; landQualifierIndice=1; gl>0.01). (enleve jh, jan 2001)
    ! NOT doNE ANYMORE 
    
    ! 9) test 9: Uncorrected Tb check (single)
    ! Uncorrected datum (flag bit #6 off). In this case switch bit 11 ON.
    call amsuABTest9UncorrectedTbCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                        RESETQC, obsFlags2D, qcIndicator) 

    ! 11) test 11: Radiance observation "Gross" check (single) 
    !  Change this test from full to single. jh nov 2000.
    call amsuABTest11RadianceGrossValueCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                              obsTb2D, GROSSMIN, GROSSMAX, obsFlags2D, qcIndicator)

    ! 12) test 12: Grody cloud liquid water check (partial)
    ! For Cloud Liquid Water > clwQcThreshold, reject AMSUA-A channels 1-5 and 15.
    call amsuaTest12GrodyClwCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                   cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, landQualifierIndice, &
                                   ICLWREJ, obsFlags2D, qcIndicator)

    ! 13) test 13: Grody scattering index check (partial)
    ! For Scattering Index > 9, reject AMSUA-A channels 1-6 and 15.
    call amsuaTest13GrodyScatteringIndexCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                               scatIndexOverWaterObs, landQualifierIndice, terrainTypeIndice, &
                                               ISCATREJ, obsFlags2D, qcIndicator)

    ! 14) test 14: "Rogue check" for (O-P) Tb residuals out of range. (single/full)
    ! Les observations, dont le residu (O-P) depasse par un facteur (roguefac) l'erreur totale des TOVS.
    ! N.B.: a reject by any of the 3 surface channels produces the rejection of AMSUA-A channels 1-5 and 15. 
    call amsuaTest14RogueCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, ROGUEFAC, &
                                landQualifierIndice, ompTb2D, cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                                ISFCREJ, obsFlags2D, qcIndicator)

    ! 15) test 15: Channel Selection using array oer_tovutil(chan,sat)
    !  oer_tovutil = 0 (blacklisted)
    !                1 (assmilate)
    !                2 (assimilate over open water only)
    !
    !  We also set QC flag bits 7 and 9 ON for channels with oer_tovutil=2 
    !  over land or sea-ice
    !    and 
    !  we set QC flag bits 7 and 9 ON for channels 1-3,15 over land
    !  or sea-ice REGARDLESS of oer_tovutil value (but oer_tovutil=0 always for
    !  these unassimilated channels).
    call amsuABTest15ChannelSelectionWithTovutil (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                                  landQualifierIndice, terrainTypeIndice, modelInterpSeaIce, &
                                                  ISFCREJ2, obsFlags2D, qcIndicator)

    ! 16) test 16: exclude radiances affected by extreme scattering in deep convective region in all-sky mode.
    call amsuaTest16ExcludeExtremeScattering(obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                             landQualifierIndice, obsTb2D, btClear2D, ompTb2D,  obsFlags2D, qcIndicator) 

    !  Synthese de la controle de qualite au niveau de chaque point
    !  d'observation. Code:
    !            =0, aucun rejet,
    !            >0, au moins un canal rejete.

    do JJ = 1, numObsToProcess
      KCHKPRF(JJ) = 0
      do JI = 1, actualNumChannel
        KCHKPRF(JJ) = MAX(KCHKPRF(JJ),qcIndicator(JI,JJ))
      end do
    end do

    if ( mwbg_debug ) then
      write(*,*)'KCHKPRF = ',(KCHKPRF(JJ),JJ=1,numObsToProcess)
    end if

    ! Copy the modified FLAG to the 1D array, used outside this s/r.
    do JJ = 1, numObsToProcess
      do JI = 1, actualNumChannel
        INDX = (JJ-1) * actualNumChannel + JI
        obsFlags(INDX) = obsFlags2D(JI,JJ)
      end do
    end do

    ! reset global marker flag (55200) and mark it if observtions are rejected
    call resetQcCases(RESETQC, KCHKPRF, obsGlobalMarker)

    !###############################################################################
    ! FINAL STEP: set terrain type to sea ice given certain conditions
    !###############################################################################
    call setTerrainTypeToSeaIce(modelInterpSeaIce, landQualifierIndice, terrainTypeIndice)

  end subroutine mwbg_tovCheckAmsua

  !--------------------------------------------------------------------------
  ! mwbg_tovCheckAmsub
  !--------------------------------------------------------------------------
  subroutine mwbg_tovCheckAmsub(landQualifierIndice, obsChannels, obsTb, btClear, obsTbBiasCorr, &
                                ompTb, qcIndicator, actualNumChannel, numObsToProcess, sensorIndex, &
                                satScanPosition, modelInterpGroundIce, modelInterpTerrain, &
                                modelInterpSeaIce, terrainTypeIndice, satZenithAngle, &
                                obsGlobalMarker, obsFlags, newInformationFlag, &
                                cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                                scatIndexOverWaterObs, scatIndexOverWaterFG, stnId, RESETQC)
  
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
    integer, allocatable, intent(inout)    :: obsGlobalMarker(:)          !Marqueurs globaux  
    integer, intent(in)                    :: landQualifierIndice(:)      ! indicateur terre/mer
    integer, intent(in)                    :: satScanPosition(:)          ! position sur le "scan"
    integer, intent(in)                    :: obsChannels(:)              ! canaux des observations
    integer, intent(inout)                 :: terrainTypeIndice(:)        ! indicateur du type de terrain
    integer, intent(in)                    :: actualNumChannel            ! nombre de canaux des observations 
    integer, intent(in)                    :: numObsToProcess             ! nombre de tovs
    integer, intent(in)                    :: sensorIndex                 ! numero de satellite (i.e. indice)
    integer, intent(inout)                 :: obsFlags(:)                 ! marqueurs des radiances
    real, intent(in)                       :: obsTb(:)                    ! radiances from observation
    real, intent(in)                       :: btClear(:)                  ! clear-sky radiances from background
    real, intent(in)                       :: obsTbBiasCorr(:)            ! correction aux radiances
    real, intent(in)                       :: ompTb(:)                    ! residus (o-p)
    real, intent(in)                       :: modelInterpGroundIce(:)     ! masque terre/mer du modele
    real, intent(in)                       :: modelInterpTerrain(:)       ! topographie du modele
    real, intent(in)                       :: modelInterpSeaIce(:)        ! etendue de glace du modele
    real, intent(in)                       :: satZenithAngle(:)           ! angle zenith du satellite (deg.)
    character *9, intent(in)               :: stnId                       ! identificateur du satellite
    logical, intent(in)                    :: RESETQC                     ! reset du controle de qualite?
    integer, allocatable, intent(out)      :: qcIndicator(:,:)            ! indicateur controle de qualite tovs par canal 
    !                                                                       =0 ok, >0 rejet,
    real, allocatable, intent(out)         :: cloudLiquidWaterPathObs(:)  ! retrieved cloud liquid water from observation 
    real, allocatable, intent(out)         :: cloudLiquidWaterPathFG(:)   ! retrieved cloud liquid water from background 
    real, allocatable, intent(out)         :: scatIndexOverWaterObs(:)    ! scattering index over water from observation
    real, allocatable, intent(out)         :: scatIndexOverWaterFG(:)     ! scattering index over water from background

    integer, allocatable, intent(out)      :: newInformationFlag(:)       !ATMS Information flag (newInformationFlag) values 
                                                                          ! (new BURP element 025174 in header)
                                                                          ! FOR AMSUA just fill with zeros

    !locals
    integer, parameter                     :: mwbg_maxScanAngleHIRS= 56 
    integer, parameter                     :: maxScanAngleAMSU= 90 
    real, parameter                        :: ZANGL =  117.6/maxScanAngleAMSU
    
    integer                                :: obsFlags2D(actualNumChannel,numObsToProcess)
    integer                                :: obsChannels2D(actualNumChannel,numObsToProcess)
    real                                   :: obsTb2D(actualNumChannel,numObsToProcess)
    real                                   :: btClear2D(actualNumChannel,numObsToProcess)
    real                                   :: obsTbBiasCorr2D(actualNumChannel,numObsToProcess)
    real                                   :: ompTb2D(actualNumChannel,numObsToProcess)
    integer, allocatable                   :: KCHKPRF(:)          
    integer                                :: JI
    integer                                :: JJ
    integer                                :: INDX
    integer                                :: ISFCREJ(2)
    integer                                :: ICH2OMPREJ(4)
    integer                                :: ISFCREJ2(1)
    integer                                :: chanIgnoreInAllskyGenCoeff(5)
    real                                   :: EPSILON
    real, allocatable                      :: GROSSMIN(:)
    real, allocatable                      :: GROSSMAX(:) 
    real, allocatable                      :: ROGUEFAC(:)
    real                                   :: tb89 (numObsToProcess)
    real                                   :: tb150 (numObsToProcess)
    real                                   :: tb1831 (numObsToProcess)
    real                                   :: tb1832 (numObsToProcess)
    real                                   :: tb1833 (numObsToProcess)
    real                                   :: tb89FG (numObsToProcess)
    real                                   :: tb150FG(numObsToProcess)
    real                                   :: tb89FgClear(numObsToProcess)
    real                                   :: tb150FgClear(numObsToProcess)    
    real                                   :: scatIndexOverLandObs(numObsToProcess)
    integer                                :: err (numObsToProcess)
    integer                                :: channelForTopoFilter(3)
    real                                   :: altitudeForTopoFilter(3)
    logical, save                          :: LLFIRST = .true.

    EPSILON = 0.01

    call utl_reAllocate(ROGUEFAC, actualNumChannel+tvs_channelOffset(sensorIndex))
    ROGUEFAC(:) =(/ 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
                    4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
                    4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 2.0, 2.0, 2.0, &
                    3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, &
                    4.0, 2.0, 2.0, 2.0, 4.0, 4.0, 4.0/)

    ISFCREJ(:) = (/ 43, 44 /)
    ISFCREJ2(:) = (/ 43 /)
    ICH2OMPREJ(:) = (/ 44, 45, 46, 47 /)
    call utl_reAllocate(GROSSMIN, actualNumChannel+tvs_channelOffset(sensorIndex))
    GROSSMIN(:) = (/ 200., 190., 190., 180., 180., 180., 170., &
                    170., 180., 170., 170., 170., 180., 180., &
                    180., 180., 170., 180., 180., 000., 120., &
                    190., 180., 180., 180., 190., 200., 120., &
                    120., 160., 190., 190., 200., 190., 180., &
                    180., 180., 180., 190., 190., 200., 130., &
                    130., 130., 130., 130., 130./)
    call utl_reAllocate(GROSSMAX, actualNumChannel+tvs_channelOffset(sensorIndex))
    GROSSMAX(:) = (/ 270., 250., 250., 250., 260., 280., 290., &
                    320., 300., 320., 300., 280., 320., 300., &
                    290., 280., 330., 350., 350., 000., 310., &
                    300., 250., 250., 270., 280., 290., 310., &
                    310., 310., 300., 300., 260., 250., 250., &
                    250., 260., 260., 270., 280., 290., 330., &
                    330., 330., 330., 330., 330./)  

    channelForTopoFilter(:) = (/ 45, 46, 47 /)
    altitudeForTopoFilter(:) = (/ 2500., 2000., 1000./)

    ! Channels excluded from genCoeff in all-sky mode
    chanIgnoreInAllskyGenCoeff(:) = (/43, 44, 45, 46, 47/)

    ! Allocation
    call utl_reAllocate(scatIndexOverWaterObs, numObsToProcess)
    call utl_reAllocate(scatIndexOverWaterFG,  numObsToProcess)
    call utl_reAllocate(cloudLiquidWaterPathObs,   numObsToProcess)
    call utl_reAllocate(cloudLiquidWaterPathFG,    numObsToProcess)

    call utl_reAllocate(kchkprf, numObsToProcess)
    call utl_reAllocate(newInformationFlag, numObsToProcess)
    call utl_reAllocate(qcIndicator, actualNumChannel, numObsToProcess)

    ! copy the original input 1D array to 2D array. The 2D arrays are used in this s/r.
    call copy1Dimto2DimIntegerArray(obsChannels, actualNumChannel, numObsToProcess, obsChannels2D)
    call copy1Dimto2DimIntegerArray(obsFlags, actualNumChannel, numObsToProcess, obsFlags2D)
    call copy1Dimto2DimRealArray(obsTbBiasCorr, actualNumChannel, numObsToProcess, obsTbBiasCorr2D)
    call copy1Dimto2DimRealArray(obsTb, actualNumChannel, numObsToProcess, obsTb2D)
    call copy1Dimto2DimRealArray(btClear, actualNumChannel, numObsToProcess, btClear2D)
    call copy1Dimto2DimRealArray(ompTb, actualNumChannel, numObsToProcess, ompTb2D)

    ! Initialisation, la premiere fois seulement!
    if (LLFIRST) then
      rejectionCodArray(:,:,:) = 0
      LLFIRST = .FALSE.
    end if
    ! fill newInformationFlag with zeros ONLY for consistency with ATMS
    newInformationFlag(:) = 0
    qcIndicator(:,:) = 0
    if ( RESETQC ) obsFlags2D(:,:) = 0

    ! Bennartz parameters are   extract required channels:
    call extractParamForBennartzRun (obsChannels2D, obsTb2D, btClear2D, ompTb2D, obsTbBiasCorr2D, numObsToProcess, actualNumChannel, &
                                     tb89, tb150, tb1831, tb1832, tb1833, &
                                     tb89FG, tb150FG, tb89FgClear, tb150FgClear)
    
    !  Run Bennartz AMSU-B algorithms.
    call bennartz (err, numObsToProcess, tb89, tb150, tb89FG, tb150FG, tb89FgClear, tb150FgClear, &
                   satZenithAngle, landQualifierIndice, scatIndexOverLandObs, &
                   scatIndexOverWaterObs, scatIndexOverWaterFG, &
                   cloudLiquidWaterPathObs, cloudLiquidWaterPathFG)

    ! 10) test 10: RTTOV reject check (single)
    ! Rejected datum flag has bit #9 on.
    call amsuABTest10RttovRejectCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, RESETQC, &
                                       stnId, obsFlags2D, qcIndicator)

    ! 1) test 1: Topography check (partial)
    ! Channel 3- 45 is rejected for topography >  2500m.
    ! Channel 4 - 46 is rejected for topography > 2000m.
    ! Channel 5 - 47 is rejected for topography > 1000m.
    call amsuABTest1TopographyCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                     modelInterpTerrain, channelForTopoFilter, altitudeForTopoFilter, obsFlags2D, &
                                     qcIndicator)
 
    ! 2) test 2: "Land/sea qualifier" code check (full)
    ! allowed values are: 0, land,
    !                     1, sea,
    !                     2, coast.
    call amsuABTest2LandSeaQualifierCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                           landQualifierIndice, obsFlags2D, qcIndicator)

    ! 3) test 3: "Terrain type" code check (full)
    ! allowed values are: -1 missing, 0 sea-ice, 1 snow on land.
    call amsuABTest3TerrainTypeCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                      terrainTypeIndice, obsFlags2D, qcIndicator)
 
    ! 4) test 4: Field of view number check (full)
    !
    ! Field of view acceptable range is [1,maxScanAngleAMSU]  for AMSU footprints.
    call amsuABTest4FieldOfViewCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                      satScanPosition, maxScanAngleAMSU, obsFlags2D, qcIndicator)

    ! 5) test 5: Satellite zenith angle check (full)
    ! Satellite zenith angle acceptable range is [0.,60.].
    call amsuABTest5ZenithAngleCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                      satZenithAngle, obsFlags2D, qcIndicator)

    ! 6) test 6: "Sat. zenith angle"/"field of view" consistency check.  (full)
    ! Acceptable difference between "Satellite zenith angle"  and
    ! "approximate angle computed from field of view number" is 1.8 degrees.
    call amsuABTest6ZenAngleAndFovConsistencyCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, &
                                                    stnId, satZenithAngle, ZANGL, satScanPosition, &
                                                    maxScanAngleAMSU, obsFlags2D, qcIndicator) 

    ! 7) test 7: "Land/sea qual."/"model land/sea" consistency check.    (full)
    ! Acceptable conditions are:
    !       a) both over ocean (landQualifierIndice=1; mg<0.01), new threshold 0.20, jh dec 2000,
    !       b) both over land  (landQualifierIndice=0; mg>0.80), new threshold 0.50, jh dec 2000.
    ! Other conditions are unacceptable.
    call amsuABTest7landSeaQualifyerAndModelLandSeaConsistencyCheck (obsChannels2D, sensorIndex, actualNumChannel, &
                                                                     numObsToProcess, stnId, modelInterpGroundIce, &
                                                                     landQualifierIndice, obsFlags2D, qcIndicator)

    ! 8) test 8: "Terrain type"/"Land/sea qual."/"model ice" consistency check. (full)
    ! Unacceptable conditions are:
    !        a) terrain is sea-ice and model has no ice(terrainTypeIndice=0; gl<0.01).
    !        b) terrain is sea-ice and land/sea qualifier is land (terrainTypeIndice=0; landQualifierIndice=0).
    !        c) terrain is snow on land and land/sea qualifier is sea (terrainTypeIndice=1; landQualifierIndice=1).
    !        d) terrain is missing, land/sea qualifier is sea and model has ice(terrainTypeIndice=-1; landQualifierIndice=1; gl>0.01). (enleve jh, jan 2001)
    ! NOT doNE ANYMORE 
    
    ! 9) test 9: Uncorrected Tb check (single) SKIP FOR NOW
    ! Uncorrected datum (flag bit #6 off). In this case switch bit 11 ON.
    ! call amsuABTest9UncorrectedTbCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, RESETQC, obsFlags2D, qcIndicator) 

    ! 11) test 11: Radiance observation "Gross" check (single) 
    !  Change this test from full to single. jh nov 2000.
    call amsuABTest11RadianceGrossValueCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                              obsTb2D, GROSSMIN, GROSSMAX, obsFlags2D, qcIndicator)

    ! 12) test 12:  Dryness index check 
    !The difference between channels AMSUB-3 and AMSUB-5 is used as an indicator
    !of "dryness" of the atmosphere. In extreme dry conditions, channels AMSUB-3 4 and 5
    ! are sensitive to the surface.
    ! Therefore, various thresholds are used to reject channels AMSUB-3 4 and 5
    !  over land and ice
    call amsubTest12DrynessIndexCheck (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                       tb1831, tb1833, landQualifierIndice, modelInterpSeaIce, obsFlags2D, qcIndicator, &
                                       skipTestArr_opt=skipTestArr(:))

    ! 13) test 13: Bennartz scattering index check (full)
    call amsubTest13BennartzScatteringIndexCheck(obsChannels2D, sensorIndex, numObsToProcess, actualNumChannel, stnId, &
                                                 scatIndexOverWaterObs, scatIndexOverWaterFG, scatIndexOverLandObs, &
                                                 landQualifierIndice, modelInterpSeaIce, &
                                                 obsFlags2D, qcIndicator, chanIgnoreInAllskyGenCoeff, &
                                                 skipTestArr_opt=skipTestArr(:))

    ! 14) test 14: "Rogue check" for (O-P) Tb residuals out of range. (single/full)
    ! Les observations, dont le residu (O-P) depasse par un facteur (roguefac) l'erreur totale des TOVS.
    ! N.B.: a reject by any of the 3 surface channels produces the rejection of AMSUA-A channels 1-5 and 15. 
    call amsubTest14RogueCheck(obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, ROGUEFAC, &
                               terrainTypeIndice, landQualifierIndice, ompTb2D, scatIndexOverWaterObs, scatIndexOverWaterFG, &
                               ICH2OMPREJ, obsFlags2D, qcIndicator, skipTestArr_opt=skipTestArr(:))

    ! 15) test 15: Channel Selection using array oer_tovutil(chan,sat)
    !  oer_tovutil = 0 (blacklisted)
    !                1 (assmilate)
    !                2 (assimilate over open water only)
    !
    !  We also set QC flag bits 7 and 9 ON for channels with oer_tovutil=2 
    !  over land or sea-ice
    !    and 
    !  we set QC flag bits 7 and 9 ON for channels 1-3,15 over land
    !  or sea-ice REGARDLESS of oer_tovutil value (but oer_tovutil=0 always for
    !  these unassimilated channels).

    call amsuABTest15ChannelSelectionWithTovutil (obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                                  landQualifierIndice, terrainTypeIndice, modelInterpSeaIce, &
                                                  ISFCREJ2, obsFlags2D, qcIndicator)

    !  Synthese de la controle de qualite au niveau de chaque point
    !  d'observation. Code:
    !            =0, aucun rejet,
    !            >0, au moins un canal rejete.

    do JJ = 1, numObsToProcess
      KCHKPRF(JJ) = 0
      do JI = 1, actualNumChannel
        KCHKPRF(JJ) = MAX(KCHKPRF(JJ),qcIndicator(JI,JJ))
      end do
    end do

    if ( mwbg_debug ) then
      write(*,*)'KCHKPRF = ',(KCHKPRF(JJ),JJ=1,numObsToProcess)
    end if

    ! Copy the modified FLAG to the 1D array, used outside this s/r.
    do JJ = 1, numObsToProcess
      do JI = 1, actualNumChannel
        INDX = (JJ-1) * actualNumChannel + JI
        obsFlags(INDX) = obsFlags2D(JI,JJ)
      end do
    end do

    ! reset global marker flag (55200) and mark it if observtions are rejected
    call resetQcCases(RESETQC, KCHKPRF, obsGlobalMarker)

    !###############################################################################
    ! FINAL STEP: set terrain type to sea ice given certain conditions
    !###############################################################################
    call setTerrainTypeToSeaIce(modelInterpSeaIce, landQualifierIndice, terrainTypeIndice)

  end subroutine mwbg_tovCheckAmsub

  !--------------------------------------------------------------------------
  ! mwbg_qcStats
  !--------------------------------------------------------------------------
  subroutine mwbg_qcStats(instName, qcIndicator, obsChannels, sensorIndex, &
                          actualNumChannel, numObsToProcess, satelliteId, &
                          LDprint)

    !:Purpose:          Cumuler ou imprimer des statistiques decriptives
    !                   des rejets tovs.
    implicit none 
    !Arguments:
    character(*), intent(in)               :: instName                                        ! Instrument Name
    integer, intent(in)                    :: qcIndicator(:,:)                                ! indicateur controle de qualite tovs par canal 
    !                                                                                           =0 ok, >0 rejet,
    integer, intent(in)                    :: actualNumChannel                                ! nombre de canaux des observations 
    integer, intent(in)                    :: numObsToProcess                                 ! nombre de tovs
    integer, intent(in)                    :: obsChannels(actualNumChannel*numObsToProcess)   ! canaux des observations
    integer, intent(in)                    :: sensorIndex                                     ! numero d'identificateur du satellite
    character(len=15), intent(in)          :: satelliteId(:)                                  ! identificateur du satellite
    logical, intent(in)                    :: LDprint                                         ! mode: imprimer ou cumuler?
    !Locals
    integer                                :: numSats
    integer                                :: JI
    integer                                :: JJ
    integer                                :: JK
    integer                                :: INTOTOBS
    integer                                :: INTOTACC
    integer, allocatable, save             :: INTOT(:)                                        ! INTOT(mwbg_maxNumSat)
    integer, allocatable, save             :: INTOTRJF(:)                                     ! INTOTRJF(mwbg_maxNumSat)
    integer, allocatable, save             :: INTOTRJP(:)                                     ! INTOTRJP(mwbg_maxNumSat)
    integer                                :: obsChannels2D(actualNumChannel,numObsToProcess) ! canaux des observations

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
      call copy1Dimto2DimIntegerArray(obsChannels, actualNumChannel, numObsToProcess, obsChannels2D)
      ! Accumulate statistics on rejects
      do JJ = 1, numObsToProcess

        INTOT(sensorIndex) = INTOT(sensorIndex) + 1
        ! Fully accepted, fully rejected or partially rejected?
        FULLREJCT = .true.
        FULLACCPT = .true.
        if (instName == "AMSUA") then 
          do JI = 1, actualNumChannel
            if ( obsChannels2D(JI,JJ) /= 20 ) then
              if ( qcIndicator(JI,JJ) /= 0 ) then
                FULLACCPT = .false.
              else
                FULLREJCT = .false.
              end if
            end if
          end do
          if ( FULLREJCT ) then
            INTOTRJF(sensorIndex) = INTOTRJF(sensorIndex) + 1
          end if
          if ( .not. FULLREJCT .and. .not.FULLACCPT ) then
            INTOTRJP(sensorIndex) = INTOTRJP(sensorIndex) + 1
          end if
        else if  (instName == "ATMS") then 
          do JI = 1, actualNumChannel
            if ( qcIndicator(JI,JJ) /= 0 ) then
              FULLACCPT = .false.
            else
              FULLREJCT = .false.
            end if
          end do
          if ( FULLREJCT ) then
            INTOTRJF(sensorIndex) = INTOTRJF(sensorIndex) + 1
          end if
          if ( .not. FULLREJCT .and. .not.FULLACCPT ) then
            INTOTRJP(sensorIndex) = INTOTRJP(sensorIndex) + 1
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
          do JJ = 1, actualNumChannel 
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
          do JJ = 1, actualNumChannel 
            write(*,'(3X,I2,t10,"|",5I7)') JJ,(rejectionCodArray(JI,JJ,JK), &
                                        JI=1,mwbg_maxNumTest)
          end do
          write(*,'(1x,59("-"))')
          write(*,'(//,1x,59("-"))')
          write(*,'(t10,"|",t19,"2. QC2 REJECT CATEGORIES")')
          write(*,'(" CHANNEL",t10,"|",50("-"))') 
          write(*,'(t10,"|",5i7)') (JI,JI=1,mwbg_maxNumTest)
          write(*,'(1x,"--------|",50("-"))')
          do JJ = 1, actualNumChannel
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
  subroutine resetQcCases(RESETQC, KCHKPRF, obsGlobalMarker)
    !:Purpose:        allumer la bit (6) indiquant que l'observation a un element
    !                 rejete par le controle de qualite de l'AO.
    !                 N.B.: si on est en mode resetqc, on remet le marqueur global a
    !                 sa valeur de defaut, soit 1024,  avant de faire la mise a jour.
    implicit none
    !Arguments:
    logical,              intent(in)     :: RESETQC               !reset the quality control flags before adding the new ones ?
    integer,              intent(in)     :: KCHKPRF(:)            !indicateur global controle de qualite tovs. Code:
    integer,              intent(inout)  :: obsGlobalMarker(:)    !Marqueurs globaux  

    !Locals
    integer                              :: dataNum 
    integer                              :: dataIndex

    dataNum = size(obsGlobalMarker)
    do dataIndex = 1, dataNum
      if (RESETQC) then
        obsGlobalMarker(dataIndex) = 1024  
      end if
      if ( KCHKPRF(dataIndex) /= 0  ) then
        obsGlobalMarker(dataIndex) = OR (obsGlobalMarker(dataIndex),2**6)
      end if
    end do
    if (mwbg_debug) then
      write(*,*) ' KCHKPRF   = ', (KCHKPRF(dataIndex),dataIndex=1,dataNum)
      write(*,*) ' NEW FLAGS = ', (obsGlobalMarker(dataIndex),dataIndex=1,dataNum)
    end if

  end  subroutine resetQcCases

  !--------------------------------------------------------------------------
  ! setTerrainTypeToSeaIce
  !--------------------------------------------------------------------------
  subroutine setTerrainTypeToSeaIce(modelInterpSeaIce, landQualifierIndice, terrainTypeIndice)

    !:Purpose:       Dans les conditions suivantes:
    !                1) l'indicateur terre/mer indique l'ocean (landQualifierIndice=1),
    !                2) le "terrain type" est manquant (terrainTypeIndice=-1),
    !                3) le modele indique de la glace (gl >= 0.01),
    !                on specifie "sea ice" pour le "terrain type" (terrainTypeIndice=0).
    implicit none 
    !Arguments:
    real,                 intent(in)     :: modelInterpSeaIce(:)    ! sea ice
    integer,              intent(in)     :: landQualifierIndice(:)  ! land sea qualifier
    integer,              intent(inout)  :: terrainTypeIndice(:)    ! terrain type

    !Locals
    integer                              :: dataNum 
    integer                              :: dataIndex

    dataNum = size(terrainTypeIndice)

    if ( mwbg_debug ) then
      write(*,*) ' OLD TERRAIN type = ', (terrainTypeIndice(dataIndex),dataIndex=1,dataNum)
      write(*,*) ' landQualifierIndice = ', (landQualifierIndice(dataIndex),dataIndex=1,dataNum)
      write(*,*) ' modelInterpSeaIce = ', (modelInterpSeaIce(dataIndex),dataIndex=1,dataNum)
    end if
    do dataIndex = 1, dataNum
      if (landQualifierIndice (dataIndex) == 1 .and. terrainTypeIndice(dataIndex) == -1 .and. &
          modelInterpSeaIce (dataIndex) >= 0.01) then
        terrainTypeIndice(dataIndex) = 0
      end if
    end do
    if ( mwbg_debug ) then
      write(*,*) ' NEW TERRAIN type = ', (terrainTypeIndice(dataIndex),dataIndex=1,dataNum)
    end if
    
  end  subroutine setTerrainTypeToSeaIce

  !--------------------------------------------------------------------------
  ! GRODY
  !--------------------------------------------------------------------------
  subroutine GRODY (ier, numObsToProcess, tb23, tb31, tb50, tb53, tb89, tb23FG, tb31FG, &
                    satZenithAngle, obsLat, landQualifierIndice, ice, tpw, &
                    cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                    rain, snow, scatIndexOverLandObs, scatIndexOverWaterObs)
    !:Purpose: Compute the following parameters using 5 AMSU-A channels:
    !          - sea ice, 
    !          - total precipitable water, 
    !          - cloud liquid water, 
    !          - ocean/land rain, 
    !          - snow cover/glacial ice,
    !          - scattering index (sur la terre et sur l'eau).
    !          The four channels used are: 23Ghz, 31Ghz, 50Ghz and 89Ghz.
    !          REGERENCES     N. Grody, NOAA/NESDIS, ....

    implicit none

    ! Arguments:
    integer, intent(out) :: ier(:)                      ! error return code: 0 ok, 1 input parameter out of range.
    integer,  intent(in) :: numObsToProcess             ! number of points to process
    real,     intent(in) :: tb23  (:)                   ! 23Ghz brightness temperature (K)
    real,     intent(in) :: tb31  (:)                   ! 31Ghz brightness temperature (K)
    real,     intent(in) :: tb50  (:)                   ! 50Ghz brightness temperature (K)
    real,     intent(in) :: tb53  (:)                   ! 53Ghz brightness temperature (K)
    real,     intent(in) :: tb89  (:)                   ! 89Ghz brightness temperature (K)
    real,     intent(in) :: tb23FG(:)                   ! 23Ghz brightness temperature from background (K)
    real,     intent(in) :: tb31FG(:)                   ! 31Ghz brightness temperature from background (K)
    real,     intent(in) :: satZenithAngle(:)           ! satellite zenith angle (deg.)
    real,     intent(in) :: obsLat(:)                   ! latitude (deg.)
    integer,  intent(in) :: landQualifierIndice(:)      ! land/sea indicator (0=land;1=ocean)
    real,    intent(out) :: ice (:)                     ! sea ice concentration (0-100%)
    real,    intent(out) :: tpw (:)                     ! total precipitable water (0-70mm)
    real,    intent(out) :: cloudLiquidWaterPathObs(:)  ! retrieved cloud liquid water from observation (0-3mm)
    real,    intent(out) :: cloudLiquidWaterPathFG (:)  ! retrieved cloud liquid water from background (0-3mm)
    integer, intent(out) :: rain(:)                     ! rain identification (0=no rain; 1=rain)
    integer, intent(out) :: snow(:)                     ! snow cover and glacial ice identification:
                                                        ! (0=no snow; 1=snow; 2=glacial ice)
    real,    intent(out) :: scatIndexOverLandObs (:)    ! scattering index over land
    real,    intent(out) :: scatIndexOverWaterObs (:)   ! scattering index over water

    ! Locals:
    real :: siw, sil, df1, df2, df3, a, b, c, d, e23
    real :: ei, cosz, tt, scat, sc31, abslat, t23, t31, t50, t89
    real :: sc50, par, t53
    real :: dif285t23, dif285t31, epsilon
    real :: dif285t23FG, dif285t31FG
    integer :: i

    data epsilon / 1.E-30 /

    logical skipLoopChan15Missing 
    !
    ! Notes: In the case where an output parameter cannot be calculated, the
    !        value of this parameter is to to the missing value, i.e. -99.

    ! 1) Initialise output parameters:
    do i = 1, numObsToProcess
      ice  (i) = mwbg_realMissing
      tpw  (i) = mwbg_realMissing
      cloudLiquidWaterPathObs(i) = mwbg_realMissing
      cloudLiquidWaterPathFG(i) = mwbg_realMissing
      scatIndexOverLandObs(i) = mwbg_realMissing
      scatIndexOverWaterObs(i) = mwbg_realMissing
      rain (i) = nint(mwbg_realMissing)
      snow (i) = nint(mwbg_realMissing)
    end do

    ! 2) Validate input parameters:
    do i = 1, numObsToProcess
      if ( tb23(i)    < 120.  .or. &
           tb23(i)    > 350.  .or. &
           tb31(i)    < 120.  .or. &
           tb31(i)    > 350.  .or. &
           tb50(i)    < 120.  .or. &
           tb50(i)    > 350.  .or. &
           tb53(i)    < 120.  .or. &
           tb53(i)    > 350.  .or. &
           satZenithAngle(i)   < -90.  .or. &
           satZenithAngle(i)   > 90.   .or. &
           obsLat(i)  < -90.  .or. &
           obsLat(i)  > 90.   .or. &
           landQualifierIndice(i) < 0 .or. &
           landQualifierIndice(i) > 1 ) then
        ier(i) = 1
      else
        ier(i) = 0
      end if
    end do

    !3) Compute parameters:
    loopObsGrody: do i = 1, numObsToProcess
      if ( ier(i) == 0 ) then
        abslat = abs(obsLat(i))
        cosz   = cosd(satZenithAngle(i))
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

          scatIndexOverLandObs (i) = sil
          scatIndexOverWaterObs (i) = siw
        end if

        ! discriminate functions:
        df1 =  2.85 + 0.020*t23 - 0.028*t50 ! used to identify (also remove) sea ice
        df2 =  5.10 + 0.078*t23 - 0.096*t50 ! used to identify (also remove) warm deserts
        df3 = 10.20 + 0.036*t23 - 0.074*t50 ! used to identify (also remove) cold deserts

        if ( landQualifierIndice(i)== 1 ) then

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
            tpw(i) = mwbg_realMissing
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

          !3.3) Cloud liquid water from obs (cloudLiquidWaterPathObs) and background state (cloudLiquidWaterPathFG):
          ! identify and remove sea ice
          if ( abslat > 50.  .and. &
              df1     >  0.0        ) then  
            cloudLiquidWaterPathObs(i) = mwbg_realMissing
            cloudLiquidWaterPathFG(i) = mwbg_realMissing
          else
            a =  8.240 - (2.622 - 1.846*cosz)*cosz
            b =  0.754
            c = -2.265
            cloudLiquidWaterPathObs(i) = a + b*log(dif285t23) & 
                      + c*log(dif285t31)
            cloudLiquidWaterPathObs(i) = cloudLiquidWaterPathObs(i)*cosz           ! theoretical cloud liquid water (0-3mm)
            cloudLiquidWaterPathObs(i) = cloudLiquidWaterPathObs(i) - 0.03         ! corrected   cloud liquid water 
            cloudLiquidWaterPathObs(i) = min(3.,max(0.,cloudLiquidWaterPathObs(i)))   ! jh       

            cloudLiquidWaterPathFG(i) = a + b*log(dif285t23FG) & 
                      + c*log(dif285t31FG)
            cloudLiquidWaterPathFG(i) = cloudLiquidWaterPathFG(i)*cosz           ! theoretical cloud liquid water (0-3mm)
            cloudLiquidWaterPathFG(i) = cloudLiquidWaterPathFG(i) - 0.03         ! corrected   cloud liquid water 
            cloudLiquidWaterPathFG(i) = min(3.,max(0.,cloudLiquidWaterPathFG(i)))   ! jh       
          end if

          if ( skipLoopChan15Missing ) cycle loopObsGrody

          !3.4) Ocean rain: 0=no rain; 1=rain.
          ! identify and remove sea ice
          if ( abslat > 50.  .and. &
              df1    >  0.0        ) then  
            rain(i) = nint(mwbg_realMissing)
          else                                   ! remove non-precipitating clouds
            if ( cloudLiquidWaterPathObs(i) > 0.3 .or. &
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
        print *, 'GRODY: i,tb23(i),tb31(i),tb50(i),tb89(i),satZenithAngle(i),obsLat(i), &
                  landQualifierIndice(i) = ', &
                  i,tb23(i),tb31(i),tb50(i),tb89(i),satZenithAngle(i),obsLat(i), &
                  landQualifierIndice(i)
        print *, 'GRODY: ier(i),ice(i),tpw(i),cloudLiquidWaterPathObs(i),cloudLiquidWaterPathFG(i),rain(i),snow(i)=', &
                  ier(i),ice(i),tpw(i),cloudLiquidWaterPathObs(i),cloudLiquidWaterPathFG(i),rain(i),snow(i)
      end if

    end do loopObsGrody

  end subroutine GRODY

  !------------------------------------------------------------------------------------
  ! bennartz
  !------------------------------------------------------------------------------------
  subroutine bennartz (ier, numObsToProcess, tb89, tb150, tb89FG, tb150FG, tb89FgClear, tb150FgClear, &
                       satZenithAngle, landQualifierIndice, scatIndexOverLandObs, &
                       scatIndexOverWaterObs, scatIndexOverWaterFG, &
                       cloudLiquidWaterPathObs, cloudLiquidWaterPathFG)

    !:Purpose: Compute the following parameters using 2 AMSU-B channels:
    !          - scattering index (over land and ocean).*
    !          The two channels used are: 89Ghz, 150Ghz.
    !          REGERENCES     Bennartz, R., A. Thoss, A. Dybbroe and D. B. Michelson, 
    !                         1999: Precipitation Analysis from AMSU, Nowcasting SAF, 
    !                         Swedish Meteorologicali and Hydrological Institute, 
    !                         Visiting Scientist Report, November 1999.
    implicit none
    ! arguments: 
    integer, intent(out) :: ier(numObsToProcess)        ! error return code:
    integer, intent(in)  :: numObsToProcess             ! number of points to process
    real,    intent(in)  :: tb89(:)                     ! 89Ghz AMSU-B brightness temperature (K)
    real,    intent(in)  :: tb150(:)                    ! 150Ghz AMSU-B brightness temperature (K)
    real,    intent(in)  :: tb89FG(:)                   ! 89Ghz AMSU-B brightness temperature from background (K)
    real,    intent(in)  :: tb150FG(:)                  ! 150Ghz AMSU-B brightness temperature from background (K)
    real,    intent(in)  :: tb89FgClear(:)              ! 89Ghz clear-sky brightness temperature from background (K)
    real,    intent(in)  :: tb150FgClear(:)             ! 150Ghz clear-sky brightness temperature from background (K)
    real,    intent(in)  :: satZenithAngle(:)           !  satellite zenith angle (deg.)
    integer, intent(in)  :: landQualifierIndice(:)      ! land/sea indicator (0=land;1=ocean)
    real,    intent(out) :: scatIndexOverLandObs(:)     ! scattering index over land
    real,    intent(out) :: scatIndexOverWaterObs(:)    ! scattering index over water from observation
    real,    intent(out) :: scatIndexOverWaterFG(:)     ! scattering index over water from background
    real,    intent(out) :: cloudLiquidWaterPathObs(:)  ! obs cloud liquid water content (not computed for NOW)
    real,    intent(out) :: cloudLiquidWaterPathFG(:)   ! first guess cloud liquid water content (not computed for NOW)

    ! Notes: In the case where an output parameter cannot be calculated, the
    !        value of this parameter is to to the missing value, i.e. -99.
    ! Locals: 
    integer :: i

    ! 1) Initialise output parameters
    do i = 1, numObsToProcess 
      scatIndexOverLandObs(i) = mwbg_realMissing
      scatIndexOverWaterObs(i) = mwbg_realMissing
      scatIndexOverWaterFG(i) = mwbg_realMissing
      cloudLiquidWaterPathObs(i) = mwbg_realMissing
      cloudLiquidWaterPathFG(i) = mwbg_realMissing
    end do

    ! 2) Validate input parameters
    do i = 1, numObsToProcess
      if (tb89(i)      < 120. .or. &
          tb89(i)     > 350.  .or. &
          tb150(i)    < 120.  .or. &
          tb150(i)    > 350.  .or. & 
          satZenithAngle(i)    < -90.  .or. &
          satZenithAngle(i)    >  90.  .or. & 
          landQualifierIndice(i) < 0 .or. &
          landQualifierIndice(i) > 1) then
          ier(i) = 1        
      else
          ier(i) = 0      
      end if 
    end do

    ! 3) Compute parameters
    do i = 1, numObsToProcess 
      if (ier(i) == 0) then
        if (landQualifierIndice(i) == 1) then
          if (tvs_mwAllskyAssim) then
            scatIndexOverWaterObs(i) = (tb89(i) - tb150(i)) - (tb89FgClear(i) - tb150FgClear(i))
            scatIndexOverWaterFG(i) = (tb89FG(i) - tb150FG(i)) - (tb89FgClear(i) - tb150FgClear(i))
          else
            scatIndexOverWaterObs(i) = (tb89(i) - tb150(i)) - (-39.2010 + 0.1104 * satZenithAngle(i))
          end if
        else
          scatIndexOverLandObs(i) = (tb89(i) - tb150(i)) - (0.158 + 0.0163 * satZenithAngle(i))
        end if ! if (landQualifierIndice(i) == 1)
      else if ((ier(i) /= 0) .and. (i <= 100 )) then 
        print *, 'bennartz: input Parameters are not all valid: '
        print *, 'bennartz: i, tb89(i), tb150(i), satZenithAngle(i), landQualifierIndice(i) = ', &
                            i, tb89(i), tb150(i), satZenithAngle(i), landQualifierIndice(i)
        print *, 'bennartz: ier(i), scatIndexOverLandObs(i), scatIndexOverWaterObs(i), scatIndexOverWaterFG(i)=', &
                            ier(i), scatIndexOverLandObs(i), scatIndexOverWaterObs(i), scatIndexOverWaterFG(i)
      end if ! if (ier(i) == 0)
    end do 

  end subroutine bennartz

  !--------------------------------------------------------------------------
  ! atmsMwhs2Test1Flagbit7Check
  !--------------------------------------------------------------------------
  subroutine atmsMwhs2Test1Flagbit7Check (itest, obsChannels2D, obsFlags2D, sensorIndex, qcIndicator, &
                                          actualNumChannel, numObsToProcess, stnId,  B7CHCK)

    !:Purpose:               1) test 1: Check flag bit 7 on from the first bgckAtms/bgckMwhs2 program
    !                        Includes observations flagged for cloud liquid water, scattering index,
    !                        dryness index plus failure of several QC checks.


    ! Arguments
    integer,     intent(in)                                   :: itest(:)                                         ! test number
    integer,     intent(in)                                   :: actualNumChannel                                 ! nombre de canaux des observations 
    integer,     intent(in)                                   :: numObsToProcess                                  ! nombre de tovs    
    integer,     intent(in)                                   :: obsChannels2D(actualNumChannel,numObsToProcess)  ! observations channels
    integer,     intent(inout)                                :: obsFlags2D(actualNumChannel,numObsToProcess)     ! observations channels
    integer,     intent(in)                                   :: sensorIndex                                      ! numero de satellite (i.e. indice) 
    character *9,intent(in)                                   :: stnId                                            ! identificateur du satellite
    integer,     intent(inout)                                :: B7CHCK(actualNumChannel,numObsToProcess) 
    integer,     intent(inout)                                :: qcIndicator(actualNumChannel,numObsToProcess)    ! indicateur du QC par canal
    ! Locals
    integer                                                   :: nDataIndex
    integer                                                   :: nChannelIndex
    integer                                                   :: testIndex
    integer                                                   :: IBIT 

    testIndex = 1
    if ( itest(testIndex) == 1 ) then
      do nDataIndex=1,numObsToProcess
        do nChannelIndex=1,actualNumChannel
          IBIT = AND(obsFlags2D(nChannelIndex,nDataIndex), 2**7)
          if ( IBIT /= 0  ) then
            qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
            B7CHCK(nChannelIndex,nDataIndex) = 1
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
            rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                 rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1
            if ( mwbg_debug ) then
              write(*,*)stnId(2:9),' first bgckAtms/bgckMwhs2 program REJECT.', &
                        'CHANNEL=', obsChannels2D(nChannelIndex,nDataIndex), &
                        ' obsFlags2D= ',obsFlags2D(nChannelIndex,nDataIndex)
            end if
          end if
        end do
      end do
    end if

  end subroutine atmsMwhs2Test1Flagbit7Check

  !--------------------------------------------------------------------------
  ! atmsMwhs2Test2TopographyCheck
  !--------------------------------------------------------------------------
  subroutine atmsMwhs2Test2TopographyCheck(itest, obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                           modelInterpTerrain, obsFlags2D, ICHTOPO, ZCRIT, B7CHCK, qcIndicator)

    !:Purpose:               1) test 2: Topography check (partial)

    ! Arguments
    integer,     intent(in)                                   :: itest(:)                                         ! test number
    integer,     intent(in)                                   :: actualNumChannel                                 ! nombre de canaux des observations 
    integer,     intent(in)                                   :: numObsToProcess                                  ! nombre de tovs    
    integer,     intent(in)                                   :: obsChannels2D(actualNumChannel,numObsToProcess)  ! observations channels
    integer,     intent(in)                                   :: sensorIndex                                      ! numero de satellite (i.e. indice) 
    character *9,intent(in)                                   :: stnId                                            ! identificateur du satellite
    real,        intent(in)                                   :: modelInterpTerrain(numObsToProcess)              ! topo aux point d'obs
    integer,     intent(inout)                                :: obsFlags2D(actualNumChannel,numObsToProcess)     ! marqueur de radiance 
    integer,     intent(in)                                   :: ICHTOPO(:) 
    real ,       intent(in)                                   :: ZCRIT(:)
    integer,     intent(inout)                                :: qcIndicator(actualNumChannel,numObsToProcess)    ! indicateur du QC par canal
    integer,     intent(inout)                                :: B7CHCK(actualNumChannel,numObsToProcess)
    ! Locals
    integer                                                   :: nDataIndex
    integer                                                   :: nChannelIndex
    integer                                                   :: testIndex
    integer                                                   :: INDXTOPO 

    testIndex = 2
    if ( itest(testIndex) == 1 ) then
      do nDataIndex=1,numObsToProcess
        do nChannelIndex=1,actualNumChannel
          INDXTOPO = ISRCHEQI(ICHTOPO,obsChannels2D(nChannelIndex,nDataIndex))
          if ( INDXTOPO > 0 ) then
            if (modelInterpTerrain(nDataIndex) >= ZCRIT(INDXTOPO)) then
              qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
              obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
              obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**18)
              rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                   rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1
              if ( B7CHCK(nChannelIndex,nDataIndex) == 0 ) then
                rejectionCodArray2(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                   rejectionCodArray2(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1                 
              end if
              if ( mwbg_debug ) then
                write(*,*)stnId(2:9),' TOPOGRAPHY REJECT.', &
                          'CHANNEL=', obsChannels2D(nChannelIndex,nDataIndex), &
                          ' TOPO= ',modelInterpTerrain(nDataIndex)
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
  subroutine atmsMwhs2Test3UncorrectedTbCheck(itest, obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                              RESETQC, obsFlags2D, B7CHCK, qcIndicator)

    !:Purpose:                       Test 3: Uncorrected Tb check (single)
    !                                Uncorrected datum (flag bit #6 off). 
    !                                In this case switch bit 11 ON.

    ! Arguments
    integer,     intent(in)                                   :: itest(:)                                         ! test number
    integer,     intent(in)                                   :: actualNumChannel                                 ! nombre de canaux des observations 
    integer,     intent(in)                                   :: numObsToProcess                                  ! nombre de tovs    
    integer,     intent(in)                                   :: obsChannels2D(actualNumChannel,numObsToProcess)  ! observations channels
    integer,     intent(in)                                   :: sensorIndex                                      ! numero de satellite (i.e. indice) 
    character *9,intent(in)                                   :: stnId                                            ! identificateur du satellite
    logical,     intent(in)                                   :: RESETQC                                          ! resetqc logical
    integer,     intent(inout)                                :: obsFlags2D(actualNumChannel,numObsToProcess)     ! marqueur de radiance 
    integer,     intent(inout)                                :: qcIndicator(actualNumChannel,numObsToProcess)    ! indicateur du QC par canal
    integer,     intent(inout)                                :: B7CHCK(actualNumChannel,numObsToProcess)
    ! Locals
    integer                                                   :: nDataIndex
    integer                                                   :: nChannelIndex
    integer                                                   :: testIndex
    integer                                                   :: max 
    integer                                                   :: IBIT 

 
    if (.not. RESETQC) then
      testIndex = 3
      if ( itest(testIndex) == 1 ) then
        do nDataIndex=1,numObsToProcess
          do nChannelIndex=1,actualNumChannel
            IBIT = AND(obsFlags2D(nChannelIndex,nDataIndex), 2**6)
            if ( IBIT == 0  ) then
              qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
              obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**11)
              rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                 rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1
              if ( B7CHCK(nChannelIndex,nDataIndex) == 0 ) then
                rejectionCodArray2(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                    rejectionCodArray2(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1                 
              end if
              if ( mwbg_debug ) then
                write(*,*)stnId(2:9),' UNCORRECTED TB REJECT.', &
                          'CHANNEL=', obsChannels2D(nChannelIndex,nDataIndex), &
                          ' obsFlags= ',obsFlags2D(nChannelIndex,nDataIndex)
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
  subroutine atmsTest4RogueCheck(itest, obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, ROGUEFAC, &
                                 waterobs, ompTb2D, cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                                 newInformationFlag, ISFCREJ, ICH2OMPREJ, &
                                 obsFlags2D, B7CHCK, qcIndicator)

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
    integer,     intent(in)              :: itest(:)                                        ! test number
    integer,     intent(in)              :: actualNumChannel                                ! nombre de canaux des observations 
    integer,     intent(in)              :: numObsToProcess                                 ! nombre de tovs
    integer,     intent(in)              :: obsChannels2D(actualNumChannel,numObsToProcess) ! observations channels
    integer,     intent(in)              :: sensorIndex                                     ! numero de satellite (i.e. indice) 
    character *9,intent(in)              :: stnId                                           ! identificateur du satellite
    real,        intent(in)              :: ROGUEFAC(:)                                     ! rogue factor 
    logical,     intent(in)              :: waterobs(:)                                     ! open water obs
    real,        intent(in)              :: ompTb2D(actualNumChannel,numObsToProcess)       ! radiance o-p 
    real,        intent(in)              :: cloudLiquidWaterPathObs(:)
    real,        intent(in)              :: cloudLiquidWaterPathFG(:)
    integer,     intent(in)              :: newInformationFlag(numObsToProcess)             ! data flag newInformationFlag  
    integer,     intent(in)              :: ISFCREJ(:)
    integer,     intent(in)              :: ICH2OMPREJ(:)
    integer,     intent(inout)           :: obsFlags2D(actualNumChannel,numObsToProcess)    ! marqueur de radiance 
    integer,     intent(inout)           :: qcIndicator(actualNumChannel,numObsToProcess)   ! indicateur du QC par canal
    integer,     intent(inout)           :: B7CHCK(actualNumChannel,numObsToProcess)

    ! Locals
    integer                              :: channelval
    integer                              :: nDataIndex
    integer                              :: nChannelIndex
    integer                              :: testIndex
    integer                              :: INDXCAN 
    real                                 :: XCHECKVAL
    real                                 :: clwThresh1 
    real                                 :: clwThresh2
    real                                 :: errThresh1 
    real                                 :: errThresh2
    real                                 :: sigmaObsErrUsed  
    real                                 :: clwObsFGaveraged 
    logical                              :: SFCREJCT
    logical                              :: CH2OMPREJCT
    logical                              :: IBIT 

    testIndex = 4
    if ( itest(testIndex) == 1 ) then
      do nDataIndex=1,numObsToProcess
        SFCREJCT = .FALSE.
        CH2OMPREJCT = .FALSE.
        do nChannelIndex=1,actualNumChannel
          channelval = obsChannels2D(nChannelIndex,nDataIndex)
          ! using state-dependent obs error only over water.
          ! obs over sea-ice will be rejected in test 15.
          if ( tvs_mwAllskyAssim .and. oer_useStateDepSigmaObs(channelval,sensorIndex) &
                .and. waterobs(nDataIndex) ) then
            clwThresh1 = oer_cldPredThresh(channelval,sensorIndex,1)
            clwThresh2 = oer_cldPredThresh(channelval,sensorIndex,2)
            errThresh1 = oer_errThreshAllsky(channelval,sensorIndex,1)
            errThresh2 = oer_errThreshAllsky(channelval,sensorIndex,2)
            clwObsFGaveraged = 0.5 * (cloudLiquidWaterPathObs(nDataIndex) + cloudLiquidWaterPathFG(nDataIndex))
            if (cloudLiquidWaterPathObs(nDataIndex) == mwbg_realMissing .or. &
                cloudLiquidWaterPathFG(nDataIndex) == mwbg_realMissing) then
              sigmaObsErrUsed = MPC_missingValue_R4
            else
              sigmaObsErrUsed = calcStateDepObsErr_r4(clwThresh1,clwThresh2,errThresh1, &
                                                        errThresh2,clwObsFGaveraged)
            end if
          else
            sigmaObsErrUsed = oer_toverrst(channelval,sensorIndex)
          end if
          ! For sigmaObsErrUsed=MPC_missingValue_R4 (cloudLiquidWaterPathObs[FG]=mwbg_realMissing
          ! in all-sky mode), the observation is flagged for rejection in 
          ! mwbg_reviewAllCritforFinalFlagsAtms.
          XCHECKVAL = ROGUEFAC(channelval) * sigmaObsErrUsed
          if ( ompTb2D(nChannelIndex,nDataIndex)      /= mwbg_realMissing    .and. &
              ABS(ompTb2D(nChannelIndex,nDataIndex))  >= XCHECKVAL .and. &
              sigmaObsErrUsed /= MPC_missingValue_R4 ) then
            qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**16)
            rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) =  &
               rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) + 1
            if ( B7CHCK(nChannelIndex,nDataIndex) == 0 ) then
              rejectionCodArray2(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                 rejectionCodArray2(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1                 
            end if
            if ( mwbg_debug ) then
              write(*,*)stnId(2:9),'ROGUE CHECK REJECT.NO.', &
                     ' OBS = ',nDataIndex, &
                     ' CHANNEL= ',obsChannels2D(nChannelIndex,nDataIndex), &
                     ' CHECK VALUE= ',XCHECKVAL, &
                     ' TBOMP= ',ompTb2D(nChannelIndex,nDataIndex), &
                     ' TOVERRST= ',oer_toverrst(channelval,sensorIndex)
            end if
            if ( channelval == 1 .or. &
                 channelval == 2 .or. &
                 channelval == 3    ) then
              SFCREJCT = .TRUE.
            end if
          end if
          if ( channelval == 17 .and. ompTb2D(nChannelIndex,nDataIndex) /= mwbg_realMissing .and. &
              ABS(ompTb2D(nChannelIndex,nDataIndex)) > 5.0 ) then
            CH2OMPREJCT = .TRUE.
          end if
        end do

        if ( SFCREJCT ) then
          do nChannelIndex=1,actualNumChannel
            INDXCAN = ISRCHEQI(ISFCREJ,obsChannels2D(nChannelIndex,nDataIndex))
            if ( INDXCAN /= 0 )  then
              if ( qcIndicator(nChannelIndex,nDataIndex) /= testIndex ) then
                qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
                obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
                obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**16)
                rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                        rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1
                if ( B7CHCK(nChannelIndex,nDataIndex) == 0 ) then
                  rejectionCodArray2(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                     rejectionCodArray2(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1                 
                end if
              end if
            end if
          end do
        end if

        !  amsub channels 17-22 obs are rejected if, for ch17 ABS(O-P) > 5K
        !    Apply over open water only (bit 0 ON in QC integer newInformationFlag).
        !    Only apply if obs not rejected in this test already.
        IBIT = AND(newInformationFlag(nDataIndex), 2**0)
        if ( CH2OMPREJCT .and. (IBIT /= 0) ) then
          do nChannelIndex=1,actualNumChannel
            INDXCAN = ISRCHEQI (ICH2OMPREJ,obsChannels2D(nChannelIndex,nDataIndex))
            if ( INDXCAN /= 0 )  then
              if ( qcIndicator(nChannelIndex,nDataIndex) /= testIndex ) then
                qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
                obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
                obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**16)
                rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                        rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1
                if ( B7CHCK(nChannelIndex,nDataIndex) == 0 ) then
                  rejectionCodArray2(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                     rejectionCodArray2(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1                 
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
  subroutine Mwhs2Test4RogueCheck(itest, obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, ROGUEFAC, &
                                  waterobs, ompTb2D, cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                                  newInformationFlag, ICH2OMPREJ, obsFlags2D, B7CHCK, qcIndicator)

    !:Purpose:                         test 4: "Rogue check" for (O-P) Tb residuals out of range.
    !                                  (single/full).
    !                                  Also, over WATER remove CH.10-15 if CH.10 |O-P|>5K (full)
    !                                  Les observations, dont le residu (O-P)
    !                                  depasse par un facteur (roguefac) l'erreur totale des TOVS.
    !                                  OVER OPEN WATER
    !                                  ch. 10 Abs(O-P) > 5K produces rejection of all ATMS amsub channels 10-15.

    ! Arguments
    integer,     intent(in)              :: itest(:)                                        ! test number
    integer,     intent(in)              :: actualNumChannel                                ! nombre de canaux des observations
    integer,     intent(in)              :: numObsToProcess                                 ! nombre de tovs    
    integer,     intent(in)              :: obsChannels2D(actualNumChannel,numObsToProcess) ! observations channels
    integer,     intent(in)              :: sensorIndex                                     ! numero de satellite (i.e. indice)
    character *9,intent(in)              :: stnId                                           ! identificateur du satellite
    real,        intent(in)              :: ROGUEFAC(:)                                     ! rogue factor
    logical,     intent(in)              :: waterobs(:)                                     ! open water obs
    real,        intent(in)              :: ompTb2D(actualNumChannel,numObsToProcess)       ! radiance o-p
    real,        intent(in)              :: cloudLiquidWaterPathObs(:)
    real,        intent(in)              :: cloudLiquidWaterPathFG(:)
    integer,     intent(in)              :: newInformationFlag(numObsToProcess)             ! data flag newInformationFlag
    integer,     intent(in)              :: ICH2OMPREJ(:)
    integer,     intent(inout)           :: obsFlags2D(actualNumChannel,numObsToProcess)    ! marqueur de radiance
    integer,     intent(inout)           :: qcIndicator(actualNumChannel,numObsToProcess)   ! indicateur du QC par canal
    integer,     intent(inout)           :: B7CHCK(actualNumChannel,numObsToProcess)

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
      do nDataIndex=1,numObsToProcess
        CH2OMPREJCT = .FALSE.
        do nChannelIndex=1,actualNumChannel
          channelval = obsChannels2D(nChannelIndex,nDataIndex)
          ! using state-dependent obs error only over water.
          ! obs over sea-ice will be rejected in test 15.
          if ( tvs_mwAllskyAssim .and. oer_useStateDepSigmaObs(channelval,sensorIndex) &
                .and. waterobs(nDataIndex) ) then
            clwThresh1 = oer_cldPredThresh(channelval,sensorIndex,1)
            clwThresh2 = oer_cldPredThresh(channelval,sensorIndex,2)
            sigmaThresh1 = oer_errThreshAllsky(channelval,sensorIndex,1)
            sigmaThresh2 = oer_errThreshAllsky(channelval,sensorIndex,2)
            clwObsFGaveraged = 0.5 * (cloudLiquidWaterPathObs(nDataIndex) + cloudLiquidWaterPathFG(nDataIndex))
            if ( cloudLiquidWaterPathObs(nDataIndex) == mwbg_realMissing ) then
              sigmaObsErrUsed = MPC_missingValue_R4
            else
              sigmaObsErrUsed = calcStateDepObsErr_r4(clwThresh1,clwThresh2,sigmaThresh1, &
                                                        sigmaThresh2,clwObsFGaveraged)
            end if
          else
            sigmaObsErrUsed = oer_toverrst(channelval,sensorIndex)
          end if
          ! For sigmaObsErrUsed=MPC_missingValue_R4 (cloudLiquidWaterPathObs=mwbg_realMissing
          ! in all-sky mode), the observation is flagged for rejection in
          ! mwbg_reviewAllCritforFinalFlagsMwhs2.
          XCHECKVAL = ROGUEFAC(channelval) * sigmaObsErrUsed
          if ( ompTb2D(nChannelIndex,nDataIndex)      /= mwbg_realMissing    .and. &
              ABS(ompTb2D(nChannelIndex,nDataIndex))  >= XCHECKVAL .and. &
              sigmaObsErrUsed /= MPC_missingValue_R4 ) then
            qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
            obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**16)
            rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) =  &
               rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) + 1
            if ( B7CHCK(nChannelIndex,nDataIndex) == 0 ) then
              rejectionCodArray2(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                 rejectionCodArray2(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1
            end if
            if ( mwbg_debug ) then
              write(*,*)stnId(2:9),'ROGUE CHECK REJECT.NO.', &
                     ' OBS = ',nDataIndex, &
                     ' CHANNEL= ',obsChannels2D(nChannelIndex,nDataIndex), &
                     ' CHECK VALUE= ',XCHECKVAL, &
                     ' TBOMP= ',ompTb2D(nChannelIndex,nDataIndex), &
                     ' TOVERRST= ',oer_toverrst(channelval,sensorIndex)
            end if
          end if
          if ( channelval == 10 .and. ompTb2D(nChannelIndex,nDataIndex) /= mwbg_realMissing .and. &
              ABS(ompTb2D(nChannelIndex,nDataIndex)) > 5.0 ) then
            CH2OMPREJCT = .TRUE.
          end if
        end do

        !    Channels 10-15 are rejected if, for ch10 ABS(O-P) > 5K
        !    Apply over open water only (bit 0 ON in QC integer newInformationFlag).
        !    Only apply if obs not rejected in this test already.
        IBIT = AND(newInformationFlag(nDataIndex), 2**0)
        if ( CH2OMPREJCT .and. (IBIT /= 0) ) then
          do nChannelIndex=1,actualNumChannel
            INDXCAN = ISRCHEQI(ICH2OMPREJ,obsChannels2D(nChannelIndex,nDataIndex))
            if ( INDXCAN /= 0 )  then
              if ( qcIndicator(nChannelIndex,nDataIndex) /= testIndex ) then
                qcIndicator(nChannelIndex,nDataIndex) = MAX(qcIndicator(nChannelIndex,nDataIndex),testIndex)
                obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**9)
                obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**16)
                rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                        rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1
                if ( B7CHCK(nChannelIndex,nDataIndex) == 0 ) then
                  rejectionCodArray2(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                     rejectionCodArray2(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1
                end if
              end if
            end if
          end do
        end if

      end do
    end if

  end subroutine Mwhs2Test4RogueCheck

  !--------------------------------------------------------------------------
  ! atmsMwhs2Test5ChannelSelectionUsingTovutil
  !--------------------------------------------------------------------------
  subroutine atmsMwhs2Test5ChannelSelectionUsingTovutil(itest, obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, &
                                                        stnId, obsFlags2D, qcIndicator)

    !:Purpose: test 5: Channel selection using array oer_tovutil(chan,sat)
    !          oer_tovutil = 0 (blacklisted), 1 (assmilate)

    ! Arguments
    integer,     intent(in)               :: itest(:)                                         ! test number
    integer,     intent(in)               :: actualNumChannel                                 ! nombre de canaux des observations 
    integer,     intent(in)               :: numObsToProcess                                  ! nombre de tovs    
    integer,     intent(in)               :: obsChannels2D(actualNumChannel,numObsToProcess)  ! observations channels
    integer,     intent(in)               :: sensorIndex                                      ! numero de satellite (i.e. indice) 
    character *9,intent(in)               :: stnId                                            ! identificateur du satellite
    integer,     intent(in)               :: qcIndicator(actualNumChannel,numObsToProcess)    ! indicateur du QC par canal
    integer,     intent(inout)            :: obsFlags2D(actualNumChannel,numObsToProcess)     ! marqueur de radiance 

    ! Locals
    integer                               :: channelval
    integer                               :: nDataIndex
    integer                               :: nChannelIndex
    integer                               :: testIndex

    testIndex = 5
    if ( itest(testIndex) == 1 ) then
      do nDataIndex=1,numObsToProcess
        do nChannelIndex=1,actualNumChannel
           channelval = obsChannels2D(nChannelIndex,nDataIndex)
           if ( oer_tovutil(channelval,sensorIndex) == 0 ) then
             obsFlags2D(nChannelIndex,nDataIndex) = OR(obsFlags2D(nChannelIndex,nDataIndex),2**8)
             rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex) = &
                   rejectionCodArray(testIndex,obsChannels2D(nChannelIndex,nDataIndex),sensorIndex)+ 1
             if ( mwbg_debug ) then
               write(*,*)stnId(2:9),'CHANNEL REJECT: ', &
                      ' OBS = ',nDataIndex, &
                      ' CHANNEL= ',channelval                  
             end if
           end if
        end do
      end do      
    end if

    if ( mwbg_debug ) then
      write(*,*) 'qcIndicator = ',((qcIndicator(nChannelIndex,nDataIndex),nChannelIndex=1,actualNumChannel),nDataIndex=1,numObsToProcess)
    end if

  end subroutine atmsMwhs2Test5ChannelSelectionUsingTovutil

  !--------------------------------------------------------------------------
  ! mwbg_tovCheckAtms 
  !--------------------------------------------------------------------------
  subroutine mwbg_tovCheckAtms(obsLat, obsLon, landQualifierIndice, terrainTypeIndice, satZenithAngle, &
                               obsQcFlag2, obsQcFlag1, obsChannels, obsTb, obsTbBiasCorr, ompTb, qcIndicator, &
                               actualNumChannel, numObsToProcess, sensorIndex, newInformationFlag, &
                               satScanPosition, modelInterpTerrain, obsGlobalMarker, obsFlags, &
                               cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                               scatIndexOverWaterObs, scatIndexOverWaterFG, &
                               stnId, RESETQC)

    !:Purpose:                   Effectuer le controle de qualite des radiances tovs.
    !

    implicit none
    !Arguments
    integer, intent(in)              :: actualNumChannel                    ! nombre de canaux des observations
    integer, intent(in)              :: numObsToProcess                     ! nombre de tovs
    real,    intent(in)              :: obsLat(:)
    real,    intent(in)              :: obsLon(:)
    integer, intent(in)              :: landQualifierIndice(:)
    integer, intent(in)              :: terrainTypeIndice(:)
    real,    intent(inout)           :: satZenithAngle(:)
    integer, intent(in)              :: obsQcFlag2(:)
    integer, intent(in)              :: obsQcFlag1(:,:)
    integer, intent(inout)           :: obsGlobalMarker(:)
    integer, intent(in)              :: satScanPosition(numObsToProcess)    ! position sur le "scan"
    integer, intent(in)              :: obsChannels(:)                      ! canaux des observations
    integer, intent(in)              :: sensorIndex                         ! numero de satellite (i.e. indice)
    real, intent(inout)              :: obsTb(:)                            ! radiances
    real, intent(in)                 :: obsTbBiasCorr(:)                    ! correction aux radiances
    real, intent(in)                 :: ompTb(:)                            ! residus (o-p)
    real, intent(in)                 :: modelInterpTerrain(numObsToProcess) ! topographie du modele
    integer, allocatable, intent(out):: newInformationFlag(:)               ! flag to identify all obs pts in report
    !                                                                         as being over land/ice, cloudy, bad IWV
    character *9, intent(in)         :: stnId                               ! identificateur du satellite
    logical, intent(in)              :: RESETQC                             ! reset du controle de qualite?
    integer,allocatable, intent(out) :: qcIndicator(:,:)                    ! indicateur controle de qualite tovs par canal 
    !                                                                         =0 ok, >0 rejet,
    integer, intent(inout)           :: obsFlags(:)                         ! marqueurs des radiances
                                                                            ! satellite, critere et par canal
                                                                            ! (chech n2) par satellite, critere et par canal
    real, allocatable, intent(out)   :: cloudLiquidWaterPathObs(:)
    real, allocatable, intent(out)   :: cloudLiquidWaterPathFG(:)
    real, allocatable, intent(out)   :: scatIndexOverWaterObs(:)            ! scattering index over water from observation
    real, allocatable, intent(out)   :: scatIndexOverWaterFG(:)             ! scattering index over water from background

    !locals
    real                             :: ompTb2D(actualNumChannel,numObsToProcess)
    integer                          :: obsChannels2D(actualNumChannel,numObsToProcess)
    integer                          :: obsFlags2D(actualNumChannel,numObsToProcess)
    integer, allocatable             :: calcLandQualifierIndice(:)
    integer, allocatable             :: calcTerrainTypeIndice(:)
    integer,allocatable              :: KCHKPRF(:)

    logical, allocatable             :: waterobs(:)
    logical, allocatable             :: grossrej(:)
    logical                          :: reportHasMissingTb
    logical, allocatable             :: qcRejectLogic(:,:)
    logical, allocatable             :: cloudobs(:)
    logical, allocatable             :: iwvreject(:)
    logical, allocatable             :: precipobs(:)
    real, allocatable                :: zdi(:)
    real, allocatable                :: scatec(:)
    real, allocatable                :: scatbg(:)
    real, allocatable                :: SeaIce(:)
    real, allocatable                :: riwv(:)

    integer, parameter               :: maxScanAngleAMSU = 96
    integer, parameter               :: ilsmOpt = 1   ! OPTION for values of MG (land/sea mask) and LG (ice) 
                                                      ! at each observation point using values on 5x5 mesh 
                                                      ! centered at each point.
                                                      ! ilsmOpt = 1 --> use MAX value from all 25 mesh points
                                                      ! ilsmOpt = 2 --> use value at central mesh point (obs location)
                                                      ! ilsmOpt = 3 --> use AVG value from all 25 mesh points

    integer                          :: iRej
    integer                          :: iNumSeaIce
    integer                          :: JI
    integer                          :: JJ
    integer                          :: kk
    integer                          :: INDX
    integer                          :: ISFCREJ(8)
    integer                          :: ICH2OMPREJ(6)
    integer                          :: B7CHCK(actualNumChannel,numObsToProcess)
    real, allocatable                :: ROGUEFAC(:)
    real                             :: ZCRIT(5)
    integer                          :: ITEST(mwbg_maxNumTest)
    integer                          :: chanIgnoreInAllskyGenCoeff(6)
    integer                          :: ICHTOPO(5)
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

    call utl_reAllocate(ROGUEFAC, actualNumChannel+tvs_channelOffset(sensorIndex))
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
    chanIgnoreInAllskyGenCoeff(:) = (/ 1, 2, 3, 4, 5, 6/)

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
    call atmsMwhs2landIceMask(numObsToProcess, obsLat, obsLon, calcLandQualifierIndice, &
                              calcTerrainTypeIndice, waterobs, ilsmOpt)

    !###############################################################################
    ! STEP 2 ) Check for values of TB that are missing or outside physical limits.
    !###############################################################################
    call mwbg_grossValueCheck(numObsToProcess, actualNumChannel, obsTb, obsTbBiasCorr, 50., 380., grossrej)

    !###############################################################################
    ! STEP 3 ) Preliminary QC checks --> set qcRejectLogic(numObsToProcess,actualNumChannel)=.true.
    !          for data that fail QC
    !###############################################################################
    call mwbg_firstQcCheckAtms(satZenithAngle, landQualifierIndice, terrainTypeIndice, obsLat, obsLon, obsTb, satScanPosition, &
                               actualNumChannel, numObsToProcess, qcRejectLogic, grossrej, calcLandQualifierIndice, calcTerrainTypeIndice, &
                               obsQcFlag1, obsQcFlag2, obsChannels, reportHasMissingTb)

    if ( reportHasMissingTb ) numReportWithMissingTb = numReportWithMissingTb + 1
    !  Exclude problem points from further calculations
    do kk = 1,numObsToProcess
      if ( COUNT(qcRejectLogic(kk,:)) == actualNumChannel ) grossrej(kk) = .true.
    end do

    !###############################################################################
    ! STEP 4 ) mwbg_nrlFilterAtms returns cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, scatec, scatbg and also does sea-ice
    !          detection missing value for  cloudLiquidWaterPathObs, scatec, scatbg  is -99.0 (e.g. over
    !          land or sea-ice).Sets calcTerrainTypeIndice=0 (sea ice) for points where retrieved SeaIce
    !          >=0.55. Does nothing if calcTerrainTypeIndice=0 (sea ice) and retrieved SeaIce<0.55.
    !###############################################################################
    call mwbg_nrlFilterAtms(numObsToProcess, actualNumChannel, obsTb, ompTb, obsTbBiasCorr, satZenithAngle, obsLat, &
                            calcLandQualifierIndice, calcTerrainTypeIndice, waterobs, grossrej, &
                            cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                            scatec, scatbg, iNumSeaIce, iRej, SeaIce)

    seaIcePointNum = seaIcePointNum + iNumSeaIce
    clwMissingPointNum = clwMissingPointNum + iRej

    !###############################################################################
    ! STEP 5 ) Apply NRL cloud filter, scattering index and sea-ice detection algorithms
    !          to OPEN WATER (waterobs=true) points.
    ! Points with SeaIce>0.55 are set to sea-ice points (waterobs --> false)
    !###############################################################################
    call mwbg_flagDataUsingNrlCritAtms(numObsToProcess, actualNumChannel, obsTb, obsTbBiasCorr, cloudLiquidWaterPathObs, &
                                       scatec, scatbg, SeaIce, grossrej, waterobs, mwbg_useUnbiasedObsForClw, &
                                       iwvreject, cloudobs, precipobs, cldcnt , newInformationFlag, riwv, zdi)

    !###############################################################################
    ! STEP 6 ) ! Review all the checks previously made to determine which obs are to be
    !            accepted for assimilation and which are to be flagged for exclusion
    !            (obsFlags).
    !            grossrej()  = .true. if any channel had a gross error at the point
    !            cloudobs()  = .true. if CLW > clw_atms_nrl_LTrej (0.175) or precipobs
    !            precipobs() = .true. if precip. detected through NRL scattering indices
    !            waterobs()  = .true. if open water point
    !            iwvreject() = .true. if Mean 183 Ghz [ch. 18-22] Tb < 240K (too dry
    !            for ch.20-22 over land)
    !###############################################################################
    call mwbg_reviewAllCritforFinalFlagsAtms(numObsToProcess, actualNumChannel, qcRejectLogic, grossrej, waterobs, &
                                             precipobs, cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, scatec, scatbg, &
                                             scatIndexOverWaterObs, scatIndexOverWaterFG, iwvreject, riwv, obsFlags, &
                                             obsGlobalMarker, zdi, newInformationFlag, drycnt, landcnt, &
                                             rejcnt, iwvcnt, pcpcnt, flgcnt, &
                                             chanIgnoreInAllskyGenCoeff, obsChannels)

    !###############################################################################
    ! PART 2 TESTS:
    !###############################################################################

    ! copy the original input 1D array to 2D array. The 2D arrays are used in this s/r.
    call copy1Dimto2DimRealArray(ompTb, actualNumChannel, numObsToProcess, ompTb2D)
    call copy1Dimto2DimIntegerArray(obsChannels, actualNumChannel, numObsToProcess, obsChannels2D)
    call copy1Dimto2DimIntegerArray(obsFlags, actualNumChannel, numObsToProcess, obsFlags2D)
    ! allocations
    call utl_reAllocate(kchkprf, numObsToProcess)
    call utl_reAllocate(qcIndicator, actualNumChannel, numObsToProcess)
    !  Initialisations
    qcIndicator(:,:) = 0
    B7CHCK(:,:) = 0

    if ( RESETQC ) obsFlags2D(:,:) = 0

    ! 1) test 1: Check flag bit 7 on from the first bgckAtms program
    !  Includes observations flagged for cloud liquid water, scattering index,
    !  dryness index plus failure of several QC checks.
    call atmsMwhs2Test1Flagbit7Check (itest, obsChannels2D, obsFlags2D, sensorIndex, qcIndicator, &
                                      actualNumChannel, numObsToProcess, stnId,  B7CHCK)

    ! 2) test 2: Topography check (partial)
    call atmsMwhs2Test2TopographyCheck (itest, obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                        modelInterpTerrain, obsFlags2D, ICHTOPO, ZCRIT, B7CHCK, qcIndicator)

    ! 3) test 3: Uncorrected Tb check (single)
    !  Uncorrected datum (flag bit #6 off). In this case switch bit 11 ON.
    call atmsMwhs2Test3UncorrectedTbCheck (itest, obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                           RESETQC, obsFlags2D, B7CHCK, qcIndicator)

    ! 4) test 4: "Rogue check" for (O-P) Tb residuals out of range. (single/full)
    !             Also, over WATER remove CH.17-22 if CH.17 |O-P|>5K (partial)
    !  Les observations, dont le residu (O-P) depasse par un facteur (roguefac)
    !   l'erreur totale des TOVS.
    !  N.B.: a reject by any of the 3 amsua surface channels 1-3 produces the
    !           rejection of ATMS sfc/tropospheric channels 1-6 and 16-17.
    !  OVER OPEN WATER
    !    ch. 17 Abs(O-P) > 5K produces rejection of all ATMS amsub channels 17-22.
    call atmsTest4RogueCheck (itest, obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, ROGUEFAC, &
                              waterobs, ompTb2D, cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                              newInformationFlag, ISFCREJ, ICH2OMPREJ, &
                              obsFlags2D, B7CHCK, qcIndicator)

    ! 5) test 5: Channel selection using array oer_tovutil(chan,sat)
    !  oer_tovutil = 0 (blacklisted)
    !                1 (assmilate)
    call atmsMwhs2Test5ChannelSelectionUsingTovutil(itest, obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, &
                                                    stnId, obsFlags2D, qcIndicator)

    !  Synthese de la controle de qualite au niveau de chaque point
    !  d'observation. Code:
    !            =0 aucun rejet, >0 au moins un canal rejete.
    do JJ = 1, numObsToProcess
      KCHKPRF(JJ) = 0
      do JI = 1, actualNumChannel
        KCHKPRF(JJ) = MAX(KCHKPRF(JJ),qcIndicator(JI,JJ))
      end do
    end do

    if ( mwbg_debug ) then
      write(*,*)'KCHKPRF = ',(KCHKPRF(JJ),JJ=1,numObsToProcess)
    end if

    ! Copy the modified FLAG to the 1D array, used outside this s/r.
    do JJ = 1, numObsToProcess
      do JI = 1, actualNumChannel
        INDX = (JJ-1) * actualNumChannel + JI
        obsFlags(INDX) = obsFlags2D(JI,JJ)
      end do
    end do

    ! reset global marker flag (55200) and mark it if observtions are rejected
    call resetQcCases(RESETQC, KCHKPRF, obsGlobalMarker)

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

      write(*,*) '   Meaning of newInformationFlag flag bits: '
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
      write(*,*) '   New Element 25174 in BURP file = newInformationFlag flag'
      write(*,*) ' '
    end if

  end subroutine mwbg_tovCheckAtms

  !--------------------------------------------------------------------------
  ! mwbg_tovCheckMwhs2
  !--------------------------------------------------------------------------
  subroutine mwbg_tovCheckMwhs2(obsLat, obsLon, landQualifierIndice, terrainTypeIndice, satZenithAngle, &
                                obsChannels, obsTb, obsTbBiasCorr, ompTb, qcIndicator, &
                                actualNumChannel, numObsToProcess, sensorIndex, &
                                newInformationFlag, satScanPosition, &
                                modelInterpTerrain, obsGlobalMarker, obsFlags, &
                                cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                                scatIndexOverWaterObs, scatIndexOverWaterFG, &
                                stnId, RESETQC, modLSQ, lastHeader)

    !:Purpose:                   Effectuer le controle de qualite des radiances tovs.

    implicit none
    !Arguments
    integer, intent(in)              :: actualNumChannel                    ! nombre de canaux des observations
    integer, intent(in)              :: numObsToProcess                     ! nombre de tovs
    real,    intent(in)              :: obsLat(:)
    real,    intent(in)              :: obsLon(:)
    integer, intent(in)              :: landQualifierIndice(:)
    integer, intent(in)              :: terrainTypeIndice(:)
    real,    intent(inout)           :: satZenithAngle(:)
    integer, intent(inout)           :: obsGlobalMarker(:)
    integer, intent(in)              :: satScanPosition(numObsToProcess)    ! position sur le "scan"
    integer, intent(in)              :: obsChannels(:)                      ! canaux des observations
    integer, intent(in)              :: sensorIndex                         ! numero de satellite (i.e. indice)
    real, intent(inout)              :: obsTb(:)                            ! radiances
    real, intent(in)                 :: obsTbBiasCorr(:)                    ! correction aux radiances
    real, intent(in)                 :: ompTb(:)                            ! residus (o-p)
    real, intent(in)                 :: modelInterpTerrain(numObsToProcess) ! topographie du modele
    integer, allocatable, intent(out):: newInformationFlag(:)               ! flag to identify all obs pts in report
    !                                                                         as being over land/ice, cloudy, bad IWV
    character *9, intent(in)         :: stnId                               ! identificateur du satellite
    logical, intent(in)              :: RESETQC                             ! reset du controle de qualite?
    logical, intent(in)              :: modLSQ                              ! If active, recalculate values for land/sea
                                                                            ! qualifier and terrain type based on LG/MG
    logical, intent(in)              :: lastHeader                          ! active if last header
    integer,allocatable, intent(out) :: qcIndicator(:,:)                    ! indicateur controle de qualite tovs par canal
    !                                                                         =0 ok, >0 rejet,
    integer, intent(inout)           :: obsFlags(:)                         ! marqueurs des radiances
                                                                            ! satellite, critere et par canal
                                                                            ! (chech n2) par satellite, critere et par canal
    real, allocatable, intent(out)   :: cloudLiquidWaterPathObs(:)
    real, allocatable, intent(out)   :: cloudLiquidWaterPathFG(:)
    real, allocatable, intent(out)   :: scatIndexOverWaterObs(:)            ! scattering index over water from observation
    real, allocatable, intent(out)   :: scatIndexOverWaterFG(:)             ! scattering index over water from background

    !locals
    real                             :: ompTb2D(actualNumChannel,numObsToProcess)
    integer                          :: obsChannels2D(actualNumChannel,numObsToProcess)
    integer                          :: obsFlags2D(actualNumChannel,numObsToProcess)
    integer, allocatable             :: calcLandQualifierIndice(:)
    integer, allocatable             :: calcTerrainTypeIndice(:)
    integer,allocatable              :: KCHKPRF(:)

    logical, allocatable             :: waterobs(:)
    logical, allocatable             :: grossrej(:)
    logical                          :: reportHasMissingTb
    logical, allocatable             :: qcRejectLogic(:,:)
    logical, allocatable             :: cloudobs(:)
    logical, allocatable             :: iwvreject(:)
    logical, allocatable             :: precipobs(:)
    real, allocatable                :: zdi(:)
    real, allocatable                :: scatec(:)
    real, allocatable                :: scatbg(:)
    real, allocatable                :: SeaIce(:)
    real, allocatable                :: riwv(:)

    integer, parameter               :: maxScanAngleAMSU = 98
    integer, parameter               :: ilsmOpt = 2   ! OPTION for values of MG (land/sea mask) and LG (ice) 
                                                      ! at each observation point using values on 5x5 mesh 
                                                      ! centered at each point.
                                                      ! ilsmOpt = 1 --> use MAX value from all 25 mesh points
                                                      ! ilsmOpt = 2 --> use value at central mesh point (obs location)
                                                      ! ilsmOpt = 3 --> use AVG value from all 25 mesh points
    
    integer                          :: iRej
    integer                          :: iNumSeaIce
    integer                          :: JI
    integer                          :: JJ
    integer                          :: kk
    integer                          :: INDX
    integer                          :: ICH2OMPREJ(6)
    integer                          :: B7CHCK(actualNumChannel,numObsToProcess)
    real, allocatable                :: ROGUEFAC(:)
    real                             :: ZCRIT(3)
    integer                          :: ITEST(mwbg_maxNumTest)
    integer                          :: chanIgnoreInAllskyGenCoeff(6)
    integer                          :: ICHTOPO(3)
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

    call utl_reAllocate(ROGUEFAC, actualNumChannel+tvs_channelOffset(sensorIndex))
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
    chanIgnoreInAllskyGenCoeff(:) = (/ 10, 11, 12, 13, 14, 15/)

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
    call atmsMwhs2landIceMask(numObsToProcess, obsLat, obsLon, calcLandQualifierIndice, &
                              calcTerrainTypeIndice, waterobs, ilsmOpt)

    !###############################################################################
    ! STEP 2 ) Check for values of TB that are missing or outside physical limits.
    !###############################################################################
    call mwbg_grossValueCheck(numObsToProcess, actualNumChannel, obsTb, obsTbBiasCorr, 50., 380., grossrej)

    !###############################################################################
    ! STEP 3 ) Preliminary QC checks --> set qcRejectLogic(numObsToProcess,actualNumChannel)=.true.
    !          for data that fail QC
    !###############################################################################
    call mwbg_firstQcCheckMwhs2(satZenithAngle, landQualifierIndice, terrainTypeIndice, obsLat, obsLon, obsTb, satScanPosition, &
                                actualNumChannel, numObsToProcess, qcRejectLogic, calcLandQualifierIndice, calcTerrainTypeIndice, &
                                obsChannels, reportHasMissingTb, modLSQ)

    if ( reportHasMissingTb ) numReportWithMissingTb = numReportWithMissingTb + 1
    !  Exclude problem points from further calculations
    do kk = 1,numObsToProcess
      if ( COUNT(qcRejectLogic(kk,:)) == actualNumChannel ) grossrej(kk) = .true.
    end do

    !###############################################################################
    ! STEP 4 ) mwbg_nrlFilterMwhs2 returns cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, scatec, scatbg and also does sea-ice
    !          detection missing value for  cloudLiquidWaterPathObs, scatec, scatbg  is -99.0 (e.g. over
    !          land or sea-ice).Sets calcTerrainTypeIndice=0 (sea ice) for points where retrieved SeaIce
    !          >=0.55. Does nothing if calcTerrainTypeIndice=0 (sea ice) and retrieved SeaIce<0.55.
    !###############################################################################
    call mwbg_nrlFilterMwhs2(numObsToProcess, actualNumChannel, obsTb, obsTbBiasCorr, satZenithAngle, obsLat, &
                             calcLandQualifierIndice, calcTerrainTypeIndice, waterobs, grossrej, &
                             cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                             scatec, scatbg, iNumSeaIce, iRej, SeaIce)

    seaIcePointNum = seaIcePointNum + iNumSeaIce
    clwMissingPointNum = clwMissingPointNum + iRej

    !###############################################################################
    ! STEP 5 ) Apply NRL cloud filter, scattering index and sea-ice detection algorithms
    !          to OPEN WATER (waterobs=true) points.
    ! Points with SeaIce>0.55 are set to sea-ice points (waterobs --> false)
    !###############################################################################
    call mwbg_flagDataUsingNrlCritMwhs2(numObsToProcess, actualNumChannel, obsTb, obsTbBiasCorr, cloudLiquidWaterPathObs, &
                                        scatec, SeaIce, grossrej, waterobs, mwbg_useUnbiasedObsForClw, &
                                        iwvreject, cloudobs, precipobs, cldcnt , newInformationFlag, riwv, zdi)

    !###############################################################################
    ! STEP 6 ) ! Review all the checks previously made to determine which obs are to be
    !            accepted for assimilation and which are to be flagged for exclusion
    !            (obsFlags).
    !            grossrej()  = .true. if any channel had a gross error at the point
    !            cloudobs()  = .true. if CLW > clw_mwhs2_nrl_LTrej (0.175) or precipobs
    !            precipobs() = .true. if precip. detected through NRL scattering indices
    !            waterobs()  = .true. if open water point
    !            iwvreject() = .true. if Mean 183 Ghz [ch. 11-15] Tb < 240K (too dry
    !            for ch.11-13 over land)
    !###############################################################################
    call mwbg_reviewAllCritforFinalFlagsMwhs2(numObsToProcess, actualNumChannel, qcRejectLogic, grossrej, calcTerrainTypeIndice, waterobs, &
                                              precipobs, cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, scatec, scatbg, &
                                              scatIndexOverWaterObs, scatIndexOverWaterFG, iwvreject, riwv, obsFlags, &
                                              obsGlobalMarker, zdi, newInformationFlag, allcnt, drycnt, landcnt, &
                                              rejcnt, iwvcnt, pcpcnt, flgcnt, &
                                              chanIgnoreInAllskyGenCoeff, obsChannels)

    !###############################################################################
    ! PART 2 TESTS:
    !###############################################################################

    ! copy the original input 1D array to 2D array. The 2D arrays are used in this s/r.
    call copy1Dimto2DimRealArray(ompTb, actualNumChannel, numObsToProcess, ompTb2D)
    call copy1Dimto2DimIntegerArray(obsChannels, actualNumChannel, numObsToProcess, obsChannels2D)
    call copy1Dimto2DimIntegerArray(obsFlags, actualNumChannel, numObsToProcess, obsFlags2D)
    ! allocations
    call utl_reAllocate(kchkprf, numObsToProcess)
    call utl_reAllocate(qcIndicator, actualNumChannel, numObsToProcess)
    !  Initialisations
    qcIndicator(:,:) = 0
    B7CHCK(:,:) = 0

    if ( RESETQC ) obsFlags2D(:,:) = 0

    ! 1) test 1: Check flag bit 7 on from the first bgckMwhs2 program
    !  Includes observations flagged for cloud liquid water, scattering index,
    !  dryness index plus failure of several QC checks.
    call atmsMwhs2Test1Flagbit7Check (itest, obsChannels2D, obsFlags2D, sensorIndex, qcIndicator, &
                                      actualNumChannel, numObsToProcess, stnId,  B7CHCK)

    ! 2) test 2: Topography check (partial)
    call atmsMwhs2Test2TopographyCheck (itest, obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                        modelInterpTerrain, obsFlags2D, ICHTOPO, ZCRIT, B7CHCK, qcIndicator)

    ! 3) test 3: Uncorrected Tb check (single)
    !  Uncorrected datum (flag bit #6 off). In this case switch bit 11 ON.
    call atmsMwhs2Test3UncorrectedTbCheck (itest, obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, &
                                           RESETQC, obsFlags2D, B7CHCK, qcIndicator)

    ! 4) test 4: "Rogue check" for (O-P) Tb residuals out of range. (single/full)
    !             Also, over WATER remove CH.10-15 if CH.10 |O-P|>5K (full)
    !  Les observations, dont le residu (O-P) depasse par un facteur (roguefac)
    !   l'erreur totale des TOVS.
    !  OVER OPEN WATER
    !    ch. 10 Abs(O-P) > 5K produces rejection of all ATMS amsub channels 10-15.
    call Mwhs2Test4RogueCheck (itest, obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, stnId, ROGUEFAC, &
                              waterobs, ompTb2D, cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                              newInformationFlag, ICH2OMPREJ, obsFlags2D, B7CHCK, qcIndicator)

    ! 5) test 5: Channel selection using array oer_tovutil(chan,sat)
    !  oer_tovutil = 0 (blacklisted)
    !                1 (assmilate)
    call atmsMwhs2Test5ChannelSelectionUsingTovutil(itest, obsChannels2D, sensorIndex, actualNumChannel, numObsToProcess, &
                                                    stnId, obsFlags2D, qcIndicator)

    !  Synthese de la controle de qualite au niveau de chaque point
    !  d'observation. Code:
    !            =0 aucun rejet, >0 au moins un canal rejete.
    do JJ = 1, numObsToProcess
      KCHKPRF(JJ) = 0
      do JI = 1, actualNumChannel
        KCHKPRF(JJ) = MAX(KCHKPRF(JJ),qcIndicator(JI,JJ))
      end do
    end do

    if ( mwbg_debug ) then
      write(*,*)'KCHKPRF = ',(KCHKPRF(JJ),JJ=1,numObsToProcess)
    end if

    ! Copy the modified FLAG to the 1D array, used outside this s/r.
    do JJ = 1, numObsToProcess
      do JI = 1, actualNumChannel
        INDX = (JJ-1) * actualNumChannel + JI
        obsFlags(INDX) = obsFlags2D(JI,JJ)
      end do
    end do

    ! reset global marker flag (55200) and mark it if observtions are rejected
    call resetQcCases(RESETQC, KCHKPRF, obsGlobalMarker)

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

      write(*,*) '   Meaning of newInformationFlag flag bits: '
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
      write(*,*) '   New Element 25174 in BURP file = newInformationFlag flag'
      write(*,*) ' '
    end if

  end subroutine mwbg_tovCheckMwhs2

  !--------------------------------------------------------------------------
  ! mwbg_readGeophysicFieldsAndInterpolate
  !--------------------------------------------------------------------------
  subroutine mwbg_readGeophysicFieldsAndInterpolate(instName, obsLat, obsLon, modelInterpTerrain, &
                                                    modelInterpGroundIce, modelInterpSeaIce)

    implicit none

    !:Purpose: Reads Modele Geophysical variables and save for the first time
    !         TOPOGRAPHIE (MF ou MX):
    !             MF est la topographie filtree avec unites en metres (filtered ME).
    !             MX est la topographie filtree avec unites en m2/s2  (geopotential topography).
    !         Glace de Mer (GL)
    !         Masque Terre-Mer (MG)
    !         Then Interpolate Those variables to observation location
    !Arguments:
    character(*),       intent(in)   :: instName                ! Instrument Name
    real,               intent(in)   :: obsLat(:)               ! Obseravtion Lats
    real,               intent(in)   :: obsLon(:)               ! Observation Lons
    real, allocatable,  intent(out)  :: modelInterpGroundIce(:) ! Glace de mer interpolees au pt d'obs.
    real, allocatable,  intent(out)  :: modelInterpTerrain(:)   ! topographie filtree (en metres) et interpolees
    real ,allocatable,  intent(out)  :: modelInterpSeaIce(:)    ! Glace de mer interpolees au pt d'obs.

    ! Locals:
    real, allocatable, save  :: GL(:)                           ! Modele Glace de Mer (GL)
    real, allocatable, save  :: MG(:)                           ! Modele Masque Terre-Mer (MG)
    real, allocatable, save  :: MT(:)                           ! Modele Topographie (MT)
    real,              save  :: TOPOFACT                        ! Facteur x topo pour avoir des unites en metre
    logical,           save  :: ifFirstCall = .True.            ! If .True. we read GL, MT and MG
    integer,           save  ::  gdmt                           ! topo interpolation param
    integer,           save  ::  gdmg                           ! mask terre-mer interpolation param
    integer,           save  ::  gdgl                           ! glace interpolation param
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

    ! STEP 0: CHECK if obsLat AND obsLon ARE SAME DIMENSION
    zlatNum = size(obsLat)
    zlonNum = size(obsLon)
    if (zlatNum /= zlonNum) then
      call utl_abort ('bgckMicrowave_mod: ERREUR: OBSERVATION obsLat and obsLon should have SAME LENGTH')
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
        XLAT = obsLat(dataIndex) +latIndex*DLAT
        XLAT = MAX(-90.0,MIN(90.0,XLAT))
        do lonIndex = -NLON, NLON
          boxPointIndex = boxPointIndex + 1
          XLON = obsLon(dataIndex) +lonIndex*DLON
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

    if(allocated(modelInterpTerrain)) deallocate(modelInterpTerrain)
    allocate (modelInterpTerrain(dataNum) , STAT=ier)
    if(allocated(modelInterpGroundIce)) deallocate(modelInterpGroundIce)
    allocate (modelInterpGroundIce(dataNum) , STAT=ier)
    if(allocated(modelInterpSeaIce)) deallocate(modelInterpSeaIce)
    allocate (modelInterpSeaIce(dataNum) , STAT=ier)
    do dataIndex = 1, dataNum
      if (mwbg_debug) then
        print *, ' ------------------  '
        print *, ' dataIndex = ', dataIndex
        print *, '   '
        print *, ' obsLat,obsLon = ', obsLat(dataIndex), obsLon(dataIndex)
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
      modelInterpGroundIce(dataIndex) = 0.0
      modelInterpTerrain(dataIndex) = 0.0
      modelInterpSeaIce(dataIndex) = 0.0
      do boxPointIndex=1,MXLAT*MXLON
        modelInterpTerrain(dataIndex) = MAX(modelInterpTerrain(dataIndex),MTINTBOX(boxPointIndex,dataIndex)/TOPOFACT)
        if(readGlaceMask) then
          modelInterpGroundIce(dataIndex) = MAX(modelInterpGroundIce(dataIndex),MGINTBOX(boxPointIndex,dataIndex))
          modelInterpSeaIce(dataIndex) = MAX(modelInterpSeaIce(dataIndex),GLINTBOX(boxPointIndex,dataIndex))
        end if
      end do
      if (mwbg_debug) then
        print *, ' modelInterpGroundIce = ', modelInterpGroundIce(dataIndex)
        print *, ' modelInterpTerrain = ', modelInterpTerrain(dataIndex)
        print *, ' modelInterpSeaIce = ', modelInterpSeaIce(dataIndex)
      end if
    end do
  end subroutine mwbg_readGeophysicFieldsAndInterpolate

  !--------------------------------------------------------------------------
  ! atmsMwhs2landIceMask
  !--------------------------------------------------------------------------
  subroutine atmsMwhs2landIceMask(numObsToProcess, obsLat, obsLon, calcLandQualifierIndice, &
                                  calcTerrainTypeIndice, waterobs, ilsmOpt)
    ! Adapted from: land_ice_mask_ssmis.ftn90 of mwbg_ssmis (D. Anselmo, S. Macpherson)
    !
    ! Object:   This routine sets waterobs array by performing a land/ice proximity check using
    !           using analysis MG and LG (or GL) fields used by the model which produces the trial field.
    !           The purpose of this check is to remove obs that reside close to coasts or ice,
    !           and so whose TBs may be contaminated.
    !           The GEM Global (glbhyb2) analysis contains MG and LG fields (on different grids).
    !
    !           NOTE: The 0.1 deg binary ice field check from land_ice_mask_ssmis.ftn90
    !           was removed. The land/sea qualifier (calcLandQualifierIndice) and terrain type (calcTerrainTypeIndice) 
    !           are modified to indicate proximity to land and sea-ice but are NOT changed in output BURP file.
    !
    !           In the application of this check, a 5x5 mesh, with spacing defined by rlat_km and
    !           rlon_km, is positioned with its center over an obs pt (2 grid pts on either side
    !           of the obs pt; size of mesh is equal to 4*rlat_km x 4*rlon_km). The values of MG
    !           and LG are evaluated at the grid points of this mesh. For ilsmOpt=1 (or 3), the maximum
    !           (or average) determines whether the obs pt is too close to ice or land to be retained.
    !           For ilsmOpt=2, the value at the central mesh point is used.
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
    ! - numObsToProcess : input  -  number of input obs pts in report
    ! - obsLat     : input  -  array holding lat values for all obs pts in report
    ! - obsLon     : input  -  array holding lon values for all obs pts in report
    ! - calcLandQualifierIndice  : in/out -  array holding land/sea qualifier values for all obs
    !                        pts of report (0 = land, 1 = sea)
    ! - calcTerrainTypeIndice  : in/out -  array holding terrain-type values for all obs pts
    !                        of current report (-1 land/open water, 0 = ice)
    ! - waterobs   : output -  logical array identifying for each obs in current report
    !                        whether it is over open water, far from coast/ice
    ! - ilsmOpt    : input  -  option for "interpolated" value of MG, LG at each location
    !                        1 = use MAX value taken from all mesh grid points
    !                        2 = use CENTRAL mesh point value (value at obs location)
    !                        3 = use AVG value of all mesh grid points    
    ! - mxlat      : internal-  number of grid pts in lat. direction for mesh
    ! - mxlon      : internal-  number of grid pts in lon. direction for mesh
    ! - rlat_km    : internal-  spacing desired between mesh grid points in km
    !                           along lat. direction
    ! - rlon_km    : internal-  spacing desired between mesh grid points in km
    !                           along lon. direction
    ! - dlat       : internal-  spacing between mesh grid points along lon. direction
    !                           in degrees computed from rlat_km
    ! - dlon       : internal-  spacing between mesh grid points along lon. direction
    !                           in degrees computed from rlon_km
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
    integer, intent(in)                   :: numObsToProcess
    integer, intent(in)                   :: ilsmOpt    
    real,    intent(in)                   :: obsLat(:)
    real,    intent(in)                   :: obsLon(:)
    integer, intent(out), allocatable     :: calcLandQualifierIndice(:)
    integer, intent(out), allocatable     :: calcTerrainTypeIndice(:)
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
    call utl_reAllocate(zlatbox, mxlat*mxlon, numObsToProcess)
    call utl_reAllocate(zlonbox, mxlat*mxlon, numObsToProcess)
    call utl_reAllocate(calcLandQualifierIndice, numObsToProcess)
    call utl_reAllocate(calcTerrainTypeIndice, numObsToProcess)
    call utl_reAllocate(waterobs, numObsToProcess)

    ! Open FST file.
    iungeo = 0
    ier = fnom( iungeo,fileMgLg,'STD+RND+R/O',0 )
    ier = fstouv( iungeo,'RND' )

    if (firstCall) then
      firstCall = .false.

      ! Read MG field.
      key = fstinf(iungeo,ni,nj,nk,-1,' ',-1,-1,-1,' ' ,'MG')
      if ( key <  0 ) then
        call utl_abort('atmsMwhs2landIceMask: The MG field is MISSING')
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
          call utl_abort('atmsMwhs2landIceMask: The LG or GL field is MISSING')
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
    do kk = 1, numObsToProcess
      indx = 0

      do ii = -nlat, nlat
        rii = float(ii)
        xlat = obsLat(kk) + rii*dlat
        xlat = max( -90.0, min(90.0,xlat) )
        xlatrad = xlat*pi/180.0

        do jj = -nlon, nlon
          dlon = rlon_km / ( rkm_per_deg*cos(xlatrad) )
          rjj = float(jj)
          indx = indx + 1
          xlon = obsLon(kk) + rjj*dlon
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

    call utl_reAllocate(mgintrp, numObsToProcess)
    call utl_reAllocate(lgintrp, numObsToProcess)

    mgintrp(:) = 0.0
    lgintrp(:) = 0.0
    do kk = 1, numObsToProcess

      latmesh = zlatbox(:,kk)
      lonmesh = zlonbox(:,kk)

      ier  = gdllsval(gdid,mgintob,mg,latmesh,lonmesh,mxlat*mxlon)
      ier  = gdllsval(gdidlg,lgintob,lg,latmesh,lonmesh,mxlat*mxlon)

      if ( ilsmOpt == 1 ) then
        mgintrp(kk) = maxval(mgintob(:))
        lgintrp(kk) = maxval(lgintob(:))
      elseif ( ilsmOpt == 2) then
        mgintrp(kk) = mgintob(ii_obsloc)
        lgintrp(kk) = lgintob(ii_obsloc)
      else
        mgintrp(kk) = sum(mgintob(:))/real((mxlat*mxlon))
        lgintrp(kk) = sum(lgintob(:))/real((mxlat*mxlon))
      end if      

    end do

    !  Initialize all obs as being over land and free of ice or snow.
    !  Determine which obs are over open water.
    waterobs(:) = .false.   ! not over open water
    calcTerrainTypeIndice(:) = -1             ! no ice (reset terain type)
    calcLandQualifierIndice(:) = 0 ! land   (reset land/sea qualifier)

    do kk = 1, numObsToProcess
      if ( mgintrp(kk) < MGthresh ) calcLandQualifierIndice(kk) = 1  ! ocean point away from coast
      if ( lgintrp(kk) >= LGthresh .and. calcLandQualifierIndice(kk) == 1 ) calcTerrainTypeIndice(kk) = 0  ! sea-ice affected point
      if ( lgintrp(kk)  < LGthresh .and. calcLandQualifierIndice(kk) == 1 ) then
        waterobs(kk) = .true.  ! water point not in close proximity to land or sea-ice
      end if
    end do

    ier = fstfrm(iungeo)
    ier = fclos(iungeo)

  end subroutine atmsMwhs2landIceMask

  !--------------------------------------------------------------------------
  ! mwbg_computeMwhs2SurfaceType
  !--------------------------------------------------------------------------

  subroutine mwbg_computeMwhs2SurfaceType(obsSpaceData)
    ! :Purpose: Compute surface type element and update obsSpaceData.

    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsSpaceData           ! ObsSpaceData object
    ! Locals:
    integer, allocatable :: calcLandQualifierIndice(:)
    integer, allocatable :: calcTerrainTypeIndice(:)
    logical, allocatable :: waterobs(:)
    integer              :: codtyp
    integer              :: headerIndex

    integer, parameter :: ilsmOpt = 2   ! OPTION for values of MG (land/sea mask) and LG (ice) 
                                        ! at each observation point using values on 5x5 mesh 
                                        ! centered at each point.
                                        ! ilsmOpt = 1 --> use MAX value from all 25 mesh points
                                        ! ilsmOpt = 2 --> use value at central mesh point (obs location)
                                        ! ilsmOpt = 3 --> use AVG value from all 25 mesh points
    logical              :: mwhs2DataPresent

    real                 :: obsLat(1)
    real                 :: obsLon(1)

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
      obsLat(1)  = obs_headElem_r( obsSpaceData, OBS_LAT, headerIndex )
      obsLat(1)  = obsLat(1) *MPC_DEGREES_PER_RADIAN_R8
      obsLon(1) = obs_headElem_r( obsSpaceData, OBS_LON, headerIndex )
      obsLon(1) = obsLon(1)*MPC_DEGREES_PER_RADIAN_R8
      if( obsLon(1) > 180. ) obsLon(1) = obsLon(1) - 360.
      call atmsMwhs2landIceMask(1, obsLat, obsLon, calcLandQualifierIndice, &
                                calcTerrainTypeIndice, waterobs, ilsmOpt)
      call obs_headSet_i(obsSpaceData, OBS_STYP, headerIndex, calcLandQualifierIndice(1))
      call obs_headSet_i(obsSpaceData, OBS_TTYP, headerIndex, calcTerrainTypeIndice(1))

    end do HEADER1

    write(*,*) 'ssbg_computeMwhs2SurfaceType: Finished'

  end subroutine mwbg_computeMwhs2SurfaceType

  !--------------------------------------------------------------------------
  ! mwbg_grossValueCheck  
  !--------------------------------------------------------------------------

  subroutine mwbg_grossValueCheck(numObsToProcess, actualNumChannel, obsTb, obsTbBiasCorr, ztbThresholdMin, ztbThresholdMax, grossrej)

    !:Purpose: Check Tbs for values that are missing or outside physical limits.
    !          **NOTE: REJECT ALL CHANNELS OF ONE IS FOUND TO BE BAD.
    implicit none

    ! Arguments
    integer, intent(in)               :: numObsToProcess  ! number of obs pts to process
    integer, intent(in)               :: actualNumChannel ! number of ichannels
    real,    intent(in)               :: obsTb(:)         ! radiances
    real,    intent(in)               :: obsTbBiasCorr(:) ! Bias correction
    real,    intent(in)               :: ztbThresholdMin  ! ztb threshold for rejection
    real,    intent(in)               :: ztbThresholdMax  ! ztb threshold for rejection
    logical, intent(out), allocatable :: grossrej(:)      ! logical array defining which obs are to be rejected

    ! Locals
    integer :: pt, indx, channel

    real, allocatable                 :: ztb(:)           ! biased or unbiased radiances

    call utl_reAllocate(ztb, actualNumChannel)
    call utl_reAllocate(grossrej, numObsToProcess)
    
    grossrej(1:numObsToProcess) = .true.
    indx = 0
    do pt = 1, numObsToProcess
      if ( mwbg_useUnbiasedObsForClw ) then
        do channel = 1, actualNumChannel
          ztb(channel) = obsTb(indx+channel)
        end do
      else
        do channel = 1, actualNumChannel
          if (obsTbBiasCorr(indx+channel) /= mwbg_realMissing) then
            ztb(channel) = obsTb(indx+channel) - obsTbBiasCorr(indx+channel)
          else
            ztb(channel) = obsTb(indx+channel)
          end if
        end do
      end if
      if ( all( ztb > ztbThresholdMin ) .and. all( ztb < ztbThresholdMax ) ) then
        grossrej(pt) = .false.
      end if
      indx = pt*actualNumChannel

    end do

  end subroutine mwbg_grossValueCheck

  !--------------------------------------------------------------------------
  ! mwbg_firstQcCheckAtms
  !--------------------------------------------------------------------------
  subroutine mwbg_firstQcCheckAtms(satZenithAngle, landQualifierIndice, terrainTypeIndice, obsLat, obsLon, obsTb, satScanPosition, &
                                   actualNumChannel, numObsToProcess, qcRejectLogic, grossrej, calcLandQualifierIndice, calcTerrainTypeIndice, &
                                   obsQcFlag1, obsQcFlag2, obsChannels, reportHasMissingTb)
    !  :Purpose: This routine performs basic quality control checks on the data. It sets array
    !            qcRejectLogic(numObsToProcess,actualNumChannel) elements to .true. to flag data with failed checks.
    !
    !  The 6 QC checks are:
    !                 - 1) Invalid land/sea qualifier or terrain type,
    !                 - 2) Invalid field of view number,
    !                 - 3) Satellite zenith angle missing or out of range, (> 75 deg),
    !                 - 4) lat,lon check (lat,lon = O(-90.), 0(-180.))
    !                 - 5) Change in (computed) calcLandQualifierIndice,calcTerrainTypeIndice from (input) 
    !                      landQualifierIndice,terrainTypeIndice (from MG,LG fields).
    !                      landQualifierIndice= 0,1 (from hi-res land/sea mask interpolated to obs point [CMDA])
    !                      terrainTypeIndice=-1,0 (from hi-res ice analysis  interpolated to obs point [CMDA])
    !                      calcLandQualifierIndice= 0,1 (from max interp MG (0.0 to 1.0) in box surrounding obs point)
    !                      calcTerrainTypeIndice=-1,0 (from max interp LG (0.0 to 1.0) in box surrounding obs point)
    !                 - 6) ATMS quality flag check (qual. flag elements 33078,33079,33080,33081)

    !
    !  In most cases, qcRejectLogic(ii,actualNumChannel) is set to .true. for all channels at point ii
    !  if the check detects a problem. In addition, Tb (obsTb) is set to missing_value
    !  for checks 3 and 4 fails.
    implicit none

    ! Arguments
    integer,              intent(in)                :: landQualifierIndice(:)
    integer,              intent(in)                :: terrainTypeIndice(:)
    integer,              intent(in)                :: satScanPosition(:)
    integer,              intent(in)                :: obsChannels(:)
    integer,              intent(in)                :: obsQcFlag2(:)
    integer,              intent(in)                :: obsQcFlag1(:,:)
    integer,              intent(in)                :: numObsToProcess
    integer,              intent(in)                :: actualNumChannel
    integer,              intent(in)                :: calcLandQualifierIndice(:)
    integer,              intent(in)                :: calcTerrainTypeIndice(:)
    logical,              intent(in)                :: grossrej(:)        ! dim(numObsToProcess), true if 1 or more Tb fail gross error check
    real,                 intent(in)                :: obsLat(:)
    real,                 intent(in)                :: obsLon(:)
    real,                 intent(inout)             :: obsTb(:)
    real,                 intent(inout)             :: satZenithAngle(:)
    logical,              intent(out)               :: reportHasMissingTb ! true if Tb(obsTb) are set to missing_value
    logical, allocatable, intent(out)               :: qcRejectLogic(:,:) ! dim(numObsToProcess,actualNumChannel) 
                                                                          ! qcRejectLogic = .false. on input

    ! Locals
    integer :: ii, jj, indx1, icount
    logical :: fail, fail1, fail2

    reportHasMissingTb = .false.
    call utl_reAllocate(qcRejectLogic, numObsToProcess, actualNumChannel)
    qcRejectLogic(:,:) = .false.  ! Flag for preliminary QC checks
    ! Global rejection checks

    ! Check if number of channels is correct
    !if ( actualNumChannel /= mwbg_maxNumChan ) then
    !  write(*,*) 'WARNING: Number of channels (',actualNumChannel, ') is not equal to mwbg_maxNumChan (', mwbg_maxNumChan,')'
    !  write(*,*) '         All data flagged as bad and returning to calling routine!'
    !  qcRejectLogic(:,:) = .true.  ! flag all data in report as bad
    !  return
    !end if

    ! Check for errors in channel numbers (should be 1-22 for each location ii)
    indx1 = 1
    fail = .false.
    do ii = 1,numObsToProcess
      do jj = 1, actualNumChannel
        if ( obsChannels(indx1+jj-1) /= jj ) fail = .true.
      end do
      indx1 = indx1 + actualNumChannel
    end do
    if ( fail ) then
      write(*,*) 'WARNING: Bad channel number(s) detected!'
      write(*,*) '         All data flagged as bad and returning to calling routine!'
      write(*,*) '  obsChannels(numObsToProcess*actualNumChannel) array = ', obsChannels(:)
      qcRejectLogic(:,:) = .true.  ! flag all data in report as bad
      return
    end if

    ! 1) invalid land/sea qualifier or terrain type
    !  landQualifierIndice = 0 (land),     1 (sea)
    !  terrainTypeIndice = 0 (sea-ice), -1 otherwise
    !  calcLandQualifierIndice = 1 (sea, away from land/coast [MG]),      0 otherwise
    !  calcTerrainTypeIndice = 0 (over or near analyzed sea-ice [LG]), -1 otherwise
    do ii = 1,numObsToProcess
      fail = .false.
      if ( landQualifierIndice(ii) < 0  .or. landQualifierIndice(ii) > 2 ) fail = .true.
      if ( terrainTypeIndice(ii) < -1 .or. terrainTypeIndice(ii) > 1 ) fail = .true.
      if ( fail ) then
        write(*,*) 'WARNING: Invalid land/sea qualifier or terrain type!'
        write(*,*) '  landQualifierIndice, terrainTypeIndice, (lat, lon) = ', landQualifierIndice(ii), terrainTypeIndice(ii), '(',obsLat(ii), obsLon(ii),')'
      end if

      if ( landQualifierIndice(ii) == 0 .and. terrainTypeIndice(ii) == 0 ) then
        fail = .true.
        write(*,*) 'WARNING: Sea ice point (terrainTypeIndice=0) at land point (landQualifierIndice=0)!'
        write(*,*) ' lat, lon =  ', obsLat(ii), obsLon(ii)
      end if
      if ( fail ) qcRejectLogic(ii,:) = .true.
    end do

    do ii = 1,numObsToProcess
      fail = .false.
      if ( calcLandQualifierIndice(ii) < 0  .or. calcLandQualifierIndice(ii) > 2 ) fail = .true.
      if ( calcTerrainTypeIndice(ii) < -1 .or. calcTerrainTypeIndice(ii) > 1 ) fail = .true.
      if ( fail ) then
        write(*,*) 'WARNING: Invalid model-based (MG/LG) land/sea qualifier or terrain type!'
        write(*,*) '  calcLandQualifierIndice, calcTerrainTypeIndice, (lat, lon) = ', &
                   calcLandQualifierIndice(ii), calcTerrainTypeIndice(ii), '(',obsLat(ii), obsLon(ii),')'
      end if
      if ( fail ) qcRejectLogic(ii,:) = .true.
    end do

    ! 2) invalid field of view number
    do ii = 1,numObsToProcess
      fail = .false.
      if ( satScanPosition(ii) < 1  .or. satScanPosition(ii) > mwbg_maxScanAngle ) then
        fail = .true.
        write(*,*) 'WARNING: Invalid field of view! satScanPosition, lat, lon = ', satScanPosition(ii), obsLat(ii), obsLon(ii)
      end if
      if ( fail ) qcRejectLogic(ii,:) = .true.
    end do

    ! 3) satellite zenith angle missing or out of range (> 75 deg)
    !  If bad zenith, then set Tb (and zenith) = missing value
    indx1 = 1
    do ii = 1,numObsToProcess
      fail = .false.
      if ( satZenithAngle(ii) > 75.0 .or. satZenithAngle(ii) < 0. ) then
        fail = .true.
        write(*,*) 'WARNING: Bad or missing zenith angle! zenith, lat, lon = ', satZenithAngle(ii), obsLat(ii), obsLon(ii)
        satZenithAngle(ii) = mwbg_realMissing
        reportHasMissingTb = .true.
      end if
      do jj = 1, actualNumChannel
        if ( fail ) then
          qcRejectLogic(ii,jj) = .true.
          obsTb(indx1+jj-1) = mwbg_realMissing
        end if
      end do
      indx1 = indx1 + actualNumChannel
    end do

    ! 4) Lat,lon check
    ! Check for undecoded BURP file integer values of lat,lon = 0,0
    ! (usually associated with missing zenith angle and erroneous Tb=330K)

    icount = 0
    indx1 = 1
    do ii = 1,numObsToProcess
      fail = .false.
      if ( obsLat(ii) == -90.0  .and. obsLon(ii) == -180.0 ) then
        fail = .true.
        icount =  icount + 1
        reportHasMissingTb = .true.
      end if
      do jj = 1, actualNumChannel
        if ( fail ) then
          qcRejectLogic(ii,jj) = .true.
          obsTb(indx1+jj-1) = mwbg_realMissing
        end if
      end do
      indx1 = indx1 + actualNumChannel
    end do
    if ( icount > 0 ) write(*,*) 'WARNING: Bad lat,lon pair(s) detected. Number of locations = ', icount

    icount = 0
    indx1 = 1
    do ii = 1,numObsToProcess
      fail = .false.
      if ( abs(obsLat(ii)) > 90.0  .or. abs(obsLon(ii)) > 180.0 ) then
        fail = .true.
        icount =  icount + 1
        reportHasMissingTb = .true.
      end if
      do jj = 1, actualNumChannel
        if ( fail ) then
          qcRejectLogic(ii,jj) = .true.
          obsTb(indx1+jj-1) = mwbg_realMissing
        end if
      end do
      indx1 = indx1 + actualNumChannel
    end do
    if ( icount > 0 ) write(*,*) 'WARNING: Lat or lon out of range! Number of locations = ', icount

    !  5) Change in land/sea qualifier or terrain-type based on MG,LG fields
    icount = 0
    do ii = 1,numObsToProcess
      fail = .false.
      if ( (landQualifierIndice(ii) /= calcLandQualifierIndice(ii)) .or. &
           (terrainTypeIndice(ii) /= calcTerrainTypeIndice(ii)) ) then
        fail = .true.
      end if
      if ( fail ) then
        icount =  icount + 1
      end if
    end do

    ! 6) ATMS quality flag check (qual. flag elements 33078,33079,33080,33081)

    !  33078 Geolocation quality code     obsQcFlag1(ii,1)  code value = 0-15 (0= OK, 15=misg)
    !  33079 Granule level quality flags  obsQcFlag1(ii,2)  16 bit flag  (start bit 6(2^5)=32) (misg=2^16-1 = 65535)
    !  33080 Scan level quality flags     obsQcFlag1(ii,3)  20 bit flag  (start bit 7(2^6)=64) (misg=2^20-1) 
    !  33081 Channel data quality flags   obsQcFlag2        12 bit flag  (start bit 3(2^2)=4)  (misg=2^12-1)
    !
    !  See http://www.wmo.int/pages/prog/www/WMOCodes/WMO306_vI2/2010edition/BUFRver16/BUFR_16_0_0_TableD.pdf

    indx1 = 1
    do ii = 1,numObsToProcess
      fail1 = .false.
      fail = .false.
      if ( (obsQcFlag1(ii,1) > 0) .or. (obsQcFlag1(ii,2) >= 32) .or. (obsQcFlag1(ii,3) >= 64) ) then
        write(*,*) 'WARNING: INFO BLOCK QC flag(s) indicate problem with data'
        write(*,*) ' ele33078 = ',obsQcFlag1(ii,1),' ele33079 = ',obsQcFlag1(ii,2),' ele33080 = ', obsQcFlag1(ii,3)
        write(*,*) ' lat, lon = ', obsLat(ii), obsLon(ii)
        fail1 = .true.
        if ( grossrej(ii) ) write(*,*) ' NOTE: grossrej is also true for this point!'
      end if
      do jj = 1, actualNumChannel
        fail2 = .false.
        if ( obsQcFlag2(indx1+jj-1) >= 4 ) then
          !write(*,*) 'WARNING: DATA BLOCK QC flag ele33081 = ', obsQcFlag2(indx1+jj-1)
          !write(*,*) '    Lat, lon, channel = ', obsLat(ii), obsLon(ii), obsChannels(indx1+jj-1)
          fail2 = .true.
          fail = .true.
          !if ( (.not. fail1) .and. grossrej(ii) ) write(*,*) ' NOTE: grossrej is also true for this point!'
        end if
        if ( fail2 .or. fail1 ) qcRejectLogic(ii,jj) = .true.
      end do
      if ( fail ) write(*,*) 'WARNING: DATA BLOCK QC flag ele33081 >= 4 for one or more channels! lat, lon = ', obsLat(ii), obsLon(ii)
      indx1 = indx1 + actualNumChannel
    end do

    !write(*,*) 'mwbg_firstQcCheckAtms: Number of data processed and flagged = ', &
    !           numObsToProcess*actualNumChannel, count(qcRejectLogic)

  end subroutine mwbg_firstQcCheckAtms

  !--------------------------------------------------------------------------
  ! mwbg_firstQcCheckMwhs2
  !--------------------------------------------------------------------------
  subroutine mwbg_firstQcCheckMwhs2(satZenithAngle, landQualifierIndice, terrainTypeIndice, obsLat, obsLon, obsTb, satScanPosition, &
                                    actualNumChannel, numObsToProcess, qcRejectLogic, calcLandQualifierIndice, calcTerrainTypeIndice, &
                                    obsChannels, reportHasMissingTb, modLSQ)
    !  :Purpose: This routine performs basic quality control checks on the data. It sets array
    !            qcRejectLogic(numObsToProcess,actualNumChannel) elements to .true. to flag data with failed checks. Check 1
    !            (for landQualifierIndice,terrainTypeIndice) and check 5 are skipped if modlsqtt=.true., 
    !            as the original values will be replaced in output file by calcLandQualifierIndice,calcTerrainTypeIndice.
    !
    !  The 5 QC checks are:
    !                 - 1) Invalid land/sea qualifier or terrain type,
    !                 - 2) Invalid field of view number,
    !                 - 3) Satellite zenith angle missing or out of range, (> 75 deg),
    !                 - 4) lat,lon check (lat,lon = O(-90.), 0(-180.))
    !                 - 5) Change in (computed) calcLandQualifierIndice,calcTerrainTypeIndice 
    !                      from (input) landQualifierIndice,terrainTypeIndice (from MG,LG fields).
    !                      landQualifierIndice= 0,1 (from hi-res land/sea mask interpolated to obs point [CMDA])
    !                      terrainTypeIndice=-1,0 (from hi-res ice analysis  interpolated to obs point [CMDA])
    !                      calcLandQualifierIndice= 0,1 (from max interp MG (0.0 to 1.0) in box surrounding obs point)
    !                      calcTerrainTypeIndice=-1,0 (from max interp LG (0.0 to 1.0) in box surrounding obs point)

    !
    !  In most cases, qcRejectLogic(ii,actualNumChannel) is set to .true. for all channels at point ii
    !  if the check detects a problem. In addition, Tb (obsTb) is set to missing_value
    !  for checks 3 and 4 fails.
    implicit none

    ! Arguments
    integer,              intent(in)                :: landQualifierIndice(:)
    integer,              intent(in)                :: terrainTypeIndice(:)
    integer,              intent(in)                :: satScanPosition(:)
    integer,              intent(in)                :: obsChannels(:)
    integer,              intent(in)                :: numObsToProcess
    integer,              intent(in)                :: actualNumChannel
    integer,              intent(in)                :: calcLandQualifierIndice(:)
    integer,              intent(in)                :: calcTerrainTypeIndice(:)
    real,                 intent(in)                :: obsLat(:)
    real,                 intent(in)                :: obsLon(:)
    real,                 intent(inout)             :: obsTb(:)
    real,                 intent(inout)             :: satZenithAngle(:)
    logical,              intent(out)               :: reportHasMissingTb ! true if Tb(obsTb) are set to missing_value
    logical,              intent(in)                :: modLSQ
    logical, allocatable, intent(out)               :: qcRejectLogic(:,:) ! dim(numObsToProcess,actualNumChannel)
                                                                          ! qcRejectLogic = .false. on input

    ! Locals
    integer :: ii, jj, indx1, icount
    logical :: fail

    reportHasMissingTb = .false.
    call utl_reAllocate(qcRejectLogic, numObsToProcess, actualNumChannel)
    qcRejectLogic(:,:) = .false.  ! Flag for preliminary QC checks
    ! Global rejection checks

    ! Check if number of channels is correct
    !if ( actualNumChannel /= mwbg_maxNumChan ) then
    !  write(*,*) 'WARNING: Number of channels (',actualNumChannel, ') is not equal to mwbg_maxNumChan (', mwbg_maxNumChan,')'
    !  write(*,*) '         All data flagged as bad and returning to calling routine!'
    !  qcRejectLogic(:,:) = .true.  ! flag all data in report as bad
    !  return
    !end if

    ! Check for errors in channel numbers (should be 1-15 for each location ii)
    indx1 = 1
    fail = .false.
    do ii = 1,numObsToProcess
      do jj = 1, actualNumChannel
        if ( obsChannels(indx1+jj-1) /= jj ) fail = .true.
      end do
      indx1 = indx1 + actualNumChannel
    end do
    if ( fail ) then
      write(*,*) 'WARNING: Bad channel number(s) detected!'
      write(*,*) '         All data flagged as bad and returning to calling routine!'
      write(*,*) '  obsChannels(numObsToProcess*actualNumChannel) array = ', obsChannels(:)
      qcRejectLogic(:,:) = .true.  ! flag all data in report as bad
      return
    end if

    ! 1) invalid land/sea qualifier or terrain type
    !  landQualifierIndice = 0 (land),     1 (sea)
    !  terrainTypeIndice = 0 (sea-ice), -1 otherwise
    !  calcLandQualifierIndice = 1 (sea, away from land/coast [MG]),      0 otherwise
    !  calcTerrainTypeIndice = 0 (over or near analyzed sea-ice [LG]), -1 otherwise

    ! Checks on landQualifierIndice,terrainTypeIndice are not done if values are to be replaced in output file.

    if ( .not. modLSQ ) then
      do ii = 1,numObsToProcess
        fail = .false.
        if ( landQualifierIndice(ii) < 0  .or. landQualifierIndice(ii) > 2 ) fail = .true.
        if ( terrainTypeIndice(ii) < -1 .or. terrainTypeIndice(ii) > 1 ) fail = .true.
        if ( fail ) then
          write(*,*) 'WARNING: Invalid land/sea qualifier or terrain type!'
          write(*,*) '  landQualifierIndice, terrainTypeIndice, (lat, lon) = ', landQualifierIndice(ii), terrainTypeIndice(ii), '(',obsLat(ii), obsLon(ii),')'
        end if

        if ( landQualifierIndice(ii) == 0 .and. terrainTypeIndice(ii) == 0 ) then
          fail = .true.
          write(*,*) 'WARNING: Sea ice point (terrainTypeIndice=0) at land point (landQualifierIndice=0)!'
          write(*,*) ' lat, lon =  ', obsLat(ii), obsLon(ii)
        end if
        if ( fail ) qcRejectLogic(ii,:) = .true.
      end do
    end if

    do ii = 1,numObsToProcess
      fail = .false.
      if ( calcLandQualifierIndice(ii) < 0  .or. calcLandQualifierIndice(ii) > 2 ) fail = .true.
      if ( calcTerrainTypeIndice(ii) < -1 .or. calcTerrainTypeIndice(ii) > 1 ) fail = .true.
      if ( fail ) then
        write(*,*) 'WARNING: Invalid model-based (MG/LG) land/sea qualifier or terrain type!'
        write(*,*) '  calcLandQualifierIndice, calcTerrainTypeIndice, (lat, lon) = ', &
                   calcLandQualifierIndice(ii), calcTerrainTypeIndice(ii), '(',obsLat(ii), obsLon(ii),')'
      end if
      if ( fail ) qcRejectLogic(ii,:) = .true.
    end do

    ! 2) invalid field of view number
    do ii = 1,numObsToProcess
      fail = .false.
      if ( satScanPosition(ii) < 1  .or. satScanPosition(ii) > mwbg_maxScanAngle ) then
        fail = .true.
        write(*,*) 'WARNING: Invalid field of view! satScanPosition, lat, lon = ', satScanPosition(ii), obsLat(ii), obsLon(ii)
      end if
      if ( fail ) qcRejectLogic(ii,:) = .true.
    end do

    ! 3) satellite zenith angle missing or out of range (> 75 deg)
    !  If bad zenith, then set Tb (and zenith) = missing value
    indx1 = 1
    do ii = 1,numObsToProcess
      fail = .false.
      if ( satZenithAngle(ii) > 75.0 .or. satZenithAngle(ii) < 0. ) then
        fail = .true.
        write(*,*) 'WARNING: Bad or missing zenith angle! zenith, lat, lon = ', satZenithAngle(ii), obsLat(ii), obsLon(ii)
        satZenithAngle(ii) = mwbg_realMissing
        reportHasMissingTb = .true.
      end if
      do jj = 1, actualNumChannel
        if ( fail ) then
          qcRejectLogic(ii,jj) = .true.
          obsTb(indx1+jj-1) = mwbg_realMissing
        end if
      end do
      indx1 = indx1 + actualNumChannel
    end do

    ! 4) Lat,lon check
    ! Check for undecoded BURP file integer values of lat,lon = 0,0
    ! (usually associated with missing zenith angle and erroneous Tb=330K)

    icount = 0
    indx1 = 1
    do ii = 1,numObsToProcess
      fail = .false.
      if ( obsLat(ii) == -90.0  .and. obsLon(ii) == -180.0 ) then
        fail = .true.
        icount =  icount + 1
        reportHasMissingTb = .true.
      end if
      do jj = 1, actualNumChannel
        if ( fail ) then
          qcRejectLogic(ii,jj) = .true.
          obsTb(indx1+jj-1) = mwbg_realMissing
        end if
      end do
      indx1 = indx1 + actualNumChannel
    end do
    if ( icount > 0 ) write(*,*) 'WARNING: Bad lat,lon pair(s) detected. Number of locations = ', icount

    icount = 0
    indx1 = 1
    do ii = 1,numObsToProcess
      fail = .false.
      if ( abs(obsLat(ii)) > 90.0  .or. abs(obsLon(ii)) > 180.0 ) then
        fail = .true.
        icount =  icount + 1
        reportHasMissingTb = .true.
      end if
      do jj = 1, actualNumChannel
        if ( fail ) then
          qcRejectLogic(ii,jj) = .true.
          obsTb(indx1+jj-1) = mwbg_realMissing
        end if
      end do
      indx1 = indx1 + actualNumChannel
    end do
    if ( icount > 0 ) write(*,*) 'WARNING: Lat or lon out of range! Number of locations = ', icount

    !  5) Change in land/sea qualifier or terrain-type based on MG,LG fields
    if ( .not. modLSQ ) then
      icount = 0
      do ii = 1,numObsToProcess
        fail = .false.
        if ( (landQualifierIndice(ii) /= calcLandQualifierIndice(ii)) .or. (terrainTypeIndice(ii) /= calcTerrainTypeIndice(ii)) ) then
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
  subroutine mwbg_nrlFilterAtms(numObsToProcess, actualNumChannel, obsTb, ompTb, obsTbBiasCorr, satZenithAngle, obsLat, &
                                calcLandQualifierIndice, calcTerrainTypeIndice, waterobs, grossrej, &
                                cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                                si_ecmwf, si_bg, iNumSeaIce, iRej, SeaIce)
    !OBJET          Compute the following parameters using 5 ATMS channels:
    !                  - sea ice,
    !                  - cloud liquid water from observation (cloudLiquidWaterPathObs),
    !                  - cloud liquid water from first guess (cloudLiquidWaterPathFG),
    !                  - 2 scattering indices (si) (ECMWF, Bennartz-Grody)
    !               The five channels used are: 23Ghz, 31Ghz, 50Ghz, 89Ghz, and 165Ghz.
    !
    !NOTES*
    !                - open water points are converted to sea-ice points if sea ice concentration >= 0.55
    !                   and calcTerrainTypeIndice (terrainTypeIndice or terrain type) is changed accordingly
    !                - cloudLiquidWaterPathObs are missing when out-of-range parameters/Tb detected or grossrej = .true.
    !                - cloudLiquidWaterPathObs and si only computed over open water away from coasts and sea-ice
    !                - cloudLiquidWaterPathObs and si = -99.0 where value cannot be computed.
    !
    !REFERENCES     Ben Ruston, NRL Monterey
    !                  JCSDA Seminar 12/12/12: Impact of NPP Satellite Assimilation in the U.S. Navy Global Modeling System
    !
    !
    !ARGUMENTS      - ier         - output - error return code for each location:
    !                                        0, ok,
    !                                        1, input parameter out of range or grossrej=.true.
    !               - numObsToProcess - input  -  number of points to process (= numObsToProcess)
    !               - tb23            - input  -  23Ghz brightness temperature (K) -- ch. 1
    !               - tb23FG          - input  -  23Ghz brightness temperature (K) from first guess -- ch. 1
    !               - tb31            - input  -  31Ghz brightness temperature (K) -- ch. 2
    !               - tb31FG          - input  -  31Ghz brightness temperature (K) from first guess -- ch. 2
    !               - tb50            - input  -  50Ghz brightness temperature (K) -- ch. 3
    !               - tb89            - input  -  89Ghz brightness temperature (K) -- ch. 16
    !               - tb165           - input  -  165Ghz brightness temperature (K) -- ch. 17
    !               - satZenithAngle  - input  -  satellite zenith angle (deg.)
    !               - obsLat          - input  -  latitude (deg.)
    !               - calcLandQualifierIndice - input  -  land/sea indicator (0=land, 1=ocean)
    !               - calcTerrainTypeIndice   - in/out -  terrain type (0=ice, -1 otherwise)
    !               - waterobs    - in/out -  .true. if open water point (away from coasts and sea-ice)
    !               - grossrej    - input  -  .true. if any channel had a gross error from mwbg_grossValueCheck
    !               - cloudLiquidWaterPathObs - output -  cloud liquid water from observation (kg/m**2) from tb23 & tb31
    !               - cloudLiquidWaterPathFG  - output -  cloud liquid water from first guess (kg/m**2) from tb23FG & tb31FG
    !               - si_ecmwf    - output -  ECMWF scattering index from tb89 & tb165
    !               - si_bg       - output -  Bennartz-Grody scattering index from tb89 & tb165
    !               - iNumSeaIce  - in/out -  running counter for number of open water points
    !                                         with sea-ice detected (from algorithm)
    !               - iRej        - in/out -  running counter for number of locations with bad
    !                                         satZenithAngle, obsLat, calcLandQualifierIndice, or with grossrej=true
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

    integer, intent(in)                   ::  numObsToProcess
    integer, intent(in)                   ::  actualNumChannel
    integer, intent(out)                  ::  iNumSeaIce
    integer, intent(in)                   ::  calcLandQualifierIndice(:)
    integer, intent(inout)                ::  calcTerrainTypeIndice(:)
    integer, intent(out)                  ::  iRej


    logical, intent(in)                   ::  grossrej(:)
    logical, intent(inout)                ::  waterobs(:)

    real, intent(in)                      ::  obsTb(:)
    real, intent(in)                      ::  ompTb(:)
    real, intent(in)                      ::  obsTbBiasCorr(:)
    real, intent(in)                      ::  satZenithAngle(:)
    real, intent(in)                      ::  obsLat(:)
    real, allocatable, intent(out)        ::  cloudLiquidWaterPathObs(:)
    real, allocatable, intent(out)        ::  cloudLiquidWaterPathFG(:)
    real, allocatable, intent(out)        ::  si_ecmwf(:)
    real, allocatable, intent(out)        ::  si_bg(:)
    real, allocatable, intent(out)        ::  SeaIce(:)

    ! Locals
    integer                               :: ier(numObsToProcess)
    real                                  ::  ice(numObsToProcess)
    real                                  :: tb23(numObsToProcess)
    real                                  :: tb23FG(numObsToProcess)
    real                                  :: tb31(numObsToProcess)
    real                                  :: tb31FG(numObsToProcess)
    real                                  :: tb50(numObsToProcess)
    real                                  :: tb89(numObsToProcess)
    real                                  :: tb165(numObsToProcess)
    real                                  :: bcor23(numObsToProcess)
    real                                  :: bcor31(numObsToProcess)
    real                                  :: bcor50(numObsToProcess)
    real                                  :: bcor89(numObsToProcess)
    real                                  :: bcor165(numObsToProcess)
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
    call utl_reAllocate(cloudLiquidWaterPathObs,numObsToProcess)
    call utl_reAllocate(cloudLiquidWaterPathFG,numObsToProcess)
    call utl_reAllocate(si_ecmwf,numObsToProcess)
    call utl_reAllocate(si_bg,numObsToProcess)
    call utl_reAllocate(SeaIce,numObsToProcess)

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
    do ii = 1, numObsToProcess
      indx2 = ii*actualNumChannel
      tb23(ii)      = obsTb(indx1)
      tb23FG(ii)    = obsTb(indx1) - ompTb(indx1)
      bcor23(ii)    = obsTbBiasCorr(indx1)
      tb31(ii)      = obsTb(indx1+1)
      tb31FG(ii)    = obsTb(indx1+1) - ompTb(indx1+1)
      bcor31(ii)    = obsTbBiasCorr(indx1+1)
      tb50(ii)      = obsTb(indx1+2)
      bcor50(ii)    = obsTbBiasCorr(indx1+2)
      tb89(ii)      = obsTb(indx1+15)
      bcor89(ii)    = obsTbBiasCorr(indx1+15)
      tb165(ii)    = obsTb(indx1+16)
      bcor165(ii)    = obsTbBiasCorr(indx1+16)
      indx1 = indx2 + 1
    end do

    ier = 0

    ! 1) Initialise parameters:
    do i = 1, numObsToProcess
      ice(i)      = mwbg_realMissing
      cloudLiquidWaterPathObs(i)   = mwbg_realMissing
      cloudLiquidWaterPathFG(i)    = mwbg_realMissing
      si_ecmwf(i) = mwbg_realMissing
      si_bg(i)    = mwbg_realMissing
      SeaIce(i)   = 0.0
    end do

    ! 2) Validate input parameters:
    do i = 1, numObsToProcess
      if ( satZenithAngle(i) < 0.  .or. &
           satZenithAngle(i) > 70. .or. &
           obsLat(i)  < -90.  .or. &
           obsLat(i)  >  90.  .or. &
           calcLandQualifierIndice(i) < 0 .or. &
           calcLandQualifierIndice(i) > 1 ) then
         ier(i) = 1
      end if

      ! Skip computations for points where all data are rejected  (bad Tb ANY channel)
      if ( grossrej(i) ) then
        ier(i) = 1
      end if
    end do

    ! 3) Compute parameters:
    do i = 1, numObsToProcess

      if ( ier(i) == 0 ) then

        abslat = abs(obsLat(i))
        cosz   = cosd(satZenithAngle(i))

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
        if ( calcLandQualifierIndice(i) == 1 ) then  ! water point

          if ( abslat < 50. ) then
            ice(i) = 0.0
          else
            ice(i) = 2.85 + 0.020*t23 - 0.028*t50
          end if

          SeaIce(i) = ice(i)

          if ( ice(i) >= 0.55 .and. waterobs(i) ) then
            iNumSeaIce = iNumSeaIce + 1
            waterobs(i) = .false.
            calcTerrainTypeIndice(i) = 0
          end if

        end if

        ! Compute cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, and Scattering Indices (over open water only)
        if ( waterobs(i) ) then
          if ( t23 < 284. .and. t31 < 284. ) then
            aa = 8.24 - (2.622 - 1.846 * cosz) * cosz
            cloudLiquidWaterPathObs(i) = aa + 0.754 * alog(285.0 - t23) - 2.265 * alog(285.0 - t31)
            cloudLiquidWaterPathObs(i) = cloudLiquidWaterPathObs(i) * cosz
            if ( cloudLiquidWaterPathObs(i) < 0.0 ) cloudLiquidWaterPathObs(i) = 0.0

            cloudLiquidWaterPathFG(i) = aa + 0.754 * alog(285.0 - t23FG) - 2.265 * alog(285.0 - t31FG)
            cloudLiquidWaterPathFG(i) = cloudLiquidWaterPathFG(i) * cosz
            if ( cloudLiquidWaterPathFG(i) < 0.0 ) cloudLiquidWaterPathFG(i) = 0.0
          end if
          si_ecmwf(i) = deltb - (-46.94 + 0.248 * satZenithAngle(i))
          si_bg(i)    = deltb - (-39.201 + 0.1104 * satZenithAngle(i))
        end if

      else  ! ier(i) == 1 case
         iRej = iRej + 1

      end if ! if ( ier(i) == 0 )

      if ( mwbg_debug .and. (i <= 100) ) then
        write(*,*) ' '
        write(*,*) ' i,tb23(i),tb23FG(i),tb31(i),tb31FG(i),tb50(i),tb89(i),tb165(i),satZenithAngle(i),obsLat(i), calcLandQualifierIndice(i) = ', &
     &             i,tb23(i),tb23FG(i),tb31(i),tb31FG(i),tb50(i),tb89(i),tb165(i),satZenithAngle(i),obsLat(i), calcLandQualifierIndice(i)
        write(*,*) ' ier(i),ice(i),cloudLiquidWaterPathObs(i),cloudLiquidWaterPathFG(i),si_ecmwf(i),si_bg(i),calcTerrainTypeIndice(i),waterobs(i) =',ier(i),ice(i),&
     &             cloudLiquidWaterPathObs(i),cloudLiquidWaterPathFG(i),si_ecmwf(i),si_bg(i),calcTerrainTypeIndice(i),waterobs(i)
      end if

    end do   ! i loop over numObsToProcess points

  end subroutine mwbg_nrlFilterAtms

  !--------------------------------------------------------------------------
  ! mwbg_nrlFilterMwhs2
  !--------------------------------------------------------------------------
  subroutine mwbg_nrlFilterMwhs2(numObsToProcess, actualNumChannel, obsTb, obsTbBiasCorr, satZenithAngle, obsLat, &
                                 calcLandQualifierIndice, calcTerrainTypeIndice, waterobs, grossrej, &
                                 cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                                 si_ecmwf, si_bg, iNumSeaIce, iRej, SeaIce)
    !OBJET          Compute the following parameters using 2 MWHS2 channels:
    !                  - sea ice,
    !                  - cloud liquid water from observation (cloudLiquidWaterPathObs),
    !                  - cloud liquid water from first guess (cloudLiquidWaterPathFG),
    !                  - 2 scattering indices (si) (ECMWF, Bennartz-Grody)
    !               The two channels used are: 89Ghz, and 165Ghz.
    !
    !NOTES*
    !                - open water points are converted to sea-ice points if sea ice concentration >= 0.55
    !                   and calcTerrainTypeIndice (terrainTypeIndice or terrain type) is changed accordingly
    !                - cloudLiquidWaterPathObs are missing when out-of-range parameters/Tb detected or grossrej = .true.
    !                - cloudLiquidWaterPathObs and si_ecmwf only computed over open water away from coasts and sea-ice
    !                - si_bg is computed for all points
    !                - cloudLiquidWaterPathObs and si = -99.0 where value cannot be computed.
    !
    !REFERENCES     Ben Ruston, NRL Monterey
    !                  JCSDA Seminar 12/12/12: Impact of NPP Satellite Assimilation in the U.S. Navy Global Modeling System
    !
    !
    !ARGUMENTS      - ier         - output - error return code for each location:
    !                                        0, ok,
    !                                        1, input parameter out of range or grossrej=.true.
    !               - numObsToProcess - input  -  number of points to process (= numObsToProcess)
    !               - tb23            - input  -  23Ghz brightness temperature (K) -- ch. 1 [missing for MWHS-2]
    !               - tb23FG          - input  -  23Ghz brightness temperature (K) from first guess -- ch. 1
    !               - tb31            - input  -  31Ghz brightness temperature (K) -- ch. 2 [missing for MWHS-2]
    !               - tb31FG          - input  -  31Ghz brightness temperature (K) from first guess -- ch. 2
    !               - tb50            - input  -  50Ghz brightness temperature (K) -- ch. 3 [missing for MWHS-2]
    !               - tb89            - input  -  89Ghz brightness temperature (K) -- ch. 16
    !               - tb165           - input  -  165Ghz brightness temperature (K) -- ch. 17
    !               - satZenithAngle  - input  -  satellite zenith angle (deg.)
    !               - obsLat          - input  -  latitude (deg.)
    !               - calcLandQualifierIndice   - input  -  land/sea indicator (0=land, 1=ocean)
    !               - calcTerrainTypeIndice     - in/out -  terrain type (0=ice, -1 otherwise)
    !               - waterobs    - in/out -  .true. if open water point (away from coasts and sea-ice)
    !               - grossrej    - input  -  .true. if any channel had a gross error from mwbg_grossValueCheck
    !               - cloudLiquidWaterPathObs - output -  cloud liquid water from observation (kg/m**2) from tb23 & tb31
    !               - cloudLiquidWaterPathFG  - output -  cloud liquid water from first guess (kg/m**2) from tb23FG & tb31FG
    !               - si_ecmwf    - output -  ECMWF scattering index from tb89 & tb165
    !               - si_bg       - output -  Bennartz-Grody scattering index from tb89 & tb165
    !               - iNumSeaIce  - in/out -  running counter for number of open water points
    !                                         with sea-ice detected (from algorithm)
    !               - iRej        - in/out -  running counter for number of locations with bad
    !                                         satZenithAngle, obsLat, calcLandQualifierIndice, or with grossrej=true
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

    integer, intent(in)                   ::  numObsToProcess
    integer, intent(in)                   ::  actualNumChannel
    integer, intent(out)                  ::  iNumSeaIce
    integer, intent(in)                   ::  calcLandQualifierIndice(:)
    integer, intent(inout)                ::  calcTerrainTypeIndice(:)
    integer, intent(out)                  ::  iRej


    logical, intent(in)                   ::  grossrej(:)
    logical, intent(inout)                ::  waterobs(:)

    real, intent(in)                      ::  obsTb(:)
    real, intent(in)                      ::  obsTbBiasCorr(:)
    real, intent(in)                      ::  satZenithAngle(:)
    real, intent(in)                      ::  obsLat(:)
    real, allocatable, intent(out)        ::  cloudLiquidWaterPathObs(:)
    real, allocatable, intent(out)        ::  cloudLiquidWaterPathFG(:)
    real, allocatable, intent(out)        ::  si_ecmwf(:)
    real, allocatable, intent(out)        ::  si_bg(:)
    real, allocatable, intent(out)        ::  SeaIce(:)

    ! Locals
    integer                               :: ier(numObsToProcess)
    real                                  ::  ice(numObsToProcess)
    real                                  :: tb23(numObsToProcess)
    real                                  :: tb23FG(numObsToProcess)
    real                                  :: tb31(numObsToProcess)
    real                                  :: tb31FG(numObsToProcess)
    real                                  :: tb50(numObsToProcess)
    real                                  :: tb89(numObsToProcess)
    real                                  :: tb165(numObsToProcess)
    real                                  :: bcor23(numObsToProcess)
    real                                  :: bcor31(numObsToProcess)
    real                                  :: bcor50(numObsToProcess)
    real                                  :: bcor89(numObsToProcess)
    real                                  :: bcor165(numObsToProcess)
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
    call utl_reAllocate(cloudLiquidWaterPathObs,numObsToProcess)
    call utl_reAllocate(cloudLiquidWaterPathFG,numObsToProcess)
    call utl_reAllocate(si_ecmwf,numObsToProcess)
    call utl_reAllocate(si_bg,numObsToProcess)
    call utl_reAllocate(SeaIce,numObsToProcess)

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
    do ii = 1, numObsToProcess
      indx2 = ii*actualNumChannel
      tb23(ii)      = mwbg_realMissing
      tb23FG(ii)    = mwbg_realMissing
      bcor23(ii)    = mwbg_realMissing
      tb31(ii)      = mwbg_realMissing
      tb31FG(ii)    = mwbg_realMissing
      bcor31(ii)    = mwbg_realMissing
      tb50(ii)      = mwbg_realMissing
      bcor50(ii)    = mwbg_realMissing
      tb89(ii)      = obsTb(indx1)
      bcor89(ii)    = obsTbBiasCorr(indx1)
      tb165(ii)    = obsTb(indx1+9)
      bcor165(ii)    = obsTbBiasCorr(indx1+9)
      indx1 = indx2 + 1
    end do

    ier = 0

    ! 1) Initialise parameters:
    do i = 1, numObsToProcess
      ice(i)      = mwbg_realMissing
      cloudLiquidWaterPathObs(i) = mwbg_realMissing
      cloudLiquidWaterPathFG(i)    = mwbg_realMissing
      si_ecmwf(i) = mwbg_realMissing
      si_bg(i)    = mwbg_realMissing
      SeaIce(i)   = 0.0
    end do

    ! 2) Validate input parameters:
    do i = 1, numObsToProcess
      if ( satZenithAngle(i) <  0. .or. &
           satZenithAngle(i) > 70. .or. &
           obsLat(i)  < -90.  .or. &
           obsLat(i)  >  90.  .or. &
           calcLandQualifierIndice(i) < 0 .or. &
           calcLandQualifierIndice(i) > 1 ) then
         ier(i) = 1
      end if

      ! Skip computations for points where all data are rejected  (bad Tb ANY channel)
      if ( grossrej(i) ) then
        ier(i) = 1
      end if
    end do

    ! 3) Compute parameters:
    do i = 1, numObsToProcess

      if ( ier(i) == 0 ) then

        abslat = abs(obsLat(i))
        cosz   = cosd(satZenithAngle(i))

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
        if ( calcLandQualifierIndice(i) == 1 .and. t23 /= mwbg_realMissing ) then  ! water point

          if ( abslat < 50. ) then
            ice(i) = 0.0
          else
            ice(i) = 2.85 + 0.020*t23 - 0.028*t50
          end if

          SeaIce(i) = ice(i)
          if ( ice(i) >= 0.55 .and. waterobs(i) ) then
            iNumSeaIce = iNumSeaIce + 1
            waterobs(i) = .false.
            calcTerrainTypeIndice(i) = 0
          end if

        end if

        ! Compute cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, and Scattering Indices (over open water only)
        if ( waterobs(i) ) then
          if ( t23 /= mwbg_realMissing ) then
            if ( t23 < 284. .and. t31 < 284. ) then
              aa = 8.24 - (2.622 - 1.846 * cosz) * cosz
              cloudLiquidWaterPathObs(i) = aa + 0.754 * alog(285.0 - t23) - 2.265 * alog(285.0 - t31)
              cloudLiquidWaterPathObs(i) = cloudLiquidWaterPathObs(i) * cosz
              if ( cloudLiquidWaterPathObs(i) < 0.0 ) cloudLiquidWaterPathObs(i) = 0.0

              cloudLiquidWaterPathFG(i) = aa + 0.754 * alog(285.0 - t23FG) - 2.265 * alog(285.0 - t31FG)
              cloudLiquidWaterPathFG(i) = cloudLiquidWaterPathFG(i) * cosz
              if ( cloudLiquidWaterPathFG(i) < 0.0 ) cloudLiquidWaterPathFG(i) = 0.0
            end if
          end if
          si_ecmwf(i) = deltb - (-46.94 + 0.248 * satZenithAngle(i))
          si_bg(i)    = deltb - (-39.201 + 0.1104 * satZenithAngle(i))
        else
          si_bg(i)    = deltb - (0.158 + 0.0163 * satZenithAngle(i))
        end if

      else  ! ier(i) == 1 case
         iRej = iRej + 1

      end if ! if ( ier(i) == 0 )

      if ( mwbg_debug .and. (i <= 100) ) then
        write(*,*) ' '
        write(*,*) ' i,tb23(i),tb23FG(i),tb31(i),tb31FG(i),tb50(i),tb89(i),tb165(i),satZenithAngle(i),obsLat(i), calcLandQualifierIndice(i) = ', &
     &             i,tb23(i),tb23FG(i),tb31(i),tb31FG(i),tb50(i),tb89(i),tb165(i),satZenithAngle(i),obsLat(i), calcLandQualifierIndice(i)
        write(*,*) ' ier(i),ice(i),cloudLiquidWaterPathObs(i),cloudLiquidWaterPathFG(i),si_ecmwf(i),si_bg(i),calcTerrainTypeIndice(i),waterobs(i) =',ier(i),ice(i),&
     &             cloudLiquidWaterPathObs(i),cloudLiquidWaterPathFG(i),si_ecmwf(i),si_bg(i),calcTerrainTypeIndice(i),waterobs(i)
      end if

    end do   ! i loop over numObsToProcess points

  end subroutine mwbg_nrlFilterMwhs2

  !--------------------------------------------------------------------------
  ! mwbg_flagDataUsingNrlCritAtms
  !--------------------------------------------------------------------------
  subroutine mwbg_flagDataUsingNrlCritAtms(numObsToProcess, actualNumChannel, obsTb, obsTbBiasCorr, cloudLiquidWaterPathObs, &
                                           scatec, scatbg, SeaIce, grossrej, waterobs, useUnbiasedObsForClw, &
                                           iwvreject, cloudobs, precipobs,  cldcnt, newInformationFlag, riwv, zdi)

    !:Purpose:                       Set the  Information flag (newInformationFlag) values (new BURP element 025174 in header)
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
    integer, intent(in)                        :: numObsToProcess
    integer, intent(in)                        :: actualNumChannel
    real, intent(in)                           :: obsTb(:)
    real, intent(in)                           :: obsTbBiasCorr(:)
    real, intent(in)                           :: cloudLiquidWaterPathObs (:)
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
    integer, allocatable, intent(out)          :: newInformationFlag(:)
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


    call utl_reAllocate(cloudobs, numObsToProcess)
    call utl_reAllocate(iwvreject, numObsToProcess)
    call utl_reAllocate(newInformationFlag, numObsToProcess)
    call utl_reAllocate(precipobs, numObsToProcess)
    call utl_reAllocate(riwv, numObsToProcess)
    call utl_reAllocate(ztb_amsub3, numObsToProcess)
    call utl_reAllocate(bcor_amsub3, numObsToProcess)
    call utl_reAllocate(ztb_amsub5, numObsToProcess)
    call utl_reAllocate(bcor_amsub5, numObsToProcess)

    ! To begin, assume that all obs are good.
    newInformationFlag(:) = 0
    cloudobs(:)  = .false.
    iwvreject(:) = .false.
    precipobs(:) = .false.

    ! Extract Tb for channels 16 (AMSU-B 1) and 17 (AMSU-B 2) for Bennartz SI
    ! Extract Tb for channels 22 (AMSU-B 3) and 18 (AMSU-B 5) for Dryness Index (DI)

    indx1 = 1
    do ii = 1, numObsToProcess
      indx2 = ii * actualNumChannel
      ztb_amsub3(ii) = obsTb(indx1+21)
      bcor_amsub3(ii) = obsTbBiasCorr(indx1+21)
      ztb_amsub5(ii) = obsTb(indx1+17)
      bcor_amsub5(ii) = obsTbBiasCorr(indx1+17)
      indx1 = indx2 + 1
    end do


    ! Flag data using NRL criteria

    ! Compute Mean 183 Ghz [ch. 18-22] Tb (riwv)
    riwv = -99.0
    indx1 = 1
    do ii = 1, numObsToProcess
      indx2 = ii * actualNumChannel
      if (.not.grossrej(ii)) then
        do indx = 1, 5
          if (obsTbBiasCorr(indx1+indx+9) == mwbg_realMissing .or. useUnbiasedObsForClw) then
            ztb183(indx) = obsTb(indx1+indx+16)
          else
            ztb183(indx) = obsTb(indx1+indx+16) - obsTbBiasCorr(indx1+indx+16)
          end if
        end do
        riwv(ii)  = sum(ztb183)/5.0
        if ( riwv(ii) < mean_Tb_183Ghz_min ) iwvreject(ii) = .true.
      else
        iwvreject(ii) = .true.
      end if
      indx1 = indx2 + 1
    end do

    !  Set bits in newInformationFlag flag to identify where various data selection criteria are met
    !     precipobs = .true  where ECMWF or BG scattering index > min_threshold (LT)
    !     cloudobs  = .true. where CLW > min_threshold (LT) or if precipobs = .true

    where ( grossrej ) newInformationFlag = IBSET(newInformationFlag,11)
    where ( scatec > scatec_atms_nrl_LTrej .or. scatbg > scatbg_atms_nrl_LTrej ) precipobs = .true.
    n_cld = count(cloudLiquidWaterPathObs > clw_atms_nrl_LTrej)
    cldcnt  = cldcnt  + n_cld
    where ( (cloudLiquidWaterPathObs > clw_atms_nrl_LTrej) .or. precipobs ) cloudobs = .true.
    where ( waterobs )  newInformationFlag = IBSET(newInformationFlag,0)
    where ( iwvreject ) newInformationFlag = IBSET(newInformationFlag,5)
    where ( precipobs ) newInformationFlag = IBSET(newInformationFlag,4)
    where ( cloudLiquidWaterPathObs > clw_atms_nrl_LTrej) newInformationFlag = IBSET(newInformationFlag,3)
    where ( cloudLiquidWaterPathObs > clw_atms_nrl_UTrej) newInformationFlag = IBSET(newInformationFlag,6)
    where ( scatec > scatec_atms_nrl_UTrej .or. scatbg > scatbg_atms_nrl_UTrej ) newInformationFlag = IBSET(newInformationFlag,8)
    where ( SeaIce >= 0.55 ) newInformationFlag = IBSET(newInformationFlag,10)

    where ( waterobs .and. (cloudLiquidWaterPathObs == -99.) ) newInformationFlag = IBSET(newInformationFlag,2)
    where ( riwv == -99.)                   newInformationFlag = IBSET(newInformationFlag,1)

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
  subroutine mwbg_flagDataUsingNrlCritMwhs2(numObsToProcess, actualNumChannel, obsTb, obsTbBiasCorr, cloudLiquidWaterPathObs, &
                                            scatec, SeaIce, grossrej, waterobs, useUnbiasedObsForClw, &
                                            iwvreject, cloudobs, precipobs,  cldcnt, newInformationFlag, riwv, zdi)

    !:Purpose:                       Set the  Information flag (newInformationFlag) values (new BURP element 025174 in header)
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
    integer, intent(in)                        :: numObsToProcess
    integer, intent(in)                        :: actualNumChannel
    real, intent(in)                           :: obsTb(:)
    real, intent(in)                           :: obsTbBiasCorr(:)
    real, intent(in)                           :: cloudLiquidWaterPathObs (:)
    real, intent(in)                           :: scatec(:)
    real, intent(in)                           :: SeaIce (:)

    logical, intent(in)                        :: useUnbiasedObsForClw
    logical, intent(in)                        :: grossrej(:)
    logical, intent(in)                        :: waterobs(:)
    integer, intent(inout)                     :: cldcnt
    logical, allocatable, intent(out)          :: cloudobs(:)
    logical, allocatable, intent(out)          :: iwvreject(:)
    logical, allocatable, intent(out)          :: precipobs(:)
    integer, allocatable, intent(out)          :: newInformationFlag(:)
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


    call utl_reAllocate(cloudobs, numObsToProcess)
    call utl_reAllocate(iwvreject, numObsToProcess)
    call utl_reAllocate(newInformationFlag, numObsToProcess)
    call utl_reAllocate(precipobs, numObsToProcess)
    call utl_reAllocate(riwv, numObsToProcess)
    call utl_reAllocate(ztb_amsub3, numObsToProcess)
    call utl_reAllocate(bcor_amsub3, numObsToProcess)
    call utl_reAllocate(ztb_amsub5, numObsToProcess)
    call utl_reAllocate(bcor_amsub5, numObsToProcess)

    ! To begin, assume that all obs are good.
    newInformationFlag(:) = 0
    cloudobs(:)  = .false.
    iwvreject(:) = .false.
    precipobs(:) = .false.

    ! Extract Tb for channels 1 (AMSU-B 1) and 10 (AMSU-B 2) for Bennartz SI
    ! Extract Tb for channels 11 (AMSU-B 3) and 15 (AMSU-B 5) for Dryness Index (DI)

    indx1 = 1
    do ii = 1, numObsToProcess
      indx2 = ii * actualNumChannel
      ztb_amsub3(ii) = obsTb(indx1+10)
      bcor_amsub3(ii) = obsTbBiasCorr(indx1+10)
      ztb_amsub5(ii) = obsTb(indx1+14)
      bcor_amsub5(ii) = obsTbBiasCorr(indx1+14)
      indx1 = indx2 + 1
    end do


    ! Flag data using NRL criteria

    ! Compute Mean 183 Ghz [ch. 11-15] Tb (riwv)
    riwv = -99.0
    indx1 = 1
    do ii = 1, numObsToProcess
      indx2 = ii * actualNumChannel
      if (.not.grossrej(ii)) then
        do indx = 1, 5
          if (obsTbBiasCorr(indx1+indx+9) == mwbg_realMissing .or. useUnbiasedObsForClw) then
            ztb183(indx) = obsTb(indx1+indx+9)
          else
            ztb183(indx) = obsTb(indx1+indx+9) - obsTbBiasCorr(indx1+indx+9)
          end if
        end do
        riwv(ii)  = sum(ztb183)/5.0
        if ( riwv(ii) < mean_Tb_183Ghz_min ) iwvreject(ii) = .true.
      else
        iwvreject(ii) = .true.
      end if
      indx1 = indx2 + 1
    end do

    !  Set bits in newInformationFlag flag to identify where various data selection criteria are met
    !     precipobs = .true  where ECMWF or BG scattering index > min_threshold (LT)
    !     cloudobs  = .true. where CLW > min_threshold (LT) or if precipobs = .true

    where ( grossrej ) newInformationFlag = IBSET(newInformationFlag,11)
    where ( scatec > scatec_mwhs2_nrl_LTrej ) precipobs = .true.
    n_cld = count(cloudLiquidWaterPathObs > clw_mwhs2_nrl_LTrej)
    cldcnt  = cldcnt  + n_cld
    where ( (cloudLiquidWaterPathObs > clw_mwhs2_nrl_LTrej) .or. precipobs ) cloudobs = .true.
    where ( waterobs )  newInformationFlag = IBSET(newInformationFlag,0)
    where ( iwvreject ) newInformationFlag = IBSET(newInformationFlag,5)
    where ( precipobs ) newInformationFlag = IBSET(newInformationFlag,4)
    where ( cloudLiquidWaterPathObs > clw_mwhs2_nrl_LTrej) newInformationFlag = IBSET(newInformationFlag,3)
    where ( cloudLiquidWaterPathObs > clw_mwhs2_nrl_UTrej) newInformationFlag = IBSET(newInformationFlag,6)
    where ( SeaIce >= 0.55 ) newInformationFlag = IBSET(newInformationFlag,10)

    where ( waterobs .and. (cloudLiquidWaterPathObs == -99.) ) newInformationFlag = IBSET(newInformationFlag,2)
    where ( riwv == -99.)                   newInformationFlag = IBSET(newInformationFlag,1)

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
  subroutine mwbg_reviewAllCritforFinalFlagsAtms(numObsToProcess, actualNumChannel, qcRejectLogic, grossrej, waterobs, &
                                                 precipobs, cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, scatec, scatbg, &
                                                 scatIndexOverWaterObs, scatIndexOverWaterFG, iwvreject, riwv, obsFlags, &
                                                 obsGlobalMarker, zdi, newInformationFlag, drycnt, landcnt, &
                                                 rejcnt, iwvcnt, pcpcnt, flgcnt, &
                                                 chanIgnoreInAllskyGenCoeff, obsChannels)

    !:Purpose:                   Review all the checks previously made to determine which obs are to be accepted
    !                            for assimilation and which are to be flagged for exclusion (lflagchn).
    !                            - grossrej()  = .true. if any channel had a gross error at the point
    !                            - cloudobs()  = .true. if CLW > clw_atms_nrl_LTrej (0.175) or precipobs
    !                            - precipobs() = .true. if precip. detected through NRL scattering indices
    !                            - waterobs()  = .true. if open water point
    !                            - iwvreject() = .true. if Mean 183 Ghz [ch. 18-22] Tb < 240K (too dry for ch.20-22 over land)
    ! Arguments
    integer, intent(in)                        :: numObsToProcess
    integer, intent(in)                        :: actualNumChannel
    logical, intent(in)                        :: qcRejectLogic(:,:)
    real, intent(inout)                        :: cloudLiquidWaterPathObs(:)
    real, intent(inout)                        :: cloudLiquidWaterPathFG(:)
    real, intent(in)                           :: scatec(:)
    real, intent(in)                           :: scatbg(:)
    real, allocatable, intent(out)             :: scatIndexOverWaterObs(:)
    real, allocatable, intent(out)             :: scatIndexOverWaterFG(:)    
    logical, intent(in)                        :: grossrej(:)
    logical, intent(in)                        :: waterobs(:)
    logical, intent(in)                        :: iwvreject(:)
    logical, intent(in)                        :: precipobs(:)
    integer, intent(inout)                     :: newInformationFlag(:)
    real, intent(in)                           :: zdi(:)
    real, intent(in)                           :: riwv(:)
    integer, intent(inout)                     :: obsFlags(:)
    integer, intent(inout)                     :: obsGlobalMarker(:)
    integer, intent(inout)                     :: drycnt
    integer, intent(inout)                     :: landcnt
    integer, intent(inout)                     :: rejcnt
    integer, intent(inout)                     :: iwvcnt
    integer, intent(inout)                     :: pcpcnt
    integer, intent(inout)                     :: flgcnt
    integer, intent(in)                        :: chanIgnoreInAllskyGenCoeff(:)
    integer, intent(in)                        :: obsChannels(:)

    ! Locals
    real                                       :: clwObsFGaveraged
    logical, allocatable                       :: lflagchn(:,:)
    integer                                    :: kk, j, ipos, INDXCAN


    ! Allocation
    call utl_reAllocate(lflagchn,numObsToProcess, actualNumChannel)
    call utl_reAllocate(scatIndexOverWaterObs,numObsToProcess)
    call utl_reAllocate(scatIndexOverWaterFG, numObsToProcess)

    lflagchn(:,:) = qcRejectLogic(:,:)  ! initialize with flags set in mwbg_firstQcCheckAtms
    do kk = 1, numObsToProcess
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
            newInformationFlag(kk) = IBSET(newInformationFlag(kk),7)
          end if
          if ( zdi(kk) > -5.0 ) then
            lflagchn(kk,21) = .true.
            newInformationFlag(kk) = IBSET(newInformationFlag(kk),9)
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

          if ( cloudLiquidWaterPathObs(kk)   >  clw_atms_nrl_LTrej )  then
            if ( tvs_mwAllskyAssim ) then
              lflagchn(kk,1:4) = .true.
              clwObsFGaveraged = 0.5 * (cloudLiquidWaterPathObs(kk) + cloudLiquidWaterPathFG(kk))
              if ( clwObsFGaveraged > mwbg_clwQcThreshold ) lflagchn(kk,5:6) = .true.
            else
              lflagchn(kk,1:mwbg_atmsNumSfcSensitiveChannel) = .true.
            end if
            lflagchn(kk,16:20) = .true.
          end if
          if ( cloudLiquidWaterPathObs(kk) > clw_atms_nrl_UTrej )  then
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
          if ( cloudLiquidWaterPathObs(kk) == -99. ) then
            newInformationFlag(kk) = IBSET(newInformationFlag(kk),2)
            lflagchn(kk,1:9)   = .true.
            lflagchn(kk,16:22) = .true.
          end if
          if ( riwv(kk) == -99. ) then     ! riwv = mean_Tb_183Ghz
            newInformationFlag(kk) = IBSET(newInformationFlag(kk),1)
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

    ! RESET scatIndexOverWaterObs array to ECMWF scattering index for output to BURP file
    scatIndexOverWaterObs(:) = scatec(:)
    ! Set missing cloudLiquidWaterPathObs and scatIndexOverWaterObs to BURP missing value (mwbg_realMissing)
    where (cloudLiquidWaterPathObs == -99. ) cloudLiquidWaterPathObs = mwbg_realMissing
    where (cloudLiquidWaterPathObs == -99. ) cloudLiquidWaterPathFG = mwbg_realMissing
    where (scatIndexOverWaterObs == -99. ) scatIndexOverWaterObs = mwbg_realMissing
    scatIndexOverWaterFG(:) = mwbg_realMissing

    ! Modify data flag values (set bit 7) for rejected data
    ! In all-sky mode, turn on bit=23 for channels in chanIgnoreInAllskyGenCoeff(:)
    ! as cloud-affected radiances over sea when there is mismatch between 
    ! cloudLiquidWaterPathObs and cloudLiquidWaterPathFG (to be used in gen_bias_corr)
    ipos=0
    do kk =1, numObsToProcess
      clwObsFGaveraged = 0.5 * (cloudLiquidWaterPathObs(kk) + cloudLiquidWaterPathFG(kk))

      do j = 1, actualNumChannel
        ipos = ipos + 1
        if (lflagchn(kk,j)) then
          obsFlags(ipos) = IBSET(obsFlags(ipos),7)
        end if

        INDXCAN = ISRCHEQI(chanIgnoreInAllskyGenCoeff,obsChannels(ipos))
        if (tvs_mwAllskyAssim .and. waterobs(kk) .and. INDXCAN /= 0 .and. &
            (clwObsFGaveraged > mwbg_cloudyClwThresholdBcorr .or. &
             cloudLiquidWaterPathObs(kk) == mwbg_realMissing .or. cloudLiquidWaterPathFG(kk) == mwbg_realMissing)) then
          obsFlags(ipos) = IBSET(obsFlags(ipos),23)
        end if
      end do
    end do


    ! Set bit 6 in 24-bit global flags if any data rejected
    do kk =1, numObsToProcess
      if ( ANY(lflagchn(kk,:)) ) obsGlobalMarker(kk) = IBSET(obsGlobalMarker(kk),6)
    end do

  end subroutine mwbg_reviewAllCritforFinalFlagsAtms

  !--------------------------------------------------------------------------
  ! mwbg_reviewAllCritforFinalFlagsMwhs2
  !--------------------------------------------------------------------------
  subroutine mwbg_reviewAllCritforFinalFlagsMwhs2(numObsToProcess, actualNumChannel, qcRejectLogic, grossrej, calcTerrainTypeIndice, waterobs, &
                                                  precipobs, cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, scatec, scatbg, &
                                                  scatIndexOverWaterObs, scatIndexOverWaterFG, iwvreject, riwv, obsFlags, &
                                                  obsGlobalMarker, zdi, newInformationFlag, allcnt, drycnt, landcnt, &
                                                  rejcnt, iwvcnt, pcpcnt, flgcnt, &
                                                  chanIgnoreInAllskyGenCoeff, obsChannels)

    !:Purpose:                   Review all the checks previously made to determine which obs are to be accepted
    !                            for assimilation and which are to be flagged for exclusion (lflagchn).
    !                            - grossrej()  = .true. if any channel had a gross error at the point
    !                            - cloudobs()  = .true. if CLW > clw_mwhs2_nrl_LTrej (0.175) or precipobs
    !                            - precipobs() = .true. if precip. detected through NRL scattering indices
    !                            - waterobs()  = .true. if open water point
    !                            - iwvreject() = .true. if Mean 183 Ghz [ch. 18-22] Tb < 240K (too dry for ch.20-22 over land)
    ! Arguments
    integer, intent(in)                        :: numObsToProcess
    integer, intent(in)                        :: actualNumChannel
    logical, intent(in)                        :: qcRejectLogic(:,:)
    real, intent(inout)                        :: cloudLiquidWaterPathObs(:)
    real, intent(inout)                        :: cloudLiquidWaterPathFG(:)
    real, intent(in)                           :: scatec(:)
    real, intent(in)                           :: scatbg(:)
    real, allocatable, intent(out)             :: scatIndexOverWaterObs(:)
    real, allocatable, intent(out)             :: scatIndexOverWaterFG(:)
    logical, intent(in)                        :: grossrej(:)
    logical, intent(in)                        :: waterobs(:)
    logical, intent(in)                        :: iwvreject(:)
    logical, intent(in)                        :: precipobs(:)
    integer, intent(inout)                     :: newInformationFlag(:)
    real, intent(in)                           :: zdi(:)
    real, intent(in)                           :: riwv(:)
    integer, intent(inout)                     :: obsFlags(:)
    integer, intent(inout)                     :: obsGlobalMarker(:)
    integer, intent(inout)                     :: allcnt
    integer, intent(inout)                     :: drycnt
    integer, intent(inout)                     :: landcnt
    integer, intent(inout)                     :: rejcnt
    integer, intent(inout)                     :: iwvcnt
    integer, intent(inout)                     :: pcpcnt
    integer, intent(inout)                     :: flgcnt
    integer, intent(inout)                     :: calcTerrainTypeIndice(:)
    integer, intent(in)                        :: chanIgnoreInAllskyGenCoeff(:)
    integer, intent(in)                        :: obsChannels(:)

    ! Locals
    real                                       :: clwObsFGaveraged
    real                                       :: scatbg_rej
    logical, allocatable                       :: lflagchn(:,:)
    integer                                    :: kk, j, ipos, INDXCAN


    ! Allocation
    call utl_reAllocate(lflagchn,numObsToProcess, actualNumChannel)
    call utl_reAllocate(scatIndexOverWaterObs,numObsToProcess)
    call utl_reAllocate(scatIndexOverWaterFG, numObsToProcess)

    lflagchn(:,:) = qcRejectLogic(:,:)  ! initialize with flags set in mwbg_firstQcCheckMwhs2
    do kk = 1, numObsToProcess
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
            newInformationFlag(kk) = IBSET(newInformationFlag(kk),7)
          end if
          if ( zdi(kk) > -5.0 ) then
            lflagchn(kk,12) = .true.
            newInformationFlag(kk) = IBSET(newInformationFlag(kk),9)
            drycnt = drycnt + 1
          end if
          if ( zdi(kk) > -8.0 ) then
            lflagchn(kk,13) = .true.
          end if

          ! Bennartz -Grody SI check thresholds (same as for QC of AMSU-B/MHS)
          if ( calcTerrainTypeIndice(kk) == 0 ) then ! sea-ice
            scatbg_rej = scatbg_mwhs2_cmc_ICErej
          else                     ! land
            scatbg_rej = scatbg_mwhs2_cmc_LANDrej
          end if
          if ( scatbg(kk) > scatbg_rej ) then
            lflagchn(kk,:) = .true.
            newInformationFlag(kk) = IBSET(newInformationFlag(kk),8)
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

          if ( cloudLiquidWaterPathObs(kk) > clw_mwhs2_nrl_LTrej )  then
            if ( tvs_mwAllskyAssim ) then ! NEVER TRUE SINCE NOT IMPLEMENTED YET
              clwObsFGaveraged = 0.5 * (cloudLiquidWaterPathObs(kk) + cloudLiquidWaterPathFG(kk))
              if ( clwObsFGaveraged > mwbg_clwQcThreshold ) lflagchn(kk,1) = .true.
            else
              lflagchn(kk,1) = .true.
            end if
            lflagchn(kk,(/ 10,13,14,15 /)) = .true.
          end if
          if ( cloudLiquidWaterPathObs(kk) > clw_mwhs2_nrl_UTrej )  then
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
            newInformationFlag(kk) = IBSET(newInformationFlag(kk),1)
            lflagchn(kk,1) = .true.
            lflagchn(kk,10:15) = .true.
          end if
          ! Bennartz-Grody SI check thresholds (same as for QC of AMSU-B/MHS)
          if ( scatbg(kk) > scatbg_mwhs2_cmc_SEA ) then
            lflagchn(kk,:) = .true.
            newInformationFlag(kk) = IBSET(newInformationFlag(kk),8)
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

    ! RESET scatIndexOverWaterObs array to Bennartz-Grody scattering index for output to BURP file
    scatIndexOverWaterObs(:) = scatbg(:)
    ! Set missing cloudLiquidWaterPathObs and scatIndexOverWaterObs to BURP missing value (mwbg_realMissing)
    where (cloudLiquidWaterPathObs == -99. ) cloudLiquidWaterPathObs = mwbg_realMissing
    where (cloudLiquidWaterPathObs == -99. ) cloudLiquidWaterPathFG = mwbg_realMissing
    where (scatIndexOverWaterObs == -99. ) scatIndexOverWaterObs = mwbg_realMissing
    scatIndexOverWaterFG(:) = mwbg_realMissing

    ! Modify data flag values (set bit 7) for rejected data
    ! In all-sky mode, turn on bit=23 for channels in chanIgnoreInAllskyGenCoeff(:)
    ! as cloud-affected radiances over sea when there is mismatch between 
    ! cloudLiquidWaterPathObs and cloudLiquidWaterPathFG (to be used in gen_bias_corr)
    ipos=0
    do kk =1, numObsToProcess
      clwObsFGaveraged = 0.5 * (cloudLiquidWaterPathObs(kk) + cloudLiquidWaterPathFG(kk))

      do j = 1, actualNumChannel
        ipos = ipos + 1
        if (lflagchn(kk,j)) then
          obsFlags(ipos) = IBSET(obsFlags(ipos),7)
        end if

        INDXCAN = ISRCHEQI(chanIgnoreInAllskyGenCoeff,obsChannels(ipos))
        if (tvs_mwAllskyAssim .and. waterobs(kk) .and. INDXCAN /= 0 .and. &
            (clwObsFGaveraged > mwbg_cloudyClwThresholdBcorr .or. &
             cloudLiquidWaterPathObs(kk) == mwbg_realMissing .or. cloudLiquidWaterPathFG(kk) == mwbg_realMissing)) then
          obsFlags(ipos) = IBSET(obsFlags(ipos),23)
        end if
      end do
    end do


    ! Set bit 6 in 24-bit global flags if any data rejected
    do kk =1, numObsToProcess
      if ( ANY(lflagchn(kk,:)) ) obsGlobalMarker(kk) = IBSET(obsGlobalMarker(kk),6)
    end do

  end subroutine mwbg_reviewAllCritforFinalFlagsMwhs2

  function calcStateDepObsErr_r4(cldPredThresh1, cldPredThresh2, &
                                 errThresh1, errThresh2, cldPredUsed) result(sigmaObsErrUsed)
    !
    ! :Purpose: Calculate single-precision state-dependent observation error.
    !                                 
    implicit none

    ! Arguments:
    real, intent(in) :: cldPredThresh1  ! first cloud predictor threshold
    real, intent(in) :: cldPredThresh2  ! second cloud predictor threshold
    real, intent(in) :: errThresh1      ! sigmaO corresponding to first cloud predictor threshold
    real, intent(in) :: errThresh2      ! sigmaO corresponding to second cloud predictor threshold
    real, intent(in) :: cldPredUsed     ! cloud predictor for the obs
    real :: sigmaObsErrUsed             ! estimated sigmaO for the obs

    if (cldPredUsed <= cldPredThresh1) then
      sigmaObsErrUsed = errThresh1
    else if (cldPredUsed >  cldPredThresh1 .and. & 
             cldPredUsed <= cldPredThresh2) then
      sigmaObsErrUsed = errThresh1 + &
                        (errThresh2 - errThresh1) / &
                        (cldPredThresh2 - cldPredThresh1) * &
                        (cldPredUsed - cldPredThresh1) 
    else
      sigmaObsErrUsed = errThresh2
    end if

  end function calcStateDepObsErr_r4

  !--------------------------------------------------------------------------
  ! mwbg_updateObsSpaceAfterQc
  !--------------------------------------------------------------------------
  subroutine mwbg_updateObsSpaceAfterQc(obsSpaceData, sensorIndex, headerIndex, obsTb, obsFlags, &
                                        cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                                        scatIndexOverWaterObs, scatIndexOverWaterFG, &
                                        obsGlobalMarker, newInformationFlag)

    !:Purpose:      Update obspacedata variables (obstTB and obs flags) after QC
    implicit None

    !Arguments
    type(struct_obs),     intent(inout)     :: obsSpaceData               ! obspaceData Object
    integer,              intent(in)        :: sensorIndex                ! tvs_sensorIndex
    integer,              intent(in)        :: headerIndex                ! current header index
    integer,              intent(in)        :: obsFlags(:)                ! data flags
    real,                 intent(in)        :: obsTb(:)                   ! obs Tb
    real,                 intent(in)        :: cloudLiquidWaterPathObs(:) ! obs CLW
    real,                 intent(in)        :: cloudLiquidWaterPathFG(:)  ! trial CLW
    real,                 intent(in)        :: scatIndexOverWaterObs(:)   ! atmospheric scatering index from observation
    real,                 intent(in)        :: scatIndexOverWaterFG(:)    ! atmospheric scatering index from background
    integer,              intent(in)        :: newInformationFlag(:)      ! information flag used with satplot
    integer,              intent(in)        :: obsGlobalMarker(:)         ! information flag used with satplot
    ! Locals
    integer                                 :: bodyIndex
    integer                                 :: obsNumCurrentLoc
    integer                                 :: bodyIndexbeg
    integer                                 :: currentChannelNumber
    integer                                 :: codtyp

    codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
    call obs_headSet_r(obsSpaceData, OBS_CLWO, headerIndex, cloudLiquidWaterPathObs(1))

    if (tvs_isInstrumAllskyTtAssim(tvs_getInstrumentId(codtyp_get_name(codtyp)))) then
      call obs_headSet_r(obsSpaceData, OBS_CLWB, headerIndex, cloudLiquidWaterPathFG(1))
    end if

    if (scatIndexOverWaterObs(1) /= mwbg_realMissing) then
      call obs_headSet_r(obsSpaceData, OBS_SIO, headerIndex, scatIndexOverWaterObs(1))
    else
      call obs_headSet_r(obsSpaceData, OBS_SIO, headerIndex, MPC_missingValue_R4)
    end if

    if (tvs_isInstrumAllskyHuAssim(tvs_getInstrumentId(codtyp_get_name(codtyp)))) then
      if (scatIndexOverWaterFG(1) /= mwbg_realMissing) then
        call obs_headSet_r(obsSpaceData, OBS_SIB, headerIndex, scatIndexOverWaterFG(1))
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
                                      terrainTypeIndice, obsLat, obsLon, &
                                      satScanPosition, obsQcFlag1, satOrbit, & 
                                      obsGlobalMarker, stnId, obsTb, btClear, &
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
    real   , allocatable, intent(out)    :: obsLat(:)              ! latitude values (btyp=5120,ele=5002)
    real   , allocatable, intent(out)    :: obsLon(:)              ! longitude values (btyp=5120,ele=6002)
    integer, allocatable, intent(out)    :: satScanPosition(:)     ! scan position (fov)    (btyp=3072,ele=5043)
    integer, allocatable, intent(out)    :: obsQcFlag1(:,:)        ! flag values for btyp=3072 block ele 033078, 033079, 033080
    integer, allocatable, intent(out)    :: satOrbit(:)            ! orbit number
    integer, allocatable, intent(out)    :: obsGlobalMarker(:)     ! global Marqueur Data
    character(*),intent(out)             :: stnId                  ! Platform Name
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
    call utl_reAllocate(obsLat, numObsToProcess)
    call utl_reAllocate(obsLon, numObsToProcess)
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

        
    stnId                              = obs_elem_c    (obsSpaceData, 'STID' , headerIndex) 
    satIdentifier(headerCompt)         = obs_headElem_i(obsSpaceData, OBS_SAT, headerIndex) 
    satZenithAngle(headerCompt)        = obs_headElem_r(obsSpaceData, OBS_SZA, headerIndex) 
    landQualifierIndice(headerCompt)   = obs_headElem_i(obsSpaceData, OBS_STYP, headerIndex) 
    terrainTypeIndice(headerCompt)     = obs_headElem_i(obsSpaceData, OBS_TTYP, headerIndex) 
    ! If terrain type is missing, set it to -1 for the QC programs
    if (terrainTypeIndice(headerCompt) ==  99) terrainTypeIndice(headerCompt) = -1
    obsLat (headerCompt)               = obs_headElem_r(obsSpaceData, OBS_LAT, headerIndex) 
    obsLon(headerCompt)                = obs_headElem_r(obsSpaceData, OBS_LON, headerIndex) 
    ! Convert lat/lon to degrees
    obsLon(headerCompt) = obsLon(headerCompt)*MPC_DEGREES_PER_RADIAN_R8
    if (obsLon(headerCompt) > 180.) obsLon(headerCompt) = obsLon(headerCompt) - 360.
    obsLat(headerCompt)       = obsLat(headerCompt) *MPC_DEGREES_PER_RADIAN_R8
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
          tvs_isInstrumAllskyHuAssim(tvs_getInstrumentId(codtyp_get_name(codtyp)))) then
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
    character(len=9)              :: stnId                         ! station id in burp file
    real, allocatable             :: modelInterpTerrain(:)         ! topo in standard file interpolated to obs point
    real, allocatable             :: modelInterpSeaIce(:)          ! Glace de mer " "
    real, allocatable             :: modelInterpGroundIce(:)       ! Glace de continent " "
    real,    allocatable          :: obsLat(:)                     ! obs. point latitudes
    real,    allocatable          :: obsLon(:)                     ! obs. point longitude
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
    real,    allocatable          :: scatIndexOverWaterObs(:)      ! scattering index from observation.
    real,    allocatable          :: scatIndexOverWaterFG(:)       ! scattering index from background.
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
                                   terrainTypeIndice, obsLat, obsLon,      &
                                   satScanPosition, obsQcFlag1, satOrbit,             &
                                   obsGlobalMarker, stnId, obsTb, btClear,    &
                                   obsTbBiasCorr, ompTb, obsQcFlag2, obsChannels,     &
                                   obsFlags, sensorIndex, actualNumChannel, obsSpaceData)

      !###############################################################################
      ! STEP 3) Interpolation de le champ MX(topogrpahy), MG et GL aux pts TOVS.
      !###############################################################################
      call mwbg_readGeophysicFieldsAndInterpolate(instName, obsLat, obsLon, modelInterpTerrain,     &
                                                  modelInterpGroundIce, modelInterpSeaIce)
      !###############################################################################
      ! STEP 4) Controle de qualite des TOVS. Data QC flags (obsFlags) are modified here!
      !###############################################################################

      if (instName == 'AMSUA') then
        call mwbg_tovCheckAmsua(landQualifierIndice, obsChannels, obsTb, btClear, obsTbBiasCorr, &
                                ompTb, qcIndicator, actualNumChannel, numObsToProcess, sensorIndex, &
                                satScanPosition, modelInterpGroundIce, modelInterpTerrain,&
                                modelInterpSeaIce, terrainTypeIndice, satZenithAngle,     &
                                obsGlobalMarker, obsFlags, newInformationFlag, &
                                cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                                scatIndexOverWaterObs, scatIndexOverWaterFG, &
                                stnId, RESETQC, obsLat)
      else if (instName == 'AMSUB') then
        call mwbg_tovCheckAmsub(landQualifierIndice, obsChannels, obsTb, btClear, obsTbBiasCorr, &
                                ompTb, qcIndicator, actualNumChannel, numObsToProcess, sensorIndex, &
                                satScanPosition, modelInterpGroundIce, modelInterpTerrain,&
                                modelInterpSeaIce, terrainTypeIndice, satZenithAngle,     &
                                obsGlobalMarker, obsFlags, newInformationFlag,        & 
                                cloudLiquidWaterPathObs, cloudLiquidWaterPathFG,       &
                                scatIndexOverWaterObs, scatIndexOverWaterFG, stnId, RESETQC)
      else if (instName == 'ATMS') then
        call mwbg_tovCheckAtms(obsLat, obsLon, landQualifierIndice, terrainTypeIndice, satZenithAngle,   &
                               obsQcFlag2, obsQcFlag1, obsChannels, obsTb, obsTbBiasCorr, ompTb, qcIndicator, &
                               actualNumChannel, numObsToProcess, sensorIndex, newInformationFlag, &
                               satScanPosition, modelInterpTerrain, obsGlobalMarker, obsFlags, &
                               cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                               scatIndexOverWaterObs, scatIndexOverWaterFG, &
                               stnId, RESETQC)
      else if (instName == 'MWHS2') then
        call mwbg_tovCheckMwhs2(obsLat, obsLon, landQualifierIndice, terrainTypeIndice, satZenithAngle,  &
                                obsChannels, obsTb, obsTbBiasCorr, ompTb, qcIndicator,   &
                                actualNumChannel, numObsToProcess, sensorIndex,          &
                                newInformationFlag, satScanPosition,   &
                                modelInterpTerrain, obsGlobalMarker, obsFlags,            &
                                cloudLiquidWaterPathObs, cloudLiquidWaterPathFG, &
                                scatIndexOverWaterObs, scatIndexOverWaterFG, &
                                stnId, RESETQC, modLSQ, lastHeader)
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
                                      scatIndexOverWaterObs, scatIndexOverWaterFG, &
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

