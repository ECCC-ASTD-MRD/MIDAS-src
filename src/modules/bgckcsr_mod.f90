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

module bgckcsr_mod
  ! MODULE bgckcsr_mod(prefix='csrbg' category='1. High-level functionality')
  !
  ! :Purpose: To perform CSR data background Check
  !
  use midasMpi_mod
  use utilities_mod
  use obsSpaceData_mod
  use tovs_nl_mod
  use obsErrors_mod

  implicit none
  save
  private

  real,    parameter :: csrbg_realMissing=-99.
  
  ! Public functions/subroutines

  public :: csrbg_bgCheckCSR
  
  real, parameter :: csrbg_ompThreshold = 4.2
  integer, parameter :: maxNumsat  = 20     ! nb max de satellites
  integer, parameter :: maxNumchan = 15     ! nb max de canaux

  ! namelist variables
  character (len=9)  :: burpSatName(maxNumSat)
  integer            :: satCloudCoverLimit(maxNumSat,maxNumchan)

  namelist /namcsr/burpSatName,satCloudCoverLimit

contains

  !----------------------------------------------------------------------------------------
  ! csrbg_init
  !----------------------------------------------------------------------------------------
  subroutine csrbg_init()
    !
    !:Purpose: This subroutine reads the namelist section NAMCSR
    !          for the module.
    implicit none

    ! Locals:
    integer :: nulnam, ierr
    integer, external :: fnom, fclos

    nulnam = 0
    ierr = fnom(nulnam, './flnml','FTN+SEQ+R/O', 0)
    read(nulnam, nml=namcsr, iostat=ierr)
    if (ierr /= 0) call utl_abort('csrbg_init: Error reading namelist')
    if (mmpi_myid == 0) write(*, nml=namcsr)
    ierr = fclos(nulnam)

  end subroutine csrbg_init
  
  !----------------------------------------------------------------------------------------
  ! csrbg_bgCheckCSR
  !----------------------------------------------------------------------------------------

  subroutine csrbg_bgCheckCSR (obsSpaceData)

    !: Purpose: Effectuer le controle que qualite sur les donnees CSR.  
    !           Modifier les marqueurs de donnees selon le type de rejet.
    implicit none

    !argument:
    type(struct_obs),     intent(inout)  :: obsSpaceData           ! obspaceData Object

    ! Locals

    real   , allocatable  :: obsTb(:)            ! brightness temperature (btyp=9248/9264,ele=12163)
    real   , allocatable  :: ompTb(:)            ! OMP values
    real   , allocatable  :: satZenithAngle(:)   ! satellite zenith angle (btyp=3072,ele=7024)
    integer, allocatable  :: obsFlags(:)         ! data flags
    integer, allocatable  :: cloudAmount(:)      ! data flags
    integer, allocatable  :: obsChannels(:)      ! Tb Channels 
    integer, allocatable  :: obsDate(:)          ! date YYYYMMDD
    integer, allocatable  :: obsHour(:)          ! Hour HHMM
    integer               :: sensorIndex         ! find tvs_sensor index corresponding to current obs
    character(len=9)          :: burpFileSatId       ! Platform Name

    ! Data derived from elements read in obspace data

    logical, allocatable  :: maxAngleReached(:)     ! satellite angle exceed max angle at obs
    logical, allocatable  :: topographicData(:)  ! data flagged as topo data
    logical, allocatable  :: nonCorrectedData(:) ! data non corrected by bias corr
    logical, allocatable  :: isTbPresent(:)      ! non missing data
    logical, allocatable  :: isClearSky(:)       ! clear sky obs
    logical, allocatable  :: strayLight(:)       !
    logical, allocatable  :: goesMidi(:)         ! goes noon
    logical, allocatable  :: isToAssim(:)            ! is channel assimilable
    logical, allocatable  :: ompOutOfRange(:)        
    integer               :: categorieRejet(7)        
    integer               :: headerIndex 
    integer               :: codtyp

    logical               :: csrDataPresent


    

    categorieRejet(:) = 0
    csrDataPresent = .false.
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER0: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER0
      codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if (tvs_isIdBurpInst(codtyp,'radianceclear'  )) csrDataPresent = .true.
    end do HEADER0

    if ( .not. csrDataPresent ) then
      write(*,*) 'WARNING: WILL NOT RUN csrbg_bgCheckCSR since no CSR Data'
      return
    end if

    call utl_tmg_start(116,'--BgckCSR')
    write(*,*) ' CSRBG QC PROGRAM STARTS ....'
    ! Read Namelist
    call csrbg_init()

    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER
      codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)
      if ( .not. tvs_isIdBurpInst(codtyp,'radianceclear') ) then 
        write(*,*) 'WARNING: Observation with codtyp = ', codtyp, ' is not a CSR data'
        cycle HEADER
      end if

      !###############################################################################
      ! STEP 1) read obs from obsSpacedata to start QC                               !
      !###############################################################################

      call csrbg_readObsFromObsSpace(obsSpaceData, headerIndex, obsTb, ompTb, satZenithAngle, obsFlags, cloudAmount, &
                                     obsChannels, sensorIndex, obsDate, obsHour, burpFileSatId, &
                                     maxAngleReached, topographicData, nonCorrectedData, isTbPresent, &
                                     isClearSky, strayLight, goesMidi, isToAssim, ompOutOfRange) 

      !###############################################################################
      ! STEP 2) Controle de qualite des CSR. Data QC flags (obsFlags) are modified here!
      !###############################################################################

      call csrbg_csrCheckQc(obsFlags, categorierejet, maxAngleReached, topographicData, &
                            nonCorrectedData, isTbPresent, isClearSky, strayLight, goesMidi, isToAssim, &
                            ompOutOfRange)

      !###############################################################################
      ! STEP 3) Update Flags and obs in obsspace data
      !###############################################################################

      call csrbg_updateObsSpaceAfterQc(obsSpaceData, obsFlags, headerIndex, sensorIndex)

    end do HEADER 
    write(*,*) "Nombre de donnees rejetees"
    write(*,*) "Attention, une donnee peut etre rejetee pour plusieurs raisons"
    write(*,*) "Maxangle, straylight ou goesmid    " , categorieRejet(1)
    write(*,*) "Topographie                        " , categorieRejet(2)
    write(*,*) "TB sans correction de biais        " , categorieRejet(3)
    write(*,*) "Pas TB                             " , categorieRejet(4)
    write(*,*) "std(O-P)*sigma trop grand          " , categorieRejet(5)
    write(*,*) "canal non assimile                 " , categorieRejet(6)
    write(*,*) "Ciel non clair                     " , categorieRejet(7)
    write(*,*) "*******"

    call utl_tmg_stop(116)
  end subroutine csrbg_bgCheckCSR

  !--------------------------------------------------------------------------
  !  csrbg_readObsFromObsSpace
  !--------------------------------------------------------------------------

  subroutine csrbg_readObsFromObsSpace(obsSpaceData, headerIndex, obsTb, ompTb, satZenithAngle, obsFlags, cloudAmount, &
                                       obsChannels, sensorIndex, obsDate, obsHour, burpFileSatId, &
                                       maxAngleReached, topographicData, nonCorrectedData, isTbPresent, &
                                       isClearSky, strayLight, goesMidi, isToAssim, ompOutOfRange)

    !:Purpose:        copy headers and bodies from obsSpaceData object to arrays
    !                 compute some parameters from the read variables

    implicit None

    
    integer , intent(in)                 :: headerIndex         ! current header index
    real   , allocatable, intent(out)    :: obsTb(:)            ! brightness temperature (btyp=9248/9264,ele=12163)
    real   , allocatable, intent(out)    :: ompTb(:)            ! OMP values
    real   , allocatable, intent(out)    :: satZenithAngle(:)   ! satellite zenith angle (btyp=3072,ele=7024)
    integer, allocatable, intent(out)    :: obsFlags(:)         ! data flags
    integer, allocatable, intent(out)    :: cloudAmount(:)      ! data flags
    integer, allocatable, intent(out)    :: obsChannels(:)      ! Tb Channels 
    integer, allocatable, intent(out)    :: obsDate(:)          ! date YYYYMMDD
    integer, allocatable, intent(out)    :: obsHour(:)          ! Hour HHMM
    integer,              intent(out)    :: sensorIndex         ! find tvs_sensor index corresponding to current obs
    character(*),intent(out)             :: burpFileSatId       ! Platform Name
    type(struct_obs),     intent(inout)  :: obsSpaceData        ! obspaceData Object

    ! Data derived from elements read in obspace data

    logical, allocatable, intent(out)    :: maxAngleReached(:)     ! satellite angle exceed max angle at obs
    logical, allocatable, intent(out)    :: topographicData(:)  ! data flagged as topo data
    logical, allocatable, intent(out)    :: nonCorrectedData(:) ! data non corrected by bias corr
    logical, allocatable, intent(out)    :: isTbPresent(:)      ! non missing data
    logical, allocatable, intent(out)    :: isClearSky(:)       ! clear sky obs
    logical, allocatable, intent(out)    :: strayLight(:)       !
    logical, allocatable, intent(out)    :: goesMidi(:)         ! goes noon
    logical, allocatable, intent(out)    :: isToAssim(:)        ! is channel assimilable
    logical, allocatable, intent(out)    :: ompOutOfRange(:)    ! 
  
    ! Locals
    integer                              :: bodyIndex
    integer                              :: bodyIndexbeg
    integer                              :: headerCompt
    integer                              :: currentChannelNumber
    integer                              :: channelIndex
    integer                              :: actualnumChannel
    integer                              :: numObsToProcess
    integer                              :: iplatform
    integer                              :: instrum
    integer                              :: isat, iplatf
    integer                              :: instr
    integer                              :: mois
    integer                              :: jour
    integer                              :: heure
    integer                              :: midnight
    integer                              :: dataIndex
    integer                              :: indexSat
    logical                              :: sensorIndexFound
    logical                              :: indexSatFound
    real                                 :: errorThreshold



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
    
    if ( .not. sensorIndexFound ) call utl_abort('csrbg_readObsFromObsSpace: sensor Index not found')

    ! find actual Number of channels
    actualNumChannel = tvs_coefs(sensorIndex)%coef%fmv_ori_nchn
    

    headerCompt = 1
    numObsToProcess = 1
    ! Allocate Header elements
    call utl_reAllocate(obsDate, numObsToProcess)
    call utl_reAllocate(obsHour, numObsToProcess)
    call utl_reAllocate(satZenithAngle, numObsToProcess)

    ! Allocate Body elements
    call utl_reAllocate(obsTb, numObsToProcess*actualNumChannel)
    call utl_reAllocate(ompTb, numObsToProcess*actualNumChannel)
    call utl_reAllocate(obsFlags, numObsToProcess*actualNumChannel)
    call utl_reAllocate(cloudAmount, numObsToProcess*actualNumChannel)
    call utl_reAllocate(obsChannels, numObsToProcess*actualNumChannel)

    !initialization
    obsTb(:) = csrbg_realMissing
    ompTb(:) = csrbg_realMissing
    obsChannels(:) = csrbg_realMissing
    cloudAmount(:) = csrbg_realMissing
    !
    ! Header elements
    !
    burpFileSatId                      = obs_elem_c    ( obsSpaceData, 'STID' , headerIndex )
    satZenithAngle(headerCompt)        = obs_headElem_r( obsSpaceData, OBS_SZA, headerIndex )
    obsDate(headerCompt)               = obs_headElem_i( obsSpaceData, OBS_DAT, headerIndex )
    obsHour(headerCompt)               = obs_headElem_i( obsSpaceData, OBS_ETM, headerIndex )
    !
    ! Body elements
    !
    bodyIndexbeg        = obs_headElem_i( obsSpaceData, OBS_RLN, headerIndex )

    BODY: do bodyIndex =  bodyIndexbeg, bodyIndexbeg + obs_headElem_i( obsSpaceData, OBS_NLV, headerIndex ) - 1
      currentChannelNumber = nint(obs_bodyElem_r( obsSpaceData,  OBS_PPP, bodyIndex )) - & 
                             tvs_channelOffset(sensorIndex)
      obsTb(currentChannelNumber)          = obs_bodyElem_r( obsSpaceData,  OBS_VAR, bodyIndex )
      ompTb(currentChannelNumber)          = obs_bodyElem_r( obsSpaceData,  OBS_OMP, bodyIndex )
      obsFlags(currentChannelNumber)       = obs_bodyElem_i( obsSpaceData,  OBS_FLG, bodyIndex )
      cloudAmount(currentChannelNumber)    = obs_bodyElem_i( obsSpaceData, OBS_CLA, bodyIndex)

    end do BODY

    do channelIndex=1, actualNumChannel
      obsChannels(channelIndex)    = channelIndex+tvs_channelOffset(sensorIndex)
    end do


    ! Get index for satCloudCoverLimit for current current sat  
    indexSatFound = .false.
    do indexSat = 1, maxNumChan
      if (trim(burpSatName(indexSat)) == trim(burpFileSatId)) then
        indexSatFound = .true.
        exit
      end if 
    end do
    if ( .not. indexSatFound ) call utl_abort('csrbg_readObsFromObsSpace: Cloud Cover Limit Not & 
                                               Found for ' // trim(burpFileSatId))
    ! compute data for QC
    call utl_reAllocate(topographicData, numObsToProcess)
    call utl_reAllocate(maxAngleReached, numObsToProcess)
    call utl_reAllocate(strayLight, numObsToProcess)
    call utl_reAllocate(goesMidi, numObsToProcess)
    call utl_reAllocate(isTbPresent, numObsToProcess*actualNumChannel)
    call utl_reAllocate(nonCorrectedData, numObsToProcess*actualNumChannel)
    call utl_reAllocate(isClearSky, numObsToProcess*actualNumChannel)
    call utl_reAllocate(ompOutOfRange, numObsToProcess*actualNumChannel)
    call utl_reAllocate(isToAssim, numObsToProcess*actualNumChannel)
    
    isTbPresent(:) = .false.
    isClearSky(:) = .false.
    nonCorrectedData = .false.
    
    do channelIndex=1, numObsToProcess*actualNumChannel
      if (obsTb(channelIndex) /= csrbg_realMissing) isTbPresent(channelIndex) = .true.
      if (cloudAmount(channelIndex) /= csrbg_realMissing .and. & 
          cloudAmount(channelIndex) < satCloudCoverlimit(indexSat,channelIndex)) then
        isClearSky(channelIndex) = .true.
      end if 
      if (.not. btest(obsFlags(channelIndex), 6)) nonCorrectedData(channelIndex) = .true.
    end do

    topographicData(:) = .false.
    maxAngleReached(:) = .false.
    do dataIndex=1, numObsToProcess
      if (btest(obsFlags(dataIndex), 18)) topographicData(dataIndex) = .true.
      if (satZenithAngle(dataIndex) > 15250) maxAngleReached(dataIndex) = .true.
    end do

    !!! check O-P of Tb and set assim flag for each channel of satellite numsat using stat_iutilst
    !!! If channel not found in stats file, NO ASSIMILATE flag is set and O-P check is skipped
    ompOutOfRange(:) = .false.
    isToAssim(:) = .false.
    do  channelIndex=1, numObsToProcess*actualNumChannel
      if ( .not. oer_tovutil(obsChannels(channelIndex), sensorIndex) == 0) then
        isToAssim(channelIndex) = .true.
      end if
      errorThreshold = oer_toverrst(obsChannels(channelIndex), sensorIndex)*csrbg_ompThreshold
      if (abs(ompTb(channelIndex)) > errorThreshold) ompOutOfRange(channelIndex) = .true.
    end do
    
    strayLight(:) = .false.
    goesMidi(:) = .false.
    do  dataIndex = 1, numObsToProcess
      mois  = (obsDate(dataIndex)/100) - (obsDate(dataIndex)/10000)*100
      jour  = (obsDate(dataIndex)/100) -  (obsDate(dataIndex)/100)*100
      heure = obsHour(dataIndex)/100 
      if (burpFileSatId == '^METSAT7') then
        if ( heure == 21 ) straylight(dataIndex) = 1
        if ( (heure == 20) .or. (heure == 22) ) then
          if ( (jour + (mois-1)*30) > 20 .and. (jour + (mois-1)*30) < 135 ) straylight(dataIndex) = .true.
          if ( (jour + (mois-1)*30) > 210 .and. (jour + (mois-1)*30) < 300 ) straylight(dataIndex) = .true.
        end if
      end if
      if (burpFileSatId == '^GOES11' .or.  burpFileSatId == '^GOES15') midnight = 9
      if (burpFileSatId == '^GOES13' .or.  burpFileSatId == '^GOES14') midnight = 5
      if (burpFileSatId == '^GOES11' .or.  burpFileSatId == '^GOES15' .or. &
          burpFileSatId == '^GOES13' .or.  burpFileSatId == '^GOES14') then
        if ( (heure == (midnight - 1)) .or. (heure == midnight) ) goesMidi(dataIndex) = .true.
      end if 
    end do
    
  end subroutine csrbg_readObsFromObsSpace


  !--------------------------------------------------------------------------
  !  csrbg_csrCheckQc
  !--------------------------------------------------------------------------

  subroutine csrbg_csrCheckQc(obsFlags, categorieRejet, maxAngleReached, topographicData, &
                              nonCorrectedData, isTbPresent, isClearSky, strayLight, goesMidi, isToAssim, &
                              ompOutOfRange)

    !:Purpose:        Modify obsFlags

    implicit None


    integer,              intent(inout)     :: obsFlags(:)         ! obs Flags to update
    integer,              intent(inout)     :: categorieRejet(7)   ! the 7 categories of rejections
    logical,              intent(in)        :: maxAngleReached(:)  ! satellite angle exceed max angle at obs
    logical,              intent(in)        :: topographicData(:)  ! data flagged as topo data
    logical,              intent(in)        :: nonCorrectedData(:) ! data non corrected by bias corr
    logical,              intent(in)        :: isTbPresent(:)      ! non missing data
    logical,              intent(in)        :: isClearSky(:)       ! clear sky obs
    logical,              intent(in)        :: strayLight(:)       !
    logical,              intent(in)        :: goesMidi(:)         ! goes noon
    logical,              intent(in)        :: isToAssim(:)        ! is channel assimilable
    logical,              intent(in)        :: ompOutOfRange(:)    ! abs of omp greater than threshold 
    
    ! Locals
    integer                                 :: currentChannelNumber
    integer                                 :: channelIndex
    integer                                 :: numObsToProcess
    integer                                 :: numData
    integer                                 :: dataIndex


    numObsToProcess = 1
    currentChannelNumber = size(obsFlags)/numObsToProcess

    ! tests de controle de qualite sur les profils et canaux
    numData = 0
    do dataIndex = 1, numObsToProcess
      do channelIndex = 1, currentChannelNumber
        numData = numData+1
        !Maxangle & Straylight & Goesmid : Bit 7 et 9 pour tous les canaux
        if (maxAngleReached(dataIndex) == .true. .or. & 
            strayLight(dataIndex)      == .true. .or. & 
            goesMidi(dataIndex)        == .true. ) then 
          obsFlags(numData) = ibset(obsFlags(numData),7) 
          obsFlags(numData) = ibset(obsFlags(numData),9) 
          categorieRejet(1) =  categorieRejet(1) + 1
        end if

        !Topo : bit 9 et 18 sont déjà là provenant du EnVar
        if (topographicData(dataIndex) == .true.) then
          obsFlags(numData) = ibset(obsFlags(numData),9) 
          obsFlags(numData) = ibset(obsFlags(numData),18) 
          categorieRejet(2) =  categorieRejet(2) + 1
        end if 

        !Pas corrige : bit 11
        if (nonCorrectedData(numData) == .true.) then
          obsFlags(numData) = ibset(obsFlags(numData),11) 
          categorieRejet(3) =  categorieRejet(3) + 1
        end if 

        !Tbyela: bit  7 et 9 pour les canaux concernés
        if (isTbPresent(numData) == .false.) then 
          obsFlags(numData) = ibset(obsFlags(numData),7) 
          obsFlags(numData) = ibset(obsFlags(numData),9) 
          categorieRejet(4) =  categorieRejet(4) + 1
        end if

        !omp Out of range: Bit 9 et 16 pour les canaux concernés
        if (ompOutOfRange(numData) == .true.) then 
          obsFlags(numData) = ibset(obsFlags(numData),9)
          obsFlags(numData) = ibset(obsFlags(numData),16)
          categorieRejet(5) =  categorieRejet(5) + 1
        end if

        !Assim: Bit 11 pour les canaux concernés
        if (isToAssim(numData) == .false.) then 
          obsFlags(numData) = ibset(obsFlags(numData),11)
          categorieRejet(6) =  categorieRejet(6) + 1
        end if 

        !Clearsky : Bit 7 & Bit 9
        if (isClearSky(numData) == .false.) then 
          obsFlags(numData) = ibset(obsFlags(numData),7)
          obsFlags(numData) = ibset(obsFlags(numData),9)
          categorieRejet(7) =  categorieRejet(7) + 1
        end if  
      end do
    end do    

  end subroutine csrbg_csrCheckQc

  !--------------------------------------------------------------------------
  ! csrbg_updateObsSpaceAfterQc
  !--------------------------------------------------------------------------
  subroutine csrbg_updateObsSpaceAfterQc(obsSpaceData, obsFlags, headerIndex, sensorIndex)

    !:Purpose:      Update obspacedata variables (obstTB and obs flags) after QC
    implicit None

    !Arguments
    type(struct_obs),     intent(inout)     :: obsSpaceData           ! obspaceData Object
    integer,              intent(in)        :: obsFlags(:)            ! data flags
    integer,              intent(in)        :: sensorIndex            ! sensor Index 
    integer,              intent(in)        :: headerIndex            ! sensor Index 
    ! Locals
    integer                                 :: bodyIndex
    integer                                 :: bodyIndexbeg
    integer                                 :: currentChannelNumber

    bodyIndexbeg        = obs_headElem_i( obsSpaceData, OBS_RLN, headerIndex )

    BODY: do bodyIndex =  bodyIndexbeg, bodyIndexbeg + obs_headElem_i( obsSpaceData, OBS_NLV, headerIndex ) - 1
      currentChannelNumber=nint(obs_bodyElem_r( obsSpaceData,  OBS_PPP, bodyIndex ))-tvs_channelOffset(sensorIndex)
      call obs_bodySet_i(obsSpaceData, OBS_FLG, bodyIndex, obsFlags(currentChannelNumber))
    end do BODY

  end subroutine csrbg_updateObsSpaceAfterQc


end module bgckcsr_mod


















