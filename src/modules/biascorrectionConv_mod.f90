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

MODULE biasCorrectionConv_mod
  ! MODULE biasCorrectionConv_mod (prefix="bcc" category='1. High-level functionality')
  !
  ! :Purpose: Performs bias correction for conventional observations
  !
  use utilities_mod
  use obsSpaceData_mod
  use MathPhysConstants_mod
  use midasMpi_mod
  use codePrecision_mod
  use bufr_mod
  use codtyp_mod

  implicit none
  save
  private
  public               :: bcc_applyAIBcor, bcc_applyGPBcor
  public               :: bcc_biasActive

  integer, parameter :: nPhases=3, nLevels=5, nAircraftMax=100000 
  integer, parameter :: nStationMaxGP=10000
  integer, parameter :: phaseLevel   = 3
  integer, parameter :: phaseAscent  = 5
  integer, parameter :: phaseDescent = 6
  integer, parameter :: phaseLevelIndex   = 1
  integer, parameter :: phaseAscentIndex  = 2
  integer, parameter :: phaseDescentIndex = 3


  integer  :: nbAircrafts, nbGpStations
  real(8), allocatable  :: ttCorrections(:,:,:)
  real(8), allocatable  :: ztdCorrections(:)
  character(len=9), allocatable :: aircraftIds(:), gpsStations(:)
  
  logical :: bcc_aiBiasActive, bcc_gpBiasActive
  logical :: initialized = .false.
  
  ! Bias correction files (must be in program working directory)
  character(len=8), parameter :: aiBcFile = "ai_bcors", gpBcFile = "gp_bcors"

  integer, external    :: fnom, fclos
  
  logical :: aiBiasActive ! Control if bias correction is applied to aircraft data
  logical :: gpBiasActive ! Control if bias correction is applied to ground bases GPS data
  logical :: aiRevOnly    ! Don't apply new correction but simply reverse any old corrections for AI
  logical :: gpRevOnly    ! Don't apply new correction but simply reverse any old corrections for GP
  namelist /nambiasconv/ aiBiasActive, gpBiasActive, aiRevOnly, gpRevOnly
  
CONTAINS
 
  !-----------------------------------------------------------------------
  ! bcc_readConfig
  !-----------------------------------------------------------------------
  subroutine bcc_readConfig()
    !
    ! :Purpose: Read nambiasconv namelist section
    !
    implicit none
    !Locals:
    integer  :: ierr,nulnam
  
    ! set default values for namelist variables
    aiBiasActive = .false.  ! bias correct AI data (TT)
    gpBiasActive = .false.  ! bias correct GP data (ZTD)
    aiRevOnly    = .false.  ! AI: don't apply new correction but simply reverse any old corrections
    gpRevOnly    = .false.  ! GP: don't apply new correction but simply reverse any old corrections
    ! read in the namelist NAMBIASCONV
    if ( utl_isNamelistPresent('nambiasconv','./flnml') ) then
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam,nml=nambiasconv,iostat=ierr)
      if ( ierr /= 0 )  call utl_abort('bcc_readConfig: Error reading namelist section nambiasconv')
      if ( mmpi_myid == 0 ) write(*,nml=nambiasconv)
      ierr = fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'bcc_readconfig: nambiasconv is missing in the namelist. The default value will be taken.'
    end if
    
    bcc_aiBiasActive = aiBiasActive
    bcc_gpBiasActive = gpBiasActive
    
    initialized = .true.
    
    
  end subroutine bcc_readConfig

  !-----------------------------------------------------------------------
  ! bcc_readAIBiases
  !-----------------------------------------------------------------------
  subroutine bcc_readAIBiases(biasEstimateFile)
    !
    ! :Purpose: Read aircraft (AI) TT bias estimates from bias file and fill bias correction array ttCorrections. 
    !           The first line of the file is the number of aircraft plus one.
    !           The rest of the file gives 15 values of Mean O-A for each aircraft, with each (AC,value) line written in format "a9,1x,f6.2".
    !           The order is the same as what is written by genbiascorr script genbc.aircraft_bcor.py.
    !           The first "aircraft" (AC name = BULKBCORS) values are the bulk biases by layer for All-AC (first 5 values), AIREP/ADS 
    !           (second 5 values) and AMDAR/BUFR (last 5 values).
    !           Missing value = 99.0.
    !
    implicit none
    !Arguments:
    character(len=*), intent(in) :: biasEstimateFile

    !Locals:
    integer :: ierr, nulcoeff
    integer :: stationIndex, phaseIndex, levelIndex
    real(8) :: biasEstimate, correctionValue
    character(len=9) :: stationId

    if (.not.initialized) call bcc_readConfig()
    
    if ( aiRevOnly ) return
    
    nulcoeff = 0
    ierr = fnom(nulcoeff, biasEstimateFile, 'FMT+R/O', 0)
    if ( ierr /= 0 ) then
      call utl_abort('bcc_readAIBiases: unable to open airplanes bias correction file ' // biasEstimateFile )
    end if
    read (nulcoeff, '(i5)', iostat=ierr ) nbAircrafts
    if ( ierr /= 0 ) then
      call utl_abort('bcc_readAIBiases: error 1 while reading airplanes bias correction file ' // biasEstimateFile )
    end if
    
    allocate( ttCorrections(nAircraftMax,nPhases,nLevels) )
    allocate( aircraftIds(nAircraftMax) )

    ttCorrections(:,:,:) =  MPC_missingValue_R8
   
    do stationIndex=1,nbAircrafts
      do phaseIndex=1,3
        do levelIndex=1,5
          read (nulcoeff, *, iostat=ierr) stationId, biasEstimate
          if ( ierr /= 0 ) then
            call utl_abort('bcc_readAIBiases: error 2 while reading airplanes bias correction file ' // biasEstimateFile )
          end if
          if ( biasEstimate == 99.0D0 ) then
            correctionValue = MPC_missingValue_R8
          else
            correctionValue = -1.0D0*biasEstimate
          end if
          ttCorrections(stationIndex,phaseIndex,levelIndex) = correctionValue
          aircraftIds(stationIndex)                         = stationId
          !print*, stationIndex, phaseIndex, levelIndex, aircraftIds(stationIndex), ttCorrections(stationIndex,phaseIndex,levelIndex)
        end do
      end do
    end do
    ierr = fclos(nulcoeff)

    ! Check for bulk bias corrections at start of file
    if ( aircraftIds(1) /= "BULKBCORS" ) then
      call utl_abort('bcc_readAIBiases: ERROR: Bulk bias corrections are missing in bias correction file!' )
    end if

  end subroutine bcc_readAIBiases

  !-----------------------------------------------------------------------
  ! bcc_applyAIBcor
  !-----------------------------------------------------------------------
  subroutine bcc_applyAIBcor(obsSpaceData)
    !
    ! :Purpose:  to apply aircraft (AI) bias corrections to observations in ObsSpaceData
    !
    implicit none
    !Arguments:
    type(struct_obs)        :: obsSpaceData
    !Locals:
    integer  :: headerIndex, bodyIndex, codtyp
    integer  :: flag, phase, bufrCode
    integer  :: phaseIndex, levelIndex, stationIndex, stationNumber
    integer  :: countTailCorrections,  countBulkCorrections
    integer  :: headerFlag
    real(8)  :: corr, tt, oldCorr, pressure
    character(len=9) :: stnid, stnId1, stnId2

    if (.not.initialized) call bcc_readConfig()
    
    if ( .not. aiBiasActive ) return

    write(*,*) "bcc_applyAIBcor: start"

    if ( .not.aiRevOnly ) call bcc_readAIBiases(aiBcFile)
    
    countTailCorrections = 0
    countBulkCorrections = 0

    call obs_set_current_header_list(obsSpaceData,'AI')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if ( headerIndex < 0 ) exit HEADER
      
      headerFlag = obs_headElem_i(obsSpaceData, OBS_ST1, headerIndex )
      
      call obs_set_current_body_list(obsSpaceData, headerIndex)

      BODY: do

        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if ( bodyIndex < 0 ) exit BODY

        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) cycle BODY 

        bufrCode = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex )

        if ( bufrCode == BUFR_NETT ) then
          tt      = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
          flag    = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex )
          oldCorr = obs_bodyElem_r(obsSpaceData, OBS_BCOR, bodyIndex )
          corr = MPC_missingValue_R8
          
          if ( tt /= MPC_missingValue_R8 ) then
          
            if ( btest(flag, 6) .and. oldCorr /= MPC_missingValue_R8 ) then
              if ( btest(headerFlag, 15) ) then
                tt = tt - oldCorr
              else
                tt = tt + oldCorr
              end if
              flag = ibclr(flag, 6)
            end if
            if ( aiRevOnly ) corr = 0.0D0
             
            if ( .not.aiRevOnly ) then
                
              pressure = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex ) * MPC_MBAR_PER_PA_R8

              ! Get level index and current (mar 2020) static bulk corrections applied to AI TT data at derivate stage
              if ( (pressure <= 1100.d0) .and. (pressure > 700.d0) ) then
                corr = 0.0D0
                levelIndex = 5
              else if ( (pressure <= 700.d0)  .and. (pressure > 500.d0) ) then
                corr = -0.1D0
                levelIndex = 4
              else if ( (pressure <= 500.d0)  .and. (pressure > 400.d0) ) then
                corr = -0.2D0
                levelIndex = 3
              else if ( (pressure <= 400.d0)  .and. (pressure > 300.d0) ) then
                corr = -0.3D0
                levelIndex = 2
              else if ( (pressure <= 300.d0)  .and. (pressure > 100.d0) ) then
                corr = -0.5D0
                levelIndex = 1
              else 
                levelIndex = 0
                corr = 0.0D0
              end if
 
              codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)

              ! Default bulk corrections read from bcor file (applied if dynamic corrections are not availble for the aircraft)
              select case(  trim( codtyp_get_name(codtyp) ) )
                case('airep','ads')
                  phaseIndex = phaseAscentIndex
                case('amdar','acars')
                  phaseIndex = phaseDescentIndex
                case default
                  write(*,*) 'bcc_applyAIBcor: codtyp=', codtyp
                  call utl_abort('bcc_applyAIBcor: unknown codtyp') 
              end select

              if ( levelIndex /= 0 ) then
                if ( ttCorrections(1,phaseIndex,levelIndex) /= MPC_missingValue_R8 ) corr = ttCorrections(1,phaseIndex,levelIndex)
                countBulkCorrections = countBulkCorrections + 1
              end if

              headerIndex = obs_bodyElem_i(obsSpaceData, OBS_HIND, bodyIndex  )
              stnid = trim( obs_elem_c(obsSpaceData,'STID',headerIndex) )

              ! on verifie si la station est dans le dictionnaire du fichier de correction de biais
              !---------------------------------------------------------------------------------
              stationNumber = 0
              stnId2 = trim(stnid)
              do stationIndex = 1, nbAircrafts
                stnId1 = trim(aircraftIds(stationIndex))
                if ( stnId2(2:9) == stnId1(1:8) ) stationNumber = stationIndex
              end do
            
              phase =  obs_headElem_i(obsSpaceData, OBS_PHAS, headerIndex  )

              ! If the aircraft is in the bias correction file, get the new correction
              ! ttCorrections(stationNumber,phaseIndex,levelIndex) where
              !     stationNumber = index for this AC in bias correction file (0 if not found)
              !     phaseIndex   = index for the 3 phases of flight (level, asc, desc)
              !     levelIndex   = index for the 5 layers (100-300, 300-400,400-500,500-700,700-1100)
              !  and use it instead of bulk value (if it is not missing value).
              if ( stationNumber /= 0 ) then 
                phaseIndex = 0
                if ( phase == phaseLevel   ) phaseIndex = phaseLevelIndex
                if ( phase == phaseAscent  ) phaseIndex = phaseAscentIndex
                if ( phase == phaseDescent ) phaseIndex = phaseDescentIndex
                if ( levelIndex /= 0 .and. phaseIndex /= 0 ) then
                  if ( ttCorrections(stationNumber,phaseIndex,levelIndex) /= MPC_missingValue_R8 ) then
                    corr = ttCorrections(stationNumber,phaseIndex,levelIndex)
                    countTailCorrections = countTailCorrections + 1
                    countBulkCorrections = countBulkCorrections - 1
                  end if
                end if
              end if
            
              ! Apply the bias correction (bulk or new) and set the "bias corrected" bit in TT data flag ON
              tt = tt + corr
              flag = ibset(flag, 6)
            
            end if
            
          end if

          call obs_bodySet_r( obsSpaceData, OBS_BCOR, bodyIndex, corr )
          call obs_bodySet_r( obsSpaceData, OBS_VAR , bodyIndex, tt   )
          call obs_bodySet_i( obsSpaceData, OBS_FLG , bodyIndex, flag )        

        end if
        
      end do BODY
      
      headerFlag = ibset(headerFlag, 15)
      call obs_headSet_i( obsSpaceData, OBS_ST1, headerIndex, headerFlag )
      
    end do HEADER
    
    if ( countBulkCorrections + countTailCorrections /= 0 ) then
      write (*, '(a50, i10)' ) "bcc_applyAIBcor: Number of obs with TT bulk correction  = ", countBulkCorrections
      write (*, '(a50, i10)' ) "bcc_applyAIBcor: Number of obs with TT tail correction  = ", countTailCorrections
    else
      write(*,*) "bcc_applyAIBcor: No AI data found"
    end if
    
    if ( allocated(ttCorrections) ) deallocate(ttCorrections)
    if ( allocated(aircraftIds)   ) deallocate(aircraftIds)
    
    write(*,*) "bcc_applyAIBcor: end"
    
  end subroutine bcc_applyAIBcor

  !-----------------------------------------------------------------------
  ! bcc_readGPBiases
  !-----------------------------------------------------------------------
  subroutine bcc_readGPBiases(biasEstimateFile)
    !
    ! :Purpose: Read GB-GPS bias estimates (mean ZTD O-A [mm] by station) and fill bias correction array ztdCorrections.
    !           Missing value = -999.00
    !
    implicit none
    !Arguments:
    character(len=*), intent(in) :: biasEstimateFile
    !Locals:
    integer :: ierr, nulcoeff
    integer :: stationIndex
    real(8) :: biasEstimate
    character(len=9) :: stationId

    if (.not.initialized) call bcc_readConfig()
    
    if ( gpRevOnly ) return
    
    nulcoeff = 0
    ierr = fnom(nulcoeff, biasEstimateFile, 'FMT+R/O', 0)
    if ( ierr /= 0 ) then
      call utl_abort('bcc_readGPBiases: unable to open GB-GPS bias correction file ' // biasEstimateFile )
    end if
    read (nulcoeff, '(i5)', iostat=ierr ) nbGpStations
    if ( ierr /= 0 ) then
      call utl_abort('bcc_readGPBiases: error 1 while reading GB-GPS bias correction file ' // biasEstimateFile )
    end if

    allocate( ztdCorrections(nStationMaxGP) )
    allocate( gpsStations(nStationMaxGP)  )
    
    ztdCorrections(:) =  MPC_missingValue_R8
    
    do stationIndex=1,nbGpStations
       read (nulcoeff, *, iostat=ierr) stationId, biasEstimate
       if ( ierr /= 0 ) then
          call utl_abort('bcc_readGPBiases: error 2 while reading GB-GPS bias correction file ' // biasEstimateFile )
       end if
       if ( biasEstimate /= MPC_missingValue_R8 ) ztdCorrections(stationIndex) = -1.0D0*(biasEstimate/1000.0D0)  ! mm to m
       gpsStations(stationIndex) = stationId
    end do
    ierr = fclos(nulcoeff)
    
  end subroutine bcc_readGPBiases

  !-----------------------------------------------------------------------
  ! bcc_applyGPBcor
  !-----------------------------------------------------------------------
  subroutine bcc_applyGPBcor(obsSpaceData)
    !
    ! :Purpose:  to apply GB-GPS (GP) ZTD bias corrections to ZTD observations in ObsSpaceData
    !
    implicit none
    !Arguments:
    type(struct_obs)  :: obsSpaceData
    !Locals:
    integer  :: headerIndex, bodyIndex
    integer  :: flag, bufrCode
    integer  :: stationIndex, stationNumber
    integer  :: nbCorrected
    real(8)  :: corr, ztd, oldCorr
    character(len=9) :: stnid, stnId1, stnId2

    if (.not.initialized) call bcc_readConfig()
    
    if ( .not. gpBiasActive ) return

    write(*,*) "bcc_applyGPBcor: start"

    if ( .not.gpRevOnly ) call bcc_readGPBiases(gpBcFile)
    
    nbCorrected = 0

    call obs_set_current_header_list(obsSpaceData,'GP')
    
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if ( headerIndex < 0 ) exit HEADER
      
      call obs_set_current_body_list(obsSpaceData, headerIndex)

      BODY: do

        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if ( bodyIndex < 0 ) exit BODY

        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) cycle BODY 

        bufrCode = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex )

        if ( bufrCode == BUFR_NEZD ) then
          
          ztd     = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
          oldCorr = obs_bodyElem_r(obsSpaceData, OBS_BCOR, bodyIndex )
          flag    = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex )
          
          corr = MPC_missingValue_R8
          
          if ( ztd /= MPC_missingValue_R8 ) then  
            
            ! Remove any previous bias correction
            if ( btest(flag, 6) .and. oldCorr /= MPC_missingValue_R8 ) then
              ztd = ztd - oldCorr
              flag = ibclr(flag, 6)
            end if

            if ( .not. gpRevOnly ) then
              headerIndex = obs_bodyElem_i(obsSpaceData, OBS_HIND, bodyIndex  )
              stnid = trim( obs_elem_c(obsSpaceData,'STID',headerIndex) )

              ! on verifie si la station est dans le dictionnaire du fichier de correction de biais
              ! ---------------------------------------------------------------------------------
              stationNumber = 0
              stnId2 = trim(stnid)
              do stationIndex = 1, nbGpStations
                stnId1 = trim(gpsStations(stationIndex))
                if ( stnId2 == stnId1 ) stationNumber = stationIndex
              end do
               
              if (stationNumber /= 0) then 
                corr = ztdCorrections(stationNumber)
              end if
               
              ! Apply the bias correction and set the "bias corrected" bit in ZTD data flag ON
              if ( corr /= MPC_missingValue_R8 ) then
                ztd = ztd + corr
                nbCorrected = nbCorrected + 1
                flag = ibset(flag, 6)
              else 
                corr = 0.0D0
              end if
            else
              corr = 0.0D0
            end if

          end if

          call obs_bodySet_r( obsSpaceData, OBS_BCOR, bodyIndex, corr )
          call obs_bodySet_r( obsSpaceData, OBS_VAR , bodyIndex, ztd  )
          call obs_bodySet_i( obsSpaceData, OBS_FLG , bodyIndex, flag )
           
        end if
        
      end do BODY
    end do HEADER
    
    if ( nbCorrected /= 0 ) then
      write (*, '(a50, i10)' ) "bcc_applyGPBcor: Number of ZTD observations corrected  = ", nbCorrected
    else 
      write(*,*) "bcc_applyGPBcor: No GP data bias corrections made"
    end if
    
    if ( allocated(ztdCorrections) ) deallocate( ztdCorrections )
    if ( allocated(gpsStations) )    deallocate( gpsStations  )
    
    write(*,*) "bcc_applyGPBcor: end"
    
  end subroutine bcc_applyGPBcor


  !-----------------------------------------------------------------------
  ! bcc_biasActive
  !-----------------------------------------------------------------------
  logical function bcc_biasActive(obsFam)
    !
    ! :Purpose: returns if bias correction is active for the given conventional observation family
    !
    implicit none
    character(len=*),intent(in) :: obsFam

    if (.not.initialized) call bcc_readConfig()
    
    bcc_biasActive = (bcc_gpBiasActive .and. trim(obsFam) == 'GP') .or. (bcc_aiBiasActive .and. trim(obsFam) == 'AI')

  end function bcc_biasActive


end MODULE biasCorrectionConv_mod

