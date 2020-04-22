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
  use mpi_mod
  use codePrecision_mod
  use bufr_mod

  implicit none
  save
  private
  integer, parameter :: nPhases=3, nLevels=5, nStationMax=100000 
  integer, parameter :: nStationMaxGP=10000
  integer  :: nb_aircraft_bias, nb_gps_bias
  real, allocatable  :: corrects_TT(:,:,:)
  real, allocatable  :: corrects_ZTD(:)
  character(len=9), allocatable :: aircraft_ID(:), gps_stn(:)
  
  ! Bias correction files (must be in program working directory)
  character(len=8), parameter :: ai_bcfile = "ai_bcors", gp_bcfile = "gp_bcors"

  public               :: bcc_readConfig, bcc_applyAIBcor, bcc_applyGPBcor
  
  integer, external    :: fnom, fclos
  
  logical :: aiBiasActive, gpBiasActive, aiRevOnly, gpRevOnly, gpRejBit11
  namelist /nambiasconv/ aiBiasActive, gpBiasActive, aiRevOnly, gpRevOnly, gpRejBit11
  
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
    gpRejBit11   = .true.   ! GP: Set data flag bit 11 (reject for assim) for any uncorrected data
    ! read in the namelist NAMBIASCONV
    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=nambiasconv,iostat=ierr)
    if ( ierr /= 0 .and. mpi_myid == 0 )  &
         write(*,*) 'bcc_readConfig: WARNING: Error reading namelist, assume it will not be used!'
    if ( mpi_myid == 0 ) write(*,nml=nambiasconv)
    ierr = fclos(nulnam)
    
  end subroutine bcc_readConfig

  !-----------------------------------------------------------------------
  ! bcc_readAIBcor
  !-----------------------------------------------------------------------
  subroutine bcc_readAIBcor(file_cor)
    !
    ! :Purpose: Read TT bias corrections 
    !     The first line of the file is the number of aircraft plus one.
    !     The rest of the file gives 15 values of Mean O-A for each aircraft, with each (AC,value) line written in format "a9,1x,f6.2".
    !     The order is the same as what is written by genbiascorr script genbc.aircraft_bcor.py.
    !     The first "aircraft" (AC name = BULKBCORS) values are the bulk corrections by layer for All-AC (first 5 values), AIREP/ADS 
    !     (second 5 values) and AMDAR/BUFR (last 5 values).
    !     Missing value = 99.0.
    !
    implicit none
    !Arguments:
    character(len=*), intent(in) :: file_cor

    !Locals:
    integer :: ierr, nulcoeff
    integer :: stationIndex, phaseIndex, levelIndex
    real    :: corr_ligne
    character(len=9) :: id_ligne

    if (aiRevOnly) return
    
    nulcoeff = 0
    ierr = fnom(nulcoeff, file_cor, 'FMT+R/O', 0)
    if ( ierr /= 0 ) then
      call utl_abort('bcc_readAIBcor: unable to open airplanes bias correction file ' // file_cor )
    end if
    read (nulcoeff, '(i5)', iostat=ierr ) nb_aircraft_bias
    if ( ierr /= 0 ) then
      call utl_abort('bcc_readAIBcor: error 1 while reading airplanes bias correction file ' // file_cor )
    end if
    
    allocate( corrects_TT(nStationMax,nPhases,nLevels) )
    allocate( aircraft_ID(nStationMax) )

    corrects_TT(:,:,:) =  MPC_missingValue_R8
   
    do stationIndex=1,nb_aircraft_bias
      do phaseIndex=1,3
        do levelIndex=1,5
          read (nulcoeff, *, iostat=ierr) id_ligne,corr_ligne
          if ( ierr /= 0 ) then
            call utl_abort('bcc_readAIBcor: error 2 while reading airplanes bias correction file ' // file_cor )
          end if
          if (corr_ligne == 99.) corr_ligne = MPC_missingValue_R8
          corrects_TT(stationIndex,phaseIndex,levelIndex) = corr_ligne
          aircraft_ID(stationIndex)                       = id_ligne
          !print*, stationIndex, phaseIndex, levelIndex, aircraft_ID(stationIndex), corrects_TT(stationIndex,phaseIndex,levelIndex)
        end do
      end do
    end do
    ierr = fclos(nulcoeff)

    ! Check for bulk bias corrections at start of file
    if ( aircraft_ID(1) /= "BULKBCORS" ) then
      call utl_abort('bcc_readAIBcor: ERROR: Bulk bias corrections are missing in bias correction file!' )
    end if

  end subroutine bcc_readAIBcor

  !-----------------------------------------------------------------------
  ! bcc_applyAIBcor
  !-----------------------------------------------------------------------
  subroutine bcc_applyAIBcor(obsSpaceData)
    !
    ! :Purpose:  to fill OBS_BCOR column of ObsSpaceData body with aircraft TT bias correction 
    !
    implicit none
    !Arguments:
    type(struct_obs)        :: obsSpaceData
    !Locals:
    integer  :: headerIndex, bodyIndex, codtyp
    integer  :: flag, phase, bufrCode
    integer  :: phaseIndex, levelIndex, stationIndex, stationNumber
    integer  :: n_cor_ac,  n_cor_bk
    real(8)  :: corr, tt, oldCorr, pressure
    character(len=9) :: stnid, tempo1, tempo2

    if ( .not. aiBiasActive ) return

    write(*,*) "bcc_applyAIBcor: start"

    if ( .not.aiRevOnly ) call bcc_readAIBcor(ai_bcfile)
    
    n_cor_ac = 0
    n_cor_bk = 0

    call obs_set_current_header_list(obsSpaceData,'AI')
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if ( headerIndex < 0 ) exit HEADER
      
      call obs_set_current_body_list(obsSpaceData, headerIndex)

      BODY: do

        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if ( bodyIndex < 0 ) exit BODY

        if ( obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) /= obs_assimilated ) cycle BODY 

        bufrCode = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex )

        if ( bufrCode == BUFR_NETT) then
          tt      = obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex )
          flag    = obs_bodyElem_i(obsSpaceData, OBS_FLG, bodyIndex )
          oldCorr = obs_bodyElem_r(obsSpaceData, OBS_BCOR, bodyIndex )
          corr = MPC_missingValue_R8
          
          if ( tt /= real(MPC_missingValue_R8,OBS_REAL) ) then
          
            if ( btest(flag, 6) .and. oldCorr /= real(MPC_missingValue_R8,OBS_REAL) ) then
               tt = tt + oldCorr
               flag = ibclr(flag, 6)
            end if
            if (aiRevOnly) corr = 0.0
             
            IF ( .not.aiRevOnly ) THEN
                
            pressure = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex ) * MPC_MBAR_PER_PA_R8

            ! Get level index and current (mar 2020) static bulk corrections applied to AI TT data at derivate stage
            if ( (pressure <= 1100.d0) .and. (pressure > 700.d0) ) then
              corr = 0.0
              levelIndex = 5
            else if ( (pressure <= 700.d0)  .and. (pressure > 500.d0) ) then
              corr = 0.1
              levelIndex = 4
            else if ( (pressure <= 500.d0)  .and. (pressure > 400.d0) ) then
              corr = 0.2
              levelIndex = 3
            else if ( (pressure <= 400.d0)  .and. (pressure > 300.d0) ) then
              corr = 0.3
              levelIndex = 2
            else if ( (pressure <= 300.d0)  .and. (pressure > 100.d0) ) then
              corr = 0.5
              levelIndex = 1
            else 
              levelIndex = 0
              corr = 0.0
            end if       
 
            codtyp = obs_headElem_i(obsSpaceData, OBS_ITY, headerIndex)

            ! Default bulk corrections read from bcor file (applied if dynamic corrections are not availble for the aircraft)
            if ( codtyp == 128 .or. codtyp == 177 ) then  ! AIREP/ADS
              phaseIndex = 2
            else  ! AMDAR/BUFR
              phaseIndex = 3
            end if
            if ( levelIndex /= 0 ) then
              if ( corrects_TT(1,phaseIndex,levelIndex) /= MPC_missingValue_R8) corr = corrects_TT(1,phaseIndex,levelIndex)
              n_cor_bk = n_cor_bk + 1
            end if

            headerIndex = obs_bodyElem_i(obsSpaceData, OBS_HIND, bodyIndex  )
            stnid = trim( obs_elem_c(obsSpaceData,'STID',headerIndex) )

            ! on verifie si la station est dans le dictionnaire du fichier de correction de biais
            !---------------------------------------------------------------------------------
            stationNumber = 0
            do stationIndex = 1, nb_aircraft_bias
              tempo1 = trim(aircraft_ID(stationIndex))
              tempo2 = trim(stnid)
              if ( tempo2(2:9) == tempo1(1:8) ) stationNumber = stationIndex
            end do
            
            phase =  obs_headElem_i(obsSpaceData, OBS_PHAS, headerIndex  )

            ! If the aircraft is in the bias correction file, get the new correction
            ! corrects_TT(stationNumber,phaseIndex,levelIndex) where
            !     stationNumber = index for this AC in bias correction file (0 if not found)
            !     phaseIndex   = index for the 3 phases of flight (level, asc, desc)
            !     levelIndex   = index for the 5 layers (100-300, 300-400,400-500,500-700,700-1100)
            !  and use it instead of bulk value (if it is not missing value).
            if (stationNumber /= 0) then 
              phaseIndex = 0
              if ( phase == 3 ) phaseIndex = 1 ! level
              if ( phase == 5 ) phaseIndex = 2 ! ascent
              if ( phase == 6 ) phaseIndex = 3 ! descent
              if (levelIndex /= 0 .and. phaseIndex /= 0) then
                if ( corrects_TT(stationNumber,phaseIndex,levelIndex) /= MPC_missingValue_R8 ) then
                   corr = corrects_TT(stationNumber,phaseIndex,levelIndex)
                   n_cor_ac = n_cor_ac + 1
                   n_cor_bk = n_cor_bk - 1
                end if
              end if
            end if
            
            ! Apply the bias correction (bulk or new) and set the "bias corrected" bit in TT data flag ON
            tt = tt - corr
            flag = ibset(flag, 6)
            
            END IF
            
          end if

          call obs_bodySet_r( obsSpaceData, OBS_BCOR, bodyIndex, corr )
          call obs_bodySet_r( obsSpaceData, OBS_VAR , bodyIndex, tt   )
          call obs_bodySet_i( obsSpaceData, OBS_FLG , bodyIndex, flag )        

        end if
        
      end do BODY
    end do HEADER
    
    if ( n_cor_bk+n_cor_ac /= 0 ) then
       write (*, '(a50, i10)' ) "bcc_applyAIBcor: Number of obs with TT bulk correction  = ", n_cor_bk
       write (*, '(a50, i10)' ) "bcc_applyAIBcor: Number of obs with TT tail correction  = ", n_cor_ac
    else
       write(*,*) "No AI data found"
    end if
    
    if (allocated(corrects_TT)) deallocate(corrects_TT)
    if (allocated(aircraft_ID)) deallocate(aircraft_ID)
    
    write(*,*) "bcc_applyAIBcor: end"
    
  end subroutine bcc_applyAIBcor

  !-----------------------------------------------------------------------
  ! bcc_readGPBcor
  !-----------------------------------------------------------------------
  subroutine bcc_readGPBcor(file_cor)
    !
    ! :Purpose: Read GB-GPS bias corrections (mean O-A by station)
    !     Missing value = -999.00
    !
    implicit none
    !Arguments:
    character(len=*), intent(in) :: file_cor
    !Locals:
    integer :: ierr, nulcoeff
    integer :: stationIndex
    real    :: corr_ligne
    character(len=9) :: id_ligne

    if (gpRevOnly) return
    
    nulcoeff = 0
    ierr = fnom(nulcoeff, file_cor, 'FMT+R/O', 0)
    if ( ierr /= 0 ) then
      call utl_abort('bcc_readGPBcor: unable to open GB-GPS bias correction file ' // file_cor )
    end if
    read (nulcoeff, '(i5)', iostat=ierr ) nb_gps_bias
    if ( ierr /= 0 ) then
      call utl_abort('bcc_readGPBcor: error 1 while reading GB-GPS bias correction file ' // file_cor )
    end if

    allocate( corrects_ZTD(nStationMaxGP) )
    allocate( gps_stn(nStationMaxGP)  )
    
    corrects_ZTD(:) =  MPC_missingValue_R8
    
    do stationIndex=1,nb_gps_bias
       read (nulcoeff, *, iostat=ierr) id_ligne,corr_ligne
       if ( ierr /= 0 ) then
          call utl_abort('bcc_readGPBcor: error 2 while reading GB-GPS bias correction file ' // file_cor )
       end if
       if (corr_ligne /= -999.00) corrects_ZTD(stationIndex) = -corr_ligne
       gps_stn(stationIndex) = id_ligne
    end do
    ierr = fclos(nulcoeff)
    
  end subroutine bcc_readGPBcor

  !-----------------------------------------------------------------------
  ! bcc_applyGPBcor
  !-----------------------------------------------------------------------
  subroutine bcc_applyGPBcor(obsSpaceData)
    !
    ! :Purpose:  to fill OBS_BCOR column of ObsSpaceData body with GB-GPS ZTD bias correction 
    !
    implicit none
    !Arguments:
    type(struct_obs)  :: obsSpaceData
    !Locals:
    integer  :: headerIndex, bodyIndex
    integer  :: flag, bufrCode
    integer  :: stationIndex, stationNumber
    integer  :: n_cor
    real(8)  :: corr, ztd, oldCorr
    character(len=9) :: stnid, tempo1, tempo2

    if ( .not. gpBiasActive ) return

    write(*,*) "bcc_applyGPBcor: start"

    if (.not.gpRevOnly) call bcc_readGPBcor(gp_bcfile)
    
    n_cor = 0

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
          
          if ( ztd /= real(MPC_missingValue_R8,OBS_REAL) ) then  
            
            ! Remove any previous bias correction
            if ( btest(flag, 6) .and. oldCorr /= real(MPC_missingValue_R8,OBS_REAL) ) then
              ztd = ztd - oldCorr
              flag = ibclr(flag, 6)
            end if

            IF (.not. gpRevOnly ) THEN
               headerIndex = obs_bodyElem_i(obsSpaceData, OBS_HIND, bodyIndex  )
               stnid = trim( obs_elem_c(obsSpaceData,'STID',headerIndex) )

               ! on verifie si la station est dans le dictionnaire du fichier de correction de biais
               !---------------------------------------------------------------------------------
               stationNumber = 0
               do stationIndex = 1, nb_gps_bias
                 tempo1 = trim(gps_stn(stationIndex))
                 tempo2 = trim(stnid)
                 if ( tempo2 == tempo1 ) stationNumber = stationIndex
               end do
               
               if (stationNumber /= 0) then 
                 corr = corrects_ZTD(stationNumber)
               end if
               
               ! Apply the bias correction and set the "bias corrected" bit in ZTD data flag ON
               if ( corr /= MPC_missingValue_R8 ) then
                 ztd = ztd + corr
                 n_cor = n_cor + 1
                 flag = ibset(flag, 6)
               else 
                 corr = 0.0
               end if
            END IF

          end if

          call obs_bodySet_r( obsSpaceData, OBS_BCOR, bodyIndex, corr )
          call obs_bodySet_r( obsSpaceData, OBS_VAR , bodyIndex, ztd  )
          ! Flag uncorrected data for rejection
          if (.not. btest(flag, 6) .and. gpRejBit11) flag = ibset(flag, 11)
          call obs_bodySet_i( obsSpaceData, OBS_FLG , bodyIndex, flag )
           
        end if
        
      end do BODY
    end do HEADER
    
    if ( n_cor /= 0 ) then
      write (*, '(a50, i10)' ) "bcc_applyGPBcor: Number of ZTD observations corrected  = ", n_cor
    else 
      write(*,*) "No GP data bias corrections made"
    end if
    
    if (allocated(corrects_ZTD)) deallocate( corrects_ZTD )
    if (allocated(gps_stn))      deallocate( gps_stn  )
    
    write(*,*) "bcc_applyGPBcor: end"
    
  end subroutine bcc_applyGPBcor
    

end MODULE biasCorrectionConv_mod

