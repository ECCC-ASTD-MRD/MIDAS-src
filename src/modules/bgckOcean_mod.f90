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

module bgckOcean_mod
  ! MODULE bgckOcean_mod (prefix='ocebg' category='1. High-level functionality')
  !
  ! :Purpose: to perform ocean data background Check
  !
  use midasMpi_mod
  use utilities_mod
  use obsSpaceData_mod
  use columnData_mod
  use codtyp_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use statetocolumn_mod
  use bufr_mod
  use mathPhysConstants_mod
  use timeCoord_mod

  implicit none

  save
  private

  ! Public functions/subroutines
  public :: ocebg_bgCheckSST, ocebg_bgCheckSeaIce

  ! External functions
  integer, external :: fnom, fclos

  ! mpi topology
  integer           :: myLatBeg, myLatEnd
  integer           :: myLonBeg, myLonEnd
  integer           :: latPerPE, latPerPEmax, lonPerPE, lonPerPEmax

 ! namelist variables with setting default values
  character(len=20) :: timeInterpType_nl = 'NEAREST' ! 'NEAREST' or 'LINEAR'
  integer           :: numObsBatches     = 20        ! number of batches for calling interp setup
  logical           :: checkWinds        = .false.   ! if .true., check the winds for the last four
                                                     ! days to amplify the error in the zone of maximum wind speed
  integer           :: ndaysWinds        = 4         ! number of days in the 'winds' file to detect tropical storm (TS)
  integer           :: timeStepWinds     = 6         ! in hours, winds are available every timeStepWinds-hours
  integer           :: windForecastLeadtime = 6      ! in hours, lead time of wind forecast in the input file
  real(4)           :: minLatNH = 10.                ! min lat of Northern hemisphere latutude band where TS is detected
  real(4)           :: maxLatNH = 40.                ! max lat of Northern hemisphere latutude band where TS is detected
  real(4)           :: maxLatExceptionNH = 45.       ! max lat of Northern hemisphere latutude band
                                                     ! allows for TS to penetrate further North in some months 
  integer           :: nmonthsExceptionNH  = 0.      ! number of months where TS penetrates exceptionnally North
  character(len=3)  :: monthExceptionNH(12) = '   '  ! exceptional months where TS allowed to penetrated further North 
  real(4)           :: minLatSH = -35.               ! min lat of Southern hemisphere latutude band where TS is detected
  real(4)           :: maxLatSH = -10.               ! max lat of Southern hemisphere latutude band where TS is detected
  real(8)           :: smoothLenghtScale = 50000.    ! lenght scale. in m, to smooth the amplification error field
  real(8)           :: globalSelectCriteria(3) = (/5.d0, 25.d0, 30.d0/) ! global selection criteria
  logical           :: separateSelectCriteria = .false.! to apply a separate selection criteria on sea/inland waters 
                                                       ! and insitu/satellite data 
  real(8)           :: inlandWaterSelectCriteriaSatData(3) = (/5.d0, 25.d0, 30.d0/)
  real(8)           :: inlandWaterSelectCriteriaInsitu(3)  = (/5.d0, 25.d0, 30.d0/)
  real(8)           :: seaWaterSelectCriteriaSatData(3)    = (/5.d0, 25.d0, 30.d0/)
  real(8)           :: seaWaterSelectCriteriaInsitu(3)     = (/5.d0, 25.d0, 30.d0/) 
  real(4)           :: seaWaterThreshold = 0.1       ! threshold to distinguish inland water from sea water
  namelist /namOceanBGcheck/ timeInterpType_nl, numObsBatches, checkWinds, ndaysWinds, timeStepWinds, &
                             windForecastLeadtime, minLatNH, maxLatNH, maxLatExceptionNH, nmonthsExceptionNH, &
                             monthExceptionNH, minLatSH, maxLatSH, smoothLenghtScale, &
                             globalSelectCriteria, separateSelectCriteria, inlandWaterSelectCriteriaSatData, &
                             inlandWaterSelectCriteriaInsitu, seaWaterSelectCriteriaSatData, &
                             seaWaterSelectCriteriaInsitu, seaWaterThreshold 

  character(len=3), parameter :: months(12) = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

  contains

  !----------------------------------------------------------------------------------------
  ! ocebg_bgCheckSST
  !----------------------------------------------------------------------------------------
  subroutine ocebg_bgCheckSST(obsData, dateStamp, columnTrlOnTrlLev, hco)
    !
    ! :Purpose: to compute SST data background Check
    !

    implicit none

    ! Arguments:
    type(struct_obs)       , intent(inout)       :: obsData           ! obsSpaceData object
    integer                , intent(in)          :: dateStamp         ! date stamp
    type(struct_columnData), intent(inout)       :: columnTrlOnTrlLev ! column data on trl levels
    type(struct_hco)       , intent(in), pointer :: hco               ! horizontal trl grid

    ! Locals:
    type(struct_gsv)            :: stateVectorFGE        ! state vector containing std B estimation field
    type(struct_gsv)            :: stateVectorAmplFactor ! state vector for error amplification field
    real(4), pointer            :: stateVectorAmplFactor_ptr(:,:,:)
    type(struct_gsv)            :: stateVectorSeaWaterFraction ! statevector for sea water fraction
    integer                     :: nulnam, ierr, headerIndex, bodyIndex, obsFlag, obsVarno
    integer                     :: numberObs, numberObsRejected
    integer                     :: numberObsInsitu, numberObsInsituRejected, codeType  
    real(8)                     :: OER, OmP, FGE, bgCheck, seaWaterFraction
    type(struct_columnData)     :: columnFGE, columnSeaWaterFraction
    real(4), pointer            :: stateVectorFGE_ptr(:,:,:)
    logical                     :: checkMonth, llok
    integer                     :: lonIndex, latIndex, monthIndex, exceptMonthIndex

    call msg('ocebg_bgCheckSST', 'performing background check for the SST data...')

    ! get mpi topology
    call mmpi_setup_lonbands(hco%ni, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd)
    call mmpi_setup_latbands(hco%nj, latPerPE, latPerPEmax, myLatBeg, myLatEnd)
    
    ! Read the namelist
    if (.not. utl_isNamelistPresent('namOceanBGcheck','./flnml')) then
      if (mmpi_myid == 0) then
        call msg('ocebg_bgCheckSST', 'namOceanBGcheck is missing in the namelist.')
        call msg('ocebg_bgCheckSST', 'the default values will be taken.')
      end if
    else
      ! reading namelist variables
      nulnam = 0
      ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
      read(nulnam, nml = namOceanBGcheck, iostat = ierr)
      if (ierr /= 0) call utl_abort('ocebg_bgCheckSST: Error reading namelist')
      ierr = fclos(nulnam)
    end if
    call msg('ocebg_bgCheckSST', 'interpolation type: '//timeInterpType_nl)
    call msg('ocebg_bgCheckSST', 'number obs batches: '//numObsBatches)
    call msg('ocebg_bgCheckSST', 'check winds to detect tropical storms (TS): '//checkWinds)

    if (separateSelectCriteria) then
      call msg('ocebg_bgCheckSST', 'a separate selection criteria will be applied...')
      write(*,'(a,3f7.2)') 'ocebg_bgCheckSST: inland water selection criteria for INSITU data rejection set to: ', &
                           inlandWaterSelectCriteriaInsitu(:)
      write(*,'(a,3f7.2)') 'ocebg_bgCheckSST: inland water selection criteria for SATELLITE data rejection set to: ', &
                           inlandWaterSelectCriteriaSatData(:)
      write(*,'(a,3f7.2)') 'ocebg_bgCheckSST: sea water selection criteria for INSITU data rejection set to: ', &
                           seaWaterSelectCriteriaInsitu(:)
      write(*,'(a,3f7.2)') 'ocebg_bgCheckSST: sea water selection criteria for SATELLITE data rejection set to: ', &
                           seaWaterSelectCriteriaSatData(:)
      call msg('ocebg_bgCheckSST', 'sea water fraction threshold is set to: '//seaWaterThreshold)
      ! read sea water fraction
      call gsv_allocate(stateVectorSeaWaterFraction, 1, hco, columnTrlOnTrlLev%vco, dataKind_opt = 4, &
                        datestamp_opt = -1, mpi_local_opt = .true., &
                        varNames_opt = (/'VF'/), hInterpolateDegree_opt ='LINEAR')
      call gio_readFromFile(stateVectorSeaWaterFraction, './seaice_analysis', ' ','A', &
                            unitConversion_opt=.false., containsFullField_opt=.true.)
      ! Convert sea water fraction stateVector to column object
      call col_setVco(columnSeaWaterFraction, col_getVco(columnTrlOnTrlLev))
      call col_allocate(columnSeaWaterFraction, col_getNumCol(columnTrlOnTrlLev), varNames_opt = (/'VF'/))
      call s2c_nl(stateVectorSeaWaterFraction, obsData, columnSeaWaterFraction, hco, varName_opt = 'VF', &
                  timeInterpType = timeInterpType_nl, moveObsAtPole_opt = .true., &
                  numObsBatches_opt = numObsBatches, dealloc_opt = .true.)
    else
      write(*,'(a,3f7.2)') 'ocebg_bgCheckSST: global selection criteria for data rejection is set to: ', globalSelectCriteria(:)
    end if

    ! Read First Guess Error (FGE) and put it into stateVector
    call gsv_allocate(stateVectorFGE, 1, hco, columnTrlOnTrlLev%vco, dataKind_opt = 4, &
                      hInterpolateDegree_opt = 'NEAREST', &
                      datestamp_opt = -1, mpi_local_opt = .true., varNames_opt = (/'TM'/))
    call gio_readFromFile(stateVectorFGE, './bgstddev', 'STDDEV', 'X', &
                          unitConversion_opt=.false., containsFullField_opt=.true.)
    if (checkWinds) then
      call utl_tmg_start(123, '--checkWindsForSST') 
      call msg('ocebg_bgCheckSST', 'looking for tropical storms...')
      call msg('ocebg_bgCheckSST', 'number of days with available winds in the input winds file: '//ndaysWinds)
      call msg('ocebg_bgCheckSST', 'winds are provided every: '//timeStepWinds//' hours')
      call msg('ocebg_bgCheckSST', 'wind forecast lead time: '//windForecastLeadtime//' hours')
      do exceptMonthIndex = 1, nmonthsExceptionNH
        checkMonth = .False.
        loop_month: do monthIndex = 1, 12
          if (monthExceptionNH(exceptMonthIndex) == months(monthIndex)) then
            checkMonth = .True.
            exit loop_month
          end if
        end do loop_month
        if (.not. checkMonth) then
          call msg('ocebg_bgCheckSST', 'month should be one of these: '//months(:))
          call utl_abort('ocebg_bgCheckSST: unknown month '//monthExceptionNH(exceptMonthIndex))
        end if
      end do
      
      ! amplification error field state vector  
      call gsv_allocate(stateVectorAmplFactor, 1, hco, columnTrlOnTrlLev%vco, dataKind_opt = 4, &
                        hInterpolateDegree_opt = 'LINEAR', datestamp_opt = dateStamp, &
                        mpi_local_opt = .true., varNames_opt = (/'TM'/))
      call gsv_getField(stateVectorAmplFactor, stateVectorAmplFactor_ptr)
      stateVectorAmplFactor_ptr(myLonBeg:myLonEnd,myLatBeg:myLatEnd,1) = 1.0d0

      call ocebg_getFGEamplification(stateVectorAmplFactor, dateStamp, hco)
      
      ! FGE state vector
      call gsv_getField(stateVectorFGE, stateVectorFGE_ptr)

      ! Apply tropical storm correction to the FGE field:
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          stateVectorFGE_ptr(lonIndex, latIndex, 1) = stateVectorFGE_ptr(lonIndex, latIndex, 1) * &
                                                      stateVectorAmplFactor_ptr(lonIndex, latIndex, 1)
        end do
      end do      
      call gsv_deallocate(stateVectorAmplFactor)
      call utl_tmg_stop(123)
    end if

    ! Convert FGE stateVector to column object
    call col_setVco(columnFGE, col_getVco(columnTrlOnTrlLev))
    call col_allocate(columnFGE, col_getNumCol(columnTrlOnTrlLev), varNames_opt = (/'TM'/))
    call s2c_nl(stateVectorFGE, obsData, columnFGE, hco, varName_opt = 'TM', &
                timeInterpType = timeInterpType_nl, moveObsAtPole_opt = .true., &
                numObsBatches_opt = numObsBatches, dealloc_opt = .true.)

    numberObs = 0
    numberObsRejected = 0
    numberObsInsitu = 0
    numberObsInsituRejected = 0

    do headerIndex = 1, obs_numheader(obsData)

      bodyIndex = obs_headElem_i(obsData, obs_rln, headerIndex)
      obsVarno  = obs_bodyElem_i(obsData, obs_vnm, bodyIndex)
      llok = (obs_bodyElem_i(obsData, obs_ass, bodyIndex) == obs_assimilated)
      if (llok) then
        if (obsVarno == bufr_sst) then

	  FGE = col_getElem(columnFGE, 1, headerIndex, 'TM')
	  OmP = obs_bodyElem_r(obsData, OBS_OMP , bodyIndex)
          OER = obs_bodyElem_r(obsData, OBS_OER , bodyIndex)
	  codeType = obs_headElem_i(obsData, obs_ity, headerIndex)

	  if (FGE /= MPC_missingValue_R8 .and. OmP /= MPC_missingValue_R8) then

	    numberObs = numberObs + 1
	    if (codeType /= codtyp_get_codtyp('satob')) numberObsInsitu = numberObsInsitu + 1
	    call obs_bodySet_r(obsData, OBS_HPHT, bodyIndex, FGE)
	    bgCheck = (OmP)**2 / (FGE**2 + OER**2)

            if (separateSelectCriteria) then
              seaWaterFraction = col_getElem(columnSeaWaterFraction, 1, headerIndex, 'VF')
              if (seaWaterFraction <= seaWaterThreshold) then
                if(codeType == codtyp_get_codtyp('satob')) then
                  obsFlag = ocebg_setFlag(obsVarno, bgCheck, inlandWaterSelectCriteriaSatData)
                else
                  obsFlag = ocebg_setFlag(obsVarno, bgCheck, inlandWaterSelectCriteriaInsitu)
                end if
              else
                if(codeType == codtyp_get_codtyp('satob')) then
                  obsFlag = ocebg_setFlag(obsVarno, bgCheck, seaWaterSelectCriteriaSatData)
                else
                  obsFlag = ocebg_setFlag(obsVarno, bgCheck, seaWaterSelectCriteriaInsitu)
                end if 
              end if
            else
              obsFlag = ocebg_setFlag(obsVarno, bgCheck, globalSelectCriteria)
            end if
	
            if (obsFlag >= 2) then
              numberObsRejected = numberObsRejected + 1
	      if (codeType /= codtyp_get_codtyp('satob')) numberObsInsituRejected = numberObsInsituRejected + 1
	      write(*,'(i10,a,i5,4f10.4,i5)') numberObsRejected, &
                                              ', sensor, codtype, lon, lat, obs.value, OmP, flag: '&
                                              //obs_elem_c(obsData, 'STID' , headerIndex)//' data: ', codeType, &
                                              obs_headElem_r(obsData, obs_lon, headerIndex) * MPC_DEGREES_PER_RADIAN_R8,&
                                              obs_headElem_r(obsData, obs_lat, headerIndex) * MPC_DEGREES_PER_RADIAN_R8,&
                                              obs_bodyElem_r(obsData, obs_var, bodyIndex), OmP, obsFlag 
            end if

	    ! update background check flags based on bgCheck
            ! (element flags + global header flags)
	    if (obsFlag == 1) then
              call obs_bodySet_i(obsData, obs_flg, bodyIndex  , ibset(obs_bodyElem_i(obsData, obs_flg, bodyIndex)  , 13))
            else if (obsFlag == 2) then
              call obs_bodySet_i(obsData, obs_flg, bodyIndex  , ibset(obs_bodyElem_i(obsData, obs_flg, bodyIndex)  , 14))
              call obs_bodySet_i(obsData, obs_flg, bodyIndex  , ibset(obs_bodyElem_i(obsData, obs_flg, bodyIndex)  , 16))
              call obs_bodySet_i(obsData, obs_flg, bodyIndex  , ibset(obs_bodyElem_i(obsData, obs_flg, bodyIndex)  , 09))
              call obs_headSet_i(obsData, obs_st1, headerIndex, ibset(obs_headElem_i(obsData, obs_st1, headerIndex), 06))
            else if (obsFlag == 3) then
              call obs_bodySet_i(obsData, obs_flg, bodyIndex  , ibset(obs_bodyElem_i(obsData, obs_flg, bodyIndex)  , 15))
              call obs_bodySet_i(obsData, obs_flg, bodyIndex  , ibset(obs_bodyElem_i(obsData, obs_flg, bodyIndex)  , 16))
              call obs_bodySet_i(obsData, obs_flg, bodyIndex  , ibset(obs_bodyElem_i(obsData, obs_flg, bodyIndex)  , 09))
              call obs_headSet_i(obsData, obs_st1, headerIndex, ibset(obs_headElem_i(obsData, obs_st1, headerIndex), 06))
            end if

          end if
	end if
      end if

    end do

    if (numberObs > 0) then
      write(*,*)
      call msg('ocebg_bgCheckSST', 'background check of SST (TM) data is computed')
      write(*,*) '***************************************************************************************'
      write(*,'(a, i7,a,i7,a)') 'ocebg_bgCheckSST: total ', numberObsRejected, ' observations out of (ALL) ', numberObs,' rejected'
      write(*,'(a, i7,a,i7,a)') 'ocebg_bgCheckSST: where ', numberObsInsituRejected, ' insitu observations out of ', &
                                numberObsInsitu,' insitu obs. rejected'
      write(*,*) '***************************************************************************************'
      write(*,*)
    end if

    call gsv_deallocate(stateVectorFGE)
    call col_deallocate(columnFGE)
    if (separateSelectCriteria) then
      call col_deallocate(columnSeaWaterFraction)
      call gsv_deallocate(stateVectorSeaWaterFraction)
    end if      

  end subroutine ocebg_bgCheckSST

  !----------------------------------------------------------------------------------------
  ! ocebg_bgCheckSeaIce
  !----------------------------------------------------------------------------------------
  subroutine ocebg_bgCheckSeaIce(obsData)
    !
    ! :Purpose: Compute sea ice data background check
    !

    implicit none

    ! Arguments:
    type(struct_obs), intent(inout) :: obsData           ! obsSpaceData object

    ! Locals:
    integer            :: nulnam, ierr
    integer            :: headerIndex, bodyIndex, stationIndex, bodyCount
    integer            :: obsChid, obsDate, obsTime, obsVarno, obsFlag
    integer            :: obsDateStamp
    integer, parameter :: maxSwath = 10, maxPerSwath = 100000
    integer            :: numberObs(maxSwath), bodyIndexList(maxPerSwath,maxSwath)
    integer            :: minDateStamp(maxSwath), maxDateStamp(maxSwath), swathID(maxSwath)
    real               :: rmsDiff(maxSwath)
    integer            :: swathIndex, nSwath
    real(8)            :: OmP, deltaHours
    integer            :: numberObsRejected
    character(len=12)  :: cstnid
    integer, external  :: newdate
    integer            :: imode, prntdate, prnttime

    ! Namelist variables: (local)
    integer           :: numStation = 11
    character(len=12) :: idStation(11) = 'null'
    real              :: OmpRmsdThresh(11) = 0.0
    namelist /namIceBGcheck/ numStation, idStation, OmpRmsdThresh

    call msg('ocebg_bgCheckSeaIce', 'performing background check for the SeaIce data...')
 
    ! Read the namelist
    if (.not. utl_isNamelistPresent('namIceBGcheck','./flnml')) then
      if (mpi_myid == 0) then
        call msg('ocebg_bgCheckSeaIce', 'namIceBGcheck is missing in the namelist.')
        call msg('ocebg_bgCheckSeaIce', 'the default values will be taken.')
      end if
    else
      ! reading namelist variables
      nulnam = 0
      ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
      read(nulnam, nml = namIceBGcheck, iostat = ierr)
      if (ierr /= 0) call utl_abort('ocebg_bgCheckSeaIce: Error reading namelist')
      ierr = fclos(nulnam)
    end if
    write(*, nml = namIceBGcheck)

    STATION: do stationIndex = 1, numStation

        call msg('ocebg_bgCheckSeaIce', 'doing idStation: '//idStation(stationIndex))

        if (idStation(stationIndex) /= 'null') then

        numberObs(:) = 0
        rmsDiff(:) = 0.0
        minDateStamp(:) = tim_getDatestamp()
        maxDateStamp(:) = tim_getDatestamp()
        nSwath = 0

        HEADER: do headerIndex = 1, obs_numheader(obsData)

          cstnid = obs_elem_c(obsData, 'STID' , headerIndex )
          bodyIndex = obs_headElem_i(obsData, obs_rln, headerIndex)

!!$          obsVarno = obs_bodyElem_i(obsData, obs_vnm, bodyIndex)
!!$
!!$          if (obsVarno == BUFR_ICEC) then

          obsFlag = obs_bodyElem_i(obsData, OBS_FLG, bodyIndex)

          if ( (cstnid == idStation(stationIndex)) .and. &
               (obsFlag == 0 .or. cstnid == 'CIS_REGIONAL') ) then

            obsChid = obs_headElem_i(obsData, obs_chid, headerIndex)
            obsDate = obs_headElem_i(obsData, OBS_DAT, headerIndex)
            obsTime = obs_headElem_i(obsData, OBS_ETM, headerIndex)

            imode = 3 ! printable to stamp
            ierr = newdate(obsDateStamp, obsDate, obsTime*10000, imode)

            do swathIndex=1, nSwath

              if(obsDateStamp <= maxDateStamp(swathIndex) .and. &
                 obsDateStamp >= minDateStamp(swathIndex) .and. &
                 obsChid == swathID(swathIndex)) exit

              call difdatr(obsDateStamp, minDateStamp(swathIndex), deltaHours)
            
              if(abs(deltaHours) <= 0.5d0 .and. obsChid == swathID(swathIndex)) exit

            end do

            nSwath = max(nSwath, swathIndex)

            if(swathIndex <= maxSwath) then

              if(numberObs(swathIndex) == 0) then
                minDateStamp(swathIndex) = obsDateStamp
                maxDateStamp(swathIndex) = obsDateStamp
              else
                minDateStamp(swathIndex) = min(minDateStamp(swathIndex), obsDateStamp)
                maxDateStamp(swathIndex) = max(maxDateStamp(swathIndex), obsDateStamp)
              end if

              swathID(swathIndex) = obsChid

              OmP = obs_bodyElem_r(obsData, OBS_OMP , bodyIndex)

              if (OmP /= MPC_missingValue_R8) then

                numberObs(swathIndex) = numberObs(swathIndex) + 1
                rmsDiff(swathIndex) = rmsDiff(swathIndex) + (OmP)**2
                if (numberObs(swathIndex) <= maxPerSwath) then
                  bodyIndexList(numberObs(swathIndex), swathIndex) = bodyIndex
                else
                  call utl_abort('ocebg_bgCheckSeaIce: Increase maxPerSwath')
                end if

              end if

            else

              call utl_abort('ocebg_bgCheckSeaIce: Increase maxSwath')

            end if

          end if

!!$            end if

        end do HEADER

        do swathIndex = 1, maxSwath

          if (numberObs(swathIndex) > 0) then

            rmsDiff(swathIndex) = sqrt(rmsDiff(swathIndex)/numberObs(swathIndex))

            call msg('ocebg_bgCheckSeaIce', 'swathIndex = '//swathIndex)
            call msg('ocebg_bgCheckSeaIce', 'Timespan:')
            imode = -3 ! stamp to printable
            ierr = newdate(minDateStamp(swathIndex), prntdate, prnttime, imode)
            call msg('ocebg_bgCheckSeaIce', 'From: '//prntdate//prnttime/100)
            ierr = newdate(maxDateStamp(swathIndex), prntdate, prnttime, imode)
            call msg('ocebg_bgCheckSeaIce', 'To  : '//prntdate//prnttime/100)

            write(*,'(a,f10.4)') 'ocebg_bgCheckSeaIce: RMSD = ', rmsDiff(swathIndex)
            write(*,'(a,i10)') 'ocebg_bgCheckSeaIce: nobs = ', numberObs(swathIndex)

            if (rmsDiff(swathIndex) > OmpRmsdThresh(stationIndex)) then

              numberObsRejected = 0
              do bodyCount = 1, numberObs(swathIndex)
                numberObsRejected = numberObsRejected + 1
                ! update background check flag
                call obs_bodySet_i(obsData, obs_flg, bodyIndexList(bodyCount,swathIndex), ibset(obs_bodyElem_i(obsData, obs_flg, bodyIndexList(bodyCount,swathIndex)), 10))
              end do

	      write(*,'(a,i7,a)')'ocebg_bgCheckSeaIce: ********** reject: ', numberObsRejected, &
                                    ' '//idStation(stationIndex)//' data'

            end if

          end if

        end do

      end if

    end do STATION

  end subroutine ocebg_bgCheckSeaIce

  !--------------------------------------------------------------------------
  ! ocebg_setFlag
  !--------------------------------------------------------------------------
  function ocebg_setFlag(obsVarno, bgCheck, selectCriteria) result(obsFlag)
    !
    ! :Purpose: Set background-check flags according to values set in a table.
    !           Original values in table come from ECMWF.
    !

    implicit none

    integer             :: obsFlag  ! obs flag

    ! Arguments:
    integer, intent(in) :: obsVarno          ! obsVarno, Universal Field-Identity Numbers defined in bufr_mod
    real(8), intent(in) :: bgCheck           ! normalized background departure
    real(8), intent(in) :: selectCriteria(:) ! selection criteria for three levels 

    obsFlag = 0
 
    if (obsVarno == bufr_sst) then
      if (bgCheck >= selectCriteria(1) .and. bgCheck < selectCriteria(2)) then
        obsFlag = 1
      else if (bgCheck >= selectCriteria(2) .and. bgCheck < selectCriteria(3)) then
        obsFlag = 2
      else if (bgCheck >= selectCriteria(3)) then
        obsFlag = 3
      end if
    end if

  end function ocebg_setFlag

  !--------------------------------------------------------------------------
  ! ocebg_getFGEamplification
  !--------------------------------------------------------------------------
  subroutine ocebg_getFGEamplification(stateVectorAmplFactor, dateStamp, hco)
    !
    ! :Purpose: Read wind speed fields for the last four days.
    !           In the operations: 
    !           The background error used during the background check is then
    !           amplified in those regions by a factor that varies from 1,
    !           where the maximum wind speed is 21m/s or less, to 12,
    !           where the maximum wind speed is 24m/s or more.
    !           The factor is then filtered to produce a smoothly varying
    !           field. This amplified background error is used only to
    !           perform the background check. 

    implicit none

    ! Arguments
    type(struct_gsv), intent(inout)       :: stateVectorAmplFactor ! state vector to save amplification factor
    integer         , intent(in)          :: dateStamp             ! date stamp
    type(struct_hco), intent(in), pointer :: hco                   ! horizontal trl grid

    ! locals
    type(struct_gsv)          :: stateVector         ! state vector for surface winds
    type(struct_vco), pointer :: vco_winds           ! vertical grid structure for winds
    real(4)         , pointer :: uu_ptr4d(:,:,:,:)   ! pointer to get UU wind component
    real(4)         , pointer :: vv_ptr4d(:,:,:,:)   ! pointer to get VV wind component
    integer                   :: dataStampList(ndaysWinds * 24 / timeStepWinds) ! datastamp list for wind fields
    real(4)                   :: windSpeed
    integer                   :: hour, day, monthNumber
    integer                   :: yyyy, ndays, timeStepIndex
    real(8)                   :: deltaT, lat
    integer                   :: lonIndex, latIndex, monthIndex
    real(4)         , pointer :: stateVectorAmplFactor_ptr(:,:,:)
    real(8)                   :: amplFactor

    nullify(vco_winds)
    
    ! looking for the earliest valid time:
    deltaT = - real(ndaysWinds * 24 - windForecastLeadtime)
    call incdatr(dataStampList(1), dateStamp, deltaT)

    ! scanning the datastamps up from the earliest to the last valid time
    deltaT = real(timeStepWinds)
    do timeStepIndex = 2, ndaysWinds * 24 / timeStepWinds
      call incdatr(dataStampList(timeStepIndex), dataStampList(timeStepIndex - 1), deltaT)
    end do

    call vco_SetupFromFile(vco_winds,'./winds')
    call gsv_allocate(stateVector, ndaysWinds * 24 / timeStepWinds, hco, vco_winds, &
                      dateStampList_opt = dataStampList, dataKind_opt = 4, &
                      varNames_opt=(/'UU','VV'/), mpi_local_opt=.true., &
                      hInterpolateDegree_opt='LINEAR')
      
    call msg('ocebg_getFGEamplification', 'reading wind speed fields...')
    do timeStepIndex = 1, ndaysWinds * 24 / timeStepWinds
      call tim_dateStampToYYYYMMDDHH(dataStampList(timeStepIndex), hour, day, monthNumber, &
                                     ndays, yyyy, verbose_opt = .False.)
      call msg('ocebg_getFGEamplification', ' '//timeStepIndex//dataStampList(timeStepIndex)//yyyy//monthNumber//day//hour)
      call gio_readFromFile(stateVector, './winds', ' ', ' ', stepIndex_opt = timeStepIndex, &
                            containsFullField_opt=.true.)
    end do

    call gsv_getField(stateVector, uu_ptr4d, 'UU')
    call gsv_getField(stateVector, vv_ptr4d, 'VV')
    call gsv_getField(stateVectorAmplFactor, stateVectorAmplFactor_ptr)

    write(*,*)
    call msg('ocebg_getFGEamplification', 'detecting tropical storms (TS)...')

    time_loop: do timeStepIndex = 1, ndaysWinds * 24 / timeStepWinds

      call tim_dateStampToYYYYMMDDHH(dataStampList(timeStepIndex), hour, day, monthNumber, &
                                     ndays, yyyy, verbose_opt = .False.)
      do monthIndex = 1, nmonthsExceptionNH
        if (months(monthNumber) == monthExceptionNH(monthIndex)) then
          maxLatNH = maxLatExceptionNH
          call msg('ocebg_getFGEamplification', 'TS is allowed to penetrate up to '//maxLatNH &
                     //' degrees North; current month: '//months(monthNumber))
        end if
      end do

      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd

          lat = real(hco%lat2d_4(lonIndex, latIndex), 8) * MPC_DEGREES_PER_RADIAN_R8
          if ((lat > minLatNH .and. lat < maxLatNH) .or. (lat > minLatSH .and. lat < maxLatSH)) then
            windSpeed = sqrt(uu_ptr4d(lonIndex, latIndex, 1, timeStepIndex) * &
                             uu_ptr4d(lonIndex, latIndex, 1, timeStepIndex) + &
                             vv_ptr4d(lonIndex, latIndex, 1, timeStepIndex) * &
                             vv_ptr4d(lonIndex, latIndex, 1, timeStepIndex))
            amplFactor = max(1.0d0, min(4.0d0 * max(0.0d0,(windSpeed - 20.75d0)), 12.0d0))
            ! stateVectorAmplFactor_ptr may already be assigned with a value > 1.0 in a previous state 
            if(amplFactor > stateVectorAmplFactor_ptr(lonIndex, latIndex, 1)) then
              stateVectorAmplFactor_ptr(lonIndex, latIndex, 1) = amplFactor
            end if
          endif

        end do
      end do

    end do time_loop

    call gsv_deallocate(stateVector)

    call gio_writeToFile(stateVectorAmplFactor, './amplification', 'ORIG')      
    call gsv_smoothHorizontal(stateVectorAmplFactor, smoothLenghtScale)
    call gio_writeToFile(stateVectorAmplFactor, './amplification', 'SMOOTH')        

  end subroutine ocebg_getFGEamplification

end module bgckOcean_mod
