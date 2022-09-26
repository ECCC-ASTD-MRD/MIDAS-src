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

module ensPostProcess_mod
  ! MODULE ensPostProcess_mod (prefix='epp' category='1. High-level functionality')
  !
  ! :Purpose: Various routines that are used to modify or process
  !           ensembles, usually produced by the LETKF.
  !
  use midasMpi_mod
  use utilities_mod
  use mathPhysConstants_mod
  use timeCoord_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use ensembleStateVector_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use interpolation_mod
  use randomNumber_mod
  use controlVector_mod
  use gridVariableTransforms_mod
  use bMatrix_mod
  use humidityLimits_mod
  use localizationFunction_mod
  use varNameList_mod
  use fileNames_mod
  use clib_interfaces_mod
  use physicsFunctions_mod
  implicit none
  save
  private

  ! public procedures
  public :: epp_postProcess

contains

  !----------------------------------------------------------------------
  ! epp_postProcess
  !----------------------------------------------------------------------
  subroutine epp_postProcess(ensembleTrl, ensembleAnl, &
                             stateVectorHeightSfc, stateVectorCtrlTrl, &
                             writeTrlEnsemble, outputOnlyEnsMean_opt)
    !
    !:Purpose:  Perform numerous post-processing steps to the ensemble
    !           produced by the LETKF algorithm.
    !
    implicit none

    ! Arguments
    type(struct_ens), intent(inout) :: ensembleTrl
    type(struct_ens), intent(inout) :: ensembleAnl
    type(struct_gsv), intent(in)    :: stateVectorHeightSfc
    type(struct_gsv), intent(inout) :: stateVectorCtrlTrl
    logical, intent(in)             :: writeTrlEnsemble
    logical, optional, intent(in)   :: outputOnlyEnsMean_opt

    ! Locals
    integer                   :: ierr, nEns, dateStamp, datePrint, timePrint, imode, randomSeedRandomPert
    integer                   :: stepIndex, middleStepIndex, nulnam
    integer, allocatable      :: dateStampListInc(:)
    type(struct_hco), pointer :: hco_ens
    type(struct_vco), pointer :: vco_ens
    type(struct_gsv)          :: stateVectorMeanAnl, stateVectorMeanAnlRaw, stateVectorMeanTrl
    type(struct_gsv)          :: stateVectorMeanInc
    type(struct_gsv)          :: stateVectorAnalIncMask
    type(struct_gsv)          :: stateVectorStdDevAnl, stateVectorStdDevAnlRaw, stateVectorStdDevAnlPert, stateVectorStdDevTrl
    type(struct_gsv)          :: stateVectorMeanIncSubSample
    type(struct_gsv)          :: stateVectorMeanAnlSubSample
    type(struct_gsv)          :: stateVectorMeanAnlSfcPres
    type(struct_gsv)          :: stateVectorMeanAnlSfcPresMpiGlb
    type(struct_ens)          :: ensembleTrlSubSample
    type(struct_ens)          :: ensembleAnlSubSample
    type(struct_ens)          :: ensembleAnlSubSampleUnPert
    character(len=12)         :: etiket
    type(struct_ens)          :: ensembleAnlInc
    type(struct_ens)          :: ensembleAnlIncSubSample
    character(len=256)        :: outFileName
    character(len=4), pointer :: varNames(:)
    character(len=12)         :: hInterpolationDegree = 'LINEAR'
    integer, external         :: fnom, fclos, newdate
    logical                   :: outputOnlyEnsMean

    ! Namelist variables
    integer  :: randomSeed           ! seed used for random perturbation additive inflation
    logical  :: includeYearInSeed    ! switch for choosing to include year in default random seed
    logical  :: writeSubSample       ! write sub-sample members for initializing medium-range fcsts
    logical  :: writeSubSampleUnPert ! write unperturbed sub-sample members for medium-range fcsts
    real(8)  :: alphaRTPS            ! RTPS coefficient (between 0 and 1; 0 means no relaxation)
    real(8)  :: alphaRTPP            ! RTPP coefficient (between 0 and 1; 0 means no relaxation)
    real(8)  :: alphaRandomPert      ! Random perturbation additive inflation coeff (0->1)
    real(8)  :: alphaRandomPertSubSample ! Random pert. additive inflation coeff for medium-range fcsts
    logical  :: huLimitsBeforeRecenter   ! Choose to apply humidity limits before recentering
    logical  :: imposeSaturationLimit  ! switch for choosing to impose saturation limit of humidity
    logical  :: imposeRttovHuLimits    ! switch for choosing to impose the RTTOV limits on humidity
    real(8)  :: weightRecenter         ! weight applied to recentering increment
    real(8)  :: weightRecenterLand     ! weight applied to recentering increment for land variables
    integer  :: numMembersToRecenter   ! number of members that get recentered on supplied analysis
    logical  :: useOptionTableRecenter ! use values in the optiontable file
    character(len=8)  :: etiket_anl ! etiket for ensemble output files (member number will be appended)
    character(len=8)  :: etiket_inc ! etiket for ensemble output files (member number will be appended)
    character(len=8)  :: etiket_trl ! etiket for ensemble output files (member number will be appended)
    character(len=12) :: etiket_anlmean ! etiket for mean of analyses and mean of increments files
    character(len=12) :: etiket_anlrms  ! etiket for rms of analyses files
    character(len=12) :: etiket_anlmean_raw ! etiket for mean of raw analyses
    character(len=12) :: etiket_anlrms_raw  ! etiket for rms of raw analyses
    character(len=12) :: etiket_anlmeanpert ! etiket for mean of perturbed analyses
    character(len=12) :: etiket_anlrmspert  ! etiket for rms of perturbed analyses
    character(len=12) :: etiket_trlmean     ! etiket for mean of trials
    character(len=12) :: etiket_trlrms      ! etiket for rms of trials
    integer  :: numBits ! number of bits when writing ensemble mean and spread
    logical  :: useAnalIncMask        ! mask out the increment on the pilot zone
    logical  :: writeRawAnalStats     ! write mean and standard deviation of the raw analysis ensemble
    logical  :: useMemberAsHuRefState ! use each member as reference state for variable transforms 

    NAMELIST /namEnsPostProcModule/randomSeed, includeYearInSeed, writeSubSample, writeSubSampleUnPert,  &
                                   alphaRTPS, alphaRTPP, alphaRandomPert, alphaRandomPertSubSample,      &
                                   huLimitsBeforeRecenter, imposeSaturationLimit, imposeRttovHuLimits,   &
                                   weightRecenter, weightRecenterLand, numMembersToRecenter,  &
                                   useOptionTableRecenter,  &
                                   etiket_anl, etiket_inc, etiket_trl, etiket_anlmean, etiket_anlrms,    &
                                   etiket_anlmeanpert, etiket_anlrmspert,  &
                                   etiket_anlmean_raw, etiket_anlrms_raw,  &
                                   etiket_trlmean, etiket_trlrms, numBits, useAnalIncMask,  &
                                   writeRawAnalStats, useMemberAsHuRefState

    if (present(outputOnlyEnsMean_opt)) then
      outputOnlyEnsMean = outputOnlyEnsMean_opt
    else
      outputOnlyEnsMean = .false.
    end if

    !- Extract the grid definitions and ensemble size
    if (ens_isAllocated(ensembleTrl)) then
      hco_ens => ens_getHco(ensembleTrl)
      vco_ens => ens_getVco(ensembleTrl)
      nEns = ens_getNumMembers(ensembleTrl)
    else
      hco_ens => ens_getHco(ensembleAnl)
      vco_ens => ens_getVco(ensembleAnl)
      nEns = ens_getNumMembers(ensembleAnl)
    end if

    !- Setting default namelist variable values
    randomSeed            =  -999
    includeYearInSeed     = .false.
    writeSubSample        = .false.
    writeSubSampleUnPert  = .false.
    alphaRTPS             =  0.0D0
    alphaRTPP             =  0.0D0
    alphaRandomPert       =  0.0D0
    alphaRandomPertSubSample =  -1.0D0
    huLimitsBeforeRecenter = .true.
    imposeSaturationLimit = .false.
    imposeRttovHuLimits   = .false.
    weightRecenter        = 0.0D0  ! means no recentering applied
    weightRecenterLand    = -1.0D0 ! means same recentering for land as other variables
    numMembersToRecenter  = -1     ! means all members recentered by default
    useOptionTableRecenter = .false.
    ! For the next 3 variables, the member number will be appended to this string
    ! for files '${trialdate}_006_${member}'
    etiket_trl = 'E27_0_0P'
    ! for files '${analdate}_000_${member}', 'subspace/${analdate}_000_${member}' and 'subspace_unpert/${analdate}_000_${member}'
    ! the member=0000 will contain the mean analysis
    etiket_anl = 'E27_0_0P'
    ! for files '${analdate}_000_inc_${member}' and 'subspace/${analdate}_000_inc_${member}'
    ! the member=0000 will contain the mean increment
    etiket_inc = 'E27_0_0P'
    !
    etiket_anlmean = 'E27_0_0PAVG' ! for file '${analdate}_000_analmean'
    etiket_anlrms = 'E27_0_0PRMS'  ! for file '${analdate}_000_analrms'
    etiket_anlmeanpert = 'E27_0_0PAVGP'  ! for file '${analdate}_000_analpertmean'
    etiket_anlrmspert = 'E27_0_0PRMSP'  ! for file '${analdate}_000_analpertrms'
    etiket_anlmean_raw = 'E27_0_0PAVGR'  ! for file '${analdate}_000_analmean_raw'
    etiket_anlrms_raw = 'E27_0_0PRMSR'  ! for file '${analdate}_000_analrms_raw'
    etiket_trlmean = 'E27_0_0PAVG' ! for file '${trialdate}_006_trialmean'
    etiket_trlrms = 'E27_0_0PRMS' ! for file '${trialdate}_006_trialrms'
    !
    numBits = 16
    useAnalIncMask = .false.
    writeRawAnalStats = .false.
    useMemberAsHuRefState = .false.

    !- Read the namelist
    nulnam = 0
    ierr = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0)
    read(nulnam, nml=namEnsPostProcModule, iostat=ierr)
    if ( ierr /= 0) call utl_abort('epp_postProc: Error reading namelist')
    if ( mmpi_myid == 0 ) write(*,nml=namEnsPostProcModule)
    ierr = fclos(nulnam)

    if (alphaRTPS < 0.0D0) alphaRTPS = 0.0D0
    if (alphaRTPP < 0.0D0) alphaRTPP = 0.0D0
    if (alphaRandomPert < 0.0D0) alphaRandomPert = 0.0D0
    if (alphaRandomPertSubSample < 0.0D0) alphaRandomPertSubSample = 0.0D0
    if (numMembersToRecenter == -1) numMembersToRecenter = nEns ! default behaviour

    if (writeSubSample) then
      if (.not.(ens_isAllocated(ensembleTrl).and.ens_isAllocated(ensembleAnl))) then
        call utl_abort('epp_postProc: subSample can only be produced if both Anl and Trl ensembles available')
      end if
    end if

    if (writeSubSampleUnPert) then
      if (.not.ens_isAllocated(ensembleAnl)) then
        call utl_abort('epp_postProc: subSampleUnPert can only be produced if Anl ensemble available')
      end if
    end if

    if (writeRawAnalStats) then
      if (.not.ens_isAllocated(ensembleAnl)) then
        call utl_abort('epp_postProc: RawAnalStats can only be produced if Anl ensemble available')
      end if
    end if

    if (ens_isAllocated(ensembleTrl)) then
      !- Allocate and compute ensemble mean Trl
      call gsv_allocate( stateVectorMeanTrl, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                         mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                         dataKind_opt=4, allocHeightSfc_opt=.true., &
                         hInterpolateDegree_opt = hInterpolationDegree, &
                         allocHeight_opt=.false., allocPressure_opt=.false. )
      call gsv_zero(stateVectorMeanTrl)
      call ens_computeMean(ensembleTrl)
      call ens_copyEnsMean(ensembleTrl, stateVectorMeanTrl)
      
      !- Allocate and compute ensemble spread stddev Trl
      call gsv_allocate( stateVectorStdDevTrl, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                         mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                         hInterpolateDegree_opt = hInterpolationDegree, &
                         dataKind_opt=4, allocHeight_opt=.false., allocPressure_opt=.false. )
      call gsv_zero(stateVectorStdDevTrl)
      call ens_computeStdDev(ensembleTrl)
      call ens_copyEnsStdDev(ensembleTrl, stateVectorStdDevTrl)
    end if

    if (ens_isAllocated(ensembleAnl)) then
      !- Allocate and compute ensemble mean Anl
      call gsv_allocate( stateVectorMeanAnl, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                         mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                         dataKind_opt=4, allocHeightSfc_opt=.true., &
                         hInterpolateDegree_opt = hInterpolationDegree, &
                         allocHeight_opt=.false., allocPressure_opt=.false. )
      call gsv_zero(stateVectorMeanAnl)
      call ens_computeMean(ensembleAnl)
      call ens_copyEnsMean(ensembleAnl, stateVectorMeanAnl)

      !- Allocate and compute ensemble spread stddev Anl
      call gsv_allocate( stateVectorStdDevAnl, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                         mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                         hInterpolateDegree_opt = hInterpolationDegree, &
                         dataKind_opt=4, allocHeight_opt=.false., allocPressure_opt=.false. )
      call gsv_zero(stateVectorStdDevAnl)
      call ens_computeStdDev(ensembleAnl)
      call ens_copyEnsStdDev(ensembleAnl, stateVectorStdDevAnl)

      !- Keep the mean and standard deviation of the raw analysis for outputs, if requested
      if (writeRawAnalStats) then
        call gsv_allocate( stateVectorMeanAnlRaw, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                           mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                           dataKind_opt=4, allocHeightSfc_opt=.true., &
                           hInterpolateDegree_opt = hInterpolationDegree, &
                           allocHeight_opt=.false., allocPressure_opt=.false. )
        call gsv_allocate( stateVectorStdDevAnlRaw, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                           mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                           hInterpolateDegree_opt = hInterpolationDegree, &
                           dataKind_opt=4, allocHeight_opt=.false., allocPressure_opt=.false. )
        call gsv_copy(stateVectorMeanAnl, stateVectorMeanAnlRaw)
        call gsv_copy(stateVectorStdDevAnl, stateVectorStdDevAnlRaw)
      end if

      if (ens_isAllocated(ensembleTrl)) then
        !- Apply RTPP, if requested
        if (alphaRTPP > 0.0D0) then
          call epp_RTPP(ensembleAnl, ensembleTrl, stateVectorMeanAnl, &
                        stateVectorMeanTrl, alphaRTPP)
          ! recompute the analysis spread stddev
          call ens_computeStdDev(ensembleAnl)
          call ens_copyEnsStdDev(ensembleAnl, stateVectorStdDevAnl)
        end if

        !- Apply RTPS, if requested
        if (alphaRTPS > 0.0D0) then
          call epp_RTPS(ensembleAnl, stateVectorStdDevAnl, stateVectorStdDevTrl, &
                        stateVectorMeanAnl, alphaRTPS)
          ! recompute the analysis spread stddev
          call ens_computeStdDev(ensembleAnl)
          call ens_copyEnsStdDev(ensembleAnl, stateVectorStdDevAnl)
        end if
      end if

      !- Impose limits on humidity *before* recentering, if requested
      if (huLimitsBeforeRecenter) then
        if (imposeSaturationLimit .or. imposeRttovHuLimits) then
          if (mmpi_myid == 0) write(*,*) ''
          if (mmpi_myid == 0) write(*,*) 'epp_postProcess: limits will be imposed on the humidity of analysis ensemble'
          if (mmpi_myid == 0 .and. imposeSaturationLimit ) write(*,*) '              -> Saturation Limit'
          if (mmpi_myid == 0 .and. imposeRttovHuLimits   ) write(*,*) '              -> Rttov Limit'
          if ( imposeSaturationLimit ) call qlim_saturationLimit(ensembleAnl)
          if ( imposeRttovHuLimits   ) call qlim_rttovLimit     (ensembleAnl)
          ! And recompute analysis mean
          call ens_computeMean(ensembleAnl)
          call ens_copyEnsMean(ensembleAnl, stateVectorMeanAnl)
        end if
      end if

      !- Recenter analysis ensemble on supplied analysis
      if (weightRecenter > 0.0D0 .or. useOptionTableRecenter) then
        write(*,*) 'epp_postProcess: Recenter analyses on supplied analysis'
        call epp_hybridRecentering(ensembleAnl, weightRecenter, weightRecenterLand, &
                                   useOptionTableRecenter, numMembersToRecenter)
        ! And recompute analysis mean
        call ens_computeMean(ensembleAnl)
        call ens_copyEnsMean(ensembleAnl, stateVectorMeanAnl)
        ! And recompute the analysis spread stddev
        call ens_computeStdDev(ensembleAnl)
        call ens_copyEnsStdDev(ensembleAnl, stateVectorStdDevAnl)
      end if

      !- Impose limits on humidity *after* recentering, if requested
      if (.not.huLimitsBeforeRecenter) then
        if (imposeSaturationLimit .or. imposeRttovHuLimits) then
          if (mmpi_myid == 0) write(*,*) ''
          if (mmpi_myid == 0) write(*,*) 'epp_postProcess: limits will be imposed on the humidity of analysis ensemble'
          if (mmpi_myid == 0 .and. imposeSaturationLimit ) write(*,*) '              -> Saturation Limit'
          if (mmpi_myid == 0 .and. imposeRttovHuLimits   ) write(*,*) '              -> Rttov Limit'
          if ( imposeSaturationLimit ) call qlim_saturationLimit(ensembleAnl)
          if ( imposeRttovHuLimits   ) call qlim_rttovLimit     (ensembleAnl)
          ! And recompute analysis mean
          call ens_computeMean(ensembleAnl)
          call ens_copyEnsMean(ensembleAnl, stateVectorMeanAnl)
        end if
      end if

      !- If SubSample requested, copy sub-sample of analysis and trial members
      if (writeSubSample) then
        ! Copy sub-sampled analysis and trial ensemble members
        call epp_selectSubSample(ensembleAnl, ensembleAnlSubSample,  &
                                 ensembleTrl, ensembleTrlSubSample)

        ! Create subdirectory for outputting sub sample increments
        ierr = clib_mkdir_r('subspace')

        ! Allocate stateVectors to store and output sub-sampled ensemble mean analysis
        call gsv_allocate( stateVectorMeanAnlSubSample, tim_nstepobsinc,  &
                           hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                           mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                           dataKind_opt=4, allocHeightSfc_opt=.true., &
                           hInterpolateDegree_opt = hInterpolationDegree, &
                           allocHeight_opt=.false., allocPressure_opt=.false. )
        call gsv_zero(stateVectorMeanAnlSubSample)

      end if

      !- If unperturbed SubSample requested, copy sub-sample of analysis members
      if (writeSubSampleUnPert) then
        ! Copy sub-sampled analysis ensemble members
        call epp_selectSubSample(ensembleAnl, ensembleAnlSubSampleUnPert)

        ! Create subdirectory for outputting sub sample members without perturbations
        ierr = clib_mkdir_r('subspace_unpert')

      end if

      !- Apply random additive inflation, if requested
      if (alphaRandomPert > 0.0D0) then
        ! If namelist value is -999, set random seed using the date (as in standard EnKF)
        if (randomSeed == -999) then
          imode = -3 ! stamp to printable date and time: YYYYMMDD, HHMMSShh
          dateStamp = tim_getDateStamp()
          ierr = newdate(dateStamp, datePrint, timePrint, imode)
          timePrint = timePrint/1000000
          datePrint =  datePrint*100 + timePrint
          if (includeYearInSeed) then
            ! Remove the century, keeping 2 digits of the year
            randomSeedRandomPert = datePrint - 100000000*(datePrint/100000000)
          else
            ! Remove the year and add 9
            randomSeedRandomPert = 9 + datePrint - 1000000*(datePrint/1000000)
          end if
        else
          randomSeedRandomPert = randomSeed
        end if
        write(*,*) 'epp_postProcess: randomSeed for additive inflation set to ', &
             randomSeedRandomPert
        if (ens_isAllocated(ensembleTrl)) then
          call epp_addRandomPert(ensembleAnl, stateVectorMeanTrl, alphaRandomPert, &
               randomSeedRandomPert, useMemberAsHuRefState)
        else
          call epp_addRandomPert(ensembleAnl, stateVectorMeanAnl, alphaRandomPert, &
               randomSeedRandomPert, useMemberAsHuRefState)
        end if
      end if

      !- Recompute the analysis spread stddev after inflation and humidity limits
      call gsv_allocate( stateVectorStdDevAnlPert, tim_nstepobsinc, hco_ens, vco_ens, &
                         dateStamp_opt=tim_getDateStamp(),  &
                         mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                         hInterpolateDegree_opt = hInterpolationDegree, &
                         dataKind_opt=4, allocHeight_opt=.false., &
                         allocPressure_opt=.false. )
      call ens_computeStdDev(ensembleAnl)
      call ens_copyEnsStdDev(ensembleAnl, stateVectorStdDevAnlPert)

      !- If SubSample requested, do remaining processing and output of sub-sampled members
      if (writeSubSample) then

        ! Apply random additive inflation to sub-sampled ensemble, if requested
        if (alphaRandomPertSubSample > 0.0D0) then
          ! If namelist value is -999, set random seed using the date (as in standard EnKF)
          if (randomSeed == -999) then
            imode = -3 ! stamp to printable date and time: YYYYMMDD, HHMMSShh
            dateStamp = tim_getDateStamp()
            ierr = newdate(dateStamp, datePrint, timePrint, imode)
            timePrint = timePrint/1000000
            datePrint =  datePrint*100 + timePrint
            if (includeYearInSeed) then
              ! Remove the century, keeping 2 digits of the year
              randomSeedRandomPert = datePrint - 100000000*(datePrint/100000000)
            else
              ! Remove the year and add 9
              randomSeedRandomPert = 9 + datePrint - 1000000*(datePrint/1000000)
            end if
          else
            randomSeedRandomPert = randomSeed
          end if
          write(*,*) 'epp_postProcess: randomSeed for additive inflation set to ', &
               randomSeedRandomPert
          call epp_addRandomPert(ensembleAnlSubSample, stateVectorMeanTrl,  &
                                 alphaRandomPertSubSample, randomSeedRandomPert, &
                                 useMemberAsHuRefState)
        end if

        ! Compute analysis mean of sub-sampled ensemble
        call ens_computeMean(ensembleAnlSubSample)

        ! Shift members to have same mean as full ensemble
        call ens_recenter(ensembleAnlSubSample, stateVectorMeanAnl,  &
                          recenteringCoeff_opt=1.0D0)

        ! Re-compute analysis mean of sub-sampled ensemble
        call ens_computeMean(ensembleAnlSubSample)
        call ens_copyEnsMean(ensembleAnlSubSample, stateVectorMeanAnlSubSample)
        
      end if ! writeSubsample

      !- If SubSample requested, do remaining processing and output of sub-sampled members
      if (writeSubSampleUnPert) then

        ! Compute analysis mean of sub-sampled ensemble
        call ens_computeMean(ensembleAnlSubSampleUnPert)

        ! Shift members to have same mean as full ensemble
        call ens_recenter(ensembleAnlSubSampleUnPert, stateVectorMeanAnl,  &
                          recenteringCoeff_opt=1.0D0)

      end if

    end if ! ens_isAllocated(ensembleAnl)

    !
    !- Transform data before computing the increments and writing to files.
    !
    if (ens_isAllocated(ensembleAnl)) then
      call gvt_transform(ensembleAnl,'AllTransformedToModel',allowOverWrite_opt=.true.)
      call gvt_transform(stateVectorMeanAnl,'AllTransformedToModel',allowOverWrite_opt=.true.)
      if (writeSubSample) then
        call gvt_transform(ensembleAnlSubSample,'AllTransformedToModel',allowOverWrite_opt=.true.)
        call gvt_transform(stateVectorMeanAnlSubSample,'AllTransformedToModel',allowOverWrite_opt=.true.)
      end if
      if (writeSubSampleUnPert) then
        call gvt_transform(ensembleAnlSubSampleUnPert,'AllTransformedToModel',allowOverWrite_opt=.true.)
      end if
    end if

    ! When we read ensemble trials we always need to transform them either for incremnets or for writing
    if (ens_isAllocated(ensembleTrl)) then
      call gvt_transform(ensembleTrl,'AllTransformedToModel',allowOverWrite_opt=.true.)
      call gvt_transform(stateVectorCtrlTrl,'AllTransformedToModel',allowOverWrite_opt=.true.)
      call gvt_transform(stateVectorMeanTrl,'AllTransformedToModel',allowOverWrite_opt=.true.)
      if (writeSubSample) then
        call gvt_transform(ensembleTrlSubSample,'AllTransformedToModel',allowOverWrite_opt=.true.)
      end if
    end if

    !
    !- Compute increments
    !
    if (ens_isAllocated(ensembleAnl).and.ens_isAllocated(ensembleTrl)) then

      !- Read the analysis mask (in LAM mode only) - N.B. different from land/sea mask!!!
      if (.not. hco_ens%global .and. useAnalIncMask) then
        call gio_getMaskLAM(stateVectorAnalIncMask, hco_ens, vco_ens, hInterpolationDegree)
      end if

      ! initialize the vector for the ensemble increment
      allocate(dateStampListInc(tim_nstepobsinc))
      call tim_getstamplist(dateStampListInc,tim_nstepobsinc,tim_getDatestamp())
      nullify(varNames)
      call ens_varNamesList(varNames,ensembleTrl)
      call ens_allocate(ensembleAnlInc, nEns, tim_nstepobsinc, hco_ens, vco_ens, &
                        dateStampListInc, varNames_opt=varNames)
      call ens_copy(ensembleTrl,ensembleAnlInc)
      deallocate(varNames)

      ! Compute the ensemble increments
      call ens_add(ensembleAnl, ensembleAnlInc, scaleFactorInOut_opt=-1.0D0)
      
      !- Mask the ensemble increments for LAM grid and recompute ensemble analysis
      if (.not. hco_ens%global .and. useAnalIncMask) then
        call ens_applyMaskLAM(ensembleAnlInc, stateVectorAnalIncMask)
        call ens_copy(ensembleAnlInc,ensembleAnl)
        call ens_add(ensembleTrl,ensembleAnl)
      end if

      ! Compute mean increment for converted model variables (e.g. VIS and PR)
      nullify(varNames)
      call gsv_varNamesList(varNames, stateVectorMeanAnl)
      call gsv_allocate( stateVectorMeanInc, tim_nstepobsinc, hco_ens, vco_ens, &
                         dateStamp_opt=tim_getDateStamp(),  &
                         mpi_local_opt=.true., mpi_distribution_opt='Tiles',  &
                         hInterpolateDegree_opt = hInterpolationDegree, &
                         dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=varNames )
      call gsv_copy(stateVectorMeanAnl, stateVectorMeanInc)
      deallocate(varNames)

      call gsv_add(stateVectorCtrlTrl, stateVectorMeanInc, scaleFactor_opt=-1.0D0)

      !- Mask the mean increment for LAM grid and recompute the mean analysis
      if (.not. hco_ens%global .and. useAnalIncMask) then
        call gsv_applyMaskLAM(stateVectorMeanInc, stateVectorAnalIncMask)
        call gsv_copy(stateVectorMeanInc,stateVectorMeanAnl)
        call gsv_add(stateVectorCtrlTrl,stateVectorMeanAnl)
      end if

      ! If subsample reqested compute the subsample increments
      if (writeSubSample) then
        ! Compute ensemble increments for subsample
        nullify(varNames)
        call ens_varNamesList(varNames,ensembleAnlSubSample)
        call ens_allocate(ensembleAnlIncSubSample, ens_getNumMembers(ensembleAnlSubSample), &
                          tim_nstepobsinc, hco_ens, vco_ens, dateStampListInc, varNames_opt=varNames)
        call ens_copy(ensembleTrlSubSample,ensembleAnlIncSubSample)
        deallocate(varNames)
        
        call ens_add(ensembleAnlSubSample,ensembleAnlIncSubSample,scaleFactorInOut_opt=-1.0D0)

        !- Mask the ensemble increment for LAM grid and recompute the mean analysis
        if (.not. hco_ens%global .and. useAnalIncMask) then
          call ens_applyMaskLAM(ensembleAnlIncSubSample, stateVectorAnalIncMask)
          call ens_copy(ensembleAnlIncSubSample,ensembleAnlSubSample)
          call ens_add(ensembleTrlSubSample,ensembleAnlSubSample)
        end if

        ! Compute mean increment with respect to mean of full trial ensemble
        nullify(varNames)
        call gsv_varNamesList(varNames, stateVectorMeanAnlSubSample)
        call gsv_allocate( stateVectorMeanIncSubSample, tim_nstepobsinc,  &
                           hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                           mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                           dataKind_opt=4, allocHeightSfc_opt=.true., &
                           hInterpolateDegree_opt = hInterpolationDegree, &
!                           allocHeight_opt=.false., allocPressure_opt=.false., &
                           varNames_opt=varNames )
        call gsv_copy(stateVectorMeanAnlSubSample, stateVectorMeanIncSubSample)
        deallocate(varNames)

        call gsv_add(stateVectorCtrlTrl, stateVectorMeanIncSubSample, scaleFactor_opt=-1.0D0)

        !- Mask the mean increment for LAM grid and recompute the mean analysis
        if (.not. hco_ens%global .and. useAnalIncMask) then
          call gsv_applyMaskLAM(stateVectorMeanIncSubSample, stateVectorAnalIncMask)
          call gsv_copy(stateVectorMeanIncSubSample,stateVectorMeanAnlSubSample)
          call gsv_add(stateVectorCtrlTrl,stateVectorMeanAnlSubSample)
        end if
      end if

    end if !- end of increment computation

    !
    !- Output everything
    !

    !- Output ens stddev and mean in trialrms, analrms and analpertrms files

    ! determine middle timestep for output of these files
    middleStepIndex = (tim_nstepobsinc + 1) / 2

    if (ens_isAllocated(ensembleTrl)) then
      ! output trialmean, trialrms
      call utl_tmg_start(5,'--WriteEnsMeanRms')
      call fln_ensTrlFileName(outFileName, '.', tim_getDateStamp())
      outFileName = trim(outFileName) // '_trialmean'
      call ens_copyMaskToGsv(ensembleTrl, stateVectorMeanTrl)
      do stepIndex = 1, tim_nstepobsinc
        call gio_writeToFile(stateVectorMeanTrl, outFileName, etiket_trlmean,  &
                             typvar_opt='P', writeHeightSfc_opt=.false., numBits_opt=numBits,  &
                             stepIndex_opt=stepIndex, containsFullField_opt=.true.)
      end do
      call fln_ensTrlFileName(outFileName, '.', tim_getDateStamp())
      outFileName = trim(outFileName) // '_trialrms'
      call ens_copyMaskToGsv(ensembleTrl, stateVectorStdDevTrl)
      call gio_writeToFile(stateVectorStdDevTrl, outFileName, etiket_trlrms,  &
                           typvar_opt='P', writeHeightSfc_opt=.false., numBits_opt=numBits, &
                           stepIndex_opt=middleStepIndex, containsFullField_opt=.false.)
      outFileName = trim(outFileName) // '_ascii'
      call epp_printRmsStats(stateVectorStdDevTrl,outFileName,elapsed=0.0D0,ftype='F',nEns=nEns)
      call utl_tmg_stop(5)

      ! output the trial ensemble if requested (because it was interpolated)
      if (writeTrlEnsemble) then
        call utl_tmg_start(3,'--WriteEnsemble')
        if (.not. outputOnlyEnsMean) then
          call ens_writeEnsemble(ensembleTrl, '.', '', etiket_trl, 'P',  &
                                 numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                                 containsFullField_opt=.true.)
        end if
        call utl_tmg_stop(3)
      end if
    end if

    ! all outputs related to analysis ensemble
    if (ens_isAllocated(ensembleAnl)) then

      !- Prepare stateVector with only MeanAnl surface pressure and surface height
      if (gsv_isAllocated(stateVectorHeightSfc)) then
        call gsv_allocate( stateVectorMeanAnlSfcPres, tim_nstepobsinc, hco_ens, vco_ens,   &
                           dateStamp_opt=tim_getDateStamp(),  &
                           mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                           hInterpolateDegree_opt = hInterpolationDegree, &
                           dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0'/) )
        call gsv_zero(stateVectorMeanAnlSfcPres)
        if (mmpi_myid <= (nEns-1)) then
          call gsv_allocate( stateVectorMeanAnlSfcPresMpiGlb, tim_nstepobsinc, hco_ens, vco_ens,   &
                             dateStamp_opt=tim_getDateStamp(),  &
                             mpi_local_opt=.false., &
                             hInterpolateDegree_opt = hInterpolationDegree, &
                             dataKind_opt=4, allocHeightSfc_opt=.true., varNames_opt=(/'P0'/) )
          call gsv_zero(stateVectorMeanAnlSfcPresMpiGlb)
        end if
        call gsv_copy(stateVectorMeanAnl, stateVectorMeanAnlSfcPres, allowVarMismatch_opt=.true.)
        call gsv_copyHeightSfc(stateVectorHeightSfc, stateVectorMeanAnlSfcPres)
        call gsv_transposeTilesToMpiGlobal(stateVectorMeanAnlSfcPresMpiGlb, stateVectorMeanAnlSfcPres)
      end if
      
      ! output analmean, analrms
      call utl_tmg_start(5,'--WriteEnsMeanRms')
      call fln_ensAnlFileName(outFileName, '.', tim_getDateStamp())
      outFileName = trim(outFileName) // '_analmean'
      call ens_copyMaskToGsv(ensembleAnl, stateVectorMeanAnl)
      do stepIndex = 1, tim_nstepobsinc
        call gio_writeToFile(stateVectorMeanAnl, outFileName, etiket_anlmean,  &
                             typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=numBits, &
                             stepIndex_opt=stepIndex, containsFullField_opt=.true.)
      end do
      call fln_ensAnlFileName(outFileName, '.', tim_getDateStamp())
      outFileName = trim(outFileName) // '_analrms'
      call ens_copyMaskToGsv(ensembleAnl, stateVectorStdDevAnl)
      call gio_writeToFile(stateVectorStdDevAnl, outFileName, etiket_anlrms,  &
                           typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=numBits, &
                           stepIndex_opt=middleStepIndex, containsFullField_opt=.false.)
      outFileName = trim(outFileName) // '_ascii'
      call epp_printRmsStats(stateVectorStdDevAnl,outFileName,elapsed=0.0D0,ftype='A',nEns=nEns)

      ! output analmean_raw and analrms_raw, if requested
      if (writeRawAnalStats) then
        call fln_ensAnlFileName(outFileName, '.', tim_getDateStamp())
        outFileName = trim(outFileName) // '_analmean_raw'
        call ens_copyMaskToGsv(ensembleAnl, stateVectorMeanAnlRaw)
        do stepIndex = 1, tim_nstepobsinc
          call gio_writeToFile(stateVectorMeanAnlRaw, outFileName, etiket_anlmean_raw,  &
                               typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=numBits, &
                               stepIndex_opt=stepIndex, containsFullField_opt=.true.)
        end do
        call fln_ensAnlFileName(outFileName, '.', tim_getDateStamp())
        outFileName = trim(outFileName) // '_analrms_raw'
        call ens_copyMaskToGsv(ensembleAnl, stateVectorStdDevAnl)
        call gio_writeToFile(stateVectorStdDevAnlRaw, outFileName, etiket_anlrms_raw,  &
                             typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=numBits, &
                             stepIndex_opt=middleStepIndex, containsFullField_opt=.false.)
        outFileName = trim(outFileName) // '_ascii'
        call epp_printRmsStats(stateVectorStdDevAnlRaw,outFileName,elapsed=0.0D0,ftype='A',nEns=nEns)
      end if  ! writeRawAnalStats

      if (alphaRandomPert > 0.0D0) then
        ! output analpertmean, analpertrms
        call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp() )
        outFileName = trim(outFileName) // '_analpertmean'
        call ens_copyMaskToGsv(ensembleAnl, stateVectorMeanAnl)
        do stepIndex = 1, tim_nstepobsinc
          call gio_writeToFile(stateVectorMeanAnl, outFileName, etiket_anlmeanpert,  &
                               typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=numBits, &
                               stepIndex_opt=stepIndex, containsFullField_opt=.true.)
        end do
        call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp() )
        outFileName = trim(outFileName) // '_analpertrms'
        call ens_copyMaskToGsv(ensembleAnl, stateVectorStdDevAnlPert)
        call gio_writeToFile(stateVectorStdDevAnlPert, outFileName, etiket_anlrmspert,  &
                             typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=numBits, &
                             stepIndex_opt=middleStepIndex, containsFullField_opt=.false.)
        outFileName = trim(outFileName) // '_ascii'
        call epp_printRmsStats(stateVectorStdDevAnlPert,outFileName,elapsed=0.0D0,ftype='P',nEns=nEns)
      end if
      call utl_tmg_stop(5)

      !- Output the ensemble mean increment (include MeanAnl Psfc) and ensemble increments
      if (ens_isAllocated(ensembleTrl)) then
        call utl_tmg_start(5,'--WriteEnsMeanRms')
        ! output ensemble mean increment
        call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp(), 0, ensFileNameSuffix_opt='inc' )
        call ens_copyMaskToGsv(ensembleAnl, stateVectorMeanInc)
        ! here we assume 4 digits for the ensemble member!!!!
        etiket = trim(etiket_inc) // '0000'
        do stepIndex = 1, tim_nstepobsinc
          call gio_writeToFile(stateVectorMeanInc, outFileName, etiket,  &
                               typvar_opt='R', writeHeightSfc_opt=.false., numBits_opt=numBits, &
                               stepIndex_opt=stepIndex, containsFullField_opt=.false.)
          if (gsv_isAllocated(stateVectorMeanAnlSfcPres)) then
            call gio_writeToFile(stateVectorMeanAnlSfcPres, outFileName, etiket,  &
                                 typvar_opt='A', writeHeightSfc_opt=.true., &
                                 stepIndex_opt=stepIndex, containsFullField_opt=.true.)
          end if
        end do
        call utl_tmg_stop(5)

        !- Output all ensemble member increments
        call utl_tmg_start(3,'--WriteEnsemble')
        if (.not. outputOnlyEnsMean) then
          call ens_writeEnsemble(ensembleAnlInc, '.', '', etiket_inc, 'R',  &
                                 numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                                 containsFullField_opt=.false., resetTimeParams_opt=.true.)
          if (gsv_isAllocated(stateVectorMeanAnlSfcPresMpiGlb)) then
            ! Also write the reference (analysis) surface pressure to increment files
            call epp_writeToAllMembers(stateVectorMeanAnlSfcPresMpiGlb, nEns,  &
                                       etiket='ENS_INC', typvar='A', fileNameSuffix='inc',  &
                                       ensPath='.')
          end if
        end if
        call utl_tmg_stop(3)

      end if ! allocated(ensembleTrl)

      ! output ensemble mean analysis state
      call utl_tmg_start(5,'--WriteEnsMeanRms')
      call fln_ensAnlFileName( outFileName, '.', tim_getDateStamp(), 0 )
      call ens_copyMaskToGsv(ensembleAnl, stateVectorMeanAnl)
      ! here we assume 4 digits for the ensemble member!!!!
      etiket = trim(etiket_anl) // '0000'
      do stepIndex = 1, tim_nstepobsinc
        call gio_writeToFile(stateVectorMeanAnl, outFileName, etiket,  &
                             typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=numBits, &
                             stepIndex_opt=stepIndex, containsFullField_opt=.true.)
      end do
      call utl_tmg_stop(5)

      !- Output all ensemble member analyses
      ! convert transformed to model variables for analysis and trial ensembles
      call utl_tmg_start(3,'--WriteEnsemble')
      if (.not. outputOnlyEnsMean) then
        call ens_writeEnsemble(ensembleAnl, '.', '', etiket_anl, 'A',  &
                               numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                               containsFullField_opt=.true.)
      end if
      call utl_tmg_stop(3)

      !- Output the sub-sampled ensemble analyses and increments
      if (writeSubSample) then
        call utl_tmg_start(5,'--WriteEnsMeanRms')
        ! Output the ensemble mean increment (include MeanAnl Psfc)
        call fln_ensAnlFileName( outFileName, 'subspace', tim_getDateStamp(), 0, ensFileNameSuffix_opt='inc' )
        ! here we assume 4 digits for the ensemble member!!!!
        etiket = trim(etiket_inc) // '0000'
        do stepIndex = 1, tim_nstepobsinc
          call gio_writeToFile(stateVectorMeanIncSubSample, outFileName, etiket,  &
                               typvar_opt='R', writeHeightSfc_opt=.false., numBits_opt=numBits, &
                               stepIndex_opt=stepIndex, containsFullField_opt=.false.)
          if (gsv_isAllocated(stateVectorMeanAnlSfcPres)) then
            call gio_writeToFile(stateVectorMeanAnlSfcPres, outFileName, etiket,  &
                                 typvar_opt='A', writeHeightSfc_opt=.true., &
                                 stepIndex_opt=stepIndex, containsFullField_opt=.true.)
          end if
        end do

        ! Output the ensemble mean analysis state
        call fln_ensAnlFileName( outFileName, 'subspace', tim_getDateStamp(), 0 )
        ! here we assume 4 digits for the ensemble member!!!!
        etiket = trim(etiket_anl) // '0000'
        do stepIndex = 1, tim_nstepobsinc
          call gio_writeToFile(stateVectorMeanAnlSubSample, outFileName, etiket,  &
                               typvar_opt='A', writeHeightSfc_opt=.false., numBits_opt=numBits, &
                               stepIndex_opt=stepIndex, containsFullField_opt=.true.)
        end do
        call utl_tmg_stop(3)

        ! Output the sub-sampled analysis ensemble members
        call utl_tmg_start(3,'--WriteEnsemble')
        if (.not. outputOnlyEnsMean) then
          call ens_writeEnsemble(ensembleAnlSubSample, 'subspace', '', etiket_anl, 'A',  &
                                 numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                                 containsFullField_opt=.true.)
        end if
        call utl_tmg_stop(3)

        ! Output the sub-sampled ensemble increments (include MeanAnl Psfc)
        call utl_tmg_start(3,'--WriteEnsemble')
        if (.not. outputOnlyEnsMean) then
          call ens_writeEnsemble(ensembleAnlIncSubSample, 'subspace', '', etiket_inc, 'R',  &
                                 numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                                 containsFullField_opt=.false., resetTimeParams_opt=.true.)
          ! Also write the reference (analysis) surface pressure to increment files
          if (gsv_isAllocated(stateVectorMeanAnlSfcPresMpiGlb)) then
            call epp_writeToAllMembers(stateVectorMeanAnlSfcPresMpiGlb,  &
                                       ens_getNumMembers(ensembleAnlSubSample),  &
                                       etiket=etiket_inc, typvar='A', fileNameSuffix='inc',  &
                                       ensPath='subspace')
          end if
        end if
        call utl_tmg_stop(3)

      end if ! writeSubSample

      !- Output the unperturbed sub-sampled ensemble analyses
      if (writeSubSampleUnPert) then

        ! Output the sub-sampled analysis ensemble members
        call utl_tmg_start(3,'--WriteEnsemble')
        if (.not. outputOnlyEnsMean) then
          call ens_writeEnsemble(ensembleAnlSubSampleUnPert, 'subspace_unpert', '', etiket_anl, 'A',  &
                                 numBits_opt=16, etiketAppendMemberNumber_opt=.true.,  &
                                 containsFullField_opt=.true.)
        end if
        call utl_tmg_stop(3)

      end if

    end if ! ens_isAllocated(ensembleAnl)

  end subroutine epp_postProcess

  !----------------------------------------------------------------------
  ! epp_writeToAllMembers (private subroutine)
  !----------------------------------------------------------------------
  subroutine epp_writeToAllMembers(stateVector, nEns, etiket, typvar,  &
                                    fileNameSuffix, ensPath)
    !
    !:Purpose:   Write the contents of the supplied stateVector to all
    !            ensemble member files in an efficient parallel way.
    !
    implicit none

    ! arguments:
    type(struct_gsv) :: stateVector
    integer          :: nEns
    character(len=*) :: etiket
    character(len=*) :: typvar
    character(len=*) :: fileNameSuffix
    character(len=*) :: ensPath

    ! locals:
    integer            :: memberIndex, stepIndex, writeFilePE(nEns)
    character(len=4)   :: memberIndexStr
    character(len=256) :: outFileName

    do memberIndex = 1, nEns
      writeFilePE(memberIndex) = mod(memberIndex-1, mmpi_nprocs)
    end do

    do memberIndex = 1, nEns

      if (mmpi_myid == writeFilePE(memberIndex)) then

        call fln_ensAnlFileName( outFileName, ensPath, tim_getDateStamp(),  &
                                 memberIndex, ensFileNameSuffix_opt=fileNameSuffix )
        write(memberIndexStr,'(I4.4)') memberIndex

        do stepIndex = 1, tim_nstepobsinc
          call gio_writeToFile(stateVector, outFileName,  &
                               trim(etiket) // memberIndexStr,  &
                               typvar_opt=trim(typvar), writeHeightSfc_opt=.true., &
                               stepIndex_opt=stepIndex, containsFullField_opt=.true., &
                               numBits_opt=16)
        end do

      end if

    end do

  end subroutine epp_writeToAllMembers

  !--------------------------------------------------------------------------
  ! epp_RTPS
  !--------------------------------------------------------------------------
  subroutine epp_RTPS(ensembleAnl, stateVectorStdDevAnl, stateVectorStdDevTrl, &
                      stateVectorMeanAnl, alphaRTPS)
    ! :Purpose: Apply Relaxation To Prior Spread ensemble inflation according
    !           to the factor alphaRTPS (usually between 0 and 1).
    implicit none

    ! Arguments
    type(struct_ens) :: ensembleAnl
    type(struct_gsv) :: stateVectorStdDevAnl
    type(struct_gsv) :: stateVectorStdDevTrl
    type(struct_gsv) :: stateVectorMeanAnl
    real(8)          :: alphaRTPS

    ! Locals
    integer :: varLevIndex, latIndex, lonIndex, stepIndex, memberIndex
    integer :: nEns, numVarLev, myLonBeg, myLonEnd, myLatBeg, myLatEnd
    real(8) :: factorRTPS
    real(4), pointer     :: stdDevTrl_ptr_r4(:,:,:,:), stdDevAnl_ptr_r4(:,:,:,:)
    real(4), pointer     :: meanAnl_ptr_r4(:,:,:,:), memberAnl_ptr_r4(:,:,:,:)

    write(*,*) 'epp_RTPS: Starting'

    call gsv_getField(stateVectorStdDevTrl,stdDevTrl_ptr_r4)
    call gsv_getField(stateVectorStdDevAnl,stdDevAnl_ptr_r4)
    call gsv_getField(stateVectorMeanAnl,meanAnl_ptr_r4)

    nEns = ens_getNumMembers(ensembleAnl)
    numVarLev = ens_getNumK(ensembleAnl)
    call ens_getLatLonBounds(ensembleAnl, myLonBeg, myLonEnd, myLatBeg, myLatEnd)
    do varLevIndex = 1, numVarLev
      memberAnl_ptr_r4 => ens_getOneLev_r4(ensembleAnl,varLevIndex)
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          do stepIndex = 1, tim_nstepobsinc
            ! compute the inflation factor for RTPS
            if ( stdDevAnl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) > 0.0 ) then
              factorRTPS = 1.0D0 + alphaRTPS *  &
                           ( stdDevTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) -  &
                             stdDevAnl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) ) /  &
                           stdDevAnl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex)
            else
              factorRTPS = 0.0D0
            end if
            ! apply the inflation factor to all Anl members (in place)
            if (factorRTPS > 0.0D0) then
              do memberIndex = 1, nEns
                memberAnl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) =  &
                     factorRTPS *  &
                     ( memberAnl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) -  &
                       meanAnl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) ) +  &
                     meanAnl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex)
              end do ! memberIndex
            end if ! factorRTPS > 0
          end do ! stepIndex
        end do ! lonIndex
      end do ! latIndex
    end do ! varLevIndex

    write(*,*) 'epp_RTPS: Finished'

  end subroutine epp_RTPS

  !--------------------------------------------------------------------------
  ! epp_RTPP
  !--------------------------------------------------------------------------
  subroutine epp_RTPP(ensembleAnl, ensembleTrl, stateVectorMeanAnl, &
                      stateVectorMeanTrl, alphaRTPP)
    ! :Purpose: Apply Relaxation To Prior Perturbation ensemble inflation according
    !           to the factor alphaRTPP (usually between 0 and 1).
    implicit none

    ! Arguments
    type(struct_ens) :: ensembleAnl
    type(struct_ens) :: ensembleTrl
    type(struct_gsv) :: stateVectorMeanAnl
    type(struct_gsv) :: stateVectorMeanTrl
    real(8)          :: alphaRTPP

    ! Locals
    integer :: varLevIndex, latIndex, lonIndex, stepIndex, memberIndex
    integer :: nEns, numVarLev, myLonBeg, myLonEnd, myLatBeg, myLatEnd
    real(4), pointer     :: meanAnl_ptr_r4(:,:,:,:), meanTrl_ptr_r4(:,:,:,:)
    real(4), pointer     :: memberAnl_ptr_r4(:,:,:,:), memberTrl_ptr_r4(:,:,:,:)

    write(*,*) 'epp_RTPP: Starting'

    call gsv_getField(stateVectorMeanAnl, meanAnl_ptr_r4)
    call gsv_getField(stateVectorMeanTrl, meanTrl_ptr_r4)

    nEns = ens_getNumMembers(ensembleAnl)
    numVarLev = ens_getNumK(ensembleAnl)
    call ens_getLatLonBounds(ensembleAnl, myLonBeg, myLonEnd, myLatBeg, myLatEnd)

    do varLevIndex = 1, numVarLev
      memberAnl_ptr_r4 => ens_getOneLev_r4(ensembleAnl,varLevIndex)
      memberTrl_ptr_r4 => ens_getOneLev_r4(ensembleTrl,varLevIndex)
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          do stepIndex = 1, tim_nstepobsinc
            ! apply RTPP to all Anl members (in place)
            ! member_anl = mean_anl + (1-a)*(member_anl-mean_anl) + a*(member_trl-mean_trl)
            do memberIndex = 1, nEns
              memberAnl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) =  &
                   meanAnl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) + &
                   (1.0D0 - alphaRTPP) * &
                   ( memberAnl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) -  &
                     meanAnl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) ) +  &
                   alphaRTPP * &
                   ( memberTrl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) -  &
                     meanTrl_ptr_r4(lonIndex,latIndex,varLevIndex,stepIndex) )
            end do ! memberIndex
          end do ! stepIndex
        end do ! lonIndex
      end do ! latIndex
    end do ! varLevIndex

    write(*,*) 'epp_RTPP: Finished'

  end subroutine epp_RTPP

  !--------------------------------------------------------------------------
  ! epp_addRandomPert
  !--------------------------------------------------------------------------
  subroutine epp_addRandomPert(ensembleAnl, stateVectorRefState, alphaRandomPert, &
                               randomSeed, useMemberAsHuRefState)
    ! :Purpose: Apply additive inflation using random perturbations from sampling
    !           the B matrix as defined by the regular namelist block NAMBHI, NAMBEN, etc.
    !           The scale factor alphaRandomPert (usually between 0 and 1) is used
    !           to simply multiply the resulting perturbations before adding to the original
    !           members. The perturbations have zero ensemble mean.
    implicit none

    ! Arguments
    type(struct_ens), intent(inout) :: ensembleAnl
    type(struct_gsv), intent(in)    :: stateVectorRefState
    real(8)         , intent(in)    :: alphaRandomPert
    integer         , intent(in)    :: randomSeed
    logical         , intent(in)    :: useMemberAsHuRefState

    ! Locals
    type(struct_gsv)         :: stateVectorPerturbation
    type(struct_gsv)         :: stateVectorPerturbationInterp
    type(struct_gsv)         :: stateVectorHuRefState
    type(struct_gsv)         :: stateVectorHuRefStateInterp
    type(struct_vco), pointer :: vco_randomPert, vco_ens
    type(struct_hco), pointer :: hco_randomPert, hco_ens, hco_core
    character(len=12)  :: etiket
    real(8), allocatable :: controlVector_mpiglobal(:), controlVector(:)
    real(8), allocatable :: perturbationMean(:,:,:)
    real(8), allocatable :: PsfcReference(:,:,:)
    real(8), pointer     :: perturbation_ptr(:,:,:)
    real(4), pointer     :: memberAnl_ptr_r4(:,:,:,:)
    integer :: cvIndex, memberIndex, varLevIndex, lonIndex, latIndex, stepIndex
    integer :: nEns, numVarLev, myLonBeg, myLonEnd, myLatBeg, myLatEnd, varIndex
    integer :: middleStepIndex
    logical, save :: firstCall = .true.
    character(len=4), pointer :: varNamesWithLQ(:)
    character(len=4) :: varName

    call utl_tmg_start(4,'--AddEnsRandomPert')

    ! Determine middle timestep
    middleStepIndex = (tim_nstepobsinc + 1) / 2

    ! Get ensemble dimensions
    nEns = ens_getNumMembers(ensembleAnl)
    numVarLev = ens_getNumK(ensembleAnl)
    call ens_getLatLonBounds(ensembleAnl, myLonBeg, myLonEnd, myLatBeg, myLatEnd)
    vco_ens => ens_getVco(ensembleAnl)
    hco_ens => ens_getHco(ensembleAnl)

    ! Define the horiz/vertical coordinate for perturbation calculation
    nullify(vco_randomPert)
    nullify(hco_randomPert)
    call hco_setupFromFile(hco_randomPert, './analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN
    if ( hco_randomPert % global ) then
      etiket = 'BGCK_STDDEV'
    else
      etiket = 'STDDEV'
    end if
    call vco_setupFromFile(vco_randomPert, './bgcov', etiket)

    ! Special modification to vco_randomPert if a soil variables is present
    if ( vco_getNumLev(vco_ens,'OT','I0') > 0 .or. &
         vco_getNumLev(vco_ens,'OT','I1') > 0 ) then
      write(*,*)
      write(*,*) 'epp_addRandomPert: copying number of levels for land surface variables!' 
      write(*,*)
      vco_randomPert%nlev_Other(:) = vco_ens%nlev_Other(:)
    end if

    ! Get list of variable names in the ensemble and modify if necessary
    nullify(varNamesWithLQ)
    call ens_varNamesList(varNamesWithLQ,ensembleAnl)
    if (useMemberAsHuRefState .and. ens_varExist(ensembleAnl,'HU')) then
      ! Replace HU with LQ so we can apply the transform directly in this routine
      do varIndex = 1, size(varNamesWithLQ)
        if (varNamesWithLQ(varIndex) == 'HU') varNamesWithLQ(varIndex) = 'LQ'
      end do
      if (mmpi_myid == 0) write(*,*) 'epp_addRandomPert: varNamesWithLQ = ', varNamesWithLQ(:)
    end if

    hco_core => hco_randomPert
    if (firstCall) then
      call utl_tmg_stop(4) ! stop counter, since Bmat has it's own counters
      call bmat_setup(hco_randomPert, hco_core, vco_randomPert)
      firstCall = .false.
      call utl_tmg_start(4,'--AddEnsRandomPert')
    end if
    call gvt_setup(hco_randomPert, hco_core, vco_randomPert)

    call rng_setup(abs(randomSeed))

    allocate(controlVector(cvm_nvadim))
    allocate(controlVector_mpiglobal(cvm_nvadim_mpiglobal))
    allocate(perturbationMean(myLonBeg:myLonEnd,myLatBeg:myLatEnd,numVarLev))
    perturbationMean(:,:,:) = 0.0d0

    ! NOTE: The following stateVectors will include LQ instead of HU when
    !       useMemberAsHuRefState=true. This results in the perturbations
    !       being provided by bmat_sqrtB in terms of LQ, allowing the conversion
    !       from LQ to HU to be done within this subroutine using the HU field
    !       of each ensemble member.
    call gsv_allocate(stateVectorPerturbation, 1, hco_randomPert, vco_randomPert, &
                      dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                      hInterpolateDegree_opt='LINEAR', &
                      varNames_opt=varNamesWithLQ)
    call gsv_allocate(stateVectorPerturbationInterp, 1, hco_ens, vco_ens, &
                      dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                      hInterpolateDegree_opt='LINEAR', &
                      varNames_opt=varNamesWithLQ)
    call gsv_getField(stateVectorPerturbationInterp,perturbation_ptr)
    allocate(PsfcReference(myLonBeg:myLonEnd,myLatBeg:myLatEnd,1))
    PsfcReference(:,:,:) = 100000.0D0

    ! prepare the reference state HU field for transforming LQ to HU perturbations
    if (ens_varExist(ensembleAnl,'HU')) then
      call gsv_allocate(stateVectorHuRefState, 1, hco_ens, vco_ens,   &
                        dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                        allocHeightSfc_opt=.true., hInterpolateDegree_opt='LINEAR', &
                        hExtrapolateDegree_opt='MINIMUM', &
                        varNames_opt=(/'HU','P0'/) )
      if (.not. useMemberAsHuRefState) then
        ! use the single provided state as the reference state and interpolate to perturbation grid
        call gsv_copy(stateVectorRefState, stateVectorHuRefState, allowVarMismatch_opt=.true.)
        call gsv_allocate(stateVectorHuRefStateInterp, 1, hco_randomPert, vco_randomPert,   &
                          dateStamp_opt=tim_getDateStamp(), mpi_local_opt=.true., &
                          allocHeightSfc_opt=.true., hInterpolateDegree_opt='LINEAR', &
                          hExtrapolateDegree_opt='MINIMUM', &
                          varNames_opt=(/'HU','P0'/) )
        call int_interp_gsv(stateVectorHuRefState, stateVectorHuRefStateInterp)
        call gsv_deallocate(stateVectorHuRefState)
      end if
    end if

    do memberIndex = 1, nEns

      if( mmpi_myid == 0 ) then
        write(*,*) 
        write(*,*) 'Computing random perturbation number= ', memberIndex
      end if

      ! global vector random control vector (independent of mpi topology)
      do cvIndex = 1, cvm_nvadim_mpiglobal
        controlVector_mpiglobal(cvIndex) = rng_gaussian()
      end do
      call bmat_reduceToMPILocal( controlVector, controlVector_mpiglobal )

      call utl_tmg_stop(4) ! stop counter, since Bmat has it's own counters
      if (ens_varExist(ensembleAnl,'HU') .and. .not.useMemberAsHuRefState) then
        ! Use supplied reference state for LQ to HU conversion
        call bmat_sqrtB(controlVector, cvm_nvadim, &       ! IN
                        stateVectorPerturbation,   &       ! OUT
                        stateVectorRef_opt=stateVectorHuRefStateInterp) ! IN
      else
        ! No conversion from LQ to HU done in bmat_sqrtB
        call bmat_sqrtB(controlVector, cvm_nvadim, &   ! IN
                        stateVectorPerturbation)       ! OUT
      end if
      call utl_tmg_start(4,'--AddEnsRandomPert')

      call int_interp_gsv(stateVectorPerturbation, stateVectorPerturbationInterp, &
                          PsfcReference_opt=PsfcReference)

      ! If desired, use member itself as reference state for LQ to HU conversion
      if (ens_varExist(ensembleAnl,'HU') .and. useMemberAsHuRefState) then
        call ens_copyMember(ensembleAnl, stateVectorHuRefState, memberIndex)
        call gvt_transform( stateVectorPerturbationInterp,  &          ! INOUT
                            'LQtoHU_tlm', &                            ! IN
                            stateVectorRef_opt=stateVectorHuRefState ) ! IN
      end if

      ! scale the perturbation by the specified factor
      call gsv_scale(stateVectorPerturbationInterp, alphaRandomPert)

      write(*,*) 'epp_addRandomPert: member ', memberIndex, ', perturbation min/maxval = ',  &
                 minval(perturbation_ptr), maxval(perturbation_ptr)

      do varLevIndex = 1, numVarLev
        varName = ens_getVarNameFromK(ensembleAnl,varLevIndex)
        memberAnl_ptr_r4 => ens_getOneLev_r4(ensembleAnl,varLevIndex)
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd

            perturbationMean(lonIndex, latIndex, varLevIndex) =   &
                 perturbationMean(lonIndex, latIndex, varLevIndex) +  &
                 perturbation_ptr(lonIndex, latIndex, varLevIndex) / real(nEns, 8)

            ! Add perturbation to member
            do stepIndex = 1, tim_nstepobsinc
              memberAnl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) =  &
                   memberAnl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) +  &
                   perturbation_ptr(lonIndex, latIndex, varLevIndex)
            end do ! stepIndex

          end do ! lonIndex
        end do ! latIndex
      end do ! varLevIndex

    end do ! memberIndex

    ! remove the ensemble mean of the perturbations
    do varLevIndex = 1, numVarLev
      memberAnl_ptr_r4 => ens_getOneLev_r4(ensembleAnl,varLevIndex)
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          do stepIndex = 1, tim_nstepobsinc
            do memberIndex = 1, nEns
              memberAnl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) =  &
                   memberAnl_ptr_r4(memberIndex,stepIndex,lonIndex,latIndex) -  &
                   perturbationMean(lonIndex, latIndex, varLevIndex)
            end do ! memberIndex
          end do ! stepIndex
        end do ! lonIndex
      end do ! latIndex
    end do ! varLevIndex

    deallocate(controlVector)
    deallocate(controlVector_mpiglobal)
    deallocate(perturbationMean)
    call gsv_deallocate(stateVectorPerturbation)
    call gsv_deallocate(stateVectorPerturbationInterp)
    deallocate(PsfcReference)
    if (gsv_isAllocated(stateVectorHuRefStateInterp)) then
      call gsv_deallocate(stateVectorHuRefStateInterp)
    end if

    call utl_tmg_stop(4)

  end subroutine epp_addRandomPert

  !-----------------------------------------------------------------
  ! epp_selectSubSample
  !-----------------------------------------------------------------
  subroutine epp_selectSubSample(ensembleAnl, ensembleAnlSubSample,  &
                                  ensembleTrl_opt, ensembleTrlSubSample_opt)
    ! :Purpose: Create sub-sampled ensembles of analyses and trials based on
    !           the contents of the ascii files 'sampletable' which lists
    !           the member indices for the subsample.
    implicit none

    ! Arguments
    type(struct_ens)           :: ensembleAnl
    type(struct_ens)           :: ensembleAnlSubSample
    type(struct_ens), optional :: ensembleTrl_opt
    type(struct_ens), optional :: ensembleTrlSubSample_opt

    ! Locals
    type(struct_gsv) :: stateVectorMember
    integer :: nulFile, ierr, status, numSubSample
    integer :: memberIndex, memberIndexSubSample, memberIndexFull
    integer :: memberIndexesSubSample(1000), memberIndexesFull(1000)
    integer, allocatable :: dateStampListInc(:)
    integer, external :: fnom, fclos

    numSubSample = 0

    nulFile = 0
    ierr = fnom(nulFile, './sampletable', 'FTN+SEQ+R/O', 0)
    do
      read(nulFile,*, IOSTAT=status) memberIndexSubSample, memberIndexFull
      if (status < 0) exit

      numSubSample = numSubSample + 1
      write(*,*) 'epp_selectSubSample: ', memberIndexSubSample, memberIndexFull
      memberIndexesSubSample(numSubSample) = memberIndexSubSample
      memberIndexesFull(numSubSample) = memberIndexFull
    end do
    ierr = fclos(nulFile)

    write(*,*) 'epp_selectSubSample: number of subSample members = ', numSubSample

    allocate(dateStampListInc(tim_nstepobsinc))
    call tim_getstamplist(dateStampListInc,tim_nstepobsinc,tim_getDatestamp())

    call ens_allocate(ensembleAnlSubSample, numSubSample, tim_nstepobsinc,  &
                      ens_getHco(ensembleAnl), ens_getVco(ensembleAnl), dateStampListInc)
    if (present(ensembleTrlSubSample_opt)) then
      call ens_allocate(ensembleTrlSubSample_opt, numSubSample, tim_nstepobsinc,  &
                        ens_getHco(ensembleAnl), ens_getVco(ensembleAnl), dateStampListInc)
    end if

    call gsv_allocate(stateVectorMember, tim_nstepobsinc,  &
                      ens_getHco(ensembleAnl), ens_getVco(ensembleAnl),  &
                      dateStamp_opt=tim_getDateStamp(),  &
                      mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                      dataKind_opt=4, allocHeightSfc_opt=.false., &
                      allocHeight_opt=.false., allocPressure_opt=.false. )

    do memberIndex = 1, numSubSample
      ! copy analysis ensemble member
      call ens_copyMember(ensembleAnl, stateVectorMember,  &
                          memberIndexesFull(memberIndex))
      call ens_insertMember(ensembleAnlSubSample, stateVectorMember,  &
                            memberIndexesSubSample(memberIndex))

      if (present(ensembleTrl_opt) .and. present(ensembleTrlSubSample_opt)) then
        ! copy trial ensemble member
        call ens_copyMember(ensembleTrl_opt, stateVectorMember,  &
                            memberIndexesFull(memberIndex))
        call ens_insertMember(ensembleTrlSubSample_opt, stateVectorMember,  &
                              memberIndexesSubSample(memberIndex))
      end if
    end do

    call gsv_deallocate(stateVectorMember)
    deallocate(dateStampListInc)

  end subroutine epp_selectSubSample

  !-----------------------------------------------------------------
  ! epp_hybridRecentering
  !-----------------------------------------------------------------
  subroutine epp_hybridRecentering(ensembleAnl, weightRecenter, weightRecenterLand, &
                                   useOptionTableRecenter, numMembersToRecenter)
    ! :Purpose: Modify an ensemble by recentering the members on a state provided
    !           in the file "recentering_analysis".
    !           The "weightRecenter" and "numMembersToRecenter" are used in the calculation
    !           to determine the amount of recentering and how many members it is
    !           applied to. Alternatively the information in the "optiontable" file
    !           can be used to perform a different amount of recentering on each
    !           member.
    implicit none

    ! Arguments
    type(struct_ens) :: ensembleAnl
    real(8)          :: weightRecenter
    real(8)          :: weightRecenterLand
    logical          :: useOptionTableRecenter
    integer          :: numMembersToRecenter

    ! Locals
    type(struct_gsv) :: stateVectorRecenterAnl
    type(struct_hco), pointer :: hco_ens => null()
    type(struct_vco), pointer :: vco_ens => null()
    character(len=30)    :: recenterAnlFileName = 'recentering_analysis'
    character(len=20)    :: stringArray(100)
    character(len=1000)  :: textLine
    integer              :: stepIndex, memberIndex, columnIndex
    integer              :: numMembers, numColumns, nulFile, status
    logical              :: recenterAnlFileExists
    real(8), allocatable :: weightArray(:)
    integer, external    :: fnom, fclos

    ! check if recentering analysis file exists
    inquire(file=recenterAnlFileName, exist=recenterAnlFileExists)
    if (.not. recenterAnlFileExists) then
      write(*,*) 'epp_hybridRecentering: RecenterAnlFileName = ', recenterAnlFileName
      call utl_abort('epp_hybridRecentering: The recentering analysis file does not exist')
    end if

    numMembers = ens_getNumMembers(ensembleAnl)

    ! read the optiontable file, if requested
    if (useOptionTableRecenter) then
      write(*,*) 'epp_hybridRecentering: using optiontable file to specify recentering weights.'
      nulFile = 0
      status = fnom(nulFile, './optiontable', 'FMT+SEQ+R/O', 0)
      read(nulFile,'(a)', IOSTAT=status) textLine
      if (status /= 0) then
        call utl_abort('epp_hybridRecentering: unable to read optiontable file')
      end if
      call utl_parseColumns(textLine, numColumns)
      if (mmpi_myid==0) write(*,*) 'epp_hybridRecentering: optiontable file has ', numColumns, ' columns.'
      allocate( weightArray(0:numMembers) )
      rewind(nulFile)
      do memberIndex = 0, numMembers
        read(nulFile,'(a)') textLine
        call utl_parseColumns(textLine, numColumns, stringArray_opt=stringArray)
        if (mmpi_myid==0) write(*,*) memberIndex, (stringArray(columnIndex),columnIndex=1,numColumns)
        read(stringArray(numColumns),'(f6.3)') weightArray(memberIndex)
        if (mmpi_myid==0) write(*,*) 'weightArray = ', weightArray(memberIndex)
      end do
      status = fclos(nulFile)
    end if

    ! allocate and read in recentering analysis state
    hco_ens => ens_getHco(ensembleAnl)
    vco_ens => ens_getVco(ensembleAnl)
    call gsv_allocate( stateVectorRecenterAnl, tim_nstepobsinc, hco_ens, vco_ens, dateStamp_opt=tim_getDateStamp(),  &
                       mpi_local_opt=.true., mpi_distribution_opt='Tiles', &
                       dataKind_opt=4, allocHeightSfc_opt=.false., &
                       hInterpolateDegree_opt = 'LINEAR', &
                       allocHeight_opt=.false., allocPressure_opt=.false. )
    call gsv_zero(stateVectorRecenterAnl)

    do stepIndex = 1, tim_nstepobsinc
      call gio_readFromFile( stateVectorRecenterAnl, recenterAnlFileName, ' ', ' ',  &
                             stepIndex_opt=stepIndex, containsFullField_opt=.true., &
                             readHeightSfc_opt=.false. )
    end do

    ! apply recentering
    if (useOptionTableRecenter) then
      call ens_recenter(ensembleAnl, stateVectorRecenterAnl, &
                        recenteringCoeffArray_opt=weightArray(1:numMembers),  &
                        numMembersToRecenter_opt=numMembersToRecenter, &
                        recenteringCoeffLand_opt=weightRecenterLand)
    else
      call ens_recenter(ensembleAnl, stateVectorRecenterAnl,  &
                        recenteringCoeff_opt=weightRecenter,  &
                        numMembersToRecenter_opt=numMembersToRecenter, &
                        recenteringCoeffLand_opt=weightRecenterLand)
    end if

    call gsv_deallocate(stateVectorRecenterAnl)
    if ( allocated(weightArray) ) deallocate(weightArray)

  end subroutine epp_hybridRecentering

  !-----------------------------------------------------------------
  ! epp_printRmsStats
  !-----------------------------------------------------------------
  subroutine epp_printRmsStats(stateVectorStdDev,fileName,elapsed,ftype,nEns)
    !
    ! :Purpose: Print statistics of a field to an ASCII output file
    !
    implicit none

    ! arguments
    type(struct_gsv)             :: stateVectorStdDev
    character(len=*)             :: fileName
    real(8), intent(in)          :: elapsed
    character(len=1), intent(in) :: ftype
    integer, intent(in)          :: nEns

    ! locals
    real(8), allocatable          :: rmsvalue(:) 
    type(struct_vco), pointer     :: vco
    type(struct_hco), pointer     :: hco
    character(len=4), allocatable :: nomvar_v(:)
    character(len=4)              :: varLevel
    real(8), allocatable          :: weight(:,:), scaleFactor(:)
    real(4), allocatable          :: pressureOrDepth(:)
    integer :: ierr, latIndex, lonIndex, nulFile
    integer :: kIndex, kIndexCount, levIndex, numK, nLev_M, kIndexUU, kIndexVV
    real(4), pointer              :: stdDev_ptr_r4(:,:,:)
    real(8)                       :: pSfc(1,1)
    real(8), pointer              :: pressures_T(:,:,:), pressures_M(:,:,:)
    integer, external             :: fnom, fclos

    vco => gsv_getVco(stateVectorStdDev)
    nLev_M = vco_getNumLev(vco,'MM')

    numK = gsv_getNumK(stateVectorStdDev)
    allocate(nomvar_v(numK))
    allocate(pressureOrDepth(numK))
    allocate(rmsvalue(numK))
    allocate(scaleFactor(numK))
    allocate(weight(stateVectorStdDev%ni,stateVectorStdDev%nj))
    weight(:,:) = 0.0d0

    ! compute a 2D weight field used for horizontal averaging
    hco => gsv_getHco(stateVectorStdDev)
    ! weights are reduced in the overlap region in case of a Yin-Yang grid
    call hco_weight(hco, weight)
    ! compute global mean variance accounting for weights
    call gsv_getField(stateVectorStdDev,stdDev_ptr_r4)
    rmsvalue(:) = 0.0D0
    do kIndex = 1, numK
      do latIndex = stateVectorStdDev%myLatBeg, stateVectorStdDev%myLatEnd
        do lonIndex = stateVectorStdDev%myLonBeg, stateVectorStdDev%myLonEnd
          rmsvalue(kIndex) = rmsvalue(kIndex) +  &
               (stdDev_ptr_r4(lonIndex,latIndex,kIndex)**2) *  &
               weight(lonIndex,latIndex)
        end do
      end do
      call mmpi_allreduce_sumreal8scalar(rmsvalue(kIndex),'GRID')
      rmsvalue(kIndex) = rmsvalue(kIndex)**0.5
    end do

    if (vco%vgridPresent) then
      ! compute pressure for a column where Psfc=1000hPa
      pSfc(1,1) = 1000.0D2 !1000 hPa
      ! pressure on momentum levels
      nullify(pressures_M)
      ierr = vgd_levels(vco%vgrid,           &
                        ip1_list=vco%ip1_M,  &
                        levels=pressures_M,  &
                        sfc_field=pSfc,      &
                        in_log=.false.)
      if ( ierr /= VGD_OK ) call utl_abort('epp_printRmsStats: ERROR with vgd_levels')
      ! pressure on thermodynamic levels
      nullify(pressures_T)
      ierr = vgd_levels(vco%vgrid,           &
                        ip1_list=vco%ip1_T,  &
                        levels=pressures_T,  &
                        sfc_field=pSfc,      &
                        in_log=.false.)
      if ( ierr /= VGD_OK ) call utl_abort('epp_printRmsStats: ERROR with vgd_levels')

      ! set the variable name and pressure for each element of column
      do kIndex = 1, numK
        levIndex = gsv_getLevFromK(stateVectorStdDev, kIndex)
        nomvar_v(kIndex) = gsv_getVarNameFromK(stateVectorStdDev,kIndex)
        varLevel = vnl_varLevelFromVarname(nomvar_v(kIndex))
        if (varLevel == 'MM') then
          pressureOrDepth(kIndex) = pressures_M(1,1,levIndex)/Psfc(1,1)
        else if (varLevel == 'TH') then
          pressureOrDepth(kIndex) = pressures_T(1,1,levIndex)/Psfc(1,1)
        else
          pressureOrDepth(kIndex) = 1.0
        end if
      end do
      deallocate(pressures_M)
      deallocate(pressures_T)

    else if (vco%nLev_depth > 0) then

      ! set the variable name and depth for each element of column
      do kIndex = 1, numK
        nomvar_v(kIndex) = gsv_getVarNameFromK(stateVectorStdDev,kIndex)
        if (vnl_varLevelFromVarName(nomvar_v(kIndex)) == 'SSDP') then
          pressureOrDepth(kIndex) = 0.0
        else
          levIndex = gsv_getLevFromK(stateVectorStdDev, kIndex)
          pressureOrDepth(kIndex) = vco%depths(levIndex)
        end if
      end do

    else

      call utl_abort('epp_printRmsStats: Unknown type of levels')

    end if

    ! set the scaleFactor and rmsvalue for each element of column
    do kIndex = 1, numK
      if ( (nomvar_v(kIndex) == 'UU') .or. (nomvar_v(kIndex) == 'VV') ) then
        scaleFactor(kIndex) = MPC_KNOTS_PER_M_PER_S_R8
      else if (nomvar_v(kIndex) == 'P0') then
        scaleFactor(kIndex) = MPC_MBAR_PER_PA_R8
      else
        scaleFactor(kIndex) = 1.0d0
      end if
      rmsvalue(kIndex) = scaleFactor(kIndex) * rmsvalue(kIndex)
    end do

    ! write to file
    if (mmpi_myid == 0) then
      write(*,*) 'epp_printRmsStats: Opening ascii output file: ', trim(fileName)
      nulFile = 0
      ierr = fnom (nulFile, fileName, 'SEQ+R/W', 0)
      if (ierr /= 0) then
        call utl_abort('epp_printRmsStats: Cannot open ascii output file')
      end if

      kIndexCount = 0
      if ( (nomvar_v(1) == 'UU') .and. (nomvar_v(1+nLev_M) == 'VV') ) then        
        do levIndex = 1, nLev_M
          kIndexCount = kIndexCount + 2
          kIndexUU = levIndex
          kIndexVV = levIndex + nLev_M
          write(nulFile,100) &
               elapsed,ftype,nEns,nomvar_v(kIndexUU),pressureOrDepth(kIndexUU),rmsvalue(kIndexUU)
          write(nulFile,100) &
               elapsed,ftype,nEns,nomvar_v(kIndexVV),pressureOrDepth(kIndexVV),rmsvalue(kIndexVV)
        end do
      end if
      do kIndex = kIndexCount+1, numK
        write(nulFile,100) &
             elapsed,ftype,nEns,nomvar_v(kIndex),pressureOrDepth(kIndex),rmsvalue(kIndex)   
      end do
      ierr = fclos(nulFile)
    end if

100 format(f7.2,1x,A1,1x,I5,1x,A4,1x,f12.7,1x,(2E12.5))

    deallocate(weight)
    deallocate(nomvar_v)
    deallocate(pressureOrDepth)
    deallocate(scaleFactor)
    deallocate(rmsvalue)

  end subroutine epp_printRmsStats

end module ensPostProcess_mod
