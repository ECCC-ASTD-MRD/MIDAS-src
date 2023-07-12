
module scaleDecomposition_mod
  ! MODULE scaleDecomposition_mod (prefix='scd' category='4. Data Object transformations')
  !
  ! :Purpose: To perform horizontal and vertical scale decomposition of an ensemble state vector
  !
  use utilities_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use ensembleStatevector_mod
  use gridStateVector_mod
  use midasMpi_mod
  use globalSpectralTransform_mod
  use lamSpectralTransform_mod
  use verticalModes_mod
  use varNameList_mod
  use timeCoord_mod
  use earthConstants_mod
  use mathPhysConstants_mod
  implicit none
  save
  private

  ! Public subroutines
  public :: scd_horizontal, scd_vertical
  
contains

  !--------------------------------------------------------------------------
  ! scd_horizontal
  !-------------------------------------------------------------------------- 
  subroutine scd_horizontal(ensembleStateVector, nEnsOverDimension, nHorizWaveBand, horizWaveBandPeaks, &
                            decompositionMode, filterResponseFunctionMode, &
                            horizWaveBandIndexSelected_opt, writeResponseFunction_opt)
    !
    ! :Purpose: Perform a horizontal scale decomposition of an ensemble state vector
    !
    ! --> Mode Split (e.g. for SDL in bMatrixEnsemble_mod)
    !
    ! --- Ensemble Data at the Start  ---
    ! ensembleStateVector(1               ) contains the full data
    ! ensembleStateVector(2:nHorizWaveBand) already allocated but empty
    !
    ! --- Ensemble Data at the End    ---
    ! ensembleStateVector(nHorizWaveBand  ) contains the largest scales
    ! ...
    ! ensembleStateVector(1               ) contains the smallest scales
    !
    ! --> Mode Select (e.g. for SDLwSL in bMatrixEnsemble_mod)
    !
    ! --- Ensemble Data at the Start  ---
    ! ensembleStateVector(1) contains the full data
    !
    ! --- Ensemble Data at the End    ---
    ! ensembleStateVector(1) contains the selected scales
    !
    implicit none

    ! Arguments:
    type(struct_ens), intent(inout) :: ensembleStateVector(:)
    integer, intent(in) :: nEnsOverDimension
    integer, intent(in) :: nHorizWaveBand
    integer, intent(in) :: horizWaveBandPeaks(:)
    character(len=*), intent(in) :: decompositionMode
    character(len=*), intent(in) :: filterResponseFunctionMode
    integer, optional, intent(in) :: horizWaveBandIndexSelected_opt
    logical, optional, intent(in) :: writeResponseFunction_opt
    
    ! Locals:
    type(struct_hco), pointer :: hco
    type(struct_vco), pointer :: vco
    integer :: nEns, nTrunc, horizWaveBandIndexSelected
    integer :: horizWaveBandIndex, memberindex, stepIndex, levIndex, latIndex, lonIndex
    integer :: ila_filter, p, nla_filter, nphase_filter
    real(8), allocatable :: ResponseFunction(:,:)
    real(8), allocatable :: bandSum(:,:)
    real(8) :: totwvnb_r8, temp_r8, waveLength
    real(8), allocatable :: ensPertSP(:,:,:)
    real(8), allocatable :: ensPertSPfiltered(:,:,:)
    real(8), allocatable :: ensPertGD(:,:,:)
    real(4), pointer     :: ptr4d_r4(:,:,:,:)
    integer, allocatable :: nIndex_vec(:)
    integer :: totwvnb, totwvnbMax
    integer :: myLonBeg, myLonEnd, myLatBeg, myLatEnd
    integer :: lonPerPE, latPerPE, lonPerPEmax, latPerPEmax
    integer :: gstFilterID, mIndex, nIndex, mymBeg, mymEnd, mynBeg, mynEnd, mymSkip, mynSkip
    integer :: mymCount, mynCount
    integer :: horizWaveBandIndexStart, horizWaveBandIndexEnd
    integer :: horizWaveBandIndexLoopStart, horizWaveBandIndexLoopEnd, horizWaveBandIndexLoopDirection
    type(struct_lst)   :: lst_ben_filter ! Spectral transform Parameters for filtering
    character(len=128) :: outfilename
    character(len=2)   :: wbnum
    character(len=19)  :: kind
    logical :: writeResponseFunction

    if ( mmpi_myid == 0 ) then
      write(*,*)
      write(*,*) 'Horizontal scale decomposition of ensemble'
      write(*,*) '   number of WaveBands for filtering = ', nHorizWaveBand
      write(*,*) '   WaveBand Peaks (total wavenumber)...'
      do horizWaveBandIndex = 1, nHorizWaveBand
        write(*,*) horizWaveBandIndex, horizWaveBandPeaks(horizWaveBandIndex)
      end do
    end if

    if (present(writeResponseFunction_opt)) then
      writeResponseFunction = writeResponseFunction_opt
    else
      writeResponseFunction = .false.
    end if
    
    nEns = ens_getNumMembers(ensembleStateVector(1))
    
    hco => ens_getHco(ensembleStateVector(1))
    vco => ens_getVco(ensembleStateVector(1))

    call mmpi_setup_latbands(hco%nj, latPerPE, latPerPEmax, myLatBeg, myLatEnd)
    call mmpi_setup_lonbands(hco%ni, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd)

    if (present(horizWaveBandIndexSelected_opt)) then
      horizWaveBandIndexSelected=horizWaveBandIndexSelected_opt
    else
      horizWaveBandIndexSelected=-1
    end if
    
    !
    !- Setup a spectral transform for filtering (nk = nEnsOverDimension)
    !
    if (trim(decompositionMode) == 'Split') then
      nTrunc = horizWaveBandPeaks(1)
    else if (trim(decompositionMode) == 'Select') then
      if (horizWaveBandIndexSelected == 1) then
        nTrunc = -1 ! no truncation needed to extract the smallest scales
      else if (horizWaveBandIndexSelected == -1) then
        call utl_abort('scd_horizontal: horizWaveBandIndexSelected_opt must be provided in Select mode')
      else
        nTrunc = horizWaveBandPeaks(horizWaveBandIndexSelected-1)
      end if
    else
      write(*,*)
      write(*,*) 'decomposition mode = ', trim(decompositionMode)
      call utl_abort('scd_horizontal: unknown decomposition mode')
    end if

    if (trim(filterResponseFunctionMode) /= 'SumToOne'      .and. &
        trim(filterResponseFunctionMode) /= 'SquareSumToOne') then
      write(*,*)
      write(*,*) 'filter response function mode = ', trim(filterResponseFunctionMode)
      call utl_abort('scd_horizontal: unknown filter response function mode')
    end if
    
    if (hco%global) then
      ! Global mode
      gstFilterID = gst_setup(hco%ni,hco%nj,nTrunc,nEnsOverDimension)
      if (mmpi_myid == 0) write(*,*) 'ben : returned value of gstFilterID = ',gstFilterID

      nla_filter = gst_getNla(gstFilterID)
      nphase_filter = 2

      allocate(nIndex_vec(nla_filter))
      call mmpi_setup_m(gst_getNtrunc(gstFilterID),mymBeg,mymEnd,mymSkip,mymCount)
      call mmpi_setup_n(gst_getNtrunc(gstFilterID),mynBeg,mynEnd,mynSkip,mynCount)
      ila_filter = 0
      do mIndex = mymBeg, mymEnd, mymSkip
        do nIndex = mynBeg, mynEnd, mynSkip
          if (mIndex.le.nIndex) then
            ila_filter = ila_filter + 1
            nIndex_vec(ila_filter) = nIndex
          end if
        end do
      end do

    else
      ! LAM mode
      call lst_Setup(lst_ben_filter,                                                      & ! OUT
                     hco%ni, hco%nj, hco%dlon, nTrunc,                                    & ! IN
                     'LatLonMN', maxlevels_opt=nEnsOverDimension, gridDataOrder_opt='kij')  ! IN

      nla_filter    = lst_ben_filter%nla
      nphase_filter = lst_ben_filter%nphase
    end if

    !
    !- 1.  Scale decomposition for every wave band except for wave band #1
    !
    if (trim(decompositionMode) == 'Split') then
      horizWaveBandIndexStart     = 2 ! Skip the smallest scales
      horizWaveBandIndexEnd       = nHorizWaveBand
      horizWaveBandIndexLoopStart = horizWaveBandIndexEnd ! Start with the largest scales
      horizWaveBandIndexLoopEnd   = horizWaveBandIndexStart
      horizWaveBandIndexLoopDirection = -1 
    else ! Select
      horizWaveBandIndexStart     = horizWaveBandIndexSelected
      horizWaveBandIndexEnd       = horizWaveBandIndexSelected
      horizWaveBandIndexLoopStart = horizWaveBandIndexStart
      horizWaveBandIndexLoopEnd   = horizWaveBandIndexEnd 
      horizWaveBandIndexLoopDirection = 1 
    end if

    allocate(ResponseFunction(nla_filter,horizWaveBandIndexStart:horizWaveBandIndexEnd))
    allocate(ensPertSP(nla_filter,nphase_filter,nEnsOverDimension))
    allocate(ensPertSPfiltered(nla_filter,nphase_filter,nEnsOverDimension))
    allocate(ensPertGD(nEnsOverDimension,myLonBeg:myLonEnd,myLatBeg:myLatEnd))

    ensPertSP        (:,:,:) = 0.0d0
    ensPertSPfiltered(:,:,:) = 0.0d0

    !- 1.1 Pre-compute the response function
    do horizWaveBandIndex = horizWaveBandIndexLoopStart, horizWaveBandIndexLoopEnd, horizWaveBandIndexLoopDirection

      totwvnbMax=-1
      do ila_filter = 1, nla_filter
        if (hco%global) then
          totwvnb_r8 = real(nIndex_vec(ila_filter),8)
        else
          totwvnb_r8 = lst_ben_filter%k_r8(ila_filter)
        end if
        totwvnbMax=max(totwvnbMax,nint(totwvnb_r8))

        responseFunction(ila_filter,horizWaveBandIndex) = &
             scd_filterResponseFunction(totwvnb_r8,horizWaveBandIndex, horizWaveBandPeaks, &
                                        nHorizWaveBand)
        if (trim(filterResponseFunctionMode) == 'SquareSumToOne') then
          responseFunction(ila_filter,horizWaveBandIndex) = sqrt(responseFunction(ila_filter,horizWaveBandIndex)) 
        end if

        write(*,*) totwvnb_r8, ResponseFunction(ila_filter,horizWaveBandIndex)
      end do

      ! For visualization only: Compute and output the response function for each integer total wave number
      if (writeResponseFunction .and. mmpi_myid == 0) then
        
        write(wbnum,'(I2.2)') horizWaveBandIndex
        outfilename = "./ResponseFunction_"//wbnum//".txt"
        open (unit=99,file=outfilename,action="write",status="new")

        do totwvnb = 0, totwvnbMax
          temp_r8 = scd_filterResponseFunction(real(totwvnb,8),horizWaveBandIndex, horizWaveBandPeaks, nHorizWaveBand)
          if (trim(filterResponseFunctionMode) == 'SquareSumToOne') then
            temp_r8 = sqrt(temp_r8)
          end if
          if (hco%global) then
            if (totwvnb /= 0) then
              waveLength=2.d0*MPC_PI_R8*ec_ra/real(totwvnb,8)
            else
              waveLength=0.d0
            end if
            write(99,'(I4,2X,F7.1,2X,F5.3)') totwvnb, waveLength/1000.d0, temp_r8
          else
            write(99,'(I4,2X,F5.3)') totwvnb, temp_r8
          end if
        end do
        
        close(unit=99)
      end if

    end do ! horizWaveBandIndex

    if (hco%global) deallocate(nIndex_vec)

    do stepIndex = 1, ens_getNumStep(ensembleStateVector(1)) ! Loop on ensemble time bin
      do levIndex = 1, ens_getNumK(ensembleStateVector(1)) ! Loop on variables and vertical levels

        ptr4d_r4 => ens_getOneLev_r4(ensembleStateVector(1),levIndex)
          
        !- 1.2 GridPoint space -> Spectral Space
        !$OMP PARALLEL DO PRIVATE (latIndex)
        do latIndex = myLatBeg, myLatEnd
          ensPertGD(:,:,latIndex) = 0.0d0
        end do
        !$OMP END PARALLEL DO
        !$OMP PARALLEL DO PRIVATE (memberIndex,latIndex,lonIndex)
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            do memberIndex = 1, nEns
              ensPertGD(memberIndex,lonIndex,latIndex) = dble(ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex))
            end do
          end do
        end do
        !$OMP END PARALLEL DO
        if (hco%global) then
          ! Global Mode
          call gst_setID(gstFilterID) ! IN
          call gst_reespe_kij(ensPertSP, & ! OUT
                              ensPertGD)   ! IN
        else
          ! LAM mode
          kind = 'GridPointToSpectral'
          call lst_VarTransform(lst_ben_filter,         & ! IN
                                ensPertSP,              & ! OUT
                                ensPertGD,              & ! IN 
                                kind, nEnsOverDimension ) ! IN
        end if

        !- 1.3 Filtering and transformation back to grid point space 
        do horizWaveBandIndex = horizWaveBandIndexLoopStart, horizWaveBandIndexLoopEnd, horizWaveBandIndexLoopDirection

          ! Filtering
          !$OMP PARALLEL DO PRIVATE (memberIndex,p,ila_filter)
          do memberIndex = 1, nEns
            do p = 1, nphase_filter
              do ila_filter = 1, nla_filter
                ensPertSPfiltered(ila_filter,p,memberIndex) = &
                     ensPertSP(ila_filter,p,memberIndex) * responseFunction(ila_filter,horizWaveBandIndex)
              end do
            end do
          end do
          !$OMP END PARALLEL DO

          ! Spectral Space -> GridPoint space
          if (hco%global) then
            ! Global Mode
            call gst_setID(gstFilterID) ! IN
            call gst_speree_kij(ensPertSPfiltered, & ! IN
                                ensPertGD)           ! OUT
          else
            ! LAM mode
            kind = 'SpectralToGridPoint'
            call lst_VarTransform(lst_ben_filter,         & ! IN
                                  ensPertSPfiltered,      & ! IN
                                  ensPertGD,              & ! OUT
                                  kind, nEnsOverDimension)  ! IN
          end if

          if (trim(decompositionMode) == 'Split') then
            ptr4d_r4 => ens_getOneLev_r4(ensembleStateVector(horizWaveBandIndex),levIndex)
          else
            ptr4d_r4 => ens_getOneLev_r4(ensembleStateVector(1),levIndex)
          end if
          !$OMP PARALLEL DO PRIVATE (memberIndex,latIndex,lonIndex)
          do latIndex = myLatBeg, myLatEnd
            do lonIndex = myLonBeg, myLonEnd
              do memberIndex = 1, nEns
                ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex) = sngl(ensPertGD(memberIndex,lonIndex,latIndex))
              end do
            end do
          end do
          !$OMP END PARALLEL DO

        end do ! horizWaveBandIndex
      end do ! time bins
    end do ! variables&levels

    deallocate(ensPertGD)
    deallocate(responseFunction)
    deallocate(ensPertSP)
    deallocate(ensPertSPfiltered)

    !
    !- 2.  In split mode, isolate the smallest scales in horizWaveBandIndex = 1 by difference in grid point space
    !
    if (trim(decompositionMode) == 'Split') then

      allocate(bandSum(myLonBeg:myLonEnd,myLatBeg:myLatEnd))
      do stepIndex = 1, ens_getNumStep(ensembleStateVector(1))
        !$OMP PARALLEL DO PRIVATE (memberIndex,levIndex,latIndex,lonIndex,horizWaveBandIndex,bandsum,ptr4d_r4)
        do levIndex = 1, ens_getNumK(ensembleStateVector(1))
          do memberIndex = 1, nEns
            bandSum(:,:) = 0.d0
            do horizWaveBandIndex = 2, nHorizWaveBand
              ptr4d_r4 => ens_getOneLev_r4(ensembleStateVector(horizWaveBandIndex),levIndex)
              do latIndex = myLatBeg, myLatEnd
                do lonIndex = myLonBeg, myLonEnd
                  bandSum(lonIndex,latIndex) = bandSum(lonIndex,latIndex) + dble(ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex))
                end do
              end do
            end do
            ptr4d_r4 => ens_getOneLev_r4(ensembleStateVector(1),levIndex)
            do latIndex = myLatBeg, myLatEnd
              do lonIndex = myLonBeg, myLonEnd
                ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex) = sngl(dble(ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex)) - bandSum(lonIndex,latIndex))
              end do
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do
      deallocate(bandSum)

    end if

  end subroutine scd_horizontal

  !--------------------------------------------------------------------------
  ! scd_vertical
  !--------------------------------------------------------------------------
  subroutine scd_vertical(ensembleStateVector, nVertWaveBand,                &
                          vertWaveBandPeaks, vertmodesLengthScale,           &
                          decompositionMode, vertWaveBandIndexSelected_opt,  &
                          writeResponseFunction_opt, writeTransformInfo_opt, &
                          TGhandling_opt)
    !
    ! :Purpose: Perform a vertical scale decomposition of an ensemble state vector
    !
    ! --> Mode Split (e.g. for SDL in bMatrixEnsemble_mod)
    !
    ! --- Ensemble Data at the Start  ---
    ! ensembleStateVector(1               ) contains all the vertical scales
    ! ensembleStateVector(2:nVertWaveBand) already allocated but empty
    !
    ! --- Ensemble Data at the End    ---
    ! ensembleStateVector(nVertWaveBand  ) contains the deepest scales
    ! ...
    ! ensembleStateVector(1               ) contains the shallowest scales
    !
    ! --> Mode Select (e.g. for calcStatsGlb_mod)
    !
    ! --- Ensemble Data at the Start  ---
    ! ensembleStateVector(1) contains all the vertical scales
    !
    ! --- Ensemble Data at the End    ---
    ! ensembleStateVector(1) contains the selected scales
    !
    implicit none

    ! Arguments:
    type(struct_ens), intent(inout) :: ensembleStateVector(:)
    integer, intent(in) :: nVertWaveBand
    integer, intent(in) :: vertWaveBandPeaks(:)
    real(8), intent(in) :: vertModesLengthScale(2)
    character(len=*), intent(in) :: decompositionMode
    integer, optional, intent(in) :: vertWaveBandIndexSelected_opt
    logical, optional, intent(in) :: writeResponseFunction_opt
    logical, optional, intent(in) :: writeTransformInfo_opt
    character(len=*), optional, intent(in) :: TGhandling_opt

    ! Locals:
    type(struct_gsv), allocatable :: gridStateVector_oneMember(:)
    type(struct_hco), pointer :: hco
    real(8), allocatable :: vertModesState4d(:,:,:,:)
    real(8), allocatable :: vertModesState4dFiltered(:,:,:,:)
    real(8), allocatable, target :: gridState4d(:,:,:,:)
    real(8), pointer     :: ptr4d_r8(:,:,:,:)
    real(8), allocatable :: responseFunction(:,:)
    real(8), allocatable :: bandSum(:,:)
    real(4), pointer     :: ptr4d_r4(:,:,:,:)
    character(len=4), pointer :: varNamesList(:)
    character(len=128) :: outfilename
    character(len=2)   :: wbnum
    character(len=4)   :: varNameForTransform
    character(len=64)  :: TGhandling    
    integer :: vertWaveBandIndexSelected
    integer :: nMode, nModeMax, modeIndex, modeBeg, modeEnd, mTrunc
    integer :: myLonBeg, myLonEnd, myLatBeg, myLatEnd
    integer :: lonPerPE, latPerPE, lonPerPEmax, latPerPEmax
    integer :: numVar, varIndex, memberIndex
    integer :: numStep, stepIndex
    integer :: nLev, nLev_M, nLev_T
    integer :: levIndex,latIndex,lonIndex
    integer :: vertWaveBandIndex, vertWaveBandIndexStart, vertWaveBandIndexEnd
    integer :: vertWaveBandIndexLoopStart, vertWaveBandIndexLoopEnd, vertWaveBandIndexLoopDirection
    logical :: writeResponseFunction
    logical :: writeTransformInfo
    type(struct_vms), save :: vModes
    logical,          save :: firstTime = .true.

    if ( mmpi_myid == 0 ) then
      write(*,*)
      write(*,*) 'Vertical scale decomposition of ensemble'
      write(*,*) '   number of WaveBands for filtering = ', nVertWaveBand
      write(*,*) '   WaveBand Peaks (total wavenumber)...'
      do vertWaveBandIndex = 1, nVertWaveBand
        write(*,*) vertWaveBandIndex, vertWaveBandPeaks(vertWaveBandIndex)
      end do
    end if
    
    !
    !- 1.  Setup some variables
    !
    if (present(vertWaveBandIndexSelected_opt)) then
      vertWaveBandIndexSelected = vertWaveBandIndexSelected_opt
    else
      vertWaveBandIndexSelected = -1
    end if

    if (present(TGhandling_opt)) then
      TGhandling=TGhandling_opt
    else
      TGhandling='fullPertsInMediumScales'
    end if

    if (trim(decompositionMode) == 'Split') then
      mTrunc = vertWaveBandPeaks(1)
      vertWaveBandIndexStart     = 2 ! Skip the shallowest scales
      vertWaveBandIndexEnd       = nVertWaveBand
      vertWaveBandIndexLoopStart = vertWaveBandIndexEnd ! Start with the deepest scales
      vertWaveBandIndexLoopEnd   = vertWaveBandIndexStart
      vertWaveBandIndexLoopDirection = -1
    else if (trim(decompositionMode) == 'Select') then
      if (vertWaveBandIndexSelected == -1) then
        call utl_abort('scd_vertical: vertWaveBandIndexSelected_opt must be provided in Select mode')
      else if (vertWaveBandIndexSelected == 1) then
        mTrunc = -1 ! no truncation needed to extract the shallowest scales
      else
        mTrunc = vertWaveBandPeaks(vertWaveBandIndexSelected-1)
      end if
      
      vertWaveBandIndexStart     = vertWaveBandIndexSelected
      vertWaveBandIndexEnd       = vertWaveBandIndexSelected
      vertWaveBandIndexLoopStart = vertWaveBandIndexStart
      vertWaveBandIndexLoopEnd   = vertWaveBandIndexEnd 
      vertWaveBandIndexLoopDirection = 1 
    else
      write(*,*)
      write(*,*) 'decomposition mode = ', trim(decompositionMode)
      call utl_abort('scd_vertical: unknown decomposition mode')
    end if

    nLev_M = ens_getNumLev(ensembleStateVector(1),'MM')
    nLev_T = ens_getNumLev(ensembleStateVector(1),'TH')
    nModeMax=max(nLev_M,nLev_T)

    numStep = ens_getNumStep(ensembleStateVector(1))

    hco => ens_getHco(ensembleStateVector(1))
    call mmpi_setup_latbands(hco%nj, latPerPE, latPerPEmax, myLatBeg, myLatEnd)
    call mmpi_setup_lonbands(hco%ni, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd)
    
    !- Pre-compute the response function
    if (present(writeResponseFunction_opt)) then
      writeResponseFunction = writeResponseFunction_opt
    else
      writeResponseFunction = .false.
    end if
    if (present(writeTransformInfo_opt)) then
      writeTransformInfo = writeTransformInfo_opt
    else
      writeTransformInfo = .false.
    end if

    allocate(responseFunction(nModeMax,vertWaveBandIndexStart:vertWaveBandIndexEnd))

    do vertWaveBandIndex = vertWaveBandIndexLoopStart, vertWaveBandIndexLoopEnd, vertWaveBandIndexLoopDirection
      
      if (writeResponseFunction .and. mmpi_myid == 0) then
        write(wbnum,'(I2.2)') vertWaveBandIndex
        outfilename = "./vertResponseFunction_"//wbnum//".txt"
        open (unit=99,file=outfilename,action="write",status="new")
      end if
      
      do modeIndex = 1, nModeMax
        responseFunction(modeIndex,vertWaveBandIndex) = &
             scd_filterResponseFunction(dble(modeIndex), vertWaveBandIndex, &
                                        vertWaveBandPeaks, nVertWaveBand)
        if (mmpi_myid == 0) then
          write(* ,'(I4,2X,F5.3)') modeIndex, responseFunction(modeIndex,vertWaveBandIndex)
          if (writeResponseFunction) then
            write(99,'(I4,2X,F5.3)') modeIndex, responseFunction(modeIndex,vertWaveBandIndex)
          end if
        end if
      end do

      if (writeResponseFunction .and. mmpi_myid == 0) then
        close(unit=99)
      end if
      
    end do

    !
    !- 2. Setup transform
    !
    if (firstTime) then
      call vms_computeModesFromFunction(ens_getVco(ensembleStateVector(1)),               & ! IN
                                        vertModesLengthScale(1), vertModesLengthScale(2), & ! IN
                                        vModes)                                             ! OUT
      if (writeTransformInfo) then
        call vms_writeModes(vModes)
      end if
      firstTime = .false.
    end if

    !
    !- 3.  Vertical scale decomposition
    !
    nullify(varNamesList)
    call ens_varNamesList(varNamesList,ensembleStateVector(1))
    numVar = size(varNamesList)

    if (trim(decompositionMode) == 'Split') then
      allocate(gridStateVector_oneMember(nVertWaveBand))
      do vertWaveBandIndex = 1, nVertWaveBand
        call gsv_allocate(gridStateVector_oneMember(vertWaveBandIndex), numStep,                  &
                          ens_getHco(ensembleStateVector(1)), ens_getVco(ensembleStateVector(1)), &
                          varNames_opt=varNamesList, datestamp_opt=tim_getDatestamp(),            &
                          mpi_local_opt=.true., mpi_distribution_opt='Tiles', dataKind_opt=8)
      end do
    else ! Select
      allocate(gridStateVector_oneMember(1))
      call gsv_allocate(gridStateVector_oneMember(1), numStep,                                  &
                        ens_getHco(ensembleStateVector(1)), ens_getVco(ensembleStateVector(1)), &
                        varNames_opt=varNamesList, datestamp_opt=tim_getDatestamp(),            &
                        mpi_local_opt=.true., mpi_distribution_opt='Tiles', dataKind_opt=8)
    end if

    do memberIndex = 1, ens_getNumMembers(ensembleStateVector(1))
      call ens_copyMember(ensembleStateVector(1), gridStateVector_oneMember(1), memberIndex) ! extract all the vertical scales
      
      do varIndex = 1, numVar

        varNameForTransform = varNamesList(varIndex)
        nullify(ptr4d_r8)
        if (vnl_varLevelFromVarName(trim(varNamesList(varIndex))).eq.'MM') then
          nLev  = nLev_M
          call gsv_getField(gridStateVector_oneMember(1),ptr4d_r8,varName_opt=varNamesList(varIndex))
          
        else if (vnl_varLevelFromVarName(trim(varNamesList(varIndex))).eq.'TH') then
          nLev  = nLev_T
          call gsv_getField(gridStateVector_oneMember(1),ptr4d_r8,varName_opt=varNamesList(varIndex))
        else
          if (trim(varNamesList(varIndex)) == 'TG') then
            if (trim(TGhandling) == 'expandWithTT') then
              if (allocated(gridState4d)) deallocate(gridState4d)
              call expandTG(gridStateVector_oneMember(1), gridState4d, nLev)
              ptr4d_r8 => gridState4d
              varNameForTransform = 'TT'
            else
              cycle ! varIndex
            end if
          else if (trim(varNamesList(varIndex)) == 'P0') then
            ! Temporary measures
            call gsv_getField(gridStateVector_oneMember(1),ptr4d_r8,varName_opt=varNamesList(varIndex))
            ptr4d_r8(:,:,:,:) = 0.d0
            cycle ! varIndex
          else
            write(*,*)
            write(*,*) 'varname  = ', trim(varNamesList(varIndex))
            write(*,*) 'varlevel = ', vnl_varLevelFromVarName(trim(varNamesList(varIndex)))
            call utl_abort('scd_vertical: variable not handle yet')
          end if
        end if
        nMode = nLev
        
        if (allocated(vertModesState4d)) deallocate(vertModesState4d)
        allocate(vertModesState4d(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nMode,numStep))

        !- GridPoint space -> Vertical modes space
        do stepIndex = 1, numStep ! Loop on ensemble time bin
          call vms_transform(vModes, vertModesState4d(:,:,:,stepIndex),         &
                             ptr4d_r8(:,:,:,stepIndex), 'GridPointToVertModes', &
                             myLonBeg, myLonEnd, myLatBeg, myLatEnd, nLev,      &
                             varNameForTransform, modeEnd_opt=mTrunc)
        end do

        !- Filtering
        if (allocated(vertModesState4dFiltered)) deallocate(vertModesState4dFiltered)
        allocate(vertModesState4dFiltered(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nMode,numStep))
        
        do vertWaveBandIndex = vertWaveBandIndexLoopStart, vertWaveBandIndexLoopEnd, vertWaveBandIndexLoopDirection

          if (trim(decompositionMode) == 'Split' .and. trim(varNamesList(varIndex)) /= 'TG' ) then
            ! Note that this step will only be done if TGhandling = 'expandWithTT'
            nullify(ptr4d_r8)
            call gsv_getField(gridStateVector_oneMember(vertWaveBandIndex),ptr4d_r8,varName_opt=varNamesList(varIndex))
          end if

          ! Select only the modes needed to make the transform faster
          if (vertWaveBandIndex == nVertWaveBand ) then
            modeBeg = 1
            modeEnd = vertWaveBandPeaks(vertWaveBandIndex-1)
          else if ( vertWaveBandIndex /= 1 ) then
            modeBeg = vertWaveBandPeaks(vertWaveBandIndex+1)
            modeEnd = vertWaveBandPeaks(vertWaveBandIndex-1)
          else
            modeBeg = vertWaveBandPeaks(vertWaveBandIndex+1)
            modeEnd = nMode
          end if
          
          do stepIndex = 1, numStep ! Loop on ensemble time bin
            
            !$OMP PARALLEL DO PRIVATE (modeIndex,latIndex,lonIndex)
            do modeIndex = modeBeg, modeEnd
              do latIndex = myLatBeg, myLatEnd
                do lonIndex = myLonBeg, myLonEnd
                  vertModesState4dFiltered(lonIndex,latIndex,modeIndex,stepIndex) =   &
                       responseFunction(modeIndex,vertWaveBandIndex) *        &
                       vertModesState4d(lonIndex,latIndex,modeIndex,stepIndex)
                end do
              end do
            end do
            !$OMP END PARALLEL DO

            !- Vertical modes space -> GridPoint space
            call vms_transform(vModes, vertModesState4dFiltered(:,:,:,stepIndex), &
                               ptr4d_r8(:,:,:,stepIndex), 'VertModesToGridPoint', &
                               myLonBeg, myLonEnd, myLatBeg, myLatEnd, nLev,      &
                               varNameForTransform, modeBeg_opt=modeBeg,          &
                               modeEnd_opt=modeEnd)

          end do ! stepIndex

          if (trim(varNamesList(varIndex)) == 'TG') then
            ! Note that this step will only be done if TGhandling = 'expandWithTT'
            if (trim(decompositionMode) == 'Split') then
              call extractTG(gridStateVector_oneMember(vertWaveBandIndex), gridState4d, nLev)
            else
              call extractTG(gridStateVector_oneMember(1), gridState4d, nLev)
            end if
          end if
          
        end do ! vertWaveBandIndex
      end do ! variables

      if (ens_varExist(ensembleStateVector(1),'TG') .and. trim(TGhandling) /= 'expandWithTT') then
        call adhocTGdecomposition(gridStateVector_oneMember, decompositionMode, TGhandling, &
                                  vertWaveBandIndexLoopStart, vertWaveBandIndexLoopEnd, &
                                  vertWaveBandIndexLoopDirection, nVertWaveBand)
      end if

      if (trim(decompositionMode) == 'Split') then
        do vertWaveBandIndex = vertWaveBandIndexLoopStart, vertWaveBandIndexLoopEnd, vertWaveBandIndexLoopDirection
          call ens_insertMember(ensembleStateVector(vertWaveBandIndex),                  &
                                gridStateVector_oneMember(vertWaveBandIndex), memberIndex)
        end do
      else ! Select
        call ens_insertMember(ensembleStateVector(1), gridStateVector_oneMember(1), memberIndex)
      end if
 
    end do ! memberIndex

    deallocate(vertModesState4dFiltered)
    deallocate(vertModesState4d)
    if (trim(decompositionMode) == 'Split') then
      do vertWaveBandIndex = 1, nVertWaveBand
        call gsv_deallocate(gridStateVector_oneMember(vertWaveBandIndex))
      end do
    else
      call gsv_deallocate(gridStateVector_oneMember(1))
    end if
    deallocate(responseFunction)
    
    !
    !- 3.  In split mode, isolate the shallowest scales in vertWaveBandIndex = 1 by difference in grid point space
    !
    if (trim(decompositionMode) == 'Split') then

      allocate(bandSum(myLonBeg:myLonEnd,myLatBeg:myLatEnd))
      do stepIndex = 1, ens_getNumStep(ensembleStateVector(1))
        !$OMP PARALLEL DO PRIVATE (memberIndex,levIndex,latIndex,lonIndex,vertWaveBandIndex,bandsum,ptr4d_r4)
        do levIndex = 1, ens_getNumK(ensembleStateVector(1))
          do memberIndex = 1, ens_getNumMembers(ensembleStateVector(1))
            bandSum(:,:) = 0.d0
            do vertWaveBandIndex = 2, nVertWaveBand
              ptr4d_r4 => ens_getOneLev_r4(ensembleStateVector(vertWaveBandIndex),levIndex)
              do latIndex = myLatBeg, myLatEnd
                do lonIndex = myLonBeg, myLonEnd
                  bandSum(lonIndex,latIndex) = bandSum(lonIndex,latIndex) + dble(ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex))
                end do
              end do
            end do
            ptr4d_r4 => ens_getOneLev_r4(ensembleStateVector(1),levIndex)
            do latIndex = myLatBeg, myLatEnd
              do lonIndex = myLonBeg, myLonEnd
                ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex) = sngl(dble(ptr4d_r4(memberIndex,stepIndex,lonIndex,latIndex)) - bandSum(lonIndex,latIndex))
              end do
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end do
      deallocate(bandSum)

    end if
    
  end subroutine scd_vertical
  
  !--------------------------------------------------------------------------
  ! scd_filterResponseFunction
  !--------------------------------------------------------------------------
  function scd_filterResponseFunction(totalWaveNumber, waveBandIndex, waveBandPeaks, &
                                      nWaveBand) result(ResponseFunction) 
    implicit none

    ! Arguments:
    real(8), intent(in) :: totalWaveNumber
    integer, intent(in) :: waveBandIndex
    integer, intent(in) :: nWaveBand
    integer, intent(in) :: waveBandPeaks(:)
    ! Result:
    real(8) :: ResponseFunction 

    ! Locals:
    real(8) :: linearResponse, lowerLimit, center, upperLimit
    real(8), parameter :: pi = 2.d0*asin(1.d0)

    if (waveBandIndex == nWaveBand ) then
       ! This wave band contains the largest scales.
       !
       ! The response function is 1 total wave number <= waveBandPeaks(nWaveBand)
       ! and decreases to 0 at waveBandPeaks(nWaveBand-1)
       !
       !                    response=1 |---
       !                               |   \
       !                               |    \
       !                    response=0 |------------
       !       waveBandPeaks(nWaveBand) <-|  |-> waveBandPeaks(nWaveBand-1)
       !
       lowerlimit = real(waveBandPeaks(waveBandIndex  ),8)
       upperlimit = real(waveBandPeaks(waveBandIndex-1),8)

       if ( totalWaveNumber < lowerlimit ) then
          ResponseFunction = 1.d0
       else if ( totalWaveNumber <= upperlimit ) then
          linearResponse = (upperlimit-totalWaveNumber) / (upperlimit-lowerlimit)
          ResponseFunction = sin( (pi/2.d0) * linearResponse)**2
       else
          ResponseFunction = 0.d0
       end if

    else if ( waveBandIndex /= 1 ) then
       ! This wave band contains intermediate scales (i.e., not the largest or the smallest).
       !
       ! The response function is 1 (only) for the total wave number = waveBandPeaks(waveBandIndex)
       ! and decreases to 0 at both waveBandPeaks(waveBandIndex+1) and waveBandPeaks(waveBandIndex-1)
       !
       !                    response=1 |      -
       !                               |     / \
       !                               |    /   \
       !                    response=0 |------------
       !  waveBandPeaks(waveBandIndex+1) <-|     |-> waveBandPeaks(waveBandIndex-1)
       !                                      |-> waveBandPeaks(waveBandIndex)
       !
       center     = real(waveBandPeaks(waveBandIndex  ),8)
       upperlimit = real(waveBandPeaks(waveBandIndex-1),8)
       lowerlimit = real(waveBandPeaks(waveBandIndex+1),8)

       if (      totalWaveNumber >  lowerlimit .and. &
                 totalWaveNumber <= center     ) then
          linearResponse = (totalWaveNumber-lowerlimit) / (center-lowerlimit)
          ResponseFunction = sin( (pi/2.d0) * linearResponse)**2
       else if ( totalWaveNumber >  center      .and. &
                 totalWaveNumber <  upperlimit  ) then
          linearResponse = (upperlimit-totalWaveNumber) / (upperlimit-center)
          ResponseFunction = sin( (pi/2.d0) * linearResponse)**2
       else
          ResponseFunction = 0.d0
       end if

    else
       !
       ! This wave band contains the smallest scales.
       !
       ! The response function is 1 total wave number >= waveBandPeaks(nWaveBand)
       ! and decreases to 0 at waveBandPeaks(nWaveBand-1)
       !
       !                    response=1 |    ---
       !                               |   /
       !                               |  / 
       !                    response=0 |------------
       !waveBandPeaks(waveBandIndex+1) <-|  |-> waveBandPeaks(1)
       !
       upperlimit = real(waveBandPeaks(waveBandIndex  ),8)
       lowerlimit = real(waveBandPeaks(waveBandIndex+1),8)

       if      ( totalWaveNumber > upperlimit ) then
          ResponseFunction = 1.d0
       else if ( totalWaveNumber > lowerlimit ) then
          linearResponse = (totalWaveNumber-lowerlimit) / (upperlimit-lowerlimit)
          ResponseFunction = sin( (pi/2.d0) * linearResponse)**2
       else
          ResponseFunction = 0.d0
       end if

    end if

  end function scd_filterResponseFunction

  !--------------------------------------------------------------------------
  ! expandTG
  !--------------------------------------------------------------------------
   subroutine expandTG(statevector, gridState4d, nLev)
    !
    ! :Purpose: Combine TT and TG to create a new 4D field with TG at the lower boundary
    !
    implicit none

    ! Arguments:
    type(struct_gsv)    , intent(in)  :: statevector
    real(8), allocatable, intent(out) :: gridState4d(:,:,:,:)
    integer             , intent(out) :: nLev

    ! Locals:
    real(8), pointer     :: TTptr4d_r8(:,:,:,:)
    real(8), pointer     :: TGptr4d_r8(:,:,:,:)
    integer :: stepIndex, levIndex, lonIndex, latIndex

    call gsv_getField(stateVector,TTptr4d_r8,varName_opt='TT')
    call gsv_getField(stateVector,TGptr4d_r8,varName_opt='TG')

    nLev = gsv_getNumLevFromVarName(stateVector,'TT')

    allocate(gridState4d(statevector%myLonBeg:statevector%myLonEnd,statevector%myLatBeg:statevector%myLatEnd,nLev,statevector%numStep))

    !$OMP PARALLEL DO PRIVATE (stepIndex,levIndex,latIndex,lonIndex)
    do stepIndex = 1, statevector%numStep
      do levIndex = 1, nLev
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do lonIndex = statevector%myLonBeg, statevector%myLonEnd
            if (levIndex == nLev) then
              gridState4d(lonIndex,latIndex,levIndex,stepIndex) = TGptr4d_r8(lonIndex,latIndex,1,stepIndex)
            else
              gridState4d(lonIndex,latIndex,levIndex,stepIndex) = TTptr4d_r8(lonIndex,latIndex,levIndex,stepIndex)
            end if
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine expandTG

  !--------------------------------------------------------------------------
  ! extractTG
  !--------------------------------------------------------------------------
  subroutine extractTG(statevector, gridState4d, nLev)
    !
    ! :Purpose: Copy the TG field from gridState4d into the input stateVector
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: statevector
    integer,          intent(in)    :: nLev
    real(8),          intent(in)    :: gridState4d(statevector%myLonBeg:statevector%myLonEnd,statevector%myLatBeg:statevector%myLatEnd,nLev,statevector%numStep)

    ! Locals:
    real(8), pointer :: TGptr4d_r8(:,:,:,:)
    integer :: stepIndex, lonIndex, latIndex

    call gsv_getField(stateVector,TGptr4d_r8,varName_opt='TG')
    
    !$OMP PARALLEL DO PRIVATE (stepIndex,latIndex,lonIndex)
    do stepIndex = 1, stateVector%numStep
      do latIndex = stateVector%myLatBeg, stateVector%myLatEnd
        do lonIndex = stateVector%myLonBeg, stateVector%myLonEnd
          TGptr4d_r8(lonIndex,latIndex,1,stepIndex) = gridState4d(lonIndex,latIndex,nLev,stepIndex)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine extractTG

  !--------------------------------------------------------------------------
  ! adhocTGdecomposition
  !--------------------------------------------------------------------------
  subroutine adhocTGdecomposition(statevector, decompositionMode, TGhandling, &
                                  vertWaveBandIndexLoopStart, vertWaveBandIndexLoopEnd, &
                                  vertWaveBandIndexLoopDirection, nVertWaveBand)
    !
    ! :Purpose: Decompose TG with an ad hoc (i.e. crude) procedure
    !
    implicit none

    ! Arguments:
    type(struct_gsv), intent(inout) :: statevector(:)
    character(len=*), intent(in)    :: decompositionMode
    character(len=*), intent(in)    :: TGhandling
    integer,          intent(in)    :: vertWaveBandIndexLoopStart
    integer,          intent(in)    :: vertWaveBandIndexLoopEnd
    integer,          intent(in)    :: vertWaveBandIndexLoopDirection
    integer,          intent(in)    :: nVertWaveBand
    
    ! Locals:
    real(8), pointer :: TGfullPtr4d_r8(:,:,:,:)
    real(8), pointer :: TGdecompPtr4d_r8(:,:,:,:)
    real(8) :: factor
    integer :: vertWaveBandIndex, fullPertsInThisWaveBandIndex

    nullify(TGfullPtr4d_r8)
    call gsv_getField(stateVector(1),TGfullPtr4d_r8,varName_opt='TG')

    select case (trim(TGhandling))
    case ('fullPertsInShallowScales')
      fullPertsInThisWaveBandIndex = 1
    case ('fullPertsInMediumScales')
      fullPertsInThisWaveBandIndex = nint(real(nVertWaveBand,8) / 2.d0)
    case ('splitEvenly')
      factor = 1.d0 / real(nVertWaveBand,8)
    case default
      write(*,*)
      write(*,*) 'Error:  Unknown TGhandling ', trim(TGhandling)
      call utl_abort('scaleDecomposition : adhocTGdecomposition')
    end select

    if (trim(decompositionMode) == 'Split') then
      do vertWaveBandIndex = vertWaveBandIndexLoopStart, vertWaveBandIndexLoopEnd, vertWaveBandIndexLoopDirection
        nullify(TGdecompPtr4d_r8)
        call gsv_getField(stateVector(vertWaveBandIndex),TGdecompPtr4d_r8,varName_opt='TG')
        if (trim(TGhandling) == 'splitEvenly') then
          TGdecompPtr4d_r8(:,:,:,:) = TGfullPtr4d_r8(:,:,:,:) * factor
        else
          if (vertWaveBandIndex == fullPertsInThisWaveBandIndex) then
            TGdecompPtr4d_r8(:,:,:,:) = TGfullPtr4d_r8(:,:,:,:)
          else
            TGdecompPtr4d_r8(:,:,:,:) = 0.d0
          end if
        end if
      end do
    else ! select
      if (trim(TGhandling) == 'splitEvenly') then
          TGfullPtr4d_r8(:,:,:,:) = TGfullPtr4d_r8(:,:,:,:) / factor
        else
          if (vertWaveBandIndexLoopStart == fullPertsInThisWaveBandIndex) then
            TGfullPtr4d_r8(:,:,:,:) = TGfullPtr4d_r8(:,:,:,:)
          else
            TGfullPtr4d_r8(:,:,:,:) = 0.d0
          end if
        end if
    end if

  end subroutine adhocTGdecomposition
  
end module scaleDecomposition_mod
