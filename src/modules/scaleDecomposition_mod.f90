
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
                          writeResponseFunction_opt, writeTransformInfo_opt)
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

    ! Locals:
    type(struct_gsv), allocatable :: gridStateVector_oneMember(:)
    type(struct_hco), pointer :: hco
    real(8), allocatable :: vertModesState4d(:,:,:,:)
    real(8), allocatable :: vertModesState4dFiltered(:,:,:,:)
    real(8), pointer     :: gridState4d(:,:,:,:)
    real(8), allocatable :: responseFunction(:,:)
    character(len=128) :: outfilename
    character(len=2)   :: wbnum
    character(len=4), pointer :: varNamesList(:)
    integer :: vertWaveBandIndexSelected
    integer :: nMode, nModeMax, modeIndex, latIndex, lonIndex
    integer :: myLonBeg, myLonEnd, myLatBeg, myLatEnd
    integer :: lonPerPE, latPerPE, lonPerPEmax, latPerPEmax
    integer :: numVar, varIndex, memberIndex
    integer :: numStep, stepIndex
    integer :: nLev, nLev_M, nLev_T
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
    
    if (trim(decompositionMode) == 'Split') then
      vertWaveBandIndexStart     = 1
      vertWaveBandIndexEnd       = nVertWaveBand
      vertWaveBandIndexLoopStart = vertWaveBandIndexEnd ! Start with the largest scales
      vertWaveBandIndexLoopEnd   = vertWaveBandIndexStart
      vertWaveBandIndexLoopDirection = -1
    else if (trim(decompositionMode) == 'Select') then
      if (vertWaveBandIndexSelected == -1) then
        call utl_abort('scd_vertical: vertWaveBandIndexSelected_opt must be provided in Select mode')
      end if
      vertWaveBandIndexStart     = vertWaveBandIndexSelected
      vertWaveBandIndexEnd       = vertWaveBandIndexSelected
      vertWaveBandIndexLoopStart = vertWaveBandIndexStart
      vertWaveBandIndexLoopEnd   = vertWaveBandIndexEnd 
      vertWaveBandIndexLoopDirection = 1 
    else
      write(*,*)
      write(*,*) 'decomposition mode = ', trim(decompositionMode)
      call utl_abort('scd_vertocal: unknown decomposition mode')
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
      call gsv_allocate(gridStateVector_oneMember(1), numStep,                                    &
                        ens_getHco(ensembleStateVector(1)), ens_getVco(ensembleStateVector(1)), &
                        varNames_opt=varNamesList, datestamp_opt=tim_getDatestamp(),            &
                        mpi_local_opt=.true., mpi_distribution_opt='Tiles', dataKind_opt=8)
    end if

    do memberIndex = 1, ens_getNumMembers(ensembleStateVector(1))
      call ens_copyMember(ensembleStateVector(1), gridStateVector_oneMember(1), memberIndex) ! extract all the vertical scales

      do varIndex = 1, numVar

        nullify(gridState4d)
        call gsv_getField(gridStateVector_oneMember(1),gridState4d,varName_opt=varNamesList(varIndex))
        
        if (vnl_varLevelFromVarName(trim(varNamesList(varIndex))).eq.'MM') then
          nLev  = nLev_M
          nMode = nLev
        else if (vnl_varLevelFromVarName(trim(varNamesList(varIndex))).eq.'TH') then
          nLev  = nLev_T
          nMode = nLev
        else
          gridState4d(:,:,:,:) = 0.d0
          cycle
        end if

        if (allocated(vertModesState4d)) deallocate(vertModesState4d)
        allocate(vertModesState4d(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nMode,numStep))

        !- GridPoint space -> Vertical modes space
        do stepIndex = 1, numStep ! Loop on ensemble time bin
          call vms_transform(vModes, vertModesState4d(:,:,:,stepIndex),            &
                             gridState4d(:,:,:,stepIndex), 'GridPointToVertModes', &
                             myLonBeg, myLonEnd, myLatBeg, myLatEnd, nLev,         &
                             varNamesList(varIndex))
        end do

        !- Filtering
        if (allocated(vertModesState4dFiltered)) deallocate(vertModesState4dFiltered)
        allocate(vertModesState4dFiltered(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nMode,numStep))
        
        do vertWaveBandIndex = vertWaveBandIndexLoopStart, vertWaveBandIndexLoopEnd, vertWaveBandIndexLoopDirection

          if (trim(decompositionMode) == 'Split') then
            nullify(gridState4d)
            call gsv_getField(gridStateVector_oneMember(vertWaveBandIndex),gridState4d,varName_opt=varNamesList(varIndex))
          end if

          do stepIndex = 1, numStep ! Loop on ensemble time bin
            
            !$OMP PARALLEL DO PRIVATE (modeIndex,latIndex,lonIndex)
            do modeIndex = 1, nMode
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
            call vms_transform(vModes, vertModesState4dFiltered(:,:,:,stepIndex),    &
                               gridState4d(:,:,:,stepIndex), 'VertModesToGridPoint', &
                               myLonBeg, myLonEnd, myLatBeg, myLatEnd, nLev,         &
                               varNamesList(varIndex))

          end do ! stepIndex
        end do ! vertWaveBandIndex
      end do ! variables

      if (trim(decompositionMode) == 'Split') then
        do vertWaveBandIndex = vertWaveBandIndexLoopStart, vertWaveBandIndexLoopEnd, vertWaveBandIndexLoopDirection
          call ens_insertMember(ensembleStateVector(vertWaveBandIndex),                  &
                                gridStateVector_oneMember(vertWaveBandIndex), memberIndex)
        end do
      else ! Select
        call ens_insertMember(ensembleStateVector(1), gridStateVector_oneMember(1), memberIndex)
      end if
 
    end do ! memberIndex

    !
    !- 3.  Ending
    !
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
  
end module scaleDecomposition_mod
