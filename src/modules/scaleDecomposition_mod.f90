module scaleDecomposition_mod
  ! MODULE scaleDecomposition_mod (prefix='scd' category='4. Data Object transformations')
  !
  ! :Purpose: To
  !
  !use earthConstants_mod
  use utilities_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use ensembleStatevector_mod
  use midasMpi_mod
  use globalSpectralTransform_mod
  use lamSpectralTransform_mod
  implicit none
  save
  private

  ! Public subroutines
  public :: scd_horizontal !, scd_vertical
  public :: scd_filterResponseFunction
  
contains

  !--------------------------------------------------------------------------
  ! scd_horizontal
  !-------------------------------------------------------------------------- 
  subroutine scd_horizontal(ensembleStateVector, nEnsOverDimension, nHorizWaveBand, horizWaveBandPeaks, &
                            decompositionMode, filterResponseFunctionMode, &
                            horizWaveBandIndexSelected_opt)
    implicit none

    ! Arguments:
    type(struct_ens), intent(inout) :: ensembleStateVector(:)
    integer, intent(in) :: nEnsOverDimension
    integer, intent(in) :: nHorizWaveBand
    integer, intent(in) :: horizWaveBandPeaks(:)
    character(len=*), intent(in) :: decompositionMode
    character(len=*), intent(in) :: filterResponseFunctionMode
    integer, optional, intent(in) :: horizWaveBandIndexSelected_opt
    
    ! Locals:
    type(struct_hco), pointer :: hco
    type(struct_vco), pointer :: vco
    integer :: nEns, nTrunc, horizWaveBandIndexSelected
    integer :: horizWaveBandIndex, memberindex, stepIndex, levIndex, latIndex, lonIndex
    integer :: ila_filter, p, nla_filter, nphase_filter
    real(8), allocatable :: ResponseFunction(:,:)
    real(8), allocatable :: bandSum(:,:)
    real(8) :: totwvnb_r8
    real(8), allocatable :: ensPertSP(:,:,:)
    real(8), allocatable :: ensPertSPfiltered(:,:,:)
    real(8), allocatable :: ensPertGD(:,:,:)
    real(4), pointer     :: ptr4d_r4(:,:,:,:)
    integer, allocatable :: nIndex_vec(:)
    integer :: myLonBeg, myLonEnd, myLatBeg, myLatEnd
    integer        :: lonPerPE, latPerPE, lonPerPEmax, latPerPEmax
    integer :: gstFilterID, mIndex, nIndex, mymBeg, mymEnd, mynBeg, mynEnd, mymSkip, mynSkip
    integer :: mymCount, mynCount
    integer :: horizWaveBandIndexStart, horizWaveBandIndexEnd
    integer :: horizWaveBandIndexLoopStart, horizWaveBandIndexLoopEnd, horizWaveBandIndexLoopDirection
    type(struct_lst)    :: lst_ben_filter ! Spectral transform Parameters for filtering
    character(len=19)   :: kind

    !
    ! --> Mode Split (e.g. for SDL)
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

    !
    ! --> Mode Select (e.g. for SDLwSL)
    !
    ! --- Ensemble Data at the Start  ---
    ! ensembleStateVector(1) contains the full data
    !
    ! --- Ensemble Data at the End    ---
    ! ensembleStateVector(1) contains the selected scales
    !

    if ( mmpi_myid == 0 ) then
      write(*,*)
      write(*,*) 'Scale decomposition of the ensemble perturbations'
      write(*,*) '   number of WaveBands for filtering = ', nHorizWaveBand
      write(*,*) '   WaveBand Peaks (total wavenumber)...'
      do horizWaveBandIndex = 1, nHorizWaveBand
        write(*,*) horizWaveBandIndex, horizWaveBandPeaks(horizWaveBandIndex)
      end do
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
      call mmpi_setup_m(nTrunc,mymBeg,mymEnd,mymSkip,mymCount)
      call mmpi_setup_n(nTrunc,mynBeg,mynEnd,mynSkip,mynCount)
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
      do ila_filter = 1, nla_filter
        if (hco%global) then
          totwvnb_r8 = real(nIndex_vec(ila_filter),8)
        else
          totwvnb_r8 = lst_ben_filter%k_r8(ila_filter)
        end if
        ResponseFunction(ila_filter,horizWaveBandIndex) = &
             scd_filterResponseFunction(totwvnb_r8,horizWaveBandIndex, horizWaveBandPeaks, &
                                        nHorizWaveBand)
        if (trim(filterResponseFunctionMode) == 'SquareSumToOne') then
          responseFunction(ila_filter,horizWaveBandIndex) = sqrt(responseFunction(ila_filter,horizWaveBandIndex)) 
        end if
        write(*,*) totwvnb_r8, ResponseFunction(ila_filter,horizWaveBandIndex)
      end do
    end do
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
