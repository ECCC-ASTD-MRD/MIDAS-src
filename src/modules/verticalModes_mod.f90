
module verticalModes_mod
  ! MODULE verticalModes_mod (prefix='vms' category='4. Data Object transformations')
  !
  !:Purpose:  To 1) compute empirical orthogonal functions (EOFs) from either ensemble-derived
  !           vertical background-error covariances matrices or a prescribed vertical correlation
  !           function (i.e., the so-called vertical modes) and to 2) project back or forth
  !           ensemble pertubations onto these modes.
  !           Therefore, capablity #2 behaves like a spectral transform but in the vertical
  !           dimension.
  !
  use utilities_mod
  use varNameList_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use localizationFunction_mod
  use calcHeightAndPressure_mod
  use ensembleStatevector_mod
  use midasMpi_mod
  implicit none
  save
  private
  
  ! Public derived type
  public :: struct_vms
  ! Public subroutines
  public :: vms_computeModesFromEns, vms_computeModesFromFunction
  public :: vms_writeModes, vms_transform

  type :: struct_oneVar3d
    integer :: nLev
    character(len=4) :: varName
    real(8), allocatable :: autoCovariance(:,:)
    real(8), allocatable :: eigenVectors(:,:)
    real(8), allocatable :: eigenValues(:)
    real(8), allocatable :: eigenVectorsInv(:,:)
  end type struct_oneVar3d
  
  type :: struct_vms
    private
    integer :: nVar3d
    type(struct_oneVar3d), allocatable :: allVar3d(:)
    logical :: initialized = .false.
  end type struct_vms

contains

  !--------------------------------------------------------------------------
  ! vms_setup
  !--------------------------------------------------------------------------
  subroutine vms_setup(varNamesList, vco, vModes)
    !
    ! :Purpose: To setup the vModes structure
    !
    implicit none

    ! Arguments:
    character(len=4),          intent(in)    :: varNamesList(:)
    type(struct_vco), pointer, intent(in)    :: vco
    type(struct_vms),          intent(inout) :: vModes

    ! Locals:
    integer :: nLev, varNameIndex, var3dIndex
    
    if (vModes%initialized) return
    
    !
    !- Set parameters 
    !
    vModes%nVar3d = 0
    do varNameIndex = 1, size(varNamesList)
      if (vnl_varLevelFromVarname(varNamesList(varNameIndex)) /= 'SF' ) then
        vModes%nVar3d = vModes%nVar3d + 1
      end if
    end do

    allocate(vModes%allVar3d(vModes%nVar3d))

    var3dIndex = 0
    do varNameIndex = 1, size(varNamesList)
      if (vnl_varLevelFromVarname(varNamesList(varNameIndex)) /= 'SF' ) then
        var3dIndex = var3dIndex + 1
        vModes%allVar3d(var3dIndex)%varName = varNamesList(varNameIndex)
        if (vnl_varLevelFromVarName(vModes%allVar3d(var3dIndex)%varName) == 'MM') then
          nLev = vco_getNumLev(vco,'MM')
        else
          nLev = vco_getNumLev(vco,'TH')
        end if
        allocate(vModes%allVar3d(var3dIndex)%autoCovariance(nLev,nLev))
        allocate(vModes%allVar3d(var3dIndex)%eigenVectors(nLev,nLev))
        allocate(vModes%allVar3d(var3dIndex)%eigenValues(nLev))
        allocate(vModes%allVar3d(var3dIndex)%eigenVectorsInv(nLev,nLev))
        vModes%allVar3d(var3dIndex)%nLev = nLev
      end if
    end do

    !
    !- Ending
    !
    vModes%initialized = .true.

  end subroutine vms_setup
  
  !--------------------------------------------------------------------------
  ! vms_computeModesFromEns
  !--------------------------------------------------------------------------
  subroutine vms_computeModesFromEns(ensPerts,vModes)
    !
    ! :Purpose: To compute vertical modes from ensemble-derived
    !           vertical background-error covariances matrices
    !
    implicit none

    ! Arguments:
    type(struct_ens), intent(inout) :: ensPerts
    type(struct_vms), intent(inout) :: vModes

    ! Locals:
    type(struct_hco), pointer :: hco_ens
    type(struct_vco), pointer :: vco_ens
    real(4),          pointer :: ptr4d_1_r4(:,:,:,:)
    real(4),          pointer :: ptr4d_2_r4(:,:,:,:)
    real(8),      allocatable :: latWeight(:) ! Weight given to grid point in the statistic computation
    real(8),      allocatable :: sumWeight(:,:)
    integer :: varLevIndex1, varLevIndex2, levIndex1, levIndex2
    integer :: lonIndex, latIndex, nEns, memberIndex
    integer :: var3dIndex, myLonBeg, myLonEnd, myLatBeg, myLatEnd
    character(len=4), pointer :: varNamesList(:)
    
    !
    !- Structure initialization
    !
    nullify(varNamesList)
    call ens_varNamesList(varNamesList, ensPerts)

    vco_ens => ens_getVco(ensPerts)
    
    call vms_setup(varNamesList, vco_ens, & ! IN
                   vModes)                  ! OUT

    !
    !- Compute the vertical auto-covariance matrices
    !
    nEns = ens_getNumMembers(ensPerts)
    
    hco_ens => ens_getHco(ensPerts)
    allocate(latWeight(hco_ens%nj))
    do latIndex = 1, hco_ens%nj
      latWeight(latIndex) = cos(hco_ens%lat(latIndex))
      if (mmpi_myid == 0) then
        write(*,*) latIndex, hco_ens%lat(latIndex), cos(hco_ens%lat(latIndex))
      end if
    end do 
    
    call ens_getLatLonBounds(ensPerts, myLonBeg, myLonEnd, myLatBeg, myLatEnd)
    
    do var3dIndex = 1, vModes%nVar3d

      allocate(sumWeight(vModes%allVar3d(var3dIndex)%nLev,vModes%allVar3d(var3dIndex)%nLev))
      sumWeight(:,:)      = 0.d0
      
      vModes%allVar3d(var3dIndex)%autoCovariance(:,:) = 0.d0

      do memberIndex = 1, nEns

        !$OMP PARALLEL DO PRIVATE (varLevIndex1,varLevIndex2,levIndex2,levIndex1,ptr4d_1_r4,ptr4d_2_r4,latIndex,lonIndex)
        do levIndex2 = 1, vModes%allVar3d(var3dIndex)%nLev
          varLevIndex2 = levIndex2 + ens_getOffsetFromVarName(ensPerts,vModes%allVar3d(var3dIndex)%varName)
          do levIndex1 = 1, vModes%allVar3d(var3dIndex)%nLev
            varLevIndex1 = levIndex1 + ens_getOffsetFromVarName(ensPerts,vModes%allVar3d(var3dIndex)%varName)

            ptr4d_1_r4 => ens_getOneLev_r4(ensPerts,varLevIndex1)
            ptr4d_2_r4 => ens_getOneLev_r4(ensPerts,varLevIndex2)
          
            do latIndex = myLatBeg, myLatEnd
              do lonIndex = myLonBeg, myLonEnd

                vModes%allVar3d(var3dIndex)%autoCovariance(levIndex1,levIndex2) =                     &
                                     vModes%allVar3d(var3dIndex)%autoCovariance(levIndex1,levIndex2)  &
                                     + real(ptr4d_1_r4(memberIndex,1,lonIndex,latIndex),8)    &
                                     * real(ptr4d_2_r4(memberIndex,1,lonIndex,latIndex),8)    &
                                     * latWeight(latIndex)

                sumWeight(levIndex1,levIndex2) = sumWeight(levIndex1,levIndex2) + latWeight(latIndex)
                
              end do
            end do
            
          end do
        end do
        !$OMP END PARALLEL DO

      end do ! Loop on Ensemble

      !- Communication
      call mmpi_allreduce_sumR8_2d(vModes%allVar3d(var3dIndex)%autoCovariance, "GRID")
      call mmpi_allreduce_sumR8_2d(sumWeight, "GRID")

      !- Apply the appropriate scaling
      do levIndex2 = 1, vModes%allVar3d(var3dIndex)%nLev
        do levIndex1 = 1, vModes%allVar3d(var3dIndex)%nLev
          vModes%allVar3d(var3dIndex)%autoCovariance(levIndex1,levIndex2) =      &
               vModes%allVar3d(var3dIndex)%autoCovariance(levIndex1,levIndex2)   &
               / (sumWeight(levIndex1,levIndex2) / dble(nEns) * dble(nEns-1))
        end do
      end do

      deallocate(sumWeight)

    end do ! Loop on 3d variables

    !
    !- Compute the eigenvectors and eigenvalue
    !
    call vms_computeModes(vModes)

  end subroutine vms_computeModesFromEns

  !--------------------------------------------------------------------------
  ! vms_computeModesFromFunction
  !--------------------------------------------------------------------------
  subroutine vms_computeModesFromFunction(vco, lengthScaleTop, lengthScaleBot, &
                                          vModes)
    !
    ! :Purpose: To compute vertical modes from a prescribed correlation function
    !
    implicit none

    ! Arguments:
    type(struct_vco), pointer, intent(in)    :: vco
    real(8)         ,          intent(in)    :: lengthScaleTop
    real(8)         ,          intent(in)    :: lengthScaleBot
    type(struct_vms),          intent(inout) :: vModes

    ! Locals:
    real(8), pointer :: vertLocation_MM(:)
    real(8), pointer :: vertLocation_TH(:)
    real(8), pointer :: vertLocation(:)
    real(8) :: pSurfRef, distance, correl
    real(8) :: lengthScaleGradient
    real(8) :: lengthScaleLev1, lengthScaleLev2, lengthScaleAvg
    integer :: levIndex, levIndex1, levIndex2
    integer :: var3dIndex
    character(len=4) :: varNamesList(2)

    if (lengthScaleTop <= 0.d0 .or. lengthScaleBot <= 0.d0) then
      write(*,*) 'lengthScaleTop =', lengthScaleTop
      write(*,*) 'lengthScaleBot =', lengthScaleBot
      call utl_abort('vms_computeModes: problem with the provided lengthScales')
    end if

    !
    !- Structure initialization
    !
    varNamesList = (/'P_M','P_T'/)
    call vms_setup(varNamesList, vco, & ! IN
                   vModes)              ! OUT

    !
    !- Compute the vertical auto-covariance matrices
    !
    pSurfRef = 101000.D0
    call czp_fetch1DLevels(vco, pSurfRef,            & ! IN
                           profM_opt=vertLocation_MM)  ! OUT
    call czp_fetch1DLevels(vco, pSurfRef,            & ! IN
                           profT_opt=vertLocation_TH)  ! OUT
      
    do levIndex = 1, vco%nLev_M
      vertLocation_MM(levIndex) = log(vertLocation_MM(levIndex))
    end do
    do levIndex = 1, vco%nLev_T
      vertLocation_TH(levIndex) = log(vertLocation_TH(levIndex))
    end do

    call lfn_setup('FifthOrder') ! IN

    do var3dIndex = 1, vModes%nVar3d

      if (vnl_varLevelFromVarName(vModes%allVar3d(var3dIndex)%varName) == 'MM') then
        vertLocation => vertLocation_MM
      else
        vertLocation => vertLocation_TH
      end if

      lengthScaleGradient = (lengthScaleTop-lengthScaleBot) / &
                            (vertLocation(1)-vertLocation(vModes%allVar3d(var3dIndex)%nLev))
      
      do levIndex1 = 1, vModes%allVar3d(var3dIndex)%nLev
        do levIndex2 = 1, vModes%allVar3d(var3dIndex)%nLev
          distance = abs(vertLocation(levIndex2) - vertLocation(levIndex1))

          lengthScaleLev1 = lengthScaleGradient * &
               (vertLocation(levIndex1)-vertLocation(vModes%allVar3d(var3dIndex)%nLev)) + &
               lengthScaleBot
          lengthScaleLev2 = lengthScaleGradient * &
               (vertLocation(levIndex2)-vertLocation(vModes%allVar3d(var3dIndex)%nLev)) + &
               lengthScaleBot
          lengthScaleAvg = (lengthScaleLev1+lengthScaleLev2)/2.d0

          if (levIndex1 == 1 .and. mmpi_myid == 0) then
            write(*,*) 'vms_computeModesFromFunction: function length scale (', &
                        levIndex2,') = ', lengthScaleAvg
          end if
 
          correl = lfn_response(distance,lengthScaleAvg)
          vModes%allVar3d(var3dIndex)%autoCovariance(levIndex1,levIndex2) = correl
        end do
      end do

    end do

    !
    !- Compute the eigenvectors and eigenvalue
    !
    call vms_computeModes(vModes)

  end subroutine vms_computeModesFromFunction

  !--------------------------------------------------------------------------
  ! vms_computeModes
  !--------------------------------------------------------------------------
  subroutine vms_computeModes(vModes)
    !
    ! :Purpose: Compute vertical modes from the vertical covariances matrices
    !           contained in the vModes structure.
    !
    implicit none

    ! Arguments:
    type(struct_vms), intent(inout) :: vModes

    ! Locals:
    real(8) :: tolerance    
    integer :: levIndex1, levIndex2, var3dIndex
    integer :: matrixRank
    
    !
    !- Compute the eigenvectors and eigenvalue
    !
    tolerance = 1.0D-50
    do var3dIndex = 1, vModes%nVar3d

      !- Compute the eigenvectors and eigenvalue
      call utl_eigenDecomp(vModes%allVar3d(var3dIndex)%autoCovariance, & ! IN
                           vModes%allVar3d(var3dIndex)%eigenValues,    & ! OUT
                           vModes%allVar3d(var3dIndex)%eigenVectors,   & ! OUT
                           tolerance,                                  & ! IN
                           matrixRank)                                   ! OUT
      if ( matrixRank < vModes%allVar3d(var3dIndex)%nLev ) then
        write(*,*) 'varName    =', vModes%allVar3d(var3dIndex)%varName
        write(*,*) 'matrixRank =', matrixRank
        write(*,*) 'nLev       =', vModes%allVar3d(var3dIndex)%nLev
        call utl_abort('vms_computeModes: vertical matrix is rank deficient')
      end if

      !- Compute the inverse of the eigenVector matrix to go from grid point space
      !  to vertical mode space
      !    Note: The inverse of an orthogonal matrix is simply its transpose
      do levIndex2 = 1, vModes%allVar3d(var3dIndex)%nLev
        do levIndex1 = 1, vModes%allVar3d(var3dIndex)%nLev
          vModes%allVar3d(var3dIndex)%eigenVectorsInv(levIndex2,levIndex1) = &
               vModes%allVar3d(var3dIndex)%eigenVectors(levIndex1,levIndex2)
        end do
      end do

    end do
    
  end subroutine vms_computeModes

  !--------------------------------------------------------------------------
  ! vms_transform
  !--------------------------------------------------------------------------
  subroutine vms_transform(vModes, vertModesState, gridState,             &
                           TransformDirection, lonBeg, lonEnd, latBeg,    &
                           latEnd, nLev, varName, modeBeg_opt, modeEnd_opt)
    !
    ! :Purpose: To project back or forth ensemble pertubations onto the vertical
    !           modes contained in the vModes structure.
    !
    implicit none

    ! Arguments:
    type(struct_vms), intent(in)    :: vModes
    integer,          intent(in)    :: lonBeg
    integer,          intent(in)    :: lonEnd
    integer,          intent(in)    :: latBeg
    integer,          intent(in)    :: latEnd
    integer,          intent(in)    :: nLev
    real(8),          intent(inout) :: vertModesState(lonBeg:lonEnd,latBeg:latEnd,nLev) ! 3D vertical modes coefficients
    real(8),          intent(inout) :: gridState(lonBeg:lonEnd,latBeg:latEnd,nLev) ! 3D field in grid point space
    character(len=*), intent(in)    :: TransformDirection ! VertModesToGridPoint or GridPointToVertModes
    character(len=*), intent(in)    :: varName
    integer, optional,intent(in)    :: modeBeg_opt
    integer, optional,intent(in)    :: modeEnd_opt

    ! Locals:
    integer :: latIndex, lonIndex
    integer :: varIndexAssociated, varIndex
    integer :: nMode, modeBeg, modeEnd 
    
    !
    !- 1.  Find the appropriate modes for the provided varName
    !
    varIndexAssociated = -1

    ! Search for a match in varName
    do varIndex = 1, vModes%nVar3d
      if (vModes%allVar3d(varIndex)%varName == varName) then
        varIndexAssociated = varIndex
        exit
      end if
    end do

    ! Search for a match in level type (MM or TH)
    if (varIndexAssociated == -1) then
      if (vnl_varLevelFromVarName(varName) == 'MM') then
        do varIndex = 1, vModes%nVar3d
          if (vModes%allVar3d(varIndex)%varName == 'P_M') then
            varIndexAssociated = varIndex
            exit
          end if
        end do
      else if (vnl_varLevelFromVarName(varName) == 'TH') then
        do varIndex = 1, vModes%nVar3d
          if (vModes%allVar3d(varIndex)%varName == 'P_T') then
            varIndexAssociated = varIndex
            exit
          end if
        end do
      else
        write(*,*) trim(vModes%allVar3d(varIndex)%varName), vnl_varLevelFromVarName(varName)
        call utl_abort('vms_transform: varName is not on momentum or thermodynamic levels')
      end if
    end if

    if (varIndexAssociated == -1) then
      write(*,*)
      write(*,*) 'varName in input', trim(varName)
      write(*,*) 'varNames found in vModes'
      do varIndex = 1, vModes%nVar3d
        write(*,*) ' ... ', trim(vModes%allVar3d(varIndex)%varName)
      end do
      call utl_abort('vms_transform: could not match varName with any modes')
    end if

    if (nLev /= vModes%allVar3d(varIndexAssociated)%nLev) then
      call utl_abort('vms_transform: the number of levels are not consistent')
    end if
    nMode = nLev

    if (present(modeBeg_opt)) then
      modeBeg = modeBeg_opt
    else
      modeBeg = 1
    end if
    if (present(modeEnd_opt)) then
      if (modeEnd_opt /= -1) then
        modeEnd = modeEnd_opt
      else
        modeEnd = nMode
      end if
    else
      modeEnd = nMode
    end if
    
    !
    !- 2.  Do the transform
    !
    select case (trim(TransformDirection))
    case ('GridPointToVertModes')

      !$OMP PARALLEL DO PRIVATE (latIndex,lonIndex)
      do latIndex = latBeg, latEnd
        do lonIndex = lonBeg, lonEnd
          vertModesState(lonIndex,latIndex,modeBeg:modeEnd) = &
               matmul(vModes%allVar3d(varIndexAssociated)%eigenVectorsInv(modeBeg:modeEnd,1:nLev), &
                      gridState(lonIndex,latIndex,1:nLev))
        end do
      end do
      !$OMP END PARALLEL DO
      
    case ('VertModesToGridPoint')

      !$OMP PARALLEL DO PRIVATE (latIndex,lonIndex)
      do latIndex = latBeg, latEnd
        do lonIndex = lonBeg, lonEnd
          gridState(lonIndex,latIndex,1:nLev) = &
               matmul(vModes%allVar3d(varIndexAssociated)%eigenVectors(1:nLev,modeBeg:modeEnd), &
                      vertModesState(lonIndex,latIndex,modeBeg:modeEnd))
        end do
      end do
      !$OMP END PARALLEL DO

    case default
      write(*,*)
      write(*,*) 'Error: TranformDirection Unknown ', trim(TransformDirection)
      call utl_abort('vms_transform')
    end select
    
  end subroutine vms_transform
  
  !--------------------------------------------------------------------------
  ! vms_writeModes
  !--------------------------------------------------------------------------
  subroutine vms_writeModes(vModes)
    !
    ! :Purpose: To write the content of the provided vModes structure
    !
    implicit none

    ! Arguments:
    type(struct_vms), intent(in) :: vModes

    ! Locals:
    integer :: var3dIndex, levIndex1, levIndex2
    integer :: fstouv, fnom, fstfrm, fclos, iunstats, errorID
    character(len=4), allocatable :: varName3d(:)
    character(len=128) :: outfilename
    
    if (.not. vModes%initialized) then
      call utl_abort('vms_writeModes: The vModes structure is not initialized')
    end if

    if (mmpi_myid == 0) then
      !
      !- Write in FST file format
      !
      allocate(varName3d(vModes%nVar3d))
      iunstats = 0
      errorID = fnom(iunstats,'./vertAutoCov.fst','RND',0)
      errorID = fstouv(iunstats,'RND')      
      do var3dIndex = 1, vModes%nVar3d
        call writeMatrix2d_r8(vModes%allVar3d(var3dIndex)%autoCovariance, vModes%allVar3d(var3dIndex)%nlev, &
                              iunstats, vModes%allVar3d(var3dIndex)%varName,'AUTOCOVAR')
        call writeMatrix2d_r8(vModes%allVar3d(var3dIndex)%eigenVectors, vModes%allVar3d(var3dIndex)%nlev, &
                              iunstats, vModes%allVar3d(var3dIndex)%varName,'EIGENVECTORS')
        call writeArray1d_r8(vModes%allVar3d(var3dIndex)%eigenValues, vModes%allVar3d(var3dIndex)%nlev, &
                             iunstats, vModes%allVar3d(var3dIndex)%varName,'EIGENVALUES')
        varName3d(var3dIndex) = vModes%allVar3d(var3dIndex)%varName
      end do
      call writeArray1d_c4(varName3d, vModes%nVar3d, iunstats, 'VAR', 'VARNAME')
      errorID = fstfrm(iunstats)
      errorID = fclos (iunstats)

      !
      !- Write in txt format (for plotting purposes)
      !
      do var3dIndex = 1, vModes%nVar3d
        outfilename = "./vModes_autoCovar_"//trim(vModes%allVar3d(var3dIndex)%varName)//".txt"
        open (unit=99,file=outfilename,action="write",status="new")
        do levIndex2 = 1, vModes%allVar3d(var3dIndex)%nLev
          do levIndex1 = 1, vModes%allVar3d(var3dIndex)%nLev
            if (levIndex1 == vModes%allVar3d(var3dIndex)%nLev) then
              write(99,'(2X,E10.4)') vModes%allVar3d(var3dIndex)%autoCovariance(levIndex1,levIndex2)
            else
              write(99,'(2X,E10.4,$)') vModes%allVar3d(var3dIndex)%autoCovariance(levIndex1,levIndex2)
            end if
          end do
        end do
        close(unit=99)
        outfilename = "./vModes_eigenVectors_"//trim(vModes%allVar3d(var3dIndex)%varName)//".txt"
        open (unit=99,file=outfilename,action="write",status="new")
        do levIndex2 = 1, vModes%allVar3d(var3dIndex)%nLev
          do levIndex1 = 1, vModes%allVar3d(var3dIndex)%nLev
            if (levIndex1 == vModes%allVar3d(var3dIndex)%nLev) then
              write(99,'(2X,E10.4)') vModes%allVar3d(var3dIndex)%eigenVectors(levIndex1,levIndex2)
            else
              write(99,'(2X,E10.4,$)') vModes%allVar3d(var3dIndex)%eigenVectors(levIndex1,levIndex2)
            end if
          end do
        end do
        close(unit=99)
        outfilename = "./vModes_eigenValues_"//trim(vModes%allVar3d(var3dIndex)%varName)//".txt"
        open (unit=99,file=outfilename,action="write",status="new")
        do levIndex1 = 1, vModes%allVar3d(var3dIndex)%nLev
          write(99,'(2X,E10.4)') vModes%allVar3d(var3dIndex)%eigenValues(levIndex1)
        end do
        close(unit=99)
      end do
    end if

  end subroutine vms_writeModes

  !--------------------------------------------------------------------------
  ! writeMatrix2d_r8
  !--------------------------------------------------------------------------
  subroutine writeMatrix2d_r8(matrix2d,rank,iun,nomvar_in,etiket_in)
    !
    ! :Purpose: To write a 2D matrix rank x rank
    !
    implicit none

    ! Arguments:
    integer,          intent(in) :: rank
    real(8),          intent(in) :: matrix2d(rank,rank)
    integer,          intent(in) :: iun
    character(len=*), intent(in) :: nomvar_in
    character(len=*), intent(in) :: etiket_in

    ! Locals:
    real(4), allocatable :: workecr(:,:)
    real(4) :: work
    integer :: errorID, fstecr
    integer :: dateo, npak, ni, nj, nk
    integer :: ip1, ip2, ip3, deet, npas, datyp
    integer :: ig1 ,ig2 ,ig3 ,ig4
    character(len=1)  :: grtyp
    character(len=4)  :: nomvar
    character(len=2)  :: typvar
    character(len=12) :: etiket

    allocate(workecr(rank, rank))

    npak   = -32
    dateo  = 0
    deet   = 0
    npas   = 0
    ni     = rank
    nj     = rank
    nk     = 1
    ip1    = 0
    ip2    = 0
    ip3    = 0
    typvar = 'XX'
    nomvar = nomvar_in
    etiket = etiket_in
    grtyp  = 'X'
    ig1    = 0
    ig2    = 0
    ig3    = 0
    ig4    = 0
    datyp  = 5

    !- Convert to real 4
    workecr(:,:) = real(matrix2d(:,:),4)

    !- Writing 
    errorID = fstecr(workecr, work, npak, iun, dateo, deet, npas, ni, nj, &
                     nk, ip1, ip2, ip3, typvar, nomvar, etiket, grtyp,    &
                     ig1, ig2, ig3, ig4, datyp, .true.)

    deallocate(workecr)

  end subroutine writeMatrix2d_r8

  !--------------------------------------------------------------------------
  ! writeArray1d_r8
  !--------------------------------------------------------------------------
  subroutine writeArray1d_r8(array,size,iun,nomvar_in,etiket_in)
    !
    ! :Purpose: To write a 1D array in real 8
    !
    implicit none

    ! Arguments:
    integer,          intent(in) :: size
    real(8),          intent(in) :: array(size)
    integer,          intent(in) :: iun
    character(len=*), intent(in) :: nomvar_in
    character(len=*), intent(in) :: etiket_in

    ! Locals:
    real(4), allocatable :: workecr(:)
    real(4) :: work
    integer :: errorID, fstecr
    integer :: dateo, npak, ni, nj, nk
    integer :: ip1, ip2, ip3, deet, npas, datyp
    integer :: ig1 ,ig2 ,ig3 ,ig4
    character(len=1)  :: grtyp
    character(len=4)  :: nomvar
    character(len=2)  :: typvar
    character(len=12) :: etiket

    allocate(workecr(size))

    npak   = -32
    dateo  = 0
    deet   = 0
    npas   = 0
    ni     = size
    nj     = 1
    nk     = 1
    ip1    = 0
    ip2    = 0
    ip3    = 0
    typvar = 'XX'
    nomvar = nomvar_in
    etiket = etiket_in
    grtyp  = 'X'
    ig1    = 0
    ig2    = 0
    ig3    = 0
    ig4    = 0
    datyp  = 5

    !- Convert to real 4
    workecr(:) = real(array(:),4)

    !- Writing 
    errorID = fstecr(workecr, work, npak, iun, dateo, deet, npas, ni, nj, &
                     nk, ip1, ip2, ip3, typvar, nomvar, etiket, grtyp,    &
                     ig1, ig2, ig3, ig4, datyp, .true.)

    deallocate(workecr)
    
  end subroutine writeArray1d_r8

  !--------------------------------------------------------------------------
  ! writeArray1d_c4
  !--------------------------------------------------------------------------
  subroutine writeArray1d_c4(array,size,iun,nomvar_in,etiket_in)
    !
    ! :Purpose: To write a 1D array in character with len = 4
    !
    implicit none

    ! Arguments:
    integer,          intent(in) :: size
    character(len=4), intent(in) :: array(size)
    integer,          intent(in) :: iun
    character(len=*), intent(in) :: nomvar_in
    character(len=*), intent(in) :: etiket_in

    ! Locals:
    real(4) :: work
    integer :: errorID, fstecr_s
    integer :: dateo, npak, ni, nj, nk
    integer :: ip1, ip2, ip3, deet, npas, datyp
    integer :: ig1 ,ig2 ,ig3 ,ig4
    character(len=1)  :: grtyp
    character(len=4)  :: nomvar
    character(len=2)  :: typvar
    character(len=12) :: etiket

    npak   = -32
    dateo  = 0
    deet   = 0
    npas   = 0
    ni     = 4 ! 4 Characters
    nj     = size
    nk     = 1
    ip1    = 0
    ip2    = 0
    ip3    = 0
    typvar = 'XX'
    nomvar = nomvar_in
    etiket = etiket_in
    grtyp  = 'X'
    ig1    = 0
    ig2    = 0
    ig3    = 0
    ig4    = 0
    datyp  = 7 ! Character

    !- Writing
    errorID = fstecr_s(array, work, npak, iun, dateo, deet, npas, ni, nj, &
                      nk, ip1, ip2, ip3, typvar, nomvar, etiket, grtyp,    &
                      ig1, ig2, ig3, ig4, datyp, .true.)

  end subroutine writeArray1d_c4
  
end module verticalModes_mod
