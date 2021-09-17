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

module menetrierDiag_mod
  ! MODULE menetrierDiag_mod (prefix='bmd' category='1. High-level functionality')
  !
  ! :Purpose: To compute optimal localization radii according to the theory 
  !           developed by Benjamin Menetrier (Meteo-France) and reported
  !           in Menetrier, Michel, Montmerle and Berre, 2015, Parts 1 and 2.
  !
  use earthConstants_mod
  use utilities_mod
  use localizationFunction_mod
  use varNameList_mod
  use horizontalCoord_mod
  use gridStatevector_mod
  use ensembleStatevector_mod
  use timeCoord_mod
  use mpi_mod
  implicit none
  save
  private

  real(8), pointer :: pressureProfile_M(:), pressureProfile_T(:)
  
  type(struct_hco), pointer :: hco_ens ! Ensemble horizontal grid parameters

  integer :: nens, ni, nj, nLevEns_M, nLevEns_T, nkgdimEns
  integer :: nvar3d, nvar2d, nWaveBand

  integer, allocatable :: varLevOffset(:)

  character(len=4), allocatable :: nomvar3d(:), nomvar2d(:)

  integer :: strideForHLoc, strideForVloc, horizPadding
  logical :: hLoc, vLoc, global

  logical :: initialized = .false.

  ! Public Subroutines
  public :: bmd_setup, bmd_localizationRadii

contains

  !--------------------------------------------------------------------------
  ! bmd_setup
  !--------------------------------------------------------------------------
  subroutine bmd_setup(statevector_template, hco_core_in, nens_in, pressureProfile_M_in, &
                       pressureProfile_T_in, nWaveBand_in)
    implicit none

    type(struct_gsv) :: statevector_template
    type(struct_hco), pointer :: hco_core_in

    integer, intent(in) :: nens_in
    integer, intent(in) :: nWaveBand_in

    real(8), intent(in), pointer :: pressureProfile_M_in(:), pressureProfile_T_in(:)

    integer :: nulnam, ierr, fclos, fnom
    integer :: nVar, varNameIndex, var2dIndex, var3dIndex

    character(len=4), pointer :: varNamesList(:)

    NAMELIST /NAMLOCALIZATIONRADII/strideForHLoc,strideForVloc,hLoc,vLoc,horizPadding

    !
    !- 1.  Input parameters 
    !
    hco_ens   => hco_core_in
    nens      = nens_in
    ni        = hco_ens%ni
    nj        = hco_ens%nj
    nLevEns_M = gsv_getNumLev(statevector_template,'MM') !nLevEns_M_in
    nLevEns_T = gsv_getNumLev(statevector_template,'TH') !nLevEns_T_in
    nkgdimEns = statevector_template%nk
    pressureProfile_M => pressureProfile_M_in
    pressureProfile_T => pressureProfile_T_in
    nWaveBand = nWaveBand_in
    global    = hco_ens%global

    nullify(varNamesList)
    call gsv_varNamesList(varNamesList, statevector_template)
    
    nVar = size(varNamesList)
    allocate(varLevOffset(nVar))
    nVar3d = 0
    nVar2d = 0
    do varNameIndex = 1, size(varNamesList)
      varLevOffset(varNameIndex) = gsv_getOffsetFromVarName(statevector_template,varNamesList(varNameIndex))
      if (vnl_varLevelFromVarname(varNamesList(varNameIndex)) == 'SF' ) then
        nVar2d = nVar2d + 1
      else
        nVar3d = nVar3d + 1
      end if
    end do

    allocate(nomvar3d(nvar3d))
    allocate(nomvar2d(nvar2d))

    var2dIndex = 0
    var3dIndex = 0
    do varNameIndex = 1, size(varNamesList)
      if (vnl_varLevelFromVarname(varNamesList(varNameIndex)) == 'SF' ) then
        var2dIndex = var2dIndex + 1
        nomvar2d(var2dIndex) =  varNamesList(varNameIndex)
      else
        var3dIndex = var3dIndex + 1
        nomvar3d(var3dIndex) =  varNamesList(varNameIndex)
      end if
    end do

    !
    !- 2.  Read Namelist
    !
    hLoc          = .true.
    vLoc          = .true.
    strideForHLoc = 100 ! Horizontal correlations will be computed every "stride" gridpoint in x and y
    strideForVLoc = 50  ! Vertical   correlations will be computed every "stride" gridpoint in x and y
    horizPadding  = 0   ! Number of grid point to discard along the horizontal edges (only for LAM)

    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=NAMLOCALIZATIONRADII)
    write(*,nml=NAMLOCALIZATIONRADII)
    ierr = fclos(nulnam)

    !
    !- 3.  Ending
    !
    initialized = .true.

  end subroutine bmd_setup

  !--------------------------------------------------------------------------
  ! bmd_localizationRadii
  !--------------------------------------------------------------------------
  subroutine bmd_localizationRadii(ensPerts,waveBandIndex_opt)
    implicit none

    type(struct_ens) :: ensPerts
    integer,optional, intent(in) :: waveBandIndex_opt

    !
    !- Estimate the horionzontal and vertical localization radii
    !
    if ( hLoc ) then
      call calcHorizLocalizationRadii(ensPerts, strideForHLoc,           & ! IN
                                      waveBandIndex_opt=waveBandIndex_opt) ! IN
    end if
    if ( vLoc ) then
      call calcVertLocalizationRadii(ensPerts, strideForVloc,           &  ! IN
                                     waveBandIndex_opt=waveBandIndex_opt)  ! IN
    end if

  end subroutine bmd_localizationRadii

  !--------------------------------------------------------------------------
  ! CALCHORIZLOCALIZATIONRADII
  !--------------------------------------------------------------------------
  subroutine calcHorizLocalizationRadii(ensPerts,stride,waveBandIndex_opt)
    implicit none

    type(struct_ens)    :: ensPerts
    integer, intent(in) :: stride
    integer,optional, intent(in) :: waveBandIndex_opt

    type(struct_gsv) :: statevector_ensStdDev
    type(struct_gsv) :: statevector_ensStdDev_tiles
    type(struct_gsv) :: statevector_oneMemberTiles
    type(struct_gsv) :: statevector_oneMember(nens)

    real(8), pointer :: ensStdDev(:,:,:)
    real(4), pointer :: ptr3d_r4(:,:,:)

    real(4), allocatable :: ensPert_local(:,:,:)

    real(8), allocatable :: meanCorrel(:,:)
    real(8), allocatable :: meanCorrelSquare(:,:)
    real(8), allocatable :: meanVarianceProduct(:,:)
    real(8), allocatable :: meanCovarianceSquare(:,:)
    real(8), allocatable :: meanFourthMoment(:,:)

    real(8), allocatable :: meanCorrel_local(:,:)
    real(8), allocatable :: meanCorrelSquare_local(:,:)
    real(8), allocatable :: meanVarianceProduct_local(:,:)
    real(8), allocatable :: meanCovarianceSquare_local(:,:)
    real(8), allocatable :: meanFourthMoment_local(:,:)

    real(8), allocatable :: localizationFunctions(:,:,:) ! Eq. 19-21 in MMMB 2015 Part 2

    real(8), allocatable :: localizationRadii(:,:)

    real(8), allocatable :: distanceBinThresholds(:) ! Maximum distance for each distance-bin
    real(8), allocatable :: distanceBinMean(:)       ! Mean distance for each distance-bin
    real(8), allocatable :: distanceBinWeight(:)     ! Weight given to each bin in the curve fitting step

    real(8), allocatable :: gridPointWeight(:,:,:)   ! Weight given to grid point in the statistic computation

    real(8), allocatable :: sumWeight(:,:)    ! Sample size for each distance-bin
    real(8), allocatable :: sumWeight_local(:,:)

    real(8), pointer :: PressureProfile(:)

    logical, allocatable :: gridPointAlreadyUsed(:,:)

    real(8) :: dnens, correlation, covariance, fourthMoment, distance, maxDistance, weight
    real(8) :: t1, t2, t3, rmse

    integer :: i, j, k, f, ens, bin, numbins, numFunctions, nSize
    integer :: iref, jref, ier
    integer :: nLevEns, jvar, mykBeg, mykEnd

    character(len=128) :: outfilename
    character(len=2)   :: wbnum
    character(len=4), pointer :: varNamesList(:)

    !
    !- 1.  Setup
    !
    call ens_copyEnsStdDev(ensPerts, statevector_ensStdDev_tiles)

    nullify(varNamesList)
    call ens_varNamesList(varNamesList,ensPerts)
    call gsv_allocate(statevector_ensStdDev, ens_getNumStep(ensPerts),                       &
                      ens_getHco(ensPerts), ens_getVco(ensPerts), varNames_opt=varNamesList, &
                      datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true.,                &
                      mpi_distribution_opt='VarsLevs', dataKind_opt=8 )
    call gsv_transposeTilesToVarsLevs(statevector_ensStdDev_tiles, statevector_ensStdDev)
    call gsv_getField(statevector_ensStdDev,ensStdDev)

    call gsv_deallocate(statevector_ensStdDev_tiles)

    mykBeg = statevector_ensStdDev%mykBeg
    mykEnd = statevector_ensStdDev%mykEnd

    numFunctions = 3

    if ( global ) then
      numBins = ni/4 ! 1/4 of the Earth
      if ( horizPadding /= 0 ) then
        write(*,*)
        write(*,*) 'WARNING: horizPadding /= 0. Force it to zero since we are in global mode.'
        horizPadding = 0
      end if
    else
      numBins = ni-(2*horizPadding)
      if ( horizPadding == 0 ) then
        write(*,*)
        write(*,*) 'WARNING: horizPadding == 0. The rim and blending zone WILL HAVE an impact on the results.'
      end if
    end if

    allocate(meanCorrelSquare(numBins,nkgdimEns))
    allocate(meanCorrel(numBins,nkgdimEns))
    allocate(meanVarianceProduct(numBins,nkgdimEns))
    allocate(meanCovarianceSquare(numBins,nkgdimEns))
    allocate(meanFourthMoment(numBins,nkgdimEns))
    allocate(sumWeight(numBins,nkgdimEns))

    allocate(meanCorrelSquare_local(numBins,nkgdimEns))
    allocate(meanCorrel_local(numBins,nkgdimEns))
    allocate(meanVarianceProduct_local(numBins,nkgdimEns))
    allocate(meanCovarianceSquare_local(numBins,nkgdimEns))
    allocate(meanFourthMoment_local(numBins,nkgdimEns))
    allocate(sumWeight_local(numBins,nkgdimEns))

    allocate(gridPointWeight(ni,nj,nkgdimEns))

    allocate(distanceBinThresholds(numBins))

    ! Assign the (upper) threshold to each separation-distance-bin
    write(*,*)
    write(*,*) 'Separation distance bin'
    do bin = 1, numbins
      distanceBinThresholds(bin) = calcDistance(hco_ens%lat(nj/2),hco_ens%lon(1+horizPadding),hco_ens%lat(nj/2),hco_ens%lon(bin+horizPadding))
      write(*,*) bin, hco_ens%lat(nj/2),hco_ens%lon(1+horizPadding),hco_ens%lon(bin+horizPadding), distanceBinThresholds(bin)/1000.d0
    end do

    maxDistance=distanceBinThresholds(numBins)

    dnens = 1.0d0/dble(nens-1)

    ! Grid point Weight
    write(*,*)
    write(*,*) 'Grid point Weight'
    do j = 1, nj
      gridPointWeight(:,j,:) = cos(hco_ens%lat(j))
      write(*,*) j, hco_ens%lat(j), cos(hco_ens%lat(j))
    end do

    !
    !- 2.  Estimation of localization functions
    !

    !- 2.1  Computation of various statistics for different separation distances
    meanCorrelSquare_local(:,:)     = 0.d0
    meanCorrel_local(:,:)           = 0.d0
    meanCovarianceSquare_local(:,:) = 0.d0
    meanVarianceProduct_local(:,:)  = 0.d0
    meanFourthMoment_local(:,:)     = 0.d0
    sumWeight_local(:,:) = 0.d0

    allocate(ensPert_local(nens,ni,nj))
    allocate(gridPointAlreadyUsed(ni,nj))

    call gsv_allocate(statevector_oneMemberTiles, ens_getNumStep(ensPerts),                  &
                      ens_getHco(ensPerts), ens_getVco(ensPerts), varNames_opt=varNamesList, &
                      datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true.,                &
                      mpi_distribution_opt='Tiles', dataKind_opt=4 )

    do ens = 1, nens
      call gsv_allocate(statevector_oneMember(ens), ens_getNumStep(ensPerts),      &
                        ens_getHco(ensPerts), ens_getVco(ensPerts), varNames_opt=varNamesList, &
                        datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true.,                &
                        mpi_distribution_opt='VarsLevs', dataKind_opt=4 )
      call ens_copyMember(ensPerts, statevector_oneMemberTiles, ens)
      call gsv_transposeTilesToVarsLevs(statevector_oneMemberTiles, statevector_oneMember(ens))
    end do

    call gsv_deallocate(statevector_oneMemberTiles)

    !$OMP PARALLEL DO PRIVATE (ptr3d_r4,k,iref,jref,j,i,ens,correlation,covariance,fourthMoment,distance,bin,weight,ensPert_local,gridPointAlreadyUsed)
    do k = mykBeg, mykEnd
      write(*,*) 'Computing distance-bin statistics for ensemble level: ', k

      !- Select data needed to speed up the process (ensemble member index must come first in ensPert_Local 
      !  because ensemble member is the last loop index below)
      do ens = 1, nens
        call gsv_getField(statevector_oneMember(ens),ptr3d_r4)
        do j = 1, nj
          do i = 1, ni
            ensPert_local(ens,i,j) = ptr3d_r4(i,j,k)
          end do
        end do
      end do
      
      gridPointAlreadyUsed(:,:) = .false.

      do jref = nint(stride/2.0)+horizPadding, nj-horizPadding, stride    ! Pick every stride point to save cost.
        do iref = nint(stride/2.0)+horizPadding, ni-horizPadding, stride  ! Pick every stride point to save cost.

          if (k == 1) then
            write(*,*) 'grid point info', iref, jref
          end if

          do j = 1+horizPadding, nj-horizPadding
            do i = 1+horizPadding, ni-horizPadding
              
              if ( gridPointAlreadyUsed(i,j) ) cycle ! prevent using the same pair of points more than once

              distance=calcDistance(hco_ens%lat(jref),hco_ens%lon(iref),hco_ens%lat(j),hco_ens%lon(i))
              if (distance <= maxDistance .and. gridPointWeight(i,j,k) > 0.d0) then
                covariance = 0.d0
                fourthMoment = 0.d0
                do ens = 1, nens
                  covariance = covariance + &
                       ensPert_local(ens,i,j) * ensPert_local(ens,iref,jref)
                  fourthMoment = fourthMoment + &
                       (ensPert_local(ens,i,j) * ensPert_local(ens,iref,jref))**2
                end do
                covariance   = covariance * dnens
                fourthMoment = fourthMoment / real(nens,8)
                if ( ensStdDev(i,j,k) > 0.d0 .and. ensStdDev(iref,jref,k) > 0.d0 ) then
                  correlation  = covariance / (ensStdDev(i,j,k)*ensStdDev(iref,jref,k))
                else
                  correlation  = 0.d0
                end if

                bin=findBinIndex(distance,distanceBinThresholds,numBins)

                weight = sqrt(gridPointWeight(iref,jref,k)) * sqrt(gridPointWeight(i,j,k))

                sumWeight_local(bin,k) = sumWeight_local(bin,k) + weight

                meanCorrel_local(bin,k)           = meanCorrel_local(bin,k)           + &
                                                    correlation    * weight
                meanCorrelSquare_local(bin,k)     = meanCorrelSquare_local(bin,k)     + &
                                                    correlation**2 * weight
                meanCovarianceSquare_local(bin,k) = meanCovarianceSquare_local(bin,k) + &
                                                    covariance**2  * weight
                meanVarianceProduct_local(bin,k)  = meanVarianceProduct_local(bin,k)  + &
                                                    ensStdDev(i,j,k)**2 * ensStdDev(iref,jref,k)**2 * weight
                meanFourthMoment_local(bin,k)     = meanFourthMoment_local(bin,k)     + &
                                                    fourthMoment   * weight
              end if
            end do
          end do

          ! From now on, omit the current reference point
          gridPointAlreadyUsed(iref,jref) = .true.

        end do ! iref
      end do   ! jref

    end do ! nkgdimEns
    !$OMP END PARALLEL DO

    deallocate(ensPert_local)
    deallocate(gridPointAlreadyUsed)
    do ens = 1, nens
      call gsv_deallocate(statevector_oneMember(ens))
    end do

    !- 2.2 Gather the all the info in processor 0
    nSize = nkgdimEns * numbins
    call rpn_comm_reduce(sumWeight_local           ,sumWeight           ,nSize, &
         "mpi_double_precision","mpi_sum",0,"GRID",ier)
    call rpn_comm_reduce(meanCorrel_local          ,meanCorrel          ,nSize, &
         "mpi_double_precision","mpi_sum",0,"GRID",ier)
    call rpn_comm_reduce(meanCorrelSquare_local    ,meanCorrelSquare    ,nSize, &
         "mpi_double_precision","mpi_sum",0,"GRID",ier)
    call rpn_comm_reduce(meanVarianceProduct_local ,meanVarianceProduct ,nSize, &
         "mpi_double_precision","mpi_sum",0,"GRID",ier)
    call rpn_comm_reduce(meanFourthMoment_local    ,meanFourthMoment    ,nSize, &
         "mpi_double_precision","mpi_sum",0,"GRID",ier)
    call rpn_comm_reduce(meanCovarianceSquare_local,meanCovarianceSquare,nSize, &
         "mpi_double_precision","mpi_sum",0,"GRID",ier)

    deallocate(sumWeight_local)
    deallocate(meanCorrel_local)
    deallocate(meanCorrelSquare_local)
    deallocate(meanVarianceProduct_local)
    deallocate(meanFourthMoment_local)
    deallocate(meanCovarianceSquare_local)

    if (mpi_myid == 0) then
      !- 2.3  Computation of the localization functions
      allocate(localizationFunctions(numFunctions,numBins,nkgdimEns))

      t1=dble((nens-1)**2)/dble(nens*(nens-3))
      t2=dble(nens)/dble((nens-2)*(nens-3))
      t3=dble(nens-1)/dble(nens*(nens-2)*(nens-3))
      !$OMP PARALLEL DO PRIVATE (k,bin)
      do k = 1, nkgdimEns
        do bin = 1, numbins
          
          meanCorrel(bin,k)           = meanCorrel(bin,k)           / sumWeight(bin,k)
          meanCorrelSquare(bin,k)     = meanCorrelSquare(bin,k)     / sumWeight(bin,k)
          meanCovarianceSquare(bin,k) = meanCovarianceSquare(bin,k) / sumWeight(bin,k)
          meanVarianceProduct(bin,k)  = meanVarianceProduct(bin,k)  / sumWeight(bin,k)
          meanFourthMoment(bin,k)     = meanFourthMoment(bin,k)     / sumWeight(bin,k)
          
          if ( meanCovarianceSquare(bin,k) /= 0.d0 ) then
            ! Form 1: General formulation (Eq. 19 in MMMB 2015 Part 2)
            localizationFunctions(1,bin,k) = t1 - t2*meanFourthMoment(bin,k)/meanCovarianceSquare(bin,k) + &
                 t3*meanVarianceProduct(bin,k)/meanCovarianceSquare(bin,k)
            ! Form 2: Gaussian sample distribution (Eq. 20 in MMMB 2015 Part 2)
            localizationFunctions(2,bin,k) = dble(nens-1)/dble((nens+1)*(nens-2)) * &
                 (dble(nens-1)-meanVarianceProduct(bin,k)/meanCovarianceSquare(bin,k))
          else
            write(*,*) ' !!! Warning !!! meanCovarianceSquare = 0 in bin, level = ', bin, k
            localizationFunctions(1,bin,k) = 0.d0
            localizationFunctions(2,bin,k) = 0.d0
          end if
          ! Form 3: Gaussian sample distribution and correlation-based formulation (Eq. 21 in MMMB 2015 Part 2)
          if ( meanCorrelSquare(bin,k) /= 0.d0 ) then
            localizationFunctions(3,bin,k) = dble(nens-1)/dble((nens+1)*(nens-2)) * &
                 (dble(nens-1)-1.d0/meanCorrelSquare(bin,k))
          else
            write(*,*) ' !!! Warning !!! meanCorrelSquare = 0 in bin, level = ', bin, k
            localizationFunctions(3,bin,k) = 0.d0
          end if
        end do
      end do
      !$OMP END PARALLEL DO

      !
      !- 3.  Estimation of localization radii (curve fitting)
      !
      allocate(localizationRadii(numFunctions,nkgdimEns))
      allocate(distanceBinMean(numBins))
      allocate(distanceBinWeight(numBins))
      
      call lfn_setup('FifthOrder') ! IN
      
      localizationRadii(:,:) = 2000.d0*1000.d0 ! First Guess (meter)
      distanceBinWeight(:)   = 1.d0            ! Even weight
      do bin = 1, numBins
        if (bin == 1) then
          distanceBinMean(bin) = distanceBinThresholds(bin) ! = 0 m
        else
          distanceBinMean(bin) = 0.5d0*(distanceBinThresholds(bin)+distanceBinThresholds(bin-1))
        end if
      end do
      
      do f = 1, numFunctions
        write(*,*)
        write(*,*) 'Localization Function : ', f
        do k =  1, nkgdimEns
          write(*,*) '         ----- EnsLev : ', k
          call lfn_lengthscale( localizationRadii(f,k),       & ! INOUT
               rmse,                         & ! OUT
               localizationFunctions(f,:,k), & ! IN
               distanceBinMean,              & ! IN
               distanceBinWeight,            & ! IN
               numbins )                       ! IN
        end do
      end do
      
      !
      !- 4.  Write to file
      !
      if ( nWaveBand /= 1 ) then
        if (.not. present(waveBandIndex_opt)) then
          write(*,*) 'calcLocalizationRadii: No waveBandIndex was supplied!!!'
          call utl_abort('calbmatrix_glb')
        end if
        write(wbnum,'(I2.2)') waveBandIndex_opt
      end if
      
      !- 4.1 Localization functions in txt format (for plotting purposes)
      do jvar = 1, nvar3d
        if ( nWaveBand == 1 ) then
          outfilename = "./horizLocalizationFunctions_"//trim(nomvar3d(jvar))//".txt"
        else
          outfilename = "./horizLocalizationFunctions_"//trim(nomvar3d(jvar))//"_"//wbnum//".txt"
        end if
        open (unit=99,file=outfilename,action="write",status="new")
        
        if(vnl_varLevelFromVarName(nomvar3d(jvar)).eq.'MM') then
          nLevEns = nLevEns_M
        else
          nLevEns = nLevEns_T
        endif
        do k=1,nlevEns
          do bin = 1, numbins
            write(99,'(I3,2X,I3,2X,F7.1,2X,I7,2X,F6.4,2X,F6.4,2X,F6.4,2X,F6.4,2X,F6.4,2X,E10.3,2X,E10.3,2X,E10.3)') k, bin, &
                 distanceBinThresholds(bin)/1000.d0, nint(sumWeight(bin,varLevOffset(jvar)+k)), &
                 meanCorrel(bin,varLevOffset(jvar)+k), &
                 localizationFunctions(3,bin,varLevOffset(jvar)+k), &
                 localizationFunctions(2,bin,varLevOffset(jvar)+k), &
                 localizationFunctions(1,bin,varLevOffset(jvar)+k), &
                 meanCorrelSquare(bin,varLevOffset(jvar)+k), &
                 meanCovarianceSquare(bin,varLevOffset(jvar)+k), &
                 meanVarianceProduct(bin,varLevOffset(jvar)+k), &
                 meanFourthMoment(bin,varLevOffset(jvar)+k)
          end do
        end do
        close(unit=99)
      end do
      
      do jvar = 1, nvar2d
        k = varLevOffset(nvar3d+1)+jvar
        if ( nWaveBand == 1 ) then
          outfilename = "./horizLocalizationFunctions_"//trim(nomvar2d(jvar))//".txt"
        else
          outfilename = "./horizLocalizationFunctions_"//trim(nomvar2d(jvar))//"_"//wbnum//".txt"
        end if
        open (unit=99,file=outfilename,action="write",status="new")
        do bin = 1, numbins
          write(99,'(I3,2X,I3,2X,F7.1,2X,I7,2X,F6.4,2X,F6.4,2X,F6.4,2X,F6.4,2X,F6.4,2X,E10.3,2X,E10.3,2X,E10.3)') 1, bin, &
               distanceBinThresholds(bin)/1000.d0, nint(sumWeight(bin,k)), &
               meanCorrel(bin,k), &
               localizationFunctions(3, bin,k), &
               localizationFunctions(2, bin,k), &
               localizationFunctions(1, bin,k), &
               meanCorrelSquare(bin,k), &
               meanCovarianceSquare(bin,k), &
               meanVarianceProduct(bin,k), &
               meanFourthMoment(bin,k)
        end do
        close(unit=99)
      end do
      
      !- 4.2 Localization radii in txt format (for plotting purposes)
      do jvar = 1, nvar3d
        if ( nWaveBand == 1 ) then
          outfilename = "./horizLocalizationRadii_"//trim(nomvar3d(jvar))//".txt"
        else
          outfilename = "./horizLocalizationRadii_"//trim(nomvar3d(jvar))//"_"//wbnum//".txt"
        end if
        open (unit=99,file=outfilename,action="write",status="new")
        
        if(vnl_varLevelFromVarName(nomvar3d(jvar)).eq.'MM') then
          nLevEns = nLevEns_M
          PressureProfile => pressureProfile_M
        else
          nLevEns = nLevEns_T
          PressureProfile => pressureProfile_T
        endif
        do k=1,nlevEns
          write(99,'(I3,2X,F7.2,2X,F7.1,2X,F7.1,2X,F7.1)') k, PressureProfile(k)/100.d0, &
               min(localizationRadii(3,varLevOffset(jvar)+k)/1000.d0,99999.9d0), &
               min(localizationRadii(2,varLevOffset(jvar)+k)/1000.d0,99999.9d0), &
               min(localizationRadii(1,varLevOffset(jvar)+k)/1000.d0,99999.9d0)
        end do
        close(unit=99)
      end do
      
      do jvar = 1, nvar2d
        k = varLevOffset(nvar3d+1)+jvar
        if ( nWaveBand == 1 ) then
          outfilename = "./horizLocalizationRadii_"//trim(nomvar2d(jvar))//".txt"
        else
          outfilename = "./horizLocalizationRadii_"//trim(nomvar2d(jvar))//"_"//wbnum//".txt"
        end if
        open (unit=99,file=outfilename,action="write",status="new")
        write(99,'(I3,2X,F7.2,2X,F7.1,2X,F7.1,2X,F7.1)') 1, 1010.0, &
             min(localizationRadii(3,k)/1000.d0,99999.9d0), &
             min(localizationRadii(2,k)/1000.d0,99999.9d0), &
             min(localizationRadii(1,k)/1000.d0,99999.9d0)
        close(unit=99)
      end do

      deallocate(localizationFunctions)
      deallocate(localizationRadii)
      deallocate(distanceBinMean)
      deallocate(distanceBinWeight)

    end if ! mpi_myid == 0

    deallocate(meanCorrelSquare)
    deallocate(meanCorrel)
    deallocate(meanCovarianceSquare)
    deallocate(meanVarianceProduct)
    deallocate(meanFourthMoment)
    deallocate(distanceBinThresholds)
    deallocate(sumWeight)
    deallocate(gridPointWeight)

    write(*,*)
    write(*,*) 'finished estimating the horizontal localization radii...'

  end subroutine calcHorizLocalizationRadii

  !--------------------------------------------------------------------------
  ! FINDBININDEX
  !--------------------------------------------------------------------------
  function findBinIndex(distance,distanceBinThresholds,numBins) result(binIndex)
    implicit none

    integer :: numBins, binIndex
    real(8) :: distance
    real(8) :: distanceBinThresholds(numBins)

    integer :: bin

    binIndex = -1
    do bin = 1, numbins
      if ( distance <= distanceBinThresholds(bin) ) then
        binIndex = bin
        exit
      end if
    end do

    if (binIndex == -1) then
      write(*,*) 'findBinIndex: No match found! ABORTING'
      call utl_abort('findBinIndex')
    end if

  end function findBinIndex

  !--------------------------------------------------------------------------
  ! DISTANCE
  !--------------------------------------------------------------------------
  function calcDistance(lat2, lon2, lat1, lon1) result(distanceInM)
    !:Purpose: To compute the distance between two points on Earth: (lat1,lon1)
    !          and (lat2,lon2). Calcul utilisant la Formule d'Haversine
    !          Reference: R.W. Sinnott,'Virtues of Haversine',Sky and Telescope,
    !          vol.68, no.2, 1984, p.159)
    implicit none
    real(8) :: distanceInM

    ! Arguments:
    real(8) :: lat1, lon1, lat2, lon2

    ! Locals:
    real(8) :: dlat, dlon, a, c

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = (sin(dlat/2.d0))**2 + cos(lat1)*cos(lat2)*(sin(dlon/2.d0))**2
    c = 2.d0 * atan2(sqrt(a),sqrt(1.d0-a))
    distanceInM = EC_RA * c

  end function calcDistance

  !--------------------------------------------------------------------------
  ! CALCVERTLOCALIZATIONRADII
  !--------------------------------------------------------------------------
  subroutine calcVertLocalizationRadii(ensPerts,stride,waveBandIndex_opt)
    implicit none

    type(struct_ens)    :: ensPerts
    integer, intent(in) :: stride
    integer,optional, intent(in) :: waveBandIndex_opt

    type(struct_gsv) :: statevector_ensStdDev

    real(8), pointer :: ensStdDev(:,:,:)
    real(4), pointer :: ptr4d_r4(:,:,:,:)

    real*4, allocatable :: ensPert_local(:,:)

    real(8), allocatable :: meanCorrel(:,:)
    real(8), allocatable :: meanCorrelSquare(:,:)
    real(8), allocatable :: meanVarianceProduct(:,:)
    real(8), allocatable :: meanCovarianceSquare(:,:)
    real(8), allocatable :: meanFourthMoment(:,:)

    real(8), allocatable :: meanCorrel_local(:,:)
    real(8), allocatable :: meanCorrelSquare_local(:,:)
    real(8), allocatable :: meanVarianceProduct_local(:,:)
    real(8), allocatable :: meanCovarianceSquare_local(:,:)
    real(8), allocatable :: meanFourthMoment_local(:,:)

    real(8), allocatable :: localizationFunctions(:,:,:) ! Eq. 19-21 in MMMB 2015 Part 2

    real(8), allocatable :: localizationRadii(:,:)

    real(8), allocatable, target :: distanceBinInLnP_T(:,:) ! Distance between each pair of thermo vertical levels in ln(Pressure)
    real(8), allocatable, target :: distanceBinInLnP_M(:,:) ! Distance between each pair of momentum vertical levels in ln(Pressure)
    real(8), allocatable :: distanceBinWeight(:)    ! Weight given to each bin in the curve fitting step

    real(8), allocatable :: gridPointWeight(:,:)   ! Weight given to grid point in the statistic computation

    real(8), allocatable :: sumWeight(:,:)    ! Sample size for each distance-bin
    real(8), allocatable :: sumWeight_local(:,:)

    real(8), pointer :: PressureProfile(:)
    real(8), pointer :: distanceBinInLnP(:,:)

    real(8) :: dnens, correlation, covariance, fourthMoment, weight
    real(8) :: t1, t2, t3, rmse

    integer :: k, k2, kens, f, ens, bin, numbins, numFunctions
    integer :: iref, jref
    integer :: nLevEns, nLevStart, nLevEnd, jvar
    integer :: myLonBeg, myLonEnd, myLatBeg, myLatEnd, nSize, ier

    character(len=128) :: outfilename
    character(len=2)   :: wbnum

    !
    !- 1.  Setup
    !
    call ens_copyEnsStdDev(ensPerts, statevector_ensStdDev)
    call gsv_getField(statevector_ensStdDev,ensStdDev)

    numFunctions = 3

    numBins=max(nLevEns_M,nLevEns_T)

    call ens_getLatLonBounds(ensPerts, myLonBeg, myLonEnd, myLatBeg, myLatEnd)

    allocate(meanCorrelSquare(numBins,nkgdimEns))
    allocate(meanCorrel(numBins,nkgdimEns))
    allocate(meanVarianceProduct(numBins,nkgdimEns))
    allocate(meanCovarianceSquare(numBins,nkgdimEns))
    allocate(meanFourthMoment(numBins,nkgdimEns))
    allocate(sumWeight(numBins,nkgdimEns))
    
    allocate(meanCorrelSquare_local(numBins,nkgdimEns))
    allocate(meanCorrel_local(numBins,nkgdimEns))
    allocate(meanVarianceProduct_local(numBins,nkgdimEns))
    allocate(meanCovarianceSquare_local(numBins,nkgdimEns))
    allocate(meanFourthMoment_local(numBins,nkgdimEns))
    allocate(sumWeight_local(numBins,nkgdimEns))

    allocate(gridPointWeight(ni,nj))

    allocate(distanceBinInLnP_M(nLevEns_M,nLevEns_M))
    allocate(distanceBinInLnP_T(nLevEns_T,nLevEns_T))

    ! Compute the separation-distance in ln(p) between vertical levels
    do k = 1, nLevEns_M
      do k2 = 1, nLevEns_M
        distanceBinInLnP_M(k,k2) = abs(log(pressureProfile_M(k))-log(pressureProfile_M(k2)))
      end do
    end do
    do k = 1, nLevEns_T
      do k2 = 1, nLevEns_T
        distanceBinInLnP_T(k,k2) = abs(log(pressureProfile_T(k))-log(pressureProfile_T(k2)))
      end do
    end do

    dnens = 1.0d0/dble(nens-1)

    ! Grid point Weight
    do jref = 1, nj
      gridPointWeight(:,jref) = cos(hco_ens%lat(jref))
    end do

    !
    !- 2.  Estimation of localization functions
    !

    !- 2.1  Computation of various statistics for different separation distances
    meanCorrelSquare_local(:,:)     = 0.d0
    meanCorrel_local(:,:)           = 0.d0
    meanCovarianceSquare_local(:,:) = 0.d0
    meanVarianceProduct_local(:,:)  = 0.d0
    meanFourthMoment_local(:,:)     = 0.d0
    sumWeight_local(:,:)            = 0.d0

    allocate(ensPert_local(nens,nkgdimEns))

    do jref = nint(stride/2.0)+horizPadding, nj-horizPadding, stride    ! Pick every stride point to save cost.
      do iref = nint(stride/2.0)+horizPadding, ni-horizPadding, stride  ! Pick every stride point to save cost.

        if (iref < myLonBeg .or. iref > myLonEnd .or. &
            jref < myLatBeg .or. jref > myLatEnd ) cycle

        !- Select data needed to speed up the process (ensemble member index must come first in ensPert_Local 
        !  because ensemble member is the last loop index below)
        do ens = 1, nens
          do k = 1, nkgdimEns
            ptr4d_r4 => ens_getOneLev_r4(ensPerts,k)
            ensPert_local(ens,k) = ptr4d_r4(ens,1,iref,jref)
          end do
        end do

        ! Loop on all 3D variables
        do jvar = 1, nvar3d

          if (vnl_varLevelFromVarName(nomvar3d(jvar)).eq.'MM') then
            nLevEns = nLevEns_M
          else
            nLevEns = nLevEns_T
          endif
          nLevStart = varLevOffset(jvar)+ 1
          nLevEnd   = varLevOffset(jvar)+ nLevEns

          !$OMP PARALLEL DO PRIVATE (k,k2,ens,correlation,covariance,fourthMoment,bin,weight)
          do k = nLevStart, nLevEnd
            do k2 = nLevStart, nLevEnd

              covariance = 0.d0
              fourthMoment = 0.d0
              do ens = 1, nens
                covariance = covariance + &
                     ensPert_local(ens,k2) * ensPert_local(ens,k)
                fourthMoment = fourthMoment + &
                     (ensPert_local(ens,k2) * ensPert_local(ens,k))**2
              end do
              covariance   = covariance * dnens
              fourthMoment = fourthMoment / real(nens,8)
              if ( ensStdDev(iref,jref,k2) > 0.d0 .and. ensStdDev(iref,jref,k) > 0.d0 ) then
                correlation  = covariance / (ensStdDev(iref,jref,k2)*ensStdDev(iref,jref,k))
              else
                correlation  = 0.d0
              end if

              bin=k2-nLevStart+1

              weight = gridPointWeight(iref,jref)

              sumWeight_local(bin,k) = sumWeight_local(bin,k) + weight

              meanCorrel_local(bin,k)           = meanCorrel_local(bin,k)           + &
                                                  correlation    * weight
              meanCorrelSquare_local(bin,k)     = meanCorrelSquare_local(bin,k)     + &
                                                  correlation**2 * weight
              meanCovarianceSquare_local(bin,k) = meanCovarianceSquare_local(bin,k) + &
                                                  covariance**2  * weight
              meanVarianceProduct_local(bin,k)  = meanVarianceProduct_local(bin,k)  + &
                                                  ensStdDev(iref,jref,k2)**2 * ensStdDev(iref,jref,k)**2 * weight
              meanFourthMoment_local(bin,k)     = meanFourthMoment_local(bin,k)     + &
                                                  fourthMoment   * weight
            end do
          end do
          !$OMP END PARALLEL DO
        end do ! var3D

      end do ! iref
    end do ! jref

    deallocate(ensPert_local)

    !- 2.2 Gather the all the info in processor 0
    nSize = nkgdimEns * numbins
    call rpn_comm_reduce(sumWeight_local           ,sumWeight           ,nSize, &
         "mpi_double_precision","mpi_sum",0,"GRID",ier)
    call rpn_comm_reduce(meanCorrel_local          ,meanCorrel          ,nSize, &
         "mpi_double_precision","mpi_sum",0,"GRID",ier)
    call rpn_comm_reduce(meanCorrelSquare_local    ,meanCorrelSquare    ,nSize, &
         "mpi_double_precision","mpi_sum",0,"GRID",ier)
    call rpn_comm_reduce(meanVarianceProduct_local ,meanVarianceProduct ,nSize, &
         "mpi_double_precision","mpi_sum",0,"GRID",ier)
    call rpn_comm_reduce(meanFourthMoment_local    ,meanFourthMoment    ,nSize, &
         "mpi_double_precision","mpi_sum",0,"GRID",ier)
    call rpn_comm_reduce(meanCovarianceSquare_local,meanCovarianceSquare,nSize, &
         "mpi_double_precision","mpi_sum",0,"GRID",ier)

    deallocate(sumWeight_local)
    deallocate(meanCorrel_local)
    deallocate(meanCorrelSquare_local)
    deallocate(meanVarianceProduct_local)
    deallocate(meanFourthMoment_local)
    deallocate(meanCovarianceSquare_local)

    if (mpi_myid == 0) then
      !- 2.3  Computation of the localization functions
      allocate(localizationFunctions(numFunctions,numBins,nkgdimEns))
      
      t1=dble((nens-1)**2)/dble(nens*(nens-3))
      t2=dble(nens)/dble((nens-2)*(nens-3))
      t3=dble(nens-1)/dble(nens*(nens-2)*(nens-3))

      do jvar = 1, nvar3d
        if (vnl_varLevelFromVarName(nomvar3d(jvar)).eq.'MM') then
          nLevEns = nLevEns_M
        else
          nLevEns = nLevEns_T
        endif
        nLevStart = varLevOffset(jvar)+ 1
        nLevEnd   = varLevOffset(jvar)+ nLevEns
        !$OMP PARALLEL DO PRIVATE (k,bin)       
        do k = nLevStart, nLevEnd
          do bin = 1, nLevEns
            meanCorrel(bin,k)           = meanCorrel(bin,k)           / sumWeight(bin,k)
            meanCorrelSquare(bin,k)     = meanCorrelSquare(bin,k)     / sumWeight(bin,k)
            meanCovarianceSquare(bin,k) = meanCovarianceSquare(bin,k) / sumWeight(bin,k)
            meanVarianceProduct(bin,k)  = meanVarianceProduct(bin,k)  / sumWeight(bin,k)
            meanFourthMoment(bin,k)     = meanFourthMoment(bin,k)     / sumWeight(bin,k)

            if ( meanCovarianceSquare(bin,k) /= 0.d0 ) then
              ! Form 1: General formulation (Eq. 19 in MMMB 2015 Part 2)
              localizationFunctions(1,bin,k) = t1 - t2*meanFourthMoment(bin,k)/meanCovarianceSquare(bin,k) + &
                   t3*meanVarianceProduct(bin,k)/meanCovarianceSquare(bin,k)
              ! Form 2: Gaussian sample distribution (Eq. 20 in MMMB 2015 Part 2)
              localizationFunctions(2,bin,k) = dble(nens-1)/dble((nens+1)*(nens-2)) * &
                   (dble(nens-1)-meanVarianceProduct(bin,k)/meanCovarianceSquare(bin,k))
            else
              write(*,*) ' !!! Warning !!! meanCovarianceSquare = 0 in bin, level = ', bin, k
              localizationFunctions(1,bin,k) = 0.d0
              localizationFunctions(2,bin,k) = 0.d0
            end if
            ! Form 3: Gaussian sample distribution and correlation-based formulation (Eq. 21 in MMMB 2015 Part 2)
            if ( meanCorrelSquare(bin,k) /= 0.d0 ) then
              localizationFunctions(3,bin,k) = dble(nens-1)/dble((nens+1)*(nens-2)) * &
                   (dble(nens-1)-1.d0/meanCorrelSquare(bin,k))
            else
              write(*,*) ' !!! Warning !!! meanCorrelSquare = 0 in bin, level = ', bin, k
              localizationFunctions(3,bin,k) = 0.d0
            end if
          end do
        end do
        !$OMP END PARALLEL DO
      end do

      !
      !- 3.  Estimation of localization radii (curve fitting)
      !
      allocate(localizationRadii(numFunctions,nkgdimEns))
      allocate(distanceBinWeight(numBins))
      
      call lfn_setup('FifthOrder') ! IN
      
      localizationRadii(:,:) = 2.d0 ! First Guess (in ln(p) distance)
      distanceBinWeight(:)   = 1.d0 ! Even weight

      do jvar = 1, nvar3d
        if (vnl_varLevelFromVarName(nomvar3d(jvar)).eq.'MM') then
          nLevEns = nLevEns_M
          distanceBinInLnP => distanceBinInLnP_M
        else
          nLevEns = nLevEns_T
          distanceBinInLnP => distanceBinInLnP_T
        endif
        nLevStart = varLevOffset(jvar)+ 1
        nLevEnd   = varLevOffset(jvar)+ nLevEns
        
        write(*,*)
        write(*,*) nomvar3d(jvar)
        
        do f = 1, numFunctions
          
          write(*,*)
          write(*,*) 'Localization Function : ', f
          do k =  nLevStart, nLevEnd
            kens = k-nLevStart+1
            write(*,*) '         ----- EnsLev : ', k
            call lfn_lengthscale( localizationRadii(f,k),       & ! INOUT
                                  rmse,                         & ! OUT
                                  localizationFunctions(f,1:nLevEns,k), & ! IN
                                  distanceBinInLnP(kens,:),     & ! IN
                                  distanceBinWeight(1:nLevEns), & ! IN
                                  nLevEns )                       ! IN
          end do
        end do
      end do

      !
      !- 4.  Write to file
      !
      if ( nWaveBand /= 1 ) then
        if (.not. present(waveBandIndex_opt)) then
          write(*,*) 'calcLocalizationRadii: No waveBandIndex was supplied!!!'
          call utl_abort('calbmatrix_glb')
        end if
        write(wbnum,'(I2.2)') waveBandIndex_opt
      end if
      
      !- 4.1 Localization functions in txt format (for plotting purposes)
      do jvar = 1, nvar3d
        if ( nWaveBand == 1 ) then
          outfilename = "./vertLocalizationFunctions_"//trim(nomvar3d(jvar))//".txt"
        else
          outfilename = "./vertLocalizationFunctions_"//trim(nomvar3d(jvar))//"_"//wbnum//".txt"
        end if
        open (unit=99,file=outfilename,action="write",status="new")
        
        if(vnl_varLevelFromVarName(nomvar3d(jvar)).eq.'MM') then
          nLevEns = nLevEns_M
          PressureProfile => pressureProfile_M
          distanceBinInLnP => distanceBinInLnP_M 
        else
          nLevEns = nLevEns_T
          PressureProfile => pressureProfile_T
          distanceBinInLnP => distanceBinInLnP_T
        endif
        
        do k=1,nlevEns
          do bin = 1, nlevEns
            write(99,'(I3,2X,I3,2X,F7.2,2X,F7.3,2X,I7,2X,F6.4,2X,F6.4,2X,F6.4,2X,F6.4,2X,F6.4,2X,E10.3,2X,E10.3,2X,E10.3)') k, bin, &
                 PressureProfile(bin)/100.d0, distanceBinInLnP(k,bin), nint(sumWeight(bin,varLevOffset(jvar)+k)), &
                 meanCorrel(bin,varLevOffset(jvar)+k), &
                 localizationFunctions(3,bin,varLevOffset(jvar)+k), &
                 localizationFunctions(2,bin,varLevOffset(jvar)+k), &
                 localizationFunctions(1,bin,varLevOffset(jvar)+k), &
                 meanCorrelSquare(bin,varLevOffset(jvar)+k), &
                 meanCovarianceSquare(bin,varLevOffset(jvar)+k), &
                 meanVarianceProduct(bin,varLevOffset(jvar)+k), &
                 meanFourthMoment(bin,varLevOffset(jvar)+k)
          end do
        end do
        close(unit=99)
      end do
      
      !- 4.2 Localization radii in txt format (for plotting purposes)
      do jvar = 1, nvar3d
        if ( nWaveBand == 1 ) then
          outfilename = "./vertLocalizationRadii_"//trim(nomvar3d(jvar))//".txt"
        else
          outfilename = "./vertLocalizationRadii_"//trim(nomvar3d(jvar))//"_"//wbnum//".txt"
        end if
        open (unit=99,file=outfilename,action="write",status="new")
        
        if(vnl_varLevelFromVarName(nomvar3d(jvar)).eq.'MM') then
          nLevEns = nLevEns_M
          PressureProfile => pressureProfile_M
        else
          nLevEns = nLevEns_T
          PressureProfile => pressureProfile_T
        endif
        do k=1,nlevEns
          write(99,'(I3,2X,F7.2,2X,F7.2,2X,F7.2,2X,F7.2)') k, PressureProfile(k)/100.d0, &
               localizationRadii(3,varLevOffset(jvar)+k), &
               localizationRadii(2,varLevOffset(jvar)+k), &
               localizationRadii(1,varLevOffset(jvar)+k)
        end do
        close(unit=99)
      end do

      deallocate(localizationFunctions)
      deallocate(localizationRadii)
      deallocate(distanceBinWeight)

    end if ! mpi_myid == 0

    deallocate(meanCorrelSquare)
    deallocate(meanCorrel)
    deallocate(meanCovarianceSquare)
    deallocate(meanVarianceProduct)
    deallocate(meanFourthMoment)
    deallocate(distanceBinInLnP_M)
    deallocate(distanceBinInLnP_T)
    deallocate(sumWeight)
    deallocate(gridPointWeight)

    write(*,*) 'finished estimating the vertical localization radii...'

  end subroutine calcVertLocalizationRadii

end module menetrierDiag_mod
