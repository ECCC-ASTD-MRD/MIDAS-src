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

!--------------------------------------------------------------------------
!!
!! *Purpose*: Main program for computing diagnostics of the B and L matrices
!!
!--------------------------------------------------------------------------
program main_diagBmatrix
  use mpivar_mod
  use MathPhysConstants_mod
  use controlVector_mod
  use variableTransforms_mod
  use gridStateVector_mod
  use bmatrix_mod
  use bmatrixEnsemble_mod
  use localization_mod
  use horizontalCoord_mod
  use analysisGrid_mod
  use verticalCoord_mod
  use timeCoord_mod
  use randomNumber_mod
  use utilities_mod
  use ramDisk_mod
  use globalSpectralTransform_mod
  IMPLICIT NONE

  type(struct_gsv) :: statevector, statevectorEnsAmp
  type(struct_hco), pointer :: hco_anl => null()
  type(struct_hco), pointer :: hco_core => null()
  type(struct_vco), pointer :: vco_anl => null()
  type(struct_loc), pointer :: locInfo => null()

  real(8), pointer :: field(:,:,:)
  real(8), allocatable :: ensAmplitude(:,:,:,:)
  real(4), allocatable :: randomEns(:,:,:,:)
  real(8), allocatable :: mean(:,:,:)
  real(8), allocatable :: stddev(:,:,:)
  real(8), allocatable :: stddev_zm(:,:),stddev_zm2(:,:)
  real(8), allocatable :: stddev_dm(:,:),stddev_dm2(:,:)
  real(8), allocatable :: zonalMeanStddev(:)
  real(8), allocatable :: controlVector(:), controlVector_global(:)

  integer :: fclos, fnom, fstopc, newdate, get_max_rss
  integer :: ierr, nsize, iseed, cvDim_local
  integer :: ensIndex, index, lonIndex, latIndex, kIndex, nkgdim, levIndex
  integer :: idate, itime, nulnam, nultxt, dateStamp, numLoc
  integer :: nlons, nlats, nlevs, nlevs2, jvar, ip3
  integer :: latIndex2, lonIndex2, levIndex2

  integer :: latPerPE, lonPerPE
  integer :: myLatBeg, myLatEnd
  integer :: myLonBeg, myLonEnd

  character(len=128) :: filename, filenameEnsAmp
  character(len=10)  :: datestr
  character(len=12)  :: etiket
  character(len=4)   :: varName

  ! namelist variables
  integer :: numperturbations, nrandseed, diagdate
  integer :: oneobs_levs(100),oneobs_lons(100),oneobs_lats(100)
  logical :: writeEnsAmplitude
  logical :: writeTextStddev
  logical :: writePsiChiStddev

  namelist /namdiag/numperturbations, nrandseed, diagdate, oneobs_levs, oneobs_lons, oneobs_lats, &
                    writeEnsAmplitude, writeTextStddev, writePsiChiStddev

  write(*,*) " -------------------------------------------"
  write(*,*) " --- START OF MAIN PROGRAM diagBmatrix   ---"
  write(*,*) " --- Diagnositcs of the B matrix         ---"
  write(*,*) " -------------------------------------------"

  ! MPI, tmg initialization
  call mpi_initialize 
  call tmg_init(mpi_myid, 'TMG_DIAGBMATRIX' )
  call tmg_start(1,'MAIN')
  ierr = fstopc('MSGLVL','ERRORS',0)

  ! Set default values for namelist NAMDIAG parameters
  diagdate = 2011020100
  numperturbations = -1
  nrandseed = 1
  oneobs_levs(:)=-1
  oneobs_lons(:)=-1
  oneobs_lats(:)=-1
  writeEnsAmplitude = .false.
  writeTextStddev = .false.
  writePsiChiStddev = .false.

  ! Read the parameters from NAMDIAG
  nulnam=0
  ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
  read(nulnam,nml=namdiag,iostat=ierr)
  if(ierr.ne.0) call utl_abort('diagBmatrix: Error reading namelist')
  write(*,nml=namdiag)
  ierr=fclos(nulnam)

  nlevs=0
  do index = 1, size(oneobs_levs)
    if(oneobs_levs(index).ge.1) nlevs=nlevs+1  
  end do
  nlons=0
  do index = 1, size(oneobs_lons)
    if(oneobs_lons(index).ge.1) nlons=nlons+1  
  end do
  nlats=0
  do index = 1, size(oneobs_lats)
    if(oneobs_lats(index).ge.1) nlats=nlats+1  
  end do

  ! Decompose diagdate(yyyymmddhh) into idate(YYYYMMDD) itime(HHMMSShh)
  ! and calculate date-time stamp
  idate = diagdate/100
  itime = (diagdate-idate*100)*1000000
  ierr = newdate(dateStamp,idate,itime,3)
  write(datestr,'(i10.10)') diagdate
  write(*,*)' idate= ',idate,' time= ',itime
  write(*,*)' date= ',diagdate,' stamp= ',dateStamp

  ! Top Level Control setup
  call ram_setup

  !- Initialize the Temporal grid
  call tim_setup
  call tim_setDatestamp(dateStamp)

  ! Initialize variables of the model states
  call gsv_setup

  ! Initialize the Analysis horizontal grid
  call hco_SetupFromFile( hco_anl,'./analysisgrid', 'ANALYSIS', 'Analysis' ) ! IN

  if ( hco_anl % global ) then
    call agd_SetupFromHCO( hco_anl ) ! IN
  else
    !- Iniatilized the core (Non-Exteded) analysis grid
    call hco_SetupFromFile( hco_core, './analysisgrid', 'COREGRID', 'AnalysisCore' ) ! IN
    !- Setup the LAM analysis grid metrics
    call agd_SetupFromHCO( hco_anl, hco_core ) ! IN
  end if

  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! Initialize the vertical coordinate from the statistics file
  call vco_SetupFromFile( vco_anl,        & ! OUT
                          './analysisgrid') ! IN

  ! Allocate the statevector
  call gsv_allocate(statevector, tim_nstepobsinc, hco_anl, vco_anl, &
                    datestamp=tim_getDatestamp(), mpi_local=.true.)
  call gsv_zero(statevector)
  nkgdim = statevector%nk

  ! Setup the B matrix
  call bmat_setup(hco_anl,vco_anl)

  ! Setup of the L matrix done in bmat_setup

  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !==============================================
  !- Compute columns of B and L matrices
  !==============================================
  !
  if ( nlevs.ge.1 .and. nlons.ge.1 .and. nlats.ge.1 ) then

    !
    !- Compute columns of the B matrix
    !
    allocate(controlVector(cvm_nvadim))

    write(*,*) '********************************************'
    write(*,*) 'Compute columns of B matrix'
    write(*,*) '********************************************'

    write(*,*) 'number of levels     =',nlevs
    write(*,*) 'number of longitudes =',nlons
    write(*,*) 'number of latitudes  =',nlats

    do jvar = 1, vnl_numvarmax

      if (.not. gsv_varExist(varName=vnl_varNameList(jvar)) ) cycle

      filename       = 'columnB_'      // trim(vnl_varNameList(jvar)) // '_' // datestr // '.fst'
      filenameEnsAmp = 'ensAmplitude_' // trim(vnl_varNameList(jvar)) // '_' // datestr

      if(vnl_varLevelFromVarname(vnl_varNameList(jvar)).eq.'SF') then
        nlevs2 = 1
      else
        nlevs2 = nlevs
      end if

      ip3 = 0
      do levIndex = 1, nlevs2
        do lonIndex = 1, nlons
          do latIndex = 1, nlats

            call gsv_zero(statevector)                   
            field => gsv_getField3d_r8(statevector,vnl_varNameList(jvar))

            if(oneobs_lats(latIndex).ge.statevector%myLatBeg .and. oneobs_lats(latIndex).le.statevector%myLatEnd .and.  &
                 oneobs_lons(lonIndex).ge.statevector%myLonBeg .and. oneobs_lons(lonIndex).le.statevector%myLonEnd) then
              if(vnl_varLevelFromVarname(vnl_varNameList(jvar)).eq.'SF') then
                field(oneobs_lons(lonIndex),oneobs_lats(latIndex),1) = 1.0D0
              else
                field(oneobs_lons(lonIndex),oneobs_lats(latIndex),oneobs_levs(levIndex)) = 1.0D0
              end if
            end if

            controlVector(:)=0.0d0
            call bmat_sqrtBT(controlVector,cvm_nvadim,statevector)
            call bmat_sqrtB (controlVector,cvm_nvadim,statevector)

            write(*,*)'diagBmatrix: writing out the column of B, levIndex,lonIndex,latIndex=',levIndex,lonIndex,latIndex
            call flush(6)

            ip3 = ip3 + 1

            call gsv_writeToFile(statevector,filename,'ONEOBS_'//trim(vnl_varNameList(jvar)),  &
                 ip3_in=ip3,HUcontainsLQ=.true.,unitConversion=.true.)

            ! Write the ensemble amplitude fields (i.e., the alphas) when Bens is active
            if (writeEnsAmplitude) call ben_writeAmplitude('./',filenameEnsAmp, ip3) ! IN

          end do
        end do
      end do

    end do

    deallocate(controlVector)

    !
    !- Compute columns of the L matrix
    !
    numLoc = loc_getNumLocActive()

    if (numLoc /= 0) then
      locInfo => loc_getLocInfo(1) ! Grab the first one...

      call mpivar_setup_latbands(locInfo%hco%nj,latPerPE,myLatBeg,myLatEnd)
      call mpivar_setup_lonbands(locInfo%hco%ni,lonPerPE,myLonBeg,myLonEnd)

      call gsv_allocate(statevectorEnsAmp, 1, locInfo%hco, locInfo%vco, &
                        datestamp=tim_getDatestamp(),mpi_local=.true.,varName='ALFA')

      allocate(ensAmplitude(locInfo%nEnsOverDimension,myLonBeg:myLonEnd,myLatBeg:myLatEnd,locInfo%vco%nLev_M))

      allocate(controlVector(locInfo%cvDim))

      write(*,*) '********************************************'
      write(*,*) 'Compute columns of L matrix'
      write(*,*) '********************************************'
      
      write(*,*) 'number of levels     =',nlevs
      write(*,*) 'number of longitudes =',nlons
      write(*,*) 'number of latitudes  =',nlats
      
      ip3 = 0
      filename = 'columnL_' // trim(locInfo%locType) // '_' // datestr // '.fst'

      do levIndex = 1, nlevs2
        do lonIndex = 1, nlons
          do latIndex = 1, nlats

            ensAmplitude(:,:,:,:) = 0.d0
            if(oneobs_lats(latIndex).ge.myLatBeg .and. oneobs_lats(latIndex).le.myLatEnd .and.  &
                 oneobs_lons(lonIndex).ge.myLonBeg .and. oneobs_lons(lonIndex).le.myLonEnd) then
              if ( locInfo%hco%global ) then
                ensAmplitude(1,oneobs_lons(lonIndex),oneobs_lats(latIndex),oneobs_levs(levIndex)) = &
                     1.d0 * real(locInfo%hco%ni,8) / gst_getRWT(oneobs_lats(latIndex)) ! This normalization should not be done here
              else
                ensAmplitude(1,oneobs_lons(lonIndex),oneobs_lats(latIndex),oneobs_levs(levIndex)) = 1.d0
              end if
            end if
            controlVector(:)=0.0d0

            call loc_LsqrtAd(1,             & ! IN
                             ensAmplitude,  & ! IN
                             controlVector)   ! OUT
            call loc_Lsqrt  (1,             & ! IN
                             controlVector, & ! IN
                             ensAmplitude)    ! OUT

            write(*,*)'diagBmatrix: writing out the column of L, levIndex,lonIndex,latIndex=',levIndex,lonIndex,latIndex
            call flush(6)
            
            ip3 = ip3 + 1
            
            nullify(field)
            field => gsv_getField3D_r8(statevectorEnsAmp,'ALFA')
!$OMP PARALLEL DO PRIVATE(latIndex2,lonIndex2,levIndex2)
            do levIndex2 = 1, locInfo%vco%nLev_M
              do latIndex2 = myLatBeg, myLatEnd
                do lonIndex2 = myLonBeg, myLonEnd
                  field(lonIndex2,latIndex2,levIndex2) = &
                      ensAmplitude(1,lonIndex2,latIndex2,levIndex2)
                end do
              end do
            end do
!$OMP END PARALLEL DO

            call gsv_writeToFile(statevectorEnsAmp,filename,'ONEOBS',  &
                 ip3_in=ip3,unitConversion=.false.)

          end do
        end do
      end do

      deallocate(controlVector)
      deallocate(ensAmplitude)
      call gsv_deallocate(statevectorEnsAmp)

    end if ! if localization is active in B

  end if ! if any oneobs selected

  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  !
  !==============================================
  !- Compute the stddev from random perturbations
  !==============================================
  !
  if ( numperturbations > 1 ) then

    write(*,*) '********************************************'
    write(*,*) 'Compute the stddev from random perturbations'
    write(*,*) '********************************************'

    allocate(controlVector(cvm_nvadim))
    allocate(controlVector_global(cvm_nvadim_mpiglobal))

    call gsv_zero(statevector)

    ! Allocate the randomEns, mean and stddev
    allocate(randomEns(statevector%myLonBeg:statevector%myLonEnd,statevector%myLatBeg:statevector%myLatEnd,nkgdim,numperturbations))
    allocate(mean(statevector%myLonBeg:statevector%myLonEnd,statevector%myLatBeg:statevector%myLatEnd,nkgdim))
    allocate(stddev(statevector%myLonBeg:statevector%myLonEnd,statevector%myLatBeg:statevector%myLatEnd,nkgdim))

    iseed = abs(nrandseed)
    call rng_setup(iseed)

    field => gsv_getField3d_r8(statevector)

    !
    !- Compute the ensemble of random perturbations
    !
    do ensIndex = 1, numperturbations
      write(*,*) ' computing member number= ',ensIndex
      call flush(6)

      !- Global vector (same for each processors)
      do index = 1, cvm_nvadim_mpiglobal
        controlVector_global(index) = rng_gaussian()
      end do

      !- Extract only the subvector for this processor
      call bmat_reduceToMPILocal(controlVector,        & ! OUT
           controlVector_global, & ! IN
           cvDim_local )           ! OUT

      !- Transform to control variables in physical space
      call bmat_sqrtB(controlVector,cvm_nvadim,statevector)

      if ( writePsiChiStddev ) call vtr_transform(statevector,'UVtoPsiChi')

      !$OMP PARALLEL DO PRIVATE (lonIndex,latIndex,kIndex)    
      do kIndex = 1, nkgdim
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do lonIndex = statevector%myLonBeg, statevector%myLonEnd
            randomEns(lonIndex,latIndex,kIndex,ensIndex) = field(lonIndex,latIndex,kIndex)
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    end do ! Loop on ens member

    !
    !- Compute the ensemble mean
    !
    mean(:,:,:) = 0.0d0
    do ensIndex = 1, numperturbations
      !$OMP PARALLEL DO PRIVATE (lonIndex,latIndex,kIndex)
      do kIndex = 1, nkgdim
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do lonIndex = statevector%myLonBeg, statevector%myLonEnd
            mean(lonIndex,latIndex,kIndex) = mean(lonIndex,latIndex,kIndex) + randomEns(lonIndex,latIndex,kIndex,ensIndex)
          end do
        end do
      end do
      !$OMP END PARALLEL DO
    end do

    !$OMP PARALLEL DO PRIVATE (lonIndex,latIndex,kIndex)    
    do kIndex = 1, nkgdim
      do latIndex = statevector%myLatBeg, statevector%myLatEnd
        do lonIndex = statevector%myLonBeg, statevector%myLonEnd
          mean(lonIndex,latIndex,kIndex) = mean(lonIndex,latIndex,kIndex)/real(numperturbations,8)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    !
    !- Remove the ensemble mean from the ensemble
    !
    !$OMP PARALLEL DO PRIVATE (lonIndex,ensIndex,latIndex,kIndex)    
    do ensIndex = 1, numperturbations
      do kIndex = 1, nkgdim
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do lonIndex = statevector%myLonBeg, statevector%myLonEnd
            randomEns(lonIndex,latIndex,kIndex,ensIndex) = randomEns(lonIndex,latIndex,kIndex,ensIndex) - mean(lonIndex,latIndex,kIndex)
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    deallocate(mean)

    !
    !- Compute the ensemble stddev
    !
    stddev(:,:,:) = 0.0d0

    do ensIndex = 1, numperturbations
      !$OMP PARALLEL DO PRIVATE (lonIndex,latIndex,kIndex)
      do lonIndex = statevector%myLonBeg, statevector%myLonEnd
        do latIndex = statevector%myLatBeg, statevector%myLatEnd
          do kIndex = 1, nkgdim
            stddev(lonIndex,latIndex,kIndex) = stddev(lonIndex,latIndex,kIndex) + &
                 (randomEns(lonIndex,latIndex,kIndex,ensIndex)**2)/real(numperturbations,8)
          end do
        end do
      end do
      !$OMP END PARALLEL DO
    end do
    deallocate(randomEns)

    !$OMP PARALLEL DO PRIVATE (lonIndex,latIndex,kIndex)
    do kIndex = 1, nkgdim
      do latIndex = statevector%myLatBeg, statevector%myLatEnd
        do lonIndex = statevector%myLonBeg, statevector%myLonEnd
          stddev(lonIndex,latIndex,kIndex) = sqrt(stddev(lonIndex,latIndex,kIndex))
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    !- Insert results in statevector
    !$OMP PARALLEL DO PRIVATE (lonIndex,latIndex,kIndex)    
    do kIndex = 1, nkgdim
      do latIndex = statevector%myLatBeg, statevector%myLatEnd
        do lonIndex = statevector%myLonBeg, statevector%myLonEnd
          field(lonIndex,latIndex,kIndex) = stddev(lonIndex,latIndex,kIndex)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    deallocate(stddev)

    !- Write to file
    call gsv_writeToFile(statevector,'stddev_' // datestr // '.fst','GD_STDDEV',  &
         HUcontainsLQ=.true.,unitConversion=.true.)

    !
    !- Compute the zonal mean std dev
    !
    write(*,*) 'Compute the zonal mean stddev'
    call flush(6)

    allocate(stddev_zm(hco_anl%nj,nkgdim))
    allocate(stddev_zm2(hco_anl%nj,nkgdim))
    stddev_zm(:,:) = 0.d0
    stddev_zm2(:,:) = 0.d0

    !$OMP PARALLEL DO PRIVATE (lonIndex,latIndex,kIndex)
    do kIndex = 1, nkgdim
      do latIndex = statevector%myLatBeg, statevector%myLatEnd
        do lonIndex = statevector%myLonBeg, statevector%myLonEnd
          stddev_zm(latIndex,kIndex) = stddev_zm(latIndex,kIndex) + (field(lonIndex,latIndex,kIndex)**2)/real(hco_anl%ni,8)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    nsize = statevector%nj*nkgdim
    call rpn_comm_allreduce(stddev_zm,stddev_zm2,nsize,  &
         "MPI_DOUBLE_PRECISION","MPI_SUM","GRID",ierr)

    !- Insert results in statevector
    !$OMP PARALLEL DO PRIVATE (lonIndex,latIndex,kIndex)    
    do kIndex = 1, nkgdim
      do latIndex = statevector%myLatBeg, statevector%myLatEnd
        do lonIndex = statevector%myLonBeg, statevector%myLonEnd
          if(stddev_zm2(latIndex,kIndex)> 0.0d0) then
            field(lonIndex,latIndex,kIndex) = sqrt(stddev_zm2(latIndex,kIndex))
          else
            field(lonIndex,latIndex,kIndex) = 0.0d0
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    deallocate(stddev_zm)
    deallocate(stddev_zm2)

    call gsv_writeToFile(statevector,'stddev_' // datestr // '.fst','ZM_STDDEV',  &
                         HUcontainsLQ=.true.,unitConversion=.true.)

    ! Write the zonal mean stddev to a text file, if requested
    if ( writeTextStddev ) then
      allocate(zonalMeanStddev(statevector%latPerPE * mpi_npey))
      do jvar = 1, vnl_numvarmax
        if (.not. gsv_varExist(varName=vnl_varNameList(jvar)) ) cycle

        write(*,*) ' writing zonal mean stddev to text file for variable: ', vnl_varNameList(jvar)
        field => gsv_getField3d_r8(statevector,vnl_varNameList(jvar))

        varName = vnl_varNameList(jvar)
        if ( varName == 'HU' ) varName = 'LQ'
        filename = 'stddev_ZM_' // trim(varName) // '_' // datestr // '.txt'
        nultxt = 0
        if ( mpi_myid == 0 ) ierr = fnom(nultxt,trim(filename),'FTN',0)

        do levIndex = 1, gsv_getNumLevFromVarName(statevector,vnl_varNameList(jvar))
          nsize = statevector%latPerPE
          call rpn_comm_gather(field(1,:,levIndex), nsize, 'mpi_double_precision',  &
                               zonalMeanStddev,     nsize, 'mpi_double_precision', 0, 'NS', ierr )
          if ( mpi_myid == 0 ) then
            do latIndex = 1, statevector%nj
              write(nultxt,*) field(1,latIndex,levIndex)
            end do
          end if
        end do

        if ( mpi_myid == 0 ) ierr = fclos(nulnam)

      end do
      deallocate(zonalMeanStddev)

    end if

    !
    !- Compute the domain mean std dev
    !
    write(*,*) 'Compute the domain mean stddev'
    call flush(6)

    allocate(stddev_dm(hco_anl%ni,nkgdim))
    allocate(stddev_dm2(hco_anl%ni,nkgdim))
    stddev_dm(:,:) = 0.d0
    stddev_dm2(:,:) = 0.d0

    !$OMP PARALLEL DO PRIVATE (latIndex,kIndex)
    do kIndex = 1, nkgdim
      do latIndex = statevector%myLatBeg, statevector%myLatEnd
        do lonIndex = statevector%myLonBeg, statevector%myLonEnd
          stddev_dm(lonIndex,kIndex) = stddev_dm(lonIndex,kIndex) + (field(lonIndex,latIndex,kIndex)**2)/real(hco_anl%nj,8)
        end do
      end do
    end do
    !$OMP END PARALLEL DO 

    nsize = statevector%ni*nkgdim
    call rpn_comm_allreduce(stddev_dm,stddev_dm2,nsize,  &
         "MPI_DOUBLE_PRECISION","MPI_SUM","GRID",ierr)

    !- Insert results in statevector
    !$OMP PARALLEL DO PRIVATE (lonIndex,latIndex,kIndex)
    do kIndex = 1, nkgdim
      do latIndex = statevector%myLatBeg, statevector%myLatEnd
        do lonIndex = statevector%myLonBeg, statevector%myLonEnd
          if(stddev_dm2(lonIndex,kIndex)> 0.0d0) then
            field(lonIndex,latIndex,kIndex) = sqrt(stddev_dm2(lonIndex,kIndex))
          else
            field(lonIndex,latIndex,kIndex) = 0.0d0
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    deallocate(stddev_dm)
    deallocate(stddev_dm2)
    deallocate(controlVector)
    deallocate(controlVector_global)

    call gsv_writeToFile(statevector,'stddev_' // datestr // '.fst','DM_STDDEV',  &
         HUcontainsLQ=.true.,unitConversion=.true.)

  end if ! if numperturbations.gt.1

  call gsv_deallocate(statevector)

  write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

  ! MPI, tmg finalize
  call tmg_stop(1)
  call tmg_terminate(mpi_myid, 'TMG_DIAGBMATRIX' )
  call rpn_comm_finalize(ierr) 

  write(*,*) ' --------------------------------'
  write(*,*) ' diagBmatrix ENDS'
  write(*,*) ' --------------------------------'

end program main_diagBmatrix
