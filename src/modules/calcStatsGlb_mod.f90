
module calcStatsGlb_mod
  ! MODULE calcStatsGlb_mod (prefix='csg' category='1. High-level functionality')
  !
  !:Purpose:  To compute homogeneous and isotropic background error covariances
  !           from forecast error estimate in model variable space (global
  !           version).
  !
  use codePrecision_mod
  use midasMpi_mod
  use gridStateVector_mod
  use gridStateVectorFileIO_mod
  use ensembleStateVector_mod
  use globalSpectralTransform_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use varNameList_mod
  use calcHeightAndPressure_mod
  use earthConstants_mod
  use utilities_mod
  use scaleDecomposition_mod
  use menetrierDiag_mod
  use fileNames_mod
  use timeCoord_mod
  use gridVariableTransforms_mod
  use gridBinning_mod
  use verticalModes_mod
  implicit none
  save
  private

  ! Public Subroutines
  public :: csg_setup, csg_computeBhi, csg_toolbox

  type(struct_hco), pointer :: hco_ens ! Ensemble horizontal grid parameters
  type(struct_vco), pointer :: vco_ens ! Ensemble horizontal grid parameters

  integer :: nens,ni,nj,nLevEns_M,nLevEns_T,nLevPtoT,nkgdimEns,varLevOffset(6),nla
  integer :: myLatBeg, myLatEnd, myLonBeg, myLonEnd, latPerPE, latPerPEmax, lonPerPE, lonPerPEmax
  integer :: mymBeg, mymEnd, mymSkip, mymCount, mynBeg, mynEnd, mynSkip, mynCount
  integer :: myMemberBeg, myMemberEnd, myMemberCount, maxMyMemberCount
  integer :: nEnsOverDimension
  integer :: nla_mpilocal, maxMyNla
  integer, pointer    :: ilaList_mpiglobal(:), ilaList_mpilocal(:)
  
  character(len=256), allocatable :: cflensin(:)
  integer :: gstID_nkgdimEns, gstID_nLevEns_M, gstID_nLevEns_T_P1
  integer, allocatable :: nip1_M(:),nip1_T(:)
  real(8), pointer :: pressureProfile_M(:), pressureProfile_T(:)
  integer,external    :: get_max_rss

  integer,parameter  :: nvar3d=4, nvar2d=1, nvar=nvar3d+nvar2d
  character(len=4) :: nomvar3d(nvar3d,3), nomvar2d(nvar2d,3), nomvar(nvar,3)
  integer, parameter :: modelSpace   = 1
  integer, parameter :: cvSpace      = 2
  integer, parameter :: cvUnbalSpace = 3

  real(8) :: gridSpacingInKm

  ! For wave band decomposition
  integer, parameter  :: maxNumLocalLength = 20
  integer             :: nHorizWaveBand
  integer             :: nVertWaveBand

  logical :: initialized = .false.

  ! Namelist variables
  integer :: ntrunc
  integer :: horizWaveBandPeaks(maxNumLocalLength) ! For horizontal wave band decomposition
  integer :: vertWaveBandPeaks(maxNumLocalLength) ! For vertical wave band decomposition
  
  contains

  !--------------------------------------------------------------------------
  ! CSG_SETUP
  !--------------------------------------------------------------------------
  subroutine csg_setup(nens_in, hco_in, vco_in)
    !
    !:Purpose: Main setup routine for this module
    !
    implicit none

    ! Arguments:
    integer,                   intent(in) :: nens_in
    type(struct_vco), pointer, intent(in) :: vco_in
    type(struct_hco), pointer, intent(in) :: hco_in

    ! Locals:
    integer :: nulnam, ierr, horizWaveBandIndex, vertWaveBandIndex, memberIndex
    integer :: fclos, fnom
    real(8) :: zps

    NAMELIST /NAMCALCSTATS_GLB/ntrunc,horizWaveBandPeaks,vertWaveBandPeaks

    write(*,*)
    write(*,*) 'csg_setup: Starting...'

    nens=nens_in
    allocate(cflensin(nens))
    do memberIndex = 1, nEns
      call fln_ensfileName(cflensin(memberIndex), 'ensemble', memberIndex_opt=memberIndex)
    end do
    
    if ( mmpi_myid == 0 ) then
      call mpc_printConstants(6)
      call pre_printPrecisions
    end if

    ! parameters from namelist (date in filename should come directly from sequencer?)
    ntrunc=108
    horizWaveBandPeaks(:) = -1.0d0
    vertWaveBandPeaks(:) = -1.0d0
    
    nulnam=0
    ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=NAMCALCSTATS_GLB)
    if (ierr /= 0) call utl_abort('csg_setup: Error reading namelist NAMCALCSTATS_GLB')
    if (mmpi_myid == 0) write(*,nml=NAMCALCSTATS_GLB)
    ierr=fclos(nulnam)

    !- Setup horizontal grid
    hco_ens => hco_in
    ni=hco_in%ni
    nj=hco_in%nj

    gridSpacingInKm = ec_ra * hco_in%dlon / 1000.d0
    write(*,*) 'Grid Spacing in Km = ', gridSpacingInKm

    ! setup mpi local grid parameters
    call mmpi_setup_latbands(nj, latPerPE, latPerPEmax, myLatBeg, myLatEnd)
    call mmpi_setup_lonbands(ni, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd)
    
    !- Setup vertical levels
    vco_ens => vco_in
    nLevEns_M=vco_in%nlev_M
    nLevEns_T=vco_in%nlev_T
    nLevPtot=nLevEns_M-1 ! ignore streamfunction at hyb=1, since highly correlated with next level
    varLevOffset(1) = 0
    varLevOffset(2) = 1*nLevEns_M
    varLevOffset(3) = 2*nLevEns_M
    varLevOffset(4) = 2*nLevEns_M+1*nLevEns_T
    varLevOffset(5) = 2*nLevEns_M+2*nLevEns_T
    nkgdimEns=nLevEns_M*2+nLevEns_T*2+1 ! NO TG !!!
    nla=(ntrunc+1)*(ntrunc+2)/2

    !- Setup the global spectral transform 
    gstID_nkgdimEns = gst_setup(ni,nj,ntrunc,nkgdimEns)
    gstID_nLevEns_M = gst_setup(ni,nj,ntrunc,nLevEns_M)
    gstID_nLevEns_T_P1 = gst_setup(ni,nj,ntrunc,nLevEns_T+1)

    ! setup mpi local spectral vector parameters
    call mmpi_setup_m(ntrunc, mymBeg, mymEnd, mymSkip, mymCount)
    call mmpi_setup_n(ntrunc, mynBeg, mynEnd, mynSkip, mynCount)
    call gst_ilaList_mpiglobal(ilaList_mpiglobal, nla_mpilocal, maxMyNla,  &
         gstID_nkgdimEns, mymBeg, mymEnd, mymSkip, mynBeg, mynEnd, mynSkip)
    call gst_ilaList_mpilocal(ilaList_mpilocal,  &
         gstID_nkgdimEns, mymBeg, mymEnd, mymSkip, mynBeg, mynEnd, mynSkip)

    ! setup ensemble members mpi partinionning (when working with struct_ens)
    call mmpi_setup_levels(nEns,myMemberBeg,myMemberEnd,myMemberCount)
    call rpn_comm_allreduce(myMemberCount, maxMyMemberCount, &
                            1,"MPI_INTEGER","MPI_MAX","GRID",ierr)
    nEnsOverDimension = mmpi_npex * maxMyMemberCount
    
    !- Setup ip1s
    allocate(nip1_M(nLevEns_M))
    nip1_M(:)=vco_in%ip1_M(:)
    allocate(nip1_T(nLevEns_T))
    nip1_T(:)=vco_in%ip1_T(:)

    !- Estimate the pressure profile for each vertical grid    
    zps = 101000.D0
    call czp_fetch1DLevels(vco_in, zps, &
                           profM_opt=pressureProfile_M, profT_opt=pressureProfile_T)

    !
    !- Setup variable names
    !
    nomvar3d(1,modelSpace)='UU'
    nomvar3d(2,modelSpace)='VV'
    nomvar3d(3,modelSpace)='TT'
    nomvar3d(4,modelSpace)='LQ'
    nomvar2d(1,modelSpace)='P0'
    
    nomvar3d(1,cvSpace)='PP'
    nomvar3d(2,cvSpace)='CC'
    nomvar3d(3,cvSpace)='TT'
    nomvar3d(4,cvSpace)='LQ'
    nomvar2d(1,cvSpace)='P0'
    
    nomvar3d(1,cvUnbalSpace)='PP'
    nomvar3d(2,cvUnbalSpace)='UC'
    nomvar3d(3,cvUnbalSpace)='UT'
    nomvar3d(4,cvUnbalSpace)='LQ'
    nomvar2d(1,cvUnbalSpace)='UP'

    nomvar(1:4,:) = nomvar3d(1:4,:)
    nomvar(5,:)   = nomvar2d(1,:)
    
    !
    !- Horizontal wave band decomposition option
    !
    nHorizWaveBand = count(horizWaveBandPeaks >= 0)
    if ( nHorizWaveBand < 1 ) then
       nHorizWaveBand = 1
    else if (nHorizWaveBand == 1) then
       write(*,*) 'You have specified only ONE horizWaveBandPeaks'
       call utl_abort('calbmatrix_glb')
    else
       write(*,*)
       write(*,*) 'Horizontal waveBand decomposition is ACTIVATED'
    end if
    
    ! Make sure that the wavenumbers are in the correct (decreasing) order
    do horizWaveBandIndex = 1, nHorizWaveBand-1
       if ( horizWaveBandPeaks(horizWaveBandIndex)-horizWaveBandPeaks(horizWaveBandIndex+1) <= 0 ) then
          write(*,*) 'csg_setup: horizWaveBandPeaks are not in decreasing wavenumber order'
          call utl_abort('calbmatrix_glb')
       end if
    end do

    ! Make sure the truncation is compatible with the horizWaveBandPeaks
    if ( ntrunc < nj-1 .and. nHorizWaveBand > 1 ) then
       write(*,*) 'csg_setup: The truncation is not compatible with wave band decomposition'
       write(*,*) '                 ntrunc should = ', nj-1
       call utl_abort('calbmatrix_glb')
    end if

    !
    !- Vertical wave band decomposition option
    !
    nVertWaveBand = count(vertWaveBandPeaks >= 0)
    if ( nVertWaveBand < 1 ) then
       nVertWaveBand = 1
    else if (nVertWaveBand == 1) then
       write(*,*) 'You have specified only ONE horizWaveBandPeaks'
       call utl_abort('calbmatrix_glb')
    else
       write(*,*)
       write(*,*) 'Vertical waveBand decomposition is ACTIVATED'
    end if
    
    ! Make sure that the wavenumbers are in the correct (decreasing) order
    do vertWaveBandIndex = 1, nVertWaveBand-1
       if ( vertWaveBandPeaks(vertWaveBandIndex)-vertWaveBandPeaks(vertWaveBandIndex+1) <= 0 ) then
          write(*,*) 'csg_setup: vertWaveBandPeaks are not in decreasing mode order'
          call utl_abort('calbmatrix_glb')
       end if
    end do
    
    !
    !- Ending
    !
    initialized = .true.

    write(*,*)
    write(*,*) 'csg_setup: Done!'

  end subroutine csg_setup

  !--------------------------------------------------------------------------
  ! csg_computeBhi
  !--------------------------------------------------------------------------
  subroutine csg_computeBhi()
    !
    !:Purpose: Master routine for Bhi computation in global mode
    !
    implicit none

    ! Locals:
    integer :: ierr, nulnam
    integer :: fclos, fnom
    character(len=12) :: formulation ! Bhi formulation

    NAMELIST /NAMCOMPUTEBHI/formulation

    ! parameters from namelist
    formulation='legacy'

    nulnam=0
    ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=NAMCOMPUTEBHI)
    if (ierr /= 0) call utl_abort('csg_computeBhi: Error reading namelist NAMCOMPUTEBHI')
    if (mmpi_myid == 0) write(*,nml=NAMCOMPUTEBHI)
    ierr=fclos(nulnam)

    select case(trim(formulation))
    case ('legacy')
      call csg_computeBhiLegacy
    case ('latbands')
      call csg_computeBhiLatBands
    case default
      write(*,*)
      write(*,*) 'csg_computeBhi: Unknown value of FORMULATION for Bhi computation: ',formulation
      write(*,*) 'Please select legacy or latbands'
      call utl_abort('csg_computeBhi')
    end select

  end subroutine csg_computeBhi
  
  !--------------------------------------------------------------------------
  ! csg_computeBhiLegacy
  !--------------------------------------------------------------------------
  subroutine csg_computeBhiLegacy
    !
    !:Purpose: Computation of Bhi using the legacy formulation
    !
    implicit none

    ! Locals:
    real(4), pointer :: ensPerturbations(:,:,:,:), ens_ptr(:,:,:,:)
    real(4), pointer :: ensBalPerturbations(:,:,:,:)
    real(8), pointer :: stddev3d(:,:,:), stddev3dBal(:,:,:), stddev3dUnbal(:,:,:)
    real(8), pointer :: stddev3d_ptr(:,:,:)
    real(8), pointer :: stddevZonAvg(:,:), stddevZonAvgBal(:,:), stddevZonAvgUnbal(:,:)
    real(8), allocatable :: PtoT(:,:,:),theta1(:,:),theta2(:,:)
    real(8), allocatable :: corns(:,:,:),rstddev(:,:)
    integer :: variableType 

    allocate(ensPerturbations(myLonBeg:myLonEnd, myLatBeg:myLatEnd, nkgdimEns, nens))
    allocate(ensBalPerturbations(myLonBeg:myLonEnd, myLatBeg:myLatEnd, nLevEns_T+1, nens))
    allocate(theta1(nlevEns_M,nj))
    allocate(stddev3d(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdimEns))
    allocate(stddevZonAvg(nkgdimEns,nj))
    allocate(PtoT(nlevEns_T+1,nlevEns_M,nj))
    allocate(theta2(nlevEns_M,nj))
    allocate(stddev3dBal(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLevEns_T+1))
    allocate(stddev3dUnbal(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdimEns))
    allocate(stddevZonAvgBal(nLevEns_T+1,nj))
    allocate(stddevZonAvgUnbal(nkgdimEns,nj))
    allocate(corns(nkgdimEns,nkgdimEns,0:ntrunc))
    allocate(rstddev(nkgdimEns,0:ntrunc))

    write(*,*) 'Initializing ensemble arrays to claim memory'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    ensPerturbations(:,:,:,:) = 0.0
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    ensBalPerturbations(:,:,:,:) = 0.0
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    call readEnsemble(ensPerturbations)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    call removeMean(ensPerturbations)

    call uv_to_psichi(ensPerturbations)

    call calcStddev3d(ensPerturbations,stddev3d,nkgdimens)

    call calcZonAvg(stddevZonAvg,stddev3d,nkgdimens)

    call calcTheta(ensPerturbations,theta1) ! theta1 is put in glbcov and used for analysis!
    if (mmpi_myid == 0) write(301,*) theta1

    call removeBalancedChi(ensPerturbations,theta1)

    call normalize3d(ensPerturbations,stddev3d)

    call calcPtoT(ensPerturbations,PtoT)
    if (mmpi_myid == 0) write(303,*) PTOT(:,:,1)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

!    call calcTheta(ensPerturbations,theta2) ! theta2 is used previously for computing unbalanced Chi!
!    if (mmpi_myid == 0) write(302,*) theta2

    call removeBalancedT_Ps(ensPerturbations,ensBalPerturbations,PtoT)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

!    call removeBalancedChi(ensPerturbations,theta2)

    call multiply3d(ensPerturbations,stddev3d,nkgdimens)

    ens_ptr(myLonBeg:,myLatBeg:,1:,1:) => ensBalPerturbations(:,:,1:nLevEns_T,:)
    stddev3d_ptr(myLonBeg:,myLatBeg:,1:) => stddev3d(:,:,(2*nLevEns_M+1):(2*nLevEns_M+nLevEns_T))
    call multiply3d(ens_ptr, stddev3d_ptr, nLevEns_T)

    ens_ptr(myLonBeg:,myLatBeg:,1:,1:) => ensBalPerturbations(:,:,(nLevEns_T+1):(nLevEns_T+1),:)
    stddev3d_ptr(myLonBeg:,myLatBeg:,1:) => stddev3d(:,:,(2*nLevEns_M+2*nLevEns_T+1):(2*nLevEns_M+2*nLevEns_T+1))
    call multiply3d(ens_ptr, stddev3d_ptr, 1)

    call spectralFilterLegacy(ensPerturbations,nkgdimens)

    call spectralFilterLegacy(ensBalPerturbations,nLevEns_T+1)

    call calcStddev3d(ensPerturbations,stddev3dUnbal,nkgdimens)

    call calcStddev3d(ensBalPerturbations,stddev3dBal,nLevEns_T+1)

    call calcZonAvg(stddevZonAvgUnbal,stddev3dUnbal,nkgdimens)

    call calcZonAvg(stddevZonAvgBal,stddev3dBal,nLevEns_T+1)

    call normalize3d(ensPerturbations,stddev3dUnbal)

    call removeGlobalMean(ensPerturbations)

    call calcCorrelations(ensPerturbations,corns,rstddev)

    variableType = cvUnbalSpace
    call writeStats(corns,rstddev,ptot,theta1)

    call writeStddev(stddevZonAvg,stddev3d,stddevZonAvgUnbal,stddev3dUnbal)

    call writeStddevBal(stddevZonAvgBal,stddev3dBal)

    call writeSpStats(ptot,theta1)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if (mmpi_myid == 0) then
      write(200,*) stddevZonAvg(1:nlevEns_M,:)
      write(201,*) stddevZonAvg((1+1*nlevEns_M):(2*nlevEns_M),:)
      write(202,*) stddevZonAvg((1+2*nlevEns_M):(3*nlevEns_T),:)
      write(203,*) stddevZonAvg((1+2*nlevEns_M+1*nlevEns_T):(2*nlevEns_M+2*nlevEns_T),:)
      write(204,*) stddevZonAvg((1+2*nlevEns_M+2*nlevEns_T),:)/1.0d2

      write(400,*) stddevZonAvgUnbal(1:nlevEns_M,:)
      write(401,*) stddevZonAvgUnbal((1+1*nlevEns_M):(2*nlevEns_M),:)
      write(402,*) stddevZonAvgUnbal((1+2*nlevEns_M):(3*nlevEns_T),:)
      write(403,*) stddevZonAvgUnbal((1+2*nlevEns_M+1*nlevEns_T):(2*nlevEns_M+2*nlevEns_T),:)
      write(404,*) stddevZonAvgUnbal((1+2*nlevEns_M+2*nlevEns_T),:)/1.0d2
    end if

  end subroutine csg_computeBhiLegacy

  !--------------------------------------------------------------------------
  ! csg_computeBhiLatBands
  !--------------------------------------------------------------------------
  subroutine csg_computeBhiLatBands
    !
    !:Purpose: Computation of Bhi on a set of latitude bands
    !
    implicit none

    ! Locals:
    integer :: variableType, latIndex, jlatband, lat1, lat2, lat3
    real(4), pointer     :: ensPerturbations(:,:,:,:)
    real(8), pointer     :: stddev3d(:,:,:)
    real(8), pointer     :: stddevZonAvg(:,:)
    real(8), allocatable :: corns(:,:,:),rstddev(:,:)
    real(8) :: latMask(nj)

    allocate(ensPerturbations(ni,nj,nkgdimEns,nens))
    allocate(stddev3d(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdimEns))
    allocate(stddevZonAvg(nkgdimEns,nj))
    allocate(corns(nkgdimEns,nkgdimEns,0:ntrunc))
    allocate(rstddev(nkgdimEns,0:ntrunc))

    call readEnsemble(ensPerturbations)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    call removeMean(ensPerturbations)

    call uv_to_psichi(ensPerturbations)

    call calcStddev3d(ensPerturbations,stddev3d,nkgdimens)

    call calcZonAvg(stddevZonAvg,stddev3d,nkgdimens)

    call normalize3d(ensPerturbations,stddev3d)

    call removeGlobalMean(ensPerturbations)

    do jlatband = 1, 3
      write(*,*) 'csg_computeBhiLatBands: selected LATBAND = ',jlatband
      lat1=nj/4
      lat2=nj/2
      lat3=3*nj/4
      write(*,*) 'lat1,2,3=',lat1,lat2,lat3
      if(jlatband==1) then
        ! Southern extratropics
        latMask(1:lat1) = 1.0d0
        do latIndex = lat1, lat2
          !latMask(latIndex) = sqrt(dble((lat2-latIndex)*4)/dble(nj))
          latMask(latIndex) = sqrt(0.5d0*(1.0d0+cos(dble((latIndex-lat1)*4)*MPC_PI_R8/dble(nj))))
        end do
        latMask(lat2:nj) = 0.0d0
      else if(jlatband==2) then
        ! Tropics
        !latMask(1:lat1) = 0.0d0
        !do latIndex = lat1, lat2
        !  latMask(latIndex) = sqrt(dble((latIndex-lat1)*4)/dble(nj))
        !end do
        !do latIndex = lat2,lat3
        !  latMask(latIndex) = sqrt(dble((lat3-latIndex)*4)/dble(nj))
        !end do
        !latMask(lat3:nj) = 0.0d0

        ! NOTE: use much broader band for tropics to avoid shortening horizontal correlations
        ! ok, since the masks do not have to sum to one for calculation, but they do
        ! when used in bmatrixhi_mod
        ! Tropics
        do latIndex = 1, lat1
          !latMask(latIndex) = sqrt(dble((latIndex-1)*4)/dble(nj))
          latMask(latIndex) = sqrt(0.5d0*(1.0d0+cos(dble((lat1-latIndex)*4)*MPC_PI_R8/dble(nj))))
        end do
        latMask(lat1:lat3) = 1.0d0
        do latIndex = lat3,nj
          !latMask(latIndex) = sqrt(dble((nj-latIndex)*4)/dble(nj))
          latMask(latIndex) = sqrt(0.5d0*(1.0d0+cos(dble((latIndex-lat3)*4)*MPC_PI_R8/dble(nj))))
        end do
      else if(jlatband==3) then
        ! Northern extratropics
        latMask(1:lat2) = 0.0d0
        do latIndex = lat2, lat3
          !latMask(latIndex) = sqrt(dble((latIndex-lat2)*4)/dble(nj))
          latMask(latIndex) = sqrt(0.5d0*(1.0d0+cos(dble((lat3-latIndex)*4)*MPC_PI_R8/dble(nj))))
        end do
        latMask(lat3:nj) = 1.0d0
      end if
      write(*,*) 'latMask = ',latMask(:)
      call calcCorrelations(ensPerturbations,corns,rstddev,latMask_opt=latMask)

      variableType = cvSpace
      call writeStats(corns,rstddev,latBand_opt=jlatBand)
    end do

    call writeStddev(stddevZonAvg,stddev3d)

    if (mmpi_myid == 0) then
      write(200,*) stddevZonAvg(1:nlevEns_M,:)
      write(201,*) stddevZonAvg((1+1*nlevEns_M):(2*nlevEns_M),:)
      write(202,*) stddevZonAvg((1+2*nlevEns_M):(3*nlevEns_T),:)
      write(203,*) stddevZonAvg((1+2*nlevEns_M+1*nlevEns_T):(2*nlevEns_M+2*nlevEns_T),:)
      write(204,*) stddevZonAvg((1+2*nlevEns_M+2*nlevEns_T),:)/1.0d2
    end if

  end subroutine csg_computeBhiLatBands

  !--------------------------------------------------------------------------
  ! CSG_TOOLBOOX
  !--------------------------------------------------------------------------
  subroutine csg_toolbox
    !
    !:Purpose: High-level routine to do a variety of diagnostic operations
    !          on the ensemble of error samples
    !
    implicit none

    ! NOTE: The diagnostic computed here are in model variable space 
    !       (no variable transform!!!)

    ! Locals:
    integer :: waveBandIndex, vertWaveBandIndex
    integer :: nulnam, ierr, fclos, fnom, numStep
    real(8), allocatable :: corns(:,:,:), rstddev(:,:), powerSpec(:,:)
    integer, allocatable :: dateStampList(:)    
    type(struct_ens), target  :: ensPerts
    type(struct_ens), target  :: ensPertsScaleDecomp(1)
    type(struct_gsv) :: statevector_template
    type(struct_gbi) :: gbi_zonalMean
    type(struct_gbi) :: gbi_globalMean    
    type(struct_vms) :: vModes    
    integer :: variableType 
    logical :: ensContainsFullField
    logical :: makeBiPeriodic
    logical :: doSpectralFilter
    real(8) :: vertModesLengthScale(2)
    character(len=60) :: tool
    character(len=2)  :: wbnum
    character(len=2)  :: ctrlVarHumidity
    
    NAMELIST /NAMTOOLBOX/tool, ensContainsFullField, doSpectralFilter, vertModesLengthScale, &
                         ctrlVarHumidity

    write(*,*)
    write(*,*) 'csg_toolbox'
    write(*,*)

    !
    !- Set options
    !
    variableType             = modelSpace ! hardwired
    ensContainsFullField     = .true.     ! default value
    doSpectralFilter         = .true.
    vertModesLengthScale(1) = 6.d0
    vertModesLengthScale(2) = -1.d0
    ctrlVarHumidity          = 'HU'
    
    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=NAMTOOLBOX)
    if (ierr /= 0) call utl_abort('csg_toolbox: Error reading namelist NAMTOOLBOX')
    if (mmpi_myid == 0) write(*,nml=NAMTOOLBOX)
    ierr = fclos(nulnam)

     if (vertModesLengthScale(2) == -1.d0) then
      vertModesLengthScale(2) = vertModesLengthScale(1)
     end if
    
    !
    !- Read ensemble
    !
    numStep = 1
    allocate(dateStampList(numStep))
    dateStampList(:)  = -1
    call ens_allocate(ensPerts, nEns, numStep, hco_ens, vco_ens, dateStampList)

    makeBiPeriodic = .false.
    call ens_readEnsemble(ensPerts, './ensemble', makeBiPeriodic, &
                          containsFullField_opt=ensContainsFullField)

    if ( ctrlVarHumidity == 'LQ' .and. ensContainsFullField ) then
      call gvt_transform(ensPerts,'HUtoLQ')
      call ens_modifyVarName(ensPerts, 'HU', 'LQ')
    end if

    !
    !- Compute and remove the ensemble mean; compute the stdDev
    !
    call ens_computeMean(ensPerts)
    call ens_removeMean (ensPerts)
    
    !
    !- Tool selection
    !
    select case(trim(tool))
    case ('HVCORREL_HI')
      write(*,*)
      write(*,*) 'Computing Homogeneous and Isotropic Correlation'

      if (mmpi_nprocs > 1) then
        call utl_abort('csg_toolbox: this tool is not yet MPI capable') ! only due to horizCorrelFunction
      end if
       
      call spectralFilter(ensPerts) ! INOUT
      call ens_computeStdDev(ensPerts)
      call ens_normalize(ensPerts)
      call ens_removeGlobalMean(ensPerts)

      allocate(corns(nkgdimEns,nkgdimEns,0:ntrunc))
      allocate(rstddev(nkgdimEns,0:ntrunc))
       
      call calcCorrelations2(ensPerts, & ! IN
                             corns,        & ! OUT (vertical correlation in spectral space)
                             rstddev)        ! OUT ( sqrt(normalized power spectrum) )

      call writeStats(corns,rstddev,waveBandIndex_opt=1) ! IN
      call calcHorizScale(rstddev,variableType,waveBandIndex_opt=1) ! IN
      call horizCorrelFunction(rstddev,variableType,waveBandIndex_opt=1) ! IN

      deallocate(rstddev)
      deallocate(corns)
       
    case ('HVCORREL_LOCAL')
      write(*,*)
      write(*,*) 'Computing Local Correlation'

      call ens_removeGlobalMean(ensPerts)
      call spectralFilter(ensPerts) ! INOUT
      call ens_computeStdDev(ensPerts)
      call ens_normalize(ensPerts)
      call calcLocalCorrelations(ensPerts) ! IN

    case ('VCORRMATRIX')
      write(*,*)
      write(*,*) 'Computing Local Correlation'

      if (doSpectralFilter) then
        call ens_removeGlobalMean(ensPerts)
        call spectralFilter(ensPerts) ! INOUT
      end if
      call ens_computeStdDev(ensPerts)
      call ens_normalize(ensPerts)
      call calcLocalVertCorrMatrix(ensPerts) ! IN

    case ('LOCALIZATIONRADII')
      write(*,*)
      write(*,*) 'Estimating the optimal covariance localization radii'

      if (nHorizWaveBand > 1 .and. nVertWaveBand > 1) then
        call utl_abort('csg_toolbox: cannot do horizontal AND vertical scale-decomposition')
      end if

      call ens_allocate(ensPertsScaleDecomp(1), nEns, numStep, hco_ens, vco_ens, dateStampList)
      if ( ctrlVarHumidity == 'LQ' .and. ensContainsFullField ) then
        call ens_modifyVarName(ensPertsScaleDecomp(1), 'HU', 'LQ')
      end if
      
      call ens_removeGlobalMean(ensPerts)
      call spectralFilter(ensPerts)

      do waveBandIndex = 1, max(nHorizWaveBand,nVertWaveBand)

        call ens_copy(ensPerts,                & ! IN
                      ensPertsScaleDecomp(1))    ! OUT
        
        if (nVertWaveBand > 1) then
          call scd_vertical(ensPertsScaleDecomp,                                    & ! INOUT
                            nVertWaveBand, vertWaveBandPeaks, vertModesLengthScale, & ! IN
                            'Select', vertWaveBandIndexSelected_opt=waveBandIndex,  & ! IN
                            writeResponseFunction_opt=.true.)                         ! IN
        else if (nHorizWaveBand > 1) then
          call scd_horizontal(ensPertsScaleDecomp,                                   & ! INOUT
                              nEnsOverDimension, nHorizWaveBand, horizWaveBandPeaks, & ! IN
                              'Select', 'SumToOne',                                  & ! IN
                              horizWaveBandIndexSelected_opt=waveBandIndex,          & ! IN
                              writeResponseFunction_opt=.true.)                        ! IN
        end if

        call ens_computeStdDev(ensPertsScaleDecomp(1))
       
        if (waveBandIndex == 1) then
          call ens_copyEnsStdDev(ensPertsScaleDecomp(1), statevector_template) ! IN
          call bmd_setup(statevector_template, hco_ens, nEns, pressureProfile_M, & ! IN
                         pressureProfile_T, max(nHorizWaveBand,nVertWaveBand))     ! IN
        end if

        call bmd_localizationRadii(ensPertsScaleDecomp(1), waveBandIndex_opt=waveBandIndex) ! IN

      end do

    case ('STDDEV')
      write(*,*)
      write(*,*) 'Computing standard deviations using PSI/CHI for winds'

      ! u,v to psi,chi 
      call gvt_setup(hco_ens, hco_ens, vco_ens)
      call gvt_transform(ensPerts, 'UVtoPsiChi')
      if (.not. ens_varExist(ensPerts,'PP') .and. &
          .not. ens_varExist(ensPerts,'CC') ) then
        call ens_modifyVarName(ensPerts, 'UU', 'PP')
        call ens_modifyVarName(ensPerts, 'VV', 'CC')
      end if

      ! Compute the grid point std dev
      call ens_computeStdDev(ensPerts)
      call ens_copyEnsStdDev(ensPerts, statevector_template) ! IN
      call gio_writeToFile(statevector_template, './stddev.fst', 'STDDEV_GRIDP', &
                           typvar_opt = 'E', numBits_opt = 32)

      ! Compute the zonal std dev
      call gbi_setup(gbi_zonalMean, 'YrowBand', statevector_template, hco_ens)
      call gbi_stdDev(gbi_zonalMean, ensPerts, & ! IN
                      statevector_template)          ! OUT
      call gio_writeToFile(statevector_template, './stddev.fst', 'STDDEV_ZONAL', &
                           typvar_opt = 'E', numBits_opt = 32)

    case ('POWERSPEC')
      write(*,*)
      write(*,*) 'Computing power spectra'

      call calcPowerSpec(ensPerts, & ! IN
                         powerSpec)      ! OUT
      call writePowerSpec(powerSpec, modelSpace) ! IN

    case ('VERTMODES_SPEC')
      write(*,*)
      write(*,*) 'Computing vertical modes spectra'

      if (nHorizWaveBand > 1 .or. nVertWaveBand > 1) then
        call utl_abort('csg_toolbox: waveband decomposition cannot be use when TOOLBOX=VERTMODES_SPEC')
      end if

      call vms_computeModesFromFunction(vco_ens, vertModesLengthScale(1), & ! IN
                                        vertModesLengthScale(2),          & ! IN
                                        vModes)                             ! OUT
      call vms_writeModes(vModes)

      call calcVertModesSpec(ensPerts,vModes) ! IN

    case ('VERTMODES_WAVEBAND')
      write(*,*)
      write(*,*) 'Computing vertical-scale decomposed ensemble perturbations'

      if (nVertWaveBand == 1) then
        call utl_abort('csg_toolbox: please specify a vertical waveband decomposition when TOOLBOX=VERTMODES_WAVEBAND')
      end if

      call ens_allocate(ensPertsScaleDecomp(1), nEns, numStep, hco_ens, vco_ens, dateStampList)
      if ( ctrlVarHumidity == 'LQ' .and. ensContainsFullField ) then
        call ens_modifyVarName(ensPertsScaleDecomp(1), 'HU', 'LQ')
      end if
      
      do vertWaveBandIndex = 1, nVertWaveBand
        
        call ens_copy(ensPerts,                & ! IN
                      ensPertsScaleDecomp(1))    ! OUT
        
        call scd_vertical(ensPertsScaleDecomp,                                           & ! INOUT
                          nVertWaveBand, vertWaveBandPeaks, vertModesLengthScale,        & ! IN
                          'Select', vertWaveBandIndexSelected_opt=vertWaveBandIndex,     & ! IN
                          writeResponseFunction_opt=.true., writeTransformInfo_opt=.true.) ! IN
        
        write(wbnum,'(I2.2)') vertWaveBandIndex
        
        ! Compute the grid point std dev
        call ens_computeStdDev(ensPertsScaleDecomp(1))
        call ens_copyEnsStdDev(ensPertsScaleDecomp(1), & ! IN
                               statevector_template)     ! OUT
        call gio_writeToFile(statevector_template, './stddev_'//trim(wbnum)//'.fst', 'STD_GRIDP_'//trim(wbnum), &
                             typvar_opt = 'E', numBits_opt = 32)

        ! Compute the global std dev
        if (vertWaveBandIndex == 1) then
          call gbi_setup(gbi_globalMean, 'HorizontalMean', statevector_template, hco_ens)
        end if
        call gbi_stdDev(gbi_globalMean, ensPertsScaleDecomp(1), & ! IN
                        statevector_template)           ! OUT
        call gio_writeToFile(statevector_template, './stddev_'//trim(wbnum)//'.fst', 'STD_GLB_'//trim(wbnum), &
                             typvar_opt = 'E', numBits_opt = 32)

        ! Write the scale-decomposed ensemble perturbations for member=1
        call ens_copyMember(ensPertsScaleDecomp(1), stateVector_template, 1)
        call gio_writeToFile(statevector_template, './ensPert_0001_'//trim(wbnum)//'.fst', 'ENSPERT_'//trim(wbnum), &
                             numBits_opt = 32)
        
      end do

      ! Write the full ensemble perturbations for member=1
      call ens_copyMember(ensPerts, stateVector_template, 1)
      call gio_writeToFile(statevector_template, './ensPert_0001.fst', 'ENSPERT', &
                           numBits_opt = 32)

      ! Compute the grid point std dev for the full perturbations
      call ens_computeStdDev(ensPerts)
      call ens_copyEnsStdDev(ensPerts,           & ! IN
                             statevector_template) ! OUT
      call gio_writeToFile(statevector_template, './stddev.fst', 'STD_GRIDP', &
                           typvar_opt = 'E', numBits_opt = 32)

      ! Compute the global std dev for the full perturbations
      call gbi_stdDev(gbi_globalMean, ensPerts, & ! IN
                      statevector_template)       ! OUT
      call gio_writeToFile(statevector_template, './stddev.fst', 'STD_GLB', &
                           typvar_opt = 'E', numBits_opt = 32)

    case default
      write(*,*)
      write(*,*) 'Unknown TOOL in csg_toolbox : ', trim(tool)
      call utl_abort('csg_toolbox')
    end select

    !
    !- Write the estimated pressure profiles
    !
    if (vco_ens%vgridPresent) then
      call writePressureProfiles
    end if

    !
    !- Ending
    !
    call ens_deallocate(ensPerts)

    write(*,*)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    write(*,*) 'csl_toolbox: Done!'

  end subroutine csg_toolbox

  !--------------------------------------------------------------------------
  ! WRITESPSTATS
  !--------------------------------------------------------------------------
  subroutine writeSpStats(ptot,theta)
    !
    !:Purpose: Write the spectral representation of PtoT and THETA to files
    !
    implicit none

    ! Arguments:
    real(8), intent(in) :: PtoT(:,:,:)
    real(8), intent(in) :: theta(:,:)

    ! Locals:
    integer jn,ierr,ipak,latIndex,levIndex1,levIndex2,nlev
    integer fstouv,fnom,fstfrm,fclos
    integer ip1,ip3,kni,knj,idatyp,idateo
    integer :: nulstats
    real(8) :: bufz(nLevEns_M),bufyz(nj,nLevEns_M),zsp(0:ntrunc,nLevEns_M)
    real(8) :: bufptot(nj,(nLevEns_T+1)*nLevEns_M),spptot(0:ntrunc,(nLevEns_T+1)*nLevEns_M)
    real(8) :: zspptot(nLevEns_T+1,nLevEns_M)

    if (mmpi_myid /= 0) return

    nulstats=0
    ierr =  fnom  (nulstats,'./stats_sp.fst','RND',0)
    ierr =  fstouv(nulstats,'RND')

    ipak = -32
    idatyp = 5
    ip1 = 0
    ip3 = nens
    idateo = 0

    ! write out SP_THETA

    do latIndex = 1, nj
      do levIndex1 = 1, nLevEns_M
        bufyz(latIndex,levIndex1) = theta(levIndex1,latIndex)
      end do
    end do

    call gst_zlegdir(gstID_nkgdimEns,bufyz,zsp,nLevEns_M)

    do jn = 0, ntrunc
      do levIndex1=1, nLevEns_M
        bufz(levIndex1) = zsp(jn,levIndex1)
      end do

      ierr = utl_fstecr(bufz,ipak,nulstats,idateo,0,0,nlevEns_M,1,1,   &
                        ip1,jn,ip3,'X','ZZ','SP_THETA','X',0,0,0,0,idatyp,.true.)

    end do

    ! write out SP_PTOT

    do latIndex = 1, nj
      do levIndex1 = 1, (nLevEns_T+1)
        do levIndex2 = 1, nLevEns_M
          bufptot(latIndex,(levIndex2-1)*(nLevEns_T+1)+levIndex1) = PtoT(levIndex1,levIndex2,latIndex)
        end do
      end do
    end do

    nlev=(nLevEns_T+1)*nLevEns_M
    call gst_zlegdir(gstID_nkgdimEns,bufptot,spptot,nLev)

    do jn = 0, ntrunc
      do levIndex1 = 1, (nLevEns_T+1)
        do levIndex2 = 1, nLevEns_M
          zspptot(levIndex1,levIndex2) = spptot(jn,(levIndex2-1)*(nLevEns_T+1)+levIndex1)
        end do
      end do

      kni=nLevEns_T+1
      knj=nLevEns_M
      ierr = utl_fstecr(zspptot,ipak,nulstats,idateo,0,0,kni,knj,1,  &
                        ip1,jn,ip3,'X','ZZ','SP_PTOT ','X',0,0,0,0,idatyp,.true.)
    end do

    ierr =  fstfrm(nulstats)
    ierr =  fclos (nulstats)

    write(*,*) 'finished writing statistics...'

  end subroutine writeSpStats

  !--------------------------------------------------------------------------
  ! REMOVEBALANCEDCHI
  !--------------------------------------------------------------------------
  subroutine removeBalancedChi(ensPerturbations,theta)
    !
    !:Purpose: Subtract the balanced components of velocity potential
    !          from the full variable
    !
    implicit none

    ! Arguments:
    real(4), pointer, intent(inout) :: ensPerturbations(:,:,:,:)
    real(8),          intent(in)    :: theta(:,:)

    ! Locals:
    real(4), pointer :: psi_ptr(:,:,:), chi_ptr(:,:,:)
    integer :: ensIndex,latIndex,levIndex,lonIndex

    do ensIndex = 1,nens
      psi_ptr(myLonBeg:,myLatBeg:,1:) => ensPerturbations(:,:,1:nlevEns_M,ensIndex)
      chi_ptr(myLonBeg:,myLatBeg:,1:) => ensPerturbations(:,:,(nlevEns_M+1):(2*nlevEns_M),ensIndex)

      do levIndex = 1, nLevEns_M
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            chi_ptr(lonIndex,latIndex,levIndex) = chi_ptr(lonIndex,latIndex,levIndex) +  &
                 tan(theta(levIndex,latIndex))*psi_ptr(lonIndex,latIndex,levIndex)
          end do
        end do
      end do

    end do

    write(*,*) 'finished removing balanced chi...'

  end subroutine removeBalancedChi

  !--------------------------------------------------------------------------
  ! REMOVEBALANCEDT_PS
  !--------------------------------------------------------------------------
  subroutine removeBalancedT_Ps(ensPerturbations,ensBalPerturbations,PtoT)
    !
    !:Purpose: Subtract the balanced components of temperature and surface
    !          pressure from the full variables
    !
    implicit none

    ! Arguments:
    real(4), pointer, intent(inout) :: ensPerturbations(:,:,:,:)
    real(4), pointer, intent(inout) :: ensBalPerturbations(:,:,:,:)
    real(8),          intent(in)    :: PtoT(:,:,:)

    ! Locals:
    real(4),pointer :: tt_ptr(:,:,:), ps_ptr(:,:,:), ttb_ptr(:,:,:), psb_ptr(:,:,:)
    real(8) :: spectralState(nla_mpilocal,2,nLevEns_M), spBalancedP(nla_mpilocal,2,nlevEns_M)
    real(8) :: balancedP(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nlevEns_M)
    real(8) :: psi(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLevEns_M)
    integer :: ensIndex, latIndex, lonIndex, jk1, jk2

    do ensIndex=1,nens

      write(*,*) 'removing balanced T and Ps for member ', ensIndex

      psi(:,:,:) = ensPerturbations(:,:,1:nlevEns_M,ensIndex)
      call gst_setID(gstID_nLevEns_M)
      call gst_reespe(spectralState,psi)
      call calcBalancedP(spectralState,spBalancedP)
      call gst_speree(spBalancedP,balancedP)

      tt_ptr(myLonBeg:,myLatBeg:,1:)  => ensPerturbations(:,:,(1+2*nLevEns_M):(2*nLevEns_M+1*nLevEns_T),ensIndex)
      ps_ptr(myLonBeg:,myLatBeg:,1:)  => ensPerturbations(:,:,(1+2*nLevEns_M+2*nLevEns_T):(1+2*nLevEns_M+2*nLevEns_T),ensIndex)
      ttb_ptr(myLonBeg:,myLatBeg:,1:) => ensBalPerturbations(:,:,1:nLevEns_T,ensIndex)
      psb_ptr(myLonBeg:,myLatBeg:,1:) => ensBalPerturbations(:,:,(1+nLevEns_T):(1+nLevEns_T),ensIndex)

      ttb_ptr(:,:,:)=0.0d0
      psb_ptr(:,:,:)=0.0d0

      !$OMP PARALLEL DO PRIVATE (latIndex,lonIndex,jk1,jk2)
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          do jk1 = 1, nLevEns_T
            do jk2 = 1, nlevptot
              ttb_ptr(lonIndex,latIndex,jk1) = ttb_ptr(lonIndex,latIndex,jk1) + PtoT(jk1,jk2,latIndex)*balancedP(lonIndex,latIndex,jk2)
            end do
            tt_ptr(lonIndex,latIndex,jk1) = tt_ptr(lonIndex,latIndex,jk1) - ttb_ptr(lonIndex,latIndex,jk1)
          end do
          do jk2 = 1, nlevptot
            psb_ptr(lonIndex,latIndex,1) = psb_ptr(lonIndex,latIndex,1) + PtoT(nLevEns_T+1,jk2,latIndex)*balancedP(lonIndex,latIndex,jk2)
          end do
          ps_ptr(lonIndex,latIndex,1) = ps_ptr(lonIndex,latIndex,1) - psb_ptr(lonIndex,latIndex,1)
        end do
      end do
      !$OMP END PARALLEL DO

    end do

    write(*,*) 'finished removing balanced T and Ps...'

  end subroutine removeBalancedT_Ps

  !--------------------------------------------------------------------------
  ! CALCCORRELATIONS
  !--------------------------------------------------------------------------
  subroutine calcCorrelations(ensPerturbations,corns,rstddev,latMask_opt)
    !
    !:Purpose: Calculate the homogeneous and isotropic correlations in spectral space
    !
    implicit none

    ! Arguments:
    real(4), pointer,  intent(in)  :: ensPerturbations(:,:,:,:)
    real(8),           intent(out) :: corns(nkgdimEns,nkgdimEns,0:ntrunc)
    real(8),           intent(out) :: rstddev(nkgdimEns,0:ntrunc)
    real(8), optional, intent(in)  :: latMask_opt(:)

    ! Locals:
    real(8) :: spectralState(nla_mpilocal,2,nkgdimEns)
    real(8) :: corns_mpiglobal(nkgdimEns,nkgdimEns,0:ntrunc)
    real(8) :: gridState(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdimEns)
    real(8) :: dfact,dfact2,dsummed
    integer :: ensIndex,ila_mpilocal,ila_mpiglobal,jn,jm,jk1,jk2,latIndex,nsize,ierr

    call utl_tmg_start(120,'--CalcStats_Corr')

    corns(:,:,:) = 0.0d0
    do ensIndex = 1, nens

      write(*,*) 'calcCorrelations: processing member ',ensIndex

      gridState(:,:,:) = ensPerturbations(:,:,:,ensIndex)
      if(present(latMask_opt)) then
        do latIndex = myLatBeg, myLatEnd
          gridState(:,latIndex,:) = latMask_opt(latIndex)*gridState(:,latIndex,:)
        end do
      end if
      call gst_setID(gstID_nkgdimEns)
      call gst_reespe(spectralState,gridState)

      !$OMP PARALLEL DO PRIVATE (jn,jm,dfact,ila_mpilocal,ila_mpiglobal,jk1,jk2)
      do jn = mynBeg, mynEnd, mynSkip
        do jm = mymBeg, mymEnd, mymSkip
          if(jm.le.jn) then
            dfact = 2.0d0
            if (jm.eq.0) dfact = 1.0d0
            ila_mpiglobal = gst_getNind(jm,gstID_nkgdimEns) + jn - jm
            ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
            do jk1 = 1, nkgdimEns
              do jk2 = 1, nkgdimEns
                corns(jk1,jk2,jn) = corns(jk1,jk2,jn) +     &
                     dfact*( spectralState(ila_mpilocal,1,jk1)*spectralState(ila_mpilocal,1,jk2) +   &
                     spectralState(ila_mpilocal,2,jk1)*spectralState(ila_mpilocal,2,jk2)  )
              end do
            end do
          end if
        end do
      end do
      !$OMP END PARALLEL DO

    end do

    ! communicate between all tasks
    nsize = nkgdimEns*nkgdimEns*(1+ntrunc)
    call rpn_comm_allreduce(corns,corns_mpiglobal,nsize,"mpi_double_precision","mpi_sum","GRID",ierr)
    corns(:,:,:) = corns_mpiglobal(:,:,:)
    
    !$OMP PARALLEL DO PRIVATE (jn,jk1)
    do jn = 0, ntrunc
      do jk1 = 1, nkgdimEns
        if(abs(corns(jk1,jk1,jn)).gt.0.0d0) then
          rstddev(jk1,jn) = dsqrt(abs(corns(jk1,jk1,jn)))
        else
          rstddev(jk1,jn) = 0.0d0
        end if
      end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE (jn,jk1,jk2)
    do jn = 0, ntrunc
      do jk1 = 1, nkgdimEns
        do jk2 = 1, nkgdimEns
          if(rstddev(jk1,jn).ne.0..and.rstddev(jk2,jn).ne.0.) then
            corns(jk1,jk2,jn) =  corns(jk1,jk2,jn)/(rstddev(jk1,jn)*rstddev(jk2,jn))
          else
            corns(jk1,jk2,jn) = 0.0d0
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    dfact2 = 1.0d0/sqrt(dble(nens-1))
    do jn = 0, ntrunc
      dfact = 1.0d0/sqrt(2.0d0*dble(jn) + 1.0d0)
      do jk1 = 1, nkgdimEns
        rstddev(jk1,jn) = rstddev(jk1,jn)*dfact2*dfact
      end do
    end do

    ! Normalize to ensure correlations in horizontal and Multiply by sqrt(0.5) to make valid for m.ne.0
    !$OMP PARALLEL DO PRIVATE (jk1,jn,dsummed)
    do jk1 = 1, nkgdimEns
      dsummed=0.0d0
      do jn = 0, ntrunc
        dsummed=dsummed + (rstddev(jk1,jn)**2)*((2.0d0*dble(jn))+1.0d0)/2.0d0
      end do
      do jn = 0, ntrunc
        if(dsummed.gt.0.0d0) rstddev(jk1,jn)=rstddev(jk1,jn)*sqrt(0.5d0/dsummed)
      end do
    end do
    !$OMP END PARALLEL DO

    call utl_tmg_stop(120)
    write(*,*) 'finished computing correlations...'

  end subroutine calcCorrelations

  !--------------------------------------------------------------------------
  ! CALCCORRELATIONS2
  !--------------------------------------------------------------------------
  subroutine calcCorrelations2(ensPerts,corns,rstddev,latMask_opt)
    !
    !:Purpose: Calculate the homogeneous and isotropic correlations in spectral space
    !
    implicit none

    ! Arguments:
    type(struct_ens),  intent(inout) :: ensPerts
    real(8),           intent(out)   :: corns(nkgdimEns,nkgdimEns,0:ntrunc)
    real(8),           intent(out)   :: rstddev(nkgdimEns,0:ntrunc)
    real(8), optional, intent(in)    :: latMask_opt(:)

    ! Locals:
    real(4), pointer :: ptr4d_r4(:,:,:,:)
    real(8) :: corns_mpiglobal(nkgdimEns,nkgdimEns,0:ntrunc)
    real(8) :: spectralState(nla_mpilocal,2,nkgdimEns)
    real(8) :: gridState(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdimEns)
    real(8) :: dfact, dfact2, dsummed
    integer :: ensIndex, ila_mpilocal, ila_mpiglobal, jn, jm, jk1, jk2
    integer :: levIndex, latIndex, nsize, ierr

    call utl_tmg_start(120,'--CalcStats_Corr')

    corns(:,:,:) = 0.0d0
    do ensIndex = 1, nens

      write(*,*) 'calcCorrelations: processing member ',ensIndex

      !- 2.1 Extract fields from ensPerturbations
      do levIndex = 1, ens_getNumK(ensPerts)
        ptr4d_r4 => ens_getOneLev_r4(ensPerts,levIndex)
        gridState(:,:,levIndex) = real(ptr4d_r4(ensIndex,1,:,:),8)
      end do
      if ( present(latMask_opt) ) then
        do latIndex = myLatBeg, myLatEnd
          gridState(:,latIndex,:) = latMask_opt(latIndex)*gridState(:,latIndex,:)
        end do
      end if
      call gst_setID(gstID_nkgdimEns)
      call gst_reespe(spectralState,gridState)

      !$OMP PARALLEL DO PRIVATE (jn,jm,dfact,ila_mpilocal,ila_mpiglobal,jk1,jk2)
      do jn = mynBeg, mynEnd, mynSkip
        do jm = mymBeg, mymEnd, mymSkip
          if(jm.le.jn) then
            dfact = 2.0d0
            if (jm.eq.0) dfact = 1.0d0
            ila_mpiglobal = gst_getNind(jm,gstID_nkgdimEns) + jn - jm
            ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
            do jk1 = 1, nkgdimEns
              do jk2 = 1, nkgdimEns
                corns(jk1,jk2,jn) = corns(jk1,jk2,jn) +     &
                     dfact*( spectralState(ila_mpilocal,1,jk1)*spectralState(ila_mpilocal,1,jk2) +   &
                     spectralState(ila_mpilocal,2,jk1)*spectralState(ila_mpilocal,2,jk2)  )
              end do
            end do
          end if
        end do
      end do
      !$OMP END PARALLEL DO

    end do

    ! communicate between all tasks
    nsize = nkgdimEns*nkgdimEns*(1+ntrunc)
    call rpn_comm_allreduce(corns,corns_mpiglobal,nsize,"mpi_double_precision","mpi_sum","GRID",ierr)
    corns(:,:,:) = corns_mpiglobal(:,:,:)
    
    !$OMP PARALLEL DO PRIVATE (jn,jk1)
    do jn = 0, ntrunc
      do jk1 = 1, nkgdimEns
        if(abs(corns(jk1,jk1,jn)).gt.0.0d0) then
          rstddev(jk1,jn) = dsqrt(abs(corns(jk1,jk1,jn)))
        else
          rstddev(jk1,jn) = 0.0d0
        end if
      end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE (jn,jk1,jk2)
    do jn = 0, ntrunc
      do jk1 = 1, nkgdimEns
        do jk2 = 1, nkgdimEns
          if(rstddev(jk1,jn).ne.0..and.rstddev(jk2,jn).ne.0.) then
            corns(jk1,jk2,jn) =  corns(jk1,jk2,jn)/(rstddev(jk1,jn)*rstddev(jk2,jn))
          else
            corns(jk1,jk2,jn) = 0.0d0
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    dfact2 = 1.0d0/sqrt(dble(nens-1))
    do jn = 0, ntrunc
      dfact = 1.0d0/sqrt(2.0d0*dble(jn) + 1.0d0)
      do jk1 = 1, nkgdimEns
        rstddev(jk1,jn) = rstddev(jk1,jn)*dfact2*dfact
      end do
    end do

    ! Normalize to ensure correlations in horizontal and Multiply by sqrt(0.5) to make valid for m.ne.0
    !$OMP PARALLEL DO PRIVATE (jk1,jn,dsummed)
    do jk1 = 1, nkgdimEns
      dsummed=0.0d0
      do jn = 0, ntrunc
        dsummed=dsummed + (rstddev(jk1,jn)**2)*((2.0d0*dble(jn))+1.0d0)/2.0d0
      end do
      do jn = 0, ntrunc
        if(dsummed.gt.0.0d0) rstddev(jk1,jn)=rstddev(jk1,jn)*sqrt(0.5d0/dsummed)
      end do
    end do
    !$OMP END PARALLEL DO

    call utl_tmg_stop(120)
    write(*,*) 'finished computing correlations...'

  end subroutine calcCorrelations2

  !--------------------------------------------------------------------------
  ! CALCPOWERSPEC
  !--------------------------------------------------------------------------
  subroutine calcPowerSpec(ensPerts,powerSpec)
    !
    !:Purpose: Calculate the horizontal power spectrum of the ensemble
    !          of error samples
    !
    implicit none

    ! Arguments:
    type(struct_ens),     intent(inout) :: ensPerts
    real(8), allocatable, intent(out)   :: powerSpec(:,:)

    ! Locals:
    real(8), allocatable :: ensPertSP(:,:,:)
    real(8), allocatable :: ensPertGD(:,:,:)
    real(4), pointer     :: ptr4d_r4(:,:,:,:)
    real(8) :: dfact, dfact2
    integer :: gstPowerSpecID
    integer :: memberIndex, levIndex, latIndex, lonIndex
    integer :: jn, jm, ila_mpilocal, ila_mpiglobal

    allocate(powerSpec          (ens_getNumK(ensPerts),0:ntrunc))

    !
    !- Spectral decomposition and spectral coefficient summation
    !
    gstPowerSpecID = gst_setup(ni,nj,nTrunc,nEnsOverDimension)
    
    allocate(ensPertSP(nla_mpilocal,2,nEnsOverDimension))
    allocate(ensPertGD(nEnsOverDimension,myLonBeg:myLonEnd,myLatBeg:myLatEnd))

    powerSpec(:,:)=0.0d0
    
    do levIndex = 1, ens_getNumK(ensPerts) ! Loop on variables and vertical levels
      
      ptr4d_r4 => ens_getOneLev_r4(ensPerts,levIndex)

      do latIndex = myLatBeg, myLatEnd
        ensPertGD(:,:,latIndex) = 0.0d0
      end do

      !$OMP PARALLEL DO PRIVATE (memberIndex,latIndex,lonIndex)
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          do memberIndex = 1, nEns
            ensPertGD(memberIndex,lonIndex,latIndex) = dble(ptr4d_r4(memberIndex,1,lonIndex,latIndex))
          end do
        end do
      end do
      !$OMP END PARALLEL DO
      
      !- GridPoint space -> Spectral Space
      call gst_setID(gstPowerSpecID) ! IN
      call gst_reespe_kij(ensPertSP, & ! OUT
                          ensPertGD)   ! IN

      !$OMP PARALLEL DO PRIVATE (memberIndex,jn,jm,dfact,ila_mpilocal,ila_mpiglobal)
      do jn = mynBeg, mynEnd, mynSkip
        do jm = mymBeg, mymEnd, mymSkip
          if (jm .le. jn) then
            dfact = 2.0d0
            if (jm.eq.0) dfact = 1.0d0
            ila_mpiglobal = gst_getNind(jm,gstPowerSpecID) + jn - jm
            ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
            do memberIndex = 1, nEns
              powerSpec(levIndex,jn) = powerSpec(levIndex,jn) +     &
                   dfact*( ensPertSP(ila_mpilocal,1,memberIndex)**2 + ensPertSP(ila_mpilocal,2,memberIndex)**2 )
            end do
          end if
        end do
      end do
      !$OMP END PARALLEL DO

    end do

    !
    !- Communicate between all tasks
    !
    call mmpi_allreduce_sumR8_2d(powerSpec, "GRID")

    !
    !- Apply the appropriate scaling
    !
    dfact2 = 1.0d0/sqrt(dble(nens-1))
    do jn = 0, ntrunc
      dfact = 1.0d0/sqrt(2.0d0*dble(jn) + 1.0d0)
      do levIndex = 1, ens_getNumK(ensPerts)
        powerSpec(levIndex,jn) = powerSpec(levIndex,jn)*dfact2*dfact
      end do
    end do

    write(*,*) 'finished computing power spectrum...'

  end subroutine calcPowerSpec
  
  !--------------------------------------------------------------------------
  ! WRITESTATS
  !--------------------------------------------------------------------------
  subroutine writeStats(corns, rstddev, ptot_opt, theta_opt, waveBandIndex_opt, latBand_opt)
    !
    !:Purpose: Write several components of the BHI matrix to a file
    !
    implicit none

    ! Arguments:
    real(8),           intent(in) :: corns(nkgdimEns,nkgdimEns,0:ntrunc)
    real(8),           intent(in) :: rstddev(nkgdimEns,0:ntrunc)
    real(8), optional, intent(in) :: PtoT_opt(:,:,:)
    real(8), optional, intent(in) :: theta_opt(:,:)
    integer, optional, intent(in) :: waveBandIndex_opt
    integer, optional, intent(in) :: latBand_opt

    ! Locals:
    real(8) :: prcor(nkgdimEns,nkgdimEns)
    integer :: jn,ierr,ipak,jk,jl
    integer :: fstouv,fnom,fstfrm,fclos
    integer :: ip1,ip2,ip3,idatyp,idateo
    integer :: nulstats
    character(len=128) :: outfilename
    character(len=2) :: wbnum

    if (mmpi_myid /= 0) return
    
    nulstats=0
    if ( nHorizWaveBand == 1 ) then
       outfilename='./bgcov.fst'
    else
       if (.not. present(waveBandIndex_opt)) then
          write(*,*) 'writeStats: No waveBandIndex was supplied!!!'
          call utl_abort('calbmatrix_glb')
       end if
       write(wbnum,'(I2.2)') waveBandIndex_opt
       outfilename='./bgcov_'//trim(wbnum)//'.fst'
    end if
    ierr =  fnom  (nulstats,trim(outfilename),'RND',0)
    ierr =  fstouv(nulstats,'RND')

    ipak = -32
    idatyp = 5
    ip1 = 0
    ip2 = 0
    ip3 = nens
    idateo = 0

    if(present(latBand_opt)) ip1 = latBand_opt

    if (present(ptot_opt)) then
       ierr = utl_fstecr(ptot_opt,ipak,nulstats,idateo,0,0,nlevEns_T+1,nlevEns_M,nj,  &
                         ip1,ip2,ip3,'X','ZZ','P_to_T  ','X',0,0,0,0,idatyp,.true.)
    end if
    if (present(theta_opt)) then
       ierr = utl_fstecr(theta_opt,ipak,nulstats,idateo,0,0,nlevEns_M,nj,1,   &
                         ip1,ip2,ip3,'X','ZZ','THETA   ','X',0,0,0,0,idatyp,.true.)
    end if

    do jn = 0, ntrunc
      ip2 = jn
      ierr = utl_fstecr(corns(:,:,jn),ipak,nulstats,idateo,0,0,nkgdimEns,nkgdimEns,1,  &
                        ip1,ip2,ip3,'X','ZZ','CORRNS  ','X',0,0,0,0,idatyp,.true.)
    end do

    do jn = 0, ntrunc
      ip2 = jn
      ierr = utl_fstecr(rstddev(:,jn),ipak,nulstats,idateo,0,0,nkgdimEns,1,1,   &
                        ip1,ip2,ip3,'X','SS','RSTDDEV ','X',0,0,0,0,idatyp,.true.)
    end do

    
    ! Computing the total vertical correlation matrix
    do jk = 1, nkgdimEns
      do jl = 1, nkgdimEns
        prcor(jk,jl) = 0.0d0
        do jn = 0, ntrunc
          prcor(jk,jl) = prcor(jk,jl) + ((2*jn+1)*rstddev(jk,jn)*rstddev(jl,jn)*corns(jk,jl,jn))
        end do
      end do
    end do

    do jk = 1, nkgdimEns
      do jl = 1, nkgdimEns
        if(prcor(jk,jk)*prcor(jl,jl) .gt. 0.0d0) then
          prcor(jk,jl) = prcor(jk,jl) / (sqrt(prcor(jk,jk)*prcor(jl,jl)))
        else
          prcor(jk,jl) = 0.0d0
        end if
      end do
    end do
    ip2 =0
    ierr = utl_fstecr(prcor(:,:),ipak,nulstats,idateo,0,0,nkgdimEns,nkgdimEns,1,   &
                      ip1,ip2,ip3,'X','ZV','CORVERT ','X',0,0,0,0,idatyp,.true.)

    ierr =  fstfrm(nulstats)
    ierr =  fclos (nulstats)

    write(*,*) 'finished writing statistics...'

  end subroutine writeStats

  !--------------------------------------------------------------------------
  ! WRITEPRESSUREPROFILES
  !--------------------------------------------------------------------------
  subroutine writePressureProfiles
    !
    !:Purpose: Write the profiles of pressure to ascii files
    !
    implicit none

    ! Locals:
    character(len=128) :: outfilename
    integer :: jk

    if (mmpi_myid /= 0) return

    outfilename = "./pressureProfile_M.txt"
    open (unit=99,file=outfilename,action="write",status="new")
    do jk = 1, nLevEns_M
       write(99,'(I3,2X,F7.2)') jk, pressureProfile_M(jk)/100.d0
    end do
    close(unit=99)
       
    outfilename = "./pressureProfile_T.txt"
    open (unit=99,file=outfilename,action="write",status="new")
    do jk = 1, nLevEns_T
       write(99,'(I3,2X,F7.2)') jk, pressureProfile_T(jk)/100.d0
    end do
    close(unit=99)

    write(*,*) 'finished writing pressure profiles...'

  end subroutine writePressureProfiles

  !--------------------------------------------------------------------------
  ! WRITESTDDEV
  !--------------------------------------------------------------------------
  subroutine writeStddev(stddevZonAvg,stddev3d,stddevZonAvgUnbal_opt,stddev3dUnbal_opt)
    !
    !:Purpose: Write the stddev to a file
    !
    implicit none

    ! Arguments:
    real(8), pointer,           intent(in) :: stddevZonAvg(:,:)
    real(8), pointer,           intent(in) :: stddev3d(:,:,:)
    real(8), pointer, optional, intent(in) :: stddevZonAvgUnbal_opt(:,:)
    real(8), pointer, optional, intent(in) :: stddev3dUnbal_opt(:,:,:)

    ! Locals:
    type(struct_gsv) :: stateVector
    real(8) :: dfact, zbufyz(nj,max(nLevEns_M,nLevens_T)), zbufy(nj)
    integer :: latIndex, levIndex, ierr, varIndex, varIndexStddev, nLevEns, numVarToWrite
    integer :: ip1, ip2, ip3, idatyp, idateo, numBits
    integer :: nulstats
    real(8), pointer :: field(:,:,:)
    integer, allocatable :: dateStampList(:)
    character(len=4) :: nomVarToWrite(1:20)
    integer :: fstouv, fnom, fstfrm, fclos
    
    numBits = 32
    idatyp = 5
    ip1 = 0
    ip2 = 0
    ip3 = nens
    idateo = 0

    ! figure out full list of variables to write
    numVarToWrite = size(nomvar,1)
    nomVarToWrite(1:numVarToWrite) = nomvar(:,cvSpace)
    if (present(stddevZonAvgUnbal_opt) .and. present(stddev3dUnbal_opt)) then
      do varIndex = 1, nvar
        if( all(nomvar(:,cvSpace) /= nomvar(varIndex,cvUnbalSpace)) ) then
          numVarToWrite = numVarToWrite + 1
          nomVarToWrite(numVarToWrite) = nomvar(varIndex,cvUnbalSpace)
        end if
      end do
    end if

    ! first, write stddev3d
    
    ! allocate stateVector used for writing stddev3d
    allocate(dateStampList(1))
    dateStampList(:)  = 0
    call gsv_allocate( stateVector, 1, hco_ens, vco_ens,  &
                       mpi_local_opt=.true., dateStampList_opt=dateStampList,   &
                       varNames_opt=nomVarToWrite(1:numVarToWrite) )
    do varIndex = 1, numVarToWrite
      nLevEns = gsv_getNumLevFromVarName(stateVector,nomVarToWrite(varIndex))
      call gsv_getField(stateVector,field,nomVarToWrite(varIndex))
      if ( any(nomVarToWrite(varIndex) == nomvar(:,cvSpace)) ) then
        field(:,:,:) = stddev3d(:, :,      (varLevOffset(varIndex)+1):(varLevOffset(varIndex)+nLevEns) )
      else if(present(stddevZonAvgUnbal_opt) .and. present(stddev3dUnbal_opt)) then
        varIndexStddev = findloc(nomvar(:,cvUnbalSpace), value=nomVarToWrite(varIndex), dim=1)
        field(:,:,:) = stddev3dUnbal_opt(:, :, (varLevOffset(varIndexStddev)+1):(varLevOffset(varIndexStddev)+nLevEns) )
      end if
    end do
    call gio_writeToFile(stateVector, './stddev.fst', etiket_in='STDDEV3D', ip3_opt=ip3, &
                         typvar_opt='E', numBits_opt=numBits, containsFullField_opt=.false.)

    ! second, write stddev (zonal average)

    if (mmpi_myid == 0) then
      nulstats=0
      ierr =  fnom  (nulstats,'./stddev.fst','RND',0)
      ierr =  fstouv(nulstats,'RND')

      ! do 3d variables
      do varIndex=1,nvar3d
        nLevEns = gsv_getNumLevFromVarName(stateVector,nomvar(varIndex,cvSpace))

        dfact=1.0d0
        do levIndex = 1, nlevEns
          do latIndex = 1, nj
            zbufyz(latIndex,levIndex)=dfact*stddevZonAvg(varLevOffset(varIndex)+levIndex,latIndex)
          end do
        end do
        ierr = utl_fstecr(zbufyz(:,1:nLevEns),-numBits,nulstats,idateo,0,0,1,nj,nlevEns,ip1,ip2,ip3,   &
                          'E',nomvar3d(varIndex,cvSpace),'STDDEV  ','X',0,0,0,0,idatyp,.true.)

        if (present(stddevZonAvgUnbal_opt) .and. present(stddev3dUnbal_opt)) then
          if(nomvar3d(varIndex,cvSpace).ne.nomvar3d(varIndex,cvUnbalSpace)) then
            dfact=1.0d0
            do levIndex = 1, nlevEns
              do latIndex = 1, nj
                zbufyz(latIndex,levIndex)=dfact*stddevZonAvgUnbal_opt(varLevOffset(varIndex)+levIndex,latIndex)
              end do
            end do
            ierr = utl_fstecr(zbufyz(:,1:nLevEns),-numBits,nulstats,idateo,0,0,1,nj,nlevEns,ip1,ip2,ip3,   &
                              'E',nomvar3d(varIndex,cvUnbalSpace),'STDDEV  ','X',0,0,0,0,idatyp,.true.)
          end if
        end if

      end do

      ! now do 2D variables
      do varIndex=1,nvar2d
        if(nomvar2d(varIndex,cvSpace).eq.'P0') then
          dfact=1.0d0/1.0d2
        else
          dfact=1.0d0
        end if

        do latIndex = 1, nj
          zbufy(latIndex)=dfact*stddevZonAvg(varLevOffset(nvar3d+1)+varIndex,latIndex)
        end do
        ierr = utl_fstecr(zbufy,-numBits,nulstats,idateo,0,0,1,nj,1,ip1,ip2,ip3,   &
                          'E',nomvar2d(varIndex,cvSpace),'STDDEV  ','X',0,0,0,0,idatyp,.true.)

        if (present(stddevZonAvgUnbal_opt) .and. present(stddev3dUnbal_opt)) then
          if(nomvar2d(varIndex,cvSpace).ne.nomvar2d(varIndex,cvUnbalSpace)) then
            if(nomvar2d(varIndex,cvUnbalSpace).eq.'UP') then
              dfact=1.0d0/1.0d2
            else
              dfact=1.0d0
            end if

            do latIndex = 1, nj
              zbufy(latIndex)=dfact*stddevZonAvgUnbal_opt(varLevOffset(nvar3d+1)+varIndex,latIndex)
            end do
            ierr = utl_fstecr(zbufy,-numBits,nulstats,idateo,0,0,1,nj,1,ip1,ip2,ip3,   &
                              'E',nomvar2d(varIndex,cvUnbalSpace),'STDDEV  ','X',0,0,0,0,idatyp,.true.)
          end if
        end if

      end do

      ierr =  fstfrm(nulstats)
      ierr =  fclos (nulstats)
    end if

    call gsv_deallocate(stateVector)

    write(*,*) 'finished writing stddev...'

  end subroutine writeStddev

  !--------------------------------------------------------------------------
  ! WRITESTDDEVBAL
  !--------------------------------------------------------------------------
  subroutine writeStddevBal(stddevZonAvgBal,stddev3dBal)
    !
    !:Purpose: Write the stddev of the balanced variables to a file
    !
    implicit none

    ! Arguments:
    real(8), intent(in) :: stddevZonAvgBal(:,:)
    real(8), intent(in) :: stddev3dBal(:,:,:)

    ! Locals:
    type(struct_gsv) :: stateVector
    real(8) :: dfact, zbufyz(nj,max(nLevEns_M,nLevens_T)), zbufy(nj)
    integer :: latIndex, levIndex, ierr, varIndex, nLevEns
    integer :: fstouv, fnom, fstfrm, fclos
    integer :: ip1, ip2, ip3, idatyp, idateo, numBits
    integer :: nulstats
    real(8), pointer :: field(:,:,:)
    integer, allocatable :: dateStampList(:)
    integer, parameter :: nvar3d=1, nvar2d=1, nvar=nvar3d+nvar2d
    character(len=4) :: nomVarToWrite(nvar)
    character(len=4) :: nomvar3dBal(nvar3d), nomvar2dBal(nvar2d)
    integer          :: varLevOffsetBal(nvar)

    nomvar3dBal(1)='TB'
    nomvar2dBal(1)='PB'
    varLevOffsetBal(1) = 0
    varLevOffsetBal(2) = 1*nLevEns_T

    numBits = 32
    idatyp = 5
    ip1 = 0
    ip2 = 0
    ip3 = nens
    idateo = 0

    nomVarToWrite(1:nvar3d)        = nomvar3dBal(:)
    nomVarToWrite((1+nvar3d):nvar) = nomvar2dBal(:)

    ! first, write stddev3d
    
    ! allocate stateVector used for writing stddev3d
    allocate(dateStampList(1))
    dateStampList(:)  = 0
    call gsv_allocate( stateVector, 1, hco_ens, vco_ens,  &
                       mpi_local_opt=.true., dateStampList_opt=dateStampList,   &
                       varNames_opt=nomVarToWrite(1:nvar) )
    do varIndex = 1, nvar
      nLevEns = gsv_getNumLevFromVarName(stateVector,nomVarToWrite(varIndex))
      call gsv_getField(stateVector,field,nomVarToWrite(varIndex))
      field(:,:,:) = stddev3dBal(:, :, (varLevOffsetBal(varIndex)+1):(varLevOffsetBal(varIndex)+nLevEns) )
    end do
    call gio_writeToFile(stateVector, './stddev_balanced.fst', etiket_in='STDDEV3D', ip3_opt=ip3, &
                         typvar_opt='E', numBits_opt=numBits, containsFullField_opt=.false.)

    ! second, write stddev (zonal average)

    if (mmpi_myid == 0) then
      nulstats=0
      ierr =  fnom  (nulstats,'./stddev_balanced.fst','RND',0)
      ierr =  fstouv(nulstats,'RND')

      ! do 3d variables
      do varIndex=1,nvar3d
        nLevEns = gsv_getNumLevFromVarName(stateVector,nomVar3dBal(varIndex))
        !nip1_l(1:nLevEns_T)=nip1_T(1:nLevEns_T)
        dfact=1.0d0
        do levIndex = 1, nlevEns
          do latIndex = 1, nj
            zbufyz(latIndex,levIndex)=dfact*stddevZonAvgBal(varLevOffsetBal(varIndex)+levIndex,latIndex)
          end do
        end do
        ierr = utl_fstecr(zbufyz(:,1:nLevEns),-numBits,nulstats,idateo,0,0,1,nj,nlevEns,ip1,ip2,ip3,   &
                          'E',nomvar3dBal(varIndex),'STDDEV  ','X',0,0,0,0,idatyp,.true.)
      end do

      ! now do 2D variables
      do varIndex=1,nvar2d
        dfact=1.0d0/1.0d2
        do latIndex = 1, nj
          zbufy(latIndex)=dfact*stddevZonAvgBal(varLevOffsetBal(nvar3d+1)+varIndex,latIndex)
        end do
        ierr = utl_fstecr(zbufy,-numBits,nulstats,idateo,0,0,1,nj,1,ip1,ip2,ip3,   &
                          'E',nomvar2dBal(varIndex),'STDDEV  ','X',0,0,0,0,idatyp,.true.)
      end do

      ierr =  fstfrm(nulstats)
      ierr =  fclos (nulstats)
    end if

    call gsv_deallocate(stateVector)

    write(*,*) 'finished writing stddev...'

  end subroutine writeStddevBal

  !--------------------------------------------------------------------------
  ! spectralFilterLegacy
  !--------------------------------------------------------------------------
  subroutine spectralFilterLegacy(ensPerturbations,nlev)
    !
    !:Purpose: Apply a spectral filter
    !
    implicit none

    ! Arguments:
    real(4), pointer,  intent(inout) :: ensPerturbations(:,:,:,:)
    integer,           intent(in)    :: nlev

    ! Locals:
    real(8) :: spectralState(nla_mpilocal,2,nlev)
    real(8) :: member(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nlev)
    integer :: ensIndex

    do ensIndex=1,nens
      member(:,:,:)=dble(ensPerturbations(:,:,:,ensIndex))     
      if(nlev.eq.nkgdimEns) then
        call gst_setID(gstID_nkgdimEns)
      else if(nlev.eq.nLevEns_T+1) then
        call gst_setID(gstID_nLevEns_T_P1)
      else
        write(*,*) 'spectralFilter: nlev = ',nlev
        call utl_abort('spectralFilter: spectral transform not initialized for this number of levels')
      end if
      call gst_reespe(spectralState,member)
      call gst_speree(spectralState,member)
      ensPerturbations(:,:,:,ensIndex)=sngl(member(:,:,:))
    end do

    write(*,*) 'finished applying spectral filter...'

  end subroutine spectralFilterLegacy

  !--------------------------------------------------------------------------
  ! spectralFilter
  !--------------------------------------------------------------------------
  subroutine spectralFilter(ensPerts)
    !
    !:Purpose: Apply a spectral filter
    !
    implicit none

    ! Arguments:
    type(struct_ens),  intent(inout) :: ensPerts

    ! Locals:
    real(8), allocatable :: ensPertSP(:,:,:)
    real(8), allocatable :: ensPertGD(:,:,:)
    real(4), pointer     :: ptr4d_r4(:,:,:,:)    
    integer :: memberIndex, levIndex, latIndex, lonIndex
    integer :: gstFilterID

    !
    !- 2.  Apply Filter
    !
    gstFilterID = gst_setup(ni,nj,nTrunc,nEnsOverDimension)

    allocate(ensPertSP(nla_mpilocal,2,nEnsOverDimension))
    allocate(ensPertGD(nEnsOverDimension,myLonBeg:myLonEnd,myLatBeg:myLatEnd))
    ensPertSP(:,:,:) = 0.0d0
    
    do levIndex = 1, ens_getNumK(ensPerts) ! Loop on variables and vertical levels

      ptr4d_r4 => ens_getOneLev_r4(ensPerts,levIndex)

      do latIndex = myLatBeg, myLatEnd
        ensPertGD(:,:,latIndex) = 0.0d0
      end do

      !$OMP PARALLEL DO PRIVATE (memberIndex,latIndex,lonIndex)
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          do memberIndex = 1, nEns
            ensPertGD(memberIndex,lonIndex,latIndex) = dble(ptr4d_r4(memberIndex,1,lonIndex,latIndex))
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      !- GridPoint space -> Spectral Space
      call gst_setID(gstFilterID) ! IN
      call gst_reespe_kij(ensPertSP, & ! OUT
                          ensPertGD)   ! IN

      ! Spectral Space -> GridPoint space
      call gst_setID(gstFilterID) ! IN
      call gst_speree_kij(ensPertSP, & ! IN
                          ensPertGD)   ! OUT

      !$OMP PARALLEL DO PRIVATE (memberIndex,latIndex,lonIndex)
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          do memberIndex = 1, nEns
            ptr4d_r4(memberIndex,1,lonIndex,latIndex) = sngl(ensPertGD(memberIndex,lonIndex,latIndex))
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    end do

    deallocate(ensPertGD)
    deallocate(ensPertSP)
    
    write(*,*) 'finished applying spectral filter...'

  end subroutine spectralFilter
  
  !--------------------------------------------------------------------------
  ! CALCTHETA
  !--------------------------------------------------------------------------
  subroutine calcTheta(ensPerturbations,theta)
    !
    !:Purpose: Calculate the Theta turning angle according to Ekman balance
    !
    implicit none

    ! Arguments:
    real(4), pointer, intent(in)  :: ensPerturbations(:,:,:,:)
    real(8),          intent(out) :: theta(:,:)

    ! Locals:
    real(8) :: zchipsi(nLevEns_M,nj), zpsipsi(nLevEns_M,nj)
    real(8) :: zchipsi_mpiglobal(nLevEns_M,nj), zpsipsi_mpiglobal(nLevEns_M,nj)
    real(4), pointer :: psi_ptr(:,:,:), chi_ptr(:,:,:)
    integer :: latIndex,lonIndex,levIndex,ensIndex, ierr, nsize

    theta(:,:) = 0.0d0
    zchipsi(:,:) = 0.0d0
    zpsipsi(:,:) = 0.0d0
    zchipsi_mpiglobal(:,:) = 0.0d0
    zpsipsi_mpiglobal(:,:) = 0.0d0

    do ensIndex = 1,nens
      psi_ptr(myLonBeg:,myLatBeg:,1:) => ensPerturbations(:,:,1:nlevEns_M,ensIndex)
      chi_ptr(myLonBeg:,myLatBeg:,1:) => ensPerturbations(:,:,(nlevEns_M+1):(2*nlevEns_M),ensIndex)

      ! update zchipsi and zpsipsi covariances
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          do levIndex = 1, nLevEns_M
            zpsipsi(levIndex,latIndex) = zpsipsi(levIndex,latIndex) +  &
                 psi_ptr(lonIndex,latIndex,levIndex) * psi_ptr(lonIndex,latIndex,levIndex)
            zchipsi(levIndex,latIndex) = zchipsi(levIndex,latIndex) +  &
                 chi_ptr(lonIndex,latIndex,levIndex) * psi_ptr(lonIndex,latIndex,levIndex)
          end do
        end do
      end do
    end do

    nsize = nLevEns_M*nj
    call rpn_comm_allreduce(zchipsi,zchipsi_mpiglobal,nsize,"mpi_double_precision","mpi_sum","GRID",ierr)
    call rpn_comm_allreduce(zpsipsi,zpsipsi_mpiglobal,nsize,"mpi_double_precision","mpi_sum","GRID",ierr)

    !  calculate THETA
    do latIndex = 1, nj
      do levIndex = 1, nLevEns_M
        theta(levIndex,latIndex) = atan(-zchipsi_mpiglobal(levIndex,latIndex) /  &
                                         zpsipsi_mpiglobal(levIndex,latIndex))
      end do
    end do

    write(*,*) 'finished computing theta...'

  end subroutine calcTheta

  !--------------------------------------------------------------------------
  ! CALCPTOT
  !--------------------------------------------------------------------------
  subroutine calcPtoT(ensPerturbations,PtoT)
    !
    !:Purpose: Calculate the "P" to Temperature transform matrix using
    !          a regression analysis of the "P" and temperature samples
    !
    implicit none

    ! Arguments:
    real(4), pointer, intent(in)  :: ensPerturbations(:,:,:,:)
    real(8),          intent(out) :: PtoT(:,:,:)

    ! Locals:
    real(8) :: spectralState(nla_mpilocal,2,nLevEns_M), spBalancedP(nla_mpilocal,2,nlevEns_M)
    real(8) :: balancedP(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nlevEns_M)
    real(8) :: psi(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nLevEns_M)
    real(4), pointer :: tt_ptr(:,:,:),ps_ptr(:,:)
    INTEGER :: ensIndex, JK1, JK2, nsize
    INTEGER :: IERR, JK, latIndex, lonIndex, JB, JPNLATBND
    PARAMETER (JPNLATBND = 3)
    REAL(8) :: ZFACT,zlat(nj)
    REAL(8) :: ZFACTTOT
    REAL(8) :: ZM1(NLEVENS_T+1,NLEVENS_M,JPNLATBND), ZM2(NLEVPTOT,NLEVPTOT,JPNLATBND)
    REAL(8) :: ZM1_mpiglobal(NLEVENS_T+1,NLEVENS_M,JPNLATBND), ZM2_mpiglobal(NLEVPTOT,NLEVPTOT,JPNLATBND)
    REAL(8) :: ZPTOTBND(NLEVENS_T+1,NLEVENS_M)
    REAL(8) :: ZM2INV(NLEVPTOT,NLEVPTOT,JPNLATBND)
    REAL(8) :: DLA2, DL1SA2
    REAL(8) :: DLLATMIN(JPNLATBND), DLLATMAX(JPNLATBND)
    real(8) :: zeigwrk(4*nlevPtoT),zeigen(nlevPtoT,nlevPtoT),zeigenv(nlevPtoT)
    real(8) :: zeigenvi(nlevPtoT)
    integer :: iwork,info

    DATA DLLATMIN / -60.0D0, -30.0D0, 30.0D0 /
    DATA DLLATMAX / -30.0D0,  30.0D0, 60.0D0 /

    DLA2 = ec_ra * ec_ra
    DL1SA2 = 1.D0/DLA2

    ! 1. Initialize P_to_T, ZM1, ZM2

    ZFACTTOT = 0.0D0
    DO latIndex = 1, nj
      ZFACTTOT = ZFACTTOT + cos(GST_GETRLATI(latIndex))
    END DO
    ZFACTTOT = NJ/ZFACTTOT

    PtoT(:,:,:) = 0.0d0
    ZM1(:,:,:) = 0.0d0
    ZM2(:,:,:) = 0.0d0
    ZPTOTBND(:,:) = 0.0d0
    ZM1_mpiglobal(:,:,:) = 0.0d0
    ZM2_mpiglobal(:,:,:) = 0.0d0

    do ensIndex = 1, nens

      write(*,*) 'calcPtoT: processing member ',ensIndex

      psi(:,:,:) = ensPerturbations(:,:,1:nlevEns_M,ensIndex)
      call gst_setID(gstID_nLevEns_M)
      call gst_reespe(spectralState,psi)
      call calcBalancedP(spectralState,spBalancedP)
      call gst_speree(spBalancedP,balancedP)

      tt_ptr(myLonBeg:,myLatBeg:,1:) => ensPerturbations(:,:,(2*nLevEns_M+1):(2*nLevEns_M+nLevEns_T),ensIndex)
      ps_ptr(myLonBeg:,myLatBeg:)    => ensPerturbations(:,:,2*nLevEns_M+2*nLevEns_T+1,ensIndex)

      do latIndex = myLatBeg, myLatEnd
        zlat(latIndex)=GST_GETRLATI(latIndex)
      end do

      !$OMP PARALLEL DO PRIVATE (JK1,JB,latIndex,ZFACT,lonIndex,JK2)
      DO JK1 = 1, (nLevEns_T+1)
        DO JB=1,JPNLATBND
          DO latIndex = myLatBeg, myLatEnd
            if ((ZLAT(latIndex) .gt. 2.D0*MPC_PI_R8*DLLATMIN(JB)/360.D0) .and.  &
                (ZLAT(latIndex) .le. 2.D0*MPC_PI_R8*DLLATMAX(JB)/360.D0)) then
              ZFACT = cos(ZLAT(latIndex))*ZFACTTOT
              DO lonIndex = myLonBeg, myLonEnd

                ! update ZM1 = sum_over_t_x_y[vec(T lnPs) vec(P_b)^T]
                DO JK2 = 1, nLevEns_M
                  IF(JK1.LE.nLevEns_T) THEN
                    zm1(jk1,jk2,jb) = zm1(jk1,jk2,jb) + zfact * tt_ptr(lonIndex,latIndex,jk1) * balancedP(lonIndex,latIndex,jk2)
                  ELSE
                    zm1(jk1,jk2,jb) = zm1(jk1,jk2,jb) + zfact * ps_ptr(lonIndex,latIndex) * balancedP(lonIndex,latIndex,jk2)
                  END IF
                END DO

                ! update ZM2 = sum_over_t_x_y[vec(P_b) vec(P_b)^T]
                IF(JK1.LE.NLEVPTOT) THEN
                  DO JK2 = 1, NLEVPTOT
                    zm2(jk1,jk2,jb) = zm2(jk1,jk2,jb) + zfact * balancedP(lonIndex,latIndex,jk1) * balancedP(lonIndex,latIndex,jk2)
                  END DO
                END IF

              END DO
            end if
          END DO
        END DO ! Loop on JPNLATBND
      END DO  ! Loop on JK1
      !$OMP END PARALLEL DO

    END DO

    ! communicate matrices to have global result on all tasks
    nsize = (NLEVENS_T+1)*NLEVENS_M*JPNLATBND
    call rpn_comm_allreduce(zm1,zm1_mpiglobal,nsize,"mpi_double_precision","mpi_sum","GRID",ierr)
    nsize = NLEVPTOT*NLEVPTOT*JPNLATBND
    call rpn_comm_allreduce(zm2,zm2_mpiglobal,nsize,"mpi_double_precision","mpi_sum","GRID",ierr)
    
    ! SET ZM1_MPIGLOBAL, ZM2_MPIGLOBAL EQUAL FOR ALL THREE REGIONS
    DO JK1 = 1, NLEVPTOT
      DO JK2 = 1, NLEVPTOT
        ZM2_MPIGLOBAL(JK1,JK2,1)=ZM2_MPIGLOBAL(JK1,JK2,1)+ZM2_MPIGLOBAL(JK1,JK2,3)
        ZM2_MPIGLOBAL(JK1,JK2,2)=ZM2_MPIGLOBAL(JK1,JK2,1)
        ZM2_MPIGLOBAL(JK1,JK2,3)=ZM2_MPIGLOBAL(JK1,JK2,1)
      END DO
    END DO
    DO JK1 = 1, (nLevEns_T+1)
      DO JK2 = 1, NLEVPTOT
        ZM1_MPIGLOBAL(JK1,JK2,1)=ZM1_MPIGLOBAL(JK1,JK2,1)+ZM1_MPIGLOBAL(JK1,JK2,3)
        ZM1_MPIGLOBAL(JK1,JK2,2)=ZM1_MPIGLOBAL(JK1,JK2,1)
        ZM1_MPIGLOBAL(JK1,JK2,3)=ZM1_MPIGLOBAL(JK1,JK2,1)
      END DO
    END DO

    DO JK1=1,NLEVPTOT
      DO JK2=1,NLEVPTOT
        ZEIGEN(JK1,JK2)=ZM2_MPIGLOBAL(JK1,JK2,1)
      END DO
    END DO
    IWORK=4*NLEVPTOT
    CALL DSYEV('V','U',NLEVPTOT,ZEIGEN,NLEVPTOT,ZEIGENV,ZEIGWRK,IWORK,INFO)

    write(*,*) 'calcPtot: info=',info
    write(*,*) 'calcPtot: eigen values=',zeigenv(:)

    do JK1=1,NLEVPTOT
      if (ZEIGENV(JK1).gt.0.0d0) then
        ZEIGENVI(JK1)=1.0d0/ZEIGENV(JK1)
      else
        ZEIGENVI(JK1)=0.0d0
      end if
    end do

    DO JK1=1,NLEVPTOT
      DO JK2=1,NLEVPTOT
        ZM2INV(JK1,JK2,1)=0.0d0
        DO JK=1,NLEVPTOT
          ZM2INV(JK1,JK2,1)=ZM2INV(JK1,JK2,1)+ZEIGEN(JK1,JK)*ZEIGENVI(JK)*ZEIGEN(JK2,JK)
        END DO
      END DO
    END DO


    ! Calculate A = ZM1_MPIGLOBAL*inv(ZM2)
    DO JK1 = 1, (nLevEns_T+1)
      DO JK2 = 1, NLEVPTOT
        DO JK = 1, NLEVPTOT
          ZPTOTBND(JK1,JK2) = ZPTOTBND(JK1,JK2) + ZM1_MPIGLOBAL(JK1,JK,1) * ZM2INV(JK,JK2,1)
        END DO
      END DO
    END DO

    DO JK1 = 1, nLevEns_T+1
      DO JK2 = 1, NLEVPTOT
        DO latIndex = 1, nj
          PTOT(JK1,JK2,latIndex) = ZPTOTBND(JK1,JK2)
        END DO
      END DO
    END DO

    write(*,*) 'finished computing PtoT...'

  end subroutine calcPtoT

  !--------------------------------------------------------------------------
  ! REMOVEGLOBALMEAN
  !--------------------------------------------------------------------------
  subroutine removeGlobalMean(ensPerturbations)
    !
    !:Purpose: Calculate and subtract the horizontal mean value from the ensemble
    !          of samples for each member, vertical level and variable
    !
    implicit none

    ! Arguments:
    real(4), pointer, intent(inout) :: ensPerturbations(:,:,:,:)

    ! Locals:
    integer :: lonIndex,latIndex,levIndex,ensIndex,ierr
    real(8)  :: dmean, dmean_mpiglobal

    do ensIndex = 1, nens
      do levIndex = 1, nkgdimEns
        dmean = 0.0d0
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            dmean = dmean + ensPerturbations(lonIndex,latIndex,levIndex,ensIndex)
          end do
        end do
        call rpn_comm_allreduce(dmean, dmean_mpiglobal,1,"mpi_double_precision","mpi_sum","GRID",ierr)
        dmean_mpiglobal = dmean_mpiglobal/(dble(ni)*dble(nj))
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            ensPerturbations(lonIndex,latIndex,levIndex,ensIndex) =   &
                 ensPerturbations(lonIndex,latIndex,levIndex,ensIndex) - dmean_mpiglobal
          end do
        end do
      end do
    end do

    write(*,*) 'finished removing global mean...'

  end subroutine removeGlobalMean

  !--------------------------------------------------------------------------
  ! CALCZONAVG
  !--------------------------------------------------------------------------
  subroutine calcZonAvg(fieldsZonAvg_mpiglobal,fields3D,nlev)
    !
    !:Purpose: Calculate the zonal average of the supplied 3D fields
    !
    implicit none

    ! Arguments:
    real(8), pointer, intent(inout) :: fieldsZonAvg_mpiglobal(:,:)
    real(8), pointer, intent(in)    :: fields3D(:,:,:)
    integer,          intent(in)    :: nlev

    ! Locals:
    integer :: lonIndex, latIndex, levIndex, ierr, nsize
    real(8) :: dfact
    real(8), allocatable :: fieldsZonAvg(:,:) 

    allocate(fieldsZonAvg(nlev,nj))
    fieldsZonAvg(:,:)=0.0d0

    dfact=1.0d0/dble(ni)
    !$OMP PARALLEL DO PRIVATE (levIndex,latIndex,lonIndex)
    do levIndex = 1, nlev
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          fieldsZonAvg(levIndex,latIndex) = fieldsZonAvg(levIndex,latIndex) +  &
                                            dfact*fields3D(lonIndex,latIndex,levIndex)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    ! combine info from all mpi tasks
    nsize = nlev*nj
    call rpn_comm_allreduce(fieldsZonAvg,fieldsZonAvg_mpiglobal,nsize,"mpi_double_precision","mpi_sum","GRID",ierr)

    deallocate(fieldsZonAvg)

    write(*,*) 'finished computing the zonal average...'

  end subroutine calcZonAvg

  !--------------------------------------------------------------------------
  ! CALCSTDDEV3D
  !--------------------------------------------------------------------------
  subroutine calcStddev3d(ensPerturbations,stddev3d,nlev)
    !
    !:Purpose: Calculate the 3d stddev field
    !
    implicit none

    ! Arguments:
    real(8), pointer, intent(inout) :: stddev3d(:,:,:)
    real(4), pointer, intent(in)    :: ensPerturbations(:,:,:,:)
    integer,          intent(in)    :: nlev

    ! Locals:
    integer :: lonIndex,latIndex,levIndex,ensIndex
    real(8) :: dnens

    write(*,*) 'started computing the stddev...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    stddev3d(:,:,:) = 0.0d0
    dnens = 1.0d0/dble(nens-1)
    !$OMP PARALLEL DO PRIVATE (levIndex,ensIndex,latIndex,lonIndex)
    do levIndex = 1, nlev
      do ensIndex = 1, nens
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            stddev3d(lonIndex,latIndex,levIndex)=stddev3d(lonIndex,latIndex,levIndex)+ensPerturbations(lonIndex,latIndex,levIndex,ensIndex)**2
          end do
        end do
      end do
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          if(stddev3d(lonIndex,latIndex,levIndex).gt.0.0d0) then
            stddev3d(lonIndex,latIndex,levIndex)=sqrt(stddev3d(lonIndex,latIndex,levIndex)*dnens)
          else
            stddev3d(lonIndex,latIndex,levIndex)=0.0d0
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    write(*,*) 'finished computing the stddev...'
  
  end subroutine calcStddev3d

  !--------------------------------------------------------------------------
  ! CALCBALANCEDP
  !--------------------------------------------------------------------------
  subroutine calcBalancedP(sppsi,spgz)
    !
    !:Purpose: Calculate the balanced "P" variable using geostrophy
    !
    implicit none

    ! Arguments:
    real(8), intent(in)  :: sppsi(:,:,:)
    real(8), intent(out) :: spgz(:,:,:)

    ! Locals:
    real(8) :: spvor_mpiglobal(nla,2,nlevEns_M)
    real(8) :: spvor_mpiglobal2(nla,2,nlevEns_M)
    real(8) :: spgz_mpiglobal(nla,2,nlevEns_M)
    integer :: ia, ib, ji, jm, levIndex, jla_mpilocal, ila_mpiglobal
    integer :: ierr, nsize
    real(8) :: zn, zm, zenm, zenmp1, zcon, dl1sa2
    
    ! convert PSI to vorticity 
    dl1sa2   = 1.0d0/(ec_ra*ec_ra)
    spvor_mpiglobal(:,:,:) = 0.0d0
    spvor_mpiglobal2(:,:,:) = 0.0d0
    do levIndex = 1, nlevEns_M
      do jla_mpilocal = 1, nla_mpilocal
        ila_mpiglobal = ilaList_mpiglobal(jla_mpilocal)
        spvor_mpiglobal(ila_mpiglobal,1,levIndex) = sppsi(jla_mpilocal,1,levIndex)*dl1sa2*gst_getRnnp1(ila_mpiglobal)
        spvor_mpiglobal(ila_mpiglobal,2,levIndex) = sppsi(jla_mpilocal,2,levIndex)*dl1sa2*gst_getRnnp1(ila_mpiglobal)
      end do
      ! ensure input field is zero for spectral component (0,0)
      spvor_mpiglobal(1,1,levIndex) = 0.0D0
      spvor_mpiglobal(1,2,levIndex) = 0.0D0
    end do

    nsize = nla*2*nlevEns_M
    call rpn_comm_allreduce(spvor_mpiglobal,spvor_mpiglobal2,nsize,"mpi_double_precision","mpi_sum","GRID",ierr)
    
    ! initialize output field to zero
    spgz_mpiglobal(:,:,:)=0.0d0

    ! loop over levels and zonal wavenumbers
    ! n.b.: at the tip of the triangle, no contributions
    
    zcon = -2.D0*ec_ROmega*ec_ra**2
    do levIndex = 1, nlevEns_M

      ! the base address ia will point to the spherical harmonic
      ! coefficient (m,m), in the input field
      ia = 1
      do jm = 0, ntrunc-1
        ib = ia + ntrunc - jm
        zm = dble(jm)

        ! at the base, contributions from n+1 coeff only
        zn = zm
        zenmp1 = sqrt ( ((zn+1)**2-zm**2)/(4.D0*(zn+1)**2-1.D0) )
        spgz_mpiglobal(ia,1,levIndex)=zcon*spvor_mpiglobal2(ia+1,1,levIndex)*zenmp1/((zn+1.0D0)**2)
        spgz_mpiglobal(ia,2,levIndex)=zcon*spvor_mpiglobal2(ia+1,2,levIndex)*zenmp1/((zn+1.0D0)**2)

        zn = zn+1
        do ji = ia+1, ib-1
          zenm = sqrt ( (zn**2-zm**2)/(4.D0*zn**2-1.D0) )
          zenmp1 = sqrt ( ((zn+1)**2-zm**2)/(4.D0*(zn+1)**2-1.D0) )
          spgz_mpiglobal(ji,1,levIndex)=spvor_mpiglobal2(ji-1,1,levIndex)*zenm/(zn**2)
          spgz_mpiglobal(ji,2,levIndex)=spvor_mpiglobal2(ji-1,2,levIndex)*zenm/(zn**2)
          spgz_mpiglobal(ji,1,levIndex)=zcon*(spgz_mpiglobal(ji,1,levIndex)+spvor_mpiglobal2(ji+1,1,levIndex)*zenmp1/((zn+1.0D0)**2))
          spgz_mpiglobal(ji,2,levIndex)=zcon*(spgz_mpiglobal(ji,2,levIndex)+spvor_mpiglobal2(ji+1,2,levIndex)*zenmp1/((zn+1.0D0)**2))
          zn = zn + 1.0D0
        end do

        ! at the top, contributions from n-1 coeff only
        zenm = sqrt ( (zn**2-zm**2)/(4.D0*zn**2-1.D0) )
        spgz_mpiglobal(ib,1,levIndex) = zcon*spvor_mpiglobal2(ib-1,1,levIndex)*zenm/(zn**2)
        spgz_mpiglobal(ib,2,levIndex) = zcon*spvor_mpiglobal2(ib-1,2,levIndex)*zenm/(zn**2)
        ia = ib + 1
      end do
    end do

    ! ensure correct value for mass spectral-coefficient for m=n=0
    do levIndex = 1, nlevens_M
      spgz_mpiglobal(1,1,levIndex) = 0.0D0
      spgz_mpiglobal(1,2,levIndex) = 0.0D0
    end do

    ! copy to mpilocal output array
    do levIndex = 1, nlevEns_M
      do jla_mpilocal = 1, nla_mpilocal
        ila_mpiglobal = ilaList_mpiglobal(jla_mpilocal)
        spgz(jla_mpilocal,1,levIndex) = spgz_mpiglobal(ila_mpiglobal,1,levIndex)
        spgz(jla_mpilocal,2,levIndex) = spgz_mpiglobal(ila_mpiglobal,2,levIndex)
      end do
    end do

  end subroutine calcBalancedP

  !--------------------------------------------------------------------------
  ! NORMALIZED3D
  !--------------------------------------------------------------------------
  subroutine normalize3d(ensPerturbations,stddev3d)
    !
    !:Purpose: Divide the ensemble perturbations by the supplied 3d stddev field
    !
    implicit none

    ! Arguments:
    real(8), pointer, intent(in)    :: stddev3d(:,:,:)
    real(4), pointer, intent(inout) :: ensPerturbations(:,:,:,:)

    ! Locals:
    integer :: lonIndex,latIndex,levIndex,ensIndex
    real(8) :: dfact

    !$OMP PARALLEL DO PRIVATE (levIndex,ensIndex,latIndex,lonIndex,DFACT)
    do levIndex = 1, nkgdimEns
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          if(stddev3d(lonIndex,latIndex,levIndex).gt.0.0d0) then
            dfact=1.0d0/stddev3d(lonIndex,latIndex,levIndex)
          else
            dfact=0.0d0
          end if
          do ensIndex = 1, nens
            ensPerturbations(lonIndex,latIndex,levIndex,ensIndex)=ensPerturbations(lonIndex,latIndex,levIndex,ensIndex)*dfact
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    write(*,*) 'finished normalizing by stddev3D...'
  
  end subroutine normalize3d

  !--------------------------------------------------------------------------
  ! MULTIPLY3D
  !--------------------------------------------------------------------------
  subroutine multiply3d(ensPerturbations,stddev3d,nlev)
    !
    !:Purpose: Multiply the ensemble perturbations by the supplied 3d stddev field
    !
    implicit none

    ! Arguments:
    real(8), pointer, intent(in)    :: stddev3d(:,:,:)
    real(4), pointer, intent(inout) :: ensPerturbations(:,:,:,:)
    integer,          intent(in)    :: nlev

    ! Locals:
    integer :: lonIndex,latIndex,levIndex,ensIndex

    !$OMP PARALLEL DO PRIVATE (levIndex,ensIndex,latIndex,lonIndex)
    do ensIndex = 1, nens
      do levIndex = 1, nlev
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            ensPerturbations(lonIndex,latIndex,levIndex,ensIndex)=ensPerturbations(lonIndex,latIndex,levIndex,ensIndex)*stddev3d(lonIndex,latIndex,levIndex)
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    write(*,*) 'finished multiplying by stddev3D...'
  
  end subroutine multiply3d

  !--------------------------------------------------------------------------
  ! READENSEMBLE
  !--------------------------------------------------------------------------
  subroutine readEnsemble(ensPerturbations)
    !
    !:Purpose: Read the ensemble of error samples from files
    !
    implicit none

    ! Arguments:
    real(4), pointer, intent(inout) :: ensPerturbations(:,:,:,:)

    ! Locals:
    integer :: lonIndex, latIndex, levIndex, ensIndex, numStep
    integer, allocatable :: dateStampList(:)
    real(4), pointer :: field_r4(:,:,:)
    logical :: makeBiPeriodic
    type(struct_ens) :: ensPerts
    type(struct_gsv) :: stateVector

    write(*,*) 'Before reading the ensemble:'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    numStep = 1
    allocate(dateStampList(numStep))
    dateStampList(:)  = -1
    call ens_allocate(ensPerts, nEns, numStep, hco_ens, vco_ens, dateStampList)

    makeBiPeriodic = .false.
    call ens_readEnsemble(ensPerts, './ensemble', makeBiPeriodic, &
                          containsFullField_opt=.false.)

    call gsv_allocate(stateVector, 1, hco_ens, vco_ens, dateStampList_opt=dateStampList,  &
                      mpi_local_opt=.true., dataKind_opt=4)

    do ensIndex = 1, nEns
      write(*,*) 'readEnsemble: copying over member ', ensIndex
      call ens_copyMember(ensPerts, stateVector, ensIndex)

      call gsv_getField(stateVector,field_r4,'UU')
      do levIndex = 1, nLevEns_M
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            ensPerturbations(lonIndex,latIndex,levIndex+varLevOffset(1),ensIndex)= field_r4(lonIndex,latIndex,levIndex)
          end do
        end do
      end do

      call gsv_getField(stateVector,field_r4,'VV')
      do levIndex = 1, nLevEns_M
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            ensPerturbations(lonIndex,latIndex,levIndex+varLevOffset(2),ensIndex)= field_r4(lonIndex,latIndex,levIndex)
          end do
        end do
      end do

      call gsv_getField(stateVector,field_r4,'TT')
      do levIndex = 1, nLevEns_T
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            ensPerturbations(lonIndex,latIndex,levIndex+varLevOffset(3),ensIndex)= field_r4(lonIndex,latIndex,levIndex)
          end do
        end do
      end do
    
      call gsv_getField(stateVector,field_r4,'LQ')
      do levIndex=1,nLevEns_T
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            ensPerturbations(lonIndex,latIndex,levIndex+varLevOffset(4),ensIndex) = field_r4(lonIndex,latIndex,levIndex)
          end do
        end do
      end do

      call gsv_getField(stateVector,field_r4,'P0')
      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          ensPerturbations(lonIndex,latIndex,1+varLevOffset(5),ensIndex)= field_r4(lonIndex,latIndex,1)
        end do
      end do

    end do

    call gsv_deallocate(stateVector)
    call ens_deallocate(ensPerts)
    
    write(*,*) 'After reading the ensemble:'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    write(*,*) 'finished reading ensemble members...'

  end subroutine readEnsemble

  !--------------------------------------------------------------------------
  ! UV_TO_PSICHI
  !--------------------------------------------------------------------------
  subroutine uv_to_psichi(ensPerturbations)
    !
    !:Purpose: Transform wind components to Psi and Chi
    !
    implicit none

    ! Arguments:
    real(4), intent(inout) :: ensPerturbations(:,:,:,:)

    ! Locals:
    integer :: ensIndex, levIndex, jla_mpilocal, ila_mpiglobal
    real(8) :: dla2
    real(8) :: spectralState(nla_mpilocal,2,nkgdimEns)
    real(8) :: member(myLonBeg:myLonEnd,myLatBeg:myLatend,nkgdimens)

    ! Convert from U/V to PSI/CHI and spectrally filter all fields
    call utl_tmg_start(121,'--CalcStats_UVtoPsiChi')
    dla2   = ec_ra * ec_ra
    do ensIndex=1,nens
      write(*,*) '  doing u/v -> psi/chi and spectral filter for member ', ensIndex
      member(:,:,:)=dble(ensPerturbations(:,:,:,ensIndex))
      call gst_setID(gstID_nkgdimEns)
      call gst_gdsp(spectralState,member,nlevEns_M)
      do levIndex = 1, nlevEns_M
        do jla_mpilocal = 1, nla_mpilocal
          ila_mpiglobal = ilaList_mpiglobal(jla_mpilocal)
          spectralState(jla_mpilocal,1,levIndex)           = spectralState(jla_mpilocal,1,levIndex)           * dla2*gst_getR1snp1(ila_mpiglobal)
          spectralState(jla_mpilocal,2,levIndex)           = spectralState(jla_mpilocal,2,levIndex)           * dla2*gst_getR1snp1(ila_mpiglobal)
          spectralState(jla_mpilocal,1,levIndex+nlevEns_M) = spectralState(jla_mpilocal,1,levIndex+nlevEns_M) * dla2*gst_getR1snp1(ila_mpiglobal)
          spectralState(jla_mpilocal,2,levIndex+nlevEns_M) = spectralState(jla_mpilocal,2,levIndex+nlevEns_M) * dla2*gst_getR1snp1(ila_mpiglobal)
        end do
      end do
      call gst_speree(spectralState,member)
      ensPerturbations(:,:,:,ensIndex)=sngl(member(:,:,:))
    end do

    call utl_tmg_stop(121)
    write(*,*) 'finished doing u/v -> psi/chi and spectral filter...'
    
  end subroutine uv_to_psichi

  !--------------------------------------------------------------------------
  ! REMOVEMEAN
  !--------------------------------------------------------------------------
  subroutine removeMean(ensPerturbations)
    !
    !:Purpose: Compute and subtract the ensemble mean 
    !
    implicit none

    ! Arguments:
    real(4), pointer, intent(inout) :: ensPerturbations(:,:,:,:)

    ! Locals:
    integer :: ensIndex, levIndex, latIndex, lonIndex
    real(8) :: dnens, gd2d(myLonBeg:myLonEnd,myLatBeg:myLatEnd)

    ! remove mean and divide by sqrt(2*(NENS-1)) - extra 2 is needed?
    dnens=1.0d0/dble(nens)
      !$OMP PARALLEL DO PRIVATE (levIndex,gd2d,ensIndex,latIndex,lonIndex)
      do levIndex = 1, nkgdimEns
        gd2d(:,:)=0.0d0
        do ensIndex = 1, nens
          do latIndex = myLatBeg, myLatEnd
            do lonIndex = myLonBeg, myLonEnd
              gd2d(lonIndex,latIndex)=gd2d(lonIndex,latIndex)+ensPerturbations(lonIndex,latIndex,levIndex,ensIndex)
            end do
          end do
        end do
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            gd2d(lonIndex,latIndex)=gd2d(lonIndex,latIndex)*dnens
          end do
        end do
        do ensIndex = 1, nens
          do latIndex = myLatBeg, myLatEnd
            do lonIndex = myLonBeg, myLonEnd
              ensPerturbations(lonIndex,latIndex,levIndex,ensIndex)=     &
                ensPerturbations(lonIndex,latIndex,levIndex,ensIndex)-gd2d(lonIndex,latIndex)
            end do
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    write(*,*) 'finished removing the ensemble mean...'

  end subroutine removeMean

  !--------------------------------------------------------------------------
  ! HORIZCORRELFUNCTION
  !--------------------------------------------------------------------------
  subroutine horizCorrelFunction(rstddev,variableType, waveBandIndex_opt)
    !
    !:Purpose: Compute homogeneous-isotropic horizontal correlation
    !          function from spectral variances
    !
    implicit none

    ! Arguments:
    real(8),          intent(in) :: rstddev(nkgdimEns,0:ntrunc)
    integer,          intent(in) :: variableType
    integer,optional, intent(in) :: waveBandIndex_opt

    ! Locals:
    real(8)  :: spectralState(nla,2,nkgdimEns)
    real(8)  :: gridState(ni,nj,nkgdimEns)
    integer :: ji, jk, jn, jm, ila, iref, jref
    integer :: nLevEns, nLevStart, nLevEnd, varIndex, iStart, iEnd
    character(len=128) :: outfilename
    character(len=2) :: wbnum
    
    write(*,*)
    write(*,*) 'Computing horizontal correlation functions'

    !
    !- 1.  Spectral transform of a delta function (at the center of the domain)
    !

    !- 1.1 Create the delta function
    iref = ni/2
    jref = nj/2

    GridState(:,:,:)       = 0.d0
    GridState(iref,jref,:) = 1.d0

    !- 1.2 Adjoint of the identity (change of norm)
    GridState(iref,jref,:) = GridState(iref,jref,:) * real(ni,8) / gst_getRWT(jref)

    !- 1.3 Move to spectral space
    call gst_setID(gstID_nkgdimEns)
    call gst_reespe(spectralState,  & ! OUT
                    GridState)        ! IN

    !
    !- 2.  Apply the horizontal correlation function
    !
    !$OMP PARALLEL DO PRIVATE (jk,jn,jm,ila)
    do jk = 1, nkgdimEns
       do jn = 0, ntrunc
          do jm = 0, jn
             ila = gst_getNind(jm) + jn - jm
             if (jm.eq.0) then
                spectralState(ila,1,jk) = spectralState(ila,1,jk) * rstddev(jk,jn)**2 * 2.0d0
                spectralState(ila,2,jk) = 0.d0
             else
                spectralState(ila,1,jk) = spectralState(ila,1,jk) * rstddev(jk,jn)**2 * 2.0d0
                spectralState(ila,2,jk) = spectralState(ila,2,jk) * rstddev(jk,jn)**2 * 2.0d0
             end if
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !
    !- 3.  Move back to grid point space
    !
    call gst_setID(gstID_nkgdimEns)
    call gst_speree(spectralState,     & ! IN
                    GridState)           ! OUT

    !
    !- 4.  Write to file
    !
    if (mmpi_myid == 0) then
      if ( nHorizWaveBand /= 1 ) then
        if (.not. present(waveBandIndex_opt)) then
          write(*,*) 'horizCorrelFunction: No waveBandIndex was supplied!!!'
          call utl_abort('calbmatrix_glb')
        end if
        write(wbnum,'(I2.2)') waveBandIndex_opt
      end if
      
      !- 4.1 2D correlation function in fst format
      if ( nHorizWaveBand == 1 ) then
        outfilename = "./horizCorrel.fst"
      else
        outfilename = "./horizCorrel_"//wbnum//".fst"
      end if
      call write3d(GridState,outfilename,'HORIZCORFUNC',variableType)
      
      !- 4.2 1D correlation function in txt format (for plotting purposes)
      iStart=iref
      iEnd=3*ni/4 ! About 10 000 km away from the center of the domain
      
      do varIndex = 1, nvar3d
        if ( nHorizWaveBand == 1 ) then
          outfilename = "./horizCorrel_"//trim(nomvar3d(varIndex,variableType))//".txt"
        else
          outfilename = "./horizCorrel_"//trim(nomvar3d(varIndex,variableType))//"_"//wbnum//".txt"
        end if
        open (unit=99,file=outfilename,action="write",status="new")
        
        if(vnl_varLevelFromVarName(nomvar3d(varIndex,variableType)).eq.'MM') then
          nLevEns = nLevEns_M
        else
          nLevEns = nLevEns_T
        end if
        nLevStart = varLevOffset(varIndex)+ 1 
        nLevEnd   = varLevOffset(varIndex)+ nLevEns
        
        do ji=iStart,iEnd
          do jk = nLevStart,nLevEnd
            if ( jk == nLevStart  ) then
              write(99,'(I7,2X,F7.1,2X,F6.4,$)')  ji-iStart, (ji-iStart)*gridSpacingInKm, GridState(ji,jref,jk)
            else if ( jk == nLevEnd ) then 
              write(99,'(2X,F6.4)')  GridState(ji,jref,jk) ! Saut de ligne
            else
              write(99,'(2X,F6.4,$)')  GridState(ji,jref,jk)
            end if
          end do
        end do
        close(unit=99)
      end do
      
      do varIndex = 1, nvar2d
        if ( nHorizWaveBand == 1 ) then
          outfilename = "./horizCorrel_"//trim(nomvar2d(varIndex,variableType))//".txt"
        else
          outfilename = "./horizCorrel_"//trim(nomvar2d(varIndex,variableType))//"_"//wbnum//".txt"
        end if
        open (unit=99,file=outfilename,action="write",status="new")
        do ji=iStart,iEnd
          write(99,'(I7,2X,F7.1,2X,F6.4)')  ji-iStart, (ji-iStart)*gridSpacingInKm, GridState(ji,jref,varLevOffset(nvar3d+1)+varIndex)
        end do
        close(unit=99)
      end do

    end if

  end subroutine horizCorrelFunction

  !--------------------------------------------------------------------------
  ! WRITE3D
  !--------------------------------------------------------------------------
  subroutine write3d(gridpoint3d,filename,etiket_in,variableType)
    !
    !:Purpose: Write the 3D stddev fields for all variables
    !
    implicit none

    ! Arguments:
    real(8),          intent(in) :: gridpoint3d(:,:,:)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: etiket_in
    integer,          intent(in) :: variableType

    ! Locals:
    real(8) :: dfact,zbuf(ni,nj)
    integer latIndex,lonIndex,levIndex,ierr,varIndex,nLevEns
    integer fstouv,fnom,fstfrm,fclos
    integer ip1,ip2,ip3,idatyp,idateo,ipak,nip1_l(max(nLevEns_M,nLevens_T))
    integer :: nulstats
    character(len=12) :: etiket

    etiket=trim(etiket_in)

    nulstats=0
    ierr =  fnom  (nulstats,filename,'RND',0)
    ierr =  fstouv(nulstats,'RND')

    ipak = -32
    idatyp = 5
    ip1 = 0
    ip2 = 0
    ip3 = nens
    idateo = 0

    ! do 3d variables
    do varIndex = 1, nvar3d

      if(vnl_varLevelFromVarName(nomvar3d(varIndex,variableType)).eq.'MM') then
        nLevEns = nLevEns_M
        nip1_l(1:nLevEns_M)=nip1_M(1:nLevEns_M)
      else
        nLevEns = nLevEns_T
        nip1_l(1:nLevEns_T)=nip1_T(1:nLevEns_T)
      end if
      dfact=1.0d0

      do levIndex = 1, nlevEns
        do latIndex = myLatBeg, myLatEnd
          do lonIndex = myLonBeg, myLonEnd
            zbuf(lonIndex,latIndex) = dfact*gridpoint3d(lonIndex,latIndex,varLevOffset(varIndex)+levIndex)
          end do
        end do
        ierr = utl_fstecr(zbuf,ipak,nulstats,idateo,0,0,ni,nj,1,nip1_l(levIndex),ip2,ip3,   &
                          'E',nomvar3d(varIndex,variableType),etiket,'G',0,0,0,0,idatyp,.true.)
      end do

    end do

    ! now do 2D variables
    do varIndex = 1, nvar2d
      !if(nomvar2d(varIndex,variableType).eq.'P0') then
      !  dfact=1.0d0/1.0d2
      !else
        dfact = 1.0d0
      !end if

      do latIndex = myLatBeg, myLatEnd
        do lonIndex = myLonBeg, myLonEnd
          zbuf(lonIndex,latIndex) = dfact*gridpoint3d(lonIndex,latIndex,varLevOffset(nvar3d+1)+varIndex)
        end do
      end do
      ierr = utl_fstecr(zbuf,ipak,nulstats,idateo,0,0,ni,nj,1,0,ip2,ip3,   &
                        'E',nomvar2d(varIndex,variableType),etiket,'G',0,0,0,0,idatyp,.true.)

    end do

    ierr =  fstfrm(nulstats)
    ierr =  fclos (nulstats)

    write(*,*) 'finished writing 3d array...'

  end subroutine write3d

  !--------------------------------------------------------------------------
  ! calcHorizScale
  !--------------------------------------------------------------------------
  subroutine calcHorizScale(rstddev,variableType,waveBandIndex_opt)
    !
    !:Purpose: Calculate the horizontal correlation length scale
    !
    implicit none

    ! Based on subroutine corrlength.ftn in the "old" var code

    ! Arguments:
    real(8),           intent(in) :: rstddev(nkgdimEns,0:ntrunc)
    integer,           intent(in) :: variableType
    integer, optional, intent(in) :: waveBandIndex_opt

    ! Locals:
    real(8) :: HorizScale(nkgdimEns)
    real(8), pointer :: PressureProfile(:)
    integer :: jk, jn, nLevEns, varIndex
    real(8) :: rjn, fact, temp, a, b
    character(len=128) :: outfilename
    character(len=2) :: wbnum

    if (mmpi_myid /= 0) return
    
    write(*,*)
    write(*,*) 'Computing horizontal correlation lengthscales'

    !
    !- 1.  Compute the horizontal correlation lengthscales
    !
    do jk = 1, nkgdimEns
       a = 0.d0
       b = 0.d0   
       do jn = 0, ntrunc
          rjn = dble(jn)
          fact = (2.0d0*rjn + 1.0d0)/2.0d0
          temp  = (rstddev(jk,jn)**2) * fact
          a = a + temp
          if (jn /= 0) then
             b = b - temp*rjn*(rjn+1.d0)
          end if
       end do
       if ( a > 0.d0 .and. b /= 0.d0 ) then
          HorizScale(jk) = ec_ra * sqrt(-2.0d0*a/b)
       else
          HorizScale(jk) = 0.d0
       end if
    end do

    !
    !- 2. Write the results
    !
    if ( nHorizWaveBand /= 1 ) then
       if (.not. present(waveBandIndex_opt)) then
          write(*,*) 'CalcHorizScale: No waveBandIndex was supplied!!!'
          call utl_abort('calbmatrix_glb')
       end if
       write(wbnum,'(I2.2)') waveBandIndex_opt
    end if

    do varIndex = 1, nvar3d
       write(*,*)
       write(*,*) nomvar3d(varIndex,variableType)

       if ( nHorizWaveBand == 1 ) then
          outfilename = "./horizScale_"//trim(nomvar3d(varIndex,variableType))//".txt"
       else
          outfilename = "./horizScale_"//trim(nomvar3d(varIndex,variableType))//"_"//wbnum//".txt"
       end if
       open (unit=99,file=outfilename,action="write",status="new")

       if(vnl_varLevelFromVarName(nomvar3d(varIndex,variableType)).eq.'MM') then
          nLevEns = nLevEns_M
          PressureProfile => pressureProfile_M
       else
          nLevEns = nLevEns_T
          PressureProfile => pressureProfile_T
       end if
       do jk=1,nlevEns
          write(* ,'(I3,2X,F6.1,2X,F6.1)')  jk, PressureProfile(jk)/100.d0, HorizScale(varLevOffset(varIndex)+jk)/1000.d0
          write(99,'(I3,2X,F6.1,2X,F6.1)')  jk, PressureProfile(jk)/100.d0, HorizScale(varLevOffset(varIndex)+jk)/1000.d0
       end do

       close(unit=99)
    end do

    do varIndex = 1, nvar2d
       write(*,*)
       write(*,*) nomvar2d(varIndex,variableType)
       
       if ( nHorizWaveBand == 1 ) then
          outfilename = "./horizScale_"//trim(nomvar2d(varIndex,variableType))//".txt"
       else
          outfilename = "./horizScale_"//trim(nomvar2d(varIndex,variableType))//"_"//wbnum//".txt"
       end if
       open (unit=99,file=outfilename,action="write",status="new")

       write(* ,'(I3,2X,F6.1,2X,F6.1)') 1, 1010.0, HorizScale(varLevOffset(nvar3d+1)+varIndex)/1000.d0
       write(99,'(I3,2X,F6.1,2X,F6.1)') 1, 1010.0, HorizScale(varLevOffset(nvar3d+1)+varIndex)/1000.d0

       close(unit=99)
    end do

  end subroutine CalcHorizScale

  !--------------------------------------------------------------------------
  ! WRITEPOWERSPEC
  !--------------------------------------------------------------------------
  subroutine writePowerSpec(powerSpec,variableType)
    !
    !:Purpose: Write the computed power spectrum to an ascii file
    !
    implicit none

    ! Arguments:
    real(8), intent(in) :: powerSpec(nkgdimEns,0:ntrunc)
    integer, intent(in) :: variableType

    ! Locals:
    integer :: jk, nLevEns, nLevStart, nLevEnd, varIndex, jn
    real(8) :: waveLength 
    character(len=128) :: outfilename

    if (mmpi_myid /= 0) return
    
    !- Write to txt files
    do varIndex = 1, nvar3d

       outfilename = "./PowerSpec_"//trim(nomvar3d(varIndex,variableType))//".txt"
       open (unit=99,file=outfilename,action="write",status="new")

       if(vnl_varLevelFromVarName(nomvar3d(varIndex,variableType)).eq.'MM') then
          nLevEns = nLevEns_M
       else
          nLevEns = nLevEns_T
       end if
       nLevStart = varLevOffset(varIndex)+ 1 
       nLevEnd   = varLevOffset(varIndex)+ nLevEns

       do jn = 0, ntrunc
          if ( jn /= 0) then
             waveLength=2*MPC_PI_R8*ec_ra/dble(jn)
          else
             waveLength=0.d0
          end if
          do jk = nLevStart,nLevEnd
             if ( jk == nLevStart  ) then
                write(99,'(I4,2X,F7.1,2X,e10.3,$)')  jn, waveLength/1000.d0, sngl(powerSpec(jk,jn))
             else if ( jk == nLevEnd ) then 
                write(99,'(2X,e10.3)')  sngl(powerSpec(jk,jn)) ! Saut de ligne
             else
                write(99,'(2X,e10.3,$)')  sngl(powerSpec(jk,jn))
             end if
          end do
       end do
       close(unit=99)
    end do

    do varIndex = 1, nvar2d

       outfilename = "./PowerSpec_"//trim(nomvar2d(varIndex,variableType))//".txt"
       open (unit=99,file=outfilename,action="write",status="new")
       do jn = 0, ntrunc
          if ( jn /= 0) then
             waveLength=2*MPC_PI_R8*ec_ra/dble(jn)
          else
             waveLength=0.d0
          end if
          write(99,'(I4,2X,F7.1,2X,e10.3)')  jn, waveLength/1000.d0, sngl(powerSpec(varLevOffset(nvar3d+1)+varIndex,jn))
       end do
       close(unit=99)
    end do

  end subroutine writePowerSpec

  !--------------------------------------------------------------------------
  ! calcLocalCorrelations (identical to the routine with the same name in calcstatslam)
  !--------------------------------------------------------------------------
  subroutine calcLocalCorrelations(ensPerts)
    !
    !:Purpose: Compute local horizontal correlations
    !
    implicit none

    ! Arguments:
    type(struct_ens), intent(in) :: ensPerts

    ! Locals:
    type(struct_gsv) :: statevector_locHorizCor
    type(struct_gsv) :: statevector_oneMember
    type(struct_gsv) :: statevector_oneMemberTiles
    real(8), pointer :: ptr3d_r8(:,:,:)
    real(8), pointer :: ptr3d_r8_oneMember(:,:,:)
    real(8) :: dnEns
    integer :: i, j, k, ens
    integer :: blocklength_x, blocklength_y
    integer :: iref_id, jref_id, iref, jref
    integer :: imin, imax, jmin, jmax
    character(len=4), pointer :: varNamesList(:)
    integer :: ierr, fclos, fnom, nulnam

    ! Namelist variables
    integer :: blockpadding
    integer :: nirefpoint
    integer :: njrefpoint

    NAMELIST /NAMHVCORREL_LOCAL/nirefpoint, njrefpoint, blockpadding

    !
    ! To compute the local horizontal correlation for some 'reference' grid point
    ! ... we assume that the ensemble grid point mean was removed and that
    !     the ensemble values were divided by the grid point std dev.
    !

    nirefpoint = 4 ! Number of reference grid point in x
    njrefpoint = 2 ! Number of reference grid point in y
    blockpadding = 4  ! Number of grid point padding between blocks (to set correlation to 0 between each block)

    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=NAMHVCORREL_LOCAL)
    if (ierr /= 0) call utl_abort('calcLocalCorrelations: Error reading namelist NAMHVCORREL_LOCAL')
    if (mmpi_myid == 0) write(*,nml=NAMHVCORREL_LOCAL)
    ierr = fclos(nulnam)

    blocklength_x = hco_ens%ni / nirefpoint ! Horizontal correlation will be compute blocklength x blocklength gridpoint
    ! around each reference point
    blocklength_y = hco_ens%nj / njrefpoint ! Horizontal correlation will be compute blocklength x blocklength gridpoint
    ! around each reference point

    nullify(varNamesList)
    call ens_varNamesList(varNamesList,ensPerts) 

    call gsv_allocate(statevector_locHorizCor, ens_getNumStep(ensPerts),                     &
                      ens_getHco(ensPerts), ens_getVco(ensPerts), varNames_opt=varNamesList, &
                      datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true.,                &
                      mpi_distribution_opt='VarsLevs', dataKind_opt=8 )

    call gsv_allocate(statevector_oneMemberTiles, ens_getNumStep(ensPerts),                  &
                      ens_getHco(ensPerts), ens_getVco(ensPerts), varNames_opt=varNamesList, &
                      datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true.,                &
                      mpi_distribution_opt='Tiles', dataKind_opt=8 )

    call gsv_allocate(statevector_oneMember, ens_getNumStep(ensPerts),                       &
                      ens_getHco(ensPerts), ens_getVco(ensPerts), varNames_opt=varNamesList, &
                      datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true.,                &
                      mpi_distribution_opt='VarsLevs', dataKind_opt=8 )

    call gsv_zero(statevector_locHorizCor)

    dnEns = 1.0d0/dble(nEns-1)

    call gsv_getField(statevector_locHorizCor,ptr3d_r8)

    do ens = 1, nEns
      call ens_copyMember(ensPerts, statevector_oneMemberTiles, ens)
      call gsv_transposeTilesToVarsLevs(statevector_oneMemberTiles, statevector_oneMember)
      call gsv_getField(statevector_oneMember,ptr3d_r8_oneMember)

      do k = statevector_locHorizCor%mykBeg, statevector_locHorizCor%mykEnd
        do jref_id = 1, njrefpoint
          do iref_id = 1, nirefpoint
            iref = (2*iref_id-1)*blocklength_x/2
            jref = (2*jref_id-1)*blocklength_y/2
            jmin = max(jref-(blocklength_y-blockpadding)/2,1)
            jmax = min(jref+(blocklength_y-blockpadding)/2,hco_ens%nj)
            imin = max(iref-(blocklength_x-blockpadding)/2,1)
            imax = min(iref+(blocklength_x-blockpadding)/2,hco_ens%ni)
            do j = jmin, jmax
              do i = imin, imax
                ptr3d_r8(i,j,k) = ptr3d_r8(i,j,k) + &
                     ptr3d_r8_oneMember(i,j,k)*ptr3d_r8_oneMember(iref,jref,k)
              end do
            end do
          end do
        end do
      end do

    end do

    call gsv_scale(statevector_locHorizCor,dnEns)

    write(*,*) 'finished computing the local horizontal correlations...'

    !
    !- 4.  Write to file
    !
    call gio_writeToFile(statevector_locHorizCor, './horizCorrelLocal.fst', 'HCORREL_LOC', &
                         typvar_opt = 'E', numBits_opt = 32)

    call gsv_deallocate(statevector_locHorizCor)
    call gsv_deallocate(statevector_oneMember)
    call gsv_deallocate(statevector_oneMemberTiles)

  end subroutine calcLocalCorrelations
  
  !--------------------------------------------------------------------------
  ! calcLocalVertCorrMatrix
  !--------------------------------------------------------------------------
  subroutine calcLocalVertCorrMatrix(ensPerts)
    !
    !:Purpose: Compute all vertical and between-variable local correlations
    !
    implicit none

    ! Arguments:
    type(struct_ens), intent(in) :: ensPerts

    ! Locals:
    type(struct_gsv) :: statevector_vertCorr
    type(struct_gsv) :: statevector_oneMember
    real(8), pointer :: ptr3d_r8(:,:,:)
    real(8), pointer :: ptr3d_r8_oneMember(:,:,:)
    real(8) :: dnEns
    integer :: lonIndex, latIndex, varLevIndex1, varLevIndex2, memberIndex
    integer :: levIndex1, dateStamp
    character(len=4), pointer :: varNamesList(:)
    character(len=12)  :: etiket
    character(len=4)   :: varName
    character(len=3)   :: levIndexStr
    character(len=256) :: ensFileName, outFileName

    ! Set the dateStamp using the first ensemble member
    call fln_ensfileName(ensFileName, ens_getPathName(ensPerts), memberIndex_opt=1)
    dateStamp = tim_getDatestampFromFile(ensFileName)

    ! Get list of variable names in the ensemble
    nullify(varNamesList)
    call ens_varNamesList(varNamesList,ensPerts) 

    call gsv_allocate(statevector_vertCorr, ens_getNumStep(ensPerts),                        &
                      ens_getHco(ensPerts), ens_getVco(ensPerts), varNames_opt=varNamesList, &
                      datestamp_opt=dateStamp, mpi_local_opt=.true.,                         &
                      mpi_distribution_opt='Tiles', dataKind_opt=8 )

    call gsv_allocate(statevector_oneMember, ens_getNumStep(ensPerts),                       &
                      ens_getHco(ensPerts), ens_getVco(ensPerts), varNames_opt=varNamesList, &
                      datestamp_opt=dateStamp, mpi_local_opt=.true.,                         &
                      mpi_distribution_opt='Tiles', dataKind_opt=8 )

    dnEns = 1.0d0/dble(nEns-1)

    call gsv_getField(statevector_vertCorr,ptr3d_r8)

    ! Loop over all vertical levels and variables
    varLev1: do varLevIndex1 = statevector_vertCorr%mykBeg, statevector_vertCorr%mykEnd

      varName = ens_getVarNameFromK(ensPerts,varLevIndex1)
      levIndex1 = ens_getLevFromK(ensPerts,varLevIndex1)

      ! Compute vertical correlations relative to varLev1
      call gsv_zero(statevector_vertCorr)
      member: do memberIndex = 1, nEns
        call ens_copyMember(ensPerts, statevector_oneMember, memberIndex)
        call gsv_getField(statevector_oneMember,ptr3d_r8_oneMember)
        varLev2: do varLevIndex2 = statevector_vertCorr%mykBeg, statevector_vertCorr%mykEnd
          do latIndex = statevector_vertCorr%myLatBeg, statevector_vertCorr%myLatEnd
            do lonIndex = statevector_vertCorr%myLonBeg, statevector_vertCorr%myLonEnd
              ptr3d_r8(lonIndex,latIndex,varLevIndex2) = &
                     ptr3d_r8(lonIndex,latIndex,varLevIndex2) + &
                     ptr3d_r8_oneMember(lonIndex,latIndex,varLevIndex1)* &
                     ptr3d_r8_oneMember(lonIndex,latIndex,varLevIndex2)
            end do
          end do
        end do varLev2
      end do member
      call gsv_scale(statevector_vertCorr,dnEns)

      ! Write to file the correlation matrix 'row' for this value of varLev1
      write(levIndexStr,'(i3.3)') levIndex1
      etiket = 'VCOR_' // trim(varName) // levIndexStr
      outFileName = './vertCorr_' // trim(varName) // levIndexStr // '.fst'
      call gio_writeToFile(statevector_vertCorr, &
                           trim(outFileName), etiket_in = etiket, &
                           typvar_opt = 'E', numBits_opt = 32)
      write(*,*) 'calcLocalVertCorrMatrix: finished variable/level =', varName, levIndex1
    end do varLev1

    write(*,*) 'calcLocalVertCorrMatrix: finished computing the local vertical correlation matrix...'

    call gsv_deallocate(statevector_vertCorr)
    call gsv_deallocate(statevector_oneMember)

  end subroutine calcLocalVertCorrMatrix

  !--------------------------------------------------------------------------
  ! calcVertModesSpec
  !--------------------------------------------------------------------------
  subroutine calcVertModesSpec(ensPerts, vModes)
    !
    ! :Purpose: Compute the amplitude of the ensemble perturbations after their
    !           projection onto the vertical modes. This is the vertical equivalent
    !           of power spectra in the horizontal.
    !
    implicit none

    ! Arguments:
    type(struct_ens), intent(inout)  :: ensPerts
    type(struct_vms), intent(in)     :: vModes

    ! Locals:
    type(struct_gsv)     :: gridStateVector_oneMember    
    real(8), allocatable :: vertModesState3d(:,:,:)
    real(8), pointer     :: gridState3d(:,:,:)
    real(8), allocatable :: powerSpec(:,:)
    real(8), allocatable :: latWeight(:) ! Weight given to grid point in the statistic computation
    real(8) :: sumWeight
    character(len=4), pointer :: varNamesList(:)
    character(len=128) :: outfilename
    integer :: numVar, nLev, varIndex, memberIndex
    integer :: nMode, modeIndex, latIndex, lonIndex

    !
    !- Setup
    !
    nullify(varNamesList)
    call ens_varNamesList(varNamesList,ensPerts)
    numVar = size(varNamesList)

    allocate(powerSpec(numVar,max(nLevEns_M,nLevEns_T)))
    powerSpec(:,:)=0.0d0
    
    call gsv_allocate(gridStateVector_oneMember, ens_getNumStep(ensPerts), &
                      ens_getHco(ensPerts), ens_getVco(ensPerts), varNames_opt=varNamesList, &
                      datestamp_opt=tim_getDatestamp(), mpi_local_opt=.true.,                         &
                      mpi_distribution_opt='Tiles', dataKind_opt=8 )

    allocate(latWeight(hco_ens%nj))
    do latIndex = 1, hco_ens%nj
      latWeight(latIndex) = cos(hco_ens%lat(latIndex))
      if (mmpi_myid == 0) then
        write(*,*) latIndex, hco_ens%lat(latIndex), cos(hco_ens%lat(latIndex))
      end if
    end do

    sumWeight = 0.d0
    do latIndex = myLatBeg, myLatEnd
      do lonIndex = myLonBeg, myLonEnd
        sumWeight = sumWeight + latWeight(latIndex)
      end do
    end do
    
    !
    !- Vertical modes decomposition and coefficient summation
    !    
    do memberIndex = 1, nEns
      write(*,*) 'calcVertModesSpec, member index = ', memberIndex
      call ens_copyMember(ensPerts, gridStateVector_oneMember, memberIndex)

      do varIndex = 1, numVar
        write(*,*) 'calcVertModesSpec, varName = ', varIndex, varNamesList(varIndex)
        if (vnl_varLevelFromVarName(trim(varNamesList(varIndex))).eq.'MM') then
          nLev  = nLevEns_M
          nMode = nLev
        else if (vnl_varLevelFromVarName(trim(varNamesList(varIndex))).eq.'TH') then
          nLev  = nLevEns_T
          nMode = nLev
        else
          if (memberIndex == 1) then
            write(*,*) '   Skipping this variable not on momentum or thermodynamic levels'
          end if
          cycle  
        end if

        nullify(gridState3d)
        call gsv_getField(gridStateVector_oneMember,gridState3d,varName_opt=varNamesList(varIndex))

        if (allocated(vertModesState3d)) deallocate(vertModesState3d)
        allocate(vertModesState3d(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nMode))
        
        call vms_transform(vModes, vertModesState3d, gridState3d, &
                          'GridPointToVertModes', myLonBeg, myLonEnd, &
                           myLatBeg, myLatEnd, nLev, varNamesList(varIndex))

        !$OMP PARALLEL DO PRIVATE (modeIndex,latIndex,lonIndex)
        do modeIndex = 1, nMode
          do latIndex = myLatBeg, myLatEnd
            do lonIndex = myLonBeg, myLonEnd
              powerSpec(varIndex,modeIndex) = powerSpec(varIndex,modeIndex) &
                                            + vertModesState3d(lonIndex,latIndex,modeIndex)**2 &
                                            * latWeight(latIndex)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
        
      end do

    end do
    
    deallocate(vertModesState3d)
    call gsv_deallocate(gridStateVector_oneMember)

    !
    !- Communicate between all tasks
    !
    call mmpi_allreduce_sumR8_2d(powerSpec, "GRID")
    call mmpi_allreduce_sumreal8scalar(sumWeight, "GRID")
    
    !
    !- Apply the appropriate scaling
    !
    powerSpec(:,:) = powerSpec(:,:)/(dble(nEns-1)*sumWeight)

    !
    !- Write to file
    !
    if (mmpi_myid == 0) then
      do varIndex = 1, numVar

        if (vnl_varLevelFromVarName(trim(varNamesList(varIndex))) /= 'MM' .and. &
            vnl_varLevelFromVarName(trim(varNamesList(varIndex))) /= 'TH') cycle
        
        outfilename = "./vertModesSpec_"//trim(varNamesList(varIndex))//".txt"
        write(*,*) 'calcVertModesSpec, writing ', trim(outfilename)
        open (unit=99,file=outfilename,action="write",status="new")

        if (vnl_varLevelFromVarName(trim(varNamesList(varIndex))) == 'MM') then
          nMode = nLevEns_M
        else
          nMode = nLevEns_T
        end if
        
        do modeIndex = 1, nMode
          write(99,'(I4,2X,E10.4)') modeIndex, powerSpec(varIndex,modeIndex)
        end do
        close(unit=99)
      end do
    end if
    
  end subroutine calcVertModesSpec
  
end module calcStatsGlb_mod
