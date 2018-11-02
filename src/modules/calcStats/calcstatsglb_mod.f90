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
!! MODULE CalcStatsGlb (prefix="csg" category='1')
!!
!! *Purpose*: Compute homogeneous and isotropic background error covariances 
!!            from forecast error estimate in model variable space (global version).
!!
!--------------------------------------------------------------------------
module calcstatsglb_mod
  use gridStateVector_mod
  use globalSpectralTransform_mod
  use mathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use varNameList_mod
  use earthConstants_mod, only: RA
  use utilities_mod
  use spectralFilter_mod
  use menetrierDiag_mod
  implicit none
  save
  private

  ! Public Subroutines
  public :: csg_setup, csg_computeStats, csg_computeStatsLatBands, csg_stddev
  public :: csg_toolbox, csg_powerspec

  type(struct_hco), pointer :: hco_ens ! Ensemble horizontal grid parameters

  integer :: nens,ntrunc,ni,nj,nLevEns_M,nLevEns_T,nLevPtoT,nkgdimEns,varLevOffset(6),nla
  character(len=256), allocatable :: cflensin(:)
  integer :: gstID_nkgdimEns, gstID_nLevEns_M, gstID_nLevEns_T_P1
  integer, allocatable :: nip1_M(:),nip1_T(:)
  real(8), pointer :: pressureProfile_M(:), pressureProfile_T(:)
  integer,external    :: get_max_rss

  integer,parameter  :: nvar3d=4,nvar2d=1
  character*4 :: nomvar3d(nvar3d,3),nomvar2d(nvar2d,3)
  integer, parameter :: modelSpace   = 1
  integer, parameter :: cvSpace      = 2
  integer, parameter :: cvUnbalSpace = 3
!  integer :: variableType

  real(8) :: gridSpacingInKm

  ! For wave band decomposition
  integer, parameter  :: maxNumLocalLength = 20
  integer             :: waveBandPeaks(maxNumLocalLength)
  integer             :: nWaveBand

  logical :: initialized = .false.

  contains

!--------------------------------------------------------------------------
! CSG_SETUP
!--------------------------------------------------------------------------
  subroutine csg_setup( nens_in, cflens_in, hco_in, vco_in)
    implicit none

    integer, intent(in)            :: nens_in

    character(len=*), intent(in)   :: cflens_in(:)

    type(struct_vco), pointer, intent(in)   :: vco_in
    type(struct_hco), pointer, intent(in)   :: hco_in

    integer :: nulnam,ierr,status, waveBandIndex
    integer :: fclos,fnom,fstouv,fstfrm

    real(8) :: zps

    NAMELIST /NAMCALCSTATS_GLB/ntrunc,waveBandPeaks

    write(*,*)
    write(*,*) 'csg_setup: Starting...'

    nens=nens_in
    allocate(cflensin(nens))
    cflensin(:)=cflens_in(1:nens)
    call mpc_printConstants(6)

    ! parameters from namelist (date in filename should come directly from sequencer?)
    ntrunc=108
    waveBandPeaks(:) = -1.0d0

    nulnam=0
    ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=NAMCALCSTATS_GLB)
    write(*,nml=NAMCALCSTATS_GLB)
    ierr=fclos(nulnam)

    !- Setup horizontal grid
    hco_ens => hco_in
    ni=hco_in%ni
    nj=hco_in%nj
    gridSpacingInKm = ra * hco_in%dlon / 1000.d0
    
    write(*,*) 'Grid Spacing in Km = ', gridSpacingInKm

    !- Setup vertical levels
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
    !  NOTE: code only run with mpi 1x1 until now!!!
    gstID_nkgdimEns = gst_setup(ni,nj,ntrunc,nkgdimEns)
    gstID_nLevEns_M = gst_setup(ni,nj,ntrunc,nLevEns_M)
    gstID_nLevEns_T_P1 = gst_setup(ni,nj,ntrunc,nLevEns_T+1)

    !- Setup ip1s
    allocate(nip1_M(nLevEns_M))
    nip1_M(:)=vco_in%ip1_M(:)
    allocate(nip1_T(nLevEns_T))
    nip1_T(:)=vco_in%ip1_T(:)

    !- Estimate the pressure profile for each vertical grid    
    zps = 101000.D0
    status = vgd_levels( vco_in%vgrid, ip1_list=vco_in%ip1_M, levels=pressureProfile_M, &
                         sfc_field=zps, in_log=.false.)

    status = vgd_levels( vco_in%vgrid, ip1_list=vco_in%ip1_T, levels=pressureProfile_T, &
                         sfc_field=zps, in_log=.false.)

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

    !
    !- Wave band decomposition option (available in )
    !
    nWaveBand = count(waveBandPeaks .ge. 0)
    if ( nWaveBand < 1 ) then
       nWaveBand = 1
    else if (nWaveBand == 1) then
       write(*,*) 'You have specified only ONE waveBandPeaks'
       call utl_abort('calbmatrix_glb')
    else
       write(*,*)
       write(*,*) 'WaveBand decomposition is ACTIVATED'
    end if
    
    ! Make sure that the wavenumbers are in the correct (decreasing) order
    do waveBandIndex = 1, nWaveBand-1
       if ( waveBandPeaks(waveBandIndex)-waveBandPeaks(waveBandIndex+1) <= 0 ) then
          write(*,*) 'csg_setup: waveBandPeaks are not in decreasing wavenumber order'
          call utl_abort('calbmatrix_glb')
       end if
    end do

    ! Make sure the truncation is compatible with the waveBandPeaks
    if ( ntrunc < nj-1 .and. nWaveBand > 1 ) then
       write(*,*) 'csg_setup: The truncation is not compatible with wave band decomposition'
       write(*,*) '                 ntrunc should = ', nj-1
       call utl_abort('calbmatrix_glb')
    end if

    !
    !- Ending
    !
    initialized = .true.

    write(*,*)
    write(*,*) 'csg_setup: Done!'

  end subroutine csg_setup

!--------------------------------------------------------------------------
! CSG_COMPUTESTATS
!--------------------------------------------------------------------------
  subroutine csg_computeStats
    implicit none
    integer :: ierr
    real*4,pointer  :: ensPerturbations(:,:,:,:)
    real*4,pointer  :: ensBalPerturbations(:,:,:,:)
    real*8,allocatable :: stddev3d(:,:,:),stddev3dBal(:,:,:),stddev3dUnbal(:,:,:)
    real*8,allocatable :: stddevZonAvg(:,:),stddevZonAvgBal(:,:),stddevZonAvgUnbal(:,:)
    real*8,allocatable :: PtoT(:,:,:),theta1(:,:),theta2(:,:)
    real*8,allocatable :: corns(:,:,:),rstddev(:,:)

    integer :: variableType 

    allocate(ensPerturbations(ni,nj,nkgdimEns,nens),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 1')
    allocate(ensBalPerturbations(ni,nj,nLevEns_T+1,nens),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 1')
    allocate(theta1(nlevEns_M,nj),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 1')
    allocate(stddev3d(ni,nj,nkgdimEns),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 1')
    allocate(stddevZonAvg(nkgdimEns,nj),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 1')
    allocate(PtoT(nlevEns_T+1,nlevEns_M,nj),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 1')
    allocate(theta2(nlevEns_M,nj),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 1')
    allocate(stddev3dBal(ni,nj,nLevEns_T+1),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 1')
    allocate(stddev3dUnbal(ni,nj,nkgdimEns),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 1')
    allocate(stddevZonAvgBal(nLevEns_T+1,nj),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 1')
    allocate(stddevZonAvgUnbal(nkgdimEns,nj),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 1')
    allocate(corns(nkgdimEns,nkgdimEns,0:ntrunc),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 1')
    allocate(rstddev(nkgdimEns,0:ntrunc),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 1')

    call readEnsemble(ensPerturbations)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    call removeMean(ensPerturbations)

    call uv_to_psichi(ensPerturbations)

    call calcStddev3d(ensPerturbations,stddev3d,nkgdimens)

    call calcZonAvg(stddevZonAvg,stddev3d,nkgdimens)

    call calcTheta(ensPerturbations,theta1) ! theta1 is put in glbcov and used for analysis!
    write(301,*) theta1

    call removeBalancedChi(ensPerturbations,theta1)

    call normalize3d(ensPerturbations,stddev3d)

    call calcPtoT(ensPerturbations,PtoT)
    write(303,*) PTOT(:,:,1)
    call flush(303)

!    call calcTheta(ensPerturbations,theta2) ! theta2 is used previously for computing unbalanced Chi!
!    write(302,*) theta2

    call removeBalancedT_Ps(ensPerturbations,ensBalPerturbations,PtoT)

!    call removeBalancedChi(ensPerturbations,theta2)

    call multiply3d(ensPerturbations,stddev3d,nkgdimens)

    call multiply3d(ensBalPerturbations(:,:,1:nLevEns_T,:),   &
                    stddev3d(:,:,(2*nLevEns_M+1):(2*nLevEns_M+nLevEns_T)),nLevEns_T)

    call multiply3d(ensBalPerturbations(:,:,(nLevEns_T+1):(nLevEns_T+1),:),  &
                    stddev3d(:,:,(2*nLevEns_M+2*nLevEns_T+1):(2*nLevEns_M+2*nLevEns_T+1)),1)

    call spectralFilter(ensPerturbations,nkgdimens)

    call spectralFilter(ensBalPerturbations,nLevEns_T+1)

    call calcStddev3d(ensPerturbations,stddev3dUnbal,nkgdimens)

    call calcStddev3d(ensBalPerturbations,stddev3dBal,nLevEns_T+1)

    call calcZonAvg(stddevZonAvgUnbal,stddev3dUnbal,nkgdimens)

    call calcZonAvg(stddevZonAvgBal,stddev3dBal,nLevEns_T+1)

    call normalize3d(ensPerturbations,stddev3dUnbal)

    call removeGlobalMean(ensPerturbations)

    call calcCorrelations(ensPerturbations,corns,rstddev)

    variableType = cvUnbalSpace
    call writeStats(corns,rstddev,variableType,ptot,theta1)

    call writeStddev(stddevZonAvg,stddevZonAvgUnbal,stddev3d,stddev3dUnbal)

    call writeStddevBal(stddevZonAvgBal,stddev3dBal)

    call writeSpStats(ptot,theta1)

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

  end subroutine csg_computeStats

!--------------------------------------------------------------------------
! CSG_COMPUTESTATSLATBANDS
!--------------------------------------------------------------------------
  subroutine csg_computeStatsLatBands
    integer :: ierr, variableType, jlat, jlatband, lat1, lat2, lat3
    real*4,pointer     :: ensPerturbations(:,:,:,:)
    real*8,allocatable :: stddev3d(:,:,:)
    real*8,allocatable :: stddevZonAvg(:,:),stddevZonAvgBal(:,:)
    real*8,allocatable :: corns(:,:,:),rstddev(:,:)
    real*8 :: latMask(nj)

    allocate(ensPerturbations(ni,nj,nkgdimEns,nens),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 1')
    allocate(stddev3d(ni,nj,nkgdimEns),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 2')
    allocate(stddevZonAvg(nkgdimEns,nj),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 3')
    allocate(corns(nkgdimEns,nkgdimEns,0:ntrunc),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 4')
    allocate(rstddev(nkgdimEns,0:ntrunc),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 5')

    call readEnsemble(ensPerturbations)
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    call removeMean(ensPerturbations)

    call uv_to_psichi(ensPerturbations)

    call calcStddev3d(ensPerturbations,stddev3d,nkgdimens)

    call calcZonAvg(stddevZonAvg,stddev3d,nkgdimens)

    call normalize3d(ensPerturbations,stddev3d)

    call removeGlobalMean(ensPerturbations)

    do jlatband = 1, 3
      write(*,*) 'calcb_glb2_computeStats: selected LATBAND = ',jlatband
      lat1=nj/4
      lat2=nj/2
      lat3=3*nj/4
      write(*,*) 'lat1,2,3=',lat1,lat2,lat3
      if(jlatband==1) then
        ! Southern extratropics
        latMask(1:lat1) = 1.0d0
        do jlat = lat1, lat2
          !latMask(jlat) = sqrt(dble((lat2-jlat)*4)/dble(nj))
          latMask(jlat) = sqrt(0.5d0*(1.0d0+cos(dble((jlat-lat1)*4)*MPC_PI_R8/dble(nj))))
        enddo
        latMask(lat2:nj) = 0.0d0
      elseif(jlatband==2) then
        ! Tropics
        !latMask(1:lat1) = 0.0d0
        !do jlat = lat1, lat2
        !  latMask(jlat) = sqrt(dble((jlat-lat1)*4)/dble(nj))
        !enddo
        !do jlat = lat2,lat3
        !  latMask(jlat) = sqrt(dble((lat3-jlat)*4)/dble(nj))
        !enddo
        !latMask(lat3:nj) = 0.0d0

        ! NOTE: use much broader band for tropics to avoid shortening horizontal correlations
        ! ok, since the masks do not have to sum to one for calculation, but they do
        ! when used in bmatrixhi_mod
        ! Tropics
        do jlat = 1, lat1
          !latMask(jlat) = sqrt(dble((jlat-1)*4)/dble(nj))
          latMask(jlat) = sqrt(0.5d0*(1.0d0+cos(dble((lat1-jlat)*4)*MPC_PI_R8/dble(nj))))
        enddo
        latMask(lat1:lat3) = 1.0d0
        do jlat = lat3,nj
          !latMask(jlat) = sqrt(dble((nj-jlat)*4)/dble(nj))
          latMask(jlat) = sqrt(0.5d0*(1.0d0+cos(dble((jlat-lat3)*4)*MPC_PI_R8/dble(nj))))
        enddo
      elseif(jlatband==3) then
        ! Northern extratropics
        latMask(1:lat2) = 0.0d0
        do jlat = lat2, lat3
          !latMask(jlat) = sqrt(dble((jlat-lat2)*4)/dble(nj))
          latMask(jlat) = sqrt(0.5d0*(1.0d0+cos(dble((lat3-jlat)*4)*MPC_PI_R8/dble(nj))))
        enddo
        latMask(lat3:nj) = 1.0d0
      endif
      write(*,*) 'latMask = ',latMask(:)
      call calcCorrelations(ensPerturbations,corns,rstddev,latMask_opt=latMask)

      variableType = cvSpace
      call writeStats(corns,rstddev,variableType,latBand_opt=jlatBand)
    enddo

    call writeStddev2(stddevZonAvg,stddev3d)

    write(200,*) stddevZonAvg(1:nlevEns_M,:)
    write(201,*) stddevZonAvg((1+1*nlevEns_M):(2*nlevEns_M),:)
    write(202,*) stddevZonAvg((1+2*nlevEns_M):(3*nlevEns_T),:)
    write(203,*) stddevZonAvg((1+2*nlevEns_M+1*nlevEns_T):(2*nlevEns_M+2*nlevEns_T),:)
    write(204,*) stddevZonAvg((1+2*nlevEns_M+2*nlevEns_T),:)/1.0d2

  end subroutine csg_computeStatsLatBands

!--------------------------------------------------------------------------
! CSG_TOOLBOOX
!--------------------------------------------------------------------------
  subroutine csg_toolbox
    implicit none

    ! NOTE: The diagnostic computed here are in model variable space 
    !       (no variable transform!!!)

    integer :: waveBandIndex
    integer :: nulnam, ierr, fclos, fnom

    real(4),pointer     :: ensPerturbations(:,:,:,:)
    real(8),allocatable :: stddev3d(:,:,:)
    real(8),allocatable :: corns(:,:,:), rstddev(:,:)

    integer :: variableType

    character(len=60) :: tool

    NAMELIST /NAMTOOLBOX/tool

    write(*,*)
    write(*,*) 'csg_toolbox'
    write(*,*)

    variableType = modelSpace

    !
    !- Tool selection
    !
    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=NAMTOOLBOX)
    write(*,nml=NAMTOOLBOX)
    ierr = fclos(nulnam)

    select case(trim(tool))
    case ('HVCORREL_HI')
       write(*,*)
       write(*,*) 'Computing Homogeneous and Isotropic Correlation'
    case ('HVCORREL_LOCAL')
       write(*,*)
       write(*,*) 'Computing Local Correlation'
    case ('LOCALIZATIONRADII')
       write(*,*)
       write(*,*) 'Estimating the optimal covariance localization radii'
       call bmd_setup( hco_ens, nens, nLevEns_M, nLevEns_T,              & ! IN
                       nkgdimEns, pressureProfile_M, pressureProfile_T,  & ! IN
                       nvar3d, nvar2d, varLevOffset, nomvar3d, nomvar2d, & ! IN
                       nWaveBand)                                          ! IN
    case default
       write(*,*)
       write(*,*) 'Unknown TOOL in csg_toolbox : ', trim(tool)
       call utl_abort('calbmatrix_glb')
    end select

    !
    !- Horizontal and vertical correlation diagnostics
    !
    allocate(ensPerturbations(ni,nj,nkgdimEns,nens))
    allocate(stddev3d(ni,nj,nkgdimEns))
    if (trim(tool) == 'HVCORREL_HI') then
       allocate(corns(nkgdimEns,nkgdimEns,0:ntrunc))
       allocate(rstddev(nkgdimEns,0:ntrunc))
    end if

    do waveBandIndex = 1, nWaveBand

       if ( nWaveBand /= 1 ) then
          write(*,*)
          write(*,*) ' ********* Processing WaveBand #',waveBandIndex
          write(*,*)
       end if

       call readEnsemble(ensPerturbations) ! OUT
       write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

       call removeMean(ensPerturbations) ! INOUT
       
       if (trim(tool) /= 'HVCORREL_HI') then
          call removeGlobalMean(ensPerturbations) ! INOUT
       end if

       call spectralFilter(ensPerturbations,nkgdimens,waveBandIndex_opt=waveBandIndex) ! INOUT, IN, IN

       call calcStddev3d(ensPerturbations,stddev3d,nkgdimens) ! IN, OUT, IN

       if (trim(tool) == 'HVCORREL_HI') then
          call normalize3d(ensPerturbations,stddev3d) ! INOUT, IN

          call removeGlobalMean(ensPerturbations) ! INOUT

          call calcCorrelations(ensPerturbations, & ! IN
                                corns,            & ! OUT (vertical correlation in spectral space)
                                rstddev)            ! OUT ( sqrt(normalized power spectrum) )

          call writeStats(corns,rstddev,variableType,waveBandIndex_opt=waveBandIndex) ! IN

          call calcHorizScale(rstddev,variableType,waveBandIndex_opt=waveBandIndex) ! IN

          call horizCorrelFunction(rstddev,variableType,waveBandIndex_opt=waveBandIndex) ! IN
       else if (trim(tool) == 'HVCORREL_LOCAL') then
          call normalize3d(ensPerturbations,stddev3d) ! INOUT, IN
          call calcLocalCorrelations(ensPerturbations, variableType, waveBandIndex_opt=waveBandIndex) ! IN
       else
          call bmd_localizationRadii(ensPerturbations, stddev3d, variableType, waveBandIndex_opt=waveBandIndex) ! IN
       end if

    end do

    !
    !- Write the estimated pressure profiles
    !
    call writePressureProfiles

  end subroutine csg_toolbox

!--------------------------------------------------------------------------
! CSG_POWERSPEC
!--------------------------------------------------------------------------
  subroutine csg_powerspec
    implicit none

    ! NOTE: The diagnostic computed here are in model variable space 
    !       (no variable transform!!!)

    integer :: variableType

    integer :: ierr,  waveBandIndex

    real*4,pointer     :: ensPerturbations(:,:,:,:)
    real*8,allocatable :: powerSpec(:,:)

    write(*,*)
    write(*,*) 'csg_powerspec'
    write(*,*)
    
    variableType = modelSpace

    allocate(ensPerturbations(ni,nj,nkgdimEns,nens))
    allocate(powerspec(nkgdimEns,0:ntrunc))

    call readEnsemble(ensPerturbations) ! OUT, IN
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    call removeMean(ensPerturbations) ! INOUT

    call calcPowerSpec(ensPerturbations, & ! IN
                       powerSpec)          ! OUT

    call writePowerSpec(powerSpec, variableType) ! IN

  end subroutine csg_powerspec

!--------------------------------------------------------------------------
! CSG_STDDEV
!--------------------------------------------------------------------------
  subroutine csg_stddev
    implicit none
    integer :: ierr
    real*4,pointer  :: ensPerturbations(:,:,:,:)
    real*8,allocatable :: stddev3d(:,:,:),stddevZonAvg(:,:)

    allocate(ensPerturbations(ni,nj,nkgdimEns,nens),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 1')
    allocate(stddev3d(ni,nj,nkgdimEns),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 1')
    allocate(stddevZonAvg(nkgdimEns,nj),stat=ierr)
    if(ierr.ne.0) call utl_abort('Problem allocating memory 1')

    call readEnsemble(ensPerturbations)

    call removeMean(ensPerturbations)

    call uv_to_psichi(ensPerturbations)

    call calcStddev3d(ensPerturbations,stddev3d,nkgdimens)

    call calcZonAvg(stddevZonAvg,stddev3d,nkgdimens)

    call writeStddev2(stddevZonAvg,stddev3d)

    write(200,*) stddevZonAvg(1:nlevEns_M,:)
    write(201,*) stddevZonAvg((1+1*nlevEns_M):(2*nlevEns_M),:)
    write(202,*) stddevZonAvg((1+2*nlevEns_M):(3*nlevEns_T),:)
    write(203,*) stddevZonAvg((1+2*nlevEns_M+1*nlevEns_T):(2*nlevEns_M+2*nlevEns_T),:)
    write(204,*) stddevZonAvg((1+2*nlevEns_M+2*nlevEns_T),:)/1.0d2

  end subroutine csg_stddev

!--------------------------------------------------------------------------
! WRITESPSTATS
!--------------------------------------------------------------------------
  subroutine writeSpStats(ptot,theta)
    implicit none
    real*8 :: PtoT(:,:,:),theta(:,:)
    integer jn,ierr,ipak,jlat,jlev1,jlev2,nlev
    integer fstouv,fnom,fstfrm,fclos
    integer ip1,ip2,ip3,kni,knj,idatyp,idateo
    integer :: nulstats
    real*8 :: bufz(nLevEns_M),bufyz(nj,nLevEns_M),zsp(0:ntrunc,nLevEns_M)
    real*8 :: bufptot(nj,(nLevEns_T+1)*nLevEns_M),spptot(0:ntrunc,(nLevEns_T+1)*nLevEns_M)
    real*8 :: zspptot(nLevEns_T+1,nLevEns_M)

    nulstats=0
    ierr =  fnom  (nulstats,'./stats_sp.fst','RND',0)
    ierr =  fstouv(nulstats,'RND')

    ipak = -32
    idatyp = 5
    ip1 = 0
    ip3 = nens
    idateo = 0

    ! write out SP_THETA

    do jlat = 1, nj
      do jlev1 = 1, nLevEns_M
        bufyz(jlat,jlev1) = theta(jlev1,jlat)
      enddo
    enddo

    call gst_zlegdir(gstID_nkgdimEns,bufyz,zsp,nLevEns_M)

    do jn = 0, ntrunc
      do jlev1=1, nLevEns_M
        bufz(jlev1) = zsp(jn,jlev1)
      enddo

      ierr = utl_fstecr(bufz,ipak,nulstats,idateo,0,0,nlevEns_M,1,1,   &
                        ip1,jn,ip3,'X','ZZ','SP_THETA','X',0,0,0,0,idatyp,.true.)

    enddo

    ! write out SP_PTOT

    do jlat = 1, nj
      do jlev1 = 1, (nLevEns_T+1)
        do jlev2 = 1, nLevEns_M
          bufptot(jlat,(jlev2-1)*(nLevEns_T+1)+jlev1) = PtoT(jlev1,jlev2,jlat)
        enddo
      enddo
    enddo

    nlev=(nLevEns_T+1)*nLevEns_M
    call gst_zlegdir(gstID_nkgdimEns,bufptot,spptot,nLev)

    do jn = 0, ntrunc
      do jlev1 = 1, (nLevEns_T+1)
        do jlev2 = 1, nLevEns_M
          zspptot(jlev1,jlev2) = spptot(jn,(jlev2-1)*(nLevEns_T+1)+jlev1)
        enddo
      enddo

      kni=nLevEns_T+1
      knj=nLevEns_M
      ierr = utl_fstecr(zspptot,ipak,nulstats,idateo,0,0,kni,knj,1,  &
                        ip1,jn,ip3,'X','ZZ','SP_PTOT ','X',0,0,0,0,idatyp,.true.)
    enddo




    ierr =  fstfrm(nulstats)
    ierr =  fclos (nulstats)

    write(*,*) 'finished writing statistics...'
    call flush(6)

  end subroutine writeSpStats

!--------------------------------------------------------------------------
! REMOVEBALANCEDCHI
!--------------------------------------------------------------------------
  subroutine removeBalancedChi(ensPerturbations,theta)
    implicit none
    real*4,pointer :: ensPerturbations(:,:,:,:)
    real*8 :: theta(:,:)
    real*4,pointer :: psi_ptr(:,:,:),chi_ptr(:,:,:)
    integer :: jens,jlat,jlev,jlon

    do jens = 1,nens
      psi_ptr => ensPerturbations(:,:,1:nlevEns_M,jens)
      chi_ptr => ensPerturbations(:,:,(nlevEns_M+1):(2*nlevEns_M),jens)

      do jlev = 1, nLevEns_M
        do jlat = 1, nj
          do jlon = 1, ni
            chi_ptr(jlon,jlat,jlev) = chi_ptr(jlon,jlat,jlev) + tan(theta(jlev,jlat))*psi_ptr(jlon,jlat,jlev)
          enddo
        enddo
      enddo

    enddo

    write(*,*) 'finished removing balanced chi...'
    call flush(6)

  end subroutine removeBalancedChi

!--------------------------------------------------------------------------
! REMOVEBALANCEDT_PS
!--------------------------------------------------------------------------
  subroutine removeBalancedT_Ps(ensPerturbations,ensBalPerturbations,PtoT)
    implicit none
    real*4,pointer :: ensPerturbations(:,:,:,:)
    real*4,pointer :: ensBalPerturbations(:,:,:,:)
    real*8 :: PtoT(:,:,:)

    real*4,pointer :: tt_ptr(:,:,:),ps_ptr(:,:,:),ttb_ptr(:,:,:),psb_ptr(:,:,:)
    real*8  :: spectralState(nla,2,nLevEns_M),spBalancedP(nla,2,nlevEns_M),balancedP(ni,nj,nlevEns_M),psi(ni,nj,nLevEns_M)
    integer :: jens,jlat,jlon,jk1,jk2

    do jens=1,nens

      psi(:,:,:)=ensPerturbations(:,:,1:nlevEns_M,jens)
      call gst_setID(gstID_nLevEns_M)
      call gst_reespe(spectralState,psi)
      call calcBalancedP(spectralState,spBalancedP)
      call gst_speree(spBalancedP,balancedP)

      tt_ptr => ensPerturbations(:,:,(1+2*nLevEns_M):(2*nLevEns_M+1*nLevEns_T),jens)
      ps_ptr => ensPerturbations(:,:,(1+2*nLevEns_M+2*nLevEns_T):(1+2*nLevEns_M+2*nLevEns_T),jens)
      ttb_ptr => ensBalPerturbations(:,:,1:nLevEns_T,jens)
      psb_ptr => ensBalPerturbations(:,:,(1+nLevEns_T):(1+nLevEns_T),jens)

      ttb_ptr(:,:,:)=0.0d0
      psb_ptr(:,:,:)=0.0d0

!$OMP PARALLEL DO PRIVATE (jlat,jlon,jk1,jk2)
      do jlat = 1, nj
        do jlon = 1, ni
          do jk1 = 1, nLevEns_T
            do jk2 = 1, nlevptot
              ttb_ptr(jlon,jlat,jk1) = ttb_ptr(jlon,jlat,jk1) + PtoT(jk1,jk2,jlat)*balancedP(jlon,jlat,jk2)
            enddo
            tt_ptr(jlon,jlat,jk1) = tt_ptr(jlon,jlat,jk1) - ttb_ptr(jlon,jlat,jk1)
          enddo
          do jk2 = 1, nlevptot
            psb_ptr(jlon,jlat,1) = psb_ptr(jlon,jlat,1) + PtoT(nLevEns_T+1,jk2,jlat)*balancedP(jlon,jlat,jk2)
          enddo
          ps_ptr(jlon,jlat,1) = ps_ptr(jlon,jlat,1) - psb_ptr(jlon,jlat,1)
        enddo
      enddo
!$OMP END PARALLEL DO

    enddo

    write(*,*) 'finished removing balanced T and Ps...'
    call flush(6)

  end subroutine removeBalancedT_Ps

!--------------------------------------------------------------------------
! CALCCORRELATIONS
!--------------------------------------------------------------------------
  subroutine calcCorrelations(ensPerturbations,corns,rstddev,latMask_opt)
    implicit none
    real*4,pointer :: ensPerturbations(:,:,:,:)
    real*8 :: corns(nkgdimEns,nkgdimEns,0:ntrunc),rstddev(nkgdimEns,0:ntrunc)
    real*8,  optional :: latMask_opt(:)

    real*8  :: spectralState(nla,2,nkgdimEns),gridState(ni,nj,nkgdimEns)
    real*8  :: dfact,dfact2,dsummed
    integer :: jens,ila,jn,jm,jk1,jk2,jlat

    call tmg_start(3,'CALCCORRELATIONS')

    corns(:,:,:)=0.0d0
    do jens=1,nens

      write(*,*) 'calcCorrelations: processing member ',jens
      call flush(6)

      gridState(:,:,:)=ensPerturbations(:,:,:,jens)
      if(present(latMask_opt)) then
        do jlat = 1, nj
          gridState(:,jlat,:) = latMask_opt(jlat)*gridState(:,jlat,:)
        enddo
      endif
      call gst_setID(gstID_nkgdimEns)
      call gst_reespe(spectralState,gridState)

!$OMP PARALLEL DO PRIVATE (jn,jm,dfact,ila,jk1,jk2)
      do jn = 0, ntrunc
        do jm = 0, jn
          dfact = 2.0d0
          if (jm.eq.0) dfact = 1.0d0
          ila = gst_getNind(jm) + jn - jm
          do jk1 = 1, nkgdimEns
            do jk2 = 1, nkgdimEns
              corns(jk1,jk2,jn) = corns(jk1,jk2,jn) +     &
                 dfact*( spectralState(ila,1,jk1)*spectralState(ila,1,jk2) +   &
                         spectralState(ila,2,jk1)*spectralState(ila,2,jk2)  )
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

    enddo

!$OMP PARALLEL DO PRIVATE (jn,jk1)
    do jn = 0, ntrunc
      do jk1 = 1, nkgdimEns
        if(abs(corns(jk1,jk1,jn)).gt.0.0d0) then
          rstddev(jk1,jn) = dsqrt(abs(corns(jk1,jk1,jn)))
        else
          rstddev(jk1,jn) = 0.0d0
        endif
      enddo
    enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE (jn,jk1,jk2)
    do jn = 0, ntrunc
      do jk1 = 1, nkgdimEns
        do jk2 = 1, nkgdimEns
          if(rstddev(jk1,jn).ne.0..and.rstddev(jk2,jn).ne.0.) then
            corns(jk1,jk2,jn) =  corns(jk1,jk2,jn)/(rstddev(jk1,jn)*rstddev(jk2,jn))
          else
            corns(jk1,jk2,jn) = 0.0d0
          endif
        enddo
      enddo
    enddo
!$OMP END PARALLEL DO

    dfact2 = 1.0d0/sqrt(dble(nens-1))
    do jn = 0, ntrunc
      dfact = 1.0d0/sqrt(2.0d0*dble(jn) + 1.0d0)
      do jk1 = 1, nkgdimEns
        rstddev(jk1,jn) = rstddev(jk1,jn)*dfact2*dfact
      enddo
    enddo

    ! Normalize to ensure correlations in horizontal and Multiply by sqrt(0.5) to make valid for m.ne.0
!$OMP PARALLEL DO PRIVATE (jk1,jn,dsummed)
    do jk1 = 1, nkgdimEns
      dsummed=0.0d0
      do jn = 0, ntrunc
        dsummed=dsummed + (rstddev(jk1,jn)**2)*((2.0d0*dble(jn))+1.0d0)/2.0d0
      enddo
      do jn = 0, ntrunc
        if(dsummed.gt.0.0d0) rstddev(jk1,jn)=rstddev(jk1,jn)*sqrt(0.5d0/dsummed)
      enddo
    enddo
!$OMP END PARALLEL DO

    call tmg_stop(3)
    write(*,*) 'finished computing correlations...'
    call flush(6)

  end subroutine calcCorrelations

!--------------------------------------------------------------------------
! CALCPOWERSPEC
!--------------------------------------------------------------------------
  subroutine calcPowerSpec(ensPerturbations,powerSpec)
    implicit none
    real*4,pointer, intent(in)  :: ensPerturbations(:,:,:,:)
    real*8,         intent(out) :: powerSpec(nkgdimEns,0:ntrunc)

    real*8  :: spectralState(nla,2,nkgdimEns),gridState(ni,nj,nkgdimEns)
    real*8  :: dfact,dfact2,dsummed

    integer :: jens,ila,jn,jm,jk

    powerSpec(:,:)=0.0d0

    do jens = 1, nens

      write(*,*) 'calcPowerSpec: processing member ',jens
      call flush(6)

      gridState(:,:,:)=ensPerturbations(:,:,:,jens)
      call gst_setID(gstID_nkgdimEns)
      call gst_reespe(spectralState,gridState)

!$OMP PARALLEL DO PRIVATE (jn,jm,dfact,ila,jk)
      do jn = 0, ntrunc
        do jm = 0, jn
          dfact = 2.0d0
          if (jm.eq.0) dfact = 1.0d0
          ila = gst_getNind(jm) + jn - jm
          do jk = 1, nkgdimEns
              powerSpec(jk,jn) = powerSpec(jk,jn) +     &
                 dfact*( spectralState(ila,1,jk)**2 + spectralState(ila,2,jk)**2 )
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

    enddo

    dfact2 = 1.0d0/sqrt(dble(nens-1))
    do jn = 0, ntrunc
      dfact = 1.0d0/sqrt(2.0d0*dble(jn) + 1.0d0)
      do jk = 1, nkgdimEns
        powerSpec(jk,jn) = powerSpec(jk,jn)*dfact2*dfact
      enddo
    enddo

    write(*,*) 'finished computing power spectrum...'
    call flush(6)

  end subroutine calcPowerSpec

!--------------------------------------------------------------------------
! WRITESTATS
!--------------------------------------------------------------------------
  subroutine writeStats(corns, rstddev, variableType, ptot_opt, &
                        theta_opt, waveBandIndex_opt, latBand_opt)
    implicit none

    real*8 :: corns(nkgdimEns,nkgdimEns,0:ntrunc),rstddev(nkgdimEns,0:ntrunc)
    integer, intent(in) :: variableType
    real*8, optional :: PtoT_opt(:,:,:),theta_opt(:,:)
    integer, optional :: waveBandIndex_opt
    integer, optional :: latBand_opt

    real*8 prcor(nkgdimEns,nkgdimEns)

    integer :: jn,ierr,ipak,jk,jl
    integer :: fstouv,fnom,fstfrm,fclos
    integer :: ip1,ip2,ip3,idatyp,idateo
    integer :: nulstats
    integer :: jvar1, jvar2
    integer :: nLevEns1, nLevEns2, nlevstart1, nlevend1, nlevstart2, nlevend2

    character(len=128) :: outfilename
    character(len=2) :: wbnum

    nulstats=0
    if ( nWaveBand == 1 ) then
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
    enddo

    do jn = 0, ntrunc
      ip2 = jn
      ierr = utl_fstecr(rstddev(:,jn),ipak,nulstats,idateo,0,0,nkgdimEns,1,1,   &
                        ip1,ip2,ip3,'X','SS','RSTDDEV ','X',0,0,0,0,idatyp,.true.)
    enddo

    
    ! Computing the total vertical correlation matrix
    do jk = 1, nkgdimEns
      do jl = 1, nkgdimEns
        prcor(jk,jl) = 0
        do jn = 0, ntrunc
          prcor(jk,jl) = prcor(jk,jl) + ((2*jn+1)*rstddev(jk,jn)*rstddev(jl,jn)*corns(jk,jl,jn))
        enddo
      enddo
    enddo

    do jk = 1, nkgdimEns
      do jl = 1, nkgdimEns
        if(prcor(jk,jk)*prcor(jl,jl) .gt. 0.0d0) then
          prcor(jk,jl) = prcor(jk,jl) / (sqrt(prcor(jk,jk)*prcor(jl,jl)))
        else
          prcor(jk,jl) = 0.0d0
        endif
      enddo
    enddo
    ip2 =0
    ierr = utl_fstecr(prcor(:,:),ipak,nulstats,idateo,0,0,nkgdimEns,nkgdimEns,1,   &
                      ip1,ip2,ip3,'X','ZV','CORVERT ','X',0,0,0,0,idatyp,.true.)

    ierr =  fstfrm(nulstats)
    ierr =  fclos (nulstats)

    !- For plotting purposes...
    if(.false.) then
    do jvar1 = 1, nvar3d
       do jvar2 = jvar1, nvar3d

          if ( nWaveBand == 1 ) then
             outfilename = "./vertCorrel_"//trim(nomvar3d(jvar1,variableType))//"_"//trim(nomvar3d(jvar2,variableType))//".txt"
          else
             outfilename = "./vertCorrel_"//trim(nomvar3d(jvar1,variableType))//"_"//trim(nomvar3d(jvar2,variableType))//"_"//wbnum//".txt"
          end if
          open (unit=99,file=outfilename,action="write",status="new")

          if (vnl_varLevelFromVarName(nomvar3d(jvar1,variableType)).eq.'MM') then
             nLevEns1 = nLevEns_M
          else
             nLevEns1 = nLevEns_T
          endif
          nLevStart1 = varLevOffset(jvar1)+ 1 
          nLevEnd1   = varLevOffset(jvar1)+ nLevEns1

          if (vnl_varLevelFromVarName(nomvar3d(jvar2,variableType)).eq.'MM') then
             nLevEns2 = nLevEns_M
          else
             nLevEns2 = nLevEns_T
          endif
          nLevStart2 = varLevOffset(jvar2)+ 1 
          nLevEnd2   = varLevOffset(jvar2)+ nLevEns2

          do jk = nLevStart1, nLevEnd1
             do jl = nLevStart2, nLevEnd2
                if ( jl == nLevEnd2 ) then 
                   write(99,'(2X,F6.4)')    prcor(jk,jl) ! Saut de ligne
                else
                   write(99,'(2X,F6.4,$)')  prcor(jk,jl)
                end if
             end do
          end do

          close(unit=99)
          
       end do
    end do
    end if

    write(*,*) 'finished writing statistics...'
    call flush(6)

  end subroutine writeStats

!--------------------------------------------------------------------------
! WRITEPRESSUREPROFILES
!--------------------------------------------------------------------------
  subroutine writePressureProfiles
    implicit none

    character(len=128) :: outfilename

    integer :: jk

    outfilename = "./pressureProfile_M.txt"
    open (unit=99,file=outfilename,action="write",status="new")
    do jk = 1, nLevEns_M
       write(99,'(I3,2X,F6.1)') jk, pressureProfile_M(jk)/100.d0
    end do
    close(unit=99)
       
    outfilename = "./pressureProfile_T.txt"
    open (unit=99,file=outfilename,action="write",status="new")
    do jk = 1, nLevEns_T
       write(99,'(I3,2X,F6.1)') jk, pressureProfile_T(jk)/100.d0
    end do
    close(unit=99)

    write(*,*) 'finished writing pressure profiles...'
    call flush(6)

  end subroutine writePressureProfiles

!--------------------------------------------------------------------------
! WRITESTDDEV
!--------------------------------------------------------------------------
  subroutine writeStddev(stddevZonAvg,stddevZonAvgUnbal,stddev3d,stddev3dUnbal)
    implicit none
    real*8 :: stddevZonAvg(:,:),stddevZonAvgUnbal(:,:),stddev3d(:,:,:),stddev3dUnbal(:,:,:)
    real*8 :: dfact,zbuf(ni,nj),zbufyz(nj,max(nLevEns_M,nLevens_T)),zbufy(nj)
    integer jlat,jlon,jlev,ierr,jvar,nLevEns
    integer fstouv,fnom,fstfrm,fclos
    integer ip1,ip2,ip3,idatyp,idateo,ipak,nip1_l(max(nLevEns_M,nLevens_T))
    integer :: nulstats

    nulstats=0
    ierr =  fnom  (nulstats,'./stddev.fst','RND',0)
    ierr =  fstouv(nulstats,'RND')

    ipak = -32
    idatyp = 5
    ip1 = 0
    ip2 = 0
    ip3 = nens
    idateo = 0

    ! do 3d variables
    do jvar=1,nvar3d
      if(vnl_varLevelFromVarName(nomvar3d(jvar,cvSpace)).eq.'MM') then
        nLevEns = nLevEns_M
        nip1_l(1:nLevEns_M)=nip1_M(1:nLevEns_M)
      else
        nLevEns = nLevEns_T
        nip1_l(1:nLevEns_T)=nip1_T(1:nLevEns_T)
      endif
      dfact=1.0d0

      do jlev=1,nlevEns
        do jlat=1,nj
          do jlon=1,ni
            zbuf(jlon,jlat)=dfact*stddev3d(jlon,jlat,varLevOffset(jvar)+jlev)
          enddo
        enddo
        ierr = utl_fstecr(zbuf,ipak,nulstats,idateo,0,0,ni,nj,1,nip1_l(jlev),ip2,ip3,   &
                          'E',nomvar3d(jvar,cvSpace),'STDDEV3D','G',0,0,0,0,idatyp,.true.)
      enddo

      do jlev=1,nlevEns
        do jlat=1,nj
          zbufyz(jlat,jlev)=dfact*stddevZonAvg(varLevOffset(jvar)+jlev,jlat)
        enddo
      enddo
      ierr = utl_fstecr(zbufyz(:,1:nLevEns),ipak,nulstats,idateo,0,0,1,nj,nlevEns,ip1,ip2,ip3,   &
                        'E',nomvar3d(jvar,cvSpace),'STDDEV  ','X',0,0,0,0,idatyp,.true.)

      if(nomvar3d(jvar,cvSpace).ne.nomvar3d(jvar,cvUnbalSpace)) then
        dfact=1.0d0

        do jlev=1,nlevEns
          do jlat=1,nj
            do jlon=1,ni
              zbuf(jlon,jlat)=dfact*stddev3dUnbal(jlon,jlat,varLevOffset(jvar)+jlev)
            enddo 
          enddo
          ierr = utl_fstecr(zbuf,ipak,nulstats,idateo,0,0,ni,nj,1,nip1_l(jlev),ip2,ip3,   &
                            'E',nomvar3d(jvar,cvUnbalSpace),'STDDEV3D','G',0,0,0,0,idatyp,.true.)
        enddo

        do jlev=1,nlevEns
          do jlat=1,nj
            zbufyz(jlat,jlev)=dfact*stddevZonAvgUnbal(varLevOffset(jvar)+jlev,jlat)
          enddo
        enddo
        ierr = utl_fstecr(zbufyz(:,1:nLevEns),ipak,nulstats,idateo,0,0,1,nj,nlevEns,ip1,ip2,ip3,   &
                          'E',nomvar3d(jvar,cvUnbalSpace),'STDDEV  ','X',0,0,0,0,idatyp,.true.)
      endif

    enddo

    ! now do 2D variables
    do jvar=1,nvar2d
      if(nomvar2d(jvar,cvSpace).eq.'P0') then
        dfact=1.0d0/1.0d2
      else
        dfact=1.0d0
      endif

      do jlat=1,nj
        do jlon=1,ni
          zbuf(jlon,jlat)=dfact*stddev3d(jlon,jlat,varLevOffset(nvar3d+1)+jvar)
        enddo
      enddo
      ierr = utl_fstecr(zbuf,ipak,nulstats,idateo,0,0,ni,nj,1,0,ip2,ip3,   &
                        'E',nomvar2d(jvar,cvSpace),'STDDEV3D','G',0,0,0,0,idatyp,.true.)

      do jlat=1,nj
        zbufy(jlat)=dfact*stddevZonAvg(varLevOffset(nvar3d+1)+jvar,jlat)
      enddo
      ierr = utl_fstecr(zbufy,ipak,nulstats,idateo,0,0,1,nj,1,ip1,ip2,ip3,   &
                        'E',nomvar2d(jvar,cvSpace),'STDDEV  ','X',0,0,0,0,idatyp,.true.)

      if(nomvar2d(jvar,cvSpace).ne.nomvar2d(jvar,cvUnbalSpace)) then
        if(nomvar2d(jvar,cvUnbalSpace).eq.'UP') then
          dfact=1.0d0/1.0d2
        else
          dfact=1.0d0
        endif

        do jlat=1,nj
          do jlon=1,ni
            zbuf(jlon,jlat)=dfact*stddev3dUnbal(jlon,jlat,varLevOffset(nvar3d+1)+jvar)
          enddo 
        enddo
        ierr = utl_fstecr(zbuf,ipak,nulstats,idateo,0,0,ni,nj,1,0,ip2,ip3,   &
                          'E',nomvar2d(jvar,cvUnbalSpace),'STDDEV3D','G',0,0,0,0,idatyp,.true.)

        do jlat=1,nj
          zbufy(jlat)=dfact*stddevZonAvgUnbal(varLevOffset(nvar3d+1)+jvar,jlat)
        enddo
        ierr = utl_fstecr(zbufy,ipak,nulstats,idateo,0,0,1,nj,1,ip1,ip2,ip3,   &
                          'E',nomvar2d(jvar,cvUnbalSpace),'STDDEV  ','X',0,0,0,0,idatyp,.true.)
      endif

    enddo

    ierr =  fstfrm(nulstats)
    ierr =  fclos (nulstats)

    write(*,*) 'finished writing stddev...'
    call flush(6)

  end subroutine writeStddev

!--------------------------------------------------------------------------
! WRITESTDDEV2
!--------------------------------------------------------------------------
  subroutine writeStddev2(stddevZonAvg,stddev3d)
    implicit none
    real*8 :: stddevZonAvg(:,:),stddev3d(:,:,:)
    real*8 :: dfact,zbuf(ni,nj),zbufyz(nj,max(nLevEns_M,nLevens_T)),zbufy(nj)
    integer jlat,jlon,jlev,ierr,jvar,nLevEns
    integer fstouv,fnom,fstfrm,fclos
    integer ip1,ip2,ip3,idatyp,idateo,ipak,nip1_l(max(nLevEns_M,nLevens_T))
    integer :: nulstats

    nulstats=0
    ierr =  fnom  (nulstats,'./stddev.fst','RND',0)
    ierr =  fstouv(nulstats,'RND')

    ipak = -32
    idatyp = 5
    ip1 = 0
    ip2 = 0
    ip3 = nens
    idateo = 0

    ! do 3d variables
    do jvar=1,nvar3d
      if(vnl_varLevelFromVarName(nomvar3d(jvar,cvSpace)).eq.'MM') then
        nLevEns = nLevEns_M
        nip1_l(1:nLevEns_M)=nip1_M(1:nLevEns_M)
      else
        nLevEns = nLevEns_T
        nip1_l(1:nLevEns_T)=nip1_T(1:nLevEns_T)
      endif

      do jlev=1,nlevEns
        do jlat=1,nj
          dfact=1.0d0
          do jlon=1,ni
            zbuf(jlon,jlat)=dfact*stddev3d(jlon,jlat,varLevOffset(jvar)+jlev)
          enddo
        enddo
        ierr = utl_fstecr(zbuf,ipak,nulstats,idateo,0,0,ni,nj,1,nip1_l(jlev),ip2,ip3,   &
                          'E',nomvar3d(jvar,cvSpace),'STDDEV3D','G',0,0,0,0,idatyp,.true.)
      enddo

      do jlev=1,nlevEns
        do jlat=1,nj
          dfact=1.0d0
          zbufyz(jlat,jlev)=dfact*stddevZonAvg(varLevOffset(jvar)+jlev,jlat)
        enddo
      enddo
      ierr = utl_fstecr(zbufyz(:,1:nLevEns),ipak,nulstats,idateo,0,0,1,nj,nlevEns,ip1,ip2,ip3,   &
                        'E',nomvar3d(jvar,cvSpace),'STDDEV  ','X',0,0,0,0,idatyp,.true.)

    enddo

    ! now do 2D variables
    do jvar=1,nvar2d
      if(nomvar2d(jvar,cvSpace).eq.'P0') then
        dfact=1.0d0/1.0d2
      else
        dfact=1.0d0
      endif

      do jlat=1,nj
        do jlon=1,ni
          zbuf(jlon,jlat)=dfact*stddev3d(jlon,jlat,varLevOffset(nvar3d+1)+jvar)
        enddo
      enddo
      ierr = utl_fstecr(zbuf,ipak,nulstats,idateo,0,0,ni,nj,1,0,ip2,ip3,   &
                        'E',nomvar2d(jvar,cvSpace),'STDDEV3D','G',0,0,0,0,idatyp,.true.)

      do jlat=1,nj
        zbufy(jlat)=dfact*stddevZonAvg(varLevOffset(nvar3d+1)+jvar,jlat)
      enddo
      ierr = utl_fstecr(zbufy,ipak,nulstats,idateo,0,0,1,nj,1,ip1,ip2,ip3,   &
                        'E',nomvar2d(jvar,cvSpace),'STDDEV  ','X',0,0,0,0,idatyp,.true.)

    enddo

    ierr =  fstfrm(nulstats)
    ierr =  fclos (nulstats)

    write(*,*) 'finished writing stddev...'
    call flush(6)

  end subroutine writeStddev2

!--------------------------------------------------------------------------
! WRITESTDDEVBAL
!--------------------------------------------------------------------------
  subroutine writeStddevBal(stddevZonAvgBal,stddev3dBal)
    implicit none
    real*8 :: stddevZonAvgBal(:,:),stddev3dBal(:,:,:)
    real*8 :: dfact,zbuf(ni,nj),zbufyz(nj,max(nLevEns_M,nLevens_T)),zbufy(nj)
    integer jlat,jlon,jlev,ierr,jvar,nLevEns
    integer fstouv,fnom,fstfrm,fclos
    integer ip1,ip2,ip3,idatyp,idateo,ipak,nip1_l(max(nLevEns_M,nLevens_T))
    integer :: nulstats
    integer,parameter :: nvar3d=1,nvar2d=1
    character*4 :: nomvar3dBal(nvar3d),nomvar2dBal(nvar2d)

    nomvar3dBal(1)='TB'
    nomvar2dBal(1)='PB'

    nulstats=0
    ierr =  fnom  (nulstats,'./stddev_balanced.fst','RND',0)
    ierr =  fstouv(nulstats,'RND')

    ipak = -32
    idatyp = 5
    ip1 = 0
    ip2 = 0
    ip3 = nens
    idateo = 0

    ! do 3d variables
    do jvar=1,nvar3d
      nLevEns = nLevEns_T
      nip1_l(1:nLevEns_T)=nip1_T(1:nLevEns_T)
      dfact=1.0d0

      do jlev=1,nlevEns
        do jlat=1,nj
          do jlon=1,ni
            zbuf(jlon,jlat)=dfact*stddev3dBal(jlon,jlat,varLevOffset(jvar)+jlev)
          enddo
        enddo
        ierr = utl_fstecr(zbuf,ipak,nulstats,idateo,0,0,ni,nj,1,nip1_l(jlev),ip2,ip3,   &
                          'E',nomvar3dBal(jvar),'STDDEV3D','G',0,0,0,0,idatyp,.true.)
      enddo

      do jlev=1,nlevEns
        do jlat=1,nj
          zbufyz(jlat,jlev)=dfact*stddevZonAvgBal(varLevOffset(jvar)+jlev,jlat)
        enddo
      enddo
      ierr = utl_fstecr(zbufyz(:,1:nLevEns),ipak,nulstats,idateo,0,0,1,nj,nlevEns,ip1,ip2,ip3,   &
                        'E',nomvar3dBal(jvar),'STDDEV  ','X',0,0,0,0,idatyp,.true.)

    enddo

    ! now do 2D variables
    do jvar=1,nvar2d
      dfact=1.0d0/1.0d2

      do jlat=1,nj
        do jlon=1,ni
          zbuf(jlon,jlat)=dfact*stddev3dBal(jlon,jlat,varLevOffset(nvar3d+1)+jvar)
        enddo
      enddo
      ierr = utl_fstecr(zbuf,ipak,nulstats,idateo,0,0,ni,nj,1,0,ip2,ip3,   &
                        'E',nomvar2dBal(jvar),'STDDEV3D','G',0,0,0,0,idatyp,.true.)

      do jlat=1,nj
        zbufy(jlat)=dfact*stddevZonAvgBal(varLevOffset(nvar3d+1)+jvar,jlat)
      enddo
      ierr = utl_fstecr(zbufy,ipak,nulstats,idateo,0,0,1,nj,1,ip1,ip2,ip3,   &
                        'E',nomvar2dBal(jvar),'STDDEV  ','X',0,0,0,0,idatyp,.true.)

    enddo

    ierr =  fstfrm(nulstats)
    ierr =  fclos (nulstats)

    write(*,*) 'finished writing stddev...'
    call flush(6)

  end subroutine writeStddevBal

!--------------------------------------------------------------------------
! SPECTRALFILTER
!--------------------------------------------------------------------------
  subroutine spectralFilter(ensPerturbations,nlev,waveBandIndex_opt)
    implicit none
    real*4,pointer :: ensPerturbations(:,:,:,:)
    integer, intent(in) :: nlev
    integer, optional, intent(in) :: waveBandIndex_opt

    real*8  :: spectralState(nla,2,nlev)
    real*8  :: member(ni,nj,nlev)
    integer :: jens, jk, jn, jm, ila

    real(8), allocatable :: ResponseFunction(:)
    real(8) :: waveLength

    character(len=128) :: outfilename
    character(len=2) :: wbnum

    if ( nWaveBand /= 1 ) then
       write(*,*) 'Bandpass filtering step'
       if (.not. present(waveBandIndex_opt)) then
          write(*,*) 'Error: No waveBandIndex was supplied!!!'
          call utl_abort('calbmatrix_glb')
       end if
       allocate(ResponseFunction(0:ntrunc))
       write(wbnum,'(I2.2)') waveBandIndex_opt
       outfilename = "./ResponseFunction_"//wbnum//".txt"
       open (unit=99,file=outfilename,action="write",status="new")
       do jn = 0, ntrunc
          ResponseFunction(jn) = spf_filterResponseFunction(dble(jn),waveBandIndex_opt, waveBandPeaks, nWaveBand)
          if ( jn /= 0) then
             waveLength=4.d0*asin(1.d0)*ra/dble(jn)
          else
             waveLength=0.d0
          end if
          write(* ,'(I4,2X,F7.1,2X,F5.3)') jn, waveLength/1000.d0, ResponseFunction(jn)
          write(99,'(I4,2X,F7.1,2X,F5.3)') jn, waveLength/1000.d0, ResponseFunction(jn)
       end do
       close(unit=99)
    end if

    do jens=1,nens
      member(:,:,:)=dble(ensPerturbations(:,:,:,jens))     
      if(nlev.eq.nkgdimEns) then
        call gst_setID(gstID_nkgdimEns)
      elseif(nlev.eq.nLevEns_T+1) then
        call gst_setID(gstID_nLevEns_T_P1)
      else
        write(*,*) 'spectralFilter: nlev = ',nlev
        call utl_abort('spectralFilter: spectral transform not initialized for this number of levels')
      endif
      call gst_reespe(spectralState,member)
      if ( nWaveBand /= 1 ) then
!$OMP PARALLEL DO PRIVATE (jk,jn,jm,ila)
         do jk = 1, nlev
            do jn = 0, ntrunc
               do jm = 0, jn
                  ila = gst_getnind(jm,gstID_nkgdimEns)+jn-jm
                  spectralState(ila,1,jk) = spectralState(ila,1,jk) * ResponseFunction(jn)
                  spectralState(ila,2,jk) = spectralState(ila,2,jk) * ResponseFunction(jn)
               end do
            end do
         end do
!$OMP END PARALLEL DO
      end if
      call gst_speree(spectralState,member)
      ensPerturbations(:,:,:,jens)=sngl(member(:,:,:))
    enddo

    if ( nWaveBand /= 1 ) then
       deallocate(ResponseFunction)
    end if

    write(*,*) 'finished applying spectral filter...'
    call flush(6)

  end subroutine spectralFilter

!--------------------------------------------------------------------------
! CALCTHETA
!--------------------------------------------------------------------------
  subroutine calcTheta(ensPerturbations,theta)
    implicit none
    real*4,pointer :: ensPerturbations(:,:,:,:)
    real*8  :: theta(:,:)
    real*8 zchipsi(nLevEns_M,nj), zpsipsi(nLevEns_M,nj)
    real*4, pointer :: psi_ptr(:,:,:),chi_ptr(:,:,:)
    integer :: jlat,jlon,jlev,jens

    theta(:,:) = 0.0d0
    zchipsi(:,:) = 0.0d0
    zpsipsi(:,:) = 0.0d0

    do jens = 1,nens
      psi_ptr => ensPerturbations(:,:,1:nlevEns_M,jens)
      chi_ptr => ensPerturbations(:,:,(nlevEns_M+1):(2*nlevEns_M),jens)

      ! update zchipsi and zpsipsi covariances
      do jlat = 1, nj
        do jlon = 1, ni
          do jlev = 1, nLevEns_M
            zpsipsi(jlev,jlat) = zpsipsi(jlev,jlat) + psi_ptr(jlon,jlat,jlev) * psi_ptr(jlon,jlat,jlev)
            zchipsi(jlev,jlat) = zchipsi(jlev,jlat) + chi_ptr(jlon,jlat,jlev) * psi_ptr(jlon,jlat,jlev)
          enddo
        enddo
      enddo
    enddo

    !  calculate THETA
    do jlat = 1, nj
      do jlev = 1, nLevEns_M
        theta(jlev,jlat) = atan(-zchipsi(jlev,jlat) / zpsipsi(jlev,jlat))
      enddo
    enddo

    write(*,*) 'finished computing theta...'
    call flush(6)

  end subroutine calcTheta

!--------------------------------------------------------------------------
! CALCPTOT
!--------------------------------------------------------------------------
  subroutine calcPtoT(ensPerturbations,PtoT)
    implicit none
    real*4,pointer :: ensPerturbations(:,:,:,:)
    real*8  :: PtoT(:,:,:)

    real*8  :: spectralState(nla,2,nLevEns_M),spBalancedP(nla,2,nlevEns_M),balancedP(ni,nj,nlevEns_M),psi(ni,nj,nLevEns_M)
    real*4, pointer :: tt_ptr(:,:,:),ps_ptr(:,:)
    INTEGER JENS, IENS, JK1, JK2, JLA, JN, JM, ILA, JLEV
    INTEGER IERR, JFILE, JK, JLAT, ILON, JLON, JB, NLATBAND
    INTEGER IBND1,IBND2,JPNLATBND,ILAT
    PARAMETER (JPNLATBND = 3)
    REAL*8 ZFACT,ZMAXI,ZWT,ZPS,zlat(nj)
    REAL*8 ZFACT2,ZFACTTOT
    REAL*8 ZM1(NLEVENS_T+1,NLEVENS_M,JPNLATBND), ZM2(NLEVPTOT,NLEVPTOT,JPNLATBND)
    REAL*8 ZPTOTBND(NLEVENS_T+1,NLEVENS_M)
    REAL*8 ZM2INV(NLEVPTOT,NLEVPTOT,JPNLATBND),ZWORK(NLEVPTOT*NLEVPTOT),ZDET,ZEPS
    REAL*8  DLA2, DL1SA2
    REAL*8  DLLATMIN(JPNLATBND), DLLATMAX(JPNLATBND)
    REAL*8  DLLATMID(JPNLATBND)
    REAL*8  ZLC,ZTLEN,ZR,ZCORR,ZPRES1,ZPRES2
    real*8 zeigwrk(4*nlevPtoT),zeigen(nlevPtoT,nlevPtoT),zeigenv(nlevPtoT)
    real*8 zeigenvi(nlevPtoT)
    real   zfix
    integer iwork,info

    DATA DLLATMIN / -60.0D0, -30.0D0, 30.0D0 /
    DATA DLLATMAX / -30.0D0,  30.0D0, 60.0D0 /
    DATA DLLATMID / -45.0D0,  00.0D0, 45.0D0 /

    DLA2 = DBLE(RA)*DBLE(RA)
    DL1SA2 = 1.D0/DLA2

!   1. Initialize P_to_T, ZM1, ZM2

    ZFACTTOT = 0.0D0
    DO JLAT = 1, NJ
      ZFACTTOT = ZFACTTOT + cos(GST_GETRLATI(JLAT))
    ENDDO
    ZFACTTOT = NJ/ZFACTTOT

    PtoT(:,:,:) = 0.0d0
    ZM1(:,:,:) = 0.0d0
    ZPTOTBND(:,:) = 0.0d0
    ZM2(:,:,:) = 0.0d0

    do jens = 1,nens

      write(*,*) 'calcPtoT: processing member ',jens
      call flush(6)

      psi(:,:,:)=ensPerturbations(:,:,1:nlevEns_M,jens)
      call gst_setID(gstID_nLevEns_M)
      call gst_reespe(spectralState,psi)
      CALL calcBalancedP(spectralState,spBalancedP)
      call gst_speree(spBalancedP,balancedP)

      tt_ptr  => ensPerturbations(:,:,(2*nLevEns_M+1):(2*nLevEns_M+nLevEns_T),jens)
      ps_ptr  => ensPerturbations(:,:,2*nLevEns_M+2*nLevEns_T+1,jens)

      DO JLAT = 1, NJ
        zlat(jlat)=GST_GETRLATI(JLAT)
      enddo

!$OMP PARALLEL DO PRIVATE (JK1,JB,JLAT,ZFACT,JLON,JK2)
      DO JK1 = 1, (nLevEns_T+1)
        DO JB=1,JPNLATBND
          DO JLAT = 1, NJ
            if ((ZLAT(JLAT) .gt. 2.D0*MPC_PI_R8*DLLATMIN(JB)/360.D0)   &
          .and. (ZLAT(JLAT) .le. 2.D0*MPC_PI_R8*DLLATMAX(JB)/360.D0)) then
              ZFACT = cos(ZLAT(JLAT))*ZFACTTOT
              DO JLON = 1, NI

                ! update ZM1 = sum_over_t_x_y[vec(T lnPs) vec(P_b)^T]
                DO JK2 = 1, nLevEns_M
                  IF(JK1.LE.nLevEns_T) THEN
                    zm1(jk1,jk2,jb) = zm1(jk1,jk2,jb) + zfact * tt_ptr(jlon,jlat,jk1) * balancedP(jlon,jlat,jk2)
                  ELSE
                    zm1(jk1,jk2,jb) = zm1(jk1,jk2,jb) + zfact * ps_ptr(jlon,jlat) * balancedP(jlon,jlat,jk2)
                  ENDIF
                ENDDO

                ! update ZM2 = sum_over_t_x_y[vec(P_b) vec(P_b)^T]
                IF(JK1.LE.NLEVPTOT) THEN
                  DO JK2 = 1, NLEVPTOT
                    zm2(jk1,jk2,jb) = zm2(jk1,jk2,jb) + zfact * balancedP(jlon,jlat,jk1) * balancedP(jlon,jlat,jk2)
                  ENDDO
                ENDIF

              ENDDO
            endif
          ENDDO
        ENDDO ! Loop on JPNLATBND
      END DO  ! Loop on JK1
!$OMP END PARALLEL DO

    ENDDO

! SET ZM1,ZM2 EQUAL FOR ALL THREE REGIONS
    DO JK1 = 1, NLEVPTOT
      DO JK2 = 1, NLEVPTOT
        ZM2(JK1,JK2,1)=ZM2(JK1,JK2,1)+ZM2(JK1,JK2,3)
        ZM2(JK1,JK2,2)=ZM2(JK1,JK2,1)
        ZM2(JK1,JK2,3)=ZM2(JK1,JK2,1)
      ENDDO
    ENDDO
    DO JK1 = 1, (nLevEns_T+1)
      DO JK2 = 1, NLEVPTOT
        ZM1(JK1,JK2,1)=ZM1(JK1,JK2,1)+ZM1(JK1,JK2,3)
        ZM1(JK1,JK2,2)=ZM1(JK1,JK2,1)
        ZM1(JK1,JK2,3)=ZM1(JK1,JK2,1)
      ENDDO
    ENDDO

    DO JK1=1,NLEVPTOT
      DO JK2=1,NLEVPTOT
        ZEIGEN(JK1,JK2)=ZM2(JK1,JK2,1)
      ENDDO
    ENDDO
    IWORK=4*NLEVPTOT
    CALL DSYEV('V','U',NLEVPTOT,ZEIGEN,NLEVPTOT,ZEIGENV,ZEIGWRK,IWORK,INFO)

    write(*,*) 'calcPtot: info=',info
    write(*,*) 'calcPtot: eigen values=',zeigenv(:)

    do JK1=1,NLEVPTOT
      if (ZEIGENV(JK1).gt.0.0d0) then
        ZEIGENVI(JK1)=1.0d0/ZEIGENV(JK1)
      else
        ZEIGENVI(JK1)=0.0d0
      endif
    enddo

    DO JK1=1,NLEVPTOT
      DO JK2=1,NLEVPTOT
        ZM2INV(JK1,JK2,1)=0.0d0
        DO JK=1,NLEVPTOT
          ZM2INV(JK1,JK2,1)=ZM2INV(JK1,JK2,1)+ZEIGEN(JK1,JK)*ZEIGENVI(JK)*ZEIGEN(JK2,JK)
        ENDDO
      ENDDO
    ENDDO

!    write(*,*) 'zm1=',zm1(:,:,1)
!    write(*,*) 'zm2=',zm2(:,:,1)
!    write(*,*) 'zm2inv=',zm2inv(:,:,1)

! Calculate A = ZM1*inv(ZM2)
!
    DO JK1 = 1, (nLevEns_T+1)
      DO JK2 = 1, NLEVPTOT
        DO JK = 1, NLEVPTOT
          ZPTOTBND(JK1,JK2) = ZPTOTBND(JK1,JK2) + ZM1(JK1,JK,1) * ZM2INV(JK,JK2,1)
        ENDDO
      ENDDO
    ENDDO
!
    DO JK1 = 1, nLevEns_T+1
      DO JK2 = 1, NLEVPTOT
        DO JLAT = 1,NJ
          PTOT(JK1,JK2,JLAT) = ZPTOTBND(JK1,JK2)
        ENDDO
      ENDDO
    ENDDO

    write(*,*) 'finished computing PtoT...'
    call flush(6)

  end subroutine calcPtoT

!--------------------------------------------------------------------------
! REMOVEGLOBALMEAN
!--------------------------------------------------------------------------
  subroutine removeGlobalMean(ensPerturbations)
    implicit none
    integer :: jlon,jlat,jlev,jens
    real*4  :: ensPerturbations(:,:,:,:)
    real*8  :: dmean

!$OMP PARALLEL DO PRIVATE (JENS,JLEV,JLAT,JLON,DMEAN)
    do jens=1,nens
      do jlev=1,nkgdimEns
        dmean=0.0d0
        do jlat=1,nj
          do jlon=1,ni
            dmean=dmean+ensPerturbations(jlon,jlat,jlev,jens)
          enddo
        enddo
        dmean=dmean/(dble(ni)*dble(nj))
        do jlat=1,nj
          do jlon=1,ni
            ensPerturbations(jlon,jlat,jlev,jens)=ensPerturbations(jlon,jlat,jlev,jens)-dmean
          enddo
        enddo
      enddo
    enddo
!$OMP END PARALLEL DO

    write(*,*) 'finished removing global mean...'
    call flush(6)

  end subroutine removeGlobalMean

!--------------------------------------------------------------------------
! CALCZONAVG
!--------------------------------------------------------------------------
  subroutine calcZonAvg(fieldsZonAvg,fields3D,nlev)
    implicit none

    integer :: jlon,jlat,jlev,nlev
    real*8  :: fieldsZonAvg(:,:),fields3D(:,:,:),dfact

    fieldsZonAvg(:,:)=0.0d0
    dfact=1.0d0/dble(ni)
!$OMP PARALLEL DO PRIVATE (JLEV,JLAT,JLON)
      do jlev=1,nlev
        do jlat=1,nj
          do jlon=1,ni
            fieldsZonAvg(jlev,jlat)=fieldsZonAvg(jlev,jlat)+dfact*fields3D(jlon,jlat,jlev)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

    write(*,*) 'finished computing the zonal average...'
    call flush(6)

  end subroutine calcZonAvg

!--------------------------------------------------------------------------
! CALCSTDDEV3D
!--------------------------------------------------------------------------
  subroutine calcStddev3d(ensPerturbations,stddev3d,nlev)
    implicit none

    integer :: jlon,jlat,jlev,jens,nlev
    real*8  :: dnens,stddev3d(:,:,:)
    real*4  :: ensPerturbations(:,:,:,:)

    write(*,*) 'started computing the stddev...'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call flush(6)

    stddev3d(:,:,:)=0.0d0
    dnens=1.0d0/dble(nens-1)
!$OMP PARALLEL DO PRIVATE (JLEV,JENS,JLAT,JLON)
      do jlev=1,nlev
        do jens=1,nens
          do jlat=1,nj
            do jlon=1,ni
              stddev3d(jlon,jlat,jlev)=stddev3d(jlon,jlat,jlev)+ensPerturbations(jlon,jlat,jlev,jens)**2
            enddo
          enddo
        enddo
        do jlat=1,nj
          do jlon=1,ni
            if(stddev3d(jlon,jlat,jlev).gt.0.0d0) then
              stddev3d(jlon,jlat,jlev)=sqrt(stddev3d(jlon,jlat,jlev)*dnens)
            else
              stddev3d(jlon,jlat,jlev)=0.0d0
            endif
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

    write(*,*) 'finished computing the stddev...'
    call flush(6)
  
  end subroutine calcStddev3d

!--------------------------------------------------------------------------
! CALCBALANCEDP
!--------------------------------------------------------------------------
  subroutine calcBalancedP(sppsi,spgz)
    implicit none

    real*8 :: sppsi(:,:,:),spgz(:,:,:)
    real*8 :: spvor(nla,2,nlevEns_M)
    integer  ia, ib, ji, jm, jlev,jlat
    real*8 :: zn,zm,zenm,zenmp1,zcon,dl1sa2
    ! constants
    real*8             :: rday
    real*8             :: rsiyea
    real*8             :: rsiday
    real*8             :: romega

    ! some constants
    RDAY=86400.D0
    RSIYEA=365.25D0*RDAY*2.*MPC_PI_R8/6.283076D0
    RSIDAY=RDAY/(1.D0+RDAY/RSIYEA)
    ROMEGA=2.D0*MPC_PI_R8/RSIDAY

    ! convert PSI to vorticity 
    dl1sa2   = 1.0d0/(dble(ra)*dble(ra))
    do jlev = 1, nlevEns_M
      do jlat = 1, nla
        spvor(jlat,1,jlev) = sppsi(jlat,1,jlev)*dl1sa2*gst_getRnnp1(jlat)
        spvor(jlat,2,jlev) = sppsi(jlat,2,jlev)*dl1sa2*gst_getRnnp1(jlat)
      enddo
    enddo

    ! ensure input field is zero for spectral component (0,0)
    do jlev = 1, nlevEns_M
      if(spvor(1,1,jlev).ne.0.D0) then
        spvor(1,1,jlev) = 0.0D0
      endif
      if(spvor(1,2,jlev).ne.0.D0) then
        spvor(1,2,jlev) = 0.0D0
      endif
    enddo

    ! initialize outout field to zero
    spgz(:,:,:)=0.0d0

    ! loop over levels and zonal wavenumbers
    ! n.b.: at the tip of the triangle, no contributions
    
    zcon = -2.D0*romega*ra**2
    do jlev = 1, nlevEns_M

      ! the base address ia will point to the spherical harmonic
      ! coefficient (m,m), in the input field
      ia = 1
      do jm = 0, ntrunc-1
        ib = ia + ntrunc - jm
        zm = dble(jm)

        ! at the base, contributions from n+1 coeff only
        zn = zm
        zenmp1 = sqrt ( ((zn+1)**2-zm**2)/(4.D0*(zn+1)**2-1.D0) )
        spgz(ia,1,jlev)=zcon*spvor(ia+1,1,jlev)*zenmp1/((zn+1.0D0)**2)
        spgz(ia,2,jlev)=zcon*spvor(ia+1,2,jlev)*zenmp1/((zn+1.0D0)**2)

        zn = zn+1
        do ji = ia+1, ib-1
          zenm = sqrt ( (zn**2-zm**2)/(4.D0*zn**2-1.D0) )
          zenmp1 = sqrt ( ((zn+1)**2-zm**2)/(4.D0*(zn+1)**2-1.D0) )
          spgz(ji,1,jlev)=spvor(ji-1,1,jlev)*zenm/(zn**2)
          spgz(ji,2,jlev)=spvor(ji-1,2,jlev)*zenm/(zn**2)
          spgz(ji,1,jlev)=zcon*(spgz(ji,1,jlev)+spvor(ji+1,1,jlev)*zenmp1/((zn+1.0D0)**2))
          spgz(ji,2,jlev)=zcon*(spgz(ji,2,jlev)+spvor(ji+1,2,jlev)*zenmp1/((zn+1.0D0)**2))
          zn = zn + 1.0D0
        enddo

        ! at the top, contributions from n-1 coeff only
        zenm = sqrt ( (zn**2-zm**2)/(4.D0*zn**2-1.D0) )
        spgz(ib,1,jlev) = zcon*spvor(ib-1,1,jlev)*zenm/(zn**2)
        spgz(ib,2,jlev) = zcon*spvor(ib-1,2,jlev)*zenm/(zn**2)
        ia = ib + 1
      enddo
    enddo

    ! ensure correct value for mass spectral-coefficient for m=n=0
    do jlev = 1, nlevens_M
      spgz(1,1,jlev) = 0.0D0
      spgz(1,2,jlev) = 0.0D0
    enddo

  end subroutine calcBalancedP

!--------------------------------------------------------------------------
! NORMALIZED3D
!--------------------------------------------------------------------------
  subroutine normalize3d(ensPerturbations,stddev3d)
    implicit none

    integer :: jlon,jlat,jlev,jens
    real*8  :: dfact,stddev3d(:,:,:)
    real*4  :: ensPerturbations(:,:,:,:)

!$OMP PARALLEL DO PRIVATE (JLEV,JENS,JLAT,JLON,DFACT)
      do jlev=1,nkgdimEns
        do jlat=1,nj
          do jlon=1,ni
            if(stddev3d(jlon,jlat,jlev).gt.0.0d0) then
              dfact=1.0d0/stddev3d(jlon,jlat,jlev)
            else
              dfact=0.0d0
            endif
            do jens=1,nens
              ensPerturbations(jlon,jlat,jlev,jens)=ensPerturbations(jlon,jlat,jlev,jens)*dfact
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

    write(*,*) 'finished normalizing by stddev3D...'
    call flush(6)
  
  end subroutine normalize3d

!--------------------------------------------------------------------------
! MULTIPLY3D
!--------------------------------------------------------------------------
  subroutine multiply3d(ensPerturbations,stddev3d,nlev)
    implicit none

    integer :: jlon,jlat,jlev,jens,nlev
    real*8  :: stddev3d(:,:,:)
    real*4  :: ensPerturbations(:,:,:,:)

!$OMP PARALLEL DO PRIVATE (JLEV,JENS,JLAT,JLON)
      do jens=1,nens
        do jlev=1,nlev
          do jlat=1,nj
            do jlon=1,ni
              ensPerturbations(jlon,jlat,jlev,jens)=ensPerturbations(jlon,jlat,jlev,jens)*stddev3d(jlon,jlat,jlev)
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

    write(*,*) 'finished multiplying by stddev3D...'
    call flush(6)
  
  end subroutine multiply3d

!--------------------------------------------------------------------------
! READENSEMBLE
!--------------------------------------------------------------------------
  subroutine readEnsemble(ensPerturbations)
    implicit none

    real*4 :: ensPerturbations(:,:,:,:)

    integer :: jens, fclos, fnom, fstfrm, fstouv, ierr
    integer :: nulens(nens)

    nulens(:)=0
    do jens = 1, nens
       ierr = fnom(nulens(jens),cflensin(jens),'RND+OLD+R/O',0)
       ierr = fstouv(nulens(jens),'RND+OLD')
    end do

    write(*,*) 'Before reading the ensemble:'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call flush(6)

!$OMP PARALLEL DO PRIVATE (jens)
    do jens = 1, nens
        write(*,*) 'Reading ensemble member:',trim(cflensin(jens))
        call flush(6)
        call readEnsembleMember(ensPerturbations,nulens,jens)
        write(*,*) 'done reading member ',jens
        call flush(6)
    end do
!$OMP END PARALLEL DO

    write(*,*) 'After reading the ensemble:'
    write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    call flush(6)

    do jens = 1, nens
       ierr =  fstfrm(nulens(jens))
       ierr =  fclos (nulens(jens))
    end do

    write(*,*) 'finished reading ensemble members...'
    call flush(6)

  end subroutine readEnsemble

!--------------------------------------------------------------------------
! READENSEMBLEMEMBER
!--------------------------------------------------------------------------
  subroutine readEnsembleMember(ensPerturbations,nulens,jens)
    implicit none

    real*4 :: ensPerturbations(:,:,:,:)
    integer, intent(in) :: jens
    integer, intent(in) :: nulens(nens)
    
    integer :: stamp_in
    real*4 :: gd2d(ni,nj)
    real*8 :: rmsknt,rmbtpa,r1sa
    real*8 :: rhumin = 2.5d-6
    integer :: jlon,jlat,jlev
    integer :: fstlir
    integer :: ngposituu,ngpositvv,ngposittt,ngpositq,ngpositps,ngposittg

    ! standard file variables
    integer ini,inj,ink,ip1,ip2,ip3,ierr,idateo,ikey
    character(len=2)   :: cltypvar
    character(len=1)   :: clgrtyp
    character(len=4)   :: clnomvar
    character(len=12)  :: cletiket

    ! this should come from state vector object
    ngposituu=1
    ngpositvv=1+1*nLevEns_M
    ngposittt=1+2*nLevEns_M
    ngpositq =1+2*nLevEns_M+1*nLevEns_T
    ngpositps=1+2*nLevEns_M+2*nLevEns_T
    ngposittg=2+2*nLevEns_M+2*nLevEns_T

    ! some physical constants
    rmsknt = 1.d0/1.94246d0
    rmbtpa = 1.0d2
    r1sa=1.d0/6371229.d0

    ! read in raw ensemble (UU,VV,TT,P0,LQ (convert HU to LQ) - covariances)
    ip2 = 6
!    ip2 = -1
    ip3=-1
    idateo = -1
    cltypvar = ' '
    cletiket = ' '

    clnomvar = 'P0' 
    ikey=fstlir(gd2d,nulens(jens),ini,inj,ink,idateo,cletiket,-1,ip2,ip3, cltypvar,clnomvar)
    if(ikey.lt.0) call utl_abort('readEnsembleMember: Problem with P0 ENS')
    do jlat=1,nj
       do jlon=1,ni
          ensPerturbations(jlon,jlat,ngpositps,jens)= gd2d(jlon,jlat)*rmbtpa
       enddo
    enddo
    
    do jlev=1,nLevEns_T
       clnomvar = 'TT'
       ikey=fstlir(gd2d,nulens(jens),ini,inj,ink,idateo,cletiket,nip1_T(jlev),ip2,ip3,cltypvar,clnomvar)
       if(ikey.lt.0) then
          write(*,*) idateo,cletiket,nip1_T(jlev),ip2,ip3,cltypvar,clnomvar
          call utl_abort('readEnsembleMember: Problem with TT ENS')
       endif
       call flush(6)
       do jlat=1,nj
          do jlon=1,ni
             ensPerturbations(jlon,jlat,jlev-1+ngposittt,jens)= gd2d(jlon,jlat)
          enddo
       enddo
    enddo
    
    do jlev=1,nLevEns_T
       clnomvar = 'HU' 
       ikey=fstlir(gd2d,nulens(jens),ini,inj,ink,idateo,cletiket,nip1_T(jlev),ip2,ip3,cltypvar,clnomvar)
       if(ikey.lt.0) then
          clnomvar = 'LQ' 
          ikey=fstlir(gd2d,nulens(jens),ini,inj,ink,idateo,cletiket,nip1_T(jlev),ip2,ip3,cltypvar,clnomvar)
          if(ikey.lt.0) then
             write(*,*) idateo,cletiket,nip1_T(jlev),ip2,ip3,cltypvar,clnomvar
             call utl_abort('readEnsembleMember: Problem with HU and LQ ENS')
          else
             do jlat=1,nj
                do jlon=1,ni
                   ensPerturbations(jlon,jlat,jlev-1+ngpositq,jens)= gd2d(jlon,jlat)
                enddo
             enddo
          endif
       else
          do jlat=1,nj
             do jlon=1,ni
                ensPerturbations(jlon,jlat,jlev-1+ngpositq,jens)= log(max(gd2d(jlon,jlat),real(rhumin,4)))
             enddo
          enddo
       endif
    enddo
    
    do jlev=1,nLevEns_M
       clnomvar = 'UU' 
       ikey=fstlir(gd2d,nulens(jens),ini,inj,ink,idateo,cletiket,nip1_M(jlev),ip2,ip3,cltypvar,clnomvar)
       if(ikey.lt.0) then
          write(*,*) idateo,cletiket,nip1_M(jlev),ip2,ip3,cltypvar,clnomvar
          call utl_abort('readEnsembleMember: Problem with UU ENS')
       endif
       do jlat=1,nj
          do jlon=1,ni
             ensPerturbations(jlon,jlat,jlev-1+ngposituu,jens)= gd2d(jlon,jlat)*rmsknt
          enddo
       enddo
    enddo
    
    do jlev=1,nLevEns_M
       clnomvar = 'VV' 
       ikey=fstlir(gd2d,nulens(jens),ini,inj,ink,idateo,cletiket,nip1_M(jlev),ip2,ip3,cltypvar,clnomvar)
       if(ikey.lt.0) then
          write(*,*) idateo,cletiket,nip1_M(jlev),ip2,ip3,cltypvar,clnomvar
          call utl_abort('readEnsembleMember: Problem with VV ENS')
       endif
       do jlat=1,nj
          do jlon=1,ni
             ensPerturbations(jlon,jlat,jlev-1+ngpositvv,jens)= gd2d(jlon,jlat)*rmsknt
          enddo
       enddo
    enddo

!          clnomvar = 'TG' 
!          ikey = utl_fstlir(gd2d,nulens,ini,inj,ink,idateo,cletiket,-1,ip2,ip3,cltypvar,clnomvar)
!          if(ikey.lt.0)  then
!            write(*,*) idateo,cletiket,ip2,ip3,cltypvar,clnomvar
!            call utl_abort('readEnsembleMember: Problem with TG ENS')
!          else
!            do jlat=1,nj
!              do jlon=1,ni
!                ensPerturbations(jlon,jlat,ngposittg,jens)= gd2d(jlon,jlat)
!              enddo
!            enddo
!          endif

    
  end subroutine readEnsembleMember

!--------------------------------------------------------------------------
! UV_TO_PSICHI
!--------------------------------------------------------------------------
  subroutine uv_to_psichi(ensPerturbations)
    implicit none

    integer :: jens,jlev,jla
    real*8  :: dla2
    real*8  :: spectralState(nla,2,nkgdimEns)
    real*4  :: ensPerturbations(:,:,:,:)
    real*8  :: member(ni,nj,nkgdimens)
!
! Convert from U/V to PSI/CHI and spectrally filter all fields
!
    call tmg_start(2,'UV_TO_PSICHI')
    dla2   = dble(ra)*dble(ra)
    do jens=1,nens
      member(:,:,:)=dble(ensPerturbations(:,:,:,jens))
      call gst_setID(gstID_nkgdimEns)
      call gst_gdsp(spectralState,member,nlevEns_M)
      do jlev = 1, nlevEns_M
        do jla = 1, nla
          spectralState(jla,1,jlev)           = spectralState(jla,1,jlev)           * dla2*gst_getR1snp1(jla)
          spectralState(jla,2,jlev)           = spectralState(jla,2,jlev)           * dla2*gst_getR1snp1(jla)
          spectralState(jla,1,jlev+nlevEns_M) = spectralState(jla,1,jlev+nlevEns_M) * dla2*gst_getR1snp1(jla)
          spectralState(jla,2,jlev+nlevEns_M) = spectralState(jla,2,jlev+nlevEns_M) * dla2*gst_getR1snp1(jla)
        enddo
      enddo
      call gst_speree(spectralState,member)
      ensPerturbations(:,:,:,jens)=sngl(member(:,:,:))
    enddo

    call tmg_stop(2)
    write(*,*) 'finished doing u/v -> psi/chi and spectral filter...'
    call flush(6)
    
  end subroutine uv_to_psichi

!--------------------------------------------------------------------------
! REMOVEMEAN
!--------------------------------------------------------------------------
  subroutine removeMean(ensPerturbations)
    implicit none

    integer :: jens,jlev,jlat,jlon
    real*8  :: dnens,gd2d(ni,nj)
    real*4  :: ensPerturbations(:,:,:,:)

    ! remove mean and divide by sqrt(2*(NENS-1)) - extra 2 is needed?
    dnens=1.0d0/dble(nens)
!$OMP PARALLEL DO PRIVATE (JLEV,GD2D,JENS,JLAT,JLON)
      do jlev=1,nkgdimEns
        gd2d(:,:)=0.0d0
        do jens=1,nens
          do jlat=1,nj
            do jlon=1,ni
              gd2d(jlon,jlat)=gd2d(jlon,jlat)+ensPerturbations(jlon,jlat,jlev,jens)
            enddo
          enddo
        enddo
        do jlat=1,nj
          do jlon=1,ni
            gd2d(jlon,jlat)=gd2d(jlon,jlat)*dnens
          enddo
        enddo
        do jens=1,nens
          do jlat=1,nj
            do jlon=1,ni
              ensPerturbations(jlon,jlat,jlev,jens)=     &
                ensPerturbations(jlon,jlat,jlev,jens)-gd2d(jlon,jlat)
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

    write(*,*) 'finished removing the ensemble mean...'
    call flush(6)

  end subroutine removeMean


!--------------------------------------------------------------------------
! HORIZCORRELFUNCTION
!--------------------------------------------------------------------------
  subroutine horizCorrelFunction(rstddev,variableType, waveBandIndex_opt)
    implicit none

    real*8,  intent(in) :: rstddev(nkgdimEns,0:ntrunc)
    integer, intent(in) :: variableType
    integer,optional, intent(in) :: waveBandIndex_opt

    real*8  :: spectralState(nla,2,nkgdimEns)
    real*8  :: gridState(ni,nj,nkgdimEns)

    integer :: ji, jj, jk, jn, jm, ila, iref, jref
    integer :: nLevEns, nLevStart, nLevEnd, jvar, iStart, iEnd

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
             endif
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
    if ( nWaveBand /= 1 ) then
       if (.not. present(waveBandIndex_opt)) then
          write(*,*) 'horizCorrelFunction: No waveBandIndex was supplied!!!'
          call utl_abort('calbmatrix_glb')
       end if
       write(wbnum,'(I2.2)') waveBandIndex_opt
    end if

    !- 4.1 2D correlation function in fst format
    if ( nWaveBand == 1 ) then
       outfilename = "./horizCorrel.fst"
    else
       outfilename = "./horizCorrel_"//wbnum//".fst"
    end if
    call write3d(GridState,outfilename,'HORIZCORFUNC',variableType)

    !- 4.2 1D correlation function in txt format (for plotting purposes)
    iStart=iref
    iEnd=3*ni/4 ! About 10 000 km away from the center of the domain

    do jvar = 1, nvar3d
       if ( nWaveBand == 1 ) then
          outfilename = "./horizCorrel_"//trim(nomvar3d(jvar,variableType))//".txt"
       else
          outfilename = "./horizCorrel_"//trim(nomvar3d(jvar,variableType))//"_"//wbnum//".txt"
       end if
       open (unit=99,file=outfilename,action="write",status="new")

       if(vnl_varLevelFromVarName(nomvar3d(jvar,variableType)).eq.'MM') then
          nLevEns = nLevEns_M
       else
          nLevEns = nLevEns_T
       endif
       nLevStart = varLevOffset(jvar)+ 1 
       nLevEnd   = varLevOffset(jvar)+ nLevEns

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

    do jvar = 1, nvar2d
       if ( nWaveBand == 1 ) then
          outfilename = "./horizCorrel_"//trim(nomvar2d(jvar,variableType))//".txt"
       else
          outfilename = "./horizCorrel_"//trim(nomvar2d(jvar,variableType))//"_"//wbnum//".txt"
       end if
       open (unit=99,file=outfilename,action="write",status="new")
       do ji=iStart,iEnd
          write(99,'(I7,2X,F7.1,2X,F6.4)')  ji-iStart, (ji-iStart)*gridSpacingInKm, GridState(ji,jref,varLevOffset(nvar3d+1)+jvar)
       end do
       close(unit=99)
    end do

  end subroutine horizCorrelFunction

!--------------------------------------------------------------------------
! WRITE3D
!--------------------------------------------------------------------------
  subroutine write3d(gridpoint3d,filename,etiket_in,variableType)
    implicit none

    real*8, intent(in) :: gridpoint3d(:,:,:)
    character(len=*), intent(in) :: filename,etiket_in
    integer, intent(in) :: variableType

    real*8 :: dfact,zbuf(ni,nj)
    integer jlat,jlon,jlev,ierr,jvar,nLevEns
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
    do jvar = 1, nvar3d

      if(vnl_varLevelFromVarName(nomvar3d(jvar,variableType)).eq.'MM') then
        nLevEns = nLevEns_M
        nip1_l(1:nLevEns_M)=nip1_M(1:nLevEns_M)
      else
        nLevEns = nLevEns_T
        nip1_l(1:nLevEns_T)=nip1_T(1:nLevEns_T)
      endif
      dfact=1.0d0

      do jlev=1,nlevEns
        do jlat=1,nj
          do jlon=1,ni
            zbuf(jlon,jlat)=dfact*gridpoint3d(jlon,jlat,varLevOffset(jvar)+jlev)
          enddo
        enddo
        ierr = utl_fstecr(zbuf,ipak,nulstats,idateo,0,0,ni,nj,1,nip1_l(jlev),ip2,ip3,   &
                          'E',nomvar3d(jvar,variableType),etiket,'G',0,0,0,0,idatyp,.true.)
      enddo

    enddo

    ! now do 2D variables
    do jvar = 1, nvar2d
      !if(nomvar2d(jvar,variableType).eq.'P0') then
      !  dfact=1.0d0/1.0d2
      !else
        dfact=1.0d0
      !endif

      do jlat=1,nj
        do jlon=1,ni
          zbuf(jlon,jlat)=dfact*gridpoint3d(jlon,jlat,varLevOffset(nvar3d+1)+jvar)
        enddo
      enddo
      ierr = utl_fstecr(zbuf,ipak,nulstats,idateo,0,0,ni,nj,1,0,ip2,ip3,   &
                        'E',nomvar2d(jvar,variableType),etiket,'G',0,0,0,0,idatyp,.true.)

    enddo

    ierr =  fstfrm(nulstats)
    ierr =  fclos (nulstats)

    write(*,*) 'finished writing 3d array...'
    call flush(6)

  end subroutine write3d

!--------------------------------------------------------------------------
! CALCHORIZSCALE
!--------------------------------------------------------------------------
  subroutine CalcHorizScale(rstddev,variableType,waveBandIndex_opt)
    implicit none
    
    ! Based on subroutine corrlength.ftn in the "old" var code
    real(8),intent(in) :: rstddev(nkgdimEns,0:ntrunc)
    integer,intent(in) :: variableType
    integer,optional, intent(in) :: waveBandIndex_opt

    real(8) :: HorizScale(nkgdimEns)
    real(8), pointer :: PressureProfile(:)

    integer :: jk, jn, nLevEns, jvar
    real(8) :: rjn, fact, temp, a, b

    character(len=128) :: outfilename
    character(len=2) :: wbnum

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
          HorizScale(jk) = ra * sqrt(-2.0d0*a/b)
       else
          HorizScale(jk) = 0.d0
       end if
    end do

    !
    !- 2. Write the results
    !
    if ( nWaveBand /= 1 ) then
       if (.not. present(waveBandIndex_opt)) then
          write(*,*) 'CalcHorizScale: No waveBandIndex was supplied!!!'
          call utl_abort('calbmatrix_glb')
       end if
       write(wbnum,'(I2.2)') waveBandIndex_opt
    end if

    do jvar = 1, nvar3d
       write(*,*)
       write(*,*) nomvar3d(jvar,variableType)

       if ( nWaveBand == 1 ) then
          outfilename = "./horizScale_"//trim(nomvar3d(jvar,variableType))//".txt"
       else
          outfilename = "./horizScale_"//trim(nomvar3d(jvar,variableType))//"_"//wbnum//".txt"
       end if
       open (unit=99,file=outfilename,action="write",status="new")

       if(vnl_varLevelFromVarName(nomvar3d(jvar,variableType)).eq.'MM') then
          nLevEns = nLevEns_M
          PressureProfile => pressureProfile_M
       else
          nLevEns = nLevEns_T
          PressureProfile => pressureProfile_T
       endif
       do jk=1,nlevEns
          write(* ,'(I3,2X,F6.1,2X,F6.1)')  jk, PressureProfile(jk)/100.d0, HorizScale(varLevOffset(jvar)+jk)/1000.d0
          write(99,'(I3,2X,F6.1,2X,F6.1)')  jk, PressureProfile(jk)/100.d0, HorizScale(varLevOffset(jvar)+jk)/1000.d0
       end do

       close(unit=99)
    end do

    do jvar = 1, nvar2d
       write(*,*)
       write(*,*) nomvar2d(jvar,variableType)
       
       if ( nWaveBand == 1 ) then
          outfilename = "./horizScale_"//trim(nomvar2d(jvar,variableType))//".txt"
       else
          outfilename = "./horizScale_"//trim(nomvar2d(jvar,variableType))//"_"//wbnum//".txt"
       end if
       open (unit=99,file=outfilename,action="write",status="new")

       write(* ,'(I3,2X,F6.1,2X,F6.1)') 1, 1010.0, HorizScale(varLevOffset(nvar3d+1)+jvar)/1000.d0
       write(99,'(I3,2X,F6.1,2X,F6.1)') 1, 1010.0, HorizScale(varLevOffset(nvar3d+1)+jvar)/1000.d0

       close(unit=99)
    end do

  end subroutine CalcHorizScale

!--------------------------------------------------------------------------
! WRITEPOWERSPEC
!--------------------------------------------------------------------------
  subroutine writePowerSpec(powerSpec,variableType)
    implicit none

    real*8,intent(in) :: powerSpec(nkgdimEns,0:ntrunc)
    integer,intent(in) :: variableType

    integer :: jk, nLevEns, nLevStart, nLevEnd, jvar, jn

    real(8) :: waveLength 

    character(len=128) :: outfilename

    !- Write to txt files
    do jvar = 1, nvar3d

       outfilename = "./PowerSpec_"//trim(nomvar3d(jvar,variableType))//".txt"
       open (unit=99,file=outfilename,action="write",status="new")

       if(vnl_varLevelFromVarName(nomvar3d(jvar,variableType)).eq.'MM') then
          nLevEns = nLevEns_M
       else
          nLevEns = nLevEns_T
       endif
       nLevStart = varLevOffset(jvar)+ 1 
       nLevEnd   = varLevOffset(jvar)+ nLevEns

       do jn = 0, ntrunc
          if ( jn /= 0) then
             waveLength=4.d0*asin(1.d0)*ra/dble(jn)
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

    do jvar = 1, nvar2d

       outfilename = "./PowerSpec_"//trim(nomvar2d(jvar,variableType))//".txt"
       open (unit=99,file=outfilename,action="write",status="new")
       do jn = 0, ntrunc
          if ( jn /= 0) then
             waveLength=4.d0*asin(1.d0)*ra/dble(jn)
          else
             waveLength=0.d0
          end if
          write(99,'(I4,2X,F7.1,2X,e10.3)')  jn, waveLength/1000.d0, sngl(powerSpec(varLevOffset(nvar3d+1)+jvar,jn))
       end do
       close(unit=99)
    end do

  end subroutine writePowerSpec

  !--------------------------------------------------------------------------
  ! CALCLOCALCORRELATIONS
  !--------------------------------------------------------------------------
  subroutine calcLocalCorrelations(ensPerturbations,variableType,waveBandIndex_opt)
    implicit none

    real(4), intent(in) :: ensPerturbations(:,:,:,:)
    integer, intent(in) :: variableType
    integer,optional, intent(in) :: waveBandIndex_opt

    real(8), allocatable :: localHorizCorrel(:,:,:)

    real(8) :: dnens

    integer :: i, j, k, ens
    integer :: blocklength, blockpadding, nirefpoint, njrefpoint
    integer :: iref_id, jref_id, iref, jref
    integer :: imin, imax, jmin, jmax

    integer :: nulnam, ierr, fclos, fnom
    
    character(len=128) :: outfilename
    character(len=2)   :: wbnum

    NAMELIST /NAMHVCORREL_LOCAL/blocklength, blockpadding

    !
    ! To compute the local horizontal correlation for some 'reference' grid point
    ! ... we assume that the ensemble grid point mean was removed and that
    !     the ensemble values were divided by the grid point std dev.
    !

    blocklength = 100 ! Horizontal correlation will be compute blocklength x blocklength gridpoint
                      ! around each reference point
    blockpadding = 4  ! Number of grid point padding between blocks (to set correlation to 0 between each block)

    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=NAMHVCORREL_LOCAL)
    write(*,nml=NAMHVCORREL_LOCAL)
    ierr = fclos(nulnam)

    nirefpoint = ni/blocklength ! Number of reference grid point in x
    njrefpoint = nj/blocklength ! Number of reference grid point in y

    allocate(localHorizCorrel(ni,nkgdimEns,nj))

    localHorizCorrel(:,:,:)=0.0d0

    dnens = 1.0d0/dble(nens-1)

!$OMP PARALLEL DO PRIVATE (k,jref_id,iref_id,iref,jref,jmin,jmax,imin,imax,j,i,ens)
    do k = 1, nkgdimEns

       do ens = 1, nens
          do jref_id = 1, njrefpoint
             do iref_id = 1, nirefpoint
                iref = (2*iref_id-1)*blocklength/2
                jref = (2*jref_id-1)*blocklength/2
                jmin = max(jref-(blocklength-blockpadding)/2,1)
                jmax = min(jref+(blocklength-blockpadding)/2,nj)
                imin = max(iref-(blocklength-blockpadding)/2,1)
                imax = min(iref+(blocklength-blockpadding)/2,ni)
                do j = jmin, jmax
                   do i = imin, imax
                      localHorizCorrel(i,k,j)=localHorizCorrel(i,k,j) + &
                           ensPerturbations(i,j,k,ens) * ensPerturbations(iref,jref,k,ens)
                   end do
                end do
             end do
          end do
       end do

       do j = 1, nj
          do i = 1, ni
             localHorizCorrel(i,k,j) = localHorizCorrel(i,k,j)*dnens
          end do
       end do

    end do
!$OMP END PARALLEL DO

    write(*,*) 'finished computing the local horizontal correlations...'
    call flush(6)

    !
    !- 4.  Write to file
    !
    if ( nWaveBand /= 1 ) then
       if (.not. present(waveBandIndex_opt)) then
          write(*,*) 'calcLocalCorrelations: No waveBandIndex was supplied!!!'
          call utl_abort('calbmatrix_glb')
       end if
       write(wbnum,'(I2.2)') waveBandIndex_opt
    end if

    !- 4.1 2D correlation function in fst format
    if ( nWaveBand == 1 ) then
       outfilename = "./horizCorrelLocal.fst"
    else
       outfilename = "./horizCorrelLocal_"//wbnum//".fst"
    end if
    call write3d(localHorizCorrel,outfilename,'HCORREL_LOC',variableType)

    localHorizCorrel(:,:,:) = localHorizCorrel(:,:,:)**2
    call write3d(localHorizCorrel,outfilename,'HCORREL2_LOC',variableType)

    deallocate(localHorizCorrel)

  end subroutine calcLocalCorrelations

end module calcstatsglb_mod
