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

module bMatrixLatBands_mod
  ! MODULE bMatrixLatBands_mod (prefix='blb' category='2. B and R matrices')
  !
  ! :Purpose: Performs transformation from control vector to analysis increment 
  !           using the background-error covariance matrix based on homogeneous
  !           and isotropic correlations. This is the Global version. A separate 
  !           module exists for limited-area applications.
  !
  use midasMpi_mod
  use mathPhysConstants_mod
  use earthConstants_mod
  use gridStateVector_mod
  use globalSpectralTransform_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use varNameList_mod
  use utilities_mod
  use interpolation_mod
  implicit none
  save
  private

  ! public procedures
  public :: blb_setup, blb_Bsqrt, blb_BsqrtAd, blb_finalize
  public :: blb_expandToMPIglobal, blb_expandToMPIglobal_r4, blb_reduceToMPIlocal, blb_reduceToMPIlocal_r4
  public :: blb_getScaleFactor, blb_truncateCV


  logical             :: initialized = .false.
  integer             :: nj,ni
  integer             :: AnalGridID ! EZscintID
  integer             :: nlev_M,nlev_T,nkgdim
  integer             :: ntrunc,nla_mpiglobal,nla_mpilocal
  integer             :: cvDim_mpilocal,cvDim_mpiglobal
  integer             :: gstID
  integer             :: nlev_bdl
  type(struct_vco),pointer :: vco_anl

  real(8),pointer     :: rgsig(:,:)
  real(8),pointer     :: rgsiguu(:,:),rgsigvv(:,:),rgsigtt(:,:),rgsigq(:,:)
  real(8),pointer     :: rgsigps(:)
  real(8),allocatable :: tgstdbg(:,:)

  real(8),allocatable :: corns(:,:,:,:)
  real(8),allocatable :: rstddev(:,:,:)

  ! originally from common blocks and possibly from the namelist:
  integer,parameter   :: maxNumLevels=200
  real(8)             :: scaleFactor(maxNumLevels)
  real(8)             :: scaleFactorLQ(maxNumLevels)
  logical             :: scaleTG
  real(8)             :: rcscltg(1)=100000.d0
  real(8)             :: rfacthum=1.0d0
  real(8)             :: rlimsuptg=3.0d0
  logical             :: llimtg=.true.
  logical             :: TweakTG
  integer             :: nulbgst=0
  real(8)             :: rvlocpsichittps
  real(8)             :: rvloclq
  real(8)             :: rlimlv_bdl  = 85000.0d0
  character(len=4)    :: stddevMode
  integer             :: filterStddev
  real(8)             :: blendMeanStddev
  integer             :: numLatBand = 3
  real(8),allocatable :: latMask(:,:)
  logical             :: zeroTropicsCrossCorr

  ! this should come from state vector object
  integer             :: numvar3d
  integer             :: numvar2d
  integer             :: nspositVO 
  integer             :: nspositDI 
  integer             :: nspositTT 
  integer             :: nspositQ
  integer             :: nspositPS 
  integer             :: nspositTG

  real(8), pointer    :: pressureProfile_M(:),pressureProfile_T(:)

  integer             :: mymBeg,mymEnd,mymSkip,mymCount
  integer             :: mynBeg,mynEnd,mynSkip,mynCount
  integer,allocatable :: mynIndex_fromn(:)
  integer             :: maxMyNla
  integer             :: myLatBeg,myLatEnd
  integer             :: myLonBeg,myLonEnd
  integer, pointer    :: ilaList_mpiglobal(:)
  integer, pointer    :: ilaList_mpilocal(:)

  integer,external    :: get_max_rss


contains

  subroutine blb_setup(hco_in,vco_in,CVDIM_OUT, mode_opt)
    implicit none

    type(struct_hco),pointer :: hco_in
    type(struct_vco),pointer :: vco_in
    integer                  :: cvDim_out
    character(len=*), intent(in), optional :: mode_opt

    character(len=15) :: bhi_mode
    integer :: jlev, nulnam, ierr, fnom, fclos, fstouv, fstfrm
    integer :: jm, jn, mynIndex, status, latPerPE, latPerPEmax, lonPerPE, lonPerPEmax, Vcode_anl
    integer :: jlat, jLatBand, lat1, lat2, lat3
    logical :: llfound, lExists
    real(8) :: zps
    type(struct_vco),pointer :: vco_file => null()
    character(len=8) :: bFileName = './bgcov'
    
    NAMELIST /NAMBLB/ntrunc, scaleFactor, scaleFactorLQ, scaleTG, TweakTG,  &
         stddevMode, filterStddev, blendMeanStddev, zeroTropicsCrossCorr, rvlocpsichittps, rvloclq

    if ( mmpi_myid == 0 ) write(*,*) 'blb_setup: starting'
    if ( mmpi_myid == 0 ) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    ! Default values for namelist variables
    ntrunc = 108
    scaleFactor(:) = 0.0d0
    scaleFactorLQ(:) = 1.0d0
    scaleTG = .true.
    TweakTG = .false.
    stddevMode = 'GD2D'
    filterStddev = -1
    blendMeanStddev = -1.0d0
    zeroTropicsCrossCorr = .true.
    rvlocpsichittps = 6.0d0
    rvloclq = 4.0d0

    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=namblb,iostat=ierr)
    if ( ierr /= 0 ) then
      if ( mmpi_myid == 0 ) write(*,*) 'WARNING: blb_setup: Error reading namelist, ' //  &
                                      'assume it will not be used!'
      cvdim_out = 0
      return
    end if
    if ( mmpi_myid == 0 ) write(*,nml=namblb)
    ierr = fclos(nulnam)

    if ( present(mode_opt) ) then
      if ( trim(mode_opt) == 'Analysis' .or. trim(mode_opt) == 'BackgroundCheck') then
        bhi_mode = trim(mode_opt)
        if ( mmpi_myid == 0 ) write(*,*)
        if ( mmpi_myid == 0 ) write(*,*) 'blb_setup: Mode activated = ', trim(bhi_mode)
      else
        write(*,*)
        write(*,*) 'mode = ', trim(mode_opt)
        call utl_abort('blb_setup: unknown mode')
      end if
    else
      bhi_mode = 'Analysis'
      if ( mmpi_myid == 0 ) write(*,*)
      if ( mmpi_myid == 0 ) write(*,*) 'blb_setup: Analysis mode activated (by default)'
    end if

    vco_anl => vco_in
    nLev_M = vco_anl%nlev_M
    nLev_T = vco_anl%nlev_T
    if ( mmpi_myid == 0 ) write(*,*) 'blb_setup: nLev_M, nLev_T =',nLev_M, nLev_T

    do jlev = 1, max(nLev_M,nLev_T)
      if ( scaleFactor(jlev) > 0.0d0 ) then 
        scaleFactor(jlev) = sqrt(scaleFactor(jlev))
      else
        scaleFactor(jlev) = 0.0d0
      end if
    end do

    if ( sum(scaleFactor(1:max(nLev_M,nLev_T))) == 0.0d0 ) then
      if ( mmpi_myid == 0 ) write(*,*) 'blb_setup: scaleFactor=0, skipping rest of setup'
      cvdim_out = 0
      return
    end if

    do jlev = 1, max(nLev_M,nLev_T)
      if ( scaleFactorLQ(jlev) > 0.0d0 ) then 
        scaleFactorLQ(jlev) = sqrt(scaleFactorLQ(jlev))
      else
        scaleFactorLQ(jlev) = 0.0d0
      end if
    end do

    ! check if analysisgrid and covariance file have the same vertical levels
    call vco_SetupFromFile( vco_file,  & ! OUT
                            bFileName )  ! IN
    if ( .not. vco_equal(vco_anl,vco_file) ) then
      call utl_abort('bmatrixLatBands: vco from analysisgrid and cov file do not match')
    end if

    status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)
    if ( Vcode_anl /= 5002 .and. Vcode_anl /= 5005 ) then
      write(*,*) 'Vcode_anl = ',Vcode_anl
      call utl_abort('blb_setup: unknown vertical coordinate type!')
    end if

    if ( .not. (gsv_varExist(varName='TT').and.gsv_varExist(varName='UU').and.gsv_varExist(varName='VV').and. &
               gsv_varExist(varName='HU').and.gsv_varExist(varName='P0').and.gsv_varExist(varName='TG')) ) then
      call utl_abort('bmatrixLatBands: Some or all weather fields are missing. If it is desired to deactivate the weather assimilation, then all entries of the array SCALEFACTOR in the namelist NAMBLB should be set to zero.')
    end if

    if ( trim(bhi_mode) == 'BackgroundCheck' ) then
      cvDim_out = 9999 ! Dummy value > 0 to indicate to the background check (s/r ose_compute_HBHT_ensemble) 
                       ! that Bhi is used
      return
    end if

    numvar3d = 4
    numvar2d = 2

    nspositVO = 1
    nspositDI = 1*nLev_M+1
    nspositTT = 2*nLev_M+1
    nspositQ  = 2*nLev_M+1*nLev_T+1
    nspositPS = 2*nLev_M+2*nLev_T+1
    nspositTG = 2*nLev_M+2*nLev_T+2
    nkgdim = nLev_M*2 + nLev_T*2 + numvar2d
    nla_mpiglobal = (ntrunc+1)*(ntrunc+2)/2
    
    ni = hco_in%ni
    nj = hco_in%nj
    AnalGridID = hco_in%EZscintID

    gstID  = gst_setup(ni,nj,ntrunc,nkgdim)
    if ( mmpi_myid == 0 ) write(*,*) 'blb_setup: returned value of gstID    =',gstID

    call mmpi_setup_latbands(nj, latPerPE, latPerPEmax, myLatBeg, myLatEnd)
    call mmpi_setup_lonbands(ni, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd)

    call mmpi_setup_m(ntrunc,mymBeg,mymEnd,mymSkip,mymCount)
    call mmpi_setup_n(ntrunc,mynBeg,mynEnd,mynSkip,mynCount)
    allocate(mynIndex_fromn(0:ntrunc))
    mynIndex = 0
    do jn = mynBeg, mynEnd, mynSkip
      mynIndex = mynIndex + 1
      mynIndex_fromn(jn) = mynIndex
    end do

    call gst_ilaList_mpiglobal(ilaList_mpiglobal,nla_mpilocal,maxMyNla,gstID,mymBeg,mymEnd,mymSkip,mynBeg,mynEnd,mynSkip)
    call gst_ilaList_mpilocal(ilaList_mpilocal,gstID,mymBeg,mymEnd,mymSkip,mynBeg,mynEnd,mynSkip)

    ! setup latitude masks
    allocate(latMask(nj,numLatBand))
    do jLatBand = 1, numLatBand
      write(*,*) 'blb_setup: selected LATBAND = ',jlatband
      lat1=nj/4
      lat2=nj/2
      lat3=3*nj/4
      write(*,*) 'lat1,2,3=',lat1,lat2,lat3
      if ( jlatband == 1 ) then
        ! Southern extratropics
        latMask(1:lat1,jLatBand) = 1.0d0
        do jlat = lat1, lat2
          latMask(jlat,jLatBand) = sqrt(0.5d0*(1.0d0+cos(dble((jlat-lat1)*4)*MPC_PI_R8/dble(nj))))
        end do
        latMask(lat2:nj,jLatBand) = 0.0d0
      else if ( jlatband == 2 ) then
        ! Tropics
        latMask(1:lat1,jLatBand) = 0.0d0
        do jlat = lat1, lat2
          latMask(jlat,jLatBand) = sqrt(0.5d0*(1.0d0+cos(dble((lat2-jlat)*4)*MPC_PI_R8/dble(nj))))
        end do
        do jlat = lat2,lat3
          latMask(jlat,jLatBand) = sqrt(0.5d0*(1.0d0+cos(dble((jlat-lat2)*4)*MPC_PI_R8/dble(nj))))
        end do
        latMask(lat3:nj,jLatBand) = 0.0d0
      else if ( jlatband == 3 ) then
        ! Northern extratropics
        latMask(1:lat2,jLatBand) = 0.0d0
        do jlat = lat2, lat3
          latMask(jlat,jLatBand) = sqrt(0.5d0*(1.0d0+cos(dble((lat3-jlat)*4)*MPC_PI_R8/dble(nj))))
        end do
        latMask(lat3:nj,jLatBand) = 1.0d0
      end if
      write(*,*) 'latMask = ',latMask(:,jLatBand)
    end do

    ! compute mpilocal control vector size
    cvDim_mpilocal = 0
    do jLatBand = 1, numLatBand
      do jm = mymBeg, mymEnd, mymSkip
        do jn = mynBeg, mynEnd, mynSkip
          if ( jm <= jn ) then
            if ( jm == 0 ) then
              ! only real component for jm=0
              cvDim_mpilocal = cvDim_mpilocal + 1*nkgdim
            else
              ! both real and imaginary components for jm>0
              cvDim_mpilocal = cvDim_mpilocal + 2*nkgdim
            end if
          end if
        end do
      end do
    end do
    cvDim_out = cvDim_mpilocal

    ! also compute mpiglobal control vector dimension
    call rpn_comm_allreduce(cvDim_mpilocal,cvDim_mpiglobal,1,"mpi_integer","mpi_sum","GRID",ierr)

    allocate(rgsig(nj,nkgdim))
    allocate(tgstdbg(ni,nj))
    rgsig(:,:) = 0.0d0
    rgsiguu => rgsig(1:nj,nspositVO:nspositVO+nlev_M-1)
    rgsigvv => rgsig(1:nj,nspositDI:nspositDI+nlev_M-1)
    rgsigtt => rgsig(1:nj,nspositTT:nspositTT+nlev_T-1)
    rgsigq  => rgsig(1:nj,nspositQ :nspositQ +nlev_T-1)
    rgsigps => rgsig(1:nj,nspositPS)
    allocate(corns(nkgdim,nkgdim,mynCount,numLatBand))
    allocate(rstddev(nkgdim,0:ntrunc,numLatBand))
    rstddev(:,:,:) = 0.0d0

    if ( mmpi_myid == 0 ) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    zps = 101000.D0
    status = vgd_levels( vco_anl%vgrid, ip1_list=vco_anl%ip1_M, levels=pressureProfile_M, &
                         sfc_field=zps, in_log=.false.)
    status = vgd_levels( vco_anl%vgrid, ip1_list=vco_anl%ip1_T, levels=pressureProfile_T, &
                         sfc_field=zps, in_log=.false.)

    llfound = .false.
    nlev_bdl = 0
    do jlev = 1, nlev_M
      if ( .not. llfound .and. (pressureProfile_M(jlev) >= rlimlv_bdl) ) then
        nlev_bdl = jlev
        llfound = .true.
      end if
    end do

    inquire(file=bFileName,exist=lExists)
    if ( lexists ) then
      ierr = fnom(nulbgst,bFileName,'RND+OLD+R/O',0)
      if ( ierr == 0 ) then
        ierr =  fstouv(nulbgst,'RND+OLD')
      else
        call utl_abort('blb_setup:NO BACKGROUND STAT FILE!!')
      end if
    end if

    call blb_readCorns

    call blb_setCrossCorr

    call blb_setupTg

    call blb_setupCorns

    if ( stddevmode == 'GD2D' ) then
      call blb_readStd
    else if ( stddevmode == 'GD3D' ) then
      call blb_readStd3d
    else
      call utl_abort('blb_setup: unknown stddevMode')
    end if

    call blb_scalestd

    ierr = fstfrm(nulbgst)
    ierr = fclos(nulbgst)

    if ( mmpi_myid == 0 ) write(*,*) 'blb_setup: finished'

    initialized = .true.

  end subroutine blb_setup


  subroutine blb_getScaleFactor(scaleFactor_out)
    implicit none
    real(8) :: scaleFactor_out(:)
    integer :: jlev

    do jlev = 1, max(nLev_M,nLev_T)
      scaleFactor_out(jlev) = scaleFactor(jlev)
    end do

  end subroutine blb_getScaleFactor


  subroutine blb_scaleStd
    implicit none

    integer :: jlev, jlon, jlat, shift_level, Vcode_anl, status

    status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)
    if ( Vcode_anl == 5002 ) then
      shift_level = 1
    else
      shift_level = 0
    end if

    do jlev = 1, nlev_M
      do jlat = 1, nj
        rgsiguu(jlat,jlev) = scaleFactor(jlev+shift_level)*rgsiguu(jlat,jlev)
        rgsigvv(jlat,jlev) = scaleFactor(jlev+shift_level)*rgsigvv(jlat,jlev)
      end do
    end do
    do jlev = 1, nlev_T
      do jlat = 1, nj
        rgsigtt(jlat,jlev) = scaleFactor(jlev)*rgsigtt(jlat,jlev)
        rgsigq(jlat,jlev)  = scaleFactorLQ(jlev)*scaleFactor(jlev)*rgsigq(jlat,jlev)
      end do
    end do
    do jlat = 1, nj
      rgsigps(jlat)  = scaleFactor(max(nLev_M,nLev_T))*rgsigps(jlat)
    end do
    ! User has the option to not scale down the STDDEV of TG (because underestimated in Benkf)
    if ( scaleTG ) then
      do jlat = 1, nj
        do jlon = 1, ni
          tgstdbg(jlon,jlat) = scaleFactor(max(nLev_M,nLev_T))*tgstdbg(jlon,jlat)
        end do
      end do
    end if

  end subroutine blb_scaleStd


  subroutine blb_setupCorns
    implicit none

    real(8) :: eigenval(nkgdim), eigenvec(nkgdim,nkgdim), result(nkgdim,nkgdim)
    real(8) :: eigenvalsqrt(nkgdim)

    integer :: mynIndex,jk1,jk2,jk3,jLatBand
    integer :: ilwork
    integer :: iulcorvert, ikey, nsize

    real(8) :: zwork(2*4*nkgdim)
    real(8) :: ztlen,zcorr,zr,zpres1,zpres2
    real(8) :: corvert(nkgdim,nkgdim), corvert_temp(nkgdim,nkgdim)

    ! standard file variables
    integer :: ni_file,nj_file,nk_file, inpas, inbits, idatyp, ideet 
    integer :: ip1,ip2,ip3,ig1,ig2,ig3,ig4,iswa,ilength,idltf
    integer :: iubc,iextr1,iextr2,iextr3,ierr
    integer :: idateo
    character(len=2)  :: cltypvar
    character(len=1)  :: clgrtyp
    character(len=4)  :: varName
    character(len=12) :: cletiket
    integer :: fstprm,fstinf
    integer :: fnom,fstouv,fstfrm,fclos

    ! Apply vertical localization to corrns

    ! streamfunction 
    ztlen = rvlocpsichittps    ! specify length scale (in units of ln(Pressure))
    if ( ztlen > 0.0d0 ) then
      ! calculate 5'th order function (from Gaspari and Cohn)
      do jk1 = 1, nlev_M
        zpres1 = log(pressureProfile_M(jk1))
        do jk2 = 1, nlev_M
          zpres2 = log(pressureProfile_M(jk2))
          zr = abs(zpres2 - zpres1)
          zcorr = gasparicohn(ztlen,zr)
          do mynIndex = 1, mynCount
            corns(jk1,jk2,mynIndex,:) = corns(jk1,jk2,mynIndex,:)*zcorr
          end do
        end do
      end do
    end if

    ! velocity potential
    ztlen = rvlocpsichittps    ! specify length scale (in units of ln(Pressure))
    if ( ztlen > 0.0d0 ) then
      ! calculate 5'th order function (from Gaspari and Cohn)
      do jk1 = 1, nlev_M
        zpres1 = log(pressureProfile_M(jk1))
        do jk2 = 1, nlev_M
          zpres2 = log(pressureProfile_M(jk2))
          zr = abs(zpres2 - zpres1)
          zcorr = gasparicohn(ztlen,zr)
          do mynIndex = 1, mynCount
            corns(jk1+nlev_M,jk2+nlev_M,mynIndex,:) = corns(jk1+nlev_M,jk2+nlev_M,mynIndex,:)*zcorr
          end do
        end do
      end do
    end if

    ! psi-chi cross-correlations
    ztlen = rvlocpsichittps    ! specify length scale (in units of ln(Pressure))
    if ( ztlen > 0.0d0 ) then
      ! calculate 5'th order function (from Gaspari and Cohn)
      do jk1 = 1, nlev_M
        zpres1 = log(pressureProfile_M(jk1))
        do jk2 = 1, nlev_M
          zpres2 = log(pressureProfile_M(jk2))
          zr = abs(zpres2 - zpres1)
          zcorr = gasparicohn(ztlen,zr)
          do mynIndex = 1, mynCount
            corns(jk1,jk2+nlev_M,mynIndex,:) = corns(jk1,jk2+nlev_M,mynIndex,:)*zcorr
            corns(jk2+nlev_M,jk1,mynIndex,:) = corns(jk2+nlev_M,jk1,mynIndex,:)*zcorr
          end do
        end do
      end do
    end if

    ! temperature
    ztlen = rvlocpsichittps
    if ( ztlen > 0.0d0 ) then
      do jk1 = 1, nlev_T
        zpres1 = log(pressureProfile_T(jk1))
        do jk2 = 1, nlev_T
          zpres2 = log(pressureProfile_T(jk2))
          zr = abs(zpres2 - zpres1)
          zcorr = gasparicohn(ztlen,zr)
          do mynIndex = 1, mynCount
            corns(jk1+2*nlev_M,jk2+2*nlev_M,mynIndex,:)  =   &
                 corns(jk1+2*nlev_M,jk2+2*nlev_M,mynIndex,:)*zcorr
          end do
        end do
      end do
    end if

    ! temp-psi cross-correlations
    ztlen = rvlocpsichittps    ! specify length scale (in units of ln(Pressure))
    if ( ztlen > 0.0d0 ) then
      ! calculate 5'th order function (from Gaspari and Cohn)
      do jk1 = 1, nlev_M
        zpres1 = log(pressureProfile_M(jk1))
        do jk2 = 1, nlev_T
          zpres2 = log(pressureProfile_T(jk2))
          zr = abs(zpres2 - zpres1)
          zcorr = gasparicohn(ztlen,zr)
          do mynIndex = 1, mynCount
            corns(jk1,jk2+2*nlev_M,mynIndex,:) = corns(jk1,jk2+2*nlev_M,mynIndex,:)*zcorr
            corns(jk2+2*nlev_M,jk1,mynIndex,:) = corns(jk2+2*nlev_M,jk1,mynIndex,:)*zcorr
          end do
        end do
      end do
    end if

    ! temp-chi cross-correlations
    ztlen = rvlocpsichittps    ! specify length scale (in units of ln(Pressure))
    if ( ztlen > 0.0d0 ) then
      ! calculate 5'th order function (from Gaspari and Cohn)
      do jk1 = 1, nlev_M
        zpres1 = log(pressureProfile_M(jk1))
        do jk2 = 1, nlev_T
          zpres2 = log(pressureProfile_T(jk2))
          zr = abs(zpres2 - zpres1)
          zcorr = gasparicohn(ztlen,zr)
          do mynIndex = 1, mynCount
            corns(jk1+nlev_M,jk2+2*nlev_M,mynIndex,:) = corns(jk1+nlev_M,jk2+2*nlev_M,mynIndex,:)*zcorr
            corns(jk2+2*nlev_M,jk1+nlev_M,mynIndex,:) = corns(jk2+2*nlev_M,jk1+nlev_M,mynIndex,:)*zcorr
          end do
        end do
      end do
    end if

    ! cross-correlation psi-ps
    ztlen = rvlocpsichittps    ! specify length scale (in units of ln(Pressure))
    if ( ztlen > 0.0d0 ) then
      ! calculate 5'th order function (from Gaspari and Cohn)
      zpres1 = log(pressureProfile_T(nlev_T))
      do jk2 = 1, nlev_M
        zpres2 = log(pressureProfile_M(jk2))
        zr = abs(zpres2 - zpres1)
        zcorr = gasparicohn(ztlen,zr)
        do mynIndex = 1, mynCount
          corns(1+2*nlev_M+2*nlev_T,jk2,mynIndex,:)  =       &
               corns(1+2*nlev_M+2*nlev_T,jk2,mynIndex,:)*zcorr
          corns(jk2,1+2*nlev_M+2*nlev_T,mynIndex,:)  =       &
               corns(jk2,1+2*nlev_M+2*nlev_T,mynIndex,:)*zcorr
        end do
      end do
    end if

    ! cross-correlation chi-ps
    ztlen = rvlocpsichittps    ! specify length scale (in units of ln(Pressure))
    if ( ztlen > 0.0d0 ) then
      ! calculate 5'th order function (from Gaspari and Cohn)
      zpres1 = log(pressureProfile_T(nlev_T))
      do jk2 = 1, nlev_M
        zpres2 = log(pressureProfile_M(jk2))
        zr = abs(zpres2 - zpres1)
        zcorr = gasparicohn(ztlen,zr)
        do mynIndex = 1, mynCount
          corns(1+2*nlev_M+2*nlev_T,jk2+nlev_M,mynIndex,:)  =       &
               corns(1+2*nlev_M+2*nlev_T,jk2+nlev_M,mynIndex,:)*zcorr
          corns(jk2+nlev_M,1+2*nlev_M+2*nlev_T,mynIndex,:)  =       &
               corns(jk2+nlev_M,1+2*nlev_M+2*nlev_T,mynIndex,:)*zcorr
        end do
      end do
    end if

    ! cross-correlation temp-ps
    ztlen = rvlocpsichittps    ! specify length scale (in units of ln(Pressure))
    if ( ztlen > 0.0d0 ) then
      ! calculate 5'th order function (from Gaspari and Cohn)
      zpres1 = log(pressureProfile_T(nlev_T))
      do jk2 = 1, nlev_T
        zpres2 = log(pressureProfile_T(jk2))
        zr = abs(zpres2 - zpres1)
        zcorr = gasparicohn(ztlen,zr)
        do mynIndex = 1, mynCount
          corns(1+2*nlev_M+2*nlev_T,jk2+2*nlev_M,mynIndex,:)  =       &
               corns(1+2*nlev_M+2*nlev_T,jk2+2*nlev_M,mynIndex,:)*zcorr
          corns(jk2+2*nlev_M,1+2*nlev_M+2*nlev_T,mynIndex,:)  =       &
               corns(jk2+2*nlev_M,1+2*nlev_M+2*nlev_T,mynIndex,:)*zcorr
        end do
      end do
    end if

    ! humidity
    ztlen = rvloclq    ! specify length scale (in units of ln(Pressure))
    if ( ztlen > 0.0d0 ) then
      ! calculate 5'th order function (from Gaspari and Cohn)
      do jk1 = 1, nlev_T
        zpres1 = log(pressureProfile_T(jk1))
        do jk2 = 1, nlev_T
          zpres2 = log(pressureProfile_T(jk2))
          zr = abs(zpres2 - zpres1)
          zcorr = gasparicohn(ztlen,zr)
          do mynIndex = 1, mynCount
            corns(jk1+2*nlev_M+nlev_T,jk2+2*nlev_M+nlev_T,mynIndex,:)  =       &
                 corns(jk1+2*nlev_M+nlev_T,jk2+2*nlev_M+nlev_T,mynIndex,:)*zcorr
          end do
        end do
      end do
    end if

    ! compute total vertical correlations (including for balanced temperature)
    if ( .true. ) then

      iulcorvert = 0
      if ( mmpi_myid == 0 ) then
        ierr = fnom(iulcorvert,'corvert_localized.fst','RND',0)
        ierr = fstouv(iulcorvert,'RND')
      end if

      do jLatBand = 1, numLatBand

        do jk2 = 1, nkgdim
          do jk1 = 1, nkgdim
            corvert_temp(jk1,jk2) = 0.0d0
            do mynIndex = mynBeg, mynEnd, mynSkip
              corvert_temp(jk1,jk2) = corvert_temp(jk1,jk2)+((2*mynIndex+1)*corns(jk1,jk2,mynIndex_fromn(mynIndex),jLatBand))
            end do
          end do
        end do
        corvert(:,:) = 0.0d0
        nsize = nkgdim*nkgdim
        call rpn_comm_allreduce(corvert_temp,corvert,nsize,"mpi_double_precision","mpi_sum","EW",ierr)

        if ( mmpi_myid == 0 ) then
          ikey = fstinf(NULBGST,ni_file,nj_file,nk_file,-1,'CORRNS',-1,0,-1,' ','ZZ')
          ierr = fstprm(ikey,idateo,ideet,inpas,ni_file,nj_file,nk_file, inbits        &
               ,idatyp,ip1,ip2,ip3,cltypvar,varName,cletiket,clgrtyp      &
               ,ig1,ig2,ig3,ig4,iswa,ilength,idltf,iubc,iextr1,iextr2      &
               ,iextr3)

          ni_file = nkgdim
          nj_file = nkgdim
          nk_file = 1
          ip1 = jLatBand
          ip2 = ntrunc
          ip3 = 0
          varName = 'ZV'
          cletiket = 'CORVERT'
          idatyp = 5

          ierr = utl_fstecr(corvert, -inbits, iulcorvert, idateo,    &
                            ideet,inpas, ni_file, nj_file, nk_file, ip1, ip2, ip3, cltypvar,     &
                            varName,cletiket,clgrtyp,ig1, ig2, ig3, ig4, idatyp,      &
                            .true.)

        end if

      end do ! jLatBand

      if ( mmpi_myid == 0 ) then
        ierr = fstfrm(iulcorvert)
        ierr = fclos(iulcorvert)
      end if

    end if

    ! compute square-root of corns for each total wavenumber
    do jLatBand = 1, numLatBand

      do mynIndex = 1, mynCount

        do jk1 = 1, nkgdim
          do jk2 = 1, nkgdim
            eigenvec(jk2,jk1) = corns(jk2,jk1,mynIndex,jLatBand)
          end do
        end do

        ! CALCULATE EIGENVALUES AND EIGENVECTORS.
        ilwork = 4*nkgdim*2
        call dsyev('V','U',nkgdim,eigenvec,nkgdim,eigenval,zwork,ilwork,ierr)
        if ( ierr /= 0 ) then
          write(*,*) 'blb_setupcorns: non-zero value of ierr for dsyev =',ierr,' returned by dsyev for wavenumber index ',mynIndex
          call utl_abort('blb_SUCORNS')
        end if

        do jk1 = 1, nkgdim
          if ( eigenval(jk1) < 1.0d-15 ) then
            eigenvalsqrt(jk1) = 0.0d0
          else
            eigenvalsqrt(jk1) = sqrt(eigenval(jk1))
          end if
        end do
 
        ! compute E * lambda^1/2
        result(:,:) = 0.0d0
        do jk1 = 1, nkgdim
          do jk2 = 1, nkgdim
            result(jk2,jk1) = eigenvec(jk2,jk1)*eigenvalsqrt(jk1)
          end do
        end do

        ! compute (E * lambda^1/2) * E^T if new formulation
        corns(:,:,mynIndex,jLatBand) = 0.0d0
        do jk1 = 1, nkgdim
          do jk2 = 1, nkgdim
            do jk3 = 1, nkgdim
              corns(jk2,jk1,mynIndex,jLatBand) = corns(jk2,jk1,mynIndex,jLatBand) + result(jk2,jk3)*eigenvec(jk1,jk3)
            end do
          end do
        end do

      end do ! mynIndex

    end do ! jLatBand

  end subroutine blb_setupCorns


  function gaspariCohn(ztlen,zr)

    real(8)  :: gasparicohn
    real(8)  :: ztlen,zr,zlc

    zlc = ztlen/2.0d0
    if ( zr <= zlc ) then
      gasparicohn = -0.250d0*(zr/zlc)**5+0.5d0*(zr/zlc)**4             &
                  +0.625d0*(zr/zlc)**3-(5.0d0/3.0d0)*(zr/zlc)**2+1.0d0
    else if ( zr <= (2.0d0*zlc) ) then
      gasparicohn = (1.0d0/12.0d0)*(zr/zlc)**5-0.5d0*(zr/zlc)**4         &
                  +0.625d0*(zr/zlc)**3+(5.0d0/3.0d0)*(zr/zlc)**2       &
                  -5.0d0*(zr/zlc)+4.0d0-(2.0d0/3.0d0)*(zlc/zr)
    else
      gasparicohn = 0.0d0
    end if
    if ( gasparicohn < 0.0d0 ) gasparicohn = 0.0d0

  end function gaspariCohn


  subroutine blb_calcCorr(zgd,pcscl,klev)
    implicit none
    integer :: klev
    real(8) :: zgd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,klev)
    real(8) :: pcscl(klev)

    integer :: jlev, jlat, jlon
    real(8) :: zr, dlfac, dltemp, dln, dlcsurn, dlc, dlcorr

    ! parameters that define the correlation function
    integer :: ntoar = 3
    real(8) :: dlalpha = 0.2d0
    integer :: kcorrtyp = 1

    dlfac   = 1.d0/(1.d0+dlalpha)
    dln     = 1.d0*real(ntoar,8)
    dltemp  = (3.d0*(1.d0 + dlalpha))/(1.d0 + dlalpha/(dln*dln))
    dltemp  = dsqrt(dltemp)

    if ( kcorrtyp == 1 ) then
      ! Gaussian correlation
      do  jlev = 1, klev
        dlc = 1.d0/dble(pcscl(jlev))
        dlc = 0.5d0*dlc*dlc
        do  jlat = myLatBeg, myLatEnd
          zr = ec_ra * acos(gst_getRmu(jlat,gstID))
          dlcorr = dexp(-(zr**2)*dlc)
          do  jlon = myLonBeg, myLonEnd
            zgd(jlon,jlat,jlev) = dlcorr
          end do
        end do
      end do
    else if ( kcorrtyp == 2 ) then
      ! Autoregressive (SOAR) correlation
      do jlev = 1, klev
        dlc = dltemp/dble(pcscl(jlev))
        dlcsurn = dlc/dln
        do jlat = myLatBeg, myLatEnd
          zr = ec_ra * acos(gst_getRmu(jlat,gstID))
          dlcorr = (1.d0 + dlc*zr + zr*dlc*zr*dlc/3.d0)*dexp(-zr*dlc)    &
            + dlalpha*(1.d0 + dlcsurn*zr + zr*dlcsurn*zr*dlcsurn/3.d0)*dexp(-zr*dlcsurn)
          dlcorr = dlcorr*dlfac
          do jlon = myLonBeg, myLonEnd
            zgd(jlon,jlat,jlev) = dlcorr
          end do
        end do
      end do
    else
      call utl_abort('CALCCORR- Undefined correlation type')
    end if

  end subroutine blb_calcCorr


  subroutine blb_setupTG
    use timeCoord_mod
    implicit none

    logical :: llpb
    integer :: ikey, jlat, jlon, jla, ezgprm, ezqkdef
    integer :: jn, jm, ila_mpilocal, ila_mpiglobal, inlev, itggid, inmxlev, iset, nsize
    integer :: ezdefset
    integer :: ip1style,ip1kind
    integer :: koutmpg
    real(8), allocatable :: dltg(:,:), tgstdbg_tmp(:,:)
    real(8) :: cortgg(nla_mpiglobal,2)
    real(8) :: rcscltg_vec(nkgdim)
    real(8) :: zabs, zpole, dlfac
    real(8) :: zsp_mpilocal(nla_mpilocal,2,nkgdim)
    real(8) :: zgd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)
    real(8) :: zsp_mpiglobal(nla_mpiglobal,2,1)

    real(8),allocatable :: my_zsp_mpiglobal(:,:,:)

    real(4), allocatable :: TrialLandSeaMask(:,:), TrialSeaIceMask(:,:)
    real(4), allocatable :: AnalLandSeaMask(:,:), AnalSeaIceMask(:,:)

    ! standard file variables
    integer :: ni_file,nj_file,nk_file
    integer :: ip1,ip2,ip3,ig1,ig2,ig3,ig4,idatyp
    integer :: ierr,ntrials
    integer :: idateo
    integer :: fstprm,fstinf,iultg,fnom,fclos,fstouv,fstfrm
    integer :: ip1s(1), nulbgsts(1)

    integer :: TrlmNumberWanted
    integer :: fstlir, key, nultrl, ni_trial, nj_trial
    integer :: deet, npas, nbits, datyp
    integer :: swa, lng, dltf, ubc
    integer :: extra1, extra2, extra3
    
    integer :: TrialGridID

    character(len=2)  :: cltypvar
    character(len=1)  :: clgrtyp
    character(len=4)  :: varName
    character(len=12) :: cletiket

    character(len=2) :: flnum
    character(len=128) :: trialfile

    logical :: trialExists

    !
    !- 1.  Reading and processing TG standard deviations
    !

    !- 1.1 Read the Std. Dev. field from ./bgcov file
    varName = 'TG'
    idateo = -1
    inmxlev = 1
    ntrials = 1
    nulbgsts(1) = nulbgst

    call utl_getfldprm(IP1S, IP2, IP3, INLEV, CLETIKET, CLTYPVAR, ITGGID,  &
                       varName, idateo, inmxlev, nulbgsts, ip1style, ip1kind,  &
                       ntrials, koutmpg)
    ip1 = ip1s(1)

    ierr = ezgprm(itggid,CLGRTYP,NI_FILE,NJ_FILE,IG1,IG2,IG3,IG4)
    allocate(dltg(ni_file,nj_file))

    ikey = utl_fstlir(dltg,koutmpg,ni_file,nj_file,nk_file,idateo,cletiket,ip1,   &
           ip2, ip3, cltypvar, varName)

    !- 1.2 Rearrange the data according to the analysis grid (if necessary)
    if ( clgrtyp == 'G' .and. ni == ni_file .and. nj == nj_file .and. ig1 == 0  &
          .and. ig2 ==0 .and. ig3 == 0 .and.ig4 == 0 ) then
      !- 1.2.1 The std. dev. on the analysis grid 
      do jlat = 1, nj
        do jlon = 1,ni
          tgstdbg(jlon,jlat) = dltg(jlon,jlat)
        end do
      end do

    else if ( clgrtyp == 'G' .and. ni == ni_file .and. nj == nj_file .and.   &
             ig1 == 0 .and. ig2 ==1 .and. ig3 == 0 .and.ig4 == 0 ) then
       !- 1.2.2 flipped Gaussian grid no longer supported
       call utl_abort('blb_setupTg: The flipped Gaussian grid is no longer supported!')

    else
       !- 1.2.3 The std. dev. are NOT on the analysis grid. Interpolation is needed
       iset = ezdefset(AnalGridID,itggid)
       if ( TweakTG ) then
          ierr = int_ezsint(tgstdbg,dltg,interpDegree='NEAREST')
       else
          ierr = int_ezsint(tgstdbg,dltg,interpDegree='CUBIC')
       end if

    end if

    !- 1.3 Tweaking the Std Dev.

    !- 1.3.1 If specified at the top of the module, do not accept TG errors of more than value specified above
    if ( llimtg ) then
       if ( mmpi_myid == 0 ) write(*,*)
       if ( mmpi_myid == 0 ) write(*,*) 'Capping TG Std. Dev. using a max value (K) = ', rlimsuptg
       where ( tgstdbg > rlimsuptg) tgstdbg = rlimsuptg
    end if

    !- 1.3.2 Take into account the Land-Sea mask and the Sea-Ice mask of the day
    if ( TweakTG ) then

      if ( mmpi_myid == 0 ) write(*,*)
      if ( mmpi_myid == 0 ) write(*,*) 'Adjusting TG Std Dev based on LandSea and SeaIce masks'

      !- Read MG and GL in the middle of the assimilation time window
      if ( tim_nStepObs == 1 ) then
         TrlmNumberWanted = 1
      else
         TrlmNumberWanted = nint( (tim_nStepObs + 1.d0) / 2.d0)
      end if

      write(flnum,'(I2.2)') TrlmNumberWanted
      trialfile='./trlm_'//trim(flnum)
      inquire(file=trim(trialfile),exist=trialExists)

      if ( .not. trialExists ) then
        if ( mmpi_myid == 0 ) write(*,*)
        if ( mmpi_myid == 0 ) write(*,*) 'Trial file not found = ', trialfile
        if ( mmpi_myid == 0 ) write(*,*) 'Look for an ensemble of trial files '

        trialfile='./trlm_'//trim(flnum)//'_0001'
        inquire(file=trim(trialfile),exist=trialExists)
        if ( .not. trialExists ) then
           if ( mmpi_myid == 0 ) write(*,*) 'Ensemble trial file not found = ', trialfile
           call utl_abort('blb_setupTg : DID NOT FIND A TRIAL FIELD FILE')
        end if
      end if

      nultrl = 0
      ierr = fnom(nultrl,trim(trialfile),'RND+OLD+R/O',0)
      ierr = fstouv(nultrl,'RND+OLD')

      !- Determine grid size and EZSCINT ID
      idateo    = -1
      cletiket = ' '
      ip1      = -1
      ip2      = -1
      ip3      = -1
      cltypvar = ' '
      varName = 'MG'

      key = fstinf( nultrl,                                             & ! IN
                    ni_trial, nj_trial, nk_file,                            & ! OUT
                    idateo, cletiket, ip1, ip2, ip3, cltypvar, varName ) ! IN

      if ( key < 0 ) then
         write(*,*)
         write(*,*) 'blb_setupTg: Unable to find trial field = ',varName
         call utl_abort('blb_setupTg')
      end if

      ierr = fstprm( key,                                                 & ! IN
                     idateo, deet, npas, ni_trial, nj_trial, nk_file, nbits,  & ! OUT
                     datyp, ip1, ip2, ip3, cltypvar, varName, cletiket,  & ! OUT
                     clgrtyp, ig1, ig2, ig3,                              & ! OUT
                     ig4, swa, lng, dltf, ubc, extra1, extra2, extra3 )     ! OUT

      allocate(TrialLandSeaMask(ni_trial, nj_trial))
      allocate(TrialSeaIceMask(ni_trial, nj_trial))

      idateo   = -1
      cletiket = ' '
      ip1      = -1
      ip2      = -1
      ip3      = -1
      cltypvar = ' '
      varName = 'MG'
      ierr = fstlir(TrialLandSeaMask, nultrl, ni_file, nj_file, nk_file,  &
                    idateo ,cletiket, ip1, ip2, ip3, cltypvar, varName)
      if ( ierr < 0 ) then
         write(*,*)
         write(*,*) 'blb_setupTg: Unable to read trial field = ',varName
         call utl_abort('BMatrixLatBands : fstlir failed')
      end if

      if ( ni_file /= ni_trial .or. nj_file /= nj_trial ) then
          write(*,*)
          write(*,*) 'blb_setupTg: Invalid dimensions for ...'
          write(*,*) 'nomvar      =', trim(varName)
          write(*,*) 'etiket      =', trim(cletiket)
          write(*,*) 'ip1         =', ip1
          write(*,*) 'Found ni,nj =', ni_file, nj_file 
          write(*,*) 'Should be   =', ni_trial, nj_trial
          call utl_abort('blb_setupTg')
        end if

      varName = 'GL'
      ierr = fstlir(TrialSeaIceMask, nultrl, ni_file, nj_file, nk_file,  &
                    idateo ,cletiket, ip1, ip2, ip3, cltypvar, varName)
      if ( ierr < 0 ) then
         write(*,*)
         write(*,*) 'blb_setupTg: Unable to read trial field = ',varName
         call utl_abort('blb_setupTg: fstlir failed')
      end if

      if (ni_file /= ni_trial .or. nj_file /= nj_trial) then
          write(*,*)
          write(*,*) 'blb_setupTg: Invalid dimensions for ...'
          write(*,*) 'nomvar      =', trim(varName)
          write(*,*) 'etiket      =', trim(cletiket)
          write(*,*) 'ip1         =', ip1
          write(*,*) 'Found ni,nj =', ni_file, nj_file 
          write(*,*) 'Should be   =', ni_trial, nj_trial
          call utl_abort('blb_setupTg')
      end if

      TrialGridID  = ezqkdef( ni_trial, nj_trial, clgrtyp, ig1, ig2, ig3, ig4, nultrl )   ! IN
 
      ierr = fstfrm(nultrl)  
      ierr = fclos(nultrl)

      !- Interpolate to the Analysis Grid
      allocate(AnalLandSeaMask(ni, nj))
      allocate(AnalSeaIceMask(ni, nj))

      ierr = ezdefset(AnalGridID     , TrialGridID     ) ! IN,  IN

      ! Nearest-neighbor interpolation
      ierr = int_ezsint(AnalLandSeaMask, TrialLandSeaMask, interpDegree='NEAREST') ! OUT, IN
      ierr = int_ezsint(AnalSeaIceMask , TrialSeaIceMask, interpDegree='NEAREST')  ! OUT, IN

      deallocate(TrialLandSeaMask)
      deallocate(TrialSeaIceMask)

      !- Modify the input/regridded TG Std. Dev.
      do jlat = 1, nj
         do jlon = 1,ni

           if ( AnalLandSeaMask(jlon,jlat) > 0.1 ) then
             ! We take this as a land point.
             ! Force std. dev. to capping value
             tgstdbg(jlon,jlat) = rlimsuptg
           else if ( AnalSeaIceMask(jlon,jlat) > 0.2 ) then
             ! We have significant sea ice on this sea point.
             ! Force std. dev. to capping value
             tgstdbg(jlon,jlat) = rlimsuptg
           else
             ! We have an open water point. Make sure that the std. dev. is realistic
             ! 1.55 is slightly above the max value over the ocean in the legacy TG Std. Dev. field (in 2014)
             if ( tgstdbg(jlon,jlat) > 1.55d0 ) then
                tgstdbg(jlon,jlat) = 0.9d0
             end if
           end if

         end do
      end do

      !- Write the modified Std. Dev.
      if ( mmpi_myid == 0 ) then
        allocate(tgstdbg_tmp(ni,nj))
        do jlat = 1, nj
          do jlon = 1,ni
             tgstdbg_tmp(jlon,jlat) = tgstdbg(jlon,jlat)
          end do
        end do

        iultg = 0
        ierr = fnom(iultg,'tg_stddev_of_the_day.fst','RND',0)
        ierr = fstouv(iultg,'RND')

        ni_file = ni
        nj_file = nj
        nk_file = 1
        ip1 = 0
        ip2 = 0
        ip3 = 0
        idateo = 0
        cltypvar = 'E'
        varName = 'TG'
        cletiket = 'TWEAK_STDDEV'
        clgrtyp = 'G'
        ig1 = 0
        ig2 = 0
        ig3 = 0
        ig4 = 0
        idatyp = 1

        ierr = utl_fstecr(tgstdbg_tmp, -32, iultg, idateo,         &
                       0, 0, ni_file, nj_file, nk_file, ip1, ip2, ip3, cltypvar,         &
                       varName,cletiket,clgrtyp,ig1, ig2, ig3, ig4, idatyp, &
                       .true.)

        do jlat = 1, nj
          do jlon = 1,ni
             tgstdbg_tmp(jlon,jlat) = AnalLandSeaMask(jlon,jlat)
          end do
        end do
        cltypvar = 'P'
        varName = 'MG'
        cletiket = 'TRIAL2ANAL'
        ierr = utl_fstecr(tgstdbg_tmp, -32, iultg, idateo,         &
                       0, 0, ni_file, nj_file, nk_file, ip1, ip2, ip3, cltypvar,         &
                       varName,cletiket,clgrtyp,ig1, ig2, ig3, ig4, idatyp, &
                       .true.)

        do jlat = 1, nj
          do jlon = 1,ni
             tgstdbg_tmp(jlon,jlat) = AnalSeaIceMask(jlon,jlat)
          end do
        end do       
        cltypvar = 'P'
        varName = 'GL'
        cletiket = 'TRIAL2ANAL'
        ierr = utl_fstecr(tgstdbg_tmp, -32, iultg, idateo,         &
                       0, 0, ni_file, nj_file, nk_file, ip1, ip2, ip3, cltypvar,         &
                       varName,cletiket,clgrtyp,ig1, ig2, ig3, ig4, idatyp, &
                       .true.)

        ierr = fstfrm(iultg)
        ierr = fclos(iultg)

        deallocate(tgstdbg_tmp)
      end if

      deallocate(AnalLandSeaMask)
      deallocate(AnalSeaIceMask)

    end if ! if ( TweakTG )

    !
    !- 2. Compute correlations in spectral space (CORNS)
    !
    zgd(:,:,:) = 0.0d0
    zsp_mpilocal(:,:,:) = 0.0d0
    allocate(my_zsp_mpiglobal(nla_mpiglobal,2,1)) 
    my_zsp_mpiglobal(:,:,:) = 0.0d0

    do jla = 1, nla_mpiglobal
       cortgg(jla,1) = 0.0d0
       cortgg(jla,2) = 0.0d0
    end do

    ! 2.4.2  Compute correlations in physical space
    rcscltg_vec(:) = rcscltg(1)
    call blb_calccorr(zgd,rcscltg_vec,nkgdim)

    ! 2.4.3  Bring back the result in spectral space
    call gst_setID(gstID)
    call gst_reespe(zsp_mpilocal,zgd)

    ! and make the result mpiglobal
    do jm = mymBeg, mymEnd, mymSkip
      do jn = mynBeg, mynEnd, mynSkip
        if ( jm <= jn ) then
          ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
          ila_mpilocal  = ilaList_mpilocal(ila_mpiglobal)
          my_zsp_mpiglobal(ila_mpiglobal,:,1) = zsp_mpilocal(ila_mpilocal,:,1)
        end if
      end do
    end do
    nsize = 2*nla_mpiglobal
    call rpn_comm_allreduce(my_zsp_mpiglobal(:,:,1),zsp_mpiglobal(:,:,1),nsize,"mpi_double_precision","mpi_sum","GRID",ierr)
    deallocate(my_zsp_mpiglobal) 
    ! 2.4.4  Check positiveness
    llpb = .false.
    do jla = 1, ntrunc+1
      zabs = abs(zsp_mpiglobal(jla,1,1))
      llpb = llpb.or.((zsp_mpiglobal(jla,1,1) < 0.).and.(zabs > epsilon(zabs)))
    end do
    if ( llpb ) then
      call utl_abort(' AUTOCORRELATION  NEGATIVES')
    end if
    do jla = 1, ntrunc+1
      zsp_mpiglobal(jla,1,1) = abs(zsp_mpiglobal(jla,1,1))
    end do

    zpole = 0.d0
    do  jla = 1, ntrunc+1
      jn = jla-1
      zpole = zpole + zsp_mpiglobal(jla,1,1)*sqrt((2.d0*jn+1.d0)/2.d0)
    end do
    if ( zpole <= 0.d0 ) then
      call utl_abort('POLE VALUE NEGATIVE IN SETUPTG')
    end if
    do jla = 1, ntrunc+1
      zsp_mpiglobal(jla,1,1) = zsp_mpiglobal(jla,1,1)/zpole
      zsp_mpiglobal(jla,2,1) = zsp_mpiglobal(jla,2,1)/zpole
    end do

    !  2.4.5  Correlation
    do jm = 0, ntrunc
      do jn = jm, ntrunc
        jla = gst_getNIND(jm,gstID) + jn - jm
        dlfac = 0.5d0/dsqrt((2*jn+1.d0)/2.d0)
        cortgg(jla,1) = dlfac * zsp_mpiglobal(jn+1,1,1)
        cortgg(jla,2) = dlfac * zsp_mpiglobal(jn+1,1,1)
      end do
    end do

    ! 2.5. For zonal modes : set to zero the imaginary part and set the correct factor 1.0 for the real part
    do jla = 1, ntrunc + 1
      cortgg(jla,1) = 0.5d0*cortgg(jla,1)
      cortgg(jla,2) = 0.0d0
    end do

    ! 2.6. Result in corns array
    do jn = mynBeg, mynEnd, mynSkip
      ila_mpiglobal = jn + 1
      corns(nspositTG,nspositTG,mynIndex_fromn(jn),:) = 2.d0*cortgg(ila_mpiglobal,1)
    end do

    deallocate(dltg)

  end subroutine blb_setupTG


  subroutine blb_convol
    implicit none

    real(8) dlfact2,dlc,dsummed
    real(8) dtlen,zr,dlfact
    integer jn,jlat,jk,jLatBand
    real(8) zsp(0:ntrunc,nkgdim),zgr(nj,nkgdim)
    real(8) dlwti(nj),zrmu(nj)

    real(8)         :: RPORVO   = 6000.D3
    real(8)         :: RPORDI   = 6000.D3
    real(8)         :: RPORTT   = 3000.D3
    real(8)         :: RPORQ    = 3000.D3
    real(8)         :: RPORPS   = 3000.D3

    do jlat = 1, nj
       dlwti(jlat) = gst_getrwt(jlat,gstID)
       zrmu(jlat)  = gst_getrmu(jlat,gstID)
    end do

!     1.2 CONVERT THE CORRELATIONS IN SPECTRAL SPACE INTO SPECTRAL
!         COEFFICIENTS OF THE CORRELATION FUNCTION AND FUNCTION TO BE
!         SELF-CONVOLVED
    do jLatBand = 1, numLatBand

      do jn = 0, ntrunc
        dlfact = ((2.0d0*jn +1.0d0)/2.0d0)**(0.25d0)
        dlfact2= ((2.0d0*jn +1.0d0)/2.0d0)**(0.25d0)
        do jk = 1, nkgdim
          zsp(jn,jk) = rstddev(jk,jn,jLatBand)*dlfact*dlfact2
        end do
      end do

      ! Transform to physical space
      call gst_zleginv(gstID,zgr,zsp,nkgdim)

      ! Truncate in horizontal extent with Gaussian window
      do jk = 1, nkgdim
        if ( jk >= nspositVO .and. jk < nspositVO+nlev_M ) then
          dtlen = rporvo
        else if ( jk >= nspositDI .and. jk < nspositDI+nlev_M ) then
          dtlen = rpordi
        else if ( jk >= nspositTT .and. jk < nspositTT+nlev_T ) then
          dtlen = rportt
        else if ( jk >= nspositQ .and. jk < nspositQ+nlev_T ) then
          dtlen = rporq
        else if ( jk == nspositPS ) then
          dtlen = rporps
        end if

        if ( dtlen > 0.0d0 ) then
          dlc = 1.d0/dble(dtlen)
          dlc = 0.5d0*dlc*dlc
          do jlat = 1, nj
            zr = ec_ra * acos(zrmu(jlat))
            dlfact = dexp(-(zr**2)*dlc)
            zgr(jlat,jk) = dlfact*zgr(jlat,jk)
          end do
        end if

      end do

      ! Transform back to spectral space
      call gst_zlegdir(gstID,zgr,zsp,nkgdim)

      ! Convert back to correlations
      do jk = 1, nkgdim
        do jn = 0, ntrunc
          zsp(jn,jk) = zsp(jn,jk)*(2.0d0/(2.0d0*jn+1.0d0))**(0.25d0)
        end do
      end do

      ! PUT BACK INTO RSTDDEV
      do jn = 0, ntrunc
        do jk = 1, nkgdim
          rstddev(jk,jn,jLatBand) = zsp(jn,jk)
        end do
      end do
 
      ! Re-normalize to ensure correlations
      do jk = 1, nkgdim
        dsummed = 0.d0
        do jn = 0, ntrunc
          dsummed = dsummed+ dble(rstddev(jk,jn,jLatBand)**2)*sqrt(((2.d0*jn)+1.d0)/2.d0)
        end do
        if ( dsummed > 0.0d0 ) dsummed = sqrt(dsummed)
        do jn = 0, ntrunc
          if ( dsummed > 1.d-30 ) rstddev(jk,jn,jLatBand) = rstddev(jk,jn,jLatBand)/dsummed
        end do
      end do

      !     CONVERT THE SPECTRAL COEFFICIENTS OF THE CORRELATION FUNCTION
      !     .  BACK INTO CORRELATIONS OF SPECTRAL COMPONENTS
      do jn = 0, ntrunc
        dlfact = sqrt(0.5d0)*(1.0d0/((2.0d0*jn+1)/2.0d0))**0.25d0
        do jk = 1, nkgdim
          rstddev(jk,jn,jLatBand) = rstddev(jk,jn,jLatBand)*dlfact
        end do
      end do

    end do ! jLatBand

  end subroutine blb_convol


  subroutine blb_setCrossCorr
    implicit none

    integer :: jblock1, inbrblock, jblock2, jLatBand
    integer :: jk1, jk2, nlev_all(numvar3d+numvar2d), levOffset(numvar3d+numvar2d)

    inbrblock = numvar3d+numvar2d
    nlev_all(1) = nLev_M
    nlev_all(2) = nLev_M
    nlev_all(3) = nLev_T
    nlev_all(4) = nLev_T
    nlev_all(5) = 1
    nlev_all(6) = 1
    levOffset(1) = 0
    levOffset(2) = 1*nLev_M
    levOffset(3) = 2*nLev_M
    levOffset(4) = 2*nLev_M+1*nLev_T
    levOffset(5) = 2*nLev_M+2*nLev_T
    levOffset(6) = 2*nLev_M+2*nLev_T+1

    do jLatBand = 1, numLatBand

      if ( jLatBand == 2 .and. zeroTropicsCrossCorr ) then
        ! Set all cross-variable correlations to zero for tropics
        do jblock1 = 1, inbrblock
          do jblock2 = 1, jblock1
            if ( jblock1 /= jblock2 ) then
              do jk2 = 1, nlev_all(jblock2)
                do jk1 = 1, nlev_all(jblock1)
                  corns(jk1 + levOffset(jblock1),jk2 + levOffset(jblock2),:,jLatBand) = 0.0d0
                  corns(jk2 + levOffset(jblock2),jk1 + levOffset(jblock1),:,jLatBand) = 0.0d0
                end do
              end do
            end if
          end do
        end do
      else
        ! Only set cross-variable correlations with humidity and TG to zero in extra-tropics
        do jblock1 = 1, inbrblock
          do jblock2 = 1, jblock1-1
            if ( jblock1==4 .or. jblock1==6 .or. jblock2==4 .or. jblock2==6 ) then
              do jk2 = 1, nlev_all(jblock2)
                do jk1 = 1, nlev_all(jblock1)
                  corns(jk1 + levOffset(jblock1),jk2 + levOffset(jblock2),:,jLatBand) = 0.0d0
                  corns(jk2 + levOffset(jblock2),jk1 + levOffset(jblock1),:,jLatBand) = 0.0d0
                end do
              end do
            end if
          end do
        end do
      end if

    end do ! jLatBand

  end subroutine blb_setCrossCorr


  subroutine blb_readCorns
    implicit none

    integer :: jn, istdkey,icornskey
    integer :: iksdim,jcol,jrow,jLatBand
    real(8), allocatable, dimension(:) :: zstdsrc
    real(8), allocatable, dimension(:,:) :: zcornssrc

    ! standard file variables
    integer :: ni_file,nj_file,nk_file
    integer :: ip1,ip2,ip3
    integer :: idateo
    character(len=2)  :: cltypvar
    character(len=4)  :: varName
    character(len=12) :: cletiket

    iksdim = 2*nlev_M+2*nlev_T+1    ! assume 4 3d variables and 1 2d variable (TG not included)
    allocate(zcornssrc(iksdim,iksdim))
    allocate(zstdsrc(iksdim))

    do jlatBand = 1, numLatBand

      do jn = 0, ntrunc

        ! Looking for FST record parameters..

        idateo = -1
        cletiket = 'RSTDDEV'
        ip1 = jLatBand
        ip2 = jn
        ip3 = -1
        cltypvar = 'X'
        varName = 'SS'

        istdkey = utl_fstlir(ZSTDSRC,nulbgst,NI_FILE,NJ_FILE,NK_FILE,idateo,cletiket,ip1,ip2,ip3,cltypvar,varName)

        if ( istdkey < 0 ) then
          call utl_abort('READCORNS: Problem with background stat file (RSTDDEV)')
        end if

        if ( ni_file /= iksdim ) then
          call utl_abort('READCORNS: BG stat levels inconsitencies')
        end if

        do jrow = 1, iksdim
          rstddev(jrow,jn,jLatBand) = zstdsrc(jrow)
        end do

      end do

      do jn = mynBeg, mynEnd, mynSkip

        ! Looking for FST record parameters..

        idateo = -1
        cletiket = 'CORRNS'
        ip1 = jLatBand
        IP2 = jn
        ip3 = -1
        cltypvar = 'X'
        varName = 'ZZ'
        icornskey = utl_fstlir(ZCORNSSRC,nulbgst,NI_FILE,NJ_FILE,NK_FILE,idateo,cletiket,ip1,ip2,ip3,cltypvar,varName)

        if ( icornskey < 0 ) then
          call utl_abort('READCORNS: Problem with background stat file (CORRNS)')
        end if

        if ( ni_file /= iksdim .or. nj_file /= iksdim ) then
          call utl_abort('READCORNS: BG stat levels inconsitencies')
        end if

        do jcol = 1, iksdim
          do jrow = 1, iksdim
            corns(jrow,jcol,mynIndex_fromn(jn),jLatBand) = zcornssrc(jrow,jcol)
          end do
        end do

      end do

    end do

    ! Apply convolution to RSTDDEV correlations
    call blb_convol

    ! Re-build of correlation matrix: factorization of corns with convoluted RSTDDEV
    do jLatBand = 1, numLatBand
      do jn = mynBeg, mynEnd, mynSkip
        do jcol = 1, nkgdim
          do jrow = 1, nkgdim
            corns(jrow,jcol,mynIndex_fromn(jn),jLatBand) =  &
              rstddev(jrow,jn,jLatBand) * corns(jrow,jcol,mynIndex_fromn(jn),jLatBand)* rstddev(jcol,jn,jLatBand)
          end do
        end do
      end do
    end do

    deallocate(zcornssrc)
    deallocate(zstdsrc)

  end subroutine blb_readCorns


  subroutine blb_readStd
    implicit none

    integer, parameter  :: inbrvar3d=4
    integer, parameter  :: inbrvar2d=1
    integer :: jvar,jfilt,count
    integer :: ikey, jlev, jlat, jlat_file
    real(8), allocatable :: zgr(:,:)
    real(8) :: zgr_interp(nj,max(nlev_M,nlev_T))
    character(len=4) :: varName3d(inbrvar3d),varName2d(inbrvar2d)
    real(8), allocatable :: rgsig_filter(:,:)
    real(8) :: globalmean

    ! standard file variables
    integer :: ni_file,nj_file,nk_file  
    integer :: ip1,ip2,ip3
    integer :: idate,nlev_MT
    character(len=2)  :: cltypvar
    character(len=4)  :: varName
    character(len=12) :: cletiket
    integer :: fstinf

    data varName3d/'PP  ','CC  ','TT  ','LQ  '/   
    data varName2d/'P0  '/

    rgsig(:,:) = 0.0d0

!   2. Reading the data

    idate = -1
    ip1   = -1
    ip2   = -1
    ip3   = -1

    cletiket = 'STDDEV'
    cltypvar = 'E'

    ! allocate array used to read 2D stddev
    varName = varName3d(1)
    ikey = fstinf(nulbgst,ni_file,nj_file,nk_file,idate,cletiket,ip1,ip2,ip3,cltypvar,varName)
    allocate(zgr(nj_file,max(nlev_M,nlev_T)))

    do jvar = 1, inbrvar3d
      varName = varName3d(jvar)
      if ( vnl_varLevelFromVarName(varName) == 'MM' ) then
        nlev_MT = nlev_M
      else
        nlev_MT = nlev_T
      end if

      ikey = fstinf(nulbgst,ni_file,nj_file,nk_file,idate,cletiket,ip1,ip2,ip3,cltypvar,varName)

      if ( nk_file /= nlev_MT ) then
        write(*,*) 'nk_file, nlev_MT=', nk_file, nlev_MT
        call utl_abort('blb_readStd: BG stat levels inconsitencies')
      end if

      if ( ikey >= 0 ) then
        ikey = utl_fstlir(zgr(:,1:nlev_MT),nulbgst,ni_file,nj_file,nk_file,idate,cletiket,ip1,ip2,ip3,cltypvar,varName)
      else
        write(*,*) 'blb_readStd: could not read varName=',varName
        call utl_abort('blb_readStd') 
      end if

      if ( nj_file == nj ) then
        zgr_interp(:,1:nlev_MT) = zgr(:,1:nlev_MT)
      else
        do jlat = 1, nj
          jlat_file = (real(nj_file-1)/real(nj-1))*(jlat-1) + 1
          zgr_interp(jlat,1:nlev_MT) = zgr(jlat_file,1:nlev_MT)
        end do
      end if

      if ( varName == 'PP' ) then
        rgsiguu(:,:) = zgr_interp(:,1:nlev_MT)
      else if ( varName == 'UC' .or. varName == 'CC' ) then
        rgsigvv(:,:) = zgr_interp(:,1:nlev_MT)
      else if ( varName == 'TT' ) then
        rgsigtt(:,:) = zgr_interp(:,1:nlev_MT)
      else if ( varName == 'LQ' ) then
        rgsigq(:,:) = max(0.10d0,zgr_interp(:,1:nlev_MT)*rfacthum)
      end if

    end do

    nlev_MT = 1
    do jvar = 1, inbrvar2d
      varName = varName2d(jvar)

      ikey = fstinf(nulbgst,ni_file,nj_file,nk_file,idate,cletiket,ip1,ip2,ip3,cltypvar,varName)

      if ( ikey >= 0 ) then
        ikey = utl_fstlir(zgr,nulbgst,ni_file,nj_file,nk_file,idate,cletiket,ip1,ip2,ip3,cltypvar,varName)
      else
        write(*,*) 'blb_readStd: could not read varName=',varName
        call utl_abort('blb_readStd') 
      end if

      if ( nj_file == nj ) then
        zgr_interp(:,1) = zgr(:,1)
      else
        do jlat = 1, nj
          jlat_file = (real(nj_file-1)/real(nj-1))*(jlat-1) + 1
          zgr_interp(jlat,1) = zgr(jlat_file,1)
        end do
      end if

      if ( varName == 'P0' ) then
        rgsigps(:) = zgr_interp(:,1)*100.0d0
      end if

    end do

    if ( filterStddev > 0 ) then

      allocate(rgsig_filter(nj,nkgdim))
      rgsig_filter(:,:) = 0.0d0
      do jlat = 1, nj
        count = 0
        do jfilt = max(1,jlat-filterStddev), min(nj,jlat+filterStddev)
          count=count+1
          rgsig_filter(jlat,:) = rgsig_filter(jlat,:) + rgsig(jfilt,:)
        end do
        rgsig_filter(jlat,:) = rgsig_filter(jlat,:)/count
      end do
      rgsig(:,:) = rgsig_filter(:,:)
      deallocate(rgsig_filter)

    end if

    if ( blendMeanStddev > 0.0d0 ) then

      do jlev = 1, nkgdim
        globalmean = 0.0d0
        do jlat = 1, nj
          globalmean = globalmean + rgsig(jlat,jlev)
        end do
        rgsig(:,jlev) = blendMeanStddev*globalmean/dble(nj) + (1.0d0-blendMeanStddev)*rgsig(:,jlev)
      end do

    end if

    deallocate(zgr)

  end subroutine blb_readStd


  subroutine blb_readStd3d
    implicit none

    integer, parameter  :: inbrvar=5
    integer :: varIndex, jfilt, count
    integer :: ikey, ierr, levIndex, jlat
    real(8) :: zgr(ni,nj)
    character(len=4) :: varNames(inbrvar)
    real(8), allocatable :: rgsig_filter(:,:)
    real(8) :: globalmean
    character(len=4) :: varLevel

    ! standard file variables
    integer :: ni_file,nj_file,nk_file  
    integer :: ip1,ip2,ip3
    integer :: idate,nlev_MT
    character(len=2)  :: cltypvar
    character(len=4)  :: varName
    character(len=12) :: cletiket
    integer :: fstinf

    data varNames/'PP  ','CC  ','TT  ','LQ  ','P0  '/   

    if ( mmpi_myid == 0 ) write(*,*) 'blb_readStd3d: starting'
    rgsig(:,:) = 0.0d0

!   2. Reading the data

    idate = -1
    ip2   = -1
    ip3   = -1

    cletiket = 'ZM_STDDEV' ! zonal mean STDDEV from diagbmatrix
    cltypvar = ' '

    do varIndex = 1, inbrvar
      varName = varNames(varIndex)
      varLevel = vnl_varLevelFromVarname(varName)
      if ( mmpi_myid == 0 ) write(*,*) 'blb_readStd3d: reading stddev for variable ', trim(varName)

      if ( varLevel == 'MM' ) then
        nlev_MT = nlev_M
      else if ( varLevel == 'TH' ) then
        nlev_MT = nlev_T
      else if ( varLevel == 'SF' ) then
        nlev_MT = 1
      else
        call utl_abort('blb_readStd3d: unknown varLevel')
      end if

      do levIndex = 1, nlev_MT
        if ( varLevel == 'MM' ) then
          ip1 = vco_anl%ip1_M(levIndex)
        else if ( varLevel == 'TH' ) then
          ip1 = vco_anl%ip1_T(levIndex)
        else if ( varLevel == 'SF' ) then
          ip1 = -1
        end if

        ikey = fstinf(nulbgst,ni_file,nj_file,nk_file,idate,cletiket,ip1,ip2,ip3,cltypvar,varName)

        if ( nj_file /= nj ) then
          write(*,*) 'blb_readStd3d: number of latitudes not consistent between '
          write(*,*) '             stats file and analysis grid: ', nj_file, nj
          call utl_abort('blb_readStd3d')
        end if

        if ( ikey >= 0 ) then
          ierr = utl_fstlir(zgr,nulbgst,ni_file,nj_file,nk_file,idate,cletiket,ip1,ip2,ip3,cltypvar,varName)
        else
          write(*,*) 'blb_readStd3d: could not read varName, ip1=',varName,ip1
          call utl_abort('blb_readStd3d') 
        end if

        if ( varName == 'PP' ) then
          rgsiguu(:,levIndex) = zgr(1,:)
        else if ( varName == 'UC' .or. varName == 'CC' ) then
          rgsigvv(:,levIndex) = zgr(1,:)
        else if ( varName == 'TT' ) then
          rgsigtt(:,levIndex) = zgr(1,:)
        else if ( varName == 'LQ' ) then
          rgsigq(:,levIndex) = max(0.10d0, zgr(1,:)*rfacthum)
        else if ( varName == 'P0' ) then
          rgsigps(:) = zgr(1,:)*100.0d0
        end if

      end do ! levIndex

    end do ! varIndex

    if ( filterStddev > 0 ) then

      allocate( rgsig_filter(nj,nkgdim) )
      rgsig_filter(:,:) = 0.0d0
      do jlat = 1, nj
        count = 0
        do jfilt = max(1,jlat-filterStddev), min(nj,jlat+filterStddev)
          count=count+1
          rgsig_filter(jlat,:) = rgsig_filter(jlat,:) + rgsig(jfilt,:)
        end do
        rgsig_filter(jlat,:) = rgsig_filter(jlat,:)/count
      end do
      rgsig(:,:) = rgsig_filter(:,:)
      deallocate( rgsig_filter )

    end if

    if ( blendMeanStddev > 0.0d0 ) then

      do levIndex = 1, nkgdim
        globalmean = 0.0d0
        do jlat = 1, nj
          globalmean = globalmean + rgsig(jlat,levIndex)
        end do
        rgsig(:,levIndex) = blendMeanStddev*globalmean/dble(nj) + (1.0d0-blendMeanStddev)*rgsig(:,levIndex)
      end do

    end if

    if ( mmpi_myid == 0 ) write(*,*) 'blb_readStd3d: finished'

  end subroutine blb_readStd3d


  subroutine blb_truncateCV(controlVector_inout,ntruncCut)
    !
    !:Purpose: To set to zero all coefficients with total wavenumber > ntruncCut
    implicit none

    real(8), pointer :: controlVector_inout(:)
    integer          :: ntruncCut
    integer          :: jn, jm, ila_mpiglobal, ila_mpilocal, jlev, jdim

    if ( .not. initialized ) then
      if ( mmpi_myid == 0 ) write(*,*) 'blb_truncateCV: bMatrixLatBands not initialized'
      return
    end if

    if ( ntruncCut > ntrunc ) then
      write(*,*) ntruncCut, ntrunc
      call utl_abort('blb_truncateCV: ntruncCut is greater than ntrunc!')
    end if

    jdim = 0
    do jlev = 1, nkgdim
      do jm = mymBeg, mymEnd, mymSkip
        do jn = mynBeg, mynEnd, mynSkip
          if ( jm <= jn ) then
            ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
            ila_mpilocal  = ilaList_mpilocal(ila_mpiglobal)
            if ( jm == 0 ) then
              ! only real component for jm=0
              jdim = jdim + 1
              if ( jn > ntruncCut ) controlVector_inout(jdim) = 0.0d0
            else
              ! both real and imaginary components for jm>0
              jdim = jdim + 1
              if ( jn > ntruncCut ) controlVector_inout(jdim) = 0.0d0
              jdim = jdim + 1
              if ( jn > ntruncCut ) controlVector_inout(jdim) = 0.0d0
            end if
          end if
        end do
      end do
    end do

  end subroutine blb_truncateCV


  subroutine blb_bSqrt(controlVector_in,statevector)
    implicit none

    real(8)   :: controlVector_in(cvDim_mpilocal)
    type(struct_gsv) :: statevector
    real(8),allocatable :: gd_out(:,:,:)
    real(8)   :: hiControlVector(nla_mpilocal,2,nkgdim,numLatBand)

    if ( mmpi_myid == 0 ) write(*,*) 'blb_bsqrt: starting'
    if ( mmpi_myid == 0 ) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if ( .not. initialized ) then
      if ( mmpi_myid == 0 ) write(*,*) 'bMatrixLatBands not initialized'
      return
    end if

    allocate(gd_out(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim))

    call blb_cain(controlVector_in,hiControlVector)
    call blb_spa2gd(hiControlVector,gd_out)

    call copyToStatevector(statevector,gd_out)

    deallocate(gd_out)

    if ( mmpi_myid == 0 ) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if ( mmpi_myid == 0 ) write(*,*) 'blb_bsqrt: done'

  end subroutine blb_bSqrt


  subroutine blb_bSqrtAd(statevector,controlVector_out)
    implicit none

    real(8)   :: controlVector_out(cvDim_mpilocal)
    type(struct_gsv) :: statevector
    real(8), allocatable :: gd_in(:,:,:)
    real(8)   :: hiControlVector(nla_mpilocal,2,nkgdim,numLatBand)

    if ( .not. initialized ) then
      if ( mmpi_myid == 0 ) write(*,*) 'bMatrixLatBands not initialized'
      return
    end if

    if ( mmpi_myid == 0 ) write(*,*) 'blb_bsqrtad: starting'
    if ( mmpi_myid == 0 ) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    allocate(gd_in(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim))

    call copyFromStatevector(statevector,gd_in)

    call blb_spa2gdad(gd_in,hiControlVector)

    call blb_cainad(hiControlVector,controlVector_out)

    deallocate(gd_in)

    if ( mmpi_myid == 0 ) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if ( mmpi_myid == 0 ) write(*,*) 'blb_bsqrtad: done'

  end subroutine blb_bSqrtAd


  subroutine copyToStatevector(statevector,gd)
    implicit none
    type(struct_gsv) :: statevector
    real(8) :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)
    integer :: jlon, jlev, jlev2, jlat, jvar, ilev1, ilev2
    real(8), pointer :: field(:,:,:)

    do jvar = 1, vnl_numvarmax 
      if ( gsv_varExist(statevector,vnl_varNameList(jvar)) ) then
        call gsv_getField(statevector,field,vnl_varNameList(jvar))
        if ( vnl_varNameList(jvar) == 'UU  ' ) then
          ilev1 = nspositVO
        else if ( vnl_varNameList(jvar) == 'VV  ' ) then
          ilev1 = nspositDI
        else if ( vnl_varNameList(jvar) == 'TT  ' ) then
          ilev1 = nspositTT
        else if ( vnl_varNameList(jvar) == 'HU  ' ) then
          ilev1 = nspositQ
        else if ( vnl_varNameList(jvar) == 'P0  ' ) then
          ilev1 = nspositPS
        else if ( vnl_varNameList(jvar) == 'TG  ' ) then
          ilev1 = nspositTG
        else
          ! Cycle (instead of abort) to allow for non-NWP assimilation (e.g. chemical data assimilation)
!          call utl_abort('bmatrixlatbands_mod: copyToStatevector: No covariances available for variable:' // vnl_varNameList(jvar))
          cycle
        end if
        ilev2 = ilev1 - 1 + gsv_getNumLev(statevector,vnl_varLevelFromVarname(vnl_varNameList(jvar)))
        do jlev = ilev1, ilev2
          jlev2 = jlev-ilev1+1
          do jlat = myLatBeg, myLatEnd
            do jlon = myLonBeg, myLonEnd
              field(jlon,jlat,jlev2) = gd(jlon,jlat,jlev)
            end do
          end do
        end do
      end if
    end do

  end subroutine copyToStatevector


  subroutine copyFromStatevector(statevector,gd)
    implicit none
    type(struct_gsv) :: statevector
    real(8)          :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)
    integer :: jlon, jlev, jlev2, jlat, jvar, ilev1, ilev2
    real(8), pointer :: field(:,:,:)

    do jvar = 1, vnl_numvarmax 
      if ( gsv_varExist(statevector,vnl_varNameList(jvar)) ) then
        call gsv_getField(statevector,field,vnl_varNameList(jvar))
        if ( vnl_varNameList(jvar) == 'UU  ' ) then
          ilev1 = nspositVO
        else if ( vnl_varNameList(jvar) == 'VV  ' ) then
          ilev1 = nspositDI
        else if ( vnl_varNameList(jvar) == 'TT  ' ) then
          ilev1 = nspositTT
        else if ( vnl_varNameList(jvar) == 'HU  ' ) then
          ilev1 = nspositQ
        else if ( vnl_varNameList(jvar) == 'P0  ' ) then
          ilev1 = nspositPS
        else if ( vnl_varNameList(jvar) == 'TG  ' ) then
          ilev1 = nspositTG
        else
          ! Cycle (instead of abort) to allow for non-NWP assimilation (e.g. chemical data assimilation)
!          call utl_abort('bmatrixlatbands_mod: copyFromStatevector: No covariances available for variable:' // vnl_varNameList(jvar))
          cycle
        end if
        ilev2 = ilev1 - 1 + gsv_getNumLev(statevector,vnl_varLevelFromVarname(vnl_varNameList(jvar)))
        do jlev = ilev1, ilev2
          jlev2 = jlev-ilev1+1
          do jlat = myLatBeg, myLatEnd
            do jlon = myLonBeg, myLonEnd
              gd(jlon,jlat,jlev) = field(jlon,jlat,jlev2)
            end do
          end do
        end do
      end if
    end do

  end subroutine copyFromStatevector


  subroutine blb_reduceToMPILocal(cv_mpilocal,cv_mpiglobal)
    implicit none
    real(8), intent(out) :: cv_mpilocal(cvDim_mpilocal)
    real(8), intent(in)  :: cv_mpiglobal(:)

    real(8), allocatable :: cv_allMaxMpilocal(:,:)

    integer, allocatable :: cvDim_allMpilocal(:), displs(:)
    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)

    integer :: jproc,jlatBand,cvDim_maxmpilocal,ierr
    integer :: jlev,jn,jm,ila_mpiglobal,jdim_mpilocal,jdim_mpiglobal

    call rpn_comm_allreduce(cvDim_mpilocal, cvDim_maxmpilocal, &
                            1,"MPI_INTEGER","MPI_MAX","GRID",ierr)

    if ( mmpi_myid == 0 ) then
      allocate(cvDim_allMpiLocal(mmpi_nprocs))
    else
      allocate(cvDim_allMpiLocal(1))
    end if

    call rpn_comm_gather(cvDim_mpiLocal   ,1,"mpi_integer",       &
                         cvDim_allMpiLocal,1,"mpi_integer",0,"GRID",ierr)

    if ( mmpi_myid == 0 ) then
      allocate(allnBeg(mmpi_nprocs))
      allocate(allnEnd(mmpi_nprocs))
      allocate(allnSkip(mmpi_nprocs))
      allocate(allmBeg(mmpi_nprocs))
      allocate(allmEnd(mmpi_nprocs))
      allocate(allmSkip(mmpi_nprocs))
    else
      allocate(allnBeg(1))
      allocate(allnEnd(1))
      allocate(allnSkip(1))
      allocate(allmBeg(1))
      allocate(allmEnd(1))
      allocate(allmSkip(1))
    end if

    call rpn_comm_gather(mynBeg  ,1,"mpi_integer",       &
                         allnBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynEnd  ,1,"mpi_integer",       &
                         allnEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynSkip ,1,"mpi_integer",       &
                         allnSkip,1,"mpi_integer",0,"GRID",ierr)

    call rpn_comm_gather(mymBeg  ,1,"mpi_integer",       &
                         allmBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymEnd  ,1,"mpi_integer",       &
                         allmEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymSkip ,1,"mpi_integer",       &
                         allmSkip,1,"mpi_integer",0,"GRID",ierr)

    ! Prepare data to be distributed
    if ( mmpi_myid == 0 ) then

      allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mmpi_nprocs))

!$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,jlatBand,jlev,jm,jn,ila_mpiglobal,jdim_mpiglobal)
      do jproc = 0, (mmpi_nprocs-1)
        cv_allmaxmpilocal(:,jproc+1) = 0.d0
        jdim_mpilocal = 0

        do jlatBand = 1, numLatBand
          do jlev = 1, nkgdim
            do jm = allmBeg(jproc+1), allmEnd(jproc+1), allmSkip(jproc+1)
              do jn = allnBeg(jproc+1), allnEnd(jproc+1), allnSkip(jproc+1)

                if ( jm <= jn ) then
                      
                  ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
                      
                  ! figure out index into global control vector
                  if ( jm == 0 ) then
                    ! for jm=0 only real part
                    jdim_mpiglobal = ila_mpiglobal
                  else
                    ! for jm>0 both real and imaginary part
                    jdim_mpiglobal = 2*ila_mpiglobal-1 - (ntrunc+1)
                  end if
                  ! add offset for level
                  jdim_mpiglobal = jdim_mpiglobal + (jlev-1) * (ntrunc+1)*(ntrunc+1)
                  ! add offset for latBand
                  jdim_mpiglobal = jdim_mpiglobal + (jlatBand-1) * nkgdim * (ntrunc+1)*(ntrunc+1)
                      
                  ! index into local control vector computer as in cain
                  if ( jm == 0 ) then
                    ! only real component for jm=0
                    jdim_mpilocal = jdim_mpilocal + 1
                    cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal)
                  else
                    ! both real and imaginary components for jm>0
                    jdim_mpilocal = jdim_mpilocal + 1
                    cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal)
                    jdim_mpilocal = jdim_mpilocal + 1
                    cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal+1)
                  end if
                      
                  if ( jdim_mpilocal > cvDim_allMpiLocal(jproc+1) ) then
                    write(*,*)
                    write(*,*) 'ERROR: jdim_mpilocal > cvDim_allMpiLocal(jproc+1)', jdim_mpilocal, cvDim_mpilocal
                    write(*,*) '       proc, jlev, jn, jm = ',jproc, jlev, jn, jm
                    call utl_abort('blb_reduceToMPILocal')
                  end if
                  if ( jdim_mpiglobal > cvDim_mpiglobal ) then
                    write(*,*)
                    write(*,*) 'ERROR: jdim_mpiglobal > cvDim_mpiglobal', jdim_mpiglobal, cvDim_mpiglobal
                    write(*,*) '       proc, jlev, jn, jm = ',jproc, jlev, jn, jm
                    call utl_abort('blb_reduceToMPILocal')
                  end if

                end if
              end do
            end do
          end do
        end do
 
      end do ! jproc
!$OMP END PARALLEL DO

    else
      allocate(cv_allmaxmpilocal(1,1))
    end if

    !- Distribute
    allocate(displs(mmpi_nprocs))
    do jproc = 0, (mmpi_nprocs-1)
      displs(jproc+1) = jproc*cvDim_maxMpiLocal ! displacement wrt cv_allMaxMpiLocal from which
                                                ! to take the outgoing data to process jproc
    end do

    call rpn_comm_scatterv(cv_allMaxMpiLocal, cvDim_allMpiLocal, displs, "mpi_double_precision", &
                           cv_mpiLocal, cvDim_mpiLocal, "mpi_double_precision", &
                           0, "GRID", ierr)

    !- End
    deallocate(displs)
    deallocate(cv_allMaxMpiLocal)
    deallocate(cvDim_allMpiLocal)
    deallocate(allnBeg)
    deallocate(allnEnd)
    deallocate(allnSkip)
    deallocate(allmBeg)
    deallocate(allmEnd)
    deallocate(allmSkip)

  end subroutine blb_reduceToMPILocal


  subroutine blb_reduceToMPILocal_r4(cv_mpilocal,cv_mpiglobal)
    implicit none
    real(4), intent(out) :: cv_mpilocal(cvDim_mpilocal)
    real(4), intent(in)  :: cv_mpiglobal(:)

    real(4), allocatable :: cv_allMaxMpilocal(:,:)

    integer, allocatable :: cvDim_allMpilocal(:), displs(:)
    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)

    integer :: jproc,jlatBand,cvDim_maxmpilocal,ierr
    integer :: jlev,jn,jm,ila_mpiglobal,jdim_mpilocal,jdim_mpiglobal

    call rpn_comm_allreduce(cvDim_mpilocal, cvDim_maxmpilocal, &
                            1,"MPI_INTEGER","MPI_MAX","GRID",ierr)

    if ( mmpi_myid == 0 ) then
      allocate(cvDim_allMpiLocal(mmpi_nprocs))
    else
      allocate(cvDim_allMpiLocal(1))
    end if

    call rpn_comm_gather(cvDim_mpiLocal   ,1,"mpi_integer",       &
                         cvDim_allMpiLocal,1,"mpi_integer",0,"GRID",ierr)

    if ( mmpi_myid == 0 ) then
      allocate(allnBeg(mmpi_nprocs))
      allocate(allnEnd(mmpi_nprocs))
      allocate(allnSkip(mmpi_nprocs))
      allocate(allmBeg(mmpi_nprocs))
      allocate(allmEnd(mmpi_nprocs))
      allocate(allmSkip(mmpi_nprocs))
    else
      allocate(allnBeg(1))
      allocate(allnEnd(1))
      allocate(allnSkip(1))
      allocate(allmBeg(1))
      allocate(allmEnd(1))
      allocate(allmSkip(1))
    end if

    call rpn_comm_gather(mynBeg  ,1,"mpi_integer",       &
                         allnBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynEnd  ,1,"mpi_integer",       &
                         allnEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynSkip ,1,"mpi_integer",       &
                         allnSkip,1,"mpi_integer",0,"GRID",ierr)

    call rpn_comm_gather(mymBeg  ,1,"mpi_integer",       &
                         allmBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymEnd  ,1,"mpi_integer",       &
                         allmEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymSkip ,1,"mpi_integer",       &
                         allmSkip,1,"mpi_integer",0,"GRID",ierr)

    ! Prepare data to be distributed
    if (mmpi_myid == 0 ) then

      allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mmpi_nprocs))

!$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,jlatBand,jlev,jm,jn,ila_mpiglobal,jdim_mpiglobal)
      do jproc = 0, (mmpi_nprocs-1)
        cv_allmaxmpilocal(:,jproc+1) = 0.d0
        jdim_mpilocal = 0

        do jlatBand = 1, numLatBand
          do jlev = 1, nkgdim
            do jm = allmBeg(jproc+1), allmEnd(jproc+1), allmSkip(jproc+1)
              do jn = allnBeg(jproc+1), allnEnd(jproc+1), allnSkip(jproc+1)

                if ( jm <= jn ) then
                      
                  ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
                      
                  ! figure out index into global control vector
                  if ( jm == 0 ) then
                    ! for jm=0 only real part
                    jdim_mpiglobal = ila_mpiglobal
                  else
                    ! for jm>0 both real and imaginary part
                    jdim_mpiglobal = 2*ila_mpiglobal-1 - (ntrunc+1)
                  end if
                  ! add offset for level
                  jdim_mpiglobal = jdim_mpiglobal + (jlev-1) * (ntrunc+1)*(ntrunc+1)
                  ! add offset for latBand
                  jdim_mpiglobal = jdim_mpiglobal + (jlatBand-1) * nkgdim * (ntrunc+1)*(ntrunc+1)
                        
                  ! index into local control vector computer as in cain
                  if ( jm == 0 ) then
                    ! only real component for jm=0
                    jdim_mpilocal = jdim_mpilocal + 1
                    cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal)
                  else
                    ! both real and imaginary components for jm>0
                    jdim_mpilocal = jdim_mpilocal + 1
                    cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal)
                    jdim_mpilocal = jdim_mpilocal + 1
                    cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal+1)
                  end if
                      
                  if ( jdim_mpilocal > cvDim_allMpiLocal(jproc+1) ) then
                    write(*,*)
                    write(*,*) 'ERROR: jdim_mpilocal > cvDim_allMpiLocal(jproc+1)', jdim_mpilocal, cvDim_mpilocal
                    write(*,*) '       proc, jlev, jn, jm = ',jproc, jlev, jn, jm
                    call utl_abort('blb_reduceToMPILocal')
                  end if
                  if ( jdim_mpiglobal > cvDim_mpiglobal ) then
                    write(*,*)
                    write(*,*) 'ERROR: jdim_mpiglobal > cvDim_mpiglobal', jdim_mpiglobal, cvDim_mpiglobal
                    write(*,*) '       proc, jlev, jn, jm = ',jproc, jlev, jn, jm
                    call utl_abort('blb_reduceToMPILocal')
                  end if

                end if ! jm <= jn
              end do ! jn
            end do ! jm
          end do
        end do
 
      end do ! jproc
!$OMP END PARALLEL DO

    else
      allocate(cv_allmaxmpilocal(1,1))
    end if

    !- Distribute
    allocate(displs(mmpi_nprocs))
    do jproc = 0, (mmpi_nprocs-1)
      displs(jproc+1) = jproc*cvDim_maxMpiLocal ! displacement wrt cv_allMaxMpiLocal from which
                                                ! to take the outgoing data to process jproc
    end do
    call rpn_comm_scatterv(cv_allMaxMpiLocal, cvDim_allMpiLocal, displs, "mpi_real4", &
                           cv_mpiLocal, cvDim_mpiLocal, "mpi_real4", &
                           0, "GRID", ierr)

    deallocate(displs)
    deallocate(cv_allMaxMpiLocal)
    deallocate(cvDim_allMpiLocal)
    deallocate(allnBeg)
    deallocate(allnEnd)
    deallocate(allnSkip)
    deallocate(allmBeg)
    deallocate(allmEnd)
    deallocate(allmSkip)

  end subroutine blb_reduceToMPILocal_r4


  subroutine blb_expandToMPIGlobal(cv_mpilocal,cv_mpiglobal)
    implicit none
    real(8), intent(in)  :: cv_mpilocal(cvDim_mpilocal)
    real(8), intent(out) :: cv_mpiglobal(:)

    real(8), allocatable :: cv_maxmpilocal(:)
    real(8), pointer :: cv_allmaxmpilocal(:,:) => null()
    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)
    integer :: jlev, jn, jm, jproc, jlatBand, ila_mpiglobal, jdim_mpilocal, jdim_mpiglobal, ierr, cvDim_maxmpilocal

    !
    !- 1.  Gather all local control vectors onto mpi task 0
    !
    call rpn_comm_allreduce(cvDim_mpilocal,cvDim_maxmpilocal,1,"mpi_integer","mpi_max","GRID",ierr)

    allocate(cv_maxmpilocal(cvDim_maxmpilocal))

    nullify(cv_allmaxmpilocal)
    if ( mmpi_myid == 0 ) then
       allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mmpi_nprocs))
    else
       allocate(cv_allmaxmpilocal(1,1))
    end if

    cv_maxmpilocal(:) = 0.0d0
    cv_maxmpilocal(1:cvDim_mpilocal) = cv_mpilocal(1:cvDim_mpilocal)

    call rpn_comm_gather(cv_maxmpilocal,    cvDim_maxmpilocal, "mpi_double_precision",  &
                         cv_allmaxmpilocal, cvDim_maxmpilocal, "mpi_double_precision", 0, "GRID", ierr )

    deallocate(cv_maxmpilocal)

    !
    !- 2.  Reorganize gathered mpilocal control vectors into the mpiglobal control vector
    !
    if ( mmpi_myid == 0 ) then
       allocate(allnBeg(mmpi_nprocs))
       allocate(allnEnd(mmpi_nprocs))
       allocate(allnSkip(mmpi_nprocs))
       allocate(allmBeg(mmpi_nprocs))
       allocate(allmEnd(mmpi_nprocs))
       allocate(allmSkip(mmpi_nprocs))
    else
       allocate(allnBeg(1))
       allocate(allnEnd(1))
       allocate(allnSkip(1))
       allocate(allmBeg(1))
       allocate(allmEnd(1))
       allocate(allmSkip(1))
    end if

    call rpn_comm_gather(mynBeg  ,1,"mpi_integer",       &
                         allnBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynEnd  ,1,"mpi_integer",       &
                         allnEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynSkip ,1,"mpi_integer",       &
                         allnSkip,1,"mpi_integer",0,"GRID",ierr)

    call rpn_comm_gather(mymBeg  ,1,"mpi_integer",       &
                         allmBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymEnd  ,1,"mpi_integer",       &
                         allmEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymSkip ,1,"mpi_integer",       &
                         allmSkip,1,"mpi_integer",0,"GRID",ierr)

    if ( mmpi_myid == 0 ) then
      cv_mpiglobal(:) = 0.0d0

!$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,jlatBand,jlev,jm,jn,ila_mpiglobal,jdim_mpiglobal)
      do jproc = 0, (mmpi_nprocs-1)
        jdim_mpilocal = 0

        do jlatBand = 1, numLatBand
          do jlev = 1, nkgdim
            do jm = allmBeg(jproc+1), allmEnd(jproc+1), allmSkip(jproc+1)
              do jn = allnBeg(jproc+1), allnEnd(jproc+1), allnSkip(jproc+1)
                if ( jm <= jn ) then

                  ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm

                  ! figure out index into global control vector
                  if ( jm == 0 ) then
                    ! for jm=0 only real part
                    jdim_mpiglobal = ila_mpiglobal
                  else
                    ! for jm>0 both real and imaginary part
                    jdim_mpiglobal = 2*ila_mpiglobal-1 - (ntrunc+1)
                  end if
                  ! add offset for level
                  jdim_mpiglobal = jdim_mpiglobal + (jlev-1) * (ntrunc+1)*(ntrunc+1)
                  ! add offset for latBand
                  jdim_mpiglobal = jdim_mpiglobal + (jlatBand-1) * nkgdim * (ntrunc+1)*(ntrunc+1)

                  ! index into local control vector
                  if ( jm == 0 ) then
                    ! only real component for jm=0
                    jdim_mpilocal = jdim_mpilocal + 1
                    cv_mpiglobal(jdim_mpiglobal) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                  else
                    ! both real and imaginary components for jm>0
                    jdim_mpilocal = jdim_mpilocal + 1
                    cv_mpiglobal(jdim_mpiglobal) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                    jdim_mpilocal = jdim_mpilocal + 1
                    cv_mpiglobal(jdim_mpiglobal+1) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                  end if

                  if ( jdim_mpiglobal > cvDim_mpiglobal )   &
                    write(*,*) 'ERROR: jdim,cvDim,mpiglobal=',jdim_mpiglobal,cvDim_mpiglobal,jlev,jn,jm

                end if
              end do
            end do
          end do
        end do
      end do ! jproc
!$OMP END PARALLEL DO

    end if ! myid == 0 

    deallocate(allnBeg)
    deallocate(allnEnd)
    deallocate(allnSkip)
    deallocate(allmBeg)
    deallocate(allmEnd)
    deallocate(allmSkip)
    deallocate(cv_allmaxmpilocal)

  end subroutine blb_expandToMPIGlobal


  subroutine blb_expandToMPIGlobal_r4(cv_mpilocal,cv_mpiglobal)
    implicit none
    real(4), intent(in)  :: cv_mpilocal(cvDim_mpilocal)
    real(4), intent(out) :: cv_mpiglobal(:)

    real(4), allocatable :: cv_maxmpilocal(:)
    real(4), pointer :: cv_allmaxmpilocal(:,:) => null()
    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)
    integer :: jlev, jn, jm, jproc, jlatBand, ila_mpiglobal, jdim_mpilocal, jdim_mpiglobal, ierr, cvDim_maxmpilocal

    !
    !- 1.  Gather all local control vectors onto mpi task 0
    !
    call rpn_comm_allreduce(cvDim_mpilocal,cvDim_maxmpilocal,1,"mpi_integer","mpi_max","GRID",ierr)

    allocate(cv_maxmpilocal(cvDim_maxmpilocal))

    nullify(cv_allmaxmpilocal)
    if ( mmpi_myid == 0 ) then
       allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mmpi_nprocs))
    else
       allocate(cv_allmaxmpilocal(1,1))
    end if

    cv_maxmpilocal(:) = 0.0d0
    cv_maxmpilocal(1:cvDim_mpilocal) = cv_mpilocal(1:cvDim_mpilocal)

    call rpn_comm_gather(cv_maxmpilocal,    cvDim_maxmpilocal, "mpi_real4",  &
                         cv_allmaxmpilocal, cvDim_maxmpilocal, "mpi_real4", 0, "GRID", ierr )

    deallocate(cv_maxmpilocal)

    !
    !- 2.  Reorganize gathered mpilocal control vectors into the mpiglobal control vector
    !
    if ( mmpi_myid == 0 ) then
       allocate(allnBeg(mmpi_nprocs))
       allocate(allnEnd(mmpi_nprocs))
       allocate(allnSkip(mmpi_nprocs))
       allocate(allmBeg(mmpi_nprocs))
       allocate(allmEnd(mmpi_nprocs))
       allocate(allmSkip(mmpi_nprocs))
    else
       allocate(allnBeg(1))
       allocate(allnEnd(1))
       allocate(allnSkip(1))
       allocate(allmBeg(1))
       allocate(allmEnd(1))
       allocate(allmSkip(1))
    end if

    call rpn_comm_gather(mynBeg  ,1,"mpi_integer",       &
                         allnBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynEnd  ,1,"mpi_integer",       &
                         allnEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mynSkip ,1,"mpi_integer",       &
                         allnSkip,1,"mpi_integer",0,"GRID",ierr)

    call rpn_comm_gather(mymBeg  ,1,"mpi_integer",       &
                         allmBeg ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymEnd  ,1,"mpi_integer",       &
                         allmEnd ,1,"mpi_integer",0,"GRID",ierr)
    call rpn_comm_gather(mymSkip ,1,"mpi_integer",       &
                         allmSkip,1,"mpi_integer",0,"GRID",ierr)

    if ( mmpi_myid == 0 ) then
      cv_mpiglobal(:) = 0.0d0

!$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,jlatBand,jlev,jm,jn,ila_mpiglobal,jdim_mpiglobal)
      do jproc = 0, (mmpi_nprocs-1)
        jdim_mpilocal = 0

        do jlatBand = 1, numLatBand
          do jlev = 1, nkgdim
            do jm = allmBeg(jproc+1), allmEnd(jproc+1), allmSkip(jproc+1)
              do jn = allnBeg(jproc+1), allnEnd(jproc+1), allnSkip(jproc+1)
                if ( jm <= jn ) then

                  ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm

                  ! figure out index into global control vector
                  if ( jm == 0 ) then
                    ! for jm=0 only real part
                    jdim_mpiglobal = ila_mpiglobal
                  else
                    ! for jm>0 both real and imaginary part
                    jdim_mpiglobal = 2*ila_mpiglobal-1 - (ntrunc+1)
                  end if
                  ! add offset for level
                  jdim_mpiglobal = jdim_mpiglobal + (jlev-1) * (ntrunc+1)*(ntrunc+1)
                  ! add offset for latBand
                  jdim_mpiglobal = jdim_mpiglobal + (jlatBand-1) * nkgdim * (ntrunc+1)*(ntrunc+1)

                  ! index into local control vector
                  if ( jm == 0 ) then
                    ! only real component for jm=0
                    jdim_mpilocal = jdim_mpilocal + 1
                    cv_mpiglobal(jdim_mpiglobal) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                  else
                    ! both real and imaginary components for jm>0
                    jdim_mpilocal = jdim_mpilocal + 1
                    cv_mpiglobal(jdim_mpiglobal) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                    jdim_mpilocal = jdim_mpilocal + 1
                    cv_mpiglobal(jdim_mpiglobal+1) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                  end if

                  if ( jdim_mpiglobal > cvDim_mpiglobal )   &
                    write(*,*) 'ERROR: jdim,cvDim,mpiglobal=',jdim_mpiglobal,cvDim_mpiglobal,jlev,jn,jm

                end if ! jm <= jn
              end do ! jn
            end do ! jm
          end do ! jlev
        end do ! jlatBand
      end do ! jproc
!$OMP END PARALLEL DO

    end if ! myid == 0 

    deallocate(allnBeg)
    deallocate(allnEnd)
    deallocate(allnSkip)
    deallocate(allmBeg)
    deallocate(allmEnd)
    deallocate(allmSkip)
    deallocate(cv_allmaxmpilocal)

  end subroutine blb_expandToMPIGlobal_r4


  subroutine blb_cain(controlVector_in,hiControlVector_out)
    implicit none

    real(8) :: controlVector_in(cvDim_mpilocal)
    real(8) :: hiControlVector_out(nla_mpilocal,2,nkgdim,numLatBand)

    integer :: jdim, jlev, jm, jn, jLatBand, ila_mpilocal, ila_mpiglobal

    jdim = 0
    hiControlVector_out(:,:,:,:) = 0.0d0
    do jLatBand = 1, numLatBand
      do jlev = 1, nkgdim
        do jm = mymBeg, mymEnd, mymSkip
          do jn = mynBeg, mynEnd, mynSkip
            if ( jm <= jn ) then
              ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
              ila_mpilocal  = ilaList_mpilocal(ila_mpiglobal)
              if ( jm == 0 ) then
                ! only real component for jm=0
                jdim = jdim + 1
                hiControlVector_out(ila_mpilocal,1,jlev,jLatBand) = controlVector_in(jdim)
              else
                ! both real and imaginary components for jm>0
                jdim = jdim + 1
                hiControlVector_out(ila_mpilocal,1,jlev,jLatBand) = controlVector_in(jdim)
                jdim = jdim + 1
                hiControlVector_out(ila_mpilocal,2,jlev,jLatBand) = controlVector_in(jdim)
              end if
            end if
          end do
        end do
      end do
    end do

  end subroutine blb_cain


  subroutine blb_cainAd(hiControlVector_in,controlVector_out)
    implicit none

    real(8) :: controlVector_out(cvDim_mpilocal)
    real(8) :: hiControlVector_in(nla_mpilocal,2,nkgdim,numLatBand)

    integer :: jdim, jlev, jm, jn, jLatBand, ila_mpilocal, ila_mpiglobal

    jdim = 0
    do jLatBand = 1, numLatBand
      do jlev = 1, nkgdim
        do jm = mymBeg, mymEnd, mymSkip
          do jn = mynBeg, mynEnd, mynSkip
            if ( jm <= jn ) then
              ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
              ila_mpilocal  = ilaList_mpilocal(ila_mpiglobal)
              if ( jm == 0 ) then
                ! only real component for jm=0
                jdim = jdim + 1
                controlVector_out(jdim) = controlVector_out(jdim) + hiControlVector_in(ila_mpilocal,1,jlev,jLatBand)
              else
                ! both real and imaginary components for jm>0
                jdim = jdim + 1
                controlVector_out(jdim) = controlVector_out(jdim) + hiControlVector_in(ila_mpilocal,1,jlev,jLatBand)*2.0d0
                jdim = jdim + 1
                controlVector_out(jdim) = controlVector_out(jdim) + hiControlVector_in(ila_mpilocal,2,jlev,jLatBand)*2.0d0
              end if
            end if
          end do
        end do
      end do
    end do

  end subroutine blb_cainAd


  subroutine blb_spa2gd(hiControlVector_in,gd_out)
    implicit none

    real(8) :: hiControlVector_in(nla_mpilocal,2,nkgdim,numLatBand)
    real(8) :: gd_out(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)

    integer :: jlev, jlon, jlat, jla_mpilocal, jLatBand, jn, jm
    integer :: ila_mpiglobal, ila_mpilocal, icount
    real(8) :: sq2, dl1sa2, dla2
    real(8) :: sp(nla_mpilocal,2,nkgdim)
    real(8), allocatable :: zsp(:,:,:), zsp2(:,:,:)
    real(8) :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)

    sq2 = sqrt(2.0d0)
    allocate(zsp(nkgdim,2,mymCount))
    allocate(zsp2(nkgdim,2,mymCount))

!$OMP PARALLEL DO PRIVATE(JLEV)
    do jlev = 1, nkgdim
      gd_out(:,:,jlev) = 0.0d0
    end do
!$OMP END PARALLEL DO

    do jLatBand = 1, numLatBand

      sp(:,:,:) = 0.0d0
!$OMP PARALLEL DO PRIVATE(jn,jm,jlev,ila_mpiglobal,ila_mpilocal,zsp2,zsp,icount)
      do jn = mynBeg, mynEnd, mynSkip

        icount = 0
        do jm = mymBeg, mymEnd, mymSkip
          if ( jm <= jn ) then
            icount = icount+1
            ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
            ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
            do jlev = 1, nkgdim
              zsp(jlev,1,icount) = hiControlVector_in(ila_mpilocal,1,jlev,jLatBand)
              zsp(jlev,2,icount) = hiControlVector_in(ila_mpilocal,2,jlev,jLatBand)
            end do
          end if
        end do

        if ( icount > 0 ) then

          CALL DGEMM('N','N',nkgdim,2*icount,nkgdim,1.0d0,corns(1,1,mynIndex_fromn(jn),jLatBand),nkgdim,zsp(1,1,1),nkgdim,0.0d0,zsp2(1,1,1),nkgdim)

          icount = 0
          do jm = mymBeg, mymEnd, mymSkip
            if ( jm <= jn ) then
              icount = icount+1
              ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
              ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
              do jlev = 1, nkgdim
                sp(ila_mpilocal,1,jlev) = zsp2(jlev,1,icount)
                sp(ila_mpilocal,2,jlev) = zsp2(jlev,2,icount)
              end do
            end if
          end do

        end if

        ! make adjustments for jm=0
        if ( mymBeg == 0 ) then
          ila_mpiglobal = gst_getNind(0,gstID) + jn
          ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
          do jlev = 1, nkgdim
            sp(ila_mpilocal,1,jlev) = sp(ila_mpilocal,1,jlev)*sq2
            sp(ila_mpilocal,2,jlev) = 0.0d0
          end do
        end if

      end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(JLEV)
      do jlev = 1, nkgdim
        gd(:,:,jlev) = 0.0d0
      end do
!$OMP END PARALLEL DO
      call gst_setID(gstID)
      call gst_speree(sp,gd)

!$OMP PARALLEL DO PRIVATE(jlat,jlev,jlon)
      do jlev = 1, nkgdim
        if ( jlev == nspositTG ) then
          do jlat = myLatBeg, myLatEnd
            do jlon = myLonBeg, myLonEnd
              gd(jlon,jlat,jlev) = gd(jlon,jlat,jlev)*tgstdbg(jlon,jlat)*latMask(jlat,jLatBand)
            end do
          end do
        else
          do jlat = myLatBeg, myLatEnd
            do jlon = myLonBeg, myLonEnd
              gd(jlon,jlat,jlev) = gd(jlon,jlat,jlev)*rgsig(jlat,jlev)*latMask(jlat,jLatBand)
            end do
          end do
        end if

        do jlat = myLatBeg, myLatEnd
          do jlon = myLonBeg, myLonEnd
            gd_out(jlon,jlat,jlev) = gd_out(jlon,jlat,jlev) + gd(jlon,jlat,jlev)
          end do
        end do
      end do
!$OMP END PARALLEL DO

    end do ! jLatBand

    deallocate(zsp)
    deallocate(zsp2)

    ! convert final result from PSI/CHI to U/V
    call gst_setID(gstID)
    call gst_reespe(sp,gd_out)

    dla2   = ec_ra * ec_ra
    dl1sa2 = 1.d0/dla2
!$OMP PARALLEL DO PRIVATE(JLEV,JLA_MPILOCAL,ILA_MPIGLOBAL)
    do jlev = 1, nlev_M
      do jla_mpilocal = 1, nla_mpilocal
        ila_mpiglobal = ilaList_mpiglobal(jla_mpilocal)
        sp(jla_mpilocal,1,nspositVO+jlev-1) = sp(jla_mpilocal,1,nspositVO+jlev-1)*dl1sa2*gst_getrnnp1(ila_mpiglobal,gstID)
        sp(jla_mpilocal,2,nspositVO+jlev-1) = sp(jla_mpilocal,2,nspositVO+jlev-1)*dl1sa2*gst_getrnnp1(ila_mpiglobal,gstID)
        sp(jla_mpilocal,1,nspositDI+jlev-1) = sp(jla_mpilocal,1,nspositDI+jlev-1)*dl1sa2*gst_getrnnp1(ila_mpiglobal,gstID)
        sp(jla_mpilocal,2,nspositDI+jlev-1) = sp(jla_mpilocal,2,nspositDI+jlev-1)*dl1sa2*gst_getrnnp1(ila_mpiglobal,gstID)
      end do
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(JLEV)
    do jlev = 1, 2*nlev_M
      gd_out(:,:,jlev) = 0.0d0
    end do
!$OMP END PARALLEL DO
    call gst_setID(gstID)
    call gst_spgd(sp,gd_out,nLev_M)

  end subroutine blb_spa2gd


  subroutine blb_spa2gdad(gd_in,hiControlVector_out)
    implicit none

    real(8) :: hiControlVector_out(nla_mpilocal,2,nkgdim,numLatBand)
    real(8) :: gd_in(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)

    integer :: jn, jm, ila_mpilocal, ila_mpiglobal, icount, jLatBand
    integer :: jlev, jlon, jlat, jla_mpilocal
    real(8) :: sq2, dl1sa2, dla2
    real(8) :: sp(nla_mpilocal,2,nkgdim)
    real(8) ,allocatable :: zsp(:,:,:), zsp2(:,:,:)
    real(8) :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)
    real(8) :: gd_in2(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)

    allocate(zsp(nkgdim,2,mymCount))
    allocate(zsp2(nkgdim,2,mymCount))

!$OMP PARALLEL DO PRIVATE(JLEV)
    do jlev = 1, nkgdim
      gd_in2(:,:,jlev) = gd_in(:,:,jlev)
      sp(:,:,jlev) = 0.0d0
    end do
!$OMP END PARALLEL DO

    ! adjoint of convert final result from PSI/CHI to U/V

    call gst_setID(gstID)
    call gst_spgda(sp,gd_in2,nLev_M)

    dla2   = ec_ra * ec_ra
    dl1sa2 = 1.d0/dla2
!$OMP PARALLEL DO PRIVATE(JLEV,JLA_MPILOCAL,ILA_MPIGLOBAL)
    do jlev = 1, nlev_M
      do jla_mpilocal = 1, nla_mpilocal
        ila_mpiglobal = ilaList_mpiglobal(jla_mpilocal)
        sp(jla_mpilocal,:,nspositVO+jlev-1) = sp(jla_mpilocal,:,nspositVO+jlev-1)*dl1sa2*gst_getrnnp1(ila_mpiglobal,gstID)
        sp(jla_mpilocal,:,nspositDI+jlev-1) = sp(jla_mpilocal,:,nspositDI+jlev-1)*dl1sa2*gst_getrnnp1(ila_mpiglobal,gstID)
      end do
    end do
!$OMP END PARALLEL DO

    call gst_setID(gstID)
    call gst_speree(sp,gd_in2)

    do jLatBand = 1, numLatBand

!$OMP PARALLEL DO PRIVATE(JLAT,JLEV,JLON)
      do jlev = 1, nkgdim
        do jlat = myLatBeg, myLatEnd
          do jlon = myLonBeg, myLonEnd
            gd(jlon,jlat,jlev) = gd_in2(jlon,jlat,jlev)
          end do
        end do

        if ( jlev == nspositTG ) then
          do jlat = myLatBeg, myLatEnd
            do jlon = myLonBeg, myLonEnd
              gd(jlon,jlat,jlev) = gd(jlon,jlat,jlev)*tgstdbg(jlon,jlat)*latMask(jlat,jLatBand)
            end do
          end do
        else
          do jlat = myLatBeg, myLatEnd
            do jlon = myLonBeg, myLonEnd
              gd(jlon,jlat,jlev) = gd(jlon,jlat,jlev)*rgsig(jlat,jlev)*latMask(jlat,jLatBand)
            end do
          end do
        end if
      end do
!$OMP END PARALLEL DO 

      call gst_setID(gstID)
      call gst_reespe(sp,gd)

      hiControlVector_out(:,:,:,jLatBand) = 0.0d0
      sq2 = sqrt(2.0d0)
!$OMP PARALLEL DO PRIVATE(JN,JM,JLEV,ILA_MPILOCAL,ILA_MPIGLOBAL,zsp,zsp2,icount)
      do jn = mynBeg, mynEnd, mynSkip

        icount = 0
        do jm = mymBeg, mymEnd, mymSkip
          if ( jm <= jn ) then
            icount = icount+1
            ila_mpiglobal = gst_getNind(jm,gstID) + jn - jm
            ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
            do jlev = 1, nkgdim
              zsp2(jlev,1,icount) = sp(ila_mpilocal,1,jlev)
              zsp2(jlev,2,icount) = sp(ila_mpilocal,2,jlev)
            end do
          end if
        end do

        if ( icount > 0 ) then

          call dgemm('T', 'N', nkgdim, 2*icount, nkgdim, 1.0d0,  &
              corns(1,1,mynIndex_fromn(jn),jLatBand), nkgdim,  &
              zsp2(1,1,1), nkgdim, 0.0d0, zsp(1,1,1), nkgdim)

          icount = 0
          do jm = mymBeg, jn, mymSkip
            icount=icount+1
            ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
            ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
            do jlev = 1, nkgdim
              hiControlVector_out(ila_mpilocal,1,jlev,jLatBand) = zsp(jlev,1,icount)
              hiControlVector_out(ila_mpilocal,2,jlev,jLatBand) = zsp(jlev,2,icount)
            end do
          end do

        end if

        ! make adjustments for jm=0
        if ( mymBeg == 0 ) then
          ila_mpiglobal = gst_getNIND(0,gstID) + jn
          ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
          do jlev = 1, nkgdim
            hiControlVector_out(ila_mpilocal,1,jlev,jLatBand) = hiControlVector_out(ila_mpilocal,1,jlev,jLatBand)*sq2
            hiControlVector_out(ila_mpilocal,2,jlev,jLatBand) = hiControlVector_out(ila_mpilocal,2,jlev,jLatBand)*sq2
          end do
        end if

      end do
!$OMP END PARALLEL DO

    end do ! jLatBand

    deallocate(zsp)
    deallocate(zsp2)

  end subroutine blb_spa2gdad


  subroutine blb_Finalize()
    implicit none

    if ( initialized ) then
      deallocate(pressureProfile_M)
      deallocate(pressureProfile_T)
      deallocate(rgsig)
      deallocate(tgstdbg)
      deallocate(corns)
      deallocate(rstddev)
    end if

  end subroutine blb_Finalize


end module bMatrixLatBands_mod
