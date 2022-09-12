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

MODULE BmatrixHI_mod
  ! MODULE BmatrixHI_mod (prefix='bhi' category='2. B and R matrices')
  !
  ! :Purpose: Performs transformation from control vector to analysis increment 
  !           using the background-error covariance matrix based on homogeneous
  !           and isotropic correlations. This is the Global version. A separate 
  !           module exists for limited-area applications.
  !
  use midasMpi_mod
  use earthConstants_mod
  use gridStateVector_mod
  use globalSpectralTransform_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use varNameList_mod
  use utilities_mod
  use gridVariableTransforms_mod
  implicit none
  save
  private

  ! public procedures
  public :: bhi_Setup,bhi_BSqrt,bhi_BSqrtAd,bhi_Finalize,bhi_expandToMPIglobal,bhi_expandToMPIglobal_r4,bhi_reduceToMPIlocal,bhi_reduceToMPIlocal_r4
  public :: bhi_getScaleFactor,bhi_truncateCV


  logical             :: initialized = .false.
  integer             :: nj_l,ni_l
  integer             :: AnalGridID ! EZscintID
  integer             :: nlev_M,nlev_T,nlev_T_even,nkgdim,nkgdim2,nkgdimSqrt
  integer             :: ntrunc,nla_mpiglobal,nla_mpilocal
  integer             :: cvDim_mpilocal,cvDim_mpiglobal
  logical             :: squareSqrt
  integer             :: gstID, gstID2
  integer             :: nlev_bdl
  type(struct_vco),pointer :: vco_anl

  real(8),allocatable :: tantheta(:,:)
  real(8),allocatable :: PtoT(:,:,:)

  real(8),pointer     :: rgsig(:,:)
  real(8),pointer     :: rgsiguu(:,:),rgsigvv(:,:),rgsigtt(:,:),rgsigtb(:,:),rgsigq(:,:)
  real(8),pointer     :: rgsigps(:),rgsigpsb(:)
  real(8),allocatable :: tgstdbg(:,:)

  real(8),allocatable :: corns(:,:,:)
  real(8),allocatable :: rstddev(:,:)

  ! originally from common blocks and possibly from the namelist:
  integer,parameter   :: maxNumLevels=200
  real(8)             :: scaleFactor(maxNumLevels)
  real(8)             :: scaleFactorLQ(maxNumLevels)
  real(8)             :: scaleFactorCC(maxNumLevels)
  logical             :: scaleTG
  real(8)             :: rcscltg(1)=100000.d0
  real(8)             :: rlimsuptg=3.0d0
  logical             :: llimtg=.true.
  logical             :: TweakTG
  logical             :: ReadWrite_sqrt
  integer             :: nulbgst=0
  integer             :: nLevPtoT
  real(8)             :: rvlocbalt   = 6.0d0
  real(8)             :: rvlocpsi    = 6.0d0
  real(8)             :: rvlocchi    = 6.0d0
  real(8)             :: rvlocpsitt  = 6.0d0
  real(8)             :: rvlocunbalt = 4.0d0
  real(8)             :: rvloclq     = 4.0d0
  real(8)             :: rlimlv_bdl  = 85000.0d0
  integer             :: numModeZero  ! number of eigenmodes to set to zero
  character(len=4)    :: stddevMode

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
  integer             :: maxMyNla
  integer             :: myLatBeg,myLatEnd
  integer             :: myLonBeg,myLonEnd
  integer, pointer    :: ilaList_mpiglobal(:)
  integer, pointer    :: ilaList_mpilocal(:)

  integer,external    :: get_max_rss


CONTAINS

  SUBROUTINE BHI_setup(hco_in,vco_in,CVDIM_OUT, mode_opt)
    implicit none

    type(struct_hco),pointer :: hco_in
    type(struct_vco),pointer :: vco_in
    integer                  :: cvDim_out
    character(len=*), intent(in), optional :: mode_opt

    character(len=15) :: bhi_mode

    integer :: jlev, nulnam, ierr, fnom, fclos, fstouv, fstfrm
    integer :: jm, jn, status, latPerPE, lonPerPE, latPerPEmax, lonPerPEmax, Vcode_anl
    logical :: llfound, lExists
    real(8) :: zps
    type(struct_vco),pointer :: vco_file => null()
    character(len=8) :: bFileName = './bgcov'

    NAMELIST /NAMBHI/ntrunc,scaleFactor,scaleFactorLQ,scaleFactorCC,scaleTG,numModeZero,squareSqrt,TweakTG,ReadWrite_sqrt,stddevMode

    if(mmpi_myid == 0) write(*,*) 'bhi_setup: starting'
    if(mmpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if ( present(mode_opt) ) then
       if ( trim(mode_opt) == 'Analysis' .or. trim(mode_opt) == 'BackgroundCheck') then
         bhi_mode = trim(mode_opt)
         if(mmpi_myid == 0) write(*,*)
         if(mmpi_myid == 0) write(*,*) 'bmatrixHI: Mode activated = ', trim(bhi_mode)
       else
          write(*,*)
          write(*,*) 'mode = ', trim(mode_opt)
          call utl_abort('bmatrixHI: unknown mode')
       end if
    else
       bhi_mode = 'Analysis'
       if(mmpi_myid == 0) write(*,*)
       if(mmpi_myid == 0) write(*,*) 'bmatrixHI: Analysis mode activated (by default)'
    end if

    ! default values for namelist variables
    ntrunc = 108
    scaleFactor(:) = 0.0d0
    scaleFactorLQ(:) = 1.0d0
    scaleFactorCC(:) = 1.0d0
    scaleTG = .true.
    numModeZero = 0
    squareSqrt = .false.
    TweakTG = .false.
    ReadWrite_sqrt = .false.
    stddevMode = 'SP2D'

    nulnam = 0
    ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
    read(nulnam,nml=nambhi,iostat=ierr)
    if ( ierr /= 0 ) call utl_abort( 'bhi_setup: Error reading namelist' )
    if ( mmpi_myid == 0 ) write( *, nml = nambhi )
    ierr = fclos( nulnam )

    do jlev = 1, maxNumLevels
      if( scaleFactor( jlev ) > 0.0d0 ) then 
        scaleFactor( jlev ) = sqrt( scaleFactor( jlev ))
      else
        scaleFactor( jlev ) = 0.0d0
      endif
    enddo

    if ( sum( scaleFactor( 1 : maxNumLevels ) ) == 0.0d0 ) then
      if ( mmpi_myid == 0 ) write(*,*) 'bmatrixHI: scaleFactor=0, skipping rest of setup'
      cvdim_out = 0
      return
    end if

    vco_anl => vco_in
    nLev_M = vco_anl%nlev_M
    nLev_T = vco_anl%nlev_T
    ! need an even number of levels for spectral transform (gstID2)
    if(mod(nLev_T,2).ne.0) then
      nLev_T_even = nLev_T+1
    else
      nLev_T_even = nLev_T
    endif
    if(mmpi_myid == 0) write(*,*) 'BHI_setup: nLev_M, nLev_T, nLev_T_even=',nLev_M, nLev_T, nLev_T_even

    ! check if analysisgrid and covariance file have the same vertical levels
    call vco_SetupFromFile( vco_file,  & ! OUT
                            bFileName )  ! IN
    if (.not. vco_equal(vco_anl,vco_file)) then
      call utl_abort('bmatrixHI: vco from analysisgrid and cov file do not match')
    end if

    status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)
    if(Vcode_anl .ne. 5002 .and. Vcode_anl .ne. 5005) then
      write(*,*) 'Vcode_anl = ',Vcode_anl
      call utl_abort('bmatrixHI: unknown vertical coordinate type!')
    endif

    if (.not. (gsv_varExist(varName='TT') .and.  &
               gsv_varExist(varName='UU') .and.  &
               gsv_varExist(varName='VV') .and.  &
               (gsv_varExist(varName='HU').or.gsv_varExist(varName='LQ')) .and.  &
               gsv_varExist(varName='P0')) ) then
      call utl_abort('bmatrixHI: Some or all weather fields are missing. If it is desired to deactivate the weather assimilation, then all entries of the array SCALEFACTOR in the namelist NAMBHI should be set to zero.')
    end if
    if (.not. gsv_varExist(varName='TG')) then
      write(*,*) 'bmatrixHI: WARNING: The TG field is missing. This must be present when assimilating'
      write(*,*) '                    radiance observations'
    end if

    do jlev = 1, max(nLev_M,nLev_T)
      if(scaleFactorLQ(jlev).gt.0.0d0) then 
        scaleFactorLQ(jlev) = sqrt(scaleFactorLQ(jlev))
      else
        scaleFactorLQ(jlev) = 0.0d0
      endif
    enddo

    do jlev = 1, max(nLev_M,nLev_T)
      if(scaleFactorCC(jlev).gt.0.0d0) then 
        scaleFactorCC(jlev) = sqrt(scaleFactorCC(jlev))
      else
        scaleFactorCC(jlev) = 0.0d0
      endif
    enddo

    if ( trim(bhi_mode) == 'BackgroundCheck' ) then
       cvDim_out = 9999 ! Dummy value > 0 to indicate to the background check (s/r ose_compute_HBHT_ensemble) 
                        ! that Bhi is used
       return
    end if

    numvar3d = 4
    numvar2d = 2

    nLevPtot = nLev_M-1 ! ignore streamfunction at hyb=1, since highly correlated with next level
    nspositVO = 1
    nspositDI = 1*nLev_M+1
    nspositTT = 2*nLev_M+1
    nspositQ  = 2*nLev_M+1*nLev_T+1
    nspositPS = 2*nLev_M+2*nLev_T+1
    nspositTG = 2*nLev_M+2*nLev_T+2
    nkgdim = nLev_M*2 + nLev_T*2 + numvar2d
    nkgdim2 = nkgdim + nLev_T
    if(squareSqrt) then
      nkgdimSqrt = nkgdim2
    else
      nkgdimSqrt = nkgdim
    endif
    nla_mpiglobal = (ntrunc+1)*(ntrunc+2)/2
    
    ni_l = hco_in%ni
    nj_l = hco_in%nj
    AnalGridID = hco_in%EZscintID

    gstID  = gst_setup(ni_l,nj_l,ntrunc,nkgdim)
    gstID2 = gst_setup(ni_l,nj_l,ntrunc,nlev_T_even)
    if(mmpi_myid == 0) write(*,*) 'BHI:returned value of gstID =',gstID
    if(mmpi_myid == 0) write(*,*) 'BHI:returned value of gstID2=',gstID2

    call mmpi_setup_latbands(nj_l, latPerPE, latPerPEmax, myLatBeg, myLatEnd)
    call mmpi_setup_lonbands(ni_l, lonPerPE, lonPerPEmax, myLonBeg, myLonEnd)

    call mmpi_setup_m(ntrunc,mymBeg,mymEnd,mymSkip,mymCount)
    call mmpi_setup_n(ntrunc,mynBeg,mynEnd,mynSkip,mynCount)

    call gst_ilaList_mpiglobal(ilaList_mpiglobal,nla_mpilocal,maxMyNla,gstID,mymBeg,mymEnd,mymSkip,mynBeg,mynEnd,mynSkip)
    call gst_ilaList_mpilocal(ilaList_mpilocal,gstID,mymBeg,mymEnd,mymSkip,mynBeg,mynEnd,mynSkip)

    ! compute mpilocal control vector size
    cvDim_mpilocal = 0
    do jm = mymBeg, mymEnd, mymSkip
      do jn = mynBeg, mynEnd, mynSkip
        if(jm.le.jn) then
          if(jm == 0) then
            ! only real component for jm=0
            cvDim_mpilocal = cvDim_mpilocal + 1*nkgdimSqrt
          else
            ! both real and imaginary components for jm>0
            cvDim_mpilocal = cvDim_mpilocal + 2*nkgdimSqrt
          endif
        endif
      enddo
    enddo
    cvDim_out = cvDim_mpilocal

    ! also compute mpiglobal control vector dimension
    call rpn_comm_allreduce(cvDim_mpilocal,cvDim_mpiglobal,1,"mpi_integer","mpi_sum","GRID",ierr)

    allocate(PtoT(nlev_T+1,nlev_M,nj_l))
    allocate(tantheta(nlev_M,nj_l))
    allocate(rgsig(nj_l,nkgdim))
    allocate(tgstdbg(ni_l,nj_l))
    rgsiguu => rgsig(1:nj_l,nspositVO:nspositVO+nlev_M-1)
    rgsigvv => rgsig(1:nj_l,nspositDI:nspositDI+nlev_M-1)
    rgsigtt => rgsig(1:nj_l,nspositTT:nspositTT+nlev_T-1)
    rgsigq  => rgsig(1:nj_l,nspositQ :nspositQ +nlev_T-1)
    rgsigps => rgsig(1:nj_l,nspositPS)
    allocate(rgsigtb(nj_l,nlev_T))
    allocate(rgsigpsb(nj_l))
    allocate(corns(nkgdim2,nkgdim2,0:ntrunc))
    allocate(rstddev(nkgdim2,0:ntrunc))

    if(mmpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    zps = 101000.D0
    status = vgd_levels( vco_anl%vgrid, ip1_list=vco_anl%ip1_M, levels=pressureProfile_M, &
                         sfc_field=zps, in_log=.false.)
    status = vgd_levels( vco_anl%vgrid, ip1_list=vco_anl%ip1_T, levels=pressureProfile_T, &
                         sfc_field=zps, in_log=.false.)

    llfound = .false.
    nlev_bdl = 0
    do jlev = 1, nlev_M
      if(.not.llfound .and. (pressureProfile_M(jlev) .ge. rlimlv_bdl  )) then
        nlev_bdl = jlev
        llfound = .true.
      endif
    enddo

    inquire(file=bFileName,exist=lExists)
    IF ( lexists )then
      ierr = fnom(nulbgst,bFileName,'RND+OLD+R/O',0)
      if ( ierr == 0 ) then
        ierr =  fstouv(nulbgst,'RND+OLD')
      else
        call utl_abort('BHI_setup:NO BACKGROUND STAT FILE!!')
      endif
    endif

    call BHI_rdspPtoT

    call BHI_readcorns2
    if(mmpi_myid == 0) write(*,*) 'Memory Used (after readcorns2): ',get_max_rss()/1024,'Mb'

    call BHI_sutg

    if(stddevmode == 'GD2D') then
      call BHI_rdstd
    elseif(stddevMode == 'SP2D') then
      call BHI_rdspstd_newfmt
    else
      call utl_abort('BHI_setup: unknown stddevMode')
    endif

    call BHI_scalestd

    call BHI_sucorns2
    if(mmpi_myid == 0) write(*,*) 'Memory Used (after sucorns2): ',get_max_rss()/1024,'Mb'

    ierr = fstfrm(nulbgst)
    ierr = fclos(nulbgst)

    if(mmpi_myid == 0) write(*,*) 'END OF BHI_SETUP'

    initialized = .true.

  END SUBROUTINE BHI_setup


  subroutine bhi_getScaleFactor(scaleFactor_out)
    implicit none
    real(8) :: scaleFactor_out(:)
    integer :: jlev

    do jlev = 1, max(nLev_M,nLev_T)
      scaleFactor_out(jlev) = scaleFactor(jlev)
    enddo

  end subroutine bhi_getScaleFactor


  SUBROUTINE BHI_scalestd
    implicit none

    integer :: jlev, jlon, jlat, shift_level, Vcode_anl, status

    status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)
    if(Vcode_anl == 5002) then
      shift_level = 1
    else
      shift_level = 0
    endif

    do jlev = 1, nlev_M
      do jlat = 1, nj_l
        rgsiguu(jlat,jlev) =                                 scaleFactor(jlev+shift_level)*rgsiguu(jlat,jlev)
        rgsigvv(jlat,jlev) = scaleFactorCC(jlev+shift_level)*scaleFactor(jlev+shift_level)*rgsigvv(jlat,jlev)
      enddo
    enddo
    do jlev = 1, nlev_T
      do jlat = 1, nj_l
        rgsigtt(jlat,jlev) =                     scaleFactor(jlev)*rgsigtt(jlat,jlev)
        rgsigq(jlat,jlev)  = scaleFactorLQ(jlev)*scaleFactor(jlev)*rgsigq(jlat,jlev)
        rgsigtb(jlat,jlev) =                     scaleFactor(jlev)*rgsigtb(jlat,jlev)
      enddo
    enddo
    do jlat = 1, nj_l
      rgsigpsb(jlat) = scaleFactor(max(nLev_M,nLev_T))*rgsigpsb(jlat)
      rgsigps(jlat)  = scaleFactor(max(nLev_M,nLev_T))*rgsigps(jlat)
    enddo
    ! User has the option to not scale down the STDDEV of TG (because underestimated in Benkf)
    if(scaleTG) then
    do jlat = 1, nj_l
      do jlon = 1, ni_l
        tgstdbg(jlon,jlat) = scaleFactor(max(nLev_M,nLev_T))*tgstdbg(jlon,jlat)
      enddo
    enddo
    endif

  END SUBROUTINE BHI_scalestd


  SUBROUTINE BHI_SUCORNS2
    implicit none

    real(8) :: eigenval(nkgdim2), eigenvec(nkgdim2,nkgdim2), result(nkgdim2,nkgdim2)
    real(8) :: eigenvalsqrt(nkgdim2), eigenvec2(nkgdim2,nkgdim2), eigenvalsqrt2(nkgdim2)

    integer :: jlat,jn,jk1,jk2,jk3
    integer :: ilwork,info,klatPtoT
    integer :: iulcorvert, ikey, nsize

    real(8) :: zwork(2*4*nkgdim2)
    real(8) :: ztt(nlev_T,nlev_T,(ntrunc+1)),ztpsi(nlev_T,nlev_M,(ntrunc+1))
    real(8) :: ztlen,zcorr,zr,zpres1,zpres2
    real(8) :: zfact,zfact2,zcoriolis,zpsips(nLevPtoT)
    real(8) :: zpsi(nlev_M,nlev_M),zfacttb(nj_l,nlev_T),zfactpsb(nj_l)
    real(8) :: corvert(nkgdim2,nkgdim2)
    real(8),allocatable :: corns_temp(:,:,:)
    logical :: lldebug, lfound_sqrt

    ! standard file variables
    integer :: ini,inj,ink, inpas, inbits, idatyp, ideet 
    integer :: ip1,ip2,ip3,ig1,ig2,ig3,ig4,iswa,ilength,idltf
    integer :: iubc,iextr1,iextr2,iextr3,ierr
    integer :: idateo
    character(len=2)  :: cltypvar
    character(len=1)  :: clgrtyp
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fstprm,fstinf
    integer :: fnom,fstouv,fstfrm,fclos

    lldebug = .false.

    iulcorvert = 0
    if(mmpi_myid==0) then
      ierr = fnom(iulcorvert,'corvert_modular.fst','RND',0)
      ierr = fstouv(iulcorvert,'RND')
    endif

    klatPtoT = 1
    zfactpsb(:) = 0.0d0
    zfacttb(:,:) = 0.0d0

    if(lldebug) then
      do jk1 = 1, nlev_T
        do jk2 = 1, nlevPtoT
          write(622,*) jk1,jk2,klatPtoT,PtoT(jk1,jk2,klatPtoT)
        enddo
      enddo
    endif

    ! explicitly compute the balanced temperature and temperature-psi correlations

    do jn = 0, ntrunc

      ztpsi(:,:,jn+1) = 0.0d0
      ztt(:,:,jn+1) = 0.0d0
      do jk1 = 1, nlevPtoT
        do jk2 = 1, nlev_T
          do jk3 = 1, nlevPtoT
            ztpsi(jk2,jk1,jn+1) = ztpsi(jk2,jk1,jn+1)+PtoT(jk2,jk3,klatPtoT)*corns(jk3,jk1,jn)
          enddo
        enddo
      enddo
      if(nlevPtoT.lt.nlev_M) then
        do jk1 = (nlevPtoT+1), nlev_M
          do jk2 = 1, nlev_T
            ztpsi(jk2,jk1,jn+1) = ztpsi(jk2,nlevPtoT,jn+1)
          enddo
        enddo
      endif
      do jk1 = 1, nlev_T
        do jk2 = 1, nlev_T
          do jk3 = 1, nlevPtoT
            ztt(jk2,jk1,jn+1) = ztt(jk2,jk1,jn+1)+ztpsi(jk2,jk3,jn+1)*PtoT(jk1,jk3,klatPtoT)
          enddo
        enddo
      enddo
    enddo

    if(lldebug) then
      write(620,*) ztt
      write(621,*) ztpsi
    endif

    ! fill in blocks for balance temperature

    do jn = 0, ntrunc
      do jk1 = 1, nlev_T
        do jk2 = 1, nlev_T
          corns(nkgdim+jk2,nkgdim+jk1,jn) = ztt(jk2,jk1,jn+1)
        enddo
      enddo
      do jk1 = 1, nlev_M
        do jk2 = 1, nlev_T
          corns(       jk1,nkgdim+jk2,jn) = ztpsi(jk2,jk1,jn+1)
          corns(nkgdim+jk2,       jk1,jn) = ztpsi(jk2,jk1,jn+1)
        enddo
      enddo
    enddo

    ! Save un-localized PSI correlations
    do jk2 = 1, nlev_M
      do jk1 = 1, nlev_M
        zpsi(jk1,jk2) = 0.0d0
        do jn = 0, ntrunc
          zpsi(jk1,jk2) = zpsi(jk1,jk2)+((2*jn+1)*corns(jk1,jk2,jn))
        enddo
      enddo
    enddo

    ! Apply vertical localization to corrns

    ! unbalanced temperature
    ztlen = rvlocunbalt
    if(ztlen.gt.0.0d0) then
      ! calculate 5'th order function (from Gaspari and Cohn)
      do jk1 = 1, nlev_T
        zpres1 = log(pressureProfile_T(jk1))
        do jk2 = 1, nlev_T
          zpres2 = log(pressureProfile_T(jk2))
          zr = abs(zpres2 - zpres1)
          zcorr = gasparicohn(ztlen,zr)
          do jn = 0, ntrunc
            corns(jk1+2*nlev_M,jk2+2*nlev_M,jn)  =   &
                 corns(jk1+2*nlev_M,jk2+2*nlev_M,jn)*zcorr
          enddo
        enddo
      enddo
    endif

    ! balanced temperature
    ztlen = rvlocbalt
    if(ztlen.gt.0.0d0) then
      ! calculate 5'th order function (from Gaspari and Cohn)
      do jk1 = 1, nlev_T
        zpres1 = log(pressureProfile_T(jk1))
        do jk2 = 1, nlev_T
          zpres2 = log(pressureProfile_T(jk2))
          zr = abs(zpres2 - zpres1)
          zcorr = gasparicohn(ztlen,zr)
          do jn = 0, ntrunc
            corns(jk1+nkgdim,jk2+nkgdim,jn)  =        &
                 corns(jk1+nkgdim,jk2+nkgdim,jn)*zcorr
          enddo
        enddo
      enddo
    endif

    ! streamfunction 
    ztlen = rvlocpsi    ! specify length scale (in units of ln(Pressure))
    if(ztlen.gt.0.0d0) then
      ! calculate 5'th order function (from Gaspari and Cohn)
      do jk1 = 1, nlev_M
        zpres1 = log(pressureProfile_M(jk1))
        do jk2 = 1, nlev_M
          zpres2 = log(pressureProfile_M(jk2))
          zr = abs(zpres2 - zpres1)
          zcorr = gasparicohn(ztlen,zr)
          do jn = 0, ntrunc
            corns(jk1,jk2,jn) = corns(jk1,jk2,jn)*zcorr
          enddo
        enddo
      enddo
    endif

    ! temp-psi cross-correlations
    ztlen = rvlocpsitt    ! specify length scale (in units of ln(Pressure))
    if(ztlen.gt.0.0d0) then
      ! calculate 5'th order function (from Gaspari and Cohn)
      do jk1 = 1, nlev_M
        zpres1 = log(pressureProfile_M(jk1))
        do jk2 = 1, nlev_T
          zpres2 = log(pressureProfile_T(jk2))
          zr = abs(zpres2 - zpres1)
          zcorr = gasparicohn(ztlen,zr)
          do jn = 0, ntrunc
            corns(jk1,jk2+nkgdim,jn) = corns(jk1,jk2+nkgdim,jn)*zcorr
            corns(jk2+nkgdim,jk1,jn) = corns(jk2+nkgdim,jk1,jn)*zcorr
          enddo
        enddo
      enddo
    endif

    ! velocity potential (unbalanced)
    ztlen = rvlocchi    ! specify length scale (in units of ln(Pressure))
    if(ztlen.gt.0.0d0) then
      ! calculate 5'th order function (from Gaspari and Cohn)
      do jk1 = 1, nlev_M
        zpres1 = log(pressureProfile_M(jk1))
        do jk2 = 1, nlev_M
          zpres2 = log(pressureProfile_M(jk2))
          zr = abs(zpres2 - zpres1)
          zcorr = gasparicohn(ztlen,zr)
          do jn = 0, ntrunc
            corns(jk1+nlev_M,jk2+nlev_M,jn) = corns(jk1+nlev_M,jk2+nlev_M,jn)*zcorr
          enddo
        enddo
      enddo
    endif

    ! cross-correlation t'-ps'
    if(.true.) then
    ztlen = rvlocunbalt    ! specify length scale (in units of ln(Pressure))
    if(ztlen.gt.0.0d0) then
      ! calculate 5'th order function (from Gaspari and Cohn)
      zpres1 = log(pressureProfile_T(nlev_T))
      do jk2 = 1, nlev_T
        zpres2 = log(pressureProfile_T(jk2))
        zr = abs(zpres2 - zpres1)
        zcorr = gasparicohn(ztlen,zr)
        do jn = 0, ntrunc
          corns(1+2*nlev_M+2*nlev_T,jk2+2*nlev_M,jn)  =       &
               corns(1+2*nlev_M+2*nlev_T,jk2+2*nlev_M,jn)*zcorr
          corns(jk2+2*nlev_M,1+2*nlev_M+2*nlev_T,jn)  =       &
               corns(jk2+2*nlev_M,1+2*nlev_M+2*nlev_T,jn)*zcorr
        enddo
      enddo
    endif
    endif

    ! humidity
    ztlen = rvloclq    ! specify length scale (in units of ln(Pressure))
    if(ztlen.gt.0.0d0) then
      ! calculate 5'th order function (from Gaspari and Cohn)
      do jk1 = 1, nlev_T
        zpres1 = log(pressureProfile_T(jk1))
        do jk2 = 1, nlev_T
          zpres2 = log(pressureProfile_T(jk2))
          zr = abs(zpres2 - zpres1)
          zcorr = gasparicohn(ztlen,zr)
          do jn = 0, ntrunc
            corns(jk1+2*nlev_M+nlev_T,jk2+2*nlev_M+nlev_T,jn)  =       &
                 corns(jk1+2*nlev_M+nlev_T,jk2+2*nlev_M+nlev_T,jn)*zcorr
          enddo
        enddo
      enddo
    endif

    ! compute total vertical correlations (including for balanced temperature)
    if(.true.) then
      do jk2 = 1, nkgdim2
        do jk1 = 1, nkgdim2
          corvert(jk1,jk2) = 0.0d0
          do jn = 0, ntrunc
            corvert(jk1,jk2) = corvert(jk1,jk2)+((2*jn+1)*corns(jk1,jk2,jn))
          enddo
        enddo
      enddo

      if(lldebug) then
        write(701,*) corvert
        write(702,*) zpsi
      endif

      if(mmpi_myid == 0) then
        ikey = fstinf(NULBGST,ini,inj,ink,-1,'CORRNS',-1,0,-1,' ','ZZ')
        ierr = fstprm(ikey,idateo,ideet,inpas,ini,inj,ink, inbits        &
             ,idatyp,ip1,ip2,ip3,cltypvar,clnomvar,cletiket,clgrtyp      &
             ,ig1,ig2,ig3,ig4,iswa,ilength,idltf,iubc,iextr1,iextr2      &
             ,iextr3)

        ini = nkgdim2
        inj = nkgdim2
        ink = 1
        ip1 = 0
        ip2 = ntrunc
        ip3 = 0
        clnomvar = 'ZV'
        cletiket = 'CORVERT'
        idatyp = 5

        ierr = utl_fstecr(corvert, -inbits, iulcorvert, idateo    &
             , ideet,inpas, ini, inj, ink, ip1, ip2, ip3, cltypvar,     &
             clnomvar,cletiket,clgrtyp,ig1, ig2, ig3, ig4, idatyp,      &
             .true.)

      endif

      ! Modify RGSIGTB to obtain correct sigma_Tb
      do jk1 = 1, nlev_T
        zfact = corvert(jk1+nkgdim,jk1+nkgdim)
        do jlat = 1, nj_l
          zcoriolis = abs(2.d0*ec_Omega*gst_getrmu(jlat,gstID))
          if(zfact.gt.0.0d0.and.zcoriolis.ne.0.0d0) then
            zfact2 = 1.0d0/(zfact*zcoriolis*zcoriolis)
          else 
            zfact2 = 0.0d0
          endif
          zfacttb(jlat,jk1) = zfacttb(jlat,jk1)+zfact2
        enddo
      enddo

      ! Modify RGSIGPSB to obtain correct sigma_PSb
      do jlat = 1, nj_l
        do jk2 = 1, nlevPtoT
          zpsips(jk2) = 0.0d0
          do jk1 = 1, nlevPtoT
            zpsips(jk2) = zpsips(jk2)+PtoT(nlev_T+1,jk1,klatPtoT)*zpsi(jk1,jk2)
          enddo
        enddo
        zfact = 0.0d0
        do jk1 = 1, nlevPtoT
          zfact = zfact+PtoT(nlev_T+1,jk1,klatPtoT)*zpsips(jk1)
        enddo
        zcoriolis = abs(2.d0*ec_Omega*gst_getrmu(jlat,gstID))
        if(zfact.gt.0.0d0.and.zcoriolis.ne.0.0d0) then
          zfact2 = 1.0d0/(zfact*zcoriolis*zcoriolis)
        else 
          zfact2 = 0.0d0
        endif
        zfactpsb(jlat) = zfactpsb(jlat)+zfact2
      enddo
    endif

    ! Modify RGSIGTB and RGSIGPSB to obtain correct sigma_Tb and sigma_Psb
    do jlat = 1, nj_l
      if(zfactpsb(jlat).gt.0.0d0) then
        rgsigpsb(jlat) = rgsigpsb(jlat)*sqrt(zfactpsb(jlat))
      else
        rgsigpsb(jlat) = 0.0d0
      endif          
      do jk1 = 1, nlev_T
        if(zfacttb(jlat,jk1).gt.0.0d0) then
          rgsigtb(jlat,jk1) = rgsigtb(jlat,jk1)*sqrt(zfacttb(jlat,jk1))
        else
          rgsigtb(jlat,jk1) = 0.0d0
        endif
      enddo
    enddo

    ! compute square-root of corns for each total wavenumber
    allocate(corns_temp(nkgdim2,nkgdim2,0:ntrunc))
    corns_temp(:,:,:)=0.0d0
    do jn = mmpi_myid, ntrunc, mmpi_nprocs

      do jk1 = 1, nkgdim2
         do jk2 = 1, nkgdim2
            eigenvec(jk2,jk1) = corns(jk2,jk1,jn)
         enddo
      enddo

      ! CALCULATE EIGENVALUES AND EIGENVECTORS.
      ilwork = 4*nkgdim2*2
      call dsyev('V','U',nkgdim2,eigenvec,nkgdim2,eigenval,zwork,ilwork,info)
      if(info.ne.0) then
        write(*,*) 'bhi_sucorns2: non-zero value of info =',info,' returned by dsyev for wavenumber ',jn
        call utl_abort('BHI_SUCORNS')
      endif

      ! set selected number of eigenmodes to zero
      if(numModeZero.gt.0) then
        write(*,*) 'bhi_sucorns2: setting ',numModeZero,' eigenvalues to zero for wavenumber n=',jn
        write(*,*) 'bhi_sucorns2: original eigenvalues=',eigenval(:)
        do jk1 = 1, numModeZero
          eigenval(jk1) = 0.0d0
        enddo
        write(*,*) 'bhi_sucorns2: modified eigenvalues=',eigenval(:)
      endif

      do jk1 = 1, nkgdim2
        if(eigenval(jk1).lt.1.0d-15) then
          eigenvalsqrt(jk1) = 0.0d0
        else
          eigenvalsqrt(jk1) = sqrt(eigenval(jk1))
        endif
      enddo
 
      ! Reverse the order of E-Values if old formulation (for compatibility)
      if(.not.squareSqrt) then
        do jk1 = 1, nkgdim2
          eigenvalsqrt2(jk1) = eigenvalsqrt(nkgdim2-jk1+1)
          do jk2 = 1, nkgdim2
            eigenvec2(jk2,jk1) = eigenvec(jk2,nkgdim2-jk1+1)
          enddo
        enddo
        eigenvalsqrt(:) = eigenvalsqrt2(:)
        eigenvec(:,:) = eigenvec2(:,:)
      endif

      ! compute E * lambda^1/2
      result(:,:) = 0.0d0
      do jk1 = 1, nkgdimSqrt
         do jk2 = 1, nkgdim2
            result(jk2,jk1) = eigenvec(jk2,jk1)*eigenvalsqrt(jk1)
         enddo
      enddo

      ! compute (E * lambda^1/2) * E^T if new formulation
      if(squareSqrt) then
        do jk1 = 1, nkgdim2
          do jk2 = 1, nkgdim2
            do jk3 = 1, nkgdim2
              corns_temp(jk2,jk1,jn) = corns_temp(jk2,jk1,jn) + result(jk2,jk3)*eigenvec(jk1,jk3)
            enddo
          enddo
        enddo
      else
        corns_temp(:,:,jn) = result(:,:)
      endif

      !if(jn == 30) then
      !  write(200,*) corns(:,:,jn)
      !  write(201,*) corns_temp(:,:,jn)
      !  write(202,*) eigenval(:)
      !  write(203,*) eigenvec(:,:)
      !  call flush(200)
      !  call flush(201)
      !  call flush(202)
      !  call flush(203)
      !endif

    enddo ! jn

    nsize = nkgdim2*nkgdim2*(ntrunc+1)
    call rpn_comm_allreduce(corns_temp,corns,nsize,"mpi_double_precision","mpi_sum","GRID",ierr)
    deallocate(corns_temp)

    if(mmpi_myid==0) then
      ierr = fstfrm(iulcorvert)
      ierr = fclos(iulcorvert)
    endif

    if(ReadWrite_sqrt) then
      ! if desired, read precomputed sqrt of corns which overwrites computed matrix
      call readcorns_sqrt(lfound_sqrt)
      if(.not.lfound_sqrt) then
        ! if precomputed sqrt does not exist in stats, then write it out to a separate file
        call writecorns_sqrt
      endif
    endif

  END SUBROUTINE BHI_SUCORNS2


  SUBROUTINE WRITECORNS_SQRT
    implicit none

    integer :: jn, nulcorns_sqrt, ierr

    ! standard file variables
    integer :: ip1,ip2,ip3
    integer :: idateo, ipak, idatyp
    integer :: fnom, fstouv,  fstfrm, fclos

    write(*,*) 'WRITECORNS_SQRT: CORNS_SQRT will be written to file corns_sqrt.fst for NTRUNC =', ntrunc

    if(mmpi_myid==0) then
      nulcorns_sqrt = 0
      ierr = fnom(nulcorns_sqrt,'corns_sqrt.fst','RND',0)
      ierr = fstouv(nulcorns_sqrt,'RND')

      ipak = -32
      idatyp = 5
      ip1 = 0
      ip3 = ntrunc
      idateo = 0

      do jn = 0, ntrunc
        ip2 = jn
        ierr = utl_fstecr(corns(1,1,jn),ipak,nulcorns_sqrt,idateo,0,0,nkgdim2,nkgdim2,1,  &
                       ip1,ip2,ip3,'X','ZZ','CORNS_SQRT','X',0,0,0,0,idatyp,.true.)
      enddo

      ierr = fstfrm(nulcorns_sqrt)  
      ierr = fclos(nulcorns_sqrt)
    endif

  END SUBROUTINE WRITECORNS_SQRT


  SUBROUTINE READCORNS_SQRT(lfound_sqrt)
    implicit none
    logical :: lfound_sqrt

    integer :: jn, icornskey
    integer :: jcol,jrow
    real(8), allocatable :: zcornssrc(:,:)

    ! standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3
    integer :: idateo
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket

    allocate(zcornssrc(nkgdim2,nkgdim2))

    write(*,*) 'READCORNS_SQRT: looking for CORNS_SQRT with NTRUNC =', ntrunc
    do jn = 0, ntrunc

      ! Looking for FST record parameters..
      idateo = -1
      cletiket = 'CORNS_SQRT'
      ip1 = -1
      ip2 = jn
      ip3 = ntrunc
      cltypvar = 'X'
      clnomvar = 'ZZ'
      icornskey = utl_fstlir(ZCORNSSRC,nulbgst,INI,INJ,INK,idateo,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)

      if( jn == 0 ) then
        if(icornskey .lt.0 ) then
          write(*,*) 'READCORNS_SQRT: CORNS_SQRT not found in stats file, use computed sqrt'
          lfound_sqrt = .false.
          return
        else
          write(*,*) 'READCORNS_SQRT: CORNS_SQRT found in stats file, will use it instead of computed sqrt'
          lfound_sqrt = .true.
        endif
      endif

      if (ini .ne. nkgdim2 .or. inj .ne. nkgdim2) then
        call utl_abort('READCORNS2: BG stat levels inconsitencies')
      endif

      do jcol = 1, nkgdim2
        do jrow = 1, nkgdim2
          corns(jrow,jcol,jn) = zcornssrc(jrow,jcol)
        enddo
      enddo

    enddo

    deallocate(zcornssrc)

  END SUBROUTINE READCORNS_SQRT


  FUNCTION GASPARICOHN(ztlen,zr)

    real(8)  :: gasparicohn
    real(8)  :: ztlen,zr,zlc

    zlc = ztlen/2.0d0
    if(zr.le.zlc) then
      gasparicohn = -0.250d0*(zr/zlc)**5+0.5d0*(zr/zlc)**4             &
                  +0.625d0*(zr/zlc)**3-(5.0d0/3.0d0)*(zr/zlc)**2+1.0d0
    elseif(zr.le.(2.0d0*zlc)) then
      gasparicohn = (1.0d0/12.0d0)*(zr/zlc)**5-0.5d0*(zr/zlc)**4         &
                  +0.625d0*(zr/zlc)**3+(5.0d0/3.0d0)*(zr/zlc)**2       &
                  -5.0d0*(zr/zlc)+4.0d0-(2.0d0/3.0d0)*(zlc/zr)
    else
      gasparicohn = 0.0d0
    endif
    if(gasparicohn.lt.0.0d0) gasparicohn = 0.0d0

  END FUNCTION GASPARICOHN


  SUBROUTINE BHI_CALCCORR(zgd,pcscl,klev)
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

    if (kcorrtyp == 1) then
      ! Gaussian correlation
      do  jlev = 1, klev
        dlc = 1.d0/dble(pcscl(jlev))
        dlc = 0.5d0*dlc*dlc
        do  jlat = myLatBeg, myLatEnd
          zr = ec_ra * acos(gst_getRmu(jlat,gstID))
          dlcorr = dexp(-(zr**2)*dlc)
          do  jlon = myLonBeg, myLonEnd
            zgd(jlon,jlat,jlev) = dlcorr
          enddo
        enddo
      enddo
    elseif (kcorrtyp == 2) then
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
          enddo
        enddo
      enddo
    else
      call utl_abort('CALCCORR- Undefined correlation type')
    endif

  END SUBROUTINE BHI_calcCorr


  SUBROUTINE BHI_SUTG
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
    real(8) :: rcscltg_vec(nlev_T_even)
    real(8) :: zabs, zpole, dlfac
    real(8) :: zsp_mpilocal(nla_mpilocal,2,nlev_T_even)
    real(8) :: zgd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nlev_T_even)
    real(8) :: zsp_mpiglobal(nla_mpiglobal,2,1)

    real(8),allocatable :: my_zsp_mpiglobal(:,:,:)

    real(4), allocatable :: TrialLandSeaMask(:,:), TrialSeaIceMask(:,:)
    real(4), allocatable :: AnalLandSeaMask(:,:), AnalSeaIceMask(:,:)

    ! standard file variables
    integer :: ini,inj,ink, idatyp
    integer :: ip1,ip2,ip3
    integer :: ierr,ntrials
    integer :: idateo
    integer :: fstprm,fstinf,iultg,fnom,fclos,fstouv,fstfrm
    integer :: ip1s(1), nulbgsts(1)

    integer :: TrlmNumberWanted
    integer :: fstlir, key, nultrl, ni_trial, nj_trial
    integer :: deet, npas, nbits, datyp
    integer :: ig1,ig2,ig3,ig4,swa, lng, dltf, ubc
    integer :: extra1, extra2, extra3
    
    integer :: TrialGridID

    character(len=2)  :: cltypvar
    character(len=1)  :: clgrtyp
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket

    character(len=2) :: flnum
    character(len=128) :: trialfile

    logical :: trialExists

    !
    !- 1.  Reading and processing TG standard deviations
    !

    !- 1.1 Read the Std. Dev. field from ./bgcov file
    clnomvar = 'TG'
    idateo = -1
    inmxlev = 1
    ntrials = 1
    nulbgsts(1) = nulbgst

    call utl_getfldprm(IP1S, IP2, IP3, INLEV, CLETIKET, CLTYPVAR, ITGGID,  &
                       clnomvar, idateo, inmxlev, nulbgsts, ip1style, ip1kind,  &
                       ntrials, koutmpg)
    ip1 = ip1s(1)

    ierr = ezgprm(itggid,CLGRTYP,INI,INJ,IG1,IG2,IG3,IG4)
    allocate(dltg(ini,inj))

    ikey = utl_fstlir(dltg,koutmpg,ini,inj,ink,idateo,cletiket,ip1,   &
           ip2, ip3, cltypvar, clnomvar)

    !- 1.2 Rearrange the data according to the analysis grid (if necessary)
    if(clgrtyp == 'G' .and. ni_l == ini .and. nj_l == inj .and. ig1 == 0  &
          .and. ig2 ==0 .and. ig3 == 0 .and.ig4 == 0) then
      !- 1.2.1 The std. dev. on the analysis grid 
      do jlat = 1, nj_l
        do jlon = 1,ni_l
          tgstdbg(jlon,jlat) = dltg(jlon,jlat)
        enddo
      enddo

    elseif(clgrtyp == 'G' .and. ni_l == ini .and. nj_l == inj .and. ig1 ==   &
            0 .and. ig2 ==1 .and. ig3 == 0 .and.ig4 == 0) then
       !- 1.2.2 flipped Gaussian grid no longer supported
       call utl_abort('bhi_sutg: The flipped Gaussian grid is no longer supported!')

    else
       !- 1.2.3 The std. dev. are NOT on the analysis grid. Interpolation is needed
       iset = ezdefset(AnalGridID,itggid)
       if ( TweakTG ) then
          ierr = utl_ezsint(tgstdbg,dltg,interpDegree='NEAREST')
       else
          ierr = utl_ezsint(tgstdbg,dltg,interpDegree='CUBIC')
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
           call utl_abort('BMatrixHI : DID NOT FIND A TRIAL FIELD FILE')
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
      clnomvar = 'MG'

      key = fstinf( nultrl,                                             & ! IN
                    ni_trial, nj_trial, ink,                            & ! OUT
                    idateo, cletiket, ip1, ip2, ip3, cltypvar, clnomvar ) ! IN

      if (key < 0) then
         write(*,*)
         write(*,*) 'bMatrixHI: Unable to find trial field = ',clnomvar
         call utl_abort('BMatrixHI')
      end if

      ierr = fstprm( key,                                                 & ! IN
                     idateo, deet, npas, ni_trial, nj_trial, ink, nbits,  & ! OUT
                     datyp, ip1, ip2, ip3, cltypvar, clnomvar, cletiket,  & ! OUT
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
      clnomvar = 'MG'
      ierr = fstlir(TrialLandSeaMask, nultrl, ini, inj, ink,  &
                    idateo ,cletiket, ip1, ip2, ip3, cltypvar, clnomvar)
      if ( ierr < 0 ) then
         write(*,*)
         write(*,*) 'bMatrixHI: Unable to read trial field = ',clnomvar
         call utl_abort('BMatrixHI : fstlir failed')
      end if

      if (ini /= ni_trial .or. inj /= nj_trial) then
          write(*,*)
          write(*,*) 'bMatrixHI: Invalid dimensions for ...'
          write(*,*) 'nomvar      =', trim(clnomvar)
          write(*,*) 'etiket      =', trim(cletiket)
          write(*,*) 'ip1         =', ip1
          write(*,*) 'Found ni,nj =', ini, inj 
          write(*,*) 'Should be   =', ni_trial, nj_trial
          call utl_abort('bMatrixHI')
        end if

      clnomvar = 'GL'
      ierr = fstlir(TrialSeaIceMask, nultrl, ini, inj, ink,  &
                    idateo ,cletiket, ip1, ip2, ip3, cltypvar, clnomvar)
      if ( ierr < 0 ) then
         write(*,*)
         write(*,*) 'bMatrixHI: Unable to read trial field = ',clnomvar
         call utl_abort('BMatrixHI : fstlir failed')
      end if

      if (ini /= ni_trial .or. inj /= nj_trial) then
          write(*,*)
          write(*,*) 'bMatrixHI: Invalid dimensions for ...'
          write(*,*) 'nomvar      =', trim(clnomvar)
          write(*,*) 'etiket      =', trim(cletiket)
          write(*,*) 'ip1         =', ip1
          write(*,*) 'Found ni,nj =', ini, inj 
          write(*,*) 'Should be   =', ni_trial, nj_trial
          call utl_abort('bMatrixHI')
      end if

      TrialGridID  = ezqkdef( ni_trial, nj_trial, clgrtyp, ig1, ig2, ig3, ig4, nultrl )   ! IN
 
      ierr = fstfrm(nultrl)  
      ierr = fclos(nultrl)

      !- Interpolate to the Analysis Grid
      allocate(AnalLandSeaMask(ni_l, nj_l))
      allocate(AnalSeaIceMask(ni_l, nj_l))

      ierr = ezdefset(AnalGridID     , TrialGridID     ) ! IN,  IN

      ! Nearest-neighbor interpolation
      ierr = utl_ezsint(AnalLandSeaMask, TrialLandSeaMask, interpDegree='NEAREST') ! OUT, IN
      ierr = utl_ezsint(AnalSeaIceMask , TrialSeaIceMask, interpDegree='NEAREST' ) ! OUT, IN

      deallocate(TrialLandSeaMask)
      deallocate(TrialSeaIceMask)

      !- Modify the input/regridded TG Std. Dev.
      do jlat = 1, nj_l
         do jlon = 1,ni_l

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
        allocate(tgstdbg_tmp(ni_l,nj_l))
        do jlat = 1, nj_l
          do jlon = 1,ni_l
             tgstdbg_tmp(jlon,jlat) = tgstdbg(jlon,jlat)
          end do
        end do

        iultg = 0
        ierr = fnom(iultg,'tg_stddev_of_the_day.fst','RND',0)
        ierr = fstouv(iultg,'RND')

        ini = ni_l
        inj = nj_l
        ink = 1
        ip1 = 0
        ip2 = 0
        ip3 = 0
        idateo = 0
        cltypvar = 'E'
        clnomvar = 'TG'
        cletiket = 'TWEAK_STDDEV'
        clgrtyp = 'G'
        ig1 = 0
        ig2 = 0
        ig3 = 0
        ig4 = 0
        idatyp = 1

        ierr = utl_fstecr(tgstdbg_tmp, -32, iultg, idateo,         &
                       0, 0, ini, inj, ink, ip1, ip2, ip3, cltypvar,         &
                       clnomvar,cletiket,clgrtyp,ig1, ig2, ig3, ig4, idatyp, &
                       .true.)

        do jlat = 1, nj_l
          do jlon = 1,ni_l
             tgstdbg_tmp(jlon,jlat) = AnalLandSeaMask(jlon,jlat)
          end do
        end do
        cltypvar = 'P'
        clnomvar = 'MG'
        cletiket = 'TRIAL2ANAL'
        ierr = utl_fstecr(tgstdbg_tmp, -32, iultg, idateo,         &
                       0, 0, ini, inj, ink, ip1, ip2, ip3, cltypvar,         &
                       clnomvar,cletiket,clgrtyp,ig1, ig2, ig3, ig4, idatyp, &
                       .true.)

        do jlat = 1, nj_l
          do jlon = 1,ni_l
             tgstdbg_tmp(jlon,jlat) = AnalSeaIceMask(jlon,jlat)
          end do
        end do       
        cltypvar = 'P'
        clnomvar = 'GL'
        cletiket = 'TRIAL2ANAL'
        ierr = utl_fstecr(tgstdbg_tmp, -32, iultg, idateo,         &
                       0, 0, ini, inj, ink, ip1, ip2, ip3, cltypvar,         &
                       clnomvar,cletiket,clgrtyp,ig1, ig2, ig3, ig4, idatyp, &
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
    enddo

    ! 2.4.2  Compute correlations in physical space
    rcscltg_vec(:) = rcscltg(1)
    call BHI_calccorr(zgd,rcscltg_vec,nlev_T_even)

    ! 2.4.3  Bring back the result in spectral space
    call gst_setID(gstID2)
    call gst_reespe(zsp_mpilocal,zgd)

    ! and make the result mpiglobal
    do jm = mymBeg, mymEnd, mymSkip
      do jn = mynBeg, mynEnd, mynSkip
        if(jm.le.jn) then
          ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
          ila_mpilocal  = ilaList_mpilocal(ila_mpiglobal)
          my_zsp_mpiglobal(ila_mpiglobal,:,1) = zsp_mpilocal(ila_mpilocal,:,1)
        endif
      enddo
    enddo
    nsize = 2*nla_mpiglobal
    call rpn_comm_allreduce(my_zsp_mpiglobal(:,:,1),zsp_mpiglobal(:,:,1),nsize,"mpi_double_precision","mpi_sum","GRID",ierr)
    deallocate(my_zsp_mpiglobal) 
    ! 2.4.4  Check positiveness
    llpb = .false.
    do jla = 1, ntrunc+1
      zabs = abs(zsp_mpiglobal(jla,1,1))
      llpb = llpb.or.((zsp_mpiglobal(jla,1,1).lt.0.).and.(zabs.gt.epsilon(zabs)))
    enddo
    if(llpb) then
      call utl_abort(' AUTOCORRELATION  NEGATIVES')
    endif
    do jla = 1, ntrunc+1
      zsp_mpiglobal(jla,1,1) = abs(zsp_mpiglobal(jla,1,1))
    enddo

    zpole = 0.d0
    do  jla = 1, ntrunc+1
      jn = jla-1
      zpole = zpole + zsp_mpiglobal(jla,1,1)*sqrt((2.d0*jn+1.d0)/2.d0)
    enddo
    if(zpole.le.0.d0) then
      call utl_abort('POLE VALUE NEGATIVE IN SUTG')
    endif
    do jla = 1, ntrunc+1
      zsp_mpiglobal(jla,1,1) = zsp_mpiglobal(jla,1,1)/zpole
      zsp_mpiglobal(jla,2,1) = zsp_mpiglobal(jla,2,1)/zpole
    enddo

    !  2.4.5  Correlation
    do jm = 0, ntrunc
      do jn = jm, ntrunc
        jla = gst_getNIND(jm,gstID) + jn - jm
        dlfac = 0.5d0/dsqrt((2*jn+1.d0)/2.d0)
        cortgg(jla,1) = dlfac * zsp_mpiglobal(jn+1,1,1)
        cortgg(jla,2) = dlfac * zsp_mpiglobal(jn+1,1,1)
      enddo
    enddo

    ! 2.5. For zonal modes : set to zero the imaginary part and set the correct factor 1.0 for the real part
    do jla = 1, ntrunc + 1
      cortgg(jla,1) = 0.5d0*cortgg(jla,1)
      cortgg(jla,2) = 0.0d0
    enddo

    ! 2.6. Result in corns array
    do jn = 0, ntrunc
      ila_mpiglobal = jn + 1
      corns(nspositTG,nspositTG,jn) = 2.d0*cortgg(ila_mpiglobal,1)
    enddo

    deallocate(dltg)
    !write(*,*)'DONE in SUTG'

  END SUBROUTINE BHI_sutg


  SUBROUTINE BHI_convol
    implicit none

    real(8) dlfact2,dlc,dsummed
    real(8) dtlen,zr,dlfact
    integer jn,jlat,jk
    real(8) zsp(0:ntrunc,nkgdim),zgr(nj_l,nkgdim)
    real(8) dlwti(nj_l),zrmu(nj_l)

    real(8)         :: RPORVO   = 6000.D3
    real(8)         :: RPORDI   = 6000.D3
    real(8)         :: RPORTT   = 3000.D3
    real(8)         :: RPORQ    = 3000.D3
    real(8)         :: RPORPS   = 3000.D3

    do jlat = 1, nj_l
       dlwti(jlat) = gst_getrwt(jlat,gstID)
       zrmu(jlat)  = gst_getrmu(jlat,gstID)
    end do

!     1.2 CONVERT THE CORRELATIONS IN SPECTRAL SPACE INTO SPECTRAL
!         COEFFICIENTS OF THE CORRELATION FUNCTION AND FUNCTION TO BE
!         SELF-CONVOLVED
    do jn = 0, ntrunc
      dlfact = ((2.0d0*jn +1.0d0)/2.0d0)**(0.25d0)
      dlfact2= ((2.0d0*jn +1.0d0)/2.0d0)**(0.25d0)
      do jk = 1, nkgdim
        zsp(jn,jk) = rstddev(jk,jn)*dlfact*dlfact2
      enddo
    enddo

    ! Transform to physical space
    call gst_zleginv(gstID,zgr,zsp,nkgdim)

    ! Truncate in horizontal extent with Gaussian window
    do jk = 1, nkgdim
      if (jk.ge.nspositVO.and.jk.lt.nspositVO+nlev_M) then
        dtlen = rporvo
      elseif (jk.ge.nspositDI.and.jk.lt.nspositDI+nlev_M) then
        dtlen = rpordi
      elseif (jk.ge.nspositTT.and.jk.lt.nspositTT+nlev_T) then
        dtlen = rportt
      elseif (jk.ge.nspositQ.and.jk.lt.nspositQ+nlev_T) then
        dtlen = rporq
      elseif (jk == nspositPS) then
        dtlen = rporps
      endif

      if(dtlen.gt.0.0d0) then
        dlc = 1.d0/dble(dtlen)
        dlc = 0.5d0*dlc*dlc
        do jlat = 1, nj_l
          zr = ec_ra * acos(zrmu(jlat))
          dlfact = dexp(-(zr**2)*dlc)
          zgr(jlat,jk) = dlfact*zgr(jlat,jk)
        enddo
      endif

      !write(*,*) 'zeroing length (km)=',jk,dtlen/1000.0
    enddo

    ! Transform back to spectral space
    call gst_zlegdir(gstID,zgr,zsp,nkgdim)

    ! Convert back to correlations
    do jk = 1, nkgdim
      do jn = 0, ntrunc
         zsp(jn,jk) = zsp(jn,jk)*(2.0d0/(2.0d0*jn+1.0d0))**(0.25d0)
      enddo
    enddo

    ! PUT BACK INTO RSTDDEV
    do jn = 0, ntrunc
      do jk = 1, nkgdim
         rstddev(jk,jn) = zsp(jn,jk)
      enddo
    enddo
 
    ! Re-normalize to ensure correlations
    do jk = 1, nkgdim
      dsummed = 0.d0
      do jn = 0, ntrunc
        dsummed = dsummed+ dble(rstddev(jk,jn)**2)*sqrt(((2.d0*jn)+1.d0)/2.d0)
      enddo
      dsummed = sqrt(dsummed)
      do jn = 0, ntrunc
        if(dsummed.gt.1.d-30) rstddev(jk,jn) = rstddev(jk,jn)/dsummed
      enddo
    enddo

    !     CONVERT THE SPECTRAL COEFFICIENTS OF THE CORRELATION FUNCTION
    !     .  BACK INTO CORRELATIONS OF SPECTRAL COMPONENTS
    do jn = 0, ntrunc
      dlfact = sqrt(0.5d0)*(1.0d0/((2.0d0*jn+1)/2.0d0))**0.25d0
      do jk = 1, nkgdim
        rstddev(jk,jn) = rstddev(jk,jn)*dlfact
      enddo
    enddo

  END SUBROUTINE BHI_convol


  SUBROUTINE BHI_setCrossCorr(kn)
    implicit none

    integer :: kn, jblock1, inbrblock, jblock2
    integer :: jk1, jk2, nlev_all(numvar3d), levOffset(numvar3d+1)

    inbrblock = numvar3d
    nlev_all(1) = nLev_M
    nlev_all(2) = nLev_M
    nlev_all(3) = nLev_T
    nlev_all(4) = nLev_T
    levOffset(1) = 0
    levOffset(2) = 1*nLev_M
    levOffset(3) = 2*nLev_M
    levOffset(4) = 2*nLev_M+1*nLev_T
    levOffset(5) = 2*nLev_M+2*nLev_T

    ! Set cross-variable correlations to 0 ...
    do jblock1 = 1, inbrblock
      do jblock2 = 1, inbrblock
        if (jblock1.ne.jblock2) then
          do jk2 = 1, nlev_all(jblock2)
            do jk1 = 1, nlev_all(jblock1)
              corns(jk1 + levOffset(jblock1),jk2 + levOffset(jblock2),kn) = 0.0d0
            enddo
          enddo
        endif
      enddo
    enddo

    ! ... but T'ln(ps') correlations
    do jk2 = 1, nkgdim
      do jk1 = levOffset(5)+1, levOffset(5)+numvar2d
        if ((jk1.ne.nspositPS.or.jk2.lt.nspositTT.or.    &
             jk2.ge.(nspositTT+nlev_T)).and.(jk1.ne.jk2)) then
          corns(jk1,jk2,kn) = 0.0d0
        endif
      enddo
    enddo

    do jk2 = levOffset(5)+1, levOffset(5)+numvar2d
      do jk1 = 1, nkgdim
        if ((jk2.ne.nspositPS.or.jk1.lt.nspositTT.or.   &
             jk1.ge.(nspositTT+nlev_T)) .and.(jk1.ne.jk2)) then
          corns(jk1,jk2,kn) = 0.0d0
        endif
      enddo
    enddo

  END SUBROUTINE BHI_setCrossCorr


  SUBROUTINE BHI_READCORNS2
    implicit none

    integer :: kip1
    integer :: jn, istdkey,icornskey
    integer :: iksdim,jcol,jrow
    real(8), allocatable, dimension(:) :: zstdsrc
    real(8), allocatable, dimension(:,:) :: zcornssrc

    ! standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3
    integer :: idateo
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket

    iksdim = 2*nlev_M+2*nlev_T+1    ! assume 4 3d variables and 1 2d variable (TG not included)
    allocate(zcornssrc(iksdim,iksdim))
    allocate(zstdsrc(iksdim))

    kip1 = -1

    do jn = 0, ntrunc

      ! Looking for FST record parameters..

      idateo = -1
      cletiket = 'RSTDDEV'
      ip1 = kip1
      ip2 = jn
      ip3 = -1
      cltypvar = 'X'
      clnomvar = 'SS'

      istdkey = utl_fstlir(ZSTDSRC,nulbgst,INI,INJ,INK,idateo,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)

      if(istdkey .lt.0 ) then
        call utl_abort('READCORNS2: Problem with background stat file')
      endif

      if (ini .ne. iksdim) then
        call utl_abort('READCORNS2: BG stat levels inconsitencies')
      endif

      ! Looking for FST record parameters..

      idateo = -1
      cletiket = 'CORRNS'
      ip1 = kip1
      IP2 = JN
      ip3 = -1
      cltypvar = 'X'
      clnomvar = 'ZZ'
      icornskey = utl_fstlir(ZCORNSSRC,nulbgst,INI,INJ,INK,idateo,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)

      if(icornskey .lt.0 ) then
        call utl_abort('READCORNS2: Problem with background stat file')
      endif

      if (ini .ne. iksdim .or. inj .ne. iksdim) then
        call utl_abort('READCORNS2: BG stat levels inconsitencies')
      endif

      do jcol = 1, nkgdim2
        rstddev(jcol,jn) = 0.0d0
        do jrow = 1, nkgdim2
          corns(jrow,jcol,jn) = 0.0d0
        enddo
      enddo

      do jcol = 1, iksdim
        do jrow = 1, iksdim
          corns(jrow,jcol,jn) = zcornssrc(jrow,jcol)
        enddo
      enddo

      ! Set cross-variable correlations to zero except between T' and ln(ps')
      call BHI_setcrosscorr(jn)

      do jrow = 1, iksdim
        rstddev(jrow,jn) = zstdsrc(jrow)
      enddo

    enddo

    ! Apply convolution to RSTDDEV correlations

    call BHI_convol 

    do jn = 0, ntrunc

      ! Re-build of correlation matrix: factorization of corns with convoluted RSTDDEV
      do jcol = 1, nkgdim
        do jrow = 1, nkgdim
          corns(jrow,jcol,jn) = rstddev(jrow,jn) * corns(jrow,jcol,jn)* rstddev(jcol,jn)
        enddo
      enddo

    enddo

    deallocate(zcornssrc)
    deallocate(zstdsrc)

    !write(*,*) 'Done in READCORNS2'
  END SUBROUTINE BHI_READCORNS2


  SUBROUTINE BHI_RDSPSTD
    implicit none

    integer, parameter  :: inbrvar3d=5
    integer, parameter  :: inbrvar2d=2
    integer :: jvar,jn,inix,injx,inkx
    integer :: ikey, jlevo, jlat,firstn,lastn
    real(8) :: zsp(0:ntrunc,max(nlev_M,nlev_T)),zspbuf(max(nlev_M,nlev_T))
    real(8) :: zgr(nj_l,max(nlev_M,nlev_T))
    character(len=4) :: varName3d(inbrvar3d),varName2d(inbrvar2d)

    ! standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3
    integer :: idate(100),nlev_MT
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fstinf

    data varName3d/'PP  ','UC  ','UT  ','LQ  ','TB  '/
    data varName2d/'UP  ','PB  '/

    rgsig(:,:) = 0.0d0
    rgsigtb(:,:) = 0.0d0
    rgsigpsb(:) = 0.0d0

!   2. Reading the data

    idate(1) = -1
    ip1      = -1
    ip2      = -1
    ip3      = -1

    cletiket = 'SPSTDDEV'
    cltypvar = 'X'

    do jvar = 1, inbrvar3d
      clnomvar = varName3d(jvar)
      if(vnl_varLevelFromVarName(clnomvar) == 'MM') then
        nlev_MT = nlev_M
      else
        nlev_MT = nlev_T
      endif
      firstn = -1
      do jn = 0, ntrunc
        ip2 = jn
        ikey = fstinf(nulbgst,inix,injx,inkx,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)

        if(ikey .ge.0 ) then
          ikey = utl_fstlir(zspbuf(1:nlev_MT),nulbgst,ini,inj,ink,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
        else
          if(firstn == -1) firstn = jn
          lastn = jn
          zspbuf(:) = 0.0d0
        endif

        if (ini .ne. nlev_MT) then
          call utl_abort('RDSPSTD: BG stat levels inconsitencies')
        endif

        do jlevo = 1, nlev_MT
          zsp(jn,jlevo) = zspbuf(jlevo)
        enddo
      enddo
      if(mmpi_myid == 0.and.firstn.ne.-1) then
        write(*,*) 'WARNING: CANNOT FIND SPSTD FOR ',clnomvar, &
                     ' AT N BETWEEN ',firstn,' AND ',lastn,', SETTING TO ZERO!!!'
      endif

      call gst_zleginv(gstID,zgr(:,1:nlev_MT),zsp(:,1:nlev_MT),nlev_MT)

      if(clnomvar == 'PP') then
        do jlat = 1, nj_l
          do jlevo = 1, nlev_M
            rgsiguu(jlat,jlevo) = zgr(jlat,jlevo)
          enddo
        enddo
      elseif(clnomvar == 'UC' .or. clnomvar == 'CC') then
        do jlat = 1, nj_l
          do jlevo = 1, nlev_M
            rgsigvv(jlat,jlevo) = zgr(jlat,jlevo)
          enddo
        enddo
      elseif(clnomvar == 'UT') then
        do jlat = 1, nj_l
          do jlevo = 1, nlev_T
            rgsigtt(jlat,jlevo) = zgr(jlat,jlevo)
          enddo
        enddo
      elseif(clnomvar == 'TB') then
        do jlat = 1, nj_l
          do jlevo = 1, nlev_T
            rgsigtb(jlat,jlevo) = zgr(jlat,jlevo)
          enddo
        enddo
      elseif(clnomvar == 'LQ') then
        do jlat = 1, nj_l
          do jlevo = 1, nlev_T
            rgsigq(jlat,jlevo) = zgr(jlat,jlevo)
          enddo
        enddo
      endif

    enddo

    nlev_MT = 1
    do jvar = 1, inbrvar2d
      clnomvar = varName2d(jvar)
      firstn = -1
      do jn = 0, ntrunc
        ip2 = jn
        ikey = fstinf(nulbgst,inix,injx,inkx,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)

        if(ikey .ge.0 ) then
          ikey = utl_fstlir(zspbuf,nulbgst,ini,inj,ink,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
        else
          if(firstn == -1) firstn = jn
          lastn = jn
          zspbuf(:) = 0.0d0
        endif

        zsp(jn,1) = zspbuf(1)

      enddo
      if(mmpi_myid == 0.and.firstn.ne.-1) then
        write(*,*) 'WARNING: CANNOT FIND SPSTD FOR ',clnomvar, &
                     ' AT N BETWEEN ',firstn,' AND ',lastn,', SETTING TO ZERO!!!'

      endif

      call gst_zleginv(gstID,zgr(:,1:nlev_MT),zsp(:,1:nlev_MT),nlev_MT)

      if(clnomvar == 'UP') then
        do jlat = 1, nj_l
          rgsigps(jlat) = zgr(jlat,1)*100.0d0
        enddo
      endif
      if(clnomvar == 'PB') then
        do jlat = 1, nj_l
          rgsigpsb(jlat) = zgr(jlat,1)*100.0d0
        enddo
      endif

    enddo

  END SUBROUTINE BHI_RDSPSTD


  SUBROUTINE BHI_RDSTD
    implicit none

    integer, parameter  :: inbrvar3d=5
    integer, parameter  :: inbrvar2d=2
    integer :: jvar,inix,injx,inkx
    integer :: ikey
    real(8) :: zgr(nj_l,max(nlev_M,nlev_T))
    character(len=4) :: varName3d(inbrvar3d),varName2d(inbrvar2d)

    ! standard file variables
    integer :: ini,inj,ink  
    integer :: ip1,ip2,ip3
    integer :: idate(100),nlev_MT
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fstinf

    data varName3d/'PP  ','UC  ','UT  ','LQ  ','TB  '/
    data varName2d/'UP  ','PB  '/

    rgsig(:,:) = 0.0d0
    rgsigtb(:,:) = 0.0d0
    rgsigpsb(:) = 0.0d0

!   2. Reading the data

    idate(1) = -1
    ip1      = -1
    ip2      = -1
    ip3      = -1

    cletiket = 'STDDEV'
    cltypvar = 'E'

    do jvar = 1, inbrvar3d
      clnomvar = varName3d(jvar)
      if(vnl_varLevelFromVarName(clnomvar) == 'MM') then
        nlev_MT = nlev_M
      else
        nlev_MT = nlev_T
      endif

      ikey = fstinf(nulbgst,inix,injx,inkx,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)

      if(ikey .ge.0 ) then
        ikey = utl_fstlir(zgr(:,1:nlev_MT),nulbgst,ini,inj,ink,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
      else
        write(*,*) 'RDSTD: could not read varName=',clnomvar
        call utl_abort('RDSTD') 
      endif

      if (ink .ne. nlev_MT) then
        write(*,*) 'ink, nlev_MT=', ink, nlev_MT
        call utl_abort('RDSPSTD: BG stat levels inconsitencies')
      endif

      if(clnomvar == 'PP') then
        rgsiguu(:,:) = zgr(:,1:nlev_M)
      elseif(clnomvar == 'UC' .or. clnomvar == 'CC') then
        rgsigvv(:,:) = zgr(:,1:nlev_M)
      elseif(clnomvar == 'UT') then
        rgsigtt(:,:) = zgr(:,1:nlev_T)
      elseif(clnomvar == 'TB') then
        rgsigtb(:,:) = zgr(:,1:nlev_T)
      elseif(clnomvar == 'LQ') then
        rgsigq(:,:) = zgr(:,1:nlev_T)
      endif

    enddo

    nlev_MT = 1
    do jvar = 1, inbrvar2d
      clnomvar = varName2d(jvar)

      ikey = fstinf(nulbgst,inix,injx,inkx,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)

      if(ikey .ge.0 ) then
        ikey = utl_fstlir(zgr,nulbgst,ini,inj,ink,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
      else
        write(*,*) 'RDSTD: could not read varName=',clnomvar
        call utl_abort('RDSTD') 
      endif

      if(clnomvar == 'UP') then
        rgsigps(:) = zgr(:,1)*100.0d0
      endif
      if(clnomvar == 'PB') then
        rgsigpsb(:) = zgr(:,1)*100.0d0
      endif

    enddo

  END SUBROUTINE BHI_RDSTD


  SUBROUTINE BHI_RDSPSTD_NEWFMT
    implicit none

    integer, parameter  :: inbrvar3d=5
    integer, parameter  :: inbrvar2d=2
    integer :: jvar,jn,inix,injx,inkx,ntrunc_file
    integer :: ikey,jlevo,jlat
    real(8) :: zsp(0:ntrunc,max(nlev_M,nlev_T))
    real(8), allocatable :: zspbuf(:)
    real(8) :: zgr(nj_l,max(nlev_M,nlev_T))
    character(len=4) :: varName3d(inbrvar3d),varName2d(inbrvar2d)

    ! standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3
    integer :: idate(100),nlev_MT
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fstinf

    data varName3d/'PP  ','UC  ','UT  ','LQ  ','TB  '/
    data varName2d/'UP  ','PB  '/

    rgsig(:,:) = 0.0d0
    rgsigtb(:,:) = 0.0d0
    rgsigpsb(:) = 0.0d0

!   2. Reading the data

    idate(1) = -1
    ip2      = -1
    ip3      = -1

    cletiket = 'SPSTDDEV'
    cltypvar = 'X'

    ! check if file is old format
    ip1 = -1
    clnomvar = varName3d(1)
    ikey = fstinf(nulbgst,inix,injx,inkx,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
    write(*,*) 'ini,inj,ink=',inix,injx,inkx
    if(inix.gt.1) then
      write(*,*) 'BHI_RDSPSTD_NEWFMT: ini>1, SPSTDDEV is in old format, calling BHI_RDSPSTD...'
      call bhi_rdspstd
      return
    endif

    !write(*,*) 'Reading 3D variables'
    do jvar = 1, inbrvar3d
      clnomvar = varName3d(jvar)
      if(vnl_varLevelFromVarName(clnomvar) == 'MM') then
        nlev_MT = nlev_M
      else
        nlev_MT = nlev_T
      endif
      !write(*,*)'Reading ',clnomvar
      do jlevo = 1, nlev_MT
        if(vnl_varLevelFromVarName(clnomvar) == 'MM') then
          ip1 = vco_anl%ip1_M(jlevo)
        else
          ip1 = vco_anl%ip1_T(jlevo)
        endif
        ikey = fstinf(nulbgst,inix,ntrunc_file,inkx,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
        ntrunc_file = ntrunc_file-1

        allocate(zspbuf(0:ntrunc_file))
        if(ikey .ge.0 ) then
          ikey = utl_fstlir(zspbuf(0:ntrunc_file),nulbgst,ini,inj,ink,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
        else
          write(*,*) 'RDSPSTD_NEWFMT: ',jvar,clnomvar,nlev_MT,jlevo,ikey,ntrunc,ntrunc_file
          call utl_abort('RDSPSTD_NEWFMT: SPSTDDEV record not found')
        endif

        zsp(:,jlevo) = 0.0d0
        do jn = 0, min(ntrunc,ntrunc_file)
          zsp(jn,jlevo) = zspbuf(jn)
        enddo
        deallocate(zspbuf)

      enddo

      call gst_zleginv(gstID,zgr(:,1:nlev_MT),zsp(:,1:nlev_MT),nlev_MT)

      if(clnomvar == 'PP') then
        do jlat = 1, nj_l
          do jlevo = 1, nlev_M
            rgsiguu(jlat,jlevo) = zgr(jlat,jlevo)
          enddo
        enddo
      elseif(clnomvar == 'UC' .or. clnomvar == 'CC') then
        do jlat = 1, nj_l
          do jlevo = 1, nlev_M
            rgsigvv(jlat,jlevo) = zgr(jlat,jlevo)
          enddo
        enddo
      elseif(clnomvar == 'UT') then
        do jlat = 1, nj_l
          do jlevo = 1, nlev_T
            rgsigtt(jlat,jlevo) = zgr(jlat,jlevo)
          enddo
        enddo
      elseif(clnomvar == 'TB') then
        do jlat = 1, nj_l
          do jlevo = 1, nlev_T
            rgsigtb(jlat,jlevo) = zgr(jlat,jlevo)
          enddo
        enddo
      elseif(clnomvar == 'LQ') then
        do jlat = 1, nj_l
          do jlevo = 1, nlev_T
            rgsigq(jlat,jlevo) = zgr(jlat,jlevo)
          enddo
        enddo
      endif

    enddo

    nlev_MT = 1
    do jvar = 1, inbrvar2d
      clnomvar = varName2d(jvar)
      ip1 = -1
      ikey = fstinf(nulbgst,inix,ntrunc_file,inkx,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
      ntrunc_file = ntrunc_file-1

      allocate(zspbuf(0:ntrunc_file))

      if(ikey .ge.0 ) then
        ikey = utl_fstlir(zspbuf(0:ntrunc_file),nulbgst,ini,inj,ink,idate(1),cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
      else
        write(*,*) 'RDSPSTD_NEWFMT: ',jvar,clnomvar,nlev_MT,jlevo,ikey,ntrunc,ntrunc_file
        call utl_abort('RDSPSTD_NEWFMT: SPSTDDEV record not found')
      endif

      zsp(:,1) = 0.0d0
      do jn = 0, min(ntrunc,ntrunc_file)
        zsp(jn,1) = zspbuf(jn)
      enddo
      deallocate(zspbuf)

      call gst_zleginv(gstID,zgr(:,1:nlev_MT),zsp(:,1:nlev_MT),nlev_MT)

      if(clnomvar == 'UP') then
        do jlat = 1, nj_l
          rgsigps(jlat) = zgr(jlat,1)*100.0d0
        enddo
      endif
      if(clnomvar == 'PB') then
        do jlat = 1, nj_l
          rgsigpsb(jlat) = zgr(jlat,1)*100.0d0
        enddo
      endif

    enddo

  END SUBROUTINE BHI_RDSPSTD_NEWFMT


  SUBROUTINE BHI_RDSPPTOT
    IMPLICIT NONE

    integer :: jn, jk1, jk2, ikey, ilen,jlat,inix,injx,inkx
    real(8) :: zsptheta(0:ntrunc,nlev_M)
    real(8) :: zgrtheta(nj_l,nlev_M)
    real(8) :: zPtoTsrc(nlev_T+1,nlev_M)
    real(8) :: zspPtoT(0:ntrunc,nlev_T+1,nlev_M)
    real(8) :: zgrPtoT(nj_l,nlev_T+1,nlev_M)
    real(8) :: ztheta(nlev_M)
    
    ! standard file variables
    integer :: ini,inj,ink
    integer :: ip1,ip2,ip3
    integer :: idateo
    character(len=2)  :: cltypvar
    character(len=4)  :: clnomvar
    character(len=12) :: cletiket
    integer :: fstinf
    
    ip1 = -1
    ip3 = -1
    idateo = -1
    cletiket = 'SP_THETA'
    cltypvar = 'X'
    clnomvar = 'ZZ'

    ! read spectral coefficients for theta

    do jn = 0, ntrunc
      ip2 = jn
      ikey = fstinf(nulbgst,inix,injx,inkx,idateo,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)

      if(ikey .ge.0 ) then
        ikey = utl_fstlir(ztheta,nulbgst,ini,inj,ink,idateo,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
      else
        if(mmpi_myid == 0) write(*,*) 'WARNING: CANNOT FIND THETA FOR ',jn,' SETTING TO ZERO!!!'
        ztheta(:) = 0.0d0
      endif

      do jk1 = 1, nlev_M
        zsptheta(jn,jk1) = ztheta(jk1)
      enddo

    enddo

    ! converting theta in physical space

    call gst_zleginv(gstID,zgrtheta,zsptheta,nlev_M)

    do jlat = 1, nj_l
      do jk1 = 1, nlev_M
        tantheta(jk1,jlat) = tan(zgrtheta(jlat,jk1))
      end do
    end do

    ip1 = -1
    ip2 = -1
    ip3 = -1
    idateo = -1
    cletiket = 'SP_PtoT'
    cltypvar = 'X'
    clnomvar = 'ZZ'

    ! read of spectral coefficients for P to T operator

    do jn = 0, ntrunc
      ip2 = jn
      ikey = fstinf(nulbgst,inix,injx,inkx,idateo,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)

      if(ikey .ge.0 ) then
        ikey = utl_fstlir(zPtoTsrc,nulbgst,ini,inj,ink,idateo,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
      else
        if(mmpi_myid == 0) write(*,*) 'WARNING: CANNOT FIND P_to_T FOR ',jn,' SETTING TO ZERO!!!'
        zPtoTsrc(:,:) = 0.0d0
      endif

      do jk2 = 1, nlev_M
        do jk1 = 1, nlev_T+1
          zspPtoT(jn,jk1,jk2) = zPtoTsrc(jk1,jk2)
        enddo
      enddo

    enddo

    ilen = nlev_M*(nlev_T+1)
    call gst_zleginv(gstID,zgrPtoT,zspPtoT,ilen)

    do jlat = 1, nj_l
      do jk2 = 1, nlev_M
        do jk1 = 1, nlev_T+1
          PtoT(jk1,jk2,jlat) = zgrPtoT(jlat,jk1,jk2)
        end do
      end do
    enddo

  END SUBROUTINE BHI_RDSPPTOT


  SUBROUTINE BHI_truncateCV(controlVector_inout,ntruncCut)
    implicit none
    ! set to zero all coefficients with total wavenumber > ntruncCut

    real(8), pointer :: controlVector_inout(:)
    integer          :: ntruncCut
    integer          :: jn, jm, ila_mpiglobal, ila_mpilocal, jlev, jdim

    if(.not. initialized) then
      if(mmpi_myid == 0) write(*,*) 'bhi_truncateCV: bMatrixHI not initialized'
      return
    endif

    if(ntruncCut.gt.ntrunc) then
      write(*,*) ntruncCut, ntrunc
      call utl_abort('bhi_truncateCV: ntruncCut is greater than ntrunc!')
    endif

    jdim = 0
    do jlev = 1, nkgdimSqrt
      do jm = mymBeg, mymEnd, mymSkip
        do jn = mynBeg, mynEnd, mynSkip
          if(jm.le.jn) then
            ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
            ila_mpilocal  = ilaList_mpilocal(ila_mpiglobal)
            if(jm == 0) then
              ! only real component for jm=0
              jdim = jdim + 1
              if(jn.gt.ntruncCut) controlVector_inout(jdim) = 0.0d0
            else
              ! both real and imaginary components for jm>0
              jdim = jdim + 1
              if(jn.gt.ntruncCut) controlVector_inout(jdim) = 0.0d0
              jdim = jdim + 1
              if(jn.gt.ntruncCut) controlVector_inout(jdim) = 0.0d0
            endif
          endif
        enddo
      enddo
    enddo

  END SUBROUTINE BHI_truncateCV


  SUBROUTINE BHI_bSqrt(controlVector_in, statevector, stateVectorRef_opt)
    implicit none

    ! Arguments
    real(8)                    :: controlVector_in(cvDim_mpilocal)
    type(struct_gsv)           :: statevector
    type(struct_gsv), optional :: statevectorRef_opt

    ! Locals
    real(8),allocatable :: gd_out(:,:,:)
    real(8)   :: hiControlVector(nla_mpilocal,2,nkgdimSqrt)

    if(mmpi_myid == 0) write(*,*) 'bhi_bsqrt: starting'
    if(mmpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    if(.not. initialized) then
      if(mmpi_myid == 0) write(*,*) 'bMatrixHI not initialized'
      return
    endif

    allocate(gd_out(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim))

    call bhi_cain(controlVector_in,hiControlVector)
    call bhi_spa2gd(hiControlVector,gd_out)

    call copyToStatevector(statevector,gd_out)

    if (gsv_varExist(statevector,'HU')) then
      call gvt_transform( statevector,  &                         ! INOUT
                          'LQtoHU_tlm', &                         ! IN
                          stateVectorRef_opt=stateVectorRef_opt ) ! IN
    end if

    deallocate(gd_out)

    if(mmpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if(mmpi_myid == 0) write(*,*) 'bhi_bsqrt: done'

  END SUBROUTINE BHI_bSqrt


  SUBROUTINE BHI_bSqrtAd(statevector, controlVector_out, stateVectorRef_opt)
    implicit none

    ! Arguments
    real(8)                    :: controlVector_out(cvDim_mpilocal)
    type(struct_gsv)           :: statevector
    type(struct_gsv), optional :: statevectorRef_opt

    ! Locals
    real(8), allocatable :: gd_in(:,:,:)
    real(8)   :: hiControlVector(nla_mpilocal,2,nkgdimSqrt)

    if(.not. initialized) then
      if(mmpi_myid == 0) write(*,*) 'bMatrixHI not initialized'
      return
    endif

    if(mmpi_myid == 0) write(*,*) 'bhi_bsqrtad: starting'
    if(mmpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'

    allocate(gd_in(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim))
    gd_in(:,:,:) = 0.d0
    if (gsv_varExist(statevector,'HU')) then
      call gvt_transform( statevector, &                          ! INOUT
                          'LQtoHU_ad', &                          ! IN
                          stateVectorRef_opt=stateVectorRef_opt ) ! IN
    end if

    call copyFromStatevector(statevector,gd_in)

    call bhi_spa2gdad(gd_in,hiControlVector)

    call bhi_cainad(hiControlVector,controlVector_out)

    deallocate(gd_in)

    if(mmpi_myid == 0) write(*,*) 'Memory Used: ',get_max_rss()/1024,'Mb'
    if(mmpi_myid == 0) write(*,*) 'bhi_bsqrtad: done'

  END SUBROUTINE BHI_bSqrtAd


  SUBROUTINE copyToStatevector(statevector,gd)
    implicit none
    type(struct_gsv) :: statevector
    real(8) :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)
    integer :: jlon, jlev, jlev2, jlat, jvar, ilev1, ilev2
    real(8), pointer :: field_r8(:,:,:)
    real(4), pointer :: field_r4(:,:,:)

    do jvar = 1, vnl_numvarmax 
      if(gsv_varExist(statevector,vnl_varNameList(jvar))) then
        if (gsv_getDataKind(statevector) == 8) then
          call gsv_getField(statevector,field_r8,vnl_varNameList(jvar))
        else
          call gsv_getField(statevector,field_r4,vnl_varNameList(jvar))
        end if
        if(vnl_varNameList(jvar) == 'UU  ') then
          ilev1 = nspositVO
        elseif(vnl_varNameList(jvar) == 'VV  ') then
          ilev1 = nspositDI
        elseif(vnl_varNameList(jvar) == 'TT  ') then
          ilev1 = nspositTT
        elseif(vnl_varNameList(jvar) == 'HU  ' .or. vnl_varNameList(jvar) == 'LQ  ') then
          ilev1 = nspositQ
        elseif(vnl_varNameList(jvar) == 'P0  ') then
          ilev1 = nspositPS
        elseif(vnl_varNameList(jvar) == 'TG  ') then
          ilev1 = nspositTG
        else
          ! Cycle (instead of abort) to allow for non-NWP assimilation (e.g. chemical data assimilation)
!          call utl_abort('bmatrixhi_mod: copyToStatevector: No covariances available for variable:' // vnl_varNameList(jvar))
          cycle
        endif
        ilev2 = ilev1 - 1 + gsv_getNumLev(statevector,vnl_varLevelFromVarname(vnl_varNameList(jvar)))
        if (gsv_getDataKind(statevector) == 8) then
          do jlev = ilev1, ilev2
            jlev2 = jlev-ilev1+1
            do jlat = myLatBeg, myLatEnd
              do jlon = myLonBeg, myLonEnd
                field_r8(jlon,jlat,jlev2) = gd(jlon,jlat,jlev)
              enddo
            enddo
          enddo
        else
          do jlev = ilev1, ilev2
            jlev2 = jlev-ilev1+1
            do jlat = myLatBeg, myLatEnd
              do jlon = myLonBeg, myLonEnd
                field_r4(jlon,jlat,jlev2) = gd(jlon,jlat,jlev)
              enddo
            enddo
          enddo
        end if
      endif
    enddo

  END SUBROUTINE copyToStatevector


  SUBROUTINE copyFromStatevector(statevector,gd)
    implicit none
    type(struct_gsv) :: statevector
    real(8)          :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)
    integer :: jlon, jlev, jlev2, jlat, jvar, ilev1, ilev2
    real(8), pointer :: field_r8(:,:,:)
    real(4), pointer :: field_r4(:,:,:)

    do jvar = 1, vnl_numvarmax 
      if(gsv_varExist(statevector,vnl_varNameList(jvar))) then
        if (gsv_getDataKind(statevector) == 8) then
          call gsv_getField(statevector,field_r8,vnl_varNameList(jvar))
        else
          call gsv_getField(statevector,field_r4,vnl_varNameList(jvar))
        end if
        if(vnl_varNameList(jvar) == 'UU  ') then
          ilev1 = nspositVO
        elseif(vnl_varNameList(jvar) == 'VV  ') then
          ilev1 = nspositDI
        elseif(vnl_varNameList(jvar) == 'TT  ') then
          ilev1 = nspositTT
        elseif(vnl_varNameList(jvar) == 'HU  ') then
          ilev1 = nspositQ
        elseif(vnl_varNameList(jvar) == 'P0  ') then
          ilev1 = nspositPS
        elseif(vnl_varNameList(jvar) == 'TG  ') then
          ilev1 = nspositTG
        else
          ! Cycle (instead of abort) to allow for non-NWP assimilation (e.g. chemical data assimilation)
          cycle
        endif
        ilev2 = ilev1 - 1 + gsv_getNumLev(statevector,vnl_varLevelFromVarname(vnl_varNameList(jvar)))
        if (gsv_getDataKind(statevector) == 8) then
          do jlev = ilev1, ilev2
            jlev2 = jlev-ilev1+1
            do jlat = myLatBeg, myLatEnd
              do jlon = myLonBeg, myLonEnd
                gd(jlon,jlat,jlev) = field_r8(jlon,jlat,jlev2)
              enddo
            enddo
          enddo
        else
          do jlev = ilev1, ilev2
            jlev2 = jlev-ilev1+1
            do jlat = myLatBeg, myLatEnd
              do jlon = myLonBeg, myLonEnd
                gd(jlon,jlat,jlev) = field_r4(jlon,jlat,jlev2)
              enddo
            enddo
          enddo
        end if
      endif
    enddo

  END SUBROUTINE copyFromStatevector


  SUBROUTINE BHI_reduceToMPILocal(cv_mpilocal,cv_mpiglobal)
    implicit none
    real(8), intent(out) :: cv_mpilocal(cvDim_mpilocal)
    real(8), intent(in)  :: cv_mpiglobal(:)

    real(8), allocatable :: cv_allMaxMpilocal(:,:)

    integer, allocatable :: cvDim_allMpilocal(:), displs(:)
    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)

    integer :: jproc,cvDim_maxmpilocal,ierr
    integer :: jlev,jn,jm,ila_mpiglobal,jdim_mpilocal,jdim_mpiglobal

    call rpn_comm_allreduce(cvDim_mpilocal, cvDim_maxmpilocal, &
                            1,"MPI_INTEGER","MPI_MAX","GRID",ierr)

    if(mmpi_myid == 0) then
       allocate(cvDim_allMpiLocal(mmpi_nprocs))
    else
       allocate(cvDim_allMpiLocal(1))
    end if

    call rpn_comm_gather(cvDim_mpiLocal   ,1,"mpi_integer",       &
                         cvDim_allMpiLocal,1,"mpi_integer",0,"GRID",ierr)

    if(mmpi_myid == 0) then
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

    ! Prepare to data to be distributed
    if (mmpi_myid == 0) then

       allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mmpi_nprocs))

       !$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,jlev,jm,jn,ila_mpiglobal,jdim_mpiglobal)
       do jproc = 0, (mmpi_nprocs-1)
          cv_allmaxmpilocal(:,jproc+1) = 0.d0
          jdim_mpilocal = 0

          do jlev = 1, nkgdimSqrt
             do jm = allmBeg(jproc+1), allmEnd(jproc+1), allmSkip(jproc+1)
                do jn = allnBeg(jproc+1), allnEnd(jproc+1), allnSkip(jproc+1)

                   if(jm.le.jn) then
                      
                      ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
                      
                      ! figure out index into global control vector
                      if(jm == 0) then
                         ! for jm=0 only real part
                         jdim_mpiglobal = ila_mpiglobal
                      else
                         ! for jm>0 both real and imaginary part
                         jdim_mpiglobal = 2*ila_mpiglobal-1 - (ntrunc+1)
                      endif
                      ! add offset for level
                      jdim_mpiglobal = jdim_mpiglobal + (jlev-1) * (ntrunc+1)*(ntrunc+1)
                      
                      ! index into local control vector computer as in cain
                      if(jm == 0) then
                         ! only real component for jm=0
                         jdim_mpilocal = jdim_mpilocal + 1
                         cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal)
                      else
                         ! both real and imaginary components for jm>0
                         jdim_mpilocal = jdim_mpilocal + 1
                         cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal)
                         jdim_mpilocal = jdim_mpilocal + 1
                         cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal+1)
                      endif
                      
                      if (jdim_mpilocal > cvDim_allMpiLocal(jproc+1)) then
                         write(*,*)
                         write(*,*) 'ERROR: jdim_mpilocal > cvDim_allMpiLocal(jproc+1)', jdim_mpilocal, cvDim_mpilocal
                         write(*,*) '       proc, jlev, jn, jm = ',jproc, jlev, jn, jm
                         call utl_abort('bhi_reduceToMPILocal')
                      end if
                      if (jdim_mpiglobal > cvDim_mpiglobal) then
                         write(*,*)
                         write(*,*) 'ERROR: jdim_mpiglobal > cvDim_mpiglobal', jdim_mpiglobal, cvDim_mpiglobal
                         write(*,*) '       proc, jlev, jn, jm = ',jproc, jlev, jn, jm
                         call utl_abort('bhi_reduceToMPILocal')
                      end if

                   endif
                enddo
             enddo
          enddo
 
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

  END SUBROUTINE BHI_reduceToMPILocal


  SUBROUTINE BHI_reduceToMPILocal_r4(cv_mpilocal,cv_mpiglobal)
    implicit none
    real(4), intent(out) :: cv_mpilocal(cvDim_mpilocal)
    real(4), intent(in)  :: cv_mpiglobal(:)

    real(4), allocatable :: cv_allMaxMpilocal(:,:)

    integer, allocatable :: cvDim_allMpilocal(:), displs(:)
    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)

    integer :: jproc,cvDim_maxmpilocal,ierr
    integer :: jlev,jn,jm,ila_mpiglobal,jdim_mpilocal,jdim_mpiglobal

    call rpn_comm_allreduce(cvDim_mpilocal, cvDim_maxmpilocal, &
                            1,"MPI_INTEGER","MPI_MAX","GRID",ierr)

    if(mmpi_myid == 0) then
       allocate(cvDim_allMpiLocal(mmpi_nprocs))
    else
       allocate(cvDim_allMpiLocal(1))
    end if

    call rpn_comm_gather(cvDim_mpiLocal   ,1,"mpi_integer",       &
                         cvDim_allMpiLocal,1,"mpi_integer",0,"GRID",ierr)

    if(mmpi_myid == 0) then
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

    ! Prepare to data to be distributed
    if (mmpi_myid == 0) then

       allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mmpi_nprocs))

       !$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,jlev,jm,jn,ila_mpiglobal,jdim_mpiglobal)
       do jproc = 0, (mmpi_nprocs-1)
          cv_allmaxmpilocal(:,jproc+1) = 0.d0
          jdim_mpilocal = 0

          do jlev = 1, nkgdimSqrt
             do jm = allmBeg(jproc+1), allmEnd(jproc+1), allmSkip(jproc+1)
                do jn = allnBeg(jproc+1), allnEnd(jproc+1), allnSkip(jproc+1)

                   if(jm.le.jn) then
                      
                      ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
                      
                      ! figure out index into global control vector
                      if(jm == 0) then
                         ! for jm=0 only real part
                         jdim_mpiglobal = ila_mpiglobal
                      else
                         ! for jm>0 both real and imaginary part
                         jdim_mpiglobal = 2*ila_mpiglobal-1 - (ntrunc+1)
                      endif
                      ! add offset for level
                      jdim_mpiglobal = jdim_mpiglobal + (jlev-1) * (ntrunc+1)*(ntrunc+1)
                      
                      ! index into local control vector computer as in cain
                      if(jm == 0) then
                         ! only real component for jm=0
                         jdim_mpilocal = jdim_mpilocal + 1
                         cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal)
                      else
                         ! both real and imaginary components for jm>0
                         jdim_mpilocal = jdim_mpilocal + 1
                         cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal)
                         jdim_mpilocal = jdim_mpilocal + 1
                         cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal+1)
                      endif
                      
                      if (jdim_mpilocal > cvDim_allMpiLocal(jproc+1)) then
                         write(*,*)
                         write(*,*) 'ERROR: jdim_mpilocal > cvDim_allMpiLocal(jproc+1)', jdim_mpilocal, cvDim_mpilocal
                         write(*,*) '       proc, jlev, jn, jm = ',jproc, jlev, jn, jm
                         call utl_abort('bhi_reduceToMPILocal')
                      end if
                      if (jdim_mpiglobal > cvDim_mpiglobal) then
                         write(*,*)
                         write(*,*) 'ERROR: jdim_mpiglobal > cvDim_mpiglobal', jdim_mpiglobal, cvDim_mpiglobal
                         write(*,*) '       proc, jlev, jn, jm = ',jproc, jlev, jn, jm
                         call utl_abort('bhi_reduceToMPILocal')
                      end if

                   endif
                enddo
             enddo
          enddo
 
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

  END SUBROUTINE BHI_reduceToMPILocal_r4


  SUBROUTINE BHI_expandToMPIGlobal(cv_mpilocal,cv_mpiglobal)
    implicit none
    real(8), intent(in)  :: cv_mpilocal(cvDim_mpilocal)
    real(8), intent(out) :: cv_mpiglobal(:)

    real(8), allocatable :: cv_maxmpilocal(:)
    real(8), pointer :: cv_allmaxmpilocal(:,:) => null()
    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)
    integer :: jlev, jn, jm, jproc, ila_mpiglobal, jdim_mpilocal, jdim_mpiglobal, ierr, cvDim_maxmpilocal

    !
    !- 1.  Gather all local control vectors onto mpi task 0
    !
    call rpn_comm_allreduce(cvDim_mpilocal,cvDim_maxmpilocal,1,"mpi_integer","mpi_max","GRID",ierr)

    allocate(cv_maxmpilocal(cvDim_maxmpilocal))

    nullify(cv_allmaxmpilocal)
    if(mmpi_myid == 0) then
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
    if(mmpi_myid == 0) then
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

    if(mmpi_myid == 0) then
      cv_mpiglobal(:) = 0.0d0

      !$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,jlev,jm,jn,ila_mpiglobal,jdim_mpiglobal)
      do jproc = 0, (mmpi_nprocs-1)
        jdim_mpilocal = 0

        do jlev = 1, nkgdimSqrt
          do jm = allmBeg(jproc+1), allmEnd(jproc+1), allmSkip(jproc+1)
            do jn = allnBeg(jproc+1), allnEnd(jproc+1), allnSkip(jproc+1)
              if(jm.le.jn) then

                ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm

                ! figure out index into global control vector
                if(jm == 0) then
                  ! for jm=0 only real part
                  jdim_mpiglobal = ila_mpiglobal
                else
                  ! for jm>0 both real and imaginary part
                  jdim_mpiglobal = 2*ila_mpiglobal-1 - (ntrunc+1)
                endif
                ! add offset for level
                jdim_mpiglobal = jdim_mpiglobal + (jlev-1) * (ntrunc+1)*(ntrunc+1)

                ! index into local control vector
                if(jm == 0) then
                  ! only real component for jm=0
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_mpiglobal(jdim_mpiglobal) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                else
                  ! both real and imaginary components for jm>0
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_mpiglobal(jdim_mpiglobal) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_mpiglobal(jdim_mpiglobal+1) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                endif

                if(jdim_mpiglobal.gt.cvDim_mpiglobal)   &
                  write(*,*) 'ERROR: jdim,cvDim,mpiglobal=',jdim_mpiglobal,cvDim_mpiglobal,jlev,jn,jm

              endif
            enddo
          enddo
        enddo
      enddo ! jproc
      !$OMP END PARALLEL DO

    endif ! myid == 0 

    deallocate(allnBeg)
    deallocate(allnEnd)
    deallocate(allnSkip)
    deallocate(allmBeg)
    deallocate(allmEnd)
    deallocate(allmSkip)
    deallocate(cv_allmaxmpilocal)

  end SUBROUTINE BHI_expandToMPIGlobal


  SUBROUTINE BHI_expandToMPIGlobal_r4(cv_mpilocal,cv_mpiglobal)
    implicit none
    real(4), intent(in)  :: cv_mpilocal(cvDim_mpilocal)
    real(4), intent(out) :: cv_mpiglobal(:)

    real(4), allocatable :: cv_maxmpilocal(:)
    real(4), pointer :: cv_allmaxmpilocal(:,:) => null()
    integer, allocatable :: allnBeg(:),allnEnd(:),allnSkip(:)
    integer, allocatable :: allmBeg(:),allmEnd(:),allmSkip(:)
    integer :: jlev, jn, jm, jproc, ila_mpiglobal, jdim_mpilocal, jdim_mpiglobal, ierr, cvDim_maxmpilocal

    !
    !- 1.  Gather all local control vectors onto mpi task 0
    !
    call rpn_comm_allreduce(cvDim_mpilocal,cvDim_maxmpilocal,1,"mpi_integer","mpi_max","GRID",ierr)

    allocate(cv_maxmpilocal(cvDim_maxmpilocal))

    nullify(cv_allmaxmpilocal)
    if(mmpi_myid == 0) then
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
    if(mmpi_myid == 0) then
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

    if(mmpi_myid == 0) then
      cv_mpiglobal(:) = 0.0d0

      !$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,jlev,jm,jn,ila_mpiglobal,jdim_mpiglobal)
      do jproc = 0, (mmpi_nprocs-1)
        jdim_mpilocal = 0

        do jlev = 1, nkgdimSqrt
          do jm = allmBeg(jproc+1), allmEnd(jproc+1), allmSkip(jproc+1)
            do jn = allnBeg(jproc+1), allnEnd(jproc+1), allnSkip(jproc+1)
              if(jm.le.jn) then

                ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm

                ! figure out index into global control vector
                if(jm == 0) then
                  ! for jm=0 only real part
                  jdim_mpiglobal = ila_mpiglobal
                else
                  ! for jm>0 both real and imaginary part
                  jdim_mpiglobal = 2*ila_mpiglobal-1 - (ntrunc+1)
                endif
                ! add offset for level
                jdim_mpiglobal = jdim_mpiglobal + (jlev-1) * (ntrunc+1)*(ntrunc+1)

                ! index into local control vector
                if(jm == 0) then
                  ! only real component for jm=0
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_mpiglobal(jdim_mpiglobal) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                else
                  ! both real and imaginary components for jm>0
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_mpiglobal(jdim_mpiglobal) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                  jdim_mpilocal = jdim_mpilocal + 1
                  cv_mpiglobal(jdim_mpiglobal+1) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)
                endif

                if(jdim_mpiglobal.gt.cvDim_mpiglobal)   &
                  write(*,*) 'ERROR: jdim,cvDim,mpiglobal=',jdim_mpiglobal,cvDim_mpiglobal,jlev,jn,jm

              endif
            enddo
          enddo
        enddo
      enddo ! jproc
      !$OMP END PARALLEL DO

    endif ! myid == 0 

    deallocate(allnBeg)
    deallocate(allnEnd)
    deallocate(allnSkip)
    deallocate(allmBeg)
    deallocate(allmEnd)
    deallocate(allmSkip)
    deallocate(cv_allmaxmpilocal)

  end SUBROUTINE BHI_expandToMPIGlobal_r4


  SUBROUTINE BHI_cain(controlVector_in,hiControlVector_out)
    implicit none

    real(8) :: controlVector_in(cvDim_mpilocal)
    real(8) :: hiControlVector_out(nla_mpilocal,2,nkgdimSqrt)

    integer :: jdim, jlev, jm, jn, ila_mpilocal, ila_mpiglobal

    jdim = 0
    hiControlVector_out(:,:,:) = 0.0d0
    do jlev = 1, nkgdimSqrt
      do jm = mymBeg, mymEnd, mymSkip
        do jn = mynBeg, mynEnd, mynSkip
          if(jm.le.jn) then
            ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
            ila_mpilocal  = ilaList_mpilocal(ila_mpiglobal)
            if(jm == 0) then
              ! only real component for jm=0
              jdim = jdim + 1
              hiControlVector_out(ila_mpilocal,1,jlev) = controlVector_in(jdim)
            else
              ! both real and imaginary components for jm>0
              jdim = jdim + 1
              hiControlVector_out(ila_mpilocal,1,jlev) = controlVector_in(jdim)
              jdim = jdim + 1
              hiControlVector_out(ila_mpilocal,2,jlev) = controlVector_in(jdim)
            endif
          endif
        enddo
      enddo
    enddo

  end SUBROUTINE BHI_cain


  SUBROUTINE BHI_cainAd(hiControlVector_in,controlVector_out)
    IMPLICIT NONE

    real(8) :: controlVector_out(cvDim_mpilocal)
    real(8) :: hiControlVector_in(nla_mpilocal,2,nkgdimSqrt)

    integer :: jdim, jlev, jm, jn, ila_mpilocal, ila_mpiglobal

    jdim = 0
    do jlev = 1, nkgdimSqrt
      do jm = mymBeg, mymEnd, mymSkip
        do jn = mynBeg, mynEnd, mynSkip
          if(jm.le.jn) then
            ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
            ila_mpilocal  = ilaList_mpilocal(ila_mpiglobal)
            if(jm == 0) then
              ! only real component for jm=0
              jdim = jdim + 1
              controlVector_out(jdim) = controlVector_out(jdim) + hiControlVector_in(ila_mpilocal,1,jlev)
            else
              ! both real and imaginary components for jm>0
              jdim = jdim + 1
              controlVector_out(jdim) = controlVector_out(jdim) + hiControlVector_in(ila_mpilocal,1,jlev)*2.0d0
              jdim = jdim + 1
              controlVector_out(jdim) = controlVector_out(jdim) + hiControlVector_in(ila_mpilocal,2,jlev)*2.0d0
            endif
          endif
        enddo
      enddo
    enddo

  END SUBROUTINE BHI_cainAd


  SUBROUTINE BHI_SPA2GD(hiControlVector_in,gd_out)
    IMPLICIT NONE

    real(8) :: hiControlVector_in(nla_mpilocal,2,nkgdimSqrt)
    real(8) :: gd_out(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)

    real(8) :: sptb(nla_mpilocal,2,nlev_T_even),sp(nla_mpilocal,2,nkgdim)
    real(8) :: tb0(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nlev_T_even)

    integer :: jn,jm,ila_mpilocal,ila_mpiglobal,icount
    real(8) :: sq2, zp
    real(8) , allocatable :: zsp(:,:,:), zsp2(:,:,:)
    integer :: jlev, jlon, jlat, jla_mpilocal, klatPtoT
    real(8), pointer :: zgdpsi(:,:,:),zgdchi(:,:,:)
    real(8), target  :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)
    real(8) :: dla2, dl1sa2, zcoriolis, zpsb(myLonBeg:myLonEnd,myLatBeg:myLatEnd)

    klatPtoT = 1

    ! maybe not needed:
    sp(:,:,:) = 0.0d0
    sptb(:,:,:) = 0.0d0

    sq2 = sqrt(2.0d0)
    allocate(zsp(nkgdimSqrt,2,mymCount))
    allocate(zsp2(nkgdim2,2,mymCount))
    !$OMP PARALLEL DO PRIVATE(jn,jm,jlev,ila_mpiglobal,ila_mpilocal,zsp2,zsp,icount)
    do jn = mynBeg, mynEnd, mynSkip

      icount = 0
      do jm = mymBeg, mymEnd, mymSkip
        if(jm.le.jn) then
          icount = icount+1
          ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
          ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
          do jlev = 1, nkgdimSqrt
            zsp(jlev,1,icount) = hiControlVector_in(ila_mpilocal,1,jlev)
            zsp(jlev,2,icount) = hiControlVector_in(ila_mpilocal,2,jlev)
          enddo
        endif
      enddo

      if(icount.gt.0) then

        CALL DGEMM('N','N',nkgdim2,2*icount,nkgdimSqrt,1.0d0,corns(1,1,jn),nkgdim2,zsp(1,1,1),nkgdimSqrt,0.0d0,zsp2(1,1,1),nkgdim2)

        icount = 0
        do jm = mymBeg, mymEnd, mymSkip
          if(jm.le.jn) then
            icount = icount+1
            ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
            ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
            do jlev = 1, nkgdim
              sp(ila_mpilocal,1,jlev) = zsp2(jlev,1,icount)
              sp(ila_mpilocal,2,jlev) = zsp2(jlev,2,icount)
            enddo
            do jlev = 1, nlev_T
              sptb(ila_mpilocal,1,jlev) = zsp2(jlev+nkgdim,1,icount)
              sptb(ila_mpilocal,2,jlev) = zsp2(jlev+nkgdim,2,icount)
            enddo
          endif
        enddo

      endif

      ! make adjustments for jm=0
      if(mymBeg == 0) then

        ila_mpiglobal = gst_getNind(0,gstID) + jn
        ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)

        do jlev = 1, nkgdim
          sp(ila_mpilocal,1,jlev) = sp(ila_mpilocal,1,jlev)*sq2
          sp(ila_mpilocal,2,jlev) = 0.0d0
        enddo
        do jlev = 1, nlev_T
          sptb(ila_mpilocal,1,jlev) = sptb(ila_mpilocal,1,jlev)*sq2
          sptb(ila_mpilocal,2,jlev) = 0.0d0
        enddo

      endif

    enddo
    !$OMP END PARALLEL DO
    deallocate(zsp)
    deallocate(zsp2)

    !$OMP PARALLEL DO PRIVATE(JLAT,JLEV,JLON)
    do jlev = 1, nkgdim
      do jlat = myLatBeg, myLatEnd
        do jlon = myLonBeg, myLonEnd
          gd(jlon,jlat,jlev) = 0.0d0
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(JLAT,JLEV,JLON)
    do jlev = 1, nlev_T_even
      do jlat = myLatBeg, myLatEnd
        do jlon = myLonBeg, myLonEnd
          tb0(jlon,jlat,jlev) = 0.0d0
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    call gst_setID(gstID)
    call gst_speree(sp,gd)
    call gst_setID(gstID2)
    call gst_speree(sptb,tb0)

    !$OMP PARALLEL DO PRIVATE(jlat,zcoriolis,jlev,jlon,zp)
    do jlat = myLatBeg, myLatEnd
      zcoriolis = 2.d0*ec_Omega*gst_getRmu(jlat,gstID)
      do jlon = myLonBeg, myLonEnd
        zpsb(jlon,jlat) = 0.0d0
        do jlev = 1, nlevPtoT
         zp = zcoriolis*gd(jlon,jlat,nspositVO+jlev-1)
         zpsb(jlon,jlat) = zpsb(jlon,jlat) + PtoT(nlev_T+1,jlev,klatPtoT)*zp
        enddo
      enddo

      do jlev = 1, nlev_T
        do jlon = myLonBeg, myLonEnd
          tb0(jlon,jlat,jlev) = zcoriolis*tb0(jlon,jlat,jlev)
        enddo
      enddo

      do jlev = 1, nkgdim
        do jlon = myLonBeg, myLonEnd
          if(jlev.ne.nspositTG) then
            gd(jlon,jlat,jlev) = gd(jlon,jlat,jlev)*rgsig(jlat,jlev)
          else
            gd(jlon,jlat,jlev) = gd(jlon,jlat,jlev)*tgstdbg(jlon,jlat)
          endif
        enddo
      enddo

      do jlev = 1, nlev_T
        do jlon = myLonBeg, myLonEnd
          tb0(jlon,jlat,jlev) = tb0(jlon,jlat,jlev)*rgsigtb(jlat,jlev)
          gd(jlon,jlat,nspositTT+jlev-1) = gd(jlon,jlat,nspositTT+jlev-1)+tb0(jlon,jlat,jlev)
        enddo
      enddo
      do jlon = myLonBeg, myLonEnd
        zpsb(jlon,jlat) = zpsb(jlon,jlat)*rgsigpsb(jlat)
        gd(jlon,jlat,nspositPS) = gd(jlon,jlat,nspositPS)+zpsb(jlon,jlat)
      enddo
    enddo  ! jlat
    !$OMP END PARALLEL DO

    zgdpsi(myLonBeg:,myLatBeg:,1:) => gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nspositVO:(nspositVO+nlev_M-1))
    zgdchi(myLonBeg:,myLatBeg:,1:) => gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nspositDI:(nspositDI+nlev_M-1))
    !$OMP PARALLEL DO PRIVATE(jlat,jlev,jlon)
    do jlev = nlev_bdl, nlev_M
      do jlat = myLatBeg, myLatEnd
        do jlon = myLonBeg, myLonEnd
          zgdchi(jlon,jlat,jlev) = zgdchi(jlon,jlat,jlev) - tantheta(jlev,jlat)*zgdpsi(jlon,jlat,jlev)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    sp(:,:,:) = 0.0d0

    call gst_setID(gstID)
    call gst_reespe(sp,gd)

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
      enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(JLAT,JLEV,JLON)
    do jlev = 1, nkgdim
      do jlat = myLatBeg, myLatEnd
        do jlon = myLonBeg, myLonEnd
          gd(jlon,jlat,jlev) = 0.0d0
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    call gst_setID(gstID)
    call gst_spgd(sp,gd,nlev_M)

    !$OMP PARALLEL DO PRIVATE(JLAT,JLEV,JLON)
    do jlev = 1, nkgdim
      do jlat = myLatBeg, myLatEnd
        do jlon = myLonBeg, myLonEnd
          gd_out(jlon,jlat,jlev) = gd(jlon,jlat,jlev)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

  END SUBROUTINE BHI_SPA2GD


  SUBROUTINE BHI_SPA2GDAD(gd_in,hiControlVector_out)
    implicit none

    real(8) :: hiControlVector_out(nla_mpilocal,2,nkgdimSqrt)
    real(8) :: gd_in(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)

    real(8) :: sptb(nla_mpilocal,2,nlev_T_even)
    real(8) :: sp(nla_mpilocal,2,nkgdim)
    real(8) :: tb0(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nlev_T_even)

    integer :: jn, jm, ila_mpilocal, ila_mpiglobal, icount
    real(8) :: sq2, zp
    real(8) ,allocatable :: zsp(:,:,:), zsp2(:,:,:)

    integer :: jlev, jlon, jlat, jla_mpilocal, klatPtoT
    real(8) :: dl1sa2, dla2, zcoriolis, zpsb(myLonBeg:myLonEnd,myLatBeg:myLatEnd)
    real(8),pointer :: zgdpsi(:,:,:) ,zgdchi(:,:,:)
    real(8), target :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nkgdim)

    klatPtoT = 1

    !$OMP PARALLEL DO PRIVATE(JLAT,JLEV,JLON)
    do jlev = 1, nkgdim
      do jlat = myLatBeg, myLatEnd
        do jlon = myLonBeg, myLonEnd
          gd(jlon,jlat,jlev) = gd_in(jlon,jlat,jlev)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    call gst_setID(gstID)
    call gst_spgda(sp,gd,nlev_M)

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
      enddo
    enddo
    !$OMP END PARALLEL DO

    call gst_setID(gstID)
    call gst_speree(sp,gd)

    zgdpsi(myLonBeg:,myLatBeg:,1:) => gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nspositVO:(nspositVO+nlev_M-1))
    zgdchi(myLonBeg:,myLatBeg:,1:) => gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nspositDI:(nspositDI+nlev_M-1))
    !$OMP PARALLEL DO PRIVATE(jlat,jlev,jlon)
    do jlev = nlev_bdl, nlev_M
      do jlat = myLatBeg, myLatEnd
        do jlon = myLonBeg, myLonEnd
          zgdpsi(jlon,jlat,jlev) = zgdpsi(jlon,jlat,jlev) - tantheta(jlev,jlat)*zgdchi(jlon,jlat,jlev)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(jlat,zcoriolis,jlev,jlon,zp)
    do jlat = myLatBeg, myLatEnd
      zcoriolis = 2.d0*ec_Omega*gst_getRMU(jlat,gstID)
      tb0(:,jlat,:) = 0.0d0
      do jlev = 1, nlev_T
        do jlon = myLonBeg, myLonEnd
          tb0(jlon,jlat,jlev) = gd(jlon,jlat,nspositTT+jlev-1)
          tb0(jlon,jlat,jlev) = tb0(jlon,jlat,jlev)*rgsigtb(jlat,jlev)
        enddo
      enddo
      do jlon = myLonBeg, myLonEnd
        zpsb(jlon,jlat) = gd(jlon,jlat,nspositPS)
        zpsb(jlon,jlat) = zpsb(jlon,jlat)*rgsigpsb(jlat)
      enddo

      do jlev = 1, nkgdim
        do jlon = myLonBeg, myLonEnd
          if(jlev.ne.nspositTG) then
            gd(jlon,jlat,jlev) = gd(jlon,jlat,jlev)*rgsig(jlat,jlev)
          else
            gd(jlon,jlat,nspositTG) = gd(jlon,jlat,nspositTG)*tgstdbg(jlon,jlat)
          endif
        enddo
      enddo

      do jlev = 1, nlev_T
        do jlon = myLonBeg, myLonEnd
          tb0(jlon,jlat,jlev) = zcoriolis*tb0(jlon,jlat,jlev)
        enddo
      enddo

      do jlev = 1, nlevPtoT
        do jlon = myLonBeg, myLonEnd
          zp = PtoT(nlev_T+1,jlev,klatPtoT)*zpsb(jlon,jlat)
          gd(jlon,jlat,nspositVO+jlev-1) = zcoriolis*zp+gd(jlon,jlat,nspositVO+jlev-1)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO 

    call gst_setID(gstID)
    call gst_reespe(sp,gd)
    call gst_setID(gstID2)
    call gst_reespe(sptb,tb0)

    hiControlVector_out(:,:,:) = 0.0d0
    sq2 = sqrt(2.0d0)
    allocate(zsp(nkgdimSqrt,2,mymCount))
    allocate(zsp2(nkgdim2,2,mymCount))
    !$OMP PARALLEL DO PRIVATE(JN,JM,JLEV,ILA_MPILOCAL,ILA_MPIGLOBAL,zsp,zsp2,icount)
    do jn = mynBeg, mynEnd, mynSkip

      icount = 0
      do jm = mymBeg, mymEnd, mymSkip
        if(jm.le.jn) then
          icount = icount+1
          ila_mpiglobal = gst_getNind(jm,gstID) + jn - jm
          ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
          do jlev = 1, nkgdim
            zsp2(jlev,1,icount) = sp(ila_mpilocal,1,jlev)
            zsp2(jlev,2,icount) = sp(ila_mpilocal,2,jlev)
          enddo
          do jlev = 1, nlev_T
            zsp2(jlev+nkgdim,1,icount) = sptb(ila_mpilocal,1,jlev)
            zsp2(jlev+nkgdim,2,icount) = sptb(ila_mpilocal,2,jlev)
          enddo
        endif
      enddo

      if(icount.gt.0) then

        CALL DGEMM('T','N',nkgdimSqrt,2*icount,nkgdim2,1.0d0,corns(1,1,jn),nkgdim2,zsp2(1,1,1),nkgdim2,0.0d0,zsp(1,1,1),nkgdimSqrt)

        icount = 0
        do jm = mymBeg, jn, mymSkip
          icount=icount+1
          ila_mpiglobal = gst_getNIND(jm,gstID) + jn - jm
          ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)
          do jlev = 1, nkgdimSqrt
            hiControlVector_out(ila_mpilocal,1,jlev) = zsp(jlev,1,icount)
            hiControlVector_out(ila_mpilocal,2,jlev) = zsp(jlev,2,icount)
          enddo
        enddo

      endif

      ! make adjustments for jm=0
      if(mymBeg == 0) then

        ila_mpiglobal = gst_getNIND(0,gstID) + jn
        ila_mpilocal = ilaList_mpilocal(ila_mpiglobal)

        do jlev = 1, nkgdimSqrt
          hiControlVector_out(ila_mpilocal,1,jlev) = hiControlVector_out(ila_mpilocal,1,jlev)*sq2
          hiControlVector_out(ila_mpilocal,2,jlev) = hiControlVector_out(ila_mpilocal,2,jlev)*sq2
        enddo

      endif

    enddo
    !$OMP END PARALLEL DO
    deallocate(zsp)
    deallocate(zsp2)

  END SUBROUTINE BHI_SPA2GDAD


  SUBROUTINE BHI_Finalize()
    implicit none

    if (initialized) then
       deallocate(pressureProfile_M)
       deallocate(pressureProfile_T)
       deallocate(PtoT)
       deallocate(tantheta)
       deallocate(rgsig)
       deallocate(tgstdbg)
       deallocate(rgsigtb)
       deallocate(rgsigpsb)
       deallocate(corns)
       deallocate(rstddev)
    end if

  END SUBROUTINE BHI_Finalize


END MODULE BmatrixHI_mod
