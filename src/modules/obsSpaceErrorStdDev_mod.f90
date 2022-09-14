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

module obsSpaceErrorStdDev_mod
  ! MODULE obsSpaceErrorStdDev_mod (prefix='ose' category='1. High-level functionality')
  !
  ! :Purpose: Contains subroutines for computing background-error and OmP-error
  !           standard deviations in observation space
  !
  use midasMpi_mod
  use obsSpaceData_mod
  use columnData_mod
  use bufr_mod
  use utilities_mod
  use earthConstants_mod
  use MathPhysConstants_mod
  use stateToColumn_mod
  use gridStateVector_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use bmatrixhi_mod
  use obsOperators_mod
  use gps_mod
  use codtyp_mod
  use bCovarSetupChem_mod
  use timeCoord_mod
  use obsTimeInterp_mod
  use bMatrixEnsemble_mod
  use varNameList_mod
  use obsOperatorsChem_mod
  use obsFamilyList_mod

  implicit none
  private

  ! public procedures
  public :: ose_computeStddev

  ! module structures
  ! -----------------

  type :: struct_OmPStdDevCH
     !
     ! Structure containing information retrieved from auxiliary file for holding 
     ! OmP std dev information
     !
     !  Variable               Description
     !  --------               -----------
     !  n_stnid                Number of sub-families (identified via STNIDs)
     !  stnids                 Sub-families (STNIDs; * are wild cards)
     !  element                BUFR element in data block 
     !  source                 0: Set entirely from the auxiliary file being read 
     !                            as a function latitude, vertical level, and month.
     !                         1: Values to be calculated from OmP differences 
     !                            (requires a large data set size for each assim)
     !                            using read in latitudes and vertical levels 
     !                         2: Values to be calculated from OBS_OER and estimated 
     !                            OBS_HPHT (default in absence of stnid and element)
     !  std_type               Index of std dev being read for source=0 or calculated for source=1
     !                         0: absolute value
     !                         1: fraction of measurement value (% / 100)
     !  ibegin                 Position index of start of data for given
     !                         sub-family in the arrays OmPstd,levels,lat,month
     !  n_lvl                  Number of vertical levels (max number when source=2)
     !  levels                 Vertical levels (in coordinate of sub-family data)
     !  n_lat                  Number of latitudes
     !  lat                    Latitudes (degrees; ordered in increasing size)
     !  n_month                Number of months (requires n_lat>1 to be >1)
     !  month                  Months numbered between 1 and 12.  
     !  std                    OmP error std dev (or fraction)

     integer ::  n_stnid
     character(len=12), allocatable :: stnids(:)
     integer, allocatable :: element(:),n_lat(:),n_month(:),std_type(:)
     integer, allocatable :: source(:),ibegin(:),n_lvl(:)
     real(8), allocatable :: std(:)
     real(8), allocatable :: levels(:),lat(:),month(:)

  end type struct_OmPStdDevCH

  type(struct_OmPStdDevCH)  :: OmPstdCH
  type(struct_hco), pointer :: hco_anl => null()

  contains

  !--------------------------------------------------------------------------
  ! ose_computeStddev
  !--------------------------------------------------------------------------
  subroutine ose_computeStddev(columnTrlOnAnlIncLev,hco_anl_in,obsSpaceData)
    !
    !:Purpose: To set OmP-error std dev when possible. Otherwise 
    !          compute background-error stddev in observation space to 
    !          estimate OmP-error std dev.
    !
    implicit none

    ! Arguments:
    type(struct_columnData) :: columnTrlOnAnlIncLev ! Columns of the background interpolated to analysis levels and to obs horizontal locations
    type(struct_hco), pointer :: hco_anl_in
    type(struct_obs) :: obsSpaceData         ! Observation-related data
    
    ! Locals:
    real(8) :: HBHT_static, HBHT_ensemble, HBHT_hybrid

    logical :: staticHBHT = .false.
    logical :: staticOMPE = .false.
    logical :: staticHBHT_ch   = .false.
    logical :: staticOMPE_ch   = .false.
    logical :: ensemble = .false.

    integer :: index_body
    integer :: fnom, fclos, ierr, nulnam

    character(len=12) :: hybrid_mode

    !namelist
    NAMELIST /NAMHBHT/hybrid_mode

    ! set the module variable hco_anl
    hco_anl => hco_anl_in

    ! return if there is no observation.
    if ( obs_numheader(obsSpaceData) == 0 ) return

    !- 1.  Compute HBHT (sigma_b in observation space) or OmP error std dev

    !- 1.1 Compute error std dev from static statistics
    !      OmP error std dev (OBS_OMPE) when possible and required. 
    !      Otherwise, compute HBHT
    !      obsSpaceData - INOUT ( OmP error std dev outputted in OBS_OMPE or/and
    !                            sqrt(diag(H*B*H^T)) with B_static_chm outputted in OBS_HPHT )
  
    call ose_setStaticErrorStddev( columnTrlOnAnlIncLev, obsSpaceData, staticHBHT, staticHBHT_ch, staticOMPE_ch )

    !- 1.2 HBHT from the Bens
    !      obsSpaceData - INOUT (HBensHT std. dev. outputted in OBS_WORK)
    
    call ose_compute_hbht_ensemble( columnTrlOnAnlIncLev,      & 
                                    obsSpaceData, &
                                    ensemble )                               
    
    if ( .not. staticHBHT .and. .not. staticOMPE .and. .not. ensemble .and. &
         .not. staticHBHT_ch .and. .not. staticOMPE_ch ) &
         call utl_abort('ose_computeStddev: no OmP std dev or sqrt(HBHT) was initialized')
    
    !- 2. Select/Blend HBHT.
   
    if ( staticHBHT .and. .not. ensemble ) then
      ! Bnmc only
      write(*,*)
      write(*,*) 'ose_computeStddev: Using B_static ONLY'
      ! HBnmcHT std. dev. already in OBS_HPHT
    else if ( .not. staticHBHT .and. ensemble ) then
      write(*,*)
      write(*,*) 'ose_computeStddev: Using B_ensemble ONLY'
      ! Transfer HBensHT std. dev. values in OBS_WORK to OBS_HPHT
      do index_body = 1, obs_numBody(obsSpaceData)
        HBHT_ensemble = obs_bodyElem_r(obsSpaceData,OBS_WORK,index_body)
        call obs_bodySet_r(obsSpaceData,OBS_HPHT,index_body, HBHT_ensemble)
      end do
    else if ( staticHBHT .and. ensemble ) then
      ! Read Namelist first
      hybrid_mode = 'WEIGHTED_SUM' ! default value
      nulnam = 0
      ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam,nml=namhbht,iostat=ierr)
      if ( ierr /= 0) call utl_abort('ose_computeStddev: Error reading namelist')
      if ( mmpi_myid == 0 ) write(*,nml=namhbht)
      ierr = fclos(nulnam)

      write(*,*)
      write(*,*) 'ose_computeStddev: Using hybrid approach (blend of B_static and B_ensemble) in mode = ', trim(hybrid_mode)

      do index_body = 1, obs_numBody(obsSpaceData)
        HBHT_static = obs_bodyElem_r(obsSpaceData,OBS_HPHT,index_body)
        HBHT_ensemble = obs_bodyElem_r(obsSpaceData,OBS_WORK,index_body)
        if ( HBHT_static <= 0.0d0 ) then
          call obs_bodySet_r(obsSpaceData,OBS_HPHT,index_body, HBHT_ensemble)
        else      
          select case ( trim(hybrid_mode) )
          case ('WEIGHTED_SUM')
            HBHT_hybrid = sqrt(HBHT_static**2 + HBHT_ensemble**2)
          case ('MAX_VALUE')
            HBHT_hybrid = max(HBHT_static,HBHT_ensemble)
          case default
            write(*,*)
            write(*,*) 'ose_computeStddev: Unknown hybrid_mode ', trim(hybrid_mode)
            call utl_abort('ose_compute_HBHT')
          end select
          call obs_bodySet_r(obsSpaceData,OBS_HPHT,index_body, HBHT_hybrid)
        end if
      end do

    end if

  end subroutine ose_computeStddev
  
  !--------------------------------------------------------------------------
  ! ose_setStaticErrorStddev
  !--------------------------------------------------------------------------
  subroutine ose_setStaticErrorStddev( columnTrlOnAnlIncLev, obsSpaceData, statusHBHT, statusHBHT_ch, statusOMPE_ch )
    !
    !:Purpose: To assign or compute the OmP error standard deviations in
    !          observation space where requested. If not possible or available,
    !          compute background-error standard deviation.
    !
    !:Note: OmP error std dev assigment currently only available for the CH obs family. 
    ! 
    !-------------------------------------------------------------------------
    implicit none
  
    ! Arguments:
    type(struct_columnData) :: columnTrlOnAnlIncLev
    type(struct_obs)        :: obsSpaceData  ! observation-space data, output saved in OBS_HPHT column
    logical, intent(inout)  :: statusHBHT, statusHBHT_ch, statusOMPE_ch
    
    ! Locals:
    integer :: famIndex
    character(len=4), allocatable :: availableOMPE(:)
    
    allocate(availableOMPE(ofl_numFamily)) 
    availableOMPE(:) = '    '
           
    ! Assignment of OBS_OMPE in obsSpaceData according to obs family where possible and required
    
    do famIndex=1,ofl_numFamily  
      if ( .not. obs_famExist(obsSpaceData,ofl_familyList(famIndex),localMPI_opt=.true.) ) cycle
      
      if ( ofl_familyList(famIndex) == 'CH' ) then
        write(*,*)
        write(*,*) 'ose_setStaticErrorStddev: Setting of OmP error std dev begins for family ',ofl_familyList(famIndex)
        availableOMPE(famIndex) = ose_setOmPstddevCH(obsSpaceData)
      else
        ! Not available 
        availableOMPE(famIndex) = 'None'
      end if
      if ( trim(availableOMPE(famIndex)) == 'Some' .or. &
           trim(availableOMPE(famIndex))  == 'All' ) statusOMPE_ch = .true.

      if ( trim(availableOMPE(famIndex)) == '' ) then
        write(*,*) 'ose_setStaticErrorStddev: No ',ofl_familyList(famIndex),' obs. Setting of error std dev not required.'
      else if ( trim(availableOMPE(famIndex)) == 'None' ) then
        if ( ofl_familyList(famIndex) == 'CH' ) then
          write(*,*) 'ose_setStaticErrorStddev: Setting of ',ofl_familyList(famIndex),' OmP error std dev to be estimated via HBHT calc for all obs'
        end if        
      else if ( trim(availableOMPE(famIndex)) == 'Some' ) then
        write(*,*) 'ose_setStaticErrorStddev: Setting of ',ofl_familyList(famIndex),' OmP error std dev partially completed (some HBHT calc needed to complete)'       
      else
        write(*,*) 'ose_setStaticErrorStddev: Setting of ',ofl_familyList(famIndex),' OmP error std dev completed (no HBHT calc needed)'
      end if
      
    end do
    
    ! Computation of background-error standard deviation for obs families or groups of families when required

    ! HBHT from the Bnmc matrix for weather 
    ! obsSpaceData - INOUT (HBnmcHT std. dev. output in OBS_HPHT)
     
    if ( any(ofl_familyList /= 'CH' .and. ofl_familyList /= 'TO' .and. &
       ( availableOMPE == 'Some' .or. availableOMPE == 'None' ) ) ) &
       call ose_compute_hbht_static( columnTrlOnAnlIncLev, obsSpaceData, statusHBHT )

    ! HBHT from the B matrix for constituents
    
    if ( any(ofl_familyList == 'CH' .and. ( availableOMPE == 'Some' .or. availableOMPE == 'None' ) ) ) &
       call ose_compute_hbht_static_chem( columnTrlOnAnlIncLev, obsSpaceData, statusHBHT_ch )
    
    deallocate(availableOMPE)
    
  end subroutine ose_setStaticErrorStddev

  !--------------------------------------------------------------------------
  ! ose_compute_hbht_static
  !--------------------------------------------------------------------------
  subroutine ose_compute_hbht_static(columnTrlOnAnlIncLev,lobsSpaceData,active)
    !
    !:Purpose: To compute background-error stddev in observation space using
    !          fixed statistics specific in stats file.
    !
    implicit none

      ! Arguments:
      type(struct_obs)        :: lobsSpaceData
      type(struct_columnData) :: columnTrlOnAnlIncLev
      logical                 :: active
      
      ! Locals:
      type(struct_vco), pointer        :: vco_anl
      type(struct_columnData) :: column
      type(struct_gsv)        :: statevector

      INTEGER JLAT, JLON, JLEV, JOBS

      CHARACTER*12 CLETIKET
      CHARACTER*2 CLTYPVAR
      CHARACTER*1 CLGRTYP
      CHARACTER*4 CLNOMVAR

      INTEGER IULSSF,IDATEO
      INTEGER FSTPRM,FNOM,FSTOUV
      INTEGER IKEY,IERR,IDATE

      REAL*8, allocatable :: ZBUFFER(:,:)
      real*8, pointer     :: height_column(:), tt_column(:), field_ptr(:,:,:)

      INTEGER INI,INJ,INK, INPAS, INBITS, IDATYP, IDEET
      INTEGER IP1,IP2,IP3,IG1,IG2,IG3,IG4,ISWA,ILENGTH,IDLTF
      INTEGER IUBC,IEXTR1,IEXTR2,IEXTR3

      integer :: nLev_M,nLev_T,status,shift_level,Vcode_anl

      integer :: cvdim

      real(8), allocatable  :: scaleFactor(:)

      !- Get the appropriate Vertical Coordinate
      vco_anl => col_getVco(columnTrlOnAnlIncLev)

      ! Note : Here we can use the global B_hi even if we are in LAM mode since, 
      !        in BackgroundCheck mode, the only purpose of bhi_setup is to read 
      !        the scaleFactor and to check the consistency between the vco from 
      !        analysisgrid and cov file
      call bhi_setup( hco_anl,vco_anl,           &  ! IN
                      cvdim,                     &  ! OUT
                      mode_opt='BackgroundCheck' )  ! IN

      if ( cvdim > 0 ) then
         write(*,*)
         write(*,*) 'Computing HBHT from Bnmc - Start'
         active = .true.
      else
         if ( mmpi_myid == 0 ) write(*,*) 'ose_compute_hbht_static: option NOT ACTIVATED'
         active = .false.
         return
      end if

      nlev_T = col_getNumLev(columnTrlOnAnlIncLev,'TH')
      nlev_M = col_getNumLev(columnTrlOnAnlIncLev,'MM')

      allocate(scaleFactor(max(nLev_M,nLev_T)))
      call bhi_getScaleFactor(scaleFactor)

      status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)
      if ( Vcode_anl == 5001 ) then
        shift_level = 0
      else if ( Vcode_anl == 5002) then
        shift_level = 1
      else if ( Vcode_anl == 5005 ) then
        shift_level = 0
      else
        write(*,*) 'Vcode_anl = ',Vcode_anl
        call utl_abort('ose_compute_hbht_static: unknown vertical coordinate type!')
      end if

      allocate(ZBUFFER(HCO_ANL%NI,HCO_ANL%NJ))

      call gsv_allocate(statevector, 1, hco_anl, vco_anl, &
                        allocHeight_opt=.false., allocPressure_opt=.false.)
      call gsv_zero(statevector)

      call col_setVco(column,col_getVco(columnTrlOnAnlIncLev))
      call col_allocate(column,col_getNumCol(columnTrlOnAnlIncLev))

      ! Set the value of OBS_LYR required by setfge routines

      call oop_vobslyrs(columnTrlOnAnlIncLev,lobsSpaceData,beSilent=.false.)

      ! 1. Opening the statistics file

      IULSSF=0
      IERR=FNOM(iulssf,'./bgcov','RND+OLD+R/O',0)
      IF ( IERR .EQ. 0 ) THEN
        write(*,*) 'IBGST - File : ./bgcov'
        write(*,*) ' opened as unit file ',iulssf
        ierr =  fstouv(iulssf,'RND+OLD')
      ELSE
        CALL utl_abort('HBHT_static:NO BACKGROUND STAT FILE!!')
      ENDIF

      ! 2.1 Background error standard deviations

      CLETIKET = 'BGCK_STDDEV'
      write(*,*) 'HBHT_static: CLETIKET = ',CLETIKET
      IDATE    = -1
      IP2      = -1
      IP3      = -1
      CLTYPVAR =' '

      ! READ IN STANDARD DEVIATION FOR EACH OBSERVATION TYPE

      clnomvar = 'UU'
      write(*,*) clnomvar 
      call gsv_getField(statevector,field_ptr,'UU')       
      do jlev = 1, nlev_M
         ip1 = vco_anl%ip1_M(jlev)
         ikey = utl_fstlir(zbuffer,iulssf,ini,inj,ink,idate,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
         ierr = fstprm(ikey,idateo,ideet,inpas   &
              ,ini,inj,ink, inbits, idatyp       &
              ,ip1,ip2,ip3,cltypvar,clnomvar,cletiket,clgrtyp   &
              ,ig1,ig2,ig3,ig4,iswa,ilength,idltf   &
              ,iubc,iextr1,iextr2,iextr3)
         if (ini /= hco_anl%ni .or. inj /= hco_anl%nj .or. ink /= 1) then
            write(*,*)
            write(*,*) 'HBHT_static: Invalid dimensions for...'
            write(*,*) 'nomvar         =', trim(clnomvar)
            write(*,*) 'etiket         =', trim(cletiket)
            write(*,*) 'Found ni,nj,nk =', ini, inj, ink
            write(*,*) 'Should be      =', hco_anl%ni, hco_anl%nj, nlev_M
            call utl_abort('ose_compute_hbht_static')
         end if
         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scaleFactor(jlev+shift_level)*zbuffer(jlon,jlat)*MPC_M_PER_S_PER_KNOT_R8
             end do
          end do
      end do

      clnomvar = 'VV'
      write(*,*) clnomvar
      call gsv_getField(statevector,field_ptr,'VV')
      do jlev = 1, nlev_M
         ip1 = vco_anl%ip1_M(jlev)
         ikey = utl_fstlir(zbuffer,iulssf,ini,inj,ink,idate,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
         if (ini /= hco_anl%ni .or. inj /= hco_anl%nj .or. ink /= 1) then
            write(*,*)
            write(*,*) 'HBHT_static: Invalid dimensions for...'
            write(*,*) 'nomvar         =', trim(clnomvar)
            write(*,*) 'etiket         =', trim(cletiket)
            write(*,*) 'Found ni,nj,nk =', ini, inj, ink
            write(*,*) 'Should be      =', hco_anl%ni, hco_anl%nj, nlev_M
            call utl_abort('ose_compute_hbht_static')
         end if

         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scalefactor(jlev+shift_level)*zbuffer(jlon,jlat)*MPC_M_PER_S_PER_KNOT_R8
            end do
         end do
      end do

      clnomvar = 'ES'
      write(*,*) clnomvar
      call gsv_getField(statevector,field_ptr,'HU')
      do jlev = 1, nlev_T
         ip1 = vco_anl%ip1_T(jlev)
         ikey = utl_fstlir(zbuffer,iulssf,ini,inj,ink,idate,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
         if (ini /= hco_anl%ni .or. inj /= hco_anl%nj .or. ink /= 1) then
            write(*,*)
            write(*,*) 'HBHT_static: Invalid dimensions for...'
            write(*,*) 'nomvar         =', trim(clnomvar)
            write(*,*) 'etiket         =', trim(cletiket)
            write(*,*) 'Found ni,nj,nk =', ini, inj, ink
            write(*,*) 'Should be      =', hco_anl%ni, hco_anl%nj, nlev_T
            call utl_abort('ose_compute_hbht_static')
         end if
         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scaleFactor(jlev)*zbuffer(jlon,jlat)
            end do
         end do
      end do

      ! height is put into TT slot in gridStateVector
      clnomvar = 'GZ'
      write(*,*) clnomvar
      call gsv_getField(statevector,field_ptr,'TT')
      do jlev = 1, nlev_T
         ip1 = vco_anl%ip1_T(jlev)
         ikey = utl_fstlir(zbuffer,iulssf,ini,inj,ink,idate,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
         if (ini /= hco_anl%ni .or. inj /= hco_anl%nj .or. ink /= 1) then
            write(*,*)
            write(*,*) 'HBHT_static: Invalid dimensions for...'
            write(*,*) 'nomvar         =', trim(clnomvar)
            write(*,*) 'etiket         =', trim(cletiket)
            write(*,*) 'Found ni,nj,nk =', ini, inj, ink
            write(*,*) 'Should be      =', hco_anl%ni, hco_anl%nj, nlev_T
            call utl_abort('ose_compute_hbht_static')
         end if
         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scaleFactor(jlev)*zbuffer(jlon,jlat)*ec_rg*10.d0
            end do
         end do
      end do

      call s2c_bgcheck_bilin(column,statevector,lobsSpaceData)

      ! copy height data from TT to height slot in columnData
      do jobs= 1, col_getNumCol(column)
         height_column => col_getColumn(column,jobs,'Z_T')
         tt_column => col_getColumn(column,jobs,'TT')
         do jlev = 1,col_getNumLev(column,'TH')
            height_column(jlev)=tt_column(jlev)
         enddo
      enddo

      ! SET THE FIRST-GUESS ERRORS FOR CONVENTIONAL DATA ON PRESSURE LEVELS
      ! --------------------------------------------------------------------

      call setfgefam('AI',column,columnTrlOnAnlIncLev,lobsSpaceData)
      call setfgefam('SW',column,columnTrlOnAnlIncLev,lobsSpaceData)
      call setfgefam('UA',column,columnTrlOnAnlIncLev,lobsSpaceData)
      call setfgefam('SF',column,columnTrlOnAnlIncLev,lobsSpaceData)
      call setfgefam('HU',column,columnTrlOnAnlIncLev,lobsSpaceData)
      call setfgefamz('PR',column,columnTrlOnAnlIncLev,lobsSpaceData)
      call setfgefamz('AL',column,columnTrlOnAnlIncLev,lobsSpaceData)
      call setfgefamz('RA',column,columnTrlOnAnlIncLev,lobsSpaceData)

      ! SET THE FIRST-GUESS ERRORS FOR RADIO OCCULTATION DATA
      ! -----------------------------------------------------

      call setfgedif('RO',columnTrlOnAnlIncLev,lobsSpaceData)

      ! DO TEMPERATURE FIRST-GUESS ERROR
      ! ---------------------------------

      clnomvar = 'TT'
      write(*,*) clnomvar
      call gsv_getField(statevector,field_ptr,'TT')
      do jlev = 1, nlev_T
         ip1 = vco_anl%ip1_T(jlev)
         ikey = utl_fstlir(zbuffer,iulssf,ini,inj,ink,idate,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
         if (ini /= hco_anl%ni .or. inj /= hco_anl%nj .or. ink /= 1) then
            write(*,*)
            write(*,*) 'HBHT_static: Invalid dimensions for...'
            write(*,*) 'nomvar         =', trim(clnomvar)
            write(*,*) 'etiket         =', trim(cletiket)
            write(*,*) 'Found ni,nj,nk =', ini, inj, ink
            write(*,*) 'Should be      =', hco_anl%ni, hco_anl%nj, nlev_T
            call utl_abort('ose_compute_hbht_static')
         end if
         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scaleFactor(jlev)*zbuffer(jlon,jlat)
            end do
         end do
      end do

      call s2c_bgcheck_bilin(column,statevector,lobsSpaceData)
      call setfgett(column,columnTrlOnAnlIncLev,lobsSpaceData)

      ! RELOAD DATA TO DO SURFACE FIRST-GUESS ERRORS
      ! ---------------------------------------------

      clnomvar = 'P0'
      write(*,*) clnomvar
      ip1 = -1
      ikey = utl_fstlir(zbuffer,iulssf,ini,inj,ink,idate,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
      if (ini /= hco_anl%ni .or. inj /= hco_anl%nj .or. ink /= 1) then
         write(*,*)
         write(*,*) 'HBHT_static: Invalid dimensions for...'
         write(*,*) 'nomvar         =', trim(CLNOMVAR)
         write(*,*) 'etiket         =', trim(CLETIKET)
         write(*,*) 'Found ni,nj,nk =', ini, inj, ink
         write(*,*) 'Should be      =', hco_anl%ni, hco_anl%nj, 1
         call utl_abort('ose_compute_hbht_static')
      end if
      call gsv_getField(statevector,field_ptr,'P0')
      do jlev = 1, ink
         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scaleFactor(max(nLev_M,nLev_T))*zbuffer(jlon,jlat)*MPC_PA_PER_MBAR_R8
            end do
         end do
      end do

      clnomvar = 'UU'
      write(*,*) clnomvar
      call gsv_getField(statevector,field_ptr,'UU')
      do jlev = 1, nlev_M
         ip1 = vco_anl%ip1_M(jlev)
         ikey = utl_fstlir(zbuffer,iulssf,ini,inj,ink,idate,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
         if (ini /= hco_anl%ni .or. inj /= hco_anl%nj .or. ink /= 1) then
            write(*,*)
            write(*,*) 'HBHT_static: Invalid dimensions for...'
            write(*,*) 'nomvar         =', trim(clnomvar)
            write(*,*) 'etiket         =', trim(cletiket)
            write(*,*) 'Found ni,nj,nk =', ini, inj, ink
            write(*,*) 'Should be      =', hco_anl%ni, hco_anl%nj, nlev_M
            call utl_abort('ose_compute_hbht_static')
         end if
         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scaleFactor(jlev+shift_level)*zbuffer(jlon,jlat)*MPC_M_PER_S_PER_KNOT_R8
            end do
         end do
      end do

      clnomvar = 'VV'
      write(*,*) clnomvar
      call gsv_getField(statevector,field_ptr,'VV')
      do jlev = 1, nlev_M
         ip1 = vco_anl%ip1_M(jlev)
         ikey = utl_fstlir(zbuffer,iulssf,ini,inj,ink,idate,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
         if (ini /= hco_anl%ni .or. inj /= hco_anl%nj .or. ink /= 1) then
            write(*,*)
            write(*,*) 'HBHT_static: Invalid dimensions for...'
            write(*,*) 'nomvar         =', trim(clnomvar)
            write(*,*) 'etiket         =', trim(cletiket)
            write(*,*) 'Found ni,nj,nk =', ini, inj, ink
            write(*,*) 'Should be      =', hco_anl%ni, hco_anl%nj, nlev_M
            call utl_abort('ose_compute_hbht_static')
         end if
         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scaleFactor(jlev+shift_level)*zbuffer(jlon,jlat)*MPC_M_PER_S_PER_KNOT_R8
            end do
         end do
      end do

      clnomvar = 'TT'
      write(*,*) clnomvar
      call gsv_getField(statevector,field_ptr,'TT')
      do jlev = 1, nlev_T
         ip1 = vco_anl%ip1_T(jlev)
         ikey = utl_fstlir(zbuffer,iulssf,ini,inj,ink,idate,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
         if (ini /= hco_anl%ni .or. inj /= hco_anl%nj .or. ink /= 1) then
            write(*,*)
            write(*,*) 'HBHT_static: Invalid dimensions for...'
            write(*,*) 'nomvar         =', trim(clnomvar)
            write(*,*) 'etiket         =', trim(cletiket)
            write(*,*) 'Found ni,nj,nk =', ini, inj, ink
            write(*,*) 'Should be      =', hco_anl%ni, hco_anl%nj, nlev_T
            call utl_abort('ose_compute_hbht_static')
         end if
         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scaleFactor(jlev)*zbuffer(jlon,jlat)
            end do
         end do
      end do

      clnomvar = 'ES'
      write(*,*) clnomvar
      call gsv_getField(statevector,field_ptr,'HU')
      do jlev = 1, nlev_T
         ip1 = vco_anl%ip1_T(jlev)
         ikey = utl_fstlir(zbuffer,iulssf,ini,inj,ink,idate,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
         if (ini /= hco_anl%ni .or. inj /= hco_anl%nj .or. ink /= 1) then
            write(*,*)
            write(*,*) 'HBHT_static: Invalid dimensions for...'
            write(*,*) 'nomvar         =', trim(clnomvar)
            write(*,*) 'etiket         =', trim(cletiket)
            write(*,*) 'Found ni,nj,nk =', ini, inj, ink
            write(*,*) 'Should be      =', hco_anl%ni, hco_anl%nj, nlev_T
            call utl_abort('ose_compute_hbht_static')
         end if
         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scaleFactor(jlev)*zbuffer(jlon,jlat)
            end do
         end do
      end do

      clnomvar = 'LVIS'
      if (gsv_varExist(statevector,clnomvar)) then
        write(*,*) clnomvar
        call gsv_getField(statevector,field_ptr,clnomvar)
        do jlev = 1, nlev_T
          ip1 = vco_anl%ip1_T(jlev)
          ikey = utl_fstlir(zbuffer,iulssf,ini,inj,ink,idate,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
          if (ini /= hco_anl%ni .or. inj /= hco_anl%nj .or. ink /= 1) then
            write(*,*)
            write(*,*) 'HBHT_static: Invalid dimensions for...'
            write(*,*) 'nomvar         =', trim(clnomvar)
            write(*,*) 'etiket         =', trim(cletiket)
            write(*,*) 'Found ni,nj,nk =', ini, inj, ink
            write(*,*) 'Should be      =', hco_anl%ni, hco_anl%nj, nlev_T
            call utl_abort('ose_compute_hbht_static')
         end if
          do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
              field_ptr(jlon,jlat,jlev) = scaleFactor(jlev)*zbuffer(jlon,jlat)
            end do
          end do
        end do
      end if

      clnomvar = 'WGE'
      if (gsv_varExist(statevector,clnomvar)) then
        write(*,*) clnomvar
        ip1 = -1
        ikey = utl_fstlir(zbuffer,iulssf,ini,inj,ink,idate,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
        if (ini /= hco_anl%ni .or. inj /= hco_anl%nj .or. ink /= 1) then
          write(*,*)
          write(*,*) 'HBHT_static: Invalid dimensions for...'
          write(*,*) 'nomvar         =', trim(CLNOMVAR)
          write(*,*) 'etiket         =', trim(CLETIKET)
          write(*,*) 'Found ni,nj,nk =', ini, inj, ink
          write(*,*) 'Should be      =', hco_anl%ni, hco_anl%nj, 1
          call utl_abort('ose_compute_hbht_static')
        end if
        call gsv_getField(statevector,field_ptr,clnomvar)
        do jlev = 1, ink
          do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
              field_ptr(jlon,jlat,jlev) = scaleFactor(max(nLev_M,nLev_T))*zbuffer(jlon,jlat)
            end do
          end do
        end do
      end if

      call s2c_bgcheck_bilin(column,statevector,lobsSpaceData)

      ! SET THE FIRST-GUESS ERRORS FOR THE SURFACE DATA
      ! ------------------------------------------------

      call setfgesurf(column,columnTrlOnAnlIncLev,lobsSpaceData)

      ! READ IN LN Q FIRST-GUESS ERRORS FOR SETFGEGPS
      ! ---------------------------------------------
      
      clnomvar = 'LQ'
      write(*,*) clnomvar
      call gsv_getField(statevector,field_ptr,'HU')
      do jlev = 1, nlev_T
         ip1 = vco_anl%ip1_T(jlev)
         ikey = utl_fstlir(zbuffer,iulssf,ini,inj,ink,idate,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
         if (ini /= hco_anl%ni .or. inj /= hco_anl%nj .or. ink /= 1) then
            write(*,*)
            write(*,*) 'HBHT_static: Invalid dimensions for...'
            write(*,*) 'nomvar         =', trim(clnomvar)
            write(*,*) 'etiket         =', trim(cletiket)
            write(*,*) 'Found ni,nj,nk =', ini, inj, ink
            write(*,*) 'Should be      =', hco_anl%ni, hco_anl%nj, nlev_T
            call utl_abort('ose_compute_hbht_static')
         end if
         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scaleFactor(jlev)*zbuffer(jlon,jlat)
            end do
         end do
      end do

      call s2c_bgcheck_bilin(column,statevector,lobsSpaceData)

      ! OPTIONAL TEST OF THE GB-GPS ZTD OPERATOR JACOBIAN
      ! -------------------------------------------------

      if (ltestop) call setfgegps(column,columnTrlOnAnlIncLev,lobsSpaceData)

      deallocate(scaleFactor)
      call col_deallocate(column)
      deallocate(zbuffer)

      write(*,*)
      write(*,*) 'Computing HBHT from Bnmc - END'
      
  end subroutine ose_compute_hbht_static
  
  !--------------------------------------------------------------------------
  ! ose_compute_hbht_static_chem
  !--------------------------------------------------------------------------
  subroutine ose_compute_hbht_static_chem(columnTrlOnAnlIncLev,obsSpaceData,active)
    !
    !:Purpose: To compute the background error standard deviations in
    !          observation space, sqrt(diag(H*B_static*H^T)).
    !
    implicit none
  
    ! Arguments:
    type(struct_columnData) :: columnTrlOnAnlIncLev      ! column at observation location
    type(struct_obs)        :: obsSpaceData ! observation-space data, output saved in OBS_HPHT column
    logical                 :: active        ! flag to indicate if chemical constituents are to be used

    ! Locals:
    type(struct_vco), pointer :: vco_anl
 
    !- Get the appropriate Vertical Coordinate
    vco_anl => col_getVco(columnTrlOnAnlIncLev)
  
    call bcsc_setupCH( hco_anl,vco_anl,active,'BackgroundCheck' )
  
    if (active) then
      write(*,*)
      write(*,*) 'Computing H*B*H^T using B_static_chm - Start'
    else
      if ( mmpi_myid == 0 ) write(*,*) 'ose_compute_HBHT_static_chem: option NOT ACTIVATED'
      return
    end if
          
    call oopc_CHobsoperators(columnTrlOnAnlIncLev,obsSpaceData,'HBHT')
  
    write(*,*)
    write(*,*) 'Computing H*B*H^T using B_static_chm - End'
  
    RETURN
  end subroutine ose_compute_hbht_static_chem

  !--------------------------------------------------------------------------
  ! ose_compute_hbht_ensemble
  !--------------------------------------------------------------------------
  subroutine ose_compute_hbht_ensemble(columnTrlOnAnlIncLev,obsSpaceData,active)
    !
    !:Purpose: To compute background-error stddev in observation space using
    !          ensemble-based statistics.
    !
    implicit none

    ! Arguments:
    type(struct_columnData) :: columnTrlOnAnlIncLev      ! Columns of the background interpolated to analysis levels and to obs horizontal locations
    type(struct_obs)        :: obsSpaceData ! Observation-related data
    logical                 :: active

    ! Locals:
    type(struct_columnData) :: column
    type(struct_gsv)        :: statevector

    type(struct_vco), pointer :: vco_anl

    real(8), allocatable :: HBHT_ens(:)

    integer :: memberIndex, index_body
    integer, allocatable :: cvdim(:)

    !
    !- 1.  Initialization
    !

    !- 1.1 Get vertical analysis grid attributes
    vco_anl => col_getVco(columnTrlOnAnlIncLev)

    !- 1.2 Initialize/Read the flow-dependent ensemble perturbations
    call ben_Setup( hco_anl, hco_anl,    & ! IN
                    vco_anl,             & ! IN
                    cvdim,               & ! OUT
                    'BackgroundCheck' )    ! IN

    if ( cvdim(1) > 0 ) then
      write(*,*)
      write(*,*) 'Computing HBHT from ensemble perturbations - START'
      active = .true.
    else
      if ( mmpi_myid == 0 ) write(*,*) 'ose_compute_hbht_ensemble: option NOT ACTIVATED'
      active = .false.
      return
    end if

    !- 1.3 Create a gridstatevector to store the ensemble perturbations
    call gsv_allocate(statevector, tim_nstepobsinc, hco_anl, vco_anl, &
      mpi_local_opt=.true., allocHeight_opt=.false., allocPressure_opt=.false.)

    !- 1.4 Create column vectors to store the ens perturbation interpolated to obs horizontal locations
    call col_setVco(column,vco_anl)
    call col_allocate(column,col_getNumCol(columnTrlOnAnlIncLev),mpiLocal_opt=.true.)

    !- 1.5 Create a working a array to sum H ensPert HT
    allocate(HBHT_ens(obs_numBody(obsSpaceData)))
    HBHT_ens(:) = 0.d0

    !- 1.6
    call oti_timeBinning(obsSpaceData,tim_nstepobsinc)

    !
    !- 2.  Compute HBHT from the ensemble perturbations
    !
    do memberIndex = 1, ben_getnEns()

      !- 2.1 Extract perturbations from the current memberIndex
      write(*,*)
      write(*,*) 'Reading ensemble perturbation from member = ', memberIndex
      call ben_getPerturbation( statevector,    &
                                memberIndex,    &
                                'ConstantValue' )

      !- 2.2 Interpolation to the observation horizontal locations
      !       column - OUT (H_horiz EnsPert)
      call s2c_tl( statevector,           &
                   column,                &
                   columnTrlOnAnlIncLev, obsSpaceData )
                   
      !- 2.3 Interpolation to observation space
      !         obsSpaceData - OUT (Save as OBS_WORK: H_vert H_horiz EnsPert = H EnsPert)
      call oop_Htl( column, columnTrlOnAnlIncLev, &
                    obsSpaceData,    &
                    1 )

      !- 2.4 alpha * HBH^T = sum(OBS_WORK^2)
      do index_body = 1, obs_numBody(obsSpaceData)
        HBHT_ens(index_body) = HBHT_ens(index_body) + &
                               (obs_bodyElem_r(obsSpaceData,OBS_WORK,index_body))**2
      end do

    end do

    !- 2.5 Insert the standard deviations in OBS_WORK
    do index_body = 1, obs_numBody(obsSpaceData)
      call obs_bodySet_r(obsSpaceData,OBS_WORK,index_body,sqrt(HBHT_ens(index_body)))
    end do

    !
    !- 3.  Ending/Deallocation
    !
    deallocate(HBHT_ens)
    call col_deallocate(column)
    call gsv_deallocate(statevector)

    write(*,*)
    write(*,*) 'Computing HBHT from ensemble perturbations - END'

  end subroutine ose_compute_hbht_ensemble

  !-------------------- Weather obs FGE std dev routines --------------------
  
  !--------------------------------------------------------------------------
  ! setfgefam
  !--------------------------------------------------------------------------
  subroutine setfgefam(cdfam,column,columnTrlOnAnlIncLev,lobsSpaceData)
    !
    !:Purpose: To interpolate vertically the contents of "column" to
    !          the pressure levels of the observations. Then to compute
    !          THE FIRST GUESS ERROR VARIANCES. A linear interpolation in ln(p)
    !          is performed.
    !
    implicit none

    ! Arguments:
    character*2 cdfam
    type(struct_columnData) :: column
    type(struct_columnData) :: columnTrlOnAnlIncLev
    type(struct_obs) :: lobsSpaceData

    ! Locals:
      INTEGER IPB,IPT
      INTEGER INDEX_HEADER,ITYP,IK
      INTEGER INDEX_BODY
      REAL*8 ZWB,ZWT, columnElem
      REAL*8 ZLEV,ZPB,ZPT
      character(len=4) :: varLevel

      ! loop over all header indices of the CDFAM family
      call obs_set_current_header_list(lobsSpaceData,CDFAM)
      HEADER: do
         index_header = obs_getHeaderIndex(lobsSpaceData)
         if (index_header < 0) exit HEADER

         ! loop over all body indices for this index_header
         call obs_set_current_body_list(lobsSpaceData, index_header)
         BODY: do 
            index_body = obs_getBodyIndex(lobsSpaceData)
            if (index_body < 0) exit BODY
            !
            !*    1. Computation of sigmap
            !     .  -----------------------------
            !
            IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,index_body) == obs_assimilated .AND.   &
                 obs_bodyElem_i(lobsSpaceData,OBS_VCO,index_body) .EQ. 2      ) then
            IF  (obs_bodyElem_i(lobsSpaceData,OBS_XTR,index_body) .NE. 0) THEN
               ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body)
               varLevel = vnl_varLevelFromVarnum(ityp)
               IK   = col_getNumLev(columnTrlOnAnlIncLev,varLevel)
               IPB  = IK + col_getOffsetFromVarno(columnTrlOnAnlIncLev,ityp)
               if(ITYP .ne. BUFR_NEGZ) then
                 call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,col_getElem(column,IPB,INDEX_HEADER))
               else
                 call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,col_getHeight(column,IK,INDEX_HEADER,'TH'))
               endif
            ELSE
               ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body)
               varLevel = vnl_varLevelFromVarnum(ityp)
               ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,index_body)
               IK   = obs_bodyElem_i(lobsSpaceData,OBS_LYR,index_body)
               IPT  = IK + col_getOffsetFromVarno(columnTrlOnAnlIncLev,ityp)
               IPB  = IPT+1
               ZPT  = col_getPressure(columnTrlOnAnlIncLev,IK,INDEX_HEADER,varLevel)
               ZPB  = col_getPressure(columnTrlOnAnlIncLev,IK+1,INDEX_HEADER,varLevel)
               ZWB  = LOG(ZLEV/ZPT)/LOG(ZPB/ZPT)
               ZWT  = 1.0D0 - ZWB

               ! FIRST GUESS ERROR VARIANCE

               if(ITYP .ne. BUFR_NEGZ) then
                 call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,   &
                      (ZWB*col_getElem(column,IPB,INDEX_HEADER) + ZWT*col_getElem(column,IPT,INDEX_HEADER)))
               else
                 call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,   &
                      (ZWB*col_getHeight(column,IK+1,INDEX_HEADER,'TH') + ZWT*col_getHeight(column,IK,INDEX_HEADER,'TH')))
               endif
               if(obs_bodyElem_r(lobsSpaceData,OBS_HPHT,index_body).le.0.d0) then
                 write(*,*) 'SETFGEFAM: CDFAM = ',CDFAM
                 write(*,*) 'SETFGEFAM: IPB,IPT,ZWB,ZWT,ITYP,ZLEV=',IPB,IPT,ZWB,ZWT,ITYP,ZLEV
                 columnElem = col_getElem(column,IPB,INDEX_HEADER)
                 write(*,*) 'SETFGEFAM: column_all(IPB,INDEX_HEADER)=',columnElem
                 columnElem = col_getElem(column,IPT,INDEX_HEADER)
                 write(*,*) 'SETFGEFAM: column_all(IPT,INDEX_HEADER)=',columnElem
                 columnElem = col_getHeight(column,IK+1,INDEX_HEADER,'TH')
                 write(*,*) 'SETFGEFAM: get_height(IK+1,INDEX_HEADER)=',columnElem
                 columnElem = col_getHeight(column,IK  ,INDEX_HEADER,'TH')
                 write(*,*) 'SETFGEFAM: get_height(IK  ,INDEX_HEADER)=',columnElem
                 CALL utl_abort('SETFGEFAM: First-guess stdev bad value')
               endif
            ENDIF
            ENDIF

         END DO BODY

      END DO HEADER

      RETURN
  end subroutine setfgefam

  !--------------------------------------------------------------------------
  ! setfgefamz
  !--------------------------------------------------------------------------
  subroutine setfgefamz(cdfam, column, columnTrlOnAnlIncLev, obsSpaceData)
    !
    !:Purpose: To interpolate vertically the contents of "column" to the levels
    !          of the observations (in meters). Then to compute THE FIRST GUESS
    !          ERROR VARIANCES. A linear interpolation in z is performed.
    !
    implicit none

    ! Arguments:
    character(len=2),        intent(in)    :: cdfam
    type(struct_columnData), intent(in)    :: column
    type(struct_columnData), intent(in)    :: columnTrlOnAnlIncLev
    type(struct_obs),        intent(inout) :: obsSpaceData

    ! Locals:
    integer :: IPB, IPT
    integer :: headerIndex, ITYP, IK
    integer :: bodyIndex
    integer :: bodyIndexStart, bodyIndexEnd, bodyIndex2
    real*8  :: ZWB, ZWT, sigmap_uuT, sigmap_vvT , sigmap_uuB, sigmap_vvB 
    real*8  :: ZLEV, ZPB, ZPT
    real(8) :: azimuth, fge_uu, fge_vv, fge_fam
    character(len=4) :: varLevel

    ! loop over all header indices of the CDFAM family
    call obs_set_current_header_list(obsSpaceData, CDFAM)
    HEADER: do
      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if ( headerIndex < 0 ) exit HEADER

      if  ( cdfam == 'RA' ) then
        ! Azimuth of the radar beam
         azimuth = obs_headElem_r(obsSpaceData, OBS_RZAM, headerIndex )
      end if

      ! loop over all body indices for this headerIndex
      call obs_set_current_body_list(obsSpaceData, headerIndex)
      BODY: do 
        bodyIndex = obs_getBodyIndex(obsSpaceData)
        if ( bodyIndex < 0 ) exit BODY

        !*    1. Computation of sigmap
        !     .  -----------------------------
        if ( obs_bodyElem_i(obsSpaceData, OBS_ASS, bodyIndex) == obs_assimilated .and. &
             obs_bodyElem_i(obsSpaceData, OBS_VCO, bodyIndex) == 1 )then
          ITYP = obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex)
          if ( ITYP == bufr_radarPrecip ) cycle BODY
          varLevel = vnl_varLevelFromVarnum(ityp)

          ! Interpolate the background-covariance statistics
          if ( obs_bodyElem_i(obsSpaceData, OBS_XTR, bodyIndex) /= 0 ) then
            IK=col_getNumLev(columnTrlOnAnlIncLev, varLevel)-1
            IPT  = IK + col_getOffsetFromVarno(columnTrlOnAnlIncLev, ityp)
            IPB  = IPT +1
            fge_uu = col_getElem(column, IPB, headerIndex, 'UU')
            fge_vv = col_getElem(column, IPB, headerIndex, 'VV')

          else
            ZLEV = obs_bodyElem_r(obsSpaceData, OBS_PPP, bodyIndex)
            IK   = obs_bodyElem_i(obsSpaceData, OBS_LYR, bodyIndex)
            IPT  = IK + col_getOffsetFromVarno(columnTrlOnAnlIncLev, ityp)
            IPB  = IPT+1
            ZPT  = col_getHeight(columnTrlOnAnlIncLev, IK, headerIndex, varLevel)
            ZPB  = col_getHeight(columnTrlOnAnlIncLev, IK+1, headerIndex, varLevel)
            ZWB  = (ZPT-ZLEV)/(ZPT-ZPB)
            ZWT  = 1.d0 - ZWB
            fge_uu = ZWB*col_getElem(column, IPB, headerIndex, 'UU') &
                   + ZWT*col_getElem(column, IPT, headerIndex, 'UU')
            fge_vv = ZWB*col_getElem(column, IPB, headerIndex, 'VV') &
                   + ZWT*col_getElem(column, IPT, headerIndex, 'VV')
          end if

          ! First-Guess Error Variance
          if ( cdfam == 'AL' )then
            ! Scan body indices for the azimuth
            bodyIndexStart = obs_headElem_i(obsSpaceData, OBS_RLN, headerIndex)
            bodyIndexEnd   = obs_headElem_i(obsSpaceData, OBS_NLV, headerIndex) &
                                  + bodyIndexStart - 1
            BODY_SUPP: do bodyIndex2 = bodyIndexStart, bodyIndexEnd
              if(obs_bodyElem_i(obsSpaceData, OBS_VNM, bodyIndex2) == 5021)then
                azimuth=obs_bodyElem_r(obsSpaceData, OBS_VAR, bodyIndex2) * MPC_RADIANS_PER_DEGREE_R8
                    exit BODY_SUPP
              end if
            end do BODY_SUPP

            fge_fam = sqrt((fge_vv*cos(azimuth))**2 + (fge_uu*sin(azimuth))**2)

          else if( cdfam == 'PR' )then
            fge_fam = ZWB*col_getElem(column, IPB, headerIndex) &
                    + ZWT*col_getElem(column, IPT, headerIndex)
                  
          else if( cdfam == 'RA' .and. ITYP == bufr_radvel ) then    

            ! Calculation of sigmap^2 = diag(H*P*H^t) to save sigmap in OBS_HPHT  
            ! 
            ! H includes vertical interpolation  
            !   and projection of U and V wind components along the direction of the beam
            !
            ! P forecast error covariance
            !  
            ! H = [ ZWT*sin(az) ZWT*cos(az) ZWB*sin(az) ZWB*cos(az) ] 
            ! diag(P) = [ sigmap_uuT^2  sigmap_vvT^2  sigmap_uuB^2  sigmap_vvB^2  ]
            sigmap_uuT = col_getElem(column, IPT, headerIndex, 'UU')
            sigmap_uuB = col_getElem(column, IPB, headerIndex, 'UU')
            sigmap_vvT = col_getElem(column, IPT, headerIndex, 'VV')
            sigmap_vvB = col_getElem(column, IPB, headerIndex, 'VV')

            fge_fam = sqrt((sigmap_uuT*ZWT*sin(azimuth))**2 + (sigmap_vvT*ZWT*cos(azimuth))**2 + &
                           (sigmap_uuB*ZWB*sin(azimuth))**2 + (sigmap_vvB*ZWB*cos(azimuth))**2)
                 
          else

            write(*,*)"ERROR:  The family", cdfam, " is not supported by setfgefamz"
            call utl_abort('setfgefamz')

          end if

          ! Store fge_fam in OBS_HPHT
          call obs_bodySet_r(obsSpaceData, OBS_HPHT, bodyIndex, fge_fam)
        
        end if

      end do BODY

    end do HEADER

  end subroutine setfgefamz

  !--------------------------------------------------------------------------
  ! setfgett
  !--------------------------------------------------------------------------
  subroutine setfgett(column,columnTrlOnAnlIncLev,lobsSpaceData)
    !
    !:Purpose: To interpolate vertically the contents of "column" to the
    !          pressure levels of the observations. Then to compute THE FIRST
    !          GUESS ERROR VARIANCES. A linear interpolation in ln(p) is
    !          performed.
    !
    implicit none

    ! Arguments:
    type(struct_columnData) :: column
    type(struct_columnData) :: columnTrlOnAnlIncLev
    type(struct_obs) :: lobsSpaceData

    ! Locals:
      INTEGER IPB,IPT
      INTEGER INDEX_HEADER,ITYP,IK
      INTEGER INDEX_BODY
      REAL*8 ZWB,ZWT
      REAL*8 ZLEV,ZPB,ZPT
      character(len=4) :: varLevel

      ! loop over all body rows
      BODY: do index_body=1,obs_numbody(lobsSpaceData)

         ! 1. Computation of sigmap
         !    ---------------------
         
         IF ( (obs_bodyElem_i(lobsSpaceData,OBS_ASS,index_body) == obs_assimilated) .and.    &
              (obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body).EQ. BUFR_NETT) ) THEN

            IF ( (obs_bodyElem_i(lobsSpaceData,OBS_XTR,index_body) .NE. 0) .and.    &
                 (obs_bodyElem_i(lobsSpaceData,OBS_VCO,index_body) .EQ. 2) ) THEN
               ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body)
               varLevel = vnl_varLevelFromVarnum(ityp)
               IK=col_getNumLev(columnTrlOnAnlIncLev,varLevel)-1
               INDEX_HEADER = obs_bodyElem_i(lobsSpaceData,OBS_HIND,index_body)
               IPT  = IK + col_getOffsetFromVarno(columnTrlOnAnlIncLev,ityp)
               IPB  = IPT +1
               call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,col_getElem(column,IPB,INDEX_HEADER))
            ELSE
               INDEX_HEADER = obs_bodyElem_i(lobsSpaceData,OBS_HIND,index_body)
               ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,index_body)
               IK   = obs_bodyElem_i(lobsSpaceData,OBS_LYR,index_body)
               IPT  = IK
               IPB  = IPT+1
               ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body)
               varLevel = vnl_varLevelFromVarnum(ityp)
               ZPT  = col_getPressure(columnTrlOnAnlIncLev,IK,INDEX_HEADER,varLevel)
               ZPB  = col_getPressure(columnTrlOnAnlIncLev,IK+1,INDEX_HEADER,varLevel)
               ZWB  = LOG(ZLEV/ZPT)/LOG(ZPB/ZPT)
               ZWT  = 1.0D0 - ZWB

               ! FIRST GUESS ERROR VARIANCE

               call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,   &
                  (ZWB*col_getElem(column,IPB,INDEX_HEADER,'TT') + ZWT*col_getElem(column,IPT,INDEX_HEADER,'TT')))
            ENDIF

         ENDIF

      END DO BODY

      RETURN
  end subroutine setfgett

  !--------------------------------------------------------------------------
  ! setfgeSurf
  !--------------------------------------------------------------------------
  subroutine setfgeSurf(column, columnTrlOnAnlIncLev, lobsSpaceData)
    !
    !:Purpose: To interpolate vertically the contents of "column" to the
    !          pressure levels of the observations. A linear interpolation in
    !          ln(p) is performed.
    !
    implicit none

    ! Arguments
    type(struct_columnData) :: column
    type(struct_columnData) :: columnTrlOnAnlIncLev
    type(struct_obs)        :: lobsSpaceData

    ! Locals
    integer          :: ipb, ipt, idim, headerIndex, ik, bodyIndex, ityp, bodyElem_i
    real(8)          :: zwb, zwt, zlev, zpt, zpb, zhhh, bodyElem_r, colElem1, colElem2
    character(len=2) :: cfam
    character(len=4) :: varLevel
    character(len=12) :: stnid
    logical          :: ok

    ! loop over all body rows
    BODY: do bodyIndex = 1, obs_numbody( lobsSpaceData )

      cfam = obs_getFamily( lobsSpaceData, bodyIndex = bodyIndex )
      if( cfam == 'SF'.or. cfam == 'TM' .or. cfam == 'UA' .or. cfam  == 'SC' .or. cfam == 'GP' .or. cfam == 'GL' ) then

        ! Process all data within the domain of the model (excluding GB-GPS ZTD data)
        ok = .false.

        if ( obs_bodyElem_i( lobsSpaceData, OBS_VCO, bodyIndex ) == 1 ) then

          ityp = obs_bodyElem_i( lobsSpaceData, OBS_VNM, bodyIndex )
          if ( ityp == BUFR_NETS .or. ityp == BUFR_NEPS .or. ityp == BUFR_NEPN .or. ityp == BUFR_NESS .or. &
             ityp == BUFR_NEUS .or. ityp == BUFR_NEVS .or. ityp == BUFR_NEFS .or. ityp == BUFR_NEDS .or. &
             ityp == bufr_sst  .or. ityp == BUFR_ICEC .or. ityp == bufr_logVis  .or. ityp == bufr_gust .or. &
             ityp == bufr_riverFlow) then

            ok = ( obs_bodyElem_i( lobsSpaceData, OBS_ASS, bodyIndex ) == obs_assimilated )

          else if ( ityp == BUFR_NEZD ) then

            ! make sure total zenith delay (from ground-based GPS) not treated
            ok=.false.

          else

            ok=(obs_bodyElem_i( lobsSpaceData, OBS_ASS, bodyIndex ) == obs_assimilated .and. &
                obs_bodyElem_i( lobsSpaceData, OBS_XTR, bodyIndex ) >= 0)
            if ( ok ) write(*,*) 'setfgesurf: WARNING!!! unknown obs seen'
            if ( ok ) write(*,*) 'setfgesurf: ityp=',ityp,', cfam=',cfam

          end if

          if ( ok ) then

            headerIndex = obs_bodyElem_i( lobsSpaceData, OBS_HIND, bodyIndex )
            ityp         = obs_bodyElem_i( lobsSpaceData, OBS_VNM , bodyIndex )
            varLevel     = vnl_varLevelFromVarnum( ityp )
            idim = 1
            if ( varLevel == 'SF') idim = 0
            ik   = obs_bodyElem_i( lobsSpaceData, OBS_LYR, bodyIndex )
            zlev = obs_bodyElem_r( lobsSpaceData, OBS_PPP, bodyIndex )
            zhhh = zlev

            if ( ityp == BUFR_NETS .or. ityp == BUFR_NEPS .or. ityp == BUFR_NEPN .or. &
              ityp == BUFR_NESS .or. ityp == BUFR_NEUS .or. ityp == BUFR_NEVS .or. &
              ityp == bufr_logVis  .or. ityp == bufr_gust) then

              ipt  = ik + col_getOffsetFromVarno(columnTrlOnAnlIncLev,ityp)
              ipb  = ipt+1
              call obs_bodySet_r( lobsSpaceData, OBS_HPHT, bodyIndex, col_getElem( column, ipb, headerIndex ) )

            else

              ipt  = ik + col_getOffsetFromVarno(columnTrlOnAnlIncLev,ityp)
              ipb  = ipt+1
              zpt  = col_getHeight(columnTrlOnAnlIncLev,ik,headerIndex,varLevel)
              zpb  = col_getHeight(columnTrlOnAnlIncLev,ik+1,headerIndex,varLevel)
              zwb  = idim*(zpt-zhhh)/(zpt-zpb)
              zwt  = 1.d0 - zwb

              if ( obs_bodyElem_i( lobsSpaceData, OBS_XTR, bodyIndex ) == 0 ) then

                call obs_bodySet_r( lobsSpaceData, OBS_HPHT, bodyIndex,   &
                zwb * col_getElem( column, ipb, headerIndex ) + zwt * col_getElem( column, ipt, headerIndex ))

              else

                call obs_bodySet_r( lobsSpaceData, OBS_HPHT, bodyIndex,   &
                  col_getElem( column, ik + col_getOffsetFromVarno( columnTrlOnAnlIncLev, ityp ), headerIndex ))

              end if

              stnid = obs_elem_c( lobsSpaceData, 'STID', headerIndex )
              if(stnid == '99999999' ) then

                bodyElem_i = obs_bodyElem_i( lobsSpaceData, OBS_XTR, bodyIndex )
                write(*,*) 'setfgesurf: stn, ityp, xtr, ipt, ipb, zwt, zwb',  &
                   stnid, ityp, &
                   bodyElem_i, ipt, ipb, zwt, zwb
                bodyElem_r = obs_bodyElem_i( lobsSpaceData, OBS_HPHT, bodyIndex )
                colElem1 = col_getElem( column, ipb, headerIndex )
                colElem2 = col_getElem( column, ipt, headerIndex )
                write(*,*) 'setfgesurf: gobs(ipb), gobs(ipt), fge',   &
                    colElem1, colElem2, bodyElem_r

              endif

            end if
          end if
        end if

      end if

    end do BODY

  end subroutine setfgeSurf

  !--------------------------------------------------------------------------
  ! setfgedif
  !--------------------------------------------------------------------------
  subroutine setfgedif(cdfam,columnTrlOnAnlIncLev,lobsSpaceData)
    !
    !:Purpose: To construct the FIRST GUESS ERROR VARIANCES from the
    !          diff-calculated dependencies and the primary errors.
    !
    implicit none

    ! Arguments:
    character*2 cdfam
    type(struct_columnData) :: columnTrlOnAnlIncLev
    type(struct_obs)        :: lobsSpaceData

    ! Locals:
      INTEGER INDEX_HEADER, IDATYP, INDEX_BODY, iProfile, varNum
      REAL*8 zLat, Lat, sLat
      REAL*8 zLon, Lon
      REAL*8 zAzm !, Azm
      INTEGER ISAT
      REAL*8 Rad, Geo
      REAL*8, allocatable :: zPP(:)
      REAL*8, allocatable :: zDP(:)
      REAL*8, allocatable :: zTT(:)
      REAL*8, allocatable :: zHU(:)
      REAL*8, allocatable :: zUU(:)
      REAL*8, allocatable :: zVV(:)
      INTEGER status
      INTEGER JL
      REAL*8 ZP0, ZMT
      REAL*8 ZFGE, ZERR
      INTEGER JV, NGPSLEV, NWNDLEV
      LOGICAL  ASSIM, LFIRST, FIRSTHEADER

      INTEGER NH, NH1

      ! REAL*8 JAC(ngpscvmx)
      REAL*8 DV (ngpscvmx)
      TYPE(GPS_PROFILE)           :: PRF
      REAL*8       , allocatable :: H   (:),AZMV(:)
      TYPE(GPS_DIFF), allocatable :: RSTV(:)
      type(struct_vco), pointer  :: vco_anl
      real*8, dimension(:), pointer :: dPdPs

      WRITE(*,*)'ENTER SETFGEDIFF'

      ! Initializations

      nullify(dPdPs)

      NGPSLEV=col_getNumLev(columnTrlOnAnlIncLev,'TH')
      NWNDLEV=col_getNumLev(columnTrlOnAnlIncLev,'MM')
      LFIRST=.FALSE.
      if ( .NOT.allocated(gps_vRO_Jacobian) ) then
         LFIRST = .TRUE.
         allocate(zPP (NGPSLEV))
         allocate(zDP (NGPSLEV))
         allocate(zTT (NGPSLEV))
         allocate(zHU (NGPSLEV))
         allocate(zUU (NGPSLEV))
         allocate(zVV (NGPSLEV))

         allocate(gps_vRO_Jacobian(gps_numROProfiles,GPSRO_MAXPRFSIZE,2*NGPSLEV+1))
         allocate(gps_vRO_lJac    (gps_numROProfiles))
         gps_vRO_lJac=.false.

         allocate( H    (GPSRO_MAXPRFSIZE) )
         allocate( AZMV (GPSRO_MAXPRFSIZE) )
         allocate( RSTV (GPSRO_MAXPRFSIZE) )
      endif

      vco_anl => col_getVco(columnTrlOnAnlIncLev)


      ! Loop over all header indices of the 'RO' family:

      call obs_set_current_header_list(lobsSpaceData,CDFAM)
      FIRSTHEADER=.TRUE.
      HEADER: do
         index_header = obs_getHeaderIndex(lobsSpaceData)
         if (index_header < 0) exit HEADER

         ! Process only refractivity data (codtyp 169)

         IDATYP = obs_headElem_i(lobsSpaceData,OBS_ITY,INDEX_HEADER)
         IF ( IDATYP .EQ. 169 ) THEN
            iProfile = gps_iprofile_from_index(index_header)
            varNum = gps_vRO_IndexPrf(iProfile, 2)

            ! Scan for requested data values of the profile, and count them
            
            ASSIM = .FALSE.
            NH = 0
            call obs_set_current_body_list(lobsSpaceData, INDEX_HEADER)
            BODY: do
               INDEX_BODY = obs_getBodyIndex(lobsSpaceData)
               if (INDEX_BODY < 0) exit BODY
               IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY) == obs_assimilated ) THEN
                  ASSIM = .TRUE.
                  NH = NH + 1
               ENDIF
            ENDDO BODY

            ! If assimilations are requested, prepare and apply the observation operator

            IF (ASSIM) THEN
               iProfile=gps_iprofile_from_index(INDEX_HEADER)

               ! Profile at the observation location:

               if (.not.gps_vRO_lJac(iProfile)) then

                  ! Basic geometric variables of the profile:

                  zLat = obs_headElem_r(lobsSpaceData,OBS_LAT,INDEX_HEADER)
                  zLon = obs_headElem_r(lobsSpaceData,OBS_LON,INDEX_HEADER)
                  ISAT = obs_headElem_i(lobsSpaceData,OBS_SAT,INDEX_HEADER)
                  Rad  = obs_headElem_r(lobsSpaceData,OBS_TRAD,INDEX_HEADER)
                  Geo  = obs_headElem_r(lobsSpaceData,OBS_GEOI,INDEX_HEADER)
                  zAzm = obs_headElem_r(lobsSpaceData,OBS_AZA,INDEX_HEADER) / MPC_DEGREES_PER_RADIAN_R8
                  zMT  = col_getHeight(columnTrlOnAnlIncLev,NGPSLEV,INDEX_HEADER,'TH')
                  Lat  = zLat * MPC_DEGREES_PER_RADIAN_R8
                  Lon  = zLon * MPC_DEGREES_PER_RADIAN_R8
                  !Azm  = zAzm * MPC_DEGREES_PER_RADIAN_R8
                  sLat = sin(zLat)
                  zMT  = zMT * ec_rg / gpsgravitysrf(sLat)
                  zP0  = col_getElem(columnTrlOnAnlIncLev,1,INDEX_HEADER,'P0')

                  ! approximation for dPdPs               
                  if (associated(dPdPs)) then
                    deallocate(dPdPs,stat=status)
                    nullify(dPdPs)
                  end if
                  status = vgd_dpidpis(vco_anl%vgrid,vco_anl%ip1_T,dPdPs,zP0)
                  zDP(1:NGPSLEV) = dPdPs(1:NGPSLEV)

                  DO JL = 1, NGPSLEV

                     ! Profile x

                     zPP(JL) = col_getPressure(columnTrlOnAnlIncLev,JL,INDEX_HEADER,'TH')
                     zTT(JL) = col_getElem(columnTrlOnAnlIncLev,JL,INDEX_HEADER,'TT') - p_TC
                     zHU(JL) = col_getElem(columnTrlOnAnlIncLev,JL,INDEX_HEADER,'HU')
                     zUU(JL) = 0.d0
                     zVV(JL) = 0.d0
                  ENDDO
                  DO JL = 1, NWNDLEV
                     zUU(JL) = col_getElem(columnTrlOnAnlIncLev,JL,INDEX_HEADER,'UU')
                     zVV(JL) = col_getElem(columnTrlOnAnlIncLev,JL,INDEX_HEADER,'VV')
                  ENDDO
                  zUU(NGPSLEV) = zUU(NWNDLEV)
                  zVV(NGPSLEV) = zUU(NWNDLEV)
  
                  ! GPS profile structure:

                  call gps_struct1sw(ngpslev,zLat,zLon,zAzm,zMT,Rad,geo,zP0,zPP,zDP,zTT,zHU,zUU,zVV,prf)

                  ! Prepare the vector of all the observations:

                  NH1 = 0
                  call obs_set_current_body_list(lobsSpaceData, INDEX_HEADER)
                  BODY_2: do
                     INDEX_BODY = obs_getBodyIndex(lobsSpaceData)
                     if (INDEX_BODY < 0) exit BODY_2
                     IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY) == obs_assimilated ) THEN
                        NH1      = NH1 + 1
                        H(NH1)   = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)
                        AZMV(NH1)= zAzm
                     ENDIF
                  ENDDO BODY_2

                  ! Apply the observation operator:

                  IF (varNum == bufr_nebd) THEN
                     CALL GPS_BNDOPV1(H, AZMV, NH, PRF, RSTV)
                  ELSE
                     CALL GPS_REFOPV (H,       NH, PRF, RSTV)
                  ENDIF
                  DO NH1=1,NH
                     gps_vRO_Jacobian(iProfile,NH1,:)= RSTV(NH1)%DVAR(1:2*NGPSLEV+1)
                  ENDDO
                  gps_vRO_lJac(iProfile)=.true.
               endif

               ! Local error

               DO JL = 1, NGPSLEV
                  DV (        JL) = 1.d0
                  DV (NGPSLEV+JL) = 1.d0
               ENDDO
               DV (2*NGPSLEV+1)   = 2.d0

               ! Perform the H(xb)DV operation:

               NH1 = 0
               call obs_set_current_body_list(lobsSpaceData, INDEX_HEADER)
               BODY_3: do
                  INDEX_BODY = obs_getBodyIndex(lobsSpaceData)
                  if (INDEX_BODY < 0) exit BODY_3
                  IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY) == obs_assimilated ) THEN
                     NH1 = NH1 + 1

                     ! Observation jacobian

!                     JAC = RSTV(NH1)%DVAR

                     ! Evaluate sqrt( H(xb)DV **2 )

                     ZFGE = 0.d0
                     DO JV = 1, 2*PRF%NGPSLEV+1
                        ZFGE = ZFGE + (gps_vRO_Jacobian(iProfile,NH1,JV) * DV(JV))**2
                     ENDDO
                     ZFGE = SQRT(ZFGE)
                     ZERR = obs_bodyElem_r(lobsSpaceData,OBS_OER,INDEX_BODY)
    
                     ! FIRST GUESS ERROR VARIANCE

                     call obs_bodySet_r(lobsSpaceData,OBS_HPHT,INDEX_BODY,ZFGE)
                     IF (FIRSTHEADER) THEN
11                      FORMAT(A12,2I5,F12.2,3F16.8)
                        WRITE(*,11)'SETFGEDIFFGE',NH1,NH,H(NH1),RSTV(NH1)%VAR,ZFGE,ZERR
                     ENDIF
                  ENDIF
               ENDDO BODY_3
            ENDIF
         ENDIF
         FIRSTHEADER = .FALSE.
      ENDDO HEADER

      IF (LFIRST) THEN
         deallocate( RSTV )
         deallocate( AZMV )
         deallocate( H    )

         deallocate(zVV)
         deallocate(zUU)
         deallocate(zHU)
         deallocate(zTT)
         deallocate(zDP)
         deallocate(zPP)
         deallocate(gps_vRO_Jacobian)
      ENDIF

      WRITE(*,*)'EXIT SETFGEDIFF'
      RETURN
  end subroutine setfgedif

  !--------------------------------------------------------------------------
  ! setfgegps
  !--------------------------------------------------------------------------
  subroutine setfgegps(column,columnTrlOnAnlIncLev,lobsSpaceData)
    !
    !:Purpose: To set FGE for all GPS ZTD observations using Jacobians from ZTD
    !          observation operator
    !
    !:Option: Test ZTD operators (compares H(x+dx)-H(x) with (dH/dx)*dx
    !         when LTESTOP = .true.)
    !
    !:Note:
    !      _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    !                             9 October 2015                               
    !                                                                          
    !          :NOTE: Effective Rev644M, this routine is no longer used!       
    !                 FGE for ZTD is no longer needed for background check.    
    !                 Routine is called only when LTESTOP=.true., in which     
    !                 case the operator test only is done.                     
    !                                                                          
    !      _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    !
    implicit none

    ! Arguments:
    type(struct_columnData) :: column
    type(struct_columnData) :: columnTrlOnAnlIncLev
    type(struct_obs) :: lobsSpaceData

    ! Locals:
    ! column  contains background errors for control variables on model levels
    ! columnTrlOnAnlIncLev contains lo-res first guess profiles at obs locations
      type(struct_vco), pointer :: vco_anl
      REAL*8 ZLAT, Lat
      REAL*8 ZLON, Lon
      REAL*8, allocatable :: ZPP(:)
      REAL*8, allocatable :: ZDP(:)
      REAL*8, allocatable :: ZTT(:)
      REAL*8, allocatable :: ZHU(:)
      REAL*8, allocatable :: zHeight(:)
      REAL*8, allocatable :: zHeight2(:)
      REAL*8, allocatable :: ZTTB(:)
      REAL*8, allocatable :: ZHUB(:)
      REAL*8, allocatable :: ZQQB(:)
      REAL*8, allocatable :: ZQQ(:)
      REAL*8, allocatable :: ZTTB_P(:)
      REAL*8, allocatable :: ZQQB_P(:)
      REAL*8, allocatable :: zHeight_P(:)
      REAL*8, allocatable :: ZPP_P(:)
      
      REAL*8 ZP0
      REAL*8 ZP0B, ZP0B_P
      REAL*8 ZMT, ZTOP, ZBOT
 
      REAL*8 JAC(ngpscvmx)
      REAL*8 DX (ngpscvmx)

      REAL*8 ZLEV, ZTDOBS, ZPSMOD
      REAL*8 ZLSUM
      REAL*8 DELTAH_NL, DELTAH_TL
      REAL*8 PERTFAC, ZTDM
      REAL*8 ZDZMIN, ZSUMTEST

      INTEGER INDEX_HEADER, FIRST_HEADER
      INTEGER IDATYP, ITYP
      INTEGER IDATA, IDATEND, INDEX_BODY
      INTEGER JL, NFLEV_T, ILYR, IOBS
      INTEGER INOBS_OPT, INOBS_JAC, icount, status, iversion

      LOGICAL  ASSIM, OK, LSTAG
      CHARACTER*9  STN_JAC
      character(len=12) :: stnid
      
      CHARACTER(len=4) :: varLevel
      
      TYPE(GPS_PROFILEZD)    :: PRF, PRFP
      TYPE(GPS_DIFF)         :: ZTDopv, ZTDopvP

      real*8, dimension(:), pointer :: dPdPs

      IF (numGPSZTD .EQ. 0) RETURN

      ! Initializations

      nullify(dPdPs)

      NFLEV_T = col_getNumLev(columnTrlOnAnlIncLev,'TH')
      allocate(ZPP(NFLEV_T))
      allocate(ZDP(NFLEV_T))
      allocate(ZTT(NFLEV_T))
      allocate(ZHU(NFLEV_T))
      allocate(zHeight(NFLEV_T))
      allocate(ZTTB(NFLEV_T))
      allocate(ZHUB(NFLEV_T))
      allocate(ZQQB(NFLEV_T))
      allocate(ZQQ(NFLEV_T))

      ! Number of locations/sites for observation operator test
      INOBS_OPT = 50
      ! Number of locations/sites for Jacobian printout
      INOBS_JAC  = 5
      ! Factor to multiply background errors for perturbation vector
      PERTFAC = 0.75d0

      STN_JAC = 'FSL_BRFT '
      
      ZDZMIN = DZMIN

      vco_anl => col_getVco(columnTrlOnAnlIncLev)
      status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=iversion)
      if (iversion .eq. 5002) then
         LSTAG = .TRUE. 
         WRITE(*,*)'VERTICAL COORD OF ANALYSIS FIELDS IS STAGGERED'
         WRITE(*,*)'VCODE= ',iversion,' LSTAG= ',LSTAG
      else
         LSTAG = .FALSE.
         WRITE(*,*)'VERTICAL COORD OF ANALYSIS FIELDS IS NOT STAGGERED'
         WRITE(*,*)'VCODE= ',iversion,' LSTAG= ',LSTAG
      endif

      IF ( .NOT.LTESTOP ) THEN

      first_header=-1
      icount = 0

      ! loop over all header indices of the 'GP' family
      call obs_set_current_header_list(lobsSpaceData,'GP')
      HEADER: do
         index_header = obs_getHeaderIndex(lobsSpaceData)
         if (index_header < 0) exit HEADER
         if (first_header .eq. -1) first_header = index_header
     
               ! Process only zenith delay data (codtyp 189 and BUFR_NEZD)

               IDATYP = obs_headElem_i(lobsSpaceData,OBS_ITY,INDEX_HEADER)
               IF ( IDATYP .EQ. 189 ) THEN

                  ! Loop over data in the observations

                  IDATA   = obs_headElem_i(lobsSpaceData,OBS_RLN,INDEX_HEADER)
                  IDATEND = obs_headElem_i(lobsSpaceData,OBS_NLV,INDEX_HEADER) + IDATA - 1
                  ASSIM = .FALSE.

                  ! Scan for requested assimilations, and count them.

                  DO INDEX_BODY= IDATA, IDATEND
                     ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,INDEX_BODY)
                     OK = ( (ITYP .EQ. BUFR_NEZD) .AND. (obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY) == obs_assimilated) )
                     IF ( OK ) THEN
                        ASSIM = .TRUE.
                        ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)
                        icount = icount + 1
                     ENDIF
                  ENDDO

                  ! If assimilations are requested, apply the AD observation operator

                  IF (ASSIM) THEN
    
                     ! LR background profile and background errors at the observation location x :

                     Lat  = obs_headElem_r(lobsSpaceData,OBS_LAT,INDEX_HEADER)
                     Lon  = obs_headElem_r(lobsSpaceData,OBS_LON,INDEX_HEADER)
                     ZLAT = Lat * MPC_DEGREES_PER_RADIAN_R8
                     ZLON = Lon * MPC_DEGREES_PER_RADIAN_R8
                     ZP0B = col_getElem(columnTrlOnAnlIncLev,1,INDEX_HEADER,'P0')
                     DO JL = 1, NFLEV_T
                       ZPP(JL)  = col_getPressure(columnTrlOnAnlIncLev,JL,INDEX_HEADER,'TH')
                       ZTTB(JL) = col_getElem(columnTrlOnAnlIncLev,JL,INDEX_HEADER,'TT')- 273.15d0
                       ZTT(JL)  = col_getElem(column,JL,INDEX_HEADER,'TT')
                       DX(JL)   = ZTT(JL)
                       ZHUB(JL) = col_getElem(columnTrlOnAnlIncLev,JL,INDEX_HEADER,'HU')
                       ZQQB(JL) = ZHUB(JL)
                       ZHU(JL)  = col_getElem(column,JL,INDEX_HEADER,'HU')
                       DX(NFLEV_T+JL) = ZHU(JL)
                       zHeight(JL)  = col_getHeight(columnTrlOnAnlIncLev,JL,INDEX_HEADER,'TH')
                       DX(2*NFLEV_T+JL) = col_getHeight(column,JL,INDEX_HEADER,'TH')
                     ENDDO
                     ZP0  = col_getElem(column,1,INDEX_HEADER,'P0')
                     DX(3*NFLEV_T+1) = ZP0
                     ZMT  = zHeight(NFLEV_T)
                     CALL gps_structztd_v2(NFLEV_T,Lat,Lon,ZMT,ZP0B,ZPP,ZTTB,ZHUB,zHeight,LBEVIS,IREFOPT,PRF)
                     CALL gps_ztdopv(ZLEV,PRF,LBEVIS,ZDZMIN,ZTDopv,ZPSMOD,IZTDOP)
                     JAC = ZTDopv%DVar

                     DO INDEX_BODY= IDATA, IDATEND
                        ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,INDEX_BODY)
                        IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY) == obs_assimilated .AND. ITYP.EQ.BUFR_NEZD ) THEN

                           ! Observation error    SDERR
!                           ZOER = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)

                           ! Observation height (m)
                           ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)

                           ZLSUM  = 0.0d0

                           DO JL = 1, 3*NFLEV_T+1
                             ZLSUM = ZLSUM + (JAC(JL)*DX(JL))**2
                           ENDDO
                           call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,SQRT(ZLSUM))

                           IF (icount .LE. INOBS_JAC) THEN
!                           IF ( obs_elem_c(lobsSpaceData,'STID',INDEX_HEADER) .EQ. STN_JAC ) THEN
                             stnid = obs_elem_c(lobsSpaceData,'STID',INDEX_HEADER)
                             WRITE(*,'(A11,A9)') 'SETFGEGPS: ', stnid
                             WRITE(*,*) '  ZTD, ZTD FGE = ', ZTDopv%Var, SQRT(ZLSUM)
                             WRITE(*,'(A11,A9,3(1x,f7.2))')   &
                               'SETFGEGPS: ',stnid,ZLAT,ZLON,ZLEV
                             WRITE(*,*) 'JL JACT JACQ FGE_T FGE_LQ QQ'
                             DO JL = 1, NFLEV_T
                               WRITE(*,'(1X,I2,5(1x,E13.6))') JL,JAC(JL),JAC(JL+NFLEV_T)/ZQQB(JL),ZTT(JL),ZHU(JL),ZQQB(JL)
                             ENDDO                         
                             WRITE(*,*) 'JACPS FGE_PS'
                             WRITE(*,'(2(1x,E13.6))') JAC(3*NFLEV_T+1), ZP0
                           ENDIF

                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF

      ENDDO HEADER
      
      ENDIF

      !-------------------------------------------------------------------------

      IF ( LTESTOP ) THEN
      
      allocate(ZTTB_P(NFLEV_T))
      allocate(ZQQB_P(NFLEV_T))
      allocate(zHeight2(NFLEV_T))
      allocate(zHeight_P(NFLEV_T))
      allocate(ZPP_P(NFLEV_T))

      icount = 0
      ZSUMTEST = 0
      
      ! loop over all header indices of the 'GP' family
      call obs_set_current_header_list(lobsSpaceData,'GP')
      HEADER2: DO
         index_header = obs_getHeaderIndex(lobsSpaceData)
         if (index_header < 0) exit HEADER2
         if (icount > INOBS_OPT ) exit HEADER2

         ! Loop over data in the observations

         IDATA   = obs_headElem_i(lobsSpaceData,OBS_RLN,INDEX_HEADER)
         IDATEND = obs_headElem_i(lobsSpaceData,OBS_NLV,INDEX_HEADER) + IDATA - 1

         ! LR background profile and background errors at the observation location x :

         Lat  = obs_headElem_r(lobsSpaceData,OBS_LAT,INDEX_HEADER)
         Lon  = obs_headElem_r(lobsSpaceData,OBS_LON,INDEX_HEADER)
         ZLAT = Lat * MPC_DEGREES_PER_RADIAN_R8
         ZLON = Lon * MPC_DEGREES_PER_RADIAN_R8
         ZP0B = col_getElem(columnTrlOnAnlIncLev,1,INDEX_HEADER,'P0')

         ! approximation for dPdPs               
         if (associated(dPdPs)) then
           deallocate(dPdPs,stat=status)
           nullify(dPdPs)
         end if
         status = vgd_dpidpis(vco_anl%vgrid,vco_anl%ip1_T,dPdPs,ZP0B)
         zDP(1:NFLEV_T) = dPdPs(1:NFLEV_T)

         DO JL = 1, NFLEV_T
            ZPP(JL)  = col_getPressure(columnTrlOnAnlIncLev,JL,INDEX_HEADER,'TH')
            ZTTB(JL) = col_getElem(columnTrlOnAnlIncLev,JL,INDEX_HEADER,'TT')- 273.15d0
            ZTT(JL)  = col_getElem(column,JL,INDEX_HEADER,'TT') * PERTFAC
            ZQQB(JL) = col_getElem(columnTrlOnAnlIncLev,JL,INDEX_HEADER,'HU')
            ZQQ(JL)  = col_getElem(column,JL,INDEX_HEADER,'HU') * PERTFAC
            zHeight(JL)  = col_getHeight(columnTrlOnAnlIncLev,JL,INDEX_HEADER,'TH')
            zHeight2(JL)  = col_getHeight(column,JL,INDEX_HEADER,'TH') * PERTFAC
         ENDDO
         ZP0  = col_getElem(column,1,INDEX_HEADER,'P0') * PERTFAC
         ZMT  = zHeight(NFLEV_T)

         DO JL = 1, NFLEV_T
             DX (      JL) = ZTT(JL)
             DX (NFLEV_T+JL) = ZQQ(JL)
             DX (2*NFLEV_T+JL) = zHeight2(JL)
         ENDDO
         DX (3*NFLEV_T+1) = ZP0

         ZTDOBS = -1.0d0
         DO INDEX_BODY = IDATA, IDATEND
           ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,INDEX_BODY)
           IOBS = obs_bodyElem_i(lobsSpaceData,OBS_HIND,INDEX_BODY)
           IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY) == obs_assimilated .AND. ITYP .EQ. BUFR_NEZD ) THEN
             varLevel = vnl_varLevelFromVarnum(ITYP)
             ZTDOBS  = obs_bodyElem_r(lobsSpaceData,OBS_VAR,INDEX_BODY)
             ZLEV    = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)
             ILYR    = obs_bodyElem_i(lobsSpaceData,OBS_LYR,INDEX_BODY)
             ZTOP    = col_getHeight(columnTrlOnAnlIncLev,ILYR,IOBS,varLevel)
             if ( ILYR .LT. NFLEV_T ) then
               ZBOT    = col_getHeight(columnTrlOnAnlIncLev,ILYR+1,IOBS,varLevel)
             else
               ZBOT    = ZTOP
             endif
             icount  = icount + 1
           ENDIF
         ENDDO

         IF ( ZTDOBS .GT. 0.d0 ) THEN
           ! Create the pertubation control vector
           DO JL = 1, NFLEV_T
             ZPP_P(JL)  = ZPP(JL)  + ZDP(JL)*ZP0
             ZTTB_P(JL) = ZTTB(JL) + ZTT(JL)
             ZQQB_P(JL) = ZQQB(JL) + ZQQ(JL)
             zHeight_P(JL) = zHeight(JL) + zHeight2(JL)
           ENDDO
           ZP0B_P = ZP0B + ZP0

           ! Non-linear observation operator --> delta_H = H(x+delta_x) - H(x)

           CALL gps_structztd_v2(NFLEV_T,Lat,Lon,ZMT,ZP0B,ZPP,ZTTB,ZQQB,zHeight,LBEVIS,IREFOPT,PRF)
           CALL gps_structztd_v2(NFLEV_T,Lat,Lon,ZMT,ZP0B_P,ZPP_P,ZTTB_P,ZQQB_P,zHeight_P,LBEVIS,IREFOPT,PRFP)
           CALL gps_ztdopv(ZLEV,PRF,LBEVIS,ZDZMIN,ZTDopv,ZPSMOD,IZTDOP)
           JAC  = ZTDopv%DVar
           ZTDM = ZTDopv%Var
           CALL gps_ztdopv(ZLEV,PRFP,LBEVIS,ZDZMIN,ZTDopvP,ZPSMOD,IZTDOP)
           DELTAH_NL = ZTDopvP%Var - ZTDopv%Var

           ! Linear  --> delta_H = dH/dx * delta_x

           DELTAH_TL = 0.0d0
           DO JL = 1, 3*NFLEV_T+1
             DELTAH_TL = DELTAH_TL + JAC(JL)*DX(JL)
           ENDDO

           stnid = obs_elem_c(lobsSpaceData,'STID',INDEX_HEADER)
           WRITE(*,*) 'SETFGEGPS: GPS ZTD OBSOP TEST FOR SITE ', stnid
           WRITE(*,*) ' '
           WRITE(*,*) '  DZ (M), MODEL LEVEL ABOVE = ', ZLEV-ZMT, ILYR
           WRITE(*,*) '  ZLEV (M), ZTOP (M), ZBOT (M) = ', ZLEV, ZTOP, ZBOT
           WRITE(*,*) '  ZTD OBS (MM)            = ', ZTDOBS*1000.d0
           WRITE(*,*) '  ZTD_MOD                 = ', ZTDM*1000.d0
           WRITE(*,*) '  DELTAH_NL, DELTAH_TL = ', DELTAH_NL*1000.d0, DELTAH_TL*1000.d0
           WRITE(*,*) ' '
           WRITE(*,*) '  DELTAH_TL/DELTAH_NL = ', DELTAH_TL/DELTAH_NL
           WRITE(*,*) ' '  
           
           ZSUMTEST = ZSUMTEST + (DELTAH_TL/DELTAH_NL)
           
         ENDIF

      ENDDO HEADER2
      
      WRITE(*,*) ' '
      WRITE(*,*) 'SETFGEGPS: ----- GPS ZTD OBSOP TEST SUMMARY -----'
      WRITE(*,*) '           NUMBER OF TESTS (sites) = ', icount
      WRITE(*,*) '           AVG DELTAH_TL/DELTAH_NL = ', ZSUMTEST/FLOAT(icount)
      WRITE(*,*) ' '  

      deallocate(ZTTB_P)
      deallocate(ZQQB_P)
      deallocate(zHeight2)
      deallocate(zHeight_P)
      deallocate(ZPP_P)

      ENDIF
      !-------------------------------------------------------------------------

      deallocate(ZPP)
      deallocate(ZDP)
      deallocate(ZTT)
      deallocate(ZHU)
      deallocate(zHeight)
      deallocate(ZTTB)
      deallocate(ZHUB)
      deallocate(ZQQB)
      deallocate(ZQQ)

      RETURN
  end subroutine setfgegps

  !------------------ CH obs family OmP error std dev routines --------------
  
  !--------------------------------------------------------------------------
  ! ose_setOmPstddevCH
  !--------------------------------------------------------------------------
  character(len=4) function ose_setOmPstddevCH(obsSpaceData) result(availableOMPE)
    !
    !:Purpose: To read OmP error std dev from auxiliary file or calculate from OmP.
    !    
    implicit none

    ! Arguments:
    type(struct_obs) :: obsSpaceData ! observation-space data; output saved in OBS_OMPE column

    ! Externals:
    ! Local:
    integer, parameter :: ndim=1

    ! Check for the presence of CH observations
    availableOMPE = '    '
    if ( .not.obs_famExist(obsSpaceData,'CH',localMPI_opt=.true.) ) return
            
    ! read the OmP error std. dev. information from the auxiliary file
    call ose_readOmPstddev_auxfileCH
    
    availableOMPE='None'
    if ( OmPstdCH%n_stnid == 0 ) return ! All CH family OBS_OMPE to be estimated elsewhere via use of OBS_HPHT and OBS_OER

    ! Calc from the OmP dataset if requested (and there is enough data)
    if ( any(OmPstdCH%source(1:OmPstdCH%n_stnid) == 1) )  call ose_calcOmPstddevCH(obsSpaceData)

    ! Assign OmP error std dev values to OBS_OMPE for CH family
    if ( any(OmPstdCH%source(1:OmPstdCH%n_stnid) <= 1) ) call ose_fillOmPstddevCH(obsSpaceData)
    
    ! Deallocate OmPstdCH space
    call ose_deallocOmPstddevCH
    
    ! Check if ALL CH family obs have been assigned usable values in OBS_OMPE
    if ( .not. ose_OmPstddevExistsForAllCH(ObsSpaceData) ) then
      ! For some obs, need to calc OBS_HPHT to get OBS_OMPE
      availableOMPE='Some'
    else
      ! All OBS_OMPE are available
      availableOMPE='All '
    end if
    
  end function ose_setOmPstddevCH
  
  !--------------------------------------------------------------------------
  ! ose_readOmPstddev_auxfileCH
  !--------------------------------------------------------------------------
  subroutine ose_readOmPstddev_auxfileCH
    !
    !:Purpose:  To read and store OmP error std. dev. as needed for CH
    !           family obs - if/when available.
    !
    implicit none

    integer, external :: fnom, fclos
    integer :: ierr, nulstat
    logical :: LnewExists
  
    character (len=128) :: ligne
    character(len=11) :: AuxObsDataFileCH = 'obsinfo_chm'
    integer :: ipos, stnidIndex, monthIndex, levIndex, ios, isize, icount
    character(len=20) :: abortText

    ! Initialization

    OmPstdCH%n_stnid=0

    ! Check the existence of the text file with statistics

    INQUIRE(FILE=trim(AuxObsDataFileCH),EXIST=LnewExists)
    if (.not.LnewExists) then
      WRITE(*,*) '-----------------------------------------------------------------------------'
      WRITE(*,*) 'WARNING! ose_readOmPstddev_auxfileCH: auxiliary file ' // trim(AuxObsDataFileCH)
      WRITE(*,*) 'WARNING! not available. Default CH family OmP stddev to be applied if needed.'
      WRITE(*,*) '-----------------------------------------------------------------------------'
      return
    end if

    ! Read error std dev. from file auxiliary file for constituent data

    nulstat=0
    ierr=fnom(nulstat,trim(AuxObsDataFileCH),'SEQ',0)
    if ( ierr == 0 ) then
      open(unit=nulstat, file=trim(AuxObsDataFileCH), status='OLD')
    else
      call utl_abort('ose_readOmPstddev_auxfileCH: COULD NOT OPEN AUXILIARY FILE ' //  trim(AuxObsDataFileCH) )
    end if
  
    ! Read OmP error standard deviations for constituents or related directives if available.
    
    ios=0
    read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne
    do while (trim(adjustl(ligne(1:12))) /= 'SECTION V:') 
        read(nulstat,'(A)',iostat=ios,err=10,end=15) ligne
    end do 
  
    ! Read number of observation set sub-families (STNIDs and ...) and allocate space

    read(nulstat,*,iostat=ios,err=10,end=10) OmPstdCH%n_stnid
    read(nulstat,*,iostat=ios,err=10,end=10) isize
  
    allocate(OmPstdCH%stnids(OmPstdCH%n_stnid),OmPstdCH%std_type(OmPstdCH%n_stnid))
    allocate(OmPstdCH%n_month(OmPstdCH%n_stnid),OmPstdCH%n_lat(OmPstdCH%n_stnid))
    allocate(OmPstdCH%source(OmPstdCH%n_stnid),OmPstdCH%ibegin(OmPstdCH%n_stnid))
    allocate(OmPstdCH%element(OmPstdCH%n_stnid),OmPstdCH%n_lvl(OmPstdCH%n_stnid))
    allocate(OmPstdCH%std(isize))
    allocate(OmPstdCH%levels(isize),OmPstdCH%lat(isize),OmPstdCH%month(isize))
 
    OmPstdCH%element(:)=0
    OmPstdCH%source(:)=0
    OmPstdCH%std_type(:)=0
    OmPstdCH%n_lvl(:)=1
    OmPstdCH%n_lat(:)=1
    OmPstdCH%n_month(:)=1

    ! Begin reading for each sub-family
    ! Important: Combination of STNID, BUFR element and number of vertical levels
    !            to determine association to the observations.

    icount=0
    STNIDLOOP: do stnidIndex=1,OmPstdCH%n_stnid

      ! disregard line of dashes
      read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne

      ! Read STNID (* as wildcard)    
      read(nulstat,'(2X,A9)',iostat=ios,err=10,end=10) OmPstdCH%stnids(stnidIndex) 

      !   Read (1) BUFR element
      !        (2) Flag indication if OmP error std dev provided from this auxiliary file or
      !            to be calculated via OmP differences or HPHT
      !        (3) Type of input/set std dev
      !        (4) Number of vertical levels
      !        (5) Number of latitudes
      !        (6) Number of months
      !
      !   Important: Combination of STNID and BUFR element
      !              to determine association to the observations.
    
      read(nulstat,*,iostat=ios,err=10,end=10) &
        OmPstdCH%element(stnidIndex),OmPstdCH%source(stnidIndex),  &
        OmPstdCH%std_type(stnidIndex),OmPstdCH%n_lvl(stnidIndex),  &
        OmPstdCH%n_lat(stnidIndex),OmPstdCH%n_month(stnidIndex)

      if ( OmPstdCH%n_lvl(stnidIndex) < 1 ) OmPstdCH%n_lvl(stnidIndex)=1
      if ( OmPstdCH%n_lat(stnidIndex) < 1 ) OmPstdCH%n_lat(stnidIndex)=1
      if ( OmPstdCH%n_month(stnidIndex) < 1 ) OmPstdCH%n_month(stnidIndex)=1

      if ( icount+OmPstdCH%n_lvl(stnidIndex)*OmPstdCH%n_lat(stnidIndex)*OmPstdCH%n_month(stnidIndex) > isize ) then
         write(*,'(10X,"Max array size exceeded: ",I6)') isize
         CALL utl_abort('ose_readOmPstddev_auxfileHCHin: PROBLEM READING OBSERR STD DEV.')  
      else if ( OmPstdCH%n_lat(stnidIndex) == 1 .and. OmPstdCH%n_month(stnidIndex) > 1 ) then
         write(*,'(10X,"Fails for stnid number: ",I6)') stnidIndex
         CALL utl_abort('ose_readOmPstddev_auxfileHCHin: Cannot depend on month if not dependent on latitude')  
      else if ( OmPstdCH%n_month(stnidIndex) /= 1 .and. OmPstdCH%n_month(stnidIndex) /= 12 ) then
         write(*,'(10X,"Fails for stnid number: ",I6)') stnidIndex
         CALL utl_abort('ose_readOmPstddev_auxfileHCHin: Number of months must be 1 or 12')  
      end if

      ! disregard line of dashes
      read(nulstat,'(A)',iostat=ios,err=10,end=10) ligne

      if ( OmPstdCH%source(stnidIndex) > 1 ) then
         ! Disregard data section
         OmPstdCH%ibegin(stnidIndex)=icount
         cycle STNIDLOOP 
      else if ( OmPstdCH%source(stnidIndex) == 1 ) then

        OmPstdCH%ibegin(stnidIndex)=icount+1
        
        ! Only reference vertical levels and latitudes to be determined 
        ! For use in determining error std dev from OmPs of current processing window

        OmPstdCH%n_month(stnidIndex) = 1          
        if ( OmPstdCH%n_lvl(stnidIndex) == 1 .and. OmPstdCH%n_lat(stnidIndex) == 1 ) then
           
          icount=icount+1
       
        else if ( OmPstdCH%n_lvl(stnidIndex) == 1 .and. OmPstdCH%n_lat(stnidIndex) > 1 ) then
           
          ! Read reference latitudes (must be in order of increasing size)
       
          read(nulstat,*,iostat=ios,err=10,end=10)                      &
               OmPstdCH%lat(icount+1:icount+OmPstdCH%n_lat(stnidIndex))
          icount=icount+OmPstdCH%n_lat(stnidIndex)
               
        else if ( OmPstdCH%n_lvl(stnidIndex) > 1 .and. OmPstdCH%n_lat(stnidIndex) == 1 ) then
        
          ! Read reference vertical levels (must be in order of increasing size)
       
          read(nulstat,*,iostat=ios,err=10,end=10)                      &
               OmPstdCH%levels(icount+1:icount+OmPstdCH%n_lvl(stnidIndex))
          icount=icount+OmPstdCH%n_lvl(stnidIndex)
                 
        else 
        
          ! Read reference latitudes (must be in order of increasing size)
          read(nulstat,*,iostat=ios,err=10,end=10)                      &
               OmPstdCH%lat(icount+1:icount+OmPstdCH%n_lat(stnidIndex))

          ! Read reference vertical levels (must be in order of increasing size)
          read(nulstat,*,iostat=ios,err=10,end=10)                      &
               OmPstdCH%levels(icount+1:icount+OmPstdCH%n_lvl(stnidIndex))
               
          icount=icount+OmPstdCH%n_lvl(stnidIndex)*OmPstdCH%n_lat(stnidIndex)

        end if
        
        cycle STNIDLOOP
         
      end if
      
      ! For OmPstdCH%source(stnidIndex) == 0
      
      OmPstdCH%ibegin(stnidIndex)=icount+1
      if ( OmPstdCH%n_lvl(stnidIndex) == 1 .and. OmPstdCH%n_lat(stnidIndex) == 1 ) then
    
        ! Read one value only (independent of level and latitude)
        ! Assumes one month only.
        
        read(nulstat,*,iostat=ios,err=10,end=10) OmPstdCH%std(icount+1)
        icount=icount+1

      else if ( OmPstdCH%n_lvl(stnidIndex) == 1 .and. OmPstdCH%n_lat(stnidIndex) > 1 ) then
    
        ! Value dependent on latitude 
       
        ! Read reference latitudes (must be in order of increasing size)
       
        read(nulstat,*,iostat=ios,err=10,end=10)                      &
               OmPstdCH%lat(icount+1:icount+OmPstdCH%n_lat(stnidIndex))
      
        ! Read OMP error std dev related values
  
        if ( OmPstdCH%n_month(stnidIndex) == 1 ) then
        
          read(nulstat,*,iostat=ios,err=10,end=10)                 &
                   OmPstdCH%std(icount+1:icount+OmPstdCH%n_lat(stnidIndex))
        
        else
        
          do monthIndex=1,OmPstdCH%n_month(stnidIndex)
          
            ! Read reference month (must be in order of increasing size)     
            read(nulstat,*,iostat=ios,err=10,end=10) OmPstdCH%month(icount+monthIndex)

            ipos=icount+(monthIndex-1)*OmPstdCH%n_lat(stnidIndex)
                       
            read(nulstat,*,iostat=ios,err=10,end=10)                 &
                   OmPstdCH%std(ipos+1:ipos+OmPstdCH%n_lat(stnidIndex))
          end do
                               
        end if
        
        icount = icount + OmPstdCH%n_lat(stnidIndex)*OmPstdCH%n_month(stnidIndex)
        
      else if ( OmPstdCH%n_lvl(stnidIndex) > 1 .and. OmPstdCH%n_lat(stnidIndex) == 1 ) then
    
       ! Value dependent on vertical level and not latitude
       
        if ( OmPstdCH%n_month(stnidIndex) == 1 ) then
        
          do levIndex=1,OmPstdCH%n_lvl(stnidIndex)
            icount=icount+1
            
            ! Read vertical level and OMP error std dev related value.
          
            read(nulstat,*,iostat=ios,err=10,end=10)                 &
                 OmPstdCH%levels(icount),OmPstdCH%std(icount)

          end do
          
        else
        
          do monthIndex=1,OmPstdCH%n_month(stnidIndex)
          
            ! Read reference month (must be in order of increasing size)     
            read(nulstat,*,iostat=ios,err=10,end=10) OmPstdCH%month(icount+monthIndex)

            ipos=icount+(monthIndex-1)*OmPstdCH%n_lvl(stnidIndex)
                       
            read(nulstat,*,iostat=ios,err=10,end=10)                 &
                   OmPstdCH%std(ipos+1:ipos+OmPstdCH%n_lvl(stnidIndex))
          end do
          
        end if
         
        icount = icount + OmPstdCH%n_lvl(stnidIndex)*OmPstdCH%n_month(stnidIndex)
   
      else if ( OmPstdCH%n_lvl(stnidIndex) > 1 .and. OmPstdCH%n_lat(stnidIndex) > 1 ) then
    
        ! Value dependent on vertical level and latitude 
       
        ! Read reference latitudes (must be in order of increasing size)
        read(nulstat,*,iostat=ios,err=10,end=10)                      &
               OmPstdCH%lat(icount+1:icount+OmPstdCH%n_lat(stnidIndex))

        if ( OmPstdCH%n_month(stnidIndex) == 1 ) then
            
          do levIndex=1,OmPstdCH%n_lvl(stnidIndex)
          
            ! Read vertical level and OMP error std dev related lat-dependent values.
          
            read(nulstat,*,iostat=ios,err=10,end=10)    &
                 OmPstdCH%levels(icount+levIndex),    &
                 OmPstdCH%std(icount+(levIndex-1)*    &
                 OmPstdCH%n_lat(stnidIndex)+1:icount+levIndex*OmPstdCH%n_lat(stnidIndex))

          end do

        else
        
          do monthIndex=1,OmPstdCH%n_month(stnidIndex)
          
            ! Read reference month (must be in order of increasing size)     
            read(nulstat,*,iostat=ios,err=10,end=10) OmPstdCH%month(icount+monthIndex)

            ipos=icount+(monthIndex-1)*OmPstdCH%n_lat(stnidIndex)*OmPstdCH%n_lvl(stnidIndex)
                             
            do levIndex=1,OmPstdCH%n_lvl(stnidIndex)
          
              ! Read vertical level and OMP error std dev related lat-dependent values.
          
              read(nulstat,*,iostat=ios,err=10,end=10)     &
                 OmPstdCH%levels(icount+levIndex),       &
                 OmPstdCH%std(ipos+(levIndex-1)* &
                 OmPstdCH%n_lat(stnidIndex)+1:ipos+levIndex*OmPstdCH%n_lat(stnidIndex))

            end do

          end do
                
        end if
        
        icount=icount+OmPstdCH%n_lat(stnidIndex)*OmPstdCH%n_lvl(stnidIndex)*OmPstdCH%n_month(stnidIndex)
      end if
    end do STNIDLOOP
     
 10 if (ios.gt.0) then
      write(abortText,*) ios
      call utl_abort('ose_readOmPstddev_auxfileCHin: PROBLEM READING OMP ERR STD DEV. - ' // &
                     'File read error message number: ' // trim(abortText) )    
    end if
    
    return
    
    ! Reached end of file and no related section found.
    
 15 OmPstdCH%n_stnid = 0

    close(unit=nulstat)
    ierr=fclos(nulstat)    

  end subroutine ose_readOmPstddev_auxfileCH
  
  !--------------------------------------------------------------------------
  ! ose_calcOmPstddevCH
  !--------------------------------------------------------------------------
  subroutine ose_calcOmPstddevCH(obsSpaceData)
    ! 
    !:Purpose: To calc OmP error std dev for some obs sets of the CH family
    !
    implicit none
    
    ! Arguments:
    type(struct_obs) :: obsSpaceData ! observation-space data; output saved in OBS_OMPE column
    
    ! Local:
    logical :: availableOmP
    integer :: stnidIndex, headerIndex, bodyIndex, bodyIndex_start, bodyIndex_end, icodtyp, ierr
    integer :: idate, itime, iass, latIndex, levIndex, monthIndex, ibegin, loopIndex, posIndex
    real(8) :: zlat, zval, zlev, lat, sumOmP, sumSqrOmP, varOmP, maxOmP,meanOmP,medianOmP
    character(len=12) :: stnid
    real(8), allocatable :: series(:,:,:),sumOmP2d(:,:),sumSqrOmP2d(:,:)
    real(8), allocatable :: sumOmP2dt(:,:),sumSqrOmP2dt(:,:)
    integer, allocatable :: nSeries(:,:),nSeriest(:,:)
    integer, parameter :: maxCount = 10000
    integer, parameter :: minCount = 5
    real(8), parameter :: stdScale = 2.0
    real(8), parameter :: minVal = 1.0d-20
    real(4) :: rseries(maxCount)
    integer :: ip(maxCount)
    
    ! Loop over all obs types
 
    do stnidIndex=1,OmPstdCH%n_stnid
      ! Check if OBS_OMPE is to be estimated elsewhere via use of OBS_HPHT and OBS_OER
      ! instead of via ose_fillOmPstddevCH
      if ( OmPstdCH%source(stnidIndex) /= 1 ) cycle

      write(*,*) 
      write(*,*) 'ose_calcOmPstddevCH: Online calculation of OmP error std dev of CH family stnid ',OmPstdCH%stnids(stnidIndex)
      write(*,*)

      ! Identify start position index (-1) in OmPstdCH  
      ibegin=OmPstdCH%ibegin(stnidIndex)-1
      
      allocate(series(OmPstdCH%n_lvl(stnidIndex),OmPstdCH%n_lat(stnidIndex),maxCount))
      allocate(nSeries(OmPstdCH%n_lvl(stnidIndex),OmPstdCH%n_lat(stnidIndex)))        
      allocate(sumOmP2d(OmPstdCH%n_lvl(stnidIndex),OmPstdCH%n_lat(stnidIndex)))
      allocate(sumSqrOmP2d(OmPstdCH%n_lvl(stnidIndex),OmPstdCH%n_lat(stnidIndex)))
      nSeries(:,:)=0
      sumOmP2d(:,:)=0.0d0
      sumSqrOmP2d(:,:)=0.0d0
      
      allocate(nSeriest(OmPstdCH%n_lvl(stnidIndex),OmPstdCH%n_lat(stnidIndex)))        
      allocate(sumOmP2dt(OmPstdCH%n_lvl(stnidIndex),OmPstdCH%n_lat(stnidIndex)))
      allocate(sumSqrOmP2dt(OmPstdCH%n_lvl(stnidIndex),OmPstdCH%n_lat(stnidIndex)))
      nSeriest(:,:)=0
      
      ! Loop over all header indices of the 'CH' family:
      
      call obs_set_current_header_list(obsSpaceData,'CH')
      HEADER: do

        headerIndex = obs_getHeaderIndex(obsSpaceData)
        if (headerIndex < 0) exit HEADER
  
        icodtyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)        
        if (icodtyp.ne.codtyp_get_codtyp('CHEMREMOTE').and.icodtyp.ne.codtyp_get_codtyp('CHEMINSITU')) cycle HEADER
      
        stnid = obs_elem_c(obsSpaceData,'STID',headerIndex)
        if ( .not. utl_stnid_equal(OmPstdCH%stnids(stnidIndex),stnid) ) cycle HEADER

        bodyIndex_start = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
        bodyIndex_end = bodyIndex_start+obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex)-1

        ! Check OBS_VNM value.
        do bodyIndex=bodyIndex_start,bodyIndex_end
          if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) /= BUFR_SCALE_EXPONENT) then
             if ( obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) /= OmPstdCH%element(stnidIndex) ) cycle HEADER
          end if
        end do
           
        zlat   = obs_headElem_r( obsSpaceData, OBS_LAT, headerIndex )
        idate  = obs_headElem_i( obsSpaceData, OBS_DAT, headerIndex ) 
        itime  = obs_headElem_i( obsSpaceData, OBS_ETM, headerIndex )

        ! Identify lat and month index
        if (OmPstdCH%n_lat(stnidIndex) == 1) then      
           latIndex = 1
        else

          ! Find latitude index for specifying the latitude bin
          ! Assuming increasing latitudes in OmPstdCH%lat

          lat = zlat / MPC_RADIANS_PER_DEGREE_R8  ! radians to degrees

          if (lat >= 0.5*(OmPstdCH%lat(ibegin+OmPstdCH%n_lat(stnidIndex)) &
                          +OmPstdCH%lat(ibegin+OmPstdCH%n_lat(stnidIndex)-1)) ) then
            latIndex = OmPstdCH%n_lat(stnidIndex)
          else
            do latIndex=1,OmPstdCH%n_lat(stnidIndex)-1
              if (lat <= 0.5*(OmPstdCH%lat(ibegin+latIndex)+OmPstdCH%lat(ibegin+latIndex+1)) ) exit
            end do
          end if
        end if
        
        if (OmPstdCH%n_month(stnidIndex) == 1) then
          monthIndex = 1 
        else
          ! Find month index
          monthIndex=(idate-(idate/10000)*10000)/100
        end if

        do bodyIndex=bodyIndex_start,bodyIndex_end
          if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) == BUFR_SCALE_EXPONENT) cycle

          iass  = obs_bodyElem_i( obsSpaceData, OBS_ASS, bodyIndex )
          if ( iass /= obs_assimilated ) cycle
           
          ! Identify vertical level 
          zlev  = obs_bodyElem_r( obsSpaceData, OBS_PPP, bodyIndex )
          if ( OmPstdCH%n_lvl(stnidIndex) == 1 ) then
            levIndex = 1
          else if ( OmPstdCH%levels(ibegin+1) < OmPstdCH%levels(ibegin+2) ) then             
            if (zlev <= 0.5*(OmPstdCH%levels(ibegin+1)+OmPstdCH%levels(ibegin+2)) ) then
              levIndex = 1
            else
              do levIndex=OmPstdCH%n_lvl(stnidIndex),2,-1
                if ( zlev >= 0.5*(OmPstdCH%levels(ibegin+levIndex)+OmPstdCH%levels(ibegin+levIndex-1)) ) exit
              end do
            end if
          else
            if (zlev <= 0.5*(OmPstdCH%levels(ibegin+OmPstdCH%n_lvl(stnidIndex)) &
                            +OmPstdCH%levels(ibegin+OmPstdCH%n_lvl(stnidIndex)-1)) ) then
              levIndex = OmPstdCH%n_lvl(stnidIndex)
            else
              do levIndex=1,OmPstdCH%n_lvl(stnidIndex)-1
                if ( zlev >= 0.5*(OmPstdCH%levels(ibegin+levIndex)+OmPstdCH%levels(ibegin+levIndex+1)) ) exit
              end do
            end if
          end if
          
          ! Store OmP in work array
          
          if ( OmPstdCH%std_type(stnidIndex) == 1 ) then
            zval  = obs_bodyElem_r( obsSpaceData, OBS_VAR, bodyIndex )
            if ( zval > minVal ) then
              nSeries(levIndex,latIndex) = nSeries(levIndex,latIndex) + 1
              series(levIndex,latIndex,nSeries(levIndex,latIndex)) = obs_bodyElem_r( obsSpaceData, OBS_OMP, bodyIndex )/zval
            end if
          else
            nSeries(levIndex,latIndex) = nSeries(levIndex,latIndex) + 1
            series(levIndex,latIndex,nSeries(levIndex,latIndex)) =  obs_bodyElem_r( obsSpaceData, OBS_OMP, bodyIndex )
          end if
          ! Check if enough data
          if ( nSeries(levIndex,latIndex) == maxCount ) exit
          
        end do
        
      end do HEADER
     
      if ( any(nSeries > 0) ) then
      
        ! Calc OmP error std dev
      
        !write(*,*) 'latbin levbin Npts      AvgOmP    OmPStddev'
        do levIndex=1,OmPstdCH%n_lvl(stnidIndex)
        do latIndex=1,OmPstdCH%n_lat(stnidIndex)
          availableOmP=.true.
          if (nSeries(levIndex,latIndex) > minCount ) then
            rseries(1:nSeries(levIndex,latIndex))=series(levIndex,latIndex,1:nSeries(levIndex,latIndex))
            call ipsort(ip,rseries(1:nSeries(levIndex,latIndex)),nSeries(levIndex,latIndex))
            medianOmP=series(levIndex,latIndex,ip(nSeries(levIndex,latIndex)/2))
               
            sumOmP = sum(series(levIndex,latIndex,1:nSeries(levIndex,latIndex)))
            meanOmP=sumOmP/nSeries(levIndex,latIndex)            
            sumSqrOmP = sum(series(levIndex,latIndex,1:nSeries(levIndex,latIndex))*series(levIndex,latIndex,1:nSeries(levIndex,latIndex)))
            varOmP = sumSqrOmP/nSeries(levIndex,latIndex) - (sumOmP/nSeries(levIndex,latIndex))**2
            if (varOmP  > minVal*minVal ) then
              do loopIndex=1,nSeries(levIndex,latIndex)
                maxOmP = stdScale*sqrt(varOmP) 
                ! if ( abs(series(levIndex,latIndex,loopIndex)-meanOmP) > maxOmP ) then
                if ( abs(series(levIndex,latIndex,loopIndex)-medianOmP) > maxOmP ) then
                  sumOmP = sumOmP - series(levIndex,latIndex,loopIndex)
                  sumSqrOmP = sumSqrOmP - series(levIndex,latIndex,loopIndex)*series(levIndex,latIndex,loopIndex)
                  nSeries(levIndex,latIndex) = nSeries(levIndex,latIndex) - 1
                end if
              end do
              varOmP = sumSqrOmP/nSeries(levIndex,latIndex) - (sumOmP/nSeries(levIndex,latIndex))**2
              !write(*,110) latIndex,levIndex,nSeries(levIndex,latIndex),sumOmP/nSeries(levIndex,latIndex),sqrt(varOmP)
              sumOmP2d(levIndex,latIndex)=sumOmP
              sumSqrOmP2d(levIndex,latIndex)=sumSqrOmP   
              if (varOmP  < minVal*minVal ) availableOmP = .false.
            else
              availableOmP = .false.
            end if
          else
            availableOmP = .false.
          end if
          if (.not.availableOmP) nSeries(levIndex,latIndex) = 0
        end do
        end do
      end if
                  
      ! Combine from all processors

      call rpn_comm_allreduce(nSeries,nSeriest,OmPstdCH%n_lvl(stnidIndex)*OmPstdCH%n_lat(stnidIndex),"MPI_INTEGER","MPI_SUM","GRID",ierr)
      call rpn_comm_allreduce(sumOmP2d,sumOmP2dt,OmPstdCH%n_lvl(stnidIndex)*OmPstdCH%n_lat(stnidIndex),"MPI_DOUBLE_PRECISION","MPI_SUM","GRID",ierr)
      call rpn_comm_allreduce(sumSqrOmP2d,sumSqrOmP2dt,OmPstdCH%n_lvl(stnidIndex)*OmPstdCH%n_lat(stnidIndex),"MPI_DOUBLE_PRECISION","MPI_SUM","GRID",ierr)
      
      if (any(nSeriest > 5*minCount)) then
      
        where ( nSeriest > 5*minCount ) sumSqrOmP2dt = sumSqrOmP2dt/nSeriest - (sumOmP2dt/nSeriest)**2
        
        if (mmpi_myid == 0) then
           write(*,*) 'Calculated OmP error std dev'
           write(*,*) 'latbin levbin Npts      AvgOmP    OmPStddev'
        end if
        do levIndex=1,OmPstdCH%n_lvl(stnidIndex)
        do latIndex=1,OmPstdCH%n_lat(stnidIndex)
          availableOmP=.true.
          if (nSeriest(levIndex,latIndex) > 5*minCount ) then
            varOmP =  sumSqrOmP2dt(levIndex,latIndex)
            if (varOmP  < minVal*minVal ) then
               availableOmP = .false.
            else     
              if (mmpi_myid == 0 ) write(*,110) latIndex,levIndex,nSeriest(levIndex,latIndex),sumOmP2dt(levIndex,latIndex)/nSeriest(levIndex,latIndex),sqrt(varOmP)
            end if
          else
            availableOmP = .false.
          end if
  
          ! Store in structure
   
          ! Indentify position in structure
          if (OmPstdCH%n_lvl(stnidIndex) > 1) then
            posIndex=ibegin+(levIndex-1)*OmPstdCH%n_lat(stnidIndex)
          else
            posIndex=ibegin
          end if
      
          if ( OmPstdCH%n_month(stnidIndex) == 1 ) then
            posIndex = posIndex + latIndex
          else
            posIndex = posIndex + (monthIndex-1)*OmPstdCH%n_lvl(stnidIndex)*OmPstdCH%n_lat(stnidIndex) + latIndex
          end if
        
          if (availableOmP) then
            OmPstdCH%std(posIndex)  = sqrt(varOmP)
          else
             OmPstdCH%std(posIndex)  = MPC_missingValue_R8
          end if
        
        end do
        end do
      end if  
110   format(I5,I7,I7,3x,G11.3,G11.3)
                 
      deallocate(series,nSeries,sumOmP2d,sumSqrOmP2d)
      deallocate(nSeriest,sumOmP2dt,sumSqrOmP2dt)
    end do
    
  end subroutine ose_calcOmPstddevCH

  !--------------------------------------------------------------------------
  ! ose_fillOmPstddevCH
  !--------------------------------------------------------------------------
  subroutine ose_fillOmPstddevCH(obsSpaceData)
    ! 
    !:Purpose: To assign the Omp error std dev where possible for the obs 
    !          of the CH obs family.
    !
    implicit none

    ! Arguments:
    type(struct_obs) :: obsSpaceData ! observation-space data; output saved in OBS_OMPE column
    
    ! Local:
    integer :: stnidIndex, headerIndex, bodyIndex, bodyIndex_start, bodyIndex_end, icodtyp
    integer :: idate, itime, iass, latIndex, monthIndex, ibegin
    real(8) :: zlat, zlev, OmP_err_stddev, lat
    character(len=12) :: stnid

    ! Loop over all obs types
    
    do stnidIndex=1,OmPstdCH%n_stnid
      ! Check if OBS_OMPE is to be estimated elsewhere via use of OBS_HPHT and OBS_OER
      ! instead of via ose_fillOmPstddevCH
      if ( OmPstdCH%source(stnidIndex) > 1 ) cycle

      ! Loop over all header indices of the 'CH' family:
       
      call obs_set_current_header_list(obsSpaceData,'CH')
      HEADER: do
      
        headerIndex = obs_getHeaderIndex(obsSpaceData)
        if (headerIndex < 0) exit HEADER
  
        icodtyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
        if (icodtyp.ne.codtyp_get_codtyp('CHEMREMOTE').and.icodtyp.ne.codtyp_get_codtyp('CHEMINSITU')) cycle HEADER
      
        stnid = obs_elem_c(obsSpaceData,'STID',headerIndex)
        if ( .not. utl_stnid_equal(OmPstdCH%stnids(stnidIndex),stnid) ) cycle HEADER
        
        bodyIndex_start = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
        bodyIndex_end = bodyIndex_start+obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex)-1
        
        ! Check OBS_VNM value.
        do bodyIndex=bodyIndex_start,bodyIndex_end
          if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) /= BUFR_SCALE_EXPONENT) then
            if ( obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) /= OmPstdCH%element(stnidIndex) ) cycle HEADER
          end if
        end do

        zlat   = obs_headElem_r( obsSpaceData, OBS_LAT, headerIndex )
        idate  = obs_headElem_i( obsSpaceData, OBS_DAT, headerIndex )
        itime  = obs_headElem_i( obsSpaceData, OBS_ETM, headerIndex )

        ! Identify lat and month index
        if (OmPstdCH%n_lat(stnidIndex) == 1) then      
           latIndex = 0 ! no interpolation
        else

          ! Find latitude index for interpolation.
          ! Assuming increasing latitudes in OmPstdCH%lat
          ! Interpolation to between latIndex and latIndex-1
          
          lat = zlat / MPC_RADIANS_PER_DEGREE_R8  ! radians to degrees

          ibegin=OmPstdCH%ibegin(stnidIndex)-1
          if ( lat <= OmPstdCH%lat(ibegin+1) ) then
            latIndex=2
          else if ( lat >= OmPstdCH%lat(ibegin+OmPstdCH%n_lat(stnidIndex)) ) then
            latIndex=OmPstdCH%n_lat(stnidIndex)
          else
            do latIndex=2,OmPstdCH%n_lat(stnidIndex)
              if (lat >= OmPstdCH%lat(ibegin+latIndex-1) .and. &
                  lat <= OmPstdCH%lat(ibegin+latIndex)) exit
            end do
          end if
        end if
        
        if (OmPstdCH%n_month(stnidIndex) == 1) then
          monthIndex = 1 
        else
          ! Find month index
          monthIndex=(idate-(idate/10000)*10000)/100
        end if

        do bodyIndex=bodyIndex_start,bodyIndex_end
          if (obs_bodyElem_i(obsSpaceData,OBS_VNM,bodyIndex) == BUFR_SCALE_EXPONENT) cycle

          zlev  = obs_bodyElem_r( obsSpaceData, OBS_PPP, bodyIndex )
                    
          ! Get OmP error std dev 
          OmP_err_stddev = ose_getOmPstddevCH( lat, zlev, stnidIndex, latIndex, monthIndex )
          if ( OmPstdCH%std_type(stnidIndex) == 1 ) then
            iass  = obs_bodyElem_i( obsSpaceData, OBS_ASS, bodyIndex )
            if ( iass == obs_assimilated ) then
              OmP_err_stddev = obs_bodyElem_r( obsSpaceData, OBS_VAR, bodyIndex )*OmP_err_stddev
              call obs_bodySet_r( obsSpaceData, OBS_OMPE, bodyIndex, OmP_err_stddev )
            end if
          else
             call obs_bodySet_r( obsSpaceData, OBS_OMPE, bodyIndex, OmP_err_stddev )             
          end if

        end do
      end do HEADER
    end do
   
  end subroutine ose_fillOmPstddevCH
  
  !--------------------------------------------------------------------------
  ! ose_getOmPstddevCH
  !--------------------------------------------------------------------------
  real(8) function ose_getOmPstddevCH(zlat,zlev,stnidIndex,latIndex,monthIndex) result(OmP_err_stddev) 
    ! 
    !:Purpose: To return the OmP error std dev for a CH family measurement
    !
    implicit none

    ! Arguments:
    real(8), intent(in) :: zlat  ! latitude (radians)
    real(8), intent(in) :: zlev  ! vertical coordinate value
    integer, intent(in) :: stnidIndex ! station and obs type index
    integer, intent(in) :: latIndex   ! reference lat for interpolation
    integer, intent(in) :: monthIndex ! month index

    ! Locals:
    real(8) :: levDiff
    integer :: ibegin,levIndex,loopIndex,posIndex
    real(8), parameter :: minDiff = 0.01*abs(MPC_missingValue_R8)

    OmP_err_stddev = MPC_missingValue_R8

    ! Get OmP error std. dev.

    ibegin=OmPstdCH%ibegin(stnidIndex)-1
    
    if (OmPstdCH%n_lvl(stnidIndex) > 1) then
                 
      ! Find nearest vertical level (no vertical interpolation in this version)
                 
      levDiff=1.0d10
      do loopIndex=1,OmPstdCH%n_lvl(stnidIndex)
         if ( levDiff > abs(zlev-OmPstdCH%levels(ibegin+loopIndex)) ) THEN
            levIndex=loopIndex
            levDiff=abs(zlev-OmPstdCH%levels(ibegin+loopIndex))
         END IF
      END DO
      posIndex=ibegin+(levIndex-1)*OmPstdCH%n_lat(stnidIndex)
    else
      posIndex=ibegin
    end if
      
    if ( OmPstdCH%n_month(stnidIndex) == 1 ) then
       posIndex = posIndex + latIndex
    else
       posIndex = posIndex + (monthIndex-1)*OmPstdCH%n_lat(stnidIndex)*OmPstdCH%n_lvl(stnidIndex) + latIndex
    end if   
    
    if ( OmPstdCH%n_lat(stnidIndex) > 1 ) then
      ! Apply interpolation in latitude
  
      if ( latIndex == 1 .or. latIndex > OmPstdCH%n_lat(stnidIndex)) then
        if ( abs(OmPstdCH%std(posIndex) - MPC_missingValue_R8) > minDiff ) then
          OmP_err_stddev = MPC_missingValue_R8    
        else     
          OmP_err_stddev = OmPstdCH%std(posIndex)
        end if
      else        
        if ( abs(OmPstdCH%std(posIndex-1) - MPC_missingValue_R8) < minDiff .and. &
             abs(OmPstdCH%std(posIndex) - MPC_missingValue_R8) < minDiff ) then
          OmP_err_stddev = MPC_missingValue_R8  
        else if ( abs(OmPstdCH%std(posIndex-1) - MPC_missingValue_R8) < minDiff ) then
          OmP_err_stddev = OmPstdCH%std(posIndex)
        else if ( abs(OmPstdCH%std(posIndex) - MPC_missingValue_R8) < minDiff ) then
          OmP_err_stddev = OmPstdCH%std(posIndex-1)
        else
          OmP_err_stddev = (OmPstdCH%std(posIndex-1)*(OmPstdCH%lat(ibegin+latIndex)-zlat)+ &
            OmPstdCH%std(posIndex)*(zlat-OmPstdCH%lat(ibegin+latIndex-1)))/                 &
            (OmPstdCH%lat(ibegin+latIndex)-OmPstdCH%lat(ibegin+latIndex-1))
        end if
      end if       
    else
      if ( abs(OmPstdCH%std(posIndex) - MPC_missingValue_R8) > minDiff ) then
        OmP_err_stddev = MPC_missingValue_R8    
      else     
        OmP_err_stddev = OmPstdCH%std(posIndex)  
      end if           
    end if 

  end function ose_getOmPstddevCH

  !--------------------------------------------------------------------------
  ! ose_OmPstddevExistsForAllCH
  !--------------------------------------------------------------------------
  logical function ose_OmPstddevExistsForAllCH(obsSpaceData) result(allOMPE)
    ! 
    !:Purpose: To determine if all obs to be processed have usable OBS_OMPE values 
    !          for the CH obs family.
    !
    implicit none

    ! Arguments:
    type(struct_obs) :: obsSpaceData ! observation-space data; output saved in OBS_OMPE column
    
    ! Local:
    integer :: headerIndex, bodyIndex, bodyIndex_start, bodyIndex_end, icodtyp
     
    ! Loop over all header indices of the 'CH' family:
       
    allOMPE = .true.
    
    call obs_set_current_header_list(obsSpaceData,'CH')
    HEADER: do

      headerIndex = obs_getHeaderIndex(obsSpaceData)
      if (headerIndex < 0) exit HEADER
  
      icodtyp = obs_headElem_i(obsSpaceData,OBS_ITY,headerIndex)
      if (icodtyp.ne.codtyp_get_codtyp('CHEMREMOTE').and.icodtyp.ne.codtyp_get_codtyp('CHEMINSITU')) cycle HEADER
      
      bodyIndex_start = obs_headElem_i(obsSpaceData,OBS_RLN,headerIndex)
      bodyIndex_end = bodyIndex_start+obs_headElem_i(obsSpaceData,OBS_NLV,headerIndex)-1

      ! Check for cases were OmP error std dev is not available
      
      do bodyIndex=bodyIndex_start,bodyIndex_end
         if ( obs_bodyElem_r(obsSpaceData,OBS_OMPE,bodyIndex) <= 0.0d0 .and. &
              obs_bodyElem_i(obsSpaceData,OBS_ASS,bodyIndex) == obs_assimilated ) then
            allOMPE = .false.
            write(*,*)
            write(*,*) 'ose_OmPstddevExistsForAllCH: Not all CH obs to be processed have available OmP error std dev.'
            write(*,*) 'HPHT will be calculated for the obs that do not have OmP error std dev available.'
            write(*,*)
            return
         end if
      end do
    end do HEADER
   
  end function ose_OmPstddevExistsForAllCH

  !--------------------------------------------------------------------------
  ! ose_deallocOmPstddevCH
  !--------------------------------------------------------------------------
  subroutine ose_deallocOmPstddevCH
    ! 
    !:Purpose: To deallocate temporary storage space used for OmP error std dev
    !          for the CH family.
    !
    implicit none

    if (OmPstdCH%n_stnid.eq.0) return

    if (allocated(OmPstdCH%stnids))   deallocate(OmPstdCH%stnids)
    if (allocated(OmPstdCH%n_lvl))    deallocate(OmPstdCH%n_lvl)
    if (allocated(OmPstdCH%ibegin))   deallocate(OmPstdCH%ibegin)
    if (allocated(OmPstdCH%element))  deallocate(OmPstdCH%element)
    if (allocated(OmPstdCH%source))   deallocate(OmPstdCH%source)
    if (allocated(OmPstdCH%std_type)) deallocate(OmPstdCH%std_type)
    if (allocated(OmPstdCH%n_lat))    deallocate(OmPstdCH%n_lat)
    if (allocated(OmPstdCH%std))      deallocate(OmPstdCH%std)
    if (allocated(OmPstdCH%levels))   deallocate(OmPstdCH%levels)
    if (allocated(OmPstdCH%lat))      deallocate(OmPstdCH%lat)
    if (allocated(OmPstdCH%n_month))  deallocate(OmPstdCH%n_month)
    if (allocated(OmPstdCH%month))    deallocate(OmPstdCH%month)

  end subroutine ose_deallocOmPstddevCH

end module obsSpaceErrorStdDev_mod
