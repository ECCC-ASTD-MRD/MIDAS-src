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
!! MODULE multi_ir_bgcheck (prefix='irbg' category='1. High-level functionality')
!!
!! *Purpose*: VARIABLES FOR MULTISPECTRAL INFRARED BACKGROUND CHECK
!!            AND QUALITY CONTROL.
!!
!! @author A. BEAULNE (CMDA/SMC) February 2006
!!
!! REVISION: adapted to IASI and CrIS by S. Heilliette
!!
!--------------------------------------------------------------------------
module multi_ir_bgck_mod
  use rttov_interfaces_mod
  use tovs_nl_mod
  use rttov_const, only : inst_id_iasi
  use utilities_mod
  use obsSpaceData_mod
  use mpi_mod
  use columnData_mod
  use mpivar_mod
  use obsFiles_mod
  use burp_module
  use EarthConstants_mod
  use MathPhysConstants_mod
  use verticalCoord_mod
  use presProfileOperators_mod

  implicit none
  save
  private
! Public functions (methods)
  public :: irbg_setup, irbg_bgCheckIR

  integer ,parameter :: nClassAVHRR = 7
  integer ,parameter :: NIR = 3, NVIS = 3
  integer ,parameter :: nChanAVHRR = NIR + NVIS
  integer ,parameter :: NMAXINST = 10

  character(len=6) :: INST(NMAXINST)
  integer :: NINST
  ! Reference (and alternate) window channel for clear / cloudy profile detection
  ! (subroutine cloud_height)

  integer :: IWINDOW(NMAXINST), IWINDOW_ALT(NMAXINST)
 
  ! Number of channels (and their values) to use for cloud top height detection
  ! with the "background profile matching" method (subroutine cloud_top)

  integer, parameter        :: NCH_HE = 4

  integer :: ILIST1(NMAXINST,NCH_HE) 

  ! Number of channels (and their values) to use for cloud top height detection
  ! with the CO2-slicing method. IREFR is the reference channel number (and alternate).
  ! (subroutine co2_slicing)


  integer, parameter  :: NCO2 = 13

  integer :: ILIST2(NMAXINST,NCO2), ILIST2_PAIR(NMAXINST,NCO2)

  ! Cloud top units : (1) mb, (2) meters
  ! (subroutines cloud_height (IOPT1) and cloud_top (IOPT2))

  integer, parameter        :: IOPT1 = 2   ! verify subr input if iopt1 changes
  integer, parameter        :: IOPT2 = 1

  ! Cloud top based on which background profile matching (subroutine cloud_top)
  ! (0) brightness temperature, (1) radiance, (2) both

  integer, parameter        :: IHGT = 2

  ! Maximum delta temperature allowed between guess and true skin temperature
  ! over water (DTW) and land (DTL)   (subroutine airsqc)

  real(8) :: DTW,DTL 

  ! Minimum and maximum RTTOV levels for LEV_START variable entering CO2 slicing
  ! In mb, between 50mb and 325mb (subroutine co2_slicing)

  real(8) :: PCO2MIN, PCO2MAX

  ! First channel affected by sun (for channels used only at night)
  ! (subroutine airsqc)

  integer :: ICHN_SUN(NMAXINST)
  
  ! Minimum solar zenith angle for night (between 90 and 180)
  ! (subroutine airsqc)

  real(8) :: NIGHT_ANG

  ! Highest flag in post files (value of N in 2^N)
  ! Currently 21

  integer, parameter :: BITFLAG = 29

  real(8),parameter :: seuilalb_static(NIR,0:2)= reshape( (/ 70.0,67.0,50.0, &
                                                             40.0,37.0,37.0, &
                                                             70.0,57.0,40. /),(/3,3/) ) 
  real(8),parameter :: seuilalb_homog(NIR,0:2)= reshape( (/ 15.0,18.0,13.0, &
                                                            9.0,10.0,10.0, &
                                                            18.0,16.0,10.0 /),(/3,3/) )
  
  real(8) :: seuilbt_homog(NVIS+1:NVIS+NIR,0:2,1:2)= reshape( (/5.d0, 4.d0, 4.d0, 4.d0, 3.d0, 3.d0, &
                                                                5.d0, 4.d0, 4.d0, 5.d0, 5.d0, 5.d0, &
                                                                4.d0, 3.d0, 3.d0, 5.d0, 5.d0, 5.d0/), (/3,3,2/) )

  namelist /NAMBGCKIR/ NINST, INST, IWINDOW, IWINDOW_ALT, ILIST1, ILIST2, ILIST2_PAIR, ICHN_SUN
  namelist /NAMBGCKIR/ DTW, DTL, PCO2MIN, PCO2MAX, NIGHT_ANG

  type( rttov_coefs ) :: coefs_avhrr

  type avhrr_bgck_iasi
     real(8)              :: RADMOY(nClassAVHRR,nChanAVHRR)
     real(8)              :: RADSTD(nClassAVHRR,nChanAVHRR)
     real(8)              :: CFRAC(nClassAVHRR)
     real(8)              :: TBMOY(nClassAVHRR,NVIS+1:NVIS+NIR)
     real(8)              :: TBSTD(nClassAVHRR,NVIS+1:NVIS+NIR)
     real(8)              :: ALBEDMOY(nClassAVHRR,1:NVIS)
     real(8)              :: ALBEDSTD(nClassAVHRR,1:NVIS)
     real(8)              :: TBSTD_PIXELIASI(NVIS+1:NVIS+NIR)
     real(8)              :: ALBSTD_PIXELIASI(1:NVIS)
     real(8)              :: RADCLEARCALC(NVIS+1:NVIS+NIR)
     real(8)              :: TBCLEARCALC(NVIS+1:NVIS+NIR)
     real(8)              :: RADOVCALC( tvs_nlevels-1,NVIS+1:NVIS+NIR)
     real(8)              :: TRANSMCALC( tvs_nlevels,NVIS+1:NVIS+NIR)
     real(8)              :: TRANSMSURF(NVIS+1:NVIS+NIR)
     real(8)              :: EMISS(NVIS+1:NVIS+NIR)
  end type avhrr_bgck_iasi

  type(avhrr_bgck_iasi)  , allocatable :: avhrr_bgck(:)      ! avhrr parameters for IASI quality control

contains

  subroutine  irbg_init()
   
    implicit none
    integer :: nulnam, ierr
    logical, save :: lfirst = .true.
    integer ,external :: fnom, fclos

    if (lfirst) then
      nulnam = 0
      ierr = fnom(nulnam, './flnml','FTN+SEQ+R/O', 0)
      read(nulnam, nml=NAMBGCKIR, iostat=ierr)
      if (ierr /= 0) call utl_abort('irbg_init: Error reading namelist')
      if (mpi_myid == 0) write(*, nml=NAMBGCKIR)
      ierr = fclos(nulnam)
      lfirst = .false.
    end if
  end subroutine irbg_init


  subroutine irbg_setup(lobsSpaceData)
!
!  s/r irbg_setup : Memory allocation for the Hyperspectral Infrared
!                background check variables
!          (original name of routine: sutovalo)
!
! Revision:

!           S.  Heilliette
!            - creation from tovs_setup_allo  December 2013

    implicit none
!implicits

    type(struct_obs) :: lobsSpaceData

    integer :: alloc_status(2)

    integer :: KRTID
    integer ::  JO,NCMAX
    integer ::  ISENS, NC, NL


!     Memory allocation for background check related variables
!     .  -----------------------------------------------------
    alloc_status(:) = 0
    allocate( tvs_surfaceParameters(tvs_nobtov), stat=alloc_status(1))
    call utl_checkAllocationStatus(alloc_status(1:1), " irbg_setup tvs_surfaceParameters")

!___ radiance by profile
    do jo = 1, tvs_nobtov
      isens = tvs_lsensor(jo)
      nc = tvs_nchan(isens)
      nl = tvs_coefs(isens) % coef % nlevels
 ! allocate clear sky radiance output
      allocate( tvs_radiance(jo)  % clear  ( nc ) ,stat= alloc_status(2) )
      tvs_radiance(jo)  % clear  ( : ) = 0.d0
 !  allocate overcast black cloud sky radiance output
      allocate( tvs_radiance(jo)  % overcast  (nl - 1,nc), stat=alloc_status(1))
      call utl_checkAllocationStatus(alloc_status(1:1), " irbg_setup")
      tvs_radiance(jo)  % overcast  (:,:) = 0.d0
    end do


!___ transmission by profile

    alloc_status(:) = 0
    allocate( tvs_transmission(tvs_nobtov), stat=alloc_status(1))
    call utl_checkAllocationStatus(alloc_status(1:1), " irbg_setup tvs_transmission")

    do jo = 1, tvs_nobtov
      isens = tvs_lsensor(jo)
      nc = tvs_nchan(isens)
      nl = tvs_coefs(isens) % coef % nlevels
      ! allocate transmittance from surface and from pressure levels
      allocate( tvs_transmission(jo) % tau_total ( nc ), stat= alloc_status(1))
      allocate( tvs_transmission(jo) % tau_levels(nl,nc), stat= alloc_status(2))
      call utl_checkAllocationStatus(alloc_status, " irbg_setup")
    end do

!___ emissivity by profile

    ncmax = 1
    do jo = 1, tvs_nobtov
      isens = tvs_lsensor(jo)
      nc = tvs_nchan(isens)
      if (nc > ncmax) ncmax=nc
    end do

    allocate( tvs_emissivity (ncmax,tvs_nobtov), stat=alloc_status(1))
    call utl_checkAllocationStatus(alloc_status(1:1), " irbg_setup tvs_emissivity")
    
    do KRTID = 1, tvs_nsensors

      if ( tvs_instruments(KRTID) == inst_id_iasi ) then
        allocate ( avhrr_bgck(tvs_nobtov), stat=alloc_status(1))
        call utl_checkAllocationStatus(alloc_status(1:1), " irbg_setup avhrr_bgck")
        exit
      end if

    end do

  end subroutine irbg_setup

!--------------------------------------------------------------------------
!! *Purpose*: Do background check on all hyperspectral infrared observations
!!
!! @author P. Koclas *CMC/CMDA  Nov 1998
!!
!--------------------------------------------------------------------------
  subroutine irbg_bgCheckIR(columnhr,obsSpaceData)
    IMPLICIT NONE

    type(struct_obs) :: obsSpaceData
    type(struct_columnData) :: columnhr
    INTEGER J
    INTEGER,allocatable :: nobir(:)
    INTEGER :: index_header,idatyp,krtid

    call tmg_start(3,'BGCHECKIR')

    write(*,'(A)') " ****************"
    write(*,'(A)') " BEGIN IR BACKGROUND CHECK"
    write(*,'(A)') " **************** **************** ****************"

    !
    !     Preliminary initializations
    !     ---------------------------
    !

    tvs_nobtov = 0
    call irbg_init()
    allocate (nobir(ninst))
    nobir = 0
  

    ! loop over all header indices of the 'TO' family
    ! Set the header list (and start at the beginning of the list)
    call obs_set_current_header_list(obsSpaceData,'TO')
    HEADER: do
      index_header = obs_getHeaderIndex(obsSpaceData)
      if (index_header < 0) exit HEADER
      
      IDATYP = obs_headElem_i(obsSpaceData,OBS_ITY,index_header)

      if ( .not. tvs_isIdBurpTovs(IDATYP) ) cycle HEADER   ! Proceed to the next header_index

      tvs_nobtov = tvs_nobtov + 1
      do j=1, ninst
        if ( tvs_isIdBurpInst(IDATYP,INST(j)) ) then
          nobir(j) = nobir(j) + 1
          exit
        end if
      end do
    end do HEADER
    
    do j=1, ninst
      if (nobir(j) > 0) then
        do krtid=1,tvs_nsensors
          if (tvs_instruments(krtid) ==  tvs_getInstrumentId(inst(j)) ) then
            call irbg_doQualityControl (columnhr, obsSpaceData, INST(j), krtid)
          end if
        end do
      end if
    end do
    deallocate (nobir)

  !     Write out contents of obsSpaceData into BURP files
  !
    call obsf_writeFiles(obsSpaceData)
    ! add cloud parameter data to burp files (AIRS,IASI,CrIS,...)
    call ADD_CLOUDPRMS(obsSpaceData)
    do j =1, min(1,obs_numHeader(obsSpaceData))
      call obs_prnthdr(obsSpaceData,j)
      call obs_prntbdy(obsSpaceData,j)
    end do

    ! deallocate obsSpaceData
    call obs_finalize(obsSpaceData)

    call tmg_stop(3)

  END subroutine IRBG_BGCHECKIR

  subroutine add_cloudprms(lobsSpaceData)
    implicit none

    type(struct_obs) :: lobsSpaceData
    integer          :: fileIndex
    character(len=10) :: obsFileType
    
    ! If obs files not split and I am not task 0, then return
    if ( .not.obsf_filesSplit() .and. mpi_myid /= 0 ) return

    do fileIndex = 1, obsf_nfiles

      write(*,*) 'INPUT FILE TO  hir_cldprm_to_brp= ', trim( obsf_cfilnam(fileIndex) )
      call obsf_determineSplitFileType( obsFileType, obsf_cfilnam(fileIndex) )
      if ( trim(obsFileType) /= 'BURP' ) then
        write(*,*) 'obsFileType = ',obsFileType
        call utl_abort('add_cloudprms: this s/r is currently only compatible with BURP files')
      else
        call hir_cldprm_to_brp( lobsspacedata, obsf_cfilnam(fileIndex) )
      end if
    end do

  end subroutine add_cloudprms


  subroutine hir_cldprm_to_brp(lobsspacedata,brp_file)
    IMPLICIT NONE
    !implicits
    CHARACTER(LEN=128),intent(in)     :: BRP_FILE
    type(struct_obs),intent(inout)    :: lobsSpaceData

    TYPE(BURP_FILE)        :: FILE_IN
    TYPE(BURP_RPT)         :: RPT_IN,CP_RPT
    TYPE(BURP_BLOCK)       :: BLOCK_IN
      
    CHARACTER(LEN=9)       :: OPT_MISSING
    INTEGER                :: NEW_BTYP
    INTEGER                :: BTYP10
    INTEGER                :: BTYP10DES,BTYP10INF,BTYP10OBS,BTYP10FLG,BTYP10OMP

    INTEGER                :: NB_RPTS,REF_RPT,REF_BLK,COUNT
    INTEGER, ALLOCATABLE   :: ADDRESS(:), GOODPROF(:)
    REAL(8), ALLOCATABLE   :: BTOBS(:,:)
    REAL(8)                :: ETOP,VTOP,ECF,VCF,HE,ZTS,emisfc
    INTEGER                :: NBELE,NVALE,NTE
    INTEGER, ALLOCATABLE   :: GLBFLAG(:)

    INTEGER                :: I,J,K,KK,L,BTYP,BFAM,ERROR
    INTEGER                :: IND008012,IND012163,IND055200,INDCHAN,ICHN,ICHNB
    INTEGER                :: IDATA2,IDATA3,IDATA,IDATEND
    INTEGER                :: FLAG_PASSAGE1,FLAG_PASSAGE2,FLAG_PASSAGE3
    INTEGER                :: FLAG_PASSAGE4,FLAG_PASSAGE5
    INTEGER                :: IDATYP
    REAL                   :: VAL_OPTION_R4
    CHARACTER(LEN=9)       :: station_id
    
    write(*,*) '---------------------------------------'
    write(*,*) '------- BEGIN hir_cldprm_to_brp -------'
    write(*,*) '---------------------------------------'


    ! initialisation
    ! --------------

    flag_passage1 = 0
    flag_passage2 = 0
    flag_passage3 = 0
    flag_passage4 = 0
    flag_passage5 = 0
    
    opt_missing = 'MISSING'
    val_option_r4  = -7777.77

    call BURP_Set_Options( &
         REAL_OPTNAME       = opt_missing, &
         REAL_OPTNAME_VALUE = val_option_r4, &
         IOSTAT             = error )

    call BURP_Init(File_in,IOSTAT=error)
    call BURP_Init(Rpt_in,Cp_rpt,IOSTAT=error)
    call BURP_Init(Block_in,IOSTAT=error)


    ! opening file
    ! ------------

    write(*,*) 'OPENED FILE = ', trim(brp_file)

    call BURP_New(File_in, &
         FILENAME = brp_file, &
         MODE     = FILE_ACC_APPEND, &
         IOSTAT   = error )


    ! obtain input burp file number of reports
    ! ----------------------------------------

    call BURP_Get_Property(File_in, NRPTS=nb_rpts)


    ! scan input burp file to get all reports address
    ! -----------------------------------------------

    allocate(address(nb_rpts))
    address(:) = 0
    count = 0
    ref_rpt = 0

    do
      ref_rpt = BURP_Find_Report(File_in, &
           REPORT      = Rpt_in,  &
           SEARCH_FROM = ref_rpt, &
           IOSTAT      = error)
      if (ref_rpt < 0) Exit

      call BURP_Get_Property(Rpt_in,STNID =station_id)
      if (station_id(1:2)==">>") cycle

      count = count + 1
      address(count) = ref_rpt
    end do

    write(*,*) 
    write(*,*) 'NUMBER OF REPORTS WITH OBSERVATIONS = ',count
    write(*,*) 
    
    if ( count > 0 ) then

      ! create a new report
      ! ------------------
      
      !pik  call BURP_New(Cp_rpt, ALLOC_SPACE=10000000, IOSTAT=error)
      call BURP_New(Cp_rpt, ALLOC_SPACE=20000000, IOSTAT=error)
      if (error/=burp_noerr) then
        Write(*,*) "Error creating new directory ",error 
        call handle_error('hir_cldprm_to_brp')
      end if

      ! LOOP ON REPORTS
      ! ---------------

      REPORTS: do kk = 1, count

        call BURP_Get_Report(File_in, &
             REPORT    = Rpt_in, &
             REF       = address(kk), &
             IOSTAT    = error)
        
        if (kk == 1) then
          call BURP_Get_Property(Rpt_in, IDTYP = IDATYP)
          Write(*,*) "hir_cldprm_to_brp idatyp ", IDATYP
          idata2 = -1
          call obs_set_current_header_list(lobsSpaceData,'TO')
          HEADER: do
            i = obs_getHeaderIndex(lobsSpaceData)
            if (i < 0) exit HEADER  
            if  ( obs_headElem_i(lobsSpaceData,OBS_ITY,i) == idatyp) then
              idata2 = i
              exit HEADER
            end if
          end do HEADER
          if (idata2 == -1) then
            Write(*,*) "datyp ",idatyp," not found in input file !"
            call utl_abort("hir_cldprm_to_brp")
          end if
          idata3 = idata2
        end if

        ! FIRST LOOP ON BLOCKS
        ! --------------------

        ! find bad profiles not in CMA. This occurs if :
        !  - all observations are -1 and/or have a quality flag not zero


        ref_blk = 0

        BLOCKS1: do

          ref_blk = BURP_Find_Block(Rpt_in, &
               BLOCK       = Block_in, &
               SEARCH_FROM = ref_blk, &
               IOSTAT      = error)

          if (ref_blk < 0) EXIT BLOCKS1

          call BURP_Get_Property(Block_in, &
               NELE   = nbele, &
               NVAL   = nvale, &
               NT     = nte,   &
               BFAM   = bfam,  &
               BTYP   = btyp,  &
               IOSTAT = error)

          ! observation block (btyp = 0100 100011X XXXX)
          ! 0100 1000110 0000 = 9312
          btyp10    = ishft(btyp,-5)
          btyp10obs = 291
           
          if ( btyp10 - btyp10obs == 0 .and. bfam == 0 ) then

            ALLOCATE(goodprof(nte),btobs(nvale,nte))

            goodprof(:) = 0
            btobs(:,:)  = 0.

            ind012163  = BURP_Find_Element(Block_in, ELEMENT=012163, IOSTAT=error)

            do k=1,nte
              do j=1,nvale
                btobs(j,k) =  BURP_Get_Rval(Block_in, &
                     NELE_IND = ind012163, &
                     NVAL_IND = j, &
                     NT_IND   = k )
                if ( btobs(j,k) > 0. ) goodprof(k) = 1
              end do
            end do
            
          end if

        end do BLOCKS1


        call BURP_Copy_Header(TO=Cp_rpt,FROM=Rpt_in)
        IF (error /= BURP_NOERR) then
          Write(*,*) "Error= ",error
          call handle_error("Erreur dans BURP_Copy_Header")
        end if

        call BURP_Init_Report_Write(File_in,Cp_Rpt, IOSTAT=error)
        IF (error /= BURP_NOERR) then
          Write(*,*) "Error= ",error
          call handle_error("Erreur dans BURP_Init_Report_Write")
        end if

        ! SECOND LOOP ON BLOCKS
        ! ---------------------

        ! add new informations


        ref_blk = 0
        
        BLOCKS2: do

          if ( .not. allocated(goodprof) ) then
            write(*,*)
            write(*,*) 'Resume report is position # ',kk
            EXIT BLOCKS2
          end if

          ref_blk = BURP_Find_Block(Rpt_in, &
               BLOCK       = Block_in, &
               SEARCH_FROM = ref_blk, &
               IOSTAT      = error)
          
          if (ref_blk < 0) EXIT BLOCKS2

          call BURP_Get_Property(Block_in, &
               NELE   = nbele, &
               NVAL   = nvale, &
               NT     = nte, &
               BFAM   = bfam, &
               BTYP   = btyp, &
               IOSTAT = error)
          

          ! descriptor block (btyp = 0010 100000X XXXX) 
          ! 0010 1000000 0000==5120 )
          !    if profile contains rejected observations (apart from blacklisted channels),
          !     set bit 6 in global flags.

          btyp10    = ishft(btyp,-5)
          btyp10des = 160

          if ( btyp10 - btyp10des == 0 ) then

            flag_passage1 = 1

            ALLOCATE(glbflag(nte))

            ind055200  = BURP_Find_Element(Block_in, ELEMENT=055200, IOSTAT=error)
            do k = 1, nte
              glbflag(k) =  BURP_Get_Tblval(Block_in, &
                   NELE_IND = ind055200, &
                   NVAL_IND = 1, &
                   NT_IND   = k )
            end do

            do k = 1, nte
              if (goodprof(k)/= 1) glbflag(k) = ibset(glbflag(k),6)            
            end do

            do k = 1, nte
              call BURP_Set_Tblval(Block_in, &
                   NELE_IND = ind055200, &
                   NVAL_IND = 1, &
                   NT_IND   = k, &
                   TBLVAL   = glbflag(k), &
                   IOSTAT   = error)
            end do
              
            DEALLOCATE(glbflag)

          end if


          ! info block (btyp = 0001 100000X XXXX) 
          ! 0001 100000X XXXX = 3072
          btyp10    = ishft(btyp,-5)
          btyp10inf = 96

          if ( btyp10 - btyp10inf == 0 ) then

            flag_passage2 = 1

            call BURP_Resize_Block(Block_in, ADD_NELE=11, IOSTAT=error)
            if (error/=burp_noerr) then
              call handle_error("Erreur dans BURP_Resize_Block info")
            end if
            call BURP_Set_Element(Block_in, NELE_IND=nbele+ 1, ELEMENT=014213, IOSTAT=error)
            call BURP_Set_Element(Block_in, NELE_IND=nbele+ 2, ELEMENT=014214, IOSTAT=error)
            call BURP_Set_Element(Block_in, NELE_IND=nbele+ 3, ELEMENT=014215, IOSTAT=error)
            call BURP_Set_Element(Block_in, NELE_IND=nbele+ 4, ELEMENT=014216, IOSTAT=error)
            call BURP_Set_Element(Block_in, NELE_IND=nbele+ 5, ELEMENT=014217, IOSTAT=error)
            call BURP_Set_Element(Block_in, NELE_IND=nbele+ 6, ELEMENT=014218, IOSTAT=error)
            call BURP_Set_Element(Block_in, NELE_IND=nbele+ 7, ELEMENT=014219, IOSTAT=error)
            call BURP_Set_Element(Block_in, NELE_IND=nbele+ 8, ELEMENT=014220, IOSTAT=error)
            call BURP_Set_Element(Block_in, NELE_IND=nbele+ 9, ELEMENT=014221, IOSTAT=error)
            call BURP_Set_Element(Block_in, NELE_IND=nbele+10, ELEMENT=013214, IOSTAT=error)
            call BURP_Set_Element(Block_in, NELE_IND=nbele+11, ELEMENT=059182, IOSTAT=error)
            
            ind008012 = BURP_Find_Element(Block_in, &
                 ELEMENT  = 008012, &
                 IOSTAT   = error)
            
            do k = 1, nte
              
              if ( goodprof(k) == 1 ) then

                call Insert_into_burp_r8(obs_headElem_r(lobsSpaceData,OBS_ETOP,idata2),nbele+1,1,k)

                call Insert_into_burp_r8(obs_headElem_r(lobsSpaceData,OBS_VTOP,idata2),nbele+2,1,k)

                call Insert_into_burp_r8(obs_headElem_r(lobsSpaceData,OBS_ECF,idata2),nbele+3,1,k)
                
                call Insert_into_burp_r8(obs_headElem_r(lobsSpaceData,OBS_VCF,idata2),nbele+4,1,k)
                
                call Insert_into_burp_r8(obs_headElem_r(lobsSpaceData,OBS_HE,idata2),nbele+5,1,k)
                
                call Insert_into_burp_r8(obs_headElem_r(lobsSpaceData,OBS_ZTSR,idata2),nbele+6,1,k)
                
                call Insert_into_burp_i(obs_headElem_i(lobsSpaceData,OBS_NCO2,idata2),nbele+7,1,k)
                
                call Insert_into_burp_r8(obs_headElem_r(lobsSpaceData,OBS_ZTM,idata2),nbele+8,1,k)
                
                call Insert_into_burp_r8(obs_headElem_r(lobsSpaceData,OBS_ZTGM,idata2),nbele+9,1,k)
                
                call Insert_into_burp_r8(obs_headElem_r(lobsSpaceData,OBS_ZLQM,idata2),nbele+10,1,k)
                
                call Insert_into_burp_r8(obs_headElem_r(lobsSpaceData,OBS_ZPS,idata2),nbele+11,1,k)
                
                call Insert_into_burp_i(obs_headElem_i(lobsSpaceData,OBS_STYP,idata2),ind008012,1,k)
                                
                idata2 = idata2 + 1

              else

                do i = 1, 11
                  call Insert_into_burp_r8(-1.d0,nbele + i,1,k)
                end do

                call Insert_into_burp_i(-1,ind008012,1,k)
                                
              end if
                             
            end do

          end if


          ! observation block (btyp = 0100 100011X XXXX)
          ! 0100 1000110 0000 = 9312
          btyp10    = ishft(btyp,-5)
          btyp10obs = 291

          if ( btyp10 - btyp10obs == 0 .and. bfam == 0 ) then
            flag_passage3 = 1

            call BURP_Resize_Block(Block_in, ADD_NELE=1, IOSTAT=error)
            if (error/=burp_noerr) then
              call handle_error("Erreur dans BURP_Resize_Block data")
            end if
            call BURP_Set_Element(Block_in, NELE_IND=nbele+1, ELEMENT=055043, IOSTAT=error)
            indchan  = BURP_Find_Element(Block_in, ELEMENT=005042, IOSTAT=error)
            do k = 1, nte
              do j = 1, nvale
                call Insert_into_burp_i(-1,nbele+1,j,k)
              end do
                 
              if ( goodprof(k) == 1 ) then

                IDATA   = obs_headElem_i(lobsSpaceData,OBS_RLN,idata3)
                IDATEND = obs_headElem_i(lobsSpaceData,OBS_NLV,idata3) + IDATA - 1
                do j = IDATA,IDATEND
                  emisfc=100.d0*obs_bodyElem_r(lobsspacedata,OBS_SEM,j)
                  ICHN = NINT(obs_bodyElem_r(lobsSpaceData,OBS_PPP,j))
                  ICHN = MAX(0,MIN(ICHN,tvs_maxChannelNumber+1))
                  bl: do l=1,nvale
                    ichnb=BURP_Get_Tblval(Block_in, &
                         NELE_IND = indchan, &
                         NVAL_IND = l, &
                         NT_IND   = k)
                    if (ichn==ichnb) then
                      call Insert_into_burp_r8(emisfc,nbele+1,l,k)
                      exit bl
                    end if
                  end do bl
                  
                end do
                    
                idata3 = idata3 + 1
                    
              end if
                       
            end do

          end if


          ! flag block (btyp = 0111 100011X XXXX)
          ! 0111 1000110 0000 = 15456
          btyp10    = ishft(btyp,-5)
          btyp10flg = 483
              
          if ( btyp10 - btyp10flg == 0 ) then
            flag_passage4 = 1
            
            call BURP_Resize_Block(Block_in, ADD_NELE=1, IOSTAT=error)
            if (error/=burp_noerr) then
              call handle_error("Erreur dans BURP_Resize_Block marqueur")
            end if
            call BURP_Set_Element(Block_in, NELE_IND=nbele+1, ELEMENT=255043, IOSTAT=error)
            
            do k = 1, nte
              do j = 1, nvale
                call BURP_Set_Tblval(Block_in, &
                     NELE_IND = nbele+1, &
                     NVAL_IND = j, &
                     NT_IND   = k, &
                     TBLVAL   = 0, &
                     IOSTAT   = error)
              end do
            end do
          end if
              

          ! O-P block (btyp = 0100 100011X XXXX)
          ! 0100 1000110 0000 = 9312
          btyp10    = ishft(btyp,-5)
          btyp10omp = 291
          
          if ( btyp10 - btyp10omp == 0 .and. bfam == 14 ) then
            flag_passage5 = 1
              
            call BURP_Resize_Block(Block_in, ADD_NELE=1, IOSTAT=error)
            if (error/=burp_noerr) then
              call handle_error("Erreur dans BURP_Resize_Block O-P")
            end if
            call BURP_Set_Element(Block_in, NELE_IND=nbele+1, ELEMENT=055043, IOSTAT=error)
                
            do k = 1, nte
              do j = 1, nvale
                call Insert_into_burp_i(-1,nbele+1,j,k)
              end do
            end do
                
          end if

          ! add block into new report
          ! -------------------------

          if ( btyp == 5120 ) then
            call BURP_Write_Block(Cp_rpt, Block_in, &
                 ENCODE_BLOCK  = .true., &
                 IOSTAT        = error)
          else
            call BURP_Write_Block(Cp_rpt, Block_in, &
                 ENCODE_BLOCK  = .true., &
                 CONVERT_BLOCK = .true., &
                 IOSTAT        = error)
          end if
          if (error/=burp_noerr) then
            write(*,*)"Btyp= ",btyp
            call handle_error("Erreur dans BURP_Write_Block")
          end if
        end do BLOCKS2


        if ( allocated(goodprof) ) then
          DEALLOCATE (goodprof,btobs)
        end if


        ! write new report into file
        ! --------------------------
        
        call BURP_Delete_Report(File_in,Rpt_in, IOSTAT=error)
        call BURP_Write_Report(File_in,Cp_rpt, IOSTAT=error)
      end do REPORTS

      if ( flag_passage1 == 0 ) then
        write(*,*)
        write(*,*) 'ERROR - descriptor block not seen ? Verify btyp'
      end if
      if ( flag_passage2 == 0 ) then
        write(*,*)
        write(*,*) 'ERROR - info block not seen ? Verify btyp'
      end if
      if ( flag_passage3 == 0 ) then
        write(*,*)
        write(*,*) 'ERROR - observation block not seen ? Verify btyp'
      end if
      if ( flag_passage4 == 0 ) then
        write(*,*)
        write(*,*) 'ERROR - flag block not seen ? Verify btyp'
      end if
      if ( flag_passage5 == 0 ) then
        write(*,*)
        write(*,*) 'ERROR - O-P block not seen ? Verify btyp'
      end if
          
    end if !! End of 'if ( count > 0 )'

    deallocate(address)

    call BURP_Free(File_in,IOSTAT=error)
    call BURP_Free(Rpt_in,Cp_rpt,IOSTAT=error)
    call BURP_Free(Block_in,IOSTAT=error)

  contains

    !------------------------------------- HANDLE_ERROR -----
  
    subroutine handle_error(errormessage)
      implicit none
      character (len=*) :: errormessage
      write(*,*) BURP_STR_ERROR()
      write(*,*) "history"
      call BURP_STR_ERROR_HISTORY()
      if (allocated(address)) deallocate(address)
      call BURP_Free(File_in)
      call BURP_Free(Rpt_in,Cp_rpt)
      call BURP_Free(Block_in)
      call utl_abort(trim(errormessage))
    end subroutine handle_error

    subroutine Insert_into_burp_r8(r8val,pele,pval,pt)
      implicit none
      real (8), intent(in):: r8val
      integer, intent(in) :: pele,pval,pt
      integer :: error
      
      if ( r8val >= 0.d0 ) then
        call BURP_Set_Rval(Block_in, &
             NELE_IND = pele, &
             NVAL_IND = pval, &
             NT_IND   = pt, &
             RVAL     = sngl(r8val), &
             IOSTAT   = error)
      else
        call BURP_Set_Rval(Block_in, &
             NELE_IND = pele, &
             NVAL_IND = pval, &
             NT_IND   = pt, &
             RVAL     = val_option_r4, &
             IOSTAT   = error)
      end if
      if (error/=burp_noerr) then
        Write(*,*) "r8val,pele,pval,pt",r8val,pele,pval,pt
        call handle_error("Insert_into_burp_r8")
      end if

    end subroutine Insert_into_burp_r8

    subroutine Insert_into_burp_i(ival,pele,pval,pt)
      implicit none
      integer, intent(in) :: ival
      integer, intent(in) :: pele,pval,pt
      integer :: error
      
      if ( ival >= 0 ) then
        call BURP_Set_Rval(Block_in, &
             NELE_IND = pele, &
             NVAL_IND = pval, &
             NT_IND   = pt, &
             RVAL   = real(ival), &
             IOSTAT   = error)
      else
        call BURP_Set_Rval(Block_in, &
             NELE_IND = pele, &
             NVAL_IND = pval, &
             NT_IND   = pt, &
             RVAL   = val_option_r4, &
             IOSTAT   = error)
      end if
      
      if (error/=burp_noerr) then
        Write(*,*) "ival,pele,pval,pt",ival,pele,pval,pt
        call handle_error("Insert_into_burp_i")
      end if
    end subroutine Insert_into_burp_i


  END subroutine HIR_CLDPRM_TO_BRP


  subroutine BGCK_GET_QCID(CINSTR,QCID)
    implicit none
    character (len=*),intent(in) :: CINSTR
    integer ,intent (out) :: QCID
    !**********
    integer :: i 

    QCID = -1

    do i=1, NINST
      if (trim(CINSTR) == trim(INST(i))) then
        QCID = i
        exit
      end if
    end do

    if (QCID == -1) then
      Write(*,*) "Unknown instrument ",CINSTR
      call utl_abort('BGCK_GET_QCID')
    end if

  end subroutine BGCK_GET_QCID

  subroutine irbg_doQualityControl ( lcolumnhr, lobsSpaceData,CINST,id_opt)
!
!**ID irbg_doQualityControl -- QUALITY CONTROL OF HYPERSPECTRAL INFRARED OBSERVATIONS
!
!       SCIENCE:  L. GARAND
!       AUTHOR:   A. BEAULNE (CMDA/SMC) August 2004
!                 A. BEAULNE (CMDA/SMC)   June 2006  (ADAPT TO 3DVAR)
!                 S. HEILLIETTE           February 2008 (adaptation to IASI)
!                 S. MACPHERSON, S.HEILLIETTE (ARMA) February 2013 
!                   -- modify test pour detecter le isatzen manquant ou anormal
!
!       REVISION:
!
!       OBJECT: ASSIGN ASSIMILATION FLAGS TO OBSERVATIONS 
!
!       ARGUMENTS:
!          INPUT:
!            -LOOP_DONE : NUMBER OF PREVIOUS CALLS TO irbg_doQualityControl
!
!          OUTPUT:
!            -LEND       : AT THE END OF THIS CALL TO irbg_doQualityControl, DO ALL 
!                               PROFILES BEEN TREATED (true) OR NOT (false)
!
    implicit none
    integer,intent(in),optional :: id_opt
    type(struct_columnData),intent(in) :: lcolumnhr
    type(struct_obs),intent(inout) :: lobsSpaceData
    character (len=*),intent(in) :: CINST
!******************************************************************
    integer       :: JC,NCHN,JCH,JF,JL,NLEV,NLEVB,iextr,NPRF,NFLG,ICHN
    integer       :: IWINDO,IWINDO_ALT
    integer       :: INDEX_BODY,IDATA,IDATEND,INDEX_HEADER
    integer       :: IDATYP
    real(8)       :: DIFFTOP_MIN
    integer       :: IMODTOP
    integer       :: count
    real(8)       :: T_EFFECTIVE
    integer       :: alloc_status(28)

    real(8) :: ZTG,ZPS,ZTS,ptop_T,ZLQS
    real(8), allocatable :: ZT(:),ZHT(:,:),ZVLEV(:)
    real(8), allocatable :: ZLEVMOD(:,:)
    real(8), allocatable :: BTOBSERR(:),BTOBS(:),BTCALC(:),RCAL_CLR(:),SFCTAU(:)
    real(8), allocatable :: ROBS(:),RCLD(:,:),TRANSM(:,:),EMI_SFC(:) 
    real(8), allocatable :: TOEXT(:),ZHOEXT(:,:)
    real(8), allocatable :: PTOP_BT(:),PTOP_RD(:)
    real(8), allocatable :: PMIN(:),DTAUDP1(:),MAXWF(:)
    real(8), allocatable :: RCLD_AVHRR(:,:)
    integer, allocatable :: REJFLAG(:,:) 
    integer, allocatable :: NTOP_BT(:),NTOP_RD(:)
    integer, allocatable :: MINP(:),FATE(:), channelIndex(:)
    real(8), allocatable :: xpres(:)

    real(8) :: CLFR,SUNZA,SATAZIM,SATZEN,SUNAZIM
    real(8) :: ALBEDO,ICE,PCNT_WAT,PCNT_REG
    real(8) :: PTOP_EQ,PTOP_MB
    real(8) :: PTOP_CO2(NCO2),FCLOUD_CO2(NCO2)
    real(8) :: ETOP,VTOP,ECF,VCF,HEFF
    real(8) :: TAMPON,CFSUB
    real(8) :: ZTS_AVHRR(nClassAVHRR),SFCTAU_AVHRR(NIR),EMI_SFC_AVHRR(NIR),RCAL_CLR_AVHRR(NIR)
    real(8) :: PTOP_BT_AVHRR(NIR,nClassAVHRR),PTOP_RD_AVHRR(NIR,nClassAVHRR)
    real(8) :: BTOBS_AVHRR(NIR,nClassAVHRR),ROBS_AVHRR(NIR,nClassAVHRR),PTOP_EQ_AVHRR(nClassAVHRR)
    real(8) :: CFRAC_AVHRR
    real(8) :: avhrr_surfem1(NIR)
    real(8) :: seuil_albed(NIR)

    integer :: KSURF,LTYPE
    integer :: CLDFLAG,LEV_START   
    integer :: GNCLDFLAG
    integer :: ICHREF,INDX(1)
    integer :: NTOP_EQ,NTOP_MB
    integer :: NGOOD
    integer :: NTOP_CO2(NCO2)
    integer :: CLDFLAG_AVHRR(nClassAVHRR),LEV_START_AVHRR(nClassAVHRR),ICHREF_AVHRR(nClassAVHRR),NTOP_RD_AVHRR(NIR,nClassAVHRR)
    integer :: NTOP_BT_AVHRR(NIR,nClassAVHRR),NTOP_EQ_AVHRR(nClassAVHRR)
    integer :: ICL

    logical :: ASSIM_ALL
  
    integer ,parameter :: nn=2
    integer ,parameter :: ilist_avhrr(nn)=(/ 2 ,3 /)
    integer :: cpt,iclass
    logical :: bad
    real(8),parameter :: sunzenmax=87.12d0
    real(8) :: minpavhrr(2:3)
    real(8) :: anisot,zlamb,zcloud,scos,del,deltaphi
    integer :: ier,ijour,iloc(2:3),co2min(1),co2max(1),iobs
    integer :: isatzen
    integer :: chan_indx,ILIST_SUN,ilist_co2(NCO2),ilist_co2_pair(NCO2),ilist_he(NCH_HE)
!***************************************************************************************
    integer :: nlv_T,id,KRTID, QCID, nchannels
    logical :: liasi,lairs,lcris
!****************************************
    write (*,*) "Entering irbg_doQualityControl"

    liasi= ( trim(cinst) == "IASI" .or.  trim(cinst) == "iasi")
    lairs= ( trim(cinst) == "AIRS" .or.  trim(cinst) == "airs")
    lcris= ( trim(cinst) == "CRIS" .or.  trim(cinst) == "cris")

    call BGCK_GET_QCID(cinst,QCID)
    
    if (present(id_opt)) then
      id = id_opt
    else
! ** find sensor number corresponding to the desired instrument
      ID = -1
      do KRTID = 1, tvs_nsensors
        if ( trim(tvs_instrumentName(KRTID)) == TRIM(CINST)) then
          ID = KRTID
          exit
        end if
      end do
      if (ID < 0) call utl_abort("irbg_doQualityControl: should not happen !")
    end if

    ! ** find number of profiles 
    count = 0

    ! loop over all header indices of the 'TO' family
    call obs_set_current_header_list(lobsSpaceData,'TO')
    HEADER: do
      index_header = obs_getHeaderIndex(lobsSpaceData)
      if (index_header < 0) exit HEADER
       
      IDATYP = obs_headElem_i(lobsSpaceData,OBS_ITY,INDEX_HEADER)
      if ( tvs_isIdBurpInst(IDATYP,CINST) .and. tvs_lsensor(tvs_ltovsno (index_header)) == id ) then
        count = count + 1
      end if
    end do HEADER

    if ( count == 0 ) return
    ! ** find number of channels and RTTOV levels

    NCHN = tvs_coefs(id)%coef%fmv_chn

    NLEV = tvs_coefs(id) % coef % nlevels 
    allocate (xpres(NLEV))
    xpres(1:NLEV) = tvs_coefs(id) % coef % ref_prfl_p(1:NLEV)

    select case(nlev)
    case(43)
      iextr = 0
    case(44)
      iextr = 1
    case(51)
      iextr = 2
    case(54)
      iextr = 4
    case(101)
      iextr = 4
    case default
      Write(*,*) "Attention: modification necessaire dans irbg_doQualityControl", nlev
      call utl_abort('irbg_doQualityControl')
    end select
   
    NLEVB = NLEV - iextr

    write(*,*) ' irbg_doQualityControl - nchn ', nchn
  
    nlv_T = col_getNumLev(lcolumnhr,'TH')
   

! information to extract (transvidage)
! ------------------------------------
!
! ZTG -- guess skin temperatures (deg K)
! ZPS(NPRF) -- surface pressure (hPa)
! ZT(nlv_T) -- temperature profiles on NWP model levels (deg K)
! ZHT(nlv_T,1) -- height profiles on NWP model levels (m)
! ZLQS -- surface specific humidity in ln q (kg/kg)
! BTOBSERR(nchn) -- observation error standard deviation
! BTOBS(nchn) -- observed brightness temperatures (deg K)
! BTCALC(nchn) -- computed brightness temperatures (deg K)
! RCAL_CLR(nchn) -- computed clear radiances (mw/m2/sr/cm-1)
! SFCTAU(nchn) -- surface to space transmittances (0-1)
! RCLD(nchn,NLEV) -- overcast cloudy radiances (mw/m2/sr/cm-1)
! TRANSM(nchn,NLEV) -- layer to space transmittances (0-1)
! EMI_SFC(nchn) -- surface emissivities (0-1)
! KSURF -- surface type in obs file (0, 1)
! CLFR -- cloud fraction (%)
! TOEXT(NLEV) -- temperature profiles on RT model levels (deg K)
! ZHOEXT(NLEV) -- height profiles on RT model levels (m)
! SUNZA -- sun zenith angle (deg)
! SATAZIM -- satellite azimuth angle (deg)
! SATZEN -- satellite zenith angle (deg)
! ALBEDO -- surface albedo (0-1)
! ICE -- ice fraction (0-1)
! LTYPE -- surface type (1,...,20)
! PCNT_WAT -- water fraction (0-1)
! PCNT_REG -- water fraction in the area (0-1)
! ROBS(nchn) -- observed radiances (mW/m2/sr/cm-1)

    alloc_status(:) = 0
 
    allocate ( BTOBSERR(nchn),                   stat= alloc_status(1))
    allocate ( BTOBS(nchn),                      stat= alloc_status(2))
    allocate ( BTCALC(nchn),                     stat= alloc_status(3))
    allocate ( RCAL_CLR(nchn),                   stat= alloc_status(4))
    allocate ( SFCTAU(nchn),                     stat= alloc_status(5))
    allocate ( RCLD(nchn,NLEVB),                 stat= alloc_status(6))
    allocate ( TRANSM(nchn,NLEVB),               stat= alloc_status(7))
    allocate ( EMI_SFC(nchn),                    stat= alloc_status(8))
    allocate ( TOEXT(NLEVB),                     stat= alloc_status(9))
    allocate ( ZHOEXT(NLEVB,1),                  stat= alloc_status(10))
    allocate ( ROBS(nchn),                       stat= alloc_status(11))
    allocate ( REJFLAG(nchn,0:BITFLAG),          stat= alloc_status(12))
    allocate ( NTOP_BT(nchn),                    stat= alloc_status(13))
    allocate ( NTOP_RD(nchn),                    stat= alloc_status(14))
    allocate ( PTOP_BT(nchn),                    stat= alloc_status(15))
    allocate ( PTOP_RD(nchn),                    stat= alloc_status(16))
    allocate ( MINP(nchn),                       stat= alloc_status(17))
    allocate ( PMIN(nchn),                       stat= alloc_status(18))
    allocate ( DTAUDP1(nchn),                    stat= alloc_status(19))
    allocate ( FATE(nchn),                       stat= alloc_status(20))
    if (liasi) allocate ( RCLD_AVHRR(NIR,NLEVB), stat= alloc_status(21))
    allocate ( maxwf(nchn),                      stat= alloc_status(22))
    allocate ( ZVLEV(NLEVB),                     stat= alloc_status(23))
    allocate ( ZLEVMOD(nlv_T,1),                 stat= alloc_status(24))
    allocate ( ZT(nlv_T),                        stat= alloc_status(25))
    allocate ( ZHT(nlv_T,1),                     stat= alloc_status(26))
    allocate ( channelIndex(nchn),               stat= alloc_status(27))
    call utl_checkAllocationStatus(alloc_status, " irbg_doQualityControl 1")

    do JL = 1, NLEVB
      ZVLEV(JL) = XPRES(JL + iextr)
    end do

  
    DIFFTOP_MIN = 100000.d0
    IMODTOP = 1

    ptop_T = col_getPressure(lcolumnhr,1,1,'TH')
    do JL = 1, NLEVB
      if ( abs(ptop_T - 100.d0 * ZVLEV(JL)) < DIFFTOP_MIN ) then
        DIFFTOP_MIN = abs(ptop_T - 100.d0 * ZVLEV(JL))
        IMODTOP = JL
      end if
    end do
!* -- FIND RADIATIVE TRANSFER MODEL LEVEL NEAREST TO TRIAL TOP (only compute one time)
    write(*,*) 'TOIT DU MODELE (MB)'
    write(*,*) 0.01d0 * ptop_T
    write(*,*) 'NIVEAU DU MODELE DE TRANSFERT RADIATIF LE PLUS PRES DU TOIT DU MODELE'
    write(*,*) IMODTOP

    CO2MIN = minloc( abs( ZVLEV(:) - pco2min ) )
    CO2MAX = minloc( abs( ZVLEV(:) - pco2max ) )

    tvs_nobtov = 0

    ! loop over all header indices of the 'TO' family
    call obs_set_current_header_list(lobsSpaceData, 'TO')
    HEADER_2: do
      index_header = obs_getHeaderIndex(lobsSpaceData)
      if (index_header < 0) exit HEADER_2

      IDATYP = obs_headElem_i(lobsSpaceData,OBS_ITY,INDEX_HEADER)

      if ( tvs_isIdBurpTovs(idatyp) ) tvs_nobtov = tvs_nobtov + 1

      if ( tvs_isIdBurpInst(IDATYP,CINST) .and. tvs_lsensor(tvs_ltovsno (index_header)) == id) then
        BTOBS(:)    = -1.d0
        BTCALC(:)   = -1.d0
        BTOBSERR(:) = -1.d0
        RCAL_CLR(:) = -1.d0
        SFCTAU(:)   = -1.d0
        RCLD(:,:)   = -1.d0
        TRANSM(:,:) = -1.d0
        EMI_SFC(:)  = -1.d0
        REJFLAG(:,:) = 0
        channelIndex(:) = -1

        if (liasi) then
          INDX = index_header
          iclass = 1
          do iobs=OBS_CF1, OBS_CF7
            avhrr_bgck(INDEX_HEADER)%CFRAC(iclass) = obs_headElem_i(lobsSpaceData,iobs,index_header)
            iclass = iclass + 1
          end do
          iclass = 1
          ichn = 1
          do iobs = OBS_M1C1, OBS_M7C6
            avhrr_bgck(INDEX_HEADER)%radmoy(iclass,ichn) = obs_headElem_r(lobsSpaceData,iobs,index_header)
            ichn = ichn + 1
            if (ichn > nChanAVHRR) then
              ichn = 1
              iclass = iclass + 1
            end if
          end do
          iclass = 1
          ichn = 1
          do iobs=OBS_S1C1, OBS_S7C6
            avhrr_bgck(INDEX_HEADER)%radstd(iclass,ichn) = obs_headElem_r(lobsSpaceData,iobs,index_header)
            ichn = ichn + 1
            if (ichn > nChanAVHRR) then
              ichn = 1
              iclass = iclass + 1
            end if
          end do
          SUNAZIM = 0.01d0 * obs_headElem_i(lobsSpaceData,OBS_SAZ,index_header)
        end if

        ZTG = col_getElem(lcolumnhr,1,INDEX_HEADER,'TG')
        ZPS = col_getElem(lcolumnhr,1,INDEX_HEADER,'P0') * MPC_MBAR_PER_PA_R8

        do JL = 1, nlv_T
          ZT(JL) = col_getElem(lcolumnhr,JL,INDEX_HEADER,'TT')
          ZHT(JL,1) = col_getHeight(lcolumnhr,JL,INDEX_HEADER,'TH') / RG
          ZLEVMOD(JL,1)= col_getPressure(lcolumnhr,JL,INDEX_HEADER,'TH') * MPC_MBAR_PER_PA_R8
        end do
        ZLQS = col_getElem(lcolumnhr,nlv_T,INDEX_HEADER,'HU')

        call ppo_lintv (zlevmod(:,1:1),zht(:,1:1),nlv_T,nlv_T,1, &
             nlevb,zvlev,zhoext(:,1:1))

        IDATA   = obs_headElem_i(lobsSpaceData,OBS_RLN,index_header)
        IDATEND = obs_headElem_i(lobsSpaceData,OBS_NLV,index_header) + IDATA - 1
        BAD = .false.
        if (lcris) BAD=( obs_headElem_i(lobsSpaceData,OBS_GQF,index_header)/=0 .or. &
             obs_headElem_i(lobsSpaceData,OBS_GQL,index_header) /=0)
        if (liasi) BAD=( obs_headElem_i(lobsSpaceData,OBS_GQF,index_header)/=0 .or. &
             obs_headElem_i(lobsSpaceData,OBS_GQL,index_header) >1) 


        nchannels = 0 ! number of channels available at that observation point
        do INDEX_BODY= IDATA, IDATEND
          if ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY)==1 ) then
            ICHN = nint(obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY))
            ICHN = max( 0,min( ICHN,tvs_maxChannelNumber + 1 ) )
            call tvs_getChannelIndexFromChannelNumber(id,chan_indx,ichn)
            nchannels = nchannels + 1
            channelIndex(nchannels) = chan_indx
!            channelNumber(nchannels) = ICHN
            BTOBSERR(chan_indx) = obs_bodyElem_r(lobsSpaceData,OBS_OER,INDEX_BODY)
            BTOBS(chan_indx) = obs_bodyElem_r(lobsSpaceData,OBS_VAR,INDEX_BODY)

            ! *** Flag check on observed BTs ***
            if (.not.liasi .and. btest(obs_bodyElem_i(lobsSpaceData,OBS_FLG,INDEX_BODY),2)) REJFLAG(chan_indx,9) = 1
            if (BAD) REJFLAG(chan_indx,9) = 1

              ! *** Gross check on observed BTs ***
            if (BTOBS(chan_indx)<150.d0) REJFLAG(chan_indx,9) = 1
            if (BTOBS(chan_indx)>350.d0) REJFLAG(chan_indx,9) = 1
          end if
        end do

        if (nchannels==0) cycle HEADER_2
        do JC = 1, nchannels
          chan_indx = channelIndex(jc)
          BTCALC(chan_indx) = tvs_radiance(tvs_nobtov) % bt(chan_indx)
          RCAL_CLR(chan_indx) = tvs_radiance(tvs_nobtov) % clear(chan_indx)
          SFCTAU(chan_indx) = tvs_transmission(tvs_nobtov) % tau_total(chan_indx)
          do JL = 1, NLEVB
            RCLD(chan_indx,JL) = tvs_radiance(tvs_nobtov) % overcast(jl + iextr - 1,chan_indx)
            TRANSM(chan_indx,JL) = tvs_transmission(tvs_nobtov) % tau_levels(jl + iextr,chan_indx)
          end do
          EMI_SFC(chan_indx) = tvs_emissivity(chan_indx,tvs_nobtov)
! *** Gross check on computed BTs ***
          if (BTCALC(chan_indx) < 150.d0) REJFLAG(chan_indx,9) = 1
          if (BTCALC(chan_indx) > 350.d0) REJFLAG(chan_indx,9) = 1
        end do

        KSURF = tvs_profiles(tvs_nobtov) % skin % surftype
!Test pour detecter le isatzen manquant (-1) ou anormal
! (angle negatif ou superieur a 75 degres )
        isatzen= obs_headElem_i(lobsSpaceData,OBS_SZA,INDEX_HEADER)
        if ( isatzen < 9000 .or. &
             isatzen > 16500 ) then
          REJFLAG(:,9) = 1
        end if
!**************************************************************
        CLFR = 0.
        if (lairs) CLFR = obs_headElem_i(lobsSpaceData,OBS_CLF,INDEX_HEADER)

        do JL = 1, NLEVB
          TOEXT(JL) = tvs_profiles(tvs_nobtov)%t(jl + iextr)
        end do

        SUNZA = tvs_profiles(tvs_nobtov) % sunzenangle
        if (liasi) then
          SATAZIM = tvs_profiles(tvs_nobtov) % azangle 
          SATZEN = tvs_profiles(tvs_nobtov) % zenangle
        end if
        ALBEDO =  tvs_surfaceParameters(tvs_nobtov) % albedo
        ICE =  tvs_surfaceParameters(tvs_nobtov) % ice
        LTYPE =  tvs_surfaceParameters(tvs_nobtov) % ltype
        if (LTYPE == 20) KSURF = 2
        PCNT_WAT =  tvs_surfaceParameters(tvs_nobtov) % pcnt_wat
        PCNT_REG =  tvs_surfaceParameters(tvs_nobtov) % pcnt_reg
           
! ** find TOA radiances converted from observed BT's

        ROBS(:) = -1.d0
        
        channels: do JC = 1, nchannels
          chan_indx=channelIndex(JC)
          if ( REJFLAG(chan_indx,9) == 1 ) cycle channels
          t_effective =  tvs_coefs(id) % coef % ff_bco(chan_indx) &
               + tvs_coefs(id) % coef % ff_bcs(chan_indx) * BTOBS(chan_indx)
          ROBS(chan_indx) =  tvs_coefs(id) % coef % planck1(chan_indx) / &
               ( exp( tvs_coefs(id) % coef % planck2(chan_indx) / t_effective ) - 1.d0 )
        end do channels

! ** set height fields to 'height above ground' fields
        do JL = 1, NLEVB
          ZHOEXT(JL,1) = ZHOEXT(JL,1) - ZHT(nlv_T,1)
        end do
        do JL = 1, nlv_T
          ZHT(JL,1) = ZHT(JL,1) - ZHT(nlv_T,1)
        end do

!**********************************************************************************************
!* ///// ---------------------------------------------------- /////
!* ///// DETERMINATION OF THE CLEAR/CLOUDY PROFILES (CLDFLAG) /////
!* ///// ---------------------------------------------------- /////           
        CLDFLAG = 0
        
!* -- REFERENCE FOR WINDOW CHANNEL
        call tvs_getChannelIndexFromChannelNumber(id, IWINDO, IWINDOW(QCID) )
        call tvs_getChannelIndexFromChannelNumber(id, IWINDO_ALT, IWINDOW_ALT(QCID) )
        ICHREF = IWINDO
           
        if ( REJFLAG(IWINDO,9) == 1 ) then
          ICHREF = IWINDO_ALT
          if ( REJFLAG(IWINDO_ALT,9) == 1 ) then
            ICHREF = -1
            CLDFLAG = -1
            REJFLAG(:,9) = 1
            write(*,*) 'WARNING'
            write(*,*) 'WINDOW AND ALTERNATE WINDOW CHANNEL OBSERVATIONS'
            write(*,*) 'HAVE BEEN REJECTED.                             '
            write(*,*) 'ALL '//cinst//' OBSERVATIONS FROM THIS PROFILE REJECTED'
          end if
        end if

!* -- CLOUD TOP BASED ON MATCHING OBSERVED BRIGHTNESS TEMPERATURE 
!* -- AT A REFERENCE SURFACE CHANNEL WITH BACKGROUND TEMPERATURE PROFILE (PTOP_EQ)
!* -- ON GUESS VERTICAL LEVELS.

        LEV_START = 0

!iopt2=1 : calcul de la hauteur en hPa PTOP_MB et du NTOP_MB correspondant
        call CLOUD_HEIGHT (PTOP_MB,NTOP_MB, btobs,cldflag,zt, &
             zht(:,1),zps,zlevmod(:,1),nlv_T,nchn,ichref,lev_start,iopt2)

!iopt1=2 : calcul de la hauteur em metres PTOP_EQ et du NTOP_EQ correspondant
        call CLOUD_HEIGHT (PTOP_EQ,NTOP_EQ, btobs,cldflag,zt, &
             zht(:,1),zps,zlevmod(:,1),nlv_T,nchn,ichref,lev_start,iopt1)

        if (liasi) then
! appel de RTTOV pour calculer les radiances des 3 canaux IR (3b, 4 et 5) de AVHRR 3
           
          call get_avhrr_emiss(emi_sfc(channelIndex(1:nchannels)),tvs_coefs(id) % coef % ff_cwn(channelIndex(1:nchannels)), &
               nchannels,avhrr_surfem1)

          call tovs_rttov_AVHRR_for_IASI(indx,avhrr_surfem1,tvs_satellites(id))
                 
          IOBS = INDX(1)
          call convert_avhrr(sunza, avhrr_bgck(IOBS) )
          call stat_avhrr(avhrr_bgck(IOBS))
          
          LEV_START_AVHRR(:) = 0
          cldflag_avhrr(:) = 0
          do JC=1,nClassAVHRR
            btobs_avhrr(:,JC) = avhrr_bgck(IOBS) % TBMOY(JC,:)
            robs_avhrr(1:NIR,JC) = avhrr_bgck(IOBS) % RADMOY(JC,NVIS + 1:NIR + NVIS)
            RCAL_CLR_AVHRR(:) = avhrr_bgck(IOBS) % RADCLEARCALC(:)
            EMI_SFC_AVHRR(:) = avhrr_bgck(IOBS) % EMISS(:)
            SFCTAU_AVHRR(:) = avhrr_bgck(IOBS) % TRANSMSURF(:)
            
            do JL=1,NLEVB
              RCLD_AVHRR(:,JL) = avhrr_bgck(IOBS) % RADOVCALC(JL + iextr - 1,:)
            end do
           
            if (btobs_avhrr(2,JC) > 100.d0 ) then
              ichref_avhrr(JC) = 2
            else if (btobs_avhrr(3,JC) > 100.d0 ) then
              ichref_avhrr(JC) = 3
            else
              ichref_avhrr(JC) = -1
              cldflag_avhrr(JC) = -1
            end if
            
            call CLOUD_HEIGHT (PTOP_EQ_AVHRR(JC),NTOP_EQ_AVHRR(JC), btobs_avhrr(:,JC),cldflag_avhrr(JC),zt, &
                 zht(:,1),zps,zvlev,nlv_T,NIR,ichref_avhrr(JC),lev_start_avhrr(JC),iopt1)
          end do
          
        end if

!* -- CLEAR/CLOUDY PROFILE DETECTION USING THE GARAND & NADON ALGORITHM

        call GARAND1998NADON (CLDFLAG, btobs,ztg,zt, &
             zht(:,1),nlv_T,nchn,ptop_eq,ntop_eq,ichref)

        if (liasi) then
          do JC=1,nClassAVHRR
            call GARAND1998NADON (CLDFLAG_AVHRR(jC), btobs_avhrr(:,JC),ztg,zt, &
                 zht(:,1),nlv_T,NIR,ptop_eq_avhrr(JC),ntop_eq_avhrr(JC),ichref_avhrr(JC))
          end do
        end if
        
!* -- FURTHER TESTS TO REMOVE POTENTIAL CLOUDY PROFILES
! *** TEST # A ***
! *** In daytime, set cloudy if cloud fraction over 5% ***
        CFSUB = -1.d0
        if (lairs) then
          if ( CLDFLAG == 0 .and. CLFR > 5.d0 .and. SUNZA < 90.d0 ) then
            CLDFLAG = 1
            CFSUB = 0.01d0 * CLFR !conversion % -> 0-1
          end if
        end if
! *** TEST # B ***
! *** Set cloudy if temperature difference between guess (ZTG)     ***
! *** and estimated true (ZTS) skin temperatures is over threshold ***

        call ESTIM_TS(ZTS, ztg, emi_sfc, rcal_clr, robs, &
             sfctau, cldflag, ichref, nchn, tvs_coefs(id) )

        if ( CLDFLAG == 0 .and. KSURF == 1 &
             .and. abs(ZTS-ZTG) > DTW ) CLDFLAG = 1 

        if ( CLDFLAG == 0 .and. KSURF /= 1 &
             .and. abs(ZTS-ZTG) > DTL ) CLDFLAG = 1

        if (liasi) then

          do JC=1,nClassAVHRR
            call ESTIM_TS(ZTS_AVHRR(JC), ztg,emi_sfc_avhrr,rcal_clr_avhrr,robs_avhrr(:,JC), &
                 sfctau_avhrr,CLDFLAG_AVHRR(JC),ichref_avhrr(JC),NIR, coefs_avhrr)
          end do

          do JC=1,nClassAVHRR
            if ( CLDFLAG_AVHRR(JC) == 0 .and. KSURF == 1 &
                 .and. abs(ZTS_AVHRR(JC)-ZTG) > DTW ) CLDFLAG_AVHRR(JC) = 1
              
            if ( CLDFLAG_AVHRR(JC) == 0 .and. KSURF /= 1 &
                 .and. abs(ZTS_AVHRR(JC)-ZTG) > DTL ) CLDFLAG_AVHRR(JC) = 1
              
          end do

!criteres AVHRR utilisant les canaux visibles (de jour seulement)
          if (sunza < sunzenmax) then 
            ANISOT = 1.d0
            deltaphi = abs(SATAZIM - SUNAZIM )
           
            if (deltaphi > 180.d0) deltaphi = 360.d0 - deltaphi
            
            if (ALBEDO < 0.17d0) then               
              call VISOCN(sunza,satzen,deltaphi,ANISOT,ZLAMB,ZCLOUD,IER)
              SEUIL_ALBED = 10.d0 * max(1.d0,ANISOT) 
            else
              SEUIL_ALBED = 100.d0 * ALBEDO + 10.d0
            end if
              
            if (ANISOT < 1.5d0) then !to avoid sun glint
              SCOS = cos ( sunza * MPC_DEGREES_PER_RADIAN_R8 )
              call  cor_albedo ( DEL, SCOS )
              SEUIL_ALBED = SEUIL_ALBED * DEL
              do JC=1,nClassAVHRR
                if (avhrr_bgck(IOBS)%ALBEDMOY(JC,1) > SEUIL_ALBED(1) ) then
                  CLDFLAG_AVHRR(JC) = 1
                end if
                  !static AVHRR thresholds v3
                do JL=1,NVIS
                  if (avhrr_bgck(IOBS)%ALBEDMOY(JC,JL) > seuilalb_static(JL,KSURF) ) then
                    CLDFLAG_AVHRR(JC) = 1
                  end if
                end do
              end do
             
            end if
          end if

!Calcul de la pseudo fraction nuageuse AVHRR

          CFRAC_AVHRR = 0.d0
          do JC=1,nClassAVHRR
            if (CLDFLAG_AVHRR(JC) == 1) CFRAC_AVHRR = CFRAC_AVHRR + avhrr_bgck(IOBS) % CFRAC(JC)
          end do

          CFSUB = -1.0d0
          if ( CLDFLAG == 0 .and. CFRAC_AVHRR > 5.d0 ) then
            CLDFLAG = 1
            CFSUB = 0.01d0 * min(CFRAC_AVHRR,100.d0) !conversion % -> 0-1 avec seuil car parfois CFRAC_AVHRR=101
          end if

!AVHRR Homogeneity criteria
          if (CLDFLAG == 0) then
            IJOUR = 1
            if (SUNZA < 90.d0) IJOUR=2
            ! 1 NUIT
            ! 2 JOUR
            if (IJOUR == 2) then
              do JC=1,NVIS
                if (avhrr_bgck(IOBS)%ALBSTD_PIXELIASI(JC) > seuilalb_homog(JC,KSURF) ) CLDFLAG = 1
              end do
            end if
            do JC=NVIS+1,NVIS+NIR
              if (avhrr_bgck(IOBS)%TBSTD_PIXELIASI(JC) > seuilbt_homog(JC,KSURF,IJOUR)) CLDFLAG = 1
            end do
          end if
        end if

        GNCLDFLAG = CLDFLAG

!* ///// ------------------------------------------------------- /////
!* ///// DETERMINATION OF THE ASSIMILABLE OBSERVATIONS (REJFLAG) /////
!* ///// ------------------------------------------------------- /////


!* -- FIRST TESTS TO REJECT OBSERVATIONS


! *** TEST # 1 ***
! *** Do not assimilate where cloudy ***

        if ( CLDFLAG == 1 ) then
          REJFLAG(:,11) = 1
          REJFLAG(:,23) = 1
        end if

! *** TEST # 2 ***
! *** Gross check on valid BTs ***

!     already done


!* -- CLOUD TOP BASED ON MATCHING 
!* -- OBSERVED BRIGHTNESS TEMPERATURE WITH BACKGROUND TEMPERATURE PROFILES (PTOP_BT)
!* -- OR COMPUTED OBSERVED RADIANCES WITH BACKGROUND RADIANCE PROFILES (PTOP_RD)
!* -- ON RTTOV VERTICAL LEVELS

        LEV_START = 0

        do JCH = 1, NCH_HE
          call tvs_getChannelIndexFromChannelNumber(id,ILIST_HE(JCH),ILIST1(QCID,JCH))
        end do

        call CLOUD_TOP ( PTOP_BT,PTOP_RD,NTOP_BT,NTOP_RD, &
             btobs,toext,zhoext(:,1),rcal_clr,zps,robs,rcld,zvlev,nlevb, &
             nchn,cldflag,rejflag,lev_start,iopt2,ihgt,ichref,nch_he,ilist_he)

        if (liasi) then
          LEV_START_AVHRR(:) = 0
          
          do JC=1,nClassAVHRR
            call CLOUD_TOP_AVHRR ( PTOP_BT_AVHRR(:,JC),PTOP_RD_AVHRR(:,JC),NTOP_BT_AVHRR(:,JC),NTOP_RD_AVHRR(:,JC), &
                 btobs_avhrr(:,JC),toext,zhoext(:,1),rcal_clr_avhrr,zps,robs_avhrr(:,JC),rcld_avhrr,zvlev,nlevb, &
                 NIR,cldflag_avhrr(jc),lev_start_avhrr(JC),iopt2,ihgt,nn,ilist_avhrr)
          end do
        end if

!* -- REFERENCE CHANNEL FOR CO2-SLICING

        do JCH = 1, NCO2
          call tvs_getChannelIndexFromChannelNumber(id, ILIST_CO2(JCH), ILIST2(QCID,JCH)  )
          call tvs_getChannelIndexFromChannelNumber(id, ILIST_CO2_PAIR(JCH), ILIST2_PAIR(QCID,JCH)  )
        end do

        cpt = 0
        do JCH=1,NCO2
          if ( REJFLAG(ILIST_CO2(JCH),9) == 1 .or. &
               REJFLAG(ILIST_CO2_PAIR(JCH),9) == 1 ) cpt = cpt + 1
        end do
         
        if (cpt == nco2) then
          CLDFLAG = -1
          REJFLAG(:,9) = 1
          write(*,*) 'WARNING'
          write(*,*) 'CO2 REFERENCE AND ALTERNATE CHANNEL OBSERVATIONS'
          write(*,*) 'HAVE BEEN REJECTED.                             '
          write(*,*) 'ALL '//CINST//' OBSERVATIONS FROM THIS PROFILE REJECTED'
        end if

!* -- EQUIVALENT HEIGHT OF SELECTED WINDOW CHANNEL
        call tvs_getChannelIndexFromChannelNumber(id,chan_indx,ILIST1(QCID,2))
        HEFF = PTOP_RD( chan_indx )

              
        if (ICHREF == IWINDO_ALT) then
          call tvs_getChannelIndexFromChannelNumber(id,chan_indx,ILIST1(QCID,3))
          HEFF = PTOP_RD( chan_indx )
        end if
!* -- CLOUD TOP BASED ON CO2 SLICING 

        
        LEV_START = max( min(LEV_START,CO2MAX(1)), CO2MIN(1) )

        call CO2_SLICING ( PTOP_CO2,NTOP_CO2,FCLOUD_CO2, &
             rcal_clr,rcld,robs,zps,zvlev,nlevb,nchn,cldflag,rejflag, &
             lev_start,ichref,ilist_co2,ilist_co2_pair)

!* -- FIND CONSENSUS CLOUD TOP AND FRACTION
 
        call SELTOP ( ETOP,VTOP,ECF,VCF,NGOOD, heff,ptop_co2,fcloud_co2, &
             CFSUB,PTOP_MB,zps,cldflag,gncldflag )

        if (liasi) then
! Correction pour les nuages trop bas:
! en principe Pco2 < Heff.
! on cherche les cas pathologiques avec Pco2>Min(Heff(AVHRR))
          minpavhrr(2:3) = 12200
          ILOC(2:3) = -1      ! pour eviter les catastrophes...
          do JC=1,nClassAVHRR
            if (avhrr_bgck(IOBS)%CFRAC(JC) > 0.d0) then
              if (PTOP_RD_AVHRR(2,JC) < minpavhrr(2)) then
                ILOC(2) = JC
                minpavhrr(2) = PTOP_RD_AVHRR(2,JC)
              end if
              if (PTOP_RD_AVHRR(3,JC) < minpavhrr(3)) then
                ILOC(3) = JC
                minpavhrr(3) = PTOP_RD_AVHRR(3,JC)
              end if
            end if
          end do
          if ( ILOC(2) /= -1 .and. ILOC(3) /= -1) then ! pour eviter les catastrophes...
            ! on se limite aux cas "surs" ou les deux hauteurs effectives sont > a Pco2
            ! et ou un accord raisonnable existe entre les deux hauteurs effectives
            if ( ILOC(2) == ILOC(3) .and. &
                 minpavhrr(2) < ETOP .and. &
                 minpavhrr(3) < ETOP .and. &
                 abs(minpavhrr(2)- minpavhrr(3)) < 25.d0 .and. &
                 CLDFLAG_AVHRR(ILOC(2)) /= -1 .and. CLDFLAG_AVHRR(ILOC(3)) /= -1) then
              
              if (ECF == 0.d0 .and. CLDFLAG == 1) then
                ! cas predetermine nuageux mais ramene a clair 
                ECF = 0.01d0 * min(100.d0,CFRAC_AVHRR)
                ! cette ligne peut generer des fractions nuageuses inferieures a 20 %.
                ETOP = 0.5d0 * (minpavhrr(2) + minpavhrr(3))
              end if

              if (ECF > 0.d0 .and. CLDFLAG == 1) then
                !cas predetermine nuageux pas ramene clair (==normal)
                ETOP = 0.5d0 * ( minpavhrr(2) + minpavhrr(3))
              end if

              if (CLDFLAG == 0) then
                !cas predetermine clair ... que faire
                CLDFLAG = 1
                ETOP = 0.5d0 * (minpavhrr(2) + minpavhrr(3))
                ECF = 0.01d0 * min(100.d0,CFRAC_AVHRR)
              end if
            end if
          end if
        end if

        !* -- FIND MINIMUM LEVEL OF SENSITIVITY FOR CHANNEL ASSIMILATION NOT SENSIBLE TO CLOUDS        
        call MIN_PRES_new (MAXWF, MINP,PMIN,DTAUDP1, zps,transm,zvlev,cldflag,nlevb,nchn,imodtop )
!* -- ASSIMILATION OF OBSERVATIONS WHEN CLOUDY PROFILES

! *** TEST # 3 ***
! *** Assimilation above clouds (refinement of test 1)             ***
! *** Set security margin to 2x the std on height from CO2-slicing *** 

        TAMPON = max(50.d0, 2.d0*VTOP)                                                          

        do JC = 1, nchn        
          if ( REJFLAG(JC,11) == 1 .and. REJFLAG(JC,23) == 1 .and. ETOP - TAMPON > PMIN(JC) ) then
            REJFLAG(JC,11) = 0
            REJFLAG(JC,23) = 0
          end if
        end do

!     LOOK AT THE FATE OF THE OBSERVATIONS
        FATE(:) = sum(REJFLAG(:,:), DIM=2)            

!     FURTHER REASONS TO REJECT OBSERVATIONS

        call  tvs_getChannelIndexFromChannelNumber(id,ILIST_SUN,ICHN_SUN(QCID))

        do JC = 1, nchn

          if ( FATE(JC) == 0 ) then

! *** TEST # 4 ***
! *** Background check, do not assimilate if O-P > 3sigma ***

            if ( abs(BTOBS(JC) - BTCALC(JC)) > 3.d0 * BTOBSERR(JC) ) then
              REJFLAG(JC,9) = 1
              REJFLAG(JC,16) = 1
            end if

! *** TEST # 5 ***
! *** Do not assimilate shortwave channels during the day ***

            if ( JC >= ILIST_SUN .and. SUNZA < NIGHT_ANG ) then
              REJFLAG(JC,11) = 1
              REJFLAG(JC,7)  = 1
            end if

! *** TEST # 6 ***
! *** Do not assimilate surface channels over land ***

            if ( MINP(JC) == NLEVB .or. ZPS-PMIN(JC) < 100.d0 ) then
              if ( KSURF == 0 ) then
                REJFLAG(JC,11) = 1    !!! comment this line if assimilation under conditions
                REJFLAG(JC,19) = 1    !!! comment this line if assimilation under conditions
                if ( PCNT_WAT > 0.01d0 .or. PCNT_REG > 0.1d0 .or. EMI_SFC(JC) < 0.97d0 ) then
                  REJFLAG(JC,11) = 1
                  REJFLAG(JC,19) = 1
                end if

! *** TEST # 7 ***
! *** Do not assimilate surface channels over water under conditions ***

              else if ( KSURF == 1 ) then
                if ( PCNT_WAT < 0.99d0 .or. PCNT_REG < 0.97d0 .or. &
                     ICE > 0.001d0 .or. ALBEDO >= 0.17d0 .or. EMI_SFC(JC) < 0.9d0 ) then
                  REJFLAG(JC,11) = 1   
                  REJFLAG(JC,19) = 1   
                end if
                
! *** TEST # 8 ***
! *** Do not assimilate surface channels over sea ice ***
                          
              else if ( KSURF == 2 ) then
                REJFLAG(JC,11) = 1
                REJFLAG(JC,19) = 1   
              end if
            end if
            
          end if

! *** TEST # 9 ***
! *** Do not assimilate if jacobian has a significant contribution over model top ***

! Condition valid if model top at 10mb or lower only
          if ( nint(ptop_T) >= 1000 ) then
            if ( REJFLAG(JC,9) /= 1 .and. DTAUDP1(JC) > 0.50d0 ) then
              REJFLAG(JC,11) = 1
              REJFLAG(JC,21) = 1
            end if
          end if
        
! Condition valid if model top at 10mb or lower only
          if ( nint(ptop_T) >= 1000 ) then
            if ( REJFLAG(JC,9) /= 1 .and. TRANSM(JC,1) < 0.99d0 ) then
              REJFLAG(JC,11) = 1
              REJFLAG(JC,21) = 1 
            end if
          end if

! Condition valid if model top is higher than 10 mb
          if ( nint(ptop_T) < 1000 ) then
            if ( REJFLAG(JC,9) /= 1 .and. TRANSM(JC,1) < 0.95d0 ) then
              REJFLAG(JC,11) = 1
              REJFLAG(JC,21) = 1 
            end if
          end if

        end do

        nchannels =0 
        do INDEX_BODY= IDATA, IDATEND
          if ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY)==1 ) then
            nchannels =  nchannels + 1
            if (btest(obs_bodyElem_i(lobsSpaceData,OBS_FLG,INDEX_BODY),8)) REJFLAG(channelIndex(nchannels),8) = 1
          end if
        end do

!* -- FOR EACH PROFILE, ARE ALL NON-BLACKLISTED CHANNELS ASSIMILATED

        ASSIM_ALL = .true.
        FATE(:) = sum(REJFLAG(:,:),DIM=2)            
        
        chn: do JC = 1, nchn
          if ( REJFLAG(JC,8) == 0 ) then
            if ( FATE(JC) /= 0 ) then
              ASSIM_ALL = .false.
              exit chn
            end if
          end if
        end do chn

        if  (.not.ASSIM_ALL) then
          call obs_headSet_i(lobsSpaceData, OBS_ST1, index_header,ibset(obs_headElem_i(lobsSpaceData,OBS_ST1,INDEX_HEADER),6) )
        end if
!* -- ADDITION OF BACKGROUND CHECK PARAMETERS TO BURP FILE
!* ------------------------------------------------

        call obs_headSet_r(lobsSpaceData, OBS_ETOP, index_header, ETOP )
        call obs_headSet_r(lobsSpaceData, OBS_VTOP, index_header, VTOP )
        call obs_headSet_r(lobsSpaceData, OBS_ECF,  index_header, 100.d0 * ECF )
        call obs_headSet_r(lobsSpaceData, OBS_VCF,  index_header, 100.d0 * VCF )
        call obs_headSet_r(lobsSpaceData, OBS_HE,   index_header, HEFF )
        call obs_headSet_r(lobsSpaceData, OBS_ZTSR, index_header, ZTS )
        call obs_headSet_i(lobsSpaceData, OBS_NCO2, index_header, NGOOD)
        call obs_headSet_r(lobsSpaceData, OBS_ZTM,  index_header, ZT(nlv_T) )
        call obs_headSet_r(lobsSpaceData, OBS_ZTGM, index_header, ZTG )
        call obs_headSet_r(lobsSpaceData, OBS_ZLQM, index_header, exp(ZLQS) )
        call obs_headSet_r(lobsSpaceData, OBS_ZPS,  index_header, 100.d0 * ZPS )
        call obs_headSet_i(lobsSpaceData, OBS_STYP, index_header, KSURF )

        do INDEX_BODY= IDATA, IDATEND
          ICHN = nint(obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY))
          ICHN = max(0, min(ICHN, tvs_maxChannelNumber + 1))
          call tvs_getChannelIndexFromChannelNumber(id,chan_indx,ICHN)
          call obs_bodySet_r(lobsSpaceData,OBS_SEM,INDEX_BODY,EMI_SFC(chan_indx))
          do NFLG = 0, BITFLAG
            if ( REJFLAG(chan_indx,NFLG) == 1 ) &
                 call obs_bodySet_i(lobsSpaceData,OBS_FLG,INDEX_BODY,IBSET(obs_bodyElem_i(lobsSpaceData,OBS_FLG,INDEX_BODY),NFLG))
          end do
        end do
          
      end if

    end do HEADER_2

    deallocate ( channelIndex, stat= alloc_status(1))
    deallocate ( ZHT,       stat= alloc_status(2))
    deallocate ( ZT,        stat= alloc_status(3))
    deallocate ( ZLEVMOD,   stat= alloc_status(4))
    deallocate ( ZVLEV,     stat= alloc_status(5))
    deallocate ( maxwf,     stat= alloc_status(6))
    if (liasi) deallocate ( RCLD_AVHRR , stat= alloc_status(7))
    deallocate ( FATE,      stat= alloc_status(8))
    deallocate ( DTAUDP1,   stat= alloc_status(9))
    deallocate ( PMIN,      stat= alloc_status(10))
    deallocate ( MINP,      stat= alloc_status(11))
    deallocate ( PTOP_RD,   stat= alloc_status(12))
    deallocate ( PTOP_BT,   stat= alloc_status(13))
    deallocate ( NTOP_RD,   stat= alloc_status(14))
    deallocate ( NTOP_BT,   stat= alloc_status(15))
    deallocate ( REJFLAG,   stat= alloc_status(16))
    deallocate ( ROBS,      stat= alloc_status(17))
    deallocate ( ZHOEXT,    stat= alloc_status(18))
    deallocate ( TOEXT,     stat= alloc_status(19))
    deallocate ( EMI_SFC,   stat= alloc_status(20))
    deallocate ( TRANSM,    stat= alloc_status(21))
    deallocate ( RCLD,      stat= alloc_status(22))
    deallocate ( SFCTAU,    stat= alloc_status(23))
    deallocate ( RCAL_CLR,  stat= alloc_status(24))
    deallocate ( BTCALC,    stat= alloc_status(25))
    deallocate ( BTOBS,     stat= alloc_status(26))
    deallocate ( BTOBSERR,  stat= alloc_status(27))
    deallocate ( XPRES,     stat= alloc_status(28))
    call utl_checkAllocationStatus(alloc_status, " irbg_doQualityControl", .false.)
        
  end subroutine irbg_doQualityControl

  subroutine convert_avhrr(sunzen,avhrr)
! conversion des radiance IR en temperatures de brillance
! et des radiances visibles en "albedo"
  
    implicit none
    real(8) ,intent(in) :: sunzen
    type (avhrr_bgck_iasi) ,intent(inout) :: avhrr

    integer :: ICL
    real (8) :: tb(NIR),dtbsdrad(NIR)
    real (8) :: FREQ(NIR),OFFSET(NIR),SLOPE(NIR)


    freq = coefs_avhrr%coef%ff_cwn (:)
    offset = coefs_avhrr%coef%ff_bco(:)
    slope = coefs_avhrr%coef%ff_bcs(:)

    do ICL=1,nClassAVHRR
      call calcbt(avhrr % radmoy(ICL,4:6), tb, dtbsdrad,freq,offset,slope)
      avhrr % tbmoy(ICL,4:6) = tb(1:3)
      avhrr % tbstd(ICL,4:6) = avhrr % radstd(ICL,4:6) * dtbsdrad(1:3)
      call calcreflect(avhrr % radmoy(ICL,1:3) ,sunzen,avhrr % ALBEDMOY(ICL,1:3) )
      call calcreflect(avhrr % radstd(ICL,1:3) ,sunzen,avhrr % ALBEDSTD(ICL,1:3) )
    end do

  end subroutine convert_avhrr

  subroutine calcreflect(rad,sunzen,reflect)
    implicit none

    real (8) , intent(in) :: rad(nvis)
    real (8) , intent(in) :: sunzen
    real (8) , intent(out):: reflect(nvis) ! reflectivite en %
    !************
    real (8) :: SOLAR_FILTERED_IRRADIANCE(nvis)
    data SOLAR_FILTERED_IRRADIANCE /139.873215d0,232.919556d0,14.016470d0/
!# equivalent widths, integrated solar irradiance,  effective central wavelength
!0.084877,139.873215,0.632815
!0.229421,232.919556,0.841679
!0.056998,14.016470,1.606119
    ! pour la definition de l'albedo voir http://calval.cr.usgs.gov/PDF/Rao.CRN_IJRS.24.9.2003_Chander.pdf
    real (8) :: RADB ! radiance en W/m2/str
    integer :: i
    !**************************************************************

    do i = 1, nvis
      if (rad(i) >= 0.0d0 ) then
        radb = rad(i) / 1000.0d0
        reflect(i) = (MPC_PI_R8 * radb) / SOLAR_FILTERED_IRRADIANCE(I)
        if (sunzen < 90.0d0 ) reflect(i) = reflect(i) / cos(sunzen * MPC_RADIANS_PER_DEGREE_R8)
      else
        reflect(i) = -1
      end if
    end do
  
  end subroutine calcreflect

  subroutine calcbt(rad,tb,dtbsdrad,freq,offset,slope)
    implicit none
    integer,parameter  :: nchan=3
    real(8) ,parameter :: c1= 1.19106590D-05   ! first planck constant
    real(8) ,parameter :: c2= 1.438833d0     ! second planck constant 
    real (8) , intent(in) :: rad(nchan), freq(nchan), offset(nchan), slope(nchan)
    real (8) , intent(out):: tb(nchan), dtbsdrad(nchan)
    !************
    integer :: i
    real (8) ::  radtotal,tstore,planck1,planck2

    do i = 1, nchan
      if (rad(i) > 1.d-20) then
        planck2 = c2 * freq(I)
        planck1 = c1 * ( freq(I) ** 3 ) 
        tstore = planck2 / log( 1.0d0 + planck1 / rad(i) )
        tb(i) = ( tstore - offset(i) ) / slope(i)
        
        radtotal = rad(i)
        
        dtbsdrad(i) = planck1 * tstore ** 2 / ( planck2 * radtotal * ( radtotal + planck1 ) )
        
        dtbsdrad(i) = dtbsdrad(i) / slope(i)
        
      else
        tb(i) = 0.d0
        dtbsdrad(i) = 0.d0
      end if
      
    end do

  end subroutine calcbt

  subroutine stat_avhrr(avhrr)
    ! calcul de statistiques
    ! sur l'information sous-pixel AVHRR
    implicit none
    type (avhrr_bgck_iasi) ,intent(inout) :: avhrr
    integer :: ICL,ICH
    real (8) :: SUMFRAC(NVIS+NIR),TBMIN(NVIS+1:NVIS+NIR),TBMAX(NVIS+1:NVIS+NIR),SUMTB(NVIS+1:NVIS+NIR),SUMTB2(NVIS+1:NVIS+NIR)
    real (8) :: SUMALB(1:NVIS),SUMALB2(1:NVIS)
    !******************************************
    
    SUMFRAC(:) = 0.d0
    SUMTB(:) = 0.d0
    SUMTB2(:) = 0.d0
    SUMALB(:) = 0.d0
    SUMALB2(:) = 0.d0
    
    do ICL=1,nClassAVHRR
      if (avhrr%CFRAC(ICL) > 0.d0 ) then
        do ICH=1,NVIS
          if (avhrr%ALBEDMOY(ICL,ICH) >= 0.d0 ) then
            SUMFRAC(ICH) = SUMFRAC(ICH) + avhrr%CFRAC(ICL)
            SUMALB(ICH) = SUMALB(ICH) + avhrr%CFRAC(ICL) * avhrr%ALBEDMOY(ICL,ICH)
            SUMALB2(ICH) = SUMALB2(ICH) + avhrr%CFRAC(ICL) * ( avhrr%ALBEDMOY(ICL,ICH)**2 + avhrr%ALBEDSTD(ICL,ICH)**2)
          end if
        end do
        do ICH=1+NVIS,NVIS+NIR
          if (avhrr%TBMOY(ICL,ICH) > 0.d0 ) then
            SUMFRAC(ICH) = SUMFRAC(ICH) + avhrr%CFRAC(ICL)
            SUMTB(ICH) = SUMTB(ICH) + avhrr%CFRAC(ICL) * avhrr%TBMOY(ICL,ICH)
            SUMTB2(ICH) = SUMTB2(ICH) + avhrr%CFRAC(ICL) * (avhrr%TBMOY(ICL,ICH)**2 + avhrr%TBSTD(ICL,ICH)**2 )
          end if
        end do
      end if
    end do
      
    do ICH=1,NVIS
      if (SUMFRAC(ICH) > 0.d0 ) then
        SUMALB(ICH) = SUMALB(ICH) / SUMFRAC(ICH)
        SUMALB2(ICH) = SUMALB2(ICH)/SUMFRAC(ICH) - SUMALB(ICH)**2
        if (SUMALB2(ICH) > 0.d0) then
          SUMALB2(ICH) = sqrt( SUMALB2(ICH) )
        else
          SUMALB2(ICH) = 0.d0
        end if
      end if
    end do
      
    do ICH=NVIS+1,NVIS+NIR
      if (SUMFRAC(ICH) > 0.d0 ) then
        SUMTB(ICH) = SUMTB(ICH) / SUMFRAC(ICH)
        SUMTB2(ICH) = SUMTB2(ICH)/SUMFRAC(ICH) - SUMTB(ICH)**2
        if (SUMTB2(ICH) > 0.d0) then
          SUMTB2(ICH)= sqrt ( SUMTB2(ICH) )
        else
          SUMTB2(ICH) = 0.d0
        end if
      end if
    end do
      
    avhrr % TBSTD_PIXELIASI = SUMTB2
    avhrr % ALBSTD_PIXELIASI = SUMALB2
      
  end subroutine stat_avhrr

  subroutine CO2_SLICING ( PTOP,NTOP,FCLOUD,    &
       rcal,rcld,robs,ps,plev,nlev,nchn,cldflag,rejflag, &
       lev_start,ichref,ilist,ilist_pair)
!
!**ID CO2_SLICING -- CLOUD TOP HEIGHT COMPUTATION
!
!       AUTHOR:   L. GARAND               July 2004
!                 A. BEAULNE (CMDA/SMC)  March 2006  (ADAPT TO 3DVAR)                 
!
!       REVISION: 001 O. Pancrati various improvements
!
!       OBJECT:   CLOUD TOP FROM CO2 SLICING AND CLOUD FRACTION ESTIMATE
!
!       ARGUMENTS:
!          INPUT:
!            -RCAL(NCHN)      : COMPUTED CLEAR RADIANCES (MW/M2/SR/CM-1)
!            -RCLD(NCHN,NLEV) : COMPUTED CLOUD RADIANCES FROM EACH LEVEL (")
!            -ROBS(NCHN)      : COMPUTED OBSERVED RADIANCES (")
!            -PS             : SURFACE PRESSURE (HPA)
!            -PLEV(NLEV)           : PRESSURE LEVELS (HPA)
!            -NLEV                 : NUMBER OF VERTICAL LEVELS
!            -NCHN                 : NUMBER OF CHANNELS
!            -CLDFLAG        : (0) CLEAR, (1) CLOUDY, (-1) UNDEFINED PROFILE
!            -REJFLAG(NCHN,0:BITFLAG) : FLAGS FOR REJECTED OBSERVATIONS
!            -BITFLAG              : HIGHEST FLAG IN POST FILES (VALUE OF N IN 2^N)
!            -ICHREF         : WINDOW CHANNEL TO PREDETERMINE CLEAR
!            -NCO2                 : NUMBER OF CHANNELS TO GET ESTIMATES IN
!                                     COMBINATION WITH ICHREF_CO2 (NOT INCLUDED)
!            -ILIST(NCO2)          : LIST OF THE CHANNEL NUMBERS, ICHREF_CO2 NOT INCLUDED
!                                     (SUBSET VALUES)
!
!          INPUT/OUTPUT:
!            -LEV_START      : LEVEL TO START ITERATION (IDEALLY TROPOPAUSE)
!
!          OUTPUT:
!            -PTOP(NCO2)      : CLOUD TOP (HPA)
!            -FCLOUD(NCO2)    : CLOUD FRACTION
!            -NTOP(NCO2)      : NEAREST PRESSURE LEVEL CORRESPONDING TO PTOP
!                                     (PTOP <= PS)
!
    implicit none
    integer ,intent (in) :: NLEV,NCHN
    real(8) ,intent (in) :: RCAL(NCHN),RCLD(NCHN,NLEV),ROBS(NCHN)
    real(8) ,intent (in) :: PLEV(NLEV),PS
    integer ,intent (in) :: ICHREF,CLDFLAG,REJFLAG(NCHN,0:BITFLAG)
    integer ,intent (in) :: ILIST(NCO2),ILIST_PAIR(NCO2)
    integer ,intent (inout) :: LEV_START
    real(8) ,intent (out) :: PTOP(NCO2),FCLOUD(NCO2)
    integer ,intent (out) :: NTOP(NCO2)
    !*********************************************************************************
    integer     :: J,JCH,JC,JPMAX,JMAX
    integer     :: SUMREJ
    real(8)     :: EPS
    real(8)     :: FC(NCHN,NLEV),RAPG,RADP
    real(8)     :: DRAP(NCO2,NLEV),A_DRAP(NLEV)
    real(8)     :: VAL,VAL1,VAL2,VAL3,FCINT
    real(8)     :: EMI_RATIO
    integer     :: JC_PAIR
    integer     :: ITER,NITER
      
    EPS = 1.D-12
    
    PTOP(:) = -1.d0
    NTOP(:) = -1
    FCLOUD(:) = -1.d0


    !**     profile not assimilated if data from 2 windows channels bad
    !**     and/or if data from 2 reference co2 channels bad

    if ( CLDFLAG == -1 ) return

    !**     define closest level jpmax to surface pressure ps

    JPMAX = NLEV
    
    do J = LEV_START, NLEV
      if ( PLEV(J) > PS ) then
        JPMAX = J
        exit
      end if
    end do
    
    !**     define jmax as last level for co2-slicing calculations
  
    JMAX = JPMAX - 1
    
    !**     predetermined clear window channel, all nco2 estimates clear

    SUMREJ = sum(REJFLAG(ICHREF,:))

    if ( SUMREJ == 0 ) then
      PTOP(:) = PS
      NTOP(:) = JPMAX
      FCLOUD(:) = 0.d0
      return
    end if

    channels: do JCH = 1, NCO2
      
      JC = ILIST(JCH)
      JC_PAIR = ILIST_PAIR(JCH)
      FC(JC_PAIR,:) = RCAL(JC_PAIR) - RCLD(JC_PAIR,:)
      NITER = 1
      if ( JCH > 13) NITER = 2 
     
      iteration: do ITER = 1, NITER
        DRAP(JCH,:)   = 9999.d0
        NTOP(JCH) = -1
        !-------------------------------------------------------------------------------
        !         calcul EMI_RATIO
        if (JCH > 13) then       
          if ( ITER == 1 ) then
            EMI_RATIO = 1.0376d0
          else
            EMI_RATIO = 1.09961d0 - 0.09082d0 * FCLOUD(JCH)
          end if
        else
          EMI_RATIO = 1.0d0
        end if
!-------------------------------------------------------------------------------

        FC(JC,:) = RCAL(JC) - RCLD(JC,:)

!**       gross check failure

        if ( REJFLAG(JC,9) == 1 ) cycle channels
        if ( REJFLAG(JC_PAIR,9) == 1 ) cycle channels

        if ( abs( RCAL(JC_PAIR) - ROBS(JC_PAIR) ) > EPS ) then
          RAPG = (RCAL(JC) - ROBS(JC)) / (RCAL(JC_PAIR) - ROBS(JC_PAIR))
        else
          RAPG = 0.0d0
        end if

        do J = LEV_START, JPMAX
          if ( FC(JC,J) > 0.d0 .and. FC(JC_PAIR,J) > 0.d0 )  &
               DRAP(JCH,J) = RAPG - (FC(JC,J) / FC(JC_PAIR,J)) * EMI_RATIO
        end do

        A_DRAP(:) = abs( DRAP(JCH,:) )

        levels: do J = LEV_START + 1, JMAX

          !**         do not allow fc negative (i.e. drap(jch,j) = 9999.)

          if ( DRAP(JCH,J) > 9000.d0 .and. &
               A_DRAP(J-1) < EPS .and. &
               A_DRAP(J+1) < EPS ) cycle channels

          VAL = DRAP(JCH,J) / DRAP(JCH,J - 1)

!**         find first, hopefully unique, zero crossing

          if ( VAL < 0.d0 ) then

!**         conditions near zero crossing of isolated minimum need monotonically
!**         decreasing drap from j-3 to j-1 as well increasing from j to j+1

            VAL1 = DRAP(JCH,J - 2) / DRAP(JCH,J - 1)
            VAL2 = DRAP(JCH,J - 3) / DRAP(JCH,J - 1)
            VAL3 = DRAP(JCH,J) / DRAP(JCH,J + 1)

            if ( VAL1 > 0.d0 .and.  & 
                 VAL2 > 0.d0 .and.  & 
                 VAL3 > 0.d0 .and.  &
                 A_DRAP(J-2) > A_DRAP(J-1) .and.  &
                 A_DRAP(J-3) > A_DRAP(J-2) .and.  &
                 A_DRAP(J)   < 9000.d0     .and.  &
                 A_DRAP(J+1) > A_DRAP(J) )        &
                 then
              PTOP(JCH) = PLEV(J)
              NTOP(JCH) = J
            end if
            
            exit levels
                      
          end if
              
        end do levels

        J = NTOP(JCH)

!**       special cases of no determination

        if ( J < 1) then
          PTOP(JCH)   = -1.d0
          NTOP(JCH)   = -1
          FCLOUD(JCH) = -1.d0
          cycle channels
        end if
      
        if ( J <= LEV_START .or. DRAP(JCH,J) > 9000.d0 ) then
          !if ( ITER == 1) then
          PTOP(JCH) = -1.d0
          NTOP(JCH) = -1
          FCLOUD(JCH) = -1.d0
          !end if
          cycle channels
        end if
        
        if ( abs( RCLD(JC,J) - RCAL(JC) ) > 0.d0 )  &
             FCLOUD(JCH) = (ROBS(JC) - RCAL(JC)) /    &
             (RCLD(JC,J) - RCAL(JC))

        !**       find passage to zero if it exists and interpolate to exact pressure

        PTOP(JCH) = PLEV(J - 1) - DRAP(JCH,J - 1) /    &
             ( DRAP(JCH,J) - DRAP(JCH,J-1) ) * ( PLEV(J) - PLEV(J-1) )
        !**       find cloud radiance at zero crossing to use to get cloud fraction

        FCINT = FC(JC,J - 1) + ( FC(JC,J) - FC(JC,J - 1) ) /   &
             ( PLEV(J) - PLEV(J - 1) ) * ( PTOP(JCH) - PLEV(J - 1) )

        !**       find cloud fraction based on exact cloud top

        if ( abs(FCINT) > 0.d0 )             &
             FCLOUD(JCH) = ( RCAL(JC) - ROBS(JC) ) / FCINT

        FCLOUD(JCH) = min ( FCLOUD(JCH),  1.5d0 )
        FCLOUD(JCH) = max ( FCLOUD(JCH), -0.5d0 )

        if (FCLOUD(JCH) < 0.0d0 .or. FCLOUD(JCH) > 1.0d0 )  cycle channels
      
      end do iteration
     
    end do channels
      
  end subroutine CO2_SLICING

  subroutine SELTOP ( ETOP,VTOP,ECF,VCF,NGOOD, he,ht,cf,cfsub,ptop_mb,ps,cldflag,gncldflag )
!
!**ID SELTOP -- SELECT CLOUD TOP
!
!       AUTHOR:   L. GARAND                  July 2004
!                 A. BEAULNE (CMDA/SMC)     March 2006  (ADAPT TO 3DVAR)                 
!
!       REVISION:
!
!       OBJECT:   SELECT CLOUD TOP BY AVERAGING CO2-SLICING RESULTS
!          JUDGED CORRECT. ALL MISSING VALUES ARE -1.
!
!       ARGUMENTS:
!          INPUT:
!            -HE(NPRF)      : EQUIVALENT CLOUD TOP HEIGHTS 
!                              FROM A WINDOW CHANNEL (HPA)
!            -HT(NCO2,NPRF) : CLOUD TOPS FROM CO2-SLICING (HPA)
!            -CF(NCO2,NPRF) : EFFECTIVE CLOUD FRACTION FOR CO2-SLICING
!            -CFSUB(NPRF)   : visible ("subpixel") cloud fraction
!            -PTOP_MB(NPRF) : height (mb) from cloud_height subroutine           
!            -PS(NPRF)      : SURFACE PRESSURE IN (HPA)
!            -CLDFLAG(NPRF) : (0) CLEAR, (1) CLOUDY, (-1) UNDEFINED PROFILE
!            -NPRF          : NUMBER OF PROFILES
!
!          OUTPUT:
!            -ETOP(NPRF)    : CONSENSUS CLOUD TOP (HPA)
!            -VTOP(NPRF)    : CORRESPONDING VARIANCE ON ETOP (HPA)
!            -ECF(NPRF)     : CONSENSUS EFFECTIVE CLOUD FRACTION
!            -VCF(NPRF)     : CORRESPONDING VARIANCE ON ECF
!            -NGOOD(NPRF)   : NUMBER OF GOOD ESTIMATES
!
    implicit none
    real(8) ,intent (in) :: HE,HT(NCO2),CF(NCO2),PS,CFSUB
    integer ,intent (in) :: CLDFLAG, GNCLDFLAG
    real(8) ,intent (out):: ETOP,VTOP,ECF,VCF
    integer ,intent (out):: NGOOD
    !***********************************************************************************
    integer    :: N,JCH
    real(8)    :: PTOP_MB
    real(8)    :: H(NCO2),F(NCO2)


    ETOP = -1.d0
    VTOP = -1.d0
    ECF  = -1.d0
    VCF  = -1.d0
    NGOOD = 0

    !**     profile not assimilated if data from 2 windows channels bad
    !**     and/or if data from 2 reference co2 channels bad    
    if ( CLDFLAG == -1 ) return

    N = 0
    H(:) = 0.d0
    F(:) = 0.d0

    do JCH = 1, NCO2

      !*        CHECK FOR ZERO CLOUD FRACTION

      if ( CF(JCH) > -0.9d0 .and. CF(JCH) < 1.D-6 ) then
        N = N + 1
        H(N) = PS
        F(N) = 0.d0
      else


        !*        CONSIDER ONLY VALID VALUES OF CLOUD FRACTION ABOVE SOME THRESHOLD
        
        !         IMPORTANT LOGIC: FOR VALUES ABOVE 1.0 OF CO2-SLICING CLOUD FRACTION,
        !         SET IT TO 1.0 AND FORCE THE TOP EQUAL TO THE EFFECTIVE HEIGHT HE.
        !         CO2-SLICING NOT ALLOWED TO GIVE ESTIMATES BELOW HE, WHICH HAPPENS
        !         FOR CLOUD FRACTION CF > 1.0.

        if ( HT(JCH) > 0.0d0 ) then
          N = N + 1
          H(N) = HT(JCH)
          F(N) = min(CF(JCH), 1.0d0)
          F(N) = max(F(N), 0.d0)
          if ( CF(JCH) > 1.0d0 ) H(N) = HE
        end if
      end if

    end do


    NGOOD = N

    !*      COMPUTE MEAN AND VARIANCE

    if ( N >= 1 ) then
         
      !         ETOP = SUM(H(1:N)) / N
      !         ECF  = SUM(F(1:N)) / N

      call calcul_median_fast(N,NCO2,H,F,ETOP,ECF)
    
      VTOP = sqrt ( sum((H(1:N) - ETOP)**2) / N )
      VCF  = sqrt ( sum((F(1:N) - ECF)**2) / N )         

      if ( N == 1 ) then
        VTOP = 50.d0
        VCF  = 0.20d0
      end if
       
    else

      !*      IF NO SOLUTION FROM CO2-SLICING, AND NOT PREDETERMINED CLEAR, 
      !*      ASSUME CLOUDY WITH TOP EQUAL TO EFFECTIVE HEIGHT HE;
      !*      HOWEVER IF HE IS VERY CLOSE TO SURFACE PRESSURE PS, ASSUME CLEAR.

      ETOP = HE
      ECF  = 1.0d0
      if (CFSUB >= 0.05d0) then
        ECF = CFSUB
        ETOP = min( min(HE,PTOP_MB) , PS - 50.0d0)
      end if
      VTOP = 50.d0
      VCF  = 0.30d0
      if ( HE > (PS - 10.d0) ) ECF = 0.d0
      if ( GNCLDFLAG == 0 ) then
        ECF = 0.0d0
        ETOP = PS
      end if
    end if

    if ( ECF < 0.05d0 ) then
      ECF = 0.0d0
      ETOP = PS
    end if
  
  end subroutine SELTOP
  

  subroutine calcul_median_fast(NN,Nmax,Hin,Fin,CTP,CFR)
! 
    implicit none
    integer ,intent (in) :: NN
    integer ,intent (in) :: Nmax
    real (8) ,intent (in):: Hin(Nmax),Fin(Nmax)
    real (8) ,intent (out):: CTP,CFR
!*********************************************
    integer    :: index(NN)
    real (4) :: H(NN)
!*******
    integer :: i

    if (NN == 1) then
      CTP = Hin(NN)
      CFR = Fin(NN)
    else
      
      H(1:NN) = Hin(1:NN)
        
      call IPSORT(index,H,NN)

      if (mod(NN,2) == 0) then ! N - pair
        i = index(NN / 2)
        CTP = Hin(i)
        CFR = Fin(i)
      else                     ! N - impair
        i = index(1 + NN / 2)
        CTP = Hin(i)
        CFR = Fin(i)
      end if
    
    end if

  end subroutine calcul_median_fast

  subroutine MIN_PRES_new(MAXHEIGHT,MINP,PMIN,DT1, ps,tau,plev,cldflag,nlev,nchn,imodtop)
!
!**ID MIN_PRES -- FIND MINIMUM HEIGHT LEVEL OF SENSITIVITY
!
!       AUTHOR:   L. GARAND                   May 2004
!                 A. BEAULNE (CMDA/SMC)     March 2006  (ADAPT TO 3DVAR)                 
!
!       REVISION:
!
!       OBJECT:   FROM TOTAL TRANSMITTANCE ARRAY, FIND MINIMUM HEIGHT 
!          LEVEL OF SENSITIVITY FOR A NUMBER OF PROFILES AND CHANNELS.
!          THIS MAY BE USED TO SELECT FOR ASSIMILATION ONLY THE
!          OBSERVATIONS WITHOUT SENSITIVITY TO CLOUDS, THAT IS THE
!          RESPONSE FUNCTION SIGNIFICANT ONLY ABOVE CLOUD LEVEL.
!          THE CRITERION IS THAT dTAU/dPLEV > 0.01 FOR A 100 MB LAYER.
!
!       ARGUMENTS:
!          INPUT:
!            -PS            : SURFACE PRESSURE (HPA)
!            -TAU(NCHN,NLEV) : LAYER TO SPACE TRANSMITTANCES (0.-1.)
!            -PLEV(NLEV)          : PRESSURE LEVELS (HPA)
!            -CLDFLAG       : (0) CLEAR, (1) CLOUDY, (-1) UNDEFINED PROFILE
!            -NLEV                : NUMBER OF VERTICAL LEVELS
!            -NCHN                : NUMBER OF CHANNELS
!            -IMODTOP             : RT MODEL LEVEL NEAREST TO MODEL TOP
!
!          OUTPUT:
!            -PMIN(NCHN)     : MINIMUM HEIGHT OF SENSITIVITY (HPA)
!            -MINP(NCHN)     : VERTICAL LEVEL CORRESPONDING TO PMIN
!            -DT1(NCHN)      : VALUE OF 'DTAU/DLOGP' AT MODEL TOP
!            -MAXHEIGHT(NCHN): Height (hPa) of the maximum of the weighting function
!
    implicit none
    integer ,intent(in)   :: NCHN,NLEV,IMODTOP,CLDFLAG
    real(8), intent(in)   :: PLEV(NLEV),PS,TAU(NCHN,NLEV)
    integer, intent (out) :: MINP(NCHN)
    real(8), intent(out)  :: PMIN(NCHN), DT1(NCHN),MAXHEIGHT(NCHN)
    !*******************************************************************************

    real(8) :: MAXWF
    integer   :: J,JC,ipos(1)
    real(8)   :: WFUNC(NLEV-1),RAP(NLEV-1)

    MINP(:) = -1
    PMIN(:) = -1.d0
    DT1(:)  = -1.d0

    if ( CLDFLAG == -1 ) return

    do J = 1, NLEV - 1
      RAP(J) = log( PLEV(J + 1) / PLEV(J) )
    end do

    channels: do JC = 1, NCHN

!**       profile not assimilated if data from 2 windows channels bad
!**       and/or if data from 2 reference co2 channels bad
    
      do J = 1, NLEV
        if ( TAU(JC,J) < 0.d0) cycle channels
      end do

      MINP(JC) = NLEV
      PMIN(JC) = min(PLEV(NLEV),PS)

!*        COMPUTE ENTIRE ARRAY OF dTAU/dlog(P)
          
      do J = 1, NLEV - 1
        WFUNC(J) = (TAU(JC,J) - TAU(JC,J + 1)) / RAP(J) 
      end do
       
      DT1(JC) = WFUNC(IMODTOP)

!*        IF CHANNEL SEES THE SURFACE, DON'T RECALCULATE MINP AND PMIN

      if ( TAU(JC,NLEV) > 0.01d0 ) cycle channels

      ! Recherche du maximum
      IPOS = maxloc( WFUNC(:) )
      ! Calcul de la valeur du maximum
      MAXWF = WFUNC(IPOS(1))
      ! maximum entre les 2 niveaux puisque WF calculee pour une couche finie ( discutable ?)
      MAXHEIGHT(JC)= 0.5d0 * ( PLEV(IPOS(1)) +  PLEV(IPOS(1) + 1)  )

      !*        IF CHANNEL DOESN'T SEE THE SURFACE, SEE WHERE dTAU/dlog(PLEV) BECOMES IMPORTANT
      !*        FOR RECOMPUTATION OF MINP AND PMIN.

      do J = NLEV - 1, IPOS(1), -1
        if ( ( WFUNC(J)/ MAXWF ) > 0.01d0) then
          MINP(JC) = J + 1
          PMIN(JC) = min(PLEV(J + 1),PS)
          exit
        end if
      end do
     
    end do channels

  end subroutine MIN_PRES_NEW

  subroutine CLOUD_HEIGHT (PTOP,NTOP, &
       btobs,cldflag,tt,gz,ps,plev,nlev, &
       nchn,ichref,lev_start,iopt)
!
!**ID CLOUD_HEIGHT -- CLOUD TOP HEIGHT COMPUTATION
!
!       SCIENCE:  L. GARAND
!       AUTHOR:   A. BEAULNE (CMDA/SMC)   August 2004
!                 A. BEAULNE (CMDA/SMC) February 2006  (ADAPT TO 3DVAR)                 
!
!       REVISION:
!
!       OBJECT:   COMPUTATION OF CLOUD TOP HEIGHT (ABOVE THE GROUND)
!          BASED ON MATCHING OBSERVED BRIGHTNESS TEMPERATURE AT A 
!          REFERENCE SURFACE CHANNEL WITH BACKGROUND TEMPERATURE PROFILE.
!          TO USE WITH ONE REFERENCE CHANNEL. USED HERE ON MODEL LEVELS.
!
!       ARGUMENTS:
!          INPUT:
!            -BTOBS(NCHN) : OBSERVED BRIGHTNESS TEMPERATURE (DEG K)
!            -CLDFLAG    : CLEAR(0), CLOUDY(1), UNDEFINED(-1) PROFILES
!            -TT(NLEV)    : TEMPERATURE PROFILES (DEG K)
!            -GZ(NLEV)    : HEIGHT PROFILES ABOVE GROUND (M)
!            -PS(NPRF)         : SURFACE PRESSURE (HPA)
!            -PLEV(NLEV)  : PRESSURE LEVELS (HPA)
!            -NLEV             : NUMBER OF VERTICAL LEVELS
!            -NCHN             : NUMBER OF CHANNELS
!            -ICHREF     : CHOSEN REFERENCE SURFACE CHANNEL
!            -IOPT             : LEVELS USING PLEV (1) OR GZ (2)
!
!
!          INPUT/OUTPUT:
!            -LEV_START : LEVEL TO START ITERATION (IDEALLY TROPOPAUSE)
!
!          OUTPUT:
!            -PTOP    : CHOSEN EQUIVALENT CLOUD TOPS 
!                             (IN HPA|M WITH IOPT = 1|2) 
!            -NTOP    : NUMBER OF POSSIBLE PTOP SOLUTIONS
!
!
    implicit none
    integer ,intent (in) :: NCHN,NLEV,IOPT,ICHREF,CLDFLAG
    real(8) ,intent (in) :: BTOBS(NCHN),TT(NLEV),GZ(NLEV),PS,PLEV(NLEV)
    integer ,intent (inout) :: LEV_START
    real(8) ,intent (out) :: PTOP
    integer ,intent (out) :: NTOP
    !**********************************************************************************************

    integer     :: JN 
    integer     :: ITOP
    integer     :: NHT
    real(8)     :: HT(NLEV)
 
    if ( IOPT == 1 ) then
     
      PTOP = PS
      NTOP = 1      

      if ( CLDFLAG == -1 ) return
      
      call GET_TOP ( HT,NHT, btobs(ichref),tt,plev,nlev,lev_start,iopt ) 

      ITOP = 1
      if ( NHT >= 2 ) ITOP = 2
      PTOP = min ( HT(ITOP), PS )
      NTOP = NHT

    else if ( IOPT == 2 ) then
      
      PTOP = 0.d0
      NTOP = 1      

      if ( CLDFLAG == -1 ) return

      call GET_TOP ( HT,NHT, btobs(ichref),tt,gz,nlev,lev_start,iopt )

      ITOP = 1
      if ( NHT >= 2 ) ITOP = 2
      PTOP = max ( HT(ITOP), 0.d0 )
      NTOP = NHT
       
    end if

  end subroutine CLOUD_HEIGHT
 
  subroutine GARAND1998NADON (CLDFLAG, btobs,tg,tt,gz,nlev, &
       nchn,ptop_eq,ntop_eq,ichref)
!
!**ID GARAND1998NADON -- DETERMINE IF PROFILES ARE CLEAR OR CLOUDY
!
!       SCIENCE:  L. GARAND AND S. NADON
!       AUTHOR:   A. BEAULNE (CMDA/SMC)      June 2004
!                 A. BEAULNE (CMDA/SMC)     March 2006  (ADAPT TO 3DVAR)                 
!
!       REVISION:
!
!       OBJECT:   DETERMINE IF THE PROFILES ARE CLEAR OR CLOUDY BASED ON
!          THE ALGORITHM OF GARAND & NADON 98 J.CLIM V11 PP.1976-1996
!          WITH CHANNEL IREF
!
!       ARGUMENTS:
!          INPUT:
!            -BTOBS(NCHN) : OBSERVED BRIGHTNESS TEMPERATURES (DEG K)
!            -TG         : GUESS SKIN TEMPERATURES (DEG K)
!            -TT(NLEV)    : GUESS TEMPERATURE PROFILES (DEG K)
!            -GZ(NLEV)    : GUESS HEIGHT PROFILE ABOVE GROUND (M)
!            -NLEV             : NUMBER OF VERTICAL LEVELS
!            -NCHN             : NUMBER OF CHANNELS
!            -PTOP_EQ    : CHOSEN EQUIVALENT CLOUD TOPS (M)
!            -NTOP_EQ    : NUMBER OF POSSIBLE PTOP_EQ SOLUTIONS
!            -ICHREF     : CHOSEN REFERENCE SURFACE CHANNEL
!
!          INPUT/OUTPUT:
!            -CLDFLAG(NPRF)  : CLEAR(0), CLOUDY(1), UNDEFINED(-1) PROFILES
!
    implicit none
    integer ,intent (in) :: NLEV,NCHN
    real(8) ,intent (in) :: BTOBS(NCHN),TG,GZ(NLEV),TT(NLEV),PTOP_EQ
    integer ,intent (in) :: NTOP_EQ,ICHREF
    integer ,intent (inout) :: CLDFLAG
!*********************************************************************************************
    integer    :: NINV
    real(8)    :: LEV(2)

      
    LEV(1) = 222.d0
    LEV(2) = 428.d0


    if ( CLDFLAG == -1 ) return

    if ( BTOBS(ICHREF) >= TG - 3.d0 .and. BTOBS(ICHREF) <= TG + 3.d0 ) then
      CLDFLAG = 0
      return
    end if

    if ( BTOBS(ICHREF) >= TG - 4.d0 .and. BTOBS(ICHREF) <= TG - 3.d0 ) then
      if ( PTOP_EQ > 1100.d0 ) then
        CLDFLAG = 1
        return
      else
        CLDFLAG = 0
        return
      end if
    end if
    
    if ( PTOP_EQ > 728.d0 ) then
      CLDFLAG = 1
      return
    end if

    if ( TG - BTOBS(ICHREF) > 8.d0 ) then 
      if ( NTOP_EQ >= 3 ) then
        if ( PTOP_EQ > 73.d0 ) then
          CLDFLAG=1
          return
        else
          CLDFLAG=0
          return
        end if
      else
        call MONOTONIC_INVERSION (NINV, tg,tt,gz,nlev,lev(1))
        if ( NINV == 1 ) then
          if ( PTOP_EQ > 222.d0 ) then
            CLDFLAG = 1
            return
          else
            CLDFLAG = 0 
            return
          end if
        else
          CLDFLAG = 0
          return
        end if
      end if
    end if
    
    if ( TG - BTOBS(ICHREF) > 5.d0 ) then
      if ( NTOP_EQ >= 3 ) then
        if ( PTOP_EQ > 222.d0 ) then
          CLDFLAG = 1
          return
        else
          CLDFLAG = 0
          return
        end if
      else
        call MONOTONIC_INVERSION (NINV, tg,tt,gz,nlev,lev(2))
        if ( NINV == 1) then
          if( PTOP_EQ > 428.d0 ) then
            CLDFLAG = 1
            return
          else
            CLDFLAG = 0
            return
          end if
        else
          CLDFLAG = 0
        end if
      end if
    else
      CLDFLAG = 0
    end if
    
  end subroutine GARAND1998NADON

  subroutine MONOTONIC_INVERSION (NINVR, ptg,ptt,pgz,npr,lvl)

!***********************************************************************
!
!**ID MONOTONIC_INVERSION -- DETECT TEMPERATURE INVERSION
!
!       SCIENCE:  L. GARAND AND S. NADON
!       AUTHOR:   A. BEAULNE (CMDA/SMC)      June 2004
!                 A. BEAULNE (CMDA/SMC)     March 2006  (ADAPT TO 3DVAR)                 
!
!       REVISION:
!
!       OBJECT:   DETERMINE IF THERE IS A PRESENCE (NINVR=1) OR NOT (NINVR=0)
!           OF A TEMPERATURE INVERSION GOING FROM THE SURFACE UP TO THE
!           HEIGHT LVL
!
!       ARGUMENTS:
!          INPUT:
!            -PTG       : SKIN TEMPERATURE (DEG K)
!            -PTT(NPR) : TEMPERATURE PROFILE (DEG K)
!            -PGZ(NPR) : HEIGHT PROFILE ABOVE GROUND (M)
!            -NPR     : NUMBER OF VERTICAL LVLELS
!            -LVL      : HEIGHT TO SEARCH FOR TEMPERATURE INVERSION (M)
!
!          OUTPUT:
!            -NINVR     : PRESENCE (1) OR NOT (0) OF A TEMPERATURE INVERSION
!                         FROM THE SURFACE TO HEIGHT LVL
!
!
!***********************************************************************

    implicit none
    integer ,intent (in) :: npr
    real(8),intent (in)  :: PTT(NPR),PGZ(NPR),PTG,LVL
    integer ,intent (out):: ninvr
!**************************************************
    integer   :: NL

    NINVR = 0
    if ( PTG - PTT(NPR) < 0.d0 ) then
      NINVR = 1
      do NL = NPR - 1, 1, -1
        if ( PGZ(NL) > LVL ) exit
        if ( PTT(NL+1) - PTT(NL) > 0.d0 ) then
          NINVR = 0
          exit
        end if
      end do
    end if

  end subroutine MONOTONIC_INVERSION

  subroutine ESTIM_TS(TS, tg,emi,rcal,radobs,sfctau,cldflag, &
       ichref,nchnkept,myCoefs)

!
!**ID ESTIM_TS -- GET AN ESTIMATED SKIN TEMPERATURE
!
!       AUTHOR:   L. GARAND                   May 2004
!                 A. BEAULNE (CMDA/SMC)     March 2006  (ADAPT TO 3DVAR)                 
!
!       REVISION:
!
!       OBJECT:   GET AN ESTIMATED SKIN TEMPERATURE BY INVERSION OF
!          RADIATIVE TRANSFER EQUATION ASSUMING GUESS T AND Q PROFILES
!          ARE PERFECT. DESIGNED FOR A SINGLE CHANNEL ICHREF AND NPRF
!          PROFILES. ASSUMES A REAL TG (GUESS) OVER OCEANS AND A TG 
!          WITH HYPOTHESIS OF UNITY EMISSIVITY OVER LAND.
!      
!          USES:  RCAL = B(TG)*EMI*SFCTAU + ATMOS_PART
!             TS = B(TS)*EMI*SFCTAU + ATMOS_PART
!          SOLVES FOR TS
!
!       ARGUMENTS:
!          INPUT:
!            -TG          : GUESS SKIN TEMPERATURE (DEG K)
!            -EMI(NCHNKEPT)    : SURFACE EMISSIVITIES FROM WINDOW CHANNEL (0.-1.)
!            -RCAL(NCHNKEPT)   : COMPUTED CLEAR RADIANCES (MW/M2/SR/CM-1)
!            -RADOBS(NCHNKEPT) : OBSERVED RADIANCES (")
!            -SFCTAU(NCHNKEPT) : SURFACE TO SPACE TRANSMITTANCES (0.-1.)
!            -CLDFLAG     : CLEAR(0), CLOUDY(1) OR UNDEFINED(-1) PROFILES
!            -ICHREF      : REFERENCE SURFACE CHANNEL (SUBSET VALUES)
!            -NCHNKEPT          : NUMBER OF CHANNELS KEPT IN CMA
!
!          OUTPUT:
!            -TS          : RETRIEVED SKIN TEMPERATURE (-1. FOR MISSING)
!
    implicit none
    integer ,intent(in) :: NCHNKEPT
    integer ,intent(in) :: ICHREF,CLDFLAG
    real(8) ,intent(in) :: TG,EMI(NCHNKEPT),RCAL(NCHNKEPT),RADOBS(NCHNKEPT)
    real(8) ,intent(in) :: SFCTAU(NCHNKEPT)
    real(8) ,intent(out):: TS
    type( rttov_coefs ) ,intent(in) :: myCoefs
!************************************************************************************
    real(8)    :: RTG,RADTG
    real(8)    :: RADTS,tstore,t_effective
  
    TS = -1.d0

    if ( CLDFLAG /= 0 ) return
    if ( ichref == -1 ) return


!*    transform guess skin temperature to plank radiances 

    t_effective =  myCoefs % coef % ff_bco(ichref) + myCoefs % coef % ff_bcs(ichref) * TG

    RADTG =  myCoefs % coef % planck1(ichref) / &
         ( exp( myCoefs % coef % planck2(ichref) / t_effective ) - 1.0d0 )


    !*   compute TOA planck radiances due to guess skin planck radiances

    RTG = RADTG * EMI(ICHREF) * SFCTAU(ICHREF)


    !*   compute true skin planck radiances due to TOA true planck radiances
    
    RADTS = ( RADOBS(ICHREF) + RTG - RCAL(ICHREF) ) / &
         ( EMI(ICHREF) * SFCTAU(ICHREF) )

    if (RADTS <= 0.d0) then
      Write(*,'(A25,1x,8e14.6)') "Warning, negative radts", RADOBS(ICHREF), RTG, RCAL(ICHREF), EMI(ICHREF), SFCTAU(ICHREF), &
           ( RADOBS(ICHREF) + RTG - RCAL(ICHREF) ), ( EMI(ICHREF) * SFCTAU(ICHREF) ), RADTS
      Write(*,*) "Skipping tskin retrieval."
      return
    end if

    
    !*   transform true skin planck radiances to true skin temperatures

    tstore = myCoefs % coef % planck2(ichref) / log( 1.0d0 + myCoefs % coef % planck1(ichref) / RADTS )

    TS = ( tstore - myCoefs % coef % ff_bco(ichref) ) / myCoefs % coef % ff_bcs(ichref)
    

  end subroutine ESTIM_TS


  subroutine CLOUD_TOP ( PTOP_BT,PTOP_RD,NTOP_BT,NTOP_RD,  &
       btobs,tt,gz,rcal,ps,robs,rcld,plev,nlev,nchn, &
       cldflag,rejflag,lev_start,iopt,ihgt,ichref,nch,ilist)
!
!**ID CLOUD_TOP -- CLOUD TOP HEIGHT COMPUTATION
!
!       AUTHOR:   L. GARAND             August 2004
!                 A. BEAULNE (CMDA/SMC)  March 2006  (ADAPT TO 3DVAR)      
!                
!       REVISION:  001 S. Heilliette: removal of hard-coded rttov level
!
!       OBJECT:   COMPUTATION OF CLOUD TOP HEIGHT (ABOVE THE GROUND)
!          BASED ON MATCHING OBSERVED BRIGHTNESS TEMPERATURE WITH 
!          BACKGROUND TEMPERATURE PROFILES AND/OR COMPUTED OBSERVED
!          RADIANCES WITH BACKGROUND RADIANCE PROFILES.
!          TO USE WITH MORE THAN ONE CHANNEL. USED HERE ON RTTOV LEVELS.
!
!       ARGUMENTS:
!          INPUT:
!            -BTOBS(NCHN)     : OBSERVED BRIGHTNESS TEMPERAUTRES (DEG K)
!            -TT(NLEV)        : TEMPERATURE PROFILES (DEG K)
!            -GZ(NLEV)        : HEIGHT PROFILES ABOVE GROUND (M)
!            -RCAL(NCHN)      : COMPUTED CLEAR RADIANCES (MW/M2/SR/CM-1)
!            -PS            : SURFACE PRESSURE (HPA)
!            -ROBS(NCHN)      : COMPUTED OBSERVED RADIANCES (MW/M2/SR/CM-1)
!            -RCLD(NCHN,NLEV) : COMPUTED CLOUD RADIANCES FROM EACH LEVEL (")
!            -PLEV(NLEV)           : PRESSURE LEVELS (HPA)
!            -NLEV                 : NUMBER OF VERTICAL LEVELS
!            -NCHN                 : NUMBER OF CHANNELS
!            -CLDFLAG        : CLEAR(0), CLOUDY(1), UNDEFINED(-1) PROFILES
!            -REJFLAG(NCHN,0:BITFLAG) : FLAGS FOR REJECTED OBSERVATIONS
!            -IOPT                 : LEVELS USING PLEV (1) OR GZ (2)
!            -IHGT                 : GET *_BT* ONLY (0), *_RD* ONLY (1), BOTH (2)
!            -ICHREF         : REFERENCE SURFACE CHANNEL (SUBSET VALUE)
!            -NCH                  : NUMBER OF CHANNELS WE WANT OUTPUTS
!            -ILIST(NCH )          : LIST OF THE CHANNEL NUMBERS (SUBSET VALUES) 
!
!          INPUT/OUTPUT:
!            -LEV_START      : LEVEL TO START ITERATION (IDEALLY TROPOPAUSE)
!
!          OUTPUT:
!            -PTOP_BT(NCHN)  : CHOSEN EQUIVALENT CLOUD TOPS BASED ON 
!                                    BRIGHTNESS TEMPERATURES (IN HPA|M WITH IOPT = 1|2)
!            -PTOP_RD(NCHN)  : CHOSEN EQUIVALENT CLOUD TOPS BASED ON 
!                                    RADIANCES (IN HPA|M WITH IOPT = 1|2)
!            -NTOP_BT       : NUMBER OF POSSIBLE PTOP_BT SOLUTIONS
!            -NTOP_RD       : NUMBER OF POSSIBLE PTOP_RD SOLUTIONS
!
    implicit none
    integer, intent (in) :: NCHN,NCH,NLEV,IOPT,IHGT
    real(8), intent (in) :: BTOBS(NCHN),RCLD(NCHN,NLEV)
    real(8), intent (in) :: ROBS(NCHN),RCAL(NCHN)
    real(8), intent (in) :: TT(NLEV),GZ(NLEV),PLEV(NLEV),PS
    integer, intent (in) :: REJFLAG(NCHN,0:BITFLAG),ILIST(NCH),CLDFLAG,ICHREF
    integer, intent (inout) :: LEV_START
    real(8), intent (out) ::  PTOP_BT(NCHN),PTOP_RD(NCHN)
    integer, intent (out) ::  NTOP_BT(NCHN),NTOP_RD(NCHN)
    !******************************************************************
    integer      ::  JCH,JC,ITOP,NHT,i10,i
    integer      ::  SUMREJ
    real(8)      ::  HT(NLEV)

      
    i10=1
    do I=2,NLEV
      if (plev(i - 1) <= 100.d0 .and. plev(i) > 100.d0) then
        I10 = I
        exit
      end if
    end do

    PTOP_BT(:) = -10.d0
    PTOP_RD(:) = -10.d0

    NTOP_BT(:) = 0.d0
    NTOP_RD(:) = 0.d0

    !**     profile not assimilated if data from 2 windows channels bad

    if ( CLDFLAG == -1 ) return

    !**     predetermined clear

    SUMREJ = sum( REJFLAG(ICHREF,:) )

    if ( SUMREJ == 0 ) then
      
      if ( IOPT == 1 ) then
        PTOP_BT(:) = min ( PLEV(NLEV), PS )
        PTOP_RD(:) = min ( PLEV(NLEV), PS )
      else if ( IOPT == 2 ) then
        PTOP_BT(:) = 0.d0
        PTOP_RD(:) = 0.d0
      end if
     
      NTOP_BT(:) = 1
      NTOP_RD(:) = 1
     
      LEV_START = max ( LEV_START , i10 )
    
      return

    end if


    channels: do JCH = 1, NCH
       
      JC = ILIST(JCH)
       
      !**       gross check failure

      if ( REJFLAG(JC,9) == 1 ) cycle channels

      !**       no clouds if observed radiance warmer than clear estimate

      if ( ROBS(JC) > RCAL(JC) ) then

        if ( IOPT == 1 ) then
          PTOP_BT(JC) = min ( PLEV(NLEV), PS )
          PTOP_RD(JC) = min ( PLEV(NLEV), PS )
        else if ( IOPT == 2 ) then
          PTOP_BT(JC) = 0.d0
          PTOP_RD(JC) = 0.d0
        end if
      
        NTOP_BT(JC) = 1
        NTOP_RD(JC) = 1

        cycle channels
           
      end if

      !**       cloudy
      
      if ( REJFLAG(JC,11) == 1 .and. REJFLAG(JC,23) == 1 ) then
        
        if ( IOPT == 1 ) then
          
          if ( IHGT == 0 .or. IHGT == 2 ) then
            call GET_TOP ( HT,NHT, btobs(jc),tt,plev,nlev,lev_start,iopt) 
            ITOP = 1
            if ( NHT >= 2 ) ITOP = 2
            PTOP_BT(JC) = min ( HT(ITOP), PS )
            NTOP_BT(JC) = NHT
          end if
          
          if ( IHGT == 1 .or. IHGT == 2 ) then
            call GET_TOP ( HT,NHT, robs(jc),rcld(jc,:),plev,nlev,lev_start,iopt)
            ITOP = 1
            if ( NHT >= 2 ) ITOP = 2
            PTOP_RD(JC) = min ( HT(ITOP), PS )
            NTOP_RD(JC) = NHT
          end if
          
        else if ( IOPT == 2 ) then 
          
          if ( IHGT == 0 .or. IHGT == 2 ) then
            call GET_TOP ( HT,NHT, btobs(jc),tt,gz,nlev,lev_start,iopt) 
            ITOP = 1
            if ( NHT >= 2 ) ITOP = 2
            PTOP_BT(JC) = max ( HT(ITOP), 0.d0 )
            NTOP_BT(JC) = NHT
          end if
          
          if ( IHGT == 1 .or. IHGT == 2 ) then
            call GET_TOP ( HT,NHT, robs(jc),rcld(jc,:),gz,nlev,lev_start,iopt)
            ITOP = 1
            if ( NHT >= 2 ) ITOP = 2
            PTOP_RD(JC) = max ( HT(ITOP), 0.d0 )
            NTOP_RD(JC) = NHT
          end if
          
        end if
        
      end if
      
    end do channels

  end subroutine CLOUD_TOP

  subroutine CLOUD_TOP_AVHRR ( PTOP_BT,PTOP_RD,NTOP_BT,NTOP_RD,  &
       btobs,tt,gz,rcal,ps,robs,rcld,plev,nlev,nchn, &
       cldflag,lev_start,iopt,ihgt,nch,ilist)
!
!**ID CLOUD_TOP -- CLOUD TOP HEIGHT COMPUTATION
!
!       AUTHOR:   L. GARAND             August 2004
!                 A. BEAULNE (CMDA/SMC)  March 2006  (ADAPT TO 3DVAR)                 
!
!       REVISION: 001 S. Heilliette:
!                     -to remove hard-coded rttov pressure level numbers
!
!       OBJECT:   COMPUTATION OF CLOUD TOP HEIGHT (ABOVE THE GROUND)
!          BASED ON MATCHING OBSERVED BRIGHTNESS TEMPERATURE WITH 
!          BACKGROUND TEMPERATURE PROFILES AND/OR COMPUTED OBSERVED
!          RADIANCES WITH BACKGROUND RADIANCE PROFILES.
!          TO USE WITH MORE THAN ONE CHANNEL. USED HERE ON RTTOV LEVELS.
!
!       ARGUMENTS:
!          INPUT:
!            -BTOBS(NCHN)     : OBSERVED BRIGHTNESS TEMPERAUTRES (DEG K)
!            -TT(NLEV)        : TEMPERATURE PROFILES (DEG K)
!            -GZ(NLEV)        : HEIGHT PROFILES ABOVE GROUND (M)
!            -RCAL(NCHN)      : COMPUTED CLEAR RADIANCES (MW/M2/SR/CM-1)
!            -PS             : SURFACE PRESSURE (HPA)
!            -ROBS(NCHN)      : COMPUTED OBSERVED RADIANCES (MW/M2/SR/CM-1)
!            -RCLD(NCHN,NLEV) : COMPUTED CLOUD RADIANCES FROM EACH LEVEL (")
!            -PLEV(NLEV)           : PRESSURE LEVELS (HPA)
!            -NLEV                 : NUMBER OF VERTICAL LEVELS
!            -NCHN                 : NUMBER OF CHANNELS
!            -CLDFLAG        : CLEAR(0), CLOUDY(1), UNDEFINED(-1) PROFILES
!            -IOPT                 : LEVELS USING PLEV (1) OR GZ (2)
!            -IHGT                 : GET *_BT* ONLY (0), *_RD* ONLY (1), BOTH (2)
!            -NCH                  : NUMBER OF CHANNELS WE WANT OUTPUTS
!            -ILIST(NCH)           : LIST OF THE CHANNEL NUMBERS (SUBSET VALUES) 
!
!          INPUT/OUTPUT:
!            -LEV_START      : LEVEL TO START ITERATION (IDEALLY TROPOPAUSE)
!
!          OUTPUT:
!            -PTOP_BT(NCHN)  : CHOSEN EQUIVALENT CLOUD TOPS BASED ON 
!                                    BRIGHTNESS TEMPERATURES (IN HPA|M WITH IOPT = 1|2)
!            -PTOP_RD(NCHN)  : CHOSEN EQUIVALENT CLOUD TOPS BASED ON 
!                                    RADIANCES (IN HPA|M WITH IOPT = 1|2)
!            -NTOP_BT(NCHN)  : NUMBER OF POSSIBLE PTOP_BT SOLUTIONS
!            -NTOP_RD(NCHN)  : NUMBER OF POSSIBLE PTOP_RD SOLUTIONS
!
    implicit none
    integer ,intent(in) :: NCH,IOPT,IHGT,NLEV,NCHN
    integer ,intent(in) :: ILIST(NCH),CLDFLAG
    real(8) ,intent(in) ::  PLEV(NLEV),PS
    real(8) ,intent(in) ::  ROBS(NCHN),RCAL(NCHN)
    real(8) ,intent(in) ::  BTOBS(NCHN),RCLD(NCHN,NLEV)
    real(8) ,intent(in) ::  TT(NLEV),GZ(NLEV)
    integer ,intent(inout) :: LEV_START
    real(8) ,intent(out) ::  PTOP_BT(NCHN),PTOP_RD(NCHN)
    integer ,intent(out) ::  NTOP_BT(NCHN),NTOP_RD(NCHN)
    !*********************************************************************
    integer      ::  JCH,JC,ITOP,NHT,i10,i
    real(8)      ::  HT(NLEV)
     
    i10 = 1
    do I=2,NLEV
      if (plev(i - 1) <= 100.d0 .and. plev(i) > 100.d0) then
        I10 = I
        exit
      end if
    end do
    
    PTOP_BT(:) = -10.d0
    PTOP_RD(:) = -10.d0

    NTOP_BT(:) = 0.d0
    NTOP_RD(:) = 0.d0


    !**     profile not assimilated if data from 2 windows channels bad

    if ( CLDFLAG == -1 ) return

    !**     predetermined clear

        
    if ( CLDFLAG == 0 ) then
        
      if ( IOPT == 1 ) then
        PTOP_BT(:) = min ( PLEV(NLEV), PS )
        PTOP_RD(:) = min ( PLEV(NLEV), PS )
      else if ( IOPT == 2 ) then
        PTOP_BT(:) = 0.d0
        PTOP_RD(:) = 0.d0
      end if
      
      NTOP_BT(:) = 1
      NTOP_RD(:) = 1
      
      LEV_START = max ( LEV_START , i10 )
       
      return

    end if

    channels: do JCH = 1, NCH

      JC = ILIST(JCH)

      !**       gross check failure
      
      if ( BTOBS(JC) < 150.d0 .or. BTOBS(JC) > 350.d0) cycle channels

      !**       no clouds if observed radiance warmer than clear estimate
      
      if ( ROBS(JC) > RCAL(JC) ) then
          
        if ( IOPT == 1 ) then
          PTOP_BT(JC) = min ( PLEV(NLEV), PS )
          PTOP_RD(JC) = min ( PLEV(NLEV), PS )
        else if ( IOPT == 2 ) then
          PTOP_BT(JC) = 0.d0
          PTOP_RD(JC) = 0.d0
        end if
        
        NTOP_BT(JC) = 1
        NTOP_RD(JC) = 1
        
        cycle channels
        
      end if

!**       cloudy

      if ( CLDFLAG == 1 ) then
        
        if ( IOPT == 1 ) then

          if ( IHGT == 0 .or. IHGT == 2 ) then
            call GET_TOP ( HT,NHT, btobs(jc),tt,plev,nlev,lev_start,iopt) 
            ITOP = 1
            if ( NHT >= 2 ) ITOP = 2
            PTOP_BT(JC) = min ( HT(ITOP), PS )
            NTOP_BT(JC) = NHT
          end if
              
          if ( IHGT == 1 .or. IHGT == 2 ) then
            call GET_TOP ( HT,NHT, robs(jc),rcld(jc,:),plev,nlev,lev_start,iopt)
            ITOP = 1
            if ( NHT >= 2 ) ITOP = 2
            PTOP_RD(JC) = min ( HT(ITOP), PS )
            NTOP_RD(JC) = NHT
          end if
          
        else if ( IOPT == 2 ) then 
          
          if ( IHGT == 0 .or. IHGT == 2 ) then
            call GET_TOP ( HT,NHT, btobs(jc),tt,gz,nlev,lev_start,iopt) 
            ITOP = 1
            if ( NHT >= 2 ) ITOP = 2
            PTOP_BT(JC) = max ( HT(ITOP), 0.d0 )
            NTOP_BT(JC) = NHT
          end if
          
          if ( IHGT == 1 .or. IHGT == 2 ) then
            call GET_TOP ( HT,NHT, robs(jc),rcld(jc,:),gz,nlev,lev_start,iopt)
            ITOP = 1
            if ( NHT >= 2 ) ITOP = 2
            PTOP_RD(JC) = max ( HT(ITOP), 0.d0 )
            NTOP_RD(JC) = NHT
          end if
              
        end if
           
      end if

    end do channels

  end subroutine CLOUD_TOP_AVHRR

  subroutine GET_TOP (HT,NHT, bt,tt,pp,nlev,lev_start,iopt)

!**************************************************************************
!
!**ID GET_TOP -- CLOUD TOP HEIGHT COMPUTATION
!
!       AUTHOR:   L. GARAND                       2004
!                 A. BEAULNE (CMDA/SMC)  February 2006  (ADAPT TO 3DVAR)                 
!
!       REVISION:
!
!       OBJECT:   COMPUTATION OF CLOUD TOP HEIGHT AND NUMBER OF POSSIBLE HEIGHTS
!
!       ARGUMENTS:
!          INPUT:
!            -BT        : OBSERVED BRIGHTNESS TEMPERATURES (DEG K)
!                          OR COMPUTED OBSERVED RADIANCES (MW/M2/SR/CM-1)
!            -TT(NLEV)  : TEMPERATURE PROFILE (DEG K)
!                          OR COMPUTED CLOUD RADIANCE FROM EACH LEVEL TO TOP (")
!            -PP(NLEV)  : PRESSURE (HPA) OR HEIGHTS (M) PROFILE (IOPT=1 OR 2)
!            -NLEV      : NUMBER OF VERTICAL LEVELS
!            -IOPT      : HEIGHT UNITS IN HPA (1) OR IN METERS (2)
!
!          INPUT/OUTPUT:
!            -LEV_START : LEVEL TO START ITERATION (IDEALLY TROPOPAUSE)
!                         (IF <= 0, SEARCH & START AT COLDEST LEVEL)
!
!          OUTPUT:
!            -HT(NLEV)  : CLOUD TOP HEIGHT IN HPA OR METERS (IOPT = 1 OR 2)
!            -NHT       : NUMBER OF POSSIBLE CLOUD HEIGHT SOLUTIONS
!
!***************************************************************************

    implicit none
    integer ,intent   (in) :: NLEV,IOPT
    integer ,intent (inout) :: lev_start
    real(8)  ,intent (in) :: BT, TT(NLEV), PP(NLEV)
    real(8)  ,intent (out) :: HT(NLEV)
    integer ,intent (out) :: NHT
!*******************************************************
     
    integer   :: I, IM, i10
    real(8)   :: P(NLEV)
    real(8)   :: DT, A, AA, B
    
    HT(:) = -1.
    
    if (IOPT == 1) P(:) = log(PP(:))


    IM = LEV_START


    if ( LEV_START <= 0 ) then

      !*      SEARCH INDEX IM WHERE TT IS MINIMUM

      call FMIN (IM, tt, nlev)
  
      i10 = -1
      do I=2,NLEV
        if (pp(i-1) <= 100.d0 .and. pp(i) > 100.d0) then
          I10 = I
          exit
        end if
      end do
       
      LEV_START = IM
     
      if ( IM == NLEV ) then
        LEV_START = max(LEV_START,i10)
        NHT = 1
        HT(1) = PP(NLEV)
        return
      end if
       
    end if


    NHT = 0        
    
    do I = IM, NLEV - 1
      DT = TT(I + 1) - TT(I) + 1.D-12
      if ( BT > TT(I) .and. BT <= TT(I + 1) ) then
        
        NHT = NHT + 1
        
        if (IOPT == 1) then
          A = P(I) + (P(I + 1) - P(I)) / DT * ( BT - TT(I))
          HT(NHT) = exp(A)
        end if

        if (IOPT == 2) then
          B  = PP(I) + (PP(I+1) - PP(I)) / DT * (BT - TT(I))
          HT(NHT) = B
        end if
         
      else if ( BT >= TT(I+1) .and. BT < TT(I) ) then
      
        NHT = NHT + 1
      
        if (IOPT == 1) then
          A  = P(I + 1)- (P(I + 1)-P(I)) / DT * (TT(I + 1) - BT)
          HT(NHT) = exp(A)
        end if

        if(IOPT == 2) then
          B = PP(I + 1)- (PP(I + 1) - PP(I)) / DT * (TT(I + 1) - BT)
          HT(NHT) = B
        end if
       
      end if
    end do
    
    
    if ( NHT == 0 .and. BT < TT(IM) )  then
      NHT  = 1
      HT(1) = PP(IM)
    else if ( NHT == 0 .and. BT > TT(NLEV) )  then
      NHT   = 1
      HT(1) = PP(NLEV)
    end if
    
  end subroutine GET_TOP

  subroutine FMIN ( IMIN, f, ndim )

!***********************************************************************
!
!**ID FMIN -- SEARCH MINIMUM VALUE OF VECTOR
!
!       AUTHOR:   L. GARAND
!                 A. BEAULNE (CMDA/SMC)  February 2006  (ADAPT TO 3DVAR)                 
!
!       REVISION:
!
!       OBJECT:   SEARCH THE POSITION IN VECTOR F WHERE VALUE IF MINIMUM
!
!       ARGUMENTS:
!          INPUT:
!            -F     : 1D VECTOR
!            -NDIM  : VECTOR DIMENSION
!
!          OUTPUT:
!            -IM    : INDEX OF VECTOR F WHERE VALUE IS MINIMUM
!
!
!***********************************************************************
    implicit none
    integer,intent (in) :: ndim
    real(8),intent (in) :: f(ndim)
    integer,intent (out):: imin
    !*********************************
    integer  :: I
    real(8)  :: X


    IMIN = 1
    X = F(1)
    
    do I=2,NDIM
      if (F(I) < X) then
        X = F(I)
        IMIN = I
      end if
    end do
    

  end subroutine FMIN

  subroutine get_avhrr_emiss(iasi_surfem1,freqiasi,nchaniasi,avhrr_surfem1)
    ! choisi l'emissivite d'un canal IASI proche pour AVHRR
    ! a raffiner pour prendre en  compte la largeur  des canaux AVHRR ??
    implicit none
    integer ,intent(in) :: nchaniasi
    real (8) ,intent (in) :: iasi_surfem1 ( nchaniasi )
    real (8) ,intent (in) :: freqiasi( nchaniasi )
    real (8) ,intent (out):: avhrr_surfem1( NIR )
    !****************************
    real (8),parameter :: freqavhrr(NIR)= (/0.2687000000D+04 , 0.9272000000D+03 , 0.8377000000D+03/)
    integer,save :: indxavhrr(NIR)
    integer :: i,pos(1)
    !*************************************************************8
    do I=1,NIR
      pos = minloc ( abs (freqiasi(:) - freqavhrr(I)) )
      indxavhrr(i) = pos(1)
    end do

    do I=1,NIR
      avhrr_surfem1(i) = iasi_surfem1(indxavhrr(i))
    end do
  
  end subroutine get_avhrr_emiss
 
  subroutine tovs_rttov_AVHRR_for_IASI (iptobs,surfem1_avhrr,idiasi)

!
!**s/r tovs_rttov_AVHRR_for_IASI  - Computation of forward radiance with rttov_direct
!                   (for AVHRR)
!
!
!author        : S. Heilliette
!
!revision 001  : s. heilliette october 2010
!                  - adaptation to rttov 10.0
!revision 002  : s. heilliette may 2017
!                  - adaptation to rttov 12.1
!    -------------------
!     purpose:
!
!arguments
!

! appel de RTTOV pour le calcul des radiances AVHRR
! (non assimilees mais necessaires au background check IASI)

    implicit none
    integer ,intent(in) :: idiasi
    integer ,intent (in) :: iptobs(1)
    real (8) , intent (in) :: surfem1_avhrr(3)
!*********************************************************************
    type (rttov_chanprof)  :: chanprof(3)
    logical :: calcemis  (3)
    integer ::  list_sensor (3),errorstatus,alloc_status(2)
    integer, save :: idiasi_old=-1
    integer :: ich,i,j,jn
    integer :: ichan_avhrr (NIR)
    type ( rttov_transmission )  :: transmission
    type ( rttov_radiance )      :: radiancedata_d
    type ( rttov_emissivity )    :: emissivity(3)
    integer :: nchannels
    integer :: asw,nlevels,io
!***********************************************
  
    if (IDIASI_OLD /= IDIASI) then
      LIST_SENSOR(1) = 10
      LIST_SENSOR(2) = idiasi
      LIST_SENSOR(3) = 5
      do ICH=1,NIR
        ICHAN_AVHRR(ICH)=ICH
      end do
    
      errorstatus = 0

      if (IDIASI_OLD > 0) then
        call rttov_dealloc_coefs(errorstatus, coefs_avhrr )
        if ( errorstatus /= 0) then
          write(*,*) "Probleme dans rttov_dealloc_coefs !"
          call utl_abort("tovs_rttov_AVHRR_for_IASI")
        end if
      end if

      call rttov_read_coefs ( errorstatus, &! out
           coefs_avhrr,     &! out
           tvs_opts(1),     &! in
           channels=ichan_avhrr,     &! in
           instrument=list_sensor )     ! in

       
      if ( errorstatus /= 0) then
        write(*,*) "Probleme dans rttov_read_coefs !"
        call utl_abort("tovs_rttov_AVHRR_for_IASI")
      end if
     
      IDIASI_OLD = IDIASI
   
    end if


    nlevels = coefs_avhrr % coef % nlevels
  
    nchannels = NIR

    calcemis(:) = .false.
    emissivity(1:3)%emis_in = surfem1_avhrr(1:3)
    ! Build the list of channels/profiles indices

    do  ich = 1, nchannels
      chanprof(ich) % prof = 1
      chanprof(ich) % chan = ich
    end do

    ! allocate transmittance structure
    alloc_status = 0
    allocate( transmission % tau_levels     ( nlevels, nchannels ) ,stat= alloc_status(1))
    allocate( transmission % tau_total      ( nchannels )          ,stat= alloc_status(2))
    call utl_checkAllocationStatus(alloc_status, " tovs_rttov_AVHRR_for_IASI transmission")
    transmission % tau_levels (:,:) = 0.0D0
    transmission % tau_total (:) = 0.0D0
    ! allocate radiance structure

    asw = 1 ! 1 to allocate,0 to deallocate
    call rttov_alloc_rad (alloc_status(1),nchannels,radiancedata_d,nlevels,asw)
    call utl_checkAllocationStatus(alloc_status(1:1), " tovs_rttov_AVHRR_for_IASI radiances")
 
    call rttov_direct(            &
         errorstatus,             & ! out
         chanprof,                & ! in
         tvs_opts(1),             & ! in
         tvs_profiles(iptobs(:)), & ! in
         coefs_avhrr,             & ! in
         transmission,            & ! inout
         radiancedata_d,          & ! out
         calcemis=calcemis,       & ! in
         emissivity=emissivity)     ! inout
    
    io = iptobs(1)
    avhrr_bgck(io)% RADCLEARCALC(NVIS+1:NVIS+NIR) = radiancedata_d % clear(1:NIR)
    avhrr_bgck(io)% TBCLEARCALC(NVIS+1:NVIS+NIR)  = radiancedata_d % bt(1:NIR)
    avhrr_bgck(io)% RADOVCALC(1:nlevels-1,NVIS+1:NVIS+NIR) = radiancedata_d % overcast(1:nlevels-1,1:NIR)
    avhrr_bgck(io)% TRANSMCALC(1:nlevels,NVIS+1:NVIS+NIR) =  transmission % tau_levels(1:nlevels,1:NIR)
    avhrr_bgck(io)% EMISS(NVIS+1:NVIS+NIR) = emissivity(1:NIR)%emis_out
    avhrr_bgck(io)% TRANSMSURF(NVIS+1:NVIS+NIR) = transmission% tau_total(1:NIR)


    deallocate( transmission % tau_total    ,stat= alloc_status(1))
    deallocate( transmission % tau_levels   ,stat= alloc_status(2))
    call utl_checkAllocationStatus(alloc_status, " tovs_rttov_AVHRR_for_IASI transmission", .false.)

    asw = 0 ! 1 to allocate,0 to deallocate
    call rttov_alloc_rad (alloc_status(1),nchannels,radiancedata_d,nlevels,asw)
    call utl_checkAllocationStatus(alloc_status(1:1), " tovs_rttov_AVHRR_for_IASI radiances", .false.)
  
  end subroutine tovs_rttov_AVHRR_for_IASI

  subroutine  COR_ALBEDO  ( DEL, SCOS )
!***subroutine     COR_ALBEDO
!*
!*auteur           Louis Garand  - rpn - dorval
!*
!*revision 001     Jacques Halle - ddo - dorval - 421-4660
!*                                 fev 1991
!*                 adapter au systeme operationel GOES.
!*
!*REVISION 002     JACQUES HALLE - DDO - DORVAL - 421-4660
!*                                 Decembre 1995
!*                 Generaliser pour toutes les plateformes satellitaires.
!*
!*objet            ce sous-programme calcule un facteur de correction
!*                 pour l'albedo a partir du cosinus de l'angle solaire. 
!*
!*appel            CALL COR_ALBEDO  ( DEL, SCOS )
!*
!*arguments        del   - output - facteur de correction
!*                 scos  - input  - cosinus de l'angle solaire
!**
    implicit  none
    real(8),intent(in)  ::  scos
    real(8),intent(out) ::  del
    !************************************
    integer  i1, i2
    real(8)  x1, x2, g1, g2, a, b
    real(8)  S(11)
 
    data  S / 00.00d0, 18.19d0, 31.79d0, 41.41d0, 49.46d0, 56.63d0, 63.26d0, 69.51d0, 75.52d0, 81.37d0, 87.13d0 /
 
    I1  = 12 -( SCOS + 0.05d0) * 10.d0 
    I2  = I1 + 1 
    I1  = min(I1,11)
    I2  = min(I2,11)
    X1  = cos ( S(I1) * MPC_RADIANS_PER_DEGREE_R8 )  
    X2  = cos ( S(I2) * MPC_RADIANS_PER_DEGREE_R8 ) 
    G1  = DRCLD(I1)
    G2  = DRCLD(I2)
    if (I1 == I2) then
      DEL =G1
    else
      call  SOLU ( G1, X1, G2 ,X2, A, B )
      DEL = A * SCOS + B
    end if
  
  end subroutine COR_ALBEDO


  subroutine  SOLU ( YY1, XX1, YY2, XX2, AA, BB )
!**subroutine     SOLU
!
!auteur           Louis Garand  - rpn - dorval
!
!revision 001     Jacques Halle - ddo - dorval - 421-4660
!                                 fev 1991
!                 adapter au systeme operationel GOES.
!
!REVISION 002     JACQUES HALLE - DDO - DORVAL - 421-4660
!                                 Decembre 1995
!                 Generaliser pour toutes les plateformes satellitaires.
!
!langage          fortran 90
!
!objet            ce sous-programme calcule la pente et l'intercept
!                 a partir de deux couples de donnees.
!
!appel            CALL SOLU ( Y1, X1, Y2, X2, A, B )
!
!arguments        XY1    - input - coordonnee Y du point 1
!                 XX1    - input - coordonnee X du point 1
!                 YY2    - input - coordonnee Y du point 2
!                 YX2    - input - coordonnee X du point 2
!                 AA     - output- pente
!                 BB     - output- intercept
!*
    implicit none
    real(8),intent (in)  ::  YY1, XX1, YY2, XX2
    real(8),intent (out) ::  AA, BB
    ! 
!  DROITE PASSANT PAR DEUX POINTS PENTE A ET INTERCEPT B
!


    AA = (YY1 - YY2) / (XX1 - XX2)
    BB = YY1 - AA * XX1

  end subroutine SOLU

 
  real(8) function  DRCLD ( IZ ) 
!**fonction       DRCLD
!
!auteur           Louis Garand  - rpn - dorval
!
!revision 001     Jacques Halle - ddo - dorval - 421-4660
!                                 fev 1991
!                 adapter au systeme operationel GOES.
!
!REVISION 002     JACQUES HALLE - DDO - DORVAL - 421-4660
!                                 Decembre 1995
!                 Generaliser pour toutes les plateformes satellitaires.
!
!langage          fortran 5
!
!objet            ce sous-programme calcule la normalisation due
!                 a l'angle zenith solaire selon 
!                 MINNIS-HARRISSON (COURBE FIG 7), P1038,JCAM 84.  
!
!appel            xnorm = DRCLD ( IZ )
!
!arguments        xnorm - output - facteur de normalisation
!                 iz    - input  - cosinus de l'angle solaire
!*
    implicit  none

    integer,intent (in) ::  iz
    
    real(8)  DRF(11) 
    
    data  DRF / 1.000d0, 1.002d0, 1.042d0, 1.092d0, 1.178d0, 1.286d0, &
         1.420d0, 1.546d0, 1.710d0, 1.870d0, 2.050d0  / 

    DRCLD = DRF (IZ)
    
  end function DRCLD

  subroutine VISOCN(SZ,SATZ,RZ,ANISOT,ZLAMB,ZCLOUD,IER)
!***subroutine     VISOCN
!*
!*auteur           LOUIS GARAND 1985
!*
!*REVISION 001     JACQUES HALLE - DDO - DORVAL - 421-4660
!*                                 Decembre 1995
!*                 Generaliser pour toutes les plateformes satellitaires.
!*
!*objet            THIS ROUTINE PROVIDES THE CORRECTIVE FACTORS FOR THE ANISOTROPY
!*                 OF REFLECTANCE OVER CLEAR OCEAN.
!*                 
!*
!*appel            CALL VISOCN(SZ,SATZ,RZ,ANISOT,ZLAMB,ZCLOUD,IER)
!*
!*arguments        sz     - input  - SUN ZENITH ANGLE IN DEGREES (0 TO 90)
!*                 satz   - input  - SATELLITE ZENITH ANGLE (0 TO 90)
!*                 rz     - input  - RELATIVE   ANGLE IN DEGREES (0 TO 180) WITH
!*                                   0 AS BACKSCATTERING AND 
!*                                   180 AS FORWARD SCATTERING
!*                 anisot - output - ANISOTROPIC CORRECTIVE FACTOR 
!*                                  (KHI IN MINNIS-HARRISSON)
!*                 zlamb  - output - CORRECTIVE FACTOR FOR LAMBERTIAN REFLECTANCE
!*                                   (DELTA """") ZLAMB IS A FUNCTION OF SZ ONLY.
!*                                   THIS IS FOR OCEAN SURFACE.
!*                 zcloud - output - SAME AS ZLAMB BUT FOR CLOUD SURFACE
!*                 ier    - output - error code (0=ok; -1=problem with interpolation)
!*
!*notes            OBTAINED FROM DR PAT MINNIS,LANGLEY , AND BASED ON THE WORK
!*                 OF MINNIS AND HARRISSON,JCAM 1984,P993.
!*                 THE ROUTINE IS A LOOK UP TABLE ALONG WITH INTERPOLATION ON THE 
!*                 THREE ANGLES. 
!**
    implicit  none
    real (8),intent(in) :: SZ,SATZ,rz
    real (8),intent(out):: ANISOT,ZLAMB,ZCLOUD
    integer ,intent(out) :: ier
    !********************************************************
    integer  i1, i2, j1, j2, k1, k2, l, i, n, m, j, k
    real(8) cc, d1, d2, slop, cept, x1, x2
    real(8) g1, g2
    real(8) VNORM(11,10,13),S(11),V(10),R(13),DA(2),DD(2) 

    data S/0.0d0,18.19d0,31.79d0,41.41d0,49.46d0,56.63d0,63.26d0,69.51d0,75.52d0,81.37d0,87.13d0/ 
    
    data R/0.0d0,15.0d0,30.0d0,45.0d0,60.0d0,75.0d0,90.0d0,105.0d0,120.0d0,135.0d0,150.0d0,165.0d0,180.0d0/ 

    data V/0.0d0,10.0d0,20.0d0,30.0d0,40.0d0,50.0d0,60.0d0,70.0d0,80.0d0,90.0d0/

    data ((VNORM(1,J,K),J=1,10),K=1,13)/  &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0, &
         2.668d0,2.210d0,1.105d0,0.979d0,0.810d0,0.735d0,0.785d0,0.979d0,1.092d0,1.174d0/ 
    
    data ((VNORM(2,J,K),J=1,10),K=1,13)/  &
         1.154d0, .960d0, .896d0, .818d0, .748d0, .825d0, .922d0,1.018d0,1.179d0,1.334d0, &
         1.154d0, .954d0, .838d0, .799d0, .735d0, .786d0, .883d0, .960d0,1.128d0,1.250d0, &
         1.514d0, .973d0, .825d0, .786d0, .722d0, .754d0,0.838d0,0.922d0,1.063d0,1.160d0, &
         1.514d0,0.967d0,0.864d0,0.818d0,0.715d0,0.728d0,0.793d0,0.876d0,1.005d0,1.102d0, &
         1.514d0,0.967d0,0.896d0,0.889d0,0.702d0,0.696d0,0.773d0,0.851d0,0.954d0,1.038d0, &
         1.514d0,1.070d0,0.986d0,0.922d0,0.677d0,0.696d0,0.754d0,0.838d0,0.922d0,1.012d0, &
         1.514d0,1.270d0,0.967d0,0.870d0,0.677d0,0.664d0,0.709d0,0.773d0,0.857d0,0.954d0, &
         1.514d0,1.495d0,1.166d0,0.960d0,0.683d0,0.690d0,0.728d0,0.806d0,0.896d0,0.999d0, &
         1.514d0,1.959d0,1.534d0,1.025d0,0.973d0,0.709d0,0.754d0,0.857d0,0.954d0,1.050d0, &
         1.514d0,2.165d0,2.165d0,1.270d0,1.038d0,0.760d0,0.812d0,0.902d0,1.012d0,1.115d0, &
         1.514d0,2.275d0,2.262d0,1.688d0,1.115d0,0.780d0,0.857d0,0.954d0,1.070d0,1.173d0, &
         1.514d0,2.326d0,2.520d0,2.172d0,1.257d0,0.812d0,0.883d0,1.005d0,1.108d0,1.212d0, &
         1.514d0,2.359d0,2.951d0,2.255d0,1.411d0,0.980d0,0.915d0,1.050d0,1.160d0,1.295d0/ 

    data ((VNORM(3,J,K),J=1,10),K=1,13)/   &
         0.897d0,0.792d0,0.765d0,0.765d0,0.778d0,0.897d0,0.996d0,1.095d0,1.306d0,1.431d0, &
         0.897d0,0.712d0,0.739d0,0.745d0,0.765d0,0.891d0,0.970d0,1.069d0,1.214d0,1.359d0, &
         0.897d0,0.666d0,0.699d0,0.745d0,0.759d0,0.811d0,0.917d0,1.042d0,1.148d0,1.306d0, &
         0.897d0,0.646d0,0.693d0,0.739d0,0.693d0,0.752d0,0.858d0,0.989d0,1.102d0,1.234d0, &
         0.897d0,0.686d0,0.679d0,0.726d0,0.679d0,0.693d0,0.792d0,0.924d0,1.049d0,1.154d0, &
         0.897d0,0.660d0,0.673d0,0.693d0,0.646d0,0.660d0,0.759d0,0.858d0,1.003d0,1.102d0, &
         0.897d0,0.673d0,0.765d0,0.792d0,0.712d0,0.600d0,0.699d0,0.811d0,0.963d0,1.055d0, &
         0.897d0,0.706d0,0.772d0,0.917d0,0.904d0,0.613d0,0.726d0,0.858d0,1.055d0,1.121d0, &
         0.897d0,0.825d0,0.924d0,0.996d0,0.989d0,0.686d0,0.778d0,0.937d0,1.115d0,1.181d0, &
         0.897d0,1.036d0,1.253d0,1.286d0,1.260d0,0.778d0,0.858d0,0.996d0,1.181d0,1.260d0, &
         0.897d0,1.201d0,1.788d0,1.986d0,1.827d0,0.884d0,0.851d0,1.062d0,1.227d0,1.333d0, &
         0.897d0,1.530d0,2.249d0,2.546d0,2.381d0,1.352d0,0.891d0,1.108d0,1.286d0,1.405d0, &
         0.897d0,1.854d0,2.401d0,3.325d0,2.559d0,1.590d0,0.937d0,1.168d0,1.214d0,1.425d0/ 

    data ((VNORM(4,J,K),J=1,10),K=1,13)/  &
         0.752d0,0.800d0,0.745d0,0.717d0,0.759d0,0.891d0,1.149d0,1.309d0,1.469d0,1.650d0, &
         0.752d0,0.773d0,0.717d0,0.703d0,0.752d0,0.835d0,1.065d0,1.246d0,1.406d0,1.552d0, &
         0.752d0,0.731d0,0.689d0,0.703d0,0.745d0,0.814d0,0.988d0,1.176d0,1.323d0,1.476d0, &
         0.752d0,0.689d0,0.675d0,0.654d0,0.696d0,0.752d0,0.940d0,1.100d0,1.246d0,1.378d0, &
         0.752d0,0.675d0,0.661d0,0.633d0,0.668d0,0.717d0,0.877d0,1.030d0,1.176d0,1.309d0, &
         0.752d0,0.647d0,0.640d0,0.620d0,0.613d0,0.682d0,0.814d0,0.947d0,1.107d0,1.232d0, &
         0.752d0,0.633d0,0.620d0,0.613d0,0.606d0,0.640d0,0.773d0,0.898d0,1.044d0,1.162d0, &
         0.752d0,0.626d0,0.626d0,0.626d0,0.620d0,0.654d0,0.821d0,0.947d0,1.128d0,1.225d0, &
         0.752d0,0.633d0,0.633d0,0.633d0,0.647d0,0.675d0,0.877d0,1.009d0,1.183d0,1.274d0, &
         0.752d0,0.682d0,0.717d0,0.961d0,1.023d0,0.968d0,0.940d0,1.142d0,1.274d0,1.413d0, &
         0.752d0,0.856d0,1.037d0,1.434d0,1.594d0,1.441d0,1.044d0,1.225d0,1.323d0,1.545d0, &
         0.752d0,1.044d0,1.295d0,2.207d0,1.610d0,2.311d0,1.385d0,1.274d0,1.441d0,1.636d0, &
         0.752d0,1.079d0,1.524d0,2.541d0,3.564d0,3.014d0,1.942d0,1.462d0,1.552d0,1.726d0/ 

    data ((VNORM(5,J,K),J=1,10),K=1,13)/  &
         0.552d0,0.588d0,0.617d0,0.638d0,0.724d0,0.860d0,1.133d0,1.362d0,1.556d0,1.678d0, &
         0.552d0,0.581d0,0.602d0,0.617d0,0.652d0,0.803d0,1.075d0,1.326d0,1.484d0,1.592d0, &
         0.552d0,0.559d0,0.588d0,0.595d0,0.617d0,0.731d0,1.018d0,1.283d0,1.412d0,1.527d0, &
         0.552d0,0.531d0,0.538d0,0.574d0,0.595d0,0.710d0,0.946d0,1.240d0,1.341d0,1.463d0, &
         0.552d0,0.516d0,0.523d0,0.552d0,0.559d0,0.695d0,0.911d0,1.226d0,1.291d0,1.412d0, &
         0.552d0,0.516d0,0.523d0,0.538d0,0.538d0,0.652d0,0.882d0,1.154d0,1.240d0,1.348d0, &
         0.552d0,0.516d0,0.523d0,0.538d0,0.523d0,0.595d0,0.774d0,1.075d0,1.169d0,1.269d0, &
         0.552d0,0.531d0,0.545d0,0.552d0,0.566d0,0.609d0,0.817d0,1.140d0,1.248d0,1.369d0, &
         0.552d0,0.538d0,0.545d0,0.566d0,0.581d0,0.645d0,0.911d0,1.240d0,1.319d0,1.441d0, &
         0.552d0,0.566d0,0.552d0,0.574d0,0.710d0,0.839d0,0.982d0,1.298d0,1.391d0,2.323d0, &
         0.552d0,0.566d0,0.559d0,0.710d0,1.147d0,1.176d0,1.040d0,1.348d0,1.671d0,2.674d0, &
         0.552d0,0.588d0,1.133d0,1.355d0,2.194d0,2.803d0,2.201d0,2.459d0,2.904d0,3.126d0, &
         0.552d0,0.710d0,1.341d0,1.757d0,3.026d0,3.900d0,4.445d0,4.503d0,4.445d0,4.503d0/ 

    data ((VNORM(6,J,K),J=1,10),K=1,13)/  &
         0.551d0,0.627d0,0.665d0,0.734d0,0.826d0,0.971d0,1.231d0,1.537d0,1.721d0,1.866d0, &
         0.551d0,0.604d0,0.619d0,0.665d0,0.765d0,0.895d0,1.185d0,1.476d0,1.568d0,1.652d0, &
         0.551d0,0.597d0,0.604d0,0.619d0,0.734d0,0.849d0,1.101d0,1.346d0,1.453d0,1.568d0, &
         0.551d0,0.581d0,0.589d0,0.597d0,0.665d0,0.795d0,1.032d0,1.262d0,1.346d0,1.445d0, &
         0.551d0,0.558d0,0.558d0,0.566d0,0.612d0,0.727d0,0.987d0,1.201d0,1.262d0,1.399d0, &
         0.551d0,0.505d0,0.505d0,0.512d0,0.566d0,0.696d0,0.925d0,1.117d0,1.185d0,1.308d0, &
         0.551d0,0.474d0,0.497d0,0.512d0,0.535d0,0.673d0,0.864d0,1.048d0,1.124d0,1.216d0, &
         0.551d0,0.497d0,0.505d0,0.520d0,0.551d0,0.681d0,0.902d0,1.124d0,1.201d0,1.323d0, &
         0.551d0,0.535d0,0.535d0,0.551d0,0.566d0,0.711d0,1.017d0,1.201d0,1.269d0,1.422d0, &
         0.551d0,0.535d0,0.543d0,0.558d0,0.704d0,1.193d0,1.247d0,1.285d0,1.346d0,1.950d0, &
         0.551d0,0.543d0,0.551d0,0.581d0,0.994d0,1.545d0,1.583d0,1.354d0,2.019d0,2.883d0, &
         0.551d0,0.566d0,0.612d0,0.788d0,1.468d0,2.233d0,2.340d0,2.531d0,2.983d0,3.365d0, &
         0.551d0,0.658d0,0.665d0,1.101d0,2.134d0,3.120d0,4.221d0,4.856d0,4.956d0,5.613d0/ 

    data ((VNORM(7,J,K),J=1,10),K=1,13)/  &
         0.545d0,0.606d0,0.683d0,0.744d0,0.798d0,0.990d0,1.228d0,1.704d0,1.850d0,2.049d0, &
         0.545d0,0.576d0,0.583d0,0.714d0,0.783d0,0.952d0,1.144d0,1.573d0,1.758d0,1.888d0, &
         0.545d0,0.560d0,0.568d0,0.629d0,0.744d0,0.875d0,1.105d0,1.504d0,1.642d0,1.788d0, &
         0.545d0,0.553d0,0.560d0,0.599d0,0.629d0,0.791d0,1.028d0,1.420d0,1.527d0,1.696d0, &
         0.545d0,0.545d0,0.553d0,0.599d0,0.606d0,0.714d0,0.990d0,1.335d0,1.451d0,1.581d0, &
         0.545d0,0.530d0,0.537d0,0.568d0,0.583d0,0.683d0,0.890d0,1.243d0,1.351d0,1.489d0, &
         0.545d0,0.491d0,0.499d0,0.507d0,0.576d0,0.622d0,0.791d0,1.182d0,1.282d0,1.389d0, &
         0.545d0,0.507d0,0.514d0,0.507d0,0.576d0,0.675d0,0.890d0,1.197d0,1.328d0,1.451d0, &
         0.545d0,0.522d0,0.537d0,0.522d0,0.591d0,0.760d0,0.944d0,1.259d0,1.389d0,1.527d0, &
         0.545d0,0.537d0,0.545d0,0.553d0,0.614d0,0.906d0,1.028d0,1.389d0,1.504d0,2.533d0, &
         0.545d0,0.553d0,0.553d0,0.576d0,0.637d0,1.036d0,1.550d0,1.658d0,1.934d0,3.277d0, &
         0.545d0,0.560d0,0.568d0,0.606d0,1.174d0,1.781d0,2.563d0,3.170d0,3.791d0,4.966d0, &
         0.545d0,0.591d0,0.614d0,1.259d0,2.065d0,2.824d0,3.761d0,4.498d0,5.902d0,6.148d0/ 

    data ((VNORM(8,J,K),J=1,10),K=1,13)/  &
         0.514d0,0.539d0,0.596d0,0.694d0,0.832d0,1.004d0,1.444d0,1.869d0,2.203d0,2.538d0, &
         0.514d0,0.539d0,0.571d0,0.645d0,0.751d0,0.906d0,1.387d0,1.779d0,2.056d0,2.317d0, &
         0.514d0,0.547d0,0.555d0,0.612d0,0.702d0,0.824d0,1.281d0,1.681d0,1.934d0,2.203d0, &
         0.514d0,0.539d0,0.555d0,0.588d0,0.653d0,0.743d0,1.028d0,1.404d0,1.624d0,2.024d0, &
         0.514d0,0.539d0,0.547d0,0.555d0,0.588d0,0.710d0,0.889d0,1.191d0,1.420d0,1.820d0, &
         0.514d0,0.522d0,0.522d0,0.539d0,0.563d0,0.710d0,0.849d0,1.044d0,1.208d0,1.534d0, &
         0.514d0,0.481d0,0.506d0,0.514d0,0.539d0,0.694d0,0.824d0,1.028d0,1.200d0,1.371d0, &
         0.514d0,0.481d0,0.514d0,0.547d0,0.563d0,0.702d0,0.898d0,1.134d0,1.297d0,1.501d0, &
         0.514d0,0.490d0,0.514d0,0.555d0,0.588d0,0.726d0,0.955d0,1.265d0,1.379d0,1.648d0, &
         0.514d0,0.547d0,0.547d0,0.571d0,0.604d0,0.767d0,1.036d0,1.355d0,1.550d0,3.142d0, &
         0.514d0,0.563d0,0.579d0,0.604d0,0.612d0,0.832d0,1.909d0,2.848d0,3.917d0,4.790d0, &
         0.514d0,0.522d0,0.563d0,0.677d0,0.767d0,1.420d0,2.040d0,3.158d0,4.863d0,6.291d0, &
         0.514d0,0.588d0,0.588d0,0.612d0,0.824d0,2.032d0,3.109d0,4.969d0,6.846d0,7.695d0/ 

    data ((VNORM(9,J,K),J=1,10),K=1,13)/  &
         0.572d0,0.608d0,0.679d0,0.751d0,0.831d0,1.001d0,1.377d0,1.913d0,2.512d0,2.879d0, &
         0.572d0,0.572d0,0.608d0,0.679d0,0.760d0,0.930d0,1.243d0,1.707d0,2.369d0,2.700d0, &
         0.572d0,0.563d0,0.590d0,0.644d0,0.706d0,0.831d0,1.171d0,1.618d0,2.190d0,2.378d0, &
         0.572d0,0.554d0,0.563d0,0.599d0,0.662d0,0.760d0,1.010d0,1.502d0,2.011d0,2.235d0, &
         0.572d0,0.545d0,0.563d0,0.590d0,0.626d0,0.715d0,0.885d0,1.323d0,1.815d0,2.119d0, &
         0.572d0,0.527d0,0.554d0,0.572d0,0.608d0,0.670d0,0.724d0,1.144d0,1.618d0,1.868d0, &
         0.572d0,0.545d0,0.572d0,0.572d0,0.599d0,0.662d0,0.724d0,1.117d0,1.484d0,1.761d0, &
         0.572d0,0.554d0,0.590d0,0.599d0,0.608d0,0.679d0,0.760d0,1.216d0,1.582d0,1.922d0, &
         0.572d0,0.572d0,0.599d0,0.608d0,0.635d0,0.715d0,0.822d0,1.377d0,1.707d0,2.056d0, &
         0.572d0,0.590d0,0.608d0,0.635d0,0.662d0,0.742d0,0.912d0,1.529d0,3.075d0,4.693d0, &
         0.572d0,0.590d0,0.626d0,0.644d0,0.670d0,0.760d0,1.109d0,1.564d0,3.111d0,4.702d0, &
         0.572d0,0.599d0,0.644d0,0.662d0,0.688d0,0.822d0,1.788d0,2.816d0,5.346d0,7.295d0, &
         0.572d0,0.608d0,0.662d0,0.670d0,0.715d0,1.851d0,3.227d0,4.810d0,6.669d0,9.557d0/ 
    
    data ((VNORM(10,J,K),J=1,10),K=1,13)/   &
         0.552d0,0.606d0,0.639d0,0.671d0,0.704d0,0.899d0,1.223d0,2.479d0,3.194d0,3.573d0, &
         0.552d0,0.574d0,0.606d0,0.628d0,0.682d0,0.855d0,1.148d0,2.339d0,2.642d0,3.378d0, &
         0.552d0,0.563d0,0.552d0,0.595d0,0.639d0,0.834d0,1.061d0,2.014d0,2.404d0,2.891d0, &
         0.552d0,0.563d0,0.509d0,0.552d0,0.628d0,0.801d0,0.985d0,1.689d0,2.176d0,2.653d0, &
         0.552d0,0.574d0,0.509d0,0.520d0,0.585d0,0.747d0,0.888d0,1.332d0,1.970d0,2.458d0, &
         0.552d0,0.531d0,0.509d0,0.509d0,0.531d0,0.682d0,0.801d0,1.191d0,1.819d0,2.425d0, &
         0.552d0,0.498d0,0.498d0,0.498d0,0.520d0,0.639d0,0.747d0,1.126d0,1.711d0,2.317d0, &
         0.552d0,0.498d0,0.509d0,0.509d0,0.541d0,0.671d0,0.780d0,1.278d0,1.862d0,2.598d0, &
         0.552d0,0.498d0,0.509d0,0.520d0,0.574d0,0.693d0,0.812d0,1.602d0,2.035d0,2.793d0, &
         0.552d0,0.520d0,0.520d0,0.531d0,0.595d0,0.725d0,0.844d0,1.916d0,2.588d0,3.768d0, &
         0.552d0,0.531d0,0.541d0,0.574d0,0.628d0,0.780d0,1.039d0,2.349d0,3.313d0,5.652d0, &
         0.552d0,0.574d0,0.563d0,0.606d0,0.660d0,0.812d0,1.797d0,3.010d0,5.478d0,7.492d0, &
         0.552d0,0.650d0,0.671d0,0.704d0,0.801d0,1.029d0,2.436d0,3.465d0,7.828d0,10.578d0/

    data ((VNORM(11,J,K),J=1,10),K=1,13)/   &
         0.518d0,0.576d0,0.605d0,0.633d0,0.662d0,0.864d0,1.238d0,2.620d0,3.455d0,3.887d0, &
         0.518d0,0.547d0,0.576d0,0.576d0,0.633d0,0.835d0,1.123d0,2.447d0,2.821d0,3.656d0, &
         0.518d0,0.518d0,0.518d0,0.547d0,0.605d0,0.806d0,1.036d0,2.102d0,2.533d0,3.080d0, &
         0.518d0,0.518d0,0.461d0,0.518d0,0.576d0,0.777d0,0.950d0,1.727d0,2.274d0,2.821d0, &
         0.518d0,0.547d0,0.461d0,0.489d0,0.547d0,0.720d0,0.864d0,1.353d0,2.044d0,2.591d0, &
         0.518d0,0.489d0,0.461d0,0.461d0,0.489d0,0.662d0,0.777d0,1.180d0,1.871d0,2.562d0, &
         0.518d0,0.461d0,0.461d0,0.461d0,0.489d0,0.605d0,0.720d0,1.123d0,1.756d0,2.418d0, &
         0.518d0,0.461d0,0.461d0,0.461d0,0.518d0,0.633d0,0.749d0,1.296d0,1.929d0,2.764d0, &
         0.518d0,0.461d0,0.461d0,0.489d0,0.547d0,0.662d0,0.777d0,1.641d0,2.130d0,2.994d0, &
         0.518d0,0.489d0,0.489d0,0.489d0,0.547d0,0.691d0,0.806d0,1.986d0,2.735d0,4.117d0, &
         0.518d0,0.489d0,0.489d0,0.547d0,0.576d0,0.749d0,1.008d0,2.476d0,3.599d0,6.334d0, &
         0.518d0,0.547d0,0.518d0,0.576d0,0.633d0,0.777d0,1.842d0,3.224d0,6.132d0,8.550d0, &
         0.518d0,0.605d0,0.633d0,0.662d0,0.777d0,1.008d0,2.562d0,3.771d0,8.953d0,12.293d0/
 
    !   COMPUTE SUN ZENITH BIN
    CC  = cos( SZ * MPC_RADIANS_PER_DEGREE_R8)
    I1  = 12.d0 - (CC + 0.05d0) * 10.d0
    I2  = I1 + 1 
    if (I1 >= 11) I1 = 11 
    if (I1 == 11) I2 = I1 

    !  COMPUTE SAT ZENITH BIN 
    J1  = int(SATZ / 10.d0) + 1 
    J2  = J1 + 1 
    if (J1 == 10) J2 = J1 

    !  COMPUTE RELATIVE AZIMUTH BIN 
    K1  = RZ / 15.d0 + 1.d0
    K2  = K1 + 1 
    if (K1 == 13) K2 = K1 

    !  INTERPOLATE
    IER = 0 
    do L=I1,I2  
      I = L -I1 + 1
       
    !     BETWEEN R'S FOR CONSTANT S
      do N=K1,K2 

!        BETWEEN V'S FOR CONSTANT R AND S 
        M  = N - K1 + 1
        D1 = VNORM(L,J1,N)
        D2 = VNORM(L,J2,N)
        if (D1 == D2) then
          DA(M) = D1
        else
          call LINEQ(V(J1),V(J2),D1,D2,SLOP,CEPT,IER) 
          DA(M) = SLOP * SATZ + CEPT
        end if
      end do
      if(K1 == K2) then 
        DD(I)  = DA(1) 
      else 
        call LINEQ(R(K1),R(K2),DA(1),DA(2),SLOP,CEPT,IER) 
        DD(I) = SLOP * RZ + CEPT
      end if
    end do

    !  BETWEEN S'S USING RESULT OF OTHER INTERPOLATIONS 
    if(I1 == I2) then
      ZLAMB  = DRM(I1) 
      ZCLOUD = DRCLD(I1)
      ANISOT = DD(1)
    else
      X1 = cos(S(I1) * MPC_RADIANS_PER_DEGREE_R8) 
      X2 = cos(S(I2) * MPC_RADIANS_PER_DEGREE_R8) 
      call LINEQ(X1,X2,DD(1),DD(2),SLOP,CEPT,IER) 
      ANISOT = SLOP * CC + CEPT 
      G1 = DRM(I1)
      G2 = DRM(I2)
      call LINEQ(X1,X2,G1,G2,SLOP,CEPT,IER) 
      ZLAMB  = SLOP * CC + CEPT
      G1 = DRCLD(I1)
      G2 = DRCLD(I2)
      call LINEQ(X1,X2,G1,G2,SLOP,CEPT,IER) 
      ZCLOUD = SLOP * CC + CEPT 
    end if
    
    if (ANISOT < 0.) then 
      IER = -1
      ANISOT = 1.d0 
      ZLAMB  = DRM(I1) 
      ZCLOUD = DRCLD(I1)
    end if
    
  end subroutine VISOCN

  
  subroutine LINEQ(XX1,XX2,YY1,YY2,AA,BB,IERR) 
!**subroutine     LINEQ
!
!auteur           Louis Garand  - rpn - dorval
!
!REVISION 001     JACQUES HALLE - DDO - DORVAL - 421-4660
!                                 Decembre 1995
!                Generaliser pour toutes les plateformes satellitaires.
!
!langage          fortran 90
!
!objet            calculate slope and intercept of a line.
!
!appel            CALL  LINEQ(X1,X2,Y1,Y2,A,B,IER)
!
!arguments        x1  - input  - coordinate x of point 1
!                 x2  - input  - coordinate x of point 2
!                 y1  - input  - coordinate y of point 1
!                 y2  - input  - coordinate y of point 2
!                 a   - output - slope
!                 b   - output - intercept
!                 ier - output - error code (0=ok)
!
    real(8) ,intent(in)     :: XX1,XX2,YY1,YY2
    real(8) ,intent(out)    :: AA,BB
    integer ,intent(out) :: ierr
!****************************************
     
    ierr = 0
    
    if ( (XX2 - XX1) == 0.d0) then 
      IERR = -1
      return
    end if

    AA = ( YY2 - YY1) / (XX2 - XX1) 
    BB = YY1 - AA * XX1 
    
  end subroutine LINEQ

  real(8) function DRM(IZ) 
!**function       DRM
!
!auteur           Louis Garand  - rpn - dorval
!
!REVISION 001     JACQUES HALLE - DDO - DORVAL - 421-4660
!                                 Decembre 1995
!                 Generaliser pour toutes les plateformes satellitaires.
!
!langage          fortran 90
!
!objet            NORMALIZATION FOR SUN ZENITH ANGLE (LAMBERTIAN)
!                 FOR OCEAN.
!
!appel            val = DRM(IZ)
!
!arguments        iz  - input  - index
!                 val - output - normalization factor
!*

    integer,intent (in) ::  iz

    real(8)  DRF(11)

    data DRF /1.d0,1.0255d0,1.1197d0,1.2026d0,1.3472d0,1.4926d0,1.8180d0,2.1980d0, &
         2.8180d0,3.8615d0,4.3555d0/

    DRM = DRF(IZ) 
  
  end function DRM
      

end module MULTI_IR_BGCK_MOD
