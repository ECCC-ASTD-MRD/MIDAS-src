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
!! MODULE computeHBHT_mod (prefix='hbht' category='1. High-level functionality')
!!
!! *Purpose*: Contains subroutines for computing the background error 
!!            variance in observation space
!!
!--------------------------------------------------------------------------
module computeHBHT_mod
  use mpi_mod
  use mpivar_mod
  use obsSpaceData_mod
  use columnData_mod
  use bufr_mod
  use utilities_mod
  use EarthConstants_mod
  use MathPhysConstants_mod
  use stateToColumn_mod
  use gridStateVector_mod
  use verticalCoord_mod
  use horizontalCoord_mod
  use analysisGrid_mod
  use bmatrixhi_mod
  use obsOperators_mod
  use gps_mod
  use codtyp_mod
  use bmatrixchem_mod
  use chem_obsoperators_mod
  use timeCoord_mod
  use obsTimeInterp_mod
  use bMatrixEnsemble_mod
  use varNameList_mod
  implicit none
  private

  ! public procedures
  public :: hbht_compute, hbht_compute_static, hbht_compute_static_chem, hbht_compute_ensemble


  contains

!--------------------------------------------------------------------------
!! *Purpose*: Compute background error stddev in observation space.
!!
!! Revision:   
!!v            M. Sitwell (ARQI/AQRD)  May 2015
!!v               - Added call to compute_HBHT_static_chem for computation
!!v                 with B_chm for chemical constituents
!!v            Y.J. Rochon (ARQI/AQDR) March 2016
!!v               - Allowed static=ensemble=.false. when static_chm=.true.
!!v               - Allowed possibility of static_chm=.true. and ensemble=.true.
!!v                 with and without static=.true.
!--------------------------------------------------------------------------
subroutine hbht_compute(columng,columnhr,obsSpaceData)

  implicit none

  type(struct_obs)        :: obsSpaceData ! Observation-related data
  type(struct_columnData) :: columng      ! Columns of the background interpolated 
                                          ! to analysis levels and to obs horizontal locations
  type(struct_columnData) :: columnhr     ! Columns of the background interpolated 
                                          ! to obs horizontal locations

  real(8) :: HBHT_static, HBHT_ensemble, HBHT_hybrid

  logical :: static   = .false.
  logical :: ensemble = .false.
  logical :: static_chm = .false.

  integer :: index_body
  integer :: fnom, fclos, ierr, nulnam

  character(len=12) :: hybrid_mode
  character(len=2)  :: fam

  !namelist
  NAMELIST /NAMHBHT/hybrid_mode

  !
  !- 1.  Compute HBHT (sigma_b in observation space)
  !

  !- 1.1 HBHT from the Bnmc matrix
  call hbht_compute_static( columng,      & ! IN
                            columnhr,     & ! IN
                            obsSpaceData, & ! INOUT (HBnmcHT std. dev. outputted in OBS_HPHT)
                            static        ) ! OUT   (Active if TRUE)

  !- 1.2 HBHT from the Bens
  call hbht_compute_ensemble( columng,      & ! IN
                              columnhr,     & ! IN
                              obsSpaceData, & ! INOUT (HBensHT std. dev. outputted in OBS_WORK)
                              ensemble )      ! OUT   (Active if TRUE)

  !- 1.3 HBHT from the B_static matrix for chemistry
  call hbht_compute_static_chem( columng,      & ! IN
                                 obsSpaceData, & ! INOUT ( sqrt(diag(H*B*H^T)) with B_static_chm outputted in OBS_HPHT )
                                 static_chm )    ! OUT   (Active if TRUE)

  !
  !- 2. Select/Blend HBHT
  !
  if ( (static.or.static_chm) .and. .not. ensemble ) then
     ! Bnmc only
     write(*,*)
     write(*,*) 'compute_HBHT: Using B_static ONLY'
     ! HBnmcHT std. dev. already in OBS_HPHT
  else if ( .not. static .and. ensemble .and. .not. static_chm ) then
     write(*,*)
     write(*,*) 'compute_HBHT: Using B_ensemble ONLY'
     ! Transfer HBensHT std. dev. values in OBS_WORK to OBS_HPHT
     do index_body = 1, obs_numBody(obsSpaceData)
        HBHT_ensemble = obs_bodyElem_r(obsSpaceData,OBS_WORK,index_body)
        call obs_bodySet_r(obsSpaceData,OBS_HPHT,index_body, HBHT_ensemble)
     end do
  else if ( (static.or.static_chm) .and. ensemble ) then
     ! Read Namelist first
     hybrid_mode = 'WEIGHTED_SUM' ! default value
     nulnam = 0
     ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
     read(nulnam,nml=namhbht,iostat=ierr)
     if ( ierr /= 0) call utl_abort('compute_HBHT: Error reading namelist')
     if ( mpi_myid == 0 ) write(*,nml=namhbht)
     ierr = fclos(nulnam)

     write(*,*)
     write(*,*) 'compute_HBHT: Using hybrid approach (blend of B_static and B_ensemble) in mode = ', trim(hybrid_mode)

     do index_body = 1, obs_numBody(obsSpaceData)
        fam = obs_getFamily(obsSpaceData,obs_bodyElem_i(obsSpaceData,OBS_HIND,index_body))
        if ( (trim(fam).eq.'CH'.and..not.static_chm) .or. (trim(fam).ne.'CH'.and..not.static) ) then
           HBHT_static = 0.0D0
        else
           HBHT_static = obs_bodyElem_r(obsSpaceData,OBS_HPHT,index_body)
        end if
        HBHT_ensemble = obs_bodyElem_r(obsSpaceData,OBS_WORK,index_body)
        select case ( trim(hybrid_mode) )
        case ('WEIGHTED_SUM')
           HBHT_hybrid = sqrt(HBHT_static**2 + HBHT_ensemble**2)
        case ('MAX_VALUE')
           HBHT_hybrid = max(HBHT_static,HBHT_ensemble)
        case default
           write(*,*)
           write(*,*) 'compute_HBHT: Unknown hybrid_mode ', trim(hybrid_mode)
           call utl_abort('compute_HBHT')
        end select
        call obs_bodySet_r(obsSpaceData,OBS_HPHT,index_body, HBHT_hybrid)
     end do

  else
     call utl_abort('compute_HBHT: no B matrix was initialized')
  end if

end subroutine hbht_compute


!--------------------------------------------------------------------------
!! *Purpose*: Compute background error stddev in observation space using
!!            fixed statistics specific in stats file.
!!
!--------------------------------------------------------------------------
SUBROUTINE hbht_compute_static(lcolumng,lcolumnhr,lobsSpaceData,active)
      IMPLICIT NONE

      type(struct_hco), pointer        :: hco_anl
      type(struct_vco), pointer        :: vco_anl
      type(struct_obs)        :: lobsSpaceData
      type(struct_columnData) :: lcolumn,lcolumng,lcolumnhr
      type(struct_gsv)        :: statevector

      logical                 :: active

      INTEGER JLAT, JLON, JLEV, JOBS

      CHARACTER*12 CLETIKET
      CHARACTER*2 CLTYPVAR
      CHARACTER*1 CLGRTYP
      CHARACTER*4 CLNOMVAR

      INTEGER IULSSF,IDATEO
      INTEGER FSTPRM,FNOM,FSTOUV,FCLOS,FSTFRM
      INTEGER IKEY,ILEN,IERR,IDATE

      REAL*8, allocatable :: ZBUFFER(:,:)
      real*8, pointer     :: gz_column(:), tt_column(:), field_ptr(:,:,:)

      INTEGER INI,INJ,INK, INPAS, INBITS, IDATYP, IDEET
      INTEGER IP1,IP2,IP3,IG1,IG2,IG3,IG4,ISWA,ILENGTH,IDLTF
      INTEGER IUBC,IEXTR1,IEXTR2,IEXTR3

      integer :: nLev_M,nLev_T,status,shift_level,Vcode_anl,cvdim

      real(8), allocatable  :: scaleFactor(:)

      logical             :: is_staggered

      !- Get the appropriate Horizontal and Vertical Coordinate
      hco_anl => agd_getHco('ComputationalGrid')
      vco_anl => col_getVco(lcolumng)

      call bhi_setup( hco_anl,vco_anl, &  ! IN
                      cvdim,           &  ! OUT
                      'BackgroundCheck' ) ! IN

      if ( cvdim > 0 ) then
         write(*,*)
         write(*,*) 'Computing HBHT from Bnmc - Start'
         active = .true.
      else
         if ( mpi_myid == 0 ) write(*,*) 'compute_HBHT_static: option NOT ACTIVATED'
         active = .false.
         return
      end if

      nlev_T = col_getNumLev(LCOLUMNG,'TH')
      nlev_M = col_getNumLev(LCOLUMNG,'MM')

      allocate(scaleFactor(max(nLev_M,nLev_T)))
      call bhi_getScaleFactor(scaleFactor)

      status = vgd_get(vco_anl%vgrid,key='ig_1 - vertical coord code',value=Vcode_anl)
      if ( Vcode_anl == 5001 ) then
         is_staggered = .false.
      else if ( Vcode_anl == 5002 ) then
         is_staggered = .true.
      else
         write(*,*) 'Vcode_anl = ',Vcode_anl
         call utl_abort('compute_HBHT_static: unknown vertical coordinate type!')
      end if
      if ( mpi_myid == 0 ) write(*,*) 'compute_HBHT: vertical coord is_staggered = ',is_staggered

      if ( is_staggered ) then
         shift_level = 1
      else
         shift_level = 0
      endif

      allocate(ZBUFFER(HCO_ANL%NI,HCO_ANL%NJ))

      call gsv_allocate(statevector, 1, hco_anl, vco_anl)
      call gsv_zero(statevector)

      call col_setVco(lcolumn,col_getVco(lcolumng))
      call col_allocate(lcolumn,col_getNumCol(lcolumng))
!
!     Set the value of OBS_LYR required by setfge routines
!
      call oop_vobslyrs(lcolumng,lobsSpaceData)
!
!     1. Opening the statistics file
!
      IULSSF=0
      IERR=FNOM(iulssf,'./bgcov','RND+OLD+R/O',0)
      IF ( IERR .EQ. 0 ) THEN
        write(*,*) 'IBGST - File : ./bgcov'
        write(*,*) ' opened as unit file ',iulssf
        ierr =  fstouv(iulssf,'RND+OLD')
      ELSE
        CALL utl_abort('HBHT_static:NO BACKGROUND STAT FILE!!')
      ENDIF
!
!     2.1 Background error standard deviations
!
      CLETIKET = 'BGCK_STDDEV'
      write(*,*) 'HBHT_static: CLETIKET = ',CLETIKET
      IDATE    = -1
      IP2      = -1
      IP3      = -1
      CLTYPVAR =' '
!
!     READ IN STANDARD DEVIATION FOR EACH OBSERVATION TYPE
!
      clnomvar = 'UU'
      write(*,*) clnomvar 
      field_ptr => gsv_getField3D_r8(statevector,'UU')       
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
            call utl_abort('compute_HBHT_static')
         end if
         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scaleFactor(jlev+shift_level)*zbuffer(jlon,jlat)*MPC_M_PER_S_PER_KNOT_R8
             end do
          end do
      end do

      clnomvar = 'VV'
      write(*,*) clnomvar
      field_ptr => gsv_getField3D_r8(statevector,'VV')
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
            call utl_abort('compute_HBHT_static')
         end if

         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scalefactor(jlev+shift_level)*zbuffer(jlon,jlat)*MPC_M_PER_S_PER_KNOT_R8
            end do
         end do
      end do

      clnomvar = 'ES'
      write(*,*) clnomvar
      field_ptr => gsv_getField3D_r8(statevector,'HU')
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
            call utl_abort('compute_HBHT_static')
         end if
         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scaleFactor(jlev)*zbuffer(jlon,jlat)
            end do
         end do
      end do

      ! GZ is put into TT slot in gridStateVector
      clnomvar = 'GZ'
      write(*,*) clnomvar
      field_ptr => gsv_getField3D_r8(statevector,'TT')
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
            call utl_abort('compute_HBHT_static')
         end if
         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scaleFactor(jlev)*zbuffer(jlon,jlat)*RG*10.d0
            end do
         end do
      end do

      call s2c_bgcheck_bilin(lcolumn,statevector,lobsSpaceData)

      ! copy GZ data from TT to GZ slot in columnData
      do jobs= 1, col_getNumCol(lcolumn)
         gz_column => col_getColumn(lcolumn,jobs,'GZ','TH')
         tt_column => col_getColumn(lcolumn,jobs,'TT')
         do jlev = 1,col_getNumLev(lcolumn,'TH')
            gz_column(jlev)=tt_column(jlev)
         enddo
      enddo

!
!     SET THE FIRST-GUESS ERRORS FOR CONVENTIONAL DATA ON PRESSURE LEVELS
!     --------------------------------------------------------------------
!
      call setfgefam('AI',lcolumn,lcolumng,lobsSpaceData)
      call setfgefam('SW',lcolumn,lcolumng,lobsSpaceData)
      call setfgefam('UA',lcolumn,lcolumng,lobsSpaceData)
      call setfgefam('SF',lcolumn,lcolumng,lobsSpaceData)
      call setfgefam('HU',lcolumn,lcolumng,lobsSpaceData)
      call setfgefamz('PR',lcolumn,lcolumng,lobsSpaceData)
      call setfgefamz('AL',lcolumn,lcolumng,lobsSpaceData)
!
!     SET THE FIRST-GUESS ERRORS FOR RADIO OCCULTATION DATA
!     -----------------------------------------------------
!
      call setfgedif('RO',lcolumng,lobsSpaceData)

!
!     DO TEMPERATURE FIRST-GUESS ERROR
!     ---------------------------------
!
      clnomvar = 'TT'
      write(*,*) clnomvar
      field_ptr => gsv_getField3D_r8(statevector,'TT')
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
            call utl_abort('compute_HBHT_static')
         end if
         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scaleFactor(jlev)*zbuffer(jlon,jlat)
            end do
         end do
      end do

      call s2c_bgcheck_bilin(lcolumn,statevector,lobsSpaceData)
      call setfgett(lcolumn,lcolumng,lobsSpaceData)

!
!     RELOAD DATA TO DO SURFACE FIRST-GUESS ERRORS
!     ---------------------------------------------
!
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
         call utl_abort('compute_HBHT_static')
      end if
      field_ptr => gsv_getField3D_r8(statevector,'P0')
      do jlev = 1, ink
         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scaleFactor(max(nLev_M,nLev_T))*zbuffer(jlon,jlat)*MPC_PA_PER_MBAR_R8
            end do
         end do
      end do

      clnomvar = 'UU'
      write(*,*) clnomvar
      field_ptr => gsv_getField3D_r8(statevector,'UU')
      do jlev = 1, nlev_M
         ip1 = vco_anl%ip1_M(jlev)
         ikey = utl_fstlir(zbuffer,iulssf,ini,inj,ink,idate,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scaleFactor(jlev+shift_level)*zbuffer(jlon,jlat)*MPC_M_PER_S_PER_KNOT_R8
            end do
         end do
      end do

      clnomvar = 'VV'
      write(*,*) clnomvar
      field_ptr => gsv_getField3D_r8(statevector,'VV')
      do jlev = 1, nlev_M
         ip1 = vco_anl%ip1_M(jlev)
         ikey = utl_fstlir(zbuffer,iulssf,ini,inj,ink,idate,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scaleFactor(jlev+shift_level)*zbuffer(jlon,jlat)*MPC_M_PER_S_PER_KNOT_R8
            end do
         end do
      end do

      clnomvar = 'TT'
      write(*,*) clnomvar
      field_ptr => gsv_getField3D_r8(statevector,'TT')
      do jlev = 1, nlev_T
         ip1 = vco_anl%ip1_T(jlev)
         ikey = utl_fstlir(zbuffer,iulssf,ini,inj,ink,idate,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scaleFactor(jlev)*zbuffer(jlon,jlat)
            end do
         end do
      end do

      clnomvar = 'ES'
      write(*,*) clnomvar
      field_ptr => gsv_getField3D_r8(statevector,'HU')
      do jlev = 1, nlev_T
         ip1 = vco_anl%ip1_T(jlev)
         ikey = utl_fstlir(zbuffer,iulssf,ini,inj,ink,idate,cletiket,ip1,ip2,ip3,cltypvar,clnomvar)
         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scaleFactor(jlev)*zbuffer(jlon,jlat)
            end do
         end do
      end do

      call s2c_bgcheck_bilin(lcolumn,statevector,lobsSpaceData)
!
!     SET THE FIRST-GUESS ERRORS FOR THE SURFACE DATA
!     ------------------------------------------------
!
      call setfgesurf(lcolumn,lcolumng,lobsSpaceData)

!     READ IN LN Q FIRST-GUESS ERRORS FOR SETFGEGPS
!     ---------------------------------------------
!
      clnomvar = 'LQ'
      write(*,*) clnomvar
      field_ptr => gsv_getField3D_r8(statevector,'HU')
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
            call utl_abort('compute_HBHT_static')
         end if
         do jlat = 1, hco_anl%nj
            do jlon = 1, hco_anl%ni
               field_ptr(jlon,jlat,jlev) = scaleFactor(jlev)*zbuffer(jlon,jlat)
            end do
         end do
      end do

      call s2c_bgcheck_bilin(lcolumn,statevector,lobsSpaceData)
!
!     OPTIONAL TEST OF THE GB-GPS ZTD OPERATOR JACOBIAN
!     -------------------------------------------------
!
      if (ltestop) call setfgegps(lcolumn,lcolumng,lobsSpaceData)

      deallocate(scaleFactor)
      call col_deallocate(lcolumn)
      deallocate(zbuffer)

      write(*,*)
      write(*,*) 'Computing HBHT from Bnmc - END'
      
    END SUBROUTINE hbht_compute_static


SUBROUTINE hbht_compute_static_chem(lcolumng,lobsSpaceData,active)
!
! Author  : M. Sitwell, May 2015
!
! Purpose: To compute the background error standard deviations in
!          observation space, sqrt(diag(H*B_static*H^T)).
!
! Arguments:
!
!  Input
!
!           lcolumng             column at observation location
!
!  Inout:
!
!           lobsSpaceData        observation space data, output saved in OBS_HPHT column
!
!  Output:
!           active               flag to indicate if chemical consituents are to be used
!
! Revision:
!
!-----------------------------------------------------------------------------------------

  implicit none
  
  type(struct_hco), pointer :: hco_anl
  type(struct_vco), pointer :: vco_anl
  type(struct_obs)        :: lobsSpaceData
  type(struct_columnData) :: lcolumng
  logical                 :: active
      
  integer :: cvdim

  
  !- Get the appropriate Horizontal and Vertical Coordinate
  hco_anl => agd_getHco('ComputationalGrid')
  vco_anl => col_getVco(lcolumng)
  
  call bchm_setup( hco_anl,vco_anl, &  ! IN
                   cvdim, &            ! OUT
                  'BackgroundCheck' )  ! IN

  active = bchm_is_initialized()
  
  if (active) then
     write(*,*)
     write(*,*) 'Computing H*B*H^T using B_static_chm - Start'
  else
     if ( mpi_myid == 0 ) write(*,*) 'compute_HBHT_static_chem: option NOT ACTIVATED'
     return
  end if
          
  call chm_observation_operators(lcolumng,lobsSpaceData,kmode=1) ! kmode=1 for background check to compute HBH^T
  
  write(*,*)
  write(*,*) 'Computing H*B*H^T using B_static_chm - End'
  
  RETURN
END SUBROUTINE hbht_compute_static_chem

!--------------------------------------------------------------------------
!! *Purpose*: Compute background error stddev in observation space using
!!            ensemble-based statistics.
!!
!--------------------------------------------------------------------------
subroutine hbht_compute_ensemble(columng,columnhr,obsSpaceData,active)
  implicit none

  type(struct_obs)        :: obsSpaceData ! Observation-related data
  type(struct_columnData) :: columng      ! Columns of the background interpolated 
                                          ! to analysis levels and to obs horizontal locations
  type(struct_columnData) :: columnhr     ! Columns of the background interpolated 
                                          ! to obs horizontal locations
  logical                 :: active

  type(struct_columnData) :: column
  type(struct_gsv)        :: statevector

  type(struct_hco), pointer :: hco_anl
  type(struct_vco), pointer :: vco_anl

  real(8), allocatable :: HBHT_ens(:)

  integer :: memberIndex, index_body, cvdim

  !
  !- 1.  Initialization
  !

  !- 1.1 Get vertical and horizontal analysis grid attributes
  vco_anl => col_getVco(columng)
  hco_anl => agd_getHco('ComputationalGrid')

  !- 1.2 Initialize/Read the flow-dependent ensemble perturbations
  call ben_Setup( hco_anl,             & ! IN
                  vco_anl,             & ! IN
                  cvdim,               & ! OUT
                  'BackgroundCheck' )    ! IN

  if ( cvdim > 0 ) then
     write(*,*)
     write(*,*) 'Computing HBHT from ensemble perturbations - START'
     active = .true.
  else
     if ( mpi_myid == 0 ) write(*,*) 'compute_HBHT_ensemble: option NOT ACTIVATED'
     active = .false.
     return
  end if

  !- 1.3 Create a gridstatevector to store the ensemble perturbations
  call gsv_allocate(statevector, tim_nstepobsinc, hco_anl, vco_anl, &
       mpi_local_opt=.true.)

  !- 1.4 Create column vectors to store the ens perturbation interpolated to obs horizontal locations
  call col_setVco(column,vco_anl)
  call col_allocate(column,col_getNumCol(columng),mpiLocal_opt=.true.)

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
     call ben_getPerturbation( statevector,    & ! OUT
                               memberIndex,    & ! IN
                               'ConstantValue' ) ! IN

     !- 2.2 Interpolation to the observation horizontal locations
     call s2c_tl( statevector,           & ! IN
                  column,                & ! OUT (H_horiz EnsPert)
                  columng, obsSpaceData )  ! IN

     !- 2.3 Interpolation to observation space
     call oop_Htl( column, columng, & ! IN
                   obsSpaceData,    & ! OUT (Save as OBS_WORK: H_vert H_horiz EnsPert = H EnsPert)
                   1 )                ! IN

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

end subroutine hbht_compute_ensemble


!--------------------------------------------------------------------------
!! *Purpose*: Interpolate vertically the contents of "column" to
!!            the pressure levels of the observations. Then
!!            compute THE FIRST GUESS ERROR VARIANCES
!!            A linear interpolation in ln(p) is performed.
!!
!! @author P. Koclas *CMC/CMSV November 1998
!--------------------------------------------------------------------------
      SUBROUTINE setfgefam(CDFAM,lcolumn,lcolumng,lobsSpaceData)
      IMPLICIT NONE
      type(struct_columnData) :: lcolumn,lcolumng
      type(struct_obs) :: lobsSpaceData
      CHARACTER*2 CDFAM
      INTEGER IPB,IPT
      INTEGER INDEX_HEADER,ITYP,IK
      INTEGER INDEX_BODY
      REAL*8 ZWB,ZWT
      REAL*8 ZLEV,ZPB,ZPT
      character(len=2) :: varLevel

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
            IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,index_body) .EQ. 1 .AND.   &
                 obs_bodyElem_i(lobsSpaceData,OBS_VCO,index_body) .EQ. 2      ) then
            IF  (obs_bodyElem_i(lobsSpaceData,OBS_XTR,index_body) .NE. 0) THEN
               ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body)
               varLevel = vnl_varLevelFromVarnum(ityp)
               IK   = col_getNumLev(LCOLUMNG,varLevel)
               IPB  = IK + col_getOffsetFromVarno(lcolumng,ityp)
               if(ITYP .ne. BUFR_NEGZ) then
                 call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,col_getElem(lcolumn,IPB,INDEX_HEADER))
               else
                 call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,col_getHeight(lcolumn,IK,INDEX_HEADER,'TH'))
               endif
            ELSE
               ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body)
               varLevel = vnl_varLevelFromVarnum(ityp)
               ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,index_body)
               IK   = obs_bodyElem_i(lobsSpaceData,OBS_LYR,index_body)
               IPT  = IK + col_getOffsetFromVarno(lcolumng,ityp)
               IPB  = IPT+1
               ZPT  = col_getPressure(lcolumng,IK,INDEX_HEADER,varLevel)
               ZPB  = col_getPressure(lcolumng,IK+1,INDEX_HEADER,varLevel)
               ZWB  = LOG(ZLEV/ZPT)/LOG(ZPB/ZPT)
               ZWT  = 1.0D0 - ZWB
!
!              FIRST GUESS ERROR VARIANCE
!
               if(ITYP .ne. BUFR_NEGZ) then
                 call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,   &
                      (ZWB*col_getElem(lcolumn,IPB,INDEX_HEADER) + ZWT*col_getElem(lcolumn,IPT,INDEX_HEADER)))
               else
                 call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,   &
                      (ZWB*col_getHeight(lcolumn,IK+1,INDEX_HEADER,'TH') + ZWT*col_getHeight(lcolumn,IK,INDEX_HEADER,'TH')))
               endif
               if(obs_bodyElem_r(lobsSpaceData,OBS_HPHT,index_body).le.0.d0) then
                 write(*,*) 'SETFGEFAM: CDFAM = ',CDFAM
                 write(*,*) 'SETFGEFAM: IPB,IPT,ZWB,ZWT,ITYP,ZLEV=',IPB,IPT,ZWB,ZWT,ITYP,ZLEV
                 write(*,*) 'SETFGEFAM: lcolumn_all(IPB,INDEX_HEADER)=',col_getElem(lcolumn,IPB,INDEX_HEADER)
                 write(*,*) 'SETFGEFAM: lcolumn_all(IPT,INDEX_HEADER)=',col_getElem(lcolumn,IPT,INDEX_HEADER)
                 write(*,*) 'SETFGEFAM: get_height(IK+1,INDEX_HEADER)=',col_getHeight(lcolumn,IK+1,INDEX_HEADER,'TH')
                 write(*,*) 'SETFGEFAM: get_height(IK  ,INDEX_HEADER)=',col_getHeight(lcolumn,IK  ,INDEX_HEADER,'TH')
                 CALL utl_abort('SETFGEFAM: First-guess stdev bad value')
               endif
            ENDIF
            ENDIF

         END DO BODY

      END DO HEADER

      RETURN
      END SUBROUTINE setfgefam


!--------------------------------------------------------------------------
!!    Purpose: Interpolate vertically the contents of "column" to
!!             the levels of the observations (in meters). Then
!!             compute THE FIRST GUESS ERROR VARIANCES
!!             A linear interpolation in z is performed.
!!
!! @author J. St-James, CMDA/SMC November 2002
!!
!--------------------------------------------------------------------------
      SUBROUTINE setfgefamz(CDFAM,lcolumn,lcolumng,lobsSpaceData)
      IMPLICIT NONE
      type(struct_columnData) :: lcolumn,lcolumng
      type(struct_obs) :: lobsSpaceData
      CHARACTER*2 CDFAM
      INTEGER IPB,IPT
      INTEGER INDEX_HEADER,ITYP,IK,IBEGIN,ILAST
      INTEGER J,INDEX_BODY
      integer :: bodyIndexStart, bodyIndexEnd, bodyIndex2
      REAL*8 ZWB,ZWT
      REAL*8 ZLEV,ZPB,ZPT
      real(8) :: azimuth, fge_uu, fge_vv, fge_fam
      character(len=2) :: varLevel

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
!C
!C*    1. Computation of sigmap
!C     .  -----------------------------
!C
               IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,index_body) == 1 .AND. &
                    obs_bodyElem_i(lobsSpaceData,OBS_VCO,index_body) == 1 )then
                  ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body)
                  varLevel = vnl_varLevelFromVarnum(ityp)

                  ! Interpolate the background-covariance statistics
                  IF  (obs_bodyElem_i(lobsSpaceData,OBS_XTR,index_body) /= 0)THEN
                     IK=col_getNumLev(LCOLUMNG,varLevel)-1
                     IPT  = IK + col_getOffsetFromVarno(lcolumng,ityp)
                     IPB  = IPT +1
                     fge_uu = col_getElem(lcolumn,IPB,INDEX_HEADER,'UU')
                     fge_vv = col_getElem(lcolumn,IPB,INDEX_HEADER,'VV')

                  ELSE
                     ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,index_body)
                     IK   = obs_bodyElem_i(lobsSpaceData,OBS_LYR,index_body)
                     IPT  = IK + col_getOffsetFromVarno(lcolumng,ityp)
                     IPB  = IPT+1
                     ZPT  = col_getHeight(lcolumng,IK  ,INDEX_HEADER,varLevel)/RG
                     ZPB  = col_getHeight(lcolumng,IK+1,INDEX_HEADER,varLevel)/RG
                     ZWB  = (ZPT-ZLEV)/(ZPT-ZPB)
                     ZWT  = 1.d0 - ZWB
                     fge_uu =   ZWB*col_getElem(lcolumn,IPB,INDEX_HEADER,'UU') &
                              + ZWT*col_getElem(lcolumn,IPT,INDEX_HEADER,'UU')
                     fge_vv =   ZWB*col_getElem(lcolumn,IPB,INDEX_HEADER,'VV') &
                              + ZWT*col_getElem(lcolumn,IPT,INDEX_HEADER,'VV')
                  ENDIF

                  ! First-Guess Error Variance
                  if(cdfam == 'AL')then
                    ! Scan body indices for the azimuth
                    bodyIndexStart= obs_headElem_i(lobsSpaceData, OBS_RLN, &
                                                   index_header)
                    bodyIndexEnd  = obs_headElem_i(lobsSpaceData, OBS_NLV, &
                                                   index_header)&
                                  + bodyIndexStart - 1
                    BODY_SUPP: do bodyIndex2 = bodyIndexStart, bodyIndexEnd
                      if(obs_bodyElem_i(lobsSpaceData, OBS_VNM, bodyIndex2) &
                         == 5021)then
                        azimuth=obs_bodyElem_r(lobsSpaceData,OBS_VAR,bodyIndex2)&
                                * MPC_RADIANS_PER_DEGREE_R8
                        exit BODY_SUPP
                      end if
                    end do BODY_SUPP

                    fge_fam =sqrt((fge_vv*cos(azimuth))**2 + &
                                  (fge_uu*sin(azimuth))**2)

                  else if(cdfam == 'PR')then
                     fge_fam =   ZWB*col_getElem(lcolumn,IPB,INDEX_HEADER) &
                               + ZWT*col_getElem(lcolumn,IPT,INDEX_HEADER)
                  else
                     write(*,*)"ERROR:  The family, ", cdfam, &
                               ", is not supported by setfgefamz"
                     return
                  end if
                  call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,fge_fam)
               ENDIF

         END DO BODY

      END DO HEADER

      RETURN
      END SUBROUTINE setfgefamz


!--------------------------------------------------------------------------
!! *Purpose*: Interpolate vertically the contents of "column" to
!!            the pressure levels of the observations. Then
!!            compute THE FIRST GUESS ERROR VARIANCES
!!            A linear interpolation in ln(p) is performed.
!!
!! @author P. Koclas *CMC/CMSV November 1998
!!
!--------------------------------------------------------------------------
      SUBROUTINE setfgett(lcolumn,lcolumng,lobsSpaceData)
      IMPLICIT NONE
      type(struct_columnData) :: lcolumn,lcolumng
      type(struct_obs) :: lobsSpaceData
      INTEGER IPB,IPT
      INTEGER INDEX_HEADER,ITYP,IK,IBEGIN,ILAST
      INTEGER J,INDEX_BODY
      REAL*8 ZWB,ZWT
      REAL*8 ZLEV,ZPB,ZPT
      character(len=2) :: varLevel

      ! loop over all body rows
      BODY: do index_body=1,obs_numbody(lobsSpaceData)
!
!*    1. Computation of sigmap
!     .  -----------------------------
!
         IF ( (obs_bodyElem_i(lobsSpaceData,OBS_ASS,index_body) .EQ. 1) .and.    &
              (obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body).EQ. BUFR_NETT) ) THEN

            IF ( (obs_bodyElem_i(lobsSpaceData,OBS_XTR,index_body) .NE. 0) .and.    &
                 (obs_bodyElem_i(lobsSpaceData,OBS_VCO,index_body) .EQ. 2) ) THEN
               ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body)
               varLevel = vnl_varLevelFromVarnum(ityp)
               IK=col_getNumLev(lcolumng,varLevel)-1
               INDEX_HEADER = obs_bodyElem_i(lobsSpaceData,OBS_HIND,index_body)
               IPT  = IK + col_getOffsetFromVarno(lcolumng,ityp)
               IPB  = IPT +1
               call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,col_getElem(lcolumn,IPB,INDEX_HEADER))
            ELSE
               INDEX_HEADER = obs_bodyElem_i(lobsSpaceData,OBS_HIND,index_body)
               ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,index_body)
               IK   = obs_bodyElem_i(lobsSpaceData,OBS_LYR,index_body)
               IPT  = IK
               IPB  = IPT+1
               ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,index_body)
               varLevel = vnl_varLevelFromVarnum(ityp)
               ZPT  = col_getPressure(lcolumng,IK,INDEX_HEADER,varLevel)
               ZPB  = col_getPressure(lcolumng,IK+1,INDEX_HEADER,varLevel)
               ZWB  = LOG(ZLEV/ZPT)/LOG(ZPB/ZPT)
               ZWT  = 1.0D0 - ZWB
!
!              FIRST GUESS ERROR VARIANCE
!
               call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,   &
                  (ZWB*col_getElem(lcolumn,IPB,INDEX_HEADER,'TT') + ZWT*col_getElem(lcolumn,IPT,INDEX_HEADER,'TT')))
            ENDIF

         ENDIF

      END DO BODY

      RETURN
      END SUBROUTINE setfgett


!--------------------------------------------------------------------------
!! *Purpose*: Interpolate vertically the contents of "column" to
!!            the pressure levels of the observations.
!!            A linear interpolation in ln(p) is performed.
!!
!! @author P. Koclas *CMC/AES  September 2000
!!
!--------------------------------------------------------------------------
  subroutine setfgeSurf( lcolumn, lcolumng, lobsSpaceData )
  
  implicit none
  ! arguments
  type(struct_columnData) :: lcolumn, lcolumng
  type(struct_obs)        :: lobsSpaceData
  ! locals
  integer          :: ipb, ipt, idim, headerIndex, ik, bodyIndex, ityp
  real(8)          :: zwb, zwt, zlev, zpt, zpb, zhhh
  character(len=2) :: cfam, varLevel
  logical          :: llok

  ! loop over all body rows
  BODY: do bodyIndex = 1, obs_numbody( lobsSpaceData )

    cfam = obs_getFamily( lobsSpaceData, bodyIndex = bodyIndex )
    if( cfam == 'SF'.or. cfam == 'TM' .or. cfam == 'UA' .or. cfam  == 'SC' .or. cfam == 'GP' ) then

      ! Process all data within the domain of the model (excluding GB-GPS ZTD data)
      llok = .false.

      if ( obs_bodyElem_i( lobsSpaceData, OBS_VCO, bodyIndex ) == 1 ) then

        ityp = obs_bodyElem_i( lobsSpaceData, OBS_VNM, bodyIndex )
        if ( ityp == BUFR_NETS .or. ityp == BUFR_NEPS .or. ityp == BUFR_NEPN .or. ityp == BUFR_NESS .or. &
             ityp == BUFR_NEUS .or. ityp == BUFR_NEVS .or. ityp == BUFR_NEFS .or. ityp == BUFR_NEDS .or. &
             ityp == bufr_sst ) then

          llok = ( obs_bodyElem_i( lobsSpaceData, OBS_ASS, bodyIndex ) == 1 )

        else if ( ityp == BUFR_NEZD ) then

          ! make sure total zenith delay (from ground-based GPS) not treated
          llok=.false.

        else

          llok=(obs_bodyElem_i( lobsSpaceData, OBS_ASS, bodyIndex ) == 1 .and. &
                obs_bodyElem_i( lobsSpaceData, OBS_XTR, bodyIndex ) >= 0)
          if ( llok ) write(*,*) 'setfgesurf: WARNING!!! unknown obs seen'
          if ( llok ) write(*,*) 'setfgesurf: ityp=',ityp,', cfam=',cfam

        end if

        if ( llok ) then

          headerIndex = obs_bodyElem_i( lobsSpaceData, OBS_HIND, bodyIndex )
          ityp         = obs_bodyElem_i( lobsSpaceData, OBS_VNM , bodyIndex )
          varLevel     = vnl_varLevelFromVarnum( ityp )
          idim = 1
          if ( varLevel == 'SF') idim = 0
          ik   = obs_bodyElem_i( lobsSpaceData, OBS_LYR, bodyIndex )
          zlev = obs_bodyElem_r( lobsSpaceData, OBS_PPP, bodyIndex )
          zhhh = zlev * grav

          if ( ityp == BUFR_NETS .or. ityp == BUFR_NEPS .or. ityp == BUFR_NEPN .or. &
               ityp == BUFR_NESS .or. ityp == BUFR_NEUS .or. ityp == BUFR_NEVS ) then

            ipt  = ik + col_getOffsetFromVarno(lcolumng,ityp)
            ipb  = ipt+1
            call obs_bodySet_r( lobsSpaceData, OBS_HPHT, bodyIndex, col_getElem( lcolumn, ipb, headerIndex ) )

          else

            ipt  = ik + col_getOffsetFromVarno(lcolumng,ityp)
            ipb  = ipt+1
            zpt  = col_getHeight(lcolumng,ik,headerIndex,varLevel)
            zpb  = col_getHeight(lcolumng,ik+1,headerIndex,varLevel)
            zwb  = idim*(zpt-zhhh)/(zpt-zpb)
            zwt  = 1.d0 - zwb
    
            if ( obs_bodyElem_i( lobsSpaceData, OBS_XTR, bodyIndex ) == 0 ) then

              call obs_bodySet_r( lobsSpaceData, OBS_HPHT, bodyIndex,   &
              zwb * col_getElem( lcolumn, ipb, headerIndex ) + zwt * col_getElem( lcolumn, ipt, headerIndex ))

            else

              call obs_bodySet_r( lobsSpaceData, OBS_HPHT, bodyIndex,   &
                col_getElem( lcolumn, ik + col_getOffsetFromVarno( lcolumng, ityp ), headerIndex ))

            end if

            if(obs_elem_c( lobsSpaceData, 'STID', headerIndex ) == '99999999' ) then

              write(*,*) 'setfgesurf: stn, ityp, xtr, ipt, ipb, zwt, zwb',  &
                   obs_elem_c( lobsSpaceData, 'STID', headerIndex ), ityp, &
                   obs_bodyElem_i( lobsSpaceData, OBS_XTR, bodyIndex ), ipt, ipb, zwt, zwb
              write(*,*) 'setfgesurf: gobs(ipb), gobs(ipt), fge',   &
                    col_getElem( lcolumn, ipb, headerIndex ), col_getElem( lcolumn, ipt, headerIndex ), &
                    obs_bodyElem_r( lobsSpaceData, OBS_HPHT, bodyIndex )

            endif

          end if
        end if
      end if

    end if

  end do BODY

  end subroutine setfgeSurf


!--------------------------------------------------------------------------
!! *Purpose*: Construct the FIRST GUESS ERROR VARIANCES from the
!!            diff-calculated dependencies and the primary errors.
!!
!! @author J.M. Aparicio *MSC/ARMA Nov 2004
!!
!!          Adapted Nov 2012 for both refractivity and bending angle data
!!
!--------------------------------------------------------------------------
      SUBROUTINE setfgedif(CDFAM,lcolumng,lobsSpaceData)
      IMPLICIT NONE
!C
      type(struct_columnData) :: lcolumng
      type(struct_obs)        :: lobsSpaceData
!C
      INTEGER INDEX_HEADER, IDATYP, INDEX_BODY, iProfile
      CHARACTER*2 CDFAM
      REAL*8 zLat, Lat, sLat
      REAL*8 zLon, Lon
      REAL*8 zAzm, Azm
      INTEGER IAZM, ISAT
      REAL*8 Rad, Geo, WFGPS
      REAL*8, allocatable :: zPP(:)
      REAL*8, allocatable :: zDP(:)
      REAL*8, allocatable :: zTT(:)
      REAL*8, allocatable :: zHU(:)
      REAL*8, allocatable :: zUU(:)
      REAL*8, allocatable :: zVV(:)
      INTEGER JF, stat
      INTEGER JL, JJ
      REAL*8 ZP0, ZMT
      REAL*8 HNH1, ZFGE, ZERR
      INTEGER JV, NGPSLEV, NWNDLEV
      LOGICAL  ASSIM, LFIRST, FIRSTHEADER

      INTEGER NH, NH1

!      REAL*8 JAC(ngpscvmx)
      REAL*8 DV (ngpscvmx)
      TYPE(GPS_PROFILE)           :: PRF
      REAL*8       , allocatable :: H   (:),AZMV(:)
      TYPE(GPS_DIFF), allocatable :: RSTV(:),RSTVP(:),RSTVM(:)
      type(struct_vco), pointer  :: vco_anl

      WRITE(*,*)'ENTER SETFGEDIFF'
!C
!C     * 1.  Initializations
!C     *     ---------------
!C
      NGPSLEV=col_getNumLev(lcolumng,'TH')
      NWNDLEV=col_getNumLev(lcolumng,'MM')
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
!C         IF (LEVELGPSRO.EQ.1) THEN
!C            allocate( RSTVP(GPSRO_MAXPRFSIZE) )
!C            allocate( RSTVM(GPSRO_MAXPRFSIZE) )
!C         ENDIF
      endif

      vco_anl => col_getVco(lcolumng)
!C
!C    Loop over all header indices of the 'RO' family:
!C
      call obs_set_current_header_list(lobsSpaceData,CDFAM)
      FIRSTHEADER=.TRUE.
      HEADER: do
         index_header = obs_getHeaderIndex(lobsSpaceData)
         if (index_header < 0) exit HEADER
!C
!C     * Process only refractivity data (codtyp 169)
!C
         IDATYP = obs_headElem_i(lobsSpaceData,OBS_ITY,INDEX_HEADER)
         IF ( IDATYP .EQ. 169 ) THEN
!C
!C     *    Scan for requested data values of the profile, and count them
!C
            ASSIM = .FALSE.
            NH = 0
            call obs_set_current_body_list(lobsSpaceData, INDEX_HEADER)
            BODY: do
               INDEX_BODY = obs_getBodyIndex(lobsSpaceData)
               if (INDEX_BODY < 0) exit BODY
               IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY).EQ.1 ) THEN
                  ASSIM = .TRUE.
                  NH = NH + 1
               ENDIF
            ENDDO BODY
!C
!C     *    If assimilations are requested, prepare and apply the observation operator
!C
            IF (ASSIM) THEN
               iProfile=gps_iprofile_from_index(INDEX_HEADER)
!C
!C     *       Profile at the observation location:
!C
               if (.not.gps_vRO_lJac(iProfile)) then
!C
!C     *          Basic geometric variables of the profile:
!C
                  zLat = obs_headElem_r(lobsSpaceData,OBS_LAT,INDEX_HEADER)
                  zLon = obs_headElem_r(lobsSpaceData,OBS_LON,INDEX_HEADER)
                  IAZM = obs_headElem_i(lobsSpaceData,OBS_AZA,INDEX_HEADER)
                  ISAT = obs_headElem_i(lobsSpaceData,OBS_SAT,INDEX_HEADER)
                  Rad  = obs_headElem_r(lobsSpaceData,OBS_TRAD,INDEX_HEADER)
                  Geo  = obs_headElem_r(lobsSpaceData,OBS_GEOI,INDEX_HEADER)
                  zAzm = 0.01d0*IAZM / MPC_DEGREES_PER_RADIAN_R8
                  zMT  = col_getHeight(lcolumng,NGPSLEV,INDEX_HEADER,'TH')/RG
                  WFGPS= 0.d0
                  DO JJ=1,NUMGPSSATS
                     IF (ISAT.EQ.IGPSSAT(JJ)) WFGPS=WGPS(JJ)
                  ENDDO
                  Lat  = zLat * MPC_DEGREES_PER_RADIAN_R8
                  Lon  = zLon * MPC_DEGREES_PER_RADIAN_R8
                  Azm  = zAzm * MPC_DEGREES_PER_RADIAN_R8
                  sLat = sin(zLat)
                  zMT  = zMT * RG / gpsgravitysrf(sLat)
                  zP0  = col_getElem(lcolumng,1,INDEX_HEADER,'P0')
                  DO JL = 1, NGPSLEV
!C
!C     *             Profile x
!C
                     zPP(JL) = col_getPressure(lcolumng,JL,INDEX_HEADER,'TH')
!C     *             True implementation of zDP (dP/dP0)
                     zDP(JL) = col_getPressureDeriv(lcolumng,JL,INDEX_HEADER,'TH')
                     zTT(JL) = col_getElem(lcolumng,JL,INDEX_HEADER,'TT') - p_TC
                     zHU(JL) = col_getElem(lcolumng,JL,INDEX_HEADER,'HU')
                     zUU(JL) = 0.d0
                     zVV(JL) = 0.d0
                  ENDDO
                  DO JL = 1, NWNDLEV
                     zUU(JL) = col_getElem(lcolumng,JL,INDEX_HEADER,'UU') * p_knot
                     zVV(JL) = col_getElem(lcolumng,JL,INDEX_HEADER,'VV') * p_knot
                  ENDDO
                  zUU(NGPSLEV) = zUU(NWNDLEV)
                  zVV(NGPSLEV) = zUU(NWNDLEV)
!C     
!C     *          GPS profile structure:
!C
                  call gps_struct1sw(ngpslev,zLat,zLon,zAzm,zMT,Rad,geo,zP0,zPP,zDP,zTT,zHU,zUU,zVV,prf)
!C
!C     *          Prepare the vector of all the observations:
!C
                  NH1 = 0
                  call obs_set_current_body_list(lobsSpaceData, INDEX_HEADER)
                  BODY_2: do
                     INDEX_BODY = obs_getBodyIndex(lobsSpaceData)
                     if (INDEX_BODY < 0) exit BODY_2
                     IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY).EQ.1 ) THEN
                        NH1      = NH1 + 1
                        H(NH1)   = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)
                        AZMV(NH1)= zAzm
                     ENDIF
                  ENDDO BODY_2
!C
!C     *          Apply the observation operator:
!C
                  IF (LEVELGPSRO.EQ.1) THEN
                     CALL GPS_BNDOPV1(H      , AZMV, NH, PRF, RSTV)
!C                     CALL GPS_BNDOPV1(H+WFGPS, AZMV, NH, PRF, RSTVP)
!C                     CALL GPS_BNDOPV1(H-WFGPS, AZMV, NH, PRF, RSTVM)
!C                     do nh1 = 1, nh
!C                        RSTV(nh1)=(RSTVP(nh1)+RSTV(nh1)+RSTVM(nh1))/3.d0
!C                     enddo
                  ELSE
                     CALL GPS_REFOPV (H,       NH, PRF, RSTV)
                  ENDIF
                  DO NH1=1,NH
                     gps_vRO_Jacobian(iProfile,NH1,:)= RSTV(NH1)%DVAR(1:2*NGPSLEV+1)
                  ENDDO
                  gps_vRO_lJac(iProfile)=.true.
               endif
!C
!C     *       Local error
!C
               DO JL = 1, NGPSLEV
                  DV (        JL) = 1.d0
                  DV (NGPSLEV+JL) = 1.d0
               ENDDO
               DV (2*NGPSLEV+1)   = 2.d0
!C
!C     *       Perform the H(xb)DV operation:
!C
               NH1 = 0
               call obs_set_current_body_list(lobsSpaceData, INDEX_HEADER)
               BODY_3: do
                  INDEX_BODY = obs_getBodyIndex(lobsSpaceData)
                  if (INDEX_BODY < 0) exit BODY_3
                  IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY).EQ.1 ) THEN
                     NH1 = NH1 + 1
!C
!C     *             Observation jacobian
!C
!                     JAC = RSTV(NH1)%DVAR
!C
!C     *             Evaluate sqrt( H(xb)DV **2 )
!C
                     ZFGE = 0.d0
                     DO JV = 1, 2*PRF%NGPSLEV+1
                        ZFGE = ZFGE + (gps_vRO_Jacobian(iProfile,NH1,JV) * DV(JV))**2
                     ENDDO
                     ZFGE = SQRT(ZFGE)
                     ZERR = obs_bodyElem_r(lobsSpaceData,OBS_OER,INDEX_BODY)
!C     
!C     *             FIRST GUESS ERROR VARIANCE
!C
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
!C         IF (LEVELGPSRO.EQ.1) THEN
!C            deallocate( RSTVM )
!C            deallocate( RSTVP )
!C         ENDIF
         deallocate( RSTV )
         deallocate( AZMV )
         deallocate( H    )

         deallocate(zVV)
         deallocate(zUU)
         deallocate(zHU)
         deallocate(zTT)
         deallocate(zDP)
         deallocate(zPP)
      ENDIF

      WRITE(*,*)'EXIT SETFGEDIFF'
      RETURN
      END SUBROUTINE setfgedif


!--------------------------------------------------------------------------
!! *Purpose*: Set FGE for all GPS ZTD observations using
!!            Jacobians from ZTD observation operator
!!
!! OPTION: Test ZTD operators (compares H(x+dx)-H(x) with (dH/dx)*dx
!!         when LTESTOP = .true.)
!!
!! @author S. Macpherson *ARMA/MSC  December 2004
!!
!! Revisions:
!!
!!        -S. Macpherson *ARMA/MSC  18 March 2010
!!           - add optional NL, TL and AD operator tests
!!        -S. Macpherson *ARMA/MSC   August 2010
!!           - use new GPS ZTD observation operator (from GPS-RO modules)
!!        -S. Macpherson *ARMA/MSC   December 2012
!!           - update from Rev189 to Rev213
!!           - use new ZTD-specific GPS modules modgps04profilezd, modgps08ztdop
!!           - LTESTOP option now set in 3dvar namelist
!!           - if numGPSZTD=0, does nothing and returns
!!        -S. Macpherson *ARMA/MSC   November 2014
!!           - add surface pressure (P0) argument to call gps_structztd()
!!
!!v     *********************************************************************
!!v     ****                   9 October 2015                            ****
!!v     ****                                                             ****
!!v     **** NOTE: Effective Rev644M, this routine is no longer used!    ****
!!v     ****       FGE for ZTD is no longer needed for background check. ****
!!v     ****       Routine is only called when LTESTOP=.true., in which  ****
!!v     ****       case the operator test only is done.                  ****
!!v     ****                                                             ****
!!v     *********************************************************************
!!
!--------------------------------------------------------------------------
      SUBROUTINE setfgegps(lcolumn,lcolumng,lobsSpaceData)
      IMPLICIT NONE
!   lcolumn  contains background errors for control variables on model levels
!   lcolumng contains lo-res first guess profiles at obs locations
      type(struct_columnData) :: lcolumn, lcolumng
      type(struct_obs) :: lobsSpaceData
      type(struct_vco), pointer :: vco_anl
      REAL*8 ZLAT, Lat
      REAL*8 ZLON, Lon
      REAL*8, allocatable :: ZPP(:)
      REAL*8, allocatable :: ZDP(:)
      REAL*8, allocatable :: ZTT(:)
      REAL*8, allocatable :: ZHU(:)
      REAL*8, allocatable :: ZGZ(:)
      REAL*8, allocatable :: ZTTB(:)
      REAL*8, allocatable :: ZHUB(:)
      REAL*8, allocatable :: ZQQB(:)
      REAL*8, allocatable :: ZQQ(:)
      REAL*8, allocatable :: ZTTB_P(:)
      REAL*8, allocatable :: ZQQB_P(:)
      REAL*8, allocatable :: RZHUB_P(:)
      REAL*8, allocatable :: ZPP_P(:)
      
      REAL*8 ZP0
      REAL*8 ZP0B, ZP0B_P
      REAL*8 ZMT, ZTOP, ZBOT
 
      REAL*8 JAC(ngpscvmx)
      REAL*8 DX (ngpscvmx)

      REAL*8 ZOER, ZLEV, ZTDOBS, ZVAR, ZPSMOD
      REAL*8 ZJP0, ZLSUM
      REAL*8 DELTAH_NL, DELTAH_TL
      REAL*8 PERTFAC, ZTDM
      REAL*8 ZDZMIN, ZSUMTEST

      INTEGER INDEX_HEADER, FIRST_HEADER
      INTEGER IDATYP, ITYP
      INTEGER IDATA, IDATEND, INDEX_BODY
      INTEGER JL, JK, NFLEV_T, ILYR, IOBS
      INTEGER INOBS_OPT, INOBS_JAC, icount, status, iversion

      LOGICAL  ASSIM, LLOK, LSTAG
      CHARACTER*9  STN_JAC
      
      CHARACTER(len=2) :: varLevel
      
      TYPE(GPS_PROFILEZD)    :: PRF, PRFP
      TYPE(GPS_DIFF)         :: ZTDopv, ZTDopvP

      IF (numGPSZTD .EQ. 0) RETURN

!C
!C     * 1.  Initializations
!C     *     ---------------
!C
      NFLEV_T = col_getNumLev(lcolumng,'TH')
      allocate(ZPP(NFLEV_T))
      allocate(ZDP(NFLEV_T))
      allocate(ZTT(NFLEV_T))
      allocate(ZHU(NFLEV_T))
      allocate(ZGZ(NFLEV_T))
      allocate(ZTTB(NFLEV_T))
      allocate(ZHUB(NFLEV_T))
      allocate(ZQQB(NFLEV_T))
      allocate(ZQQ(NFLEV_T))

!c     Number of locations/sites for observation operator test
      INOBS_OPT = 50
!c     Number of locations/sites for Jacobian printout
      INOBS_JAC  = 5
!c     Factor to multiply background errors for perturbation vector
      PERTFAC = 0.75d0
!C
      STN_JAC = 'FSL_BRFT '
!c      
      ZDZMIN = DZMIN
!c
      vco_anl => col_getVco(lcolumng)
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

!C
      IF ( .NOT.LTESTOP ) THEN

      first_header=-1
      icount = 0
!C
      ! loop over all header indices of the 'GP' family
      call obs_set_current_header_list(lobsSpaceData,'GP')
      HEADER: do
         index_header = obs_getHeaderIndex(lobsSpaceData)
         if (index_header < 0) exit HEADER
         if (first_header .eq. -1) first_header = index_header
!C     
!C     *    .   Process only zenith delay data (codtyp 189 and BUFR_NEZD)
!C
               IDATYP = obs_headElem_i(lobsSpaceData,OBS_ITY,INDEX_HEADER)
               IF ( IDATYP .EQ. 189 ) THEN
!C
!C                 Loop over data in the observations
!C
                  IDATA   = obs_headElem_i(lobsSpaceData,OBS_RLN,INDEX_HEADER)
                  IDATEND = obs_headElem_i(lobsSpaceData,OBS_NLV,INDEX_HEADER) + IDATA - 1
                  ASSIM = .FALSE.
!C
!C                 Scan for requested assimilations, and count them.
!C
                  DO INDEX_BODY= IDATA, IDATEND
                     ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,INDEX_BODY)
                     LLOK = ( (ITYP .EQ. BUFR_NEZD) .AND. (obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY) .EQ. 1) )
                     IF ( LLOK ) THEN
                        ASSIM = .TRUE.
                        ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)
                        icount = icount + 1
                     ENDIF
                  ENDDO
!C
!C     *           If assimilations are requested, apply the AD observation operator
!C
                  IF (ASSIM) THEN
!C     
!C     *        LR background profile and background errors at the observation location x :
!C
                     Lat  = obs_headElem_r(lobsSpaceData,OBS_LAT,INDEX_HEADER)
                     Lon  = obs_headElem_r(lobsSpaceData,OBS_LON,INDEX_HEADER)
                     ZLAT = Lat * MPC_DEGREES_PER_RADIAN_R8
                     ZLON = Lon * MPC_DEGREES_PER_RADIAN_R8
                     ZP0B = col_getElem(lcolumng,1,INDEX_HEADER,'P0')
                     DO JL = 1, NFLEV_T
                       ZPP(JL)  = col_getPressure(lcolumng,JL,INDEX_HEADER,'TH')
!C                     Get ZDP = dP/dP0
                       ZDP(JL)  = col_getPressureDeriv(lcolumng,JL,INDEX_HEADER,'TH')
                       ZTTB(JL) = col_getElem(lcolumng,JL,INDEX_HEADER,'TT')- 273.15d0
                       ZTT(JL)  = col_getElem(lcolumn,JL,INDEX_HEADER,'TT')
                       DX(JL)   = ZTT(JL)
                       ZHUB(JL) = col_getElem(lcolumng,JL,INDEX_HEADER,'HU')
                       ZQQB(JL) = ZHUB(JL)
                       ZHU(JL)  = col_getElem(lcolumn,JL,INDEX_HEADER,'HU')
                       DX(NFLEV_T+JL) = ZHU(JL)
                       ZGZ(JL)  = col_getHeight(lcolumng,JL,INDEX_HEADER,'TH')
                     ENDDO
                     ZP0  = col_getElem(lcolumn,1,INDEX_HEADER,'P0')
                     DX(2*NFLEV_T+1) = ZP0
                     ZMT  = ZGZ(NFLEV_T)/GRAV
                     CALL gps_structztd(NFLEV_T,Lat,Lon,ZMT,ZP0B,ZPP,ZDP,ZTTB,ZHUB,LBEVIS,IREFOPT,PRF)
                     CALL gps_ztdopv(ZLEV,PRF,LBEVIS,ZDZMIN,ZTDopv,ZPSMOD,IZTDOP)
                     JAC = ZTDopv%DVar
!c
                     DO INDEX_BODY= IDATA, IDATEND
                        ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,INDEX_BODY)
                        IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY).EQ.1 .AND. ITYP.EQ.BUFR_NEZD ) THEN
!C
!C     *                    Observation error    SDERR
!c                           ZOER = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)

!C     *                    Observation height (m)
                           ZLEV = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)

                           ZLSUM  = 0.0d0
!C
                           DO JL = 1, 2*NFLEV_T+1
                             ZLSUM = ZLSUM + (JAC(JL)*DX(JL))**2
                           ENDDO
                           call obs_bodySet_r(lobsSpaceData,OBS_HPHT,index_body,SQRT(ZLSUM))

                           IF (icount .LE. INOBS_JAC) THEN
!                           IF ( obs_elem_c(lobsSpaceData,'STID',INDEX_HEADER) .EQ. STN_JAC ) THEN
                             WRITE(*,'(A11,A9)') 'SETFGEGPS: ',obs_elem_c(lobsSpaceData,'STID',INDEX_HEADER)
                             WRITE(*,*) '  ZTD, ZTD FGE = ', ZTDopv%Var, SQRT(ZLSUM)
                             WRITE(*,'(A11,A9,3(1x,f7.2))')   &
                               'SETFGEGPS: ',obs_elem_c(lobsSpaceData,'STID',INDEX_HEADER),ZLAT,ZLON,ZLEV
                             WRITE(*,*) 'JL JACT JACQ FGE_T FGE_LQ QQ'
                             DO JL = 1, NFLEV_T
                               WRITE(*,'(1X,I2,5(1x,E13.6))') JL,JAC(JL),JAC(JL+NFLEV_T)/ZQQB(JL),ZTT(JL),ZHU(JL),ZQQB(JL)
                             ENDDO                         
                             WRITE(*,*) 'JACPS FGE_PS'
                             WRITE(*,'(2(1x,E13.6))') JAC(2*NFLEV_T+1), ZP0
                           ENDIF

                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF

      ENDDO HEADER
      
      ENDIF

!c---------------------------------------------------------------------------------------------------------------

      IF ( LTESTOP ) THEN
      
      allocate(ZTTB_P(NFLEV_T))
      allocate(ZQQB_P(NFLEV_T))
      allocate(ZPP_P(NFLEV_T))

      icount = 0
      ZSUMTEST = 0
!C
      ! loop over all header indices of the 'GP' family
      call obs_set_current_header_list(lobsSpaceData,'GP')
      HEADER2: DO
         index_header = obs_getHeaderIndex(lobsSpaceData)
         if (index_header < 0) exit HEADER2
         if (icount > INOBS_OPT ) exit HEADER2

!C       Loop over data in the observations
!C
         IDATA   = obs_headElem_i(lobsSpaceData,OBS_RLN,INDEX_HEADER)
         IDATEND = obs_headElem_i(lobsSpaceData,OBS_NLV,INDEX_HEADER) + IDATA - 1
!C     
!C       LR background profile and background errors at the observation location x :
!C
         Lat  = obs_headElem_r(lobsSpaceData,OBS_LAT,INDEX_HEADER)
         Lon  = obs_headElem_r(lobsSpaceData,OBS_LON,INDEX_HEADER)
         ZLAT = Lat * MPC_DEGREES_PER_RADIAN_R8
         ZLON = Lon * MPC_DEGREES_PER_RADIAN_R8
         ZP0B = col_getElem(lcolumng,1,INDEX_HEADER,'P0')
         DO JL = 1, NFLEV_T
            ZPP(JL)  = col_getPressure(lcolumng,JL,INDEX_HEADER,'TH')
!C          Get ZDP = dP/dP0
            ZDP(JL)  = col_getPressureDeriv(lcolumng,JL,INDEX_HEADER,'TH')
            ZTTB(JL) = col_getElem(lcolumng,JL,INDEX_HEADER,'TT')- 273.15d0
            ZTT(JL)  = col_getElem(lcolumn,JL,INDEX_HEADER,'TT') * PERTFAC
            ZQQB(JL) = col_getElem(lcolumng,JL,INDEX_HEADER,'HU')
            ZQQ(JL)  = col_getElem(lcolumn,JL,INDEX_HEADER,'HU') * PERTFAC
            ZGZ(JL)  = col_getHeight(lcolumng,JL,INDEX_HEADER,'TH')
         ENDDO
         ZP0  = col_getElem(lcolumn,1,INDEX_HEADER,'P0') * PERTFAC
         ZMT  = ZGZ(NFLEV_T)/GRAV

         DO JL = 1, NFLEV_T
             DX (      JL) = ZTT(JL)
             DX (NFLEV_T+JL) = ZQQ(JL)
         ENDDO
         DX (2*NFLEV_T+1) = ZP0

         ZTDOBS = -1.0d0
         DO INDEX_BODY = IDATA, IDATEND
           ITYP = obs_bodyElem_i(lobsSpaceData,OBS_VNM,INDEX_BODY)
           IOBS = obs_bodyElem_i(lobsSpaceData,OBS_HIND,INDEX_BODY)
           IF ( obs_bodyElem_i(lobsSpaceData,OBS_ASS,INDEX_BODY).EQ.1 .AND. ITYP .EQ. BUFR_NEZD ) THEN
             varLevel = vnl_varLevelFromVarnum(ITYP)
             ZTDOBS  = obs_bodyElem_r(lobsSpaceData,OBS_VAR,INDEX_BODY)
             ZLEV    = obs_bodyElem_r(lobsSpaceData,OBS_PPP,INDEX_BODY)
             ILYR    = obs_bodyElem_i(lobsSpaceData,OBS_LYR,INDEX_BODY)
             ZTOP    = col_getHeight(lcolumng,ILYR,IOBS,varLevel)/GRAV
             if ( ILYR .LT. NFLEV_T ) then
               ZBOT    = col_getHeight(lcolumng,ILYR+1,IOBS,varLevel)/GRAV
             else
               ZBOT    = ZTOP
             endif
             icount  = icount + 1
           ENDIF
         ENDDO

         IF ( ZTDOBS .GT. 0.d0 ) THEN
!c         Create the pertubation control vector
           DO JL = 1, NFLEV_T
             ZPP_P(JL)  = ZPP(JL)  + ZDP(JL)*ZP0
             ZTTB_P(JL) = ZTTB(JL) + ZTT(JL)
             ZQQB_P(JL) = ZQQB(JL) + ZQQ(JL)
           ENDDO
           ZP0B_P = ZP0B + ZP0
!C
!C         Non-linear observation operator --> delta_H = H(x+delta_x) - H(x)
!c
           CALL gps_structztd(NFLEV_T,Lat,Lon,ZMT,ZP0B,ZPP,ZDP,ZTTB,ZQQB,LBEVIS,IREFOPT,PRF)
           CALL gps_structztd(NFLEV_T,Lat,Lon,ZMT,ZP0B_P,ZPP_P,ZDP,ZTTB_P,ZQQB_P,LBEVIS,IREFOPT,PRFP)
           CALL gps_ztdopv(ZLEV,PRF,LBEVIS,ZDZMIN,ZTDopv,ZPSMOD,IZTDOP)
           JAC  = ZTDopv%DVar
           ZTDM = ZTDopv%Var
           CALL gps_ztdopv(ZLEV,PRFP,LBEVIS,ZDZMIN,ZTDopvP,ZPSMOD,IZTDOP)
           DELTAH_NL = ZTDopvP%Var - ZTDopv%Var
!c
!c         Linear  --> delta_H = dH/dx * delta_x
!c
           DELTAH_TL = 0.0d0
           DO JL = 1, 2*NFLEV_T+1
             DELTAH_TL = DELTAH_TL + JAC(JL)*DX(JL)
           ENDDO
!c
           WRITE(*,*) 'SETFGEGPS: GPS ZTD OBSOP TEST FOR SITE ', obs_elem_c(lobsSpaceData,'STID',INDEX_HEADER)
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
      deallocate(ZPP_P)

      ENDIF
!-----------------------------------------------------------------------------------------------------------

      deallocate(ZPP)
      deallocate(ZDP)
      deallocate(ZTT)
      deallocate(ZHU)
      deallocate(ZGZ)
      deallocate(ZTTB)
      deallocate(ZHUB)
      deallocate(ZQQB)
      deallocate(ZQQ)

      RETURN
      END SUBROUTINE setfgegps

end module computeHBHT_mod
