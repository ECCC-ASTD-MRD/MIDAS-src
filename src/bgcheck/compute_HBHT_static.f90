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
!! *Purpose*: Compute background error stddev in observation space using
!!            fixed statistics specific in stats file.
!!
!--------------------------------------------------------------------------
SUBROUTINE COMPUTE_HBHT_STATIC(lcolumng,lcolumnhr,lobsSpaceData,active)
      use mpivar_mod
      use EarthConstants_mod
      use MathPhysConstants_mod
      use obsSpaceData_mod
      use columnData_mod
      use gridStateVector_mod
      use verticalCoord_mod
      use horizontalCoord_mod
      use analysisGrid_mod
      use bmatrixhi_mod
      use obsOperators_mod
      use gps_mod
      use utilities_mod
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
      call col_copyLatLon(lcolumng,lcolumn)
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

      call bilin(lcolumn,statevector,lobsSpaceData)

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

      call bilin(lcolumn,statevector,lobsSpaceData)
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

      call bilin(lcolumn,statevector,lobsSpaceData)
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

      call bilin(lcolumn,statevector,lobsSpaceData)
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
      
    END SUBROUTINE COMPUTE_HBHT_STATIC