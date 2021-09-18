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

module analysisGrid_mod
  ! MODULE analysisGrid_mod (prefix='agd' category='7. Low-level data objects and utilities')
  !
  ! :Purpose: Performs horizontal grid-point variable transforms 
  !           for the limited-area computational analysis grids (extended and
  !           non-extended).
  !
  use earthConstants_mod
  use MathPhysConstants_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use mpi_mod
  use mpivar_mod
  use utilities_mod
  implicit none
  save
  private

  ! public procedures
  public :: agd_SetupFromHCO, agd_mach, agd_mach_r4
  public :: agd_PsiChiToUV, agd_PsiChiToUVAdj, agd_UVToVortDiv
  public :: agd_createLamTemplateGrids

  ! Definition of some parameters characterizing the geometry of
  ! the Limited-Area (LA) analysis grid and associated metric factors
  type :: struct_glmf
     real(8), allocatable :: rlat   (:) ! Latitudes of Scalar gridpoints
     real(8), allocatable :: rlon   (:) ! Longitudes of Scalar gridpoints
     real(8)              :: rdlon      ! Latitude difference of gridpoints
     real(8)              :: rdlat      ! Longitude differences of gridpoints
     real(8), allocatable :: cos2 (:)   ! Grid metric for Psi-Chi to U-V
     real(8), allocatable :: cos2h(:)   ! Grid metric for Psi-Chi to U-V
     real(8), allocatable :: cos2vd (:) ! Grid metric for U-V to Vort-Div
     real(8), allocatable :: cos2hvd(:) ! Grid metric for U-V to Vort-Div
     real(8), allocatable :: idmuh(:)   ! Grid metric for U-V to Vort-Div
     real(8), allocatable :: idmu (:)   ! Grid metric for U-V to Vort-Div
     real(8)              :: dx         ! Grid-point spacing (uniform)
     real(8), allocatable :: conphy (:) ! to go from Wind Images to true winds
     real(8), allocatable :: conima (:) ! to go from true winds to Wind Images
  end type struct_glmf

  type(struct_hco), pointer :: hco_core => null()
  type(struct_hco), pointer :: hco_ext => null()
  type(struct_glmf):: glmf

  integer :: ni_ext, nj_ext   ! With    Extension for Bi-Fourrier
  integer :: ni_core, nj_core ! Without Extension for Bi-Fourrier
  integer :: ext_i , ext_j    ! Extension gridpoints

  integer :: ni_ext_per, nj_ext_per

  integer :: istart, iend, jstart, jend

  integer :: LonPerPE, LonPerPEmax, myLonBeg, myLonEnd
  integer :: LatPerPE, LatPerPEmax, myLatBeg, myLatEnd

  logical :: initialized = .false.

contains

  !--------------------------------------------------------------------------
  ! agd_SetupFromHCO
  !--------------------------------------------------------------------------
  subroutine agd_SetupFromHCO( hco_ext_in, hco_core_opt )
    implicit none

    type(struct_hco), pointer :: hco_ext_in
    type(struct_hco), pointer, optional :: hco_core_opt

    real(8), allocatable :: rlath  (:) ! Latitudes of half grid-points of gridpoints in lat-direction
    real(8), allocatable :: rrcos  (:) ! 1.0/Cos(Latitudes of gridpoints)
    !real(8), allocatable :: rrcosh (:) ! 1.0/Cos(Half-Latitudes of gridpoints)
    real(8), allocatable :: rdmu   (:) ! Differences of mu=sin(lat)
    real(8), allocatable :: rdmuh  (:) ! Differences of muh=sin(lath)
    real(8), allocatable :: r1mmu2 (:) ! (1.-mu**2)
    real(8), allocatable :: r1mmu2h(:) ! (1.-muh**2)

    real(8) :: dlon_test, dlat_test, dlon_ref, dlat_ref

    integer :: i, j

    logical, save :: firstCall = .true.

    ! Ensure subroutine only runs one time during program execution
    if (firstCall) then
      firstCall = .false.
    else
      return
    end if

    write(*,*)
    write(*,*) 'agd_SetupFromHCO: Starting...'

    initialized = .true.

    !
    !- 1.  Check the hco structures
    !
    hco_ext  => hco_ext_in
    if( present(hco_core_opt) ) then
      hco_core => hco_core_opt
    else
      hco_core => hco_ext_in
    end if

    if ( (.not. hco_core%initialized) .or.  & 
         (.not. hco_ext%initialized) ) then
      write(*,*)
      write(*,*) 'agd_SetupFromHCO: At least one hco structure was not initilzed'
      write(*,*) 'hco_core = ', hco_core%initialized
      write(*,*) 'hco_ext = ', hco_ext%initialized
      call utl_abort('agd_SetupFromHCO: abort')
    end if

    if ( (hco_core%global) .or.  & 
         (hco_ext%global) ) then
      write(*,*)
      write(*,*) 'agd_SetupFromHCO: At least one hco structure is from a global grid, skipping rest of setup'
      write(*,*) 'hco_core = ', hco_core%global
      write(*,*) 'hco_ext = ', hco_ext%global
      return
    end if

    !
    !- 2.  Dimension settings and Memory allocation
    !

    !- 2.1 Dimensions
    ni_core = hco_core%ni
    nj_core = hco_core%nj

    ni_ext = hco_ext%ni
    nj_ext = hco_ext%nj

    ext_i  = ni_ext - ni_core
    ext_j  = nj_ext - nj_core

    if ( ext_i == 0 .and. ext_j == 0 ) then
      write(*,*)
      write(*,*) 'agd_SetupFromHCO: LAM core and extended grids are identical'
    else if ( ext_i < 10 .or. ext_j < 10 ) then
      write(*,*)
      write(*,*) 'agd_SetupFromHCO: LAM domain extension is less than 10 gridpoints'
      write(*,*) ' ext_i = ', ext_i,' ext_j = ', ext_j
    end if

    ni_ext_per = ni_ext + 1  ! Fields will be periodic (i = 1 repeated) at ni_ext+1
    nj_ext_per = nj_ext + 1  ! Fields will be periodic (j = 1 repeated) at nj_ext+1

    !- 2.2 First Make sure we have a uniform grid in x and y direction

    !- 2.2.1 Analysis core grid
    dlon_ref = hco_core%lon(2) - hco_core%lon(1)
    do i = 2, ni_core
      dlon_test = hco_core%lon(i) - hco_core%lon(i-1)
      if ( (dlon_test - dlon_ref) > dlon_ref/100.0d0 ) then
        write(*,*)
        write(*,*) 'agd_SetupFromHCO: Core grid spacing is not uniform in x-direction'
        write(*,*) ' i         = ', i
        write(*,*) ' dlon      = ', dlon_test
        write(*,*) ' dlon ref  = ', dlon_ref
        call utl_abort('agd_SetupFromHCO')
      end if
    end do

    dlat_ref = hco_core%lat(2) - hco_core%lat(1)
    do j = 2, nj_core
      dlat_test = hco_core%lat(j) - hco_core%lat(j-1) 
      if ( (dlat_test - dlat_ref) > dlat_ref/100.0d0 ) then
        write(*,*)
        write(*,*) 'agd_SetupFromHCO: Core grid spacing is not uniform in x-direction'
        write(*,*) ' j         = ', j
        write(*,*) ' dlat      = ', dlon_test
        write(*,*) ' dlat ref  = ', dlon_ref
        call utl_abort('agd_SetupFromHCO') 
      end if
    end do

    !- 2.2.2 Extended Analysis grid
    dlon_ref = hco_ext%lon(2) - hco_ext%lon(1)
    do i = 2, ni_ext
      dlon_test = hco_ext%lon(i) - hco_ext%lon(i-1)
      if ( (dlon_test - dlon_ref) > dlon_ref/100.0d0 ) then
        write(*,*)
        write(*,*) 'agd_SetupFromHCO: Extended grid spacing is not uniform in x-direction'
        write(*,*) ' i         = ', i
        write(*,*) ' dlon      = ', dlon_test
        write(*,*) ' dlon ref  = ', dlon_ref
        call utl_abort('agd_SetupFromHCO')
      end if
    end do

    dlat_ref = hco_ext%lat(2) - hco_ext%lat(1)
    do j = 2, nj_ext
      dlat_test = hco_ext%lat(j) - hco_ext%lat(j-1)
      if ( (dlat_test - dlat_ref) > dlat_ref/100.0d0 ) then
        write(*,*)
        write(*,*) 'agd_SetupFromHCO: Extended grid spacing is not uniform in x-direction'
        write(*,*) ' j         = ', j
        write(*,*) ' dlat      = ', dlon_test
        write(*,*) ' dlat ref  = ', dlon_ref
        call utl_abort('agd_SetupFromHCO') 
      end if
    end do

    !- 2.3 Dimensions for variables needed to be symmetric
    istart = -4          ! 5 gridpoints padding in West direction
    iend   = ni_ext + 4  ! 4 gridpoints padding in East direction
    jstart = -4          ! 5 gridpoints padding in South direction
    jend   = nj_ext + 4  ! 4 gridpoints padding in North direction

    !- 2.4 Allocations
    allocate(glmf%rlon(istart:iend))
    allocate(glmf%rlat(jstart:jend))
    allocate(rlath    (jstart:jend))
    allocate(rrcos    (jstart:jend))
    !allocate(rrcosh   (jstart:jend))
    allocate(rdmu     (jstart:jend))
    allocate(rdmuh    (jstart:jend))
    allocate(r1mmu2   (jstart:jend))
    allocate(r1mmu2h  (jstart:jend))

    !
    !- 3.  Set Lat-Lon of the computational grid
    !

    !  3.1 Set (constant) Grid spacing
    glmf%rdlon = (hco_ext%lon(2) - hco_ext%lon(1))
    glmf%rdlat = (hco_ext%lat(2) - hco_ext%lat(1))

    !- 3.2 Lat-Lon

    !  3.2.1 Extract the lat-lon from the extended grid
    glmf%rlon(1:ni_ext) = hco_ext%lon(1:ni_ext)
    glmf%rlat(1:nj_ext) = hco_ext%lat(1:nj_ext)

    !  3.2.2 Extend to the full computational grid (larger than the extended grid!)

    ! West
    do i = istart, 0
      glmf%rlon(i) = glmf%rlon(1) + (i-1) * glmf%rdlon
    end do
    ! East
    do i = ni_ext+1, iend
      glmf%rlon(i) = glmf%rlon(ni_ext) + (i-ni_ext) * glmf%rdlon
    end do
    ! North
    do j = nj_ext+1, jend
      glmf%rlat(j) = glmf%rlat(nj_ext) + (j-nj_ext) * glmf%rdlat
    end do
    ! South
    do j = jstart, 0
      glmf%rlat(j) = glmf%rlat(1) + (j-1) * glmf%rdlat
    end do

    !- 3.2 Half Lat-Lon
    do j = jstart+1, jend-1
      rlath(j) = ( glmf%rlat(j) + glmf%rlat(j+1) ) / 2.0d0
    end do

    !
    !- 4. Set Metric Factors
    !

    !- 4.1  Compute local factors
    do j = jstart+1, jend-2
      rrcos  (j) = 1.0d0 / cos(glmf%rlat (j))
      !rrcosh (j) = 1.0d0 / cos(rlath(j))
      rdmu   (j) = sin(glmf%rlat (j+1)) - sin(glmf%rlat (j))
      rdmuh  (j) = sin(rlath(j+1)) - sin(rlath(j))
      r1mmu2 (j) = (cos(glmf%rlat (j)))**2
      r1mmu2h(j) = (cos(rlath(j)))**2
    end do

    !- 4.2  Bi-periodize and symmetrize Metric coefficients
    call agd_mach(rdmu   (1:nj_ext), & ! INOUT
                   1, nj_ext,1)        ! IN
    call agd_mach(rdmuh  (1:nj_ext), & ! INOUT
                   1, nj_ext,1)        ! IN 
    call agd_mach(r1mmu2 (1:nj_ext), & ! INOUT
                   1, nj_ext,1)        ! IN
    call agd_mach(r1mmu2h(1:nj_ext), & ! INOUT
                   1, nj_ext,1)        ! IN

    call symmetrize_coef(rdmu   ) ! INOUT
    call symmetrize_coef(rdmuh  ) ! INOUT
    call symmetrize_coef(r1mmu2 ) ! INOUT
    call symmetrize_coef(r1mmu2h) ! INOUT

    !- 4.3 Compute global factors
    glmf%dx  = 1.d0 / (EC_RA * glmf%rdlon)

    allocate(glmf%cos2   ( 0:nj_ext+1))
    allocate(glmf%cos2h  ( 0:nj_ext+1))
    do j = 0, nj_ext+1
      glmf%cos2 (j) = r1mmu2 (j) / (EC_RA * rdmuh(j-1))
      glmf%cos2h(j) = r1mmu2h(j) / (EC_RA * rdmu (j  ))
    end do

    allocate(glmf%idmuh  (0:nj_ext-1))
    allocate(glmf%idmu   (1:nj_ext  ))
    allocate(glmf%cos2vd (1:nj_ext  ))
    allocate(glmf%cos2hvd(1:nj_ext  ))

    do j = 1, nj_ext
      glmf%idmuh  (j-1) = 1.d0 / rdmuh(j-1)
      glmf%idmu   (j)   = 1.d0 / rdmu   (j)
      glmf%cos2vd (j)   = 1.d0 / r1mmu2 (j)
      glmf%cos2hvd(j)   = 1.d0 / r1mmu2h(j)
    end do

    !
    ! Conversion of wind images to physical winds and vice-versa
    ! N.B.: Those are geometrical factors of the COMPUTATIONAL grid 
    !        ==> only computational latitude variation...
    !
    allocate(glmf%conphy(nj_ext))
    allocate(glmf%conima(nj_ext))

    call agd_mach(rrcos(1:nj_ext),  & ! INOUT
                  1, nj_ext,1)        ! IN

    do j = 1, nj_ext
      glmf%conphy(j) = rrcos(j)                ! to go from Wind Images to true winds
      glmf%conima(j) = 1.0d0 / glmf%conphy(j)  ! to go from true winds to Wind Images
    end do

    !
    !- 5. MPI partitionning
    !
    if ( mpi_nprocs /= 0 ) then
       call mpivar_setup_lonbands(ni_ext,                      & ! IN
                                  lonPerPE, lonPerPEmax, myLonBeg, myLonEnd ) ! OUT

       call mpivar_setup_latbands(nj_ext,                      & ! IN
                                  latPerPE, latPerPEmax, myLatBeg, myLatEnd ) ! OUT
    else
       ! This option is needed for the biper program
       lonPerPE = ni_ext
       myLonBeg = 1
       myLonEnd = ni_ext
       latPerPE = nj_ext
       myLatBeg = 1
       myLatEnd = nj_ext
       write(*,*)
       write(*,*) 'WARNING: the module will be executed in NO-MPI MODE!'
    end if

    !
    !- 6.  Ending
    !
    deallocate(rlath  )
    deallocate(rrcos  )
    !deallocate(rrcosh )
    deallocate(rdmu   )
    deallocate(rdmuh  )
    deallocate(r1mmu2 )
    deallocate(r1mmu2h)

    write(*,*)
    write(*,*) 'agd_SetupFromHCO: Done!'

  end subroutine agd_SetupFromHCO

  !--------------------------------------------------------------------------
  ! symmetrize_coef
  !--------------------------------------------------------------------------
  subroutine symmetrize_coef(coef_inout)
    !
    !:Purpose: Extend symmetrically Metric coefficients.
    !
    implicit none

    real(8), intent(inout) :: coef_inout(jstart:jend)

    integer j
    
    do j = jstart, 0
      coef_inout(j) = coef_inout(nj_ext+j)
    end do

    do j = nj_ext+1, jend
      coef_inout(j) = coef_inout(j-nj_ext) 
    end do

  end subroutine symmetrize_coef

  !--------------------------------------------------------------------------
  ! agd_PsiChiToUV
  !--------------------------------------------------------------------------
  subroutine agd_PsiChiToUV(psi, chi, uphy, vphy, nk)
    implicit none

    integer,          intent(in)  :: nk

    real(8),          intent(in)  :: psi(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nk)
    real(8),          intent(in)  :: chi(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nk)

    real(8),          intent(out) :: uphy(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nk)
    real(8),          intent(out) :: vphy(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nk)

    real(8), allocatable :: psi_ext(:,:,:)
    real(8), allocatable :: chi_ext(:,:,:)
    real(8), allocatable :: uimg(:,:,:)
    real(8), allocatable :: vimg(:,:,:)
    real(8), allocatable :: uimgs(:,:,:)
    real(8), allocatable :: vimgs(:,:,:)

    integer :: i,j,k

    if ( hco_ext%global ) then
      call utl_abort('agd_PsiChiToUV: Not compatible with global grid')
    endif

    if ( .not. initialized ) then
      call utl_abort('agd_PsiChiToUV: AnalysisGrid not initialized')
    endif

    allocate(psi_ext( (myLonBeg-1):(myLonEnd+1), (myLatBeg-1):(myLatEnd+1), nk))
    allocate(chi_ext( (myLonBeg-1):(myLonEnd+1), (myLatBeg-1):(myLatEnd+1), nk))
    allocate(uimg   ( myLonBeg:myLonEnd        , myLatBeg:myLatEnd        , nk))
    allocate(vimg   ( myLonBeg:myLonEnd        , myLatBeg:myLatEnd        , nk))
    allocate(uimgs  ( (myLonBeg-1):myLonEnd    , myLatBeg:(myLatEnd+1)    , nk))
    allocate(vimgs  ( myLonBeg:(myLonEnd+1)    , (myLatBeg-1):myLatEnd    , nk))

    !
    !- 1.  Symmetrize
    !
    call symmetrize( psi_ext,                                        & ! OUT
                     psi, myLonBeg, myLonEnd, myLatBeg, myLatEnd, nk ) ! IN
    call symmetrize( chi_ext,                                        & ! OUT
                     chi, myLonBeg, myLonEnd, myLatBeg, myLatEnd, nk ) ! IN

    !
    !- 2.  Compute Wind on staggered grid
    !

    !$OMP PARALLEL DO PRIVATE (k,j,i)
    do k = 1, nk
      !- 2.1 u-wind component
       do j = myLatBeg, myLatEnd+1
        do i = myLonBeg-1, myLonEnd
           uimgs(i,j,k) =   glmf%dx      * ( chi_ext(i+1,j,k) - chi_ext(i,j,k)   ) &
                          - glmf%cos2(j) * ( psi_ext(i,j,k)   - psi_ext(i,j-1,k) )
        end do
      end do
      !- 2.2 v-wind component
      do j = myLatBeg-1, myLatEnd
        do i = myLonBeg, myLonEnd+1
           vimgs(i,j,k) =   glmf%cos2h(j) * ( chi_ext(i,j+1,k) - chi_ext(i,j,k)   ) &
                          + glmf%dx       * ( psi_ext(i,j,k)   - psi_ext(i-1,j,k) )
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    !
    !- 3.  Move to collocated (scalar) grid
    !
    call uvStagToColloc( uimgs, vimgs,                              & ! IN
                         uimg , vimg ,                              & ! OUT
                         myLonBeg, myLonEnd, myLatBeg, myLatEnd, nk ) ! IN

    !
    !- 4.  Convert Wind images to Physical (true) winds
    !
    !$OMP PARALLEL DO PRIVATE (j)
    do j = myLatBeg, myLatEnd
      uphy(:,j,:) =  glmf%conphy(j) * uimg(:,j,:)
      vphy(:,j,:) =  glmf%conphy(j) * vimg(:,j,:)
    end do
    !$OMP END PARALLEL DO

    deallocate(psi_ext)
    deallocate(chi_ext)
    deallocate(uimg)
    deallocate(vimg)
    deallocate(uimgs)
    deallocate(vimgs)

  end subroutine agd_PsiChiToUV

  !--------------------------------------------------------------------------
  ! agd_PsiChiToUVAdj
  !--------------------------------------------------------------------------
  subroutine agd_PsiChiToUVAdj(psi, chi, uphy, vphy, nk)
    implicit none

    integer,          intent(in)  :: nk

    real(8),          intent(out) :: psi(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nk)
    real(8),          intent(out) :: chi(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nk)
    
    real(8),          intent(in)  :: uphy(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nk)
    real(8),          intent(in)  :: vphy(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nk)

    real(8), allocatable :: psi_ext(:,:,:)
    real(8), allocatable :: chi_ext(:,:,:)
    real(8), allocatable :: uimg(:,:,:)
    real(8), allocatable :: vimg(:,:,:)
    real(8), allocatable :: uimgs(:,:,:)
    real(8), allocatable :: vimgs(:,:,:)

    integer :: i,j,k

    if ( hco_ext%global ) then
      call utl_abort('agd_PsiChiToUVAdj: Not compatible with global grid')
    endif

    if ( .not. initialized ) then
      call utl_abort('agd_PsiChiToUV: AnalysisGrid not initialized')
    endif

    allocate(psi_ext( (myLonBeg-1):(myLonEnd+1), (myLatBeg-1):(myLatEnd+1), nk))
    allocate(chi_ext( (myLonBeg-1):(myLonEnd+1), (myLatBeg-1):(myLatEnd+1), nk))
    allocate(uimg   ( myLonBeg:myLonEnd        , myLatBeg:myLatEnd        , nk))
    allocate(vimg   ( myLonBeg:myLonEnd        , myLatBeg:myLatEnd        , nk))
    allocate(uimgs  ( (myLonBeg-1):myLonEnd    , myLatBeg:(myLatEnd+1)    , nk))
    allocate(vimgs  ( myLonBeg:(myLonEnd+1)    , (myLatBeg-1):myLatEnd    , nk))

    !
    !- 4.  Convert Physical (true) winds to Wind images
    !
    !$OMP PARALLEL DO PRIVATE (j)
    do j = myLatBeg, myLatEnd
      uimg(:,j,:) =  glmf%conphy(j) * uphy(:,j,:)
      vimg(:,j,:) =  glmf%conphy(j) * vphy(:,j,:)
    end do
    !$OMP END PARALLEL DO

    !
    !- 3.  Move to stagerred grid
    !
    uimgs(:,:,:) = 0.d0
    vimgs(:,:,:) = 0.d0

    call uvStagToCollocAdj( uimgs, vimgs,                              & ! OUT
                            uimg , vimg ,                              & ! IN
                            myLonBeg, myLonEnd, myLatBeg, myLatEnd, nk ) ! IN

    !
    !- 2.  Compute Psi-Chi
    !
    chi_ext(:,:,:) = 0.d0
    psi_ext(:,:,:) = 0.d0

    !$OMP PARALLEL DO PRIVATE (k,j,i)
    do k = 1, nk
      !- 2.2 from v-wind component
      do j = myLatEnd, myLatBeg-1, -1
        do i = myLonEnd+1, myLonBeg, -1
          chi_ext(i  ,j+1,k) = chi_ext(i  ,j+1,k) + vimgs(i,j,k) * glmf%cos2h(j)
          chi_ext(i  ,j  ,k) = chi_ext(i  ,j  ,k) - vimgs(i,j,k) * glmf%cos2h(j)
          psi_ext(i  ,j  ,k) = psi_ext(i  ,j  ,k) + vimgs(i,j,k) * glmf%dx
          psi_ext(i-1,j  ,k) = psi_ext(i-1,j  ,k) - vimgs(i,j,k) * glmf%dx
        end do
      end do
      !- 2.1 from u-wind component
      do j = myLatEnd+1, myLatBeg, -1
        do i = myLonEnd, myLonBeg-1, -1
          chi_ext(i+1,j  ,k) = chi_ext(i+1,j  ,k) + uimgs(i,j,k) * glmf%dx
          chi_ext(i  ,j  ,k) = chi_ext(i  ,j  ,k) - uimgs(i,j,k) * glmf%dx
          psi_ext(i  ,j  ,k) = psi_ext(i  ,j  ,k) - uimgs(i,j,k) * glmf%cos2(j)
          psi_ext(i  ,j-1,k) = psi_ext(i  ,j-1,k) + uimgs(i,j,k) * glmf%cos2(j)
        end do
     end do
    end do
    !$OMP END PARALLEL DO

    !
    !- 1.  De-Symmetrize
    !
    chi(:,:,:) = 0.d0
    psi(:,:,:) = 0.d0

    call symmetrizeAdj( psi_ext,                                   & ! IN
                        psi,                                       & ! OUT
                        myLonBeg, myLonEnd, myLatBeg, myLatEnd, nk ) ! IN
    call symmetrizeAdj( chi_ext,                                   & ! IN
                        chi,                                       & ! OUT
                        myLonBeg, myLonEnd, myLatBeg, myLatEnd, nk ) ! IN

    deallocate(psi_ext)
    deallocate(chi_ext)
    deallocate(uimg)
    deallocate(vimg)
    deallocate(uimgs)
    deallocate(vimgs)

  end subroutine agd_PsiChiToUVAdj

  !--------------------------------------------------------------------------
  ! uvStagToColloc
  !--------------------------------------------------------------------------
  subroutine uvStagToColloc(uStag, vStag, uColloc, vColloc, iBeg, iEnd, jBeg, jEnd , nk)
    implicit none

    integer,          intent(in)  :: iBeg, iEnd, jBeg, JEnd, nk
    real(8),          intent(out) :: uColloc(iBeg:iEnd  ,jBeg  :jEnd  ,nk)
    real(8),          intent(out) :: vColloc(iBeg:iEnd  ,jBeg  :jEnd  ,nk)
    real(8),          intent(in)  :: uStag  (iBeg-1:iEnd,jBeg  :jEnd+1,nk)
    real(8),          intent(in)  :: vStag  (iBeg:iEnd+1,jBeg-1:jEnd  ,nk)

    integer :: i,j,k

    !$OMP PARALLEL DO PRIVATE (k,j,i)
    do k = 1, nk
      do j = jBeg, jEnd
        do i = iBeg, iEnd
          uColloc(i,j,k) = ( uStag(i-1,j  ,k) + uStag(i,j,k) ) / 2.0d0
          vColloc(i,j,k) = ( vStag(i  ,j-1,k) + vStag(i,j,k) ) / 2.0d0
        end do
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine uvStagToColloc

  !--------------------------------------------------------------------------
  ! uvStagToCollocAdj
  !--------------------------------------------------------------------------
  subroutine uvStagToCollocAdj(uStag, vStag, uColloc, vColloc, iBeg, iEnd, jBeg, jEnd , nk)
    implicit none

    integer,          intent(in)  :: iBeg, iEnd, jBeg, jEnd, nk
    real(8),          intent(in)  :: uColloc(iBeg:iEnd  ,jBeg  :jEnd  ,nk)
    real(8),          intent(in)  :: vColloc(iBeg:iEnd  ,jBeg  :jEnd  ,nk)
    real(8),          intent(out) :: uStag  (iBeg-1:iEnd,jBeg  :jEnd+1,nk)
    real(8),          intent(out) :: vStag  (iBeg:iEnd+1,jBeg-1:jEnd  ,nk)

    integer :: i, j, k

    !$OMP PARALLEL DO PRIVATE (k,j,i)
    do k = 1, nk
      do j = jEnd, jBeg, -1
        do i = iEnd, iBeg, -1
          uStag(i-1,j  ,k) = uStag(i-1,j  ,k) + uColloc(i,j,k) / 2.0d0
          uStag(i  ,j  ,k) = uStag(i  ,j  ,k) + uColloc(i,j,k) / 2.0d0
          vStag(i  ,j-1,k) = vStag(i  ,j-1,k) + vColloc(i,j,k) / 2.0d0
          vStag(i  ,j  ,k) = vStag(i  ,j  ,k) + vColloc(i,j,k) / 2.0d0
        end do
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine uvStagToCollocAdj

  !--------------------------------------------------------------------------
  ! Symmetrize
  !--------------------------------------------------------------------------
  subroutine symmetrize(field_out, field_in, iBeg, iEnd, jBeg, jEnd, nk)
    !
    !:Purpose: Extend symmetrically outside 1 grid point all around LAM-boundary
    !          ready for finite differences
    !
    implicit none
    integer, intent(in) :: iBeg, iEnd, jBeg, jEnd, nk
    real(8), intent(out):: field_out(iBeg-1:iEnd+1,jBeg-1:jEnd+1, nk)
    real(8), intent(in) :: field_in(iBeg:iEnd,jBeg:jEnd,nk)

    real(8), allocatable :: field_8(:,:,:)

    integer :: ni,nj

    ni = iEnd-iBeg+1
    nj = jEnd-jBeg+1

    allocate(field_8(0:ni+1,0:nj+1, nk))
    field_8(:,:,:) = 0.d0

    field_8(1:ni,1:nj,:) = field_in(iBeg:iEnd,jBeg:jEnd,:)

    call RPN_COMM_xch_halo_8(field_8,                & ! INOUT
                             0,ni+1,0,nj+1,ni,nj,nk, & ! IN
                             1,1,.true.,.true.,ni,0)   ! IN

    field_out(iBeg-1:iEnd+1,jBeg-1:jEnd+1,:) = field_8(0:ni+1,0:nj+1,:)
    deallocate(field_8)

  end subroutine symmetrize

  !--------------------------------------------------------------------------
  ! SymmetrizeAdj
  !--------------------------------------------------------------------------
  subroutine symmetrizeAdj(field_in, field_out, iBeg, iEnd, jBeg, jEnd, nk)
    !
    !:Purpose: Adjoint of sub. symmetrize.
    !
    implicit none
    integer, intent(in)  :: iBeg, iEnd, jBeg, jEnd, nk
    real(8), intent(in)  :: field_in(iBeg-1:iEnd+1,jBeg-1:jEnd+1, nk)
    real(8), intent(out) :: field_out(iBeg:iEnd,jBeg:jEnd,nk)

    real(8), allocatable :: field_8(:,:,:)

    integer :: i,j,k,ni,nj

    ni = iEnd-iBeg+1
    nj = jEnd-jBeg+1

    allocate(field_8(0:ni+1,0:nj+1, nk))
    field_8(0:ni+1,0:nj+1,:) = field_in(iBeg-1:iEnd+1,jBeg-1:jEnd+1,:)

    call RPN_COMM_adj_halo8(field_8,                 & ! INOUT
                            0,ni+1,0,nj+1,ni,nj,nk,  & ! IN
                            1,1,.true.,.true.,ni,0)    ! IN

    !$OMP PARALLEL DO PRIVATE (k,j,i)
    do k = 1, nk
       do j = jBeg, jEnd
          do i = iBeg, iEnd
             field_out(i,j,k) = field_out(i,j,k) + field_8(i-iBeg+1,j-jBeg+1,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    deallocate(field_8)

  end subroutine symmetrizeAdj

  !--------------------------------------------------------------------------
  ! agd_Mach
  !--------------------------------------------------------------------------
  subroutine agd_mach(gd,ni,nj,nk)
    !
    !:Purpose:  [to be completed]
    !
    !:Arguments:
    !       :ni: Maximum I-dimension where the input array is assumed to carry
    !            information.  Will be used as I-limit where backward
    !            derivatives will be evaluated
    !
    !       :nj: Maximum J-dimension where the input array is assumed to carry
    !            information. Will be used as J-limit where backward derivatives
    !            will be evaluated
    !
    implicit none

    integer,          intent(in)    :: ni, nj, nk
    real(8),          intent(inout) :: gd(ni,nj,nk)

    integer :: istart, jstart
    integer :: i,j,k

    real(8) :: con, xp, yp, a0, a1, b1, b2
    real(8) :: deriv_istart, deriv_jstart, deriv_i0, deriv_j0, del

    if ( hco_ext%global ) then
      call utl_abort('agd_Mach: Not compatible with global grid')
    endif

    if ( .not. initialized ) then
      call utl_abort('agd_Mach: AnalysisGrid not initialized')
    endif

    if ( (ni /= 1 .and. ni /= ni_ext) .or. &
         (nj /= 1 .and. nj /= nj_ext) ) then
      call utl_abort('agd_Mach: Invalid Dimensions')
    end if

    !$OMP PARALLEL
    !$OMP DO PRIVATE (k,j,i,istart,jstart,con,a0,a1,b1,b2,del,deriv_istart,deriv_i0,deriv_jstart,deriv_j0,xp,yp)
    do k = 1, nk
    !
    !- 1.  Periodicized in x-direction from ni_core to ni
    !    
    if ( ni > 1 ) then
      istart = ni_core - 2  ! I-limit where backward derivatives will be evaluated
      con = 1.d0 / ( glmf%rdlon * MPC_DEGREES_PER_RADIAN_R8 * 111.d+3)
      do j = 1, nj
        a0           = 0.5d0 * ( gd(istart,j,k)  + gd(1,j,k)       )
        a1           = 0.5d0 * ( gd(istart,j,k)  - gd(1,j,k)       )
        deriv_istart = con   * ( gd(istart,j,k)  - gd(istart-1,j,k))
        deriv_i0     = con   * ( gd(2,j,k)       - gd(1,j,k)       )
        b1           = 0.5d0 * ( deriv_istart  - deriv_i0       )
        b2           = 0.25d0* ( deriv_istart  + deriv_i0       )
        del          = real(ni_ext_per-istart,8)
        do i = istart, ni
          xp = MPC_PI_R8 * real(i-istart,8) / del
          gd(i,j,k) = a0 + a1*cos(xp) + b1*sin(xp) + b2*sin(2.d0*xp)
        end do
      end do
    end if
    !
    !- 2.  Periodicized in y-direction from nj_core to nj
    !
    if ( nj > 1 ) then
      jstart = nj_core - 2
      con = 1.d0 / (glmf%rdlat * MPC_DEGREES_PER_RADIAN_R8 * 111.d+3)
      do i = 1, ni
        a0           = 0.5d0 * ( gd(i,jstart,k) + gd(i,1,k)       )
        a1           = 0.5d0 * ( gd(i,jstart,k) - gd(i,1,k)       )
        deriv_jstart = con   * ( gd(i,jstart,k) - gd(i,jstart-1,k))
        deriv_j0     = con   * ( gd(i,2,k)      - gd(i,1,k)       )
        b1           = 0.5d0 * ( deriv_jstart - deriv_j0     )
        b2           = 0.25d0* ( deriv_jstart + deriv_j0     )
        del          = real(nj_ext_per-jstart,8)
        do j = jstart , nj
          yp = MPC_PI_R8 * real(j-jstart,8) / del
          gd(i,j,k) = a0 + a1*cos(yp) + b1*sin(yp) + b2*sin(2.d0*yp)
        end do
      end do
    end if

    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine agd_mach

  !--------------------------------------------------------------------------
  ! agd_Mach_r4
  !--------------------------------------------------------------------------
  subroutine agd_mach_r4(gd,ni,nj,nk)
    !:Purpose:  [to be completed]
    !
    !:Arguments:
    !       :ni:  Maximum I-dimension where the input array is assumed to carry
    !             information. Will be used as I-limit where backward
    !             derivatives will be evaluated
    !       :nj:  Maximum J-dimension where the input array is assumed to carry
    !             information. Will be used as J-limit where backward
    !             derivatives will be evaluated
    !
    implicit none

    integer,          intent(in)    :: ni, nj, nk
    real(4),          intent(inout) :: gd(ni,nj,nk)

    integer :: istart, jstart
    integer :: i,j,k

    real(4) :: con, xp, yp, a0, a1, b1, b2
    real(4) :: deriv_istart, deriv_jstart, deriv_i0, deriv_j0, del

    if ( hco_ext%global ) then
      call utl_abort('agd_Mach_r4: Not compatible with global grid')
    endif

    if ( .not. initialized ) then
      call utl_abort('agd_Mach_r4: AnalysisGrid not initialized')
    endif

    if ( (ni /= 1 .and. ni /= ni_ext) .or. &
         (nj /= 1 .and. nj /= nj_ext) ) then
      call utl_abort('agd_Mach_r4 : Invalid Dimensions')
    end if

    !$OMP PARALLEL
    !$OMP DO PRIVATE (k,j,i,istart,jstart,con,a0,a1,b1,b2,del,deriv_istart,deriv_i0,deriv_jstart,deriv_j0,xp,yp)
    do k = 1, nk
    !
    !- 1.  Periodicized in x-direction from ni_core to ni
    !    
    if ( ni > 1 ) then
      istart = ni_core - 2  ! I-limit where backward derivatives will be evaluated
      con = 1.0 / ( glmf%rdlon * MPC_DEGREES_PER_RADIAN_R4 * 111.e+3)
      do j = 1, nj
        a0           = 0.50 * ( gd(istart,j,k)  + gd(1,j,k)       )
        a1           = 0.50 * ( gd(istart,j,k)  - gd(1,j,k)       )
        deriv_istart = con   * ( gd(istart,j,k)  - gd(istart-1,j,k))
        deriv_i0     = con   * ( gd(2,j,k)       - gd(1,j,k)       )
        b1           = 0.50 * ( deriv_istart  - deriv_i0       )
        b2           = 0.250* ( deriv_istart  + deriv_i0       )
        del          = real(ni_ext_per-istart,8)
        do i = istart, ni
          xp = MPC_PI_R4 * real(i-istart,8) / del
          gd(i,j,k) = a0 + a1*cos(xp) + b1*sin(xp) + b2*sin(2.0*xp)
        end do
      end do
    end if
    !
    !- 2.  Periodicized in y-direction from nj_core to nj
    !
    if ( nj > 1 ) then
      jstart = nj_core - 2
      con = 1.0 / (glmf%rdlat * MPC_DEGREES_PER_RADIAN_R4 * 111.e+3)
      do i = 1, ni
        a0           = 0.50 * ( gd(i,jstart,k) + gd(i,1,k)       )
        a1           = 0.50 * ( gd(i,jstart,k) - gd(i,1,k)       )
        deriv_jstart = con   * ( gd(i,jstart,k) - gd(i,jstart-1,k))
        deriv_j0     = con   * ( gd(i,2,k)      - gd(i,1,k)       )
        b1           = 0.50 * ( deriv_jstart - deriv_j0     )
        b2           = 0.250* ( deriv_jstart + deriv_j0     )
        del          = real(nj_ext_per-jstart,8)
        do j = jstart , nj
          yp = MPC_PI_R4 * real(j-jstart,8) / del
          gd(i,j,k) = a0 + a1*cos(yp) + b1*sin(yp) + b2*sin(2.0*yp)
        end do
      end do
    end if

    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine agd_mach_r4

  !--------------------------------------------------------------------------
  ! agd_UVToVortDiv
  !--------------------------------------------------------------------------
  subroutine agd_UVToVortDiv(Vorticity, Divergence, uphy, vphy, nk)
    implicit none

    integer,          intent(in)  :: nk

    real(8),          intent(out)  :: Vorticity (myLonBeg:myLonEnd,myLatBeg:myLatEnd,nk)
    real(8),          intent(out)  :: Divergence(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nk)

    real(8),          intent(in) :: uphy(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nk)
    real(8),          intent(in) :: vphy(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nk)

    real(8), allocatable :: uimg(:,:,:)
    real(8), allocatable :: vimg(:,:,:)
    real(8), allocatable :: uimg_sym(:,:,:)
    real(8), allocatable :: vimg_sym(:,:,:)

    integer :: i,j,k

    if ( hco_ext%global ) then
      call utl_abort('agd_UVToVortDiv: Not compatible with global grid')
    endif

    if ( .not. initialized ) then
      call utl_abort('agd_UVToVortDiv: AnalysisGrid not initialized')
    endif

    allocate(uimg_sym( (myLonBeg-1):(myLonEnd+1), (myLatBeg-1):(myLatEnd+1), nk))
    allocate(vimg_sym( (myLonBeg-1):(myLonEnd+1), (myLatBeg-1):(myLatEnd+1), nk))
    allocate(uimg    ( myLonBeg:myLonEnd        , myLatBeg:myLatEnd        , nk))
    allocate(vimg    ( myLonBeg:myLonEnd        , myLatBeg:myLatEnd        , nk))

    !
    !- 2.  Convert Physical (true) winds to Wind images
    !
    !$OMP PARALLEL DO PRIVATE (j)
    do j = myLatBeg, myLatEnd
          uimg(:,j,:) =  glmf%conima(j) * uphy(:,j,:)
          vimg(:,j,:) =  glmf%conima(j) * vphy(:,j,:)
    end do
    !$OMP END PARALLEL DO

    !
    !- 3.  Symmetrize
    !
    call symmetrize( uimg_sym,                                        & ! OUT
                     uimg, myLonBeg, myLonEnd, myLatBeg, myLatEnd, nk ) ! IN
    call symmetrize( vimg_sym,                                        & ! OUT
                     vimg, myLonBeg, myLonEnd, myLatBeg, myLatEnd, nk ) ! IN
    deallocate(uimg,vimg)

    !
    !- 4.  Compute Vorticity and Divergence
    !

    !$OMP PARALLEL DO PRIVATE (k,j,i)
     do k = 1, nk
       do j = myLatBeg, myLatEnd
         do i = myLonBeg, myLonEnd
           vorticity (i,j,k) = &
                glmf%cos2hvd(j) * glmf%dx       * (vimg_sym(i+1,j,k)-vimg_sym(i,j,k)) - &
                       EC_R1SA  * glmf%idmu(j)  * (uimg_sym(i,j+1,k)-uimg_sym(i,j,k))

           divergence(i,j,k) = &
                glmf%cos2vd(j) * glmf%dx        * (uimg_sym(i,j,k)-uimg_sym(i-1,j,k)) + &
                       EC_R1SA * glmf%idmuh(j-1)* (vimg_sym(i,j,k)-vimg_sym(i,j-1,k))
         end do
       end do
     end do
     !$OMP END PARALLEL DO

    deallocate(uimg_sym)
    deallocate(vimg_sym)

  end subroutine agd_UVToVortDiv

  !--------------------------------------------------------------------------
  ! agd_createLamTemplateGrids
  !--------------------------------------------------------------------------
  subroutine agd_createLamTemplateGrids(templateFileName, hco_core, vco, &
                                        grd_ext_x, grd_ext_y)
    implicit none

    character(len=*), intent(in) :: templateFileName
    type(struct_hco) :: hco_core
    type(struct_vco) :: vco
    integer         , intent(in) :: grd_ext_x
    integer         , intent(in) :: grd_ext_y

    integer :: ni_ext, nj_ext, i, j, lev, ni, nj, nk
    integer :: iun = 0
    integer :: ier, fnom, fstouv, fstfrm, fclos, fstecr

    real(8), allocatable :: Field2d(:,:)
    real(8), allocatable :: lat_ext(:)
    real(8), allocatable :: lon_ext(:)

    real(4), allocatable :: dummy2D(:,:)

    real(8) :: dlat, dlon
    real(4) :: work

    integer :: dateo,npak,status
    integer :: ip1,ip2,ip3,deet,npas,datyp,ig1,ig2,ig3,ig4
    integer :: ig1_tictac,ig2_tictac,ig3_tictac,ig4_tictac

    character(len=1)  :: grtyp
    character(len=2)  :: typvar
    character(len=12) :: etiket

    !
    !- 1.  Opening the output template file
    !
    ier = fnom(iun, trim(templateFileName), 'RND', 0)
    ier = fstouv(iun, 'RND')

    npak     = -32

    !
    !- 2.  Writing the core grid (Ensemble) template
    !

    !- 2.1 Tic-Tac
    deet     =  0
    ip1      =  hco_core%ig1
    ip2      =  hco_core%ig2
    ip3      =  hco_core%ig3
    npas     =  0
    datyp    =  1
    grtyp    =  hco_core%grtypTicTac
    typvar   = 'X'
    etiket   = 'COREGRID'
    dateo =  0

    call cxgaig ( grtyp,                                          & ! IN
         ig1_tictac, ig2_tictac, ig3_tictac, ig4_tictac, & ! OUT
         real(hco_core%xlat1), real(hco_core%xlon1),   & ! IN
         real(hco_core%xlat2), real(hco_core%xlon2)  )   ! IN

    ig1      =  ig1_tictac
    ig2      =  ig2_tictac
    ig3      =  ig3_tictac
    ig4      =  ig4_tictac

    ier = utl_fstecr(hco_core%lon*MPC_DEGREES_PER_RADIAN_R8, npak, &
         iun, dateo, deet, npas, hco_core%ni, 1, 1, ip1,    &
         ip2, ip3, typvar, '>>', etiket, grtyp, ig1,          &
         ig2, ig3, ig4, datyp, .true.)

    ier = utl_fstecr(hco_core%lat*MPC_DEGREES_PER_RADIAN_R8, npak, &
         iun, dateo, deet, npas, 1, hco_core%nj, 1, ip1,    &
         ip2, ip3, typvar, '^^', etiket, grtyp, ig1,          &
         ig2, ig3, ig4, datyp, .true.)

    !- 2.2 2D Field
    allocate(Field2d(hco_core%ni,hco_core%nj))
    Field2d(:,:) = 10.d0

    deet      =  0
    ip1       =  0
    ip2       =  0
    ip3       =  0
    npas      =  0
    datyp     =  1
    grtyp     =  hco_core%grtyp
    typvar    = 'A'
    etiket    = 'COREGRID'
    dateo     =  0
    ig1       =  hco_core%ig1
    ig2       =  hco_core%ig2
    ig3       =  hco_core%ig3
    ig4       =  hco_core%ig4

    ier = utl_fstecr(Field2d, npak,                                    &
         iun, dateo, deet, npas, hco_core%ni, hco_core%nj, 1, ip1, &
         ip2, ip3, typvar, 'P0', etiket, grtyp, ig1,                &
         ig2, ig3, ig4, datyp, .true.)

    deallocate(Field2d)

    !
    !- 3.  Create and Write the extended grid (Analysis) template
    !
    ni_ext = hco_core%ni + grd_ext_x
    nj_ext = hco_core%nj + grd_ext_y

    !- 3.1 Tic-Tac
    allocate(lon_ext(ni_ext))
    allocate(lat_ext(nj_ext))

    !- Copy core grid info
    lon_ext(1:hco_core%ni) = hco_core%lon(:) 
    lat_ext(1:hco_core%nj) = hco_core%lat(:)

    !- Extend the lat lon
    dlon = hco_core%lon(2) - hco_core%lon(1) 
    do i = hco_core%ni + 1, ni_ext
      lon_ext(i) = lon_ext(hco_core%ni) + (i - hco_core%ni) * dlon
    end do

    dlat = hco_core%lat(2) - hco_core%lat(1) 
    do j = hco_core%nj + 1, nj_ext
      lat_ext(j) = lat_ext(hco_core%nj) + (j - hco_core%nj) * dlat
    end do

    !- Write
    deet     =  0
    ip1      =  hco_core%ig1 + 100 ! Must be different from the core grid
    ip2      =  hco_core%ig2 + 100 ! Must be different from the core grid
    if (hco_core%ig3 > 0) then
      ip3      =  hco_core%ig3 + 100 ! Must be different from the core grid
    else
      ip3      =  0
    end if
    npas     =  0
    datyp    =  1
    grtyp    =  hco_core%grtypTicTac
    typvar   = 'X'
    etiket   = 'ANALYSIS'
    dateo    =  0
    ig1      =  ig1_tictac
    ig2      =  ig2_tictac
    ig3      =  ig3_tictac
    ig4      =  ig4_tictac

    ier = utl_fstecr(lon_ext*MPC_DEGREES_PER_RADIAN_R8, npak, &
         iun, dateo, deet, npas, ni_ext, 1, 1, ip1,  &
         ip2, ip3, typvar, '>>', etiket, grtyp, ig1,    &
         ig2, ig3, ig4, datyp, .true.)

    ier = utl_fstecr(lat_ext*MPC_DEGREES_PER_RADIAN_R8, npak, &
         iun, dateo, deet, npas, 1, nj_ext, 1, ip1,  &
         ip2, ip3, typvar, '^^', etiket, grtyp, ig1,    &
         ig2, ig3, ig4, datyp, .true.)

    deallocate(lon_ext)
    deallocate(lat_ext)

    !- 3.2 2D Field
    allocate(Field2d(ni_ext,nj_ext))
    Field2d(:,:) = 10.d0

    deet      =  0
    ip1       =  0
    ip2       =  0
    ip3       =  0
    npas      =  0
    datyp     =  1
    grtyp     =  hco_core%grtyp
    typvar    = 'A'
    etiket    = 'ANALYSIS'
    dateo     =  0
    ig1       =  hco_core%ig1 + 100 ! Must be different from the core grid
    ig2       =  hco_core%ig2 + 100 ! Must be different from the core grid
    if (hco_core%ig3 > 0) then
      ig3      =  hco_core%ig3 + 100 ! Must be different from the core grid
    else
      ig3      =  0
    end if
    ig4       =  0

    ier = utl_fstecr(Field2d, npak,                                  &
         iun, dateo, deet, npas, ni_ext, nj_ext, 1, ip1, &
         ip2, ip3, typvar, 'P0', etiket, grtyp, ig1,     &
         ig2, ig3, ig4, datyp, .true.)

    deallocate(Field2d)

    !
    !- 4. Write the vertical grid description
    !

    if (vco%vgridPresent) then
      !- 4.1 Write the toc-toc
      status = vgd_write(vco%vgrid,iun,'fst')
      
      if ( status /= VGD_OK ) then
        call utl_abort('createLamTemplateGrids: ERROR with vgd_write')
      end if
      
      !- 4.2 Write a dummy 2D field for each MM and TH levels
      npak   = -12
      dateo  = 0
      deet   = 0
      npas   = 0
      ni     = 4
      nj     = 2
      nk     = 1
      ip2    = 0
      ip3    = 0
      typvar = 'A'
      etiket = 'VERTICALGRID'
      grtyp  = 'G'
      ig1    = 0
      ig2    = 0
      ig3    = 0
      ig4    = 0
      datyp  = 1
      
      allocate(dummy2D(ni,nj))
      dummy2D(:,:) = 0.0
      
      do lev = 1, vco%nlev_M
        ip1 = vco%ip1_M(lev)
        ier = fstecr(dummy2D, work, npak, iun, dateo, deet, npas, ni, nj, &
                     nk, ip1, ip2, ip3, typvar, 'MM', etiket, grtyp,              &
                     ig1, ig2, ig3, ig4, datyp, .true.)
      end do
      do lev = 1, vco%nlev_T
        ip1 = vco%ip1_T(lev)
        ier = fstecr(dummy2D, work, npak, iun, dateo, deet, npas, ni, nj, &
                     nk, ip1, ip2, ip3, typvar, 'TH', etiket, grtyp,              &
                     ig1, ig2, ig3, ig4, datyp, .true.)
      end do

      deallocate(dummy2D)
    end if ! vco%vgridPresent

    !
    !- 5.  Closing the output template file
    !
    ier = fstfrm(iun)
    ier = fclos (iun)

  end subroutine agd_createLamTemplateGrids

end module analysisGrid_mod
