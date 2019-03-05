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
!! MODULE HorizontalCoord (prefix="hco" category='7. Low-level data objects and utilities')
!!
!! *Purpose*: Derived type and procedures related to the horizontal grid 
!!            coordinate for various grids (global and limited area).
!!
!--------------------------------------------------------------------------
module HorizontalCoord_mod
  use mpi_mod
  use mpivar_mod
  use mathPhysConstants_mod
  use utilities_mod
  implicit none
  save
  private

  ! Public derived type
  public :: struct_hco

  ! Public Subroutines
  public :: hco_SetupFromFile, hco_equal, hco_deallocate, hco_mpiBcast

  integer, parameter :: maxNumSubGrid = 2

  type :: struct_hco
     character(len=32)    :: gridname
     logical              :: initialized = .false.
     integer              :: ni
     integer              :: nj
     character(len=1)     :: grtyp
     integer              :: ig1
     integer              :: ig2
     integer              :: ig3
     integer              :: ig4
     integer              :: EZscintID
     integer              :: numSubGrid
     integer              :: EZscintIDsubGrids(maxNumSubGrid)
     real(8), allocatable :: lat(:) ! in radians
     real(8), allocatable :: lon(:) ! in radians
     real(4), allocatable :: lat2d_4(:,:) ! in radians
     real(4), allocatable :: lon2d_4(:,:) ! in radians
     real(8)              :: dlat   ! in radians
     real(8)              :: dlon   ! in radians
     logical              :: global
     logical              :: rotated
     real(8)              :: xlat1, xlat1_yan
     real(8)              :: xlon1, xlon1_yan
     real(8)              :: xlat2, xlat2_yan
     real(8)              :: xlon2, xlon2_yan
     real(4), allocatable :: tictacU(:)
     integer, allocatable :: mask(:,:)
  end type struct_hco

  contains

!--------------------------------------------------------------------------
! hco_SetupFromFile
!--------------------------------------------------------------------------
  subroutine hco_SetupFromFile(hco, TemplateFile, EtiketName, GridName_opt, &
       varName_opt)
    implicit none

    type(struct_hco), pointer    :: hco
    character(len=*), intent(in) :: TemplateFile
    character(len=*), intent(in) :: EtiketName
    character(len=*), intent(in), optional :: GridName_opt
    character(len=*), intent(in), optional :: varName_opt

    real(8), allocatable :: lat_8(:)
    real(8), allocatable :: lon_8(:)

    real    :: xlat1_4, xlon1_4, xlat2_4, xlon2_4
    real    :: xlat1_yan_4, xlon1_yan_4, xlat2_yan_4, xlon2_yan_4

    integer :: iu_template, numSubGrid
    integer :: fnom, fstlir, fstouv, fstfrm, fclos, fstluk
    integer :: ezqkdef, ezget_nsubgrids, ezget_subgridids, ezgprm
    integer :: key, fstinf, fstprm, ier, fstinl, EZscintID, EZscintIDsubGrids(maxNumSubGrid)
    integer :: ni, nj, ni_tictacU, ni_t, nj_t, nlev_t, i, j, gdll
    integer :: dateo, deet, npas, nk, nbits, datyp
    integer :: ip1, ip2, ip3, swa, lng, dltf, ubc
    integer :: extra1, extra2, extra3
    integer :: ig1, ig2, ig3, ig4
    integer :: ig1_tictac, ig2_tictac, ig3_tictac, ig4_tictac
    integer :: ni_yy, nj_yy,  ig1_yy, ig2_yy, ig3_yy, ig4_yy

    logical :: FileExist, global, rotated

    character(len=4 ) :: nomvar
    character(len=2 ) :: typvar
    character(len=1 ) :: grtyp, grtyp_tictac, grtyp_yy
    character(len=12) :: etiket

    if( .not.associated(hco) ) then
      allocate(hco)
    else
      call utl_abort('hco_setupFromFile: supplied hco must be null')
    endif

    !
    !- 1.  Open/Check template file
    !
    inquire(file=trim(TemplateFile), exist=FileExist)

    if ( FileExist ) then
      iu_template = 0
      ier = fnom(iu_template,trim(TemplateFile),'RND+OLD+R/O',0)
      if ( ier == 0 ) then
        write(*,*)
        write(*,*) 'hco_setupFromFile: Template File =', trim(TemplateFile)
        ier = fstouv(iu_template,'RND+OLD')
      else
        write(*,*)
        write(*,*) 'hco_SetupFromFile: Error in opening the template grid file'
        write(*,*) trim(TemplateFile)
        call utl_abort('hco_SetupFromFile')
      end if
    else
      write(*,*)
      write(*,*) 'hco_SetupFromFile: template grid file DOES NOT EXIST'
      write(*,*) trim(TemplateFile)
      call utl_abort('hco_SetupFromFile')
    end if

    !
    !- 2.  Get Horizontal grid info
    !

    !- 2.1 Grid size and grid projection info
    dateo  = -1
    etiket = EtiketName
    ip1    = -1
    ip2    = -1
    ip3    = -1
    typvar = ' '
    if (present(varName_opt)) then
      nomvar = varName_opt
    else
      nomvar = 'P0'
    end if

    key = fstinf( iu_template,                                & ! IN
                  ni, nj, nk,                                 & ! OUT
                  dateo, etiket, ip1, ip2, ip3, typvar, nomvar )! IN

    if (key < 0) then
      write(*,*)
      write(*,*) 'hco_SetupFromFile: Unable to find output horiz grid info using = ',nomvar
      write(*,*) 'hco_SetupFromFile: with etiket = ',trim(EtiketName)
      call utl_abort('hco_setupFromFile')
    end if

    ier = fstprm( key,                                             & ! IN
                  dateo, deet, npas, ni, nj, nk, nbits,            & ! OUT
                  datyp, ip1, ip2, ip3, typvar, nomvar, etiket,    & ! OUT
                  grtyp, ig1, ig2, ig3,                            & ! OUT
                  ig4, swa, lng, dltf, ubc, extra1, extra2, extra3 ) ! OUT

    if ( trim(grtyp) == 'G' .and. ig2 == 1 ) then
       call utl_abort('hco_setupFromFile: ERROR: due to bug in ezsint, Gaussian grid with ig2=1 no longer supported')
    end if

    EZscintID  = ezqkdef( ni, nj, grtyp, ig1, ig2, ig3, ig4, iu_template )   ! IN
    numSubGrid = 1
    EZscintIDsubGrids(:) = -999

    allocate(lat_8(1:nj))
    allocate(lon_8(1:ni))

    allocate(hco%lat2d_4(1:ni,1:nj))
    allocate(hco%lon2d_4(1:ni,1:nj))

    ier = gdll( EZscintID,       & ! IN
                hco%lat2d_4, hco%lon2d_4 ) ! OUT

    xlat1_yan_4 = -999.9
    xlon1_yan_4 = -999.9
    xlat2_yan_4 = -999.9
    xlon2_yan_4 = -999.9

    if (mpi_myid == 0) write(*,*) 'hco_setupFromFile: grtyp, ni, nj = ', grtyp, ni, nj

    !- 2.2 Rotated lat-lon grid
    if ( trim(grtyp) == 'Z' ) then

       !-  2.2.1 Read the Longitudes
       dateo  = -1
       etiket = EtiketName
       ip1    = ig1
       ip2    = ig2
       ip3    = -1
       typvar = 'X'
       nomvar = '>>'

       ier = utl_fstlir( lon_8,                                      & ! OUT 
                         iu_template,                                & ! IN
                         ni_t, nj_t, nlev_t,                         & ! OUT
                         dateo, etiket, ip1, ip2, ip3, typvar,nomvar)  ! IN

       if (ier < 0) then
         write(*,*)
         write(*,*) 'hco_SetupFromFile: Unable to find >> grid descriptors'
         call utl_abort('hco_setupFromFile')
       end if

       !  Test if the dimensions are compatible with the grid
       if ( ni_t /= ni .or. nj_t /= 1 ) then
         write(*,*)
         write(*,*) 'hco_SetupFromFile: Incompatible >> grid descriptors !'
         write(*,*) 'Found     :', ni_t, nj_t
         write(*,*) 'Should be :', ni, 1
         call utl_abort('hco_setupFromFile')
       end if

       !-  2.2.2 Read the latitudes
       dateo  = -1
       etiket = EtiketName
       ip1    = ig1
       ip2    = ig2
       ip3    = -1
       typvar = 'X'
       nomvar = '^^'

       ier = utl_fstlir( lat_8,                                      & ! OUT 
                         iu_template,                                & ! IN
                         ni_t, nj_t, nlev_t,                         & ! OUT
                         dateo, etiket, ip1, ip2, ip3, typvar,nomvar)  ! IN

       if (ier < 0) then
         write(*,*)
         write(*,*) 'hco_SetupFromFile: Unable to find ^^ grid descriptors'
         call utl_abort('hco_setupFromFile')
       end if

       !  Test if the dimensions are compatible with the grid
       if ( ni_t /= 1 .or. nj_t /= nj ) then
         write(*,*)
         write(*,*) 'hco_SetupFromFile: Incompatible ^^ grid descriptors !'
         write(*,*) 'Found     :', ni_t, nj_t
         write(*,*) 'Should be :', 1, nj
         call utl_abort('hco_setupFromFile')
       end if

       !- 2.2.3 Do we have a rotated grid ?
       dateo  = -1
       etiket = EtiketName
       ip1    = ig1
       ip2    = ig2
       ip3    = -1
       typvar = 'X'
       nomvar = '^^'

       key = fstinf( iu_template,                                   & ! IN
                     ni_t, nj_t, nk,                                & ! OUT
                     dateo, etiket, ip1, ip2, ip3, typvar, nomvar )   ! IN

       ier = fstprm( key,                                           & ! IN
                     dateo, deet, npas, ni_t, nj_t, nk, nbits,      & ! OUT
                     datyp, ip1, ip2, ip3, typvar, nomvar, etiket,  & ! OUT
                     grtyp_tictac, ig1_tictac, ig2_tictac,          & ! OUT
                     ig3_tictac, ig4_tictac, swa, lng, dltf,        & ! OUT
                     ubc, extra1, extra2, extra3 )                    ! OUT

       call cigaxg ( grtyp_tictac,                                  & ! IN
                     xlat1_4, xlon1_4, xlat2_4, xlon2_4,            & ! OUT
                     ig1_tictac, ig2_tictac, ig3_tictac, ig4_tictac ) ! IN

       if ( xlat1_4 == 0.0 .and. xlat2_4 == 0.0 ) then
          rotated = .false.
       else
          rotated = .true.
       end if
       if (mpi_myid == 0) then
         write(*,*) 'hco_setupFromFile: xlat1/2, xlon1/2, rotated = ',  &
              xlat1_4, xlat2_4, xlon1_4, xlon2_4, rotated
       end if

       !- 2.2.4 Is this a global or a LAM domain ?
       call global_or_lam( global,     & ! OUT
                           lon_8, ni )   ! IN

    !- 2.3 Gaussian Grid
    elseif ( trim(grtyp) == 'G' ) then

      !-  2.3.1 Find the latitudes and longitudes
      lon_8(:) = real(hco%lon2d_4(:,nj/2),8)
      lat_8(:) = real(hco%lat2d_4(1,:),8)

      !- 2.3.2 This grid type is not rotated
      rotated = .false.
      xlat1_4 =   0.0
      xlon1_4 = 180.0
      xlat2_4 =   0.0
      xlon2_4 = 180.0

      !- 2.3.3 We know this is a global grid
      global = .true.

    !- 2.4 Universal Grid (Yin-Yang) - not fully supported: use at own risk!
    elseif ( trim(grtyp) == 'U' ) then

      !-  2.4.1 Read the tic-tac vector
      dateo  = -1
      etiket = ' '
      ip1    = ig1
      ip2    = ig2
      ip3    = ig3
      typvar = 'X'
      nomvar = '^>'

      ni_tictacU = 5 + 2 * (10 + ni + nj/2)
      allocate(hco%tictacU(ni_tictacU))
      ier = fstlir( hco%tictacU,                                 & ! OUT 
                    iu_template,                                 & ! IN
                    ni_t, nj_t, nlev_t,                          & ! OUT
                    dateo, etiket, ip1, ip2, ip3, typvar, nomvar)  ! IN

      if (ier < 0) then
        write(*,*)
        write(*,*) 'hco_SetupFromFile: Unable to find ^> grid descriptors'
        call utl_abort('hco_setupFromFile')
      end if

      !  Test if the dimensions are compatible with the grid
      if ( ni_t /= ni_tictacU .or. nj_t /= 1 ) then
        write(*,*)
        write(*,*) 'hco_SetupFromFile: Incompatible ^> grid descriptors !'
        write(*,*) 'Found     :', ni_t, nj_t
        write(*,*) 'Should be :', ni_tictacU, 1
        call utl_abort('hco_setupFromFile')
      end if

      !-  2.4.1 Initialize latitudes and longitudes to dummy values - should not be used!
      lon_8(:) = -999.999d0
      lat_8(:) = -999.999d0

      !-  2.4.2 Yin-Yan subgrid IDs
      numSubGrid = ezget_nsubgrids(EZscintID)
      if ( numSubGrid /= 2 ) then
        call utl_abort('hco_setupFromFile: CONFUSED! number of sub grids must be 2')
      end if
      ier = ezget_subgridids(EZscintID, EZscintIDsubGrids)

      !-  2.4.3 Determine parameters related to Yin and Yan grid rotations
      rotated = .true.  ! since Yin-Yan is made up of 2 grids with different rotations

      ier = ezgprm( EZscintIDsubGrids(1), grtyp_yy, ni_yy, nj_yy, ig1_yy, ig2_yy, ig3_yy, ig4_yy )
      grtyp_yy = 'E' ! needed since ezgprm returns 'Z', but grtyp for tictac should be 'E'
      call cigaxg ( grtyp_yy,                           & ! IN
                    xlat1_4, xlon1_4, xlat2_4, xlon2_4, & ! OUT
                    ig1_yy, ig2_yy, ig3_yy, ig4_yy )      ! IN

      ier = ezgprm( EZscintIDsubGrids(2), grtyp_yy, ni_yy, nj_yy, ig1_yy, ig2_yy, ig3_yy, ig4_yy )
      grtyp_yy = 'E' ! needed since ezgprm returns 'Z', but grtyp for tictac should be 'E'
      call cigaxg ( grtyp_yy,                                           & ! IN
                    xlat1_yan_4, xlon1_yan_4, xlat2_yan_4, xlon2_yan_4, & ! OUT
                    ig1_yy, ig2_yy, ig3_yy, ig4_yy )                      ! IN

      rotated = .true.  ! since Yin-Yan is made up of 2 grids with different rotations

      !-  2.4.3 We know this is a global grid
      global = .true.

    else
      write(*,*)
      write(*,*) 'hco_SetupFromFile: Only grtyp = Z or G or U are supported !, grtyp = ', trim(grtyp)
      call utl_abort('hco_setupFromFile')
    end if

    !
    !- 3.  Initialized Horizontal Grid Structure
    !
    allocate(hco % lat(1:nj))
    allocate(hco % lon(1:ni))

    if ( present(gridName_opt) ) then
      hco % gridname     = trim(gridName_opt)
    else
      hco % gridname     = 'UNDEFINED'
    end if
    hco % ni                   = ni
    hco % nj                   = nj
    hco % grtyp                = trim(grtyp) 
    hco % ig1                  = ig1
    hco % ig2                  = ig2
    hco % ig3                  = ig3
    hco % ig4                  = ig4
    hco % EZscintID            = EZscintID
    hco % numSubGrid           = numSubGrid
    hco % EZscintIDsubGrids(:) = EZscintIDsubGrids(:)
    hco % lon(:)               = lon_8(:) * MPC_RADIANS_PER_DEGREE_R8
    hco % lat(:)               = lat_8(:) * MPC_RADIANS_PER_DEGREE_R8
    hco % dlon                 = (lon_8(2) - lon_8(1)) * MPC_RADIANS_PER_DEGREE_R8
    hco % dlat                 = (lat_8(2) - lat_8(1)) * MPC_RADIANS_PER_DEGREE_R8
    hco % global               = global
    hco % rotated              = rotated
    hco % xlat1                = real(xlat1_4,8)
    hco % xlon1                = real(xlon1_4,8)
    hco % xlat2                = real(xlat2_4,8)
    hco % xlon2                = real(xlon2_4,8)
    hco % xlat1_yan            = real(xlat1_yan_4,8)
    hco % xlon1_yan            = real(xlon1_yan_4,8)
    hco % xlat2_yan            = real(xlat2_yan_4,8)
    hco % xlon2_yan            = real(xlon2_yan_4,8)
    hco % initialized          = .true.

    hco%lat2d_4(:,:) = hco%lat2d_4(:,:) * MPC_RADIANS_PER_DEGREE_R8
    hco%lon2d_4(:,:) = hco%lon2d_4(:,:) * MPC_RADIANS_PER_DEGREE_R8

    deallocate(lat_8)
    deallocate(lon_8)

    !- 3.1  Read the mask if present

    dateo  = -1
    etiket = EtiketName
    ip1    = -1
    ip2    = -1
    ip3    = -1
    typvar = '@@'
    nomvar = ' '

    key = fstinf( iu_template,                                   & ! IN
                  ni_t, nj_t, nk,                                & ! OUT
                  dateo, etiket, ip1, ip2, ip3, typvar, nomvar )   ! IN

    if (key >= 0) then

      !  Test if the dimensions are compatible with the grid
      if ( ni_t /= ni .or. nj_t /= nj ) then
        write(*,*)
        write(*,*) 'hco_SetupFromFile: Incompatible mask grid descriptors !'
        write(*,*) 'Found     :', ni_t, nj_t
        write(*,*) 'Should be :', ni, nj
        call utl_abort('hco_setupFromFile')
      end if

      allocate(hco % mask(ni,nj))

      key = fstluk(hco%mask, key, ni_t, nj_t, nk)

    end if

    !
    !- 4.  Close the input file
    !
    ier = fstfrm(iu_template)
    ier = fclos (iu_template)

  end subroutine hco_SetupFromFile

!--------------------------------------------------------------------------
! Global_or_lam
!--------------------------------------------------------------------------
  subroutine global_or_lam(global, lon, ni)
    implicit none

    integer, intent(in)  :: ni
    real(8), intent(in)  :: lon(ni)
    logical, intent(out) :: global

    real(8) :: dx, next_lon

    dx       = lon(2) - lon(1)
    next_lon = lon(ni) + 1.5d0 * dx

    write(*,*)
    write(*,*) 'dx       = ',dx
    write(*,*) 'lon(ni)  = ',lon(ni)
    write(*,*) 'next_lon = ',next_lon
    write(*,*) 'lon(1)   = ',lon(1)

    if ( next_lon - lon(1) > 360.0d0 ) then

      global = .true.
      if ( lon(1) == lon(ni) ) then
        write(*,*)
        write(*,*) ' *** Global Grid where i = ni (repetition) '
      else  
        write(*,*)
        write(*,*) ' *** Global Grid where i /= ni '
      end if

    else

      global = .false.
      write(*,*)
      write(*,*) ' *** Limited-Area Grid '

    end if

  end subroutine global_or_lam

!--------------------------------------------------------------------------
! mpiBcast
!--------------------------------------------------------------------------
  subroutine hco_mpiBcast(hco)
    implicit none
    type(struct_hco), pointer :: hco
    integer :: ierr
    integer, external :: ezqkdef
    logical           :: maskAllocated

    write(*,*) 'hco_mpiBcast: starting'

    if ( mpi_myid > 0 ) then
      if( .not.associated(hco) ) then
        allocate(hco)
      else
        call utl_abort('hco_mpiBcast: hco must be nullified for mpi task id > 0')
      endif
    endif

    call rpn_comm_bcastc(hco%gridname, len(hco%gridname), 'MPI_CHARACTER', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%initialized, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%ni, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%nj, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcastc(hco%grtyp, len(hco%grtyp), 'MPI_CHARACTER', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%ig1, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%ig2, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%ig3, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%ig4, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    if ( mpi_myid > 0 ) then
      allocate(hco % lat(1:hco%nj))
      allocate(hco % lon(1:hco%ni))
    endif
    call rpn_comm_bcast(hco%lat, size(hco%lat), 'MPI_REAL8', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%lon, size(hco%lon), 'MPI_REAL8', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%dlat, 1, 'MPI_REAL8', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%dlon, 1, 'MPI_REAL8', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%global, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%rotated, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%xlat1, 1, 'MPI_REAL8', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%xlon1, 1, 'MPI_REAL8', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%xlat2, 1, 'MPI_REAL8', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%xlon2, 1, 'MPI_REAL8', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%xlat1_yan, 1, 'MPI_REAL8', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%xlon1_yan, 1, 'MPI_REAL8', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%xlat2_yan, 1, 'MPI_REAL8', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%xlon2_yan, 1, 'MPI_REAL8', 0, 'GRID', ierr)
    if ( hco%grtyp == 'U' ) then
      if ( mpi_myid > 0 ) then
        allocate(hco%tictacU(5 + 2 * (10 + hco%ni + hco%nj/2)))
      endif
      call rpn_comm_bcast(hco%tictacU, size(hco%tictacU), 'MPI_REAL4', 0, 'GRID', ierr)
    end if

    if ( mpi_myid > 0 ) then
      if ( hco%grtyp == 'G' ) then
        hco%EZscintID  = ezqkdef( hco%ni, hco%nj, hco%grtyp, hco%ig1, hco%ig2, hco%ig3, hco%ig4, 0 )
      else
        write(*,*) 'hco_mpiBcast: Warning! Grid ID for EZSCINT not set for grtyp = ',hco%grtyp
        hco%EZscintID  = -1
      endif
    endif

    maskAllocated = allocated(hco%mask)
    call rpn_comm_bcast(maskAllocated, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)
    if ( maskAllocated ) then
       if ( mpi_myid > 0 ) then
          allocate(hco % mask(hco%ni,hco%nj))
       endif
       call rpn_comm_bcast(hco%mask, size(hco%mask), 'MPI_INTEGER', 0, 'GRID', ierr)
    endif

    write(*,*) 'hco_mpiBcast: done'

  end subroutine hco_mpiBcast

!--------------------------------------------------------------------------
! Equal ?
!--------------------------------------------------------------------------
  function hco_equal(hco1, hco2) result(equal)
    implicit none
    type(struct_hco), pointer :: hco1, hco2
    logical                   :: equal

    equal = .true.
    equal = equal .and. (hco1%ni == hco2%ni)
    equal = equal .and. (hco1%nj == hco2%nj)
    if (.not. equal) then
      write(*,*) 'hco_equal: dimensions not equal ', hco1%ni, hco1%nj, hco2%ni, hco2%nj
      return
    endif

    equal = equal .and. (hco1%grtyp   ==    hco2%grtyp)
    if (.not. equal) then
      write(*,*) 'hco_equal: grid type not equal'
      return
    endif

    equal = equal .and. (hco1%dlat == hco2%dlat)
    equal = equal .and. (hco1%dlon == hco2%dlon)
    if (.not. equal) then
      write(*,*) 'hco_equal: grid spacing not equal'
      return
    endif

    if(hco1%grtyp .eq. 'G') then
      equal = equal .and. (hco1%ig2 == hco2%ig2)
      if (.not. equal) then
        write(*,*) 'hco_equal: Gaussian grid ig2 not equal'
        return
      endif
    endif    

    equal = equal .and. (hco1%rotated .eqv. hco2%rotated)
    equal = equal .and. (hco1%xlat1   ==    hco2%xlat1)
    equal = equal .and. (hco1%xlon1   ==    hco2%xlon1)
    equal = equal .and. (hco1%xlat2   ==    hco2%xlat2)
    equal = equal .and. (hco1%xlon2   ==    hco2%xlon2)
    equal = equal .and. (hco1%xlat1_yan ==  hco2%xlat1_yan)
    equal = equal .and. (hco1%xlon1_yan ==  hco2%xlon1_yan)
    equal = equal .and. (hco1%xlat2_yan ==  hco2%xlat2_yan)
    equal = equal .and. (hco1%xlon2_yan ==  hco2%xlon2_yan)
    if (.not. equal) then
      write(*,*) 'hco_equal: rotation not equal'
      return
    endif

    equal = equal .and. all(hco1%lat(:) == hco2%lat(:))
    equal = equal .and. all(hco1%lon(:) == hco2%lon(:))
    if (.not. equal) then
      write(*,*) 'hco_equal: lat/lon not equal'
      return
    endif

  end function hco_equal


  subroutine hco_deallocate( hco )
    implicit none
    type(struct_hco), pointer :: hco

    deallocate(hco % lat)
    deallocate(hco % lon)
    deallocate(hco % lat2d_4)
    deallocate(hco % lon2d_4)
    deallocate(hco)
    nullify(hco)

  end subroutine hco_deallocate

end module HorizontalCoord_mod
