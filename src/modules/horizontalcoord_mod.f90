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

module HorizontalCoord_mod
  ! MODULE HorizontalCoord_mod (prefix='hco' category='7. Low-level data objects and utilities')
  !
  ! :Purpose: Derived type and procedures related to the horizontal grid
  !           coordinate for various grids (global and limited area).
  
  use mpi_mod
  use mpivar_mod
  use earthConstants_mod
  use mathPhysConstants_mod
  use utilities_mod
  use varNameList_mod
  use physicsFunctions_mod
  implicit none
  save
  private

  ! Public derived type
  public :: struct_hco

  ! Public Subroutines
  public :: hco_SetupFromFile, hco_equal, hco_deallocate, hco_mpiBcast, hco_weight, hco_setupYgrid

  integer, parameter :: maxNumSubGrid = 2

  type :: struct_hco
    character(len=32)    :: gridname
    logical              :: initialized = .false.
    integer              :: ni
    integer              :: nj
    character(len=1)     :: grtyp
    character(len=1)     :: grtypTicTac
    integer              :: ig1
    integer              :: ig2
    integer              :: ig3
    integer              :: ig4
    integer              :: EZscintID = -1
    integer              :: numSubGrid
    integer              :: EZscintIDsubGrids(maxNumSubGrid)
    real(8), allocatable :: lat(:) ! in radians
    real(8), allocatable :: lon(:) ! in radians
    real(4), allocatable :: lat2d_4(:,:) ! in radians
    real(4), allocatable :: lon2d_4(:,:) ! in radians
    real(8)              :: dlat   ! in radians
    real(8)              :: dlon   ! in radians
    real(8)              :: maxGridSpacing ! in meter
    logical              :: global
    logical              :: rotated
    real(8)              :: xlat1, xlat1_yan
    real(8)              :: xlon1, xlon1_yan
    real(8)              :: xlat2, xlat2_yan
    real(8)              :: xlon2, xlon2_yan
    real(4), allocatable :: tictacU(:)
  end type struct_hco

contains

  !--------------------------------------------------------------------------
  ! hco_SetupFromFile
  !--------------------------------------------------------------------------
  subroutine hco_SetupFromFile(hco, TemplateFile, EtiketName, GridName_opt, &
                               varName_opt)
    !
    ! :Purpose: to initialize hco structure from a template file
    !           
    implicit none
    ! arguments:
    type(struct_hco), pointer    :: hco
    character(len=*), intent(in) :: TemplateFile
    character(len=*), intent(in) :: EtiketName
    character(len=*), intent(in), optional :: GridName_opt
    character(len=*), intent(in), optional :: varName_opt
    ! locals:
    real(8), allocatable :: lat_8(:)
    real(8), allocatable :: lon_8(:)

    real(8) :: maxDeltaLat, deltaLon, maxDeltaLon, maxGridSpacing 
    real(8), save :: maxGridSpacingPrevious = -1.0d0
    real(4) :: xlat1_4, xlon1_4, xlat2_4, xlon2_4
    real(4) :: xlat1_yan_4, xlon1_yan_4, xlat2_yan_4, xlon2_yan_4
    
    integer :: iu_template, numSubGrid, varIndex
    integer :: fnom, fstlir, fstouv, fstfrm, fclos
    integer :: ezqkdef, ezget_nsubgrids, ezget_subgridids, ezgprm
    integer :: key, fstinf, fstprm, ier, EZscintID, EZscintIDsubGrids(maxNumSubGrid)
    integer :: ni, nj, ni_tictacU, ni_t, nj_t, nlev_t, gdll
    integer :: dateo, deet, npas, nk, nbits, datyp
    integer :: ip1, ip2, ip3, swa, lng, dltf, ubc
    integer :: extra1, extra2, extra3
    integer :: ig1, ig2, ig3, ig4
    integer :: ig1_tictac, ig2_tictac, ig3_tictac, ig4_tictac
    integer :: ni_yy, nj_yy,  ig1_yy, ig2_yy, ig3_yy, ig4_yy
    integer :: latIndex, lonIndex, latIndexBeg, latIndexEnd  

    logical :: FileExist, global, rotated, foundVarNameInFile

    character(len=4 ) :: nomvar
    character(len=2 ) :: typvar
    character(len=1 ) :: grtyp, grtypTicTac
    character(len=12) :: etiket
    
    if( .not.associated(hco) ) then
      allocate(hco)
    else
      call utl_abort('hco_setupFromFile: supplied hco must be null')
    end if

    !
    !- 1.1  Determine which variable to use for defining the grid
    !
    if (present(varName_opt)) then
      
      ! User specified variable name
      nomvar = varName_opt
      
    else

      ! First try to use P0
      nomvar = 'P0'
      
      if ( .not. utl_varNamePresentInFile(nomvar,fileName_opt=trim(TemplateFile)) ) then
        ! P0 not present, look for another suitable variable in the file
        
        foundVarNameInFile = .false.
        do varIndex = 1, vnl_numvarmax
          nomvar = vnl_varNameList(varIndex)

          ! check if variable is in the file
          if ( .not. utl_varNamePresentInFile(nomvar,fileName_opt=trim(TemplateFile)) ) cycle

          foundVarNameInFile = .true.
          exit
          
        end do

        if ( .not. foundVarNameInFile) call utl_abort('hco_SetupFromFile: NO variables found in the file!!!')

      end if

    end if

    write(*,*) 'hco_SetupFromFile: defining hco by varname= ', nomvar


    !
    !- 1.2  Open/Check template file
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
    
    key = fstinf( iu_template,                                & ! IN
                  ni, nj, nk,                                 & ! OUT
                  dateo, etiket, ip1, ip2, ip3, typvar, nomvar )! IN

    if (key < 0) then
      write(*,*)
      write(*,*) 'hco_SetupFromFile: Unable to find output horiz grid info using = ',nomvar
      write(*,*) '                   with etiket = ',trim(EtiketName)
      call utl_abort('hco_setupFromFile: unable to setup the structure')
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
    
    ier = gdll( EZscintID,               & ! IN
                hco%lat2d_4, hco%lon2d_4 ) ! OUT
    
    xlat1_yan_4 = -999.9
    xlon1_yan_4 = -999.9
    xlat2_yan_4 = -999.9
    xlon2_yan_4 = -999.9

    grtypTicTac = 'X'

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
                    grtypTicTac, ig1_tictac, ig2_tictac,           & ! OUT
                    ig3_tictac, ig4_tictac, swa, lng, dltf,        & ! OUT
                    ubc, extra1, extra2, extra3 )                    ! OUT

      call cigaxg ( grtypTicTac,                                   & ! IN
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
    else if ( trim(grtyp) == 'G' ) then
    
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
    else if ( trim(grtyp) == 'U' ) then
      
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
      
      ier = ezgprm( EZscintIDsubGrids(1), grtypTicTac, ni_yy, nj_yy, ig1_yy, ig2_yy, ig3_yy, ig4_yy )
      grtypTicTac = 'E' ! needed since ezgprm returns 'Z', but grtyp for tictac should be 'E'
      call cigaxg ( grtypTicTac,                        & ! IN
                    xlat1_4, xlon1_4, xlat2_4, xlon2_4, & ! OUT
                    ig1_yy, ig2_yy, ig3_yy, ig4_yy )      ! IN
      
      ier = ezgprm( EZscintIDsubGrids(2), grtypTicTac, ni_yy, nj_yy, ig1_yy, ig2_yy, ig3_yy, ig4_yy )
      grtypTicTac = 'E' ! needed since ezgprm returns 'Z', but grtyp for tictac should be 'E'
      call cigaxg ( grtypTicTac,                                        & ! IN
                    xlat1_yan_4, xlon1_yan_4, xlat2_yan_4, xlon2_yan_4, & ! OUT
                    ig1_yy, ig2_yy, ig3_yy, ig4_yy )                      ! IN
      
      rotated = .true.  ! since Yin-Yan is made up of 2 grids with different rotations
      
      !-  2.4.3 We know this is a global grid
      global = .true.
      
      !- 2.5 Irregular structure
    else if ( trim(grtyp) == 'Y' ) then
      
      !- 2.5.1 This grid type is not rotated
      rotated = .false.
      xlat1_4 = 0.0
      xlon1_4 = 0.0
      xlat2_4 = 1.0
      xlon2_4 = 1.0
      
      grtypTicTac = 'L'
      
      !- 2.5.2 Test using first row of longitudes (should work for ORCA grids)
      lon_8(:) = hco%lon2d_4(:,1)
      call global_or_lam( global,     & ! OUT
                          lon_8, ni )   ! IN
      
      !-  2.5.3 Initialize latitudes and longitudes to dummy values - should not be used!
      lon_8(:) = -999.999d0
      lat_8(:) = -999.999d0
      
    else
      write(*,*)
      write(*,*) 'hco_SetupFromFile: Only grtyp = Z or G or U or Y are supported !, grtyp = ', trim(grtyp)
      call utl_abort('hco_setupFromFile')
    end if
    
    !
    !- 3.  Initialized Horizontal Grid Structure
    !
    allocate(hco%lat(1:nj))
    allocate(hco%lon(1:ni))
    
    if ( present(gridName_opt) ) then
      hco%gridname     = trim(gridName_opt)
    else
      hco%gridname     = 'UNDEFINED'
    end if
    hco%ni                   = ni
    hco%nj                   = nj
    hco%grtyp                = trim(grtyp) 
    hco%grtypTicTac          = trim(grtypTicTac)
    hco%ig1                  = ig1
    hco%ig2                  = ig2
    hco%ig3                  = ig3
    hco%ig4                  = ig4
    hco%EZscintID            = EZscintID
    hco%numSubGrid           = numSubGrid
    hco%EZscintIDsubGrids(:) = EZscintIDsubGrids(:)
    hco%lon(:)               = lon_8(:) * MPC_RADIANS_PER_DEGREE_R8
    hco%lat(:)               = lat_8(:) * MPC_RADIANS_PER_DEGREE_R8
    hco%dlon                 = (lon_8(2) - lon_8(1)) * MPC_RADIANS_PER_DEGREE_R8
    hco%dlat                 = (lat_8(2) - lat_8(1)) * MPC_RADIANS_PER_DEGREE_R8
    hco%global               = global
    hco%rotated              = rotated
    hco%xlat1                = real(xlat1_4,8)
    hco%xlon1                = real(xlon1_4,8)
    hco%xlat2                = real(xlat2_4,8)
    hco%xlon2                = real(xlon2_4,8)
    hco%xlat1_yan            = real(xlat1_yan_4,8)
    hco%xlon1_yan            = real(xlon1_yan_4,8)
    hco%xlat2_yan            = real(xlat2_yan_4,8)
    hco%xlon2_yan            = real(xlon2_yan_4,8)
    hco%initialized          = .true.
  
    hco%lat2d_4(:,:) = hco%lat2d_4(:,:) * MPC_RADIANS_PER_DEGREE_R8
    hco%lon2d_4(:,:) = hco%lon2d_4(:,:) * MPC_RADIANS_PER_DEGREE_R8
  
    deallocate(lat_8)
    deallocate(lon_8)
    
    !- 3.1 Compute maxGridSpacing 
    if ( trim(grtyp) == 'U' ) then
      latIndexBeg = 1
      latIndexEnd = nj / 2
    else
      latIndexBeg = 1
      latIndexEnd = nj
    end if
    
    maxDeltaLat = maxval( abs(hco%lat2d_4(2:ni,(latIndexBeg+1):latIndexEnd) - &
         hco%lat2d_4(1:(ni-1),latIndexBeg:(latIndexEnd-1))) )
    maxDeltaLon = 0.0d0
    do lonIndex = 1, ni - 1
      do latIndex = latIndexBeg, latIndexEnd - 1
        deltaLon = abs(hco%lon2d_4(lonIndex+1,latIndex+1) - hco%lon2d_4(lonIndex,latIndex))
        
        if ( deltaLon > MPC_PI_R8 ) deltaLon = deltaLon - 2.0d0 * MPC_PI_R8 
        
        deltaLon = deltaLon * cos(hco%lat2d_4(lonIndex,latIndex))
      
        if ( deltaLon > maxDeltaLon ) maxDeltaLon = deltaLon
      end do
    end do

    maxGridSpacing = EC_RA * sqrt(2.0d0) * max(maxDeltaLon,maxDeltaLat)
  
    if ( mpi_myid == 0 .and. maxGridSpacing /= maxGridSpacingPrevious ) then
      maxGridSpacingPrevious = maxGridSpacing
      write(*,*) 'hco_setupFromFile: maxDeltaLat=', maxDeltaLat * MPC_DEGREES_PER_RADIAN_R8, ' deg'
      write(*,*) 'hco_setupFromFile: maxDeltaLon=', maxDeltaLon * MPC_DEGREES_PER_RADIAN_R8, ' deg'
      write(*,*) 'hco_setupFromFile: maxGridSpacing=', maxGridSpacing, ' m'
    end if
  
    if ( maxGridSpacing > 1.0d6 ) then
      call utl_abort('hco_setupFromFile: maxGridSpacing is greater than 1000 km.')
    end if
    
    hco%maxGridSpacing = maxGridSpacing
  
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
    !
    ! :Purpose: to decide if a given grid is global or lam from input longitude array
    !           
    implicit none
    ! arguments:
    integer, intent(in)  :: ni
    real(8), intent(in)  :: lon(ni)
    logical, intent(out) :: global
    ! locals:
    real(8) :: dx, next_lon
    
    dx       = lon(2) - lon(1)
    next_lon = lon(ni) + 1.5d0 * dx
    
    write(*,*)
    write(*,*) 'dx       = ',dx
    write(*,*) 'lon(ni)  = ',lon(ni)
    write(*,*) 'next_lon = ',next_lon
    write(*,*) 'lon(1)   = ',lon(1)
    
    if ( next_lon - lon(1) > 360.0d0 .or. &
         next_lon - lon(1) < 3.0*dx ) then
      
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
    !
    ! :Purpose: to broadcast hco strucure from MPI task 0 to other tasks 
    !        
    implicit none
    ! arguments:
    type(struct_hco), pointer :: hco
    ! locals:
    integer :: ierr
    integer, external :: ezqkdef
    
    write(*,*) 'hco_mpiBcast: starting'
    
    if ( mpi_myid > 0 ) then
      if( .not.associated(hco) ) then
        allocate(hco)
      else
        call utl_abort('hco_mpiBcast: hco must be nullified for mpi task id > 0')
      end if
    end if
    
    call rpn_comm_bcastc(hco%gridname, len(hco%gridname), 'MPI_CHARACTER', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%initialized, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%ni, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%nj, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcastc(hco%grtyp, len(hco%grtyp), 'MPI_CHARACTER', 0, 'GRID', ierr)
    call rpn_comm_bcastc(hco%grtypTicTac, len(hco%grtypTicTac), 'MPI_CHARACTER', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%ig1, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%ig2, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%ig3, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%ig4, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    if ( mpi_myid > 0 ) then
      allocate(hco%lat(hco%nj))
      allocate(hco%lon(hco%ni))
      allocate(hco%lat2d_4(hco%ni,hco%nj))
      allocate(hco%lon2d_4(hco%ni,hco%nj))
    end if
    call rpn_comm_bcast(hco%lat, size(hco%lat), 'MPI_REAL8', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%lon, size(hco%lon), 'MPI_REAL8', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%lat2d_4, size(hco%lat2d_4), 'MPI_REAL4', 0, 'GRID', ierr)
    call rpn_comm_bcast(hco%lon2d_4, size(hco%lon2d_4), 'MPI_REAL4', 0, 'GRID', ierr)
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
      end if
      call rpn_comm_bcast(hco%tictacU, size(hco%tictacU), 'MPI_REAL4', 0, 'GRID', ierr)
    end if
    
    if ( mpi_myid > 0 ) then
      if ( hco%grtyp == 'G' ) then
        hco%EZscintID  = ezqkdef( hco%ni, hco%nj, hco%grtyp, hco%ig1, hco%ig2, hco%ig3, hco%ig4, 0 )
      else
        ! special treatment since EZscintID not properly communicated: keep as is
        write(*,*) 'hco_mpiBcast: Warning! Grid ID for EZSCINT not communicated for grtyp = ', hco%grtyp
        write(*,*) 'hco_mpiBcast: Warning! Grid ID for EZSCINT equals = ', hco%EZscintID
      end if
    end if

    write(*,*) 'hco_mpiBcast: done'

  end subroutine hco_mpiBcast

  !--------------------------------------------------------------------------
  ! Equal ?
  !--------------------------------------------------------------------------
  function hco_equal(hco1, hco2) result(equal)
    !
    ! :Purpose: to check if two given hco strucures are equal or not
    !        
    implicit none
    ! arguments:
    type(struct_hco), pointer :: hco1, hco2
    logical                   :: equal

    equal = .true.
    equal = equal .and. (hco1%ni == hco2%ni)
    equal = equal .and. (hco1%nj == hco2%nj)
    if (.not. equal) then
      write(*,*) 'hco_equal: dimensions not equal ', hco1%ni, hco1%nj, hco2%ni, hco2%nj
      return
    end if

    equal = equal .and. (hco1%grtyp   ==    hco2%grtyp)
    if (.not. equal) then
      write(*,*) 'hco_equal: grid type not equal'
      return
    end if

    equal = equal .and. (hco1%dlat == hco2%dlat)
    equal = equal .and. (hco1%dlon == hco2%dlon)
    if (.not. equal) then
      write(*,*) 'hco_equal: grid spacing not equal'
      return
    end if
    
    if(hco1%grtyp == 'G') then
      equal = equal .and. (hco1%ig2 == hco2%ig2)
      if (.not. equal) then
        write(*,*) 'hco_equal: Gaussian grid ig2 not equal'
        return
      end if
    end if
    
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
    end if
    
    equal = equal .and. all(hco1%lat(:) == hco2%lat(:))
    equal = equal .and. all(hco1%lon(:) == hco2%lon(:))
    if (.not. equal) then
      write(*,*) 'hco_equal: lat/lon not equal'
      return
    end if
    
  end function hco_equal

  !--------------------------------------------------------------------------
  ! hco_deallocate
  !--------------------------------------------------------------------------
  subroutine hco_deallocate( hco )
    implicit none
    type(struct_hco), pointer :: hco

    if (allocated(hco % lat)) deallocate(hco % lat)
    if (allocated(hco % lon)) deallocate(hco % lon)
    if (allocated(hco % lat2d_4)) deallocate(hco % lat2d_4)
    if (allocated(hco % lon2d_4)) deallocate(hco % lon2d_4)
    if (allocated(hco % tictacU)) deallocate(hco % tictacU)

    if (associated(hco)) deallocate(hco)
    nullify(hco)

  end subroutine hco_deallocate

  !--------------------------------------------------------------------------
  ! grid_mask
  !--------------------------------------------------------------------------
  subroutine grid_mask (F_mask_8,dx,dy,xg,yg,ni,nj)
    !
    ! :Purpose: 1) Find out where YIN lat lon points are in (YAN) grid with call to smat.
    !           2) If they are not outside of Yin grid, put area to zero for those points.
    !
    ! Author Qaddouri
    !
    implicit none

    ! arguments:
    integer,  intent(in)  :: Ni,Nj   
    real(8) , intent(out) :: F_mask_8(Ni,Nj)
    real(8) :: dx, dy
    real    :: xg(ni), yg(nj)

    ! locals: 
    integer :: lonIndex,latIndex,np_subd
    real(8)  :: poids(ni,nj),x_a_4,y_a_4,sp,sf,sp1,sf1
    real     :: area_4(ni,nj)

    np_subd = 4*ni

    sp    = 0.d0
    sf    = 0.d0

    do latIndex = 1, nj
      y_a_4 = yg(latIndex)
      do lonIndex = 1, ni

        x_a_4 = xg(lonIndex)-acos(-1.d0)

        area_4(lonIndex,latIndex) = dx*dy*cos(yg(latIndex))
        poids (lonIndex,latIndex) = yyg_weight (x_a_4,y_a_4,dx,dy,np_subd)

        !Check if poids <0
        if (poids(lonIndex,latIndex)*(1.d0-poids(lonIndex,latIndex)) > 0.d0) then
          sp = sp + poids(lonIndex,latIndex)*area_4(lonIndex,latIndex)
        else if (abs(poids(lonIndex,latIndex)-1.d0) < 1.d-14) then
          sf = sf + poids(lonIndex,latIndex)*area_4(lonIndex,latIndex)
        end if

      end do
    end do

    !Correct and scale poids
    !-----------------------
    sp1 = 0.d0
    sf1 = 0.d0

    do latIndex = 1, nj
      do lonIndex = 1, ni

        x_a_4 = poids(lonIndex,latIndex)*(2.d0*acos(-1.d0) - sf)/sp

        if (poids(lonIndex,latIndex)*(1.d0-poids(lonIndex,latIndex)) > 0.d0) then
          poids(lonIndex,latIndex) = min( 1.0d0, x_a_4 )
        end if
        if (poids(lonIndex,latIndex)*(1.0-poids(lonIndex,latIndex)) > 0.d0) then
          sp1 = sp1 + poids(lonIndex,latIndex)*area_4(lonIndex,latIndex)
        else if (abs(poids(lonIndex,latIndex)-1.d0) < 1.d-14) then
          sf1 = sf1 + poids(lonIndex,latIndex)*area_4(lonIndex,latIndex)
        end if

      end do
    end do

    !Correct
    !-------
    do latIndex = 1, nj
      do lonIndex = 1, ni
        x_a_4 = poids(lonIndex,latIndex)*(2.d0*acos(-1.d0) - sf1)/sp1

        if (poids(lonIndex,latIndex)*(1.d0-poids(lonIndex,latIndex)) > 0.d0) then
          poids(lonIndex,latIndex) = min( 1.d0, x_a_4 )
        end if
 
      end do
    end do

    F_mask_8 = 0.d0
    do latIndex=1,nj
      do lonIndex = 1,ni
        F_mask_8(lonIndex,latIndex) = poids(lonIndex,latIndex)
      end do
    end do

  end subroutine grid_mask

  !--------------------------------------------------------------------------
  ! inter_curve_boundary_yy
  !--------------------------------------------------------------------------
  subroutine inter_curve_boundary_yy (x,y,xi,yi,np)
    ! 
    ! :Purpose: compute the intersections between a line and the panel
    !         (yin or yang) boundary. The line passes through the panel
    !         center point (0, 0) and the cell center point (x, y).
    !
    ! Author A. Qaddouri. October 2016
    !
    ! Note: this routine has been taken and adjusted from a routine
    !       with the same name in the GEM model.
    !
    ! :Arguments: input:  (x, y):   longitude, latitude of the cell 
    !                              center point,
    !                     np:      workspace,
    !            output: (xi, yi): longitude, latitude of the 
    !                              intersection point.

    implicit none
    ! arguments:
    real(8) :: x,y,xi,yi
    integer :: np

    ! locals:
    real(8) :: tol, pi, xmin, ymin, xb
    real(8) :: xc, yc, s1, s2, x1, test, dxs
    real(8) :: xp1, xp2, xr1, yr1, xr2, yr2
    integer :: i

    tol = 1.0d-16
    pi  = MPC_PI_R8 
    xmin = -3.d0*pi/4.d0
    ymin = -pi/4.d0
    xb  = -0.5d0*pi

    xc = x
    yc = y

    if ( x > 0.d0 ) xc = - x
    if ( y > 0.d0 ) yc = - y

    if ( abs(xc) < tol ) then
      xi = xc
      yi = ymin
    else

      s1 = yc / xc
      s2 = ymin/(xb-xmin)

      x1 = s2*xmin/(s2-s1)
      if (x1 > xb ) then
        xi = ymin/s1
        yi = ymin
      else
        test = -1.d0
        dxs  = -xb/(np-1)
        i = 1
        do while (test < 0.d0 )
          xp1 = (i-1)*dxs + xb
          xp2 =   (i)*dxs + xb
          xr1 = atan2(sin(ymin),-cos(ymin)*cos(xp1))
          yr1 = asin(cos(ymin)*sin(xp1))
          xr2 = atan2(sin(ymin),-cos(ymin)*cos(xp2))
          yr2 = asin(cos(ymin)*sin(xp2))
          s2 = (yr1-yr2)/(xr1-xr2)
          xi = (s2*xr2-yr2)/(s2-s1)
          yi = s1*xi
          test=(xi-xr1)*(xr2-xi)
          i = i+1
        end do
      end if
    end if

    if ( x > 0.d0 ) xi = - xi
    if ( y > 0.d0 ) yi = - yi

  end subroutine inter_curve_boundary_yy

  !--------------------------------------------------------------------------
  ! hco_weight
  !--------------------------------------------------------------------------
  subroutine hco_weight(hco, weight)
    ! 
    ! :Purpose: given the horizontal grid definition of the grid,
    !         return appropriate weights for individual points (avoiding 
    !         double counting in the overlap regions in the case of a Yin-Yang grid). 
    !
    ! author: Abdessamad Qaddouri and Peter Houtekamer
    !         October 2016
    !
    ! Revision: imported code for the Yin-Yang grid on May 2021 from the EnKF library and 
    !           combined with code for other grid types from the MIDAS library.
    !
    ! :Arguments: 
    !    input:
    !        hco: structure with the specification of the horizontal grid 
    !    output:        
    !        weight: weight to be given when computing a horizontal average
    !
    implicit none

    ! arguments:
    type(struct_hco), intent(in) :: hco
    real(8), intent(out) :: weight(:,:)
    
    ! locals:
    integer :: sindx
    integer :: ni,nj
    integer :: lonIndex,latIndex,lonIndexP1,latIndexP1

    real(8),  allocatable :: F_mask_8(:,:), F_mask(:,:)

    real(8)  :: deg2rad,dx,dy,sum_weight
    real(8)  :: lon1,lon2,lon3,lat1,lat2,lat3
    real , allocatable :: xg(:),yg(:)
               
    deg2rad= MPC_RADIANS_PER_DEGREE_R8 
    sindx  = 6

    if (trim(hco%grtyp) == 'U') then ! case of a Yin-Yang grid
      write(*,*) 'compute weights for Yin_Yang grid'      
      ni = nint(hco%tictacU(sindx  ))
      nj = nint(hco%tictacU(sindx+1))
      allocate (F_mask_8 (ni,nj))
      allocate (F_mask(ni,2*nj))
      allocate (xg(ni))
      allocate (yg(nj))
         
      dx = hco%tictacU(sindx+10+1) -hco%tictacU(sindx+10)
      dy=  hco%tictacU(sindx+10+ni+1)-hco%tictacU(sindx+10+ni)
      dx=  deg2rad* dx
      dy=  deg2rad* dy
      do lonIndex=1,ni
        xg(lonIndex)=deg2rad* hco%tictacU(sindx+10+lonIndex-1)
      end do
      do latIndex=1,nj
        yg(latIndex)=  deg2rad*hco%tictacU(sindx+10+ni+latIndex-1)
        weight (:,latIndex) = cos( deg2rad* hco%tictacU(sindx+10+ni+latIndex-1))
        weight (:,nj+latIndex)= weight (:,latIndex)
      end do
      call  grid_mask (F_mask_8,dx,dy,xg,yg,ni,nj)
      do latIndex=1,nj
        F_mask (:,latIndex) = F_mask_8(:,latIndex)
        F_mask (:,nj+latIndex)= F_mask (:,latIndex)
      end do
      do latIndex=1,nj*2
        weight(:,latIndex) = weight(:, latIndex) * F_mask(:, latIndex)
      end do
      sum_weight=sum(weight)
      weight = weight/sum_weight
      deallocate(F_mask_8)
      deallocate(F_mask)
      deallocate(xg)
      deallocate(yg)
    else
      write(*,*) 'compute weights for grid type: ',hco%grtyp
      do latIndex=1,hco%nj
        latIndexP1 = min(hco%nj,latIndex+1)
        do lonIndex=1,hco%ni
          lonIndexP1 = min(hco%ni,lonIndex+1)
          lon1 = hco%lon2d_4(lonIndex,latIndex)
          lon2 = hco%lon2d_4(lonIndexP1,latIndex)
          lon3 = hco%lon2d_4(lonIndex,latIndexP1)
          lat1 = hco%lat2d_4(lonIndex,latIndex)
          lat2 = hco%lat2d_4(lonIndexP1,latIndex)
          lat3 = hco%lat2d_4(lonIndex,latIndexP1)
          dx = phf_calcDistance(lat1, lon1, lat2, lon2)/1000.0D0
          dy = phf_calcDistance(lat1, lon1, lat3, lon3)/1000.0D0
          weight(lonIndex,latIndex) = dx * dy
        end do
      end do
      sum_weight=sum(weight)
      weight = weight/sum_weight      
    end if

  end subroutine hco_weight

  !--------------------------------------------------------------------------
  ! yyg_weight
  !--------------------------------------------------------------------------
  real(8) function yyg_weight (x,y,dx,dy,np)
    Implicit none

    ! arguments:
    real(8) :: x,y,dx,dy
    integer :: np

    !
    ! :Purpose:
    !    Based on a Draft from Zerroukat (2013) - Evaluates weight
    !    for each cell, so that the overlap is computed once.
    !    
    !author
    !     Author Abdessamad Qaddouri -- Summer 2014
    !
    ! Note: this routine has been taken and adjusted from a routine
    !       with the same name in the GEM model.
    ! GEM revision:
    ! v4_70 - Qaddouri A.     - initial version
    !
    ! :Arguments: input: x,y: (longitude, latitude) in radians of
    !                           the cell center point.
    !                   dx, dy: horizontal grid resolution in radiancs.
    !                   np:  working dimension
    !            output: cell weight. 

    ! locals:
    real(8) :: pi, xmin, xmax, ymin, ymax, xb1, xb2
    real(8) :: t1x, t2x, t1y, dcell, xi, yi, di, dp, df, d

    pi   = MPC_PI_R8 
    xmin = -3.d0*pi/4.d0 ;  xmax = 3.d0*pi/4.d0
    ymin = -pi/4.d0       ;  ymax = pi/4.d0
    xb1  = -0.5d0*pi       ;  xb2  = 0.5d0*pi

    t1x = (x-xmin)*(xmax-x)
    t2x = (x-xb1)*(xb2-x)
    t1y = (y-ymin)*(ymax-y)

    if ( t1x < 0.d0 .or. t1y < 0.d0 ) then
      yyg_weight = 0.0
    else if ( t2x > 0.d0 .and. t1y > 0.d0 ) then
      yyg_weight = 1.d0
    else
      dcell = 0.5d0*dsqrt(dx**2 + dy**2)

      call inter_curve_boundary_yy (x, y, xi, yi, np)

      di = sqrt( xi**2 + yi**2 )
      dp = sqrt(  x**2 +  y**2 )
      df = dp - di
      d  = min(max(-dcell,df),dcell)
      yyg_weight = 0.5d0*(1.d0 - (d/dcell))
    end if

  end function yyg_weight

  !--------------------------------------------------------------------------
  ! hco_setupYgrid
  !--------------------------------------------------------------------------
  subroutine hco_setupYgrid( hco, ni, nj)
    !
    ! :Purpose: to initialize hco structure for a Y grid
    !           
    implicit none
    ! arguments:
    type(struct_hco), pointer :: hco
    integer, intent(in)       :: ni, nj

    allocate(hco)
    if (mpi_myId == 0 ) then
      hco%initialized = .true.
      hco%ni = ni
      hco%nj = nj
      hco%grtyp = 'Y'
      hco%grtypTicTac = 'L'
      if (allocated(hco%lat2d_4) ) then
        deallocate(hco%lat2d_4)
        deallocate(hco%lon2d_4) 
      end if
      allocate(hco%lat2d_4(ni, nj))
      allocate(hco%lon2d_4(ni, nj)) 
      hco%xlat1 = 0.d0
      hco%xlon1 = 0.d0
      hco%xlat2 = 1.d0 
      hco%xlon2 = 1.d0
    end if

  end subroutine hco_setupYgrid

end module HorizontalCoord_mod
