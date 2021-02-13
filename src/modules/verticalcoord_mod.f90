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

module verticalCoord_mod
  ! MODULE verticalcoord (prefix='vco' category='7. Low-level data objects and utilities')
  !
  ! :Purpose: Derived type and procedures related to the vertical levels.
  !           The derived type includes a pointer to the associated VGRID
  !           descriptor.
  !
  use mpi_mod
  use MathPhysConstants_mod
  use Vgrid_Descriptors
  use varNameList_mod
  use utilities_mod
  implicit none
  private

  ! public derived type
  public :: struct_vco
  ! public variables
  public :: vco_ip1_other
  ! public procedures
  public :: vco_setupFromFile, vco_getNumLev, vco_equal, vco_deallocate, vco_mpiBcast
  public :: vco_ensureCompatibleTops
  public :: vco_subsetOrNot, vco_levelMatchingList

  ! public entities accessed through inheritance (from module vgrid_descriptors)
  public :: vgd_get,vgd_levels,vgd_ok,vgd_dpidpis,vgd_write

  integer, parameter :: maxNumOtherLevels = 20
  integer :: vco_ip1_other(maxNumOtherLevels)

  integer, parameter :: maxNumDepthLevels = 200

  type struct_vco
     logical :: initialized=.false.
     integer :: Vcode = -1
     integer :: nlev_T = 0
     integer :: nlev_M = 0
     integer :: nlev_Other(vnl_numvarmaxOther) = 0
     integer :: nlev_Depth = 0
     integer :: ip1_sfc   ! ip1 value for the surface (hybrid = 1)
     integer :: ip1_T_2m  ! ip1 value for the 2m thermodynamic level
     integer :: ip1_M_10m ! ip1 value for the 10m momentum level
     integer,pointer :: ip1_T(:) => null()
     integer,pointer :: ip1_M(:) => null()  ! encoded IP1 levels (Thermo/Moment)
     integer :: ip1_depth(maxNumDepthLevels) = 0  ! encoded IP1 levels (Ocean depth levels)
     type(vgrid_descriptor) :: vgrid
     logical :: vgridPresent
     real(8) :: depths(maxNumDepthLevels)
  end type struct_vco

contains
  
  !--------------------------------------------------------------------------
  ! vco_allocateIp1
  !--------------------------------------------------------------------------
  subroutine vco_allocateIp1(vco)
    !
    ! :Purpose: Allocate the ip1 arrays of a vertical coordinate object.
    !
    implicit none

    ! Arguments:
    type(struct_vco), pointer :: vco ! Vertical coordinate object

    ! locals
    integer :: numLev

    numLev = vco_getNumLev(vco,'MM')
    if (numLev > 0) allocate (vco%ip1_M(numLev))

    numLev = vco_getNumLev(vco,'TH')
    if (numLev > 0) allocate (vco%ip1_T(numLev))

  end subroutine vco_allocateIp1

  !--------------------------------------------------------------------------
  ! vco_SetupFromFile
  !--------------------------------------------------------------------------
  subroutine vco_SetupFromFile(vco,templatefile,etiket_opt,beSilent_opt)
    ! 
    ! :Purpose: Initialize vertical coordinate object with information from 
    !           a standard file.
    !
    implicit none
    type(struct_vco),pointer   :: vco          ! Vertical coordinate object 
    character(len=*)           :: templatefile ! Template file
    character(len=*), optional :: etiket_opt   ! Optional argument etiket
    logical, optional          :: beSilent_opt ! Optional argument beSilent

    logical           :: beSilent
    character(len=12) :: etiket
    integer :: Vcode,jlev,nlevMatched,stat,nultemplate,ierr,ikey
    integer, parameter :: maxnumRecords = 500
    integer :: recordIndex, numRecords, ikeys(maxnumRecords)
    integer :: fnom,fstouv,fstfrm,fclos,fstinf,fstprm,fstinl
    integer :: vgd_nlev_M, vgd_nlev_T
    integer,   pointer :: vgd_ip1_M(:), vgd_ip1_T(:)
    integer :: ip1_sfc, ip1SeaLevel
    character(len=10) :: blk_S
    logical :: isExist_L, ip1_found, sfcFieldFound, oceanFieldFound
    integer :: ni, nj, nk, varListIndex, IP1kind
    character(len=4) :: nomvar_T, nomvar_M, nomvar_Other
    character(len=10) :: IP1string
    real :: otherVertCoordValue, vertCoordValue

    integer :: ideet, inpas, dateStamp_origin, ini, inj, ink, inbits, idatyp
    integer :: ip1, ip2, ip3, ig1, ig2, ig3, ig4, iswa, ilng, idltf, iubc
    integer :: iextra1, iextra2, iextra3
    character(len=2)  :: typvar
    character(len=4)  :: nomvar
    character(len=1)  :: grtyp

    nullify(vgd_ip1_M, vgd_ip1_T)

    if ( associated(vco) ) then
      call utl_abort('vco_setupFromFile: the supplied vco pointer is not null!')
    end if

    allocate(vco)

    if (present(etiket_opt)) then
      etiket = etiket_opt
    else
      etiket = ' '
    end if

    if (present(beSilent_opt)) then
      beSilent = beSilent_opt
    else
      beSilent = .false.
    end if

    if (mpi_myid == 0 .and. .not. beSilent) then
      write(*,*) 'vco_setupFromFile: Template File = ', trim(templatefile)
    end if
    inquire(file=templatefile,exist=isExist_L)
    if ( isExist_L )then
      nultemplate=0
      ierr=fnom(nultemplate,templatefile,'RND+OLD+R/O',0)
      if ( ierr == 0 ) then
        ierr =  fstouv(nultemplate,'RND+OLD')
      else
        call utl_abort('vco_setupFromFile: CANNOT OPEN TEMPLATE FILE!')
      end if
    else
      call utl_abort('vco_setupFromFile: CANNOT FIND TEMPLATE FILE!')
    end if

    !==========================================================================
    ! Get vertical coordinate descriptors from standard file (vgd_new reads "!!" record)

    stat = vgd_new(vco%vgrid,unit=nultemplate,format="fst",ip1=-1,ip2=-1)
    if (stat == VGD_OK) then
      vco%vgridPresent = .true.
    else
      write(*,*) 'vco_setupFromFile: Problem with vgd_new, check if surface-only file'

      ierr = fstinl(nultemplate,ini,inj,ink,-1,etiket,-1,-1,-1,' ', &
                    ' ',ikeys,numRecords,maxnumRecords)
      if ( ikeys(1) <= 0 ) then
        call utl_abort('vco_setupFromFile: Could not find any records ' //  &
                       'in the supplied file')
      end if
      write(*,*) 'vco_setupFromFile: number of records found = ', numRecords
      sfcFieldFound   = .false.
      oceanFieldFound = .false.
      record_loop: do recordIndex = 1, numRecords
        ierr = fstprm(ikeys(recordIndex), dateStamp_origin, ideet, inpas, ini, inj, &
                      ink, inbits, idatyp, ip1, ip2, ip3, &
                      typvar, nomvar, etiket, grtyp, ig1, ig2, ig3, ig4, &
                      iswa, ilng, idltf, iubc, iextra1, iextra2, iextra3)

        ! abort if !! record present, even though vgd_new failed
        if (trim(nomvar) == '!!') call utl_abort('vco_setupFromFile: found !! record')

        ! ignore horizontal coordinate records
        if (trim(nomvar) == '^^' .or. trim(nomvar) == '>>' .or.  &
            trim(nomvar) == '^>') cycle record_loop

        ! ignore any mask records
        if (typvar == '@@') cycle record_loop

        ! check for record with surface data
        call convip(ip1_sfc, 1.0, 5, 2, blk_s, .false.) 
        call convip(ip1SeaLevel, 0.0, 0, 2, blk_s, .false.) 
        if (ip1 == 0 .or. ip1 == ip1_sfc .or. ip1 == ip1SeaLevel) then
          sfcFieldFound = .true.
          cycle record_loop
        end if

        ! check for record with ocean data
        call convip(ip1, vertCoordValue, Ip1Kind, -1, blk_s, .false.) 
        if (Ip1Kind == 0 .and. vertCoordValue > 0.0 .and. &
            vnl_varKindFromVarname(trim(nomvar)) == 'OC' ) then
          oceanFieldFound = .true.
          ! check if we've NOT already recorded this depth level
          if ( .not. any(vertCoordValue == vco%depths(1:vco%nLev_depth)) ) then
            vco%nLev_depth = vco%nLev_depth + 1
            vco%depths(vco%nLev_depth) = vertCoordValue
            vco%ip1_depth(vco%nLev_depth) = ip1
            if (mpi_myid == 0) then
              write(*,*) 'vco_setupFromFile: found ocean record on height levels'
              write(*,*) 'vco_setupFromFile: found ocean record: nLev_depth = ', &
                   vco%nLev_depth, 'varName = ', trim(nomvar), &
                   ', value = ', vco%depths(vco%nLev_depth)
            end if
          end if
          cycle record_loop
        end if

        ! something else was found, abort
        write(*,*) 'vco_setupFromFile: found non-surface record'
        write(*,*) 'varName = ', trim(nomvar), ' typvar = ', typvar, ' ip1 = ', ip1
        call utl_abort('vco_setupFromFile: found a non-surface field')

      end do record_loop

      ! Check if ocean depth levels are in correct order (ascending in value)
      if ( oceanFieldFound .and. vco%nLev_depth > 1 ) then
        if ( any(vco%depths(2:vco%nLev_depth)-vco%depths(1:(vco%nLev_depth-1)) < 0.0) ) then
          call utl_abort('vco_setupFromFile: some depth levels not in ascending order')
        end if
      end if

      if (.not.sfcFieldFound .and. .not.oceanFieldFound) then
        call utl_abort('vco_setupFromFile: no surface or ocean fields found')
      end if
      write(*,*) 'vco_setupFromFile: found only surface or ocean fields, proceed without vgrid descriptor'
      vco%vgridPresent = .false.
      vco%nlev_T = 0
      vco%nlev_M = 0
      vco%Vcode  = 0
      vco%initialized=.true.

      ierr =  fstfrm(nultemplate)
      ierr =  fclos (nultemplate)

      return ! skip the rest
    end if

    ! Print out vertical structure 
    if (mpi_myid == 0 .and. .not. beSilent) then
      call flush(6) ! possibly needed so vgd_print output appears correctly in listing
      stat = vgd_print(vco%vgrid)
      if ( stat /= VGD_OK )then
        call utl_abort('vco_setupFromFile: ERROR with vgd_print')
      end if
    end if

    !==========================================================================
    ! Get version of the vertical coordinate

    stat = 0
    stat = vgd_get(vco%vgrid,key='ig_1 - vertical coord code',value=Vcode)
    if ( stat /= VGD_OK ) then
      call utl_abort('vco_setupFromFile: problem with vgd_get: key= ig_1 - vertical coord code')
    end if
    if (Vcode /= 5002 .and. Vcode /= 5005) then
      call utl_abort('vco_setupFromFile: Invalid Vcode. Currently only 5002 and 5005 supported.')
    end if
    vco%Vcode = Vcode

    ! Get vgrid values for ip1
    stat = 0

    stat = vgd_get(vco%vgrid,key='vipm - vertical levels (m)',value=vgd_ip1_m)
    stat = stat + VGD_OK

    stat = vgd_get(vco%vgrid,key='vipt - vertical ip1 levels (t)',value=vgd_ip1_t)
    stat = stat + VGD_OK

    if (stat /= 0) then
      call utl_abort('vco_setupFromFile: problem with vgd_get')
    end if

    vgd_nlev_M = size(vgd_ip1_M)
    vgd_nlev_T = size(vgd_ip1_T)

    !==========================================================================
    ! Set the number of vertical levels and allocate vco arrays

    vco%nlev_T = 0
    nomvar_T = 'TH  '
    do jlev = 1, vgd_nlev_T
      ikey = fstinf(nultemplate, ni, nj, nk, -1 ,etiket, vgd_ip1_T(jlev), -1, -1, ' ', nomvar_T)
      if (ikey > 0) vco%nlev_T = vco%nlev_T + 1
    end do
    if (vco%nlev_T == 0) then
      if (mpi_myid == 0 .and. .not. beSilent) then
        write(*,*) 'vco_setupFromFile: TH not found looking for TT to get nlev_T'
      end if
      nomvar_T = 'TT  '
      do jlev = 1, vgd_nlev_T
        ikey = fstinf(nultemplate, ni, nj, nk, -1 ,etiket, vgd_ip1_T(jlev), -1, -1, ' ', nomvar_T)
        if (ikey > 0) vco%nlev_T = vco%nlev_T + 1
      end do
    end if
    if ( vco%nlev_T == 0 ) then
      write(*,*) 
      write(*,*) 'vco_setupfromfile: Could not find a valid thermodynamic variable in the template file!'
    end if

    vco%nlev_M = 0
    nomvar_M = 'MM  '
    do jlev = 1, vgd_nlev_M
      ikey = fstinf(nultemplate, ni, nj, nk, -1 ,etiket, vgd_ip1_M(jlev), -1, -1, ' ', nomvar_M)
      if (ikey > 0) vco%nlev_M = vco%nlev_M + 1
    end do
    if (vco%nlev_M == 0) then
      if (mpi_myid == 0 .and. .not. beSilent) then
        write(*,*) 'vco_setupFromFile: MM not found looking for UU to get nlev_M'
      end if
      nomvar_M = 'UU  '
      do jlev = 1, vgd_nlev_M
        ikey = fstinf(nultemplate, ni, nj, nk, -1 ,etiket, vgd_ip1_M(jlev), -1, -1, ' ', nomvar_M)
        if (ikey.gt.0) vco%nlev_M = vco%nlev_M + 1
      end do
    end if
    if (vco%nlev_M == 0) then
      if (mpi_myid == 0 .and. .not. beSilent) then
        write(*,*) 'vco_setupFromFile: UU not found looking for PP to get nlev_M'
      end if
      nomvar_M = 'PP  '
      do jlev = 1, vgd_nlev_M
        ikey = fstinf(nultemplate, ni, nj, nk, -1 ,etiket, vgd_ip1_M(jlev), -1, -1, ' ', nomvar_M)
        if ( ikey > 0 ) vco%nlev_M = vco%nlev_M + 1
      end do
    end if
    if (vco%nlev_M == 0) then
      write(*,*) 
      write(*,*) 'vco_setupfromfile: Could not find a valid momentum variable in the template file!'
    end if

    do jlev = 1, maxNumOtherLevels
      otherVertCoordValue = real(jlev)
      IP1kind = 3
      call convip_plus(vco_ip1_other(jlev), otherVertCoordValue, IP1kind, +2, IP1string, .false.)
    end do
    vco%nlev_Other(:) = 0
    do varListIndex = 1, vnl_numvarmaxOther
      nomvar_Other = vnl_varNameListOther(varListIndex)
      do jlev = 1, maxNumOtherLevels
        ikey = fstinf(nultemplate, ni, nj, nk, -1 ,etiket, vco_ip1_other(jlev), -1, -1, ' ', nomvar_Other)
        if (ikey > 0) vco%nlev_Other(varListIndex) = vco%nlev_Other(varListIndex) + 1
      end do
      if (mpi_myid == 0 .and. .not. beSilent) then
        if (vco%nlev_Other(varListIndex) == 0) then
          write(*,*) 
          write(*,*) 'vco_setupfromfile: Found no levels in template file for OTHER type variable ', nomvar_Other
        else
          write(*,*) 'vco_setupfromfile: Found ', vco%nlev_Other(varListIndex),  &
               ' levels in template file for OTHER type variable ', nomvar_Other
        end if
      end if
    end do

    if (vco%nlev_M == 0 .and. vco%nlev_T == 0) then
      call utl_abort('vco_setupfromfile: they were no valid momentum and thermodynamic variables in the template file!')
    end  if

    if (mpi_myid == 0 .and. .not.beSilent) then
      write(*,*) 'vco_setupFromFile: nlev_M, nlev_T=',vco%nlev_M,vco%nlev_T
    end if
    
    call vco_allocateIp1(vco)

    !==========================================================================
    ! Define levels ip1 for momentum levels

    ! Match up ip1 values from file and vgrid
    nlevMatched = 0
    do jlev = 1, vgd_nlev_M
      ikey = fstinf(nultemplate, ni, nj, nk, -1 ,etiket, vgd_ip1_M(jlev), -1, -1, ' ', nomvar_M)
      if ( ikey > 0 ) then
        nlevMatched = nlevMatched + 1
        if (nlevMatched.gt.vco%nlev_M) then
          call utl_abort('vco_setupFromFile: Problem with consistency between vgrid descriptor and template file (momentum)')
        end if
        vco%ip1_M(nlevMatched) = vgd_ip1_M(jlev)
      end if
    end do
    if (nlevMatched /= vco%nlev_M) then
      write(*,*) 'vco_setupFromFile: nlevMatched = ', nlevMatched, ', nlev_M = ', vco%nlev_M
      call utl_abort('vco_setupFromFile: Problem with consistency between vgrid descriptor and template file (momentum)')
    end if

    !==========================================================================
    ! Define levels ip1 for thermo levels

    ! Match up ip1 values from file and vgrid
    nlevMatched = 0
    do jlev = 1, vgd_nlev_T
      ikey = fstinf(nultemplate, ni, nj, nk, -1 ,etiket, vgd_ip1_T(jlev), -1, -1, ' ', nomvar_T)
      if ( ikey > 0 ) then
        nlevMatched = nlevMatched + 1
        if (nlevMatched > vco%nlev_T) then
          call utl_abort('vco_setupFromFile: Problem with consistency between vgrid descriptor and template file (thermo)')
        end if
        vco%ip1_T(nlevMatched) = vgd_ip1_T(jlev)
      end if
    end do
    if (nlevMatched /= vco%nlev_T) then
      write(*,*) 'vco_setupFromFile: nlevMatched = ', nlevMatched, ', nlev_T = ', vco%nlev_T
      call utl_abort('vco_setupFromFile: Problem with consistency between vgrid descriptor and template file (thermo)')
    end if

    !
    !- Define level ip1 for surface (only used for Vcode=5005)
    !

    ! determine IP1 of sfc (hyb=1.0)
    call convip(ip1_sfc, 1.0, 5, 2, blk_s, .false.) 
    ip1_found = .false.
    do jlev = 1, vgd_nlev_T
      if (ip1_sfc .eq. vgd_ip1_T(jlev)) then
        ip1_found = .true.
        vco%ip1_sfc = vgd_ip1_T(jlev)
      end if
    end do
    if (.not.ip1_found) then
      write(*,*) 'vco_setupFromFile: Could not find IP1=',ip1_sfc
      call utl_abort('vco_setupFromFile: No surface level found in Vgrid!!!')
    else
      if ( mpi_myid == 0 .and. .not. beSilent ) write(*,*) 'vco_setupFromFile: Set surface level IP1=',vco%ip1_sfc
    end if

    !
    !- determine IP1s of 2m and 10m levels
    !
    call set_2m_10m_levels(vco)

    !
    !- Ending
    !
    vco%initialized=.true.

    ierr =  fstfrm(nultemplate)
    ierr =  fclos (nultemplate)

  end subroutine vco_setupFromFile


  subroutine vco_deallocate(vco)
    !
    ! :Purpose: Deallocate vertical coordinate object
    !
    implicit none

    type(struct_vco), pointer :: vco ! Vertical coordinate object

    integer :: stat

    if ( vco%vgridPresent ) then
      deallocate (vco%ip1_M)
      deallocate (vco%ip1_T)
      stat = vgd_free(vco%vgrid)
    end if

    nullify(vco)

  end subroutine vco_deallocate


  function vco_getNumLev(vco,varLevel,varName_opt) result(nlev)
    ! 
    ! :Purpose: get number of vertical levels
    ! 
    implicit none

    ! arguments:
    type(struct_vco), pointer              :: vco         ! Vertical coordinate object
    character(len=*), intent(in)           :: varLevel    ! 'TH', 'MM', 'SF', 'SFMM',
                                                          ! 'SFTH', 'DP', 'SFDP' or 'OT'
    character(len=*), optional, intent(in) :: varName_opt ! only needed for varLevel='OT'
    ! locals:
    integer :: nlev, varListIndex

    if (varLevel == 'MM') then
      nlev = vco%nlev_M
    else if (varLevel == 'TH') then
      nlev = vco%nlev_T
    else if (varLevel == 'SF'   .or. varLevel == 'SFTH' .or. &
             varLevel == 'SFMM' .or. varLevel == 'SFDP') then
      nlev = 1
    else if (varLevel == 'OT') then
      if (.not. present(varName_opt)) then
        call utl_abort('vco_getNumLev: varName must be specified for varLevel=OT')
      end if
      varListIndex = vnl_varListIndexOther(varName_opt)
      nlev = vco%nlev_Other(varListIndex)
    else if (varLevel == 'DP') then
      nlev = vco%nlev_depth
    else
      call utl_abort('vco_getNumLev: Unknown variable type! ' // varLevel)
    end if

  end function vco_getNumLev


  subroutine vco_mpiBcast(vco)
    !
    ! :Purpose: MPI broadcast of vertical coordinate object
    !
    implicit none
    ! Arguments:
    type(struct_vco), pointer :: vco ! vertical coordinate object
    ! locals
    integer :: ierr, vgd_nlev_M, vgd_nlev_T
    integer :: vgdig1, vgdig2, vgdig3, vgdig4, vgdip1, vgdip2, vgdip3, vgddate
    integer :: vgdtable_dim1, vgdtable_dim2, vgdtable_dim3
    character(len=12) :: vgdetik
    real(8), pointer :: vgdtable(:,:,:)

    nullify(vgdtable)

    write(*,*) 'vco_mpiBcast: starting'

    if ( mpi_myid > 0 ) then
      if ( .not.associated(vco) ) then
        allocate(vco)
      else 
        call utl_abort('vco_mpiBcast: vco must be nullified for mpi task id > 0')
      end if
    end if

    call rpn_comm_bcast(vco%initialized , 1, 'MPI_LOGICAL', 0, 'GRID', ierr)
    call rpn_comm_bcast(vco%vgridPresent, 1, 'MPI_LOGICAL', 0, 'GRID', ierr)
    call rpn_comm_bcast(vco%nlev_T      , 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(vco%nlev_M      , 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(vco%ip1_sfc     , 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(vco%ip1_T_2m    , 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(vco%ip1_M_10m   , 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    call rpn_comm_bcast(vco%Vcode       , 1, 'MPI_INTEGER', 0, 'GRID', ierr)
    if (vco%vgridPresent) then
      if ( mpi_myid == 0 ) then
        vgd_nlev_M = size(vco%ip1_M)
        vgd_nlev_T = size(vco%ip1_T)
      end if
      call rpn_comm_bcast(vgd_nlev_M, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
      call rpn_comm_bcast(vgd_nlev_T, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
      if ( mpi_myid > 0 ) then
        allocate(vco%ip1_M(vgd_nlev_M))
        allocate(vco%ip1_T(vgd_nlev_T))
      end if
      call rpn_comm_bcast(vco%ip1_M, vgd_nlev_M, 'MPI_INTEGER', 0, 'GRID', ierr)
      call rpn_comm_bcast(vco%ip1_T, vgd_nlev_T, 'MPI_INTEGER', 0, 'GRID', ierr)
    end if

    ! now do bcast for vgrid object
    if (vco%vgridPresent) then

      if ( mpi_myid == 0 ) then
        ierr = vgd_get(vco%vgrid,'VTBL',vgdtable)
        vgdtable_dim1 = size(vgdtable,1)
        vgdtable_dim2 = size(vgdtable,2)
        vgdtable_dim3 = size(vgdtable,3)
        ierr = vgd_get(vco%vgrid,'DATE',vgddate)
        ierr = vgd_get(vco%vgrid,'ETIK',vgdetik)
        ierr = vgd_get(vco%vgrid,'IG_1',vgdig1)
        ierr = vgd_get(vco%vgrid,'IG_2',vgdig2)
        ierr = vgd_get(vco%vgrid,'IG_3',vgdig3)
        ierr = vgd_get(vco%vgrid,'IG_4',vgdig4)
        ierr = vgd_get(vco%vgrid,'IP_1',vgdip1)
        ierr = vgd_get(vco%vgrid,'IP_2',vgdip2)
        ierr = vgd_get(vco%vgrid,'IP_3',vgdip3)
      end if

      ! 3D table of real*8
      call rpn_comm_bcast(vgdtable_dim1, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
      call rpn_comm_bcast(vgdtable_dim2, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
      call rpn_comm_bcast(vgdtable_dim3, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
      if ( mpi_myid > 0 ) allocate(vgdtable(vgdtable_dim1, vgdtable_dim2, vgdtable_dim3)) 
      call rpn_comm_bcast(vgdtable, size(vgdtable), 'MPI_REAL8', 0, 'GRID', ierr)
      ! others
      call rpn_comm_bcast(vgddate, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
      call rpn_comm_bcastc(vgdetik, len(vgdetik), 'MPI_CHARACTER', 0, 'GRID', ierr)
      call rpn_comm_bcast(vgdig1, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
      call rpn_comm_bcast(vgdig2, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
      call rpn_comm_bcast(vgdig3, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
      call rpn_comm_bcast(vgdig4, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
      call rpn_comm_bcast(vgdip1, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
      call rpn_comm_bcast(vgdip2, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
      call rpn_comm_bcast(vgdip3, 1, 'MPI_INTEGER', 0, 'GRID', ierr)
      
      if ( mpi_myid > 0 ) then
        ierr = vgd_new(vco%vgrid,vgdtable)
        ierr = vgd_put(vco%vgrid,'DATE',vgddate)
        ierr = vgd_put(vco%vgrid,'ETIK',vgdetik)
        ierr = vgd_put(vco%vgrid,'IG_1',vgdig1)
        ierr = vgd_put(vco%vgrid,'IG_2',vgdig2)
        ierr = vgd_put(vco%vgrid,'IG_3',vgdig3)
        ierr = vgd_put(vco%vgrid,'IG_4',vgdig4)
        ierr = vgd_put(vco%vgrid,'IP_1',vgdip1)
        ierr = vgd_put(vco%vgrid,'IP_2',vgdip2)
        ierr = vgd_put(vco%vgrid,'IP_3',vgdip3)
      end if
      
      deallocate(vgdtable) 

    end if

    write(*,*) 'vco_mpiBcast: done'

  end subroutine vco_mpiBcast


  subroutine vco_ensureCompatibleTops(vco_sourceGrid,vco_destGrid) 
    !
    ! :Purpose: This function checks if the top of a destination grid
    !           is ~equal or lower than the top of the source grid
    !           the code aborts if this is not the case
    !
    implicit none

    ! Arguments:
    type(struct_vco), pointer :: vco_sourceGrid ! vertical coordinate source grid
    type(struct_vco), pointer :: vco_destGrid   ! vertical coordinate destination grid

    ! internal
    integer :: status, nAbove, numLevSource, numLevDest
    real(kind=8) :: sourceModelTop
    real(kind=8), dimension(1,1) :: pSfc
    real(kind=8), pointer, dimension(:,:,:) :: sourcePressureLevels
    real(kind=8), pointer, dimension(:,:,:) :: destPressureLevels

    nullify(sourcePressureLevels)
    nullify(destPressureLevels)

    numLevSource = vco_getNumLev(vco_sourceGrid,'MM')
    numLevDest   = vco_getNumLev(vco_destGrid,'MM')
    if (numLevSource == 1 .and. numLevDest == 1) then
      write(*,*) 'vco_ensureCompatibleTops: both grids only have 1 level, skip the test'
      return
    end if

    !dummy pressure value
    pSfc(1,1) = 100.0D3 !100 kPa

    !pressure on momentum levels of source grid
    status=vgd_levels(vco_sourceGrid%vgrid,          &
                      ip1_list=vco_sourceGrid%ip1_M, &
                      levels=sourcePressureLevels,   &
                      sfc_field=pSfc,                &
                      in_log=.false.)
    if ( status /= VGD_OK ) call utl_abort('vco_ensureCompatibleTops: ERROR with vgd_levels')

    !pressure on momentum levels of destination grid
    status=vgd_levels(vco_destGrid%vgrid,            &
                      ip1_list=vco_destGrid%ip1_M,   &
                      levels=destPressureLevels,     &
                      sfc_field=pSfc,                &
                      in_log=.false.)
    if ( status /= VGD_OK ) call utl_abort('vco_ensureCompatibleTops: ERROR with vgd_levels')

    !count number of levels where output grid is higher than input grid
    sourceModelTop = sourcePressureLevels(1,1,1)
    nAbove=0
    do while (sourceModelTop > destPressureLevels(1,1,nAbove+1))
      nAbove = nAbove + 1
    end do

    !Destination grid has "nAbove" levels above source grid;  tolerate one
    if ( nAbove > 1 ) then
      call utl_abort('vco_ensureCompatibleTops: top of destination grid is more than one level higher than the top of source grid')
    end if

  end subroutine vco_ensureCompatibleTops


  function vco_equal(vco1,vco2) result(equal)
    !
    ! :Purpose: Compare two vertical grid object and provide a logical result if they are equal or not
    !
    implicit none
  
    !Arguments:
    type(struct_vco), pointer :: vco1 ! vertical coordinate object one
    type(struct_vco), pointer :: vco2 ! vertical coordinate object two

    logical                   :: equal

    equal = .true.

    equal = equal .and. (vco1%Vcode == vco2%Vcode)
    if (.not. equal) then
      write(*,*) 'vco_equal: Vcode not equal'
      return
    end if
    if ( vco1%vgridPresent .and. vco2%vgridPresent ) then
       equal = equal .and. (vco1%vgrid == vco2%vgrid)
       if (.not. equal) then
          write(*,*) 'vco_equal: vgrid not equal'
          return
       end if
    end if

    ! Even if vgrid defined, not enough just to compare vgrid, must compare everything
    equal = equal .and. (vco1%nlev_T == vco2%nlev_T)
    if (.not. equal) then
       write(*,*) 'vco_equal: nlev_T not equal', vco1%nlev_T, vco2%nlev_T
       return
    end if
    equal = equal .and. (vco1%nlev_M == vco2%nlev_M)
    if (.not. equal) then
       write(*,*) 'vco_equal: nlev_M not equal', vco1%nlev_M, vco2%nlev_M
       return
    end if
    if ( vco1%vgridPresent .and. vco2%vgridPresent ) then
      equal = equal .and. all(vco1%ip1_T(:) == vco2%ip1_T(:))
      if (.not. equal) then
        write(*,*) 'vco_equal: ip1_T not equal'
        return
      end if
      equal = equal .and. all(vco1%ip1_M(:) == vco2%ip1_M(:))
      if (.not. equal) then
        write(*,*) 'vco_equal: ip1_M not equal'
        return
      end if
      equal = equal .and. (vco1%ip1_sfc == vco2%ip1_sfc)
      if (.not. equal) then
        write(*,*) 'vco_equal: ip1_sfc not equal'
        return
      end if
      if (vco1%Vcode == 5002 .or. vco1%Vcode == 5005) then
        equal = equal .and. hybridCoefEqualOrNot(vco1, vco2)
        if (.not. equal) then
          write(*,*) 'vco_equal: hybrid parameters are not equal'
          return
        end if
      end if
    end if

  end function vco_equal


  function vco_subsetOrNot(vco_template, vco_full) result(subset)
    !
    ! :Purpose: This function determines if vco_template is a subset of vco_full.
    !
    implicit none
    
    ! Arguments:
    type(struct_vco), pointer, intent(in)  :: vco_full     ! vertical coordinate object full
    type(struct_vco), pointer, intent(in)  :: vco_template ! vertical coordinate object template
    ! Result:
    logical :: subset

    integer, allocatable :: THlevelWanted(:), MMlevelWanted(:)

    !
    !- Compare the vCode
    !
    if (vco_template%Vcode /= vco_full%Vcode) then
      subset = .false.
      return
    end if

    !
    !- Compare the IP1s
    !
    allocate(THlevelWanted(vco_template%nlev_T))
    allocate(MMlevelWanted(vco_template%nlev_M))

    call vco_levelMatchingList( THlevelWanted, MMlevelWanted, & ! OUT
                                vco_template, vco_full )        ! IN

    if (any(THlevelWanted == -1) .or. any(MMlevelWanted == -1)) then
      subset = .false.
      return
    end if

    deallocate(MMlevelWanted)
    deallocate(THlevelWanted)

    !
    !- For hybrid coordinates, compare additional grid parameters
    !
    if ( .not. hybridCoefEqualOrNot(vco_template, vco_full) ) then
      subset = .false.
      return
    end if

    !
    !- When reaching this point, we assume that we have a subset
    !
    subset = .true.

  end function vco_subsetOrNot


  function hybridCoefEqualOrNot(vco1, vco2) result(equal)
    !
    ! :Purpose: To compare two vertical coordinate hybrid coefficient object 
    !
    implicit none

    ! Arguments:
    type(struct_vco), pointer, intent(in)  :: vco1 ! vertical coordinate object one
    type(struct_vco), pointer, intent(in)  :: vco2 ! vertical coordinate object two
    ! Result:
    logical :: equal

    real(8) :: ptop1, ptop2

    real(8), pointer :: coefA1(:), coefA2(:)
    real(8), pointer :: coefB1(:), coefB2(:)

    real :: coefR11, coefR12
    real :: coefR21, coefR22

    integer :: stat, levIndex

    if (vco1%Vcode == 5002) then
      !- Ptop
      stat = vgd_get(vco1%vgrid,key='PTOP - top level pressure',value=ptop1)
      stat = vgd_get(vco2%vgrid,key='PTOP - top level pressure',value=ptop2)
      if (ptop1 /= ptop2) then
        equal = .false.
        return
      end if
    end if

    if (vco1%Vcode == 5002 .or. vco1%Vcode == 5005) then

      !- Pref
      stat = vgd_get(vco1%vgrid,key='PREF - reference pressure',value=ptop1)
      stat = vgd_get(vco2%vgrid,key='PREF - reference pressure',value=ptop2)
      if (ptop1 /= ptop2) then
        equal = .false.
        return
      end if
      !- R-coef 1
      stat = vgd_get(vco1%vgrid,key='RC_1 - first R-coef value',value=coefR11)
      stat = vgd_get(vco2%vgrid,key='RC_1 - first R-coef value',value=coefR12)
      if (coefR11 /= coefR12) then
        equal = .false.
        return
      end if
      !- R-coef 2
      stat = vgd_get(vco1%vgrid,key='RC_2 - second R-coef value',value=coefR21)
      stat = vgd_get(vco2%vgrid,key='RC_2 - second R-coef value',value=coefR22)
      if (coefR21 /= coefR22) then
        equal = .false.
        return
      end if
      !- A
      nullify(coefA1)
      nullify(coefA2)
      stat = vgd_get(vco1%vgrid,key='CA_M - vertical A coefficient (m)',value=coefA1)
      stat = vgd_get(vco2%vgrid,key='CA_M - vertical A coefficient (m)',value=coefA2)
      if ( size(coefA1) /= size(coefA2) ) then
        equal = .false.
        return
      end if
      do levIndex = 1, size(coefA1)
        if (coefA1(levIndex) /= coefA2(levIndex)) then
          equal = .false.
          return
        end if
      end do
      !- B
      nullify(coefB1)
      nullify(coefB2)
      stat = vgd_get(vco1%vgrid,key='CB_M - vertical B coefficient (m)',value=coefB1)
      stat = vgd_get(vco2%vgrid,key='CB_M - vertical B coefficient (m)',value=coefB2)
      if ( size(coefB1) /= size(coefB2) ) then
        equal = .false.
        return
      end if
      do levIndex = 1, size(coefB1)
        if (coefB1(levIndex) /= coefB2(levIndex)) then
          equal = .false.
          return
        end if
      end do

    end if

    !
    !- When reaching this point, we assume that they are equal
    !
    equal = .true.

  end function hybridCoefEqualOrNot


  subroutine vco_levelMatchingList(THmatchingList, MMmatchingList, vco1, vco2)
    !
    ! :Purpose: This subroutine returns arrays of array indices of the levels (ip1s) in vco2 
    !           corresponding with the levels (ip1s) in vco1
    !
    implicit none

    ! Arguments:
    type(struct_vco), pointer, intent(in) :: vco1       ! vertical coordinate object one
    type(struct_vco), pointer, intent(in) :: vco2       ! vertical coordinate object two
    integer, intent(out) :: THmatchingList(vco1%nlev_T) ! TH matching list
    integer, intent(out) :: MMmatchingList(vco1%nlev_M) ! MM matching list

    ! locals
    integer :: levIndex1, levIndex2

    !
    !- Do momentum levels...
    !
    MMmatchingList(:) = -1
    do levIndex1 = 1, vco1%nlev_M
      do levIndex2 =  1, vco2%nlev_M
        if ( (vco2%ip1_M(levIndex2) == vco1%ip1_M(levIndex1)) .or. &
             (vco2%ip1_M(levIndex2) == vco2%ip1_M_10m .and.        &
              vco1%ip1_M(levIndex1) == vco1%ip1_M_10m) ) then
          MMmatchingList(levIndex1) = levIndex2
          exit
        end if
      end do
    end do
    
    !
    !- Do thermo levels...
    !
    THmatchingList(:) = -1
    do levIndex1 = 1, vco1%nlev_T
      do levIndex2 = 1, vco2%nlev_T
        if ( (vco2%ip1_T(levIndex2) == vco1%ip1_T(levIndex1)) .or. &
             (vco2%ip1_T(levIndex2) == vco2%ip1_T_2m .and.        &
              vco1%ip1_T(levIndex1) == vco1%ip1_T_2m) ) then 
          THmatchingList(levIndex1) = levIndex2
          exit
        end if
      end do
    end do

  end subroutine vco_levelMatchingList


  subroutine set_2m_10m_levels(vco)
    !
    ! :Purpose: To set 2-m and 10-m levels
    !
    implicit none

    type(struct_vco), pointer, intent(in) :: vco ! vertical coordinate object

    character(len=10) :: blk_s

    if      (vco%Vcode == 5002) then
      vco%ip1_T_2m  = vco%ip1_sfc 
      vco%ip1_M_10m = vco%ip1_sfc
    else if (vco%Vcode == 5005) then
      call convip(vco%ip1_T_2m ,  1.5, 4, 2, blk_s, .false.)
      call convip(vco%ip1_M_10m, 10.0, 4, 2, blk_s, .false.)
    else
      vco%ip1_T_2m  = -1
      vco%ip1_M_10m = -1
    end if

  end subroutine set_2m_10m_levels

end module verticalCoord_mod
