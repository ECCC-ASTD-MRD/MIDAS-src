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

module lamBMatrixHI_mod
  ! MODULE lamBMatrixHI_mod (prefix='lbhi' category='2. B and R matrices')
  !
  ! :Purpose: Performs transformation from control vector to analysis increment 
  !           using the homogeneous and isotropic background error covariance 
  !           matrix.
  !
  use mpi_mod
  use horizontalCoord_mod
  use verticalCoord_mod
  use LamSpectralTransform_mod
  use gridStateVector_mod
  use analysisGrid_mod
  use utilities_mod
  use gridVariableTransforms_mod
  use varNameList_mod
  implicit none
  save
  private

  ! public procedures
  public :: lbhi_Setup, lbhi_bSqrt, lbhi_bSqrtAdj, lbhi_Finalize
  public :: lbhi_expandToMPIglobal, lbhi_expandToMPIglobal_r4, lbhi_reduceToMPIlocal, lbhi_reduceToMPIlocal_r4

  integer, parameter :: cv_model = 1
  integer, parameter :: cv_bhi   = 2
  type  :: lbhi_cv
    character(len=4)      :: NomVar(2)
    character(len=2)      :: GridType ! TT=Thermo, MM=Momentum, NS=Non-staggered
    integer               :: nlev
    integer               :: kDimStart
    integer               :: kDimEnd
    real(8), allocatable  :: GpStdDev(:,:,:)
    integer, allocatable  :: ip1(:)
  end type lbhi_cv

  integer,parameter    :: nMaxControlVar = 10
  type(lbhi_cv)        :: ControlVariable(nMaxControlVar)

  integer :: UWindID = -1
  integer :: VWindID = -1

  character(len=12) :: WindTransform

  type(struct_hco), pointer :: hco_bhi    ! Analysis horizontal grid parameters
  type(struct_hco), pointer :: hco_bstats ! Grid-point std dev horizontal grid parameters

  type(struct_lst)     :: lst_bhi    ! Spectral transform Parameters for B
  type(struct_lst)     :: lst_wind   ! Spectral transform Parameters for Vort/Div -> Psi/Chi

  real(8), allocatable :: bsqrt  (:,:,:)  ! B^1/2

  integer              :: nControlVariable
  integer              :: trunc
  integer              :: nksdim
  integer              :: nkgdim
  integer              :: cvDim
  integer              :: cvDim_mpiglobal
  
  integer              :: nlev_M
  integer              :: nlev_T

  logical              :: regrid
  logical              :: initialized = .false.

  integer              :: LatPerPE, LatPerPEmax, myLatBeg, myLatEnd
  integer              :: LonPerPE, LonPerPEmax, myLonBeg, myLonEnd

  integer,parameter    :: maxNumLevels=200
  real(8)              :: scaleFactor(maxNumLevels)

  character(len=8), parameter     :: BStatsFilename = './bgcov'

contains

  !--------------------------------------------------------------------------
  ! lbhi_SETUP
  !--------------------------------------------------------------------------
  subroutine lbhi_Setup( hco_anl_in, hco_core_in, vco_anl_in, cvDim_out )
    implicit none

    type(struct_hco), pointer, intent(in)    :: hco_anl_in
    type(struct_hco), pointer, intent(in)    :: hco_core_in
    type(struct_vco), pointer, intent(in)    :: vco_anl_in
    integer,          intent(out)   :: cvDim_out

    integer  :: var
    integer  :: ntrunc

    integer  :: iu_bstats = 0
    integer  :: iu_flnml = 0

    integer  :: ier, fnom, fstouv, fstfrm, fclos, levIndex, nLev

    logical  :: FileExist

    type(struct_vco), pointer :: vco_file => null()

    !namelist
    NAMELIST /NAMBHI/ntrunc,scaleFactor

    write(*,*)
    write(*,*) 'lbhi_Setup: Starting...'

    !
    !- 0.  Read namelist options
    !
    ntrunc         = 75     ! default values
    scaleFactor(:) =  0.0d0 ! default values

    ier = fnom(iu_flnml,'./flnml','FTN+SEQ+R/O',0)
    write(*,*)
    write(*,*) 'lbhi_setup: Reading namelist, ier = ',ier
    read(iu_flnml,nml=nambhi)
    write(*,nml=nambhi)
    ier = fclos(iu_flnml)

    nLev = max(max(vco_anl_in%nlev_M,vco_anl_in%nlev_T),1)

    do levIndex = 1, nLev
      if ( scaleFactor(levIndex) > 0.0d0 ) then 
        scaleFactor(levIndex) = sqrt(scaleFactor(levIndex))
      else
        scaleFactor(levIndex) = 0.0d0
      end if
    end do

    write(*,*) ' sum(scaleFactor) : ',sum(scaleFactor(1:nLev))
    if ( sum(scaleFactor(1:nLev)) == 0.0d0 ) then
      write(*,*) 'lambmatrixHI: scaleFactor=0, skipping rest of setup'
      cvDim_out   = 0
      return
    end if

    !- Setup the LAM analysis grid metrics
    call agd_SetupFromHCO( hco_anl_in, hco_core_in ) ! IN
    
    trunc = ntrunc
    write(*,*)
    write(*,*) 'Spectral TRUNCATION = ', trunc

    !
    !- 1.  Open the background stats file
    !
    inquire(file=trim(BStatsFilename), exist=FileExist)

    if ( FileExist ) then
      ier = fnom(iu_bstats,trim(BStatsFilename),'RND+OLD+R/O',0)
      if ( ier == 0 ) then
        write(*,*)
        write(*,*) 'Background Stats File :', trim(BStatsFilename)
        write(*,*) 'opened as unit file ',iu_bstats
        ier = fstouv(iu_bstats,'RND+OLD')
      else
        write(*,*)
        write(*,*) 'lbhi_Setup: Error in opening the background stats file'
        write(*,*) trim(BStatsFilename)
        call utl_abort('lbhi_Setup')
      end if
    else
      write(*,*)
      write(*,*) 'lbhi_Setup: The background stats file DOES NOT EXIST'
      write(*,*) trim(BStatsFilename)
      call utl_abort('lbhi_Setup')
    end if

    ! Check if analysisgrid and covariance file have the same vertical levels
    call vco_SetupFromFile( vco_file,      & ! OUT
                            BStatsFilename ) ! IN
    if (.not. vco_equal(vco_anl_in,vco_file)) then
      call utl_abort('lamBmatrixHI: vco from analysisgrid and cov file do not match')
    end if

    !
    !- 2.  Set some variables
    !
    hco_bhi => hco_anl_in
    nlev_M  = vco_anl_in%nlev_M
    nlev_T  = vco_anl_in%nlev_T

    !- Read variables and vertical grid info from the background stats file
    call lbhi_GetControlVariableInfo( iu_bstats ) ! IN
    call lbhi_GetHorizGridInfo()

    nkgdim = 0
    do var = 1, nControlVariable
      allocate( ControlVariable(var)%GpStdDev (1:hco_bhi%ni, 1:hco_bhi%nj, 1:ControlVariable(var)%nlev) )
      allocate( ControlVariable(var)%ip1 (1:ControlVariable(var)%nlev) )

      if (ControlVariable(var)%GridType == 'TH') then
         ControlVariable(var)%ip1(:) = vco_anl_in%ip1_T(:)
      elseif (ControlVariable(var)%GridType == 'MM') then
         ControlVariable(var)%ip1(:) = vco_anl_in%ip1_M(:)
      else
        ControlVariable(var)%ip1(:) = 0
      end if

      ControlVariable(var)%kDimStart = nkgdim + 1
      nkgdim = nkgdim + ControlVariable(var)%nlev
      ControlVariable(var)%kDimEnd    = nkgdim
    end do

    nksdim = nkgdim ! + nlev

    allocate( bsqrt  (1:nksdim, 1:nksdim ,0:trunc) )

    !- 2.2 Initialized the LAM spectral transform
    call mmpi_setup_lonbands(hco_bhi%ni,                  & ! IN
                               lonPerPE, lonPerPEmax, myLonBeg, myLonEnd ) ! OUT

    call mmpi_setup_latbands(hco_bhi%nj,                  & ! IN
                               latPerPE, latPerPEmax, myLatBeg, myLatEnd ) ! OUT

    call lst_Setup(lst_bhi,                         & ! OUT
                   hco_bhi%ni, hco_bhi%nj,          & ! IN
                   hco_bhi%dlon, trunc,             & ! IN
                   'LatLonMN', maxlevels_opt=nksdim)  ! IN

    cvDim     = nkgdim * lst_bhi%nla * lst_bhi%nphase
    cvDim_out = cvDim

    ! also compute mpiglobal control vector dimension
    call rpn_comm_allreduce(cvDim,cvDim_mpiglobal,1,"mpi_integer","mpi_sum","GRID",ier)

    !- 2.3 Initialized the Wind spectral transform
    if ( trim(WindTransform) == 'VortDiv' ) then
       call lst_Setup(lst_wind,                       & ! OUT
                      hco_bhi%ni, hco_bhi%nj,         & ! IN
                      hco_bhi%dlon, trunc,            & ! IN
                      'LatLonMN', maxlevels_opt=nlev_M) ! IN
    end if

    !
    !- 3.  Read info from the background error statistics file
    !
    call lbhi_ReadStats( iu_bstats )  ! IN

    !
    !- 4.  Close the background stats files
    !
    ier = fstfrm(iu_bstats)
    ier = fclos (iu_bstats)

    !
    !- 6.  Ending
    !
    initialized = .true.

    write(*,*)
    write(*,*) 'lbhi_Setup: Done!'

  end subroutine lbhi_Setup

  !--------------------------------------------------------------------------
  ! lbhi_GetControlVariableInfo
  !--------------------------------------------------------------------------
  subroutine lbhi_GetControlVariableInfo( iu_bstats )
    implicit none

    integer, intent(in) :: iu_bstats

    integer :: key, fstinf, fstlir, fstlir_s
    integer :: ni, nj, nlev
    integer :: dateo, nk
    integer :: ip1, ip2, ip3
    integer :: var

    character(len=4 )      :: nomvar
    character(len=2 )      :: typvar
    character(len=12)      :: etiket

    character(len=4), allocatable :: ControlModelVarnameList(:)
    character(len=4), allocatable :: ControlBhiVarnameList  (:)
    character(len=2), allocatable :: ControlVarGridTypeList (:)
    integer, allocatable          :: ControlVarNlevList     (:)

    !
    !- 1.  How Many Control Variables do we have?
    !
    dateo  = -1
    etiket = 'NLEV'
    ip1    = -1
    ip2    = -1
    ip3    = -1
    typvar = ' '
    nomvar = 'CVL'

    key = fstinf( iu_bstats,                                  & ! IN
                  ni, nj, nk,                                 & ! OUT
                  dateo, etiket, ip1, ip2, ip3, typvar, nomvar )! IN

    if (key < 0) then
      write(*,*)
      write(*,*) 'lbhi_GetControlVariableInfo: Unable to find variable =',nomvar
      call utl_abort('lbhi_GetControlVariableInfo')
    end if

    nControlVariable = ni
    write(*,*)
    write(*,*) 'Number of Control Variables found = ', nControlVariable

    allocate(ControlModelVarnameList(nControlVariable))
    allocate(ControlBhiVarnameList  (nControlVariable))
    allocate(ControlVarGridTypeList (nControlVariable))
    allocate(ControlVarNlevList     (nControlVariable))

    !
    !- 2. Read the info from the input file
    !
    nomvar = 'CVN'

    etiket = 'MODEL'
    key = fstlir_s(ControlModelVarnameList,                    & ! OUT 
                   iu_bstats,                                  & ! IN
                   ni, nj, nlev,                               & ! OUT
                   dateo, etiket, ip1, ip2, ip3, typvar,nomvar)  ! IN
    if (key < 0) then
      write(*,*)
      write(*,*) 'lbhi_GetControlVariableInfo: Cannot find variable ', nomvar 
      call utl_abort('lbhi_GetControlVariableInfo') 
    end if

    etiket = 'B_HI'
    key = fstlir_s(ControlBhiVarnameList,                      & ! OUT 
                   iu_bstats,                                  & ! IN
                   ni, nj, nlev,                               & ! OUT
                   dateo, etiket, ip1, ip2, ip3, typvar,nomvar)  ! IN
    if (key < 0) then
      write(*,*)
      write(*,*) 'lbhi_GetControlVariableInfo: Cannot find variable ', nomvar
      call utl_abort('lbhi_GetControlVariableInfo') 
    end if
    

    nomvar = 'CVL'

    etiket = 'NLEV'
    key = fstlir  (ControlVarNlevList,                         & ! OUT 
                   iu_bstats,                                  & ! IN
                   ni, nj, nlev,                               & ! OUT
                   dateo, etiket, ip1, ip2, ip3, typvar,nomvar)  ! IN
    if (key < 0) then
      write(*,*)
      write(*,*) 'lbhi_GetControlVariableInfo: Cannot find variable ', nomvar
      call utl_abort('lbhi_GetControlVariableInfo') 
    end if

    etiket = 'LEVTYPE'
    key = fstlir_s(ControlVarGridTypeList,                     & ! OUT 
                   iu_bstats,                                  & ! IN
                   ni, nj, nlev,                               & ! OUT
                   dateo, etiket, ip1, ip2, ip3, typvar,nomvar)  ! IN
    if (key < 0) then
      write(*,*)
      write(*,*) 'lbhi_GetControlVariableInfo: Cannot find variable ', nomvar
      call utl_abort('lbhi_GetControlVariableInfo') 
    end if

    !
    !- 3. Introduce the info in the ControlVariable structure
    !
    do var = 1, nControlVariable
       ControlVariable(var)%nomvar(cv_model)= trim(ControlModelVarnameList(var))
       ControlVariable(var)%nomvar(cv_bhi)  = trim(ControlBhiVarnameList(var))
       ControlVariable(var)%nlev            = ControlVarNlevList(var)
       ControlVariable(var)%GridType        = trim(ControlVarGridTypeList(var))

       if      (trim(ControlVariable(var)%nomvar(cv_model)) == 'UU' ) then 
          UWindID = var
       else if (trim(ControlVariable(var)%nomvar(cv_model)) == 'VV' ) then
          VWindID = var
       end if

       if ( trim(ControlVariable(var)%nomvar(cv_model)) == 'LQ' ) then
         ControlVariable(var)%nomvar(cv_model) = 'HU'  ! PATCH because gridStateVector uses HU
       end if

       write(*,*)
       write(*,*) 'nomvar(cv_model) = ', ControlVariable(var)%nomvar(cv_model)
       write(*,*) 'nomvar(cv_bhi)   = ', ControlVariable(var)%nomvar(cv_bhi)
       write(*,*) 'nlev             = ', ControlVariable(var)%nlev
       write(*,*) 'gridtype         = ', ControlVariable(var)%GridType
    end do

    deallocate(ControlModelVarnameList)
    deallocate(ControlBhiVarnameList)
    deallocate(ControlVarGridTypeList)
    deallocate(ControlVarNlevList)

    !
    !- 4. Error traps
    !

    !- Make sure that all the control variables are present in GridStateVector
    do var = 1, nControlVariable
       if ( .not. gsv_varExist(varName=ControlVariable(var)%nomvar(cv_model)) ) then
          write(*,*)
          write(*,*) 'lbhi_GetControlVariableInfo: The following variable is MISSING in GridStateVector'
          write(*,*) trim(ControlVariable(var)%nomvar(cv_model))
          call utl_abort('lbhi_GetControlVariableInfo')
       end if
    end do

    !
    !- 5.  Set the type of momentum control variables
    !
    if (UWindID /= -1 .and. VWindID /= -1) then
       if ( trim(ControlVariable(UWindID)%nomvar(cv_model)) == 'UU' .and. &
            trim(ControlVariable(UWindID)%nomvar(cv_bhi))   == 'PP' .and. &
            trim(ControlVariable(VWindID)%nomvar(cv_model)) == 'VV' .and. &
            trim(ControlVariable(VWindID)%nomvar(cv_bhi))   == 'CC' ) then
          write(*,*)
          write(*,*) 'Momentum Control Variables = Psi-Chi'
          WindTransform = 'PsiChi'
       else if ( trim(ControlVariable(UWindID)%nomvar(cv_model)) == 'UU' .and. &
            trim(ControlVariable(UWindID)%nomvar(cv_bhi))   == 'QR' .and. &
            trim(ControlVariable(VWindID)%nomvar(cv_model)) == 'VV' .and. &
            trim(ControlVariable(VWindID)%nomvar(cv_bhi))   == 'DD' ) then
          write(*,*)
          write(*,*) 'Momentum Control Variables = Vort-Div'
          WindTransform = 'VortDiv'
       else if ( trim(ControlVariable(UWindID)%nomvar(cv_model)) == 'UU' .and. &
            trim(ControlVariable(UWindID)%nomvar(cv_bhi))   == 'UU' .and. &
            trim(ControlVariable(VWindID)%nomvar(cv_model)) == 'VV' .and. &
            trim(ControlVariable(VWindID)%nomvar(cv_bhi))   == 'VV' ) then
          write(*,*)
          write(*,*) 'Momentum Control Variables = U-V '
          WindTransform = 'UV'
       else
          call utl_abort('lbhi_GetControlVariableInfo: Unkown Wind Tranform')
       end if
    end if

  end subroutine lbhi_GetControlVariableInfo

  !--------------------------------------------------------------------------
  ! lbhi_GetHorizGridInfo
  !--------------------------------------------------------------------------
  subroutine lbhi_GetHorizGridInfo()
    implicit none

    character(len=4)          :: varName

    !
    !- 1.  Get horizontal grid parameters
    !
    varName = ControlVariable(1)%nomvar(cv_bhi)

    call hco_setupFromFile(hco_bstats, BStatsFilename, ' ', 'bstats', &  ! IN
                           varName_opt=varName)                          ! IN

    !- 1.3 Regridding needed ?
    if ( hco_equal(hco_bstats,hco_bhi) ) then
      Regrid = .false.
      write(*,*)
      write(*,*) 'lbhi_GetHorizGridInfo: No Horizontal regridding needed'
    else
      Regrid = .true.
      write(*,*)
      write(*,*) 'lbhi_GetHorizGridInfo: Horizontal regridding is needed'
    end if

  end subroutine lbhi_GetHorizGridInfo

  !--------------------------------------------------------------------------
  ! lbhi_READSTATS
  !--------------------------------------------------------------------------
  subroutine lbhi_ReadStats( iu_bstats )
    implicit none

    integer, intent(in) :: iu_bstats

    !
    !- 1.  Read the background error statistics
    !

    !- 1.1 Verical correlations of control variables in spectral space
    call lbhi_ReadBSqrt( iu_bstats )        ! IN

    !- 1.2 Mass - Rotational wind statistical linear balance operator

    ! JFC: Pas encore code
    !if ( usePtoT ) then
    !  call lbhi_ReadPtoT( iu_bstats, PtoT_Type )
    !end if

    !- 1.3 Read grid-point standard deviations of control variables
    call lbhi_ReadGridPointStdDev(iu_bstats ) ! IN

  end subroutine lbhi_ReadStats

  !--------------------------------------------------------------------------
  ! lbhi_ReadBSqrt
  !--------------------------------------------------------------------------
  subroutine lbhi_ReadBSqrt( iu_bstats )
    implicit none

    integer, intent(in) :: iu_bstats

    real(8), allocatable :: bsqrt2d  (:,:)

    integer :: key, fstinf, fstinl, totwvnb, infon
    integer, parameter :: nmax=2000
    integer :: liste(nmax)

    integer                     :: ip1, ip2, ip3
    integer                     :: ni_t, nj_t, nlev_t, dateo  
    character(len=4 )           :: nomvar
    character(len=2 )           :: typvar
    character(len=12)           :: etiket

    dateo  = -1
    etiket = 'B_SQUAREROOT'
    ip1    = -1
    ip3    = -1
    typvar = ' '
    nomvar = 'ZN'

    !
    !- 1.  Find the truncation in the stats file
    !
    ip2    = -1

    key = fstinl(iu_bstats,                                   & ! IN
                ni_t, nj_t, nlev_t,                           & ! OUT
                dateo, etiket, ip1, ip2, ip3, typvar, nomvar, & ! IN
                liste, infon,                                 & ! OUT
                nmax )                                          ! IN

    if (key >= 0) then
      !- 1.2 Ensure spectral trunctation are the same
      if ( infon - 1  /= trunc ) then
        write(*,*)
        write(*,*) 'lbhi_ReadBSqrt: Truncation here and on stats file different'
        write(*,*) 'VAR truncation        = ', trunc
        write(*,*) 'Stats file truncation = ', infon-1
        call utl_abort('lbhi_ReadBSqrt')
      end if
    else
      write(*,*)
      write(*,*) 'lbhi_ReadBSqrt: Cannot find B square-root ', nomvar 
      call utl_abort('lbhi_ReadBSqrt')
    end if

    !
    !- 2.   Read B^0.5
    !
    allocate( bsqrt2d  (1:nksdim, 1:nksdim) )

    do totwvnb = 0, trunc

      ip2    = totwvnb

      !- 2.1 Check if field exists and its dimensions
      key = fstinf( iu_bstats,                                  & ! IN
                    ni_t, nj_t, nlev_t,                         & ! OUT
                    dateo, etiket, ip1, ip2, ip3, typvar,nomvar)  ! IN

      if (key >= 0) then
        !- 2.2 Ensure that the number of vertical levels are compatible
        if ( ni_t /= nksdim .or. nj_t /= nksdim  ) then
          write(*,*)
          write(*,*) 'lbhi_ReadBSqrt: BG stat levels inconsitencies'
          write(*,*) 'for BSQRT: ni_t, nj_t, nksdim =', ni_t, nj_t, nksdim
          call utl_abort('lbhi_ReadBSqrt')
        endif

        !- 2.3 Reading
        key = utl_fstlir( bsqrt2d,                                    & ! OUT 
                          iu_bstats,                                  & ! IN
                          ni_t, nj_t, nlev_t,                         & ! OUT
                          dateo, etiket, ip1, ip2, ip3, typvar,nomvar)  ! IN
      else
        write(*,*)
        write(*,*) 'lbhi_ReadBSqrt: Cannot find BSQRT for totwvnb = ', totwvnb
        call utl_abort('lbhi_ReadBSqrt')
      end if

      !- 2.4 Transfer to a 3D array
      bsqrt(:,:,totwvnb) = bsqrt2d(:,:)

    end do

    deallocate( bsqrt2d )

  end subroutine lbhi_ReadBSqrt

  !--------------------------------------------------------------------------
  ! lbhi_ReadGridPointStdDev
  !--------------------------------------------------------------------------
  subroutine lbhi_ReadGridPointStdDev(iu_bstats)
    implicit none

    integer, intent(in) :: iu_bstats

    real(8), allocatable :: StdDev2D(:,:)
    real(8), allocatable :: StdDev2D_Regrid(:,:)

    integer :: ezdefset, ier
    integer :: ni_t, nj_t, nlev_t, var, k
    integer :: dateo, ip1,ip2,ip3

    character(len=4 )      :: nomvar
    character(len=2 )      :: typvar
    character(len=12)      :: etiket

    real(8) :: UnitConv

    !
    !- 1.  Read grid point standard deviations
    !
    allocate( StdDev2D(1:hco_bstats%ni,1:hco_bstats%nj) )
    if (Regrid) then
      allocate( StdDev2D_Regrid(1:hco_bhi%ni, 1:hco_bhi%nj) )
      ier = ezdefset(hco_bhi%EZscintID, hco_bstats%EZscintID)             ! IN
    end if

    !- 1.1 Loop over Control Variables
    do var = 1, nControlVariable

      !- 1.2 Loop over vertical Levels
      do k = 1, ControlVariable(var)%nlev
        dateo  = -1
        ip1    = ControlVariable(var)%ip1(k)
        ip2    = -1
        ip3    = -1
        typvar = ' '
        nomvar = trim(ControlVariable(var)%nomvar(cv_bhi))
        etiket = 'STDDEV'
        if ( trim(nomvar) == 'P0') then
          UnitConv = 100.0d0 ! hPa -> Pa
        else
          UnitConv = 1.0d0
        end if

        !- 1.2.1 Reading
        ier = utl_fstlir( StdDev2D,                                   & ! OUT 
                          iu_bstats,                                  & ! IN
                          ni_t, nj_t, nlev_t,                         & ! OUT
                          dateo, etiket, ip1, ip2, ip3, typvar,nomvar)  ! IN 

        if (ier < 0) then
          write(*,*)
          write(*,*) 'lbhi_ReadGridPointStdDev: Cannot find Std Deviations'
          write(*,*) 'nomvar =', trim(ControlVariable(var)%nomvar(cv_bhi))
          write(*,*) 'etiket =', trim(etiket)
          write(*,*) 'ip1    =', ControlVariable(var)%ip1(k)
          call utl_abort('lbhi_ReadGridPointStdDev')
        end if

        if (ni_t /= hco_bstats%ni .or. nj_t /= hco_bstats%nj) then
          write(*,*)
          write(*,*) 'lbhi_ReadGridPointStdDev: Invalid dimensions for ...'
          write(*,*) 'nomvar      =', trim(ControlVariable(var)%nomvar(cv_bhi))
          write(*,*) 'etiket      =', trim(etiket)
          write(*,*) 'ip1         =', ControlVariable(var)%ip1(k)
          write(*,*) 'Found ni,nj =', ni_t, nj_t 
          write(*,*) 'Should be   =', hco_bstats%ni, hco_bstats%nj
          call utl_abort('lbhi_ReadGridPointStdDev')
        end if

        !- 1.2.2 Regrid (if necessary) and transfer to 3D array
        if ( .not. Regrid) then
           ControlVariable(var)%GpStdDev(:,:,k) = StdDev2D(:,:)
        else
           ! Note: EZSCINT setup was done above
           ier = utl_ezsint(StdDev2D_Regrid, StdDev2D, interpDegree='LINEAR')
           ControlVariable(var)%GpStdDev(:,:,k) = StdDev2D_Regrid(:,:)
        end if

        !- 1.3 Scaling
        ControlVariable(var)%GpStdDev(:,:,k) = ControlVariable(var)%GpStdDev(:,:,k) * &
                                                 UnitConv * scaleFactor(k)

      end do

    end do

    deallocate( StdDev2D )
    if (Regrid) then
      deallocate( StdDev2D_Regrid )
    end if

  end subroutine lbhi_ReadGridPointStdDev

  !--------------------------------------------------------------------------
  ! lbhi_bSqrt
  !--------------------------------------------------------------------------
  subroutine lbhi_bSqrt(controlVector_in, statevector, stateVectorRef_opt)
    implicit none

    ! Arguments
    real(8),          intent(in)           :: controlVector_in(cvDim)
    type(struct_gsv), intent(inout)        :: statevector
    type(struct_gsv), intent(in), optional :: statevectorRef_opt

    ! Locals
    real(8), allocatable :: gd_out(:,:,:)
    real(8), allocatable :: hiControlVector(:,:,:)

    if ( .not. initialized ) then
      if(mpi_myid == 0) write(*,*) 'lbhi_bSqrt: LAM_bMatrixHI not initialized'
      return
    endif

    write(*,*)
    write(*,*) 'lbhi_bSqrt: Starting ...'

    !
    !-  1.  Extract data from the 1D controlVector array
    !
    allocate( hiControlVector(lst_bhi%nla, lst_bhi%nphase, nksdim) )

    call lbhi_cain( controlVector_in,  & ! IN
                    hiControlVector )    ! OUT

    !
    !-  2.  Move from control variables space to model variables space
    !
    allocate( gd_out  (myLonBeg:myLonEnd, myLatBeg:myLatEnd, nksdim) )

    call lbhi_cv2gd( hiControlVector,   & ! IN
                     gd_out           )   ! OUT
    
    deallocate(hiControlVector)

    !
    !-  3.  Transfer results to statevector structure
    !
    call StatevectorInterface( statevector,   & ! INOUT
                               gd_out,        & ! IN
                              'ToStateVector' ) ! IN
    deallocate(gd_out)

    !
    !-  4. Convert LQ_inc to HU_inc
    !
    if ( gsv_varExist(varName='HU') ) then
      call gvt_transform( statevector,   &                        ! INOUT
                          'LQtoHU_tlm',  &                        ! IN
                          stateVectorRef_opt=stateVectorRef_opt ) ! IN
    end if

    write(*,*)
    write(*,*) 'lbhi_bSqrt: Done'

  end subroutine lbhi_bSqrt

  !--------------------------------------------------------------------------
  ! lbhi_bSqrtAdj
  !--------------------------------------------------------------------------
  subroutine lbhi_bSqrtAdj(statevector, controlVector_out, stateVectorRef_opt)
    implicit none

    ! Arguments
    real(8),          intent(out)          :: controlVector_out(cvDim)
    type(struct_gsv), intent(inout)        :: statevector
    type(struct_gsv), intent(in), optional :: statevectorRef_opt

    ! Locals
    real(8), allocatable :: gd_in(:,:,:)
    real(8), allocatable :: hiControlVector(:,:,:)

    if ( .not. initialized ) then
      if(mpi_myid == 0) write(*,*) 'lbhi_bSqrtAdj: LAM_bMatrixHI not initialized'
      return
    endif

    write(*,*)
    write(*,*) 'lbhi_bSqrtAdj: Starting ...'

    !
    !-  4. Convert LQ_inc to HU_inc
    !
    if ( gsv_varExist(varName='HU') ) then
      call gvt_transform( statevector,  &                         ! INOUT
                          'LQtoHU_ad',  &                         ! IN
                          stateVectorRef_opt=stateVectorRef_opt ) ! IN
    end if

    !
    !-  3.  Extract data from the StateVector
    !
    allocate( gd_in(myLonBeg:myLonEnd, myLatBeg:myLatEnd, nksdim) )

    call StatevectorInterface ( statevector,      & ! IN
                                gd_in,            & ! OUT
                               'FromStateVector' )  ! IN

    !
    !-  2.  Move from model variables space to control variables space
    !
    allocate( hiControlVector(lst_bhi%nla, lst_bhi%nphase, nksdim) )
    hiControlVector(:,:,:) = 0.d0

    call lbhi_cv2gdAdj( hiControlVector, & ! OUT
                        gd_in          )   ! IN

    !
    !-  1.  Put data into the 1D controlVector array
    !
    controlVector_out(:) = 0.d0
    call lbhi_cainAdj(controlVector_out, hiControlVector)

    deallocate(gd_in)
    deallocate(hiControlVector)

    write(*,*)
    write(*,*) 'lbhi_bSqrtAdj: Done'

  end subroutine lbhi_bSqrtAdj

  !--------------------------------------------------------------------------
  ! lbhi_cv2gd
  !--------------------------------------------------------------------------
  subroutine lbhi_cv2gd(hiControlVector_in, gd_out)
    implicit none

    real(8), intent(inout) :: hiControlVector_in(lst_bhi%nla, lst_bhi%nphase, nksdim)
    real(8), intent(out)   :: gd_out(myLonBeg:myLonEnd  ,myLatBeg:myLatEnd  ,1:nksdim)

    real(8), allocatable :: uphy(:,:,:)
    real(8), allocatable :: vphy(:,:,:)
    real(8), allocatable :: psi(:,:,:)
    real(8), allocatable :: chi(:,:,:)

    integer :: kstart, kend, var

    character(len=19)   :: kind

    !
    !- 1. B^1/2 * xi (in spectral space)
    !
    call lbhi_bSqrtXi(hiControlVector_in)    ! INOUT

    !
    !- 2. Spectral Space -> Grid Point Space
    !
    kind = 'SpectralToGridPoint'
    call lst_VarTransform(lst_bhi,            & ! IN
                          hiControlVector_in, & ! IN
                          gd_out,             & ! OUT
                          kind, nksdim )        ! IN

    !
    !- 3.  Multiply by the grid point standard deviations
    !
    !$OMP PARALLEL DO PRIVATE(var,kstart,kend)
    do var = 1, nControlVariable
      kstart = ControlVariable(var)%kDimStart
      kend   = ControlVariable(var)%kDimEnd
      gd_out(:,:,kstart:kend) = gd_out(:,:,kstart:kend) * ControlVariable(var)%GpStdDev(myLonBeg:myLonEnd,myLatBeg:myLatEnd,:)
    end do
    !$OMP END PARALLEL DO

    !
    !- 4.  Momentum variables transform
    !
    if ( UWindID /= -1 .and. VWindID /= -1) then
       if ( ControlVariable(UWindID)%nlev /= nlev_M .or. &
            ControlVariable(VWindID)%nlev /= nlev_M  ) then
          call utl_abort('lbhi_cv2gd: Error in Wind related parameters')
       end if
    
       if ( trim(WindTransform) /= 'UV') then

          allocate(uphy(myLonBeg:myLonEnd,myLatBeg:myLatEnd,1:nlev_M))
          allocate(vphy(myLonBeg:myLonEnd,myLatBeg:myLatEnd,1:nlev_M))
          allocate(psi (myLonBeg:myLonEnd,myLatBeg:myLatEnd,1:nlev_M))
          allocate(chi (myLonBeg:myLonEnd,myLatBeg:myLatEnd,1:nlev_M))

          !- 4.1 Vort / Div -> Psi / Chi
          if ( trim(WindTransform) == 'VortDiv') then

             !  4.1.1 Putting vorticity in psi array and divergence in chi array
             psi(:,:,:) = gd_out(:,:,ControlVariable(UWindID)%kDimStart:ControlVariable(UWindID)%kDimEnd)
             chi(:,:,:) = gd_out(:,:,ControlVariable(VWindID)%kDimStart:ControlVariable(VWindID)%kDimEnd)
             
             !  4.1.2 Vort -> Psi
             call lst_Laplacian(lst_wind,           & ! IN
                                Psi,                & ! INOUT
                                'Inverse', nlev_M)    ! IN    

             !  4.1.3 Div -> Chi
             call lst_Laplacian(lst_wind,          & ! IN
                                Chi,               & ! INOUT
                                'Inverse', nlev_M)   ! IN

          end if

          !- 4.2 Psi / Chi -> U-wind / V-wind
          if ( trim(WindTransform) == 'PsiChi') then
             psi(:,:,:) = gd_out(:,:,ControlVariable(UWindID)%kDimStart:ControlVariable(UWindID)%kDimEnd)
             chi(:,:,:) = gd_out(:,:,ControlVariable(VWindID)%kDimStart:ControlVariable(VWindID)%kDimEnd)
          end if

          !  4.2.1 Do Transform
          call agd_PsiChiToUV(psi, chi,           & ! IN
                              uphy, vphy,         & ! OUT
                              nlev_M)               ! IN

          !  4.2.2 Insert results in gd_out and deallocate memories
          gd_out(:,:,1       :  nlev_M) = uphy(:,:,:)
          gd_out(:,:,nlev_M+1:2*nlev_M) = vphy(:,:,:)
          
          deallocate(chi)
          deallocate(psi)
          deallocate(vphy)
          deallocate(uphy)
          
       end if

    end if

  end subroutine lbhi_cv2gd

  !--------------------------------------------------------------------------
  ! lbhi_cv2gdAdj
  !--------------------------------------------------------------------------
  subroutine lbhi_cv2gdAdj(hiControlVector_out, gd_in)
    implicit none

    real(8), intent(out)   :: hiControlVector_out(lst_bhi%nla, lst_bhi%nphase, nksdim)
    real(8), intent(inout) :: gd_in(myLonBeg:myLonEnd, myLatBeg:myLatEnd ,1:nksdim)

    real(8), allocatable :: uphy(:,:,:)
    real(8), allocatable :: vphy(:,:,:)
    real(8), allocatable :: psi(:,:,:)
    real(8), allocatable :: chi(:,:,:)

    integer :: kstart, kend, var

    character(len=19)   :: kind

    !
    !- 4.  Momentum variables transform
    !
    if ( UWindID /= -1 .and. VWindID /= -1) then
       if ( ControlVariable(UWindID)%nlev /= nlev_M .or. &
            ControlVariable(VWindID)%nlev /= nlev_M  ) then
          call utl_abort('lbhi_cv2gdadj: Error in Wind related parameters')
       end if
       
       if ( trim(WindTransform) /= 'UV' ) then

          allocate(uphy(myLonBeg:myLonEnd,myLatBeg:myLatEnd,1:nlev_M))
          allocate(vphy(myLonBeg:myLonEnd,myLatBeg:myLatEnd,1:nlev_M))
          allocate(psi (myLonBeg:myLonEnd,myLatBeg:myLatEnd,1:nlev_M))
          allocate(chi (myLonBeg:myLonEnd,myLatBeg:myLatEnd,1:nlev_M))
          
          !- 4.2  U-wind / V-wind -> Psi / Chi

          !  4.2.2 Extract winds
          uphy(:,:,:) = gd_in(:,:,1       :  nlev_M)
          vphy(:,:,:) = gd_in(:,:,nlev_M+1:2*nlev_M)

          !  4.2.1 Do Transform
          call agd_PsiChiToUVAdj( psi, chi,           & ! OUT
                                  uphy, vphy,         & ! IN
                                  nlev_M)               ! IN

          !- 4.1 Psi / Chi -> Vort / Div (Laplacian and inverse Laplacian are auto-adjoint operators)
          if ( trim(WindTransform) == 'VortDiv' ) then

             !  4.1.2 Psi -> Vort
             call lst_Laplacian(lst_wind,          & ! IN
                                Psi,               & ! INOUT
                                'Inverse', nlev_M)   ! IN    

             !  4.1.3 Chi -> Vort
             call lst_Laplacian(lst_wind,          & ! IN
                                Chi,               & ! INOUT
                                'Inverse', nlev_M)   ! IN

          end if

          !- Insert results in gd and deallocate moemories
          gd_in(:,:,ControlVariable(UWindID)%kDimStart:ControlVariable(UWindID)%kDimEnd) = psi(:,:,:)
          gd_in(:,:,ControlVariable(VWindID)%kDimStart:ControlVariable(VWindID)%kDimEnd) = chi(:,:,:)

          deallocate(chi)
          deallocate(psi)
          deallocate(vphy)
          deallocate(uphy)

       end if

    end if

    !
    !- 3.  Multiply by the grid point standard deviations
    !

    !$OMP PARALLEL DO PRIVATE(var,kstart,kend)
    do var = 1, nControlVariable
      kstart = ControlVariable(var)%kDimStart
      kend   = ControlVariable(var)%kDimEnd
      gd_in(:,:,kstart:kend) = gd_in(:,:,kstart:kend) * ControlVariable(var)%GpStdDev(myLonBeg:myLonEnd,myLatBeg:myLatEnd,:)
    end do
    !$OMP END PARALLEL DO

    !
    !- 2. Grid Point Space -> Spectral Space
    !
    kind = 'GridPointToSpectral'
    call lst_VarTransform(lst_bhi,                & ! IN
                          hiControlVector_out,    & ! OUT
                          gd_in,                  & ! IN
                          kind, nksdim )            ! IN

    !
    !- 1. B^1/2 * xi (in spectral space)
    !
    call lbhi_bSqrtXi( hiControlVector_out )    ! INOUT

  end subroutine lbhi_cv2gdAdj

  !--------------------------------------------------------------------------
  ! lbhi_bSqrtXi
  !--------------------------------------------------------------------------
  subroutine lbhi_bSqrtXi(hiControlVector_in)
    implicit none

    real(8), intent(inout) :: hiControlVector_in(lst_bhi%nla, lst_bhi%nphase, nksdim)

    real(8), allocatable :: sp_in (:,:,:)
    real(8), allocatable :: sp_out(:,:,:)

    integer :: totwvnb, e, k, ila
    integer :: m, n, lda, ldb, ldc

    !
    !- 1. B^1/2 * xi (in spectral space)
    !
    do totwvnb = 0, trunc

      if ( lst_bhi%nePerK(totwvnb) == 0 ) then
         cycle
      end if

      allocate( sp_in (nksdim,lst_bhi%nphase,lst_bhi%nePerK(totwvnb)) )
      allocate( sp_out(nksdim,lst_bhi%nphase,lst_bhi%nePerK(totwvnb)) )

      !- 1.1 Select spectral elements associated with the total wavenumber
      !$OMP PARALLEL DO PRIVATE(e,ila,k)
      do e = 1, lst_bhi%nePerK(totwvnb)
        ila = lst_bhi%ilaFromEK(e,totwvnb)
        do k = 1, nksdim
          sp_in(k,1:lst_bhi%nphase,e) = hiControlVector_in(ila,1:lst_bhi%nphase,k)
        end do
      end do
      !$OMP END PARALLEL DO

      !- 1.2 Compute bsqrt * sp_in using DGEMM 

      ! For documentation on dgemm, see: http://www.netlib.org/blas/dgemm.f
      ! Matrix A = BSQRT(:,:,totwvnb)
      ! Matrix B = SP_IN
      ! Matrix C = SP_OUT
      m   = nksdim
      n   = lst_bhi%nphase * lst_bhi%nePerK(totwvnb)
      k   = nksdim
      lda = nksdim
      ldb = nksdim
      ldc = nksdim

      call dgemm( 'N', 'N', m, n, k, 1.d0,                   &  ! IN
                  bsqrt(:,:,totwvnb), lda, sp_in, ldb, 0.d0, &  ! IN
                  sp_out,                                    &  ! OUT
                  ldc )                                         ! IN

      !- 1.3 Replace sp values with output matrix
      !$OMP PARALLEL DO PRIVATE(e,ila,k)
      do e = 1, lst_bhi%nePerK(totwvnb)
        ila = lst_bhi%ilaFromEK(e,totwvnb)
        do k = 1, nksdim
          hiControlVector_in(ila,1:lst_bhi%nphase,k) = sp_out(k,1:lst_bhi%nphase,e)
        end do
      end do
      !$OMP END PARALLEL DO

      deallocate(sp_in)
      deallocate(sp_out)

    end do ! Total Wavenumber

  end subroutine lbhi_bSqrtXi

  !--------------------------------------------------------------------------
  ! lbhi_cain
  !--------------------------------------------------------------------------
  SUBROUTINE lbhi_cain(controlVector_in, hiControlVector_out)
    implicit none

    real(8), intent(in)    :: controlVector_in(cvDim)
    real(8), intent(out)   :: hiControlVector_out(lst_bhi%nla,lst_bhi%nphase,nksdim)

    integer :: dim, k, ila, p

    dim = 0
    hiControlVector_out(:,:,:) = 0.0d0
    do k = 1, nksdim
      do ila = 1, lst_bhi%nla
        do p = 1, lst_bhi%nphase
          dim = dim + 1
          hiControlVector_out(ila,p,k) = controlVector_in(dim) * lst_bhi%NormFactor(ila,p)
        end do
      end do
    end do

  end SUBROUTINE lbhi_cain

  !--------------------------------------------------------------------------
  ! lbhi_cainAdj
  !--------------------------------------------------------------------------
  SUBROUTINE lbhi_cainAdj(controlVector_out, hiControlVector_in)
    implicit none

    real(8), intent(out)   :: controlVector_out(cvDim)
    real(8), intent(in )   :: hiControlVector_in(lst_bhi%nla,lst_bhi%nphase,nksdim)

    integer :: dim, k, ila, p

    dim = 0
    do k = 1, nksdim
      do ila = 1, lst_bhi%nla
        do p = 1, lst_bhi%nphase
          dim = dim + 1
          controlVector_out(dim) = controlVector_out(dim) + &
                                   hiControlVector_in(ila,p,k) * lst_bhi%NormFactorAd(ila,p)
        end do
      end do
    end do

  end SUBROUTINE lbhi_cainAdj

  !--------------------------------------------------------------------------
  ! StatevectorInterface
  !--------------------------------------------------------------------------
  subroutine StatevectorInterface(statevector, gd, Direction)
    implicit none

    type(struct_gsv), intent(inout) :: statevector
    real(8),          intent(inout) :: gd(myLonBeg:myLonEnd,myLatBeg:myLatEnd,nksdim)
    character(len=*), intent(in)    :: Direction

    integer :: var
    integer :: kgdStart, kgdEnd, i, j, k, kgd, nlev

    real(4), pointer :: field_r4(:,:,:)
    real(8), pointer :: field_r8(:,:,:)

    character(len=4 )      :: varname

    logical :: ToStateVector

    select case ( trim(Direction) )
    case ('ToStateVector')
      ToStateVector = .true.
    case ('FromStateVector')
      ToStateVector = .false.
    case default
      write(*,*)
      write(*,*) 'StatevectorInterface: Unknown Direction ', trim(Direction)
      call utl_abort('StatevectorInterface')
    end select

    do var = 1, nControlVariable

      varname = ControlVariable(var)%nomvar(cv_model)

      if (.not. gsv_varExist(statevector,varname) ) then
         write(*,*)
         write(*,*) 'StatevectorInterface: The following variable is MISSING in GridStateVector'
         write(*,*) varname
         call utl_abort('StatevectorInterface')
      end if

      if (gsv_getDataKind(statevector) == 4) then
        call gsv_getField(statevector,field_r4,varname)
      else
        call gsv_getField(statevector,field_r8,varname)
      end if

      kgdStart = ControlVariable(var)%kDimStart
      kgdEnd   = ControlVariable(var)%kDimEnd
   
      nlev = gsv_getNumLev(statevector,vnl_varLevelFromVarname(varname))
      if ( kgdEnd - kgdStart + 1  /= nlev ) then
         write(*,*)
         write(*,*) 'StatevectorInterface: Number of vertical level mismatch'
         write(*,*) kgdEnd - kgdStart + 1, nlev
         call utl_abort('StatevectorInterface')
      end if

      if ( ToStateVector ) then
        if (gsv_getDataKind(statevector) == 4) then
          !$OMP PARALLEL DO PRIVATE(j,kgd,k,i)
          do kgd = kgdStart, kgdEnd
            k = kgd - kgdStart + 1
            do j = myLatBeg, myLatEnd
              do i = myLonBeg, myLonEnd
                field_r4(i,j,k) = gd(i,j,kgd)
              end do
            end do
          end do
          !$OMP END PARALLEL DO
        else
          !$OMP PARALLEL DO PRIVATE(j,kgd,k,i)
          do kgd = kgdStart, kgdEnd
            k = kgd - kgdStart + 1
            do j = myLatBeg, myLatEnd
              do i = myLonBeg, myLonEnd
                field_r8(i,j,k) = gd(i,j,kgd)
              end do
            end do
          end do
          !$OMP END PARALLEL DO
        end if
      else
        if (gsv_getDataKind(statevector) == 4) then
          !$OMP PARALLEL DO PRIVATE(j,kgd,k,i)
          do kgd = kgdStart, kgdEnd
            k = kgd - kgdStart + 1
            do j = myLatBeg, myLatEnd
              do i = myLonBeg, myLonEnd
                gd(i,j,kgd)  = field_r4(i,j,k)
              end do
            end do
          end do
          !$OMP END PARALLEL DO
        else
          !$OMP PARALLEL DO PRIVATE(j,kgd,k,i)
          do kgd = kgdStart, kgdEnd
            k = kgd - kgdStart + 1
            do j = myLatBeg, myLatEnd
              do i = myLonBeg, myLonEnd
                gd(i,j,kgd)  = field_r8(i,j,k)
              end do
            end do
          end do
          !$OMP END PARALLEL DO
        end if
      end if
   end do

  end subroutine StatevectorInterface

  !--------------------------------------------------------------------------
  ! lbhi_reduceToMPILocal
  !--------------------------------------------------------------------------
  SUBROUTINE lbhi_reduceToMPILocal(cv_mpilocal,cv_mpiglobal)
    implicit none
    real(8), intent(out) :: cv_mpilocal(cvDim)
    real(8), intent(in)  :: cv_mpiglobal(:)

    real(8), allocatable :: cv_allmaxmpilocal(:,:)

    integer, allocatable :: cvDim_allMpilocal(:), displs(:)    
    integer, allocatable :: ilaGlobal(:), allnlaLocal(:)
    integer, allocatable :: allilaGlobal(:,:)

    integer :: k, ila, p, ila_mpiglobal, jdim_mpilocal, jdim_mpiglobal
    integer :: ier, nlaMax, cvDim_maxmpilocal, jproc

    call rpn_comm_allreduce(cvDim, cvDim_maxmpilocal, &
         1,"MPI_INTEGER","MPI_MAX","GRID",ier)

    allocate(cvDim_allMpiLocal(mpi_nprocs))

    call rpn_comm_allgather(cvDim   ,1,"mpi_integer",       &
                            cvDim_allMpiLocal,1,"mpi_integer","GRID",ier)

    call rpn_comm_allreduce(lst_bhi%nla,nlaMax,1,"mpi_integer","mpi_max","GRID",ier)

    if (mpi_myid == 0) then
       allocate(allnlaLocal(mpi_nprocs))
       allocate(allilaGlobal(nlaMax,mpi_nprocs))
    else
       allocate(allnlaLocal(1))
       allocate(allilaGlobal(1,1))
    end if
    
    allocate(ilaGlobal(nlaMax))
    ilaGlobal(:)             = -1
    ilaGlobal(1:lst_bhi%nla) = lst_bhi%ilaGlobal(:)
    
    call rpn_comm_gather(lst_bhi%nla, 1, "mpi_integer",       &
                         allnlaLocal, 1, "mpi_integer", 0, "GRID", ier)
    call rpn_comm_gather(ilaGlobal   , nlaMax, "mpi_integer",       &
                         allilaGlobal, nlaMax, "mpi_integer",0 ,"GRID", ier)

    deallocate(ilaGlobal)

    ! assign part of mpiglobal vector from current mpi process
    if (mpi_myid == 0) then

       allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mpi_nprocs))

       !$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,k,ila,p,ila_mpiglobal,jdim_mpiglobal)
       do jproc = 0, (mpi_nprocs-1)
          cv_allmaxmpilocal(:,jproc+1) = 0.d0
          
          do k = 1, nksdim
             do ila = 1, allnlaLocal(jproc+1)
                do p = 1, lst_bhi%nphase

                   jdim_mpilocal = ( (k-1) * allnlaLocal(jproc+1) * lst_bhi%nphase ) + &
                                                        ( (ila-1) * lst_bhi%nphase ) + p

                   ila_mpiglobal = allilaGlobal(ila,jproc+1)
                   if ( ila_mpiglobal <= 0 ) then 
                      write(*,*) 'lbhi_reduceToMPILocal: invalid ila_mpiglobal index ', ila_mpiglobal
                      call utl_abort('lbhi_reduceToMPILocal')
                   end if

                   jdim_mpiglobal = ( (k-1) * lst_bhi%nlaGlobal * lst_bhi%nphase ) + &
                                            ( (ila_mpiglobal-1) * lst_bhi%nphase ) + p
  
                   if (jdim_mpilocal > cvDim_allMpiLocal(jproc+1)) then
                      write(*,*)
                      write(*,*) 'ERROR: jdim_mpilocal > cvDim_allMpiLocal(jproc+1)', jdim_mpilocal, cvDim_allMpiLocal(jproc+1) 
                      write(*,*) '       proc, k, ila, p = ',jproc,k,ila,p
                      call utl_abort('lbhi_reduceToMPILocal')
                   end if
                   if (jdim_mpiglobal > cvDim_mpiglobal) then
                      write(*,*)
                      write(*,*) 'ERROR: jdim_mpiglobal > cvDim_mpiglobal', jdim_mpiglobal, cvDim_mpiglobal
                      write(*,*) '       proc, k, ila, p = ',jproc,k,ila,p
                      call utl_abort('lbhi_reduceToMPILocal')
                   end if
                   
                   cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal)
                   
                end do
             end do
          end do
       end do
       !$OMP END PARALLEL DO

    else
       allocate(cv_allmaxmpilocal(1,1))
    end if

    deallocate(allnlaLocal)
    deallocate(allilaGlobal)

    !- Distribute
    allocate(displs(mpi_nprocs))
    !$OMP PARALLEL DO PRIVATE(jproc)
    do jproc = 0, (mpi_nprocs-1)
       displs(jproc+1) = jproc*cvDim_maxMpiLocal ! displacement wrt cv_allMaxMpiLocal from which
                                                 ! to take the outgoing data to process jproc
    end do
    !$OMP END PARALLEL DO

    call rpn_comm_scatterv(cv_allMaxMpiLocal, cvDim_allMpiLocal, displs, "mpi_double_precision", &
                           cv_mpiLocal      , cvDim , "mpi_double_precision", &
                           0, "GRID", ier)

   deallocate(displs) 
   deallocate(cv_allMaxMpiLocal)
   deallocate(cvDim_allMpiLocal)


  END SUBROUTINE lbhi_reduceToMPILocal

  !--------------------------------------------------------------------------
  ! lbhi_reduceToMPILocal_r4
  !--------------------------------------------------------------------------
  SUBROUTINE lbhi_reduceToMPILocal_r4(cv_mpilocal,cv_mpiglobal)
    implicit none
    real(4), intent(out) :: cv_mpilocal(cvDim)
    real(4), intent(in)  :: cv_mpiglobal(:)

    real(4), allocatable :: cv_allmaxmpilocal(:,:)

    integer, allocatable :: cvDim_allMpilocal(:), displs(:)    
    integer, allocatable :: ilaGlobal(:), allnlaLocal(:)
    integer, allocatable :: allilaGlobal(:,:)

    integer :: k, ila, p, ila_mpiglobal, jdim_mpilocal, jdim_mpiglobal
    integer :: ier, nlaMax, cvDim_maxmpilocal, jproc

    call rpn_comm_allreduce(cvDim, cvDim_maxmpilocal, &
         1,"MPI_INTEGER","MPI_MAX","GRID",ier)

    allocate(cvDim_allMpiLocal(mpi_nprocs))

    call rpn_comm_allgather(cvDim   ,1,"mpi_integer",       &
                            cvDim_allMpiLocal,1,"mpi_integer","GRID",ier)

    call rpn_comm_allreduce(lst_bhi%nla,nlaMax,1,"mpi_integer","mpi_max","GRID",ier)

    if (mpi_myid == 0) then
       allocate(allnlaLocal(mpi_nprocs))
       allocate(allilaGlobal(nlaMax,mpi_nprocs))
    else
       allocate(allnlaLocal(1))
       allocate(allilaGlobal(1,1))
    end if
    
    allocate(ilaGlobal(nlaMax))
    ilaGlobal(:)             = -1
    ilaGlobal(1:lst_bhi%nla) = lst_bhi%ilaGlobal(:)
    
    call rpn_comm_gather(lst_bhi%nla, 1, "mpi_integer",       &
                         allnlaLocal, 1, "mpi_integer", 0, "GRID", ier)
    call rpn_comm_gather(ilaGlobal   , nlaMax, "mpi_integer",       &
                         allilaGlobal, nlaMax, "mpi_integer",0 ,"GRID", ier)

    deallocate(ilaGlobal)

    ! assign part of mpiglobal vector from current mpi process
    if (mpi_myid == 0) then

       allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mpi_nprocs))

       !$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,k,ila,p,ila_mpiglobal,jdim_mpiglobal)
       do jproc = 0, (mpi_nprocs-1)
          cv_allmaxmpilocal(:,jproc+1) = 0.d0
          
          do k = 1, nksdim
             do ila = 1, allnlaLocal(jproc+1)
                do p = 1, lst_bhi%nphase

                   jdim_mpilocal = ( (k-1) * allnlaLocal(jproc+1) * lst_bhi%nphase ) + &
                                                        ( (ila-1) * lst_bhi%nphase ) + p

                   ila_mpiglobal = allilaGlobal(ila,jproc+1)
                   if ( ila_mpiglobal <= 0 ) then 
                      write(*,*) 'lbhi_reduceToMPILocal: invalid ila_mpiglobal index ', ila_mpiglobal
                      call utl_abort('lbhi_reduceToMPILocal')
                   end if

                   jdim_mpiglobal = ( (k-1) * lst_bhi%nlaGlobal * lst_bhi%nphase ) + &
                                            ( (ila_mpiglobal-1) * lst_bhi%nphase ) + p
  
                   if (jdim_mpilocal > cvDim_allMpiLocal(jproc+1)) then
                      write(*,*)
                      write(*,*) 'ERROR: jdim_mpilocal > cvDim_allMpiLocal(jproc+1)', jdim_mpilocal, cvDim_allMpiLocal(jproc+1) 
                      write(*,*) '       proc, k, ila, p = ',jproc,k,ila,p
                      call utl_abort('lbhi_reduceToMPILocal')
                   end if
                   if (jdim_mpiglobal > cvDim_mpiglobal) then
                      write(*,*)
                      write(*,*) 'ERROR: jdim_mpiglobal > cvDim_mpiglobal', jdim_mpiglobal, cvDim_mpiglobal
                      write(*,*) '       proc, k, ila, p = ',jproc,k,ila,p
                      call utl_abort('lbhi_reduceToMPILocal')
                   end if
                   
                   cv_allmaxmpilocal(jdim_mpilocal,jproc+1) = cv_mpiglobal(jdim_mpiglobal)
                   
                end do
             end do
          end do
       end do
       !$OMP END PARALLEL DO

    else
       allocate(cv_allmaxmpilocal(1,1))
    end if

    deallocate(allnlaLocal)
    deallocate(allilaGlobal)

    !- Distribute
    allocate(displs(mpi_nprocs))
    !$OMP PARALLEL DO PRIVATE(jproc)
    do jproc = 0, (mpi_nprocs-1)
       displs(jproc+1) = jproc*cvDim_maxMpiLocal ! displacement wrt cv_allMaxMpiLocal from which
                                                 ! to take the outgoing data to process jproc
    end do
    !$OMP END PARALLEL DO

    call rpn_comm_scatterv(cv_allMaxMpiLocal, cvDim_allMpiLocal, displs, "mpi_real4", &
                           cv_mpiLocal      , cvDim , "mpi_real4", &
                           0, "GRID", ier)

   deallocate(displs) 
   deallocate(cv_allMaxMpiLocal)
   deallocate(cvDim_allMpiLocal)


  END SUBROUTINE lbhi_reduceToMPILocal_r4

  !--------------------------------------------------------------------------
  ! lbhi_expandToMPIGlobal
  !--------------------------------------------------------------------------
  SUBROUTINE lbhi_expandToMPIGlobal(cv_mpilocal,cv_mpiglobal)
    implicit none
    real(8), intent(in)  :: cv_mpilocal(cvDim)
    real(8), intent(out) :: cv_mpiglobal(:)

    real(8), allocatable :: cv_maxmpilocal(:)
    real(8), pointer     :: cv_allmaxmpilocal(:,:) => null()

    integer, allocatable :: cvDim_allMpilocal(:)

    integer, allocatable :: ilaGlobal(:), allnlaLocal(:)
    integer, allocatable :: allilaGlobal(:,:)

    integer :: k, ila, p, ila_mpiglobal, jdim_mpilocal, jdim_mpiglobal
    integer :: ier, nlaMax, cvDim_maxmpilocal, jproc

    !
    !- 1.  Gather all local control vectors onto mpi task 0
    !
    allocate(cvDim_allMpiLocal(mpi_nprocs))
    call rpn_comm_allgather(cvDim            ,1,"mpi_integer",       &
                            cvDim_allMpiLocal,1,"mpi_integer","GRID",ier)

    call rpn_comm_allreduce(cvDim,cvDim_maxmpilocal,1,"mpi_integer","mpi_max","GRID",ier)

    allocate(cv_maxmpilocal(cvDim_maxmpilocal))

    cv_maxmpilocal(:) = 0.0d0
    cv_maxmpilocal(1:cvDim) = cv_mpilocal(1:cvDim)

    nullify(cv_allmaxmpilocal)
    if (mpi_myid == 0) then
       allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mpi_nprocs))
    else
       allocate(cv_allmaxmpilocal(1,1))
    end if

    call rpn_comm_gather(cv_maxmpilocal,    cvDim_maxmpilocal, "mpi_double_precision",  &
                         cv_allmaxmpilocal, cvDim_maxmpilocal, "mpi_double_precision", 0, "GRID", ier )

    deallocate(cv_maxmpilocal)

    !
    !- 2.  Reorganize gathered mpilocal control vectors into the mpiglobal control vector
    !

    call rpn_comm_allreduce(lst_bhi%nla,nlaMax,1,"mpi_integer","mpi_max","GRID",ier)

    if (mpi_myid == 0) then
       allocate(allnlaLocal(mpi_nprocs))
       allocate(allilaGlobal(nlaMax,mpi_nprocs))
    else
       allocate(allnlaLocal(1))
       allocate(allilaGlobal(1,1))
    end if
    
    allocate(ilaGlobal(nlaMax))
    ilaGlobal(:)             = -1
    ilaGlobal(1:lst_bhi%nla) = lst_bhi%ilaGlobal(:)

    call rpn_comm_gather(lst_bhi%nla, 1, "mpi_integer",       &
                         allnlaLocal, 1, "mpi_integer", 0, "GRID", ier)
    call rpn_comm_gather(ilaGlobal   , nlaMax, "mpi_integer",       &
                         allilaGlobal, nlaMax, "mpi_integer",0 ,"GRID", ier)

    deallocate(ilaGlobal)

    if (mpi_myid == 0) then
       cv_mpiglobal(:) = 0.0d0

       !$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,k,ila,p,ila_mpiglobal,jdim_mpiglobal)
       do jproc = 0, (mpi_nprocs-1)
          do k = 1, nksdim
             do ila = 1, allnlaLocal(jproc+1)
                do p = 1, lst_bhi%nphase

                   jdim_mpilocal = ( (k-1) * allnlaLocal(jproc+1) * lst_bhi%nphase ) + &
                                                        ( (ila-1) * lst_bhi%nphase ) + p

                   ila_mpiglobal = allilaGlobal(ila,jproc+1)
                   if ( ila_mpiglobal <= 0 ) then 
                      write(*,*) 'lbhi_expandToMPIGlobal: invalid ila_mpiglobal index ', ila_mpiglobal
                      call utl_abort('lbhi_expandToMPIGlobal')
                   end if

                   jdim_mpiglobal = ( (k-1) * lst_bhi%nlaGlobal * lst_bhi%nphase ) + &
                                            ( (ila_mpiglobal-1) * lst_bhi%nphase ) + p
  
                   if (jdim_mpilocal > cvDim_allMpiLocal(jproc+1)) then
                      write(*,*)
                      write(*,*) 'ERROR: jdim_mpilocal > cvDim_allMpiLocal(jproc+1)', jdim_mpilocal, cvDim_allMpiLocal(jproc+1) 
                      write(*,*) '       proc, k, ila, p = ',jproc,k,ila,p
                      call utl_abort('lbhi_expandToMPIGlobal')
                   end if
                   if (jdim_mpiglobal > cvDim_mpiglobal) then
                      write(*,*)
                      write(*,*) 'ERROR: jdim_mpiglobal > cvDim_mpiglobal', jdim_mpiglobal, cvDim_mpiglobal
                      write(*,*) '       proc, k, ila, p = ',jproc,k,ila,p
                      call utl_abort('lbhi_expandToMPIGlobal')
                   end if
                   
                   cv_mpiglobal(jdim_mpiglobal) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)

                end do
             end do
          end do
       end do
       !$OMP END PARALLEL DO

    end if

    deallocate(allnlaLocal)
    deallocate(allilaGlobal)
    deallocate(cv_allmaxmpilocal)
    deallocate(cvDim_allMpiLocal)

  end SUBROUTINE lbhi_expandToMPIGlobal

  !--------------------------------------------------------------------------
  ! lbhi_expandToMPIGlobal_r4
  !--------------------------------------------------------------------------
  SUBROUTINE lbhi_expandToMPIGlobal_r4(cv_mpilocal,cv_mpiglobal)
    implicit none
    real(4), intent(in)  :: cv_mpilocal(cvDim)
    real(4), intent(out) :: cv_mpiglobal(:)

    real(4), allocatable :: cv_maxmpilocal(:)
    real(4), pointer     :: cv_allmaxmpilocal(:,:) => null()

    integer, allocatable :: cvDim_allMpilocal(:)

    integer, allocatable :: ilaGlobal(:), allnlaLocal(:)
    integer, allocatable :: allilaGlobal(:,:)

    integer :: k, ila, p, ila_mpiglobal, jdim_mpilocal, jdim_mpiglobal
    integer :: ier, nlaMax, cvDim_maxmpilocal, jproc

    !
    !- 1.  Gather all local control vectors onto mpi task 0
    !
    allocate(cvDim_allMpiLocal(mpi_nprocs))
    call rpn_comm_allgather(cvDim            ,1,"mpi_integer",       &
                            cvDim_allMpiLocal,1,"mpi_integer","GRID",ier)

    call rpn_comm_allreduce(cvDim,cvDim_maxmpilocal,1,"mpi_integer","mpi_max","GRID",ier)

    allocate(cv_maxmpilocal(cvDim_maxmpilocal))

    cv_maxmpilocal(:) = 0.0d0
    cv_maxmpilocal(1:cvDim) = cv_mpilocal(1:cvDim)

    nullify(cv_allmaxmpilocal)
    if (mpi_myid == 0) then
       allocate(cv_allmaxmpilocal(cvDim_maxmpilocal,mpi_nprocs))
    else
       allocate(cv_allmaxmpilocal(1,1))
    end if

    call rpn_comm_gather(cv_maxmpilocal,    cvDim_maxmpilocal, "mpi_real4",  &
                         cv_allmaxmpilocal, cvDim_maxmpilocal, "mpi_real4", 0, "GRID", ier )

    deallocate(cv_maxmpilocal)

    !
    !- 2.  Reorganize gathered mpilocal control vectors into the mpiglobal control vector
    !

    call rpn_comm_allreduce(lst_bhi%nla,nlaMax,1,"mpi_integer","mpi_max","GRID",ier)

    if (mpi_myid == 0) then
       allocate(allnlaLocal(mpi_nprocs))
       allocate(allilaGlobal(nlaMax,mpi_nprocs))
    else
       allocate(allnlaLocal(1))
       allocate(allilaGlobal(1,1))
    end if
    
    allocate(ilaGlobal(nlaMax))
    ilaGlobal(:)             = -1
    ilaGlobal(1:lst_bhi%nla) = lst_bhi%ilaGlobal(:)

    call rpn_comm_gather(lst_bhi%nla, 1, "mpi_integer",       &
                         allnlaLocal, 1, "mpi_integer", 0, "GRID", ier)
    call rpn_comm_gather(ilaGlobal   , nlaMax, "mpi_integer",       &
                         allilaGlobal, nlaMax, "mpi_integer",0 ,"GRID", ier)

    deallocate(ilaGlobal)

    if (mpi_myid == 0) then
       cv_mpiglobal(:) = 0.0d0

       !$OMP PARALLEL DO PRIVATE(jproc,jdim_mpilocal,k,ila,p,ila_mpiglobal,jdim_mpiglobal)
       do jproc = 0, (mpi_nprocs-1)
          do k = 1, nksdim
             do ila = 1, allnlaLocal(jproc+1)
                do p = 1, lst_bhi%nphase

                   jdim_mpilocal = ( (k-1) * allnlaLocal(jproc+1) * lst_bhi%nphase ) + &
                                                        ( (ila-1) * lst_bhi%nphase ) + p

                   ila_mpiglobal = allilaGlobal(ila,jproc+1)
                   if ( ila_mpiglobal <= 0 ) then 
                      write(*,*) 'lbhi_expandToMPIGlobal: invalid ila_mpiglobal index ', ila_mpiglobal
                      call utl_abort('lbhi_expandToMPIGlobal')
                   end if

                   jdim_mpiglobal = ( (k-1) * lst_bhi%nlaGlobal * lst_bhi%nphase ) + &
                                            ( (ila_mpiglobal-1) * lst_bhi%nphase ) + p
  
                   if (jdim_mpilocal > cvDim_allMpiLocal(jproc+1)) then
                      write(*,*)
                      write(*,*) 'ERROR: jdim_mpilocal > cvDim_allMpiLocal(jproc+1)', jdim_mpilocal, cvDim_allMpiLocal(jproc+1) 
                      write(*,*) '       proc, k, ila, p = ',jproc,k,ila,p
                      call utl_abort('lbhi_expandToMPIGlobal')
                   end if
                   if (jdim_mpiglobal > cvDim_mpiglobal) then
                      write(*,*)
                      write(*,*) 'ERROR: jdim_mpiglobal > cvDim_mpiglobal', jdim_mpiglobal, cvDim_mpiglobal
                      write(*,*) '       proc, k, ila, p = ',jproc,k,ila,p
                      call utl_abort('lbhi_expandToMPIGlobal')
                   end if
                   
                   cv_mpiglobal(jdim_mpiglobal) = cv_allmaxmpilocal(jdim_mpilocal,jproc+1)

                end do
             end do
          end do
       end do
       !$OMP END PARALLEL DO

    end if

    deallocate(allnlaLocal)
    deallocate(allilaGlobal)
    deallocate(cv_allmaxmpilocal)
    deallocate(cvDim_allMpiLocal)

  end SUBROUTINE lbhi_expandToMPIGlobal_r4

  !--------------------------------------------------------------------------
  ! lbhi_Finalize
  !--------------------------------------------------------------------------
  subroutine lbhi_Finalize
    implicit none

    integer :: var

    if (initialized) then
       deallocate(bsqrt)
       do var = 1, nControlVariable
          deallocate(ControlVariable(var)%GpStdDev)
          deallocate(ControlVariable(var)%ip1     )
       end do
    end if

  end subroutine lbhi_Finalize

end module lamBMatrixHI_mod
