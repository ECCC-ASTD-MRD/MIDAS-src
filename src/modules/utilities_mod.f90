
module utilities_mod
  ! MODULE utilities_mod (prefix='utl' category='8. Low-level utilities and constants')
  !
  !:Purpose: A place to collect numerous simple utility routines.
  !
  use mpi_f08 ! this is the Fortran 2008 MPI library module
  use omp_lib
  use netcdf
  use rmn_fst98
  use clibInterfaces_mod
  use randomNumber_mod
  use mathPhysConstants_mod
  use omp_lib

  implicit none
  save
  private

  ! Public procedures
  public :: utl_readNml, utl_flnml, utl_flnml_static
  public :: utl_fstlir,  utl_fstlir_r4, utl_fstecr
  public :: utl_unitConvMultFactor_r8, utl_unitConvMultFactor_r4
  public :: utl_writeStatus, utl_getfldprm, utl_abort
  public :: utl_printTime
  public :: utl_open_asciifile, utl_stnid_equal, utl_resize, utl_str
  public :: utl_get_stringId, utl_get_Id, utl_isNamelistPresent
  public :: utl_readFstField
  public :: utl_reAllocate
  public :: utl_heapsort2d
  public :: utl_heapsort1d
  public :: utl_isEqual
  public :: utl_combineString, utl_splitString, utl_removeEmptyStrings
  public :: utl_stringArrayToIntegerArray, utl_parseColumns
  public :: utl_copyFile, utl_findloc, utl_findlocs
  public :: utl_randomOrderInt, utl_cosDegrees
  public :: utl_tmg_start, utl_tmg_stop, utl_medianIndex
  public :: utl_fileType, utl_checkNetCDFstatus
  public :: utl_varPresentInNetcdfFile

  ! module interfaces
  ! -----------------

  ! interface for resizing arrays
  interface utl_resize
    module procedure utl_resize_1d_real
    module procedure utl_resize_1d_int
    module procedure utl_resize_1d_str
    module procedure utl_resize_2d_real
    module procedure utl_resize_3d_real
  end interface utl_resize

  ! interface for conversion to a left-justified string (useful for calls to utl_abort)
  interface utl_str
    module procedure utl_int2str
    module procedure utl_float2str
  end interface utl_str

  interface utl_reAllocate
    module procedure utl_reAllocate_char_1d
    module procedure utl_reAllocate_char_2d
    module procedure utl_reAllocate_char_3d
    module procedure utl_reAllocate_log_1d
    module procedure utl_reAllocate_log_2d
    module procedure utl_reAllocate_log_3d
    module procedure utl_reAllocate_int_1d
    module procedure utl_reAllocate_int_2d
    module procedure utl_reAllocate_int_3d
    module procedure utl_reAllocate_r4_1d
    module procedure utl_reAllocate_r8_1d
    module procedure utl_reAllocate_r4_2d
    module procedure utl_reAllocate_r8_2d
    module procedure utl_reAllocate_r4_3d
    module procedure utl_reAllocate_r8_3d
    module procedure utl_reAllocate_r4_4d
    module procedure utl_reAllocate_r8_4d
    module procedure utl_reAllocate_r4_5d
    module procedure utl_reAllocate_r8_5d
  end interface utl_reAllocate

  interface utl_findloc
    module procedure utl_findloc_logical
    module procedure utl_findloc_char
    module procedure utl_findloc_int
  end interface utl_findloc

  interface utl_findlocs
    module procedure utl_findlocs_char
    module procedure utl_findlocs_int
  end interface utl_findlocs

  interface utl_cosDegrees
    module procedure utl_cosDegrees_real4
    module procedure utl_cosDegrees_real8
  end interface utl_cosDegrees

  interface utl_isEqual
    module procedure utl_isEqual_real4
    module procedure utl_isEqual_real8
    module procedure utl_isEqual_real4Arrays
    module procedure utl_isEqual_real8Arrays
    module procedure utl_isEqual_real4ArraysScalar
    module procedure utl_isEqual_real8ArraysScalar
  end interface utl_isEqual

  ! For namelist reading
  character(len=:), target, allocatable :: utl_flnml, utl_flnml_static

contains

  subroutine utl_readNml()
    !
    !:Purpose: Read the namelist files into strings for quicker access.
    !
    !:Note:    It currently does not work correctly for "NAMBEN" which
    !          may have multiple instance of the namelist block within the
    !          same file.
    !
    implicit none

    ! Locals:
    integer :: myid, fileSize, ierr, nulnam, positionBeg, positionEnd
    logical :: fileExists

    write(*,*)
    write(*,*) 'utl_readNML: reading namelist files into strings for later use'

    ! We cannot use 'midasMPI_mod' modules variables 'mmpi_myid',
    ! 'mmpi_myidx' and 'mmpi_myidy' because depending on that module
    ! would introduce a circular dependency.
    call mpi_comm_rank(mpi_comm_world, myid, ierr)
    if ( ierr /= 0 ) then
      call utl_abort('MPI error raised in mpi_comm_rank called from utl_readNml')
    end if

    ! First read the file flnml which must exist
    inquire(file='./flnml', exist=fileExists, size=fileSize)
    if (fileExists) then
      allocate( character(len=filesize) :: utl_flnml )

      ! read flnml into a string
      nulnam = 0
      open(newunit=nulnam,file='./flnml',status='OLD',&
           form='UNFORMATTED',access='STREAM',iostat=ierr)
      read(nulnam, pos=1, iostat=ierr) utl_flnml
      close(nulnam, iostat=ierr)

      ! print to the listing
      if (myid == 0) then
        write(*,*)
        write(*,*) '============BEGIN CONTENTS OF FLNML================'
        write(*,'(A)') utl_flnml
        write(*,*) '============END   CONTENTS OF FLNML================'
      end if
    else
      call utl_abort('utl_readNml: The file "flnml" is not accessible')
    end if

    ! Check for and remove comments
    do
      positionBeg = index(utl_flnml, '!')
      if (positionBeg < 1) then
        ! No (more) comments found
        exit
      else
        ! Found a comment, replace it with blank space
        do positionEnd = positionBeg, positionBeg+1000
          if (utl_flnml(positionEnd:positionEnd) == new_line('A')) exit
        end do
        positionEnd = positionEnd - 1
        utl_flnml(positionBeg:PositionEnd) = ' '
      end if
    end do

    ! Second read the file flnml_static which may exist
    inquire(file='./flnml_static', exist=fileExists, size=fileSize)
    if (fileExists) then

      allocate( character(len=filesize) :: utl_flnml_static )

      ! read flnml_static into a string
      nulnam = 0
      open(newunit=nulnam,file='./flnml_static',status='OLD',&
           form='UNFORMATTED',access='STREAM',iostat=ierr)
      read(nulnam, pos=1, iostat=ierr) utl_flnml_static
      close(nulnam, iostat=ierr)

      ! print to the listing
      if (myid == 0) then
        write(*,*)
        write(*,*) '============BEGIN CONTENTS OF FLNML_STATIC================'
        write(*,'(A)') utl_flnml_static
        write(*,*) '============END   CONTENTS OF FLNML_STATIC================'
      end if

      ! Check for and remove comments
      do
        positionBeg = index(utl_flnml_static, '!')
        if (positionBeg < 1) then
          ! No (more) comments found
          exit
        else
          ! Found a comment, replace it with blank space
          do positionEnd = positionBeg, positionBeg+1000
            if (utl_flnml_static(positionEnd:positionEnd) == new_line('A')) exit
          end do
          positionEnd = positionEnd - 1
          utl_flnml_static(positionBeg:PositionEnd) = ' '
        end if
      end do

    else
      allocate( character(len=1) :: utl_flnml_static )
      utl_flnml_static = ' '
      write(*,*) 'utl_readNml: The file "flnml_static" is not accessible, will use default values'
    end if

    write(*,*)
    write(*,*) 'utl_readNML: finished'
    write(*,*)

  end subroutine utl_readNml

  function utl_fstlir(fld8, iun, ni, nj, nk, datev, etiket, &
                      ip1, ip2, ip3, typvar, nomvar) result(vfstlir)
    implicit none

    ! Arguments:
    real(8),          intent(inout) :: fld8(*)
    integer,          intent(in)    :: iun
    integer,          intent(inout) :: ni
    integer,          intent(inout) :: nj
    integer,          intent(inout) :: nk
    integer,          intent(in)    :: datev
    integer,          intent(in)    :: ip1
    integer,          intent(in)    :: ip2
    integer,          intent(in)    :: ip3
    character(len=*), intent(in)    :: etiket
    character(len=*), intent(in)    :: nomvar
    character(len=*), intent(in)    :: typvar
    ! Result:
    integer :: vfstlir

    ! Locals:
    integer :: key1,key2, ilen, jk1, jk2, jk3, la
    real(4), allocatable :: buffer4(:)

    !     Get field dimensions and allow memory for REAL copy of fld8.
    key1 = fstinf(iun, ni, nj, nk, datev, etiket, &
         ip1, ip2, ip3, typvar, nomvar)

    if(key1 >= 0) then
       ilen = ni*nj*nk
       allocate(buffer4(ilen))
       !     Read field
       key2 = fstluk(buffer4, key1, ni, nj, nk)
       if(key2 >= 0) then
          do jk3 = 1,nk
             do jk2 = 1,nj
                do jk1 = 1,ni
                   la=jk1+(jk2-1)*ni+(jk3-1)*ni*nj
                   fld8(la) = buffer4(la)
                end do
             end do
          end do
       end if

       deallocate(buffer4)
    end if

    vfstlir=key1

  end function utl_fstlir

  function utl_fstlir_r4(fld_r4, iun, ni, nj, nk, datev, etiket, &
                         ip1, ip2, ip3, typvar, nomvar) result(vfstlir)
    implicit none

    ! Arguments:
    real(4),          intent(inout) :: fld_r4(*)
    integer,          intent(in)    :: iun
    integer,          intent(inout) :: ni
    integer,          intent(inout) :: nj
    integer,          intent(inout) :: nk
    integer,          intent(in)    :: datev
    integer,          intent(in)    :: ip1
    integer,          intent(in)    :: ip2
    integer,          intent(in)    :: ip3
    character(len=*), intent(in)    :: etiket
    character(len=*), intent(in)    :: nomvar
    character(len=*), intent(in)    :: typvar
    ! Result:
    integer :: vfstlir

    ! Locals:
    integer :: key1,key2, ilen, jk1, jk2, jk3, la
    real(4), allocatable :: buffer_r4(:)

    !     Get field dimensions.
    key1 = fstinf(iun, ni, nj, nk, datev, etiket, &
         ip1, ip2, ip3, typvar, nomvar)

    if(key1 >= 0) then
       ilen = ni*nj*nk
       allocate(buffer_r4(ilen))
       !     Read field
       key2 = fstluk(buffer_r4, key1, ni, nj, nk)
       if(key2 >= 0) then
          do jk3 = 1,nk
             do jk2 = 1,nj
                do jk1 = 1,ni
                   la=jk1+(jk2-1)*ni+(jk3-1)*ni*nj
                   fld_r4(la) = buffer_r4(la)
                end do
             end do
          end do
       end if

       deallocate(buffer_r4)
    end if

    vfstlir=key1

  end function utl_fstlir_r4

  function utl_fstecr(fld8, npak, iun, dateo, deet, &
                      npas, ni, nj, nk, ip1, ip2, ip3, typvar, &
                      nomvar, etiket, grtyp, ig1, ig2, ig3, ig4, &
                      datyp, rewrit) result(vfstecr)
    implicit none

    ! Arguments:
    integer,          intent(in) :: ni
    integer,          intent(in) :: nj
    integer,          intent(in) :: nk
    real(8),          intent(in) :: fld8(ni,nj,nk)
    integer,          intent(in) :: iun
    integer,          intent(in) :: ip1
    integer,          intent(in) :: ip2
    integer,          intent(in) :: ip3
    integer,          intent(in) :: ig1
    integer,          intent(in) :: ig2
    integer,          intent(in) :: ig3
    integer,          intent(in) :: ig4
    integer,          intent(in) :: npak
    integer,          intent(in) :: dateo
    integer,          intent(in) :: deet
    integer,          intent(in) :: npas
    integer,          intent(in) :: datyp
    logical,          intent(in) :: rewrit
    character(len=*), intent(in) :: etiket
    character(len=*), intent(in) :: typvar
    character(len=*), intent(in) :: grtyp
    character(len=*), intent(in) :: nomvar
    ! Result:
    integer :: vfstecr

    ! Locals:
    real(4) :: work
    integer :: ikey, jk1, jk2, jk3
    real(4), allocatable :: buffer4(:,:,:)

    allocate(buffer4(ni,nj,nk))

    do jk3 = 1,nk
      do jk2 = 1,nj
        do jk1 = 1,ni
          buffer4(jk1,jk2,jk3) = real(fld8(jk1,jk2,jk3),4)
        end do
      end do
    end do

    ikey = fstecr(buffer4, work, npak, iun, dateo, deet, &
         npas, ni, nj, nk, ip1, ip2, ip3, typvar, nomvar, &
         etiket, grtyp, ig1, ig2, ig3, ig4, datyp, rewrit)

    deallocate(buffer4)

    vfstecr=ikey

  end function utl_fstecr

  !--------------------------------------------------------------------------
  ! utl_unitConvMultFactor_r4
  !--------------------------------------------------------------------------
  function utl_unitConvMultFactor_r4(varName, direction) result(multFactor)
    !
    !:Purpose: Return the multiplicative factor for unit conversion between
    !          FST files and MIDAS (real 4 version)
    !
    implicit none

    ! arguments
    character(len=*), intent(in) :: varName
    character(len=*), intent(in) :: direction ! toFSTfile (writing) or fromFSTfile (reading)
    ! output
    real(4) :: multFactor

    if (trim(direction) /= 'toFSTfile' .and. trim(direction) /= 'fromFSTfile' ) then
      call utl_abort('utl_unitConvMultFactor: invalid direction ' // direction )
    end if

    ! Multiplicative factor for data conversion
    select case(trim(varName))
    case ('P0','PB','UP')
      if (trim(direction) == 'fromFSTfile' ) then
        multFactor = mpc_pa_per_mbar_r4 ! hPa -> Pa
      else
        multFactor = 0.01 ! Pa -> hPa
      end if
    case ('UU','VV','UV')
      if (trim(direction) == 'fromFSTfile' ) then
        multFactor = mpc_m_per_s_per_knot_r4 ! knots -> m/s
      else
        multFactor = mpc_knots_per_m_per_s_r4 ! m/s -> knots
      end if
    case ('GZ')
      if (trim(direction) == 'fromFSTfile' ) then
        multFactor = 10.0 ! dam -> m
      else
        multFactor = real(1.0d0 / 10.0d0,4) ! m -> dam
      end if
    case ('TO3','O3L')
      if (trim(direction) == 'fromFSTfile' ) then
        multFactor = 1.0E9 * mpc_molar_mass_O3_r4 / &
                     mpc_molar_mass_dry_air_r4 ! vmr -> micrograms/kg
      else
        multFactor = 1.0E-9 * mpc_molar_mass_dry_air_r4 / &
                     mpc_molar_mass_O3_r4 ! micrograms/kg -> vmr
      end if
    case default
      multFactor = 1.0
    end select

  end function utl_unitConvMultFactor_r4

  !--------------------------------------------------------------------------
  ! utl_unitConvMultFactor_r8
  !--------------------------------------------------------------------------
  function utl_unitConvMultFactor_r8(varName, direction) result(multFactor)
    !
    !:Purpose: Return the multiplicative factor for unit conversion between
    !          FST files and MIDAS (real 8 version)
    !
    implicit none

    ! arguments
    character(len=*), intent(in) :: varName
    character(len=*), intent(in) :: direction ! toFSTfile (writing) or fromFSTfile (reading)
    ! output
    real(8) :: multFactor

    if (trim(direction) /= 'toFSTfile' .and. trim(direction) /= 'fromFSTfile' ) then
      call utl_abort('utl_unitConvMultFactor: invalid direction ' // direction )
    end if

    ! Multiplicative factor for data conversion
    select case(trim(varName))
    case ('P0','UP','PB')
      if (trim(direction) == 'fromFSTfile' ) then
        multFactor = mpc_pa_per_mbar_r8 ! hPa -> Pa
      else
        multFactor = 0.01d0 ! Pa -> hPa
      end if
    case ('UU','VV','UV')
      if (trim(direction) == 'fromFSTfile' ) then
        multFactor = mpc_m_per_s_per_knot_r8 ! knots -> m/s
      else
        multFactor = mpc_knots_per_m_per_s_r8 ! m/s -> knots
      end if
    case ('GZ')
      if (trim(direction) == 'fromFSTfile' ) then
        multFactor = 10.0d0 ! dam -> m
      else
        multFactor = 1.0d0 / 10.0d0 ! m -> dam
      end if
    case ('TO3','O3L')
      if (trim(direction) == 'fromFSTfile' ) then
        multFactor = 1.0d9 * mpc_molar_mass_O3_r8 / &
                     mpc_molar_mass_dry_air_r8 ! vmr -> micrograms/kg
      else
        multFactor = 1.0d-9 * mpc_molar_mass_dry_air_r8 / &
                     mpc_molar_mass_O3_r8 ! micrograms/kg -> vmr
      end if
    case default
      multFactor = 1.d0
    end select

  end function utl_unitConvMultFactor_r8

  !--------------------------------------------------------------------------
  ! utl_writeStatus
  !--------------------------------------------------------------------------
  subroutine utl_writeStatus(cmsg)
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: cmsg

    ! Locals:
    INTEGER :: iulstatus,fnom,fclos, ierr
    character(len=22):: clmsg

    clmsg='VAR3D_STATUS='//cmsg
    iulstatus = 0
    IERR =  FNOM(iulstatus,'VAR3D_STATUS.dot','SEQ+FMT',0)
    rewind (iulstatus)
    WRITE(iulstatus,'(a22)') clmsg
    ierr = fclos(iulstatus)

  end subroutine utl_writeStatus

  !--------------------------------------------------------------------------
  ! utl_printTime
  !--------------------------------------------------------------------------
  subroutine utl_printTime(reset_opt)
    !
    !:Purpose: Print the elapsed time in the listing. Use of the optional
    !          argument `reset_opt=.true.` resets the accumulator to zero.
    !
    implicit none

    ! Arguments:
    logical, optional, intent(in) :: reset_opt ! Allow user to reset the accumulator

    ! Locals:
    real(8), save :: startTime = -1.0d0
    real(8), save :: accumulatedStart = -1.0d0
    real(8), save :: previousTime = -1.0d0
    real(8)       :: currentTime
    logical, save :: firstCall = .true.
    logical       :: reset
    character(len=8)  :: dateString
    character(len=10) :: timeString

    if (present(reset_opt)) then
      reset = reset_opt
    else
      reset = .false.
    end if

    currentTime = omp_get_wtime()

    if (startTime < 0.0d0) then
      startTime = currentTime
    end if

    if (previousTime < 0.0d0) then
      previousTime = currentTime
    end if

    if (accumulatedStart < 0.0d0 .or. reset) then
      accumulatedStart = currentTime
    end if

    ! Also get the actual date and time
    call date_and_time(dateString, timeString)

    if (firstCall) then
      write(*,'(A,A)') &
           ' utl_printTime: '//dateString//' '//timeString(1:2)//'h '//timeString(3:4)//'m '//timeString(5:10)//'s, ', &
           'First call, counters initialized'
    end if

    if (reset .and. .not.firstCall) then
      write(*,'(A,A)') &
           ' utl_printTime: '//dateString//' '//timeString(1:2)//'h '//timeString(3:4)//'m '//timeString(5:10)//'s, ', &
           'Accumulator reset'
    end if

    if (.not. firstCall) then
      write(*,'(A,A,f10.4,A,A,f10.4,A,A,f10.4,A)') &
           ' utl_printTime: '//dateString//' '//timeString(1:2)//'h '//timeString(3:4)//'m '//timeString(5:10)//'s, ', &
                          'deltaT = ', (currentTime - previousTime), ' s, ', &
                          'accumT = ', (currentTime - accumulatedStart), ' s, ', &
                          'totalT = ', (currentTime - startTime), ' s'
    end if

    previousTime = currentTime
    firstCall = .false.

  end subroutine utl_printTime


  subroutine utl_getfldprm(kip1s,kip2,kip3,knlev,cdetiket,cdtypvar,kgid, &
                           cdvar,kstampv,knmaxlev,kinmpg,kip1style,kip1kind, &
                           ktrials,koutmpg)
    !
    !:Purpose:  Get 3D grid parameters for a specific trial field
    !           and check for consitancies between grid parameters
    !           of the levels.
    !
    implicit none

    ! Arguments:
    integer,          intent(out) :: kstampv
    integer,          intent(in)  :: knmaxlev
    integer,          intent(out) :: knlev
    integer,          intent(out) :: kgid
    integer,          intent(out) :: kip1s(knmaxlev)
    integer,          intent(out) :: kip1style
    integer,          intent(out) :: kip1kind
    integer,          intent(out) :: kip2
    integer,          intent(out) :: kip3
    integer,          intent(in)  :: ktrials
    integer,          intent(out) :: koutmpg
    integer,          intent(in)  :: kinmpg(ktrials)
    character(len=*), intent(out) :: cdtypvar
    character(len=*), intent(in)  :: cdvar
    character(len=*), intent(out) :: cdetiket

    ! Locals:
    integer :: ini,inj,ink,jlev,ier
    integer :: idateo, idateo2, idatyp, idatyp2, ideet, ideet2, idltf
    integer :: iextra1, iextra2, iextra3, iig12, iig22
    integer :: iig32, iig42, ilng, inbits,iig1,iig2,iig3,iig4
    integer :: inpas,inpas2, iswa, iubc, iip2, iip3
    integer :: ipmode,idate2,idate3,idatefull
    integer :: k,ier1
    real(4) :: zlev_r4
    character(len=12) :: cletiket
    character(len=4) :: clnomvar
    character(len=3) :: clnomvar_3
    character(len=2) :: cltypvar
    character(len=1) :: clgrtyp2,clgrtyp,clstring
    logical :: llflag
    integer :: ikeys(knmaxlev)
    ! external definitions
    integer, external :: newdate, ezqkdef

    knlev = 0

    do k=1,ktrials
       if(cdvar.eq.'U1') then
          clnomvar_3='UT1'
          ier = fstinl(kinmpg(k),INI,INJ, INK, -1, ' ', -1, -1, -1, &
               ' ',clnomvar_3,ikeys, knlev, knmaxlev)
       else if(cdvar.eq.'V1') then
          clnomvar_3='VT1'
          ier = fstinl(kinmpg(k),ini,inj, ink, -1, ' ', -1, -1, -1, &
               ' ',clnomvar_3,ikeys, knlev, knmaxlev)
       else
          ier = fstinl(kinmpg(k),INI,INJ, INK, -1, ' ', -1, -1, -1, &
               ' ',cdvar,IKEYS, KNLEV, knmaxlev)
       end if
       !
       if(knlev > 0 ) then
          ier1   = newdate(kstampv,idate2,idate3,-3)

          idatefull = idate2*100 + idate3/1000000
          idateo = MPC_missingValue_INT
          ideet = MPC_missingValue_INT
          inpas = MPC_missingValue_INT
          cdetiket = '-9999999'
          clgrtyp = '-'
          kip2 = MPC_missingValue_INT
          kip3 = MPC_missingValue_INT
          cdtypvar = '-'
          idatyp = MPC_missingValue_INT
          iig1 = MPC_missingValue_INT
          iig2 = MPC_missingValue_INT
          iig3 = MPC_missingValue_INT
          iig4 = MPC_missingValue_INT
          llflag = .true.
          koutmpg = kinmpg(k)
          exit
       end if
    end do ! End of loop k
    !
    if (knlev.gt.0) then
       do jlev = 1, knlev
          ier = fstprm(ikeys(jlev), idateo2, ideet2, inpas2, ini, inj, &
               ink,inbits,idatyp2, kip1s(jlev),iip2, iip3, &
               cltypvar,clnomvar,cletiket,clgrtyp2, iig12, iig22,iig32 &
               ,iig42,iswa,ilng,idltf,iubc,iextra1, iextra2, iextra3)
          llflag = (llflag.and.(idateo.eq.idateo2.or.idateo.eq.MPC_missingValue_INT))
          llflag = (llflag.and.(ideet.eq.ideet2.or.ideet.eq.MPC_missingValue_INT))
          llflag = (llflag.and.(inpas.eq.inpas2.or.inpas.eq.MPC_missingValue_INT))
          !          llflag = (llflag.and.(cdetiket.eq.cletiket.or.cdetiket.eq.
          !     &         '-9999999'))
          llflag = (llflag.and.(clgrtyp.eq.clgrtyp2.or.clgrtyp.eq.'-'))
          llflag = (llflag.and.(kip2.eq.iip2.or.kip2.eq.MPC_missingValue_INT))
          llflag = (llflag.and.(kip3.eq.iip3.or.kip3.eq.MPC_missingValue_INT))
          llflag = (llflag.and.(cdtypvar.eq.cltypvar.or.cdtypvar.eq.'-'))
          llflag = (llflag.and.(idatyp.eq.idatyp2.or.idatyp.eq.MPC_missingValue_INT))
          llflag = (llflag.and.(iig1.eq.iig12.or.iig1.eq.MPC_missingValue_INT))
          llflag = (llflag.and.(iig2.eq.iig22.or.iig2.eq.MPC_missingValue_INT))
          llflag = (llflag.and.(iig3.eq.iig32.or.iig3.eq.MPC_missingValue_INT))
          llflag = (llflag.and.(iig4.eq.iig42.or.iig4.eq.MPC_missingValue_INT))
          if (llflag) then
             idateo = idateo2
             ideet = ideet2
             inpas = inpas2
             cdetiket = cletiket
             clgrtyp = clgrtyp2
             kip2 = iip2
             kip3 = iip3
             cdtypvar = cltypvar
             idatyp = idatyp2
             iig1 = iig12
             iig2 = iig22
             iig3 = iig32
             iig4 = iig42
          else
             write(*,*) &
                  '****** Unit ', kinmpg &
                  ,' contains mixed dateo,deet,npas,etiket,grtyp,ip2,ip3' &
                  ,',typvar,datyp,ig1,ig2,ig3 and/or ig4 ' &
                  ,'for variable ',cdvar,' and datev, ',kstampv
             call utl_abort('GETFLDPRM2')
          end if
       end do
       !
       kgid = ezqkdef(ini,inj,clgrtyp,iig1,iig2,iig3,iig4,koutmpg)
       !
       !-------Determine the style in which ip1 is encoded (15bits or 31 bits)
       !       A value <= 32767 (2**16 -1)  means that ip1 is compacted in 15 bits
       !       Determine the type of P which was encoded in IP1
       !
       if(kip1s(1) .le. 32767) then
          kip1style = 3
       else
          kip1style = 2
       end if
       !
       !-------Determine the type of P  (see doc. of convip)
       !
       ipmode = -1
       call CONVIP(kip1s(1),zlev_r4,KIP1KIND, &
            ipmode,clstring, .false. )
    else
       do k=1,ktrials
          ier = fstinl(kinmpg(k),ini,inj, ink, -1, ' ', -1, -1, -1, &
               ' ',cdvar,ikeys, knlev, knmaxlev)
       end do
       write(*,*) 'Error - getfldprm2: no record found at time ' &
            ,idatefull,' for field ',cdvar,' but',knlev, &
            ' records found in unit ',kinmpg(k)
       call utl_abort('GETFLDPRM2')
    end if
    !
  end subroutine utl_getfldprm


  subroutine utl_abort(message)
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: message

    ! Locals:
    integer :: ierr

    write(6,9000) message
9000 format(//,4X,"!!!---ABORT---!!!",/,8X,"MIDAS stopped in ",A)
    flush(6)

    call mpi_abort(mpi_comm_world, 1, ierr)

  end subroutine utl_abort

  subroutine utl_open_asciifile(filename,unit)
    !
    !:Purpose: Opens an ascii file for output
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in)  :: filename
    integer,          intent(out) :: unit

    ! Locals:
    logical :: file_exists
    integer :: ier
    character(len=20) :: mode

    inquire(file=trim(filename), exist=file_exists)

    if (file_exists) then
       mode = 'FTN+APPEND+R/W'
    else
       mode = 'FTN+R/W'
    end if

    unit=0

    ier = utl_open_file(unit,trim(filename),trim(mode))

    if (ier.ne.0) call utl_abort('utl_open_messagefile: Error associating unit number')

  end subroutine utl_open_asciifile


  function utl_open_file(unit,filename,mode) result(ier)
    !
    !:Purpose: This is a temporary subroutine to open a file with fnom that is needed due to
    !          a bug in fnom that does not allow an ascii file to be opened in 'APPEND' mode.
    !
    implicit none

    ! Arguments:
    integer,          intent(inout) :: unit
    character(len=*), intent(in)    :: filename
    character(len=*), intent(in)    :: mode
    ! Result:
    integer :: ier

    ! Locals:
    character(len=10) :: position,action
    integer :: fnom

    if (index(mode,'APPEND').gt.0) then
       position = 'APPEND'
    else
       position = 'ASIS'
    end if

    if (index(mode,'R/W').gt.0) then
       action = 'READWRITE'
    else
       action = 'READ'
    end if

    ier = fnom(unit,filename,mode,0)

    close(unit=unit)
    open(unit=unit, file=filename, position=position, action=action)

  end function utl_open_file


  function utl_stnid_equal(id1,id2) result(same)
    !
    !:Purpose: Compares STNID values allowing for * as wildcards and trailing blanks
    !
    !:Arguments:
    !           :id1: reference stnid
    !           :id2: stnid being verified
    !           :same: logical indicating if id1 and id2 match
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: id1
    character(len=*), intent(in) :: id2
    ! Result:
    logical :: same

    ! Locals:
    integer :: ilen1,ilen2,ji

    same=.true.
    ilen1=len_trim(id1)
    ilen2=len_trim(id2)

    do ji=1,min(ilen1,ilen2)
       if ( id1(ji:ji).ne.'*' .and. id2(ji:ji).ne.'*' .and. id2(ji:ji).ne.id1(ji:ji) ) then
          same = .false.
          exit
       end if
    end do

    if (same.and.ilen1.gt.ilen2) then
       do ji=ilen2+1,ilen1
          if (id1(ji:ji).ne.'*') then
              same=.false.
              exit
          end if
       end do
    else if (same.and.ilen2.gt.ilen1) then
       do ji=ilen1+1,ilen2
          if (id2(ji:ji).ne.'*') then
              same=.false.
              exit
          end if
       end do
    end if

  end function utl_stnid_equal


  character(len=20) function utl_int2str(i)
    !
    !:Purpose: Function for integer to string conversion. Helpful when calling subroutine utl_abort.
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: i

    write(utl_int2str,*) i
    utl_int2str = adjustl(utl_int2str)

  end function utl_int2str


  character(len=20) function utl_float2str(x)
    !
    !:Purpose: Function for integer to string conversion. Helpful when calling subroutine utl_abort.
    !
    implicit none

    ! Arguments:
    real(8), intent(in) :: x

    write(utl_float2str,*) x
    utl_float2str = adjustl(utl_float2str)

  end function utl_float2str


  subroutine utl_resize_1d_real(arr,dim1)
    !
    !:Purpose: Resize 1D array
    !
    implicit none

    ! Arguments:
    real(8), pointer, intent(inout) :: arr(:)
    integer,          intent(in)    :: dim1

    ! Locals:
    real(8), pointer :: tmp(:)
    integer :: dim1_in,d1

    dim1_in = size(arr)
    d1 = min(dim1_in, dim1)

    allocate(tmp(dim1))
    tmp(1:d1) = arr(1:d1)

    if (dim1.gt.dim1_in) tmp(d1+1:dim1) = 0.0D0

    deallocate(arr)

    arr => tmp

    nullify(tmp)

  end subroutine utl_resize_1d_real


  subroutine utl_resize_1d_int(arr,dim1)
    !
    !:Purpose: Resize 1D array
    !
    implicit none

    ! Arguments:
    integer, pointer, intent(inout) :: arr(:)
    integer,          intent(in)    :: dim1

    ! Locals:
    integer, pointer :: tmp(:)
    integer :: dim1_in,d1

    dim1_in = size(arr)
    d1 = min(dim1_in, dim1)

    allocate(tmp(dim1))
    tmp(1:d1) = arr(1:d1)

    if (dim1.gt.dim1_in) tmp(d1+1:dim1) = 0

    deallocate(arr)

    arr => tmp

    nullify(tmp)

  end subroutine utl_resize_1d_int


  subroutine utl_resize_1d_str(arr,dim1)
    !
    !:Purpose: Resize 1D array
    !
    implicit none

    ! Arguments:
    character(len=*), pointer, intent(inout) :: arr(:)
    integer,                   intent(in)    :: dim1

    ! Locals:
    character(len=len(arr(1))), pointer :: tmp(:)
    integer :: dim1_in,d1

    dim1_in = size(arr)
    d1 = min(dim1_in, dim1)

    allocate(tmp(dim1))
    tmp(1:d1) = arr(1:d1)

    if (dim1.gt.dim1_in) tmp(d1+1:dim1) = ""

    deallocate(arr)
    arr => tmp
    nullify(tmp)

  end subroutine utl_resize_1d_str


  subroutine utl_resize_2d_real(arr,dim1,dim2)
    !
    !:Purpose: Resize 2D array
    !
    implicit none

    ! Arguments:
    real(8), pointer, intent(inout) :: arr(:,:)
    integer,          intent(in)    :: dim1
    integer,          intent(in)    :: dim2

    ! Locals:
    real(8), pointer :: tmp(:,:)
    integer :: dim1_in,dim2_in,d1,d2

    dim1_in = size(arr,dim=1)
    dim2_in = size(arr,dim=2)
    d1 = min(dim1_in, dim1)
    d2 = min(dim2_in, dim2)

    allocate(tmp(dim1,dim2))
    tmp(1:d1,1:d2) = arr(1:d1,1:d2)

    if (dim1.gt.dim1_in) tmp(d1+1:dim1,:) = 0.0D0
    if (dim2.gt.dim2_in) tmp(:,d2+1:dim2) = 0.0D0

    deallocate(arr)

    arr => tmp

    nullify(tmp)

  end subroutine utl_resize_2d_real


  subroutine utl_resize_3d_real(arr,dim1,dim2,dim3)
    !
    !:Purpose: Resize 3D array
    !
    implicit none

    ! Arguments:
    real(8), pointer, intent(inout) :: arr(:,:,:)
    integer,          intent(in)    :: dim1
    integer,          intent(in)    :: dim2
    integer,          intent(in)    :: dim3

    ! Locals:
    real(8), pointer :: tmp(:,:,:)
    integer :: dim1_in,dim2_in,dim3_in,d1,d2,d3

    dim1_in = size(arr,dim=1)
    dim2_in = size(arr,dim=2)
    dim3_in = size(arr,dim=3)
    d1 = min(dim1_in, dim1)
    d2 = min(dim2_in, dim2)
    d3 = min(dim3_in, dim3)

    allocate(tmp(dim1,dim2,dim3))
    tmp(1:d1,1:d2,1:d3) = arr(1:d1,1:d2,1:d3)

    if (dim1.gt.dim1_in) tmp(d1+1:dim1,:,:) = 0.0D0
    if (dim2.gt.dim2_in) tmp(:,d2+1:dim2,:) = 0.0D0
    if (dim3.gt.dim3_in) tmp(:,:,d3+1:dim3) = 0.0D0

    deallocate(arr)

    arr => tmp

    nullify(tmp)

  end subroutine utl_resize_3d_real


  subroutine utl_get_stringId(cstringin,nobslev,CList,NListSize,Nmax,elemId)
    !
    !:Purpose: Get element ID from a list of accumulating character strings (e.g. stnids).
    !          Called by filt_topoChm in filterobs_mod.ftn90
    !
    implicit none

    ! Arguments:
    integer,           intent(in)    :: Nmax
    integer,           intent(in)    :: nobslev
    integer,           intent(inout) :: NListSize
    integer,           intent(out)   :: elemId
    character(len=*),  intent(in)    :: cstringin
    character(len=*),  intent(inout) :: CList(Nmax)

    ! Locals:
    integer :: i
    character(len=120) :: cstring

    elemId=0
    if (NListSize.gt.Nmax-1) then
       write(*,*) 'utl_get_stringId: NListSize > Nmax-1 (', NListSize, '>', Nmax-1, ')'
       call utl_abort('utl_get_stringId: Dimension error, NListSize > Nmax-1.')
    else if (NListSize.gt.0) then
       if (nobslev.eq.1) then
          cstring=trim(cstringin)//'U'
          do i=1,NListSize
             if (trim(cstring).eq.trim(CList(i))) then
                 elemId=i
                 exit
             end if
          end do
       else
          cstring=trim(cstringin)
          do i=1,NListSize
             if (trim(cstring).eq.trim(CList(i))) then
                 elemId=i
                 exit
             end if
          end do
       end if

       if (elemId.eq.0) then
          do i=1,NListSize
             if (utl_stnid_equal(trim(CList(i)),trim(cstring))) then
                elemId=i
                exit
             end if
          end do
       end if
    end if

    if (elemID.eq.0) then
        NListSize=NListSize+1
        elemId=NListSize
        if (nobslev.eq.1) then
           CList(NListSize)=trim(cstringin)//'U'
        else
           CList(NListSize)=trim(cstringin)
        end if
    end if

  end subroutine utl_get_stringId


  subroutine utl_get_Id(id,IdList,NListSize,Nmax,elemId)
    !
    !:Purpose: Get element ID from list of accumulating integer IDs.
    !
    implicit none

    ! Arguments:
    integer, intent(in)    :: Nmax
    integer, intent(in)    :: id
    integer, intent(inout) :: NListSize
    integer, intent(inout) :: IdList(Nmax)
    integer, intent(out)   :: elemId

    ! Locals:
    integer :: i

    elemId=0
    if (NListSize.gt.Nmax-1) then
       call utl_abort('utl_get_Id: Dimension error, NListSize > Nmax-1.')
    else if (NListSize.gt.0) then
       do i=1,NListSize
          if (id.eq.IdList(i)) then
              elemId=i
              exit
          end if
       end do
    end if

    if (elemID.eq.0) then
        NListSize=NListSize+1
        elemId=NListSize
        IdList(NListSize)=id
    end if


  end subroutine utl_get_Id


  subroutine utl_readFstField( fname, varName, iip1, iip2, iip3, etiketi, &
                               ni, nj, nkeys, array, xlat_opt, xlong_opt, lvls_opt, kind_opt )
    !
    !:Purpose:  Read specified field from standard RPN/fst file. Could be one
    !           to all levels depending on the input iip1,iip2,iip3 values.
    !
    !           Currently assumes lat/long (or Gaussian) type grids.
    !           See hco_SetupFromFile for example toward future generalizations.
    !           Generalization would require having xlat and xlong being 2D.
    !
    !:Arguments:
    !                 :fname: input filename
    !                 :varName:  search nomvar
    !                 :iip1: search ip1
    !                 :iip2: search ip2
    !                 :iip3: search ip3
    !                 :etiketi: search etiket
    !                 :ni: ni values
    !                 :nj: nj values
    !                 :nkeys: number of records satisfying search criteria
    !                 :array: data arrray
    !                 :xlat_opt: 1D latitude array (optional)
    !                 :xlong_opt: 1D longitude array (optional)
    !                 :lvls_opt: 1D vertical coordinate array (optional)
    !                 :kind_opt: vertical coordinate type according to convip (optional)
    !
    implicit none

    ! Arguments:
    integer,                        intent(in)  :: iip1
    integer,                        intent(in)  :: iip2
    integer,                        intent(in)  :: iip3
    character(len=*),               intent(in)  :: varName
    character(len=*),               intent(in)  :: fname
    character(len=*),               intent(in)  :: etiketi
    integer,                        intent(out) :: ni
    integer,                        intent(out) :: nj
    integer,                        intent(out) :: nkeys
    integer, optional,              intent(out) :: kind_opt
    real(8), allocatable,           intent(out) :: array(:,:,:)
    real(8), allocatable, optional, intent(out) :: lvls_opt(:)
    real(8), allocatable, optional, intent(out) :: xlat_opt(:)
    real(8), allocatable, optional, intent(out) :: xlong_opt(:)

    ! Locals:
    real(4) :: lvl_r4
    logical :: Exists
    character(len=1) :: string
    integer :: iun
    integer :: i,ier, kindi
    integer, parameter :: maxkeys=1000
    integer :: keys(maxkeys),ini,inj,nk
    integer :: dateo, deet, npas, nbits, datyp
    integer :: ip1, ip2, ip3, swa, lng, dltf, ubc
    integer :: extra1, extra2, extra3
    integer :: ig1, ig2, ig3, ig4
    character(1) clgrtyp
    character(2) cltypvar
    character(4) nomvar
    character(12) cletiket
    real(4), allocatable :: buffer(:,:)
    real(4), allocatable :: buffer3D(:,:,:)
    real :: xlat1_4, xlon1_4, xlat2_4, xlon2_4, dincr
    ! external definitions
    integer, external :: fnom, fclos

    ! Open file
    iun = 0
    inquire(file=trim(fname),exist=Exists)
    if(.not.Exists) then
      call utl_abort('utl_readFstField: Did not find file ' // trim(fname))
    else
      ier=fnom(iun,trim(fname),'RND+OLD+R/O',0)
      ier=fstouv(iun,'RND+OLD')
    end if

    ! Find reports in file for specified varName and iip*.
    ier = fstinl(iun,ni,nj,nk,-1,etiketi,iip1,iip2,iip3,'',varName,keys,nkeys,maxkeys)

    if(ier.lt.0.or.nkeys.eq.0) then
      call utl_abort('utl_readFstField: Search field missing ' // trim(varName) // &
           ' from file ' // trim(fname) // '. IPs and etiket: ' // &
           utl_str(iip1) // ', ' // utl_str(iip2) // ', ' // utl_str(iip3) // &
           ',  ' // trim(etiketi) // '.')
    else if (nk.gt.1) then
      if (nkeys > 1 .or. present(kind_opt) .or. present(lvls_opt) ) then
        call utl_abort('utl_readFstField: Unexpected size nk ' // trim(utl_str(nk)) // &
             ' for ' // trim(varName) // ' of file ' // trim(fname))
      end if
    end if

    if (present(xlat_opt).and.present(xlong_opt)) then

      !  Get lat and long if available.

      if (allocated(xlat_opt)) deallocate(xlat_opt,xlong_opt)
      allocate(xlat_opt(nj),xlong_opt(ni),buffer(ni*nj,1))
      xlat_opt(:)=MPC_missingValue_R8
      xlong_opt(:)=MPC_missingValue_R8

      ier = fstprm(keys(1),dateo, deet, npas, ni, nj, nk, nbits,    &
           datyp, ip1, ip2, ip3, cltypvar, nomvar, cletiket, &
           clgrtyp, ig1, ig2, ig3,                           &
           ig4, swa, lng, dltf, ubc, extra1, extra2, extra3)

      if (trim(clgrtyp) /= 'B') then
        if (ni.gt.1) then
          ier=fstlir(buffer,iun,ni,inj,nk,-1,'',ig1,ig2,ig3,'','>>')
          if (ier.ge.0) xlong_opt(:)=buffer(1:ni,1)
        end if
        if (nj.gt.1) then
          ier=fstlir(buffer,iun,ini,nj,nk,-1,'',ig1,ig2,ig3,'','^^')
          if (ier.ge.0) xlat_opt(:)=buffer(1:nj,1)
        end if
        deallocate(buffer)
      end if

      if ( trim(clgrtyp) == 'Z' ) then

        ! Check for rotated grid

        ier = fstprm(ier,                                         & ! IN
             dateo, deet, npas, ni, nj, nk, nbits,             & ! OUT
             datyp, ip1, ip2, ip3, cltypvar, nomvar, cletiket, & ! OUT
             clgrtyp, ig1, ig2, ig3,                           & ! OUT
             ig4, swa, lng, dltf, ubc, extra1, extra2, extra3 )  ! OUT

        call cigaxg (clgrtyp,                                     & ! IN
             xlat1_4, xlon1_4, xlat2_4, xlon2_4,          & ! OUT
             ig1, ig2, ig3, ig4 ) ! IN

        if ( .not. utl_isEqual(xlat1_4, xlat2_4) .or. .not. utl_isEqual(xlon1_4,xlon2_4) ) &
             call utl_abort('utl_readFstField: Cannot currently handle rotated grid')

      else if (trim(clgrtyp) == 'B') then

        ! Set B type lat-long grid

        dincr=360.0e0/(ni-1)
        do i=1,ni
          xlong_opt(i) = (i-1)*dincr
        end do

        if (ig1 == 0) then
          ! Global
          dincr=180.0e0/(nj-1)
          if (ig2 == 0) then
            do i=1,nj
              xlat_opt(i) = -90.0 + (i-1)*dincr
            end do
          else
            do i=1,nj
              xlat_opt(i) = 90.0 - (i-1)*dincr
            end do
          end if
        else if (ig1 == 1) then
          ! Northern hemispheric
          dincr=90.0e0/(nj-1)
          if (ig2 == 0) then
            do i=1,nj
              xlat_opt(i) = 0.0 + (i-1)*dincr
            end do
          else
            do i=1,nj
              xlat_opt(i) = 90.0 - (i-1)*dincr
            end do
          end if
        else
          ! Southern hemispheric
          dincr=90.0e0/(nj-1)
          if (ig2 == 0) then
            do i=1,nj
              xlat_opt(i) = -90.0 + (i-1)*dincr
            end do
          else
            do i=1,nj
              xlat_opt(i) = 0 - (i-1)*dincr
            end do
          end if
        end if

      else if (trim(clgrtyp) /= 'G') then

        call utl_abort('utl_readFstField: Cannot currently handle grid type ' // trim(clgrtyp) )

      end if

    end if

    ! Get vertical coordinate

    if (present(lvls_opt)) then
      if (allocated(lvls_opt)) deallocate(lvls_opt)
      allocate(lvls_opt(nkeys))

      do i=1,nkeys
        ier = fstprm(keys(i),dateo, deet, npas, ni, nj, nk, nbits,    &
             datyp, ip1, ip2, ip3, cltypvar, nomvar, cletiket, &
             clgrtyp, ig1, ig2, ig3,                           &
             ig4, swa, lng, dltf, ubc, extra1, extra2, extra3)
        call convip(ip1,lvl_r4,kindi,-1,string,.false.)
        lvls_opt(i)=lvl_r4
      end do
    end if

    if (present(kind_opt)) then
      if (present(lvls_opt)) then
        kind_opt=kindi
      else
        kind_opt=-1
      end if
    end if

    ! Get field

    if ( nk == 1 ) then
      if (allocated(buffer)) deallocate(buffer)
      allocate(array(ni,nj,nkeys),buffer(ni,nj))
      do i=1,nkeys
        ier=fstluk(buffer,keys(i),ni,nj,nk)
        array(:,:,i)=buffer(:,:)
      end do
      deallocate(buffer)
    else
      if (allocated(buffer3D)) deallocate(buffer3D)
      allocate(array(ni,nj,nk),buffer3D(ni,nj,nk))
      ier=fstluk(buffer,keys(1),ni,nj,nk)
      array(:,:,:)=buffer3D(:,:,:)
      deallocate(buffer3D)
      nkeys=nk
    end if

    ier=fstfrm(iun)
    ier=fclos(iun)

  end subroutine utl_readFstField

  function utl_fileType(fileName_opt) result(fileType)
    !
    !:Purpose: Return a string that identifies the file type, relying mostly
    !          on the `rmnlib` function `wkoffit`.
    !
    implicit none

    ! arguments
    character(len=*), optional, intent(in) :: fileName_opt
    ! output
    character(len=20) :: fileType

    ! locals
    integer :: wkoffit, typeCode

    ! Temporary fix for testing purposes
    if (.not. present(fileName_opt)) then
      fileType = 'FST'
      return
    end if

    typeCode = wkoffit(trim(fileName_opt))
    select case(typeCode)
    case (1,2,3,33,34,39)
      fileType = 'FST'
    case (6)
      fileType = 'BURP'
    case (35,38)
      fileType = 'NetCDF'
    case (41)
      fileType = 'sqliteOrObsdb'
    case default
      ! check if filename contain '.nc' in case it is a new netCDF version
      if (index(trim(fileName_opt),'.nc') /= 0) then
        write(*,*) 'utl_fileType: assume NetCDF based on file name extension'
        fileType = 'NetCDF'
      else
        write(*,*) 'utl_fileType: fileName     = ', trim(fileName_opt)
        write(*,*) 'utl_fileType: wkoffit code = ', typeCode
        call utl_abort('utl_fileType: unknown file type')
      end if
    end select

  end function utl_fileType


  subroutine utl_reAllocate_char_1d(array,dim1)
    implicit none

    ! Arguments:
    character(len=128), allocatable, intent(inout) :: array(:)
    integer,                         intent(in) :: dim1

    if( allocated(array) ) then
      if ( size(array) == dim1 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1))

  end subroutine utl_reAllocate_char_1d

  subroutine utl_reAllocate_char_2d(array,dim1,dim2)
    implicit none

    ! Arguments:
    character(len=128), allocatable, intent(inout) :: array(:,:)
    integer,                         intent(in)    :: dim1
    integer,                         intent(in)    :: dim2

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2))

  end subroutine utl_reAllocate_char_2d

  subroutine utl_reAllocate_char_3d(array,dim1,dim2,dim3)
    implicit none

    ! Arguments:
    character(len=128), allocatable, intent(inout) :: array(:,:,:)
    integer,                         intent(in)    :: dim1
    integer,                         intent(in)    :: dim2
    integer,                         intent(in)    :: dim3

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2*dim3 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2,dim3))

  end subroutine utl_reAllocate_char_3d

  subroutine utl_reAllocate_log_1d(array,dim1)
    implicit none

    ! Arguments:
    logical, allocatable, intent(inout) :: array(:)
    integer,              intent(in)    :: dim1

    if( allocated(array) ) then
      if ( size(array) == dim1 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1))
    array(:) = .true.

  end subroutine utl_reAllocate_log_1d

  subroutine utl_reAllocate_log_2d(array,dim1,dim2)
    implicit none

    ! Arguments:
    logical, allocatable, intent(inout) :: array(:,:)
    integer,              intent(in)    :: dim1
    integer,              intent(in)    :: dim2

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2))
    array(:,:) = .true.

  end subroutine utl_reAllocate_log_2d

  subroutine utl_reAllocate_log_3d(array,dim1,dim2,dim3)
    implicit none

    ! Arguments:
    logical, allocatable, intent(inout) :: array(:,:,:)
    integer,              intent(in)    :: dim1
    integer,              intent(in)    :: dim2
    integer,              intent(in)    :: dim3

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2*dim3 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2,dim3))
    array(:,:,:) = .true.

  end subroutine utl_reAllocate_log_3d

  subroutine utl_reAllocate_int_1d(array,dim1)
    implicit none

    ! Arguments:
    integer, allocatable, intent(inout) :: array(:)
    integer,              intent(in)    :: dim1

    if( allocated(array) ) then
      if ( size(array) == dim1 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1))
    array(:) = 0d0

  end subroutine utl_reAllocate_int_1d

  subroutine utl_reAllocate_int_2d(array,dim1,dim2)
    implicit none

    ! Arguments:
    integer, allocatable, intent(inout) :: array(:,:)
    integer,              intent(in)    :: dim1
    integer,              intent(in)    :: dim2

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2))
    array(:,:) = 0d0

  end subroutine utl_reAllocate_int_2d

  subroutine utl_reAllocate_int_3d(array,dim1,dim2,dim3)
    implicit none

    ! Arguments:
    integer, allocatable, intent(inout) :: array(:,:,:)
    integer,              intent(in)    :: dim1
    integer,              intent(in)    :: dim2
    integer,              intent(in)    :: dim3

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2*dim3 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2,dim3))
    array(:,:,:) = 0d0

  end subroutine utl_reAllocate_int_3d

  subroutine utl_reAllocate_r4_1d(array,dim1)
    implicit none

    ! Arguments:
    real(4), allocatable, intent(inout) :: array(:)
    integer,              intent(in)    :: dim1

    if( allocated(array) ) then
      if ( size(array) == dim1 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1))
    array(:) = 0.0d0

  end subroutine utl_reAllocate_r4_1d

  subroutine utl_reAllocate_r8_1d(array,dim1)
    implicit none

    ! Arguments:
    real(8), allocatable, intent(inout) :: array(:)
    integer,              intent(in)    :: dim1

    if( allocated(array) ) then
      if ( size(array) == dim1 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1))
    array(:) = 0.0d0

  end subroutine utl_reAllocate_r8_1d

  subroutine utl_reAllocate_r4_2d(array,dim1,dim2)
    implicit none

    ! Arguments:
    real(4), allocatable, intent(inout) :: array(:,:)
    integer,              intent(in)    :: dim1
    integer,              intent(in)    :: dim2

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2))
    array(:,:) = 0.0d0

  end subroutine utl_reAllocate_r4_2d

  subroutine utl_reAllocate_r8_2d(array,dim1,dim2)
    implicit none

    ! Arguments:
    real(8), allocatable, intent(inout) :: array(:,:)
    integer,              intent(in)    :: dim1
    integer,              intent(in)    :: dim2

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2))
    array(:,:) = 0.0d0

  end subroutine utl_reAllocate_r8_2d

  subroutine utl_reAllocate_r4_3d(array,dim1,dim2,dim3)
    implicit none

    ! Arguments:
    real(4), allocatable, intent(inout) :: array(:,:,:)
    integer,              intent(in)    :: dim1
    integer,              intent(in)    :: dim2
    integer,              intent(in)    :: dim3

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2*dim3 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2,dim3))
    array(:,:,:) = 0.0d0

  end subroutine utl_reAllocate_r4_3d


  subroutine utl_reAllocate_r8_3d(array,dim1,dim2,dim3)
    implicit none

    ! Arguments:
    real(8), allocatable, intent(inout) :: array(:,:,:)
    integer,              intent(in)    :: dim1
    integer,              intent(in)    :: dim2
    integer,              intent(in)    :: dim3

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2*dim3 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2,dim3))
    array(:,:,:) = 0.0d0

  end subroutine utl_reAllocate_r8_3d


  subroutine utl_reAllocate_r4_4d(array,dim1,dim2,dim3,dim4)
    implicit none

    ! Arguments:
    real(4), allocatable, intent(inout) :: array(:,:,:,:)
    integer,              intent(in)    :: dim1
    integer,              intent(in)    :: dim2
    integer,              intent(in)    :: dim3
    integer,              intent(in)    :: dim4

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2*dim3*dim4 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2,dim3,dim4))
    array(:,:,:,:) = 0.0d0

  end subroutine utl_reAllocate_r4_4d


  subroutine utl_reAllocate_r8_4d(array,dim1,dim2,dim3,dim4)
    implicit none

    ! Arguments:
    real(8), allocatable, intent(inout) :: array(:,:,:,:)
    integer,              intent(in)    :: dim1
    integer,              intent(in)    :: dim2
    integer,              intent(in)    :: dim3
    integer,              intent(in)    :: dim4

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2*dim3*dim4 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2,dim3,dim4))
    array(:,:,:,:) = 0.0d0

  end subroutine utl_reAllocate_r8_4d


  subroutine utl_reAllocate_r4_5d(array,dim1,dim2,dim3,dim4,dim5)
    implicit none

    ! Arguments:
    real(4), allocatable, intent(inout) :: array(:,:,:,:,:)
    integer,              intent(in)    :: dim1
    integer,              intent(in)    :: dim2
    integer,              intent(in)    :: dim3
    integer,              intent(in)    :: dim4
    integer,              intent(in)    :: dim5

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2*dim3*dim4*dim5 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2,dim3,dim4,dim5))
    array(:,:,:,:,:) = 0.0d0

  end subroutine utl_reAllocate_r4_5d


  subroutine utl_reAllocate_r8_5d(array,dim1,dim2,dim3,dim4,dim5)
    implicit none

    ! Arguments:
    real(8), allocatable, intent(inout) :: array(:,:,:,:,:)
    integer,              intent(in)    :: dim1
    integer,              intent(in)    :: dim2
    integer,              intent(in)    :: dim3
    integer,              intent(in)    :: dim4
    integer,              intent(in)    :: dim5

    if( allocated(array) ) then
      if ( size(array) == dim1*dim2*dim3*dim4*dim5 ) then
        return
      else
        deallocate(array)
      end if
    end if

    allocate(array(dim1,dim2,dim3,dim4,dim5))
    array(:,:,:,:,:) = 0.0d0

  end subroutine utl_reAllocate_r8_5d


  subroutine utl_heapsort2d(array)
    !
    !:Purpose: Sort a real 2D array in ascending order according
    !          to the first column
    !
    implicit none

    ! Arguments:
    real(4), intent(inout) :: array(:,:)

    ! Locals:
    real(4) :: values(2) ! temporary value
    integer :: i,j,nsize
    integer :: ileft,iright

    nsize  = size(array,1)
    ileft  = nsize/2+1
    iright = nsize

    if (nsize == 1) return

    do
      if(ileft > 1)then
        ileft = ileft-1
        values(:) = array(ileft,:)
      else
        values(:) = array(iright,:)
        array(iright,:) = array(1,:)
        iright = iright-1
        if (iright == 1) then
          array(1,:) = values(:)
          return
        end if
      end if
      i = ileft
      j = 2*ileft
      do while (j <= iright)
        if (j < iright) then
          if (array(j,1) < array(j+1,1)) j = j+1
        endif
        if (values(1) < array(j,1)) then
          array(i,:) = array(j,:)
          i = j
          j = j+j
        else
          j = iright+1
        end if
      end do
      array(i,:) = values(:)
    end do

  end subroutine utl_heapsort2d


  subroutine utl_heapsort1d(rvalues,indices)
    !
    !:Purpose: Sort a real 1D array in ascending order and give its indexes
    !
    implicit none

    ! Arguments:
    real(8), intent(inout) :: rvalues(:) ! 1D array of real values to be sorted
    integer, intent(inout) :: indices(:) ! indexes of the 1D array

    ! Locals:
    real(8) :: tmpval ! temporary value
    integer :: tmpind ! temporary value
    integer :: i,j,nsize
    integer :: ileft,iright

    if (size(rvalues) /= size(rvalues)) then
      call utl_abort('utl_heapsort1d: input arrays have different sizes.')
    endif
    nsize  = size(rvalues)
    ileft  = nsize/2+1
    iright = nsize

    if (nsize == 1) return

    do
      if(ileft > 1)then
        ileft=ileft-1
        tmpval = rvalues(ileft)
        tmpind = indices(ileft)
      else
        tmpval = rvalues(iright)
        tmpind = indices(iright)
        rvalues(iright) = rvalues(1)
        indices(iright) = indices(1)
        iright = iright-1
        if (iright == 1) then
          rvalues(1) = tmpval
          indices(1) = tmpind
          return
        end if
      end if
      i = ileft
      j = 2*ileft
      do while (j <= iright)
        if (j < iright) then
          if (rvalues(j) < rvalues(j+1)) j = j+1
        endif
        if (tmpval < rvalues(j)) then
          rvalues(i) = rvalues(j)
          indices(i) = indices(j)
          i = j
          j = j+j
        else
          j = iright+1
        end if
      end do
      rvalues(i) = tmpval
      indices(i) = tmpind
    end do

  end subroutine utl_heapsort1d


  subroutine utl_splitString(string, separator, stringArray)
    !
    !:Purpose: Divide a string into several parts using a specified separator
    !
    implicit none

    ! Arguments:
    character(len=*),              intent(in)    :: string         ! input string
    character(len=*),              intent(in)    :: separator      ! separator
    character(len=*), allocatable, intent(inout) :: stringArray(:) ! seperated strings

    ! Locals:
    integer :: stringArraySize, start, sepPos, substringIndex

    ! Calculate the number of substrings
    stringArraySize = count(transfer(string, 'a', len(string)) == separator) + 1
    allocate(stringArray(stringArraySize))

    start = 1
    substringIndex = 1
    do while (start <= len(string))
      sepPos = index(string(start:), separator)
      if (sepPos == 0) then
        stringArray(substringIndex) = string(start:)
        exit
      else
        stringArray(substringIndex) = string(start:start+sepPos-2)
        start = start + sepPos
        substringIndex = substringIndex + 1
      end if
    end do
  end subroutine utl_splitString

  subroutine utl_combineString(string,separator,stringArray)
    implicit none

    ! Arguments:
    character(len=*), intent(out) :: string
    character(len=*), intent(in)  :: separator
    character(len=*), intent(in)  :: stringArray(:)

    ! Locals:
    integer :: stringArraySize, stringIndex, stringCount, stringCountTotal

    stringArraySize = size(stringArray)

    stringCountTotal = 0
    do stringIndex = 1, size(stringArray)
      if (trim(stringArray(stringIndex)) == '') cycle
      stringCountTotal = stringCountTotal + 1
    end do

    string = ''
    stringCount = 0
    do stringIndex = 1, size(stringArray)
      if (trim(stringArray(stringIndex)) == '') cycle
      stringCount = stringCount + 1
      if (stringCount < stringCountTotal) then
        string = trim(string) // ' ' // trim(stringArray(stringIndex)) // trim(separator)
      else
        string = trim(string) // ' ' // trim(stringArray(stringIndex))
      end if
    end do

    write(*,*) 'utl_combineString: stringCountTotal = ', stringCountTotal
    write(*,*) 'utl_combineString: string     = ', trim(string)

  end subroutine utl_combineString


  subroutine utl_removeEmptyStrings(stringArray)
    implicit none

    ! Arguments:
    character(len=*), allocatable, intent(inout) :: stringArray(:)

    ! Locals:
    integer :: stringArraySize, stringIndex, stringCount, stringCountTotal
    character(len=len(stringArray(1))), allocatable :: newStringArray(:)

    stringArraySize = size(stringArray)

    stringCountTotal = 0
    do stringIndex = 1, size(stringArray)
      if (trim(stringArray(stringIndex)) == '') cycle
      stringCountTotal = stringCountTotal + 1
    end do

    allocate(newStringArray(stringCountTotal))

    stringCount = 0
    do stringIndex = 1, size(stringArray)
      if (trim(stringArray(stringIndex)) == '') cycle
      stringCount = stringCount + 1
      newStringArray(stringCount) = stringArray(stringIndex)
    end do

    deallocate(stringArray)
    allocate(stringArray(stringCountTotal))
    stringArray(:) = newStringArray(:)
    deallocate(newStringArray)

    write(*,*) 'utl_removeEmptyStrings: stringCountTotal = ', stringCountTotal
    write(*,*) 'utl_removeEmptyStrings: stringArray      = ', &
               (trim(stringArray(stringIndex))//' ', stringIndex = 1, size(stringArray))

  end subroutine utl_removeEmptyStrings


  subroutine utl_stringArrayToIntegerArray(stringArray,integerArray)
    implicit none

    ! Arguments:
    character(len=256),   intent(in)  :: stringArray(:)
    integer, allocatable, intent(out) :: integerArray(:)

    ! Locals:
    integer :: arraySize, arrayIndex

    arraySize = size(stringArray)

    allocate(integerArray(arraySize))

    do arrayIndex = 1, arraySize
      read(stringArray(arrayIndex),'(i5)')  integerArray(arrayIndex)
    end do

    write(*,*) 'utl_stringArrayToIntegerArray: integerArray = ', integerArray(:)

  end subroutine utl_stringArrayToIntegerArray

  !--------------------------------------------------------------------------
  ! utl_isNamelistPresent
  !--------------------------------------------------------------------------
  function utl_isNamelistPresent(namelistSectionName, namelistFileName, failMode_opt) result(found)
    !
    !:Purpose: To find if a namelist name tag is present in a namelist file
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: namelistSectionName
    character(len=*), intent(in) :: namelistFileName
    character(len=*), optional, intent(in) :: failMode_opt
    ! Result:
    logical :: found

    ! Locals:
    integer :: unit, fnom, fclos, ierr, positionBeg, positionEnd
    character(len=1000) :: line, failMode
    character(len=100)  :: word, namelistSectionNameUpper
    logical :: namelistExist
    character(len=:), pointer :: nameListStr_ptr

    ! Set action if namelist is missing
    if (present(failMode_opt)) then
      failMode = failMode_opt
    else
      failMode = 'ABORT'
    end if
    
    ! Check if namelistFileName is present
    inquire(file=namelistFileName,exist=namelistExist)
    if (.not. namelistExist) then
      if (trim(failMode) == 'ABORT') then
        call utl_abort('utl_isNamelistPresent: namelist file is missing : '// namelistFileName)
      else
        write(*,*)
        write(*,*) 'utl_isNamelistPresent: WARNING, namelist file is missing : ' // namelistFileName
        found = .false.
        return
      end if
    end if
    
    ! Ensure namelist section name is all in upper case
    namelistSectionNameUpper = namelistSectionName
    ierr = clib_toUpper(namelistSectionNameUpper)

    ! Check if namelist file is already read into a character string
    if (index(namelistFileName, 'flnml_static') > 0) then
      nameListStr_ptr => utl_flnml_static
    else if (index(namelistFileName, 'flnml') > 0) then
      nameListStr_ptr => utl_flnml
    else
      nullify(nameListStr_ptr)
    end if

    if (associated(nameListStr_ptr)) then

      ! Search for namelistSectionName in string
      found = .false.
      positionBeg = 1
      namelistLoop1 : do
        ! read next line from string
        positionEnd = positionBeg - 1 + index(nameListStr_ptr(positionBeg:),new_line('A'))
        read (nameListStr_ptr(positionBeg:(positionEnd-1)),"(a)",iostat=ierr) line ! read line into character variable
        positionBeg = positionEnd + 1
        if (positionBeg >= len(nameListStr_ptr)) exit namelistLoop1 ! reached end of the file
        if (trim(line) == '') cycle namelistLoop1                   ! skip empty lines

        read (line,*) word ! read first word of line
        ierr = clib_toUpper(word)
        if (trim(word) == '&'//trim(namelistSectionNameUpper)) then ! case insensitive
          ! found search string at beginning of line
          found = .true.
          exit
        end if
      end do namelistLoop1

    else

      write(*,*) 'utl_isNamelistPresent: Namelist file is not already in a string, checking the file instead'

      ! Open the namelist file
      unit=0
      ierr=fnom(unit,namelistFileName,'FTN+SEQ+R/O',0)

      ! Search for namelistSectionName in the file
      found = .false.
      namelistLoop2 : do
        read (unit,"(a)",iostat=ierr) line ! read line into character variable
        if (ierr /= 0) exit namelistLoop2
        if (trim(line) == "") cycle namelistLoop2 ! skip empty lines
        read (line,*) word          ! read first word of line
        ierr = clib_toUpper(word)
        if (trim(word) == '&'//trim(namelistSectionNameUpper)) then ! case insensitive
          ! found search string at beginning of line
          found = .true.
          exit namelistLoop2
        end if
      end do namelistLoop2

      ! Close the namelist file
      ierr=fclos(unit)

    end if

  end function utl_isNamelistPresent

  !-----------------------------------------------------------------
  ! utl_parseColumns
  !-----------------------------------------------------------------
  subroutine utl_parseColumns(line, numColumns, stringArray_opt)
    !
    !:Purpose: To return column values in array of strings and
    !          the number of space-delimited columns in a string
    !
    implicit none

    ! Arguments:
    character(len=*),           intent(in)  :: line
    integer,                    intent(out) :: numColumns
    character(len=*), optional, intent(out) :: stringArray_opt(:)

    ! Locals:
    integer :: linePosition, wordPosition, lineLength

    linePosition = 1
    lineLength = len_trim(line)
    numColumns = 0

    do while(linePosition <= lineLength)

      do while(line(linePosition:linePosition) == ' ')
        linePosition = linePosition + 1
        if (lineLength < linePosition) return
      end do

      numColumns = numColumns + 1
      wordPosition = 0
      if (present(stringArray_opt)) then
        stringArray_opt(numColumns) = ''
      end if

      do
        if (linePosition > lineLength) return
        if (line(linePosition:linePosition) == ' ') exit
        if (present(stringArray_opt)) then
          wordPosition = wordPosition + 1
          stringArray_opt(numColumns)(wordPosition:wordPosition) = line(linePosition:linePosition)
        end if
        linePosition = linePosition + 1
      end do

    end do

  end subroutine utl_parseColumns

  !--------------------------------------------------------------------------
  ! utl_copyFile
  !--------------------------------------------------------------------------
  function utl_copyFile(filein, fileout, concatenate_opt) result(status)
    !
    !:Purpose: Copy the specified file to the new location and/or name
    !          This function is very general, but was initially written to
    !          copy files from the disk to the ram disk
    !
    !
    implicit none

    ! Arguments:
    character(len=*),  intent(in) :: filein
    character(len=*),  intent(in) :: fileout
    logical, optional, intent(in) :: concatenate_opt ! default is .false.
    ! Result:
    integer :: status

    ! Locals:
    logical :: concatenate, exists
    integer :: ierr, unitin, unitout
    integer(8) :: numChar
    character :: bufferB
    integer, parameter :: bufferSizeKB = 1024
    character :: bufferKB(bufferSizeKB)
    integer, parameter :: bufferSizeMB = 1024*1024
    character :: bufferMB(bufferSizeMB)
    character(len=7) :: open_status ! Will contain 'NEW', 'OLD' or 'REPLACE'
    character(len=6) :: position    ! Will contain 'ASIS' or 'APPEND'

    write(*,*) 'utl_copyFile: copy from ', trim(filein), ' to ', trim(fileout)

    call utl_tmg_start(175,'low-level--utl_copyFile')

    unitin=10
    open(unit=unitin, file=trim(filein), status='OLD', form='UNFORMATTED', &
         action='READ', access='STREAM')

    concatenate = .false.
    if ( present(concatenate_opt) ) then
      concatenate = concatenate_opt
    end if

    unitout=11
    if ( concatenate ) then
      inquire(file=trim(fileout), exist=exists)
      if ( exists ) then
        open_status = 'OLD'
      else
        open_status = 'NEW'
      end if

      position = 'APPEND'
    else
      open_status = 'REPLACE'
      position = 'ASIS' ! this is the default
    end if

    open(unit=unitout, file=trim(fileout), status=open_status, form='UNFORMATTED', &
         action='WRITE', access='STREAM', position=position)

    numChar = 0
    do
      read(unitin,iostat=ierr) bufferMB
      if (ierr < 0) exit
      numChar = numChar + bufferSizeMB
      write(unitout) bufferMB
    end do

    do
      read(unitin,iostat=ierr,pos=numChar+1) bufferKB
      if (ierr < 0) exit
      numChar = numChar + bufferSizeKB
      write(unitout) bufferKB
    end do

    do
      read(unitin,iostat=ierr,pos=numChar+1) bufferB
      if (ierr < 0) exit
      numChar = numChar + 1
      write(unitout) bufferB
    end do

    write(*,*) 'utl_copyFile: copied ', numChar, ' bytes'

    close(unit=unitin)
    close(unit=unitout)

    if (numChar > 0) then
      status = 0
    else
      status = -1
      if (numChar == 0) then
        call utl_abort('utl_copyFile: ERROR, zero bytes copied')
      else
        ! Note: If 'numChar' becomes negative then it means it got bigger
        !       than the maximum integer the 'integer' type and so the
        !       variable 'numChar' wraps around and becomes negative.
        call utl_abort('utl_copyFile: ERROR, overflow detected since number of bytes copied is negative!')
      end if
    end if

    call utl_tmg_stop(175)

  end function utl_copyFile

  !--------------------------------------------------------------------------
  ! utl_findloc_logical
  !--------------------------------------------------------------------------
  function utl_findloc_logical(logicalArray, value_opt) result(location)
    !
    !:Purpose: A modified version of the fortran function `findloc`.
    !          If multiple matches are found in the array, a warning
    !          message is printed to the listing.
    !
    implicit none

    ! Arguments:
    logical,           intent(in) :: logicalArray(:)
    logical, optional, intent(in) :: value_opt
    ! Result:
    integer             :: location

    ! Locals:
    integer :: numFound, arrayIndex
    logical :: value

    if (present(value_opt)) then
      value = value_opt
    else
      value = .true.
    end if

    numFound = 0
    LOOP: do arrayIndex = 1, size(logicalArray)
      if (logicalArray(arrayIndex) .eqv. value) then
        numFound = numFound + 1
        ! return the first location found
        if (numFound == 1) location = arrayIndex
      end if
    end do LOOP

    ! give warning if more than 1 found
    if (numFound > 1) then
      write(*,*) 'utl_findloc_logical: found multiple locations of ', value
      write(*,*) 'utl_findloc_logical: number locations found =  ', numFound
    end if

    ! return zero if not found
    if (numFound == 0) then
      location = 0
    end if

  end function utl_findloc_logical

  !--------------------------------------------------------------------------
  ! utl_findloc_char
  !--------------------------------------------------------------------------
  function utl_findloc_char(charArray, value) result(location)
    !
    !:Purpose: A modified version of the fortran function `findloc`.
    !          If multiple matches are found in the array, a warning
    !          message is printed to the listing.
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: charArray(:)
    character(len=*), intent(in) :: value
    ! Result:
    integer                      :: location

    ! Locals:
    integer :: numFound, arrayIndex

    location = 0
    numFound = 0
    LOOP: do arrayIndex = 1, size(charArray)
      if (trim(charArray(arrayIndex)) == trim(value)) then
        numFound = numFound + 1
        ! return the first location found
        if (numFound == 1) location = arrayIndex
      end if
    end do LOOP

    ! give warning if more than 1 found
    if (numFound > 1) then
      write(*,*) 'utl_findloc_char: found multiple locations of ', trim(value)
      write(*,*) 'utl_findloc_char: number locations found =  ', numFound
    end if

    ! return zero if not found
    if (numFound == 0) then
      location = 0
    end if

  end function utl_findloc_char

  !--------------------------------------------------------------------------
  ! utl_findloc_int
  !--------------------------------------------------------------------------
  function utl_findloc_int(intArray, intValue) result(location)
    !
    !:Purpose: A modified version of the fortran function `findloc`.
    !          If multiple matches are found in the array, a warning
    !          message is printed to the listing.
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: intArray(:)
    integer, intent(in) :: intValue
    ! Result:
    integer             :: location

    ! Locals:
    integer :: numFound, arrayIndex

    numFound = 0
    LOOP: do arrayIndex = 1, size(intArray)
      if (intArray(arrayIndex) == intValue) then
        numFound = numFound + 1
        ! return the first location found
        if (numFound == 1) location = arrayIndex
      end if
    end do LOOP

    ! give warning if more than 1 found
    if (numFound > 1) then
      write(*,*) 'utl_findloc_int: found multiple locations of ', intValue
      write(*,*) 'utl_findloc_int: number locations found =  ', numFound
    end if

    ! return zero if not found
    if (numFound == 0) then
      location = 0
    end if

  end function utl_findloc_int

  !--------------------------------------------------------------------------
  ! utl_findlocs_char
  !--------------------------------------------------------------------------
  subroutine utl_findlocs_char(charArray, charValue, locations)
    !
    !:Purpose: A modified version of the fortran function `findloc`.
    !          Returns an array of all matches found in the array.
    !
    implicit none

    ! Arguments:
    character(len=*),     intent(in)  :: charArray(:)
    character(len=*),     intent(in)  :: charValue
    integer, allocatable, intent(out) :: locations(:)

    ! Locals:
    integer :: numFound, arrayIndex

    if (allocated(locations)) deallocate(locations)

    ! count number of matches found
    numFound = 0
    do arrayIndex = 1, size(charArray)
      if (trim(charArray(arrayIndex)) == trim(charValue)) numFound = numFound + 1
    end do

    if (numFound > 0) then

      ! return all found locations
      allocate(locations(numFound))
      numFound = 0
      do arrayIndex = 1, size(charArray)
        if (trim(charArray(arrayIndex)) == trim(charValue)) then
          numFound = numFound + 1
          locations(numFound) = arrayIndex
        end if
      end do

    else

      ! return zero if not found
      allocate(locations(1))
      locations(1) = 0

    end if

  end subroutine utl_findlocs_char

  !--------------------------------------------------------------------------
  ! utl_findlocs_int
  !--------------------------------------------------------------------------
  subroutine utl_findlocs_int(intArray, intValue, locations)
    !
    !:Purpose: A modified version of the fortran function `findloc`.
    !          Returns an array of all matches found in the array.
    !
    implicit none

    ! Arguments:
    integer,              intent(in)  :: intArray(:)
    integer,              intent(in)  :: intValue
    integer, allocatable, intent(out) :: locations(:)

    ! Locals:
    integer :: numFound, arrayIndex

    if (allocated(locations)) deallocate(locations)

    ! count number of matches found
    numFound = 0
    do arrayIndex = 1, size(intArray)
      if (intArray(arrayIndex) == intValue) numFound = numFound + 1
    end do

    if (numFound > 0) then

      ! return all found locations
      allocate(locations(numFound))
      numFound = 0
      do arrayIndex = 1, size(intArray)
        if (intArray(arrayIndex) == intValue) then
          numFound = numFound + 1
          locations(numFound) = arrayIndex
        end if
      end do

    else

      ! return zero if not found
      allocate(locations(1))
      locations(1) = 0

    end if

  end subroutine utl_findlocs_int

  !--------------------------------------------------------------------------
  ! utl_randomOrderInt
  !--------------------------------------------------------------------------
  subroutine utl_randomOrderInt(intArray,randomSeed)
    !
    !:Purpose: Randomly shuffle the order of the integer array elements.
    !
    implicit none

    ! Arguments:
    integer, intent(inout) :: intArray(:)
    integer, intent(in)    :: randomSeed

    ! Locals:
    integer              :: arraySize, arrayIndex, arrayIndexMin
    integer, allocatable :: intArrayOut(:)
    real(8), allocatable :: realRandomArray(:)

    arraySize = size(intArray)
    allocate(realRandomArray(arraySize))
    allocate(intArrayOut(arraySize))

    call rng_setup(randomSeed, beSilent_opt=.true.)
    do arrayIndex = 1, arraySize
      realRandomArray(arrayIndex) = rng_uniform()
    end do

    do arrayIndex = 1, arraySize
      arrayIndexMin = minloc(realRandomArray,dim=1)
      realRandomArray(arrayIndexMin) = huge(1.0D0)
      intArrayOut(arrayIndex) = intArray(arrayIndexMin)
    end do

    intArray(:) = intArrayOut(:)

    deallocate(intArrayOut)
    deallocate(realRandomArray)

  end subroutine utl_randomOrderInt

  !--------------------------------------------------------------------------
  ! utl_tmg_start
  !--------------------------------------------------------------------------
  subroutine utl_tmg_start(blockIndex, blockLabel)
    !
    !:Purpose: Wrapper for rpnlib subroutine tmg_start
    !
    implicit none

    ! Arguments:
    integer,          intent(in) :: blockIndex
    character(len=*), intent(in) :: blockLabel

    ! Locals:
    integer            :: labelLength
    integer, parameter :: labelPaddedLength = 40
    character(len=labelPaddedLength) :: blockLabelPadded

    ! only the first thread does the timing
    if (omp_get_thread_num() > 0) return

    blockLabelPadded = '........................................'
    labelLength = min(len_trim(blockLabel), labelPaddedLength)
    blockLabelPadded(1:labelLength) = blockLabel(1:labelLength)

    call tmg_start(blockIndex, blockLabelPadded)

  end subroutine utl_tmg_start

  !--------------------------------------------------------------------------
  ! utl_tmg_stop
  !--------------------------------------------------------------------------
  subroutine utl_tmg_stop(blockIndex)
    !
    !:Purpose: Wrapper for rpnlib subroutine tmg_stop
    !
    implicit none

    ! Arguments:
    integer,          intent(in) :: blockIndex

    ! only the first thread does the timing
    if (omp_get_thread_num() > 0) return

    call tmg_stop(blockIndex)

  end subroutine utl_tmg_stop

  !--------------------------------------------------------------------------
  ! utl_medianIndex
  !--------------------------------------------------------------------------
  function utl_medianIndex(inputVector) result(medianIndex)
    !
    !:Purpose: to find the median index of an input vector
    !
    implicit none

    ! Arguments:
    real(4), intent(in) :: inputVector(:)
    ! Result:
    integer             :: medianIndex

    ! Locals:
    integer :: vectorIndex, vectorDim
    logical :: maskVector(size(inputVector))
    real(4) :: sortedArray(size(inputVector))
    real(4) :: median

    vectorDim = size(inputVector)

    ! sorting array:
    maskVector(:) = .true.
    do vectorIndex = 1, vectorDim
      sortedArray(vectorIndex) = minval(inputVector, maskVector)
      maskVector(minloc(inputVector, maskVector)) = .false.
    end do

    if (mod(size(inputVector), 2) == 0) then
      median = sortedArray(vectorDim / 2)
    else
      median = sortedArray((vectorDim + 1) / 2)
    end if

    medianIndex =  MPC_missingValue_INT
    do vectorIndex = 1, vectorDim
      if ( utl_isEqual(inputVector(vectorIndex),median) ) then
        medianIndex = vectorIndex
        exit
      end if
    end do

  end function utl_medianIndex

  !--------------------------------------------------------------------------
  ! utl_checkNetCDFstatus
  !--------------------------------------------------------------------------
  subroutine utl_checkNetCDFstatus(status, subroutineName_opt, ncCommand_opt, addInfo_opt)
    !
    ! :Purpose: to check the status of a netCDF command. This subroutine may be usefull
    !           when adding a new netCDF realted code.
    !
    implicit none

    ! Arguments:
    integer                   , intent(in) :: status
    character(len=*), optional, intent(in) :: subroutineName_opt
    character(len=*), optional, intent(in) :: ncCommand_opt
    character(len=*), optional, intent(in) :: addInfo_opt

    ! Locals:
    character(len=256) :: errorMessage

    if(status /= nf90_noerr) then

      errorMessage = 'netCDF error'
      if (present(subroutineName_opt)) then
        errorMessage = trim(subroutineName_opt)//': '//trim(errorMessage)
      end if
      if (present(ncCommand_opt)) then
        errorMessage = trim(errorMessage)//' when calling: '//trim(ncCommand_opt)
      end if
      if (present(addInfo_opt)) then
        errorMessage=trim(errorMessage)//' for '//trim(addInfo_opt)
      end if

      write(*,*) 'nf90 error status: ', trim(nf90_strerror(status))
      call utl_abort(trim(errorMessage))

    end if

  end subroutine utl_checkNetCDFstatus

  !--------------------------------------------------------------------------
  ! utl_varPresentInNetcdfFile
  !--------------------------------------------------------------------------
  function utl_varPresentInNetcdfFile(varName, fileName) result(variableFound)
    !
    ! :Purpose: to verify if the given varName is present in netCDF file.
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in)  :: varName  ! variable name
    character(len=*), intent(in)  :: fileName ! netCDF filename

    ! Result:
    logical :: variableFound

    ! Locals:
    integer :: ncid, varID, ierr

    ! Open the template file
    call utl_checkNetCDFstatus(nf90_open(trim(fileName), nf90_nowrite, ncid))

    ! Inquire variable name
    ierr = nf90_inq_varid(ncid, trim(varName), varID)
    if (ierr == nf90_noerr) then
      write(*,*) 'utl_varPresentInNetcdfFile: variable: ', &
                 trim(varName), ' is found in file: ', trim(fileName)
      variableFound = .true.
    else
      write(*,*) 'utl_varPresentInNetcdfFile: variable: ', &
                   trim(varName), ' is missing from file: ', trim(fileName)
      variableFound = .false.
    end if

    ! Close the file
    call utl_checkNetCDFstatus(nf90_close(ncid))

  end function utl_varPresentInNetcdfFile

  !--------------------------------------------------------------------------
  ! utl_cosDegrees_real4
  !--------------------------------------------------------------------------
  function utl_cosDegrees_real4(degrees) result(cosinus)
    !
    ! :Purpose: Computes the cosinus of the angle where the angle is
    !           specified in degrees.
    !           All arguments are in single precision floating point numbers, real(4).
    !
    implicit none

    ! Arguments:
    real(4), intent(in)  :: degrees ! angle in degrees
    ! Result:
    real(4) :: cosinus

    ! Locals:
    real(4) :: radians

    radians = MPC_RADIANS_PER_DEGREE_R4 * degrees

    cosinus = cos(radians)
  end function utl_cosDegrees_real4

  !--------------------------------------------------------------------------
  ! utl_cosDegrees_real8
  !--------------------------------------------------------------------------
  function utl_cosDegrees_real8(degrees) result(cosinus)
    !
    ! :Purpose: Computes the cosinus of the angle where the angle is
    !           specified in degrees.
    !           All arguments are in double precision floating point numbers, real(8).
    !
    implicit none

    ! Arguments:
    real(8), intent(in)  :: degrees ! angle in degrees
    ! Result:
    real(8) :: cosinus

    ! Locals:
    real(8) :: radians

    radians = MPC_RADIANS_PER_DEGREE_R8 * degrees

    cosinus = cos(radians)
  end function utl_cosDegrees_real8

  !--------------------------------------------------------------------------
  ! utl_isEqual_real4
  !--------------------------------------------------------------------------
  function utl_isEqual_real4(firstValue, secondValue) result(areTheyEqual)
    !
    ! :Purpose: Checks if two real(4) values are equal according to the machine precision
    !           All arguments are in single precision floating point numbers, real(4).
    !
    implicit none

    ! Arguments:
    real(4), intent(in) :: firstValue  ! First  real(4) value to compare with the second value
    real(4), intent(in) :: secondValue ! Second real(4) value to compare with the first value
    ! Result:
    logical :: areTheyEqual

    ! tiny(X) returns the smallest positive (non zero) number in the model of the type of X.
    areTheyEqual = abs(firstValue-secondValue) < tiny(firstValue)

  end function utl_isEqual_real4

  !--------------------------------------------------------------------------
  ! utl_isEqual_real8
  !--------------------------------------------------------------------------
  function utl_isEqual_real8(firstValue, secondValue) result(areTheyEqual)
    !
    ! :Purpose: Checks if two real(8) values are equal according to the machine precision
    !           All arguments are double precision floating point numbers, real(8).
    !
    implicit none

    ! Arguments:
    real(8), intent(in) :: firstValue  ! First  real(8) value to compare with the second value
    real(8), intent(in) :: secondValue ! Second real(8) value to compare with the first value
    ! Result:
    logical :: areTheyEqual

    ! tiny(X) returns the smallest positive (non zero) number in the model of the type of X.
    areTheyEqual = abs(firstValue-secondValue) < tiny(firstValue)

  end function utl_isEqual_real8

  !--------------------------------------------------------------------------
  ! utl_isEqual_real4Arrays
  !--------------------------------------------------------------------------
  function utl_isEqual_real4Arrays(firstArray, secondArray) result(areTheyEqual)
    !
    ! :Purpose: Checks if two arrays of real(4) values are pair-wise
    !           equal according to the machine precision.
    !           It returns an array of logical values comparing the two arrays index by index.
    !           If the arrays are not of the same size, it will return an array of '.false.'
    !           of the same size as 'firstArray'.
    !           All arguments are single precision floating point numbers, real(4).
    !
    implicit none

    ! Arguments:
    real(4), intent(in) :: firstArray(:)  ! First  array of real(4) values to compare with the second value
    real(4), intent(in) :: secondArray(:) ! Second array of real(4) values to compare with the first value
    ! Result:
    logical :: areTheyEqual(size(firstArray)) ! Array of logical values

    ! Locals:
    integer :: arrayIndex

    ! If the arrays does not have the same sizes, they can't be equal
    if ( size(firstArray) /= size(secondArray) ) then
      areTheyEqual(:) = .false.
      return
    end if

    do arrayIndex = 1, size(firstArray)
      areTheyEqual(arrayIndex) = utl_isEqual(firstArray(arrayIndex), secondArray(arrayIndex))
    end do

    ! If we didn't catch any different value in the array, then the
    ! arrays are equal.
    areTheyEqual = .true.

  end function utl_isEqual_real4Arrays

  !--------------------------------------------------------------------------
  ! utl_isEqual_real8Arrays
  !--------------------------------------------------------------------------
  function utl_isEqual_real8Arrays(firstArray, secondArray) result(areTheyEqual)
    !
    ! :Purpose: Checks if two arrays of real(8) values are pair-wise
    !           equal according to the machine precision.
    !           It returns an array of logical values comparing the two arrays index by index.
    !           If the arrays are not of the same size, it will return an array of '.false.'
    !           of the same size as 'firstArray'.
    !           All arguments are double precision floating point numbers, real(8).
    !
    implicit none

    ! Arguments:
    real(8), intent(in) :: firstArray(:)  ! First  array of real(8) values to compare with the second value
    real(8), intent(in) :: secondArray(:) ! Second array of real(8) values to compare with the first value
    ! Result:
    logical :: areTheyEqual(size(firstArray))

    ! Locals:
    integer :: arrayIndex

    ! If the arrays does not have the same sizes, they can't be equal
    if ( size(firstArray) /= size(secondArray) ) then
      areTheyEqual(:) = .false.
      return
    end if

    do arrayIndex = 1, size(firstArray)
      areTheyEqual(arrayIndex) = utl_isEqual(firstArray(arrayIndex), secondArray(arrayIndex))
    end do

  end function utl_isEqual_real8Arrays

  !--------------------------------------------------------------------------
  ! utl_isEqual_real4ArraysScalar
  !--------------------------------------------------------------------------
  function utl_isEqual_real4ArraysScalar(array, scalar) result(areTheyEqual)
    !
    ! :Purpose: Checks if values of an array of real(4) are equal to a scalar
    !           equal according to the machine precision.
    !           It returns an array of logical values comparing each value of the array to the scalar.
    !           All arguments are single precision floating point numbers, real(4).
    !
    implicit none

    ! Arguments:
    real(4), intent(in) :: array(:) ! Array of real(4) values to compare with the second value
    real(4), intent(in) :: scalar   ! A real(4)' scalar to compare with the values in the array
    ! Result:
    logical :: areTheyEqual(size(array))

    ! Locals:
    integer :: arrayIndex

    do arrayIndex = 1, size(array)
      areTheyEqual(arrayIndex) = utl_isEqual(array(arrayIndex), scalar)
    end do

  end function utl_isEqual_real4ArraysScalar

  !--------------------------------------------------------------------------
  ! utl_isEqual_real8ArraysScalar
  !--------------------------------------------------------------------------
  function utl_isEqual_real8ArraysScalar(array, scalar) result(areTheyEqual)
    !
    ! :Purpose: Checks if values of an array of real(8) are equal to a scalar
    !           equal according to the machine precision.
    !           It returns an array of logical values comparing each value of the array to the scalar.
    !           All arguments are double precision floating point numbers, real(8).
    !
    implicit none

    ! Arguments:
    real(8), intent(in) :: array(:) ! Array of real(8) values to compare with the second value
    real(8), intent(in) :: scalar   ! A real(8)' scalar to compare with the values in the array
    ! Result:
    logical :: areTheyEqual(size(array))

    ! Locals:
    integer :: arrayIndex

    do arrayIndex = 1, size(array)
      areTheyEqual(arrayIndex) = utl_isEqual(array(arrayIndex), scalar)
    end do

  end function utl_isEqual_real8ArraysScalar

end module utilities_mod
