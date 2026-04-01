
module climatologies_mod
  ! MODULE climatologies_mod (prefix='clm' category='5. Observation operators')
  !
  ! :Purpose: Access to climatologies
  !
  use midasMpi_mod
  use bufr_mod
  use mathPhysConstants_mod
  use physicsFunctions_mod
  use presProfileOperators_mod
  use obsSubSpaceData_mod
  use utilities_mod
  use varNameList_mod
  use timeCoord_mod

  implicit none
  save
  private

  ! Public procedures
  public :: clm_readFields, clm_setColumn, clm_getColumn

  ! Public derived type
  public :: struct_clm

  !-------------------------------------------------------------------------
  ! Declaration of structures and parameters

  ! Arrays containing input climatology/reference fields and fields
  ! interpolated to obs locations

  type :: struct_clm

    !  Structure for storing reference (climatological) fields)
    real(8), allocatable :: field(:,:,:) ! Gridded 3D field (lon,lat,vlev) or 2D field (1,lat,vlev)
    real(8), allocatable :: lat(:)       ! lat grid in radians
    real(8), allocatable :: lon(:)       ! lon grid in radians
    real(8), allocatable :: vlev(:)      ! vertical levels
    integer              :: nlev         ! number of vertical levels
    integer              :: nlon         ! number of longitudes
    integer              :: nlat         ! number of latitudes
    integer              :: vertCoordKind! Index of vertical coordinate type. Defintion may vary according to source.
                                         ! For fields read from RPN files and use of convip:
                                         !     0: P is in height [m] (metres) with respect to sea level
                                         !     1: P is in stddev [sg] (0.0 -111.0)
                                         !     2: P is in pressure [mb] (millibars)
                                         !     3: P is in an arbitrary code
                                         !     4: P is in height [M] (metres) with respect to ground level
                                         !     5: P is in hybrid coordinates [hy]
                                         !     6: P is in theta [th]
                                         ! For use with obs
  end type struct_clm

  ! clm_readFields climatology fields and parameters needed externally to clm_readFields
  type(struct_clm), allocatable :: climatFields(:,:) ! First dimension will be 0:maxNumConstituents
  integer :: maxNumTypes  ! Max of climatFields second dimension
  integer :: maxNumFields ! Number of variables with climatologies
  integer, parameter :: maxNumConstituents = BUFR_NECH_maxValue
  logical :: nearestNeighbourInterp(0:maxNumConstituents)

contains

  !--------------------------------------------------------------------------
  ! clm_readFields
  !--------------------------------------------------------------------------
  subroutine clm_readFields(modelName_opt)
    !
    !:Purpose:  To read climatrology (reference) fields as directed by input
    !
    !:Comments:
    !      - Fields are provided in RPN/fst files
    !      - Fields set/assumed to be of the same units as those of the
    !        corresponding input trial fields (micrograms/kg in most cases)
    !      - Case sourceIndex=2 for numFields=2 to combine separate strato
    !        and tropo climatologies is missing the setting of etiket below.
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in), optional :: modelName_opt ! Model name

    ! Locals:
    character(len=256) :: fname
    character(len=4) :: varName, climVarName
    character(len=12) :: etiket
    character(len=12), parameter :: climFields = 'climatFields'
    character(len=11), parameter :: ozoneBaseline = 'ozoneclim98'
    integer :: varIndex, constituentId, SourceIndex
    integer :: ijour, imonth, iday, itime, iyear
    real(8) :: day, scaleFactor
    integer :: datestamp
    logical :: initialized = .false.
    logical :: timeInterp
    integer :: ierr
    logical :: fileExists
    integer :: ni, nj, nkeys, vertCoordKind, loopIndex
    real(8), allocatable :: array1(:,:,:), array2(:,:,:), lvls(:), xlat(:), xlong(:)
    integer, parameter :: baselineLevelNum = 28
    real(8) :: climatLevelsDefault(baselineLevelNum)  !  Fortuin-Kelder 1998 pressure levels
    data climatLevelsDefault / &
        0.010d0, 0.015d0, 0.022d0, 0.032d0, 0.046d0, 0.068d0, 0.100d0,  &
        0.150d0, 0.200d0, 0.300d0, 0.500d0, 1.000d0, 2.000d0, 3.000d0,  &
        5.000d0, 7.000d0, 10.00d0, 20.00d0, 30.00d0, 50.00d0, 70.00d0,  &
        100.0d0, 150.0d0, 200.0d0, 300.0d0, 500.0d0, 700.0d0, 1000.d0 /

    ! Namelist variables (local):
    character(len=4) :: requiredConstituents(maxNumConstituents+1) ! List of constituent NOMVARs for which climatology fields are required
    character(len=4) :: inputClimVarNames(maxNumConstituents+1)    ! List of input climatology NOMVARs associated to requiredConstituents
    character(len=256) :: climatSourceFileDefault                  ! Default climatology source file
    character(len=256) :: climatSourceFiles(maxNumConstituents+1)  ! Climatology source file names
    integer :: fieldDimension(maxNumConstituents+1)                ! Dimension of input field
    integer :: numFields(maxNumConstituents+1)                     ! Number of climatology fields per constituent
    logical :: timeInterpolation                                   ! Indicating if interpolation between mid-months
    logical :: nearNeighbourInterp(maxNumConstituents+1)           ! Nearest neighbour or linear interpolation
    real(8) :: climatLevels(maxNumConstituents+1,100)              ! Pressure levels (hPa) needed if fieldDimension=1 or 2.

    ! NAMCLIMATOLOGY namelist parameters
    namelist /NAMCLIMATOLOGY/ climatSourceFileDefault
    namelist /NAMCLIMATOLOGY/ climatSourceFiles
    namelist /NAMCLIMATOLOGY/ requiredConstituents
    namelist /NAMCLIMATOLOGY/ inputClimVarNames
    namelist /NAMCLIMATOLOGY/ nearNeighbourInterp
    namelist /NAMCLIMATOLOGY/ numFields
    namelist /NAMCLIMATOLOGY/ timeInterpolation
    namelist /NAMCLIMATOLOGY/ fieldDimension
    namelist /NAMCLIMATOLOGY/ climatLevels

    if (initialized) return

    ! Set defaults
    requiredConstituents(:) = ''
    requiredConstituents(1) = 'O3L'
    inputClimVarNames(:) = ''
    climatSourceFileDefault = climFields
    climatSourceFiles(:) = ''
    climatSourceFiles(1) = ozoneBaseline
    numFields(:) = 0
    numFields(1) = 1
    fieldDimension(:) = 3
    fieldDimension(1) = 2
    nearNeighbourInterp(:) = .false.
    nearNeighbourInterp(1) = .true.
    timeInterpolation = .false.
    climatLevels(:,:) = -1.0

    ! Check for presence or read NAMCLIMATOLOGY
    if ( .not. utl_isNamelistPresent('NAMCLIMATOLOGY','./flnml') ) then
      write(*,*)
      write(*,*) 'clm_readFields: Namelist block NAMCLIMATOLOGY is missing. Using defaults'
      write(*,*)
    else
      call utl_tmg_start(181,'low-level--readNML')
      read(utl_flnml, nml=namclimatology, iostat=ierr)
      if (ierr /= 0) call utl_abort('clm_readFields: Error reading namelist NAMCLIMATOLOGY')
      call utl_tmg_stop(181)
    end if
    if (mmpi_myid == 0) write(*,nml=namclimatology)

    ! Set dimensions
    maxNumFields = 1
    maxNumTypes = 1
    do while (trim(requiredConstituents(maxNumFields)) /= '')
      ! The call to 'vnl_varListIndex3d' will abort if it does not find the variable
      ! so we just need to call it from here.
      loopIndex = vnl_varListIndex3d(trim(requiredConstituents(maxNumFields)))
      if (numFields(maxNumFields) < 1) then
        numFields(maxNumFields) = 1
      else if (numFields(maxNumFields) > maxNumTypes) then
        maxNumTypes = numFields(maxNumFields)
      end if
      maxNumFields = maxNumFields + 1
    end do
    maxNumFields = maxNumFields - 1

    if ( maxNumTypes > 2 ) then
      call utl_abort('clm_readFields: Allowed max number of fields per constituent is 2')
    end if

    allocate(climatFields(0:maxNumConstituents,maxNumTypes))

    ! Initialize dimensions

    climatFields(:,:)%nlon = 0
    climatFields(:,:)%nlat = 0
    climatFields(:,:)%nlev = 0

    ! Initialization of work parameters

    datestamp = tim_getDateStamp()
    call tim_dateStampToYYYYMMDDHH(datestamp,ijour,itime)
    iyear = ijour/10000
    imonth = MOD(ijour/100,100)
    iday = MOD(ijour,100)
    day = iday + itime*1.0D-8
    if (day > 15.) then
      day = day - 15.0
    else
      day = day + 15.0
    end if

    ! Get needed fields for each varIndex
    do varIndex = 1, maxNumFields
      ! Identify input file
      if (trim(climatSourceFiles(varIndex)) /= '') then
        inquire(file=trim(climatSourceFiles(varIndex)),exist=fileExists)
        if (fileExists) then
          fname = trim(climatSourceFiles(varIndex))
        else
          write(*,*)
          write(*,*) 'clm_readFields: Climatologies file does not exist. Filename ', &
              trim(climatSourceFiles(varIndex))
          write(*,*) 'Specifying file ',trim(climatSourceFileDefault)
          write(*,*)
          fname = trim(climatSourceFileDefault)
          inquire(file=trim(fname),exist=fileExists)
          if (.not.fileExists) then
            call utl_abort('clm_readFields: Climatologies file ' // trim(fname) // ' is missing.')
          end if
        end if
      else
        write(*,*)
        write(*,*) 'clm_readFields: Climatologies file name not set.'
        write(*,*) 'Specifying file ', trim(climatSourceFileDefault)
        write(*,*)
        fname = trim(climatSourceFileDefault)
        inquire(file=trim(fname),exist=fileExists)
        if (.not.fileExists) then
          call utl_abort('clm_readFields: Climatologies file ' // trim(fname) // ' is missing.')
        end if
      end if

      varName = trim(requiredConstituents(varIndex))
      if (trim(fname) == ozoneBaseline) then
        if (trim(varName) == &
            vnl_varnameFromVarnum(0,varNumberChm_opt=00,modelName_opt=modelName_opt)) then
          climVarName = 'O3'
          if (fieldDimension(varIndex) /= 2) then
            call utl_abort('clm_readFields: Invalid field dimension for ' // trim(varName) // &
                           ' in ' // trim(fname))
          end if
        else
          call utl_abort('clm_readFields: Invalid climatology file ' // trim(fname) // &
                         ' for ' // trim(varName))
        end if
      else
        climVarName = trim(inputClimVarNames(varIndex))
        if (trim(fname) == climFields) then
          if (trim(varName) == &
              vnl_varnameFromVarnum(0,varNumberChm_opt=00,modelName_opt=modelName_opt)) then
            if (trim(climVarName) == '') climVarName = 'O3CE'
          else if (trim(varName) == &
              vnl_varnameFromVarnum(0,varNumberChm_opt=02,modelName_opt=modelName_opt)) then
            if (trim(climVarName) == '') climVarName = 'CH4C'
          else if (trim(varName) == &
              vnl_varnameFromVarnum(0,varNumberChm_opt=06,modelName_opt=modelName_opt)) then
            if (trim(climVarName) == '') climVarName = 'N2OC'
          else
            if (trim(climVarName) == '') climVarName = varName
          end if
          if (trim(climVarName) == 'O3CE' .or. trim(climVarName) == 'CH4C' .or. &
              trim(climVarName) == 'N2OC') fieldDimension(varIndex) = 3
        else
          if (trim(climVarName) == '') climVarName = varName
        end if
      end if

      constituentId = vnl_varnumFromVarName(varName,'CH')
      if (constituentId /=0 .and. numFields(varIndex) > 1) then
        call utl_abort('clm_readFields: numFields > 1 only allowed for ozone currently')
      end if
      ! Set interpolation approach
      nearestNeighbourInterp(constituentId) = nearNeighbourInterp(varIndex)

      if ( fieldDimension(varIndex) <= 1 .and. numFields(varIndex) > 1 ) then
        call utl_abort('clm_readFields: Separate strato and tropo input ' &
            // 'climatologies only for 2D and 3D fields')
      end if

      timeInterp = timeInterpolation

      do sourceIndex = 1, numFields(varIndex)

        ! Read climatology for specified field

        if ( sourceIndex == 2 ) then
          ! Not currently applied. Standin for eventual use.
          etiket = '            '
          call utl_abort('clm_readFields: numFields=2 case tbc')
        else
          etiket = '            '
        end if

        if ( fieldDimension(varIndex) == 0 ) then

          ! Mixing ratio set later by the constant in climatScaling
          ! Currently assumed to be relevant only for GHGs

          if ( .not.clm_isGHG(trim(varName)) ) call utl_abort('clm_readFields: Needs to be a GHG')

          if (trim(climVarName) == '') climVarName = varName
          allocate(array1(1,1,1),lvls(1),xlat(1),xlong(1))
          vertCoordKind = 2
          ni = 1
          nj = 1
          nkeys = 1
          xlat(1) = 0.0
          xlong(1) = 0.0
          lvls(1) = 0.0d0
          array1(1:ni,1:nj,1:nkeys) = 1.0
          timeInterp =  .false.

        else if ( fieldDimension(varIndex) <= 2) then

          ! Read 1D or 2D climatology for specified month
          call clm_readSub3DFields(imonth,array1)

          ! Set lat and long
          allocate(xlong(ni),xlat(nj),lvls(nkeys))
          xlong(1) = 0.0d0
          if (nj > 1) then
            do loopIndex = 1, nj
              xlat(loopIndex) = (loopIndex-1)*180.0d0/real(nj-1,8) - 90.0d0
            end do
          else
            xlat(1) = 0.0d0
          end if

          ! Set pressures (hPa)
          if (nkeys == baselineLevelNum .and. &
              any(climatLevels(varIndex,1:nkeys) < 0.0) ) then
            ! Default pressure levels
            lvls(1:nkeys) = climatLevelsDefault(1:nkeys) ! Store in hPa for consistency with kind=2
            vertCoordKind = 2
          else
            if (all(climatLevels(varIndex,1:nkeys) > 0.0)) then
              lvls(1:nkeys) = climatLevels(varIndex,1:nkeys)
              vertCoordKind = 2
            else
              call utl_abort('clm_readFields: Missing set up for levels')
            end if
          end if

        else

          call utl_readFstField(trim(fname),trim(climVarName),-1,imonth,-1,etiket, &
              ni,nj,nkeys,array1,xlat_opt=xlat,xlong_opt=xlong,         &
              lvls_opt=lvls,kind_opt=vertCoordKind)

        end if

        climatFields(constituentId,sourceIndex)%nlon=ni
        climatFields(constituentId,sourceIndex)%nlat=nj
        climatFields(constituentId,sourceIndex)%nlev=nkeys
        climatFields(constituentId,sourceIndex)%vertCoordKind=vertCoordKind

        allocate(climatFields(constituentId,sourceIndex)%field(ni,nj,nkeys))
        allocate(climatFields(constituentId,sourceIndex)%vlev(nkeys))
        allocate(climatFields(constituentId,sourceIndex)%lon(ni))
        allocate(climatFields(constituentId,sourceIndex)%lat(nj))

        climatFields(constituentId,sourceIndex)%lat(1:nj) = xlat(1:nj)
        climatFields(constituentId,sourceIndex)%lon(1:ni) = xlong(1:ni)
        where (climatFields(constituentId,sourceIndex)%lon(1:ni) < 0.0)
          climatFields(constituentId,sourceIndex)%lon(1:ni) = 2.0*MPC_PI_R8 &
              + climatFields(constituentId,sourceIndex)%lon(1:ni)
        end where
        climatFields(constituentId,sourceIndex)%vlev(1:nkeys) = lvls(1:nkeys)

        if (.not.timeInterp) then
          climatFields(constituentId,sourceIndex)%field(:,:,:) = array1(:,:,:)
        else

          ! Following for interpolation as a function of days from mid-months.
          ! Assumes same grid and level
          if (iday > 15) then
            if ( fieldDimension(varIndex) <= 2 ) then
              if (imonth == 12) then
                call clm_readSub3DFields(1,array2)
              else
                call clm_readSub3DFields(imonth+1,array2)
              end if
            else
              if (imonth == 12) then
                call utl_readFstField(trim(fname),trim(climVarName),-1,1,-1,etiket, &
                    ni,nj,nkeys,array2)
              else
                call utl_readFstField(trim(fname),trim(climVarName),-1,imonth+1,-1, &
                    etiket,ni,nj,nkeys,array2)
              end if
            end if

            ! Linearly interpolate in time
            ! (approximation - assumes 30 day months)

            climatFields(constituentId,sourceIndex)%field(:,:,:) = &
                (array1(:,:,:)*(30.0-day) + array2(:,:,:)*day)/30.0

          else if (iday <= 15) then
            if ( fieldDimension(varIndex) <= 2 ) then
              if (imonth == 1) then
                call clm_readSub3DFields(12,array2)
              else
                call clm_readSub3DFields(imonth-1,array2)
              end if
            else
              if (imonth == 1) then
                call utl_readFstField(trim(fname),trim(climVarName),-1,12,-1,etiket, &
                    ni,nj,nkeys,array2)
              else
                call utl_readFstField(trim(fname),trim(climVarName),-1,imonth-1,-1, &
                    etiket,ni,nj,nkeys,array2)
              end if
            end if

            ! Linearly interpolate in time
            ! (approximation - applies 30 day months)

            climatFields(constituentId,sourceIndex)%field(:,:,:) = &
                (array2(:,:,:)*(30.0-day)+array1(:,:,:)*day)/30.0
          end if
        end if

        if (clm_isGHG(trim(varName)) ) then
          ! Apply units scaling factor to convert from volumetric units (e.g. ppmv or ppbv)
          ! to micrograms/kg as needed. Assumes in micrograms/kg otherwise.
          scaleFactor = clm_getUnitsScaling(varName,iyear)
          climatFields(constituentId,sourceIndex)%field(:,:,:)= &
              climatFields(constituentId,sourceIndex)%field(:,:,:) &
              *scaleFactor*vnl_varMassFromVarName(trim(varName)) &
              /MPC_MOLAR_MASS_DRY_AIR_R8
        else if (trim(climVarName) == 'O3CE') then
          ! Convert kg/kg to ug/kg for ozone
          climatFields(constituentId,sourceIndex)%field(:,:,:)= 1.0D9 &
              *climatFields(constituentId,sourceIndex)%field(:,:,:)
        end if

        if (allocated(array1)) deallocate(array1,lvls,xlat,xlong)
        if (allocated(array2)) deallocate(array2)

      end do
    end do

    initialized = .true.

    contains

    !--------------------------------------------------------------------------
    ! clm_readSub3DFields
    !--------------------------------------------------------------------------
    subroutine clm_readSub3DFields(month,outarray)
      !
      !:Purpose:  Read 1D or 2D climatology (reference) field
      !
      implicit none

      ! Arguments:
      integer, intent(in)  :: month                        ! Month for input field
      real(8), allocatable, intent(out) :: outarray(:,:,:) ! Final output field

      ! Locals:
      real(8), allocatable :: inarray(:,:,:)               ! Input field

      call utl_readFstField(trim(fname),climVarName,-1,-1,month,etiket, &
                            ni,nj,nkeys,inarray)

      if (ni == 1) then
        ! nj for latitudes and nkeys for levels
        allocate(outarray(ni,nj,nkeys))
        outarray(1,1:nj,1:nkeys) =  inarray(1,1:nj,1:nkeys)
      else if (nkeys == 1) then
        if (nj == 1) then
          nkeys = ni
          ni = 1
          allocate(outarray(ni,nj,nkeys))
          outarray(1,1,1:nkeys) = inarray(1:nkeys,1,1)
        else
          nkeys = nj
          nj = ni
          ni = 1
          allocate(outarray(ni,nj,nkeys))
          if (trim(fname) == ozoneBaseline) then
            ! Convert from ppmv to micrograms/kg
            outarray(1,1:nj,1:nkeys) = inarray(1:nj,1:nkeys,1) * &
              1.0d3 * MPC_MOLAR_MASS_O3_R8 / MPC_MOLAR_MASS_DRY_AIR_R8 ! Same as above
          else
            outarray(1,1:nj,1:nkeys) = inarray(1:nj,1:nkeys,1)
          end if
        end if
      end if
      deallocate(inarray)

    end subroutine clm_readSub3DFields

  end subroutine clm_readFields

  !--------------------------------------------------------------------------
  ! clm_setColumn
  !--------------------------------------------------------------------------
  subroutine clm_setColumn(numModelLevs,modelPressLevs, &
                           modelHeightLevs,obsLat,obsLong,obsIndex, &
                           maxsize,constituentId,tt_opt,hu_opt,     &
                           climatProfileSet_opt,climatProfile_opt)

    !:Purpose: To determine and store a column profile at obs location.
    !          Optional storage as part of a cumulative profile set for a
    !          specific variable.
    !
    !:Comments:
    !
    !          - tt_opt and hu_opt required only when maxNumTypes = 2.
    !

    implicit none

    ! Arguments:
    integer,                            intent(in)    :: obsIndex       ! Unique measurement identifier
    integer,                            intent(in)    :: numModelLevs   ! Number of model levels
    integer,                            intent(in)    :: maxsize        ! Max number of obs for which climatProfileSet will be used
    real(8),                            intent(in)    :: modelPressLevs(numModelLevs)  ! Model pressure array (Pa)
    real(8),                            intent(in)    :: modelHeightLevs(numModelLevs) ! Model height (m)
    real(8),                            intent(in)    :: obsLat         ! Latitude (degrees)
    real(8),                            intent(in)    :: obsLong        ! Longitude (degrees)
    integer,                            intent(in)    :: constituentId  ! Constituent id
    real(8),                  optional, intent(in)    :: tt_opt(:) ! Model temperature (Kelvin)
    real(8),                  optional, intent(in)    :: hu_opt(:) ! Model specific humidity
    type(struct_oss_obsdata), optional, intent(inout) :: climatProfileSet_opt ! I/O profile set
    real(8),                  optional, intent(out)   :: climatProfile_opt(numModelLevs)  ! Output profile

    ! Locals;
    integer :: level, start
    real(8) :: tropo_press, refprof(numModelLevs), refprof2(numModelLevs), dt
    real(8), allocatable :: pressrefin(:)
    logical, allocatable :: success(:)

    if (.not.present(climatProfileSet_opt) .and. .not.present(climatProfile_opt)) then
      call utl_abort('clm_setColumn: Missing output array argument.')
    end if

    if (.not.allocated(climatFields)) then
      call utl_abort('clm_setColumn: Climatologies not set.')
    end if
    if (climatFields(constituentId,1)%nlat == 0 .or. &
        climatFields(constituentId,1)%nlon == 0 .or. &
        climatFields(constituentId,1)%nlev == 0 ) return

    if (climatFields(constituentId,1)%nlat == 1 .and. &
        climatFields(constituentId,1)%nlon == 1 .and. &
        climatFields(constituentId,1)%nlev == 1) then

      refprof(:) = climatFields(constituentId,1)%field(1,1,1)

      ! The following scaling accounts for mixing ratios being values in dry air.
      ! This would be needed when integrating vertically in pressure.
      ! It should be applied in the vertical integration operator
      ! and not in the setting of mixing ratios here (hence commented below).

      ! Account for relative difference in pressure of dry air and humid air.

      !if ( .not. present(hu_opt) ) then
      !  call utl_abort('clm_setColumn: Missing HU for determining ' // &
      !               'mixing ratio in dry air ')
      !end if
      !if ( any(hu_opt <= 0.0d0) ) then
      !  call utl_abort('clm_setColumn: Invalid HU for determining ' // &
      !                 'mixing ratio in dry air density')
      !end if

      !do level=1,numModelLevs
      !  refprof(level)=1.0D0-phf_FOEFQ8(hu_opt(level),modelPressLevs(level))&
      !                  /modelPressLevs(level))
      !end do
      !refprof(:)=climatFields(constituentId,1)%field(1,1,1)*refprof(:)

    else

      ! Spatial interpolation

      ! Set vertical levels of reference.
      ! Convert to pressure coordinate if needed.

      if (allocated(pressrefin)) deallocate(pressrefin)
      allocate(pressrefin(climatFields(constituentId,1)%nlev))
      pressrefin(:) = climatFields(constituentId,1)%vlev(1:climatFields(constituentId,1)%nlev)

      if (allocated(success)) deallocate(success)
      allocate(success(climatFields(constituentId,1)%nlev))
      success(:) = .true.

      if (climatFields(constituentId,1)%vertCoordKind == 2) then
        pressrefin(:) = pressrefin(:) * MPC_PA_PER_MBAR_R8
      else if (climatFields(constituentId,1)%vertCoordKind == 0) then
        where (pressrefin < modelHeightLevs(numModelLevs))
          pressrefin = modelHeightLevs(numModelLevs)
        end where
        pressrefin(:) = phf_convertZtoPressure(pressrefin,modelHeightLevs,modelPressLevs, &
            climatFields(constituentId,1)%nlev,numModelLevs, &
            obsLat*MPC_RADIANS_PER_DEGREE_R8,success)
      else if (climatFields(constituentId,1)%vertCoordKind == 4) then
        pressrefin(:) = pressrefin(:) + modelHeightLevs(numModelLevs)
        pressrefin(:) = phf_convertZtoPressure(pressrefin,modelHeightLevs,modelPressLevs, &
            climatFields(constituentId,1)%nlev,numModelLevs, &
            obsLat*MPC_RADIANS_PER_DEGREE_R8,success)
      else if (climatFields(constituentId,1)%vertCoordKind == 1) then
        pressrefin(:) = pressrefin(:)*modelPressLevs(numModelLevs) ! Convert from sigma to Pa
      else
        call utl_abort('clm_setColumn: Cannot handle vertical coordinate of kind ' &
            // trim(utl_str(climatFields(constituentId,1)%vertCoordKind)))
      end if

      ! Interpolate to obs lat/long (or lat) location and model level

      ! Latitudes provided in degrees to ensure backward compatibility of results
      ! when the ozone Fortuin-Kelder climatology is used.

      call ppo_getColumn(climatFields(constituentId,1)%field, &
          climatFields(constituentId,1)%nlon,climatFields(constituentId,1)%nlat, &
          climatFields(constituentId,1)%nlev,climatFields(constituentId,1)%lon*MPC_RADIANS_PER_DEGREE_R8, &
          climatFields(constituentId,1)%lat,pressrefin, &
          obsLong*MPC_RADIANS_PER_DEGREE_R8,obsLat,numModelLevs,modelPressLevs, &
          nearestNeighbourInterp(constituentId),refprof)

      if (maxNumTypes == 1) then

        ! No second climatology

      else if (climatFields(constituentId,2)%nlat > 0 .and. &
          climatFields(constituentId,2)%nlon > 0 .and. &
          climatFields(constituentId,2)%nlev > 0) then

        ! Get second climatology field
        ! Assumes to be for replacing/filling the troposphere and
        ! requires at least tt_opt
        ! (i.e. for combining distinct lower and middle atmosphere sources)

        if ( .not. present(tt_opt) ) then
          call utl_abort('clm_setColumn: Missing TT for determining ' // &
              'tropopause pressure')
        end if
        if ( any(tt_opt <= 0.0d0) ) then
          call utl_abort('clm_setColumn: Invalid TT for determining ' // &
              'tropopause pressure')
        end if

        tropo_press = -1.0

        if ( present(hu_opt) .and. present(tt_opt) ) then
          if (all(hu_opt >= 0.0D0)) then
            tropo_press = phf_calcTropopause(numModelLevs,modelPressLevs, &
                tt_opt,modelHeightLevs,hu_opt=hu_opt)
          else
            tropo_press = phf_calcTropopause(numModelLevs,modelPressLevs, &
                tt_opt,modelHeightLevs)
          end if
        else
          tropo_press = phf_calcTropopause(numModelLevs,modelPressLevs,tt_opt,modelHeightLevs)
        end if

        if (tropo_press <= 0) then
          call utl_abort('clm_setColumn: Invalid tropopause level')
        end if
        ! Set vertical levels of reference.
        ! Convert to pressure coordinate if needed

        if (allocated(pressrefin)) deallocate(pressrefin)
        allocate(pressrefin(climatFields(constituentId,2)%nlev))
        pressrefin(:) = climatFields(constituentId,2)%vlev(1:climatFields(constituentId,2)%nlev)

        if (allocated(success)) deallocate(success)
        allocate(success(climatFields(constituentId,2)%nlev))
        success(:) = .true.

        if (climatFields(constituentId,2)%vertCoordKind == 2) then
          pressrefin(:) = pressrefin(:)*100. ! Conversion from hPa to Pa.
        else if (climatFields(constituentId,2)%vertCoordKind == 0) then
          where (pressrefin < modelHeightLevs(numModelLevs))
            pressrefin = modelHeightLevs(numModelLevs)
          end where
          pressrefin(:) = phf_convertZtoPressure(pressrefin, &
              modelHeightLevs,modelPressLevs, &
              climatFields(constituentId,2)%nlev,numModelLevs, &
              obsLat*MPC_RADIANS_PER_DEGREE_R8,success)
        else if (climatFields(constituentId,2)%vertCoordKind == 4) then
          pressrefin(:) = pressrefin(:) + modelHeightLevs(numModelLevs)
          pressrefin(:) = phf_convertZtoPressure(pressrefin, &
              modelHeightLevs,modelPressLevs, &
              climatFields(constituentId,2)%nlev,numModelLevs, &
              obsLat*MPC_RADIANS_PER_DEGREE_R8,success)
        else if (climatFields(constituentId,2)%vertCoordKind == 1) then
          pressrefin(:) = pressrefin(:)*modelPressLevs(numModelLevs) ! Convert from sigma to Pa
        else
          call utl_abort('clm_setColumn: Cannot handle vertical ' // &
              'coordinate of kind ' // trim(utl_str(climatFields(constituentId,2)%vertCoordKind)))
        end if

        ! Interpolate to obs lat/long (or lat) and model levels

        ! Latitudes provided in degrees to ensure backgward compatibility of results
        ! when the ozone Fortuin-Kelder climatology is used.

        call ppo_getColumn(climatFields(constituentId,2)%field, &
            climatFields(constituentId,2)%nlon,climatFields(constituentId,2)%nlat, &
            climatFields(constituentId,2)%nlev,climatFields(constituentId,2)%lon*MPC_RADIANS_PER_DEGREE_R8, &
            climatFields(constituentId,2)%lat,pressrefin, &
            obsLong*MPC_RADIANS_PER_DEGREE_R8,obsLat,numModelLevs,modelPressLevs, &
            nearestNeighbourInterp(constituentId),refprof2)

        ! Combine with upper level profile

        do level = numModelLevs, 3, -1
          if (modelPressLevs(level) < tropo_press) exit
          refprof(level) = refprof2(level)
        end do
        start = level

        ! Apply linear combination of four levels just above the tropopause

        do level = start, max(2,start-3), -1
          dt = (start+1.0-level)/5.0
          refprof(level) = dt*refprof2(level) + (1.0-dt)*refprof(level)
        end do

      end if

      if (allocated(pressrefin)) deallocate(pressrefin)
      if (allocated(success)) deallocate(success)

    end if

    ! Profile for output
    if (present(climatProfile_opt)) climatProfile_opt(:) = refprof(:)

    ! ------- Save in climatProfileSet_opt ---------

    if (.not.present(climatProfileSet_opt)) return

    if (.not.associated(climatProfileSet_opt%data1d)) then
      call oss_obsdata_alloc(climatProfileSet_opt, maxsize, dim1=numModelLevs)
      climatProfileSet_opt%nrep = 0
    end if

    ! Here, nrep will count the number of filled elements in the data arrays
    climatProfileSet_opt%nrep = climatProfileSet_opt%nrep+1

    if (climatProfileSet_opt%nrep > maxsize) then
      call utl_abort('clm_setColumn: Reach max size of array ' // &
          trim(utl_str(maxsize)) )
    end if

    ! obsIndex serves as the unique locator code
    write(climatProfileSet_opt%code(climatProfileSet_opt%nrep),'(I22)') obsIndex

    ! Save profile in climatProfileSet_opt

    climatProfileSet_opt%data1d(:,climatProfileSet_opt%nrep) = refprof(:)

  end subroutine clm_setColumn

  !--------------------------------------------------------------------------
  ! clm_getColumn
  !--------------------------------------------------------------------------
  function clm_getColumn(climatProfileSet,code) result(profile)
    !
    !:Purpose: To extract and provide column profile from climatProfileSet
    !          according to code value.
    !
    implicit none

    ! Arguments:
    type(struct_oss_obsdata), intent(inout) :: climatProfileSet  ! Profile set
    character(len=*),         intent(in)    :: code              ! unique obs identifying code
    ! Result:
    real(8) :: profile(climatProfileSet%dim1) ! retrieved array from obsdata%data1d of dimension obsdata%dim1

    ! Locals:
    integer :: status ! search success (0 = found; 1 = no data; 2 = not found)

    profile = oss_obsdata_get_array1d(climatProfileSet,code,status)
    if (status > 0) then
      call utl_abort("clm_getColumn: Code not found - " // trim(code))
    end if

  end function clm_getColumn

  !--------------------------------------------------------------------------
  ! clm_isGHG
  !--------------------------------------------------------------------------
  function clm_isGHG(varName) result(ghg_identified)

    !:Purpose: Identify if this is a greenhouse gas in the 'climScaling' file
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: varName   ! Variable name
    ! Result:
    logical :: ghg_identified                 ! Indicates if varName is a GHG

    if (trim(varName) == 'TCO2' .or. &
        ! trim(varName) == 'TCO' .or. &
        trim(varName) == 'CH4L' .or. trim(varName) == 'TCH4' .or. &
        trim(varName) == 'CF1L' .or. trim(varName) == 'CF2L' .or. &
        trim(varName) == 'N2OL' .or. trim(varName) == 'TN2O') then
      ghg_identified = .true.
    else
      ghg_identified = .false.
    end if

  end function clm_isGHG

  !--------------------------------------------------------------------------
  ! clm_getUnitsScaling
  !--------------------------------------------------------------------------
  function clm_getUnitsScaling(varName,iyear) result(scalefactor)

    !:Purpose: Identify unit scaling factor from 'climScaling' file.
    !
    !:Comments:
    !
    !    Scaling factors in 'climScaling' are dry air molar fractions and
    !    therefore the same as volumetric fractions (numbers of constituant
    !    molecules per number of dry air molecules such as ppmv or ppbv)
    !
    implicit none

    ! Arguments:
    character(len=*), intent(in) :: varName   ! Variable name
    integer, intent(in) :: iyear  ! Year for required scaling factor
    ! Result:
    real(8) :: scaleFactor        ! Units scaling factor

    ! Locals:
    integer, parameter :: numSpecies=5
    integer :: MRYear, countSpecies, varNum
    character(len=6) :: lineOffset,speciesNames(numSpecies),speciesUnits(numSpecies)
    character(len=13), parameter :: climScaling  = 'climatScaling'
    real(8) :: speciesMR(numSpecies)
    integer, external :: fnom, fclos
    integer :: ierr, nulunScaling, iosScaling
    logical :: fileExists
    character (len=128) :: ligne

    inquire(file=climScaling,exist=fileExists)
    if ( .not. fileExists ) then
      call utl_abort('clm_getUnitsScaling: Did not find file climatScaling' )
    end if

    nulunScaling = 0
    ierr = 0
    ierr = fnom(nulunScaling,trim(climScaling),'SEQ',0)
    open(unit=nulunScaling, file=trim(climScaling), status='OLD')
    iosScaling = 0
    read(nulunScaling,'(A)',iostat=iosScaling,err=20,end=20) ligne
    do while (trim(adjustl(ligne(1:7))) /= '#     c')
      read(nulunScaling,'(A)',iostat=iosScaling,err=20,end=20) ligne
    end do
    read(ligne,*) lineOffset,speciesNames(1:numSpecies)
    read(nulunScaling,'(A)',iostat=iosScaling,err=20,end=20) ligne
    read(ligne,*) lineOffset,speciesUnits(1:numSpecies)

    ! Use last year if there is no matching year
    read(nulunScaling,*) MRYear,speciesMR(1:numSpecies)
    do while (MRYear /= iyear)
      read(nulunScaling,*) MRYear,speciesMR(1:numSpecies)
    end do

    close(unit=nulunScaling)
    ierr = fclos(nulunScaling)

    ! Get unit scaling factors to get micrograms/kg

    countSpecies = 0
    scaleFactor = 1.0d0

    varNum = vnl_varnumFromVarName(varName)
    select case (varNum)
    case (BUFR_NECH_CH4)
      do countSpecies = 1, numSpecies
        if (trim(speciesNames(countSpecies)) == 'ch4') then
          scaleFactor=speciesMR(countSpecies)
          exit
        end if
      end do
    case (BUFR_NECH_CO2)
      do countSpecies = 1, numSpecies
        if (trim(speciesNames(countSpecies)) == 'co2') then
          scaleFactor=speciesMR(countSpecies)
          exit
        end if
      end do
    case (BUFR_NECH_CO)
      do countSpecies = 1, numSpecies
        if (trim(speciesNames(countSpecies)) == 'co') then
          scaleFactor=speciesMR(countSpecies)
          exit
        end if
      end do
    case (BUFR_NECH_N2O)
      do countSpecies = 1, numSpecies
        if (trim(speciesNames(countSpecies)) == 'n2o') then
          scaleFactor = speciesMR(countSpecies)
          exit
        end if
      end do
    case default
      call utl_abort('clm_getUnitsScaling: GHG scaling factor not found for ' &
          // trim(varName) )
    end select

    if (countSpecies > 0 .and. countSpecies <= numSpecies) then
      ! conversion factor to go to micrograms/kg
      if (trim(speciesUnits(countSpecies)) == 'ppm') then
        scaleFactor=scaleFactor * 1.0d03
      else if (trim(speciesUnits(countSpecies)) == 'ppb') then
        scaleFactor=scaleFactor * 1.0d06
      end if
    end if

20  if (iosScaling > 0) then
      call utl_abort('clm_getUnitsScaling: READING PROBLEM.' // &
          ' File read error message number: ' // trim(utl_str(iosScaling)))
    end if

  end function clm_getUnitsScaling

end module climatologies_mod
