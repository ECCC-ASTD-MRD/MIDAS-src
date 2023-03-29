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

module climatologies_mod
  ! MODULE climatologies_mod (prefix='clm' category='5. Observation operators')
  !
  ! :Purpose: Access to climatologies
  !
  use mathPhysConstants_mod
  use physicsFunctions_mod
  use obsSpaceData_mod
  use obsSubSpaceData_mod
  use presProfileOperators_mod
  use utilities_mod
  use varNameList_mod
  use timeCoord_mod

  implicit none
  save
  private

  ! public procedures  
  public :: clm_getOzoneProfile, clm_readOzoneClimatology
  public :: clm_readFields, clm_addToProfileSet, clm_getProfile

  ! public structures
  public :: struct_clm_field
  
  ! Number of latitudes and vertical levels in ozoneclim98 climatology file
  INTEGER, PARAMETER    :: NLATO3=19, NLEVO3=28

  ! Climatological ozone field (ppmv) and total ozone (ozoneclim98)
  REAL                  :: FOZO_r4(NLATO3,NLEVO3), TOTOZO_r4(NLATO3,12)

  ! Pressure height of ozoneclim98 climatology file vertical levels (mb)
  REAL(8)               :: PO3(NLEVO3)

  DATA PO3 /    0.010D0, 0.015D0, 0.022D0, 0.032D0, 0.046D0, 0.068D0, 0.100D0,   &
       0.150D0, 0.200D0, 0.300D0, 0.500D0, 1.000D0, 2.000D0, 3.000D0,   &
       5.000D0, 7.000D0, 10.00D0, 20.00D0, 30.00D0, 50.00D0, 70.00D0,   &
       100.0D0, 150.0D0, 200.0D0, 300.0D0, 500.0D0, 700.0D0, 1000.D0 / 

  ! Arrays containing input reference fields and fields interpolated 
  ! to obs locations
  
  type :: struct_clm_field

    !  Structure for storing reference (climatological) fields)
    !     
    !  Variable               Description
    !  --------               -----------
    !  field                  Gridded 3D field (lon,lat,vlev) or 2D field (1,lat,vlev)
    !  nlat                   number of latitudes
    !  nlon                   number of longitudes
    !  nlev                   number of vertical levels
    !  lat,lon                lat,lon grid in radians
    !  vlev                   vertical levels
    !  ivkind                 Index of vertical coordinate type. Defintion may vary according to source.
    !                         For fields read from RPN files and use of convip:
    !                             0: P is in height [m] (metres) with respect to sea level 
    !                             1: P is in stddev [sg] (0.0 -> 1.0) 
    !                             2: P is in pressure [mb] (millibars) 
    !                             3: P is in an arbitrary code 
    !                             4: P is in height [M] (metres) with respect to ground level 
    !                             5: P is in hybrid coordinates [hy] 
    !                             6: P is in theta [th] 
    !                         For use with obs                      
    
    real(8), allocatable :: field(:,:,:),lat(:),lon(:),vlev(:)
    integer :: nlev,nlon,nlat,ivkind
  
  end type struct_clm_field

contains

  !--------------------------------------------------------------------------
  ! clm_getOzoneProfile
  !--------------------------------------------------------------------------
  subroutine clm_getOzoneProfile(o3p,zlat,plev,nlev,nprf)
    !
    !:Purpose: Get ozone profile from climatology interpolated to desired P levels
    !
    IMPLICIT NONE       

    integer ,intent(in) :: nlev            ! NUMBER OF VERTICAL LEVELS
    integer ,intent(in) :: nprf            ! NUMBER OF PROFILES
    REAL(8),intent(in)  :: ZLAT(NPRF)      ! ARRAY OF LATITUDE (-90S TO 90N)
    REAL(8),intent(in)  :: PLEV(NLEV,NPRF) ! PRESSURE LEVELS (HPA)
    REAL(8),intent(out) :: O3P(NLEV,NPRF)  ! OZONE PROFILES (PPMV)

    INTEGER   :: JN, K, NUMLAT
    REAL(8)   :: QO3B(NLEVO3,NPRF)
    REAL(8)   :: PRO3(NLEVO3,NPRF)


    !* assign default qgas values if need be

    DO JN = 1, NPRF
       NUMLAT = NINT( (ZLAT(JN)+90.D0) / (180.D0/(REAL(NLATO3-1,8))) ) + 1
       DO K = 1, NLEVO3
          QO3B(K,JN) = FOZO_r4(NUMLAT,K)
       END DO
    END DO

    !* interpolation of field QO3B at NLEVO3 levels of height PO3mbb
    !* into field O3P at NLEV levels of height PLEV

    FORALL(K=1:NLEVO3) PRO3(K,:) = PO3(K)

    CALL ppo_LINTV(pro3,qo3b,nlevo3,nprf,nlev,plev,O3P)

  end subroutine clm_getOzoneProfile

  !--------------------------------------------------------------------------
  ! clm_readOzoneClimatology
  !--------------------------------------------------------------------------
  subroutine clm_readOzoneClimatology(datestamp,nlat_opt,nlev_opt,press_opt,ozone_opt)
    !
    !:Purpose: READ OZONE CLIMATOLOGICAL FIELDS
    !
    IMPLICIT NONE
    
    !Arguments
    integer            :: datestamp            ! Datestamp
    integer, intent(out), optional :: nlat_opt ! Number of latitudes
    integer, intent(out), optional :: nlev_opt ! Number of vertical levels
    real(8), allocatable, intent(out), optional :: ozone_opt(:,:) ! Ozone field
    real(8), allocatable, intent(out), optional :: press_opt(:)   ! Pressure levels

    !Locals
    INTEGER            :: IJOUR,ITIME,IMONTH,IJ,IER
    CHARACTER(len=100) :: CFILE
    INTEGER            :: NIOZO,NJOZO,NKOZO
    INTEGER, EXTERNAL  :: FNOM,FSTOUV,FSTLIR,FSTFRM,FCLOS,NEWDATE

    integer            :: IOZTEST
    integer            :: iv1,iv2,iv3,iv4,iv5,iv6

    ier = newdate(datestamp,ijour,itime,-3)

    IJ= IJOUR/100
    IMONTH = IJ - (IJ/100)*100

    ioztest=0

    CFILE='ozoneclim98'
    IV1=FNOM(IOZTEST,CFILE,'RND+R/O',0)
    IV2=FSTOUV(IOZTEST,'RND')
    IV3=FSTLIR(FOZO_r4,IOZTEST,NIOZO,NJOZO,NKOZO,-1,' ',-1,-1,IMONTH,' ','O3')
    IV4=FSTLIR(TOTOZO_r4,IOZTEST,NIOZO,NJOZO,NKOZO,-1,' ',-1,-1,-1,' ','TO')
    IV5=FSTFRM(IOZTEST)
    IV6=FCLOS(IOZTEST)

    if(iv1.lt.0.or.iv2.lt.0.or.iv3.lt.0.or.iv4.lt.0.or.iv5.lt.0.or.iv6.lt.0) then
       write(*,*) 'LES IV DE OZO_READ_CLIMATOLOGY ',iv1,iv2,iv3,iv4,iv5,iv6
       write(*,*) 'THESE NUMBERS SHOULD NOT BE NEGATIVE'
       write(*,*) 'datestamp,ijour,itime,imonth = ',datestamp,ijour,itime,imonth
       call utl_abort('Problem with file in ozo_read_climatology (ozoneclim_mod)')
    endif

    if (present(nlat_opt)) then
       nlat_opt=NLATO3
       nlev_opt=NLEVO3
       allocate(press_opt(nlevo3),ozone_opt(nlato3,nlevo3))
       press_opt(1:nlevo3)=PO3(1:nlevo3)
       ozone_opt(1:nlato3,1:nlevo3)=FOZO_r4(1:nlato3,1:nlevo3)
    endif
    
  end subroutine clm_readOzoneClimatology

  !--------------------------------------------------------------------------
  ! clm_readFields
  !--------------------------------------------------------------------------
  subroutine clm_readFields(climatFields,filename,variable,    &
                             maxNumFields,maxNumTypes,         &
                             fieldTypeRequired,success,filetype_opt)
    !
    !:Purpose:  To read climatrology (reference) fields as directed by input
    !
    ! Comments:
    !      - Fields are provided in RPN/fst files 
    !      - Reference fields can be in a separate RPN file with name provided
    !        within 'filename' if filetype='TXT' or provided as 'filename' if it
    !        refers to an RPN standaard file.
    !      - Fields assumed to be of the same units as those of the
    !        corresponding input trial fields
    !
    !:Arguments:
    !
    !    :fieldTypeRequired:
    !                   Reference profile type required. Accepted values:
    !                   ======================================================
    !                   'Climat' Climatology (or other reference) field 
    !                   'Diff'   Climatology (or other reference) for possible
    !                            use in generating differences from trial
    !                   ======================================================
    !
    !  
    implicit none

    ! Arguments:    
    type(struct_clm_field),     intent(out) :: climatFields(0:maxNumFields,maxNumTypes)
    character(len=*),           intent(in)  :: filename
    integer,                    intent(in)  :: maxNumFields,maxNumTypes
    logical,                    intent(out) :: success
    character(len=*),           intent(in)  :: variable
    character(len=*),           intent(in)  :: fieldTypeRequired(0:maxNumFields) ! Indicates
    character(len=*), optional, intent(in)  :: filetype_opt

    ! Locals:
    character(len=3) :: filetype
    character(len=256) :: fname
    character(len=4) :: varName
    character(len=12) :: etiket
    integer :: varIndex,id,nd,j,numvar,ijour,imonth,iday,itime,latIndex
    real(8) :: day
    integer :: datestamp
    integer, external :: newdate   
    integer, external :: fnom, fclos
    integer :: ierr, nulun, ios
    logical :: fileExists
    logical :: timeInterp
    integer :: ni, nj, nkeys, kind
    real(8), allocatable :: array1(:,:,:),array2(:,:,:),lvls(:),xlat(:),xlong(:) 
    real(8), allocatable :: pressclim(:),ozoneclim(:,:)
    character (len=128) :: ligne

    ! Initialize dimensions to zero
    
    climatFields(:,:)%nlon=0
    climatFields(:,:)%nlat=0
    climatFields(:,:)%nlev=1
 
    if ( trim(variable) == 'CH' ) then
      if ( all(fieldTypeRequired(:) /= 'Diff') .and. &
           all(fieldTypeRequired(:) /= 'Climat') ) then
        ! Not needed
        success=.true.
	return
      end if
    end if
   
    inquire(file=trim(filename),exist=fileExists)
    if ( .not.fileExists ) then
      write(*,*)  '----------------------------------------------------'
      write(*,*)  'clm_readFields: COULD NOT FIND file ' // trim(filename)
      write(*,*)  '----------------------------------------------------'
      success = .false.
      return
    else
      success = .true.
    end if

    ! Check for file names containing climatological fields or input directives

    if ( present(filetype_opt) ) then
      filetype = trim(fileType_opt) 
    else
      filetype = 'RPN'      
    end if
    
    nulun=0
    ierr=0
    if ( filetype == 'TXT' ) then
      ierr=fnom(nulun,trim(filename),'SEQ',0)
      if ( ierr == 0 ) then
        open(unit=nulun, file=trim(filename), status='OLD')
        ios=0
        
        if ( trim(variable) == 'CH' ) then 
        
          ! CH variable kind (for constituent fields)
          
          read(nulun,'(A)',iostat=ios,err=10,end=10) ligne
          do while (trim(adjustl(ligne(1:14))) /= 'SECTION IV:') 
            read(nulun,'(A)',iostat=ios,err=10,end=11) ligne
          end do    
          
          ! Read number of constituents with associated input file(s)
   
          read(nulun,*,iostat=ios,err=10,end=10) numvar
          if (numvar <= 0) go to 10
        else
          numvar=1
          nd=1
          call utl_abort('clm_readFields: Variable kind or name ' // &
	                 trim(variable) // ' not taken into account')  
        end if
      end if
    else if ( filetype == 'RPN' ) then
      numvar=1
      nd=1
    else if ( filetype /= 'RPN' ) then
      call utl_abort('clm_readFields: File type ' // trim(filetype) // &
                     ' not recognized') 
    end if

    if ( ierr /= 0 ) then
      call utl_abort('clm_readFields: COULD NOT OPEN file ' // trim(filename))
    end if
       
    ! Initialization

    timeInterp = .true.
    datestamp=tim_getDateStamp()
    ierr = newdate(datestamp,ijour,itime,-3)
    if ( ierr < 0 ) then
      call utl_abort('clm_readFields: Invalid datestamp ' // &
                     trim(utl_str(datestamp)) )
    end if
    imonth = MOD(ijour/100,100)
    iday = MOD(ijour,100)
    day=iday+itime*1.0D-8
    if (day > 15.) then
      day=day-15.0
    else
      day=day+15.0
    end if
    
    ! Get needed fields for each file/varIndex

    do varIndex=1,numvar

      if ( trim(variable) == 'CH' ) then 
       
        ! Read id,nd
        ! id: constituent identifier code; (0 for ozone, ...)
        ! nd: number of sets; 1 or 2 (nd=2 required when different profile 
        !       sets need to be merged according to the tropopause height 
	!       when the first set referring to strato files and teh second 
	!       to tropo fields)
       
        read(nulun,*,iostat=ios,err=10,end=10)
        read(nulun,*,iostat=ios,err=10,end=10) id,nd   
        varName=vnl_varnameFromVarnum(0,id)

        read(nulun,*,iostat=ios,err=10,end=10) fname
        inquire(file=trim(fname),exist=fileExists)
        if ( .not. fileExists ) then
          call utl_abort('clm_readFields: Did not find file ' // trim(fname))
        end if
      else
        id=varIndex
        ! Currently assumes nunmar = 1 and fname = filename. Could be extended
        fname = filename
        varname = trim(variable)
      end if
              
      do j=1,nd
       
        if ( trim(fname) ==  'ozoneclim98' ) then
          timeInterp = .false.
          call clm_readOzoneClimatology(datestamp,nlat_opt=nj,nlev_opt=nkeys, &
	                            press_opt=pressclim,ozone_opt=ozoneclim) 
          id=0
          ni=1
          allocate(array1(1,nj,nkeys),lvls(nkeys),xlat(nj),xlong(1))
          ! Convert from ppmv to microgram/kg
          array1(1,1:nj,1:nkeys) =  ozoneclim(1:nj,1:nkeys) * &
	         MPC_MOLAR_MASS_O3_R8 / (1.0d-3 * MPC_MOLAR_MASS_DRY_AIR_R8)
          lvls(1:nkeys) = pressclim(1:nkeys)
	  deallocate(ozoneclim,pressclim)
          kind = 2
          xlong(1)=0.0d0
          do latIndex = 1, nj
            xlat(latIndex) = (latIndex-1)*180.0d0/(nj-1) - 90.0d0
          end do
          etiket = '            '       
        else   
          if ( nd == 2 ) then
            read(nulun,*,iostat=ios,err=10,end=10) etiket    
          else
            etiket = '            '       
          end if                   
          call utl_readFstField(trim(fname),varName,-1,imonth,-1,etiket, &
	       ni,nj,nkeys,array1,xlat_opt=xlat,xlong_opt=xlong,         &
               lvls_opt=lvls,kind_opt=kind)
        end if      

        climatFields(id,j)%nlon=ni
        climatFields(id,j)%nlat=nj
        climatFields(id,j)%nlev=nkeys
        climatFields(id,j)%ivkind=kind   
                         
        allocate(climatFields(id,j)%field(ni,nj,nkeys))
        allocate(climatFields(id,j)%vlev(nkeys),climatFields(id,j)%lon(ni))
        allocate(climatFields(id,j)%lat(nj))
              
        climatFields(id,j)%lat(1:nj)=xlat(1:nj)*MPC_RADIANS_PER_DEGREE_R8
        climatFields(id,j)%lon(1:ni)=xlong(1:ni)*MPC_RADIANS_PER_DEGREE_R8
        where (climatFields(id,j)%lon(1:ni) < 0.0) 
	  climatFields(id,j)%lon(1:ni)=2.0*MPC_PI_R8 + climatFields(id,j)%lon(1:ni)
	end where
        climatFields(id,j)%vlev(1:nkeys)=lvls(1:nkeys)              

        if (.not.timeInterp) then

          climatFields(id,j)%field(:,:,:) = array1(:,:,:)

        else

          ! Following for interpolation as a function of days from mid-months.
             
          if (iday > 15) then
            if (imonth == 12) then
              call utl_readFstField(trim(fname),varName,-1,1,-1,etiket, &
                   ni,nj,nkeys,array2,lvls_opt=lvls,kind_opt=kind)
            else
              call utl_readFstField(trim(fname),varName,-1,imonth+1,-1, &
	           etiket,ni,nj,nkeys,array2,lvls_opt=lvls,kind_opt=kind)
            end if
          
            ! Linearly interpolate in time 
	    ! (approximately - assumes 30 day months)

            climatFields(id,j)%field(:,:,:) = (array1(:,:,:)*(30.0-day)+array2(:,:,:)*day)/30.0
             
          else if (iday <= 15) then
            if (imonth == 1) then
              call utl_readFstField(trim(fname),varName,-1,12,-1,etiket, &
                   ni,nj,nkeys,array2,lvls_opt=lvls,kind_opt=kind)
            else
              call utl_readFstField(trim(fname),varName,-1,imonth-1,-1, &
                   etiket,ni,nj,nkeys,array2,lvls_opt=lvls,kind_opt=kind)
            end if

            ! Linearly interpolate in time 
            ! (approximately - assumes 30 day months)

            climatFields(id,j)%field(:,:,:) = (array2(:,:,:)* &
		                             (30.0-day)+array1(:,:,:)*day)/30.0             
          end if          
        end if
 
        if (allocated(array1)) deallocate(array1,lvls,xlat,xlong)
        if (allocated(array2)) deallocate(array2)   
                 
      end do
    end do 
     
 10 if (ios > 0) then
      call utl_abort('clm_readFields: READING PROBLEM.' // &
                     ' File read error message number: ' // trim(utl_str(ios))) 
    end if   
    close(unit=nulun)
    ierr = fclos(nulun)
    if ( (any(fieldTypeRequired(:) == 'Diff') .or. &
          any(fieldTypeRequired(:) == 'Climat') ) .and.trim(variable) == 'CH' ) then
      do j=0,maxNumFields
        if ( climatFields(j,1)%nlon == 0 .and.           &
	     ( trim(fieldTypeRequired(j)) == 'Diff' .or. &
	       trim(fieldTypeRequired(j)) == 'Climat' ) ) then
          call utl_abort('clm_readFields: READING PROBLEM. Did not' // &
	                 ' find SECTION IV required for constituent ID ' // &
			 trim(utl_str(j)))
	end if
      end do
    end if 
	
    return    

 11 close(unit=nulun)
    ierr = fclos(nulun)
    if ( trim(variable) == 'CH' ) then
      call utl_abort('clm_readFields: READING PROBLEM. Did not find ' // &
                      'SECTION IV.') 
    end if 
	     
  end subroutine clm_readFields

  !--------------------------------------------------------------------------
  ! clm_addToProfileSet
  !--------------------------------------------------------------------------
  subroutine clm_addToProfileSet(climatFields,climatProfileSet,maxNumFields,maxNumTypes, &
                                  numModelLevs,modelPressLevs,modelHeightLevs,obsLat, &
                                  obsLong,obsIndex,maxsize,varKind_opt,varNumber_opt,tt_opt,hu_opt)
    
    !:Purpose: To determine and to store a profile at obs location as part of a cumulative
    !          profile set for a specific variable
    !
    !:Input:
    !
    !    :climatFields:           Input fields from which interpolations are done
    !    :climatProfileSet:       Input profile set
    !    :maxNumFields:           Size of first dimension for climatFields
    !    :maxNumTypes:            Size of second dimension for climatFields
    !    :numModelLevs:           Number of model levels
    !    :modelPressLevs          Model pressure array (Pa)
    !    :modelHeightLevs:        Model height (m)
    !    :obsLat:                 Latitude (rad)
    !    :obsLong:                Longitude (rad)
    !    :obsIndex:               Unique measurement identifier    
    !    :varKind_opt:            variable kind (currently only relevant for 'CH')
    !    :varNumber_opt:          Constituent id
    !    :tt_opt:                 Model temperature (Kelvin)
    !    :hu_opt:                 Specific humidity 
    !    :maxsize:                Max number of obs for which climatProfileSet will be used
    !
    !:Output:
    ! 
    !    :climatProfileSet:       Updated profile set (with one profile added for (obs_long,obs_lat))
    !
    implicit none

    ! Arguments:
    type(struct_clm_field),     intent(in)    :: climatFields(0:maxNumFields,maxNumTypes)
    type(struct_oss_obsdata),   intent(inout) :: climatProfileSet
    integer,                    intent(in)    :: maxNumFields
    integer,                    intent(in)    :: maxNumTypes
    integer,                    intent(in)    :: obsIndex
    integer,                    intent(in)    :: numModelLevs
    integer,                    intent(in)    :: maxsize
    real(8),                    intent(in)    :: modelPressLevs(numModelLevs)
    real(8),                    intent(in)    :: modelHeightLevs(numModelLevs)
    real(8),                    intent(in)    :: obsLat
    real(8),                    intent(in)    :: obsLong
    integer,          optional, intent(in)    :: varNumber_opt
    real(8),          optional, intent(in)    :: tt_opt(:)
    real(8),          optional, intent(in)    :: hu_opt(:)
    character(len=*), optional, intent(in)    :: varKind_opt
    
    ! Locals
    integer :: level,start,id
    real(8) :: tropo_press, refprof(numModelLevs),refprof2(numModelLevs),dt
    real(8), allocatable :: pressrefin(:)
    logical, allocatable :: success(:)

    if ( present(varNumber_opt) ) then
      if ( varNumber_opt  < 0 ) return
      id = varNumber_opt
    else
      id = 0
    end if

    if ( present(varKind_opt) ) then
      ! Not currently used
    end if
    
    if (climatFields(id,1)%nlat == 0) return
    
    ! Set vertical levels of reference.
    ! Convert to pressure coordinate if needed.
    
    if (allocated(pressrefin)) deallocate(pressrefin)
    allocate(pressrefin(climatFields(id,1)%nlev))
    pressrefin(:) = climatFields(id,1)%vlev(1:climatFields(id,1)%nlev)

    if (allocated(success)) deallocate(success)
    allocate(success(climatFields(id,1)%nlev))
    success(:)=.true.
    
    if (climatFields(id,1)%ivkind == 2) then
      pressrefin(:)=pressrefin(:)*100. ! Conversion from hPa to Pa.
    else if (climatFields(id,1)%ivkind == 0) then
      where (pressrefin < modelHeightLevs(numModelLevs))
        pressrefin=modelHeightLevs(numModelLevs)
      end where
      pressrefin(:) = phf_convert_z_to_pressure(pressrefin,modelHeightLevs,modelPressLevs, &
                      climatFields(id,1)%nlev,numModelLevs,obsLat,success)
    else if (climatFields(id,1)%ivkind == 4) then
      pressrefin(:)=pressrefin(:) + modelHeightLevs(numModelLevs)
      pressrefin(:) = phf_convert_z_to_pressure(pressrefin,modelHeightLevs,modelPressLevs, &
                      climatFields(id,1)%nlev,numModelLevs,obsLat,success)
    else if (climatFields(id,1)%ivkind == 1) then
      pressrefin(:)=pressrefin(:)*modelPressLevs(numModelLevs) ! Convert from sigma to Pa   
    else
       call utl_abort('clm_addToProfileSet: Cannot handle vertical coordinate of kind ' // trim(utl_str(climatFields(id,1)%ivkind)))
    end if
    
    ! Interpolate to obs lat/long (or lat) location and model level

    call clm_columnHBilin(climatFields(id,1)%field,pressrefin, &
                  climatFields(id,1)%nlon,climatFields(id,1)%nlat,climatFields(id,1)%nlev, &
                  climatFields(id,1)%lon,climatFields(id,1)%lat,obsLong,obsLat, &
                  refprof,modelPressLevs,numModelLevs)

    if (climatFields(id,2)%nlat > 0 .and. climatFields(id,2)%nlon > 0 &
        .and. climatFields(id,2)%nlev > 0) then
        
      if ( .not. present(tt_opt) ) then
        call utl_abort('clm_addToProfileSet: Missing TT for determining ' // &
	               'tropopause pressure')
      end if
      if ( any(tt_opt <= 0.0d0) ) then
        call utl_abort('clm_addToProfileSet: Invalid TT for determining ' // &
	               'tropopause pressure')
      end if
        
      ! Get second reference field (for troposphere)
        
      tropo_press=-1.0
        
      if ( present(hu_opt) ) then
        if (all(hu_opt >= 0.0D0)) then
          tropo_press=phf_get_tropopause(numModelLevs,modelPressLevs, &
	              tt_opt,modelHeightLevs,hu_opt=hu_opt)
        else
          tropo_press=phf_get_tropopause(numModelLevs,modelPressLevs, &
	              tt_opt,modelHeightLevs)
        end if
      else
        tropo_press=phf_get_tropopause(numModelLevs,modelPressLevs,tt_opt,modelHeightLevs)
      end if
	
      if (tropo_press > 0) then
          
        ! Set vertical levels of reference.
        ! Convert to pressure coordinate if needed
 
        if (allocated(pressrefin)) deallocate(pressrefin)
        allocate(pressrefin(climatFields(id,2)%nlev))    
        pressrefin(:)= climatFields(id,2)%vlev(1:climatFields(id,2)%nlev)

        if (allocated(success)) deallocate(success)
        allocate(success(climatFields(id,2)%nlev))
        success(:)=.true.

        if (climatFields(id,2)%ivkind == 2) then
          pressrefin(:)=pressrefin(:)*100. ! Conversion from hPa to Pa.
        else if (climatFields(id,2)%ivkind == 0) then
          where (pressrefin < modelHeightLevs(numModelLevs)) 
	    pressrefin=modelHeightLevs(numModelLevs)
	  end where 
          pressrefin(:) = phf_convert_z_to_pressure(pressrefin, &
	                  modelHeightLevs,modelPressLevs, &
                          climatFields(id,2)%nlev,numModelLevs, &
	                  obsLat,success)
        else if (climatFields(id,2)%ivkind == 4) then
          pressrefin(:)=pressrefin(:) + modelHeightLevs(numModelLevs)
          pressrefin(:) = phf_convert_z_to_pressure(pressrefin, &
	                   modelHeightLevs,modelPressLevs, &
                           climatFields(id,2)%nlev,numModelLevs,obsLat, &
	                   success)
        else if (climatFields(id,2)%ivkind == 1) then
          pressrefin(:)=pressrefin(:)*modelPressLevs(numModelLevs) ! Convert from sigma to Pa   
        else
          call utl_abort('clm_addToProfileSet: Cannot handle vertical ' // &
	      'coordinate of kind ' // trim(utl_str(climatFields(id,2)%ivkind)))
        end if
            
        ! Interpolate to obs lat/long (or lat) and model levels
            
        call clm_columnHBilin(climatFields(id,2)%field,pressrefin, &
             climatFields(id,2)%nlon,climatFields(id,2)%nlat,climatFields(id,2)%nlev, &
             climatFields(id,2)%lon,climatFields(id,2)%lat,obsLong,obsLat, &
             refprof2,modelPressLevs,numModelLevs)
    
      end if

       ! Combine with upper level profile
       
       do level=numModelLevs,3,-1
         if (modelPressLevs(level) < tropo_press) exit
         refprof(level)=refprof2(level)            
       end do
       start=level
            
       ! Apply linear combination of four levels just above the tropopause
        
       do level=start,max(2,start-3),-1
         dt=(start+1.0-level)/5.0
         refprof(level)=dt*refprof2(level) + (1.0-dt)*refprof(level)
      end do
                    
    end if 

    if (allocated(pressrefin)) deallocate(pressrefin)
    if (allocated(success)) deallocate(success) 

    ! ------- Save in climatProfileSet ---------
       
    if (.not.associated(climatProfileSet%data1d)) then
      call oss_obsdata_alloc(climatProfileSet, maxsize, dim1=numModelLevs)
      climatProfileSet%nrep = 0
    end if

    ! Here, nrep will count the number of filled elements in the data arrays
    climatProfileSet%nrep = climatProfileSet%nrep+1 

    if (climatProfileSet%nrep > maxsize) then
      call utl_abort('clm_addToProfilesSet: Reach max size of array ' // &
	             trim(utl_str(maxsize)) )
    end if
    
    ! obsIndex serves as the unique locator code 
    write(climatProfileSet%code(climatProfileSet%nrep),'(I22)') obsIndex
    
    ! Save profile in climatProfileSet
    
    climatProfileSet%data1d(:,climatProfileSet%nrep) = refprof(:)

  end subroutine clm_addToProfileSet
  
  !--------------------------------------------------------------------------
  ! clm_getProfile
  !--------------------------------------------------------------------------
  function clm_getProfile(climatProfileSet,code) result(profile)
    !
    !:Purpose: To extract and provide profile from climatProfileSet according to 
    !          code value.     
    !  
    implicit none
  
    ! Arguments
    type(struct_oss_obsdata), intent(inout) :: climatProfileSet  ! Profile set
    character(len=*),         intent(in)    :: code              ! unique obs identifying code    
    ! Result:
    real(8) :: profile(climatProfileSet%dim1) ! retrieved array from obsdata%data1d of dimension obsdata%dim1

    ! Locals:
    integer :: status ! search success (0 = found; 1 = no data; 2 = not found)

    profile = oss_obsdata_get_array1d(climatProfileSet,code,status)
    if (status > 0) then
      call utl_abort("clm_getProfile: Code not found - " // trim(code))
    end if
    
  end function clm_getProfile

  !--------------------------------------------------------------------------
  ! clm_columnHBilin
  !--------------------------------------------------------------------------
  subroutine clm_columnHBilin(field,vlev,nlong,nlat,nlev,xlong,xlat, &
                               plong,plat,vprof,vlevout,nlevout)
    !
    ! :Purpose: Horizontal bilinear interpolation from a 3D field to a profile at (plong,plat).
    !           Assumes vertical interpolation not needed or already done.
    !
    !           This version can be used with fields that are not part of the background state,
    !           such as climatologies.
    !
    !           This version does not depend in column_data and gridstatevector modules.
    !
    implicit none

    ! Arguments:
    integer, intent(in)  :: nlong            ! number or longitudes
    integer, intent(in)  :: nlat             ! number or latitudes
    integer, intent(in)  :: nlev             ! number of vertical levels
    integer, intent(in)  :: nlevout          ! number of target vertical levels
    real(8), intent(in)  :: field(nlong,nlat,nlev) ! 3D field
    real(8), intent(in)  :: vlev(nlev)       ! vertical levels of input field (in pressure)
    real(8), intent(in)  :: xlong(nlong)     ! longitudes (radians)
    real(8), intent(in)  :: xlat(nlat)       ! latitudes (radians)
    real(8), intent(in)  :: plong            ! target longitude (radians)
    real(8), intent(in)  :: plat             ! target latitude (radian)
    real(8), intent(in)  :: vlevout(nlevout) ! target vertical levels (in pressure)
    real(8), intent(out) :: vprof(nlevout)   ! profile at (plong,plat)
    
    ! Locals:
    real(8) :: lnvlev(nlev),lnvlevout(nlevout),plong2
    integer :: ilev,lonIndex,latIndex,i,j
    real(8) :: DLDX, DLDY, DLDP, DLW1, DLW2, DLW3, DLW4

    call utl_tmg_start(30,'--StateToColumn')

    ! Find near lat/long grid points

    if ( nlong > 1 ) then
      plong2 = plong
      if (plong2 < 0.0) plong2 = 2.D0*MPC_PI_R8 + plong2
      do lonIndex = 2, nlong
        if  (xlong(lonIndex-1) < xlong(lonIndex)) then
          if (plong2 >= xlong(lonIndex-1) .and. plong2 <= xlong(lonIndex)) exit
        else 
          ! Assumes this is a transition between 360 to 0 (if it exists). Skip over.
        end if
      end do
      lonIndex = lonIndex-1
    else
      lonIndex=0
    end if

    do latIndex = 2, nlat
      if (plat <= xlat(latIndex)) exit
    end do
    latIndex = latIndex-1

    if ( lonIndex == 0 ) then
    
      ! Set lat interpolation weights

      DLDY = (plat - xlat(latIndex))/(xlat(latIndex+1)-xlat(latIndex))

      DLW1 = (1.d0-DLDY)
      DLW2 = DLDY

      ! Set vertical interpolation weights (assumes pressure vertical coordinate)

      lnvlevout(:) = log(vlevout(:))    
      lnvlev(:) = log(vlev(:))    

      ilev = 1
      do i = 1, nlevout
        do j = ilev, nlev          
          if (lnvlevout(i) < lnvlev(j)) exit ! assumes lnvlevout and lnvlev increase with index
        end do
        ilev = j-1
        if (ilev < 1) then
          ilev = 1
        else if (ilev >= nlev) then
          ilev = nlev-1
        end if

        DLDP = (lnvlev(ilev+1)-lnvlevout(i))/(lnvlev(ilev+1)-lnvlev(ilev))
          
        vprof(i) = DLDP* (DLW1 * field(lonIndex,latIndex,ilev)      &
                        + DLW2 * field(lonIndex,latIndex+1,ilev))   & 
          + (1.d0-DLDP)* (DLW1 * field(lonIndex,latIndex,ilev+1)    &
                        + DLW2 * field(lonIndex,latIndex+1,ilev+1))  
      end do
      
    else
    
      ! Set lat/long interpolation weights

      DLDX = (plong - xlong(lonIndex))/(xlong(lonIndex+1)-xlong(lonIndex))
      DLDY = (plat - xlat(latIndex))/(xlat(latIndex+1)-xlat(latIndex))

      DLW1 = (1.d0-DLDX) * (1.d0-DLDY)
      DLW2 =       DLDX  * (1.d0-DLDY)
      DLW3 = (1.d0-DLDX) *       DLDY
      DLW4 =       DLDX  *       DLDY

      ! Set vertical interpolation weights (assumes pressure vertical coordinate)

      lnvlevout(:) = log(vlevout(:))    
      lnvlev(:) = log(vlev(:))    

      ilev = 1
      do i = 1, nlevout
        do j = ilev, nlev          
          if (lnvlevout(i) < lnvlev(j)) exit ! assumes lnvlevout and lnvlev increase with index
        end do
        ilev = j-1
        if (ilev < 1) then
          ilev = 1
        else if (ilev >= nlev) then
          ilev = nlev-1
        end if

        DLDP = (lnvlev(ilev+1)-lnvlevout(i))/(lnvlev(ilev+1)-lnvlev(ilev))
          
        vprof(i) = DLDP* (DLW1 * field(lonIndex,latIndex,ilev)      &
                        + DLW2 * field(lonIndex+1,latIndex,ilev)    &
                        + DLW3 * field(lonIndex,latIndex+1,ilev)    &
                        + DLW4 * field(lonIndex+1,latIndex+1,ilev)) &
          + (1.d0-DLDP)* (DLW1 * field(lonIndex,latIndex,ilev+1)    &
                        + DLW2 * field(lonIndex+1,latIndex,ilev+1)  &
                        + DLW3 * field(lonIndex,latIndex+1,ilev+1)  &
                        + DLW4 * field(lonIndex+1,latIndex+1,ilev+1))                               
      end do
    end if

    call utl_tmg_stop(30)

  end subroutine clm_columnHBilin

end module climatologies_mod
