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

module sqliteRead_mod
  ! MODULE sqliteRead (prefix='sqlr' category='6. Observation input/output')
  !
  ! :Purpose: To read and update SQLITE observation files. Data is stored in 
  !           obsSpaceData object.
  !

use codePrecision_mod
use obsSpaceData_mod
use mpi_mod
use fSQLite
use mathPhysConstants_mod
use mpiVar_mod
use obsUtil_mod
use utilities_mod
use bufr_mod
use ramDisk_mod
use tovs_nl_mod
use codtyp_mod
use obsVariableTransforms_mod
use obsFilter_mod

implicit none   
 
save

private
public :: sqlr_insertSqlite, sqlr_updateSqlite, sqlr_readSqlite, sqlr_query
public :: sqlr_cleanSqlite, sqlr_writeAllSqlDiagFiles

contains
  
  subroutine sqlr_initData(obsdat, vertCoord, obsValue, obsVarno, obsFlag, vertCoordType, numberData )
    !
    ! :Purpose: Initialize data values for an observation object.
    !
    implicit none

    ! arguments
    type (struct_obs), intent(inout) :: obsdat
    real(pre_obsReal), intent(in)    :: vertCoord
    real(pre_obsReal), intent(in)    :: obsValue
    integer          , intent(in)    :: obsVarno
    integer          , intent(in)    :: obsFlag 
    integer          , intent(in)    :: vertCoordType
    integer          , intent(in)    :: numberData

    call obs_bodySet_r(obsdat, OBS_PPP, numberData, vertCoord     )
    call obs_bodySet_r(obsdat, OBS_VAR, numberData, obsValue      )
    call obs_bodySet_i(obsdat, OBS_VNM, numberData, obsVarno      )
    call obs_bodySet_i(obsdat, OBS_FLG, numberData, obsFlag       )
    call obs_bodySet_i(obsdat, OBS_VCO, numberData, vertCoordType )

  end subroutine sqlr_initData


  subroutine sqlr_initHeader( obsdat, rdbSchema, familyType, headerIndex, elev, obsSat, azimuth, geoidUndulation, &
                              earthLocRadCurv, roQcFlag, instrument, zenith, cloudCover, solarZenith, &
                              solarAzimuth, terrainType, landSea, obsIdo, obsLat, obsLon, codeType, obsDate, obsTime, &
                              obsStatus, idStation, idProf, trackCellNum, modelWindSpeed)
    !
    ! :Purpose: To initialize the header information when SQLite files are read.
    !
    implicit none

    ! arguments
    type (struct_obs), intent(inout) :: obsdat
    character(len=*) , intent(in)    :: rdbSchema
    character(len=*) , intent(in)    :: idStation
    character(len=2) , intent(in)    :: familyType
    integer          , intent(in)    :: headerIndex
    integer          , intent(in)    :: obsIdo
    integer          , intent(in)    :: codeType
    integer          , intent(in)    :: obsDate
    integer          , intent(in)    :: obsTime
    integer          , intent(in)    :: obsStatus
    integer          , intent(in)    :: terrainType
    integer          , intent(in)    :: landSea
    integer          , intent(in)    :: roQcFlag
    integer          , intent(in)    :: obsSat
    integer          , intent(in)    :: instrument
    integer          , intent(in)    :: idProf
    integer          , intent(in)    :: trackCellNum
    real(pre_obsReal), intent(in)    :: geoidUndulation
    real(pre_obsReal), intent(in)    :: earthLocRadCurv
    real(pre_obsReal), intent(in)    :: elev
    real(pre_obsReal), intent(in)    :: obsLat
    real(pre_obsReal), intent(in)    :: obsLon
    real(pre_obsReal), intent(in)    :: solarAzimuth
    real(pre_obsReal), intent(in)    :: cloudCover
    real(pre_obsReal), intent(in)    :: solarZenith
    real(pre_obsReal), intent(in)    :: zenith
    real(pre_obsReal), intent(in)    :: azimuth
    real(pre_obsReal), intent(in)    :: modelWindSpeed

    call obs_setFamily( obsdat, trim(familyType), headerIndex       )
    call obs_headSet_i( obsdat, OBS_IDO, headerIndex, obsIdo        )       
    call obs_headSet_i( obsdat, OBS_ONM, headerIndex, headerIndex   )
    call obs_headSet_i( obsdat, OBS_ITY, headerIndex, codeType      )
    call obs_headSet_r( obsdat, OBS_LAT, headerIndex, obsLat        )
    call obs_headSet_r( obsdat, OBS_LON, headerIndex, obsLon        )
    call obs_headSet_i( obsdat, OBS_DAT, headerIndex, obsDate       )
    call obs_headSet_i( obsdat, OBS_ETM, headerIndex, obsTime       )
    call obs_headSet_i( obsdat, OBS_ST1, headerIndex, obsStatus     )
    call     obs_set_c( obsdat, 'STID' , headerIndex, trim(idStation) )

    if ( trim(familyType) == 'TO' ) then

      call obs_headSet_i( obsdat, OBS_SAT , headerIndex, obsSat      )
      call obs_headSet_i( obsdat, OBS_TTYP, headerIndex, terrainType     )
      call obs_headSet_i( obsdat, OBS_STYP, headerIndex, landSea     )
      call obs_headSet_i( obsdat, OBS_INS , headerIndex, instrument  )
      call obs_headSet_r( obsdat, OBS_SZA , headerIndex, zenith      )
      call obs_headSet_r( obsdat, OBS_SUN , headerIndex, solarZenith )

      if ( trim(rdbSchema) /= 'csr' ) then
        call obs_headSet_r( obsdat, OBS_AZA , headerIndex, azimuth )
      end if

      if ( trim(rdbSchema) == 'airs' .or. trim(rdbSchema) == 'iasi' .or. trim(rdbSchema) == 'cris' ) then
        call obs_headSet_r( obsdat, OBS_CLF , headerIndex, cloudCover   )
        call obs_headSet_r( obsdat, OBS_SAZ , headerIndex, solarAzimuth )  
      else if ( trim(rdbSchema) == 'amsua' .or. trim(rdbSchema) == 'amsub' .or. trim(rdbSchema) == 'atms' ) then   
        call obs_headSet_r( obsdat, OBS_SAZ , headerIndex, solarAzimuth )
      end if

    else

      call obs_headSet_r( obsdat, OBS_ALT ,headerIndex , elev )
      if ( trim(rdbSchema) == 'ro' ) then
        call obs_headSet_i( obsdat, OBS_ROQF, headerIndex, roQcFlag        )
        call obs_headSet_r( obsdat, OBS_GEOI, headerIndex, geoidUndulation )
        call obs_headSet_r( obsdat, OBS_TRAD, headerIndex, earthLocRadCurv )
        call obs_headSet_i( obsdat, OBS_SAT , headerIndex, obsSat          )
        call obs_headSet_r( obsdat, OBS_AZA , headerIndex, azimuth         )
      else if ( trim(rdbSchema) == 'al' ) then
        call obs_headSet_i( obsdat, OBS_PRFL, headerIndex, idProf          )
      else if ( trim(rdbSchema) == 'gl_ascat' ) then
         call obs_headSet_i( obsdat, OBS_FOV, headerIndex, trackCellNum   )
         call obs_headSet_r( obsdat, OBS_MWS, headerIndex, modelWindSpeed    )
      end if

    end if

  end subroutine sqlr_initHeader


  subroutine sqlr_readSqlite(obsdat, familyType, fileName )
    !
    ! :Purpose: To read SQLite namelist and files.
    ! 
    implicit none

    ! Arguments:
    type (struct_obs), intent(inout) :: obsdat     ! ObsSpaceData Structure
    character(len=*) , intent(in)    :: familyType ! Family Type may be TOVS or CONV
    character(len=*) , intent(in)    :: fileName   ! SQLite filename

    ! Locals:
    character(len=9)         :: rdbSchema
    character(len=12)        :: idStation
    integer                  :: obsIdo, obsIdd, codeType, obsDate, obsTime, obsStatus, obsFlag, obsVarno
    real(pre_obsReal)        :: elevReal, xlat, xlon, vertCoord
    real                     :: elev    , obsLat, obsLon, elevFact
    integer                  :: vertCoordType, vertCoordFact, fnom, fclos, nulnam, ierr, idProf
    real                     :: zenithReal, solarZenithReal, CloudCoverReal, solarAzimuthReal
    integer                  :: roQcFlag
    real(pre_obsReal)        :: geoidUndulation, earthLocRadCurv, obsValue, surfEmiss, biasCorrection
    real(8)                  :: geoidUndulation_R8, earthLocRadCurv_R8, azimuthReal_R8
    integer                  :: trackCellNum
    real(pre_obsReal)        :: modelWindSpeed
    real(8)                  :: modelWindSpeed_R8
    integer                  :: obsSat, landSea, terrainType, instrument, sensor, numberElem
    integer                  :: i, rowIndex, obsNlv, headerIndex, headerIndexStart, bodyIndex
    integer                  :: bitsFlagOn, bitsFlagOff, reportLocation, numBody, numHeader
    real(pre_obsReal), parameter :: zemFact = 0.01
    character(len=512)       :: query, queryData, queryHeader
    character(len=256)       :: cfgSqlite, csqlcrit, columnsHeader, columnsData
    logical                  :: finished
    real(8), allocatable     :: matdata(:,:)
    type(fSQL_DATABASE)      :: db         ! type for SQLIte  file handle
    type(fSQL_STATEMENT)     :: stmt,stmt2 ! type for precompiled SQLite statements
    type(fSQL_STATUS)        :: stat,stat2 !type for error status
    character(len=256)       :: listElem, sqlExtraDat, sqlExtraHeader, sqlNull, sqlLimit
    character(len=256),allocatable :: listElemArray(:)
    integer,allocatable            :: listElemArrayInteger(:)
    integer                  :: numberBitsOff, numberBitsOn, bitsOff(15), bitsOn(15), numberRows, numberColumns, lastId
    character(len=*), parameter :: myName = 'sqlr_readSqlite:'
    character(len=*), parameter :: myError   = myName //' ERROR: '
    namelist /NAMSQLamsua/numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLamsub/numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLairs/ numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLiasi/ numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLcris/ numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLcrisfsr/ numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLssmi/ numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLgo/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLcsr/  numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLatms/ numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLua/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLai/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLsw/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLro/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLsfc/  numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLsc/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLpr/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLal/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLgl/   numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn
    namelist /NAMSQLradar/numberElem,listElem,sqlExtraDat,sqlExtraHeader,sqlNull,sqlLimit,numberBitsOff,bitsOff,numberBitsOn,bitsOn

    write(*,*) 'Subroutine '//myName
    write(*,*) myName//': fileName   : ', trim(fileName)
    write(*,*) myName//': familyType : ', trim(familyType)
    call fSQL_open( db, trim(fileName) ,stat )
    if ( fSQL_error(stat) /= FSQL_OK ) then
      write(*,*) myError//'fSQL_open: ', fSQL_errmsg(stat)
      call utl_abort( myError//'fSQL_open' )
    end if

    query = "select schema from rdb4_schema ;"
    rdbSchema = sqlr_query(db,trim(query))
    write(*,'(4a)') myName//' rdbSchema is ---> ', trim(rdbSchema)

    sqlExtraHeader = ''
    sqlExtraDat    = ''
    sqlLimit       = ''
    sqlNull        = ''
    listElem       = ''
    numberElem = 0
    numberBitsOff = 0
    bitsOff =  0
    numberBitsOn = 0
    bitsOn = 0

    if ( trim(familyType) == 'TO' ) then
      write(listElem,*) bufr_nbt3
    else
      write(listElem,'(i5.5,",",i5.5)') bufr_nedd, bufr_neff
    end if

    vertCoordfact  = 1
    columnsData = " id_data, id_obs, vcoord, varno, obsvalue, flag "
    if ( trim(familyType) == 'TO' ) then      
      columnsHeader = " id_obs, lat, lon, codtyp, date, time, id_stn, status, id_sat, land_sea, instrument, zenith, solar_zenith "
      vertCoordType = 3
    else if ( trim(familyType) == 'GL' ) then
      columnsHeader = " id_obs, lat, lon, codtyp, date, time, id_stn"
    else if ( trim(familyType) == 'RA' ) then
      columnsHeader = " id_obs, lat, lon, codtyp, date, time, id_stn"
      vertCoordType = 1
    else
      columnsHeader = " id_obs, lat, lon, codtyp, date, time, id_stn, status, elev"  
    end if

    nulnam = 0
    ierr=fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0 )
    select case(trim(rdbSchema))
      case ( 'ua' )
        vertCoordFact = 1
        vertCoordType = 2
        read(nulnam, nml = NAMSQLua, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLua )
      case ( 'ai' )
        vertCoordFact = 1
        vertCoordType = 2
        read(nulnam, nml = NAMSQLai, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLai )
      case ( 'sw' )
        vertCoordFact = 1
        vertCoordType = 2
        read(nulnam, nml = NAMSQLsw, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLsw )
      case ( 'pr' )
        vertCoordFact = 1
        vertCoordType = 1
        read(nulnam, nml = NAMSQLpr, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
        if (mpi_myid == 0) write(*, nml =  NAMSQLpr )
      case ( 'al' )  
        columnsHeader = trim(columnsHeader)//", id_prof"
        vertCoordFact = 1
        vertCoordType = 1
        read(nulnam, nml = NAMSQLal, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
        if (mpi_myid == 0) write(*, nml =  NAMSQLal )
      case ( 'ro' )     
        columnsHeader = trim(columnsHeader)//",ro_qc_flag, geoid_undulation, earth_local_rad_curv, id_sat, azimuth"
        vertCoordFact = 1
        vertCoordType = 1
        read(nulnam, nml = NAMSQLro, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLro )
      case ( 'sf' )
        vertCoordFact = 0
        vertCoordType = 1
        read(nulnam, nml = NAMSQLsfc, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLsfc )
      case ( 'scat' )
        vertCoordFact = 0
        vertCoordType = 1
        read(nulnam, nml = NAMSQLsc, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLsc )
      case( 'airs' )
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type, cloud_cover, solar_azimuth "
        columnsData = trim(columnsData)//", surf_emiss, bias_corr "
        read(nulnam, nml = NAMSQLairs, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLairs )
      case( 'iasi' )
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type, cloud_cover, solar_azimuth "
        columnsData = trim(columnsData)//", surf_emiss, bias_corr "
        read(nulnam, nml = NAMSQLiasi, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLiasi )
      case( 'cris' )
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type, cloud_cover, solar_azimuth "
        columnsData = trim(columnsData)//", surf_emiss, bias_corr "
        read(nulnam, nml = NAMSQLcris, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLcris )
      case( 'crisfsr' )
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type, cloud_cover, solar_azimuth "
        columnsData = trim(columnsData)//", surf_emiss "
        read(nulnam, nml = NAMSQLcrisfsr, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLcrisfsr )
      case( 'amsua' )
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type, sensor, solar_azimuth "
        columnsData = trim(columnsData)//", bias_corr "
        read(nulnam, nml = NAMSQLamsua, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLamsua )
      case( 'amsub' )
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type, sensor, solar_azimuth "
        columnsData = trim(columnsData)//", bias_corr "
        read(nulnam, nml = NAMSQLamsub, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLamsub )
      case( 'atms')
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type, sensor, solar_azimuth "
        columnsData = trim(columnsData)//", bias_corr "
        read(nulnam, nml = NAMSQLatms, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLatms )
      case( 'ssmi' )
        columnsHeader = trim(columnsHeader)//", azimuth, terrain_type "
        columnsData = trim(columnsData)//", bias_corr "
        read(nulnam, nml = NAMSQLssmi, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
        if (mpi_myid == 0) write(*, nml = NAMSQLssmi )
      case( 'csr' )
        columnsData = trim(columnsData)//", bias_corr "
        read(nulnam, nml = NAMSQLcsr, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
        if (mpi_myid == 0) write(*, nml =  NAMSQLcsr )
      case( 'gl' )
        read(nulnam, nml = NAMSQLgl, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
        if (mpi_myid == 0) write(*, nml =  NAMSQLgl )
      case( 'gl_ascat' )
        columnsHeader = trim(columnsHeader)//", track_cell_no, mod_wind_spd "
        read(nulnam, nml = NAMSQLgl, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
        if (mpi_myid == 0) write(*, nml =  NAMSQLgl )
      case( 'ra' )
        read(nulnam, nml = NAMSQLradar, iostat = ierr )
        if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
        if (mpi_myid == 0) write(*, nml =  NAMSQLradar ) 
      case DEFAULT
        write(*,*) myError//'Unsupported  SCHEMA ---> ',trim(rdbSchema), ' ABORT!!! '
        call utl_abort( myError//'Unsupported  SCHEMA in SQLITE file!' )
    end select
    ierr=fclos( nulnam )

    ! Set the observation variable transforms
    if (numberElem > 0) then
      call utl_splitString(listElem,',', & ! IN
                           listElemArray)  ! OUT
      call utl_stringArrayToIntegerArray(listElemArray,      & ! IN
                                         listElemArrayInteger) ! OUT
      call ovt_setup(listElemArrayInteger) ! IN
      deallocate(listElemArrayInteger)
      deallocate(listElemArray)
    end if

    ! Compose SQL queries
    bitsFlagOff=0
    do i = 1, numberBitsOff
      bitsFlagOff = ibset ( bitsFlagOff, 13 - bitsOff(i) )
    end do
    bitsFlagOn=0
    do i = 1, numberBitsOn
      bitsFlagOn = ibset ( bitsFlagOn, 13 - bitsOn(i) )
    end do

    write(cfgSqlite, '(i6)' ) bitsFlagOff
    csqlcrit=trim(" (flag & "//cfgSqlite)//" = 0 "

    if ( numberBitsOn > 0 ) then
      write(cfgSqlite, '(i6)' ) bitsFlagOn
      csqlcrit = trim(csqlcrit)//trim(" and flag & "//cfgSqlite)
      csqlcrit = trim(csqlcrit)//" = "//cfgSqlite
    end if

    csqlcrit = trim(csqlcrit)//" or flag is null) and varno in ( "//trim(listElem)//" )"//trim(SQLNull)//trim(sqlExtraDat)
    queryData= "select "//columnsData
    queryData = trim(queryData)//trim(" from data where ")//trim(csqlcrit)//trim(sqlLimit)//";"
    queryHeader="select "//trim(columnsHeader)//" from header "//trim(sqlExtraHeader)//" order by id_obs;"
    write(*,'(4a)') myName//': ',trim(rdbSchema),' queryData    --> ', trim(queryData)
    write(*,'(4a)') myName//': ',trim(rdbSchema),' queryHeader --> ', trim(queryHeader)
    write(*,*)' ========================================== '

    if ( trim(rdbSchema)=='pr' .or. trim(rdbSchema)=='sf' ) then
      elevFact=1.
    else
      elevFact=0.
    end if

    call fSQL_prepare( db, trim(queryHeader) , stmt, stat )
    call fSQL_prepare( db, trim(queryData), stmt2, stat2 )
    if ( fSQL_error(stat)  /= FSQL_OK ) call sqlr_handleError(stat ,'fSQL_prepare hdr: ')
    if ( fSQL_error(stat2) /= FSQL_OK ) call sqlr_handleError(stat2,'fSQL_prepare dat: ')

    headerIndex  = obs_numHeader(obsdat)
    bodyIndex = obs_numBody(obsdat)
    headerIndexStart = headerIndex
    numHeader = obs_numHeader(obsdat)
    numBody   = obs_numBody(obsdat)
    write(*,*) myName//' DEBUT numheader  =', numHeader
    write(*,*) myName//' DEBUT numbody    =', numBody
    call fSQL_get_many( stmt2, nrows=numberRows, ncols=numberColumns, mode=FSQL_REAL8 )
    write(*,*) myName//'  numberRows numberColumns =', numberRows, numberColumns
    write(*,*) myName//'  rdbSchema = ', rdbSchema
    write(*,*)' ========================================== '
    allocate( matdata(numberRows, numberColumns) )
    matdata = 0.0d0
    call fSQL_fill_matrix( stmt2, matdata )
    call fSQL_free_mem( stmt2 )
    call fSQL_finalize( stmt2 )

    obsNlv = 0
    lastId = 1

    HEADER: do

      call fSQL_get_row( stmt, finished )   ! Fetch the next row
      if (finished) then                    ! exit LOOP WHEN LAST ROW HAS BEEN FETCHED
        if ( headerIndex > 1 .and. obsNlv > 0 ) &
          call sqlr_initHeader( obsdat, rdbSchema, familyType, headerIndex, elevReal, obsSat, real(azimuthReal_R8,kind=pre_obsReal), geoidUndulation, &
                                earthLocRadCurv, roQcFlag, instrument, real(zenithReal,kind=pre_obsReal), real(cloudCoverReal,kind=pre_obsReal), real(solarZenithReal,kind=pre_obsReal), &
                                real(solarAzimuthReal,kind=pre_obsReal), terrainType, landSea, obsIdo, xlat, xlon, codeType, obsDate, obsTime/100, obsStatus, idStation, idProf, &
                                trackCellNum, modelWindSpeed )
        exit HEADER
      end if

      ! The query result is inserted into variables
      terrainType = MPC_missingValue_INT
      landSea     = MPC_missingValue_INT
      sensor      = MPC_missingValue_INT
      cloudCoverReal = MPC_missingValue_R4
      elev = 0.; elevReal = 0.
      solarAzimuthReal = MPC_missingValue_R4
      solarZenithReal  = MPC_missingValue_R4
      zenithReal = MPC_missingValue_R4
      instrument = MPC_missingValue_INT
      azimuthReal_R8 = MPC_missingValue_R8

      call fSQL_get_column( stmt, COL_INDEX = 1, INT_VAR   = obsIdo    )
      call fSQL_get_column( stmt, COL_INDEX = 2, REAL_VAR  = obsLat    )
      call fSQL_get_column( stmt, COL_INDEX = 3, REAL_VAR  = obsLon    )
      call fSQL_get_column( stmt, COL_INDEX = 4, INT_VAR   = codeType  )
      call fSQL_get_column( stmt, COL_INDEX = 5, INT_VAR   = obsDate   )
      call fSQL_get_column( stmt, COL_INDEX = 6, INT_VAR   = obsTime   )
      call fSQL_get_column( stmt, COL_INDEX = 7, CHAR_VAR  = idStation )

      if ( trim(familyType) == 'TO' ) then

        call fSQL_get_column( stmt, COL_INDEX = 8,  INT_VAR   = obsStatus )
        call fSQL_get_column( stmt, COL_INDEX = 9,  INT_VAR   = obsSat )
        call fSQL_get_column( stmt, COL_INDEX = 10, INT_VAR   = landSea ,   INT_MISSING=MPC_missingValue_INT  )
        call fSQL_get_column( stmt, COL_INDEX = 11, INT_VAR   = instrument, INT_MISSING=MPC_missingValue_INT  )
        call fSQL_get_column( stmt, COL_INDEX = 12, REAL_VAR  = zenithReal, REAL_MISSING=MPC_missingValue_R4  )
        call fSQL_get_column( stmt, COL_INDEX = 13, REAL_VAR  = solarZenithReal, REAL_MISSING=MPC_missingValue_R4 )

        if ( trim(rdbSchema) /= 'csr' ) then

          call fSQL_get_column( stmt, COL_INDEX = 14, REAL8_VAR  = azimuthReal_R8 )
         
          call fSQL_get_column( stmt, COL_INDEX = 15, INT_VAR   = terrainType, INT_MISSING=MPC_missingValue_INT )

        end if

        if ( trim(rdbSchema) == 'airs' .or. trim(rdbSchema) == 'iasi' .or. trim(rdbSchema) == 'cris' ) then

          call fSQL_get_column( stmt, COL_INDEX = 16, REAL_VAR = cloudCoverReal, REAL_MISSING=MPC_missingValue_R4   )
          call fSQL_get_column( stmt, COL_INDEX = 17, REAL_VAR = solarAzimuthReal, REAL_MISSING=MPC_missingValue_R4 )

        else if ( trim(rdbSchema) == 'amsua' .or. trim(rdbSchema) == 'amsub' .or. trim(rdbSchema) == 'atms' ) then

          call fSQL_get_column( stmt, COL_INDEX = 16, INT_VAR  = sensor, INT_MISSING=MPC_missingValue_INT )
          call fSQL_get_column( stmt, COL_INDEX = 17, REAL_VAR = solarAzimuthReal, REAL_MISSING=MPC_missingValue_R4 )

        end if

        if ( instrument == 420 ) obsSat  = 784
        if ( codeType == 202 .and. instrument == 620 ) instrument = 2046
        if ( sensor == MPC_missingValue_INT ) then
          sensor = 0
          if (instrument == MPC_missingValue_INT ) instrument = 0
        else
          instrument = obsu_cvt_obs_instrum(sensor)
        end if

      else if ( trim(familyType) == 'GL' ) then

        ! It does not have the obsStatus column.

        if ( idStation(1:6) == 'METOP-') then
          call fSQL_get_column( stmt, COL_INDEX = 8,   INT_VAR  = trackCellNum )
          if (trackCellNum > 21) trackCellNum = 43 - trackCellNum
          call fSQL_get_column( stmt, COL_INDEX = 9, REAL8_VAR  = modelWindSpeed_R8 )
          modelWindSpeed = modelWindSpeed_R8
        end if

      else if ( trim(familyType) == 'RA' ) then

        ! Nothing more to read for RA now.
        ! It does not have the obsStatus column.
      else  ! familyType = CONV

        call fSQL_get_column( stmt, COL_INDEX = 8,  INT_VAR  = obsStatus )
        call fSQL_get_column( stmt, COL_INDEX = 9, REAL_VAR  = elev, REAL_MISSING=MPC_missingValue_R4 )
        elevReal=elev

        if ( trim(rdbSchema)=='ro' ) then

          call fSQL_get_column( stmt, COL_INDEX = 10, INT_VAR   = roQcFlag, INT_MISSING=MPC_missingValue_INT )
          call fSQL_get_column( stmt, COL_INDEX = 11, REAL8_VAR = geoidUndulation_R8 )
          geoidUndulation = geoidUndulation_R8
          call fSQL_get_column( stmt, COL_INDEX = 12, REAL8_VAR = earthLocRadCurv_R8 )
          earthLocRadCurv = earthLocRadCurv_R8
          call fSQL_get_column( stmt, COL_INDEX = 13, INT_VAR   = obsSat, INT_MISSING=MPC_missingValue_INT )
          call fSQL_get_column( stmt, COL_INDEX = 14, REAL8_VAR = azimuthReal_R8 )

        else if ( trim(rdbSchema)=='al' ) then
          call fSQL_get_column( stmt, COL_INDEX = 10, INT_VAR   = idProf )
        end if

      end if  ! TOVS or CONV

      if ( obsLon < 0. ) obsLon = obsLon + 360.
      xlat = obsLat * MPC_RADIANS_PER_DEGREE_R8
      xlon = obsLon * MPC_RADIANS_PER_DEGREE_R8
      headerIndex = headerIndex + 1 

      obsNlv = 0
      DATA: do rowIndex = lastId,numberRows
        ! ---------------------------------------------------
        if ( int(matdata(rowIndex,2)) > obsIdo ) then 

          call sqlr_initHeader( obsdat, rdbSchema, familyType, headerIndex, elevReal, obsSat,   &
               real(azimuthReal_R8,kind=pre_obsReal), geoidUndulation, earthLocRadCurv,         &
               roQcFlag, instrument, real(zenithReal,kind=pre_obsReal),                         &
               real(cloudCoverReal,kind=pre_obsReal), real(solarZenithReal,kind=pre_obsReal),   &
               real(solarAzimuthReal,kind=pre_obsReal), terrainType, landSea, obsIdo, xlat,     &
               xlon, codeType, obsDate, obsTime/100, obsStatus, idStation, idProf,              &
               trackCellNum, modelWindSpeed )
          exit DATA

        else if ( int(matdata(rowIndex,2)) == obsIdo ) then

          if ( headerIndex == headerIndexStart .and. obsNlv == 0 ) then
            call sqlr_initHeader( obsdat, rdbSchema, familyType, headerIndex, elevReal, obsSat,   &
                 real(azimuthReal_R8,kind=pre_obsReal), geoidUndulation, earthLocRadCurv,         &
                 roQcFlag, instrument, real(zenithReal,kind=pre_obsReal),                         &
                 real(cloudCoverReal,kind=pre_obsReal), real(solarZenithReal,kind=pre_obsReal),   &
                 real(solarAzimuthReal,kind=pre_obsReal), terrainType, landSea, obsIdo, xlat,     &
                 xlon, codeType, obsDate, obsTime/100, obsStatus, idStation, idProf,              &
                 trackCellNum, modelWindSpeed )
          end if

          lastId = rowIndex + 1
          obsIdd = int(matdata(rowIndex,1))
          vertCoord = matdata(rowIndex,3)
          obsVarno = int(matdata(rowIndex,4))
          obsValue = matdata(rowIndex,5)
          obsFlag = int(matdata(rowIndex,6))

          bodyIndex = bodyIndex + 1
          obsNlv   = obsNlv + 1
          call obs_bodySet_i(obsdat, OBS_IDD, bodyIndex, obsIdd)

          if (trim(rdbSchema) == 'airs' .or. trim(rdbSchema) == 'iasi' .or. trim(rdbSchema) == 'cris' ) then

            surfEmiss = matdata(rowIndex,7)
            call obs_bodySet_r(obsdat, OBS_SEM, bodyIndex, surfEmiss * zemFact)

            biasCorrection = matdata(rowIndex,8)

            if ( obs_columnActive_RB(obsdat,OBS_BCOR) ) then
              call obs_bodySet_r(obsdat, OBS_BCOR, bodyIndex, biasCorrection )
            end if

          end if

          if ( trim(rdbSchema) == 'amsua' .or. trim(rdbSchema) == 'amsub' .or. &
               trim(rdbSchema) == 'atms'  .or.   trim(rdbSchema) =="ssmi" .or. &
               trim(rdbSchema) =="csr" ) then

            biasCorrection = matdata(rowIndex,7)

            if ( obs_columnActive_RB(obsdat,OBS_BCOR) ) &
                 call obs_bodySet_r(obsdat, OBS_BCOR, bodyIndex, biasCorrection )

          end if

          if ( trim(familyType) == 'TO' ) then

            call sqlr_initData(obsdat, vertCoord, obsValue, obsVarno, obsFlag, vertCoordType, bodyIndex)

          else ! CONV
 
           call sqlr_initData( obsdat, vertCoord * vertCoordFact + elevReal * elevFact, &
                               obsValue, obsVarno, obsFlag, vertCoordType, bodyIndex)

           if (.not. filt_bufrCodeAssimilated(obsVarno) .and. &
               .not. ovt_bufrCodeSkipped(obsVarno)) then
             call obs_bodySet_i( obsdat, OBS_IDD, bodyIndex + 1, -1)
             call sqlr_initData( obsdat, vertCoord * vertCoordFact + elevReal * elevFact, &
                                 obs_missingValue_R, ovt_getDestinationBufrCode(obsVarno), &
                                 0, vertCoordType, bodyIndex + 1 )
             bodyIndex = bodyIndex + 1
             obsNlv = obsNlv + 1
             if (ovt_isWindObs(obsVarno)) then
               ! Add an extra row for the other wind component
               call obs_bodySet_i( obsdat, OBS_IDD, bodyIndex + 1, -1)
               call sqlr_initData( obsdat, vertCoord * vertCoordFact + elevReal * elevFact, &
                                   obs_missingValue_R, ovt_getDestinationBufrCode(obsVarno,extra_opt=.true.), &
                                   0, vertCoordType, bodyIndex + 1 )
               bodyIndex = bodyIndex + 1
               obsNlv = obsNlv + 1
             end if
           end if
           
         end if       ! TOVS or CONV

       end if          !  obsIdo   

     end do DATA  ! END OF DATA LOOP

     if ( obsNlv > 0 ) then

        if ( headerIndex == 1 ) call obs_headSet_i(obsdat, OBS_RLN, headerIndex, 1 )

        call obs_headSet_i(obsdat, OBS_NLV, headerIndex, obsNlv )

        if ( headerIndex > 1 ) then
          reportLocation = obs_headElem_i(obsdat, OBS_RLN, headerIndex - 1 ) +  obs_headElem_i(obsdat, OBS_NLV, headerIndex - 1 )
          call obs_headSet_i(obsdat, OBS_RLN, headerIndex, reportLocation )
        end if   

        if ( lastId > numberRows ) &
          call sqlr_initHeader( obsdat, rdbSchema, familyType, headerIndex, elevReal, obsSat, real(azimuthReal_R8,kind=pre_obsReal), geoidUndulation, &
                                earthLocRadCurv, roQcFlag, instrument, real(zenithReal,kind=pre_obsReal), real(cloudCoverReal,kind=pre_obsReal), real(solarZenithReal,kind=pre_obsReal), &
                                real(solarAzimuthReal,kind=pre_obsReal), terrainType, landSea, obsIdo, xlat, xlon, codeType,obsDate, obsTime/100, obsStatus, idStation, idProf, trackCellNum, modelWindSpeed )

      else

        headerIndex = headerIndex - 1

      end if

    end do HEADER ! HEADER

    deallocate(matdata)
    numHeader = obs_numHeader(obsdat)
    numBody   = obs_numBody(obsdat)
    write(*,*)  myName//' FIN numheader  =', numHeader
    write(*,*)  myName//' FIN numbody    =', numBody
    write(*,*)  myName//' fin header '
    call fSQL_finalize( stmt )
    call fSQL_close( db, stat ) 
    write(*,*) 'end subroutine: ', myName

  end subroutine sqlr_readSqlite


  function sqlr_query(db,query)
    !
    ! :Purpose: To create a query to read an SQLite file
    !
    implicit none

    ! arguments
    type(fSQL_DATABASE)  :: db    ! type handle for  SQLIte file
    character(len = *)   :: query
    ! locals
    character(len = 256) :: sqlr_query, result
    logical finished
    type(fSQL_STATEMENT) :: stmt !  prepared statement for  SQLite
    type(fSQL_STATUS)    :: stat !type error status

    result=''
    call fSQL_prepare( db, trim(query), stmt, stat)
    if ( fSQL_error(stat) /= FSQL_OK ) call sqlr_handleError(stat,'fSQL_prepare: ')
    finished=.false.
    call fSQL_get_row( stmt, finished )

    ! Put result of query into variable
    call fSQL_get_column( stmt, COL_INDEX = 1, CHAR_VAR = result )
    call fSQL_get_row( stmt, finished )
    if ( .not. finished ) write(*,*)' SQL QUERY ---> QUERY RETURNS MORE THAN ONE ROW...  '
    call fSQL_finalize( stmt )
    sqlr_query = trim(result)

  end function sqlr_query


  subroutine sqlr_handleError(stat, message)
    implicit none

    type(FSQL_STATUS)  :: stat
    character(len = *) :: message

    write(*,*) message, fSQL_errmsg(stat)
    call utl_abort( trim(message) )

  end subroutine sqlr_handleError

  
  subroutine sqlr_updateSqlite(db, obsdat, familyType, fileName, fileNumber )
    implicit none
    ! arguments
    type(fSQL_DATABASE),intent(inout):: db
    type (struct_obs),  intent(inout) :: obsdat
    character(len=*),   intent(in)    :: fileName   
    character(len=*),   intent(in)    :: familyType ! Observation Family Type 
    integer         ,   intent(in)    :: fileNumber ! FILE NUMBER ASSOCIATED WITH db
    !  locals
    type(fSQL_STATEMENT)             :: stmt ! prepared statement for  SQLite
    type(fSQL_STATUS)                :: stat ! type error status
    integer                          :: obsRln, obsNlv, obsIdf, obsIdd, obsFlag
    integer                          :: obsIdo, obsStatus, last_question
    integer                          :: itemId
    integer                          :: headerIndex, bodyIndex, numberUpdateItems
    character(len =   3)             :: item, itemUpdateList(15)
    integer                          :: updateList(20), fnom, fclos, nulnam, ierr
    character(len =   9)             :: item2
    character(len = 128)             :: query
    character(len = 356)             :: itemChar,item2Char
    logical                          :: back
    real                             :: romp, obsValue
    character(len=*), parameter      :: myName = 'sqlr_updateSqlite:'
    character(len=*), parameter      :: myError = myName //' ERROR: '
    namelist/namSQLUpdate/ numberUpdateItems, itemUpdateList

    write(*,*) myName//' Starting ===================  '

    ! set default values of namelist variables
    itemUpdateList(:) = ''
    numberUpdateItems = 0

    ! Read the namelist for directives
    nulnam = 0
    ierr   = fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0 )
    read(nulnam,nml = namSQLUpdate, iostat = ierr )
    if ( ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
    if ( mpi_myid == 0 ) write(*, nml = namSQLUpdate )
    ierr = fclos( nulnam )

    write(*,*) myName//': Family Type   = ', trim(familyType)
    write(*,*) myName//': Number of items to update: ', numberUpdateItems
    write(*,*) myName//': File Name     = ', trim(fileName)
    write(*,*) myName//': Missing Value = ', MPC_missingValue_R8    

    ! CREATE QUERY
    itemChar='  '

    do itemId = 1, numberUpdateItems
      item = itemUpdateList(itemId)
      write(*,*)'Updating ', itemId, item
      select case(item)
        case('OMA')
          updateList(itemId) = OBS_OMA
          item2='oma'
        case('OMP')
          updateList(itemId) = OBS_OMP
          item2='omp'
        case('VAR')
          updateList(itemId) = OBS_VAR
          item2='obsvalue'
        case('OER')
          updateList(itemId) = OBS_OER
          item2='obs_error'
        case('FGE')
          updateList(itemId) = OBS_HPHT
          item2='fg_error'
        case DEFAULT
          write(*,*)'invalid item: ', item2,' EXIT sqlr_updateSQL!!!'
          call utl_abort( myError//'invalid item ' )
      end select
      itemChar = trim(itemChar)//','//trim(item2)//trim(' = ? ')
    end do

    back=.true.
    last_question  = scan(itemChar, '?', back)
    item2Char   = itemChar(1:last_question)
    itemChar    = item2Char
    query = ' update data set flag = ? '//trim(itemChar)
    query = trim(query)//' where id_data = ?  ;'
    write(*,*) ' Update query --->  ', query
    call fSQL_prepare( db, query , stmt, stat )
    if ( fSQL_error(stat) /= FSQL_OK ) call sqlr_handleError(stat, 'fSQL_prepare : ')
    call fSQL_begin(db)

    HEADER: do headerIndex = 1, obs_numHeader(obsdat)
 
      obsIdf = obs_headElem_i( obsdat,OBS_IDF, headerIndex )
 
      if ( obsIdf /= fileNumber) cycle HEADER
      obsIdo = obs_headElem_i( obsdat, OBS_IDO, headerIndex )
      obsRln = obs_headElem_i( obsdat, OBS_RLN, headerIndex )
      obsNlv = obs_headElem_i( obsdat, OBS_NLV, headerIndex )

      BODY: do bodyIndex = obsRln, obsNlv + obsRln - 1

        obsFlag = obs_bodyElem_i( obsdat, OBS_FLG, bodyIndex )
        obsIdd  = obs_bodyElem_i( obsdat, OBS_IDD, bodyIndex )

        call fSQL_bind_param(stmt, PARAM_INDEX = 1,   INT_VAR  = obsFlag  )
        ITEMS: do itemId = 1, numberUpdateItems

          obsValue = obs_bodyElem_r(obsdat, OBS_VAR, bodyIndex)
          if ( obsValue /= obs_missingValue_R ) then  
            romp = obs_bodyElem_r(obsdat, updateList(itemId), bodyIndex )
            if ( romp == obs_missingValue_R ) then
              call fSQL_bind_param(stmt, PARAM_INDEX = itemId + 1                  )  ! sql null values
            else
              call fSQL_bind_param(stmt, PARAM_INDEX = itemId + 1, REAL_VAR = romp )
            end if
          end if

        end do ITEMS

        call fSQL_bind_param(stmt, PARAM_INDEX = numberUpdateItems + 2, INT_VAR  = obsIdd )
        call fSQL_exec_stmt (stmt)

      end do BODY

    end do HEADER

    call fSQL_finalize( stmt )

    if ( trim(familyType) /= 'GL'.and. trim(familyType) /= 'RA' ) then

       ! UPDATES FOR THE STATUS FLAGS IN THE HEADER TABLE
       query = ' update header set status  = ? where id_obs = ? '
       call fSQL_prepare( db, query , stmt, stat)
       if ( fSQL_error(stat) /= FSQL_OK ) call sqlr_handleError(stat,'fSQL_prepare : ')

       HEADER2: do headerIndex = 1,obs_numHeader(obsdat)

          obsIdf = obs_headElem_i(obsdat, OBS_IDF, headerIndex )
          if ( obsIdf /= fileNumber ) cycle HEADER2
          obsIdo    = obs_headElem_i(obsdat, OBS_IDO, headerIndex )
          obsStatus = obs_headElem_i(obsdat, OBS_ST1, headerIndex )
          call fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR  = obsStatus )
          call fSQL_bind_param( stmt, PARAM_INDEX = 2, INT_VAR  = obsIdo )
          call fSQL_exec_stmt ( stmt)

       end do HEADER2
    
       call fSQL_finalize( stmt )

    end if

    call fSQL_commit(db)
    write(*,*) myName//' End ===================  ', trim(familyType)

  end subroutine sqlr_updateSqlite


  subroutine sqlr_insertSqlite( db, obsdat, familyType, fileName, fileNumber )
    implicit none
    ! arguments
    type(fSQL_DATABASE)    :: db   ! type for SQLIte  file handle
    type(struct_obs)       :: obsdat
    character(len=*)       :: familyType
    character(len=*)       :: fileName
    integer                :: fileNumber
    ! locals
    integer                :: itemInsertList(15), numberInsertItems
    type(fSQL_STATEMENT)   :: stmt ! type for precompiled SQLite statements
    type(fSQL_STATUS)      :: stat !type for error status
    integer                :: obsVarno, obsFlag, vertCoordType, fnom, fclos, nulnam, ierr 
    real                   :: obsValue, OMA, OMP, OER, FGE, PPP
    integer                :: numberInsert, headerIndex, bodyIndex, numHeader
    integer                :: obsNlv, obsRln, obsIdd, obsIdo, obsIdf, insertItem
    character(len = 256)   :: query
    logical                :: llok    
    character(len=*), parameter :: myName = 'sqlr_insertSqlite:'
    character(len=*), parameter :: myError   = myName //' ERROR: '
    namelist/namSQLInsert/ numberInsertItems, itemInsertList

    write(*,*)  myName//' --- Starting ---   '
    numHeader = obs_numHeader(obsdat)
    write(*,*)' FAMILY ---> ', trim(familyType), '  headerIndex  ----> ', numHeader
    write(*,*)' fileName -> ', trim(fileName)   

    ! set default values of namelist variables
    numberInsertItems = 0
    itemInsertList(:) = 0

    nulnam = 0
    ierr=fnom(nulnam, './flnml', 'FTN+SEQ+R/O', 0 )
    read(nulnam, nml = namSQLInsert, iostat = ierr )
    if (ierr /= 0 ) call utl_abort( myError//'Error reading namelist' )
    if (mpi_myid == 0) write(*, nml = namSQLInsert )
    ierr=fclos( nulnam )

    write(*,*)' INSERT INTO SQLITE FILE ELEMENTS :--> ',( itemInsertList(insertItem), insertItem = 1, numberInsertItems )

    select case( trim( familyType ) )

      case( 'SF', 'SC', 'GP' )

        query = 'insert into data (id_obs,varno,vcoord,obsvalue,flag,oma,omp,fg_error,obs_error) values(?,?,?,?,?,?,?,?,?);'

      case DEFAULT

        query = 'insert into data (id_obs,varno,vcoord,vcoord_type,obsvalue,flag,oma,omp,fg_error,obs_error) values(?,?,?,?,?,?,?,?,?,?);'

    end select

    query=trim(query)
    write(*,*) ' === Family Type === ',trim(familyType)
    write(*,*) ' Insert query = ', trim(query)

    call fSQL_begin(db)
    call fSQL_prepare( db, query, stmt, stat )
    if ( fSQL_error(stat) /= FSQL_OK ) call sqlr_handleError(stat, 'fSQL_prepare : ')

    numberInsert=0
    HEADER: do headerIndex = 1, obs_numHeader(obsdat)

      obsIdf = obs_headElem_i(obsdat, OBS_IDF, headerIndex )
      if ( obsIdf /= fileNumber ) cycle HEADER
      obsIdo = obs_headElem_i(obsdat, OBS_IDO, headerIndex )
      obsRln = obs_headElem_i(obsdat, OBS_RLN, headerIndex )
      obsNlv = obs_headElem_i(obsdat, OBS_NLV, headerIndex )

      BODY: do bodyIndex = obsRln, obsNlv + obsRln -1

        obsIdd        = obs_bodyElem_i(obsdat, OBS_IDD , bodyIndex )
        obsVarno      = obs_bodyElem_i(obsdat, OBS_VNM , bodyIndex )
        obsFlag       = obs_bodyElem_i(obsdat, OBS_FLG , bodyIndex )
        vertCoordType = obs_bodyElem_i(obsdat, OBS_VCO , bodyIndex )
        obsValue      = obs_bodyElem_r(obsdat, OBS_VAR , bodyIndex )
        OMA           = obs_bodyElem_r(obsdat, OBS_OMA , bodyIndex )
        OMP           = obs_bodyElem_r(obsdat, OBS_OMP , bodyIndex )
        OER           = obs_bodyElem_r(obsdat, OBS_OER , bodyIndex )
        FGE           = obs_bodyElem_r(obsdat, OBS_HPHT, bodyIndex )
        PPP           = obs_bodyElem_r(obsdat, OBS_PPP , bodyIndex )
        llok = .false.

        do insertItem = 1, numberInsertItems
          if ( obsVarno == itemInsertList(insertItem) ) llok=.true.
        end do

        if ( obsIdd == -1 ) then
          if ( llok ) then
            select case(trim(familyType))
              case( 'SF', 'SC', 'GP' )
                call fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR  = obsIdo   )
                call fSQL_bind_param( stmt, PARAM_INDEX = 2, INT_VAR  = obsVarno )
                call fSQL_bind_param( stmt, PARAM_INDEX = 3, REAL_VAR = PPP      )
                if ( obsValue == obs_missingValue_R ) then          ! sql null values
                  call fSQL_bind_param( stmt, PARAM_INDEX = 4                      )
                else
                  call fSQL_bind_param( stmt, PARAM_INDEX = 4, REAL_VAR = obsValue )
                end if
                call fSQL_bind_param( stmt, PARAM_INDEX = 5, INT_VAR  = obsFlag  )
                if ( OMA == obs_missingValue_R ) then
                  call fSQL_bind_param( stmt, PARAM_INDEX = 6                    ) 
                else 
                  call fSQL_bind_param( stmt, PARAM_INDEX = 6, REAL_VAR = OMA    ) 
                end if
                if ( OMP == obs_missingValue_R ) then
                  call fSQL_bind_param( stmt, PARAM_INDEX = 7                    ) 
                else
                  call fSQL_bind_param( stmt, PARAM_INDEX = 7, REAL_VAR = OMP    ) 
                end if
                if ( FGE == obs_missingValue_R ) then
                  call fSQL_bind_param( stmt, PARAM_INDEX = 8                    ) 
                else
                  call fSQL_bind_param( stmt, PARAM_INDEX = 8, REAL_VAR = FGE    )
                end if
                if ( OER == obs_missingValue_R ) then
                  call fSQL_bind_param( stmt, PARAM_INDEX = 9                    ) 
                else
                  call fSQL_bind_param( stmt, PARAM_INDEX = 9, REAL_VAR = OER    ) 
                end if
              case DEFAULT
                call fSQL_bind_param( stmt, PARAM_INDEX = 1, INT_VAR  = obsIdo        )
                call fSQL_bind_param( stmt, PARAM_INDEX = 2, INT_VAR  = obsVarno      )
                call fSQL_bind_param( stmt, PARAM_INDEX = 3, REAL_VAR = PPP           ) 
                call fSQL_bind_param( stmt, PARAM_INDEX = 4, INT_VAR  = vertCoordType )
                if ( obsValue == obs_missingValue_R ) then
                  call fSQL_bind_param( stmt, PARAM_INDEX = 5                         )
                else
                  call fSQL_bind_param( stmt, PARAM_INDEX = 5, REAL_VAR = obsValue    )
                end if 
                call fSQL_bind_param( stmt, PARAM_INDEX = 6, INT_VAR  = obsFlag       )
                if ( OMA == obs_missingValue_R ) then
                  call fSQL_bind_param( stmt, PARAM_INDEX = 7                         ) 
                else 
                  call fSQL_bind_param( stmt, PARAM_INDEX = 7, REAL_VAR = OMA         )
                end if
                if ( OMP == obs_missingValue_R ) then
                  call fSQL_bind_param( stmt, PARAM_INDEX = 8                         ) 
                else
                  call fSQL_bind_param( stmt, PARAM_INDEX = 8, REAL_VAR = OMP         ) 
                end if
                if ( FGE == obs_missingValue_R ) then
                  call fSQL_bind_param( stmt, PARAM_INDEX = 9                         ) 
                else
                  call fSQL_bind_param( stmt, PARAM_INDEX = 9, REAL_VAR = FGE         ) 
                end if
                if ( OER == obs_missingValue_R ) then
                  call fSQL_bind_param( stmt, PARAM_INDEX = 10                        ) 
                else
                  call fSQL_bind_param( stmt, PARAM_INDEX = 10, REAL_VAR = OER        )
                end if 
            end select
            call fSQL_exec_stmt ( stmt )

            numberInsert = numberInsert + 1
          end if    ! llok
        end if

      end do BODY

    end do HEADER

    call fSQL_finalize( stmt )
    call fSQL_commit(db)
    write(*,'(3a,i8)') myName//' FAMILY ---> ' ,trim(familyType), '  NUMBER OF INSERTIONS ----> ', numberInsert

  end subroutine sqlr_insertSqlite


  subroutine sqlr_cleanSqlite(db, fileName)
    !
    ! :Purpose: Remove flagged (bit 11 set) observations in an SQLite file
    !
    implicit none

    ! arguments
    type(fSQL_DATABASE), intent(inout) :: db   ! SQLite file handle
    character(len=*),    intent(in) :: fileName

    ! locals
    character(len=*), parameter :: myName = 'sqlr_cleanSqlite:'
    character(len=*), parameter :: myError = myName //' ERROR: '

    character(len = 128) :: query
    type(fSQL_STATEMENT) :: statement ! prepared statement for SQLite
    type(fSQL_STATUS)    :: status

    call fSQL_open(db, fileName, status)
    if (fSQL_error(status) /= FSQL_OK) then
      write(*,*) myError, fSQL_errmsg(status)
    end if

    ! Mark for deletion all records with bit 11 (2048) set
    query = ' delete from data where flag & 2048;'
    call fSQL_prepare(db, query, statement, status)
    if (fSQL_error(status) /= FSQL_OK) &
      call sqlr_handleError(status, 'thinning fSQL_prepare : ')
    call fSQL_begin(db)
    call fSQL_exec_stmt(statement)
    call fSQL_finalize(statement)
    call fSQL_commit(db)

    write(*,*)'  closed database -->', trim(FileName)
    call fSQL_close( db, status )
  end subroutine sqlr_cleanSqlite


  function getObsFileName(obsFamily, sfFileName_opt, codetype_opt) result(fileName)
    !
    ! :Purpose: Return the part of the observation file name associated
    !           with the type of observation it contains.
    !
    implicit none

    ! arguments:
    character(len=*)           :: obsFamily
    character(len=*), optional :: sfFileName_opt ! fileName acronym used for surface obs file
    integer, optional          :: codetype_opt
    character(len=20) :: fileName

    if ( obsFamily == 'TO' ) then
      if (.not. present(codetype_opt)) then
        call utl_abort('getObsFileName: codetype_opt must be specified for TO family')
      end if

      if ( codtyp_get_name( codeType_opt ) == 'radianceclear' ) then
        fileName  = 'csr'
      else if ( codtyp_get_name( codeType_opt ) == 'mhs' .or. codtyp_get_name( codeType_opt ) == 'amsub' ) then
        fileName = 'to_amsub'
      else if ( codtyp_get_name( codeType_opt ) == 'amsua' ) then
        fileName = 'to_amsua'
      else if ( codtyp_get_name( codeType_opt ) == 'ssmi' ) then
        fileName = 'ssmis'
      else if ( codtyp_get_name( codeType_opt ) == 'crisfsr' ) then
        fileName = 'cris'
      else
        fileName = codtyp_get_name( codeType_opt )
      end if
    else
      if (.not. present(sfFileName_opt)) then
        call utl_abort('getObsFileName: sfFileName_opt must be specified')
      end if
      call up2low( obsFamily, fileName )
      if ( fileName == 'ra' ) fileName = 'radar'
      if ( fileName == 'sf' ) then
        ! use either 'sf' or 'sfc' for filename with surface obs
        fileName = sfFileName_opt
      end if
    end if

  end function getObsFileName


  subroutine sqlr_writeAllSqlDiagFiles( obsdat, sfFileName, onlyAssimObs )
    !
    ! :Purpose: To prepare the writing of obsSpaceData content into SQLite format files
    !  
    implicit none

    ! arguments:
    type(struct_obs)       :: obsdat       ! obsSpaceData object
    character(len=*)       :: sfFileName   ! fileName acronym used for surface obs file
    logical                :: onlyAssimObs ! only write assimilated obs

    ! locals:
    integer                :: familyIndex, codeTypeIndex, fileIndex
    character(len=2)       :: obsFamilyList(50)
    integer                :: obsFamilyListSize
    integer                :: tovsAllCodeTypeListSize, tovsAllCodeTypeList(30)
    integer                :: tovsCodeTypeListSize, tovsCodeTypeList(10)
    integer                :: tovsFileNameListSize
    character(len=20)      :: tovsFileNameList(30)
    character(len=20)      :: fileName

    ! ensure all mpi tasks have same list of common obs family names
    call getObsFamilyListMpiGlobal(obsdat, obsFamilyListSize, obsFamilyList)

    ! get list of all possible tovs codetype values and unique list of corresponding filenames
    call tvs_getAllIdBurpTovs(tovsAllCodeTypeListSize, tovsAllCodeTypeList)
    write(*,*) 'tovsAllCodeTypeListSize = ', tovsAllCodeTypeListSize
    write(*,*) 'tovsAllCodeTypeList = ', tovsAllCodeTypeList(1:tovsAllCodeTypeListSize)
    
    tovsFileNameListSize = 0
    tovsFileNameList(:) = 'XXXXX'
    do codeTypeIndex = 1, tovsAllCodeTypeListSize
      fileName = getObsFileName('TO', codeType_opt=tovsAllCodeTypeList(codeTypeIndex))
      if ( all(tovsFileNameList(:) /= fileName) ) then
        tovsFileNameListSize = tovsFileNameListSize + 1
        tovsFileNameList(tovsFileNameListSize) = fileName
      end if
    end do
    write(*,*) 'tovsFileNameListSize = ', tovsFileNameListSize
    write(*,*) 'tovsFileNameList = ', tovsFileNameList(1:tovsFileNameListSize)
    
    do familyIndex = 1, obsFamilyListSize

      write(*,*) 'sqlr_writeAllSqlDiagFiles: Family = ', familyIndex, obsFamilyList(familyIndex)

      if ( obsFamilyList(familyIndex) == 'TO' ) then

        do fileIndex = 1, tovsFileNameListSize
          fileName = tovsFileNameList(fileIndex)
          write(*,*) 'tovs filename = ', fileName

          ! get list of codetypes associated with this filename
          tovsCodeTypeListSize = 0
          tovsCodeTypeList(:) = MPC_missingValue_INT
          do codeTypeIndex = 1, tovsAllCodeTypeListSize
            if (fileName == getObsFileName('TO', codeType_opt=tovsAllCodeTypeList(codeTypeIndex))) then
              tovsCodeTypeListSize = tovsCodeTypeListSize + 1
              tovsCodeTypeList(tovsCodeTypeListSize) = tovsAllCodeTypeList(codeTypeIndex)
            end if
          end do

          write(*,*) 'tovsCodeTypeListSize = ', tovsCodeTypeListSize
          write(*,*) 'tovsCodeTypeList = ', tovsCodeTypeList(1:tovsCodeTypeListSize) 
          call sqlr_writeSqlDiagFile(obsdat, 'TO', onlyAssimObs, &
                                     tovsFileNameList(fileIndex), &
                                     tovsCodeTypeList(1:tovsCodeTypeListSize)) 
        end do

      else

        fileName = getObsFileName(obsFamilyList(familyIndex), sfFileName_opt=sfFileName)
        call sqlr_writeSqlDiagFile(obsdat, obsFamilyList(familyIndex),  &
                                   onlyAssimObs, fileName) 

      end if   
      
    end do

  end subroutine sqlr_writeAllSqlDiagFiles


  subroutine getObsFamilyListMpiGlobal(obsdat, obsFamilyListSizeCommon,  &
                                       obsFamilyListCommon)
    !
    ! :Purpose: Obtain a common set of obs family names over all mpi tasks
    !
    implicit none
      
    ! arguments:
    type(struct_obs) :: obsdat
    integer          :: obsFamilyListSizeCommon
    character(len=*) :: obsFamilyListCommon(:)

    ! locals:
    integer                       :: headerIndex, familyIndex, charIndex, procIndex, nsize, ierr
    integer                       :: obsFamilyListSizeMpiLocal, obsFamilyListSizeMaxMpiLocal, obsFamilyListSizeMax
    character(len=2), allocatable :: obsFamilyListMpiLocal(:)
    character(len=2), allocatable :: obsFamilyListMpiGlobal(:,:)
    character(len=2)              :: currentObsFamily
    integer, allocatable          :: intObsFamilyListMpiLocal(:,:)
    integer, allocatable          :: intObsFamilyListMpiGlobal(:,:,:)
    integer, allocatable          :: allObsFamilyListSizeMpiLocal(:)

    obsFamilyListSizeMax = size(obsFamilyListCommon)
    write(*,*) 'obsFamilyListSizeMax =', obsFamilyListSizeMax

    ! get family list for this mpi task
    obsFamilyListSizeMpiLocal = 0
    allocate(obsFamilyListMpiLocal(obsFamilyListSizeMax))
    obsFamilyListMpiLocal(:) = 'XX'
    HEADER: do headerIndex = 1, obs_numHeader( obsdat )
      currentObsFamily = obs_getFamily( obsdat, headerIndex ) 
      if ( any( obsFamilyListMpiLocal(:) == currentObsFamily ) ) cycle HEADER
      obsFamilyListSizeMpiLocal = obsFamilyListSizeMpiLocal + 1
      obsFamilyListMpiLocal( obsFamilyListSizeMpiLocal ) = currentObsFamily
      write(*,*) 'add the family: ', currentObsFamily
    end do HEADER
    write(*,*) 'obsFamilyListSizeMpiLocal =', obsFamilyListSizeMpiLocal
    write(*,*) 'obsFamilyListMpiLocal = ', obsFamilyListMpiLocal(1:obsFamilyListSizeMpiLocal)

    allocate(allObsFamilyListSizeMpiLocal(mpi_nprocs))
    call rpn_comm_allgather(obsFamilyListSizeMpiLocal,    1, 'mpi_integer',  &
                            allObsFamilyListSizeMpiLocal, 1, 'mpi_integer', 'GRID', ierr)
    call rpn_comm_allreduce(obsFamilyListSizeMpiLocal, obsFamilyListSizeMaxMpiLocal,1,'mpi_integer','mpi_max','GRID',ierr)

    ! convert local family list from characters to integers
    allocate(intObsFamilyListMpiLocal(len(currentObsFamily),obsFamilyListSizeMaxMpiLocal))
    intObsFamilyListMpiLocal(:,:)=0
    do familyIndex = 1, obsFamilyListSizeMpiLocal
      do charIndex = 1, len(currentObsFamily)
        intObsFamilyListMpiLocal(charIndex,familyIndex) =  &
             iachar(obsFamilyListMpiLocal(familyIndex)(charIndex:charIndex))
      end do
    end do

    ! communicate obs family list to all mpi tasks as integers
    allocate(intObsFamilyListMpiGlobal(len(currentObsFamily),obsFamilyListSizeMaxMpiLocal,mpi_nprocs))
    nsize = size(intObsFamilyListMpiLocal)
    call rpn_comm_allgather(intObsFamilyListMpiLocal,  nsize, 'mpi_integer',  &
                            intObsFamilyListMpiGlobal, nsize, 'mpi_integer', 'GRID', ierr)

    ! convert global family lists from integers to characters
    allocate(obsFamilyListMpiGlobal(obsFamilyListSizeMaxMpiLocal,mpi_nprocs))
    obsFamilyListMpiGlobal(:,:) = 'XX'
    do procIndex = 1, mpi_nprocs
      do familyIndex = 1, allObsFamilyListSizeMpiLocal(procIndex)
        do charIndex=1,len(currentObsFamily)
          obsFamilyListMpiGlobal(familyIndex,procIndex)(charIndex:charIndex) =  &
               achar(intObsFamilyListMpiGlobal(charIndex,familyIndex,procIndex))
        end do
      end do
      write(*,*) 'obsFamilyListMpiGlobal = ', procIndex,  &
           obsFamilyListMpiGlobal(1:allObsFamilyListSizeMpiLocal(procIndex),procIndex)
    end do

    ! construct single common list of families to be used for all mpi tasks
    obsFamilyListCommon(:) = 'YY'
    obsFamilyListSizeCommon = 0
    do procIndex = 1, mpi_nprocs
      FAMILY: do familyIndex = 1, obsFamilyListSizeMaxMpiLocal
        if (obsFamilyListMpiGlobal(familyIndex,procIndex) == 'XX') cycle FAMILY
        if (any( obsFamilyListCommon(:) == obsFamilyListMpiGlobal(familyIndex,procIndex) )) cycle FAMILY
        obsFamilyListSizeCommon = obsFamilyListSizeCommon + 1
        obsFamilyListCommon(obsFamilyListSizeCommon) = obsFamilyListMpiGlobal(familyIndex,procIndex)
      end do FAMILY
    end do
    write(*,*) 'obsFamilyListSizeCommon = ', obsFamilyListSizeCommon
    write(*,*) 'obsFamilyListCommon = ', obsFamilyListCommon(1:obsFamilyListSizeCommon)

    deallocate(allObsFamilyListSizeMpiLocal)
    deallocate(obsFamilyListMpiGlobal)
    deallocate(intObsFamilyListMpiGlobal)
    deallocate(intObsFamilyListMpiLocal)
    deallocate(obsFamilyListMpiLocal)

  end subroutine getObsFamilyListMpiGlobal


  subroutine sqlr_writeSqlDiagFile( obsdat, obsFamily, onlyAssimObs, &
                                    instrumentFileName, codeTypeList_opt )
    !
    ! :Purpose: To write the obsSpaceData content into SQLite format files
    !
    implicit none

    ! arguments
    type(struct_obs)  :: obsdat
    character(len=*)  :: obsFamily
    logical           :: onlyAssimObs
    character(len=*)  :: instrumentFileName
    integer, optional :: codeTypeList_opt(:)

    ! locals
    type(fSQL_DATABASE)    :: db                   ! type for SQLIte  file handle
    type(fSQL_STATEMENT)   :: stmtData, stmtHeader ! type for precompiled SQLite statements
    type(fSQL_STATUS)      :: stat                 ! type for error status
    integer                :: obsVarno, obsFlag, ASS, vertCoordType, codeType, date, time, idObs, idData 
    real                   :: obsValue, OMA, OMP, OER, FGE, PPP, lon, lat, altitude
    real                   :: ensInnovStdDev, ensObsErrStdDev, zhad
    integer                :: numberInsertions, numHeaders, headerIndex, bodyIndex, obsNlv, obsRln
    character(len = 512)   :: queryData, queryHeader, queryCreate 
    character(len = 12 )   :: idStation
    character(len=*), parameter :: myName = 'sqlr_writeSqlDiagFile:'
    character(len=*), parameter :: myWarning = myName //' WARNING: '
    character(len=*), parameter :: myError   = myName //' ERROR: '
    character(len=30)      :: fileNameExtention
    character(len=256)     :: fileName, fileNameDir
    character(len=4)       :: cmyidx, cmyidy
    logical                :: writeHeader

    ! determine initial idData,idObs to ensure unique values across mpi tasks
    call getInitialIdObsData(obsDat, obsFamily, idObs, idData, codeTypeList_opt)

    ! return if this mpi task does not have any observations of this family
    if (trim(instrumentFileName) == 'XXXXX') return

    ! check if any obs exist for this file, if not return
    numHeaders = 0
    call obs_set_current_header_list( obsdat, obsFamily )
    HEADERCOUNT: do
      headerIndex = obs_getHeaderIndex( obsdat )
      if ( headerIndex < 0 ) exit HEADERCOUNT
      if ( present( codeTypeList_opt ) ) then
        codeType  = obs_headElem_i( obsdat, OBS_ITY, headerIndex )
        if ( all( codeTypeList_opt(:) /= codeType ) ) cycle HEADERCOUNT
      end if
      numHeaders = numHeaders + 1
    end do HEADERCOUNT
    if (numHeaders == 0) return

    fileNameDir = trim(ram_getRamDiskDir())
    if ( fileNameDir == ' ' ) &
    write(*,*) myWarning//' The program may be slow creating many sqlite files in the same directory.'
    write(*,*) myWarning//' Please, use the ram disk option prior to MIDAS run!'

    if ( obs_mpiLocal( obsdat ) ) then
      write(cmyidy,'(I4.4)') ( mpi_myidy + 1 )
      write(cmyidx,'(I4.4)') ( mpi_myidx + 1 )
      fileNameExtention  = trim(cmyidx) // '_' // trim( cmyidy )
    else
      if ( mpi_myid > 0 ) return
      fileNameExtention = ' '
    end if

    fileName = trim(fileNameDir) // 'obs/dia' // trim(instrumentFileName) // '_' // trim( fileNameExtention )

    write(*,*) myName//' Creating file: ', trim(fileName)
    call fSQL_open( db, fileName, stat )
    if ( fSQL_error( stat ) /= FSQL_OK ) write(*,*) myError//'fSQL_open: ', fSQL_errmsg( stat ),' filename: '//trim(fileName)

    ! Create the tables HEADER and DATA
    queryCreate = 'create table header (id_obs integer primary key, id_stn varchar(50), lat real, lon real, &
                   &codtyp integer, date integer, time integer, elev real); &
                   &create table data (id_data integer primary key, id_obs integer, varno integer, vcoord real, &
                   &vcoord_type integer, obsvalue real, flag integer, oma real, ompt real, oma0 real, omp real, &
                   &an_error real, fg_error real, obs_error real, sigi real, sigo real, zhad real);'
    call fSQL_do_many( db, queryCreate, stat )
    if ( fSQL_error(stat) /= FSQL_OK ) call sqlr_handleError( stat, 'fSQL_do_many with query: '//trim(queryCreate) )

    queryData = 'insert into data (id_data, id_obs, varno, vcoord, vcoord_type, obsvalue, flag, oma, oma0, ompt, fg_error, obs_error, sigi, sigo, zhad) values(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);'
    queryHeader = ' insert into header (id_obs, id_stn, lat, lon, date, time, codtyp, elev ) values(?,?,?,?,?,?,?,?); '

    write(*,*) myName//' Insert query Data   = ', trim( queryData )
    write(*,*) myName//' Insert query Header = ', trim( queryHeader )

    call fSQL_begin(db)
    call fSQL_prepare( db, queryData, stmtData, stat )
    if ( fSQL_error(stat) /= FSQL_OK ) call sqlr_handleError(stat, 'fSQL_prepare : ')
    call fSQL_prepare( db, queryHeader, stmtHeader, stat )
    if ( fSQL_error(stat) /= FSQL_OK ) call sqlr_handleError(stat, 'fSQL_prepare : ')

    numberInsertions = 0

    call obs_set_current_header_list( obsdat, obsFamily )
    HEADER: do

      headerIndex = obs_getHeaderIndex( obsdat )
      if ( headerIndex < 0 ) exit HEADER
        
      codeType  = obs_headElem_i( obsdat, OBS_ITY, headerIndex )
      if ( present( codeTypeList_opt ) ) then
        if ( all( codeTypeList_opt(:) /= codeType ) ) cycle HEADER
      end if

      obsRln    = obs_headElem_i( obsdat, OBS_RLN, headerIndex )
      obsNlv    = obs_headElem_i( obsdat, OBS_NLV, headerIndex )
      idStation = obs_elem_c    ( obsdat, 'STID' , headerIndex ) 
      altitude  = obs_headElem_r( obsdat, OBS_ALT, headerIndex )      
      lon       = obs_headElem_r( obsdat, OBS_LON, headerIndex ) * MPC_DEGREES_PER_RADIAN_R8
      lat       = obs_headElem_r( obsdat, OBS_LAT, headerIndex ) * MPC_DEGREES_PER_RADIAN_R8
      if (  lon > 180. ) lon = lon - 360.
      date      = obs_headElem_i( obsdat, OBS_DAT, headerIndex )
      time      = obs_headElem_i( obsdat, OBS_ETM, headerIndex ) * 100.

      ! check if at least one observation will be written, if not skip this header
      if (onlyAssimObs) then
        writeHeader = .false.
        BODYCHECK: do bodyIndex = obsRln, obsNlv + obsRln -1
          OMP = obs_bodyElem_r( obsdat, OBS_OMP , bodyIndex )
          FGE = obs_bodyElem_r( obsdat, OBS_HPHT, bodyIndex )
          ASS = obs_bodyElem_i( obsdat, OBS_ASS , bodyIndex )
          if ( (ASS == obs_assimilated) .and.     &
               (OMP /= obs_missingValue_R) .and.  &
               (FGE /= obs_missingValue_R) ) writeHeader = .true.
        end do BODYCHECK
        if (.not. writeHeader) cycle HEADER
      end if

      idObs = idObs + 1
      call fSQL_bind_param( stmtHeader, PARAM_INDEX = 1, INT_VAR  = idObs     )
      call fSQL_bind_param( stmtHeader, PARAM_INDEX = 2, CHAR_VAR = idStation )
      call fSQL_bind_param( stmtHeader, PARAM_INDEX = 3, REAL_VAR = lat       ) 
      call fSQL_bind_param( stmtHeader, PARAM_INDEX = 4, REAL_VAR = lon       ) 
      call fSQL_bind_param( stmtHeader, PARAM_INDEX = 5, INT_VAR  = date      ) 
      call fSQL_bind_param( stmtHeader, PARAM_INDEX = 6, INT_VAR  = time      ) 
      call fSQL_bind_param( stmtHeader, PARAM_INDEX = 7, INT_VAR  = codeType  ) 
      call fSQL_bind_param( stmtHeader, PARAM_INDEX = 8, REAL_VAR = altitude  ) 
      call fSQL_exec_stmt ( stmtHeader )

      BODY: do bodyIndex = obsRln, obsNlv + obsRln -1
         
        obsVarno      = obs_bodyElem_i( obsdat, OBS_VNM , bodyIndex )
        obsFlag       = obs_bodyElem_i( obsdat, OBS_FLG , bodyIndex )
        vertCoordType = obs_bodyElem_i( obsdat, OBS_VCO , bodyIndex )
        obsValue      = obs_bodyElem_r( obsdat, OBS_VAR , bodyIndex )
        OMA           = obs_bodyElem_r( obsdat, OBS_OMA , bodyIndex )
        OMP           = obs_bodyElem_r( obsdat, OBS_OMP , bodyIndex )
        OER           = obs_bodyElem_r( obsdat, OBS_OER , bodyIndex )
        FGE           = obs_bodyElem_r( obsdat, OBS_HPHT, bodyIndex )
        PPP           = obs_bodyElem_r( obsdat, OBS_PPP , bodyIndex )
        ASS           = obs_bodyElem_i( obsdat, OBS_ASS , bodyIndex )

        ! skip obs if it was not assimilated
        if (onlyAssimObs) then
          if ( (ASS /= obs_assimilated) .or.     &
               (OMP == obs_missingValue_R) .or.  &
               (FGE == obs_missingValue_R) ) cycle BODY
        end if

        if ( obs_columnActive_RB(obsdat, OBS_SIGI) ) then
          ensInnovStdDev = obs_bodyElem_r(obsdat, OBS_SIGI, bodyIndex)
        else
          ensInnovStdDev = obs_missingValue_R
        end if
        if ( obs_columnActive_RB(obsdat, OBS_SIGO) ) then
          ensObsErrStdDev = obs_bodyElem_r(obsdat, OBS_SIGO, bodyIndex)
        else
          ensObsErrStdDev = obs_missingValue_R
        end if
        if ( obs_columnActive_RB(obsdat, OBS_ZHA) ) then
          zhad = obs_bodyElem_r(obsdat, OBS_ZHA, bodyIndex)
        else
          zhad = obs_missingValue_R
        end if

        select case( obsFamily )
          case ( 'UA', 'AI', 'SW' )
            if ( vertCoordType == 2 ) vertCoordType = 7004
          case ( 'RO' )
            vertCoordType = 7007
          case ( 'PR' )
            vertCoordType = 7006
            PPP = PPP - altitude
          case ( 'TO' )
            vertCoordType = 5042
            if( codeType == 164 .or. codeType == 181 .or. codeType == 182 ) vertCoordType = 2150
          case ( 'SF', 'SC', 'GP' )
            vertCoordType = MPC_missingValue_INT 
        end select

        ! insert order: id_obs,varno,vcoord,vcoord_type,obsvalue,flag,oma,oma0,ompt,fg_error,obs_error,sigi,sigo
        idData = idData + 1
        call fSQL_bind_param( stmtData, PARAM_INDEX = 1, INT_VAR  = idData        )
        call fSQL_bind_param( stmtData, PARAM_INDEX = 2, INT_VAR  = idObs         )
        call fSQL_bind_param( stmtData, PARAM_INDEX = 3, INT_VAR  = obsVarno      )
        call fSQL_bind_param( stmtData, PARAM_INDEX = 4, REAL_VAR = PPP           )
        if ( vertCoordType == MPC_missingValue_INT ) then
          call fSQL_bind_param( stmtData, PARAM_INDEX = 5                         ) 
        else
          call fSQL_bind_param( stmtData, PARAM_INDEX = 5, INT_VAR  = vertCoordType ) 
        end if
        call fSQL_bind_param( stmtData, PARAM_INDEX = 6, REAL_VAR = obsValue      ) 
        call fSQL_bind_param( stmtData, PARAM_INDEX = 7, INT_VAR  = obsFlag       )
        if ( OMA == obs_missingValue_R ) then
          call fSQL_bind_param( stmtData, PARAM_INDEX = 8                         ) 
          call fSQL_bind_param( stmtData, PARAM_INDEX = 9                         ) 
        else
          call fSQL_bind_param( stmtData, PARAM_INDEX = 8, REAL_VAR = OMA         )
          call fSQL_bind_param( stmtData, PARAM_INDEX = 9, REAL_VAR = OMA         )
        end if
        if ( OMP == obs_missingValue_R ) then
          call fSQL_bind_param( stmtData, PARAM_INDEX = 10                         ) 
        else
          call fSQL_bind_param( stmtData, PARAM_INDEX = 10, REAL_VAR = OMP         )
        end if
        if ( FGE == obs_missingValue_R ) then
          call fSQL_bind_param( stmtData, PARAM_INDEX = 11                         ) 
        else
          call fSQL_bind_param( stmtData, PARAM_INDEX = 11, REAL_VAR = FGE         )
        end if
        if ( OER == obs_missingValue_R ) then
          call fSQL_bind_param( stmtData, PARAM_INDEX = 12                        ) 
        else
          call fSQL_bind_param( stmtData, PARAM_INDEX = 12, REAL_VAR = OER        )
        end if 
        if ( ensInnovStdDev == obs_missingValue_R ) then
          call fSQL_bind_param( stmtData, PARAM_INDEX = 13                        ) 
        else
          call fSQL_bind_param( stmtData, PARAM_INDEX = 13, REAL_VAR = ensInnovStdDev )
        end if 
        if ( ensObsErrStdDev == obs_missingValue_R ) then
          call fSQL_bind_param( stmtData, PARAM_INDEX = 14                        ) 
        else
          call fSQL_bind_param( stmtData, PARAM_INDEX = 14, REAL_VAR = ensObsErrStdDev )
        end if 
        if ( zhad == obs_missingValue_R ) then
          call fSQL_bind_param( stmtData, PARAM_INDEX = 15                        ) 
        else
          call fSQL_bind_param( stmtData, PARAM_INDEX = 15, REAL_VAR = zhad )
        end if 

        call fSQL_exec_stmt ( stmtData )

        numberInsertions = numberInsertions + 1

      end do BODY
     
    end do HEADER

    write(*,*) myName// ' Observation Family: ', obsFamily,', number of insertions: ', numberInsertions
    call fSQL_finalize( stmtData )
    call fSQL_commit(db)
    call fSQL_close( db, stat )

  end subroutine sqlr_writeSqlDiagFile


  subroutine getInitialIdObsData(obsDat, obsFamily, idObs, idData, codeTypeList_opt)
    !
    ! :Purpose: Compute initial value for idObs and idData that will ensure
    !           unique values over all mpi tasks
    !
    implicit none

    ! arguments:
    type(struct_obs)  :: obsdat
    character(len=*)  :: obsFamily    
    integer           :: idObs, idData
    integer, optional :: codeTypeList_opt(:)

    ! locals:
    integer                :: headerIndex, numHeader, numBody, codeType, ierr
    integer, allocatable   :: allNumHeader(:), allNumBody(:)

    numHeader = 0
    numBody = 0
    call obs_set_current_header_list( obsdat, obsFamily )
    HEADERCOUNT: do
      headerIndex = obs_getHeaderIndex( obsdat )
      if ( headerIndex < 0 ) exit HEADERCOUNT
      if ( present( codeTypeList_opt ) ) then
        codeType  = obs_headElem_i( obsdat, OBS_ITY, headerIndex )
        if ( all( codeTypeList_opt(:) /= codeType ) ) cycle HEADERCOUNT
      end if
      numHeader = numHeader + 1
      numBody = numBody + obs_headElem_i( obsdat, OBS_NLV, headerIndex )
    end do HEADERCOUNT
    allocate(allNumHeader(mpi_nprocs))
    allocate(allNumBody(mpi_nprocs))
    call rpn_comm_allgather(numHeader,1,'mpi_integer',       &
                            allNumHeader,1,'mpi_integer','GRID',ierr)
    call rpn_comm_allgather(numBody,1,'mpi_integer',       &
                            allNumBody,1,'mpi_integer','GRID',ierr)
    if (mpi_myid > 0) then
      idObs = sum(allNumHeader(1:mpi_myid))
      idData = sum(allNumBody(1:mpi_myid))
    else
      idObs = 0
      idData = 0
    end if
    deallocate(allNumHeader)
    deallocate(allNumBody)
  end subroutine getInitialIdObsData

end module sqliteRead_mod
