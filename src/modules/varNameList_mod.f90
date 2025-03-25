
module varNameList_mod
  ! MODULE varNameList_mod (prefix='vnl' category='7. Low-level data objects')
  !
  !:Purpose:  Contains a list of all possible variable names that can be used
  !           as analysis variables along with additional information for each
  !           and procedures for accessing this information
  !
  use bufr_mod
  use midasMpi_mod
  use utilities_mod
  use MathPhysConstants_mod
  use netcdf

  implicit none
  save
  private

  ! Public variables (parameters)
  public :: vnl_numvarmax3D, vnl_numvarmax2D, vnl_numvarmaxOther, vnl_numvarmax
  public :: vnl_varNameList3D, vnl_varNameList2D, vnl_varNameListOther, vnl_varNameList
  public :: vnl_numvarmaxCloud, vnl_varNameListCloud

  ! Public procedures
  public :: vnl_varListIndex3d, vnl_varListIndex2d, vnl_varListIndexOther
  public :: vnl_varListIndex, vnl_varnameFromVarnum, vnl_varnameIsValid
  public :: vnl_varLevelFromVarname, vnl_varLevelFromVarnum
  public :: vnl_varKindFromVarname, vnl_varnumFromVarname
  public :: vnl_varNamesFromExistList, vnl_varMassFromVarNum, vnl_varMassFromVarName
  public :: vnl_isPhysicsVar, vnl_isCloudVar
  public :: vnl_addToVarNames
  public :: vnl_varNamePresentInFile, vnl_varNameNetCDF

  ! These private parameters permit side-stepping a conflict with the Sphinx documenter,
  ! and an infinite loop
  integer, parameter          :: VNLnumvarmax3D    = 52
  integer, parameter          :: VNLnumvarmax2D    = 39
  integer, parameter          :: VNLnumvarmaxOther =  7
  integer, parameter          :: VNLnumvarmaxCloud =  5

  integer, parameter          :: vnl_numvarmax3D    = VNLnumvarmax3D
  integer, parameter          :: vnl_numvarmax2D    = VNLnumvarmax2D
  integer, parameter          :: vnl_numvarmaxOther = VNLnumvarmaxOther
  integer, parameter          :: vnl_numvarmaxCloud = VNLnumvarmaxCloud

  character(len=4), parameter :: vnl_varNameList3D(vnl_numvarmax3D) = (/                         &
                                 'UU  ','VV  ','Z_T ','Z_M ','P_T ','P_M ',                      &
                                 'TT  ','HU  ','LQ  ','ES  ','VT  ',                             &
                                 'PP  ','CC  ','UC  ','UT  ','TB  ','DW  ','QR  ','DD  ',        &
                                 'TO3 ','O3L ','TCH4','TCO2','TCO ','TNO2','TN2O','THCH',        &
                                 'TSO2','TNH3','AF  ','AC  ','TNO ','ALFA','VIS ','LVIS',        &
                                 'HR  ','TD  ','ALFT','UV  ','LWCR','IWCR','QC  ','CH4L',        &
                                 'N2OL','UUW ','VVW ','TM  ','SALW','ALFO','RF  ','SF  ',        &
                                 'CLDR' /)

  character(len=4), parameter :: varLevelList3D(vnl_numvarmax3D)     = (/                        &
                                 'MM',  'MM',  'TH',  'MM',  'TH',  'MM',                        &
                                 'TH',  'TH',  'TH',  'TH',  'TH',                               &
                                 'MM',  'MM',  'MM',  'TH',  'TH',  'TH',  'MM',  'MM',          &
                                 'TH',  'TH',  'TH',  'TH',  'TH',  'TH',  'TH',  'TH',          &
                                 'TH',  'TH',  'TH',  'TH',  'TH',  'MM',  'TH',  'TH',          &
                                 'TH',  'TH',  'TH',  'MM',  'TH',  'TH',  'TH',  'TH',          &
                                 'TH',  'DP',  'DP',  'DP',  'DP',  'DP',  'TH',  'TH',          &
                                 'TH' /)

  character(len=2), parameter :: varKindList3D(vnl_numvarmax3D)     = (/                         &
                                 'MT',  'MT',  'MT',  'MT',  'MT',  'MT',                        &
                                 'MT',  'MT',  'MT',  'MT',  'MT',                               &
                                 'MT',  'MT',  'MT',  'MT',  'MT',  'MT',  'MT',  'MT',          &
                                 'CH',  'CH',  'CH',  'CH',  'CH',  'CH',  'CH',  'CH',          &
                                 'CH',  'CH',  'CH',  'CH',  'CH',  'MT',  'MT',  'MT',          &
                                 'MT',  'MT',  'MT',  'MT',  'MT',  'MT',  'MT',  'CH',          &
                                 'CH',  'OC',  'OC',  'OC',  'OC',  'OC',  'MT',  'MT',          &
                                 'MT' /)

  character(len=4), parameter :: vnl_varNameList2D(vnl_numvarmax2D) = (/ &
                                 'P0  ','TG  ','UP  ','PB  ','ECO ','ENO2','EHCH','ESO2','ENH3', &
                                 'GL  ','WGE ','BIN ','MG  ','SSH ','QI1 ','QO1 ','STOR','ALFS', &
                                 'PN  ','PR  ','LPR ','I2  ','I3  ','I4  ','I5  ','I6  ','I8  ', &
                                 'DN  ','FB  ','FI  ','MSKC','LZS ','WT  ','LG  ','VF  ','DSLO', &
                                 'P0LS','MELS','SSS '/)

  character(len=4), parameter :: varLevelList2D(vnl_numvarmax2D) = (/    &
                                 'SF',  'SF',  'SF',  'SF',  'SF',  'SF',  'SF',  'SF',  'SF',  &
                                 'SS',  'SF',  'SF',  'SF',  'SS',  'SF',  'SF',  'SF',  'SF',  &
                                 'SF',  'SF',  'SF',  'SF',  'SF',  'SF',  'SF',  'SF',  'SF',  &
                                 'SF',  'SF',  'SF',  'SF',  'SF',  'SF',  'SS',  'SS',  'SS',  &
                                 'SF',  'SF',  'SS'/)

  character(len=2), parameter :: varKindList2D(vnl_numvarmax2D) = (/     &
                                 'MT',  'MT',  'MT',  'MT',  'CH',  'CH',  'CH',  'CH',  'CH', &
                                 'OC',  'MT',  'MT',  'MT',  'OC',  'HY',  'HY',  'HY',  'HY', &
                                 'MT',  'MT',  'MT',  'MT',  'MT',  'MT',  'MT',  'MT',  'MT', &
                                 'MT',  'MT',  'MT',  'MT',  'HY',  'MT',  'OC',  'OC',  'OC', &
                                 'MT',  'MT',  'OC'/)

  character(len=4), parameter :: vnl_varNameListOther(vnl_numvarmaxOther) = (/ &
                                 'I0  ','I1  ','I7  ','I9  ','SD  ','AL  ','EMMW'/)

  character(len=4), parameter :: varLevelListOther(vnl_numvarmaxOther) = (/    &
                                 'OT',  'OT',  'OT',  'OT',  'OT',  'OT',  'OT'  /)

  character(len=2), parameter :: varKindListOther(vnl_numvarmaxOther) = (/     &
                                 'LD',  'LD',  'LD',  'LD',  'LD',  'LD',  'LD'  /) ! LD = Land

  character(len=4), parameter :: vnl_varNameListCloud(vnl_numvarmaxCloud) = (/ &
                                 'LWCR', 'IWCR', 'RF  ', 'SF  ', 'CLDR' /)

  integer, parameter          :: vnl_numvarmax = VNLnumvarmax3D + VNLnumvarmax2D + VNLnumvarmaxOther

  character(len=4), parameter :: vnl_varNameList(vnl_numvarmax) =  &
       (/ vnl_varNameList3D, vnl_varNameList2D, vnl_varNameListOther /)
  character(len=4), parameter :: varLevelList   (vnl_numvarmax) =  &
       (/ varLevelList3D   , varLevelList2D   , varLevelListOther    /)
  character(len=2), parameter :: varKindList    (vnl_numvarmax) =  &
       (/ varKindList3D    , varKindList2D    , varKindListOther     /)

  contains

   !--------------------------------------------------------------------------
   ! vnl_varListIndex3d
   !--------------------------------------------------------------------------
    function vnl_varListIndex3d(varName) result(listIndex)
      !
      ! :Purpose: To get the 3d list index from the variable name
      !
      implicit none

      ! Arguments:
      character(len=*), intent(in) :: varName
      ! Result:
      integer                      :: listIndex
      
      ! Locals:
      integer                      :: jvar

      listIndex=-1
      do jvar=1,vnl_numvarmax3D
        if( varName == vnl_varNameList3d(jvar)) then
          listIndex=jvar
          exit
        end if
      end do

      if( listIndex <= 0 ) then
        call utl_abort('vnl_varListIndex3D: Unknown variable name! ' // varName)
      end if

    end function vnl_varListIndex3d

   !--------------------------------------------------------------------------
   ! vnl_varListIndex2d
   !--------------------------------------------------------------------------
    function vnl_varListIndex2d(varName) result(listIndex)
      !
      ! :Purpose: To get the 2d list index from the variable name
      !
      implicit none

      ! Arguments:
      character(len=*), intent(in) :: varName
      ! Result:
      integer                      :: listIndex

      ! Locals:
      integer                      :: jvar

      listIndex=-1
      do jvar = 1, vnl_numvarmax2D
        if( varName == vnl_varNameList2d(jvar) ) then 
          listIndex=jvar
          exit
        end if
      end do

      if( listIndex <= 0 ) then
        call utl_abort('vnl_varListIndex2D: Unknown variable name! ' // varName)
      end if

    end function vnl_varListIndex2d

   !--------------------------------------------------------------------------
   ! vnl_varListIndexOther
   !--------------------------------------------------------------------------
    function vnl_varListIndexOther(varName) result(listIndex)
      !
      ! :Purpose: To get the "Other" list index from the variable name
      !
      implicit none

      ! Arguments:
      character(len=*), intent(in) :: varName
      ! Result:
      integer                      :: listIndex

      ! Locals:
      integer                      :: jvar

      listIndex=-1
      do jvar = 1, vnl_numvarmaxOther
        if( varName == vnl_varNameListOther(jvar) ) then 
          listIndex=jvar
          exit
        end if
      end do

      if(listIndex <= 0) then
        call utl_abort('vnl_varListIndexOther: Unknown variable name! ' // varName)
      end if

    end function vnl_varListIndexOther

   !--------------------------------------------------------------------------
   ! vnl_varListIndex
   !--------------------------------------------------------------------------
    function vnl_varListIndex(varName) result(listIndex)
      !
      ! :Purpose: To get the varlist index from the variable name
      !

      implicit none

      ! Arguments:
      character(len=*), intent(in) :: varName
      ! Result:
      integer                      :: listIndex

      ! Local:
      integer                      :: jvar

      listIndex=-1
      do jvar=1,vnl_numvarmax
        if(varName == vnl_varNameList(jvar)) then 
          listIndex=jvar
          exit
        end if
      end do

      if(listIndex <= 0) then
        call utl_abort('vnl_varListIndex: Unknown variable name! ' // varName)
      end if

    end function vnl_varListIndex

    !--------------------------------------------------------------------------
    ! vnl_varnameIsValid
    !--------------------------------------------------------------------------
    function vnl_varnameIsValid(varName) result(isValid)
      !
      ! :Purpose: Check if the supplied variable name is known by MIDAS.
      !
      implicit none
      
      ! Arguments:
      character(len=*), intent(in) :: varName
      ! Result:
      logical                      :: isValid
      
      ! Local::
      integer                      :: varIndex

      isValid = .false.
      do varIndex = 1, vnl_numvarmax
        if(varName == vnl_varNameList(varIndex)) then 
          isValid = .true.
          exit
        end if
      end do

    end function vnl_varnameIsValid

    !--------------------------------------------------------------------------
    ! vnl_varnameFromVarnum
    !--------------------------------------------------------------------------
    function vnl_varnameFromVarnum(varNumber, varNumberChm_opt, modelName_opt) result(varName)
      !
      ! :Purpose: To get the variable name from the variable number
      !
      implicit none

      ! Arguments:
      integer,                    intent(in) :: varNumber
      integer,          optional, intent(in) :: varNumberChm_opt
      character(len=*), optional, intent(in) :: modelName_opt
      ! Result:
      character(len=4)    :: varName
      character(len=8)    :: modelName

      if (present(modelName_opt)) then
        modelName = trim(modelName_opt)
      else
        modelName = 'GEM'
      end if

      varName = '    '
      select case (varNumber)
        case ( BUFR_NEUU, BUFR_NEUS, BUFR_NEAL )
          varName='UU'
        case( BUFR_NEVV, BUFR_NEVS )
          varName='VV'
        case( BUFR_NETT, BUFR_NETS )
          varName='TT'
        case( BUFR_NEDZ, BUFR_NEGZ )
          varName='Z_T'
        case( BUFR_NEHU, BUFR_NEHS, BUFR_NEES, BUFR_NESS )
          varName='HU'
        case( BUFR_NEPS, BUFR_NEPN )
          varName='P0'
        case ( BUFR_NERF, BUFR_NEBD, BUFR_NEZD )
          varName='TT'   ! temporarily associate refractivity and ZTD with temperature
        case ( BUFR_NEDW )
          varName='DW'
        case ( BUFR_SST )
          varname='TG'
        case ( BUFR_ICEC, BUFR_ICEP, BUFR_ICEV, BUFR_ICES )
          varname='GL'
        case ( bufr_vis )
          varname='VIS'
        case ( bufr_logVis )
          varname='LVIS'
        case ( bufr_radarPrecip )
          varname='PR'
        case ( bufr_logRadarPrecip )
          varname='LPR'
        case ( bufr_gust )
          varname='WGE'
        case ( bufr_riverFlow )
          varname='QO1'
        case ( BUFR_NEFS, bufr_radvel)
          varname='UV'
        case default
          !
          ! Search for constituents. Identification depends on value and presence of second parameter.
          !
          if (present(varNumberChm_opt)) then
            if (trim(modelName) == 'GEM-MACH') then
              select case (varNumberChm_opt)
                case(BUFR_NECH_O3)
                  varname='TO3'
                case(BUFR_NECH_H2O)
                  varname='HU'
                case(BUFR_NECH_CH4)
                  varname='TCH4'
                case(BUFR_NECH_CO2)
                  varname='TCO2'
                case(BUFR_NECH_CO)
                  varname='TCO'
                case(BUFR_NECH_NO2)
                  varname='TNO2'
                case(BUFR_NECH_N2O)
                  varname='TN2O'
                case(BUFR_NECH_NO)
                  varname='TNO'
                case(BUFR_NECH_HCHO)
                  varname='THCH'
                case(BUFR_NECH_SO2)
                  varname='TSO2'
                case(BUFR_NECH_NH3)
                  varname='TNH3'
                case(BUFR_NECH_PM25)
                  varname='AF'
                case(BUFR_NECH_PM10)
                  varname='AC'
                case default
                  call utl_abort('vnl_varnameFromVarnum: Unknown variable number! ' // &
                    utl_str(varNumber) // ', ' // utl_str(varNumberChm_opt))
              end select
            else if (trim(modelName) == 'GEM') then
              select case (varNumberChm_opt)
                case(BUFR_NECH_O3)
                  varname='O3L'
                case(BUFR_NECH_CH4)
                  varname='CH4L'
                case(BUFR_NECH_N2O)
                  varname='N2OL'
                case default
                  call utl_abort('vnl_varnameFromVarnum: Unknown variable number! ' // &
                    utl_str(varNumber) // ', ' // utl_str(varNumberChm_opt))
              end select
            else
              call utl_abort('vnl_varnameFromVarnum: Unknown model! ' // trim(modelName))
            end if
          else
            call utl_abort('vnl_varnameFromVarnum: Unknown variable number! ' // utl_str(varNumber))
          end if
      end select

    end function vnl_varnameFromVarnum

   !--------------------------------------------------------------------------
   ! vnl_varnumFromVarName
   !--------------------------------------------------------------------------
    function vnl_varnumFromVarName(varName,varKind_opt) result(varNumber)
      !
      ! :Purpose: Identifies varNumber from varName for use in assimilating
      !           obs in the CH family.   
      !           Here, for weather variables, there is a 1-1 association between
      !           a variable name and an observation unit.
      !           So one must provide the name directly associated to a single
      !           BUFR code.
      !           As such, weather variable varNames may not necessarily be a
      !           member of the vnl_varNameList for this routine only.
      !
      !           For constituents, the varNumber refers only to the field/
      !           variable and not units. As consequence, there is a unique
      !           pairing of varNumbers with the varNames from vnl_VarNameList.
      !
      implicit none

      ! Arguments:
      character(len=*),            intent(in) :: varName
      character(len=*), optional,  intent(in) :: varKind_opt
      ! Result:
      integer    :: varNumber

      varNumber=0
      select case (varName)

      ! Weather variables. Must provide name directly associated to a single
      ! BUFR code. As such, the varName may not necessarily be a member of the
      ! vnl_varNameList for this routine only.

      case('UU')
        varNumber=BUFR_NEUU
      case('US')
        varNumber=BUFR_NEUS
      case('VV')
        varNumber=BUFR_NEVV
      case('VS')
        varNumber=BUFR_NEVS
      case('TT')
        varNumber=BUFR_NETT
      case('TS')
        varNumber=BUFR_NETS
      case('RF')
        varNumber=BUFR_NERF
      case('BD')
        varNumber=BUFR_NEBD
      case('ZD')
        varNumber=BUFR_NEZD
      case('GZ')
        varNumber=BUFR_NEGZ
      case('DZ')
        varNumber=BUFR_NEDZ
      case('UV')
        varNumber=BUFR_NEFS 
      case('HU')
        if (present(varKind_opt)) then
           if (varKind_opt == 'CH') then
              varNumber=BUFR_NECH_H2O
           else
              varNumber=BUFR_NEHU
           end if
         else
            varNumber=BUFR_NEHU
         end if
      case('HS')
        varNumber=BUFR_NEHS
      case('SS')
        varNumber=BUFR_NESS
      case('P0','PS')
        varNumber=BUFR_NEPS
      case('PN')
        varNumber=BUFR_NEPN
      case('DW')
        varNumber=BUFR_NEDW
      case('WGE')
        varNumber=bufr_gust
      case('LVIS')
        varNumber=bufr_logVis
      case('VIS')
        varNumber=bufr_vis

      ! Atmospheric constituents other than HU
      case('TO3','O3L')
        varNumber=BUFR_NECH_O3
      case('TH2O')
        varNumber=BUFR_NECH_H2O
      case('TCH4','CH4L')
        varNumber=BUFR_NECH_CH4
      case('TCO2')
        varNumber=BUFR_NECH_CO2
      case('TCO','ECO')
        varNumber=BUFR_NECH_CO
      case('TNO2','ENO2')
        varNumber=BUFR_NECH_NO2
      case('TN2O','N2OL')
        varNumber=BUFR_NECH_N2O
      case('TNO')
        varNumber=BUFR_NECH_NO
      case('THCH','EHCH')
        varNumber=BUFR_NECH_HCHO
      case('TSO2','ESO2')
        varNumber=BUFR_NECH_SO2
      case('TNH3','ENH3')
        varNumber=BUFR_NECH_NH3
      case('AF')
        varNumber=BUFR_NECH_PM25
      case('AC')
        varNumber=BUFR_NECH_PM10

      case default
         call utl_abort('vnl_varnumFromVarName: Unknown variable name ' // trim(varName) )
      end select

    end function vnl_varnumFromVarname

   !--------------------------------------------------------------------------
   ! vnl_varLevelFromVarname
   !--------------------------------------------------------------------------
    function vnl_varLevelFromVarname(varName) result(varLevel)
      !
      ! :Purpose: To get variable level list from variable name 
      !
      implicit none

      ! Arguments:
      character(len=*), intent(in)   :: varName
      ! Result:
      character(len=4)               :: varLevel

      ! Locals:
      integer                :: ierr
      logical, save          :: firstTime = .true.

      ! Namelist variables
      character(len=4), save :: forceSfcOnly(vnl_numVarMax) ! List of 3D variable names only allocated at the surface

      NAMELIST /namvnl/forceSfcOnly

      if (firstTime) then
        firstTime = .false.
        ! default values (not a valid variable name)
        forceSfcOnly(:) = 'XXXX'

        if (utl_isNamelistPresent('namvnl','./flnml')) then
          call utl_tmg_start(181,'low-level--readNML')
          read (utl_flnml, nml = NAMVNL, iostat = ierr)
          if ( ierr /= 0 ) call utl_abort('vnl_varLevelFromVarname: Error reading namelist')
          if ( mmpi_myid == 0 ) write(*,nml=namvnl)
          call utl_tmg_stop(181)
        else
          write(*,*)
          write(*,*) 'vnl_varLevelFromVarname: namvnl is missing in the namelist. The default value will be taken.'
        end if
      end if

      varLevel = varLevelList(vnl_varListIndex(varName))
      if (any(forceSfcOnly(:) == varName)) then
        if (varLevel == 'TH') then
          varLevel = 'SFTH'
        else if (varLevel == 'MM') then
          varLevel = 'SFMM'
        else
          call utl_abort('vnl_varLevelFromVarname: something is wrong')
        end if
      end if

    end function vnl_varLevelFromVarname

   !--------------------------------------------------------------------------
   ! vnl_varLevelFromVarnum
   !--------------------------------------------------------------------------
    function vnl_varLevelFromVarnum(varNumber) result(varLevel)
      !
      ! :Purpose: To get variable level list from the variable number 
      !
      implicit none

      ! Arguments:
      integer,                    intent(in) :: varNumber
      ! Result:
      character(len=4)              :: varLevel

      ! Locals:
      character(len=4)              :: varName

      if (bufr_isAtmosConstituent(varNumber)) then
        varLevel = 'TH'
      else
        varName = vnl_varnameFromVarnum(varNumber)
        varLevel = varLevelList(vnl_varListIndex(varName))
      end if

    end function vnl_varLevelFromVarnum

   !--------------------------------------------------------------------------
   ! vnl_varKindFromVarname
   !--------------------------------------------------------------------------
    function vnl_varKindFromVarname(varName) result(varKind)
      !
      ! :Purpose: To get variable kind list from the variable number 
      !
      implicit none

      ! Arguments:
      character(len=*), intent(in) :: varName
      ! Result:
      character(len=2) :: varKind
      
      varKind = varKindList(vnl_varListIndex(varName))

    end function vnl_varKindFromVarname

   !--------------------------------------------------------------------------
   ! vnl_varNamesFromExistList
   !--------------------------------------------------------------------------
    subroutine vnl_varNamesFromExistList(varNames, varExistList)
      !
      ! :Purpose: To get variable names from the variable existList 
      !

      implicit none

      ! Arguments:
      logical,                   intent(in)  :: varExistList(:) ! a logical switch for the current variable 
      character(len=4), pointer, intent(out) :: varNames(:)     ! variable names

      ! Local:
      integer :: varIndex, numFound

      if (associated(varNames)) then
        call utl_abort('vnl_varNamesFromExistList: varNames must be NULL pointer on input')
      end if

      numFound = 0
      do varIndex = 1, vnl_numvarmax
        if (varExistList(varIndex)) numFound = numFound + 1
      end do
      allocate(varNames(numFound))

      numFound = 0
      do varIndex = 1, vnl_numvarmax
        if (varExistList(varIndex)) then
          numFound = numFound + 1
          varNames(numFound) = vnl_varNameList(varIndex)
        end if
      end do

    end subroutine vnl_varNamesFromExistList
 
   !--------------------------------------------------------------------------
   ! vnl_varMassFromVarNum
   !--------------------------------------------------------------------------
    function vnl_varMassFromVarNum(varNumber) result(varMass)
      !
      ! :Purpose: Identifies constituent molar mass from varNum for use in conversions for the CH family.   
      !
      implicit none

      ! Arguments:
      integer, intent(in) :: varNumber
      ! Result:
      real(8)             :: varMass

      if ( varNumber == BUFR_NECH_O3 ) then
        varMass = MPC_MOLAR_MASS_O3_R8
      else if ( varNumber == BUFR_NECH_H2O ) then
        varMass = MPC_MOLAR_MASS_VAPOUR_R8
      else if ( varNumber == BUFR_NECH_CH4 ) then
        varMass = MPC_MOLAR_MASS_CH4_R8
      else if ( varNumber == BUFR_NECH_CO2 ) then
        varMass = MPC_MOLAR_MASS_CO2_R8
      else if ( varNumber == BUFR_NECH_CO ) then
        varMass = MPC_MOLAR_MASS_CO_R8
      else if ( varNumber == BUFR_NECH_NO2 ) then
        varMass = MPC_MOLAR_MASS_NO2_R8
      else if ( varNumber == BUFR_NECH_N2O ) then
        varMass = MPC_MOLAR_MASS_N2O_R8
      else if ( varNumber == BUFR_NECH_HCHO ) then
        varMass = MPC_MOLAR_MASS_HCHO_R8
      else if ( varNumber == BUFR_NECH_SO2 ) then
        varMass = MPC_MOLAR_MASS_SO2_R8
      else if ( varNumber == BUFR_NECH_NH3 ) then
        varMass = MPC_MOLAR_MASS_NH3_R8
      else if ( varNumber == BUFR_NECH_NO ) then
        varMass = MPC_MOLAR_MASS_NO_R8
      else if ( varNumber == BUFR_NECH_PM25 ) then
        varMass = 1.0d0 ! no scaling
      else if ( varNumber == BUFR_NECH_PM10 ) then
        varMass = 1.0d0 ! no scaling
      else
        call utl_abort('vnl_varMassFromVarNum: Constituent id number ' // &
                       utl_str(varNumber) // ' not recognized' )
      end if
      
    end function vnl_varMassFromVarNum

   !--------------------------------------------------------------------------
   ! vnl_varMassFromVarName
   !--------------------------------------------------------------------------
    function vnl_varMassFromVarName(varName) result(varMass)
      !
      ! :Purpose: Identifies constituent molar mass from varName for use in conversions for the CH family.   
      !
      implicit none

      ! Arguments:
      character(len=*),  intent(in) :: varName
      ! Result:
      real(8)                       :: varMass

      if ( varName == 'TO3' .or. varName == 'O3L'  ) then
        varMass = MPC_MOLAR_MASS_O3_R8
      else if ( varName == 'LQ' .or.  varName == 'HU' ) then
        varMass = MPC_MOLAR_MASS_VAPOUR_R8
      else if ( varName == 'TCH4' .or. varName == 'CH4L'   ) then
        varMass = MPC_MOLAR_MASS_CH4_R8
      else if ( varName == 'TCO2' ) then
        varMass = MPC_MOLAR_MASS_CO2_R8
      else if ( varName == 'TCO' ) then
        varMass = MPC_MOLAR_MASS_CO_R8
      else if ( varName == 'TNO2' ) then
        varMass = MPC_MOLAR_MASS_NO2_R8
      else if ( varName == 'TN2O' .or. varName == 'N2OL'   ) then
        varMass = MPC_MOLAR_MASS_N2O_R8
      else if ( varName == 'THCH' ) then
        varMass = MPC_MOLAR_MASS_HCHO_R8
      else if ( varName == 'TSO2' ) then
        varMass = MPC_MOLAR_MASS_SO2_R8
      else if ( varName == 'TNH3' ) then
        varMass = MPC_MOLAR_MASS_NH3_R8
      else if ( varName == 'TNO' ) then
        varMass = MPC_MOLAR_MASS_NO_R8
      else if ( varName == 'AF' ) then
        varMass = 1.0d0 ! no scaling
      else if ( varName == 'AC' ) then
        varMass = 1.0d0 ! no scaling
      else
        call utl_abort('vnl_varMassFromVarName: Molar mass not found for varName ' // &
                       trim(varName) )
      end if
      
    end function vnl_varMassFromVarName

    !--------------------------------------------------------------------------
    ! vnl_isPhysicsVar
    !--------------------------------------------------------------------------
    function vnl_isPhysicsVar(varName) result(isPhysicsVar)
      !
      ! :Purpose: Signals if variable is expected to be on the "physics" grid.
      !
      implicit none

      ! Arguments:
      character(len=*),  intent(in) :: varName
      ! Result:
      logical                       :: isPhysicsVar

      select case (trim(varName))
      case ( 'I0','I1','I2','I3','I4','I5','I6','I7','I8','I9', &
             'DN','FB','FI','PR','LPR' )
        isPhysicsVar = .true.
      case default
        isPhysicsVar = .false.
      end select

    end function vnl_isPhysicsVar

    !-----------------------------------------------------------------------
    ! vnl_isCloudVar
    !----------------------------------------------------------------------
    function vnl_isCloudVar(varName) result(isCloud)
      !
      ! :Purpose: determine if varName is cloud variable.
      !
      implicit none
  
      ! Arguments:
      character(len=*), intent(in) :: varName
      ! Result:
      logical                      :: isCloud
  
      ! Locals:
      integer :: varNameIndex
  
      isCloud = .false.
      do varNameIndex = 1, vnl_numvarmaxCloud
        if (trim(varName) == trim(vnl_varNameListCloud(varNameIndex))) then
          isCloud = .true.
          return
        end if
      end do
  
    end function vnl_isCloudVar

    !--------------------------------------------------------------------------
    ! vnl_addToVarNames
    !--------------------------------------------------------------------------
    subroutine vnl_addToVarNames(varNamesInOut, varNameToAdd, imposeVnlOrder_opt)
      !
      ! :Purpose: Add an additional varName to an existing list of varNames.
      !           But only add if it is not already in the list.
      !
      implicit none

      ! Arguments:
      character(len=*), pointer, intent(inout) :: varNamesInOut(:) ! input/output variable names
      character(len=*),          intent(in)    :: varNameToAdd     ! variable name to add to the list of existing variables
      logical,         optional, intent(in)    :: imposeVnlOrder_opt ! choose to sort output list by order in vnl list

      ! Locals:
      integer :: lenVarNames, nameIndex, varIndex
      logical :: alreadyInList, imposeVnlOrder
      character(len=4), pointer :: varNamesOut(:)
      character(len=4), pointer :: varNamesSorted(:)
      character(len=4)          :: varName

      lenVarNames = size(varNamesInOut)

      if (present(imposeVnlOrder_opt)) then
        imposeVnlOrder = imposeVnlOrder_opt
      else
        imposeVnlOrder = .false.
      end if

      ! Check if the name is already in the list
      alreadyInList = .false.
      CheckLoop: do nameIndex = 1, lenVarNames
        if (trim(varNamesInOut(nameIndex)) == trim(varNameToAdd)) then
          alreadyInList = .true.
          exit CheckLoop
        end if
      end do CheckLoop

      if (alreadyInList) then
        ! If already in list, just return leaving InOut list unchanged
        return
      end if

      ! Create output list 1 element longer than input list
      allocate(varNamesOut(lenVarNames+1))
      varNamesOut(1:lenVarNames) = varNamesInOut(:)
      varNamesOut(lenVarNames+1) = trim(varNameToAdd)

      lenVarNames = size(varNamesOut)

      if (imposeVnlOrder) then
        ! Ensure order follows that in varNameList_mod
        allocate(varNamesSorted(lenVarNames))
        nameIndex = 0
        do varIndex = 1, vnl_numvarmax
          varName = vnl_varNameList(varIndex)
          if (ANY(varNamesOut(:) == varName)) then
            nameIndex = nameIndex + 1
            varNamesSorted(nameIndex) = varName
          end if
        end do
        varNamesOut(:) = varNamesSorted(:)
        deallocate(varNamesSorted)
      end if

      ! Replace the input list with output list
      deallocate(varNamesInOut)
      varNamesInOut => varNamesOut

    end subroutine vnl_addToVarNames

    !-----------------------------------------------------------------------
    ! vnl_varNamePresentInFile
    !----------------------------------------------------------------------
    function vnl_varNamePresentInFile(varName, fileName, typvar_opt) result(found)
      !
      !:Purpose: Determine if a given variable name is present within a file.
      !          This function supports both "standard" files and NetCDF files.
      !
      implicit none

      ! Arguments:
      character(len=*), intent(in)           :: varName    ! variable name
      character(len=*), intent(in)           :: fileName   ! file name
      character(len=*), intent(in), optional :: typvar_opt ! typvar used for RPN standard files

      ! Result:
      logical                                :: found      ! true if variable name is found, false, the converse

      ! Locals:
      integer :: fnom, fstouv, fstfrm, fclos, fstinf
      integer :: ni, nj, nk, key, ierr
      integer :: unit, ncid, varID
      character(len=2)   :: typvar
      character(len=10)  :: varNameNetCDF

      unit = 0
      
      if (present(typvar_opt)) then
        typvar = trim(typvar_opt)
      else
        typvar = ' '
      end if

      if (trim(utl_fileType(trim(fileName))) == 'FST') then

        ierr = fnom(unit, fileName, 'RND+OLD+R/O', 0)
        ierr = fstouv(unit, 'RND+OLD')
      
        key = fstinf(unit, ni, nj, nk, -1 ,' ', -1, -1, -1, typvar, trim(varName))
  
        if (key > 0)  then
          found = .true.
        else
          found = .false.
        end if

        ierr =  fstfrm(unit)
        ierr =  fclos (unit)

      else if (trim(utl_fileType(trim(fileName))) == 'NetCDF') then

        varNameNetCDF = vnl_varNameNetCDF(varName, trim(fileName))

        if (varNameNetCDF == 'notFound') then

          found = .false.

        else

          call utl_checkNetCDFstatus(nf90_open(trim(fileName), nf90_nowrite, ncid))
          ierr = nf90_inq_varid(ncid, trim(varNameNetCDF), varID)      
          if (ierr == nf90_noerr) then
            found = .true.
          else
            found = .false.
          end if
          call utl_checkNetCDFstatus(nf90_close(ncid))

        end if

      else

        call utl_abort('vnl_varNamePresentInFile: unknown input file type: '//&
                       trim(utl_fileType(trim(fileName))))
      end if

    end function vnl_varNamePresentInFile

    !-----------------------------------------------------------------------
    ! vnl_varNameNetCDF
    !----------------------------------------------------------------------
    function vnl_varNameNetCDF(varName, fileName) result(varNameNetCDF)
      !
      ! :Purpose: Return the equivalent variable name used for netCDF files
      !           in use by the NEMO ocean model.
      implicit none
  
      ! Arguments:
      character(len=*), intent(in) :: varName  ! input MIDAS variable name
      character(len=*), intent(in) :: fileName ! NEMO trial file   
          
      ! Result:
      character(len=20) :: varNameNetCDF ! variable name used in NEMO netCDF

      select case(trim(varName))
      
      case('SSH')
        
        if (utl_varPresentInNetcdfFile('zos', trim(fileName))) then
          varNameNetCDF = 'zos'
        else if (utl_varPresentInNetcdfFile('sshn', trim(fileName))) then
          varNameNetCDF = 'sshn'
        else        
          call utl_abort('vnl_varNameNetCDF: no equivalent name varName: '//trim(varName)//&
                         'found in '//trim(fileName))
        end if   

      case('TM') 
   
        if (utl_varPresentInNetcdfFile('toce', trim(fileName))) then
          varNameNetCDF = 'toce'    ! NEMO 3D ocean temperature field, trial
        else if (utl_varPresentInNetcdfFile('tos', trim(fileName))) then
          varNameNetCDF = 'tos'     ! NEMO 2D SST field
        else if (utl_varPresentInNetcdfFile('tn', trim(fileName))) then
          varNameNetCDF = 'tn'      ! NEMO 3D ocean temperature field, restart
        else
          call utl_abort('vnl_varNameNetCDF: no equivalent name varName: '//trim(varName)//&
                         'found in '//trim(fileName))
        end if

      case('SALW')

        if (utl_varPresentInNetcdfFile('soce', trim(fileName))) then
          varNameNetCDF = 'soce'
        else if (utl_varPresentInNetcdfFile('sn', trim(fileName))) then
          varNameNetCDF = 'sn'
        else
          call utl_abort('vnl_varNameNetCDF: no equivalent name varName: '//trim(varName)//&
                         'found in '//trim(fileName))
        end if
   
      case('UUW')

        if (utl_varPresentInNetcdfFile('uo', trim(fileName))) then
          varNameNetCDF = 'uo'
        else if (utl_varPresentInNetcdfFile('un', trim(fileName))) then
          varNameNetCDF = 'un'
        else 
          call utl_abort('vnl_varNameNetCDF: no equivalent name varName: '//trim(varName)//&
                         'found in '//trim(fileName))
        end if

      case('VVW')

        if (utl_varPresentInNetcdfFile('vo', trim(fileName))) then
          varNameNetCDF = 'vo'
        else if (utl_varPresentInNetcdfFile('vn', trim(fileName))) then
          varNameNetCDF = 'vn'
        else 
          call utl_abort('vnl_varNameNetCDF: no equivalent name varName: '//trim(varName)//&
                         'found in '//trim(fileName))
        end if

      case('SSS')

        if (utl_varPresentInNetcdfFile('sss', trim(fileName))) then
          varNameNetCDF = 'sss'
        else
          call utl_abort('vnl_varNameNetCDF: no equivalent name varName: '//trim(varName)//&
                         'found in '//trim(fileName))
        end if
         
      case default

        varNameNetCDF = 'notFound'
        write(*,*) 'vnl_varNameNetCDF: WARNING: no equivalent name for NetCDF files for varName: '//trim(varName)

      end select
      
      write(*,*) 'vnl_varNameNetCDF: ', trim(varName),' is ', trim(varNameNetCDF), ' in netCDF file.'

    end function vnl_varNameNetCDF

end module varNameList_mod
