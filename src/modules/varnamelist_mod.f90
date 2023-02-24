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

module varNameList_mod
  ! MODULE varNameList (prefix='vnl' category='7. Low-level data objects')
  !
  ! :Purpose: Contains a list of all possible variable names that can be used
  !           as analysis variables along with additional information for each
  !           and procedures for accessing this information
  !
  use bufr_mod
  use midasMpi_mod
  use utilities_mod
  use MathPhysConstants_mod

  implicit none
  save
  private

  ! public variables (parameters)
  public :: vnl_numvarmax3D, vnl_numvarmax2D, vnl_numvarmaxOther, vnl_numvarmax
  public :: vnl_varNameList3D, vnl_varNameList2D, vnl_varNameListOther, vnl_varNameList
  public :: vnl_numvarmaxCloud, vnl_varNameListCloud

  ! public procedures
  public :: vnl_varListIndex3d, vnl_varListIndex2d, vnl_varListIndexOther
  public :: vnl_varListIndex, vnl_varnameFromVarnum, vnl_varnameIsValid
  public :: vnl_varLevelFromVarname, vnl_varLevelFromVarnum
  public :: vnl_varKindFromVarname, vnl_varnumFromVarname
  public :: vnl_varNamesFromExistList, vnl_varMassFromVarNum, vnl_varMassFromVarName
  public :: vnl_isPhysicsVar, vnl_isCloudVar
  public :: vnl_addToVarNames

  ! These private parameters permit side-stepping a conflict with the Sphinx documenter,
  ! and an infinite loop
  integer, parameter          :: VNLnumvarmax3D    = 52
  integer, parameter          :: VNLnumvarmax2D    = 37
  integer, parameter          :: VNLnumvarmaxOther =  6
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
                                 'P0LS'/)

  character(len=4), parameter :: varLevelList2D(vnl_numvarmax2D) = (/    &
                                 'SF',  'SF',  'SF',  'SF',  'SF',  'SF',  'SF',  'SF',  'SF',  &
                                 'SS',  'SF',  'SF',  'SF',  'SS',  'SF',  'SF',  'SF',  'SF',  &
                                 'SF',  'SF',  'SF',  'SF',  'SF',  'SF',  'SF',  'SF',  'SF',  &
                                 'SF',  'SF',  'SF',  'SF',  'SF',  'SF',  'SS',  'SS',  'SS',  &
                                 'SF'/)

  character(len=2), parameter :: varKindList2D(vnl_numvarmax2D) = (/     &
                                 'MT',  'MT',  'MT',  'MT',  'CH',  'CH',  'CH',  'CH',  'CH', &
                                 'OC',  'MT',  'MT',  'MT',  'OC',  'HY',  'HY',  'HY',  'HY', &
                                 'MT',  'MT',  'MT',  'MT',  'MT',  'MT',  'MT',  'MT',  'MT', &
                                 'MT',  'MT',  'MT',  'MT',  'HY',  'MT',  'OC',  'OC',  'OC', &
                                 'MT'/)

  character(len=4), parameter :: vnl_varNameListOther(vnl_numvarmaxOther) = (/ &
                                 'I0  ','I1  ','I7  ','I9  ','SD  ','AL  '/)

  character(len=4), parameter :: varLevelListOther(vnl_numvarmaxOther) = (/    &
                                 'OT',  'OT',  'OT',  'OT',  'OT',  'OT'  /)

  character(len=2), parameter :: varKindListOther(vnl_numvarmaxOther) = (/     &
                                 'LD',  'LD',  'LD',  'LD',  'LD',  'LD'  /) ! LD = Land

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
      integer                      :: listIndex
      
      !Local:
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
      integer                      :: listIndex

      !Local:
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
      integer                      :: listIndex

      !Local:
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

      !Arguments:
      character(len=*), intent(in) :: varName
      integer                      :: listIndex

      !Local:
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
      logical                      :: isValid
      
      ! Local:
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
    function vnl_varnameFromVarnum( varNumber, varNumberChm_opt, modelName_opt ) result(varName)
      !
      ! :Purpose: To get the variable name from the variable number
      !

      implicit none

      !Arguments:
      integer, intent(in) :: varNumber
      integer, intent(in), optional :: varNumberChm_opt
      character(len=*), intent(in), optional :: modelName_opt
      character(len=4)    :: varName

      varName='    '
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
                 call utl_abort( 'vnl_varnameFromVarnum: Unknown variable number! ' // &
                   utl_str(varNumber) // ', ' // utl_str(varNumberChm_opt) )
           end select
           if (present(modelName_opt)) then
             if (trim(modelName_opt) == 'GEM') then
               select case (varNumberChm_opt)
                 case(BUFR_NECH_O3)
                   varname='O3L'  
                 case(BUFR_NECH_CH4)
                   varname='CH4L'
                 case(BUFR_NECH_N2O)
                   varname='N2OL'
                 case default
                   call utl_abort( 'vnl_varnameFromVarnum: Unknown variable number or model! ' // &
                      utl_str(varNumber) // ', ' // utl_str(varNumberChm_opt) // ', ' // trim(modelName_opt) )
               end select
             end if
           end if      
        else
           write(*,*) 'vnl_varnameFromVarnum: Unknown variable number! ',varNumber
           call utl_abort('vnl_varnameFromVarnum')
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

      !Arguments:
      character(len=*),  intent(in) :: varName
      character(len=*),  intent(in), optional :: varKind_opt
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

      !Arguments:
      character(len=*), intent(in)   :: varName
      character(len=4)               :: varLevel

      !Locals:
      integer                :: nulnam, ierr
      integer, external      :: fnom, fclos
      logical, save          :: firstTime = .true.

      ! Namelist variables
      character(len=4), save :: forceSfcOnly(vnl_numVarMax) ! List of 3D variable names only allocated at the surface

      NAMELIST /namvnl/forceSfcOnly

      if (firstTime) then
        firstTime = .false.
        ! default values (not a valid variable name)
        forceSfcOnly(:) = 'XXXX'

        if (utl_isNamelistPresent('namvnl','./flnml')) then
          nulnam = 0
          ierr = fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
          read (nulnam, nml = NAMVNL, iostat = ierr)
          if ( ierr /= 0 ) call utl_abort('vnl_varLevelFromVarname: Error reading namelist')
          if ( mmpi_myid == 0 ) write(*,nml=namvnl)
          ierr = fclos(nulnam)
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
    function vnl_varLevelFromVarnum(varNumber,varNumberChm_opt,modelName_opt) result(varLevel)
      !
      ! :Purpose: To get variable level list from the variable number 
      !
      implicit none

      !Arguments:
      integer, intent(in)           :: varNumber
      integer, intent(in), optional :: varNumberChm_opt
      character(len=*), intent(in), optional :: modelName_opt
      character(len=4)              :: varLevel

      !Local:
      character(len=4)              :: varName

      varName = vnl_varnameFromVarnum(varNumber,varNumberChm_opt=varNumberChm_opt,modelName_opt=modelName_opt)
      varLevel = varLevelList(vnl_varListIndex(varName))

    end function vnl_varLevelFromVarnum

   !--------------------------------------------------------------------------
   ! vnl_varKindFromVarname
   !--------------------------------------------------------------------------
    function vnl_varKindFromVarname(varName) result(varKind)
      !
      ! :Purpose: To get variable kind list from the variable number 
      !

      implicit none

      !Arguments:
      character(len=*), intent(in) :: varName
      character(len=2) :: varKind
      
      varKind = varKindList(vnl_varListIndex(varName))

    end function vnl_varKindFromVarname

   !--------------------------------------------------------------------------
   ! vnl_varNamesFromExistList
   !--------------------------------------------------------------------------
    subroutine vnl_varNamesFromExistList(varNames,varExistList)
      !
      ! :Purpose: To get variable names from the variable existList 
      !

      implicit none

      !Arguments:
      logical :: varExistList(:)
      character(len=4), pointer :: varNames(:)

      !Local:
      integer :: varIndex, numFound

      if (associated(varNames)) then
        call utl_abort('vnl_varNamesFromExistList: varNames must be NULL pointer on input')
      end if

      numFound = 0
      do varIndex = 1, vnl_numvarmax
        if ( varExistList(varIndex) ) numFound = numFound + 1
      end do
      allocate(varNames(numFound))

      numFound = 0
      do varIndex = 1, vnl_numvarmax
        if ( varExistList(varIndex) ) then
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

      !Arguments:
      integer, intent(in) :: varNumber
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

      !Arguments:
      character(len=*),  intent(in) :: varName
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
    function vnl_addToVarNames(varNamesIn,varNameToAdd) result(varNamesOut)
      !
      ! :Purpose: Add an additional varName to an existing list of varNames
      !
      implicit none

      ! Arguments:
      character(len=*),  intent(in) :: varNamesIn(:)
      character(len=*),  intent(in) :: varNameToAdd
      character(len=4), pointer     :: varNamesOut(:)

      ! Locals:
      integer :: lenVarNames

      lenVarNames = size(varNamesIn)
      allocate(varNamesOut(lenVarNames+1))

      varNamesOut(1:lenVarNames) = varNamesIn(:)
      varNamesOut(lenVarNames+1) = varNameToAdd

    end function vnl_addToVarNames

end module varNameList_mod
