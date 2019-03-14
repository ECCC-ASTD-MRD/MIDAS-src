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
!! MODULE varNameList (prefix="vnl" category='7. Low-level data objects and utilities')
!!
!! *Purpose*: Contains a list of all possible variable names that can be used 
!!            as analysis variables along with additional information for each and
!!            procedures for accessing this information
!!
!--------------------------------------------------------------------------
module varNameList_mod
  use bufr_mod
  use utilities_mod
  implicit none
  save
  private

  ! public variables (parameters)
  public :: vnl_numvarmax3D, vnl_numvarmax2D, vnl_numvarmax
  public :: vnl_varNameList3D, vnl_varNameList2D, vnl_varNameList

  ! public procedures
  public :: vnl_varListIndex3d, vnl_varListIndex2d, vnl_varListIndex, vnl_varnameFromVarnum
  public :: vnl_varLevelFromVarname, vnl_varLevelFromVarnum
  public :: vnl_varKindFromVarname, vnl_varnumFromVarname
  public :: vnl_varNamesFromExistList

  integer, parameter          :: vnl_numvarmax3D = 31, vnl_numvarmax2D = 13

  character(len=4), parameter :: vnl_varNameList3D(vnl_numvarmax3D) = (/                         &
                                 'UU  ','VV  ','GZ  ','TT  ','HU  ','LQ  ','ES  ','VT  ',        &
                                 'PP  ','CC  ','UC  ','UT  ','TB  ','DW  ','QR  ','DD  ',        &
                                 'TO3 ','TCH4','TCO2','TCO ','TNO2','TN2O','THCH','TSO2',        &
                                 'TNH3','AF  ','AC  ','TNO ','ALFA','VIS ','LVIS'/)

  character(len=2), parameter :: varLevelList3D(vnl_numvarmax3D)     = (/                        &
                                 'MM',  'MM',  'TH',  'TH',  'TH',  'TH',  'TH',  'TH',          &
                                 'MM',  'MM',  'MM',  'TH',  'TH',  'TH',  'MM',  'MM',          &
                                 'TH',  'TH',  'TH',  'TH',  'TH',  'TH',  'TH',  'TH',          &
                                 'TH',  'TH',  'TH',  'TH',  'MM',  'TH',  'TH'/)

  character(len=2), parameter :: varKindList3D(vnl_numvarmax3D)     = (/                         &
                                 'MT',  'MT',  'MT',  'MT',  'MT',  'MT',  'MT',  'MT',          &
                                 'MT',  'MT',  'MT',  'MT',  'MT',  'MT',  'MT',  'MT',          &
                                 'CH',  'CH',  'CH',  'CH',  'CH',  'CH',  'CH',  'CH',          &
                                 'CH',  'CH',  'CH',  'CH',  'MT',  'MT',  'MT'/)

  character(len=4), parameter :: vnl_varNameList2D(vnl_numvarmax2D) = (/ &
                                 'P0  ','TG  ','UP  ','PB  ','ECO ', 'ENO2', 'EHCH', 'ESO2', 'ENH3' , &
                                 'GL  ','WGE ','BIN ','MG  '/)

  character(len=2), parameter :: varLevelList2D(vnl_numvarmax2D) = (/    &
                                 'SF',  'SF',  'SF',  'SF', 'SF',  'SF',  'SF',  'SF',  'SF',  &
                                 'SF',  'SF',  'SF',  'SF'/)

  character(len=2), parameter :: varKindList2D(vnl_numvarmax2D) = (/     &
                                 'MT', 'MT', 'MT', 'MT', 'CH', 'CH', 'CH', 'CH', 'CH', &
                                 'MT', 'MT', 'MT', 'MT'/)

  integer, parameter          :: vnl_numvarmax = vnl_numvarmax3D + vnl_numvarmax2D

  character(len=4), parameter :: vnl_varNameList(vnl_numvarmax) = (/ vnl_varNameList3D, vnl_varNameList2D /)
  character(len=2), parameter :: varLevelList   (vnl_numvarmax) = (/ varLevelList3D   , varLevelList2D    /)
  character(len=2), parameter :: varKindList    (vnl_numvarmax) = (/ varKindList3D    , varKindList2D     /)

  contains

    !--------------------------------------------------------------------------
    ! vnl_varListIndex3d
    !--------------------------------------------------------------------------
    function vnl_varListIndex3d(varName) result(listIndex)
      implicit none
      character(len=*), intent(in) :: varName
      integer                      :: jvar,listIndex

      listIndex=-1
      do jvar=1,vnl_numvarmax3D
        if(varName.eq.vnl_varNameList3d(jvar)) then
          listIndex=jvar
          exit
        endif
      enddo

      if(listIndex.le.0) then
        call utl_abort('vnl_varListIndex3D: Unknown variable name! ' // varName)
      endif

    end function vnl_varListIndex3d

    !--------------------------------------------------------------------------
    ! vnl_varListIndex2d
    !--------------------------------------------------------------------------
    function vnl_varListIndex2d(varName) result(listIndex)
      implicit none
      character(len=*), intent(in) :: varName
      integer                      :: jvar,listIndex

      listIndex=-1
      do jvar = 1, vnl_numvarmax2D
        if( varName == vnl_varNameList2d(jvar) ) then 
          listIndex=jvar
          exit
        endif
      enddo

      if(listIndex <= 0) then
        call utl_abort('vnl_varListIndex2D: Unknown variable name! ' // varName)
      endif

    end function vnl_varListIndex2d

    !--------------------------------------------------------------------------
    ! vnl_varListIndex
    !--------------------------------------------------------------------------
    function vnl_varListIndex(varName) result(listIndex)
      implicit none
      character(len=*), intent(in) :: varName
      integer                      :: jvar,listIndex

      listIndex=-1
      do jvar=1,vnl_numvarmax
        if(varName == vnl_varNameList(jvar)) then 
          listIndex=jvar
          exit
        endif
      enddo

      if(listIndex <= 0) then
        call utl_abort('vnl_varListIndex: Unknown variable name! ' // varName)
      endif

    end function vnl_varListIndex

    !--------------------------------------------------------------------------
    ! vnl_varnameFromVarnum
    !--------------------------------------------------------------------------
    function vnl_varnameFromVarnum( varNumber, varNumberChm_opt ) result(varName)
      implicit none
      integer, intent(in) :: varNumber
      integer, intent(in), optional :: varNumberChm_opt
      character(len=4)    :: varName
      integer             :: i
      
      varName='    '
      select case (varNumber)
      case ( BUFR_NEUU, BUFR_NEUS, BUFR_NEAL )
        varName='UU'
      case( BUFR_NEVV, BUFR_NEVS )
        varName='VV'
      case( BUFR_NETT, BUFR_NETS )
        varName='TT'
      case( BUFR_NEDZ, BUFR_NEGZ )
        varName='GZ'
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
      case ( BUFR_ICEC  )
        varname='GL'
      case ( bufr_vis )
        varname='LVIS'
      case ( bufr_gust )
        varname='WGE'
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
                 write(*,*) 'vnl_varnameFromVarnum: Unknown variable number! ',varNumber, varNumberChm_opt
                 call utl_abort('vnl_varnameFromVarnum')
           end select
        else
           write(*,*) 'vnl_varnameFromVarnum: Unknown variable number! ',varNumber
           call utl_abort('vnl_varnameFromVarnum')
        endif 
      end select

    end function vnl_varnameFromVarnum

    !--------------------------------------------------------------------------
    ! vnl_varnumFromVarname
    !-------------------------------------------------------------------------- 
    function vnl_varnumFromVarName(varName,varKind_opt) result(varNumber)
      !   Purpose: Identifies varNumber from varName for use in assimilating obs in the CH family.   
      !  
      !   Here, for weather variables, there is a 1-1 association between a variable name and an observation unit.
      !   So one must provide the name directly associated to a single BUFR code.
      !   As such, weather variable varNames may not necessarily be a member of the vnl_varNameList for this routine only.
      !   
      !   For constituents, the varNumber refers only to the field/variable and not units. As consequence,
      !   there is a unique pairing of varNumbers with the varNames from vnl_VarNameList.
      implicit none
      character(len=*),  intent(in) :: varName
      character(len=*),  intent(in), optional :: varKind_opt
      integer    :: varNumber
      
      varNumber=0
      select case (varName)
      
      ! Weather variables. Must provide name directly associated to a single BUFR code.
      ! As such, the varName may not necessarily be a member of the vnl_varNameList for this routine only.

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
        varNumber=bufr_vis

      ! Atmospheric constituents other than HU
      case('TO3')
        varNumber=BUFR_NECH_O3
      case('TH2O')
        varNumber=BUFR_NECH_H2O
      case('TCH4')
        varNumber=BUFR_NECH_CH4
      case('TCO2')
        varNumber=BUFR_NECH_CO2
      case('TCO','ECO')
        varNumber=BUFR_NECH_CO
      case('TNO2','ENO2')
        varNumber=BUFR_NECH_NO2
      case('TN2O')
        varNumber=BUFR_NECH_N2O
      case('TNO')
        varNumber=BUFR_NECH_NO
      case('HCHO','THCH','EHCH')
        varNumber=BUFR_NECH_HCHO
      case('TSO2','ESO2')
        varNumber=BUFR_NECH_HCHO
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
      implicit none

      character(len=*), intent(in)   :: varName
      character(len=2)               :: varLevel

      varLevel = varLevelList(vnl_varListIndex(varName))

    end function vnl_varLevelFromVarname


    !--------------------------------------------------------------------------
    ! vnl_varLevelFromVarnum
    !--------------------------------------------------------------------------
    function vnl_varLevelFromVarnum(varNumber,varNumberChm_opt) result(varLevel)
      implicit none

      integer, intent(in)           :: varNumber
      integer, intent(in), optional :: varNumberChm_opt
      character(len=2)              :: varLevel
      character(len=4)              :: varName

      varName = vnl_varnameFromVarnum(varNumber,varNumberChm_opt=varNumberChm_opt)
      varLevel = varLevelList(vnl_varListIndex(varName))

    end function vnl_varLevelFromVarnum

    !--------------------------------------------------------------------------
    ! vnl_varKindFromVarname
    !--------------------------------------------------------------------------
    function vnl_varKindFromVarname(varName) result(varKind)
      implicit none

      character(len=*), intent(in) :: varName
      character(len=2) :: varKind
      
      varKind = varKindList(vnl_varListIndex(varName))

    end function vnl_varKindFromVarname

    !--------------------------------------------------------------------------
    ! vnl_varNamesFromExistList
    !--------------------------------------------------------------------------
    subroutine vnl_varNamesFromExistList(varNames,varExistList)
      implicit none

      ! arguments
      logical :: varExistList(:)
      character(len=4), pointer :: varNames(:)

      ! locals
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

end module varNameList_mod
