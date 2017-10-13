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
!! MODULE varNameList (prefix="vnl")
!!
!! *Purpose*: Contains a list of all possible variable names that can be used 
!!            as analysis variables along with additional information for each and
!!            procedures for accessing this information
!!
!--------------------------------------------------------------------------
module varNameList_mod
!
! Revisions: 
!            M. Sitwell (ARQI/AQRD) Oct 2015
!            - Added varKindList and vnl_varKindFromVarname to identify the kind of field
!              where 'MT'=meteorological and 'CH'=chemical constituent
!
  use bufr
  use utilities_mod
  implicit none
  save
  private

  ! public variables (parameters)
  public :: vnl_numvarmax3D, vnl_numvarmax2D, vnl_numvarmax
  public :: vnl_varNameList3D, vnl_varNameList2D, vnl_varNameList

  ! public procedures
  public :: vnl_varListIndex3d, vnl_varListIndex2d, vnl_varListIndex, vnl_varnameFromVarnum
  public :: vnl_varLevelFromVarname, vnl_varLevelFromVarnum, vnl_varTypeFromVarname
  public :: vnl_varKindFromVarname, vnl_varnumFromVarname

  integer, parameter          :: vnl_numvarmax3D = 29, vnl_numvarmax2D = 10

  character(len=4), parameter :: vnl_varNameList3D(vnl_numvarmax3D) = (/                         &
                                 'UU  ','VV  ','GZ  ','TT  ','HU  ','LQ  ','ES  ','VT  ',        &
                                 'PP  ','CC  ','UC  ','UT  ','TB  ','DW  ','QR  ','DD  ',        &
                                 'TO3 ','TCH4','TCO2','TCO ','TNO2','TN2O','THCH','TSO2',        &
                                 'TNH3','AF  ','AC  ','TNO ','ALFA'/)

  character(len=2), parameter :: varLevelList3D(vnl_numvarmax3D)     = (/                        &
                                 'MM',  'MM',  'TH',  'TH',  'TH',  'TH',  'TH',  'TH',          &
                                 'MM',  'MM',  'MM',  'TH',  'TH',  'TH',  'MM',  'MM',          &
                                 'TH',  'TH',  'TH',  'TH',  'TH',  'TH',  'TH',  'TH',          &
                                 'TH',  'TH',  'TH',  'TH',  'MM'/)

  character(len=5), parameter :: varTypeList3D(vnl_numvarmax3D)     = (/                                  &
                                 'MODEL', 'MODEL', 'MODEL', 'MODEL', 'MODEL', 'DIAG ', 'DIAG ', 'DIAG ',  &
                                 'DIAG ', 'DIAG ', 'DIAG ', 'DIAG ', 'DIAG ', 'DIAG ', 'DIAG ', 'DIAG ',  &
                                 'MODEL', 'MODEL', 'MODEL', 'MODEL', 'MODEL', 'MODEL', 'MODEL', 'MODEL',  &
                                 'MODEL', 'MODEL', 'MODEL', 'MODEL', 'OTHER'/)

  character(len=2), parameter :: varKindList3D(vnl_numvarmax3D)     = (/                         &
                                 'MT',  'MT',  'MT',  'MT',  'MT',  'MT',  'MT',  'MT',          &
                                 'MT',  'MT',  'MT',  'MT',  'MT',  'MT',  'MT',  'MT',          &
                                 'CH',  'CH',  'CH',  'CH',  'CH',  'CH',  'CH',  'CH',          &
                                 'CH',  'CH',  'CH',  'CH',  'MT'/)

  character(len=4), parameter :: vnl_varNameList2D(vnl_numvarmax2D) = (/ &
                                 'P0  ','TG  ','UP  ','PB  ','ECO ', 'ENO2', 'EHCH', 'ESO2', 'ENH3' , 'GL  '/)

  character(len=2), parameter :: varLevelList2D(vnl_numvarmax2D) = (/    &
                                 'SF',  'SF',  'SF',  'SF', 'SF',  'SF',  'SF',  'SF',  'SF',  'SF'/)

  character(len=5), parameter :: varTypeList2D(vnl_numvarmax2D) = (/     &
                                 'MODEL', 'MODEL', 'DIAG ', 'DIAG ', 'MODEL', 'MODEL', 'MODEL', 'MODEL', 'MODEL', 'MODEL'/)

  character(len=2), parameter :: varKindList2D(vnl_numvarmax2D) = (/     &
                                 'MT', 'MT', 'MT', 'MT', 'CH', 'CH', 'CH', 'CH', 'CH', 'MT'/)

  integer, parameter          :: vnl_numvarmax = vnl_numvarmax3D + vnl_numvarmax2D

  character(len=4), parameter :: vnl_varNameList(vnl_numvarmax) = (/ vnl_varNameList3D, vnl_varNameList2D /)
  character(len=2), parameter :: varLevelList   (vnl_numvarmax) = (/ varLevelList3D   , varLevelList2D    /)
  character(len=5), parameter :: varTypeList    (vnl_numvarmax) = (/ varTypeList3D    , varTypeList2D     /)
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
      do jvar=1,vnl_numvarmax2D
        if(varName.eq.vnl_varNameList2d(jvar)) then 
          listIndex=jvar
          exit
        endif
      enddo

      if(listIndex.le.0) then
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
        if(varName.eq.vnl_varNameList(jvar)) then 
          listIndex=jvar
          exit
        endif
      enddo

      if(listIndex.le.0) then
        call utl_abort('vnl_varListIndex: Unknown variable name! ' // varName)
      endif

    end function vnl_varListIndex


    !--------------------------------------------------------------------------
    ! vnl_varnameFromVarnum
    !
    !   Revisions:
    !             Y.J. Rochon (ARQI), Jan. 2015
    !             - Modified search for chemical constituents 
    !          
    !--------------------------------------------------------------------------
    function vnl_varnameFromVarnum(varNumber,varNumber_chm) result(varName)
      implicit none
      integer, intent(in) :: varNumber
      integer, intent(in), optional :: varNumber_chm
      character(len=4)    :: varName
      integer             :: i
      
      varName='    '
      select case (varNumber)
      case(BUFR_NEUU,BUFR_NEUS)
        varName='UU'
      case(BUFR_NEVV,BUFR_NEVS)
        varName='VV'
      case(BUFR_NETT,BUFR_NETS)
        varName='TT'
      case(BUFR_NEDZ,BUFR_NEGZ)
        varName='GZ'
      case(BUFR_NEHU,BUFR_NEHS,BUFR_NEES,BUFR_NESS)
        varName='HU'
      case(BUFR_NEPS,BUFR_NEPN)
        varName='P0'
      case(BUFR_NERF,BUFR_NEBD,BUFR_NEZD)
        varName='TT'   ! temporarily associate refractivity and ZTD with temperature
      case(BUFR_NEDW)
        varName='DW'
      case default
        !
        ! Search for constituents. Identification depends on value and presence of second parameter.
        !
        if (present(varNumber_chm)) then 
           select case (varnumber_chm)
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
                 write(*,*) 'vnl_varnameFromVarnum: Unknown variable number! ',varNumber, varNumber_chm
                 call utl_abort('aborting in vnl_varnameFromVarnum')
           end select
        else
           write(*,*) 'vnl_varnameFromVarnum: Unknown variable number! ',varNumber
           call utl_abort('aborting in vnl_varnameFromVarnum')
        endif 
      end select

    end function vnl_varnameFromVarnum

    !--------------------------------------------------------------------------
    ! vnl_varnumFromVarname
    !
    !   Author: Y.J. Rochon (ARQI), Jan. 2016
    !
    !   Revisions:
    !          
    !   Purpose: Identifies varNumber from varName for use in assimilating obs in the CH family.   
    !  
    !   Here, for weather variables, there is a 1-1 association between a variable name and an observation unit.
    !   So one must provide the name directly associated to a single BUFR code.
    !   As such, weather variable varNames may not necessarily be a member of the vnl_varNameList for this routine only.
    !   
    !   For constituents, the varNumber refers only to the field/variable and not units. As consequence,
    !   there is a unique pairing of varNumbers with the varNames from vnl_VarNameList.
    !   
    !-------------------------------------------------------------------------- 
    function vnl_varnumFromVarName(varName,kind) result(varNumber)

      implicit none
      character(len=*),  intent(in) :: varName
      character(len=*),  intent(in), optional :: kind
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
        if (present(kind)) then
           if (kind.eq.'CH') then
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
    function vnl_varLevelFromVarnum(varNumber,varNumber_chm) result(varLevel)
      implicit none

      integer, intent(in)           :: varNumber
      integer, intent(in), optional :: varNumber_chm
      character(len=2)              :: varLevel
      character(len=4)              :: varName

      if (present(varNumber_chm)) then      
         varName = vnl_varnameFromVarnum(varNumber,varNumber_chm)
      else 
         varName = vnl_varnameFromVarnum(varNumber)
      endif 
      varLevel = varLevelList(vnl_varListIndex(varName))

    end function vnl_varLevelFromVarnum

    !--------------------------------------------------------------------------
    ! vnl_varTypeFromVarname
    !--------------------------------------------------------------------------
    function vnl_varTypeFromVarname(varName) result(varType)
      implicit none

      character(len=*), intent(in)   :: varName
      character(len=5)               :: varType

      varType = varTypeList(vnl_varListIndex(varName))

    end function vnl_varTypeFromVarname


    function vnl_varKindFromVarname(varName) result(varKind)
!
!   Purpose: Get the variable kind from the variable name
!            as set in the array vnl_varKindList
!
!   Author: M. Sitwell (ARQI/AQRD) Oct 2015
!           Following recommendation by M. Buehner.
!
!---------------------------------------------------------------

      implicit none

      character(len=*), intent(in) :: varName
      character(len=2) :: varKind
      
      varKind = varKindList(vnl_varListIndex(varName))

    end function vnl_varKindFromVarname

end module varNameList_mod
