
module bufr_mod
  ! MODULE bufr_mod (prefix='bufr' category='8. Low-level utilities and constants')
  !
  !:Purpose:  This module is a holder for information that is needed from *the
  !           Binary Universal Form for the Representation of meteorological
  !           data* (BUFR), maintained by the World Meteorological Organization.
  !
  !
  !:Routines:
  !
  !       - bufr_IsAtmosConstituent: determines if observation is
  !         constituent/chemistry data
  !
  !       - bufr_IsIntegral: determines if vertical integral measurement
  !
  !       - bufr_isWindComponent: determines if a wind component variable
  !
  !       - bufr_isOceanObs: determines if it is an ocean observation
  !
  public

 ! Table B:
 !          Universal Field-Identity Numbers
 !          (At the CMC, these have also been known as burp_id.)
 !
 ! Use the program codbuf to find out the meaning of these IDs.
 !
 integer, parameter :: bufr_suWindSpeed = 08194 ! source and units of wind speed IW, code 1855, code table
 integer, parameter :: BUFR_NETT        = 12001
 integer, parameter :: BUFR_NEUU        = 11003
 integer, parameter :: BUFR_NEVV        = 11004
 integer, parameter :: BUFR_NEGZ        = 10194
 integer, parameter :: BUFR_NEES        = 12192
 integer, parameter :: BUFR_NEFF        = 11002
 integer, parameter :: BUFR_NEDD        = 11001
 integer, parameter :: BUFR_NEUS        = 11215
 integer, parameter :: BUFR_NEVS        = 11216
 integer, parameter :: BUFR_NETS        = 12004
 integer, parameter :: BUFR_NESS        = 12203
 integer, parameter :: BUFR_NEFS        = 11012
 integer, parameter :: BUFR_NEDS        = 11011
 integer, parameter :: BUFR_NEDZ        = 10192
 integer, parameter :: BUFR_NEPP        = 07004
 integer, parameter :: BUFR_NEPS        = 10004
 integer, parameter :: BUFR_NEPN        = 10051
 integer, parameter :: bufr_dewPoint2m  = 12006 ! Dewpoint temperature at 2m
 integer, parameter :: BUFR_NBT1        = 12062
 integer, parameter :: BUFR_NBT2        = 12063
 integer, parameter :: BUFR_NBT3        = 12163
 integer, parameter :: BUFR_NEBD        = 15037
 integer, parameter :: BUFR_NERF        = 15036
 integer, parameter :: BUFR_NEHU        = 13210
 integer, parameter :: BUFR_NEHS        = 13220
 integer, parameter :: BUFR_NEZD        = 15031
 integer, parameter :: BUFR_NEFE        = 15032
 integer, parameter :: BUFR_NEZW        = 15035
 integer, parameter :: BUFR_ZTDSCORE    = 40026
 integer, parameter :: BUFR_NEAZ        = 05021
 integer, parameter :: BUFR_NEAL        = 40030 ! aladin HLOS wind
 integer, parameter :: BUFR_NEDWDP      = 40032 ! derivative of HLOS wrt P
 integer, parameter :: BUFR_NEDWDT      = 40033 ! derivative of HLOS wrt T
 integer, parameter :: BUFR_NEDW        = 11200 ! Doppler wind
 integer, parameter :: BUFR_radarPrecip = 21036 ! radar precipitation
 integer, parameter :: BUFR_logRadarPrecip = 51036 ! radar precipitation
 integer, parameter :: bufr_sst         = 22042 ! sea/water temperature
 integer, parameter :: bufr_soz         = 07025
 integer, parameter :: BUFR_ICEC        = 20237 ! concentration (%) from ice charts
 integer, parameter :: BUFR_ICEP        = 20222 ! concentration (%) from passive microwave retrievals
 integer, parameter :: BUFR_ICEV        = 21169 ! presence of ice retrieval from Vis/IR
 integer, parameter :: BUFR_ICES        = 21199 ! backscatter anisotropy from scatterometer
 integer, parameter :: bufr_vis         = 20001 ! horizontal visibility
 integer, parameter :: bufr_logVis      = 50001 ! log(horizontal visibility)
 integer, parameter :: bufr_gust        = 11041
 integer, parameter :: bufr_riverFlow   = 23040
 integer, parameter :: bufr_cloudInSeg  = 20081
 integer, parameter :: bufr_radvel      = 21014 ! Doppler velocity (Radial Wind)

 !
 ! Table C-14: (for Code Table 08046)
 !             Atmospheric chemical or physical constituent type
 !             - Includes local values as well!
 integer, parameter :: BUFR_NECH_O3    =0   
 integer, parameter :: BUFR_NECH_H2O   =1   
 integer, parameter :: BUFR_NECH_CH4   =2   
 integer, parameter :: BUFR_NECH_CO2   =3   
 integer, parameter :: BUFR_NECH_CO    =4   
 integer, parameter :: BUFR_NECH_NO2   =5   
 integer, parameter :: BUFR_NECH_N2O   =6   
 integer, parameter :: BUFR_NECH_HCHO  =7
 integer, parameter :: BUFR_NECH_SO2   =8  
 integer, parameter :: BUFR_NECH_NH3   =9  
 integer, parameter :: BUFR_NECH_NO    =11
 integer, parameter :: BUFR_NECH_PM25  =26 
 integer, parameter :: BUFR_NECH_PM10  =27  

 ! Table B elements associated to constituents/chemistry
 ! -----------------------------------------------------
 !
 ! Element denoting exponents accompanying obs bufr element
 integer, parameter :: BUFR_SCALE_EXPONENT      = 8090
 
 ! Element used for denoting averaging kernel elements 
 integer, parameter :: BUFR_AVGKERN      = 25143
 
 ! Element denoting error correlation matrix elements
 integer, parameter :: BUFR_CORREL      = 33205
 
 ! Elements for units of constituent observations
 !
 ! When needed, used in tandem with BUFR_NECH_* identifying the consituent
 ! and BUFR_SCALE_EXPONENT when the input file contains an exponent for
 ! a power ten scale factor. 
 !
 ! For some units, multiple BUFR are available, some originally defined to
 ! be specific to certain constituents and some being local elements devised
 ! prior to official elements being assigned. While a single set would be
 ! sufficient, all are included for completeness. A sufficient set might be
 ! 15008,15026,13002,15022,15207,15027,15003,15029,15198,15201,15021,15202,15024,15028.

 integer, parameter :: BUFR_UNIT_VMR            = 15008   ! Volume mixing ratio (vmr)
 integer, parameter :: BUFR_UNIT_VMR2           = 15208   ! Volume mixing ratio  
 integer, parameter :: BUFR_UNIT_MolePerMole    = 15026   ! Pollutant concentration (mole/mole)
 integer, parameter :: BUFR_UNIT_MolePerMole2   = 15197   ! Mixing ratio (mole/mole)
 integer, parameter :: BUFR_UNIT_MMR            = 13002   ! Mass mixing ratio (kg/kg)
 integer, parameter :: BUFR_UNIT_MMR2           = 13001   ! Humidity mass mixing ratio (kg/kg) - same as above
 integer, parameter :: BUFR_UNIT_NumberDensity  = 15207   ! Number density (1/m^3) 
 integer, parameter :: BUFR_UNIT_MolarDensity   = 15230   ! Molar density (mole/m^3)  
 integer, parameter :: BUFR_UNIT_Density        = 15027   ! Concentration in kg/m^3
 integer, parameter :: BUFR_UNIT_Density2       = 15223   ! Concentration = Density (kg/m^3) 
 integer, parameter :: BUFR_UNIT_AirDensity     = 15194   ! Air density (kg/m^3) 
 integer, parameter :: BUFR_UNIT_PMDensity      = 15195   ! Density of PM2.5 (kg/m^3)  
 integer, parameter :: BUFR_UNIT_PartPress      = 15003   ! Partial pressure in Pa (not just for ozone)
 integer, parameter :: BUFR_UNIT_PartPress2     = 15199   ! Partial pressure in Pa (same as above)
 integer, parameter :: BUFR_UNIT_MR_NVaerosol   = 15055   ! Non-volatile aerosol mixing ratio (unitless)  
 integer, parameter :: BUFR_UNIT_ExtinctCoef    = 15029   ! Extinction coefficient (1/m)
 
 ! The following elements are all vertically integrated quantities

 integer, parameter :: BUFR_UNIT_DU             = 15198   ! Integrate number density in Dobson units DU
 integer, parameter :: BUFR_UNIT_DU2            = 15001   ! Total ozone in DU (same as above) - applicable for all
 integer, parameter :: BUFR_UNIT_DU3            = 15005   ! Partial column for ozone in DU (same as above) - applicable for all
 integer, parameter :: BUFR_UNIT_DU4            = 15045   ! Partial column for SO2 in DU (same as above) - applicable for all
 integer, parameter :: BUFR_UNIT_IntegND        = 15201   ! Integrated number density (1/m^2)
 integer, parameter :: BUFR_UNIT_IntegND2       = 15012   ! Electron density per m^2 (1/m^2) 
 integer, parameter :: BUFR_UNIT_IntegDens      = 15021   ! Integrated density (kg/m^2)
 integer, parameter :: BUFR_UNIT_IntegDens2     = 15020   ! Integrated density for ozone (kg/m^2; same as above) - applicable to all
 integer, parameter :: BUFR_UNIT_IntegDens3     = 15200   ! Integrated density (kg/m^2)
 integer, parameter :: BUFR_UNIT_IntegMolarDens = 15202   ! Integrated molar density (mole/m^2)

 integer, parameter :: BUFR_UNIT_OptDepth       = 15024   ! Optical depth (unitless)
 integer, parameter :: BUFR_UNIT_OptDepth2      = 15196   ! Optical depth (unitless)
 integer, parameter :: BUFR_UNIT_OptDepth3      = 15062   ! Aerosol Optical depth (unitless)

 ! Additional elements associated to chemistry
 
 integer, parameter :: BUFR_UNIT_PhotoDissoc    = 15028   ! Photodissociate rate (1/sec)
 
 ! ----------------------------------------------
  
contains

  function bufr_IsAtmosConstituent(varNumber) result(var_chm)
    !
    !:Purpose: To determine whether 'varNumber' refers to constituent data from
    !          the CH family with recognized data units.
    !
    implicit none

    ! Arguments:
    integer, intent(in)           :: varNumber ! BUFR element number
    ! Result:
    logical                       :: var_chm
      
    if (any(varNumber.eq. (/ BUFR_UNIT_VMR, BUFR_UNIT_VMR2, BUFR_UNIT_MolePerMole, BUFR_UNIT_MolePerMole2, &
                             BUFR_UNIT_MMR, BUFR_UNIT_MMR2, BUFR_UNIT_NumberDensity, BUFR_UNIT_MolarDensity,  &
                             BUFR_UNIT_Density, BUFR_UNIT_Density2, &
                             BUFR_UNIT_AirDensity, BUFR_UNIT_PMDensity, &
                             BUFR_UNIT_OptDepth, BUFR_UNIT_OptDepth2, BUFR_UNIT_OptDepth3, BUFR_UNIT_MR_NVaerosol, &
                             BUFR_UNIT_PartPress, BUFR_UNIT_PartPress2, &
                             BUFR_UNIT_DU, BUFR_UNIT_DU2, BUFR_UNIT_DU3, BUFR_UNIT_DU4, &
                             BUFR_UNIT_IntegND, BUFR_UNIT_IntegND2, &
                             BUFR_UNIT_IntegDens, BUFR_UNIT_IntegDens2, BUFR_UNIT_IntegDens3, &
                             BUFR_UNIT_IntegMolarDens, BUFR_UNIT_ExtinctCoef, BUFR_UNIT_PhotoDissoc /) )) then          
       var_chm=.true.
    else         
       var_chm=.false.
    end if
      
  end function bufr_IsAtmosConstituent


  logical function bufr_IsIntegral(varNumber)
    !
    !:Purpose: To identify whether obs is a vertically integrated constituent
    !          measurement.
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: varNumber ! BUFR element number
 
    if (any(varNumber .eq. (/ BUFR_UNIT_DU, BUFR_UNIT_DU2, BUFR_UNIT_DU3, BUFR_UNIT_DU4, &
                              BUFR_UNIT_IntegND, BUFR_UNIT_IntegND2, BUFR_UNIT_IntegDens, &
                              BUFR_UNIT_IntegDens2, BUFR_UNIT_IntegDens3, BUFR_UNIT_IntegMolarDens /) )) then      
      bufr_IsIntegral=.true.     
    else
      bufr_IsIntegral=.false.
    end if
  
  end function bufr_IsIntegral


  logical function bufr_isWindComponent(varNumber)
    !
    !:Purpose: True if the variable is a wind component
    !
    implicit none
    
    ! Arguments:
    integer, intent(in) :: varNumber ! BUFR element number

    select case(varNumber)
    case(BUFR_NEUU, BUFR_NEVV, BUFR_NEUS, BUFR_NEVS)
      bufr_isWindComponent=.true.
    case default
      bufr_isWindComponent=.false.
    end select

  end function bufr_isWindComponent

  
  logical function bufr_isOceanObs(varNumber)
    !
    !:Purpose: True if the variable is an ocean observation
    !
    implicit none
    
    ! Arguments:
    integer, intent(in) :: varNumber ! BUFR element number

    select case(varNumber)
    case(BUFR_SST)
      bufr_isOceanObs = .true.
    case default
      bufr_isOceanObs = .false.
    end select

  end function bufr_isOceanObs

end module bufr_mod
