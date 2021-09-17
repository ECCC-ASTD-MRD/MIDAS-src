  
! This file is considered to be the source of the contained module for the
! assimilation code.  It can be noted, however, that it was produced with the aid
! of the script, thermoconsts_to_MathPhysConstants.py.  That script is simply a
! convenient tool that may or may not be used to update this file.   This is the
! command that produces this file:
!          ./thermoconsts_to_MathPhysConstants.py

MODULE MathPhysConstants_mod
  ! MODULE MathPhysConstants_mod (prefix='mpc' category='8. Global constants and interfaces')  
  !
  ! :Purpose: To supply mathematical and physical constants in a universal,
  !           reliable fashion.
  !
  implicit none
  public

! <<<<<< F U N D A M E N T A L   C O N S T A N T S >>>>>>
   real(8), parameter ::MPC_PI_R8                         =  3.141592653589793D+00 ! (unitless)
   real(8), parameter ::MPC_SPEED_OF_LIGHT_R8             =  2.997924580000000D+08 ! m/s - speed of light
   real(8), parameter ::MPC_PLANCK_R8                     =  6.626075500000000D-34 ! J/s - Planck's constant
   real(8), parameter ::MPC_BOLTZMANN_R8                  =  1.380658000000000D-23 ! J/K - Boltzmann constant
   real(8), parameter ::MPC_AVOGADRO_R8                   =  6.022136700000000D+23 ! mol-1 - Avogadro's number
!
! *** UNIT-CONVERSION RATIOS ***
   real(8), parameter ::MPC_KNOTS_PER_M_PER_S_R8          =  1.942535902919826D+00 ! knots/(m/s) - conversion: m/s to knots
! N.B.:  The correct value for MPC_M_PER_S_PER_KNOT is 0.5144444....
!        because 1 knot == 1.852 km/hr; 1 m/s = 3.6 km/hr
! However, this value must match that used by others
   real(8), parameter ::MPC_M_PER_S_PER_KNOT_R8           =  5.147910000000000D-01 ! (m/s)/knot - conversion: knots to m/s
   real(8), parameter ::MPC_PA_PER_MBAR_R8                =  1.000000000000000D+02 ! Pa/mbar - conversion: mbars to Pascals
   real(8), parameter ::MPC_MBAR_PER_PA_R8                =  1.000000000000000D-02 ! mbar/Pa - conversion: Pascals to mbars
   real(8), parameter ::MPC_RADIANS_PER_DEGREE_R8         =  1.745329251994330D-02 ! rad/deg - conversion: degrees to radians
   real(8), parameter ::MPC_DEGREES_PER_RADIAN_R8         =  5.729577951308232D+01 ! deg/rad - conversion: radians to degrees
   real(8), parameter ::MPC_K_C_DEGREE_OFFSET_R8          =  2.731500000000000D+02 ! K - offset between degrees K and C
!
! <<<<<< V A L U E S   T H A T   A R E   C O N S T A N T  ...  F O R   N O W >>>>>>
!
! *** ASTRONOMICAL CONSTANTS ***
   real(8)            ::MPC_ASTRONOMICAL_UNIT_R8          =  1.495978700000000D+11 ! m - half g. axis
   real(8)            ::MPC_MEAN_ANOMALY_R8               =  4.090930000000000D-01 ! (unitless) - mean anomaly
   real(8)            ::MPC_SIDEREAL_YEAR_R8              =  3.155815000000000D+07 ! s - sidereal year
!
! *** EARTH, THE HEAVENLY BODY ***
!     These are best obtained directly from:
!          module modgps02wgs84const (constants), or
!          module modgps06gravity (functions)
!
! *** RADIATION CONSTANTS ***
   real(8)            ::MPC_STEFAN_BOLTZMANN_R8           =  5.669800000000000D-08 ! W m-2 K-4 - Stefan-Boltzmann constant
   real(8)            ::MPC_SOLAR_CONSTANT_R8             =  1.367000000000000D+03 ! W/m2 - solar constant
!
! *** THERMODYNAMIC CONSTANTS (GAS PHASE) ***
   real(8)            ::MPC_RGAS_IDEAL_R8                 =  8.314511211948600D+00 ! J mol-1 K-1 - ideal-gas constant
   real(8)            ::MPC_RGAS_DRY_AIR_R8               =  2.870500000000000D+02 ! J kg-1 K-1 - gas constant, dry air
   real(8)            ::MPC_RGAS_VAPOUR_R8                =  4.615100000000000D+02 ! J kg-1 K-1 - gas constant, water vapour
   real(8)            ::MPC_EPS1_R8                       =  6.219800221014000D-01 ! (unitless) - R(air) / R(vapour)
   real(8)            ::MPC_EPS2_R8                       =  3.780199778986000D-01 ! (unitless) - 1 - eps1
   real(8)            ::MPC_DELTA_R8                      =  6.077686814144000D-01 ! (unitless) - 1/eps1 - 1 (unitless) \u2013 [R(vapour) / R(air)] - 1
!
! *** THERMODYNAMIC CONSTANTS (MOLAR MASSES) ***
   real(8)            ::MPC_MOLAR_MASS_DRY_AIR_R8         =  2.896440000000000D+01 ! g/mol - dry-air molar mass
   real(8)            ::MPC_MOLAR_MASS_VAPOUR_R8          =  1.801530000000000D+01 ! g/mol - water-vapour molar mass
   real(8)            ::MPC_MOLAR_MASS_O3_R8              =  4.799820000000000D+01 ! g/mol - ozone molar mass
   real(8)            ::MPC_MOLAR_MASS_CH4_R8             =  1.604246000000000D+01 ! g/mol - methane molar mass
   real(8)            ::MPC_MOLAR_MASS_CO2_R8             =  4.400950000000000D+01 ! g/mol - CO2 molar mass
   real(8)            ::MPC_MOLAR_MASS_CO_R8              =  2.801010000000000D+01 ! g/mol - CO molar mass
   real(8)            ::MPC_MOLAR_MASS_NO2_R8             =  4.600550000000000D+01 ! g/mol - NO2 molar mass
   real(8)            ::MPC_MOLAR_MASS_NO_R8              =  3.000610000000000D+01 ! g/mol - NO molar mass
   real(8)            ::MPC_MOLAR_MASS_N2O_R8             =  4.401280000000000D+01 ! g/mol - N2O molar mass
   real(8)            ::MPC_MOLAR_MASS_HCHO_R8            =  3.002598000000000D+01 ! g/mol - Formaldehyde molar mass
   real(8)            ::MPC_MOLAR_MASS_SO2_R8             =  6.406380000000000D+01 ! g/mol - SO2 molar mass
   real(8)            ::MPC_MOLAR_MASS_NH3_R8             =  1.703052000000000D+01 ! g/mol - NH3 molar mass
!
! *** THERMODYNAMIC CONSTANTS (SPECIFIC HEATS) ***
   real(8)            ::MPC_CV_DRY_AIR_R8                 =  7.184100000000001D+02 ! J kg-1 K-1 - sp. heat(V) of dry air
   real(8)            ::MPC_CV_VAPOUR_R8                  =  1.407950000000000D+03 ! J kg-1 K-1 - sp. heat(V) of water vapour
   real(8)            ::MPC_CP_DRY_AIR_R8                 =  1.005460000000000D+03 ! J kg-1 K-1 - sp. heat(P) of dry air
   real(8)            ::MPC_CP_VAPOUR_R8                  =  1.869460000000000D+03 ! J kg-1 K-1 - sp. heat(P) of water vapour
   real(8)            ::MPC_CP_ICE_R8                     =  2.115300000000000D+03 ! J kg-1 K-1 - sp. heat(P?) of ice
   real(8)            ::MPC_KAPPA_R8                      =  2.854912179500000D-01 ! (unitless) - for dry air: Rgas / Cp
!
! *** THERMODYNAMIC CONSTANTS (ENTHALPIES) ***
   real(8)            ::MPC_HEAT_CONDENS_WATER_R8         =  2.501000000000000D+06 ! J/kg - heat of condensation at 0C (water)
   real(8)            ::MPC_HEAT_FUSION_WATER_R8          =  3.340000000000000D+05 ! J/kg - heat of fusion (water)
   real(8)            ::MPC_HEAT_SUBL_WATER_R8            =  2.834000000000000D+06 ! J/kg - heat of sublimation (water)
!
! *** THERMODYNAMIC CONSTANTS (FLUID DYNAMICS) ***
   real(8)            ::MPC_KARMAN_R8                     =  4.000000000000000D-01 ! (unitless) - von Karman constant
   real(8)            ::MPC_CRITICAL_RICHARDSON_R8        =  2.000000000000000D-01 ! (unitless) -critical Richardson number
   real(8)            ::MPC_DENSITY_WATER_R8              =  1.000000000000000D+03 ! kg m-3 - density of (liquid) water
   real(8)            ::MPC_SCHUMANN_NEWELL_LAPSE_RATE_R8 =  6.628486583943000D-04 ! K s2 m-2 - Schuman-Newell lapse rate
!
! *** THERMODYNAMIC CONSTANTS (OTHER CHARACTERISTICS OF WATER) ***
   real(8)            ::MPC_T_ICE_R8                      =  2.731600000000000D+02 ! K - ice temperature in the atmosphere
   real(8)            ::MPC_TRIPLE_POINT_R8               =  2.731600000000000D+02 ! K - triple point of water
!
! *** USED TO CALCULATE L/CP IN FUNC HTVOCP ***
! Consider removing these 7 variables that are 'USED TO CALCULATE L/CP IN FUNC HTVOCP'
   real(8)            ::MPC_AI_R8                         =  2.864887713087000D+03 !  
   real(8)            ::MPC_AW_R8                         =  3.135012829948000D+03 !  
   real(8)            ::MPC_BI_R8                         =  1.660931315020000D-01 !  
   real(8)            ::MPC_BW_R8                         =  2.367075766316000D+00 !  
   real(8)            ::MPC_SLP_R8                        =  6.666666666667000D-02 !  
   real(8)            ::MPC_T1S_R8                        =  2.731600000000000D+02 ! K
   real(8)            ::MPC_T2S_R8                        =  2.581600000000000D+02 ! K
!
! <<<<<< C O N S T A N T S   C O N C E R N I N G   T H E   L I M I T A T I O N S   O F   D I G I T A L   C A L C U L A T I O N >>>>>>
   real(8)            ::MPC_MINIMUM_HU_R8                 =  2.500000000000000D-06 !  
   real(8)            ::MPC_MAXIMUM_ES_R8                 =  3.000000000000000D+01 !  
   real(8)            ::MPC_MINIMUM_VIS_R8                =  1.000000000000000D+00 !
   real(8)            ::MPC_MAXIMUM_VIS_R8                =  1.500000000000000D+04 ! 15km
   real(8)            ::MPC_MINIMUM_PR_R8                 =  1.000000000000000D-04 ! 0.1 mm/h = 0.0001 m/h
   real(8)            ::MPC_MINIMUM_CH_R8                 =  1.000000000000000D-01 ! 0.1 micrograms/kg (for constituents)
   real(8)            ::MPC_MINIMUM_PM_R8                 =  0.000000000000000D+00 ! 0 micrograms/m^3 (for particulate matter)

! <<<<<< F U N D A M E N T A L   C O N S T A N T S >>>>>>
   real(4), parameter ::MPC_PI_R4                         =  3.141592653589793D+00 ! (unitless)
   real(4), parameter ::MPC_SPEED_OF_LIGHT_R4             =  2.997924580000000D+08 ! m/s - speed of light
   real(4), parameter ::MPC_PLANCK_R4                     =  6.626075500000000D-34 ! J/s - Planck's constant
   real(4), parameter ::MPC_BOLTZMANN_R4                  =  1.380658000000000D-23 ! J/K - Boltzmann constant
   real(4), parameter ::MPC_AVOGADRO_R4                   =  6.022136700000000D+23 ! mol-1 - Avogadro's number
!
! *** UNIT-CONVERSION RATIOS ***
   real(4), parameter ::MPC_KNOTS_PER_M_PER_S_R4          =  1.942535902919826D+00 ! knots/(m/s) - conversion: m/s to knots
! N.B.:  The correct value for MPC_M_PER_S_PER_KNOT is 0.5144444....
!        because 1 knot == 1.852 km/hr; 1 m/s = 3.6 km/hr
! However, this value must match that used by others
   real(4), parameter ::MPC_M_PER_S_PER_KNOT_R4           =  5.147910000000000D-01 ! (m/s)/knot - conversion: knots to m/s
   real(4), parameter ::MPC_PA_PER_MBAR_R4                =  1.000000000000000D+02 ! Pa/mbar - conversion: mbars to Pascals
   real(4), parameter ::MPC_MBAR_PER_PA_R4                =  1.000000000000000D-02 ! mbar/Pa - conversion: Pascals to mbars
   real(4), parameter ::MPC_RADIANS_PER_DEGREE_R4         =  1.745329251994330D-02 ! rad/deg - conversion: degrees to radians
   real(4), parameter ::MPC_DEGREES_PER_RADIAN_R4         =  5.729577951308232D+01 ! deg/rad - conversion: radians to degrees
   real(4), parameter ::MPC_K_C_DEGREE_OFFSET_R4          =  2.731500000000000D+02 ! K - offset between degrees K and C
!
! <<<<<< V A L U E S   T H A T   A R E   C O N S T A N T  ...  F O R   N O W >>>>>>
!
! *** ASTRONOMICAL CONSTANTS ***
   real(4)            ::MPC_ASTRONOMICAL_UNIT_R4          =  1.495978700000000D+11 ! m - half g. axis
   real(4)            ::MPC_MEAN_ANOMALY_R4               =  4.090930000000000D-01 ! (unitless) - mean anomaly
   real(4)            ::MPC_SIDEREAL_YEAR_R4              =  3.155815000000000D+07 ! s - sidereal year
!
! *** EARTH, THE HEAVENLY BODY ***
!     These are best obtained directly from:
!          module modgps02wgs84const (constants), or
!          module modgps06gravity (functions)
!
! *** RADIATION CONSTANTS ***
   real(4)            ::MPC_STEFAN_BOLTZMANN_R4           =  5.669800000000000D-08 ! W m-2 K-4 - Stefan-Boltzmann constant
   real(4)            ::MPC_SOLAR_CONSTANT_R4             =  1.367000000000000D+03 ! W/m2 - solar constant
!
! *** THERMODYNAMIC CONSTANTS (GAS PHASE) ***
   real(4)            ::MPC_RGAS_IDEAL_R4                 =  8.314511211948600D+00 ! J mol-1 K-1 - ideal-gas constant
   real(4)            ::MPC_RGAS_DRY_AIR_R4               =  2.870500000000000D+02 ! J kg-1 K-1 - gas constant, dry air
   real(4)            ::MPC_RGAS_VAPOUR_R4                =  4.615100000000000D+02 ! J kg-1 K-1 - gas constant, water vapour
   real(4)            ::MPC_EPS1_R4                       =  6.219800221014000D-01 ! (unitless) - R(air) / R(vapour)
   real(4)            ::MPC_EPS2_R4                       =  3.780199778986000D-01 ! (unitless) - 1 - eps1
   real(4)            ::MPC_DELTA_R4                      =  6.077686814144000D-01 ! (unitless) - 1/eps1 - 1 (unitless) \u2013 [R(vapour) / R(air)] - 1
!
! *** THERMODYNAMIC CONSTANTS (MOLAR MASSES) ***
   real(4)            ::MPC_MOLAR_MASS_DRY_AIR_R4         =  2.896440000000000D+01 ! g/mol - dry-air molar mass
   real(4)            ::MPC_MOLAR_MASS_VAPOUR_R4          =  1.801530000000000D+01 ! g/mol - water-vapour molar mass
   real(4)            ::MPC_MOLAR_MASS_O3_R4              =  4.799820000000000D+01 ! g/mol - ozone molar mass
   real(4)            ::MPC_MOLAR_MASS_CH4_R4             =  1.604246000000000D+01 ! g/mol - methane molar mass
   real(4)            ::MPC_MOLAR_MASS_CO2_R4             =  4.400950000000000D+01 ! g/mol - CO2 molar mass
   real(4)            ::MPC_MOLAR_MASS_CO_R4              =  2.801010000000000D+01 ! g/mol - CO molar mass
   real(4)            ::MPC_MOLAR_MASS_NO2_R4             =  4.600550000000000D+01 ! g/mol - NO2 molar mass
   real(4)            ::MPC_MOLAR_MASS_NO_R4              =  3.000610000000000D+01 ! g/mol - NO molar mass
   real(4)            ::MPC_MOLAR_MASS_N2O_R4             =  4.401280000000000D+01 ! g/mol - N2O molar mass
   real(4)            ::MPC_MOLAR_MASS_HCHO_R4            =  3.002598000000000D+01 ! g/mol - Formaldehyde molar mass
   real(4)            ::MPC_MOLAR_MASS_SO2_R4             =  6.406380000000000D+01 ! g/mol - SO2 molar mass
   real(4)            ::MPC_MOLAR_MASS_NH3_R4             =  1.703052000000000D+01 ! g/mol - NH3 molar mass
!
! *** THERMODYNAMIC CONSTANTS (SPECIFIC HEATS) ***
   real(4)            ::MPC_CV_DRY_AIR_R4                 =  7.184100000000001D+02 ! J kg-1 K-1 - sp. heat(V) of dry air
   real(4)            ::MPC_CV_VAPOUR_R4                  =  1.407950000000000D+03 ! J kg-1 K-1 - sp. heat(V) of water vapour
   real(4)            ::MPC_CP_DRY_AIR_R4                 =  1.005460000000000D+03 ! J kg-1 K-1 - sp. heat(P) of dry air
   real(4)            ::MPC_CP_VAPOUR_R4                  =  1.869460000000000D+03 ! J kg-1 K-1 - sp. heat(P) of water vapour
   real(4)            ::MPC_CP_ICE_R4                     =  2.115300000000000D+03 ! J kg-1 K-1 - sp. heat(P?) of ice
   real(4)            ::MPC_KAPPA_R4                      =  2.854912179500000D-01 ! (unitless) - for dry air: Rgas / Cp
!
! *** THERMODYNAMIC CONSTANTS (ENTHALPIES) ***
   real(4)            ::MPC_HEAT_CONDENS_WATER_R4         =  2.501000000000000D+06 ! J/kg - heat of condensation at 0C (water)
   real(4)            ::MPC_HEAT_FUSION_WATER_R4          =  3.340000000000000D+05 ! J/kg - heat of fusion (water)
   real(4)            ::MPC_HEAT_SUBL_WATER_R4            =  2.834000000000000D+06 ! J/kg - heat of sublimation (water)
!
! *** THERMODYNAMIC CONSTANTS (FLUID DYNAMICS) ***
   real(4)            ::MPC_KARMAN_R4                     =  4.000000000000000D-01 ! (unitless) - von Karman constant
   real(4)            ::MPC_CRITICAL_RICHARDSON_R4        =  2.000000000000000D-01 ! (unitless) -critical Richardson number
   real(4)            ::MPC_DENSITY_WATER_R4              =  1.000000000000000D+03 ! kg m-3 - density of (liquid) water
   real(4)            ::MPC_SCHUMANN_NEWELL_LAPSE_RATE_R4 =  6.628486583943000D-04 ! K s2 m-2 - Schuman-Newell lapse rate
!
! *** THERMODYNAMIC CONSTANTS (OTHER CHARACTERISTICS OF WATER) ***
   real(4)            ::MPC_T_ICE_R4                      =  2.731600000000000D+02 ! K - ice temperature in the atmosphere
   real(4)            ::MPC_TRIPLE_POINT_R4               =  2.731600000000000D+02 ! K - triple point of water
!
! *** USED TO CALCULATE L/CP IN FUNC HTVOCP ***
! Consider removing these 7 variables that are 'USED TO CALCULATE L/CP IN FUNC HTVOCP'
   real(4)            ::MPC_AI_R4                         =  2.864887713087000D+03 !  
   real(4)            ::MPC_AW_R4                         =  3.135012829948000D+03 !  
   real(4)            ::MPC_BI_R4                         =  1.660931315020000D-01 !  
   real(4)            ::MPC_BW_R4                         =  2.367075766316000D+00 !  
   real(4)            ::MPC_SLP_R4                        =  6.666666666667000D-02 !  
   real(4)            ::MPC_T1S_R4                        =  2.731600000000000D+02 ! K
   real(4)            ::MPC_T2S_R4                        =  2.581600000000000D+02 ! K
!
! <<<<<< C O N S T A N T S   C O N C E R N I N G   T H E   L I M I T A T I O N S   O F   D I G I T A L   C A L C U L A T I O N >>>>>>
   real(4)            ::MPC_MINIMUM_HU_R4                 =  2.500000000000000D-06 !  
   real(4)            ::MPC_MAXIMUM_ES_R4                 =  3.000000000000000D+01 !  
   real(4)            ::MPC_MINIMUM_VIS_R4                =  1.000000000000000D+00 !
   real(4)            ::MPC_MAXIMUM_VIS_R4                =  1.500000000000000D+04 ! 15km
   real(4)            ::MPC_MINIMUM_PR_R4                 =  1.000000000000000D-04 ! 0.1 mm/h = 0.0001 m/h
   real(4)            ::MPC_MINIMUM_CH_R4                 =  1.000000000000000D-01 ! 0.1 micrograms/kg (for constituents)
   real(4)            ::MPC_MINIMUM_PM_R4                 =  0.000000000000000D+00 ! 0 micrograms/m^3 (for particulate matter)

! <<<<<< OBS FILE CONSTANTS
   real(4), parameter :: MPC_missingValue_R4 = -999.
   real(8), parameter :: MPC_missingValue_R8 = -999.D0
   integer, parameter :: MPC_missingValue_INT = -999

contains

  subroutine mpc_setValue(name, var_r4, var_r8, value)
     !
     !:Purpose: To provide a means to change a (non-parameter) value
     implicit none

     ! Arguments:
     character(len=*) :: name
     real(4) :: var_r4
     real(8) :: var_r8, value

     var_r4 = value
     var_r8 = value

     write(6,*)'****************************************************'
     write(6,*)'*'
     write(6,*)'*    THE VALUE OF THE CONSTANT ', name
     write(6,*)'*'
     write(6,'(A25, D22.15)')' *    HAS BEEN CHANGED TO ', var_r8
     write(6,*)'*'
     write(6,*)'****************************************************'
  end subroutine mpc_setValue

  subroutine mpc_printConstants(kulout)
     !
     !:Purpose: To print all of the constants that are provided by this module.
     !          The intent is to make it clear in a program listing which values
     !          were used.
     implicit none

     ! Arguments:
     integer, intent(in) :: kulout     ! unit number for printing

     write(kulout,FMT='(//,4x ,"*** mpc_printConstants: definition of Mathematical and Physical constants  ***",/)')
     write(kulout,*) " <<<<<< F U N D A M E N T A L   C O N S T A N T S >>>>>>"
     write(kulout,"(A36, D22.15, A)") "                         MPC_PI_R8= ", MPC_PI_R8   , " (unitless)"
     write(kulout,"(A36, D22.15, A)") "             MPC_SPEED_OF_LIGHT_R8= ", MPC_SPEED_OF_LIGHT_R8   , " m/s - speed of light"
     write(kulout,"(A36, D22.15, A)") "                     MPC_PLANCK_R8= ", MPC_PLANCK_R8   , " J/s - Planck's constant"
     write(kulout,"(A36, D22.15, A)") "                  MPC_BOLTZMANN_R8= ", MPC_BOLTZMANN_R8   , " J/K - Boltzmann constant"
     write(kulout,"(A36, D22.15, A)") "                   MPC_AVOGADRO_R8= ", MPC_AVOGADRO_R8   , " mol-1 - Avogadro's number"
     write(kulout,*) ""
     write(kulout,*) " *** UNIT-CONVERSION RATIOS ***"
     write(kulout,"(A36, D22.15, A)") "          MPC_KNOTS_PER_M_PER_S_R8= ", MPC_KNOTS_PER_M_PER_S_R8   , " knots/(m/s) - conversion: m/s to knots"
     write(kulout,*) " N.B.:  The correct value for MPC_M_PER_S_PER_KNOT is 0.5144444...."
     write(kulout,*) "        because 1 knot == 1.852 km/hr; 1 m/s = 3.6 km/hr"
     write(kulout,*) " However, this value must match that used by others"
     write(kulout,"(A36, D22.15, A)") "           MPC_M_PER_S_PER_KNOT_R8= ", MPC_M_PER_S_PER_KNOT_R8   , " (m/s)/knot - conversion: knots to m/s"
     write(kulout,"(A36, D22.15, A)") "                MPC_PA_PER_MBAR_R8= ", MPC_PA_PER_MBAR_R8   , " Pa/mbar - conversion: mbars to Pascals"
     write(kulout,"(A36, D22.15, A)") "                MPC_MBAR_PER_PA_R8= ", MPC_MBAR_PER_PA_R8   , " mbar/Pa - conversion: Pascals to mbars"
     write(kulout,"(A36, D22.15, A)") "         MPC_RADIANS_PER_DEGREE_R8= ", MPC_RADIANS_PER_DEGREE_R8   , " rad/deg - conversion: degrees to radians"
     write(kulout,"(A36, D22.15, A)") "         MPC_DEGREES_PER_RADIAN_R8= ", MPC_DEGREES_PER_RADIAN_R8   , " deg/rad - conversion: radians to degrees"
     write(kulout,"(A36, D22.15, A)") "          MPC_K_C_DEGREE_OFFSET_R8= ", MPC_K_C_DEGREE_OFFSET_R8   , " K - offset between degrees K and C"
     write(kulout,*) ""
     write(kulout,*) " <<<<<< V A L U E S   T H A T   A R E   C O N S T A N T  ...  F O R   N O W >>>>>>"
     write(kulout,*) ""
     write(kulout,*) " *** ASTRONOMICAL CONSTANTS ***"
     write(kulout,"(A36, D22.15, A)") "          MPC_ASTRONOMICAL_UNIT_R8= ", MPC_ASTRONOMICAL_UNIT_R8   , " m - half g. axis"
     write(kulout,"(A36, D22.15, A)") "               MPC_MEAN_ANOMALY_R8= ", MPC_MEAN_ANOMALY_R8   , " (unitless) - mean anomaly"
     write(kulout,"(A36, D22.15, A)") "              MPC_SIDEREAL_YEAR_R8= ", MPC_SIDEREAL_YEAR_R8   , " s - sidereal year"
     write(kulout,*) ""
     write(kulout,*) " *** EARTH, THE HEAVENLY BODY ***"
     write(kulout,*) "     These are best obtained directly from:"
     write(kulout,*) "          module modgps02wgs84const (constants), or"
     write(kulout,*) "          module modgps06gravity (functions)"
     write(kulout,*) ""
     write(kulout,*) " *** RADIATION CONSTANTS ***"
     write(kulout,"(A36, D22.15, A)") "           MPC_STEFAN_BOLTZMANN_R8= ", MPC_STEFAN_BOLTZMANN_R8   , " W m-2 K-4 - Stefan-Boltzmann constant"
     write(kulout,"(A36, D22.15, A)") "             MPC_SOLAR_CONSTANT_R8= ", MPC_SOLAR_CONSTANT_R8   , " W/m2 - solar constant"
     write(kulout,*) ""
     write(kulout,*) " *** THERMODYNAMIC CONSTANTS (GAS PHASE) ***"
     write(kulout,"(A36, D22.15, A)") "                 MPC_RGAS_IDEAL_R8= ", MPC_RGAS_IDEAL_R8   , " J mol-1 K-1 - ideal-gas constant"
     write(kulout,"(A36, D22.15, A)") "               MPC_RGAS_DRY_AIR_R8= ", MPC_RGAS_DRY_AIR_R8   , " J kg-1 K-1 - gas constant, dry air"
     write(kulout,"(A36, D22.15, A)") "                MPC_RGAS_VAPOUR_R8= ", MPC_RGAS_VAPOUR_R8   , " J kg-1 K-1 - gas constant, water vapour"
     write(kulout,"(A36, D22.15, A)") "                       MPC_EPS1_R8= ", MPC_EPS1_R8   , " (unitless) - R(air) / R(vapour)"
     write(kulout,"(A36, D22.15, A)") "                       MPC_EPS2_R8= ", MPC_EPS2_R8   , " (unitless) - 1 - eps1"
     write(kulout,"(A36, D22.15, A)") "                      MPC_DELTA_R8= ", MPC_DELTA_R8   , " (unitless) - 1/eps1 - 1 (unitless) \u2013 [R(vapour) / R(air)] - 1"
     write(kulout,*) ""
     write(kulout,*) " *** THERMODYNAMIC CONSTANTS (MOLAR MASSES) ***"
     write(kulout,"(A36, D22.15, A)") "         MPC_MOLAR_MASS_DRY_AIR_R8= ", MPC_MOLAR_MASS_DRY_AIR_R8   , " g/mol - dry-air molar mass"
     write(kulout,"(A36, D22.15, A)") "          MPC_MOLAR_MASS_VAPOUR_R8= ", MPC_MOLAR_MASS_VAPOUR_R8   , " g/mol - water-vapour molar mass"
     write(kulout,"(A36, D22.15, A)") "              MPC_MOLAR_MASS_O3_R8= ", MPC_MOLAR_MASS_O3_R8   , " g/mol - ozone molar mass"
     write(kulout,"(A36, D22.15, A)") "             MPC_MOLAR_MASS_CH4_R8= ", MPC_MOLAR_MASS_CH4_R8   , " g/mol - methane molar mass"
     write(kulout,"(A36, D22.15, A)") "             MPC_MOLAR_MASS_CO2_R8= ", MPC_MOLAR_MASS_CO2_R8   , " g/mol - CO2 molar mass"
     write(kulout,"(A36, D22.15, A)") "              MPC_MOLAR_MASS_CO_R8= ", MPC_MOLAR_MASS_CO_R8   , " g/mol - CO molar mass"
     write(kulout,"(A36, D22.15, A)") "             MPC_MOLAR_MASS_NO2_R8= ", MPC_MOLAR_MASS_NO2_R8   , " g/mol - NO2 molar mass"
     write(kulout,"(A36, D22.15, A)") "              MPC_MOLAR_MASS_NO_R8= ", MPC_MOLAR_MASS_NO_R8   , " g/mol - NO molar mass"
     write(kulout,"(A36, D22.15, A)") "             MPC_MOLAR_MASS_N2O_R8= ", MPC_MOLAR_MASS_N2O_R8   , " g/mol - N2O molar mass"
     write(kulout,"(A36, D22.15, A)") "            MPC_MOLAR_MASS_HCHO_R8= ", MPC_MOLAR_MASS_HCHO_R8   , " g/mol - Formaldehyde molar mass"
     write(kulout,"(A36, D22.15, A)") "             MPC_MOLAR_MASS_SO2_R8= ", MPC_MOLAR_MASS_SO2_R8   , " g/mol - SO2 molar mass"
     write(kulout,"(A36, D22.15, A)") "             MPC_MOLAR_MASS_NH3_R8= ", MPC_MOLAR_MASS_NH3_R8   , " g/mol - NH3 molar mass"
     write(kulout,*) ""
     write(kulout,*) " *** THERMODYNAMIC CONSTANTS (SPECIFIC HEATS) ***"
     write(kulout,"(A36, D22.15, A)") "                 MPC_CV_DRY_AIR_R8= ", MPC_CV_DRY_AIR_R8   , " J kg-1 K-1 - sp. heat(V) of dry air"
     write(kulout,"(A36, D22.15, A)") "                  MPC_CV_VAPOUR_R8= ", MPC_CV_VAPOUR_R8   , " J kg-1 K-1 - sp. heat(V) of water vapour"
     write(kulout,"(A36, D22.15, A)") "                 MPC_CP_DRY_AIR_R8= ", MPC_CP_DRY_AIR_R8   , " J kg-1 K-1 - sp. heat(P) of dry air"
     write(kulout,"(A36, D22.15, A)") "                  MPC_CP_VAPOUR_R8= ", MPC_CP_VAPOUR_R8   , " J kg-1 K-1 - sp. heat(P) of water vapour"
     write(kulout,"(A36, D22.15, A)") "                     MPC_CP_ICE_R8= ", MPC_CP_ICE_R8   , " J kg-1 K-1 - sp. heat(P?) of ice"
     write(kulout,"(A36, D22.15, A)") "                      MPC_KAPPA_R8= ", MPC_KAPPA_R8   , " (unitless) - for dry air: Rgas / Cp"
     write(kulout,*) ""
     write(kulout,*) " *** THERMODYNAMIC CONSTANTS (ENTHALPIES) ***"
     write(kulout,"(A36, D22.15, A)") "         MPC_HEAT_CONDENS_WATER_R8= ", MPC_HEAT_CONDENS_WATER_R8   , " J/kg - heat of condensation at 0C (water)"
     write(kulout,"(A36, D22.15, A)") "          MPC_HEAT_FUSION_WATER_R8= ", MPC_HEAT_FUSION_WATER_R8   , " J/kg - heat of fusion (water)"
     write(kulout,"(A36, D22.15, A)") "            MPC_HEAT_SUBL_WATER_R8= ", MPC_HEAT_SUBL_WATER_R8   , " J/kg - heat of sublimation (water)"
     write(kulout,*) ""
     write(kulout,*) " *** THERMODYNAMIC CONSTANTS (FLUID DYNAMICS) ***"
     write(kulout,"(A36, D22.15, A)") "                     MPC_KARMAN_R8= ", MPC_KARMAN_R8   , " (unitless) - von Karman constant"
     write(kulout,"(A36, D22.15, A)") "        MPC_CRITICAL_RICHARDSON_R8= ", MPC_CRITICAL_RICHARDSON_R8   , " (unitless) -critical Richardson number"
     write(kulout,"(A36, D22.15, A)") "              MPC_DENSITY_WATER_R8= ", MPC_DENSITY_WATER_R8   , " kg m-3 - density of (liquid) water"
     write(kulout,"(A36, D22.15, A)") " MPC_SCHUMANN_NEWELL_LAPSE_RATE_R8= ", MPC_SCHUMANN_NEWELL_LAPSE_RATE_R8   , " K s2 m-2 - Schuman-Newell lapse rate"
     write(kulout,*) ""
     write(kulout,*) " *** THERMODYNAMIC CONSTANTS (OTHER CHARACTERISTICS OF WATER) ***"
     write(kulout,"(A36, D22.15, A)") "                      MPC_T_ICE_R8= ", MPC_T_ICE_R8   , " K - ice temperature in the atmosphere"
     write(kulout,"(A36, D22.15, A)") "               MPC_TRIPLE_POINT_R8= ", MPC_TRIPLE_POINT_R8   , " K - triple point of water"
     write(kulout,*) ""
     write(kulout,*) " *** USED TO CALCULATE L/CP IN FUNC HTVOCP ***"
     write(kulout,*) " Consider removing these 7 variables that are 'USED TO CALCULATE L/CP IN FUNC HTVOCP'"
     write(kulout,"(A36, D22.15, A)") "                         MPC_AI_R8= ", MPC_AI_R8   , "  "
     write(kulout,"(A36, D22.15, A)") "                         MPC_AW_R8= ", MPC_AW_R8   , "  "
     write(kulout,"(A36, D22.15, A)") "                         MPC_BI_R8= ", MPC_BI_R8   , "  "
     write(kulout,"(A36, D22.15, A)") "                         MPC_BW_R8= ", MPC_BW_R8   , "  "
     write(kulout,"(A36, D22.15, A)") "                        MPC_SLP_R8= ", MPC_SLP_R8   , "  "
     write(kulout,"(A36, D22.15, A)") "                        MPC_T1S_R8= ", MPC_T1S_R8   , " K"
     write(kulout,"(A36, D22.15, A)") "                        MPC_T2S_R8= ", MPC_T2S_R8   , " K"
     write(kulout,*) ""
     write(kulout,*) " <<<<<< C O N S T A N T S   C O N C E R N I N G   T H E   L I M I T A T I O N S   O F   D I G I T A L   C A L C U L A T I O N >>>>>>"
     write(kulout,"(A36, D22.15, A)") "                 MPC_MINIMUM_HU_R8= ", MPC_MINIMUM_HU_R8   , "  "
     write(kulout,"(A36, D22.15, A)") "                 MPC_MAXIMUM_ES_R8= ", MPC_MAXIMUM_ES_R8   , "  "
     write(kulout,"(A36, D22.15, A)") "                 MPC_MINIMUM_CH_R8= ", MPC_MINIMUM_CH_R8   , "  "
     write(kulout,"(A36, D22.15, A)") "                 MPC_MINIMUM_PM_R8= ", MPC_MINIMUM_PM_R8   , "  "
     write(kulout,*) "\n\n"
  end subroutine mpc_printConstants
end MODULE MathPhysConstants_mod
