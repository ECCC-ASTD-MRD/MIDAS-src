#!/usr/bin/env python2

"""This utility is NOT part of the MathPhysConstants source code.  It is merely
a convenient tool for producing the file, MathPhysConstants_mod.ftn90, which is part of the
MathPhysConstants source code.

The purpose of this utility is to provide a simple means of ensuring that the
file, MathPhysConstants_mod.ftn90, is consistent with the file thermoconsts.  The utility reads the
latter file and incorporates that information into MathPhysConstants_mod.ftn90.  After executing this
utility, MathPhysConstants_mod.ftn90 should be saved in the revision-control data base of the source
code."""

"It can be noted that, if Python 3.1 or greater is used, almost all platforms map Python floats to IEEE-754 'double precision', the same representation that is used by fortran. IEEE-754 doubles contain 53 bits of precision.  Thus, values that are calculated here should yield double precision, without explicitly requesting it.  It can be further noted that 53 base-two digits is equivalent to 15.95 base-ten digits.  Thus, if the constants were printed with 16 digits of precision, the non-perfect base-two representation of the numbers would be evident; e.g. 4.0000000000000002D-01 instead of 4.0000000000000000D-01.  Therefore, this utility writes only 15 digits of precision, for the sake of clarity to humans."

from sys    import stdin,stdout

VALUE=1                                 # index to acces the value in the PhysDict


def init_phys_values():
    """Create a dictionary containing the initial values of the mathematical and
    physical constants / variables."""

                                        # In this dictionary:
                                        #    key   - the constant/variable name in assimilation
                                        #    value - a list containg a type, the value, and a comment (units - description)
                                        #          - type indicates whether this value should be constant or variable
                                        #          - for entries to be read from thermoconsts, the initial value is a message
                                        #          - for calculated values, the initial value is 'TBD'
    PhysDict={
              # <<<<<< F U N D A M E N T A L   C O N S T A N T S >>>>>>
              "MPC_PI"                         : ('const', 3.1415926535897932384, "(unitless)"),
              "MPC_SPEED_OF_LIGHT"             : ('const', 299792458.E+0, "m/s - speed of light"),
              "MPC_PLANCK"                     : ('const', 6.6260755E-34, "J/s - Planck's constant"),
              "MPC_BOLTZMANN"                  : ('const', 1.380658E-23, "J/K - Boltzmann constant"),
              "MPC_AVOGADRO"                   : ('const', 6.0221367E+23, "mol-1 - Avogadro's number"),

              # *** UNIT-CONVERSION RATIOS ***
              "MPC_KNOTS_PER_M_PER_S"          : ('const', 'TBD', "knots/(m/s) - conversion: m/s to knots"),
                                        # N.B.:  The correct value for MPC_M_PER_S_PER_KNOT is 0.5144444....
                                        #        because 1 knot == 1.852 km/hr; 1 m/s = 3.6 km/hr
                                        # However, this value must match that used by others
              "MPC_M_PER_S_PER_KNOT"           : ('const', 'reading thermoconsts FAILED', "(m/s)/knot - conversion: knots to m/s"),
              "MPC_PA_PER_MBAR"                : ('const', 1.0E+02, "Pa/mbar - conversion: mbars to Pascals"),
              "MPC_MBAR_PER_PA"                : ('const', 1.0E-02, "mbar/Pa - conversion: Pascals to mbars"),
              "MPC_RADIANS_PER_DEGREE"         : ('const', 'TBD', "rad/deg - conversion: degrees to radians"),
              "MPC_DEGREES_PER_RADIAN"         : ('const', 'TBD', "deg/rad - conversion: radians to degrees"),
              "MPC_K_C_DEGREE_OFFSET"          : ('const', 'reading thermoconsts FAILED', "K - offset between degrees K and C"),

              # <<<<<< V A L U E S   T H A T   A R E   C O N S T A N T  ...  F O R   N O W >>>>>>

              # *** ASTRONOMICAL CONSTANTS ***
              "MPC_ASTRONOMICAL_UNIT"          : ('var', 149597870000.E+0, "m - half g. axis"),
              "MPC_MEAN_ANOMALY"               : ('var', 0.409093E+0, "(unitless) - mean anomaly"),
              "MPC_SIDEREAL_YEAR"              : ('var', 0.3155815E+08, "s - sidereal year"),

              # *** EARTH, THE HEAVENLY BODY ***
              #     These are best obtained directly from:
              #          module modgps02wgs84const (constants), or
              #          module modgps06gravity (functions)

              # *** RADIATION CONSTANTS ***
              "MPC_STEFAN_BOLTZMANN"           : ('var', 'reading thermoconsts FAILED',
                                                         "W m-2 K-4 - Stefan-Boltzmann constant"),
              "MPC_SOLAR_CONSTANT"             : ('var', 'reading thermoconsts FAILED',
                                                         "W/m2 - solar constant"),

              # *** THERMODYNAMIC CONSTANTS (GAS PHASE) ***
              "MPC_RGAS_IDEAL"                 : ('var', 'TBD',
                                                         "J mol-1 K-1 - ideal-gas constant"),
              "MPC_RGAS_DRY_AIR"               : ('var', 'reading thermoconsts FAILED',
                                                         "J kg-1 K-1 - gas constant, dry air"),
              "MPC_RGAS_VAPOUR"                : ('var', 'reading thermoconsts FAILED',
                                                         "J kg-1 K-1 - gas constant, water vapour"),
              "MPC_EPS1"                       : ('var', 'reading thermoconsts FAILED',
                                                         "(unitless) - R(air) / R(vapour)"),
              "MPC_EPS2"                       : ('var', 'reading thermoconsts FAILED',
                                                         "(unitless) - 1 - eps1"),
              "MPC_DELTA"                      : ('var', 'reading thermoconsts FAILED',
                                                         "(unitless) - 1/eps1 - 1 (unitless) \u2013 [R(vapour) / R(air)] - 1"),

              # *** THERMODYNAMIC CONSTANTS (MOLAR MASSES) ***
              "MPC_MOLAR_MASS_DRY_AIR"         : ('var', 28.9644E+0,
                                                         "g/mol - dry-air molar mass"),
              "MPC_MOLAR_MASS_VAPOUR"          : ('var', 18.0153E+0,
                                                         "g/mol - water-vapour molar mass"),

              # *** THERMODYNAMIC CONSTANTS (SPECIFIC HEATS) ***
              "MPC_CV_DRY_AIR"                 : ('var', 'TBD',
                                                         "J kg-1 K-1 - sp. heat(V) of dry air"),
              "MPC_CV_VAPOUR"                  : ('var', 'TBD',
                                                         "J kg-1 K-1 - sp. heat(V) of water vapour"),
              "MPC_CP_DRY_AIR"                 : ('var', 'reading thermoconsts FAILED',
                                                         "J kg-1 K-1 - sp. heat(P) of dry air"),
              "MPC_CP_VAPOUR"                  : ('var', 'reading thermoconsts FAILED',
                                                         "J kg-1 K-1 - sp. heat(P) of water vapour"),
              "MPC_CP_ICE"                     : ('var', 0.2115300000000E+04,
                                                         "J kg-1 K-1 - sp. heat(P?) of ice"),
              "MPC_KAPPA"                      : ('var', 'reading thermoconsts FAILED',
                                                         "(unitless) - for dry air: Rgas / Cp"),

              # *** THERMODYNAMIC CONSTANTS (ENTHALPIES) ***
              "MPC_HEAT_CONDENS_WATER"         : ('var', 'reading thermoconsts FAILED',
                                                         "J/kg - heat of condensation at 0C (water)"),
              "MPC_HEAT_FUSION_WATER"          : ('var', 'reading thermoconsts FAILED',
                                                         "J/kg - heat of fusion (water)"),
              "MPC_HEAT_SUBL_WATER"            : ('var', 0.28340000000000E+07,
                                                         "J/kg - heat of sublimation (water)"),

              # *** THERMODYNAMIC CONSTANTS (FLUID DYNAMICS) ***
              "MPC_KARMAN"                     : ('var', 0.4000000000000E+00,
                                                         "(unitless) - von Karman constant"),
              "MPC_CRITICAL_RICHARDSON"        : ('var', 'reading thermoconsts FAILED',
                                                         "(unitless) -critical Richardson number"),
              "MPC_DENSITY_WATER"              : ('var', 'reading thermoconsts FAILED',
                                                         "kg m-3 - density of (liquid) water"),
              "MPC_SCHUMANN_NEWELL_LAPSE_RATE" : ('var', 'reading thermoconsts FAILED',
                                                         "K s2 m-2 - Schuman-Newell lapse rate"),

              # *** THERMODYNAMIC CONSTANTS (OTHER CHARACTERISTICS OF WATER) ***
              "MPC_T_ICE"                      : ('var', 'reading thermoconsts FAILED',
                                                         "K - ice temperature in the atmosphere"),
              "MPC_TRIPLE_POINT"               : ('var', 'reading thermoconsts FAILED',
                                                         "K - triple point of water"), 

              # *** USED TO CALCULATE L/CP IN FUNC HTVOCP ***
              "MPC_AI"                         : ('var', 'reading thermoconsts FAILED', " "),
              "MPC_AW"                         : ('var', 'reading thermoconsts FAILED', " "),
              "MPC_BI"                         : ('var', 'reading thermoconsts FAILED', " "),
              "MPC_BW"                         : ('var', 'reading thermoconsts FAILED', " "),
              "MPC_SLP"                        : ('var', 'reading thermoconsts FAILED', " "),
              "MPC_T1S"                        : ('var', 'reading thermoconsts FAILED', "K"),
              "MPC_T2S"                        : ('var', 'reading thermoconsts FAILED', "K"),

              # <<<<<< C O N S T A N T S   C O N C E R N I N G   T H E   L I M I T A T I O N S   O F   D I G I T A L   C A L C U L A T I O N >>>>>>
              "MPC_MINIMUM_HU"                 : ('var', 2.50E-06, " "),
              "MPC_MAXIMUM_ES"                 : ('var', 30.0E+00, " ")
             }

    # Make a list of the 'nice' order of entries in PhysDict
    PhysList=(
          "# <<<<<< F U N D A M E N T A L   C O N S T A N T S >>>>>>",
          "MPC_PI",
          "MPC_SPEED_OF_LIGHT",
          "MPC_PLANCK",
          "MPC_BOLTZMANN",
          "MPC_AVOGADRO",
          "#",
          "# *** UNIT-CONVERSION RATIOS ***",
          "MPC_KNOTS_PER_M_PER_S",
                                    "# N.B.:  The correct value for MPC_M_PER_S_PER_KNOT is 0.5144444....",
                                    "#        because 1 knot == 1.852 km/hr; 1 m/s = 3.6 km/hr",
                                    "# However, this value must match that used by others",
          "MPC_M_PER_S_PER_KNOT",
          "MPC_PA_PER_MBAR",
          "MPC_MBAR_PER_PA",
          "MPC_RADIANS_PER_DEGREE",
          "MPC_DEGREES_PER_RADIAN",
          "MPC_K_C_DEGREE_OFFSET",
          "#",
          "# <<<<<< V A L U E S   T H A T   A R E   C O N S T A N T  ...  F O R   N O W >>>>>>",
          "#",
          "# *** ASTRONOMICAL CONSTANTS ***",
          "MPC_ASTRONOMICAL_UNIT",
          "MPC_MEAN_ANOMALY",
          "MPC_SIDEREAL_YEAR",
          "#",
          "# *** EARTH, THE HEAVENLY BODY ***",
          "#     These are best obtained directly from:",
          "#          module modgps02wgs84const (constants), or",
          "#          module modgps06gravity (functions)",
          "#",
          "# *** RADIATION CONSTANTS ***",
          "MPC_STEFAN_BOLTZMANN",
          "MPC_SOLAR_CONSTANT",
          "#",
          "# *** THERMODYNAMIC CONSTANTS (GAS PHASE) ***",
          "MPC_RGAS_IDEAL",
          "MPC_RGAS_DRY_AIR",
          "MPC_RGAS_VAPOUR",
          "MPC_EPS1",
          "MPC_EPS2",
          "MPC_DELTA",
          "#",
          "# *** THERMODYNAMIC CONSTANTS (MOLAR MASSES) ***",
          "MPC_MOLAR_MASS_DRY_AIR",
          "MPC_MOLAR_MASS_VAPOUR",
          "#",
          "# *** THERMODYNAMIC CONSTANTS (SPECIFIC HEATS) ***",
          "MPC_CV_DRY_AIR",
          "MPC_CV_VAPOUR",
          "MPC_CP_DRY_AIR",
          "MPC_CP_VAPOUR",
          "MPC_CP_ICE",
          "MPC_KAPPA",
          "#",
          "# *** THERMODYNAMIC CONSTANTS (ENTHALPIES) ***",
          "MPC_HEAT_CONDENS_WATER",
          "MPC_HEAT_FUSION_WATER",
          "MPC_HEAT_SUBL_WATER",
          "#",
          "# *** THERMODYNAMIC CONSTANTS (FLUID DYNAMICS) ***",
          "MPC_KARMAN" ,
          "MPC_CRITICAL_RICHARDSON",
          "MPC_DENSITY_WATER",
          "MPC_SCHUMANN_NEWELL_LAPSE_RATE",
          "#",
          "# *** THERMODYNAMIC CONSTANTS (OTHER CHARACTERISTICS OF WATER) ***",
          "MPC_T_ICE",
          "MPC_TRIPLE_POINT",
          "#",
          "# *** USED TO CALCULATE L/CP IN FUNC HTVOCP ***",
          "# Consider removing these 7 variables that are 'USED TO CALCULATE L/CP IN FUNC HTVOCP'",
          "MPC_AI",
          "MPC_AW",
          "MPC_BI",
          "MPC_BW",
          "MPC_SLP",
          "MPC_T1S",
          "MPC_T2S",
          "#",
          "# <<<<<< C O N S T A N T S   C O N C E R N I N G   T H E   L I M I T A T I O N S   O F   D I G I T A L   C A L C U L A T I O N >>>>>>",
          "MPC_MINIMUM_HU",
          "MPC_MAXIMUM_ES"
        )
    return(PhysDict, PhysList)


def calculate_values(PhysDict):
    """Some values are calculated, based on values set in this script and on those read from outside."""
    
    (item_type,item_value,item_units) = PhysDict["MPC_RADIANS_PER_DEGREE"]
    item_value = PhysDict["MPC_PI"][VALUE] / 180.
    PhysDict["MPC_RADIANS_PER_DEGREE"] = (item_type,item_value,item_units)

    (item_type,item_value,item_units) = PhysDict["MPC_DEGREES_PER_RADIAN"]
    item_value = 1.0 / PhysDict["MPC_RADIANS_PER_DEGREE"][VALUE]
    PhysDict["MPC_DEGREES_PER_RADIAN"] = (item_type,item_value,item_units)

    (item_type,item_value,item_units) = PhysDict["MPC_KNOTS_PER_M_PER_S"]
    item_value = 1.0 / PhysDict["MPC_M_PER_S_PER_KNOT"][VALUE]
    PhysDict["MPC_KNOTS_PER_M_PER_S"] = (item_type,item_value,item_units)

    (item_type,item_value,item_units) = PhysDict["MPC_RGAS_IDEAL"]
    item_value = PhysDict["MPC_AVOGADRO"][VALUE] * PhysDict["MPC_BOLTZMANN"][VALUE]
    PhysDict["MPC_RGAS_IDEAL"] = (item_type,item_value,item_units)

    (item_type,item_value,item_units) = PhysDict["MPC_CV_DRY_AIR"]
    item_value = PhysDict["MPC_CP_DRY_AIR"][VALUE] - PhysDict["MPC_RGAS_DRY_AIR"][VALUE]
    PhysDict["MPC_CV_DRY_AIR"] = (item_type,item_value,item_units)

    (item_type,item_value,item_units) = PhysDict["MPC_CV_VAPOUR"]
    item_value = PhysDict["MPC_CP_VAPOUR"][VALUE] - PhysDict["MPC_RGAS_VAPOUR"][VALUE]
    PhysDict["MPC_CV_VAPOUR"] = (item_type,item_value,item_units)

    return(PhysDict)


def read_thermoconsts(PhysDict):
    """Read in the content of the file, thermoconsts."""

    stdout.write( ''.join(["\nReading values from the physics library: ", phys_lib_in_file.name, "\n\n"]))

                                        # In this dictionary:
                                        #    key   - the name in thermoconsts
                                        #    value - the corresponding constant/variable name in assimilation
    # It can be noted that 'KARMAN' is ignored from the file, thermoconsts.
    # That is because it is set here to the value used by most groups.
    localNameDict={'CPD'    : 'MPC_CP_DRY_AIR',
                   'CPV'    : 'MPC_CP_VAPOUR',
                   'RGASD'  : 'MPC_RGAS_DRY_AIR',
                   'RGASV'  : 'MPC_RGAS_VAPOUR',
                   'TRPL'   : 'MPC_TRIPLE_POINT',
                   'TCDK'   : 'MPC_K_C_DEGREE_OFFSET',
                   'RAUW'   : 'MPC_DENSITY_WATER',
                   'EPS1'   : 'MPC_EPS1',
                   'EPS2'   : 'MPC_EPS2',
                   'DELTA'  : 'MPC_DELTA',
                   'CAPPA'  : 'MPC_KAPPA',
                   'TGL'    : 'MPC_T_ICE',
                   'CONSOL' : 'MPC_SOLAR_CONSTANT',
                   'GRAV'   : 'NOT_HERE',
                   'RAYT'   : 'NOT_HERE',
                   'STEFAN' : 'MPC_STEFAN_BOLTZMANN',
                   'OMEGA'  : 'NOT_HERE',
                   'KNAMS'  : 'MPC_M_PER_S_PER_KNOT',
                   'STLO'   : 'MPC_SCHUMANN_NEWELL_LAPSE_RATE',
                   'RIC'    : 'MPC_CRITICAL_RICHARDSON',
                   'CHLC'   : 'MPC_HEAT_CONDENS_WATER',
                   'CHLF'   : 'MPC_HEAT_FUSION_WATER',
                   'T1S'    : 'MPC_T1S',
                   'T2S'    : 'MPC_T2S',
                   'AW'     : 'MPC_AW',
                   'BW'     : 'MPC_BW',
                   'AI'     : 'MPC_AI',
                   'BI'     : 'MPC_BI',
                   'SLP'    : 'MPC_SLP'
                  }
    
    while 1:
                                        # Read a line
        s_InLine = phys_lib_in_file.readline()
        if len(s_InLine) == 0: break

        # Split the line into words
        o_InWords = s_InLine.split()
        name_thermoconsts  = o_InWords[0]
        value_thermoconsts = float(o_InWords[1])

        # Massage the format of the value
##                                        # Convert E+ and e+ into D+
##        value_thermoconsts = 'D+'.join(value_thermoconsts.split('E+'))
##        value_thermoconsts = 'D+'.join(value_thermoconsts.split('e+'))

##                                        # Convert E- and e- into E-
##        value_thermoconsts = 'D-'.join(value_thermoconsts.split('E-'))
##        value_thermoconsts = 'D-'.join(value_thermoconsts.split('e-'))

        if name_thermoconsts in localNameDict:
            name_local = localNameDict[name_thermoconsts]
            if name_local in PhysDict:
                item_type  = PhysDict[name_local][0]
                item_units = PhysDict[name_local][2]
                PhysDict[name_local] = (item_type, value_thermoconsts, item_units)

        else:
            stdout.write( ''.join([name_thermoconsts, " read from thermoconsts, but will be ignored.\n"]))

    file_out.write("  \n")

    return(PhysDict)


def write_fortran_MPCmodule(PhysDict, PhysList):
    """Write a fortran file that contains the mathematical and physical constants/variables."""

    stdout.write( ''.join(["\nCreating the file, ", file_out.name, "\n"]))

    #
    # Create the module structure (beginning) around the data
    #
    file_out.write("! This file is considered to be the source of the contained module for the\n")
    file_out.write("! assimilation code.  It can be noted, however, that it was produced with the aid\n")
    file_out.write("! of the script, thermoconsts_to_MathPhysConstants.py.  That script is simply a\n")
    file_out.write("! convenient tool that may or may not be used to update this file.   This is the\n")
    file_out.write("! command that produces this file:\n")
    file_out.write("!          ./thermoconsts_to_MathPhysConstants.py\n")
    file_out.write("\n")
    file_out.write("MODULE MathPhysConstants_mod\n")
    file_out.write("  ! Mathematical and Physical Constants\n")
    file_out.write("  implicit none\n")
    file_out.write("  public\n")

    file_out.write("  !-----------------------------------------------------------------------------\n")
    file_out.write("  ! MODULE MathPhysConstants_mod (prefix='mpc')  \n")
    file_out.write("  !\n")
    file_out.write("  ! Purpose: To supply mathematical and physical constants in a universal,\n")
    file_out.write("  !          reliable fashion.\n")
    file_out.write("  !-----------------------------------------------------------------------------\n")
    file_out.write("  ! Feb  9, 2012 Jeff Blezius\n")


    #
    # Write out the values as real(8) and then real(4) constants / variables
    #
    for real_size in ('8', '4'):
        file_out.write('\n')
        for item_name in PhysList:
            if item_name.find('#') != -1:
                file_out.write(''.join(['!', item_name[1:], '\n']))
            else:
                (item_type,item_value,item_units) = PhysDict[item_name]
                # Determine constant or variable
                if item_type.find('const') != -1:
                    # This value is to be a constant
                    s_parameter = ', parameter'
                else:
                    # This value is to be a variable
                    s_parameter = '           '

                # Prepare the string to be written
                print_string = "   real({size}){parameter} ::{name:34}= {value:22.15E} ! {units}\n".format(size=real_size, parameter=s_parameter, name=''.join([item_name, '_R', real_size]), value=item_value, units=item_units)
                # Write the value in double precision
                print_string = 'D+'.join(print_string.split('E+'))
                print_string = 'D-'.join(print_string.split('E-'))
                file_out.write(print_string)

    #
    # Create the module structure (ending) around the data
    #
    file_out.write("\n")
    file_out.write("contains\n")
    file_out.write("\n")
    file_out.write("  subroutine mpc_setValue(name, var_r4, var_r8, value)\n")
    file_out.write("     ! A means to change a (non-parameter) value\n")
    file_out.write("     character(len=*) :: name\n")
    file_out.write("     real(4) :: var_r4\n")
    file_out.write("     real(8) :: var_r8, value\n")
    file_out.write("\n")
    file_out.write("     var_r4 = value\n")
    file_out.write("     var_r8 = value\n")
    file_out.write("\n")
    file_out.write("     write(6,*)'****************************************************'\n")
    file_out.write("     write(6,*)'*'\n")
    file_out.write("     write(6,*)'*    THE VALUE OF THE CONSTANT ', name\n")
    file_out.write("     write(6,*)'*'\n")
    file_out.write("     write(6,'(A25, D22.15)')' *    HAS BEEN CHANGED TO ', var_r8\n")
    file_out.write("     write(6,*)'*'\n")
    file_out.write("     write(6,*)'****************************************************'\n")
    file_out.write("  end subroutine mpc_setValue\n")
    file_out.write("\n")
    file_out.write("  subroutine mpc_printConstants(kulout)\n")
    file_out.write("     integer, intent(in) :: kulout     ! unit number for printing\n")
    file_out.write("\n")
    file_out.write("     write(kulout,FMT='(//,4x ,\"*** mpc_printConstants: definition of Mathematical and Physical constants  ***\",/)')\n")

    for item_name in PhysList:
        if item_name.find('#') != -1:
            file_out.write(''.join(["     write(kulout,*) \"", item_name[1:], "\"\n"]))
        else:
            (item_type,item_value,item_units) = PhysDict[item_name]

            var_name = ''.join([item_name, '_R8'])
            file_out.write( ''.join(['     write(kulout,"(A36, D22.15, A)") "', var_name.rjust(34) + '= ", ', var_name, '   ', ', "',' '+item_units, '\"\n']))
    file_out.write( '     write(kulout,*) "\\n\\n"')


    file_out.write("\n")
    file_out.write("  end subroutine mpc_printConstants\n")
    file_out.write("end MODULE MathPhysConstants_mod\n")


# Execution starts here

phys_lib_in_file = open('/home/binops/afsi/sio/datafiles/constants/thermoconsts', 'r')
## this path on the 'science' network if
##     /home/smco502/datafiles/constants/thermoconsts
file_out  = open('mathphysconstants_mod.ftn90', 'w+')

(PhysDict, PhysList) = init_phys_values()
PhysDict = read_thermoconsts(PhysDict)
PhysDict = calculate_values(PhysDict)
write_fortran_MPCmodule(PhysDict, PhysList)

file_out.close()
phys_lib_in_file.close()
