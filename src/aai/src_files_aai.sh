#!/bin/bash

set -e

LIBAPPL="descrip $MPILIB"

SRC_FILES="utilities_mod.f90 bufr_mod.f90 ramdisk_mod.ftn90"
SRC_FILES="$SRC_FILES mathphysconstants_mod.f90 earthconstants_mod.f90 mpi_mod.f90 mpivar_mod.f90"
SRC_FILES="$SRC_FILES physicsfunctions_mod.f90 horizontalcoord_mod.f90 analysisgrid_mod.f90"
SRC_FILES="$SRC_FILES varnamelist_mod.f90 obsspacedata_mod.ftn90 timecoord_mod.f90"
SRC_FILES="$SRC_FILES verticalcoord_mod.f90"
SRC_FILES="$SRC_FILES gridstatevector_mod.f90 humiditylimits_mod.f90"
