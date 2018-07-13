#!/bin/bash

set -e

LIBAPPL="burp_module descrip $MPILIB f90sqlite udfsqlite"

SRC_FILES="utilities_mod.f90 ramdisk_mod.f90 mpi_mod.f90 mpivar_mod.f90 timecoord_mod.f90"
SRC_FILES="$SRC_FILES mathphysconstants_mod.f90 earthconstants_mod.f90 bufr_mod.f90 codtyp_mod.f90"
SRC_FILES="$SRC_FILES codeprecision_mod.f90 obsspacedata_mod.f90 obssubspacedata_mod.f90 obsUtil_mod.f90"
SRC_FILES="$SRC_FILES sqlite_read_mod.f90 burpread_mod.f90"
SRC_FILES="$SRC_FILES sqlitefiles_mod.f90 burpfiles_mod.f90 cmafiles_mod.f90 obsfiles_mod.f90"
