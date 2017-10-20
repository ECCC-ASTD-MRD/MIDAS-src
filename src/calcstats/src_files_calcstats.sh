#!/bin/bash

set -e

LIBAPPL="descrip $MPILIB"

SRC_FILES="randomnumber_mod.f90 utilities_mod.f90 bufr_mod.f90"
SRC_FILES="$SRC_FILES mathphysconstants_mod.f90 earthconstants_mod.f90 mpi_mod.f90 mpivar_mod.f90"
SRC_FILES="$SRC_FILES physicsfunctions_mod.f90 horizontalcoord_mod.f90 obsspacedata_mod.ftn90"
SRC_FILES="$SRC_FILES controlvector_mod.f90 varnamelist_mod.f90 timecoord_mod.f90 localizationfunction_mod.f90"
SRC_FILES="$SRC_FILES globalspectraltransform_mod.f90 verticalcoord_mod.f90"
SRC_FILES="$SRC_FILES lamspectraltransform_mod.f90 gridstatevector_mod.f90 spectralfilter_mod.f90"
SRC_FILES="$SRC_FILES columndata_mod.f90 analysisgrid_mod.f90 lambmatrixhi_mod.f90"
SRC_FILES="$SRC_FILES menetrierdiag_mod.f90 calcstatsglb_mod.f90 calcstatslam_mod.f90"
