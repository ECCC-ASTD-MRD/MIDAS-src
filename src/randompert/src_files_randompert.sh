#!/bin/bash

set -e

LIBAPPL="descrip $MPILIB"

SRC_FILES="randomnumber_mod.f90 utilities_mod.f90 ramdisk_mod.ftn90 bufr_mod.f90"
SRC_FILES="$SRC_FILES mathphysconstants_mod.f90 earthconstants_mod.f90 mpi_mod.f90 mpivar_mod.f90"
SRC_FILES="$SRC_FILES physicsfunctions_mod.f90 horizontalcoord_mod.f90 obsspacedata_mod.ftn90 analysisgrid_mod.f90"
SRC_FILES="$SRC_FILES localizationfunction_mod.f90 controlvector_mod.f90 varnamelist_mod.f90 timecoord_mod.f90"
SRC_FILES="$SRC_FILES globalspectraltransform_mod.f90 verticalcoord_mod.f90"
SRC_FILES="$SRC_FILES lamspectraltransform_mod.f90 gridstatevector_mod.f90 ensemblestatevector_mod.f90"
SRC_FILES="$SRC_FILES columndata_mod.f90 spectralfilter_mod.f90 lambmatrixhi_mod.f90"
SRC_FILES="$SRC_FILES localizationspectral_mod.f90 localization_mod.f90 variabletransforms_mod.f90 bmatrixensemble_mod.f90 bmatrixhi_mod.f90 bmatrixchem_mod.f90"
SRC_FILES="$SRC_FILES bmatrix_mod.f90"
