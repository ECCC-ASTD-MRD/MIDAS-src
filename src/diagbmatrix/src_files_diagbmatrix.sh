#!/bin/bash

set -e

LIBAPPL="descrip $MPILIB"

SRC_FILES="randomnumber_mod.ftn90 utilities_mod.ftn90 ramdisk_mod.ftn90 bufr_mod.ftn90"
SRC_FILES="$SRC_FILES mathphysconstants_mod.ftn90 earthconstants_mod.ftn90 mpi_mod.ftn90 mpivar_mod.ftn90"
SRC_FILES="$SRC_FILES  physicsfunctions_mod.ftn90 horizontalcoord_mod.ftn90 obsspacedata_mod.ftn90 analysisgrid_mod.ftn90"
SRC_FILES="$SRC_FILES localizationfunction_mod.ftn90 controlvector_mod.ftn90 varnamelist_mod.ftn90 timecoord_mod.ftn90"
SRC_FILES="$SRC_FILES globalspectraltransform_mod.ftn90 verticalcoord_mod.ftn90"
SRC_FILES="$SRC_FILES lamspectraltransform_mod.ftn90 gridstatevector_mod.ftn90 ensemblestatevector_mod.ftn90"
SRC_FILES="$SRC_FILES columndata_mod.ftn90 spectralfilter_mod.ftn90 lambmatrixhi_mod.ftn90"
SRC_FILES="$SRC_FILES localizationspectral_mod.ftn90 localization_mod.ftn90"
SRC_FILES="$SRC_FILES variabletransforms_mod.ftn90 bmatrixensemble_mod.ftn90 bmatrixhi_mod.ftn90 bmatrixchem_mod.ftn90"
SRC_FILES="$SRC_FILES bmatrix_mod.ftn90"
