#!/bin/bash

set -e

LIBAPPL="descrip $MPILIB"

SRC_FILES="codeprecision_mod.f90 clib_interfaces_mod.ftn90 randomnumber_mod.f90 utilities_mod.f90 ramdisk_mod.f90 bufr_mod.f90"
SRC_FILES="$SRC_FILES mathphysconstants_mod.f90 earthconstants_mod.f90 mpi_mod.f90 mpivar_mod.f90 codtyp_mod.f90"
SRC_FILES="$SRC_FILES physicsfunctions_mod.f90 horizontalcoord_mod.f90 obsspacedata_mod.f90 analysisgrid_mod.f90"
SRC_FILES="$SRC_FILES filenames_mod.f90"
SRC_FILES="$SRC_FILES localizationfunction_mod.f90 controlvector_mod.f90 varnamelist_mod.f90 timecoord_mod.f90"
SRC_FILES="$SRC_FILES globalspectraltransform_mod.f90 diffusion_mod.f90 verticalcoord_mod.f90"
SRC_FILES="$SRC_FILES lamspectraltransform_mod.f90 gridstatevector_mod.f90 humiditylimits_mod.f90 ensemblestatevector_mod.f90"
SRC_FILES="$SRC_FILES spectralfilter_mod.f90 variabletransforms_mod.f90 lambmatrixhi_mod.f90"
SRC_FILES="$SRC_FILES localizationspectral_mod.f90 localization_mod.f90 advection_mod.f90"
SRC_FILES="$SRC_FILES bmatrixensemble_mod.f90 bmatrixhi_mod.f90 bmatrixlatbands_mod.f90 bmatrixdiff_mod.f90 bmatrixchem_mod.f90"
SRC_FILES="$SRC_FILES bmatrix_mod.f90"
