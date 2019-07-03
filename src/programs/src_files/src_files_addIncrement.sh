#!/bin/bash

set -e

LIBAPPL="f90sqlite udfsqlite burp_module vgrid $MPILIB"

SRC_FILES="clib_interfaces_mod.ftn90 codeprecision_mod.ftn90 utilities_mod.f90 bufr_mod.f90 ramdisk_mod.f90 randomnumber_mod.f90"
SRC_FILES="$SRC_FILES mathphysconstants_mod.f90 earthconstants_mod.f90 mpi_mod.f90 mpivar_mod.f90 codtyp_mod.f90"
SRC_FILES="$SRC_FILES physicsfunctions_mod.f90 horizontalcoord_mod.f90 analysisgrid_mod.f90 windrotation_mod.f90 "
SRC_FILES="$SRC_FILES obsspacedata_mod.f90 obssubspacedata_mod.f90 obsUtil_mod.f90 burpread_mod.f90 burpfiles_mod.f90 sqlite_read_mod.f90 sqlitefiles_mod.f90 cmafiles_mod.f90 obsfiles_mod.f90"
SRC_FILES="$SRC_FILES timecoord_mod.f90 obstimeinterp_mod.f90 varnamelist_mod.f90 diffusion_mod.f90 filenames_mod.f90"
SRC_FILES="$SRC_FILES verticalcoord_mod.f90 columndata_mod.f90 gridstatevector_mod.f90"
SRC_FILES="$SRC_FILES globalspectraltransform_mod.f90 lamspectraltransform_mod.f90 variabletransforms_mod.f90 controlvector_mod.f90"
SRC_FILES="$SRC_FILES humiditylimits_mod.f90 ensemblestatevector_mod.f90 advection_mod.f90 localizationfunction_mod.f90 localizationspectral_mod.f90 localization_mod.f90 spectralfilter_mod.f90"
SRC_FILES="$SRC_FILES lambmatrixhi_mod.f90 bmatrixdiff_mod.f90 bmatrixlatbands_mod.f90 bmatrixhi_mod.f90 bmatrixensemble_mod.f90 bmatrixchem_mod.f90 bmatrix_mod.f90"
SRC_FILES="$SRC_FILES tt2phi_mod.f90 statetocolumn_mod.f90 chem_setup_mod.f90 chem_postproc_mod.f90 increment_mod.f90"
