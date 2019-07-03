#!/bin/bash

set -e

if [ "${ORDENV_PLAT}" = sles-11-amd64-64 -o "${ORDENV_PLAT}" = sles-11-broadwell-64-xc40 ]; then
    echo "... loading cray-netcdf"
    module load cray-netcdf
fi

LIBAPPL="netcdff burp_module vgrid $MPILIB"

SRC_FILES="codeprecision_mod.ftn90 clib_interfaces_mod.ftn90 utilities_mod.f90"
SRC_FILES="$SRC_FILES ramdisk_mod.f90 filenames_mod.f90 mpi_mod.f90 mpivar_mod.f90 bufr_mod.f90 codtyp_mod.f90"
SRC_FILES="$SRC_FILES mathphysconstants_mod.f90 earthconstants_mod.f90"
SRC_FILES="$SRC_FILES randomnumber_mod.f90"
SRC_FILES="$SRC_FILES controlvector_mod.f90"
SRC_FILES="$SRC_FILES physicsfunctions_mod.f90"
SRC_FILES="$SRC_FILES varnamelist_mod.f90"
SRC_FILES="$SRC_FILES obsspacedata_mod.f90"
SRC_FILES="$SRC_FILES horizontalcoord_mod.f90 verticalcoord_mod.f90 timecoord_mod.f90"
SRC_FILES="$SRC_FILES analysisgrid_mod.f90"
SRC_FILES="$SRC_FILES gridstatevector_mod.f90 humiditylimits_mod.f90 ensemblestatevector_mod.f90"
SRC_FILES="$SRC_FILES globalspectraltransform_mod.f90 lamspectraltransform_mod.f90 variabletransforms_mod.f90 lambmatrixhi_mod.f90"
SRC_FILES="$SRC_FILES bmatrixhi_mod.f90"
SRC_FILES="$SRC_FILES localizationfunction_mod.f90 localizationspectral_mod.f90 localization_mod.f90"
SRC_FILES="$SRC_FILES spectralfilter_mod.f90 bmatrixlatbands_mod.f90"
SRC_FILES="$SRC_FILES advection_mod.f90 bmatrixensemble_mod.f90"
SRC_FILES="$SRC_FILES diffusion_mod.f90 bmatrixdiff_mod.f90"
SRC_FILES="$SRC_FILES bmatrixchem_mod.f90 bmatrix_mod.f90"
