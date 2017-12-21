#!/bin/bash

set -e

LIBAPPL="netcdff rttov_coef_io rttov_hdf rttov_parallel rttov_main rttov_emis_atlas rttov_other ${HDF5_LIBS} burp_module descrip $MPILIB"

SRC_FILES="rttov_interfaces_mod.ftn90 codeprecision_mod.f90 utilities_mod.f90"
SRC_FILES="$SRC_FILES ramdisk_mod.f90 mpi_mod.f90 mpivar_mod.f90 bufr_mod.f90 codtyp_mod.f90"
SRC_FILES="$SRC_FILES mathphysconstants_mod.f90 earthconstants_mod.f90"
SRC_FILES="$SRC_FILES randomnumber_mod.f90"
SRC_FILES="$SRC_FILES controlvector_mod.f90"
SRC_FILES="$SRC_FILES physicsfunctions_mod.f90"
SRC_FILES="$SRC_FILES varnamelist_mod.f90"
SRC_FILES="$SRC_FILES obsspacedata_mod.f90"
SRC_FILES="$SRC_FILES presprofileoperators_mod.f90 tovs_extrap_mod.f90"
SRC_FILES="$SRC_FILES horizontalcoord_mod.f90 verticalcoord_mod.f90 timecoord_mod.f90"
SRC_FILES="$SRC_FILES rmatrix_mod.f90 hirchannels_mod.f90 columndata_mod.f90 ozoneclim_mod.f90 tovs_nl_mod.f90"
SRC_FILES="$SRC_FILES analysisgrid_mod.f90"
SRC_FILES="$SRC_FILES gridstatevector_mod.f90"
SRC_FILES="$SRC_FILES lamspectraltransform_mod.f90 lambmatrixhi_mod.f90"
SRC_FILES="$SRC_FILES globalspectraltransform_mod.f90 bmatrixhi_mod.f90"
SRC_FILES="$SRC_FILES localizationfunction_mod.f90 localizationspectral_mod.f90 localization_mod.f90"
SRC_FILES="$SRC_FILES variabletransforms_mod.f90"
SRC_FILES="$SRC_FILES spectralfilter_mod.f90 bmatrixlatbands_mod.f90"
SRC_FILES="$SRC_FILES ensemblestatevector_mod.f90 bmatrixensemble_mod.f90"
SRC_FILES="$SRC_FILES diffusion_mod.f90 bmatrixdiff_mod.f90"
SRC_FILES="$SRC_FILES bmatrixchem_mod.f90 biascorrection_mod.f90 bmatrix_mod.f90"
