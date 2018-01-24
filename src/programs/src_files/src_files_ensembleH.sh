#!/bin/bash

set -e

LIBAPPL="rttov_coef_io rttov_hdf rttov_parallel rttov_main rttov_emis_atlas rttov_other ${HDF5_LIBS} burp_module descrip $MPILIB"

SRC_FILES="rttov_interfaces_mod.ftn90 codeprecision_mod.f90 utilities_mod.f90 bufr_mod.f90 ramdisk_mod.f90 filenames_mod.f90"
SRC_FILES="$SRC_FILES mathphysconstants_mod.f90 earthconstants_mod.f90 mpi_mod.f90 mpivar_mod.f90 codtyp_mod.f90 randomnumber_mod.f90"
SRC_FILES="$SRC_FILES physicsfunctions_mod.f90 horizontalcoord_mod.f90 analysisgrid_mod.f90"
SRC_FILES="$SRC_FILES varnamelist_mod.f90 obsspacedata_mod.f90 timecoord_mod.f90"
SRC_FILES="$SRC_FILES lqtoes_mod.f90 presprofileoperators_mod.f90"
SRC_FILES="$SRC_FILES tovs_extrap_mod.f90 rmatrix_mod.f90 hirchannels_mod.f90 "
SRC_FILES="$SRC_FILES verticalcoord_mod.f90 columndata_mod.f90 ozoneclim_mod.f90 tovs_nl_mod.f90"
SRC_FILES="$SRC_FILES gridstatevector_mod.f90"
SRC_FILES="$SRC_FILES lamspectraltransform_mod.f90 globalspectraltransform_mod.f90"
SRC_FILES="$SRC_FILES variabletransforms_mod.f90"
SRC_FILES="$SRC_FILES burpread_mod.f90 obssubspacedata_mod.f90"
SRC_FILES="$SRC_FILES gps_mod.f90 tt2phi_mod.f90 windrotation_mod.f90 statetocolumn_mod.f90"
SRC_FILES="$SRC_FILES burpfiles_mod.f90 multi_ir_bgck_mod.f90 chem_setup_mod.f90 chem_obserrors_mod.f90"
SRC_FILES="$SRC_FILES bmatrixchem_mod.f90 chem_obsoperators_mod.f90 chem_postproc_mod.f90"
SRC_FILES="$SRC_FILES obserrors_mod.f90 varqc_mod.f90 obsfilter_mod.f90 tovs_lin_mod.f90 obsoperators_mod.f90"
SRC_FILES="$SRC_FILES innovation_mod.f90 enkf_mod.f90"