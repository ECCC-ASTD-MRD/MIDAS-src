#!/bin/bash

set -e

LIBAPPL="rttov_coef_io rttov_hdf rttov_parallel  rttov_main rttov_emis_atlas rttov_other ${HDF5_LIBS} burp_module descrip $MPILIB"
SRC_FILES="utilities_mod.f90 bufr_mod.f90 ramdisk_mod.ftn90"
SRC_FILES="$SRC_FILES mathphysconstants_mod.f90 earthconstants_mod.f90 mpi_mod.f90 mpivar_mod.f90"
SRC_FILES="$SRC_FILES physicsfunctions_mod.f90 horizontalcoord_mod.f90 analysisgrid_mod.f90"
SRC_FILES="$SRC_FILES varnamelist_mod.f90 obsspacedata_mod.ftn90 timecoord_mod.f90"
SRC_FILES="$SRC_FILES verticalcoord_mod.f90"
SRC_FILES="$SRC_FILES gridstatevector_mod.f90 ensemblestatevector_mod.f90"
SRC_FILES="$SRC_FILES lamspectraltransform_mod.f90 globalspectraltransform_mod.f90"
SRC_FILES="$SRC_FILES variabletransforms_mod.f90"
