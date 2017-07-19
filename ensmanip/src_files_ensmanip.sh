#!/bin/bash

set -e

SRC_FILES="utilities_mod.ftn90 toplevelcontrol_mod.ftn90 bufr_mod.ftn90"
SRC_FILES="$SRC_FILES mathphysconstants_mod.ftn90 earthconstants_mod.ftn90 mpi_mod.ftn90 mpivar_mod.ftn90"
SRC_FILES="$SRC_FILES physicsfunctions_mod.ftn90 horizontalcoord_mod.ftn90 analysisgrid_mod.ftn90"
SRC_FILES="$SRC_FILES varnamelist_mod.ftn90 obsspacedata_mod.ftn90 timecoord_mod.ftn90"
SRC_FILES="$SRC_FILES verticalcoord_mod.ftn90"
SRC_FILES="$SRC_FILES gridstatevector_mod.ftn90 ensemblestatevector_mod.ftn90"
SRC_FILES="$SRC_FILES lamspectraltransform_mod.ftn90"
SRC_FILES="$SRC_FILES variabletransforms_mod.ftn90"
