#!/bin/bash

set -e

LIBAPPL="f90sqlite udfsqlite rttov_coef_io rttov_hdf rttov_parallel rttov_main rttov_emis_atlas rttov_other ${HDF5_LIBS} burp_module descrip $MPILIB"

