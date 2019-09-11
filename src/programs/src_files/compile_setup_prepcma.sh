#!/bin/bash

set -e

LIBAPPL="rttov_coef_io rttov_hdf rttov_parallel  rttov_main rttov_emis_atlas rttov_other ${HDF5_LIBS} f90sqlite udfsqlite burp_module vgrid irc $MPILIB"

COMPF="${COMPF} -DCODEPRECISION_OBS_REAL_SINGLE"
