#!/bin/bash

set -e

if [ "${ORDENV_PLAT}" = sles-11-amd64-64 -o "${ORDENV_PLAT}" = sles-11-broadwell-64-xc40 ]; then
    echo "... loading cray-netcdf"
    module load cray-netcdf
fi

LIBAPPL="netcdff burp_module descrip $MPILIB"

