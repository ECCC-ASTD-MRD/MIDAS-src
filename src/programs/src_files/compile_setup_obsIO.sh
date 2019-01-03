#!/bin/bash

set -e

LIBAPPL="burp_module descrip $MPILIB f90sqlite udfsqlite"

COMPF="${COMPF} -defines =-DCODEPRECISION_OBS_REAL_SINGLE"
