#!/bin/bash

set -e

LIBAPPL="f90sqlite udfsqlite burp_module descrip $MPILIB"

COMPF="${COMPF} -DCODEPRECISION_OBS_REAL_SINGLE"
