#!/bin/bash

set -e

## this script verifies that all the programs have been compiled for a platform
platform=${1:-${ORDENV_PLAT}}
midas_abs=${2}

compiledir_main=${COMPILEDIR_MIDAS_MAIN:-"../../compiledir"}
absdir=${compiledir_main}/midas_abs
revnum=$(git describe --abbrev=7 --always --dirty=_M 2>/dev/null || ssh eccc-ppp4 "cd $PWD; git describe --abbrev=7 --always --dirty=_M" 2>/dev/null || echo unkown revision)

program_missing=0

[ -n "${midas_abs}" ] && [ ! -d "${midas_abs}" ] && mkdir -p ${midas_abs}

for file in *.f90; do
    program=$(basename ${file} .f90)
    midasAbs=midas-${program}_${platform}-${revnum}.Abs
    if [ -f "${absdir}/${midasAbs}" ]; then
        if [ -n "${midas_abs}" ]; then
            echo "The program '${absdir}/${midasAbs}' exists and is copied to '${midas_abs}'"
            cp ${absdir}/${midasAbs} ${midas_abs}
        else
            echo "The program '${absdir}/${midasAbs}' exists"
        fi
    else
        echo "The program '${absdir}/${midasAbs}' does not exist" >&2
        program_missing=1
    fi
done

if [ "${program_missing}" -ne 0 ]; then
    echo "Some program was missing... Abort!" >&2
    exit 1
fi
