#!/bin/bash

set -e

## this script verifies that all the programs have been compiled for a platform
platform=${1:-${ORDENV_PLAT}}
midas_abs=${2}

SEQ_MAESTRO_SHORTCUT=${SEQ_MAESTRO_SHORTCUT:-". ssmuse-sh -d eccc/cmo/isst/maestro/1.5.3.3"}
which getdef 1>/dev/null 2>&1 || ${SEQ_MAESTRO_SHORTCUT}

suite=$(git rev-parse --show-toplevel)/maestro/suites/midas_system_tests
if [ -z "${COMPILING_MACHINE_PPP}" ]; then
    COMPILING_MACHINE_PPP=$(cd ${suite}; getdef --exp ${suite} resources/resources.def FRONTEND)
fi

absdir=${MIDAS_ABS:-${COMPILEDIR_MIDAS_MAIN:-"../../compiledir"}/midas_abs}
revnum=$(git describe --abbrev=7 --always --dirty=_M 2>/dev/null || ssh ${COMPILING_MACHINE_PPP} "cd $PWD; git describe --abbrev=7 --always --dirty=_M" 2>/dev/null || echo unkown revision)

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
