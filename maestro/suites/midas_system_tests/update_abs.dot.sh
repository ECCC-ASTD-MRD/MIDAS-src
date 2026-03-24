#!/bin/bash

set -euo pipefail

typeset -r toplevel=${1:-$(git rev-parse --show-toplevel)}

echo "MIDAS_toplevel=${toplevel}" > ${TMPDIR}/abs.dot
MIDAS_ABS_LEAFDIR=${MIDAS_ABS_LEAFDIR:-midas_abs}
if [ "${MIDAS_COMPILE_DIR_MAIN:-}" = build_directory_local_to_the_repository -o -z "${MIDAS_COMPILE_DIR_MAIN:-}" ]; then
    echo "ABS_DIR=${toplevel}/compiledir/${MIDAS_ABS_LEAFDIR}" >> ${TMPDIR}/abs.dot
else
    echo "ABS_DIR=${MIDAS_COMPILE_DIR_MAIN}/${MIDAS_ABS_LEAFDIR}" >> ${TMPDIR}/abs.dot
fi
echo "MIDAS_version=\$(cd ${toplevel}; ./midas.version)" >> ${TMPDIR}/abs.dot

typeset -r absdotfile=${toplevel}/maestro/suites/midas_system_tests/abs.dot

if [ -f "${absdotfile}" ]; then
    status_diff=0
    cmp --silent ${absdotfile} ${TMPDIR}/abs.dot || status_diff=1

    if [ "${status_diff}" -ne 0 ]; then
        echo "Updating '${absdotfile}'"
    else
        rm ${TMPDIR}/abs.dot
        exit
    fi
else
    echo "Creating '${absdotfile}'"
fi

mv -f ${TMPDIR}/abs.dot ${absdotfile}
