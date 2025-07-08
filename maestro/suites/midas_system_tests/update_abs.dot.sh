#!/bin/bash

set -euo pipefail

typeset -r toplevel=${1:-$(git rev-parse --show-toplevel)}

cat > ${TMPDIR}/abs.dot <<EOF
MIDAS_toplevel=${toplevel}
ABS_DIR=${MIDAS_COMPILE_DIR_MAIN:-${toplevel}/compiledir}/midas_abs >> abs.dot
MIDAS_version=\$(cd ${toplevel}; ./midas.version)
EOF

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
