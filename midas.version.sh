#!/bin/bash

set -e

if [ $# -gt 1 ]; then
    echo "$0: Cannot give more than 1 arguments.  You gave '${*}'." >&2
    exit 1
fi

abbrev=7
if [ $# -eq 1 ]; then
    abbrev=${1}
fi

command="git describe --abbrev=${abbrev} --always --dirty=_M"

status=0
${command} || status=1

if [ "${status}" -ne 0 ]; then
    toplevel=$(git rev-parse --show-toplevel)
    ${toplevel}/set_resources_def.sh
    suite=${toplevel}/maestro/suites/midas_system_tests
    . ${suite}/set_machine_list.dot ${suite}

    ssh ${MACHINE_PPP} "cd ${toplevel}; ${command}"  || echo unkown revision
fi
