#!/bin/bash

set -e

SEQ_MAESTRO_SHORTCUT=${SEQ_MAESTRO_SHORTCUT:-". ssmuse-sh -d eccc/cmo/isst/maestro/1.5.3"}
which nodehistory 1>/dev/null 2>&1 || ${SEQ_MAESTRO_SHORTCUT}

## This scrits looks in the logs of the maestro suite:
##  - exit with 1 if it finds an 'abortx?' message
##  - exit with 0 if it find the message 'endx'

if [ $# -ne 1 ]; then
    echo "extractRunTime.sh: this scripts accepts only and only one argument which is the maestro suite." >&2
    exit 1
fi

suite=$1
if [[ "${suite}" != /* ]]; then
    suite=~/.suites/${suite}
fi

if [ ! -d "${suite}" ]; then
    echo "The suite given '${suite}' does not exist." >&2
    exit 1
fi

echo "Extract run times for tests in ${suite}"

logs=$(/bin/ls -t ${suite}/logs/*_nodelog 2>/dev/null | head -1)
if [ -z "${logs}" ]; then
    echo "The suite '${suite}' does not contain any logs in '${suite}/logs'." >&2
    exit 1
fi

logdate=$(basename ${logs} _nodelog)

findRunTime () {
    set -e

    findRunTime_node=$1
    findRunTime_nodes="$(nodeinfo -n ${findRunTime_node} | grep '^node\.submit=' | cut -d= -f2)"
    if [[ "${findRunTime_nodes}" = /*/UnitTest ]]; then
        echo ${findRunTime_nodes%/*}
        __findRunTime_line__=$(nodehistory -n ${findRunTime_nodes}/run -history 0 -edate ${logdate} | grep 'The runtime was [.0-9][.0-9]* seconds' | head -1)
        if [ -n "${__findRunTime_line__}" ]; then
            echo -e "\t${__findRunTime_line__}"
        else
            echo -e "\tNo run time was available for that test"
        fi
        unset __findRunTime_line__
    else
        for __node__ in ${findRunTime_nodes}; do
            findRunTime ${__node__}
        done
        unset __node__
    fi

    unset findRunTime_node findRunTime_nodes
}  ## End of function 'findRunTime'

export SEQ_EXP_HOME=${suite}
findRunTime /Tests

