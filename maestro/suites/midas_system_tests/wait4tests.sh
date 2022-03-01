#!/bin/bash

set -e

SEQ_MAESTRO_SHORTCUT=${SEQ_MAESTRO_SHORTCUT:-". ssmuse-sh -d eccc/cmo/isst/maestro/1.8.0"}
which nodehistory 1>/dev/null 2>&1 || ${SEQ_MAESTRO_SHORTCUT}

## This scrits looks in the logs of the maestro suite:
##  - exit with 1 if it finds an 'abortx?' message
##  - exit with 0 if it find the message 'endx'

if [ $# -ne 1 ]; then
    echo "wait4tests.sh: this scripts accepts only and only one argument which is the maestro suite." >&2
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

echo "Waiting for suite ${suite}"
echo "To monitor the suite:"
printf "\tSEQ_EXP_HOME=${suite} xflow\n"

logs=$(/bin/ls -t ${suite}/logs/*_nodelog 2>/dev/null | head -1)
if [ -z "${logs}" ]; then
    echo "The suite '${suite}' does not contain any logs in '${suite}/logs'." >&2
    exit 1
fi

logdate=$(basename ${logs} _nodelog)

findTestsWhoAborted () {
    set -e

    local __findTestsWhoAborted_abort__=0

    __findTestsWhoAborted_known_aborts__=${__findTestsWhoAborted_known_aborts__}

    findTestsWhoAborted_node=$1
    findTestsWhoAborted_nodes="$(nodeinfo -n ${findTestsWhoAborted_node} | grep '^node\.submit=' | cut -d= -f2)"
    if [[ "${findTestsWhoAborted_nodes}" = /*/UnitTest ]]; then
        __abort__=1
        nodehistory -n ${findTestsWhoAborted_nodes} -history 0 -edate ${logdate} | grep -E '^TIMESTAMP=.* MESSAGE=abortx?' 1>/dev/null || __abort__=0
        if [ "${__abort__}" -ne 0 ]; then
            findJobWhoAborted="get $(nodeinfo -n ${findTestsWhoAborted_nodes}/get | grep '^node\.sibling=' | cut -d= -f2)"
            for __node__ in ${findJobWhoAborted}; do
                __abortn__=1
                nodehistory -n ${findTestsWhoAborted_nodes}/${__node__} -history 0 -edate ${logdate} | grep -v '^[[:space:]]*$' | tail -1 | grep -E '^TIMESTAMP=.* MESSAGE=abortx?' 1>/dev/null || __abortn__=0
                if [ "${__abortn__}" -ne 0 ]; then
                    __known_aborted_node_found__=no
                    for __known_aborted_node__ in ${__findTestsWhoAborted_known_aborts__}; do
                        if [ "${findTestsWhoAborted_nodes}/${__node__}" = "${__known_aborted_node__}" ]; then
                            __known_aborted_node_found__=yes
                            break
                        fi
                    done
                    unset __known_aborted_node__
                    if [ "${__known_aborted_node_found__}" != yes ]; then
                        echo "The test '${findTestsWhoAborted_node}' aborted"
                        printf "\tSEQ_EXP_HOME=${suite} nodelister -n ${findTestsWhoAborted_nodes}/${__node__} -type abort\n"
                        if [ -n "${__findTestsWhoAborted_known_aborts__}" ]; then
                            __findTestsWhoAborted_known_aborts__="${__findTestsWhoAborted_known_aborts__} ${findTestsWhoAborted_nodes}/${__node__}"
                        else
                            __findTestsWhoAborted_known_aborts__=${findTestsWhoAborted_nodes}/${__node__}
                        fi
                        __findTestsWhoAborted_abort__=1
                    fi
                    unset __known_aborted_node_found__
                fi
                unset __abortn__
            done
            unset __node__
        fi
        unset __abort__
    else
        for __node__ in ${findTestsWhoAborted_nodes}; do
            findTestsWhoAborted ${__node__} || return 1
        done
        unset __node__
    fi

    unset findTestsWhoAborted_node findTestsWhoAborted_nodes
    if [ "${__findTestsWhoAborted_abort__}" -ne 0 ]; then
        #unset __findTestsWhoAborted_abort__
        return 1
    else
        #unset __findTestsWhoAborted_abort__
        return 0
    fi
}

export SEQ_EXP_HOME=${suite}
while true; do
    status=0
    #nodehistory -n /Tests -history 0 -edate ${logdate} | grep -E '^TIMESTAMP=.* MESSAGE=(abortx?|endx)' 1>/dev/null || status=1
    nodehistory -n /Tests -history 0 -edate ${logdate} | grep -E '^TIMESTAMP=.* MESSAGE=endx' 1>/dev/null || status=1
    if [ "${status}" -eq 0 ]; then
        break
    fi
    abort=0
    findTestsWhoAborted /Tests || abort=1
    sleep ${WAIT4TESTS_SLEEP:-60}
done

abort=1
nodehistory -n /Tests -history 0 -edate ${logdate} | grep -v '^[[:space:]]*$' | tail -1 | grep -E '^TIMESTAMP=.* MESSAGE=abortx?' 1>/dev/null || abort=0
if [ "${abort}" -ne 0 ]; then
    echo "An abort has been found in ${logdate} in suite ${suite}"
    exit 1
fi

echo "The tests has run normally according to ${logdate}"

