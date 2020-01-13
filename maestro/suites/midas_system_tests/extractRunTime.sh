#!/bin/bash

set -e

SEQ_MAESTRO_SHORTCUT=${SEQ_MAESTRO_SHORTCUT:-". ssmuse-sh -d eccc/cmo/isst/maestro/1.5.3.3"}
which nodehistory 1>/dev/null 2>&1 || ${SEQ_MAESTRO_SHORTCUT}

## This scrits looks in the logs of the maestro suite:
##  - exit with 1 if it finds an 'abortx?' message
##  - exit with 0 if it find the message 'endx'

if [ $# -ne 1 -a $# -ne 2 -a $# -ne 3 ]; then
    echo "extractRunTime.sh: this scripts accepts only one, two or three argument which are the maestro suite, if statistics are computed (default is not) and if we search for outliers." >&2
    exit 1
fi

suite=$1
if [[ "${suite}" != /* ]]; then
    suite=~/.suites/${suite}
fi

computeStats=${2:-no}
if [ "${computeStats}" = yes ]; then
    echo "The statistics are given like this:"
    echo "Mean, Stddev, Mean/Stddev, min, max, Number of cases"
fi

findOutliers=${3:-no}

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
        if [ "${findOutliers}" = yes ]; then
            nodehistory -n ${findRunTime_nodes}/run -history 0 -edate ${logdate} | grep  'The runtime was [.0-9][.0-9]* seconds which is greater than the maximum allowed' || true
        fi
        __findRunTime_runtime__=$(nodehistory -n ${findRunTime_nodes}/run -history 0 -edate ${logdate} | grep 'The runtime was [.0-9][.0-9]* seconds' | sed 's/%/%%/g')
        if [ "${computeStats}" = yes ]; then
            __findRunTime_stats__=$(printf "${__findRunTime_runtime__}" | awk '
BEGIN {
   number=0
   sum=0
   sum2=0
   max=0
   min=10000
}

{
   timing=$6
   sum+=timing
   sum2+=timing**2
   if (timing<min) min=timing
   if (timing>max) max=timing
   number++
}

END {
   mean=sum/number
   var=sum2/number-mean**2
   print mean, sqrt(var), sqrt(var)/mean, min, max, number
}')
            printf "\t${__findRunTime_stats__}\n"
            unset __findRunTime_stats__
        else
            if [ -n "${__findRunTime_runtime__}" ]; then
                echo -e "\t$(printf "${__findRunTime_runtime__}" | head -1)"
            else
                echo -e "\tNo run time was available for that test"
            fi
        fi
        unset __findRunTime_runtime__
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

