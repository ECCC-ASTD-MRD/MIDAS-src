#!/bin/bash

set -e

SEQ_MAESTRO_SHORTCUT=${SEQ_MAESTRO_SHORTCUT:-". ssmuse-sh -d eccc/cmo/isst/maestro/1.8.0-beta"}
which nodehistory 1>/dev/null 2>&1 || ${SEQ_MAESTRO_SHORTCUT}

eval $(cclargs_lite -D ' ' $0 "[ extract timings of 'run' task from MIDAS test suite ]" \
 -suite "" "" "[ suite to extract timings from (default to the same suite as this script is to) ]" \
 -all "no" "yes" "[ extract all the timings ]" \
 -computeStats "no" "yes" "[ compute mean and other statistics from several execution of the tests for the same maestro date (default if 'no')]" \
 -findOutliers "no" "yes" "[ check for outliers in the listings if 'notify' then send an email to a list of emails when outliers are found (see '-emails' argument) (default is 'no') ]" \
 -emails "" "" "[ List of emails to send a message when we find outliers for the execution time ]" \
 -log    "latest" "latest" "[ log date to use to find timings (default: latest log file) ]" \
 ++ $*)

## Note on '-log' option:
##    This option allows the user to choose which log to look.  By default, it
##    is the latest one.  It must a date or a part of a date.  If several log
##    files qualify for the value given, then it takes the latest one.

if [ -z "${suite}" ]; then
    suite=$(dirname $(true_path $0))
elif [[ "${suite}" != /* ]]; then
    suite=${HOME}/.suites/${suite}
fi

if [ ! -d "${suite}" ]; then
    echo "The suite given '${suite}' does not exist." >&2
    exit 1
fi

if [ "${all}" = yes -o "${all}" = no ]; then
    ## transform 'all' into 'extractAll' for a more significant variable in the script
    extractAll=${all}
else
    echo "$0: The '-all' argument myst be 'yes' or 'no' and not '${all}'" >&2
    exit 1
fi

if [ "${computeStats}" != no -a "${computeStats}" != yes ]; then
    echo "$0: The '-computeStats' argument must be 'yes' or 'no' and not '${computeStats}" >&2
    exit 1
fi

if [ "${findOutliers}" != no -a "${findOutliers}" != yes ]; then
    echo "$0: The '-findOutliers' argument must be 'yes' or 'no' and not '${findOutliers}" >&2
    exit 1
fi

if [ "${findOutliers}" = no -a -n "${emails}" ]; then
    echo "$0: WARNING: Since '-findOutliers' argument is 'no', no email will be sent to the list you gave." >&2
fi

echo "Extract run times for tests in ${suite}"

if [ "${log}" = latest ]; then
    logs=$(/bin/ls -t ${suite}/logs/*_nodelog 2>/dev/null | head -1)
    if [ -z "${logs}" ]; then
        echo "The suite '${suite}' does not contain any logs in '${suite}/logs'." >&2
        exit 1
    fi
else
    logs=$(/bin/ls -t ${suite}/logs/*${log}*_nodelog 2>/dev/null | head -1)
    if [ -z "${logs}" ]; then
        echo "The suite '${suite}' does not contain any logs matching '${suite}/logs/*${log}*_nodelog'." >&2
        exit 1
    fi
fi

if [ ! -f "${logs}" ]; then
    echo "Cannot access log file '${logs}'." >&2
    exit 1
fi

logdate=$(basename ${logs} _nodelog)
echo "Using the logdate ${logdate}"

if [ "${computeStats}" = yes ]; then
    echo
    echo "The statistics are given like this:"
    echo -e "\tMean, Stddev, Stddev/Mean, Min, Max, Number of cases"
    echo
fi

findRunTime () {
    set -e

    findRunTime_node=$1
    findRunTime_nodes="$(nodeinfo -n ${findRunTime_node} | grep '^node\.submit=' | cut -d= -f2)"
    if [[ "${findRunTime_nodes}" = /*/UnitTest ]]; then
        echo ${findRunTime_nodes%/*}
        __findRunTime_runtime__=$(nodehistory -n ${findRunTime_nodes}/run -history 0 -edate ${logdate} | grep 'The runtime was [.0-9][.0-9]* seconds' | sed 's/%/%%/g')
        if [ "${computeStats}" = yes ]; then
            __findRunTime_stats__=$(printf "${__findRunTime_runtime__}" | awk '
BEGIN {
   min=10000
   max=0
   sum=0
   sum2=0
   number=0
}

{
   match($0, /The runtime was ([.0-9]+) seconds/, array_timing)
   timing=array_timing[1]
   sum+=timing
   sum2+=timing**2
   if (timing<min) min=timing
   if (timing>max) max=timing
   number++
}

END {
    if (number>0) {
        mean=sum/number
        var=sum2/number-mean**2
        print mean, sqrt(var), sqrt(var)/mean, min, max, number
    }
}')
            printf "\t${__findRunTime_stats__}\n"
            unset __findRunTime_stats__
        fi
        if [ "${extractAll}" = yes ]; then
            printf "${__findRunTime_runtime__}\n" | sed 's/^/\t/'
        else
            if [ -n "${__findRunTime_runtime__}" ]; then
                noOutlier=$(printf "${__findRunTime_runtime__}" | head -1 | grep -v 'The runtime was [.0-9][.0-9]* seconds which is greater than the maximum allowed' || true)
                if [ -n "${noOutlier}" ]; then
                    echo -e "\t${noOutlier}"
                fi
            else
                echo -e "\tNo run time was available for that test"
            fi
        fi
        unset __findRunTime_runtime__
        if [ "${findOutliers}" = yes ]; then
            outlier=$(nodehistory -n ${findRunTime_nodes}/run -history 0 -edate ${logdate} | grep 'The runtime was [.0-9][.0-9]* seconds which is greater than the maximum allowed' | sed 's/%/%%/g')
            if [ -n "${outlier}" ]; then
                printf "${outlier}\n" | sed 's/^/\t/'
                line=$(printf "${outlier}" | sed 's/^/\t/' | sed 's/%/%%/g')
                outliers="${outliers}"${findRunTime_nodes%/*}"\n"${line}"\n"
            fi
        fi
    else
        for __node__ in ${findRunTime_nodes}; do
            findRunTime ${__node__}
        done
        unset __node__
    fi

    unset findRunTime_node findRunTime_nodes
}  ## End of function 'findRunTime'

## Initialize 'outliers' variable used in 'findRunTime'
[ "${findOutliers}" = yes ] && outliers=

export SEQ_EXP_HOME=${suite}
findRunTime /Tests

if [ "${findOutliers}" = yes ]; then
    if [ -n "${outliers}" ]; then
        echo "Some timing outliers were found"
        if [ -n "${emails}" ]; then
            echo "Sending a notification to '${emails}'"
            printf "We found some timing outliers in the timing in MIDAS test suite '${suite}':\n\n${outliers}\n" | mail -s "Timing outliers found in MIDAS test suite '${suite}'" ${emails}
        fi
    else
        echo "No timing outliers found"
    fi
fi
