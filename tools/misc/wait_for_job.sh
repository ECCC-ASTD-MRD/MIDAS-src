#!/bin/bash

set -e

if [ $# -lt 3 ]; then
    echo "wait_for_job.sh: This scripts needs at least 3 arguments"
    exit 1
fi
if [ $# -gt 4 ]; then
    echo "wait_for_job.sh: This scripts takes at most 4 arguments"
    exit 1
fi

jobname=${1}
jobid=${2}
host=${3}
sleep_interval=${4:-5}

listing=$(/bin/ls -t ${jobname}.${host}-*-$(hostname)-*.out | head -1)

while true; do
    status=0
    jobchk -c ${host} ${jobid} || status=$?
    if [ "${status}" -ne 0 ]; then
        echo "The compilation on host '${host}' with job '${jobname}' has finished."
        if [ -f "${listing}" ]; then
            echo "======================================"
            echo "Begin of listing of '${jobname}'"
            cat ${listing}
            rm ${listing}
            echo "End of listing of '${jobname}'"
            echo "======================================"
        else
            echo "Cannot find listing: ${listing}"
        fi
        break
    fi
    sleep ${sleep_interval}
done
