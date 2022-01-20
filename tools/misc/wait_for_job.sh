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

logLine() {
    set -e
    printf "\t${*}\n"
} ## End of function 'logLine'

logLine
logLine "Waiting for job '${jobid}' on host '${host}' with job '${jobname}' to finish"
logLine

findListing() {
    set -e
    __findListing_jobname__=${1}
    __findListing_host__=${2}

    /bin/ls -t ${__findListing_jobname__}.${__findListing_host__}-*-$(hostname)-*.out 2>/dev/null | head -1

    unset __findListing_jobname__ __findListing_host__
} ## End of function 'findListing'

listing=$(findListing ${jobname} ${host} || true)

linesListed=0
while true; do
    status=0
    if [ "${ORDENV_PLAT}" = rhel-8-icelake-64 ]; then
        echo qstat ${jobid} | ssh ${host} bash --login 1>/dev/null 2>&1 || status=$?
    else
        jobchk -c ${host} ${jobid} || status=$?
    fi

    if [ -f "${listing}" ]; then
        if [ "${linesListed}" -eq 0 ]; then
            logLine ======================================
            logLine "Begin of listing of '${jobname}'"
            logLine
        fi
        let oldLines=linesListed+1
        linesListed=$(cat ${listing} | wc -l)
        ## Print from line '${oldLines}' to line '${linesListed}' in file '${listing}'
        sed -n -e "${oldLines},${linesListed}p" ${listing}
    else
        listing=$(findListing ${jobname} ${host} || true)
    fi

    if [ "${status}" -ne 0 ]; then
        if [ -f "${listing}" ]; then
            let oldLines=linesListed+1
            tail -n +${oldLines} ${listing}
            rm ${listing}
            logLine
            logLine "End of listing of '${jobname}'"
            logLine ======================================
            logLine
            logLine "The job '${jobid}' on host '${host}' with job '${jobname}' has finished."
        else
            echo "Cannot find listing: ${listing}"
        fi
        break
    fi
    sleep ${sleep_interval}
done
