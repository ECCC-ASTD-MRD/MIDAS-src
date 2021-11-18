#!/bin/bash

set -e

arguments=$*
eval `${CCLARGS:-cclargs} $0 \
  -exp   "not defined" "not defined" "[base maestro experiment" \
  -node  "not defined" "not defined" "[test to launch interactively]" \
  -date  "not defined" "not defined" "[date of the working directory prepared with maestro]" \
  ++ $arguments`

if [ "${exp}" = 'not defined' ]; then
    echo "You must specify the maestro experiment associated with your code"
    exit 1
fi

if [ "${test}" = 'not defined' ]; then
    echo "You must specify the test you want to launch"
    exit 1
fi

if [ "${date}" = 'not defined' ]; then
    lastlog=$(/bin/ls -1t ${exp}/logs/*_nodelog | head -1)
    date=$(basename ${lastlog} _nodelog)
fi

date_number_of_characters=14
if [ "${#date}" -gt "${date_number_of_characters}" ]; then
    echo "The date '${date}' contains more than ${date_number_of_characters} characters and should not!"
    exit 1
fi

if [ "${#date}" -lt "${date_number_of_characters}" ]; then
    ## padding with zeros at the right
    ##  trick from https://stackoverflow.com/a/5349842
    padding=$(printf "0%.0s" $(seq $((${#date}+1)) ${date_number_of_characters}))
    date=${date}${padding}
fi


host=$(nodeinfo -e ${exp} -n ${node} | grep '^node\.machine=' | cut -d= -f2)
working_directory=${exp}/hub/${host}/work/${date}/${node}/work

extract_from_XML() {
    set -e

    if [ $# -ne 2 ]; then
        echo "extract_from_XML: needs exactly two arguments"
        exit 1
    fi

    property=${1}
    xmlfile=${2}

    xmllint --xpath "string(/NODE_RESOURCES/BATCH/@${property})" ${xmlfile}
}

resource=$(nodeinfo -e ${exp} -n ${node} | grep '^node\.resourcepath=' | cut -d= -f2 | sed 's|^${SEQ_EXP_HOME}|'"${exp}|")


cpus=$(extract_from_XML cpu ${resource})
isMPI=$(extract_from_XML mpi ${resource})
soumet_args=$(extract_from_XML soumet_args ${resource})
memory=$(extract_from_XML memory ${resource})

if [ "${isMPI}" = 1 ]; then
    mpi='-mpi'
else
    mpi=
fi

# The following lines can be used for debugging:
#echo cpus=${cpus}
#echo mpi=${mpi}
#echo memory=${memory}
#echo soumet_args=${soumet_args}

unittestname=$(echo ${node} | sed 's|^/Tests/||' | sed 's|/UnitTest/run$||' | sed 's|/|.|')
jobname=MIDAS.${unittestname}

# The following lines can be used for debugging:
#echo unittestname=${unittestname}
#echo jobname=${jobname}

sleep_job=${TMPDIR}/sleep_forever.sh
if [ -f "${sleep_job}" ]; then
    echo "The file ${sleep_job} exists"
    echo "Erase it or move it if you want to keep it"
    exit 1
fi

cat > ${sleep_job} <<EOF
#!/bin/bash

sleep $((360*60))

EOF

echo
echo "Submitting job ${jobname} on ${host} with cpus=${cpus} memory=${memory} ${mpi:+with mpi} and ${soumet_args}"
echo

jobid=$(ord_soumet ${sleep_job} -jn ${jobname} -mach ${host} -listing ${PWD} -w 360 -cpus ${cpus} -m ${memory} ${soumet_args})
rm ${sleep_job}

wait_for_job () {
    set -e

    while true; do
        jobstate=$(jobst -j ${jobid} -f | grep '^job_state=' | cut -d= -f2)
        if [ "${jobstate}" = R ]; then
            echo "Job ${jobid} just started"
            echo
            break
        elif [ -z "${jobstate}" ]; then
            echo "Did not find the job in the queue"
            exit 1
        fi
        echo
        echo "Waiting job ${jobid} on ${host} to start ..."
        echo
        sleep 2
    done
}

wait_for_job

cat <<EOF
Now you can log in the interactive job ${jobid} with
   sshj -j ${jobid}
   cd ${working_directory}
   . ./load_env.sh
   ./launch_programs.sh \${path_to_program}

EOF
