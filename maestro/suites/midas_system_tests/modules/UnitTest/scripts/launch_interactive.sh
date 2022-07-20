#!/bin/bash

set -e

## 3 hours in minutes
wallclock_default=$((3*60))

if ! which qsub; then
    echo
    echo "This command only works when 'qsub' is available"
    echo "exit!"
    exit 1
fi

echo 'loading MAESTRO SSM ...'
which maestro 1>/dev/null 2>&1 || . ssmuse-sh -d eccc/cmo/isst/maestro/1.8.2

arguments=$*
eval `${CCLARGS:-cclargs} $0 \
  -exp   "not defined" "not defined" "[base maestro experiment" \
  -node  "not defined" "not defined" "[test to launch interactively]" \
  -date  "not defined" "not defined" "[date of the working directory prepared with maestro]" \
  -wallclock  "${wallclock_default}" "${wallclock_default}" "[wallclock time in minutes for the interactive job (default: 180 for 3 hours)]" \
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
        echo "extract_from_XML: needs exactly two arguments" >&2
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

unittestname=$(echo ${node} | sed 's|^/Tests/||' | sed 's|/UnitTest/run$||' | sed 's|/|.|g')
jobname=MIDAS.${unittestname}

# The following lines can be used for debugging:
#echo unittestname=${unittestname}
#echo jobname=${jobname}

echo
echo "Submitting an interactive job on ${host} with cpus=${cpus} memory=${memory} ${mpi:+with mpi} and ${soumet_args}"
echo

if [ "${TRUE_HOST}" != "${host}" ]; then
    echo "To launch on the interactive job, you must be on the same cluster ${host} as the targeted one" >&2
    exit 1
fi

## This function is checking if a file does already exist.
## If it does, then the function returns 1.
check_file() {
    set -e
    if [ $# -ne 2 ]; then
        echo "The function 'check_file' does accept only two arguments which are:" >&2
        echo "       function name" >&2
        echo "       file name" >&2
        echo "The function 'check_file' has been called with" >&2
        echo "       check_file $*" >&2
        return 1
    fi

    __check_file_function__=${1}
    __check_file_file__=${2}
    if [ -f "${__check_file_file__}" ]; then
        echo "${__check_file_function__}: The file ${__check_file_file__} exists" >&2
        echo "${__check_file_function__}: Erase it or move it if you want to keep it" >&2
        rm -i ${__check_file_file__}
    fi

    unset __check_file_file__ __check_file_function__
}

## This function is using 'ord_soumet' called with '-nosubmit' to
## extract the ressources computed by 'ord_soumet'.  We can then use
## those resources when calling 'jobsubi'.
find_resources () {
    set -e
    if [ $# -ne 0 ]; then
        echo "The function 'find_resources' does not accept any argument" >&2
        echo "The function 'find_resources' has been called with" >&2
        echo "       find_resources $*" >&2
        return 1
    fi

    __sleep_job__=${jobname}.sleep_forever.sh
    check_file find_resources ${__sleep_job__}

    __lajobtar__=lajob.tar
    check_file find_resources ${__lajobtar__}

    cat > ${__sleep_job__} <<EOF
#!/bin/bash

sleep $((360*60))

EOF
    ord_soumet_cmd="ord_soumet ${__sleep_job__} -jn ${jobname} -mach ${host} -listing ${PWD} -w 360 -cpus ${cpus} -m ${memory} ${soumet_args} -nosubmit 2>/dev/null"
    echo ord_soumet command: >&2
    echo ${ord_soumet_cmd} >&2

    ${ord_soumet_cmd}
    if [ $? -ne 0 ]; then
      echo Problem with call to ord_soumet
      exit 1
    fi

    rm ${__sleep_job__}

    ## The following command will output
    ##     #PBS -l select=10:ncpus=40:mem=100000M:res_tmpfs=5120:res_image=eccc/eccc_all_ppp_ubuntu-18.04-amd64_latest -l place=free -r n
    tar --wildcards -xOf ${__lajobtar__} "${jobname}.${host}*" | grep '^#PBS' | grep -v '^#PBS -[jo]' | grep -v '^#PBS -l walltime='

    rm ${__lajobtar__}

    unset __sleep_job__ __lajobtar__
} ## End of function 'find_resources'

pbsdirectives=${PWD}/${jobname}.pbsdirectives
check_file launch_interactive ${pbsdirectives}

find_resources > ${pbsdirectives}
echo "#PBS -l walltime=$((wallclock*60))" >> ${pbsdirectives}

rcfile=${PWD}/${jobname}.rcfile
check_file launch_interactive ${rcfile}

cat > ${rcfile} <<EOF
echo
echo Changing directory for ${working_directory}
cd ${working_directory}

echo
echo Loading execution environment in file load_env.sh
. ./load_env.sh

echo
echo You can now run your program interactively with:
echo
echo "   ./launch_program.sh pgm"
echo
echo for the latest executable used by 'run.tsk' or
echo
echo "   ./launch_program.sh /path/to/midas-bld/midas_abs/midas-\\\${pgm}_\\\${plat}_\\\${version}_\\\${commit}_\\\${M}.Abs"
echo
echo for any executable of the midas program used in this test.
EOF

#TODO mpiprocs should not be hard coded.
#TODO figure out how to launch rcfile and let interactive session open
qsub_cmd="qsub -X -I ${pbsdirectives}"
echo "Launching the interactive job with"
echo "   ${qsub_cmd}"
echo

echo "When the interactive job is started, do"
echo
echo "    source ${rcfile}"
echo
echo "to load the environment and start working"

${qsub_cmd}

rm ${rcfile} ${pbsdirectives}
