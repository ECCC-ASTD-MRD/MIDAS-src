#!/bin/bash

set -e

## 6 hours in seconds
wallclock_default=$((6*60*60))

arguments=$*
eval `${CCLARGS:-cclargs} $0 \
  -exp   "not defined" "not defined" "[base maestro experiment" \
  -node  "not defined" "not defined" "[test to launch interactively]" \
  -date  "not defined" "not defined" "[date of the working directory prepared with maestro]" \
  -wallclock  "${wallclock_default}" "${wallclock_default}" "[wallclock time for the interactive job (default: 6 hours)]" \
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

which nodeinfo 1>/dev/null 2>&1 || . ssmuse-sh -d eccc/cmo/isst/maestro/1.7.1

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

unittestname=$(echo ${node} | sed 's|^/Tests/||' | sed 's|/UnitTest/run$||' | sed 's|/|.|')
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
        return 1
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

    __sleep_job__=sleep_forever.sh
    check_file find_resources ${__sleep_job__}

    __lajobtar__=lajob.tar
    check_file find_resources ${__lajobtar__}

    cat > ${__sleep_job__} <<EOF
#!/bin/bash

sleep $((360*60))

EOF
    ord_soumet ${__sleep_job__} -jn ${jobname} -mach ${host} -listing ${PWD} -w 360 -cpus ${cpus} -m ${memory} ${soumet_args} -nosubmit 2>/dev/null
    rm ${__sleep_job__}

    ## The following command will output
    ##     #PBS -l select=10:ncpus=40:mem=100000M:res_tmpfs=5120:res_image=eccc/eccc_all_ppp_ubuntu-18.04-amd64_latest -l place=free -r n
    __PBSresources__=$(tar --wildcards -xOf ${__lajobtar__} "${jobname}.${host}*" | grep '^#PBS -l select=')

    rm ${__lajobtar__}

    # The following line can be used for debugging:
    # echo __PBSresources__=${__PBSresources__} >&2

    ## extract the memory setting
    echo ${__PBSresources__} | grep -oE "mem=[0-9]+[KMG]"
    ## extract the select setting (nslots)
    echo ${__PBSresources__} | grep -oE "select=[0-9]+"
    ## extract the npus setting
    echo ${__PBSresources__} | grep -oE "ncpus=[0-9]+"

    unset __sleep_job__ __lajobtar__
} ## End of function 'find_resources'

extract_resource() {
    set -e
    if [ $# -ne 2 ]; then
        echo "The function 'extract_resource' accepts only two arguments." >&2
        echo "    It has been called with:" >&2
        echo "            $0 $*" >&2
        return 1
    fi
    __element__=${1}
    __resources__=${2}

    printf "${__resources__}\n" | grep "^${__element__}=" | cut -d= -f2

    unset __element__ __resources__
}

resources=$(find_resources)
# The following line can be used for debugging:
# printf "resources=${resources}\n"

ncores=$(extract_resource ncpus "${resources}")
nslots=$(extract_resource select "${resources}")
memory_per_node=$(extract_resource mem "${resources}")

## Parsing 'soumet_args' to find the '-tmpfs ${TMPFSSIZE}'
set -- ${soumet_args}
other_resources=
while [ $# -ne 0 ]; do
    if [ "${1}" = -tmpfs ]; then
        if [ $# -lt 2 ]; then
            echo "The argument '-tmpfs' under 'soumet_args' needs an argument" >&2
            exit 1
        else
            other_resources="${other_resources} -r tmpfs=${2}"
            shift 2
        fi
    else
        shift
    fi
done

rcfile=${PWD}/rcfile
check_file launch_interactive ${rcfile}

cat > ${rcfile} <<EOF
. /etc/profile
. $HOME/.profile

if [ -z "\${PBS_JOBNAME}" ]; then
    profile=\$(/bin/ls /var/tmp/pbs.*/profile.sh)
    if [ -f "\${profile}" ]; then
        echo "Sourcing PBS setup in file \${profile}"
        source \${profile}
    else
        echo "Cannot find the PBS setup \${profile}"
    fi
fi

echo
echo Changing directory for ${working_directory}
cd ${working_directory}

echo
echo Loading execution environment in file load_env.sh
. ./load_env.sh

echo
echo You can now run interactively your program with
echo "   ./launch_program.sh \\\${path_to_program}"
echo
EOF

submit_host=$(echo ${host} | cut -d- -f2) ## remove 'eccc-' from 'eccc-ppp3'

which jobsubi 1> /dev/null 2>&1 || . ssmuse-sh -x hpco/exp/jobsubi/jobsubi-0.3

jobsubi_cmd="jobsubi -r memory=${memory_per_node} -r nslots=${nslots} -r ncores=${ncores} -r wallclock=$((6*60*60)) ${other_resources} ${submit_host} -- bash --rcfile ${rcfile} -i"
echo "Launching the interactive job with"
echo "   ${jobsubi_cmd}"
echo

${jobsubi_cmd}

rm ${rcfile}
