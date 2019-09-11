#!/bin/bash

set -e

if [ $# -eq 0 ]; then
    MIDAS_ABS=
    codedir=${PWD}
elif [ $# -eq 1 ]; then
    MIDAS_ABS=${1}
    codedir=${PWD}
elif [ $# -eq 2 ]; then
    MIDAS_ABS=${1}
    codedir=${2}
else
    echo "$0 only takes 0, 1 or 2 arguments"
    exit 1
fi

SEQ_MAESTRO_SHORTCUT=${SEQ_MAESTRO_SHORTCUT:-". ssmuse-sh -d eccc/cmo/isst/maestro/1.5.3.3"}
which getdef 1>/dev/null 2>&1 || ${SEQ_MAESTRO_SHORTCUT}

suite=$(git rev-parse --show-toplevel)/maestro/suites/midas_system_tests
if [ -z "${COMPILING_MACHINE_PPP}" ]; then
    COMPILING_MACHINE_PPP=$(cd ${suite}; getdef --exp ${suite} resources/resources.def FRONTEND)
fi
if [ -z "${COMPILING_MACHINE_SUPER}" ]; then
    COMPILING_MACHINE_SUPER=$(cd ${suite}; getdef --exp ${suite} resources/resources.def BACKEND)
fi

rev=${CI_BUILD_REF:-$(git describe)}
jobname=${rev}_midasCompile

## get the number of programs
number_of_programs=$(/bin/ls -1 *.f90 | wc -l)

cat > compile_job <<EOF
#!/bin/bash

set -ex

export COMPILE_MIDAS_ADD_DEBUG_OPTIONS=${COMPILE_MIDAS_ADD_DEBUG_OPTIONS:-no}

cd ${codedir}
echo Launching compilation on '\${TRUE_HOST}' for platform '\${ORDENV_PLAT}'
yes '' | head -n ${number_of_programs} | ./compile_all.sh
EOF

#pbs_extra1='-Wblock=true'
#ord_soumet compile_job -jn ${jobname} -mach daley    -listing ${PWD} -w 60 -cpus 36

## Using as many cpus as there are programs to compile
jobid=$(ord_soumet compile_job -jn ${jobname} -mach eccc-ppp4 -listing ${PWD} -w 60 -cpus ${number_of_programs}  -m 8G)

## On evite d'attendre en queue en faisant un 'ssh' directement sur 'daley'
cat compile_job | ssh daley bash --login
rm compile_job

function is_compilation_done {
    set -e
    __is_compilation_done_host__=${1}

    if [ "${__is_compilation_done_host__}" = daley -o "${__is_compilation_done_host__}" = eccc-ppp4 ]; then
        # the jobname is cut with 15 characters by 'jobst'
        jobstname=$(echo ${jobname} | cut -c-15)
    else
        jobstname=${jobname}
    fi

    while true; do
        status=0
        ## jobst -c ${__is_compilation_done_host__} | grep ${USER} | grep "${jobstname}" || status=1
        jobchk -c ${__is_compilation_done_host__} ${jobid} || status=$?
        if [ "${status}" -ne 0 ]; then
            echo "The compilation on __is_compilation_done_host__ '${__is_compilation_done_host__}' with job '${jobname}' has finished."
            listing=$(/bin/ls -t ${jobname}.${__is_compilation_done_host__}-*-$(hostname)-*.out | head -1)
            cat ${listing}
            rm ${listing}
            break
        fi
        sleep 60
    done
    unset __is_compilation_done_host__
}

for host in eccc-ppp4; do
    is_compilation_done ${host}
done

status=0
echo "Checking if all programs have been compiled on '${TRUE_HOST}' for platform '${ORDENV_PLAT}'"
./check_if_all_programs_compiled.sh ${ORDENV_PLAT}          ${MIDAS_ABS} || status=1
echo "Checking if all programs have been compiled on '${host}' for platform 'sles-15-skylake-64-xc50'"
./check_if_all_programs_compiled.sh sles-15-skylake-64-xc50 ${MIDAS_ABS} || status=1

if [ "${status}" -eq 0 ]; then
    echo "All programs have been compiled correctly!"
else
    echo "Some programs could not compile!"
    exit 1
fi
