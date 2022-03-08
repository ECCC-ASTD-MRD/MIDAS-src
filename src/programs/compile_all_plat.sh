#!/bin/bash

set -e
__toplevel=$(git rev-parse --show-toplevel)

## env. variable new naming convention + retrocompatibility
source ${__toplevel}/src/programs/commons/retroComp_warning.sh

# set the resources.def file, which depends on the TRUE_HOST name
../../set_resources_def.sh

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

suite=${__toplevel}/maestro/suites/midas_system_tests
compiledir_main=${MIDAS_COMPILE_DIR_MAIN:-"../../compiledir"}

if [ -z "${COMPILING_MACHINE_PPP}" -o -z "${COMPILING_MACHINE_SUPER}" ]; then
    ${__toplevel}/set_resources_def.sh
    . ${suite}/set_machine_list.dot

    if [ -z "${COMPILING_MACHINE_PPP}" ]; then
        COMPILING_MACHINE_PPP=${MACHINE_PPP}
    fi
    if [ -z "${COMPILING_MACHINE_SUPER}" ]; then
        COMPILING_MACHINE_SUPER=${MACHINE_SUPER}
    fi
fi

if [ "${COMPILING_MACHINE_SUPER}" = daley -o "${COMPILING_MACHINE_SUPER}" = banting ]; then
    PLAT_SUPER=sles-15-skylake-64-xc50
fi

if [ "${ORDENV_PLAT}" = rhel-8-icelake-64 ]; then
    ## On 'rhel-8-icelake-64', we can compile only on PPP
    unset COMPILING_MACHINE_SUPER
fi

rev=${CI_BUILD_REF:-$(${__toplevel}/midas.version.sh)}
jobname=${rev}_midasCompile

## get the number of programs
number_of_programs=$(/bin/ls -1 *.f90 | wc -l)

cat > compile_job <<EOF
#!/bin/bash

set -ex

export MIDAS_COMPILE_ADD_DEBUG_OPTIONS=${MIDAS_COMPILE_ADD_DEBUG_OPTIONS:-no}

cd ${codedir}
echo Launching compilation on "\${TRUE_HOST}" for platform "\${ORDENV_PLAT}"
yes '' | head -n ${number_of_programs} | ./compile_all.sh
EOF

cat > compile_splitobs <<EOF
#!/bin/bash

set -ex

cd ${__toplevel}/tools/splitobs
make splitobs_\${ORDENV_PLAT}

splitobs_pgm=midas.splitobs_\${ORDENV_PLAT}-\$(../../midas.version.sh).Abs
make install PGM=${MIDAS_ABS:-${compiledir_main}/midas_abs}/\${splitobs_pgm}

EOF

#pbs_extra1='-Wblock=true'
#ord_soumet compile_job -jn ${jobname} -mach ${COMPILING_MACHINE_SUPER} -listing ${PWD} -w 60 -cpus 36

set -x
jobname=${rev}_midasCompile
## Using as many cpus as there are programs to compile
jobid=$(ord_soumet compile_job -jn ${jobname} -mach ${COMPILING_MACHINE_PPP} -listing ${PWD} -w 60 -cpus ${number_of_programs} -m 24G -tmpfs 2G)
## Compile 'midas.splitobs'
jobname_splitobs=${rev}_splitobsCompile
jobid_splitobs=$(ord_soumet compile_splitobs -jn ${jobname_splitobs} -mach ${COMPILING_MACHINE_PPP} -listing ${PWD} -w 60 -cpus 1 -m 4G -tmpfs 1G)

status=0
if [ -n "${COMPILING_MACHINE_SUPER}" ]; then
    ## On evite d'attendre en queue en faisant un 'ssh' directement sur '${COMPILING_MACHINE_SUPER}'
    cat compile_job | ssh ${COMPILING_MACHINE_SUPER} bash --login || status=1
fi
rm compile_job

## if previous compilation aborted, then kill the compilation job
if [ "${status}" -ne 0 ]; then
    if [ "${ORDENV_PLAT}" = rhel-8-icelake-64 ]; then
        echo qdel ${jobid} | ssh ${COMPILING_MACHINE_PPP} bash --login
    else
        jobdel -c ${COMPILING_MACHINE_PPP} ${jobid}
    fi
    echo "Compilation aborted!"
    exit 1
fi

## then wait for all other programs to compile
${__toplevel}/tools/misc/wait_for_job.sh ${jobname} ${jobid} ${COMPILING_MACHINE_PPP}

status=0
echo "Checking if all programs have been compiled on '${TRUE_HOST}' for platform '${ORDENV_PLAT}' in ${MIDAS_ABS}"
./check_if_all_programs_compiled.sh ${ORDENV_PLAT} ${MIDAS_ABS} || status=1
if [ -n "${COMPILING_MACHINE_SUPER}" ]; then
    echo "Checking if all programs have been compiled on '${host}' for platform '${PLAT_SUPER}' in ${MIDAS_ABS}"
    ./check_if_all_programs_compiled.sh ${PLAT_SUPER} ${MIDAS_ABS} || status=1
fi

if [ "${status}" -eq 0 ]; then
    echo "All programs have been compiled correctly!"
else
    echo "Some programs could not compile!"
    exit 1
fi
