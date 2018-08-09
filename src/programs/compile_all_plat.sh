#!/bin/bash

set -ex

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

rev=${CI_BUILD_REF:-$(git describe)}
jobname=${rev}_midasCompile

## get the number of programs
number_of_programs=$(/bin/ls -1 *.f90 | wc -l)

cat > compile_job <<EOF
#!/bin/bash

set -ex

cd ${codedir}
echo Launching compilation on '\${TRUE_HOST}' for platform '\${ORDENV_PLAT}'
yes '' | head -n ${number_of_programs} | ./compile_all.sh
EOF

#pbs_extra1='-Wblock=true'
#ord_soumet compile_job -jn ${jobname} -mach brooks    -listing ${PWD} -w 60 -cpus 36

## Using as many cpus as there are programs to compile
jobid=$(ord_soumet compile_job -jn ${jobname} -mach eccc-ppp1 -listing ${PWD} -w 60 -cpus ${number_of_programs}  -m 8G)

## On evite d'attendre en queue en faisant un 'ssh' directement sur 'brooks'
cat compile_job | ssh brooks bash --login
rm compile_job

function is_compilation_done {
    set -e
    __is_compilation_done_host__=${1}

    if [ "${__is_compilation_done_host__}" = brooks -o "${__is_compilation_done_host__}" = eccc-ppp1 ]; then
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

for host in eccc-ppp1; do
    is_compilation_done ${host}
done

echo "Checking if all programs have been compiled on '${TRUE_HOST}' for platform '${ORDENV_PLAT}'"
./check_if_all_programs_compiled.sh ${ORDENV_PLAT}            ${MIDAS_ABS}
echo "Checking if all programs have been compiled on '${host}' for platform 'sles-11-broadwell-64-xc40'"
./check_if_all_programs_compiled.sh sles-11-broadwell-64-xc40 ${MIDAS_ABS}
