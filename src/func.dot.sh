#! /bin/sh

##-----------------------------------------------
function is_compilation_done_frontend {
    set -e
    __jobid=$1

    while true; do
        __status=0

        jobchk -c ${FRONTEND} ${__jobid} || __status=$?
        if [ "${__status}" -ne 0 ]; then
            echo "Compilation on '${FRONTEND}' '${JOBNAME}' has finished."
            __listing=$(/bin/ls -t ${JOBNAME}.${FRONTEND}-*-$(hostname)-*.out | head -1)
            cat ${__listing}
            rm ${__listing}
            break
        fi
        sleep 5
    done
}

##-----------------------------------------------
function copy_depend {
    __target_machine=$1

    __revnum=$(../midas.version.sh)
    __arch=$( \
        (echo 'echo EC_ARCH=${EC_ARCH}' | ssh ${__target_machine} bash --login )\
        2>/dev/null | grep EC_ARCH | cut -d= -f2)
    __target_dir=${DIR_BLD_ROOT}/${__revnum}/${__arch}/
    __src_dir=${DIR_BLD_ROOT}/${__revnum}/${EC_ARCH}
    echo "Copy dependecies from ${__src_dir} to ${__target_machine}:${__target_dir}"
    [ -d ${__target_dir} ] || mkdir -p ${__target_dir}
    cp -p ${__src_dir}/dep.* ${__target_dir}
    cp -p ${__src_dir}/*.f90 ${__target_dir}
}
