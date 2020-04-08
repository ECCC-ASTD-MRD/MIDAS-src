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
            rm ${__liisting}
            break
        fi
        sleep 5
    done
}

##-----------------------------------------------
function copy_depend {
    __target_machine=$1

    __arch=$( \
        (echo 'echo EC_ARCH=${EC_ARCH}' | ssh ${__target_machine} bash --login )\
        2>/dev/null | grep EC_ARCH | cut -d= -f2)
    for __prec in real4 real8
    do
        __target_dir=${DIR_BLD_ROOT}/${__arch}/${__prec}
        [ -d ${__target_dir} ] || mkdir -p ${DIR_BLD_ROOT}/${__arch}/${__prec}
        cp -p ${DIR_BLD_ROOT}/${EC_ARCH}/${__prec}/dep.${__prec}.* ${__target_dir}
        cp -p ${DIR_BLD_ROOT}/${EC_ARCH}/${__prec}/*.f90 ${__target_dir}
    done
}
