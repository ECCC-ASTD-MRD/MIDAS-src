#! /bin/sh


##=========================================================
##  functions
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
##=========================================================


##  sourcing user configuration
source ./config.dot.sh


##=========================================================
##  Compilation on backend (in the background)
here=$(pwd)
cat > .compile_job << EOF
#! /bin/bash
set -ex
cd ${here}
source ${DOT_CONFIG} 
make all -j ${NCORES} -O  DIR_BLD_ROOT=${DIR_BLD_ROOT} VERBOSE=${VERBOSE} && \
make install  -j ${NCORES} -O DIR_BLD_ROOT=${DIR_BLD_ROOT} VERBOSE=${VERBOSE}
EOF

echo "#####################################"
echo "... Launching compilation on BACKEND"
echo "    > listing_${BACKEND}" 
(cat .compile_job | ssh ${BACKEND} bash --login > listing_${BACKEND} 2>&1) &  
PID_BACKEND=$!
echo "BACKEND compilation process: ${PID_BACKEND}"

##=========================================================
##  Compilation on frontend
echo "######################################"
echo "... Launching compilation on FRONTEND"
echo "    > listing_ppp"
if ${DIRECT_FRONTEND_COMPILE}
then
    ## compile directly on head node
    make all -j ${NCORES} \
        DIR_BLD_ROOT=${DIR_BLD_ROOT} VERBOSE=${VERBOSE} > listing_ppp 2>&1 
    make install -j ${NCORES} \
        DIR_BLD_ROOT=${DIR_BLD_ROOT} VERBOSE=${VERBOSE} >> listing_ppp 2>&1
else
    ## use ord_soumet
    JOBID_FRONTEND=$(ord_soumet .compile_job -jn ${JOBNAME} -mach ${FRONTEND} \
                    -listing ${PWD} -w 60 -cpus ${NCORES}  -m 8G -shell /bin/bash)

    ## waiting for frontend compilation to terminate
    is_compilation_done_frontend ${JOBID_FRONTEND} || \
        echo "Something went wrong on ${FRONTEND} with ${JOBID_FRONTEND}:" ;\
        cat .compile_job
fi

## waiting for backend compilation to terminate
echo "waiting for BACKEND (${PID_BACKEND}) compilation to finish"
wait ${PID_BACKEND}

if ${CLEAN}
then
    make cleanabs DIR_BLD_ROOT=${DIR_BLD_ROOT} VERBOSE=${VERBOSE}
    rm -rf .compile_job
fi


echo "######################################"
echo "#"
echo "#  MIDAS COMPILATION COMPLETED"
echo "#"
echo "#  > ${DIR_BLD_ROOT}"
echo "#"
echo "######################################"

