#! /bin/sh

##  sourcing user configuration
source ./config.dot.sh

##  sourcing utilitary functions
source ./func.dot.sh


##=========================================================
##  Compilation on backend (in the background)
here=$(pwd)
cat > .compile_job << EOF
#! /bin/sh
set -ex
cd ${here}
source ${DOT_CONFIG} ; \
make all -j ${NCORES} \
    DIR_BLD_ROOT=${DIR_BLD_ROOT} VERBOSE=${VERBOSE} && \
make install  -j ${NCORES} \
    DIR_BLD_ROOT=${DIR_BLD_ROOT} VERBOSE=${VERBOSE}
EOF

echo "#####################################"
echo "... Launching compilation on BACKEND"
echo "    > listing_${BACKEND}" 
(cat .compile_job | ssh ${BACKEND} bash --login > listing_${BACKEND} 2>&1) &
PID_BACKEND=$!

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
                    -listing ${PWD} -w 60 -cpus ${NCORES}  -m 8G)

    ## waiting for frontend compilation to terminate
    is_compilation_done_frontend ${JOBID_FRONTEND} || \
        echo "Something went wrong on ${FRONTEND} with ${JOBID_FRONTEND}:" ;\
        cat .compile_job
fi

## waiting for backend compilation to terminate
wait ${PID_BACKEND}

if ${CLEAN}
then
    make cleanabs DIR_BLD_ROOT=${DIR_BLD_ROOT} VERBOSE=${VERBOSE}
    #rm -rf .compile_job
fi


echo "######################################"
echo "#"
echo "#  MIDAS COMPILATION COMPLETED"
echo "#"
echo "#  > ${DIR_BLD_ROOT}"
echo "#"
echo "######################################"

