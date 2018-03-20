#!/bin/ksh

set -ex

fasttmp=${1}
ramdiskpath=${2:-/tmp/${USER}}

## the observations are splitted into many files
splitobs=yes

npex=${SEQ_NPEX}
npey=${SEQ_NPEY}
let nprocs=npex*npey

typeset -Z4 burpfile_x burpfile_y

burpfile_y=$((npey-MP_CHILD/npex))
burpfile_x=$((MP_CHILD+1-npex*(MP_CHILD/npex)))

if [ "${fasttmp}" = yes ]; then
    FASTTMPDIR=${ramdiskpath}/oavar_${MP_CHILD}
    /bin/mkdir -p ${FASTTMPDIR}

    ## move observations in RAMDisk
    /bin/mkdir ${FASTTMPDIR}/obs
    if [ "${splitobs}" = no ]; then
        /bin/cp brp* ${FASTTMPDIR}/obs
    else
        for file in burpfiles_*/brp*_${burpfile_x}_${burpfile_y}; do
            [ -f "${file}" ] && /bin/cp ${file} ${FASTTMPDIR}/obs
        done
    fi
    export OAVAR_BURP_DIRECTORY=${FASTTMPDIR}/obs

    ## Move Ensemble files in RAMDisk
    if [ -d ensemble ]; then
        /bin/mkdir ${FASTTMPDIR}/ensemble
        for file in ensemble/2*; do
        ## the file ensemble/2018020500_006_0128
        ##    -> 2018020500_006_0128
        ##      -> 8210_600_0050208102
        ##        -> 8210
        ##          -> 0128 (which is 'number')
            number=$(/usr/bin/basename ${file} | /usr/bin/rev | /usr/bin/cut -d_ -f1 | /usr/bin/rev)
            if [ "$(((number-1)%nprocs))" -eq "${MP_CHILD}" ]; then
                /bin/cp ensemble/*_${number} ${FASTTMPDIR}/ensemble
            fi
        done
        export OAVAR_ENSPATHNAME=${FASTTMPDIR}/ensemble
    fi

    ## Move Trial files in RAMDisk
    typeset -Z2 counter=1
    for file in trlm_* ; do
        if [ "$(((counter-1)%nprocs))" -eq "${MP_CHILD}" ]; then
            /bin/cp ${file} ${FASTTMPDIR}/
        fi
	let counter=counter+1
    done

    /bin/df -hP ${FASTTMPDIR}

    export OAVAR_RAM_DISK_DIR=${FASTTMPDIR}
else
    if [ "${splitobs}" != no ]; then
       [ "${MP_CHILD}" -eq 0 ] && /bin/mkdir ./obs
       while [ ! -d ./obs ]; do
           /bin/sleep 1
       done
       /bin/cp burpfiles_*/brp*_${burpfile_x}_${burpfile_y} ./obs
       export OAVAR_BURP_DIRECTORY=./obs
    fi
fi

${RUN_PGM}

if [ "${splitobs}" = no ]; then
    if [ "${fasttmp}" = yes ]; then
        if [ "${MP_CHILD}" -eq 0 ]; then
            /bin/cp ${FASTTMPDIR}/obs/brp* .
        fi
    fi
else
    if [ "${fasttmp}" = yes ]; then
        for file in ${FASTTMPDIR}/obs/brp*; do
            if [ -f "${file}" ]; then
                fam=$(/usr/bin/basename ${file} | /usr/bin/rev | /usr/bin/cut -d_ -f3- | /usr/bin/rev | /usr/bin/cut -c4-)
                while [ ! -d burpfiles_${fam}.updated ]; do
                    if [ "${MP_CHILD}" -eq 0 ]; then
                        /bin/mkdir burpfiles_${fam}.updated
                        break
                    fi
                    sleep 1
                done
                /bin/cp ${file} burpfiles_${fam}.updated/$(/usr/bin/basename ${file})
            fi
        done
    else
        for file in ./obs/brp*_${burpfile_x}_${burpfile_y}; do
            fam=$(/usr/bin/basename ${file} | /usr/bin/rev | /usr/bin/cut -d_ -f3- | /usr/bin/rev | /usr/bin/cut -c4-)
            while [ ! -d burpfiles_${fam}.updated ]; do
                if [ "${MP_CHILD}" -eq 0 ]; then
                    /bin/mkdir burpfiles_${fam}.updated
                    break
                fi
                sleep 1
            done
            /bin/mv ${file} burpfiles_${fam}.updated/$(/usr/bin/basename ${file})
        done
    fi
fi

if [ "${fasttmp}" = yes ]; then
    /bin/df -hP ${FASTTMPDIR}
    /bin/rm -rf ${FASTTMPDIR}
fi

