#!/bin/ksh

set -ex

SECONDS=0

fasttmp=${1}
ramdiskpath=${2:-/tmp/${USER}}

## the observations are split into many files
splitobs=yes

npex=${SEQ_NPEX}
npey=${SEQ_NPEY}
let nprocs=npex*npey

typeset -Z4 obsfile_x obsfile_y

obsfile_y=$((npey-MP_CHILD/npex))
obsfile_x=$((MP_CHILD+1-npex*(MP_CHILD/npex)))

if [ "${fasttmp}" = yes ]; then
    FASTTMPDIR=${ramdiskpath}/oavar_${MP_CHILD}
    /bin/mkdir -p ${FASTTMPDIR}

    ## move observations in RAMDisk
    /bin/mkdir ${FASTTMPDIR}/obs
    if [ "${splitobs}" = no ]; then
        /bin/cp brp* ${FASTTMPDIR}/obs
    else
        for file in burpfiles_*/brp*_${obsfile_x}_${obsfile_y} \
                    sqlfiles_*/obs*_${obsfile_x}_${obsfile_y}; do
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
       if [ "$(echo burpfiles_*/brp*)" != "burpfiles_*/brp*" ]; then
          /bin/cp burpfiles_*/brp*_${obsfile_x}_${obsfile_y} ./obs
          export OAVAR_BURP_DIRECTORY=./obs
       fi
       if [ "$(echo sqlfiles_*/obs*)" != "sqlfiles_*/obs*" ]; then
          /bin/cp sqlfiles_*/obs*_${obsfile_x}_${obsfile_y} ./obs
       fi
    fi
fi

echo "The preparation of the working directory took ${SECONDS} seconds"

SECONDS=0
${RUN_PGM}
echo "The program itself took ${SECONDS} seconds"

SECONDS=0

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
                while [ ! -d burpfiles_${fam}.thinned ]; do
                    if [ "${MP_CHILD}" -eq 0 ]; then
                        /bin/mkdir burpfiles_${fam}.thinned
                        break
                    fi
                    sleep 1
                done
                /bin/cp ${file} burpfiles_${fam}.thinned/$(/usr/bin/basename ${file})
            fi
        done
    else
        if [ "$(echo burpfiles_*/brp*)" != "burpfiles_*/brp*" ]; then
            for file in ./obs/brp*_${obsfile_x}_${obsfile_y}; do
                fam=$(/usr/bin/basename ${file} | /usr/bin/rev | /usr/bin/cut -d_ -f3- | /usr/bin/rev | /usr/bin/cut -c4-)
                while [ ! -d burpfiles_${fam}.thinned ]; do
                    if [ "${MP_CHILD}" -eq 0 ]; then
                        /bin/mkdir burpfiles_${fam}.thinned
                        break
                    fi
                    sleep 1
                done
                /bin/mv ${file} burpfiles_${fam}.thinned/$(/usr/bin/basename ${file})
            done
        fi

        if [ "$(echo sqlfiles_*/obs*)" != "sqlfiles_*/obs*" ]; then
            for file in ./obs/obs*_${obsfile_x}_${obsfile_y}; do
                fam=$(/usr/bin/basename ${file} | /usr/bin/rev | /usr/bin/cut -d_ -f3- | /usr/bin/rev | /usr/bin/cut -c4-)
                while [ ! -d sqlfiles_${fam}.thinned ]; do
                    if [ "${MP_CHILD}" -eq 0 ]; then
                        /bin/mkdir sqlfiles_${fam}.thinned
                        break
                    fi
                    sleep 1
                done

#                # Clean deleted records from the SQL file
#                sqlite3 ${file} 'vacuum'

                /bin/mv ${file} sqlfiles_${fam}.thinned/$(/usr/bin/basename ${file})
            done
        fi
    fi
fi

if [ "${fasttmp}" = yes ]; then
    /bin/df -hP ${FASTTMPDIR}
    /bin/rm -rf ${FASTTMPDIR}
fi

echo "The finalization took ${SECONDS} seconds"
