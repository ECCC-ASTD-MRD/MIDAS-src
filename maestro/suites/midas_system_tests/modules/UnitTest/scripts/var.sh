#!/bin/ksh

set -ex

SECONDS=0

fasttmp=${1}
ramdiskpath=${2:-/tmp/${USER}}

## assume the observations are splitted into many files

npex=${SEQ_NPEX}
npey=${SEQ_NPEY}
let nprocs=npex*npey

typeset -Z4 burpfile_x burpfile_y

burpfile_y=$((npey-MP_CHILD/npex))
burpfile_x=$((MP_CHILD+1-npex*(MP_CHILD/npex)))

if [ "${fasttmp}" = yes ]; then
    FASTTMPDIR=${ramdiskpath}/midas_${MP_CHILD}
    /bin/mkdir -p ${FASTTMPDIR}
    export MIDAS_RAMDISKDIR=${FASTTMPDIR}
fi

[ "${MP_CHILD}" -eq 0 ] && /bin/mkdir ./obs
while [ ! -d ./obs ]; do
    /bin/sleep 1
done
/bin/cp burpfiles_*/brp*_${burpfile_x}_${burpfile_y} ./obs

echo "The preparation of the working directory took ${SECONDS} seconds"

SECONDS=0

${RUN_PGM}
echo "The program itself took ${SECONDS} seconds"

SECONDS=0

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
            /bin/rm -f obs/$(/usr/bin/basename ${file})
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

if [ "${fasttmp}" = yes ]; then
    /bin/df -hP ${FASTTMPDIR}
    /bin/rm -rf ${FASTTMPDIR}
fi

echo "The finalization took ${SECONDS} seconds"
