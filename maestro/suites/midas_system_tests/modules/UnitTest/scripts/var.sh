#!/bin/bash

set -ex

SECONDS=0

fasttmp=${1}
ramdiskpath=${2:-/tmp/${USER}}

## assume the observations are splitted into many files

if [ "${fasttmp}" = yes ]; then
    FASTTMPDIR=${ramdiskpath}/midas_${MP_CHILD}
    /bin/mkdir -p ${FASTTMPDIR}
    export MIDAS_RAMDISKDIR=${FASTTMPDIR}
    /bin/df -hP ${FASTTMPDIR}
fi

echo "The preparation of the working directory took ${SECONDS} seconds"

SECONDS=0
${RUN_PGM}
echo "The program itself took ${SECONDS} seconds"

SECONDS=0

extract_family () {
    set -ex
    file=$1
    ## on prend le 'basename'
    file1=${file##*/}
    ## On enleve la partie '_0014_0008'
    file2=${file1%_*_*}
    ## On enleve la partie 'brp' ce qui nous donne la famille
    echo ${file2:3}
}

if [ "${fasttmp}" = yes ]; then
    for file in ${FASTTMPDIR}/obs/brp*; do
        if [ -f "${file}" ]; then
            if [[ "${file}" = *.num_headers ]]; then
                continue
            fi
            bfile=${file##*/}
            fam=$(extract_family ${file})
            while [ ! -d burpfiles_${fam}.updated ]; do
                if [ "${MP_CHILD}" -eq 0 ]; then
                    /bin/mkdir burpfiles_${fam}.updated
                    break
                fi
                /bin/sleep 1
            done
            /bin/cp ${file} burpfiles_${fam}.updated/${bfile}
            /bin/rm -f obs/${bfile}
        fi
    done

    for file in ${FASTTMPDIR}/obs/obs*; do
        if [ -f "${file}" ]; then
            if [[ "${file}" = *.num_headers ]]; then
                continue
            fi
            bfile=${file##*/}
            fam=$(extract_family ${file})
            while [ ! -d sqlfiles_${fam}.updated ]; do
                if [ "${MP_CHILD}" -eq 0 ]; then
                    /bin/mkdir sqlfiles_${fam}.updated
                    break
                fi
                sleep 1
            done
            /bin/cp ${file} sqlfiles_${fam}.updated/${bfile}
            /bin/rm -f obs/${bfile}
        fi
    done

else
    obsfile_y=$((SEQ_NPEY-MP_CHILD/SEQ_NPEX))
    obsfile_x=$((MP_CHILD+1-SEQ_NPEX*(MP_CHILD/SEQ_NPEX)))

    ## The numbering for files must contain 4 characters padded with '0'.
    obsfile_y=$(/usr/bin/printf "%0.4d" ${obsfile_y})
    obsfile_x=$(/usr/bin/printf "%0.4d" ${obsfile_x})

    for file in ./obs/brp*_${obsfile_x}_${obsfile_y}; do
        if [ -f "${file}" ]; then
            if [[ "${file}" = *.num_headers ]]; then
                /bin/rm ${file}
                continue
            fi
            bfile=${file##*/}
            fam=$(extract_family ${file})
            while [ ! -d burpfiles_${fam}.updated ]; do
                if [ "${MP_CHILD}" -eq 0 ]; then
                    /bin/mkdir burpfiles_${fam}.updated
                    break
                fi
                /bin/sleep 1
            done
            /bin/mv ${file} burpfiles_${fam}.updated/${bfile}
        fi
    done

    for file in ./obs/obs*_${obsfile_x}_${obsfile_y}; do
        if [ -f "${file}" ]; then
            if [[ "${file}" = *.num_headers ]]; then
                /bin/rm ${file}
                continue
            fi
            bfile=${file##*/}
            fam=$(extract_family ${file})
            while [ ! -d sqlfiles_${fam}.updated ]; do
                if [ "${MP_CHILD}" -eq 0 ]; then
                    /bin/mkdir sqlfiles_${fam}.updated
                    break
                fi
                /bin/sleep 1
            done
            /bin/mv ${file} sqlfiles_${fam}.updated/${bfile}
        fi
    done
fi

if [ "${fasttmp}" = yes ]; then
    /bin/df -hP ${FASTTMPDIR}
    /bin/rm -rf ${FASTTMPDIR}
fi

echo "The finalization took ${SECONDS} seconds"
echo "Ending var.sh at $(/bin/date +%Y%m%d:%H:%M:%S.%N)"
