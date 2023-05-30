#!/bin/bash

set -ex

date_cmd=/bin/date

echo starting var.sh at $(${date_cmd} +%Y%m%d:%H:%M:%S.%N)

SECONDS=0

fasttmp=${1}
ramdiskpath=${PBSTMPFSDIR-/tmp/${USER}}

## assume the observations are splitted into many files

if [ "${fasttmp}" = yes ]; then
    FASTTMPDIR=${ramdiskpath}/midas_${MP_CHILD}
    [ -d "${FASTTMPDIR}" ] && /bin/rm -rf ${FASTTMPDIR}
    /bin/mkdir -p ${FASTTMPDIR}
    export MIDAS_RAMDISKDIR=${FASTTMPDIR}
    echo "Start of the job $(hostname) ${FASTTMPDIR} $(df -h ${FASTTMPDIR})"
fi

echo "The preparation of the working directory took ${SECONDS} seconds"

ulimit -c unlimited

echo starting ${RUN_PGM} at $(${date_cmd} +%Y%m%d:%H:%M:%S.%N)
SECONDS=0

#status=0
#/usr/bin/time --format "Total Memory: %M Mb" ${RUN_PGM} || status=1
${RUN_PGM}

echo "After program of the job $(hostname) ${FASTTMPDIR} $(df -h ${FASTTMPDIR})"

echo ending ${RUN_PGM} at $(${date_cmd} +%Y%m%d:%H:%M:%S.%N)
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

obsfile_y=$((1+MP_CHILD/SEQ_NPEX))
obsfile_x=$((MP_CHILD+1-SEQ_NPEX*(MP_CHILD/SEQ_NPEX)))

## The numbering for files must contain 4 characters padded with '0'.
obsfile_y=$(/usr/bin/printf "%0.4d" ${obsfile_y})
obsfile_x=$(/usr/bin/printf "%0.4d" ${obsfile_x})

if [ -d obsBeforeThinning ]; then
    for prefix in brp obs sql dia; do
        for file in ./obsBeforeThinning/${prefix}*_${obsfile_x}_${obsfile_y}; do
            if [ -f "${file}" ]; then
                if [[ "${file}" = *.num_headers ]]; then
                    /bin/rm ${file}
                    continue
                fi
                bfile=${file##*/}
                fam=$(extract_family ${file})
                if [ "${prefix}" = brp ]; then
                    updated_dir=burpfiles_${fam}.beforeThinning
                else
                    updated_dir=${prefix}files_${fam}.beforeThinning
                fi
                while [ ! -d "${updated_dir}" ]; do
                    if [ "${MP_CHILD}" -eq 0 ]; then
                        /bin/mkdir "${updated_dir}"
                        break
                    fi
                    /bin/sleep 1
                done
                /bin/mv ${file} ${updated_dir}/${bfile}
            fi
        done
    done
fi

if [ "${fasttmp}" = yes ]; then
    for obsDir in obs obsDB; do
        if [ ! -d "${FASTTMPDIR}/${obsDir}" ]; then
            continue
        fi
        for prefix in brp obs sql dia bcr; do
            for file in ${FASTTMPDIR}/${obsDir}/${prefix}*; do
                if [ -f "${file}" ]; then
                    if [[ "${file}" = *.num_headers ]]; then
                        continue
                    fi
                    bfile=${file##*/}
                    fam=$(extract_family ${file})
                    if [ "${obsDir}" == "obsDB" ]; then
                        bfile2=$(echo ${bfile} | sed "s/^${prefix}/odb/g")
                    else
                        bfile2=${bfile}
                    fi

                    if [ "${prefix}" = brp ]; then
                        updated_dir=burpfiles_${fam}.updated
                    else
                        updated_dir=${prefix}files_${fam}.updated
                    fi
                    while [ ! -d "${updated_dir}" ]; do
                        if [ "${MP_CHILD}" -eq 0 ]; then
                            /bin/mkdir ${updated_dir}
                            break
                        fi
                        /bin/sleep 1
                    done
                    /bin/cp ${file} ${updated_dir}/${bfile2}
                    /bin/rm -f ${obsDir}/${bfile}
                fi
            done
        done
    done
else
    for obsDir in obs obsDB; do
        if [ ! -d "./${obsDir}" ]; then
            continue
        fi
        for prefix in brp obs sql dia bcr; do
            for file in ./${obsDir}/${prefix}*_${obsfile_x}_${obsfile_y}; do
                if [ -f "${file}" ]; then
                    if [[ "${file}" = *.num_headers ]]; then
                        /bin/rm ${file}
                        continue
                    fi
                    bfile=${file##*/}
                    fam=$(extract_family ${file})
                    if [ "${obsDir}" == "obsDB" ]; then
                        bfile2=$(echo ${bfile} | sed "s/^${prefix}/odb/g")
                    else
                        bfile2=${bfile}
                    fi                    
                    
                    if [ "${prefix}" = brp ]; then
                        updated_dir=burpfiles_${fam}.updated
                    else
                        updated_dir=${prefix}files_${fam}.updated
                    fi
                    while [ ! -d "${updated_dir}" ]; do
                        if [ "${MP_CHILD}" -eq 0 ]; then
                            /bin/mkdir "${updated_dir}"
                            break
                        fi
                        /bin/sleep 1
                    done
                    /bin/mv ${file} ${updated_dir}/${bfile2}
                fi
            done
        done
    done
    
    if [ -d ./obsDB -a -z "$(ls -A ./obsDB)" ]; then
        /bin/rm -rf ./obsDB
    fi
fi

if [ "${fasttmp}" = yes ]; then
    /bin/df -hP ${FASTTMPDIR}
    /bin/rm -rf ${FASTTMPDIR}
fi

echo "The finalization took ${SECONDS} seconds"
echo "ending var.sh at $(/bin/date +%Y%m%d:%H:%M:%S.%N)"
