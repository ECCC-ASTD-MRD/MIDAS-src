#!/bin/ksh

## Ce script gère l'interpolation des ensembles sur la grille horizontale
## et verticale de l'analyse.  Ce script appelle
## `midas.interpEnsTrials.ksh`.

set -ex

INTERPENSTRIALS_SSM_DOMAIN=${INTERPENSTRIALS_SSM_DOMAIN:-arma/envar/run/2.0.0}

INTERPENSTRIALS_MPIRUN=${INTERPENSTRIALS_MPIRUN:-midas.mpirun}
INTERPENSTRIALS_INTERP_SCRIPT=${INTERPENSTRIALS_INTERP_SCRIPT:-$(which midas.interpEnsTrials.ksh)}
INTERPENSTRIALS_RAMDISKPATH=${INTERPENSTRIALS_RAMDISKPATH-${TMPDIR}}

## Cette variable permet de decider si on fait l'interpolation verticale des membres.
export INTERPENSTRIALS_DO_INTERPOLATION=${INTERPENSTRIALS_DO_INTERPOLATION:-'yes'}
export INTERPENSTRIALS_INTERP2GAUSS=${INTERPENSTRIALS_INTERP2GAUSS:-'no'} ## si different de 'no' alors ce doit etre un chiffre qui indique la troncature

export INTERPENSTRIALS_VGRID_PACKAGE=${INTERPENSTRIALS_VGRID_PACKAGE:-eccc/cmd/cmdn/utils/20190812/09/default}
export INTERPENSTRIALS_CMDN_PACKAGE=${INTERPENSTRIALS_CMDN_PACKAGE:-eccc/cmd/cmdn/pxs2pxt/3.17.4/default}

#
# launch using: ord_soumet stg2nonstg_cub_enkf_driver_${VERSION_SCRIPT}.ksh -mach hadar -t 1800 -cpus 192x1 -mpi -args "script arguments"
#

if [ $# -lt 5 -o $# -gt 8 ]; then
    set +x
    echo "You must specify a source for the staggered EnKF files and a working directory and"
    echo "a prefix for the input files and the topology of the next MIDAS program execution (for separation into latitude and longitude bands)"
    echo "and possibly a prototype file defining the vertical coordinate (defaut: ${proto_file})"
    echo "and there is also the possibility to add the option 'menage=oui' as the last argument to remove the original staggered files"
    set -x
    exit 1
fi
stg_trial_dir=${1}
WORKDIR=${2}
prefix=${3}
npex=${4}
npey=${5}
proto_file=${6}
mpi_ssm=${7}
omp_stacksize=${8}
if [ $# -eq 9 ]; then
    if [ "${8}" = 'menage=oui' ]; then
	menage='oui'
    fi
fi

which midas.interpEnsTrials.ksh || . ssmuse-sh -d ${INTERPENSTRIALS_SSM_DOMAIN}
. ssmuse-sh -d ${INTERPENSTRIALS_VGRID_PACKAGE} ## get d.add_toctoc
. ssmuse-sh -d ${INTERPENSTRIALS_CMDN_PACKAGE}  ## get d.pxs2pxt

if [ -d "${WORKDIR}" ]; then
    rm -rf ${WORKDIR}
fi
mkdir ${WORKDIR}
cd ${WORKDIR}

let num_procs=${SEQ_NPEX:-1}*${SEQ_NPEY:-1}
echo 'num_procs= ' ${num_procs}

if [ -d "${stg_trial_dir##*:}" ]; then
    stg_trial_dir=${stg_trial_dir##*:}
    liste=$(/bin/ls ${stg_trial_dir}/${prefix}*)
    prefix_machine=
else
    liste=$(ssh ${stg_trial_dir%%:*} /bin/ls ${stg_trial_dir##*:}/${prefix}*)
    prefix_machine=${stg_trial_dir%%:*}:
fi

integer index
typeset -Z4 index
index=0

cp_cmd=$(which cp)
mkdir_cmd=$(which mkdir)
rm_cmd=$(which rm)
basename_cmd=$(which basename)
date_cmd=$(which date)
touch_cmd=$(which touch)
true_cmd=$(which true)
mv_cmd=$(which mv)
echo_cmd=$(which echo)
cut_cmd=$(which cut)
tail_cmd=$(which tail)
printf_cmd=$(which printf)

for jfile in ${liste}; do
    echo JOBINDEX=$index
    rm -fR ${WORKDIR}/workdir_$index
    mkdir -p ${WORKDIR}/workdir_$index
    cd ${WORKDIR}/workdir_$index
    ## This file is erased in 'task.sh' if everything has been executed correctly.
    touch abort

    cat << EOF > lance.sh
#!/bin/ksh
set -ex
${touch_cmd} abort.lance
status=0
time ${INTERPENSTRIALS_INTERP_SCRIPT} ${prefix_machine}${jfile} ${proto_file} ${WORKDIR} $(basename ${jfile}) ${npex} ${npey}
status=\$?
[ "\${status}" -eq 0 ] && ${rm_cmd} abort.lance
[ "\${status}" -eq 0 ] || exit \${status}
EOF
    cat lance.sh
    chmod +x ./lance.sh
    cp $proto_file ./proto_file_in.fst

    cd ../

    let index=index+1
done

if [ -n "${INTERPENSTRIALS_RAMDISKPATH}" ]; then
    cat > task.sh <<EOF
#!/bin/ksh
set -ex
typeset -Z4 job
let job=\${base}+\${MP_CHILD} || ${true_cmd} ## Pour eviter que ca retourne une erreur lorsque job=0+0
echo MP_CHILD = \${MP_CHILD} \$(${date_cmd}) ${WORKDIR}/workdir_\${job}
if [ -d ${WORKDIR}/workdir_\${job} ]; then
    RAMTMPDIR=${INTERPENSTRIALS_RAMDISKPATH}/interpEnsTrials_\${MP_CHILD}_\$\$
    ${mkdir_cmd} -pv \${RAMTMPDIR}
    cd \${RAMTMPDIR}
    abort=0
    ${WORKDIR}/workdir_\${job}/lance.sh > ${WORKDIR}/listing.job_\${job} 2>&1 || abort=1
    if [ "\${abort}" -ne 0 ]; then
        ${mv_cmd} -f * ${WORKDIR}/workdir_\${job} || true
    fi
    cd ${WORKDIR}
    ${rm_cmd} -rf \${RAMTMPDIR}
    if [ "\${abort}" -ne 0 ]; then
        exit 1
    fi
    rm ${WORKDIR}/workdir_\${job}/abort
fi
EOF
else
    cat > task.sh <<EOF
#!/bin/ksh
set -ex
typeset -Z4 job
let job=\${base}+\${MP_CHILD} || true ## Pour eviter que ca retourne une erreur lorsque job=0+0
${echo_cmd} MP_CHILD = \${MP_CHILD} \$(${date_cmd}) ${WORKDIR}/workdir_\${job}
if [ -d ${WORKDIR}/workdir_\${job} ]; then
    cd ${WORKDIR}/workdir_\${job}
    ${WORKDIR}/workdir_\${job}/lance.sh > ${WORKDIR}/listing.job_\${job} 2>&1
    rm ${WORKDIR}/workdir_\${job}/abort
fi
EOF
fi  ## Fin du 'else' relie au 'if [ -n "${INTERPENSTRIALS_RAMDISKPATH}" ]'

chmod +x task.sh

export REFLEX=$(which reflex)
export EDITFST=$(which editfst)
export FSTLISTE=$(which r.fstliste)
export RDATE=$(which r.date)
export PGSM=$(which pgsm)
export add_toctoc=$(which d.add_toctoc)
export pxs2pxt=$(which d.pxs2pxt)

export cp_cmd rm_cmd basename_cmd echo_cmd touch_cmd cut_cmd tail_cmd mv_cmd printf_cmd

LOCAL_BIN=${PWD}/local_bin
mkdir ${LOCAL_BIN}
cp $(which cclargs) $(which r.fstinfo) ${LOCAL_BIN}

typeset -Z4 job number
export base=0
while [ "${base}" -lt "${index}" ]; do
    echo CALLING POE for
    echo base=$base
    echo num_procs=$num_procs

    abort=0
    ${INTERPENSTRIALS_MPIRUN} ./task.sh "Lauching interpolation with ${INTERPENSTRIALS_INTERP_SCRIPT}" "${mpi_ssm}" ${omp_stacksize} ${TASK_BIN}/r.run_in_parallel -npex ${num_procs} -nompi pseudo_mpi -extrapath ${LOCAL_BIN} -tmpdir ${PWD}/interp_mpi_tmpdir || abort=1

    job=${base}
    while [ "${job}" -lt "$((base+num_procs))" ]; do
	if [ -d ${WORKDIR}/workdir_${job} ]; then
	    number=$((job+1))
	    if [ -f ${WORKDIR}/workdir_${job}/abort -o -f ${WORKDIR}/workdir_${job}/abort.lance ]; then
		echo "The job '${job}' processing the member ${number} did not end properly" >> aborted_jobs
	    fi
	    if [ -f ${WORKDIR}/workdir_${job}/damaged ]; then
		${SEQ_BIN}/nodelogger -n $SEQ_NODE -s info -m "Ensemble member ${number} is damaged"
	    fi
	fi
	let job=job+1
    done
    let base=${base}+${num_procs}
done

if [ "${menage}" = oui ]; then
    rm -fR ${WORKDIR}/workdir_*
    rm -f ${WORKDIR}/listing.job*
    if [ -d "${stg_trial_dir##*:}" ]; then
	stg_trial_dir=${stg_trial_dir##*:}
	rm -f ${stg_trial_dir}/${prefix}*
    else
	ssh ${stg_trial_dir%%:*} rm -f ${stg_trial_dir##*:}/${prefix}*
    fi
fi

if [ -s aborted_jobs ]; then
    echo "Problem detected while processing the ensemble members.  Here are the members that aborted"
    cat aborted_jobs
    exit 1
fi
