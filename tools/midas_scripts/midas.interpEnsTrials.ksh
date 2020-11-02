#!/bin/ksh

## Ce script fait l'interpolation pour un membre d'ensemble en
## particulier.  On y utilise amplement l'outil `pxs2pxt`.

set -ex

if [ $# -ne 6 ]; then
    echo "The script '$0' takes six arguments"
    exit 1
fi

jfile=$1
vertical_prototype=$2
dir_out=$3
file_out=$4
npex=$5
npey=$6

reflex=${REFLEX:-$(which reflex)}
editfst=${EDITFST:-$(which editfst)}
fstliste=${FSTLISTE:-$(which r.fstliste)}
rdate=${RDATE:-$(which r.date)}
pgsm=${PGSM:-$(which pgsm)}

cp_cmd=${cp_cmd:-$(which cp)}
mkdir_cmd=${mkdir_cmd:-$(which mkdir)}
rm_cmd=${rm_cmd:-$(which rm)}
basename_cmd=${basename_cmd:-$(which basename)}
date_cmd=${date_cmd:-$(which date)}
touch_cmd=${touch_cmd:-$(which touch)}
true_cmd=${true_cmd:-$(which true)}
mv_cmd=${mv_cmd:-$(which mv)}
echo_cmd=${echo_cmd:-$(which echo)}
printf_cmd=${printf_cmd:-$(which printf)}

if [ -f "${jfile}" ]; then
    [ -f "$(${basename_cmd} ${jfile})" ] && ${rm_cmd} -f $(${basename_cmd} ${jfile})
    ${cp_cmd} ${jfile} .
else
    ${scp_cmd:-r.srcp} ${jfile} .
fi

jfile=`${basename_cmd} ${jfile##*:}`

if [ "${INTERPENSTRIALS_DO_INTERPOLATION:-yes}" = yes ]; then

__damaged__=0
${reflex} -ixent ${jfile} -stats -errexit || __damaged__=1
if [ "${__damaged__}" -ne 0 ]; then
    number=$(${echo_cmd} ${jfile} | ${cut_cmd} -d_ -f3)
    ## indicate that this file is damaged
    ${touch_cmd} damaged
    exit 2
fi

OLD_TMPDIR=${TMPDIR}
export TMPDIR="."
list_datev_file=$(${fstliste} -izfst ${jfile} -col 11 -nomvar TT)
${echo_cmd} "LIST OF VALID DATES IN FILE: ${list_datev_file}"

trial_date=$(${echo_cmd} ${jfile} | cut -d_ -f1)
if [[ "${trial_date}" = [0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9] ]]; then
    begin_date=$(${rdate} ${trial_date} +3)
    end_date=$(${rdate} ${trial_date} +9)
    list_datev=""
    for trldate in ${list_datev_file}; do
	if [[ "${begin_date}" -le "${trldate}" && "${trldate}" -le "${end_date}" ]]; then
	    list_datev="${list_datev} ${trldate}"
	fi
    done
    ${echo_cmd} "LIST OF VALID DATES IN the 6 hours assimilation window: ${list_datev}"
fi

if [ -z "${add_toctoc}" ]; then
    if [ -z "$(which d.add_toctoc || true)" ]; then
	. s.ssmuse.dot ${INTERPENSTRIALS_VGRID_PACKAGE:-cmdn/201502/00/base}
	if [ -n "$(which d.add_toctoc || true)" ]; then
	    ${echo_cmd} "Using add_toctoc=$(which d.add_toctoc)"
	else
	    ${echo_cmd} "Cannot find the program 'd.add_toctoc' in the PATH"
	    exit 1
	fi
    fi
    add_toctoc=$(which d.add_toctoc)
fi
${echo_cmd} "Using add_toctoc=${add_toctoc}"

if [ -z "${pxs2pxt}" ]; then
    if [ -z "$(which d.pxs2pxt || true)" ]; then
	. ssmuse-sh -d ${INTERPENSTRIALS_CMDN_PACKAGE:-cmdn/pxs2pxt/3.5.0/base}
	if [ -n "$(which d.pxs2pxt || true)" ]; then
	    ${echo_cmd} "Using pxs2pxt=$(which d.pxs2pxt)"
	else
	    ${echo_cmd} "Cannot find the program 'd.pxs2pxt' in the PATH"
	    exit 1
	fi
    fi
    pxs2pxt=$(which d.pxs2pxt)
fi
${echo_cmd} "Using d.pxs2pxt=${pxs2pxt}"

[ -f vertical_prototype ] && ${rm_cmd} -f vertical_prototype
${cp_cmd} ${vertical_prototype} vertical_prototype

## Look into the prototype to see if it is the 'stg2nonstg' or 'stg2stg' mode
TOCTOC=$(${fstliste} -nomvar '!!' -col 1 -izfst vertical_prototype)
if [ "${TOCTOC}" = "!!" ]; then
    ig1_toctoc=$(${fstliste} -nomvar '!!' -col 15 -izfst vertical_prototype)
    if [ "${ig1_toctoc}" = 5001 ]; then
	__INTERPENSTRIALS_MODE__=stg2nonstg
    else
	__INTERPENSTRIALS_MODE__=stg2stg
    fi
    ${echo_cmd} "  desire(-1,['!!'])" > edit.dir
    ${editfst} -s vertical_prototype -d target_prototype -i edit.dir
    IS_P0_IN_VERTICAL_PROTOTYPE=$(${fstliste} -nomvar P0 -col 1 -izfst vertical_prototype)
    if [ -n "${IS_P0_IN_VERTICAL_PROTOTYPE}" ]; then
        ${echo_cmd}  "  desire(-1,['>>','^^','^>'])" > edit.dir
        ${editfst} -s vertical_prototype -d target_prototype_p0 -i edit.dir
        ${echo_cmd}  "  desire(-1,'P0')" > edit.dir
        ${echo_cmd}  "  zap(-1,'P0TA')" >> edit.dir
        ${editfst} -s vertical_prototype -d target_prototype_p0 -i edit.dir
    fi
else
    HYFIELD=$(${fstliste} -nomvar 'HY' -col 1 -izfst vertical_prototype)
    if [ "${HYFIELD}" = "HY" ]; then
	__INTERPENSTRIALS_MODE__=stg2nonstg
	${add_toctoc} -s vertical_prototype -d target_prototype
    else
	${echo_cmd} "Cannot find any '!!' or 'HY' in the vertical coordinate prototype file '${vertical_prototype}'"
	exit 1
    fi
fi

process_jdate () {
    set -ex

    jdate=${1}
    MAIN_DIR=${PWD}

    mkdir ${jdate}_work
    cd ${jdate}_work

    echo start jdate=${jdate} SECONDS=${SECONDS}

    date_encode=$(${rdate} ${jdate})
    ${rm_cmd} -f file_datev.fst file_datev_int.fst
    ${echo_cmd} " desire(-1,['>>','^^','^>','!!'])" > edit.dir
    ${echo_cmd} " desire(-1,-1,-1,${date_encode})" >> edit.dir
    ${echo_cmd} " end" >> edit.dir
    ${editfst} -s ${MAIN_DIR}/${jfile} -d file_datev.fst -n -i edit.dir

#====================================================================

# Get etiket for the output file
    CETIK=UNDEF
    CETIK=$(${fstliste} -nomvar P0 -col 9  -izfst ./file_datev.fst | ${tail_cmd} -1 | ${cut_cmd} -c1-12)
# Get date time stamp of the starting date of the trial
    DATEO=$(${fstliste} -nomvar P0 -col 10 -izfst ./file_datev.fst | ${tail_cmd} -1 | ${cut_cmd} -c1-12)
    STAMPO=$(${rdate} -Sn ${DATEO})

# Get date time stamp of the valid date of the trial
    DATEV=$(${fstliste} -nomvar P0 -col 11 -izfst ./file_datev.fst | ${tail_cmd} -1 | ${cut_cmd} -c1-12)
    STAMPV=$(${rdate} -Sn ${DATEV})

#====================================================================
# Prepare px_target

    ${cp_cmd} ${MAIN_DIR}/target_prototype px_target
    if [ -f ${MAIN_DIR}/target_prototype_p0 ]; then
        ## if there is a P0 in the vertical prototype then we must interpolate the P0 in 'file_datev.fst' onto that grid
        ${echo_cmd} "  SORTIE(STD,400,R)" > edit.dir
        ${echo_cmd} "  GRILLE(COMME, FENTREE, 'P0TA', '  ',-1, -1, -1, -1, '            ')" >> edit.dir
        ${echo_cmd} "  HEURE(TOUT)" >> edit.dir
        ${echo_cmd} "  CHAMP('P0',TOUT)" >> edit.dir
        ${echo_cmd} "  CHAMP('TG',TOUT)" >> edit.dir
        rm -f px_target_p0
        ${pgsm} -iment ./file_datev.fst ${MAIN_DIR}/target_prototype_p0 -ozsrt px_target_p0 -i edit.dir
        ## no need to have the '!!' copied from 'file_datev.fst' into 'px_target_p0' by P
        ${echo_cmd} "  desire(-1,['P0','>>','^^','^>'])" > edit.dir
        ${editfst} -s px_target_p0 -d px_target -i edit.dir
        ${echo_cmd} "  desire(-1,['TG','>>','^^','^>'])" > edit.dir
        ${editfst} -s px_target_p0 -d target.fst -i edit.dir
        PXS2PXT_CUB_UV=CUB_UV
    else
        ${echo_cmd} "  desire(-1,['P0','>>','^^','^>'])" > edit.dir
        ${editfst} -s ./file_datev.fst -d px_target  -i edit.dir
        ${echo_cmd} "  desire(-1,['TG','>>','^^','^>'])" > edit.dir
        ${editfst} -s ./file_datev.fst -d target.fst -i edit.dir
        PXS2PXT_CUB_UV="CUB_UU CUB_VV"
    fi

#====================================================================
# Do the interpolation

    PXT_LEVELS_M=MOMENTUM
    PXT_LEVELS_T=MOMENTUM
    if [ "${__INTERPENSTRIALS_MODE__}" = stg2stg ];then
        PXT_LEVELS_T=THERMO
    fi

    vars_momentum=${MIDAS_INTERPENSTRIALS_VARS_MOMENTUM:-${PXS2PXT_CUB_UV}}
    vars_thermo=${MIDAS_INTERPENSTRIALS_VARS_THERMO:-"CUB_TT LIN_HU"}

    ${rm_cmd} -f target_MM target_TH
    ${pxs2pxt} \
	-s   file_datev.fst \
	-pxs file_datev.fst \
	-pxt px_target \
	-d target_MM \
	-etiket ${CETIK} \
	-var ${vars_momentum} \
        -pxt_levels ${PXT_LEVELS_M} \
	-above 0.0

    ${pxs2pxt} \
	-s   file_datev.fst \
	-pxs file_datev.fst \
	-pxt px_target \
	-d target_TH \
	-etiket ${CETIK} \
	-var ${vars_thermo} \
        -pxt_levels ${PXT_LEVELS_T} \
	-above 0.0

    ${editfst} -s target_MM -d target.fst -i 0
    ${echo_cmd} " exclure(-1,['P0'])" > edit.dir
    ${echo_cmd} " end" >> edit.dir
    ${editfst} -s target_TH -d target.fst -i edit.dir

    ${rm_cmd} -f target_MM target_TH

    if [ "${__INTERPENSTRIALS_MODE__}" = stg2nonstg ]; then
        # Add HY from prototype to the output file
	${echo_cmd} " desire(-1,'HY')" > edit.dir
	${echo_cmd} " zap(-1,'HY','${CETIK}',${STAMPO})" >> edit.dir
    elif [ "${__INTERPENSTRIALS_MODE__}" = stg2stg ]; then
        # Add toc-toc from prototype to the output file
	${echo_cmd} " desire(-1,'!!')" > edit.dir
    else
	${echo_cmd} "The interpolation mode variable '\${__INTERPENSTRIALS_MODE__}' can only be 'stg2nonstg' or 'stg2stg' and not '${__INTERPENSTRIALS_MODE__}'"
	return 1
    fi
    ${echo_cmd} " end" >> edit.dir
    ${editfst} -s ${MAIN_DIR}/vertical_prototype -d target.fst -i edit.dir

    ${echo_cmd} "    ***FIN DE ${__INTERPENSTRIALS_MODE__}***"

    ${echo_cmd} "   ***FIN de la datev " ${jdate}
    rm -f file_datev.fst px_target edit.dir
    echo end jdate=${jdate} SECONDS=${SECONDS}
    cd ${MAIN_DIR}
}  ## End of function 'process_jdate'

if [ -n "${INTERPENSTRIALS_THREADS_PER_MEMBER}" -a -n "${MIDAS_INTERPENSTRIALS_THREADS_PER_MEMBER}" ]; then
    echo "The variable 'INTERPENSTRIALS_THREADS_PER_MEMBER' which is equal to '${INTERPENSTRIALS_THREADS_PER_MEMBER}' is no longer supported."
    echo "It has been replaced by variable 'MIDAS_INTERPENSTRIALS_THREADS_PER_MEMBER'."
    exit 1
fi

SECONDS=0
counter=0
for jdate in ${list_datev}; do
    process_jdate ${jdate} || ${touch_cmd} abort.${jfile}.${jdate} &
    let counter=counter+1
    if [ "${counter}" -eq "${MIDAS_INTERPENSTRIALS_THREADS_PER_MEMBER:-1}" ]; then
        wait
        counter=0
    fi
done

## wait for any other thread launched in background
wait
echo final time for jdate loop SECONDS=${SECONDS} with INTERPENSTRIALS_THREADS_PER_MEMBER=${INTERPENSTRIALS_THREADS_PER_MEMBER}

abort=0
for file in *_work/abort.${jfile}.*; do
    if [ -f "${file}" ]; then
        jdate=$(${basename_cmd} ${file} | ${cut_cmd} -d. -f3)
        echo "The thread for jdate=${jdate} aborted!"
        abort=1
    fi
done
if [ "${abort}" -ne 0 ]; then
    echo "At least one thread aborted"
    exit 1
fi

${rm_cmd} -f target.fst
for jdate in ${list_datev}; do
    ${editfst} -s ${jdate}_work/target.fst -d target.fst -i 0
done

## Si on ne separe pas en bandes de latitudes alors on doit faire ce traitement
## pour garantir la validation avec les executions precedentes
if [ "${INTERPENSTRIALS_INTERP2GAUSS:-no}" != no ]; then
    if [ "${INTERPENSTRIALS_INTERP2GAUSS}" = yes ]; then
	INTERPENSTRIALS_INTERP2GAUSS=300
    fi
    ${rm_cmd} -f file_datev_int_gg.fst
    ${echo_cmd} "Interpolated to a gauss grid defined by grille(GAUSS,$((INTERPENSTRIALS_INTERP2GAUSS*2)),${INTERPENSTRIALS_INTERP2GAUSS},GLOBAL)"
    ${pgsm} -iment target.fst -ozsrt file_datev_int_gg.fst -i <<EOF
 sortie(STD,1440,A)
 grille(GAUSS,$((INTERPENSTRIALS_INTERP2GAUSS*2)),${INTERPENSTRIALS_INTERP2GAUSS},GLOBAL)
 compress=oui
 heure(tout)
 champ(TT,TOUT)
 champ(HU,TOUT)
 champ(P0,TOUT)
 champ(TG,TOUT)
 champ(UV,TOUT)
 end
EOF
    ${editfst} -s file_datev_int_gg.fst -d target_gg.fst -n -i <<EOF
 exclure(-1,'!!')
 end
EOF
    ${mv_cmd} target.fst target.fst_no_pgsm
    ${mv_cmd} target_gg.fst target.fst

fi ## Fin du 'if [ "${npex}" = 1 -a "${npey}" = 1 ]'

else ## else relie au 'if [ "${INTERPENSTRIALS_DO_INTERPOLATION}" = yes ]'
    ${cp_cmd} ${jfile} target.fst
fi ## Fin du else relie 'if [ "${INTERPENSTRIALS_DO_INTERPOLATION}" = yes ]'

if [ "${npex}" = 1 -a "${npey}" = 1 ]; then
    ${cp_cmd} target.fst ${dir_out}/${file_out}
else
    # Now split horizontally to a set of latitude bands and possibly longitude bands
    ${rm_cmd} -f subdomain_*
    time ${INTERPENSTRIALS_WRITE_SUBDOMAINS_PROGRAM:-$(which write_subdomains.Abs)} target.fst ${npex} ${npey} subdomain_

    ## Les sorties auront le nom 'subdomain_0012_0005.fst'
    typeset -i latband

    typeset -i lonband=1
    while [[ "${lonband}" -le "${npex}" ]]; do
	latband=1
	lonband_str=`${printf_cmd} "%.4d" ${lonband}`
	while [[ "${latband}" -le "${npey}" ]]; do
	    latband_str=`${printf_cmd} "%.4d" ${latband}`
	    subdomain=${lonband_str}_${latband_str}

	    ${mkdir_cmd} -p ${dir_out}/subdomain_${subdomain}

	    if [ "${__INTERPENSTRIALS_MODE__}" = stg2nonstg ]; then
	    # Add toc-toc from prototype to the output file
		${echo_cmd} " exclure(-1,'!!')" > edit.dir
		${echo_cmd} "  end" >> edit.dir
		${editfst} -s subdomain_${subdomain}.fst -d subdomain_${subdomain}_final.fst -n -i edit.dir
		${rm_cmd} -f subdomain_${subdomain}.fst
	    elif [ "${__INTERPENSTRIALS_MODE__}" = stg2stg ]; then
		${mv_cmd} subdomain_${subdomain}.fst subdomain_${subdomain}_final.fst
	    else
		${echo_cmd} "The interpolation mode variable '\${__INTERPENSTRIALS_MODE__}' can only be 'stg2nonstg' or 'stg2stg' and not '${__INTERPENSTRIALS_MODE__}'"
		exit 1
	    fi

	    ${mv_cmd} subdomain_${subdomain}_final.fst ${dir_out}/subdomain_${subdomain}/${file_out}

	    let latband=latband+1
	done
	let lonband=lonband+1
    done
fi ## Fin du 'if [ "${npex}" != 1 -o "${npey}" != 1 ]'

[ -f file_datev.fst ]        && ${rm_cmd} -f file_datev.fst
[ -f file_datev_int_gg.fst ] && ${rm_cmd} -f file_datev_int_gg.fst
[ -f target.fst ]            && ${rm_cmd} -f target.fst
[ -f target_MM.fst ]         && ${rm_cmd} -f target_MM.fst
[ -f target_TH.fst ]         && ${rm_cmd} -f target_TH.fst
[ -f px_source_momentum ]    && ${rm_cmd} -f px_source_momentum
[ -f px_source_thermo ]      && ${rm_cmd} -f px_source_thermo
[ -f px_target ]             && ${rm_cmd} -f px_target

${echo_cmd} "   ***FIN INTERPOLATION***"
