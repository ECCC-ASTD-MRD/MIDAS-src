#!/bin/ksh

set -ex

if [ $# -ne 4 ]; then
    echo "The script '$0' takes 4 arguments"
    exit 1
fi

jfile=$1
grid_prototype=$2
dir_out=$3
file_out=$4

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

[ -f grid_prototype ] && ${rm_cmd} -f grid_prototype
${cp_cmd} ${grid_prototype} grid_prototype

## Look into the prototype to see if it is the 'stg2nonstg' or 'stg2stg' mode
TOCTOC=$(${fstliste} -nomvar '!!' -col 1 -izfst grid_prototype)
if [ "${TOCTOC}" = "!!" ]; then
    ig1_toctoc=$(${fstliste} -nomvar '!!' -col 15 -izfst grid_prototype)
    if [ "${ig1_toctoc}" = 5001 ]; then
	__INTERPENSTRIALS_MODE__=stg2nonstg
    else
	__INTERPENSTRIALS_MODE__=stg2stg
    fi
    ${echo_cmd} "  desire(-1,['!!'])" > edit.dir
    ${editfst} -s grid_prototype -d target_prototype -i edit.dir
else
    HYFIELD=$(${fstliste} -nomvar 'HY' -col 1 -izfst grid_prototype)
    if [ "${HYFIELD}" = "HY" ]; then
	__INTERPENSTRIALS_MODE__=stg2nonstg
	${add_toctoc} -s grid_prototype -d target_prototype
    else
	${echo_cmd} "Cannot find any '!!' or 'HY' in the vertical coordinate prototype file '${grid_prototype}'"
	exit 1
    fi
fi

## Look into the prototype to see if horionzal interpolation is demanded
if [ -n "$(${fstliste} -izfst grid_prototype -nomvar 'P0' -col 1)" ]; then
    ${echo_cmd} "Horizontal Interpolation WILL BE PERFORMED"
    DO_HORIZ_INTERP=yes
    # Is it a LAM grid or a GLOBAL GRID
    if [ -n "$(${fstliste} -etiket "COREGRID" -col 9 -izfst grid_prototype)" ] ; then
	GRID_DOMAIN=lam
    else
	GRID_DOMAIN=global
    fi
else
    ${echo_cmd} "NO horizontal interpolation will be performed"
    DO_HORIZ_INTERP=no
fi

if [ "$DO_HORIZ_INTERP" = yes ] ; then
    ${echo_cmd} "Horizontal Interpolation WILL BE PERFORMED using tic-tic & tac-tac"
    if [ "$GRID_DOMAIN" = lam ] ; then
	${echo_cmd} "  desire(-1,['>>','^^'],'COREGRID')" > edit.dir
    else
	${echo_cmd} "  desire(-1,['>>','^^'],'ANALYSIS')" > edit.dir
    fi   
    ${editfst} -s grid_prototype -d target_prototype -i edit.dir
elif [ -n "$(${fstliste} -izfst grid_prototype -nomvar '^>' -col 1)" ]; then
    ${echo_cmd} ""
    ${echo_cmd} " !!!! Horizontal interpolation toward a grid using tic-tac (^>) IS NOT POSSIBLE !!!"
    exit 1
fi

rm -f target.fst
for jdate in ${list_datev}; do
    date_encode=$(${rdate} ${jdate})
    ${rm_cmd} -f file_datev.fst file_datev_int.fst
    ${echo_cmd} " desire(-1,'^^')" > edit.dir
    ${echo_cmd} " desire(-1,'>>')" >> edit.dir
    ${echo_cmd} " desire(-1,'!!')" >> edit.dir
    ${echo_cmd} " desire(-1,-1,-1,${date_encode})" >> edit.dir
    ${echo_cmd} " end" >> edit.dir
    ${editfst} -s ${jfile} -d file_datev.fst -n -i edit.dir

#====================================================================

# Get etiket for the output file
    CETIK=UNDEF
    CETIK=$(${fstliste} -nomvar P0 -col 9 -izfst ./file_datev.fst | ${tail_cmd} -1 | ${cut_cmd} -c1-12)
# Get date time stamp of the starting date of the trial
    DATEO=$(${fstliste} -nomvar P0 -col 10 -izfst ./file_datev.fst | ${tail_cmd} -1 | ${cut_cmd} -c1-12)
    STAMPO=$(${rdate} -Sn ${DATEO})

# Get date time stamp of the valid date of the trial
    DATEV=$(${fstliste} -nomvar P0 -col 11 -izfst ./file_datev.fst | ${tail_cmd} -1 | ${cut_cmd} -c1-12)
    STAMPV=$(${rdate} -Sn ${DATEV})

#====================================================================
# Prepare px_target

    cp target_prototype px_target

    if [ "$DO_HORIZ_INTERP" = yes ] ; then

	# Extract P0 without the !! because PGSM automatically copy the !!
	rm -f P0_source.fst
	${echo_cmd} "  desire(-1,['P0','>>','^^'])" > edit.dir
	${editfst} -s ./file_datev.fst -d P0_source.fst -i edit.dir

	# Inperpolate P0 (the horizontal grid descriptors are already present in px_target)
	ip1=$(${fstliste} -izfst px_target -nomvar ">>" -col 3)
	ip2=$(${fstliste} -izfst px_target -nomvar ">>" -col 4)
	ip3=$(${fstliste} -izfst px_target -nomvar ">>" -col 5)

	${pgsm} -iment ./P0_source.fst -ozsrt px_target -i <<EOF
 SORTIE(STD,4000,R)
 GRILLE(TAPE2,${ip1},${ip2},${ip3})
 SETINTX(LINEAIR)
 HEURE(TOUT)
 CHAMP('P0')
 END
EOF

    else
	# Copy P0 as is
	${echo_cmd} "  desire(-1,['P0','>>','^^','^>'])" > edit.dir
	${editfst} -s ./file_datev.fst -d px_target -i edit.dir
    fi

#====================================================================
# Prepare px_source

    ${rm_cmd} -f px_source
    ${echo_cmd} "  desire(-1,['P0','>>','^^','^>','!!'])" > edit.dir
    ${editfst} -s ./file_datev.fst -d px_source -i edit.dir

#====================================================================
# Do the interpolation

    PXT_LEVELS_M=MOMENTUM
    PXT_LEVELS_T=MOMENTUM
    if [ "${__INTERPENSTRIALS_MODE__}" = stg2stg ];then
       PXT_LEVELS_T=THERMO
    fi

    if [ "$DO_HORIZ_INTERP" = yes ] ; then
	MOMENTUM_INTERP="LIN_UV"        # VECTOR interpolation
	THERMO_INTERP="LIN_TT LIN_HU"
    else
        MOMENTUM_INTERP="CUB_UU CUB_VV"	# SCALAR interpolation
	THERMO_INTERP="CUB_TT LIN_HU"
    fi

    rm -f target_MM target_TH
    ${pxs2pxt} \
	-s ./file_datev.fst \
	-pxs px_source \
	-pxt px_target \
	-d target_MM \
	-etiket ${CETIK} \
	-var ${MOMENTUM_INTERP} \
        -pxt_levels ${PXT_LEVELS_M} \
	-above 0.0

    ${pxs2pxt} \
	-s ./file_datev.fst \
	-pxs px_source \
	-pxt px_target \
	-d target_TH \
	-etiket ${CETIK} \
	-var ${THERMO_INTERP} \
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
	exit 1
    fi
    ${echo_cmd} " end" >> edit.dir
    ${editfst} -s grid_prototype -d target.fst  -i edit.dir

# Add surf and geophys fields  to the output file
    ${rm_cmd} -f sfc_geo.fst
#   Exclude GZ P0 and all the fields just interpolated by pxs2pxt
#   EXclude VT PX if present ( not interpolated)
#   Exclude !! >> ^^
    ${echo_cmd} " exclure(-1,['P0','GZ','TT','HU','ES','UU','VV'])" > edit.dir
    ${echo_cmd} " exclure(-1,['!!','>>','^^','HY'])" >> edit.dir
    ${echo_cmd} " exclure(-1,['VT','PX'])" >> edit.dir
    ${echo_cmd} " end" >> edit.dir
    
    ${editfst} -s ./file_datev.fst -d sfc_geo.fst -i edit.dir
#  Add sfc_geo.fst to target
    if [ "$DO_HORIZ_INTERP" = yes ] ; then

	${echo_cmd} " desire(-1,['>>','^^'])" > edit.dir
	${editfst} -s ./file_datev.fst -d sfc_geo.fst -i edit.dir

	rm -f sfc_geo_target.fst
	${echo_cmd} " desire(-1,['>>','^^'])" > edit.dir
	${editfst} -s target_prototype -d sfc_geo_target.fst -i edit.dir

	ip1=$(${fstliste} -izfst sfc_geo_target.fst -nomvar ">>" -col 3)
        ip2=$(${fstliste} -izfst sfc_geo_target.fst -nomvar ">>" -col 4)
	ip3=$(${fstliste} -izfst sfc_geo_target.fst -nomvar ">>" -col 5)

        ${pgsm} -iment sfc_geo.fst -ozsrt sfc_geo_target.fst -i<<EOF
 SORTIE(STD,4000,R)
 GRILLE(TAPE2,${ip1},${ip2},${ip3})
 SETINTX(LINEAIR)
 HEURE(TOUT)
 CHAMP(TOUT)
 END
EOF

	${echo_cmd} " exclure(-1,['>>','^^','!!'])" > edit.dir
	${editfst} -s sfc_geo_target.fst -d target.fst -i 0

    else
	${editfst} -s sfc_geo.fst -d target.fst -i 0
    fi

    ${echo_cmd} "    ***FIN DE ${__INTERPENSTRIALS_MODE__}***"

    ${echo_cmd} "   ***FIN de la datev " ${jdate}
done

## En mode global legacy, on doit faire ce traitement
## pour garantir la validation avec les executions precedentes
if [ "${INTERPENSTRIALS_INTERP2GAUSS:-no}" != no -a "$DO_HORIZ_INTERP" = no ]; then
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

fi ## Global legacy

else ## else relie au 'if [ "${INTERPENSTRIALS_DO_INTERPOLATION}" = yes ]'

    ln -s ${jfile} target.fst

fi ## Fin du else relie 'if [ "${INTERPENSTRIALS_DO_INTERPOLATION}" = yes ]'

# Final copy to the output directory
cp target.fst ${dir_out}/${file_out}

[ -f file_datev.fst ]        && ${rm_cmd} -f file_datev.fst
[ -f file_datev_int_gg.fst ] && ${rm_cmd} -f file_datev_int_gg.fst
[ -f target.fst ]            && ${rm_cmd} -f target.fst
[ -f target_MM.fst ]         && ${rm_cmd} -f target_MM.fst
[ -f target_TH.fst ]         && ${rm_cmd} -f target_TH.fst
[ -f px_source_momentum ]    && ${rm_cmd} -f px_source_momentum
[ -f px_source_thermo ]      && ${rm_cmd} -f px_source_thermo
[ -f px_target ]             && ${rm_cmd} -f px_target

${echo_cmd} "   ***FIN INTERPOLATION***"  
