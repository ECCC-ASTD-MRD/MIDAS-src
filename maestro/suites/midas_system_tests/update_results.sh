#!/bin/ksh

set -e

UnitTest_name=${1:-all}
version=${2}

SEQ_EXP_HOME=${SEQ_EXP_HOME:-$(dirname $(true_path ${0}))}
export SEQ_EXP_HOME

toplevel=$(git rev-parse --show-toplevel)

SEQ_MAESTRO_SHORTCUT=${SEQ_MAESTRO_SHORTCUT:-". ssmuse-sh -d eccc/cmo/isst/maestro/1.8.0-beta"}
which getdef 1>/dev/null 2>&1 || ${SEQ_MAESTRO_SHORTCUT}

${toplevel}/set_resources_def.sh

if [ -z "${version}" ]; then
    version=$(${toplevel}/midas.version.sh)
fi

if [ "${version}" = 'unknown revision' ]; then
    echo "Cannot determine the version (version='${version}')"
    exit 1
fi

if [ "${UnitTest_name}" = all ]; then
    echo "Call '${0}' on each test"
    export __oavar_tests_update_results_do_commit__=no
    for unittest in $(grep SUBMITS ${SEQ_EXP_HOME}/modules/Tests/flow.xml | grep -v '<SUBMITS sub_name="UnitTest"/>' | awk -F\" '{print $2}'); do
	${0} ${unittest} ${version}
    done
    git --git-dir=${toplevel}/.git --work-tree=${level} commit -F - <<EOF
Update the results for all Unit Tests with version=${version}
EOF
else
    __oavar_tests_update_results_do_commit__=${__oavar_tests_update_results_do_commit__:-yes}

    echo "Call '${0}' on '${UnitTest_name}' with version=${version} and SEQ_EXP_HOME=${SEQ_EXP_HOME}"

    cfg=modules/Tests/${UnitTest_name}/container.cfg
    testname=$(awk -F= '/UnitTest_name=/ {print $2}' ${SEQ_EXP_HOME}/${cfg})
    [ -z "${testname}" ] && testname=${UnitTest_name}
    FRONTEND=$(getdef --exp ${SEQ_EXP_HOME} resources FRONTEND)
    ENVAR_output=${SEQ_EXP_HOME}/hub/${FRONTEND}/${testname}
    if [ ! -d "${ENVAR_output}" ]; then
	echo "The source directory '${ENVAR_output}' is not accessible!"
	exit 1
    fi

    eval $(grep ENVAR_DATE ${SEQ_EXP_HOME}/${cfg})
    if [ -z "${ENVAR_DATE}" ]; then
	echo "Cannot get the date in '${SEQ_EXP_HOME}/${cfg}'"
	exit 1
    fi

    output_results=$(dirname $(true_path ${SEQ_EXP_HOME}/hub/${FRONTEND}/..))/UnitTests/oavar_UnitTests/${UnitTest_name}/${version}
    if [ ! -d "${output_results}" ]; then

	mkdir -p ${output_results}
	cd ${output_results}

	if [ -f ${ENVAR_output}/gridpt/anal/hyb/${ENVAR_DATE}_000_nosfc ]; then
	    ln -si ${ENVAR_output}/gridpt/anal/hyb/${ENVAR_DATE}_000_nosfc anal.anl.model.${ENVAR_DATE}_000_nosfc
	    cmcarc -f ${ENVAR_DATE}.ca --dereference -a anal.anl.model.${ENVAR_DATE}_000_nosfc
	    rm anal.anl.model.${ENVAR_DATE}_000_nosfc
	fi

	if [ -f ${ENVAR_output}/gridpt/anal/hyb/${ENVAR_DATE}_inc_lo ]; then
	    ln -si ${ENVAR_output}/gridpt/anal/hyb/${ENVAR_DATE}_inc_lo anal.inc.lo.model.${ENVAR_DATE}_000
	    cmcarc -f ${ENVAR_DATE}.ca --dereference -a anal.inc.lo.model.${ENVAR_DATE}_000
	    rm anal.inc.lo.model.${ENVAR_DATE}_000
	fi

	for file in ${ENVAR_output}/gridpt/anal/hyb/${ENVAR_DATE}*_inc_hi; do
	    if [ -f "${file}" ]; then
		ln -si ${file} anal.inc.hi.model.$(basename $file)
		cmcarc -f ${ENVAR_DATE}.ca --dereference -a anal.inc.hi.model.$(basename $file)
		rm anal.inc.hi.model.$(basename $file)
	    fi
	done

	for file in ${ENVAR_output}/banco/postalt/${ENVAR_DATE}*; do
	    if [ -f "${file}" ]; then
		ln -si ${file} banco.postalt.$(basename ${file})
		cmcarc -f ${ENVAR_DATE}.ca --dereference -a banco.postalt.$(basename ${file})
		rm banco.postalt.$(basename ${file})
	    fi
	done

	chmod -wx ${ENVAR_DATE}.ca

	echo "Update UnitTest configuration to compare with these results"
	echo "UnitTest_results=${output_results}" >> ${SEQ_EXP_HOME}/${cfg}
	git --git-dir=${toplevel}/.git --work-tree=${toplevel} add "*/${cfg}"
	if [ "${__oavar_tests_update_results_do_commit__}" = yes ]; then
	    git --git-dir=${toplevel}/.git --work-tree=${toplevel} commit -F - <<EOF
Update the results for unittest '${Unittest_name}'
EOF
	fi
    else
	echo "The directory 'output_results=${output_results}' already exists!"
	cd ${output_results}
    fi

    echo "Collect the listings for '${UnitTest_name}'"
    if [ -f listings_${version}.ca ]; then
	echo "The file 'listings_${version}.ca' already exists so we do not update them"
	exit 0
    fi

    tasks="get check envar/getObservations envar/getTrials envar/getStats envar/getEnsTrials envar/Transfer envar/Transfer_anlm envar/interpEnsTrials envar/VAR envar/wait4increments envar/AddAnalInc envar/postalt envar/getPostalt"

    for task in ${tasks}; do
	nodename=/Tests/${UnitTest_name}/UnitTest/${task}
	filename=$(echo ${nodename} | cut -c2- | sed 's!/!\.!g')
	nodelister -n ${nodename} > ${filename}
	status=0
	cat <<EOF | diff ${filename} - 1> /dev/null 2>&1 || status=1
${nodename} listing not available
EOF
	if [ "${status}" -eq 0 ]; then
	    echo "The listing for '${nodename}' does not exist"
	    rm ${filename}
	    continue
	fi
	cmcarc -f listings_${version}.ca -a ${filename}
	rm ${filename}
    done
    ## cmcarc -f listings_${version}.ca -t -v

    chmod -wx listings_${version}.ca
fi ## Fin du 'else' relie au 'if [ "${UnitTest_name}" = all ]'
