#!/bin/bash

set -e

MIDAS_SUITE_LAUNCH_DIRECTORY=$(dirname $(true_path $0))
SEQ_MAESTRO_SHORTCUT=${SEQ_MAESTRO_SHORTCUT:-". ssmuse-sh -d eccc/cmo/isst/maestro/1.5.1-rc21"}

which clone_suite 1>/dev/null 2>&1 || . ssmuse-sh -d eccc/cmd/cmdi/utils/2.0
which maestro     1>/dev/null 2>&1 || ${SEQ_MAESTRO_SHORTCUT}
which r.date      1>/dev/null 2>&1 || . ssmuse-sh -d eccc/mrd/rpn/utils/16.2.3

DEFAULT_SUITE_NAME=midas-$(git rev-parse --abbrev-ref HEAD | cut -d- -f1)

if [ $# -ne 0 ]; then
    MIDAS_TESTS_SUITE=$1
    MIDAS_SUITE_RUN=$2
    MIDAS_ABS=$3
else
    while [ -z "${REPONSE}" ]; do
        echo "Please choose a suite name to install of maestro midas suite of tests (default: ${DEFAULT_SUITE_NAME})"
        echo "It will be installed under ${HOME}/.suites"
        read REPONSE
        if [ -z "${REPONSE}" ]; then
	    echo "Nothing was entered... so we use the default '${DEFAULT_SUITE_NAME}'."
	    REPONSE=${DEFAULT_SUITE_NAME}
        fi
        if [ -a ~/.suites/${REPONSE} ]; then

	    trpwd=$(true_path ${MIDAS_SUITE_LAUNCH_DIRECTORY})
	    trsuite=$(true_path ~/.suites/${REPONSE})
	    if [ "${trpwd}" = "${trsuite}" ]; then
	        echo ' '
	        echo "The suite '~/.suites/${REPONSE}' is already existing."
	        echo "Do you want to use that suite? (Y,y,yes,YES,...)"
	        read REPONSE_YES_OR_NO
	        if [[ "${REPONSE_YES_OR_NO}" = [yY] ]] || [[ "${REPONSE_YES_OR_NO}" = [yY][eE][sS] || \
		    "${REPONSE_YES_OR_NO}" = [oO] ]] || [[ "${REPONSE_YES_OR_NO}" = [oO][uU][iI] ]]; then
		    echo "Then do:"
		    echo "   cd ~/.suites/${REPONSE}"
		    echo "   xflow"
		    exit
	        else
		    echo "Then use another name..."
		    echo
		    REPONSE=
		    continue
	        fi
	    else
	        echo ' '
	        echo "The suite '~/.suites/${REPONSE}' already exists but does not point to this directory"
	        echo "   MIDAS_SUITE_LAUNCH_DIRECTORY=${MIDAS_SUITE_LAUNCH_DIRECTORY}"
	        echo "   ~/.suites/${REPONSE} -> ${trsuite}"
	        echo "Do you want to override this suite?"
	        read REPONSE_YES_OR_NO
	        if [[ "${REPONSE_YES_OR_NO}" = [yY] ]] || [[ "${REPONSE_YES_OR_NO}" = [yY][eE][sS] || \
		    "${REPONSE_YES_OR_NO}" = [oO] ]] || [[ "${REPONSE_YES_OR_NO}" = [oO][uU][iI] ]]; then
		    rm ~/.suites/${REPONSE}
		    break
	        else
		    echo "Then use another name..."
		    echo
		    REPONSE=
		    continue
	        fi
	    fi

	    echo "Do you want to use that suite? (Y,y,yes,YES,...)"
	    read REPONSE_YES_OR_NO
	    if [[ "${REPONSE_YES_OR_NO}" = [yY] ]] || [[ "${REPONSE_YES_OR_NO}" = [yY][eE][sS] || \
                "${REPONSE_YES_OR_NO}" = [oO] ]] || [[ "${REPONSE_YES_OR_NO}" = [oO][uU][iI] ]]; then
	        echo "Then do:"
	        echo "   cd ~/.suites/${REPONSE}"
	        echo "   xflow"
	        exit
	    else
	        echo "Then use another name..."
	        echo
	        REPONSE=
	        continue
	    fi
        else
	    echo "The suite '~/.suites/${REPONSE}' will be created.  Please confirm (Y,y,yes,YES,...)"
	    read REPONSE_YES_OR_NO
	    if [[ "${REPONSE_YES_OR_NO}" = [yY] ]] || [[ "${REPONSE_YES_OR_NO}" = [yY][eE][sS] || \
                "${REPONSE_YES_OR_NO}" = [oO] ]] || [[ "${REPONSE_YES_OR_NO}" = [oO][uU][iI] ]]; then
	        break
	    else
	        echo
	        REPONSE=
	        REPONSE_YES_OR_NO=
	        continue
	    fi
        fi
    done

    MIDAS_TESTS_SUITE=${REPONSE}
fi  ## Fin du 'else' relie au 'if [ $# -ne 0 ]'

cd ${MIDAS_SUITE_LAUNCH_DIRECTORY}

echo "Installing suite ~/.suites/${MIDAS_TESTS_SUITE}"

if [[ ${MIDAS_TESTS_SUITE} = */* ]]; then
    mkdir -p $(dirname ${MIDAS_TESTS_SUITE})
fi

[ -L ~/.suites/${MIDAS_TESTS_SUITE} ] && rm ~/.suites/${MIDAS_TESTS_SUITE}
ln -s $PWD ~/.suites/${MIDAS_TESTS_SUITE}
export MAKE_LINKS_START_DATE=$(date +%Y%m%d000000)
make_links ${MIDAS_TESTS_SUITE}

echo "ABS_DIR=${COMPILEDIR_MIDAS_MAIN:-$(dirname $(dirname $(dirname ${PWD})))/compiledir}/midas_abs" > abs.dot
echo "MIDAS_version=\$(cd ${PWD}/..; git describe --abbrev=7 --always --dirty=_M 2>/dev/null || ssh eccc-ppp1 'cd ${PWD}/..; git describe --abbrev=7 --always --dirty=_M' 2>/dev/null || echo unkown revision)" >> abs.dot

## Ajouter la creation pour chaque usager de repertoires de reference pour les tests
##    test_results

export SEQ_EXP_HOME=~/.suites/${MIDAS_TESTS_SUITE}
if tty -s && [ "${MIDAS_INSTALL_SUITE_START_FLOW:-yes}" = yes ]; then
    xflow &
fi

if [ -n "${MIDAS_ABS}" ]; then
    . ./abs.dot
    mkdir -p ${ABS_DIR}
    cp ${MIDAS_ABS}/midas-*-${MIDAS_version}.Abs ${ABS_DIR}
fi

if [ "${MIDAS_SUITE_RUN}" = run ]; then
    maestro -s submit -n /Tests -d ${MAKE_LINKS_START_DATE}
fi
