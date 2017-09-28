#!/bin/bash

set -e

OAVAR_SUITE_LAUNCH_DIRECTORY=$(dirname $(true_path $0))
SEQ_MAESTRO_SHORTCUT=${SEQ_MAESTRO_SHORTCUT:-". ssmuse-sh -d eccc/cmo/isst/maestro/1.5.1-rc11"}

which clone_suite 1>/dev/null 2>&1 || . ssmuse-sh -d eccc/cmd/cmda/maestro/dev/2.10
which maestro     1>/dev/null 2>&1 || ${SEQ_MAESTRO_SHORTCUT}

while [ -z "${REPONSE}" ]; do
    echo "Please choose a suite name to install of maestro oavar suite of tests"
    echo "It will be installed under ${HOME}/.suites"
    read REPONSE
    if [ -z "${REPONSE}" ]; then
	echo "Nothing was entered... so we use the default 'oavar_system_tests_science'."
	REPONSE=oavar_system_tests_science
    fi
    if [ -a ~/.suites/${REPONSE} ]; then

	trpwd=$(true_path ${OAVAR_SUITE_LAUNCH_DIRECTORY})
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
	    echo "   OAVAR_SUITE_LAUNCH_DIRECTORY=${OAVAR_SUITE_LAUNCH_DIRECTORY}"
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

OAVAR_TESTS_SUITE=${REPONSE}

cd ${OAVAR_SUITE_LAUNCH_DIRECTORY}

if [[ ${OAVAR_TESTS_SUITE} = */* ]]; then
    mkdir -p $(dirname ${OAVAR_TESTS_SUITE})
fi

ln -s $PWD ~/.suites/${OAVAR_TESTS_SUITE}
MAKE_LINKS_START_DATE=$(date +%Y%m%d000000) ~erv000/maestro_utilities/make_links/make_links ${OAVAR_TESTS_SUITE}

echo "ABS_DIR=${COMPILEDIR_OAVAR_MAIN:-$(dirname $(dirname $(dirname ${PWD})))/compiledir}/oavar_abs" > abs.dot
echo "OAVAR_version=\$(cd ${PWD}/..; git describe --abbrev=7 --always --dirty=_M 2>/dev/null || ssh eccc-ppp1 'cd ${PWD}/..; git describe --abbrev=7 --always --dirty=_M' 2>/dev/null || echo unkown revision)" >> abs.dot

## Ajouter la creation pour chaque usager de repertoires de reference pour les tests
##    test_results

export SEQ_EXP_HOME=~/.suites/${OAVAR_TESTS_SUITE}
xflow &
