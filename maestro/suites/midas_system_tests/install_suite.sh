#!/bin/bash

set -e

__toplevel=$(git rev-parse --show-toplevel)

MIDAS_SUITE_LAUNCH_DIRECTORY=${__toplevel}/maestro/suites/midas_system_tests

# If it does not yet exist, set the resources.def file, which depends
# on '${TRUE_HOST}'
if [ ! -f ${MIDAS_SUITE_LAUNCH_DIRECTORY}/resources/resources.def ]; then
    ## Run this command as a script since we only need the file
    ## 'maestro/suites/midas_system_tests/resources/resources.def' to
    ## be generated.
    status=0
    cmake_listing=${TMPDIR}/cmake_install_suite_resources.listing
    ${__toplevel}/src/config.dot.sh --build-id resources --fresh --no-cd-build --no-show-instructions > ${cmake_listing} 2>&1 || status=1
    if [[ "${status}" -ne 0 ]]; then
        echo "install_suite.sh: There was an error while updating the resources file (see the listing ${cmake_listing})" >&2
        exit 1
    else
        rm ${cmake_listing}
    fi
fi
## Initialize the hosts list for the test suite
. ${MIDAS_SUITE_LAUNCH_DIRECTORY}/set_machine_list.dot

which maestro 1>/dev/null 2>&1 || {
    SEQ_MAESTRO_SHORTCUT=${SEQ_MAESTRO_SHORTCUT:-". ssmuse-sh -d eccc/cmo/isst/maestro/20250808"}
    ${SEQ_MAESTRO_SHORTCUT}
}

__rmnlib_version__=20260202
if [[ "${ORDENV_PLAT}" = rhel-8-icelake-64 ]]; then
    __cmdiutils_version__=2.9
elif [[ "${ORDENV_PLAT}" = rhel-9-graniterapids-64 ]]; then
    __cmdiutils_version__=2.12
else
    echo "The platform '${ORDENV_PLAT}' is not supported.  Only 'rhel-8-icelake-64' and 'rhel-9-graniterapids-64' are!" >&2
    exit 1
fi
which r.date      1>/dev/null 2>&1 || . r.load.dot    eccc/mrd/rpn/utils/${__rmnlib_version__}
which clone_suite 1>/dev/null 2>&1 || . ssmuse-sh -d eccc/cmd/cmdi/utils/${__cmdiutils_version__}

unset __rmnlib_version__ __cmdiutils_version__

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

                    echo "Do you want to erase everything and reinstall the suite? (Y,y,yes,YES,...)"
	            read REPONSE_YES_OR_NO
	            if [[ "${REPONSE_YES_OR_NO}" = [yY] ]] || [[ "${REPONSE_YES_OR_NO}" = [yY][eE][sS] || \
		        "${REPONSE_YES_OR_NO}" = [oO] ]] || [[ "${REPONSE_YES_OR_NO}" = [oO][uU][iI] ]]; then
                        echo "The suite will be reinstalled"
                        break
                    else
		        echo "Then do:"
		        echo "   xflow -exp ~/.suites/${REPONSE}"
		        exit
                    fi
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

if [ ! -f ${MIDAS_SUITE_LAUNCH_DIRECTORY}/abs.dot ]; then
    ## Run this command as a script since we only need the file
    ## 'maestro/suites/midas_system_tests/abs.dot' to
    ## be generated.
    status=0
    cmake_listing=${TMPDIR}/cmake_install_suite_absdot.listing
    ${__toplevel}/src/config.dot.sh --build-id resources --fresh --no-cd-build --no-show-instructions > ${cmake_listing} 2>&1 || status=1
    if [[ "${status}" -ne 0 ]]; then
        echo "install_suite.sh: There was an error while updating the 'abs.dot' file (see the listing ${cmake_listing})" >&2
        exit 1
    else
        rm ${cmake_listing}
    fi
fi

if [ -n "${MIDAS_ABS}" ]; then
    . ./abs.dot
    mkdir -p ${ABS_DIR}
    cp ${MIDAS_ABS}/midas-*-${MIDAS_version}.Abs ${ABS_DIR}
fi

[ -L ~/.suites/${MIDAS_TESTS_SUITE} ] && rm ~/.suites/${MIDAS_TESTS_SUITE}
ln -s $PWD ~/.suites/${MIDAS_TESTS_SUITE}

## if 'MIDAS_MAKE_LINKS_MACHINE_LIST' exists then use that list
if [ -n "${MIDAS_MAKE_LINKS_MACHINE_LIST}" ]; then
    export MAKE_LINKS_MACHINE_LIST=${MIDAS_MAKE_LINKS_MACHINE_LIST}
fi

export MAKE_LINKS_START_DATE=$(date +%Y%m%d000000)
make_links ${MIDAS_TESTS_SUITE}

## Ajouter la creation pour chaque usager de repertoires de reference pour les tests
##    test_results

export SEQ_EXP_HOME=~/.suites/${MIDAS_TESTS_SUITE}
if tty -s && [ "${MIDAS_INSTALL_SUITE_START_FLOW:-yes}" = yes ]; then
    xflow &
fi

if [ "${MIDAS_SUITE_RUN}" = run ]; then
    maestro -s submit -n /Tests -d ${MAKE_LINKS_START_DATE}
fi
