#!/bin/bash

set -euo pipefail

typeset -r toplevel=$(git rev-parse --show-toplevel)

## Check if any of 'MIDAS_COMPILE_ADD_DEBUG_OPTIONS' or
## 'MIDAS_COMPILE_CODECOVERAGE_DATAPATH' is set.
## If not set, just do nothing!
if [ "${MIDAS_COMPILE_ADD_DEBUG_OPTIONS:-}" != yes -a -z "${MIDAS_COMPILE_CODECOVERAGE_DATAPATH:-}" ]; then
    exit
fi

## If any of these variables is set, we avoid running the 'check'
## and 'clean' tasks in the 'UnitTest' module by setting
## 'CHECK_RESULTS_CATCHUP' and 'CHECK_RESULTS_CATCHUP' to '9'

typeset -r resources_file=${toplevel}/maestro/suites/midas_system_tests/resources/resources.def
## If the resources file does not exist, then the user has not yet
## run 'maestro/suites/midas_system_tests/install_suite.sh' so we
## just cannot proceed furter here.
if [ ! -f "${resources_file}" ]; then
    exit
fi

## First, check if there is one and only one definition of the
## variable 'CHECK_RESULTS_CATCHUP' in the resources file
if [ "$(grep -c '^CHECK_RESULTS_CATCHUP=' ${resources_file})" -ne 1 ]; then
    echo "Found none or several definitions of the variable 'CHECK_RESULTS_CATCHUP' in '${resources_file}'"
    echo "We will not modify the resources file"
    exit
fi
## Same thing with variable 'CLEAN_UNITTEST_CATCHUP'
if [ "$(grep -c '^CLEAN_UNITTEST_CATCHUP=' ${resources_file})" -ne 1 ]; then
    echo "Found none or several definitions of the variable 'CLEAN_UNITTEST_CATCHUP' in '${resources_file}'"
    echo "We will not modify the resources file"
    exit
fi

## If the variables 'CHECK_RESULTS_CATCHUP' and
## 'CLEAN_UNITTEST_CATCHUP' both have already a value of '9', then
## we don't need to change the file
typeset -r check_result_catchup=$(awk -F= '/^CHECK_RESULTS_CATCHUP=/ {print $2}' ${resources_file})
typeset -r clean_unittest_catchup=$(awk -F= '/^CLEAN_UNITTEST_CATCHUP=/ {print $2}' ${resources_file})
if [ "${check_result_catchup}"  -eq 9 -a "${clean_unittest_catchup}" -eq 9 ]; then
    echo "Both variables 'CHECK_RESULTS_CATCHUP' and 'CLEAN_UNITTEST_CATCHUP' are already set to 9 in '${resources_file}'"
    echo "We will not modify the resources file"
    exit
fi

## If 'CHECK_RESULTS_CATCHUP' is not set to 9, then update it
if [ "${check_result_catchup}" -ne 9 ]; then
    echo "Set 'CHECK_RESULTS_CATCHUP=9' in ${resources_file}"
    sed -i 's/^CHECK_RESULTS_CATCHUP=[1-8]$/CHECK_RESULTS_CATCHUP=9/' ${resources_file}
fi

## If 'CLEAN_UNITTEST_CATCHUP' is not set to 9, then update it
if [ "${clean_unittest_catchup}" -ne 9 ]; then
    echo "Set 'CLEAN_UNITTEST_CATCHUP=9' in ${resources_file}"
    sed -i 's/^CLEAN_UNITTEST_CATCHUP=[1-8]$/CLEAN_UNITTEST_CATCHUP=9/' ${resources_file}
fi
