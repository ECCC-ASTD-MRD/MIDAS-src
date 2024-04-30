#!/bin/bash

set -euo pipefail

typeset -r toplevel=${1:-$(git rev-parse --show-toplevel)}

## Check if any of 'MIDAS_COMPILE_ADD_DEBUG_OPTIONS' or
## 'MIDAS_COMPILE_CODECOVERAGE_DATAPATH' is set.
if [ "${MIDAS_COMPILE_ADD_DEBUG_OPTIONS:-}" = yes -o -n "${MIDAS_COMPILE_CODECOVERAGE_DATAPATH:-}" ]; then
    ## If any of these environment variables above is set, we want
    ## avoid running the 'check' and 'clean' tasks in the 'UnitTest'
    ## module by setting 'CHECK_RESULTS_CATCHUP' and
    ## 'CHECK_RESULTS_CATCHUP' to '9' in
    ## 'maestro/suites/midas_system_tests/resources/resources.def'
    typeset -r check_result_catchup_expected_value=9
    typeset -r clean_unittest_catchup_expected_value=9
else
    ## If neither of these environment variable above are set, we want
    ## to run the 'check' and 'clean' tasks in the 'UnitTest' module
    ## so we set 'CHECK_RESULTS_CATCHUP' and 'CHECK_RESULTS_CATCHUP'
    ## to '2' in
    ## 'maestro/suites/midas_system_tests/resources/resources.def'
    typeset -r check_result_catchup_expected_value=2
    typeset -r clean_unittest_catchup_expected_value=2
fi


typeset -r resources_file=${toplevel}/maestro/suites/midas_system_tests/resources/resources.def
## If the resources file does not exist, then the user has not yet
## run 'set_resources_def.sh' so we will call it from here
if [ ! -f "${resources_file}" ]; then
    ${toplevel}/set_resources_def.sh
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
## 'CLEAN_UNITTEST_CATCHUP' both have already the expected value, then
## we don't need to change the file
typeset -r check_result_catchup=$(awk -F= '/^CHECK_RESULTS_CATCHUP=/ {print $2}' ${resources_file})
typeset -r clean_unittest_catchup=$(awk -F= '/^CLEAN_UNITTEST_CATCHUP=/ {print $2}' ${resources_file})
if [ "${check_result_catchup}"   -eq "${check_result_catchup_expected_value}" -a \
     "${clean_unittest_catchup}" -eq "${clean_unittest_catchup_expected_value}" ]; then
    ## echo "Both variables 'CHECK_RESULTS_CATCHUP' and 'CLEAN_UNITTEST_CATCHUP' are already set"
    ## echo "to ${check_result_catchup_expected_value} and ${clean_unittest_catchup_expected_value}"
    ## echo "respectively in '${resources_file}'."
    ## echo "We will not modify the resources file"
    exit
fi

## If 'CHECK_RESULTS_CATCHUP' is not set to 9, then update it
if [ "${check_result_catchup}" -ne "${check_result_catchup_expected_value}" ]; then
    echo "Set 'CHECK_RESULTS_CATCHUP=${check_result_catchup_expected_value}' in ${resources_file}"
    sed -i "s/^CHECK_RESULTS_CATCHUP=[1-9]$/CHECK_RESULTS_CATCHUP=${check_result_catchup_expected_value}/" ${resources_file}
fi

## If 'CLEAN_UNITTEST_CATCHUP' is not set to 9, then update it
if [ "${clean_unittest_catchup}" -ne "${clean_unittest_catchup_expected_value}" ]; then
    echo "Set 'CLEAN_UNITTEST_CATCHUP=${clean_unittest_catchup_expected_value}' in ${resources_file}"
    sed -i "s/^CLEAN_UNITTEST_CATCHUP=[1-9]$/CLEAN_UNITTEST_CATCHUP=${clean_unittest_catchup_expected_value}/" ${resources_file}
fi
