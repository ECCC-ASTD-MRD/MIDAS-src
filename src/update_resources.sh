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
    # and we want to use 2 OpenMP threads and 11G of memory for task '/Tests/letkf/glb_15km/UnitTest/run'
    typeset -r letkf_glb15km_num_threads_expected_value=4
    typeset -r letkf_glb15km_memory_expected_value=22G
else
    ## If neither of these environment variable above are set, we want
    ## to run the 'check' and 'clean' tasks in the 'UnitTest' module
    ## so we set 'CHECK_RESULTS_CATCHUP' and 'CHECK_RESULTS_CATCHUP'
    ## to '2' in
    ## 'maestro/suites/midas_system_tests/resources/resources.def'
    typeset -r check_result_catchup_expected_value=2
    typeset -r clean_unittest_catchup_expected_value=2
    # and we want to use 2 OpenMP threads and 11G of memory for task '/Tests/letkf/glb_15km/UnitTest/run'
    typeset -r letkf_glb15km_num_threads_expected_value=2
    typeset -r letkf_glb15km_memory_expected_value=11G
fi

typeset -r resources_file=${toplevel}/maestro/suites/midas_system_tests/resources/resources.def

## Rebuild the 'resources.def' file by calling 'set_resources_def.sh'
${toplevel}/set_resources_def.sh

## First, check if there is one and only one definition of the
## variable 'CHECK_RESULTS_CATCHUP' in the resources file
if [ "$(grep -c '^CHECK_RESULTS_CATCHUP=' ${resources_file})" -ne 1 ]; then
    echo "Found none or several definitions of the variable 'CHECK_RESULTS_CATCHUP' in '${resources_file}'"
    echo "We will not modify the resources file"
    exit
fi
# Same thing with variable 'CLEAN_UNITTEST_CATCHUP'
if [ "$(grep -c '^CLEAN_UNITTEST_CATCHUP=' ${resources_file})" -ne 1 ]; then
    echo "Found none or several definitions of the variable 'CLEAN_UNITTEST_CATCHUP' in '${resources_file}'"
    echo "We will not modify the resources file"
    exit
fi

echo "Set 'CHECK_RESULTS_CATCHUP=${check_result_catchup_expected_value}' in ${resources_file}"
sed -i "s/^CHECK_RESULTS_CATCHUP=[1-9]$/CHECK_RESULTS_CATCHUP=${check_result_catchup_expected_value}/" ${resources_file}

echo "Set 'CLEAN_UNITTEST_CATCHUP=${clean_unittest_catchup_expected_value}' in ${resources_file}"
sed -i "s/^CLEAN_UNITTEST_CATCHUP=[1-9]$/CLEAN_UNITTEST_CATCHUP=${clean_unittest_catchup_expected_value}/" ${resources_file}

echo "Set 'LETKF_GLB15KM_NUM_THREADS=${letkf_glb15km_num_threads_expected_value}' in ${resources_file}"
sed -i "s/^LETKF_GLB15KM_NUM_THREADS=[1-9]$/LETKF_GLB15KM_NUM_THREADS=${letkf_glb15km_num_threads_expected_value}/" ${resources_file}

echo "Set 'LETKF_GLB15KM_MEMORY=${letkf_glb15km_memory_expected_value}' in ${resources_file}"
sed -i "s/^LETKF_GLB15KM_MEMORY=[0-9]\+[MGT]$/LETKF_GLB15KM_MEMORY=${letkf_glb15km_memory_expected_value}/" ${resources_file}
