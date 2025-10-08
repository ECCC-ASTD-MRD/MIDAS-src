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

    ## In normal mode, the tests '/Tests/var/EnVar/gdps' and
    ## '/Tests/var/EnVar/hrdps' run very fine with the operational
    ## number of iterations.
    typeset -r var_nitermax_expected_value=15
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

    ## In DEBUG mode, the tests '/Tests/var/EnVar/gdps' and
    ## '/Tests/var/EnVar/hrdps' are so long that we cannot run them
    ## with 70 iterations within the maximum wall clock time.
    typeset -r var_nitermax_expected_value=70
fi

## For the moment, the test '/Tests/letkf/glb_10km' is always run.
typeset -r letkf_glb10km_catchup_expected_value=2

typeset -r resources_file=${toplevel}/maestro/suites/midas_system_tests/resources/resources.def

## Rebuild the 'resources.def' file by calling 'set_resources_def.sh'
${toplevel}/set_resources_def.sh ${toplevel}

## First, check if there is one and only one definition of at least one of the
## variables in the resources file
found_several_definitions_of_the_same_variable=false
for variable in CHECK_RESULTS_CATCHUP CLEAN_UNITTEST_CATCHUP LETKF_GLB15KM_NUM_THREADS LETKF_GLB15KM_MEMORY LETKF_GLB10KM_CATCHUP VAR_NITERMAX; do
    if [ "$(grep -c "^${variable}=" ${resources_file})" -ne 1 ]; then
        echo "Found none or several definitions of the variable '${variable}' in '${resources_file}'" >&2
        found_several_definitions_of_the_same_variable=true
    fi
done

if [[ "${found_several_definitions_of_the_same_variable}" = true ]]; then
    echo "This is not expected so We will not modify the resources file" >&2
    exit 1
fi

## If the variables 'CHECK_RESULTS_CATCHUP', 'CLEAN_UNITTEST_CATCHUP',
## 'LETKF_GLB15KM_NUM_THREADS' and 'LETKF_GLB15KM_MEMORY' all have
## already the expected value, then we don't need to change the file.
typeset -r check_result_catchup=$(awk -F= '/^CHECK_RESULTS_CATCHUP=/ {print $2}' ${resources_file})
typeset -r clean_unittest_catchup=$(awk -F= '/^CLEAN_UNITTEST_CATCHUP=/ {print $2}' ${resources_file})
typeset -r letkf_glb15km_num_threads=$(awk -F= '/^LETKF_GLB15KM_NUM_THREADS=/ {print $2}' ${resources_file})
typeset -r letkf_glb15km_memory=$(awk -F= '/^LETKF_GLB15KM_MEMORY=/ {print $2}' ${resources_file})
typeset -r letkf_glb10km_catchup=$(awk -F= '/^LETKF_GLB10KM_CATCHUP=/ {print $2}' ${resources_file})
typeset -r var_nitermax=$(awk -F= '/^VAR_NITERMAX=/ {print $2}' ${resources_file})
if [ "${check_result_catchup}"      -eq "${check_result_catchup_expected_value}"      -a \
     "${clean_unittest_catchup}"    -eq "${clean_unittest_catchup_expected_value}"    -a \
     "${letkf_glb15km_num_threads}" -eq "${letkf_glb15km_num_threads_expected_value}" -a \
     "${letkf_glb15km_memory}"       =  "${letkf_glb15km_memory_expected_value}"      -a \
     "${letkf_glb10km_catchup}"     -eq "${letkf_glb10km_catchup_expected_value}"     -a \
     "${var_nitermax}"              -eq "${var_nitermax_expected_value}" ]; then
    ## All variables 'CHECK_RESULTS_CATCHUP',
    ## 'CLEAN_UNITTEST_CATCHUP', 'LETKF_GLB15KM_NUM_THREADS',
    ## 'LETKF_GLB15KM_MEMORY', 'COMPILER', 'LETKF_GLB10KM_CATCHUP' and
    ## 'VAR_NITERMAX' are all already set to ${check_result_catchup},
    ## ${clean_unittest_catchup}, ${letkf_glb15km_num_threads},
    ## ${letkf_glb15km_memory}, ${letkf_glb10km_catchup} and ${var_nitermax}
    ## respectively in '${resources_file}'.
    ## We will not modify the resources file
    exit
fi

if [ "${check_result_catchup}" -ne "${check_result_catchup_expected_value}" ]; then
    echo "Set 'CHECK_RESULTS_CATCHUP=${check_result_catchup_expected_value}' in ${resources_file}"
    sed -i "s/^CHECK_RESULTS_CATCHUP=[1-9]$/CHECK_RESULTS_CATCHUP=${check_result_catchup_expected_value}/" ${resources_file}
fi

if [ "${clean_unittest_catchup}" -ne "${clean_unittest_catchup_expected_value}" ]; then
    echo "Set 'CLEAN_UNITTEST_CATCHUP=${clean_unittest_catchup_expected_value}' in ${resources_file}"
    sed -i "s/^CLEAN_UNITTEST_CATCHUP=[1-9]$/CLEAN_UNITTEST_CATCHUP=${clean_unittest_catchup_expected_value}/" ${resources_file}
fi

if [ "${letkf_glb15km_num_threads}" -ne "${letkf_glb15km_num_threads_expected_value}" ]; then
    echo "Set 'LETKF_GLB15KM_NUM_THREADS=${letkf_glb15km_num_threads_expected_value}' in ${resources_file}"
    sed -i "s/^LETKF_GLB15KM_NUM_THREADS=[1-9]$/LETKF_GLB15KM_NUM_THREADS=${letkf_glb15km_num_threads_expected_value}/" ${resources_file}
fi

if [ "${letkf_glb15km_memory}" !=  "${letkf_glb15km_memory_expected_value}" ]; then
    echo "Set 'LETKF_GLB15KM_MEMORY=${letkf_glb15km_memory_expected_value}' in ${resources_file}"
    sed -i "s/^LETKF_GLB15KM_MEMORY=[0-9]\+[MGT]$/LETKF_GLB15KM_MEMORY=${letkf_glb15km_memory_expected_value}/" ${resources_file}
fi

if [ "${letkf_glb10km_catchup}" !=  "${letkf_glb10km_catchup_expected_value}" ]; then
    echo "Set 'LETKF_GLB10KM_CATCHUP=${letkf_glb10km_catchup_expected_value}' in ${resources_file}"
    sed -i "s/^LETKF_GLB10KM_CATCHUP=[^ ]\+$/LETKF_GLB10KM_CATCHUP=${letkf_glb10km_catchup_expected_value}/" ${resources_file}
fi

if [ "${var_nitermax}" !=  "${var_nitermax_expected_value}" ]; then
    echo "Set 'VAR_NITERMAX=${var_nitermax_expected_value}' in ${resources_file}"
    sed -i "s/^VAR_NITERMAX=[^ ]\+$/VAR_NITERMAX=${var_nitermax_expected_value}/" ${resources_file}
fi

#if [ "${}" !=  "${_expected_value}" ]; then
#    echo "Set '=${_expected_value}' in ${resources_file}"
#    sed -i "s/^=[^ ]\+$/=${_expected_value}/" ${resources_file}
#fi
