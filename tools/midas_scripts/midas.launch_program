#!/bin/bash

set -e

running_mode=run
while [ $# -ne 0 ]; do
    if [ "${1}" = --ddt ]; then
        running_mode=ddt
        which ddt 1>/dev/null 2>&1 || . ssmuse-sh -x main/opt/forge/21.1.3
    elif [ "${1}" = --gdb ]; then
        running_mode=gdb
    elif [[ "${1}" = --vtune* ]]; then
        ## Can be '--vtune', '--vtune=${vtune options}'
        ## See https://www.intel.com/content/www/us/en/docs/vtune-profiler/user-guide/2024-0/collect.html
        running_mode=${1}
        ## Load the full compiler suite
        if ! which vtune 1>/dev/null 2>&1; then
            if [ -n "${VTUNE_SSM}" ]; then
                . ssmuse-sh -x ${VTUNE_SSM}
            else
                echo "${0}: aborts because it cannot find 'vtune' and 'VTUNE_SSM' environment variable does not exist"
                exit 1
            fi
        fi
    elif [ "${1}" = -h -o "${1}" = --help ]; then
        echo "${0} [--ddt] [--gdb] [--vtune] \${path_to_program}"
        exit
    else
        run_pgm=${1}
    fi
    shift
done

if [ -z "${run_pgm}" ]; then
    echo 'You must provide a program to launch' >&2
    exit 1
fi

if [[ ${run_pgm} != */* ]]; then
    status=0
    which ${run_pgm} 1> /dev/null 2>&1 || status=1
    if [ "${status}" -ne 0 ]; then
        if [ -x ${run_pgm} ]; then
            run_pgm=./${run_pgm}
        fi
    fi
fi

if [ ! -f "${run_pgm}" ]; then
    echo "The program specified '${run_pgm}' does not exist!" >&2
    exit 1
fi

if [ "${running_mode}" = ddt ]; then
    true ## no script to prepare when running in DDT mode
elif [ "${running_mode}" = gdb ]; then
    cat > launch_cmd <<EOFLAUNCH
#!/bin/bash

set -ex

gdb -ex run -ex where -ex quit ${run_pgm}

EOFLAUNCH
    chmod +x launch_cmd

elif [[ "${running_mode}" = --vtune* ]]; then
     cat > launch_cmd <<EOFLAUNCH
#!/bin/bash

set -ex

## this will generate directories named 'vtune.\${hostname}'
EOFLAUNCH

    if [ "${running_mode}" = --vtune ]; then
        vtune_collect_mode=hotspots
    else
        ## The script 'launch_program.sh' can be called with '--vtune=hotspots' or '--vtune=hpc-performance'
        vtune_collect_mode=$(echo ${running_mode} | sed 's/^--vtune=//')
    fi

    cat >> launch_cmd <<EOFLAUNCH
## this will generate directories named 'vtune.\${hostname}'
vtune -collect ${vtune_collect_mode} –r vtune -- ${run_pgm}

EOFLAUNCH
    chmod +x launch_cmd

elif [ "${running_mode}" = run ]; then
    cat > launch_cmd <<EOFLAUNCH
#!/bin/bash

set -e

/usr/bin/time --format="Total Memory %M Kb" ${run_pgm}

EOFLAUNCH
    chmod +x launch_cmd

else
    echo "The 'running_mode' can only be 'run', 'ddt', 'gdb', 'vtune' and not '${running_mode}'!" >&2
    exit 1
fi

## This may be needed for 'mpiscript'
export RUN_PGM=${PWD}/launch_cmd

echo "$(date +%Y%m%d:%H:%M:%S.%N): Copy of observations starts"
[ -d obs ] && rm -rf obs
sscp -r obs.original obs

echo "$(date +%Y%m%d:%H:%M:%S.%N): Cleaning working directory"

for file in *; do
    if ./isFileInList.sh ${file} not in obs $(cat localFilesNotToErase); then
        rm -r ${file}
    fi
done

if [ "${running_mode}" = ddt ]; then
    echo "$(date +%Y%m%d:%H:%M:%S.%N): Launching ${run_pgm} with DDT"
    ## DDT cannot run a script, it must launch the program itself
    ddt mpirun -n $((SEQ_NPEX*SEQ_NPEY)) ${run_pgm}
else
    if [ "${running_mode}" = --vtune ]; then
        info_mode=" with vtune"
    elif [ "${running_mode}" = gdb ]; then
        info_mode=" with gdb"
    fi

    echo "$(date +%Y%m%d:%H:%M:%S.%N): Launching ${run_pgm}${info_mode}"
    echo Using r.run_in_parallel from $(which r.run_in_parallel)

    SECONDS=0
    r.run_in_parallel -pgm ${LAUNCH_CMD} -npex ${SEQ_NPEX} -npey ${SEQ_NPEY} -processorder -verbose -tag -nocleanup -tmpdir ${PWD}/mpitmpdir ${run_in_parallel_extra_args} -args ${UnitTest_run_pgm_args}
    echo RUNTIME=${SECONDS}
fi

echo "End of ${pgmpath} at $(date +%Y%m%d:%H:%M:%S.%N)"

if [[ "${running_mode}" = --vtune* ]]; then
    echo "You can now vizualize the profiling data using the commands"
    if ! which vtune-gui 1>/dev/null 2>&1; then
        echo "    . ssmuse-sh -x ${VTUNE_SSM}"
    fi
    echo "    vtune-gui ./vtune.*"
    echo "Select only one of the 'vtune.*' directories"
fi
