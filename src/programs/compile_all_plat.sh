#!/bin/bash

set -e

MIDAS_ABS=${1}
codedir=${2:-${PWD}}

host=brooks
echo "Launching compilation on '${host}' for platform 'sles-11-broadwell-64-xc40'"
echo "cd ${codedir}; yes '' | ./compile_all.sh" | ssh brooks bash --login > compiling.listing.${host} &

echo "Launching compilation on '${TRUE_HOST}' for platform '${ORDENV_PLAT}'"
cd ${codedir}
yes '' | ./compile_all.sh || true

## wait for compilation on 'brooks' to finish
wait
cat compiling.listing.${host}
rm compiling.listing.${host}

echo "Checking if all programs have been compiled on '${TRUE_HOST}' for platform '${ORDENV_PLAT}'"
./check_if_all_programs_compiled.sh ${ORDENV_PLAT}            ${MIDAS_ABS}
echo "Checking if all programs have been compiled on '${host}' for platform 'sles-11-broadwell-64-xc40'"
./check_if_all_programs_compiled.sh sles-11-broadwell-64-xc40 ${MIDAS_ABS}
