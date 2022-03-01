#!/bin/bash

set -e

resourcesDir=$(git rev-parse --show-toplevel)/maestro/suites/midas_system_tests/resources

# immediately return if resources.def already defined
if [ -f "${resourcesDir}/resources.def" ]; then
    exit
fi
    
# select which file to use for resources.def
if [ "${TRUE_HOST}" = eccc-ppp1 -o "${TRUE_HOST}" = eccc-ppp2 -o "${TRUE_HOST}" = hare -o "${TRUE_HOST}" = brooks ]; then
    resources_file=resources.def_ppp12_xc40
elif [ "${TRUE_HOST}" = eccc-gpsc1 -o "${TRUE_HOST}" = eccc-ppp3 -o "${TRUE_HOST}" = eccc-ppp4 -o "${TRUE_HOST}" = daley -o "${TRUE_HOST}" = banting -o "${TRUE_HOST}" = xc3 -o "${TRUE_HOST}" = xc4 ]; then
    resources_file=resources.def_ppp34_xc50
elif [ "${TRUE_HOST}" = ppp5 -o "${TRUE_HOST}" = ppp6 -o "${TRUE_HOST}" = robert -o "${TRUE_HOST}" = underhill ]; then
    resources_file=resources.def_hpcr-u2
else
    echo "Unknown TRUE_HOST: ${TRUE_HOST}"
    exit
fi

# set symbolic link
echo
echo "Setting symbolic link: resources.def ===> ${resources_file}"
echo
cd ${resourcesDir}
rm -f resources.def
ln -s ${resources_file} resources.def
cd -
