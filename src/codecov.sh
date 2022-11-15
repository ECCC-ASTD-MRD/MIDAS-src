#!/bin/bash

## This script is using the variable
## 'MIDAS_COMPILE_ADD_CODECOVERAGE_OPTIONS' as a directory to find the
## code coverage information.

## Those steps are inspired from this discussion:
##    https://github.com/Goddard-Fortran-Ecosystem/pFUnit/issues/309#issuecomment-874886020

set -euo pipefail

typeset -r scriptName=codecov.sh
typeset -r defautCovDir=${MIDAS_COMPILE_ADD_CODECOVERAGE_OPTIONS}
typeset -r toplevelCmd=$"git rev-parse --show-toplevel"
typeset -r toplevel=$(eval ${toplevelCmd})
typeset -r defautSrcDir=${toplevel}/src
typeset -r revstring=$(${toplevel}/midas.version.sh)
typeset -r defaultWebDir=~/public_html/midas/codecoverage-${revstring}

set +u

eval `cclargs_lite                                                  \
 -cov_dir   ${defautCovDir} ${defautCovDir} "[code coverage directory for files generated by the execution (default: value of env variable 'MIDAS_COMPILE_ADD_CODECOVERAGE_OPTIONS']"      \
 -src_dir   ${defautSrcDir} ${defautSrcDir} "[ (default: result of '${toplevelCmd}')]" \
 -web_dir   ${defaultWebDir} ${defaultWebDir} "[directory for storing html pages (default: ${defaultWebDir})]"   \
 -setmx      no yes "[Insert 'set -x' in the script (default: no)]" \
 ++ $*`

## load Intel compiler so gain access to 'codecov' and 'profmerge':
. r.load.dot eccc/mrd/rpn/code-tools/ENV/cdt-1.5.7-inteloneapi-2022.1.2
## Or, we could just source the 'config.dot.sh'

## On met 'set -u' apres le 'cclargs' et le 'r.load.dot' parce que
## sinon le script se plaint qu'il y a des variable non-initialisees.
## Mais a ce point-ci, ca n'est pas trop grave puisqu'on fait
## confiance a 'cclargs' qu'il gere ses trucs comme il faut.
set -u

function say {
    printf "${scriptName}: ${*}\n"
} ## End of 'function say'

function die {
    say $@ >&2
    exit 1
} ## End of 'function die'

typeset -r cov_dir
typeset -r src_dir
typeset -r web_dir
typeset -r setmx

if [ "${setmx}" != yes -a "${setmx}" != no ]; then
    die "Only 'yes' or 'no' is accepted for argument '-setmx'.  You gave '${setmx}'."
fi
if [ "${setmx}" = yes ]; then
    set -x
fi

if [ -z "${cov_dir}" -o "${cov_dir}" = no -o ! -d "${cov_dir}"  ]; then
    die "The argument '-cov_dir' must be set to an existing directory containing the code coverage information.  We have '${cov_dir}'."
fi
#if [ ! -d "${cov_dir}/codecov" ]; then
#    die "The environment variable 'cov_dir' must be set to an existing directory in which the code coverage information is in the subdirectory 'codecov'."
#fi

cd ${cov_dir}

## This tool is documented at https://www.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top/optimization-and-programming/tools/pgo-tools/profmerge-and-proforder-tools.html
say "Calling profmerge"
profmerge -cov_dir . -prof_dir . -src_root ${src_dir} -verbose
#profmerge -cov_dir ./codecov -prof_dir ./codecov -src_root ${src_dir} -verbose

say "Calling codecov"
## This tool is documented at https://www.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top/optimization-and-programming/tools/pgo-tools/code-coverage-tool.html
codecov -ccolor '#d7fad2'

say "Copy the html report to ${web_dir}"
mkdir -p ${web_dir}
rsync -a CodeCoverage CODE_COVERAGE.HTML ${web_dir}

if [[ "${web_dir}" = ${HOME}/public_html/* ]]; then
    typeset -r url=${web_dir/${HOME}\/public_html/http:\/\/goc-dx.science.gc.ca\/~${USER}}
    echo "You can look at the code coverage report at this URL:"
    echo "    ${url}/CODE_COVERAGE.HTML"
else
    echo "You can find the code coverage report in:"
    echo "    ${web_dir}/CODE_COVERAGE.HTML"
fi
