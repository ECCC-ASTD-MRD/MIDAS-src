#! /bin/sh

__toplevel=$(git rev-parse --show-toplevel)
__revstring=$(${__toplevel}/midas.version.sh)
__revnum=$(echo ${__revstring} | sed -e 's/v_\([^-]*\)-.*/\1/')

set -x
###########################################################
##
##  USER CONFIGURATION
##
###########################################################
MIDAS_COMPILE_ADD_DEBUG_OPTIONS=${MIDAS_COMPILE_ADD_DEBUG_OPTIONS:-no}
MIDAS_COMPILE_BACKEND=${MIDAS_COMPILE_BACKEND:-daley}
MIDAS_COMPILE_DO_BACKEND=${MIDAS_COMPILE_DO_BACKEND:-true}
MIDAS_COMPILE_BACKEND_ARCH=${MIDAS_COMPILE_BACKEND_ARCH:-sles-15-skylake-64-xc50}
MIDAS_COMPILE_CLEAN=${MIDAS_COMPILE_CLEAN:-true}
MIDAS_COMPILE_COMPF_GLOBAL=${MIDAS_COMPILE_COMPF_GLOBAL:-}
MIDAS_COMPILE_DIR_MAIN=${MIDAS_COMPILE_DIR_MAIN:-${HOME}/data_maestro/ords/midas-bld}
MIDAS_COMPILE_FRONTEND=${MIDAS_COMPILE_FRONTEND:-eccc-ppp4}
MIDAS_COMPILE_HEADNODE_FRONTEND=${MIDAS_COMPILE_HEADNODE_FRONTEND:-false}
MIDAS_COMPILE_JOBNAME=${MIDAS_COMPILE_JOBNAME:-midasCompilation}
MIDAS_COMPILE_KEEP_LISTING=${MIDAS_COMPILE_KEEP_LISTING:-false}
MIDAS_COMPILE_NCORES=${MIDAS_COMPILE_NCORES:-8}
MIDAS_COMPILE_VERBOSE=${MIDAS_COMPILE_VERBOSE:-2}

###########################################################
##  SSM Packaging configuration 
##
MIDAS_SSM_DESCRIPTION=${MIDAS_SSM_DESCRIPTION:-"The Modular and Integrated Data Assimilation System"}
MIDAS_SSM_GITREPO=${MIDAS_SSM_GITREPO:-https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas.git}
MIDAS_SSM_MAINTAINER=${MIDAS_SSM_MAINTAINER:-ervig.lapalme@canada.ca}
MIDAS_SSM_PKGNAME=${MIDAS_SSM_PKGNAME:-midas}
MIDAS_SSM_TARGET=${MIDAS_SSM_TARGET:-/fs/ssm/eccc/mrd/rpn/anl/midas}
MIDAS_SSM_VERSION=${MIDAS_SSM_VERSION:-${__revnum}}


###########################################################
##  LESS-USER-FRIENDLY CONFIGURATION
##
##  these should not be changed unless you know what you're doing
##  it can impact the maestro testing suite or the cleaning targets
##  in unwated ways

__exec_leafdir_midas=midas_abs
__install_always_midas=true
__compiledir_link=${__compiledir_link:-${__toplevel}/compiledir}
__build_dir_version=${MIDAS_COMPILE_DIR_MAIN}/midas_bld-${__revstring}
__keep_jobsubmit_ofile=false
if [ "${MIDAS_COMPILE_DO_BACKEND}" = true ]; then 
  __plat_super=${MIDAS_COMPILE_BACKEND_ARCH}
else
  __plat_super=''
fi

###########################################################
##  compilation and SSM needed for compilation
##
## -- should not change that
__midas_dot_cfg=${__toplevel}/src/programs/commons/compile_setup.sh
set +x

##  linking the build directory where it used to be
([ -d  ${__compiledir_link} ]||[ -L ${__compiledir_link} ]) \
&& echo "${__compiledir_link}  already exists: not creating link." \
|| ln -s ${MIDAS_COMPILE_DIR_MAIN} ${__compiledir_link}

##  sourcing compilation configuration and SSM packages
source ${__midas_dot_cfg} || false


## check makedepf90 install
if ! which makedepf90
then 
    echo "<!> makedepf90 unavailable on the system."
    echo "Check ${__midas_dot_cfg} for proper SSM package"
    false 
fi

export MIDAS_COMPILE_BACKEND
export MIDAS_COMPILE_DO_BACKEND
export MIDAS_COMPILE_FRONTEND
export MIDAS_COMPILE_JOBNAME
export __exec_leafdir_midas
export __install_always_midas 
export MIDAS_COMPILE_NCORES
export MIDAS_COMPILE_VERBOSE
export MIDAS_COMPILE_CLEAN
export MIDAS_COMPILE_HEADNODE_FRONTEND

export MIDAS_SSM_TARGET
export MIDAS_SSM_PKGNAME
export MIDAS_SSM_MAINTAINER
export MIDAS_SSM_DESCRIPTION
export MIDAS_SSM_GITREPO
export MIDAS_SSM_VERSION
