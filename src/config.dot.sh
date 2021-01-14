#! /bin/sh

###########################################################
##
##  USER CONFIGURATION
##
###########################################################
set -x
BACKEND=${BACKEND:-daley}
FRONTEND=${FRONTEND:-eccc-ppp4}
JOBNAME=${JOBNAME:-midasCompilation}
COMPILEDIR_MIDAS_MAIN=${COMPILEDIR_MIDAS_MAIN:-${HOME}/data_maestro/ords/midas-bld}
NCORES=${NCORES:-8}
VERBOSE=${VERBOSE:-2}
CLEAN=${CLEAN:-false}
DIRECT_FRONTEND_COMPILE=${DIRECT_FRONTEND_COMPILE:-false}
FRONTEND_PLAT=${FRONTEND_PLAT:-ubuntu-18.04-skylake-64}


###########################################################
##  SSM Packaging configuration 
##
SSM_TARGET=${SSM_TARGET:-/fs/ssm/eccc/mrd/rpn/anl/midas}
SSM_PKGNAME=${SSM_PKGNAME:-midas}
SSM_MAINTAINER=${SSM_MAINTAINER:-ervig.lapalme@canada.ca}
SSM_DESCRIPTION=${SSM_DESCRIPTION:-"The Modular and Integrated Data Assimilation System"}
SSM_GITREPO=${SSM_GITREPO:-https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas.git}
SSM_VERSION=${SSM_VERSION:-3.6.0}


###########################################################
##  LESS-USER-FRIENDLY CONFIGURATION
##
##  these should not be changed unless you know what you're doing
##  it can impact the maestro testing suite or the cleaning targets
##  in unwated ways
__root=$(git rev-parse --show-toplevel)
__exec_leafdir_midas=${__exec_leafdir_midas:-midas_abs}
__install_always_midas=${__install_always_midas:-true}
__compiledir_link=${__compiledir_link:-${__root}/compiledir}


###########################################################
##  compilation and SSM needed for compilation
##
## -- should not change that
__midas_dot_cfg=${PWD}/programs/commons/compile_setup.sh
set +x

##  linking the build directory where it used to be
([ -d  ${__compiledir_link} ]||[ -L ${__compiledir_link} ]) \
&& echo "${__compiledir_link}  already exists: not creating link." \
|| ln -s ${COMPILEDIR_MIDAS_MAIN} ${__compiledir_link}

##  sourcing compilation configuration and SSM packages
source ${__midas_dot_cfg} || false


## check makedepf90 install
if ! which makedepf90
then 
    echo "<!> makedepf90 unavailable on the system."
    echo "Check ${__midas_dot_cfg} for proper SSM package"
    false 
fi

export BACKEND
export FRONTEND
export JOBNAME
export __exec_leafdir_midas
export __install_always_midas 
export NCORES
export VERBOSE
export CLEAN
export DIRECT_FRONTEND_COMPILE
export FRONTEND_PLAT

export SSM_TARGET
export SSM_PKGNAME
export SSM_MAINTAINER
export SSM_DESCRIPTION
export SSM_GITREPO
export SSM_VERSION
