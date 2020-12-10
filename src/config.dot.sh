#! /bin/sh

###########################################################
##
##  USER CONFIGURATION
##
###########################################################
set -x
BACKEND=daley
FRONTEND=eccc-ppp4
JOBNAME=midasCompilation
DIR_BLD_ROOT=${HOME}/data_maestro/ords/midas-bld
DIR_BLD_LINK=../compiledir
NCORES=8
VERBOSE=2
CLEAN=false
DIRECT_FRONTEND_COMPILE=false
FRONTEND_PLAT=ubuntu-18.04-skylake-64

###########################################################
##  compilation and SSM needed for compilation
##
## -- should not change that
DOT_CONFIG=./programs/commons/compile_setup.sh
set +x

##  linking the build directory where it used to be
([ -d  ${DIR_BLD_LINK} ]||[ -L ${DIR_BLD_LINK} ]) \
&& echo "${DIR_BLD_LINK}  already exists: not creating link." \
|| ln -s ${DIR_BLD_ROOT} ${DIR_BLD_LINK}

##  sourcing compilation configuration and SSM packages
source ${DOT_CONFIG} || false


## check makedepf90 install
if ! which makedepf90
then 
    echo "<!> makedepf90 unavailable on the system."
    echo "Check ${DOT_CONFIG} for proper SSM package"
    false 
fi

export BACKEND
export FRONTEND
export JOBNAME
export DIR_BLD_ROOT
export NCORES
export VERBOSE
export CLEAN
export DIRECT_FRONTEND_COMPILE
export FRONTEND_PLAT
