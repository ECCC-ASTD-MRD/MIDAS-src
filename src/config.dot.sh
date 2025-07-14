#! /bin/sh
##
## DO NOT MODIFY THIS FILE
##
## as it is part of the versioned repository and contains the default configuration
## for the build environment.
##
## If you need to change the default configuration provided by this file, either
## export the desired variables in your environment before sourcing it or define
## them in your profile.
##

__toplevel=$(git rev-parse --show-toplevel)
__revstring=$(${__toplevel}/midas.version)
__revnum=$(echo ${__revstring} | sed -e 's/v_\([^-]*\)-.*/\1/')
__status=true

set -x
###########################################################
##
##  USER CONFIGURATION
##
###########################################################

## The variable 'MIDAS_COMPILE_DIR_MAIN' can be defined as
## 'build_directory_local_to_the_repository' and so the build
## directory will be in '${__toplevel}/compiledir' (see variable '${__compiledir_link}' below).
## If it is not defined, we will use a default which appends the
## basename of the toplevel directory to '${HOME}/data_maestro/ords/midas-bld'.

MIDAS_COMPILE_ADD_DEBUG_OPTIONS=${MIDAS_COMPILE_ADD_DEBUG_OPTIONS:-no}
MIDAS_COMPILE_APPEND_VERSION_ID_BUILDDIR=${MIDAS_COMPILE_APPEND_VERSION_ID_BUILDDIR:-true}
MIDAS_COMPILE_CODECOVERAGE_DATAPATH=${MIDAS_COMPILE_CODECOVERAGE_DATAPATH:-}
MIDAS_COMPILE_FRONTEND=${MIDAS_COMPILE_FRONTEND:-ppp5}
MIDAS_COMPILE_CLEAN=${MIDAS_COMPILE_CLEAN:-true}
MIDAS_COMPILE_COMPF_GLOBAL=${MIDAS_COMPILE_COMPF_GLOBAL:-}
MIDAS_COMPILE_HEADNODE_FRONTEND=${MIDAS_COMPILE_HEADNODE_FRONTEND:-false}
MIDAS_COMPILE_JOBNAME=${MIDAS_COMPILE_JOBNAME:-midasCompilation}
MIDAS_COMPILE_KEEP_LISTING=${MIDAS_COMPILE_KEEP_LISTING:-false}
MIDAS_COMPILE_NCORES=${MIDAS_COMPILE_NCORES:-8}
MIDAS_COMPILE_VERBOSE=${MIDAS_COMPILE_VERBOSE:-1}
MIDAS_COMPILE_OPTIMIZE_REPORT=${MIDAS_COMPILE_OPTIMIZE_REPORT:-no}

set +x
__compiledir_link=${__compiledir_link:-${__toplevel}/compiledir}

## If 'MIDAS_COMPILE_DIR_MAIN' is equal to the special value
## 'build_directory_local_to_the_repository', then we keep the build
## directory local to the Git repository.
if [ "${MIDAS_COMPILE_DIR_MAIN}" = build_directory_local_to_the_repository ]; then
    echo "Creating '${__compiledir_link}' since 'MIDAS_COMPILE_DIR_MAIN' is defined but empty"
    [ ! -d "${__compiledir_link}" ] && mkdir ${__compiledir_link}
    __midas_compile_dir_main=${__compiledir_link}
else
    if [ -n "${MIDAS_COMPILE_DIR_MAIN}" ]; then
        if [[ "${MIDAS_COMPILE_DIR_MAIN}" != /* ]]; then
            echo "Please provide of value of MIDAS_COMPILE_DIR_MAIN which is an absolute path"
            echo "   MIDAS_COMPILE_DIR_MAIN=${MIDAS_COMPILE_DIR_MAIN}"
            echo "was given"
            return 1
        fi
        __midas_compile_dir_main=${MIDAS_COMPILE_DIR_MAIN}
    else
        ## If the variable 'MIDAS_COMPILE_DIR_MAIN' is not defined,
        ## we add the leaf part of the toplevel directory to
        ## '${HOME}/data_maestro/ords/midas-bld'.
        __toplevel_leaf=$(basename ${__toplevel})
        __midas_compile_dir_main=${HOME}/data_maestro/ords/midas-bld/${__toplevel_leaf}
    fi

    if [ ! -d "${__midas_compile_dir_main}" ]; then
        mkdir -p ${__midas_compile_dir_main}
    fi
    ##  linking the build directory where it used to be
    if [ -d "${__compiledir_link}" -o -L "${__compiledir_link}" ]; then
        echo "${__compiledir_link} already exists: not creating link."
    else
        ln -s ${__midas_compile_dir_main} ${__compiledir_link}
    fi
fi
set -x
export MIDAS_SOURCE_DIR=${__toplevel}
export MIDAS_COMPILE_DIR_MAIN=${__midas_compile_dir_main}

###########################################################
##  LESS-USER-FRIENDLY CONFIGURATION
##
##  these should not be changed unless you know what you're doing
##  it can impact the maestro testing suite or the cleaning targets
##  in unwated ways
MIDAS_ABS_LEAFDIR=${MIDAS_ABS_LEAFDIR:-midas_abs}
__install_always_midas=true

if [ "${MIDAS_COMPILE_APPEND_VERSION_ID_BUILDDIR}" = true ]; then
    export MIDAS_COMPILE_DIR_BUILD=${MIDAS_COMPILE_DIR_MAIN}/midas_bld-${__revstring}
else
    export MIDAS_COMPILE_DIR_BUILD=${MIDAS_COMPILE_DIR_MAIN}/midas_bld
fi
if [ ! -d "${MIDAS_COMPILE_DIR_BUILD}" ]; then
  mkdir ${MIDAS_COMPILE_DIR_BUILD}
fi
__keep_jobsubmit_ofile=false
__ordsoumet_wallclock=${__ordsoumet_wallclock:-20}

###########################################################
##  compilation needed for compilation
##
## -- should not change that
set +x

# User-specified compilation options
#export MIDAS_COMPILE_COMPF_GLOBAL="-DCODEPRECISION_INCR_REAL_SINGLE"
#export MIDAS_COMPILE_COMPF_GLOBAL="-DCODEPRECISION_SPECTRANS_REAL_SINGLE"
if [ -n "${MIDAS_COMPILE_COMPF_GLOBAL}" ];then
     echo "..."
     echo "... Additional user-specified compilation options = ${MIDAS_COMPILE_COMPF_GLOBAL}"
     echo "..."
fi

# Set the optimization level
if [ "${ORDENV_PLAT}" = rhel-8-icelake-64 ];then
    FOPTMIZ="-O3 -fast-transcendentals -no-prec-div -fpic -ip -no-prec-sqrt"
else
    echo "... This platform 'ORDENV_PLAT=${ORDENV_PLAT}' is not supported."
    return 1
fi

#----------------------------------------------------------------
#  Set up dependent librarys and tools.
#---------------------------------------------------------------
echo "... loading rpn/code-tools/20250521/env/inteloneapi-2022.1.2"
. r.load.dot rpn/code-tools/20250521/env/inteloneapi-2022.1.2

echo "... loading eccc/mrd/rpn/libs/20250604-beta"
. r.load.dot eccc/mrd/rpn/libs/20250604-beta
echo "... loading hdf5"
. ssmuse-sh -d main/opt/hdf5-netcdf4/serial/static/${COMP_ARCH}/01

echo "... loading eccc/mrd/rpn/utils/20250604-beta/burp-tools_20.0.8.2-${COMP_ARCH}_${ORDENV_PLAT}"
. r.load.dot eccc/mrd/rpn/utils/20250604-beta/burp-tools_20.0.8.2-${COMP_ARCH}_${ORDENV_PLAT}

echo "... loading eccc/cmd/cmda/libs/20250604-beta/${COMP_ARCH}"
. ssmuse-sh -x eccc/cmd/cmda/libs/20250604-beta/${COMP_ARCH}

echo "... loading main/opt/perftools/perftools-2.0/${COMP_ARCH}"
. ssmuse-sh -x main/opt/perftools/perftools-2.0/${COMP_ARCH}

echo "... loading eccc/mrd/rpn/anl/rttov/13v2.0/${COMP_ARCH}"
. r.load.dot eccc/mrd/rpn/anl/rttov/13v2.0/${COMP_ARCH}

COMPF_GLOBAL="${MIDAS_COMPILE_COMPF_GLOBAL}"
OPTF="-stand f08 -diag-disable=5268 -check noarg_temp_created -no-wrap-margin -warn all -warn errors"
OPTF="-qmkl ${OPTF} -warn noexternal"

# add compiler option to produce reports on code optimization and deactivate cleaning
if [ "${MIDAS_COMPILE_OPTIMIZE_REPORT:-no}" = yes ]; then
    OPTF="${OPTF} -qopt-report=5"
    echo "... > !WARNING! Compiler optimization reports will be produced in the compile directory."
    echo "... >           To be able to see them, we ensure cleaning is not activated."
    MIDAS_COMPILE_CLEAN=false
fi

if [ "${MIDAS_COMPILE_ADD_DEBUG_OPTIONS:-no}" = yes ]; then
    FOPTMIZ=-O0
    COMPF_NOC="${COMPF_GLOBAL} ${OPTF} -g -ftrapuv"
    COMPF="${COMPF_NOC} -check all -fp-speculation=safe -init=snan,arrays"
    echo "... > !WARNING! You are compiling in DEBUG MODE: '${COMPF}'"
else
    COMPF="${COMPF_GLOBAL} ${OPTF}"
    COMPF_NOC=${COMPF}
fi

if [ -n "${MIDAS_COMPILE_CODECOVERAGE_DATAPATH}" ]; then
    echo "... > !WARNING! You are compiling in CODE COVERAGE MODE: '${COMPF}'"
    FOPTMIZ="-O0"
    [[ "${MIDAS_COMPILE_CODECOVERAGE_DATAPATH}" != /* ]] && {
        echo "Please provide an absolute path to variable 'MIDAS_COMPILE_CODECOVERAGE_DATAPATH'"
        echo "This value was given: ${MIDAS_COMPILE_CODECOVERAGE_DATAPATH}"
        __status=false
    }
    [ ! -d "${MIDAS_COMPILE_CODECOVERAGE_DATAPATH}" ] && mkdir -p ${MIDAS_COMPILE_CODECOVERAGE_DATAPATH}
    [ ! -d "${MIDAS_COMPILE_CODECOVERAGE_DATAPATH}" ] && {
        echo "Could not create the directory ${MIDAS_COMPILE_CODECOVERAGE_DATAPATH}"
        __status=false
    }
    COMPF="${COMPF} -prof-gen=srcpos -prof-dir=${MIDAS_COMPILE_CODECOVERAGE_DATAPATH}"
    COMPF_NOC="${COMPF_NOC} -prof-gen=srcpos -prof-dir=${MIDAS_COMPILE_CODECOVERAGE_DATAPATH}"
fi

## loading docopt for analyzeDep.py
## https://gitlab.science.gc.ca/hpc/hpcr_upgrade_2/issues/252
. ssmuse-sh -x comm/eccc/arqi/modules-python/1.0

# Shortcut commands for cmake and make
alias cado=${PWD}/cado

export COMPF
export FOPTMIZ

export MIDAS_COMPILE_FRONTEND
export MIDAS_COMPILE_JOBNAME
export MIDAS_ABS_LEAFDIR
export MIDAS_COMPILE_VERBOSE

# config return status
${__status}
