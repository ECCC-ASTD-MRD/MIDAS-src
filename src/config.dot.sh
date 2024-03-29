#! /bin/sh

__toplevel=$(git rev-parse --show-toplevel)
__revstring=$(${__toplevel}/midas.version.sh)
__revnum=$(echo ${__revstring} | sed -e 's/v_\([^-]*\)-.*/\1/')
__status=true

set -x
###########################################################
##
##  USER CONFIGURATION
##
###########################################################
MIDAS_COMPILE_DIR_MAIN=${MIDAS_COMPILE_DIR_MAIN:-${HOME}/data_maestro/ords/midas-bld}
MIDAS_COMPILE_ADD_DEBUG_OPTIONS=${MIDAS_COMPILE_ADD_DEBUG_OPTIONS:-no}
MIDAS_COMPILE_CODECOVERAGE_DATAPATH=${MIDAS_COMPILE_CODECOVERAGE_DATAPATH:-}
MIDAS_COMPILE_FRONTEND=${MIDAS_COMPILE_FRONTEND:-ppp5}
MIDAS_COMPILE_CLEAN=${MIDAS_COMPILE_CLEAN:-true}
MIDAS_COMPILE_COMPF_GLOBAL=${MIDAS_COMPILE_COMPF_GLOBAL:-}
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
MIDAS_ABS_LEAFDIR=${MIDAS_ABS_LEAFDIR:-midas_abs}
MIDAS_MAKEDEP_TIMEOUT=${MIDAS_MAKEDEP_TIMEOUT:-5s}
__install_always_midas=true
__compiledir_link=${__compiledir_link:-${__toplevel}/compiledir}
__build_dir_version=${MIDAS_COMPILE_DIR_MAIN}/midas_bld-${__revstring}
__keep_jobsubmit_ofile=false
__ordsoumet_wallclock=${__ordsoumet_wallclock:-20}

##  linking the build directory where it used to be
([ -d  ${__compiledir_link} ]||[ -L ${__compiledir_link} ]) \
&& echo "${__compiledir_link}  already exists: not creating link." \
|| ln -s ${MIDAS_COMPILE_DIR_MAIN} ${__compiledir_link}

###########################################################
##  compilation and SSM needed for compilation
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
    FOPTMIZ=4
else
    echo "... This platform 'ORDENV_PLAT=${ORDENV_PLAT}' is not supported."
    return 1
fi

## https://stackoverflow.com/a/4025065
## if $1 = $2, returns '='
## if $1 < $2, returns '<'
## if $1 > $2, returns '>'
vercomp () {
    if [[ $1 == $2 ]]
    then
        echo '='
    fi
    local IFS=.
    local i ver1=($1) ver2=($2)
    # fill empty fields in ver1 with zeros
    for ((i=${#ver1[@]}; i<${#ver2[@]}; i++))
    do
        ver1[i]=0
    done
    for ((i=0; i<${#ver1[@]}; i++))
    do
        if [[ -z ${ver2[i]} ]]
        then
            # fill empty fields in ver2 with zeros
            ver2[i]=0
        fi
        if ((10#${ver1[i]} > 10#${ver2[i]}))
        then
            echo '>'
            return
        fi
        if ((10#${ver1[i]} < 10#${ver2[i]}))
        then
            echo '<'
            return
        fi
    done
}

check_ec_atomic_profile_version () {
    if [ "$(vercomp 1.28.0 ${EC_ATOMIC_PROFILE_VERSION})" = '>' ]; then
	echo
	echo
        echo "+---| ERROR ERROR ERROR |----------------------------------------------------------------"
        echo "  EC_ATOMIC_PROFILE_VERSION=${EC_ATOMIC_PROFILE_VERSION} but should be greater or equal to 1.28.0"
        echo "  Please use login profile greater of equal to /fs/ssm/eccc/mrd/ordenv/profile/1.28.0"
        echo "+---| ERROR ERROR ERROR |----------------------------------------------------------------"
	echo
	echo
        return 1
    fi
}


#----------------------------------------------------------------
#  Set up dependent librarys and tools. 
#---------------------------------------------------------------
echo "... loading eccc/mrd/rpn/code-tools/ENV/cdt-1.6.2/SCIENCE/inteloneapi-2022.1.2"
. r.load.dot eccc/mrd/rpn/code-tools/ENV/cdt-1.6.2/SCIENCE/inteloneapi-2022.1.2

## for hdf5
HDF5_LIBS="hdf5hl_fortran hdf5_hl hdf5_fortran hdf5 z"

## for rmn, vgrid, rpncomm
VGRID_LIBNAME="vgrid"
echo "... loading eccc/mrd/rpn/libs/20230531"
. r.load.dot eccc/mrd/rpn/libs/20230531
echo "... loading hdf5"
. ssmuse-sh -d main/opt/hdf5-netcdf4/serial/static/${COMP_ARCH}/01

echo "... loading eccc/cmd/cmda/libs/20230531/${COMP_ARCH}"
. ssmuse-sh -d eccc/cmd/cmda/libs/20230531/${COMP_ARCH}

echo "... loading main/opt/perftools/perftools-2.0/${COMP_ARCH}"
. ssmuse-sh -x main/opt/perftools/perftools-2.0/${COMP_ARCH}

echo "... loading eccc/mrd/rpn/anl/rttov/13v1.3/${COMP_ARCH}"
. r.load.dot eccc/mrd/rpn/anl/rttov/13v1.3/${COMP_ARCH}

## loading makedep90
echo "... loading makedepf90"
. ssmuse-sh -d eccc/mrd/rpn/anl/makedepf90/2.8.9

COMPF_GLOBAL="-openmp -mpi ${MIDAS_COMPILE_COMPF_GLOBAL}"
OPTF="-check noarg_temp_created -no-wrap-margin -warn all -warn errors"
OPTF="-qmkl ${OPTF} -warn noexternal"

if [ "${MIDAS_COMPILE_ADD_DEBUG_OPTIONS:-no}" = yes ]; then
    FOPTMIZ=0
    COMPF_NOC="${COMPF_GLOBAL} ${OPTF} -debug"
    COMPF="${COMPF_NOC} -check all -fp-speculation=safe -init=snan,arrays"
    echo "... > !WARNING! You are compiling in DEBUG MODE: '${COMPF}'"
else
    COMPF="${COMPF_GLOBAL} ${OPTF}"
    COMPF_NOC=${COMPF}
fi

if [ -n "${MIDAS_COMPILE_CODECOVERAGE_DATAPATH}" ]; then
    echo "... > !WARNING! You are compiling in CODE COVERAGE MODE: '${COMPF}'"
    FOPTMIZ=0
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

GPP_INCLUDE_PATH="$(s.prefix -I $(s.generate_ec_path --include))"
GPP_OPTS="-lang-f90+ -chop_bang -gpp -F ${GPP_INCLUDE_PATH} -D__FILE__=\"#file\" -D__LINE__=#line"

## check makedepf90 install
if ! which makedepf90
then 
    echo "<!> makedepf90 unavailable on the system."
    __status=false
fi

## loading docopt for analyzeDep.py
## https://gitlab.science.gc.ca/hpc/hpcr_upgrade_2/issues/252
. ssmuse-sh -x comm/eccc/arqi/modules-python/1.0


export COMPF
export FOPTMIZ
export GPP_OPTS

export LIBIRC
export HDF5_LIBS
export VGRID_LIBNAME

export MIDAS_COMPILE_FRONTEND
export MIDAS_COMPILE_JOBNAME
export MIDAS_ABS_LEAFDIR
export MIDAS_COMPILE_VERBOSE
export MIDAS_MAKEDEP_TIMEOUT

export MIDAS_SSM_TARGET
export MIDAS_SSM_PKGNAME
export MIDAS_SSM_MAINTAINER
export MIDAS_SSM_DESCRIPTION
export MIDAS_SSM_GITREPO
export MIDAS_SSM_VERSION

# config return status
${__status}

# check version of profile
check_ec_atomic_profile_version
