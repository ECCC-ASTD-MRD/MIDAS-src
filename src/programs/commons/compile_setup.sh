#!/bin/bash

# This script contains general configuration and should be dotted in other compile scripts

echo
echo "... compile_setup.sh script"

# User-specified compilation options
#export COMPILE_MIDAS_COMPF_GLOBAL="-DCODEPRECISION_INCR_REAL_SINGLE"
#export COMPILE_MIDAS_COMPF_GLOBAL="-DCODEPRECISION_SPECTRANS_REAL_SINGLE"
if [ -n "${COMPILE_MIDAS_COMPF_GLOBAL}" ];then
     echo "..."
     echo "... Additional user-specified compilation options = ${COMPILE_MIDAS_COMPF_GLOBAL}"
     echo "..."
fi

# Set the optimization level
if [ "${ORDENV_PLAT}" = ubuntu-18.04-skylake-64 -o "${ORDENV_PLAT}" = sles-15-skylake-64-xc50 -o "${ORDENV_PLAT}" = rhel-8-icelake-64 ];then
    FOPTMIZ=4
else
    echo "... This platform 'ORDENV_PLAT=${ORDENV_PLAT}' is not supported."
    exit 1
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
    set -e
    if [ "$(vercomp 1.11.0 ${EC_ATOMIC_PROFILE_VERSION})" = '>' ]; then
        echo "EC_ATOMIC_PROFILE_VERSION=${EC_ATOMIC_PROFILE_VERSION} but should be greater or equal to 1.11.0"
        echo "Please use login profile greater of equal to /fs/ssm/eccc/mrd/ordenv/profile/1.11.0"
        exit 1
    fi
}

#----------------------------------------------------------------
#  Set up dependent librarys and tools. 
#---------------------------------------------------------------
if [ "${ORDENV_PLAT}" = ubuntu-18.04-skylake-64 ]; then
    check_ec_atomic_profile_version
    echo "... loading compiler hpco/exp/intelpsxe-cluster-19.0.3.199"
    . ssmuse-sh -d hpco/exp/intelpsxe-cluster-19.0.3.199
    LIBIRC=irc
    echo "... loading code tools rpn/code-tools/1.5.0"
    . ssmuse-sh -d eccc/mrd/rpn/code-tools/1.5.0
elif [ "${ORDENV_PLAT}" = sles-15-skylake-64-xc50 ]; then
    check_ec_atomic_profile_version
    echo "... loading Intel compiler"
    . r.env.dot --comp 19.0.5
    echo "... loaded compiler ${COMP_ARCH}"
    echo "... loading craype-hugepages16M"
    module load craype-hugepages16M
elif [ "${ORDENV_PLAT}" =  rhel-8-icelake-64 ]; then
    echo "... loading rpn/code-tools/ENV/cdt-1.5.3-b-inteloneapi-2021.4.0"
    . r.load.dot rpn/code-tools/ENV/cdt-1.5.3-b-inteloneapi-2021.4.0
else
    echo "... This platform 'ORDENV_PLAT=${ORDENV_PLAT}' is not supported."
    exit 1
fi

## for hdf5
HDF5_LIBS="hdf5hl_fortran hdf5_hl hdf5_fortran hdf5 z"

if [ "${ORDENV_PLAT}" = ubuntu-18.04-skylake-64 ]; then
    ## for rmn, rpncomm
    echo "... loading eccc/mrd/rpn/libs/19.6.0"
    . r.load.dot eccc/mrd/rpn/libs/19.6.0
    ## for openmpi
    echo "... loading hpco/exp/openmpi/openmpi-3.1.2--hpcx-2.2.0--ofed-4.4.2--intel-2019.0.045"
    . ssmuse-sh -d hpco/exp/openmpi/openmpi-3.1.2--hpcx-2.2.0--ofed-4.4.2--intel-2019.0.045
    echo "... loading hdf5"
    . ssmuse-sh -x hpco/exp/hdf5-netcdf4/serial/static/intel-19.0.3.199/01
elif [ "${ORDENV_PLAT}" = sles-15-skylake-64-xc50 ]; then
    ## for rmn, rpncomm
    echo "... loading eccc/mrd/rpn/libs/19.6.0"
    . r.load.dot eccc/mrd/rpn/libs/19.6.0
    echo "... loading cray-hdf5"
    module load cray-hdf5
    echo "... loading cray-netcdf"
    module load cray-netcdf
    echo "... loading hpco/exp/sqlite/3.29.0"
    . ssmuse-sh -d hpco/exp/sqlite/3.29.0
elif [ "${ORDENV_PLAT}" =  rhel-8-icelake-64 ]; then
    ## for rmn, rpncomm
    echo "... loading eccc/mrd/rpn/libs/19.7.1"
    . r.load.dot eccc/mrd/rpn/libs/19.7.1
    echo "... loading hdf5"
    . ssmuse-sh -d main/opt/hdf5-netcdf4/serial/static/inteloneapi-2021.3.0/01
fi

if [ "${ORDENV_PLAT}" = ubuntu-18.04-skylake-64 -o "${ORDENV_PLAT}" = sles-15-skylake-64-xc50 ]; then
    ## for 'vgrid'
    echo "... loading eccc/mrd/rpn/vgrid/6.5.0"
    . ssmuse-sh -d eccc/mrd/rpn/vgrid/6.5.0
    VGRID_LIBNAME="vgrid"

    echo "... loading eccc/cmd/cmda/libs/19.6.0-2/${COMP_ARCH}"
    . ssmuse-sh -d eccc/cmd/cmda/libs/19.6.0-2/${COMP_ARCH}

    ## For 'perftools' needed for TMG timings
    if [ "${ORDENV_PLAT}" = sles-15-skylake-64-xc50 ]; then
        echo "... loading main/opt/perftools/perftools-2.0/PrgEnv-intel-6.0.5"
        . ssmuse-sh -x main/opt/perftools/perftools-2.0/PrgEnv-intel-6.0.5
    else
        echo "... loading main/opt/perftools/perftools-2.0/${COMP_ARCH}"
        . ssmuse-sh -x main/opt/perftools/perftools-2.0/${COMP_ARCH}
    fi

    echo "... loading eccc/mrd/rpn/anl/rttov/12v1.5.1/${COMP_ARCH}"
    . r.load.dot eccc/mrd/rpn/anl/rttov/12v1.5.1/${COMP_ARCH}

    ## for 'random_tools'
    echo "... loading eccc/mrd/rpn/anl/random_tools/Release_1.0.0-HPCRU1"
    . ssmuse-sh -d eccc/mrd/rpn/anl/random_tools/Release_1.0.0-HPCRU1
elif [ "${ORDENV_PLAT}" =  rhel-8-icelake-64 ]; then
    ## for 'vgrid'
    echo "... loading eccc/mrd/rpn/vgrid/6.6.0"
    . ssmuse-sh -d eccc/mrd/rpn/vgrid/6.6.0
    VGRID_LIBNAME="vgrid"

    echo "... loading eccc/cmd/cmda/libs/19.7.1/${COMP_ARCH}"
    . ssmuse-sh -d eccc/cmd/cmda/libs/19.7.1/${COMP_ARCH}

    echo "... loading main/opt/perftools/perftools-2.0/${COMP_ARCH}"
    . ssmuse-sh -x main/opt/perftools/perftools-2.0/${COMP_ARCH}

    echo "... /home/erv000/SSM/rttov/12v1.6"
    . r.load.dot /home/erv000/SSM/rttov/12v1.6

    ## for 'random_tools'
    echo "... loading /home/erv000/SSM/randomtools/hpcr-u2"
    . ssmuse-sh -d /home/erv000/SSM/randomtools/hpcr-u2
fi

COMPF_GLOBAL="-openmp -mpi ${COMPILE_MIDAS_COMPF_GLOBAL}"
OPTF="-check noarg_temp_created -no-wrap-margin -warn all -warn errors"
if [ "${ORDENV_PLAT}" = ubuntu-18.04-skylake-64 -o "${ORDENV_PLAT}" = rhel-8-icelake-64 ]; then
    OPTF="-mkl ${OPTF}"
elif [ "${ORDENV_PLAT}" = sles-15-skylake-64-xc50 ]; then
    OPTF="${OPTF}"
else
    echo "... This platform 'ORDENV_PLAT=${ORDENV_PLAT}' is not supported."
    exit 1
fi

if [ "${COMPILE_MIDAS_ADD_DEBUG_OPTIONS:-no}" = yes ]; then
    FOPTMIZ=0
    COMPF_NOC="${COMPF_GLOBAL} ${OPTF} -debug"
    COMPF="${COMPF_NOC} -check all -fp-speculation=safe -init=snan,arrays"
    echo "... > !WARNING! You are compiling in DEBUG MODE: '${COMPF}'"
else
    COMPF="${COMPF_GLOBAL} ${OPTF}"
    COMPF_NOC=${COMPF}
fi

GPP_INCLUDE_PATH="$(s.prefix -I $(s.generate_ec_path --include))"
GPP_OPTS="-lang-f90+ -chop_bang -gpp -F ${GPP_INCLUDE_PATH} -D__FILE__=\"#file\" -D__LINE__=#line"
