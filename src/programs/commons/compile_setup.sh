#!/bin/bash

# This script contains general configuration and should be dotted in other compile scripts

echo
echo "... compile_setup.sh script"

if [ "${ORDENV_PLAT}" = sles-11-haswell-64-xc40 ];then
    echo "... Switching ORDENV_PLAT from '${ORDENV_PLAT}' to 'sles-11-broadwell-64-xc40'"
    . r.env.dot --arch sles-11-broadwell-64-xc40
    echo "... ORDENV_PLAT=${ORDENV_PLAT}"
fi

# Set the optimization level
if [ "${ORDENV_PLAT}" = ubuntu-14.04-amd64-64 ];then
    FOPTMIZ=2
elif [ "${ORDENV_PLAT}" = ubuntu-18.04-skylake-64 -o "${ORDENV_PLAT}" = sles-15-skylake-64-xc50 ];then
    FOPTMIZ=4
elif [ "${ORDENV_PLAT}" = sles-11-amd64-64 ];then
    FOPTMIZ=4
elif [ "${ORDENV_PLAT}" = sles-11-broadwell-64-xc40 ];then
    FOPTMIZ=4
else
    echo "... This platform 'ORDENV_PLAT=${ORDENV_PLAT}' is not supported.  Only 'ubuntu-14.04-amd64-64' and 'ubuntu-14.04-amd64-64' are."
    exit 1
fi

#----------------------------------------------------------------
#  Set up dependent librarys and tools. 
#---------------------------------------------------------------
if [ "${ORDENV_PLAT}" = ubuntu-14.04-amd64-64 ]; then
    ## for s.compile, s.f90
    echo "... loading hpco/tmp/eccc/201402/06/base"
    . ssmuse-sh -d hpco/tmp/eccc/201402/06/base
    echo "... loading compiler main/opt/intelcomp/intelcomp-2016.1.156"
    . ssmuse-sh -d main/opt/intelcomp/intelcomp-2016.1.156
elif [ "${ORDENV_PLAT}" = ubuntu-18.04-skylake-64 ]; then
    echo "... loading eccc/mrd/rpn/code-tools/01.3"
    . ssmuse-sh -d eccc/mrd/rpn/code-tools/01.3
    echo "... loading compiler hpco/exp/intelpsxe-cluster-19.0.3.199"
    . ssmuse-sh -d hpco/exp/intelpsxe-cluster-19.0.3.199
elif [ "${ORDENV_PLAT}" = sles-11-amd64-64 -o "${ORDENV_PLAT}" = sles-11-broadwell-64-xc40 ]; then
    ## for s.compile, s.f90
    echo "... loading hpco/tmp/eccc/201402/06/base"
    . ssmuse-sh -d hpco/tmp/eccc/201402/06/base
    echo "... loading compiler PrgEnv-intel-5.2.82"
    module load PrgEnv-intel/5.2.82
elif [ "${ORDENV_PLAT}" = sles-15-skylake-64-xc50 ]; then
    echo "... loading eccc/mrd/rpn/code-tools/01.3"
    . ssmuse-sh -d eccc/mrd/rpn/code-tools/01.3
    echo "... loading compiler PrgEnv-intel-6.0.5"
    module load PrgEnv-intel/6.0.5
else
    echo "... This platform 'ORDENV_PLAT=${ORDENV_PLAT}' is not supported.  Only 'ubuntu-14.04-amd64-64' and 'sles-11-amd64-64' are."
    exit 1
fi

## for hdf5
HDF5_LIBS="hdf5hl_fortran hdf5_hl hdf5_fortran hdf5 z"

if [ "${ORDENV_PLAT}" = ubuntu-14.04-amd64-64 ];then
    ## for rmn, rpncomm
    echo "... loading eccc/mrd/rpn/libs/16.2"
    . ssmuse-sh -d eccc/mrd/rpn/libs/16.2
    ## for openmpi
    echo "... loading main/opt/openmpi/openmpi-1.6.5/intelcomp-2016.1.156"
    . ssmuse-sh -d main/opt/openmpi/openmpi-1.6.5/intelcomp-2016.1.156
    echo "... loading hdf5"
    . ssmuse-sh -x comm/eccc/cmdd/rttov/rttov-1.0/serial/disable-shared/intelcomp-2016.1.156
elif [ "${ORDENV_PLAT}" = ubuntu-18.04-skylake-64 ]; then
    ## for rmn, rpncomm
    echo "... loading eccc/mrd/rpn/libs/19.5"
    . r.load.dot eccc/mrd/rpn/libs/19.5
    ## for openmpi
    echo "... loading hpco/exp/openmpi/openmpi-3.1.2--hpcx-2.2.0--ofed-4.4.2--intel-2019.0.045"
    . ssmuse-sh -d hpco/exp/openmpi/openmpi-3.1.2--hpcx-2.2.0--ofed-4.4.2--intel-2019.0.045
    echo "... loading hdf5"
    . ssmuse-sh -x hpco/exp/hdf5-netcdf4/serial/static/intel-19.0.3.199/01
elif [ "${ORDENV_PLAT}" = sles-11-amd64-64 -o "${ORDENV_PLAT}" = sles-11-broadwell-64-xc40 ];then
    ## for rmn, rpncomm
    echo "... loading eccc/mrd/rpn/libs/16.2"
    . ssmuse-sh -d eccc/mrd/rpn/libs/16.2
    echo "... loading cray-hdf5"
    module load cray-hdf5
elif [ "${ORDENV_PLAT}" = sles-15-skylake-64-xc50 ]; then
    ## for rmn, rpncomm
    echo "... loading eccc/mrd/rpn/libs/19.5"
    . r.load.dot eccc/mrd/rpn/libs/19.5
    echo "... loading cray-hdf5"
    module load cray-hdf5
    echo "... loading cray-netcdf"
    module load cray-netcdf
    echo "... loading hpco/exp/sqlite/3.29.0"
    . ssmuse-sh -d hpco/exp/sqlite/3.29.0
fi

if [ "${ORDENV_PLAT}" = ubuntu-14.04-amd64-64 -o "${ORDENV_PLAT}" = sles-11-amd64-64 -o "${ORDENV_PLAT}" = sles-11-broadwell-64-xc40 ]; then
    ## for 'vgrid'
    echo "... loading eccc/cmd/cmdn/vgrid/5.6.9/${COMP_ARCH}"
    . ssmuse-sh -d eccc/cmd/cmdn/vgrid/5.6.9/${COMP_ARCH}

    ## for 'burplib'
    echo "... loading eccc/cmd/cmda/libs/16.2-6/${COMP_ARCH}"
    . ssmuse-sh -d eccc/cmd/cmda/libs/16.2-6/${COMP_ARCH}

    ## For hpcoperf needed for TMG timings
    echo "... loading main/opt/perftools/perftools-2.0/${COMP_ARCH}"
    . ssmuse-sh -d main/opt/perftools/perftools-2.0/${COMP_ARCH}

    echo "... loading eccc/mrd/rpn/anl/rttov/12v1.2"
    . ssmuse-sh -d eccc/mrd/rpn/anl/rttov/12v1.2/${COMP_ARCH}
elif [ "${ORDENV_PLAT}" = ubuntu-18.04-skylake-64 -o "${ORDENV_PLAT}" = sles-15-skylake-64-xc50 ]; then
    ## for 'vgrid'
    echo "... loading eccc/mrd/rpn/vgrid/6.4.5"
    . ssmuse-sh -d eccc/mrd/rpn/vgrid/6.4.5

    echo "... loading eccc/cmd/cmda/libs/19.5/${COMP_ARCH}"
    . ssmuse-sh -d eccc/cmd/cmda/libs/19.5/${COMP_ARCH}

    ## For 'perftools' needed for TMG timings
    echo "... loading main/opt/perftools/perftools-2.0/${COMP_ARCH}"
    . ssmuse-sh -x main/opt/perftools/perftools-2.0/${COMP_ARCH}

    echo "... loading eccc/mrd/rpn/anl/rttov/12v1.3"
    . r.load.dot eccc/mrd/rpn/anl/rttov/12v1.3/${COMP_ARCH}
fi

COMPF_GLOBAL="-openmp -mpi"
OPTF="=-check =noarg_temp_created =-no-wrap-margin"
if [ "${ORDENV_PLAT}" = ubuntu-14.04-amd64-64 -o "${ORDENV_PLAT}" = ubuntu-18.04-skylake-64 ];then
    OPTF="=-mkl ${OPTF}"
elif [ "${ORDENV_PLAT}" = sles-11-amd64-64 -o "${ORDENV_PLAT}" = sles-11-broadwell-64-xc40 -o "${ORDENV_PLAT}" = sles-15-skylake-64-xc50 ];then
    OPTF="${OPTF}"
else
    echo "... This platform 'ORDENV_PLAT=${ORDENV_PLAT}' is not supported.  Only 'ubuntu-14.04-amd64-64' and 'sles-11-amd64-64' are."
    exit 1
fi

if [ "${COMPILE_MIDAS_ADD_DEBUG_OPTIONS:-no}" = yes ]; then
    echo "... > !WARNING! You are compiling in DEBUG MODE: '-debug -C -O 0'"
    COMPF_NOC="${COMPF_GLOBAL} -debug DEBUG -optf ${OPTF}"
    COMPF="${COMPF_NOC} =-C"
    FOPTMIZ=0
else
    COMPF="${COMPF_GLOBAL} -optf ${OPTF}"
    COMPF_NOC=${COMPF}
fi
