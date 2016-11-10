#!/bin/ksh

if [ "${BASE_ARCH}" != "AIX-powerpc7" ] ; then
    echo "Error: This code must be compiled on the IBM machine..."
    exit 1
fi

echo " Compiling for an MPI executable (even if calcstats is not MPI capable)"
echo ""
MPILIBDIR="" #"-libpath /users/dor/arma/anl/userlibs/${BASE_ARCH}/xlf13" # JFC : Mesure temporaire pour avoir acces a
MPILIB="" # "_anl_rpn_comm_4051103"                                       #       la S-R RPN_COMM_adj_halo8 de M. Valin
MPIKEY="" #"-mpi"
ABSTAG=""

calcstatsdir=$PWD
trunkdir=$PWD/../

## for s.compile
echo "loading hpcs/201402/02/base"
. ssmuse-sh -d hpcs/201402/02/base

## for the compiler
echo "loading compiler hpcs/ext/xlf_13.1.0.10"
. ssmuse-sh -d hpcs/ext/xlf_13.1.0.10

## for rmn, lapack_3.4.0, rpncomm
echo "loading rpn/libs/15.2"
. ssmuse-sh -d rpn/libs/15.2

## For hpcsperf needed for TMG timings
. ssmuse-sh -d hpcs/exp/aspgjdm/perftools

LIBAPPL="$MPILIB"

LIBSYS="essl mass"
LIBRMN="rmn"
LIBEXTRA="hpcsperf lapack-3.4.0"
COMPF_NOC="-openmp $MPIKEY -O"
COMPF="$COMPF_NOC"
#COMPF="$COMPF_NOC -debug DEBUG -optf=-C "

# Create and Move to compilation directory
compiledir="../../compiledir_optimization"
echo "COMPILEDIR IS SET TO: $compiledir"
mkdir -p $compiledir
cd $compiledir

compiledir=$PWD

cd ${trunkdir}
cp -f optimization_tests/*.ftn90 ${compiledir}
cd ${compiledir}

echo "STARTING COMPILATION AT:" 
date
s.compile $COMPF -src addEnsMember.ftn90 ${MPILIBDIR} -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -o addEnsMember_p7.abs
rm -f *.ftn* 
echo "FINISHED COMPILATION AT:"
date
