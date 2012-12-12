#!/bin/ksh
nompi=$1
if [ "$nompi" = "NOMPI" -o "$nompi" = "nompi" ] 
then
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "!!Compiling for a NON-MPI executable!!"
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  MPILIB="rpn_commstubs_40007 rpn_comm_40007"
  MPIKEY=""
  ABSTAG="_NOMPI"
else
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "!!Compiling for an MPI executable!!"
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  MPILIB="rpn_comm_40007"
  MPIKEY="-mpi"
  ABSTAG=""
fi

trunkdir=$PWD

# automatically set the global revision number in comct0.cdk by
# replacing the string XXXXX with the actual revision number
revnum=`ssh alef "cd $trunkdir ; svnversion"`
cat comct0_template.cdk |sed "s/XXXXX/${revnum}/g" > comct0.cdk

compiledir=${storage_model:-"../compiledir"}/enkf_pturb
echo "COMPILEDIR IS SET TO: $compiledir"
mkdir -p $compiledir
cd $compiledir
rm -f *.o *.f *.f90 *.mod

compiledir=$PWD

. s.ssmuse.dot Xlf13.108
. s.ssmuse.dot rmnlib-dev
. s.ssmuse.dot devtools
. /ssm/net/hpcs/shortcuts/ssmuse_ssm_v10.sh
. s.ssmuse.dot CMDN/vgrid/3.4.0
. s.ssmuse.dot rpn_comm


VAR3D_VERSION="11.2.1"
LIBAPPL="descrip $MPILIB "
LIBSYS="lapack blas mass"
LIBRMN="rmn_013_rc2"
LIBEXTRA="rtools hpm_r"
MODRTTOV="RTTOV8.7"
MODBURP="BURP1.3"
DEFINE="-DNEC=nec -DIBM=ibm"
ABI="_multi"
COMPF="-openmp $MPIKEY "
FCOMPF="-options_comp"

LIBPATH2="./ $LIBPATH"

echo "LIBPATH2="
echo $LIBPATH2
echo "INCLUDES="
echo $INCLUDES

trunkfiles="abort.ftn bmatrix_mod.ftn90 bmatrixensemble_mod.ftn90 bmatrixhi_mod.ftn90 \
            columndata_mod.ftn90 comfilt.cdk enkf_pturb.ftn controlvector_mod.ftn90 \
            enkf_pturb.ftn gasdev.ftn gaussgrid_mod.ftn90 gdout2.ftn getfldprm2.ftn \
            getstamplist.ftn90 gridstatevector_mod.ftn90 maincompileswitch.inc \
            matsqrt.ftn mpi_mod.ftn90 utils_3dvar.ftn \
            varnamelist_mod.ftn90 varout.ftn verticalcoord_mod.ftn90 comct0.cdk"

moduloptfiles="dsyev.ftn"

cd ${trunkdir}
cp -f ${trunkfiles} ${compiledir}
cp -f ${moduloptfiles} ${compiledir}
cd ${trunkdir}/shared; ls -1F | grep -v '/' | grep -v "*" | cpio -pl ${compiledir}
cd ${compiledir}

rm -f *.ftn~ *.ftn90~

echo "STARTING COMPILATION AT:"
date

echo "compiling low-level independent modules"
SRC0="mathphysconstants_mod.ftn90 earthconstants_mod.ftn90 mpi_mod.ftn90 bufr_mod.ftn90"
s.compile $INCLUDES $COMPF -O -src $SRC0 > listing0 2>&1
grep fail listing0
if [ $? = "0" ] ; then exit ; fi

echo "compiling most of the new modules"
SRC1="controlvector_mod.ftn90 fft_mod.ftn90"
SRC1="$SRC1 gaussgrid_mod.ftn90 globalspectraltransform_mod.ftn90 obsspacedata_mod.ftn90 random_mod.ftn90 varnamelist_mod.ftn90 verticalcoord_mod.ftn90"
SRC1="$SRC1 columndata_mod.ftn90 gridstatevector_mod.ftn90"
SRC1="$SRC1 bmatrixensemble_mod.ftn90 bmatrixhi_mod.ftn90"
SRC1="$SRC1 bmatrix_mod.ftn90"
s.compile $INCLUDES $COMPF -O -src $SRC1 > listing1 2>&1
grep fail listing1
if [ $? = "0" ] ; then exit ; fi

echo "compiling remaining ftn ftn90..."
filelist=""
for i in *.ftn *.ftn90
do
  xx=`echo $i |grep -v _mod.ftn` 
  filelist="$filelist $xx"
done
s.compile $INCLUDES $COMPF -O -src $filelist > listing4 2>&1
grep fail listing4
if [ $? = "0" ] ; then exit ; fi

echo "building the executable..."
s.compile -O -abi $ABI $COMPF $INCLUDES -libpriv -libpath $LIBPATH2 -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o enkf_pturb_p7.abs$ABSTAG > listing5 2>&1

grep -i ERROR listing?
if [ $? = "0" ] ; then exit; echo "ERROR found: STOP" ; fi

rm -f *.ftn* *.cdk* *.h

echo "FINISHED COMPILATION AT:"
date

