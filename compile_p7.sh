#!/bin/ksh
if [ "$1" = "NOMPI" -o "$1" = "nompi" ] 
then
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "!!Compiling for a NON-MPI executable!!"
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  MPILIB="rpn_commstubs301 rpn_comm301"
  MPIKEY=""
  ABSTAG="_NOMPI"
else
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "!!Compiling for an MPI executable!!"
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  MPILIB="rpn_comm301"
  MPIKEY="-mpi"
  ABSTAG=""
fi

trunkdir=$PWD

cd ../
mkdir -p compiledir
cd compiledir

compiledir=$PWD

. s.ssmuse.dot Xlf13.108
. s.ssmuse.dot rmnlib-dev
. s.ssmuse.dot devtools
. s.ssmuse.dot ENV/vgrid/1.0.8

VAR3D_VERSION="11.2.1"
LIBAPPL="modulopt rttov burp_module descrip $MPILIB "
LIBSYS="lapack blas mass"
LIBRMN="rmn_013_rc2"
LIBEXTRA="rtools hpm_r"
MODRTTOV="RTTOV8.7"
MODBURP="BURP1.3"
DEFINE="-DNEC=nec -DIBM=ibm"
ABI="_multi"
COMPF="-openmp $MPIKEY "
FCOMPF="-options_comp"

BASE_INCLUDE="${ARMNLIB}/modeles/ANAL/v_${VAR3D_VERSION}/include/AIX-powerpc7"
INCLUDES="-includes ${BASE_INCLUDE}/${MODBURP} ${BASE_INCLUDE}/${MODRTTOV}"
INCLUDES="$INCLUDES /home/ordenv/ssm-domains9/release/vgriddescriptors/2.0.1/multi/include/AIX-powerpc7/xlf13"

LIBPATH2="./ $LIBPATH"
LIBPATH2="${ARMNLIB}/lib/AIX/xlf13 $LIBPATH2"
LIBPATH2="${ARMNLIB}/modeles/ANAL/v_${VAR3D_VERSION}/lib/AIX-powerpc7 $LIBPATH2"
LIBPATH2="/home/ordenv/ssm-domains9/release/vgriddescriptors/2.0.1/multi/lib/AIX-powerpc7/xlf13 $LIBPATH2"
LIBPATH2="/home/ordenv/ssm-domains1/ssm-rmnlib-dev/multi/lib/AIX-powerpc7/xlf13 $LIBPATH2"

echo "LIBPATH2="
echo $LIBPATH2
echo "INCLUDES="
echo $INCLUDES

cd ${trunkdir};           ls -1F | grep -v '/' | grep -v "*" | grep -v "@" | cpio -pl $compiledir ; cd $compiledir
cd ${trunkdir}/bgcheck;           ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir
cd ${trunkdir}/shared; ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir

rm -f *.ftn~ *.ftn90~

echo "STARTING COMPILATION AT:"
date

echo "compiling low-level independent modules"
SRC0="mathphysconstants_mod.ftn90 earthconstants_mod.ftn90 mpi_mod.ftn90 bufr_mod.ftn90"
s.compile $INCLUDES $COMPF -O -src $SRC0 > listing0 2>&1
grep fail listing0
if [ $? = "0" ] ; then exit ; fi

echo "compiling most of the new modules"
SRC1="controlvector_mod.ftn90 emissivities_mod.ftn90 fft_mod.ftn90 gaussgrid_mod.ftn90 globalspectraltransform_mod.ftn90 obsspacedata_mod.ftn90 random_mod.ftn90 varnamelist_mod.ftn90 verticalcoord_mod.ftn90"
SRC1="$SRC1 columndata_mod.ftn90 gridstatevector_mod.ftn90"
SRC1="$SRC1 bmatrixensemble_mod.ftn90 bmatrixhi_mod.ftn90"
SRC1="$SRC1 bmatrix_mod.ftn90 minimization_mod.ftn90"
s.compile $INCLUDES $COMPF -O -src $SRC1 > listing1 2>&1
grep fail listing1
if [ $? = "0" ] ; then exit ; fi

echo "compiling the old modules (cdk90)..."
SRC2="airsbgcheck.cdk90 airsch.cdk90 iasibgcheck.cdk90 iasich.cdk90 mod4dv.cdk90 modgps00base.cdk90 modgps01ctmath.cdk90 modgps02wgs84const.cdk90 modgps03diff.cdk90 modgps04profile.cdk90 modgps06gravity.cdk90 out_airs.cdk90 out_iasi.cdk90 ozoneclim.cdk90 qc_profiles.cdk90"
SRC2="$SRC2 mod_tovs.cdk90 modgps05refstruct.cdk90 modgps07geostruct.cdk90 modgps08refop.cdk90 modgps09bend.cdk90"
SRC2="$SRC2 avhrr_var_mod.cdk90 common_iasi.cdk90"
s.compile $INCLUDES $COMPF -O -src $SRC2 > listing2 2>&1
grep fail listing2
if [ $? = "0" ] ; then exit ; fi

echo "compiling the remaining new modules"
SRC3="masks_mod.ftn90"
s.compile $INCLUDES $COMPF -O -src $SRC3 > listing3 2>&1
grep fail listing3
if [ $? = "0" ] ; then exit ; fi

echo "compiling remaining ftn ftn90..."
filelist=""
for i in *.ftn*
do
  xx=`echo $i |grep -v _mod.ftn` 
  filelist="$filelist $xx"
done
s.compile $INCLUDES $COMPF -O -src $filelist > listing4 2>&1
grep fail listing4
if [ $? = "0" ] ; then exit ; fi

echo "building the executable..."
s.compile -O -abi $ABI $COMPF $INCLUDES -libpriv -libpath $LIBPATH2 -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o 3dvar_p7.abs$ABSTAG > listing5 2>&1

rm -f *.ftn* *.cdk* *.h

echo "FINISHED COMPILATION AT:"
date
