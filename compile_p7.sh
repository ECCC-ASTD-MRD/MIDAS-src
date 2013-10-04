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

cd ../
mkdir -p compiledir
cd compiledir
rm -f *.o *.f *.f90 *.mod

# automatically set the global revision number in comct0.cdk by
# replacing the string XXXXX with the actual revision number
revnum=`ssh pollux "cd $trunkdir ; svnversion"`
echo Revision number= $revnum
cat ${trunkdir}/comct0_template.cdk |sed "s/XXXXX/${revnum}/g" > comct0.cdk

compiledir=$PWD

. s.ssmuse.dot Xlf13.108
. s.ssmuse.dot rmnlib-dev
. s.ssmuse.dot devtools

. /ssm/net/hpcs/shortcuts/ssmuse_ssm_v10.sh 

. s.ssmuse.dot CMDN/vgrid/3.4.0
. s.ssmuse.dot rpn_comm


VAR3D_VERSION="11.2.1"
LIBAPPL="rttov10.2.0_coef_io rttov10.2.0_main rttov10.2.0_other burp_module descrip $MPILIB "

LIBSYS="lapack blas mass"
LIBRMN="rmn_013_rc2"
LIBEXTRA="rtools hpm_r"
MODBURP="BURP1.3"
DEFINE="-DNEC=nec -DIBM=ibm"
ABI="_multi"
COMPF_NOC="-openmp $MPIKEY "
#COMPF="$COMPF_NOC"
COMPF="$COMPF_NOC -debug DEBUG -optf=-C "

BASE_INCLUDE="${ARMNLIB}/modeles/ANAL/v_${VAR3D_VERSION}/include/AIX-powerpc7"
INCLUDES="-includes ${BASE_INCLUDE}/${MODBURP} ${ARMNLIB}/modeles/ANAL_shared/rttov10/v1/AIX-powerpc7/xlf13/mod ${ARMNLIB}/modeles/ANAL_shared/rttov10/v1/AIX-powerpc7/xlf13/include"

LIBPATH2="./ $LIBPATH"
LIBPATH2="${ARMNLIB}/lib/AIX/xlf13 $LIBPATH2"
LIBPATH2="${ARMNLIB}/modeles/ANAL/v_${VAR3D_VERSION}/lib/AIX-powerpc7 $LIBPATH2"
LIBPATH2="/home/ordenv/ssm-domains1/ssm-rmnlib-dev/multi/lib/AIX-powerpc7/xlf13 ${ARMNLIB}/modeles/ANAL_shared/rttov10/v1/AIX-powerpc7/xlf13/lib  $LIBPATH2"

echo "LIBPATH2="
echo $LIBPATH2
echo "INCLUDES="
echo $INCLUDES

cd ${trunkdir};          ls -1F | grep -v '/' | grep -v "*" | grep -v "@" | cpio -pl $compiledir ; cd $compiledir
cd ${trunkdir}/bgcheck;  ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir
cd ${trunkdir}/shared;   ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir
cd ${trunkdir}/modulopt; ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir

rm -f *.ftn~ *.ftn90~

# temporarily copy the object files for a preliminary version
# of ezscint that works with the "U" grid (compiled only be AIX):
cp /users/dor/arma/erv/data/ords/yinyang/ezscint/objects_quiet3l/*.o $compiledir

echo "STARTING COMPILATION AT:"
date

# remove enkf_pturb.ftn main program from compilation directory
rm -f enkf_pturb.ftn

echo "compiling modulopt (n1qn3) [ALSO DSYEV WHICH SHOULD NOT BE HERE!]"
SRC0="dcube.ftn ddd.ftn ddds.ftn dsyev.ftn dystbl.ftn mupdts.ftn n1qn3.ftn n1qn3a.ftn nlis0.ftn"
s.compile $INCLUDES $COMPF_NOC -O -src $SRC0 > listingm 2>&1
grep fail listingm
if [ $? = "0" ] ; then exit ; fi
rm -f $SRC0

echo "compiling low-level independent modules"
SRC0="mathphysconstants_mod.ftn90 earthconstants_mod.ftn90 mpi_mod.ftn90 bufr_mod.ftn90 physicsfunctions_mod.ftn90 gaussgrid_mod.ftn90"
s.compile $INCLUDES $COMPF -O -src $SRC0 > listing0 2>&1
grep fail listing0
if [ $? = "0" ] ; then exit ; fi

echo "compiling most of the new modules"
SRC1="controlvector_mod.ftn90 hir_chans_mod.ftn90 tovs_mod.ftn90 emissivities_mod.ftn90 fft_mod.ftn90"
SRC1="$SRC1 globalspectraltransform_mod.ftn90 obsspacedata_mod.ftn90 random_mod.ftn90 varnamelist_mod.ftn90 verticalcoord_mod.ftn90"
SRC1="$SRC1 columndata_mod.ftn90 gridstatevector_mod.ftn90"
SRC1="$SRC1 bmatrixensemble_mod.ftn90 bmatrixhi_mod.ftn90"
SRC1="$SRC1 bmatrix_mod.ftn90 minimization_mod.ftn90"
SRC1="$SRC1 multi_ir_bgck_mod.ftn90 ozoneclim_mod.ftn90"


s.compile $INCLUDES $COMPF -O -src $SRC1 > listing1 2>&1
grep fail listing1
if [ $? = "0" ] ; then exit ; fi

echo "compiling burp_read module"
SRC1="burp_read_mod.ftn90 burp_functions.ftn90 selectb.ftn90 update_burpfiles.ftn90"
s.compile $INCLUDES $COMPF -O -src $SRC1 > listing_burp 2>&1
grep fail listing_burp
if [ $? = "0" ] ; then exit ; fi

echo "compiling the old modules (cdk90)..."
SRC2="modgps00base.cdk90 modgps01ctmath.cdk90 modgps01ctphys.cdk90 modgps02wgs84const.cdk90 modgps02wgs84grav.cdk90 modgps03diff.cdk90 modgps04profile.cdk90"
SRC2="$SRC2 modgps05refstruct.cdk90 modgps07geostruct.cdk90 modgps08refop.cdk90 modgps09bend.cdk90 modgpsro_mod.ftn90 modgps04profilezd.cdk90"
SRC2="$SRC2 modgps08ztdop.cdk90 modgpsztd_mod.ftn90"
s.compile $INCLUDES $COMPF -O -src $SRC2 > listing2 2>&1
grep fail listing2
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
s.compile -O -abi $ABI $COMPF $INCLUDES -libpriv -libpath $LIBPATH2 -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o 3dvar_p7.abs$ABSTAG > listing5 2>&1

grep -i ERROR listing?
if [ $? = "0" ] ; then exit; echo "ERROR found: STOP" ; fi

rm -f *.ftn* *.cdk* *.h

echo "FINISHED COMPILATION AT:"
date
