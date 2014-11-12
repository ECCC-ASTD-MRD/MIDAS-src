#!/bin/ksh

set -e

nompi=$1
if [ "$nompi" = "NOMPI" -o "$nompi" = "nompi" ] 
then
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "!!Compiling for a NON-MPI executable!!"
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  MPILIB="rpn_commstubs_40511 rpn_comm_40511"
  MPIKEY=""
  ABSTAG="_nompi"
else
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "!!Compiling for an MPI executable!!"
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  MPILIB="rpn_comm_40511"
  MPIKEY="-mpi"
  ABSTAG=""
fi

trunkdir=$PWD
if [ "$BASE_ARCH" = "AIX-powerpc7" ];then
  compiler="xlf13"
  FOPTMIZ=2 
else
  # At the moment, the defulat on Linux is "Intel13sp1" 
  compiler="intel13sp1u2" 
  FOPTMIZ=1 
fi
cd ../
mkdir -p compiledir
cd compiledir
rm -f *.o *.f *.f90 *.mod

# automatically set the global revision number in comct0.cdk by
# replacing the string XXXXX with the actual revision number
revpath=$(ssh pollux "cd $trunkdir; svn info | awk '/^URL/ {print \$2}'")
revnum=$(ssh pollux "cd $trunkdir;  svnversion")
echo "Revision number='$revnum' '$revpath'"
cat ${trunkdir}/comct0_template.cdk |sed "s!XXXXX!${revnum} ${revpath}!g" > comct0.cdk

compiledir=$PWD

#----------------------------------------------------------------
#  Set up dependent librarys and tools. 
#---------------------------------------------------------------
## for s.compile
. ssmuse-sh -d hpcs/201402/01/base
## for the compiler
if [ "$compiler" = "xlf13" ];then
    . ssmuse-sh -d hpcs/ext/xlf_13.1.0.10
else
    . ssmuse-sh -d hpcs/201402/01/${compiler}
fi 

varabs=oavar_${BASE_ARCH}${ABSTAG}

# for msg in VGRID .. 
if [ "$compiler" = "xlf13" ];then
    LIB_MODELUTILS=/home/ordenv/ssm-domains9/ENV/modelutils/modelutils_1.2.2/aix61-ppc-64/lib/AIX-powerpc7/xlf13
else
    LIB_MODELUTILS=/home/ordenv/ssm-domains9/ENV/modelutils/modelutils_1.2.2/linux26-x86-64/lib/Linux_x86-64/intel13sp1
fi 
## for rmn_015, lapack_3.4.0, rpncomm
. ssmuse-sh -d rpn/libs/15.0
## for 'vgrid'
. ssmuse-sh -d cmdn/vgrid/5.0.2/${compiler}
## for 'burplib'
. ssmuse-sh -p cmda/base/master/burplib_1.3.3-${compiler}_$(ssm platforms | cut -d' ' -f1)

## For hpcsperf needed for TMG timings
. ssmuse-sh -d hpcs/exp/aspgjdm/perftools
# For RTTOV 10v1 package... 
echo "loading arma/rttov/10v1"
if [ "$compiler" = "xlf13" ];then
    . ssmuse-sh -d arma/rttov/10v1
else
    . ssmuse-sh -d arma/rttov/10v1
fi 
#-----------------------------------------------------------------------------

LIBAPPL="rttov10.2.0_coef_io rttov10.2.0_main rttov10.2.0_other burp_module descrip $MPILIB"
if [ "$compiler" = "xlf13" ]; then
    LIBSYS="hpcsperf essl mass"
else
    LIBSYS="hpcsperf lapack blas"
fi
LIBRMN="rmn_015"
#LIBEXTRA="hpcsperf"
MODBURP="BURP1.3"
if [ "${compiler}" = "xlf13" ]; then
    DEFINE="-DNEC=nec -DIBM=ibm"
fi
COMPF_NOC="-openmp $MPIKEY "
COMPF="$COMPF_NOC"

cd ${trunkdir};          ls -1F | grep -v '/' | grep -v "*" | grep -v "@" | cpio -pl $compiledir ; cd $compiledir
cd ${trunkdir}/bgcheck;  ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir
cd ${trunkdir}/shared;   ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir
cd ${trunkdir}/modulopt; ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir

rm -f *.ftn~ *.ftn90~

echo "STARTING COMPILATION AT:"
date

# remove enkf_pturb.ftn main program from compilation directory
rm -f enkf_pturb.ftn

echo "compiling modulopt (n1qn3) [ALSO DSYEV WHICH SHOULD NOT BE HERE!]"
SRC0="dcube.ftn ddd.ftn ddds.ftn dsyev.ftn dystbl.ftn mupdts.ftn n1qn3.ftn n1qn3a.ftn nlis0.ftn"
s.compile $COMPF_NOC -O ${FOPTMIZ} -src $SRC0 > listingm 2>&1
status=0
grep fail listingm || status=1
if [ "${status}" -eq 0 ] ; then exit 1; fi
rm -f $SRC0

echo "compiling low-level independent modules"
SRC0="mathphysconstants_mod.ftn90 earthconstants_mod.ftn90 mpi_mod.ftn90 bufr_mod.ftn90 physicsfunctions_mod.ftn90 gaussgrid_mod.ftn90"
s.compile $COMPF -O ${FOPTMIZ} -src $SRC0 > listing0 2>&1
status=0
grep fail listing0 || status=1
if [ "${status}" -eq 0 ] ; then exit 1; fi

echo "compiling most of the new modules"
SRC1="controlvector_mod.ftn90 hir_chans_mod.ftn90 tovs_mod.ftn90 emissivities_mod.ftn90 fft_mod.ftn90"
SRC1="$SRC1 globalspectraltransform_mod.ftn90 obsspacedata_mod.ftn90 random_mod.ftn90 varnamelist_mod.ftn90 verticalcoord_mod.ftn90"
SRC1="$SRC1 columndata_mod.ftn90 gridstatevector_mod.ftn90"
SRC1="$SRC1 bmatrixensemble_mod.ftn90 bmatrixhi_mod.ftn90"
SRC1="$SRC1 bmatrix_mod.ftn90 minimization_mod.ftn90"
SRC1="$SRC1 multi_ir_bgck_mod.ftn90 ozoneclim_mod.ftn90"

s.compile $COMPF -O ${FOPTMIZ} -src $SRC1 > listing1 2>&1
status=0
grep fail listing1 || status=1
if [ "${status}" -eq 0 ] ; then exit 1; fi

echo "compiling burp_read module"
SRC1="burp_read_mod.ftn90 burp_functions.ftn90 selectb.ftn90 update_burpfiles.ftn90"
s.compile $COMPF -O ${FOPTMIZ} -src $SRC1 > listing_burp 2>&1
status=0
grep fail listing_burp || status=1
if [ "${status}" -eq 0 ] ; then exit 1; fi

echo "compiling the old modules (cdk90)..."
SRC2="modgps00base.cdk90 modgps01ctmath.cdk90 modgps01ctphys.cdk90 modgps02wgs84const.cdk90 modgps02wgs84grav.cdk90 modgps03diff.cdk90 modgps04profile.cdk90"
SRC2="$SRC2 modgps05refstruct.cdk90 modgps07geostruct.cdk90 modgps08refop.cdk90 modgps09bend.cdk90 modgpsro_mod.ftn90 modgps04profilezd.cdk90"
SRC2="$SRC2 modgps08ztdop.cdk90 modgpsztd_mod.ftn90"
s.compile $COMPF -O ${FOPTMIZ} -src $SRC2 > listing2 2>&1
status=0
grep fail listing2 || status=1
if [ "${status}" -eq 0 ] ; then exit 1; fi

echo "compiling remaining ftn ftn90..."
filelist=""
for i in *.ftn *.ftn90
do
  xx=`echo $i | grep -v _mod.ftn || true`
  filelist="$filelist $xx"
done
s.compile $COMPF -O ${FOPTMIZ}  -src $filelist > listing4 2>&1
status=0
grep fail listing4 || status=1
if [ "${status}" -eq 0 ] ; then exit 1; fi

echo "building the executable..."
s.compile -O ${FOPTMIZ}  $COMPF -libpath ${LIB_MODELUTILS} -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o ${varabs}.Abs > listing5 2>&1 || echo "Compilation aborted"

status=0
grep -iE 'ERROR|ERREUR' listing? || status=1
if [ "${status}" -eq 0 ] ; then
    echo "ERROR found: STOP"
    exit 1
fi

rm -f *.ftn* *.cdk* *.h

echo "FINISHED COMPILATION AT:"
date
echo "The program can be found here: ${PWD}/${varabs}.Abs"
