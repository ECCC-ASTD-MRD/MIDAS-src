#!/bin/ksh

set -e

mode=$1
nompi=$2

if [ "$mode" == "" ] ; then
  echo " "
  echo "------------------------------------------------------- "
  echo "WARNING: no compilation mode specified, assuming 'full'"
  echo "------------------------------------------------------- "
  echo " "
  mode=full
fi

if [ $mode == full ] ; then
  echo
  echo "----------------------------------------"
  echo " >>> Full Compilation"
  echo "----------------------------------------"
  echo
elif [ $mode == abs ] ; then
  echo
  echo "----------------------------------------"
  echo " >>> Building the Executable only"
  echo "----------------------------------------"
  echo
else
  echo
  echo "----------------------------------------"
  echo " >>> Compiling only routine : $1"
  echo "----------------------------------------"
  echo
fi

#nompi=$1
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
if [ "${BASE_ARCH}" = "AIX-powerpc7" ];then
    FOPTMIZ=2
elif [ "${BASE_ARCH}" = "Linux_x86-64" ];then
    FOPTMIZ=1
else
    echo "This platform 'ARCH=${ARCH}' is not supported.  Only 'AIX-powerpc7' and 'Linux_x86-64' are."
    exit 1
fi
cd ../
mkdir -p compiledir
cd compiledir
rm -f *.o *.f *.f90 *.mod

# automatically set the global revision number in comct0.cdk by
# replacing the string XXXXX with the actual revision number
revpath=$(ssh pollux "cd $trunkdir; svn info | awk '/^URL/ {print \$2}'")
revnum=$(ssh pollux "cd $trunkdir;  svnversion")
echo " " 
echo "-----------------------"
echo "Revision number='$revnum' '$revpath'"
echo "-----------------------"
echo " "

cat ${trunkdir}/toplevelcontrol_mod.ftn90_template |sed "s!XXXXX!${revnum} ${revpath}!g" > toplevelcontrol_mod.ftn90


compiledir=$PWD

#----------------------------------------------------------------
#  Set up dependent librarys and tools. 
#---------------------------------------------------------------
## for s.compile
. ssmuse-sh -d hpcs/201402/01/base
## for the compiler
if [ "${BASE_ARCH}" = "AIX-powerpc7" ];then
    echo "loading compiler hpcs/ext/xlf_13.1.0.10"
    . ssmuse-sh -d hpcs/ext/xlf_13.1.0.10
elif [ "${BASE_ARCH}" = "Linux_x86-64" ];then
    echo "loading hpcs/201402/01/intel13sp1u2"
    . ssmuse-sh -d hpcs/201402/01/intel13sp1u2
else
    echo "This platform 'ARCH=${ARCH}' is not supported.  Only 'AIX-powerpc7' and 'Linux_x86-64' are."
    exit 1
fi

varabs=oavar_${BASE_ARCH}${ABSTAG}

## for rmn_015, rpncomm
echo "loading rpn/libs/15.0"
. ssmuse-sh -d rpn/libs/15.0
## for 'vgrid'
echo "loading cmdn/vgrid/5.0.2/${COMP_ARCH}"
. ssmuse-sh -d cmdn/vgrid/5.0.2/${COMP_ARCH}
## for 'burplib'
echo "loading cmda/base/master/burplib_1.3.3-${COMP_ARCH}_$(ssm platforms | cut -d' ' -f1)"
. ssmuse-sh -p cmda/base/master/burplib_1.3.3-${COMP_ARCH}_$(ssm platforms | cut -d' ' -f1)

## For hpcsperf needed for TMG timings
. ssmuse-sh -d hpcs/exp/aspgjdm/perftools
# For RTTOV 10v1 package... 
echo "loading arma/rttov/10v1"
. ssmuse-sh -d arma/rttov/10v1
#-----------------------------------------------------------------------------

LIBAPPL="rttov10.2.0_coef_io rttov10.2.0_main rttov10.2.0_other burp_module descrip $MPILIB"
if [ "${BASE_ARCH}" = "AIX-powerpc7" ];then
    LIBSYS="hpcsperf lapack-3.4.0 essl mass"
elif [ "${BASE_ARCH}" = "Linux_x86-64" ];then
    LIBSYS="hpcsperf lapack blas"
else
    echo "This platform 'ARCH=${ARCH}' is not supported.  Only 'AIX-powerpc7' and 'Linux_x86-64' are."
    exit 1
fi
LIBRMN="rmn_015"
COMPF_NOC="-openmp $MPIKEY "
COMPF="$COMPF_NOC"
if [ $mode == full ] ; then
 rm -f  *.o *.mod *.cdk* *.h *.ftn* *.f *.f90
 cd ${trunkdir};          ls -1F | grep -v '/' | grep -v "*" | grep -v "@" | cpio -pl $compiledir ; cd $compiledir
 cd ${trunkdir}/bgcheck;  ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir
 cd ${trunkdir}/shared;   ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir
 cd ${trunkdir}/modulopt; ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir

 rm -f *.ftn~ *.ftn90~

 echo "STARTING COMPILATION AT:"
 date

# remove enkf_pturb.ftn main program from compilation directory
 rm -f enkf_pturb.ftn

# Compile the subroutines...
  echo "compiling low-level independent modules"
  SRC0="toplevelcontrol_mod.ftn90"
  SRC0="$SRC0 mathphysconstants_mod.ftn90 earthconstants_mod.ftn90 mpi_mod.ftn90 mpivar_mod.ftn90 bufr_mod.ftn90 codtyp_mod.ftn90"
  SRC0="$SRC0 physicsfunctions_mod.ftn90 obsspacedata_mod.ftn90 horizontalcoord_mod.ftn90 timecoord_mod.ftn90 verticalcoord_mod.ftn90"
  s.compile $COMPF -O ${FOPTMIZ} -src $SRC0 > listing1 2>&1
  status=1
  grep fail listing1 || status=0
  if [ "${status}" -ne 0 ] ; then exit 1; fi

  echo "compiling modulopt (n1qn3)"
  SRC0="dcube.ftn ddd.ftn ddds.ftn dystbl.ftn mupdts.ftn n1qn3.ftn n1qn3a.ftn nlis0.ftn"
  s.compile $COMPF_NOC  -O ${FOPTMIZ} -src $SRC0 > listing0 2>&1
  status=1
  grep fail listing0 || status=0
  if [ "${status}" -ne 0 ] ; then exit 1; fi
  rm -f $SRC0

  echo "compiling analysis grid modules"
  SRC0="gaussgrid_mod.ftn90 windrotation_mod.ftn90 lamanalysisgrid_mod.ftn90"
  s.compile $COMPF  -O ${FOPTMIZ} -src $SRC0 > listing2 2>&1
  status=1
  grep fail listing2 || status=0
  if [ "${status}" -ne 0 ] ; then exit 1; fi

  echo "compiling most of the new modules"
  SRC1="controlvector_mod.ftn90 rmatrix_mod.ftn90 hir_chans_mod.ftn90 tovs_nl_mod.ftn90 tovs_lin_mod.ftn90 varnamelist_mod.ftn90 columndata_mod.ftn90 multi_ir_bgck_mod.ftn90"
  SRC1="$SRC1 emissivities_mod.ftn90 fft_mod.ftn90 globalspectraltransform_mod.ftn90"
  SRC1="$SRC1 lamspectraltransform_mod.ftn90  gridstatevector_mod.ftn90"
  SRC1="$SRC1 bmatrixensemble_mod.ftn90 bmatrixhi_mod.ftn90 lambmatrixhi_mod.ftn90"
  SRC1="$SRC1 bmatrix_mod.ftn90 minimization_mod.ftn90"
  SRC1="$SRC1 ozoneclim_mod.ftn90 tovs_extrap_mod.ftn90"
  SRC1="$SRC1 burpfiles_mod.ftn90 obsspacediag_mod.ftn90 observation_erreurs_mod.ftn90"

  s.compile $COMPF  -O ${FOPTMIZ} -src $SRC1 > listing3 2>&1
  status=1
  grep fail listing3 || status=0
  if [ "${status}" -ne 0 ] ; then exit 1; fi

  echo "compiling burp_read module"
  SRC1="burp_read_mod.ftn90 burp_functions.ftn90 selectb.ftn90 update_burpfiles.ftn90"
  s.compile $COMPF  -O ${FOPTMIZ} -src $SRC1 > listing4 2>&1
  status=1
  grep fail listing4 || status=0
  if [ "${status}" -ne 0 ] ; then exit 1; fi
  
  echo "compiling the GPS modules (cdk90)..."
  SRC2="modgps00base.cdk90 modgps01ctmath.cdk90 modgps01ctphys.cdk90 modgps02wgs84const.cdk90 modgps02wgs84grav.cdk90 modgps03diff.cdk90 modgps04profile.cdk90"
  SRC2="$SRC2 modgps05refstruct.cdk90 modgps07geostruct.cdk90 modgps08refop.cdk90 modgps09bend.cdk90 modgps04profilezd.cdk90"
  SRC2="$SRC2 modgps08ztdop.cdk90"
  s.compile $COMPF  -O ${FOPTMIZ} -src $SRC2 > listing5 2>&1
  status=1
  grep fail listing5 || status=0
  if [ "${status}" -ne 0 ] ; then exit 1; fi

  echo "compiling some more modules..."
  SRC2="modgpsro_mod.ftn90 modgpsztd_mod.ftn90 filterobs_mod.ftn90 writeincrement_mod.ftn90 obsoperators_mod.ftn90"
  s.compile $COMPF  -O ${FOPTMIZ} -src $SRC2 > listing6 2>&1
  status=1
  grep fail listing6 || status=0
  if [ "${status}" -ne 0 ] ; then exit 1; fi

  echo "compiling remaining ftn ftn90..."
  filelist=""
  for i in *.ftn *.ftn90 ; do
      if [[ "${i}" != *_mod.ftn* ]]; then
	  filelist="$filelist ${i}"
      fi
  done
  s.compile $COMPF  -O ${FOPTMIZ} -src $filelist > listing7 2>&1
  status=1
  grep fail listing7 || status=0
  if [ "${status}" -ne 0 ] ; then exit 1; fi

  echo "building the executable..."
  s.compile $COMPF  -O ${FOPTMIZ} -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o ${varabs}.Abs > listing8 2>&1
  status=1
  grep -i ERROR listing? || status=0
  if [ "${status}" -ne 0 ] ; then
      echo "ERROR found: STOP"
      exit 1
  fi

  rm -f *.ftn* *.f *.f90

  echo "FINISHED COMPILATION AT:"
  date

elif [ $mode == abs ] ; then

  rm -f 3dvar_p7.abs$ABSTAG

  echo
  echo "building the executable..."
  echo
  s.compile $COMPF  -O ${FOPTMIZ} -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o ${varabs}.Abs

else

  if [ -f $trunkdir/$mode ] ; then
    file=`basename $mode`
    rm -f $file
    cp $trunkdir/$mode .
    s.compile $COMPF  -O ${FOPTMIZ} -src $file
  else
    echo "File $trunkdir/$mode does NOT exist. Stop"
    exit 1
  fi

fi

echo "FINISHED COMPILATION AT:"
date
echo "The program can be found here: ${PWD}/${varabs}.Abs"