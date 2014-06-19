#!/bin/ksh

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

if [ "$nompi" = "NOMPI" -o "$nompi" = "nompi" ] ; then
  echo " !!! Compiling for a NON-MPI executable !!! "
  echo ""
  MPILIB="rpn_commstubs_40007 rpn_comm_40007"
  MPIKEY=""
  ABSTAG="_NOMPI"
else
  echo " !!! Compiling for an MPI executable !!!"
  echo ""
  MPILIB="rpn_comm_40007"
  MPIKEY="-mpi"
  ABSTAG=""
fi

enkfpturbdir=$PWD
trunkdir=$PWD/../

# automatically set the global revision number in toplevelcontrol_mod.ftn90 by
# replacing the string XXXXX with the actual revision number
revpath=$(ssh pollux "cd $trunkdir; svn info | awk '/^URL/ {print \$2}'")
revnum=$(ssh pollux "cd $trunkdir;  svnversion")
echo " " 
echo "-----------------------"
echo "Revision number='$revnum' '$revpath'"
echo "-----------------------"
echo " "
cat ${trunkdir}/toplevelcontrol_mod.ftn90_template |sed "s!XXXXX!${revnum} ${revpath}!g" > toplevelcontrol_mod.ftn90

ARMNLIB=${ARMNLIB:-/home/dormrb02/ibmenv/armnlib}

VAR3D_VERSION="11.2.1"
LIBAPPL="rttov10.2.0_coef_io rttov10.2.0_main rttov10.2.0_other burp_module modelutils_base descrip $MPILIB "

LIBSYS="lapack blas mass"
LIBRMN="rmn_014_rc2"
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

# Create and Move to compilation directory
mkdir -p ../../compiledir_enkf_pturb
cd ../../compiledir_enkf_pturb
compiledir=$PWD

## To access 'rmn_014_rc2'
. /ssm/net/hpcs/shortcuts/ssmuse_ssm_v10.sh 
. ssmuse-sh -d /ssm/net/rpn/libs/201309/01
. s.ssmuse.dot Xlf13.110
. s.ssmuse.dot devtools
. s.ssmuse.dot ENV/d/x/modelutils/modelutils_1.1.0-a8
. s.ssmuse.dot CMDN/vgrid/4.4.1
. s.ssmuse.dot rpn_comm

if [ $mode == full ] ; then

  rm -f *.o *.mod *.cdk* *.h *.ftn* *.f *.f90

  # Create a local copy of the source code
  cd ${trunkdir};          ls -1F | grep -v '/' | grep -v "*" | grep -v "@" | cpio -pl $compiledir ; cd $compiledir
  cd ${trunkdir}/shared;   ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir
  cd ${trunkdir}/enkf_pturb; ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir
  rm -f *.ftn~ *.ftn90~

  echo "STARTING COMPILATION AT:" 
  date

  # Compile the subroutines...
  echo "compiling low-level independent modules"
  SRC0="toplevelcontrol_mod.ftn90"
  SRC0="$SRC0 mathphysconstants_mod.ftn90 earthconstants_mod.ftn90 mpi_mod.ftn90 mpivar_mod.ftn90 bufr_mod.ftn90 codtyp_mod.ftn90"
  SRC0="$SRC0 physicsfunctions_mod.ftn90 obsspacedata_mod.ftn90 horizontalcoord_mod.ftn90 timecoord_mod.ftn90 verticalcoord_mod.ftn90"
  s.compile $INCLUDES $COMPF -O -src $SRC0 > listing1 2>&1
  grep fail listing1
  if [ $? = "0" ] ; then exit ; fi

  echo "compiling analysis grid modules"
  SRC0="gaussgrid_mod.ftn90 windrotation_mod.ftn90 lamanalysisgrid_mod.ftn90"
  s.compile $INCLUDES $COMPF -O -src $SRC0 > listing2 2>&1
  grep fail listing2
  if [ $? = "0" ] ; then exit ; fi

  echo "compiling most of the new modules"
  SRC1="controlvector_mod.ftn90 fft_mod.ftn90"
  SRC1="$SRC1 globalspectraltransform_mod.ftn90 varnamelist_mod.ftn90"
  SRC1="$SRC1 lamspectraltransform_mod.ftn90 gridstatevector_mod.ftn90"
  SRC1="$SRC1 bmatrixensemble_mod.ftn90 bmatrixhi_mod.ftn90 lambmatrixhi_mod.ftn90"
  SRC1="$SRC1 bmatrix_mod.ftn90 writeincrement_mod.ftn90"
  s.compile $INCLUDES $COMPF -O -src $SRC1 > listing3 2>&1
  grep fail listing3
  if [ $? = "0" ] ; then exit ; fi

  echo "compiling remaining ftn ftn90..."
  filelist="abort.ftn utils_3dvar.ftn getstamplist.ftn90 matsqrt.ftn getfldprm.ftn getfldprm2.ftn"
  filelist="$filelist getstepobsindex.ftn90 subasic_gd.ftn matapat.ftn vintgd.ftn90 initgdg2.ftn90 gasdev.ftn"
  filelist="$filelist enkf_pturb.ftn90"
  s.compile $INCLUDES $COMPF -O -src $filelist > listing7 2>&1
  grep fail listing7
  if [ $? = "0" ] ; then exit ; fi

  echo "building the executable..."
  s.compile -O -abi $ABI $COMPF $INCLUDES -libpriv -libpath $LIBPATH2 -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o enkf_pturb_p7.abs$ABSTAG > listing8 2>&1

  grep -i ERROR listing?
  if [ $? = "0" ] ; then echo "ERROR found: STOP" ; exit ; fi

  rm -f *.ftn* *.f *.f90

  echo "FINISHED COMPILATION AT:"
  date

elif [ $mode == abs ] ; then

  rm -f enkf_pturb_p7.abs$ABSTAG

  echo
  echo "building the executable..."
  echo
  s.compile -O -abi $ABI $COMPF $INCLUDES -libpriv -libpath $LIBPATH2 -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o enkf_pturb_p7.abs$ABSTAG

else

  if [ -f $enkfpturbdir/$mode ] ; then
    file=`basename $mode`
    rm -f $file
    cp $enkfpturbdir/$mode .
    s.compile $INCLUDES $COMPF -O -src $file
  else
    echo "File $enkfpturbdir/$mode does NOT exist. Stop"
    exit 1
  fi

fi
