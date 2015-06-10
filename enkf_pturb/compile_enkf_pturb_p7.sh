#!/bin/ksh

mode=$1

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

echo " !!! Compiling for an MPI executable !!!"
echo ""
MPILIB="rpn_comm_40511"
MPIKEY="-mpi"
ABSTAG=""

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
rm -f ${trunkdir}/toplevelcontrol_mod.ftn90
cat ${trunkdir}/toplevelcontrol_mod.ftn90_template |sed "s!XXXXX!${revnum} ${revpath}!g" > toplevelcontrol_mod.ftn90

# Load the appropriate librairies
. ssmuse-sh -d hpcs/13b/04/base
## for the compiler
. ssmuse-sh -d hpcs/ext/xlf_13.1.0.10
## for rmn_014, lapack_3.4.0, rpncomm
. ssmuse-sh -d rpn/libs/15.0
. s.ssmuse.dot ENV/d/x/modelutils/modelutils_1.1.0-a8
. ssmuse-sh -d /ssm/net/cmdn/vgrid/5.3.0-a2/xlf13
. s.ssmuse.dot cmda
. ssmuse-sh -d arma/rttov/10v1

## For hpcsperf needed for TMG timings
. ssmuse-sh -d hpcs/exp/aspgjdm/perftools

VAR3D_VERSION="11.2.1"
LIBAPPL="rttov10.2.0_coef_io rttov10.2.0_main rttov10.2.0_other burp_module modelutils_base descrip $MPILIB "


LIBSYS="essl mass"
LIBRMN="rmn_015"
LIBEXTRA="hpcsperf lapack-3.4.0"
MODBURP="BURP1.3"
DEFINE="-DNEC=nec -DIBM=ibm"
COMPF_NOC="-openmp $MPIKEY -O"
#COMPF="$COMPF_NOC"
COMPF="$COMPF_NOC -debug DEBUG -optf=-C "

# Create and Move to compilation directory
mkdir -p ../../compiledir_enkf_pturb
cd ../../compiledir_enkf_pturb
compiledir=$PWD

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
  SRC0="$SRC0 physicsfunctions_mod.ftn90 obsspacedata_mod.ftn90 horizontalcoord_mod.ftn90 timecoord_mod.ftn90 verticalcoord_mod.ftn90 dsyev2.ftn"
  s.compile $COMPF -src $SRC0 > listing1 2>&1
  grep fail listing1
  if [ $? = "0" ] ; then exit ; fi

  echo "compiling analysis grid modules"
  SRC0="gaussgrid_mod.ftn90 windrotation_mod.ftn90 lamanalysisgrid_mod.ftn90"
  s.compile $COMPF -src $SRC0 > listing2 2>&1
  grep fail listing2
  if [ $? = "0" ] ; then exit ; fi

  echo "compiling most of the new modules"
  SRC1="controlvector_mod.ftn90 fft_mod.ftn90"
  SRC1="$SRC1 globalspectraltransform_mod.ftn90 varnamelist_mod.ftn90"
  SRC1="$SRC1 lamspectraltransform_mod.ftn90 gridstatevector_mod.ftn90"
  SRC1="$SRC1 bmatrixensemble_mod.ftn90 bmatrixhi_mod.ftn90 lambmatrixhi_mod.ftn90"
  SRC1="$SRC1 bmatrix_mod.ftn90 writeincrement_mod.ftn90"
  s.compile $COMPF -src $SRC1 > listing3 2>&1
  grep fail listing3
  if [ $? = "0" ] ; then exit ; fi

  echo "compiling remaining ftn ftn90..."
  filelist="abort.ftn utils_3dvar.ftn getstamplist.ftn90 matsqrt.ftn getfldprm.ftn getfldprm2.ftn"
  filelist="$filelist getstepobsindex.ftn90 matapat.ftn vintgd.ftn90 initgdg2.ftn90 gasdev.ftn"
  filelist="$filelist main_enkf_pturb.ftn90"
  s.compile $COMPF -src $filelist > listing7 2>&1
  grep fail listing7
  if [ $? = "0" ] ; then exit ; fi

  echo "building the executable..."
  s.compile $COMPF -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o enkf_pturb_p7.abs$ABSTAG > listing8 2>&1

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
  s.compile $COMPF -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o enkf_pturb_p7.abs$ABSTAG

else

  if [ -f $enkfpturbdir/$mode ] ; then
    file=`basename $mode`
    rm -f $file
    cp $enkfpturbdir/$mode .
    s.compile $COMPF -src $file
  else
    echo "File $enkfpturbdir/$mode does NOT exist. Stop"
    exit 1
  fi

fi
