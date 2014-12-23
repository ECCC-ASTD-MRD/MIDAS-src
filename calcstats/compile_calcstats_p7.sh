#!/bin/ksh

mode=$1
#nompi=$2

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

#if [ "$nompi" = "NOMPI" -o "$nompi" = "nompi" ] ; then
#  echo " !!! Compiling for a NON-MPI executable !!! "
#  echo ""
#  MPILIB="rpn_commstubs_40511 rpn_comm_40511"
#  MPIKEY=""
#  ABSTAG="_NOMPI"
#else
  echo " !!! Compiling for an MPI executable !!!"
  echo ""
  MPILIB="rpn_comm_40511"
  MPIKEY="-mpi"
  ABSTAG=""
#fi

calcstatsdir=$PWD
trunkdir=$PWD/../

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
LIBAPPL="rttov10.2.0_coef_io rttov10.2.0_main rttov10.2.0_other burp_module modelutils_base descrip $MPILIB"

LIBSYS="essl mass"
LIBRMN="rmn_015"
LIBEXTRA="hpcsperf lapack-3.4.0"
MODBURP="BURP1.3"
DEFINE="-DNEC=nec -DIBM=ibm"
COMPF_NOC="-openmp $MPIKEY -O"
#COMPF="$COMPF_NOC"
COMPF="$COMPF_NOC -debug DEBUG -optf=-C "

# Create and Move to compilation directory
compiledir="../../compiledir_calcstats"
echo "COMPILEDIR IS SET TO: $compiledir"
mkdir -p $compiledir
cd $compiledir

compiledir=$PWD

if [ $mode == full ] ; then

  rm -f *.o *.mod *.cdk* *.h *.ftn* *.f *.f90

  # Create a local copy of the source code
  trunkfiles="mpi_mod.ftn90 mpivar_mod.ftn90 abort.ftn physicsfunctions_mod.ftn90 controlvector_mod.ftn90 \
            gaussgrid_mod.ftn90 fft_mod.ftn90 globalspectraltransform_mod.ftn90 lamanalysisgrid_mod.ftn90 \
            gridstatevector_mod.ftn90 maincompileswitch.inc varnamelist_mod.ftn90 \
            utils_3dvar.ftn lamspectraltransform_mod.ftn90 \
            verticalcoord_mod.ftn90 horizontalcoord_mod.ftn90 dsyev2.ftn"

  cd ${trunkdir}
  cp -f calcstats/*.ftn90 ${compiledir}
  cp -f ${trunkfiles} ${compiledir}
  cd ${trunkdir}/shared; ls -1F | grep -v '/' | grep -v "*" | cpio -pl ${compiledir}
  cd ${compiledir}

  echo "STARTING COMPILATION AT:" 
  date

  # Remove enkf_pturb.ftn main program from compilation directory
  rm -f enkf_pturb.ftn

  # Compile the subroutines...
  echo "compiling low-level independent modules"
  SRC0="mathphysconstants_mod.ftn90 earthconstants_mod.ftn90 mpi_mod.ftn90 mpivar_mod.ftn90"
  SRC0="$SRC0 bufr_mod.ftn90 physicsfunctions_mod.ftn90 horizontalcoord_mod.ftn90"
  s.compile $COMPF -src $SRC0 > listing0a 2>&1
  grep fail listing0a
  if [ $? = "0" ] ; then exit ; fi

  echo "compiling analysis grid modules"
  SRC0="gaussgrid_mod.ftn90 lamanalysisgrid_mod.ftn90"
  s.compile $COMPF -src $SRC0 > listing0b 2>&1
  grep fail listing0b
  if [ $? = "0" ] ; then exit ; fi

  echo "compiling most of the new modules"
  SRC1="controlvector_mod.ftn90 fft_mod.ftn90 varnamelist_mod.ftn90"
  SRC1="$SRC1 globalspectraltransform_mod.ftn90 verticalcoord_mod.ftn90"
  SRC1="$SRC1 lamspectraltransform_mod.ftn90 gridstatevector_mod.ftn90"
  SRC1="$SRC1 calcbmatrix_glb_mod.ftn90 calcbmatrix_lam_mod.ftn90"
  s.compile $COMPF -src $SRC1 > listing1 2>&1
  grep fail listing1
  if [ $? = "0" ] ; then exit ; fi

  echo "compiling remaining ftn ftn90..."
  filelist=""
  for i in *.ftn *.ftn90 ; do
    xx=`echo $i |grep -v _mod.ftn` 
    filelist="$filelist $xx"
  done
  s.compile $COMPF -src $filelist > listing4 2>&1
  grep fail listing4
  if [ $? = "0" ] ; then exit ; fi

  echo "building the executable..."
  s.compile $COMPF -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o calcstats_p7.abs$ABSTAG > listing5 2>&1

  grep -i ERROR listing?
  if [ $? = "0" ] ; then exit; echo "ERROR found: STOP" ; fi

  rm -f *.ftn* 

  echo "FINISHED COMPILATION AT:"
  date

elif [ $mode == abs ] ; then

  rm -f calc_stats_p7.abs$ABSTAG

  echo
  echo "building the executable..."
  echo
  s.compile $COMPF -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o calcstats_p7.abs$ABSTAG

else

  if [ -f $calcstatsdir/$mode ] ; then
    file=`basename $mode`
    rm -f $file
    cp $calcstatsdir/$mode .
    s.compile $COMPF -src $file
  else
    echo "File $calcstatsdir/$mode does NOT exist. Stop"
    exit 1
  fi

fi
