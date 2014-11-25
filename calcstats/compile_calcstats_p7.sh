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
  echo " !!! Compiling for a NON-MPI executable !!! "
  echo ""
  MPILIB="rpn_commstubs_40007 rpn_comm_40007"
  MPIKEY=""
  ABSTAG="_NOMPI"
#else
  #echo " !!! Compiling for an MPI executable !!!"
  #echo ""
  #MPILIB="rpn_comm_40007"
  #MPIKEY="-mpi"
  #ABSTAG=""
#fi

calcstatsdir=$PWD
trunkdir=$PWD/../

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
compiledir="../../compiledir_calcstats"
echo "COMPILEDIR IS SET TO: $compiledir"
mkdir -p $compiledir
cd $compiledir

compiledir=$PWD

if [ $mode == full ] ; then

  rm -f *.o *.mod *.cdk* *.h *.ftn* *.f *.f90

  # Load the appropriate librairies
  . /ssm/net/hpcs/shortcuts/ssmuse_ssm_v10.sh 
  . ssmuse-sh -d /ssm/net/rpn/libs/201309/01
  . s.ssmuse.dot Xlf13.110
  . s.ssmuse.dot devtools
  . s.ssmuse.dot ENV/d/x/modelutils/modelutils_1.1.0-a8
  . s.ssmuse.dot CMDN/vgrid/4.4.0-a2
  . s.ssmuse.dot rpn_comm

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
  s.compile $INCLUDES $COMPF -O -src $SRC0 > listing0a 2>&1
  grep fail listing0a
  if [ $? = "0" ] ; then exit ; fi

  echo "compiling analysis grid modules"
  SRC0="gaussgrid_mod.ftn90 lamanalysisgrid_mod.ftn90"
  s.compile $INCLUDES $COMPF -O -src $SRC0 > listing0b 2>&1
  grep fail listing0b
  if [ $? = "0" ] ; then exit ; fi

  echo "compiling most of the new modules"
  SRC1="controlvector_mod.ftn90 fft_mod.ftn90 varnamelist_mod.ftn90"
  SRC1="$SRC1 globalspectraltransform_mod.ftn90 verticalcoord_mod.ftn90"
  SRC1="$SRC1 lamspectraltransform_mod.ftn90 gridstatevector_mod.ftn90"
  SRC1="$SRC1 calcbmatrix_glb_mod.ftn90 calcbmatrix_lam_mod.ftn90"
  s.compile $INCLUDES $COMPF -O -src $SRC1 > listing1 2>&1
  grep fail listing1
  if [ $? = "0" ] ; then exit ; fi

  echo "compiling remaining ftn ftn90..."
  filelist=""
  for i in *.ftn *.ftn90 ; do
    xx=`echo $i |grep -v _mod.ftn` 
    filelist="$filelist $xx"
  done
  s.compile $INCLUDES $COMPF -O -src $filelist > listing4 2>&1
  grep fail listing4
  if [ $? = "0" ] ; then exit ; fi

  echo "building the executable..."
  s.compile -O -abi $ABI $COMPF $INCLUDES -libpriv -libpath $LIBPATH2 -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o calcstats_p7.abs$ABSTAG > listing5 2>&1

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
  s.compile -O -abi $ABI $COMPF $INCLUDES -libpriv -libpath $LIBPATH2 -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o calcstats_p7.abs$ABSTAG

else

  if [ -f $calcstatsdir/$mode ] ; then
    file=`basename $mode`
    rm -f $file
    cp $calcstatsdir/$mode .
    s.compile $INCLUDES $COMPF -O -src $file
  else
    echo "File $calcstatsdir/$mode does NOT exist. Stop"
    exit 1
  fi

fi
