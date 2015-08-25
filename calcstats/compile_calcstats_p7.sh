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

echo " Compiling for an MPI executable (even if calcstats in not MPI capable)"
echo ""
MPILIBDIR="-libpath /users/dor/arma/gr3/userlibs/AIX-powerpc7/xlf13" # JFC : Mesure temporaire pour avoir acces a
MPILIB="rpn_comm_adj_halo8"                                          #       la S-R RPN_COMM_adj_halo8 de M. Valin
MPIKEY="-mpi"
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

## for 'vgrid'
echo "loading cmdn/vgrid/5.3.2/${COMP_ARCH}"
. ssmuse-sh -d cmdn/vgrid/5.3.2/${COMP_ARCH}

## for 'burplib'
echo "loading cmda/base/201411/01/${COMP_ARCH}"
. ssmuse-sh -d cmda/base/201411/01/${COMP_ARCH}

## For RTTOV 10v1 package... 
echo "loading arma/rttov/10v1"
. ssmuse-sh -d arma/rttov/10v1

## For hpcsperf needed for TMG timings
. ssmuse-sh -d hpcs/exp/aspgjdm/perftools

VAR3D_VERSION="11.2.1"
LIBAPPL="rttov10.2.0_coef_io rttov10.2.0_main rttov10.2.0_other burp_module descrip $MPILIB"

LIBSYS="essl mass"
LIBRMN="rmn"
LIBEXTRA="hpcsperf lapack-3.4.0"
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
  s.compile $COMPF ${MPILIBDIR} -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o calcstats_p7.abs$ABSTAG > listing5 2>&1

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
  s.compile $COMPF ${MPILIBDIR} -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o calcstats_p7.abs$ABSTAG

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
