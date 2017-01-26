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

echo " Compiling for an MPI executable"
echo ""
MPILIBDIR="-libpath /users/dor/arma/gr3/userlibs/AIX-powerpc7/xlf13" # JFC : Mesure temporaire pour avoir acces a
MPILIB="rpn_comm_adj_halo8"                                          #       la S-R RPN_COMM_adj_halo8 de M. Valin
MPIKEY="-mpi"
ABSTAG=""

# replacing the string XXXXX with the actual revision number
revnum=$(git describe --always --dirty=_M 2>/dev/null || ssh pollux "cd $trunkdir; git describe --always --dirty=_M" 2>/dev/null || echo unkown revision)
echo "..."
echo "... > Revision Number = '$revnum'"

enkfpturbdir=$PWD
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

## For hpcsperf needed for TMG timings
. ssmuse-sh -d hpcs/exp/aspgjdm/perftools

LIBAPPL="descrip $MPILIB"

LIBSYS="essl mass"
LIBRMN="rmn"
LIBEXTRA="hpcsperf lapack-3.4.0"
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
  sed "s!XXXXX!${revnum}!g" ${trunkdir}/toplevelcontrol_mod.ftn90_template > toplevelcontrol_mod.ftn90

  # Create a local copy of the source code
  trunkfiles="mpi_mod.ftn90 mpivar_mod.ftn90 varabort.ftn90 physicsfunctions_mod.ftn90 controlvector_mod.ftn90 \
            globalspectraltransform_mod.ftn90 lamanalysisgrid_mod.ftn90 localizationfunction_mod.ftn90 \
            gridstatevector_mod.ftn90 maincompileswitch.inc varnamelist_mod.ftn90 lambmatrixhi_mod.ftn90 \
            utils_3dvar.ftn lamspectraltransform_mod.ftn90 \
            verticalcoord_mod.ftn90 horizontalcoord_mod.ftn90 dsyev2.ftn timecoord_mod.ftn90 \
            columndata_mod.ftn90 bmatrixensemble_mod.ftn90 bmatrixhi_mod.ftn90 bmatrix_mod.ftn90 \
            matsqrt.ftn getfldprm2.ftn randomnumber_mod.ftn90 variabletransforms_mod.ftn90 filterresponsefunction_mod.ftn90"

  cd ${trunkdir}
  cp -f enkf_pturb/main_enkf_pturb.ftn90 ${compiledir}
  cp -f ${trunkfiles} ${compiledir}
  cd ${trunkdir}/shared; ls -1F | grep -v '/' | grep -v "*" | cpio -pl ${compiledir}
  cd ${compiledir}

  echo "STARTING COMPILATION AT:" 
  date

  # Compile the subroutines...
  echo "compiling low-level independent modules"
  SRC0="toplevelcontrol_mod.ftn90 randomnumber_mod.ftn90"
  SRC0="$SRC0 mathphysconstants_mod.ftn90 earthconstants_mod.ftn90 mpi_mod.ftn90 mpivar_mod.ftn90 filterresponsefunction_mod.ftn90"
  SRC0="$SRC0 bufr_mod.ftn90 physicsfunctions_mod.ftn90 horizontalcoord_mod.ftn90 obsspacedata_mod.ftn90"
  s.compile $COMPF -src $SRC0 > listing0a 2>&1
  grep fail listing0a
  if [ $? = "0" ] ; then exit ; fi

  echo "compiling analysis grid modules"
  SRC0="lamanalysisgrid_mod.ftn90"
  s.compile $COMPF -src $SRC0 > listing0b 2>&1
  grep fail listing0b
  if [ $? = "0" ] ; then exit ; fi

  echo "compiling most of the new modules"
  SRC1="localizationfunction_mod.ftn90 controlvector_mod.ftn90 varnamelist_mod.ftn90 timecoord_mod.ftn90"
  SRC1="$SRC1 globalspectraltransform_mod.ftn90 verticalcoord_mod.ftn90"
  SRC1="$SRC1 lamspectraltransform_mod.ftn90 gridstatevector_mod.ftn90"
  SRC1="$SRC1 columndata_mod.ftn90 lambmatrixhi_mod.ftn90"
  SRC1="$SRC1 variabletransforms_mod.ftn90 bmatrixensemble_mod.ftn90 bmatrixhi_mod.ftn90"
  SRC1="$SRC1 bmatrix_mod.ftn90"
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
  s.compile $COMPF ${MPILIBDIR} -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o enkf_pturb_p7.abs$ABSTAG > listing8 2>&1

  grep -i ERROR listing?
  if [ $? = "0" ] ; then echo "ERROR found: STOP" ; exit ; fi

  rm -f *.ftn* 

  echo "FINISHED COMPILATION AT:"
  date

elif [ $mode == abs ] ; then

  rm -f enkf_pturb_p7.abs$ABSTAG

  echo
  echo "building the executable..."
  echo
  s.compile $COMPF ${MPILIBDIR} -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o enkf_pturb_p7.abs$ABSTAG

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
