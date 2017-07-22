#!/bin/bash

set -e

mode=$1

. ./compile_commons.sh

echo "..."
echo "...             |========================================|"
echo "... ------------| RANDOMPERT compilation script STARTING |------------"
echo "...             |========================================|"

if [ "$mode" == "" ] ; then
  echo "... "
  echo "... !WARNING! no compilation mode specified, assuming 'full'"
  mode=full
fi

if [ $mode == full ] ; then
  echo "..."
  echo "... > Full Compilation"
elif [ $mode == abs ] ; then
  echo "..."
  echo "... > Building the Executable only"
else
  echo "..."
  echo "... > Compiling only routine : $1"
fi

echo "..."
echo "... > Compiling for an MPI executable"
MPILIB="rpn_comm"

randompertdir=$PWD
trunkdir=$PWD/../

# automatically set the global revision number in toplevelcontrol_mod.ftn90 by
# replacing the string XXXXX with the actual revision number
revnum=$(git describe --always --dirty=_M 2>/dev/null || ssh eccc-ppp1 "cd $trunkdir; git describe --always --dirty=_M" 2>/dev/null || echo unkown revision)
echo "..."
echo "... > Revision Number = '$revnum'"

# Set compiledir
compiledir_main=${COMPILEDIR_OAVAR_MAIN:-".."}
compiledir_ID=${COMPILEDIR_OAVAR_ID:-$revnum}
compiledir=${compiledir_main}/compiledir-randompert-${ORDENV_PLAT}_${compiledir_ID}
mkdir -p $compiledir
cd $compiledir
compiledir=${PWD} # needed when compiledir_main = ".."

echo "..."
echo "... > Compiledir set to ${compiledir}"
if [ ${compiledir_main} != ".." ] ; then
    if [ ! -d  ${trunkdir}/../compiledir-randompert-${ORDENV_PLAT}_${compiledir_ID} ] ; then
	ln -s ${compiledir_main}/compiledir-randompert-${ORDENV_PLAT}_${compiledir_ID} ${trunkdir}/../compiledir-randompert-${ORDENV_PLAT}_${compiledir_ID}
    fi
fi

absdir=${compiledir_main}/oavar_abs
mkdir -p ${absdir}
cd $absdir ; absdir=$PWD ; cd - >/dev/null

varabs=randompert_${ORDENV_PLAT}-${revnum}.Abs

#-----------------------------------------------------------------------------

LIBAPPL="descrip $MPILIB"
LIBSYS="hpcoperf"
LIBRMN=rmnMP

if [ $mode == full ] ; then

  rm -f *.o *.mod *.cdk* *.h *.ftn* *.f *.f90

  # Create a local copy of the source code
  sed "s!XXXXX!${revnum}!g" ${trunkdir}/toplevelcontrol_mod.ftn90_template > toplevelcontrol_mod.ftn90

  # Create a local copy of the source code
  trunkfiles="mpi_mod.ftn90 mpivar_mod.ftn90 physicsfunctions_mod.ftn90 controlvector_mod.ftn90 \
            globalspectraltransform_mod.ftn90 analysisgrid_mod.ftn90 localizationfunction_mod.ftn90 \
            gridstatevector_mod.ftn90 maincompileswitch.inc varnamelist_mod.ftn90 lambmatrixhi_mod.ftn90 \
            utilities_mod.ftn90 lamspectraltransform_mod.ftn90 spectralfilter_mod.ftn90 \
            verticalcoord_mod.ftn90 horizontalcoord_mod.ftn90 timecoord_mod.ftn90 \
            columndata_mod.ftn90 bmatrixensemble_mod.ftn90 bmatrixhi_mod.ftn90 bmatrixchem_mod.ftn90 bmatrix_mod.ftn90 \
            randomnumber_mod.ftn90 variabletransforms_mod.ftn90 ensemblestatevector_mod.ftn90"

  cd ${trunkdir}
  cp -f randompert/main_randompert.ftn90 ${compiledir}/
  cp -f ${trunkfiles} ${compiledir}/
  cp -f shared/*.ftn90 ${compiledir}/
  cp -f shared/*.inc ${compiledir}/

  cd ${compiledir}

  echo "..."
  echo "... > STARTING COMPILATION AT: $(date)"
  echo "..."

  # Compile the subroutines...
  echo "... > Compiling low-level independent modules..."
  SRC0="randomnumber_mod.ftn90 utilities_mod.ftn90 toplevelcontrol_mod.ftn90 bufr_mod.ftn90"
  SRC0="$SRC0 mathphysconstants_mod.ftn90 earthconstants_mod.ftn90 mpi_mod.ftn90 mpivar_mod.ftn90"
  SRC0="$SRC0 physicsfunctions_mod.ftn90 horizontalcoord_mod.ftn90 obsspacedata_mod.ftn90 analysisgrid_mod.ftn90"

  s.compile $COMPF -O ${FOPTMIZ} -src $SRC0 > listing0 2>&1
  status=1
  grep fail listing0 || status=0
  if [ "${status}" -ne 0 ] ; then
      echo "... !! Compilation aborted: check in ${PWD}/listing0 !!"
      exit 1
  fi

  echo "... > Compiling most of the modules..."
  SRC1="localizationfunction_mod.ftn90 controlvector_mod.ftn90 varnamelist_mod.ftn90 timecoord_mod.ftn90"
  SRC1="$SRC1 globalspectraltransform_mod.ftn90 verticalcoord_mod.ftn90"
  SRC1="$SRC1 lamspectraltransform_mod.ftn90 gridstatevector_mod.ftn90 ensemblestatevector_mod.ftn90"
  SRC1="$SRC1 columndata_mod.ftn90 spectralfilter_mod.ftn90 lambmatrixhi_mod.ftn90"
  SRC1="$SRC1 variabletransforms_mod.ftn90 bmatrixensemble_mod.ftn90 bmatrixhi_mod.ftn90 bmatrixchem_mod.ftn90"
  SRC1="$SRC1 bmatrix_mod.ftn90"
  s.compile $COMPF -O ${FOPTMIZ} -src $SRC1 > listing1 2>&1
  status=1
  grep fail listing1 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "... !! Compilation aborted: check in ${PWD}/listing1 !!"
      exit 1
  fi

  echo "... > Compiling remaining ftn90..."
  echo "...   if aborting, check in ${PWD}/listing2"
  filelist=""
  for i in *.ftn90 ; do
      if [[ "${i}" != *_mod.ftn* ]]; then
	  filelist="$filelist ${i}"
      fi
  done
  s.compile $COMPF -O ${FOPTMIZ} -src $filelist > listing2 2>&1
  status=1
  grep fail listing2 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "... !! Compilation aborted: check in ${PWD}/listing2 !!"
      exit 1
  fi

  echo "... > Building the executable ${varabs}"
  rm -f ${varabs}
  echo "...   if aborting, check in ${PWD}/listing3"
  s.compile $COMPF  -O ${FOPTMIZ} ${MPILIBDIR} -libappl $LIBAPPL -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o ${varabs} > listing3 2>&1
  status=1
  grep fail listing3 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "... !! Compilation aborted: check in ${PWD}/listing3 !!"
      exit 1
  fi

  status=1
  grep -i " ERROR " listing? || status=0
  if [ "${status}" -ne 0 ] ; then
      echo "... !! ERROR found: STOP; check listing in ${PWD} !!"
      exit 1
  fi
  cp ${varabs} ${absdir}/

elif [ $mode == abs ] ; then

  rm -f ${varabs}

  echo "..."
  echo "... > Building the executable ${varabs}"
  echo "..."
  s.compile $COMPF  -O ${FOPTMIZ} ${MPILIBDIR} -libappl $LIBAPPL -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o ${varabs}
  cp ${varabs} ${absdir}/

else

  if [ -f $randompertdir/$mode ] ; then
    file=`basename $mode`
    rm -f $file
    cp $randompertdir/$mode .
    echo "..."
    echo "... compiling $mode"
    echo "..."
    s.compile $COMPF  -O ${FOPTMIZ} -src $file
  else
    echo "..."
    echo "... !! File $randompertdir/$mode does NOT exist. Stop !!"
    exit 1
  fi

fi

echo "..."
echo "... > FINISHED COMPILATION AT: $(date)"
if [ "${mode}" == full -o "${mode}" == abs ] ; then
    echo "..."
    echo "... The program can be found here: ${absdir}/${varabs}"
    echo "..."
else
    echo "..."
fi

echo "...             |======================================|"
echo "... ------------| RANDOMPERT compilation script ENDING |------------"
echo "...             |======================================|"
