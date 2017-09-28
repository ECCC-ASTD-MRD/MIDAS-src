#!/bin/bash

set -e

mode=$1
program="ensmanip"

. ./compile_commons.sh

echo "..."
echo "...             |=====================================================|"
echo "... ------------| Compilation script STARTING for program: ${program} |------------"
echo "...             |=====================================================|"

if [ "$mode" == "" ] ; then
  echo "..."
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

codesubdir=$PWD
trunkdir=$PWD/../

# Get revision number
revnum=$(git describe --abbrev=7 --always --dirty=_M 2>/dev/null || ssh eccc-ppp1 "cd $trunkdir; git describe --abbrev=7 --always --dirty=_M" 2>/dev/null || echo unkown revision)
echo "..."
echo "... > Revision Number = '$revnum'"

# Set compiledir
compiledir_main=${COMPILEDIR_OAVAR_MAIN:-"../../compiledir"}
compiledir_ID=${COMPILEDIR_OAVAR_ID:-$revnum}
compiledir=${compiledir_main}/compiledir-${program}-${ORDENV_PLAT}_${compiledir_ID}
mkdir -p $compiledir
cd $compiledir
compiledir=${PWD} # needed when compiledir_main = ".."

echo "..."
echo "... > Compiledir set to ${compiledir}"

if [ ${compiledir_main} != "../../compiledir" ] ; then
    if [ ! -d  ${trunkdir}/../compiledir/compiledir-${program}-${ORDENV_PLAT}_${compiledir_ID} ] ; then
	if [ ! -d  ${trunkdir}/../compiledir ] ; then
	    mkdir -p ${trunkdir}/../compiledir
	fi
	ln -s ${compiledir_main}/compiledir-${program}-${ORDENV_PLAT}_${compiledir_ID} ${trunkdir}/../compiledir/compiledir-${program}-${ORDENV_PLAT}_${compiledir_ID}
    fi
fi

absdir=${compiledir_main}/oavar_abs
mkdir -p ${absdir}
cd $absdir ; absdir=$PWD ; cd - >/dev/null

varabs=${program}_${ORDENV_PLAT}-${revnum}.Abs

#-----------------------------------------------------------------------------

# LIBAPPL defined in "src_files" script
LIBSYS="hpcoperf"
LIBRMN=rmnMP
. ${codesubdir}/src_files_${program}.sh

if [ $mode == full ] ; then

  rm -f *.o *.mod *.cdk* *.h *.ftn* *.f *.f90

  # Create a local copy of the source code
  cp -f ${trunkdir}/${program}/main_${program}.f90 ${compiledir}/
  cp -f ${trunkdir}/*.f*90 ${compiledir}/
  cp -f ${trunkdir}/*.inc ${compiledir}/
  cp -f ${trunkdir}/shared/*.f*90 ${compiledir}/
  cp -f ${trunkdir}/shared/*.inc ${compiledir}/
  cp -f ${trunkdir}/bgcheck/*.f*90 ${compiledir}/

  cd ${compiledir}

  # Add revision number to the main routine
  sed -i "s|GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE|${revnum}|g" main_${program}.f90

  echo "..."
  echo "... > STARTING COMPILATION AT: $(date)"
  echo "..."

  echo "... > Compiling all modules and subroutines ..."
  echo "...   if aborting, check in ${PWD}/listing"

  s.compile $COMPF -O ${FOPTMIZ} -src $SRC_FILES > listing 2>&1
  status=1
  grep fail listing || status=0
  if [ "${status}" -ne 0 ] ; then
      echo "... !! Compilation aborted: check in ${PWD}/listing !!"
      exit 1
  fi

  echo "... > Compiling main program..."
  echo "...   if aborting, check in ${PWD}/listing_main"
  s.compile $COMPF -O ${FOPTMIZ} -src main_${program}.f90 > listing_main 2>&1
  status=1
  grep fail listing_main || status=0
  if [ "${status}" -ne 0 ]; then
      echo "... !! Compilation aborted: check in ${PWD}/listing_main !!"
      exit 1
  fi

  echo "... > Building the executable ${varabs}"
  rm -f ${varabs}
  echo "...   if aborting, check in ${PWD}/listing_abs"
  s.compile $COMPF  -O ${FOPTMIZ} ${MPILIBDIR} -libappl $LIBAPPL -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o ${varabs} > listing_abs 2>&1
  status=1
  grep fail listing_abs || status=0
  if [ "${status}" -ne 0 ]; then
      echo "... !! Compilation aborted: check in ${PWD}/listing_abs !!"
      exit 1
  fi

  status=1
  grep -i " ERROR " listing* || status=0
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

  if [ -f $codesubdir/$mode ] ; then
    file=`basename $mode`
    rm -f $file
    cp $codesubdir/$mode .
    echo "..."
    echo "... compiling $mode"
    echo "..."
    s.compile $COMPF  -O ${FOPTMIZ} -src $file
  else
    echo "..."
    echo "... !! File $codesubdir/$mode does NOT exist. Stop !!"
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

echo "...             |==================================================|"
echo "... ------------| Compilation script ENDING for program ${program} |------------"
echo "...             |==================================================|"
echo ""
