#!/bin/bash
set -e

# set the resources.def file, which depends on the TRUE_HOST name
../../set_resources_def.sh

if [ $# -lt 1 ]; then
    echo "compile_program.sh: You must give at least one argument which is the program the compile"
    exit 1
fi
program=$(basename $1 .f90)
mode=$2
deleteCompileDir=$3

export COMPILE_MIDAS_ADD_DEBUG_OPTIONS=${COMPILE_MIDAS_ADD_DEBUG_OPTIONS:-no}
. ./commons/compile_setup.sh

echo "..."
echo "...             |=====================================================|"
echo "... ------------| Compilation script STARTING for program: ${program}  "
echo "...             |=====================================================|"

if [ "$mode" == "" ] ; then
  echo "..."
  echo "... !WARNING! no compilation mode specified, assuming 'full'"
  mode=full
fi

if [ "${mode}" == full ] ; then
  echo "..."
  echo "... > Full Compilation"
elif [ "${mode}" == abs ] ; then
  echo "..."
  echo "... > Building the Executable only"
else
  echo "..."
  echo "... > Compiling only routine : $mode"
fi

if [ "$deleteCompileDir" == "" ] ; then
  echo "..."
  echo "... !WARNING! no directive given for the deletion of the compilation directory, assuming 'yes'"
  deleteCompileDir=yes
fi

if [ "${deleteCompileDir}" == yes ] ; then
  echo "..."
  echo "... > the compilation directory will be DELETED after the successful completion of this task"
else
  echo "..."
  echo "... > the compilation directory will be KEPT"
fi


echo "..."
echo "... > Compiling for an MPI executable"
MPILIB="rpn_comm"

programsDir=$PWD
modulesDir=$PWD/../modules
depotDir=$PWD/..

# Get revision number
toplevel=$(git rev-parse --show-toplevel)
revnum=$(${toplevel}/midas.version.sh)
echo "..."
echo "... > Revision Number = '$revnum'"

# Set compiledir
compiledir_main=${COMPILEDIR_MIDAS_MAIN:-"../../compiledir"}
compiledir_ID=${COMPILEDIR_MIDAS_ID:-$revnum}
compiledir=${compiledir_main}/compiledir-${program}-${ORDENV_PLAT}_${compiledir_ID}
mkdir -p $compiledir
cd $compiledir
compiledir=${PWD} # needed when compiledir_main = ".."

echo "..."
echo "... > Compiledir set to ${compiledir}"

if [ ${compiledir_main} != "../../compiledir" ] ; then
    if [ ! -d  ${depotDir}/../compiledir/compiledir-${program}-${ORDENV_PLAT}_${compiledir_ID} ] ; then
	if [ ! -d  ${depotDir}/../compiledir ] ; then
	    mkdir -p ${depotDir}/../compiledir
	fi
	ln -s ${compiledir_main}/compiledir-${program}-${ORDENV_PLAT}_${compiledir_ID} ${depotDir}/../compiledir/compiledir-${program}-${ORDENV_PLAT}_${compiledir_ID}
    fi
fi

absdir=${compiledir_main}/midas_abs
mkdir -p ${absdir}
cd $absdir ; absdir=$PWD ; cd - >/dev/null

midasAbs=midas-${program}_${ORDENV_PLAT}-${revnum}.Abs

#-----------------------------------------------------------------------------

# LIBAPPL defined in "src_files" script
LIBSYS="hpcoperf sqlite3"
LIBRMN=rmnMP
. ${programsDir}/src_files/compile_setup_${program}.sh
. ${programsDir}/src_files/src_files_${program}.sh

LINK_LIBS=
for thislib in ${LIBAPPL} ${LIBSYS} ${LIBRMN} ${LIBIRC}; do
    LINK_LIBS="${LINK_LIBS} -l${thislib}"
done

if [ $mode == full ] ; then

  rm -f *.o *.mod *.cdk* *.h *.ftn* *.f *.f90

  # Create a local copy of the source code
  cp $programsDir/${program}.f90 ${compiledir}/
  if [ -d ${modulesDir}/${program} ] ; then
      cp -f ${modulesDir}/${program}/*.f* ${compiledir}/
  fi
  cp -f ${modulesDir}/*.f* ${compiledir}/

  cd ${compiledir}

  # Add revision number to the main routine
  sed -i "s|GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE|${revnum}|g" ${program}.f90

  echo "..."
  echo "... > STARTING COMPILATION AT: $(date)"
  echo "..."

  echo "... > Preprocessing *.ftn90 ..."
  echo "...   if aborting, check in ${PWD}/listing_preproc"

  SRC_FILES_F90=
  defines=
  for file in ${SRC_FILES}; do
      if [[ "${file}" = *.ftn90 ]]; then
          if [ -z "${defines}" ]; then
              for opt in ${COMPF}; do
                  if [[ "${opt}" = -D* ]]; then
                      defines="${defines} ${opt}"
                  fi
              done
          fi
          file90=$(basename ${file} .ftn90).f90
          echo ${file} >> listing_preproc
          r.gppf -lang-f90+ ${GPP_OPTS} ${defines} ${file} > ${file90} 2>> listing_preproc
      else
          file90=${file}
      fi
      SRC_FILES_F90="${SRC_FILES_F90} ${file90}"
  done

  echo "... > Compiling all modules and subroutines ..."
  echo "...   if aborting, check in ${PWD}/listing"

  ## This variable makes 's.f90' more verbose
  export Verbose=yes

  s.f90 ${COMPF} -O ${FOPTMIZ} -c ${SRC_FILES_F90} > listing 2>&1
  status=1
  grep fail listing || status=0
  if [ "${status}" -ne 0 ] ; then
      echo "... !! Compilation aborted: check in ${PWD}/listing !!"
      exit 1
  fi

  echo "... > Compiling main program..."
  echo "...   if aborting, check in ${PWD}/listing_main"
  s.f90 ${COMPF} -O ${FOPTMIZ} -c ${program}.f90 > listing_main 2>&1
  status=1
  grep fail listing_main || status=0
  if [ "${status}" -ne 0 ]; then
      echo "... !! Compilation aborted: check in ${PWD}/listing_main !!"
      exit 1
  fi

  echo "... > Building the executable ${midasAbs}"
  rm -f ${midasAbs}
  echo "...   if aborting, check in ${PWD}/listing_abs"
  s.f90 ${COMPF} -O ${FOPTMIZ} ${LINK_LIBS} *.o -o ${midasAbs} > listing_abs 2>&1
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
  cp ${midasAbs} ${absdir}/

  if [ "${deleteCompileDir}" == yes ] ; then
      rm -rf ${compiledir}
  fi

elif [ $mode == abs ] ; then

  rm -f ${midasAbs}

  echo "..."
  echo "... > Building the executable ${midasAbs}"
  echo "..."
  s.f90 ${COMPF} -O ${FOPTMIZ} ${LINK_LIBS} *.o -o ${midasAbs}
  cp ${midasAbs} ${absdir}/

else

  if [ -f $modulesDir/$mode ] ; then
    file=`basename $mode`
    rm -f $file
    cp ${modulesDir}/$mode .
    echo "..."
    echo "... compiling $mode"
    echo "..."
    s.f90 ${COMPF} -O ${FOPTMIZ} -c ${file}
  elif [ -f $programsDir/$mode ] ; then
    file=`basename $mode`
    rm -f $file
    cp ${programsDir}/$mode .
    echo "..."
    echo "... compiling $mode"
    echo "..."
    s.f90 ${COMPF} -O ${FOPTMIZ} -c ${file}
  else
    echo "..."
    echo "... !! File $mode does NOT exist in directories ${modulesDir} and ${programsDir}. Stop !!"
    exit 1
  fi

fi

echo "..."
echo "... > FINISHED COMPILATION AT: $(date)"
if [ "${mode}" == full -o "${mode}" == abs ] ; then
    echo "..."
    echo "... The program can be found here: ${absdir}/${midasAbs}"
    echo "..."
else
    echo "..."
fi

echo "...             |==================================================|"
echo "... ------------| Compilation script ENDING for program ${program}  "
echo "...             |==================================================|"
echo ""
