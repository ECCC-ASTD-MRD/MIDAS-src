#!/bin/bash

set -e

mode=$1
nompi=$2

. ./compile_commons.sh

echo "..."
echo "...             |======================================|"
echo "... ------------|  OMINUSP compilation script STARTING |------------"
echo "...             |======================================|"

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

if [ "$nompi" = "NOMPI" -o "$nompi" = "nompi" ] ; then
  echo "..."
  echo "... > Compiling for a NON-MPI executable"
  MPILIB="rpn_commstubs rpn_comm"
  ABSTAG="_nompi"
else
  echo "..."
  echo "... > Compiling for an MPI executable"
  MPILIB="rpn_comm"
  ABSTAG=""
fi

trunkdir=$PWD/../

# Get revision number
revnum=$(git describe --abbrev=7 --always --dirty=_M 2>/dev/null || ssh eccc-ppp1 "cd $trunkdir; git describe --abbrev=7 --always --dirty=_M" 2>/dev/null || echo unkown revision)
echo "..."
echo "... > Revision Number = '$revnum'"

# Set compiledir
compiledir_main=${COMPILEDIR_OAVAR_MAIN:-"../compiledir"}
compiledir_ID=${COMPILEDIR_OAVAR_ID:-$revnum}
compiledir=${compiledir_main}/compiledir-${ORDENV_PLAT}_${compiledir_ID}
mkdir -p $compiledir
cd $compiledir
compiledir=${PWD} # needed when compiledir_main = ".."

echo "..."
echo "... > Compiledir set to $compiledir"

if [ ${compiledir_main} != "../compiledir" ] ; then
    if [ ! -d  ${trunkdir}/../compiledir/compiledir-${ORDENV_PLAT}_${compiledir_ID} ] ; then
	if [ ! -d  ${trunkdir}/../compiledir ] ; then
	    mkdir -p ${trunkdir}/../compiledir
	fi
	ln -s ${compiledir_main}/compiledir-${ORDENV_PLAT}_${compiledir_ID} ${trunkdir}/../compiledir/compiledir-${ORDENV_PLAT}_${compiledir_ID}
    fi
fi

absdir=${compiledir_main}/oavar_abs
mkdir -p ${absdir}
cd $absdir ; absdir=$PWD ; cd - >/dev/null

ominuspabs=ominusp_${ORDENV_PLAT}${ABSTAG}-${revnum}.Abs

#-----------------------------------------------------------------------------

LIBAPPL="netcdff rttov10.2.0_coef_io rttov10.2.0_main rttov10.2.0_emis_atlas rttov10.2.0_other burp_module descrip $MPILIB"
LIBSYS="hpcoperf"
LIBRMN=rmnMP

if [ "${mode}" == full ] ; then

  rm -f *.o *.mod *.cdk* *.h *.ftn* *.f *.f90

  cp -f ${trunkdir}/ominusp/main_ominusp.ftn90 ${compiledir}/
  cd ${trunkdir};          ls -1F | grep -v '/' | grep -v "*" | grep -v "@" | cpio --quiet -p $compiledir ; cd $compiledir
  cd ${trunkdir}/bgcheck;  ls -1F | grep -v '/' | grep -v "*" | cpio --quiet -p $compiledir ; cd $compiledir
  cd ${trunkdir}/shared;   ls -1F | grep -v '/' | grep -v "*" | cpio --quiet -p $compiledir ; cd $compiledir
  rm -f *.ftn~ *.ftn90~

  # Check for indented OPEN-MP directives - this is not allowed!
  status=1
  grep -i ' !$omp' *.ftn* || status=0
  if [ "${status}" -ne 0 ] ; then
      echo "... !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      echo "... !! Compilation aborted: check the code for indented OPEN-MP directives !!"
      echo "... !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      exit 1
  fi

  # Compile the subroutines...
  echo "... > Compiling low-level independent modules"
  echo "...   if aborting, check in ${PWD}/listing0"
  SRC0="mathphysconstants_mod.ftn90 earthconstants_mod.ftn90 utilities_mod.ftn90 ramdisk_mod.ftn90"
  SRC0="$SRC0 randomnumber_mod.ftn90 mpi_mod.ftn90 mpivar_mod.ftn90 bufr_mod.ftn90 codtyp_mod.ftn90"
  SRC0="$SRC0 physicsfunctions_mod.ftn90 obsspacedata_mod.ftn90 localizationfunction_mod.ftn90"
  SRC0="$SRC0 horizontalcoord_mod.ftn90 timecoord_mod.ftn90 verticalcoord_mod.ftn90"
  SRC0="$SRC0 lqtoes_mod.ftn90 presprofileoperators_mod.ftn90 spectralfilter_mod.ftn90"
  s.compile $COMPF -O ${FOPTMIZ} -src $SRC0 > listing0 2>&1
  status=1
  grep fail listing0 || status=0
  if [ "${status}" -ne 0 ] ; then
      echo "... !! Compilation aborted: check in ${PWD}/listing0 !!"
      exit 1
  fi

  echo "... > Compiling quasi-newton module"
  echo "...   if aborting, check in ${PWD}/listing1"
  SRC0="quasinewton_mod.ftn"
  s.compile $COMPF_NOC  -O ${FOPTMIZ} -src $SRC0 > listing1 2>&1
  status=1
  grep fail listing1 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "... !! Compilation aborted: check in ${PWD}/listing1 !!"
      exit 1
  fi
  rm -f $SRC0

  echo "... > Compiling analysis grid modules"
  echo "...   if aborting, check in ${PWD}/listing2"
  SRC0="windrotation_mod.ftn90 analysisgrid_mod.ftn90"
  s.compile $COMPF  -O ${FOPTMIZ} -src $SRC0 > listing2 2>&1
  status=1
  grep fail listing2 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "... !! Compilation aborted: check in ${PWD}/listing2 !!"
      exit 1
  fi

  echo "... > Compiling most of the new modules"
  echo "...   if aborting, check in ${PWD}/listing3"
  SRC1="controlvector_mod.ftn90 rmatrix_mod.ftn90 hirchannels_mod.ftn90 tovs_nl_mod.ftn90"
  SRC1="$SRC1 tovs_lin_mod.ftn90 varnamelist_mod.ftn90 columndata_mod.ftn90 multi_ir_bgck_mod.ftn90"
  SRC1="$SRC1 emissivities_mod.ftn90 globalspectraltransform_mod.ftn90 tt2phi_mod.ftn90"
  SRC1="$SRC1 lamspectraltransform_mod.ftn90 gridstatevector_mod.ftn90 ensemblestatevector_mod.ftn90 statetocolumn_mod.ftn90"
  SRC1="$SRC1 variabletransforms_mod.ftn90 localizationspectral_mod.ftn90 localization_mod.ftn90"
  SRC1="$SRC1 bmatrixensemble_mod.ftn90 bmatrixhi_mod.ftn90 lambmatrixhi_mod.ftn90"
  SRC1="$SRC1 bmatrixchem_mod.ftn90 bmatrix_mod.ftn90 residual_mod.ftn90 costfunction_mod.ftn90"
  SRC1="$SRC1 ozoneclim_mod.ftn90 tovs_extrap_mod.ftn90"

  s.compile $COMPF  -O ${FOPTMIZ} -src $SRC1 > listing3 2>&1
  status=1
  grep fail listing3 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "... !! Compilation aborted: check in ${PWD}/listing3 !!"
      exit 1
  fi

  echo "... > Compiling burpread module"
  echo "...   if aborting, check in ${PWD}/listing4"
  SRC1="burpread_mod.ftn90"
  s.compile $COMPF  -O ${FOPTMIZ} -src $SRC1 > listing4 2>&1
  status=1
  grep fail listing4 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "... !! Compilation aborted: check in ${PWD}/listing4 !!"
      exit 1
  fi
  
  echo "... > Compiling the GPS module ..."
  echo "...   if aborting, check in ${PWD}/listing5"
  SRC2="gps_mod.ftn90"
  s.compile $COMPF  -O ${FOPTMIZ} -src $SRC2 > listing5 2>&1
  status=1
  grep fail listing5 || status=0
  if [ "${status}" -ne 0 ] ; then
      echo "... !! Compilation aborted: check in ${PWD}/listing5 !!"
      exit 1
  fi

  echo "... > Compiling some more modules..."
  echo "...   if aborting, check in ${PWD}/listing6"
  SRC2="obssubspacedata_mod.ftn90 burpfiles_mod.ftn90 chem_setup_mod.ftn90 chem_obserrors_mod.ftn90"
  SRC2="$SRC2 chem_obsoperators_mod.ftn90 chem_postproc_mod.ftn90"
  SRC2="$SRC2 obserrors_mod.ftn90 varqc_mod.ftn90 obsfilter_mod.ftn90 obsoperators_mod.ftn90 obsspacediag_mod.ftn90"
  SRC2="$SRC2 innovation_mod.ftn90 minimization_mod.ftn90"
  s.compile $COMPF  -O ${FOPTMIZ} -src $SRC2 > listing6 2>&1
  status=1
  grep fail listing6 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "... !! Compilation aborted: check in ${PWD}/listing6 !!"
      exit 1
  fi

  echo "... > Compiling remaining ftn ftn90..."
  echo "...   if aborting, check in ${PWD}/listing7"
  filelist=""
  for i in *.ftn90 ; do
      if [[ "${i}" != *_mod.ftn* ]] && [[ "${i}" != main_var.ftn90 ]]; then
	  filelist="$filelist ${i}"
      fi
  done
  s.compile $COMPF  -O ${FOPTMIZ} -src $filelist > listing7 2>&1
  status=1
  grep fail listing7 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "... !! Compilation aborted: check in ${PWD}/listing7 !!"
      exit 1
  fi

  echo "... > Building the executable ${ominuspabs}"
  rm -f ${ominuspabs}
  echo "...   if aborting, check in ${PWD}/listing8"
  s.compile $COMPF  -O ${FOPTMIZ} ${MPILIBDIR} -libappl $LIBAPPL -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o ${ominuspabs} > listing8 2>&1
  status=1
  grep fail listing8 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "... !! Compilation aborted: check in ${PWD}/listing8 !!"
      exit 1
  fi

  status=1
  grep -i " ERROR " listing? || status=0
  if [ "${status}" -ne 0 ] ; then
      echo "... !! ERROR found: STOP; check listing in ${PWD} !!"
      exit 1
  fi
  cp ${ominuspabs} ${absdir}/

  #### rm -f *.ftn* *.f *.f90

elif [ "${mode}" == abs ] ; then

  rm -f ${ominuspabs}

  echo "..."
  echo "... > Building the executable ${ominuspabs}"
  echo "..."
  s.compile $COMPF  -O ${FOPTMIZ} ${MPILIBDIR} -libappl $LIBAPPL -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o ${ominuspabs}
  cp ${ominuspabs} ${absdir}/

else
    if [ -f $trunkdir/$mode ] ; then
	file=`basename $mode`
	rm -f $file
	cp $trunkdir/$mode .
	echo "..."
	echo "... compiling $mode"
	echo "..."
	s.compile $COMPF  -O ${FOPTMIZ} -src $file
    else
	echo "..."
	echo "... !! File $trunkdir/$mode does NOT exist. Stop !!"
	exit 1
    fi
fi

echo "..."
echo "... > FINISHED COMPILATION AT: $(date)"
if [ "${mode}" == full -o "${mode}" == abs ] ; then
    echo "..."
    echo "... The program can be found here: ${absdir}/${ominuspabs}"
    echo "..."
else
    echo "..."
fi

echo "...             |======================================|"
echo "... ------------|  OMINUSP compilation script ENDING   |------------"
echo "...             |======================================|"