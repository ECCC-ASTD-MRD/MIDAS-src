#!/bin/bash

set -e

mode=$1
nompi=$2

. ./compile_commons.sh

echo "..."
echo "...             |======================================|"
echo "... ------------|  OMINUSF compilation script STARTING |------------"
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
compiledir_main=${COMPILEDIR_OAVAR_MAIN:-"../../compiledir"}
compiledir_ID=${COMPILEDIR_OAVAR_ID:-$revnum}
compiledir=${compiledir_main}/compiledir-${ORDENV_PLAT}_${compiledir_ID}
mkdir -p $compiledir
cd $compiledir
compiledir=${PWD} # needed when compiledir_main = ".."

echo "..."
echo "... > Compiledir set to $compiledir"

if [ ${compiledir_main} != "../../compiledir" ] ; then
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

ominusfabs=ominusf_${ORDENV_PLAT}${ABSTAG}-${revnum}.Abs

#-----------------------------------------------------------------------------

LIBAPPL="rttov_coef_io rttov_hdf rttov_parallel  rttov_main rttov_emis_atlas rttov_other ${HDF5_LIBS} burp_module descrip $MPILIB"
LIBSYS="hpcoperf"
LIBRMN=rmnMP

if [ "${mode}" == full ] ; then

  rm -f *.o *.mod *.cdk* *.h *.ftn* *.f *.f90

  cp -f ${trunkdir}/ominusf/main_ominusf.f90 ${compiledir}/
  cd ${trunkdir};          ls -1F | grep -v '/' | grep -v "*" | grep -v "@" | cpio --quiet -p $compiledir ; cd $compiledir
  cd ${trunkdir}/bgcheck;  ls -1F | grep -v '/' | grep -v "*" | cpio --quiet -p $compiledir ; cd $compiledir
  cd ${trunkdir}/shared;   ls -1F | grep -v '/' | grep -v "*" | cpio --quiet -p $compiledir ; cd $compiledir
  rm -f *.f*~

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
  SRC0="mathphysconstants_mod.f90 earthconstants_mod.f90 utilities_mod.f90 ramdisk_mod.ftn90"
  SRC0="$SRC0 randomnumber_mod.f90 mpi_mod.f90 mpivar_mod.f90 bufr_mod.f90 codtyp_mod.f90"
  SRC0="$SRC0 physicsfunctions_mod.f90 obsspacedata_mod.ftn90 localizationfunction_mod.f90"
  SRC0="$SRC0 horizontalcoord_mod.f90 timecoord_mod.f90 verticalcoord_mod.f90"
  SRC0="$SRC0 lqtoes_mod.f90 presprofileoperators_mod.f90 spectralfilter_mod.f90 tovs_extrap_mod.f90"
  s.compile $COMPF -O ${FOPTMIZ} -src $SRC0 > listing0 2>&1
  status=1
  grep fail listing0 || status=0
  if [ "${status}" -ne 0 ] ; then
      echo "... !! Compilation aborted: check in ${PWD}/listing0 !!"
      exit 1
  fi

  echo "... > Compiling quasi-newton module"
  echo "...   if aborting, check in ${PWD}/listing1"
  SRC0="quasinewton_mod.f"
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
  SRC0="windrotation_mod.f90 analysisgrid_mod.f90"
  s.compile $COMPF  -O ${FOPTMIZ} -src $SRC0 > listing2 2>&1
  status=1
  grep fail listing2 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "... !! Compilation aborted: check in ${PWD}/listing2 !!"
      exit 1
  fi

  echo "... > Compiling most of the new modules"
  echo "...   if aborting, check in ${PWD}/listing3"
  SRC1="varnamelist_mod.f90 columndata_mod.f90 controlvector_mod.f90 rmatrix_mod.ftn90 hirchannels_mod.f90 ozoneclim_mod.f90 tovs_nl_mod.ftn90"
  SRC1="$SRC1 globalspectraltransform_mod.f90 tt2phi_mod.f90"
  SRC1="$SRC1 lamspectraltransform_mod.f90 gridstatevector_mod.f90 ensemblestatevector_mod.f90 statetocolumn_mod.f90"
  SRC1="$SRC1 variabletransforms_mod.f90 localizationspectral_mod.f90 localization_mod.f90 diffusion_mod.f90"
  SRC1="$SRC1 bmatrixensemble_mod.f90 bmatrixhi_mod.f90 bmatrixlatbands_mod.f90 lambmatrixhi_mod.f90"
  SRC1="$SRC1 bmatrixchem_mod.f90 bmatrixdiff_mod.f90 bmatrix_mod.f90 residual_mod.f90 costfunction_mod.f90"

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
  SRC2="gps_mod.f90"
  s.compile $COMPF  -O ${FOPTMIZ} -src $SRC2 > listing5 2>&1
  status=1
  grep fail listing5 || status=0
  if [ "${status}" -ne 0 ] ; then
      echo "... !! Compilation aborted: check in ${PWD}/listing5 !!"
      exit 1
  fi

  echo "... > Compiling some more modules..."
  echo "...   if aborting, check in ${PWD}/listing6"
  SRC2="obssubspacedata_mod.ftn90 burpfiles_mod.ftn90 multi_ir_bgck_mod.ftn90 chem_setup_mod.f90 chem_obserrors_mod.f90"
  SRC2="$SRC2 chem_obsoperators_mod.f90 chem_postproc_mod.f90"
  SRC2="$SRC2 obserrors_mod.f90 varqc_mod.f90 obsfilter_mod.f90 tovs_lin_mod.ftn90 obsoperators_mod.f90 obsspacediag_mod.f90"
  SRC2="$SRC2 innovation_mod.f90 minimization_mod.f90"
  s.compile $COMPF  -O ${FOPTMIZ} -src $SRC2 > listing6 2>&1
  status=1
  grep fail listing6 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "... !! Compilation aborted: check in ${PWD}/listing6 !!"
      exit 1
  fi

  echo "... > Compiling remaining .f*90..."
  echo "...   if aborting, check in ${PWD}/listing7"
  filelist=""
  for i in *.*90 ; do
      if [[ "${i}" != *_mod.f* ]] && [[ "${i}" != main_var.f90* ]] ; then
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

  echo "... > Building the executable ${ominusfabs}"
  rm -f ${ominusfabs}
  echo "...   if aborting, check in ${PWD}/listing8"
  s.compile $COMPF  -O ${FOPTMIZ} ${MPILIBDIR} -libappl $LIBAPPL -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o ${ominusfabs} > listing8 2>&1
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
  cp ${ominusfabs} ${absdir}/

  #### rm -f *.ftn* *.f *.f90

elif [ "${mode}" == abs ] ; then

  rm -f ${ominusfabs}

  echo "..."
  echo "... > Building the executable ${ominusfabs}"
  echo "..."
  s.compile $COMPF  -O ${FOPTMIZ} ${MPILIBDIR} -libappl $LIBAPPL -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o ${ominusfabs}
  cp ${ominusfabs} ${absdir}/

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
    echo "... The program can be found here: ${absdir}/${ominusfabs}"
    echo "..."
else
    echo "..."
fi

echo "...             |======================================|"
echo "... ------------|  OMINUSF compilation script ENDING   |------------"
echo "...             |======================================|"
