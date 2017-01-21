#!/bin/bash

set -e

mode=$1
nompi=$2

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

if [ "$nompi" = "NOMPI" -o "$nompi" = "nompi" ] 
then
  echo "!!Compiling for a NON-MPI executable!!"
  MPILIB="rpn_commstubs rpn_comm"
  MPIKEY=""
  ABSTAG="_nompi"
else
  echo "!!Compiling for an MPI executable!!"
  MPILIB="rpn_comm"
  MPIKEY="-mpi"
  ABSTAG=""
fi

trunkdir=$PWD

if [ "${ORDENV_PLAT}" = sles-11-haswell-64-xc40 ];then
    echo "Switching ORDENV_PLAT from '${ORDENV_PLAT}' to 'sles-11-broadwell-64-xc40'"
    . r.env.dot --arch sles-11-broadwell-64-xc40
    echo ORDENV_PLAT=${ORDENV_PLAT}
fi

if [ "${ORDENV_PLAT}" = ubuntu-14.04-amd64-64 ];then
    FOPTMIZ=2
elif [ "${ORDENV_PLAT}" = sles-11-amd64-64 ];then
    FOPTMIZ=4
elif [ "${ORDENV_PLAT}" = sles-11-broadwell-64-xc40 ];then
    FOPTMIZ=4
else
    echo "This platform 'ORDENV_PLAT=${ORDENV_PLAT}' is not supported.  Only 'ubuntu-14.04-amd64-64' and 'ubuntu-14.04-amd64-64' are."
    exit 1
fi

# automatically set the global revision number in toplevelcontrol_mod.ftn90 by
# replacing the string XXXXX with the actual revision number
revnum=$(git describe --always --dirty 2>/dev/null || ssh eccc-ppp1 "cd $trunkdir; git describe --always --dirty" 2>/dev/null || echo unkown revision)
echo " "
echo "-----------------------"
echo "Revision number='$revnum'"
echo "-----------------------"
echo " "

cd ../
mkdir -p compiledir-${ORDENV_PLAT}
cd compiledir-${ORDENV_PLAT}
#rm -f *.o *.f *.f90 *.mod
compiledir=${PWD}

#----------------------------------------------------------------
#  Set up dependent librarys and tools. 
#---------------------------------------------------------------
## for s.compile
echo "loading hpco/tmp/eccc/201402/06/base"
. ssmuse-sh -d hpco/tmp/eccc/201402/06/base
## for the compiler
if [ "${ORDENV_PLAT}" = ubuntu-14.04-amd64-64 ];then
    echo "loading compiler main/opt/intelcomp/intelcomp-2016.1.156"
    . ssmuse-sh -d main/opt/intelcomp/intelcomp-2016.1.156
elif [ "${ORDENV_PLAT}" = sles-11-amd64-64 -o "${ORDENV_PLAT}" = sles-11-broadwell-64-xc40 ];then
    echo "loading compiler PrgEnv-intel-5.2.82"
    module load PrgEnv-intel/5.2.82
else
    echo "This platform 'ORDENV_PLAT=${ORDENV_PLAT}' is not supported.  Only 'ubuntu-14.04-amd64-64' and 'sles-11-amd64-64' are."
    exit 1
fi

varabs=oavar_${ORDENV_PLAT}${ABSTAG}

## for rmn, rpncomm
echo "loading eccc/mrd/rpn/libs/16.1"
. ssmuse-sh -d eccc/mrd/rpn/libs/16.1
if [ "${ORDENV_PLAT}" = ubuntu-14.04-amd64-64 ];then
    ## for openmpi
    echo "loading main/opt/openmpi/openmpi-1.6.5/intelcomp-2016.1.156"
    . ssmuse-sh -d main/opt/openmpi/openmpi-1.6.5/intelcomp-2016.1.156
elif [ "${ORDENV_PLAT}" = sles-11-amd64-64 -o "${ORDENV_PLAT}" = sles-11-broadwell-64-xc40 ];then
    true
fi
## for 'vgrid'
echo "loading eccc/cmd/cmdn/vgrid/5.6.8/${COMP_ARCH}"
. ssmuse-sh -d eccc/cmd/cmdn/vgrid/5.6.8/${COMP_ARCH}
## for 'burplib'
echo "loading eccc/cmd/cmda/libs/16.1-3/${COMP_ARCH}"
. ssmuse-sh -d eccc/cmd/cmda/libs/16.1-3/${COMP_ARCH}
## For hpcoperf needed for TMG timings
echo "loading main/opt/perftools/perftools-2.0/${COMP_ARCH}"
. ssmuse-sh -d main/opt/perftools/perftools-2.0/${COMP_ARCH}
# For RTTOV 10v3 package...
echo "loading eccc/mrd/rpn/anl/rttov/10v3.2/${COMP_ARCH}"
. ssmuse-sh -d eccc/mrd/rpn/anl/rttov/10v3.2/${COMP_ARCH}

#-----------------------------------------------------------------------------

LIBAPPL="rttov10.2.0_coef_io rttov10.2.0_main rttov10.2.0_other burp_module descrip $MPILIB"
LIBSYS="hpcoperf"

LIBRMN=rmnMP

COMPF_GLOBAL="-openmp ${MPIKEY}"
OPTF="=-check =noarg_temp_created"
if [ "${ORDENV_PLAT}" = ubuntu-14.04-amd64-64 ];then
    OPTF="=-mkl ${OPTF}"
elif [ "${ORDENV_PLAT}" = sles-11-amd64-64 -o "${ORDENV_PLAT}" = sles-11-broadwell-64-xc40 ];then
    OPTF="${OPTF}"
else
    echo "This platform 'ORDENV_PLAT=${ORDENV_PLAT}' is not supported.  Only 'ubuntu-14.04-amd64-64' and 'sles-11-amd64-64' are."
    exit 1
fi

if [ "${COMPILE_OAVAR_ADD_DEBUG_OPTIONS:-no}" = yes ]; then
    echo "Using debug options"
    COMPF_NOC="${COMPF_GLOBAL} -debug DEBUG -optf ${OPTF}"
    COMPF="${COMPF_NOC} =-C"
else
    COMPF="${COMPF_GLOBAL} -optf ${OPTF}"
    COMPF_NOC=${COMPF}
fi

if [ "${mode}" == full ] ; then

  rm -f *.o *.mod *.cdk* *.h *.ftn* *.f *.f90

  # Create a local copy of the source code
  sed "s!XXXXX!${revnum}!g" ${trunkdir}/toplevelcontrol_mod.ftn90_template > toplevelcontrol_mod.ftn90

  cd ${trunkdir};          ls -1F | grep -v '/' | grep -v "*" | grep -v "@" | cpio -pl $compiledir ; cd $compiledir
  cd ${trunkdir}/bgcheck;  ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir
  cd ${trunkdir}/shared;   ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir
  cd ${trunkdir}/modulopt; ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir
  rm -f *.ftn~ *.ftn90~

  echo "STARTING COMPILATION AT:" 
  date

  # Remove enkf_pturb.ftn main program from compilation directory
  rm -f enkf_pturb.ftn

  # Compile the subroutines...
  echo "compiling low-level independent modules"
  echo "If aborting, check in ${PWD}/listing1"
  SRC0="toplevelcontrol_mod.ftn90"
  SRC0="$SRC0 mathphysconstants_mod.ftn90 earthconstants_mod.ftn90 mpi_mod.ftn90 mpivar_mod.ftn90 bufr_mod.ftn90 codtyp_mod.ftn90"
  SRC0="$SRC0 physicsfunctions_mod.ftn90 obsspacedata_mod.ftn90 horizontalcoord_mod.ftn90 timecoord_mod.ftn90 verticalcoord_mod.ftn90"
  s.compile $COMPF -O ${FOPTMIZ} -src $SRC0 > listing1 2>&1
  status=1
  grep fail listing1 || status=0
  if [ "${status}" -ne 0 ] ; then
      echo "Compilation aborted: check in ${PWD}/listing1"
      exit 1
  fi

  echo "compiling modulopt (n1qn3)"
  echo "If aborting, check in ${PWD}/listing0"
  SRC0="dcube.ftn ddd.ftn ddds.ftn dystbl.ftn mupdts.ftn n1qn3.ftn n1qn3a.ftn nlis0.ftn"
  s.compile $COMPF_NOC  -O ${FOPTMIZ} -src $SRC0 > listing0 2>&1
  status=1
  grep fail listing0 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "Compilation aborted: check in ${PWD}/listing0"
      exit 1
  fi
  rm -f $SRC0

  echo "compiling analysis grid modules"
  echo "If aborting, check in ${PWD}/listing2"
  SRC0="gaussgrid_mod.ftn90 windrotation_mod.ftn90 lamanalysisgrid_mod.ftn90"
  s.compile $COMPF  -O ${FOPTMIZ} -src $SRC0 > listing2 2>&1
  status=1
  grep fail listing2 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "Compilation aborted: check in ${PWD}/listing2"
      exit 1
  fi

  echo "compiling most of the new modules"
  echo "If aborting, check in ${PWD}/listing3"
  SRC1="controlvector_mod.ftn90 rmatrix_mod.ftn90 hir_chans_mod.ftn90 tovs_nl_mod.ftn90 tovs_lin_mod.ftn90 varnamelist_mod.ftn90 columndata_mod.ftn90 multi_ir_bgck_mod.ftn90"
  SRC1="$SRC1 emissivities_mod.ftn90 fft_mod.ftn90 globalspectraltransform_mod.ftn90"
  SRC1="$SRC1 lamspectraltransform_mod.ftn90  gridstatevector_mod.ftn90"
  SRC1="$SRC1 bmatrixensemble_mod.ftn90 bmatrixhi_mod.ftn90 lambmatrixhi_mod.ftn90"
  SRC1="$SRC1 bmatrix_mod.ftn90 minimization_mod.ftn90"
  SRC1="$SRC1 ozoneclim_mod.ftn90 tovs_extrap_mod.ftn90"
  SRC1="$SRC1 burpfiles_mod.ftn90 obsspacediag_mod.ftn90 observation_erreurs_mod.ftn90"

  s.compile $COMPF  -O ${FOPTMIZ} -src $SRC1 > listing3 2>&1
  status=1
  grep fail listing3 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "Compilation aborted: check in ${PWD}/listing3"
      exit 1
  fi

  echo "compiling burp_read module"
  echo "If aborting, check in ${PWD}/listing4"
  SRC1="burp_read_mod.ftn90 burp_functions.ftn90 selectb.ftn90 update_burpfiles.ftn90"
  s.compile $COMPF  -O ${FOPTMIZ} -src $SRC1 > listing4 2>&1
  status=1
  grep fail listing4 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "Compilation aborted: check in ${PWD}/listing4"
      exit 1
  fi
  
  echo "compiling the GPS modules (cdk90)..."
  echo "If aborting, check in ${PWD}/listing5"
  SRC2="modgps00base.cdk90 modgps01ctmath.cdk90 modgps01ctphys.cdk90 modgps02wgs84const.cdk90 modgps02wgs84grav.cdk90 modgps03diff.cdk90 modgps04profile.cdk90"
  SRC2="$SRC2 modgps05refstruct.cdk90 modgps07geostruct.cdk90 modgps08refop.cdk90 modgps09bend.cdk90 modgps04profilezd.cdk90"
  SRC2="$SRC2 modgps08ztdop.cdk90"
  s.compile $COMPF  -O ${FOPTMIZ} -src $SRC2 > listing5 2>&1
  status=1
  grep fail listing5 || status=0
  if [ "${status}" -ne 0 ] ; then
      echo "Compilation aborted: check in ${PWD}/listing5"
      exit 1
  fi

  echo "compiling some more modules..."
  echo "If aborting, check in ${PWD}/listing6"
  SRC2="modgpsro_mod.ftn90 modgpsztd_mod.ftn90 filterobs_mod.ftn90 writeincrement_mod.ftn90 obsoperators_mod.ftn90"
  s.compile $COMPF  -O ${FOPTMIZ} -src $SRC2 > listing6 2>&1
  status=1
  grep fail listing6 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "Compilation aborted: check in ${PWD}/listing6"
      exit 1
  fi

  echo "compiling remaining ftn ftn90..."
  echo "If aborting, check in ${PWD}/listing7"
  filelist=""
  for i in *.ftn *.ftn90 ; do
      if [[ "${i}" != *_mod.ftn* ]]; then
	  filelist="$filelist ${i}"
      fi
  done
  s.compile $COMPF  -O ${FOPTMIZ} -src $filelist > listing7 2>&1
  status=1
  grep fail listing7 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "Compilation aborted: check in ${PWD}/listing7"
      exit 1
  fi

  echo "building the executable ${varabs}.Abs"
  rm -f ${varabs}.Abs
  echo "If aborting, check in ${PWD}/listing8"
  s.compile $COMPF  -O ${FOPTMIZ} -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o ${varabs}.Abs > listing8 2>&1
  status=1
  grep fail listing8 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "Compilation aborted: check in ${PWD}/listing8"
      exit 1
  fi

  status=1
  grep -i ERROR listing? || status=0
  if [ "${status}" -ne 0 ] ; then
      echo "ERROR found: STOP; check listing in ${PWD}"
      exit 1
  fi

  #rm -f *.ftn* *.f *.f90

  echo "FINISHED COMPILATION AT:"
  date

elif [ "${mode}" == abs ] ; then

  rm -f ${varabs}.Abs

  echo
  echo "building the executable..."
  echo
  set -x
  s.compile $COMPF  -O ${FOPTMIZ} -libappl $LIBAPPL $LIBEXTRA -librmn $LIBRMN -obj *.o -o ${varabs}.Abs

else
    if [ -f $trunkdir/$mode ] ; then
	file=`basename $mode`
	rm -f $file
	cp $trunkdir/$mode .
	s.compile $COMPF  -O ${FOPTMIZ} -src $file
    else
	echo "File $trunkdir/$mode does NOT exist. Stop"
	exit 1
    fi
fi

echo "FINISHED COMPILATION AT:"
date
echo "The program can be found here: ${PWD}/${varabs}.Abs"
