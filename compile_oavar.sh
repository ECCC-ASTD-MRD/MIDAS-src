#!/bin/bash

# Revisions:
#            Mike Sitwell, Y.J. Rochon, and Ping Du (as of Dec 2014)
#            - The following are changes/additions for chemical constituents.
#            - Added compilation of bmatrixchem_mod.ftn90
#            - Addition of obssubspacedata_mod.ftn90 and chem_*_mod.ftn90 

set -e

mode=$1
nompi=$2

echo "..."
echo "...             |==================================|"
echo "... ------------|  VAR compilation script STARTING |------------"
echo "...             |==================================|"

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
  MPIKEY=""
  ABSTAG="_nompi"
else
  echo "..."
  echo "... > Compiling for an MPI executable"
  MPILIB="rpn_comm"
  MPIKEY="-mpi"
  ABSTAG=""
fi

trunkdir=$PWD

if [ "${ORDENV_PLAT}" = sles-11-haswell-64-xc40 ];then
    echo "... Switching ORDENV_PLAT from '${ORDENV_PLAT}' to 'sles-11-broadwell-64-xc40'"
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
    echo "... This platform 'ORDENV_PLAT=${ORDENV_PLAT}' is not supported.  Only 'ubuntu-14.04-amd64-64' and 'ubuntu-14.04-amd64-64' are."
    exit 1
fi

# automatically set the global revision number in toplevelcontrol_mod.ftn90 by
# replacing the string XXXXX with the actual revision number
revnum=$(git describe --always --dirty=_M 2>/dev/null || ssh eccc-ppp1 "cd $trunkdir; git describe --always --dirty=_M" 2>/dev/null || echo unkown revision)
echo "..."
echo "... > Revision Number = '$revnum'"

# Set compiledir
compiledir_main=${COMPILEDIR_OAVAR_MAIN:-".."}
compiledir_ID=${COMPILEDIR_OAVAR_ID:-$revnum}
compiledir=${compiledir_main}/compiledir-${ORDENV_PLAT}_${compiledir_ID}
mkdir -p $compiledir
cd $compiledir
compiledir=${PWD} # needed when compiledir_main = ".."

echo "..."
echo "... > Compiledir set to $compiledir"
if [ ${compiledir_main} != ".." ] ; then
    if [ ! -d  ${trunkdir}/../compiledir_${compiledir_ID} ] ; then
	ln -s ${compiledir_main}/compiledir_${compiledir_ID} ${trunkdir}/../compiledir_${compiledir_ID}
    fi
fi


#----------------------------------------------------------------
#  Set up dependent librarys and tools. 
#---------------------------------------------------------------
## for s.compile
echo "... loading hpco/tmp/eccc/201402/06/base"
. ssmuse-sh -d hpco/tmp/eccc/201402/06/base
## for the compiler
if [ "${ORDENV_PLAT}" = ubuntu-14.04-amd64-64 ];then
    echo "... loading compiler main/opt/intelcomp/intelcomp-2016.1.156"
    . ssmuse-sh -d main/opt/intelcomp/intelcomp-2016.1.156
elif [ "${ORDENV_PLAT}" = sles-11-amd64-64 -o "${ORDENV_PLAT}" = sles-11-broadwell-64-xc40 ];then
    echo "... loading compiler PrgEnv-intel-5.2.82"
    module load PrgEnv-intel/5.2.82
else
    echo "... This platform 'ORDENV_PLAT=${ORDENV_PLAT}' is not supported.  Only 'ubuntu-14.04-amd64-64' and 'sles-11-amd64-64' are."
    exit 1
fi

## for netcdf
CDF_LIBS=netcdff

varabs=oavar_${ORDENV_PLAT}${ABSTAG}-${revnum}.Abs

## for rmn, rpncomm
echo "... loading eccc/mrd/rpn/libs/16.1"
. ssmuse-sh -d eccc/mrd/rpn/libs/16.1
if [ "${ORDENV_PLAT}" = ubuntu-14.04-amd64-64 ];then
    ## for openmpi
    echo "... loading main/opt/openmpi/openmpi-1.6.5/intelcomp-2016.1.156"
    . ssmuse-sh -d main/opt/openmpi/openmpi-1.6.5/intelcomp-2016.1.156
elif [ "${ORDENV_PLAT}" = sles-11-amd64-64 -o "${ORDENV_PLAT}" = sles-11-broadwell-64-xc40 ];then
    echo "loading cray-netcdf"
    module load cray-netcdf
fi
## for 'vgrid'
echo "... loading eccc/cmd/cmdn/vgrid/5.6.9/${COMP_ARCH}"
. ssmuse-sh -d eccc/cmd/cmdn/vgrid/5.6.9/${COMP_ARCH}
## for 'burplib'
echo "... loading eccc/cmd/cmda/libs/16.1-3/${COMP_ARCH}"
. ssmuse-sh -d eccc/cmd/cmda/libs/16.1-3/${COMP_ARCH}
## For hpcoperf needed for TMG timings
echo "... loading main/opt/perftools/perftools-2.0/${COMP_ARCH}"
. ssmuse-sh -d main/opt/perftools/perftools-2.0/${COMP_ARCH}
# For RTTOV 10v3 package...
echo "... loading eccc/mrd/rpn/anl/rttov/10v3.2/${COMP_ARCH}"
. ssmuse-sh -d eccc/mrd/rpn/anl/rttov/10v3.2/${COMP_ARCH}

#-----------------------------------------------------------------------------

LIBAPPL="${CDF_LIBS} rttov10.2.0_coef_io rttov10.2.0_main rttov10.2.0_emis_atlas rttov10.2.0_other burp_module descrip $MPILIB"
LIBSYS="hpcoperf"

LIBRMN=rmnMP

COMPF_GLOBAL="-openmp ${MPIKEY}"
OPTF="=-check =noarg_temp_created"
if [ "${ORDENV_PLAT}" = ubuntu-14.04-amd64-64 ];then
    OPTF="=-mkl ${OPTF}"
elif [ "${ORDENV_PLAT}" = sles-11-amd64-64 -o "${ORDENV_PLAT}" = sles-11-broadwell-64-xc40 ];then
    OPTF="${OPTF}"
else
    echo "... This platform 'ORDENV_PLAT=${ORDENV_PLAT}' is not supported.  Only 'ubuntu-14.04-amd64-64' and 'sles-11-amd64-64' are."
    exit 1
fi

if [ "${COMPILE_OAVAR_ADD_DEBUG_OPTIONS:-no}" = yes ]; then
    echo "... > !WARNING! You are compiling in DEBUG MODE"
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
  echo "...   if aborting, check in ${PWD}/listing1"
  SRC0="mathphysconstants_mod.ftn90 earthconstants_mod.ftn90 utilities_mod.ftn90"
  SRC0="$SRC0 toplevelcontrol_mod.ftn90 randomnumber_mod.ftn90 mpi_mod.ftn90 mpivar_mod.ftn90 bufr_mod.ftn90 codtyp_mod.ftn90"
  SRC0="$SRC0 physicsfunctions_mod.ftn90 obsspacedata_mod.ftn90 localizationfunction_mod.ftn90"
  SRC0="$SRC0 horizontalcoord_mod.ftn90 timecoord_mod.ftn90 verticalcoord_mod.ftn90"
  SRC0="$SRC0 lqtoes_mod.ftn90 presprofileoperators_mod.ftn90 spectralfilter_mod.ftn90"
  s.compile $COMPF -O ${FOPTMIZ} -src $SRC0 > listing1 2>&1
  status=1
  grep fail listing1 || status=0
  if [ "${status}" -ne 0 ] ; then
      echo "... !! Compilation aborted: check in ${PWD}/listing1 !!"
      exit 1
  fi

  echo "... > Compiling quasi-newton module"
  echo "...   if aborting, check in ${PWD}/listing0"
  SRC0="quasinewton_mod.ftn"
  s.compile $COMPF_NOC  -O ${FOPTMIZ} -src $SRC0 > listing0 2>&1
  status=1
  grep fail listing0 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "... !! Compilation aborted: check in ${PWD}/listing0 !!"
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
  SRC1="$SRC1 variabletransforms_mod.ftn90 bmatrixensemble_mod.ftn90 bmatrixhi_mod.ftn90 lambmatrixhi_mod.ftn90"
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
      if [[ "${i}" != *_mod.ftn* ]]; then
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

  echo "... > Building the executable ${varabs}"
  rm -f ${varabs}
  echo "...   if aborting, check in ${PWD}/listing8"
  s.compile $COMPF  -O ${FOPTMIZ} ${MPILIBDIR} -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o ${varabs} > listing8 2>&1
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
  cp ${varabs} ~/ords/oavar_abs/

  #### rm -f *.ftn* *.f *.f90

elif [ "${mode}" == abs ] ; then

  rm -f ${varabs}

  echo "..."
  echo "... > Building the executable ${varabs}"
  echo "..."
  s.compile $COMPF  -O ${FOPTMIZ} ${MPILIBDIR} -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o ${varabs}
  cp ${varabs} ~/ords/oavar_abs/

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
    echo "... The program can be found here: ${PWD}/${varabs}"
    echo "..."
else
    echo "..."
fi

echo "...             |==================================|"
echo "... ------------|  VAR compilation script ENDING   |------------"
echo "...             |==================================|"
