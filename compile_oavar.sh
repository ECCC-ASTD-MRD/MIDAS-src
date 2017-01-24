#!/bin/ksh

# Revisions:
#            Mike Sitwell, Y.J. Rochon, and Ping Du (as of Dec 2014)
#            - The following are changes/additions for chemical constituents.
#            - Added compilation of bmatrixchem_mod.ftn90
#            - Addition of chem_mod.ftn90 and chem_interface_mod.ftn90

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
  MPILIBDIR=""
  MPILIB="rpn_commstubs rpn_comm"
  MPIKEY=""
  ABSTAG="_nompi"
else
  echo "..."
  echo "... > Compiling for an MPI executable"
  # JFC : Mesure temporaire pour avoir acces a la S-R RPN_COMM_adj_halo8 de M. Valin
  if [ "${BASE_ARCH}" = "AIX-powerpc7" ] ; then
      MPILIBDIR="-libpath /users/dor/arma/anl/userlibs/${BASE_ARCH}/xlf13"
  elif [ "${BASE_ARCH}" = "Linux_x86-64" ] ; then
      MPILIBDIR="-libpath /users/dor/arma/anl/userlibs/${BASE_ARCH}/intel13sp1u2"
  else
    echo "..."
    echo "... !! This platform 'ARCH=${ARCH}' is not supported.  Only 'AIX-powerpc7' and 'Linux_x86-64' are !!"
    exit 1
  fi
  MPILIB="_anl_rpn_comm_4051103"
  # JFC : Mesure temporaire FIN
  MPIKEY="-mpi"
  ABSTAG=""
fi

trunkdir=$PWD
if [ "${BASE_ARCH}" = "AIX-powerpc7" ];then
    FOPTMIZ=2
elif [ "${BASE_ARCH}" = "Linux_x86-64" ];then
    FOPTMIZ=1
else
    echo "..."
    echo "... !! This platform 'ARCH=${ARCH}' is not supported.  Only 'AIX-powerpc7' and 'Linux_x86-64' are !!"
    exit 1
fi

# automatically set the global revision number in toplevelcontrol_mod.ftn90 by
# replacing the string XXXXX with the actual revision number
revnum=$(git describe --always --dirty=_M 2>/dev/null || ssh pollux "cd $trunkdir; git describe --always --dirty=_M" 2>/dev/null || echo unkown revision)
echo "..."
echo "... > Revision Number = '$revnum'"

# Set compiledir
compiledir_main=${COMPILEDIR_OAVAR_MAIN:-".."}
compiledir_ID=${COMPILEDIR_OAVAR_ID:-$revnum}
compiledir=${compiledir_main}/compiledir_${compiledir_ID}

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

varabs=oavar_${BASE_ARCH}${ABSTAG}

#----------------------------------------------------------------
#  Set up dependent librarys and tools. 
#---------------------------------------------------------------
echo "..."
echo "... > Setting the environement"
# for s.compile
echo "...   loading hpcs/201402/02/base"
. ssmuse-sh -d hpcs/201402/02/base
## for the compiler
if [ "${BASE_ARCH}" = "AIX-powerpc7" ];then
    echo "...   loading compiler hpcs/ext/xlf_13.1.0.10"
    . ssmuse-sh -d hpcs/ext/xlf_13.1.0.10
    # NetCDF for the IBM:
    echo "...   loading netcdf"
    export EC_LD_LIBRARY_PATH="/ssm/net/rpn/mfv/netcdf4/lib ${EC_LD_LIBRARY_PATH}"
    export EC_INCLUDE_PATH="/ssm/net/rpn/mfv/netcdf4/include ${EC_INCLUDE_PATH}"
    CDF_LIBS="netcdf netcdff hdf5 hdf5_hl sz z"
elif [ "${BASE_ARCH}" = "Linux_x86-64" ];then
    echo "...   loading compiler hpcs/201402/02/intel13sp1u2 (this includes the netcdf library)"
    . ssmuse-sh -d hpcs/201402/02/intel13sp1u2
    . s.ssmuse.dot dot
    echo "...   EC_LD_LIBRARY_PATH=${EC_LD_LIBRARY_PATH}"
    echo "...   EC_INCLUDE_PATH=${EC_INCLUDE_PATH}"
    CDF_LIBS=netcdff
else
    echo "... !! This platform 'ARCH=${ARCH}' is not supported.  Only 'AIX-powerpc7' and 'Linux_x86-64' are !!"
    exit 1
fi
## for rmn_015, rpncomm
echo "...   loading rpn/libs/15.2"
. ssmuse-sh -d rpn/libs/15.2
## for 'vgrid'
echo "...   loading cmdn/vgrid/5.3.2/${COMP_ARCH}"
. ssmuse-sh -d cmdn/vgrid/5.3.2/${COMP_ARCH}
## for 'burplib'
echo "...   loading cmda/base/201411/01/${COMP_ARCH}"
. ssmuse-sh -d cmda/base/201411/01/${COMP_ARCH}

## For hpcsperf needed for TMG timings
echo "...   loading hpcs/exp/aspgjdm/perftools"
. ssmuse-sh -d hpcs/exp/aspgjdm/perftools
# For RTTOV package...  
echo "...   loading arma/rttov/10v4"
. ssmuse-sh -d arma/rttov/10v4

#-----------------------------------------------------------------------------

LIBAPPL="${CDF_LIBS} rttov10.2.0_coef_io rttov10.2.0_main rttov10.2.0_emis_atlas rttov10.2.0_other burp_module descrip $MPILIB"
if [ "${BASE_ARCH}" = "AIX-powerpc7" ];then
    LIBSYS="hpcsperf lapack-3.4.0 essl mass"
elif [ "${BASE_ARCH}" = "Linux_x86-64" ];then
    LIBSYS="hpcsperf lapack blas"
else
    echo "... !! This platform 'BASE_ARCH=${BASE_ARCH}' is not supported.  Only 'AIX-powerpc7' and 'Linux_x86-64' are !!"
    exit 1
fi

LIBRMN=rmn

COMPF_GLOBAL="-openmp ${MPIKEY}"
if [ "${BASE_ARCH}" = "AIX-powerpc7" ];then 
    if [ "${COMPILE_OAVAR_REMOVE_DEBUG_OPTIONS}" = yes ]; then
	COMPF=${COMPF_GLOBAL}
	COMPF_NOC=${COMPF}
	echo "..."
	echo "... > Using fully optimized compilation"
    else
	COMPF_NOC="${COMPF_GLOBAL} -debug DEBUG"
	COMPF="${COMPF_NOC} -optf =-C"
	echo "..."
	echo "... > !WARNING! You are compiling in DEBUG MODE"
    fi
elif [ "${BASE_ARCH}" = "Linux_x86-64" ];then
    OPTF="=-mkl =-fp-model =source =-check =noarg_temp_created"
    if [ "${COMPILE_OAVAR_REMOVE_DEBUG_OPTIONS}" = yes ]; then
	COMPF="${COMPF_GLOBAL} -optf ${OPTF}"
	COMPF_NOC=${COMPF}
	echo "..."
	echo "... > Using fully optimized compilation"
    else
	COMPF_NOC="${COMPF_GLOBAL} -debug DEBUG -optf ${OPTF}"
	COMPF="${COMPF_NOC} =-C"
	echo "..."
	echo "... > !WARNING! You are compiling in DEBUG MODE"
    fi
else 
    echo "..."
    echo "... !! This platform 'BASE_ARCH=${BASE_ARCH}' is not supported.  Only 'AIX-powerpc7' and 'Linux_x86-64' are !!"
    exit 1
fi

echo "..."
echo "... > STARTING COMPILATION AT: $(date)"

if [ "${mode}" == full ] ; then

  rm -f *.o *.mod *.cdk* *.h *.ftn* *.f *.f90

  # Create a local copy of the source code
  sed "s!XXXXX!${revnum}!g" ${trunkdir}/toplevelcontrol_mod.ftn90_template > toplevelcontrol_mod.ftn90

  cd ${trunkdir};          ls -1F | grep -v '/' | grep -v "*" | grep -v "@" | cpio -pl $compiledir ; cd $compiledir
  cd ${trunkdir}/bgcheck;  ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir
  cd ${trunkdir}/shared;   ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir
  cd ${trunkdir}/modulopt; ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir
  rm -f *.ftn~ *.ftn90~

  # Compile the subroutines...
  echo "... > Compiling low-level independent modules"
  echo "...   if aborting, check in ${PWD}/listing1"
  SRC0="toplevelcontrol_mod.ftn90"
  SRC0="$SRC0 mathphysconstants_mod.ftn90 earthconstants_mod.ftn90 mpi_mod.ftn90 mpivar_mod.ftn90 bufr_mod.ftn90 codtyp_mod.ftn90"
  SRC0="$SRC0 physicsfunctions_mod.ftn90 obsspacedata_mod.ftn90 localizationfunction_mod.ftn90"
  SRC0="$SRC0 horizontalcoord_mod.ftn90 timecoord_mod.ftn90 verticalcoord_mod.ftn90"
  SRC0="$SRC0 lqtoes_mod.ftn90 intavg_mod.ftn90 filterresponsefunction_mod.ftn90"
  s.compile $COMPF -O ${FOPTMIZ} -src $SRC0 > listing1 2>&1
  status=1
  grep fail listing1 || status=0
  if [ "${status}" -ne 0 ] ; then
      echo "... !! Compilation aborted: check in ${PWD}/listing1 !!"
      exit 1
  fi

  echo "... > Compiling modulopt (n1qn3)"
  echo "...   if aborting, check in ${PWD}/listing0"
  SRC0="dcube.ftn ddd.ftn ddds.ftn dystbl.ftn mupdts.ftn n1qn3.ftn n1qn3a.ftn nlis0.ftn"
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
  SRC0="windrotation_mod.ftn90 lamanalysisgrid_mod.ftn90"
  s.compile $COMPF  -O ${FOPTMIZ} -src $SRC0 > listing2 2>&1
  status=1
  grep fail listing2 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "... !! Compilation aborted: check in ${PWD}/listing2 !!"
      exit 1
  fi

  echo "... > Compiling most of the new modules"
  echo "...   if aborting, check in ${PWD}/listing3"
  SRC1="controlvector_mod.ftn90 rmatrix_mod.ftn90 hir_chans_mod.ftn90 tovs_nl_mod.ftn90"
  SRC1="$SRC1 tovs_lin_mod.ftn90 varnamelist_mod.ftn90 columndata_mod.ftn90 multi_ir_bgck_mod.ftn90"
  SRC1="$SRC1 emissivities_mod.ftn90 globalspectraltransform_mod.ftn90 tt2phi_mod.ftn90"
  SRC1="$SRC1 lamspectraltransform_mod.ftn90 gridstatevector_mod.ftn90 statetocolumn_mod.ftn90 variabletransforms_mod.ftn90"
  SRC1="$SRC1 bmatrixensemble_mod.ftn90 bmatrixhi_mod.ftn90 lambmatrixhi_mod.ftn90 bmatrixchem_mod.ftn90"
  SRC1="$SRC1 bmatrix_mod.ftn90"
  SRC1="$SRC1 chem_mod.ftn90 chem_interface_mod.ftn90"
  SRC1="$SRC1 ozoneclim_mod.ftn90 tovs_extrap_mod.ftn90"
  SRC1="$SRC1 obsspacediag_mod.ftn90 observation_erreurs_mod.ftn90"

  s.compile $COMPF  -O ${FOPTMIZ} -src $SRC1 > listing3 2>&1
  status=1
  grep fail listing3 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "... !! Compilation aborted: check in ${PWD}/listing3 !!"
      exit 1
  fi

  echo "... > Compiling burp_read module"
  echo "...   if aborting, check in ${PWD}/listing4"
  SRC1="burp_read_mod.ftn90"
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
  SRC2="varqc_mod.ftn90 filterobs_mod.ftn90 obsoperators_mod.ftn90"
  SRC2="$SRC2 minimization_mod.ftn90 burpfiles_mod.ftn90"
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
  for i in *.ftn *.ftn90 ; do
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

  echo "... > Building the executable ${varabs}.Abs"
  rm -f ${varabs}.Abs
  echo "...   if aborting, check in ${PWD}/listing8"
  s.compile $COMPF  -O ${FOPTMIZ} ${MPILIBDIR} -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o ${varabs}.Abs > listing8 2>&1
  status=1
  grep fail listing8 || status=0
  if [ "${status}" -ne 0 ]; then
      echo "... !! Compilation aborted: check in ${PWD}/listing8 !!"
      exit 1
  fi

  status=1
  grep -i ERROR listing? || status=0
  if [ "${status}" -ne 0 ] ; then
      echo "... !! ERROR found: STOP; check listing in ${PWD} !!"
      exit 1
  fi

  #### rm -f *.ftn* *.f *.f90

elif [ "${mode}" == abs ] ; then

  rm -f ${varabs}.Abs

  echo "..."
  echo "... building the executable ${varabs}.Abs"
  echo "..."
  s.compile $COMPF  -O ${FOPTMIZ} ${MPILIBDIR} -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o ${varabs}.Abs

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
    echo "... The program can be found here: ${PWD}/${varabs}.Abs"
    echo "..."
else
    echo "..."
fi

echo "...             |==================================|"
echo "... ------------|  VAR compilation script ENDING   |------------"
echo "...             |==================================|"
