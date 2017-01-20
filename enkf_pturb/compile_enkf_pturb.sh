#!/bin/bash

mode=$1
#nompi=$2

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

echo " !!! Compiling for an MPI executable !!!"
echo ""
MPILIB="rpn_comm"
MPIKEY="-mpi"
ABSTAG=""

enkfpturbdir=$PWD
trunkdir=$PWD/../

if [ "${ORDENV_PLAT}" = sles-11-haswell-64-xc40 ];then
    echo "Switching ORDENV_PLAT from '${ORDENV_PLAT}' to 'sles-11-broadwell-64-xc40'"
    . ../r.env.dot --arch sles-11-broadwell-64-xc40
    echo ORDENV_PLAT=${ORDENV_PLAT}
fi

# automatically set the global revision number in toplevelcontrol_mod.ftn90 by
# replacing the string XXXXX with the actual revision number
revnum=$(git describe --always --dirty 2>/dev/null || ssh eccc-ppp1 "cd $trunkdir; git describe --always --dirty" 2>/dev/null || echo unkown revision)
echo " "
echo "-----------------------"
echo "Revision number='$revnum'"
echo "-----------------------"
echo " "

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

LIBAPPL="rttov10.2.0_coef_io rttov10.2.0_main rttov10.2.0_other burp_module descrip $MPILIB"

LIBSYS=""
LIBRMN="rmn"
LIBEXTRA="hpcoperf"

COMPF_GLOBAL="-openmp $MPIKEY -O"

OPTF="=-check =noarg_temp_created"
if [ "${ORDENV_PLAT}" = ubuntu-14.04-amd64-64 ];then
    OPTF="=-mkl ${OPTF}"
elif [ "${ORDENV_PLAT}" = sles-11-amd64-64 -o "${ORDENV_PLAT}" = sles-11-broadwell-64-xc40 ];then
    OPTF="${OPTF}"
else
    echo "This platform 'ORDENV_PLAT=${ORDENV_PLAT}' is not supported.  Only 'ubuntu-14.04-amd64-64' and 'sles-11-amd64-64' are."
    exit 1
fi

if [ "${COMPILE_ENKF_PTURB_ADD_DEBUG_OPTIONS:-no}" = yes ]; then
    echo "Using debug options"
    COMPF_NOC="${COMPF_GLOBAL} -debug DEBUG -optf ${OPTF}"
    COMPF="${COMPF_NOC} =-C"
else
    COMPF="${COMPF_GLOBAL} -optf ${OPTF}"
    COMPF_NOC=${COMPF}
fi

# Create and Move to compilation directory
mkdir -p ../../compiledir_enkf_pturb-${ORDENV_PLAT}
cd ../../compiledir_enkf_pturb-${ORDENV_PLAT}
compiledir=${PWD}

if [ $mode == full ] ; then

  rm -f *.o *.mod *.cdk* *.h *.ftn* *.f *.f90

  # Create a local copy of the source code
  cd ${trunkdir};            ls -1F | grep -v '/' | grep -v "*" | grep -v "@" | cpio -pl $compiledir ; cd $compiledir
  cd ${trunkdir}/shared;     ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir
  cd ${trunkdir}/enkf_pturb; ls -1F | grep -v '/' | grep -v "*" | cpio -pl $compiledir ; cd $compiledir
  rm -f *.ftn~ *.ftn90~

  echo "STARTING COMPILATION AT:" 
  date

  cat ${trunkdir}/toplevelcontrol_mod.ftn90_template |sed "s!XXXXX!${revnum}!g" > toplevelcontrol_mod.ftn90

  # Compile the subroutines...
  echo "compiling low-level independent modules"
  SRC0="toplevelcontrol_mod.ftn90"
  SRC0="$SRC0 mathphysconstants_mod.ftn90 earthconstants_mod.ftn90 mpi_mod.ftn90 mpivar_mod.ftn90 bufr_mod.ftn90 codtyp_mod.ftn90"
  SRC0="$SRC0 physicsfunctions_mod.ftn90 obsspacedata_mod.ftn90 horizontalcoord_mod.ftn90 timecoord_mod.ftn90 verticalcoord_mod.ftn90"
  s.compile $COMPF -src $SRC0 > listing1 2>&1
  grep fail listing1
  if [ $? = "0" ] ; then exit ; fi

  echo "compiling analysis grid modules"
  SRC0="gaussgrid_mod.ftn90 windrotation_mod.ftn90 lamanalysisgrid_mod.ftn90"
  s.compile $COMPF -src $SRC0 > listing2 2>&1
  grep fail listing2
  if [ $? = "0" ] ; then exit ; fi

  echo "compiling most of the new modules"
  SRC1="controlvector_mod.ftn90 fft_mod.ftn90"
  SRC1="$SRC1 globalspectraltransform_mod.ftn90 varnamelist_mod.ftn90"
  SRC1="$SRC1 lamspectraltransform_mod.ftn90 gridstatevector_mod.ftn90"
  SRC1="$SRC1 bmatrixensemble_mod.ftn90 bmatrixhi_mod.ftn90 lambmatrixhi_mod.ftn90"
  SRC1="$SRC1 bmatrix_mod.ftn90 writeincrement_mod.ftn90"
  s.compile $COMPF -src $SRC1 > listing3 2>&1
  grep fail listing3
  if [ $? = "0" ] ; then exit ; fi

  echo "compiling remaining ftn ftn90..."
  filelist="abort.ftn utils_3dvar.ftn getstamplist.ftn90 matsqrt.ftn getfldprm.ftn getfldprm2.ftn"
  filelist="$filelist getstepobsindex.ftn90 matapat.ftn vintgd.ftn90 initgdg2.ftn90 gasdev.ftn"
  filelist="$filelist enkf_pturb.ftn90"
  s.compile $COMPF -src $filelist > listing7 2>&1
  grep fail listing7
  if [ $? = "0" ] ; then exit ; fi

  echo "building the executable..."
  s.compile $COMPF -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o enkf_pturb.abs$ABSTAG > listing8 2>&1

  grep -iE 'ERROR|ERREUR' listing?
  if [ $? = "0" ] ; then echo "ERROR found: STOP" ; exit ; fi

  rm -f *.ftn* *.f *.f90

elif [ $mode == abs ] ; then

  rm -f enkf_pturb.abs$ABSTAG

  echo
  echo "building the executable..."
  echo
  s.compile $COMPF -libappl $LIBAPPL $LIBEXTRA -libsys $LIBSYS -librmn $LIBRMN -obj *.o -o enkf_pturb.abs$ABSTAG

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

echo "FINISHED COMPILATION AT:"
date
echo "The program can be found here: ${PWD}/enkf_pturb.abs$ABSTAG"
