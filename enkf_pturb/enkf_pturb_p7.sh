#!/bin/ksh
#

# ----------------------------------------
# Launcher script for enkf_pturb on IBM p7
# ----------------------------------------

#
# User-defined options
#
flnml="namelist_p7.nml"
machine="hadar"
gest="/users/dor/arma/bue/power7/3dvar_modular/enkf_pturb_stag_Bhyb_enkflvls/"
anlgrid="/users/dor/arma/gr3/data_cnfs/prototype_grid/analysisgrid_glb_800x400"
#cov="spica:/home/dormrb02/ibmenv/armnlib/modeles/ANAL/stats/auto/__GEM25km_NMC_T399_stag5002_BgckStddev3d_800x400__/01"
#cov="datasvr:/users/dor/arma/bue/cnfs/test_enkf_ptb/glbcov01_stag5002"
cov="datasvr:/users/dor/arma/bue/cnfs/test_enkf_ptb/enkf_analysis_zap"
abs="/users/dor/arma/bue/home01/3dvar_latest/compiledir_enkf_pturb/enkf_pturb_p7.abs"
npex=8
npey=2
openmp=2
maxcputime=300

#
# Don't modify below ...
#
abs_basename=`basename $abs`

echo
echo "Launching ENKF_PTURB using..."
echo
echo "Using namelist file  :" $flnml
echo "Working machine      :" $machine
echo "Working directory    :" $gest
echo "Analysis grid        :" $anlgrid
echo "Executable file path :" $abs
echo "Executable file name :" $abs_basename
echo "Topology             :" ${npex}x${npey}x${openmp}
echo

ssh $machine mkdir -p $gest 
ssh $machine rm -f $gest/*
scp $flnml ${machine}:${gest}/flnml
scp $anlgrid ${machine}:${gest}/analysisgrid
scp $cov ${machine}:${gest}/bgcov
scp $abs ${machine}:${gest}/enkf_pturb.abs
ssh $machine ls -l $gest

cat << EOF > go_enkf_pturb.sh
 echo "!!STARTING SCRIPT!!"
. ssmuse-sh -d rpn/utils/15.1
 cd $gest
 export TMG_ON=YES
 export MP_STDOUTMODE=unordered
 r.run_in_parallel -pgm ./enkf_pturb.abs -npex ${npex} -npey ${npey}
EOF


cat << EOF > ptopo_nml
 &ptopo
  npex=$npex
  npey=$npey
/
EOF
scp ptopo_nml ${machine}:${gest}

ord_soumet go_enkf_pturb.sh -mach $machine -mpi -t $maxcputime -cm 1632M -cpus ${npex}x${npey}x${openmp} -listing ${gest} -jn enkf_pturb
