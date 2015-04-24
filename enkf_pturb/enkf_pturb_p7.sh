#!/bin/ksh
#

# ----------------------------------------
# Launcher script for enkf_pturb on IBM p7
# ----------------------------------------

#
# User-defined options
#
flnml="namelist_p7.nml"
machine="spica"
gest="/users/dor/arma/bue/power7/3dvar_modular/enkf_pturb_stag/"
anlgrid="/users/dor/arma/gr3/data_cnfs/prototype_grid/analysisgrid_glb_800x400"
cov="spica:/users/dor/arma/bue/power7/3dvar_modular/testnmc/glbstrato01_nmc_spstddev_80n_stag_t108_taper_uvt_fix_theta_with_toctoc_reformat.fst"
abs="/users/dor/arma/bue/home01/3dvar_latest/compiledir_enkf_pturb/enkf_pturb_p7.abs"
npex=2
npey=8
openmp=2
maxcputime=180

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
 cd $gest
 export TMG_ON=YES
 export MP_STDOUTMODE=unordered
 r.mpirun2 -pgm ./enkf_pturb.abs
EOF

cat << EOF > ptopo_nml
 &ptopo
  npex=$npex
  npey=$npey
/
EOF
scp ptopo_nml ${machine}:${gest}

ord_soumet go_enkf_pturb.sh -mach $machine -mpi -t $maxcputime -cm 1632M -q preemptable -cpus ${npex}x${npey}x${openmp} -listing ${gest} -jn enkf_pturb
