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
gest="/users/dor/arma/gr3/data_gpfs/comparaison_var/new_var/enkf_pturb/lam/atelier"
anlgrid="/users/dor/arma/gr3/data_gpfs/comparaison_var/new_var/bnmc_step/lam/analysisgrid"
cov="/users/dor/arma/gr3/data_gpfs/comparaison_var/new_var/bnmc_step/lam/120m_NormByStdDev_rev388mod_bugfix2/bgcov.fst"
abs="/users/dor/arma/gr3/home1/new_var/trunk_392_mod/compiledir_enkf_pturb/enkf_pturb_p7.abs"
npex=1
npey=8
openmp=4
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

ssh $machine rm -f $gest/*
scp $flnml ${machine}:${gest}/flnml
scp $anlgrid ${machine}:${gest}/analysisgrid
scp $cov ${machine}:${gest}/bgcov
scp $abs ${machine}:${gest}/enkf_pturb.abs
ssh $machine ls -l $gest

cat << EOF > go_enkf_pturb.sh
 echo "!!STARTING SCRIPT!!"
 cd $gest
 ./enkf_pturb.abs
EOF

cat << EOF > ptopo_nml
 &ptopo
  npex=$npex
  npey=$npey
/
EOF
scp ptopo_nml ${machine}:${gest}

ord_soumet go_enkf_pturb.sh -mach $machine -mpi -t $maxcputime -cpus ${npex}x${npey}x${openmp} -listing ${gest} -jn enkf_pturb
