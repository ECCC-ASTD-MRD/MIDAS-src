#!/bin/ksh
#

# ---------------------------------------
# Launcher script for calcstats on IBM p7
# ---------------------------------------

#
# User-defined options
#
flnml="namelist_lam_p7.nml"
machine="spica"
gest="/users/dor/arma/gr3/data_gpfs/comparaison_var/new_var/bnmc_step/lam/atelier"
ensdir="/users/dor/arma/gr3/data_gpfs/comparaison_var/new_var/bnmc_step/lam/training_data/"
abs="/users/dor/arma/gr3/home1/new_var/branches/varlam/compiledir_calcstats/calcstats_p7.abs_NOMPI"
npex=1
npey=1
openmp=32
maxcputime=180

#
# Don't modify below ...
#
abs_basename=`basename $abs`

echo
echo "Launching CALCSTATS using..."
echo
echo "Using namelist file  :" $flnml
echo "Working machine      :" $machine
echo "Working directory    :" $gest
echo "Ensemble directory   :" $ensdir
echo "Executable file path :" $abs
echo "Executable file name :" $abs_basename
echo "Topology             :" ${npex}x${npey}x${openmp}
echo

ssh $machine rm -f $gest/*
scp $flnml ${machine}:${gest}/flnml
scp $abs ${machine}:${gest}/calcb.abs
ssh $machine ln -s ${ensdir} ${gest}/ensemble
ssh $machine ls -l $gest

cat << EOF > go_calcstats.sh
 echo "!!STARTING SCRIPT!!"
 cd $gest
 ./calcb.abs
EOF

cat << EOF > ptopo_nml
 &ptopo
  npex=$npex
  npey=$npey
/
EOF
scp ptopo_nml ${machine}:${gest}

ord_soumet go_calcstats.sh -mach $machine -mpi -t $maxcputime -cpus ${npex}x${npey}x${openmp} -listing ${gest} -jn calcb
