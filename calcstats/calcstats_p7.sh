#!/bin/ksh
#

# ---------------------------------------
# Launcher script for calcstats on IBM p7
# ---------------------------------------

#
# User-defined options
#
flnml="namelist_glb_p7.nml"
machine="hadar"
gest="/users/dor/arma/gr3/data_gpfs/var/gonzalo/calcstats_hvcorrel/atelier"
ensdir="/users/dor/arma/gr3/data_gpfs/var/gonzalo/ensemble/interpEnsTrials/gaussian_grid"
abs="/users/dor/arma/gr3/home1/var/trunk_572m/compiledir_calcstats/calcstats_p7.abs"
npex=1
npey=1
openmp=32
maxcputime=500
memory=3264M

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

cat << EOF > $TMPDIR/go_calcstats.sh
 echo "!!STARTING SCRIPT!!"
 ulimit -a
 cd $gest
 export TMG_ON=YES
 ./calcb.abs
EOF

cat << EOF > $TMPDIR/ptopo_nml
 &ptopo
  npex=$npex
  npey=$npey
/
EOF
scp $TMPDIR/ptopo_nml ${machine}:${gest}

ord_soumet $TMPDIR/go_calcstats.sh -mach $machine -t $maxcputime -cpus ${npex}x${npey}x${openmp} -cm ${memory} -listing ${gest} -jn calcstats
