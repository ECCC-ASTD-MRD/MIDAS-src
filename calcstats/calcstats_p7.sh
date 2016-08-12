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
gest="/users/dor/arma/bue/power7/test_calcstats/feeaf54_M"
ensdir="/users/dor/arma/bue/power7/nmc800x400/k5_redo"
#ensdir="/users/dor/arma/gr3/data_gpfs/var/gonzalo/ensemble/interpEnsTrials/gaussian_grid"
abs="/users/dor/arma/bue/home01/3dvar_git/compiledir_calcstats/calcstats_p7.abs"
npex=1
npey=1
openmp=32
maxcputime=10800
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
 export LDR_CNTRL=TEXTPSIZE=64K@STACKPSIZE=64K@DATAPSIZE=64K@SHMPSIZE=64K@MAXDATA64=0x1900000000
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
