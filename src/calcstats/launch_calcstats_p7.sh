#!/bin/ksh
#

# ---------------------------------------
# Launcher script for calcstats on IBM p7
# ---------------------------------------

#
# User-defined options
#
flnml="namelist_lam_p7.nml"
machine="hadar"
gest="/users/dor/arma/gr3/data_gpfs/var/national_10km/test_scalingRandomPert/stddev_ens_raw/output"
ensdir="/users/dor/arma/gr3/data_gpfs/maestro/UnitTests/EnVar_LAM/work/20140702000000/Tests/EnVar_LAM/UnitTest/envar/interpEnsTrials/output"
abs="/users/dor/arma/gr3/home1/var/latest_trunk/compiledir_calcstats/calcstats_p7.abs"
npex=1
npey=1
openmp=32
maxcputime=600
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

ssh $machine rm -rf $gest
ssh $machine mkdir -p $gest  
scp $flnml ${machine}:${gest}/flnml
scp $abs ${machine}:${gest}/calcstats.abs
ssh $machine ln -s ${ensdir} ${gest}/ensemble
ssh $machine ls -l $gest

cat << EOF > $TMPDIR/go_calcstats.sh
 echo "!!STARTING SCRIPT!!"
. ssmuse-sh -d rpn/utils/15.2
 ulimit -a
 cd $gest
 export TMG_ON=YES
 export MP_STDOUTMODE=ordered
 r.run_in_parallel -pgm ./calcstats.abs -npex ${npex} -npey ${npey}
EOF

cat << EOF > $TMPDIR/ptopo_nml
 &ptopo
  npex=$npex
  npey=$npey
/
EOF
scp $TMPDIR/ptopo_nml ${machine}:${gest}

ord_soumet $TMPDIR/go_calcstats.sh -mach $machine -t $maxcputime -cpus ${npex}x${npey}x${openmp} -cm ${memory} -listing ${gest} -jn calcstats
