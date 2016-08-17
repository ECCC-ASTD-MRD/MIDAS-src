#!/bin/ksh
#

# ------------------------------------------------
# Launcher script for optimization tests on IBM p7
# ------------------------------------------------

#
# User-defined options
#
guest=$(pwd)
machine="hadar"
abs="/users/dor/arma/gr3/home1/var/test_git/compiledir_optimization/addEnsMember_p7.abs"
npex=1
npey=1
openmp=4
maxcputime=360
memory=3264M

#
# Don't modify below ...
#

abs_basename=`basename $abs`

echo
echo "Launching OPTIMIZATION TESTS using..."
echo
echo "Working machine      :" $machine
echo "Guest                :" $guest
echo "Executable file path :" $abs
echo "Executable file name :" $abs_basename
echo "Topology             :" ${npex}x${npey}x${openmp}
echo

cat << EOF > $TMPDIR/go_optimization.sh
 echo "!!STARTING SCRIPT!!"
 ulimit -a
 export TMG_ON=YES
 cd $guest
 ./addEnsMember_p7.abs
EOF

cat << EOF > $TMPDIR/ptopo_nml
 &ptopo
  npex=$npex
  npey=$npey
/
EOF
scp $TMPDIR/ptopo_nml ${machine}:${gest}

ord_soumet $TMPDIR/go_optimization.sh -mach $machine -t $maxcputime -cpus ${npex}x${npey}x${openmp} -cm ${memory} -listing ${guest} -jn optimization
