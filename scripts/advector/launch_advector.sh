#!/bin/bash
#

#
# User-defined options
#
machine=hare
abs="/home/jfc425/bin/midas/midas_abs/midas-advector_sles-11-broadwell-64-xc40-this_is_a_temporary_tag-5-g0d20ef5_M.Abs"
gest="${HOME}/data_maestro/${machine}/advector/test"
flnml="namelist.nml"
steeringFlow="/home/jfc425/tmp_big/ens_mean.fst"
fileToAdvect="/home/jfc425/resultats/advection/global/diagBmatrix/012m/4d-advecPertIncFromMiddle-obsFirstTime-012m/ens_pert1.fst"

npex=12
npey=4
openmp=3
maxcputime=3600
run_in_parallel=r.run_in_parallel

#
# Don't modify below ...
#
abs_basename=`basename $abs`

echo
echo "Launching ADVECTOR using..."
echo
echo "Using namelist file  :" $flnml
echo "Working machine      :" $machine
echo "Working directory    :" $gest
echo "Steering flow        :" $steeringFlow
echo "File to advect       :" $fileToAdvect
echo "Executable file path :" $abs
echo "Executable file name :" $abs_basename
echo "Topology             :" ${npex}x${npey}x${openmp}
echo

ssh $machine rm -rf $gest
ssh $machine mkdir -p $gest
scp $flnml ${machine}:${gest}/flnml
scp $steeringFlow ${machine}:${gest}/steeringFlow.fst
scp $fileToAdvect ${machine}:${gest}/fileToAdvect.fst
scp $abs ${machine}:${gest}/advector.abs
ssh $machine ls -l $gest

cat << EOF > $TMPDIR/go_advector.sh
#!/bin/bash
set -ex
 echo "!!STARTING SCRIPT!!"
 . ssmuse-sh -d eccc/mrd/rpn/utils/19.2
 
 # MPI SETUP FOR PPP ONLY
 if [ "\${EC_ARCH}" = ubuntu-18.04-skylake-64 ]; then
     . ssmuse-sh -d hpco/exp/intelpsxe-cluster-19.0.3.199
     . ssmuse-sh -d hpco/exp/openmpi/openmpi-3.1.2--hpcx-2.2.0--ofed-4.4.2--intel-2019.0.045
 fi

 cd $gest
 export TMG_ON=YES
 ${run_in_parallel} -pgm ./advector.abs -npex ${npex} -npey ${npey} -processorder -tag -nocleanup -verbose
EOF

cat << EOF > $TMPDIR/ptopo_nml
 &ptopo
  npex=$npex
  npey=$npey
/
EOF
scp $TMPDIR/ptopo_nml ${machine}:${gest}

ord_soumet $TMPDIR/go_advector.sh -mach $machine -mpi -t $maxcputime -cm 3000M -cpus ${npex}x${npey}x${openmp} -jn advector -waste
