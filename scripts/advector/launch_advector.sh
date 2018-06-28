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
run_in_parallel="/fs/ssm/eccc/mrd/rpn/utils/16.2/all/bin/r.run_in_parallel_1.1.28c"

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
. ssmuse-sh -d eccc/mrd/rpn/utils/16.2

# MPI SETUP FOR PPP ONLY
if [ "\${EC_ARCH}" = Linux_x86-64 ]; then
  . ssmuse-sh -d hpco/tmp/eccc/201402/05/base
  . ssmuse-sh -d main/opt/intelcomp/intelcomp-2016.1.156
  . ssmuse-sh -d main/opt/openmpi/openmpi-1.6.5/intelcomp-2016.1.156 
  export OMPI_MCA_orte_tmpdir_base=/run/shm
  export OMPI_MCA_btl_openib_if_include=mlx5_0
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
