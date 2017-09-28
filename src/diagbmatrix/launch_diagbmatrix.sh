#!/bin/bash
#

#
# User-defined options
#
flnml="namelist.nml"
machine=brooks
abs="/home/jfc425/bin/var/oavar_abs/diagbmatrix_sles-11-broadwell-64-xc40-m_2.2.2-39-ge51de4a_M.Abs"
ensdir="/home/jfc425/data_maestro/brooks/maestro/oavar_system_tests_science/work/20170724000000/Tests/EnVar_24mem_xc40/UnitTest/envar/interpEnsTrials/output/"
gest="${HOME}/data_maestro/${machine}/diagbmatrix/test/"
analysisgrid="/home/sanl000/ANAL_shared/datafiles/constants/arma/oavar/2.1.2/analysis_grid_prototype_glb_1080x540_south-to-north_80L_vcode5002"
bgcov="/home/sanl000/ANAL_shared/datafiles/constants/arma/oavar/2.1.1/__GEM25km_NMC_T399_stag5002_BgckStddev3d_corns_sqrt_800x400__/01"
npex=12
npey=12
openmp=1
maxcputime=600
run_in_parallel="/fs/ssm/eccc/mrd/rpn/utils/16.2/all/bin/r.run_in_parallel_1.1.28c"

#
# Don't modify below ...
#
abs_basename=`basename $abs`

echo
echo "Launching DIAGBMATRIX using..."
echo
echo "Using namelist file  :" $flnml
echo "Working machine      :" $machine
echo "Working directory    :" $gest
echo "Analysis grid        :" $analysisgrid
echo "Executable file path :" $abs
echo "Executable file name :" $abs_basename
echo "Topology             :" ${npex}x${npey}x${openmp}
echo

ssh $machine rm -rf $gest
ssh $machine mkdir -p $gest
ssh $machine ln -s ${ensdir} ${gest}/ensemble
scp $flnml ${machine}:${gest}/flnml
scp $analysisgrid ${machine}:${gest}/analysisgrid
scp $bgcov ${machine}:${gest}/bgcov
scp $abs ${machine}:${gest}/diagbmatrix.abs
ssh $machine ls -l $gest

cat << EOF > $TMPDIR/go_diagbmatrix.sh
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
 ${run_in_parallel} -pgm ./diagbmatrix.abs -npex ${npex} -npey ${npey} -processorder -tag -nocleanup -verbose
EOF

cat << EOF > $TMPDIR/ptopo_nml
 &ptopo
  npex=$npex
  npey=$npey
/
EOF
scp $TMPDIR/ptopo_nml ${machine}:${gest}

ord_soumet $TMPDIR/go_diagbmatrix.sh -mach $machine -mpi -t $maxcputime -cm 3000M -cpus ${npex}x${npey}x${openmp} -jn diagbmatrix -waste
