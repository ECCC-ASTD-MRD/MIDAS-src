#!/bin/bash
#

#
# User-defined options
#
flnml="namelist.nml"
machine=hare
#machine=eccc-ppp2
abs="${HOME}/data_maestro/ords/oavar_abs/randompert_sles-11-broadwell-64-xc40-m_2.2.1-6-g15400d5_M.Abs"
#abs="${HOME}/data_maestro/ords/oavar_abs/randompert_ubuntu-14.04-amd64-64-m_2.2.1-6-g15400d5_M.Abs"
gest="${HOME}/data_maestro/${machine}/randompert/test/"
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
echo "Launching RANDOMPERT using..."
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
scp $flnml ${machine}:${gest}/flnml
scp $analysisgrid ${machine}:${gest}/analysisgrid
scp $bgcov ${machine}:${gest}/bgcov
scp $abs ${machine}:${gest}/randompert.abs
ssh $machine ls -l $gest

cat << EOF > $TMPDIR/go_randompert.sh
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
 ${run_in_parallel} -pgm ./randompert.abs -npex ${npex} -npey ${npey} -processorder -tag -nocleanup -verbose
EOF

cat << EOF > $TMPDIR/ptopo_nml
 &ptopo
  npex=$npex
  npey=$npey
/
EOF
scp $TMPDIR/ptopo_nml ${machine}:${gest}

ord_soumet $TMPDIR/go_randompert.sh -mach $machine -mpi -t $maxcputime -cm 3000M -cpus ${npex}x${npey}x${openmp} -jn randompert -waste
