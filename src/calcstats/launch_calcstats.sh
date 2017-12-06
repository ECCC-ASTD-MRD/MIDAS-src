#!/bin/bash
#

#
# User-defined options
#
#flnml="namelist_hvloc_lam.nml"
#ensdir=/home/jfc425/data_maestro/brooks/ensemble/national_test #global_test/gaussian_grid
#analysisgrid="/home/jfc425/data/ords/oavarGridTemplate/analysisgrid_national10km_80L_vcode5002.fst"

flnml="namelist_bhi_glb.nml"
machine=brooks
abs=/home/mab001/data_maestro/ords/oavar_abs/calcstats_sles-11-broadwell-64-xc40-m_3.0.0-1-gaea63f4_M.Abs
expname="test_5002_latbands_m_3.0.0_enkf"
ensdir=${HOME}/data_maestro/${machine}/calcstats/ensemble_enkf
gest="${HOME}/data_maestro/${machine}/calcstats/${expname}"
analysisgrid="/home/sanl000/ANAL_shared/datafiles/constants/arma/oavar/2.1.1/analysis_grid_prototypes/analysis_grid_prototype_glb_800x400_south-to-north_80L_vcode5002"

npex=1
npey=1
openmp=36 #44
memory=3000M #220000M
maxcputime=10800

run_in_parallel="/fs/ssm/eccc/mrd/rpn/utils/16.2/all/bin/r.run_in_parallel_1.1.28c"

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
ssh $machine ln -s ${ensdir} ${gest}/ensemble
scp $flnml ${machine}:${gest}/flnml
scp $analysisgrid ${machine}:${gest}/analysisgrid
scp $abs ${machine}:${gest}/calcstats.abs
ssh $machine ls -l $gest

cat << EOF > $TMPDIR/go_calcstats.sh
#!/bin/bash
set -ex
 echo "!!STARTING SCRIPT!!"
. ssmuse-sh -d eccc/mrd/rpn/utils/16.2

export OMP_STACKSIZE=4096M

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
 ${run_in_parallel} -pgm ./calcstats.abs -npex ${npex} -npey ${npey} -processorder -tag -nocleanup -verbose
 export MP_STDOUTMODE=ordered
EOF

cat << EOF > $TMPDIR/ptopo_nml
 &ptopo
  npex=$npex
  npey=$npey
/
EOF
scp $TMPDIR/ptopo_nml ${machine}:${gest}

ord_soumet $TMPDIR/go_calcstats.sh -mach $machine -t $maxcputime -m $memory -cpus ${npex}x${npey}x${openmp} -jn calcstats -waste
