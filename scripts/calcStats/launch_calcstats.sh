#!/bin/bash
#

#
# User-defined options
#
flnml="namelist_hvloc_lam.nml"
ensdir=/home/jfc425/data_maestro/eccc-ppp2/ensemble/national_test #global_test/gaussian_grid
analysisgrid="/home/jfc425/data/ords/oavarGridTemplate/analysisgrid_national10km_80L_vcode5002.fst"

#flnml="namelist_bhi_glb.nml"
machine=eccc-ppp2
abs=/home/jfc425/bin/midas/midas_abs/midas-calcStats_ubuntu-14.04-amd64-64-m_3.1.0_M.Abs
expname="test_hvlocal"
#ensdir=${HOME}/data_maestro/${machine}/calcstats/ensemble_enkf
gest="${HOME}/data_maestro/${machine}/calcstats/${expname}"
#analysisgrid="/home/sanl000/ANAL_shared/datafiles/constants/arma/oavar/2.1.1/analysis_grid_prototypes/analysis_grid_prototype_glb_800x400_south-to-north_80L_vcode5002"

npex=1
npey=1
openmp=40 #44
memory=100000M #220000M
maxcputime=10800

run_in_parallel=r.run_in_parallel

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
 if [ "${TRUE_HOST}" = eccc-ppp1 -o "${TRUE_HOST}" = eccc-ppp2 -o "${TRUE_HOST}" = hare -o "${TRUE_HOST}" = brooks ]; then
   . ssmuse-sh -d eccc/mrd/rpn/utils/16.2.3
 elif [ "${TRUE_HOST}" = eccc-ppp3 -o "${TRUE_HOST}" = eccc-ppp4 -o "${TRUE_HOST}" = daley -o "${TRUE_HOST}" = banting -o "${TRUE_HOST}" = xc3 -o "${TRUE_HOST}" = xc4 ]; then
   . ssmuse-sh -d eccc/mrd/rpn/utils/19.2
 else
   echo "Unknown TRUE_HOST: ${TRUE_HOST}"
   exit
 fi

 # MPI SETUP FOR PPP ONLY
 if [ "\${EC_ARCH}" = Linux_x86-64 ]; then
   . ssmuse-sh -d hpco/tmp/eccc/201402/05/base
   . ssmuse-sh -d main/opt/intelcomp/intelcomp-2016.1.156
   . ssmuse-sh -d main/opt/openmpi/openmpi-1.6.5/intelcomp-2016.1.156 
   export OMPI_MCA_orte_tmpdir_base=/run/shm
   export OMPI_MCA_btl_openib_if_include=mlx5_0
 fi
 
 # MPI SETUP FOR PPP ONLY
 if [ "\${EC_ARCH}" = ubuntu-18.04-skylake-64 ]; then
   . ssmuse-sh -d hpco/exp/intelpsxe-cluster-19.0.3.199
   . ssmuse-sh -d hpco/exp/openmpi/openmpi-3.1.2--hpcx-2.2.0--ofed-4.4.2--intel-2019.0.045
   . ssmuse-sh -d hpco/exp/openmpi-setup/openmpi-setup-0.2
 fi

 export OMP_STACKSIZE=4096M

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
