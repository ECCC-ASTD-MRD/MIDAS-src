#!/bin/bash
#
set -ex

#
# User-defined options
#
#machine=hare
machine=eccc-ppp2
abs=/home/alc001/midas_issue_89/compiledir/midas_abs/midas-seaIce_ubuntu-14.04-amd64-64-m_3.1.0-99-gc5af702_M.Abs
#abs="${HOME}/data_maestro/ords/compiledir/compiledir-seaIce-ubuntu-14.04-amd64-64_m_3.1.0-97-g3f891fc_M/midas-seaIce_ubuntu-14.04-amd64-64-m_3.1.0-97-g3f891fc_M.Abs"
npex=1
npey=1
openmp=40
maxcputime=3000
run_in_parallel=r.run_in_parallel
#analysisgrid="/home/sanl000/ANAL_shared/datafiles/constants/arma/oavar/2.1.2/analysis_grid_prototype_glb_1080x540_south-to-north_80L_vcode5002"
#analysisgrid="/home/mab001/constants/analysis_grid_prototype_glb_1080x540_south-to-north_sfc_vcode5002"
analysisgrid="/home/alc001/constants/cmde/icea/v2.2.5/arctic/analgrid.std"
implicit=T
if [ ${implicit} = F ] ; then
    diff_norm_fact=/home/alc001/constants/cmde/icea/v2.2.5/arctic/diffusmod_explicit.std
elif [ ${implicit} = T ] ; then
    diff_norm_fact=/home/alc001/constants/cmde/icea/v2.2.5/arctic/diffusmod_implicit.std
else
    echo "Variable implicit has to be either T or F."
    exit 1
fi

#
# Don't modify below ...
#

gest="${HOME}/data_maestro/${machine}/seaice/test"

# build the namelist
cat << EOF > $TMPDIR/flnml
 &NAMCT0
/
 &NAMTESTSEAICE
  DATE     = 2017010100
  ONEOBS_LONS = 401,
  ONEOBS_LATS = 200,250,300,350
/
 &NAMTIME
  DSTEPOBSINC = 6.0d0
/
 &NAMSTATE
  ANLVAR(1) = 'GL'
/
 &NAMBDIFF
  corr_len = 200*10.0
  stab     = 200*0.2
  nsamp    = 200*10000
  limplicit = 200*${implicit}
  SCALEFACTOR = 200*1.0
  stddevMode = 'HOMO'
  homogeneous_std = 200*0.1
/
 &NAMBHI
  SCALEFACTOR = 200*0.0
/
 &NAMBEN
  SCALEFACTOR = 200*0.0
/
EOF

abs_basename=`basename $abs`

echo
echo "Launching SEAICE using..."
echo
echo "Working machine      :" $machine
echo "Working directory    :" $gest
echo "Executable file path :" $abs
echo "Executable file name :" $abs_basename
echo "Topology             :" ${npex}x${npey}x${openmp}
echo

ssh $machine rm -rf $gest
ssh $machine mkdir -p $gest
scp $TMPDIR/flnml ${machine}:${gest}/flnml
scp ${analysisgrid} ${machine}:${gest}/analysisgrid
if [ -e ${diff_norm_fact} ] ; then
    scp  ${diff_norm_fact} ${machine}:${gest}/diffusmod.std
fi
scp $abs ${machine}:${gest}/seaice.abs
ssh $machine ls -l $gest

cat << EOF > $TMPDIR/go_seaice.sh
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

 cd $gest
 export TMG_ON=YES
 ${run_in_parallel} -pgm ./seaice.abs -npex ${npex} -npey ${npey} -processorder -tag -nocleanup -verbose
EOF

cat << EOF > $TMPDIR/ptopo_nml
 &ptopo
  npex=$npex
  npey=$npey
/
EOF
scp $TMPDIR/ptopo_nml ${machine}:${gest}/

ord_soumet $TMPDIR/go_seaice.sh -mach $machine -mpi -t $maxcputime -cm 3000M -cpus ${npex}x${npey}x${openmp} -jn seaice -waste
