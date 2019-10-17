#!/bin/bash
#

#
# User-defined options
#
expName="4d-advecAmplitude-fromMiddle"
machine=eccc-ppp2
abs="/home/jfc425/bin/midas/midas_abs/midas-adjointTest_ubuntu-14.04-amd64-64-m_3.1.0-75-gace547f_M.Abs"
gest="${HOME}/data_maestro/${machine}/adjointTest/${expName}/"
flnml="namelist.nml"
analysisgrid="/home/sanl000/ANAL_shared/datafiles/constants/arma/oavar/2.1.3/analysis_grid_prototypes/analysis_grid_prototype_glb_1080x540_south-to-north_80L_vcode5002"
bgcov="/home/sanl000/ANAL_shared/datafiles/constants/arma/oavar/2.1.4/__GEM25km_NMC_T399_stag5002_BgckStddev3d_corns_sqrt_800x400__/06"
ensdir="/fs/site2/dev/eccc/mrd/rpndat/jfc425/adjointTest/ensemble"
trials="/home/jfc425/maestro_archives/V61C001E16/gridpt.trial.hyb.2016061500_180m /home/jfc425/maestro_archives/V61C001E16/gridpt.trial.hyb.2016061500_360m /home/jfc425/maestro_archives/V61C001E16/gridpt.trial.hyb.2016061500_540m"
#trials="/home/jfc425/maestro_archives/V61C001E16/gridpt.trial.hyb.2016061500_360m"
npey=10
npex=6
openmp=4
maxcputime=1800
run_in_parallel="r.run_in_parallel"

#
# Don't modify below ...
#
abs_basename=`basename $abs`

echo
echo "Launching adjointTest using..."
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
scp $abs ${machine}:${gest}/adjointTest.abs
scp $bgcov ${machine}:${gest}/bgcov
ssh $machine ls -l $gest

count=1
for trial in $trials ; do
    ssh $machine ln -s $trial ${gest}/trlm_0${count}
    count=$(expr $count + 1)
done

cat << EOF > $TMPDIR/go_adjointTest.sh
#!/bin/bash
set -ex
 echo "!!STARTING SCRIPT!!"
 . ssmuse-sh -d eccc/mrd/rpn/utils/19.3
 
 # MPI SETUP FOR PPP ONLY
 if [ "\${EC_ARCH}" = ubuntu-18.04-skylake-64 ]; then
     . ssmuse-sh -d hpco/exp/intelpsxe-cluster-19.0.3.199
     . ssmuse-sh -d hpco/exp/openmpi/openmpi-3.1.2--hpcx-2.2.0--ofed-4.4.2--intel-2019.0.045
     . ssmuse-sh -d hpco/exp/openmpi-setup/openmpi-setup-0.2
 fi

 cd $gest
 export TMG_ON=YES
 ${run_in_parallel} -pgm ./adjointTest.abs -npex ${npex} -npey ${npey} -processorder -tag -nocleanup -verbose
EOF

cat << EOF > $TMPDIR/ptopo_nml
 &ptopo
  npex=$npex
  npey=$npey
/
EOF
scp $TMPDIR/ptopo_nml ${machine}:${gest}

ord_soumet $TMPDIR/go_adjointTest.sh -mach $machine -mpi -t $maxcputime -cm 3000M -cpus ${npex}x${npey}x${openmp} -jn adjointTest -waste
