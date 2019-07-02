#!/bin/bash
#

#
# User-defined options
#
machine=hare
ensbase="G-NMC"
abs="/home/mab001/data_maestro/ords/midas_abs/midas-diagBmatrix_sles-11-broadwell-64-xc40-v_3.0.4-34-g3951b62_M.Abs"
gest="${HOME}/data_maestro/${machine}/diagbmatrix/${ensbase}_test_BlatBands_enkf_new/"
if [ "${ensbase}" = "G-NMC" ]; then
  flnml="namelist.nml_${ensbase}"
  #analysisgrid="/home/sanl000/ANAL_shared/datafiles/constants/arma/oavar/2.1.1/analysis_grid_prototypes/analysis_grid_prototype_glb_800x400_south-to-north_80L_vcode5002"
  analysisgrid="/home/sanl000/ANAL_shared/datafiles/constants/arma/oavar/2.1.3/analysis_grid_prototypes/analysis_grid_prototype_glb_1080x540_south-to-north_80L_vcode5002"
  #bgcov="/home/sanl000/ANAL_shared/datafiles/constants/arma/oavar/2.1.1/__GEM25km_NMC_T399_stag5002_BgckStddev3d_corns_sqrt_800x400__/01"
  #bgcov="/home/mab001/data_maestro/eccc-ppp2/calcstats/test_5002_latbands_m_3.0.0/latbands_test2_1080x540.fst"
  bgcov="/home/mab001/data_maestro/eccc-ppp2/calcstats/test_5002_latbands_m_3.0.0_enkf/latbands_2017011006_1080x540.fst"
  ensdir="/home/bed666/brooks/diagbmatrix/ensemble/G-EnKF/"
  npey=9
  npex=9
else
  if [ "${ensbase}" = "N-NMC" ]; then
    flnml="namelist.nml_${ensbase}"
    ensdir="/home/bed666/brooks/diagbmatrix/ensemble/G-EnKF/"
  else
    flnml="namelist.nml_ens"
    ensdir="/home/bed666/brooks/diagbmatrix/ensemble/${ensbase}/"
  fi
  analysisgrid="/home/jfc425/data/ords/oavarGridTemplate/analysisgrid_national10km_80L_vcode5002.fst"
  bgcov="/home/bed666/ss2/cetus/hadar/national_10km/calcstats_r640m/national_24-48_v664_h633_t250_CORR_LQ/07"
  npey=12
  npex=12
fi
openmp=1
maxcputime=1200
run_in_parallel=r.run_in_parallel

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
#ssh $machine ln -s ${ensdir} ${gest}/ensemble
scp $flnml ${machine}:${gest}/flnml
scp $analysisgrid ${machine}:${gest}/analysisgrid
scp $abs ${machine}:${gest}/diagbmatrix.abs
scp $bgcov ${machine}:${gest}/bgcov
ssh $machine ls -l $gest

cat << EOF > $TMPDIR/go_diagbmatrix.sh
#!/bin/bash
set -ex
 echo "!!STARTING SCRIPT!!"
 . ssmuse-sh -d eccc/mrd/rpn/utils/19.0
 
 # MPI SETUP FOR PPP ONLY
 if [ "\${EC_ARCH}" = ubuntu-18.04-skylake-64 ]; then
     . ssmuse-sh -d hpco/exp/intelpsxe-cluster-19.0.3.199
     . ssmuse-sh -d hpco/exp/openmpi/openmpi-3.1.2--hpcx-2.2.0--ofed-4.4.2--intel-2019.0.045
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

ord_soumet $TMPDIR/go_diagbmatrix.sh -mach $machine -mpi -t $maxcputime -cm 3000M -cpus ${npex}x${npey}x${openmp} -jn diagBmatrix -waste
