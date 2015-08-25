#!/bin/ksh
#

# ------------------------------------------
# Launcher script for diag_bmatrix on IBM p7
# ------------------------------------------

#
# User-defined options
#
flnml="namelist_glb_p7.nml"
machine="hadar"
gest="/users/dor/arma/gr3/data_gpfs/var/gonzalo/diag_bmatrix/atelier"
bgcov="cassini:/users/dor/arma/gr3/local/var/NewCovFormat4bgck/__GLBSTRATO_NMC_T380_TAPER_UVT_FIX_THETA__newFormat800x400/glbstrato10_nmc_80n_hyb_t380_taper_uvt_fix_theta_with_toctoc"
analysisgrid="datasvr:/users/dor/arma/gr3/data_cnfs/prototype_grid/analysisgrid_glb_800x400"
abs="/users/dor/arma/gr3/home1/var/trunk_572m/compiledir_diag_bmatrix/diag_bmatrix_p7.abs"
npex=2
npey=2
openmp=8
maxcputime=1800
memory=3264M

#
# Don't modify below ...
#
abs_basename=`basename $abs`

echo
echo "Launching DIAG_BMATRIX using..."
echo
echo "Using namelist file  :" $flnml
echo "Analysis grid        :" $analysisgrid
echo "Static B matrix      :" $bgcov
echo "Working machine      :" $machine
echo "Working directory    :" $gest
echo "Executable file path :" $abs
echo "Executable file name :" $abs_basename
echo "Topology             :" ${npex}x${npey}x${openmp}
echo

ssh $machine rm -f $gest/*
scp $flnml ${machine}:${gest}/flnml
scp $analysisgrid ${machine}:${gest}/analysisgrid
scp $bgcov ${machine}:${gest}/bgcov
scp $abs ${machine}:${gest}/diag_bmatrix.abs
ssh $machine ls -l $gest

cat << EOF > $TMPDIR/go_diag_bmatrix.sh
 echo "!!STARTING SCRIPT!!"
 ulimit -a
 cd $gest
 export TMG_ON=YES
 ./diag_bmatrix.abs
EOF

cat << EOF > $TMPDIR/ptopo_nml
 &ptopo
  npex=$npex
  npey=$npey
/
EOF
scp $TMPDIR/ptopo_nml ${machine}:${gest}

ord_soumet $TMPDIR/go_diag_bmatrix.sh -mach $machine -t $maxcputime -cpus ${npex}x${npey}x${openmp} -cm ${memory} -listing ${gest} -jn diag_bmatrix
