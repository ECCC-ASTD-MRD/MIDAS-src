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
gest="/users/dor/arma/bue/power7/3dvar_modular/diag"
bgcov="hadar:/home/dormrb02/modeles/ANAL_shared/stats/auto/__GEM25km_NMC_T399_stag5002_BgckStddev3d_800x400__/01"
analysisgrid="/home/dormrb02/ibmenv/armnlib/modeles/ANAL_shared/datafiles/constants/arma/oavar/2.1.1/analysis_grid_prototypes/analysis_grid_prototype_glb_800x400_south-to-north"
abs="/users/dor/arma/bue/home01/3dvar_git/compiledir_diag_bmatrix/diag_bmatrix_p7.abs"
npex=4
npey=2
openmp=2
maxcputime=300
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
. ssmuse-sh -d rpn/utils/15.2
 cd $gest
 export TMG_ON=YES
 export MP_STDOUTMODE=ordered
 r.run_in_parallel -pgm ./diag_bmatrix.abs -npex ${npex} -npey ${npey}
EOF

cat << EOF > $TMPDIR/ptopo_nml
 &ptopo
  npex=$npex
  npey=$npey
/
EOF
scp $TMPDIR/ptopo_nml ${machine}:${gest}

ord_soumet $TMPDIR/go_diag_bmatrix.sh -mach $machine -t $maxcputime -cpus ${npex}x${npey}x${openmp} -cm ${memory} -listing ${gest} -jn diag_bmatrix -waste
