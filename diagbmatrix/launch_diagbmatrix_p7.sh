#!/bin/ksh
#

# ------------------------------------------
# Launcher script for diag_bmatrix on IBM p7
# ------------------------------------------

#
# User-defined options
#
flnml="namelist_p7.nml"
machine="hadar"
gest="/users/dor/arma/gr3/data_gpfs/var/national_10km/test_scalingRandomPert/diagBmatrix/output"
bgcov=/users/dor/armn/gr4/power7/var/national_10km/calcstats_r640m/national_24-48_v664_h633_t250_CORR_chi_mass/bgcov_tapered.fst
analysisgrid=/users/dor/armn/gr4/power7/var/national_10km/calcstats_r640m/national_24-48_v664_h633_t250_CORR_chi_mass/analysisgrid
abs=/users/dor/arma/gr3/home1/var/latest_trunk/compiledir_diagbmatrix/diagbmatrix_p7.abs
npex=8
npey=4
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

ssh $machine rm -rf $gest
ssh $machine mkdir -p $gest
scp $flnml ${machine}:${gest}/flnml
scp $analysisgrid ${machine}:${gest}/analysisgrid
scp $bgcov ${machine}:${gest}/bgcov
scp $abs ${machine}:${gest}/diag_bmatrix.abs
ssh $machine ls -l $gest

cat << EOF > $TMPDIR/go_diagbmatrix.sh
 echo "!!STARTING SCRIPT!!"
. ssmuse-sh -d rpn/utils/15.2
 cd $gest
 export TMG_ON=YES
 export MP_STDOUTMODE=ordered
 r.run_in_parallel -pgm ./diag_bmatrix.abs -npex ${npex} -npey ${npey}
 ~armabue/bin/combineprocs.sh
 rm -f *_proc???
EOF

cat << EOF > $TMPDIR/ptopo_nml
 &ptopo
  npex=$npex
  npey=$npey
/
EOF
scp $TMPDIR/ptopo_nml ${machine}:${gest}

ord_soumet $TMPDIR/go_diagbmatrix.sh -mach $machine -t $maxcputime -cpus ${npex}x${npey}x${openmp} -cm ${memory} -listing ${gest} -jn diagbmatrix -waste
