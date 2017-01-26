#!/bin/ksh
#

# ----------------------------------------
# Launcher script for enkf_pturb on IBM p7
# ----------------------------------------

#
# User-defined options
#
flnml="namelist_p7.nml"
machine="hadar"
gest="/users/dor/arma/gr3/data_gpfs/var/national_10km/enkf_pturb/test/"
analysisgrid="/users/dor/armn/gr4/power7/var/national_10km/calcstats_r640m/national_24-48_v664_h633_t250_CORR_chi_mass/analysisgrid"
bgcov="/users/dor/armn/gr4/power7/var/national_10km/calcstats_r640m/national_24-48_v664_h633_t250_CORR_chi_mass/bgcov_tapered.fst"
abs="/users/dor/arma/gr3/home1/var/latest_trunk/compiledir_enkf_pturb/enkf_pturb_p7.abs"
npex=8
npey=2
openmp=2
maxcputime=300

#
# Don't modify below ...
#
abs_basename=`basename $abs`

echo
echo "Launching ENKF_PTURB using..."
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
scp $abs ${machine}:${gest}/enkf_pturb.abs
ssh $machine ls -l $gest

cat << EOF > $TMPDIR/go_enkf_pturb.sh
 echo "!!STARTING SCRIPT!!"
. ssmuse-sh -d rpn/utils/15.2
 cd $gest
 export TMG_ON=YES
 export MP_STDOUTMODE=ordered
 r.run_in_parallel -pgm ./enkf_pturb.abs -npex ${npex} -npey ${npey}
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

ord_soumet $TMPDIR/go_enkf_pturb.sh -mach $machine -mpi -t $maxcputime -cm 1632M -cpus ${npex}x${npey}x${openmp} -listing ${gest} -jn enkf_pturb -waste
