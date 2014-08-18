#!/bin/ksh
#

flnml="namelist_ens.nml"
machine="hadar"
gest="/users/dor/arma/bue/power7/3dvar_modular/calcb"
ensdir="/users/dor/arma/bue/power7/3dvar_modular/enkf_2011022200/"
abs="../../compiledir/calcb_p7.abs_NOMPI"

abs_basename=`basename $abs`

echo "Using namelist file:  " $flnml
echo "Working machine:      " $machine
echo "Working directory:    " $gest
echo "Ensemble directory:   " $ensdir
echo "Executable file path: " $abs
echo "Executable file name: " $abs_basename

scp $flnml ${machine}:${gest}/flnml
scp $abs ${machine}:${gest}/calcb.abs
ssh $machine ln -s ${ensdir} ${gest}/ensemble
ssh $machine ls -l $gest

echo "echo !!STARTING SCRIPT!!" > go_calcstats.sh
echo "cd $gest " >> go_calcstats.sh
echo "rm -f *.fst " >> go_calcstats.sh
echo "./calcb.abs " >> go_calcstats.sh

ord_soumet go_calcstats.sh -mach $machine -mpi -t 7200 -cpus 1x1x32 -listing ${gest} -jn calcb
