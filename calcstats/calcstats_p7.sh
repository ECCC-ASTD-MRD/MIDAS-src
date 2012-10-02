#!/bin/ksh
#

flnml="namelist_p7.nml"
machine="spica"
gest="/users/dor/arma/bue/power7/3dvar_modular/calcb"
abs="../../compiledir/calcb_p7.abs"

abs_basename=`basename $abs`

echo "Using namelist file:  " $flnml
echo "Working machine:      " $machine
echo "Working directory:    " $gest
echo "Executable file path: " $abs
echo "Executable file name: " $abs_basename

scp $flnml ${machine}:${gest}/flnml
scp $abs ${machine}:${gest}/calcb.abs
ssh $machine ls -l $gest

echo "echo !!STARTING SCRIPT!!" > calcb.sh
echo "cd $gest " >> calcb.sh
echo "rm -f *.fst " >> calcb.sh
echo "./calcb.abs " >> calcb.sh

ord_soumet calcb.sh -mach $machine -mpi -t 7200 -cpus 1x1x32 -listing ${gest} -jn calcb
