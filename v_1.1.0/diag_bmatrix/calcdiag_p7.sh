#!/bin/ksh
#

npex=60
npey=1
nopenmp=2

#flnml="namelist.nml_bhybrid_h9"
flnml="namelist.nml_bhi"
machine="spica"
gest="/users/dor/arma/bue/power7/3dvar_modular/diag"
abs="../../compiledir/diag_bmatrix_p7.abs"
retourmach="datasvr"
#retourdir="/users/dor/arma/bue/cnfs/diagnostics/pert100_600x300_bhybrid_h9/"
retourdir="/users/dor/arma/bue/cnfs/diagnostics/pert100_600x300_bhi/"

abs_basename=`basename $abs`

echo "Using namelist file:  " $flnml
echo "Working machine:      " $machine
echo "Working directory:    " $gest
echo "Executable file path: " $abs
echo "Executable file name: " $abs_basename

echo " &PTOPO " > ptopo_nml
echo "   npex = $npex " >> ptopo_nml
echo "   npey = $npey " >> ptopo_nml
echo "/ " >> ptopo_nml

ssh $retourmach mkdir -p $retourdir
ssh $machine rm -f ${gest}/*
ssh $machine ln -s /home/dormrb02/ibmenv/armnlib/modeles/ANAL/stats/auto/glbstrato02_nmc_80n_hyb_t380_taper_fix_theta_with_toctoc ${gest}/glbcov
scp $flnml ${machine}:${gest}/flnml
scp ptopo_nml ${machine}:${gest}/ptopo_nml
scp $abs ${machine}:${gest}/diag_bmatrix_p7.abs_NOMPI
ssh $machine ls -l $gest

echo "echo !!STARTING SCRIPT!!" > go_calcdiag_p7.sh
echo "export MP_STDOUTMODE=unordered" >> go_calcdiag_p7.sh
echo "cd $gest " >> go_calcdiag_p7.sh
echo "rm -f *.fst " >> go_calcdiag_p7.sh
echo "echo STARTING EXECUTABLE AT:" >> go_calcdiag_p7.sh
echo "date" >> go_calcdiag_p7.sh
echo "./diag_bmatrix_p7.abs_NOMPI " >> go_calcdiag_p7.sh
echo "echo FINISHED EXECUTABLE AT:" >> go_calcdiag_p7.sh
echo "date" >> go_calcdiag_p7.sh
echo "scp *.fst ${retourmach}:${retourdir}/" >> go_calcdiag_p7.sh
echo "scp *.out ${retourmach}:${retourdir}/" >> go_calcdiag_p7.sh

ord_soumet go_calcdiag_p7.sh -mach $machine -mpi -t 1800 -cpus ${npex}x${npey}x${nopenmp} -listing ${gest} -jn calcdiag -waste
