#!/bin/ksh
#

npex=60
npey=1
nopenmp=2

#flnml="namelist.nml_bhybrid_bndy3p0"
#flnml="namelist.nml_bhybrid_bndy0p9"
#flnml="namelist.nml_benkf"
#flnml="namelist.nml_bhybrid_fh"
#flnml="namelist.nml_bhi"
#flnml="namelist.nml_bhybrid50_2"
#flnml="namelist.nml_bhybrid50_3"
#flnml="namelist.nml_benkf50_2"
#flnml="namelist.nml_benkf50_3"
#flnml="namelist.nml_benkf50_2_testlatband"
#flnml="namelist.nml_benkf50_2_testlatband_ctrl"
flnml="namelist.nml_benkf50_6"

#retourdir="/users/dor/arma/bue/cnfs/diagnostics/pert100_600x300_bhybrid_bndy3p0/"
#retourdir="/users/dor/arma/bue/cnfs/diagnostics/pert100_600x300_bhybrid_bndy0p9/"
#retourdir="/users/dor/arma/bue/cnfs/diagnostics/pert100_600x300_benkf/"
#retourdir="/users/dor/arma/bue/cnfs/diagnostics/pert100_600x300_bhybrid_fh/"
#retourdir="/users/dor/arma/bue/cnfs/diagnostics/pert100_600x300_bhi/"
#retourdir="/users/dor/arma/bue/cnfs/diagnostics/pert100_600x300_bhybrid50_2/"
#retourdir="/users/dor/arma/bue/cnfs/diagnostics/pert100_600x300_bhybrid50_3/"
#retourdir="/users/dor/arma/bue/cnfs/diagnostics/pert100_600x300_benkf50_2/"
#retourdir="/users/dor/arma/bue/cnfs/diagnostics/pert100_600x300_benkf50_3/"
#retourdir="/users/dor/arma/bue/cnfs/diagnostics/pert100_800x400_benkf50_2_testlatband/"
#retourdir="/users/dor/arma/bue/cnfs/diagnostics/pert100_800x400_benkf50_2_testlatband_ctrl/"
retourdir="/users/dor/arma/bue/cnfs/diagnostics/pert100_600x300_benkf50_6/"

machine="spica"
gest="/users/dor/arma/bue/power7/3dvar_modular/diag"
abs="../../compiledir/diag_bmatrix_p7.abs"
retourmach="datasvr"

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
ssh $machine mkdir -p $gest
ssh $machine rm -f ${gest}/*
ssh $machine ln -s /home/dormrb02/ibmenv/armnlib/modeles/ANAL/stats/auto/glbstrato02_nmc_80n_hyb_t380_taper_fix_theta_with_toctoc ${gest}/glbcov
scp $flnml ${machine}:${gest}/flnml
scp ptopo_nml ${machine}:${gest}/ptopo_nml
scp $abs ${machine}:${gest}/diag_bmatrix_p7.abs
ssh $machine ls -l $gest

echo "echo !!STARTING SCRIPT!!" > go_calcdiag_p7.sh
echo "export MP_STDOUTMODE=unordered" >> go_calcdiag_p7.sh
echo "export TMG_ON=YES " >> go_calcdiag_p7.sh
echo "cd $gest " >> go_calcdiag_p7.sh
echo "rm -f *.fst " >> go_calcdiag_p7.sh
echo "echo STARTING EXECUTABLE AT:" >> go_calcdiag_p7.sh
echo "date" >> go_calcdiag_p7.sh
echo "./diag_bmatrix_p7.abs " >> go_calcdiag_p7.sh
echo "echo FINISHED EXECUTABLE AT:" >> go_calcdiag_p7.sh
echo "date" >> go_calcdiag_p7.sh
echo "srcp *.fst ${retourmach}:${retourdir}/" >> go_calcdiag_p7.sh
echo "srcp *.out ${retourmach}:${retourdir}/" >> go_calcdiag_p7.sh
echo "srcp TMG* ${retourmach}:${retourdir}/" >> go_calcdiag_p7.sh

ord_soumet go_calcdiag_p7.sh -mach $machine -mpi -t 3600 -cpus ${npex}x${npey}x${nopenmp} -listing ${gest} -jn calcdiag -waste
