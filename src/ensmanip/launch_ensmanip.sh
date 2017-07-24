#!/bin/bash
#

#
# User-defined options
#
#machine=hare
machine=eccc-ppp2
#abs="${HOME}/data_maestro/ords/oavar_abs/ensmanip_sles-11-broadwell-64-xc40-m_2.2.2-34-ga7bac3c_M.Abs"
abs="${HOME}/data_maestro/ords/oavar_abs/ensmanip_ubuntu-14.04-amd64-64-m_2.2.2-34-ga7bac3c_M.Abs"
#ensdir="/home/mab001/data_maestro/hare/ensmanip/ensemble/"
#ensdir="/home/mab001/sitestore1/maestro_archives/Tests/EnVar_small_v001/inputs/ensemble/"
#ensdir="/home/mab001/data_maestro/hare/ensmanip/kal557/"
ensdir="/home/skal001/data_maestro/eccc-ppp2/exp/archive_trial/kal557/"
npex=1
npey=267
openmp=1
maxcputime=600
run_in_parallel="/fs/ssm/eccc/mrd/rpn/utils/16.2/all/bin/r.run_in_parallel_1.1.28c"

#
# Don't modify below ...
#

ensdate=$1
if [ "${ensdate}" = "" ]; then
  echo ""
  echo " *************************************************"
  echo " ERROR: NO DATE WAS SPECIFIED ON THE COMMAND LINE!"
  echo " PLEASE GIVE THE VALID DATE FOR THE ENSEMBLES."
  echo " *************************************************"
  echo ""
  exit
fi
gest="${HOME}/data_maestro/${machine}/ensmanip/test_${ensdate}/"

# build the namelist
cat << EOF > $TMPDIR/flnml
 &NAMCT0
/
 &NAMENSMANIP
   OUTPUT_ENSEMBLE_MEAN = T
   OUTPUT_ENSEMBLE_PERTURBATIONS = T
   NENS = 8
   DATE = ${ensdate}
   WRITE_MPI = F
/
 &NAMTIME
  DSTEPOBSINC = 1.0d0
/
 &NAMSTATE
   ANLVAR(1)   ='UU'
   ANLVAR(2)   ='VV'
   ANLVAR(3)   ='TT'
   ANLVAR(4)   ='HU'
   ANLVAR(5)   ='P0'
   ANLVAR(6)   ='TG'
/
EOF

abs_basename=`basename $abs`

echo
echo "Launching ENSMANIP using..."
echo
echo "Working machine      :" $machine
echo "Working directory    :" $gest
echo "Executable file path :" $abs
echo "Executable file name :" $abs_basename
echo "Topology             :" ${npex}x${npey}x${openmp}
echo

ssh $machine rm -rf $gest
ssh $machine mkdir -p $gest
ssh $machine ln -s ${ensdir} ${gest}/ensemble
scp $TMPDIR/flnml ${machine}:${gest}/flnml
scp $abs ${machine}:${gest}/ensmanip.abs
ssh $machine ls -l $gest

cat << EOF > $TMPDIR/go_ensmanip.sh
#!/bin/bash
set -ex
 echo "!!STARTING SCRIPT!!"
. ssmuse-sh -d eccc/mrd/rpn/utils/16.2

# MPI SETUP FOR PPP ONLY
if [ "\${EC_ARCH}" = Linux_x86-64 ]; then
  . ssmuse-sh -d hpco/tmp/eccc/201402/05/base
  . ssmuse-sh -d main/opt/intelcomp/intelcomp-2016.1.156
  . ssmuse-sh -d main/opt/openmpi/openmpi-1.6.5/intelcomp-2016.1.156 
  export OMPI_MCA_orte_tmpdir_base=/run/shm
  export OMPI_MCA_btl_openib_if_include=mlx5_0
fi

 cd $gest
 export TMG_ON=YES
 export OAVAR_BURP_SPLIT=yes
 ${run_in_parallel} -pgm ./ensmanip.abs -npex ${npex} -npey ${npey} -processorder -tag -nocleanup -verbose
EOF

cat << EOF > $TMPDIR/ptopo_nml
 &ptopo
  npex=$npex
  npey=$npey
/
EOF
scp $TMPDIR/ptopo_nml ${machine}:${gest}

ord_soumet $TMPDIR/go_ensmanip.sh -mach $machine -mpi -t $maxcputime -cm 3000M -cpus ${npex}x${npey}x${openmp} -jn ensmanip_${ensdate} -waste
