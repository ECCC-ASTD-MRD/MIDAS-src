#!/bin/bash
#

#
# User-defined options
#
flnml="namelist.nml"
machine=hare
#machine=eccc-ppp1
abs="${HOME}/data_maestro/ords/oavar_abs/ensmanip_sles-11-broadwell-64-xc40-m_2.2.2-1-ge252ae0_M.Abs"
#abs="${HOME}/data_maestro/ords/oavar_abs/ensmanip_ubuntu-14.04-amd64-64-m_2.2.2-1-ge252ae0_M.Abs"
gest="${HOME}/data_maestro/${machine}/ensmanip/test/"
ensdir="/home/mab001/data_maestro/hare/ensmanip/ensemble/"
#ensdir="/home/mab001/sitestore1/maestro_archives/Tests/EnVar_small_v001/inputs/ensemble/"
npex=10
npey=10
openmp=1
maxcputime=300
run_in_parallel="/fs/ssm/eccc/mrd/rpn/utils/16.2/all/bin/r.run_in_parallel_1.1.28c"

#
# Don't modify below ...
#
abs_basename=`basename $abs`

echo
echo "Launching ENSMANIP using..."
echo
echo "Using namelist file  :" $flnml
echo "Working machine      :" $machine
echo "Working directory    :" $gest
echo "Executable file path :" $abs
echo "Executable file name :" $abs_basename
echo "Topology             :" ${npex}x${npey}x${openmp}
echo

ssh $machine rm -rf $gest
ssh $machine mkdir -p $gest
ssh $machine ln -s ${ensdir} ${gest}/ensemble
scp $flnml ${machine}:${gest}/flnml
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

ord_soumet $TMPDIR/go_ensmanip.sh -mach $machine -mpi -t $maxcputime -cm 3000M -cpus ${npex}x${npey}x${openmp} -jn ensmanip -waste
