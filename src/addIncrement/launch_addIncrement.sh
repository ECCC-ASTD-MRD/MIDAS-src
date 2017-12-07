#!/bin/bash
#

#
# User-defined options
#
machine=eccc-ppp2
abs="${HOME}/data_maestro/ords/oavar_abs/addIncrement_ubuntu-14.04-amd64-64-v_3.0.4-26-g6a1e6b2_M.Abs"
inputdir="/home/mab001/data_maestro/${machine}/addIncrement/inputs/"
npex=1
npey=139
openmp=1
maxcputime=120
run_in_parallel="/fs/ssm/eccc/mrd/rpn/utils/16.2/all/bin/r.run_in_parallel_1.1.28c"

#
# Don't modify below ...
#

gest="${HOME}/data_maestro/${machine}/addIncrement/workdir2/"

# build the namelist
cat << EOF > $TMPDIR/flnml
 &NAMCT0
/
 &NAMADDINC
/
 &NAMTIME
  dateFromTrials = .TRUE.
  DSTEPOBS = 3.0d0
  DSTEPOBSINC = 3.0d0
/
 &NAMSTATE
   ANLVAR(1)   ='UU'
   ANLVAR(2)   ='VV'
   ANLVAR(3)   ='TT'
   ANLVAR(4)   ='HU'
   ANLVAR(5)   ='P0'
   hInterpolationDegree = 'CUBIC'
/
EOF

abs_basename=`basename $abs`

echo
echo "Launching addIncrement using..."
echo
echo "Working machine      :" $machine
echo "Working directory    :" $gest
echo "Executable file path :" $abs
echo "Executable file name :" $abs_basename
echo "Topology             :" ${npex}x${npey}x${openmp}
echo

ssh $machine rm -rf $gest
ssh $machine mkdir -p $gest
ssh $machine ln -s ${inputdir} ${gest}/input
scp $TMPDIR/flnml ${machine}:${gest}/flnml
scp $abs ${machine}:${gest}/addIncrement.abs
ssh $machine ls -l $gest

cat << EOF > $TMPDIR/go_addIncrement.sh
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
 ls -l ./input/
 ln -s input/* .
 export TMG_ON=YES
 export OAVAR_BURP_SPLIT=yes
 ${run_in_parallel} -pgm ./addIncrement.abs -npex ${npex} -npey ${npey} -processorder -tag -nocleanup -verbose
EOF

cat << EOF > $TMPDIR/ptopo_nml
 &ptopo
  npex=$npex
  npey=$npey
/
EOF
scp $TMPDIR/ptopo_nml ${machine}:${gest}

ord_soumet $TMPDIR/go_addIncrement.sh -mach $machine -mpi -t $maxcputime -cm 1000M -cpus ${npex}x${npey}x${openmp} -jn addIncrement -waste
