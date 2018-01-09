#!/bin/bash
#

#
# User-defined options
#
#machine=eccc-ppp2
#abs="${HOME}/data_maestro/ords/midas_abs/midas-ensembleH_ubuntu-14.04-amd64-64-v_3.0.4-68-g19f22a7_M.Abs"

machine=brooks
#abs="${HOME}/data_maestro/ords/midas_abs/midas-ensembleH_sles-11-broadwell-64-xc40-v_3.0.4-67-ga8fcc8a_M.Abs"
abs="${HOME}/data_maestro/ords/midas_abs/midas-ensembleH_sles-11-broadwell-64-xc40-v_3.0.4-68-g19f22a7_M.Abs"

ensdir="/home/mab001/data_maestro/${machine}/kal569/with_gz"
obsdir="/home/mab001/data_maestro/${machine}/ensembleh/obssplit_16x16_noiasicris/"
coefsat="/home/scvs400/datafiles/constants/cmda/alt/v3.0.0/rtcoefsat"
statsat="/home/scvs400/datafiles/constants/cmda/alt/v3.0.0/statsat/"
obscov="/home/sanl000/ANAL_shared/datafiles/constants/arma/oavar/2.1.3/__STATOBS_CONV_201709__/"
npex=16
npey=16
openmp=2
maxcputime=600
run_in_parallel="/fs/ssm/eccc/mrd/rpn/utils/16.2/all/bin/r.run_in_parallel_1.1.28c"

#
# Don't modify below ...
#

ensdate=2017010100
gest="${HOME}/data_maestro/${machine}/ensembleh/test2_256_tinterp/"

# build the namelist
cat << EOF > $TMPDIR/flnml
 &NAMCT0
/
 &NAMENSEMBLEH
   NENS = 256
   DATE = ${ensdate}
   USETLMH = .false.
/
 &NAMTIME
  DSTEPOBS = 1.0d0
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
 &NAMDIMO
   NMXOBS=30000,
   NDATAMX=500000,
/
 &NAMRMAT
/
 &NAMFILT
  NELEMS = 16,
  NLIST(1)=10004,
  NLIST(2)=10051,
  NLIST(3)=12004,
  NLIST(4)=11215,
  NLIST(5)=11216,
  NLIST(6)=12203,
  NLIST(7)=12001,
  NLIST(8)=12192,
  NLIST(9)=11003,
  NLIST(10)=11004,
  NLIST(11)=10194,
  NLIST(12)=12062,
  NLIST(13)=12063,
  NLIST(14)=12163,
  NLIST(15)=15036,
  NLIST(16)=15031,
  NFLAGS=3,
  NLISTFLG(1)=2,
  NLISTFLG(2)=4,
  NLISTFLG(3)=5,
  LTOPOFILT=.TRUE.,
  RLIMLVHU=70.D0,
/
 &NAMTOV   
   NSENSORS =  30,
   CSATID(1)        = 'AQUA',
   CINSTRUMENTID(1) = 'AIRS',
   CSATID(2)        = 'METOP-2',
   CINSTRUMENTID(2) = 'IASI',
   CSATID(3)        = 'NOAA15',
   CINSTRUMENTID(3) = 'AMSUA',
   CSATID(4)        = 'NOAA16',
   CINSTRUMENTID(4) = 'AMSUA',
   CSATID(5)        = 'AQUA',
   CINSTRUMENTID(5) = 'AMSUA',
   CSATID(6)        = 'NOAA18',
   CINSTRUMENTID(6) = 'AMSUA',
   CSATID(7)        = 'NOAA15',
   CINSTRUMENTID(7) = 'AMSUB',
   CSATID(8)        = 'NOAA16',
   CINSTRUMENTID(8) = 'AMSUB',
   CSATID(9)        = 'NOAA17',
   CINSTRUMENTID(9) = 'AMSUB',
   CSATID(10)        = 'NOAA18',
   CINSTRUMENTID(10) = 'MHS',
   CSATID(11)        = 'DMSP16',
   CINSTRUMENTID(11) = 'SSMIS',
   CSATID(12)        = 'GOES13',
   CINSTRUMENTID(12) = 'GOESIMAGER',
   CSATID(13)        = 'MSG2',
   CINSTRUMENTID(13) = 'SEVIRI',
   CSATID(14)        = 'METEOSAT7',
   CINSTRUMENTID(14) = 'MVIRI',
   CSATID(15)        = 'METOP-2',
   CINSTRUMENTID(15) = 'AMSUA',
   CSATID(16)        = 'METOP-2',
   CINSTRUMENTID(16) = 'MHS',
   CSATID(17)        = 'NOAA19',
   CINSTRUMENTID(17) = 'AMSUA',
   CSATID(18)        = 'NOAA19',
   CINSTRUMENTID(18) = 'MHS',
   CSATID(19)        = 'METOP-1',
   CINSTRUMENTID(19) = 'IASI',
   CSATID(20)        = 'DMSP17',
   CINSTRUMENTID(20) = 'SSMIS',
   CSATID(21)        = 'DMSP18',
   CINSTRUMENTID(21) = 'SSMIS',
   CSATID(22)        = 'GOES15',
   CINSTRUMENTID(22) = 'GOESIMAGER',
   CSATID(23)        = 'MSG3',
   CINSTRUMENTID(23) = 'SEVIRI',
   CSATID(24)        = 'MTSAT1',
   CINSTRUMENTID(24) = 'GMSMTSAT',
   CSATID(25)        = 'MTSAT2',
   CINSTRUMENTID(25) = 'GMSMTSAT',
   CSATID(26)        = 'METOP-1',
   CINSTRUMENTID(26) = 'AMSUA',
   CSATID(27)        = 'METOP-1',
   CINSTRUMENTID(27) = 'MHS',
   CSATID(28)        = 'NPP',
   CINSTRUMENTID(28) = 'ATMS',
   CSATID(29)        = 'NPP',
   CINSTRUMENTID(29) = 'CRIS',
   CSATID(30)        = 'HMWARI-8',
   CINSTRUMENTID(30) = 'AHI',
   LDBGTOV    = .FALSE.,
   CRTMODL    = 'RTTOV',
/
 &NAMSAT
   LISTBURP(1)=054, ! METEOSAT 7
   LISTPLAT(1)="meteosat",
   LISTSAT(1)=007
   LISTBURP(2)=206, ! NOAA 15
   LISTPLAT(2)="noaa",
   LISTSAT(2)=015,
   LISTBURP(3)=207, ! NOAA 16
   LISTPLAT(3)="noaa",
   LISTSAT(3)=016,
   LISTBURP(4)=246, ! DMSP 13
   LISTPLAT(4)="dmsp",
   LISTSAT(4)=013,
   LISTBURP(5)=247, ! DMSP 14
   LISTPLAT(5)="dmsp",
   LISTSAT(5)=014,
   LISTBURP(6)=248, ! DMSP 15
   LISTPLAT(6)="dmsp",
   LISTSAT(6)=015,
   LISTBURP(7)=249, ! DMSP 16
   LISTPLAT(7)="dmsp",
   LISTSAT(7)=016,
   LISTBURP(8)=255, ! GOES 11
   LISTPLAT(8)="goes",
   LISTSAT(8)=011,
   LISTBURP(9)=256, ! GOES 12
   LISTPLAT(9)="goes",
   LISTSAT(9)=012,
   LISTBURP(10)=784,! AQUA (EOS 2)
   LISTPLAT(10)="eos",
   LISTSAT(10)=002,
   LISTBURP(11)=208,! NOAA 17
   LISTPLAT(11)="noaa",
   LISTSAT(11)=017,
   LISTBURP(12)=209,! NOAA 18
   LISTPLAT(12)="noaa",
   LISTSAT(12)=018,
   LISTBURP(13)=003,! METOP-1 (METOP-B)
   LISTPLAT(13)="metop",
   LISTSAT(13)=001,
   LISTBURP(14)=004,! METOP-2 (METOP-A)
   LISTPLAT(14)="metop",
   LISTSAT(14)=002,
   LISTBURP(15)=005,! METOP-3 (METOP-C)
   LISTPLAT(15)="metop",
   LISTSAT(15)=003,
   LISTBURP(16)=223,! NOAA 19
   LISTPLAT(16)="noaa",
   LISTSAT(16)=019,
   LISTBURP(17)=171,! MTSAT-1R
   LISTPLAT(17)="mtsat-1r",
   LISTSAT(17)=001,
   LISTBURP(18)=056,! METEOSAT 9 (MSG 2)
   LISTPLAT(18)="msg",
   LISTSAT(18)=002,
   LISTBURP(19)=257,! GOES 13
   LISTPLAT(19)="goes",
   LISTSAT(19)=013,
   LISTBURP(20)=285,! DMSP17
   LISTPLAT(20)="dmsp",
   LISTSAT(20)=017,
   LISTBURP(21)=286,! DMSP18
   LISTPLAT(21)="dmsp",
   LISTSAT(21)=018,
   LISTBURP(22)=224,! NPP
   LISTPLAT(22)="jpss",
   LISTSAT(22)=000,
   LISTBURP(23)=259,! GOES 15
   LISTPLAT(23)="goes",
   LISTSAT(23)=015, 
   LISTBURP(24)=172,! MTSAT-2
   LISTPLAT(24)="mtsat",
   LISTSAT(24)=002,
   LISTBURP(25)=258,! GOES 14
   LISTPLAT(25)="goes",
   LISTSAT(25)=014,
   LISTBURP(26)=057,! METEOSAT 10 alias MSG 3
   LISTPLAT(26)="msg",
   LISTSAT(26)=003,
   LISTBURP(27)=070,! METEOSAT 11 alias MSG 4
   LISTPLAT(27)="msg",
   LISTSAT(27)=004,
   LISTBURP(28)=173,! HMWARI-8
   LISTPLAT(28)="hmwari",
   LISTSAT(28)=008,
/
 &NAMINST
   LISTBURP(1)=050
   LISTINSTRUM(1)="atsr",
   LISTBURP(2)=203,
   LISTINSTRUM(2)="mhs",
   LISTBURP(3)=205,
   LISTINSTRUM(3)="mviri",
   LISTBURP(4)=207, 
   LISTINSTRUM(4)="seviri",
   LISTBURP(5)=221, 
   LISTINSTRUM(5)="iasi",
   LISTBURP(6)=295,
   LISTINSTRUM(6)="gmsim",
   LISTBURP(7)=296, 
   LISTINSTRUM(7)="gmsim",
   LISTBURP(8)=365, 
   LISTINSTRUM(8)="tmi",
   LISTBURP(9)=389, 
   LISTINSTRUM(9)="modis",
   LISTBURP(10)=420,
   LISTINSTRUM(10)="airs",
   LISTBURP(11)=570,
   LISTINSTRUM(11)="amsua",
   LISTBURP(12)=571,
   LISTINSTRUM(12)="amsua",
   LISTBURP(13)=572,
   LISTINSTRUM(13)="amsua",
   LISTBURP(14)=573,
   LISTINSTRUM(14)="amsua",
   LISTBURP(15)=574,
   LISTINSTRUM(15)="amsub",
   LISTBURP(16)=590,
   LISTINSTRUM(16)="avhrr",
   LISTBURP(17)=591,
   LISTINSTRUM(17)="avhrr",
   LISTBURP(18)=592,
   LISTINSTRUM(18)="avhrr",
   LISTBURP(19)=605,
   LISTINSTRUM(19)="hirs",
   LISTBURP(20)=606,
   LISTINSTRUM(20)="hirs",
   LISTBURP(21)=607,
   LISTINSTRUM(21)="hirs",
   LISTBURP(22)=615,
   LISTINSTRUM(22)="goesim",
   LISTBURP(23)=620,
   LISTINSTRUM(23)="cris",
   LISTBURP(24)=621,
   LISTINSTRUM(24)="atms",
   LISTBURP(25)=626,
   LISTINSTRUM(25)="goessd",
   LISTBURP(26)=623,
   LISTINSTRUM(26)="msu",
   LISTBURP(27)=627,
   LISTINSTRUM(27)="ssu",
   LISTBURP(28)=905,	
   LISTINSTRUM(28)="ssmi",
   LISTBURP(29)=908,
   LISTINSTRUM(29)="ssmis",
   LISTBURP(30)=2047,
   LISTINSTRUM(30)="airs",
   LISTBURP(31)=297,
   LISTINSTRUM(31)="ahi",
/ 
 &NAMCHANOFFSET
   LISTOFFSET(1)=27,
   LISTINSTRUM(1)="amsua",
   LISTOFFSET(2)=42,
   LISTINSTRUM(2)="amsub",
   LISTOFFSET(3)=42,
   LISTINSTRUM(3)="mhs",   
   LISTOFFSET(4)=3,
   LISTINSTRUM(4)="seviri",   
   LISTOFFSET(5)=18,
   LISTINSTRUM(5)="goesim",   
/
 &NAMTOVSINST
   inst_names(1)='amsua',
   inst_names(2)='amsub',
   inst_names(3)='mhs',
   inst_names(4)='ssmis',
   inst_names(5)='atms',
   inst_names(6)='airs',
   inst_names(7)='iasi',
   inst_names(8)='cris',
   inst_names(9)='radianceclear',
/
 &NAMHYPER
   name_inst(1)='airs',
   name_inst(2)='iasi',
   name_inst(3)='cris',
/
 &NAMGEO
   name_inst(1)='goesim',
   name_inst(2)='ahi',
   name_inst(3)='mviri',
   name_inst(4)='seviri',
   name_inst(5)='gmsim',
/
 &NAMCODTYP
/
 &NAMGPSRO
  LEVELGPSRO=2,
  SURFMIN  = 1000.D0,
  HSFMIN   = 1000.D0,
  HTPMAX   = 40000.D0,
  BGCKBAND = 0.05D0,
/
 &NAMGPSGB
   DZMIN       = 2.0D0,
   DZMAX       = 1000.0D0,
   YZTDERR     = 0.0D0,
   LASSMET     = .TRUE.,
   LLBLMET     = .TRUE.,
   YSFERRWGT   = 1.00D0,
   YZDERRWGT   = 1.00D0,
   LBEVIS      = .TRUE.,
   IREFOPT     = 1,
   L1OBS       = .FALSE.,
   LTESTOP     = .FALSE.,
   IZTDOP      = 1,
/
 &NAMBURP_FILTER_SFC
  BNBITSOFF=0,
  BBITOFF(1)=2,
  BBITOFF(2)=4,
  BBITOFF(3)=5,
  BBITOFF(4)=-5,
  NELEMS_SFC=10,
  BLISTELEMENTS_SFC=10004,12004,10051,12203,11011,11012,13220,15031,15032,15035
/
 &NAMBURP_FILTER_CONV
  UA_HIGH_PRECISION_TT_ES=.true.,
  READ_QI_GA_MT_SW=.true.,
  BNBITSOFF=0,
  BBITOFF(1)=2,
  BBITOFF(2)=4,
  BBITOFF(3)=5,
  BBITOFF(4)=-5,
  BBITOFF(5)=-8,
  NELEMS=10,
  BLISTELEMENTS=12001,11001,11002,12192,10194,15036,11011,11012,13210,13220
/
 &NAMBURP_FILTER_TOVS
  BNBITSOFF=0,
  BBITOFF(1)=2,
  BBITOFF(2)=4,
  BBITOFF(3)=5,
  BBITOFF(4)=-5,
  BBITOFF(5)=-8,
  NELEMS=1
  BLISTELEMENTS=12163
/
 &NAMPHY
  NEW_TETENS_COEFS=.true.
/
EOF

abs_basename=`basename $abs`

echo
echo "Launching ENSEMBLEH using..."
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
ssh $machine ln -s ${obsdir} ${gest}/obs
scp ${coefsat}/* ${machine}:${gest}/
scp ${statsat}/* ${machine}:${gest}/
scp ${obscov}/* ${machine}:${gest}/
scp $TMPDIR/flnml ${machine}:${gest}/flnml
scp $abs ${machine}:${gest}/ensembleh.abs
ssh $machine ls -l $gest

cat << EOF > $TMPDIR/go_ensembleh.sh
#!/bin/bash
set -ex
 echo "!!STARTING SCRIPT!!"
. ssmuse-sh -d eccc/mrd/rpn/utils/16.2
. ssmuse-sh -d eccc/mrd/rpn/anl/oavar/3.0.4

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

cat > run.sh <<EOFRUN
#!/bin/ksh

gdb -ex run -ex where ./ensembleh.abs
EOFRUN

chmod +x run.sh


 ${run_in_parallel} -pgm ./run.sh -npex ${npex} -npey ${npey} -processorder -tag -nocleanup -verbose
EOF

cat << EOF > $TMPDIR/ptopo_nml
 &ptopo
  npex=$npex
  npey=$npey
/
EOF
scp $TMPDIR/ptopo_nml ${machine}:${gest}

ord_soumet $TMPDIR/go_ensembleh.sh -mach $machine -mpi -t $maxcputime -cm 6500M -cpus ${npex}x${npey}x${openmp} -jn ensembleh -waste
