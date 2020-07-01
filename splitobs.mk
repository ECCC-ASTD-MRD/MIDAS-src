##OBJECT: splitobs
#
##LIBRARIES: standard "c", rmn
#
##MODULES: fstdlib.h fstdlib.c ok_or_notok.h splitobs.c fstdlib.c mainF.f
#
##CMOI_PLATFORM  op_linux op_ibm
#
##BUILD_LIB no
#

# export AFSISIO=${AFSISIO:-/home/binops/afsi/sio/env_ibm/afsisio}

PROGRAM = splitobs
SRC = splitobs.c fstdlib.c mainF.f

SQLITE_VERSION = 3.3.17
BURPLIB_VERSION = 1.3
VERSION_SHA1 ?= $(shell git rev-parse HEAD)
VERSION ?= $(shell git describe)


splitobs_ubuntu-18.04-skylake-64: version.h $(SRC)
	. ssmuse-sh -d eccc/mrd/rpn/code-tools/01.3; . ssmuse-sh -d hpco/exp/intelpsxe-cluster-19.0.3.199; make splitobs_libs_19

splitobs: version.h $(SRC)
	s.compile -o $(PROGRAM)_$(ORDENV_PLAT) -src $(SRC) -librmn rmn -libappl sqlite3 burp_c -O 2 $(OPT); rm version.h
splitobs_Wall: version.h $(SRC)
	s.compile -o $(PROGRAM)_$(ORDENV_PLAT) -src $(SRC) -optc =-Wall -librmn rmn -libappl sqlite3 burp_c -O 2 $(OPT); rm version.h


## on doit faire 's.use gmake as make' avant d'appeler 'make splitobs_AIX-powerpc7' sinon ca ne fonctionnera pas pour les 'make' recursifs.
splitobs_AIX-powerpc7: version.h $(SRC)
	. ssmuse-sh -d hpcs/201402/02/base -d hpcs/ext/xlf_13.1.0.10; make splitobs_libs

splitobs_Linux_x86-64: version.h $(SRC)
	. ssmuse-sh -d hpcs/201402/02/base -d hpcs/201402/02/intel13sp1u2; make splitobs_libs OPT="-optc='-fp-model precise'"

splitobs_libs: version.h $(SRC)
	. ssmuse-sh -d rpn/libs/15.2 -d cmda/libs/15.2/${COMP_ARCH}; make splitobs OPT="$(OPT)"

splitobs_ubuntu-14.04-amd64-64: version.h $(SRC)
	. ssmuse-sh -d hpco/tmp/eccc/201402/06/base -d hpco/exp/intel-2016.1.156; make splitobs_libs_16

splitobs_sles-15-skylake-64-xc50: version.h $(SRC)
	. ssmuse-sh -d eccc/mrd/rpn/code-tools/01.3; module load PrgEnv-intel/6.0.5; make splitobs_libs_19

splitobs_sles-11-amd64-64: version.h $(SRC)
	. ssmuse-sh -d hpco/tmp/eccc/201402/06/base; module load PrgEnv-intel/5.2.82; make splitobs_libs_16

splitobs_sles-11-haswell-64-xc40: version.h $(SRC)
	. ssmuse-sh -d hpco/tmp/eccc/201402/06/base; module load PrgEnv-intel/5.2.82; make splitobs_libs_16

splitobs_sles-11-broadwell-64-xc40: version.h $(SRC)
	. r.env.dot --arch sles-11-broadwell-64-xc40; . ssmuse-sh -d hpco/tmp/eccc/201402/06/base; module load PrgEnv-intel/5.2.82; make splitobs_libs_16

splitobs_libs_16: version.h $(SRC)
	. ssmuse-sh -d eccc/mrd/rpn/libs/16.2; . ssmuse-sh -d eccc/cmd/cmda/libs/16.2-4/${COMP_ARCH}; make splitobs OPT="$(OPT)"

splitobs_libs_19: version.h $(SRC)
	. r.load.dot eccc/mrd/rpn/libs/19.5; . ssmuse-sh -d eccc/cmd/cmda/libs/19.5/${COMP_ARCH}; make splitobs OPT="$(OPT)"

splitobs_linux: version.h $(SRC)
	r.compile -o $(PROGRAM)_$(BASE_ARCH) -src $(SRC) -librmn rmn -libsys m -includes /home/cmss/afsm/adm/ssm_access_points/cmda/burplib-c_$(BURPLIB_VERSION)_linux26-i386/include /data/ssm_op/sqlite/sqlite_$(SQLITE_VERSION)_linux26-i386/include -libpath /home/cmss/afsm/adm/ssm_access_points/cmda/burplib-c_$(BURPLIB_VERSION)_linux26-i386/lib /data/ssm_op/sqlite/sqlite_$(SQLITE_VERSION)_linux26-i386/lib -libappl sqlite3 burp_c -debug -O 0; rm version.h

splitobs_ibm: version.h $(SRC)
	r.compile -o $(PROGRAM)_ibm -src $(SRC) -librmn rmn -libsys m -includes /home/cmss/afsm/adm/ssm_access_points/cmda/burplib-c_$(BURPLIB_VERSION)_aix53-ppc-64/include /opt/ssm/sqlite_$(SQLITE_VERSION)_aix53-ppc-64/include -libpath /home/cmss/afsm/adm/ssm_access_points/cmda/burplib-c_$(BURPLIB_VERSION)_aix53-ppc-64/lib /opt/ssm/sqlite_$(SQLITE_VERSION)_aix53-ppc-64/lib -libappl sqlite3 burp_c; rm version.h

splitobs_ibm_ssm: version.h $(SRC)
	. s.ssmuse.dot devtools Xlf13.108 rmnlib-dev cmda; s.compile -o $(PROGRAM)_ibm -src $(SRC) -librmn rmn_013 -libappl sqlite3 burp_c -O 2; rm version.h

splitobs_ibm_nossm: version.h $(SRC)
	s.compile -o $(PROGRAM)_ibm -src $(SRC) -librmn rmn -libappl sqlite3 burp_c -O 2; rm version.h

splitobs_linux_nossm: version.h $(SRC)
	s.compile -o $(PROGRAM)_$(BASE_ARCH) -src $(SRC) -librmn rmn -libappl sqlite3 burp_c -O 2; rm version.h
        ## s.compile -o $(PROGRAM)_$(BASE_ARCH) -src $(SRC) -librmn rmn_013 -libappl sqlite3 burp_c -O 2; rm version.h

splitobs_ibm_debug: version.h $(SRC)
	s.compile -o $(PROGRAM)_$(BASE_ARCH) -src $(SRC) -librmn rmn -libappl sqlite3 burp_c -O 0 -debug; rm version.h

splitobs_ibm_efence: version.h $(SRC)
	s.compile -o $(PROGRAM)_$(BASE_ARCH) -src $(SRC) -librmn rmn -libappl sqlite3 burp_c efence_P7 -O 0 -debug; rm version.h

splitobs_linux_debug: version.h $(SRC)
	s.compile -o $(PROGRAM)_$(BASE_ARCH) -src $(SRC) -librmn rmn -libappl sqlite3 burp_c -O 0 -debug; rm version.h

splitobs_linux_efence: version.h $(SRC)
	s.compile -o $(PROGRAM)_$(BASE_ARCH) -src $(SRC) -librmn rmn -libappl efence sqlite3 burp_c -O 0 -debug; rm version.h

splitobs_linux_debug_all: version.h $(SRC)
	s.compile -o $(PROGRAM)_$(BASE_ARCH) -src $(SRC) -librmn rmn -includes $(PWD)/../burplib_c/include -libpath $(PWD)/../burplib_c/lib/linux26-x86-64 -libappl sqlite3 burp_c -O 0 -debug; rm version.h

splitobs_linux_debug_valgrind_burplib: version.h $(SRC)
	s.compile -o $(PROGRAM)_$(BASE_ARCH) -src $(SRC) -librmn rmn -includes $(PWD)/../burplib_c/include -libpath $(PWD)/../burplib_c/lib/linux26-x86-64 -libappl sqlite3 burp_c -O 0 -debug; rm version.h

splitobs_linux_debug_valgrind: version.h $(SRC)
	s.compile -o $(PROGRAM)_$(BASE_ARCH) -src $(SRC) -librmn rmn -libappl sqlite3 burp_c -O 0 -debug; rm version.h

exec_valgrind: splitobs_linux_debug_valgrind
	cd work; rm -vf 2011070112_to_amsub*; valgrind -v --show-reachable=yes --read-var-info=yes --leak-check=full --track-origins=yes ../$(PROGRAM)_$(BASE_ARCH) -burpin ../data/inputs/2011070112_to_amsub.2_enrgs -burpout 2011070112_to_amsub -fstin /users/dor/arma/erv/data/ords/programs/bgck.clipobs/grids/glbgrid -nomvar P0 -min_i 1 -max_i =1025 -min_j =0 -max_j 800 -verbose 6 > list 2> list.err

exec_valgrind_ua: splitobs_linux_debug_valgrind
	cd work; rm -vf 2011070112_to_amsub*; valgrind -v --show-reachable=yes --read-var-info=yes --leak-check=full --track-origins=yes ../$(PROGRAM)_$(BASE_ARCH) -burpin /users/dor/arma/erv/data/ords/programs/bgck.clipobs/data/inputs/2011070112_ua.north_south -burpout 2011070112_ua -fstin /users/dor/arma/erv/data/ords/programs/bgck.clipobs/grids/glbgrid -nomvar P0 -min_i 1 -max_i =1025 -min_j =0 -max_j 800 -npex 1 -npey 2 -verbose 6 > list 2> list.err

exec_valgrind_ua4d: splitobs_linux_debug
	cd work; rm -f 2011031912_ua4d.59134*; valgrind -v --show-reachable=yes --read-var-info=yes --leak-check=full --track-origins=yes ../$(PROGRAM)_$(BASE_ARCH) -verbose 6 -burpin /users/dor/arma/erv/data/ords/programs/bgck.clipobs/data/inputs/2011031912_ua4d.59134 -burpout 2011031912_ua4d.59134 -fstin /home/dormrb02/modeles/ANAL_shared/data/analysis_grid_prototypes/analysis_grid_prototype_glb_800x400 -nomvar P0 -npex 1 -npey 80 > list 2> list.err

exec_debug: splitobs_linux_debug_all
	cd work; rm -vf 2011070112_to_amsub*; ../$(PROGRAM)_$(BASE_ARCH) -burpin ../data/inputs/2011070112_to_amsub -burpout 2011070112_to_amsub -fstin /users/dor/arma/erv/data/ords/programs/bgck.clipobs/grids/glbgrid -nomvar P0 -min_i 1 -max_i =1025 -min_j =0 -max_j 800 -npex 1 -npey 2 -verbose 6 > list 2> list.err

exec_amsub: splitobs_linux_debug
	cd work; rm -vf 2011070112_to_amsub*; ../$(PROGRAM)_$(BASE_ARCH) -burpin ../data/inputs/2011070112_to_amsub.2_enrgs -burpout 2011070112_to_amsub -fstin /users/dor/arma/erv/data/ords/programs/bgck.clipobs/grids/glbgrid -nomvar P0 -min_i 1 -max_i =1025 -min_j =0 -max_j 800 -verbose 6 > list 2> list.err

exec_csr: splitobs_linux_debug
	cd work; rm -vf 2011070112_csr*; ../$(PROGRAM)_$(BASE_ARCH) -burpin ../data/inputs/2011070112_csr.1stn -burpout 2011070112_csr -fstin /users/dor/arma/erv/data/ords/programs/bgck.clipobs/grids/lamgrid -nomvar TG -verbose 6 > list 2> list.err

exec_csr_glb: splitobs_linux_debug
	cd work; rm -vf 2011070112_csr*; ../$(PROGRAM)_$(BASE_ARCH) -burpin ../data/inputs/2011070112_csr.glb -burpout 2011070112_csr -fstin /users/dor/arma/erv/data/ords/programs/bgck.clipobs/grids/lamgrid -nomvar TG -verbose 6 > list 2> list.err

exec_debug_ua: splitobs_linux_debug_all
	cd work; rm -vf 2011070112_ua*; ../$(PROGRAM)_$(BASE_ARCH) -burpin ../data/inputs/2011070112_ua.north_south.1_enrgs -burpout 2011070112_ua -fstin /users/dor/arma/erv/data/ords/programs/bgck.clipobs/grids/glbgrid -nomvar P0 -min_i 1 -max_i =1025 -min_j =0 -max_j 800 -npex 1 -npey 2 -verbose 6 > list 2> list.err

exec_debug_ua4d: splitobs_linux_debug_all
	cd work; rm -vf 2011070112_ua4d*; ../$(PROGRAM)_$(BASE_ARCH) -burpin ../data/inputs/2011070112_ua4d.1_enrgs -burpout 2011070112_ua4d -fstin /users/dor/arma/erv/data/ords/programs/bgck.clipobs/grids/glbgrid -nomvar P0 -min_i 1 -max_i =1025 -min_j =0 -max_j 800 -npex 1 -npey 2 -verbose 6 > list 2> list.err

exec_debug_ua4d_94638: splitobs_linux_debug
	cd work; rm -vf 2011020118_ua4d*; ../$(PROGRAM)_$(BASE_ARCH) -burpin ../data/inputs/2011020118_ua4d.94638 -burpout 2011020118_ua4d.94638 -fstin /users/dor/arma/erv/data/ords/programs/bgck.clipobs/grids/glbgrid -nomvar P0 -npex 1 -npey 16 -verbose 6 > list 2> list.err

exec_ua4d: splitobs_linux_debug
	cd work; rm -vf 2011070112_ua4d*; ../$(PROGRAM)_$(BASE_ARCH) -burpin ../data/inputs/2011070112_ua4d -burpout 2011070112_ua4d -fstin /users/dor/arma/erv/data/ords/programs/bgck.clipobs/grids/glbgrid -nomvar P0 -min_i 1 -max_i =1025 -min_j =0 -max_j 800 -npex 1 -npey 2 -verbose 6 > list 2>&1

exec_ua4d_1: splitobs_linux_debug
	cd work; rm -vf 2011070112_ua4d*; ../$(PROGRAM)_$(BASE_ARCH) -burpin ../data/inputs/2011070112_ua4d.1_enrgs -burpout 2011070112_ua4d -fstin /users/dor/arma/erv/data/ords/programs/bgck.clipobs/grids/glbgrid -nomvar P0 -min_i 1 -max_i =1025 -min_j =0 -max_j 800 -npex 1 -npey 2 -verbose 6 > list 2>&1

clip_ua4d: splitobs_linux_debug
	cd work; rm -vf 2011070112_ua4d*; ../$(PROGRAM)_$(BASE_ARCH) -burpin ../data/inputs/2011070112_ua4d -burpout 2011070112_ua4d -fstin /users/dor/arma/erv/data/ords/programs/bgck.clipobs/grids/reggrid -nomvar P0 -verbose 6 > list 2>&1

exec_ua4d_34: splitobs_linux_debug
	cd work; rm -vf 2011070112_ua4d*; ../$(PROGRAM)_$(BASE_ARCH) -burpin ../data/inputs/2011070112_ua4d_btyp34 -burpout 2011070112_ua4d -fstin /users/dor/arma/erv/data/ords/programs/bgck.clipobs/grids/glbgrid -nomvar P0 -min_i 1 -max_i =1025 -min_j =0 -max_j 800 -verbose 6 > list 2>&1

exec_ua4d_ASEU02: splitobs_linux_debug
	cd work; rm -vf 2011070112_ua4d*; ../$(PROGRAM)_$(BASE_ARCH) -burpin ../data/inputs/2011070112_ua4d.ASEU02 -burpout 2011070112_ua4d -fstin /users/dor/arma/erv/data/ords/programs/bgck.clipobs/grids/glbgrid -nomvar P0 -min_i 1 -max_i =1025 -min_j =0 -max_j 800 -verbose 6 > list 2>&1

exec_ua4d_resume: splitobs_linux_debug
	cd work; rm -vf 2011070112_ua4d*; ../$(PROGRAM)_$(BASE_ARCH) -burpin ../data/inputs/2011070112_ua4d.resume -burpout 2011070112_ua4d -fstin /users/dor/arma/erv/data/ords/programs/bgck.clipobs/grids/glbgrid -nomvar P0 -min_i 1 -max_i =1025 -min_j =0 -max_j 800 -verbose 6 > list 2>&1

exec_gdb: splitobs_linux_debug
	cd work; rm -vf 2011070112_ua4d*;  gdb -x ../gdb_commands ../$(PROGRAM)_$(BASE_ARCH)

version.h:
	echo "#define  VERSION   \"$(VERSION)\"" > version.h; echo "#define  VERSION_SHA1   \"$(VERSION_SHA1)\"" >> version.h
