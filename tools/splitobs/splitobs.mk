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

VERSION_SHA1 ?= $(shell git rev-parse HEAD)
VERSION ?= $(shell git describe)


splitobs_ubuntu-18.04-skylake-64: version.h $(SRC)
	. ssmuse-sh -d eccc/mrd/rpn/code-tools/01.3; . ssmuse-sh -d hpco/exp/intelpsxe-cluster-19.0.3.199; make splitobs_libs_19

splitobs_sles-15-skylake-64-xc50: version.h $(SRC)
	. ssmuse-sh -d eccc/mrd/rpn/code-tools/01.3; module load PrgEnv-intel/6.0.5; make splitobs_libs_19

splitobs: version.h $(SRC)
	s.compile -o $(PROGRAM)_$(ORDENV_PLAT) -src $(SRC) -librmn rmn -libappl sqlite3 burp_c -O 2 $(OPT); rm version.h

splitobs_Wall: version.h $(SRC)
	s.compile -o $(PROGRAM)_$(ORDENV_PLAT) -src $(SRC) -optc =-Wall -librmn rmn -libappl sqlite3 burp_c -O 2 $(OPT); rm version.h

splitobs_libs_19: version.h $(SRC)
	. r.load.dot eccc/mrd/rpn/libs/19.5; . ssmuse-sh -d eccc/cmd/cmda/libs/19.5/${COMP_ARCH}; make splitobs OPT="$(OPT)"

version.h:
	echo "#define  VERSION   \"$(VERSION)\"" > version.h; echo "#define  VERSION_SHA1   \"$(VERSION_SHA1)\"" >> version.h
