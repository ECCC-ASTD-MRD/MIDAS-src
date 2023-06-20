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

PGM ?= ./midas.splitobs_$(ORDENV_PLAT).Abs

PROGRAM = splitobs
SRC = splitobs.c fstdlib.c

VERSION_SHA1 ?= $(shell git rev-parse HEAD)
VERSION ?= $(shell git describe)

DEBUG_OPTIONS = -debug -optc =-w3 =-fp-speculation=safe
DEBUG_OPTIONS =

OPT_LEVEL = 2

install : # defining default target

splitobs_rhel-8-icelake-64: version.h $(SRC)
	. r.load.dot eccc/mrd/rpn/code-tools/ENV/cdt-1.6.2/SCIENCE/inteloneapi-2022.1.2; \
	make splitobs_libs

splitobs_libs: version.h $(SRC)
	. r.load.dot eccc/mrd/rpn/libs/20230531;  \
	. ssmuse-sh -d eccc/cmd/cmda/libs/20230531/${COMP_ARCH}; \
	make splitobs

splitobs: version.h $(SRC)
	s.cc -o $(PROGRAM)_$(ORDENV_PLAT) -lburp_c -Bstatic -lrmn -lApp -Bdynamic -static-intel -lifcore -lifport -lsqlite3 $(DEBUG_OPTIONS) -O $(OPT_LEVEL) -Wl,-rpath=/fs/ssm/main/opt/intelcomp/inteloneapi-2022.1.2/intelcomp+mpi+mkl/all/lib/ $(SRC); rm version.h


version.h:
	echo "#define  VERSION   \"$(VERSION)\"" > version.h; echo "#define  VERSION_SHA1   \"$(VERSION_SHA1)\"" >> version.h

install: $(PROGRAM)_$(ORDENV_PLAT)
	cp $(PROGRAM)_$(ORDENV_PLAT) $(PGM)

clean:
	rm -f *.o $(PROGRAM)_$(ORDENV_PLAT)
