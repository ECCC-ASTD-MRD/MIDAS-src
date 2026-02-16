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

CFLAGS = -strict -lburpc -Bstatic -lrmn -lApp -Bdynamic -static-intel -lifcore -lifport -lsqlite3

ifeq (${ORDENV_PLAT},rhel-8-icelake-64)
    COMPILER=inteloneapi-2022.1.2
    COMPILER_DYNLIBS=/fs/ssm/main/opt/intelcomp/${COMPILER}/intelcomp+mpi+mkl/all/lib
    CMDALIBS_VERSION=20251021
else ifeq (${ORDENV_PLAT},rhel-9-graniterapids-64)
    COMPILER=inteloneapi-2025.1.0
    COMPILER_DYNLIBS=/fs/ssm/main/opt/intelcomp/${COMPILER}/intelcomp+mpi+mkl/rhel-9-amd64-64/lib
    ## The version '20250202-1' only exists on HPCR-U3
    CMDALIBS_VERSION=20260202-1
else
    $(error This platform is not supported: ORDENV_PLAT=${ORDENV_PLAT})
endif

install : # defining default target

splitobs_pgm: version.h $(SRC)
	. r.load.dot rpn/code-tools/20251217/env/${COMPILER}; \
	make splitobs_libs

splitobs_libs: version.h $(SRC)
	. r.load.dot eccc/mrd/rpn/libs/20260202;                                                 \
	. r.load.dot eccc/mrd/rpn/utils/20260202/burp-tools_20.0.15-${COMP_ARCH}_${ORDENV_PLAT}; \
        . ssmuse-sh -d eccc/cmd/cmda/libs/${CMDALIBS_VERSION}/${COMP_ARCH};                      \
	make splitobs

splitobs: version.h $(SRC)
	s.cc -o $(PROGRAM)_$(ORDENV_PLAT) $(CFLAGS) $(DEBUG_OPTIONS) -O $(OPT_LEVEL) -Wl,-rpath="/fs/ssm/eccc/mrd/rpn/libs/20260202/${ORDENV_PLAT}/${COMP_ARCH}/lib:/fs/ssm/eccc/mrd/rpn/utils/20260202/${ORDENV_PLAT}/${COMP_ARCH}/lib:$(COMPILER_DYNLIBS)" $(SRC); rm version.h


version.h:
	echo "#define  VERSION   \"$(VERSION)\"" > version.h; echo "#define  VERSION_SHA1   \"$(VERSION_SHA1)\"" >> version.h

install: $(PROGRAM)_pgm
	cp $(PROGRAM)_$(ORDENV_PLAT) $(PGM)

clean:
	rm -f *.o $(PROGRAM)_$(ORDENV_PLAT)
