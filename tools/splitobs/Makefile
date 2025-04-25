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
VERSION ?= $(shell git describe --first-parent)

DEBUG_OPTIONS = -debug -optc =-w3 =-fp-speculation=safe
DEBUG_OPTIONS =

OPT_LEVEL = 2

BURPTOOLS = /fs/ssm/eccc/mrd/rpn/utils/20250604-beta/burp-tools_20.0.8.2

install : # defining default target

splitobs_rhel-8-icelake-64: version.h $(SRC)
	. r.load.dot rpn/code-tools/20250521/env/inteloneapi-2022.1.2; \
	make splitobs_libs

splitobs_libs: version.h $(SRC)
	. r.load.dot eccc/mrd/rpn/libs/20250604-beta; \
	. r.load.dot $(BURPTOOLS)-${COMP_ARCH}_${ORDENV_PLAT}; \
	make splitobs

splitobs: version.h $(SRC)
	s.cc -o $(PROGRAM)_$(ORDENV_PLAT) -lburpc -lrmn -static-intel -lifcore -lifport -lsqlite3 $(DEBUG_OPTIONS) -O $(OPT_LEVEL) -Wl,-rpath=/fs/ssm/main/opt/intelcomp/inteloneapi-2022.1.2/intelcomp+mpi+mkl/all/lib/,-rpath=${App_DIR}/lib,-rpath=$(BURPTOOLS)-${COMP_ARCH}_${ORDENV_PLAT}/lib $(SRC); rm version.h


version.h:
	echo "#define  VERSION   \"$(VERSION)\"" > version.h; echo "#define  VERSION_SHA1   \"$(VERSION_SHA1)\"" >> version.h

install: $(PROGRAM)_$(ORDENV_PLAT)
	cp $(PROGRAM)_$(ORDENV_PLAT) $(PGM)

clean:
	rm -f *.o $(PROGRAM)_$(ORDENV_PLAT)
