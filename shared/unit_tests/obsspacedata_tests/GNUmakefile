
ifndef PFUNIT
  PFUNIT =/users/dor/arma/jbl/home1/UnitTesting/pFUnit/lib/$(ARCH)/mpi
endif
PFUNIT_LIBDIR ?=$(PFUNIT)/lib
PFUNIT_INCDIR ?=$(PFUNIT)/include
PFUNIT_MODDIR ?=$(PFUNIT)/mod
PFUNIT_BINDIR ?=$(PFUNIT)/bin

MPI ?=YES

# Check for the PFUNIT environment variable
ifndef PFUNIT
   $(error The environment variable "PFUNIT" must be given a value which is the \
path for the installation of pFUnit.)
endif

.SUFFIXES :

.SUFFIXES : .o .ftn90 .o

.DEFAULT_GOAL : unit_tests

SHELL = /bin/sh

ifeq ($(DEBUG),YES)
	OPTIMIZE = -debug -O 0
else
	OPTIMIZE =
endif

COMPILER_SSM = undefined
ifeq ($(BASE_ARCH),Linux)
   COMPILER_SSM = /ssm/net/hpcs/201402/01/base /ssm/net/hpcs/201402/01/intel13sp1u2
   MPI_PKG = /ssm/net/hpcs/201402/01/intel13sp1u2
   ifeq ($(MPI),YES)
       COMPILER_SSM += /ssm/net/rpn/utils/15.0
   endif
endif
ifeq ($(BASE_ARCH),Linux_x86-64)
   COMPILER_SSM = /ssm/net/hpcs/201402/01/base /ssm/net/hpcs/201402/01/intel13sp1u2
   MPI_PKG = /ssm/net/hpcs/201402/01/intel13sp1u2
   ifeq ($(MPI),YES)
       COMPILER_SSM += /ssm/net/rpn/utils/15.0
   endif
endif
ifeq ($(BASE_ARCH),AIX-powerpc7)
   COMPILER_SSM = /ssm/net/hpcs/201402/01/base /ssm/net/hpcs/ext/xlf_13.1.0.10
   MPI_PKG = /ssm/net/hpcs/ext/xlf_13.1.0.10
endif
ifeq ($(COMPILER_SSM),undefined)
   $(error The variable BASE_ARCH=$(BASE_ARCH) is not supported)
endif

FFLAGS +=-I$(PFUNIT_MODDIR)

ifeq ($(MPI),YES)
## this should be: -libappl rpn_comm_40007
## the '-libappl' part is already in LDFLAGS
   MPI_LIB   = rpn_comm
   MPI_FLAG  = -mpi
   CPUS      = 2
else
   MPI_LIB   = 
   MPI_FLAG  = rpn_commstubs
   CPUS      = 1
endif

LIBRMN=rmn_015

LDFLAGS = -includes $(PFUNIT_INCDIR) $(PFUNIT_MODDIR) -libpath $(PFUNIT_LIBDIR) -libappl pfunit $(MPI_LIB) $(MPI_FLAG)

F90_SRC = $(wildcard *.ftn90)
F90_OBJ = $(F90_SRC:.ftn90=.o)
F90_OBJ_NO_DRIVER = $(filter-out testobsdriver.o,$(F90_OBJ))

EXECUTABLE = test_obs.x

MACH_SUBMIT = $(TRUE_HOST)
ifeq ($(MACH_SUBMIT), erg)
   MACH_SUBMIT = pollux
else ifeq ($(MACH_SUBMIT), joule)
   MACH_SUBMIT = pollux
endif

#ifeq ($(MPI),YES)
   COMMAND = ord_soumet ./runit.ksh -mach $(MACH_SUBMIT) -m 10G -t 60 $(MPI_FLAG) -cpus $(CPUS) -listing $(PWD)
#else
#   COMMAND = ./$(EXECUTABLE)
#endif

# Targets (the first target listed is the default target)

unit_tests:
	. s.ssmuse.dot $(COMPILER_SSM) /ssm/net/rpn/libs/15.0; \
	${MAKE} tests MPI=$(MPI)

.PHONY:	tests
tests: $(EXECUTABLE)
	chmod +x $(EXECUTABLE)
	@echo ". s.ssmuse.dot /ssm/net/rpn/utils/15.0" > runit.ksh
	@echo ". s.ssmuse.dot $(MPI_PKG)" >> runit.ksh
	@echo "set -ex" >> runit.ksh
	@echo "cd `pwd`" >> runit.ksh
	@echo "rm -f flnml" >> runit.ksh
	@echo "rm -f .pFUnitLog" >> runit.ksh
	@echo "export NEW_PFUNIT=test_obs.x" >> runit.ksh
	@echo "export UNIT_TEST_OBS_DATA_DIR=`pwd`/TestObsSpace_data" >> runit.ksh
	@echo "r.run_in_parallel -npex $(CPUS) -pgm ./test_obs.x" >> runit.ksh
#	@echo "export DISPLAY=feynman:0.0" >> runit.ksh
#	@echo "totalview -mpi 'poe - AIX' -procs $(CPUS) -nodes $(CPUS) -display feynman:0.0 ./test_obs.x" >> runit.ksh
	@echo "set +x" >> runit.ksh
	@echo "echo '======================='" >> runit.ksh
	@echo "echo '  T E S T S   R U N :'" >> runit.ksh
	@echo "echo '   <TEST ALSO THE OTHER GROUP,>'" >> runit.ksh
	@echo "echo '   <  IN testobsdriver.ftn90  >'" >> runit.ksh
	@echo "echo '======================='" >> runit.ksh
	@echo "cat .pFUnitLog" >> runit.ksh
	@echo " "
	@echo " "
	@echo "  - - - listing will appear in `pwd`/runit.* - - -"
	@echo " "
	@echo " "
	$(COMMAND)

$(EXECUTABLE): $(F90_OBJ)
	s.compile -obj $(F90_OBJ) $(LDFLAGS) -librmn $(LIBRMN) -o $@

.PHONY: OBJ_NO_DRIVER
OBJ_NO_DRIVER: $(F90_OBJ_NO_DRIVER)

%.o:%.ftn90
	s.compile $(OPTIMIZE) $(MPI_FLAG) -opt "= $(FFLAGS)" -defines =-DUNIT_TESTING -src $<

.PHONY: clean
clean:
	$(RM) $(F90_OBJ)
	$(RM) *.mod
	$(RM) *.f90
	$(RM) core

distclean: clean
	$(RM) $(EXECUTABLE)
	$(RM) runit.ksh
	$(RM) .pFUnitLog
	$(RM) flnml
	$(RM) *.out
	$(RM) *~

help:
	@echo 'gmake -- compile and soumet the unit tests of obsspacedata_mod'


testobs_mod.o : obsspacedata_mod.o
testobsdriver.o : testobs_mod.o
