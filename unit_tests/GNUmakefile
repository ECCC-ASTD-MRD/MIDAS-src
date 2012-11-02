
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

ifeq ($(BASE_ARCH),Linux)
   COMPILER_SSM = pgi9xx
   ifeq ($(MPI),YES)
       COMPILER_SSM += mpich2
   endif
else
   ifeq ($(BASE_ARCH),AIX-powerpc7)
      COMPILER_SSM = Xlf13.108
   else
      $(error The variable BASE_ARCH=$(BASE_ARCH) is not supported)
   endif
endif

FFLAGS +=-I$(PFUNIT_MODDIR)

ifeq ($(MPI),YES)
## this should be: -libappl rpn_comm_40007
## the '-libappl' part is already in LDFLAGS
   MPI_LIB   = rpn_comm_40007
   MPI_FLAG  = -mpi
   CPUS      = 2
else
   MPI_LIB   = 
   MPI_FLAG  = rpn_commstubs_40007
   CPUS      = 1
endif

LIBRMN=rmn_013

LDFLAGS = -includes $(PFUNIT_INCDIR) $(PFUNIT_MODDIR) -libpath $(PFUNIT_LIBDIR) -libappl pfunit $(MPI_LIB) $(MPI_FLAG)

F90_SRC = $(wildcard *.ftn90)
F90_OBJ = $(F90_SRC:.ftn90=.o)

EXECUTABLE = test_master.x

MACH_SUBMIT = $(TRUE_HOST)
ifeq ($(MACH_SUBMIT), erg)
   MACH_SUBMIT = alef
endif

#ifeq ($(MPI),YES)
   COMMAND = ord_soumet ./runit.ksh -mach $(MACH_SUBMIT) -t 60 $(MPI_FLAG) -cpus $(CPUS) -listing $(PWD)
#else
#   COMMAND = ./$(EXECUTABLE)
#endif

# Targets

.PHONY:	tests
tests: $(EXECUTABLE)
	chmod +x $(EXECUTABLE)
	@echo ". s.ssmuse.dot mpich2" > runit.ksh
	@echo "set -ex" >> runit.ksh
	@echo "cd `pwd`" >> runit.ksh
	@echo "rm -f flnml" >> runit.ksh
	@echo "export NEW_PFUNIT=test_master.x" >> runit.ksh
	@echo "export UNIT_TEST_OBS_DATA_DIR=`pwd`/obsspacedata_tests/TestObsSpace_data" >> runit.ksh
	@echo "r.mpirun2 -npex $(CPUS) -pgm ./test_master.x" >> runit.ksh
	@echo "set +x" >> runit.ksh
	@echo "echo '======================='" >> runit.ksh
	@echo "echo '  T E S T S   R U N :'" >> runit.ksh
	@echo "echo '======================='" >> runit.ksh
	@echo "cat .pFUnitLog" >> runit.ksh
	@echo " "
	@echo " "
	@echo "  - - - listing will appear in `pwd`/runit.* - - -"
	@echo " "
	@echo " "
	$(COMMAND)

$(EXECUTABLE): $(F90_OBJ)
	s.compile -obj *.o $(LDFLAGS) -librmn $(LIBRMN) -o $@

%.o:%.ftn90
	s.compile $(OPTIMIZE) $(MPI_FLAG) -opt "= $(FFLAGS)" -defines =-DUNIT_TESTING -src $<

unit_tests:
	. s.ssmuse.dot $(COMPILER_SSM) devtools rmnlib-dev rpn_comm; \
	gmake tests MPI=$(MPI)

.PHONY: clean
clean:
	$(RM) *.o
	$(RM) *.mod
	$(RM) *.f90
	cd obsspacedata_tests; $(MAKE) clean

distclean: clean
	$(RM) $(EXECUTABLE)
	$(RM) runit.ksh
	$(RM) .pFUnitLog

help:
	@echo 'gmake -- compile and soumet the unit tests of all modules'

.PHONY: testobs_mod
testobs_mod : 
	cd obsspacedata_tests; $(MAKE) OBJ_NO_DRIVER
	cp obsspacedata_tests/*.mod .
	cp obsspacedata_tests/*.o .

testmasterdriver.o : testobs_mod
