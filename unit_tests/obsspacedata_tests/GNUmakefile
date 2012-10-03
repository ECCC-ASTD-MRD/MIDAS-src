
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
	FFLAGS += -g -O0
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

OPTIMIZE =

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
F90_OBJ_NO_DRIVER = $(filter-out testobsdriver.o,$(F90_OBJ))

EXECUTABLE = test_obs.x

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
	@echo "export NEW_PFUNIT=test_obs.x" >> runit.ksh
	@echo "export UNIT_TEST_OBS_DATA_DIR=`pwd`/TestObsSpace_data" >> runit.ksh
	@echo "r.mpirun2 -npex $(CPUS) -pgm ./test_obs.x" >> runit.ksh
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
	s.compile -obj $(F90_OBJ) $(LDFLAGS) -librmn $(LIBRMN) -o $@

.PHONY: OBJ_NO_DRIVER
OBJ_NO_DRIVER: $(F90_OBJ_NO_DRIVER)

%.o:%.ftn90
	s.compile $(OPTIMIZE) $(MPI_FLAG) -opt "= $(FFLAGS)" -src $<

unit_tests:
	. s.ssmuse.dot $(COMPILER_SSM) devtools rmnlib-dev rpn_comm; \
	gmake tests MPI=$(MPI)

.PHONY: clean
clean:
	$(RM) $(F90_OBJ)
	$(RM) *.mod
	$(RM) *.f90

distclean: clean
	$(RM) $(EXECUTABLE)
	$(RM) runit.ksh
	$(RM) .pFUnitLog

help:
	@echo 'gmake -- compile and soumet the unit tests of obsspacedata_mod'


testobs_mod.o : obsspacedata_mod.o
testobsdriver.o : testobs_mod.o
