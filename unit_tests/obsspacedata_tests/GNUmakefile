# Check for the PFUNIT environment variable
ifndef PFUNIT
$(error The environment variable "PFUNIT" must be given a value which is the \
path for the installation of pFUnit.)
endif

PFUNIT_LIBDIR ?=$(PFUNIT)/lib
PFUNIT_INCDIR ?=$(PFUNIT)/include
PFUNIT_MODDIR ?=$(PFUNIT)/mod
PFUNIT_BINDIR ?=$(PFUNIT)/bin

.SUFFIXES :

.SUFFIXES : .o .ftn90 .o

SHELL = /bin/sh

ifeq ($(DEBUG),YES)
	FFLAGS += -g -O0
endif

FFLAGS +=-I$(PFUNIT_MODDIR)

OPTIMIZE =

#ifeq ($(MPI),YES)
## this should be: -libappl rpn_comm_40007
## the '-libappl' part is already in LDFLAGS
   MPI_LIB   = rpn_comm_40007
   MPI_FLAG  = -mpi
   CPUS      = -cpus 1x1x1
#else
#   MPI_LIB   = 
#  MPI_FLAG  = -mpi stubs  --> this switch is broken
#   CPUS      = 
#endif

LIBRMN=rmn_012

LDFLAGS = -includes $(PFUNIT_INCDIR) $(PFUNIT_MODDIR) -libpath $(PFUNIT_LIBDIR) -libappl pfunit $(MPI_LIB) $(MPI_FLAG)

F90_SRC = $(wildcard *.ftn90)
F90_OBJ = $(F90_SRC:.ftn90=.o)

EXECUTABLE = test_obs.x
#ifeq ($(MPI),YES)
   COMMAND = ord_soumet ./runit.ksh -mach alef -t 60 $(MPI_FLAG) $(CPUS) -listing ${HOME}/listings
#else
#   COMMAND = ./$(EXECUTABLE)
#endif

# Targets
.PHONY:	tests

DOLLAR=$$
tests: $(EXECUTABLE)
	chmod +x $(EXECUTABLE)
	echo ". s.ssmuse.dot mpich2" > runit.ksh
	echo "cd `pwd`" >> runit.ksh
	echo "export NEW_PFUNIT=test_obs.x" >> runit.ksh
	echo "export UNIT_TEST_DATA_DIR=`pwd`/TestObsSpace_data" >> runit.ksh
	echo "r.mpirun2 -npex 1 -pgm test_obs.x" >> runit.ksh
	echo "echo '======================='" >> runit.ksh
	echo "echo '  T E S T S   R U N :'" >> runit.ksh
	echo "echo '======================='" >> runit.ksh
	echo "cat .pFUnitLog" >> runit.ksh
	echo " "
	echo " "
	echo "  - - - listing will appear in $(DOLLAR)HOME/listings - - -"
	echo " "
	echo " "
	$(COMMAND)

$(EXECUTABLE): $(F90_OBJ)
	s.compile -obj $(F90_OBJ) $(LDFLAGS) -librmn $(LIBRMN) -o $@

%.o:%.ftn90
	echo "FFLAGS="$(FFLAGS)
	s.compile $(OPTIMIZE) -opt "= $(FFLAGS)" -src $<

clean:
	$(RM) $(F90_OBJ)
	$(RM) *.mod
	$(RM) *.f90

distclean: clean
	$(RM) $(EXECUTABLE)
	$(RM) runit.ksh
	$(RM) .pFUnitLog

help:
	echo 'gmake tests <MPI=YES> -- compile and soumet the test'


obsspacedata_mod.ftn90 :
	ln -sf ../../obsspacedata_mod.ftn90 .
testobs_mod.o : obsspacedata_mod.o
testobsdriver.o : testobs_mod.o
