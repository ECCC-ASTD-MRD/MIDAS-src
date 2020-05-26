
###############################################################################
## 	MIDAS Makefile
##
##	@SYNOPSIS
##		* List all targets:
##        ```
##        make help
##        ```
##		* Build all programs (on a given architecture):
##		  ```
##		  source ./config.dot.sh
##		  make [-j${NCPUS} -O] [VERBOSE=(1|2)]
##		  make install
##		  ```
##
##		* to consider what is to be done after modifying the sources
##		  without actually doing it, use `make ${TARGET} -n`.
##		  (see also https://www.gnu.org/software/make/manual/html_node/Instead-of-Execution.html)
##
##		* To build a single target (object, absolute or dependency file) 
##		  ```
##		  . ./programs/commons/compile_setup.sh 
##		  make spectralfilter_mod.o 
##		  ```
##	@SEQUENCE
##		1. preprocessing of ftn90 files
##		2. dependencies generation
##		3. object compilation
##		4. linking
##
##	@AUTHORS
##		- Martin Deshaies-Jacques (@mad001) - CMDA - January 2020
## 		- Ervig Lapalm (@erv000), prior compilation solution
##		- Philippe Blain (@phb001), technical consulting and 
##			documentation proofreading
##	
##
##	@DEV_NOTES
##		* inline comments starting with `#>` are interpreted for `make help`
##		  (search for `help:`)
##		* need `makedepf90`
##		  available on 
##        [anaconda.org](https://anaconda.org/conda-forge/makedepf90)
##		  (awaiting pending installation request on ppp)
##		* prior to launching this make, `../programs/commons/compile_setup.sh`
##		  has to be sourced 
##		  (this exports `${COMPF}` `${GPP_OPTS}` and `${FOPTMIZ}`)
##
##
##	@TODO
##		- [ ] improve `recursiveDep.py` (output file as param so that 
##				stdout is used for other stuff)
##		- [x] install target, naming of absolute, etc.
##		- [x] compilation on backends
##		- [x] out-of-tree compilation
##		- [x] move back in a single make (? versus recursive make calls)
##		- [x] separate build directories for single and double precision
##				- aggregate each Abs group and COMPF prerequisite in 
##				  a single rule in `programs/programs.mk`
##		- [x] complete `programs/programs.mk`
##		- [x] can `seaIce` specific SSM be loaded with the others?
##				`module load cray-netcdf` is already in the common script
##				 -> passes Unit Tests
##		- [?] do we really need the automatic listings?
##			  	if so: careful with trailing space between 
##				target `$@` and `>`, it can produces buggy filename
##				with trailing spaces.
##		- [?] consider Busby 2017 2-pass approach 
##				(A Note on Compiling Fortran)
##		- [?] automated configuration sourcing with a single make call 
##			  	(consider ENKF#19 strategy: it gets really messy.
##			  	Maybe a script calling make and passing it's arguments?)
##
###############################################################################

##===================================================================
##===================================================================
##===================================================================
##	Common block
##

##-----------------------------------------------
## 	default target declaration
.DEFAULT_GOAL := all 
#.DEFAULT_GOAL := help


FILE_CFG := ./config.dot.sh

##	build directory **with respect to Makefile**
## 	can only be used in target.mk
## 	can be initialized from parent shell
ifneq ($(MAKECMDGOALS),help)
   ifndef DIR_BLD_ROOT
      $(error DIR_BLD_ROOT is not set, source $(FILE_CFG) )
   endif
endif

VERSION_SCRIPT := ../midas.version.sh


## 	those can be used in both (**no absolute** path here)
DIR_MOD_RELATIVE := modules
DIR_PGM_RELATIVE := programs

##  verbosity @phb001 clever trick
##
##  (does not work as expected when called recursively
##  for preproc of ftn90 for instance)
override V := @
ifdef VERBOSE
    V :=
    export VERBOSE
    ifeq ($(VERBOSE), 2)
        export Verbose=1 ## for s.f90
    endif
endif 

##  remove all default recipes
.SUFFIXES:

##  cleaning target
## 	(independent of environment: available anytime) 
.PHONY: clean cleanabs cleanobj cleandep



## TODO
clean:	#> delete the build directory 
	@echo "rm -rf $(DIR_BLD_ROOT)"
	@read -p "<!> Sure? (<ctr+c> to cancel, <enter> to continue)" && \
	rm -rf $(DIR_BLD_ROOT) 

cleanabs: cleanobj #> delete all absolutes
	rm -rf $(DIR_BLD_ROOT)/*/*/*/*.Abs

cleanobj: cleandep #> delete all objects
	rm -rf $(DIR_BLD_ROOT)/*/*/*/*.o $(DIR_BLD_ROOT)/*/*/*/*.mod 

cleandep:	#> delete all dependency file
	rm -rf $(DIR_BLD_ROOT)/*/*/*/*.inc $(DIR_BLD_ROOT)/*/*/*/*.f90


.PHONY: help
help:	#> print this help
	@echo "USAGE:"
	@echo "    source $(FILE_CFG)"
	@echo '    make [-j NCPUS -O] [OPTIONS] [TARGETS] [VERBOSE=(1|2)]'
	@echo "OPTIONS:"
	@echo '    consult make manual: man make'
	@echo "TARGETS:"
	@grep -E '^([0-9a-zA-Z_-]+|%[^:]+) *:.*?#> .*$$' $(MAKEFILE_LIST)|sort|\
	cut -d":" -f2- | \
	awk 'BEGIN {FS = ":[^:]*#> "}; \
		{printf "    %-30s %s\n", $$1, $$2}' | cut -d":" -f2





##  Verification for some environment dependencies
##	(not needed for some targets)
ifneq ($(MAKECMDGOALS),help)
ifeq ($(filter clean%,$(MAKECMDGOALS)),)
    ##  Verify if main environment variable set when needed
    ifndef COMPF
        $(error COMPF is not set, source $(FILE_CFG))
    endif
    ifndef FOPTMIZ
        $(error FOPTMIZ is not set, source $(FILE_CFG))
    endif
    ifndef GPP_OPTS
        $(error GPP_OPTS is not set, source $(FILE_CFG))
    endif
    ifndef ORDENV_PLAT
        $(error ORDENV_PLAT is not set, source $(FILE_CFG))
    endif
    ifndef EC_ARCH
        $(error EC_ARCH is not set, source $(FILE_CFG))
    endif

## 	Verify if `makedepf90` is available
##
##  until makedepf90 is globally installed, this check stays commmented
##  because `build.sh` will first make de dependencies on
##  the frontend where `makedepf90` is loaded through `conda` 
##  but compute node won't have it.
## 
#    ifeq ($(shell which makedepf90),)
#		##  `makedepf90` not available on backend..
#		##  dependencies have to be copied there before
#	    ##	but when launched on backend we don't 
#		##  want it to fail because `makedepf90` is not available
#        ifeq ($(ORDENV_PLAT),ubuntu-18.04-skylake-64)
#            $(error makedepf90 not available)
#        endif
#    endif
endif
endif


## TODO
##	END of Common block
##
##===================================================================
##===================================================================
##===================================================================
##
## 	Multi-architecture multi-directories strategy
##	(see http://make.mad-scientist.net/papers/multi-architecture-builds)
##
##	Synopsis:
##		1. `make` is called from `./src`
##		2. gets here and `target.mk` is included
##		3. `target.mk` move us into the build directory:
##         `$(DIR_BLD_ROOT)/$(revnum)/$(EC_ARCH)`
##		4. defines the relative position of the `Makefile` with
##		   $(SRC_ROOT) used in the `else` section
##		5. it calls `make` from there (allowing clean out-of-tree compilation)
##		6. gets here again, but now the `else` section is included
##		   and compilation rules are interpreted
##
##		<!> be very careful with `target.mk` it's filled with
##		    ancestral make black magic
##
ifneq (,$(filter src%,$(notdir $(CURDIR))))
include target.mk
else
##===================================================================
##===================================================================
##===================================================================
##
##	Main make rules
##
##	(only processed when make called from within the build directory
##	see ./target.mk)
##	
##  <!> I can see potential problems with line
##      `ifneq (,$(filter src%,$(notdir $(CURDIR))))`
##      It means that make, called from anywhere except from `src` would
##      trigger this part, even out of build dir...
##      TODO: improve that condition
##

FC := s.f90

DIR_MOD = $(SRC_ROOT)/$(DIR_MOD_RELATIVE)
DIR_PGM = $(SRC_ROOT)/$(DIR_PGM_RELATIVE)

VPATH = $(SRC_ROOT) $(DIR_MOD) $(DIR_PGM)

##	Common external libraries       
##		(see also `./programs/programs.mk`)
##
MPILIB := rpn_comm
LIBSYS := hpcoperf sqlite3
LIBRMN := rmnMP


##-----------------------------------------------
##	Modules	
##
## 	1. gathering of all module and program sources
##
SRC_F90 := $(wildcard $(DIR_MOD)/*.f90)
SRC_FTN := $(wildcard $(DIR_MOD)/*.ftn90)
SRC_PGM := $(wildcard $(DIR_PGM)/*.f90)
SRC_FTN_F90 := $(notdir $(subst .ftn90,.f90,$(SRC_FTN)))
SRC_ALL := $(SRC_F90) $(SRC_FTN_F90)  $(SRC_PGM) 
OBJ := $(notdir $(subst .f90,.o,$(SRC_ALL)))

##-- ABS defined in programs/programs.mk
DEP_OBJ	= 	dep.obj.inc
DEP_ABS	= 	dep.abs.inc

##  MIDAS installed programs with version signature
#ABS_INSTALL_PRFX = $(SRC_ROOT)/$(DIR_BLD_ROOT)/midas_abs/midas-
## DIR_BLD_ROOT is now absolute path
ABS_INSTALL_DIR = $(DIR_BLD_ROOT)/midas_abs
ABS_INSTALL_PRFX = $(ABS_INSTALL_DIR)/midas-
ABS_INSTALL_SUFX = _$(ORDENV_PLAT)-$(VERSION)
ABS_INSTALL = $(addprefix $(ABS_INSTALL_PRFX),$(subst .Abs,$(ABS_INSTALL_SUFX).Abs,$(ABS)))

##  abort if there is an error with dependencies generation
##	for instance if `makedepf90` not present
.DELETE_ON_ERROR: $(DEP_OBJ) $(DEP_ABS)

##	rule for install target : simple copy and renaming
$(ABS_INSTALL_PRFX)%$(ABS_INSTALL_SUFX).Abs : $(CURDIR)/%.Abs
	+@[ -d $(@D) ] || mkdir -p $(@D)
	cp $< $@

##----------------------------------------------
##     Main ftn90 preprocessing, compilation and linking recipes
##

_D_options=$(filter -D%,$(COMPF))
%.f90: %.ftn90 #> preprocess an ftn90 file
	@echo 'Preprocessing $(@F) inplace'
	$(V) r.gppf -lang-f90+ $(GPP_OPTS) $(_D_options) $< \
	> $@ \
	2> /dev/null

%.o:	#> compile an object 
	$(FC) $(COMPF) -O $(FOPTMIZ) -c $< \
	-o $@

%.Abs: 	#>link an absolute 
	$(FC) $(COMPF) -O $(FOPTMIZ) \
	$(addprefix -l, $(LIBAPPL) $(LIBSYS) $(LIBRMN)) $^ \
	-o $@

## 	include specific environment variables and specific rules
##		(must be included after the generic rule
##		all prerequisite rules are merged 
##  	but **only the last recipe** is used
##		so for a specific rule, it must be included in 
##		the following .mk file)
-include $(SRC_ROOT)/programs/programs.mk

## 	adding dependencies prerequisite
##	(on a different recipe than with the rule for proper usage of
##	 $^ and $<)
%.o: $(DEP_OBJ) 
%.Abs: $(DEP_ABS)

##-----------------------------------------------
## 	Dynamical include of automatically generated dependencies
## 		(see `depend` target)
##		note that dependencies differ from real4 to real8
##
ifeq ($(filter clean%,$(MAKECMDGOALS)),)
ifeq ($(filter install%,$(MAKECMDGOALS)),)
ifneq ($(MAKECMDGOALS), info)
    -include $(DEP_OBJ)
    -include $(DEP_ABS)
endif
endif
endif

##-----------------------------------------------
## 	automatic generation of dependencies using `makedepf90` 
##
##	`makedepf90` can be obtained through 
##	* conda pakage at [conda-forge](https://anaconda.org/conda-forge/makedepf90)
##	* Debian or Ubuntu pakages

.PHONY: copy_depend_backend

## `makedepf90` only available on frontend
ifeq ($(ORDENV_PLAT),$(FRONTEND_PLAT))
$(DEP_OBJ): $(SRC_ALL)
	$(V) echo "Generating object dependencies > $@"	
	$(V) echo $(notdir $(shell makedepf90 $(SRC_ALL))) | \
		sed  's/\([^ ]* :\)/\n\1/g' > $@

##  ### Absolute linking prerequisites
##
##	generate prerequisites for absolutes using `recursiveDep.py`
##
$(DEP_ABS): $(DEP_OBJ)
	$(V) echo "Generating executables dependencies > $@"
	$(V) $(SRC_ROOT)/recursiveDep.py '$(notdir $(ABS))' -d $^ > $@

copy_depend_backend: $(DEP_OBJ) $(DEP_ABS)	#> copy dependency files generated on frontend
	$(V) source $(SRC_ROOT)/func.dot.sh && copy_depend ${BACKEND} &&\
	echo "Dependency copied to backend build folder" 

else
$(DEP_OBJ) $(DEP_ABS) copy_depend_backend:
	@echo "<!> Dependency analysis only available on frontend"
	@exit 1
endif

##-----------------------------------------------
##	Phony global targets
##
.PHONY: all absolutes objects info ssm doc diagrams install depend 

all: absolutes	#> compile all programs 

absolutes: $(ABS) 

objects: $(OBJ)

ssm:	#> build SSM package (not available yet)

doc:	#> build documentation (not available yet)

diagrams:	#> build diagrams (not available yet)

.PHONY: printAbs
printAbs:
	@echo $(ABS)

.PHONY: printObj
printObj:
	@echo $(OBJ)

info:	#>	print build configuration
	@echo
	@echo '===| INFO |================================='
	@echo version: $(VERSION)
	@echo build directory: $(CURDIR)
	@echo install directory: $(ABS_INSTALL_DIR)
	@echo COMPF=$(COMPF)
	@echo FOPTMIZ=$(FOPTMIZ) 
	@echo GPP_OPTS=$(GPP_OPTS) 
	@echo '============================================'
	@echo

install: $(ABS_INSTALL)			#> install all programs

depend : $(DEP_OBJ) $(DEP_ABS)	#> generate all dependency files


##
##	END of Main make rules
##
##===================================================================
##===================================================================
##===================================================================
endif
# vim: set noexpandtab noautoindent nolist:
