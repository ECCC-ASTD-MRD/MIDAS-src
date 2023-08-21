
###############################################################################
## 	MIDAS Makefile
##
##	@SYNOPSIS
##		* List all targets:
##        ```
##        make help
##        ```
##		* Clean some build elements
##        ```
##        make clean|cleanabs|cleanall|cleandep|cleanobj
##        ```
##		* Build all programs (on a given architecture):
##		  ```
##		  source ./config.dot.sh
##		  make [-j${NCORES} -O] [MIDAS_COMPILE_VERBOSE=(1|2)]
##		  make install
##		  ```
##
##		* to consider what is to be done after modifying the sources
##		  without actually doing it, use `make ${TARGET} -n`.
##		  (see also https://www.gnu.org/software/make/manual/html_node/Instead-of-Execution.html)
##
##		* To build a single target (object, absolute or dependency file) 
##		  ```
##		  source ./config.dot.sh
##		  make spectralfilter_mod.o 
##		  ```
##	@SEQUENCE
##		1. preprocessing of ftn90 files
##		2. dependencies generation
##		3. object compilation
##		4. linking
##
##  @REFERENCES 
##      * `./README.md`
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
##		* prior to launching this make, `./config.dot.sh`
##		  has to be sourced (except for help and clean* targets) 
##		  (this exports `${COMPF}` `${GPP_OPTS}` and `${FOPTMIZ}`)
##
##
##	@TODO
##		- [?] consider Busby 2017 2-pass approach 
##				(A Note on Compiling Fortran)
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


## will be redefined by target.mk
ifndef SRC_ROOT
SRC_ROOT=.
endif

FILE_CFG := $(SRC_ROOT)/config.dot.sh

VERSION_SCRIPT := $(SRC_ROOT)/../midas.version.sh

## 	those can be used in both (**no absolute** path here)
DIR_MOD_RELATIVE := modules
DIR_PGM_RELATIVE := programs

##  verbosity @phb001 clever trick
##
##  (does not work as expected when called recursively
##  for preproc of ftn90 for instance)
override V := @
ifdef MIDAS_COMPILE_VERBOSE
    V :=
    export MIDAS_COMPILE_VERBOSE
    ifeq ($(MIDAS_COMPILE_VERBOSE), 2)
        export Verbose=1 ## for s.f90
    endif
endif 

##  remove all default recipes
.SUFFIXES:

##	build directory **with respect to Makefile**
## 	can only be used in target.mk
## 	can be initialized from parent shell

## extract MIDAS_COMPILE_DIR_MAIN for simple __nopreq_goals target
## (those that can be called without sourcing config)
## extracting default value for DIR_* environment variables
ifndef VERSION
VERSION := $(shell $(VERSION_SCRIPT))
endif
__sed_env_extr='s/^.*=$${.*:-\(.*\)}/\1/'
ifndef MIDAS_COMPILE_DIR_MAIN
MIDAS_COMPILE_DIR_MAIN := $(shell cat $(FILE_CFG) | grep MIDAS_COMPILE_DIR_MAIN= | sed $(__sed_env_extr))
endif
ifndef MIDAS_ABS_LEAFDIR
MIDAS_ABS_LEAFDIR := $(shell cat $(FILE_CFG) | grep MIDAS_ABS_LEAFDIR= | sed $(__sed_env_extr))
endif
ifndef MIDAS_COMPILE_JOBNAME
MIDAS_COMPILE_JOBNAME := $(shell cat $(FILE_CFG) | grep MIDAS_COMPILE_JOBNAME= | sed $(__sed_env_extr))
endif

##  global cleaning target
## 	(independent of environment: available anytime) 
.PHONY: clean  cleanall cleanabs cleanobj cleandep cleanlst cleaninst


clean:	#> delete current build directory (keep installed programs) 
	rm -rf $(MIDAS_COMPILE_DIR_MAIN)/midas_bld-$(VERSION)

cleanall:	#> delete ALL build directories (keep installed programs) 
	@echo "rm -rf $(MIDAS_COMPILE_DIR_MAIN)/midas_bld-*"
	@read -p "<!> Sure? (<ctr+c> to cancel, <enter> to continue)" && \
		rm -rf $(MIDAS_COMPILE_DIR_MAIN)/midas_bld-*

cleaninst:	#> delete ALL installed programs
	@echo "rm -rf $(MIDAS_COMPILE_DIR_MAIN)/$(MIDAS_ABS_LEAFDIR)"
	@read -p "<!> Delete ALL installed programs? (<ctr+c> to cancel, <enter> to continue) " && \
		rm -rf $(MIDAS_COMPILE_DIR_MAIN)/$(MIDAS_ABS_LEAFDIR)

## specific cleaning target 
cleanabs:  #> delete all absolutes in current build directory
	rm -f $(MIDAS_COMPILE_DIR_MAIN)/midas_bld-$(VERSION)/*/*/*.Abs

cleanobj:  #> delete all objects in current build directory
	rm -f $(MIDAS_COMPILE_DIR_MAIN)/midas_bld-$(VERSION)/*/*/*.{o,mod}

cleandep:	#> delete all dependency files in current build directory
	rm -f $(MIDAS_COMPILE_DIR_MAIN)/midas_bld-$(VERSION)/*/*/*.{inc,f90} 

cleanlst:	#> delete all listings
	rm -f  $(MIDAS_COMPILE_JOBNAME).*.{out,o}



.PHONY: help
help:	#> print this help
	@echo "USAGE:"
	@echo '    - Cleaning targets (always available):'
	@echo "       make clean|cleanabs|cleanall|cleandep|cleanlst|cleanobj"
	@echo 
	@echo '    - For all targets but `clean*` or `help`:'
	@echo "       (source $(FILE_CFG) &&\\" 
	@echo '       make [-j $${NCORES} -O] [OPTIONS] [TARGETS] [MIDAS_COMPILE_VERBOSE=(1|2)])'
	@echo ''
	@echo '    - Complete install on current architecture:'
	@echo "       (source $(FILE_CFG) && make all && make install)"
	@echo "OPTIONS:"
	@echo '    consult make manual: man make'
	@echo "ENVIRONMENT"
	@echo "$$(cat $(FILE_CFG) | grep ^MIDAS_COMPILE.*=| sed 's/^.*{\([^:]*\):-\(.*\)}/    \1   (\2)/')"
	@echo "TARGETS:"
	@grep -E '^([0-9a-zA-Z_-]+|%[^:]+|[^:]+%) *:.*?#> .*$$' $(MAKEFILE_LIST)|sort|\
	cut -d":" -f2- | \
	awk 'BEGIN {FS = ":[^:]*#> "}; \
		{printf "    %-30s %s\n", $$1, $$2}' | cut -d":" -f2

##  Verification for some environment dependencies
##	(not needed for some targets, the __nopreq_goals)
__nopreq_goals := help clean cleanall cleaninst cleanabs cleanobj cleandep cleanlst splitobs
ifeq (,$(filter $(__nopreq_goals), $(MAKECMDGOALS)))
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
    ifeq ($(shell which makedepf90),)
        $(error makedepf90 not available)
    endif
    ifndef MIDAS_MAKEDEP_TIMEOUT
        $(error MIDAS_MAKEDEP_TIMEOUT is not set, source $(FILE_CFG))
    endif
endif

__ssm_goals := ssm  ssm_protect
ifneq (,$(filter $(__ssm_goals), $(MAKECMDGOALS)))
    ##  Verify if main environment variable set when needed
    ifndef MIDAS_SSM_TARGET
        $(error MIDAS_SSM_TARGET is not set, source $(FILE_CFG))
    endif
    ifndef MIDAS_SSM_VERSION
        $(error MIDAS_SSM_VERSION is not set, source $(FILE_CFG))
    endif
endif
MIDAS_SSM_PKGNAME ?= midas

__ssm_target_version := $(MIDAS_SSM_TARGET)/$(MIDAS_SSM_VERSION)
__ssm_workdir_base := $(MIDAS_COMPILE_DIR_MAIN)/midas_bld-$(VERSION)/ssm_build
__ssm_fullname := $(MIDAS_SSM_PKGNAME)_$(MIDAS_SSM_VERSION)_$(ORDENV_PLAT)
__ssm_workdir := $(__ssm_workdir_base)/$(__ssm_fullname)

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
##         `$(MIDAS_COMPILE_DIR_MAIN)/midas_bld-$(revnum)/$(EC_ARCH)`
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
LIBRMN := rmn-ompi


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

PGM := $(notdir $(subst .f90,,$(SRC_PGM)))
##-- program dependencies in programs/programs.mk
ABS := $(addsuffix .Abs,$(PGM)) 

DEP_OBJ	= 	dep.obj.inc
DEP_ABS	= 	dep.abs.inc

##  MIDAS installed programs with version signature
ABS_COMPILED := $(notdir $(wildcard $(CURDIR)/*.Abs))
ABS_INSTALL_DIR = $(MIDAS_COMPILE_DIR_MAIN)/$(MIDAS_ABS_LEAFDIR)
ABS_INSTALL_PRFX = $(ABS_INSTALL_DIR)/midas-
ABS_INSTALL_SUFX = _$(ORDENV_PLAT)-$(VERSION)
ABS_INSTALL = $(addprefix $(ABS_INSTALL_PRFX),\
				$(subst .Abs,$(ABS_INSTALL_SUFX).Abs,$(ABS_COMPILED)))

## splitobs as an external program has a different prefix
ABS_INSTALL := $(subst $(ABS_INSTALL_PRFX)splitobs,$(ABS_INSTALL_DIR)/midas.splitobs,$(ABS_INSTALL))

##  abort if there is an error with dependencies generation
##	for instance if `makedepf90` not present
.DELETE_ON_ERROR: $(DEP_OBJ) $(DEP_ABS)

##	rule for install target : simple copy and renaming
$(ABS_INSTALL_PRFX)%$(ABS_INSTALL_SUFX).Abs $(ABS_INSTALL_DIR)/midas.%$(ABS_INSTALL_SUFX).Abs :  $(CURDIR)/%.Abs
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

%.Abs: 	#> compile and link a program 
	@echo removing installed target if it exists
	rm -rf $(addprefix $(ABS_INSTALL_PRFX),\
		$(subst .Abs,$(ABS_INSTALL_SUFX).Abs,$@))
	$(FC) $(COMPF) -O $(FOPTMIZ)  $^ \
	$(addprefix -l, $(LIBAPPL) $(LIBSYS) $(LIBRMN)) \
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

##  External splitobs program compilation
.PHONY: splitobs splitobs.Abs
splitobs : splitobs.Abs
splitobs.Abs:	#> compile and install external splitobs program
		$(MAKE) -C $(SRC_ROOT)/../tools/splitobs \
			PGM=$(CURDIR)/splitobs.Abs && \
			$(MAKE) -C $(SRC_ROOT)/../tools/splitobs clean


## Compilation time modification to version_mod.f90
.INTERMEDIATE:  version_mod.pp.f90

version_mod.pp.f90: version_mod.f90
	sed "s|GIT-REVISION-NUMBER-WILL-BE-ADDED-HERE|$(VERSION)|g" $^ > $@ 
version_mod.o: version_mod.pp.f90
	$(FC) $(COMPF) -O $(FOPTMIZ) -c $< \
	-o $@

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
##		(in some cases, `makedepf90` will fail at analysing the dependencies and
## 		stall for a long time, a timer insure this stays under control)

$(DEP_OBJ): $(SRC_ALL) $(SRC_ROOT)/programs/programs.mk
	$(V) echo "Generating object dependencies > $@"	
	$(V) echo $(notdir $(shell timeout $(MIDAS_MAKEDEP_TIMEOUT) makedepf90 $(SRC_ALL))) | \
		sed  's/\([^ ]* :\)/\n\1/g' > $@
	$(if $(findstring $(.SHELLSTATUS),0),,$(error Error : makedepf90 timeout: problem with dependency analysis))

##  ### Absolute linking prerequisites
##
##	generate prerequisites for absolutes using `recursiveDep.py`
##
$(DEP_ABS): $(DEP_OBJ)
	$(if $(shell cat $(DEP_OBJ)),,$(error Error : empty $^))
	$(V) echo "Generating executables dependencies > $@"
	$(V) $(SRC_ROOT)/recursiveDep.py $^ '$(notdir $(ABS))'  > $@


##-----------------------------------------------
##	Phony global targets
##
.PHONY: all absolutes objects info install depend 

all: splitobs absolutes	#> compile and link all programs 

absolutes: $(ABS) 

objects: $(OBJ) #> compile all objects

##-- touch targets
##	fix for out-of-tree touch spurious file issue (#444)
TOUCH_TARGET := $(addprefix touch-,$(OBJ) $(ABS) $(DEP_OBJ) $(DEP_ABS)) \
				touch-objects touch-absolutes

.PHONY: $(TOUCH_TARGET)
$(TOUCH_TARGET) :
	make --touch $(subst touch-,,$@) -f $(SRC_ROOT)/Makefile 

touch-%:	#> mark a specific file as up-to-date
touch-objects :	#> mark all objects as up-to-date
touch-absolutes :	#> mark all absolutes as up-to-date


.PHONY: printAbs
printAbs:
	@echo $(ABS)

.PHONY: printObj
printObj:
	@echo $(OBJ)

info:	#> print build configuration
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

install: $(ABS_INSTALL)			#> install compiled programs
ifeq ($(sort $(ABS_COMPILED)), $(sort $(ABS)))
	@echo
	@echo All absolutes have been installed in $(ABS_INSTALL_DIR) from $(CURDIR)
else ifeq ($(ABS_COMPILED),)
	@echo No absolutes compiled in $(CURDIR)
else
	@echo
	@echo Compiled absolutes installed in $(ABS_INSTALL_DIR) from $(CURDIR)
	@for pgm in $(subst .Abs,,$(ABS_COMPILED)); do\
		echo " * $$pgm"; done
endif

depend : $(DEP_OBJ) $(DEP_ABS)	#> generate all dependency files



##-----------------------------------------------
##	SSM packaging
##
.PHONY: ssm ssm_protect

ssm: $(__ssm_target_version)		#> build SSM package



$(__ssm_workdir_base)/$(__ssm_fullname).ssm: $(ABS)
	test -d $(__ssm_target_version) || ssm created -d $(__ssm_target_version)
	mkdir -p $(__ssm_workdir)/.ssm.d
	mkdir -p $(__ssm_workdir)/bin
	cp $(SRC_ROOT)/ssm.d/post-install $(__ssm_workdir)/.ssm.d/
	sed -e "s|NAME|$(MIDAS_SSM_PKGNAME)|" -e "s|VERSION|$(MIDAS_SSM_VERSION)|" \
        -e "s|PLATFORM|$(ORDENV_PLAT)|" -e "s|MAINTAINER|$(MIDAS_SSM_MAINTAINER)|" \
        -e "s|BUILDINFO|$(VERSION)|" -e "s|DESCRIPTION|$(MIDAS_SSM_DESCRIPTION)|" \
        -e "s|SUMMARY|$(MIDAS_SSM_GITREPO)|" $(SRC_ROOT)/ssm.d/control.json > \
        $(__ssm_workdir)/.ssm.d/control.json
	@echo ">>> building the tarball"
	cp $^ $(__ssm_workdir)/bin/
	cd $(__ssm_workdir_base) && tar cvfzh $(notdir $@) $(__ssm_fullname)

$(__ssm_target_version): $(__ssm_workdir_base)/$(__ssm_fullname).ssm $(__ssm_target_version)
	@echo ">>> installing SSM domain"
	ssm install --debug -d $(__ssm_target_version) -f $<
	@echo ">>> publishingi SSM domain"
	ssm publish --debug -d $(__ssm_target_version) -p  $(__ssm_fullname)
	@echo ">>> SSM package published at $(__ssm_target_version)"

ssm_protect: 	#> protect SSM domain
	ssm freezed -d $(__ssm_target_version) \
        && chmod -R ugo-w $(__ssm_target_version)
	@echo ">>> SSM package $(__ssm_target_version) protected"

##
##	END of Main make rules
##
##===================================================================
##===================================================================
##===================================================================

endif
# vim: set noexpandtab noautoindent nolist:
