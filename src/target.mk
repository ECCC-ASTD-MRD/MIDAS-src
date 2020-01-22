##=========================================================
## 	Multi-real-precision/multi-architecture builds
##
##	auto-precision for programs
ifndef MIDAS_COMPILE_REAL_SINGLE
    ifneq ($(filter obsIO.%,$(MAKECMDGOALS)),)
        MIDAS_COMPILE_REAL_SINGLE = true
        MAKECMDGOALS += 'MIDAS_COMPILE_REAL_SINGLE=true'
    endif
    ifneq ($(filter prepcma.%,$(MAKECMDGOALS)),)
        MIDAS_COMPILE_REAL_SINGLE = true
        MAKECMDGOALS += 'MIDAS_COMPILE_REAL_SINGLE=true'
    endif
endif
##  Fix-hack for autocompletion for single precision programs
obsIO.Abs: 
prepcma.Abs:
##
##  **DO NOT CHANGE THESE `_DIR_SUF`**
##      (see http://make.mad-scientist.net/papers/multi-architecture-builds)
##
ifeq ($(MIDAS_COMPILE_REAL_SINGLE), true)
    _DIR_SUF := real4
else
    _DIR_SUF := real8
endif

DIR_TARGET := $(DIR_BLD_ROOT)/$(EC_ARCH)/$(_DIR_SUF)

revnum := $(shell $(VERSION_SCRIPT))

MAKETARGET = $(MAKE) --no-print-directory -C $@ -f $(CURDIR)/Makefile \
				SRC_ROOT=$(CURDIR) VERSION=$(revnum) \
				$(MAKECMDGOALS)

ifdef VERBOSE
    MAKETARGET += VERBOSE=$(VERBOSE)
endif

.PHONY: $(DIR_TARGET)

$(DIR_TARGET):
	+@[ -d $@ ] || mkdir -p $@
	+@$(MAKETARGET)

Makefile : ;

%.mk :: ;

% :: $(DIR_TARGET) ; :

# vim: set noexpandtab noautoindent:
