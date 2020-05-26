##=========================================================
## 	Multi-architecture builds
##
##      (see http://make.mad-scientist.net/papers/multi-architecture-builds)
##


revnum := $(shell $(VERSION_SCRIPT))

DIR_TARGET := $(DIR_BLD_ROOT)/$(revnum)/$(EC_ARCH)

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

# vim: set noexpandtab noautoindent nolist:
