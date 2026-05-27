#
# timing/Makefile
#

default:

clean:
	rm -f *.o *.a

install: 
	if [ "$(INSTALL_DIR)" = "" ]; then echo "error: INSTALL_DIR not set"; exit 1; fi
	mkdir -p $(INSTALL_DIR)/lib
	cp *.a $(INSTALL_DIR)/lib

lib:
	s.cc -O 2 -c -openmp hpco_fort_hpm.c weak*.c rebind.c reporter.c history_trace.c ; \
	ar rv libhpcoperf.a hpco_fort_hpm.o weak*.o rebind.o reporter.o history_trace.o
