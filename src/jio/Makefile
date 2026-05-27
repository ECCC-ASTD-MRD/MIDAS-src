#
# jio/Makefile
#

default:

clean:
	rm -f *.a *.o *.so

install:
	if [ "$(INSTALL_DIR)" = "" ]; then echo "error: INSTALL_DIR not set"; exit 1; fi
	mkdir -p $(INSTALL_DIR)/lib
	cp *.a *.so $(INSTALL_DIR)/lib

lib:	shared static

shared:
	s.cc -shared -o s.jio-prof.so -openmp jio-uni-2.c -lm || exit 0

static:
	ar rcv s.jio-prof.a s.jio-prof.so
