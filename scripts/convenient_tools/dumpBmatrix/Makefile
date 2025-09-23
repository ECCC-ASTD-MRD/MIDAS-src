FC     := s.f90
LIBRMN := rmn

ifeq (, $(shell which s.f90))
  $(error "s.f90 not in PATH: source ./config.dot.sh")
endif


dumpBmatrix: dumpBmatrix.f90
	$(FC) $< -verbose -l$(LIBRMN) -o $@

clean:
	rm dumpBmatrix
