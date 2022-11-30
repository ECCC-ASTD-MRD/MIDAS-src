
FC     := s.f90
LIBRMN := rmnMP


FC := s.f90
ifeq (, $(shell which s.f90))
  $(error "s.f90 not in PATH: source ./config.dot.sh")
endif

dumpBmatrix: dumpBmatrix.f90
	$(FC) $< -l$(LIBRMN) -o $@

clean:
	rm dumpBmatrix
