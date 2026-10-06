# Build of the benchmark driver. ALGENCAN points to the Algencan 3.1.x root.
# The Algencan callbacks are internal procedures, which gfortran implements
# with trampolines on the stack; the linker warning about an executable stack
# is therefore expected and must not be silenced with -z noexecstack.
ALGENCAN ?= $(HOME)/algencan-3.1.2
FC      = gfortran
FFLAGS ?= -O3 -w
SORT    = ../../sources/sort.f90

ovo_bench: hyperdual.o models.o sort.o ovo_bench.f90
	$(FC) $(FFLAGS) ovo_bench.f90 hyperdual.o models.o sort.o \
	    -L$(ALGENCAN)/lib -lalgencan -lhsl -llapack -lblas -o $@

hyperdual.o: hyperdual.f90
	$(FC) $(FFLAGS) -c $<

models.o: models.f90 hyperdual.o
	$(FC) $(FFLAGS) -c $<

sort.o: $(SORT)
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f *.o *.mod ovo_bench
