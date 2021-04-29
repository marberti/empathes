FC     = gfortran
FLAGS  = -g -cpp -O2 -Wall -Wunused -Wno-maybe-uninitialized

.PHONY: default
default: clean
	$(FC) $(FLAGS) -c blas.f
	$(FC) $(FLAGS) -c lbfgsb.f
	$(FC) $(FLAGS) -c linpack.f
	$(FC) $(FLAGS) -c timer.f
	ar cr liblbfgsb.a blas.o lbfgsb.o linpack.o timer.o

.PHONY: clean
clean:
	@printf "Cleaning liblbfgsb..."
	@rm -f *.o *.a
	@printf " DONE\n"

