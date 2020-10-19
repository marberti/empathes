FCC    = gfortran
CFLAGS = -cpp
EXE    = get_minima.x

.PHONY: default
default:
	$(FCC) $(CFLAGS) -c mod_utility.f90
	$(FCC) $(CFLAGS) -c main.f90
	$(FCC) $(CFLAGS) mod_utility.o main.o -o $(EXE)

.PHONY: clean
clean:
	rm -f *.o *.mod
