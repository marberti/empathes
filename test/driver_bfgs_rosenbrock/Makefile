FC    = gfortran
FLAGS = -cpp
EXE   = driver.x

.PHONY: default
default: clean link
	$(FC) $(FLAGS) -c mod_utility.f90
	$(FC) $(FLAGS) -c mod_bfgs.f90
	$(FC) $(FLAGS) -c main.f90
	$(FC) $(FLAGS) -o $(EXE) mod_utility.o mod_bfgs.o main.o

.PHONY: link
link:
	ln -s ../../mod_utility.f90 .
	ln -s ../../mod_bfgs.f90 .

.PHONY: clean
clean:
	@printf "Cleaning... "
	@rm -f *.o *.mod mod_*.f90
	@printf "DONE\n"

