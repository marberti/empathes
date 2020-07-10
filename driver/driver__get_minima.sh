gfortran -c ../mod_utility.f90
gfortran -c driver__get_minima.f90
gfortran mod_utility.o driver__get_minima.o -o driver__get_minima.x
