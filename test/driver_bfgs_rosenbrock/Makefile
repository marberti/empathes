# Copyright (C) 2020-2021  Marco Bertini
#
# This file is part of Empathes.
#
# Empathes is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Empathes is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Empathes.  If not, see <https://www.gnu.org/licenses/>.

FC    = gfortran
FLAGS = -cpp
EXE   = driver.x

.PHONY: default
default: clean link
	$(FC) $(FLAGS) -c mod_utility.f90
	$(FC) $(FLAGS) -c mod_bfgs.f90
	$(FC) $(FLAGS) -c mod_lbfgs.f90
	$(FC) $(FLAGS) -c driver.f90
	$(FC) $(FLAGS) -o $(EXE) mod_utility.o mod_bfgs.o mod_lbfgs.o driver.o

.PHONY: link
link:
	ln -s ../../mod_utility.f90 .
	ln -s ../../mod_bfgs.f90 .
	ln -s ../../mod_lbfgs.f90 .

.PHONY: clean
clean:
	@printf "Cleaning... "
	@rm -f *.o *.mod mod_*.f90
	@printf "DONE\n"

