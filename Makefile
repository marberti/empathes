CC       = gcc
FC       = gfortran
MPICC    = mpicc.openmpi
MPIFC    = mpif90.openmpi
FLAGS    = -g -cpp -O2 -Wall -Wunused -Wpedantic -Wno-maybe-uninitialized
CFLAGS   = -std=c99
FFLAGS   = -std=f2008
LPATH    = -L./lib
LIBS     = -llbfgsb

CSOURCES = c_utility.c

FSOURCES = mod_utility.f90           \
	   mod_c_utility_wrapper.f90 \
	   mod_bfgs_wrapper.f90      \
	   mod_rotation.f90          \
	   mod_geometry.f90          \
	   mod_idpp.f90              \
	   mod_pes_data.f90          \
	   mod_pes.f90               \
	   mod_elastic.f90           \
	   mod_climbing.f90          \
	   mod_output.f90            \
	   mod_optimization.f90      \
	   mod_slave.f90             \
	   mod_computation_info.f90  \
	   mod_input.f90             \
	   main.f90

COBJECTS = $(CSOURCES:.c=.o)
FOBJECTS = $(FSOURCES:.f90=.o)
OUT      = neb.x

# main compilation options --------------------------------
.PHONY: default
default: help

.PHONY: debug0
debug0: FLAGS += -DDBG0
debug0: serial

.PHONY: serial
serial: clean $(OUT)

.PHONY: parallel
parallel: CC = $(MPICC)
parallel: FC = $(MPIFC)
parallel: FLAGS += -DUSE_MPI
parallel: serial

# utility -------------------------------------------------
.PHONY: help
help:
	@echo "Usage:"
	@echo
	@echo "  make serial          Serial compilation"
	@echo "  make parallel        Parallel compilation (MPI)"

.PHONY: clean
clean:
	@printf "Cleaning..."
	@rm -f *.o *.mod
	@printf " DONE\n"

.PHONY: screenclear
screenclear:
	@clear

# core ----------------------------------------------------
$(OUT): allobjects
	$(FC) $(LPATH) $(FLAGS) $(FFLAGS) $(COBJECTS) $(FOBJECTS) $(LIBS) -o $(OUT)

.PHONY: allobjects
allobjects: $(COBJECTS) $(FOBJECTS)

$(FOBJECTS): %.o: %.f90
	$(FC) $(FLAGS) $(FFLAGS) -c $<

$(COBJECTS): %.o: %.c
	$(CC) $(FLAGS) $(CFLAGS) -c $<

