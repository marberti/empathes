# Empathes - Extensible Minimum Path Estimator

Empathes is an implementation of the CI-NEB method for finding transition states.
It works as an interface for external programs (currently Gaussian and Siesta), from which it takes the information it needs, like molecular energies and atomic forces, to get the job done.

## Installation

First clone the repository

    git clone https://github.com/marberti/empathes.git

Then enter the directory

    cd empathes

Before compiling, make sure that a C compiler, a Fortran compiler and the _make_ utility are installed on your system.
The default copilers used in the Makefile are gcc and gfortran.
If you prefer other compilers, please change the Makefile variables CC and FC accordingly.
Compile the serial version of empathes

    make serial

If, instead, you want to compile the parallel version, make sure to have an MPI implementation installed on your system.
The default MPI compilers used in the Makefile are those from OpenMPI.
If you prefer another implementation, please change the Makefile variables MPICC and MPIFC with the appropriate mpicc and mpifort executable.
Compile the parallel version of empathes

    make parallel

To confirm that the build process was successful, run empathes with the -v argument

    ./empathes -v

## References

1. Henkelman et al., Improved tangent estimate in the nudged elastic band method for finding minimum energy paths and saddle points (doi: 10.1063/1.1323224)
2. Henkelman et al., A climbing image nudged elastic band method for finding saddle points and minimum energy paths (doi: 10.1063/1.1329672)
3. Smidstrup et al., Improved initial guess for minimum energy path calculations (doi: 10.1063/1.4878664)
4. Bitzek et al., Structural Relaxation Made Simple (doi: 10.1103/PhysRevLett.97.170201)
