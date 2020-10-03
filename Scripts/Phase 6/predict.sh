#!/bin/bash
echo "Generating reference fingerprints..."
gfortran subroutines.f90 readfiles.f90 sym2rcov.f90 fp_simplex_gen.f90 -llapack -lblas -g
./a.out
echo "Generating simplex corners..."
gfortran subroutines.f90 readfiles.f90 sym2rcov.f90 ls.f90 -llapack -lblas
./a.out
echo "Compiling and running fp test generation..."
gfortran subroutines.f90 readfiles.f90 sym2rcov.f90 fp_test_gen.f90 -llapack -lblas -g
./a.out
echo "Compiling and running simplex comparison..."
gfortran simplex_comparison.f90
./a.out
