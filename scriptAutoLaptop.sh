#!/bin/bash

#set up the environment
. /usr/share/modules/init/bash
module load icc
module load impi

. mpivars.sh
. iccvars.sh intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
/usr/local/commercial/intel/xe2013/composerxe/mkl/bin/mklvars.sh intel64

cd ~/AssignmentMPI/SerieCPP/

#Compile all Serial
mpicxx -c *.cpp -lm
echo "Serial compiled"

#Link all Serial
mpicxx -o main *.o -lm
echo "Serial linked"

#Run main Serial
mpirun -n 1 ./main 
echo "Serial ran"

cd ../ParallelC/
#Compile Analytical
mpicc -o analytic analytic.c -lm
echo "Analytical compiled"

#Run Analytical on the computer to not use credits on Astral
mpirun -n 5 ./analytic
echo "Analytical ran"

#Compile Explicit
mpicc -o explicit explicitScheme.c -lm
echo "Explicit compiled"

#Run Explicit
mpirun -n 7 ./explicit
echo "Explicit ran"

#Compile Implicit
mpicc -o implicit implicitScheme.c -mkl -lm
echo "Implicit compiled"

#Run Implicit
mpirun -n 4 ./implicit
echo "Implicit ran"

#Compile Crank-Nicolson
mpicc -o crank crankScheme.cpp -lm
echo "Crank-Nicolson compiled"

#Run Crank-Nicolson
mpirun -n 5 ./crank
echo "Crank-Nicolson ran"
