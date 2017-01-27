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
qsub ~/AssignmentMPI/bashScriptSerial
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
qsub ~/AssignmentMPI/bashScriptExplicit
echo "Explicit ran"

#Compile Implicit
mpicc -o implicit implicitScheme.c -mkl -lm
echo "Implicit compiled"

#Run Implicit
qsub ~/AssignmentMPI/bashScriptImplicit
echo "Implicit ran"

#Compile Crank-Nicolson
mpicxx -o crank crankScheme.cpp -lm
echo "Crank-Nicolson compiled"

#Run Crank-Nicolson
qsub ~/AssignmentMPI/bashScriptCrank
echo "Crank-Nicolson ran"
