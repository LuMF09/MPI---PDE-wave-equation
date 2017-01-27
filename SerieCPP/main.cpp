//main.cpp
/*
Author: Lucien MAMAN
Date: 05/12/2016

Description: Main of the project
*/
#include "mpi.h"
#include "Schemes.h"
#include "Errors.h"
#include <iostream>

static const int NB_POINTS=1000;
static const double TIME_SIMULATION=5.00;

int main(int argc,char* argv[])
{
	MPI_Init(&argc, &argv);

	double step(TIME_SIMULATION), cfl(0.8);
	int grid(NB_POINTS);

	try
	{
		ExplicitUpWind explicitMPI(step, grid, cfl);
		explicitMPI.showScheme();

		ImplicitUpWind implicitMPI(step, grid, cfl);
		implicitMPI.showScheme();

		Crank_Nicolson crankMPI(step, grid, cfl);
		crankMPI.showScheme();
	}

	catch (std::exception &e)
	{
		std::cerr << e.what() << std::endl;
	}
	MPI_Finalize();
		return 0;
}
