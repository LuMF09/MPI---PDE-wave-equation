//Schemes.cpp
/*
Author: Lucien MAMAN
Date: 05/12/2016

Description: Implementation of all classes of Schemes.h file:
abstract Schemes,
Upwind Explicit, Upwind Explicit, Lax-Wendroff, Crank-Nicolson all inherited from Schemes
*/
#include "mpi.h"
#include "Schemes.h"
#include "Maths_Functions.h"
#include "Errors.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <time.h>
#include <stdlib.h>
#include <stdexcept>


void writeResult(const std::string &scheme, std::vector<double> &x, std::vector<double> &res, double time);
void lu_fact(std::vector<std::vector<double> > &a, std::vector<std::vector<double> >& l, std::vector<std::vector<double> >& u, int n);
void lu_solve(std::vector<std::vector<double> >& l, std::vector<std::vector<double> >& u, std::vector<double> &b, int n, std::vector<double> &x);

//FUNCTION TO WRITE IN THE CORRECT FILE
void writeResult(const std::string &scheme, std::vector<double> &x, std::vector<double> &res, double time)
{
	std::string fileName("../"+scheme+".csv");
	std::ofstream file(fileName.c_str(),std::ios::out | std::ios::app);
	if(file)
	{
		//Time
		file<<"Elapsed time (serial):"<<";"<<time<<std::endl;
		
		//Vector x
		file <<";";
		for(int i=0; i<res.size();i++)
		{
			file << x[i] <<";";
		}
		file << std::endl;

		//vector res
		file <<"Serial Numerical" <<";";
		for(int i=0;i<res.size();i++)
		{
			file << res[i] << ";";
		}
		file << std::endl;	
	}else std::cerr <<"File open";
}

/* LU factorisation routine
* Takes in a matrix a of size n and produces the lower (l) and
* upper (u) triangular matrices that factorise a
*/


void lu_fact(std::vector<std::vector<double> > &a, std::vector<std::vector<double> >& l, std::vector<std::vector<double> >& u, int n)
{
	double mult;

	std::vector<std::vector<double> > temp(a);

	// LU decomposition without pivoting
	for (int k = 0; k < n - 1; k++) {
		for (int i = k + 1; i < n; i++) {
			if (fabs(temp[k][k]) < 1.e-07) {
				throw std::out_of_range("zero pivot found");
			}
			mult = temp[i][k] / temp[k][k];
			temp[i][k] = mult;                      // entries of L are saved in temp
			for (int j = k + 1; j < n; j++) {
				temp[i][j] -= mult*temp[k][j];      // entries of U are saved in temp
				if (fabs(temp[i][i]) < 1.e-07) {
					throw std::out_of_range("zero pivot found");
				}
			}
		}
	}

	// create l and u from temp
	for (int i = 0; i<n; i++) {
		l[i][i] = 1.0;
	}

	for (int i = 1; i<n; i++) {
		for (int j = 0; j<i; j++) {
			l[i][j] = temp[i][j];
		}
	}

	for (int i = 0; i<n; i++) {
		for (int j = i; j<n; j++) {
			u[i][j] = temp[i][j];
		}
	}

}


/*
* Solves the equation LUx = b by performing forward and backward
* substitution. Output is the solution vector x
*/
void lu_solve(std::vector<std::vector<double> >& l, std::vector<std::vector<double> >& u, std::vector<double> &b, int n, std::vector<double> &x)
{

	std::vector<double> temp(b);

	// forward substitution for L y = b.
	for (int i = 1; i < n; i++) {
		for (int j = 0; j < i; j++) {
			temp[i] -= l[i][j] * temp[j];
		}
	}
	// back substitution for U x = y.  
	for (int i = n - 1; i >= 0; i--) {
		for (int j = i + 1; j < n; j++) {
			temp[i] -= u[i][j] * temp[j];
		}
		temp[i] /= u[i][i];
	}

	// copy solution into x
	for (int i = 0; i<n; i++) {
		x[i] = temp[i];
	}
}

double const LEFT_BOUNDARY(-50.00), RIGHT_BOUNDARY(50.00), U(1.75); //x belong to [-50;50] and u=1.75

//SCHEMES ABSTRACT CLASS

//CONSTRUCTOR

//Initialization of all protected variables
Schemes::Schemes(double step, int gridSize, double cfl1) : x(gridSize), vANAL(gridSize), vNUM(), step(step), deltaT(0), deltaX(0), cfl(0), nbIt(0)
{
	if (gridSize == 1) //Exception division by 0
	{
		throw DivideZero("in Scheme constructor (see value of gridSize)");
	}

		deltaX = fabs(RIGHT_BOUNDARY - LEFT_BOUNDARY) / (gridSize - 1);

		if (U == 0) //Exception division by 0
		{
			throw DivideZero(" in Scheme constructor (see value of U)");
		}

		deltaT = (cfl1 * deltaX)/U; 

		if (deltaT==0) //Exception division by 0
		{
			throw DivideZero(" in Scheme constructor (see value of delta T)");
		}

		nbIt = (int)round(step/deltaT);//round the value to have the closest optimised CFL

		if (nbIt==0) //Exception division by 0
		{
			throw DivideZero(" in Scheme constructor (see value of nbIt)");
		}
		deltaT = step / nbIt;

		if (deltaX == 0) //Exception division by 0
		{
			throw DivideZero(" in Scheme constructor (see value of deltaX)");
		}
		cfl=U*(deltaT / deltaX);

		for (int i = 0; i < nbIt; i++)
		{
			vNUM.push_back(std::vector<double>(gridSize));//create exactly the good size for vNUM
		}
			

	//Construction of the vector x with all positions

	x[0] = LEFT_BOUNDARY;
	x[(*this).getSizeGrid() - 1] = RIGHT_BOUNDARY;
	
	for (int i = 1; i < (*this).getSizeGrid() - 1; i++)
	{
		x[i] = x[i - 1] + deltaX;	
	}

}

//ACCESSORS

int Schemes::getSizeGrid() const
{
	return x.size();
}
double Schemes::getStep() const
{
	return step;
}

//FUNCTIONS 

void Schemes::calculAnalytical() // Compute and write in a CSV file the analytical solution
{
		//Boundary conditions

			vANAL[0] = 0.00;
			vANAL[(*this).getSizeGrid() - 1] = 0.00;
		

		// Display the last line of the analytical solution in the CSV file Results

			for (int j = 1; j < (*this).getSizeGrid() - 1; j++)
			{
				vANAL[j] = 0.5*exp(-((x[j] - U*(*this).getStep())*(x[j] - U*(nbIt - 1)*deltaT)));
			}
	}

void Schemes::calculNorms() // Compute norm 0 (max error) and norm 2 (square root of the sum of the square of errors)
{
	double norm2(0.00), norm0(0.00), sumCarre(0.00), norm2Normalized(0.00);

		for (int i = 0; i < (*this).getSizeGrid(); i++)
		{
			sumCarre += (vANAL[i] - vNUM[nbIt-1][i])*(vANAL[i] - vNUM[nbIt-1][i]); // Sum of (error)^2

			if (norm0 < fabs(vANAL[i] - vNUM[nbIt-1][i])) // compare if value of previous max is higher
			{
				norm0 = fabs(vANAL[i] - vNUM[nbIt-1][i]);
			}
		}
		norm2 = sqrt(sumCarre); // Square root of sumCarre
		norm2Normalized = norm2 / (*this).getSizeGrid(); // we divide by the number of point to normalise the norm

}

void Schemes::initNumerical(std::vector<std::vector <double> >&v) // initialize a vector of vector with boundaries depending the choice made by the user
{
	
	for (int t = 0; t < nbIt; t++)
	{
			v[t][0] = 0.00;
			v[t][(*this).getSizeGrid() - 1] = 0.00;
	}

	//Construct v from 1 to v[size of the grid - 2]
	for (int j = 1; j < (*this).getSizeGrid()-1; j++)
	{
			v[0][j] = 0.5*exp(-(x[j] * x[j]));
	}
}

void Schemes::showScheme() {

	//test if the scheme is stable
	(*this).isStable();

	//ANALYTICAL SOLYUTION
	(*this).calculAnalytical();

	//NUMERICAL SOLUTION
	(*this).calculNumerical();

	//NORMS
	(*this).calculNorms();
}

	//EXPLICIT UpWind Scheme

	//CONSTRUCTOR

ExplicitUpWind::ExplicitUpWind(double step, int size, double cfl1) : Schemes(step, size, cfl1) //The user choose f0 or f1, the time to run and the size of the grid
	{
	}

	//FUNCTIONS

void ExplicitUpWind::calculNumerical() //Compute the NUMERICAL solution of the Explicit UPWIND scheme
{

	double t1=0;
	double t2=0;
	t1=MPI_Wtime();

	(*this).initNumerical(vNUM); //Initialize the vNUM vector of vector

		for (int t = 0; t <nbIt-1; t++)
		{
			for (int space = 1; space < (*this).getSizeGrid()-1; space++)
			{
				vNUM[t+1][space] = vNUM[t][space] - cfl*(vNUM[t][space] - vNUM[t][space - 1]); //Construction of vNUM following the scheme
			}
		}
	t2=MPI_Wtime();
	double res = t2-t1;

	writeResult("Explicit", x, vNUM[5],res);
}

void ExplicitUpWind::isStable() //Test if the scheme is stable 
{
	if (cfl>1) // Stability condition
	{
		std::cout << std::endl << "This system is not stable" << std::endl;
	}
}

//IMPLICIT UPWIND SCHEME

//CONSTRUCTOR

ImplicitUpWind::ImplicitUpWind(double step, int size, double cfl1) : Schemes(step, size,cfl1) //The user choose f0 or f1, the time to run and the size of the grid and the way to compute
{
}

//FUNCTIONS

void ImplicitUpWind::calculNumerical() //Compute the NUMERICAL solution of the Implicit UPWIND scheme
{
	double t1=0;
	double t2=0;
	t1=MPI_Wtime();

	(*this).initNumerical(vNUM); //Initialize the vNUM vector of vector

	for (int t = 0; t < nbIt - 1; t++)
	{
		for (int space = 1; space < (*this).getSizeGrid() - 1; space++)
		{
			vNUM[t + 1][space] = (vNUM[t][space]) / (1 + cfl) + (cfl*(vNUM[t + 1][space - 1])) / (1 + cfl); //Construction of vNUM following the scheme
		}
	}

	t2=MPI_Wtime();
	double res = t2-t1;
	writeResult("Implicit", x, vNUM[5],res);
}

void ImplicitUpWind::isStable() //Test if the scheme is stable 
{
}

//Crank-Nicolson SCHEME

//CONSTRUCTOR

Crank_Nicolson::Crank_Nicolson(double step, int size, double cfl1) : Schemes(step, size, cfl1) 
{
}

//FUNCTIONS

void Crank_Nicolson::calculNumerical() //Compute the NUMERICAL solution of the Crank_Nicolson scheme
{
	double t1=0;
	double t2=0;
	t1=MPI_Wtime();

	(*this).initNumerical(vNUM); //Initialize the vNUM vector of vector

	double a = cfl / (2 * deltaX);

	std::vector<std::vector<double> > matrix((*this).getSizeGrid(),std::vector<double>((*this).getSizeGrid(),0)), 
					  l((*this).getSizeGrid(), std::vector<double>((*this).getSizeGrid(),0)), 
					  u((*this).getSizeGrid(), std::vector<double>((*this).getSizeGrid(), 0));
	
	for (int i = 0; i < (*this).getSizeGrid(); i++) {
		for (int j = 0; j < (*this).getSizeGrid(); j++) {
			if (j - 1 == i) {
				matrix[i][j] = -a;
			}
			if (i == j) {
				matrix[i][j] = 1.0 + 2*a;
			}
			if (j + 1 == i) {
				matrix[i][j] = -a;
			}
		}
	}

	try {
		lu_fact(matrix, l, u, (*this).getSizeGrid());
	}
	catch (std::out_of_range &e) {
		std::cerr << e.what() << std::endl;
	}
	
for(int i=0;i<5;i++)
{
	lu_solve(l,u,vNUM[i],(*this).getSizeGrid(),vNUM[i+1]);
}
	t2=MPI_Wtime();
	double res = t2-t1;
	writeResult("Crank_Nicolson", x, vNUM[5],res);
}

void Crank_Nicolson::isStable() //Test if the scheme is stable 
{
}

