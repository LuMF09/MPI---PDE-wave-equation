#include "mpi.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <stdexcept>

static const int NB_POINTS=1000;
static const int TIME_SIMULATION=5;

/* LU factorisation routine
* Takes in a matrix a of size n and produces the lower (l) and
* upper (u) triangular matrices that factorise a
*/

void lu_fact(double *a, double **l, double **u, int n);
void lu_solve(double **l, double **u, double **b, int n, double **x);

int main(int argc,char *argv[]) {

	MPI_Init(&argc,&argv);
	int npes = 0, 
	    myrank = 0,
	    tsimu=TIME_SIMULATION;

	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &npes);


	int sizeGrid = NB_POINTS;

	double t1=0.00, t2=0.00;

	double *init = new double [sizeGrid];
	double *res = new double [sizeGrid];

	double *matrix = new double [sizeGrid * sizeGrid];
	double *subMatrix = new double [(sizeGrid/npes) * sizeGrid];

	double *l = new double [sizeGrid * sizeGrid];
	double *u  = new double [sizeGrid * sizeGrid];

	t1=MPI_Wtime();
//INITIALIZATION of L and U, divided by each proc
	for(int i=0;i<sizeGrid*sizeGrid;i++)
	{
		l[i]=0.00;
		u[i]=0.00;
	}

	double  U = 1.75,
		deltaX = 100.00 / (sizeGrid - 1.00),
		deltaT = (0.8*deltaX)/U,
		cfl=U*deltaT/deltaX,
		a=cfl/(2*deltaX);

	// Matrix initialisation divided by each proc

	for (unsigned int i = 0; i < sizeGrid ; i++) {
		for (unsigned int j = 0; j < (sizeGrid/npes); j++) {
			if ((myrank * sizeGrid/npes) + j - 1  ==  i) {
				subMatrix[j * sizeGrid + i] = -a;
			}else if (i == (myrank * sizeGrid/npes) + j) {
				subMatrix[j * sizeGrid + i] = 1.0 +  2*a;
			} else	if ((myrank * sizeGrid/npes) + j + 1 ==  i) {
				subMatrix[j * sizeGrid + i] = -a;
			} else subMatrix[j*sizeGrid + i]=0.00;	
		}		
	}
	// Gathering all the subMatrix in P0
	MPI_Gather(subMatrix,(sizeGrid/npes)* sizeGrid,MPI_DOUBLE,matrix,(sizeGrid/npes)*sizeGrid,MPI_DOUBLE,0, MPI_COMM_WORLD);
	
	if(myrank == 0)
	{
		if(sizeGrid%npes!=0)
		{
			MPI_Recv(matrix+sizeGrid-sizeGrid%npes,sizeGrid%npes , MPI_DOUBLE,npes-1 , 1, MPI_COMM_WORLD, &status);
		}

		//Computing the Matrix L and U for the LU decomposition

		for (int i = 0; i < sizeGrid; i++) 
		{
			init[i] = (0.5)*exp(-pow((i-50)*deltaX, 2));
		}
		try {
			lu_fact(matrix, &l, &u, sizeGrid);
		}
		catch (std::out_of_range &e) {
			std::cerr << e.what() << std::endl;
		}

		for(int i=0; i<tsimu;i++)
		{		
			lu_solve(&l,&u,&init,sizeGrid,&res);
			
			for (int j = 0; j < sizeGrid; j++) {
				init[j]=res[j];
			}
		}
		t2=MPI_Wtime();

		//WRITE NUMERICAL
		FILE *fp;
		fp=fopen("../Crank_Nicolson.csv","a");

		fprintf(fp,"Parallel Numerical");
		fprintf(fp,";");
		for(int i=0;i<sizeGrid;i++)
		{	
			fprintf(fp,"%0.3lf ",res[i]);
			fprintf(fp,";");
		}
		fprintf(fp,"\nElapsed time (parallel):;%f\n",t2-t1);
		fclose(fp);
	}

	MPI_Finalize();

	delete[]init;
	delete[]matrix;
	delete[]l;
	delete[]u;
	delete[]subMatrix;
	delete[]res;

	return 0;
}

void lu_fact(double *a, double **l, double **u, int n)
{
	double mult;
	double *temp = new double [n*n];

	for(int i = 0;i<n*n;i++){
		temp[i] = (a[i]);
	}

	// LU decomposition without pivoting
	for (int k = 0; k < n - 1; k++) {
		if (fabs(temp[k * n + k]) < 1.e-07) {
			throw std::out_of_range("zero pivot found");
		}
			mult = temp[(k+1)*n + k]/ temp[k*n + k];
			temp[(k+1)*n + k] = mult; 
	// entries of L are saved in temp
				
			temp[(k+1)*n+k+1]-= mult*temp[k*n+k+1];      // entries of U are saved in temp
				if (fabs(temp[(k+1)*n+k+1]) < 1.e-07) {
					throw std::out_of_range("zero pivot found");
				}
	}
	// create l and u from temp
	for (int i = 0; i<n; i++) {
		(*l)[i*n +i] = 1.0;

	}

	for (int i = 1; i<n; i++) {
		for(int j=0;j<i;j++)
			{
				(*l)[i*n + j] = temp[i*n+j];
			}
	}

	for (int i = 0; i<n; i++) {
		for (int j = i; j<n; j++) {
			(*u)[i*n+j] = temp[i*n+j];
		}
	}

	delete[]temp;
}


/*
* Solves the equation LUx = b by performing forward and backward
* substitution. Output is the solution vector x
*/
void lu_solve(double **l, double **u, double **b, int n, double **x)
{

	double *temp = new double [n];
	for(int i = 0;i<n;i++){
		temp[i] = (*b)[i];
	}

	// forward substitution for L y = b.
	for (int i = 1; i < n; i++) {
		for (int j = 0; j < i; j++) {
			temp[i] -= (*l)[i*n+j] * temp[j];
		}
	}

	// back substitution for U x = y.  
	for (int i = n - 1; i >= 0; i--) {
		for (int j = i + 1; j < n; j++) {
			temp[i] -= (*u)[i*n+j] * temp[j];
		}
		temp[i] /= (*u)[i*n+i];
	}

	// copy solution into x
	for (int i = 0; i<n; i++) {
		(*x)[i] = temp[i];
	}

	delete[]temp;
}
