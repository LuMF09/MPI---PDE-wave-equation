# include "mpi.h"
#include <mkl.h>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include <limits.h>

static const int NB_POINTS=1000;
static const int TIME_SIMULATION=5;

int main ( int argc, char *argv[] )
{
int npes=0,
    myrank=0,
    sizeGrid,
    i,j, 
    tsimu=TIME_SIMULATION;

double x=0.00,
       deltaX=0.00,
       u=0.00,
       deltaT=0.00, 
       cfl=0.00,
      *coefD=NULL,
      *subCoefD=NULL,
      *res=NULL,
      *subRes=NULL,
      *temp=NULL,
       diago=0.00,
       subDiago=0.00,
       m=0.00;

float *diagoLAP=NULL,
      *upLAP=NULL,
      *lowLAP=NULL,
      *initLAP=NULL;

double t1=0.00,t2=0.00;

/* 
  Initialize MPI.
*/
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD,&npes);

  sizeGrid=NB_POINTS;
  deltaX=100.00/sizeGrid;
  u=1.75;
  deltaT=(deltaX*0.8)/u;
  cfl=u*(deltaT/deltaX);
  diago=1+cfl;

//each proc know the position in x
  x=((sizeGrid/npes)*myrank*deltaX)-50;

  subCoefD=malloc(sizeGrid/npes*sizeof(double));
  res=malloc(sizeGrid*sizeof(double));
  subRes=malloc((sizeGrid/npes)*sizeof(double));

//beginning of the time
	t1=MPI_Wtime();

	if(myrank==0)
	{
		x=-50;
		coefD=malloc(sizeGrid*sizeof(double));
		subDiago=-cfl;
		m=subDiago/diago;
		coefD[0]=0.00;

//Calculation of all d coefficients
		for(i=1;i<sizeGrid;i++)
		{
			coefD[i]=(0.5)*exp(-pow(x,2));
			x+=deltaX;
			coefD[i]-=m*coefD[i-1];
		}
	}
	//t=5
	for(i=0;i<tsimu;i++)
	{
		if(myrank==0)
		{
			if(sizeGrid%npes!=0)
			{
				MPI_Send(coefD+sizeGrid-sizeGrid%npes,sizeGrid%npes,MPI_DOUBLE,npes-1,2, MPI_COMM_WORLD);
			}
		}

			//SCATTER the d updated
		MPI_Scatter(coefD,sizeGrid/npes, MPI_DOUBLE,subCoefD, sizeGrid/npes, MPI_DOUBLE,0, MPI_COMM_WORLD);

		for(i=0;i<sizeGrid/npes;i++)
		{
			subRes[i]=subCoefD[i]/diago;
		}
		//GATHER sub results into the final result
		MPI_Gather(subRes,sizeGrid/npes,MPI_DOUBLE,res,sizeGrid/npes,MPI_DOUBLE,0, MPI_COMM_WORLD);

//If the number of proc does not divide the number of points
		if(sizeGrid%npes!=0)
		{
			 if(myrank==npes-1)
			{
				coefD=malloc(sizeGrid%npes*sizeof(double));
				MPI_Recv(coefD,sizeGrid%npes , MPI_DOUBLE,0, 2, MPI_COMM_WORLD, &status);

			  	temp=malloc((sizeGrid%npes)*sizeof(double));

				for(i=0;i<sizeGrid%npes;i++)
				{
					temp[i]=coefD[i]/diago;
				}
				MPI_Send(temp,sizeGrid%npes,MPI_DOUBLE,0,3, MPI_COMM_WORLD);
			}
			if(myrank==0)
			{
				MPI_Recv(res+sizeGrid-sizeGrid%npes,sizeGrid%npes , MPI_DOUBLE,npes-1 , 3, MPI_COMM_WORLD, &status);
			}
		}
		if(myrank==0)
		{
			for(j=0;j<sizeGrid;j++)
			{
				coefD[j]=res[j];
			}
		}
	}


t2=MPI_Wtime();

	if(myrank==0)
	{
		//WRITE NUMERICAL
		FILE *fp;
		fp=fopen("../Implicit.csv","a");

		fprintf(fp,"Parallel Numerical");
		fprintf(fp,";");
		for(i=0;i<sizeGrid;i++)
		{
			fprintf(fp,"%0.3lf ",res[i]);
			fprintf(fp,";");
		}
		fprintf(fp,"\nElapsed time (parallel):;%f\n",t2-t1);
		fclose(fp);
	}



//WITH LAPACK

	if(myrank==0)
	{
		x=-50;
		diagoLAP=malloc(sizeGrid*sizeof(double));
		upLAP=malloc((sizeGrid-1)*sizeof(double));
		lowLAP=malloc((sizeGrid-1)*sizeof(double));
		initLAP=malloc(sizeGrid*sizeof(double));

		for(i=0;i<sizeGrid-1;i++)
		{
			diagoLAP[i]=1+cfl;
			upLAP[i]=0;
			lowLAP[i]=-cfl;

			initLAP[i]=(0.5)*exp(-pow(x,2));
			x+=deltaX;
		}
		diagoLAP[sizeGrid-1]=1+cfl;
		initLAP[sizeGrid-1]=0;
		int leadB = sizeGrid;
		int info = 0;
		int colB=1;
		t1=0;
		t2=0;
		t1=MPI_Wtime();
		for(i=0;i<tsimu;i++)
		{
			//LAPACK function
			SGTSV(&sizeGrid,&colB,lowLAP,diagoLAP,upLAP,initLAP,&leadB, &info);
		}
		t2=MPI_Wtime();
		if(info==0)
		{
			FILE *fp1;
			fp1=fopen("../Implicit.csv","a");
			fprintf(fp1,"Parallel Numerical LAPACK");
			fprintf(fp1,";");
			for(i=0;i<sizeGrid;i++)
			{
				fprintf(fp1,"%0.3lf ",initLAP[i]);
				fprintf(fp1,";");
			}
			fprintf(fp1,"\n");
			fprintf(fp1,"\nElapsed time (parallel LAPACK):;%f\n",t2-t1);
			fclose(fp1);
		}
	}

/*
  Terminate MPI.
*/
  MPI_Finalize();

/*
  Free memory.
*/

free(coefD);
free(subCoefD);
free(res);
free(subRes);
free(temp);
free(diagoLAP);
free(upLAP);
free(lowLAP);
free(initLAP);

  return 0;
}
