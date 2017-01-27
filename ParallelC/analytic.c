# include "mpi.h"
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include <limits.h>

static const int NB_POINTS=1000;

int main ( int argc, char *argv[] )
{
int npes=0,
    myrank=0,
    sizeGrid,
    j;

double x=0.00,
    deltaX=0.00,
	u=0.00,
	deltaT=0.00, 
	previousValue=0.00,
	*subAnal=NULL,
	*resAnal=NULL;

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

  x=((sizeGrid/npes)*myrank*deltaX)-50;

  subAnal=malloc((sizeGrid/npes)*sizeof(double));
  resAnal=malloc(sizeGrid*sizeof(double));

//ANALYTICAL
//Work is divided by processors
		for(j=0;j<sizeGrid/npes;j++)
		{
				subAnal[j]=(0.5)*exp(-pow(x-1.75*(4*deltaT),2));
				x+=deltaX;
		}
//Data gathered in P0	
 MPI_Gather(subAnal,sizeGrid/npes,MPI_DOUBLE,resAnal,sizeGrid/npes,MPI_DOUBLE,0, MPI_COMM_WORLD);
 
//if data missing because of the number of proc
	if(sizeGrid%npes!=0)
	{
		 if(myrank==npes-1)
		{
		  MPI_Send(subAnal+(sizeGrid/npes)-sizeGrid%npes,sizeGrid%npes,MPI_DOUBLE,0,INT_MAX-1, MPI_COMM_WORLD);
		}
		if(myrank==0)
		{
			MPI_Recv(resAnal+(sizeGrid/npes)-sizeGrid%npes,sizeGrid%npes , MPI_DOUBLE,npes-1 , INT_MAX-1, MPI_COMM_WORLD, &status);
		}
	}

if(myrank==0){	
//WRITE ANALYTICAL
FILE *fp1;
fp1=fopen("../Explicit.csv","a");
fprintf(fp1,"Analytical");
fprintf(fp1,";");
for(j=0;j<sizeGrid;j++)
{
	fprintf(fp1,"%0.3lf",resAnal[j]);
	fprintf(fp1,";");
}
fprintf(fp1,"\n");
fclose(fp1);

FILE *fp2;
fp2=fopen("../Implicit.csv","a");
fprintf(fp2,"Analytical");
fprintf(fp2,";");
for(j=0;j<sizeGrid;j++)
{
	fprintf(fp2,"%0.3lf",resAnal[j]);
	fprintf(fp2,";");
}
fprintf(fp2,"\n");
fclose(fp2);

FILE *fp3;
fp3=fopen("../Crank_Nicolson.csv","a");
fprintf(fp3,"Analytical");
fprintf(fp3,";");
for(j=0;j<sizeGrid;j++)
{
	fprintf(fp3,"%0.3lf",resAnal[j]);
	fprintf(fp3,";");
}
fprintf(fp3,"\n");
fclose(fp3);
}
/*
  Terminate MPI.
*/
  MPI_Finalize();

/*
  Free memory.
*/

free(subAnal);
free(resAnal);

  return 0;
}
