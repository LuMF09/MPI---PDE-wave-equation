# include "mpi.h"
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
    sizeSubGrid,
    tsimu,
    t,
    s,
    i,
    j;

double x=0.00,
    deltaX=0.00,
	u=0.00,
	deltaT=0.00, 
	cfl=0.00,
	previousValue=0.00,
	**subNum=NULL,
        *res=NULL,
	tN1=0.00,
	tN2=0.00;

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

  tsimu=TIME_SIMULATION;

//each proc knows the position in x
  x=((sizeGrid/npes)*myrank*deltaX)-50;

  subNum=malloc(tsimu*sizeof(double *));
  res=malloc(sizeGrid*sizeof(double));

//ALlow the right size (it can changes with modulo)
	if(myrank<npes-1)
	{
	 sizeSubGrid =sizeGrid/npes;
	} else sizeSubGrid=sizeGrid/npes + sizeGrid%npes;

 for(i=0;i<tsimu;i++)
{
//we allocate a good memory space to subNum
	subNum[i]=malloc(sizeSubGrid*sizeof(double));
}

tN1=MPI_Wtime();
	//NUMERICAL
	  x=((sizeGrid/npes)*myrank*deltaX)-50;

	for(i=0;i<tsimu;i++)
	{
		if(myrank==0)
		{			
			for(j=1;j<sizeSubGrid;j++)
			{
				if(i==0)
				{
					subNum[0][j] = (0.5)*exp(-pow(x,2));
				
				}
				else
				{
					 subNum[i][j] = subNum[i-1][j] - cfl*(subNum[i-1][j]-subNum[i-1][j- 1]);
				}
			x+=deltaX;

			}
			//SEnd value to the next proc
		MPI_Send(&subNum[i][sizeSubGrid - 1], 1, MPI_DOUBLE, myrank + 1, (i+1)*(myrank+1), MPI_COMM_WORLD);
		}
		else
		{
		MPI_Recv(&previousValue, 1, MPI_DOUBLE, (myrank - 1), (i+1)*(myrank), MPI_COMM_WORLD, &status);
		}
		for(j=0;j<sizeSubGrid;j++)
		{
			if(i==0)
			{
				subNum[0][j] = (0.5)*exp(-pow(x,2));

			}
			else 
			{
				//calculate numerical solution
				subNum[i][j] = subNum[i-1][j] - cfl*(subNum[i-1][j]-previousValue);
				previousValue=subNum[i-1][j];
				
			}
		x+=deltaX;
		}

		if(myrank<npes-1)
		{
			MPI_Send(&subNum[i][sizeSubGrid - 1], 1, MPI_DOUBLE, myrank + 1, (i+1)*(myrank+1), MPI_COMM_WORLD);
		}

		if(myrank==npes-1)
		{//Boudaries
			subNum[0][sizeSubGrid-1]=0;
			subNum[i][sizeSubGrid-1]=0;
		}
	}

 MPI_Gather(subNum[tsimu-1],sizeGrid/npes,MPI_DOUBLE,res,sizeGrid/npes,MPI_DOUBLE,0, MPI_COMM_WORLD);

	if(sizeGrid%npes!=0)
	{
		 if(myrank==npes-1)
		{
		  MPI_Send(subNum[tsimu-1]+sizeSubGrid-sizeGrid%npes,sizeGrid%npes,MPI_DOUBLE,0,INT_MAX, MPI_COMM_WORLD);
		}
		if(myrank==0)
		{
			MPI_Recv(res+sizeSubGrid-sizeGrid%npes,sizeGrid%npes , MPI_DOUBLE,npes-1 , INT_MAX, MPI_COMM_WORLD, &status);
		}
	}

	//end of the time
	tN2=MPI_Wtime();

if(myrank==0){
//WRITE NUMERICAL
FILE *fp;
fp=fopen("../Explicit.csv","a");

fprintf(fp,"Parallel Numerical");
fprintf(fp,";");
for(i=0;i<sizeGrid;i++)
{
	fprintf(fp,"%0.3lf",res[i]);
	fprintf(fp,";");
}
fprintf(fp,"\n");
fprintf(fp,"Elapsed time (parallel):; %f\n",tN2-tN1);
fclose(fp);
}
/*
  Terminate MPI.
*/
  MPI_Finalize();

/*
  Free memory.
*/
 for(i=0;i<tsimu;i++)
{
	free(subNum[i]);
}
free(subNum);
free(res);
  return 0;
}
