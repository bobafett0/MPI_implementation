
#include <stdio.h>
#include <stdlib.h>
#include "derv.h"
#include <math.h>
#include <string.h>
#include <assert.h>
#include "derv.h"
#include <mpi.h>
#include <stdbool.h>


#define LEN 10
#define timeStep 5
#define deltaX 0.1
#define deltaY 0.1
#define deltaTI 0.01
#define TOL 0.00001

int main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);
  
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  printf("the size is : %d\n",size);
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  printf("the rank is %d\n", rank);	

  int rangeS = LEN/size;
  int rangeB = rangeS;
  if (rank == size-1)
  {
	rangeS = rangeS + LEN % size;
  }
  
  double* x = calloc(rangeS,sizeof(double));
  double* y = calloc(LEN,sizeof(double));

  double**f = alloc(rangeS,LEN);
	
  // initializing file to print to
  FILE *file = fopen("DATA.txt", "w");
  assert(file);
     
  // Initializing all x-y points
    
  for(int i = rank*(rangeB); i < (rank+1)*(rangeS); i++)
  {
	x[i-rank*(LEN/size)] = (i*deltaX);
  }
  
  for(int i = 0; i < LEN; i++)
  {
	y[i] = (i*deltaY);
  }
  	
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 // printf("the rank is %d\n", rank);	
  
  int counter = 0;
  // Initial value  - (X^2)*(Y^2)

  for(int i = 0; i < (rangeS); i++)
  { 
    for(int j = 0; j < LEN; j++)
	{
	  f[i][j] = x[i]*x[i]*y[j]*y[j] + 5.0;	
	  //printf("rank is %d, currently on %.10f , %d , %d\n", rank,f[i][j],i,j );
	  counter++;
	}
  }
  //free(x);free(y);
  //MPI_Request req[2];
  //MPI_Waitall(size,req,MPI_STATUSES_IGNORE );
  
  //printTo(file,f,LEN/size,LEN,rank);
  printf("%d values were initialized\n",counter );
  
  printf("the rank is %d\n",rank);
  printOut(f,rangeS,LEN);
  
  double step = deltaTI;
  double* ste = &step;
  double timePassed = 0.0;
  double** temp;
  for(int i = 0; i < timeStep; i++)
  {
	//printf("before call %d rank is %d \n",i,rank);
	temp = kutta(f,deltaX,TOL,rangeS,LEN,ste);

	for(int i = 1; i < rangeS; i++)
	{
	  free(f[i]);	
	}
	
	if (rank == 0 )
	{
	  free(f[0]);free(f);
    }
	if (rank == 2)
	  free(f[0]);
	
	f = temp;
	timePassed += *ste;
  }
  printf("Finished, the rank is %d\n",rank);
  printOut(f,LEN/rangeS,LEN);
  free(temp);
    printf("\n\n\nxfinished \n\n\n" );

  
  MPI_Finalize();
  //freePointers(f,LEN/size);
  
  
}
