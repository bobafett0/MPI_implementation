
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
#define timeStep 2000000
#define deltaX 0.1
#define deltaY 0.1
#define deltaTI 0.01
#define TOL 0.0000001

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

  double**f = alloc(rangeS,LEN);//[LEN/size][LEN];
    
  //MPI_Alloc_mem(sizeof(double*)*(LEN/size),MPI_INFO_NULL,&f); 
  // Initializing all pointers of f
	
  for(int i = rank*(rangeS); i < (rank+1)*(rangeS); i++)
  {
    //MPI_Alloc_mem(LEN*sizeof(double),MPI_INFO_NULL,&f[i] );
    f[i] = calloc(LEN,sizeof(double));
  }
  
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
  printf("the rank is %d\n", rank);	
  
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

  //MPI_Request req[2];
  //MPI_Waitall(size,req,MPI_STATUSES_IGNORE );
  
  //printTo(file,f,LEN/size,LEN,rank);
  printf("%d values were initialized\n",counter );
  
  //for(int i = 0; i < (LEN/size); i++)
  //{
    //for(int j = 0; j < LEN; j++)
	//{
	  //printf("%0.9f ",f[i][j]);	
	//}
	//printf("\n");
  //}
  printf("the rank is %d\n",rank);
  printOut(f,rangeS,LEN);
  
  double step = deltaTI;
  double* ste = &step;
  double timePassed = 0.0;
  double** temp;
  for(int i = 0; i < timeStep; i++)
  {
	//if(i % 1000 == 0)
	//{
	//printTo(file,f,LEN/size,LEN,rank);
	//if(isAlmostUni(f,LEN,LEN,0.001))
	 //break;
	//printf("calling kutta : \n");
	//printf("before call rank is %d \n",rank);
	temp = kutta(f,deltaX,TOL,rangeS,LEN,ste);
	//if(rank == 1){
	//printf("after %d call, rank is %d \n",i,rank );
	//printf("\n\n");
	//printOut(f,LEN/size,LEN);
	//printf("\n\n");
   // }
   // if (rank == 0)
	//freePointers(f,LEN/size);
	for(int i = 1; i < rangeS; i++)
	{
	  free(f[i]);	
	}
	
	//free(f[1]);free(f[2]);free(f[3]);free(f[4]);//free(f[5]);//free(f);
	if(rank == 0)
	{
	  free(f[0]);free(f);
	}
	
	//printf("pointers freed at rank %d \n",rank);
	//free(f[0]);
	//MPI_Free_mem(f);
	
	
	//for(int i = 0; i < LEN/size; i++)
	//{
	//MPI_Request req[2];
       
    //double** recv; //= calloc(LEN,sizeof(double));

    //if ( rank == 0 )
    //{
      //MPI_Isend(&f[i], LEN, MPI_DOUBLE, rank+1, 1235, MPI_COMM_WORLD, &req[0]);
      //MPI_Irecv(recv, LEN, MPI_DOUBLE, rank+1, 1235, MPI_COMM_WORLD, &req[1]);
    //}
    //else if ( rank == size-1 )
    //{
	  //MPI_Irecv(recv, LEN, MPI_DOUBLE, rank-1, 1235, MPI_COMM_WORLD, &req[1]);
      //MPI_Isend(&f[i], LEN, MPI_DOUBLE, rank-1, 1235, MPI_COMM_WORLD, &req[0]);
    //}
    //printf("rank %d waiting \n", rank);
    //MPI_Waitall(size,req,MPI_STATUSES_IGNORE );	
    ////free(recv);//free(f[i]);
    //}
    
	f = temp;
	timePassed += *ste;
	
	//int rank,size;
    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //MPI_Comm_size(MPI_COMM_WORLD, &size);
    //MPI_Request req[2];
  
    //double* recv = calloc(1,sizeof(double));
    //double* send = calloc(1,sizeof(double));
    
    //if ( rank == 0 )
    //{
      //MPI_Isend(send, 1, MPI_DOUBLE, rank+1, 1235, MPI_COMM_WORLD, &req[0]);
      //MPI_Irecv(recv, 1, MPI_DOUBLE, rank+1, 1235, MPI_COMM_WORLD, &req[1]);
    //}
    //else if ( rank == size-1 )
    //{
	  //MPI_Irecv(recv, 1, MPI_DOUBLE, rank-1, 1235, MPI_COMM_WORLD, &req[1]);
      //MPI_Isend(send, 1, MPI_DOUBLE, rank-1, 1235, MPI_COMM_WORLD, &req[0]);
    //}
    //printf("rank %d waiting \n", rank);
    //MPI_Waitall(size,req,MPI_STATUSES_IGNORE );	
    //free(recv);free(send);
    
  }
  printf("the rank is %d\n",rank);
  printOut(f,LEN/2,LEN);
  
  MPI_Finalize();
  printf("\n\n\nxfinished \n\n\n" );
  //freePointers(f,LEN/size);
  free(temp);
  
  
}
