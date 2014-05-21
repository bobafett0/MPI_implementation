

#include <stdio.h>
#include <stdlib.h>
#include "derv.h"
#include <math.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

double** change( double** u, double dx, int sizeX, int sizeY )
{
  // initializing change matrix the same size of U
	
  double** change = alloc(sizeX,sizeY);	
  double*    rec1;
  double*    rec2;
  //double*   send1;
  //double*   send2;
  
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Request req[2];
  
  
  if(rank != 0)
  {
	//send1 = calloc(sizeY,sizeof(double));
	//for (int i = 0; i < sizeY; i++)
    //{  
      //send1[i] = u[0][i];
    //}
	rec1 = calloc(sizeY,sizeof(double));
  }
  if(rank != size-1)
  {
	//send2 = calloc(sizeY,sizeof(double));
	//for (int i = 0; i < sizeY; i++)
    //{  
      //send2[i] = u[sizeX-1][i];
    //}
    rec2 = calloc(sizeY,sizeof(double));
  }
  // sending the messages 
  
  //printf("IN change, the rank is %d\n",rank);

  
  // sending the last rows to the appropriate reciever
  
  if ( rank == 0 )
  {
	MPI_Irecv(rec2, sizeY, MPI_DOUBLE, rank+1, 1244, MPI_COMM_WORLD, &req[1]);
    MPI_Isend(u[sizeX-1], sizeY, MPI_DOUBLE, rank+1, 1244, MPI_COMM_WORLD, &req[0]);
  }
  else if ( rank == size-1 )
  {
	MPI_Irecv(rec1, sizeY, MPI_DOUBLE, rank - 1, 1244, MPI_COMM_WORLD, &req[1]);
    MPI_Isend(u[0], sizeY, MPI_DOUBLE, rank - 1, 1244, MPI_COMM_WORLD, &req[0]);
  }
  else if ( rank != 0 && rank != size -1)
  {
    MPI_Irecv(rec1, sizeY, MPI_DOUBLE, rank-1, 1234, MPI_COMM_WORLD, &req[2]);
    MPI_Isend(u[0], sizeY, MPI_DOUBLE, rank-1, 1234, MPI_COMM_WORLD, &req[3]);
    
    MPI_Irecv(rec2, sizeY, MPI_DOUBLE, rank+1, 1234, MPI_COMM_WORLD, &req[2]);
    MPI_Isend(u[sizeX-1], sizeY, MPI_DOUBLE, rank+1, 1234, MPI_COMM_WORLD, &req[3]);

  }
  // setting elements of most of the points
 // printf("Setting IN change, the rank is %d\n",rank);

  int xStart = 1;
  int xBound = sizeX-1;
  
  for(int x = xStart; x < xBound; x++)
  { 
	for(int y = 1; y < sizeY-1; y++)
	{
	  change[x][y] = fpp(u[x-1][y],u[x][y],u[x+1][y],dx)
	  + fpp(u[x][y-1],u[x][y],u[x][y+1],dx);
	}
  }
  // printf("IN Waiting, rank is %d\n",rank);

  MPI_Waitall(size,req,MPI_STATUSES_IGNORE );
  
  // setting elements of the top and bottom rows
	
  for(int i = 1; i < sizeY-1; i++)
  {
   //printf("After %d iteration, the rank is %d\n",i,rank);
	if( rank != size -1 )
	{
      change[0][i] = fpp(u[0][i-1],u[0][i],u[0][i+1],dx)
	  + fpp(u[1][i],u[0][i],u[1][i],dx);
	  
	  change[sizeX-1][i] = fpp(u[sizeX-2][i],u[sizeX-1][i],rec2[i],dx)
	  + fpp(u[sizeX-1][i-1],u[sizeX-1][i],u[sizeX-1][i+1],dx);
	}	
	if( rank != 0)
	{	  
	  change[sizeX-1][i] = fpp(u[sizeX-1][i-1],u[sizeX-1][i],
	  u[sizeX-1][i+1],dx)
	  + fpp(u[sizeX-2][i],u[sizeX-1][i],u[sizeX-2][i],dx);
	  	  
	  change[0][i] = fpp(rec1[i],u[0][i],u[1][i],dx)
	  + fpp(u[0][i-1],u[0][i],u[0][i+1],dx);
    }
  }
    
 	// setting elements of the right and left cols

  for(int i = 1; i < sizeX-1; i++)
  {
	change[i][0] = fpp(u[i-1][0],u[i][0],u[i+1][0],dx)
	+ fpp(u[i][1],u[i][0],u[i][1],dx);
	
//	printf("rank is %d \n",rank);

	change[i][sizeY-1] = fpp(u[i-1][sizeY-1],u[i][sizeY-1],
	u[i+1][sizeY-1],dx)
	+ fpp(u[i][sizeY-2],u[i][sizeY-1],u[i][sizeY-2],dx);
  }
   
  if(rank == 0)
  {
    change[0][0] = fpp(u[1][0],u[0][0],u[1][0],dx)+
    fpp(u[0][1],u[0][0],u[0][1],dx);
  }
  else
  {
    change[0][0] = fpp(rec1[0],u[0][0],u[1][0],dx) + 
    fpp(u[0][1],u[0][0],u[0][1],dx);
  }

  if(rank == 0)
    change[0][sizeY-1] = fpp(u[1][sizeY-1],u[0][sizeY-1],u[1][sizeY-1],dx) + 
    fpp(u[0][sizeY-2],u[0][sizeY-1],u[0][sizeY-2],dx);
  else
    change[0][sizeY-1] = fpp(rec1[sizeY-1],u[0][sizeY-1],u[1][sizeY-1],dx) + 
    fpp(u[0][sizeY-2],u[0][sizeY-1],u[0][sizeY-2],dx);
  
  if(rank == size-1)
    change[sizeX-1][0] = fpp(u[sizeX-1][1],u[sizeX-1][0],u[sizeX-1][1],dx) + 
    fpp(u[sizeX-2][0],u[sizeX-1][0],u[sizeX-2][0],dx);
  else
    change[sizeX-1][0] = fpp(u[sizeX-1][1],u[sizeX-1][0],u[sizeX-1][1],dx) + 
    fpp(rec2[0],u[sizeX-1][0],u[sizeX-2][0],dx);
    
  if(rank == size-1)	
    change[sizeX-1][sizeY-1] = fpp(u[sizeX-1][sizeY-2],u[sizeX-1][sizeY-1],
    u[sizeX-1][sizeY-2],dx) + fpp(u[sizeX-2][sizeY-1],u[sizeX-1][sizeY-1],
    u[sizeX-2][sizeY-1],dx);
  else
    change[sizeX-1][sizeY-1] = fpp(u[sizeX-1][sizeY-2],u[sizeX-1][sizeY-1],
    u[sizeX-1][sizeY-2],dx) + fpp(rec2[sizeY-1],u[sizeX-1][sizeY-1],
    u[sizeX-2][sizeY-1],dx);

  if(rank != 0)
  free(rec1);
  
  if(rank != size-1)
  free(rec2);

  return change;
}
