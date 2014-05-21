
#include <stdio.h>
#include <stdlib.h>
#include "derv.h"
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <assert.h>

double** kutta(double** f,double dx, double tol,
int sizeX,int sizeY,double* step )
{
  //printf("printing original : \n\n");
  //printOut(f,sizeX,sizeY);
  //printf( " The average value is %0.4f ", norm(f,f,sizeX,sizeY,1.0,0.0) );

  double dt = *step;
  double*** s = calloc(6,sizeof(double**));
  double*** fn = calloc(5,sizeof(double**));
  
  s[0] = change(f ,dx,sizeX,sizeY);
  //printf ( "printing s0\n" );
 // printOut(s[0],sizeX,sizeY);  

  fn[0] = add(f,s[0],s[0],s[0],s[0],s[0],0.25*dt,
  0.0,0.0,0.0,0.0,sizeX,sizeY);
 // printf ( "printing f0\n" );
 // printOut(fn[0],sizeX,sizeY);
  
  s[1] = change(fn[0],dx,sizeX,sizeY);
 // printf ( "printing s1\n" );
 // printOut(s[1],sizeX,sizeY);  

  fn[1] = add(f,s[0],s[1],s[1],s[1],s[1],(3.0/32.0)*dt,
  (9.0/32.0)*dt,0.0,0.0,0.0,sizeX,sizeY);
 // printf ( "printing f1\n" );
 // printOut(fn[1],sizeX,sizeY);
  
  s[2] = change(fn[1],dx,sizeX,sizeY);
 // printf ( "printing s2\n" );
 // printOut(s[2],sizeX,sizeY);  
	
  fn[2] = add(f,s[0],s[1],s[2],s[2],s[2],(1932.0/2197.0)*dt,
  (-7200.0/2197.0)*dt,(7296.0/2197.0)*dt,0.0,0.0,sizeX,sizeY);
 // printf ( "printing f2\n" );
 // printOut(fn[2],sizeX,sizeY);
  
  s[3] = change(fn[2],dx,sizeX,sizeY);
 // printf ( "printing s3\n" );
 // printOut(s[3],sizeX,sizeY);  
  
  fn[3] = add(f,s[0],s[1],s[2],s[3],s[3],(439.0/216.0)*dt, -8.0*dt,
  (3680.0/513.0)*dt, (-845.0/4104.0)*dt,0.0,sizeX,sizeY);
 // printf ( "printing f3\n" );
 // printOut(fn[3],sizeX,sizeY);
  
  s[4] = change(fn[3],dx,sizeX,sizeY);
 // printf ( "printing s4\n" );
 // printOut(s[4],sizeX,sizeY);    
	
  fn[4] = add(f,s[0],s[1],s[2],s[3],s[4],(-8.0/27.0)*dt,2.0*dt,
  (-3544.0/2565.0)*dt,(1859.0/4104.0)*dt,(-11.0/40.0)*dt,sizeX,sizeY);
 // printf ( "printing f4\n" );
 // printOut(fn[4],sizeX,sizeY);
  
  s[5] = change(fn[4],dx,sizeX,sizeY);
 // printf ( "printing s5\n" );
 // printOut(s[5],sizeX,sizeY);
	
  double** est4r = add(f,s[0],s[2],s[3],s[4],s[0],dt*(25.0/216.0),dt*(1408.0/2565.0),
  dt*(2197.0/4104.0),-dt*(1.0/5.0),0.0,sizeX,sizeY);
 // printf ( "printing est4r\n" );
 // printOut(est4r,sizeX,sizeY);
  
	
  double** est5r = add(f,s[0],s[2],s[3],s[4],s[5],dt*(16.0/1351.0),dt*(6656.0/12825.0),
  dt*(28561.0/56430.0),dt*(-9.0/50.0),dt*(2.0/55.0),sizeX,sizeY);
 // printf ( "printing est5r\n" );
 // printOut(est5r,sizeX,sizeY);
  	
  deall(fn,5,sizeX);
  deall(s,6,sizeX);
	
  double error = norm(est4r,est5r,sizeX,sizeY,1.0,1.0);
  freePointers(est4r,sizeX);
  
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Request req[2];
  
  double* recv = calloc(1,sizeof(double));
  double* send = calloc(1,sizeof(double));
    
  if(error/tol > 1.2) 
    send[0] = 542545.6767676767;
  else
    send[0] = 0.0;
    
    send[1] = *step;
  if ( rank == 0 )
  {
    MPI_Isend(send, 1, MPI_DOUBLE, rank+1, 1235, MPI_COMM_WORLD, &req[0]);
    MPI_Irecv(recv, 1, MPI_DOUBLE, rank+1, 1235, MPI_COMM_WORLD, &req[1]);
  }
  else if ( rank == size-1 )
  {
	MPI_Irecv(recv, 1, MPI_DOUBLE, rank-1, 1235, MPI_COMM_WORLD, &req[1]);
    MPI_Isend(send, 1, MPI_DOUBLE, rank-1, 1235, MPI_COMM_WORLD, &req[0]);
  }
  
  MPI_Waitall(size,req,MPI_STATUSES_IGNORE );
  
  if(error/tol > 1.2 || recv[0] != 0)
  {
	free(recv); free(send);
	
	recv = calloc(1,sizeof(double));
    send = calloc(1,sizeof(double));
    
    send[0] = *step = dt = dt*pow(tol/error,1.0/5.0);
    
    if ( rank == 0 )
    {
      MPI_Isend(send, 1, MPI_DOUBLE, rank+1, 1235, MPI_COMM_WORLD, &req[0]);
      MPI_Irecv(recv, 1, MPI_DOUBLE, rank+1, 1235, MPI_COMM_WORLD, &req[1]);
    }
    else if ( rank == size-1 )
    {
	  MPI_Irecv(recv, 1, MPI_DOUBLE, rank-1, 1235, MPI_COMM_WORLD, &req[1]);
      MPI_Isend(send, 1, MPI_DOUBLE, rank-1, 1235, MPI_COMM_WORLD, &req[0]);
    }
  
    MPI_Waitall(size,req,MPI_STATUSES_IGNORE );  
    
    if(recv[0] < *step)
      *step = recv[0];
      
    
      
    //printf(" correcting, dt is %0.35f RanKK is %d\n",*step,rank );
    freePointers(est5r,sizeX);	
	free(recv); free(send);
	return kutta(f,dx,tol,sizeX,sizeY,step);
  }
  free(recv); free(send);  
  return est5r;	
}



double** add(double** f, double** sA, double** sB, double** sC, 
double** sD, double** sE, double constA, double constB, double constC,
double constD, double constE, int sizeX, int sizeY)
{
  double** ret = alloc(sizeX,sizeY);
	
  for(int i = 0; i < sizeX; i++)
	for(int u = 0; u < sizeY; u++)
	{
	  ret[i][u] = f[i][u] + constA*sA[i][u] + constB*sB[i][u] + 
	  constC*sC[i][u] + constD*sD[i][u] + constE*sE[i][u];			
	}
  return ret;	
}
