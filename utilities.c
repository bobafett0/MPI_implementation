#include <stdio.h>
#include <stdlib.h>
#include "derv.h"
#include <math.h>
#include <stdbool.h>
#include <mpi.h>

// "Prints" the arrays onto a text file

double printTo(FILE* file, double** f,int sizeX, int sizeY,int rank)
{
  fprintf(file,"The average is %0.3f \n", norm(f,f,sizeX,sizeY,1.0,0.0));
  for(int r = 0; r < sizeX; r++)
  {
    for(int u = 0; u < sizeY; u++)
	{
	  printf("%0.9f ",f[r][u]);
	  fprintf(file,"%0.9f ",f[r][u]);
	}
	printf("\n");
    fprintf(file,"\n");
  }
  fprintf(file,"\n\n");
 // printf("in printTo\n");
}

// Prints out a matrix. Made for debugging purposes.

void printOut(double** f, int sizeX, int sizeY)
{
  for(int i = 0; i < sizeX; i++)
  {
    for(int u = 0; u < sizeY; u++)
    {
	  printf(" %0.4f ",f[i][u]);
	}
	printf("\n");
  }
}

// Frees the memory refered to by a double pointer, ie, 
// the matrix

void freePointers(double** yourp, int s)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  for(int i = 0; i < s; i++)
  {
    free(yourp[i]);
  }
  free(yourp);
}

// Allocates the memory for a matrix, or a double pointer for doubles

double** alloc(int sizeO, int sizeI)
{
  double** f = calloc(sizeO,sizeof(double*));
  
  // Initializing all pointers of f
  for(int i = 0; i < sizeO; i++)
  {
	f[i] = calloc(sizeI,sizeof(double));
  }
  
  return f;
}

// Calculates the user-specified for a 2D vector

double norm(double** first, double ** second,int sizeX,int sizeY,double powE,double tCoeff)
{
  double sum = 0.0;
  for(int i = 0; i < sizeX; i++)
	for(int u = 0; u < sizeY; u++)
	{
	  sum += pow(absVal(first[i][u] - tCoeff*second[i][u]),powE);
	}
  return sum/((double)(sizeX*sizeY));
}

// Calculates the absolute value of a number and returns it

double absVal(double a)
{
  if(a < 0.0)
  return (-1.0)*a;
	
  return a;
}

// Uses deallocates memory for a collection of pointers

void deall (double*** col, int num, int out)
{
  for(int i = 0; i < num; i++)
  {
	freePointers(col[i],out);
  }
  free(col);
}

bool isAlmostUni (double** f, int sizeX, int sizeY, double tol)
{
	for(int i = 0; i < sizeX; i++)
	  for(int u = 0; u < sizeY; u++)
	    {
			if( absVal(f[i][u] - f[0][0]) > tol && (i != 0 || u != 0))
			  return false;
		}
	return true;
}


