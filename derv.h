#include <stdbool.h>


#ifndef DERV_H
#define DERV_H

double fppp(double ip, double im, double i, double ipp, double imm, double h);
double fpp (double ip, double im, double i, double h);
double fp ( double ip, double im, double h );

void freePointers(double** yourp, int s);

double printTo(FILE* file, double** f,int sizeX, int sizeY,int rank);

void printOut(double** f, int sizeX, int sizeY);

double** alloc(int sizeO, int sizeI);
double absVal(double a);
void deall (double*** col, int num, int out);

double** change( double** u, double dx, int sizeX, int sizeY );

double** add(double** f, double** sA, double** sB, double** sC, 
double** sD, double** sE, double constA, double constB, double constC,
double constD, double constE, int sizeX, int sizeY);

double** kutta(double** f,double dx, double tol,
int sizeX,int sizeY,double* step );

double norm(double** first, double ** second,int sizeX,
int sizeY,double powE,double tCoeff);

bool isAlmostUni (double** f, int sizeX, int sizeY, double tol);



#endif
