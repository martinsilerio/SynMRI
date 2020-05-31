/*We include stdlib, that is necessary for the memory allocation*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

/*This function allocates memory for a double matrix of size nr x nc */
double** doubleMatrix(int nr, int nc)
	{
	double **A;
	int i;
	A=(double**)malloc(nr*sizeof(double*));
	for(i=0;i<nr;i++)
		{
		A[i]=(double*)malloc(nc*sizeof(double));
		}
	return A;
	}

/*This function frees the memory allocated for a matrix of doubles with nr rows */
void freeDoubleMatrix(double **A,int nr)
	{
	int i;
	for(i=0;i<nr;i++)
		{
		free(A[i]);
		}
        free(A);
	}

/*This function assigns zero to all the elements of a matrix */
void zeroMatrix(double **A, int nr, int nc)
   {
   int i,j;
   for(i=0;i<nr;i++)
      {
      for(j=0;j<nc;j++)
         {
         A[i][j]=0;
         }
      }
   }

/*This function reserves memory for a vector of doubles with size n */
double* doubleVector(int n)
   {
   double *x;
   x=(double*)malloc(n*sizeof(double));
   return(x);
   }

/*This function frees the memory allocated for a vector of doubles */
void freeDoubleVector(double *x)
   {
   free(x);
   }

/*This function assigns zero to all the elements of a vector */
void zeroVector(double *x, int n)
   {
   int i;
   for(i=0;i<n;i++)
      {
      x[i]=0;
      }
   }

/*This function calculates the inner product of two vectors of doubles */
double innerProduct(double *x, double *y, int n)
   {
   double sum;
   int i;
   sum=0;
   for(i=0;i<n;i++)
      {
      sum=sum+x[i]*y[i];
      }
   return(sum);
   }

/*This function calculates the outer product of two vectors of doubles */
void outerProduct(double *x, int nr, double *y, int nc, double **A)
   {
   int i,j;
   for(i=0;i<nr;i++)
      {
      for(j=0;j<nc;j++)
         {
         A[i][j]=x[i]*y[j];
         }
      }
   }

/*This function calculates the product of a matrix and a vector of doubles */
void productMatrixVector(double **A, int nr, int nc, double *x, double *b)
   {
   int i;
   for(i=0;i<nr;i++)
      {
      b[i]=innerProduct(A[i],x,nc);
      }
   }

/*This function copies the elements of a vector of doubles to another vector of doubles*/
void assignVector(double *x, double *y,int n)
   {
   int i;
   for(i=0;i<n;i++)
      {
      x[i]=y[i];
      }
   }

/*This function copies the elements of a matrix of doubles to another matrix of doubles*/
void assignMatrix(double **A, double **B, int nr, int nc)
   {
   int i,j;
   for(i=0;i<nr;i++)
      {
      for(j=0;j<nc;j++)
         {
         A[i][j]=B[i][j];
         }
      }
   }

/*This function computes the norm of a vector of doubles*/
double norm(double *x, int n)
   {
   return(sqrt(innerProduct(x,x,n)));
   }

/*This function normalizes a vector of doubles*/
void normalize(double *x, int n)
   {
   double norm_x;
   int i;
   norm_x=norm(x,n);
   if(norm_x!=0)
      {
      for(i=0;i<n;i++)
         {
         x[i]=x[i]/norm_x;
         }
      }
   }

/*This function computes the scalar product of a scalar and a vector of doubles*/
void scalarProduct(double *x, int n, double lambda)
   {
   int i;
   for(i=0;i<n;i++)
      {
      x[i]=lambda*x[i];
      }
   }

/*This function adds a vector of doubles to another vector*/
void sumTo(double *a, double *b, int n)
   {
   int i;
   for(i=0;i<n;i++)
      {
      a[i]=a[i]+b[i];
      }
   }

/*This function adds a matrix of doubles to another matrix*/
void sumToMatrix(double **A, double **B, int nr, int nc)
   {
   int i,j;
   for(i=0;i<nr;i++)
      {
      for(j=0;j<nc;j++)
         {
         A[i][j]+=B[i][j];
         }
      }
   }

void printVector(double *x, int n)
{
int i;

for(i=0; i<n; i++)
   printf("%.4lf ", x[i]);

printf("\n");

}

void printMatrix(double **x, int nr, int nc)
{
int i, j;

for(i=0; i<nr; i++)
   {
   for(j = 0; j < nc; j++)
      printf("%.10lf ", x[i][j]);
   printf("\n");
   }

printf("\n");

}

void printMessage(char *message)
   {
   printf("---> %s\n", message);
   }

void makePause()
   {
   char c;
   scanf("%c", &c);
   }

double maxOfVector(double *x, int n)
{
int i;
double max;

max = x[0];

for(i=1; i<n; i++)
   if(x[i] > max)
      max = x[i];

return(max);

}

double minOfVector(double *x, int n)
{
int i;
double min;

min = x[0];

for(i=1; i<n; i++)
   if(x[i] < min)
      min = x[i];

return(min);

}
