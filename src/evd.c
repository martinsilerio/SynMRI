/* @file evd.c */
/* Part of the answer to Problem 1, Homework 4, STAT 580, Spring 2013 (ISU) */
/* @author Martin Silerio-Vazquez (silerio@iastate.edu) */
/* @date Sat Apr 20 20:31:00 CST 2013 */
/* */
/* This program contains the adaptation of the function "dgeev_" from LAPACK */
/* to be used in C. It was made following the example seen in class in the */
/* file svd.c, and also checking the LAPACK's documentation on */
/* http://www.netlib.org/lapack/explore-html/d9/d8e/group__double_g_eeigen.html */
/* */
/* We compile this code to construct object code with:

gcc -c evd.c -ansi -Wall -pedantic

*/

#include <stdio.h>
#include <stdlib.h>
#include "array.h"

/* Here we declare the function from LAPACK that we will use */
/* It is defined somewhere else */
void dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda, double *wr,
            double *wi, double *vl, int *ldvl, double *vr, int *ldvr, double *work,
            int *lwork, int *info);

/* We receive only a square matrix, its dimension and arrays and matrices */
/* to save the eigenvalues and the eigenvectors matrices */
/* One of those matrices will contain the eigenvectors (right eigenvectors) */
/* The another one contains the left eigenvectors. We only use them */
/* to explore what is what this function does */
/* The memory of these structures is assumed to be allocated already */
int evdd(double **a, int n, double *d, double **u, double **v)
   {
   /* We assign value to some parameters, to obtain both transformation */
   /* matrices (actually one should be the inverse of the another) */
   /* we assign the value of "lwork" only following the documentation */
   /* I do not understand completely what it does */
   double *AT, *VL, *VR;
   int lwork = 4*n;
   char jobvl = 'V';
   char jobvr = 'V';
   double *di,*work;
   int i, j, k, info;

   /* Here we reserve memory for some vectors and call the corresponding function */
   MAKE_VECTOR(AT, n*n);
   for (j=0, k=0; j<n; j++)
      for (i=0; i<n; i++)
         AT[k++] = a[i][j];

   MAKE_VECTOR(VL, n*n);
   MAKE_VECTOR(VR, n*n);

   MAKE_VECTOR(di, n);

   MAKE_VECTOR(work, lwork);

   dgeev_(&jobvl, &jobvr, &n, AT, &n, d, di, VL, &n, VR, &n, work, &lwork, &info);

   /* We save the tranformation matrices taking on count that lapack uses */
   /* a different order in the indexes of the matrices */
   for (j=0, k=0; j<n; j++)
      for (i=0; i<n; i++)
         u[i][j] = VL[k++];

   for (j=0, k=0; j<n; j++)
      for (i=0; i<n; i++)
         v[i][j] = VR[k++];

   /*Finally, we free the reserved memory and return the info variable */
   FREE_VECTOR(AT);
   FREE_VECTOR(VL);
   FREE_VECTOR(VR);
   FREE_VECTOR(work);

   return info;
   }
