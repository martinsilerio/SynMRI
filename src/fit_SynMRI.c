/* File : fit_SynMRI.c: finds the parameters needed for the syntethic MRI*/

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "lbfgsb.h"
#include "matrix.h"
#include "array.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_statistics_double.h>

/*#define MATHLIB_STANDALONE 1*/


#include "3Dimages.h"
#include "myLBFGSB.h"

/* some contants needed */
const double l_rho = 1;
const double u_rho = 400;
const double l_T1 = 300;
const double u_T1 = 12000;
const double l_T2 = 5;
const double u_T2 = 400;
const double inten_MAX = 400;

double *l_W;
double *u_W;

double *l_beta;
double *u_beta;

double *l_psi;
double *u_psi;

double *SIGMA;
double *SIGMA_subset;
double *lambda_x, *lambda_y, *lambda_z;
int flag_eigenvalues = 0;

unsigned int *index_neighbor_x_1, *index_neighbor_x_2;
unsigned int *index_neighbor_y_1, *index_neighbor_y_2;
unsigned int *index_neighbor_z_1, *index_neighbor_z_2;

unsigned int *total_neighbors_x, *total_neighbors_y, *total_neighbors_z;

const int MAX_ITE = 500;
const double TOL = 0.00000001;

/* dimension for the rosenbrock test */
const int N = 25;

double nu_i_j(double *W_i, double TE_j, double TR_j)
{

  double nu;

  nu = inten_MAX*W_i[0]*(exp(TE_j*log(W_i[2])))*(1 - exp(TR_j*log(W_i[1])));

  return nu;

}

void findLSE(int M, double *r_i, double *TE, double *TR, double *params_i,
             double *W_i)
{

  double *NU_i, *factor;
  double **der_NU_i;
  int j;
  double *x0, *g0;
  double f_max;
  double *g_aux;

  void errorAtVoxel(int p, double *x, double *error)
  {

    for(j = 0; j < M; j++)
      NU_i[j] = nu_i_j(x, TE[j], TR[j]);

    *error = 0;

    for(j = 0; j < M; j++)
      *error += (r_i[j] - NU_i[j])*(r_i[j] - NU_i[j]);

  }

  void gradErrorAtVoxel(int p, double *x, double *gradError)
  {

    for(j = 0; j < M; j++)
      NU_i[j] = nu_i_j(x, TE[j], TR[j]);

    for(j = 0; j < M; j++)
      factor[j] = 2*(NU_i[j] - r_i[j]);

    for(j = 0; j < M; j++)
    {
      der_NU_i[j][0] = NU_i[j]/x[0];
      der_NU_i[j][1] = (1/x[1])*TR[j]*(NU_i[j] - inten_MAX*x[0]*exp(log(x[2])*TE[j]));
      der_NU_i[j][2] = NU_i[j]*TE[j]/x[2];
    }

    gradError[0] = 0;
    gradError[1] = 0;
    gradError[2] = 0;

    for(j = 0; j < M; j++)
    {
      gradError[0] += der_NU_i[j][0]*factor[j];
      gradError[1] += der_NU_i[j][1]*factor[j];
      gradError[2] += der_NU_i[j][2]*factor[j];
    }

    /*printMessage("gradient ");
     printVector(gradError, 3);
     computeNumGradient(3, errorAtVoxel, x, g_aux, 0.00001);
     printMessage("num gradient ");
     printVector(g_aux, 3);
     printMessage("x ");
     printVector(x, 3);
     printMessage("lower ");
     printVector(l_W, 3);
     printMessage("upper ");
     printVector(u_W, 3);
     makePause();*/

  }

  MAKE_VECTOR(NU_i, M);

  MAKE_VECTOR(factor, M);
  MAKE_MATRIX(der_NU_i, M, 3);
  MAKE_VECTOR(x0, 3);
  MAKE_VECTOR(g_aux, 3);

  /* let set the initial value */
  x0[0] = 1.2*maxOfVector(r_i, M);
  x0[1] = 0.8*minOfVector(TR, M);
  x0[2] = 1.2*maxOfVector(TE, M);

  printVector(x0, 3);

  /* let's transform them */
  x0[0] = x0[0]/inten_MAX;
  x0[1] = exp(-1/x0[1]);
  x0[2] = exp(-1/x0[2]);

  /*printMessage("x0 ");
   printVector(x0, 3);
   printMessage("r_i ");
   printVector(r_i, 3);
   printMessage("TR ");
   printVector(TR, 3);
   printMessage(" ");
   printVector(TE, 3);
   makePause();*/

  optimLBFGSB(3, errorAtVoxel, gradErrorAtVoxel, x0, l_W, u_W, MAX_ITE, TOL,
              W_i, &f_max);

  /* printing some information */
  printf("max value %lf at ", f_max);
  printVector(W_i, 3);
  printf("\n");

  /*makePause();*/

  if(W_i[0] == 0)
  {
    scanf("%d", &j);
  }

  /* let's save the parameters on the original scale */
  params_i[0] = inten_MAX*W_i[0];
  params_i[1] = -1/log(W_i[1]);
  params_i[2] = -1/log(W_i[2]);

  FREE_VECTOR(x0);
  FREE_MATRIX(der_NU_i);
  FREE_VECTOR(factor);
  FREE_VECTOR(NU_i);
  FREE_VECTOR(g_aux);

}

void findLeastSquaresEstimates(int M, img3d *r, double *TE, double *TR,
                               img3d *params, img3d *W)
{
  /* the arrays of images params and W need to have memory already reserved */

  int n = r[0].n;
  int i, j;
  double *r_i, *params_i, *W_i;

  /* let's save n1, n2, n3 and n on the images params and w */
  for(i = 0; i < 3; i++)
  {
    params[i].n = r[0].n;
    params[i].n1 = r[0].n1;
    params[i].n2 = r[0].n2;
    params[i].n3 = r[0].n3;
    W[i].n = r[0].n;
    W[i].n1 = r[0].n1;
    W[i].n2 = r[0].n2;
    W[i].n3 = r[0].n3;
  }

  MAKE_VECTOR(r_i, M);
  MAKE_VECTOR(params_i, 3);
  MAKE_VECTOR(W_i, 3);

  /* let's find the LSE for each voxel */
  for(i = 0; i < n; i++)
  {

    for(j = 0; j < M; j++)
      r_i[j] = r[j].intensity[i];

    findLSE(M, r_i, TE, TR, params_i, W_i);

    params[0].intensity[i] = params_i[0];
    params[1].intensity[i] = params_i[1];
    params[2].intensity[i] = params_i[2];
    W[0].intensity[i] = W_i[0];
    W[1].intensity[i] = W_i[1];
    W[2].intensity[i] = W_i[2];

  }



  FREE_VECTOR(W_i);
  FREE_VECTOR(params_i);
  FREE_VECTOR(r_i);



}

void computePredictedImagesWithW(img3d *params, int M, double *TE, double *TR,
                                 img3d *predImage)
  /* predImage should have memory already reserved */
{

  int i, j;
  int n;
  double *W_i;

  n = params[0].n;
  MAKE_VECTOR(W_i, 3);

  for(j = 0; j < M; j++)
  {

    predImage[j].n = params[0].n;
    predImage[j].n1 = params[0].n1;
    predImage[j].n2 = params[0].n2;
    predImage[j].n3 = params[0].n3;

    for(i = 0; i < n; i++)
    {
      W_i[0] = params[0].intensity[i];
      W_i[1] = params[1].intensity[i];
      W_i[2] = params[2].intensity[i];
      predImage[j].intensity[i] = nu_i_j(W_i, TE[j], TR[j]);
    }
  }

  FREE_VECTOR(W_i);

}

double test = 9;

/*void matprod1 ( double *X , int * dimX , double *Y , int * dimY , double * ans ){
  double sum ;
  int ii , jj , kk ;
  int nrX = dimX [0] , ncX = dimX [1] , nrY = dimY [0] , ncY = dimY [1];

  for ( ii =0; ii < nrX ; ii ++){
    for ( jj =0; jj < ncY ; jj ++){
      sum = 0;
      for ( kk =0; kk < ncX ; kk ++){
        sum = sum + X[ ii + nrX * kk ]* Y[ kk + nrY * jj ];
      }
      ans [ ii + nrX * jj ] = sum ;
    }
  }
}*/

SEXP fit_LSE_(SEXP M_, SEXP n_, SEXP n1_, SEXP n2_, SEXP n3_, SEXP R_, SEXP TE_,
              SEXP TR_)
{
  int M = asInteger(M_);
  int n = asInteger(n_);
  int n1 = asInteger(n1_);
  int n2 = asInteger(n2_);
  int n3 = asInteger(n3_);
  int nprot = 0;
  PROTECT (R_ = AS_NUMERIC(R_)); nprot ++;
  PROTECT (TE_ = AS_NUMERIC(TE_)); nprot ++;
  PROTECT (TR_ = AS_NUMERIC(TR_)); nprot ++;

  int i, j, c;
  int i_x, i_y, i_z, index;
  double *TE, *TR;
  TR = REAL(TR_);
  TE = REAL(TE_);
  double *r;
  r = REAL(R_);

  img3d *image;
  img3d *params, *W;

  double error;

  /*set_seed(111,222);*/

  MAKE_VECTOR(image, M);

  c = 0;
  for(i = 0; i<M; i++)
  {

    MAKE_VECTOR(image[i].intensity, n);

    image[i].n = n;
    image[i].n1 = n1;
    image[i].n2 = n2;
    image[i].n3 = n3;

    /* when the intensity is 0, we subtitute it for a random value
     between 0 and 0.1  */
    for(j = 0; j<n; j++)
    {
      image[i].intensity[j] = r[c];
      /*if(image[i].intensity[j] == 0)
        image[i].intensity[j] = runif(0, 0.1);*/
      c++;
    }

  }

  saveSlices(&(image[0]), "image1.pgm");
  saveSlices(&(image[1]), "image10.pgm");
  saveSlices(&(image[2]), "image12.pgm");

  /* Now we can construct the vectors of indexes of the neighbors and total
   neighbors */

  MAKE_VECTOR(index_neighbor_x_1, n);
  MAKE_VECTOR(index_neighbor_x_2, n);
  MAKE_VECTOR(index_neighbor_y_1, n);
  MAKE_VECTOR(index_neighbor_y_2, n);
  MAKE_VECTOR(index_neighbor_z_1, n);
  MAKE_VECTOR(index_neighbor_z_2, n);
  MAKE_VECTOR(total_neighbors_x, n);
  MAKE_VECTOR(total_neighbors_y, n);
  MAKE_VECTOR(total_neighbors_z, n);

  for(i_x = 0; i_x < n1; i_x++)
    for(i_y = 0; i_y < n2; i_y++)
      for(i_z = 0; i_z < n3; i_z++)
      {

        index = getArrayIndex(&(image[0]), i_x, i_y, i_z);

        total_neighbors_x[index] = 0;
        total_neighbors_y[index] = 0;
        total_neighbors_z[index] = 0;

        /* neighbors in x */
        if(i_x != 0)
        {
          total_neighbors_x[index]++;
          index_neighbor_x_1[index] = getArrayIndex(&(image[0]), i_x - 1, i_y, i_z);
        }
        else
          index_neighbor_x_1[index] = n;

        if(i_x != n1 - 1)
        {
          total_neighbors_x[index]++;
          index_neighbor_x_2[index] = getArrayIndex(&(image[0]), i_x + 1, i_y, i_z);
        }
        else
          index_neighbor_x_2[index] = n;

        /* neighbors in y */
        if(i_y != 0)
        {
          total_neighbors_y[index]++;
          index_neighbor_y_1[index] = getArrayIndex(&(image[0]), i_x, i_y - 1, i_z);
        }
        else
          index_neighbor_y_1[index] = n;

        if(i_y != n2 - 1)
        {
          total_neighbors_y[index]++;
          index_neighbor_y_2[index] = getArrayIndex(&(image[0]), i_x, i_y + 1, i_z);
        }
        else
          index_neighbor_y_2[index] = n;

        /* neighbors in z */
        if(i_z != 0)
        {
          total_neighbors_z[index]++;
          index_neighbor_z_1[index] = getArrayIndex(&(image[0]), i_x, i_y, i_z - 1);
        }
        else
          index_neighbor_z_1[index] = n;

        if(i_z != n3 - 1)
        {
          total_neighbors_z[index]++;
          index_neighbor_z_2[index] = getArrayIndex(&(image[0]), i_x, i_y, i_z + 1);
        }
        else
          index_neighbor_z_2[index] = n;

      }

  /* Let's reserve memory and compute the bounds for the transformations of
  rho, T1 and T2 */
  MAKE_VECTOR(l_W, 3);
  MAKE_VECTOR(u_W, 3);

  l_W[0] = l_rho;
  l_W[1] = exp(-1/l_T1);
  l_W[2] = exp(-1/l_T2);
  u_W[0] = u_rho;
  u_W[1] = exp(-1/u_T1);
  u_W[2] = exp(-1/u_T2);

  /* let's find the estimates */

  MAKE_VECTOR(params, 3);
  MAKE_VECTOR(W, 3);

  for(i = 0; i < 3; i++)
  {
    MAKE_VECTOR(params[i].intensity, n);
    MAKE_VECTOR(W[i].intensity, n);
  }

  printMessage("Before Fitting the model");

   findLeastSquaresEstimates(M, image, TE, TR, params, W);

  printMessage("After Fitting the model");
  /*makePause();*/

  SEXP params_;
  PROTECT (params_ = allocVector(REALSXP, 3*n)); nprot ++;
  double *params_ptr ; params_ptr = REAL(params_);

  SEXP W_;
  PROTECT (W_ = allocVector(REALSXP, 3*n)); nprot ++;
  double *W_ptr ;   W_ptr = REAL(W_);

  c = 0;
  for(j = 0; j < 3; j++)
    for(i = 0; i < n; i++)
    {
      /*REAL(params_)[c] = params[j].intensity[i];
      REAL(W_)[c] = W[j].intensity[i];*/
      params_ptr[c] = params[j].intensity[i];
      W_ptr[c] = W[j].intensity[i];
      c++;
    }

  SEXP aux;
  PROTECT (aux = allocVector(REALSXP, 3)); nprot ++;
  double *ansptr ; ansptr = REAL(aux);

  ansptr[0] = test;
  ansptr[1] = 2*test;
  ansptr[2] = 3*test;




  SEXP OS = PROTECT(allocVector(VECSXP, 2));
  nprot++;
  SET_VECTOR_ELT(OS, 0, W_);
  SET_VECTOR_ELT(OS, 1, params_);


  for(i = 0; i<M; i++)
    FREE_VECTOR(image[i].intensity);
  FREE_VECTOR(image);

  for(i = 0; i<M; i++)
    FREE_VECTOR(params[i].intensity);
  FREE_VECTOR(params);

  for(i = 0; i<M; i++)
    FREE_VECTOR(W[i].intensity);
  FREE_VECTOR(W);

  FREE_VECTOR(SIGMA);

  FREE_VECTOR(u_W);
  FREE_VECTOR(l_W);
  FREE_VECTOR(u_beta);
  FREE_VECTOR(l_beta);
  FREE_VECTOR(u_psi);
  FREE_VECTOR(l_psi);

  FREE_VECTOR(lambda_x);
  FREE_VECTOR(lambda_y);
  FREE_VECTOR(lambda_z);

  FREE_VECTOR(index_neighbor_x_1);
  FREE_VECTOR(index_neighbor_x_2);
  FREE_VECTOR(index_neighbor_y_1);
  FREE_VECTOR(index_neighbor_y_2);
  FREE_VECTOR(index_neighbor_z_1);
  FREE_VECTOR(index_neighbor_z_2);
  FREE_VECTOR(total_neighbors_x);
  FREE_VECTOR(total_neighbors_y);
  FREE_VECTOR(total_neighbors_z);

  UNPROTECT(nprot);
  return OS;
}
