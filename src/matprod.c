/* File : matprod1 .c: Calculate the product of matrices X and Y */


#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>



void matprod1 ( double *X , int * dimX , double *Y , int * dimY , double * ans ){
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
}


SEXP matprod2_(SEXP X , SEXP Y) {
  int nprot = 0;

  PROTECT (X = AS_NUMERIC(X)); nprot ++;
  PROTECT (Y = AS_NUMERIC(Y)); nprot ++;

  double *xptr ; xptr = REAL(X);
  double *yptr ; yptr = REAL(Y);
  int *dimX ; dimX = INTEGER(GET_DIM(X));
  int *dimY ; dimY = INTEGER(GET_DIM(Y));
  SEXP ans ;
  PROTECT (ans = allocMatrix(REALSXP, dimX [0], dimY [1])); nprot ++;
  double *ansptr ; ansptr = REAL(ans);

  matprod1( xptr, dimX, yptr, dimY , ansptr );

  UNPROTECT( nprot );
  return ( ans );
}

