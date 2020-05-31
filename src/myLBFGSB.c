#include "lbfgsb.h"
#include "array.h"
#include "matrix.h"

void optimLBFGSB(int p,
                 void (*fx)(int, double*, double*),
                 void (*gx)(int, double*, double*),
                 double *x0,
                 double *l_bound,
                 double *u_bound,
                 int max_ite,
                 double tol,
                 double *x_max,
                 double *fx_max){
/* x_max should be an array of size p with space already reserved */

/* Here we declare other variable that we will use */
int DONE =  0;
int i;

/* First we need to declare/initialize a LOT of variables that the l-bfgs-b code
needs */
/* To avoid mistakes and/or confusion I will use exactly the same names
used in the documentation */
integer n = p;
integer m = 12;
double *x;
double *l;
double *u;
integer *nbd;
double f;
double *g;
double factr = 1e7;
double pgtol = tol;
double *wa;
integer *iwa;
integer taskValue;
integer *task=&taskValue;
/*char task[60];*/
integer iprint;
integer csaveValue;
integer *csave=&csaveValue;
/*char csave[60];*/
logical lsave[4];
integer isave[44];
double dsave[29];

/* Here we reserve space for all the variables that the l-bfgs-b code needs */
MAKE_VECTOR(x, n);
MAKE_VECTOR(l, n);
MAKE_VECTOR(u, n);
MAKE_VECTOR(nbd, n);
MAKE_VECTOR(g, n);
MAKE_VECTOR(wa, (2*m + 5)*n + 11*m*m + 8*m);
MAKE_VECTOR(iwa, 3*n);

/* let's initialice the necessary values */

assignVector(x, x0, n);
assignVector(l, l_bound, n);
assignVector(u, u_bound, n);
for(i = 0; i < n; i++)
   nbd[i] = 2;
iprint = -1;
factr = 1e7;
pgtol = 1e-5;
*task = (integer)START;

/* we are ready to start iterating the algortihm */

while(!DONE)
   {

   /* we call the l-bfgs-b function */
   setulb(&n, &m, x, l, u, nbd, &f, g, &factr, &pgtol, wa, iwa, task, &iprint,
          csave, lsave, isave, dsave);

   /* depending on the output we will procede on a different way */
   if(IS_FG(*task))
      {
      /* this means the function is requesting to compute both the value of
      the function and the gradient at x */

      fx(n, x, &f);
      gx(n, x, g);

      }
   else if(*task!=NEW_X)
      {
      /* if the algorithm did not return a new x and neither is as asking for
      the value of the function and gradient, the run is done */
      DONE = 1;
      }

   /* if the program reaches this line without being done, it means a new x was
   provided in the previous call */

   }

/* Here we just return both the value x that maximizes the function and the
value */
assignVector(x_max, x, n);
*fx_max = f;

/* Here we free the space reserved for all the variables used by the
l-bfgd-b */
FREE_VECTOR(x);
FREE_VECTOR(l);
FREE_VECTOR(u);
FREE_VECTOR(nbd);
FREE_VECTOR(g);
FREE_VECTOR(wa);
FREE_VECTOR(iwa);

/* Here we free the rest of the space reserved */

}
