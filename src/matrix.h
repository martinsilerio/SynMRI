

#ifndef matrix_h
#define matrix_h

/*This function allocates memory for a double matrix of size nr x nc */
double** doubleMatrix(int, int);

/*This function frees the memory allocated for a matrix of doubles with nr rows */
void freeDoubleMatrix(double**, int);

/*This function assigns zero to all the elements of a matrix */
void zeroMatrix(double**, int, int);

/*This function reserves memory for a vector of doubles with size n */
double* doubleVector(int);

/*This function frees the memory allocated for a vector of doubles */
void freeDoubleVector(double*);

/*This function assigns zero to all the elements of a vector */
void zeroVector(double*, int);

/*This function calculates the inner product of two vectors of doubles */
double innerProduct(double*, double*, int);

/*This function calculates the outer product of two vectors of doubles */
void outerProduct(double*, int, double*, int, double**);

/*This function calculates the product of a matrix and a vector of doubles */
void productMatrixVector(double **, int, int, double*, double*);

/*This function copies the elements of a vector of doubles to another vector of doubles*/
void assignVector(double*, double*, int);

/*This function copies the elements of a matrix of doubles to another matrix of doubles*/
void assignMatrix(double**, double**, int, int);

/*This function computes the norm of a vector of doubles*/
double norm(double*, int);

/*This function normalizes a vector of doubles*/
void normalize(double*, int);

/*This function computes the scalar product of a scalar and a vector of doubles*/
void scalarProduct(double*, int, double);

/*This function adds a vector of doubles to another vector*/
void sumTo(double*, double*, int);

/*This function adds a matrix of doubles to another matrix*/
void sumToMatrix(double**, double**, int, int);

void printVector(double*, int);

void printMatrix(double**, int, int);

double maxOfVector(double*, int);

double minOfVector(double*, int);

#endif /* matrix_h*/
