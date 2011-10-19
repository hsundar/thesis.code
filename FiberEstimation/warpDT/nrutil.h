#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <iostream.h>
//#include <strstream.h>
#include <strstream>
#include <fstream.h>

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


//=========== By XDR--Top  ==========

void mulMatrixABt(float **A, float **B, float **C, int row, int col);
  //matrix A (row x col), B (row x col), C (row x row)
  // C = A * B'

void mulMatrixAB(float **A, float **B, float **C, int dim);
  // C = A * B; 
  // A & B & C are (dim x dim)
void mulMatrixAB(float **A, float **B, float **C, int row, int col);
  // C=A*B, 
  // A (row x col)  B (col x row) C are (row x row)

float pythag(float a, float b);
double pythagDouble(double a, double b);

int svdcmp(float **a, int m, int n, float w[], float **v);
  // SVD decomposition of A, A is m x n , v is a vector
  // result of SVD:  a = a * matrix(w) * v
  // return 0 if result bad (does not convengent)
  // return 1 if result usable

void seeMatrix(float **m, int hh, int ww, char *st=NULL);
  // print out matrix m for checking purpose
void seeDMatrix(double **m, int hh, int ww, char *st=NULL);
  // print out matrix m for checking purpose

void seeVector(float *v, int  ll, char *st=NULL);
  // print out vector v for checking purpose
void seeDVector(double *v, int  ll, char *st=NULL);
  // print out vector v for checking purpose

// the following for EigenValue/EigenVector decomposition
void gaucof(int n, float a[], float b[], float amu0, float x[], float w[]);
void tred2(float **a, int n, float d[], float e[]);
void tred2Double(double **a, int n, double d[], double e[]);
int tqli(float d[], float e[], int n, float **z); 
int tqliDouble(double d[], double e[], int n, double **z); 
// return 1 when good result; 0 when result meaningless

void eigsrt(float d[], float **v, int n);
void eigsrtDouble(double d[], double **v, int n);

// eigOfTensor is especially for Diffusion Tensor
// return 0 if result bad
int eigOfTensor(float **inMatrix,      //input squre symmetric matrix
	 float *eigenvalue,  		//allocated 1D space 4 eigenValues
	 float **eigenvector,  		//allocated 2D space 4 eigVector
	 int dim);          		//dimension
int eigOfTensorDouble(double **inMatrix,//input squre symmetric matrix
	 double *eigenvalue,  		//allocated 1D space 4 eigenValues
	 double **eigenvector,  	//allocated 2D space 4 eigVector
	 int dim);          		//dimension

// eigOfData() will call gaucof() & tred2() & tqli() & eigsrt()
// This is copied from Songyang, revised from "C Nemuric Receipt"
void eigOfData(float **feature,   	    // matrix to be decomposed
	 float *average, 	    // return value
	 float *eigenvlaue,         // return value
	 float **eigenvector,       // return value
	 int number_of_features,    // X-dim
	 int number_of_samples);    // Y-dim

// function eig () decompose 2D matrix into eigen vectors and eigen values
// "feature" matrix:
//    +-----------------------> number_of_features, e.g. X,Y,Z coordinates
//    |				(X)
//    |
//    |
//    |
//    |
//    |
//    |
//    V
//   number_of_samples 
//   (Y)
//
//  NOTE:
//   a = matrix (1,3, 1, 5)
//  Then:
//	 a[1][1] a[1][2] a[1][3] a[1][4] a[1][5] 
//	 a[2][1] a[2][2] a[2][3] a[2][4] a[2][5] 
//	 a[3][1] a[3][2] a[3][3] a[3][4] a[3][5] 
//
//

//================ By XDR--End ===============

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

void nrerror(char error_text[],int toExit=0);
float *vector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
float **matrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);
float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_vector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);

float** matrixInverse(float **aa, int n1, int n2); // matrix inverse
  //note, when matrix is not square, the inverse result c need to be
  // used for left-multiplying: I = c*aa, but I != aa*c
void gaussj(float **a, int n, float **b, int m); //gauss approach
void matrixXvec(float **m, int h, int w, float *v, int dimv, float* nv); //nv=m*v

#else /* ANSI */
/* traditional - K&R */

void nrerror();
void matrixXvec(); //nv=m*v
float *vector();
float **matrix();
float **submatrix();
float **convert_matrix();
float ***f3tensor();
double *dvector();
double **dmatrix();
int *ivector();
int **imatrix();
unsigned char *cvector();
unsigned long *lvector();
void free_vector();
void free_dvector();
void free_ivector();
void free_cvector();
void free_lvector();
void free_matrix();
void free_submatrix();
void free_convert_matrix();
void free_dmatrix();
void free_imatrix();
void free_f3tensor();
float** matrixInverse(); // matrix inverse
  //note, when matrix is not square, the inverse result c need to be
  // used for left-multiplying: I = c*aa, but I != aa*c
void gaussj(); //gauss approach

#endif /* ANSI */

#endif /* _NR_UTILS_H_ */
