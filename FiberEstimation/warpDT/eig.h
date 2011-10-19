#ifndef _EIG_H_
#define _EIG_H_

// #### this functions have been claimed in "nrutil.h" ####

void gaucof(int n, float a[], float b[], float amu0, float x[], float w[]);
void tred2(float **a, int n, float d[], float e[]);
int tqli(float d[], float e[], int n, float **z);
void eigsrt(float d[], float **v, int n);

// eigOfTensor is especially for Diffusion Tensor
// return 0 if result bad
int eigOfTensor(float **inMatrix,      //input squre symmetric matrix
	 float *eigenvalue,  		//allocated 1D space 4 eigenValues
	 float **eigenvector,  		//allocated 2D space 4 eigVector
	 int dim);          		//dimension

// calculate FA from 3 eigenvalues
double getFA(float *eigenvalue);

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

#endif //_EIG_H_
