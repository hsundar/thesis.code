/* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
   utility file nrutil.c.  Do not confuse this file with the same-named
   file nrutil.c that is supplied in the same subdirectory or archive
   as the header file nrutil.h.  *That* file contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only ANSI C.               */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#include <string.h>
#include <math.h>

#include <iostream.h>
//#include <strstream.h>
#include <strstream>
#include <fstream.h>

#include "nrutil.h"

#define NR_END 1
#define FREE_ARG char*
static float tempr;
#define SWAP(a,b) {tempr=(a);(a)=(b);(b)=tempr;}

void nrerror(char error_text[],int toExit)
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	if (toExit)
	    exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;
	int i;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) {
	   nrerror("allocation failure in vector()",1);
	   exit(1);
	}

	for (i = nl; i<=nh;i++)
	  v[i]=0;

	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;
	int i;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) {
	   nrerror("allocation failure in ivector()",1);
	   exit(0);
	}
	for (i = nl; i<=nh;i++)
	  v[i]=0;

	return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;
	int i;

	v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) {
	   nrerror("allocation failure in cvector()",1);
	   exit(0);
	}

	for (i = nl; i<=nh;i++)
	  v[i]=0;

	return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;
	int i;

	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) {
	   nrerror("allocation failure in lvector()",1);
	   exit(0);
	}
	for (i = nl; i<=nh;i++)
	  v[i]=0;

	return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;
 	int i;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) {
	   nrerror("allocation failure in dvector()",1);
	   exit(0);
	}
	for (i = nl; i<=nh;i++)
	  v[i]=0;
	return v-nl+NR_END;
}

//  NOTE:
//   a = matrix (1,3, 1, 5)
//  Then:
//	 a[1][1] a[1][2] a[1][3] a[1][4] a[1][5] 
//	 a[2][1] a[2][2] a[2][3] a[2][4] a[2][5] 
//	 a[3][1] a[3][2] a[3][3] a[3][4] a[3][5] 
//
//

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m){
	  nrerror("allocation failure 1 in matrix()",1);
	  exit(0);
  	}
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()",1);
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        for (i = nrl; i <= nrh; i++)
        for (j = ncl; j <= nch; j++)
	   m[i][j]=0;

	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) {
	   nrerror("allocation failure 1 in matrix()",1);
	   exit(0);
	}
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()",1);
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        for (i = nrl; i <= nrh; i++)
        for (j = ncl; j <= nch; j++)
	   m[i][j]=0;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, j, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) {
	   nrerror("allocation failure 1 in matrix()",1);
	   exit(0);
	}
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()",1);
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        for (i = nrl; i <= nrh; i++)
        for (j = ncl; j <= nch; j++)
	   m[i][j]=0;

	/* return pointer to array of pointers to rows */
	return m;
}

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	float **m;

	/* allocate array of pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in submatrix()",1);
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in convert_matrix()",1);
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;

	/* allocate pointers to pointers to rows */
	t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()",1);
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()",1);
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()",1);
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}
/* (C) Copr. 1986-92 Numerical Recipes Software . */

void mulMatrixABtt(float **A, float **B, float **C, int row, int col)
{ //this is for debugging only
  //C=A*B', A is (row x col, B too
   int i,j,k;
   double t;

		 cout << "mulMatrixABtt()DDDD222: row="<<row<<"; col="<<col<<endl;
		seeMatrix(A,3,col, "matrix AA2:");
		seeMatrix(B,3,col, "matrix BB2:");
		seeMatrix(C,3,3, "matrix CC2:");

   for ( i = 1; i < row+1; i++)
   for ( j = 1; j < row+1; j++){
     t = 0;
     for (k = 1; k < 1+col; k++)
	t += A[i][k]*B[j][k];
     C[i][j] = t;
   }
}

void mulMatrixABt(float **A, float **B, float **C, int row, int col)
{//C=A*B', A is (row x col), B too
   int i,j,k;
   double t;
   if (col < 1){
     for ( i = 1; i < row+1; i++)
     for ( j = 1; j < row+1; j++)
	C[i][j]=0;
   }
   for ( i = 1; i < row+1; i++)
   for ( j = 1; j < row+1; j++){
     t = 0;
     for (k = 1; k < 1+col; k++)
	t += A[i][k]*B[j][k];
     C[i][j] = t;
   }
}
void mulMatrixAB(float **A, float **B, float **C, int dim)
{//C=A*B, A & B & C are (dim x dim)
   int i,j,k;
   double t;

   for ( i = 1; i < dim+1; i++)
   for ( j = 1; j < dim+1; j++){
     t = 0;
     for (k = 1; k < 1+dim; k++)
	t += A[i][k]*B[k][j];
     C[i][j] = t;
   }
}
void mulMatrixAB(float **A, float **B, float **C, int row, int col)
// C=A*B, 
// A (row x col)  B (col x row) C are (row x row)
{
   int i,j,k;
   double t;

   for ( i = 1; i < row+1; i++)
   for ( j = 1; j < row+1; j++){
     t = 0;
     for (k = 1; k < 1+col; k++)
	t += A[i][k]*B[k][j];
     C[i][j] = t;
   }
}
void mulMatrixAB(float **A, int r1, int c1, float **B, int r2, int c2, float **C)
// C=A*B, 
// A (r1 x c1)  B (r2 x c2) C are (r1 x c2)
{
   int i,j,k;
   double t;

   if (c1 != r2)
     nrerror("mulMatrixAB(float **A, int r1, int c1, float **B, int r2, int c2, float **C): c1 != r2");

   for ( i = 1; i < r1+1; i++)
   for ( j = 1; j < c2+1; j++){
     t = 0;
     for (k = 1; k < 1+c1; k++)
	t += A[i][k]*B[k][j];
     C[i][j] = t;
   }
}
void seeMatrix(float **m, int hh, int ww, char *st)
{// print out matrix m
   int i, j;
   cout << st<<endl;
   for (i = 1; i <= hh; i++){
     cout <<"["<<i<<"]: \n";
     for (j = 1; j <= ww; j++){
	cout << "\t"<<m[i][j];
     }
     cout << "  ;"<<endl;
   }
   cout <<"------------" << endl;
}
void seeDMatrix(double **m, int hh, int ww, char *st)
{// print out matrix m
   int i, j;
   cout << st<<endl;
   for (i = 1; i <= hh; i++){
     cout <<"["<<i<<"]: \n";
     for (j = 1; j <= ww; j++){
	cout << "\t"<<m[i][j];
     }
     cout << "  ;"<<endl;
   }
   cout <<"------------" << endl;
}

void seeVector(float *v, int  ll, char *st)
{//print out vector v
   int i;
   cout << st<<endl;
   for (i = 1; i <= ll; i++){
	cout << "\t"<<v[i];
   }
   cout << endl;
   cout <<"------------" << endl;
}

void seeDVector(double *v, int  ll, char *st)
{//print out vector v
   int i;
   cout << st<<endl;
   for (i = 1; i <= ll; i++){
	printf("\t %40.39f ", v[i]);
	//cout << "\t"<<v[i];
   }
   printf("\n---------------------------\n");
   //cout << endl;
   //cout <<"------------" << endl;
}

void gaussj(float **a, int n, float **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	float big,dum,pivinv;

	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) nrerror("gaussj: Singular Matrix-1");
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix-2");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}

float** matrixInverse(float **aa, int n1, int n2)
// **aa is a matrix for calculating inverse matrix
// **aa dim=(row: 1..n1, column: 1..n2)
// **b dim=(row: 1..n1, column: 1..m)
// result in c, return to caller, caller is responsible to free memory
{
  int    i, j, m; //, m, n1, n2 ;
  float  **a, **b,**c; //c is for result
  // a[i][j]: row i; column j
  float **At, **AtA, **AtAInv ;

  //n1 = A->height ;
  //n2 = A->width ;
  m = 1 ;

  a = matrix(1,n1,1,n2);
  /* set matrices: a <- A */
  for(i=1; i<=n1; i++)
   for(j=1; j<=n2; j++) 
      a[i][j] = aa[i][j];

  if(n1==n2){
      /* apply for  memory */
      /* Notice here the beginning index for u (v or w) is 1, not 0.  */
      c = matrix(1,n1,1,n2);
      b = matrix(1,n1,1,m);
      
      /*cheat b*/
     for(i=1; i<=n1; i++)
       for(j=1; j<=m; j++) 
	  b[i][j] = a[i][j]; // A->data[i-1][j-1] ;

     gaussj(a, n1, b, m) ; /* get inverse of a */

      /* set matrices: A <- a */
      for(i=1; i<=n1; i++)
	for(j=1; j<=n2; j++) 
	  c[i][j]=a[i][j] ;
	  //B->data[i-1][j-1] = a[i][j] ;

     free_matrix(b,1,n1,1,m);
     return c;
  }else { // n1 != n2 result in c, and c*aa ==> I (inverse to be at left)
      c = matrix(1,n2,1,n1);

      At=matrix(1,n2,1,n1); cout <<" ---1111---"<<endl;
      AtA=matrix(1,n2,1,n2);cout <<" ---2222---"<<endl;
      //AtAInv=matrix(1,n2,1,n2);

      /* At: transpose of matrix a */
      for(i=1; i<=n2; i++)
	for(j=1; j<=n1; j++) 
	  At[i][j]=aa[j][i];
      seeMatrix(At, n2, n1, "Matrix At:");

      mulMatrixAB(At,a,AtA,n2,n1); // AtA=At*A;
      seeMatrix(AtA, n2, n2, "Matrix AtA:");

      AtAInv=matrixInverse(AtA, n2, n2); //AtAInv is the inverse of AtA
      seeMatrix(AtAInv, n2, n2, "Matrix AtAInv:");

      mulMatrixAB(AtAInv,n2,n2,At,n2,n1, c); 
      seeMatrix(c, n2, n1, "Matrix c:");

      free_matrix(At,1,n2,1,n1) ;  
      free_matrix(AtA,1,n2,1,n2) ;
      free_matrix(AtAInv,1,n2,1,n2) ;
  }
  free_matrix(a,1,n1,1,n2);
  return c;
}

void matrixXvec(float **m, int h, int w, float *v, int dimv, float* nv)
//nv=m*v
{
   int i,j;

   if (w !=dimv){
     cout << "matrixXvec(): width of matrix != dim of vector "<<endl;
     exit(0);
   }
   for (i=1; i<=h; i++){
     nv[i]=0;
     for (j=1; j<=w; j++)
	nv[i] += m[i][j]*v[j];
   }
}
