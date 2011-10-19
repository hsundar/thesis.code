#include <stdio.h>
#include <math.h>

#include "Global.h"

#include "nrutil.h"
#include "eig.h"

// Decompose matrix "feature"
// to eigenvalues and eigenvectors
//
// NOTE: all arrays start from 1 rather than 0
//
//
//

void eigOfData(float **feature, 
	 float *average, 
	 float *eigenvalue, 
	 float **eigenvector, 
	 int number_of_features,    // X
	 int number_of_samples)     // Y
{// feature is a group of data, this function first 
 // compute the covariance of this group of data,
 // based on which, the eigenVectors and eigenValues
 // are calculated.
	int i,j,k;
	float *e;
	float **cov;
	float *sum;

        for (i=1;i<=number_of_features;i++){
	      average[i] =0; 
	      eigenvalue[i]=0;
   	}

	for (i=1;i<=number_of_features;i++){
	  for (j=1;j<=number_of_features;j++){
	    eigenvector[i][j]=0;
	  }
 	}

	for (i=1;i<=number_of_samples;i++)
	for(j=1;j<=number_of_features;j++)
	   average[j] += feature[i][j];

	for (i=1;i<=number_of_features;i++)
		average[i]=average[i]/number_of_samples;

 	seeVector(average, 3, "AVERAGE:");
	if (number_of_samples>=number_of_features)
	{
		for(k=1;k<=number_of_samples;k++)
		for (i=1;i<=number_of_features;i++)
		for(j=1;j<=number_of_features;j++)
		   eigenvector[i][j] += (feature[k][i]-average[i])
					*(feature[k][j]-average[j]); 



		for(i=1;i<=number_of_features;i++)
		for(j=1;j<=number_of_features;j++)
		eigenvector[i][j] = eigenvector[i][j]/(number_of_samples-1);


		e=vector(1,number_of_features);

		tred2(eigenvector,number_of_features,eigenvalue,e);
		tqli(eigenvalue,e,number_of_features,eigenvector);
		eigsrt(eigenvalue, eigenvector,number_of_features);

		free_vector(e,1,number_of_features);
	}
	else
	{
		cov=matrix(1,number_of_samples,1,number_of_samples);
		sum=vector(1,number_of_samples);
		e=vector(1,number_of_samples);

		for(k=1;k<=number_of_features;k++)
		for(i=1;i<=number_of_samples;i++)
		for(j=1;j<=number_of_samples;j++)
	 	   cov[i][j] += (feature[i][k]-average[k])
				*(feature[j][k]-average[k]);

		for(i=1;i<=number_of_samples;i++)
		for(j=1;j<=number_of_samples;j++)
		   cov[i][j] /=(number_of_samples-1);


		tred2(cov,number_of_samples,eigenvalue,e);
		tqli(eigenvalue,e,number_of_samples,cov);
		eigsrt(eigenvalue, cov,number_of_samples);
		
	
			
		for(k=1;k<=number_of_samples;k++)
		for(i=1;i<=number_of_features;i++)
		for(j=1;j<=number_of_samples;j++)
		  eigenvector[i][j] += (feature[k][i]-average[i])*cov[k][j];


		for(i=1;i<=number_of_features;i++)
		for(j=1;j<=number_of_samples;j++)
		  sum[j]+=eigenvector[i][j]*eigenvector[i][j];
		
		for(j=1;j<=number_of_samples;j++)
		  sum[j]=sqrt(sum[j]);

		for(i=1;i<=number_of_features;i++)
		for(j=1;j<=number_of_samples;j++)
		  eigenvector[i][j]=eigenvector[i][j]/sum[j];

		free_vector(e,1,number_of_samples);
		free_vector(sum,1,number_of_samples);
		free_matrix(cov,1,number_of_samples,1,number_of_samples);
	}
		
}

double getFA(float *eVal)
{
  double W123 = eVal[1]*eVal[1] + eVal[2]*eVal[2] + eVal[3]*eVal[3];
  double weight,ww;

  if (W123 < ZERO){
	weight = 0.0;
  }else {
	ww = fabs(eVal[1]*eVal[2])+fabs(eVal[2]*eVal[3])+fabs(eVal[3]*eVal[1]);
	weight =  1.0 - ww/W123;  // could be negative, due to the computer accuracy
	if (weight < 0) 
		weight = 0;
	else 
	    weight = sqrt (weight);
  }
  return weight;
}

int eigOfTensor(float **inMatrix,      //input matrix
	 float *eigenvalue, 		//allocated 1D space 4 eigenValues
	 float **eigenvector, 		//allocated 2D space 4 eigVector
	 int dim)     			//dimension
{//input must be a squre matrix
 //this function is specially for Tensor, 
 //i.e. symmetric matrix
 //return 0 if result bad
 //else return 1
	int i,j;
	float *e;
	int res = 1;
	double maxval=fabs(inMatrix[1][1]);

	for(i=1;i<=dim;i++){
	  eigenvalue[i]=0;  // clear to be 0.0
	  for(j=1;j<=dim;j++){
	    eigenvector[i][j]=0;
	    if (fabs(inMatrix[i][j])> maxval)  //pick out the largest as Scale
		maxval=fabs(inMatrix[i][j]);
	  }
	}
        if (maxval == 0){
	  for(j=1;j<=dim;j++)
	     eigenvector[j][j]=1;
	  res = 1;
	}else{
	  for(i=1;i<=dim;i++)
	  for(j=1;j<=dim;j++){
	    eigenvector[i][j]=inMatrix[i][j]/maxval; //normalize the matrix
	  }	

	  e=vector(1,dim);

	  //seeMatrix(eigenvector, dim,dim, "(eigenvector)before tred2():\n");
	  //seeVector(eigenvalue,dim," (eigenvalue) before tred2():\n");
	    //seeVector(e,dim," (e) before tred2():\n");
	  tred2(eigenvector,dim,eigenvalue,e);
	    //seeMatrix(eigenvector, dim,dim, "(eigenvector)after tred2():\n");
	    //seeVector(eigenvalue,dim," (eigenvalue) after tred2():\n");
	    //seeVector(e,dim," (e) after tred2():\n");
	  if (0==tqli(eigenvalue,e,dim,eigenvector)){
		//if (k++ <10)
		  //seeMatrix(inMatrix,dim,dim,"bad Matrix:");
		res = 0;
	  }
	  eigsrt(eigenvalue, eigenvector,dim);

	  free_vector(e,1,dim);
	}
	for(j=1;j<=dim;j++)
	   eigenvalue[j] *= maxval;

	return res;
}
int eigOfTensorDouble(double **inMatrix,      //input matrix
	 double *eigenvalue, 		//allocated 1D space 4 eigenValues
	 double **eigenvector, 		//allocated 2D space 4 eigVector
	 int dim)     			//dimension
{//input must be a squre matrix
 //this function is specially for Tensor, 
 //i.e. symmetric matrix
 //return 0 if result bad
 //else return 1
	int i,j;
	//int k=0;
	double *e;
	int res = 1;
	double maxval=fabs(inMatrix[1][1]);

	for(i=1;i<=dim;i++){
	  eigenvalue[i]=0;  // clear to be 0.0
	  for(j=1;j<=dim;j++){
	    eigenvector[i][j]=0;
	    if (fabs(inMatrix[i][j])> maxval)  //pick out the largest as Scale
		maxval=fabs(inMatrix[i][j]);
	  }
	}
        if (maxval == 0){
	  for(j=1;j<=dim;j++)
	     eigenvector[j][j]=1;
	  res = 1;
	}else{
	  for(i=1;i<=dim;i++)
	  for(j=1;j<=dim;j++){
	    eigenvector[i][j]=inMatrix[i][j]/maxval; //normalize the matrix
	  }	

	  e=dvector(1,dim);

	    //seeMatrix(eigenvector, dim,dim, "(eigenvector)before tred2():\n");
	    //seeVector(eigenvalue,dim," (eigenvalue) before tred2():\n");
	  tred2Double(eigenvector,dim,eigenvalue,e);
	    //seeDMatrix(eigenvector, dim,dim, "(eigenvector)after tred2():\n");
	    //seeDVector(eigenvalue,dim," (eigenvalue) after tred2():\n");
	  if (0==tqliDouble(eigenvalue,e,dim,eigenvector)){
		//if (k++ <10)
		  //seeMatrix(inMatrix,dim,dim,"bad Matrix:");
		res = 0;
	  }
	  eigsrtDouble(eigenvalue, eigenvector,dim);

	  free_dvector(e,1,dim);
	}
	for(j=1;j<=dim;j++)
	   eigenvalue[j] *= maxval;

	return res;
}

//--------------------------------------------------------
//   The following part was from gaucof.C
//--------------------------------------------------------
#define NRANSI

void gaucof(int n, float a[], float b[], float amu0, float x[], float w[])
{
	void eigsrt(float d[], float **v, int n);
	int tqli(float d[], float e[], int n, float **z);
	int i,j;
	float **z;

	z=matrix(1,n,1,n);
	for (i=1;i<=n;i++) {
		if (i != 1) b[i]=sqrt(b[i]);
		for (j=1;j<=n;j++) z[i][j]=(float)(i == j);
	}
	tqli(a,b,n,z);
	eigsrt(a,z,n);
	for (i=1;i<=n;i++) {
		x[i]=a[i];
		w[i]=amu0*z[1][i]*z[1][i];
	}
	free_matrix(z,1,n,1,n);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software . */


//--------------------------------------------------------
//   The following part was from tred2.C
//--------------------------------------------------------

void tred2(float **a, int n, float d[], float e[])
{
	int l,k,j,i;
	float scale,hh,h,g,f;

	for (i=n;i>=2;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 1) {
			for (k=1;k<=l;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i]=a[i][l];
			else {
				for (k=1;k<=l;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=0.0;
				for (j=1;j<=l;j++) {
					a[j][i]=a[i][j]/h;
					g=0.0;
					for (k=1;k<=j;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<=l;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h;
					f += e[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=1;j<=l;j++) {
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=1;k<=j;k++)
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}
	d[1]=0.0;
	e[1]=0.0;
	/* Contents of this loop can be omitted if eigenvectors not
			wanted except for statement d[i]=a[i][i]; */
	for (i=1;i<=n;i++) {
		l=i-1;
		if (d[i]) {
			for (j=1;j<=l;j++) {
				g=0.0;
				for (k=1;k<=l;k++)
					g += a[i][k]*a[k][j];
				for (k=1;k<=l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i];
		a[i][i]=1.0;
		for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
	}
}
void tred2Double(double **a, int n, double d[], double e[])
{
	int l,k,j,i;
	double scale,hh,h,g,f;

	for (i=n;i>=2;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 1) {
			for (k=1;k<=l;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i]=a[i][l];
			else {
				for (k=1;k<=l;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=0.0;
				for (j=1;j<=l;j++) {
					a[j][i]=a[i][j]/h;
					g=0.0;
					for (k=1;k<=j;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<=l;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h;
					f += e[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=1;j<=l;j++) {
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=1;k<=j;k++)
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}
	d[1]=0.0;
	e[1]=0.0;
	/* Contents of this loop can be omitted if eigenvectors not
			wanted except for statement d[i]=a[i][i]; */
	for (i=1;i<=n;i++) {
		l=i-1;
		if (d[i]) {
			for (j=1;j<=l;j++) {
				g=0.0;
				for (k=1;k<=l;k++)
					g += a[i][k]*a[k][j];
				for (k=1;k<=l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i];
		a[i][i]=1.0;
		for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
	}
}

/* (C) Copr. 1986-92 Numerical Recipes Software . */


//--------------------------------------------------------
//   The following part was from eigsrt.C
//--------------------------------------------------------
void eigsrt(float d[], float **v, int n)
{
	int k,j,i;
	float p;

	for (i=1;i<n;i++) {
		p=d[k=i];
		for (j=i+1;j<=n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=1;j<=n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}
void eigsrtDouble(double d[], double **v, int n)
{
	int k,j,i;
	double p;

	for (i=1;i<n;i++) {
		p=d[k=i];
		for (j=i+1;j<=n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=1;j<=n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}

//--------------------------------------------------------
//   The following part was from tqli.C
//--------------------------------------------------------
#define NRANSI
int tqli(float d[], float e[], int n, float **z)
{// IMPORTANT:
 //   the input e[] should already have been normalized to be near 1.0
 //   if the e[] is like 1e-6, this procedure would be a problem

	float pythag(float a, float b);
	int m,l,iter,i,k;
	float s,r,p,g,f,c,b;
    	//float dd;

	    //seeVector(e,n,"  $$$$$ (e)intqli():\n");

	for (i=2;i<=n;i++) e[i-1]=e[i];
	e[n]=0.0;
	for (l=1;l<=n;l++) {
	  iter=0;
	  do {
		//printf("\n     === l = %d ******", l);
	     for (m=l;m<=n-1;m++) {
		 //printf("\n   (d[m])d[%d]=%7.5f,(d[m+1]) d[%d+1]=%7.5f",
			//m,d[m],m,d[m+1]); //XDR
	  	//dd=fabs(d[m])+fabs(d[m+1]);
		  //printf("\n    dd = d[m]+d[m+1] = %27.26f",dd); //XDR
		  //printf("\n    m=%d,l=%d(fabs(e[m])+dd)= fabs(%20.19f)+%5.3f=%27.26f",m,l,e[m],dd,fabs(e[m])+dd); //XDR
		  //printf("\n    fabs(e[m])%37.36f",e[m]); //XDR
		//if ((float)(fabs(e[m])+dd) == dd){
		if (fabs(e[m]) <ZERO) { 
		   //printf("\n     EQUAL\n");
		   break;
		}//else{
		   //printf("\n       Not Equal\n");
		//}
	     }		
		//printf("\ntqli(): m=%d, l = %d",m,l);//XDR debug

	     if (m != l) {
		if (iter++ == 30) {
			nrerror("Too many iterations in tqli");
			printf("\nHey: too many iterations in tqli");
			return 0; // bad result
		}
		g=(d[l+1]-d[l])/(2.0*e[l]);
		r=pythag(g,1.0);
		g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
		s=c=1.0;
		p=0.0;
		for (i=m-1;i>=l;i--) {
			f=s*e[i];
			b=c*e[i];
			e[i+1]=(r=pythag(f,g));
			if (r == 0.0) {
			  d[i+1] -= p;
			  e[m]=0.0;
			  break;
			}
			s=f/r;
			c=g/r;
			g=d[i+1]-p;
			r=(d[i]-g)*s+2.0*c*b;
			d[i+1]=g+(p=s*r);
			g=c*r-b;
			for (k=1;k<=n;k++) {
			  f=z[k][i+1];
			  z[k][i+1]=s*z[k][i]+c*f;
			  z[k][i]=c*z[k][i]-s*f;
			}
		}
		if (r == 0.0 && i >= l) continue;
		d[l] -= p;
		e[l]=g;
		e[m]=0.0;
		}
	    } while (m != l);
	}
	return 1; //good result
}
int tqliDouble(double d[], double e[], int n, double **z)
{// IMPORTANT:
 //   the input e[] should already have been normalized to be near 1.0
 //   if the e[] is like 1e-6, this procedure would be a problem

	double pythagDouble(double a, double b);
	int m,l,iter,i,k;
	double s,r,p,g,f,c,b;
	double dd;
	    //seeDVector(e,n,"  $$$$$ (e)intqliDouble():\n");
	double volatile temp;

	for (i=2;i<=n;i++) e[i-1]=e[i];
	e[n]=0.0;
	for (l=1;l<=n;l++) {
	  iter=0;
	  do {
		//printf("\n     === l = %d ******", l);
	     for (m=l;m<=n-1;m++) {
		 //printf("\n   (d[m])d[%d]=%7.5f,(d[m+1]) d[%d+1]=%7.5f",
			//m,d[m],m,d[m+1]); //XDR
	  	dd=fabs(d[m])+fabs(d[m+1]);
		  //printf("\n    dd = d[m]+d[m+1] = %27.26f",dd); //XDR
		  //printf("\n    (fabs(e[m])+dd)= fabs(%10.9f)+%5.3f=%27.26f",e[m],dd,(float)(fabs(e[m])+dd)); //XDR
//		if ((double)(fabs(e[m])+dd) == dd) {
//		if ((fabs(e[m])+dd) == dd) { //judge if( fabs(e[m]) == 0) 

	       if( fabs(e[m]) < ZERO ) {
	       //temp = fabs(e[m]) + dd;
	       //if(temp == dd) {
		   //printf("\n     EQUAL");
		   break;
		}//else{
		   //printf("\n       Not Equal");
		//}
	     }
		//printf("\n #### tqliDouble(): m=%d, l = %d",m,l);//XDR debug
	     if (m != l) {
		if (iter++ == 30) {
		  //	nrerror("Too many iterations in tqliDouble");
		  //	printf("\nHey, too many iterations in tqliDouble");
		  //	printf("\nm=%d, l = %d\n\n",m,l);
			return 0; // bad result
		}
		g=(d[l+1]-d[l])/(2.0*e[l]);
		r=pythagDouble(g,1.0);
		g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
		s=c=1.0;
		p=0.0;
		for (i=m-1;i>=l;i--) {
			f=s*e[i];
			b=c*e[i];
			e[i+1]=(r=pythagDouble(f,g));
			if (r == 0.0) {
			  d[i+1] -= p;
			  e[m]=0.0;
			  break;
			}
			s=f/r;
			c=g/r;
			g=d[i+1]-p;
			r=(d[i]-g)*s+2.0*c*b;
			d[i+1]=g+(p=s*r);
			g=c*r-b;
			for (k=1;k<=n;k++) {
			  f=z[k][i+1];
			  z[k][i+1]=s*z[k][i]+c*f;
			  z[k][i]=c*z[k][i]-s*f;
			}
		}
		if (r == 0.0 && i >= l) continue;
		d[l] -= p;
		e[l]=g;
		e[m]=0.0;
		}
	    } while (m != l);
	}
	return 1; //good result
}
#undef NRANSI


