#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <iostream.h>
//#include <strstream.h>
#include <strstream>
#include <fstream.h>

#include "Ellipsoid.h"

Ellipsoid::Ellipsoid()
{
   a[0]=a[1]=a[2]=0;
   aSqr[0]=aSqr[1]=aSqr[2]=0;
   center.x=center.y=center.z=0;
   O2N.loadIdentity();
}
Ellipsoid::Ellipsoid(float aa, float bb, float cc, struct Pt3d orig,
		Vector ax1, Vector ax2, Vector ax3)
{
  setParameters(aa, bb, cc, orig.x,orig.y,orig.z, ax1, ax2, ax3);
}

Ellipsoid::Ellipsoid(float aa, float bb, float cc, 
		float origX, float origY, float origZ,
		Vector ax1, Vector ax2, Vector ax3)
{
  setParameters(aa, bb, cc,origX,origY, origZ,ax1, ax2, ax3);
}

Ellipsoid::~Ellipsoid()
{
}

void Ellipsoid::setParameters(float aa, float bb, float cc, 
		float origX, float origY, float origZ,
		Vector ax1, Vector ax2, Vector ax3)
{
   a[0]=aa;a[1]=bb;a[2]=cc;
   V[0]=ax1;V[1]=ax2;V[2]=ax3;
   center.x=origX; center.y=origY; center.z=origZ;

   for (int i=0;i<3;i++)
     aSqr[i]=a[i]*a[i];
   O2N.newUVN(V[0],V[1],V[2]);
}

int Ellipsoid::inside(fPoint P)
{
  float x,y,z;
  P.getXYZ(x,y,z);
  return inside(x, y, z);
}

int Ellipsoid::inside(struct Pt3d P)
{
  return inside(P.x, P.y, P.z);
}

int Ellipsoid::inside(float x, float y, float z)
{ //judge if point (x,y,z) is inside this ellipsoid
  // return 1 if inside of this ellipsoid
  // return 0 if not

   float tx = x - center.x;
   float ty = y - center.y;
   float tz = z - center.z;

   fPoint pt(tx,ty,tz);
   pt = pt * O2N;
     //pt.outputs();
   double t = tx*tx/aSqr[0]+ty*ty/aSqr[1]+tz*tz/aSqr[2];
    //printf("\n JUDGE t=%8.6f\n",t);
    //cout <<" judge t = "<<t<<endl;
   if (t<=1.0) 
	return 1;
   else 
	return 0;
}

void Ellipsoid::outputs()
{
   cout<< "The center is: ("<<center.x<<","<<center.y<<","<<center.z<<")\n";
   cout<< "The axis are:";
   V[0].outputs(); V[1].outputs(); V[2].outputs();
   cout<< "The 3 axix lengthes are:";
   cout<< "   a="<<a[0];
   cout<< "   b="<<a[1];
   cout<< "   c="<<a[2]<<endl;
   cout<< "O2N matrix:";
   O2N.outputs();
}

void Ellipsoid::getTensor(Matrix & aTensor)
{// turn the Ellipsoid into a tensor, result in aTensor
    Matrix Vect(V[0].getX(),V[1].getX(),V[2].getX(), 0.0f,
		V[0].getY(),V[1].getY(),V[2].getY(), 0.0f,
		V[0].getZ(),V[1].getZ(),V[2].getZ(), 0.0f,
		       0.0f,	   0.0f,       0.0f, 1.0f);
	
    Matrix Lamb(a[0], 0.0f, 0.0f, 0.0f,
		0.0f, a[1], 0.0f, 0.0f,
		0.0f, 0.0f, a[2], 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f);

    Matrix tVect(Vect);
    tVect.transpose();

    aTensor = Vect * Lamb * tVect;
}

void Ellipsoid::getTensor(float *tens)
{// turn the Ellipsoid into a tensor, result in tens[6]
    Matrix aTensor;

    getTensor(aTensor);

    //float  Matrix.getValue(int row, int col);
    tens[0]=aTensor.getValue(0,0);
    tens[1]=aTensor.getValue(1,1);
    tens[2]=aTensor.getValue(2,2);
    tens[3]=aTensor.getValue(0,1);
    tens[4]=aTensor.getValue(0,2);
    tens[5]=aTensor.getValue(1,2);
}

Vector Ellipsoid::getPrimaryDir()
{
   float t1 = fabs(a[0]);
   float t2 = fabs(a[1]);
   float t3 = fabs(a[2]);
   int i;

   if ( t1 >= t2 && t1 >= t3){
	i = 0;	
   }else
   if ( t2 >= t3 && t2 >= t1){
	i = 1;	
   }else
   if ( t3 >= t2 && t3 >= t1){
	i = 2;	
   }
   //V[i].getXYZ(a,b,c);
   return V[i];
}


