#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "Basics.h"
#include "GPrimitive.h"
#include <iostream.h>
//#include <strstream.h>
#include <strstream>
#include <fstream.h>

#include <sys/times.h>
#include <time.h>


//--------------------------------------------
//class fPoint members
//
//
//------------------------------------------
fPoint::fPoint()
{
  vm[0] = vm[1] = vm[2] = 0.0; 
  vm[3] = 1.0F; 
  formulated = TRUE;
}
fPoint::fPoint(struct Pt3d vox)
{
  vm[0] = vox.x;
  vm[1] = vox.y;
  vm[2] = vox.z;
  vm[3] = 1.0;
  formulated = TRUE;
}
fPoint::fPoint(int a, int b, int c, int d)
{
  vm[0] = (float)a;
  vm[1] = (float)b;
  vm[2] = (float)c; 
  vm[3] = (float)d; 
  if (fabs(d - 1.0f)<ZERO) formulated = TRUE; 
  else formulated = FALSE;
}

fPoint::fPoint(float a, float b, float c, float d)
{
  vm[0] = a;
  vm[1] = b;
  vm[2] = c;
  vm[3] = d;
 
  if (fabs(d-1.0f)<ZERO) formulated = TRUE; 
  else formulated = FALSE; 
}

fPoint::fPoint(const fPoint &pt)
{
  int i;
  for (i=0; i<4; i++) vm[i] = pt.vm[i];
  formulated = pt.formulated;
}

fPoint::~fPoint(void)
{
}

BOOL fPoint::isUndefined()
{
  if (vm[3]==0.0)
	return TRUE;
  else 
	return FALSE;
}

BOOL fPoint::equal(fPoint p)
{
  for (int i=0; i<4;i++)
    if (vm[i]!= p.vm[i]) return FALSE;
  return TRUE;
}
void fPoint::setX(int a)
{
  vm[0] = (float)a;
}

void fPoint::setX(float a)
{
  vm[0] = a;
}

void fPoint::setY(int b)
{
  vm[1] = (float)b;
}
void fPoint::setY(float b)
{
  vm[1] = b;
}

void fPoint::setZ(int c)
{
  vm[2] = (float)c;
}

void fPoint::setZ(float c)
{
  vm[2] = c;
}
 
void fPoint::setAA(int d) 
{ 
  vm[3] = (float)d; 
 
  if (fabs(d - 1.0f)<ZERO) formulated = FALSE; 
  else formulated = TRUE; 

} 
void fPoint::setAA(float d) 
{ 
  vm[3] = d; 
   
  if (fabs(d-1.0f)<ZERO) formulated = FALSE; 
  else formulated = TRUE; 
 
} 

void fPoint::setXYZ(int a,int b, int c, int d)
{
  vm[0] = (float)a; vm[1] = (float)b; vm[2] = (float)c;	 vm[3] = (float)d; 
  
 if (fabs(d-1.0f)<ZERO) formulated = FALSE; 
 else formulated = TRUE; 

}

void fPoint::setXYZ(float a,float b, float c, float d)
{
  vm[0]= a; vm[1] = b; vm[2] = c; vm[3] = d; 
 if (fabs(d-1.0f)<ZERO) formulated = FALSE; 
 else formulated = TRUE; 

}

float fPoint::getX(void)
{ 
  if (vm[3]==0) return 0;
  else return vm[0]/vm[3];
}


float fPoint::getY(void)
{ 
  if (vm[3]==0) return 0;
  else return vm[1]/vm[3];
}

float fPoint::getZ(void)
{  
  if (vm[3]==0) return 0;
  else return vm[2]/vm[3];
}
 
float fPoint::distance(fPoint pt)
{
 int i;
  float t=0.0f,tt;

  formulate();
  pt.formulate();
  for (i=0;i<3;i++){
     tt = fabs(vm[i]-pt.vm[i]);
     t += tt * tt; 
  }
  //t=_sqrt_s(t);//fabs(vm[3]);
  t=sqrt(t);
  return t;
}

float fPoint::distance(float a, float b, float c, float d)
{
  fPoint pt(a,b,c,d);
  return distance(pt);
}

fPoint fPoint::negative(void)
{
  fPoint tmp;
  int i;
  for (i=0;i<3;i++) tmp.vm[i] = -vm[i];
  tmp.vm[3] = vm[3];
  tmp.formulated = formulated;
  //printf("\nthis=");//outputs();
  //printf("tmp:"); tmp.outputs();
  return tmp;
}

fPoint fPoint::operator-(void)
{
  return negative();
}

fPoint fPoint::operator*(Matrix nM)
{
   int row,col;
   double t[4];

   for (col=0;col<dm;col++){
        t[col] = 0.0f;
        for (row=0;row<dm;row++)
                t[col] += vm[row] * nM.mm[row][col];
  }
  fPoint pp((float)(t[0]),(float)(t[1]),(float)(t[2]),(float)(t[3]));
  //printf("\n fPoint::*:\n   t=(%4.3f,%4.3f,%4.3f,%4.3f),",t[0],t[1],t[2],t[3]);
  //if (pp.formulated) printf("--Is formulated");else printf("---Not Formulated");
  pp.formulate();
  return Vector(pp);
}

fPoint fPoint::operator-(fPoint nP)
{
  return operator+(-nP);
}

fPoint fPoint::operator+(fPoint nP)
{
   fPoint tmp;
   tmp = *this;
   tmp.formulate(); nP.formulate();
   tmp.vm[0] += nP.vm[0];
   tmp.vm[1] += nP.vm[1];
   tmp.vm[2] += nP.vm[2];
   return tmp;
}
fPoint fPoint::operator*(float ss)
{
   fPoint tmp = *this;
   tmp.scale(ss,ss,ss);
   return tmp;
}

void fPoint::scale(float sx,float sy, float sz)
{
   formulate();
   vm[0] *= sx; vm[1] *= sy; vm[2] *= sz;
}
void fPoint::deScale(float sx,float sy, float sz)
{
   formulate();
   vm[0] /= sx; vm[1] /= sy; vm[2] /= sz;
}

   
void fPoint::formulate(void) 
{ 
  if (formulated) return; 
  else { 
	if (vm[3]==0) { 
	  vm[0]=vm[1]=vm[2]=0.0F; 
	  vm[3]=1.0F; 
	} 
	else {  
	  vm[0] /= vm[3]; vm[1]/= vm[3]; vm[2] /= vm[3];
	  vm[3] = 1.0f;
	} 
	formulated = TRUE; 
  } 
} 

void fPoint::getXYZ(float& a,float& b, float& c)
{ 
	formulate(); 
	a=vm[0];b=vm[1];c=vm[2]; 
}

fPoint fPoint::proportionPoint(fPoint B, float t1, float t2)
{
  fPoint tmp;

  if ((t1+t2-0)<ZERO){// r= -1, insection point at infinative place
     printf ("\n ERROR: Image::proportion: viewpoint on object.");
     exit(0);
  }
  if ((t2-0)<ZERO){//it is the same as viewpoint
     printf("\n ERROR: Image::dotOnImagePlane: viewpoint on view Plane");
     return B;
  }
  float r = t1 / t2;
 
  formulate(); B.formulate();
  int i;
  for (i=0;i<3; i++)
     tmp.vm[i] = (vm[i] + r * B.vm[i])/(1+r);
  tmp.vm[3] = 1.0f;
  return tmp;
}

void fPoint::outputs(void)
{
  cout << endl << "the point or vector is:(" << vm[0];
  cout << "," <<vm[1] << "," <<vm[2] <<"," <<vm[3]<<")" << endl;
 //  printf("\n the point or vector is:(%5.4f,%5.4f,%5.4f,%5.4f)",vm[0],\
  //	vm[1],vm[2],vm[3]);
}

//----------------------------------------------------------- 
// class Vector 
// 
// 
// 
 
Vector::Vector() 
  :fPoint() 
{} 
 
Vector::Vector(fPoint p) 
  :fPoint(p) 
{} 
 
Vector::Vector(fPoint p1, fPoint p2) 
{ 
   setValue(p1,p2);
} 
Vector::Vector(int x, int y, int z)
	:fPoint(x,y,z)
{}

Vector::Vector(float x, float y, float z)
	:fPoint(x,y,z)
{}

Vector::Vector(const Vector &V) 
{
  int i;
  for (i=0;i<4;i++) vm[i] = V.vm[i];
  formulated = V.formulated;
}

Vector::~Vector() 
{} 
 
void Vector::setValue(fPoint p)
{
  float xx,yy,zz;
  p.getXYZ(xx,yy,zz);
  vm[0] = xx;
  vm[1] = yy;
  vm[2] = zz;
  vm[3] = 1.0;
}

void Vector::setValue(fPoint p1, fPoint p2)
{
  p1.formulate();p2.formulate();
  vm[0] = p2.getX()-p1.getX();
  vm[1] = p2.getY()-p1.getY();
  vm[2] = p2.getZ()-p1.getZ();
  vm[3] = 1.0f;
  formulated = TRUE;
}

float Vector::length() 
{  
  int i;
  double t=0.0f;

  formulate();
  for (i=0;i<3;i++) t += vm[i]*vm[i]; 
  //t=_sqrt_s(t);//fabs(vm[3]); 
  t=sqrt(t);
  return t; 
} 

int Vector::normJudge()
{
  formulate();
  double t=length();

  if (fabs(t)<ZERO){
    //printf("\nERROR:Vector::norm() failed");
    setXYZ(0.0f,0.0f,0.0f);
    return -1; 
  }else{
    vm[0] /= t;
    vm[1] /= t;
    vm[2] /= t;
    return 1;
  }
}


void Vector::norm() 
{ 
  formulate(); 
  double t=length(); 
  //cout <<" $ ";
  if (fabs(t)<ZERO){ 
    cout<<"\nERROR:Vector::norm() failed,t="<<t<<endl;
    cout<<"\n (XYZ)=("<< vm[0]<<","<<vm[1]<<","<<vm[2]<<")"<<endl;
    setXYZ(0.0f,0.0f,0.0f); 
  }else{ 
    vm[0] /= t; 
    vm[1] /= t; 
    vm[2] /= t; 
  } 
} 

void Vector::normInNewSpace(float ssx,float ssy, float ssz)
{//ssx,ssy,ssz is multiplied to change Vector to the current space
 //this is to get a vector in current space equal to a normlized 
 // in the original space
 formulate();
 ssx *=ssx; ssy *= ssy; ssz *= ssz;
 double t = sqrt(vm[0]*vm[0]/ssz+vm[1]*vm[1]/ssy+vm[2]*vm[2]/ssz); 
 deScale(t,t,t);

}


float Vector::angle(const Vector &V) 
{//return value should be in [0,PI], if -1, means error
  double t1 = length();
  Vector tmpV(V);
  double t2 = tmpV.length();
  
  if ((fabs(t1)<ZERO) ||(fabs(t2)<ZERO)) {
	//printf("\n Vector undefined in Vector::angle(Vector)");//---need to be restored
        //if (t1 == 0) outputs();
	return -1.0f; 
  }
  float alpha = dotProduct(tmpV)/t1/t2;
  //alpha = facos (alpha);
  alpha = acos (alpha);
  return alpha;
}
float Vector::angle(const Plane &PL)
{
  Plane tmpPL(PL);
  float t = tmpPL.angle(*this);
  if (fabs(t+1.0f)<ZERO){
	printf("\nERROR:Vector::angle(Plane)");
	return -1.0f;
  } else return t;
}
Vector Vector::projection(const Plane &PL)
{
  Plane tmpPL(PL);
  return tmpPL.projection(*this);
}
Vector Vector::projection(Vector N, Vector PtR)
{
  Plane tmpPL(PtR,N);
  return projection(tmpPL);
}

float Vector::projectedLength(Vector B)
{//compute the length of A on B, return -1 means ERROR
  float len=B.length();
  if (fabs(len-0)<ZERO){
	printf("\nERROR: Vector::projectedLength()");
	return -1;
  }else{
	return (dotProduct(B)/len);
  }
}

Vector Vector::operator+(Vector vt)
{
  Vector tmp; 
  int i;
 
  formulate(); 
  vt.formulate(); 

  for (i=0;i<3;i++) 
     tmp.vm[i] = vm[i] + vt.vm[i]; 
  
  return tmp;
}
Vector Vector::operator-(Vector vt) 
{ 
  Vector tmp; 
  int i;
 
  formulate(); 
  vt.formulate(); 

  for (i=0;i<3;i++) 
     tmp.vm[i] = vm[i] - vt.vm[i]; 
  
  return tmp; 
} 
 
Vector Vector::operator*(float dd) 
{ 
  Vector tmp; 
 
  formulate();
  tmp.vm[0] = dd * vm[0]; 
  tmp.vm[1] = dd * vm[1]; 
  tmp.vm[2] = dd * vm[2]; 
   
  return tmp; 
} 
Vector Vector::negative(void)
{ 
  Vector tmp;
  int i;
  for (i=0;i<3;i++) tmp.vm[i] = -vm[i];
  tmp.vm[3] = vm[3];
  tmp.formulated = formulated;
  return tmp;
}

Vector Vector::operator-(void)
{ 

  return negative();
}
Vector Vector::operator/(float dd) 
{ 
  Vector tmp; 

  formulate(); 
  if (fabs(dd-0.0f)>ZERO) { 
	  tmp.setX( vm[0] / dd ); 
	  tmp.setY( vm[1] / dd ); 
	  tmp.setZ( vm[2] / dd ); 
  }   
  return tmp; 
} 

Vector Vector::operator*(Matrix nM)
{/*
   int row,col;
   double t[4];
   fPoint pp;

   formulate(); 

   for (col=0;col<dm;col++){
	t[col] = 0.0f;
	for (row=0;row<dm-1;row++) 
	   t[col] += vm[row] * nM.mm[row][col];
       // the reason to loop to dm-1 rather than dm is 
       // because this is for vector transformation, 
       //so need to remove the translation of the origin (0,0,0)
  }
  pp.setXYZ((float)t[0],(float)t[1],(float)t[2],1.0f);
  pp.formulate();
  return Vector(pp);
*/
  fPoint p1=fPoint(vm[0],vm[1],vm[2],vm[3])*nM;
  fPoint p2=fPoint(0,0,0)*nM;
  return Vector(p2,p1);
}	 

Vector Vector::rotateX(float ang)
{
   Matrix mRt;

   mRt.loadIdentity();
   mRt.rotateX(ang);
   fPoint pt(vm[0],vm[1],vm[2],vm[3]);

   pt = pt * mRt;

   return Vector(pt);   
}


Vector Vector::rotateY(float ang)
{
   Matrix mRt;

   mRt.loadIdentity();
   mRt.rotateY(ang);
   fPoint pt(vm[0],vm[1],vm[2],vm[3]);

   pt = pt * mRt;

   return Vector(pt);
}

Vector Vector::rotateZ(float ang)
{
   Matrix mRt;

   mRt.loadIdentity();
   mRt.rotateZ(ang);
   fPoint pt(vm[0],vm[1],vm[2],vm[3]);

   pt = pt * mRt;

   return Vector(pt);
}

Vector Vector::vvScale(Vector ad)
{
 formulate(); ad.formulate();
 return Vector(vm[0]*ad.vm[0], vm[1]*ad.vm[1], vm[2]*ad.vm[2]);
}

float Vector::dotProduct(Vector& Vt) 
{ 
	formulate(); 
	Vt.formulate(); 
	float t=0.0f; 
	int i;
  	for (i=0;i<3;i++)
	  t += vm[i] * Vt.vm[i]; 
	return t; 
} 
 
float Vector::innerProduct(Vector& Vt) 
{ 
	return dotProduct(Vt); 
} 
 
float Vector::scalarProduct(Vector& Vt) 
{ 
	return dotProduct(Vt); 
} 


BOOL Vector::isEqual(Vector Vt)
{
  if (sameDirection(Vt) && fabs(length()-Vt.length())<ZERO)
     return TRUE;
  else 
     return FALSE;
}

BOOL Vector::sameDirection(Vector Vt)
{
  float ang = angle(Vt);
    //printf("\nangle = %5.4f",ang); //---
  if (ang == -1.0) return FALSE;
  if ((ang<0.001)||(fabs(PI-ang)<0.001))
     return TRUE;
  else {
     printf("\nAngle of two vectors = %8.7f\n",ang);
     return FALSE;
  }
}

Vector Vector::vectorRightProduct(Vector& Vt) 
{ 
 
	Vector T; 
 
	formulate(); 
	Vt.formulate(); 
	T.vm[0] = vm[1]*Vt.vm[2] - vm[2]*Vt.vm[1]; 
	T.vm[1] = vm[2]*Vt.vm[0] - vm[0]*Vt.vm[2]; 
	T.vm[2] = vm[0]*Vt.vm[1] - vm[1]*Vt.vm[0]; 
	T.vm[3] = 1.0f; 
 
	return T; 
} 

Vector Vector::vectorRightCross(Vector& Vt)
{
  	return vectorRightProduct(Vt);
}

Vector Vector::vectorLeftProduct(Vector& Vt) 
{ 
	return Vt.vectorRightProduct(*this); 
} 

Vector Vector::vectorLeftCross(Vector& Vt)
{
	return Vt.vectorRightProduct(*this);
}

//--------------------------------------------
//class Matrix members
//
//
//------------------------------------------
Matrix::Matrix(void)
{
  reset();
}

Matrix::Matrix(float m11, float m12, float m13,float m14,
       	 float m21, float m22, float m23,float m24,
       	 float m31, float m32, float m33, float m34,
		 float m41, float m42, float m43, float m44) 
{
  setValue( m11,  m12,  m13, m14,
	  m21,  m22,  m23, m24,
       	  m31,  m32,  m33,  m34,
	  m41,  m42,  m43,  m44);
}
Matrix::Matrix(float m[16])
{
  setValue(m);
}
Matrix::Matrix(const Matrix &nM)
{
  int col,row;
  for (row=0;row<dm;row++)
  for (col=0;col<dm;col++)
     mm[row][col] = nM.mm[row][col];
}

Matrix::~Matrix(void)
{
}

void Matrix::setValue(float m11, float m12, float m13,float m14,
       	 float m21, float m22, float m23,float m24,
       	 float m31, float m32, float m33, float m34,
	 float m41, float m42, float m43, float m44)
{
  mm[0][0] = m11; mm[0][1] = m12; mm[0][2] = m13; mm[0][3] = m14;
  mm[1][0] = m21; mm[1][1] = m22; mm[1][2] = m23; mm[1][3] = m24;
  mm[2][0] = m31; mm[2][1] = m32; mm[2][2] = m33; mm[2][3] = m34;
  mm[3][0] = m41; mm[3][1] = m42; mm[3][2] = m43; mm[3][3] = m44; 
}

void Matrix::setValue(float m[16])
{
  int row,col;

  for (row=0;row<4;row++)
  for (col=0;col<4;col++){
    mm[row][col] = m [row*4+col];
  }
}

void Matrix::reset()
{ 
  int i,j;
  for (i=0;i<dm;i++)
    for (j=0;j<dm;j++){
      mm[i][j] = 0.0;
    }
}

void Matrix::loadIdentity()
{
  int i;

  reset();
  for (i=0;i<dm;i++)
    mm[i][i] =  1.0;
}
Matrix Matrix::mulRight(Matrix nM) 
{// this * nM 
  Matrix tmp; 
 
  int row,col,i; 
 
  for (row=0;row<dm;row++) 
    for(col=0;col<dm;col++){ 
      //tmp.mm[row][col]=0.0; 
      for (i=0;i<dm;i++){ 
	tmp.mm[row][col] += mm[row][i]*nM.mm[i][col]; 
      } 
    } 
  return tmp; 
} 
Matrix Matrix::mulRight3(Matrix nM) 
{// this * nM 
  Matrix tmp; 
 
  int row,col,i; 
 
  for (row=0;row<3;row++) 
    for(col=0;col<3;col++){ 
      //tmp.mm[row][col]=0.0; 
      for (i=0;i<3;i++){ 
	tmp.mm[row][col] += mm[row][i]*nM.mm[i][col]; 
      } 
    } 
  tmp.mm[0][3] = tmp.mm[1][3] = tmp.mm[2][3] = 0.0f;
  tmp.mm[3][0] = tmp.mm[3][1] = tmp.mm[3][2] = 0.0f;
  tmp.mm[3][3] = 1.0f;
  return tmp; 
} 
 
Matrix Matrix::mulLeft(Matrix nM) 
{//nM*this 
  Matrix tmp; 
 
  int row,col,i; 
 
  for (row=0;row<dm;row++) 
    for(col=0;col<dm;col++){ 
      //tmp.mm[row][col]=0.0; 
      for (i=0;i<dm;i++){ 
		tmp.mm[row][col] += nM.mm[row][i]*mm[i][col]; 
      } 
    } 
  return tmp; 
} 

Matrix Matrix::operator*(Matrix nM)
{ 
	return mulRight(nM);
}
Matrix Matrix::operator*(float ff)
{ // this[i][j]*ff
  Matrix tmp; 
 
  int row,col; 
 
  for (row=0;row<dm;row++) 
    for(col=0;col<dm;col++){            
		tmp.mm[row][col] = ff*mm[row][col];       
    } 
  return tmp;  
} 

Matrix Matrix::operator+(Matrix nM)
{ // this[i][j]+nM[i][j]
  Matrix tmp; 
 
  int row,col; 
 
  for (row=0;row<dm;row++) 
    for(col=0;col<dm;col++){            
		tmp.mm[row][col] = nM.mm[row][col]+mm[row][col];       
    } 
  return tmp;  
} 

Matrix Matrix::transpose(void)
{
  Matrix tmpM;
  int row,col;

  for(row=0;row<dm;row++)
  for(col=0;col<dm;col++)
    tmpM.mm[row][col] = mm[col][row];

  for(row=0;row<dm;row++)
  for(col=0;col<dm;col++)
    mm[row][col] = tmpM.mm[row][col];

  return tmpM;
}

Matrix  Matrix::translate(float x, float y, float z)
{
  Matrix tmpM;

  tmpM.loadIdentity();
  tmpM.mm[3][0] = -x;
  tmpM.mm[3][1] = -y;
  tmpM.mm[3][2] = -z;

  (*this) = tmpM = (*this) * tmpM;
  return tmpM;
}

Matrix Matrix::translate(fPoint pp)
{
  float XX,YY,ZZ;
  pp.getXYZ(XX,YY,ZZ);
  return translate(XX,YY,ZZ);
}
 
void Matrix::newTranslate(float x, float y, float z)
{ 
  loadIdentity();
  mm[3][0] = -x; 
  mm[3][1] = -y; 
  mm[3][2] = -z; 
}

void Matrix::newRotateX(float angle)
{
  loadIdentity();
  mm[1][1] = mm[2][2] = cosf(angle);
  mm[1][2] = sin(angle);
  mm[2][1] = -sin(angle);
}

Matrix Matrix::rotateX(float angle)
{ 
  Matrix tmpM;

  tmpM.loadIdentity();
  tmpM.mm[1][1] = tmpM.mm[2][2] = cosf(angle);
  tmpM.mm[1][2] = sin(angle);
  tmpM.mm[2][1] = -sin(angle); 
 
  (*this) = tmpM = (*this) * tmpM; 
  return tmpM;
}

Matrix Matrix::rotateY(float angle)
{
  Matrix tmpM; 
 
  tmpM.loadIdentity();
  tmpM.mm[0][0] = tmpM.mm[2][2] = cosf(angle);
  tmpM.mm[2][0] = sin(angle);
  tmpM.mm[0][2] = -sin(angle);
  
  tmpM = (*this) = (*this) * tmpM; 
  return tmpM; 
} 

void Matrix::newRotateY(float angle)
{
  loadIdentity();
  mm[0][0] = mm[2][2] = cosf(angle);
  mm[2][0] = sin(angle);
  mm[0][2] = -sin(angle);
}

Matrix Matrix::rotateZ(float angle)
{
  Matrix tmpM; 
 
  tmpM.loadIdentity();
  tmpM.mm[0][0] = tmpM.mm[1][1] = cosf(angle);
  tmpM.mm[1][0] = -sin(angle);
  tmpM.mm[0][1] = sin(angle);
  
  (*this) = tmpM = (*this) * tmpM; 
  return tmpM; 
} 

void Matrix::setValue(int row, int col, float v)
{
  mm[row][col] = v;
}
float Matrix::getValue(int row, int col)
{
  return mm[row][col];
}

void Matrix::newRotateZ(float angle)
{
  loadIdentity();
  mm[0][0] = mm[1][1] = cosf(angle);
  mm[1][0] = -sin(angle);
  mm[0][1] = sin(angle);
}

Matrix Matrix::scaleXYZ(float sa, float sb, float sc) 
{ 
  Matrix tmpM; 
 
  tmpM.loadIdentity(); 
  tmpM.mm[0][0] = sa; 
  tmpM.mm[1][1] = sb; 
  tmpM.mm[2][2] = sc; 
  
  (*this) = tmpM = (*this) * tmpM; 
  return tmpM; 
 
} 

void Matrix::newScaleXYZ(float sa, float sb, float sc)
{
  loadIdentity();
  mm[0][0] = sa;
  mm[1][1] = sb;
  mm[2][2] = sc;
}
 
Vector Matrix::mulRight(Vector vt) 
{ 
	int row, col; 
	Vector tmpV; 
 
	for (row=0;row<dm;row++){					 
		tmpV.vm[row] = 0.0f; 
		for (col =0; col <dm; col++) 
 			tmpV.vm[row] += mm[row][col] * vt.vm[row]; 
	} 
	return tmpV; 
} 

 
Vector Matrix::operator*(Vector vt) 
{ 
	return mulRight(vt); 
} 

void Matrix::newUVN(Vector U, Vector V, Vector N)
{
  U.norm();
  V.norm();
  N.norm();

  loadIdentity();
  mm[0][0] = U.vm[0];
  mm[1][0] = U.vm[1];
  mm[2][0] = U.vm[2];

  mm[0][1] = V.vm[0];
  mm[1][1] = V.vm[1];
  mm[2][1] = V.vm[2];

  mm[0][2] = N.vm[0];
  mm[1][2] = N.vm[1];
  mm[2][2] = N.vm[2];
}

void Matrix::outputs(void)
{
 int i,j;

 cout << endl; //printf("\n");
 for (i=0; i<dm;i++){
  cout << endl; printf("\n");
  for (j=0;j<dm;j++)
    cout << " "<<mm[i][j]<<" ";
    //printf("(%5.3f)-",mm[i][j]);
  }
  cout <<endl;
  //printf("\n");
}

//-------------------------------------------------
// global functions
//-------------------------------------------------
double getSeconds(void)
{
  struct tms rusage;
  times(&rusage);
  return((double)rusage.tms_utime);
}
