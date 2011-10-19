#include <stdio.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>

#include "GPrimitive.h"


//-------------------------------------------
//class Plane member functions
//
//
//-------------------------------------------

Plane::Plane(float AA, float BB, float CC, float DD)//Ax+By+Cz+D=0 
{
 	 A = AA; B = BB; C = CC; D = DD;

}

Plane::Plane(Vector PtR, Vector N) //(R-PtR).N=0
{
    setPara(PtR, N);
}

Plane::Plane(fPoint PtR, Vector N)
{
    setPara(PtR, N);
}

Plane::Plane(const Plane &P)
{
   A = P.A; B = P.B; C = P.C; D = P.D;
}

Plane::~Plane()
{

}
void Plane::outputs(void)
{
  printf("\n The plane is:");
  printf(" %5.4f X + %5.4f Y + %5.4f Z + %5.4f = 0\n", A,B,C,D);
}

void Plane::setPara(Vector PtR, Vector N)
{
   A = N.getX(); 
   B = N.getY(); 
   C = N.getZ(); 
   D = - PtR.dotProduct(N);
}
void Plane::setPara(fPoint PtR, Vector N)
{
   Vector V(PtR);
   A = N.getX(); B = N.getY(); C = N.getZ();
   D = - V.dotProduct(N);

}

Vector Plane::getNVector(void)
{
	Vector Vt;

	Vt.setXYZ(A,B,C);
	return Vt;
}

void Plane::getABCD(float& AA, float& BB, float& CC, float& DD)
{
  AA =  A; BB = B; CC = C; DD = D;
}

float Plane::distance(fPoint Pt)
{
  float t,x,y,z;

  Pt.getXYZ(x,y,z);

  //t= fabs(A*x+B*y+C*z+D)/_sqrt_s(A*A+B*B+C*C);
  t= fabs(A*x+B*y+C*z+D)/sqrt(A*A+B*B+C*C);

  return t;
}
tResult Plane::intersection(Line L, fPoint& Pt)
{
	float t;
	return intersection(L, Pt, t);
}

tResult Plane::intersection(Line L, fPoint& Pt, float &t)
{//When tResult == NORMAL: Pt contains the intersection point, t is the parameter for L to get Pt
 //     tResult == ABNORMAL: Pt and t are undefined (at this time, L and this Plane are parallel)

	fPoint R0;
	Vector V;
	int i;
	
	L.getPara(R0,V);
	R0.formulate(); V.formulate();
	float sv=A*V.vm[0]+B*V.vm[1]+C*V.vm[2];
	if (fabs(sv-0)<ZERO){
		//printf("\nERROR: Plane::intersection Abnormal");
		//printf("\n  Line: "); L.outputs();
		//printf("\n  Plane:"); outputs();
		//exit(0);
		return ABNORMAL;
	}
	else{
		sv = -(A*R0.vm[0]+B*R0.vm[1]+C*R0.vm[2]+D)/sv;
		for (i=0;i<3;i++)
		  Pt.vm[i] = R0.vm[i] + V.vm[i] * sv;
		t = sv;
		return NORMAL;
	}
}

float Plane::angle(const Plane &PL)
{//the return value is in [0,Pi/2], otherwise return -1 meaning ERROR
  Plane tmpPL(PL);
  Vector v1=getNVector();
  Vector v2=tmpPL.getNVector();
  float t1 = v1.length(), t2 = v2.length();

  if ((fabs(t1-0)<ZERO) ||(fabs(t2-0)<ZERO)){
	printf("\n ERROR: Plane::angle(Plane):Plane undefined");
	return -1.0f;
  }
  float alpha = v1.dotProduct(v2)/t1/t2;
  //printf("\nalpha = %6.5f",alpha);
  //alpha = facos (fabs(alpha));
  alpha = acos (fabs(alpha));
  //printf("\nalpha 1 = %6.5f",alpha);
  return alpha;
}

float Plane::angle(const Line &L)
{// the return value should be between 0 and 90, if -1, means error.
  Vector v1 = getNVector();
  Vector v2;
  fPoint p;
  Line LL(L);

  LL.getPara(p,v2);
  float t1 = v1.length(), t2 = v2.length();

  if ((fabs(t1-0)<ZERO) ||(fabs(t2-0)<ZERO)){
	printf("\n Error: Plane::angle(Line)");
	return -1.0f;
  }

  float alpha = v1.dotProduct(v2)/t1/t2;
  //alpha = facos (fabs(alpha));
  alpha = acos (fabs(alpha));
  return alpha;
}
float Plane::angle(const Vector &vect)
{//return value should be in [0,180], otherwise ERROR
  Vector v1 = getNVector(),v2(vect);
  float t1 = v1.length(), t2 = v2.length();

  if ((fabs(t1-0)<ZERO) ||(fabs(t2-0)<ZERO)) {
     printf("\nPlane undefined, ERROR");
     return -1.0f;
  }
  float alpha = v1.dotProduct(v2)/t1/t2;
  //alpha = facos (alpha);
  alpha = acos (alpha);
  return alpha;
}

Vector Plane::projection(const Vector &BB)
{
   Vector temp, A1 = getNVector(),B1(BB);
   float t = A1.length();

   if (fabs(t)<ZERO) {
	printf("\n the plane is undefined, ERROR");
        exit(0);
   }  
   t = t*t; 
   temp = B1 - A1*(A1.innerProduct(B1)/t);

   //A1.innerProduct(B1)/|A1| is the length of B1's projection on A1
   //A1*(A1.innerProduct(B1)/t) is the vector B1 projected on A1
   // B1 - A1*(A1.innerProduct(B1)/t) is the vector B1 on the plane
   return temp; 
}
//-----------------------------------------------------
//  class Line member functions
//
//
//
//-----------------------------------------------------
Line::Line(void)
{
	R0.setXYZ(0,0,0);
	V.setXYZ(0,0,0);
}

Line::Line(fPoint Pt, Vector Vt)
{
	R0 = Pt;
	V  = Vt;
        R0.formulate(); 
        V.formulate();
}
Line::Line(Plane P1, Plane P2)
{
	float p,q,r;
	float A1,B1,C1,D1,A2,B2,C2,D2;
	P1.getABCD(A1,B1,C1,D1);
	P2.getABCD(A2,B2,C2,D2);

	p=B1*C2-B2*C1;
	q=C1*A2-A1*C2;
	r=A1*B2-B1*A2;

	R0.setXYZ(p,q,r);

//under construction
}

Line::Line(fPoint Pt1, fPoint Pt2)
{
	Vector tVect(Pt1,Pt2);
	V = tVect;
	R0 = Pt1;
        R0.formulate(); V.formulate();
}
Line::Line(const Line &LL)
{
  V = LL.V;
  R0 = LL.R0;
  R0.formulate(); V.formulate();
}

Line::~Line(void)
{

}

fPoint Line::getPoint(float t)
{
  float Xx,Yy,Zz;
  getPoint(t, Xx, Yy, Zz);
  return fPoint(Xx,Yy,Zz);
}

void Line::getPoint(float t, float &X, float &Y, float &Z)
{
     X = R0.vm[0] + t * V.vm[0];
     Y = R0.vm[1] + t * V.vm[1];
     Z = R0.vm[2] + t * V.vm[2];
}
Vector Line::getDirection(void)
{
  return V;
}
void Line::getPara(fPoint& Pt, Vector& Vt)
{
	Pt=R0; Vt=V;
}

tResult Line::intersection(Plane P, fPoint& Pt)
{//When tResult == NORMAL, Pt contains the intersection point
 //     tResult == ABNORMAL, Pt is undefined
	float A,B,C,D;
	int i;

	P.getABCD(A,B,C,D);
	R0.formulate();V.formulate();
        float sv=A*V.vm[0]+B*V.vm[1]+C*V.vm[2];

        if (fabs(sv-0)<ZERO) return ABNORMAL;
        else{
                sv = (A*R0.vm[0]+B*R0.vm[1]+C*R0.vm[2]+D)/sv;
                for (i=0;i<3;i++)
		  Pt.vm[i] = R0.vm[i] -V.vm[i] * sv; 
                return NORMAL;
        }
}
float Line::distance(fPoint P)
{
  Vector R0P(R0,P);
  double Tv,Tdot,T;
  Tv = V.vm[0]*V.vm[0]+V.vm[1]*V.vm[1]+V.vm[2]*V.vm[2];
  if (fabs(Tv-0)<ZERO){
	printf("\nLine::distance(): the line meaningless:V = 0");
	return -1.0f;
  }else{
	Tdot = V.dotProduct(R0P);	
        //Tr0p = R0P.length();
	//T = sqrt(Tv*Tv*Tr0p*Tr0p-Tdot*Tdot)/Tv;
 	T = sqrt(((R0P.vm[0]*R0P.vm[0]+R0P.vm[1]*R0P.vm[1]+R0P.vm[2]*R0P.vm[2])
		*Tv-Tdot*Tdot) / Tv);       
	return ((float)T);
  }
}
float Line::angle(const Plane &PL)
{//return value should be in [0,90], otherwise -1 means ERROR
  Plane tmpPL(PL);
  float t = tmpPL.angle(*this);
  if(fabs(t+1.0f)<ZERO){
	printf("\nERROR: Line::angle(Plane)");
	return -1.0f;
  } else return t;
}

void Line::norm(void)
{
   V.norm();
}

void Line::normInNewSpace(float sx, float sy, float sz)
{// sx,sy,sz to be multiplied to change from old space
 // to the current one
 //this func is to get a vector equal to the normailized 
 //vector in the old space
  V.normInNewSpace(sx,sy,sz);

}


void Line::outputs(void)
{
  printf("\n Line:  X = %5.4f + %5.4f * t",R0.vm[0], V.vm[0]);
  printf("\n        Y = %5.4f + %5.4f * t",R0.vm[1], V.vm[1]);
  printf("\n        Y = %5.4f + %5.4f * t\n",R0.vm[2], V.vm[2]);
}


//-----------------------------------------------------------
//  class Coordinates
//    members: 
//	 On: 	new origin point coordinates in old system
//       Xn,
//	 Yn,
//	 Zn: 	new axis XYZ in terms of the vectors in old
//		   coordinate system
//       mO2N:  transformation matrix from old coordinate to
//		new
//       mN2O:  transformation matrix from new to old
//		mO2N * mN2O = Identity Matrix
//    --------------
//    if Po is a Vector or Point in old coordinate system,
//	 then Po * mO2N will result in a Vector or Point in 
//	 the new coordinate system
//    and vice versa
//-----------------------------------------------------------
Coordinates::Coordinates(fPoint Onew,
	Vector Xnew, Vector Ynew, Vector Znew)
{
   On = Onew;   
   Xn = Xnew;   Yn = Ynew;   Zn = Znew;
   update();
}

Coordinates::~Coordinates()
{
}

void Coordinates::selfTest()
{
  Vector vt = Xn.vectorRightProduct(Yn);
  if (vt.sameDirection(Zn)==FALSE){
    printf("\n The following 3 vectors:");
    Xn.outputs();Yn.outputs();Zn.outputs();
    printf("\n Xn * Yn = :"); vt.outputs();    
    printf("\n New Coordinates can not form a XYZ system\n");
    exit(0);
  }
}
void Coordinates::update(void)
{
   Matrix M1, M2;
   
   selfTest();
   
   M1.loadIdentity(); M1.translate(On);
   M2.newUVN(Xn,Yn,Zn);
   mO2N = M1 * M2;//Old->New
  
   Vector XX(1,0,0),YY(0,1,0),ZZ(0,0,1);
   fPoint o(0,0,0);
   o = o * mO2N; 
   M1.loadIdentity(); M1.translate(o);
   XX = XX * mO2N; YY = YY * mO2N; ZZ = ZZ * mO2N;
     /*
     printf("\n------On------\n");o.outputs();
     printf("\n------XXn------\n");XX.outputs();
     printf("\n------YYn------\n");YY.outputs();
     printf("\n------ZZn------\n");ZZ.outputs();
     */
   M2.newUVN(XX,YY,ZZ);
   mN2O = M1 * M2;//New -> Old
}

Matrix Coordinates::getMatrixO2N()
{
  return mO2N;
}
Matrix Coordinates::getMatrixN2O()
{
  return mN2O;
}

void Coordinates::reset(fPoint Onew,
	Vector Xnew, Vector Ynew, Vector Znew)
{
   On = Onew;   
   Xn = Xnew;   Yn = Ynew;   Zn = Znew;
   update();
}

void Coordinates::getOXYZ(fPoint& Onew,
	Vector& Xnew, Vector& Ynew, Vector& Znew)
{
  Onew = On;
  Xnew = Xn;
  Ynew = Yn;
  Znew = Zn;
}

void Coordinates::toNewCoor(fPoint& pt)
{//pt in Old system, resulting in pt in New system
   pt = pt * mO2N;   
}

void Coordinates::toOldCoor(fPoint& pt)
{//pt in New system, resulting in pt in Old system
   pt = pt * mN2O;   
}

void Coordinates::toNewCoor(Vector& vt)
{//vt in Old system, resulting in vt in New system
   vt = vt * mO2N;   
}

void Coordinates::toOldCoor(Vector& vt)
{//vt in New system, resulting in vt in Old system
   vt = vt * mN2O;   
}

void Coordinates::rotateAroundXn(fPoint& pt, float ang)
//pt (pt in Old Coordinates) rotates ang around Xn
//ang is [0, 360)
{
   Matrix M;
   toNewCoor(pt); //get vt in New coor-system   

   ang = PI*ang/180.0;
   M.loadIdentity();
   M.rotateX(ang);
   pt = pt * M;
   toOldCoor(pt); //cast pt to Old coor-system 
}
void Coordinates::rotateAroundYn(fPoint& pt, float ang)
//pt (pt in Old Coordinates) rotates ang around Yn
//ang is [0, 360)
{
   Matrix M;
   toNewCoor(pt); //get vt in New coor-system   

   ang = PI*ang/180.0;
   M.loadIdentity();
   M.rotateY(ang);
   pt = pt * M;
   toOldCoor(pt); //cast pt to Old coor-system 
}
void Coordinates::rotateAroundZn(fPoint& pt, float ang)
//pt (pt in Old Coordinates) rotates ang around Zn
//ang is [0, 360)
{
   Matrix M;
   toNewCoor(pt); //get vt in New coor-system   

   ang = PI*ang/180.0;
   M.loadIdentity();
   M.rotateZ(ang);
   pt = pt * M;
   toOldCoor(pt); //cast pt to Old coor-system 
}
void Coordinates::rotateAroundXn(Vector& vt, float ang)
//vt (vt in Old Coordinates) rotates ang around Zn
//ang is [0, 360)
{
   Matrix M;
   toNewCoor(vt); //get vt in New coor-system   

   ang = PI*ang/180.0;
   M.loadIdentity();
   M.rotateX(ang);
   vt = vt * M;
   toOldCoor(vt); //cast vt to Old coor-system 
}

void Coordinates::rotateAroundYn(Vector& vt, float ang)
//vt (vt in Old Coordinates) rotates ang around Zn
//ang is [0, 360)
{
   Matrix M;
   toNewCoor(vt); //get vt in New coor-system   

   ang = PI*ang/180.0;
   M.loadIdentity();
   M.rotateY(ang);
   vt = vt * M;
   toOldCoor(vt); //cast vt to Old coor-system 
}

void Coordinates::rotateAroundZn(Vector& vt, float ang)
//vt (vt in Old Coordinates) rotates ang around Zn
//ang is [0, 360)
{
   Matrix M;
   toNewCoor(vt); //get vt in New coor-system   

   ang = PI*ang/180.0;
   M.loadIdentity();
   M.rotateZ(ang);
   vt = vt * M;
   toOldCoor(vt); //cast vt to Old coor-system 
}
