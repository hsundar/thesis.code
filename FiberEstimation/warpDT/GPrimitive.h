#ifndef _GPRIMITIVE_H_
#define _GPRIMITIVE_H_

#include "Global.h"
#include "Basics.h"

class Plane
{
	float A,B,C,D;

public:
	Plane(float AA=0.0f, float BB=0.0f, float CC=0.0f, 
		float DD=0.0f);//Ax+By+Cz+D=0
	Plane(Vector PtR, Vector N); //(R-PtR).N=0
	Plane(fPoint PtR, Vector N);
	Plane(const Plane &P);
	~Plane();

        void setPara(Vector PtR, Vector N);
        void setPara(fPoint PtR, Vector N);

        void getABCD(float& AA, float& BB, float& CC, float& DD);
	Vector getNVector(void);	
	float distance(fPoint Pt);
	tResult intersection(Line L, fPoint& Pt);
	tResult intersection(Line L, fPoint& Pt, float &t);
	float angle(const Plane &PL);
	   // the return value should be between 0 and 90, if -1, means error.
	float angle(const Line  &L);
	   // the return value should be between 0 and 90, if -1, means error.
	float angle(const Vector &vect);
	   //return value should be in [0,180]otherwise ERROR with -1
	Vector projection(const Vector &BB);

//---------------------------------------------------
	void outputs(void);

};

class Line
{// Line parallel to  V(a,b,c); passing dot R0(x0,y0,z0)
 // R-R0=tV: x=x0+at, y=y0+bt; z=z0+ct
	fPoint R0;
	Vector V;
public:
	Line();
	Line(fPoint Pt, Vector Vt);
	Line(Plane P1, Plane P2);
	Line(fPoint Pt1, fPoint Pt2);//R0 = Pt1, V = Pt2 - Pt1;
	Line(const Line &LL);
	~Line();

	void 	getPara(fPoint& Pt, Vector& Vt);
        Vector  getDirection(void);
        fPoint 	getPoint(float t);
	void    getPoint(float t, float &X, float &Y, float &Z);

	tResult intersection(Plane P, fPoint& Pt);
	float	angle(const Plane &PL);//return -1 means ERROR
	float	distance(fPoint P);
	void	norm(void);
	void    normInNewSpace(float sx, float sy, float sz);

//------------------------------
	void	outputs();
};

//-----------------------------------------------------------
//  class Coordinates
//    members: 
//	 On:   	new origin point in terms of the coordinates
//		   in the old coordinate system
//       Xn,
//	 Yn,
//	 Zn: 	new axis XYZ in terms of the vectors in old
//		   coordinate system
//       mO2N:  transformation matrix from old coordinate to
//		new
//       mN2O:  transformation matrix from new to old
//		mO2N * mN2O = Identity Matrix
//
//    --------------
//    if Po is a Vector or Point in old coordinate system,
//	 then Po * mO2N will result in a Vector or Point in 
//	 the new coordinate system
//    and vice versa
//
//-----------------------------------------------------------

class Coordinates
{
  	Matrix 	mO2N, mN2O;
	Vector 	Xn,Yn,Zn;
	fPoint 	On;

	void 	update(void);
	void 	selfTest(); 
	//test if the XnYnZn legally form an orthonormal coordinate system

 public:
	Coordinates(fPoint Onew = fPoint(0,0,0),
		    Vector Xnew = Vector(1,0,0),
		    Vector Ynew = Vector(0,1,0),
		    Vector Znew = Vector(0,0,1));
	~Coordinates();

	//set new parameters to update transformation matrix
	void reset(fPoint Onew = fPoint(0,0,0),
		    Vector Xnew = Vector(1,0,0),
		    Vector Ynew = Vector(0,1,0),
		    Vector Znew = Vector(0,0,1));

	Matrix 	getMatrixO2N();//get matrix for Old->New
	Matrix 	getMatrixN2O();//get matrix for New->Old

	void 	getOXYZ(fPoint& Onew,
		  Vector& Xnew, Vector& Ynew, Vector& Znew);//get paramenters
	void 	toNewCoor(fPoint& pt);
	void 	toOldCoor(fPoint& pt);

	void 	toNewCoor(Vector& vt);//vect in old coor result in new coor
	void 	toOldCoor(Vector& vt);//vect in new coor result in old coor
	
	//point pt rotates around vector Xn/Yn/Zn
        //all coordinates in Old system
	void 	rotateAroundXn(fPoint& vt, float ang); //ang is angle [0,360)
	void 	rotateAroundYn(fPoint& vt, float ang); //ang is angle [0,360)
	void 	rotateAroundZn(fPoint& vt, float ang); //ang is angle [0,360)

	//vector vt rotates around vector Xn/Yn/Zn
        //all coordinates in Old system
	void 	rotateAroundXn(Vector& vt, float ang); //ang is angle [0,360)
	void 	rotateAroundYn(Vector& vt, float ang); //not real value [0,PI)
	void 	rotateAroundZn(Vector& vt, float ang);
        //vt in Old coor-system rotate angle 'ang' around vector Zn
};


#endif  // _GPRIMITIVE_H_
