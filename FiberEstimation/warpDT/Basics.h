#ifndef _BASICS_H_
#define _BASICS_H_

#include "Global.h"

//#define ZERO 0.00001

//enum {FALSE, TRUE}; //BOOL type
//enum {ABNORMAL, NORMAL}; // tResult type

//typedef int BOOL;
//typedef int tResult;


//---------------------------------------------------------------
//  Class fPoint:
//            a floating point coordinate (x,y,z) in 3-D space	
//     which is expressed in the form of (x*aa, y*aa, z*aa, *aa) 
//     in 4-dim space
//
//
//---------------------------------------------------------------

class fPoint
{
 protected:
    float vm[4];
    int formulated; //control if member function 'formulate' is to be done

 public:
    fPoint();
    fPoint(int a, int b,int c, int d=1);	
    fPoint(float a,float b, float c, float d=1.0f);
    fPoint(struct Pt3d vox);
    fPoint(const fPoint &pt);
    ~fPoint(void);

    
    void setX(int a);
    void setX(float a);
    void setY(int b);
    void setAA(int d);
    void setY(float b);
    void setZ(int c);
    void setZ(float c);
    void setAA(float d);

    void setXYZ(int a,int b,int c, int d=1.0f);
    void setXYZ(float a,float b, float c, float d=1.0F);

    void formulate(void);
    void scale(float sx,float sy, float sz);
    void deScale(float sx,float sy, float sz);


    float getX(void);
    float getY(void);
    float getZ(void);
    void  getXYZ(float& a,float& b, float& c);

    float distance(float a, float b, float c, float d=1.0f);
    float distance(fPoint pt);

    fPoint proportionPoint(fPoint B, float t1, float t2);

    fPoint negative(void);
    fPoint operator-(void);
    fPoint operator*(float ss);
    fPoint operator*(Matrix nM);
    fPoint operator+(fPoint nP);
    fPoint operator-(fPoint nP);
    BOOL   equal(fPoint p);
    BOOL   isUndefined();


    //friend Plane;
    //friend Line;
    //friend Image;

    friend class Plane;
    friend class Line;
    friend class Image;
//---------------------------------
    virtual void outputs(void);
};

//---------------------------------------------------------------
//  Class Vector: 
//            a floating point Vector (x,y,z) in 3-D space which 
//    is expressed in the form of (x*aa, y*aa, z*aa, *aa) in 4-d 
//    space 
// 
// 
//--------------------------------------------------------------- 

class Vector :public fPoint
{
 public:
  Vector();
  Vector(fPoint p);
  Vector(fPoint p1,fPoint p2);
  Vector(int x,int y, int z);
  Vector(float x, float y, float z);
  Vector(const Vector &V);
  ~Vector();

  void  norm();
  int   normJudge();
  void normInNewSpace(float ssx,float ssy, float ssz);
  float length(); 
  float angle(const Vector &V);
  //return value should be in [0,180], if -1, means error
  float angle(const Plane &PL);//return -1 if ERROR
  Vector projection(const Plane &PL);
  Vector projection(Vector N, Vector PtR=Vector(fPoint(0,0,0)));
  float projectedLength(Vector B);//get length of projecting *this on B, 
				  //return -1 if ERROR
  void setValue(fPoint p1, fPoint p2);
  void setValue(fPoint p);
 
  float dotProduct(Vector& Vt); 
  Vector vvScale(Vector ad);
  float innerProduct(Vector& Vt); 
  float scalarProduct(Vector& Vt); 
  Vector vectorRightProduct(Vector& Vt); 
  Vector vectorLeftProduct(Vector& Vt);
  Vector vectorRightCross(Vector& Vt);
  Vector vectorLeftCross(Vector& Vt);
 
  BOOL   sameDirection(Vector Vt);
  BOOL   isEqual(Vector Vt);

  Vector negative(void);

  Vector operator+(Vector vt);
  Vector operator-(Vector vt);
  Vector operator-(void);

  Vector operator*(float dd);
  Vector operator/(float dd);
  Vector operator*(Matrix nM);

  Vector rotateX(float ang);
  Vector rotateY(float ang);
  Vector rotateZ(float ang);

  friend class Matrix; 
  friend class Plane;
  friend class Line;
  friend class Image;
};



//---------------------------------------------------------------
//    Class Matrix:
//            define a 4-order matrix
//
//
//---------------------------------------------------------------


class Matrix
{
 protected:
#define dm 4
   float mm[dm][dm]; 
 
 public:
   Matrix(void);
   Matrix(float m[16]);
   Matrix(float m11, float m12, float m13,float m14,
       	 float m21, float m22, float m23,float m24,
       	 float m31, float m32, float m33,float m34,
       	 float m41, float m42, float m43,float m44);
   Matrix(const Matrix &nM);
   ~Matrix(void);

   void   reset();
   void   setValue(float m[16]);
   void   setValue(int row, int col, float v);
   void setValue(float m11, float m12, float m13,float m14,
       	 float m21, float m22, float m23,float m24,
       	 float m31, float m32, float m33, float m34,
	 float m41, float m42, float m43, float m44);
   float  getValue(int row, int col);
   void   loadIdentity();
   Matrix transpose(void);

   Matrix rotateX(float angle);
   Matrix rotateY(float angle);
   Matrix rotateZ(float angle);
   Matrix translate(float x, float y, float z);
   Matrix translate(fPoint pp);
   Matrix scaleXYZ(float sa, float sb, float sc); 

   void newScaleXYZ(float sa, float sb, float sc);
   void newRotateX(float angle);
   void newRotateY(float angle);   
   void newRotateZ(float angle);
   void newTranslate(float x, float y, float z);

   void newUVN(Vector U, Vector V, Vector N);

   Matrix operator*(float ff);
   Matrix operator*(Matrix nM); 
   Matrix mulRight(Matrix nM); 
   Matrix mulLeft(Matrix nM); 
   Matrix mulRight3(Matrix nM);

   Vector operator*(Vector vt); 
   Vector mulRight(Vector vt); 
 
   Matrix operator+(Matrix nM);
   void outputs(void); 
 
   //friend Vector;
   //friend fPoint;

   friend class Vector;
   friend class fPoint;
};

#endif
