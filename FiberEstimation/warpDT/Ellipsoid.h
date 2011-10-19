#ifndef _ELLIPSOID_H_
#define _ELLIPSOID_H_

#include "Global.h"
#include "Basics.h"
//===========================================================================
//  class Ellipsoid
//
//         X*X         Y*Y          Z*Z
//      --------- + ---------- + ---------- = 1
//       aSqr[0]      aSqr[1]     aSqr[2]
//
//     This class memorizes the center, the 3 axis of the ellipsoid
//     and the lengthes along the axis. A matrix O2N is then computed
//     which maps a coordinate from the original space to the ellipsoid's,
//     so that according to the above equation, it is easy to judge if
//     the point is inside of this ellipsoid
//     
//  Members:
//    a[3]    :  the 3 values along the 3 axises
//    aSqr[3] :  the squre of a[3], to easy the computation
//    V[3]    :  the 3 othorgnal axis vector
//    center  :  the center of the ellipsoid
//    O2N     :  the transformation matrix from original space
//               to the ellipsoid world
//
//  Member functions:
//    Ellipsoid(float aa, float bb, float cc, struct Pt3d orig,
//		Vector ax1, Vector ax2, Vector ax3);
//        aa,bb,cc -> a[3]; orig -> center ; ax1,ax2,ax3 ->V[3]     
//
//    int inside(struct Pt3d P);
//    int inside(float x, float y, float z);
//        judge if point (x,y,z) is inside this ellipsoid
//        return 1 if inside of this ellipsoid (including boundary)
//        return 0 if not
//
//
//==========================================================================

class Ellipsoid{
   float a[3],aSqr[3]; //aSqr[i]=a[i]*a[i]
   struct Pt3d center;
   Vector V[3];
   Matrix O2N; // to transform to the ellipsoid coordinates

public:
   Ellipsoid();
   Ellipsoid(float aa, float bb, float cc, struct Pt3d orig,
		Vector ax1, Vector ax2, Vector ax3);
   Ellipsoid(float aa, float bb, float cc, 
		float origX, float origY, float origZ,
		Vector ax1, Vector ax2, Vector ax3);

   ~Ellipsoid();
   void setParameters(float aa, float bb, float cc, 
		float origX, float origY, float origZ,
		Vector ax1, Vector ax2, Vector ax3);
   int inside(fPoint P);
   int inside(struct Pt3d P);
   int inside(float x, float y, float z);
    //judge if point (x,y,z) is inside this ellipsoid
    // return 1 if in
    // return 0 if not

   //generate a tensor (matrix) of this ellipsoid
   void getTensor(Matrix & aTensor); 
   void getTensor(float *tens); //float tens[6]

   //void getPrimaryDir(float& a, float& b, float& c);
   Vector getPrimaryDir();

   void outputs();
};



#endif // _ELLIPSOID_H_




