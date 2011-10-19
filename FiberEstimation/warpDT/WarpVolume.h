#ifndef _WARPVOLUME_H_
#define _WARPVOLUME_H_

#include "Basics.h"

//--------------------------------------------------------------
// class WarpVolume: Volume of warping vector field
//
//   warp volume dimensions:  X = 256, Y = 256, Z=?
// Members:
//    dimX,Y,Z: dimensions of the volume, file size would be
//              dimX*dimY*dimZ*(6*4): 6 float components at a voxel
//    path:     the file name to load the DT volume
//    vol :     the memory to hold the volume
//    voxSize:  size of the voxel in byte
//    planSize: number of voxel per slice
//
// Memeber functions:
//    int JacobiGrid(): calculate Jacobi matrix on voxel grid
//		return 1 if good, -1 failed
//    int JacobiRotationMatrix(): get the rotation matrix M at 
//		(xx,yy,zz) which is derived from Jacobi 
//		Matrix at the same location, 
//		return 1: succ; -1: failed
//    fPoint triInterpolationXYZ(fPoint smpt):
//    fPoint triInterpolationXYZ(struct Pt3d ptt):
//		get the warp vector at 'smpt'/'ptt' and
//		return; if undefined, the returned fPoint
//		will be undefined (the 4th component == 0)
//    void   reorganize(struct Pt3d*  &volume, int volZZ):
//		1.switch the X<->Y coordinate value, this is
//		used when the program loads Christos's vector
//		field, since he took horizontal as Y and 
//		vertical as X while we do this in the other 
//		way
//		2.the vector field defins the destination
//		while what we need is the displacement vector,
//		so the V = V - (x,y,z);
//    int    triInterpolationXYZ(struct Pt3d ptt,struct Pt3d& vii):
//              return -1 if result undefined; return 1 if good.
//              vii contains the resulting interpolated value
// 
//--------------------------------------------------------------

class WarpVolume{
    int dimX,dimY,dimZ;
    struct Pt3d *fdMap;
    float xres,yres,zres;
    char fWarp[130];
    int voxSize,planSize;

public:
    WarpVolume(char* fVF, int xx=0, int yy=0, 
	float rx=1.0,float ry=1.0, float rz=1.0 );
   ~WarpVolume();

    void   loadVolume(char* ff, Pt3d * &object, int& sampZ);
    void   reorganize(struct Pt3d*  &volume, int volZZ);

    int    JacobiGrid(int col, int row, int slice, float res[9]);
    int    triInterpolationXYZ(struct Pt3d ptt,struct Pt3d& vii);
    int    triInterpolationXYZ(fPoint smpt,    struct Pt3d& vii);
    float  triInterpolation(float* w, float v0, float v1, float v2,
		    float v3, float v4, float v5, float v6,float v7);
    int    Jacobi(float col, float row, float slice,float res[9]);
    int    JacobiRotationMatrix(float xx, float yy, float zz, Matrix& M);
    int    JacobiRotationMatrix(fPoint pt, Matrix& M);
    int    JacobiRotationMatrix(struct Pt3d pt, Matrix& M);

    struct Pt3d getVectorAt(int col, int row, int slice);

};


#endif //_WARPVOLUME_H_
