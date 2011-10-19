#ifndef _VOXMAP_H_
#define _VOXMAP_H_

#include "Global.h"
#include "GPrimitive.h"
#include "ImgVol.h"

//--------------------------------------------------------------------
//  class VoxMap:
//        get a map from the original 56-60-Z-slice space through
//        (1) double slice (2) rotate (-90,0,0) (3) slightly ac-pc
//	  adjustment (4) empty slices cropped
//        The only missing step is the final one: from the cropped
//	  volume to the 57-slice Atlas slice
//
//   members:
//        dimX,Y,Z : volume dimemsions
//        planSize : voxel number on a volume plane 
//        voxSize  : number of byte at each voxel
//	  	    ((struct Pt3d)=12B)
//        fdMap    : vox pointer, will contain the forword map
//	  fCropInfo:  file naem that contains cropping information
//		      format:   "%d %d %d" (use program:cropXDR)
//			(totalZ, 1st_nonZero_Z, last_nonZero_Z)
//	  fWarp    :  file name for a vector field warping to Atlas
//
//   member functions:
//        initVox():  initialize the mapping volume to the grid 
//		      coordinates
// 	  triInterpolationXYZ(): tri-lineal interpolation in 
//		      the given vector field
//	  loadVolume(): load a vector field volume from a given file 
//		      to the given memory space
//	  getForwardMap(): create the forword map which maps the original
//		      space to the 57 atlas space
//
//   Mapping generation procedure:
//        1. double slicing (56 slices -> 112 slices)
//	  2. get isotropic volume, and rotate (-90,0,0), then align
//	     to ac-pc (using the given paratmeters when calling 
//	     getForwardMap()), then back to the resolution of 
//	     0.9766*0.9766*1.5 at X,Y,Z directions
//
//    void averageAndSave(int nSlc, char * toStore):
// 	  average the volume at a gap of every nSlc slices
// 	  result in a volume of (dimX,dimY,nSlc)volume 
// 	  and store to file "toStore"
//    void reorgainze(int XYswitch, int ADD_REMOVE_GRID) :
//	  control to adjust the value of the vector field:
//	    XYswitch: switch X and Y coordinates in the vector field
//	    ADD_REMOVE_GRID==1 :  add
//		change the vector field from relative to absolute coordinates
//	    ADD_REMOVE_GRID== -1 :  Remove
//		change the vector field from absolute coordinates to relative
//   void  colormapRGB3Volume():
//         convert vector field components XYZ into 3 byte-volumes
//   struct Pt3d getVectorAt(int col, int row, int slice);
//             returned value Undefined when (col,row,slice) goes beyond volume
//
//   int    getVectorAt(int col, int row, int slice, struct Pt3d & P);
//             when (col,row,slice) in volume, return 1, and result in P
//             when (col,row,slice) out of volume, return 0, and P==0
//   int    getVectorAtAnyPosition(float px, float py, float pz, struct Pt3d & P);
//             same as getVectorAt(int col, int row, int slice, struct Pt3d & P)
//	       returned value based on tri-linear interpolation
//
//--------------------------------------------------------------------

class VoxMap {
     int dimX, dimY, dimZ,voxSize,planSize;
     float rX, rY, rZ; //resolution along X, Y, Z;
     struct Pt3d *fdMap;
     char fCropInfo[200],fWarp[200];
    //-- for judging point in polyhedron
     struct Pt3d ver[8], wXYZ, Box[12][2]; //Box: surrending each face
     int    faces[12][3];
     Plane  PL[12];
     float  radius;
     float  bmin[3], bmax[3]; //bounding box of the whole polyhedron
     int  dltX,dltY,dltZ;

public:
     VoxMap(); //just for debugging testPtinPolyhedron() 
     VoxMap(char* fCrop, char* fVF, int LOADEXISTING=NO,
	   int xx=0, int yy=0, int zz=0, 
	   float rx=1.0,float ry=1.0, float rz=1.0);
    ~VoxMap();

     void   initVox();
     fPoint triInterpolationXYZ(struct Pt3d* &volume,
	     const int sampZ, struct Pt3d vox);
     fPoint triInterpolationXYZ(struct Pt3d* &volume, const int sampZ, fPoint smpt);
     void   loadVolume(char* ff, Pt3d * &object, int& sampZ) ;
     void   switchXY(struct Pt3d*  &volume, int volZZ);
     void   loadCropInfo(int& infoZ, int& infoCrop1, int& infoCrop2);
     struct Pt3d getMapAt(int x1, int y1, int z1);
     void   setMapAt(int x1, int y1, int z1, struct Pt3d& vv);
     void   setMapAt(int x1, int y1, int z1, float a, float b, float c);
     void   setBlendMapAt(int x1,int y1, int z1, float& a, float& b, float& c);
     void   setAvgBlendMapAt(int x1,int y1, int z1, float a, float b, float c);

     void   createForwardMap(float rx, float ry, float rz);
     void   createSusumuMapStep1(float rx, float ry, float rz);
     struct Pt3d* &getMap();

     void   saveMapVolume(char* name);
     void   saveVecVolume(struct Pt3d * &volume, int volZ, char * where);
     int    getVolumeZ(char* ff);
     struct Pt3d getVectorAt(int col, int row, int slice);
     int    getVectorAt(int col, int row, int slice, struct Pt3d & P);
     int    getVectorAtAnyPosition(float px, float py, float pz, struct Pt3d & P);

     void   resetVFvolume(struct Pt3d* &volume, int warpZ);
     int    upper(float a); // get the upper bound integer value: upper(15.5)=16
     int    getGridCoverd(int XX,int YY, int ZZ, 
                 struct Pt3d *gridPts, float wei[][8], int limit);
     float  triInterpolation(float* w, float v0, float v1, float v2, 
		float v3, float v4, float v5, float v6,float v7);

     void   getXYZdim(int& a, int& b, int& c);

     struct Pt3d  addVec(struct Pt3d a, struct Pt3d b );
     int    inBox(struct Pt3d  q);
     void   getBox4Eachface();
     char   inPolyhedron( float x1, float y1, float z1);
     char   boxTest ( int n, struct Pt3d a, struct Pt3d b );
     void   randomRayEndpoint( struct Pt3d & endPt, float radius);
     char   segTriInt(int T[3],struct Pt3d q,struct Pt3d r, struct Pt3d & p );
     char   segPlaneInt( int  T[3], struct Pt3d q, struct Pt3d r, 
		struct Pt3d& p, int& m);
     char   inTri3D( int T[3], int m, struct Pt3d p );
     int    areaSign( float a[3], float b[3], float c[3]); 
     char   inTri2D( float Tp[3][3], float pp[3] );
     int    VolumeSign( struct Pt3d a, struct Pt3d b, 
	       struct Pt3d c, struct Pt3d d );
     char   segTriCross( int T[3], struct Pt3d q, struct Pt3d  r );
     int    planeCoeff( int T[3], struct Pt3d &N, double& D );
     void   NormalVec(struct Pt3d a, struct Pt3d b, struct Pt3d c,
		 struct Pt3d& N );
     double dot( struct Pt3d a, struct Pt3d b );
     struct Pt3d subVec( struct Pt3d a, struct Pt3d b);
     void   getWeight(int x1, int y1, int z1, float *ww);
     double distance(float a, float b, float c, struct Pt3d pp);
     double distance(struct Pt3d pp1, struct  Pt3d pp2);

     void  colormapRGB3Volume();//convert vector field components XYZ into 3 byte-volumes

     //---test purpose---
     int  testPtinPolyhedron();
     void outPt(char* ss,struct Pt3d ppp);

     void reorganize(int XYswitch, int ADD_REMOVE_GRID);

    //----for Christos's funding application
    void ChristosExample();
    void ChristosExample2();
    void averageAndSave(int nSlc, char * toStore);
    //average the volume at a gap of every nSlc slices
    // result in a volume of (dimX,dimY,nSlc)volume 
    // and store to file "toStore"

    void VF3DFiber(int Z1, float& x1, float& y1, float& z1, 
        float& deltx, float& delty, float& deltz, 
	float fx, float fy, float fz, ImgVol& NV, FILE* ff=NULL);
	//VN: an accompanying intensity image for debuging view
	//ff: an accompanying file recording all fiber vectors
    void VF3Dsimulation(FILE * binDotFile=NULL);
	// for my DT Procrustean warping experiment
	//binDotFile: an accompanying file recording all fiber vectors
    struct Pt3d getUnitVec(float x1, float y1, float z1,
	                float x2, float y2, float z2);

    void recoverFiber(struct Pt3d * dList,int dNum);

    //====dti paper experiment 1

    void DTIwarpExp1DisplacementFieldCircular();
    void DTIwarpExperiment1FiberBunch();
    void DTIwarpExp1DisplacementFieldRadiational();

    // === dti paper experiment 1 appendex
    // === explain how optimized process estimate underlying 
    // === fiber orientation
    void DTIwarpExperiment1AppFiberBunch();

    // ====DTI paper with 2 PD to be considered
    void DTIwarpExperiment2PD1FiberBunch(); //create 4 fibers 
	// 3 are the same as DTIwarpExperiment1FiberBunch();
	// another 1 running in Z direction

};


#endif //_VOXMAP_H_


