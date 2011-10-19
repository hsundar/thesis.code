
#ifndef _DTVOLUME_H_
#define _DTVOLUME_H_

#include "Basics.h"
#include "VoxMap.h"
#include "WarpVolume.h"
#include "mVolume.h"

//-------------------------------------------------------------------------
// class DTVolume: Volume of Diffusion Tensor
//    This class is to deal with a volume of diffusion tensor
//    each voxel is a 3x3 matrix :
//		[a1, a4, a5]
//		[a4, a2, a6]
//		[a5, a6, a3]
//    DT volume dimensions:  n1 = 256, n2 = 256, n3 = 52-60
// members:
//    dimX,Y,Z: dimensions of the volume, file size would be
//              dimX*dimY*dimZ*(6*4): 6 float components at a voxel
//    path:     the file name to load the DT volume
//    dt  :     the memory to hold the volume
//    voxSize:  size of the voxel in byte
//    planSize: number of voxel per slice
// member functions:
//    void forwardmapDT2atlas():input map1 contains the map from
//	       original space to the volume before to be warpped.
//	       After this step, the map1 contains the map to the 
//	       destination in atlas space, and the DT in array vol[]
//	       is updated to its oritation in the atlas space
//    void castDT2grid(VoxMap* & map1):
//	       after the forwardmapDT2atlas() step, this step is to
//             get the DT on the 57-slice grid in the atlas space
//	       by calculating the contributions of every DT on 
//	       non-grid positions
// ------------------------------------------------------------------------
//    void reorientTensor(mVolume* &mv):
//           re-orient the DT in the original space, according to 
//	     the transformation matrix volume 'mv'. but each DT is
//	     still in the original space after this re-orientation
//
//    struct DTensor * warpTensorField(VoxMap* & map1,int warpedZ=atlasZ):
//           warp original DT field to atlas space, map1 is the 
//	     warping vector field.
//
//    void averageAndSave(int nSlc, char * toStore) : average DTI
// 	  average the DT volume at a gap of every nSlc slices
// 	  result in a volume of (dimX,dimY,nSlc) DT volume 
// 	  and store to file "toStore"
//	  -in this function, it takes *dt as a stack of DTIs
//   struct Pt3d * PrimaryDirection(char *PDfilename): 
//        - extract Primary Direction:
//        if PDfilename == NULL : return the pointer to the PD vectors
//        if PDfilename != NULL : save the extracted PD vector to file
//			          and return NULL (allocated PD space freed)
//  void saveVectorFile(struct Pt3d * vf, int XX, int YY, int ZZ, char *ff)
//	  save vectors into file $ff, and create an info file $ff.info
//-------------------------------------------------------------------------
class DTVolume{
    int dimX,dimY,dimZ;
    float resX,resY,resZ;   
    DTensor *dt;
    char path[130];
    int voxSize,planSize;

    //-------
    float tmpX,tmpY,tmpZ; //temp varible for internal use

public:
    DTVolume();
    DTVolume(char* name, int dimXX=256, int dimYY=256,
		float resX=1.0,float resY=1.0,float resZ=1.0);
		//float resX=0.9766,float resY=0.9766,float resZ=3.0);
   ~DTVolume();
    void flipZComponent();
    void interlace2normal(); //change interlaced DT volume to normal order
    void normal2interlace(); //change DT volume storag format to interlaced
    void interlace2normal_slicely(); //change slicely-interlaced DT volume to normal order
    void normal2interlace_slicely(); //change normal DT volume to slicely-interlaced order

    void newVolume(int dX,int dY, int dZ, float xr=1.0f, float yr=1.0f, float zr=1.0f);
    void setResolution(float xr, float yr, float zr);

    DTensor *getTensorVol(void){ return dt;};

    void saveData(char* ff, struct DTensor* &object, int numofunit);
    void saveTensorData(char *ff); //save member dt to file 'ff'
    void loadData(char* ff, DTensor * &object, int& sampZ);
    float triInterpolation(float* w, float v0, float v1, 
	float v2, float v3, float v4, float v5, float v6,float v7);

    void  getDimXYZ(int& xx, int& yy, int& zz);
    float interTransGrid(DTensor* & dtGrid, int cn,
	     float wx,float wy, float wz);
    float interTransGrid(DTensor* & dtGrid, int cn,float ww[8]);
    int   isZEROTensor(int index); //return 1 if YES

    float weiInterpolation(float www[8], float A0,  float A1, 
	   float A2, float A3, float A4, float A5, float A6, float A7); 
    void  forwardmapDT2atlas(VoxMap* & map1, WarpVolume* & map2, 
		float rrx, float rry, float rrz);
    struct DTensor *  castDT2grid(VoxMap* & map1);
    int   fillNewDT(struct DTensor * &gridDT, struct Pt3d& grid, 
	    float *www, struct DTensor *newDT); //float www[8]

    int  isTensorWellDefined(struct DTensor * &dtGrid);//useless function
    void survey(struct DTensor * &object, int sampZ);
    void getSusumuMapLastStep(VoxMap* & map1, WarpVolume* & map2);
    double distance(float a, float b, float c, struct Pt3d pp);
    double distanceDebug(float a, float b, float c, struct Pt3d pp);

    void distributeDT(struct DTensor * &oldDT,
		float www, struct DTensor *newDT); //newDT = oldDT *www
    struct DTensor * warpTensorField(VoxMap* & map1, int warpedZ=atlasZ);
    //warp DT, similiar to castDT2grid(), but different algorithm

    struct DTensor * warpTensorFieldRev(VoxMap* & map1, int warpedZ=atlasZ);
    // similar to warpTensorField, but uses reverse vector field

    void reorientTensor(mVolume* &mv); //re-orient the DT in the original space

    void averageAndSave(int nSlc, char * toStore);
    //average DTI
    
    int averageFA_Deviation(int nSlc, float* &averFA, float* &deviation);
    // average DT's FA value and calculate standard deviation
    
    double getDT_FA(float tens[6],float **A, float *eVal, float **eVec); 
    //calculate tens[6]'s FA value

    void createSimulationDT(int XX1, int XX2, 
	int YY1, int YY2, int ZZ1, int ZZ2);// (XX1-2, YY1-2, ZZ1-2) is the subvolume

    float * createFAvolume();

    void flipVolume(int which);

    int setTensorAt(int col, int row, int slice, float vvv[6]);
    int setTensorAt(int col, int row, int slice, DTensor DDD);
    //redifine the tensor 
    void setTensorAt(int index,
	   float v1, float v2,	float v3, 
	   float v4, float v5, 	float v6);
    //add and average the value at the tensor
    void setAverageTensorAt(int index,
	   float v1, float v2,	float v3, 
	   float v4, float v5, 	float v6);
    int fgetTensorAt(float vx,float vy,float vz, DTensor& res);
      //get tensor in a trilinear way at a float position;
      // return 1 when res has result; return 0 if position out of bound
    int getTensorAt(int col, int row, int slice, DTensor &DT1);
      // retrieved DT in DT1 when returned value == 1
      // if returned 0, (col, row, slice ) goes out of volume, DT1 set to 0 
    int getTensorAt(int col, int row, int slice, float **M);
      //M is 3 x 3 matrix
      //return 1 when M is defined; otherwise return 0

   void fatterFAinRange(int ax, int ay, int az, 
        float radius,  float percentage, float acc=0.97);
   //simulate a lower FA at the given location (ax,ay,az)

   struct Pt3d * PrimaryDirection(char *PDfilename=NULL,int whichPD=1);
	//extract PD, which PD defines the 1st, 2nd or 3rd PD to be extracted
   struct Pt3d * PrimaryDirectionDouble(char *PDfilename=NULL);
   void   saveVectorFile(struct Pt3d * vf, int XX, int YY, int ZZ, char *ff);
   void outPt(char* ss,struct Pt3d ppp);
   void printTensors();
   void printTensorAt(int ind);
   void printTensorAt(DTensor *DTS, int ind);
   void printTensorInMatrixAt(int ind);
   void outputs();
};
#endif // _VECTORVOLUME_H_

