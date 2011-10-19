#ifndef _IMGVOL_H_
#define _IMGVOL_H_

#include "Global.h"
//--------------------------------------------------------------------
//  class ImgVol:
//        Image volume
//
//  members:
//     -getVoxel(int xx, int yy, int zz):
//	 get voxel value from *vol at (xx,yy,zz),
// 
//     -int saveData(char* name, char* &voxels):
//	 save voluem *voxels to file <name>, 
//	 return 1 if succeed; return -1 if failed
//
//  
//--------------------------------------------------------------------


class ImgVol {
protected:
   int dimX, dimY, dimZ,voxSize,planSize;
   char *vol;
 
public:
   ImgVol(int ddz, int ddx=256, int ddy=256);
   ImgVol(char* name,int ddx=256, int ddy=256);
  ~ImgVol();

   void  reset(char c);
   int   inVolume(int xx, int yy, int zz); //return 1 if inside;

   char*  getVolume(void);
   char*  setVolume(char * pt); //replace volume with new volume, return old

   void   loadVolume(char* ff, char* &object, int& sampZ);
   char   getVoxel(int ind); //use index instead of coordinates(x,y,z)
   char   getVoxel(int xx, int yy, int zz);//get voxel from *vol
   void   getDimXYZ(int& aa, int& bb, int& cc);
   int   getDimX(){return dimX;}
   int   getDimY(){return dimY;}
   int   getDimZ(){return dimZ;}

   void   setVoxel(int xx, int yy, int zz, char cc);
   void   setVoxel(int ind,char cc);
   float  smoothness(int xx, int yy, int zz, int nbrSize); 
         // get smoothness value; 0 the smoothest
   void   smooth(int nbrSize); //smooth image with window size "nbrSize"
   float  getVoxelRegionAt(int xx, int yy, int zz, int nbrSize);
   int    getVoxAtAnyPosition(float px, float py, float pz, float & P);
   float  getVoxelf(float vx,float vy,float vz);

   void  getZslice(int n, char* & buf);
   void  getYslice(int n, char* & buf);
   void  getXslice(int n, char* & buf);

   void   sortIntArray(int * arr, int arrLen);
   void   histogram(int * arr, int arrLen);

   int    saveVolume(char* name);
   int    saveData(char* name, char* &voxels, int xx, int yy, int zz);
};


#endif //_VOXMAP_H_
