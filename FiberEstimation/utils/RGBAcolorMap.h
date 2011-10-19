#ifndef _RGBACOLORMap_H_
#define _RGBACOLORMAP_H_

/**
 * @file RGBAcolorMap.h
 * @brief Container for RGBA colormap
 */

//#include "Global.h"


/** 
 * @brief needed by class RGBAcolorMap 
 */
struct byteRGBA{
  char rgba[4];
}; 

/**
 *  @brief class RGBAcolorMap: volume of RGBA (4-byte at each voxel)
 *
 *  class originally written by Dongrong Xu 
 *
 *  members:
 *       struct byteRGBA 
 *
 *       getVoxel(int xx, int yy, int zz):
 *	 get voxel value from *vol at (xx,yy,zz),
 * 
 *     -int saveData(char* name, byteRGBA* &voxels):
 *	 save voluem *voxels to file <name>, 
 *	 return 1 if succeed; return -1 if failed
 *
 *
*/

class RGBAcolorMap {
protected:
   int dimX, dimY, dimZ,voxSize,planSize;
   byteRGBA *vol;
 
public:

   RGBAcolorMap(int ddz, int ddx=256, int ddy=256);
   RGBAcolorMap(char* name,int ddx=256, int ddy=256);

   ~RGBAcolorMap();

   void loadVolume(char* ff, byteRGBA* &object, int& sampZ);
   void setVoxel(int xx, int yy, int zz,struct byteRGBA cc);
   void setVoxel(int xx, int yy, int zz,char rr, char gg, char bb, char aa);
   void setVoxel(int ind, struct byteRGBA cc);

   struct byteRGBA getVoxel(int xx, int yy, int zz);//use (x,y,z)
   struct byteRGBA getVoxel(int ind); //use index instead of (x,y,z)

   void   getDimXYZ(int& a, int& b, int& c);

   int   getDimX(){return dimX;}
   int   getDimY(){return dimY;}
   int   getDimZ(){return dimZ;}
   void  getZslice(int n, char* & buf);
   void  getYslice(int n, char* & buf);
   void  getXslice(int n, char* & buf);

   int  saveData(char* name, byteRGBA* &voxels,int xx, int yy, int zz);
   int  saveVolume(char* name);


};


#endif //_RGBACOLORMAP_H_
