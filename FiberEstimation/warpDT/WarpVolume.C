#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <iostream.h>
//#include <strstream.h>
#include <strstream>
#include <fstream.h>

#include "Basics.h"
#include "WarpVolume.h"

//--------------------------------------------------------------
//  class WarpVolume member functions
//--------------------------------------------------------------
WarpVolume::WarpVolume(char* fVF, int xx, int yy, 
	float rx,float ry, float rz )
  :dimX(xx),dimY(yy),fdMap(NULL),
   voxSize(sizeof(struct Pt3d)), xres(rx), yres(ry), zres(rz)
{
   strcpy(fWarp, fVF);
   planSize = dimX*dimY;
   loadVolume(fVF, fdMap, dimZ);
   reorganize(fdMap,dimZ);
}

WarpVolume::~WarpVolume()
{
  if (fdMap != NULL)
    delete []fdMap;
}

void WarpVolume::loadVolume(char* ff, Pt3d * &object, int& sampZ) 
{
   FILE *fp;
   int n,total;
   
   cout << "--opening file: " << ff << endl;
   fp = fopen(ff, "rb");
   if (fp == NULL ) {
       cout << endl << "WarpVoluem::loadVolume():file open ERROR: ";
       cout << ff << endl;
       exit (0);
   }else {       
     //fseek64(fp,0,SEEK_END);
     // total = ftell64(fp); cout << "total= " <<total<<endl;
     // sampZ= total/planSize/voxSize; //because we know slice is 256X256
     // fseek64(fp,0,SEEK_SET);

       fseek(fp,0,SEEK_END);
       total = ftell(fp); cout << "total= " <<total<<endl;
       sampZ= total/planSize/voxSize; //because we know slice is 256X256
       fseek(fp,0,SEEK_SET);

       //---Volume Data--------------
       if (object != NULL){
	 delete []object;
	 object = NULL;
       }
       object = new struct Pt3d[total/voxSize];
       if ( object == NULL ) {
           cout << endl << " WarpVoluem::loadVolume()";
	   cout << endl << "   Memory out when allocating for VOXEL";
	   cout << endl; 
           exit(0);
        }else 
	   cout << "\n Got memory for vector volume data\n" <<endl;
	n = fread(object,voxSize, total/voxSize,fp);	
	cout << endl << "\n WarpVolume size = "<< total;
	cout << "; should read = " << total/voxSize << ", while read =" << n << endl;
	fclose(fp);
   }
}
int WarpVolume::JacobiGrid(int col, int row, int slice, float res[9])
{ // return -1 if failed; return 1 if succeeded
   int ind = slice*planSize + row*dimX + col;
   struct Pt3d X0, X1, Y0, Y1, Z0, Z1;

   if (((col<0)||(col>dimX-1)) 
       || ((row<0)||(row>dimY-1)) 
       || ((slice<0)||(slice>dimZ-1)))
     return -1; //no definition
       //cout << "~~~~~~dimZ = " << dimZ << endl;
   //-- X --- Col --
   if (col == 0)
	X0 = fdMap[ind];
   else
	X0 = fdMap[ind-1];
   if (col == dimX-1)
	X1 = fdMap[ind];
   else
	X1 = fdMap[ind+1];
	 // cout << " ttttt---1" <<endl;
   //-- Y --- Row --
   if (row == 0)
	Y0 = fdMap[ind];
   else
	Y0 = fdMap[ind-dimX];
   if (row == dimY-1)
	Y1 = fdMap[ind];
   else
	Y1 = fdMap[ind+dimX];
	  //cout << " \t\t\tttttt---2" <<endl;
  
  //-- Z --- Slice --
  if (slice == 0)
	Z0 = fdMap[ind];
   else 
	Z0 = fdMap[ind-planSize];
	  //cout << " \t\t\tttttt---2.33333" <<endl;
    if (slice == dimZ-1)
	Z1 = fdMap[ind];
   else 
	Z1 = fdMap[ind+planSize];
	  //cout << " \t\t\tttttt---3" <<endl;
 
   //Question: should I use logical voxel unit 2 or (xres+xres)?
   //Answer: (xres+xres) should be used because the reoritation of the
   //tansor is related with the transformation in physical measurement
   res[0] = (X1.x-X0.x) / (xres+xres); //xres is the resolution along X;
   res[1] = (Y1.x-Y0.x) / (yres+yres); //yres is the resolution along X;
   res[2] = (Z1.x-Z0.x) / (zres+zres); //zres is the resolution along X;
	  //cout << " \t\t\tttttt---4" <<endl;

   res[3] = (X1.y-X0.y) / (xres+xres); 
   res[4] = (Y1.y-Y0.y) / (yres+yres);
   res[5] = (Z1.y-Z0.y) / (zres+zres);
	  //cout << " \t\t\tttttt---5" <<endl;

   res[6] = (X1.z-X0.z) / (xres+xres); 
   res[7] = (Y1.z-Y0.z) / (yres+yres);
   res[8] = (Z1.z-Z0.z) / (zres+zres);
	  //cout << " \t\t\tttttt---6" <<endl;

   return 1;
}
int WarpVolume::triInterpolationXYZ(fPoint smpt, struct Pt3d & vii)
{ 
  float x1,y1,z1;
  struct Pt3d ptt;

  smpt.getXYZ(x1,y1,z1);
  ptt.x= x1; ptt.y=y1;ptt.z=z1;
  return triInterpolationXYZ(ptt,vii);
}

int WarpVolume::triInterpolationXYZ(struct Pt3d ptt, struct Pt3d &vii)
{//Tri-linear interpolation in *fdMap according to the coordinate of smpt
 //This function is revised from  class VolumeData::getValue(fPoint pt)
 //  in VData.C
    float vx,vy,vz;
    int row,col,slice;
    struct Pt3d *pp, *pp1,*tpp,*tpp1;
    struct Pt3d v1,v2,v3,v4,v5,v6,v7,v8;
    float tx,ty,tz;

    vx= ptt.x; vy= ptt.y; vz=ptt.z;
    if ( ((vx<0)||(vx>dimX-1)) || ((vy<0)||(vy>dimY-1)) 
        || ((vz<0)||(vz>dimZ-1)))
       return -1; //means value undefined

    col = int(vx); //cout << "col=" << col;
    row = int(vy); //cout << ",row=" << row;
    slice=int(vz); //cout << ",slice=" << slice << endl;
   
    pp = fdMap + planSize*slice + row*dimX + col;
    if (slice == dimZ-1) pp1 = pp;
    else pp1 = pp + planSize;

    v1 = *((struct Pt3d*)pp); 
	//cout << "v1=" << v1.x <<","<< v1.y<<"," <<v1.z <<endl;
    v5 = *((struct Pt3d*)pp1);
	//cout << "v5=" << v5.x<<","<<v5.y<<","<<v5.z << endl;

    if (col == dimX-1){
	v2 = v1;
	v6 = v5;
    } else {
        v2 = *((struct Pt3d*)(pp+1));
    	v6 = *((struct Pt3d*)(pp1+1));
    }
    if (row == dimY-1) {
	v3 = v1;tpp=pp;
	v7 = v5;tpp1=pp1;
    } else {
    	tpp=pp+dimX;   v3 = *((struct Pt3d*)(tpp));
    	tpp1=pp1+dimX; v7 = *((struct Pt3d*)(tpp1));
    }
    if (col == dimX-1){
	v4 = v3;
    	v8 = v7;
    } else {
	v4 = *((struct Pt3d*)(tpp+1));
    	v8 = *((struct Pt3d*)(tpp1+1));
    } 
    //tri-linear interpolation
    tx = vx - col; ty = vy - row; tz = vz - slice;
	/*
    	cout << "v2=" << v2.x << ","<< v2.y <<","<<v2.z <<endl;
      	cout << "v6=" << v6.x << ","<< v6.y <<","<<v6.z <<endl;
	cout << "v4=" << v4.x << ","<< v4.y <<","<<v4.z <<endl;
	cout << "v8=" << v8.x << ","<< v8.y <<","<<v8.z <<endl;

	cout << "tx=" << tx << endl;
	cout << "ty=" << ty << endl;
	cout << "tz=" << tz << endl;*/
	
    //v1 & v3
    v1.x = tx*v2.x+(1-tx)*v1.x;
    v1.y = tx*v2.y+(1-tx)*v1.y;
    v1.z = tx*v2.z+(1-tx)*v1.z; 
    v3.x = tx*v4.x+(1-tx)*v3.x;
    v3.y = tx*v4.y+(1-tx)*v3.y;
    v3.z = tx*v4.z+(1-tx)*v3.z;

    // v5 & v7
    v5.x = tx*v6.x+(1-tx)*v5.x;
    v5.y = tx*v6.y+(1-tx)*v5.y;
    v5.z = tx*v6.z+(1-tx)*v5.z;
    v7.x = tx*v8.x+(1-tx)*v7.x;
    v7.y = tx*v8.y+(1-tx)*v7.y;
    v7.z = tx*v8.z+(1-tx)*v7.z;

    v1.x = ty*v3.x+(1-ty)*v1.x;
    v1.y = ty*v3.y+(1-ty)*v1.y;
    v1.z = ty*v3.z+(1-ty)*v1.z;
    v5.x = ty*v7.x+(1-ty)*v5.x;
    v5.y = ty*v7.y+(1-ty)*v5.y;
    v5.z = ty*v7.z+(1-ty)*v5.z;

    //vii contains the resulting interpolation
    vii.x = tz*v5.x+(1-tz)*v1.x;
    vii.y = tz*v5.y+(1-tz)*v1.y;
    vii.z = tz*v5.z+(1-tz)*v1.z;

    //return fPoint(v1.x,v1.y,v1.z);
    return 1;
}

/*----------------------------------------------------
         z
        /
       4----5 
      /|   /|
     0----1-+------->x
     | 7--+-6
     |/   |/ 
     3----2
     |
     Y

          (d[Ux]/d[x], d[Ux]/d[y], d[Ux]/d[z])
 Jacobi = (d[Uy]/d[x], d[Uy]/d[y], d[Uy]/d[z])
          (d[Uz]/d[x], d[Uz]/d[y], d[Uz]/d[z])
------------------------------------------------------*/
int WarpVolume::Jacobi(float col, float row, float slice,
  	float res[9])
{//calculate at any position (col, row, slice) the Jacobian
 // matrix by interpolating in the 3-d grid 
 // result in res[9] if succeeded and 1 will be returned.
 // return -1 if failed;

   int xx,yy,zz, i,j,k;
   float grid[8][9];
   float txyz[3]; //x, ty, tz;

   if (((col<0)||(col>dimX-1)) 
       || ((row<0)||(row>dimY-1)) 
       || ((slice<0)||(slice>dimZ-1)))
     return -1; //no definition

   xx = (int) col; yy = (int) row; zz = (int)slice;
   txyz[0] = col - xx;  
   txyz[1] = row - yy;  
   txyz[2] = slice - zz;

   	//cout << "xx, yy, zz - (" << xx << "," << yy << "," << zz << endl;
   if (-1 == JacobiGrid(xx  , yy  , zz  , grid[0]))
	return -1;   
    	//cout << "-111" << endl;
   if(-1 == JacobiGrid(xx+1, yy  , zz  , grid[1]))
	return -1;
    	//cout << "-2222" << endl;
   if (-1 == JacobiGrid(xx+1, yy+1, zz  , grid[2]))
	return -1;
    	//cout << "-333, wel=" << welldefined <<endl;
    if (-1 == JacobiGrid(xx  , yy+1, zz  , grid[3]))
	return -1;
     	//cout << "-444,  welldefined = " << welldefined << endl;
    if (-1 == JacobiGrid(xx  , yy  , zz+1, grid[4]))
	return -1;
    	//cout << "-555; wel=" <<  welldefined <<endl;
    if (-1 == JacobiGrid(xx+1, yy  , zz+1, grid[5]))
	return -1;
    	//cout << "-666" << endl;
    if (-1 == JacobiGrid(xx+1, yy+1, zz+1, grid[6]))
	return -1;
     	//cout << "-777" << endl;
    if (-1 == JacobiGrid(xx  , yy+1, zz+1, grid[7]))
 	return -1;
      //cout << " ####   222" << endl;
   k=0;
   for(i=0; i<3; i++){
     for (j=0;j<3;j++){ 
       res[k] = triInterpolation(txyz, 
		     grid[0][k],grid[1][k],
		     grid[2][k],grid[3][k],
		     grid[4][k],grid[5][k],
		     grid[6][k],grid[7][k]);

       k++;
     }
   }

   return 1;
}
/*----------------------------------------------------
         z
        /
       4----5 
      /|   /|
     0----1-+------->x
     | 7--|-6
     |/   |/ 
     3----2
     |
     Y
 w[0]: weightX; w[1]: weightY; w[1]: weightZ;
----------------------------------------------------*/
float WarpVolume::triInterpolation(float* w, float v0, 
	float v1, float v2, float v3, float v4, float v5, 
	float v6,float v7)
{
   float t01,t32,t0132,t45,t76,t4576,res;

   t01 = v0*w[0] + v1*(1.0-w[0]);
   t32 = v3*w[0] + v2*(1.0-w[0]); 
   t0132 = t01*w[1] + t32*(1.0-w[1]);

   t45 = v4*w[0] + v5*(1.0-w[0]);
   t76 = v7*w[0] + v6*(1.0-w[0]); 
   t4576 =t45*w[1] + t76*(1.0-w[1]);

   res = t0132*w[2] + t4576 * (1.0-w[2]);

   return res;
}
int WarpVolume::JacobiRotationMatrix(fPoint pt, Matrix& M)
{
  float x1,y1,z1;

  pt.getXYZ(x1,y1,z1);
  return JacobiRotationMatrix(x1, y1, z1, M);
}

int WarpVolume::JacobiRotationMatrix(struct Pt3d pt, Matrix& M)
{
  return JacobiRotationMatrix(pt.x, pt.y, pt.z, M);
}
int WarpVolume::JacobiRotationMatrix(float xx, float yy, float zz, Matrix& M)
{// (xx, yy, zz): the 3-D location 
 // return 1 : successful
 // return -1: failed
 // resulting matrix: M

  float Jab[9], ang, ox,oy,oz;
  Vector I(1,0,0), J(0,1,0), K(0,0,1), omega;
  	//cout <<" **111****" << endl;
  int res=Jacobi(xx,yy,zz, Jab);
  	//cout <<" **222****" << endl;

  if (res == -1)
	return -1;
 
  omega = -(I*Jab[5]) - (J*Jab[6]) - (K*Jab[1]);
	//cout <<endl <<"omega=>"; omega.outputs();
  omega.getXYZ(ox,oy,oz);
	//cout << endl << "omega <==" << endl;
  ang = omega.length(); 
     //cout << "omega & ang :" << endl;
     //omega.outputs(); cout << "ang = " << ang << endl;
  
  M.loadIdentity();
  M.setValue(0 , 0, ox*ox*(1-cos(ang)) +    cos(ang));
  M.setValue(0 , 1, ox*oy*(1-cos(ang)) - oz*sin(ang));
  M.setValue(0 , 2, ox*oz*(1-cos(ang)) + ox*sin(ang));
  M.setValue(1 , 0, ox*oy*(1-cos(ang)) + oz*sin(ang));
  M.setValue(1 , 1, oy*oy*(1-cos(ang)) +    cos(ang));
  M.setValue(1 , 2, oy*oz*(1-cos(ang)) - ox*sin(ang));
  M.setValue(2,  0, ox*oz*(1-cos(ang)) - oy*sin(ang));
  M.setValue(2 , 1, oy*oz*(1-cos(ang)) + ox*sin(ang));
  M.setValue(2 , 2, oz*oz*(1-cos(ang)) +    cos(ang));

  return 1;
}

void WarpVolume::reorganize(struct Pt3d*  &volume, int volZZ)
{//Christos took horizontal direction as Y, vert -> X;
 //and, also we need to replace the destination vector with
 //displacement vector
 //so, this step is necessary before proceeding
  int XX,YY,ZZ,k=0;
  float tp;

  cout << "reorganize() ZZ = " << volZZ;
  cout << "; XX, YY = " << dimX << "," << dimY << endl;
  for (ZZ = 0; ZZ < volZZ; ZZ ++){
   //cout << "ZZ = " << ZZ << endl;
   for (YY = 0; YY < dimY; YY ++)
     for (XX = 0; XX < dimX; XX ++){
        tp = volume[k].x;
	volume[k].x = volume[k].y; //-XX;
	volume[k].y = tp; //-YY;
	//volume[k].z -= ZZ;
	k++;
     }
   }//volume[].x/y/z was with the destination value,
    // now it is with the displacement value
  cout << " At the end of reorganize()" << endl;
}

struct Pt3d WarpVolume::getVectorAt(int col, int row, int slice)
{ //To use this function, "this" must not be a displacement field
  // but a vector field of absolute coordinates
  // i.e. 2nd part of reorganization() should be avoided.

   int ind = slice*planSize + row*dimX + col;
   struct Pt3d P;

   if (((col<0)||(col>dimX-1)) 
       || ((row<0)||(row>dimY-1)) 
       || ((slice<0)||(slice>dimZ-1))){
	P.x = P.y = P.z = UNDEFINED_PT;
	return P;
   }else
	return fdMap[ind];
} 
