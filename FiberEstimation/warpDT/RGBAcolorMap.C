#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <strstream>
#include <fstream>

using namespace std;

#include "RGBAcolorMap.h"

//--------------------------------------------------------------------
//  class RGBAcolorMap:
//        RGBA (4-byte at each voxel) 
//
//  
//  members:
//    
//--------------------------------------------------------------------
RGBAcolorMap::RGBAcolorMap(int ddz, int ddx, int ddy) //default ddx=ddy=256
 :dimZ(ddz),dimX(ddx),dimY(ddy),vol(NULL)
{
  planSize = dimX*dimY;  
  int ttt = planSize*dimZ;

  voxSize = sizeof(struct byteRGBA);
  cout <<" voxel size = "<< voxSize<<endl;
  vol = new byteRGBA [ttt];  
  if (vol == NULL) {
    cout <<" ERROR: RGBAcolorMap::RGBAcolorMap() out of memory "<< endl;
    exit(0);
  }

  for (int i=0; i<ttt;i++)
    vol[i].rgba[0] = vol[i].rgba[1] = vol[i].rgba[2] = vol[i].rgba[3] =  0;

}
RGBAcolorMap::RGBAcolorMap(char* name,int ddx, int ddy)
 :dimX(ddx),dimY(ddy),vol(NULL)
{
  planSize = dimX*dimY;  
  voxSize = sizeof(struct byteRGBA);

  loadVolume(name, vol, dimZ);
  cout << dimZ << " slices of "<<dimX<<"*"<< dimY<<" loaded." << endl;  

}

RGBAcolorMap::~RGBAcolorMap()
{
  if (vol != NULL) 
	delete []vol;
}
void RGBAcolorMap::loadVolume(char* ff, byteRGBA* &object, int& sampZ) 
{
   FILE *fp;
   int n,total;
   
   cout << "--opening file: " << ff << endl;
   fp = fopen(ff, "rb");
   if (fp == NULL ) {
       cout << endl << "file open ERROR: "<< ff << endl;
       exit (0);
   }else {       
     //       fseek64(fp,0,SEEK_END);
     fseek(fp,0,SEEK_END);
     // total = ftell64(fp); cout << "total= " <<total<<endl;
       total = ftell(fp); cout << "total= " <<total<<endl;
       sampZ= total/planSize/voxSize; 
       //fseek64(fp,0,SEEK_SET);
       fseek(fp,0,SEEK_SET);
       //---Volume Data--------------
       if (object != NULL){
	 delete []object;
	 object = NULL;
       }
       
       object = new struct byteRGBA[total/voxSize];
       if ( object == NULL ) {
           cout << endl << "   RGBAcolorMap::loadVolume()";
	   cout << endl << "   Memory out when allocating for VOXEL";
	   cout << endl; 
           exit(0);
        }else 
	   cout << "\n Got memory for RGBA volume data\n" <<endl;
	n = fread(object,voxSize, total/voxSize,fp);	
	cout << endl << "\n RGBAcolorMap size = "<< total;
	cout << "; should read = " << total/voxSize << ", while read =" << n << endl;
	fclose(fp);
   }
}
void RGBAcolorMap::setVoxel(int xx, int yy, int zz,
			    char rr, char gg, char bb, char aa)

{
  struct byteRGBA t;

  t.rgba[0]=rr; t.rgba[1]=gg;  t.rgba[2]=bb;  t.rgba[3]=aa;
  setVoxel(xx, yy, zz, t);

}
void RGBAcolorMap::setVoxel(int xx, int yy, int zz,struct byteRGBA cc)
{//set voxel to *vol

   if ( (xx<0)||(xx>dimX-1) ||
	(yy<0)||(yy>dimY-1) ||
	(zz<0)||(zz>dimZ-1) )
   	return ;
   else
 	vol[zz*planSize+yy*dimX+xx]=cc;
}

void RGBAcolorMap::setVoxel(int ind, struct byteRGBA cc)
{//set voxel to *vol

   if (( ind >= planSize*dimZ )||(ind <0)){
	cout <<"ERROR: RGBAcolorMap::setVoxel() voxel index beyond limit"<<endl;
   	return ;
   }else
 	vol[ind]=cc;
}
struct byteRGBA RGBAcolorMap::getVoxel(int ind)//use index instead of (x,y,z)
{
   struct byteRGBA t;

   if (( ind >= planSize*dimZ )||(ind <0)){
	//cout <<"ERROR: RGBAcolorMap::setVoxel() voxel index beyond limit"<<endl;
        t.rgba[0]=t.rgba[1]=t.rgba[2]=t.rgba[3]=0;
   	return t;
   }else
 	return vol[ind];
}

struct byteRGBA RGBAcolorMap::getVoxel(int xx, int yy, int zz)
{//get voxel from *vol
  struct byteRGBA t;
  
   if ( (xx<0)||(xx>dimX-1) ||
	(yy<0)||(yy>dimY-1) ||
	(zz<0)||(zz>dimZ-1) ){
        //cout <<"### OUT!!!##"<<endl;
        t.rgba[0]=t.rgba[1]=t.rgba[2]=t.rgba[3]=0;
   	return t;
   }else{/*
     cout <<" $$$ Inside!! ##"<<endl;
     cout <<" xxyyzz"<< xx<<","<<yy<<","<<zz<<endl;
     cout <<" dimxyz"<< dimX<<","<<dimY<<","<<dimZ<<endl;

          cout <<"  %%%% "<<int(vol[zz*planSize+yy*dimX+xx].rgba[0])<<endl;

     cout <<"  %%%% "<<int(vol[zz*planSize+yy*dimX+xx].rgba[1])<<endl;
     cout <<"  %%%% "<<int(vol[zz*planSize+yy*dimX+xx].rgba[2])<<endl;
     cout <<"  %%%% "<<int(vol[zz*planSize+yy*dimX+xx].rgba[3])<<endl;*/
	return vol[zz*planSize+yy*dimX+xx];
   }
}

int RGBAcolorMap::saveData(char* name, byteRGBA* &voxels, 
			   int xx, int yy, int zz)
{//this procedure can only write data in the format of uncompressed
   FILE *fp;
   int sampSize = xx*yy*zz;
 
   fp = fopen(name,"wb");
   if (fp == NULL ) 
      return -1;
   else {
      fwrite(voxels,voxSize,sampSize,fp);
      fclose(fp);
   }
   return 1;
}

int RGBAcolorMap::saveVolume(char* name)
{
   return saveData( name, vol, dimX, dimY , dimZ);
}

void  RGBAcolorMap::getDimXYZ(int& a, int& b, int& c)
{
    a=dimX; b=dimY; c=dimZ;
}

void  RGBAcolorMap::getZslice(int n, char* & buf)
{
   if ((n < 0) || (n > dimZ-1))
	return;

   byteRGBA *tpt = vol+planSize*n;
   int row, col, k, ibuf;

   ibuf=0;
   for (row = dimY-1; row >=0; row --){
      k = row * dimX; //reverse Y direction to agree with triv24 display
      for (col = 0; col < dimX; col ++){
	 /*
	 if ((n ==3)&& (row==3) &&(col==3)){
		cout << "k="<<k<<"\ntpt="<<(int)(tpt[k].rgba[0])<<","
			<<(int)(tpt[k].rgba[1])<<","
			<<(int)(tpt[k].rgba[2])<<","
			<<(int)(tpt[k].rgba[3])<<endl;
	 }*/
	 buf[ibuf]   =tpt[k].rgba[3];   //R
	 buf[ibuf+1] =tpt[k].rgba[2];   //G
	 buf[ibuf+2] =tpt[k].rgba[1];   //B
	 buf[ibuf+3] =tpt[k].rgba[0];   //A
	 k ++;
	 ibuf += 4;
      }
   }
}
void  RGBAcolorMap::getYslice(int n, char* & buf)
{// view:  Z        
 //        ^
 //        |	   
 //        |
 //        |
 //        |
 //        |
 //        |  
 //      O --------------------->X (usually 256 pixels in X)

   if ((n < 0) || (n > dimY-1))
	return;

   byteRGBA *indxZ=NULL, *indxX=NULL;
   int row, col, k, ibuf;
   /*for (k = 0; k<256*256; k++){
	    buf[k*4]=  0;
	    buf[k*4+1]= 0;
	    buf[k*4+2]=250;
	    buf[k*4+3]=0;
   }*/

   ibuf=0;k=0;
   //indxZ = vol + (dimZ - 1) * planSize + n * dimX;
   indxZ = vol + n * dimX;
   //for (row = dimZ-1; row >=0; row --, k++){
   for (row = 0; row < dimZ; row ++, k++){
      if (k >255 ){ // texture only hold 256x256 now
	 break;
      }else{
         indxX = indxZ;       
         for (col = 0; col < dimX; col ++){
	    buf[ibuf]   =indxX[0].rgba[3];   //R
	    buf[ibuf+1] =indxX[0].rgba[2];   //G
	    buf[ibuf+2] =indxX[0].rgba[1];   //B
	    buf[ibuf+3] =indxX[0].rgba[0];   //A
	    ibuf += 4;
	    indxX ++;
	 }
      }
	indxZ += planSize;
      //indxZ -= planSize;
   }
   
   for (; k < 256; k++){ // fill the rest if dimZ < 256
         for (col = 0; col < 256; col ++){
	    buf[ibuf]=buf[ibuf+1]=buf[ibuf+2]=buf[ibuf+3]=0;
	    //buf[ibuf]= 0;buf[ibuf+1]= 0; buf[ibuf+2]=0; buf[ibuf+3]=0;
	    ibuf +=4;
	 }	
   }
}

void  RGBAcolorMap::getXslice(int n, char* & buf)
{// view:  
 //      O ---------------------> x (usually 256 pixels in X)
 //        |	   
 //        |
 //        |
 //        |
 //        |
 //        | 
 //        \/ Z

   if ((n < 0) || (n > dimY-1))
	return;

   byteRGBA *indxY=NULL, *indxZ=NULL;
   int row, col, ibuf, k;

   ibuf=0; k=0;
   indxZ = vol + n;
   for (row = dimZ-1; row >=0; row --, k++){
      if (k > 255){
	break;
      }else{
        indxY = indxZ; 
        for (col = 0; col < dimY; col ++){
	  if (col>255){
	    break;
	  } else {
	    buf[ibuf]   = indxY[0].rgba[3];   //R
	    buf[ibuf+1] = indxY[0].rgba[2];   //B
	    buf[ibuf+2] = indxY[0].rgba[1];   //G
	    buf[ibuf+3] = indxY[0].rgba[0];   //A
       	    ibuf += 4;
	  }
          indxY += dimX;
        }
        for (;col<256;col++){
	    buf[ibuf]=buf[ibuf+1]=buf[ibuf+2]=buf[ibuf+3]=0;
	    ibuf +=4;
        }
      }
      indxZ += planSize;	
   }
   for (;k < 256; k++){
       for (col = 0; col < 256; col ++){
	    buf[ibuf]=buf[ibuf+1]=buf[ibuf+2]=buf[ibuf+3]=0;
	    ibuf +=4;
       }
   }
   /*
   ibuf=0;
   indxY = vol + (dimY-1)*dimX + n;
   for (row = dimY-1; row >=0; row --){
      indxZ = indxY; 
      for (col = 0; col < dimZ; col ++){
	if (col>255){
	    break;
	} else {
	    buf[ibuf]   = indxZ[0].rgba[3];   //R
	    buf[ibuf+1] = indxZ[0].rgba[2];   //B
	    buf[ibuf+2] = indxZ[0].rgba[1];   //G
	    buf[ibuf+3] = indxZ[0].rgba[0];   //A
       	    ibuf += 4;
	}
        indxZ += planSize;
      }
      for (;col<256;col++){
	    buf[ibuf]=buf[ibuf+1]=buf[ibuf+2]=buf[ibuf+3]=0;
	    ibuf +=4;
      }
      indxY -= dimX;
   }*/
}
