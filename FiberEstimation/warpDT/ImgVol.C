#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <iostream.h>
//#include <strstream.h>
#include <strstream>
#include <fstream.h>

#include "ImgVol.h"

//--------------------------------------------------------------------
//  class ImgVol:
//        load in the original image and downsample it 2x2x2
//        result in *nVol
//  
//  members:
//    
//--------------------------------------------------------------------
ImgVol::ImgVol(int ddz, int ddx, int ddy) //default ddx=ddy=256
 :dimZ(ddz),dimX(ddx),dimY(ddy),vol(NULL),voxSize(1)
{
  planSize = dimX*dimY;  
  int ttt = planSize*dimZ;
  vol = new char [ttt];  
  if (vol == NULL) {
    cout <<" ERROR: ImgVol::ImgVol() memory out"<< endl;
    exit(0);
  }
  for (int i=0; i<ttt;i++)
	vol[i] = 0;
}
ImgVol::ImgVol(char* name,int ddx, int ddy)
 :dimX(ddx),dimY(ddy),vol(NULL),voxSize(1)
{
  planSize = dimX*dimY;  

  loadVolume(name, vol, dimZ);
  cout << dimZ << " slices of plansize loaded" << endl;

}

ImgVol::~ImgVol()
{
  if (vol != NULL) 
	delete []vol;
}
void ImgVol::reset(char c)
{
  int ttt = planSize*dimZ;
  for (int i=0; i<ttt;i++)
	vol[i] = c;
}
char* ImgVol::getVolume(void)
{ 
  return vol;
}
char* ImgVol::setVolume(char * pt)
{ 
  char *tt= vol;
  vol =  pt;
  return tt;
}

void ImgVol::loadVolume(char* ff, char* &object, int& sampZ) 
{
   FILE *fp;
   int n,total;
   
   cout << "--opening file: " << ff << endl;
   fp = fopen(ff, "rb");
   if (fp == NULL ) {
       cout << endl << "file open ERROR: "<< ff << endl;
       exit (0);
   }else {       
     //fseek64(fp,0,SEEK_END);
     //total = ftell64(fp); cout << "total= " <<total<<endl;
     //sampZ= total/planSize/voxSize; //because we know slice is 256X256
     //fseek64(fp,0,SEEK_SET);

       fseek(fp,0,SEEK_END);
       total = ftell(fp); cout << "total= " <<total<<endl;
       sampZ= total/planSize/voxSize; //because we know slice is 256X256
       fseek(fp,0,SEEK_SET);
       //---Volume Data--------------
       if (object != NULL){
	 delete []object;
	 object = NULL;
       }
       
       object = new char[total/voxSize];
       if ( object == NULL ) {
           cout << endl << "VoxMap::loadVolume()";
	   cout << endl << "   Memory out when allocating for VOXEL";
	   cout << endl; 
           exit(0);
        }else 
	   cout << "\n Got memory for byte volume data\n" <<endl;
	n = fread(object,voxSize, total/voxSize,fp);	
	cout << endl << "\n WarpVolume size = "<< total;
	cout << "; should read = " << total/voxSize << ", while read =" << n << endl;
	fclose(fp);
   }
}
void ImgVol::setVoxel(int xx, int yy, int zz,char cc)
{//set voxel to *vol

   if ( (xx<0)||(xx>dimX-1) ||
	(yy<0)||(yy>dimY-1) ||
	(zz<0)||(zz>dimZ-1) )
   	return ;
   else
 	vol[zz*planSize+yy*dimX+xx]=cc;
}
void ImgVol::setVoxel(int ind,char cc)
{//set voxel to *vol

   if (( ind >= planSize*dimZ )||(ind <0)){
	cout <<"ERROR: ImgVol::setVoxel() voxel index beyond limit"<<endl;
   	return ;
   }else
 	vol[ind]=cc;
}
char ImgVol::getVoxel(int ind) //use index instead of coordinates(x,y,z)
{
  if (( ind >= planSize*dimZ )||(ind <0)){
	//cout <<"ERROR: ImgVol::getVoxel() voxel index beyond limit"<<endl;
   	return 0;
   }else
 	return vol[ind];
}
int ImgVol::inVolume(int xx, int yy, int zz)
{//return 1 if inside;
   if ( (xx<0)||(xx>dimX-1) ||
	(yy<0)||(yy>dimY-1) ||
	(zz<0)||(zz>dimZ-1) )
   	return 1;
   else 
	return 0;
}
char ImgVol::getVoxel(int xx, int yy, int zz)
{//get voxel from *vol

   if ( (xx<0)||(xx>dimX-1) ||
	(yy<0)||(yy>dimY-1) ||
	(zz<0)||(zz>dimZ-1) )
   	return char(0);
   else
	return vol[zz*planSize+yy*dimX+xx];
}

int ImgVol::saveData(char* name, char* &voxels, int xx, int yy, int zz)
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

int ImgVol::saveVolume(char* name)
{
   return saveData( name, vol, dimX, dimY , dimZ);
}
 void  ImgVol:: getDimXYZ(int& aa, int& bb, int& cc)
{
   aa=dimX;bb=dimY;cc=dimZ;
}

float ImgVol::getVoxelRegionAt(int xx, int yy, int zz, int nbrSize)
{//get an average voxel value at (xx,yy,zz)
   int i,j,k; //, nnn=0;
   float sum=0.0f; //, weit;
   
   if (nbrSize ==0)
	return (float)getVoxel(xx,yy,zz);
   for (k = zz-nbrSize; k<=zz+nbrSize; k ++)
   for (j = yy-nbrSize; j<=yy+nbrSize; j ++)
   for (i = xx-nbrSize; i<=xx+nbrSize; i ++){
       if (inVolume(i,j,k)){
         if ( i == xx && j == yy && k == zz)
	 continue;
         sum += getVoxel(i,j,k);
         //nnn ++;
       }
   }
  
   if (sum >0)
      sum += sum*0.6+ 0.4 * getVoxel(xx,yy,zz);
   else 
      sum = getVoxel(xx,yy,zz);
   return sum;
   /*
   if (nnn >0)
	return float (sum/nnn);
   else 
	return 0.0f;
   */

}
float ImgVol::smoothness(int xx, int yy, int zz, int nbrSize)
{// return 0 for good smoothness

   int i,j,k, nnn=0;
   float sum=0.0f;
   char me=getVoxel(xx,yy,zz);

   for (k = zz-nbrSize; k<=zz+nbrSize; k ++)
   for (j = yy-nbrSize; j<=yy+nbrSize; j ++)
   for (i = xx-nbrSize; i<=xx+nbrSize; i ++){
       sum += fabs(getVoxel(i,j,k)-me);
       nnn ++;
   }
   if (nnn >0)
	return float (sum/nnn);
   else 
	return 0.0f;

}
void ImgVol::smooth(int nbrSize)
{
   int xx,yy,zz, i,j,k;
   int nnn,ind;
   float sum;
   char me;

   char *tMem=   new char[dimX*dimY*dimZ];

   ind =0;
   for (zz = 0; zz< dimZ; zz++){
     for (yy = 0; yy< dimY; yy++){
       for (xx = 0; xx< dimX; xx++) {
	   sum = 0;nnn =0;
	   for (k = zz-nbrSize; k<=zz+nbrSize; k ++)
	   for (j = yy-nbrSize; j<=yy+nbrSize; j ++)
	   for (i = xx-nbrSize; i<=xx+nbrSize; i ++){
	      me = getVoxel(i,j,k);
	      if (me > 0){
	        sum += me;
	        nnn ++;
	      }
	   }
	   sum += getVoxel(xx,yy,zz);nnn++;
	   sum = sum / nnn;
	   tMem[ind]=(char)sum;
	   ind ++;
       }
     }
  }
  delete []vol;
  vol = tMem;
  saveVolume("afterSmooth.img");
}
void ImgVol::sortIntArray(int * arr, int arrLen)
{
  int t, i,j;
  
  for (i=0; i<arrLen-1; i++)
    for (j=i+1; j<arrLen; j++){
       if (arr[i]<arr[j]){
	   t = arr[i]; arr[i] = arr[j]; arr[j] = t;
       }
    }
}
void ImgVol::histogram(int * arr, int arrLen)
{
  int i, hist[256];
  int xx,yy,zz,ind;
  
  for (i=0;i<256; i++)
    hist[i]=0;

  ind =0;
  for (zz = 0; zz< dimZ; zz++){
     for (yy = 0; yy< dimY; yy++){
       for (xx = 0; xx< dimX; xx++) {
         hist[vol[ind]] ++;
         ind ++;
       }
     }
  }
  int j, max=-1,jmax=0;
  int nSS, histt[256];
  float ss;
  for (j=0;j<256; j++){
    ss = 0;nSS=0;
    for (i=j-1; i<=j+1; i++){
	if (i<0) continue;
	if (i>255) continue;
        nSS ++;
	ss += hist[i];
    }
    ss += hist[j]; nSS ++;
    histt[j] = ss / nSS;
  }
  for (i=0; i<arrLen; i++){
    for (j=0;j<256; j++){
	if (max < hist[j]){
	   max = hist[j];
	   jmax=j;
        }
    }
    arr[i] = jmax;
  }
}

/*

              |Y
              |
           v3 |___________ v4
             /|          /|
            / |         / |
         v7/__|_______ /v8|
           |  |        |  |
           |  |_______ |_ |___________X
           | / v1      | / v2
         v5|/__________|/  
           /           v6    
          /         
         /
        /
       /
      Z

*/

float ImgVol::getVoxelf(float vx,float vy,float vz)
{
    int row,col,slice;
    char *pp, *pp1,*tpp,*tpp1;

    float v1,v2,v3,v4,v5,v6,v7,v8,tx,ty,tz;

    if ((vx<0)||(vx>dimX-1)) return 0.0f;
    if ((vy<0)||(vy>dimY-1)) return 0.0f;
    if ((vz<0)||(vz>dimZ-1)) return 0.0f;

    col = int(vx);
    row = int(vy);
    slice=int(vz);
   
    pp = vol + planSize*slice + row*dimX + col;
    if (slice == dimZ-1) pp1 = pp;
    else pp1 = pp + planSize;

    v1 = unsigned (*pp);
    v5 = unsigned (*pp1);

    if (col == dimX-1){
	v2 = v1;
	v6 = v5;
    } else {
        v2 = unsigned (*(pp+1));
    	v6 = unsigned (*(pp1+1));
    }
    if (row == dimY-1) {
	v3 = v1;tpp=pp;
	v7 = v5;tpp1=pp1;
    } else {
    	tpp=pp+dimX;   v3 = unsigned(*(tpp));
    	tpp1=pp1+dimX; v7 = unsigned(*(tpp1));
    }
    if (col == dimX-1){
	v4 = v3;
    	v8 = v7;
    } else {
	v4 = unsigned(*(tpp+1));
    	v8 = unsigned(*(tpp1+1));
    } 
    tx = vx - col; ty = vy - row; tz = vz - slice;
    v1 = tx*v2+(1-tx)*v1;
    v3 = tx*v4+(1-tx)*v3;

    v5 = tx*v6+(1-tx)*v5;
    v7 = tx*v8+(1-tx)*v7;

    v1 = ty*v3+(1-ty)*v1;
    v5 = ty*v7+(1-ty)*v5;

    v1 = tz*v5+(1-tz)*v1;

    return v1;
}
int ImgVol::getVoxAtAnyPosition(float px, float py, float pz, float & P)
{//Tri-linear interpolation in *vol according to the coordinate of smpt
 //This volume have dimZ slices in Z direction
 // return 1 when P is properly assigned
 // return 0 when P is undefined (but set to be a 0-vector)

    float vx,vy,vz;
    int row,col,slice;
    char *pp, *pp1,*tpp,*tpp1;
    char v1,v2,v3,v4,v5,v6,v7,v8;
    float tx,ty,tz;
    char* volume=vol;
    int sampZ = dimZ;

    vx=px;vy=py;vz=pz;

	//cout <<"   BBB:  AnyPos: ("<<vx<<","<<vy<<","<<vz<<")"<<endl;
    if ( ((vx<0)||(vx>dimX-1)) || ((vy<0)||(vy>dimY-1)) 
        || ((vz<0)||(vz>sampZ-1))){
       P = 0;
       return  0; //means value undefined
    }

    col = int(vx); //cout << "col=" << col;
    row = int(vy); //cout << ",row=" << row;
    slice=int(vz); //cout << ",slice=" << slice << endl;
   
    pp = volume + planSize*slice + row*dimX + col;
    if (slice == sampZ-1) pp1 = pp;
    else pp1 = pp + planSize;

    v1 = *((unsigned short*)pp); //cout << "v1=" << v1.x <<","<< v1.y<<"," <<v1.z <<endl;
    v5 = *((unsigned short*)pp1);//cout << "v5=" << v5.x<<","<<v5.y<<","<<v5.z << endl;

    if (col == dimX-1){
	v2 = v1;
	v6 = v5;
    } else {
        v2 = *((unsigned short*)(pp+1));
    	v6 = *((unsigned short*)(pp1+1));
    }
    if (row == dimY-1) {
	v3 = v1;tpp=pp;
	v7 = v5;tpp1=pp1;
    } else {
    	tpp=pp+dimX;   v3 = *((unsigned short*)(tpp));
    	tpp1=pp1+dimX; v7 = *((unsigned short*)(tpp1));
    }
    if (col == dimX-1){
	v4 = v3;
    	v8 = v7;
    } else {
	v4 = *((unsigned short*)(tpp+1));
    	v8 = *((unsigned short*)(tpp1+1));
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
    v1 = tx*v2+(1-tx)*v1;
    v3 = tx*v4+(1-tx)*v3;

    // v5 & v7
    v5 = tx*v6+(1-tx)*v5;
    v7 = tx*v8+(1-tx)*v7;

    v1 = ty*v3+(1-ty)*v1;
    v5 = ty*v7+(1-ty)*v5;

    v1 = tz*v5+(1-tz)*v1;

    P=v1; 

    return 1;
}
void  ImgVol::getZslice(int n, char* & buf)
{
   if ((n < 0) || (n > dimZ-1))
	return;

   char *tpt = vol+planSize*n;
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
	 buf[ibuf]   =tpt[k];   //R
	 buf[ibuf+1] =tpt[k];   //G
	 buf[ibuf+2] =tpt[k];   //B
	 buf[ibuf+3] =255;   //A
	 k ++;
	 ibuf += 4;
      }
   }
}
void  ImgVol::getYslice(int n, char* & buf)
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

   char *indxZ=NULL, *indxX=NULL;
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
	    buf[ibuf]   =indxX[0];   //R
	    buf[ibuf+1] =indxX[0];   //G
	    buf[ibuf+2] =indxX[0];   //B
	    buf[ibuf+3] =indxX[0];   //A
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

void  ImgVol::getXslice(int n, char* & buf)
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

   char *indxY=NULL, *indxZ=NULL;
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
	    buf[ibuf]   = indxY[0];   //R
	    buf[ibuf+1] = indxY[0];   //B
	    buf[ibuf+2] = indxY[0];   //G
	    buf[ibuf+3] = indxY[0];   //A
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
}
