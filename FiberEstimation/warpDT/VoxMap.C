#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <iostream.h>
//#include <strstream.h>
#include <strstream>
#include <fstream.h>

#include <time.h>

#include "Basics.h"
#include "VoxMap.h"
#include "ImgVol.h"

//--------------------------------------------------------------
//  class VoxMap member functions
//--------------------------------------------------------------
VoxMap::VoxMap()
  :fdMap(NULL)
{
   faces[0][0]=0;   faces[0][1]=1;   faces[0][2]=3;
   faces[1][0]=2;   faces[1][1]=1;   faces[1][2]=3;
   faces[2][0]=1;   faces[2][1]=5;   faces[2][2]=2;
   faces[3][0]=2;   faces[3][1]=5;   faces[3][2]=6;
   faces[4][0]=0;   faces[4][1]=1;   faces[4][2]=5;
   faces[5][0]=0;   faces[5][1]=4;   faces[5][2]=5;
   faces[6][0]=4;   faces[6][1]=5;   faces[6][2]=7;
   faces[7][0]=5;   faces[7][1]=6;   faces[7][2]=7;
   faces[8][0]=0;   faces[8][1]=4;   faces[8][2]=7;
   faces[9][0]=0;   faces[9][1]=7;   faces[9][2]=3;
   faces[10][0]=2;  faces[10][1]=7;  faces[10][2]=3;
   faces[11][0]=2;  faces[11][1]=7;  faces[11][2]=6;
}
VoxMap::VoxMap(char* fCrop, char* fVF, int LOADEXISTING, 
	   int xx, int yy, int zz,
  	   float rx,float ry, float rz)
  :dimX(xx),dimY(yy),dimZ(zz),fdMap(NULL),
   voxSize(sizeof(struct Pt3d)), rX(rx), rY(ry), rZ(rz)
{
   int t = time(NULL);
   srand48(t); //initialize random generator  
   //int(drand48()*lenX)

   strcpy(fCropInfo, fCrop);
   strcpy(fWarp, fVF);
   planSize = dimX*dimY;
   if (LOADEXISTING){
     loadVolume(fVF,fdMap,dimZ);
     cout <<"VoxMap::loadVolume() <"<<fVF<<"> dimension=(";
     cout <<dimX<<","<<dimY<<","<<dimZ<<")"<<endl;
   }else
     initVox();
   faces[0][0]=0;   faces[0][1]=1;   faces[0][2]=3;
   faces[1][0]=2;   faces[1][1]=1;   faces[1][2]=3;
   faces[2][0]=1;   faces[2][1]=5;   faces[2][2]=2;
   faces[3][0]=2;   faces[3][1]=5;   faces[3][2]=6;
   faces[4][0]=0;   faces[4][1]=1;   faces[4][2]=5;
   faces[5][0]=0;   faces[5][1]=4;   faces[5][2]=5;
   faces[6][0]=4;   faces[6][1]=5;   faces[6][2]=7;
   faces[7][0]=5;   faces[7][1]=6;   faces[7][2]=7;
   faces[8][0]=0;   faces[8][1]=4;   faces[8][2]=7;
   faces[9][0]=0;   faces[9][1]=7;   faces[9][2]=3;
   faces[10][0]=2;  faces[10][1]=7;  faces[10][2]=3;
   faces[11][0]=2;  faces[11][1]=7;  faces[11][2]=6;
}

VoxMap::~VoxMap()
{
  if (fdMap != NULL)
    delete []fdMap;
}

void VoxMap::initVox()
{
   int XX, YY, ZZ, k;

   if (fdMap != NULL)
      delete []fdMap;
   fdMap = new struct Pt3d[planSize*dimZ];
   if (fdMap == NULL){
	cout << " VoxMap::initVox(): memory out" << endl;
	exit(1);
   }
   k=0;
   for (ZZ = 0; ZZ < dimZ; ZZ ++){
     for (YY = 0; YY < dimY; YY ++){
       for (XX = 0; XX < dimX; XX ++){
	  fdMap[k].x = XX; 
	  fdMap[k].y = YY;
	  fdMap[k++].z = ZZ; 
       }	
     }
   }
}

void VoxMap::loadCropInfo(int& infoZ, int& infoCrop1, int& infoCrop2)
{
  FILE *fp;

  fp = fopen(fCropInfo, "r");
  if (fp == NULL) {
     cout << endl << "VectorVolume::loadCropInfo(): ERROR" << endl;
     cout << "    file <" << fCropInfo  << "> open failed" << endl;
     exit(1);
  }
  fscanf(fp,"%d %d %d", &infoZ, &infoCrop1, &infoCrop2);
  // original-Z, CroppedSliceNum-inHead, CroppedSliceNum-inTail
  fclose(fp);
}

void VoxMap::createForwardMap(float rx, float ry, float rz)
{
   createSusumuMapStep1(rx,ry,rz);
   cout << "forward Map completed." << endl;
}
/*{
   Matrix Ms1,M,M90,Mx,My,Mz,Ms2;
   int XX, YY, ZZ,k;
   int mCropped, nCropped, warpZ, volZ;
   //struct Pt3d *vfWarp=NULL;

   Ms1.newScaleXYZ(1.0,1.0,1.5/0.9766);  
   M90.loadIdentity(); M90.rotateX(-90/180.0*PI);//rotateX -90
   rx = rx/180*PI;  ry = ry/180*PI;  rz = rz/180*PI;
   Mx.loadIdentity();  My.loadIdentity();  Mz.loadIdentity();
   Mx.rotateX(rx);  My.rotateY(ry);  Mz.rotateZ(rz);
   Ms2.newScaleXYZ(1.0,1.0, 0.9766/1.5);
   M = Ms1*M90*Mx*My*Mz*Ms2;
 
   cout << "fWarp = " <<fWarp<< endl;
   //loadVolume(fWarp, vfWarp, warpZ);   
   warpZ = getVolumeZ(fWarp);
   loadCropInfo(volZ, nCropped,mCropped);
      cout << endl << "--- getForwardMap: Crop Info:--- " << endl;
      cout << "(volZ, nCropped, mCropped)="<<volZ <<","<< nCropped;
      cout << "," << mCropped << endl;
   if (mCropped-nCropped+1 != warpZ){//confirm the data is correct
	cout << endl << "ERROR:VoxMap::getForwardMap()" << endl;
	cout << "     file size does not match expected,warpZ="<<warpZ << endl;
	exit(1);
   } else{
	cout << endl << "warpZ=" <<warpZ<<endl;
 	cout << "*.bottop information correct!!" << endl <<endl;
   }
   //switchXY(vfWarp,warpZ); //Christos's vector field use different X Y pair
   //resetVFvolume(vfWarp, warpZ);// for debugging purpose
   k=0;
   for (ZZ = 0; ZZ < dimZ; ZZ ++){
     for (YY = 0; YY < dimY; YY ++){
       for (XX = 0; XX < dimX; XX ++){
              //double slice: e.g. 56 -> 111 then 112-slice
	  fPoint pt(fdMap[k].x, fdMap[k].y, fdMap[k].z*2.0); //fPoint pt(XX,YY,ZZ*2.0);
		//pt.setXYZ(127.5f,138.5f,55.5f);cout << endl <<endl << "111: "; pt.outputs();
	      // move Origin to the volume center
	  pt = pt + fPoint(float((1-dimX)/2.0), float((1-dimY)/2.0), float(-dimZ+0.5));
		//cout << endl << endl <<"222: "; pt.outputs();
	      // -90 X-axel, then rx, ry, and rz rotation
	      //reslice: Y:112->256; Z: 256->112
	  pt = pt * M;  //cout << endl <<endl << "333: "; pt.outputs();
	      //move the ORIGIN point back to corner
	  pt = pt + fPoint(float((dimX-1)/2.0), float((dimY-1)/2.0), float(dimZ/2.0-0.5-nCropped));
		//cout << endl << endl <<"444: "; pt.outputs();
	      //Z deduct nCropped is because the VF file has less Z slices
	      //warp to 256*256*57 atlas space
	      //---1. vector field or destination definition in the VF file
	      //        if vector field, mathmatical addtion required
	      //	if destination coordinates, then replacement but addtion is to be used
	      //---2. how about the point of undefined value?
	      //Answer: Seems to be destination
 	      //since the real VF file has less Z slices
              //pt = pt + fPoint(0.0f, 0.0f, float(-nCropped));
	      // tri-interpolation
          //pt = triInterpolationXYZ(vfWarp,warpZ,pt); //how about pt is undefined?
	    //cout << endl << endl <<"555: "; pt.outputs();	  
          //if (pt.isUndefined()) {
	       //fdMap[k].x = fdMap[k].y = fdMap[k].z =-100.0; //let the value out of range, denoting undefined.
	       //ct_undefined ++;
	  //}else
	       pt.getXYZ(fdMap[k].x, fdMap[k].y,fdMap[k].z);
	  k++;
       }
     }
  }
  //delete []vfWarp;
  //cout << " Undefined pt="<<ct_undefined << endl;
  //cout << " total pt="<<k<<endl;
  //cout << " Forward map: undefined = " << 100.0 * ((float)(ct_undefined )/(float)k);
  //cout << " %"  << endl;
  cout << "forward Map completed." << endl;
}
*/
struct Pt3d* &VoxMap::getMap()
{
  return fdMap;
}
void VoxMap::saveVecVolume(struct Pt3d * &volume, int volZ, char * where)
{
   FILE * fp = NULL;
   int n,ss = planSize*volZ;

   cout << "VoxMap::saveVecVolume():volZ="<<volZ<<endl;
   cout << "volume * = " << volume << endl;
   fp = fopen(where, "wb");
   if(fp==NULL) {
      cout << endl <<"ERROR: VectorVolume::saveVecVolume():file open failed" << endl;
      cout << "file open failed: <" << where << "> " << endl;
      return;
   }
   n=fwrite(volume, voxSize, ss,fp); 
   if (n != ss) {
       cout << endl << " ERROR: VoxMap::saveVecVolume():writing data failed" << endl;       
       cout << "file saved size = " << n << ", while " << ss << " expected." << endl;
   }else 
       cout << "file saved correctly" << endl;
   fclose (fp);
   cout << endl << "file <" << where << "> saved " <<endl;
}
void VoxMap::saveMapVolume(char* name)
{
  cout << "fdMap * = " << fdMap <<endl;
  saveVecVolume(fdMap,dimZ,name);
}
void VoxMap::reorganize(int XYswitch, int ADD_REMOVE_GRID)
{
  int XX,YY,ZZ,k;
  float tp;

  cout << "VoxMap::reorgainze() ZZ = " << dimZ;
  //cout << "; XX, YY = " << dimX << "," << dimY << endl;
  if (XYswitch){
    cout << "  ...switching X-Y values ...."<<endl;
    k = 0;
    for (ZZ = 0; ZZ < dimZ; ZZ ++){
      cout << "ZZ = " << ZZ <<"  "<< char(13);flush(cout);
      for (YY = 0; YY < dimY; YY ++)
      for (XX = 0; XX < dimX; XX ++){
        tp = fdMap[k].x;
	fdMap[k].x = fdMap[k].y;
	fdMap[k++].y = tp;
      }
     }
     cout<<endl;
   }
   if (ADD_REMOVE_GRID==1){//relative coordinates -> absolution destination
    cout <<" ...relative displacement -> absolute coordinates"<<endl;
    k = 0;
    for (ZZ = 0; ZZ < dimZ; ZZ ++){
      cout << "ZZ = " << ZZ <<"  "<< char(13); flush(cout);
      for (YY = 0; YY < dimY; YY ++)
      for (XX = 0; XX < dimX; XX ++){
	fdMap[k].x += XX;
	fdMap[k].y += YY;
 	fdMap[k++].z += ZZ;
      }
     }
     cout <<endl;
   }else if (ADD_REMOVE_GRID== -1){// absolution coordinates -> relative coor
    cout << " ... absolute coordinates -> relative displacement....  " << char(13);
    flush(cout);
    k = 0;
    for (ZZ = 0; ZZ < dimZ; ZZ ++){
      cout << "ZZ = " << ZZ << endl;
      for (YY = 0; YY < dimY; YY ++)
      for (XX = 0; XX < dimX; XX ++){
 	fdMap[k].x -= XX;
	fdMap[k].y -= YY;
 	fdMap[k++].z -= ZZ;
      }
     }
     cout <<endl<<endl;
   }
}

void VoxMap::switchXY(struct Pt3d*  &volume, int volZZ)
{//Christos took horizontal direction as Y, vert -> X
 //so, this step is necessary before proceeding
  int XX,YY,ZZ,k=0;
  float tp;

  cout << "switchXY() ZZ = " << volZZ;
  cout << "; XX, YY = " << dimX << "," << dimY << endl;
  for (ZZ = 0; ZZ < volZZ; ZZ ++){
   //cout << "ZZ = " << ZZ << endl;
   for (YY = 0; YY < dimY; YY ++)
     for (XX = 0; XX < dimX; XX ++){
        tp = volume[k].x;
	volume[k].x = volume[k].y;
	volume[k++].y = tp;
     }
   }
  cout << " At the end of SwitchXY()" << endl;
}

int VoxMap::getVolumeZ(char* ff)
{
   FILE *fp;
   int  total,sampZ;

   fp = fopen(ff, "rb");

   if (fp == NULL ) {
       cout << endl << "file open ERROR: "<< ff << endl;
       return -1;
   }else {       
     //fseek64(fp,0,SEEK_END);
     //total = ftell64(fp); cout << "total= " <<total<<endl;
     // sampZ= total/planSize/voxSize; //because we know slice is 256X256

       fseek(fp,0,SEEK_END);
       total = ftell(fp); cout << "total= " <<total<<endl;
       sampZ= total/planSize/voxSize; //because we know slice is 256X256

       fclose(fp);
       return sampZ;
   }

}
void VoxMap::loadVolume(char* ff, Pt3d * &object, int& sampZ) 
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
       object = new struct Pt3d[total/voxSize];
       if ( object == NULL ) {
           cout << endl << "VoxMap::loadVolume()";
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

struct Pt3d VoxMap::getMapAt(int x1, int y1, int z1)
{//get the forward destination at grid (x1,y1,z1)
  //int t = z1*planSize+y1*dimX+x1;
  //return (fdMap[t]);

  struct Pt3d pt; 
  int r=getVectorAt(x1,y1,z1,pt);
  if ( r==0 )
     pt.x = pt.y = pt.z = 0; //undefined value
  return pt;
  
}
void VoxMap::setMapAt(int x1, int y1, int z1, struct Pt3d& vv)
{
  int t = z1*planSize+y1*dimX+x1;
  fdMap[t]=vv;
}
void VoxMap::setMapAt(int x1, int y1, int z1, float a, float b, float c)
{
  int t = z1*planSize+y1*dimX+x1;
  fdMap[t].x = a;
  fdMap[t].y = b; 
  fdMap[t].z = c;
}

int VoxMap::getVectorAtAnyPosition(float px, float py, float pz, struct Pt3d & P)
{//Tri-linear interpolation in *fdMap according to the coordinate of smpt
 //This volume have sampZ slices in Z direction
 //This function is revised from  class VolumeData::getValue(fPoint pt)
 //  in VData.C
 // return 1 when P is properly assigned
 // return 0 when P is undefined (but set to be a 0-vector)
    float vx,vy,vz;
    int row,col,slice;
    struct Pt3d *pp, *pp1,*tpp,*tpp1;
    struct Pt3d v1,v2,v3,v4,v5,v6,v7,v8;
    float tx,ty,tz;
    struct Pt3d* volume=fdMap;
    int sampZ = dimZ;

    vx=px;vy=py;vz=pz;

	//cout <<"   BBB:  AnyPos: ("<<vx<<","<<vy<<","<<vz<<")"<<endl;
    if ( ((vx<0)||(vx>dimX-1)) || ((vy<0)||(vy>dimY-1)) 
        || ((vz<0)||(vz>sampZ-1)))
       return  0; //means value undefined

    col = int(vx); //cout << "col=" << col;
    row = int(vy); //cout << ",row=" << row;
    slice=int(vz); //cout << ",slice=" << slice << endl;
   
    pp = volume + planSize*slice + row*dimX + col;
    if (slice == sampZ-1) pp1 = pp;
    else pp1 = pp + planSize;

    v1 = *((struct Pt3d*)pp); //cout << "v1=" << v1.x <<","<< v1.y<<"," <<v1.z <<endl;
    v5 = *((struct Pt3d*)pp1);//cout << "v5=" << v5.x<<","<<v5.y<<","<<v5.z << endl;

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

    v1.x = tz*v5.x+(1-tz)*v1.x;
    v1.y = tz*v5.y+(1-tz)*v1.y;
    v1.z = tz*v5.z+(1-tz)*v1.z;

    P.x=v1.x;  P.y=v1.y;  P.z=v1.z;

    return 1;
}

fPoint VoxMap::triInterpolationXYZ(struct Pt3d* &volume,
	const int sampZ, struct Pt3d ptt)
{
  return triInterpolationXYZ(volume, sampZ, fPoint(ptt));
}

fPoint VoxMap::triInterpolationXYZ(struct Pt3d* &volume, const int sampZ, 
	fPoint smpt)
{//Tri-linear interpolation in *volume according to the coordinate of smpt
 //This volume have sampZ slices in Z direction
 //This function is revised from  class VolumeData::getValue(fPoint pt)
 //  in VData.C
    float vx,vy,vz;
    int row,col,slice;
    struct Pt3d *pp, *pp1,*tpp,*tpp1;
    struct Pt3d v1,v2,v3,v4,v5,v6,v7,v8;
    float tx,ty,tz;

    smpt.getXYZ(vx,vy,vz);
    if ( ((vx<0)||(vx>dimX-1)) || ((vy<0)||(vy>dimY-1)) 
        || ((vz<0)||(vz>sampZ-1)))
       return fPoint(0.0f,0.0f,0.0f,0.0f); //means value undefined

    col = int(vx); //cout << "col=" << col;
    row = int(vy); //cout << ",row=" << row;
    slice=int(vz); //cout << ",slice=" << slice << endl;
   
    pp = volume + planSize*slice + row*dimX + col;
    if (slice == sampZ-1) pp1 = pp;
    else pp1 = pp + planSize;

    v1 = *((struct Pt3d*)pp); //cout << "v1=" << v1.x <<","<< v1.y<<"," <<v1.z <<endl;
    v5 = *((struct Pt3d*)pp1);//cout << "v5=" << v5.x<<","<<v5.y<<","<<v5.z << endl;

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
    	cout << "v2=" << v2.x << ","<< v2.y <<","<<v2.z <<endl;
      	cout << "v6=" << v6.x << ","<< v6.y <<","<<v6.z <<endl;
	cout << "v4=" << v4.x << ","<< v4.y <<","<<v4.z <<endl;
	cout << "v8=" << v8.x << ","<< v8.y <<","<<v8.z <<endl;

	cout << "tx=" << tx << endl;
	cout << "ty=" << ty << endl;
	cout << "tz=" << tz << endl;
	
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

    v1.x = tz*v5.x+(1-tz)*v1.x;
    v1.y = tz*v5.y+(1-tz)*v1.y;
    v1.z = tz*v5.z+(1-tz)*v1.z;

    return fPoint(v1.x,v1.y,v1.z);
}

void  VoxMap::resetVFvolume(struct Pt3d* &volume, int warpZ)
{ // this procedure is for debugging only
   int i, j,k,n;

   k=0;
   for (n = 0; n< warpZ; n++)
   for (j = 0; j< dimY; j++)
   for (i = 0; i< dimX; i++){
	volume[k].x = i;
	volume[k].y = j;
	volume[k].z = n;
	k++;
   }
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
float VoxMap::triInterpolation(float* w, float v0, 
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
int VoxMap::getGridCoverd(int XX,int YY, int ZZ, 
    struct Pt3d *gridPts, float wei[][8], int limit)
{//according to the information recorded at the grid points,
 // the regular grids are transformed to a shape in the atlas 
 // space of 57-slice, this procedure is to find out which
 // of the points in the atlas grid are insided this transformed
 // shape. If A is such a grid point stored in gridPts[i],
 // the wei[i] will be the relative weight of A to the walls of
 // the transformed cube where A locates.
 // 'limit' is the max grid number it can hold in the gridPts[] 
 // array
 // return the actual grid point number, -1 if overflow.

   float minX, minY,minZ, maxX, maxY, maxZ;
   int i; // dltX,dltY, dltZ;
   int x1, y1, z1;
   int ind = planSize*ZZ + dimX * YY + XX;
   struct Pt3d *fdMap1 = fdMap + ind;
   char rr;
   //---need to be checked -- IMPORTANT
   if (XX > dimX-1)
	dltX = 0;
   else 
	dltX = 1;
   if (YY > dimY-1)
	dltY = 0;
   else 
	dltY = dimX;
   if (ZZ > dimZ-1)
	dltZ = 0;
   else 
	dltZ = planSize;

   ver[0] = fdMap1[0];               ver[1] = fdMap1[dltX]; 
   ver[2] = fdMap1[dltY+dltX];       ver[3] = fdMap1[dltY];
   ver[4] = fdMap1[dltZ];            ver[5] = fdMap1[dltZ+dltX]; 
   ver[6] = fdMap1[dltZ+dltY+dltX];  ver[7] = fdMap1[dltZ+dltY];
   
   getBox4Eachface();
   minX=minY=minZ = 65535; maxX=maxY=maxZ=-100;

   for (i=0;i<8;i++){//********
      if (ver[i].x == UNDEFINED_PT)  //undefined; need tobe examined further
	  return 0; //undefined grid;
      if (ver[i].x < minX) minX = ver[i].x;
      if (ver[i].y < minY) minY = ver[i].y;
      if (ver[i].z < minZ) minZ = ver[i].z;
      if (ver[i].x > maxX) maxX = ver[i].x;
      if (ver[i].y > maxY) maxY = ver[i].y;
      if (ver[i].z > maxZ) maxZ = ver[i].z;
   }
   bmin[0]=minX;bmin[1]=minY;bmin[2]=minZ;
   bmax[0]=maxX;bmax[1]=maxY;bmax[2]=maxZ;
   radius = 0;
   for (i=0; i<3;i++){
	radius += fabs(bmax[i] - bmin[i]);
   }
	
   if (minX < 0) minX = 0; if (minY < 0) minY = 0; if (minZ < 0) minZ = 0;
   if (maxX > dimX-1) maxX = dimX-1; //rv
   if (maxY > dimY-1) maxY = dimY-1; //rv
   if (maxZ > dimZ-1) maxZ = dimZ-1;  //atlasZ==57 //rv
   //-------

   i=0;
   for (z1= upper(minZ); z1 < maxZ; z1 ++){
     for (y1= upper(minY); y1 < maxY; y1 ++)
     for (x1= upper(minX); x1 < maxX; x1 ++){
	 //cout << "x1y1z1= ("<< x1 <<"," <<y1<<","<<z1<<")" <<endl;
	rr = inPolyhedron(x1,y1,z1);
	if (rr != 'o'){
          //if ('o'!=inPolyhedron(x1,y1,z1)){//result weight in w
	   if (i>=limit)
	      return -1; // gridPts[] list overflow
	   gridPts[i].x = x1;
	   gridPts[i].y = y1;
	   gridPts[i].z = z1;
	   getWeight(x1,y1,z1,wei[i]); 
	   i++; 	   
       }
     }
   }
   return i;
}
void VoxMap::getWeight(int x1, int y1, int z1, float* ww)
{
   //-----need more job--calculate (x1,y1,z1) distance to ver[8]
   double dis = 0;
   int i;
   for (i=0; i< 8; i++){
	ww[i] = distance(x1,y1,z1,ver[i]);
	dis += ww[i];
   }//cout << "dis = " <<dis<<endl;
   //for (i = 0; i < 8 ; i ++){
	//ww[i] /= dis; cout << "$$$$ ww[" <<i<<"]="<<ww[i]<<endl;
   //}
     
}
double VoxMap::distance(struct Pt3d pp1, struct  Pt3d pp2)
{
   double tt = (pp1.x-pp2.x)*(pp1.x-pp2.x)+(pp1.y-pp2.y)*(pp1.y-pp2.y)+(pp1.z-pp2.z)*(pp1.z-pp2.z);
   return sqrt(tt);
}

double VoxMap::distance(float a, float b, float c, struct Pt3d pp)
{
   double tt = (a-pp.x)*(a-pp.x)+(b-pp.y)*(b-pp.y)+(c-pp.z)*(c-pp.z);
   return sqrt(tt);
}
int VoxMap::upper(float a)
{
  int t = int(a);
  if (a == t) 
      return t;
  else
  if (a>0)
      return t+1;
  else
      return t-1;
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
void VoxMap::randomRayEndpoint( struct Pt3d & endPt, float radius )
{
  double x, y, z, ww, t;

  /* Generate a random point on a sphere of radius 1. */
  /* the sphere is sliced at z, and a random point at angle t
     generated on the circle of intersection. */
  z = 2.0 * drand48() - 1.0;
  t = 2.0 * PI * drand48();
  ww = sqrt( 1 - z*z );
  x = ww * cos( t );
  y = ww * sin( t );
  endPt.x =  radius * x ;
  endPt.y =  radius * y ;
  endPt.z =  radius * z ;
  
  /*printf( "RandomRay returns %6d %6d %6d\n", endPt[X], 
	endPt[Y], endPt[Z] );*/
}
int VoxMap::inBox(struct Pt3d  q)
{
  if( ( bmin[0] <= q.x ) && ( q.x <= bmax[0] ) &&
      ( bmin[1] <= q.y ) && ( q.y <= bmax[1] ) &&
      ( bmin[2] <= q.z ) && ( q.z <= bmax[2] ) )
    return TRUE;
  return FALSE;
}
/*
  This function returns a char:
    'V': the query point a coincides with a Vertex of polyhedron P.
    'E': the query point a is in the relative interior of an Edge of polyhedron P.
    'F': the query point a is in the relative interior of a Face of polyhedron P.
    'i': the query point a is strictly interior to polyhedron P.
    'o': the query point a is strictly exterior to( or outside of) polyhedron P.
*/
char VoxMap::inPolyhedron( float x1, float y1, float z1)
{
   /* r: Ray endpoint. */
   struct Pt3d r, p;  /* p: Intersection point; not used. */
   struct Pt3d q;  /* q: the point to be judged*/
   int f, k = 0, crossings = 0;
   char code = '?';
 
		 //cout <<"ver="<<ver<<endl;
		//cout <<"---before MMM-> in inPolyhedron()" << endl;
		 //for (int i=0; i< 8; i++)
			//outPt("ver[]", ver[i]);
 
   /* If query point is outside bounding box, finished. */
   q.x = x1; q.y = y1; q.z = z1; 
   if ( !inBox(q) )
      return 'o';

   LOOP:
  		//cout <<"MMM-> in inPolyhedron()" << endl;
		// for (i=0; i< 8; i++)
			//outPt("ver[]", ver[i]);
   //cout << endl << endl << endl << "++++++++++++++++++++" <<endl << endl;
   //cout << "Loop ("<<k<<")" << endl;
   while( k++ < 24 ) {
      //cout  << "Loop ("<<k<<")" << endl;
      crossings = 0;
  
      randomRayEndpoint( r, radius ); 
      r = addVec( q, r ); 
      //cout << "Ray endpoint: (" <<r.x<<","<<r.y <<","<<r.z<< ")" <<endl;
  
      for ( f = 0; f < 12; f++ ) {  /* Begin check each face */
         if ( boxTest( f, q, r ) == '0' ) { //outside of the surrounding box
              code = '0';
              //printf("BoxTest = 0 !\n");
         }
         else code = segTriInt( faces[f], q, r, p );
         //cout << "Face ="<<f<<": BoxTest/SegTriInt returns "<<code<<endl<<endl;

         /* If ray is degenerate, then goto outer while to generate another. */
         if ( code == 'p' || code == 'v' || code == 'e' ) {
            cout << "Degenerate ray" << endl;
            goto LOOP;
         }
   
         /* If ray hits face at interior point, increment crossings. */
         else if ( code == 'f' ) {
            crossings++;
            //cout << "crossings = " <<crossings << endl;
         }

         /* If query endpoint q sits on a V/E/F, return that code. */
         else if ( code == 'V' || code == 'E' || code == 'F' )
            return( code );

         /* If ray misses triangle, do nothing. */
         else if ( code == '0' )
            ;

         else 
            fprintf( stderr, "Error, exit(EXIT_FAILURE)\n" ), exit(1);

      } /* End check each face */

      /* No degeneracies encountered: ray is generic, so finished. */
      break;

   } /* End while loop */
 
   //cout << "Crossings ="<< crossings  << endl;
   /* q strictly interior to polyhedron iff an odd number of crossings. */
   if( ( crossings % 2 ) == 1 )
      return   'i';
   else {
     if (k >= 24) 
	cout << " OH!!! test failed" << endl;
     return 'o';
   }
}
struct Pt3d VoxMap::addVec(struct Pt3d a, struct Pt3d b )
{
   struct Pt3d t;
   t.x = a.x + b.x;
   t.y = a.y + b.y;
   t.z = a.z + b.z;

   return t;
}
void VoxMap::getBox4Eachface()
{//get bounding box for each 12 faces
   int   j,k;
   float wss;
   //initialize
   for ( j=0; j < 12; j++ ) {
       Box[j][0].x = ver[ faces[j][0] ].x;
       Box[j][0].y = ver[ faces[j][0] ].y;
       Box[j][0].z = ver[ faces[j][0] ].z;
	//----
       Box[j][1].x = ver[ faces[j][0] ].x;
       Box[j][1].y = ver[ faces[j][0] ].y;
       Box[j][1].z = ver[ faces[j][0] ].z;
	
       for ( k=1; k < 3; k++ ){
          wss = ver[ faces[j][k] ].x;
          if ( wss < Box[j][0].x ) Box[j][0].x = wss;
          if ( wss > Box[j][1].x ) Box[j][1].x = wss;
          wss = ver[ faces[j][k] ].y;
          if ( wss < Box[j][0].y ) Box[j][0].y = wss;
          if ( wss > Box[j][1].y ) Box[j][1].y = wss;
          wss = ver[ faces[j][k] ].z;
          if ( wss < Box[j][0].z ) Box[j][0].z = wss;
          if ( wss > Box[j][1].z ) Box[j][1].z = wss;
      }
    }
}

char VoxMap::boxTest ( int n, struct Pt3d a, struct Pt3d b )
{//test if possible the line(a,b) intersects with faces[n]
  //n: index to the faces, a/b: start/end point of the ray
   //int i; /* Coordinate index */
   float ww;
   //--X--
   ww = Box[ n ][0].x; /* min: lower left */
   if ( (a.x < ww) && (b.x < ww) ) 
      return '0';
   ww = Box[ n ][1].x; /* max: upper right */
   if ( (a.x > ww) && (b.x > ww) ) 
      return '0';
   //--Y--
   ww = Box[ n ][0].y; /* min: lower left */
   if ( (a.y < ww) && (b.y < ww) ) 
      return '0';
   ww = Box[ n ][1].y; /* max: upper right */
   if ( (a.y > ww) && (b.y > ww) ) 
      return '0';
   //--Z--
   ww = Box[ n ][0].z; /* min: lower left */
   if ( (a.z < ww) && (b.z < ww) ) 
      return '0';
   ww = Box[ n ][1].z; /* max: upper right */
   if ( (a.z > ww) && (b.z > ww) ) 
      return '0';
   return '?';
}

/* Assumption: p lies in the plane containing T.
    Returns a char:
     'V': the query point p coincides with a Vertex of triangle T.
     'E': the query point p is in the relative interior of an Edge of triangle T.
     'F': the query point p is in the relative interior of a Face of triangle T.
     '0': the query point p does not intersect (misses) triangle T.
*/

char VoxMap::inTri3D( int T[3], int m, struct Pt3d p )
{
   int i;           /* Index for X,Y,Z           */
   int j;           /* Index for X,Y             */
   //int k;           /* Index for triangle vertex */
   float pp[3];      /* projected p */
   float Tp[3][3];   /* projected T: three new vertices */
   float p1[3];

   p1[0]=p.x; p1[1]=p.y; p1[2]=p.z;
   /* Project out coordinate m in both p and the triangular face */
   j = 0;
   for ( i = 0; i < 3; i++ ) {
     if ( i != m ) {    /* skip largest coordinate */
       pp[j] = p1[i];
       switch (i){
       case 0:
	  Tp[0][j] = ver[T[0]].x;
  	  Tp[1][j] = ver[T[1]].x;
	  Tp[2][j] = ver[T[2]].x;
	  break;
       case 1:
	  Tp[0][j] = ver[T[0]].y;
  	  Tp[1][j] = ver[T[1]].y;
	  Tp[2][j] = ver[T[2]].y;
	  break;
       case 2:
	  Tp[0][j] = ver[T[0]].z;
  	  Tp[1][j] = ver[T[1]].z;
	  Tp[2][j] = ver[T[2]].z;
	  break;
       }
       j++;
     }
   }
   return( inTri2D( Tp, pp ) );
}

int   VoxMap::areaSign( float a[3], float b[3], float c[3] )  
{
    double area2;

    area2 = ( b[0] - a[0] ) * (double)( c[1] - a[1] ) -
            ( c[0] - a[0] ) * (double)( b[1] - a[1] );

    if      ( area2 >  0.0001 ) return  1;
    else if ( area2 < -0.0001 ) return -1;
    else                     return  0;
}                           
 
char VoxMap::inTri2D( float Tp[3][3], float pp[3] )
{
   int area0, area1, area2;

   /* compute three areaSign() values for pp w.r.t. each edge of the face in 2D */
   area0 = areaSign( pp, Tp[0], Tp[1] );
   area1 = areaSign( pp, Tp[1], Tp[2] );
   area2 = areaSign( pp, Tp[2], Tp[0] );
   cout <<"area0="<<area0<<", area1=" <<area1<<", area2="<<area2<< endl;

   if ( ( area0 == 0 ) && ( area1 > 0 ) && ( area2 > 0 ) ||
        ( area1 == 0 ) && ( area0 > 0 ) && ( area2 > 0 ) ||
        ( area2 == 0 ) && ( area0 > 0 ) && ( area1 > 0 ) ) 
     return 'E';

   if ( ( area0 == 0 ) && ( area1 < 0 ) && ( area2 < 0 ) ||
        ( area1 == 0 ) && ( area0 < 0 ) && ( area2 < 0 ) ||
        ( area2 == 0 ) && ( area0 < 0 ) && ( area1 < 0 ) )
     return 'E';                 
   
   if ( ( area0 >  0 ) && ( area1 > 0 ) && ( area2 > 0 ) ||
        ( area0 <  0 ) && ( area1 < 0 ) && ( area2 < 0 ) )
     return 'F';

   if ( ( area0 == 0 ) && ( area1 == 0 ) && ( area2 == 0 ) )
     fprintf( stderr, "Error in InTriD\n" ), exit(EXIT_FAILURE);

   if ( ( area0 == 0 ) && ( area1 == 0 ) ||
        ( area0 == 0 ) && ( area2 == 0 ) ||
        ( area1 == 0 ) && ( area2 == 0 ) )
     return 'V';

   else  
     return '0';  
}
char VoxMap::segTriInt( int T[3], struct Pt3d q, 
	struct Pt3d r, struct Pt3d & p )
{//the face (plan), q: the start point, the point need to be judged if
 // it is inside the polyhedron, r: end point
 //p: the intersection
    char code = '?';
    int m = -1;

    code = segPlaneInt( T, q, r, p, m );
    //cout << "SegPlaneInt code="<< code<<", m="<<m<<"; p=("<<p.x<<")";
    //cout << "," <<p.y <<"," << p.z << endl;

    if ( code == '0') //line segment lies to one side of plane
       return '0';
    else if ( code == 'q') //only point q on plane
       return inTri3D( T, m, q ); //judge if q is inside triangle
    else if ( code == 'r') //only point r on plane
       return inTri3D( T, m, r ); //judge if r is inside triangle
    else if ( code == 'p' )//both point p & r are on plane(seg in plane)
       return 'p'; //inPlane( T, m, q, r, p );
    else if ( code == '1' )//segment intersects with the plane
       return segTriCross( T, q, r );//judge if the intersection inside triangle
    else /* Error */
       return code;
}
/*---------------------------------------------------------------------
Computes N & D and returns index m of largest component.
---------------------------------------------------------------------*/
int VoxMap::planeCoeff( int T[3], struct Pt3d &N, double& D )
{    
    double t;              /* Temp storage */
    double biggest = 0.0;  /* Largest component of normal vector. */
    int m = 0;             /* Index of largest component. */

	//cout << "T[3]=" << T[0]<<","<<T[1]<<","<<T[2]<<endl;//---
	//outPt( "ver[T[0]]",ver[T[0]]); //---XDRXDR
	//for (int i=0; i< 8; i++)
           //outPt("ver[]", ver[i]);

    NormalVec( ver[T[0]], ver[T[1]], ver[T[2]], N );
    /*printf("PlaneCoeff: N=(%lf,%lf,%lf)\n", N[X],N[Y],N[Z]);*/
    D = dot( ver[T[0]], N ); //? what is D?

    /* Find the largest component of N. */
    t = fabs( N.x );
    if ( t > biggest ) {
      biggest = t;
      m = 0;
    }
    t = fabs( N.y );
    if ( t > biggest ) {
      biggest = t;
      m = 1;
    }
    t = fabs( N.z );
    if ( t > biggest ) {
      biggest = t;
      m = 2;
    }
    return m;
}
/*---------------------------------------------------------------------
Returns the dot product of the two input vectors.
---------------------------------------------------------------------*/
double	VoxMap::dot( struct Pt3d a, struct Pt3d b )
{
    return (a.x * b.x + a.y * b.y + a.z * b.z);
}
void VoxMap::outPt(char* ss,struct Pt3d ppp)
{
   cout << ss << "=("<<ppp.x<<","<< ppp.y<<","<<ppp.z<<");"<<endl;
}
/*---------------------------------------------------------------------
 Judge the relation between of a segment of line & a plane
    'p': The segment lies wholly within the plane.
    'q': The q endpoint is on the plane (but not 'p').
    'r': The r endpoint is on the plane (but not 'p').
    '0': The segment lies strictly to one side or the other of the plane.
    '1': The segement intersects the plane, and 'p' does not hold.
---------------------------------------------------------------------*/
char VoxMap::segPlaneInt( int  T[3], struct Pt3d q, 
	struct Pt3d r, struct Pt3d &p, int& m)
{// T: the plane, q: the start point, r: the end point
 // p: the intersection, m: 0/1/2 -> X/Y/Z the biggest compont of the 
 // normal vector of plan T
    struct Pt3d N; 
    double D;
    struct Pt3d rq;
    double num, denom, t;

    //cout << "LLLLLL---->in SegPlaneInt()" << endl;
	 //for (int i=0; i< 8; i++)
           // outPt("ver[]", ver[i]);

    m = planeCoeff( T, N, D ); //outPt("N=",N); cout << "D = "<<D<<endl;
    /*printf("m=%d; plane=(%lf,%lf,%lf,%lf)\n", m, N[X],N[Y],N[Z],D);*/
    num = D - dot( q, N );
    rq = subVec( r, q);
    denom = dot( rq, N );
    //printf("SegPlaneInt: num=%lf, denom=%lf\n", num, denom );

    if ( denom == 0.0 ) {  /* Segment is parallel to plane. */
       if ( num == 0.0 )   /* q is on plane. */
           return 'p';
       else
           return '0';
    }
    else
       t = num / denom;
    /*printf("segPlaneInt: t=%lf \n", t );*/

    //for( i = 0; i < DIM; i++ )
      // p[i] = q[i] + t * ( r[i] - q[i] );
    p.x = q.x + t * (r.x - q.x);
    p.y = q.y + t * (r.y - q.y);
    p.z = q.z + t * (r.z - q.z);

    if ( (0.0 < t) && (t < 1.0) )
         return '1';
    else if ( num == 0.0 )   /* t == 0 */
         return 'q';
    else if ( num == denom ) /* t == 1 */
         return 'r';
    else return '0';
}
/*---------------------------------------------------------------------
Compute the cross product of (b-a)x(c-a) and place into N.
---------------------------------------------------------------------*/
void VoxMap::NormalVec(struct Pt3d a, struct Pt3d b, struct Pt3d c, 
	struct Pt3d &N )
{
    N.x = ( c.z - a.z ) * ( b.y - a.y ) -
           ( b.z - a.z ) * ( c.y - a.y );
    N.y = ( b.z - a.z ) * ( c.x - a.x ) -
           ( b.x - a.x ) * ( c.z - a.z );
    N.z = ( b.x - a.x ) * ( c.y - a.y ) -
           ( b.y - a.y ) * ( c.x - a.x );
    //cout << "======" << endl;
    //outPt("a=",a);outPt("b=",b);outPt("c=",c);outPt("N=",N);
    //cout << "=======" << endl;
}
/*---------------------------------------------------------------------
a - b ==> c.
---------------------------------------------------------------------*/
struct Pt3d VoxMap:: subVec( struct Pt3d a, struct Pt3d b)
{
   struct Pt3d c;

   c.x = a.x - b.x;
   c.y = a.y - b.y;
   c.z = a.z - b.z;
   return c;
}
/*---------------------------------------------------------------------
This procedure deals with the following case:
  the line segment qr intersects with the plane the triangle T[3] locates,
  which means a ray starting from q goes beyond the bounding ball intersects 
  with this plane, now the task is to tell the 3 following cases from
  each other:
    1. segment goes through the edge of the triangle: 'v' degeneracy
    2. segment goes through one of the 3 vertices of the triangle: 'e' degeneracy
    3. segment intersects interior of the triangle: 'f' good, one crossing
    4. segment does not intersect with the triangle: '0' no crossing
The signed volumes of three tetrahedra are computed, determined
by the segment qr, and each edge of the triangle.  
Returns a char:
   'v': the open segment includes a vertex of T.
   'e': the open segment includes a point in the relative interior of an edge
   of T.
   'f': the open segment includes a point in the relative interior of a face
   of T.
   '0': the open segment does not intersect triangle T.
---------------------------------------------------------------------*/
char VoxMap::segTriCross( int T[3], struct Pt3d q, struct Pt3d  r )
{
   int vol0, vol1, vol2;
   
   vol0 = VolumeSign( q, ver[ T[0] ], ver[ T[1] ], r ); 
   vol1 = VolumeSign( q, ver[ T[1] ], ver[ T[2] ], r ); 
   vol2 = VolumeSign( q, ver[ T[2] ], ver[ T[0] ], r );
 
   /*printf( "SegTriCross:  vol0 = %d; vol1 = %d; vol2 = %d\n", 
      vol0, vol1, vol2 ); */
     
   /* Same sign: segment intersects interior of triangle. */
   if ( ( ( vol0 > 0 ) && ( vol1 > 0 ) && ( vol2 > 0 ) ) || 
        ( ( vol0 < 0 ) && ( vol1 < 0 ) && ( vol2 < 0 ) ) )
      return 'f';
   
   /* Opposite sign: no intersection between segment and triangle */
   if ( ( ( vol0 > 0 ) || ( vol1 > 0 ) || ( vol2 > 0 ) ) &&
        ( ( vol0 < 0 ) || ( vol1 < 0 ) || ( vol2 < 0 ) ) )
      return '0';

   else if ( ( vol0 == 0 ) && ( vol1 == 0 ) && ( vol2 == 0 ) )
     fprintf( stderr, "Error 1 in SegTriCross\n" ), exit(EXIT_FAILURE);
   
   /* Two zeros: segment intersects vertex. */
   else if ( ( ( vol0 == 0 ) && ( vol1 == 0 ) ) || 
             ( ( vol0 == 0 ) && ( vol2 == 0 ) ) || 
             ( ( vol1 == 0 ) && ( vol2 == 0 ) ) )
      return 'v';

   /* One zero: segment intersects edge. */
   else if ( ( vol0 == 0 ) || ( vol1 == 0 ) || ( vol2 == 0 ) )
      return 'e';
   
   else
     fprintf( stderr, "Error 2 in SegTriCross\n" ), exit(EXIT_FAILURE);
   return '-'; /*error*/

}

int VoxMap::VolumeSign( struct Pt3d a, struct Pt3d b, 
	struct Pt3d c, struct Pt3d d )
{ 
   double vol;
   double ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz;
   double bxdx, bydy, bzdz, cxdx, cydy, czdz;

   ax = a.x;
   ay = a.y;
   az = a.z;
   bx = b.x;
   by = b.y;
   bz = b.z;
   cx = c.x; 
   cy = c.y;
   cz = c.z;
   dx = d.x;
   dy = d.y;
   dz = d.z;

   bxdx=bx-dx;
   bydy=by-dy;
   bzdz=bz-dz;
   cxdx=cx-dx;
   cydy=cy-dy;
   czdz=cz-dz;
   vol =   (az-dz) * (bxdx*cydy - bydy*cxdx)
         + (ay-dy) * (bzdz*cxdx - bxdx*czdz)
         + (ax-dx) * (bydy*czdz - bzdz*cydy);

   if      ( vol >  0.0001 )   return  1;
   else if ( vol < -0.0001 )  return -1;
   else                    return  0;
}
int VoxMap::testPtinPolyhedron()
{
   float minX, minY,minZ, maxX, maxY, maxZ;
   int i;
   
   ver[0].x =  0; ver[0].y = 10; ver[0].z = 10;
   ver[1].x = 10; ver[1].y = 10; ver[1].z = 10;
   ver[2].x = 10; ver[2].y =  0; ver[2].z = 10;
   ver[3].x =  0; ver[3].y =  0; ver[3].z = 10;
   ver[4].x =  0; ver[4].y = 10; ver[4].z =  0;
   ver[5].x = 10; ver[5].y = 10; ver[5].z =  0;
   ver[6].x = 10; ver[6].y =  0; ver[6].z =  0;
   ver[7].x =  0; ver[7].y =  0; ver[7].z =  0;
  
   //for (i=0; i< 8; i++)
     //outPt("ver[]", ver[i]);
   getBox4Eachface();
   minX=minY=minZ = 65535; maxX=maxY=maxZ=-100;

   cout << " 1111111" << endl;

   for (i=0;i<8;i++){//********
      if (ver[i].x == UNDEFINED_PT)  //undefined; need tobe examined further
	  return 0; //undefined grid;
      if (ver[i].x < minX) minX = ver[i].x;
      if (ver[i].y < minY) minY = ver[i].y;
      if (ver[i].z < minZ) minZ = ver[i].z;
      if (ver[i].x > maxX) maxX = ver[i].x;
      if (ver[i].y > maxY) maxY = ver[i].y;
      if (ver[i].z > maxZ) maxZ = ver[i].z;
   }
   bmin[0]=minX;bmin[1]=minY;bmin[2]=minZ;
   bmax[0]=maxX;bmax[1]=maxY;bmax[2]=maxZ;
	//cout << "minXYZ=("<<minX<<","<<minY<<","<<minZ<<")"<<endl;
	//cout << "maxXYZ=("<<maxX<<","<<maxY<<","<<maxZ<<")"<<endl;
   radius = 0;
   for (i=0; i<3;i++){
	radius += fabs(bmax[i] - bmin[i]);
   }
   cout << " 8888888" << endl;
   //for (i=0; i< 8; i++)
     //outPt("ver[]", ver[i]);

   if (minX < 0) minX = 0; if (minY < 0) minY = 0; if (minZ < 0) minZ = 0;
   if (maxX > 255) maxX = 255;  if (maxY > 255) maxY = 255;
   if (maxZ > dimZ-1) maxX = atlasZ-1;  //atlasZ==57
   //-------
   i=0;
   /*
		cout <<"AAA-> in segTriInt()" << endl;
		 for (i=0; i< 8; i++)
			outPt("ver[]", ver[i]);
		 cout <<"ver="<<ver<<endl;*/
   char cc=inPolyhedron(0,10,10.8); // judge (o,o,o)'s location
   cout << "Result =<" << cc <<">" << endl;
   return i;
}
void VoxMap::createSusumuMapStep1(float rx, float ry, float rz)
{
   Matrix Mt,Ms1,M,M90,Mx,My,Mz,Ms2,Mt1;
   int XX, YY, ZZ,k;
   int mCropped, nCropped, warpZ, volZ;

   cout << "fWarp = " <<fWarp<< endl;
   warpZ = getVolumeZ(fWarp);//get file 'fWarp' Z value for checking purpose
   loadCropInfo(volZ, nCropped,mCropped);
      cout << endl << "--- getForwardMap: Crop Info:--- " << endl;
      cout << "(volZ, nCropped, mCropped)="<<volZ <<","<< nCropped;
      cout << "," << mCropped << endl;
   if (mCropped-nCropped+1 != warpZ){//confirm the data is correct
	cout << endl << "ERROR:VoxMap::getForwardMap()" << endl;
	cout << "     file size does not match expected,warpZ="<<warpZ << endl;
	exit(1);
   } else{
	cout << endl << "warpZ=" <<warpZ<<endl;
 	cout << "*.bottop information correct!!" << endl <<endl;
   }

   Mt.newTranslate(float((1-dimX)/-2.0), float((1-dimY)/-2.0), float(-(-dimZ+0.5))); //shift origin to center
   cout <<"********************************"<<endl;
   //Ms1.newScaleXYZ(1.0,1.0,1.5/0.9766);  
   Ms1.newScaleXYZ(rX,rY,rZ);  //change to physical measurment
   M90.loadIdentity(); M90.rotateX(-90/180.0*PI);//rotateX -90
   rx = rx/180*PI;  ry = ry/180*PI;  rz = rz/180*PI;
   Mx.loadIdentity();  My.loadIdentity();  Mz.loadIdentity();
   Mx.rotateX(rx);  My.rotateY(ry);  Mz.rotateZ(rz);
   //Ms2.newScaleXYZ(1.0,1.0, 0.9766/1.5);
   Ms2.newScaleXYZ(1.0/rX,1.0/rY, 2.0/rZ);
   Mt1.newTranslate(float((dimX-1)/-2.0), float((dimY-1)/-2.0), float(-(dimZ-0.5-nCropped))); //shift origin back and do crop
   //M = Ms1*M90*Mx*My*Mz*Ms2;
   //Matrix MN=Ms1*M90*Mx*My*Mz*Ms2;
   M = Mt*Ms1*M90*Mx*My*Mz*Ms2*Mt1;
   /*
   fPoint kt(27.5f,57.5f,-8.5f), kt1;
   kt1 = kt + fPoint(float((1-dimX)/2.0), float((1-dimY)/2.0), float(-dimZ+0.5));
   kt1.outputs();
   kt1 = kt1*MN;
   kt1.outputs();
   kt1 = kt1 +fPoint(float((dimX-1)/2.0), float((dimY-1)/2.0), float(dimZ-0.5-nCropped));
   cout << "kt1:" <<endl; kt1.outputs();
    cout <<"--kt---" <<endl;
   kt = kt *M;
   cout << "kt:"<<endl;
   kt.outputs();
    //cout <<"--Mt---" <<endl;
   //Mt.outputs();
    //cout <<"--MN---" <<endl;
   //MN.outputs();
    //cout <<"--Mt1---" <<endl;
   //Mt1.outputs();

   exit(0);*/

    //switchXY(vfWarp,warpZ); //Christos's vector field use different X Y pair
   //resetVFvolume(vfWarp, warpZ);// for debugging purpose
   k=0;
   for (ZZ = 0; ZZ < dimZ; ZZ ++){
     for (YY = 0; YY < dimY; YY ++){
       for (XX = 0; XX < dimX; XX ++){
              //double slice: e.g. 56 -> 111 then 112-slice
	  //fPoint pt(fdMap[k].x, fdMap[k].y, fdMap[k].z*2.0); 
	  fPoint pt(XX,YY,ZZ*2.0); //double sliced
	      //pt.setXYZ(127.5f,138.5f,55.5f);cout << endl <<endl << "111: "; pt.outputs();
	      // move Origin to the volume center
	  //***//pt = pt + fPoint(float((1-dimX)/2.0), float((1-dimY)/2.0), float(-dimZ+0.5));//Z is double-sliced
		//cout << endl << endl <<"222: "; pt.outputs();
	      // -90 X-axel, then rx, ry, and rz rotation
	      //reslice: Y:112->256; Z: 256->112
	  pt = pt * M;  //cout << endl <<endl << "333: "; pt.outputs();
	      //move the ORIGIN point back to corner
	  //***//pt = pt + fPoint(float((dimX-1)/2.0), float((dimY-1)/2.0), float(dimZ-0.5-nCropped));
		//cout << endl << endl <<"444: "; pt.outputs();
	      //Z deduct nCropped is because the VF file has less Z slices
	      //warp to 256*256*57 atlas space
	      //---1. vector field or destination definition in the VF file
	      //        if vector field, mathmatical addtion required
	      //	if destination coordinates, then replacement but addtion is to be used
	      //---2. how about the point of undefined value?
	      //Answer: Seems to be destination
 	      //since the real VF file has less Z slices
              //pt = pt + fPoint(0.0f, 0.0f, float(-nCropped));
	      // tri-interpolation
	       pt.getXYZ(fdMap[k].x, fdMap[k].y,fdMap[k].z);
	  k++;
       }
     }
  }
  cout << "Susumu Map Step 1 completed." << endl;
}
struct Pt3d VoxMap::getVectorAt(int col, int row, int slice)
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
   }else{
	/******
	   if (( fdMap[ind].x != 0) &&  ( fdMap[ind].y != 0) && ( fdMap[ind].z != 0)
		&&( fdMap[ind].x != UNDEFINED_PT) &&  ( fdMap[ind].y != UNDEFINED_PT) && ( fdMap[ind].z != UNDEFINED_PT))
		cout <<" \n H@H ("<<fdMap[ind].x<<","<<fdMap[ind].y<<","<<fdMap[ind].z<< ")" <<endl;
	*****/
	return fdMap[ind];
   }
} 

int VoxMap::getVectorAt(int col, int row, int slice, struct Pt3d & P)
{ //To use this function, "this" must not be a displacement field
  // but a vector field of absolute coordinates
  // i.e. 2nd part of reorganization() should be avoided.

   int ind = slice*planSize + row*dimX + col;
   //struct Pt3d P;

   if (((col<0)||(col>dimX-1)) 
       || ((row<0)||(row>dimY-1)) 
       || ((slice<0)||(slice>dimZ-1))){
	  //cout <<"getVecAt(): ("<<col<<","<<row<<","<<slice<<")"<<endl;
	P.x = P.y = P.z = 0;
	return 0; //return P value undefined (meaningless)
   }else{
	/******
	   if (( fdMap[ind].x != 0) &&  ( fdMap[ind].y != 0) && ( fdMap[ind].z != 0)
		&&( fdMap[ind].x != UNDEFINED_PT) &&  ( fdMap[ind].y != UNDEFINED_PT) && ( fdMap[ind].z != UNDEFINED_PT))
		cout <<" \n H@H ("<<fdMap[ind].x<<","<<fdMap[ind].y<<","<<fdMap[ind].z<< ")" <<endl;
	*****/
	P=fdMap[ind];
	return 1; // return P value meaningful
   }
} 


void VoxMap::ChristosExample()
{ // the initial volume is 256 * 256 * 3
  int XX, YY, ZZ;
  int st=70; 
  int t = time(NULL);

  srand48(t); //initialize random generator  

  for (ZZ = 0; ZZ < dimZ; ZZ ++){
     for (YY = 20; YY < dimY-20; YY ++){
       t = ZZ*planSize + YY*dimX +st;
       for (XX = st; XX < st+100; XX ++){
	  fdMap[t].x = XX; 
	  fdMap[t].y = YY + 4.4;
	  t ++;
       }
    }
  }
  t = 0;
  for (ZZ = 0; ZZ < dimZ; ZZ ++){
     for (YY = 0; YY < dimY; YY ++){
       //t = ZZ*planSize + YY*dimX +st;
       for (XX = 0; XX < dimX; XX ++){
	  fdMap[t].x += (drand48()*2-1.0)*4.0;
	  fdMap[t].y += (drand48()*2-1.0)*4.0;
	  t ++;
       }
    }
  }
}

void VoxMap::ChristosExample2()
{ // the initial volume is 256 * 256 * 3
  int XX, YY, ZZ;
  float st=50; 
  float deltx,delty;
  int t = time(NULL);
  float ttx,tty;

  srand48(t); //initialize random generator  
  /*
  for (ZZ = 0; ZZ < dimZ; ZZ ++){
     for (YY = 20; YY < dimY-20; YY ++){
       t = ZZ*planSize + YY*dimX +st;
       for (XX = st; XX < st+30; XX ++){
	  fdMap[t].x = XX; 
	  fdMap[t].y = YY + 2.4;
	  t ++;
       }
    }
  }*/
  for (ZZ = 0; ZZ < dimZ; ZZ ++){
     deltx = 0.065;
     st = 50;
     for (YY = 20; YY < dimY-20; YY ++){
       t = ZZ*planSize + YY*dimX +st;
       ttx = sqrt(deltx*deltx + 1);
       tty = 1.0/ttx;
       ttx = deltx/ttx;
       for (XX = st; XX < st+30; XX ++){
	  fdMap[t].x = XX + ttx*2.4; 
	  fdMap[t].y = YY + tty*2.4;
	  t ++;
       }
       deltx = 1.005*deltx;
       st += deltx;
    }
  }
  for (ZZ = 0; ZZ < dimZ; ZZ ++){
     deltx = 0.42;
     st = 60;
     for (YY = 20; YY < dimY-20; YY ++){
       t = ZZ*planSize + YY*dimX +st;
       ttx = sqrt(deltx*deltx + 1);
       tty = 1.0/ttx;
       ttx = deltx/ttx;
       for (XX = st; XX < st+30; XX ++){
	  fdMap[t].x = XX + ttx*2.4; 
	  fdMap[t].y = YY + tty*2.4;
	  t ++;
       }
       deltx = 1.005*deltx;
       st += deltx;
    }
  }
  //cout <<"------------" << endl;
  for (ZZ = 0; ZZ < dimZ; ZZ ++){
     delty = 0.43;
     st = dimY-50;
     for (XX = 20; XX < dimX - 20; XX ++){
       tty = sqrt(delty*delty + 1);
       ttx = -1.0/tty;
       tty = delty/tty;
       t = ZZ*planSize + int(st)*dimX +XX;

       for (YY = st; YY >st-30; YY --){
          if ((abs(fdMap[t].x - XX) < 0.001) && (abs(fdMap[t].y - YY)<0.001)){
              fdMap[t].x = XX + ttx*2.4; 
	      fdMap[t].y = YY + tty*2.4;
	  } else {
	      fdMap[t].x = (XX + ttx*2.4+fdMap[t].x)/2.0;
	      fdMap[t].y = (YY + tty*2.4+fdMap[t].y)/2.0;
	  }
	  t -= dimX;
       }
       delty = 1.001*delty;
       st -= delty;
	//cout <<"---delty = " << delty<<"---st="<<st<<endl;
    }
  }
  t = 0;
  for (ZZ = 0; ZZ < dimZ; ZZ ++){
     for (YY = 0; YY < dimY; YY ++){
       //t = ZZ*planSize + YY*dimX +st;
       for (XX = 0; XX < dimX; XX ++){
	  fdMap[t].x += (drand48()*2-1.0)*1.5;
	  fdMap[t].y += (drand48()*2-1.0)*1.5;
	  t ++;
       }
    }
  }
}

void VoxMap::averageAndSave(int nSlc, char * toStore)
{// this function  takes fdMap[] as a stack of a number of volumes
 // nSlc is the number of slice of each volume, therefore dimZ = nSlc * n
 // the task of this function is to average the volume and get a resulting average voume
 // the result is of (dimX x dimY x nSlc) and will be store to  file "toStore"

  int i, XX, YY, ZZ, k, pk, volCt;
  int interval = nSlc *planSize;
  struct Pt3d *resVol= NULL;
  float *deviation=NULL;
  float  tx,ty,tz;

  volCt = (float)dimZ / (float)nSlc; //get the number of volume this stack has
  if ( volCt != (float)dimZ / (float)nSlc){
    cout <<"ERROR: the volume stack is not completed: Stack slice number= ";
    cout <<dimZ<<", average Slc = "<<nSlc<<endl;    
    return;
  } else{
    cout << " The stack contains " << volCt << " volumes" << endl;
  }
  resVol = new struct Pt3d [interval];
  if (resVol == NULL){
     cout <<" VoxMap::averageEveryNslice():: ERROR : memory out1"<<endl;
     return;
  } else{
     cout << " VoxMap::averageEveryNslice(): got memory1 for result"<<endl;
  }
  deviation = new float [interval];
  if (deviation == NULL){
     cout <<" VoxMap::averageEveryNslice():: ERROR : memory out2"<<endl;
     return;
  } else{
     cout << " VoxMap::averageEveryNslice(): got memory2 for result"<<endl;
  }
  //saveVecVolume(resVol,nSlc,toStore);

  k =0;
  for (ZZ = 0; ZZ < nSlc; ZZ ++){
     for (YY = 0; YY < dimY; YY ++){
       for (XX = 0; XX < dimX; XX ++){
	 pk = k;
	 tx = ty = tz = 0.0;
         for (i = 0; i < volCt ; i++, pk += interval ){
	   tx += fdMap[pk].x;
	   ty += fdMap[pk].y;
	   tz += fdMap[pk].z;
	 }
	 resVol[k].x = tx /(float)volCt;
 	 resVol[k].y = ty /(float)volCt;
 	 resVol[k].z = tz /(float)volCt;

	 Vector meanVec(resVol[k].x,resVol[k].y,resVol[k].z);
	 tx = 0; pk = k;
         for (i = 0; i < volCt ; i++, pk += interval ){
	   Vector nnv(fdMap[pk].x,fdMap[pk].y,fdMap[pk].z);
	   float ag = nnv.angle(meanVec);
	   if (ag == -1)
		ag = PI/2.0;
	   else 
	   if (ag > PI/2.0)
		ag = PI - ag;
	   tx += ag*ag;
	 }
         if (volCt > 1)
	    deviation[k] = sqrt(tx /(float)(volCt-1));
	 else 
	    deviation[k] = 0;
	 k++;
       }
     }
  }
  //---save-deviation----
    FILE * fp = NULL;
    cout << "VoxMap:  saveDeviation"<<endl;
    fp = fopen("deviationAngle.floatVolume", "wb");
    if(fp==NULL) {
      cout << endl <<"ERROR: VoxMap:  saveDeviation :file open failed" << endl;
      cout << "file open failed: <deviationAngle.floatVolume> " << endl;
    }else{
      int n=fwrite(deviation, sizeof(float), interval,fp); 
      if (n != interval) {
        cout << endl << " ERROR: VoxMap: saveDeviation :writing data failed" << endl;       
        cout << "file saved size = " << n << ", while " << interval << " expected." << endl;
      }else 
      cout << "file saved correctly" << endl;
   }
   fclose (fp);
   cout << endl << "file <deviationAngle.floatVolume> saved " <<endl;

  //---------------------
  saveVecVolume(resVol,nSlc,toStore);
  delete []resVol;
}

void VoxMap::colormapRGB3Volume()
{
   int XX, YY, ZZ, k;
   ImgVol RR(dimZ), GG(dimZ), BB(dimZ);
   char rr,gg,bb, tname[200];

   k =0;
   for (ZZ = 0; ZZ < dimZ; ZZ ++){
     for (YY = 0; YY < dimY; YY ++){
       for (XX = 0; XX < dimX; XX ++){
	   rr = char (fabs(fdMap[k].x)*255.0 + 0.5);
	   gg = char (fabs(fdMap[k].y)*255.0 + 0.5);
	   bb = char (fabs(fdMap[k].z)*255.0 + 0.5);
	   RR.setVoxel(k,rr);
	   GG.setVoxel(k,gg);
	   BB.setVoxel(k,bb);
	   k++; 
       }
     }
   }
   /*	
   char firstsection[200];
   int i = 0;
   while (fWarp[i]!='0' && fWarp[i]!='.'){
	firstsection[i]=fWarp[i];
	i++;
   }
   firstsection[i] = '\0';
   */
   //---begin saving the 3 components---
   strcpy(tname,fWarp);
   strcat(tname,".xVolR");
   cout << "\n\nSaving X-component to R-volume:<"<<tname<<">...." <<endl;
   RR.saveVolume(tname);

   strcpy(tname,fWarp);
   strcat(tname,".yVolG");
   cout << "\n\nSaving X-component to R-volume:<"<<tname<<">...." <<endl;
   GG.saveVolume(tname);

   strcpy(tname,fWarp);
   strcat(tname,".zVolB");
   cout << "\n\nSaving X-component to R-volume:<"<<tname<<">...." <<endl;
   BB.saveVolume(tname);

   cout <<"\nSaving completed."<<endl<<endl;
}
void VoxMap::setBlendMapAt(int x1,int y1, int z1, float& a, float& b, float& c)
{//average & normalize
  double a1,b1,c1, t;
  struct Pt3d ppp;

  ppp= getMapAt(x1,y1,z1);
  a1 = (ppp.x + a)/2.0 - x1;
  b1 = (ppp.y + b)/2.0 - y1;
  c1 = (ppp.z + c)/2.0 - z1;

  t = sqrt (a1*a1 + b1*b1 + c1*c1);

  if (t< 0.00001){
	a = x1; b = y1; c = z1;
  }else {
	a = a1/t + x1;
	b = b1/t + y1;
	c = c1/t + z1;
  }
  setMapAt(x1,y1,z1,a,b,c);
}


void VoxMap::setAvgBlendMapAt(int x1,int y1, int z1, float a, float b, float c)
{// simply average
  float a1,b1,c1;
  struct Pt3d ppp;

  ppp= getMapAt(x1,y1,z1);

  if ((ppp.x != x1) ||(ppp.y != y1) || (ppp.z != z1)){
    Vector V1(ppp.x-x1,ppp.y-y1,ppp.z-z1);
    Vector V2(a-x1,b-y1,c-z1);
    cout <<"\n     ##########v1 v2######################x1y1z1="<<x1<<","<<y1<<","<<z1<<endl;
      V1.outputs();
      V2.outputs();

    V1.norm(); //normalize
    V2.norm(); 
      V1.outputs();
      V2.outputs();
   
    Vector V3=(V1+V2)*0.5; //average
      cout <<" V3=";
      V3.outputs();
    V3.norm(); //normalize
      V3.outputs();
    V3.getXYZ(a1,b1,c1);
    setMapAt(x1,y1,z1,a1+x1,b1+y1,c1+z1);
  }else{
    setMapAt(x1,y1,z1,a,b,c);
  }
}


struct Pt3d VoxMap::getUnitVec(float x1, float y1, float z1,
	float x2, float y2, float z2)
{
   float tx,ty,tz, tt;
   struct Pt3d  pp;

   tx = x2-x1;
   ty = y2-y1;
   tz = z2-z1;

   tt = sqrt(tx*tx+ty*ty+tz*tz);
      
   if (tt < 0.00001)
	pp.x = pp.y = pp.z = 0;
   else {
	pp.x = tx/tt;
	pp.y = ty/tt;
	pp.z = tz/tt;
   }
   return pp;
}

void VoxMap::VF3DFiber(int Z1, float& x1, float& y1, float& z1, 
        float& deltx, float& delty, float& deltz, 
	float fx, float fy, float fz, ImgVol& NV, FILE* ff)
{/*&&&&*/
  int XX,YY,ZZ;
  struct Pt3d txyz[2];
  float MM = 2.5;
  int sz=sizeof(struct Pt3d);
  
  XX=x1; YY=y1; ZZ=z1;
  txyz[0].x = x1;  txyz[0].y = y1; txyz[0].z = z1;
  txyz[1].x=txyz[1].y =txyz[1].z = -5000;//flag of End-Of-Segment
  
  //txyz[1] = getUnitVec(XX,YY,ZZ,deltx*MM+XX, delty*MM+YY, deltz*MM+ZZ );
  //txyz[1].x += txyz[0].x; txyz[1].y += txyz[0].y; txyz[1].z += txyz[0].z;
  setAvgBlendMapAt(XX,YY,ZZ, deltx*MM+XX, delty*MM+YY, deltz*MM+ZZ);
  //setBlendMapAt(XX,YY,ZZ, txyz[1].x, txyz[1].y, txyz[1].z);
  if (ff!=NULL){
    fwrite(txyz, sz,1,ff);
  }
  NV.setVoxel(XX,YY,ZZ,235);
  while (ZZ < Z1){
    while (((int)x1 == XX) &&((int)y1 == YY) && ((int)z1 == ZZ)){
       x1 += deltx; y1 += delty; z1 += deltz;
    }
    deltx *= fx;
    delty *= fy;
    deltz *= fz;     
    XX = (int)x1; YY = int(y1); ZZ = int(z1);
    txyz[0].x = x1;  txyz[0].y = y1; txyz[0].z = z1;
      //txyz[0].x = XX;  txyz[0].y = YY; txyz[0].z = ZZ;
      //txyz[1] = getUnitVec(XX,YY,ZZ,deltx*MM+XX, delty*MM+YY, deltz*MM+ZZ );
      //txyz[1].x += txyz[0].x; txyz[1].y += txyz[0].y; txyz[1].z += txyz[0].z;
    setAvgBlendMapAt(XX,YY,ZZ, deltx*MM+XX, delty*MM+YY, deltz*MM+ZZ);
    //setBlendMapAt(XX,YY,ZZ, txyz[1].x, txyz[1].y, txyz[1].z);
    if (ff!=NULL){
	fwrite(txyz, sz,1,ff);
    }
    NV.setVoxel(XX,YY,ZZ,235);
    cout << "(XX,YY,ZZ)=("<<XX<<","<<YY<<","<<ZZ<<")"<<endl;
    cout << "    (vX,vY,vZ)=("<<deltx*MM<<","<<delty*MM<<","<<deltz*MM<<")"<<endl;
    //cout << "      -(vX,vY,vZ)=("<<txyz[1].x-XX<<","<<txyz[1].y-YY<<","<<txyz[1].z-ZZ<<")"<<endl;
    
    //deltx *= fx;
    //delty *= fy;
    //deltz *= fz;     
  }

  if (ff!=NULL){ //write a End-Of-Segment
     fwrite(txyz+1, sz,1,ff);
  }



}

void VoxMap::VF3Dsimulation(FILE * binDotFile)
{ // the initial volume is 256 * 256 * 124
  //int XX, YY, ZZ;
  float x1,y1,z1;
  //float st=50; 
  float deltx,delty, deltz;
  int t = time(NULL);
  //float ttx,tty,ttz; //starting point=(128,136,75)
  ImgVol NV(dimZ);

  srand48(t); //initialize random generator  
  /*
  for (ZZ = 0; ZZ < dimZ; ZZ ++){
     for (YY = 20; YY < dimY-20; YY ++){
       t = ZZ*planSize + YY*dimX +st;
       for (XX = st; XX < st+30; XX ++){
	  fdMap[t].x = XX; 
	  fdMap[t].y = YY + 2.4;
	  t ++;
       }
    }
  }*/
  //----bundle 1------
  
  x1=131; y1=147; z1=66;
  deltx = 0.35;
  delty = 0.018;
  deltz = 0.082;
  VF3DFiber(73, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.1, NV, binDotFile);
  VF3DFiber(101, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.05,NV, binDotFile);  
  deltx = -0.35;
  delty = 0.018;
  deltz = 0.082;
  x1=130; y1=147; z1=66;
  VF3DFiber(75, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.12, NV, binDotFile);
  VF3DFiber(104, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.06, NV, binDotFile);
  //--------
   
  x1=131; y1=146; z1=66;
  deltx = 0.35;
  delty = 0.018;
  deltz = 0.0831;
  VF3DFiber(74, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.1, NV, binDotFile);
  VF3DFiber(104, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.055,NV, binDotFile);
  deltx = -0.35;
  delty = 0.017;
  deltz = 0.0842;
  x1=130; y1=146; z1=66;
  VF3DFiber(75, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.12, NV, binDotFile);
  VF3DFiber(100, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.05, NV, binDotFile);
  //--------
  x1=131; y1=148; z1=65;
  deltx = 0.35;
  delty = 0.018;
  deltz = 0.08;
  VF3DFiber(75, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.092, NV, binDotFile);
  VF3DFiber(99, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.05,NV, binDotFile);
  deltx = -0.35;
  delty = 0.017;
  deltz = 0.083;
  x1=130; y1=148; z1=65;
  VF3DFiber(75, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.12, NV, binDotFile);
  VF3DFiber(102, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.08, NV, binDotFile);
  
  //--------
  x1=131; y1=147; z1=65; //y1=142
  deltx = 0.35;
  delty = 0.0182;
  deltz = 0.0833;
  VF3DFiber(75, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.12, NV, binDotFile);
  VF3DFiber(103, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.05,NV, binDotFile);
  deltx = -0.35;
  delty = 0.0172;
  deltz = 0.0825;
  x1=130; y1=147; z1=65;
  VF3DFiber(75, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.12, NV, binDotFile);
  VF3DFiber(105, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.06, NV, binDotFile);
  
  //--------
  x1=131; y1=148; z1=64;
  deltx = 0.35;
  delty = 0.0177;
  deltz = 0.081;
  VF3DFiber(75, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.11, NV, binDotFile);
  deltz = 0.12;
  VF3DFiber(82, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.04, NV, binDotFile);
  deltz = 0.18;
  VF3DFiber(93, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.12,NV, binDotFile);
  deltx = -0.35;
  delty = 0.018;
  deltz = 0.082;
  x1=130; y1=148; z1=64;
  VF3DFiber(84, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.11, NV, binDotFile);
  deltz = 0.18;
  VF3DFiber(95, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.11, NV, binDotFile);
   
  //--------
  x1=131; y1=147; z1=65;
  deltx = 0.35;
  delty = 0.0186;
  deltz = 0.13;
  VF3DFiber(75, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.12, NV, binDotFile);
  deltz = 0.18;
  VF3DFiber(99, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.05,NV, binDotFile);
  deltx = -0.35;
  delty = 0.021;
  deltz = 0.1;
  x1=130; y1=147; z1=65;
  VF3DFiber(84, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.12, NV, binDotFile);
  deltz = 0.18;
  VF3DFiber(99, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.05, NV, binDotFile);
  
  //--------
  x1=131; y1=148; z1=66;
  deltx = 0.35;
  delty = 0.019;
  deltz = 0.11;
  VF3DFiber(82, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.112, NV, binDotFile);
  deltz = 0.18;
  VF3DFiber(100, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.05,NV, binDotFile);
  deltx = -0.35;
  delty = 0.0192;
  deltz = 0.082;
  x1=130; y1=148; z1=66;
  VF3DFiber(76, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.112, NV, binDotFile);
  deltz = 0.1;
  VF3DFiber(90, x1, y1, z1, deltx, delty, deltz, 1.005, 1.001,1.04, NV, binDotFile);
  
  //-----bundle 2------
  //----(1)-----
  x1=131; y1=102; z1=73;
  deltx = 0.373;
  delty = -0.024;
  deltz = 0.145;
  VF3DFiber(100, x1, y1, z1, deltx, delty, deltz, 1.007, 1.003,1.10, NV, binDotFile);
  deltx = -0.368;
  delty = -0.025;
  deltz = 0.145;
  x1=130; y1=102; z1=73;
  VF3DFiber(86, x1, y1, z1, deltx, delty, deltz, 1.006, 1.003,1.11, NV, binDotFile);
  VF3DFiber(99, x1, y1, z1, deltx, delty, deltz, 1.008, 1.003,1.14, NV, binDotFile);
  //----(2)----
  x1=131; y1=103; z1=73;
  deltx = 0.373;
  delty = -0.024;
  deltz = 0.145;
  VF3DFiber(100, x1, y1, z1, deltx, delty, deltz, 1.007, 1.003,1.10, NV, binDotFile);
  deltx = -0.368;
  delty = -0.025;
  deltz = 0.145;
  x1=130; y1=103; z1=73;
  VF3DFiber(86, x1, y1, z1, deltx, delty, deltz, 1.006, 1.003,1.11, NV, binDotFile);
  VF3DFiber(99, x1, y1, z1, deltx, delty, deltz, 1.008, 1.003,1.14, NV, binDotFile);
  //----(6)----
  x1=131; y1=102; z1=74;
  deltx = 0.373;
  delty = -0.024;
  deltz = 0.145;
  VF3DFiber(100, x1, y1, z1, deltx, delty, deltz, 1.007, 1.003,1.10, NV, binDotFile);
  deltx = -0.368;
  delty = -0.025;
  deltz = 0.145;
  x1=130; y1=102; z1=74;
  VF3DFiber(86, x1, y1, z1, deltx, delty, deltz, 1.006, 1.003,1.11, NV, binDotFile);
  VF3DFiber(99, x1, y1, z1, deltx, delty, deltz, 1.008, 1.003,1.14, NV, binDotFile);
        
  //---(3)---
  x1=131; y1=104; z1=73;
  deltx = 0.373;
  delty = -0.024;
  deltz = 0.145;
  VF3DFiber(78, x1, y1, z1, deltx, delty, deltz, 1.007, 1.003,1.1, NV, binDotFile);
  deltz = 0.115;
  VF3DFiber(88, x1, y1, z1, deltx, delty, deltz, 1.007, 1.003,1.0125, NV, binDotFile);
  deltx = -0.368;
  delty = -0.025;
  deltz = 0.145;
  x1=130; y1=104; z1=73;
  VF3DFiber(80, x1, y1, z1, deltx, delty, deltz, 1.006, 1.003,1.10, NV, binDotFile);
  deltz = 0.085;
  VF3DFiber(88, x1, y1, z1, deltx, delty, deltz, 1.008, 1.003,1.04, NV, binDotFile);
  //---(4)---
  x1=131; y1=104; z1=74;
  deltx = 0.373;
  delty = -0.024;
  deltz = 0.145;
  VF3DFiber(78, x1, y1, z1, deltx, delty, deltz, 1.007, 1.003,1.1, NV, binDotFile);
  deltz = 0.115;
  VF3DFiber(88, x1, y1, z1, deltx, delty, deltz, 1.007, 1.003,1.0125, NV, binDotFile);
  deltx = -0.368;
  delty = -0.025;
  deltz = 0.145;
  x1=130; y1=104; z1=74;
  VF3DFiber(80, x1, y1, z1, deltx, delty, deltz, 1.006, 1.003,1.10, NV, binDotFile);
  deltz = 0.085;
  VF3DFiber(88, x1, y1, z1, deltx, delty, deltz, 1.008, 1.003,1.04, NV, binDotFile);
  //---(5)---
  x1=131; y1=103; z1=74;
  deltx = 0.373;
  delty = -0.024;
  deltz = 0.145;
  VF3DFiber(78, x1, y1, z1, deltx, delty, deltz, 1.007, 1.003,1.1, NV, binDotFile);
  deltz = 0.115;
  VF3DFiber(88, x1, y1, z1, deltx, delty, deltz, 1.007, 1.003,1.0125, NV, binDotFile);
  deltx = -0.368;
  delty = -0.025;
  deltz = 0.145;
  x1=130; y1=103; z1=74;
  VF3DFiber(80, x1, y1, z1, deltx, delty, deltz, 1.006, 1.003,1.10, NV, binDotFile);
  deltz = 0.085;
  VF3DFiber(88, x1, y1, z1, deltx, delty, deltz, 1.008, 1.003,1.04, NV, binDotFile);

  //----(7)-----
  x1=131; y1=103; z1=75;
  deltx = 0.373;
  delty = -0.024;
  deltz = 0.145;
  VF3DFiber(78, x1, y1, z1, deltx, delty, deltz, 1.007, 1.003,1.1, NV, binDotFile);
  deltz = 0.115;
  VF3DFiber(88, x1, y1, z1, deltx, delty, deltz, 1.007, 1.003,1.0125, NV, binDotFile);
  deltx = -0.368;
  delty = -0.025;
  deltz = 0.145;
  x1=130; y1=103; z1=74;
  VF3DFiber(80, x1, y1, z1, deltx, delty, deltz, 1.006, 1.003,1.10, NV, binDotFile);
  deltz = 0.085;
  VF3DFiber(88, x1, y1, z1, deltx, delty, deltz, 1.008, 1.003,1.04, NV, binDotFile);
  //----(8)-----
  x1=131; y1=104; z1=75;
  deltx = 0.373;
  delty = -0.024;
  deltz = 0.145;
  VF3DFiber(78, x1, y1, z1, deltx, delty, deltz, 1.007, 1.003,1.1, NV, binDotFile);
  deltz = 0.115;
  VF3DFiber(88, x1, y1, z1, deltx, delty, deltz, 1.007, 1.003,1.0125, NV, binDotFile);
  deltx = -0.368;
  delty = -0.025;
  deltz = 0.145;
  x1=130; y1=103; z1=74;
  VF3DFiber(80, x1, y1, z1, deltx, delty, deltz, 1.006, 1.003,1.10, NV, binDotFile);
  deltz = 0.085;
  VF3DFiber(88, x1, y1, z1, deltx, delty, deltz, 1.008, 1.003,1.04, NV, binDotFile);
    

  NV.saveVolume("NV.img");
}
void VoxMap::recoverFiber(struct Pt3d * dList,int dNum)
{
  int i,j;
  float a,b,c,t,dis,delt,dltX,dltY,dltZ;
  float nmx,nmy,nmz;
  int XX,YY,ZZ;
  /*
  for (i=0; i< 10; i++)
    cout <<"dot(XYZ-"<<i<<")=("<<dList[i].x<<","<<dList[i].y<<","<<dList[i].z<<")"<<endl;
  cout <<"dot(XYZ-"<<dNum-1<<")=("<<dList[dNum-1].x<<","<<dList[dNum-1].y<<","<<dList[dNum-1].z<<")"<<endl;

  */

  
  for (i=0; i< dNum-1; i++){
    if ((dList[i+1].x==-5000)
	&&(dList[i+1].y==-5000)
	&& (dList[i+1].z==-5000)){
      i ++;
      continue;

    }
    dis = distance(dList[i].x,dList[i].y,dList[i].z,dList[i+1]);
      cout << "\n-----\ni="<<i<<"; dis="<< dis<<endl;
    Vector v(dList[i+1].x-dList[i].x, dList[i+1].y-dList[i].y, dList[i+1].z-dList[i].z );
       v.outputs();
    v.norm();
       v.outputs();
    v.getXYZ(nmx,nmy,nmz);
    v.getXYZ(a,b,c);
    a = fabs(a);b=fabs(b);c=fabs(c);
    if (a>b){//sort a,b,c
      t=a;a=b;b=t;
    }
    if (a>c){
      t=a;a=c;c=t;
    }
    if (b>c){
      t=b;b=c;c=t;
    }//now a<b<c
    if (c<0.0001){ //c==0, do nothing
      //.......

    }else{
      //b = a + v*t
         cout<<" # c="<<c<<"; ";
      c = 1.0/c;
      v.scale(c,c,c); 
         v.outputs();
      v.getXYZ(dltX,dltY,dltZ);
      delt = sqrt(a*a+b*b+c*c);
      a=dList[i].x; b=dList[i].y; c=dList[i].z;
      for (j=0; j<dis;j+=delt){
	XX = (int)a; YY = (int)b; ZZ = (int)c;
	cout << "\n nmx,nmy,nmz="<<nmx<<","<<nmy<<","<<nmz<<"; ";
	setAvgBlendMapAt(XX,YY,ZZ, nmx+XX, nmy+YY, nmz+ZZ);
	  struct Pt3d tt1 = getMapAt(XX,YY,ZZ);
	  cout <<"\nXX,YY,ZZ="<<XX<<","<<YY<<","<<ZZ;
	  cout <<" ; after assign: "<<tt1.x<<","<<tt1.y<<","<<tt1.z<<"; ";
	a += dltX; b += dltY; c += dltZ;
      }
    }
  }
}

void VoxMap::getXYZdim(int& a, int& b, int& c)
{
  a=dimX;b=dimY;c=dimZ;

}

void VoxMap::DTIwarpExp1DisplacementFieldCircular()
{
  int XX,YY, ZZ;
  float a, b, c;
  struct Pt3d pp; //p1,p2;
  double d,factor;

  cout <<" The displacement field is created based on a factor of :";
  cout <<"\n    (distance^3)/factor"<<endl;
  cout <<" Smaller results in smaller vortex, larger result in larger vortex\n"; 
  cout <<" input a factor(default = 25):";
  cin >>factor;
  cout <<endl<<"You have inputted: "<<factor<<endl<<endl;

  //p1.x = 128; p1.y = 128; p1.z =0;
  //p2.x = 128; p2.y = 128; p2.z =1;
  pp.x = 128; pp.y = 128; //pp.z =1;

   for (ZZ = 0; ZZ < dimZ; ZZ ++){
     pp.z=ZZ;
     for (YY = 0; YY < dimY; YY ++){
       for (XX = 0; XX < dimX; XX ++){
	  //d = distance (XX,YY,ZZ, (ZZ==0)?p1:p2)/factor;
	  d = distance (XX,YY,ZZ, pp)/factor;
	  d = d*d*d; 
	  Vector v (fPoint(128,128,ZZ),fPoint(XX,YY,ZZ));
	  if ( XX !=128 || YY != 128 ){
	    v.norm();
            v.getXYZ(a,b,c);
	    a *= d; b *= d;
	    setMapAt(XX,YY,ZZ,XX-b,YY+a,ZZ);
	  } else
	    setMapAt(XX,YY,ZZ,XX,YY,ZZ);
       }
    }
  }
}
void VoxMap::DTIwarpExp1DisplacementFieldRadiational()
{
  int XX,YY, ZZ,i;
  float a, b, c;
  struct Pt3d p1,p2;
  double d=1.0; //scale factor
  double t, t1;

  p1.x = 128; p1.y = 128; p1.z =0;
  p2.x = 128; p2.y = 128; p2.z =1;

   for (ZZ = 0; ZZ < dimZ; ZZ ++){
     for (YY = 0; YY < dimY; YY ++){
       for (XX = 0; XX < dimX; XX ++){
	  Vector v (fPoint(128,128,ZZ),fPoint(XX,YY,ZZ));
	  v.norm();
          v.getXYZ(a,b,c);
	  t = distance (XX,YY,ZZ, (ZZ==0)?p1:p2);
          d = 1.0;
          for (i = 0 ; i < t; i++){
	    //if (t < 10)
	      //d = d*1.01;
	    //else{
	      t1 = (t+8) / 10000;
	      d = d * (1.00006 +t1);
	    //}
          }
	  d = d*d;
	  a *= d; b *= d;
	  setMapAt(XX,YY,ZZ,XX+a,YY+b,ZZ);
       }
    }
  }
}

void VoxMap::DTIwarpExperiment1FiberBunch()
{
  int XX,YY;

   for (XX = 0; XX < dimX; XX ++){
       //Z-slice 1
        setMapAt(XX,124,0,XX+1,124,0);
        setMapAt(XX,125,0,XX+1,125,0);
        setMapAt(XX,126,0,XX+1,126,0);
        setMapAt(XX,127,0,XX+1,127,0);
        setMapAt(XX,128,0,XX+1,128,0);
        setMapAt(XX,129,0,XX+1,129,0);
        setMapAt(XX,130,0,XX+1,130,0);
        setMapAt(XX,131,0,XX+1,131,0);
       //Z-slice 2
        setMapAt(XX,124,1,XX+1,124,1);
        setMapAt(XX,125,1,XX+1,125,1);
        setMapAt(XX,126,1,XX+1,126,1);
        setMapAt(XX,127,1,XX+1,127,1);
        setMapAt(XX,128,1,XX+1,128,1);
        setMapAt(XX,129,1,XX+1,129,1);
        setMapAt(XX,130,1,XX+1,130,1);
        setMapAt(XX,131,1,XX+1,131,1);
   }
   for (YY = 0; YY < dimY; YY ++){
       //Z-slice 1
        setMapAt(124,YY,0,124,YY+1,0);
        setMapAt(125,YY,0,125,YY+1,0);
        setMapAt(126,YY,0,126,YY+1,0);
        setMapAt(127,YY,0,127,YY+1,0);
        setMapAt(128,YY,0,128,YY+1,0);
        setMapAt(129,YY,0,129,YY+1,0);
        setMapAt(130,YY,0,130,YY+1,0);
        setMapAt(131,YY,0,131,YY+1,0);
       //Z-slice 2
        setMapAt(124,YY,1,124,YY+1,1);
        setMapAt(125,YY,1,125,YY+1,1);
        setMapAt(126,YY,1,126,YY+1,1);
        setMapAt(127,YY,1,127,YY+1,1);
        setMapAt(128,YY,1,128,YY+1,1);
        setMapAt(129,YY,1,129,YY+1,1);
        setMapAt(130,YY,1,130,YY+1,1);
        setMapAt(131,YY,1,131,YY+1,1);
   }
   for (YY = 205; YY < 216; YY ++){
     for (XX = 0; XX < dimX; XX ++){
       //Z-slice 1
        setMapAt(XX,YY,0,XX+1,YY,0);
       //Z-slice 2
        setMapAt(XX,YY,1,XX+1,YY,1);
     }
   }

}
void VoxMap::DTIwarpExperiment2PD1FiberBunch()
{ // simulating fibers bundles for consideration of 2PD
  // add one more fibers running Z direction
  // in addition to the 3 bundles of DTIwarpExperiment1FiberBunch()
  int XX,YY,ZZ;

   for (XX = 0; XX < dimX; XX ++){
     for (ZZ= 0; ZZ<10; ZZ ++){
        setMapAt(XX,124,ZZ,XX+1,124,ZZ);
        setMapAt(XX,125,ZZ,XX+1,125,ZZ);
        setMapAt(XX,126,ZZ,XX+1,126,ZZ);
        setMapAt(XX,127,ZZ,XX+1,127,ZZ);
        setMapAt(XX,128,ZZ,XX+1,128,ZZ);
        setMapAt(XX,129,ZZ,XX+1,129,ZZ);
        setMapAt(XX,130,ZZ,XX+1,130,ZZ);
        setMapAt(XX,131,ZZ,XX+1,131,ZZ);
     }
   }
   for (YY = 0; YY < dimY; YY ++){
     for (ZZ= 0; ZZ<10; ZZ ++){
        setMapAt(124,YY,ZZ,124,YY+1,ZZ);
        setMapAt(125,YY,ZZ,125,YY+1,ZZ);
        setMapAt(126,YY,ZZ,126,YY+1,ZZ);
        setMapAt(127,YY,ZZ,127,YY+1,ZZ);
        setMapAt(128,YY,ZZ,128,YY+1,ZZ);
        setMapAt(129,YY,ZZ,129,YY+1,ZZ);
        setMapAt(130,YY,ZZ,130,YY+1,ZZ);
        setMapAt(131,YY,ZZ,131,YY+1,ZZ);
     }
   }
   for (YY = 205; YY < 216; YY ++){
     for (XX = 0; XX < dimX; XX ++){
       for (ZZ= 0; ZZ<10; ZZ ++){
        setMapAt(XX,YY,ZZ,XX+1,YY,ZZ);
       }
     }
   }

  for (ZZ= 0; ZZ<10; ZZ ++)
   for (YY = 224; YY < 238; YY ++)
   //for (YY = 182; YY < 195; YY ++)
     //for (XX = 30; XX < 45; XX ++)
     for (XX = 105; XX < 119; XX ++)
	setMapAt(XX,YY,ZZ, XX,YY,ZZ+1);

}

void VoxMap::DTIwarpExperiment1AppFiberBunch()
{
  int ROW, XX,YY, ZZ;
  float inix, iniy, alx, aly, xdir, ydir, px, py;
  double factor=0.01;

  cout <<" This creats a simulated curved fiber left to right";
  //cout <<" input a factor (for instance: "<<factor<<"):";
  //cin >>factor;
  //cout <<endl<<"You have inputted: "<<factor<<endl<<endl;

  inix = 1; iniy = 0.001 ; 


  for (ROW = 20; ROW < 80; ROW ++){
     xdir = inix; ydir = iniy;
     alx =0.0001; aly = 0.003;
     XX= px=0; YY=py = ROW;
     Vector v (xdir, ydir, 0.0f);	  
     v.norm();
     do {
	  //cout << "-- XX="<<XX<< ";v.getX()="<<v.getX()<<"; -- XX+v.getX()="<<XX+v.getX()<<endl;
	  //cout << "-- YY="<<YY<< ";v.getY()="<<v.getY()<<"; -- YY+v.getY()="<<YY+v.getY()<<endl;
	  for (ZZ = 0; ZZ < dimZ; ZZ++)
	     setMapAt(XX,YY,ZZ,XX+v.getX(),YY+v.getY(),ZZ);
          //adjust dir
	  if (XX > 80)
	     alx *= 0.995;
	  if (XX < 120)
	     aly *= (1.0 + factor);
	  else 
	     aly *= 1.001;
	  xdir += alx; ydir += aly;
	  // adjust position         
	   v.setXYZ (xdir, ydir, 0.0f);	  
	   //v.outputs();
	  v.norm();
	   //v.outputs();
	  do {
	    px += v.getX(); py += v.getY();
	  }while ( (int)px == XX  && (int)py == YY);
	  XX = (int)px; YY = (int)py;  
	  //cout << "  ## XX="<<XX<< ";v.getX()="<<v.getX()<< endl;
	  //cout << "  ## YY="<<YY<< ";v.getY()="<<v.getY()<< endl;
     }while ( XX < 256 && YY < 256 );
	
  }

}
