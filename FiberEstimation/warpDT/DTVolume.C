#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <iostream.h>
//#include <strstream.h>
#include <strstream>
#include <fstream.h>

#include "Global.h"
#include "VoxMap.h"
#include "WarpVolume.h"
#include "DTVolume.h"

double distanceGlobal(float a, float b, float c, struct Pt3d pp, float rx=1.0f, float ry=1.0f, float rz=1.0f)
{
   double tt = (a-pp.x)*(a-pp.x)*rx*rx+ry*ry*(b-pp.y)*(b-pp.y)+rz*rz*(c-pp.z)*(c-pp.z);
   return sqrt(tt);
}

//--------------------------------------------------------------
//
// class DTVolume: Volume of Diffusion Tensor
//
//--------------------------------------------------------------
DTVolume::DTVolume()
 :dimX(0),dimY(0),dimZ(0),dt(NULL),voxSize(sizeof(DTensor))
{  
  path[0] = '\0';
  planSize = dimX*dimY;
  tmpX = tmpY = tmpZ = 1.0;
  resX=resY=resZ=1.0;
}

DTVolume::DTVolume(char* name, int dimXX, int dimYY,
	float resXX, float resYY,float resZZ)
  :voxSize(sizeof(struct DTensor)),
   resX(resXX),resY(resYY),resZ(resZZ),
   dimX(dimXX),dimY(dimYY)//dimXX,dimYY default to 256
{
  planSize = dimX*dimY;
  strcpy(path, name);
  //tmpX = tmpY = tmpZ = 1.0;
  tmpX = resX*resX;
  tmpY = resY*resY;
  tmpZ = resZ*resZ;

  dt = NULL;
  loadData(path,dt,dimZ);
    cout <<" dt Address allocated : "<< dt<<endl;

  survey(dt, dimZ);
}

DTVolume::~DTVolume()
{
  if (dt != NULL){
     delete []dt;
     cout << "~DTVolume(): dt deleted" << endl;
  }
}
//		[a1, a4, a5]	   	[XX, XY, XZ]
//		[a4, a2, a6]   i.e.	[YX, YY, YZ]
//		[a5, a6, a3]		[ZX, ZY, ZZ]
void DTVolume::flipZComponent()
{
  int xxx,yyy,zzz,n;

  n = 0;
  for (xxx = 0; xxx < dimX; xxx++){
     for (yyy = 0; yyy < dimY; yyy ++){
	for (zzz = 0; zzz < dimZ; zzz ++, n ++){
	    dt[n].a[4] = -dt[n].a[4];
	    dt[n].a[5] = -dt[n].a[5];
	}
     }
  }
}

void DTVolume::newVolume(int dX,int dY, int dZ, float xr, float yr, float zr)
{
  int i, total, j;
  dimX = dX; dimY = dY; dimZ=dZ;
  planSize = dX*dY;
  total = planSize * dZ;

  if (dt != NULL){
      delete []dt;
      dt = NULL;
  }
  dt = new struct DTensor[total];
  if ( dt == NULL ) {
      cout << endl << " DTVolume::newVolume()";
      cout << endl << "  Memory allocation failed";
      cout << endl; 
      exit(0);
  }
  for (i=0;i<total;i++)
    for (j=0; j<6;j++)
	dt[i].a[j] = 0;

  setResolution(xr, yr, zr);
}

void DTVolume::setResolution(float xr, float yr, float zr)
{
  resX=xr; resY = yr; resZ = zr;
}

void DTVolume::saveTensorData(char *ff)
{
   saveData(ff,dt,dimZ*planSize);
}
void DTVolume::saveData(char* ff, struct DTensor* &object, int numofunit)
{
   FILE *fp;
   int n;

   cout << "in SaveData(): object = " <<object<<endl;
   fp = fopen(ff, "wb");
   if (fp == NULL ) {
       cout << endl << "file open ERROR: "<< ff << endl;
       exit (0);
   }else {       
	n = fwrite(object,voxSize, numofunit,fp);	
	fclose(fp);
	if (n!= numofunit){
	   cout << "file <"<<ff<<"> writing ERROR:" << endl;
	   cout << "\t expecting to write = "<<numofunit<< ", while written="<<n <<endl;
 	} else 
	   cout << "file <"<<ff<<"> saved OK" << endl;
  }
}
void DTVolume::survey(struct DTensor * &object, int sampZ)
{//check the percentage in the DT volume of DT without undefined value
 //extend this function to inspect more information of the DT volume
   int xxx,yyy,zzz;
   int i,k,total, sum=0, zv=0;
   int DEF;
   float res;

   total = dimX * dimY * sampZ;
   i = 0;
   for (xxx = 0; xxx < dimX; xxx++){
     for (yyy = 0; yyy < dimY; yyy ++){
	for (zzz = 0; zzz < sampZ; zzz ++){
	    DEF = 0;
	    if (object[i].a[4]>0) //check if ZZ is negative
		zv++;
	    for (k = 0; k < 6; k++){
		if (object[i].a[k] != 0){
		   DEF = 1;
		   break;
	        }
	    }
	    if (DEF)
		sum ++;
	    i++;
	}
     }
   }
   res = float(sum) / float(total);
   cout << " DT volume valued at " << res * 100.0 << " %" << endl;
   cout << " ZZ value negative = " << float(zv)/float(total)*100.0 <<"%"<<endl;
}

void DTVolume::loadData(char* ff, DTensor * &object, int& sampZ) //, BOOL UseFileSize)
{
   FILE *fp;
   int total,n;
   
   cout <<"----loadData()" << endl;
   fp = fopen(ff, "rb");
   if (fp == NULL ) {
       cout << endl << "file open ERROR: "<< ff << endl;
       exit (0);
   }else {       
     //fseek64(fp,0,SEEK_END);
     //total = ftell64(fp); 
     //sampZ= total/planSize/voxSize; //because we know slice is 256X256
     //fseek64(fp,0,SEEK_SET);

       fseek(fp,0,SEEK_END);
       total = ftell(fp); 
       sampZ= total/planSize/voxSize; //because we know slice is 256X256
       fseek(fp,0,SEEK_SET);

       //---Volume Data--------------
       if (object != NULL){
	 delete []object;
	 object = NULL;
       }
       object = new struct DTensor[total/voxSize];       
       //cout <<"total tensor number = "<< total/voxSize<<endl;
       if ( object == NULL ) {
           cout << endl << " DT Volume::loadData()";
	   cout << endl << "  Memory out when allocating for VOXEL";
	   cout << endl; 
           exit(0);
        }else {
	   cout << "\n Got memory for DT volume data\n" <<endl;
	   cout <<" Address allocated : "<< object<<endl;
	}
	n = fread(object,voxSize, total/voxSize,fp);	
	cout << endl << "\n Volume size = "<< total;
	cout << "; should read = " << total/voxSize << " (tensors), while read =" << n << endl;
	fclose(fp);
      	cout <<" --- Volume loaded, dimZ = "<< sampZ<<endl<<endl;
        /*-Begin--for debug purpose, a simulation DT data---
	 int XX, YY, ZZ,ind,k;
         for (n = 0; n < total/voxSize; n ++){
	    object[n].a[0] =  object[n].a[1] = object[n].a[2] =  0;
	    object[n].a[3] =  object[n].a[4] = object[n].a[5] =  0;
	 }
         for (ZZ = 20; ZZ < 40; ZZ ++){
            ind = ZZ * planSize;
            for (YY = 120; YY < 150; YY ++){
	       k = ind + YY * dimX;
               for (XX = 70; XX < 170; XX ++){
		       object[k+XX].a[0] =  object[k+XX].a[2] = 0.0001;
		       object[k+XX].a[1] =  0.0008;
	       }
            }
         }
        //-End--for debug purpose, a simulation DT data---
	*/
   }
}

void DTVolume::createSimulationDT(int XX1, int XX2, 
	int YY1, int YY2, int ZZ1, int ZZ2)
{
	 int XX, YY, ZZ,ind,k,n;
	 int total= planSize *dimZ;

         for (n = 0; n < total; n ++){
	    dt[n].a[0] =  dt[n].a[1] = dt[n].a[2] =  0;
	    dt[n].a[3] =  dt[n].a[4] = dt[n].a[5] =  0;
	 }
         for (ZZ = XX1; ZZ < XX2; ZZ ++){
            ind = ZZ * planSize;
            for (YY = YY1; YY < YY2; YY ++){
	       k = ind + YY * dimX;
               for (XX = ZZ1; XX < ZZ2; XX ++){
		       dt[k+XX].a[0] =  dt[k+XX].a[2] = 0.0001;
		       dt[k+XX].a[1] =  0.0008;
	       }
            }
         }
}
void DTVolume::printTensors()
{
   int XX, YY, ZZ, ind, k;

   cout <<" The tensors are:"<<endl;
   
   for (ZZ = 34; ZZ < 35; ZZ ++){
     ind = ZZ * planSize;
     for (YY = 120; YY < 122; YY ++){
	k = ind + YY * dimX;
        for (XX = 120; XX < 125; XX ++)	 
	    printTensorAt(k+XX);
     }
   }
}

void DTVolume::printTensorAt(int ind)
{
  int i;
  cout <<"Tensor["<<ind<<"]: ";
  for (i=0; i< 6; i++)
    cout <<"("<<dt[ind].a[i]<<")-";
  cout << endl;
}

void DTVolume::printTensorInMatrixAt(int ind)
{
  cout <<"Tensor["<<ind<<"]: "<<endl;
  cout <<"   "<< dt[ind].a[0] <<"  "  << dt[ind].a[3] << "  " << dt[ind].a[4]<<endl;
  cout <<"   "<< dt[ind].a[3] <<"  "  << dt[ind].a[1] << "  " << dt[ind].a[5]<<endl;
  cout <<"   "<< dt[ind].a[4] <<"  "  << dt[ind].a[5] << "  " << dt[ind].a[2]<<endl;
  cout << endl;
}

void DTVolume::printTensorAt(DTensor *DTS, int ind)
{
  int i;
  cout <<"Tensor["<<ind<<"]: ";
  for (i=0; i< 6; i++)
    cout <<"("<<DTS[ind].a[i]<<")-";
  cout << endl;
}

void DTVolume::getDimXYZ(int& xx, int& yy, int& zz)
{
  xx = dimX;
  yy = dimY;
  zz = dimZ;
}
void DTVolume::outPt(char* ss,struct Pt3d ppp)
{
   cout << ss << "=("<<ppp.x<<","<< ppp.y<<","<<ppp.z<<");"<<endl;
}
void DTVolume::forwardmapDT2atlas(VoxMap* & map1, WarpVolume* & map2,
	float rrx, float rry, float rrz)
{//map1: original space -> the space before warpped
 //map2: space before warpped -> atlas space fPoint pt;(57-slice)
   Matrix Ms1,M,M90,Mx,My,Mz,Ms2, Mt, Mj,Mjt;   
   struct Pt3d pt3, ptt;
   int k,XX,YY,ZZ, meaningful;
   struct Pt3d  *vMap = map1->getMap();

   cout << " ===forwardmapDT2atlas==OOOOOO==" << endl;
   //Ms1.newScaleXYZ(1.0,1.0,1.5/0.9766);  
   M90.loadIdentity(); M90.rotateX(-90/180.0*PI);//rotateX -90
   rrx = rrx/180*PI;  rry = rry/180*PI;  rrz = rrz/180*PI;
   Mx.loadIdentity();  My.loadIdentity();  Mz.loadIdentity();
   Mx.rotateX(rrx);  My.rotateY(rry);  Mz.rotateZ(rrz);
   //Ms2.newScaleXYZ(1.0,1.0, 0.9766/1.5);
   Mt = M = M90*Mx*My*Mz;
   Mt.transpose();
   k=0;
   for (ZZ = 0; ZZ < dimZ; ZZ ++){
     for (YY = 0; YY < dimY; YY ++){
       for (XX = 0; XX < dimX; XX ++){
         //if ( XX > 190) cout <<"XXYYZZ= ("<< XX << ";"<<YY<<";"<<ZZ <<")" << endl;
	 Matrix c(dt[k].a[0], dt[k].a[3], dt[k].a[4], 0.0f,
         	  dt[k].a[3], dt[k].a[1], dt[k].a[5], 0.0f,
         	  dt[k].a[4], dt[k].a[5], dt[k].a[2], 0.0f,
                  0.0f,        0.0f,        0.0f,        1.0f);
         c = Mt*c*M; //tensor rotated to the status before warpped
	 //c.outputs();
	 pt3 = map1->getMapAt(XX,YY,ZZ);//new position of the tensor
	 meaningful = map2->triInterpolationXYZ(pt3, ptt);
	 //if ( XX > 190) cout <<" 444444---" << endl;
 	 if (meaningful == -1){ //result undefined
	   ptt.x = ptt.y = ptt.z = UNDEFINED_PT;
	   dt[k].a[0] =  dt[k].a[3] = dt[k].a[2] = UNDEFINED_PT ;
   	   dt[k].a[4] =  dt[k].a[1] = dt[k].a[5] = UNDEFINED_PT ;  
	 }else{
		//if ( XX > 190) cout <<" 444555---" << endl;
           ptt.x += pt3.x; ptt.y += pt3.y; ptt.z += pt3.z;
	     //if ( XX > 192) outPt(" => ", ptt);
	   map2->JacobiRotationMatrix(ptt, Mj); //Jacobi-rotatation Matrix at 'Pt3'
		//if ( XX > 190) cout <<" 5555---" << endl;        
	   Mjt=Mj; Mjt.transpose();
	   c = Mjt*c*Mj; // c contains the new orianted tensor
   	   dt[k].a[0] = c.getValue(0,0); dt[k].a[3] = c.getValue(0,1);
   	   dt[k].a[4] = c.getValue(0,2); 
   	   dt[k].a[1] = c.getValue(1,1); dt[k].a[5] = c.getValue(1,2);
	   dt[k].a[2] = c.getValue(2,2);	 
	 }
    	 //update grid (XX,YY,ZZ) value with the destination in Atlas space
	 map1->setMapAt(XX,YY,ZZ, ptt);
 	 k++;
	}
     } 
   }
   cout << " ===forwardmapDT2atlas==PPPPP==" << endl;

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
float DTVolume::triInterpolation(float* w, float v0, 
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

struct DTensor * DTVolume::castDT2grid(VoxMap* & map1) 
{
   int n,k,i;
   int XX,YY,ZZ;
   struct DTensor *zDT, *yDT, *xDT;
   struct DTensor *nDT= NULL;
   struct Pt3d gridPts[20];
   float wei[20][8];
   int aPlanSize = dimX * dimX ; //rv
   int sum=0,sum1=0;

   nDT = new DTensor [dimX*dimY*dimZ]; //atlas DT volume; //rv 
   if (nDT == NULL){
     cout<<"ERROR:DTVolume::castDT2grid():atlas DT volume memory out" << endl;
     return NULL;
   }else
     cout << " nDT memory obtained" << endl;
   zDT = dt;
   for (ZZ = 0; ZZ < dimZ-1 ;ZZ++){
       cout << " --Casting Slice -" <<ZZ<<"-"<<endl;
       yDT = zDT;
       for (YY = 0 ; YY < dimY-1; YY ++){
          xDT = yDT;
	  for (XX = 0; XX < dimX-1; XX++){
	        //cout <<"XXYYZZ= ("<< XX << ";"<<YY<<";"<<ZZ <<")" << endl;
	      //calculate how many grid points in Atlas space (57-slice)
	      // is covered by the transoformed grid of Map1, the nubmer
	      // will be in 'n' and actually grid points at gridPts array.
	      // And wei[20] holds the relative distance weights for 
	      // interpolation
	      // n = -1 if there are too many than gridPts[20] can hold
              n = map1->getGridCoverd(XX,YY,ZZ, gridPts,wei,20); 
		//cout << "-- @@@ -- 222--n="<< n << endl;
	      sum += n;
	      for (k = 0; k < n ; k ++){ 
		/*
		 	cout << "(2) x1y1z1["<<k<<"]=("<<gridPts[k].x<<","<<gridPts[k].y
			<< ","<<gridPts[k].z<<")"<<endl;
			 for (int L=0;L<8;L++) {
				cout << "-QQQ-"<< endl;
				cout<< "~wei("<<L<<")"<<wei[k][L] << endl;
			 }*/
		i = gridPts[k].z * aPlanSize + dimX * gridPts[k].y  //rv
		    + gridPts[k].x;  
			//cout << "--MMMM--" << endl;
		sum1 += fillNewDT(xDT, gridPts[k], wei[k],nDT+i);
			//cout <<"---MMM-222-MMM--" <<endl;
	      }
		//cout <<" -- @@@ -- 333" << endl;
	      xDT ++;	
	  }

	  yDT += dimX;
       }
       zDT += planSize;
   }
   cout << "In the new DT volume: \n Coverd percentage = ";
   cout << float(sum)/(float(dimX*dimY*dimZ))*100.0 <<"% \n\n";
   cout << "Definited percentage = ";
   cout << float(sum1)/(float(dimY*dimX*dimZ))*100.0 <<"% \n";
   cout << "struct DTensor *nDT= " <<nDT<<endl;
   return nDT;
}
int DTVolume::fillNewDT(struct DTensor * &gridDT, struct Pt3d& grid, 
	float *www, struct DTensor *newDT) //float www[8]
{// 'gridDT' is the old grid Tensor at point 'grid'
 // www[8] is the weight to the cube at point 'grid'
 // resulting in 'newDT'
 // if undefined, return 0; defined, return 1;
   int cn, dltX, dltY, dltZ;
   if (grid.x > dimX-1)
	dltX = 0;
   else 
	dltX = 1;
   if (grid.y > dimY-1)
	dltY = 0;
   else 
	dltY = dimX;
   if (grid.z > dimZ-1)
	dltZ = 0;
   else 
	dltZ = planSize;

  if ((gridDT[0].a[0] == UNDEFINED_PT) 
	|| (gridDT[dltX].a[0] == UNDEFINED_PT)
	|| (gridDT[dltY+dltX].a[0] == UNDEFINED_PT) 
	|| (gridDT[dltY].a[0] == UNDEFINED_PT)
	|| (gridDT[dltZ].a[0] == UNDEFINED_PT)
	|| (gridDT[dltZ+dltX].a[0] == UNDEFINED_PT)
	|| (gridDT[dltZ+dltY+dltX].a[0] == UNDEFINED_PT)
	|| (gridDT[dltZ+dltY].a[0] == UNDEFINED_PT))
  { //newDT undefined
	  newDT->a[0] = newDT->a[1] = newDT->a[2]
	= newDT->a[3] = newDT->a[4] = newDT->a[5] 
	= 0; //UNDEFINED_PT;
	  return 0;
  } else{
    for (cn = 0; cn < 6; cn++){
    	newDT->a[cn] = weiInterpolation( www, 
		             gridDT[0].a[cn], 
		          gridDT[dltX].a[cn],
	             gridDT[dltY+dltX].a[cn],
		          gridDT[dltY].a[cn],
		          gridDT[dltZ].a[cn],
		     gridDT[dltZ+dltX].a[cn],
		gridDT[dltZ+dltY+dltX].a[cn],
		     gridDT[dltZ+dltY].a[cn] );
    }
    return 1;
  }
}

int DTVolume::isTensorWellDefined(struct DTensor * &dtGrid)
{// this function has been abondened, useless
  if ((dtGrid[0].a[0] == UNDEFINED_PT) || (dtGrid[1].a[0] == UNDEFINED_PT)
	|| (dtGrid[dimX+1].a[0] == UNDEFINED_PT) 
	|| (dtGrid[dimX].a[0] == UNDEFINED_PT)
	|| (dtGrid[planSize].a[0] == UNDEFINED_PT)
	|| (dtGrid[planSize+1].a[0] == UNDEFINED_PT)
	|| (dtGrid[planSize+dimX+1].a[0] == UNDEFINED_PT)
	|| (dtGrid[planSize+dimX].a[0] == UNDEFINED_PT))
    		return 0;
  else
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

float DTVolume::interTransGrid(DTensor* & dtGrid, 
             int cn,float ww[8])
{//trinterpolation within the 8 tensors starting at the
 //corner (XX,YY,ZZ)which result in the index of 'dtGrid' in the
 //tensor volume. the input will garuntee the grid will not
 //go beyond the volume, so there will have no inside/outside
 //detection
 //cn: DTensor's Component Number 
 float res;
 res = weiInterpolation(ww, dtGrid[0].a[cn], dtGrid[1].a[cn],
	dtGrid[dimX+1].a[cn],dtGrid[dimX].a[cn],
	dtGrid[planSize].a[cn],dtGrid[planSize+1].a[cn],
	dtGrid[planSize+dimX+1].a[cn],
	dtGrid[planSize+dimX].a[cn]);

 return res;
}

float DTVolume::weiInterpolation(float www[8], float A0,  float A1, float A2,
		float A3, float A4, float A5, float A6, float A7) 
{
  float t = A0*www[0]+A1*www[1]+A2*www[2]+A3*www[3]+A4*www[4]+A5*www[5]+
		A6*www[6]+A7*www[7];
  return t;
}

void DTVolume::reorientTensor(mVolume* &mv)
{//map: original space  -> atlas space fPoint pt;(57-slice)
 //mv:  volume of transformation matrix
   Matrix  U, tDT, Ut;   
   int k,XX,YY,ZZ;
  
   cout << " ===DTVolume::reorientTensor()==Begin==" << endl;
   k=0;
   for (ZZ = 0; ZZ < dimZ; ZZ ++){
     for (YY = 0; YY < dimY; YY ++){
       for (XX = 0; XX < dimX; XX ++, k++){
  	  /*if (XX!=126 || YY!=210 || ZZ != 0) ;
	  else {
	    cout <<" reorient() XYZ=("<<XX<<","<<YY<<","<<ZZ<<")"<<endl;
	  }*/
         //if ( XX > 190) cout <<"XXYYZZ= ("<< XX << ";"<<YY<<";"<<ZZ <<")" << endl;
	 Matrix c(dt[k].a[0], dt[k].a[3], dt[k].a[4], 0.0f,
         	  dt[k].a[3], dt[k].a[1], dt[k].a[5], 0.0f,
         	  dt[k].a[4], dt[k].a[5], dt[k].a[2], 0.0f,
                  0.0f,        0.0f,        0.0f,        1.0f);
         Ut = U = mv->getMatrix(k);
	 Ut.transpose();
	 tDT = U*c*Ut;
	   /********
	     if ((XX == 126) && (YY == 210) &&(ZZ ==0)){
		cout <<"k="<<k<<"---dt[k].a[4]="<<dt[k].a[4]<<"---"<<endl;
		cout <<"---U-------XX="<<XX<<"----------------"<<endl;
	        U.outputs();
		cout <<"---c---"<<endl;
	        c.outputs();
		cout <<"---Ut---"<<endl;
	        Ut.outputs();
		cout <<"---tDT---"<<endl;
	        tDT.outputs();
	     }
	  ********/
	 dt[k].a[0] = tDT.getValue(0,0); dt[k].a[3] = tDT.getValue(0,1);
   	 dt[k].a[4] = tDT.getValue(0,2); 
   	 dt[k].a[1] = tDT.getValue(1,1); dt[k].a[5] = tDT.getValue(1,2);
	 dt[k].a[2] = tDT.getValue(2,2);	 
	 //k++;
      }
    }
  }
}
void DTVolume::getSusumuMapLastStep(VoxMap* & map1, WarpVolume* & map2)
{//map1: original space -> the space before warpped
 //map2: space before warpped -> atlas space fPoint pt;(57-slice)
 // final map in Map1 (from original -> atlas 57-slice)
   struct Pt3d pt3, ptt;
   int k,XX,YY,ZZ, meaningful;
   struct Pt3d  *vMap = map1->getMap();

   cout << " ==DTVolume::getSusumuMapLastStep()==BEGIN==" << endl;
   k=0;
   for (ZZ = 0; ZZ < dimZ; ZZ ++){
     for (YY = 0; YY < dimY; YY ++){
       for (XX = 0; XX < dimX; XX ++){
	 pt3 = map1->getMapAt(XX,YY,ZZ);//new position of the tensor
	 meaningful = map2->triInterpolationXYZ(pt3, ptt);
 	 if (meaningful == -1){ //result undefined
	   ptt.x = UNDEFINED_PT; //XX;
	   ptt.y = UNDEFINED_PT; //YY;
	   ptt.z = UNDEFINED_PT; //ZZ;
	   //ptt.x = XX;
	   //ptt.y = YY;
	   //ptt.z = ZZ;
	 }
    	 //update grid (XX,YY,ZZ) value with the destination in Atlas space
	 map1->setMapAt(XX,YY,ZZ, ptt);
 	 k++;
	}
     } 
   }
   cout << " ===DTVolume::getSusumuMapLastStep()==END==" << endl;
   cout << " ---Susumu Map Step 2 completed ---" << endl << endl;
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
ww[8]
----------------------------------------------------*/
struct DTensor * DTVolume::warpTensorFieldRev(VoxMap* & map1,int warpedZ) 
{

  // for each destination location, look at 8 neighbors of source
  // compute tensor at destination by linearly interpolating eight sources

   int n,i, j;
   int XX,YY,ZZ, tx,ty,tz;
   float ww[8],tww;
   struct Pt3d  tDot;

   struct DTensor *zDT, *yDT, *xDT;
   struct DTensor *nDT= NULL;

   //rv   int aTotal, aPlanSize = atlasX * atlasY;
   int aTotal, aPlanSize = dimX * dimX;
   int sum=0;

   //int NNN ; //debug

   aTotal = aPlanSize*warpedZ; //total unit numer of atlas volume
   nDT = new DTensor [aTotal]; //[atlasX*atlasY*warpedZ]; //atlas DT volume;
   if (nDT == NULL){
     cout<<"ERROR:DTVolume::castDT2grid():atlas DT volume memory out" << endl;
     return NULL;
   }else
     cout << " nDT memory obtained" << endl;

   n = 0; //initialize nDT[] volume
   for (ZZ = 0; ZZ < warpedZ ;ZZ++) //for (ZZ = 0; ZZ < warpedZ-1 ;ZZ++)
      for (YY = 0; YY <dimX; YY++)  //for (YY = 0; YY < atlasY-1; YY++) //rv
   for (XX = 0; XX < dimX; XX++){ //for (XX = 0; XX < atlasX-1; XX++){ // rv
	nDT[n].a[0]=nDT[n].a[1]=nDT[n].a[2]=nDT[n].a[3]=nDT[n].a[4]=nDT[n].a[5]= 0;
	n++;
   }


   tmpX = resX*resX;//weight for calculate DTvolume::distance()
   tmpY = resY*resY;//weight for calculate DTvolume::distance()
   tmpZ = resZ*resZ;//weight for calculate DTvolume::distance()
	   cout <<"   tmpXYZ=("<<tmpX<<","<<tmpY<<","<<tmpZ<<")"<<endl;
   cout <<"\n\n"<<endl;
   int newn= 0;
   for (ZZ = 0; ZZ < dimZ ;ZZ++){  //for (ZZ = 0; ZZ < dimZ-1 ;ZZ++)
       cout << " --Warping DT Slice -" <<ZZ<<"    ----"<<char(13);flush(cout);

       for (YY = 0 ; YY < dimY; YY ++){ //for (YY = 0 ; YY < dimY-1; YY ++){

	  for (XX = 0; XX < dimX; XX++, newn++){ //  for (XX = 0; XX < dimX-1; XX++){
	   
	    tDot = map1->getVectorAt(XX,YY,ZZ);	    
	   
 	    if (  (tDot.x>dimX-1) || (tDot.x < 0)
	        ||(tDot.y>dimX-1) || (tDot.y < 0)
	        ||(tDot.z>warpedZ-1) || (tDot.z < 0) ){
	      sum ++;
	      continue;
	    }else {
	       
		tx = (int) tDot.x; ty = (int) tDot.y; tz = (int) tDot.z; 
		n = tz * aPlanSize + dimX * ty + tx; // n is the base of the cube
		// tri-linear interpolation
		
		if(aPlanSize +n+1+dimX >= aTotal)
		  {
		    sum ++;
		    continue;
		  }

		float tl_x = tDot.x- tx;
		float tl_y = tDot.y- ty;
		float tl_z = tDot.z- tz;
				
		ww[0] = (1-tl_x)*(1-tl_y)*(1-tl_z);
		ww[1] = tl_x*(1-tl_y)*(1-tl_z);
		ww[2] = tl_x*tl_y*(1-tl_z);
		ww[3] = (1-tl_x)*tl_y*(1-tl_z);
		ww[4] = (1-tl_x)*(1-tl_y)*tl_z;
		ww[5] = tl_x*(1-tl_y)*tl_z;
		ww[6] = tl_x*tl_y*tl_z;
		ww[7] = (1-tl_x)*tl_y*tl_z;
		
		for(i=0; i<6; i++)		  
		  nDT[newn].a[i] = ww[0]*dt[n].a[i] + ww[1]*dt[n+1].a[i] + ww[2]*dt[n+1+dimX].a[i] + ww[3]*dt[n+dimX].a[i] + ww[4]*dt[aPlanSize + n].a[i] + ww[5]*dt[aPlanSize +n+1].a[i] + ww[6]*dt[aPlanSize +n+1+dimX].a[i] + ww[7]*dt[aPlanSize +n+dimX].a[i] ;

	    }

	  }

       }

   }	

   cout << "\n\nIn the old DT volume: \n Coverd percentage = ";
   cout << 100.0 - float(sum)/(float(dimX*dimX*warpedZ))*100.0 <<"% \n\n";
   cout << "struct DTensor *nDT= " <<nDT<<endl;
   cout << "NOTE: The warped NEW DT volume dimensions\n        =(";
   cout << dimX<<"x"<<dimY<<"x"<<warpedZ<<")"<<endl<<endl;
   return nDT;
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
ww[8]
----------------------------------------------------*/
struct DTensor * DTVolume::warpTensorField(VoxMap* & map1,int warpedZ) 
{// 1. take a tensor D and its destination P
 // 2. get the grid that contains P
 // 3. distribute D to the 8 grid points
   int n,i;
   int XX,YY,ZZ, tx,ty,tz;
   float ww[8],tww;
   struct Pt3d  tDot;

   struct DTensor *zDT, *yDT, *xDT;
   struct DTensor *nDT= NULL;

   //rv   int aTotal, aPlanSize = atlasX * atlasY;
   int aTotal, aPlanSize = dimX * dimX;
   int sum=0;

   //int NNN ; //debug

   aTotal = aPlanSize*warpedZ; //total unit numer of atlas volume
   nDT = new DTensor [aTotal]; //[atlasX*atlasY*warpedZ]; //atlas DT volume;
   if (nDT == NULL){
     cout<<"ERROR:DTVolume::castDT2grid():atlas DT volume memory out" << endl;
     return NULL;
   }else
     cout << " nDT memory obtained" << endl;

   n = 0; //initialize nDT[] volume
   for (ZZ = 0; ZZ < warpedZ ;ZZ++) //for (ZZ = 0; ZZ < warpedZ-1 ;ZZ++)
      for (YY = 0; YY <dimX; YY++)  //for (YY = 0; YY < atlasY-1; YY++) //rv
   for (XX = 0; XX < dimX; XX++){ //for (XX = 0; XX < atlasX-1; XX++){ // rv
	nDT[n].a[0]=nDT[n].a[1]=nDT[n].a[2]=nDT[n].a[3]=nDT[n].a[4]=nDT[n].a[5]= 0;
	n++;
   }

   zDT = dt;
   tmpX = resX*resX;//weight for calculate DTvolume::distance()
   tmpY = resY*resY;//weight for calculate DTvolume::distance()
   tmpZ = resZ*resZ;//weight for calculate DTvolume::distance()
	   cout <<"   tmpXYZ=("<<tmpX<<","<<tmpY<<","<<tmpZ<<")"<<endl;
   cout <<"\n\n"<<endl;
   for (ZZ = 0; ZZ < dimZ ;ZZ++){  //for (ZZ = 0; ZZ < dimZ-1 ;ZZ++)
       cout << " --Warping DT Slice -" <<ZZ<<"    ----"<<char(13);flush(cout);
       yDT = zDT;
       for (YY = 0 ; YY < dimY; YY ++){ //for (YY = 0 ; YY < dimY-1; YY ++){
          xDT = yDT;
	  for (XX = 0; XX < dimX; XX++,xDT++){ //  for (XX = 0; XX < dimX-1; XX++){
	    //if (ZZ==0){
		//cout <<"   XYZ=("<<XX<<","<<YY<<","<<ZZ<<")"<<endl;
	    //}
	    tDot = map1->getVectorAt(XX,YY,ZZ);	    
	    /*  if ((XX == 2) && (YY == 0) &&(ZZ ==0)){
			cout <<" #### DTVolume::warpTensorField():"<<endl;
			cout << "  AtlasXY= ("<<atlasX<<","<<atlasY<<endl;
	    		cout <<"   XYZ=("<<XX<<","<<YY<<","<<ZZ<<")"<<endl;
	    		cout <<"     tDot vec=("<<tDot.x<<","<<tDot.y<<","<<tDot.z<<")"<<endl;

			}*/
	    //if ((ZZ==0) &&(YY==00)) cout <<"  $$ 11 $$"<<endl;
 	    if (  (tDot.x>dimX-1) || (tDot.x < 0)
	        ||(tDot.y>dimX-1) || (tDot.y < 0)
	        ||(tDot.z>warpedZ-1) || (tDot.z < 0) ){
	      sum ++;
	      continue;
	    }else {
	        //if (ZZ==0) cout <<"  $$ 22 $$"<<endl;
	        /*if ((XX == 126) && (YY == 210) &&(ZZ ==0)){ //** fordebug
			cout <<" $$$$$ ---111---$$$"<<endl;
		}*/  
		tx = (int) tDot.x; ty = (int) tDot.y; tz = (int) tDot.z; 
		n = tz * aPlanSize + dimX * ty + tx; // n is the base of the cube
		/* for debug
	        if ((XX == 126) && (YY == 210) &&(ZZ ==0)){ //** fordebug
			cout <<" $$$$$ ---222---$$$"<<endl;
	    		cout <<"   txyz=("<<tx<<","<<ty<<","<<tz<<")"<<endl;
			cout <<"   n = tz * aPlanSize + atlasX * ty + tx="<<n<<endl;
		}  for debug */
		ww[0] = distance(tx  , ty  , tz  , tDot);
		ww[1] = distance(tx+1, ty  , tz  , tDot);
		ww[2] = distance(tx+1, ty+1, tz  , tDot);
		ww[3] = distance(tx  , ty+1, tz  , tDot);
		ww[4] = distance(tx  , ty  , tz+1, tDot);
		ww[5] = distance(tx+1, ty  , tz+1, tDot);
		ww[6] = distance(tx+1, ty+1, tz+1, tDot);
		ww[7] = distance(tx  , ty+1, tz+1, tDot);
		tww = 0.0;
		for (i = 0; i<8;i++)
			tww += ww[i];
	        if ((XX == 2) && (YY == 0) &&(ZZ==0)){
		    cout <<"  $$ 33 $$"<<endl;
			for (int h=0;h<8;h++)
	    			cout <<"   ww["<<h<<"]= "<<ww[h]<<endl;
			cout <<"   tww="<<tww<<";  aTotal="<<aTotal<<endl;
		}
	        /* 
		if ((XX == 126) && (YY == 210) &&(ZZ ==0)){ // fordebug
			cout <<" $$$$$ ---333---$$$"<<endl;
			for (int h=0;h<8;h++)
	    			cout <<"   ww["<<h<<"]= "<<ww[h]<<endl;
			cout <<"   tww="<<tww<<";  aTotal="<<aTotal<<endl;
		
		}// for debug 
		*/
		if (n < aTotal)
		    distributeDT(xDT, ww[0]/tww,nDT+n);
		 //if (ZZ==0)   cout <<"  $$ 44 $$"<<endl;
		if (n+1 < aTotal)
		    distributeDT(xDT, ww[1]/tww,nDT+n+1);
		 //if (ZZ==0)   cout <<"  $$ 55 $$"<<endl;
		if (n+1+dimX < aTotal)
		    distributeDT(xDT, ww[2]/tww,nDT+n+1+dimX);
		 //if (ZZ==0)   cout <<"  $$ 66 $$"<<endl;
		if (n+dimX< aTotal)
		    distributeDT(xDT, ww[3]/tww,nDT+n+dimX);
		/*
		if ((XX == 126) && (YY == 210) &&(ZZ ==0)){ // fordebug
		    cout << " ^^^ Old DT^^^^"<<endl; NNN=n;
		    printTensorAt(xDT,0);
		    cout << " ^^^ New DT^^^^n = "<<n<<endl;		    
		    printTensorAt(nDT,n);
		    printTensorAt(nDT,n+1);
		    printTensorAt(nDT,n+1+atlasX);
		    printTensorAt(nDT,n+atlasX);

		}*/
	        //if (ZZ==0) cout <<"  $$ 77 $$"<<endl;
		n += aPlanSize;
		if (n < aTotal)
		    distributeDT(xDT, ww[4]/tww,nDT+n);
		if (n+1 < aTotal)
		    distributeDT(xDT, ww[5]/tww,nDT+n+1);
		if (n+1+dimX < aTotal)
		    distributeDT(xDT, ww[6]/tww,nDT+n+1+dimX);
		if (n+dimX< aTotal)
		    distributeDT(xDT, ww[7]/tww,nDT+n+dimX);
		/*
		if ((XX == 126) && (YY == 210) &&(ZZ ==0)){ // for debug
		    cout << " ^^^ Old DT^^^^"<<endl;
		    printTensorAt(xDT,0);
		    cout << " ^^^ New DT^^^^n = "<<n<<endl;		    
		    printTensorAt(nDT,n);
		    printTensorAt(nDT,n+1);
		    printTensorAt(nDT,n+1+atlasX);
		    printTensorAt(nDT,n+atlasX);

		}*/


	    }
            //xDT ++;	
	  }
	  yDT += dimX;
       }
       zDT += planSize;
   }
		/** for debug
		    cout << " ^^^ ^^^finaly^"<<endl; NNN=n;
		    cout << " ^^^ New DT^^^^NNN = "<<NNN<<endl;		    
		    printTensorAt(nDT,NNN);
		    printTensorAt(nDT,NNN+1);
		    printTensorAt(nDT,NNN+1+atlasX);
		    printTensorAt(nDT,NNN+atlasX);
		*/
		


   cout << "\n\nIn the old DT volume: \n Coverd percentage = ";
   cout << 100.0 - float(sum)/(float(dimX*dimX*warpedZ))*100.0 <<"% \n\n";
   cout << "struct DTensor *nDT= " <<nDT<<endl;
   cout << "NOTE: The warped NEW DT volume dimensions\n        =(";
   cout << dimX<<"x"<<dimY<<"x"<<warpedZ<<")"<<endl<<endl;
   return nDT;
}

void DTVolume::distributeDT(struct DTensor * &oldDT,
		float www, struct DTensor *newDT) 
{// distribute old tensor to new according to the weight www
   int i;

   for (i=0; i<6; i++){
     newDT->a[i] = 0.5*(newDT->a[i] + oldDT->a[i]);
	// newDT->a[i] += oldDT->a[i] * www;
   }
}
double DTVolume::distance(float a, float b, float c, struct Pt3d pp)
{ // preset: tmpX = resX*resX; tmpY = resY*resY; tmpZ=resZ*resZ
   double tt = tmpX*(a-pp.x)*(a-pp.x)+tmpY*(b-pp.y)*(b-pp.y)+tmpZ*(c-pp.z)*(c-pp.z);
   return sqrt(tt);
}
double DTVolume::distanceDebug(float a, float b, float c, struct Pt3d pp)
{ // preset: tmpX = resX*resX; tmpY = resY*resY; tmpZ=resZ*resZ
	cout <<"  DTVolume::distanceDebug(): tmpXYZ=("<<tmpX<<","<<tmpY<<","<<tmpZ<<")"<<endl;
   double tt = tmpX*(a-pp.x)*(a-pp.x)+tmpY*(b-pp.y)*(b-pp.y)+tmpZ*(c-pp.z)*(c-pp.z);
	cout <<"  DTVolume::dist=("<<a-pp.x<<","<<b-pp.y<<","<<c-pp.z<<")"<<endl;

   return sqrt(tt);
}

void DTVolume::averageAndSave(int nSlc, char * toStore)
{// this function  takes dt[] as a stack of a number of volumes
 // nSlc is the number of slice of each volume, therefore dimZ = nSlc * n
 // the task of this function is to average the volume and get a resulting average volume
 // the result is of (dimX x dimY x nSlc) and will be store to  file "toStore"

  int dti,i, XX, YY, ZZ, k, pk, volCt;
  int interval = nSlc *planSize; // size of the resulting volume
  struct DTensor *resVol= NULL;
  float  tdti[6];

  volCt = (float)dimZ / (float)nSlc; //get the number of volume this stack has
  if ( volCt != (float)dimZ / (float)nSlc){
    cout <<"ERROR: the volume stack is not completed: Stack slice number= ";
    cout <<dimZ<<", average Slc = "<<nSlc<<endl;    
    return;
  } else{
    cout << " The stack contains " << volCt << " volumes" << endl;
  }
  if (volCt < 2){
    cout << " No need to do average, program quits."<<endl;
    return;
  }
  resVol = new struct DTensor[interval];
  if (resVol == NULL){
     cout <<" VoxMap::averageEveryNslice():: ERROR : memory out"<<endl;
     return;
  } else{
     cout << " VoxMap::averageEveryNslice(): got memory for result"<<endl;
  }
  //saveVecVolume(resVol,nSlc,toStore);

  k =0;
  for (ZZ = 0; ZZ < nSlc; ZZ ++){
     for (YY = 0; YY < dimY; YY ++){
       for (XX = 0; XX < dimX; XX ++){
	   pk = k;
	   tdti[0] = tdti[1] = tdti[2] = tdti[3] = tdti[4] = tdti[5] = 0.0f;
           for (i = 0; i < volCt ; i++, pk += interval ){//loop for every sub-volume
	      if (dt[pk].a[0] >=0){
	        for (dti=0; dti<6; dti++){ //DTI has 6 components
	            tdti[dti] += dt[pk].a[dti];
	        }
	      }else{
	        for (dti=0; dti<6; dti++){ //DTI has 6 components
	             tdti[dti] += -dt[pk].a[dti];
		}
	      }
	   }
           for (dti=0; dti<6; dti++){ //DTI has 6 components
	        resVol[k].a[dti] = tdti[dti] /(float)volCt;
	   }
	 k++;
       }
     }
  }
  saveData(toStore, resVol, interval);
  delete []resVol;
}

int DTVolume::averageFA_Deviation(int nSlc, float* &averFA, float* &deviation)
{// this function  takes dt[] as a stack of a number of volumes
 // nSlc is the number of slice of each volume, therefore dimZ = nSlc * n
 // the task of this function is to average the volume's FA and get a resulting average volume
 // of FA and its standard deviation. The result averFA and deviation is of (dimX x dimY x nSlc) 
 // and will be returned in the reference varibles to caller
 // The size of the volume of FA and deviation will be returned (interval).


  int i, XX, YY, ZZ, k, pk, volCt;
  int interval = nSlc *planSize; // size of the resulting volume
  double tFA, tDvi, *dtFA=NULL;

  volCt = (float)dimZ / (float)nSlc; //get the number of volume this stack has
  if ( volCt != (float)dimZ / (float)nSlc){
    cout <<"ERROR: the volume stack is not completed: Stack slice number= ";
    cout <<dimZ<<", average Slc = "<<nSlc<<endl;    
    return -1;
  } else{
    cout << " The stack contains " << volCt << " volumes" << endl;
  }
  if (volCt < 2){
    cout << " No need to do average, program quits."<<endl;
    return -1;
  }
  if (averFA!= NULL)
	delete []averFA;
  if (deviation!=NULL)
	delete []deviation;
  averFA=deviation=NULL;
  
  averFA = new float[interval];
  deviation = new float[interval];
  dtFA=new double [volCt];

  if (averFA == NULL || deviation == NULL || dtFA==NULL){
     cout <<" VoxMap::averageEveryNslice():: ERROR : memory out"<<endl;
     return -1;
  } else{
     cout << " VoxMap::averageFA_Deviation(): got memory for result"<<endl;
  }

  float **A, *eVal, **eVec;

  A = matrix (1,3, 1, 3);	// index from 1 ~3
  eVal = vector(1,3);		// index from 1 ~3
  eVec = matrix(1,3, 1, 3);	// index from 1 ~3

  k =0;
  for (ZZ = 0; ZZ < nSlc; ZZ ++){
     for (YY = 0; YY < dimY; YY ++){
       for (XX = 0; XX < dimX; XX ++){
	   pk = k;
	   tFA=0.0f;
           for (i = 0; i < volCt ; i++, pk += interval ){//loop for every sub-volume
	      dtFA[i] = getDT_FA(dt[pk].a, A, eVal, eVec);
 	      tFA += dtFA[i]; 	      
	   }
	   averFA[k]= tFA/(float)volCt;
	   tFA = averFA[k];

	   tDvi=0;
           for (i = 0; i < volCt ; i++ ){//loop for every sub-volume
	      tDvi += (tFA-dtFA[i])*(tFA-dtFA[i]); 	      
	   }
	   deviation[k] = sqrt(tDvi/((float)(volCt-1)));
	   k++;
       }
     }
  }
  //saveData(toStore, resVol, interval);
  cout << " DTVolume::averageFA_Deviation():: averFA="<<averFA<<";  deviation="<<deviation<<endl;
  delete []dtFA;

  free_matrix(A,1,3, 1, 3);
  free_vector(eVal,1,3);
  free_matrix(eVec,1,3, 1, 3);

  return interval;
}

float * DTVolume::createFAvolume()
{// this function return a volume of FA valume corresponding to 
 // the DT tensors at each voxel, 
 // return the pointer to the FA volume

  int XX, YY, ZZ, k;
  float *volFA=NULL;

  volFA = new float[dimX*dimY*dimZ];

  if (volFA == NULL){
     cout <<" VoxMap::createFAvolume():: ERROR : memory out"<<endl;
     return NULL;
  } else{
     cout << " VoxMap::createFAvolume(): got memory for result"<<endl;
  }

  float **A, *eVal, **eVec;

  A = matrix (1,3, 1, 3);	// index from 1 ~3
  eVal = vector(1,3);		// index from 1 ~3
  eVec = matrix(1,3, 1, 3);	// index from 1 ~3

  k =0;
  cout <<endl;
  for (ZZ = 0; ZZ < dimZ; ZZ ++){
     //cout << "   ...Processing Z = "<<ZZ<<endl;
     cout << "   ...Processing Z = "<<ZZ<<char(13);
     flush(cout);
     for (YY = 0; YY < dimY; YY ++){
       for (XX = 0; XX < dimX; XX ++, k++){
	   volFA[k] = getDT_FA(dt[k].a, A, eVal, eVec);
       }
     }
  }  
  //saveData(toStore, resVol, interval);
  cout << "\n\n DTVolume::createFAvolume():: volFA="<<volFA<<endl;

  free_matrix(A,1,3, 1, 3);
  free_vector(eVal,1,3);
  free_matrix(eVec,1,3, 1, 3);

  return volFA;
}
double DTVolume::getDT_FA(float tens[6],float **A, float *eVal, float **eVec)
{
  double  t = 0;

  for (int m=0;m<6;m++)
	t += tens[m];
  if (t == 0){
         A[1][1] =  A[1][2] =  A[1][3] = 0;
         A[2][1] =  A[2][2] =  A[2][3] = 0;
         A[3][1] =  A[3][2] =  A[3][3] = 0;
  }else {
         A[1][1] = tens[0]/t; A[1][2] = tens[3]/t; A[1][3] = tens[4]/t;
         A[2][1] = A[1][2]  ; A[2][2] = tens[1]/t; A[2][3] = tens[5]/t;
         A[3][1] = A[1][3]  ; A[3][2] = A[2][3]  ; A[3][3] = tens[2]/t;
  }
  int rr;
  rr = eigOfTensor(A,eVal,eVec,3);
  if (rr==0){
    cout <<"   bad Tensor"<<endl;
    //printTensorAt(n);
    //seeMatrix(A,3,3,"A matrix:");
    //bad ++;
    //if (bad > 10)
    //	exit(0);
  }

  double W123 = eVal[1]*eVal[1] + eVal[2]*eVal[2] + eVal[3]*eVal[3];
  double FA,ww;

  if (W123 < ZERO){
	FA = 0.0;
  }else {
	ww = fabs(eVal[1]*eVal[2])+fabs(eVal[2]*eVal[3])+fabs(eVal[3]*eVal[1]);
	FA =  1.0 - ww/W123;  // could be negative, due to the computer accuracy
	if (FA < 0) 
	   FA = 0;
	else 
	   FA = sqrt (FA);
  }
  return FA;
}

int DTVolume::getTensorAt(int col, int row, int slice, float **M)
{ //M is 3 x 3 matrix
  //return 1 when M is defined; otherwise return 0
  // NOTE: M subscript uses from 1 ~ N
  //
   if (( col <0 ) || (col > dimX) 
	||(row<0) || (row > dimY)
	||(slice <0) ||(slice >dimZ)) {
      M[1][1] = M[1][2] = M[1][3] = 0; 
      M[2][1] = M[2][2] = M[2][3] = 0; 
      M[3][1] = M[3][2] = M[3][3] = 0; 
      return 0;  // returned value undefined
   } else{   
      int k = planSize * slice + dimX *row + col;
      M[1][1]= dt[k].a[0]; M[1][2] = dt[k].a[3]; M[1][3] = dt[k].a[4];
      M[2][1]= dt[k].a[3]; M[2][2] = dt[k].a[1]; M[2][3] = dt[k].a[5];
      M[3][1]= dt[k].a[4]; M[3][2] = dt[k].a[5]; M[3][3] = dt[k].a[2];
      return 1; // returned value defined in DT1
   }
}
int DTVolume::setTensorAt(int col, int row, int slice, DTensor DDD)
{ //M is 3 x 3 matrix
  //return 1 when M is defined; otherwise return 0
  // NOTE: M subscript uses from 1 ~ N
  //
   if (( col <0 ) || (col > dimX) 
	||(row<0) || (row > dimY)
	||(slice <0) ||(slice >dimZ)) {
      return 0;  // returned value undefined
   } else{   
      int k = planSize * slice + dimX *row + col;
      for (int i=0; i<6; i++)
        dt[k].a[i]=DDD.a[i];
      return 1; //successful
   }
}
int DTVolume::setTensorAt(int col, int row, int slice, float vvv[6])
{ //M is 3 x 3 matrix
  //return 1 when M is defined; otherwise return 0
  // NOTE: M subscript uses from 1 ~ N
  //
   if (( col <0 ) || (col > dimX) 
	||(row<0) || (row > dimY)
	||(slice <0) ||(slice >dimZ)) {
      return 0;  // returned value undefined
   } else{   
      int k = planSize * slice + dimX *row + col;
      for (int i=0; i<6; i++)
        dt[k].a[i]=vvv[i];
      return 1; //successful
   }
}

int DTVolume::getTensorAt(int col, int row, int slice, DTensor &DT1)
{  //return 1 when DT1 is defined; otherwise return 0

   if (( col <0 ) || (col > dimX) 
	||(row<0) || (row > dimY)
	||(slice <0) ||(slice >dimZ)) {
      DT1.a[0]=0;
      DT1.a[1]=0;
      DT1.a[2]=0;
      DT1.a[3]=0;
      DT1.a[4]=0;
      DT1.a[5]=0; 
      return 0;  // returned value undefined
   } else{   
      int k = planSize * slice + dimX *row + col;
      DT1= dt[k];
      return 1; // returned value defined in DT1
   }
}
void DTVolume::setTensorAt(int index,
	float v1, float v2,	float v3, 
	float v4, float v5, 	float v6)

{
    dt[index].a[0]=v1;
    dt[index].a[1]=v2;
    dt[index].a[2]=v3;
    dt[index].a[3]=v4;
    dt[index].a[4]=v5;
    dt[index].a[5]=v6; 
}
int DTVolume::isZEROTensor(int index)
{
   if (dt[index].a[0]==0
	&& dt[index].a[1]==0
	&& dt[index].a[2]==0
	&& dt[index].a[3]==0
	&& dt[index].a[4]==0
	&& dt[index].a[5]==0)
      return 1;
   else
      return 0;
}
void DTVolume::setAverageTensorAt(int index,
	float v1, float v2,	float v3, 
	float v4, float v5, 	float v6)
{//need to average the original values
    dt[index].a[0]=0.5*(dt[index].a[0]+v1);
    dt[index].a[1]=0.5*(dt[index].a[1]+v2);
    dt[index].a[2]=0.5*(dt[index].a[2]+v3);
    dt[index].a[3]=0.5*(dt[index].a[3]+v4);
    dt[index].a[4]=0.5*(dt[index].a[4]+v5);
    dt[index].a[5]=0.5*(dt[index].a[5]+v6); 
}
void DTVolume::saveVectorFile(struct Pt3d * vf, int XX, int YY, int ZZ, char *ff)
{
   char buf[200];
   FILE *fp = fopen(ff,"wb");
   int n,len = XX*YY*ZZ;
   if (fp == NULL){
	cout <<"DTVolume::SaveVectorFile(): File <"<<ff<<"> open ERROR"<<endl;
	exit(0);
   }else
	cout <<" Saving to file <"<<ff<<">....."<<endl;

   n=fwrite(vf, sizeof(struct Pt3d), len, fp);
   cout << " Shoud write "<< len<<" unit(s), while written "<<n<<" unit(s)"<<endl;
   if (n != len)
	cout << "File saved INCORRECT!!"<<endl;
   else
	cout << "File saved !"<<endl;
   fclose (fp);	

   strcpy(buf,ff);
   strcat(buf,".info");
   fp = fopen(buf,"w");
   if (fp == NULL){
	cout <<"DTVolume::SaveVectorFile(): File <"<<buf<<"> open ERROR"<<endl;
	exit(0);
   }else
	cout <<" Writing info file <"<<buf<<">....."<<endl;
   fprintf(fp, "This is an info file for :%s (vector file of struct Pt3d)\n", ff);
   fprintf(fp, "Dimension = (%d,%d,%d)\n",XX,YY,ZZ);
   fclose(fp);
   cout <<" Writing completed"<<endl<<endl;
}
//*********************************************************
// the Primary Direction(PD)is using the following:
//
//      D123 = D(1,1)^2 +  D(2,2)^2 + D(3,3)^2;
//      t = (abs(D(1,1))-abs(D(2,2)))^2
//	  + (abs(D(3,3))-abs(D(2,2)))^2 
//	  + (abs(D(1,1))-abs(D(3,3)))^2
//      weiD = sqrt(t/2/D123)
// D1 D2 D3 are the eigen values
// A simplified version is:
// 
//    weiD = sqrt(1- x);
//    x    = y / D123
//    y    = fabs(D1*D2)+fabs(D2*D3)+fabs(D3*D1)
//
//*********************************************************
struct Pt3d * DTVolume::PrimaryDirection(char *PDfilename, int whichPD)
{//  if PDfilename == NULL : return the pointer to the PD vectors
 //  if PDfilename != NULL : save the extracted PD vector to file
 //			     and return NULL (allocated PD space freed)
 //  which PD: defines which PD to be extracted; (default == 1)
 //          == 1: the PD of the tensor  
 //	     == 2: the 2nd PD of the tensor
 //	     == 3: the 3rd (least) PD of the tensor

   int XX,YY,ZZ,n, i,iMax,bad,iMin;
   float **A, *eVal, **eVec, Max,Mini;
   //DTensor *zDT, *yDT, *xDT;
   //struct Pt3d  *xPD,*yPD,*zPD;
   struct Pt3d * PD = NULL;
   double weight,W123,ww,t;

   PD = new struct Pt3d [dimX*dimY*dimZ];
   if (PD == NULL){
	cout <<" DTVolume::PrimaryDirection(): Memory allocation failed"<<endl;
        exit(0);
   }else{
	cout <<" Allocated "<< dimX*dimY*dimZ<<" units"<<endl;
   }

   A = matrix (1,3, 1, 3);	// index from 1 ~3
   eVal = vector(1,3);		// index from 1 ~3
   eVec = matrix(1,3, 1, 3);	// index from 1 ~3

   //zDT = dt;
   //zPD = PD;

   cout <<endl;
   n=0; bad=0;
   for (ZZ = 0; ZZ < dimZ ;ZZ++){  //for (ZZ = 0; ZZ < dimZ-1 ;ZZ++)
       cout << "   Processing slice Z = " <<ZZ<<char(13);
       flush(cout);
       //yDT = zDT;
       //yPD = zPD;
       for (YY = 0 ; YY < dimY; YY ++){ //for (YY = 0 ; YY < dimY-1; YY ++){
          //xDT = yDT;
	  //xPD = yPD
	  for (XX = 0; XX < dimX; XX++,n++){ //  for (XX = 0; XX < dimX-1; XX++){
	    //if (YY>55 && YY < 70 && XX == 123){
		//cout <<"(XYZ)=("<<XX<<","<<YY<<","<<ZZ<<"):=>\n  ";
		//printTensorInMatrixAt(n);
	    //}
	     
 	     t = 0;
	     for (int m=0;m<6;m++)
		t += dt[n].a[m];
	     if (t == 0){
               A[1][1] =  A[1][2] =  A[1][3] = 0;
    	       A[2][1] =  A[2][2] =  A[2][3] = 0;
    	       A[3][1] =  A[3][2] =  A[3][3] = 0;
	     }else {
               A[1][1] = dt[n].a[0]/t; A[1][2] = dt[n].a[3]/t; A[1][3] = dt[n].a[4]/t;
    	       A[2][1] = A[1][2]     ; A[2][2] = dt[n].a[1]/t; A[2][3] = dt[n].a[5]/t;
    	       A[3][1] = A[1][3]     ; A[3][2] = A[2][3]     ; A[3][3] = dt[n].a[2]/t;
	     }
             //A[1][1] = dt[n].a[0]; A[1][2] = dt[n].a[3]; A[1][3] = dt[n].a[4];
    	     //A[2][1] = dt[n].a[3]; A[2][2] = dt[n].a[1]; A[2][3] = dt[n].a[5];
    	     //A[3][1] = dt[n].a[4]; A[3][2] = dt[n].a[5]; A[3][3] = dt[n].a[2];
	     if ((0==eigOfTensor(A,eVal,eVec,3))&&(bad<10)){
		 //cout <<"   bad Tensor["<<n<<"]"<<endl;
		 //printTensorAt(n);
		 //seeMatrix(A,3,3,"A matrix:");
		 bad ++;
		 //if (bad > 10)
		//	exit(0);
	     }
	     W123 = eVal[1]*eVal[1] + eVal[2]*eVal[2] + eVal[3]*eVal[3];
	     if (W123 < ZERO){
		weight = 0.0;
	     }else {
		ww = fabs(eVal[1]*eVal[2])+fabs(eVal[2]*eVal[3])+fabs(eVal[3]*eVal[1]);
		weight =  1.0 - ww/W123;  // could be negative, due to the computer accuracy
		if (weight < 0) 
			weight = 0;
		else 
		    weight = sqrt (weight);
	     }
	     Max = Mini = fabs(eVal[1]);
             iMax = iMin = 1;
	     for (i = 1; i<= 3; i++){  // find Primary Direction
	       if (fabs(eVal[i])>Max){
		   Max=fabs(eVal[i]);
		   iMax = i;
	       }
	       if (fabs(eVal[i]) < Mini){
		   Mini = fabs(eVal[i]);
		   iMin = i;
	       }
	     }	
	     if (iMax != iMin)		
	     switch (whichPD){
	        case 1:  // Primary Direction
		  break;
	        case 2:  // need second PD component
		  //cout <<"-----------\n";
		  //cout <<" ("<<eVal[1]<<","<<eVal[2]<<","<<eVal[3]<<")\n";
		  //cout << "iMin = "<< iMin<<"; iMax="<<iMax<<endl;
		  iMax = 6-iMin-iMax; //cout << "iMax = "<< iMax<<endl;
		  break;
	        case 3: // need the least PD component
		  iMax = iMin;
		  break;
	     }	
	     //get PD vector (weighted)
	     PD[n].x = weight*eVec[1][iMax] + XX;
	     PD[n].y = weight*eVec[2][iMax] + YY;
	     PD[n].z = weight*eVec[3][iMax] + ZZ;

	     //xPD ++;
	     //yDT ++;
	     //n++;
	  }
	  //yDT += dimX;
	  //yPD += dimX;
       }
       //zDT += planSize;
       //zPD += PlanSize;
   }
   cout <<endl<<endl;
   free_matrix(A,1,3, 1, 3);
   free_vector(eVal,1,3);
   free_matrix(eVec,1,3, 1, 3);

   if (PDfilename == NULL)
	return PD;
   else {
	saveVectorFile(PD,dimX,dimY,dimZ,PDfilename);
	delete []PD;
	return NULL;
   }
}
/*
struct Pt3d * DTVolume::PDofOneTensor()
{

   float **A, *eVal, **eVec, Max;
   eigOfTensor(A,eVal,eVec,3);


}*/

struct Pt3d * DTVolume::PrimaryDirectionDouble(char *PDfilename) //this does not work for 1e-6 DTI, why?
{//  if PDfilename == NULL : return the pointer to the PD vectors
 //  if PDfilename != NULL : save the extracted PD vector to file
 //			     and return NULL (allocated PD space freed)
 
   int XX,YY,ZZ,n, i,iMax,bad;
   double **A, *eVal, **eVec, Max;
   //DTensor *zDT, *yDT, *xDT;
   //struct Pt3d  *xPD,*yPD,*zPD;
   struct Pt3d * PD = NULL;
   double weight,W123,ww,t;

   PD = new struct Pt3d [dimX*dimY*dimZ];
   if (PD == NULL){
	cout <<" DTVolume::PrimaryDirection(): Memory allocation failed"<<endl;
        exit(0);
   }else{
	cout <<" Allocated "<< dimX*dimY*dimZ<<" units"<<endl;
   }

   A = dmatrix (1,3, 1, 3);	// index from 1 ~3
   eVal = dvector(1,3);		// index from 1 ~3
   eVec = dmatrix(1,3, 1, 3);	// index from 1 ~3

   //zDT = dt;
   //zPD = PD;

   cout <<endl;
   n=0; bad=0;
   for (ZZ = 0; ZZ < dimZ ;ZZ++){  //for (ZZ = 0; ZZ < dimZ-1 ;ZZ++)
       cout << "   Processing slice Z = " <<ZZ<<char(13);
       flush(cout);
       //yDT = zDT;
       //yPD = zPD;
       for (YY = 0 ; YY < dimY; YY ++){ //for (YY = 0 ; YY < dimY-1; YY ++){
          //xDT = yDT;
	  //xPD = yPD
	  for (XX = 0; XX < dimX; XX++,n++){ //  for (XX = 0; XX < dimX-1; XX++){
		//if (ZZ < 125) continue;
	     t = 0;
	     for (int m=0;m<6;m++)
		t += dt[n].a[m];
	     if (t == 0){
               A[1][1] =  A[1][2] =  A[1][3] = 0;
    	       A[2][1] =  A[2][2] =  A[2][3] = 0;
    	       A[3][1] =  A[3][2] =  A[3][3] = 0;
	     }else {
               A[1][1] = dt[n].a[0]/t; A[1][2] = dt[n].a[3]/t; A[1][3] = dt[n].a[4]/t;
    	       A[2][1] = A[1][2]     ; A[2][2] = dt[n].a[1]/t; A[2][3] = dt[n].a[5]/t;
    	       A[3][1] = A[1][3]     ; A[3][2] = A[2][3]     ; A[3][3] = dt[n].a[2]/t;
	     }
	     if ((0==eigOfTensorDouble(A,eVal,eVec,3))&&(bad<10)){
		 //cout <<"   bad Tensor["<<n<<"]"<<endl;
		 //printTensorAt(n);
		 //seeMatrix(A,3,3,"A matrix:");
		 bad ++;
		 //if (bad > 10)
		//	exit(0);
	     }
	     W123 = eVal[1]*eVal[1] + eVal[2]*eVal[2] + eVal[3]*eVal[3];
	     if (W123 < ZERO){
		weight = 0.0;
	     }else {
		ww = fabs(eVal[1]*eVal[2])+fabs(eVal[2]*eVal[3])+fabs(eVal[3]*eVal[1]);
		weight =  1.0 - ww/W123;  // could be negative, due to the computer accuracy
		if (weight < 0) 
			weight = 0;
		else 
		    weight = sqrt (weight);
	     }
	     Max = -1;
	     for (i = 1; i<= 3; i++){  // find Primary Direction
	       if (fabs(eVal[i])>Max){
		   Max=fabs(eVal[i]);
		   iMax = i;
	       }
	     }		
	     //get PD vector (weighted)
	     PD[n].x = weight*eVec[1][iMax] + XX;
	     PD[n].y = weight*eVec[2][iMax] + YY;
	     PD[n].z = weight*eVec[3][iMax] + ZZ;

	     //xPD ++;
	     //yDT ++;
	     //n++;
	  }
	  //yDT += dimX;
	  //yPD += dimX;
       }
       //zDT += planSize;
       //zPD += PlanSize;
   }
   cout<<endl<<endl;
   free_dmatrix(A,1,3, 1, 3);
   free_dvector(eVal,1,3);
   free_dmatrix(eVec,1,3, 1, 3);

   if (PDfilename == NULL)
	return PD;
   else {
	saveVectorFile(PD,dimX,dimY,dimZ,PDfilename);
	delete []PD;
	return NULL;
   }
}
void DTVolume::fatterFAinRange(int ax, int ay, int az, 
   float radius,  float percentage, float acc)
{ // this program is to lower FA of the tensors in a neighborhood of (ax,ay,az)
  // with the given radius (thus the neighborhood is a shpere).
  // The "percentage" indicates how much is going to be lowered
  // for example, if the original 3 eigenvalues were 1.0, 0.1, 0.2
  // the lower 20% will be : 1.0, 0.1 + (1.0-0.1)*0.2= 0.28, 0.2+(1.0-0.2)*0.2=0.36
  // a percentage of 100% will make the tensor's fA = 0.0
  // The input acc is a factor to control the percentage, the center FA will be lowered "percentage"
  // The percentage to be lowered at a little far place will be decided according to:
  //      (acc ^ distance ) * percentage
  // where "distance" is the location under consideration to the center point (ax,ay,az)
  // In this way, the simulated lower FA will be very smooth, and merges into
  // its neighborhood

  int ZZ, XX, YY, OK, m;
  DTensor DT1;
  float **A, **B, **C, *eVal, **eVec, Max, t, realpercent, dist;
  struct Pt3d pp;
  double  oFA,nFA, cFA, tcFA; //oldFA, newFA, chang of FA, total of FA
  int ct; //counter

  pp.x=ax; pp.y=ay; pp.z=az;

  A = matrix (1,3, 1, 3);	// index from 1 ~3
  eVal = vector(1,3);		// index from 1 ~3
  eVec = matrix(1,3, 1, 3);	// index from 1 ~3
  B =  matrix(1,3, 1, 3);
  C =  matrix(1,3, 1, 3);

  cout <<"radius/resZ="<<radius/resZ<<endl;
  cout <<"radius/resY="<<radius/resY<<endl;
  cout <<"radius/resX="<<radius/resX<<endl;
  cout <<"tmpX="<<tmpX<<endl;
  cout <<"tmpY="<<tmpY<<endl;
  cout <<"tmpZ="<<tmpZ<<endl;
  //exit (0);
  ct = 0; tcFA = 0;
  for (ZZ = int(az - radius/resZ+0.5); ZZ < (az+radius/resZ) ;ZZ++){
       for (YY = int(ay - radius/resY+0.5) ; YY < (ay+radius/resY); YY ++){ 
	  for (XX = int(ax - radius/resX+0.5); XX < (ax+radius/resX); XX++){ 
	        dist = distance (XX,YY,ZZ, pp);
		if (dist > radius) 
		   continue;
		OK=getTensorAt(XX, YY, ZZ, DT1);
		   //printTensor("\n\nTensor  DT1 original:\n", DT1);
		if (OK == 1){
		   ct ++;
		   realpercent = percentage;
		   for (m = 4; m < dist; m++)
			realpercent *= acc;
		   //realpercent = percentage;
	           t = 0;
	           for (int m=0;m<6;m++)
		       t += DT1.a[m];
		   //cout <<"("<<XX<<","<<YY<<","<<ZZ<<")"<<endl;
		   //cout <<"t="<<t<<endl;
	           if (t == 0){
                       //A[1][1] =  A[1][2] =  A[1][3] = 0;
    	               //A[2][1] =  A[2][2] =  A[2][3] = 0;
    	               //A[3][1] =  A[3][2] =  A[3][3] = 0;
		       continue;
	           }else {
                       A[1][1] = DT1.a[0]/t  ; A[1][2] = DT1.a[3]/t   ; A[1][3] = DT1.a[4]/t;
    	               A[2][1] = A[1][2]     ; A[2][2] = DT1.a[1]/t   ; A[2][3] = DT1.a[5]/t;
    	               A[3][1] = A[1][3]     ; A[3][2] = A[2][3]      ; A[3][3] = DT1.a[2]/t;
	           }
		   //seeMatrix(A, 3,3, " Matrix A:\n");
	           OK=eigOfTensor(A,eVal,eVec,3);
		   oFA = getFA(eVal);
		   //seeMatrix(eVec, 3,3, " Matrix eVec:\n");
		   //seeVector(eVal, 3, " Eigen values:\n");
		   Max= eVal[1];  // find the largest eigenvalue among the 3
		   if (eVal[2]>Max) 
			Max= eVal[2];
		   if (eVal[3]>Max) 
			Max= eVal[3];		   
		   //cout <<"\n"<<"Max = "<<Max<<endl;
    	           //modified eigenvalues to be a new matrix C
		   for (m=1; m<4; m++){
			B[m][m] =  eVal[m] + (Max-eVal[m])*realpercent;
		   }
                   B[1][2] = B[1][3] = 0;
   	           B[2][1] = B[2][3] = 0;
    	           B[3][1] = B[3][2] = 0;
	     	   //B[1][1] =  eVal[1];
	     	   //B[2][2] =  eVal[2];
	     	   //B[3][3] =  eVal[3];
		   //seeMatrix(B, 3,3, " Matrix B:\n");

 		   mulMatrixAB(eVec,B,C,3);
		   mulMatrixABt(C,eVec, A, 3, 3); // A = eVec * B * eVec'

		   //seeMatrix(A, 3,3, " revised Matrix A:\n");

		   //------recover to the original scale -----
		   DT1.a[0] = A[1][1] * t; 
		   DT1.a[1] = A[2][2] * t; 
		   DT1.a[2] = A[3][3] * t; 
		   DT1.a[3] = A[1][2] * t; 
		   DT1.a[4] = A[1][3] * t; 
		   DT1.a[5] = A[2][3] * t; 
		   //printTensor("\n\nTensor DT1:\n", DT1);
	           //setTensorAt(XX, YY, ZZ, DT1.a[0], DT1.a[1],DT1.a[2],DT1.a[3],DT1.a[4],DT1.a[5]);
		   setTensorAt(XX, YY, ZZ, DT1);
		   //cout <<" ~~~~~"<<endl;
	           OK=eigOfTensor(A,eVal,eVec,3);
		   nFA = getFA(eVal);
		   if (oFA == 0 ) cout<< " FA == 0 !!!! Suprise !!!"<<endl;
		   else {
 		     cFA = (oFA - nFA)/oFA; //cFA should >= 0
		     if (cFA <0) cout <<" #########!!! ";
		     cout << " change of FA = "<< cFA *100<< "% "<<endl;
		   }
		   tcFA += cFA;		   
		   //seeMatrix(eVec, 3,3, " Matrix eVec:\n");
		   //seeVector(eVal, 3, " Eigen values:\n");
		}

	  }
	}
  }
  free_matrix(A,1,3, 1, 3);
  free_vector(eVal,1,3);
  free_matrix(eVec,1,3, 1, 3);
  free_matrix(C,1,3, 1, 3);
  free_matrix(B,1,3, 1, 3);

  cout<< "\n\n FINAL REPORT :\n";
  cout <<"   Simulation at ("<<ax<<","<<ay<<","<<az<<")\n";
  cout <<"   Radius = "<<radius<<" mm\n";
  cout <<"   percentage = " <<percentage * 100<<"%; accelerate rate = "<< acc*100 <<"%\n";
  cout <<"  --> The neighborhood has: "<< ct << "  voxels"<<endl;
  if (ct > 0)
     cout <<"  --> FA has been lowered: "<< tcFA /(float)ct * 100.0f <<"% \n\n"<<endl;
  else
     cout <<"  --> FA not defined. \n\n"<<endl;
}
/************************************************
 Interlaced VOLUMELY format for *.d:
*************************************************
   Slice1 Dxx (65536 float) // the 1st component
   Slice2 Dxx (65536 float)
     .......
   SliceN Dxx (65536 float)
   Slice1 Dyy (65536 float)
   Slice2 Dyy (65536 float)  // the 2nd component
     .......

   SliceN Dyy (65536 float)
     ........
     ........
   Slice1 Dyz (65536 float)  // the 6th component
   Slice2 Dyz (65536 float)
     ......
   SliceN Dyz (65536 float)
************************************************/

void DTVolume::interlace2normal()
{// change interlaced order to normal order
 // The original interlaced format:
 // volume of Dxx + volume of Dyy + .... + volume of Dyz
  DTensor *nn= NULL;  //new normal volume
  float * oldVM = (float*)dt;  //interlaced DT volume
  long slice, row, column, ind, i, nnPt;  

  nn = new struct DTensor[dimZ * planSize];

  if (nn == NULL){
     cout << "DTVolume::interlace2normal(): Memory failed"<<endl;
     exit(0);
  }

  //-----show progress-------
  cout <<"\n|";
  for (i=0;i<50;i++)
	cout <<"_";
  cout <<"|"<<endl<<"|";
  flush(cout);
  float idd, seg = dimZ * planSize*6 / 50;
  idd = 0;  
  //--------------------------

  nnPt = 0; // point to the interlaced volume
  for (i = 0; i< 6; i++){
     ind = 0;  // pointed to the normal volume
     for (slice = 0; slice < dimZ; slice ++){   
	//-----show progress--------
        if ( nnPt > idd){
	  cout <<"=";
	  //cout <<char(172);
          flush(cout);
	  idd += seg;
        }
	//-----------------------
        for (row = 0; row < dimY; row ++)
        for (column = 0; column < dimX ; column ++){
          nn[ind].a[i] = oldVM[nnPt];  
          ind ++;
	  nnPt ++;
	}
     }
  }
  cout <<"|\n\n Interlaced ---> normal: completed.\n\n"<<endl;
  
  delete []dt;
  dt = nn;
}
void DTVolume::normal2interlace() // this procedure has not been verified
{// change normal order to interlaced order
 // after this procedure, the dt actually points to an array of "float"
 // although it is of type DTensor*, but not every 6 floats do not 
 // compose to a Tensor
  float *nn= NULL;  //new interlaced DT volume
  DTensor * oldVM = dt;  //normal volume 
  long slice, row, column, ind, i, nnPt;

  nn = new float[dimZ * planSize*6];

  if (nn == NULL){
     cout << "DTVolume::interlace2normal(): Memory failed"<<endl;
     exit(0);
  }
  
  nnPt = 0; // point to the interlaced volume
  for (i = 0; i< 6; i++){
     ind = 0;  // pointed to the normal volume
     for (slice = 0; slice < dimZ; slice ++)
     for (row = 0; row < dimY; row ++)
     for (column = 0; column < dimX ; column ++){
        nn[nnPt] = oldVM[ind].a[i];
        //nn[ind].a[i] = oldVM[nnPt];  
        ind ++;
	nnPt ++;
     }
  }
  delete []dt;
  dt = (DTensor *)nn;
}
/************************************************
 Interlaced SLICELY format for *.d:
*************************************************
   Slice1 Dxx (65536 float)
   Slice1 Dyy (65536 float)
   Slice1 Dzz (65536 float)
   Slice1 Dxy (65536 float)
   Slice1 Dxz (65536 float)
   Slice1 Dyz (65536 float)
   Slice2 Dxx (65536 float)
   Slice2 Dyy (65536 float)
   Slice2 Dzz (65536 float)
   Slice2 Dxy (65536 float)
   Slice2 Dxz (65536 float)
   Slice2 Dyz (65536 float)
     ......
   SliceN Dxx (65536 float)
   SliceN Dyy (65536 float)
   SliceN Dzz (65536 float)
   SliceN Dxy (65536 float)
   SliceN Dxz (65536 float)
   SliceN Dyz (65536 float)
************************************************/
void DTVolume::interlace2normal_slicely()
{// change interlaced order to normal order
 // The original interlaced format:
 // Slice1 of Dxx + slice1 of Dyy + .... + slice1 of Dyz
 // Slice2 of Dxx + slice2 of Dyy + .... + slice2 of Dyz
 // .....
 // SliceN of Dxx + sliceN of Dyy + .... + sliceN of Dyz

  DTensor *nn= NULL;  //new normal volume
  float * oldVM = (float*)dt;  //interlaced DT volume
  long slice, row, column, ind, i, nnPt,indSt;  

  nn = new struct DTensor[dimZ * planSize];

  if (nn == NULL){
     cout << "DTVolume::interlace2normal(): Memory failed"<<endl;
     exit(0);
  }

  //-----show progress-------
  cout <<"\n|";
  for (i=0;i<50;i++)
	cout <<"_";
  cout <<"|"<<endl<<"|";
  flush(cout);
  float idd, seg = dimZ * planSize*6 / 50;
  idd = 0;  
  //--------------------------

  nnPt = 0; // point to the interlaced volume
  indSt = 0;  // point to the normal volume, slice by slice
  for (slice = 0; slice < dimZ; slice ++){   
     for (i = 0; i< 6; i++){
	ind = indSt; //pointer to the normal volume, within a slice
	//-----show progress--------
        if ( nnPt > idd){
	  cout <<"=";
	  //cout <<char(172);
          flush(cout);
	  idd += seg;
        }
	//-----------------------
        for (row = 0; row < dimY; row ++)
        for (column = 0; column < dimX ; column ++){
          nn[ind].a[i] = oldVM[nnPt];  
          ind ++;
	  nnPt ++;
	}
     }
     indSt += planSize;
  }
  cout <<"|\n\n Interlaced ---> normal: completed.\n\n"<<endl;
  
  delete []dt;
  dt = nn;
}
void DTVolume::normal2interlace_slicely()
{
   cout <<"\n\n\n UNDER CONSTRUCTION .......\n\n\n\n"<<endl;
}
void DTVolume::outputs()
{
   cout <<endl;
   cout <<"\n Volume Name: "<<path<<endl;
   cout <<"dimension  = ("<<dimX<<","<<dimY<<","<<dimZ<<");\n";
   cout <<"resolution = ("<<resX<<","<<resY<<","<<resZ<<");\n";
   cout <<"VoxSize    = "<<voxSize<<";  planSize="<<planSize<<"\n\n"<<endl;

}
int DTVolume::fgetTensorAt(float vx,float vy,float vz, DTensor& res)
{
    int row,col,slice,i;
    DTensor *pp, *pp1,*tpp,*tpp1;  
    DTensor v1,v2,v3,v4,v5,v6,v7,v8;
    float tx,ty,tz;

    if ((vx<0)||(vx>dimX-1)) return 0.0f;
    if ((vy<0)||(vy>dimY-1)) return 0.0f;
    if ((vz<0)||(vz>dimZ-1)) return 0.0f;

    col = int(vx);
    row = int(vy);
    slice=int(vz);
   
    pp = dt + planSize*slice + row*dimX + col;
    if (slice == dimZ-1) pp1 = pp;
    else pp1 = pp + planSize;

    v1 = DTensor(*pp);
    v5 = DTensor(*pp1);

    if (col == dimX-1){
	v2 = v1;
	v6 = v5;
    } else {
        v2 = DTensor(*(pp+1));
    	v6 = DTensor(*(pp1+1));
    }
    if (row == dimY-1) {
	v3 = v1;tpp=pp;
	v7 = v5;tpp1=pp1;
    } else {
    	tpp=pp+dimX;   v3 = DTensor(*(tpp));
    	tpp1=pp1+dimX; v7 = DTensor(*(tpp1));
    }
    if (col == dimX-1){
	v4 = v3;
    	v8 = v7;
    } else {
	v4 = DTensor(*(tpp+1));
    	v8 = DTensor(*(tpp1+1));
    } 
    
    tx = vx - col; ty = vy - row; tz = vz - slice;
    for (i=0; i<6; i++){
      v1.a[i] = tx*v2.a[i]+(1.0f-tx)*v1.a[i];
      v3.a[i] = tx*v4.a[i]+(1.0f-tx)*v3.a[i];

      v5.a[i] = tx*v6.a[i]+(1.0f-tx)*v5.a[i];
      v7.a[i] = tx*v8.a[i]+(1.0f-tx)*v7.a[i];

      v1.a[i] = ty*v3.a[i]+(1.0f-ty)*v1.a[i];
      v5.a[i] = ty*v7.a[i]+(1.0f-ty)*v5.a[i];

      res.a[i] = tz*v5.a[i]+(1.0f-tz)*v1.a[i];
    }

    return 1;
}
void DTVolume::flipVolume(int which)
{// which==0:  along X dir
 // which==1:  along Y dir
 // which==2:  along Z dir

  long slice, row, column, hh;
  DTensor t1,t2;

  cout<<endl;
  switch (which){
  case 0: // X-dir
    hh = dimX /2 ;cout <<" \nNow flip DT volume along X-dir"<<endl;
    for (slice = 0; slice < dimZ; slice ++){
      cout <<" Processing Z= "<<slice<<"   "<<char(13);
      flush(cout);
      for (row = 0; row < dimY; row ++)
      for (column = 0; column < hh ; column ++){

	getTensorAt(column,row,slice,t1);
	getTensorAt(dimX-1-column,row,slice,t2);

	setTensorAt(column,row,slice,t2);
	setTensorAt(dimX-1-column,row,slice,t1);
      }
    }
    cout<<endl;
    break;
  case 1: //Y-dir
    hh = dimY /2 ;cout <<" \nNow flip DT volume along Y-dir"<<endl;
    for (slice = 0; slice < dimZ; slice ++){
      cout <<" Processing Z= "<<slice<<"   "<<char(13);
      flush(cout);
      for (row = 0; row < hh; row ++)
      for (column = 0; column < dimX ; column ++){

	getTensorAt(column,row,slice,t1);
	getTensorAt(column,dimY-1-row,slice,t2);

	setTensorAt(column,row,slice,t2);
	setTensorAt(column,dimY-1-row,slice,t1);
      }
    }
    cout<<endl;
    break;
  case 2: //Z-dir
    hh = dimZ/2 ;cout <<" \nNow flip DT volume along Z-dir"<<endl;
    for (slice = 0; slice < hh; slice ++){
      cout <<" Processing Z= "<<slice<<"   "<<char(13);
      flush(cout);
      for (row = 0; row < dimY; row ++)
      for (column = 0; column < dimX ; column ++){
	getTensorAt(column,row,slice,t1);
	getTensorAt(column,row,dimZ-1-slice,t2);

	setTensorAt(column,row,slice,t2);
	setTensorAt(column,row,dimZ-1-slice,t1);
      }
    }
    cout<<endl;
    break;
  }
}
