#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <iostream.h>
//#include <strstream.h>
#include <strstream>
#include <fstream.h>

//#include <task.h>    // for multi-cpus
//#include <ulocks.h>  // for multi-cpus

#include "mVolume.h"
#include "VoxMap.h"
#include "DTVolume.h"
#include "Ellipsoid.h"

struct Pt3d subpp(struct Pt3d P1,struct Pt3d P2)
{ 
  struct Pt3d  P;
  P.x = P2.x - P1.x;
  P.y = P2.y - P1.y;
  P.z = P2.z - P1.z;
  return P;
}

mVolume::mVolume()
  :dimX(0),dimY(0),dimZ(0),mv(NULL)
{
}

mVolume::mVolume(int xx, int yy, int zz)
  :dimX(xx),dimY(yy),dimZ(zz)
{
   getSpace();
}

mVolume::~mVolume()
{
  if (mv!=NULL)
	delete []mv;
}

void mVolume::getSpace()
{
  if (mv!=NULL)
	delete []mv;   
  mv = new Matrix[dimX*dimY*dimZ];
  if (mv == NULL){
     cout << "mVolume::getSpace(): memory allocation failed";
     exit(0);
  }else 
     cout <<"mVolume::getSpace(): Memory allocated."<<endl;
}

void mVolume::getTransMatrixOld(VoxMap* &wV)
{//wV must not be a displacement
 // it should be a vector field defining the destination coordinates
 // at each voxel. i.e. don't do wV->reorganize()
    int k,XX,YY,ZZ,kx,ky,kz,ii;//, sss=0, def=0;
    struct Pt3d p,p0;
    float **A, **B, **C, **V, **U;
    float *w;

    //getSpace();
    A = matrix(1,3, 1, 26);
    B = matrix(1,3, 1, 26);
    C = matrix(1,3, 1,3);
    V = matrix(1,3, 1,3);
    U = matrix(1,3, 1,3);
    w = vector(1,3);
    k=0;
    for (ZZ = 0; ZZ < dimZ; ZZ ++){
      for (YY = 0; YY < dimY; YY ++){
        for (XX = 0; XX < dimX; XX ++){
	  p0= wV->getVectorAt(XX  ,YY  ,ZZ);
	  if (  p0.x == UNDEFINED_PT 
	     && p0.y == UNDEFINED_PT 
	     && p0.z == UNDEFINED_PT  ){
		k++;
		continue;
	  }
 	  /*******debug ****
	  else if (  p0.x != 0 
	     && p0.y != 0 
	     && p0.z != 0) {
		cout <<" F#F"<<endl; 
		def ++;
	  }
	  *******debug *****/

	   //def =0;
	  ii = 1;
 	  for (kz = -1; kz < 2; kz++)
 	  for (ky = 1; ky > -2; ky--)
 	  for (kx = -1; kx < 2; kx++){
		if ((kx ==0 )&&(ky==0)&&(kz==0))
		   continue;
		p = wV->getVectorAt(XX+kx,YY+ky,ZZ+kz);	
		  //cout <<" -%%-" <<endl;
		if (  p.x == UNDEFINED_PT 
		   && p.y == UNDEFINED_PT 
		   && p.z == UNDEFINED_PT  ){
  	    	    B[1][ii] = B[2][ii] = B[3][ii] = 0.0f;
   	    	    A[1][ii] = A[2][ii] = A[3][ii] = 0.0f;
		} else {
		    p = subpp(p0,p);
		      //if ((p.x != 0) &&(p.x != 0) &&(p.x != 0))
		//	def ++;
   	    	    B[1][ii] = kx;  B[2][ii] = ky;  B[3][ii] = kz;
   	    	    A[1][ii] = p.x; A[2][ii] = p.y; A[3][ii] = p.z;
		}
		ii++;
	  }
	  mulMatrixABt(A,B,C, 3, 26); // C = A*B
	   /********
	   if (def >3 && sss<10){
		seeMatrix(A,3,26, "matrix AA:");
		seeMatrix(B,3,26, "matrix BB:");
		seeMatrix(C,3,3, "matrix CC:");
	
	   }
	   ****/
	  svdcmp(C, 3, 3, w, V);  //C = C * [w]*V
	  mulMatrixABt(C,V, U, 3, 3);// U = 3x3
	   /*****debug***
	   if ((def >5  )&& (sss<10 )){
		seeMatrix(C,3,3, "matrix C of C*[w]*V:");
		seeMatrix(V,3,3, "matrix V:");
		seeMatrix(U,3,3, "matrix U(final):");	
		cout <<endl <<"SSS="<<sss<<"================"<<endl;
	  	//sss ++;
	   }*/
	   /******debug**/
	  mv[k].setValue(U[1][1], U[1][2], U[1][3], 0.0f,
			 U[2][1], U[2][2], U[2][3], 0.0f,
			 U[3][1], U[3][2], U[3][3], 0.0f,
			 0.0f,    0.0f,    0.0f,    1.0f);
	  /*if ((def >5  )&& (sss<10 )){
		//seeMatrix(C,3,3, "matrix C of C*[w]*V:");
		cout << "Now mv["<<k<<"]:"<<endl;
	        mv[k].outputs();
		sss ++;
	  }*/
	  k++;
        }
      }
   }
   free_matrix(A,1,3, 1, 26);
   free_matrix(B,1,3, 1, 26);
   free_matrix(C,1,3, 1, 3);
   free_matrix(V,1,3, 1, 3);
   free_matrix(U,1,3, 1, 3);
   free_vector(w,1,3);
   //cout << "def ="<<def << ", percentage = "<<(float)def/(dimX*dimY*dimZ)*100.0<<"%."<<endl;
}

Matrix & mVolume::getMatrix(int ind)
{
   return mv[ind];
}

Matrix & mVolume::getMatrix(int xx,int yy, int zz)
{
  return getMatrix(zz*dimX*dimY+yy*dimX+xx);
}

void mVolume::getTransMatrix(VoxMap* &wV, DTVolume* TensorVol, VoxMap* &PD,VoxMap* &PD2, DTVolume* &optDTI,
	int useOptimized)
{//wV must not be a displacement
 // it should be a vector field defining the destination coordinates
 // at each voxel. i.e. don't do wV->reorganize()
    void mulMatrixABtt(float **A, float **B, float **C, int row, int col);

    int k,XX,YY,ZZ; //, sss=0, def=0;
    //struct Pt3d p,p0;
    float **A, **B, **C, **V, **U1, **U2, **UU;
    float *w;
    Ellipsoid optE;
    int neib; // bounding box of E
    int r,nAB; //dimension of **A & **B, so they will be of 3*nAB (A & B start from 1 to nAB)
    int records[1500]; //enough to hold marks for 7*7*7=343 neighbour; 11*11*11=1331
    int nCPU; //to know how many CPUs are being used
    int myID; // identify which CPU task
    int Z1, Z2; // start and end of Z slice because of multi-CPU tasks
    float **F, *eVal, **eVec, **T; // variables for subroutine: weightedShapeAt()
    int bad; // for debugging   
    long underDefined=0; // see how many places (nAB<3) is true
    //double dtScale; // scale to adjust DT to be around 1.0
    float tvalues[6]; //tensor values: 6 components
    //DTVolume optDTI;

    //optDTI.newVolume(dimX,dimY,dimZ);

    //nCPU= m_get_numprocs(); //to know how many CPUs are being used
    //myID = m_get_myid();

    nCPU=1; myID= 0;
    ZZ = dimZ / nCPU;  
	cout <<"A block contains "<<ZZ<<" Z-slice"<<endl;
    k  = ZZ * (dimX*dimY) * myID;  // k is the volume memory index
    Z1 = ZZ * myID;
    if (myID == nCPU - 1){ // the last CPU
	Z2 = dimZ;
    }else{
    	Z2 = Z1 + ZZ;
    }
    cout <<"\n CPU #"<<myID<<" working on Z slice No. "<<Z1<<".... "<<Z2-1<<endl;
    //if (myID!=2)
	  // return;
    //getSpace();
    //A = matrix(1,3, 1, 26);
    //B = matrix(1,3, 1, 26);
    A = matrix(1,3, 1, 120);
    B = matrix(1,3, 1, 120);
    C = matrix(1,3, 1,3);
    V = matrix(1,3, 1,3);
    U1 = matrix(1,3, 1,3);
    U2 = matrix(1,3, 1,3);
    UU = matrix(1,3, 1,3);
    w = vector(1,3);
    //---- for use of a subroutine: weightedShapeAt()
      F = matrix (1,3, 1, 3);
      T = matrix (1,3, 1, 3);
      eVal = vector(1,3);
      eVec = matrix(1,3, 1, 3);
    //------------------------------------
    cout <<endl<<endl;
    //k=0;
    bad =0; cout<<"\n\n"<<endl;
    for (ZZ = Z1; ZZ < Z2; ZZ ++){
      cout <<"   processing  ... Z = "<<ZZ<<"; k= "<<k<<"    "<<char(13);flush(cout);
      //if (ZZ > 8 ) exit(0);
      for (YY = 0; YY < dimY; YY ++){
	//cout <<"       processing  ... Y = "<<YY<<endl;
        for (XX = 0; XX < dimX; XX ++, k++){
	  //if (YY==195 && XX == 134)
		//cout <<"    XX="<<XX<<endl;
	  //if (XX!=126 || YY!=210 || ZZ != 0)
		//continue;
	    //cout <<"XYZ=("<<XX<<","<<YY<<","<<ZZ<<")"<<endl;
	  //if (ZZ<125)
		//continue;	
	  //r = weightedShapeAt(XX,YY,ZZ,wV,TensorVol , E, neib);//return value in E & neib

	  //return a neighborhood defined in records & neib
	  //r = weightedShapeAt(XX,YY,ZZ,wV,TensorVol , records, neib, F, eVal, eVec); //, T);	  


	  //return a neighborhood defined in records & neib; and an optimized Ellipsoid in optE
	  r = weightedShapeAt_withOptE(XX,YY,ZZ,wV,TensorVol , records, neib, F, eVal, eVec, optE);
	  if ((r == 0 ) || ( neib <1))
		continue;
  	    /*if (XX==2 && YY==2 && ZZ == 0){
		cout <<"XYZ=("<<XX<<","<<YY<<","<<ZZ<<")"<<endl;
		optE.outputs();
	    }*/
	  //rv  optE.getTensor(tvalues);//turn the optimized ellipsoid to a tensor
	    //if (XX==2 && YY==2 && ZZ == 0){
		//for (int h=0; h< 6; h++)
		  // cout <<tvalues[h]<<", ";
		//cout<<endl<<"~~~~~~~~~~~~~"<<endl;;
	    //}
	  //rv optDTI->setTensorAt(XX,YY,ZZ, tvalues ); //float tvalues[6]
	    //continue ; // this line is for debug only, MUST BE REMOVED !!!!!
  	  //if (XX==69 && YY==184 && ZZ == 7)
		//cout <<"### r = "<<r<<"; neib= "<<neib<<endl;
	  //getProcrusteanAB(XX,YY,ZZ, wV, PD, E, neib, A, B, nAB); //return value in A, B & nAB
	  //return value in A, B & nAB based on neighborhood marks in 'records'
	  getProcrusteanAB(XX,YY,ZZ, wV, PD, records, neib, A, B, nAB);
	  //if (XX==124 && YY ==0 && ZZ == 0)
		//cout <<"\n\n --$$$$$$$$$$$$ ----(124,0,0)-- nAB ="<<nAB<<" \n\n\n"<<endl; 

	  //if (( nAB<3)&&(ZZ>6)&&(ZZ<dimZ-6)){
	  if (nAB<3){
		underDefined ++;
		//cout <<" OHHH!!!! mVolume::getTransMatrix(): nAB ="<<nAB<<" (< 3)"<<endl; 
		//cout <<"XYZ=("<<XX<<","<<YY<<","<<ZZ<<")"<<endl;
			/*int ind = 0;
  			for (int kz = -neib; kz <= neib; kz++)
			for (int ky = -neib; ky <= neib; ky++)
			for (int kx = -neib; kx <= neib; kx++, ind ++){
				if (records[ind]==1)
			    	  cout <<"  -> ("<<kx<<","<<ky<<","<<kz<<")=1"<<endl;
			}*/
	  }
	  //mulMatrixABt(A,B,C, 3, 26); // C = A*B
	  mulMatrixABt(A,B,C, 3, nAB); // C = A*B

         /********
		seeMatrix(A,3,nAB, "matrix AA:");
		seeMatrix(B,3,nAB, "matrix BB:");
		seeMatrix(C,3,3, "matrix CC:");
	   ****/
	   r=svdcmp(C, 3, 3, w, V);  //C = C * [w]*V
	   if (r==0){// bad result
		float ** CC = matrix (1,3,1,3);
		cout <<"bad matrix @ ("<<XX<<","<<YY<<","<<ZZ<<")"<<endl;
		seeMatrix(A,3,nAB,"A :");
		seeMatrix(B,3,nAB,"B :");
	        mulMatrixABt(A,B,CC, 3, nAB); // C = A*B		
		seeMatrix(CC,3,3,"CC=A*B':");
		seeMatrix(C,3,3,"C of CC=C*W*V':");
		seeVector(w,3,"W of CC=C*W*V':");
		seeMatrix(V,3,3,"V of CC=C*W*V':");

		free_matrix(CC,1,3,1,3);
		bad ++;
		if (bad > 10) exit(0);
	   }
	   mulMatrixABt(C,V, U1, 3, 3);// U = 3x3

	   /*****debug***
		seeMatrix(C,3,3, "matrix C of C*[w]*V:");
		seeMatrix(V,3,3, "matrix V:");
		seeMatrix(U,3,3, "matrix U(final):");	
	   ******debug****/
	  //----------calculate PD2's rotation --------Begin-----------
	  if (PD2 != NULL){
           
	    Vector PD_b4warp;
	    if (useOptimized == 1) { // use optimized PD
		PD_b4warp=optE.getPrimaryDir();
	    } else{
		struct Pt3d ppp;
		int r1 = PD->getVectorAt(XX,YY,ZZ, ppp); //get the PD of the tensor
	  	if (r1 == 0) {
			cout<<" Impossible !!!! because r1==0 only when current location out of volume"<<endl;
			continue;
	  	}
		PD_b4warp.setXYZ(ppp.x-XX,ppp.y-YY,ppp.z-ZZ);
	    }
	    //...... warp : PD_b4warp --> PDwarped
	    Matrix trans; 
	    trans.setValue(U1[1][1], U1[1][2], U1[1][3], 0.0f,
			 U1[2][1], U1[2][2], U1[2][3], 0.0f,
			 U1[3][1], U1[3][2], U1[3][3], 0.0f,
			 0.0f,     0.0f,     0.0f,     1.0f);
	    Vector  PDwarped = trans * PD_b4warp ;
	    float q1,q2,q3;
	    PDwarped.getXYZ(q1,q2,q3);
	    if (fabs(q1) < 0.0001 && fabs(q2) < 0.0001 && fabs(q3) < 0.0001){ //PD does not exist, so use U1 enough
	      mv[k].setValue(U1[1][1], U1[1][2], U1[1][3], 0.0f,
			 U1[2][1], U1[2][2], U1[2][3], 0.0f,
			 U1[3][1], U1[3][2], U1[3][3], 0.0f,
			 0.0f,     0.0f,     0.0f,     1.0f);

	    }else {
	      int oldnAB = nAB;
	      getProcrusteanAB4PD2(XX,YY,ZZ, wV, PDwarped, PD2, records, neib, A, B, nAB);
		/* do not delete the following "if(oldnAB != nAB){...} */
		/*
		if (oldnAB != nAB){  //this only happens at the boundary regions of the volume
				     // where the neighborhood contains outside of the volume
				     // so that the sample out of the volume is omitted
				     // some places can sample PD1 but cannot sample PD2, or vice versa
		   cout <<" $";
		  //cout <<" ---!!! --- oldnAB = "<<oldnAB<<"; nAB ="<<nAB<<" ("<<XX<<","<<YY<<","<<ZZ<<")"<<endl; 
		  //exit(0);
	 	}*/

	      //if (( nAB<3)&&(ZZ>6)&&(ZZ<dimZ-6)){
	      if (nAB<3){
		underDefined ++;
		//cout <<" OHHH!!!! mVolume::getTransMatrix(): nAB ="<<nAB<<" (< 3)"<<endl; 
		//cout <<"XYZ=("<<XX<<","<<YY<<","<<ZZ<<")"<<endl;
			/*int ind = 0;
  			for (int kz = -neib; kz <= neib; kz++)
			for (int ky = -neib; ky <= neib; ky++)
			for (int kx = -neib; kx <= neib; kx++, ind ++){
				if (records[ind]==1)
			    	  cout <<"  -> ("<<kx<<","<<ky<<","<<kz<<")=1"<<endl;
			}*/
	      }
	      mulMatrixABt(A,B,C, 3, nAB); // C = A*B

              /********
		seeMatrix(A,3,nAB, "matrix AA:");
		seeMatrix(B,3,nAB, "matrix BB:");
		seeMatrix(C,3,3, "matrix CC:");
	      ****/
	      r=svdcmp(C, 3, 3, w, V);  //C = C * [w]*V
	      if (r==0){// bad result
		float ** CC = matrix (1,3,1,3);
		cout <<"bad matrix @ ("<<XX<<","<<YY<<","<<ZZ<<")"<<endl;
		seeMatrix(A,3,nAB,"A :");
		seeMatrix(B,3,nAB,"B :");
	        mulMatrixABt(A,B,CC, 3, nAB); // C = A*B		
		seeMatrix(CC,3,3,"CC=A*B':");
		seeMatrix(C,3,3,"C of CC=C*W*V':");
		seeVector(w,3,"W of CC=C*W*V':");
		seeMatrix(V,3,3,"V of CC=C*W*V':");

		free_matrix(CC,1,3,1,3);
		bad ++;
		if (bad > 10) exit(0);
	      }
	      mulMatrixABt(C,V, U2, 3, 3);// U = 3x3
	      mulMatrixAB (U2, U1, UU, 3);  //  UU = U2 * U1 because when reorient a tensor D: newD = UU * D * UU'
	      // resulting transformation matrix into the mv (Matrix Volume)
	      mv[k].setValue(UU[1][1], UU[1][2], UU[1][3], 0.0f,
			   UU[2][1], UU[2][2], UU[2][3], 0.0f,
			   UU[3][1], UU[3][2], UU[3][3], 0.0f,
			   0.0f,     0.0f,     0.0f,     1.0f);
	    }
	  }//----------calculate PD2's rotation --------End  -----------
          else { // no PD2 ---- use PD's rotation as transformation ----
	    // resulting transformation matrix into the mv (Matrix Volume)
	    mv[k].setValue(U1[1][1], U1[1][2], U1[1][3], 0.0f,
			 U1[2][1], U1[2][2], U1[2][3], 0.0f,
			 U1[3][1], U1[3][2], U1[3][3], 0.0f,
			 0.0f,     0.0f,     0.0f,     1.0f);
	  }

	  /*******debug********
	  if ((def >5  )&& (sss<10 )){
		//seeMatrix(C,3,3, "matrix C of C*[w]*V:");
		cout << "Now mv["<<k<<"]:"<<endl;
	        mv[k].outputs();
		sss ++;
	  }
	  ***** debug *******/
        }
      }      
   }
    cout <<"\n\n ***** Processor ("<<myID<<"): I am here !!! *****"<<endl;
   //free_matrix(A,1,3, 1, 26);
   //free_matrix(B,1,3, 1, 26);
   free_matrix(A,1,3, 1, 120);
   free_matrix(B,1,3, 1, 120);
   free_matrix(C,1,3, 1, 3);
   free_matrix(V,1,3, 1, 3);
   free_matrix(U1,1,3, 1, 3);
   free_matrix(U2,1,3, 1, 3);
   free_matrix(UU,1,3, 1, 3);
   free_vector(w,1,3);
   //--- free the memory for the use of a subroutine 
    free_matrix(F,1,3, 1, 3);
    free_matrix(T,1,3, 1, 3);
    free_vector(eVal,1,3);
    free_matrix(eVec,1,3, 1, 3);
   //-------------------
   //cout << "def ="<<def << ", percentage = "<<(float)def/(dimX*dimY*dimZ)*100.0<<"%."<<endl;
   if (dimX*dimY*dimZ >0)
        cout << " (Processor :"<<myID<<") "<<underDefined <<" voxels under constrained with Procrustean Analysis " \
	<< (float) underDefined / (dimX*dimY*dimZ) *100.0f <<"%.\n\n"<<endl;
   else cout<< "?????!!!! dimX*dimY*dimZ<=0 !!!\n\n"<<endl;
   //optDTI.saveTensorData("optimized.d");
}

void mVolume::getProcrusteanAB(int col, int row, int slice, VoxMap* &wV, VoxMap* &PD,
			int *records,/* Ellipsoid E, */ int neib, float **A, float **B,int & nAB)
{ //(col,row,slice) : local position
  // wV		    : the transformation displacement field
  // PD		    : the extracted PD of the DTI to be warped
  // E		    : the Ellipsoid shape within which samles are to be taken
  // neib	    : neighbor scope, the -neib +(local) to neib +(local) defines the scanning area
  // **A  	    : A of the Procrustean A = U*B
  // **B  	    : B of the Procrustean A = U*B
  // nAB	    : column dimension of **B, B is of (3 x nB)
  int kx, ky, kz, r1,r2,ind;
  struct Pt3d dtpd0, vf0,vf1, nvf;
  float x1,y1,z1,x2,y2,z2;  
  int ndt=0;
  //int flag; //for debug

  nAB=1;
  ind = 0;
	//cout <<"--PAB--"<<endl;
  /*
  if ((col == 124)&&(row==0)&&(slice==0)){
	cout << "---PAB----: neib="<<neib<<endl;
	flag = 1;
  }else 
	flag = 0;*/
  
  for (kz = -neib; kz <= neib; kz++)
  for (ky = -neib; ky <= neib; ky++)
  for (kx = -neib; kx <= neib; kx++, ind ++){
    	  //cout <<records[ind];

    //if (E.inside(col+kx,row+ky,slice+kz)){
    if (records[ind]){
	/*if (flag){
	   cout <<"==========ndt==="<<ndt<<endl;
	   cout <<"  ---OOO  AA:  XYZ="<<col<<","<<row<<","<<slice<<endl;
	   cout <<"  ---OOO  BB:  neib="<<kx<<","<<ky<<","<<kz<<endl;
	   cout <<"  ---OOO  CC:  XYZ+neib="<<col+kx<<","<<row+ky<<","<<slice+kz<<endl;
	   cout << "records[ind="<<ind<<"]="<<records[ind]<<";   nAB="<<nAB<<endl;
 
	}*/
        r1 = PD->getVectorAt(col+kx,row+ky,slice+kz, dtpd0); //get the PD of the tensor
	    //if (flag) cout <<"  ---AAA: (getPD) XYZ="<<col+kx<<","<<row+ky<<","<<slice+kz<<";  r1= "<<r1<<endl;
  	if (r1 == 0) {
	  //  cout << " r1 = 0" << endl; // rv
	  //	  if (flag) cout <<"  --- PD out: PD==0) XYZ="<<col+kx<<","<<row+ky<<","<<slice+kz<<";  r1 == 0 "<<endl; // rv reoved comment
	   continue;
	}
	   ndt ++;
        dtpd0.x -= col+kx; dtpd0.y -= row+ky; dtpd0.z -= slice+kz;
	if ((dtpd0.x==0) &&(dtpd0.y==0)&&( dtpd0.z==0)){
	   //if (flag) cout <<"  PD is ZERO !!"<<endl;
	   continue;
 	}
	Vector oVec(dtpd0.x,dtpd0.y, dtpd0.z); //old vector for B component of Procrustean's A=U*B
	oVec.norm(); 
	oVec.getXYZ(x1,y1,z1); //original unit PD
	   //if (flag){
		//cout <<"  AA--BB: oVec Unit=("<<x1<<","<<y1<<","<<z1<<endl;
	   //}
	//see how the PD is to be moved
        r1 = wV->getVectorAtAnyPosition(col+kx-0.5*x1,row+ky-0.5*y1,slice+kz-0.5*z1, vf0);
        r2 = wV->getVectorAtAnyPosition(col+kx+0.5*x1,row+ky+0.5*y1,slice+kz+0.5*z1, vf1);
	if ((r1 == 0)||(r2 == 0)) {
	  // cout << " r1 & r2 both zero " << endl; // rv comment
	     //if (flag) cout <<"  BB#BB: X1Y1Z1="<<x1<<","<<y1<<","<<z1<<endl;
	     //if (flag) cout <<"  BB#BC: XYZ="<<col+kx<<","<<row+ky<<","<<slice+kz<<";  r1 == "<<r1<<"  r2=="<<r2<<endl;
	    continue;
	}
	nvf = subpp(vf0,vf1);
	  /*if (flag){
	   cout <<"  ---  neib="<<kx<<","<<ky<<","<<kz<<endl;
	   cout <<"       ind = "<<ind<<endl;
 	   cout <<"      oVec Unit=("<<x1<<","<<y1<<","<<z1<<")"<<endl;
	   cout <<"          dtpd0="<<dtpd0.x<<","<<dtpd0.y<<","<<dtpd0.z<<endl;
	   cout <<"  nVec:  vf0="<<vf0.x<<","<<vf0.y<<","<<vf0.z<<endl;
	   cout <<"         vf1="<<vf1.x<<","<<vf1.y<<","<<vf1.z<<endl;
	   cout <<"         nvf="<<nvf.x<<","<<nvf.y<<","<<nvf.z<<endl;
	  } */
	Vector nVec(nvf.x,nvf.y,nvf.z); //new vector for A component of Procrustean's A=U*B
	if (nvf.x == 0 && nvf.y == 0 && nvf.z == 0){// current sample vecto == 0 
	   continue; // current sample vecto == 0 
	   cout << " (1) nvf === 0 !!!"<<endl;
	   cout <<"    XYZ="<<col<<","<<row<<","<<slice<<endl;
	   cout <<"    neib="<<kx<<","<<ky<<","<<kz<<endl;
	   cout <<"    XYZ+neib="<<col+kx<<","<<row+ky<<","<<slice+kz<<endl;
	   cout <<"    x1y1z1="<<x1<<","<<y1<<","<<z1<<endl;
	   cout <<"    vf0 @"<<col+kx-0.5*x1<<","<<row+ky-0.5*y1<<","<<slice+kz-0.5*z1<<endl;
	   cout <<"    vf1 @"<<col+kx+0.5*x1<<","<<row+ky+0.5*y1<<","<<slice+kz+0.5*z1<<endl;

	   cout <<"  nVec:  vf0="<<vf0.x<<","<<vf0.y<<","<<vf0.z<<endl;
	   cout <<"         vf1="<<vf1.x<<","<<vf1.y<<","<<vf1.z<<endl;
	   cout <<" ------(-)----"<<endl;
	}
  	nVec.norm(); //normailzed to unit so that transformation is pure rotation without stretching
	nVec.getXYZ(x2,y2,z2); //value after unit PD transformed 
	  /*if (flag){		
	     cout <<"      nVec Unit=("<<x2<<","<<y2<<","<<z2<<")"<<endl;
	     cout <<"      matrix index (nAB)="<<nAB<<endl;
	  }*/
 	B[1][nAB] = x1;  B[2][nAB] = y1;  B[3][nAB] = z1;
   	A[1][nAB] = x2;  A[2][nAB] = y2;  A[3][nAB] = z2;

  	nAB++;
    }
  } //cout <<"--###--"<<endl;
  nAB --; 
  // A[][1] & B[][1] are the 1st defined columns
  // A[][nAB] & B[][nAB] are the last defined columns
    //if (flag) cout <<" -- number in E = "<<ndt <<endl<<endl;
}
void mVolume::getProcrusteanAB4PD2(int col, int row, int slice, VoxMap* &wV, Vector warpedPD1, VoxMap* &PD2,
			int *records,/* Ellipsoid E, */ int neib, float **A, float **B,int & nAB)
{ //(col,row,slice) : local position
  // wV		    : the transformation displacement field
  // PD2	    : the extracted PD2 of the DTI to be warped
  // E		    : the Ellipsoid shape within which samples are to be taken
  // neib	    : neighbor scope, the -neib +(local) to neib +(local) defines the scanning area
  // **A  	    : A of the Procrustean A = U*B
  // **B  	    : B of the Procrustean A = U*B
  // nAB	    : column dimension of **B, B is of (3 x nB)
  int kx, ky, kz, r1,r2,ind;
  struct Pt3d dtpd0, vf0,vf1, nvf;
  float x1,y1,z1,x2,y2,z2;  
  int ndt=0;
  //int flag; //for debug
  Plane pl(fPoint(0,0,0), warpedPD1);
    //pl.outputs();
  nAB=1;
  ind = 0;
	//cout <<"--PAB--"<<endl;
  /*
  if ((col == 124)&&(row==0)&&(slice==0)){
	cout << "---PAB----: neib="<<neib<<endl;
	flag = 1;
  }else 
	flag = 0;
  */
  for (kz = -neib; kz <= neib; kz++)
  for (ky = -neib; ky <= neib; ky++)
  for (kx = -neib; kx <= neib; kx++, ind ++){
    	  //cout <<records[ind];

    //if (E.inside(col+kx,row+ky,slice+kz)){
    if (records[ind]){
	/*if (flag){
	   cout <<"==========ndt 4 PD 2==="<<ndt<<endl;
	   cout <<"  ---OOO  AA:  XYZ="<<col<<","<<row<<","<<slice<<endl;
	   cout <<"  ---OOO  BB:  neib="<<kx<<","<<ky<<","<<kz<<endl;
	   cout <<"  ---OOO  CC:  XYZ+neib="<<col+kx<<","<<row+ky<<","<<slice+kz<<endl;
	   cout << "records[ind="<<ind<<"]="<<records[ind]<<";   nAB="<<nAB<<endl;
 
	}*/
        r1 = PD2->getVectorAt(col+kx,row+ky,slice+kz, dtpd0); //get the PD2 of the tensor
	    //if (flag) cout <<"  ---AAA:  XYZ="<<col+kx<<","<<row+ky<<","<<slice+kz<<";  r1= "<<r1<<endl;
  	if (r1 == 0) {
	    //if (flag) cout <<"  ---AAA: PD2==0) XYZ="<<col+kx<<","<<row+ky<<","<<slice+kz<<";  r1 == 0 "<<endl;
	   continue;
	}
	   ndt ++;
        dtpd0.x -= col+kx; dtpd0.y -= row+ky; dtpd0.z -= slice+kz;
	if (fabs(dtpd0.x)< 0.0001 && fabs(dtpd0.y) <0.0001&& fabs(dtpd0.z)<0.0001){
	   //if (flag) cout <<"  PD2 is ZERO !!"<<endl;
	   continue;
 	}
	Vector oVec(dtpd0.x,dtpd0.y, dtpd0.z); //old vector for B component of Procrustean's A=U*B
	oVec.norm(); 
	oVec.getXYZ(x1,y1,z1); //normalized original PD2 sample

        float xx1, yy1, zz1;
	Vector tVec = pl.projection(oVec);
	float xt1,yt1,zt1; 
	tVec.getXYZ(xt1,yt1,zt1); // warped sample's projection on the plane perpendicular to the warped PD1
	if (fabs(xt1)< 0.0001 && fabs(yt1)< 0.0001 && fabs(zt1)< 0.0001){
	   //cout << " (1) tVec === 0 !!!"<<endl;
	   continue;
	}
	tVec.norm();    //original unit PD2's projection on the plane perpendicular to warped PD1
	tVec.getXYZ(xx1,yy1,zz1);  //xx1, yy1, zz1 form a column of B in A=U*B
	if ((xx1==0) &&(yy1==0)&&(zz1==0)){
	   continue;
 	}
	   //if (flag){
		//cout <<"  AA--BB: oVec Unit=("<<x1<<","<<y1<<","<<z1<<endl;
	   //}
	//see how the PD2 is to be moved
        r1 = wV->getVectorAtAnyPosition(col+kx-0.5*x1,row+ky-0.5*y1,slice+kz-0.5*z1, vf0);
        r2 = wV->getVectorAtAnyPosition(col+kx+0.5*x1,row+ky+0.5*y1,slice+kz+0.5*z1, vf1);
	if ((r1 == 0)||(r2 == 0)) {
	     //if (flag) cout <<"  BB#BB: X1Y1Z1="<<x1<<","<<y1<<","<<z1<<endl;
	     //if (flag) cout <<"  BB#BC: XYZ="<<col+kx<<","<<row+ky<<","<<slice+kz<<";  r1 == "<<r1<<"  r2=="<<r2<<endl;
	    continue;
	}

	nvf = subpp(vf0,vf1);
	  /*if (flag){
	   cout <<"  ---  neib="<<kx<<","<<ky<<","<<kz<<endl;
	   cout <<"       ind = "<<ind<<endl;
 	   cout <<"      oVec Unit=("<<x1<<","<<y1<<","<<z1<<")"<<endl;
	   cout <<"          dtpd0="<<dtpd0.x<<","<<dtpd0.y<<","<<dtpd0.z<<endl;
	   cout <<"  nVec:  vf0="<<vf0.x<<","<<vf0.y<<","<<vf0.z<<endl;
	   cout <<"         vf1="<<vf1.x<<","<<vf1.y<<","<<vf1.z<<endl;
	   cout <<"         nvf="<<nvf.x<<","<<nvf.y<<","<<nvf.z<<endl;
	  } */

	Vector nVec(nvf.x,nvf.y,nvf.z); //new vector for A component of Procrustean's A=U*B
        //----- project the nVec to the plane perpendular to the warped PD vector
        //r1 = PD->getVectorAt(col+kx,row+ky,slice+kz, dtpd0); //get the PD of the tensor
	// HERE .....

  	//if (r1 == 0) {
	   //cout <<" ####  THIS SHOULD NOT HAPPEN  AAA !!!!  #### "<< endl;
	   //continue;
	//}
        //dtpd0.x -= col+kx; dtpd0.y -= row+ky; dtpd0.z -= slice+kz;
	//if ((dtpd0.x==0) &&(dtpd0.y==0)&&( dtpd0.z==0)){
	   //if (flag) cout <<"  PD is ZERO !!"<<endl;
	   //cout <<" ####  THIS SHOULD NOT HAPPEN  BBB !!!!  #### "<< endl; //why cannot happen?
	   //continue;
 	//}
	//Vector PDasNorm(dtpd0.x, dtpd0.y, dtpd0.z);
	//Plane pl(fPoint(0,0,0), PDasNorm);
	Vector nnVec = pl.projection(nVec);
	//float xt1,yt1,zt1; 
        nnVec.getXYZ(xt1,yt1,zt1); // warped sample's projection on the plane perpendicular to the warped PD1
	if (xt1 == 0 && yt1 == 0 && zt1 == 0){
	   cout << " (2) nnVec === 0 !!!"<<endl;
	    cout <<"    nVec="<<nvf.x<<","<<nvf.y<<","<<nvf.z<<endl;
	    cout <<"    PD="<<dtpd0.x<<","<<dtpd0.y<<","<<dtpd0.z<<endl;
	    continue;
	}
        //-----------------
  	nnVec.norm();
	nnVec.getXYZ(x2,y2,z2); //value after unit PD2 transformed 
				//and projected to the plane perpendular to PD
	  /*if (flag){		
	     cout <<"      nVec Unit=("<<x2<<","<<y2<<","<<z2<<")"<<endl;
	     cout <<"      matrix index (nAB)="<<nAB<<endl;
	  }*/

 	B[1][nAB] = xx1;  B[2][nAB] = yy1;  B[3][nAB] = zz1;
   	A[1][nAB] = x2;  A[2][nAB] = y2;  A[3][nAB] = z2;

  	nAB++;
    }
  } //cout <<"--###--"<<endl;
  nAB --; 
  // A[][1] & B[][1] are the 1st defined columns
  // A[][nAB] & B[][nAB] are the last defined columns
    //if (flag) cout <<" -- number 4 PD 2 in E = "<<ndt <<endl;
}

int mVolume::weightedShapeAt(int col, int row, int slice,VoxMap* &wV, 
	 		DTVolume* TensorVol , int *records, /* Ellipsoid& E, */ int& neib,
			float **F, float *eVal, float **eVec) //, float **T)
{// return 1 when E is properly assigned, and 'neib' defines the bound box
 // return 0 when E is meaningless or tensor is 0
  #define RADIUS 2.0f  //the neighbor simulates a radius 2 ball
  #define MAXNEIB 4    //maximal neighbor length

  //float **F, *eVal, **eVec, **T;
  float Max, Max2;
  float axLen[3], t,tOld; //threshhold : t & tOld
  int r, ndt, Condition_Hold =1; //iMax2,iMax
  int i,kx,ky,kz,loop,ind;
  double alpha;
  DTensor dtensor,dt0, dt1;
  Ellipsoid E;

  //F = matrix (1,3, 1, 3);
  //T = matrix (1,3, 1, 3);
  //eVal = vector(1,3);
  //eVec = matrix(1,3, 1, 3);

  r = TensorVol->getTensorAt(col,row,slice,dt0);
  if (r == 0) 
	return 0; // Ellipsoid undefined
  loop = 0;
  tOld = -1;
      //printTensorAt("\n~~~~~~~~~~~~~~~\n\n  dt0:\n",dt0);
  while (Condition_Hold){
      //cout <<"============================="<<endl;
      //cout << " Loop = "<< loop << endl;
    if ((dt0.a[0]==0)&&(dt0.a[1]==0)&&(dt0.a[2]==0)&&
        (dt0.a[3]==0)&&(dt0.a[4]==0)&&(dt0.a[5]==0))
	return 0; // if this is a ZERO tensor, stop

    loop ++;  // remember iteration loop, under control if not convergent

    F[1][1]= dt0.a[0]; F[1][2] = dt0.a[3]; F[1][3] = dt0.a[4];
    F[2][1]= dt0.a[3]; F[2][2] = dt0.a[1]; F[2][3] = dt0.a[5];
    F[3][1]= dt0.a[4]; F[3][2] = dt0.a[5]; F[3][3] = dt0.a[2];

    eigOfTensor(F,eVal,eVec,3);//get eigenValue & vector

       //seeMatrix(eVec,3,3,"eigen vectors:");

    // find Primary Direction (PD)
    Max = fabs(eVal[1]) ;
    //iMax = 1;
    for (i = 1; i<= 3; i++){  // find Primary Direction
	axLen[i-1] = fabs(eVal[i]); // the axis length of the ellipsoid
	if (fabs(eVal[i])>Max){
	   Max=fabs(eVal[i]);
	   //iMax = i;
	}
    }
    if (Max == 0){// this should not happen; don't use (Max <ZERO) because some is of 1e-6
	cout << "mVolume::weightedShapeAt(): ALERT!! ALERT!! "<<endl;
	exit(0);
    }
    // normalize to 1.0 scale
    for (i = 0; i< 3; i++){  
      axLen[i] /= Max;
      if (axLen[i] < ZERO)
	  axLen[i] = 0.01;
    }    
    //calculate the size of ellipsoid according to V=4.0/3.0*PI*a*b*c = 4.0/3.0*PI *((RADIUS)^3)
    //so, alpha is for new a'=a*alpha , and (a',b',c') is of the actual ellipsoid
    alpha = RADIUS /  pow ( axLen[0] * axLen[1] * axLen[2], 0.33333333f);
    Max = -1;
    for (i = 0; i< 3; i++){  
      axLen[i] *= alpha; // actual axis length; longest to be 1.0
      if (axLen[i] > MAXNEIB){ //not too long
	  //cout << " -Axis exceed:";
	  //cout << axLen[0] <<" ,"<<axLen[1]<<" ,"<<axLen[2]<<endl;
	  axLen[i] = MAXNEIB;
      }else if (axLen[i] <0.1) //not too short
	  axLen[i] = 0.1;
      if (axLen[i] > Max) //find the largest one to decide the neighbor area 
	 Max = axLen[i];
	 //if (col==104 && row==156 && slice ==100)
 	   //cout <<"Axis "<< i+1 <<" length = "<<axLen[i]<<"  (eVal="<<eVal[i+1]<<")"<<endl;
    }    
    neib = (int)(Max + 0.5); // cout << "  Max ="<<Max<<"; neib = "<<neib<<endl;
    E.setParameters(axLen[0],axLen[1],axLen[2],col, row, slice,
		 Vector ( eVec[1][1],eVec[2][1],eVec[3][1]),
		 Vector ( eVec[1][2],eVec[2][2],eVec[3][2]),
		 Vector ( eVec[1][3],eVec[2][3],eVec[3][3]));
     
    // scan the neighbor to get an average Tensor
    for (i = 0; i< 6; i++)
	dtensor.a[i] = 0;
    ndt = 0; ind = 0;
	//cout <<"--Shape--"<<endl;
    for (kz = -neib; kz <= neib; kz++)
    for (ky = -neib; ky <= neib; ky++)
    for (kx = -neib; kx <= neib; kx++, ind ++){
        if (E.inside(col+kx,row+ky,slice+kz)){
	  records[ind] = 1;  // mark the corresponding positoin as 1 (inside)
          r = TensorVol->getTensorAt(col+kx,row+ky,slice+kz,dt1);
	  if (r == 1){
	     ndt ++;  // how many tensors in the Ellisoid shape
	     if (dt1.a[0] >= 0){
    	       for (i = 0; i< 6; i++)		
	         dtensor.a[i] += dt1.a[i];
	     } else {
    	       for (i = 0; i< 6; i++)		
	         dtensor.a[i] += -(dt1.a[i]);
	     }
	  }
	} else{
	  records[ind] = 0; // mark the corresponding positoin as 0 (outside)
	}
	  //cout <<records[ind];
    }	  //cout <<"****"<<endl;
    /*if (col==69 && row==184 && slice == 7){
	cout << "  22222 neib = "<<neib<<"; ind="<<ind<<endl;
	int ssss=0;
	for (i=0; i< ind;i++)
	   ssss += records[i];
	cout << "how many points inside Ellipsoid: "<<ssss<<"; ndt="<<ndt<<endl;
    }*/
    if (ndt > 0){
	//cout <<ndt<<" grid points in the Ellipsoid"<<endl;
      for (i = 0; i< 6; i++)
          dtensor.a[i] = dtensor.a[i]/(float)ndt;
    }
    //decide if the stop iteration:
    //step 1. scan new average tensor to scale to 1.0 so it is comparable
    //iMax = 0; 
    Max = dtensor.a[0];
    for (i = 0; i < 6 ; i ++){
	if (fabs (dtensor.a[i]) > Max){
	   Max = fabs (dtensor.a[i]);
	   //iMax = i;
	}
    }	//cout << "  333333 neib = "<<neib<<endl;

    //step 2. scan old average tensor to scale to 1.0 so it is comparable
    //iMax2 = 0;
    Max2 = dt0.a[0];
    for (i = 0; i < 6 ; i ++){
	if (fabs (dt0.a[i]) > Max2){
	   Max2 = fabs (dt0.a[i]);
	   //iMax2 = i;
	}
    }
    //step 3. compare the new and the old tensors
    t = 0;
    if ((Max != 0) && (Max2 !=0)){
      for (i = 0; i < 6 ; i ++){
	t += fabs (dtensor.a[i]/Max - dt0.a[i]/Max2);
      }
    }else if((Max != 0) && (Max2 ==0)){
      for (i = 0; i < 6 ; i ++){
	t += fabs (dtensor.a[i]/Max);
      }
    }else if ((Max == 0) && (Max2 !=0)){
      for (i = 0; i < 6 ; i ++){
	t += fabs (dt0.a[i]/Max2);
      }
    }//cout << "  4444444 neib = "<<neib<<endl;
    //step 4. decide if to stop 
    t = t/ 6.0;
    //0.01 is the threshold needs to be adjusted for different problems !!!!!
    //stop iteration if only small diff or not convergent
    if ( t < 0.01 || fabs(t-tOld) < 0.005) { 
	Condition_Hold = 0;
    }else {  //continue iteration
      if (loop > 9)
	Condition_Hold = 0;
      else
      for (i = 0; i < 6 ; i ++){
	dt0.a[i]=dtensor.a[i];
      }
    }  
	//cout << "  555555 neib = "<<neib<<endl;
    tOld = t;
     //cout << "diff 2 Tensor = "<<t<<endl;
     //printTensorAt("new dtensor = \n",dtensor);
	//cout << "  66666 neib = "<<neib<<endl;
  }
    //cout <<" ###7777###  Neib = "<< neib << endl;
  //free_matrix(F,1,3, 1, 3);
  //free_matrix(T,1,3, 1, 3);
  //free_vector(eVal,1,3);
  //free_matrix(eVec,1,3, 1, 3);
  return 1;
}
int mVolume::weightedShapeAt_withOptE(int col, int row, int slice,VoxMap* &wV, 
	 		DTVolume* TensorVol , int *records, int& neib,
			float **F, float *eVal, float **eVec, Ellipsoid& optE) //, float **T)
{// return 1 when E is properly assigned, and 'neib' defines the bound box
 // return 0 when E is meaningless or tensor is 0
 // optE: optimized Ellipsoid at (col,row,slic), hopefully it is the estimated fiber orientation
  #define RADIUS 2.0f  //the neighbor simulates a radius 2 ball
  #define MAXNEIB 4    //maximal neighbor length

  //float **F, *eVal, **eVec, **T;
  float Max, Max2;
  float axLen[3], t,tOld; //threshhold : t & tOld
  int r, ndt, Condition_Hold =1; //iMax2,iMax
  int i,kx,ky,kz,loop,ind;
  double alpha;
  DTensor dtensor,dt0, dt1;
  Ellipsoid E;

  //F = matrix (1,3, 1, 3);
  //T = matrix (1,3, 1, 3);
  //eVal = vector(1,3);
  //eVec = matrix(1,3, 1, 3);

  r = TensorVol->getTensorAt(col,row,slice,dt0);
  if (r == 0) 
	return 0; // Ellipsoid undefined
  loop = 0;
  tOld = -1;
      //printTensorAt("\n~~~~~~~~~~~~~~~\n\n  dt0:\n",dt0);
  while (Condition_Hold){
      //cout <<"============================="<<endl;
      //cout << " Loop = "<< loop << endl;
    if ((dt0.a[0]==0)&&(dt0.a[1]==0)&&(dt0.a[2]==0)&&
        (dt0.a[3]==0)&&(dt0.a[4]==0)&&(dt0.a[5]==0))
	return 0; // if this is a ZERO tensor, stop

    // this is a tensor with trace close to zero , stop: note magic number
    if ((dt0.a[0] + dt0.a[1] + dt0.a[2]) < 1e-10)
      return 0;

    loop ++;  // remember iteration loop, under control if not convergent

    F[1][1]= dt0.a[0]; F[1][2] = dt0.a[3]; F[1][3] = dt0.a[4];
    F[2][1]= dt0.a[3]; F[2][2] = dt0.a[1]; F[2][3] = dt0.a[5];
    F[3][1]= dt0.a[4]; F[3][2] = dt0.a[5]; F[3][3] = dt0.a[2];

    //eigOfTensor(F,eVal,eVec,3);//get eigenValue & vector

    //try converting to double to avoid convergence problems
    double **FD, *eValD, **eVecD;

    eValD = dvector(1,3);
    eVecD = dmatrix(1,3, 1, 3);
    FD = dmatrix(1,3, 1, 3);

    for(i=1; i<=3; i++)
      for(int i1=1; i1<=3; i1++)
	FD[i][i1] = (double) F[i][i1];

    eigOfTensorDouble(FD,eValD,eVecD,3);

    for(i=1; i<=3; i++)
      {
	eVal[i]= (double) eValD[i];
	for(int i1=1; i1<=3; i1++)
	  eVec[i][i1] = (double) eVecD[i][i1];
      }

    free_dvector(eValD, 1,3);
    free_dmatrix(eVecD, 1,3,1,3);
    free_dmatrix(FD, 1,3,1,3);

    //seeMatrix(eVec,3,3,"eigen vectors:");

    // find Primary Direction (PD)
    Max = fabs(eVal[1]) ;
    //iMax = 1;
    for (i = 1; i<= 3; i++){  // find Primary Direction
	axLen[i-1] = fabs(eVal[i]); // the axis length of the ellipsoid
	if (fabs(eVal[i])>Max){
	   Max=fabs(eVal[i]);
	   //iMax = i;
	}
    }
    if (Max == 0){// this should not happen; don't use (Max <ZERO) because some is of 1e-6
	cout << "mVolume::weightedShapeAt_withOptE(): ALERT!! ALERT!! "<<endl;
	return 0;
	//	exit(0);
    }
    // normalize to 1.0 scale
    for (i = 0; i< 3; i++){  
      axLen[i] /= Max;
      if (axLen[i] < ZERO)
	  axLen[i] = 0.01;
    }    
    //calculate the size of ellipsoid according to V=4.0/3.0*PI*a*b*c = 4.0/3.0*PI *((RADIUS)^3)
    //so, alpha is for new a'=a*alpha , and (a',b',c') is of the actual ellipsoid
    alpha = RADIUS /  pow ( axLen[0] * axLen[1] * axLen[2], 0.33333333f);
    Max = -1;
    for (i = 0; i< 3; i++){  
      axLen[i] *= alpha; // actual axis length; longest to be 1.0
    }
    optE.setParameters(axLen[0],axLen[1],axLen[2],col, row, slice,
		 Vector ( eVec[1][1],eVec[2][1],eVec[3][1]),
		 Vector ( eVec[1][2],eVec[2][2],eVec[3][2]),
		 Vector ( eVec[1][3],eVec[2][3],eVec[3][3]));
    for (i = 0; i< 3; i++){  
      if (axLen[i] > MAXNEIB){ //not too long
	  //cout << " -Axis exceed:";
	  //cout << axLen[0] <<" ,"<<axLen[1]<<" ,"<<axLen[2]<<endl;
	  axLen[i] = MAXNEIB;
      }else if (axLen[i] <0.1) //not too short
	  axLen[i] = 0.1;
      if (axLen[i] > Max) //find the largest one to decide the neighbor area 
	 Max = axLen[i];
	 //if (col==104 && row==156 && slice ==100)
 	   //cout <<"Axis "<< i+1 <<" length = "<<axLen[i]<<"  (eVal="<<eVal[i+1]<<")"<<endl;
    }    
    neib = (int)(Max + 0.5); // cout << "  Max ="<<Max<<"; neib = "<<neib<<endl;
    E.setParameters(axLen[0],axLen[1],axLen[2],col, row, slice,
		 Vector ( eVec[1][1],eVec[2][1],eVec[3][1]),
		 Vector ( eVec[1][2],eVec[2][2],eVec[3][2]),
		 Vector ( eVec[1][3],eVec[2][3],eVec[3][3]));
     
    // scan the neighbor to get an average Tensor
    for (i = 0; i< 6; i++)
	dtensor.a[i] = 0;
    ndt = 0; ind = 0;
	//cout <<"--Shape--"<<endl;
    for (kz = -neib; kz <= neib; kz++)
    for (ky = -neib; ky <= neib; ky++)
    for (kx = -neib; kx <= neib; kx++, ind ++){
        if (E.inside(col+kx,row+ky,slice+kz)){
	  records[ind] = 1;  // mark the corresponding positoin as 1 (inside)
          r = TensorVol->getTensorAt(col+kx,row+ky,slice+kz,dt1);
	  if (r == 1){
	     ndt ++;  // how many tensors in the Ellisoid shape
	     if (dt1.a[0] >= 0){
    	       for (i = 0; i< 6; i++)		
	         dtensor.a[i] += dt1.a[i];
	     } else {
    	       for (i = 0; i< 6; i++)		
	         dtensor.a[i] += -(dt1.a[i]);
	     }
	  }
	} else{
	  records[ind] = 0; // mark the corresponding positoin as 0 (outside)
	}
	  //cout <<records[ind];
    }	  //cout <<"****"<<endl;
    /*if (col==69 && row==184 && slice == 7){
	cout << "  22222 neib = "<<neib<<"; ind="<<ind<<endl;
	int ssss=0;
	for (i=0; i< ind;i++)
	   ssss += records[i];
	cout << "how many points inside Ellipsoid: "<<ssss<<"; ndt="<<ndt<<endl;
    }*/
    if (ndt > 0){
	//cout <<ndt<<" grid points in the Ellipsoid"<<endl;
      for (i = 0; i< 6; i++)
          dtensor.a[i] = dtensor.a[i]/(float)ndt;
    }
    //decide if the stop iteration:
    //step 1. scan new average tensor to scale to 1.0 so it is comparable
    //iMax = 0; 
    Max = dtensor.a[0];
    for (i = 0; i < 6 ; i ++){
	if (fabs (dtensor.a[i]) > Max){
	   Max = fabs (dtensor.a[i]);
	   //iMax = i;
	}
    }	//cout << "  333333 neib = "<<neib<<endl;

    //step 2. scan old average tensor to scale to 1.0 so it is comparable
    //iMax2 = 0;
    Max2 = dt0.a[0];
    for (i = 0; i < 6 ; i ++){
	if (fabs (dt0.a[i]) > Max2){
	   Max2 = fabs (dt0.a[i]);
	   //iMax2 = i;
	}
    }
    //step 3. compare the new and the old tensors
    t = 0;
    if ((Max != 0) && (Max2 !=0)){
      for (i = 0; i < 6 ; i ++){
	t += fabs (dtensor.a[i]/Max - dt0.a[i]/Max2);
      }
    }else if((Max != 0) && (Max2 ==0)){
      for (i = 0; i < 6 ; i ++){
	t += fabs (dtensor.a[i]/Max);
      }
    }else if ((Max == 0) && (Max2 !=0)){
      for (i = 0; i < 6 ; i ++){
	t += fabs (dt0.a[i]/Max2);
      }
    }//cout << "  4444444 neib = "<<neib<<endl;
    //step 4. decide if to stop 
    t = t/ 6.0;
    //0.01 is the threshold needs to be adjusted for different problems !!!!!
    //stop iteration if only small diff or not convergent
    if ( t < 0.01 || fabs(t-tOld) < 0.005) { 
	Condition_Hold = 0;
    }else {  //continue iteration
      if (loop > 9)
	Condition_Hold = 0;
      else
      for (i = 0; i < 6 ; i ++){
	dt0.a[i]=dtensor.a[i];
      }
    }  
	//cout << "  555555 neib = "<<neib<<endl;
    tOld = t;
     //cout << "diff 2 Tensor = "<<t<<endl;
     //printTensorAt("new dtensor = \n",dtensor);
	//cout << "  66666 neib = "<<neib<<endl;
  }
    //cout <<" ###7777###  Neib = "<< neib << endl;
  //free_matrix(F,1,3, 1, 3);
  //free_matrix(T,1,3, 1, 3);
  //free_vector(eVal,1,3);
  //free_matrix(eVec,1,3, 1, 3);
  return 1;
}


void mVolume::loadIdentity()
{// reset all matrices to Identity matrix I
  long i, nMatrix=dimX*dimY*dimZ;

  for (i=0; i<nMatrix; i++)
	mv[i].loadIdentity();

}
