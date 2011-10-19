 /*-----------------------------------------------------------
  Coordinate system:

   \Z+
    \
     \_________________X+
     |O
     |
     |
     |
     |Y+
     

  Purpose:
     Transform the diffusion tensor (DT) in the original space
  (52-136 Z-slice) to the atlas space (123 Z-slice)
     A diffusion Tensor is a matrix:
//		[a1, a4, a5]	   	[XX, XY, XZ]
//		[a4, a2, a6]   i.e.	[YX, YY, YZ]
//		[a5, a6, a3]		[ZX, ZY, ZZ]


   Make file:
     make -f mkNWDT
  
   Executable example:
     nWarpDT /nilab/meteora1/susumu/DTensor/kim_1538.d ./ori2atlas.vf -O tmp.dt -Z 57

*/


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

#include "Basics.h"
#include "GPrimitive.h"
#include "DTVolume.h"
#include "VoxMap.h"
#include "WarpVolume.h"
#include "mVolume.h"
#include "Ellipsoid.h"

//#include "DTVolume.C"
//#include "VoxMap.C"
//#include "WarpVolume.C"

float resX,resY,resZ; // resolution in X,Y,Z directions of the volume
int dimX; // dim of the entered DT image
char DTfile[200]; //DTfile file name (input)
char fPD[200];    //Primary Direction(PD) vector of DTI (input)
char fPD2[200];   //2ndPrimary Direction(PD) vector of DTI (input)
char fforwardVF[120];  //warping vector field file (input)
char freverseVF[120];  //warping vector field file (input)
char outDTfile[200];  // file name of result to be stored (output)
int destZ; // destination volume Z-slice number
int CPUs; // number of CUPs to be used, default = 1
int flipZ; // flag: whether or not to flip DT's D[x,z] and D[y,z]
int USE_OPTIMIZED_DT_VOLUME=1;  // 1: use optimized DT for warping 
				// 0: use the noisy one

VoxMap *vm=NULL;
VoxMap *vmRev=NULL;
mVolume *mSVD=NULL;
VoxMap *PD=NULL; //vector field volume of DTI primary direction
VoxMap *PD2=NULL; //vector field volume of DTI 2nd primary direction, if available
DTVolume *volDT=NULL; //pointer to a DTI volume
DTVolume *optDTI=NULL; //pointer to an optimized DTI volume

void usage(char* runname)
{
   cout << "\nUsage:\n";
   cout << "   " << runname <<" <inputDTVolume> <warpVF> <revwarpVF> <PD> -P<2nd_PDfile> -O<outfile> -X<Xdim(=Ydim)> -Z<slice> -R<rx,ry,rz> -f -N<0/1>" <<endl;
   cout << "   " << runname <<" -H" <<endl;
   cout << "   " << runname <<"  for demo" <<endl;
   
   cout << "  <inputDTVolume>: the diffusion tensor volume file to be warped" <<endl;
   cout << "  <warpVF>       : the forward map from the cropped volume to atlas\n";
   cout << "  <warpVFRev>       : the reverse map from the atlas to volume\n";
   cout << "  <PD>           : the tensor filed's Primary Direction file\n";
   cout << "  -P<2nd_PDfile> : the 2nd Primary Direction file of the DT (vector field)\n";
   cout << "                   default: NULL";
   cout << "  -O<outfile>    : resulting DT field, default: ./resultDT.dt \n";
   cout << "                   this will be in the atlas (123 slc) space"<<endl;
   cout << "  -Z<slice>      : Z-slice number of the warped volume, default: 123 \n";
   cout << "  -X<Xdim(=Ydim)>      : X or Y dim of the warped volume, default: 256 \n";								 
   cout << "  -C<int>        : number of CPU to be used, default: 1 \n";
   cout << "  -R<rx,ry,rz>   : volume resolutions in X, Y and Z, default: 1.0 \n";
   cout << "  -f             : flip the D[x,z] and D[y,z] components of DT; not flipping PD (suppose PD flipped already) \n";
   cout << "  -N<0/1>	     : use the new optimized DT for warping instead of the original noisy DT\n";
   cout << "    	     : default=1: use the new optimized DT \n";
   cout << endl <<endl;

   cout << " An example:\n";
   cout << "   " << runname <<endl;
   cout << "       to show a demo, just like running the following command line" <<endl <<endl;
   cout << "   " << runname <<" /nilab/meteora1/susumu/DTensor/kim_1538.d \n";
   cout << "        /nilab/meteora1/susumu/DTensor/kim_1538.t1b.warpped.vf \n";
   cout << "        -O dt.dt -Z 57 \n\n";

   cout << "------------------" <<endl;
   cout << " - by XU, Dongrong " <<endl;
   cout << " - last update: May 14, 2002 " <<endl;
   cout << "------------------\n\n" <<endl;
   exit(0);

}


void setDefaultValues()
{
     destZ=123; //atlasZ;
     dimX = 256;
     strcpy(DTfile,"/santorini1/Procrustean/AC01/AC01.dt0.2"); //DTfile input file name
     strcpy(fforwardVF,"/santorini1/Procrustean/AC01/AC01toAH01.vf");
     fPD[0]='\0';
     fPD2[0]='\0';
     strcpy(outDTfile,"./resultDT.dt");
     CPUs = 1;
     resX=resY=resZ = 1.0;
     USE_OPTIMIZED_DT_VOLUME=1;
     //resX=resY=0.9766; resZ = 3.0;
     //getAcPcAlignmentAdjust();
}

void getParameters(int argc, char ** argv)
{   
   extern char *optarg;
   int c;

   setDefaultValues();
   if (argc==1){
      cout << "\n\n You are in demo mode" << endl;
      cout << "  Type:"<<argv[0]<<" -H to get help" << endl <<endl;
   } else if (argc <3){
      usage(argv[0]);
      return;
   } else {
     //for (int i=0; i< argc; i++)
        // cout <<i<<": "<< argv[i] <<endl;
     strcpy(DTfile,argv[1]); //DTfile input file name
     strcpy(fforwardVF,argv[2]);
     strcpy(freverseVF,argv[3]);
     strcpy(fPD,argv[4]);

     flipZ=0;
     while ((c=getopt(argc-4,argv+4,"HP:O:Z:C:fR:N:X:")) != -1) {
      //cout << " -- --  char = " << char(c) << endl;
      switch (c) {
      case 'N': //use the new optimized DT for warping instead of the original noisy DT
		//(==1: to use optimized; ==0: to use noisy one)
	sscanf(optarg,"%d",&USE_OPTIMIZED_DT_VOLUME); 
	break;
      case 'Z': //atlas Z-number
	sscanf(optarg,"%d",&destZ); 
	//cout <<"destZ="<<destZ<<endl<<endl;
	break;
      case 'P': //input 2nd PD file (the 2nd eigen vector)
	strcpy(fPD2,optarg);
	break;
      case 'O': //define the output resulting vector field from original space to atlas space
	strcpy(outDTfile,optarg);
	break;
      case 'C': //define how many CUPs for parallel execution
	sscanf(optarg,"%d",&CPUs);  
	break;	
      case 'R': //resolution
	sscanf(optarg,"%f,%f,%f",&resX,&resY,&resZ);  
	break;	
      case 'F': //define volume resolutions
	sscanf(optarg,"%f,%f,%f",&resX,&resY,&resZ);  
	break;	
      case 'f': //flip D[x,z] D[y,z] of DT
	flipZ=1;
	break;
      case 'X':
        sscanf(optarg,"%d",&dimX);
	break;
      case 'H':
      default:
	usage(argv[0]);
	break;
      }
     }
    }
    //getAcPcAlignmentAdjust();

    cout <<"DTfile="<<DTfile<<endl;
    cout <<"fforwardVF="<<fforwardVF<<endl;
    cout <<"freverseVF="<<freverseVF<<endl;
    cout <<"PD file="<<fPD<<endl;
    cout <<"PD2 file="<<fPD2<<endl;
    cout <<"Z of warped (destZ)="<<destZ<<endl;
    cout <<"number of CPU to be used="<<CPUs<<endl;
    cout <<"outDTfile="<<outDTfile<<endl<<endl;
    cout <<"volume resolution=("<<resX<<","<<resY<<","<<resZ<<")"<<endl<<endl;
    cout <<"Z flip="<<flipZ<<endl;
    cout <<"USE_OPTIMIZED_DT_VOLUME="<<USE_OPTIMIZED_DT_VOLUME<<endl;

    //exit(0);
 
}
static void agentForGetTransMatrix()
{
  //int i=m_get_myid();
  int i=1;
  //    printf("\n agentForGetTransMatrix(): This is Processor No. %d\n", i);
    mSVD -> getTransMatrix(vm,volDT, PD, PD2, optDTI,USE_OPTIMIZED_DT_VOLUME); //generate the transformation matrix & optimized DTI
    printf("\n Task at Processor No. %d DONE! \n", i);
}

int main(int argc,char **argv) //for warp DTI
{//warp DT, this is the real Main()
   int xx,yy,zz;
   struct DTensor * nDTvol=NULL;

   getParameters(argc, argv);
 
   cout << "Step 1. loading Diffusion Tensor...." << endl << endl;
   // rv  volDT= new DTVolume(DTfile,256,256,resX,resY,resZ); //1,1,1); 
   volDT= new DTVolume(DTfile,dimX,dimX,resX,resY,resZ); //1,1,1); 
		//,256,256,0.9766,0.9766,3.0);
   volDT->getDimXYZ(xx, yy, zz);
   if (flipZ) volDT->flipZComponent();
   //volDT->printTensors();
   //return 0;

   //exit(0);
   cout << endl <<"Step 2. loading the forward warping vector field.... \n";
   cout << "         resolution of the volume not important" << endl;
   vm=new VoxMap("N/A",fforwardVF, YES, xx, yy, zz,resX,resY,resZ);
	 //1,1,1); 
	//, resX,resY,resZ); //0.9766, 0.9766, 3.0);

   cout << endl <<"Step 3. loading the reverse warping vector field.... \n";
   cout << "         resolution of the volume not important" << endl;
   vmRev=new VoxMap("N/A",freverseVF, YES, xx, yy, zz,resX,resY,resZ);

   cout << endl <<"Step 4. loading the Primary Direction (PD) vectors .... \n";
   if (fPD[0] != '\0')
     PD=new VoxMap("N/A",fPD, YES,xx, yy, zz);
   else {
     cout << " PD file not supplied. "<<endl;
     cout << " This version does not support this feature. I will now exit."<<endl;
     exit(0);
   }

   if (fPD2[0] != '\0'){
     cout << endl <<"     loading the 2nd Primary Direction (PD2) vectors .... \n";
     PD2=new VoxMap("N/A",fPD2, YES,xx, yy, zz);
   }else {
     cout << " PD2 file not supplied. "<<endl;
     PD2 = NULL;
   }

   cout << endl << "Step 5. create matrix volume mVolume: use SVD; and optimize the original DTI \n";
   mSVD = new mVolume(xx,yy,zz);

   //===============
   optDTI= new DTVolume; // create space for optimized DTI volume
   optDTI->newVolume(xx, yy, zz,resX,resY,resZ); // to be used in agentForGetTransMatrix()

   // m_set_procs(CPUs);
   // int cn= m_get_numprocs();
   int cn=1;
   // cout <<" Main(): cpu to be used = "<<cn<<endl; 
   // m_fork(agentForGetTransMatrix);  //compute transformation matrix
   //m_kill_procs();
    
     // reset all maratrices to be identity matrix I
   mSVD->loadIdentity(); //only for experiment of no reorientation, only relocation
			   // needs to be removed

   //optDTI->saveTensorData("optimized.d1"); //save the optimized DTI volume
   //exit(0);
   //================
   
   //mSVD -> getTransMatrixOld(vm,*volDT);
   //mSVD -> getTransMatrix(vm,volDT, PD);
   mSVD -> getTransMatrix(vm,volDT, PD, PD2, optDTI,USE_OPTIMIZED_DT_VOLUME);
   delete PD;

   /***** for debugging purpose ****
   Ellipsoid E;
   mSVD->weightedShapeAt(120, 17, 62,vm, volDT,E,n);
   mSVD->weightedShapeAt(126, 74, 44,vm, volDT,E,n);
   mSVD->weightedShapeAt(123, 50, 46,vm, volDT,E,n);
   cout <<" ################"<<endl<<endl;
   mSVD->weightedShapeAt(125, 122, 36,vm, volDT,E,n);
   mSVD->weightedShapeAt(125, 172, 66,vm, volDT,E,n);
   ***** for debugging purpose ****/


   //delete warp;
   cout << endl << "Step 6. warp the (optimized) DT volume from original->atlas" << endl << endl;
   //****// 

   if (USE_OPTIMIZED_DT_VOLUME){ //use optimized DT volume to replace the original noisy one
     delete volDT; //delete the old noise DT volume
     // load in the new optimized DT generated in calculating the reorientation matrix
     volDT= optDTI;
         //new DTVolume("optimized.d",256,256,resX,resY,resZ); //1,1,1); 
		//,256,256,0.9766,0.9766,3.0);
     if (flipZ) volDT->flipZComponent();
   } else{
     delete optDTI;
   }
   // mSVD->loadIdentity();
   //DTVolume *volDT;
   //getchar(); 
   volDT->reorientTensor(mSVD);
   delete mSVD;

   cout << endl << "Step 7. Cast DT to grid in Atlas" << endl << endl;
   //****//
   //nDTvol = volDT->warpTensorField(vm,destZ); 
   nDTvol = volDT->warpTensorFieldRev(vmRev,destZ); 

   delete vm;
   cout << "saving resulting DT file to <"<<outDTfile<<">...."<<endl;
   //****//
     //rv  volDT->saveData(outDTfile, nDTvol, atlasX*atlasY*destZ);
volDT->saveData(outDTfile, nDTvol, xx*yy*destZ);

   cout << "new DT volume <"<<outDTfile<<"> saved." <<endl;
   //****//
   delete []nDTvol;
   delete volDT;
   
   cout << endl << "program completed." << endl <<endl;
   cout <<endl<< "~~~~~~~~~~~~~~~~" <<endl;
   cout <<"      feed back to xdr@cbmv.jhu.edu. Thanks" <<endl;
   return 0;
}

