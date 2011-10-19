#include "Fields.h"
#include <iostream>
using namespace std;

#include "RGBAcolorMap.h"

int main(int argc, char *argv[])
{

  char str[200];
  int i;

  int SAVE_INT_RESULTS=0;

  if(argc<7 )
    {
      std::cout<<"Usage: PDcolorMap input_PDfile xsize ysize zsize output_cmpfile OffsetFlag(0/1)\n";
      exit(0);
    }
    
  int Xdim= atoi(argv[2]);
  int Ydim= atoi(argv[3]);
  int Zdim= atoi(argv[4]);
  
  float Xres = 1.72f;
  float Yres = 1.72f;
  float Zres = 3.0f;

  char input_PDfile[200];
  char output_cmpfile[200];
  
  int OffsetFlag = atoi(argv[6]);

  strcpy(input_PDfile, argv[1]);
  strcpy(output_cmpfile, argv[5]);
  
  VectorField PD;  
  PD.init(Xdim,Ydim,Zdim);
  PD.setVoxelSize(0,Xres); PD.setVoxelSize(1,Yres); PD.setVoxelSize(2,Zres);
  
  PD.importRawFieldFromFile(input_PDfile);
  std::cout<<"PD File read\n";

  if(OffsetFlag)
    PD.SubtractOffset();
  
  RGBAcolorMap *vRGBA= new RGBAcolorMap(Zdim, Xdim, Ydim);

  float *pd_val;
  struct byteRGBA aMap;
  int ix, iy, iz, n;

  for(iz=0; iz<Zdim; iz++)
    {
      cout<<"        working on slice Z = "<<iz<<endl;
      for(ix=0; ix<Xdim; ix++)
	for(iy=0; iy<Ydim; iy++, n++)
	  {
	    aMap.rgba[0]= aMap.rgba[1]= aMap.rgba[2]= aMap.rgba[3]=0;

	    pd_val = PD.getAt(ix,iy,iz);
	    
	    float d, a =fabs(pd_val[0]), b = fabs(pd_val[1]), c = fabs(pd_val[2]);

	    d = a;
	    
	    if (b>d) d = b;
	    if (c>d) d = c;
	    
	    if (d >1) 
	      {
		cout << "("<<ix<<","<<iy<<","<<iz<<"): PD component >1.0:"<<d<<endl;
	      }

	    if (d > 0)
	      {
		aMap.rgba[3]=255.0f*a;  //R
		aMap.rgba[2]=255.0f*b;  //G
		aMap.rgba[1]=255.0f*c;  //B
		
		aMap.rgba[0]= 255.0f;
	      }

	    vRGBA->setVoxel(ix, iy, iz, aMap);
	  }
    }
	
  vRGBA->saveVolume(output_cmpfile);
  cout <<" Result saved."<<endl<<endl;
  delete vRGBA;
  
  return 0;

}
