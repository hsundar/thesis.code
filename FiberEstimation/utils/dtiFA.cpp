#include "Fields.h"
#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{

  char str[200];
  int i;


  if(argc<7 )
    {
      std::cout<<"Usage: dtiFA  xsize ysize zsize input_dtfile output_fa_file switch_endian_input(0/1)\n";
      exit(0);
    }
  
 
  int Xdim= atoi(argv[1]);
  int Ydim= atoi(argv[2]);
  int Zdim= atoi(argv[3]);
  
  TensorField dti;
  
  dti.init(Xdim,Ydim,Zdim);
  
  int endian_be=0; endian_be= atoi(argv[6]);

  if(endian_be)
    {
      //std::cout<<"Reading in big endian format\n";
      dti.importRawFieldFromFile_BE(argv[4]);
    }
  else
    {
      dti.importRawFieldFromFile(argv[4]);
    }

  //std::cout<<"DTI read\n";

  ScalarField fa;
  fa.init(Xdim,Ydim,Zdim);
  // fa.setVoxelSize(0,Xres); fa.setVoxelSize(1,Yres); fa.setVoxelSize(2,Zres);
  dti.ComputeFA(&fa);
  //std::cout<<"Smooth FA computed\n";
  
  sprintf(str,"%s",argv[5]);
  fa.exportRawFieldToFile(str);
  
  return 0;

}
