#include "Fields.h"
#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{

  char str[200];
  int i;


  if(argc<7 )
    {
      std::cout<<"Usage: dtiTrace  xsize ysize zsize input_dtfile output_trace_file switch_endian_input(0/1)\n";
      exit(0);
    }
  
 
  int Xdim= atoi(argv[1]);
  int Ydim= atoi(argv[2]);
  int Zdim= atoi(argv[3]);
  
  TensorField dti;
  
  dti.init(Xdim,Ydim,Zdim);
  
  int  endian_be=0; endian_be= atoi(argv[6]);

  if(endian_be)
    dti.importRawFieldFromFile_BE(argv[4]);
  else
    dti.importRawFieldFromFile(argv[4]);
  std::cout<<"DTI read\n";

  //std::cout<<"DTI read\n";

  ScalarField trace;
  trace.init(Xdim,Ydim,Zdim);
  // trace.setVoxelSize(0,Xres); trace.setVoxelSize(1,Yres); trace.setVoxelSize(2,Zres);
  dti.ComputeTrace(&trace);
  //std::cout<<"Smooth TRACE computed\n";
  
  sprintf(str,"%s",argv[5]);
  trace.exportRawFieldToFile(str);
  
  return 0;

}
