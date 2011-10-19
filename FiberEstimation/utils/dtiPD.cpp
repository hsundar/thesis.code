#include "Fields.h"
#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{

  char str[200];
  int i;


  if(argc<10)
    {
      std::cout<<"Usage: dtiPD dirn(1/2/3) xsize ysize zsize input_dtfile output_pdfile faweighting(0/1) addoffset(0/1) switch_endian_input(0/1)\n";
      exit(0);
    }
  
  int dirn = atoi(argv[1]);
  if(dirn<1 || dirn >3)
    {
      cout<<"Dirn should be 1/2/3"<<endl;
    }
  else
    dirn -= 1;
  int Xdim= atoi(argv[2]);
  int Ydim= atoi(argv[3]);
  int Zdim= atoi(argv[4]);

  // these values do not matter
  float Xres= 1.72;
  float Yres= 1.72;
  float Zres= 3.0;  
  
  TensorField dti;
  
  dti.init(Xdim,Ydim,Zdim);
  dti.setVoxelSize(0,Xres); dti.setVoxelSize(1,Yres); dti.setVoxelSize(2,Zres);
  
  int endian_be=0; endian_be= atoi(argv[9]);

  if(endian_be)
    dti.importRawFieldFromFile_BE(argv[5]);
  else
    dti.importRawFieldFromFile(argv[5]);
  std::cout<<"DTI read\n";

  VectorField pd;
  
  pd.init(Xdim,Ydim,Zdim);
  pd.setVoxelSize(0,Xres); pd.setVoxelSize(1,Yres); pd.setVoxelSize(2,Zres);
  dti.ComputePD(&pd, dirn, atoi(argv[7]));
  std::cout<<"PD computed\n";

  if(atoi(argv[8])==1)
  {
    pd.AddOffset();


  }
  
  sprintf(str,"%s",argv[6]);
  pd.exportRawFieldToFile(str);
  std::cout<<"PD written\n";  

  return 0;

}
