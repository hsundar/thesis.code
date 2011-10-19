#include "Fields.h"
#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{

  char str[100];
  int i;


  if(argc<7)
    {
      std::cout<<"Usage: createInterleavedDTI switch_endian_input(0/1) xsize ysize zsize input_dir output_dtfile\n";
      exit(0);
    }
  
  
  int Xdim= atoi(argv[2]);
  int Ydim= atoi(argv[3]);
  int Zdim= atoi(argv[4]);
  
  /*---
  float Xres= 1.72;
  float Yres= 1.72;
  float Zres= 3.0;  
  --*/

  TensorField dti;
  
  dti.init(Xdim,Ydim,Zdim);
  //dti.setVoxelSize(0,Xres); dti.setVoxelSize(1,Yres); dti.setVoxelSize(2,Zres);
  
  int endian_be=0; endian_be= atoi(argv[1]);

  ScalarField Dxx, Dyy, Dzz, Dxy, Dxz, Dyz;
  Dxx.init(Xdim,Ydim,Zdim);
  Dyy.init(Xdim,Ydim,Zdim);
  Dzz.init(Xdim,Ydim,Zdim);
  Dxy.init(Xdim,Ydim,Zdim);
  Dxz.init(Xdim,Ydim,Zdim);
  Dyz.init(Xdim,Ydim,Zdim);

  sprintf(str,"%s/Dxx.img",argv[5]);
  if(endian_be)
    Dxx.importRawFieldFromFile_BE(str);
  else
    Dxx.importRawFieldFromFile(str);
  std::cout<<"Dxx read from "<<str<<"\n";

  sprintf(str,"%s/Dyy.img",argv[5]);
  if(endian_be)
    Dyy.importRawFieldFromFile_BE(str);
  else
    Dyy.importRawFieldFromFile(str);
  std::cout<<"Dyy read from "<<str<<"\n";

  sprintf(str,"%s/Dzz.img",argv[5]);
  if(endian_be)
    Dzz.importRawFieldFromFile_BE(str);
  else
    Dzz.importRawFieldFromFile(str);
  std::cout<<"Dzz read from "<<str<<"\n";

  sprintf(str,"%s/Dxy.img",argv[5]);
  if(endian_be)
    Dxy.importRawFieldFromFile_BE(str);
  else
    Dxy.importRawFieldFromFile(str);
  std::cout<<"Dxy read from "<<str<<"\n";

  sprintf(str,"%s/Dxz.img",argv[5]);
  if(endian_be)
    Dxz.importRawFieldFromFile_BE(str);
  else
    Dxz.importRawFieldFromFile(str);
  std::cout<<"Dxz read from "<<str<<"\n";

  sprintf(str,"%s/Dyz.img",argv[5]);
  if(endian_be)
    Dyz.importRawFieldFromFile_BE(str);
  else
    Dyz.importRawFieldFromFile(str);
  std::cout<<"Dyz read from "<<str<<"\n";
  
  dti.createInterleavedDTI(&Dxx, &Dyy, &Dzz, &Dxy, &Dxz, &Dyz);
  std::cout<<"DTI computed\n";

  dti.exportRawFieldToFile(argv[6]);
  std::cout<<"DTI written\n";

  return 0;

}
