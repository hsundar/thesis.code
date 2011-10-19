#include "Fields.h"
#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{

  char str[200];
  int i;


  if(argc<11 )
    {
      std::cout<<"Usage: smoothDTI  sigma xsize ysize zsize input_dtfile output_dtfile Xres Yres Zres switch_endian_input(0/1)\n";
      exit(0);
    }
  
  float sigma= atof(argv[1]);
  int Xdim= atoi(argv[2]);
  int Ydim= atoi(argv[3]);
  int Zdim= atoi(argv[4]);
  float Xres= atof(argv[7]);
  float Yres= atof(argv[8]);
  float Zres= atof(argv[9]);  

  int endian_be= 0;
  endian_be=atoi(argv[10]);

  TensorField dti;
  
 
  // note magic number 5.0:

  int sXdim, sYdim, sZdim;
  sXdim = 5.0/Xres*sigma;
  sYdim = 5.0/Yres*sigma;
  sZdim = 5.0/Zres*sigma;
  
  if(sXdim%2==0) sXdim +=1;
  if(sYdim%2==0) sYdim +=1;
  if(sZdim%2==0) sZdim +=1;
  
  if(sXdim > Xdim) sXdim = Xdim;
  if(sYdim > Ydim) sYdim = Ydim;
  if(sZdim > Zdim) sZdim = Zdim;

  std::cout<<"Size of smoothing filter: "<<sXdim<<"x"<<sYdim<<"x"<<sZdim<<"\n";


  dti.init(Xdim,Ydim,Zdim);
  dti.setVoxelSize(0,Xres); dti.setVoxelSize(1,Yres); dti.setVoxelSize(2,Zres);
  if(endian_be)
    dti.importRawFieldFromFile_BE(argv[5]);
  else
    dti.importRawFieldFromFile(argv[5]);
  std::cout<<"DTI read\n";

  ScalarField smoothFilter;
  
 
  smoothFilter.init(sXdim,sYdim,sZdim);
  smoothFilter.setVoxelSize(0,Xres); 
  smoothFilter.setVoxelSize(1,Yres); 
  smoothFilter.setVoxelSize(2,Zres);
 
  smoothFilter.ComputeSmoothingFilter(sigma);
    
  TensorField dti2 ;
  dti2.init(Xdim,Ydim,Zdim);
  dti2.setVoxelSize(0,Xres); dti2.setVoxelSize(1,Yres); dti2.setVoxelSize(2,Zres);
  
  dti.logESmooth(&smoothFilter, &dti2);
  
  std::cout<<"Log-Euclidean Smoothing done\n";  
  
  sprintf(str,"%s",argv[6]);
  dti2.exportRawFieldToFile(str);
  
  ScalarField fa;
  fa.init(Xdim,Ydim,Zdim);
  fa.setVoxelSize(0,Xres); fa.setVoxelSize(1,Yres); fa.setVoxelSize(2,Zres);
  dti2.ComputeFA(&fa);
  //std::cout<<"Smooth FA computed\n";
  
  sprintf(str,"%s.FA.img",argv[6]);
  fa.exportRawFieldToFile(str);
  
  return 0;

}
