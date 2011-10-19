#include "Fields.h"
#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{

  char str[200];
  int i;


  if(argc<9)
    {
      std::cout<<"Usage: dtiEigD switch_endian_input(0/1) xsize ysize zsize input_dtfile output_dir faweightingPD(0/1) addoffsetPD(0/1) \n";
      exit(0);
    }
  
  
  int Xdim= atoi(argv[2]);
  int Ydim= atoi(argv[3]);
  int Zdim= atoi(argv[4]);
  float Xres= 1.72;
  float Yres= 1.72;
  float Zres= 3.0;  
  
  TensorField dti;
  
  dti.init(Xdim,Ydim,Zdim);
  dti.setVoxelSize(0,Xres); dti.setVoxelSize(1,Yres); dti.setVoxelSize(2,Zres);
  
  int endian_be=0; endian_be= atoi(argv[1]);

  if(endian_be)
    dti.importRawFieldFromFile_BE(argv[5]);
  else
    dti.importRawFieldFromFile(argv[5]);
  std::cout<<"DTI read\n";


  ScalarField e1;
  VectorField pd1;
  ScalarField e2;
  VectorField pd2;
  ScalarField e3;
  VectorField pd3;     

  ScalarField Tr;
  ScalarField Fa;

  Tr.init(Xdim,Ydim,Zdim);
  Tr.setVoxelSize(0,Xres); Tr.setVoxelSize(1,Yres); Tr.setVoxelSize(2,Zres);
  Fa.init(Xdim,Ydim,Zdim);
  Fa.setVoxelSize(0,Xres); Fa.setVoxelSize(1,Yres); Fa.setVoxelSize(2,Zres);

  e1.init(Xdim,Ydim,Zdim);
  e1.setVoxelSize(0,Xres); e1.setVoxelSize(1,Yres); e1.setVoxelSize(2,Zres);
  e2.init(Xdim,Ydim,Zdim);
  e2.setVoxelSize(0,Xres); e2.setVoxelSize(1,Yres); e2.setVoxelSize(2,Zres);
  e3.init(Xdim,Ydim,Zdim);
  e3.setVoxelSize(0,Xres); e3.setVoxelSize(1,Yres); e3.setVoxelSize(2,Zres);
  pd1.init(Xdim,Ydim,Zdim);
  pd1.setVoxelSize(0,Xres); pd1.setVoxelSize(1,Yres); pd1.setVoxelSize(2,Zres);
  pd2.init(Xdim,Ydim,Zdim);
  pd2.setVoxelSize(0,Xres); pd2.setVoxelSize(1,Yres); pd2.setVoxelSize(2,Zres);
  pd3.init(Xdim,Ydim,Zdim);
  pd3.setVoxelSize(0,Xres); pd3.setVoxelSize(1,Yres); pd3.setVoxelSize(2,Zres);
  Tr.init(Xdim,Ydim,Zdim);
  Tr.setVoxelSize(0,Xres); Tr.setVoxelSize(1,Yres); Tr.setVoxelSize(2,Zres);
  Fa.init(Xdim,Ydim,Zdim);
  Fa.setVoxelSize(0,Xres); Fa.setVoxelSize(1,Yres); Fa.setVoxelSize(2,Zres);


  ScalarField Phi;
  ScalarField Theta;
  ScalarField Psi;
  
  Phi.init(Xdim,Ydim,Zdim);
  Phi.setVoxelSize(0,Xres); Phi.setVoxelSize(1,Yres); Phi.setVoxelSize(2,Zres);
  Theta.init(Xdim,Ydim,Zdim);
  Theta.setVoxelSize(0,Xres); Theta.setVoxelSize(1,Yres); Theta.setVoxelSize(2,Zres);

  Psi.init(Xdim,Ydim,Zdim);
  Psi.setVoxelSize(0,Xres); Psi.setVoxelSize(1,Yres); Psi.setVoxelSize(2,Zres);
  
  dti.ComputeEigD(&e1, &e2, &e3, &pd1, &pd2, &pd3, &Tr, &Fa, atoi(argv[7]), &Phi, &Theta, &Psi);
  
  //dti.ComputeEigD(&e1, &e2, &e3, &pd1, &pd2, &pd3, &Tr, &Fa, atoi(argv[7]));
  std::cout<<"Full Eigen Decomposition computed\n";

  if(atoi(argv[8])==1)
  {
    pd1.AddOffset();
    pd2.AddOffset();
    pd3.AddOffset();
  }
  
  sprintf(str,"%s/FA%d.img",argv[6],Xdim);
  Fa.exportRawFieldToFile(str);
  sprintf(str,"%s/ADC%d.img",argv[6],Xdim);
  Tr.exportRawFieldToFile(str);

  sprintf(str,"%s/e1%d.img",argv[6],Xdim);
  e1.exportRawFieldToFile(str);
  sprintf(str,"%s/e2%d.img",argv[6],Xdim);
  e2.exportRawFieldToFile(str);
  sprintf(str,"%s/e3%d.img",argv[6],Xdim);
  e3.exportRawFieldToFile(str);
  std::cout<<"Eigen values written\n";  
  
  sprintf(str,"%s/pd1%d.img",argv[6],Xdim);
  pd1.exportRawFieldToFile(str);
  sprintf(str,"%s/pd2%d.img",argv[6],Xdim);
  pd2.exportRawFieldToFile(str);
  sprintf(str,"%s/pd3%d.img",argv[6],Xdim);
  pd3.exportRawFieldToFile(str);
  std::cout<<"Eigen vectors written\n";  
  
  sprintf(str,"%s/Phi%d.img",argv[6],Xdim);
  Phi.exportRawFieldToFile(str);
  sprintf(str,"%s/Theta%d.img",argv[6],Xdim);
  Theta.exportRawFieldToFile(str);
  sprintf(str,"%s/Psi%d.img",argv[6],Xdim);
  Psi.exportRawFieldToFile(str);
  std::cout<<"Angles written\n";  
  
  return 0;

}
