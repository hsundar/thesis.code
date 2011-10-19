#include "Fields.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{

  char str[200];
  int i;

  int SUB;

  if(argc<7 )
    {
      std::cout<<"Usage: AverDTI  xsize ysize zsize input_dtfiles_txt output_dt_file switch_endian_input(0/1)\n";
      exit(0);
    }
   
  int Xdim= atoi(argv[1]);
  int Ydim= atoi(argv[2]);
  int Zdim= atoi(argv[3]);

  TensorField dti;
  
  dti.init(Xdim,Ydim,Zdim);
  
  int endian_be=0; endian_be= atoi(argv[6]);

  std::vector< std::string > DTIfilenames;
  std::string DTIfilename;

  std::ifstream DTIf;
  DTIf.open(argv[4]);
  
  int k;
  
  int ix, iy, iz;

  TensorField dti_avg;
  dti_avg.init(Xdim,Ydim,Zdim);
  
  unsigned char *mask= (unsigned char *)calloc(Xdim*Ydim*Zdim, sizeof(unsigned char));

  for(k=0; !DTIf.eof(); k++)
    {
      if(DTIf.eof())
	    break;

       DTIf>>DTIfilename;
       
       DTIfilenames.push_back(DTIfilename);

      if(endian_be)
	{
	  std::cout<<"Reading "<< DTIfilename<<" in big endian format\n";
	  dti.importRawFieldFromFile_BE(DTIfilename.c_str());
	}
      else
	{
	  std::cout<<"Reading "<< DTIfilename<<" in little endian format\n";
	  dti.importRawFieldFromFile(DTIfilename.c_str());
	}

      if(k==0)
	{
	  // mask computed from first DTI
	  // note magic number 0.001
	  dti.computeNZmask(mask,0.001);
	}

      TensorField dti_tmp;
      dti_tmp.init(Xdim,Ydim,Zdim);
  
      // compute log of tensor 
      dti.logTField(&dti_tmp, mask);

      for(ix=0; ix<Xdim; ix++)
	for(iy=0; iy<Ydim; iy++)
	  for(iz=0; iz<Zdim; iz++)
	    {
	      int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	      	     
	      if(mask[ind]==0)
		{	  
		  continue;
		}
	     
	      float *val_sub= dti_tmp.getAt(ix,iy,iz);

	      float *val_avg= dti_avg.getAt(ix,iy,iz);

	      for(i=0; i<6; i++)
		val_avg[i] += val_sub[i];
	      
	    }
    }

  DTIf.close();

  SUB = k;
  cout<<SUB<<" subjects read\n";

  for(ix=0; ix<Xdim; ix++)
    for(iy=0; iy<Ydim; iy++)
      for(iz=0; iz<Zdim; iz++)
	{
	  int ind = iz*Xdim*Ydim + iy*Xdim + ix;
	  
	  if(mask[ind]==0)
	    {	  
	      continue;
	    }
	  
	  float *val_avg= dti_avg.getAt(ix,iy,iz);
	  
	  for(i=0; i<6; i++)
	    val_avg[i] /= SUB;
	  
	}
  
  TensorField dti_tmp;
  dti_tmp.init(Xdim,Ydim,Zdim);
  
  dti_avg.expTField(&dti_tmp, mask);

  sprintf(str,"%s",argv[5]);
  dti_tmp.exportRawFieldToFile(str);

  return 0;

}
